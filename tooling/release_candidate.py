"""Retain local acquisition assertions separately from the release candidate."""

import argparse
import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path

from research_repo_tools.candidate_publication import plan_candidate_publication
from research_repo_tools.candidate_reuse import (
    CandidatePolicy,
    acquisition_assertion_digest,
    parse_candidate_policy,
    plan_candidate_archive,
    serialize_candidate_policy,
    verify_candidate_archive,
    verify_candidate_release,
)
from research_repo_tools.cli import main as shared_main
from research_repo_tools.common_measurement import validate_run_policy
from research_repo_tools.complete_runs import run_identity
from research_repo_tools.config import load as load_config
from research_repo_tools.dependency_provenance import check_release_dependencies
from research_repo_tools.evidence import parse_evidence
from research_repo_tools.files import replace_many
from research_repo_tools.measurement import load_measurement, resolve_revision
from research_repo_tools.performance_progress import ProgressReporter
from research_repo_tools.performance_release import load_performance_release, load_preparation, measure_release, serialize_preparation
from research_repo_tools.publication import publish_publication
from research_repo_tools.release_assets import (
    ReleaseTarget,
    download_release_asset,
    parse_release_target,
    preflight_release_target,
    publish_release_asset,
    revalidate_release_target,
    serialize_release_target,
)
from research_repo_tools.release_interfaces import load_evidence_policy
from research_repo_tools.release_publishing import ReviewedRelease, check_reviewed_release
from research_repo_tools.release_validation import api
from research_repo_tools.scheduled_runs import ScheduledRun
from research_repo_tools.suite_runs import parse_suite_run

SETTINGS = "tooling/performance-release.toml"
REPOSITORY = "acgetchell/la-stack"
ROLES = ("policy", "publication", "acceptance", "archive")


def owner_directory(root: Path, checkpoint: str, selected: Path | None = None) -> Path:
    """Keep the operator's immutable record outside candidate-controlled files."""
    identity = hashlib.sha256(str((root / checkpoint).resolve()).encode("utf-8")).hexdigest()
    directory = (selected or Path.home() / ".local/state/la-stack/acquisitions" / identity).expanduser().resolve()
    if directory.is_relative_to(root) or root.is_relative_to(directory):
        message = "acquisition owner directory must be separate from the checkout"
        raise ValueError(message)
    return directory


def freeze(directory: Path, files: dict[str, bytes]) -> None:
    """Atomically retain owner assertions; retries cannot replace earlier bytes."""
    directory = directory.resolve(strict=True)
    updates = {directory / name: data for name, data in files.items()}
    if any(not path.resolve().is_relative_to(directory) for path in updates):
        message = "owner files must remain within their record directory"
        raise ValueError(message)

    def unchanged() -> None:
        for path, data in updates.items():
            if path.exists() and path.read_bytes() != data:
                raise ValueError(f"immutable output differs: {path.name}")

    unchanged()
    replace_many(updates, verify=unchanged)


def capture(root: Path, owner: Path, *, retained_scratch: bool = False) -> dict[str, str]:
    """Finish the trusted owner's receipt from its already retained original bytes."""
    settings = load_performance_release(root, SETTINGS)
    preparation = load_preparation(root, settings)
    if serialize_preparation(preparation) != (owner / "preparation.json").read_bytes():
        message = "independent preparation changed since acquisition"
        raise ValueError(message)
    if retained_scratch:
        # Explicit owner recovery after a failed callback; never infer this
        # authority from an incoming artifact or silently fall back to scratch.
        preparation.verify(root)
        freeze(owner, {"run.json": (root / settings.payload).read_bytes(), "evidence.json": (root / settings.manifest).read_bytes()})
    payload, manifest = (owner / "run.json").read_bytes(), (owner / "evidence.json").read_bytes()
    # Restore only identical/missing scratch evidence for a failed capture retry.
    freeze(root, {settings.payload: payload, settings.manifest: manifest})
    evidence = parse_evidence(payload, manifest)
    run = validate_run_policy(evidence, load_measurement(root, settings.measurement), release=True)
    if not isinstance(run, ScheduledRun):
        message = "candidate acquisition requires the declared trial schedule"
        raise TypeError(message)
    policy = CandidatePolicy(
        settings.repository,
        preparation.pair.current,
        None,
        "local",
        run_identity(evidence),
        run.plan,
        dict(evidence.sources)["current"],
        tuple(
            (index, suite, phase, source) for index, trial in enumerate(run.trials, 1) for suite, phase, _sample, source in parse_suite_run(trial).observations
        ),
        report=settings.report,
    )
    acquisition = acquisition_assertion_digest(policy)
    freeze(owner, {"policy.json": serialize_candidate_policy(policy), "acquisition.sha256": (acquisition + "\n").encode("ascii")})
    receipt_path = owner / "publication.json"
    if not receipt_path.exists():
        _plan, receipt = plan_candidate_publication(root, settings.payload, settings.manifest, policy, acquisition_sha256=acquisition)
        freeze(owner, {"publication.json": receipt})
    publication = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
    freeze(owner, {"publication.sha256": (publication + "\n").encode("ascii")})
    return {"owner_directory": str(owner), "acquisition_sha256": acquisition, "publication_sha256": publication}


def acquire(root: Path, owner: Path, checkpoint: str, *, resume: bool, progress: ProgressReporter) -> dict[str, str]:
    """Acquire under owner control, freezing assertions before report promotion."""
    settings = load_performance_release(root, SETTINGS)
    preparation = load_preparation(root, settings)
    check_release_dependencies(root)
    if resume:
        if not owner.is_dir():
            message = "resume requires the original acquisition owner record"
            raise ValueError(message)
    else:
        owner.mkdir(mode=0o700, parents=True)
    freeze(owner, {"preparation.json": serialize_preparation(preparation)})

    def before_promotion() -> None:
        freeze(owner, {"run.json": (root / settings.payload).read_bytes(), "evidence.json": (root / settings.manifest).read_bytes()})
        capture(root, owner)

    measure_release(
        root,
        settings,
        allow_git_mutations=True,
        checkpoint=checkpoint,
        resume=resume,
        progress=progress,
        before_promotion=before_promotion,
        release_intent=True,
    )
    return {
        "owner_directory": str(owner),
        "acquisition_sha256": (owner / "acquisition.sha256").read_text(encoding="ascii").strip(),
        "publication_sha256": (owner / "publication.sha256").read_text(encoding="ascii").strip(),
    }


def final_directory(owner: Path, commit: str, selection: str) -> Path:
    """Name an explicit immutable review slot; never replace an earlier selection."""
    if re.fullmatch(r"[0-9a-f]{40}", commit) is None or re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", selection) is None:
        message = "final evidence requires a full commit and a short lowercase review selection name"
        raise ValueError(message)
    return owner / "final" / commit / selection


def verify(root: Path, owner: Path, selection: str = "reviewed") -> Path:
    """Bind the separately trusted record to clean HEAD and signed native CI."""
    check_release_dependencies(root)
    commit = resolve_revision(root, "HEAD")
    policy = replace(parse_candidate_policy((owner / "policy.json").read_bytes()), commit=commit)
    final = final_directory(owner, commit, selection)
    final.mkdir(parents=True, exist_ok=True)
    freeze(final, {"policy.json": serialize_candidate_policy(policy)})
    accepted = final / "acceptance.json"
    with tempfile.TemporaryDirectory(prefix="la-stack-candidate-verification-") as temporary:
        output = Path(temporary).resolve() / "acceptance.json"
        arguments = [
            "--root",
            str(root),
            "release",
            "evidence",
            "candidate",
            "tooling/release-evidence.toml",
            "--checkout",
            str(root),
            "--policy",
            str(final / "policy.json"),
            "--payload",
            str(owner / "run.json"),
            "--manifest",
            str(owner / "evidence.json"),
            "--publication-receipt",
            str(owner / "publication.json"),
            "--acquisition-sha256",
            (owner / "acquisition.sha256").read_text(encoding="ascii").strip(),
            "--publication-sha256",
            (owner / "publication.sha256").read_text(encoding="ascii").strip(),
            "--output",
            str(output),
            "--format",
            "summary",
        ]
        if accepted.exists():
            arguments.extend(("--acceptance", str(accepted)))
        status = shared_main(arguments)
        if status:
            raise RuntimeError(f"candidate verification failed with exit {status}")
        freeze(final, {"acceptance.json": output.read_bytes()})
    return accepted


def archive(owner: Path, commit: str, selection: str = "reviewed") -> Path:
    """Derive an immutable release archive without acquisition or native work."""
    final = final_directory(owner, commit, selection)
    policy = parse_candidate_policy((final / "policy.json").read_bytes())
    if policy.commit != commit:
        message = "selected candidate policy differs from the independently selected commit"
        raise ValueError(message)
    prefix = final.relative_to(owner).as_posix()
    output = f"{prefix}/la-stack-{policy.tag}-criterion-baseline.tar.gz"
    publish_publication(plan_candidate_archive(owner, "run.json", "evidence.json", f"{prefix}/policy.json", f"{prefix}/acceptance.json", output))
    verify_candidate_archive(owner / output, policy, (final / "acceptance.json").read_bytes())
    return owner / output


@dataclass(frozen=True, slots=True)
class Transfer:
    """A parsed owner selection; only an authorized invocation authenticates it."""

    target: ReleaseTarget
    acquisition: str
    publication: str
    assets: tuple[tuple[str, str, int, str], ...]


def transfer_asset(role: str, raw: object) -> tuple[str, str, int, str]:
    """Parse one exact inert asset selection."""
    if not isinstance(raw, dict) or set(raw) != {"name", "id", "sha256"}:
        message = "invalid candidate asset assertion"
        raise ValueError(message)
    name, identifier, digest = raw["name"], raw["id"], raw["sha256"]
    if not isinstance(name, str) or re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9._-]{0,199}", name) is None or type(identifier) is not int or identifier <= 0:
        message = "candidate assets require safe names and positive exact IDs"
        raise ValueError(message)
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        message = "candidate asset requires a full SHA-256 assertion"
        raise ValueError(message)
    return role, name, identifier, digest


def parse_transfer(payload: bytes) -> Transfer:
    """Reject malformed or ambiguous operator assertions before provider access."""

    def unique(pairs: list[tuple[str, object]]) -> dict[str, object]:
        if len(dict(pairs)) != len(pairs):
            message = "duplicate transfer field"
            raise ValueError(message)
        return dict(pairs)

    data = json.loads(payload, object_pairs_hook=unique)
    if (
        not isinstance(data, dict)
        or set(data) != {"schema", "target", "acquisition", "publication", "assets"}
        or type(data["schema"]) is not int
        or data["schema"] != 1
    ):
        message = "invalid candidate transfer schema"
        raise ValueError(message)
    target = parse_release_target(json.dumps(data["target"]).encode("utf-8"))
    if target.repository != REPOSITORY or target.expected_title != target.tag:
        message = "candidate transfer must select this repository and its exact release title"
        raise ValueError(message)

    def digest(value: object) -> str:
        if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
            message = "candidate transfer requires full SHA-256 assertions"
            raise ValueError(message)
        return value

    raw_assets = data["assets"]
    if not isinstance(raw_assets, dict) or set(raw_assets) not in (set(ROLES), {*ROLES, "attached"}):
        message = "candidate transfer requires four original assets and an optional attached archive"
        raise ValueError(message)
    assets = tuple(transfer_asset(role, raw_assets[role]) for role in (*ROLES, *(["attached"] if "attached" in raw_assets else [])))
    if len({item[1] for item in assets}) != len(assets) or len({item[2] for item in assets}) != len(assets):
        message = "candidate assets must have distinct names and IDs"
        raise ValueError(message)
    if "attached" in raw_assets and (raw_assets["attached"]["name"] != f"la-stack-{target.tag}-criterion-baseline.tar.gz" or assets[-1][3] != assets[3][3]):
        message = "canonical attached archive must match the originally accepted bytes"
        raise ValueError(message)
    return Transfer(target, digest(data["acquisition"]), digest(data["publication"]), assets)


def stage(root: Path, owner: Path, selection: str) -> Path:
    """Upload inert original records, retaining dispatch assertions independently."""
    accepted = verify(root, owner, selection)
    commit = resolve_revision(root, "HEAD")
    final = final_directory(owner, commit, selection)
    policy = parse_candidate_policy((final / "policy.json").read_bytes())
    artifact = archive(owner, commit, selection)
    target_path = final / "target.json"
    if not target_path.exists():
        target = preflight_release_target(root, REPOSITORY, policy.tag, commit, expected_title=policy.tag, checkout=root)
        freeze(final, {"target.json": serialize_release_target(target)})
    target = parse_release_target(target_path.read_bytes())
    if (target.repository, target.tag, target.commit) != (REPOSITORY, policy.tag, commit):
        message = "staging target differs from the accepted candidate"
        raise ValueError(message)
    revalidate_release_target(root, target, checkout=root)
    inputs = {"policy": final / "policy.json", "publication": owner / "publication.json", "acceptance": accepted, "archive": artifact}
    transport = final / "transport"
    transport.mkdir(exist_ok=True)
    assertions: dict[str, dict[str, str | int]] = {}
    for role, source in inputs.items():
        suffix = "tar.gz" if role == "archive" else "json"
        name = f"la-stack-{policy.tag}-{selection}-candidate-{role}.{suffix}"
        payload = source.read_bytes()
        freeze(transport, {name: payload})
        publish_release_asset(root, REPOSITORY, policy.tag, transport / name, target=target, checkout=root)
        release = revalidate_release_target(root, target, checkout=root)
        entry = next(item for item in release.assets if item[0] == name)
        assertions[role] = {"name": name, "id": entry[1], "sha256": hashlib.sha256(payload).hexdigest()}
    payload = (
        json.dumps(
            {
                "schema": 1,
                "target": json.loads(serialize_release_target(target)),
                "acquisition": (owner / "acquisition.sha256").read_text(encoding="ascii").strip(),
                "publication": (owner / "publication.sha256").read_text(encoding="ascii").strip(),
                "assets": assertions,
            },
            sort_keys=True,
        ).encode("utf-8")
        + b"\n"
    )
    parse_transfer(payload)
    freeze(final, {"transfer.json": payload})
    return final / "transfer.json"


def seal(root: Path, owner: Path, selection: str) -> Path:
    """Independently select the canonical attached asset before publication."""
    final = final_directory(owner, resolve_revision(root, "HEAD"), selection)
    payload = (final / "transfer.json").read_bytes()
    transfer = parse_transfer(payload)
    release = revalidate_release_target(root, transfer.target, checkout=root)
    name = f"la-stack-{transfer.target.tag}-criterion-baseline.tar.gz"
    entry = next((item for item in release.assets if item[0] == name), None)
    if entry is None:
        message = "canonical accepted benchmark archive is not yet attached"
        raise ValueError(message)
    digest = next(item[3] for item in transfer.assets if item[0] == "archive")
    with tempfile.TemporaryDirectory(prefix="la-stack-seal-") as temporary:
        download_release_asset(root, REPOSITORY, transfer.target.tag, name, Path(temporary) / "archive.tar.gz", expected_sha256=digest, allow_draft=True)
    if entry not in revalidate_release_target(root, transfer.target, checkout=root).assets:
        message = "attached archive identity changed during verification"
        raise ValueError(message)
    data = json.loads(payload)
    data["assets"]["attached"] = {"name": name, "id": entry[1], "sha256": digest}
    sealed = json.dumps(data, sort_keys=True).encode("utf-8") + b"\n"
    parse_transfer(sealed)
    freeze(final, {"publication-transfer.json": sealed})
    return final / "publication-transfer.json"


def require_dispatch(transfer: Transfer, *, published: bool) -> None:
    """Bind independent workflow input to the allowed operator, workflow and tag."""
    workflow = "publish-crates.yml" if published else "release-benchmarks.yml"
    expected = {
        "GITHUB_EVENT_NAME": "workflow_dispatch",
        "GITHUB_REPOSITORY": REPOSITORY,
        "GITHUB_REF": f"refs/tags/{transfer.target.tag}",
        "GITHUB_SHA": transfer.target.commit,
        "GITHUB_WORKFLOW_REF": f"{REPOSITORY}/.github/workflows/{workflow}@refs/tags/{transfer.target.tag}",
        "GITHUB_ACTOR": "acgetchell",
        "GITHUB_TRIGGERING_ACTOR": "acgetchell",
        "RELEASE_TAG": transfer.target.tag,
    }
    if any(os.environ.get(key) != value for key, value in expected.items()):
        message = "candidate handoff requires the maintainer's exact workflow/tag dispatch"
        raise ValueError(message)


def transfer_files(root: Path, transfer: Transfer, directory: Path, incoming: Path | None) -> dict[str, Path]:
    """Read exact dispatch-selected bytes from published assets or inert transport."""
    files: dict[str, Path] = {}
    for role, name, _identifier, digest in transfer.assets:
        path = directory / name
        if incoming is None:
            download_release_asset(root, REPOSITORY, transfer.target.tag, name, path, expected_sha256=digest)
        else:
            source = incoming / ("archive.tar.gz" if role == "archive" else f"{role}.json")
            data = source.read_bytes()
            if source.is_symlink() or hashlib.sha256(data).hexdigest() != digest:
                message = "transported candidate bytes differ from independent dispatch assertions"
                raise ValueError(message)
            freeze(directory, {name: data})
        files[role] = path
    return files


def check_published_transfer(root: Path, transfer: Transfer) -> None:
    """Require the immutable release and exact original and canonical assets."""
    release = check_reviewed_release(
        load_config(root=root),
        transfer.target.tag,
        event=ReviewedRelease(transfer.target.tag, transfer.target.commit, transfer.target.release_id),
    )
    if not release.immutable:
        message = "published candidate release must be immutable"
        raise ValueError(message)
    if (release.repository, release.tag, release.identifier, release.title) != (
        REPOSITORY,
        transfer.target.tag,
        transfer.target.release_id,
        transfer.target.expected_title,
    ):
        message = "candidate release identity changed"
        raise ValueError(message)
    for _role, name, identifier, digest in transfer.assets:
        entry = next((item for item in release.assets if item[0] == name), None)
        if entry is None or entry[1] != identifier or (entry[3] is not None and entry[3] != digest):
            message = "candidate asset identity changed"
            raise ValueError(message)


@dataclass(frozen=True, slots=True)
class ReviewedAcceptance:
    """Independent owner digest, optionally paired with the frozen approval bytes."""

    sha256: str
    expected: bytes | None = None

    def __post_init__(self) -> None:
        """Reject missing review identity or a different approval receipt."""
        if re.fullmatch(r"[0-9a-f]{64}", self.sha256) is None:
            message = "publication requires the independently reviewed full acceptance SHA-256"
            raise ValueError(message)
        if self.expected is not None and hashlib.sha256(self.expected).hexdigest() != self.sha256:
            message = "approval receipt differs from independently reviewed acceptance"
            raise ValueError(message)


def receive(root: Path, destination: Path, payload: bytes, *, source: Path | ReviewedAcceptance) -> Path:
    """Read-only hosted acceptance, pinned to original bytes and reviewed CI."""
    published = isinstance(source, ReviewedAcceptance)
    transfer = parse_transfer(payload)
    require_dispatch(transfer, published=published)
    if published and "attached" not in {item[0] for item in transfer.assets}:
        message = "publication requires the independently selected canonical archive"
        raise ValueError(message)
    check_release_dependencies(root)
    destination = owner_directory(root, "", destination)
    destination.mkdir(parents=True, exist_ok=True)

    def release_identity() -> None:
        if published:
            check_published_transfer(root, transfer)

    release_identity()
    with tempfile.TemporaryDirectory(prefix="la-stack-transfer-") as temporary:
        files = transfer_files(root, transfer, Path(temporary), source if isinstance(source, Path) else None)
        release_identity()
        policy = parse_candidate_policy(files["policy"].read_bytes())
        accepted = files["acceptance"].read_bytes()
        evidence = verify_candidate_archive(files["archive"], policy, accepted, target=transfer.target)
        actual = verify_candidate_release(
            root,
            evidence,
            policy,
            load_evidence_policy(root / "tooling/release-evidence.toml"),
            partial(api, repository=REPOSITORY),
            acquisition_sha256=transfer.acquisition,
            publication_receipt=files["publication"].read_bytes(),
            publication_sha256=transfer.publication,
            # A published release may need explicitly reviewed newer native CI.
            # Before approval, the independent owner digest authenticates the
            # newly verified result. After approval, pin those exact bytes.
            expected=source.expected if isinstance(source, ReviewedAcceptance) else accepted,
        )
        if isinstance(source, ReviewedAcceptance) and hashlib.sha256(actual).hexdigest() != source.sha256:
            message = "verified acceptance differs from the independently reviewed digest"
            raise ValueError(message)
        release_identity()
        name = f"la-stack-{transfer.target.tag}-criterion-baseline.tar.gz"
        freeze(
            destination,
            {
                "transfer.json": payload,
                "acceptance.json": actual,
                "target.json": serialize_release_target(transfer.target),
                name: files["archive"].read_bytes(),
            },
        )
    return destination / name


def dispatch_source(options: argparse.Namespace) -> Path | ReviewedAcceptance:
    """Parse mutually exclusive draft transport and published approval inputs."""
    if options.published:
        if options.incoming is not None:
            message = "published verification requires live immutable assets"
            raise ValueError(message)
        return ReviewedAcceptance(os.environ.get("CANDIDATE_ACCEPTANCE", ""), options.expected.read_bytes() if options.expected is not None else None)
    if options.incoming is None or options.expected is not None:
        message = "draft verification requires inert files and no publication approval override"
        raise ValueError(message)
    return Path(options.incoming)


def main() -> None:
    """Expose only repository orchestration; shared APIs own every evidence check."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="operation", required=True)
    measurement = commands.add_parser("acquire")
    measurement.add_argument("--checkpoint", required=True)
    measurement.add_argument("--resume", action="store_true")
    measurement.add_argument("--owner-directory", type=Path)
    measurement.add_argument("--progress", choices=("human", "json", "none"), default="human")
    measurement.add_argument("--progress-cadence", type=float, default=300)
    for name in ("capture", "verify", "archive", "stage", "seal"):
        command = commands.add_parser(name)
        command.add_argument("owner_directory", type=Path)
        if name == "capture":
            command.add_argument("--retained-scratch", action="store_true", help="explicitly retry retention of the original locally controlled acquisition")
        else:
            command.add_argument("--selection", default="reviewed", help="immutable review record; choose a new name to review a different native CI attempt")
    receiver = commands.add_parser("receive")
    receiver.add_argument("destination", type=Path)
    receiver.add_argument("--published", action="store_true")
    receiver.add_argument("--incoming", type=Path)
    receiver.add_argument("--expected", type=Path)
    options = parser.parse_args()
    root = Path.cwd().resolve()
    if options.operation == "acquire":
        owner = owner_directory(root, options.checkpoint, options.owner_directory)
        print(f"Independent acquisition record: {owner}", flush=True)
        result = acquire(
            root,
            owner,
            options.checkpoint,
            resume=options.resume,
            progress=ProgressReporter(mode=options.progress, cadence=options.progress_cadence),
        )
        print(json.dumps(result))
    elif options.operation == "receive":
        print(receive(root, options.destination, os.environ["CANDIDATE_TRANSFER"].encode("utf-8"), source=dispatch_source(options)))
    else:
        owner = owner_directory(root, "", options.owner_directory)
        if options.operation == "capture":
            print(json.dumps(capture(root, owner, retained_scratch=options.retained_scratch)))
        elif options.operation == "verify":
            accepted = verify(root, owner, options.selection)
            print(json.dumps({"acceptance": str(accepted), "sha256": hashlib.sha256(accepted.read_bytes()).hexdigest()}))
        elif options.operation == "archive":
            print(archive(owner, resolve_revision(root, "HEAD"), options.selection))
        elif options.operation == "stage":
            print(stage(root, owner, options.selection))
        else:
            print(seal(root, owner, options.selection))


if __name__ == "__main__":
    main()
