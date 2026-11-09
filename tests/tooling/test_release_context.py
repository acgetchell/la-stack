"""Resolve original release ownership and receipt selections without shell state."""

import base64
import hashlib
import json
import shlex
import subprocess
import sys
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
from research_repo_tools import release_interfaces, release_receipts
from research_repo_tools.process import run_command
from research_repo_tools.release_assets import ReleaseTarget, serialize_release_target
from research_repo_tools.release_drafts import RECEIPT_SCHEMA
from research_repo_tools.release_interfaces import load_receipt_policy

from tooling import release_candidate as candidate

ROOT = Path(__file__).resolve().parents[2]
TAG = "v0.5.42"
SHA = "a" * 40


@pytest.fixture
def owner_record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    root, owner = tmp_path / "checkout", tmp_path / "original owner"
    root.mkdir()
    owner.mkdir()
    (root / "Cargo.toml").write_text('[package]\nname = "la-stack"\nversion = "0.5.42"\n', encoding="utf-8", newline="\n")
    (owner / "policy.json").write_bytes(b"original policy")
    (owner / "publication.json").write_bytes(b"original publisher receipt\r\n")
    (owner / "acquisition.sha256").write_text("b" * 64, encoding="utf-8")
    (owner / "publication.sha256").write_text(hashlib.sha256((owner / "publication.json").read_bytes()).hexdigest(), encoding="utf-8")
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    # The complete acquisition fixture separately exercises the real policy parser/digest.
    monkeypatch.setattr(candidate, "parse_candidate_policy", lambda _payload: SimpleNamespace(tag=TAG))
    monkeypatch.setattr(candidate, "acquisition_assertion_digest", lambda _policy: "b" * 64)
    return root, owner


def test_owner_discovery_requires_the_original_step_one_registration(owner_record: tuple[Path, Path]) -> None:
    root, owner = owner_record
    with pytest.raises(ValueError, match="no original Step 1 owner selection"):
        candidate.selected_owner(root)
    assert not candidate.selection_directory(root).exists()
    candidate.remember_owner(root, owner)
    assert candidate.selected_owner(root) == owner.resolve()
    retained = (candidate.selection_directory(root) / "owner.json").read_bytes()
    other = owner.parent / "other owner"
    other.mkdir()
    for source in owner.iterdir():
        (other / source.name).write_bytes(source.read_bytes())
    assert candidate.selected_owner(root, str(other)) == other.resolve()
    with pytest.raises(ValueError, match="immutable output differs"):
        candidate.remember_owner(root, other)
    assert (candidate.selection_directory(root) / "owner.json").read_bytes() == retained


@pytest.mark.parametrize("failure", ["tag", "checkout", "owner", "schema", "digest", "publication"])
def test_owner_discovery_rejects_stale_or_changed_records(owner_record: tuple[Path, Path], failure: str) -> None:
    root, owner = owner_record
    candidate.remember_owner(root, owner)
    path = candidate.selection_directory(root) / "owner.json"
    data = json.loads(path.read_bytes())
    if failure == "digest":
        (owner / "acquisition.sha256").write_text("c" * 64, encoding="utf-8")
    elif failure == "publication":
        (owner / "publication.json").write_bytes(b"changed receipt")
    else:
        data[failure] = {"tag": "v9.9.9", "checkout": "other", "owner": str(root), "schema": True}[failure]
        path.write_text(json.dumps(data), encoding="utf-8", newline="\n")
    with pytest.raises(ValueError, match=r"selection differs|separate from the checkout|assertions changed|publication receipt differs"):
        candidate.selected_owner(root)


@pytest.fixture
def receipt_provider(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """Keep the installed tag-aware verifier and immutable CLI output transaction."""
    root, owner = tmp_path / "checkout", tmp_path / "owner"
    root.mkdir()
    owner.mkdir()
    (root / "Cargo.toml").write_text('[package]\nname = "la-stack"\nversion = "0.5.42"\n', encoding="utf-8", newline="\n")
    (root / "tooling").mkdir()
    (root / "tooling/draft-receipt.toml").write_bytes((ROOT / "tooling/draft-receipt.toml").read_bytes())
    policy = load_receipt_policy(root / "tooling/draft-receipt.toml")
    target = ReleaseTarget(candidate.REPOSITORY, TAG, 42, SHA, TAG)
    target_bytes = serialize_release_target(target)
    receipt = json.dumps(
        {
            "schema": RECEIPT_SCHEMA,
            "repository": candidate.REPOSITORY,
            "tag": TAG,
            "commit": SHA,
            "title": TAG,
            "notes_sha256": "b" * 64,
            "selection_sha256": "c" * 64,
            "state": "captured",
            "release_id": 42,
        }
    ).encode()
    repository = {"id": 11, "full_name": candidate.REPOSITORY}
    run = {
        "id": 77,
        "run_attempt": 1,
        "workflow_id": 9,
        "path": policy.workflow,
        "name": policy.workflow_name,
        "event": "push",
        "head_branch": TAG,
        "head_sha": SHA,
        "pull_requests": [],
        "repository": repository,
        "head_repository": repository,
        "status": "completed",
        "conclusion": "success",
    }
    artifact = {
        "id": 88,
        "name": "draft-receipt-77-1",
        "digest": "sha256:" + "d" * 64,
        "expired": False,
        "size_in_bytes": 512,
        "workflow_run": {"id": 77, "head_branch": TAG, "head_sha": SHA, "repository_id": 11, "head_repository_id": 11},
    }
    jobs = [
        {
            "id": index,
            "name": name,
            "run_id": 77,
            "run_attempt": 1,
            "head_branch": TAG,
            "head_sha": SHA,
            "status": "completed",
            "conclusion": "success",
            "labels": ["ubuntu-latest"],
            "steps": [{"name": step, "number": number, "status": "completed", "conclusion": "success"} for number, step in enumerate(steps, 1)],
        }
        for index, (name, steps) in enumerate(policy.required_steps, 91)
    ]
    workflow_bytes = (ROOT / policy.workflow).read_bytes()
    source = {
        "type": "file",
        "path": policy.workflow,
        "encoding": "base64",
        "size": len(workflow_bytes),
        "sha": hashlib.sha1(b"blob " + str(len(workflow_bytes)).encode() + b"\0" + workflow_bytes, usedforsecurity=False).hexdigest(),
        "content": base64.b64encode(workflow_bytes).decode(),
    }
    responses: dict[str, Any] = {
        "actions/workflows/prepare-release.yml": {"id": 9, "path": policy.workflow, "name": policy.workflow_name, "state": "active"},
        "actions/workflows/9/runs?per_page=100&page=1": {"workflow_runs": [run]},
        "actions/runs/77": run,
        "actions/runs/77/attempts/1/jobs?per_page=100&page=1": {"jobs": jobs},
        "actions/runs/77/artifacts?per_page=100&page=1": {"artifacts": [artifact]},
        "actions/artifacts/88": artifact,
        f"contents/{policy.workflow}?ref={SHA}": source,
        f"git/ref/tags/{TAG}": {"ref": f"refs/tags/{TAG}", "object": {"type": "commit", "sha": SHA}},
    }
    calls: list[str] = []

    def api(endpoint: str, *, repository: str) -> Any:
        assert repository == candidate.REPOSITORY
        calls.append(endpoint)
        assert endpoint in responses, endpoint
        return deepcopy(responses[endpoint])

    def download(repository: str, observed: dict[str, Any], destination: Path, names: tuple[str, ...], **_bounds: int) -> None:
        assert repository == candidate.REPOSITORY
        assert observed["id"] == 88
        assert names == ("draft-receipt.json", "target.json")
        destination.mkdir()
        (destination / "target.json").write_bytes(target_bytes)
        (destination / "draft-receipt.json").write_bytes(receipt)

    monkeypatch.setattr(release_interfaces, "api", api)
    monkeypatch.setattr(release_receipts, "download_artifact", download)
    monkeypatch.setattr(candidate, "resolve_revision", lambda *_args: SHA)
    monkeypatch.setattr(candidate, "require_release_checkout", lambda *_args: None)
    return {"root": root, "owner": owner, "run": run, "artifact": artifact, "jobs": jobs, "calls": calls, "responses": responses}


def test_automatic_receipt_retrieval_reuses_original_frozen_identity(receipt_provider: dict[str, Any]) -> None:
    fixture = receipt_provider
    target = candidate.preparation_target(fixture["root"], fixture["owner"], "reviewed")
    original = {path.name: path.read_bytes() for path in target.parent.iterdir()}
    assert json.loads(target.read_bytes())["release_id"] == 42
    fixture["calls"].clear()
    assert candidate.preparation_target(fixture["root"], fixture["owner"], "reviewed") == target
    assert {path.name: path.read_bytes() for path in target.parent.iterdir()} == original
    assert not any("actions/workflows/9/runs" in endpoint for endpoint in fixture["calls"])
    fixture["artifact"]["digest"] = "sha256:" + "e" * 64
    with pytest.raises(RuntimeError, match="retrieval failed"):
        candidate.preparation_target(fixture["root"], fixture["owner"], "reviewed")
    assert {path.name: path.read_bytes() for path in target.parent.iterdir()} == original


@pytest.mark.parametrize("failure", ["attempt", "step", "missing", "ambiguous"])
def test_automatic_receipt_selection_stops_without_output_on_invalid_evidence(receipt_provider: dict[str, Any], failure: str) -> None:
    fixture = receipt_provider
    if failure == "attempt":
        fixture["run"]["run_attempt"] = 2
    elif failure == "step":
        fixture["jobs"][1]["steps"][-1]["conclusion"] = "skipped"
    else:
        history = fixture["responses"]["actions/workflows/9/runs?per_page=100&page=1"]["workflow_runs"]
        history[:] = [] if failure == "missing" else [*history, {**fixture["run"], "id": 78}]
    with pytest.raises(RuntimeError, match="retrieval failed"):
        candidate.preparation_target(fixture["root"], fixture["owner"], "reviewed")
    assert not (fixture["owner"] / "final" / SHA / "reviewed/preparation").exists()


@pytest.mark.parametrize("failure", ["none", "missing-review", "changed-review", "mutable-release", "dispatch"])
def test_registry_dispatch_sends_only_the_existing_reviewed_snapshot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str) -> None:
    root, owner = tmp_path / "checkout", tmp_path / "owner"
    root.mkdir()
    owner.mkdir()
    final = candidate.final_directory(owner, SHA, "reviewed")
    final.mkdir(parents=True)
    publication = final / "publication-transfer.json"
    target = ReleaseTarget(candidate.REPOSITORY, TAG, 42, SHA, TAG)
    assets = {role: {"name": f"{role}.bin", "id": index, "sha256": "b" * 64} for index, role in enumerate((*candidate.ROLES, "attached"), 1)}
    assets["attached"]["name"] = f"la-stack-{TAG}-criterion-baseline.tar.gz"
    original = (
        json.dumps(
            {"schema": 1, "target": json.loads(serialize_release_target(target)), "acquisition": "c" * 64, "publication": "d" * 64, "assets": assets}
        ).encode()
        + b"\r\n"
    )
    publication.write_bytes(original)
    accepted = final / "acceptance.json"
    reviewed = b"original full acceptance\r\n"
    if failure != "missing-review":
        accepted.write_bytes(reviewed)
    monkeypatch.setattr(candidate, "release_tag", lambda _root: TAG)
    monkeypatch.setattr(candidate, "record_path", lambda *_args: publication)
    verified: list[Path] = []

    def verify(*_args: object) -> Path:
        verified.append(accepted)
        assert accepted.read_bytes() == reviewed
        if failure == "changed-review":
            message = "original reviewed native attempt changed"
            raise ValueError(message)
        return accepted

    def immutable(*_args: object) -> None:
        if failure == "mutable-release":
            message = "published release is not immutable"
            raise ValueError(message)

    calls: list[list[str]] = []

    def dispatch(executable: str, args: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
        assert executable == "gh"
        assert cwd == root
        assert args[:9] == ["workflow", "run", "publish-crates.yml", "--repo", candidate.REPOSITORY, "--ref", TAG, "-f", f"tag={TAG}"]
        assert args[-2:] == ["-f", "acceptance=" + hashlib.sha256(reviewed).hexdigest()]
        assert Path(args[10].removeprefix("candidate=@")).read_bytes() == original
        calls.append(args)
        if failure == "dispatch":
            raise subprocess.CalledProcessError(17, [executable, *args])
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(candidate, "verify_existing", verify)
    monkeypatch.setattr(candidate, "check_published_transfer", immutable)
    monkeypatch.setattr(candidate, "run_command", dispatch)
    assert_dispatch_outcome(root, owner, verified, failure)
    assert len(calls) == (1 if failure in {"none", "dispatch"} else 0)
    assert publication.read_bytes() == original


def assert_dispatch_outcome(root: Path, owner: Path, verified: list[Path], failure: str) -> None:
    """Require each injected boundary failure to retain its intended diagnostic."""
    if failure == "none":
        candidate.upload_dispatch(root, owner, "reviewed", TAG)
    elif failure == "dispatch":
        with pytest.raises(subprocess.CalledProcessError) as stopped:
            candidate.upload_dispatch(root, owner, "reviewed", TAG)
        assert stopped.value.returncode == 17
    elif failure == "missing-review":
        with pytest.raises(FileNotFoundError, match=r"acceptance\.json"):
            candidate.upload_dispatch(root, owner, "reviewed", TAG)
        assert verified == []
    else:
        with pytest.raises(ValueError, match=r"native attempt changed|not immutable"):
            candidate.upload_dispatch(root, owner, "reviewed", TAG)


def test_downstream_checks_never_create_missing_review(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root, owner = tmp_path / "checkout", tmp_path / "owner"
    root.mkdir()
    owner.mkdir()
    monkeypatch.setattr(candidate, "resolve_revision", lambda *_args: SHA)
    monkeypatch.setattr(candidate, "verify", lambda *_args: pytest.fail("tagging or publication must not create a new review"))
    with pytest.raises(FileNotFoundError, match=r"acceptance\.json"):
        candidate.verify_existing(root, owner, "reviewed")
    assert list(owner.iterdir()) == []


@pytest.mark.parametrize("owner_argument", [[], [""]])
def test_cli_recovers_the_registered_owner_without_shell_context(
    owner_record: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], owner_argument: list[str]
) -> None:
    root, owner = owner_record
    candidate.remember_owner(root, owner)
    accepted = owner / "acceptance.json"
    accepted.write_bytes(b"reviewed acceptance")
    monkeypatch.chdir(root)
    monkeypatch.setattr("sys.argv", ["release_candidate.py", "verify", *owner_argument])
    received: list[Path] = []

    def verify(_root: Path, selected: Path, selection: str) -> Path:
        assert _root == root
        assert selection == "reviewed"
        received.append(selected)
        return accepted

    monkeypatch.setattr(candidate, "verify", verify)
    candidate.main()
    assert received == [owner.resolve()]
    assert json.loads(capsys.readouterr().out)["sha256"] == hashlib.sha256(accepted.read_bytes()).hexdigest()


def test_conflicting_owner_selection_stops_before_acquisition(owner_record: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch) -> None:
    root, owner = owner_record
    candidate.remember_owner(root, owner)
    other = owner.parent / "new acquisition"
    monkeypatch.setattr(candidate, "load_performance_release", lambda *_args: None)
    monkeypatch.setattr(candidate, "load_preparation", lambda *_args: None)
    monkeypatch.setattr(candidate, "check_release_dependencies", lambda *_args: None)
    monkeypatch.setattr(candidate, "measure_release", lambda *_args, **_kwargs: pytest.fail("a known owner conflict must not start acquisition"))
    with pytest.raises(ValueError, match="already selects a different original owner"):
        candidate.acquire(root, other, "target/new-checkpoint", resume=False, progress=candidate.ProgressReporter(mode="none"))
    assert not other.exists()


def test_publish_recipe_requires_review_before_reaching_publication(tmp_path: Path) -> None:
    root, owner = tmp_path / "checkout", tmp_path / "owner"
    root.mkdir()
    owner.mkdir()
    source = (ROOT / "justfile").read_text(encoding="utf-8")
    header = next(line for line in source.splitlines() if line.startswith("release-publish "))
    bootstrap = root / "candidate-bootstrap.py"
    bootstrap.write_text(
        "import sys\nfrom pathlib import Path\n"
        f"sys.path.insert(0, {str(ROOT)!r})\n"
        "from tooling import release_candidate as candidate\n"
        f"candidate.resolve_revision = lambda *_args: {SHA!r}\n"
        'def unexpected(*_args):\n    raise RuntimeError("attempted to create a new review")\n'
        "candidate.verify = unexpected\ncandidate.main()\n",
        encoding="utf-8",
        newline="\n",
    )
    marker = root / "published"
    overrides = (
        "\n_release-notes-check tag:\n    @true\n\n[positional-arguments]\n_release-candidate *args:\n    "
        f'{shlex.quote(Path(sys.executable).as_posix())} -I {shlex.quote(bootstrap.as_posix())} "$@"\n\n'
        + header
        + f"\n    touch {shlex.quote(marker.as_posix())}\n"
    )
    (root / "justfile").write_text("set allow-duplicate-recipes\n" + source + overrides, encoding="utf-8", newline="\n")
    result = run_command("just", ["release-publish", TAG, str(owner)], cwd=root, check=False)
    assert result.returncode != 0
    assert "acceptance.json" in result.stderr
    assert "attempted to create a new review" not in result.stderr
    assert not marker.exists()
