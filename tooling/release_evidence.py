"""Bind Cargo publication to signed, successful main-push CI evidence.

The installed research-repo-tools owns run/job validation and safe artifact
retrieval. This adapter supplies la-stack's Cargo policy and evidence document;
it never uploads a crate, reruns checks, or accepts substitute evidence.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path
from typing import Any

from research_repo_tools.ci import export_environment
from research_repo_tools.process import format_exception_diagnostics, run_git_command
from research_repo_tools.release_validation import ValidationPolicy, api, download_artifact, gh, select_main_run, validate_jobs, validate_run

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "acgetchell/la-stack"
POLICY = ValidationPolicy(
    REPOSITORY,
    ".github/workflows/ci.yml",
    "CI",
    (
        *((f"build ({platform})", ("Run CI checks",)) for platform in ("ubuntu-latest", "macos-latest", "windows-latest")),
        ("release-package", ("Verify default and exact Cargo packages",)),
    ),
)
EVIDENCE_STEPS = ("Capture successful native validation identities", "Sign original main validation evidence", "Upload signed validation evidence")


def read(endpoint: str) -> Any:
    """Read authenticated JSON in this repository only."""
    return api(endpoint, repository=REPOSITORY)


def object_value(value: Any) -> dict[str, Any]:
    """Require a JSON object at a provider boundary."""
    if not isinstance(value, dict):
        msg = "expected an evidence object"
        raise TypeError(msg)
    return value


def positive_integer(value: Any) -> int:
    """Reject missing, boolean, or nonpositive provider identifiers."""
    if type(value) is not int or value <= 0:
        msg = "expected a positive integer identity"
        raise ValueError(msg)
    return value


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Reject duplicate keys instead of silently replacing evidence fields."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate evidence key: {key}")
        result[key] = value
    return result


def pages(endpoint: str, field: str) -> list[dict[str, Any]]:
    """Read a bounded, fully paginated provider collection."""
    result: list[dict[str, Any]] = []
    for page in range(1, 101):
        data = object_value(read(f"{endpoint}?per_page=100&page={page}"))
        entries = data.get(field)
        if not isinstance(entries, list):
            raise TypeError(f"invalid {field} collection")
        result.extend(object_value(entry) for entry in entries)
        if len(entries) < 100:
            return result
    raise ValueError(f"{field} exceeds the bounded page limit")


def checkout_commit() -> str:
    """Resolve the checkout independently of event or artifact claims."""
    return run_git_command(["--no-pager", "rev-parse", "HEAD"], cwd=ROOT, timeout=30).stdout.strip()


def current_run(commit: str, identifier: int, attempt: int, *, completed: bool = True) -> dict[str, Any]:
    """Reject reruns which replaced the selected current attempt."""
    run = object_value(read(f"actions/runs/{positive_integer(identifier)}"))
    validate_run(POLICY, run, commit, object_value(read("actions/workflows/ci.yml")), completed=completed)
    if run["run_attempt"] != positive_integer(attempt):
        msg = "selected CI attempt is stale; inspect the complete main-push run"
        raise ValueError(msg)
    return run


def manifest(run: dict[str, Any]) -> dict[str, Any]:
    """Capture native check and default/exact package job identities."""
    jobs = pages(f"actions/runs/{run['id']}/attempts/{run['run_attempt']}/jobs", "jobs")
    package = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))["package"]
    return {
        "schema": 1,
        "repository": REPOSITORY,
        "repository_id": run["repository"]["id"],
        "commit": run["head_sha"],
        "run_id": run["id"],
        "run_attempt": run["run_attempt"],
        "jobs": validate_jobs(POLICY, run, jobs),
        "package": {"name": package["name"], "version": package["version"]},
    }


def seal(output: Path) -> None:
    """Produce the document signed by this main-push CI attempt."""
    commit = checkout_commit()
    if os.environ.get("GITHUB_SHA") != commit:
        msg = "executing CI source differs from checkout"
        raise ValueError(msg)
    run = current_run(commit, int(os.environ["GITHUB_RUN_ID"]), int(os.environ["GITHUB_RUN_ATTEMPT"]), completed=False)
    document = manifest(run)
    current_run(commit, run["id"], run["run_attempt"], completed=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(document, sort_keys=True) + "\n")


def artifact_identity(data: dict[str, Any], run: dict[str, Any]) -> dict[str, Any]:
    """Validate immutable artifact association before shared ZIP retrieval."""
    name = f"release-evidence-{run['id']}-{run['run_attempt']}"
    association = object_value(data.get("workflow_run"))
    expected = {
        "id": run["id"],
        "head_sha": run["head_sha"],
        "head_branch": "main",
        "repository_id": run["repository"]["id"],
        "head_repository_id": run["head_repository"]["id"],
    }
    if data.get("name") != name or data.get("expired") is not False or any(association.get(key) != value for key, value in expected.items()):
        msg = "evidence artifact is expired or belongs to another run, attempt, repository, ref or commit"
        raise ValueError(msg)
    digest = data.get("digest")
    if not isinstance(digest, str) or re.fullmatch(r"sha256:[0-9a-f]{64}", digest) is None:
        msg = "evidence artifact lacks an authoritative SHA-256 digest"
        raise ValueError(msg)
    size = positive_integer(data.get("size_in_bytes"))
    if size > 1024 * 1024:
        msg = "evidence artifact exceeds the one-MiB limit"
        raise ValueError(msg)
    return {"id": positive_integer(data.get("id")), "name": name, "digest": digest, "size_in_bytes": size}


def find_artifact(run: dict[str, Any]) -> dict[str, Any]:
    """Require one evidence artifact for the selected attempt."""
    name = f"release-evidence-{run['id']}-{run['run_attempt']}"
    matches = [entry for entry in pages(f"actions/runs/{run['id']}/artifacts", "artifacts") if entry.get("name") == name]
    if len(matches) != 1:
        msg = "missing or ambiguous signed CI evidence; inspect retention and rerun all CI jobs if needed"
        raise ValueError(msg)
    return artifact_identity(matches[0], run)


def verify_signature(directory: Path, run: dict[str, Any]) -> None:
    """Verify Sigstore provenance, including the exact signing attempt."""
    result = json.loads(
        gh(
            "attestation",
            "verify",
            str(directory / "validation.json"),
            "--bundle",
            str(directory / "validation.sigstore.json"),
            "--repo",
            REPOSITORY,
            "--signer-workflow",
            f"{REPOSITORY}/{POLICY.workflow}",
            "--signer-digest",
            run["head_sha"],
            "--source-digest",
            run["head_sha"],
            "--source-ref",
            "refs/heads/main",
            "--deny-self-hosted-runners",
            "--predicate-type",
            "https://slsa.dev/provenance/v1",
            "--format",
            "json",
        ),
        object_pairs_hook=unique_object,
    )
    if not isinstance(result, list) or len(result) != 1:
        msg = "expected exactly one verified CI attestation"
        raise ValueError(msg)
    verification = object_value(object_value(result[0]).get("verificationResult"))
    certificate = object_value(object_value(verification.get("signature")).get("certificate"))
    invocation = f"https://github.com/{REPOSITORY}/actions/runs/{run['id']}/attempts/{run['run_attempt']}"
    if certificate.get("runInvocationURI") != invocation or certificate.get("buildTrigger") != "push":
        msg = "attestation run, attempt or push event differs"
        raise ValueError(msg)


def verify(identity: str = "") -> dict[str, Any]:
    """Reuse one exact CI attempt, or recheck the same pinned identity after approval."""
    commit = checkout_commit()
    pin = object_value(json.loads(identity, object_pairs_hook=unique_object)) if identity else None
    if pin is not None and set(pin) != {"run_id", "run_attempt", "artifact_id", "artifact_digest"}:
        msg = "invalid pinned evidence identity"
        raise ValueError(msg)
    run = current_run(commit, pin["run_id"], pin["run_attempt"]) if pin is not None else select_main_run(POLICY, commit, read)
    expected = manifest(run)
    signing_policy = ValidationPolicy(REPOSITORY, POLICY.workflow, POLICY.workflow_name, (("release-evidence", EVIDENCE_STEPS),))
    validate_jobs(signing_policy, run, pages(f"actions/runs/{run['id']}/attempts/{run['run_attempt']}/jobs", "jobs"))
    artifact = (
        artifact_identity(object_value(read(f"actions/artifacts/{positive_integer(pin['artifact_id'])}")), run) if pin is not None else find_artifact(run)
    )
    result = {"run_id": run["id"], "run_attempt": run["run_attempt"], "artifact_id": artifact["id"], "artifact_digest": artifact["digest"]}
    if pin is not None and pin != result:
        msg = "approved evidence identity changed"
        raise ValueError(msg)
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary) / "evidence"
        download_artifact(REPOSITORY, artifact, directory, ("validation.json", "validation.sigstore.json"))
        verify_signature(directory, run)
        document = json.loads((directory / "validation.json").read_text(encoding="utf-8"), object_pairs_hook=unique_object)
        if document != expected:
            msg = "signed evidence differs from current exact-commit Cargo and native validation identities"
            raise ValueError(msg)
    current_run(commit, run["id"], run["run_attempt"])
    if artifact_identity(object_value(read(f"actions/artifacts/{artifact['id']}")), run) != artifact:
        msg = "evidence artifact changed during verification"
        raise ValueError(msg)
    return result


def main() -> int:
    """Expose signing and read-only consumption without publication credentials."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("seal").add_argument("--output", type=Path, required=True)
    verification = commands.add_parser("verify")
    verification.add_argument("--identity", default="")
    verification.add_argument("--github-output", type=Path, help="append identity to an existing GitHub Actions output file instead of stdout")
    args = parser.parse_args()
    try:
        if args.command == "seal":
            seal(args.output)
        else:
            identity = json.dumps(verify(args.identity), sort_keys=True, separators=(",", ":"), allow_nan=False)
            if args.github_output is not None:
                if not args.github_output.is_file():
                    msg = "GitHub Actions output destination must be an existing regular file"
                    raise ValueError(msg)
                export_environment(args.github_output, ("identity",), environment={"identity": identity})
            else:
                print(identity)
    except (KeyError, OSError, TypeError, ValueError, subprocess.SubprocessError) as error:
        print(f"Release evidence rejected: {format_exception_diagnostics(error, single_line=True)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
