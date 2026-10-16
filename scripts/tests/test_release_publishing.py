"""Exercise la-stack's published release contract without credentials or uploads."""

import json
import os
import shlex
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml
from research_repo_tools import changelog, config, release_publishing
from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_command

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github/workflows/publish-crates.yml"
SHA = "a" * 40
POLICY = config.load(root=REPO_ROOT)
VERSION = tomllib.loads((REPO_ROOT / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
TAG = f"v{VERSION}"
CHECKS = ("build (ubuntu-latest)", "build (macos-latest)", "build (windows-latest)", "Clippy SARIF Analysis", "Repository Rule SARIF Analysis")


def test_publishing_policy_uses_the_published_pin_and_native_scientific_gates() -> None:
    assert POLICY.publishing is not None
    assert (POLICY.publishing.registry, POLICY.publishing.package, POLICY.publishing.repository) == ("crates-io", "la-stack", "acgetchell/la-stack")
    assert tuple((check.name, check.app_id) for check in POLICY.publishing.required_checks) == tuple((name, 15368) for name in CHECKS)
    assert POLICY.publishing.required_assets == ("la-stack-v{version}-criterion-baseline.tar.gz",)
    manifest = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest["dependency-groups"]["tooling"] == ["research-repo-tools==0.1.9"]
    lock = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    package = next(item for item in lock["package"] if item["name"] == "research-repo-tools")
    assert package["version"] == "0.1.9"
    assert package["source"] == {"registry": "https://pypi.org/simple"}


@pytest.fixture
def release_event(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """Supply only provider/read-only Git evidence to the installed release CLI."""
    event = {
        "action": "published",
        "repository": {"full_name": "acgetchell/la-stack"},
        "release": {"id": 42, "tag_name": TAG, "draft": False, "prerelease": False},
    }
    event_file = tmp_path / "event.json"
    event_file.write_text(json.dumps(event), encoding="utf-8")
    for key, value in {
        "GITHUB_EVENT_NAME": "release",
        "GITHUB_REPOSITORY": "acgetchell/la-stack",
        "GITHUB_SHA": SHA,
        "GITHUB_REF": f"refs/tags/{TAG}",
        "GITHUB_EVENT_PATH": str(event_file),
    }.items():
        monkeypatch.setenv(key, value)

    def git(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        assert args in (["--no-pager", "rev-parse", "HEAD"], ["--no-pager", "status", "--porcelain", "--untracked-files=all"])
        return subprocess.CompletedProcess(args, 0, stdout=SHA if "rev-parse" in args else "", stderr="")

    checks = [
        {"id": index, "name": name, "head_sha": SHA, "app": {"id": 15368}, "status": "completed", "conclusion": "success"}
        for index, name in enumerate(CHECKS, 1)
    ]
    responses: dict[str, object] = {
        "": {"full_name": "acgetchell/la-stack", "default_branch": "main"},
        "branches/main": {"name": "main", "protected": True},
        f"compare/{SHA}...main": {"status": "identical"},
        f"commits/{SHA}/check-runs?filter=latest&per_page=100": [{"total_count": len(checks), "check_runs": checks}],
    }
    monkeypatch.setattr(release_publishing, "run_git_command", git)
    monkeypatch.setattr(release_publishing, "_api", lambda _root, _repo, endpoint, **_kwargs: responses[endpoint])
    monkeypatch.setattr(release_publishing, "_remote_commit", lambda *_args: SHA)
    release = release_publishing.GitHubRelease(
        repository="acgetchell/la-stack",
        tag=TAG,
        identifier=42,
        draft=False,
        immutable=True,
        prerelease=False,
        assets=((f"la-stack-{TAG}-criterion-baseline.tar.gz", 7, 100, "digest"),),
    )
    monkeypatch.setattr(release_publishing, "lookup_release", lambda *_args: release)
    return {"event": event, "path": event_file, "checks": checks}


def test_installed_gate_accepts_actual_metadata_and_required_native_results(release_event: dict[str, Any]) -> None:
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 0


@pytest.mark.parametrize("name", CHECKS)
@pytest.mark.parametrize("state", ["missing", "failure", "pending"])
def test_every_native_gate_blocks_missing_or_unsuccessful_evidence(release_event: dict[str, Any], name: str, state: str) -> None:
    check = next(item for item in release_event["checks"] if item["name"] == name)
    if state == "missing":
        check["name"] = "unrelated check"
    else:
        check.update(status="in_progress" if state == "pending" else "completed", conclusion=None if state == "pending" else "failure")
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 1


@pytest.mark.parametrize(("field", "value"), [("draft", True), ("prerelease", True), ("tag_name", "v99.0.0")])
def test_installed_gate_rejects_mismatched_release_events(release_event: dict[str, Any], field: str, value: object) -> None:
    release_event["event"]["release"][field] = value
    release_event["path"].write_text(json.dumps(release_event["event"]), encoding="utf-8")
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 1


@pytest.mark.parametrize("mismatch", ["repository", "commit", "version", "event"])
def test_consumer_gate_rejects_wrong_identity(release_event: dict[str, Any], monkeypatch: pytest.MonkeyPatch, mismatch: str) -> None:
    tag = TAG
    if mismatch == "repository":
        monkeypatch.setenv("GITHUB_REPOSITORY", "other/repository")
    elif mismatch == "commit":
        monkeypatch.setattr(release_publishing, "_remote_commit", lambda *_args: "b" * 40)
    elif mismatch == "version":
        tag = "v99.0.0"
    else:
        monkeypatch.setenv("GITHUB_EVENT_NAME", "pull_request")
    assert main(["--root", str(REPO_ROOT), "release", "gate", tag]) == 1


@pytest.mark.parametrize("immutable", [True, False, None])
def test_recipe_requires_immutable_publication(tmp_path: Path, immutable: bool | None) -> None:
    payload = tmp_path / "event.json"
    payload.write_text(json.dumps({"release": {"id": 42}}), encoding="utf-8")
    # The shared gate has its own installed-CLI tests above. Isolate this local
    # immutable-release policy while running the actual Just recipe expansion.
    script = dry_run(REPO_ROOT, "release-gate", [TAG]).stderr
    provider = (
        f'gh() {{ [[ "$1 $2" == "api repos/acgetchell/la-stack/releases/42" ]] || return 90; echo {shlex.quote(json.dumps({"immutable": immutable}))}; }}\n'
    )
    result = run_command(
        "bash",
        ["--noprofile", "--norc", "-euc", "set -o pipefail\nuv() { return 0; }\njust() { return 0; }\n" + provider + script],
        env={
            **{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}},
            "GITHUB_EVENT_PATH": payload.as_posix(),
        },
        check=False,
    )
    assert (result.returncode == 0) == (immutable is True)


@pytest.mark.parametrize("mismatch", ["none", "title", "notes"])
def test_real_notes_recipe_checks_generated_body_and_title(tmp_path: Path, mismatch: str) -> None:
    body = changelog.notes(POLICY, TAG)[0]
    response = {"tagName": TAG, "name": "wrong" if mismatch == "title" else TAG, "body": "wrong" if mismatch == "notes" else body}
    payload = tmp_path / "release.json"
    payload.write_text(json.dumps(response), encoding="utf-8")
    # Only acquisition/provider boundaries are substituted; execute the real CLI
    # and recipe's jq comparison, preserving safe transport of Markdown bytes.
    runner = (
        'uv() { [[ "$1 $2 $3 $4" == "run --locked --no-sync research-repo-tools" ]] || return 90; shift 4; '
        f'{shlex.quote(sys.executable)} -m research_repo_tools --root {shlex.quote(str(REPO_ROOT))} "$@"; }}\n'
        f"gh() {{ cat {shlex.quote(payload.as_posix())}; }}\n"
    )
    script = dry_run(REPO_ROOT, "_release-notes-check", [TAG]).stderr
    result = run_command("bash", ["--noprofile", "--norc", "-euc", runner + script], cwd=REPO_ROOT, check=False)
    assert (result.returncode == 0) == (mismatch == "none"), result.stderr


def test_workflow_requires_approval_revalidation_and_absence_before_oidc() -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    # PyYAML's YAML 1.1 parser reads the Actions `on` key as True.
    assert workflow[True] == {"release": {"types": ["published"]}}
    assert workflow["concurrency"] == {"group": "publish-crates-la-stack", "cancel-in-progress": False}
    jobs = workflow["jobs"]
    assert jobs["publish"]["needs"] == "verify"
    assert jobs["publish"]["environment"] == "crates-io"
    assert jobs["verify"]["permissions"] == {"contents": "read", "checks": "read"}
    assert jobs["publish"]["permissions"] == {"contents": "read", "checks": "read", "id-token": "write"}
    assert "github.repository == 'acgetchell/la-stack'" in jobs["verify"]["if"]
    for job in jobs.values():
        checkout = job["steps"][0]
        assert checkout["with"] == {"ref": "${{ github.sha }}", "fetch-depth": 0, "persist-credentials": False}
        assert job["steps"][1]["with"]["cache"] == "false"
    steps = jobs["publish"]["steps"]
    gate = steps[2]["run"]
    assert 'just release-gate "$RELEASE_TAG"' in gate
    assert 'release registry "$RELEASE_TAG" --expect absent' in gate
    assert steps[3]["id"] == "auth"
    assert steps[4]["run"].endswith("cargo publish --locked --package la-stack --registry crates-io")
    assert steps[5]["if"] == "${{ !cancelled() && steps.auth.outcome == 'success' }}"
    assert steps[5]["run"] == 'just release-verify "$RELEASE_TAG"'
    assert "continue-on-error" not in WORKFLOW.read_text(encoding="utf-8")
    commands = dry_run(REPO_ROOT, "release-package").stderr
    assert commands.count("--package la-stack --registry crates-io --dry-run") == 2
    assert "--features exact" in commands
    assert "--allow-dirty" not in commands


@pytest.mark.parametrize(("failure", "expected_exit"), [("none", 0), ("evidence", 1), ("present", 1), ("unavailable", 2)])
def test_actual_preupload_shell_requires_evidence_and_proven_absence(tmp_path: Path, failure: str, expected_exit: int) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    script = workflow["jobs"]["publish"]["steps"][2]["run"]
    marker = tmp_path / "upload-reached"
    gate = "just() { return 1; }" if failure == "evidence" else "just() { return 0; }"
    registry_exit = expected_exit if failure in {"present", "unavailable"} else 0
    registry = f"uv() {{ return {registry_exit}; }}"
    script = f"{gate}\n{registry}\n{script}\ntouch {shlex.quote(marker.as_posix())}"
    result = run_command(
        "bash",
        ["--noprofile", "--norc", "-euc", script],
        env={**{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}}, "RELEASE_TAG": TAG},
        check=False,
    )
    assert result.returncode == expected_exit, result.stderr
    assert marker.exists() == (failure == "none")


def test_native_packaging_failure_stops_before_further_validation(tmp_path: Path) -> None:
    script = dry_run(REPO_ROOT, "release-package").stderr
    calls = tmp_path / "calls"
    # Cargo/dependency failure is propagated by the real recipe's strict shell.
    runner = f"uv() {{ echo invoked >> {shlex.quote(calls.as_posix())}; return 17; }}\n"
    result = run_command("bash", ["--noprofile", "--norc", "-euc", runner + script], check=False)
    assert result.returncode == 17
    assert calls.read_text(encoding="utf-8").splitlines() == ["invoked"]
