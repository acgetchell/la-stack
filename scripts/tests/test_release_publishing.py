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


def test_publishing_policy_requires_native_scientific_gates_and_the_criterion_asset() -> None:
    assert POLICY.publishing is not None
    assert (POLICY.publishing.registry, POLICY.publishing.package, POLICY.publishing.repository) == ("crates-io", "la-stack", "acgetchell/la-stack")
    assert tuple((check.name, check.app_id) for check in POLICY.publishing.required_checks) == tuple((name, 15368) for name in CHECKS)
    assert POLICY.publishing.required_assets == ("la-stack-v{version}-criterion-baseline.tar.gz",)


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
    return {"checks": checks}


def test_installed_gate_accepts_actual_metadata_and_required_native_results(release_event: dict[str, Any]) -> None:
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 0


@pytest.mark.parametrize("name", CHECKS)
def test_each_configured_native_check_is_required(release_event: dict[str, Any], capsys: pytest.CaptureFixture[str], name: str) -> None:
    check = next(item for item in release_event["checks"] if item["name"] == name)
    check["name"] = "unrelated check"
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 1
    assert f"required check {name!r} lacks one unambiguous successful result on {SHA}" in capsys.readouterr().err


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


def test_release_notes_stdout_excludes_setup_diagnostics() -> None:
    # Include the dependency commands that polluted the generated notes.
    script = dry_run(REPO_ROOT, "release-notes", [TAG]).stderr
    runner = (
        "uv() {\n"
        '  if [[ "$*" == "run --locked --no-sync --no-python-downloads research-repo-tools deps check-uv" ]]; then\n'
        "    printf '%s\\n' 'fixture setup diagnostic'\n"
        '  elif [[ "$*" == "sync --locked --group dev" ]]; then\n'
        "    return 0\n"
        "  else\n"
        '    [[ "$1 $2 $3 $4 $5" == "run --locked --group dev research-repo-tools" ]] || return 90\n'
        "    shift 5\n"
        f'    {shlex.quote(sys.executable)} -m research_repo_tools --root {shlex.quote(str(REPO_ROOT))} "$@"\n'
        "  fi\n"
        "}\n"
    )
    result = run_command("bash", ["--noprofile", "--norc", "-euc", runner + script], cwd=REPO_ROOT, check=False)
    assert result.returncode == 0, result.stderr
    assert result.stdout == changelog.notes(POLICY, TAG)[0]
    assert "fixture setup diagnostic" in result.stderr


def test_publish_recipe_passes_explicit_approval(monkeypatch: pytest.MonkeyPatch) -> None:
    requested: list[str] = []

    def publish(settings: config.Config, tag: str) -> None:
        assert settings.root == REPO_ROOT
        requested.append(tag)

    # Exercise the installed argument parser and dispatch, replacing only the
    # irreversible publication boundary with an in-memory request recorder.
    monkeypatch.setattr(release_publishing, "publish_reviewed_release", publish)
    # The publication command follows the separately tested notes precheck.
    command = shlex.split(dry_run(REPO_ROOT, "release-publish", [TAG]).stderr.splitlines()[-1])
    assert command[:7] == ["uv", "run", "--locked", "--no-sync", "research-repo-tools", "release", "publish"]
    assert main(["--root", str(REPO_ROOT), *command[5:]]) == 0
    assert requested == [TAG]


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
    assert jobs["verify"]["steps"][-1]["run"] == "just release-package"
    assert "CARGO_REGISTRY_TOKEN" not in yaml.safe_dump({"env": workflow["env"], "verify": jobs["verify"]})
    for job in jobs.values():
        checkout = job["steps"][0]
        assert checkout["with"] == {"ref": "${{ github.sha }}", "fetch-depth": 0, "persist-credentials": False}
        assert job["steps"][1]["with"]["cache"] == "false"
    steps = jobs["publish"]["steps"]
    gate = steps[2]["run"]
    assert 'just release-gate "$RELEASE_TAG"' in gate
    assert 'release registry "$RELEASE_TAG" --expect absent' in gate
    assert steps[3]["id"] == "auth"
    assert steps[4]["run"].endswith("cargo publish --locked --package la-stack --registry crates-io --no-verify")
    assert [index for index, step in enumerate(steps) if "CARGO_REGISTRY_TOKEN" in step.get("env", {})] == [4]
    assert steps[5]["if"] == "${{ !cancelled() && steps.auth.outcome == 'success' }}"
    assert steps[5]["run"] == 'just release-verify "$RELEASE_TAG"'
    assert "continue-on-error" not in WORKFLOW.read_text(encoding="utf-8")
    commands = dry_run(REPO_ROOT, "release-package").stderr
    assert commands.count("--package la-stack --registry crates-io --dry-run") == 2
    assert "--features exact" in commands
    assert "--allow-dirty" not in commands
    assert "--no-verify" not in commands


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
