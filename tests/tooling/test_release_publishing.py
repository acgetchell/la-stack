"""Exercise la-stack's published release contract without credentials or uploads."""

import json
import os
import runpy
import shlex
import subprocess
import sys
import tomllib
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
import yaml
from research_repo_tools import changelog, config, dependency_provenance, registry, release_publishing
from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_command

from tests.tooling.support_evidence import install_clean_checkout
from tests.tooling.support_toolchain import TOOLING_REQUIREMENT

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github/workflows/publish-crates.yml"
SHA = "a" * 40
POLICY = config.load(root=REPO_ROOT)
VERSION = tomllib.loads((REPO_ROOT / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
TAG = f"v{VERSION}"
CHECKS = ("build (ubuntu-latest)", "build (macos-latest)", "build (windows-latest)", "Clippy SARIF Analysis", "Repository Rule SARIF Analysis")


def test_publishing_policy_requires_native_scientific_gates_and_the_criterion_asset() -> None:
    assert POLICY.release.published_dependencies == ("research-repo-tools",)
    assert POLICY.publishing is not None
    assert (POLICY.publishing.registry, POLICY.publishing.package, POLICY.publishing.repository) == ("crates-io", "la-stack", "acgetchell/la-stack")
    assert tuple((check.name, check.app_id) for check in POLICY.publishing.required_checks) == tuple((name, 15368) for name in CHECKS)
    assert POLICY.publishing.required_assets == ("la-stack-v{version}-criterion-baseline.tar.gz",)
    assert POLICY.publishing.expected_title == "tag"
    assert POLICY.publishing.generated_notes
    assert POLICY.publishing.require_immutable


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

    # These fixtures model a published supplier while testing native release
    # identities. Separate consumer tests exercise the real development-pin gate.
    def published(root: Path, *, packages: tuple[str, ...] | None = None) -> tuple[dependency_provenance.DependencyIdentity, ...]:
        assert root == REPO_ROOT
        assert packages == ("research-repo-tools",)
        return (dependency_provenance.DependencyIdentity("research-repo-tools", "0.2.8", "https://pypi.org/simple"),)

    monkeypatch.setattr(dependency_provenance, "check_release_dependencies", published)
    install_clean_checkout(monkeypatch, REPO_ROOT, SHA)
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
        title=TAG,
        body=changelog.notes(POLICY, TAG)[0],
    )
    monkeypatch.setattr(release_publishing, "lookup_release", lambda *_args: release)
    return {"checks": checks, "release": release}


def test_installed_gate_accepts_actual_metadata_and_required_native_results(release_event: dict[str, Any]) -> None:
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 0


@pytest.mark.parametrize("name", CHECKS)
def test_each_configured_native_check_is_required(release_event: dict[str, Any], capsys: pytest.CaptureFixture[str], name: str) -> None:
    check = next(item for item in release_event["checks"] if item["name"] == name)
    check["name"] = "unrelated check"
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == 1
    diagnostic = capsys.readouterr().err
    assert f"Unmet release checks on exact commit {SHA}:" in diagnostic
    assert f"required check {name!r} (app ID 15368): missing" in diagnostic


@pytest.mark.parametrize("immutable", [True, False, None])
def test_installed_gate_requires_immutable_publication(
    release_event: dict[str, Any], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], immutable: bool | None
) -> None:
    release = replace(release_event["release"], immutable=immutable)
    monkeypatch.setattr(release_publishing, "lookup_release", lambda *_args: release)
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == (0 if immutable is True else 1)
    if immutable is not True:
        assert ("must be immutable" if immutable is False else "lifecycle fields") in capsys.readouterr().err


@pytest.mark.parametrize("mismatch", ["none", "title", "notes", "crlf"])
def test_installed_gate_checks_generated_body_and_title(
    release_event: dict[str, Any], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], mismatch: str
) -> None:
    release = release_event["release"]
    body = release.body.replace("\n", "\r\n") if mismatch == "crlf" else "wrong" if mismatch == "notes" else release.body
    release = replace(release, title="wrong" if mismatch == "title" else TAG, body=body)
    monkeypatch.setattr(release_publishing, "lookup_release", lambda *_args: release)
    assert main(["--root", str(REPO_ROOT), "release", "gate", TAG]) == (0 if mismatch in {"none", "crlf"} else 1)
    if mismatch in {"title", "notes"}:
        assert ("release title differs" if mismatch == "title" else "release body differs") in capsys.readouterr().err


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

    def publish(settings: config.Config, tag: str, *, reuse_evidence: object = None) -> None:
        assert settings.root == REPO_ROOT
        assert reuse_evidence is None
        requested.append(tag)

    # Exercise the installed argument parser and dispatch, replacing only the
    # irreversible publication boundary with an in-memory request recorder.
    monkeypatch.setattr(release_publishing, "publish_reviewed_release", publish)
    # Shared publishing policy owns the separately tested title/notes assertions.
    command = shlex.split(dry_run(REPO_ROOT, "release-publish", [TAG, "/independent-owner"]).stderr.splitlines()[-1])
    assert command[:7] == ["uv", "run", "--locked", "--no-sync", "research-repo-tools", "release", "publish"]
    assert main(["--root", str(REPO_ROOT), *command[5:]]) == 0
    assert requested == [TAG]


def test_workflow_requires_approval_revalidation_and_absence_before_oidc() -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    assert set(workflow[True]) == {"workflow_dispatch"}
    assert set(workflow[True]["workflow_dispatch"]["inputs"]) == {"tag", "candidate", "acceptance"}
    assert workflow["concurrency"] == {"group": "publish-crates-la-stack", "cancel-in-progress": False}
    jobs = workflow["jobs"]
    assert jobs["publish"]["needs"] == "verify"
    assert jobs["publish"]["environment"] == "crates-io"
    assert jobs["verify"]["permissions"] == {"contents": "read", "checks": "read", "actions": "read"}
    assert jobs["publish"]["permissions"] == {"contents": "read", "checks": "read", "actions": "read", "id-token": "write"}
    assert "github.repository == 'acgetchell/la-stack'" in jobs["verify"]["if"]
    assert "github.triggering_actor == 'acgetchell'" in jobs["verify"]["if"]
    assert "CARGO_REGISTRY_TOKEN" not in yaml.safe_dump(jobs["verify"])
    for job in jobs.values():
        assert job["steps"][0]["with"] == {"ref": "${{ github.sha }}", "fetch-depth": 0, "persist-credentials": False}
        assert job["steps"][1]["with"]["cache"] == "false"
        gate = next(step for step in job["steps"] if "just release-receive" in step.get("run", ""))
        assert gate["env"]["CANDIDATE_TRANSFER"] == "${{ inputs.candidate }}"
        assert gate["env"]["CANDIDATE_ACCEPTANCE"] == "${{ inputs.acceptance }}"
        assert 'just release-receive "$RUNNER_TEMP/accepted-candidate" --published' in gate["run"]
        assert job["steps"][1]["with"]["tool"] == "cargo"
    steps = jobs["publish"]["steps"]
    assert jobs["verify"]["steps"][3]["id"] == "approval"
    assert jobs["verify"]["steps"][3]["with"]["path"].endswith("/acceptance.json")
    assert steps[2]["env"]["ARTIFACT_DIGEST"] == "${{ needs.verify.outputs.artifact-digest }}"
    assert steps[3]["with"]["artifact-ids"] == "${{ needs.verify.outputs.artifact-id }}"
    assert '--expected "$RUNNER_TEMP/approved-acceptance/acceptance.json"' in steps[4]["run"]
    assert 'release registry "$RELEASE_TAG" --expect absent' in steps[4]["run"]
    assert steps[5]["id"] == "auth"
    assert steps[6]["run"].endswith("toolchain run --tool cargo -- cargo publish --locked --package la-stack --registry crates-io --no-verify")
    assert [index for index, step in enumerate(steps) if "CARGO_REGISTRY_TOKEN" in step.get("env", {})] == [6]
    assert steps[7]["if"] == "${{ !cancelled() && steps.auth.outcome == 'success' }}"
    assert steps[7]["run"] == 'just release-verify "$RELEASE_TAG"'
    assert "continue-on-error" not in WORKFLOW.read_text(encoding="utf-8")
    commands = dry_run(REPO_ROOT, "publish-check").stderr
    assert commands.count("--package la-stack --registry crates-io --dry-run") == 2
    assert "--features exact" in commands
    assert "--allow-dirty" not in commands
    assert "--no-verify" not in commands
    preparation = dry_run(REPO_ROOT, "release-dependency-builds").stderr
    managed = [shlex.split(line) for line in preparation.splitlines() if "toolchain run" in line]
    assert len(managed) == 2
    assert all(command[command.index("toolchain") :][:5] == ["toolchain", "run", "--tool", "cargo", "--"] for command in managed)


@pytest.mark.parametrize(
    ("failure", "expected_exit"), [("none", 0), ("dependency-builds", 1), ("evidence", 1), ("candidate", 1), ("present", 1), ("unavailable", 2)]
)
def test_actual_preupload_shell_requires_evidence_and_proven_absence(tmp_path: Path, failure: str, expected_exit: int) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    script = next(
        step["run"] for step in workflow["jobs"]["publish"]["steps"] if step.get("name") == "Recheck evidence and registry absence after deployment approval"
    )
    marker = tmp_path / "upload-reached"
    gate = "just() { return 1; }" if failure == "evidence" else "just() { return 0; }"
    if failure == "dependency-builds":
        gate = 'just() { [[ "$1" != release-dependency-builds ]]; }'
    if failure == "candidate":
        gate = 'just() { [[ "$1" != release-receive ]]; }'
    registry_exit = expected_exit if failure in {"present", "unavailable"} else 0
    registry = f"uv() {{ return {registry_exit}; }}"
    script = f"{gate}\n{registry}\n{script}\ntouch {shlex.quote(marker.as_posix())}"
    result = run_command(
        "bash",
        ["--noprofile", "--norc", "-euc", script],
        env={
            **{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}},
            "RELEASE_TAG": TAG,
            "CANDIDATE_TRANSFER": '{"schema": 1}',
            "RUNNER_TEMP": tmp_path.as_posix(),
        },
        check=False,
    )
    assert result.returncode == expected_exit, result.stderr
    assert marker.exists() == (failure == "none")


@pytest.mark.parametrize("failure", ["none", "attempt", "identifier", "digest", "expired", "different-digest", "different-run"])
def test_actual_approval_artifact_guard_rejects_drift(tmp_path: Path, failure: str) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    script = next(step["run"] for step in workflow["jobs"]["publish"]["steps"] if step.get("name") == "Require this exact approval receipt")
    metadata = {
        "expired": failure == "expired",
        "digest": "sha256:" + ("b" if failure == "different-digest" else "a") * 64,
        "workflow_run": {"id": 99 if failure == "different-run" else 42},
    }
    supplied = tmp_path / "artifact metadata.json"
    supplied.write_text(json.dumps(metadata), encoding="utf-8", newline="\n")
    marker = tmp_path / "approval-reached"
    provider = f'gh() {{ [[ "$*" == "api repos/acgetchell/la-stack/actions/artifacts/7" ]] && cat {shlex.quote(supplied.as_posix())}; }}'
    result = run_command(
        "bash",
        ["--noprofile", "--norc", "-euc", f"{provider}\n{script}\ntouch {shlex.quote(marker.as_posix())}"],
        env={
            **{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}},
            "VALIDATED_ATTEMPT": "1" if failure == "attempt" else "2",
            "GITHUB_RUN_ATTEMPT": "2",
            "ARTIFACT_ID": "07" if failure == "identifier" else "7",
            "ARTIFACT_DIGEST": "invalid" if failure == "digest" else "a" * 64,
            "GITHUB_RUN_ID": "42",
            "GITHUB_REPOSITORY": "acgetchell/la-stack",
        },
        check=False,
    )
    assert result.returncode == (0 if failure == "none" else 1), result.stderr
    assert marker.exists() == (failure == "none")


def test_main_ci_seals_only_complete_native_and_cargo_validation() -> None:
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/ci.yml").read_bytes())
    jobs = workflow["jobs"]
    platforms = ["ubuntu-latest", "macos-latest", "windows-latest"]
    for job in ("prepare", "native", "build"):
        assert jobs[job]["strategy"]["matrix"]["os"] == platforms
        assert jobs[job]["strategy"]["fail-fast"] is False
    assert jobs["prepare"]["runs-on"] == jobs["native"]["runs-on"] == "${{ matrix.os }}"
    assert jobs["prepare"]["timeout-minutes"] == 90
    assert jobs["prepare"]["steps"][-1] == {"name": "Run CI policy checks", "run": "just ci-policy"}
    assert jobs["native"]["needs"] == "prepare"
    assert jobs["native"]["name"] == "native-${{ matrix.stage }} (${{ matrix.os }})"
    assert jobs["native"]["strategy"]["matrix"]["stage"] == ["rust", "python", "candidate", "dependencies"]
    assert jobs["native"]["timeout-minutes"] == 120
    assert jobs["native"]["steps"][-1] == {
        "name": "Run native validation",
        "shell": "bash",
        "env": {"CI_STAGE": "${{ matrix.stage }}"},
        "run": 'just "ci-$CI_STAGE"',
    }
    assert jobs["build"]["needs"] == ["prepare", "native"]
    # Required checks must execute and reject cancellation, rather than skip.
    assert jobs["build"]["if"] == "${{ always() }}"
    assert jobs["build"]["runs-on"] == "ubuntu-latest"
    assert jobs["build"]["timeout-minutes"] == 5
    package = jobs["release-package"]
    assert package["steps"][1]["with"]["cache"] == "false"
    assert package["steps"][-1] == {"name": "Verify default and exact Cargo packages", "run": "just ci-dependencies\njust publish-check\n"}
    evidence = jobs["release-evidence"]
    assert evidence["needs"] == ["build", "release-package"]
    assert evidence["if"] == "github.event_name == 'push' && github.ref == 'refs/heads/main' && github.repository == 'acgetchell/la-stack'"
    assert evidence["permissions"] == {"contents": "read", "actions": "read", "attestations": "write", "id-token": "write"}
    steps = {step["name"]: step for step in evidence["steps"] if "name" in step}
    capture = steps["Capture successful native validation identities"]["run"]
    assert steps["Capture successful native validation identities"]["id"] == "evidence"
    for name in ("Sign original main validation evidence", "Retain the signature bundle", "Upload signed validation evidence"):
        assert steps[name]["if"] == "steps.evidence.outputs.eligible == 'true'"
    assert shlex.split(capture.split("<<'PY'", 1)[0].replace("\\\n", "")) == [
        "uvx",
        "--python",
        "3.15",
        "--from",
        TOOLING_REQUIREMENT,
        "python",
        "-I",
        "-",
    ]
    assert '"release", "evidence", "create", "tooling/release-evidence.toml"' in capture
    assert '"--output", "target/release-evidence/validation.json"' in capture
    assert steps["Sign original main validation evidence"]["with"]["subject-path"] == "target/release-evidence/validation.json"
    upload = steps["Upload signed validation evidence"]["with"]
    assert upload == {
        "name": "release-evidence-${{ github.run_id }}-${{ github.run_attempt }}",
        "path": "target/release-evidence/",
        "if-no-files-found": "error",
        "retention-days": 30,
    }


@pytest.mark.parametrize(
    ("prepare", "native"),
    [
        ("success", "success"),
        ("failure", "success"),
        ("skipped", "success"),
        ("cancelled", "success"),
        ("", "success"),
        ("success", "failure"),
        ("success", "skipped"),
        ("success", "cancelled"),
        ("success", ""),
    ],
)
def test_main_ci_required_checks_reject_incomplete_native_stages(prepare: str, native: str) -> None:
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/ci.yml").read_bytes())
    step = workflow["jobs"]["build"]["steps"][0]
    assert step["env"] == {"PREPARE_RESULT": "${{ needs.prepare.result }}", "NATIVE_RESULT": "${{ needs.native.result }}"}
    environment = {**os.environ, "PREPARE_RESULT": prepare, "NATIVE_RESULT": native}
    completed = run_command("bash", ["--noprofile", "--norc", "-euc", step["run"]], env=environment, check=False)
    assert completed.returncode == (0 if prepare == native == "success" else 1), completed.stderr


@pytest.mark.parametrize("source", ["development", "published", "rejected"])
def test_main_ci_creates_no_release_evidence_for_development_pins(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str) -> None:
    """Run the actual workflow adapter; only a successful release gate enables signing."""
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/ci.yml").read_bytes())
    step = next(item for item in workflow["jobs"]["release-evidence"]["steps"] if item.get("id") == "evidence")
    script = tmp_path / "capture.py"
    script.write_text(step["run"].split("<<'PY'\n", 1)[1].rsplit("\nPY", 1)[0], encoding="utf-8")
    output = tmp_path / "github output"
    original = b"prior=value\n"
    output.write_bytes(original)
    for key, value in {"GITHUB_OUTPUT": str(output), "GITHUB_SHA": SHA, "GITHUB_RUN_ID": "11", "GITHUB_RUN_ATTEMPT": "1"}.items():
        monkeypatch.setenv(key, value)
    pin = dependency_provenance.DependencyIdentity("research-repo-tools", "0.2.8", "git", "a" * 40) if source == "development" else None
    monkeypatch.setattr(dependency_provenance, "check_development_dependency", lambda _root: pin)
    calls: list[list[str]] = []

    def create(arguments: list[str]) -> int:
        calls.append(arguments)
        return 1 if source == "rejected" else 0

    monkeypatch.setattr("research_repo_tools.cli.main", create)
    with pytest.raises(SystemExit) as stopped:
        runpy.run_path(str(script), run_name="__main__")
    assert stopped.value.code == (1 if source == "rejected" else 0)
    assert len(calls) == (0 if source == "development" else 1)
    if calls:
        assert calls[0][:4] == ["release", "evidence", "create", "tooling/release-evidence.toml"]
    data = output.read_bytes()
    assert data.startswith(original)
    assert (b"eligible" in data) == (source == "published")


@pytest.mark.parametrize("outcome", ["visible", "draft", "missing-asset", "incomplete-asset", "absent-version"])
def test_release_visibility_checks_remote_state_without_local_release_preparation(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, outcome: str) -> None:
    script = tmp_path / "visibility.py"
    script.write_text(dry_run(REPO_ROOT, "release-verify", [TAG]).stderr, encoding="utf-8")
    monkeypatch.chdir(REPO_ROOT)
    monkeypatch.setattr(sys, "argv", [str(script), TAG])
    calls: list[str] = []
    release = release_publishing.GitHubRelease(
        repository="acgetchell/la-stack",
        tag=TAG,
        identifier=42,
        draft=outcome == "draft",
        immutable=outcome != "draft",
        prerelease=False,
        assets=() if outcome == "missing-asset" else ((f"la-stack-{TAG}-criterion-baseline.tar.gz", 7, 100, None),),
        title=TAG,
        asset_states=((7, "starter" if outcome == "incomplete-asset" else "uploaded"),),
    )

    def github(root: Path, repository: str, tag: str) -> release_publishing.GitHubRelease:
        assert (root, repository, tag) == (REPO_ROOT, "acgetchell/la-stack", TAG)
        calls.append("github")
        return release

    def visible(destination: str, package: str, version: str, *, attempts: int) -> registry.RegistryVersion:
        assert (destination, package, version, attempts) == ("crates-io", "la-stack", VERSION, 7)
        calls.append("registry")
        return registry.RegistryVersion("crates-io", package, version, outcome != "absent-version")

    monkeypatch.setattr("research_repo_tools.release_assets.lookup_release", github)
    monkeypatch.setattr(registry, "wait_for_version", visible)
    monkeypatch.setattr(dependency_provenance, "check_release_dependencies", lambda *_args, **_kwargs: pytest.fail("visibility repeated local provenance"))
    monkeypatch.setattr(changelog, "notes", lambda *_args: pytest.fail("visibility regenerated local notes"))
    if outcome == "visible":
        runpy.run_path(str(script), run_name="__main__")
    else:
        with pytest.raises(SystemExit, match=r"immutable|Missing completed release assets|not visible"):
            runpy.run_path(str(script), run_name="__main__")
    assert calls == (["github", "registry"] if outcome in {"visible", "absent-version"} else ["github"])


def test_native_packaging_failure_stops_before_further_validation(tmp_path: Path) -> None:
    script = dry_run(REPO_ROOT, "publish-check").stderr
    calls = tmp_path / "calls"
    # Cargo/dependency failure is propagated by the real recipe's strict shell.
    runner = f"uv() {{ echo invoked >> {shlex.quote(calls.as_posix())}; return 17; }}\n"
    result = run_command("bash", ["--noprofile", "--norc", "-euc", runner + script], check=False)
    assert result.returncode == 17
    assert calls.read_text(encoding="utf-8").splitlines() == ["invoked"]
