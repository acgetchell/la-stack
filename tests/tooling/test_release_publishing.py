"""Exercise la-stack's published release contract without credentials or uploads."""

import io
import json
import os
import shlex
import subprocess
import sys
import tomllib
from contextlib import redirect_stderr
from pathlib import Path
from typing import Any

import pytest
import yaml
from research_repo_tools import changelog, config, release_publishing
from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_command

from tests.tooling.test_release_evidence import PIN_JSON, EvidenceProvider
from tooling import release_evidence

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
    assert jobs["verify"]["permissions"] == {"contents": "read", "checks": "read", "actions": "read"}
    assert jobs["publish"]["permissions"] == {"contents": "read", "checks": "read", "actions": "read", "id-token": "write"}
    assert "github.repository == 'acgetchell/la-stack'" in jobs["verify"]["if"]
    assert jobs["verify"]["steps"][-1]["id"] == "evidence"
    assert jobs["verify"]["outputs"] == {"ci-evidence": "${{ steps.evidence.outputs.identity }}"}
    assert "CARGO_REGISTRY_TOKEN" not in yaml.safe_dump({"env": workflow["env"], "verify": jobs["verify"]})
    for job in jobs.values():
        checkout = job["steps"][0]
        assert checkout["with"] == {"ref": "${{ github.sha }}", "fetch-depth": 0, "persist-credentials": False}
        assert job["steps"][1]["with"]["cache"] == "false"
    steps = jobs["publish"]["steps"]
    gate = steps[2]["run"]
    assert 'just release-gate "$RELEASE_TAG"' in gate
    assert 'test -n "$CI_EVIDENCE"' in gate
    assert 'just release-evidence --identity "$CI_EVIDENCE"' in gate
    assert steps[2]["env"]["CI_EVIDENCE"] == "${{ needs.verify.outputs.ci-evidence }}"
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


def run_evidence_fixture() -> int:
    """Run the real verifier in a child process with simulated provider responses."""
    root = Path(os.environ["FIXTURE_ROOT"])
    provider = EvidenceProvider(release_evidence, root)
    if os.environ.get("FIXTURE_FAILURE") == "failed-job":
        provider.jobs[3]["steps"][0]["conclusion"] = "failure"
    newline = os.environ.get("FIXTURE_DIAGNOSTIC_NEWLINE", os.linesep)
    sys.stdout.buffer.write(("fixture stdout diagnostic" + newline).encode("utf-8"))
    sys.stdout.flush()
    print("fixture stderr diagnostic", file=sys.stderr)
    diagnostic = io.StringIO()
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(release_evidence, "ROOT", root)
        patch.setattr(release_evidence, "api", provider.api)
        patch.setattr(release_evidence, "gh", provider.gh)
        patch.setattr(release_evidence, "download_artifact", provider.download)
        patch.setattr(release_evidence, "run_git_command", provider.git)
        with redirect_stderr(diagnostic):
            try:
                status = release_evidence.main()
            except SystemExit as error:
                status = error.code
    assert isinstance(status, int)
    print(diagnostic.getvalue(), end="", file=sys.stderr)
    # Just's quiet mode hides stderr. Record the actual verifier outcome so an
    # import/shell failure cannot masquerade as the rejection under test.
    result = {"status": status, "diagnostic": diagnostic.getvalue(), "import_root": sys.path[0], "arguments": sys.argv[1:]}
    (root / "verifier-result.json").write_text(json.dumps(result), encoding="utf-8")
    return status


@pytest.fixture
def evidence_runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep real workflow shell, Just, uv, Python, verifier and output writer."""
    (tmp_path / "justfile").write_text(f"import '{(REPO_ROOT / 'justfile').as_posix()}'\n", encoding="utf-8")
    (tmp_path / "Cargo.toml").write_text('[package]\nname = "la-stack"\nversion = "0.5.42"\n', encoding="utf-8")
    for name in ("BASH_ENV", "SHELLOPTS", "BASHOPTS", "JUST_QUIET", "JUST_DRY_RUN", "FIXTURE_FAILURE", "FIXTURE_DIAGNOSTIC_NEWLINE", "PYTHONPATH"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("FIXTURE_ROOT", str(tmp_path))
    monkeypatch.setenv("FIXTURE_REPO_ROOT", str(REPO_ROOT))
    monkeypatch.setenv("UV_PROJECT", str(REPO_ROOT))
    monkeypatch.setenv("UV_PROJECT_ENVIRONMENT", sys.prefix)
    # Just executes the entry point relative to the fixture checkout. Supply
    # provider responses there while retaining the actual native uv/Python launch.
    bootstrap = tmp_path / "tooling/release_evidence.py"
    bootstrap.parent.mkdir()
    bootstrap.write_text(
        'import os, sys\nsys.path.insert(0, os.environ["FIXTURE_REPO_ROOT"])\n'
        "from tests.tooling.test_release_publishing import run_evidence_fixture\nsys.exit(run_evidence_fixture())\n",
        encoding="utf-8",
        newline="\n",
    )


def run_evidence_shell(script: str, root: Path) -> subprocess.CompletedProcess[str]:
    """Execute workflow source from a file, as Actions does, with explicit env."""
    source = root / "workflow.sh"
    source.write_text(script, encoding="utf-8", newline="\n")
    # Windows putenv removes empty values from the inherited native environment.
    # Build the child environment from the mapping to keep empty distinct from unset.
    return run_command("bash", ["--noprofile", "--norc", "-eu", source.as_posix()], cwd=root, env=dict(os.environ), check=False, timeout=60)


@pytest.mark.usefixtures("evidence_runner")
@pytest.mark.parametrize("quiet", [False, True], ids=["normal", "quiet"])
@pytest.mark.parametrize("outcome", ["valid", "failed-job", "missing-file", "directory", "unset-output", "empty-output"])
def test_evidence_export_uses_real_just_and_installed_writer_without_stdout_capture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, quiet: bool, outcome: str
) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    script = workflow["jobs"]["verify"]["steps"][-1]["run"]
    output = tmp_path / "github output"
    original = b"prior=value\r\n"
    if outcome == "directory":
        output.mkdir()
    elif outcome != "missing-file":
        output.write_bytes(original)
    monkeypatch.setenv("GITHUB_OUTPUT", "" if outcome == "empty-output" else output.as_posix())
    if outcome == "unset-output":
        monkeypatch.delenv("GITHUB_OUTPUT")
    monkeypatch.setenv("JUST_QUIET", str(quiet).lower())
    monkeypatch.setenv("FIXTURE_FAILURE", outcome)
    marker = tmp_path / "verification-passed"
    result = run_evidence_shell(script + f"\ntouch {shlex.quote(marker.as_posix())}", tmp_path)
    assert result.returncode == (0 if outcome == "valid" else 1), result.stderr
    assert marker.exists() == (outcome == "valid")
    report = tmp_path / "verifier-result.json"
    if outcome == "unset-output":
        assert not report.exists()
        assert "GITHUB_OUTPUT: unbound variable" in result.stderr
    else:
        assert report.exists(), result.stderr
        observed = json.loads(report.read_text(encoding="utf-8"))
        assert observed["arguments"][:2] == ["verify", "--github-output"]
        assert len(observed["arguments"]) == 3
        if outcome == "empty-output":
            assert observed["arguments"][2] == ""
        else:
            assert Path(observed["arguments"][2]) == output
        assert observed["status"] == (0 if outcome == "valid" else 1)
        rejection = (
            "did not validate required step Verify default and exact Cargo packages"
            if outcome == "failed-job"
            else "GitHub Actions output destination must be an existing regular file"
        )
        assert observed["diagnostic"] == "" if outcome == "valid" else rejection in observed["diagnostic"]
    if outcome == "valid":
        assert output.read_bytes() == original + f"identity={PIN_JSON}\n".encode()
        if not quiet:
            # Diagnostics are text; the command-file bytes above must stay LF.
            assert result.stdout.splitlines() == ["fixture stdout diagnostic"]
            assert "fixture stderr diagnostic" in result.stderr
    elif outcome == "directory":
        assert output.is_dir()
        assert list(output.iterdir()) == []
    elif outcome == "missing-file":
        assert not output.exists()
    else:
        assert output.read_bytes() == original


@pytest.mark.usefixtures("evidence_runner")
@pytest.mark.parametrize(
    ("import_root", "newline"), [(r"C:\agent\test dir\native", "\n"), (r"D:\build's root\new\repo", "\r\n")], ids=["backslashes-lf", "apostrophe-crlf"]
)
def test_evidence_bootstrap_keeps_windows_import_paths_as_environment_data(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, import_root: str, newline: str
) -> None:
    # Exercise literal Windows path data on every host. PYTHONPATH supplies the
    # real module after the intentionally nonexistent path inserted by bootstrap.
    monkeypatch.setenv("FIXTURE_REPO_ROOT", import_root)
    monkeypatch.setenv("PYTHONPATH", str(REPO_ROOT))
    monkeypatch.setenv("FIXTURE_DIAGNOSTIC_NEWLINE", newline)
    output = tmp_path / "github output"
    output.touch()
    monkeypatch.setenv("GITHUB_OUTPUT", output.as_posix())
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    script = workflow["jobs"]["verify"]["steps"][-1]["run"]
    result = run_evidence_shell(script, tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout == "fixture stdout diagnostic" + newline
    observed = json.loads((tmp_path / "verifier-result.json").read_text(encoding="utf-8"))
    assert observed["status"] == 0
    assert observed["diagnostic"] == ""
    assert observed["import_root"] == import_root
    assert observed["arguments"][:2] == ["verify", "--github-output"]
    assert len(observed["arguments"]) == 3
    assert Path(observed["arguments"][2]) == output
    assert output.read_bytes() == f"identity={PIN_JSON}\n".encode()


@pytest.mark.usefixtures("evidence_runner")
@pytest.mark.parametrize("quiet", [False, True], ids=["normal", "quiet"])
@pytest.mark.parametrize("pin", ["valid", "empty", "invalid", "changed"])
def test_exported_identity_is_revalidated_by_real_just_before_registry_and_token(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, quiet: bool, pin: str
) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    output = tmp_path / "github output"
    output.touch()
    monkeypatch.setenv("GITHUB_OUTPUT", output.as_posix())
    monkeypatch.setenv("JUST_QUIET", str(quiet).lower())
    exported = run_evidence_shell(workflow["jobs"]["verify"]["steps"][-1]["run"], tmp_path)
    assert exported.returncode == 0, exported.stderr
    report = tmp_path / "verifier-result.json"
    assert json.loads(report.read_text(encoding="utf-8"))["status"] == 0
    report.unlink()
    record = output.read_text(encoding="utf-8")
    assert record == f"identity={PIN_JSON}\n"
    identity = record.removeprefix("identity=").removesuffix("\n")
    variants = {"valid": identity, "empty": "", "invalid": "{", "changed": identity.replace('"run_attempt":2', '"run_attempt":1')}
    monkeypatch.setenv("CI_EVIDENCE", variants[pin])
    monkeypatch.setenv("RELEASE_TAG", TAG)
    registry = tmp_path / "registry-reached"
    token = tmp_path / "token-reached"
    # The other release and registry gates are independently covered below.
    # Keep evidence verification real and stop before any credentials or upload.
    runner = f'uv() {{\n  [[ "$4 $5 $6" == "research-repo-tools release registry" ]] || return 90\n  touch {shlex.quote(registry.as_posix())}\n}}\n'
    runner += 'just() { if [[ "$1" == release-gate ]]; then return 0; fi; command just "$@"; }\n'
    script = workflow["jobs"]["publish"]["steps"][2]["run"]
    result = run_evidence_shell(runner + script + f"\ntouch {shlex.quote(token.as_posix())}", tmp_path)
    if pin != "empty":
        assert report.exists(), result.stderr
        observed = json.loads(report.read_text(encoding="utf-8"))
        assert observed["arguments"] == ["verify", "--identity", variants[pin]]
    assert result.returncode == (0 if pin == "valid" else 1), result.stderr
    assert registry.exists() == (pin == "valid")
    assert token.exists() == (pin == "valid")
    if pin == "empty":
        assert not report.exists()
    else:
        observed = json.loads(report.read_text(encoding="utf-8"))
        assert observed["status"] == (0 if pin == "valid" else 1)
        rejection = "Expecting property name" if pin == "invalid" else "selected CI attempt is stale"
        assert observed["diagnostic"] == "" if pin == "valid" else rejection in observed["diagnostic"]


@pytest.mark.usefixtures("evidence_runner")
@pytest.mark.parametrize("quiet", [False, True], ids=["normal", "quiet"])
def test_explicit_empty_identity_cannot_select_fresh_evidence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, quiet: bool) -> None:
    monkeypatch.setenv("JUST_QUIET", str(quiet).lower())
    output = tmp_path / "github output"
    original = b"prior=value\n"
    output.write_bytes(original)
    monkeypatch.setenv("GITHUB_OUTPUT", output.as_posix())
    result = run_evidence_shell('just release-evidence --identity "" --github-output "$GITHUB_OUTPUT"\n', tmp_path)
    observed = json.loads((tmp_path / "verifier-result.json").read_text(encoding="utf-8"))
    assert observed["arguments"][:3] == ["verify", "--identity", ""]
    assert result.returncode == observed["status"] == 1
    assert "Expecting value" in observed["diagnostic"]
    assert output.read_bytes() == original


@pytest.mark.parametrize(("failure", "expected_exit"), [("none", 0), ("evidence", 1), ("ci-evidence", 1), ("empty-pin", 1), ("present", 1), ("unavailable", 2)])
def test_actual_preupload_shell_requires_evidence_and_proven_absence(tmp_path: Path, failure: str, expected_exit: int) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    script = workflow["jobs"]["publish"]["steps"][2]["run"]
    marker = tmp_path / "upload-reached"
    gate = "just() { return 1; }" if failure == "evidence" else "just() { return 0; }"
    if failure == "ci-evidence":
        gate = 'just() { [[ "$1" != release-evidence ]]; }'
    registry_exit = expected_exit if failure in {"present", "unavailable"} else 0
    registry = f"uv() {{ return {registry_exit}; }}"
    script = f"{gate}\n{registry}\n{script}\ntouch {shlex.quote(marker.as_posix())}"
    result = run_command(
        "bash",
        ["--noprofile", "--norc", "-euc", script],
        env={
            **{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}},
            "RELEASE_TAG": TAG,
            "CI_EVIDENCE": "" if failure == "empty-pin" else '{"run_id": 11}',
        },
        check=False,
    )
    assert result.returncode == expected_exit, result.stderr
    assert marker.exists() == (failure == "none")


def test_main_ci_seals_only_complete_native_and_cargo_validation() -> None:
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/ci.yml").read_bytes())
    jobs = workflow["jobs"]
    assert jobs["build"]["strategy"]["matrix"]["os"] == ["ubuntu-latest", "macos-latest", "windows-latest"]
    assert jobs["build"]["steps"][-1] == {"name": "Run CI checks", "run": "just ci"}
    package = jobs["release-package"]
    assert package["steps"][1]["with"]["cache"] == "false"
    assert package["steps"][-1] == {"name": "Verify default and exact Cargo packages", "run": "just release-package"}
    evidence = jobs["release-evidence"]
    assert evidence["needs"] == ["build", "release-package"]
    assert evidence["if"] == "github.event_name == 'push' && github.ref == 'refs/heads/main' && github.repository == 'acgetchell/la-stack'"
    assert evidence["permissions"] == {"contents": "read", "actions": "read", "attestations": "write", "id-token": "write"}
    steps = {step["name"]: step for step in evidence["steps"] if "name" in step}
    assert steps["Capture successful native validation identities"]["run"].endswith("seal --output target/release-evidence/validation.json")
    assert steps["Sign original main validation evidence"]["with"]["subject-path"] == "target/release-evidence/validation.json"
    upload = steps["Upload signed validation evidence"]["with"]
    assert upload == {
        "name": "release-evidence-${{ github.run_id }}-${{ github.run_attempt }}",
        "path": "target/release-evidence/",
        "if-no-files-found": "error",
        "retention-days": 30,
    }


def test_native_packaging_failure_stops_before_further_validation(tmp_path: Path) -> None:
    script = dry_run(REPO_ROOT, "release-package").stderr
    calls = tmp_path / "calls"
    # Cargo/dependency failure is propagated by the real recipe's strict shell.
    runner = f"uv() {{ echo invoked >> {shlex.quote(calls.as_posix())}; return 17; }}\n"
    result = run_command("bash", ["--noprofile", "--norc", "-euc", runner + script], check=False)
    assert result.returncode == 17
    assert calls.read_text(encoding="utf-8").splitlines() == ["invoked"]
