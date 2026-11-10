"""Exercise release command transport and the installed retained-draft handoff."""

import json
import os
import runpy
import shlex
import subprocess
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
import yaml
from research_repo_tools import release_assets, release_drafts, release_evidence, release_validation
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_command

from tests.tooling.test_release_evidence import COMMIT, POLICY, REPOSITORY, EvidenceProvider
from tests.tooling.test_toolchain_integration import TOOLING_REQUIREMENT

ROOT = Path(__file__).resolve().parents[2]
TAG = "v0.5.42"
WORKFLOW = ROOT / ".github/workflows/prepare-release.yml"


@pytest.mark.parametrize("outcome", ["success", "automatic", "dispatch-failure", "missing", "wrong-tag", "malformed"])
def test_benchmark_recipe_forwards_one_tag_and_independent_transfer(tmp_path: Path, outcome: str) -> None:
    transfer = tmp_path / "owner's transfer $literal.json"
    data = {"schema": 1, "target": {"repository": REPOSITORY, "tag": TAG, "expected_title": TAG}}
    if outcome == "wrong-tag":
        data["target"]["tag"] = "v9.9.9"
    transfer.write_text("invalid" if outcome == "malformed" else json.dumps(data), encoding="utf-8", newline="\n")
    calls = tmp_path / "calls.json"
    provider = tmp_path / "provider.py"
    provider.write_text(
        "import json, os, pathlib, sys\n"
        'pathlib.Path(os.environ["FIXTURE_CALLS"]).write_text(json.dumps(sys.argv[1:]), encoding="utf-8")\n'
        'sys.exit(int(os.environ["FIXTURE_STATUS"]))\n',
        encoding="utf-8",
        newline="\n",
    )
    script = tmp_path / "dispatch.sh"
    script.write_text(
        f'gh() {{ {shlex.quote(Path(sys.executable).as_posix())} {shlex.quote(provider.as_posix())} "$@"; }}\n'
        'just() { [[ "$*" == "--quiet _release-candidate record transfer" ]] || return 90; '
        '[[ "$FIXTURE_AUTOMATIC" == true ]] || { echo "no original Step 1 owner selection" >&2; return 2; }; printf "%s\\n" "$FIXTURE_TRANSFER"; }\n'
        + dry_run(ROOT, "release-benchmarks", [TAG]).stderr,
        encoding="utf-8",
        newline="\n",
    )
    env = {key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS", "CANDIDATE_TRANSFER_FILE"}}
    env.update(FIXTURE_CALLS=calls.as_posix(), FIXTURE_STATUS="17" if outcome == "dispatch-failure" else "0")
    env.update(FIXTURE_AUTOMATIC=str(outcome == "automatic").lower(), FIXTURE_TRANSFER=transfer.as_posix())
    if outcome not in {"missing", "automatic"}:
        env["CANDIDATE_TRANSFER_FILE"] = transfer.as_posix()
    result = run_command("bash", ["--noprofile", "--norc", "-eu", script.as_posix()], cwd=tmp_path, env=env, check=False)
    if outcome in {"success", "automatic", "dispatch-failure"}:
        assert result.returncode == (17 if outcome == "dispatch-failure" else 0), result.stderr
        assert json.loads(calls.read_text(encoding="utf-8")) == [
            "workflow",
            "run",
            "release-benchmarks.yml",
            "--repo",
            REPOSITORY,
            "--ref",
            TAG,
            "-f",
            f"tag={TAG}",
            "-F",
            f"candidate=@{transfer.as_posix()}",
        ]
    else:
        assert result.returncode != 0
        assert not calls.exists()
        if outcome == "missing":
            assert "no original Step 1 owner selection" in result.stderr


def test_benchmark_dispatch_requires_an_argument() -> None:
    result = run_command("just", ["release-benchmarks"], cwd=ROOT, check=False)
    assert result.returncode == 1
    assert "got 0 positional arguments but takes 1" in result.stderr


@pytest.mark.parametrize("tag", ["", "main", "v01.2.3", "v1.2.3; touch injected"])
def test_benchmark_dispatch_rejects_noncanonical_or_shell_tags(tmp_path: Path, tag: str) -> None:
    script = tmp_path / "invalid-tag.sh"
    script.write_text("gh() { exit 91; }\n" + dry_run(ROOT, "release-benchmarks", [tag]).stderr, encoding="utf-8", newline="\n")
    result = run_command("bash", ["--noprofile", "--norc", "-eu", script.as_posix()], cwd=tmp_path, check=False)
    assert result.returncode == 2
    assert "A canonical stable tag is required" in result.stderr
    assert not (tmp_path / "injected").exists()


def draft_provider(state: dict[str, Any], args: list[str]) -> bytes:
    """Serve release observations while retaining requests and exact notes."""
    assert args[0] == "api"
    if "POST" in args:
        request = json.loads(Path(args[-1]).read_bytes())
        state["posts"].append(request)
        state["release"] = {**request, "id": 42, "url": f"https://api.github.com/repos/{REPOSITORY}/releases/42", "immutable": False, "assets": []}
        return json.dumps(state["release"]).encode()
    endpoint = args[-1]
    if endpoint.endswith("releases?per_page=100"):
        return json.dumps([[state["release"]] if state["release"] else []]).encode()
    if endpoint.endswith("releases/42"):
        return json.dumps(state["release"]).encode()
    if endpoint.endswith("releases/42/assets?per_page=100"):
        return b"[[]]"
    if endpoint.endswith(f"git/ref/tags/{TAG}"):
        return json.dumps({"ref": f"refs/tags/{TAG}", "object": {"type": "tag", "sha": "b" * 40}}).encode()
    assert endpoint.endswith("git/tags/" + "b" * 40), args
    return json.dumps({"sha": "b" * 40, "object": {"type": "commit", "sha": state.get("remote", COMMIT)}}).encode()


@pytest.fixture
def draft_handoff(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """Keep the real evidence verifier, draft API, receipt and workflow writer."""
    provider = EvidenceProvider(tmp_path)
    provider.install(monkeypatch)
    selection = release_evidence.verify_release_evidence(POLICY, "0.5.42", COMMIT, lambda endpoint: provider.api(endpoint, repository=REPOSITORY))
    notes = b"Generated notes with `code`, $literal and CRLF\r\n"
    current = {
        "id": 77,
        "workflow_id": 9,
        "run_attempt": 1,
        "event": "push",
        "head_branch": TAG,
        "head_sha": COMMIT,
        "path": ".github/workflows/prepare-release.yml",
    }
    artifact = {
        "id": 88,
        "expired": False,
        "name": "draft-handoff-77-1",
        "digest": "sha256:" + "d" * 64,
        "workflow_run": {"id": 77, "head_sha": COMMIT},
    }
    state: dict[str, Any] = {"release": None, "posts": [], "current": current, "artifact": artifact, "history": [deepcopy(current)], "selection": selection}
    for key, value in {
        "GITHUB_REPOSITORY": REPOSITORY,
        "GITHUB_RUN_ID": "77",
        "GITHUB_RUN_ATTEMPT": "1",
        "GITHUB_EVENT_NAME": "push",
        "GITHUB_REF": f"refs/tags/{TAG}",
        "GITHUB_SHA": COMMIT,
        "RELEASE_TAG": TAG,
        "RELEASE_COMMIT": COMMIT,
        "HANDOFF_ID": "88",
        "HANDOFF_SHA256": "d" * 64,
        "VALIDATED_ATTEMPT": "1",
    }.items():
        monkeypatch.setenv(key, value)
    monkeypatch.chdir(tmp_path)

    def api(endpoint: str, *, repository: str) -> Any:
        assert repository == REPOSITORY
        if endpoint == "actions/artifacts/88":
            return deepcopy(artifact)
        if endpoint == "actions/runs/77":
            return deepcopy(current)
        if endpoint == "actions/workflows/9/runs?per_page=100&page=1":
            return {"workflow_runs": deepcopy(state["history"])}
        return provider.api(endpoint, repository=repository)

    def download(repository: str, _artifact: object, destination: Path, names: tuple[str, ...]) -> None:
        assert repository == REPOSITORY
        assert names == ("evidence.toml", "selection.json", "notes.md")
        destination.mkdir(exist_ok=True)
        (destination / "evidence.toml").write_bytes((ROOT / "tooling/release-evidence.toml").read_bytes())
        (destination / "selection.json").write_bytes(selection)
        (destination / "notes.md").write_bytes(notes)

    def gh(_root: Path, args: list[str]) -> bytes:
        return draft_provider(state, args)

    monkeypatch.setattr(release_validation, "api", api)
    monkeypatch.setattr(release_drafts, "api", api)
    monkeypatch.setattr(release_validation, "download_artifact", download)
    monkeypatch.setattr(release_assets, "_gh", gh)
    monkeypatch.setattr(release_drafts, "_gh", gh)
    # Step 2 must never launch builds, tests or acquisition, even on retries.
    monkeypatch.setattr(subprocess, "run", lambda *_args, **_kwargs: pytest.fail("draft inspection must launch no build/test/acquisition command"))
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    source = next(step["run"] for step in workflow["jobs"]["prepare"]["steps"] if "run" in step)
    script = tmp_path / "writer.py"
    script.write_text(source, encoding="utf-8", newline="\n")
    return {**state, "state": state, "provider": provider, "script": script}


@pytest.mark.parametrize("suffix", ["", "@refs/tags/" + TAG, "@" + TAG])
def test_actual_writer_creates_once_and_retains_exact_draft(draft_handoff: dict[str, Any], tmp_path: Path, suffix: str) -> None:
    draft_handoff["current"]["path"] += suffix
    script = draft_handoff["script"]
    runpy.run_path(str(script), run_name="__main__")
    original = (tmp_path / "draft-receipt.json").read_bytes()
    runpy.run_path(str(script), run_name="__main__")
    assert len(draft_handoff["posts"]) == 1
    assert (tmp_path / "draft-receipt.json").read_bytes() == original
    target = release_assets.parse_release_target((tmp_path / "target.json").read_bytes())
    assert (target.tag, target.commit, target.release_id, target.expected_title) == (TAG, COMMIT, 42, TAG)


@pytest.mark.parametrize(
    "failure", ["ref", "commit", "attempt", "attempt-spelling", "validated-attempt", "artifact", "digest", "expired", "history", "workflow", "native"]
)
def test_writer_rejects_changed_handoff_without_creating_a_draft(
    draft_handoff: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    if failure in {"attempt-spelling", "validated-attempt"}:
        monkeypatch.setenv("GITHUB_RUN_ATTEMPT" if failure == "attempt-spelling" else "VALIDATED_ATTEMPT", "01")
    elif failure in {"ref", "commit", "attempt"}:
        monkeypatch.setenv({"ref": "GITHUB_REF", "commit": "GITHUB_SHA", "attempt": "GITHUB_RUN_ATTEMPT"}[failure], "changed")
        if failure == "attempt":
            monkeypatch.setenv("GITHUB_RUN_ATTEMPT", "2")
    elif failure in {"artifact", "digest", "expired"}:
        key, value = {"artifact": ("id", 89), "digest": ("digest", "sha256:" + "e" * 64), "expired": ("expired", True)}[failure]
        draft_handoff["artifact"][key] = value
    elif failure == "history":
        draft_handoff["history"].append({**draft_handoff["current"], "id": 76})
    elif failure == "workflow":
        draft_handoff["current"]["path"] = ".github/workflows/other.yml"
    else:
        draft_handoff["provider"].jobs[0]["conclusion"] = "failure"
    diagnostic = (
        "original exact tag push"
        if failure in {"ref", "commit", "attempt", "attempt-spelling", "validated-attempt"}
        else "handoff|preparation intent|did not succeed"
    )
    with pytest.raises(ValueError, match=diagnostic):
        runpy.run_path(str(draft_handoff["script"]), run_name="__main__")
    assert draft_handoff["posts"] == []
    assert not (tmp_path / "draft-receipt.json").exists()


def test_draft_workflow_keeps_the_writer_isolated_and_recovery_explicit() -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    assert workflow[True] == {"push": {"tags": ["v*"]}}
    assert workflow["concurrency"] == {"group": "release-draft-${{ github.repository }}-${{ github.ref_name }}", "cancel-in-progress": False}
    writer = workflow["jobs"]["prepare"]
    assert writer["needs"] == "verify"
    assert writer["permissions"] == {"actions": "read", "contents": "write"}
    assert all(not step.get("uses", "").startswith("./") and "checkout@" not in step.get("uses", "") for step in writer["steps"])
    step = writer["steps"][1]
    # Custom shells use .NET argument parsing; single quotes are not delimiters.
    command = shlex.split(step["shell"], posix=False)
    assert command[:9] == ["uv", "run", "--no-config", "--no-project", "--isolated", "--no-cache", "--python", "3.15", "--with"]
    assert command[9] == f'"{TOOLING_REQUIREMENT}"'
    assert command[9].endswith('"')
    assert command[10:] == ["python", "-I", '"{0}"']
    assert step["working-directory"] == "${{ runner.temp }}"
    assert writer["steps"][-1]["if"] == "always()"
    assert "draft-receipt.json" in writer["steps"][-1]["with"]["path"]
    assert all(step.get("with", {}).get("enable-cache", False) is False for step in writer["steps"])


def test_draft_workflow_avoids_actionlint_prestart_pipe_deadlock() -> None:
    """The pinned linter writes Bash scripts before starting their pipe reader."""
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    for job in workflow["jobs"].values():
        for step in job["steps"]:
            if "run" in step and step.get("shell", "bash").split()[0] in {"bash", "sh"}:
                payload = ("set -eo pipefail\n" + step["run"] + "\n").encode("utf-8")
                assert len(payload) <= 4096, f"{step['name']} sends {len(payload)} bytes to actionlint's pre-start Windows pipe"
