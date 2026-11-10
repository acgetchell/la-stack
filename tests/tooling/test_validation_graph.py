"""Preserve every native validator while overlapping only independent lanes."""

import json
import os
import shlex
import sys
from email.message import Message
from io import BytesIO
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.error import HTTPError

import pytest
import yaml
from research_repo_tools import prebuilt_tools
from research_repo_tools.just_inspect import inspect_justfile
from research_repo_tools.process import run_command

if TYPE_CHECKING:
    from urllib.request import Request

ROOT = Path(__file__).resolve().parents[2]
POLICY = {
    "action-lint",
    "workflow-allowlist-check",
    "zizmor",
    "markdown-check",
    "spell-check",
    "_version-check",
    "changelog-check",
    "toml-parse-check",
    "toml-fmt-check",
    "toml-lint",
    "yaml-fmt-check",
    "yaml-lint",
    "citation-check",
    "validate-json",
    "justfile-fmt-check",
    "python-check",
    "cargo-lock-check",
    "fmt-check",
    "semgrep",
    "semgrep-test",
    "unused-deps",
    "shell-check",
}
RUST = {"clippy-all-targets", "doc-check", "test-rust-ci", "test-doc", "test-doc-exact", "bench-compile", "test-bench-exact", "examples"}


def reachable(recipes: dict[str, dict[str, Any]], name: str) -> set[str]:
    """Find prerequisite owners from native Just metadata."""
    return {name}.union(*(reachable(recipes, dependency["recipe"]) for dependency in recipes[name]["dependencies"]))


@pytest.fixture
def graph(tmp_path: Path) -> tuple[Path, set[str]]:
    """Keep the real dependency graph, replacing only validator command bodies."""
    recipes = inspect_justfile(ROOT).recipes
    assert {item["recipe"] for item in recipes["_ci-policy"]["dependencies"]} == POLICY
    assert {item["recipe"] for item in recipes["_ci-rust"]["dependencies"]} == RUST
    assert recipes["_ci-execution"]["attributes"] == ["parallel", "private"]
    assert [item["recipe"] for item in recipes["_ci-execution"]["dependencies"]] == ["_ci-rust", "test-python"]
    selected = reachable(recipes, "ci") - {"ci", "_ci-policy", "_ci-rust", "_ci-execution"}
    marker = tmp_path / "step.py"
    marker.write_text(
        "import json, os, pathlib, sys, time\n"
        "name = sys.argv[1]\nroot = pathlib.Path.cwd()\n"
        'done = root / (name + ".done")\n'
        'assert not done.exists(), "validator ran twice: " + name\n'
        'if name in {"clippy-all-targets", "test-python"}:\n'
        '    assert (root / "python-sync.done").exists(), "environment is not synchronized"\n'
        '    (root / (name + ".start")).touch()\n'
        '    if os.environ["FIXTURE_PARALLEL"] == "true":\n'
        '        peer = "test-python" if name == "clippy-all-targets" else "clippy-all-targets"\n'
        "        deadline = time.monotonic() + 10\n"
        '        while not (root / (peer + ".start")).exists():\n'
        '            assert time.monotonic() < deadline, "independent lane never started"\n'
        "            time.sleep(0.01)\n"
        'if os.environ["FIXTURE_FAIL"] == name:\n'
        '    print("injected validator failure: " + name, file=sys.stderr)\n    sys.exit(19)\n'
        'done.write_text(json.dumps({"name": name}), encoding="utf-8")\n',
        encoding="utf-8",
        newline="\n",
    )
    overrides = []
    for name in sorted(selected):
        recipe = recipes[name]
        parameters = " *args" if name in {"bench-exact", "zizmor", "test-python"} else ""
        assert not recipe["parameters"] or name in {"bench-exact", "zizmor", "test-python"}
        dependencies = " ".join(
            "(" + item["recipe"] + " " + " ".join(json.dumps(arg) for arg in item["arguments"]) + ")" if item["arguments"] else item["recipe"]
            for item in recipe["dependencies"]
        )
        overrides.append(f"{name}{parameters}: {dependencies}\n    {shlex.quote(Path(sys.executable).as_posix())} {shlex.quote(marker.as_posix())} {name}\n")
    source = (ROOT / "justfile").read_text(encoding="utf-8")
    (tmp_path / "justfile").write_text("set allow-duplicate-recipes\n" + source + "\n" + "\n".join(overrides), encoding="utf-8", newline="\n")
    return tmp_path, selected


@pytest.mark.parametrize("parallel", [False, True], ids=["serial", "parallel"])
@pytest.mark.parametrize("failure", ["none", "python-check", "test-python", "test-rust-ci"])
def test_real_graph_runs_each_validator_once_and_propagates_failures(graph: tuple[Path, set[str]], parallel: bool, failure: str) -> None:
    root, selected = graph
    env = {**os.environ, "FIXTURE_PARALLEL": str(parallel).lower(), "FIXTURE_FAIL": failure}
    result = run_command("just", ["--jobs", "2" if parallel else "1", "ci"], cwd=root, env=env, check=False, timeout=45)
    completed = {path.name.removesuffix(".done") for path in root.glob("*.done")}
    if failure == "none":
        assert result.returncode == 0, result.stderr
        assert completed == selected
    else:
        assert result.returncode != 0
        assert f"injected validator failure: {failure}" in result.stderr
        assert failure not in completed
        if failure == "python-check":
            assert not list(root.glob("*.start"))
    assert "python-sync" in completed


def test_managed_cache_contains_compilation_only_and_release_jobs_disable_it() -> None:
    action = yaml.safe_load((ROOT / ".github/actions/setup-tools/action.yml").read_bytes())
    steps = action["runs"]["steps"]
    cache = next(step for step in steps if step.get("uses", "").startswith("actions/cache@"))
    assert cache["if"] == "inputs.cache == 'true'"
    assert cache["with"]["path"] == "${{ runner.temp }}/rrt-cargo-build"
    assert "runner.os" in cache["with"]["key"]
    assert "runner.arch" in cache["with"]["key"]
    assert all(name in cache["with"]["key"] for name in ("rust-toolchain.toml", "pyproject.toml", "uv.lock", "Cargo.lock"))
    assert "restore-keys" not in cache["with"]
    setup = next(step for step in steps if step.get("run", "").endswith("research-repo-tools setup"))
    assert setup["env"] == {
        "GH_TOKEN": "${{ github.token }}",
        "RESEARCH_REPO_TOOLS_HOME": "${{ runner.temp }}/rrt",
        "RESEARCH_REPO_TOOLS_CARGO_BUILD_HOME": cache["with"]["path"],
    }
    for filename in ("ci.yml", "release-benchmarks.yml", "publish-crates.yml", "prepare-release.yml"):
        workflow = yaml.safe_load((ROOT / ".github/workflows" / filename).read_bytes())
        for name, job in workflow["jobs"].items():
            if filename == "ci.yml" and name in {"prepare", "native", "build"}:
                continue
            for step in job["steps"]:
                if step.get("uses") == "./.github/actions/setup-tools":
                    assert str(step["with"]["cache"]).lower() == "false"
                assert not step.get("uses", "").startswith("actions/cache@")


def test_hosted_python_jobs_partition_the_complete_collection() -> None:
    """Exercise pytest's actual marker selection without executing the slow suite."""
    recipes = inspect_justfile(ROOT).recipes
    assert "_ci-policy" in {item["recipe"] for item in recipes["ci-policy"]["dependencies"]}
    assert "_ci-rust" in {item["recipe"] for item in recipes["ci-rust"]["dependencies"]}
    selections: dict[str, set[str]] = {}
    for name in ("test-python", "ci-python", "ci-candidate"):
        arguments = [] if name == "test-python" else recipes[name]["dependencies"][0]["arguments"]
        if name != "test-python":
            assert [item["recipe"] for item in recipes[name]["dependencies"]] == ["test-python"]
        result = run_command(
            "uv",
            ["run", "--locked", "--no-sync", "pytest", "--collect-only", "-q", "--color=no", "-n", "0", *arguments],
            cwd=ROOT,
            timeout=90,
        )
        selections[name] = {line for line in result.stdout.splitlines() if line.startswith("tests/") and "::" in line}
        assert selections[name], result.stdout
    assert selections["ci-candidate"] == {"tests/tooling/test_performance_workflow.py::test_local_candidate_capture_and_final_acceptance_reuse_original_bytes"}
    assert selections["ci-python"].isdisjoint(selections["ci-candidate"])
    assert selections["ci-python"] | selections["ci-candidate"] == selections["test-python"]


@pytest.mark.parametrize("authenticated", [False, True], ids=["rate-limited-anonymous", "authenticated-setup"])
def test_setup_token_authenticates_installed_binary_release_lookup(monkeypatch: pytest.MonkeyPatch, authenticated: bool) -> None:
    """Exercise the real client's credential boundary using the action's mapping."""
    action = yaml.safe_load((ROOT / ".github/actions/setup-tools/action.yml").read_bytes())
    setup = next(step for step in action["runs"]["steps"] if step.get("run", "").endswith("research-repo-tools setup"))
    credential = next(name for name, value in setup["env"].items() if value == "${{ github.token }}")
    for name in ("GH_TOKEN", "GITHUB_TOKEN"):
        monkeypatch.delenv(name, raising=False)
    if authenticated:
        monkeypatch.setenv(credential, "fixture-lookup-token")
    observed: list[str | None] = []

    def lookup(request: Request, *, timeout: float) -> BytesIO:
        assert request.full_url == "https://api.github.com/repos/gitleaks/gitleaks/releases/tags/v8.30.1"
        assert timeout == 30
        authorization = request.get_header("Authorization")
        observed.append(authorization)
        if authorization != "Bearer fixture-lookup-token":
            raise HTTPError(request.full_url, 403, "rate limit exceeded", Message(), None)
        return BytesIO(json.dumps({"draft": False, "prerelease": False, "tag_name": "v8.30.1"}).encode())

    monkeypatch.setattr(prebuilt_tools.urllib.request, "urlopen", lookup)
    if authenticated:
        assert prebuilt_tools.release_metadata("gitleaks", "v8.30.1")["tag_name"] == "v8.30.1"
        assert observed == ["Bearer fixture-lookup-token"]
    else:
        with pytest.raises(RuntimeError, match="HTTP Error 403: rate limit exceeded"):
            prebuilt_tools.release_metadata("gitleaks", "v8.30.1")
        assert observed == [None]
