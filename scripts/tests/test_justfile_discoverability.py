"""Regression tests for the Just recipe surface."""

import json
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any

from subprocess_utils import run_cargo_command, run_safe_command

REPO_ROOT = Path(__file__).resolve().parents[2]


def run_just(
    *args: str,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run the repository's installed Just executable without a shell."""
    executable = shutil.which("just")
    assert executable is not None
    return subprocess.run(  # noqa: S603 - executable is resolved; arguments are fixed by tests.
        [executable, *args],
        cwd=REPO_ROOT,
        check=check,
        capture_output=True,
        encoding="utf-8",
        env=env,
    )


def just_recipes() -> dict[str, dict[str, Any]]:
    """Return parsed recipe metadata from the pinned Just executable."""
    result = run_just("--dump", "--dump-format", "json")
    recipes = json.loads(result.stdout)["recipes"]
    assert isinstance(recipes, dict)
    return recipes


def test_exact_benchmark_package_excludes_peer_libraries() -> None:
    # No dependency resolution: this also runs before Cargo caches exist in CI.
    metadata = json.loads(
        run_cargo_command(
            ["metadata", "--locked", "--offline", "--no-deps", "--format-version", "1"],
            cwd=REPO_ROOT,
        ).stdout
    )
    packages = {package["name"]: package for package in metadata["packages"]}
    library = packages["la-stack"]
    comparison = packages["la-stack-comparison"]
    dependencies = {dependency["name"] for dependency in library["dependencies"]}

    assert metadata["workspace_default_members"] == [library["id"]]
    assert {"criterion", "num-rational"} <= dependencies
    assert not {"nalgebra", "faer", "la-stack-comparison"} & dependencies
    assert "exact" in {target["name"] for target in library["targets"]}
    assert "vs_linalg" not in {target["name"] for target in library["targets"]}
    assert {"nalgebra", "faer"} <= {dependency["name"] for dependency in comparison["dependencies"]}
    assert "vs_linalg" in {target["name"] for target in comparison["targets"]}


def test_release_runs_each_suite_once_and_reuses_peer_measurements() -> None:
    baseline = run_just("--dry-run", "bench-save-baseline", "v0.4.5", "all")
    # Execute the real recipe's shell branching while replacing only Cargo.
    runner = run_just("--evaluate", "_run").stdout.strip()
    script = 'cargo() { printf "%s\\n" "$*"; }\n' + baseline.stderr.replace(runner + " ", "")
    baseline_run = run_safe_command("bash", ["--noprofile", "--norc", "-euc", script], cwd=REPO_ROOT)
    baseline_commands = [shlex.split(line) for line in baseline_run.stdout.splitlines()]
    recipes = just_recipes()
    current = recipes["bench-latest"]
    assert not current["body"]
    current_commands = []
    for dependency in current["dependencies"]:
        recipe = recipes[dependency["recipe"]]
        assert not recipe["dependencies"]
        commands = run_just("--dry-run", dependency["recipe"]).stderr
        current_commands.extend(shlex.split(line.removeprefix(runner + " ")) for line in commands.splitlines())

    assert len(baseline_commands) == len(current_commands) == 2
    baseline_comparison, baseline_exact = baseline_commands
    current_comparison, current_exact = current_commands
    for command in (baseline_comparison, current_comparison):
        assert command[command.index("-p") + 1] == "la-stack-comparison"
        assert command[command.index("--bench") + 1] == "vs_linalg"
    for command in (baseline_exact, current_exact):
        assert "-p" not in command
        assert command[command.index("--bench") + 1] == "exact"
    assert baseline_comparison[baseline_comparison.index("--") + 1 :] == ["--noplot", "--save-baseline", "v0.4.5"]
    assert current_comparison[current_comparison.index("--") + 1 :] == ["la_stack", "--noplot"]


def test_ci_enforces_full_python_fixture_lint_policy() -> None:
    """Canonical CI should lint fixtures without narrowing the Ruff configuration."""
    recipes = just_recipes()
    ci_dependencies = {dependency["recipe"] for dependency in recipes["ci"]["dependencies"]}
    fixture_lint_body = json.dumps(recipes["python-fixture-lint"]["body"])

    assert "python-fixture-lint" in ci_dependencies
    assert "ruff check tests/semgrep/scripts/" in fixture_lint_body
    assert "--select" not in fixture_lint_body
