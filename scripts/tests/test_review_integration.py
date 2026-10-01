"""Consumer integration with the published review CLI; no live CodeRabbit calls."""

import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest
from research_repo_tools.cli import main

from benchmark_process import run_git_command, run_safe_command

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def consumer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Import actual recipes and configuration into isolated Git history."""
    root = tmp_path / "consumer with spaces"
    root.mkdir()
    for name in ("pyproject.toml", "uv.lock", "AGENTS.md", ".coderabbit.yaml"):
        shutil.copyfile(REPO_ROOT / name, root / name)
    (root / "recipes.just").write_text(f"import '{(REPO_ROOT / 'justfile').as_posix()}'\n", encoding="utf-8")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    run_git_command(["init", "--quiet", "--initial-branch=main"], cwd=root)
    run_git_command(
        ["-c", "user.name=Review fixture", "-c", "user.email=fixture@example.invalid", "commit", "--quiet", "--allow-empty", "-m", "fixture"], cwd=root
    )
    head = run_git_command(["rev-parse", "HEAD"], cwd=root).stdout.strip()
    run_git_command(["update-ref", "refs/remotes/origin/main", head], cwd=root)
    run_git_command(["remote", "add", "origin", str(root)], cwd=root)
    return root


@pytest.fixture
def stub(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Shadow the service executable on all supported platforms."""
    directory = tmp_path / "bin"
    directory.mkdir()
    script = directory / "stub.py"
    script.write_text(
        "import json, os, sys\n"
        "print(json.dumps(sys.argv[1:]))\n"
        "print('fixture diagnostic', file=sys.stderr)\n"
        "sys.exit(int(os.environ.get('REVIEW_STUB_STATUS', '0')))\n",
        encoding="utf-8",
    )
    if os.name == "nt":
        launcher = directory / "coderabbit.cmd"
        launcher.write_text(f"@echo off\n{subprocess.list2cmdline([sys.executable, str(script)])} %*\n", encoding="utf-8")
    else:
        launcher = directory / "coderabbit"
        launcher.write_text(f"#!{sys.executable}\n{script.read_text(encoding='utf-8')}", encoding="utf-8")
        launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{directory}{os.pathsep}{os.environ['PATH']}")
    return launcher


def recipe(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    """Execute the imported adapter with the installed locked package."""
    return run_safe_command(
        "just",
        ["--justfile", str(root / "recipes.just"), "--working-directory", str(root), *args],
        cwd=root,
        env=os.environ | {"UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
        check=False,
    )


@pytest.mark.parametrize(
    ("args", "scope"),
    [(("review",), "--base=origin/main"), (("review", "HEAD"), "--base=HEAD"), (("review-uncommitted",), "--uncommitted")],
)
def test_actual_recipes_forward_scopes_and_consumer_instructions(consumer: Path, stub: Path, args: tuple[str, ...], scope: str) -> None:
    """Both scopes use structured output, untracked inputs, and consumer instructions."""
    before = run_git_command(["show-ref"], cwd=consumer).stdout
    result = recipe(consumer, *args)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [
        "review",
        "--agent",
        "--include-untracked",
        scope,
        "--config",
        str(consumer / "AGENTS.md"),
        str(consumer / ".coderabbit.yaml"),
    ]
    assert "fixture diagnostic" in result.stderr
    assert run_git_command(["show-ref"], cwd=consumer).stdout == before


@pytest.mark.parametrize("args", [("review", "HEAD"), ("review-uncommitted",)])
def test_local_scopes_do_not_require_a_remote(consumer: Path, stub: Path, args: tuple[str, ...]) -> None:
    run_git_command(["remote", "remove", "origin"], cwd=consumer)
    assert recipe(consumer, *args).returncode == 0


@pytest.mark.parametrize("status", [7, 130])
def test_service_failures_and_interruption_status_reach_just(consumer: Path, stub: Path, monkeypatch: pytest.MonkeyPatch, status: int) -> None:
    monkeypatch.setenv("REVIEW_STUB_STATUS", str(status))
    result = recipe(consumer, "review-uncommitted")
    assert result.returncode == status
    assert f"exit code {status}" in result.stderr
    assert "fixture diagnostic" in result.stderr


@pytest.mark.parametrize("problem", ["missing", "stale", "remote"])
def test_default_base_fails_closed_before_service_invocation(consumer: Path, stub: Path, problem: str) -> None:
    if problem == "missing":
        run_git_command(["update-ref", "-d", "refs/remotes/origin/main"], cwd=consumer)
    elif problem == "stale":
        run_git_command(
            ["-c", "user.name=Review fixture", "-c", "user.email=fixture@example.invalid", "commit", "--quiet", "--allow-empty", "-m", "advance"], cwd=consumer
        )
    else:
        run_git_command(["remote", "set-url", "origin", str(consumer / "missing-remote")], cwd=consumer)
    result = recipe(consumer, "review")
    assert result.returncode != 0
    assert result.stdout == ""
    assert "fixture diagnostic" not in result.stderr
    assert ("Cannot verify origin/main" if problem == "remote" else "git fetch origin") in result.stderr


def test_shell_metacharacters_are_passed_as_one_rejected_base(consumer: Path, stub: Path) -> None:
    result = recipe(consumer, "review", "HEAD'; echo INJECTION_EXECUTED; #")
    assert result.returncode != 0
    assert "review base must be" in result.stderr
    assert result.stdout == ""


@pytest.mark.parametrize("name", ["AGENTS.md", ".coderabbit.yaml"])
def test_missing_consumer_instructions_stop_review(consumer: Path, stub: Path, name: str) -> None:
    (consumer / name).unlink()
    result = recipe(consumer, "review-uncommitted")
    assert result.returncode != 0
    assert "review requires" in result.stderr
    assert result.stdout == ""


def test_missing_cli_reports_explicit_installation(consumer: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    original = shutil.which
    monkeypatch.setattr(shutil, "which", lambda command, *args, **kwargs: None if Path(command).name == "coderabbit" else original(command, *args, **kwargs))
    assert main(["--root", str(consumer), "review", "uncommitted"]) == 1
    assert "Install and authenticate it explicitly" in capsys.readouterr().err


def test_reviews_are_discoverable_and_outside_routine_gates() -> None:
    surface = json.loads(run_safe_command("just", ["--dump", "--dump-format", "json"], cwd=REPO_ROOT).stdout)["recipes"]
    help_text = run_safe_command("just", [], cwd=REPO_ROOT).stdout
    listed = run_safe_command("just", ["--list"], cwd=REPO_ROOT).stdout
    for name in ("review", "review-uncommitted"):
        assert name in help_text
        assert name in listed
        assert not surface[name]["dependencies"]
    for gate in ("check", "ci", "setup", "update"):
        pending = [gate]
        seen: set[str] = set()
        while pending:
            name = pending.pop()
            if name in seen:
                continue
            seen.add(name)
            pending.extend(dependency["recipe"] for dependency in surface[name]["dependencies"])
        assert not {"review", "review-uncommitted"} & seen
