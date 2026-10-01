"""Benchmark command adapters around the published shared process API.

The shared package owns discovery, execution, exact byte transport, CPU metadata,
and diagnostics. These adapters retain existing benchmark call signatures and
the distinction between captured commands and live measurement progress.
"""

import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict, Unpack

from research_repo_tools.process import run_command, run_command_live

if TYPE_CHECKING:
    from collections.abc import Mapping


class CommandOptions(TypedDict, total=False):
    """Options needed by captured plotting and live measurement phases."""

    env: Mapping[str, str] | None
    check: bool
    timeout: float | None
    capture_output: bool
    input: str


def run_safe_command(
    command: str,
    args: list[str],
    cwd: Path | None = None,
    **options: Unpack[CommandOptions],
) -> subprocess.CompletedProcess[str]:
    """Select the shared captured or live runner for a benchmark phase."""
    env = options.get("env")
    check = options.get("check", True)
    timeout = options.get("timeout", 300.0)
    if options.get("capture_output", True):
        return run_command(command, args, cwd=cwd, env=env, check=check, timeout=timeout, input=options.get("input"))
    if "input" in options:
        msg = "Live benchmark commands require inherited stdin"
        raise ValueError(msg)
    result = run_command_live(command, args, cwd=cwd, env=env, check=check, timeout=timeout)
    return subprocess.CompletedProcess(result.args, result.returncode, "", "")


def run_git_command(
    args: list[str],
    cwd: Path | None = None,
    *,
    check: bool = True,
    timeout: float | None = 300.0,
    env: Mapping[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Capture benchmark Git metadata without opening a pager."""
    return run_command("git", ["--no-pager", *args], cwd=cwd, check=check, timeout=timeout, env=env)


def find_project_root(start: Path | None = None) -> Path:
    """Find the benchmark consumer's nearest Cargo manifest from its working directory."""
    directory = (start or Path.cwd()).resolve()
    if directory.is_file():
        directory = directory.parent
    for candidate in (directory, *directory.parents):
        if (candidate / "Cargo.toml").is_file():
            return candidate
    msg = "Could not locate Cargo.toml to determine benchmark root"
    raise FileNotFoundError(msg)
