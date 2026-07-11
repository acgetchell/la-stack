#!/usr/bin/env python3
"""Secure subprocess utilities for Python scripts.

This module provides secure subprocess wrappers that:
- Use full executable paths instead of command names
- Validate executables exist before running
- Provide consistent error handling
- Mitigate security vulnerabilities flagged by Bandit

All scripts should use these functions instead of calling subprocess directly.

Ported from the delaunay project's scripts/subprocess_utils.py (minimal subset).
"""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

type RunKwargs = dict[str, Any]


class ExecutableNotFoundError(Exception):
    """Raised when a required executable is not found in PATH."""


def get_safe_executable(command: str) -> str:
    """Get the full path to an executable, validating it exists.

    Args:
        command: Command name to find (e.g., "git")

    Returns:
        Full path to the executable

    Raises:
        ExecutableNotFoundError: If executable is not found in PATH
    """
    full_path = shutil.which(command)
    if full_path is None:
        raise ExecutableNotFoundError(f"Required executable '{command}' not found in PATH")
    return full_path


def _build_run_kwargs(function_name: str, **kwargs: Any) -> RunKwargs:
    """Build secure kwargs for subprocess.run with consistent hardening.

    Args:
        function_name: Name of the calling function (for error messages)
        **kwargs: User-provided kwargs to validate and merge

    Returns:
        Validated and hardened kwargs dict for subprocess.run

    Raises:
        ValueError: If insecure parameters are provided
    """
    # Disallow shell=True to preserve security guarantees
    if kwargs.get("shell"):
        msg = f"shell=True is not allowed in {function_name}"
        raise ValueError(msg)
    # Disallow overriding the program to execute
    if "executable" in kwargs:
        msg = f"Overriding 'executable' is not allowed in {function_name}"
        raise ValueError(msg)
    # Enforce text mode for stable typing (CompletedProcess[str])
    kwargs.pop("text", None)
    run_kwargs = {
        "capture_output": True,
        "text": True,
        "check": True,  # Secure default
        **kwargs,  # Allow overriding other safe defaults
    }
    # Prefer deterministic UTF-8 unless caller overrides
    run_kwargs.setdefault("encoding", "utf-8")
    return run_kwargs


def run_git_command(
    args: list[str],
    cwd: Path | None = None,
    **kwargs: Any,
) -> subprocess.CompletedProcess[str]:
    """Run a git command securely using full executable path.

    Args:
        args: Git command arguments (without 'git' prefix)
        cwd: Working directory for the command
        **kwargs: Additional arguments passed to subprocess.run

    Returns:
        CompletedProcess result

    Raises:
        ExecutableNotFoundError: If git is not found
        subprocess.CalledProcessError: If command fails and check=True
        subprocess.TimeoutExpired: If command times out
    """
    git_path = get_safe_executable("git")
    run_kwargs = _build_run_kwargs("run_git_command", **kwargs)
    return subprocess.run(  # noqa: S603,PLW1510
        [git_path, *args],
        cwd=cwd,
        **run_kwargs,
    )


def run_cargo_command(
    args: list[str],
    cwd: Path | None = None,
    **kwargs: Any,
) -> subprocess.CompletedProcess[str]:
    """Run a cargo command securely using full executable path.

    Args:
        args: Cargo command arguments (without 'cargo' prefix)
        cwd: Working directory for the command
        **kwargs: Additional arguments passed to subprocess.run

    Returns:
        CompletedProcess result

    Raises:
        ExecutableNotFoundError: If cargo is not found
        subprocess.CalledProcessError: If command fails and check=True
        subprocess.TimeoutExpired: If command times out
    """
    cargo_path = get_safe_executable("cargo")
    run_kwargs = _build_run_kwargs("run_cargo_command", **kwargs)
    return subprocess.run(  # noqa: S603,PLW1510
        [cargo_path, *args],
        cwd=cwd,
        **run_kwargs,
    )


def run_safe_command(
    command: str,
    args: list[str],
    cwd: Path | None = None,
    **kwargs: Any,
) -> subprocess.CompletedProcess[str]:
    """Run any command securely using full executable path.

    Args:
        command: Command name to run (for example, "rustc" or "gnuplot")
        args: Command arguments
        cwd: Working directory for the command
        **kwargs: Additional arguments passed to subprocess.run

    Returns:
        CompletedProcess result

    Raises:
        ExecutableNotFoundError: If command is not found
        subprocess.CalledProcessError: If command fails and check=True
        subprocess.TimeoutExpired: If command times out
    """
    command_path = get_safe_executable(command)
    run_kwargs = _build_run_kwargs(f"run_safe_command for {command}", **kwargs)
    return subprocess.run(  # noqa: S603,PLW1510
        [command_path, *args],
        cwd=cwd,
        **run_kwargs,
    )


def get_git_commit_hash(cwd: Path | None = None) -> str:
    """Get the current git commit hash."""
    result = run_git_command(["rev-parse", "HEAD"], cwd=cwd)
    return result.stdout.strip()


def get_git_remote_url(remote: str = "origin", cwd: Path | None = None) -> str:
    """Get the URL of a git remote."""
    result = run_git_command(["remote", "get-url", remote], cwd=cwd)
    return result.stdout.strip()


def check_git_repo() -> bool:
    """Return true when the current directory is inside a git repository."""
    try:
        run_git_command(["rev-parse", "--git-dir"])
    except ExecutableNotFoundError, subprocess.CalledProcessError:
        return False
    else:
        return True


def check_git_history() -> bool:
    """Return true when the current git repository has at least one commit."""
    try:
        run_git_command(["log", "--oneline", "-n", "1"])
    except ExecutableNotFoundError, subprocess.CalledProcessError:
        return False
    else:
        return True


def run_git_command_with_input(
    args: list[str],
    input_data: str,
    cwd: Path | None = None,
    **kwargs: Any,
) -> subprocess.CompletedProcess[str]:
    """Run a git command securely with stdin input using full executable path.

    Args:
        args: Git command arguments (without 'git' prefix)
        input_data: Text to encode and send to stdin without newline translation
        cwd: Working directory for the command
        **kwargs: Additional arguments passed to subprocess.run

    Returns:
        CompletedProcess result

    Raises:
        ExecutableNotFoundError: If git is not found
        subprocess.CalledProcessError: If command fails and check=True
        subprocess.TimeoutExpired: If command times out
    """
    git_path = get_safe_executable("git")
    run_kwargs = _build_run_kwargs("run_git_command_with_input", **kwargs)
    encoding: str = run_kwargs.get("encoding") or "utf-8"
    errors: str = run_kwargs.get("errors") or "strict"
    with tempfile.TemporaryFile() as stdin:
        stdin.write(input_data.encode(encoding, errors))
        stdin.seek(0)
        return subprocess.run(  # noqa: S603,PLW1510
            [git_path, *args],
            cwd=cwd,
            stdin=stdin,
            **run_kwargs,
        )


class ProjectRootNotFoundError(Exception):
    """Raised when project root directory cannot be located."""


def find_project_root() -> Path:
    """Find the nearest project root by walking upward to Cargo.toml."""
    current_dir = Path.cwd()
    project_root = current_dir
    while project_root != project_root.parent:
        if (project_root / "Cargo.toml").exists():
            return project_root
        project_root = project_root.parent
    msg = "Could not locate Cargo.toml to determine project root"
    raise ProjectRootNotFoundError(msg)
