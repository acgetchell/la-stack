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
from pathlib import Path
from typing import Any


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


def _build_run_kwargs(function_name: str, **kwargs: Any) -> dict[str, Any]:
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


def run_git_command_with_input(
    args: list[str],
    input_data: str,
    cwd: Path | None = None,
    **kwargs: Any,
) -> subprocess.CompletedProcess[str]:
    """Run a git command securely with stdin input using full executable path.

    Args:
        args: Git command arguments (without 'git' prefix)
        input_data: Data to send to stdin
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
    return subprocess.run(  # noqa: S603,PLW1510
        [git_path, *args],
        cwd=cwd,
        input=input_data,
        **run_kwargs,
    )
