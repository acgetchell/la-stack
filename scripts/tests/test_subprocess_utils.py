"""Tests for subprocess_utils.py — secure subprocess wrappers."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

import subprocess_utils
from subprocess_utils import (
    ExecutableNotFoundError,
    _build_run_kwargs,
    get_safe_executable,
    run_git_command,
    run_git_command_with_input,
)

# ---------------------------------------------------------------------------
# get_safe_executable
# ---------------------------------------------------------------------------


class TestGetSafeExecutable:
    def test_finds_git(self) -> None:
        path = get_safe_executable("git")
        assert "git" in path

    def test_raises_for_nonexistent_command(self) -> None:
        with pytest.raises(ExecutableNotFoundError, match="not found in PATH"):
            get_safe_executable("definitely_not_a_real_command_12345")


# ---------------------------------------------------------------------------
# _build_run_kwargs
# ---------------------------------------------------------------------------


class TestBuildRunKwargs:
    def test_defaults(self) -> None:
        kwargs = _build_run_kwargs("test_func")
        assert kwargs["capture_output"] is True
        assert kwargs["text"] is True
        assert kwargs["check"] is True
        assert kwargs["encoding"] == "utf-8"

    def test_rejects_shell_true(self) -> None:
        with pytest.raises(ValueError, match="shell=True is not allowed"):
            _build_run_kwargs("test_func", shell=True)  # noqa: S604

    def test_rejects_executable_override(self) -> None:
        with pytest.raises(ValueError, match="Overriding 'executable' is not allowed"):
            _build_run_kwargs("test_func", executable="/bin/sh")

    def test_strips_text_kwarg(self) -> None:
        """User-provided text=False is ignored; we always enforce text=True."""
        kwargs = _build_run_kwargs("test_func", text=False)
        assert kwargs["text"] is True

    def test_allows_check_false(self) -> None:
        kwargs = _build_run_kwargs("test_func", check=False)
        assert kwargs["check"] is False

    def test_respects_custom_encoding(self) -> None:
        kwargs = _build_run_kwargs("test_func", encoding="latin-1")
        assert kwargs["encoding"] == "latin-1"


# ---------------------------------------------------------------------------
# run_git_command
# ---------------------------------------------------------------------------


class TestRunGitCommand:
    def test_runs_simple_git_command(self) -> None:
        result = run_git_command(["rev-parse", "--git-dir"])
        assert result.returncode == 0
        assert result.stdout.strip()  # should output something like ".git"

    def test_raises_on_bad_command(self) -> None:
        with pytest.raises(subprocess_utils.subprocess.CalledProcessError):
            run_git_command(["not-a-real-git-subcommand"])


# ---------------------------------------------------------------------------
# run_git_command_with_input
# ---------------------------------------------------------------------------


class TestRunGitCommandWithInput:
    def test_passes_stdin_data(self) -> None:
        """Use git hash-object --stdin to verify input piping works."""
        result = run_git_command_with_input(
            ["hash-object", "--stdin"],
            input_data="hello\n",
        )
        # git hash-object of "hello\n" is a well-known SHA
        assert result.returncode == 0
        assert result.stdout.strip()  # should be a 40-char hex hash

    @patch("subprocess_utils.get_safe_executable", return_value="/usr/bin/git")
    @patch("subprocess_utils.subprocess.run")
    def test_input_data_forwarded(self, mock_run: MagicMock, _mock_exe: MagicMock) -> None:
        """Verify input_data is passed as the 'input' kwarg to subprocess.run."""
        run_git_command_with_input(["tag", "-a", "v1.0.0", "-F", "-"], input_data="tag body")
        mock_run.assert_called_once()
        _args, kwargs = mock_run.call_args
        assert kwargs["input"] == "tag body"
