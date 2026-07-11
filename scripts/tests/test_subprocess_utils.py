"""Tests for subprocess_utils.py — secure subprocess wrappers."""

from __future__ import annotations

from typing import BinaryIO, cast
from unittest.mock import MagicMock, patch

import pytest

import subprocess_utils
from subprocess_utils import (
    ExecutableNotFoundError,
    _build_run_kwargs,
    check_git_history,
    check_git_repo,
    find_project_root,
    get_git_commit_hash,
    get_git_remote_url,
    get_safe_executable,
    run_cargo_command,
    run_git_command,
    run_git_command_with_input,
    run_safe_command,
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
        assert result.stdout.strip() == "ce013625030ba8dba906f756967f9e9ca394464a"

    def test_input_data_forwarded_as_raw_utf8(self) -> None:
        """Verify stdin bytes preserve LF even when subprocess output is text."""
        observed_input = b""

        def capture_run(*_args: object, **kwargs: object) -> subprocess_utils.subprocess.CompletedProcess[str]:
            nonlocal observed_input
            stdin = cast("BinaryIO", kwargs["stdin"])
            observed_input = stdin.read()
            return subprocess_utils.subprocess.CompletedProcess(args=["git"], returncode=0, stdout="", stderr="")

        with (
            patch("subprocess_utils.get_safe_executable", return_value="/usr/bin/git") as mock_executable,
            patch("subprocess_utils.subprocess.run", side_effect=capture_run) as mock_run,
        ):
            run_git_command_with_input(["tag", "-a", "v1.0.0", "-F", "-"], input_data="tag body\n")
        mock_executable.assert_called_once_with("git")
        mock_run.assert_called_once()
        _args, kwargs = mock_run.call_args
        assert observed_input == b"tag body\n"
        assert "input" not in kwargs
        assert kwargs["text"] is True


class TestAdditionalHelpers:
    def test_run_cargo_command_uses_safe_executable(self) -> None:
        with (
            patch("subprocess_utils.get_safe_executable", return_value="/usr/bin/cargo") as mock_executable,
            patch("subprocess_utils.subprocess.run") as mock_run,
        ):
            run_cargo_command(["--version"])
        mock_executable.assert_called_once_with("cargo")
        mock_run.assert_called_once()
        args, _kwargs = mock_run.call_args
        assert args[0] == ["/usr/bin/cargo", "--version"]

    def test_run_safe_command_uses_safe_executable(self) -> None:
        with (
            patch("subprocess_utils.get_safe_executable", return_value="/usr/bin/gnuplot") as mock_executable,
            patch("subprocess_utils.subprocess.run") as mock_run,
        ):
            run_safe_command("gnuplot", ["--version"])
        mock_executable.assert_called_once_with("gnuplot")
        mock_run.assert_called_once()
        args, _kwargs = mock_run.call_args
        assert args[0] == ["/usr/bin/gnuplot", "--version"]

    @patch("subprocess_utils.run_git_command")
    def test_git_convenience_helpers(self, mock_run_git: MagicMock) -> None:
        def fake_run_git(args: list[str], **_kwargs: object) -> subprocess_utils.subprocess.CompletedProcess[str]:
            stdout_by_args = {
                ("rev-parse", "HEAD"): "abc123def456\n",
                ("remote", "get-url", "origin"): "https://github.com/example/repo.git\n",
                ("rev-parse", "--git-dir"): ".git\n",
                ("log", "--oneline", "-n", "1"): "abc123d message\n",
            }
            return subprocess_utils.subprocess.CompletedProcess(args=["git", *args], returncode=0, stdout=stdout_by_args[tuple(args)])

        mock_run_git.side_effect = fake_run_git

        assert get_git_commit_hash() == "abc123def456"
        assert get_git_remote_url() == "https://github.com/example/repo.git"
        assert check_git_repo() is True
        assert check_git_history() is True
        assert [call_args.args[0] for call_args in mock_run_git.call_args_list] == [
            ["rev-parse", "HEAD"],
            ["remote", "get-url", "origin"],
            ["rev-parse", "--git-dir"],
            ["log", "--oneline", "-n", "1"],
        ]

    def test_find_project_root(self) -> None:
        assert (find_project_root() / "Cargo.toml").is_file()
