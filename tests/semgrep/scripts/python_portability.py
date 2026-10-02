"""Static-analysis fixtures only; never execute these example calls."""

import io
import os
import subprocess
import tempfile
from subprocess import run as run_process
from tempfile import NamedTemporaryFile as named_file
from typing import TYPE_CHECKING

from research_repo_tools import process as shared_process
from research_repo_tools.process import run_command, run_git_bytes as git_input

import benchmark_process as utils
from benchmark_process import run_git_command, run_safe_command

if TYPE_CHECKING:
    from pathlib import Path


def git_input_routing(payload: str, argv: list[str]) -> None:
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    run_safe_command("git", argv, input=payload)
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    utils.run_safe_command("git", argv, input=payload)
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    run_command("git", argv, input=payload)
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    shared_process.run_command("git", argv, input=payload)
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    subprocess.run(["git", "hash-object", "--stdin"], input=payload, text=True)
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    run_process(["git", "hash-object", "--stdin"], input=payload, text=True)
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    process = subprocess.Popen(["git", "hash-object", "--stdin"], stdin=subprocess.PIPE, text=True)
    process.communicate(payload)
    # Even binary Git input must go through the shared helper.
    # ruleid: la-stack.python.git-stdin-use-shared-helper
    process = subprocess.Popen(["git", "hash-object", "--stdin"], stdin=subprocess.PIPE, text=False)
    process.communicate(payload.encode("utf-8"))

    # ok: la-stack.python.git-stdin-use-shared-helper
    git_input(argv, input=payload.encode("utf-8"))
    # ok: la-stack.python.git-stdin-use-shared-helper
    shared_process.run_git_bytes(argv, input=payload.encode("utf-8"))
    # ok: la-stack.python.git-stdin-use-shared-helper
    run_git_command(argv)
    # ok: la-stack.python.git-stdin-use-shared-helper
    shared_process.run_command("git", argv, input=None)
    # Ordinary text input to non-Git commands is outside this rule's scope.
    # ok: la-stack.python.git-stdin-use-shared-helper
    run_safe_command("ruff", ["check", "-"], input=payload)
    # ok: la-stack.python.git-stdin-use-shared-helper
    run_command("ruff", ["check", "-"], input=payload)
    # ok: la-stack.python.git-stdin-use-shared-helper
    subprocess.run(["formatter"], input=payload, text=True, encoding="utf-8")
    # ok: la-stack.python.git-stdin-use-shared-helper
    process = subprocess.Popen(["formatter"], stdin=subprocess.PIPE, text=True, encoding="utf-8")
    process.communicate(payload)
    # Git commands without piped stdin do not need the input helper.
    # ok: la-stack.python.git-stdin-use-shared-helper
    subprocess.Popen(["git", "status"], stdout=subprocess.PIPE, text=True)


def path_text_writes(path: Path, contents: str, newline_policy: str) -> None:
    # The real extracted-notebook writer, before and after explicit LF output.
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, encoding=None, newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, encoding="utf-8", newline=None)
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, encoding="utf-8", newline="")
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(contents, encoding="utf-8", newline=newline_policy)
    # Nested reads and text transformations are not file writers.
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(contents.strip(), encoding="utf-8", newline="\n")
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(path.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(open(path, encoding="utf-8").read(), encoding="utf-8", newline="\n")
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(path.read_text(encoding=None), encoding="utf-8", newline="\n")
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_text(open(path, encoding=None, newline=None).read(), encoding="utf-8", newline="\n")


def temporary_text_writes(path: Path) -> None:
    # Production writers use both positional and keyword mode forms.
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False)
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False)
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        delete=False,
    )
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        delete=False,
    )
    # ruleid: la-stack.python.text-writes-explicit-policy
    named_file(mode="w", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    named_file(mode="w", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    named_file(mode="w", encoding="utf-8", newline=None)
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.TemporaryFile("w+", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.TemporaryFile("w+", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.TemporaryFile(mode="w+", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.TemporaryFile(mode="w+", encoding="utf-8", newline="\n")


def spooled_text_writes() -> None:
    # SpooledTemporaryFile takes max_size before its positional mode argument.
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(1024, "w+", encoding="utf-8")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(max_size=1024, mode="w+", encoding="utf-8")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(1024, "w+", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(max_size=1024, mode="w+", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(1024, "w+", encoding=None, newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(max_size=1024, mode="w+", encoding=None, newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(1024, "w+", encoding="utf-8", newline=None)
    # ruleid: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(max_size=1024, mode="w+", encoding="utf-8", newline=None)
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(1024, "w+", encoding="utf-8", newline="\n")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(max_size=1024, mode="w+", encoding="utf-8", newline="")


def other_file_writers(path: Path, descriptor: int) -> None:
    # ruleid: la-stack.python.text-writes-explicit-policy
    open(path, "w", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    open(path, "w", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    open(path, mode="a", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    open(path, mode="a", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.open("x", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    path.open("x", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.open(mode="r+", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    path.open(mode="r+", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    path.open(mode="r+", encoding=None, newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    io.open(path, "wt", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    io.open(path, "wt", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    io.open(path, mode="a+", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    io.open(path, mode="a+", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    os.fdopen(descriptor, "w", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    os.fdopen(descriptor, "w", encoding="utf-8", newline="\n")
    # ruleid: la-stack.python.text-writes-explicit-policy
    os.fdopen(descriptor, mode="w", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    os.fdopen(descriptor, mode="w", encoding="utf-8", newline="\n")


def byte_writes_and_reads_are_exempt(path: Path, descriptor: int) -> None:
    # ok: la-stack.python.text-writes-explicit-policy
    path.write_bytes(b"line\r\n")
    # ok: la-stack.python.text-writes-explicit-policy
    os.fdopen(descriptor, "wb")
    # ok: la-stack.python.text-writes-explicit-policy
    path.open("wb")
    # ok: la-stack.python.text-writes-explicit-policy
    open(path, "rb")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.NamedTemporaryFile("w+b")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.NamedTemporaryFile()
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.TemporaryFile(mode="wb")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile()
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(1024, "w+b")
    # ok: la-stack.python.text-writes-explicit-policy
    tempfile.SpooledTemporaryFile(max_size=1024, mode="w+b")
    # ok: la-stack.python.text-writes-explicit-policy
    path.open(encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    io.open(path, "r", encoding="utf-8")
    # ok: la-stack.python.text-writes-explicit-policy
    path.read_text(encoding="utf-8")
