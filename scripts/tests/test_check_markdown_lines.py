"""Regression coverage for locale-independent Markdown line-length checks."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

CHECKER = Path(__file__).resolve().parents[1] / "check_markdown_lines.py"


@pytest.mark.parametrize(
    ("content", "diagnostic"),
    [
        (("a" * 159 + "—\n").encode(), ""),
        (("a" * 159 + "🗺\r\n").encode(), ""),
        (("a" * 159 + "—").encode(), ""),
        (("a" * 160 + "—\n").encode(), ":1: line length 161 exceeds 160"),
        (b"a" * 161, ":1: line length 161 exceeds 160"),
        (b"a" * 159 + b"  \n", ":1: line length 161 exceeds 160"),
        (b"|" + b"a" * 200 + b"\n", ""),
        (b"short\r\n" + b"a" * 161 + b"\r\n", ":2: line length 161 exceeds 160"),
        (b"invalid UTF-8: \xff\n", "decode"),
    ],
    ids=["unicode-limit", "emoji-crlf", "no-final-newline", "unicode-too-long", "ascii-too-long", "trailing-spaces", "table", "line-number", "invalid-utf8"],
)
def test_raw_line_limit_under_c_locale(tmp_path: Path, content: bytes, diagnostic: str) -> None:
    """Count UTF-8 characters under the locale that exposed the Windows failure."""
    markdown = tmp_path / "document with spaces.md"
    markdown.write_bytes(content)
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    result = subprocess.run(  # noqa: S603 - fixed local checker and test-owned input.
        [sys.executable, str(CHECKER), str(markdown)],
        check=False,
        capture_output=True,
        encoding="utf-8",
        env=environment,
        timeout=30,
    )

    assert result.stdout == ""
    if diagnostic:
        assert result.returncode == 1
        assert str(markdown) in result.stderr
        assert diagnostic in result.stderr
        assert "Markdown raw line-length check failed." in result.stderr
    else:
        assert result.returncode == 0
        assert result.stderr == ""
