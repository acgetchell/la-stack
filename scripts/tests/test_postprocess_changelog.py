"""Tests for postprocess_changelog.py — trailing blank line hygiene."""

from __future__ import annotations

from typing import TYPE_CHECKING

from postprocess_changelog import postprocess

if TYPE_CHECKING:
    from pathlib import Path


class TestPostprocess:
    def test_strips_trailing_blank_lines(self, tmp_path: Path) -> None:
        f = tmp_path / "CHANGELOG.md"
        f.write_text("# Changelog\n\n- Item\n\n\n\n", encoding="utf-8")

        postprocess(f)

        assert f.read_text(encoding="utf-8") == "# Changelog\n\n- Item\n"

    def test_preserves_single_trailing_newline(self, tmp_path: Path) -> None:
        f = tmp_path / "CHANGELOG.md"
        f.write_text("# Changelog\n\n- Item\n", encoding="utf-8")

        postprocess(f)

        assert f.read_text(encoding="utf-8") == "# Changelog\n\n- Item\n"

    def test_adds_trailing_newline_if_missing(self, tmp_path: Path) -> None:
        f = tmp_path / "CHANGELOG.md"
        f.write_text("# Changelog\n\n- Item", encoding="utf-8")

        postprocess(f)

        assert f.read_text(encoding="utf-8") == "# Changelog\n\n- Item\n"

    def test_preserves_internal_blank_lines(self, tmp_path: Path) -> None:
        content = "# Changelog\n\n## [1.0.0]\n\n### Added\n\n- Item\n\n\n\n"
        f = tmp_path / "CHANGELOG.md"
        f.write_text(content, encoding="utf-8")

        postprocess(f)

        result = f.read_text(encoding="utf-8")
        # Internal blank lines preserved, only trailing ones stripped
        assert result == "# Changelog\n\n## [1.0.0]\n\n### Added\n\n- Item\n"

    def test_single_newline_file(self, tmp_path: Path) -> None:
        f = tmp_path / "CHANGELOG.md"
        f.write_text("\n", encoding="utf-8")

        postprocess(f)

        assert f.read_text(encoding="utf-8") == "\n"

    def test_empty_file(self, tmp_path: Path) -> None:
        f = tmp_path / "CHANGELOG.md"
        f.write_text("", encoding="utf-8")

        postprocess(f)

        assert f.read_text(encoding="utf-8") == "\n"
