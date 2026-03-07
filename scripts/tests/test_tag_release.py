"""Tests for tag_release.py — annotated tag creation with size-limit handling."""

from __future__ import annotations

from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

import pytest

import tag_release
from tag_release import (
    _GITHUB_TAG_ANNOTATION_LIMIT,
    _github_anchor,
    extract_changelog_section,
    find_changelog,
    parse_version,
    validate_semver,
)

if TYPE_CHECKING:
    from pathlib import Path

# ---------------------------------------------------------------------------
# SemVer validation
# ---------------------------------------------------------------------------


class TestValidateSemver:
    @pytest.mark.parametrize(
        "version",
        [
            "v0.1.0",
            "v1.0.0",
            "v12.34.56",
            "v1.2.3-rc.1",
            "v1.2.3-alpha",
            "v1.2.3+build.42",
            "v1.2.3-beta.1+build.123",
            "v1.0.0-1a",  # digit-prefixed alphanumeric prerelease
            "v1.0.0-0a",  # leading zero OK when not purely numeric
            "v1.0.0-1a.2b",  # dot-separated digit-prefixed IDs
        ],
    )
    def test_valid_versions(self, version: str) -> None:
        validate_semver(version)  # should not raise

    @pytest.mark.parametrize(
        "version",
        [
            "0.1.0",  # missing v prefix
            "v1",  # incomplete
            "v1.2",  # missing patch
            "v01.2.3",  # leading zero
            "v1.02.3",  # leading zero
            "v1.2.03",  # leading zero
            "v1.0.0-01",  # leading zero in purely numeric prerelease
            "vfoo",  # garbage
            "",  # empty
        ],
    )
    def test_invalid_versions(self, version: str) -> None:
        with pytest.raises(ValueError, match="SemVer format"):
            validate_semver(version)


class TestParseVersion:
    def test_strips_v_prefix(self) -> None:
        assert parse_version("v1.2.3") == "1.2.3"

    def test_no_prefix(self) -> None:
        assert parse_version("1.2.3") == "1.2.3"


# ---------------------------------------------------------------------------
# Changelog helpers
# ---------------------------------------------------------------------------


_SAMPLE_CHANGELOG = """\
# Changelog

## [0.2.0] - 2025-03-01

### Added

- Exact determinant sign via `det_sign_exact()`

### Changed

- Bump version to 0.2.0

## [0.1.3] - 2025-02-15

### Fixed

- Minor doc typo
"""


class TestFindChangelog:
    def test_finds_in_current_dir(self, tmp_path: Path) -> None:
        (tmp_path / "CHANGELOG.md").write_text("# Changelog\n", encoding="utf-8")
        result = find_changelog(tmp_path)
        assert result.name == "CHANGELOG.md"

    def test_finds_in_parent_dir(self, tmp_path: Path) -> None:
        (tmp_path / "CHANGELOG.md").write_text("# Changelog\n", encoding="utf-8")
        child = tmp_path / "scripts"
        child.mkdir()
        result = find_changelog(child)
        assert result.name == "CHANGELOG.md"

    def test_raises_when_missing(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match=r"CHANGELOG\.md not found"):
            find_changelog(tmp_path)


class TestExtractChangelogSection:
    def test_extracts_section(self, tmp_path: Path) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(_SAMPLE_CHANGELOG, encoding="utf-8")

        section = extract_changelog_section(changelog, "0.2.0")
        assert "det_sign_exact" in section
        assert "Bump version" in section
        # Should not include content from 0.1.3
        assert "Minor doc typo" not in section

    def test_extracts_older_section(self, tmp_path: Path) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(_SAMPLE_CHANGELOG, encoding="utf-8")

        section = extract_changelog_section(changelog, "0.1.3")
        assert "Minor doc typo" in section
        assert "det_sign_exact" not in section

    def test_raises_for_missing_version(self, tmp_path: Path) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(_SAMPLE_CHANGELOG, encoding="utf-8")

        with pytest.raises(LookupError, match="No changelog section found"):
            extract_changelog_section(changelog, "9.9.9")

    def test_raises_for_empty_section(self, tmp_path: Path) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text("# Changelog\n\n## [1.0.0] - 2025-01-01\n\n## [0.9.0] - 2024-12-01\n", encoding="utf-8")

        with pytest.raises(LookupError, match="empty"):
            extract_changelog_section(changelog, "1.0.0")


# ---------------------------------------------------------------------------
# GitHub anchor generation
# ---------------------------------------------------------------------------


class TestGitHubAnchor:
    """Verify _github_anchor matches github-slugger output."""

    def test_bracketed_heading(self, tmp_path: Path) -> None:
        """Heading ``## [1.0.0] - 2025-01-01`` should strip brackets and dots."""
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(
            "# Changelog\n\n## [1.0.0] - 2025-01-01\n\n- Item\n",
            encoding="utf-8",
        )
        assert _github_anchor(changelog, "1.0.0") == "100---2025-01-01"

    def test_plain_v_heading(self, tmp_path: Path) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(
            "# Changelog\n\n## v0.2.0\n\n- Item\n",
            encoding="utf-8",
        )
        assert _github_anchor(changelog, "0.2.0") == "v020"

    def test_fallback_when_not_found(self, tmp_path: Path) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text("# Changelog\n", encoding="utf-8")
        assert _github_anchor(changelog, "9.9.9") == "v999"

    def test_does_not_match_prerelease_heading(self, tmp_path: Path) -> None:
        """Looking for 1.0.0 must not match ## [1.0.0-rc.1]."""
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(
            "# Changelog\n\n## [1.0.0-rc.1] - 2025-01-01\n\n- Item\n",
            encoding="utf-8",
        )
        # Should fall back since no exact 1.0.0 heading exists
        assert _github_anchor(changelog, "1.0.0") == "v100"


# ---------------------------------------------------------------------------
# Tag size limit handling
# ---------------------------------------------------------------------------


class TestTagSizeLimit:
    def test_small_section_uses_full_content(self, tmp_path: Path) -> None:
        """A normal-sized changelog section should be used as the tag message."""
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(_SAMPLE_CHANGELOG, encoding="utf-8")

        section = extract_changelog_section(changelog, "0.2.0")
        assert len(section.encode("utf-8")) < _GITHUB_TAG_ANNOTATION_LIMIT

    def test_oversized_section_detected(self, tmp_path: Path) -> None:
        """Synthetic oversized changelog should exceed the limit."""
        # Build content > 125KB
        lines = [f"- Item number {i}" for i in range(20_000)]
        big_section = "\n".join(lines)

        changelog_text = f"# Changelog\n\n## [1.0.0] - 2025-01-01\n\n{big_section}\n\n## [0.9.0] - 2024-12-01\n\n- Old item\n"
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(changelog_text, encoding="utf-8")

        section = extract_changelog_section(changelog, "1.0.0")
        assert len(section.encode("utf-8")) > _GITHUB_TAG_ANNOTATION_LIMIT


# ---------------------------------------------------------------------------
# create_tag workflow (mocked git)
# ---------------------------------------------------------------------------


class TestCreateTag:
    @patch("tag_release.run_git_command_with_input")
    @patch("tag_release._tag_exists", return_value=False)
    @patch("tag_release.find_changelog")
    @patch("tag_release.extract_changelog_section", return_value="### Added\n\n- Something new")
    def test_creates_annotated_tag(
        self,
        _mock_extract: MagicMock,
        mock_find: MagicMock,
        _mock_exists: MagicMock,
        mock_git_input: MagicMock,
        tmp_path: Path,
    ) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(_SAMPLE_CHANGELOG, encoding="utf-8")
        mock_find.return_value = changelog

        tag_release.create_tag("v1.0.0")

        mock_git_input.assert_called_once()
        call_args = mock_git_input.call_args
        assert call_args[0][0] == ["tag", "-a", "v1.0.0", "-F", "-"]
        assert "Something new" in call_args[1]["input_data"]

    @patch("tag_release.run_git_command_with_input")
    @patch("tag_release._get_repo_url", return_value="https://github.com/acgetchell/la-stack")
    @patch("tag_release._tag_exists", return_value=False)
    @patch("tag_release.find_changelog")
    def test_oversized_creates_reference_tag(
        self,
        mock_find: MagicMock,
        _mock_exists: MagicMock,
        _mock_url: MagicMock,
        mock_git_input: MagicMock,
        tmp_path: Path,
    ) -> None:
        """When changelog exceeds 125KB, tag message should be a short reference."""
        lines = [f"- Item number {i}" for i in range(20_000)]
        big_section = "\n".join(lines)
        changelog_text = f"# Changelog\n\n## [1.0.0] - 2025-01-01\n\n{big_section}\n\n## [0.9.0] - 2024-12-01\n\n- Old\n"
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(changelog_text, encoding="utf-8")
        mock_find.return_value = changelog

        # Patch extract to return the real oversized content
        with patch("tag_release.extract_changelog_section", return_value=big_section):
            tag_release.create_tag("v1.0.0")

        mock_git_input.assert_called_once()
        tag_message = mock_git_input.call_args[1]["input_data"]
        assert "See full changelog" in tag_message
        assert "CHANGELOG.md" in tag_message
        assert len(tag_message) < 1000

    @patch("tag_release._tag_exists", return_value=True)
    def test_existing_tag_without_force_exits(self, _mock_exists: MagicMock) -> None:
        with pytest.raises(SystemExit):
            tag_release.create_tag("v1.0.0", force=False)

    @patch("tag_release.run_git_command_with_input")
    @patch("tag_release._delete_tag")
    @patch("tag_release._tag_exists", return_value=True)
    @patch("tag_release.find_changelog")
    @patch("tag_release.extract_changelog_section", return_value="### Fixed\n\n- Bug fix")
    def test_force_recreates_tag(
        self,
        _mock_extract: MagicMock,
        mock_find: MagicMock,
        _mock_exists: MagicMock,
        mock_delete: MagicMock,
        mock_git_input: MagicMock,
        tmp_path: Path,
    ) -> None:
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text(_SAMPLE_CHANGELOG, encoding="utf-8")
        mock_find.return_value = changelog

        tag_release.create_tag("v1.0.0", force=True)

        mock_delete.assert_called_once_with("v1.0.0")
        mock_git_input.assert_called_once()

    @patch("tag_release._tag_exists", return_value=True)
    @patch("tag_release.find_changelog")
    @patch("tag_release.extract_changelog_section", side_effect=LookupError("not found"))
    @patch("tag_release._delete_tag")
    def test_force_does_not_delete_tag_if_changelog_fails(
        self,
        mock_delete: MagicMock,
        _mock_extract: MagicMock,
        mock_find: MagicMock,
        _mock_exists: MagicMock,
        tmp_path: Path,
    ) -> None:
        """Tag must not be deleted if changelog extraction fails."""
        changelog = tmp_path / "CHANGELOG.md"
        changelog.write_text("# Changelog\n", encoding="utf-8")
        mock_find.return_value = changelog

        with pytest.raises(LookupError):
            tag_release.create_tag("v1.0.0", force=True)

        mock_delete.assert_not_called()
