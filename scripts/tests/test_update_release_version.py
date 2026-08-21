"""Tests for transactional release-version updates."""

from typing import TYPE_CHECKING

import pytest

import check_docs_version_sync
import update_release_version

if TYPE_CHECKING:
    from pathlib import Path


def _write_project(root: Path, *, metadata_version: str = "1.2.2", dependency_version: str = "1.2.2") -> None:
    files = {
        "Cargo.toml": f'[package]\nname = "other-crate"\nversion = "{metadata_version}"\n',
        "Cargo.lock": (
            f'version = 4\n\n[[package]]\nname = "either"\nversion = "1.17.0"\n\n[[package]]\nname = "other-crate"\nversion = "{metadata_version}"\n'
        ),
        "pyproject.toml": f'[project]\nname = "other-crate-scripts"\nversion = "{metadata_version}"\n',
        "uv.lock": (f'version = 1\n\n[[package]]\nname = "other-crate-scripts"\nversion = "{metadata_version}"\nsource = {{ editable = "." }}\n'),
        "CITATION.cff": (f'cff-version: 1.2.0\nversion: {metadata_version}\ndate-released: 2026-07-13\ndoi: "10.5281/zenodo.12345"\n'),
        "README.md": (
            f'other-crate = "{dependency_version}"\n'
            f'other-crate = {{ version = "{dependency_version}", features = ["exact"] }}\n'
            f"[doc](https://github.com/acgetchell/la-stack/blob/v{metadata_version}/README.md)\n"
            f"[raw](https://raw.githubusercontent.com/acgetchell/la-stack/v{metadata_version}/README.md)\n"
            f"[csv](https://github.com/acgetchell/la-stack/blob/v{metadata_version}/docs/assets/bench/result.csv)\n"
            f"[provenance](https://github.com/acgetchell/la-stack/blob/v{metadata_version}/docs/assets/bench/result.provenance.json)\n"
            f"[svg](https://raw.githubusercontent.com/acgetchell/la-stack/v{metadata_version}/docs/assets/bench/result.svg)\n"
            "<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->\n"
            "| unchanged benchmark comparison |\n"
            "<!-- BENCH_TABLE:lu_solve:median:new:END -->\n"
        ),
        "CHANGELOG.md": "# Changelog\n\n## [1.2.2] - 2026-07-13\n",
    }
    for filename, content in files.items():
        (root / filename).write_text(content, encoding="utf-8")
    docs = root / "docs"
    docs.mkdir()
    (docs / "BENCHMARKING.md").write_text(
        "just performance-release v1.2.2 v1.2.1\nThis generates a local `v1.2.1` `vs_linalg` baseline.\nHistorical v1.2.1 behavior remains documented.\n",
        encoding="utf-8",
    )


def _previous() -> update_release_version.ReleaseTag:
    return update_release_version.parse_release_tag("v1.2.2")


def test_update_release_version_updates_all_current_surfaces_without_dependency_upgrades(tmp_path: Path) -> None:
    _write_project(tmp_path)

    summary = update_release_version.update_release_version(
        tmp_path,
        "v1.2.3",
        previous=_previous(),
        release_date="2026-08-20",
    )

    assert summary.target.tag == "v1.2.3"
    assert summary.previous.tag == "v1.2.2"
    assert summary.release_date == "2026-08-20"
    assert summary.changed_paths
    assert 'name = "either"\nversion = "1.17.0"' in (tmp_path / "Cargo.lock").read_text(encoding="utf-8")
    assert 'name = "other-crate"\nversion = "1.2.3"' in (tmp_path / "Cargo.lock").read_text(encoding="utf-8")
    assert 'name = "other-crate-scripts"\nversion = "1.2.3"' in (tmp_path / "uv.lock").read_text(encoding="utf-8")
    assert 'version = "1.2.3"' in (tmp_path / "Cargo.toml").read_text(encoding="utf-8")
    assert 'version = "1.2.3"' in (tmp_path / "pyproject.toml").read_text(encoding="utf-8")
    citation = (tmp_path / "CITATION.cff").read_text(encoding="utf-8")
    assert "version: 1.2.3" in citation
    assert "date-released: 2026-08-20" in citation
    assert 'doi: "10.5281/zenodo.12345"' in citation
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert 'other-crate = "1.2.3"' in readme
    assert 'version = "1.2.3"' in readme
    assert readme.count("v1.2.3") == 2
    assert readme.count("v1.2.2/docs/assets/bench/") == 3
    assert "| unchanged benchmark comparison |" in readme
    benchmarking = (tmp_path / "docs" / "BENCHMARKING.md").read_text(encoding="utf-8")
    assert "just performance-release v1.2.3 v1.2.2" in benchmarking
    assert "This generates a local `v1.2.2` `vs_linalg` baseline." in benchmarking
    assert "Historical v1.2.1 behavior remains documented." in benchmarking
    assert "## [1.2.2] - 2026-07-13" in (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8")
    assert check_docs_version_sync.find_version_mismatches(tmp_path) == []


def test_update_release_version_is_idempotent(tmp_path: Path) -> None:
    _write_project(tmp_path)
    kwargs = {"previous": _previous(), "release_date": "2026-08-20"}

    first = update_release_version.update_release_version(tmp_path, "v1.2.3", **kwargs)
    second = update_release_version.update_release_version(tmp_path, "v1.2.3", **kwargs)

    assert first.changed_paths
    assert second.changed_paths == ()


def test_update_release_version_advances_existing_release_dates_together(tmp_path: Path) -> None:
    _write_project(tmp_path, metadata_version="1.2.3", dependency_version="1.2.3")
    citation = tmp_path / "CITATION.cff"
    citation.write_text(citation.read_text(encoding="utf-8").replace("2026-07-13", "2026-08-20"), encoding="utf-8")
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text("# Changelog\n\n## [1.2.3] - 2026-08-20\n", encoding="utf-8")
    benchmarking = tmp_path / "docs" / "BENCHMARKING.md"
    benchmarking.write_text(
        "just performance-release v1.2.3 v1.2.2\nThis generates a local `v1.2.2` `vs_linalg` baseline.\n",
        encoding="utf-8",
    )

    summary = update_release_version.update_release_version(
        tmp_path,
        "v1.2.3",
        previous=_previous(),
        release_date="2026-08-21",
    )

    assert summary.changed_paths == (changelog, citation)
    assert "date-released: 2026-08-21" in citation.read_text(encoding="utf-8")
    assert "## [1.2.3] - 2026-08-21" in changelog.read_text(encoding="utf-8")


def test_select_previous_release_tag_uses_latest_stable_published_tag() -> None:
    target = update_release_version.parse_release_tag("v1.3.0")

    previous = update_release_version.select_previous_release_tag(
        ["v1.1.9", "v1.2.0-rc.1", "not-a-release", "v1.2.0"],
        target,
    )

    assert previous.tag == "v1.2.0"


def test_select_previous_release_tag_ignores_already_published_target() -> None:
    target = update_release_version.parse_release_tag("v1.3.0")

    previous = update_release_version.select_previous_release_tag(["v1.2.0", "v1.3.0"], target)

    assert previous.tag == "v1.2.0"


@pytest.mark.parametrize("target", ["1.2.3", "v1.2", "v01.2.3", "v1.2.3-rc.1"])
def test_parse_release_tag_rejects_non_stable_tag_forms(target: str) -> None:
    with pytest.raises(ValueError, match=r"stable tag in vX\.Y\.Z form"):
        update_release_version.parse_release_tag(target)


def test_select_previous_release_tag_rejects_target_older_than_a_published_release() -> None:
    target = update_release_version.parse_release_tag("v1.2.3")

    with pytest.raises(ValueError, match="older than published"):
        update_release_version.select_previous_release_tag(["v1.2.3", "v1.3.0"], target)


def test_infer_previous_release_tag_uses_published_github_releases(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(update_release_version, "published_stable_release_tags", lambda root: ["v1.2.1", "v1.2.2"] if root == tmp_path else [])

    previous = update_release_version.infer_previous_release_tag(tmp_path, update_release_version.parse_release_tag("v1.2.3"))

    assert previous.tag == "v1.2.2"


def test_unexpected_version_fails_before_writing(tmp_path: Path) -> None:
    _write_project(tmp_path, dependency_version="1.0.0")
    originals = {path: path.read_text(encoding="utf-8") for path in tmp_path.rglob("*") if path.is_file()}

    with pytest.raises(ValueError, match="unexpected other-crate dependency version"):
        update_release_version.update_release_version(
            tmp_path,
            "v1.2.3",
            previous=_previous(),
            release_date="2026-08-20",
        )

    assert {path: path.read_text(encoding="utf-8") for path in originals} == originals


def test_validation_failure_rolls_back_every_changed_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _write_project(tmp_path)
    originals = {path: path.read_text(encoding="utf-8") for path in tmp_path.rglob("*") if path.is_file()}

    def fail_validation(*_args: object) -> None:
        msg = "simulated validation failure"
        raise ValueError(msg)

    monkeypatch.setattr(update_release_version, "_validate_updated_root", fail_validation)

    with pytest.raises(ValueError, match="simulated validation failure"):
        update_release_version.update_release_version(
            tmp_path,
            "v1.2.3",
            previous=_previous(),
            release_date="2026-08-20",
        )

    assert {path: path.read_text(encoding="utf-8") for path in originals} == originals


def test_sync_changelog_release_date_uses_citation_date(tmp_path: Path) -> None:
    _write_project(tmp_path)
    update_release_version.update_release_version(
        tmp_path,
        "v1.2.3",
        previous=_previous(),
        release_date="2026-08-20",
    )
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text("# Changelog\n\n## [1.2.3] - 2026-08-21\n", encoding="utf-8")

    changed, release_date = update_release_version.sync_changelog_release_date(
        tmp_path,
        "v1.2.3",
        previous=_previous(),
    )

    assert changed == (changelog,)
    assert release_date == "2026-08-20"
    assert "## [1.2.3] - 2026-08-20" in changelog.read_text(encoding="utf-8")
    assert check_docs_version_sync.find_version_mismatches(tmp_path) == []


def test_main_supports_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit, match="0"):
        update_release_version.main(["--help"])

    assert "Target stable release tag" in capsys.readouterr().out
