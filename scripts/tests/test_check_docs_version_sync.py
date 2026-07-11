"""Tests for documentation and package-version synchronization checks."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

import check_docs_version_sync

if TYPE_CHECKING:
    from pathlib import Path


_CARGO_TOML = '[package]\nname = "other-crate"\nversion = "1.2.3"'
_VERSION = "1.2.3"


def _write_project(
    root: Path,
    *,
    metadata_version: str = _VERSION,
    readme: str | None = None,
) -> None:
    readme_text = (
        readme
        if readme is not None
        else f'other-crate = "{_VERSION}"\n[doc](https://github.com/acgetchell/la-stack/blob/v{_VERSION}/README.md)\n[raw](https://raw.githubusercontent.com/acgetchell/la-stack/v{_VERSION}/README.md)\n'
    )
    files = {
        "Cargo.toml": f"{_CARGO_TOML}\n",
        "Cargo.lock": f'version = 4\n\n[[package]]\nname = "other-crate"\nversion = "{metadata_version}"\n',
        "pyproject.toml": f'[project]\nname = "other-crate-scripts"\nversion = "{metadata_version}"\n',
        "uv.lock": f'version = 1\n\n[[package]]\nname = "other-crate-scripts"\nversion = "{metadata_version}"\nsource = {{ editable = "." }}\n',
        "CITATION.cff": f"cff-version: 1.2.0\nversion: {metadata_version}\n",
        "README.md": readme_text,
    }
    for filename, content in files.items():
        (root / filename).write_text(content, encoding="utf-8")


def test_find_version_mismatches_accepts_matching_dependency_snippets(tmp_path: Path) -> None:
    _write_project(
        tmp_path,
        readme='other-crate = "1.2.3"\nother-crate = { version = "1.2.3", features = ["exact"] }\nla-stack = "0.4.1"',
    )

    assert check_docs_version_sync.find_version_mismatches(tmp_path) == []


def test_find_version_mismatches_reports_stale_dependency_snippets(tmp_path: Path) -> None:
    _write_project(tmp_path)
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "install.md").write_text(
        'other-crate = { version = "1.2.2", features = ["exact"] }\n',
        encoding="utf-8",
    )

    mismatches = check_docs_version_sync.find_version_mismatches(tmp_path)

    assert len(mismatches) == 1
    assert mismatches[0].reference.path == docs / "install.md"
    assert mismatches[0].reference.line == 1
    assert mismatches[0].reference.version == "1.2.2"
    assert mismatches[0].package.name == "other-crate"
    assert mismatches[0].package.version == "1.2.3"


def test_find_version_mismatches_handles_reordered_inline_table_keys(tmp_path: Path) -> None:
    _write_project(tmp_path)
    docs = tmp_path / "docs"
    docs.mkdir()
    install_doc = docs / "install.md"
    install_doc.write_text(
        'other-crate = { features = ["exact"], version = "1.2.2" }\n',
        encoding="utf-8",
    )

    mismatches = check_docs_version_sync.find_version_mismatches(tmp_path)

    assert len(mismatches) == 1
    assert mismatches[0].reference.path == install_doc
    assert mismatches[0].reference.line == 1
    assert mismatches[0].reference.version == "1.2.2"
    assert mismatches[0].package.name == "other-crate"
    assert mismatches[0].package.version == "1.2.3"


def test_find_version_mismatches_reports_all_release_metadata(tmp_path: Path) -> None:
    _write_project(
        tmp_path,
        metadata_version="1.2.2",
    )

    mismatches = check_docs_version_sync.find_version_mismatches(tmp_path)

    assert [(mismatch.reference.kind, mismatch.reference.path.name, mismatch.reference.line, mismatch.reference.version) for mismatch in mismatches] == [
        (check_docs_version_sync.ReferenceKind.CARGO_LOCK, "Cargo.lock", 5, "1.2.2"),
        (check_docs_version_sync.ReferenceKind.PYPROJECT, "pyproject.toml", 3, "1.2.2"),
        (check_docs_version_sync.ReferenceKind.UV_LOCK, "uv.lock", 5, "1.2.2"),
        (check_docs_version_sync.ReferenceKind.CITATION, "CITATION.cff", 2, "1.2.2"),
    ]


def test_find_version_mismatches_reports_readme_tag_links(tmp_path: Path) -> None:
    _write_project(
        tmp_path,
        readme=(
            "[doc](https://github.com/acgetchell/la-stack/blob/v1.2.2/README.md)\n"
            "[raw](https://raw.githubusercontent.com/acgetchell/la-stack/v1.2.1/README.md)\n"
            "[moving](https://github.com/acgetchell/la-stack/blob/main/README.md)\n"
        ),
    )

    mismatches = check_docs_version_sync.find_version_mismatches(tmp_path)

    assert [mismatch.reference.kind for mismatch in mismatches] == [check_docs_version_sync.ReferenceKind.README_TAG_LINK] * 2
    assert [mismatch.reference.line for mismatch in mismatches] == [1, 2]
    assert [mismatch.reference.version for mismatch in mismatches] == ["1.2.2", "1.2.1"]


def test_find_version_mismatches_ignores_historical_docs_and_test_fixtures(tmp_path: Path) -> None:
    _write_project(tmp_path)
    archive = tmp_path / "docs" / "archive"
    archive.mkdir(parents=True)
    fixtures = tmp_path / "tests" / "fixtures"
    fixtures.mkdir(parents=True)
    stale_snippet = 'other-crate = "0.1.0"\n'
    (tmp_path / "CHANGELOG.md").write_text(stale_snippet, encoding="utf-8")
    (archive / "old.md").write_text(stale_snippet, encoding="utf-8")
    (fixtures / "example.md").write_text(stale_snippet, encoding="utf-8")

    assert check_docs_version_sync.find_version_mismatches(tmp_path) == []


def test_find_version_mismatches_rejects_missing_editable_uv_package(tmp_path: Path) -> None:
    _write_project(tmp_path)
    (tmp_path / "uv.lock").write_text(
        'version = 1\n\n[[package]]\nname = "other-crate-scripts"\nversion = "1.2.3"\nsource = { registry = "https://pypi.org/simple" }\n',
        encoding="utf-8",
    )

    with pytest.raises(TypeError, match=r"exactly one uv\.lock editable package"):
        check_docs_version_sync.find_version_mismatches(tmp_path)


def test_find_version_mismatches_rejects_malformed_citation_version(tmp_path: Path) -> None:
    _write_project(tmp_path)
    (tmp_path / "CITATION.cff").write_text('cff-version: 1.2.0\nversion: "\n', encoding="utf-8")

    with pytest.raises(TypeError, match=r"CITATION\.cff:2: top-level version"):
        check_docs_version_sync.find_version_mismatches(tmp_path)
