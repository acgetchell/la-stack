from __future__ import annotations

from typing import TYPE_CHECKING

import check_docs_version_sync

if TYPE_CHECKING:
    from pathlib import Path


def test_find_version_mismatches_accepts_matching_dependency_snippets(tmp_path: Path) -> None:
    (tmp_path / "Cargo.toml").write_text(
        "\n".join(
            [
                "[package]",
                'name = "other-crate"',
                'version = "1.2.3"',
            ]
        ),
        encoding="utf-8",
    )
    (tmp_path / "README.md").write_text(
        "\n".join(
            [
                'other-crate = "1.2.3"',
                'other-crate = { version = "1.2.3", features = ["exact"] }',
                'la-stack = "0.4.1"',
            ]
        ),
        encoding="utf-8",
    )

    assert check_docs_version_sync.find_version_mismatches(tmp_path) == []


def test_find_version_mismatches_reports_stale_dependency_snippets(tmp_path: Path) -> None:
    (tmp_path / "Cargo.toml").write_text(
        "\n".join(
            [
                "[package]",
                'name = "other-crate"',
                'version = "1.2.3"',
            ]
        ),
        encoding="utf-8",
    )
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "install.md").write_text(
        'other-crate = { version = "1.2.2", features = ["exact"] }\n',
        encoding="utf-8",
    )

    mismatches = check_docs_version_sync.find_version_mismatches(tmp_path)

    assert len(mismatches) == 1
    assert mismatches[0].snippet.path == docs / "install.md"
    assert mismatches[0].snippet.line == 1
    assert mismatches[0].snippet.version == "1.2.2"
    assert mismatches[0].package.name == "other-crate"
    assert mismatches[0].package.version == "1.2.3"


def test_find_version_mismatches_handles_reordered_inline_table_keys(tmp_path: Path) -> None:
    (tmp_path / "Cargo.toml").write_text(
        "\n".join(
            [
                "[package]",
                'name = "other-crate"',
                'version = "1.2.3"',
            ]
        ),
        encoding="utf-8",
    )
    docs = tmp_path / "docs"
    docs.mkdir()
    install_doc = docs / "install.md"
    install_doc.write_text(
        'other-crate = { features = ["exact"], version = "1.2.2" }\n',
        encoding="utf-8",
    )

    mismatches = check_docs_version_sync.find_version_mismatches(tmp_path)

    assert len(mismatches) == 1
    assert mismatches[0].snippet.path == install_doc
    assert mismatches[0].snippet.line == 1
    assert mismatches[0].snippet.version == "1.2.2"
    assert mismatches[0].package.name == "other-crate"
    assert mismatches[0].package.version == "1.2.3"
