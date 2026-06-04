from __future__ import annotations

import check_docs_version_sync


def test_find_version_mismatches_accepts_matching_dependency_snippets(tmp_path) -> None:
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


def test_find_version_mismatches_reports_stale_dependency_snippets(tmp_path) -> None:
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
