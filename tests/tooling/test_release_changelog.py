"""Recorded citation dates make prospective changelog generation repeatable."""

from pathlib import Path

import pytest
from research_repo_tools.just_inspect import dry_run

from tooling import release_changelog as adapter

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("quotation", ["", "'", '"'])
def test_omitted_date_uses_recorded_citation_date_for_generation_and_preview(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    quotation: str,
) -> None:
    citation = tmp_path / "CITATION.cff"
    contents = f"version: 0.5.1\ndate-released: {quotation}2000-02-29{quotation}\n"
    citation.write_text(contents, encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    calls: list[list[str]] = []

    def shared_cli(argv: list[str]) -> int:
        calls.append(argv)
        return 0

    monkeypatch.setattr(adapter, "tooling_main", shared_cli)
    for extra in ([], ["--date", ""], ["--dry-run"]):
        assert adapter.main(["--tag", "v0.5.1", *extra]) == 0
        expected = ["changelog", "generate", "--tag", "v0.5.1", "--date", "2000-02-29"]
        if extra == ["--dry-run"]:
            expected.append("--dry-run")
        assert calls[-1] == expected
    assert citation.read_text(encoding="utf-8") == contents


@pytest.mark.parametrize("arguments", [["--tag", "v0.5.1", "--date", "2026-10-09"], []])
def test_explicit_date_and_untagged_preview_do_not_require_citation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    arguments: list[str],
) -> None:
    monkeypatch.chdir(tmp_path)
    calls: list[list[str]] = []

    def shared_cli(argv: list[str]) -> int:
        calls.append(argv)
        return 7

    monkeypatch.setattr(adapter, "tooling_main", shared_cli)
    assert adapter.main([*arguments, "--dry-run"]) == 7
    assert calls == [["changelog", "generate", *arguments, "--dry-run"]]
    assert not (tmp_path / "CITATION.cff").exists()


@pytest.mark.parametrize(
    ("contents", "diagnostic"),
    [
        (None, "CITATION.cff"),
        ("[", "while parsing"),
        ("- item\n", "top-level mapping"),
        ("version: 0.5.1\n", "top-level date-released"),
        ("date-released: 2000-02-29\n", "top-level version"),
        ("version: 0.5.0\ndate-released: 2000-02-29\n", "does not match v0.5.1"),
        ("version: 0.5.1\ndate-released: 2001-02-29\n", "day"),
        ("version: 0.5.1\ndate-released: 20000229\n", "YYYY-MM-DD"),
        ("version: 0.5.1\ndate-released: \n", "nonempty top-level date-released"),
        ("version: 0.5.1\ndate-released: [2000-02-29]\n", "top-level date-released scalar"),
        ("version: 0.5.1\ndate-released: 2000-02-29\ndate-released: 2000-03-01\n", "exactly one"),
        ("version: 0.5.1\nversion: 0.5.2\ndate-released: 2000-02-29\n", "exactly one"),
    ],
)
def test_invalid_default_metadata_stops_before_generation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    contents: str | None,
    diagnostic: str,
) -> None:
    citation = tmp_path / "CITATION.cff"
    if contents is not None:
        citation.write_text(contents, encoding="utf-8")
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_bytes(b"retained changelog\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(adapter, "tooling_main", lambda _argv: pytest.fail("invalid metadata must not generate anything"))
    assert adapter.main(["--tag", "v0.5.1"]) == 1
    output = capsys.readouterr()
    assert output.out == ""
    assert "Release changelog rejected:" in output.err
    assert diagnostic in output.err
    assert "just update-version" in output.err
    assert changelog.read_bytes() == b"retained changelog\n"
    assert not (tmp_path / "docs").exists()
    if contents is not None:
        assert citation.read_text(encoding="utf-8") == contents


def test_optional_recipe_argument_and_alias_delegate_to_the_date_adapter() -> None:
    implicit = dry_run(ROOT, "changelog-unreleased", ["v0.5.1"]).stderr
    explicit = dry_run(ROOT, "changelog-release", ["v0.5.1", "2000-02-29"]).stderr
    assert "tooling.release_changelog --tag 'v0.5.1' --date ''" in implicit
    assert "tooling.release_changelog --tag 'v0.5.1' --date '2000-02-29'" in explicit
    assert 'tooling.release_changelog --dry-run "$@"' in dry_run(ROOT, "changelog-preview", ["--tag", "v0.5.1"]).stderr
