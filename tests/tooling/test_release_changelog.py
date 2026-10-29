"""Consumer citation policy uses the installed deterministic date resolver."""

from pathlib import Path

import pytest
from research_repo_tools.changelog import resolve_date
from research_repo_tools.config import load
from research_repo_tools.just_inspect import dry_run

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def consumer(tmp_path: Path) -> Path:
    (tmp_path / "pyproject.toml").write_bytes((ROOT / "pyproject.toml").read_bytes())
    return tmp_path


@pytest.mark.parametrize("quotation", ["", "'", '"'])
def test_omitted_date_uses_recorded_citation_date(consumer: Path, quotation: str) -> None:
    citation = consumer / "CITATION.cff"
    contents = f"version: 0.5.1\ndate-released: {quotation}2000-02-29{quotation}\n"
    citation.write_text(contents, encoding="utf-8")
    assert resolve_date(load(root=consumer), tag="v0.5.1") == "2000-02-29"
    assert citation.read_text(encoding="utf-8") == contents


def test_explicit_date_and_untagged_preview_need_no_citation(consumer: Path) -> None:
    settings = load(root=consumer)
    assert resolve_date(settings, tag="v0.5.1", released="2026-10-09") == "2026-10-09"
    assert resolve_date(settings) is None
    assert not (consumer / "CITATION.cff").exists()


@pytest.mark.parametrize(
    "contents",
    [
        None,
        "[",
        "- item\n",
        "version: 0.5.1\n",
        "date-released: 2000-02-29\n",
        "version: 0.5.0\ndate-released: 2000-02-29\n",
        "version: 0.5.1\ndate-released: 2001-02-29\n",
        "version: 0.5.1\ndate-released: 20000229\n",
        "version: 0.5.1\ndate-released: \n",
        "version: 0.5.1\ndate-released: [2000-02-29]\n",
        "version: 0.5.1\ndate-released: 2000-02-29\ndate-released: 2000-03-01\n",
        "version: 0.5.1\nversion: 0.5.2\ndate-released: 2000-02-29\n",
    ],
)
def test_invalid_default_metadata_never_changes_release_files(consumer: Path, contents: str | None) -> None:
    citation = consumer / "CITATION.cff"
    if contents is not None:
        citation.write_text(contents, encoding="utf-8")
    changelog = consumer / "CHANGELOG.md"
    changelog.write_bytes(b"retained changelog\n")
    with pytest.raises(ValueError, match=r"citation|YYYY-MM-DD"):
        resolve_date(load(root=consumer), tag="v0.5.1")
    assert changelog.read_bytes() == b"retained changelog\n"
    if contents is not None:
        assert citation.read_text(encoding="utf-8") == contents


def test_recipe_omits_absent_date_and_delegates_to_installed_cli() -> None:
    implicit = dry_run(ROOT, "changelog-release", ["v0.5.1"]).stderr
    explicit = dry_run(ROOT, "changelog-release", ["v0.5.1", "2000-02-29"]).stderr
    assert "changelog generate --tag 'v0.5.1'" in implicit
    assert "--date" not in implicit
    assert "changelog generate --tag 'v0.5.1' --date '2000-02-29'" in explicit
    assert 'changelog generate --dry-run "$@"' in dry_run(ROOT, "changelog-preview", ["--tag", "v0.5.1"]).stderr
