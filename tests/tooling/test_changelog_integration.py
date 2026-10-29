"""Consumer policy and authored history against the installed published CLI."""

import shutil
import tomllib
from pathlib import Path
from typing import TYPE_CHECKING

from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run, inspect_justfile
from research_repo_tools.process import run_command

if TYPE_CHECKING:
    import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_consumer_changelog_policy_uses_the_shared_formatter() -> None:
    """Retain la-stack's repository, formatter and authored-body policy."""
    manifest = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest["tool"]["research-repo-tools"]["changelog"] == {
        "owner": "acgetchell",
        "repository": "la-stack",
        "formatter": "changelog-rumdl.toml",
        "dependency-bodies": "preserve",
        "date-source": "citation",
        "citation": "CITATION.cff",
    }
    assert not (ROOT / "cliff.toml").exists()


def test_consumer_generation_preserves_an_authored_dependency_note(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Exercise the installed CLI with la-stack's actual formatter and policy."""
    for name in ("pyproject.toml", "changelog-rumdl.toml", ".python-version", "rust-toolchain.toml"):
        shutil.copyfile(ROOT / name, tmp_path / name)
    run_command("git", ["--no-pager", "init", "--quiet"], cwd=tmp_path)
    for key, value in (("user.name", "Fixture"), ("user.email", "fixture@example.invalid"), ("commit.gpgsign", "false"), ("tag.gpgsign", "false")):
        run_command("git", ["--no-pager", "config", key, value], cwd=tmp_path)
    run_command("git", ["--no-pager", "commit", "--quiet", "--allow-empty", "-m", "feat: initial history"], cwd=tmp_path)
    run_command("git", ["--no-pager", "tag", "v0.4.5"], cwd=tmp_path)
    note = "Read [release notes](https://github.com/astral-sh/ruff/compare/0.16.1...0.16.2) before changing `Vector<D>` checks."
    run_command("git", ["--no-pager", "commit", "--quiet", "--allow-empty", "-m", f"chore(deps-dev): bump ruff\n\n{note}"], cwd=tmp_path)
    (tmp_path / "CHANGELOG.md").write_text("# Changelog\n\n## [0.4.5] - 2000-01-02\n\n- Historical date.\n", encoding="utf-8")
    citation = tmp_path / "CITATION.cff"
    citation.write_text('version: "0.4.6"\ndate-released: 2000-02-29\n', encoding="utf-8")
    original_citation = citation.read_bytes()
    monkeypatch.chdir(tmp_path)
    assert main(["changelog", "generate", "--tag", "v0.4.6"]) == 0
    generated = (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "### Dependencies" in generated
    assert note in generated
    assert "## [0.4.6] - 2000-02-29" in generated
    assert "## [0.4.5] - 2000-01-02" in generated
    assert main(["changelog", "generate", "--tag", "v0.4.6"]) == 0
    assert (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8") == generated
    assert citation.read_bytes() == original_citation
    assert main(["changelog", "generate", "--tag", "v0.4.6", "--dry-run"]) == 0
    assert (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8") == generated


def test_consumer_recipes_and_retained_history(capsys: pytest.CaptureFixture[str]) -> None:
    """Native recipes delegate generation and all historical notes remain readable."""
    recipes = inspect_justfile(ROOT)
    assert recipes.aliases == {}
    assert "changelog-release" in recipes.recipes
    result = dry_run(ROOT, "changelog-preview", ["--tag", "v0.4.7", "--date", "2026-10-09"])
    assert 'research-repo-tools changelog generate --dry-run "$@"' in result.stderr
    assert main(["--root", str(ROOT), "changelog", "check"]) == 0
    for tag in ("v0.1.0", "v0.2.0", "v0.3.0", "v0.4.6"):
        assert main(["--root", str(ROOT), "changelog", "notes", tag]) == 0
        assert capsys.readouterr().out.strip()
