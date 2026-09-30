"""Consumer integration with the pinned, published changelog CLI."""

import json
import re
import shutil
import tomllib
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import pytest
from research_repo_tools.cli import main

from subprocess_utils import run_git_command, run_safe_command

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def consumer(tmp_path: Path) -> Path:
    """Use real consumer configuration and policy in a disposable Git repository."""
    for name in ("pyproject.toml", "cliff.toml", "changelog-rumdl.toml", "Cargo.toml", "README.md"):
        shutil.copyfile(REPO_ROOT / name, tmp_path / name)
    run_git_command(["init", "--quiet"], cwd=tmp_path)
    run_git_command(["config", "user.name", "Changelog fixture"], cwd=tmp_path)
    run_git_command(["config", "user.email", "fixture@example.invalid"], cwd=tmp_path)
    return tmp_path


def cli(root: Path, *args: str) -> int:
    """Call only the documented installed-package entry point."""
    return main(["--root", str(root), "changelog", *args])


def markdown_bytes(root: Path) -> dict[Path, bytes]:
    """Capture changelog outputs for read-only and failed-publication assertions."""
    paths = [root / "CHANGELOG.md", *(root / "docs/archives/changelog").glob("*.md")]
    return {path.relative_to(root): path.read_bytes() for path in paths if path.is_file()}


def commit(root: Path, message: str) -> None:
    """Append one real commit to the fixture history."""
    run_git_command(["commit", "--quiet", "--allow-empty", "-m", message], cwd=root)


def test_published_package_is_exactly_pinned_and_included_in_dev() -> None:
    """Normal sync installs a PyPI release, with no sibling or editable override."""
    manifest = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    groups = manifest["dependency-groups"]
    assert groups["tooling"] == ["research-repo-tools==0.1.7"]
    assert {"include-group": "tooling"} in groups["dev"]
    assert version("research-repo-tools") == "0.1.7"
    lock = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    package = next(package for package in lock["package"] if package["name"] == "research-repo-tools")
    assert package["version"] == "0.1.7"
    assert package["source"] == {"registry": "https://pypi.org/simple"}


def test_archive_rotation_preserves_notes_references_and_relative_links(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Root and archived notes remain discoverable after a minor-series rotation."""
    changelog = consumer / "CHANGELOG.md"
    changelog.write_text(
        "# Changelog\n\n## [Unreleased]\n\n- Pending.\n\n"
        "## [0.5.0] - 2026-09-30\n\n- New minor.\n\n"
        "## [0.4.6] - 2026-09-08\n\n- Keep `Vector<D>` and [guide][guide].\n\n"
        "## [0.4.5] - 2026-08-21\n\n- Prior patch.\n\n"
        "[guide]: README.md\n"
        "[0.5.0]: https://github.com/acgetchell/la-stack/compare/v0.4.6...v0.5.0\n"
        "[0.4.6]: https://github.com/acgetchell/la-stack/compare/v0.4.5...v0.4.6\n",
        encoding="utf-8",
    )
    assert cli(consumer, "archive") == 0
    root = changelog.read_text(encoding="utf-8")
    archive = consumer / "docs/archives/changelog/0.4.md"
    archived = archive.read_text(encoding="utf-8")
    assert "## [Unreleased]" in root
    assert "## [0.5.0] - 2026-09-30" in root
    assert "## [0.4.6]" not in root
    assert "[0.4.x](docs/archives/changelog/0.4.md)" in root
    assert "## [0.4.6] - 2026-09-08" in archived
    assert "## [0.4.5] - 2026-08-21" in archived
    assert "`Vector<D>`" in archived
    assert "[guide]: ../../../README.md" in archived
    before = markdown_bytes(consumer)
    assert cli(consumer, "archive") == 0
    assert cli(consumer, "check") == 0
    assert markdown_bytes(consumer) == before
    capsys.readouterr()
    assert cli(consumer, "notes", "v0.4.6") == 0
    notes = capsys.readouterr().out
    assert "Keep `Vector<D>`" in notes
    assert "[guide]: ../../../README.md" in notes
    assert "Prior patch" not in notes
    assert "[0.5.0]:" not in notes
    assert cli(consumer, "notes", "v0.5.0") == 0
    assert "New minor" in capsys.readouterr().out
    assert markdown_bytes(consumer) == before


def test_generation_uses_consumer_policy_dates_and_transactional_preview(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Real Git/git-cliff generation retains authored dependency links and explicit dates."""
    commit(consumer, "feat: add exact arithmetic")
    run_git_command(["tag", "v0.4.5"], cwd=consumer)
    commit(
        consumer,
        "chore(deps-dev): bump ruff from 0.16.1 to 0.16.2\n\nRead [release notes](https://github.com/astral-sh/ruff/compare/0.16.1...0.16.2).",
    )
    commit(consumer, "fix!: preserve matrix invariants (#42)\n\nBREAKING CHANGE: preserve exact errors and migration instructions.")
    (consumer / "CHANGELOG.md").write_text("# Changelog\n\n## [0.4.5] - 2000-01-02\n\n- Authored date.\n", encoding="utf-8")
    metadata_before = {name: (consumer / name).read_bytes() for name in ("Cargo.toml", "pyproject.toml")}
    before = markdown_bytes(consumer)
    assert cli(consumer, "generate", "--tag", "v0.5.0", "--date", "2026-09-30", "--dry-run") == 0
    preview = capsys.readouterr().out
    assert "## [0.5.0] - 2026-09-30" in preview
    assert "https://github.com/astral-sh/ruff/compare/0.16.1...0.16.2" in preview
    assert "https://github.com/acgetchell/la-stack/pull/42" in preview
    assert "preserve exact errors and migration instructions" in preview
    assert markdown_bytes(consumer) == before
    assert not (consumer / "docs/archives/changelog").exists()
    assert cli(consumer, "generate", "--tag", "v0.5.0", "--date", "2026-09-30") == 0
    assert (consumer / "CHANGELOG.md").read_text(encoding="utf-8") == preview
    archive = consumer / "docs/archives/changelog/0.4.md"
    assert "## [0.4.5] - 2000-01-02" in archive.read_text(encoding="utf-8")
    assert {name: (consumer / name).read_bytes() for name in metadata_before} == metadata_before
    generated = markdown_bytes(consumer)
    assert cli(consumer, "generate", "--tag", "v0.5.0", "--date", "2026-09-30") == 0
    assert markdown_bytes(consumer) == generated


@pytest.mark.parametrize(
    ("args", "diagnostic"),
    [
        (("generate", "--tag", "v0.5.0"), "both --tag and --date"),
        (("generate", "--date", "2026-09-30"), "both --tag and --date"),
        (("generate", "--tag", "v0.5.0", "--date", "2026-02-30"), "day"),
        (("generate", "--tag", "v01.5.0", "--date", "2026-09-30"), "SemVer format"),
        (("notes", "0.4.6"), "SemVer format"),
        (("notes", "v9.9.9"), "not found"),
    ],
)
def test_cli_errors_are_actionable_and_preserve_outputs(consumer: Path, capsys: pytest.CaptureFixture[str], args: tuple[str, ...], diagnostic: str) -> None:
    (consumer / "CHANGELOG.md").write_text("# Changelog\n\n## [0.4.6] - 2026-09-08\n\n- Retained.\n", encoding="utf-8")
    before = markdown_bytes(consumer)
    assert cli(consumer, *args) == 1
    error = capsys.readouterr().err
    assert "research-repo-tools:" in error
    assert diagnostic.casefold() in error.casefold()
    assert "Traceback" not in error
    assert markdown_bytes(consumer) == before


@pytest.mark.parametrize(
    ("headings", "diagnostic"),
    [
        ("## [0.4.5]\n\n- Older.\n\n## [0.4.6]\n\n- Newer.\n", "out of order"),
        ("## [0.4.6]\n\n- One.\n\n## [0.4.6]\n\n- Duplicate.\n", "duplicate"),
        ("## [0.4.6] - 2026-02-30\n\n- Invalid date.\n", "date"),
    ],
)
def test_whole_history_errors_do_not_rotate_files(consumer: Path, capsys: pytest.CaptureFixture[str], headings: str, diagnostic: str) -> None:
    (consumer / "CHANGELOG.md").write_text("# Changelog\n\n" + headings, encoding="utf-8")
    before = markdown_bytes(consumer)
    for command in ("check", "archive"):
        assert cli(consumer, command) == 1
        assert diagnostic in capsys.readouterr().err.casefold()
        assert markdown_bytes(consumer) == before


def test_archive_conflicts_preserve_retained_history(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Changed notes for an already archived release fail before any replacement."""
    (consumer / "CHANGELOG.md").write_text(
        "# Changelog\n\n## [0.5.0] - 2026-09-30\n\n- Current.\n\n## [0.4.6] - 2026-09-08\n\n- Conflicting notes.\n",
        encoding="utf-8",
    )
    archive = consumer / "docs/archives/changelog/0.4.md"
    archive.parent.mkdir(parents=True)
    archive.write_text("# Changelog - 0.4.x\n\n## [0.4.6] - 2026-09-08\n\n- Authored retained notes.\n", encoding="utf-8")
    before = markdown_bytes(consumer)
    assert cli(consumer, "archive") == 1
    assert "conflicting retained release" in capsys.readouterr().err.casefold()
    assert markdown_bytes(consumer) == before


def test_tag_previews_and_force_preserve_existing_ref_on_error(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Exercise tag and force through the public CLI in a disposable repository only."""
    commit(consumer, "feat: initial fixture")
    release_date = datetime.now(UTC).date().isoformat()
    (consumer / "CHANGELOG.md").write_text(f"# Changelog\n\n## [0.4.6] - {release_date}\n\n- Exact release notes.\n", encoding="utf-8")
    before = markdown_bytes(consumer)
    assert cli(consumer, "tag", "v0.4.6", "--dry-run") == 0
    assert "Exact release notes" in capsys.readouterr().out
    assert run_git_command(["tag", "--list"], cwd=consumer).stdout == ""
    assert cli(consumer, "tag", "v0.4.6") == 0
    ref = run_git_command(["rev-parse", "refs/tags/v0.4.6"], cwd=consumer).stdout
    assert run_git_command(["cat-file", "-t", "refs/tags/v0.4.6"], cwd=consumer).stdout.strip() == "tag"
    assert cli(consumer, "tag", "v0.4.6") == 1
    assert "already exists" in capsys.readouterr().err
    assert cli(consumer, "tag", "v0.5.0", "--force") == 1
    assert "does not match package version" in capsys.readouterr().err
    (consumer / "CHANGELOG.md").write_text("# Changelog\n\n## [0.4.6] - 2026-02-30\n\n- Invalid.\n", encoding="utf-8")
    assert cli(consumer, "tag", "v0.4.6", "--force") == 1
    assert run_git_command(["rev-parse", "refs/tags/v0.4.6"], cwd=consumer).stdout == ref
    (consumer / "CHANGELOG.md").write_bytes(before[Path("CHANGELOG.md")])
    assert cli(consumer, "tag", "v0.4.6", "--force") == 0
    assert markdown_bytes(consumer) == before


def test_consumer_recipes_forward_cli_arguments_and_keep_metadata_separate() -> None:
    """Check the actual merged recipes and alias, including shell-safe argument forwarding."""
    recipes = json.loads(run_safe_command("just", ["--dump", "--dump-format", "json"], cwd=REPO_ROOT).stdout)
    assert recipes["aliases"]["changelog-unreleased"]["target"] == "changelog-release"
    commands = {
        ("changelog",): "changelog generate",
        ("changelog-preview", "--tag", "v0.5.0", "--date", "2026-09-30"): 'generate --dry-run "$@"',
        ("changelog-release", "v0.5.0", "2026-09-30"): "generate --tag 'v0.5.0' --date '2026-09-30'",
        ("changelog-unreleased", "v0.5.0", "2026-09-30"): "generate --tag 'v0.5.0' --date '2026-09-30'",
        ("changelog-archive",): "changelog archive",
        ("changelog-check",): "changelog check",
        ("release-notes", "v0.4.6"): "changelog notes 'v0.4.6'",
        ("tag", "v0.4.6"): "changelog tag 'v0.4.6'",
        ("tag-force", "v0.4.6"): "changelog tag 'v0.4.6' --force",
    }
    for args, expected in commands.items():
        result = run_safe_command("just", ["--dry-run", *args], cwd=REPO_ROOT)
        assert expected in result.stderr
        assert "uv run --locked --group dev research-repo-tools" in result.stderr
        assert "update-release-version" not in result.stderr


def test_final_repository_links_and_archived_notes_remain_valid(capsys: pytest.CaptureFixture[str]) -> None:
    """Check consumer history, including references in every completed series."""
    assert cli(REPO_ROOT, "check") == 0
    root = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    for relative in re.findall(r"\]\((docs/archives/changelog/[^)]+)\)", root):
        assert (REPO_ROOT / relative).is_file()
    for tag in ("v0.1.0", "v0.2.0", "v0.3.0", "v0.4.6"):
        assert cli(REPO_ROOT, "notes", tag) == 0
        assert capsys.readouterr().out.strip()
    policy = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["rumdl"]
    assert "MD057" not in policy["disable"]
    formatter = tomllib.loads((REPO_ROOT / "changelog-rumdl.toml").read_text(encoding="utf-8"))
    assert formatter["extends"] == "pyproject.toml"
    assert formatter["global"]["extend-disable"] == ["MD057"]


def test_tag_preview_checks_citation_date_and_handles_oversized_notes(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Tag annotations validate consumer release dates and retain a full-notes link."""
    commit(consumer, "feat: initial fixture")
    released = datetime.now(UTC).date().isoformat()
    (consumer / "CHANGELOG.md").write_text(f"# Changelog\n\n## [0.4.6] - {released}\n\n- {'Long notes. ' * 12000}\n", encoding="utf-8")
    (consumer / "CITATION.cff").write_text("version: 0.4.6\ndate-released: 2000-01-01\n", encoding="utf-8")
    assert cli(consumer, "tag", "v0.4.6", "--dry-run") == 1
    assert "citation.cff date and changelog release date differ" in capsys.readouterr().err.casefold()
    (consumer / "CITATION.cff").write_text(f"version: 0.4.6\ndate-released: {released}\n", encoding="utf-8")
    before = markdown_bytes(consumer)
    assert cli(consumer, "tag", "v0.4.6", "--dry-run") == 0
    preview = capsys.readouterr().out
    assert "https://github.com/acgetchell/la-stack/blob/v0.4.6/CHANGELOG.md#" in preview
    assert len(preview.encode("utf-8")) < 125000
    assert markdown_bytes(consumer) == before
    assert run_git_command(["tag", "--list"], cwd=consumer).stdout == ""


def test_normalization_preserves_fenced_examples_and_prerelease_notes(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Rust code, Markdown links, and fenced release examples remain literal."""
    (consumer / "CHANGELOG.md").write_text(
        "# Changelog\n\n## [0.4.6-rc.1+build.7] - 2026-09-30\n\n"
        "- Keep `Vector<D>` and [PR](https://github.com/acgetchell/la-stack/pull/42).\n\n"
        "```rust\nfn value<T>() -> Vector<T> { todo!() }\n```\n\n"
        "~~~text\n## [99.9.9]\n~~~\n",
        encoding="utf-8",
    )
    assert cli(consumer, "normalize") == 0
    assert cli(consumer, "notes", "v0.4.6-rc.1+build.7") == 0
    notes = capsys.readouterr().out
    assert "`Vector<D>`" in notes
    assert "fn value<T>() -> Vector<T>" in notes
    assert "## [99.9.9]" in notes
    assert "https://github.com/acgetchell/la-stack/pull/42" in notes
    assert cli(consumer, "check") == 0


def test_formatter_failure_preserves_every_candidate(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """A consumer configuration error cannot publish a partial root or archive."""
    commit(consumer, "feat: initial fixture")
    run_git_command(["tag", "v0.4.5"], cwd=consumer)
    commit(consumer, "fix: add an unreleased change")
    (consumer / "CHANGELOG.md").write_text("# Changelog\n\n## [0.4.5] - 2026-08-21\n\n- Retained.\n", encoding="utf-8")
    manifest = consumer / "pyproject.toml"
    manifest.write_text(manifest.read_text(encoding="utf-8").replace('formatter = "changelog-rumdl.toml"', 'formatter = "missing.toml"'), encoding="utf-8")
    before = markdown_bytes(consumer)
    assert cli(consumer, "generate", "--tag", "v0.5.0", "--date", "2026-09-30") == 1
    assert "formatter configuration not found" in capsys.readouterr().err.casefold()
    assert markdown_bytes(consumer) == before
    assert not (consumer / "docs/archives/changelog").exists()
