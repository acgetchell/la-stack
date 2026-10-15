"""Consumer policy and authored history against the installed published CLI."""

import shutil
import tomllib
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run, inspect_justfile
from research_repo_tools.process import run_command

if TYPE_CHECKING:
    import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_published_package_and_shared_policy() -> None:
    """The dependency-only environment and registry lock select the published package."""
    manifest = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest["project"]["dependencies"] == []
    assert manifest["tool"]["uv"]["package"] is False
    assert "build-system" not in manifest
    assert manifest["dependency-groups"]["tooling"] == ["research-repo-tools==0.1.9"]
    assert {"include-group": "tooling"} in manifest["dependency-groups"]["dev"]
    assert version("research-repo-tools") == "0.1.9"
    lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    package = next(item for item in lock["package"] if item["name"] == "research-repo-tools")
    assert package["version"] == "0.1.9"
    assert package["source"] == {"registry": "https://pypi.org/simple"}
    assert package["wheels"][0]["hash"] == "sha256:c4805004e9a857c072e7b735b033dd13ad256fb283c6e1a4adca0aa59cb67540"
    assert package["sdist"]["hash"] == "sha256:7bfb16d7a38e8b3cecd30ce526cf808b377b4480ebd512740198a0361c3d027b"
    assert manifest["tool"]["research-repo-tools"]["changelog"] == {
        "owner": "acgetchell",
        "repository": "la-stack",
        "formatter": "changelog-rumdl.toml",
        "dependency-bodies": "preserve",
    }
    assert not (ROOT / "cliff.toml").exists()


def test_authored_dependency_notes_survive_generation(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Exercise the consumer policy with Ruff, Ty and setuptools notes and code."""
    for name in ("pyproject.toml", "changelog-rumdl.toml", ".python-version", "rust-toolchain.toml"):
        shutil.copyfile(ROOT / name, tmp_path / name)
    run_command("git", ["--no-pager", "init", "--quiet"], cwd=tmp_path)
    for key, value in (("user.name", "Fixture"), ("user.email", "fixture@example.invalid"), ("commit.gpgsign", "false"), ("tag.gpgsign", "false")):
        run_command("git", ["--no-pager", "config", key, value], cwd=tmp_path)
    run_command("git", ["--no-pager", "commit", "--quiet", "--allow-empty", "-m", "feat: initial history"], cwd=tmp_path)
    run_command("git", ["--no-pager", "tag", "v0.4.5"], cwd=tmp_path)
    notes = (
        ("ruff", "deps-dev", "Read [release notes](https://github.com/astral-sh/ruff/compare/0.16.1...0.16.2) before changing `Vector<D>` checks."),
        ("ty", "deps-dev", "Keep [changelog](https://github.com/astral-sh/ty/releases) context and `--error all`."),
        (
            "setuptools",
            "deps",
            'See [comparison](https://github.com/pypa/setuptools/compare/v83.0.0...v84.0.0).\n\n```toml\nrequires = ["setuptools>=84.0.0"]\n```',
        ),
    )
    for dependency, scope, body in notes:
        run_command("git", ["--no-pager", "commit", "--quiet", "--allow-empty", "-m", f"chore({scope}): bump {dependency}\n\n{body}"], cwd=tmp_path)
    (tmp_path / "CHANGELOG.md").write_text("# Changelog\n\n## [0.4.5] - 2000-01-02\n\n- Historical date.\n", encoding="utf-8")
    args = ["--root", str(tmp_path), "changelog", "generate", "--tag", "v0.4.6", "--date", "2026-10-09"]
    assert main([*args, "--dry-run"]) == 0
    preview = capsys.readouterr().out
    assert "### Dependencies" in preview
    assert "## [0.4.5] - 2000-01-02" in preview
    for expected in (
        "https://github.com/astral-sh/ruff/compare/0.16.1...0.16.2",
        "https://github.com/astral-sh/ty/releases",
        "https://github.com/pypa/setuptools/compare/v83.0.0...v84.0.0",
        "before changing `Vector<D>` checks.",
        "context and `--error all`.",
        'requires = ["setuptools>=84.0.0"]',
    ):
        assert expected in preview
    assert main(args) == 0
    generated = (tmp_path / "CHANGELOG.md").read_bytes()
    assert generated == preview.encode()
    assert main(args) == 0
    assert (tmp_path / "CHANGELOG.md").read_bytes() == generated


def test_consumer_recipes_and_retained_history(capsys: pytest.CaptureFixture[str]) -> None:
    """Native recipes delegate generation and all historical notes remain readable."""
    recipes = inspect_justfile(ROOT)
    assert recipes.aliases["changelog-unreleased"] == "changelog-release"
    result = dry_run(ROOT, "changelog-preview", ["--tag", "v0.4.7", "--date", "2026-10-09"])
    assert 'research-repo-tools changelog generate --dry-run "$@"' in result.stderr
    assert main(["--root", str(ROOT), "changelog", "check"]) == 0
    for tag in ("v0.1.0", "v0.2.0", "v0.3.0", "v0.4.6"):
        assert main(["--root", str(ROOT), "changelog", "notes", tag]) == 0
        assert capsys.readouterr().out.strip()
