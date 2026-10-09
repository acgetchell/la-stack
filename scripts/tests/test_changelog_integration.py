"""Consumer policy and authored history against the installed published CLI."""

import json
import shutil
import tomllib
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

from research_repo_tools.cli import main

from benchmark_process import run_git_command, run_safe_command

if TYPE_CHECKING:
    import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_published_package_and_shared_policy() -> None:
    """Both dependency owners and the registry lock select the published package."""
    manifest = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest["project"]["dependencies"] == ["research-repo-tools==0.1.8"]
    assert manifest["dependency-groups"]["tooling"] == ["research-repo-tools==0.1.8"]
    assert {"include-group": "tooling"} in manifest["dependency-groups"]["dev"]
    assert version("research-repo-tools") == "0.1.8"
    lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    package = next(item for item in lock["package"] if item["name"] == "research-repo-tools")
    assert package["version"] == "0.1.8"
    assert package["source"] == {"registry": "https://pypi.org/simple"}
    assert package["wheels"][0]["hash"] == "sha256:9db45f09a33cac4ff32ddf010d30523e295f607e24c680bcfc3d3d9786991293"
    assert package["sdist"]["hash"] == "sha256:bc855aa930ee04301dceebc8340a22576bb25712e322b51d7afaf20f5fab1ba0"
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
    run_git_command(["init", "--quiet"], cwd=tmp_path)
    for key, value in (("user.name", "Fixture"), ("user.email", "fixture@example.invalid"), ("commit.gpgsign", "false"), ("tag.gpgsign", "false")):
        run_git_command(["config", key, value], cwd=tmp_path)
    run_git_command(["commit", "--quiet", "--allow-empty", "-m", "feat: initial history"], cwd=tmp_path)
    run_git_command(["tag", "v0.4.5"], cwd=tmp_path)
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
        run_git_command(["commit", "--quiet", "--allow-empty", "-m", f"chore({scope}): bump {dependency}\n\n{body}"], cwd=tmp_path)
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
    recipes = json.loads(run_safe_command("just", ["--dump", "--dump-format", "json"], cwd=ROOT).stdout)
    assert recipes["aliases"]["changelog-unreleased"]["target"] == "changelog-release"
    result = run_safe_command("just", ["--dry-run", "changelog-preview", "--tag", "v0.4.7", "--date", "2026-10-09"], cwd=ROOT)
    assert 'research-repo-tools changelog generate --dry-run "$@"' in result.stderr
    assert main(["--root", str(ROOT), "changelog", "check"]) == 0
    for tag in ("v0.1.0", "v0.2.0", "v0.3.0", "v0.4.6"):
        assert main(["--root", str(ROOT), "changelog", "notes", tag]) == 0
        assert capsys.readouterr().out.strip()
