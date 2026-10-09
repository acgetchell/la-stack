"""Published maintenance CLI integration with la-stack's actual consumer policy."""

import os
import re
import shutil
import sys
import tomllib
from pathlib import Path

import pytest
from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_command

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def consumer(tmp_path: Path) -> Path:
    """Copy release inputs and real documentation, without live repository state."""
    files = run_command("git", ["--no-pager", "ls-files", "-co", "--exclude-standard", "-z"], cwd=REPO_ROOT).stdout.split("\0")
    metadata = {"Cargo.toml", "Cargo.lock", "pyproject.toml", "uv.lock", "CITATION.cff", ".python-version", "rust-toolchain.toml"}
    for name in files:
        source = REPO_ROOT / name
        if not name or not source.is_file() or (source.suffix != ".md" and name not in metadata and name != "benches/comparison/Cargo.toml"):
            continue
        destination = tmp_path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    return tmp_path


def snapshot(root: Path) -> dict[str, bytes]:
    """Capture exact source bytes for previews, failures, and retained evidence."""
    return {path.relative_to(root).as_posix(): path.read_bytes() for path in root.rglob("*") if path.is_file()}


def release_args(root: Path, *args: str) -> list[str]:
    """Prepare the next patch only in a disposable consumer."""
    current = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    major, minor, patch = current.split(".")
    return ["--root", str(root), "release", "update", f"v{major}.{minor}.{int(patch) + 1}", "--previous-release", f"v{current}", "--date", "2026-10-01", *args]


def test_release_preview_and_update_preserve_scientific_evidence_and_dependency_versions(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Exercise every real release selector while retaining measured and historical data."""
    before = snapshot(consumer)
    assert main(["--root", str(consumer), "release", "check"]) == 0
    assert main(release_args(consumer, "--dry-run")) == 0
    assert snapshot(consumer) == before
    assert main(release_args(consumer)) == 0
    after = snapshot(consumer)
    for name, contents in before.items():
        if name.startswith(("docs/archive/", "docs/archives/", "docs/performance/")) or name in {"docs/performance.md", "CHANGELOG.md"}:
            assert after[name] == contents, name
    original_cargo = tomllib.loads(before["Cargo.lock"].decode())
    updated_cargo = tomllib.loads(after["Cargo.lock"].decode())
    assert [entry for entry in updated_cargo["package"] if "source" in entry] == [entry for entry in original_cargo["package"] if "source" in entry]
    original_uv = tomllib.loads(before["uv.lock"].decode())
    updated_uv = tomllib.loads(after["uv.lock"].decode())
    assert [entry for entry in updated_uv["package"] if "registry" in entry["source"]] == [
        entry for entry in original_uv["package"] if "registry" in entry["source"]
    ]
    artifact = r"https://[^\s)]+/docs/assets/bench/[^\s)]+"
    assert re.findall(artifact, after["README.md"].decode()) == re.findall(artifact, before["README.md"].decode())
    assert after["docs/BENCHMARKING.md"] == before["docs/BENCHMARKING.md"]
    assert after["benches/comparison/Cargo.toml"] == before["benches/comparison/Cargo.toml"]
    # Metadata preparation deliberately leaves old notes intact. The release gate
    # remains red until dated notes for the target are generated.
    assert main(["--root", str(consumer), "release", "check"]) == 1
    assert "latest generated changelog release" in capsys.readouterr().err
    target = tomllib.loads(after["Cargo.toml"].decode())["package"]["version"]
    changelog = consumer / "CHANGELOG.md"
    changelog.write_text(f"# Changelog\n\n## [{target}] - 2026-10-01\n\n- Fixture release notes.\n", encoding="utf-8")
    assert main(["--root", str(consumer), "release", "check", "--final-release"]) == 0


@pytest.mark.parametrize("file", ["CITATION.cff", "README.md", "REFERENCES.md"])
def test_concept_doi_policy_rejects_drift_before_publication(consumer: Path, file: str) -> None:
    path = consumer / file
    path.write_text(path.read_text(encoding="utf-8").replace("10.5281/zenodo.18158926", "10.5281/zenodo.99999999"), encoding="utf-8")
    before = snapshot(consumer)
    assert main(["--root", str(consumer), "release", "check"]) == 1
    assert main(release_args(consumer)) == 1
    assert snapshot(consumer) == before


def test_readme_release_selector_requires_the_reviewed_link_inventory(consumer: Path) -> None:
    path = consumer / "README.md"
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("/blob/main/LICENSE", "/blob/removed/LICENSE", 1), encoding="utf-8")
    before = snapshot(consumer)
    assert main(release_args(consumer, "--dry-run")) == 1
    assert snapshot(consumer) == before


def test_update_recipe_runs_the_installed_cli_with_preview_arguments(consumer: Path) -> None:
    """Execute the imported adapter with managed tools, skipping sync and gh lookup."""
    wrapper = consumer / "recipes.just"
    wrapper.write_text(
        f"set allow-duplicate-recipes\nimport '{(REPO_ROOT / 'justfile').as_posix()}'\npython-sync:\n\n_ensure-gh:\n",
        encoding="utf-8",
    )
    before = snapshot(consumer)
    args = release_args(consumer, "--dry-run")[4:]
    result = run_command(
        "just",
        ["--justfile", str(wrapper), "--working-directory", str(consumer), "update-version", *args],
        cwd=consumer,
        env=os.environ | {"UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
    )
    assert "Cargo.toml" in result.stdout
    assert snapshot(consumer) == before


def test_markdown_recipe_uses_the_published_checker() -> None:
    """Retain the consumer's tracked-file selection and quoted path forwarding."""
    result = dry_run(REPO_ROOT, "markdown-check")
    assert 'research-repo-tools docs check-lines "${files[@]}"' in result.stderr
    assert "git ls-files -co --exclude-standard -z" in result.stderr


def test_semgrep_adapter_uses_real_consumer_rules_and_fixtures() -> None:
    """The shared checker validates the supplied namespace, including hidden workflows."""
    config = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["research-repo-tools"]["semgrep"]
    assert config["config"] == "semgrep.yaml"
    assert config["fixtures"] == "tests/semgrep"
    assert config["namespace"] == "la-stack."
    result = dry_run(REPO_ROOT, "semgrep-test")
    assert "research-repo-tools semgrep check-fixtures" in result.stderr
