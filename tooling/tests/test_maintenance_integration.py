"""Published maintenance CLI integration with la-stack's actual consumer policy."""

import os
import re
import shlex
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
        if name.startswith(("docs/archive/", "docs/archives/", "docs/performance/", "docs/performance-runs/")) or name in {
            "docs/performance.md",
            "CHANGELOG.md",
        }:
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


def test_release_updates_preserve_new_complete_evidence_links_and_reports(consumer: Path) -> None:
    """New publication links retain their measured tag through later metadata updates."""
    readme = consumer / "README.md"
    links = [f"https://github.com/acgetchell/la-stack/blob/v0.5.0/docs/performance-runs/runs/frozen/{name}.json" for name in ("run", "evidence")]
    with readme.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write("\n" + "\n".join(f"[Retained {index}]({link})" for index, link in enumerate(links)) + "\n")
    report = consumer / "docs/performance-runs/frozen.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_bytes(b"# la-stack v0.5.0 measured report\n")
    assert main(release_args(consumer)) == 0
    assert all(link in readme.read_text(encoding="utf-8") for link in links)
    assert report.read_bytes() == b"# la-stack v0.5.0 measured report\n"


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


def test_markdown_recipe_preserves_active_and_historical_selection(consumer: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """The real shared selector includes new guides and preserves measured history."""
    run_command("git", ["init", "--quiet"], cwd=consumer)
    added = consumer / "docs/new guide.md"
    added.write_text("# New guide\n", encoding="utf-8")
    selection = shlex.split(run_command("just", ["--evaluate", "_markdown_inputs"], cwd=REPO_ROOT).stdout)
    assert main(["--root", str(consumer), "files", "list", *selection]) == 0
    selected = set(capsys.readouterr().out.splitlines())
    assert {"README.md", "docs/new guide.md", "docs/RELEASING.md", "tooling/README.md"} <= selected
    assert not any(name.startswith(("docs/archive/", "docs/archives/changelog/", "docs/performance-runs/")) for name in selected)
    assert not {"CHANGELOG.md", "docs/performance.md"} & selected
    checked = dry_run(REPO_ROOT, "markdown-check").stderr
    fixed = dry_run(REPO_ROOT, "markdown-fix").stderr
    assert "research-repo-tools files run" in checked
    assert "research-repo-tools files run" in fixed
    assert "--exclude README.md -- research-repo-tools docs check-lines" in checked
    assert "-- rumdl check" in checked


def test_semgrep_adapter_uses_real_consumer_rules_and_fixtures() -> None:
    """The shared checker validates the supplied namespace, including hidden workflows."""
    config = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["research-repo-tools"]["semgrep"]
    assert config["config"] == "semgrep.yaml"
    assert config["fixtures"] == "tests/semgrep"
    assert config["namespace"] == "la-stack."
    result = dry_run(REPO_ROOT, "semgrep-test")
    assert "research-repo-tools semgrep check-fixtures" in result.stderr


def test_semgrep_preserves_its_pin_and_mcp_override_in_python_314() -> None:
    """The isolated native scanner retains its exact pin and patched MCP API."""
    root = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scanner = tomllib.loads((REPO_ROOT / "tooling/semgrep/pyproject.toml").read_text(encoding="utf-8"))
    assert not any(isinstance(item, str) and item.startswith("semgrep") for item in root["dependency-groups"]["dev"])
    assert scanner["project"]["dependencies"] == ["semgrep==1.182.0"]
    assert scanner["project"]["requires-python"] == ">=3.14,<3.15"
    assert scanner["tool"]["uv"]["override-dependencies"] == ["mcp>=1.30.0,<2"]
    result = run_command("uv", ["run", "--locked", "--project", "tooling/semgrep", "semgrep", "--disable-version-check", "--version"], cwd=REPO_ROOT)
    assert result.stdout.strip() == "1.182.0"


def test_python_update_includes_the_isolated_semgrep_project() -> None:
    """The shared updater owns both locks while each project owns its constraints."""
    root = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert root["tool"]["research-repo-tools"]["deps"]["python-projects"] == ["tooling/semgrep"]
    result = dry_run(REPO_ROOT, "update-python-dependencies")
    assert "research-repo-tools deps update-python --all-projects" in result.stderr
    assert "uv lock --upgrade" not in result.stderr


@pytest.mark.parametrize("output_args", [pytest.param((), id="local"), pytest.param(("--sarif",), id="sarif")])
def test_semgrep_production_targets_include_tests_and_exclude_rule_fixtures(tmp_path: Path, output_args: tuple[str, ...]) -> None:
    """Exercise both scan formats without Git ignores masking Semgrep exclusions."""
    shutil.copyfile(REPO_ROOT / ".semgrepignore", tmp_path / ".semgrepignore")
    included = {
        "tests/proptest_exact.rs",
        "benches/comparison/tests/vs_linalg_inputs.rs",
        "tooling/tests/test_release_baseline.py",
        ".github/workflows/semgrep-sarif.yml",
        "test/control.rs",
        "testsuite/control.py",
        "src/control.go",
        "assets/control.js",
    }
    excluded = {
        f"{directory}/ignored.py"
        for directory in (
            ".svn",
            ".hg",
            "_darcs",
            "CVS",
            "build",
            "vendor",
            "dist",
            ".env",
            ".tox",
            "node_modules",
            ".npm",
            ".yarn",
            ".venv",
            "_opam",
            "_build",
            "_cargo",
            "tests/semgrep",
            "nested/dist",
            "nested/.env",
            "nested/_build",
        )
    } | {"assets/bundle.min.js", "src/example_test.go", "nested/example_test.go"}
    for name in included | excluded:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        comment = "#" if path.suffix in {".py", ".yml"} else "//"
        path.write_text(f"{comment} Target selection fixture.\n", encoding="utf-8")
    result = run_command(
        "uv",
        [
            "run",
            "--locked",
            "--project",
            str(REPO_ROOT / "tooling/semgrep"),
            "semgrep",
            "scan",
            "--metrics",
            "off",
            "--disable-version-check",
            "--error",
            "--strict",
            "--timeout",
            "30",
            "--config",
            str(REPO_ROOT / "semgrep.yaml"),
            *output_args,
            "--x-ls-long",
            ".",
        ],
        cwd=tmp_path,
    )
    selected = {Path(line.removeprefix("selected ")).as_posix() for line in result.stdout.splitlines() if line.startswith("selected ")}
    assert included <= selected
    assert not (excluded & selected)
