"""Consumer toolchain and update integration with the installed shared release."""

import json
import os
import shlex
import shutil
import sys
import tomllib
import zipfile
from importlib.metadata import distribution, version
from pathlib import Path
from typing import TYPE_CHECKING, cast

import pytest
from packaging.requirements import Requirement
from research_repo_tools.cli import main
from research_repo_tools.dependency_policy import load_policy
from research_repo_tools.dependency_provenance import check_development_dependency
from research_repo_tools.process import run_command
from research_repo_tools.python_tools import TOOLS, versions as shared_tool_versions

if TYPE_CHECKING:
    import subprocess

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLING_REVISION = "cd8594f658b5bf3ade063a685834565cc62e0cb3"
TOOLING_REPOSITORY = "https://github.com/acgetchell/research-repo-tools"
TOOLING_REQUIREMENT = f"research-repo-tools @ git+{TOOLING_REPOSITORY}@{TOOLING_REVISION}"
CHECKED = ["run", "--locked", "--no-sync", "--no-python-downloads", "research-repo-tools", "toolchain", "run", "--"]
UPDATER = ["run", "--locked", "--only-group", "tooling", "--inexact", "research-repo-tools"]
UV_UPDATE = ["run", "--no-config", "--no-sync", "--no-python-downloads", "research-repo-tools", "deps", "update-uv"]
TOOL_UPDATE = [*UPDATER, "toolchain", "upgrade"]
SETUP = ["run", "--locked", "--managed-python", "--only-group", "tooling", "research-repo-tools", "setup"]
CARGO_UPDATE = [
    [*UPDATER, "toolchain", "run", "--", "cargo", "upgrade", "--incompatible", "allow", "--exclude", "num-bigint", "--exclude", "num-rational"],
    [*UPDATER, "toolchain", "run", "--", "cargo", "update"],
]
PYTHON_UPDATE = [
    [*UPDATER, "deps", "update-python", "--all-projects"],
]
WORKFLOWS = {
    "update": [UV_UPDATE, TOOL_UPDATE, SETUP, *CARGO_UPDATE, *PYTHON_UPDATE],
    "update-tools": [UV_UPDATE, TOOL_UPDATE, SETUP],
    "update-dependencies": [*CARGO_UPDATE, *PYTHON_UPDATE],
    "update-cargo-dependencies": CARGO_UPDATE,
    "update-python-dependencies": PYTHON_UPDATE,
    "update-cargo-tools": [TOOL_UPDATE],
    "update-uv": [UV_UPDATE],
}


def test_tooling_commit_matches_manifest_lock_and_installation() -> None:
    """Distinguish the reviewed upstream commit from its unchanged version label."""
    expected_version = "0.2.6"
    manifest = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest["project"]["dependencies"] == []
    assert manifest["tool"]["uv"]["package"] is False
    assert "build-system" not in manifest
    assert manifest["dependency-groups"]["tooling"] == [f"research-repo-tools[python-tools,semgrep]=={expected_version}"]
    assert {"include-group": "tooling"} in manifest["dependency-groups"]["dev"]
    assert manifest["tool"]["research-repo-tools"]["toolchain"]["inherit-python-tools"] is True
    assert not {Requirement(item).name for item in manifest["dependency-groups"]["dev"] if isinstance(item, str)} & TOOLS
    assert version("research-repo-tools") == expected_version
    lock = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    package = next(item for item in lock["package"] if item["name"] == "research-repo-tools")
    assert package["version"] == expected_version
    assert manifest["tool"]["uv"]["sources"]["research-repo-tools"] == {"git": TOOLING_REPOSITORY, "rev": TOOLING_REVISION}
    assert package["source"] == {"git": f"{TOOLING_REPOSITORY}?rev={TOOLING_REVISION}#{TOOLING_REVISION}"}
    # Distribution.read_text decodes UTF-8 internally and has no encoding argument.
    installed = json.loads(
        # nosemgrep: la-stack.python.explicit-path-text-encoding-in-tests
        distribution("research-repo-tools").read_text("direct_url.json") or "{}"
    )
    assert installed["url"] == TOOLING_REPOSITORY
    assert installed["vcs_info"] == {"vcs": "git", "commit_id": TOOLING_REVISION, "requested_revision": TOOLING_REVISION}
    locked = {item["name"]: item["version"] for item in lock["package"]}
    for name, expected in shared_tool_versions().items():
        assert locked[name] == version(name) == expected
    identity = check_development_dependency(REPO_ROOT)
    assert identity is not None
    assert (identity.name, identity.version, identity.revision) == ("research-repo-tools", expected_version, TOOLING_REVISION)


def test_release_dependency_policy_covers_tooling_and_private_comparison_workspace() -> None:
    policy = load_policy(REPO_ROOT)
    assert policy is not None
    assert {(project.ecosystem, project.directory, project.owned) for project in policy.projects} == {
        ("python", ".", (".",)),
        ("cargo", ".", (".", "benches/comparison")),
    }
    manifest = tomllib.loads((REPO_ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    assert manifest["workspace"]["members"] == ["benches/comparison"]
    comparison = tomllib.loads((REPO_ROOT / "benches/comparison/Cargo.toml").read_text(encoding="utf-8"))
    assert comparison["package"]["publish"] is False


@pytest.mark.parametrize(
    "command",
    [
        ["release", "dependencies"],
        ["performance", "prepare", "tooling/performance-release.toml", "--release-intent"],
        ["performance", "release", "tooling/performance-release.toml", "--release-intent"],
    ],
)
def test_development_pin_cannot_enter_release_work(command: list[str], tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """The installed public gate rejects before reading plans or creating evidence."""
    for name in ("pyproject.toml", "uv.lock", "dependency-sources.toml"):
        shutil.copyfile(REPO_ROOT / name, tmp_path / name)
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    assert main(["--root", str(tmp_path), *command]) == 1
    diagnostic = capsys.readouterr().err
    assert "release requires a published registry dependency without tool.uv.sources overrides" in diagnostic
    assert "Publish the supplying upstream version" in diagnostic
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir()} == before


def test_isolated_release_runtime_uses_the_same_commit_outside_checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Exercise the actual trusted runtime on each native CI platform."""
    expected = version("research-repo-tools")
    source = (REPO_ROOT / "justfile").read_text(encoding="utf-8")
    assert source.count(f'[script("uvx", "--python", "3.15", "--from", "{TOOLING_REQUIREMENT}", "python", "-I")]') == 3
    (tmp_path / "research_repo_tools.py").write_text('raise RuntimeError("untrusted checkout import")\n', encoding="utf-8")
    (tmp_path / "tooling").mkdir()
    (tmp_path / "justfile").write_bytes((REPO_ROOT / "justfile").read_bytes())
    (tmp_path / "tooling/release_candidate.py").write_text(
        "import json, sys, research_repo_tools\nfrom importlib.metadata import distribution\n"
        'd = distribution("research-repo-tools")\n'
        'print(json.dumps([d.version, research_repo_tools.__file__, json.loads(d.read_text("direct_url.json")), sys.argv[1:]]))\n',
        encoding="utf-8",
        newline="\n",
    )
    monkeypatch.setenv("PYTHONPATH", str(tmp_path))
    result = run_command(
        "just",
        [
            "--no-deps",
            "--justfile",
            str(REPO_ROOT / "justfile"),
            "--working-directory",
            str(tmp_path),
            "release-receive",
            str(tmp_path / "output with spaces"),
            "--incoming",
            str(tmp_path / "input with spaces"),
        ],
        timeout=180,
        cwd=tmp_path,
    )
    installed_version, location, provenance, arguments = json.loads(result.stdout)
    assert installed_version == expected
    assert provenance["vcs_info"]["commit_id"] == TOOLING_REVISION
    assert not Path(location).resolve().is_relative_to(REPO_ROOT)
    assert not Path(location).resolve().is_relative_to(tmp_path)
    assert arguments == ["receive", str(tmp_path / "output with spaces"), "--incoming", str(tmp_path / "input with spaces")]


@pytest.fixture
def consumer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Use the actual recipes and declarations; record updater process boundaries."""
    root = tmp_path / "consumer with spaces"
    root.mkdir()
    for name in ("pyproject.toml", "uv.lock", "Cargo.toml", "Cargo.lock", "rust-toolchain.toml", ".python-version"):
        shutil.copyfile(REPO_ROOT / name, root / name)
    (root / "recipes.just").write_text(
        f"set allow-duplicate-recipes\nimport '{(REPO_ROOT / 'justfile').as_posix()}'\n_ensure-gh:\n\n_ensure-jq:\n",
        encoding="utf-8",
    )
    recorder = root / "record calls.py"
    recorder.write_text(
        "import json, os, pathlib, sys\n"
        "log = pathlib.Path('calls.jsonl')\n"
        "with log.open('a', encoding='utf-8') as stream: stream.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if len(log.read_text(encoding='utf-8').splitlines()) == int(os.environ.get('UPDATE_FAIL_STEP', '0')):\n"
        "    print('fixture update failed', file=sys.stderr)\n"
        "    sys.exit(23)\n",
        encoding="utf-8",
    )
    shim = root / "uv"
    shim.write_text(f'#!/bin/sh\nexec {shlex.quote(Path(sys.executable).as_posix())} {shlex.quote(recorder.as_posix())} "$@"\n', encoding="utf-8", newline="\n")
    shim.chmod(0o755)
    monkeypatch.setenv("PATH", f"{root}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.delenv("UPDATE_FAIL_STEP", raising=False)
    return root


def invoke(root: Path, recipe: str) -> tuple[subprocess.CompletedProcess[str], list[list[str]]]:
    """Run Just's actual dependency sequencing without live upgrades."""
    result = run_command("just", ["--justfile", str(root / "recipes.just"), "--working-directory", str(root), recipe], cwd=root, check=False)
    calls = [cast("list[str]", json.loads(line)) for line in (root / "calls.jsonl").read_text(encoding="utf-8").splitlines()]
    return result, calls


@pytest.mark.parametrize("recipe", WORKFLOWS)
def test_update_recipes_sequence_and_preserve_coupled_exclusions(consumer: Path, recipe: str) -> None:
    """Direct and aggregate updates share boundaries, order, exclusions, and dev sync."""
    result, calls = invoke(consumer, recipe)
    assert result.returncode == 0, result.stderr
    assert calls == WORKFLOWS[recipe]


@pytest.mark.parametrize("step", range(1, len(WORKFLOWS["update"]) + 1))
def test_update_failure_stops_later_steps(consumer: Path, monkeypatch: pytest.MonkeyPatch, step: int) -> None:
    monkeypatch.setenv("UPDATE_FAIL_STEP", str(step))
    result, calls = invoke(consumer, "update")
    assert result.returncode == 23
    assert "fixture update failed" in result.stderr
    assert calls == WORKFLOWS["update"][:step]


@pytest.mark.parametrize("workers", [None, "0", "2"])
@pytest.mark.parametrize("failure", [False, True])
def test_python_recipe_preserves_worker_selection_and_failure(consumer: Path, monkeypatch: pytest.MonkeyPatch, workers: str | None, failure: bool) -> None:
    if workers is None:
        monkeypatch.delenv("LA_STACK_TEST_WORKERS", raising=False)
    else:
        monkeypatch.setenv("LA_STACK_TEST_WORKERS", workers)
    if failure:
        monkeypatch.setenv("UPDATE_FAIL_STEP", "4")

    result, calls = invoke(consumer, "test-python")

    assert result.returncode == (23 if failure else 0), result.stderr
    assert len(calls) == 4
    invocation = calls[-1]
    assert invocation[invocation.index("pytest") :] == [
        "pytest",
        "-q",
        "-n",
        workers or "4",
        "--dist",
        "worksteal",
        "--max-worker-restart",
        "0",
        "--durations=20",
    ]
    if failure:
        assert "fixture update failed" in result.stderr


def test_declared_tools_are_verified_and_cargo_edit_runs_from_managed_store(capsys: pytest.CaptureFixture[str]) -> None:
    """Exercise the installed public CLI on each native CI platform, without installs."""
    assert main(["--root", str(REPO_ROOT), "toolchain", "check", "--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert all(status["ok"] for status in report)
    result = run_command("uv", [*CHECKED, "cargo", "upgrade", "--version"], cwd=REPO_ROOT)
    declared = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["research-repo-tools"]["toolchain"]["cargo"]
    names = {"cargo-edit": "cargo-edit-upgrade", "taplo-cli": "taplo"}
    for package, expected in declared.items():
        status = next(item for item in report if item["name"] == names.get(package, package))
        assert status["actual"] == status["required"] == expected
        directory = Path(status["path"]).parent
        assert directory.name == "bin"
        assert directory.parent.name == expected
        assert directory.parent.parent.name == package
    assert result.stdout.strip() == f"cargo-edit-upgrade {declared['cargo-edit']}"
    # The checked runner supplies the same managed Rust to subprocesses, including native Python builds.
    result = run_command("uv", [*CHECKED, "uv", "run", "--locked", "--no-sync", "python", "-c", "import shutil; print(shutil.which('cargo'))"], cwd=REPO_ROOT)
    rustup = next(status for status in report if status["name"] == "rustup")
    assert Path(result.stdout.strip()).parent == Path(rustup["path"]).parent


def fixture_wheel(registry: Path, requirement: str) -> None:
    """Supply offline resolution candidates; installed CLI code comes from the host."""
    parsed = Requirement(requirement)
    pin = next(iter(parsed.specifier))
    assert pin.operator == "=="
    name, version = parsed.name, pin.version
    distribution = name.replace("-", "_")
    info = f"{distribution}-{version}.dist-info"
    extra = ""
    if name == "research-repo-tools":
        extra = "Provides-Extra: python-tools\n" + "".join(
            f"Requires-Dist: {tool}=={expected}; extra == 'python-tools'\n" for tool, expected in shared_tool_versions().items()
        )
    with zipfile.ZipFile(registry / f"{distribution}-{version}-py3-none-any.whl", "w") as archive:
        archive.writestr(f"{info}/METADATA", f"Metadata-Version: 2.3\nName: {name}\nVersion: {version}\nRequires-Python: >=3.15\n{extra}")
        archive.writestr(f"{info}/WHEEL", "Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n")
        archive.writestr(f"{info}/RECORD", "")


def test_python_update_preserves_shared_pin_and_tools_with_explicit_dev_sync(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Run the published updater and real uv against consumer config and offline wheels."""
    root = tmp_path / "offline consumer"
    root.mkdir()
    registry = tmp_path / "wheels"
    registry.mkdir()
    original = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    # This fixture tests registry updates against inert local wheels. The
    # development Git source has separate identity coverage above and cannot
    # resolve in the fixture's deliberately empty offline cache.
    override = f'[tool.uv.sources]\nresearch-repo-tools = {{ git = "{TOOLING_REPOSITORY}", rev = "{TOOLING_REVISION}" }}\n'
    assert original.count(override) == 1
    source = original.replace(override, "").replace("package = false", "package = false\ndefault-groups = []")
    manifest = root / "pyproject.toml"
    manifest.write_text(source, encoding="utf-8")
    before = tomllib.loads(source)
    groups = before["dependency-groups"]
    shared = shared_tool_versions()
    for requirement in [
        *groups["tooling"],
        *(entry for entry in groups["dev"] if isinstance(entry, str)),
        *(f"{tool}=={expected}" for tool, expected in shared.items()),
        "yamllint==99.0.0",
        "ruff==99.0.0",
        "ty==99.0.0",
    ]:
        fixture_wheel(registry, requirement)
    for name in ("UV_PROJECT", "UV_WORKING_DIR", "UV_WORKING_DIRECTORY", "UV_CONFIG_FILE"):
        monkeypatch.delenv(name, raising=False)
    for name, value in {
        "UV_CACHE_DIR": str(tmp_path / "cache"),
        "UV_FIND_LINKS": str(registry),
        "UV_NO_INDEX": "1",
        "UV_OFFLINE": "1",
        "UV_PYTHON": sys.executable,
        "UV_PYTHON_DOWNLOADS": "never",
        "UV_PROJECT_ENVIRONMENT": str(root / ".venv"),
    }.items():
        monkeypatch.setenv(name, value)

    run_command("uv", ["lock"], cwd=root)
    assert main(["--root", str(root), "deps", "update-python"]) == 0
    after = tomllib.loads(manifest.read_text(encoding="utf-8"))
    assert "yamllint==99.0.0" in after["dependency-groups"]["dev"]
    assert after["dependency-groups"]["tooling"] == groups["tooling"]
    assert after["tool"]["research-repo-tools"] == before["tool"]["research-repo-tools"]
    assert after["tool"]["uv"] == before["tool"]["uv"]
    run_command("uv", ["lock", "--upgrade"], cwd=root)
    lock = tomllib.loads((root / "uv.lock").read_text(encoding="utf-8"))
    packages = {package["name"]: package["version"] for package in lock["package"]}
    requirement = Requirement(groups["tooling"][0])
    assert packages[requirement.name] in requirement.specifier
    assert packages["yamllint"] == "99.0.0"
    for name, expected in shared.items():
        assert packages[name] == expected
    run_command("uv", ["sync", "--locked", "--group", "dev"], cwd=root)
    installed = run_command("uv", ["run", "--locked", "--no-sync", "python", "-c", "from importlib.metadata import version; print(version('ruff'))"], cwd=root)
    assert installed.stdout.strip() == shared["ruff"]
