"""Consumer toolchain and update integration with the installed shared release."""

import json
import os
import shlex
import shutil
import sys
import tomllib
import zipfile
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, cast

import pytest
from packaging.requirements import Requirement
from research_repo_tools.cli import main
from research_repo_tools.process import run_command
from research_repo_tools.python_tools import TOOLS, versions as shared_tool_versions

if TYPE_CHECKING:
    import subprocess

REPO_ROOT = Path(__file__).resolve().parents[2]
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
    "update-python-deps": PYTHON_UPDATE,
    "update-cargo-tools": [TOOL_UPDATE],
    "update-uv": [UV_UPDATE],
}


def test_published_tooling_package_matches_manifest_lock_and_installation() -> None:
    """Keep the exact published-package contract in one consumer test."""
    expected_version = "0.2.3"
    manifest = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest["project"]["dependencies"] == []
    assert manifest["tool"]["uv"]["package"] is False
    assert "build-system" not in manifest
    assert manifest["dependency-groups"]["tooling"] == [f"research-repo-tools[python-tools]=={expected_version}"]
    assert {"include-group": "tooling"} in manifest["dependency-groups"]["dev"]
    assert manifest["tool"]["research-repo-tools"]["toolchain"]["inherit-python-tools"] is True
    assert not {Requirement(item).name for item in manifest["dependency-groups"]["dev"] if isinstance(item, str)} & TOOLS
    assert version("research-repo-tools") == expected_version
    lock = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    package = next(item for item in lock["package"] if item["name"] == "research-repo-tools")
    assert package["version"] == expected_version
    assert package["source"] == {"registry": "https://pypi.org/simple"}
    assert package["wheels"][0]["hash"] == "sha256:da5f0259027662e7706b7f2b8a38d8e1f1d8a5dce314e2c996177b2947514fc6"
    assert package["sdist"]["hash"] == "sha256:13700b174656ffad265ecb926dd66013949853b7a39c088ab6d3852aba72df95"
    locked = {item["name"]: item["version"] for item in lock["package"]}
    for name, expected in shared_tool_versions().items():
        assert locked[name] == version(name) == expected


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
    # Empty fixture wheels represent resolver candidates, not executable tools.
    source = original.replace("package = false", "package = false\ndefault-groups = []")
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
