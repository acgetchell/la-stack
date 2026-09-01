"""Regression tests for the Just recipe surface."""

import json
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Any

import update_cargo_tool_pins

REPO_ROOT = Path(__file__).resolve().parents[2]


def run_just(
    *args: str,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run the repository's installed Just executable without a shell."""
    executable = shutil.which("just")
    assert executable is not None
    return subprocess.run(  # noqa: S603 - executable is resolved; arguments are fixed by tests.
        [executable, *args],
        cwd=REPO_ROOT,
        check=check,
        capture_output=True,
        encoding="utf-8",
        env=env,
    )


def just_recipes() -> dict[str, dict[str, Any]]:
    """Return parsed recipe metadata from the pinned Just executable."""
    result = run_just("--dump", "--dump-format", "json")
    recipes = json.loads(result.stdout)["recipes"]
    assert isinstance(recipes, dict)
    return recipes


def test_uv_backed_helpers_reuse_pinned_guard() -> None:
    """Local uv consumers should share one exact-version implementation."""
    recipes = just_recipes()
    ensure_uv_body = json.dumps(recipes["_ensure-uv"]["body"])
    setup_tools_body = json.dumps(recipes["setup-tools"]["body"])

    assert "uv --version" in ensure_uv_body
    assert "uv_version" in ensure_uv_body
    assert [dependency["recipe"] for dependency in recipes["_ensure-uv"]["dependencies"]] == ["_ensure-uv-available"]
    assert "verify_tool_version uv" in setup_tools_body
    for name in ("_ensure-actionlint", "_ensure-shellcheck", "_ensure-shfmt", "_ensure-yamllint"):
        dependencies = {dependency["recipe"] for dependency in recipes[name]["dependencies"]}
        assert "_ensure-uv" in dependencies, name


def test_uv_guard_reports_expected_and_actual_versions(tmp_path: Path) -> None:
    """A mismatched uv executable should fail with actionable version details."""
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_uv = fake_bin / "uv"
    fake_uv.write_text("#!/bin/sh\nprintf '%s\\n' 'uv 9.9.9'\n", encoding="utf-8")
    fake_uv.chmod(fake_uv.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    environment = os.environ.copy()
    environment["CARGO_HOME"] = str(tmp_path / "cargo")
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"

    expected = run_just("--evaluate", "uv_version").stdout.strip()
    result = run_just("_ensure-uv", check=False, env=environment)

    assert result.returncode != 0
    assert f"version '9.9.9', expected '{expected}'" in result.stderr


def test_update_workflow_composes_scoped_dependency_and_tool_updates() -> None:
    """Update recipes should cover repo state without touching unrelated global tools."""
    recipes = just_recipes()
    update_dependencies = [dependency["recipe"] for dependency in recipes["update"]["dependencies"]]

    assert update_dependencies == ["_ensure-cargo-install-update", "update-dependencies", "update-cargo-tools"]

    aggregate_result = run_just("--dry-run", "update")
    aggregate_update = aggregate_result.stdout + aggregate_result.stderr
    cargo_upgrade = "cargo upgrade --incompatible allow --exclude num-bigint --exclude num-rational"
    assert aggregate_update.index("command -v cargo-install-update") < aggregate_update.index(cargo_upgrade)

    dependency_result = run_just("--dry-run", "update-dependencies")
    dependency_update = dependency_result.stdout + dependency_result.stderr
    dependency_preflights = [dependency["recipe"] for dependency in recipes["update-dependencies"]["dependencies"]]
    assert dependency_preflights[:2] == ["_ensure-cargo-edit", "_ensure-uv-available"]
    assert dependency_update.index("cargo upgrade --version") < dependency_update.index(cargo_upgrade)
    assert dependency_update.index("uv --version") < dependency_update.index(cargo_upgrade)
    assert cargo_upgrade in dependency_update
    assert "cargo update" in dependency_update
    assert "uv run --locked update-python-dev-pins" in dependency_update
    assert "uv lock --upgrade" in dependency_update
    assert dependency_update.index("uv run --locked update-python-dev-pins") < dependency_update.index("uv lock --upgrade")
    assert "uv sync --locked --group dev" in dependency_update
    assert "cargo install-update --all" not in dependency_update
    assert "uv tool upgrade" not in dependency_update

    tool_result = run_just("--dry-run", "update-cargo-tools")
    tool_update = tool_result.stdout + tool_result.stderr
    assert "command -v cargo-install-update" in tool_update
    assert "cargo install-update --locked" in tool_update
    assert "update-cargo-tool-pins" in tool_update
    assert "cargo install-update --all" not in tool_update
    assert "uv tool upgrade" not in tool_update
    package_block = re.search(r"packages=\(\n(?P<packages>.*?)\n\)", tool_update, re.DOTALL)
    assert package_block is not None
    updated_packages = set(re.findall(r"^\s+([a-z0-9-]+)$", package_block.group("packages"), re.MULTILINE))
    assert updated_packages == set(update_cargo_tool_pins.PIN_TO_PACKAGE.values())


def test_setup_tools_installs_and_verifies_cargo_update_provider() -> None:
    """A clean setup must provide the updater used by the update workflow."""
    body = json.dumps(just_recipes()["setup-tools"]["body"])

    assert "cargo install --locked cargo-update --version" in body
    assert "verify_tool_version cargo-install-update" in body


def test_managed_tool_pins_exist_once_in_root_justfile() -> None:
    """Every managed Cargo package and uv should map to one root Just pin."""
    justfile_text = (REPO_ROOT / "justfile").read_text(encoding="utf-8")

    for pin in update_cargo_tool_pins.PIN_TO_TOOL:
        assert len(re.findall(rf'(?m)^{re.escape(pin)}\s*:=\s*"[^"]+"\s*$', justfile_text)) == 1
