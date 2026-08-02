"""Regression tests for the Just recipe surface."""

import json
import os
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Any

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
