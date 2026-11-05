"""Consumer integration with the published Semgrep fixture checker."""

import json
import os
import shutil
from pathlib import Path

import pytest
from research_repo_tools.process import run_command

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    ("source", "diagnostic"),
    [("print('blocked')", ""), ("print('allowed')", "expected finding not reported"), ("print('blocked')\nprint('blocked')", "unexpected finding")],
)
def test_real_fixture_checker_keeps_positive_negative_gates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str, diagnostic: str) -> None:
    """Actual native scans accept a match and reject missing or extra findings."""
    for name in ("pyproject.toml", "uv.lock", ".python-version"):
        shutil.copyfile(ROOT / name, tmp_path / name)
    (tmp_path / "semgrep.yaml").write_text(
        json.dumps(
            {"rules": [{"id": "la-stack.fixture", "languages": ["python"], "message": "fixture violation", "severity": "ERROR", "pattern": "print('blocked')"}]}
        ),
        encoding="utf-8",
    )
    fixtures = tmp_path / "tests/semgrep"
    fixtures.mkdir(parents=True)
    (fixtures / "fixture.py").write_text(f"# ruleid: la-stack.fixture\n{source}\n# ok: la-stack.fixture\nprint('allowed')\n", encoding="utf-8")
    # The locked uv launch must select the project scripts ahead of ambient tools.
    ambient = tmp_path / "ambient"
    ambient.mkdir()
    for name in ("semgrep", "pysemgrep"):
        binary = ambient / (f"{name}.exe" if os.name == "nt" else name)
        binary.write_text("#!/bin/sh\nexit 97\n", encoding="utf-8", newline="\n")
        binary.chmod(0o755)
    monkeypatch.setenv("PATH", os.pathsep.join([str(ambient), os.environ["PATH"]]))
    result = run_command(
        "uv",
        ["run", "--locked", "--group", "dev", "research-repo-tools", "--root", str(tmp_path), "semgrep", "check-fixtures"],
        cwd=ROOT,
        check=False,
    )
    assert result.returncode == (1 if diagnostic else 0), result.stderr
    if diagnostic:
        assert "la-stack.fixture" in result.stderr
        assert diagnostic in result.stderr
    else:
        assert not result.stderr
