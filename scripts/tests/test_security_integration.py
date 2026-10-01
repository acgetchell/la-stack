"""Check that the consumer Gitleaks exception preserves secret detection."""

import hashlib
from pathlib import Path

import pytest

from benchmark_process import run_safe_command

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    ("relative", "credential", "status"),
    [
        ("scripts/bench_compare.py", False, 0),
        ("scripts/bench_compare.py", True, 1),
        ("scripts/other.py", False, 1),
    ],
)
def test_historical_alias_exception_is_limited_to_its_assignment_and_path(tmp_path: Path, relative: str, *, credential: bool, status: int) -> None:
    """Native scans still catch a synthetic key and the alias outside its owner."""
    source = tmp_path / relative
    source.parent.mkdir(parents=True, exist_ok=True)
    identifier = "V0_4_3_API_COMPATIBILITY"
    contents = f"_{identifier} = {identifier}\n"
    if credential:
        token = hashlib.sha256(b"public Gitleaks regression fixture").hexdigest()
        contents += "api_" + f'key = "{token}"\n'
    source.write_text(contents, encoding="utf-8")
    result = run_safe_command(
        "uv",
        [
            "run",
            "--locked",
            "--no-sync",
            "research-repo-tools",
            "toolchain",
            "run",
            "--",
            "gitleaks",
            "dir",
            "--config",
            str(REPO_ROOT / ".gitleaks.toml"),
            "--redact",
            "--ignore-gitleaks-allow",
            str(tmp_path),
        ],
        cwd=REPO_ROOT,
        check=False,
    )
    assert result.returncode == status
