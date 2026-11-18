"""Check that the consumer Gitleaks exception preserves secret detection."""

import base64
import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from research_repo_tools.process import run_command

if TYPE_CHECKING:
    import subprocess

REPO_ROOT = Path(__file__).resolve().parents[2]
FROZEN_RUN = "docs/performance-runs/runs/2e55909621d0aacbfa7670f678574ae750311ba29353e6e9ebb11461414876a2/run.json"
V054_RUN = "docs/performance-runs/runs/c71813e122e5fde26d26511b6a57ca21e5a24ed1b03bce2905029cf518ffd29e/run.json"
LEAK_FOUND = 3


@pytest.mark.parametrize(
    ("relative", "credential", "status"),
    [
        ("scripts/bench_compare.py", False, 0),
        ("scripts/bench_compare.py", True, LEAK_FOUND),
        ("scripts/other.py", False, LEAK_FOUND),
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
    result = scan_secrets(tmp_path)
    assert result.returncode == status, result.stderr


@pytest.mark.parametrize("representation", ["direct", "chunk", "payload"])
@pytest.mark.parametrize(
    ("relative", "credential", "status"),
    [
        (FROZEN_RUN, False, 0),
        (FROZEN_RUN, True, LEAK_FOUND),
        (V054_RUN, False, 0),
        (V054_RUN, True, LEAK_FOUND),
        ("docs/performance-runs/other.json", False, LEAK_FOUND),
    ],
)
def test_frozen_statistics_exception_preserves_other_secret_detection(
    tmp_path: Path, representation: str, relative: str, *, credential: bool, status: int
) -> None:
    """Only the frozen runs' encoded numerical fields are exempted."""
    fences = b"[23.01652962024839,23.323939759650607,24.143700131389842,24.451110270792054]"
    encoded = base64.b64encode(fences).decode("ascii")
    record = {}
    if credential:
        record["api_" + "key"] = hashlib.sha256(b"public Gitleaks regression fixture").hexdigest()
    record["tukey.json"] = encoded
    if representation == "payload":
        payload = json.dumps({"raw_files": record}).encode("utf-8")
        record = {"payload": base64.b64encode(payload).decode("ascii")}
    source = tmp_path / relative
    source.parent.mkdir(parents=True, exist_ok=True)
    contents = json.dumps(record, indent=2) + "\n"
    if representation == "chunk":
        contents = contents.removesuffix('"\n}\n')
    source.write_text(contents, encoding="utf-8", newline="\n")
    result = scan_secrets(tmp_path)
    assert result.returncode == status, result.stderr


@pytest.mark.parametrize(
    ("relative", "frozen", "status"),
    [(V054_RUN, True, 0), (V054_RUN, False, LEAK_FOUND), ("docs/performance-runs/other.json", True, LEAK_FOUND)],
)
def test_truncated_tukey_field_exception_is_limited_to_the_frozen_value(tmp_path: Path, relative: str, *, frozen: bool, status: int) -> None:
    """A clipped key is exempted only with its original value and report path."""
    original = b"[11.061047972370767,11.456482165636142,12.51097334767714,12.906407540942514]"
    other = b"[23.01652962024839,23.323939759650607,24.143700131389842,24.451110270792054]"
    encoded = base64.b64encode(original if frozen else other).decode("ascii")
    source = tmp_path / relative
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text(f'key.json": "{encoded}"\n', encoding="utf-8", newline="\n")
    result = scan_secrets(tmp_path)
    assert result.returncode == status, result.stderr


def scan_secrets(root: Path) -> subprocess.CompletedProcess[str]:
    """Exercise the installed scanner with the real consumer policy."""
    return run_command(
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
            "--exit-code",
            str(LEAK_FOUND),
            "--ignore-gitleaks-allow",
            str(root),
        ],
        cwd=REPO_ROOT,
        check=False,
    )
