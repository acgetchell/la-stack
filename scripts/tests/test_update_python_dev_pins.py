"""Tests for resolver-backed Python development-tool pin updates."""

import subprocess
import tomllib
from pathlib import Path

import pytest

import update_python_dev_pins


def project_text(*requirements: str) -> str:
    """Return a minimal project with exact development-tool pins."""
    rendered = "\n".join(f'    "{requirement}",' for requirement in requirements)
    return f"""[project]
name = "fixture"
version = "0.1.0"
requires-python = ">=3.14"

[dependency-groups]
dev = [
{rendered}
]
"""


def test_parse_project_accepts_exact_simple_dev_pins() -> None:
    python_version, pins = update_python_dev_pins.parse_project(project_text("ruff==0.16.2", "semgrep==1.172.0"))

    assert python_version == "3.14"
    assert pins == [
        update_python_dev_pins.DevPin("ruff", "0.16.2"),
        update_python_dev_pins.DevPin("semgrep", "1.172.0"),
    ]


def test_actual_project_keeps_included_tooling_outside_direct_dev_updates() -> None:
    text = (Path(__file__).resolve().parents[2] / "pyproject.toml").read_text(encoding="utf-8")
    python_version, pins = update_python_dev_pins.parse_project(text)
    groups = tomllib.loads(text)["dependency-groups"]
    assert python_version == "3.14"
    assert groups["tooling"] == ["research-repo-tools==0.1.7"]
    assert {"include-group": "tooling"} in groups["dev"]
    assert "research-repo-tools" not in {pin.name for pin in pins}
    assert [f"{pin.name}=={pin.version}" for pin in pins] == [entry for entry in groups["dev"] if isinstance(entry, str)]


@pytest.mark.parametrize(
    ("entry", "error", "diagnostic"),
    [
        ('{ include-group = "missing" }', ValueError, "existing dependency group"),
        ("{ include-group = 42 }", ValueError, "existing dependency group"),
        ('{ include-group = "dev" }', ValueError, "other than dev"),
        ('{ include-group = "tooling", extra = true }', TypeError, "strings or include-group tables"),
        ("42", TypeError, "strings or include-group tables"),
    ],
)
def test_parse_project_rejects_malformed_dev_entries(entry: str, error: type[Exception], diagnostic: str) -> None:
    text = project_text("ruff==0.16.2").replace("dev = [", f"dev = [\n    {entry},") + 'tooling = ["research-repo-tools==0.1.7"]\n'
    with pytest.raises(error, match=diagnostic):
        update_python_dev_pins.parse_project(text)


def test_parse_project_rejects_non_exact_dev_requirement() -> None:
    with pytest.raises(ValueError, match=r"development-tool requirements must be exact simple pins: ruff>=0\.16"):
        update_python_dev_pins.parse_project(project_text("ruff>=0.16"))


def test_parse_resolution_preserves_direct_order_and_ignores_transitives() -> None:
    pins = [
        update_python_dev_pins.DevPin("ruff", "0.16.2"),
        update_python_dev_pins.DevPin("semgrep", "1.172.0"),
    ]
    output = "packaging==26.3\nsemgrep==1.174.0\nruff==0.16.4\nmcp==1.29.0\n"

    assert update_python_dev_pins.parse_resolution(output, pins) == [
        update_python_dev_pins.DevPin("ruff", "0.16.4"),
        update_python_dev_pins.DevPin("semgrep", "1.174.0"),
    ]


def test_parse_resolution_rejects_missing_direct_tool() -> None:
    pins = [update_python_dev_pins.DevPin("semgrep", "1.172.0")]

    with pytest.raises(ValueError, match="uv resolver output omitted direct development tool: semgrep"):
        update_python_dev_pins.parse_resolution("mcp==1.29.0\n", pins)


@pytest.mark.parametrize("include_tooling", [False, True])
def test_update_dev_pins_resolves_then_applies_one_exact_transaction(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    include_tooling: bool,
) -> None:
    pyproject = tmp_path / "pyproject.toml"
    text = project_text("ruff==0.16.2", "semgrep==1.172.0")
    if include_tooling:
        text = text.replace("dev = [", 'dev = [\n    { include-group = "tooling" },') + 'tooling = ["research-repo-tools==0.1.7"]\n'
    pyproject.write_text(text, encoding="utf-8")
    calls: list[tuple[str, list[str], dict[str, object]]] = []

    def fake_run(command: str, args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append((command, args, kwargs))
        output = "ruff==0.16.4\nsemgrep==1.174.0\nmcp==1.29.0\n" if args[:2] == ["pip", "compile"] else ""
        return subprocess.CompletedProcess([command, *args], 0, stdout=output, stderr="")

    monkeypatch.setattr(update_python_dev_pins, "run_safe_command", fake_run)

    changes = update_python_dev_pins.update_dev_pins(pyproject)

    assert changes == {
        "ruff": ("0.16.2", "0.16.4"),
        "semgrep": ("1.172.0", "1.174.0"),
    }
    assert calls[0] == (
        "uv",
        [
            "pip",
            "compile",
            "-",
            "--universal",
            "--no-header",
            "--no-annotate",
            "--python-version",
            "3.14",
        ],
        {"cwd": tmp_path, "input": "ruff\nsemgrep\n"},
    )
    assert calls[1] == (
        "uv",
        ["add", "--dev", "--no-sync", "ruff==0.16.4", "semgrep==1.174.0"],
        {"cwd": tmp_path},
    )


def test_main_reports_uv_diagnostics_without_traceback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(project_text("semgrep==1.172.0"), encoding="utf-8")

    def failed_uv(_command: str, args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        raise subprocess.CalledProcessError(1, ["uv", *args], stderr="resolver conflict")

    monkeypatch.setattr(update_python_dev_pins, "run_safe_command", failed_uv)

    assert update_python_dev_pins.main(["--pyproject", str(pyproject)]) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "failed to update Python development-tool pins: resolver conflict\n"
