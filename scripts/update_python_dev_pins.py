"""Advance exact Python development-tool pins with uv's resolver."""

import argparse
import re
import subprocess
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from subprocess_utils import ExecutableNotFoundError, run_safe_command

EXACT_REQUIREMENT = re.compile(
    r"^(?P<name>[A-Za-z0-9][A-Za-z0-9._-]*)==(?P<version>[^;\s]+)$",
)
RESOLVED_REQUIREMENT = re.compile(
    r"^(?P<name>[A-Za-z0-9][A-Za-z0-9._-]*)==(?P<version>[^;\s]+)(?:\s*;\s*.+)?$",
)
PYTHON_FLOOR = re.compile(r"^>=(?P<version>[0-9]+\.[0-9]+)$")


@dataclass(frozen=True)
class DevPin:
    """One exact development-tool requirement."""

    name: str
    version: str


def canonicalize_name(name: str) -> str:
    """Return the normalized distribution name used for comparisons."""
    return re.sub(r"[-_.]+", "-", name).casefold()


def _required_table(data: dict[str, object], name: str) -> dict[str, object]:
    """Return a required TOML table with a typed failure."""
    table = data.get(name)
    if not isinstance(table, dict):
        msg = f"pyproject.toml must contain a [{name}] table"
        raise TypeError(msg)
    return cast("dict[str, object]", table)


def _python_floor(project: dict[str, object]) -> str:
    """Return the repository's single supported Python lower bound."""
    requires_python = project.get("requires-python")
    if not isinstance(requires_python, str):
        msg = "project.requires-python must be a string"
        raise TypeError(msg)
    python_match = PYTHON_FLOOR.fullmatch(requires_python)
    if python_match is None:
        msg = f"expected project.requires-python to be a single lower bound, found: {requires_python}"
        raise ValueError(msg)
    return cast("str", python_match.group("version"))


def _dev_pins(groups: dict[str, object]) -> list[DevPin]:
    """Return validated exact pins from dependency-groups.dev."""
    dev = groups.get("dev")
    if not isinstance(dev, list):
        msg = "dependency-groups.dev must be an array"
        raise TypeError(msg)
    if not dev:
        msg = "dependency-groups.dev must be a non-empty array"
        raise ValueError(msg)

    pins: list[DevPin] = []
    normalized_names: set[str] = set()
    for requirement in dev:
        if not isinstance(requirement, str):
            msg = "dependency-groups.dev entries must be strings"
            raise TypeError(msg)
        requirement_match = EXACT_REQUIREMENT.fullmatch(requirement)
        if requirement_match is None:
            msg = f"development-tool requirements must be exact simple pins: {requirement}"
            raise ValueError(msg)
        pin = DevPin(requirement_match.group("name"), requirement_match.group("version"))
        normalized = canonicalize_name(pin.name)
        if normalized in normalized_names:
            msg = f"duplicate development-tool requirement: {pin.name}"
            raise ValueError(msg)
        normalized_names.add(normalized)
        pins.append(pin)
    return pins


def parse_project(text: str) -> tuple[str, list[DevPin]]:
    """Parse the Python floor and exact development-tool pins."""
    data = tomllib.loads(text)
    return _python_floor(_required_table(data, "project")), _dev_pins(_required_table(data, "dependency-groups"))


def parse_resolution(output: str, pins: list[DevPin]) -> list[DevPin]:
    """Extract one resolver-selected version for every direct development tool."""
    requested = {canonicalize_name(pin.name): pin.name for pin in pins}
    resolved: dict[str, set[str]] = {name: set() for name in requested}

    for raw_line in output.splitlines():
        match = RESOLVED_REQUIREMENT.fullmatch(raw_line.strip())
        if match is None:
            continue
        normalized = canonicalize_name(match.group("name"))
        if normalized in resolved:
            resolved[normalized].add(match.group("version"))

    latest: list[DevPin] = []
    for pin in pins:
        versions = resolved[canonicalize_name(pin.name)]
        if not versions:
            msg = f"uv resolver output omitted direct development tool: {pin.name}"
            raise ValueError(msg)
        if len(versions) != 1:
            rendered = ", ".join(sorted(versions))
            msg = f"uv resolver selected multiple versions for {pin.name}: {rendered}"
            raise ValueError(msg)
        latest.append(DevPin(pin.name, next(iter(versions))))
    return latest


def resolve_latest_pins(pins: list[DevPin], python_version: str, project_root: Path) -> list[DevPin]:
    """Resolve the latest cross-platform set without changing repository files."""
    requirements = "".join(f"{pin.name}\n" for pin in pins)
    result = run_safe_command(
        "uv",
        [
            "pip",
            "compile",
            "-",
            "--universal",
            "--no-header",
            "--no-annotate",
            "--python-version",
            python_version,
        ],
        cwd=project_root,
        input=requirements,
    )
    return parse_resolution(result.stdout, pins)


def update_dev_pins(pyproject: Path) -> dict[str, tuple[str, str]]:
    """Resolve and apply all changed direct pins in one uv transaction."""
    python_version, current = parse_project(pyproject.read_text(encoding="utf-8"))
    latest = resolve_latest_pins(current, python_version, pyproject.parent)
    changes = {old.name: (old.version, new.version) for old, new in zip(current, latest, strict=True) if old.version != new.version}
    if not changes:
        return changes

    run_safe_command(
        "uv",
        ["add", "--dev", "--no-sync", *(f"{pin.name}=={pin.version}" for pin in latest)],
        cwd=pyproject.parent,
    )
    return changes


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pyproject",
        type=Path,
        default=Path("pyproject.toml"),
        help="project manifest containing exact dependency-groups.dev pins",
    )
    return parser.parse_args(argv)


def _subprocess_detail(error: subprocess.CalledProcessError) -> str:
    """Return captured tool diagnostics when available."""
    detail = error.stderr if isinstance(error.stderr, str) else error.stdout if isinstance(error.stdout, str) else ""
    detail = detail.strip()
    return cast("str", detail or str(error))


def main(argv: list[str] | None = None) -> int:
    """Advance exact development-tool pins and their uv lock resolution."""
    args = parse_args(argv)
    try:
        changes = update_dev_pins(args.pyproject)
    except subprocess.CalledProcessError as error:
        print(f"failed to update Python development-tool pins: {_subprocess_detail(error)}", file=sys.stderr)
        return 1
    except (ExecutableNotFoundError, OSError, subprocess.TimeoutExpired, TypeError, ValueError) as error:
        print(f"failed to update Python development-tool pins: {error}", file=sys.stderr)
        return 1

    if not changes:
        print("Python development-tool pins are already current.")
        return 0
    for name, (old_version, new_version) in changes.items():
        print(f"Updated {name}: {old_version} -> {new_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
