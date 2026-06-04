"""Check that documentation dependency snippets match Cargo.toml."""

from __future__ import annotations

import re
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SKIP_DIRS = frozenset(
    {
        ".git",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        ".tmp_pycache",
        ".venv",
        "target",
    }
)


@dataclass(frozen=True)
class PackageInfo:
    """Cargo package identity used in documented dependency snippets."""

    name: str
    version: str


@dataclass(frozen=True)
class DependencySnippet:
    """A documented dependency version snippet for the current package."""

    path: Path
    line: int
    version: str
    text: str


@dataclass(frozen=True)
class VersionMismatch:
    """A dependency snippet whose version does not match Cargo.toml."""

    snippet: DependencySnippet
    package: PackageInfo


def _read_cargo_package_info(cargo_toml: Path) -> PackageInfo:
    data: dict[str, Any] = tomllib.loads(cargo_toml.read_text(encoding="utf-8"))
    package = data.get("package")
    if not isinstance(package, dict):
        msg = f"{cargo_toml} is missing a [package] table"
        raise TypeError(msg)

    name = package.get("name")
    if not isinstance(name, str):
        msg = f"{cargo_toml} is missing a string package.name"
        raise TypeError(msg)

    version = package.get("version")
    if not isinstance(version, str):
        msg = f"{cargo_toml} is missing a string package.version"
        raise TypeError(msg)
    return PackageInfo(name=name, version=version)


def _iter_markdown_files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*.md") if path.is_file() and not (set(path.relative_to(root).parts) & SKIP_DIRS))


def _dependency_regex(package_name: str) -> re.Pattern[str]:
    escaped_name = re.escape(package_name)
    return re.compile(rf'(?<![\w.-]){escaped_name}\s*=\s*(?:"(?P<plain>[^"]+)"|\{{\s*version\s*=\s*"(?P<table>[^"]+)")')


def _dependency_snippets(path: Path, package_name: str) -> list[DependencySnippet]:
    dependency_re = _dependency_regex(package_name)
    snippets: list[DependencySnippet] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        for match in dependency_re.finditer(line):
            version = match.group("plain") or match.group("table")
            snippets.append(
                DependencySnippet(
                    path=path,
                    line=line_number,
                    version=version,
                    text=line.strip(),
                )
            )
    return snippets


def find_version_mismatches(root: Path) -> list[VersionMismatch]:
    """Return documented dependency snippets for this crate that are stale."""

    package = _read_cargo_package_info(root / "Cargo.toml")
    mismatches: list[VersionMismatch] = []
    for path in _iter_markdown_files(root):
        for snippet in _dependency_snippets(path, package.name):
            if snippet.version != package.version:
                mismatches.append(VersionMismatch(snippet=snippet, package=package))
    return mismatches


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
    try:
        mismatches = find_version_mismatches(root)
    except (OSError, TypeError, tomllib.TOMLDecodeError) as error:
        print(f"Could not check documentation dependency versions: {error}", file=sys.stderr)
        return 1

    if not mismatches:
        return 0

    print("Documentation dependency snippets are out of sync with Cargo.toml:", file=sys.stderr)
    for mismatch in mismatches:
        snippet = mismatch.snippet
        rel_path = snippet.path.relative_to(root)
        print(
            f"  {rel_path}:{snippet.line}: {mismatch.package.name} found {snippet.version}, expected {mismatch.package.version}: {snippet.text}",
            file=sys.stderr,
        )
    return 1


if __name__ == "__main__":
    sys.exit(main())
