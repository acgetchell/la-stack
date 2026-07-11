"""Validate repository-owned Semgrep fixture annotations."""

from __future__ import annotations

import collections
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TypeGuard

RULE_ANNOTATION = re.compile(r"\bruleid:\s*([A-Za-z0-9_.-]+(?:\s*,\s*[A-Za-z0-9_.-]+)*)")


type ParsedObject = dict[str, object]


@dataclass(frozen=True, slots=True)
class SemgrepResults:
    """Validated subset of Semgrep JSON needed by fixture checks."""

    results: tuple[ParsedObject, ...]


def _is_parsed_object(value: object) -> TypeGuard[ParsedObject]:
    return isinstance(value, dict) and all(isinstance(key, str) for key in value)


def _path_argument(argv: list[str]) -> Path | None:
    if len(argv) <= 1:
        print("Missing required path argument: sys.argv[1]", file=sys.stderr)
        return None

    path = Path(argv[1])
    if not path.exists():
        print(f"path argument does not exist: {path}", file=sys.stderr)
        return None
    if not path.is_file():
        print(f"path argument is not a file: {path}", file=sys.stderr)
        return None
    if not os.access(path, os.R_OK):
        print(f"path argument is not readable: {path}", file=sys.stderr)
        return None
    return path


def _semgrep_results() -> SemgrepResults | None:
    semgrep_json = os.environ.get("SEMGREP_JSON")
    if semgrep_json is None:
        print("Missing required SEMGREP_JSON environment variable", file=sys.stderr)
        return None
    try:
        data: object = json.loads(semgrep_json)
    except json.JSONDecodeError as error:
        print(f"Invalid JSON in SEMGREP_JSON: {error}", file=sys.stderr)
        return None

    if not _is_parsed_object(data):
        print("Invalid SEMGREP_JSON shape: expected a JSON object", file=sys.stderr)
        return None
    results = data.get("results")
    if not isinstance(results, list):
        print("Invalid SEMGREP_JSON shape: expected 'results' to be a list", file=sys.stderr)
        return None

    parsed_results: list[ParsedObject] = []
    malformed_results: list[str] = []
    for index, result in enumerate(results):
        if _is_parsed_object(result):
            parsed_results.append(result)
        else:
            malformed_results.append(f"result {index} is not an object")

    if malformed_results:
        print("Invalid SEMGREP_JSON shape:", file=sys.stderr)
        for malformed in malformed_results:
            print(f"  {malformed}", file=sys.stderr)
        return None

    return SemgrepResults(results=tuple(parsed_results))


def _expected_rule_counts(path: Path) -> collections.Counter[str]:
    expected: collections.Counter[str] = collections.Counter()
    for line in path.read_text(encoding="utf-8").splitlines():
        for match in RULE_ANNOTATION.finditer(line):
            expected.update(rule_id.strip() for rule_id in match.group(1).split(",") if rule_id.strip())
    return expected


def _actual_rule_counts(semgrep: SemgrepResults) -> collections.Counter[str] | None:
    actual: collections.Counter[str] = collections.Counter()
    malformed_results: list[str] = []
    for index, result in enumerate(semgrep.results):
        check_id = result.get("check_id")
        if isinstance(check_id, str):
            actual.update([check_id])
        else:
            malformed_results.append(f"result {index} is missing string field 'check_id'")

    if not malformed_results:
        return actual

    print("Invalid SEMGREP_JSON shape:", file=sys.stderr)
    for malformed in malformed_results:
        print(f"  {malformed}", file=sys.stderr)
    return None


def main() -> int:
    """Compare expected fixture annotations with the supplied Semgrep results."""
    path = _path_argument(sys.argv)
    if path is None:
        return 1

    expected = _expected_rule_counts(path)

    semgrep = _semgrep_results()
    if semgrep is None:
        return 1

    actual = _actual_rule_counts(semgrep)
    if actual is None:
        return 1

    if actual == expected:
        return 0

    print(f"Semgrep fixture mismatch in {path}")
    for rule in sorted(expected.keys() | actual.keys()):
        if expected[rule] != actual[rule]:
            print(f"  {rule}: expected {expected[rule]}, got {actual[rule]}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
