"""Validate repository-owned Semgrep fixture annotations."""

from __future__ import annotations

import collections
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

RULE_ANNOTATION = re.compile(r"\bruleid:\s*([A-Za-z0-9_.-]+(?:\s*,\s*[A-Za-z0-9_.-]+)*)")


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


def _semgrep_results() -> dict[str, Any] | None:
    semgrep_json = os.environ.get("SEMGREP_JSON")
    if semgrep_json is None:
        print("Missing required SEMGREP_JSON environment variable", file=sys.stderr)
        return None
    try:
        data = json.loads(semgrep_json)
    except json.JSONDecodeError as error:
        print(f"Invalid JSON in SEMGREP_JSON: {error}", file=sys.stderr)
        return None

    if not isinstance(data, dict):
        print("Invalid SEMGREP_JSON shape: expected a JSON object", file=sys.stderr)
        return None
    if not isinstance(data.get("results"), list):
        print("Invalid SEMGREP_JSON shape: expected 'results' to be a list", file=sys.stderr)
        return None
    return data


def main() -> int:
    path = _path_argument(sys.argv)
    if path is None:
        return 1

    expected: collections.Counter[str] = collections.Counter()
    for line in path.read_text(encoding="utf-8").splitlines():
        for match in RULE_ANNOTATION.finditer(line):
            expected.update(rule_id.strip() for rule_id in match.group(1).split(",") if rule_id.strip())

    data = _semgrep_results()
    if data is None:
        return 1

    results = data["results"]
    actual: collections.Counter[str] = collections.Counter()
    malformed_results: list[str] = []
    for index, result in enumerate(results):
        if not isinstance(result, dict):
            malformed_results.append(f"result {index} is not an object")
            continue

        check_id = result.get("check_id")
        if not isinstance(check_id, str):
            malformed_results.append(f"result {index} is missing string field 'check_id'")
            continue

        actual.update([check_id])

    if malformed_results:
        print("Invalid SEMGREP_JSON shape:", file=sys.stderr)
        for malformed in malformed_results:
            print(f"  {malformed}", file=sys.stderr)
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
