"""Inventory and validate complete release Criterion datasets before packaging."""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from bench_compare import (
    EXACT_GROUPS,
    VS_LINALG_CANONICAL_DIMS,
    VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM,
    VS_LINALG_STANDARD_BENCH_ORDER,
)
from criterion_measurements import read_object, validate_measurement
from subprocess_utils import ExecutableNotFoundError, format_exception_diagnostics, run_cargo_command

SUITES = {"vs_linalg": "bench", "exact": "bench,exact"}
RAW_FILES = ("benchmark.json", "estimates.json", "sample.json", "tukey.json")


def parse_benchmark_list(output: str) -> list[str]:
    """Read Criterion's full IDs, rejecting empty or duplicate inventories."""
    ids = [line.removesuffix(": benchmark") for line in output.splitlines() if line.endswith(": benchmark")]
    if not ids or len(ids) != len(set(ids)) or any("/" not in name for name in ids):
        msg = "Criterion inventory must contain unique, nonempty benchmark IDs"
        raise ValueError(msg)
    return sorted(ids)


def required_report_ids() -> set[str]:
    """Include every release-report row and all canonical README peer rows."""
    return {
        *(f"{group}/{bench}" for group, benches in EXACT_GROUPS.items() for bench in benches),
        *(
            f"d{dimension}/{bench}"
            for dimension in VS_LINALG_CANONICAL_DIMS
            for bench in (*VS_LINALG_STANDARD_BENCH_ORDER, *VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM.get(dimension, []))
        ),
    }


def inventory_ids(data: object) -> set[str]:
    """Validate the manifest and its coverage of the report consumers."""
    if not isinstance(data, dict) or set(data) != set(SUITES):
        msg = "inventory must contain exactly vs_linalg and exact suites"
        raise ValueError(msg)
    expected: set[str] = set()
    for suite, ids in data.items():
        if not isinstance(ids, list) or not ids or any(not isinstance(name, str) or not name for name in ids):
            raise ValueError(f"invalid benchmark inventory for {suite}")
        if len(ids) != len(set(ids)) or expected.intersection(ids):
            raise ValueError(f"duplicate benchmark IDs in {suite}")
        expected.update(ids)
    missing = required_report_ids() - expected
    if missing:
        raise ValueError(f"inventory omits report consumers: {', '.join(sorted(missing))}")
    return expected


def discover(root: Path, manifest: Path, criterion: Path) -> None:
    """Compile and list both unfiltered suites without collecting measurements."""
    if criterion.exists() and any(criterion.iterdir()):
        raise ValueError(f"release inventory requires a fresh Criterion directory: {criterion}")
    inventory: dict[str, list[str]] = {}
    for suite, features in SUITES.items():
        print(f"[release-baseline] Discovering {suite}", flush=True)
        package = ["-p", "la-stack-comparison"] if suite == "vs_linalg" else []
        result = run_cargo_command(
            ["bench", "--locked", *package, "--features", features, "--bench", suite, "--", "--list"],
            cwd=root,
            capture_output=False,
            stdout=subprocess.PIPE,
            timeout=None,  # The workflow bounds compilation and discovery together.
        )
        inventory[suite] = parse_benchmark_list(result.stdout)
        count = len(inventory[suite])
        print(
            f"[release-baseline] {suite}: {count} benchmarks; "
            f"nominal warmup + measurement {count * 8 / 60:.1f} min; "
            f"planning estimate at 12 s/benchmark {count * 12 / 60:.1f} min (not an upper bound)",
            flush=True,
        )
    inventory_ids(inventory)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8", newline="\n")


def validate(criterion: Path, manifest: Path, baseline: str) -> int:
    """Fail closed on missing, corrupt, stale, or unexpected benchmark output."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", baseline) or baseline in {"new", "base", "change", "report"}:
        msg = "baseline must be a safe, nonreserved directory name"
        raise ValueError(msg)
    expected = inventory_ids(json.loads(manifest.read_text(encoding="utf-8")))
    report_ids = required_report_ids()
    observed: set[str] = set()
    for metadata in sorted(criterion.glob("**/new/benchmark.json")):
        benchmark = read_object(metadata).get("full_id")
        if not isinstance(benchmark, str) or benchmark not in expected or benchmark in observed:
            raise ValueError(f"unexpected or duplicate benchmark ID in {metadata}: {benchmark!r}")
        directory = metadata.parent
        if benchmark in report_ids and directory.parent != criterion / benchmark:
            raise ValueError(f"report benchmark is outside its consumer path: {metadata}")
        saved = directory.parent / baseline
        for filename in RAW_FILES:
            # Criterion can log a write/copy error and still exit successfully.
            # Require every saved raw file to match the freshly measured sample.
            source = directory / filename
            destination = saved / filename
            if source.read_bytes() != destination.read_bytes():
                raise ValueError(f"saved baseline differs from new measurement: {destination}")
            json.loads(source.read_text(encoding="utf-8"))
        validate_measurement(directory)
        observed.add(benchmark)
    missing = expected - observed
    if missing:
        raise ValueError(f"incomplete Criterion dataset: {', '.join(sorted(missing))}")
    return len(observed)


def main(argv: list[str] | None = None) -> int:
    """Provide the release workflow's inventory and publication gates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("inventory", "validate"))
    parser.add_argument("--baseline")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=Path("target/release-benchmark-inventory.json"))
    parser.add_argument("--criterion-dir", type=Path, default=Path(os.environ.get("CRITERION_HOME", "target/criterion")))
    args = parser.parse_args(argv)
    try:
        if args.command == "inventory":
            discover(args.root, args.manifest, args.criterion_dir)
        else:
            if not args.baseline:
                parser.error("validate requires --baseline")
            count = validate(args.criterion_dir, args.manifest, args.baseline)
            print(f"[release-baseline] Validated {count} complete benchmarks in new and {args.baseline}")
        return 0
    except (OSError, ValueError, TypeError, ExecutableNotFoundError, subprocess.SubprocessError) as error:
        print(f"[release-baseline] {format_exception_diagnostics(error)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
