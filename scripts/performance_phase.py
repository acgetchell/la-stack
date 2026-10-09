"""Consumer Cargo commands for the captured common benchmark harness."""

import argparse
import os
from pathlib import Path

from research_repo_tools.process import run_command, run_command_live

SUITES = {"vs_linalg": ("-p", "la-stack-comparison", "--features", "bench"), "exact": ("--features", "bench,exact")}


def commands(suite: str, phase: str) -> tuple[tuple[str, ...], ...]:
    """Declare complete baseline timing and current la-stack-only timing."""
    selected = SUITES if suite == "all" else {suite: SUITES[suite]}
    return tuple(
        ("cargo", "bench", "--locked", *features, "--bench", name, "--", *(("la_stack",) if phase == "current" and name == "vs_linalg" else ()), "--noplot")
        for name, features in selected.items()
    )


def discover(root: Path, suite: str, phase: str, env: dict[str, str]) -> tuple[str, ...]:
    """Ask native Criterion binaries for semantic IDs, without timing samples."""
    expected: list[str] = []
    for command in commands(suite, phase):
        result = run_command(command[0], [*command[1:], "--list"], cwd=root, env=env, timeout=7200)
        ids = [line.removesuffix(": benchmark") for line in result.stdout.splitlines() if line.endswith(": benchmark")]
        if not ids:
            raise ValueError(f"empty Criterion inventory from {command}")
        expected.extend(ids)
    if len(expected) != len(set(expected)):
        msg = "duplicate semantic benchmark IDs"
        raise ValueError(msg)
    return tuple(sorted(expected))


def main() -> None:
    """Execute one declared correctness or timing phase with inherited live output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("gate", "measure"))
    parser.add_argument("suite", choices=("all", *SUITES))
    parser.add_argument("phase", choices=("baseline", "current"))
    args = parser.parse_args()
    selected = commands(args.suite, args.phase)
    if args.action == "gate":
        selected = (("cargo", "test", "--locked", "--workspace", "--features", "bench,exact", "--test", "vs_linalg_inputs", "--test", "exact_bench_config"),)
    for command in selected:
        print(f"[performance] {args.phase} {args.action}: {' '.join(command)}", flush=True)
        run_command_live(command[0], command[1:], cwd=Path.cwd(), env=os.environ, timeout=7200)


if __name__ == "__main__":
    main()
