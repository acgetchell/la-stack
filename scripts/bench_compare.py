#!/usr/bin/env python3
"""Compare exact-arithmetic benchmark results across Criterion baselines.

Reads Criterion output under:
  target/criterion/{group}/{bench}/{sample}/estimates.json

And writes a markdown performance table to docs/PERFORMANCE.md.

Typical workflow (see docs/RELEASING.md):

  # Save baseline at a release tag
  just bench-save-baseline v0.3.0

  # ... make optimisations ...

  # Compare current performance against the saved baseline
  just bench-compare v0.3.0

  # Or just generate a snapshot of current performance
  just bench-compare
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Benchmark group / bench discovery
# ---------------------------------------------------------------------------

# Groups and the benchmarks within each group that we track.
EXACT_GROUPS: dict[str, list[str]] = {
    "exact_d2": ["det", "det_direct", "det_exact", "det_exact_f64", "det_sign_exact", "solve_exact", "solve_exact_f64"],
    "exact_d3": ["det", "det_direct", "det_exact", "det_exact_f64", "det_sign_exact", "solve_exact", "solve_exact_f64"],
    "exact_d4": ["det", "det_direct", "det_exact", "det_exact_f64", "det_sign_exact", "solve_exact", "solve_exact_f64"],
    "exact_d5": ["det", "det_direct", "det_exact", "det_exact_f64", "det_sign_exact", "solve_exact", "solve_exact_f64"],
    "exact_near_singular_3x3": ["det_sign_exact", "det_exact"],
}


@dataclass(frozen=True, slots=True)
class BenchResult:
    """A single benchmark measurement (point estimate + confidence interval)."""

    group: str
    bench: str
    point_ns: float
    ci_lo_ns: float
    ci_hi_ns: float


@dataclass(frozen=True, slots=True)
class Comparison:
    """A comparison between baseline and current benchmark results."""

    group: str
    bench: str
    baseline_ns: float
    current_ns: float
    speedup: float  # baseline / current (>1 = faster)
    pct_change: float  # signed percent change (negative = faster)


# ---------------------------------------------------------------------------
# Criterion JSON parsing
# ---------------------------------------------------------------------------


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _read_estimate(estimates_json: Path, stat: str = "median") -> tuple[float, float, float]:
    """Read a point estimate and confidence interval from Criterion estimates.json.

    Returns (point_ns, ci_lo_ns, ci_hi_ns).
    """
    data = json.loads(estimates_json.read_text(encoding="utf-8"))

    stat_obj = data.get(stat)
    if not isinstance(stat_obj, dict):
        msg = f"stat '{stat}' not found in {estimates_json}"
        raise KeyError(msg)

    point = float(stat_obj["point_estimate"])
    ci = stat_obj.get("confidence_interval")
    if not isinstance(ci, dict):
        return (point, point, point)

    lo = float(ci.get("lower_bound", point))
    hi = float(ci.get("upper_bound", point))
    return (point, lo, hi)


def _collect_results(criterion_dir: Path, sample: str, stat: str) -> list[BenchResult]:
    """Collect all benchmark results from Criterion output."""
    results: list[BenchResult] = []

    for group, benches in EXACT_GROUPS.items():
        group_dir = criterion_dir / group
        if not group_dir.is_dir():
            continue

        for bench in benches:
            est_path = group_dir / bench / sample / "estimates.json"
            if not est_path.exists():
                continue

            point, lo, hi = _read_estimate(est_path, stat)
            results.append(BenchResult(group=group, bench=bench, point_ns=point, ci_lo_ns=lo, ci_hi_ns=hi))

    return results


def _collect_comparisons(
    criterion_dir: Path,
    baseline_name: str,
    stat: str,
) -> list[Comparison]:
    """Compare current (new) results against a named baseline."""
    comparisons: list[Comparison] = []

    for group, benches in EXACT_GROUPS.items():
        group_dir = criterion_dir / group
        if not group_dir.is_dir():
            continue

        for bench in benches:
            new_path = group_dir / bench / "new" / "estimates.json"
            base_path = group_dir / bench / baseline_name / "estimates.json"

            if not new_path.exists() or not base_path.exists():
                continue

            new_point, _, _ = _read_estimate(new_path, stat)
            base_point, _, _ = _read_estimate(base_path, stat)

            speedup = base_point / new_point if new_point > 0 else float("inf")
            pct_change = ((new_point - base_point) / base_point) * 100.0 if base_point > 0 else 0.0

            comparisons.append(
                Comparison(
                    group=group,
                    bench=bench,
                    baseline_ns=base_point,
                    current_ns=new_point,
                    speedup=speedup,
                    pct_change=pct_change,
                )
            )

    return comparisons


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def _format_time(ns: float) -> str:
    """Format nanoseconds into a human-readable string."""
    if ns < 1_000:
        return f"{ns:.1f} ns"
    if ns < 1_000_000:
        return f"{ns / 1_000:.2f} µs"
    return f"{ns / 1_000_000:.2f} ms"


def _format_pct(pct: float) -> str:
    """Format percent change with sign and colour hint."""
    if pct < -1.0:
        return f"**{pct:+.1f}%**"  # bold for improvements
    if pct > 1.0:
        return f"{pct:+.1f}%"
    return f"{pct:+.1f}%"


def _snapshot_table(results: list[BenchResult]) -> str:
    """Generate a markdown table for a single set of results (no baseline)."""
    lines = [
        "| Group | Benchmark | Median | 95% CI |",
        "|-------|-----------|-------:|-------:|",
    ]

    for r in results:
        ci_range = f"[{_format_time(r.ci_lo_ns)}, {_format_time(r.ci_hi_ns)}]"
        lines.append(f"| {r.group} | {r.bench} | {_format_time(r.point_ns)} | {ci_range} |")

    return "\n".join(lines)


def _comparison_table(comparisons: list[Comparison], baseline_name: str) -> str:
    """Generate a markdown table comparing baseline vs current."""
    lines = [
        f"| Group | Benchmark | {baseline_name} | Current | Change | Speedup |",
        "|-------|-----------|-------:|--------:|-------:|--------:|",
    ]

    for c in comparisons:
        lines.append(
            f"| {c.group} | {c.bench} "
            f"| {_format_time(c.baseline_ns)} "
            f"| {_format_time(c.current_ns)} "
            f"| {_format_pct(c.pct_change)} "
            f"| {c.speedup:.2f}x |"
        )

    return "\n".join(lines)


def _read_cargo_version(root: Path) -> str:
    cargo_toml = root / "Cargo.toml"
    if not cargo_toml.exists():
        return "unknown"
    for line in cargo_toml.read_text(encoding="utf-8").splitlines():
        m = re.match(r'^\s*version\s*=\s*"([^"]+)"', line)
        if m:
            return m.group(1)
    return "unknown"


def _get_git_info(root: Path) -> tuple[str, str]:
    """Return (short_hash, branch_or_tag)."""
    git_path = shutil.which("git")
    if git_path is None:
        return ("unknown", "unknown")

    short_hash = "unknown"
    branch = "unknown"
    try:
        result = subprocess.run(  # noqa: S603
            [git_path, "--no-pager", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=root,
        )
        short_hash = result.stdout.strip()
    except subprocess.CalledProcessError:
        pass
    try:
        result = subprocess.run(  # noqa: S603
            [git_path, "--no-pager", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=root,
        )
        branch = result.stdout.strip()
    except subprocess.CalledProcessError:
        pass
    return short_hash, branch


def _generate_markdown(
    root: Path,
    table: str,
    baseline_name: str | None,
    stat: str,
) -> str:
    """Generate the complete PERFORMANCE.md content."""
    version = _read_cargo_version(root)
    short_hash, branch = _get_git_info(root)
    now = datetime.now(tz=UTC).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = [
        "# Exact Arithmetic Performance",
        "",
        f"**la-stack** v{version} · `{short_hash}` ({branch}) · {now}",
        f"**Statistic**: {stat}",
        "",
        "## Benchmark Results",
        "",
    ]

    if baseline_name:
        lines.append(f"Comparison against baseline **{baseline_name}**:")
        lines.append("")
        lines.append("Negative change = faster. Speedup > 1.00x = improvement.")
    else:
        lines.append("Current performance snapshot (no baseline comparison).")

    lines.extend(["", table, ""])

    lines.extend(
        [
            "## How to Update",
            "",
            "```bash",
            "# Save a baseline at the current release",
            "just bench-save-baseline v0.3.0",
            "",
            "# Compare current code against a saved baseline",
            "just bench-compare v0.3.0",
            "",
            "# Generate a snapshot without comparison",
            "just bench-compare",
            "```",
            "",
            "See `docs/RELEASING.md` for the full release workflow.",
        ]
    )

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare exact-arithmetic benchmark results across Criterion baselines.",
    )
    parser.add_argument(
        "baseline",
        nargs="?",
        default=None,
        help="Baseline name to compare against (e.g. 'v0.3.0'). Omit for a snapshot.",
    )
    parser.add_argument(
        "--stat",
        default="median",
        choices=["mean", "median"],
        help="Statistic to compare (default: median).",
    )
    parser.add_argument(
        "--criterion-dir",
        default="target/criterion",
        help="Criterion output directory (default: target/criterion).",
    )
    parser.add_argument(
        "--output",
        default="docs/PERFORMANCE.md",
        help="Output markdown file (default: docs/PERFORMANCE.md).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    root = _repo_root()
    criterion_dir = root / args.criterion_dir
    output_path = Path(args.output) if Path(args.output).is_absolute() else root / args.output

    if not criterion_dir.is_dir():
        print(
            f"No Criterion results found at {criterion_dir}.\nRun benchmarks first:\n  just bench-exact\n",
            file=sys.stderr,
        )
        return 2

    if args.baseline:
        comparisons = _collect_comparisons(criterion_dir, args.baseline, args.stat)
        if not comparisons:
            print(
                f"No comparison data found for baseline '{args.baseline}'.\n"
                f"Save a baseline first:\n  just bench-save-baseline {args.baseline}\n"
                "Then run benchmarks:\n  just bench-exact\n",
                file=sys.stderr,
            )
            return 2
        table = _comparison_table(comparisons, args.baseline)
    else:
        results = _collect_results(criterion_dir, "new", args.stat)
        if not results:
            print(
                "No benchmark results found.\nRun benchmarks first:\n  just bench-exact\n",
                file=sys.stderr,
            )
            return 2
        table = _snapshot_table(results)

    md = _generate_markdown(root, table, args.baseline, args.stat)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(md, encoding="utf-8")
    print(f"📊 Wrote {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
