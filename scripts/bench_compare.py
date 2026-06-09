#!/usr/bin/env python3
"""Compare la-stack benchmark results across Criterion baselines.

Reads Criterion output under:
  target/criterion/{group}/{bench}/{sample}/estimates.json

And writes a local markdown performance report.

Typical workflow (see docs/RELEASING.md):

  # Save a full baseline for the last release
  just bench-save-last

  # ... make optimisations ...

  # Run the cheap latest-vs-last workflow
  just bench-latest-vs-last

  # Or just render the report from existing Criterion output
  just bench-compare
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tomllib
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Protocol

from criterion_dim_plot import METRICS
from subprocess_utils import ExecutableNotFoundError, run_git_command

# ---------------------------------------------------------------------------
# Benchmark group / bench discovery
# ---------------------------------------------------------------------------

# Groups and the benchmarks within each group that we track.
#
# Mirrors the structure of `benches/exact.rs`: general-case per-dimension
# groups (`exact_d{2..5}`), fixed-seed random percentile groups, plus
# adversarial/extreme-input groups that share a fixed five-bench layout
# (`det_sign_exact`, `det_exact`, `solve_exact`, `solve_exact_f64_result`,
# `solve_exact_rounded_f64`).
_EXTREME_BENCHES: list[str] = [
    "det_sign_exact",
    "det_exact",
    "solve_exact",
    "solve_exact_f64_result",
    "solve_exact_rounded_f64",
]
_RANDOM_PERCENTILE_BENCHES: list[str] = [f"{operation}_{percentile}" for operation in _EXTREME_BENCHES for percentile in ("p50", "p95", "p99")]

_EXACT_DIMENSION_BENCHES: list[str] = [
    "det",
    "det_direct",
    "det_exact",
    "det_exact_f64_result",
    "det_exact_rounded_f64",
    "det_sign_exact",
    "solve_exact",
    "solve_exact_f64_result",
    "solve_exact_rounded_f64",
]

EXACT_GROUPS: dict[str, list[str]] = {
    "exact_d2": _EXACT_DIMENSION_BENCHES,
    "exact_d3": _EXACT_DIMENSION_BENCHES,
    "exact_d4": _EXACT_DIMENSION_BENCHES,
    "exact_d5": _EXACT_DIMENSION_BENCHES,
    "exact_random_percentile_d2": _RANDOM_PERCENTILE_BENCHES,
    "exact_random_percentile_d3": _RANDOM_PERCENTILE_BENCHES,
    "exact_random_percentile_d4": _RANDOM_PERCENTILE_BENCHES,
    "exact_random_percentile_d5": _RANDOM_PERCENTILE_BENCHES,
    "exact_near_singular_3x3": _EXTREME_BENCHES,
    "exact_large_entries_3x3": _EXTREME_BENCHES,
    "exact_hilbert_4x4": _EXTREME_BENCHES,
    "exact_hilbert_5x5": _EXTREME_BENCHES,
}

EXACT_RELEASE_SIGNAL_GROUPS: frozenset[str] = frozenset(group for group in EXACT_GROUPS if not group.startswith("exact_random_percentile_d"))

# v0.4.2 and earlier named the lossy exact-to-f64 benches after the public
# `*_exact_f64` API. Current benches split that behavior into strict `*_result`
# and lossy `*_rounded_f64` variants. Use the old baseline when present so
# release reports show both compatibility-successor performance and strict
# conversion overhead instead of silently dropping the new rows.
EXACT_LEGACY_BASELINE_BENCHES: dict[str, str] = {
    "det_exact_f64_result": "det_exact_f64",
    "det_exact_rounded_f64": "det_exact_f64",
    "solve_exact_f64_result": "solve_exact_f64",
    "solve_exact_rounded_f64": "solve_exact_f64",
}

_EXACT_LEGACY_PREFIX_BASELINE_BENCHES: tuple[tuple[str, str], ...] = (
    ("solve_exact_f64_result_", "solve_exact_f64_"),
    ("solve_exact_rounded_f64_", "solve_exact_f64_"),
)

VS_LINALG_LA_STACK_ONLY_BENCHES_BY_METRIC: dict[str, list[str]] = {
    "det_via_lu": ["la_stack_det"],
}
VS_LINALG_EXTRA_PEER_BENCHES_BY_METRIC: dict[str, list[tuple[str, str, str]]] = {
    "lu": [("la_stack_ldlt", "nalgebra_cholesky", "faer_ldlt")],
    "lu_solve": [("la_stack_ldlt_solve", "nalgebra_cholesky_solve", "faer_ldlt_solve")],
    "solve_from_lu": [("la_stack_solve_from_ldlt", "nalgebra_solve_from_cholesky", "faer_solve_from_ldlt")],
    "det_from_lu": [("la_stack_det_from_ldlt", "nalgebra_det_from_cholesky", "faer_det_from_ldlt")],
}
VS_LINALG_BENCH_ORDER: list[str] = [
    bench
    for metric_key, metric in METRICS.items()
    for bench in (
        [
            metric.la_bench,
            *VS_LINALG_LA_STACK_ONLY_BENCHES_BY_METRIC.get(metric_key, []),
            metric.na_bench,
            metric.fa_bench,
            *[peer_bench for peer_group in VS_LINALG_EXTRA_PEER_BENCHES_BY_METRIC.get(metric_key, []) for peer_bench in peer_group],
        ]
    )
]

VS_LINALG_LA_STACK_BENCHES: frozenset[str] = frozenset(bench for bench in VS_LINALG_BENCH_ORDER if bench.startswith("la_stack_"))
VS_LINALG_BASELINE_PEERS: dict[str, tuple[str, str]] = {metric.la_bench: (metric.na_bench, metric.fa_bench) for metric in METRICS.values()}
VS_LINALG_BASELINE_PEERS.update(
    {
        la_stack_bench: (nalgebra_bench, faer_bench)
        for peer_groups in VS_LINALG_EXTRA_PEER_BENCHES_BY_METRIC.values()
        for la_stack_bench, nalgebra_bench, faer_bench in peer_groups
    }
)

SUITE_CHOICES: tuple[str, ...] = ("all", "exact", "vs_linalg")
SCOPE_CHOICES: tuple[str, ...] = ("release-signal", "all-benches")


@dataclass(frozen=True, slots=True)
class BenchResult:
    """A single benchmark measurement (point estimate + confidence interval)."""

    suite: str
    group: str
    bench: str
    point_ns: float
    ci_lo_ns: float
    ci_hi_ns: float


@dataclass(frozen=True, slots=True)
class Comparison:
    """A comparison between baseline and current benchmark results."""

    suite: str
    group: str
    bench: str
    baseline_ns: float
    current_ns: float
    speedup: float  # baseline / current (>1 = faster)
    pct_change: float  # signed percent change (negative = faster)
    baseline_bench: str | None = None
    baseline_nalgebra_ns: float | None = None
    baseline_faer_ns: float | None = None


@dataclass(frozen=True, slots=True)
class ReportSettings:
    """Settings rendered into the benchmark report header."""

    baseline_name: str | None
    stat: str
    suite: str
    scope: str


# ---------------------------------------------------------------------------
# Criterion JSON parsing
# ---------------------------------------------------------------------------


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _dim_from_vs_linalg_group(name: str) -> int | None:
    """Parse a vs_linalg Criterion group name such as d2 or d64."""
    if not name.startswith("d"):
        return None
    suffix = name.removeprefix("d")
    if not suffix.isdecimal():
        return None
    return int(suffix)


def _read_estimate(estimates_json: Path, stat: str = "median") -> tuple[float, float, float]:
    """Read a point estimate and confidence interval from Criterion estimates.json.

    Returns (point_ns, ci_lo_ns, ci_hi_ns).
    """
    try:
        data = json.loads(estimates_json.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        msg = f"malformed Criterion estimates JSON in {estimates_json}: {err}"
        raise ValueError(msg) from err

    if not isinstance(data, dict):
        msg = f"expected JSON object in {estimates_json}"
        raise TypeError(msg)

    stat_obj = data.get(stat)
    if not isinstance(stat_obj, dict):
        msg = f"stat '{stat}' not found in {estimates_json}"
        raise KeyError(msg)

    point = _read_numeric_field(stat_obj, "point_estimate", estimates_json, stat)
    ci = stat_obj.get("confidence_interval")
    if not isinstance(ci, dict):
        return (point, point, point)

    lo = _read_numeric_field(ci, "lower_bound", estimates_json, stat, default=point)
    hi = _read_numeric_field(ci, "upper_bound", estimates_json, stat, default=point)
    return (point, lo, hi)


def _read_numeric_field(
    obj: dict[str, object],
    field: str,
    estimates_json: Path,
    stat: str,
    *,
    default: float | None = None,
) -> float:
    """Read a numeric Criterion field with file and statistic context."""
    if field not in obj:
        if default is not None:
            return default
        msg = f"field '{field}' for stat '{stat}' not found in {estimates_json}"
        raise KeyError(msg)
    value = obj[field]
    if isinstance(value, bool) or not isinstance(value, int | float | str):
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise TypeError(msg)
    try:
        return float(value)
    except ValueError as err:
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise ValueError(msg) from err


def _collect_exact_results(criterion_dir: Path, sample: str, stat: str) -> list[BenchResult]:
    """Collect exact-arithmetic benchmark results from Criterion output."""
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
            results.append(BenchResult(suite="exact", group=group, bench=bench, point_ns=point, ci_lo_ns=lo, ci_hi_ns=hi))

    return results


def _legacy_exact_baseline_bench(bench: str) -> str | None:
    """Return the legacy exact benchmark name for renamed rows."""
    legacy_bench = EXACT_LEGACY_BASELINE_BENCHES.get(bench)
    if legacy_bench is not None:
        return legacy_bench

    for current_prefix, legacy_prefix in _EXACT_LEGACY_PREFIX_BASELINE_BENCHES:
        if bench.startswith(current_prefix):
            return f"{legacy_prefix}{bench.removeprefix(current_prefix)}"

    return None


def _exact_baseline_path(group_dir: Path, bench: str, baseline_name: str) -> tuple[str, Path]:
    """Return the exact benchmark baseline path, falling back to legacy names."""
    base_path = group_dir / bench / baseline_name / "estimates.json"
    if base_path.exists():
        return (bench, base_path)

    legacy_bench = _legacy_exact_baseline_bench(bench)
    if legacy_bench is None:
        return (bench, base_path)

    legacy_path = group_dir / legacy_bench / baseline_name / "estimates.json"
    if legacy_path.exists():
        return (legacy_bench, legacy_path)

    return (bench, base_path)


def _ordered_vs_linalg_benches(group_dir: Path, sample: str) -> list[str]:
    """Return present vs_linalg benches in a stable, metric-aware order."""
    present = {child.name for child in group_dir.iterdir() if child.is_dir() and (child / sample / "estimates.json").exists()}
    ordered = [bench for bench in VS_LINALG_BENCH_ORDER if bench in present]
    extras = sorted(present.difference(VS_LINALG_BENCH_ORDER))
    return [*ordered, *extras]


def _collect_vs_linalg_results(criterion_dir: Path, sample: str, stat: str) -> list[BenchResult]:
    """Collect vs_linalg benchmark results from Criterion output."""
    results: list[BenchResult] = []
    dim_groups: list[tuple[int, Path]] = []

    for group_dir in criterion_dir.iterdir():
        if not group_dir.is_dir():
            continue
        dim = _dim_from_vs_linalg_group(group_dir.name)
        if dim is None:
            continue
        dim_groups.append((dim, group_dir))

    for _dim, group_dir in sorted(dim_groups, key=lambda item: item[0]):
        for bench in _ordered_vs_linalg_benches(group_dir, sample):
            est_path = group_dir / bench / sample / "estimates.json"
            point, lo, hi = _read_estimate(est_path, stat)
            results.append(BenchResult(suite="vs_linalg", group=group_dir.name, bench=bench, point_ns=point, ci_lo_ns=lo, ci_hi_ns=hi))

    return results


def _collect_results(criterion_dir: Path, sample: str, stat: str, suite: str = "all") -> list[BenchResult]:
    """Collect benchmark results from Criterion output."""
    results: list[BenchResult] = []
    if suite in ("all", "exact"):
        results.extend(_collect_exact_results(criterion_dir, sample, stat))
    if suite in ("all", "vs_linalg"):
        results.extend(_collect_vs_linalg_results(criterion_dir, sample, stat))
    return results


def _collect_exact_comparisons(
    criterion_dir: Path,
    baseline_name: str,
    stat: str,
    scope: str,
) -> list[Comparison]:
    """Compare current exact-arithmetic results against a named baseline."""
    comparisons: list[Comparison] = []

    for group, benches in EXACT_GROUPS.items():
        if scope == "release-signal" and group not in EXACT_RELEASE_SIGNAL_GROUPS:
            continue

        group_dir = criterion_dir / group
        if not group_dir.is_dir():
            continue

        for bench in benches:
            new_path = group_dir / bench / "new" / "estimates.json"
            baseline_bench, base_path = _exact_baseline_path(group_dir, bench, baseline_name)

            if not new_path.exists() or not base_path.exists():
                continue

            new_point, _, _ = _read_estimate(new_path, stat)
            base_point, _, _ = _read_estimate(base_path, stat)

            speedup = base_point / new_point if new_point > 0 else float("inf")
            pct_change = ((new_point - base_point) / base_point) * 100.0 if base_point > 0 else 0.0

            comparisons.append(
                Comparison(
                    suite="exact",
                    group=group,
                    bench=bench,
                    baseline_ns=base_point,
                    current_ns=new_point,
                    speedup=speedup,
                    pct_change=pct_change,
                    baseline_bench=baseline_bench if baseline_bench != bench else None,
                )
            )

    return comparisons


def _ordered_vs_linalg_comparison_benches(group_dir: Path, baseline_name: str, scope: str) -> list[str]:
    """Return present vs_linalg benches that have both current and baseline data."""
    present = {
        child.name
        for child in group_dir.iterdir()
        if child.is_dir() and (child / "new" / "estimates.json").exists() and (child / baseline_name / "estimates.json").exists()
    }
    if scope == "release-signal":
        present = present.intersection(VS_LINALG_LA_STACK_BENCHES)
    ordered = [bench for bench in VS_LINALG_BENCH_ORDER if bench in present]
    extras = sorted(present.difference(VS_LINALG_BENCH_ORDER))
    return [*ordered, *extras]


def _read_optional_point(estimates_json: Path, stat: str) -> float | None:
    """Read an optional Criterion point estimate."""
    if not estimates_json.exists():
        return None
    point, _, _ = _read_estimate(estimates_json, stat)
    return point


def _baseline_peer_times(group_dir: Path, bench: str, baseline_name: str, stat: str) -> tuple[float | None, float | None]:
    """Return last-release nalgebra/faer context for a la-stack vs_linalg bench."""
    peers = VS_LINALG_BASELINE_PEERS.get(bench)
    if peers is None:
        return (None, None)

    nalgebra_bench, faer_bench = peers
    nalgebra_ns = _read_optional_point(group_dir / nalgebra_bench / baseline_name / "estimates.json", stat)
    faer_ns = _read_optional_point(group_dir / faer_bench / baseline_name / "estimates.json", stat)
    return (nalgebra_ns, faer_ns)


def _comparison_bench_label(comparison: Comparison) -> str:
    """Return the display label for a comparison table row."""
    if comparison.baseline_bench is None:
        return comparison.bench
    return f"{comparison.bench} (vs {comparison.baseline_bench})"


def _collect_vs_linalg_comparisons(
    criterion_dir: Path,
    baseline_name: str,
    stat: str,
    scope: str,
) -> list[Comparison]:
    """Compare current vs_linalg results against a named baseline."""
    comparisons: list[Comparison] = []
    dim_groups: list[tuple[int, Path]] = []

    for group_dir in criterion_dir.iterdir():
        if not group_dir.is_dir():
            continue
        dim = _dim_from_vs_linalg_group(group_dir.name)
        if dim is None:
            continue
        dim_groups.append((dim, group_dir))

    for _dim, group_dir in sorted(dim_groups, key=lambda item: item[0]):
        for bench in _ordered_vs_linalg_comparison_benches(group_dir, baseline_name, scope):
            new_path = group_dir / bench / "new" / "estimates.json"
            base_path = group_dir / bench / baseline_name / "estimates.json"

            new_point, _, _ = _read_estimate(new_path, stat)
            base_point, _, _ = _read_estimate(base_path, stat)
            baseline_nalgebra_ns, baseline_faer_ns = _baseline_peer_times(group_dir, bench, baseline_name, stat)

            speedup = base_point / new_point if new_point > 0 else float("inf")
            pct_change = ((new_point - base_point) / base_point) * 100.0 if base_point > 0 else 0.0

            comparisons.append(
                Comparison(
                    suite="vs_linalg",
                    group=group_dir.name,
                    bench=bench,
                    baseline_ns=base_point,
                    current_ns=new_point,
                    speedup=speedup,
                    pct_change=pct_change,
                    baseline_nalgebra_ns=baseline_nalgebra_ns,
                    baseline_faer_ns=baseline_faer_ns,
                )
            )

    return comparisons


def _collect_comparisons(
    criterion_dir: Path,
    baseline_name: str,
    stat: str,
    suite: str = "all",
    scope: str = "release-signal",
) -> list[Comparison]:
    """Compare current (new) results against a named baseline."""
    comparisons: list[Comparison] = []
    if suite in ("all", "exact"):
        comparisons.extend(_collect_exact_comparisons(criterion_dir, baseline_name, stat, scope))
    if suite in ("all", "vs_linalg"):
        comparisons.extend(_collect_vs_linalg_comparisons(criterion_dir, baseline_name, stat, scope))
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


class _GroupedItem(Protocol):
    @property
    def suite(self) -> str: ...

    @property
    def group(self) -> str: ...


def _group_by_suite[T: _GroupedItem](items: list[T]) -> dict[str, list[T]]:
    """Group results or comparisons by benchmark suite."""
    suites: dict[str, list[T]] = {}
    for item in items:
        suites.setdefault(item.suite, []).append(item)
    return suites


def _group_by_group[T: _GroupedItem](items: list[T]) -> dict[str, list[T]]:
    """Group results or comparisons by their Criterion group name."""
    groups: dict[str, list[T]] = {}
    for item in items:
        groups.setdefault(item.group, []).append(item)
    return groups


def _suite_heading(suite: str) -> str:
    """Turn an internal suite name into a readable heading."""
    if suite == "exact":
        return "Exact arithmetic"
    if suite == "vs_linalg":
        return "vs_linalg"
    return suite


def _group_heading(group: str) -> str:
    """Turn a Criterion group name into a readable heading."""
    # exact_d3 -> "D=3", exact_random_percentile_d3 ->
    # "Random percentile D=3", exact_near_singular_3x3 ->
    # "Near-singular 3x3", exact_hilbert_4x4 -> "Hilbert 4x4", etc.
    if group.startswith("exact_random_percentile_d"):
        return f"Random percentile D={group.removeprefix('exact_random_percentile_d')}"
    if group.startswith("exact_d"):
        return f"D={group.removeprefix('exact_d')}"
    if group == "exact_near_singular_3x3":
        return "Near-singular 3x3"
    if group == "exact_large_entries_3x3":
        return "Large entries 3x3"
    if group.startswith("exact_hilbert_"):
        return f"Hilbert {group.removeprefix('exact_hilbert_')}"
    return group


def _group_heading_for_suite(suite: str, group: str) -> str:
    """Turn a Criterion group name into a readable heading for its suite."""
    if suite == "vs_linalg":
        dim = _dim_from_vs_linalg_group(group)
        if dim is not None:
            return f"D={dim}"
    return _group_heading(group)


def _snapshot_tables(results: list[BenchResult], stat: str) -> str:
    """Generate per-dimension markdown tables for a single set of results."""
    stat_label = stat.capitalize()
    sections: list[str] = []

    for suite, suite_items in _group_by_suite(results).items():
        sections.append(f"## {_suite_heading(suite)}")
        for group, items in _group_by_group(suite_items).items():
            lines = [
                f"### {_group_heading_for_suite(suite, group)}",
                "",
                f"| Benchmark | {stat_label} | 95% CI |",
                "|-----------|-------:|-------:|",
            ]
            for r in items:
                ci_range = f"[{_format_time(r.ci_lo_ns)}, {_format_time(r.ci_hi_ns)}]"
                lines.append(f"| {r.bench} | {_format_time(r.point_ns)} | {ci_range} |")
            sections.append("\n".join(lines))

    return "\n\n".join(sections)


def _comparison_tables(comparisons: list[Comparison], baseline_name: str) -> str:
    """Generate per-dimension markdown tables comparing baseline vs current."""
    sections: list[str] = []

    for suite, suite_items in _group_by_suite(comparisons).items():
        sections.append(f"## {_suite_heading(suite)}")
        for group, items in _group_by_group(suite_items).items():
            has_peer_context = any(item.baseline_nalgebra_ns is not None or item.baseline_faer_ns is not None for item in items)
            lines = [
                f"### {_group_heading_for_suite(suite, group)}",
                "",
            ]
            if has_peer_context:
                lines.extend(
                    [
                        f"| Benchmark | {baseline_name} | Latest | Change | Speedup | {baseline_name} nalgebra | {baseline_name} faer |",
                        "|-----------|-------:|-------:|-------:|--------:|-------:|-------:|",
                    ]
                )
            else:
                lines.extend(
                    [
                        f"| Benchmark | {baseline_name} | Latest | Change | Speedup |",
                        "|-----------|-------:|-------:|-------:|--------:|",
                    ]
                )
            for c in items:
                cells = [
                    _comparison_bench_label(c),
                    _format_time(c.baseline_ns),
                    _format_time(c.current_ns),
                    _format_pct(c.pct_change),
                    f"{c.speedup:.2f}x",
                ]
                if has_peer_context:
                    cells.extend(
                        [
                            _format_time(c.baseline_nalgebra_ns) if c.baseline_nalgebra_ns is not None else "",
                            _format_time(c.baseline_faer_ns) if c.baseline_faer_ns is not None else "",
                        ]
                    )
                lines.append(f"| {' | '.join(cells)} |")
            sections.append("\n".join(lines))

    return "\n\n".join(sections)


def _read_cargo_version(root: Path) -> str:
    cargo_toml = root / "Cargo.toml"
    if not cargo_toml.exists():
        return "unknown"
    data = tomllib.loads(cargo_toml.read_text(encoding="utf-8"))
    package = data.get("package")
    if isinstance(package, dict):
        version = package.get("version")
        if isinstance(version, str):
            return version
    return "unknown"


def _get_git_info(root: Path) -> tuple[str, str]:
    """Return (short_hash, branch_or_tag)."""
    short_hash = "unknown"
    branch = "unknown"
    try:
        result = run_git_command(["--no-pager", "rev-parse", "--short", "HEAD"], cwd=root)
        short_hash = result.stdout.strip()
    except (ExecutableNotFoundError, subprocess.CalledProcessError):
        pass
    try:
        result = run_git_command(["--no-pager", "rev-parse", "--abbrev-ref", "HEAD"], cwd=root)
        branch = result.stdout.strip()
    except (ExecutableNotFoundError, subprocess.CalledProcessError):
        pass
    return short_hash, branch


def _generate_markdown(
    root: Path,
    table: str,
    settings: ReportSettings,
) -> str:
    """Generate the complete benchmark report content."""
    version = _read_cargo_version(root)
    short_hash, branch = _get_git_info(root)
    now = datetime.now(tz=UTC).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = [
        "# Benchmark Performance",
        "",
        f"**la-stack** v{version} · `{short_hash}` ({branch}) · {now}",
        f"**Statistic**: {settings.stat}",
        f"**Suite**: {settings.suite}",
        f"**Scope**: {settings.scope}",
        "",
        "## Benchmark Results",
        "",
    ]

    if settings.baseline_name:
        lines.append(f"Comparison against baseline **{settings.baseline_name}**:")
        lines.append("")
        lines.append("Negative change = faster. Speedup > 1.00x = improvement.")
    else:
        lines.append("Current performance snapshot (no baseline comparison).")

    lines.extend(["", table, ""])

    lines.extend(
        [
            "## How to Update",
            "",
            "Local performance reports are generated in isolated temporary worktrees:",
            "",
            "```bash",
            "# Local development: compare the current tree with the latest release",
            "just performance-local",
            "",
            "# Release PR: update docs/PERFORMANCE.md and archive the previous report",
            "just performance-release",
            "",
            "# GitHub Actions release assets",
            "just performance-github-assets",
            "",
            "# Explicit repair",
            "just performance-release <current-tag> <previous-tag>",
            "```",
            "",
            "`just performance-local` writes `target/bench-reports/performance.md`.",
            "`just performance-github-assets` writes `target/bench-reports/github-assets-performance.md`.",
            "",
            "Older curated release-to-release reports are archived in `docs/archive/performance/`.",
            "",
            "See `docs/BENCHMARKING.md` for the full comparison workflow.",
        ]
    )

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare la-stack benchmark results across Criterion baselines.",
    )
    parser.add_argument(
        "baseline",
        nargs="?",
        default="last",
        help="Baseline name to compare against (default: 'last').",
    )
    parser.add_argument(
        "--snapshot",
        action="store_true",
        help="Generate a current-performance snapshot instead of comparing against a baseline.",
    )
    parser.add_argument(
        "--stat",
        default="median",
        choices=["mean", "median"],
        help="Statistic to compare (default: median).",
    )
    parser.add_argument(
        "--suite",
        default="all",
        choices=SUITE_CHOICES,
        help="Benchmark suite to compare (default: all).",
    )
    parser.add_argument(
        "--scope",
        default="release-signal",
        choices=SCOPE_CHOICES,
        help="Comparison scope: release-signal compares la-stack latest against last with baseline peer context; all-benches compares every present bench.",
    )
    parser.add_argument(
        "--criterion-dir",
        default="target/criterion",
        help="Criterion output directory (default: target/criterion).",
    )
    parser.add_argument(
        "--output",
        default="target/bench-reports/performance.md",
        help="Output markdown file (default: target/bench-reports/performance.md).",
    )
    return parser.parse_args(argv)


def _run_bench_hint(suite: str) -> str:
    if suite == "exact":
        return "just bench-exact"
    if suite == "vs_linalg":
        return "just bench-vs-linalg-la-stack"
    return "just bench-latest"


def _save_baseline_hint(suite: str, baseline: str) -> str:
    if suite == "exact":
        return f"just bench-save-baseline {baseline} exact"
    if suite == "vs_linalg":
        return f"just bench-save-baseline {baseline} vs_linalg"
    return f"just bench-save-baseline {baseline}"


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    root = _repo_root()
    criterion_dir = root / args.criterion_dir
    output_path = Path(args.output) if Path(args.output).is_absolute() else root / args.output

    if not criterion_dir.is_dir():
        print(
            f"No Criterion results found at {criterion_dir}.\nRun benchmarks first:\n  {_run_bench_hint(args.suite)}\n",
            file=sys.stderr,
        )
        return 2

    baseline_name = None if args.snapshot else args.baseline

    if baseline_name:
        comparisons = _collect_comparisons(criterion_dir, baseline_name, args.stat, args.suite, args.scope)
        if not comparisons:
            print(
                f"No comparison data found for baseline '{baseline_name}'.\n"
                f"Save a baseline first:\n  {_save_baseline_hint(args.suite, baseline_name)}\n"
                f"Then run benchmarks:\n  {_run_bench_hint(args.suite)}\n",
                file=sys.stderr,
            )
            return 2
        table = _comparison_tables(comparisons, baseline_name)
    else:
        results = _collect_results(criterion_dir, "new", args.stat, args.suite)
        if not results:
            print(
                f"No benchmark results found.\nRun benchmarks first:\n  {_run_bench_hint(args.suite)}\n",
                file=sys.stderr,
            )
            return 2
        table = _snapshot_tables(results, args.stat)

    settings = ReportSettings(
        baseline_name=baseline_name,
        stat=args.stat,
        suite=args.suite,
        scope=args.scope,
    )
    md = _generate_markdown(root, table, settings)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(md, encoding="utf-8")
    print(f"📊 Wrote {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
