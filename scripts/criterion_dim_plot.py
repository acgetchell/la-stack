#!/usr/bin/env python3
"""Aggregate Criterion benchmark results into a time-vs-dimension chart.

Reads Criterion output under:
  target/criterion/d{D}/{benchmark}/{new|base}/estimates.json

And writes:
  docs/assets/bench/vs_linalg_{metric}_{stat}.csv
  docs/assets/bench/vs_linalg_{metric}_{stat}.svg

This is intended to create a single, README-friendly plot comparing la-stack to other
Rust linear algebra crates across dimensions.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final


@dataclass(frozen=True, slots=True)
class Metric:
    la_bench: str
    na_bench: str
    fa_bench: str
    title: str


@dataclass(frozen=True, slots=True)
class PlotRequest:
    csv_path: Path
    out_svg: Path
    title: str
    stat: str
    dims: tuple[int, ...]
    log_y: bool


@dataclass(frozen=True, slots=True)
class Row:
    dim: int
    la_time: float
    la_lo: float
    la_hi: float
    na_time: float
    na_lo: float
    na_hi: float
    fa_time: float
    fa_lo: float
    fa_hi: float


class ReadmeMarkerError(ValueError):
    """Raised when README markers are missing, duplicated, or out of order."""


METRICS: Final[dict[str, Metric]] = {
    "det_via_lu": Metric(
        la_bench="la_stack_det_via_lu",
        na_bench="nalgebra_det_via_lu",
        fa_bench="faer_det_via_lu",
        title="Determinant via LU (factor + det)",
    ),
    "lu": Metric(
        la_bench="la_stack_lu",
        na_bench="nalgebra_lu",
        fa_bench="faer_lu",
        title="LU factorization",
    ),
    "lu_solve": Metric(
        la_bench="la_stack_lu_solve",
        na_bench="nalgebra_lu_solve",
        fa_bench="faer_lu_solve",
        title="LU solve (factor + solve)",
    ),
    "solve_from_lu": Metric(
        la_bench="la_stack_solve_from_lu",
        na_bench="nalgebra_solve_from_lu",
        fa_bench="faer_solve_from_lu",
        title="Solve from precomputed LU",
    ),
    "det_from_lu": Metric(
        la_bench="la_stack_det_from_lu",
        na_bench="nalgebra_det_from_lu",
        fa_bench="faer_det_from_lu",
        title="Determinant from precomputed LU",
    ),
    "dot": Metric(
        la_bench="la_stack_dot",
        na_bench="nalgebra_dot",
        fa_bench="faer_dot",
        title="Vector dot product",
    ),
    # Different names between crates.
    "norm2_sq": Metric(
        la_bench="la_stack_norm2_sq",
        na_bench="nalgebra_norm_squared",
        fa_bench="faer_norm2_sq",
        title="Vector squared 2-norm",
    ),
    "inf_norm": Metric(
        la_bench="la_stack_inf_norm",
        na_bench="nalgebra_inf_norm",
        fa_bench="faer_inf_norm",
        title="Matrix infinity norm (max abs row sum)",
    ),
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _dim_from_group_dir(name: str) -> int | None:
    match = re.fullmatch(r"d(\d+)", name)
    if match is None:
        return None
    return int(match.group(1))


def _discover_dims(criterion_dir: Path) -> list[int]:
    dims: list[int] = []
    for child in criterion_dir.iterdir():
        if not child.is_dir():
            continue
        d = _dim_from_group_dir(child.name)
        if d is None:
            continue
        dims.append(d)
    return sorted(dims)


def _read_estimate(estimates_json: Path, stat: str) -> tuple[float, float, float]:
    data = json.loads(estimates_json.read_text(encoding="utf-8"))

    stat_obj = data.get(stat)
    if not isinstance(stat_obj, dict):
        raise KeyError(f"stat '{stat}' not found in {estimates_json}")

    point = float(stat_obj["point_estimate"])
    ci = stat_obj.get("confidence_interval")
    if not isinstance(ci, dict):
        return (point, point, point)

    lo = float(ci.get("lower_bound", point))
    hi = float(ci.get("upper_bound", point))
    return (point, lo, hi)


def _write_csv(out_csv: Path, rows: list[Row]) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8") as f:
        f.write("D,la_stack,la_lo,la_hi,nalgebra,na_lo,na_hi,faer,fa_lo,fa_hi\n")
        for row in rows:
            f.write(f"{row.dim},{row.la_time},{row.la_lo},{row.la_hi},{row.na_time},{row.na_lo},{row.na_hi},{row.fa_time},{row.fa_lo},{row.fa_hi}\n")


def _pct_reduction(baseline: float, value: float) -> str:
    """Percent time reduction relative to baseline (positive = value is faster)."""
    if baseline == 0.0:
        return "n/a"
    pct = ((baseline - value) / baseline) * 100.0
    return f"{pct:+.1f}%"


def _markdown_table(rows: list[Row], stat: str) -> str:
    lines = [
        f"| D | la-stack {stat} (ns) | nalgebra {stat} (ns) | faer {stat} (ns) | la-stack vs nalgebra | la-stack vs faer |",
        "|---:|--------------------:|--------------------:|----------------:|---------------------:|----------------:|",
    ]

    for row in rows:
        pct_vs_na = _pct_reduction(row.na_time, row.la_time)
        pct_vs_fa = _pct_reduction(row.fa_time, row.la_time)
        lines.append(f"| {row.dim} | {row.la_time:,.3f} | {row.na_time:,.3f} | {row.fa_time:,.3f} | {pct_vs_na} | {pct_vs_fa} |")

    return "\n".join(lines)


def _readme_table_markers(metric: str, stat: str, sample: str) -> tuple[str, str]:
    tag = f"BENCH_TABLE:{metric}:{stat}:{sample}"
    return (f"<!-- {tag}:BEGIN -->", f"<!-- {tag}:END -->")


def _update_readme_table(readme_path: Path, marker_begin: str, marker_end: str, table_md: str) -> bool:
    lines = readme_path.read_text(encoding="utf-8").splitlines(keepends=True)

    begin_indices = [i for i, line in enumerate(lines) if line.strip() == marker_begin]
    end_indices = [i for i, line in enumerate(lines) if line.strip() == marker_end]

    if len(begin_indices) != 1 or len(end_indices) != 1:
        raise ReadmeMarkerError(f"README markers not found or not unique. Expected exactly one of each:\n  {marker_begin}\n  {marker_end}\n")

    begin_idx = begin_indices[0]
    end_idx = end_indices[0]
    if begin_idx >= end_idx:
        msg = "README markers are out of order."
        raise ReadmeMarkerError(msg)

    table_lines = [line + "\n" for line in table_md.strip("\n").splitlines()]
    new_lines = [
        *lines[: begin_idx + 1],
        *table_lines,
        *lines[end_idx:],
    ]

    if new_lines == lines:
        return False

    readme_path.write_text("".join(new_lines), encoding="utf-8")
    return True


def _gp_quote(s: str) -> str:
    # gnuplot supports single-quoted strings; escape backslashes and single quotes.
    return "'" + s.replace("\\", "\\\\").replace("'", "\\'") + "'"


def _render_svg_with_gnuplot(req: PlotRequest) -> None:
    gnuplot_path = shutil.which("gnuplot")
    if gnuplot_path is None:
        msg = "gnuplot not found. Install it (macOS: `brew install gnuplot`) or re-run with --no-plot."
        raise FileNotFoundError(msg)

    req.out_svg.parent.mkdir(parents=True, exist_ok=True)

    xtics = ", ".join(str(d) for d in req.dims)

    gp_lines = [
        "set terminal svg size 960,540 noenhanced",
        f"set output {_gp_quote(str(req.out_svg))}",
        "set datafile separator comma",
        "set grid",
        "set key left top",
        f"set title {_gp_quote(req.title)}",
        "set xlabel 'Dimension D'",
        f"set ylabel {_gp_quote(f'{req.stat} time (ns)')}",
        f"set xtics ({xtics})",
        "set style line 1 lc rgb '#1f77b4' lt 1 lw 2 pt 7 ps 1",
        "set style line 2 lc rgb '#ff7f0e' lt 1 lw 2 pt 5 ps 1",
        "set style line 3 lc rgb '#2ca02c' lt 1 lw 2 pt 9 ps 1",
        "set style data linespoints",
        "set tics nomirror",
        "set border linewidth 1",
    ]

    if req.log_y:
        gp_lines.append("set logscale y 10")

    gp_lines.extend(
        [
            "plot \\",
            f"  {_gp_quote(str(req.csv_path))} using 1:2:3:4 with yerrorlines ls 1 title 'la-stack', \\",
            f"  {_gp_quote(str(req.csv_path))} using 1:5:6:7 with yerrorlines ls 2 title 'nalgebra', \\",
            f"  {_gp_quote(str(req.csv_path))} using 1:8:9:10 with yerrorlines ls 3 title 'faer'",
        ]
    )

    # Safe: gnuplot executable is resolved via PATH; input is a generated script with fully
    # quoted file paths.
    subprocess.run([gnuplot_path], input="\n".join(gp_lines), text=True, check=True)  # noqa: S603


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Criterion time vs dimension for la-stack vs nalgebra/faer.")

    parser.add_argument(
        "--metric",
        default="lu_solve",
        choices=sorted(METRICS.keys()),
        help="Which vs_linalg metric to plot.",
    )
    parser.add_argument(
        "--stat",
        default="median",
        choices=["mean", "median"],
        help="Statistic to plot from estimates.json.",
    )
    parser.add_argument(
        "--sample",
        default="new",
        choices=["new", "base"],
        help="Which Criterion run directory to read (new = most recent).",
    )
    parser.add_argument(
        "--criterion-dir",
        default="target/criterion",
        help="Criterion output directory (default: target/criterion).",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output SVG path (default: docs/assets/bench/vs_linalg_{metric}_{stat}.svg).",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Output CSV path (default: docs/assets/bench/vs_linalg_{metric}_{stat}.csv).",
    )
    parser.add_argument(
        "--log-y",
        action="store_true",
        help="Use a log-scale y-axis.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Only write CSV (skip gnuplot/SVG).",
    )
    parser.add_argument(
        "--update-readme",
        action="store_true",
        help="Update a Markdown table in README.md between BENCH_TABLE markers.",
    )
    parser.add_argument(
        "--readme",
        default="README.md",
        help="Path to README file to update (default: README.md at repo root).",
    )

    return parser.parse_args(argv)


def _resolve_under_root(root: Path, arg: str) -> Path:
    path = Path(arg)
    return path if path.is_absolute() else root / path


def _resolve_output_paths(root: Path, metric: str, stat: str, out_svg: str | None, out_csv: str | None) -> tuple[Path, Path]:
    svg = Path(out_svg) if out_svg is not None else Path(f"docs/assets/bench/vs_linalg_{metric}_{stat}.svg")
    csv = Path(out_csv) if out_csv is not None else Path(f"docs/assets/bench/vs_linalg_{metric}_{stat}.csv")

    if not svg.is_absolute():
        svg = root / svg
    if not csv.is_absolute():
        csv = root / csv

    return (svg, csv)


def _collect_rows(criterion_dir: Path, dims: list[int], metric: Metric, stat: str, sample: str) -> tuple[list[Row], list[str]]:
    rows: list[Row] = []
    skipped: list[str] = []

    for d in dims:
        group_dir = criterion_dir / f"d{d}"
        la_est = group_dir / metric.la_bench / sample / "estimates.json"
        na_est = group_dir / metric.na_bench / sample / "estimates.json"
        fa_est = group_dir / metric.fa_bench / sample / "estimates.json"

        if not la_est.exists() or not na_est.exists() or not fa_est.exists():
            skipped.append(f"d{d} (missing {metric.la_bench}, {metric.na_bench}, or {metric.fa_bench})")
            continue

        la, la_lo, la_hi = _read_estimate(la_est, stat)
        na, na_lo, na_hi = _read_estimate(na_est, stat)
        fa, fa_lo, fa_hi = _read_estimate(fa_est, stat)
        rows.append(
            Row(
                dim=d,
                la_time=la,
                la_lo=la_lo,
                la_hi=la_hi,
                na_time=na,
                na_lo=na_lo,
                na_hi=na_hi,
                fa_time=fa,
                fa_lo=fa_lo,
                fa_hi=fa_hi,
            )
        )

    return (rows, skipped)


def _maybe_update_readme(root: Path, args: argparse.Namespace, rows: list[Row]) -> int:
    if not args.update_readme:
        return 0

    readme_path = _resolve_under_root(root, args.readme)

    marker_begin, marker_end = _readme_table_markers(args.metric, args.stat, args.sample)
    table_md = _markdown_table(rows, args.stat)

    try:
        changed = _update_readme_table(readme_path, marker_begin, marker_end, table_md)
    except (OSError, ValueError) as e:
        print(str(e), file=sys.stderr)
        return 2

    if changed:
        print(f"Updated README table: {readme_path}")

    return 0


def _maybe_render_plot(args: argparse.Namespace, req: PlotRequest, skipped: list[str]) -> int:
    if args.no_plot:
        print(f"Wrote CSV: {req.csv_path}")
        return 0

    try:
        _render_svg_with_gnuplot(req)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(str(e), file=sys.stderr)
        print(f"Wrote CSV instead: {req.csv_path}", file=sys.stderr)
        return 1

    if skipped:
        print("Warning: some dimension groups were skipped:")
        for s in skipped:
            print(f"  - {s}")

    print(f"Wrote CSV: {req.csv_path}")
    print(f"Wrote SVG: {req.out_svg}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    root = _repo_root()

    criterion_dir = _resolve_under_root(root, args.criterion_dir)

    dims = _discover_dims(criterion_dir) if criterion_dir.exists() else []
    if not dims:
        print(
            f"No Criterion results found under {criterion_dir}.\n\nRun benchmarks first, e.g.:\n  cargo bench --bench vs_linalg\n",
            file=sys.stderr,
        )
        return 2

    metric = METRICS[args.metric]

    out_svg, out_csv = _resolve_output_paths(root, args.metric, args.stat, args.out, args.csv)

    rows, skipped = _collect_rows(criterion_dir, dims, metric, args.stat, args.sample)
    if not rows:
        print(
            "No benchmark results found to plot for the selected metric/stat.\n"
            f"Expected files like:\n  {criterion_dir}/d32/{metric.la_bench}/{args.sample}/estimates.json\n",
            file=sys.stderr,
        )
        if skipped:
            print("Skipped groups:", *skipped, sep="\n  - ", file=sys.stderr)
        return 2

    _write_csv(out_csv, rows)

    rc = _maybe_update_readme(root, args, rows)
    if rc != 0:
        return rc

    dims_present = [row.dim for row in rows]

    title = f"{metric.title}: {args.stat} time vs dimension"
    req = PlotRequest(
        csv_path=out_csv,
        out_svg=out_svg,
        title=title,
        stat=args.stat,
        dims=tuple(dims_present),
        log_y=args.log_y,
    )

    return _maybe_render_plot(args, req, skipped)


if __name__ == "__main__":
    raise SystemExit(main())
