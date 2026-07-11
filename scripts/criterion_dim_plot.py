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
import hashlib
import json
import math
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Protocol, TypeGuard, cast

from subprocess_utils import ExecutableNotFoundError, run_git_command, run_safe_command


@dataclass(frozen=True, slots=True)
class Metric:
    """Criterion benchmark names and display title for one plotted metric."""

    la_bench: str
    na_bench: str
    fa_bench: str
    title: str


@dataclass(frozen=True, slots=True)
class PlotRequest:
    """Validated inputs required to render a benchmark SVG."""

    csv_path: Path
    out_svg: Path
    title: str
    stat: str
    dims: tuple[int, ...]
    la_label: str
    na_label: str
    fa_label: str
    log_y: bool


@dataclass(frozen=True, slots=True)
class PlotCliArgs:
    """Validated command-line options for the plot generator."""

    metric: str
    stat: str
    sample: str
    criterion_dir: str
    out: str | None
    csv: str | None
    log_y: bool
    no_plot: bool
    update_readme: bool
    readme: str
    allow_partial: bool


@dataclass(frozen=True, slots=True)
class Row:
    """Validated timing estimates for one benchmark dimension."""

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

    def __post_init__(self) -> None:
        """Reject invalid dimensions, timings, and confidence intervals."""
        if self.dim <= 0:
            msg = f"dimension must be positive: {self.dim}"
            raise ValueError(msg)
        for field, value in (
            ("la_time", self.la_time),
            ("la_lo", self.la_lo),
            ("la_hi", self.la_hi),
            ("na_time", self.na_time),
            ("na_lo", self.na_lo),
            ("na_hi", self.na_hi),
            ("fa_time", self.fa_time),
            ("fa_lo", self.fa_lo),
            ("fa_hi", self.fa_hi),
        ):
            _require_nonnegative_finite_time(value, field)
        _require_confidence_interval(self.la_lo, self.la_time, self.la_hi, "la_stack row")
        _require_confidence_interval(self.na_lo, self.na_time, self.na_hi, "nalgebra row")
        _require_confidence_interval(self.fa_lo, self.fa_time, self.fa_hi, "faer row")


class ReadmeMarkerError(ValueError):
    """Base error for invalid README BENCH_TABLE markers."""


class MarkerNotFoundError(ReadmeMarkerError):
    """Raised when README markers are missing or not unique."""


class MarkerOrderError(ReadmeMarkerError):
    """Raised when README markers are out of order."""


class _ReadmeArgs(Protocol):
    @property
    def update_readme(self) -> bool: ...

    @property
    def readme(self) -> str: ...

    @property
    def metric(self) -> str: ...

    @property
    def stat(self) -> str: ...

    @property
    def sample(self) -> str: ...


class _RenderArgs(Protocol):
    @property
    def no_plot(self) -> bool: ...


type ParsedObject = dict[str, object]


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

CANONICAL_DIMS: Final[tuple[int, ...]] = (2, 3, 4, 5, 8, 16, 32, 64)
_PUBLICATION_GATE: Final[tuple[str, ...]] = ("just", "test-bench-inputs")
_PUBLICATION_BENCHMARK: Final[tuple[str, ...]] = (
    "cargo",
    "bench",
    "--locked",
    "--features",
    "bench",
    "--bench",
    "vs_linalg",
)
_COMMAND_TIMEOUT_SECONDS: Final[int] = 7200
_PROVENANCE_HARNESS_FILES: Final[tuple[str, ...]] = (
    "Cargo.toml",
    "Cargo.lock",
    "rust-toolchain.toml",
    "justfile",
    "tests/exact_bench_config.rs",
    "tests/vs_linalg_inputs.rs",
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _dim_from_group_dir(name: str) -> int | None:
    match = re.fullmatch(r"d(\d+)", name)
    if match is None:
        return None
    return int(match.group(1))


def _is_parsed_object(value: object) -> TypeGuard[ParsedObject]:
    """Return true when a parsed JSON/TOML value is an object with string keys."""
    return isinstance(value, dict) and all(isinstance(key, str) for key in value)


def _require_parsed_object(value: object, context: str) -> ParsedObject:
    if not _is_parsed_object(value):
        msg = f"expected object for {context}"
        raise TypeError(msg)
    return value


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


def _read_cargo_package_version(cargo_toml: Path) -> str | None:
    if not cargo_toml.exists():
        return None

    data = _read_cargo_toml(cargo_toml)
    package = data.get("package")
    if _is_parsed_object(package):
        version = package.get("version")
        if isinstance(version, str):
            return version
    return None


def _read_cargo_dependency_versions(cargo_toml: Path, names: set[str]) -> dict[str, str]:
    if not cargo_toml.exists():
        return {}

    data = _read_cargo_toml(cargo_toml)
    versions: dict[str, str] = {}
    for section in ("dependencies", "dev-dependencies", "build-dependencies"):
        table = data.get(section)
        if not _is_parsed_object(table):
            continue
        for name in names:
            value = table.get(name)
            if isinstance(value, str):
                versions[name] = value
            elif _is_parsed_object(value):
                version = value.get("version")
                if isinstance(version, str):
                    versions[name] = version

    return versions


def _read_cargo_toml(cargo_toml: Path) -> ParsedObject:
    data: object = tomllib.loads(cargo_toml.read_text(encoding="utf-8"))
    return _require_parsed_object(data, str(cargo_toml))


def _detect_versions(root: Path) -> dict[str, str]:
    cargo_toml = root / "Cargo.toml"
    package_version = _read_cargo_package_version(cargo_toml) or "unknown"
    dep_versions = _read_cargo_dependency_versions(cargo_toml, {"nalgebra", "faer"})

    return {
        "la-stack": package_version,
        "nalgebra": dep_versions.get("nalgebra", "unknown"),
        "faer": dep_versions.get("faer", "unknown"),
    }


def _print_versions(versions: dict[str, str]) -> None:
    order = ["la-stack", "nalgebra", "faer"]
    text = ", ".join(f"{name}={versions.get(name, 'unknown')}" for name in order)
    print(f"Detected crate versions for legend: {text}", file=sys.stderr)


def _format_legend_label(name: str, version: str) -> str:
    if version == "unknown":
        return name
    return f"{name} v{version}"


def _read_estimate(estimates_json: Path, stat: str) -> tuple[float, float, float]:
    data = _read_json_object(estimates_json)

    stat_obj = data.get(stat)
    if not _is_parsed_object(stat_obj):
        raise KeyError(f"stat '{stat}' not found in {estimates_json}")

    point = _read_numeric_field(stat_obj, "point_estimate", estimates_json, stat)
    ci = stat_obj.get("confidence_interval")
    if not _is_parsed_object(ci):
        return (point, point, point)

    lo = _read_numeric_field(ci, "lower_bound", estimates_json, stat, default=point)
    hi = _read_numeric_field(ci, "upper_bound", estimates_json, stat, default=point)
    _require_confidence_interval(lo, point, hi, f"{stat}.confidence_interval in {estimates_json}")
    return (point, lo, hi)


def _read_json_object(path: Path) -> ParsedObject:
    try:
        data: object = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        msg = f"malformed Criterion estimates JSON in {path}: {err}"
        raise ValueError(msg) from err
    return _require_parsed_object(data, str(path))


def _read_numeric_field(
    obj: ParsedObject,
    field: str,
    estimates_json: Path,
    stat: str,
    *,
    default: float | None = None,
) -> float:
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
        parsed = float(value)
    except ValueError as err:
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise ValueError(msg) from err
    return _require_nonnegative_finite_time(parsed, f"{stat}.{field} in {estimates_json}")


def _require_nonnegative_finite_time(value: float, context: str) -> float:
    if not math.isfinite(value) or value < 0.0:
        msg = f"{context} must be finite and nonnegative: {value!r}"
        raise ValueError(msg)
    return value


def _require_confidence_interval(lo: float, point: float, hi: float, context: str) -> None:
    if lo > hi:
        msg = f"{context} lower bound must be <= upper bound: {lo!r} > {hi!r}"
        raise ValueError(msg)
    if not lo <= point <= hi:
        msg = f"{context} point estimate must be inside confidence interval: {lo!r} <= {point!r} <= {hi!r}"
        raise ValueError(msg)


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
        f"| D | la-stack {stat} (ns) | nalgebra {stat} (ns) | faer {stat} (ns) | "
        "la-stack point-estimate reduction vs nalgebra | la-stack point-estimate reduction vs faer |",
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
        msg = f"README markers not found or not unique (begin={len(begin_indices)}, end={len(end_indices)})."
        raise MarkerNotFoundError(msg)

    begin_idx = begin_indices[0]
    end_idx = end_indices[0]
    if begin_idx >= end_idx:
        msg = "README markers are out of order."
        raise MarkerOrderError(msg)

    table_lines = ["\n", *[line + "\n" for line in table_md.strip("\n").splitlines()], "\n"]
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
            f"  {_gp_quote(str(req.csv_path))} using 1:2:3:4 with yerrorlines ls 1 title {_gp_quote(req.la_label)}, \\",
            f"  {_gp_quote(str(req.csv_path))} using 1:5:6:7 with yerrorlines ls 2 title {_gp_quote(req.na_label)}, \\",
            f"  {_gp_quote(str(req.csv_path))} using 1:8:9:10 with yerrorlines ls 3 title {_gp_quote(req.fa_label)}",
        ]
    )

    try:
        run_safe_command("gnuplot", [], input="\n".join(gp_lines))
    except ExecutableNotFoundError as exc:
        msg = "gnuplot not found. Install it (macOS: `brew install gnuplot`) or re-run with --no-plot."
        raise FileNotFoundError(msg) from exc


def _parse_args(argv: list[str]) -> PlotCliArgs:
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
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Allow incomplete dimensions for exploratory CSV/SVG output; incompatible with --update-readme.",
    )

    args = parser.parse_args(argv)
    return PlotCliArgs(
        metric=_required_str_attr(args, "metric"),
        stat=_required_str_attr(args, "stat"),
        sample=_required_str_attr(args, "sample"),
        criterion_dir=_required_str_attr(args, "criterion_dir"),
        out=_optional_str_attr(args, "out"),
        csv=_optional_str_attr(args, "csv"),
        log_y=_required_bool_attr(args, "log_y"),
        no_plot=_required_bool_attr(args, "no_plot"),
        update_readme=_required_bool_attr(args, "update_readme"),
        readme=_required_str_attr(args, "readme"),
        allow_partial=_required_bool_attr(args, "allow_partial"),
    )


def _required_str_attr(args: argparse.Namespace, name: str) -> str:
    value = getattr(args, name)
    if not isinstance(value, str):
        msg = f"argparse returned non-string value for {name}: {value!r}"
        raise TypeError(msg)
    return value


def _optional_str_attr(args: argparse.Namespace, name: str) -> str | None:
    value = getattr(args, name)
    if value is None or isinstance(value, str):
        return value
    msg = f"argparse returned non-string value for {name}: {value!r}"
    raise TypeError(msg)


def _required_bool_attr(args: argparse.Namespace, name: str) -> bool:
    value = getattr(args, name)
    if not isinstance(value, bool):
        msg = f"argparse returned non-bool value for {name}: {value!r}"
        raise TypeError(msg)
    return value


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

        missing = [
            bench
            for bench, path in (
                (metric.la_bench, la_est),
                (metric.na_bench, na_est),
                (metric.fa_bench, fa_est),
            )
            if not path.exists()
        ]
        if missing:
            skipped.append(f"d{d} (missing {', '.join(missing)})")
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


def _run_publication_benchmarks(root: Path) -> None:
    """Validate fixtures, then produce fresh README publication measurements."""
    _run_publication_command(root, _PUBLICATION_GATE)
    criterion_dir = root / "target" / "criterion"
    with tempfile.TemporaryDirectory(prefix="la-stack-stale-criterion-") as tmp:
        backup_root = Path(tmp)
        moved = _stage_existing_new_samples(criterion_dir, backup_root)
        try:
            _run_publication_command(root, _PUBLICATION_BENCHMARK)
        except RuntimeError:
            _remove_vs_linalg_new_samples(criterion_dir)
            for relative in moved:
                source = backup_root / relative
                destination = criterion_dir / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                source.replace(destination)
            raise


def _run_publication_command(root: Path, command: tuple[str, ...]) -> None:
    """Run one publication command with complete failure context."""
    try:
        run_safe_command(
            command[0],
            list(command[1:]),
            cwd=root,
            timeout=_COMMAND_TIMEOUT_SECONDS,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if isinstance(exc.stderr, str) else ""
        detail = f"\nstderr:\n{stderr}" if stderr else ""
        msg = f"publication command failed ({exc.returncode}): {' '.join(command)}{detail}"
        raise RuntimeError(msg) from exc
    except ExecutableNotFoundError as exc:
        msg = f"publication command could not start: {' '.join(command)}: {exc}"
        raise RuntimeError(msg) from exc
    except subprocess.TimeoutExpired as exc:
        stderr = exc.stderr.strip() if isinstance(exc.stderr, str) else ""
        detail = f"\nstderr:\n{stderr}" if stderr else ""
        msg = f"publication command timed out after {exc.timeout} seconds: {' '.join(command)}{detail}"
        raise RuntimeError(msg) from exc


def _vs_linalg_new_samples(criterion_dir: Path) -> list[Path]:
    """Return every current vs_linalg sample in deterministic order."""
    if not criterion_dir.is_dir():
        return []
    return sorted(
        (
            sample
            for group in criterion_dir.iterdir()
            if group.is_dir() and _dim_from_group_dir(group.name) is not None
            for sample in group.glob("*/new")
            if sample.is_dir()
        ),
        key=lambda path: path.relative_to(criterion_dir).as_posix(),
    )


def _stage_existing_new_samples(criterion_dir: Path, backup_root: Path) -> list[Path]:
    """Move stale `new` samples aside so only the fresh timing run can satisfy coverage."""
    moved: list[Path] = []
    for sample in _vs_linalg_new_samples(criterion_dir):
        relative = sample.relative_to(criterion_dir)
        backup = backup_root / relative
        backup.parent.mkdir(parents=True, exist_ok=True)
        sample.replace(backup)
        moved.append(relative)
    return moved


def _remove_vs_linalg_new_samples(criterion_dir: Path) -> None:
    """Remove partial current samples before restoring a failed publication run."""
    for sample in _vs_linalg_new_samples(criterion_dir):
        shutil.rmtree(sample)


def _git_value(root: Path, args: list[str]) -> str:
    """Return deterministic Git provenance or an explicit unavailable label."""
    try:
        value = run_git_command(args, cwd=root).stdout.strip()
    except ExecutableNotFoundError, subprocess.CalledProcessError:
        return "unavailable"
    return value or "unavailable"


def _git_status_metadata(root: Path) -> tuple[bool | None, str]:
    """Return checkout cleanliness and a deterministic digest of porcelain status."""
    try:
        status = run_git_command(
            ["--no-pager", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=root,
        ).stdout
    except ExecutableNotFoundError, subprocess.CalledProcessError:
        return (None, hashlib.sha256(b"unavailable").hexdigest())
    return (not status.strip(), hashlib.sha256(status.encode()).hexdigest())


def _source_state_digest(root: Path) -> tuple[str, bool]:
    """Hash the measured library source so dirty runs remain identifiable."""
    source_dir = root / "src"
    files = (
        sorted(
            (path for path in source_dir.rglob("*") if path.is_file()),
            key=lambda path: path.relative_to(root).as_posix(),
        )
        if source_dir.is_dir()
        else []
    )
    digest = hashlib.sha256()
    if not files:
        digest.update(b"MISSING:src/")
        return (digest.hexdigest(), True)
    for path in files:
        relative = path.relative_to(root).as_posix().encode()
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return (digest.hexdigest(), False)


def _provenance_harness_files(root: Path) -> tuple[list[Path], list[str]]:
    """Return stable harness files and explicitly list any missing members."""
    files: list[Path] = []
    missing: list[str] = []
    for relative in _PROVENANCE_HARNESS_FILES:
        path = root / relative
        if path.is_file():
            files.append(path)
        else:
            missing.append(relative)
    benches = root / "benches"
    if benches.is_dir():
        files.extend(path for path in benches.rglob("*") if path.is_file())
    else:
        missing.append("benches/")
    return (sorted(set(files), key=lambda path: path.relative_to(root).as_posix()), sorted(missing))


def _provenance_harness_digest(root: Path) -> tuple[str, list[str]]:
    """Hash the benchmark harness while making missing inputs visible."""
    files, missing = _provenance_harness_files(root)
    digest = hashlib.sha256()
    for path in files:
        relative = path.relative_to(root).as_posix().encode()
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    for relative in missing:
        marker = f"MISSING:{relative}".encode()
        digest.update(len(marker).to_bytes(8, "big"))
        digest.update(marker)
    return (digest.hexdigest(), missing)


def _rustc_version(root: Path) -> str:
    """Return the rustc version used by the publication workflow."""
    try:
        value = run_safe_command(
            "rustc",
            ["--version"],
            cwd=root,
            timeout=60,
        ).stdout.strip()
    except ExecutableNotFoundError, subprocess.CalledProcessError:
        return "unavailable"
    return value or "unavailable"


def _capture_provenance(
    root: Path,
    *,
    args: PlotCliArgs,
    dims: list[int],
    measurement_recorded: bool,
) -> dict[str, object]:
    """Capture deterministic provenance for CSV/SVG and README publication."""
    harness_sha256, missing_harness_files = _provenance_harness_digest(root)
    cargo_lock = root / "Cargo.lock"
    cargo_lock_sha256 = hashlib.sha256(cargo_lock.read_bytes()).hexdigest() if cargo_lock.is_file() else "unavailable"
    cpu = platform.processor().strip() or platform.machine().strip() or "unavailable"
    os_description = " ".join(part for part in (platform.system(), platform.release(), platform.machine()) if part).strip() or "unavailable"
    git_clean, git_status_sha256 = _git_status_metadata(root)
    source_state_sha256, source_missing = _source_state_digest(root)
    environment: dict[str, object] = {
        "cargo_lock_sha256": cargo_lock_sha256,
        "commit": _git_value(root, ["--no-pager", "rev-parse", "HEAD"]),
        "cpu": cpu,
        "git_clean": git_clean,
        "git_status_sha256": git_status_sha256,
        "harness_sha256": harness_sha256,
        "missing_harness_files": missing_harness_files,
        "os": os_description,
        "rustc": _rustc_version(root),
        "source_missing": source_missing,
        "source_state_sha256": source_state_sha256,
    }
    measurement: dict[str, object]
    if measurement_recorded:
        measurement = {"status": "recorded", **environment}
    else:
        measurement = {
            "reason": "the exploratory renderer did not run the benchmark command that produced these Criterion samples",
            "status": "unavailable",
        }
    criterion_version = _read_cargo_dependency_versions(root / "Cargo.toml", {"criterion"}).get("criterion", "unavailable")
    return {
        "artifact": "README vs_linalg dimension plot" if args.update_readme else "exploratory vs_linalg dimension plot",
        "criterion": {
            "benchmark_command": list(_PUBLICATION_BENCHMARK) if measurement_recorded else "unavailable",
            "criterion_dependency": criterion_version,
            "dimensions": dims,
            "log_y": args.log_y,
            "metric": args.metric,
            "sample": args.sample,
            "statistic": args.stat,
        },
        "measurement": measurement,
        "publication": {
            **environment,
            "correctness_gate": "passed" if measurement_recorded else "not-run-exploratory",
        },
        "schema": 1,
    }


def _write_provenance(path: Path, provenance: dict[str, object]) -> None:
    """Write stable, sorted JSON provenance beside generated benchmark assets."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _validate_readme_target(root: Path, args: PlotCliArgs) -> int:  # noqa: C901, PLR0911
    """Validate publication-only CLI invariants and README markers before timing."""
    if not args.update_readme:
        return 0
    if args.allow_partial:
        print("--allow-partial is exploratory-only and cannot be combined with --update-readme", file=sys.stderr)
        return 2
    if args.sample != "new":
        print("README publication requires --sample new so the gated timing run is the data being published", file=sys.stderr)
        return 2
    if args.no_plot:
        print("README publication requires SVG rendering; --no-plot is exploratory-only", file=sys.stderr)
        return 2
    criterion_dir = _resolve_under_root(root, args.criterion_dir).resolve()
    expected_criterion_dir = (root / "target" / "criterion").resolve()
    if criterion_dir != expected_criterion_dir:
        print(
            f"README publication requires Criterion output at {expected_criterion_dir}; got {criterion_dir}",
            file=sys.stderr,
        )
        return 2
    readme_path = _resolve_under_root(root, args.readme)
    canonical_readme = (root / "README.md").resolve()
    if readme_path.resolve() == canonical_readme:
        expected_svg, expected_csv = _resolve_output_paths(root, args.metric, args.stat, None, None)
        selected_svg, selected_csv = _resolve_output_paths(root, args.metric, args.stat, args.out, args.csv)
        if selected_svg.resolve() != expected_svg.resolve() or selected_csv.resolve() != expected_csv.resolve():
            print(
                f"README publication requires the canonical CSV/SVG destinations referenced by README.md; expected {expected_csv} and {expected_svg}",
                file=sys.stderr,
            )
            return 2
    marker_begin, marker_end = _readme_table_markers(args.metric, args.stat, args.sample)
    try:
        lines = readme_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    begin_count = sum(line.strip() == marker_begin for line in lines)
    end_count = sum(line.strip() == marker_end for line in lines)
    if begin_count != 1 or end_count != 1:
        print(f"README markers not found or not unique (begin={begin_count}, end={end_count}).", file=sys.stderr)
        return 2
    begin_idx = next(index for index, line in enumerate(lines) if line.strip() == marker_begin)
    end_idx = next(index for index, line in enumerate(lines) if line.strip() == marker_end)
    if begin_idx >= end_idx:
        print("README markers are out of order.", file=sys.stderr)
        return 2
    return 0


def _replace_staged_files(pairs: list[tuple[Path, Path]], backup_dir: Path) -> None:
    """Replace a group of publication files and roll back on any failure."""
    backups: dict[Path, Path | None] = {}
    for index, (_staged, destination) in enumerate(pairs):
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.is_file():
            backup = backup_dir / f"backup-{index}"
            shutil.copy2(destination, backup)
            backups[destination] = backup
        elif destination.exists():
            msg = f"publication destination is not a regular file: {destination}"
            raise ValueError(msg)
        else:
            backups[destination] = None

    replaced: list[Path] = []
    try:
        for staged, destination in pairs:
            staged.replace(destination)
            replaced.append(destination)
    except OSError:
        for destination in reversed(replaced):
            backup = backups[destination]
            if backup is None:
                destination.unlink(missing_ok=True)
            else:
                backup.replace(destination)
        raise


def _stage_and_publish_outputs(  # noqa: PLR0913
    *,
    root: Path,
    args: PlotCliArgs,
    rows: list[Row],
    req: PlotRequest,
    provenance: dict[str, object],
    skipped: list[str],
) -> int:
    """Render every output in isolation, then replace publication files together."""
    final_provenance = req.csv_path.with_suffix(".provenance.json")
    with tempfile.TemporaryDirectory(prefix=".criterion-dim-plot-", dir=root) as tmp:
        stage_dir = Path(tmp)
        staged_csv = stage_dir / "benchmark.csv"
        staged_svg = stage_dir / "benchmark.svg"
        staged_provenance = stage_dir / "benchmark.provenance.json"
        _write_csv(staged_csv, rows)
        _write_provenance(staged_provenance, provenance)

        pairs: list[tuple[Path, Path]] = [
            (staged_csv, req.csv_path),
            (staged_provenance, final_provenance),
        ]
        if not args.no_plot:
            staged_request = PlotRequest(
                csv_path=staged_csv,
                out_svg=staged_svg,
                title=req.title,
                stat=req.stat,
                dims=req.dims,
                la_label=req.la_label,
                na_label=req.na_label,
                fa_label=req.fa_label,
                log_y=req.log_y,
            )
            try:
                _render_svg_with_gnuplot(staged_request)
            except (FileNotFoundError, subprocess.CalledProcessError) as exc:
                print(str(exc), file=sys.stderr)
                print("No benchmark publication files were changed.", file=sys.stderr)
                return 1
            pairs.append((staged_svg, req.out_svg))

        if args.update_readme:
            readme_path = _resolve_under_root(root, args.readme)
            staged_readme = stage_dir / "README.md"
            shutil.copy2(readme_path, staged_readme)
            marker_begin, marker_end = _readme_table_markers(args.metric, args.stat, args.sample)
            try:
                _update_readme_table(staged_readme, marker_begin, marker_end, _markdown_table(rows, args.stat))
            except (OSError, ValueError) as exc:
                print(str(exc), file=sys.stderr)
                print("No benchmark publication files were changed.", file=sys.stderr)
                return 2
            pairs.append((staged_readme, readme_path))

        try:
            _replace_staged_files(pairs, stage_dir)
        except (OSError, ValueError) as exc:
            print(f"could not publish benchmark artifacts atomically: {exc}", file=sys.stderr)
            return 2

    if skipped:
        print("Warning: some dimension groups were skipped:")
        for item in skipped:
            print(f"  - {item}")
    print(f"Wrote CSV: {req.csv_path}")
    if not args.no_plot:
        print(f"Wrote SVG: {req.out_svg}")
    print(f"Wrote provenance: {final_provenance}")
    if args.update_readme:
        print(f"Updated README table: {_resolve_under_root(root, args.readme)}")
    return 0


def _maybe_update_readme(root: Path, args: _ReadmeArgs, rows: list[Row]) -> int:
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


def _maybe_render_plot(args: _RenderArgs, req: PlotRequest, skipped: list[str]) -> int:
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


def main(argv: list[str] | None = None) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0915
    """Generate benchmark CSV and optional SVG or README output."""
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    root = _repo_root()

    rc = _validate_readme_target(root, args)
    if rc != 0:
        return rc
    if args.update_readme:
        try:
            _run_publication_benchmarks(root)
        except (FileNotFoundError, RuntimeError) as exc:
            print(str(exc), file=sys.stderr)
            return 2

    versions = _detect_versions(root)
    _print_versions(versions)

    la_label = _format_legend_label("la-stack", versions.get("la-stack", "unknown"))
    na_label = _format_legend_label("nalgebra", versions.get("nalgebra", "unknown"))
    fa_label = _format_legend_label("faer", versions.get("faer", "unknown"))

    criterion_dir = _resolve_under_root(root, args.criterion_dir)

    discovered_dims = _discover_dims(criterion_dir) if criterion_dir.exists() else []
    dims = discovered_dims if args.allow_partial else list(CANONICAL_DIMS)
    if not args.allow_partial and not discovered_dims:
        dims = []
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

    if not args.allow_partial and skipped:
        print(
            "Canonical benchmark coverage is incomplete; no CSV, SVG, provenance, or README file was written.",
            file=sys.stderr,
        )
        print("Required dimensions: " + ", ".join(f"D={dim}" for dim in CANONICAL_DIMS), file=sys.stderr)
        print("Coverage gaps:", *skipped, sep="\n  - ", file=sys.stderr)
        return 2

    dims_present = [row.dim for row in rows]
    provenance = _capture_provenance(
        root,
        args=args,
        dims=dims_present,
        measurement_recorded=args.update_readme,
    )
    publication = provenance.get("publication")
    if not isinstance(publication, dict):
        msg = "publication provenance invariant violated"
        raise TypeError(msg)
    missing_harness_files = publication.get("missing_harness_files")
    if not isinstance(missing_harness_files, list) or not all(isinstance(path, str) for path in missing_harness_files):
        msg = "publication missing_harness_files invariant violated"
        raise AssertionError(msg)
    provenance_gaps = list(cast("list[str]", missing_harness_files))
    if publication.get("source_missing") is True:
        provenance_gaps.append("src/")
    provenance_gaps.extend(field for field in ("cargo_lock_sha256", "commit", "cpu", "os", "rustc") if publication.get(field) == "unavailable")
    if publication.get("git_clean") is None:
        provenance_gaps.append("git status")
    if args.update_readme and provenance_gaps:
        print(
            "Publication provenance is incomplete; no CSV, SVG, provenance, or README file was written because required fields are unavailable: "
            + ", ".join(provenance_gaps),
            file=sys.stderr,
        )
        return 2

    title = f"{metric.title}: {args.stat} time vs dimension"
    req = PlotRequest(
        csv_path=out_csv,
        out_svg=out_svg,
        title=title,
        stat=args.stat,
        dims=tuple(dims_present),
        la_label=la_label,
        na_label=na_label,
        fa_label=fa_label,
        log_y=args.log_y,
    )

    return _stage_and_publish_outputs(
        root=root,
        args=args,
        rows=rows,
        req=req,
        provenance=provenance,
        skipped=skipped,
    )


if __name__ == "__main__":
    raise SystemExit(main())
