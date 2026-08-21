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
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Protocol, TypeGuard, cast

from benchmark_contract import benchmark_contract_digest
from performance_artifacts import ArtifactPaths, PerformanceBundle, TimingEstimate, ensure_distinct_paths, load_bundle
from subprocess_utils import ExecutableNotFoundError, cpu_description, find_project_root, run_git_command, run_safe_command


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
    performance_csv: str
    out: str | None
    csv: str | None
    log_y: bool
    no_plot: bool
    update_readme: bool
    readme: str
    allow_partial: bool


@dataclass(frozen=True, slots=True)
class _ValidatedPerformance:
    """One retained bundle bound to the checkout used for publication."""

    bundle: PerformanceBundle
    validation: dict[str, object]
    paths: ArtifactPaths


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
            _require_positive_finite_time(value, field)
        _require_confidence_interval(self.la_lo, self.la_hi, "la_stack row")
        _require_confidence_interval(self.na_lo, self.na_hi, "nalgebra row")
        _require_confidence_interval(self.fa_lo, self.fa_hi, "faer row")


class ReadmeMarkerError(ValueError):
    """Base error for invalid README BENCH_TABLE markers."""


class MarkerNotFoundError(ReadmeMarkerError):
    """Raised when README markers are missing or not unique."""


class MarkerOrderError(ReadmeMarkerError):
    """Raised when README markers are out of order."""


class ReadmeBenchmarkLinkError(ValueError):
    """Raised when canonical README benchmark artifact links are incomplete."""


class PublicationRollbackError(RuntimeError):
    """Raised when artifact publication fails and rollback is incomplete."""


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
_DEFAULT_PERFORMANCE_CSV: Final[str] = "target/bench-reports/performance.csv"
_RELEASE_VERSION_PATTERN: Final[str] = r"[0-9]+\.[0-9]+\.[0-9]+(?:-[0-9A-Za-z.-]+)?(?:\+[0-9A-Za-z.-]+)?"
_PROVENANCE_HARNESS_FILES: Final[tuple[str, ...]] = (
    ".config/nextest.toml",
    "Cargo.toml",
    "Cargo.lock",
    "rust-toolchain.toml",
    "justfile",
    "tests/exact_bench_config.rs",
    "tests/vs_linalg_inputs.rs",
)


def _repo_root() -> Path:
    """Resolve the checkout from the caller's working tree, including wheel installs."""
    return find_project_root()


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
    if "confidence_interval" not in stat_obj:
        msg = f"field 'confidence_interval' for stat '{stat}' not found in {estimates_json}"
        raise KeyError(msg)
    ci = stat_obj["confidence_interval"]
    if not _is_parsed_object(ci):
        msg = f"field 'confidence_interval' for stat '{stat}' in {estimates_json} is not an object"
        raise TypeError(msg)

    lo = _read_numeric_field(ci, "lower_bound", estimates_json, stat)
    hi = _read_numeric_field(ci, "upper_bound", estimates_json, stat)
    _require_confidence_interval(lo, hi, f"{stat}.confidence_interval in {estimates_json}")
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
) -> float:
    if field not in obj:
        msg = f"field '{field}' for stat '{stat}' not found in {estimates_json}"
        raise KeyError(msg)

    value = obj[field]
    if isinstance(value, bool) or not isinstance(value, int | float | str):
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise TypeError(msg)

    try:
        parsed = float(value)
    except (OverflowError, ValueError) as err:
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise ValueError(msg) from err
    return _require_positive_finite_time(parsed, f"{stat}.{field} in {estimates_json}")


def _require_positive_finite_time(value: float, context: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        msg = f"{context} must be finite and positive: {value!r}"
        raise ValueError(msg)
    return value


def _require_confidence_interval(lo: float, hi: float, context: str) -> None:
    if lo > hi:
        msg = f"{context} lower bound must be <= upper bound: {lo!r} > {hi!r}"
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
        f"| D | la-stack {stat} (ns) | nalgebra {stat} (ns) | faer {stat} (ns) | reduction vs nalgebra (point est.) | reduction vs faer (point est.) |",
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


def _replace_readme_benchmark_asset_versions(text: str, *, metric: str, stat: str, version: str) -> str:
    """Update the selected README benchmark links only after its artifacts exist."""
    if re.fullmatch(_RELEASE_VERSION_PATTERN, version) is None:
        msg = f"Cargo package version is not valid for a release-pinned README link: {version!r}"
        raise ReadmeBenchmarkLinkError(msg)

    asset_stem = f"/docs/assets/bench/vs_linalg_{metric}_{stat}"
    expected_assets = sorted((f"{asset_stem}.csv", f"{asset_stem}.provenance.json", f"{asset_stem}.svg"))
    pattern = re.compile(
        r"(?P<prefix>https://(?:github\.com/acgetchell/la-stack/(?:blob|raw|tree)/|raw\.githubusercontent\.com/acgetchell/la-stack/))"
        rf"v(?P<version>{_RELEASE_VERSION_PATTERN})"
        rf"(?P<asset>{re.escape(asset_stem)}(?:\.csv|\.provenance\.json|\.svg))"
    )
    matches = list(pattern.finditer(text))
    actual_assets = sorted(match.group("asset") for match in matches)
    if actual_assets != expected_assets:
        msg = f"README must contain exactly one tag-pinned link for each published benchmark asset; expected {expected_assets}, found {actual_assets}"
        raise ReadmeBenchmarkLinkError(msg)

    return pattern.sub(lambda match: f"{match.group('prefix')}v{version}{match.group('asset')}", text)


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
        help="Criterion output directory for exploratory rendering (default: target/criterion).",
    )
    parser.add_argument(
        "--performance-csv",
        default=_DEFAULT_PERFORMANCE_CSV,
        help=(
            "Retained performance-release CSV used by --update-readme; the provenance JSON must be adjacent (default: target/bench-reports/performance.csv)."
        ),
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
        performance_csv=_required_str_attr(args, "performance_csv"),
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


def _resolve_performance_paths(root: Path, performance_csv: str) -> ArtifactPaths:
    """Resolve a retained performance CSV and its adjacent provenance JSON."""
    csv = _resolve_under_root(root, performance_csv)
    return ArtifactPaths(csv=csv, provenance=csv.with_suffix(".provenance.json"))


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


def _git_value(root: Path, args: list[str]) -> str:
    """Return deterministic Git provenance or an explicit unavailable label."""
    try:
        value = run_git_command(args, cwd=root).stdout.strip()
    except (
        ExecutableNotFoundError,
        OSError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
        return "unavailable"
    return value or "unavailable"


def _git_status_metadata(root: Path) -> tuple[bool | None, str]:
    """Return checkout cleanliness and a deterministic digest of porcelain status."""
    try:
        status = run_git_command(
            ["--no-pager", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=root,
        ).stdout
    except (
        ExecutableNotFoundError,
        OSError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
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
    except (
        ExecutableNotFoundError,
        OSError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
        return "unavailable"
    return value or "unavailable"


def _capture_environment(root: Path) -> dict[str, object]:
    """Capture the current source, harness, toolchain, host, and Git state."""
    harness_sha256, missing_harness_files = _provenance_harness_digest(root)
    cargo_lock = root / "Cargo.lock"
    cargo_lock_sha256 = hashlib.sha256(cargo_lock.read_bytes()).hexdigest() if cargo_lock.is_file() else "unavailable"
    cpu = cpu_description()
    os_description = " ".join(part for part in (platform.system(), platform.release(), platform.machine()) if part).strip() or "unavailable"
    git_clean, git_status_sha256 = _git_status_metadata(root)
    source_state_sha256, source_missing = _source_state_digest(root)
    return {
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


def _capture_exploratory_provenance(root: Path, *, args: PlotCliArgs, dims: list[int]) -> dict[str, object]:
    """Capture provenance for rendering unverified local Criterion samples."""
    environment = _capture_environment(root)
    criterion_version = _read_cargo_dependency_versions(root / "Cargo.toml", {"criterion"}).get("criterion", "unavailable")
    return {
        "artifact": "exploratory vs_linalg dimension plot",
        "criterion": {
            "benchmark_command": "unavailable",
            "criterion_dependency": criterion_version,
            "dimensions": dims,
            "log_y": args.log_y,
            "metric": args.metric,
            "sample": args.sample,
            "statistic": args.stat,
        },
        "measurement": {
            "reason": "the exploratory renderer did not run the benchmark command that produced these Criterion samples",
            "status": "unavailable",
        },
        "publication": {
            **environment,
            "correctness_gate": "not-run-exploratory",
        },
        "schema": 1,
    }


def _required_mapping(value: object, context: str) -> Mapping[str, object]:
    """Parse one already-validated frozen provenance object."""
    if not isinstance(value, Mapping) or not all(isinstance(key, str) for key in value):
        msg = f"retained performance provenance requires a {context} object"
        raise TypeError(msg)
    return cast("Mapping[str, object]", value)


def _json_compatible(value: object) -> object:
    """Convert frozen retained-provenance containers into JSON containers."""
    if isinstance(value, Mapping):
        return {str(key): _json_compatible(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_json_compatible(item) for item in value]
    return value


def _validate_performance_selection(root: Path, args: PlotCliArgs, bundle: PerformanceBundle) -> Mapping[str, object]:
    """Validate the retained report selection and return recorded measurement metadata."""
    context = bundle.context
    if context.suite not in {"all", "vs_linalg"}:
        msg = f"retained performance data does not include vs_linalg: suite={context.suite!r}"
        raise ValueError(msg)
    if context.scope != "release-signal":
        msg = f"README publication requires release-signal performance data, got {context.scope!r}"
        raise ValueError(msg)
    if context.statistic != args.stat:
        msg = f"README statistic {args.stat!r} does not match retained performance statistic {context.statistic!r}"
        raise ValueError(msg)
    if context.release.current == context.release.baseline:
        msg = "README publication requires a distinct release pair retained by `just performance-release`"
        raise ValueError(msg)

    package_version = _read_cargo_package_version(root / "Cargo.toml")
    if package_version is None:
        msg = f"{root / 'Cargo.toml'} has no string package version"
        raise ValueError(msg)
    expected_release = f"v{package_version.removeprefix('v')}"
    if context.release.current != expected_release:
        msg = (
            f"retained performance data is for {context.release.current}, but Cargo.toml is {expected_release}; "
            "run `just performance-release` for the current release"
        )
        raise ValueError(msg)

    benchmark_provenance = context.benchmark_provenance
    if benchmark_provenance.get("mode") != "shared-current-harness":
        msg = "README publication requires locally recorded shared-current-harness performance data"
        raise ValueError(msg)
    measurement = _required_mapping(benchmark_provenance.get("measurement"), "measurement")
    if measurement.get("status") != "recorded":
        msg = "README publication requires recorded performance measurements"
        raise ValueError(msg)
    return measurement


def _validate_current_measurement(root: Path, measurement: Mapping[str, object]) -> dict[str, object]:
    """Require retained measurements to match stable current checkout inputs."""
    cargo_lock = root / "Cargo.lock"
    current_source_state, source_missing = _source_state_digest(root)
    if source_missing:
        msg = "current measured source directory is missing: src/"
        raise ValueError(msg)
    current = {
        "cargo_lock_sha256": hashlib.sha256(cargo_lock.read_bytes()).hexdigest() if cargo_lock.is_file() else "unavailable",
        "commit": _git_value(root, ["--no-pager", "rev-parse", "HEAD"]),
        "source_state_sha256": current_source_state,
    }
    for retained_field, current_field in (
        ("cargo_lock_sha256", "cargo_lock_sha256"),
        ("current_commit", "commit"),
        ("current_source_state_sha256", "source_state_sha256"),
    ):
        retained = measurement.get(retained_field)
        current_value = current.get(current_field)
        if retained != current_value:
            msg = (
                f"retained performance {retained_field} {retained!r} does not match the current checkout {current_field} {current_value!r}; "
                "run `just performance-release` again"
            )
            raise ValueError(msg)

    retained_contract = measurement.get("benchmark_contract_sha256")
    contract_status = "legacy-retained-artifact"
    if retained_contract is not None:
        current_contract = benchmark_contract_digest(root)
        if retained_contract != current_contract:
            msg = (
                f"retained performance benchmark_contract_sha256 {retained_contract!r} does not match "
                f"the current checkout benchmark contract {current_contract!r}; run `just performance-release` again"
            )
            raise ValueError(msg)
        contract_status = "matched"
    return {
        **current,
        "benchmark_contract": contract_status,
    }


def _validate_performance_bundle(root: Path, args: PlotCliArgs, bundle: PerformanceBundle) -> dict[str, object]:
    """Bind retained release measurements to the current release checkout."""
    measurement = _validate_performance_selection(root, args, bundle)
    return _validate_current_measurement(root, measurement)


def _row_from_estimates(
    *,
    dim: int,
    la_stack: TimingEstimate,
    nalgebra: TimingEstimate,
    faer: TimingEstimate,
) -> Row:
    return Row(
        dim=dim,
        la_time=la_stack.median_ns,
        la_lo=la_stack.ci_lower_ns,
        la_hi=la_stack.ci_upper_ns,
        na_time=nalgebra.median_ns,
        na_lo=nalgebra.ci_lower_ns,
        na_hi=nalgebra.ci_upper_ns,
        fa_time=faer.median_ns,
        fa_lo=faer.ci_lower_ns,
        fa_hi=faer.ci_upper_ns,
    )


def _collect_performance_rows(bundle: PerformanceBundle, metric: Metric) -> list[Row]:
    """Extract current la-stack and retained same-harness peer timings."""
    rows_by_key = {(row.group, row.benchmark): row for row in bundle.rows if row.suite == "vs_linalg"}
    rows: list[Row] = []
    for dim in CANONICAL_DIMS:
        key = (f"d{dim}", metric.la_bench)
        performance_row = rows_by_key.get(key)
        if performance_row is None:
            msg = f"retained performance data is missing {key[0]}/{key[1]}"
            raise ValueError(msg)
        missing = [
            name
            for name, estimate in (
                ("current la-stack", performance_row.current),
                ("nalgebra peer", performance_row.baseline_nalgebra),
                ("faer peer", performance_row.baseline_faer),
            )
            if estimate is None
        ]
        if missing:
            msg = f"retained performance data for {key[0]}/{key[1]} is missing {', '.join(missing)} timing"
            raise ValueError(msg)
        if performance_row.current is None or performance_row.baseline_nalgebra is None or performance_row.baseline_faer is None:
            msg = "retained performance timing presence invariant violated"
            raise AssertionError(msg)
        rows.append(
            _row_from_estimates(
                dim=dim,
                la_stack=performance_row.current,
                nalgebra=performance_row.baseline_nalgebra,
                faer=performance_row.baseline_faer,
            )
        )
    return rows


def _display_path(root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _capture_performance_provenance(
    root: Path,
    *,
    args: PlotCliArgs,
    dims: list[int],
    retained: _ValidatedPerformance,
) -> dict[str, object]:
    """Derive README artifact provenance from a validated performance bundle."""
    bundle = retained.bundle
    paths = retained.paths
    benchmark_provenance = bundle.context.benchmark_provenance
    criterion = _required_mapping(benchmark_provenance.get("criterion"), "criterion")
    measurement = _required_mapping(benchmark_provenance.get("measurement"), "measurement")
    return {
        "artifact": "README vs_linalg dimension plot",
        "criterion": {
            "baseline_command": _json_compatible(criterion.get("baseline_command")),
            "criterion_dependency": criterion.get("criterion_version"),
            "current_command": _json_compatible(criterion.get("current_command")),
            "dimensions": dims,
            "log_y": args.log_y,
            "metric": args.metric,
            "sample": args.sample,
            "statistic": args.stat,
        },
        "measurement": {
            **cast("dict[str, object]", _json_compatible(measurement)),
            "la_stack_sample": "current",
            "peer_release_context": bundle.context.release.baseline,
            "peer_sample": "baseline phase under shared current harness",
            "source": "retained performance-release artifact",
        },
        "performance_artifact": {
            "csv": _display_path(root, paths.csv),
            "csv_sha256": hashlib.sha256(paths.csv.read_bytes()).hexdigest(),
            "provenance": _display_path(root, paths.provenance),
            "provenance_sha256": hashlib.sha256(paths.provenance.read_bytes()).hexdigest(),
            "release": {
                "baseline": bundle.context.release.baseline,
                "current": bundle.context.release.current,
            },
        },
        "publication": {
            **retained.validation,
            "correctness_gate": "validated-performance-release-artifact",
        },
        "schema": 2,
    }


def _write_provenance(path: Path, provenance: dict[str, object]) -> None:
    """Write stable, sorted JSON provenance beside generated benchmark assets."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _validate_readme_target(root: Path, args: PlotCliArgs) -> int:  # noqa: C901, PLR0911, PLR0912
    """Validate publication-only CLI invariants and README markers before rendering."""
    if not args.update_readme:
        return 0
    if args.allow_partial:
        print("--allow-partial is exploratory-only and cannot be combined with --update-readme", file=sys.stderr)
        return 2
    if args.sample != "new":
        print("README publication requires --sample new to identify the retained current-release measurements", file=sys.stderr)
        return 2
    if args.stat != "median":
        print("README publication requires --stat median to match performance-release artifacts", file=sys.stderr)
        return 2
    if args.no_plot:
        print("README publication requires SVG rendering; --no-plot is exploratory-only", file=sys.stderr)
        return 2
    readme_path = _resolve_under_root(root, args.readme)
    canonical_readme = (root / "README.md").resolve()
    try:
        readme_text = readme_path.read_text(encoding="utf-8")
    except OSError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if readme_path.resolve() == canonical_readme:
        expected_svg, expected_csv = _resolve_output_paths(root, args.metric, args.stat, None, None)
        selected_svg, selected_csv = _resolve_output_paths(root, args.metric, args.stat, args.out, args.csv)
        if selected_svg.resolve() != expected_svg.resolve() or selected_csv.resolve() != expected_csv.resolve():
            print(
                f"README publication requires the canonical CSV/SVG destinations referenced by README.md; expected {expected_csv} and {expected_svg}",
                file=sys.stderr,
            )
            return 2
        expected_performance = _resolve_performance_paths(root, _DEFAULT_PERFORMANCE_CSV)
        selected_performance = _resolve_performance_paths(root, args.performance_csv)
        if selected_performance.csv.resolve() != expected_performance.csv.resolve():
            print(
                f"README publication requires the canonical performance-release input {expected_performance.csv}; got {selected_performance.csv}",
                file=sys.stderr,
            )
            return 2
        try:
            package_version = _read_cargo_package_version(root / "Cargo.toml")
            if package_version is None:
                msg = f"{root / 'Cargo.toml'} has no string package version"
                raise ReadmeBenchmarkLinkError(msg)
            _replace_readme_benchmark_asset_versions(
                readme_text,
                metric=args.metric,
                stat=args.stat,
                version=package_version,
            )
        except (OSError, TypeError, ValueError, tomllib.TOMLDecodeError) as exc:
            print(str(exc), file=sys.stderr)
            return 2
    marker_begin, marker_end = _readme_table_markers(args.metric, args.stat, args.sample)
    lines = readme_text.splitlines()
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


def _validate_publication_paths(root: Path, args: PlotCliArgs, *, out_svg: Path, out_csv: Path) -> int:
    """Reject input/output aliases before publication can begin."""
    paths = {
        "CSV output": out_csv,
        "provenance output": out_csv.with_suffix(".provenance.json"),
    }
    if not args.no_plot:
        paths["SVG output"] = out_svg
    if args.update_readme:
        paths["README output"] = _resolve_under_root(root, args.readme)
        performance = _resolve_performance_paths(root, args.performance_csv)
        paths["performance input CSV"] = performance.csv
        paths["performance input provenance"] = performance.provenance
    try:
        ensure_distinct_paths(paths)
    except (OSError, ValueError) as exc:
        print(f"Invalid benchmark publication paths: {exc}", file=sys.stderr)
        return 2
    return 0


def _changed_staged_files(pairs: list[tuple[Path, Path]]) -> list[tuple[Path, Path]]:
    """Return staged files whose destination bytes differ or do not exist."""
    changed_pairs: list[tuple[Path, Path]] = []
    for staged, destination in pairs:
        if destination.is_file() and staged.read_bytes() == destination.read_bytes():
            continue
        changed_pairs.append((staged, destination))
    return changed_pairs


def _replace_staged_files(pairs: list[tuple[Path, Path]], backup_dir: Path) -> None:
    """Replace a group of publication files and roll back on any failure."""
    changed_pairs = _changed_staged_files(pairs)
    backups: dict[Path, Path | None] = {}
    for index, (_staged, destination) in enumerate(changed_pairs):
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
        for staged, destination in changed_pairs:
            staged.replace(destination)
            replaced.append(destination)
    except OSError as primary:
        rollback_errors: list[str] = []
        for destination in reversed(replaced):
            backup = backups[destination]
            try:
                if backup is None:
                    destination.unlink(missing_ok=True)
                else:
                    backup.replace(destination)
            except OSError as rollback:
                rollback_errors.append(f"could not restore {destination}: {rollback}")
        if rollback_errors:
            msg = f"artifact replacement failed ({primary}); rollback failed: {'; '.join(rollback_errors)}; backups preserved at {backup_dir}"
            raise PublicationRollbackError(msg) from primary
        raise


def _publish_staged_files(pairs: list[tuple[Path, Path]], root: Path) -> bool:
    """Publish staged files together, preserving backups after rollback failure."""
    backup_dir = Path(tempfile.mkdtemp(prefix=".criterion-dim-plot-backup-", dir=root))
    preserve_backup = False
    try:
        _replace_staged_files(pairs, backup_dir)
    except PublicationRollbackError as exc:
        preserve_backup = True
        print(f"could not publish benchmark artifacts atomically: {exc}", file=sys.stderr)
        return False
    except (OSError, ValueError) as exc:
        print(f"could not publish benchmark artifacts atomically: {exc}", file=sys.stderr)
        return False
    finally:
        if not preserve_backup:
            try:
                shutil.rmtree(backup_dir)
            except OSError as exc:
                print(f"Warning: could not remove artifact backup {backup_dir}: {exc}", file=sys.stderr)
    return True


def _update_staged_readme_publication(
    *,
    root: Path,
    args: PlotCliArgs,
    rows: list[Row],
    staged_readme: Path,
) -> None:
    """Stage the README table and canonical artifact links as one publication."""
    marker_begin, marker_end = _readme_table_markers(args.metric, args.stat, args.sample)
    _update_readme_table(staged_readme, marker_begin, marker_end, _markdown_table(rows, args.stat))

    readme_path = _resolve_under_root(root, args.readme)
    if readme_path.resolve() != (root / "README.md").resolve():
        return

    package_version = _read_cargo_package_version(root / "Cargo.toml")
    if package_version is None:
        msg = f"{root / 'Cargo.toml'} has no string package version"
        raise ReadmeBenchmarkLinkError(msg)
    updated_readme = _replace_readme_benchmark_asset_versions(
        staged_readme.read_text(encoding="utf-8"),
        metric=args.metric,
        stat=args.stat,
        version=package_version,
    )
    staged_readme.write_text(updated_readme, encoding="utf-8")


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
            try:
                _update_staged_readme_publication(root=root, args=args, rows=rows, staged_readme=staged_readme)
            except (OSError, ValueError) as exc:
                print(str(exc), file=sys.stderr)
                print("No benchmark publication files were changed.", file=sys.stderr)
                return 2
            pairs.append((staged_readme, readme_path))

        if not _publish_staged_files(pairs, root):
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
        print(f"Updated README benchmark publication: {_resolve_under_root(root, args.readme)}")
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
    out_svg, out_csv = _resolve_output_paths(root, args.metric, args.stat, args.out, args.csv)
    rc = _validate_publication_paths(root, args, out_svg=out_svg, out_csv=out_csv)
    if rc != 0:
        return rc
    versions = _detect_versions(root)
    _print_versions(versions)

    la_label = _format_legend_label("la-stack", versions.get("la-stack", "unknown"))
    na_label = _format_legend_label("nalgebra", versions.get("nalgebra", "unknown"))
    fa_label = _format_legend_label("faer", versions.get("faer", "unknown"))

    metric = METRICS[args.metric]
    if args.update_readme:
        performance_paths = _resolve_performance_paths(root, args.performance_csv)
        try:
            performance_bundle = load_bundle(performance_paths)
            validation = _validate_performance_bundle(root, args, performance_bundle)
            rows = _collect_performance_rows(performance_bundle, metric)
            retained = _ValidatedPerformance(bundle=performance_bundle, validation=validation, paths=performance_paths)
            provenance = _capture_performance_provenance(
                root,
                args=args,
                dims=[row.dim for row in rows],
                retained=retained,
            )
        except (FileNotFoundError, OSError, TypeError, ValueError, AssertionError, tomllib.TOMLDecodeError) as exc:
            print(f"Invalid retained performance-release data: {exc}", file=sys.stderr)
            return 2
        skipped: list[str] = []
    else:
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
        try:
            rows, skipped = _collect_rows(criterion_dir, dims, metric, args.stat, args.sample)
        except (OSError, KeyError, TypeError, ValueError) as exc:
            print(f"Invalid Criterion estimate data: {exc}", file=sys.stderr)
            return 2
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
        provenance = _capture_exploratory_provenance(root, args=args, dims=[row.dim for row in rows])

    dims_present = [row.dim for row in rows]
    if not args.update_readme:
        publication = provenance.get("publication")
        if not isinstance(publication, dict):
            msg = "publication provenance invariant violated"
            raise TypeError(msg)
        missing_harness_files = publication.get("missing_harness_files")
        if not isinstance(missing_harness_files, list) or not all(isinstance(path, str) for path in missing_harness_files):
            msg = "publication missing_harness_files invariant violated"
            raise AssertionError(msg)

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
