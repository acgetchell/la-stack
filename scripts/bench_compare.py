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
import math
import re
import subprocess
import sys
import tomllib
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal, Protocol, cast

from criterion_dim_plot import METRICS
from subprocess_utils import ExecutableNotFoundError, run_git_command

# ---------------------------------------------------------------------------
# Benchmark group / bench discovery
# ---------------------------------------------------------------------------

# Groups and the benchmarks within each group that we track.
#
# Mirrors the structure of `benches/exact.rs`: general-case per-dimension
# groups (`exact_d{2..5}`), fixed-seed full-corpus groups, plus
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
_RANDOM_CORPUS_BENCHES: list[str] = _EXTREME_BENCHES.copy()

_EXACT_DIMENSION_BENCHES: list[str] = [
    "det",
    "det_exact",
    "det_exact_f64_result",
    "det_exact_rounded_f64",
    "det_sign_exact",
    "solve_exact",
    "solve_exact_f64_result",
    "solve_exact_rounded_f64",
]
_EXACT_DIMENSION_BENCHES_WITH_DIRECT: list[str] = [
    "det",
    "det_direct",
    "det_direct_with_errbound",
    "det_errbound",
    *_EXACT_DIMENSION_BENCHES[1:],
]

EXACT_GROUPS: dict[str, list[str]] = {
    "exact_d2": _EXACT_DIMENSION_BENCHES_WITH_DIRECT,
    "exact_d3": _EXACT_DIMENSION_BENCHES_WITH_DIRECT,
    "exact_d4": _EXACT_DIMENSION_BENCHES_WITH_DIRECT,
    "exact_d5": _EXACT_DIMENSION_BENCHES,
    "exact_random_corpus_d2": _RANDOM_CORPUS_BENCHES,
    "exact_random_corpus_d3": _RANDOM_CORPUS_BENCHES,
    "exact_random_corpus_d4": _RANDOM_CORPUS_BENCHES,
    "exact_random_corpus_d5": _RANDOM_CORPUS_BENCHES,
    "exact_near_singular_3x3": _EXTREME_BENCHES,
    "exact_large_entries_3x3": _EXTREME_BENCHES,
    "exact_hilbert_4x4": _EXTREME_BENCHES,
    "exact_hilbert_5x5": _EXTREME_BENCHES,
}

EXACT_RELEASE_SIGNAL_GROUPS: frozenset[str] = frozenset(EXACT_GROUPS)

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
VS_LINALG_STANDARD_BENCH_ORDER: list[str] = [
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

VS_LINALG_D8_RELEASE_SIGNAL_BENCHES: list[str] = [
    "la_stack_lu_pivoting",
    "la_stack_lu_ill_conditioned",
    "la_stack_ldlt_ill_conditioned",
    "la_stack_det_from_lu_balanced_range",
    "la_stack_det_from_ldlt_balanced_range",
]
_V0_4_3_API_COMPATIBILITY = "la_stack_v0_4_3_api"
_V0_4_3_UNAVAILABLE_BASELINE_ROWS: frozenset[tuple[str, str]] = frozenset(
    {
        ("exact_d2", "det_direct_with_errbound"),
        ("exact_d3", "det_direct_with_errbound"),
        ("exact_d4", "det_direct_with_errbound"),
        ("d8", "la_stack_det_from_lu_balanced_range"),
        ("d8", "la_stack_det_from_ldlt_balanced_range"),
    }
)
_UNAVAILABLE_BASELINE_ROWS_BY_COMPATIBILITY: dict[str, frozenset[tuple[str, str]]] = {
    _V0_4_3_API_COMPATIBILITY: _V0_4_3_UNAVAILABLE_BASELINE_ROWS,
}
VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM: dict[int, list[str]] = {
    8: VS_LINALG_D8_RELEASE_SIGNAL_BENCHES,
}
VS_LINALG_CANONICAL_DIMS: tuple[int, ...] = (2, 3, 4, 5, 8, 16, 32, 64)
VS_LINALG_BENCH_ORDER: list[str] = [
    *VS_LINALG_STANDARD_BENCH_ORDER,
    *VS_LINALG_D8_RELEASE_SIGNAL_BENCHES,
]

VS_LINALG_LA_STACK_BENCHES: frozenset[str] = frozenset(bench for bench in VS_LINALG_STANDARD_BENCH_ORDER if bench.startswith("la_stack_"))
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

type ChangeAssessment = Literal["improvement", "regression", "inconclusive", "unknown"]
type BenchmarkSuite = Literal["all", "exact", "vs_linalg"]
type ComparisonScope = Literal["release-signal", "all-benches"]
type Statistic = Literal["mean", "median"]


@dataclass(frozen=True, slots=True)
class CriterionEstimate:
    """A Criterion point estimate and its optional confidence interval."""

    point_ns: float
    ci_lo_ns: float | None
    ci_hi_ns: float | None

    def __post_init__(self) -> None:
        """Keep every stored timing finite, positive, and interval-complete."""
        for field, value in (
            ("point_ns", self.point_ns),
            ("ci_lo_ns", self.ci_lo_ns),
            ("ci_hi_ns", self.ci_hi_ns),
        ):
            if value is not None and (not math.isfinite(value) or value <= 0):
                msg = f"{field} must be finite and positive: {value!r}"
                raise ValueError(msg)
        if (self.ci_lo_ns is None) != (self.ci_hi_ns is None):
            msg = "Criterion confidence interval must contain both bounds or neither"
            raise ValueError(msg)
        if self.ci_lo_ns is not None and self.ci_hi_ns is not None and self.ci_lo_ns > self.ci_hi_ns:
            msg = f"Criterion confidence interval lower bound exceeds upper bound: {self.ci_lo_ns} > {self.ci_hi_ns}"
            raise ValueError(msg)

    @property
    def has_confidence_interval(self) -> bool:
        """Return whether both confidence bounds were present."""
        return self.ci_lo_ns is not None and self.ci_hi_ns is not None


@dataclass(frozen=True, slots=True)
class BenchResult:
    """A single benchmark measurement (point estimate + confidence interval)."""

    suite: str
    group: str
    bench: str
    estimate: CriterionEstimate

    @property
    def point_ns(self) -> float:
        """Return the point estimate in nanoseconds."""
        return self.estimate.point_ns

    @property
    def ci_lo_ns(self) -> float | None:
        """Return the lower confidence bound, when Criterion recorded it."""
        return self.estimate.ci_lo_ns

    @property
    def ci_hi_ns(self) -> float | None:
        """Return the upper confidence bound, when Criterion recorded it."""
        return self.estimate.ci_hi_ns


@dataclass(frozen=True, slots=True)
class Comparison:
    """A comparison between baseline and current benchmark results."""

    suite: str
    group: str
    bench: str
    baseline: CriterionEstimate
    current: CriterionEstimate
    assessment: ChangeAssessment
    baseline_bench: str | None = None
    baseline_nalgebra: CriterionEstimate | None = None
    baseline_faer: CriterionEstimate | None = None

    @property
    def baseline_ns(self) -> float:
        """Return the baseline point estimate."""
        return self.baseline.point_ns

    @property
    def current_ns(self) -> float:
        """Return the current point estimate."""
        return self.current.point_ns

    @property
    def speedup(self) -> float:
        """Return baseline/current, where values above one are faster."""
        return self.baseline_ns / self.current_ns

    @property
    def pct_change(self) -> float:
        """Return signed point-estimate change, where negative is faster."""
        return ((self.current_ns - self.baseline_ns) / self.baseline_ns) * 100.0

    @property
    def baseline_nalgebra_ns(self) -> float | None:
        """Return the baseline nalgebra point estimate, when available."""
        return None if self.baseline_nalgebra is None else self.baseline_nalgebra.point_ns

    @property
    def baseline_faer_ns(self) -> float | None:
        """Return the baseline faer point estimate, when available."""
        return None if self.baseline_faer is None else self.baseline_faer.point_ns


@dataclass(frozen=True, slots=True)
class CoverageGap:
    """An expected comparison row missing one or both Criterion samples."""

    suite: str
    group: str
    bench: str
    baseline_bench: str
    missing_current: bool
    missing_baseline: bool


@dataclass(frozen=True, slots=True)
class ComparisonCollection:
    """Complete comparisons plus deterministic coverage gaps."""

    comparisons: list[Comparison]
    gaps: list[CoverageGap]


@dataclass(frozen=True, slots=True)
class ComparisonPolicy:
    """Coverage policy for one baseline comparison."""

    scope: str = "release-signal"
    baseline_api_compatibility: str | None = None


_DEFAULT_COMPARISON_POLICY = ComparisonPolicy()


@dataclass(frozen=True, slots=True)
class CriterionSelection:
    """Requested Criterion settings that provenance must describe exactly."""

    suite: BenchmarkSuite
    scope: ComparisonScope
    statistic: Statistic
    sample: str


@dataclass(frozen=True, slots=True)
class HarnessProvenance:
    """Validated benchmark measurement and correctness provenance."""

    schema: int
    mode: str
    sha256: str | None
    baseline: str
    measurement: dict[str, object] | None = None
    publication: dict[str, object] | None = None
    criterion: CriterionProvenance | None = None
    validation: dict[str, object] | None = None


@dataclass(frozen=True, slots=True)
class CriterionProvenance:
    """Validated Criterion settings and commands recorded for one report."""

    suite: BenchmarkSuite
    scope: ComparisonScope
    statistic: Statistic
    sample: str
    criterion_version: str
    baseline_command: tuple[str, ...]
    current_command: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ReportSettings:
    """Settings rendered into the benchmark report header."""

    baseline_name: str | None
    stat: Statistic
    suite: BenchmarkSuite
    scope: ComparisonScope
    harness_provenance: HarnessProvenance | None = None


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


def _read_estimate(estimates_json: Path, stat: str = "median") -> CriterionEstimate:
    """Read and validate a Criterion point estimate and confidence interval."""
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
    if "confidence_interval" not in stat_obj:
        return CriterionEstimate(point_ns=point, ci_lo_ns=None, ci_hi_ns=None)
    ci = stat_obj["confidence_interval"]
    if not isinstance(ci, dict):
        msg = f"field 'confidence_interval' for stat '{stat}' in {estimates_json} is not an object"
        raise TypeError(msg)

    lo = _read_numeric_field(ci, "lower_bound", estimates_json, stat)
    hi = _read_numeric_field(ci, "upper_bound", estimates_json, stat)
    if lo > hi:
        msg = f"invalid confidence interval for stat '{stat}' in {estimates_json}: lower_bound {lo} exceeds upper_bound {hi}"
        raise ValueError(msg)
    return CriterionEstimate(point_ns=point, ci_lo_ns=lo, ci_hi_ns=hi)


def _read_numeric_field(
    obj: dict[str, object],
    field: str,
    estimates_json: Path,
    stat: str,
) -> float:
    """Read a numeric Criterion field with file and statistic context."""
    if field not in obj:
        msg = f"field '{field}' for stat '{stat}' not found in {estimates_json}"
        raise KeyError(msg)
    value = obj[field]
    if isinstance(value, bool) or not isinstance(value, int | float | str):
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise TypeError(msg)
    try:
        result = float(value)
    except (OverflowError, ValueError) as err:
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} is not numeric: {value!r}"
        raise ValueError(msg) from err
    if not math.isfinite(result) or result <= 0:
        msg = f"field '{field}' for stat '{stat}' in {estimates_json} must be finite and positive: {value!r}"
        raise ValueError(msg)
    return result


def _read_harness_provenance(
    criterion_dir: Path,
    *,
    expected_baseline: str,
    expected: CriterionSelection,
) -> HarnessProvenance | None:
    """Read provenance tied to the exact Criterion samples being compared."""
    provenance_path = criterion_dir / ".la-stack-benchmark-harness.json"
    if not provenance_path.exists():
        return None

    try:
        data = json.loads(provenance_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as err:
        msg = f"malformed benchmark harness provenance JSON in {provenance_path}: {err}"
        raise ValueError(msg) from err
    if not isinstance(data, dict):
        msg = f"expected JSON object in {provenance_path}"
        raise TypeError(msg)

    schema = data.get("schema")
    mode = data.get("mode")
    baseline = _required_metadata_string(data, "baseline", provenance_path)
    if baseline != expected_baseline:
        msg = f"benchmark harness provenance baseline {baseline!r} does not match requested Criterion baseline {expected_baseline!r} in {provenance_path}"
        raise ValueError(msg)

    if not isinstance(schema, bool) and schema == 1:
        if mode != "shared-current-harness":
            msg = f"unsupported or missing mode in {provenance_path}: {mode!r}"
            raise ValueError(msg)
        sha256 = _required_sha256(data, "sha256", provenance_path)
        return HarnessProvenance(schema=1, mode=mode, sha256=sha256, baseline=baseline)

    if isinstance(schema, bool) or schema != 2:
        msg = f"unsupported or missing schema in {provenance_path}: expected 1 or 2, got {schema!r}"
        raise ValueError(msg)
    if mode not in {"shared-current-harness", "historical-assets"}:
        msg = f"unsupported or missing mode in {provenance_path}: {mode!r}"
        raise ValueError(msg)

    measurement = _required_metadata_object(data, "measurement", provenance_path)
    publication = _required_metadata_object(data, "publication", provenance_path)
    criterion_data = _required_metadata_object(data, "criterion", provenance_path)
    validation = _required_metadata_object(data, "validation", provenance_path)
    _validate_measurement_metadata(measurement, mode=mode, path=provenance_path)
    _validate_environment_metadata(publication, path=provenance_path, context="publication")
    criterion = _parse_criterion_metadata(
        criterion_data,
        path=provenance_path,
        expected=expected,
    )
    _validate_validation_metadata(validation, path=provenance_path)
    _validate_baseline_api_compatibility(validation, baseline=baseline, path=provenance_path)

    sha256: str | None = None
    if measurement.get("status") == "recorded":
        sha256 = _required_sha256(measurement, "harness_sha256", provenance_path)
    return HarnessProvenance(
        schema=2,
        mode=mode,
        sha256=sha256,
        baseline=baseline,
        measurement=measurement,
        publication=publication,
        criterion=criterion,
        validation=validation,
    )


def _required_metadata_object(data: dict[str, object], field: str, path: Path) -> dict[str, object]:
    """Return a required provenance object with contextual diagnostics."""
    value = data.get(field)
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        msg = f"invalid or missing {field} object in {path}"
        raise ValueError(msg)
    return cast("dict[str, object]", value)


def _required_metadata_string(data: dict[str, object], field: str, path: Path) -> str:
    """Return a required non-empty provenance string."""
    value = data.get(field)
    if not isinstance(value, str) or not value.strip():
        msg = f"invalid or missing {field} in {path}"
        raise ValueError(msg)
    return value


def _required_sha256(data: dict[str, object], field: str, path: Path) -> str:
    """Return a required lowercase SHA-256 digest."""
    value = data.get(field)
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        msg = f"invalid or missing {field} in {path}: expected 64 lowercase hexadecimal characters"
        raise ValueError(msg)
    return value


def _validate_environment_metadata(data: dict[str, object], *, path: Path, context: str) -> None:
    """Validate deterministic environment fields used to reproduce a run."""
    for field in ("cpu", "os", "rustc", "commit"):
        _required_metadata_string(data, field, path)
    _required_sha256(data, "cargo_lock_sha256", path)
    _required_sha256(data, "harness_sha256", path)
    _required_sha256(data, "source_state_sha256", path)
    if not isinstance(data.get("git_clean"), bool):
        msg = f"invalid or missing {context}.git_clean in {path}"
        raise TypeError(msg)
    gate = _required_metadata_string(data, "correctness_gate", path)
    if gate != "passed":
        msg = f"{context}.correctness_gate in {path} must be 'passed', got {gate!r}"
        raise ValueError(msg)


def _validate_measurement_metadata(data: dict[str, object], *, mode: object, path: Path) -> None:
    """Validate recorded or explicitly unavailable measurement provenance."""
    status = _required_metadata_string(data, "status", path)
    if status == "recorded":
        if mode != "shared-current-harness":
            msg = f"recorded measurement provenance in {path} requires shared-current-harness mode"
            raise ValueError(msg)
        for field in ("cpu", "os", "rustc", "current_commit", "baseline_commit"):
            _required_metadata_string(data, field, path)
        _required_sha256(data, "cargo_lock_sha256", path)
        _required_sha256(data, "harness_sha256", path)
        _required_sha256(data, "current_source_state_sha256", path)
        _required_sha256(data, "baseline_source_state_sha256", path)
        for field in ("current_git_clean", "baseline_git_clean"):
            if not isinstance(data.get(field), bool):
                msg = f"invalid or missing measurement.{field} in {path}"
                raise TypeError(msg)
        return
    if status == "unavailable":
        _required_metadata_string(data, "reason", path)
        return
    msg = f"unsupported measurement.status in {path}: {status!r}"
    raise ValueError(msg)


def _parse_criterion_metadata(
    data: dict[str, object],
    *,
    path: Path,
    expected: CriterionSelection,
) -> CriterionProvenance:
    """Parse Criterion metadata and bind it to the requested report settings."""
    suite = _required_metadata_string(data, "suite", path)
    scope = _required_metadata_string(data, "scope", path)
    statistic = _required_metadata_string(data, "statistic", path)
    sample = _required_metadata_string(data, "sample", path)
    criterion_version = _required_metadata_string(data, "criterion_version", path)

    expected_fields = {
        "suite": expected.suite,
        "scope": expected.scope,
        "statistic": expected.statistic,
        "sample": expected.sample,
    }
    observed_fields = {
        "suite": suite,
        "scope": scope,
        "statistic": statistic,
        "sample": sample,
    }
    for field, expected_value in expected_fields.items():
        observed = observed_fields[field]
        if observed != expected_value:
            msg = f"criterion.{field} {observed!r} in {path} does not match requested value {expected_value!r}"
            raise ValueError(msg)

    if suite not in SUITE_CHOICES:
        msg = f"unsupported criterion.suite in {path}: {suite!r}"
        raise ValueError(msg)
    if scope not in SCOPE_CHOICES:
        msg = f"unsupported criterion.scope in {path}: {scope!r}"
        raise ValueError(msg)
    if statistic not in {"mean", "median"}:
        msg = f"unsupported criterion.statistic in {path}: {statistic!r}"
        raise ValueError(msg)

    commands: dict[str, tuple[str, ...]] = {}
    for field in ("baseline_command", "current_command"):
        value = data.get(field)
        if not isinstance(value, list) or not value or not all(isinstance(part, str) and part for part in value):
            msg = f"invalid or missing criterion.{field} in {path}"
            raise ValueError(msg)
        commands[field] = tuple(cast("list[str]", value))

    return CriterionProvenance(
        suite=cast("BenchmarkSuite", suite),
        scope=cast("ComparisonScope", scope),
        statistic=cast("Statistic", statistic),
        sample=sample,
        criterion_version=criterion_version,
        baseline_command=commands["baseline_command"],
        current_command=commands["current_command"],
    )


def _validate_validation_metadata(data: dict[str, object], *, path: Path) -> None:
    """Require fixture validation for both compared revisions."""
    command = data.get("command")
    if command != ["just", "test-bench-inputs"]:
        msg = f"invalid validation.command in {path}: expected ['just', 'test-bench-inputs']"
        raise ValueError(msg)
    for field in ("current_revision", "baseline_revision"):
        value = _required_metadata_string(data, field, path)
        if value != "passed":
            msg = f"validation.{field} in {path} must be 'passed', got {value!r}"
            raise ValueError(msg)
    for field in ("current_commit", "baseline_commit"):
        _required_metadata_string(data, field, path)
    for field in ("current_source_state_sha256", "baseline_source_state_sha256"):
        _required_sha256(data, field, path)
    for field in ("current_git_clean", "baseline_git_clean"):
        if not isinstance(data.get(field), bool):
            msg = f"invalid or missing validation.{field} in {path}"
            raise TypeError(msg)
    compatibility = data.get("baseline_api_compatibility")
    if compatibility is not None and (not isinstance(compatibility, str) or not compatibility):
        msg = f"invalid validation.baseline_api_compatibility in {path}"
        raise TypeError(msg)


def _validate_baseline_api_compatibility(data: dict[str, object], *, baseline: str, path: Path) -> None:
    """Reject compatibility adapters attached to an unrelated baseline."""
    compatibility = data.get("baseline_api_compatibility")
    if compatibility == _V0_4_3_API_COMPATIBILITY and baseline != "v0.4.3":
        msg = f"validation.baseline_api_compatibility {_V0_4_3_API_COMPATIBILITY!r} in {path} is valid only for baseline 'v0.4.3', got {baseline!r}"
        raise ValueError(msg)


def _assess_change(baseline: CriterionEstimate, current: CriterionEstimate) -> ChangeAssessment:
    """Classify a change conservatively from non-overlapping Criterion intervals."""
    if not baseline.has_confidence_interval or not current.has_confidence_interval:
        return "unknown"

    if baseline.ci_lo_ns is None or baseline.ci_hi_ns is None or current.ci_lo_ns is None or current.ci_hi_ns is None:
        msg = "confidence-interval presence invariant violated"
        raise AssertionError(msg)
    if current.ci_hi_ns < baseline.ci_lo_ns:
        return "improvement"
    if current.ci_lo_ns > baseline.ci_hi_ns:
        return "regression"
    return "inconclusive"


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

            estimate = _read_estimate(est_path, stat)
            results.append(BenchResult(suite="exact", group=group, bench=bench, estimate=estimate))

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
            estimate = _read_estimate(est_path, stat)
            results.append(BenchResult(suite="vs_linalg", group=group_dir.name, bench=bench, estimate=estimate))

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
    policy: ComparisonPolicy,
) -> ComparisonCollection:
    """Compare exact results while retaining every missing expected row."""
    comparisons: list[Comparison] = []
    gaps: list[CoverageGap] = []
    unavailable_baseline_rows = _UNAVAILABLE_BASELINE_ROWS_BY_COMPATIBILITY.get(
        policy.baseline_api_compatibility or "",
        frozenset(),
    )

    for group, benches in EXACT_GROUPS.items():
        if policy.scope == "release-signal" and group not in EXACT_RELEASE_SIGNAL_GROUPS:
            continue

        group_dir = criterion_dir / group
        for bench in benches:
            new_path = group_dir / bench / "new" / "estimates.json"
            baseline_bench, base_path = _exact_baseline_path(group_dir, bench, baseline_name)
            missing_current = not new_path.exists()
            missing_baseline = not base_path.exists()
            baseline_unavailable = (group, bench) in unavailable_baseline_rows

            if baseline_unavailable:
                if missing_current:
                    gaps.append(
                        CoverageGap(
                            suite="exact",
                            group=group,
                            bench=bench,
                            baseline_bench=baseline_bench,
                            missing_current=True,
                            missing_baseline=False,
                        )
                    )
                else:
                    _read_estimate(new_path, stat)
                continue

            if missing_current or missing_baseline:
                gaps.append(
                    CoverageGap(
                        suite="exact",
                        group=group,
                        bench=bench,
                        baseline_bench=baseline_bench,
                        missing_current=missing_current,
                        missing_baseline=missing_baseline,
                    )
                )
                continue

            current = _read_estimate(new_path, stat)
            baseline = _read_estimate(base_path, stat)

            comparisons.append(
                Comparison(
                    suite="exact",
                    group=group,
                    bench=bench,
                    baseline=baseline,
                    current=current,
                    assessment=_assess_change(baseline, current),
                    baseline_bench=baseline_bench if baseline_bench != bench else None,
                )
            )

    expected_groups = [group for group in EXACT_GROUPS if policy.scope != "release-signal" or group in EXACT_RELEASE_SIGNAL_GROUPS]
    if not any((criterion_dir / group).is_dir() for group in expected_groups):
        gaps.append(_entire_suite_gap("exact"))

    return ComparisonCollection(comparisons=comparisons, gaps=gaps)


def _ordered_vs_linalg_comparison_benches(group_dir: Path, baseline_name: str, scope: str) -> list[str]:
    """Return expected or discovered comparison rows in stable order."""
    dim = _dim_from_vs_linalg_group(group_dir.name)
    if scope == "release-signal":
        present = set(VS_LINALG_LA_STACK_BENCHES)
        if dim is not None:
            present.update(VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM.get(dim, []))
    else:
        present = {
            child.name
            for child in group_dir.iterdir()
            if child.is_dir() and ((child / "new" / "estimates.json").exists() or (child / baseline_name / "estimates.json").exists())
        }

    ordered = [bench for bench in VS_LINALG_BENCH_ORDER if bench in present]
    extras = sorted(present.difference(VS_LINALG_BENCH_ORDER))
    return [*ordered, *extras]


def _vs_linalg_dimension_groups(criterion_dir: Path, scope: str) -> list[tuple[int, Path]]:
    """Return canonical or discovered dimension groups for a comparison."""
    if scope == "release-signal":
        return [(dim, criterion_dir / f"d{dim}") for dim in VS_LINALG_CANONICAL_DIMS]

    groups: list[tuple[int, Path]] = []
    for group_dir in criterion_dir.iterdir():
        if not group_dir.is_dir():
            continue
        dim = _dim_from_vs_linalg_group(group_dir.name)
        if dim is not None:
            groups.append((dim, group_dir))
    return groups


def _read_optional_estimate(estimates_json: Path, stat: str) -> CriterionEstimate | None:
    """Read an optional Criterion estimate."""
    if not estimates_json.exists():
        return None
    return _read_estimate(estimates_json, stat)


def _baseline_peer_estimates(
    group_dir: Path,
    bench: str,
    baseline_name: str,
    stat: str,
) -> tuple[CriterionEstimate | None, CriterionEstimate | None]:
    """Return last-release nalgebra/faer context for a la-stack vs_linalg bench."""
    peers = VS_LINALG_BASELINE_PEERS.get(bench)
    if peers is None:
        return (None, None)

    nalgebra_bench, faer_bench = peers
    nalgebra = _read_optional_estimate(group_dir / nalgebra_bench / baseline_name / "estimates.json", stat)
    faer = _read_optional_estimate(group_dir / faer_bench / baseline_name / "estimates.json", stat)
    return (nalgebra, faer)


def _comparison_bench_label(comparison: Comparison) -> str:
    """Return the display label for a comparison table row."""
    if comparison.baseline_bench is None:
        return comparison.bench
    return f"{comparison.bench} (vs {comparison.baseline_bench})"


def _collect_vs_linalg_comparisons(
    criterion_dir: Path,
    baseline_name: str,
    stat: str,
    policy: ComparisonPolicy,
) -> ComparisonCollection:
    """Compare vs_linalg results while retaining one-sided rows."""
    comparisons: list[Comparison] = []
    gaps: list[CoverageGap] = []
    unavailable_baseline_rows = _UNAVAILABLE_BASELINE_ROWS_BY_COMPATIBILITY.get(
        policy.baseline_api_compatibility or "",
        frozenset(),
    )
    dim_groups = _vs_linalg_dimension_groups(criterion_dir, policy.scope)

    for _dim, group_dir in sorted(dim_groups, key=lambda item: item[0]):
        expected_benches = _ordered_vs_linalg_comparison_benches(group_dir, baseline_name, policy.scope)
        for bench in expected_benches:
            new_path = group_dir / bench / "new" / "estimates.json"
            base_path = group_dir / bench / baseline_name / "estimates.json"
            missing_current = not new_path.exists()
            missing_baseline = not base_path.exists()
            baseline_unavailable = (group_dir.name, bench) in unavailable_baseline_rows

            if baseline_unavailable:
                if missing_current:
                    gaps.append(
                        CoverageGap(
                            suite="vs_linalg",
                            group=group_dir.name,
                            bench=bench,
                            baseline_bench=bench,
                            missing_current=True,
                            missing_baseline=False,
                        )
                    )
                else:
                    _read_estimate(new_path, stat)
                continue

            if missing_current or missing_baseline:
                gaps.append(
                    CoverageGap(
                        suite="vs_linalg",
                        group=group_dir.name,
                        bench=bench,
                        baseline_bench=bench,
                        missing_current=missing_current,
                        missing_baseline=missing_baseline,
                    )
                )
                continue

            current = _read_estimate(new_path, stat)
            baseline = _read_estimate(base_path, stat)
            baseline_nalgebra, baseline_faer = _baseline_peer_estimates(group_dir, bench, baseline_name, stat)

            comparisons.append(
                Comparison(
                    suite="vs_linalg",
                    group=group_dir.name,
                    bench=bench,
                    baseline=baseline,
                    current=current,
                    assessment=_assess_change(baseline, current),
                    baseline_nalgebra=baseline_nalgebra,
                    baseline_faer=baseline_faer,
                )
            )

    present_dim_groups = [group_dir for _, group_dir in dim_groups if group_dir.is_dir()]
    if not present_dim_groups:
        gaps.append(_entire_suite_gap("vs_linalg"))

    return ComparisonCollection(comparisons=comparisons, gaps=gaps)


def _entire_suite_gap(suite: str) -> CoverageGap:
    """Return an explicit marker for a wholly absent selected suite."""
    return CoverageGap(
        suite=suite,
        group="(entire suite)",
        bench="all selected rows",
        baseline_bench="all selected rows",
        missing_current=True,
        missing_baseline=True,
    )


def _snapshot_coverage_errors(
    criterion_dir: Path,
    *,
    sample: str,
    suite: str,
    scope: str,
) -> list[str]:
    """Return deterministic coverage errors for a selected snapshot scope."""
    errors: list[str] = []
    selected_suites = ("exact", "vs_linalg") if suite == "all" else (suite,)

    if "exact" in selected_suites:
        errors.extend(_exact_snapshot_coverage_errors(criterion_dir, sample=sample, scope=scope))

    if "vs_linalg" in selected_suites:
        errors.extend(_vs_linalg_snapshot_coverage_errors(criterion_dir, sample=sample, scope=scope))

    return errors


def _exact_snapshot_coverage_errors(criterion_dir: Path, *, sample: str, scope: str) -> list[str]:
    """Return exact-suite snapshot gaps."""
    if not any((criterion_dir / group).is_dir() for group in EXACT_GROUPS):
        return ["exact: entire selected suite is absent"]
    if scope != "release-signal":
        return []
    return [
        f"exact: missing {group}/{bench}/{sample}/estimates.json"
        for group, benches in EXACT_GROUPS.items()
        for bench in benches
        if not (criterion_dir / group / bench / sample / "estimates.json").is_file()
    ]


def _vs_linalg_snapshot_coverage_errors(criterion_dir: Path, *, sample: str, scope: str) -> list[str]:
    """Return vs_linalg snapshot gaps, including canonical dimensions."""
    discovered_dims = {dim for child in criterion_dir.iterdir() if child.is_dir() and (dim := _dim_from_vs_linalg_group(child.name)) is not None}
    if not discovered_dims:
        return ["vs_linalg: entire selected suite is absent"]
    if scope != "release-signal":
        return []

    errors: list[str] = []
    for dim in VS_LINALG_CANONICAL_DIMS:
        group_dir = criterion_dir / f"d{dim}"
        if not group_dir.is_dir():
            errors.append(f"vs_linalg: missing canonical dimension d{dim}")
            continue
        expected = set(VS_LINALG_LA_STACK_BENCHES)
        expected.update(VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM.get(dim, []))
        errors.extend(
            f"vs_linalg: missing d{dim}/{bench}/{sample}/estimates.json"
            for bench in VS_LINALG_BENCH_ORDER
            if bench in expected and not (group_dir / bench / sample / "estimates.json").is_file()
        )
    return errors


def _collect_comparisons(
    criterion_dir: Path,
    baseline_name: str,
    stat: str,
    suite: str = "all",
    policy: ComparisonPolicy = _DEFAULT_COMPARISON_POLICY,
) -> ComparisonCollection:
    """Compare current results and report every expected coverage gap."""
    comparisons: list[Comparison] = []
    gaps: list[CoverageGap] = []
    if suite in ("all", "exact"):
        exact = _collect_exact_comparisons(criterion_dir, baseline_name, stat, policy)
        comparisons.extend(exact.comparisons)
        gaps.extend(exact.gaps)
    if suite in ("all", "vs_linalg"):
        vs_linalg = _collect_vs_linalg_comparisons(criterion_dir, baseline_name, stat, policy)
        comparisons.extend(vs_linalg.comparisons)
        gaps.extend(vs_linalg.gaps)
    return ComparisonCollection(comparisons=comparisons, gaps=gaps)


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


def _format_pct(pct: float, assessment: ChangeAssessment) -> str:
    """Format point-estimate change without implying statistical significance."""
    del assessment
    return f"{pct:+.1f}%"


def _format_confidence_interval(estimate: CriterionEstimate) -> str:
    """Format a Criterion confidence interval without inventing missing bounds."""
    if estimate.ci_lo_ns is None or estimate.ci_hi_ns is None:
        return "unavailable"
    return f"[{_format_time(estimate.ci_lo_ns)}, {_format_time(estimate.ci_hi_ns)}]"


def _format_estimate(estimate: CriterionEstimate) -> str:
    """Format a point estimate together with its Criterion interval."""
    return f"{_format_time(estimate.point_ns)} {_format_confidence_interval(estimate)}"


def _assessment_label(assessment: ChangeAssessment) -> str:
    """Return a concise, explicit confidence-interval assessment label."""
    labels: dict[ChangeAssessment, str] = {
        "improvement": "faster point estimate; marginal CIs separated",
        "regression": "slower point estimate; marginal CIs separated",
        "inconclusive": "marginal CIs overlap",
        "unknown": "marginal CI unavailable",
    }
    return labels[assessment]


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
    # exact_d3 -> "D=3", exact_random_corpus_d3 ->
    # "Random corpus D=3", exact_near_singular_3x3 ->
    # "Near-singular 3x3", exact_hilbert_4x4 -> "Hilbert 4x4", etc.
    if group.startswith("exact_random_corpus_d"):
        return f"Random corpus D={group.removeprefix('exact_random_corpus_d')}"
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
                f"| Benchmark | {stat_label} | Criterion CI |",
                "|-----------|-------:|-------:|",
            ]
            for r in items:
                ci_range = _format_confidence_interval(r.estimate)
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
                header = (
                    f"| Benchmark | {baseline_name} (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | "
                    f"{baseline_name} nalgebra | {baseline_name} faer |"
                )
                lines.extend(
                    [
                        header,
                        "|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|",
                    ]
                )
            else:
                lines.extend(
                    [
                        f"| Benchmark | {baseline_name} (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |",
                        "|-----------|-------:|-------:|-------:|:-----------|--------:|",
                    ]
                )
            for c in items:
                cells = [
                    _comparison_bench_label(c),
                    _format_estimate(c.baseline),
                    _format_estimate(c.current),
                    _format_pct(c.pct_change, c.assessment),
                    _assessment_label(c.assessment),
                    f"{c.speedup:.2f}x",
                ]
                if has_peer_context:
                    cells.extend(
                        [
                            _format_estimate(c.baseline_nalgebra) if c.baseline_nalgebra is not None else "—",
                            _format_estimate(c.baseline_faer) if c.baseline_faer is not None else "—",
                        ]
                    )
                lines.append(f"| {' | '.join(cells)} |")
            sections.append("\n".join(lines))

    return "\n\n".join(sections)


def _coverage_table(gaps: list[CoverageGap], baseline_name: str) -> str:
    """Render missing current/baseline rows in deterministic collection order."""
    if not gaps:
        return ""

    lines = [
        "## Incomplete Comparison Coverage",
        "",
        "These expected rows were not classified because one or both Criterion samples are missing.",
        "",
        "| Suite | Group | Benchmark | Missing sample(s) |",
        "|:------|:------|:----------|:------------------|",
    ]
    for gap in gaps:
        missing: list[str] = []
        if gap.missing_current:
            missing.append("current (`new`)")
        if gap.missing_baseline:
            missing.append(f"baseline (`{baseline_name}`)")
        bench_label = gap.bench
        if gap.baseline_bench != gap.bench:
            bench_label = f"{gap.bench} (baseline row: {gap.baseline_bench})"
        lines.append(f"| {_suite_heading(gap.suite)} | {gap.group} | {bench_label} | {', '.join(missing)} |")

    return "\n".join(lines)


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
    except (
        ExecutableNotFoundError,
        OSError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
        pass
    try:
        result = run_git_command(["--no-pager", "rev-parse", "--abbrev-ref", "HEAD"], cwd=root)
        branch = result.stdout.strip()
    except (
        ExecutableNotFoundError,
        OSError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
        pass
    return short_hash, branch


def _get_git_source_date(root: Path) -> str:
    """Return the reproducible source revision timestamp normalized to UTC."""
    try:
        result = run_git_command(["--no-pager", "show", "-s", "--format=%cI", "HEAD"], cwd=root)
    except (
        ExecutableNotFoundError,
        OSError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
        return "unknown"
    value = result.stdout.strip()
    if not value:
        return "unknown"
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return "unknown"
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        return "unknown"
    return parsed.astimezone(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")


def _generate_markdown(
    root: Path,
    table: str,
    settings: ReportSettings,
) -> str:
    """Generate the complete benchmark report content."""
    version = _read_cargo_version(root)
    short_hash, branch = _get_git_info(root)
    source_date = _get_git_source_date(root)

    lines = [
        "# Benchmark Performance",
        "",
        f"**la-stack** v{version} · `{short_hash}` ({branch})",
        f"**Source revision timestamp**: {source_date} (deterministic report metadata; not the benchmark measurement time)",
        "**Benchmark measurement timestamp**: not recorded by Criterion; use the provenance below to identify the measured revisions and environment.",
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
        lines.append(
            "Negative point-estimate change means the current point estimate is smaller; "
            "a baseline/current point-estimate ratio above 1.00 has the same meaning."
        )
        lines.append(
            "The CI-relation column reports only whether the two marginal Criterion intervals overlap. "
            "These are not paired confidence intervals for the change, so the report makes no "
            "statistical-significance or performance-improvement claim from interval separation."
        )
        lines.append("")
        if settings.harness_provenance is None:
            lines.append(
                "**Reproducibility provenance**: unavailable — these Criterion samples predate provenance capture or were produced outside the "
                "publication pipeline. CPU, OS, rustc, commit, dependency lock, harness digest, Criterion configuration, and two-revision "
                "fixture validation are unknown."
            )
        else:
            lines.extend(_provenance_markdown(settings.harness_provenance))
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


def _provenance_markdown(provenance: HarnessProvenance) -> list[str]:
    """Render validated provenance without implying facts absent from metadata."""
    if provenance.schema == 1:
        return [
            "**Harness provenance**: legacy shared current harness metadata from "
            f"`.la-stack-benchmark-harness.json` (baseline source `{provenance.baseline}`, SHA-256 `{provenance.sha256}`). "
            "CPU, OS, rustc, commits, Criterion configuration, and fixture-gate results were not recorded."
        ]

    if provenance.measurement is None or provenance.publication is None or provenance.criterion is None or provenance.validation is None:
        msg = "schema-2 benchmark provenance invariant violated"
        raise AssertionError(msg)

    measurement = provenance.measurement
    publication = provenance.publication
    criterion = provenance.criterion
    validation = provenance.validation
    if provenance.mode == "shared-current-harness":
        correctness_gate = (
            "- Correctness gate: `just test-bench-inputs` passed against both the current and baseline revisions using the shared current fixture harness."
        )
    else:
        correctness_gate = (
            "- Correctness gate: `just test-bench-inputs` passed for both referenced source revisions "
            "during publication; this gate validates benchmark inputs separately from the historical timing measurements."
        )
    lines = ["### Reproducibility Provenance", ""]
    if measurement["status"] == "recorded":
        lines.extend(
            [
                "**Measurement environment**: recorded for both samples under one shared current harness.",
                "",
                f"- CPU: `{measurement['cpu']}`",
                f"- OS: `{measurement['os']}`",
                f"- rustc: `{measurement['rustc']}`",
                f"- Current commit: `{measurement['current_commit']}`",
                f"- Current Git clean: `{str(measurement['current_git_clean']).lower()}`",
                f"- Current source-state SHA-256: `{measurement['current_source_state_sha256']}`",
                f"- Baseline commit: `{measurement['baseline_commit']}`",
                f"- Baseline Git clean: `{str(measurement['baseline_git_clean']).lower()}`",
                f"- Baseline source-state SHA-256: `{measurement['baseline_source_state_sha256']}`",
                f"- Cargo.lock SHA-256: `{measurement['cargo_lock_sha256']}`",
                f"- Benchmark harness SHA-256: `{measurement['harness_sha256']}`",
            ]
        )
    else:
        lines.extend(
            [
                f"**Measurement environment**: unavailable — {measurement['reason']}",
                "",
                "The publication environment below validates report generation, not the historical timing environment.",
            ]
        )

    lines.extend(
        [
            "",
            f"- Publication CPU: `{publication['cpu']}`",
            f"- Publication OS: `{publication['os']}`",
            f"- Publication rustc: `{publication['rustc']}`",
            f"- Publication commit: `{publication['commit']}`",
            f"- Publication Git clean: `{str(publication['git_clean']).lower()}`",
            f"- Publication source-state SHA-256: `{publication['source_state_sha256']}`",
            f"- Publication Cargo.lock SHA-256: `{publication['cargo_lock_sha256']}`",
            f"- Publication harness SHA-256: `{publication['harness_sha256']}`",
            f"- Criterion suite/scope: `{criterion.suite}` / `{criterion.scope}`",
            f"- Criterion statistic/sample: `{criterion.statistic}` / `{criterion.sample}`",
            f"- Criterion dependency version: `{criterion.criterion_version}`",
            f"- Baseline command: `{' '.join(criterion.baseline_command)}`",
            f"- Current command: `{' '.join(criterion.current_command)}`",
            correctness_gate,
            (
                f"- Validated current revision: `{validation['current_commit']}` "
                f"(Git clean: `{str(validation['current_git_clean']).lower()}`, "
                f"source-state SHA-256: `{validation['current_source_state_sha256']}`)"
            ),
            (
                f"- Validated baseline revision: `{validation['baseline_commit']}` "
                f"(Git clean: `{str(validation['baseline_git_clean']).lower()}`, "
                f"source-state SHA-256: `{validation['baseline_source_state_sha256']}`)"
            ),
        ]
    )
    compatibility = validation.get("baseline_api_compatibility")
    if isinstance(compatibility, str) and compatibility != "none":
        lines.append(
            f"- Baseline API compatibility: `{compatibility}` selects only source-compatible benchmark calls; "
            "rows outside the baseline's correctness domain remain explicitly unavailable."
        )
        if compatibility == _V0_4_3_API_COMPATIBILITY and criterion.suite in {"all", "vs_linalg"}:
            lines.append(
                "- Baseline-unavailable rows: `d8/la_stack_det_from_lu_balanced_range` and "
                "`d8/la_stack_det_from_ldlt_balanced_range` were not timed because v0.4.3 returns zero for a fixture "
                "whose exact determinant is one; current samples remain required, but no speedup is claimed."
            )
        if compatibility == _V0_4_3_API_COMPATIBILITY and criterion.suite in {"all", "exact"}:
            lines.append(
                "- Baseline-unavailable rows: `exact_d2/det_direct_with_errbound`, "
                "`exact_d3/det_direct_with_errbound`, and `exact_d4/det_direct_with_errbound` were not timed because "
                "v0.4.3 predates the paired API; the comparable `det_errbound` baselines remain required."
            )
    return lines


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


def _comparison_policy(scope: str, provenance: HarnessProvenance | None) -> ComparisonPolicy:
    """Build comparison coverage policy from validated provenance."""
    if provenance is None or provenance.validation is None:
        return ComparisonPolicy(scope=scope)
    compatibility = provenance.validation.get("baseline_api_compatibility")
    return ComparisonPolicy(
        scope=scope,
        baseline_api_compatibility=compatibility if isinstance(compatibility, str) else None,
    )


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


def main(argv: list[str] | None = None) -> int:  # noqa: C901, PLR0911
    """Generate a benchmark snapshot or comparison report from CLI arguments."""
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    stat = cast("Statistic", args.stat)
    suite = cast("BenchmarkSuite", args.suite)
    scope = cast("ComparisonScope", args.scope)
    selection = CriterionSelection(suite=suite, scope=scope, statistic=stat, sample="new")

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
    harness_provenance: HarnessProvenance | None = None

    if baseline_name:
        try:
            harness_provenance = _read_harness_provenance(
                criterion_dir,
                expected_baseline=baseline_name,
                expected=selection,
            )
        except (OSError, TypeError, ValueError) as err:
            print(f"Invalid benchmark harness provenance: {err}", file=sys.stderr)
            return 2

        try:
            collection = _collect_comparisons(
                criterion_dir,
                baseline_name,
                stat,
                suite,
                _comparison_policy(scope, harness_provenance),
            )
        except (OSError, KeyError, TypeError, ValueError) as err:
            print(f"Invalid Criterion estimate data: {err}", file=sys.stderr)
            return 2
        if collection.gaps:
            print(
                f"Incomplete benchmark coverage: {len(collection.gaps)} required comparison row(s) are missing; report publication aborted.",
                file=sys.stderr,
            )
            print(_coverage_table(collection.gaps, baseline_name), file=sys.stderr)
            if not collection.comparisons:
                print(
                    f"No comparison data found for baseline '{baseline_name}'.\n"
                    f"Save a baseline first:\n  {_save_baseline_hint(args.suite, baseline_name)}\n"
                    f"Then run benchmarks:\n  {_run_bench_hint(args.suite)}\n",
                    file=sys.stderr,
                )
            return 2
        if not collection.comparisons:
            print(
                f"No comparison data found for baseline '{baseline_name}'.\n"
                f"Save a baseline first:\n  {_save_baseline_hint(args.suite, baseline_name)}\n"
                f"Then run benchmarks:\n  {_run_bench_hint(args.suite)}\n",
                file=sys.stderr,
            )
            return 2
        table = _comparison_tables(collection.comparisons, baseline_name)
    else:
        coverage_errors = _snapshot_coverage_errors(
            criterion_dir,
            sample="new",
            suite=suite,
            scope=scope,
        )
        if coverage_errors:
            print("Incomplete benchmark coverage; snapshot publication aborted:", file=sys.stderr)
            for error in coverage_errors:
                print(f"  - {error}", file=sys.stderr)
            return 2
        try:
            results = _collect_results(criterion_dir, "new", stat, suite)
        except (OSError, KeyError, TypeError, ValueError) as err:
            print(f"Invalid Criterion estimate data: {err}", file=sys.stderr)
            return 2
        if not results:
            print(
                f"No benchmark results found.\nRun benchmarks first:\n  {_run_bench_hint(args.suite)}\n",
                file=sys.stderr,
            )
            return 2
        table = _snapshot_tables(results, stat)

    settings = ReportSettings(
        baseline_name=baseline_name,
        stat=stat,
        suite=suite,
        scope=scope,
        harness_provenance=harness_provenance,
    )
    md = _generate_markdown(root, table, settings)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(md, encoding="utf-8")
    print(f"📊 Wrote {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
