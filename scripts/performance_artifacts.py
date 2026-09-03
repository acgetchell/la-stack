"""Schema-versioned performance-comparison CSV and provenance artifacts."""

import csv
import hashlib
import io
import json
import math
import os
import tempfile
from collections.abc import Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal, cast

if TYPE_CHECKING:
    from collections.abc import Iterator

SCHEMA_VERSION = 1
SUITES = ("all", "exact", "vs_linalg")
SCOPES = ("release-signal", "all-benches")
COVERAGE_STATES = ("comparable", "current-only", "baseline-only")

type CoverageState = Literal["comparable", "current-only", "baseline-only"]

CSV_COLUMNS = (
    "schema_version",
    "suite",
    "scope",
    "benchmark_id",
    "group",
    "benchmark",
    "baseline_benchmark",
    "coverage_status",
    "coverage_note",
    "baseline_median_ns",
    "baseline_ci_lower_ns",
    "baseline_ci_upper_ns",
    "current_median_ns",
    "current_ci_lower_ns",
    "current_ci_upper_ns",
    "baseline_nalgebra_median_ns",
    "baseline_nalgebra_ci_lower_ns",
    "baseline_nalgebra_ci_upper_ns",
    "baseline_faer_median_ns",
    "baseline_faer_ci_lower_ns",
    "baseline_faer_ci_upper_ns",
)


@dataclass(frozen=True, slots=True)
class TimingEstimate:
    """A finite positive median and its complete confidence interval."""

    median_ns: float
    ci_lower_ns: float
    ci_upper_ns: float

    def __post_init__(self) -> None:
        """Reject non-finite, non-positive, or reversed timing intervals."""
        for field, value in (
            ("median_ns", self.median_ns),
            ("ci_lower_ns", self.ci_lower_ns),
            ("ci_upper_ns", self.ci_upper_ns),
        ):
            if not math.isfinite(value) or value <= 0:
                msg = f"{field} must be finite and positive: {value!r}"
                raise ValueError(msg)
        if self.ci_lower_ns > self.ci_upper_ns:
            msg = f"confidence interval must be ordered: {self.ci_lower_ns} <= {self.ci_upper_ns}"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class PerformanceRow:
    """One validated row in the release-performance reporting dataset."""

    suite: str
    scope: str
    benchmark_id: str
    group: str
    benchmark: str
    baseline_benchmark: str
    coverage_status: CoverageState
    coverage_note: str
    baseline: TimingEstimate | None
    current: TimingEstimate | None
    baseline_nalgebra: TimingEstimate | None = None
    baseline_faer: TimingEstimate | None = None

    def __post_init__(self) -> None:  # noqa: C901
        """Preserve identity, coverage, and optional-timing invariants."""
        if self.suite not in SUITES:
            msg = f"unsupported suite: {self.suite!r}"
            raise ValueError(msg)
        if self.scope not in SCOPES:
            msg = f"unsupported scope: {self.scope!r}"
            raise ValueError(msg)
        for field, value in (
            ("benchmark_id", self.benchmark_id),
            ("group", self.group),
            ("benchmark", self.benchmark),
            ("baseline_benchmark", self.baseline_benchmark),
        ):
            if not value.strip():
                msg = f"{field} must not be empty"
                raise ValueError(msg)
        expected_id = f"{self.group}/{self.benchmark}"
        if self.benchmark_id != expected_id:
            msg = f"benchmark_id must be {expected_id!r}, got {self.benchmark_id!r}"
            raise ValueError(msg)
        if self.coverage_status not in COVERAGE_STATES:
            msg = f"unsupported coverage status: {self.coverage_status!r}"
            raise ValueError(msg)

        expected_presence = {
            "comparable": (True, True),
            "current-only": (False, True),
            "baseline-only": (True, False),
        }[self.coverage_status]
        observed_presence = (self.baseline is not None, self.current is not None)
        if observed_presence != expected_presence:
            msg = f"coverage status {self.coverage_status!r} requires baseline/current presence {expected_presence}, got {observed_presence}"
            raise ValueError(msg)
        if self.coverage_status == "comparable" and self.coverage_note:
            msg = "comparable rows must not contain a coverage note"
            raise ValueError(msg)
        if self.coverage_status != "comparable" and not self.coverage_note.strip():
            msg = f"{self.coverage_status} rows require a coverage note"
            raise ValueError(msg)
        if (self.baseline_nalgebra is not None or self.baseline_faer is not None) and self.baseline is None:
            msg = "baseline peer timings require a baseline la-stack timing"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ReleasePair:
    """The current and baseline package identifiers represented by a report."""

    current: str
    baseline: str

    def __post_init__(self) -> None:
        """Require two non-empty identifiers; local comparisons may match."""
        if not self.current.strip() or not self.baseline.strip():
            msg = "current and baseline releases must not be empty"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ReportSource:
    """Stable source metadata needed to reproduce the report header."""

    version: str
    commit: str
    ref: str
    revision_timestamp: str

    def __post_init__(self) -> None:
        """Reject incomplete source metadata before artifact publication."""
        for field, value in (
            ("version", self.version),
            ("commit", self.commit),
            ("ref", self.ref),
            ("revision_timestamp", self.revision_timestamp),
        ):
            if not value.strip() or value == "unknown":
                msg = f"report source {field} is incomplete: {value!r}"
                raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ArtifactContext:
    """Validated non-tabular metadata for a release-performance dataset."""

    release: ReleasePair
    statistic: Literal["median"]
    suite: str
    scope: str
    source: ReportSource
    benchmark_provenance: Mapping[str, object]

    def __post_init__(self) -> None:
        """Bind report settings and benchmark provenance to the release pair."""
        if self.statistic != "median":
            msg = f"unsupported performance-comparison statistic: {self.statistic!r}"
            raise ValueError(msg)
        if self.suite not in SUITES:
            msg = f"unsupported suite: {self.suite!r}"
            raise ValueError(msg)
        if self.scope not in SCOPES:
            msg = f"unsupported scope: {self.scope!r}"
            raise ValueError(msg)
        _validate_benchmark_provenance(
            self.benchmark_provenance,
            context=self,
        )
        object.__setattr__(self, "benchmark_provenance", freeze_mapping(self.benchmark_provenance))


@dataclass(frozen=True, slots=True)
class PerformanceBundle:
    """A complete validated release-performance dataset and its provenance."""

    context: ArtifactContext
    rows: tuple[PerformanceRow, ...]

    def __post_init__(self) -> None:
        """Require non-empty, unique rows bound to one report selection."""
        if not self.rows:
            msg = "release-performance dataset must contain at least one row"
            raise ValueError(msg)
        keys: set[tuple[str, str, str]] = set()
        for row in self.rows:
            if row.suite not in _selected_suites(self.context.suite):
                msg = f"row suite {row.suite!r} is outside selected suite {self.context.suite!r}"
                raise ValueError(msg)
            if row.scope != self.context.scope:
                msg = f"row scope {row.scope!r} does not match report scope {self.context.scope!r}"
                raise ValueError(msg)
            key = (row.suite, row.scope, row.benchmark_id)
            if key in keys:
                msg = f"duplicate benchmark key: {key!r}"
                raise ValueError(msg)
            keys.add(key)

    @property
    def sorted_rows(self) -> tuple[PerformanceRow, ...]:
        """Return rows in deterministic suite/group/benchmark order."""
        return tuple(sorted(self.rows, key=lambda row: (row.suite, row.group, row.benchmark, row.baseline_benchmark)))


@dataclass(frozen=True, slots=True)
class ArtifactPaths:
    """Adjacent CSV and JSON provenance destinations."""

    csv: Path
    provenance: Path

    def __post_init__(self) -> None:
        """Require distinct adjacent artifact files."""
        ensure_distinct_paths({"CSV": self.csv, "provenance": self.provenance})
        if self.csv.resolve(strict=False).parent != self.provenance.resolve(strict=False).parent:
            msg = "CSV and provenance sidecar must be adjacent"
            raise ValueError(msg)


def _selected_suites(suite: str) -> frozenset[str]:
    return frozenset({"exact", "vs_linalg"} if suite == "all" else {suite})


def _paths_alias(first: Path, second: Path) -> bool:
    """Return whether two paths resolve to the same filesystem target."""
    if first.exists() and second.exists() and first.samefile(second):
        return True
    return first.resolve(strict=False) == second.resolve(strict=False)


def ensure_distinct_paths(paths: Mapping[str, Path]) -> None:
    """Reject lexical, symlink, or existing-file aliases among named paths."""
    items = tuple(paths.items())
    for index, (first_name, first_path) in enumerate(items):
        for second_name, second_path in items[index + 1 :]:
            if _paths_alias(first_path, second_path):
                msg = f"{first_name} and {second_name} must use distinct paths: {first_path}"
                raise ValueError(msg)


def freeze_json(value: object) -> object:
    """Recursively detach and freeze JSON-shaped provenance data."""
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): freeze_json(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(freeze_json(item) for item in value)
    return value


def freeze_mapping(data: Mapping[str, object]) -> Mapping[str, object]:
    """Detach and freeze a JSON-shaped mapping while preserving its mapping type contract."""
    frozen = freeze_json(data)
    if not isinstance(frozen, Mapping):
        msg = "provenance mapping invariant violated"
        raise TypeError(msg)
    return cast("Mapping[str, object]", frozen)


def _thaw_json(value: object) -> object:
    """Convert frozen JSON-shaped data back to serializer-native containers."""
    if isinstance(value, Mapping):
        return {str(key): _thaw_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(item) for item in value]
    return value


def _required_provenance_object(data: Mapping[str, object], field: str, *, context: str) -> dict[str, object]:
    value = data.get(field)
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        msg = f"benchmark provenance requires a {field} object in {context}"
        raise TypeError(msg)
    return cast("dict[str, object]", value)


def _required_provenance_string(data: Mapping[str, object], field: str, *, context: str) -> str:
    value = data.get(field)
    if not isinstance(value, str) or not value.strip():
        msg = f"benchmark provenance requires a non-empty {context}.{field} string"
        raise ValueError(msg)
    return value


def _required_provenance_sha256(data: Mapping[str, object], field: str, *, context: str) -> str:
    value = _required_provenance_string(data, field, context=context)
    if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        msg = f"benchmark provenance {context}.{field} must be a lowercase SHA-256 digest"
        raise ValueError(msg)
    return value


def _required_provenance_bool(data: Mapping[str, object], field: str, *, context: str) -> bool:
    value = data.get(field)
    if not isinstance(value, bool):
        msg = f"benchmark provenance requires a boolean {context}.{field}"
        raise TypeError(msg)
    return value


def _required_provenance_command(data: Mapping[str, object], field: str, *, context: str) -> tuple[str, ...]:
    value = data.get(field)
    if not isinstance(value, list) or not value or not all(isinstance(part, str) and part for part in value):
        msg = f"benchmark provenance requires a non-empty {context}.{field} command"
        raise ValueError(msg)
    return tuple(cast("list[str]", value))


def _validate_environment_provenance(data: Mapping[str, object], *, context: str) -> None:
    for field in ("cpu", "os", "rustc", "commit"):
        _required_provenance_string(data, field, context=context)
    for field in ("cargo_lock_sha256", "harness_sha256", "source_state_sha256"):
        _required_provenance_sha256(data, field, context=context)
    if "benchmark_contract_sha256" in data:
        _required_provenance_sha256(data, "benchmark_contract_sha256", context=context)
    _required_provenance_bool(data, "git_clean", context=context)
    gate = _required_provenance_string(data, "correctness_gate", context=context)
    if gate != "passed":
        msg = f"benchmark provenance {context}.correctness_gate must be 'passed', got {gate!r}"
        raise ValueError(msg)


def _require_matching_provenance(
    first: Mapping[str, object],
    second: Mapping[str, object],
    field: str,
    *,
    first_context: str,
    second_context: str,
) -> None:
    if first.get(field) != second.get(field):
        msg = f"benchmark provenance {first_context}.{field} {first.get(field)!r} does not match {second_context}.{field} {second.get(field)!r}"
        raise ValueError(msg)


def _validate_criterion_provenance(criterion: Mapping[str, object], *, context: ArtifactContext) -> None:
    """Validate complete Criterion settings and commands."""
    expected = {
        "statistic": context.statistic,
        "suite": context.suite,
        "scope": context.scope,
        "sample": "new",
    }
    for field, expected_value in expected.items():
        if criterion.get(field) != expected_value:
            msg = f"benchmark provenance criterion.{field} {criterion.get(field)!r} does not match report value {expected_value!r}"
            raise ValueError(msg)
    _required_provenance_string(criterion, "criterion_version", context="criterion")
    _required_provenance_command(criterion, "baseline_command", context="criterion")
    _required_provenance_command(criterion, "current_command", context="criterion")


def _validate_measurement_provenance(measurement: Mapping[str, object], *, mode: object) -> str:
    """Validate recorded or explicitly unavailable measurement metadata."""
    measurement_status = _required_provenance_string(measurement, "status", context="measurement")
    if measurement_status == "recorded":
        if mode != "shared-current-harness":
            msg = "recorded benchmark measurement provenance requires shared-current-harness mode"
            raise ValueError(msg)
        for field in ("cpu", "os", "rustc", "current_commit", "baseline_commit", "baseline_api_compatibility"):
            _required_provenance_string(measurement, field, context="measurement")
        if cast("str", measurement["cpu"]).casefold() == "unavailable":
            msg = "recorded benchmark measurement provenance requires an identified CPU model"
            raise ValueError(msg)
        for field in (
            "cargo_lock_sha256",
            "harness_sha256",
            "current_source_state_sha256",
            "baseline_source_state_sha256",
        ):
            _required_provenance_sha256(measurement, field, context="measurement")
        if "benchmark_contract_sha256" in measurement:
            _required_provenance_sha256(measurement, "benchmark_contract_sha256", context="measurement")
        for field in ("current_git_clean", "baseline_git_clean"):
            _required_provenance_bool(measurement, field, context="measurement")
    elif measurement_status == "unavailable":
        _required_provenance_string(measurement, "reason", context="measurement")
    else:
        msg = f"unsupported benchmark provenance measurement.status: {measurement_status!r}"
        raise ValueError(msg)
    return measurement_status


def _expected_baseline_api_compatibility(*, current: str, baseline: str) -> str:
    """Return the compatibility adapter required by one release pair."""
    if baseline == "v0.4.3":
        return "la_stack_v0_4_3_api"
    if baseline in {"v0.4.4", "v0.4.5"} and current not in {"v0.4.4", "v0.4.5"}:
        return "la_stack_pre_rational_input_api"
    return "none"


def _validate_validation_provenance(
    validation: Mapping[str, object],
    *,
    current: str,
    baseline: str,
) -> str:
    """Validate fixture-gate evidence for both compared revisions."""
    if validation.get("command") != ["just", "test-bench-inputs"]:
        msg = "benchmark provenance validation.command must be ['just', 'test-bench-inputs']"
        raise ValueError(msg)
    for field in ("current_revision", "baseline_revision"):
        result = _required_provenance_string(validation, field, context="validation")
        if result != "passed":
            msg = f"benchmark provenance validation.{field} must be 'passed', got {result!r}"
            raise ValueError(msg)
    for field in ("current_commit", "baseline_commit"):
        _required_provenance_string(validation, field, context="validation")
    for field in ("current_source_state_sha256", "baseline_source_state_sha256"):
        _required_provenance_sha256(validation, field, context="validation")
    for field in ("current_git_clean", "baseline_git_clean"):
        _required_provenance_bool(validation, field, context="validation")
    compatibility = _required_provenance_string(validation, "baseline_api_compatibility", context="validation")
    expected_compatibility = _expected_baseline_api_compatibility(
        current=current,
        baseline=baseline,
    )
    if compatibility != expected_compatibility:
        msg = f"benchmark provenance validation.baseline_api_compatibility must be {expected_compatibility!r} for baseline {baseline!r}, got {compatibility!r}"
        raise ValueError(msg)
    harness = _required_provenance_string(validation, "harness", context="validation")
    if harness != "shared-current":
        msg = f"benchmark provenance validation.harness must be 'shared-current', got {harness!r}"
        raise ValueError(msg)
    return compatibility


def _validate_current_revision_consistency(
    publication: Mapping[str, object],
    validation: Mapping[str, object],
) -> None:
    """Bind validation of the current revision to publication metadata."""
    for field in ("commit", "git_clean", "source_state_sha256"):
        validation_field = f"current_{field}" if field != "commit" else "current_commit"
        publication_field = field
        if validation.get(validation_field) != publication.get(publication_field):
            msg = (
                f"benchmark provenance validation.{validation_field} {validation.get(validation_field)!r} does not match "
                f"publication.{publication_field} {publication.get(publication_field)!r}"
            )
            raise ValueError(msg)


def _validate_recorded_measurement_consistency(
    measurement: Mapping[str, object],
    publication: Mapping[str, object],
    validation: Mapping[str, object],
    *,
    compatibility: str,
) -> None:
    """Bind recorded timing metadata to publication and validation evidence."""
    for field in ("cpu", "os", "rustc", "cargo_lock_sha256", "harness_sha256"):
        _require_matching_provenance(
            measurement,
            publication,
            field,
            first_context="measurement",
            second_context="publication",
        )
    if "benchmark_contract_sha256" in measurement or "benchmark_contract_sha256" in publication:
        _require_matching_provenance(
            measurement,
            publication,
            "benchmark_contract_sha256",
            first_context="measurement",
            second_context="publication",
        )
    for measurement_field, publication_field in (
        ("current_commit", "commit"),
        ("current_git_clean", "git_clean"),
        ("current_source_state_sha256", "source_state_sha256"),
    ):
        if measurement.get(measurement_field) != publication.get(publication_field):
            msg = (
                f"benchmark provenance measurement.{measurement_field} {measurement.get(measurement_field)!r} does not match "
                f"publication.{publication_field} {publication.get(publication_field)!r}"
            )
            raise ValueError(msg)
    for field in ("baseline_commit", "baseline_git_clean", "baseline_source_state_sha256"):
        _require_matching_provenance(
            measurement,
            validation,
            field,
            first_context="measurement",
            second_context="validation",
        )
    if measurement.get("baseline_api_compatibility") != compatibility:
        msg = "benchmark provenance measurement.baseline_api_compatibility does not match validation.baseline_api_compatibility"
        raise ValueError(msg)


def _validate_benchmark_provenance(data: Mapping[str, object], *, context: ArtifactContext) -> None:
    """Reject incomplete or mismatched embedded benchmark provenance."""
    schema = data.get("schema")
    if type(schema) is not int or schema != 2:
        msg = f"benchmark provenance schema must be 2, got {data.get('schema')!r}"
        raise ValueError(msg)
    mode = data.get("mode")
    if mode not in {"shared-current-harness", "historical-assets"}:
        msg = f"unsupported benchmark provenance mode: {mode!r}"
        raise ValueError(msg)
    if data.get("baseline") != context.release.baseline:
        msg = f"benchmark provenance baseline {data.get('baseline')!r} does not match release baseline {context.release.baseline!r}"
        raise ValueError(msg)

    criterion = _required_provenance_object(data, "criterion", context="root")
    measurement = _required_provenance_object(data, "measurement", context="root")
    publication = _required_provenance_object(data, "publication", context="root")
    validation = _required_provenance_object(data, "validation", context="root")
    _validate_criterion_provenance(criterion, context=context)
    _validate_environment_provenance(publication, context="publication")
    source_commit = _required_provenance_string(publication, "commit", context="publication")
    if source_commit != "unavailable" and not source_commit.startswith(context.source.commit):
        msg = f"benchmark provenance publication.commit {source_commit!r} does not match report source commit {context.source.commit!r}"
        raise ValueError(msg)
    measurement_status = _validate_measurement_provenance(measurement, mode=mode)
    expected_status = "recorded" if mode == "shared-current-harness" else "unavailable"
    if measurement_status != expected_status:
        msg = f"benchmark provenance mode {mode!r} requires measurement.status {expected_status!r}, got {measurement_status!r}"
        raise ValueError(msg)
    compatibility = _validate_validation_provenance(
        validation,
        current=context.release.current,
        baseline=context.release.baseline,
    )
    _validate_current_revision_consistency(publication, validation)
    if measurement_status == "recorded":
        _validate_recorded_measurement_consistency(
            measurement,
            publication,
            validation,
            compatibility=compatibility,
        )


def _timing_fields(prefix: str, estimate: TimingEstimate | None) -> dict[str, str]:
    if estimate is None:
        return {
            f"{prefix}_median_ns": "",
            f"{prefix}_ci_lower_ns": "",
            f"{prefix}_ci_upper_ns": "",
        }
    return {
        f"{prefix}_median_ns": format(estimate.median_ns, ".17g"),
        f"{prefix}_ci_lower_ns": format(estimate.ci_lower_ns, ".17g"),
        f"{prefix}_ci_upper_ns": format(estimate.ci_upper_ns, ".17g"),
    }


def _row_to_csv(row: PerformanceRow) -> dict[str, str]:
    values = {
        "schema_version": str(SCHEMA_VERSION),
        "suite": row.suite,
        "scope": row.scope,
        "benchmark_id": row.benchmark_id,
        "group": row.group,
        "benchmark": row.benchmark,
        "baseline_benchmark": row.baseline_benchmark,
        "coverage_status": row.coverage_status,
        "coverage_note": row.coverage_note,
    }
    values.update(_timing_fields("baseline", row.baseline))
    values.update(_timing_fields("current", row.current))
    values.update(_timing_fields("baseline_nalgebra", row.baseline_nalgebra))
    values.update(_timing_fields("baseline_faer", row.baseline_faer))
    return values


def _serialize_csv(bundle: PerformanceBundle) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=CSV_COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in bundle.sorted_rows:
        writer.writerow(_row_to_csv(row))
    return output.getvalue().encode("utf-8")


def _serialize_provenance(bundle: PerformanceBundle, csv_payload: bytes) -> bytes:
    context = bundle.context
    payload = {
        "benchmark_provenance": _thaw_json(context.benchmark_provenance),
        "csv": {
            "columns": list(CSV_COLUMNS),
            "row_count": len(bundle.rows),
            "sha256": hashlib.sha256(csv_payload).hexdigest(),
        },
        "release": {
            "baseline": context.release.baseline,
            "current": context.release.current,
        },
        "report": {
            "scope": context.scope,
            "source": {
                "commit": context.source.commit,
                "ref": context.source.ref,
                "revision_timestamp": context.source.revision_timestamp,
                "version": context.source.version,
            },
            "statistic": context.statistic,
            "suite": context.suite,
        },
        "schema_version": SCHEMA_VERSION,
    }
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def serialize_bundle(bundle: PerformanceBundle) -> tuple[bytes, bytes]:
    """Serialize and self-validate one deterministic artifact pair."""
    csv_payload = _serialize_csv(bundle)
    provenance_payload = _serialize_provenance(bundle, csv_payload)
    parsed = load_bundle_bytes(csv_payload, provenance_payload, source="serialized artifact pair")
    if parsed != PerformanceBundle(context=bundle.context, rows=bundle.sorted_rows):
        msg = "release-performance artifact round trip changed the validated dataset"
        raise ValueError(msg)
    return csv_payload, provenance_payload


def _parse_required_string(data: Mapping[str, object], field: str, *, source: str) -> str:
    value = data.get(field)
    if not isinstance(value, str) or not value.strip():
        msg = f"invalid or missing {field} in {source}"
        raise ValueError(msg)
    return value


def _parse_required_object(data: Mapping[str, object], field: str, *, source: str) -> dict[str, object]:
    value = data.get(field)
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        msg = f"invalid or missing {field} object in {source}"
        raise ValueError(msg)
    return cast("dict[str, object]", value)


def _parse_float(value: str, field: str, *, row_number: int, source: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        msg = f"{field} in CSV row {row_number} of {source} is not numeric: {value!r}"
        raise ValueError(msg) from exc
    if not math.isfinite(parsed) or parsed <= 0:
        msg = f"{field} in CSV row {row_number} of {source} must be finite and positive: {value!r}"
        raise ValueError(msg)
    return parsed


def _parse_timing(row: Mapping[str, str], prefix: str, *, row_number: int, source: str) -> TimingEstimate | None:
    field_names = (
        f"{prefix}_median_ns",
        f"{prefix}_ci_lower_ns",
        f"{prefix}_ci_upper_ns",
    )
    values = tuple(row[field] for field in field_names)
    if not any(values):
        return None
    if not all(values):
        msg = f"incomplete {prefix} confidence interval in CSV row {row_number} of {source}"
        raise ValueError(msg)
    return TimingEstimate(
        median_ns=_parse_float(values[0], field_names[0], row_number=row_number, source=source),
        ci_lower_ns=_parse_float(values[1], field_names[1], row_number=row_number, source=source),
        ci_upper_ns=_parse_float(values[2], field_names[2], row_number=row_number, source=source),
    )


def _parse_rows(csv_payload: bytes, *, source: str) -> tuple[PerformanceRow, ...]:
    try:
        text = csv_payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        msg = f"release-performance CSV is not valid UTF-8 in {source}: {exc}"
        raise ValueError(msg) from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames != list(CSV_COLUMNS):
        msg = f"unsupported release-performance CSV columns in {source}: {reader.fieldnames!r}"
        raise ValueError(msg)
    rows: list[PerformanceRow] = []
    for row_number, raw in enumerate(reader, start=2):
        if None in raw or any(value is None for value in raw.values()):
            msg = f"malformed release-performance CSV row {row_number} in {source}"
            raise ValueError(msg)
        row = cast("dict[str, str]", raw)
        if row["schema_version"] != str(SCHEMA_VERSION):
            msg = f"unsupported CSV schema version in row {row_number} of {source}: {row['schema_version']!r}"
            raise ValueError(msg)
        coverage = row["coverage_status"]
        if coverage not in COVERAGE_STATES:
            msg = f"unsupported coverage status in row {row_number} of {source}: {coverage!r}"
            raise ValueError(msg)
        rows.append(
            PerformanceRow(
                suite=row["suite"],
                scope=row["scope"],
                benchmark_id=row["benchmark_id"],
                group=row["group"],
                benchmark=row["benchmark"],
                baseline_benchmark=row["baseline_benchmark"],
                coverage_status=coverage,
                coverage_note=row["coverage_note"],
                baseline=_parse_timing(row, "baseline", row_number=row_number, source=source),
                current=_parse_timing(row, "current", row_number=row_number, source=source),
                baseline_nalgebra=_parse_timing(row, "baseline_nalgebra", row_number=row_number, source=source),
                baseline_faer=_parse_timing(row, "baseline_faer", row_number=row_number, source=source),
            )
        )
    return tuple(rows)


def _parse_provenance(payload: bytes, *, source: str) -> tuple[ArtifactContext, str, int, tuple[str, ...]]:
    try:
        raw = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        msg = f"malformed release-performance provenance JSON in {source}: {exc}"
        raise ValueError(msg) from exc
    if not isinstance(raw, dict):
        msg = f"release-performance provenance must be a JSON object in {source}"
        raise TypeError(msg)
    data = cast("dict[str, object]", raw)
    schema_version = data.get("schema_version")
    if type(schema_version) is not int or schema_version != SCHEMA_VERSION:
        msg = f"unsupported release-performance provenance schema in {source}: {data.get('schema_version')!r}"
        raise ValueError(msg)

    release_data = _parse_required_object(data, "release", source=source)
    report_data = _parse_required_object(data, "report", source=source)
    source_data = _parse_required_object(report_data, "source", source=source)
    csv_data = _parse_required_object(data, "csv", source=source)
    benchmark_provenance = _parse_required_object(data, "benchmark_provenance", source=source)

    release = ReleasePair(
        current=_parse_required_string(release_data, "current", source=source),
        baseline=_parse_required_string(release_data, "baseline", source=source),
    )
    statistic = _parse_required_string(report_data, "statistic", source=source)
    if statistic != "median":
        msg = f"unsupported report statistic in {source}: {statistic!r}"
        raise ValueError(msg)
    context = ArtifactContext(
        release=release,
        statistic="median",
        suite=_parse_required_string(report_data, "suite", source=source),
        scope=_parse_required_string(report_data, "scope", source=source),
        source=ReportSource(
            version=_parse_required_string(source_data, "version", source=source),
            commit=_parse_required_string(source_data, "commit", source=source),
            ref=_parse_required_string(source_data, "ref", source=source),
            revision_timestamp=_parse_required_string(source_data, "revision_timestamp", source=source),
        ),
        benchmark_provenance=benchmark_provenance,
    )
    expected_current = f"v{context.source.version.removeprefix('v')}"
    if context.release.current != expected_current:
        msg = f"release current {context.release.current!r} does not match report source version {context.source.version!r} in {source}"
        raise ValueError(msg)

    digest = _parse_required_string(csv_data, "sha256", source=source)
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        msg = f"invalid csv.sha256 in {source}"
        raise ValueError(msg)
    row_count = csv_data.get("row_count")
    if isinstance(row_count, bool) or not isinstance(row_count, int) or row_count <= 0:
        msg = f"invalid csv.row_count in {source}: {row_count!r}"
        raise ValueError(msg)
    columns = csv_data.get("columns")
    if columns != list(CSV_COLUMNS):
        msg = f"provenance CSV columns do not match schema {SCHEMA_VERSION} in {source}"
        raise ValueError(msg)
    return context, digest, row_count, tuple(cast("list[str]", columns))


def load_bundle_bytes(csv_payload: bytes, provenance_payload: bytes, *, source: str) -> PerformanceBundle:
    """Parse and validate an artifact pair before it reaches report rendering."""
    context, expected_digest, expected_count, _columns = _parse_provenance(provenance_payload, source=source)
    observed_digest = hashlib.sha256(csv_payload).hexdigest()
    if observed_digest != expected_digest:
        msg = f"CSV digest mismatch in {source}: expected {expected_digest}, got {observed_digest}"
        raise ValueError(msg)
    rows = _parse_rows(csv_payload, source=source)
    if len(rows) != expected_count:
        msg = f"CSV row count mismatch in {source}: expected {expected_count}, got {len(rows)}"
        raise ValueError(msg)
    return PerformanceBundle(context=context, rows=rows)


def load_bundle(paths: ArtifactPaths) -> PerformanceBundle:
    """Load an adjacent CSV/provenance pair and fail closed on partial publication."""
    missing = [str(path) for path in (paths.csv, paths.provenance) if not path.is_file()]
    if missing:
        msg = f"release-performance artifact pair is incomplete; missing: {', '.join(missing)}"
        raise FileNotFoundError(msg)
    return load_bundle_bytes(
        paths.csv.read_bytes(),
        paths.provenance.read_bytes(),
        source=f"{paths.csv} and {paths.provenance}",
    )


def _stage_payload(path: Path, payload: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile("wb", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as tmp:
            staged = Path(tmp.name)
            tmp.write(payload)
            tmp.flush()
            os.fsync(tmp.fileno())
    except BaseException:
        if staged is not None:
            staged.unlink(missing_ok=True)
        raise
    if staged is None:
        msg = "temporary artifact staging completed without a path"
        raise AssertionError(msg)
    return staged


def _atomic_restore(path: Path, payload: bytes | None) -> None:
    if payload is None:
        path.unlink(missing_ok=True)
        return
    staged = _stage_payload(path, payload)
    try:
        _replace_path(staged, path)
    finally:
        staged.unlink(missing_ok=True)


def _publish_payloads(paths: ArtifactPaths, csv_payload: bytes, provenance_payload: bytes) -> None:
    staged_paths: list[Path] = []
    try:
        staged_csv = _stage_payload(paths.csv, csv_payload)
        staged_paths.append(staged_csv)
        staged_provenance = _stage_payload(paths.provenance, provenance_payload)
        staged_paths.append(staged_provenance)
        previous_csv = paths.csv.read_bytes() if paths.csv.is_file() else None
        previous_provenance = paths.provenance.read_bytes() if paths.provenance.is_file() else None
        try:
            _replace_path(staged_csv, paths.csv)
            _replace_path(staged_provenance, paths.provenance)
            load_bundle(paths)
        except BaseException as publication_error:
            rollback_errors = _restore_artifact_pair(paths, previous_csv, previous_provenance)
            if rollback_errors:
                group_message = "release-performance artifact publication and rollback failed"
                raise BaseExceptionGroup(
                    group_message,
                    [publication_error, *rollback_errors],
                ) from None
            raise
    finally:
        for staged in staged_paths:
            staged.unlink(missing_ok=True)


def _replace_path(source: Path, destination: Path) -> None:
    """Atomically replace *destination* with a staged file."""
    source.replace(destination)


def _restore_artifact_pair(
    paths: ArtifactPaths,
    previous_csv: bytes | None,
    previous_provenance: bytes | None,
) -> tuple[BaseException, ...]:
    """Attempt both artifact restorations and return every rollback failure."""
    errors: list[BaseException] = []
    for path, payload in ((paths.csv, previous_csv), (paths.provenance, previous_provenance)):
        try:
            _atomic_restore(path, payload)
        except BaseException as exc:  # noqa: BLE001
            errors.append(exc)
    return tuple(errors)


@contextmanager
def publish_bundle(paths: ArtifactPaths, bundle: PerformanceBundle) -> Iterator[None]:
    """Publish a validated pair and roll it back if a downstream promotion fails."""
    csv_payload, provenance_payload = serialize_bundle(bundle)
    previous_csv = paths.csv.read_bytes() if paths.csv.is_file() else None
    previous_provenance = paths.provenance.read_bytes() if paths.provenance.is_file() else None
    _publish_payloads(paths, csv_payload, provenance_payload)
    try:
        yield
    except BaseException as downstream_error:
        rollback_errors = _restore_artifact_pair(paths, previous_csv, previous_provenance)
        if rollback_errors:
            group_message = "release-performance downstream publication and rollback failed"
            raise BaseExceptionGroup(
                group_message,
                [downstream_error, *rollback_errors],
            ) from None
        raise


def write_bundle(paths: ArtifactPaths, bundle: PerformanceBundle) -> None:
    """Atomically publish a complete validated artifact pair."""
    csv_payload, provenance_payload = serialize_bundle(bundle)
    _publish_payloads(paths, csv_payload, provenance_payload)
