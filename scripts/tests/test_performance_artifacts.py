"""Tests for durable release-performance report artifacts."""

import csv
import hashlib
import io
import json
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import cast

import pytest

import performance_artifacts
from performance_artifacts import (
    ArtifactContext,
    ArtifactPaths,
    PerformanceBundle,
    PerformanceRow,
    ReleasePair,
    ReportSource,
    TimingEstimate,
    load_bundle,
    load_bundle_bytes,
    publish_bundle,
    serialize_bundle,
    write_bundle,
)


def _timing(value: float) -> TimingEstimate:
    return TimingEstimate(
        median_ns=value,
        ci_lower_ns=value * 0.9,
        ci_upper_ns=value * 1.1,
    )


@pytest.mark.parametrize(
    "case",
    [
        (
            "v0.4.5",
            "v0.4.4",
            False,
            "pre-rational-input",
            "pre-rational-input",
            "la_stack_pre_rational_input_api",
            "excluded",
        ),
        (
            "v0.4.4",
            "v0.4.3",
            False,
            "pre-rational-input",
            "legacy-v0.4.3",
            "la_stack_v0_4_3_api",
            "excluded",
        ),
        (
            "v0.4.5",
            "v0.4.5",
            True,
            "pre-rational-input",
            "pre-rational-input",
            "la_stack_pre_rational_input_api",
            "current-only",
        ),
        (
            "v0.4.6-rc.1",
            "v0.4.2",
            True,
            "rational-input",
            "legacy-v0.4.3",
            "la_stack_v0_4_3_api",
            "current-only",
        ),
        (
            "v1.0.0-alpha.1",
            "v0.3.9-beta.2",
            True,
            "rational-input",
            "legacy-v0.4.3",
            "la_stack_v0_4_3_api",
            "current-only",
        ),
        (
            "v0.4.6",
            "v0.4.6-rc.1",
            True,
            "rational-input",
            "rational-input",
            "none",
            "comparable",
        ),
    ],
)
def test_shared_harness_compatibility_uses_literal_release_capability_oracle(
    case: tuple[str, str, bool, str, str, str, str],
) -> None:
    current, baseline, harness_rational, expected_current, expected_baseline, adapter, coverage = case
    resolved = performance_artifacts.resolve_shared_harness_compatibility(
        current=current,
        baseline=baseline,
        shared_harness_rational_inputs=harness_rational,
    )

    assert resolved.current == expected_current
    assert resolved.baseline == expected_baseline
    assert resolved.baseline_api_compatibility == adapter
    assert resolved.rational_input_coverage == coverage


def _context(
    *,
    current: str = "v0.4.4",
    baseline: str = "v0.4.3",
    shared_harness_rational_inputs: bool = True,
) -> ArtifactContext:
    if baseline in {"v0.4.2", "v0.4.3"}:
        compatibility = "la_stack_v0_4_3_api"
    elif baseline in {"v0.4.4", "v0.4.5"}:
        compatibility = "la_stack_pre_rational_input_api"
    else:
        compatibility = "none"
    return ArtifactContext(
        release=ReleasePair(current=current, baseline=baseline),
        statistic="median",
        suite="exact",
        scope="release-signal",
        source=ReportSource(
            version=current.removeprefix("v"),
            commit="abc1234",
            ref="HEAD",
            revision_timestamp="2026-08-04 12:00:00 UTC",
        ),
        benchmark_provenance={
            "baseline": baseline,
            "current": current,
            "criterion": {
                "baseline_command": ["just", "bench-save-baseline", baseline],
                "criterion_version": "0.7.0",
                "current_command": ["just", "bench-latest"],
                "sample": "new",
                "scope": "release-signal",
                "statistic": "median",
                "suite": "exact",
            },
            "measurement": {
                "baseline_api_compatibility": compatibility,
                "baseline_commit": "def5678",
                "baseline_git_clean": True,
                "baseline_source_state_sha256": "d" * 64,
                "cargo_lock_sha256": "b" * 64,
                "cpu": "test-cpu",
                "current_commit": "abc1234",
                "current_git_clean": True,
                "current_source_state_sha256": "c" * 64,
                "harness_sha256": "a" * 64,
                "os": "TestOS 1 x86_64",
                "rustc": "rustc 1.97.1",
                "status": "recorded",
            },
            "mode": "shared-current-harness",
            "publication": {
                "cargo_lock_sha256": "b" * 64,
                "commit": "abc1234",
                "correctness_gate": "passed",
                "cpu": "test-cpu",
                "git_clean": True,
                "harness_sha256": "a" * 64,
                "os": "TestOS 1 x86_64",
                "rustc": "rustc 1.97.1",
                "source_state_sha256": "c" * 64,
            },
            "schema": 2,
            "validation": {
                "baseline_api_compatibility": compatibility,
                "baseline_commit": "def5678",
                "baseline_git_clean": True,
                "baseline_revision": "passed",
                "baseline_source_state_sha256": "d" * 64,
                "command": ["just", "test-bench-inputs"],
                "current_commit": "abc1234",
                "current_git_clean": True,
                "current_revision": "passed",
                "current_source_state_sha256": "c" * 64,
                "harness": "shared-current",
                "shared_harness_rational_inputs": shared_harness_rational_inputs,
            },
        },
    )


def _bundle(*, current_value: float = 8.0) -> PerformanceBundle:
    rows = (
        PerformanceRow(
            suite="exact",
            scope="release-signal",
            benchmark_id="exact_d2/det_exact",
            group="exact_d2",
            benchmark="det_exact",
            baseline_benchmark="det_exact",
            coverage_status="comparable",
            coverage_note="",
            baseline=_timing(10.0),
            current=_timing(current_value),
        ),
        PerformanceRow(
            suite="exact",
            scope="release-signal",
            benchmark_id="exact_d3/current_only",
            group="exact_d3",
            benchmark="current_only",
            baseline_benchmark="current_only",
            coverage_status="current-only",
            coverage_note="The baseline does not expose this correctness-compatible row.",
            baseline=None,
            current=_timing(12.0),
        ),
        PerformanceRow(
            suite="exact",
            scope="release-signal",
            benchmark_id="exact_d4/baseline_only",
            group="exact_d4",
            benchmark="baseline_only",
            baseline_benchmark="baseline_only",
            coverage_status="baseline-only",
            coverage_note="The current run does not expose this row.",
            baseline=_timing(15.0),
            current=None,
        ),
    )
    return PerformanceBundle(context=_context(), rows=rows)


def _replace_csv_and_digest(csv_payload: bytes, provenance_payload: bytes) -> tuple[bytes, bytes]:
    provenance = json.loads(provenance_payload)
    provenance["csv"]["sha256"] = hashlib.sha256(csv_payload).hexdigest()
    return csv_payload, (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode()


def _simulate_promotion_failure(paths: ArtifactPaths) -> None:
    with publish_bundle(paths, _bundle(current_value=7.0)):
        msg = "simulated promotion failure"
        raise RuntimeError(msg)


def test_artifact_round_trip_preserves_comparable_and_one_sided_rows() -> None:
    bundle = _bundle()

    csv_payload, provenance_payload = serialize_bundle(bundle)
    parsed = load_bundle_bytes(csv_payload, provenance_payload, source="round-trip fixture")

    assert parsed == PerformanceBundle(context=bundle.context, rows=bundle.sorted_rows)
    assert [row.coverage_status for row in parsed.rows] == ["comparable", "current-only", "baseline-only"]
    assert csv_payload.endswith(b"\n")
    assert provenance_payload.endswith(b"\n")


def test_bundle_rejects_rational_rows_when_historical_harness_did_not_emit_them() -> None:
    rational_row = PerformanceRow(
        suite="exact",
        scope="release-signal",
        benchmark_id="rational_input_d2/det_row_cleared_bareiss",
        group="rational_input_d2",
        benchmark="det_row_cleared_bareiss",
        baseline_benchmark="det_row_cleared_bareiss",
        coverage_status="current-only",
        coverage_note="Baseline predates rational inputs.",
        baseline=None,
        current=_timing(12.0),
    )
    context = _context(
        current="v0.4.5",
        baseline="v0.4.4",
        shared_harness_rational_inputs=False,
    )

    with pytest.raises(ValueError, match="must be excluded"):
        PerformanceBundle(context=context, rows=(rational_row,))


def test_artifact_round_trip_allows_same_version_local_comparison() -> None:
    original = _bundle()
    bundle = PerformanceBundle(
        context=_context(current="v0.4.4", baseline="v0.4.4"),
        rows=original.rows,
    )

    csv_payload, provenance_payload = serialize_bundle(bundle)
    parsed = load_bundle_bytes(csv_payload, provenance_payload, source="same-version local fixture")

    assert parsed.context.release == ReleasePair(current="v0.4.4", baseline="v0.4.4")
    assert parsed.rows == bundle.sorted_rows


def test_artifact_serialization_is_independent_of_input_row_order() -> None:
    bundle = _bundle()
    reordered = PerformanceBundle(context=bundle.context, rows=tuple(reversed(bundle.rows)))

    assert serialize_bundle(reordered) == serialize_bundle(bundle)


def test_artifact_loader_rejects_csv_digest_mismatch() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    malformed = csv_payload.replace(b"exact_d2/det_exact", b"exact_d2/det_changed", 1)

    with pytest.raises(ValueError, match="CSV digest mismatch"):
        load_bundle_bytes(malformed, provenance_payload, source="digest mismatch fixture")


def test_artifact_loader_rejects_csv_row_count_mismatch() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["csv"]["row_count"] += 1
    malformed = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode()

    with pytest.raises(ValueError, match="CSV row count mismatch"):
        load_bundle_bytes(csv_payload, malformed, source="row count mismatch fixture")


def test_artifact_loader_rejects_incomplete_confidence_interval() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    reader = csv.DictReader(io.StringIO(csv_payload.decode(), newline=""))
    rows = list(reader)
    rows[0]["current_ci_upper_ns"] = ""
    output = io.StringIO(newline="")
    fieldnames = reader.fieldnames
    assert fieldnames is not None
    writer = csv.DictWriter(output, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    malformed = output.getvalue().encode()
    malformed, provenance_payload = _replace_csv_and_digest(malformed, provenance_payload)

    with pytest.raises(ValueError, match="incomplete current confidence interval"):
        load_bundle_bytes(malformed, provenance_payload, source="incomplete interval fixture")


def test_artifact_loader_rejects_mismatched_release_pair() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["release"]["baseline"] = "v0.4.2"
    malformed = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode()

    with pytest.raises(ValueError, match="benchmark provenance baseline"):
        load_bundle_bytes(csv_payload, malformed, source="mismatched release fixture")


def test_artifact_loader_rejects_malformed_csv_schema() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    malformed = csv_payload.replace(b"benchmark_id", b"benchmark_key", 1)
    malformed, provenance_payload = _replace_csv_and_digest(malformed, provenance_payload)

    with pytest.raises(ValueError, match="unsupported release-performance CSV columns"):
        load_bundle_bytes(malformed, provenance_payload, source="malformed CSV fixture")


def test_artifact_loader_rejects_malformed_provenance_json() -> None:
    csv_payload, _provenance_payload = serialize_bundle(_bundle())

    with pytest.raises(ValueError, match="malformed release-performance provenance JSON"):
        load_bundle_bytes(csv_payload, b"{not-json}\n", source="malformed provenance fixture")


@pytest.mark.parametrize("invalid_schema", [True, 1.0])
def test_artifact_loader_rejects_non_integer_outer_schema(invalid_schema: object) -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["schema_version"] = invalid_schema

    with pytest.raises(ValueError, match="unsupported release-performance provenance schema"):
        load_bundle_bytes(
            csv_payload,
            (json.dumps(provenance) + "\n").encode(),
            source="invalid outer schema fixture",
        )


@pytest.mark.parametrize("invalid_schema", [True, 2.0])
def test_artifact_loader_rejects_non_integer_benchmark_schema(invalid_schema: object) -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["benchmark_provenance"]["schema"] = invalid_schema

    with pytest.raises(ValueError, match="benchmark provenance schema must be 2"):
        load_bundle_bytes(
            csv_payload,
            (json.dumps(provenance) + "\n").encode(),
            source="invalid benchmark schema fixture",
        )


def test_artifact_loader_rejects_incomplete_nested_provenance() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    del provenance["benchmark_provenance"]["criterion"]["criterion_version"]

    with pytest.raises(ValueError, match=r"criterion\.criterion_version"):
        load_bundle_bytes(
            csv_payload,
            (json.dumps(provenance) + "\n").encode(),
            source="incomplete nested provenance fixture",
        )


def test_artifact_loader_defaults_fields_absent_from_legacy_schema2_artifacts() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    benchmark_provenance = provenance["benchmark_provenance"]
    assert benchmark_provenance["schema"] == 2
    del benchmark_provenance["current"]
    del benchmark_provenance["validation"]["shared_harness_rational_inputs"]

    loaded = load_bundle_bytes(
        csv_payload,
        (json.dumps(provenance) + "\n").encode(),
        source="legacy schema-2 artifact fixture",
    )

    loaded_provenance = loaded.context.benchmark_provenance
    assert loaded_provenance["current"] == loaded.context.release.current
    validation = loaded_provenance["validation"]
    assert isinstance(validation, Mapping)
    assert validation["shared_harness_rational_inputs"] is False


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("current", False, r"root\.current"),
        ("shared_harness_rational_inputs", "false", "must be a boolean"),
    ],
)
def test_artifact_loader_rejects_invalid_present_legacy_default_fields(
    field: str,
    value: object,
    message: str,
) -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    benchmark_provenance = provenance["benchmark_provenance"]
    if field == "current":
        benchmark_provenance[field] = value
    else:
        benchmark_provenance["validation"][field] = value

    with pytest.raises((TypeError, ValueError), match=message):
        load_bundle_bytes(
            csv_payload,
            (json.dumps(provenance) + "\n").encode(),
            source="invalid legacy-default field fixture",
        )


def test_artifact_loader_rejects_contradictory_current_revision() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["benchmark_provenance"]["validation"]["current_commit"] = "different-commit"

    with pytest.raises(ValueError, match=r"validation\.current_commit.*publication\.commit"):
        load_bundle_bytes(
            csv_payload,
            (json.dumps(provenance) + "\n").encode(),
            source="contradictory current revision fixture",
        )


def test_artifact_loader_rejects_recorded_measurement_without_cpu_model() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["benchmark_provenance"]["measurement"]["cpu"] = "unavailable"
    provenance["benchmark_provenance"]["publication"]["cpu"] = "unavailable"

    with pytest.raises(ValueError, match="requires an identified CPU model"):
        load_bundle_bytes(
            csv_payload,
            (json.dumps(provenance) + "\n").encode(),
            source="unidentified CPU fixture",
        )


def test_bundle_rejects_duplicate_benchmark_keys() -> None:
    bundle = _bundle()

    with pytest.raises(ValueError, match="duplicate benchmark key"):
        PerformanceBundle(context=bundle.context, rows=(bundle.rows[0], bundle.rows[0]))


@pytest.mark.parametrize("value", [0.0, -1.0, float("inf"), float("nan")])
def test_timing_rejects_non_positive_or_non_finite_values(value: float) -> None:
    with pytest.raises(ValueError, match="must be finite and positive"):
        TimingEstimate(median_ns=value, ci_lower_ns=1.0, ci_upper_ns=2.0)


def test_timing_accepts_ordered_bootstrap_interval_that_excludes_point_estimate() -> None:
    estimate = TimingEstimate(median_ns=12.0, ci_lower_ns=9.0, ci_upper_ns=11.0)

    assert estimate.median_ns == 12.0


def test_timing_rejects_reversed_interval() -> None:
    with pytest.raises(ValueError, match="confidence interval must be ordered"):
        TimingEstimate(median_ns=10.0, ci_lower_ns=11.0, ci_upper_ns=9.0)


def test_artifact_context_freezes_nested_provenance() -> None:
    context = _context()

    with pytest.raises(TypeError):
        cast("dict[str, object]", context.benchmark_provenance)["schema"] = 3

    criterion = context.benchmark_provenance["criterion"]
    assert isinstance(criterion, MappingProxyType)


def test_artifact_loader_rejects_mode_measurement_contradiction() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["benchmark_provenance"]["mode"] = "historical-assets"

    with pytest.raises(ValueError, match="requires shared-current-harness mode"):
        load_bundle_bytes(csv_payload, (json.dumps(provenance) + "\n").encode(), source="contradictory mode fixture")


def test_artifact_loader_rejects_arbitrary_compatibility_adapter() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    provenance["benchmark_provenance"]["validation"]["baseline_api_compatibility"] = "custom-adapter"

    with pytest.raises(ValueError, match="baseline_api_compatibility must be"):
        load_bundle_bytes(csv_payload, (json.dumps(provenance) + "\n").encode(), source="invalid compatibility fixture")


def test_artifact_loader_accepts_pre_rational_adapter_for_v045() -> None:
    bundle = _bundle()
    context = _context(current="v0.4.6", baseline="v0.4.5")
    csv_payload, provenance_payload = serialize_bundle(PerformanceBundle(context=context, rows=bundle.rows))

    loaded = load_bundle_bytes(
        csv_payload,
        provenance_payload,
        source="pre-rational compatibility fixture",
    )

    validation = loaded.context.benchmark_provenance["validation"]
    assert isinstance(validation, Mapping)
    assert validation["baseline_api_compatibility"] == "la_stack_pre_rational_input_api"


def test_artifact_loader_requires_recorded_measurement_compatibility() -> None:
    csv_payload, provenance_payload = serialize_bundle(_bundle())
    provenance = json.loads(provenance_payload)
    del provenance["benchmark_provenance"]["measurement"]["baseline_api_compatibility"]

    with pytest.raises(ValueError, match=r"measurement\.baseline_api_compatibility"):
        load_bundle_bytes(csv_payload, (json.dumps(provenance) + "\n").encode(), source="missing compatibility fixture")


def test_artifact_loader_fails_closed_on_partial_pair(tmp_path: Path) -> None:
    paths = ArtifactPaths(
        csv=tmp_path / "performance.csv",
        provenance=tmp_path / "performance.provenance.json",
    )
    paths.csv.write_text("partial\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="artifact pair is incomplete"):
        load_bundle(paths)


def test_stage_payload_removes_temporary_file_when_fsync_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "performance.csv"

    def fail_fsync(_descriptor: int) -> None:
        msg = "simulated fsync failure"
        raise OSError(msg)

    monkeypatch.setattr(performance_artifacts.os, "fsync", fail_fsync)

    with pytest.raises(OSError, match="simulated fsync failure"):
        performance_artifacts._stage_payload(target, b"payload")

    assert not list(tmp_path.glob(".performance.csv.*.tmp"))


def test_failed_second_stage_removes_first_temporary_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = ArtifactPaths(
        csv=tmp_path / "performance.csv",
        provenance=tmp_path / "performance.provenance.json",
    )
    real_stage_payload = performance_artifacts._stage_payload

    def fail_provenance_stage(path: Path, payload: bytes) -> Path:
        if path == paths.provenance:
            msg = "simulated provenance staging failure"
            raise OSError(msg)
        return real_stage_payload(path, payload)

    monkeypatch.setattr(performance_artifacts, "_stage_payload", fail_provenance_stage)

    with pytest.raises(OSError, match="simulated provenance staging failure"):
        write_bundle(paths, _bundle())

    assert not list(tmp_path.glob(".performance.csv.*.tmp"))
    assert not paths.csv.exists()
    assert not paths.provenance.exists()


def test_failed_second_replace_restores_prior_valid_pair(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = ArtifactPaths(
        csv=tmp_path / "performance.csv",
        provenance=tmp_path / "performance.provenance.json",
    )
    original = _bundle(current_value=8.0)
    write_bundle(paths, original)
    real_replace = performance_artifacts._replace_path
    failed = False

    def fail_provenance_once(source: Path, destination: Path) -> None:
        nonlocal failed
        if Path(destination) == paths.provenance and not failed:
            failed = True
            msg = "simulated provenance replace failure"
            raise OSError(msg)
        real_replace(source, destination)

    monkeypatch.setattr(performance_artifacts, "_replace_path", fail_provenance_once)

    with pytest.raises(OSError, match="simulated provenance replace failure"):
        write_bundle(paths, _bundle(current_value=7.0))

    assert load_bundle(paths) == PerformanceBundle(context=original.context, rows=original.sorted_rows)


def test_downstream_promotion_failure_rolls_back_artifact_pair(tmp_path: Path) -> None:
    paths = ArtifactPaths(
        csv=tmp_path / "performance.csv",
        provenance=tmp_path / "performance.provenance.json",
    )
    original = _bundle(current_value=8.0)
    write_bundle(paths, original)

    with pytest.raises(RuntimeError, match="simulated promotion failure"):
        _simulate_promotion_failure(paths)

    assert load_bundle(paths) == PerformanceBundle(context=original.context, rows=original.sorted_rows)


def test_rollback_attempts_both_artifacts_when_first_restoration_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = ArtifactPaths(
        csv=tmp_path / "performance.csv",
        provenance=tmp_path / "performance.provenance.json",
    )
    write_bundle(paths, _bundle(current_value=8.0))
    real_restore = performance_artifacts._atomic_restore
    restored: list[Path] = []

    def fail_csv_restore(path: Path, payload: bytes | None) -> None:
        restored.append(path)
        if path == paths.csv:
            msg = "simulated CSV restoration failure"
            raise OSError(msg)
        real_restore(path, payload)

    monkeypatch.setattr(performance_artifacts, "_atomic_restore", fail_csv_restore)

    with pytest.raises(BaseExceptionGroup, match="downstream publication and rollback failed") as raised:
        _simulate_promotion_failure(paths)

    assert restored == [paths.csv, paths.provenance]
    assert any(isinstance(error, RuntimeError) for error in raised.value.exceptions)
    assert any(isinstance(error, OSError) for error in raised.value.exceptions)
