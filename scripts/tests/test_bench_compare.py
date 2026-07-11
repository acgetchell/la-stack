"""Tests for exact-arithmetic benchmark comparison reports."""

from __future__ import annotations

import json
import re
import subprocess
from typing import TYPE_CHECKING, cast

import pytest

import bench_compare

if TYPE_CHECKING:
    from pathlib import Path

_OVERFLOWING_TIMING = 10**400


def _write_estimates(
    path: Path,
    stat: str,
    median: float,
    *,
    lower: float | None = None,
    upper: float | None = None,
) -> None:
    """Write a minimal Criterion estimates.json."""
    lower = median * 0.9 if lower is None else lower
    upper = median * 1.1 if upper is None else upper
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                stat: {
                    "point_estimate": median,
                    "confidence_interval": {"lower_bound": lower, "upper_bound": upper},
                }
            }
        ),
        encoding="utf-8",
    )


def _write_harness_provenance(
    criterion_dir: Path,
    *,
    sha256: str = "a" * 64,
    baseline: str = "v0.4.3",
) -> None:
    """Write valid shared-harness provenance metadata."""
    criterion_dir.mkdir(parents=True, exist_ok=True)
    (criterion_dir / ".la-stack-benchmark-harness.json").write_text(
        json.dumps(
            {
                "schema": 1,
                "mode": "shared-current-harness",
                "sha256": sha256,
                "baseline": baseline,
            }
        ),
        encoding="utf-8",
    )


def _schema2_provenance_data() -> dict[str, object]:
    environment: dict[str, object] = {
        "cargo_lock_sha256": "b" * 64,
        "commit": "current-commit",
        "correctness_gate": "passed",
        "cpu": "test-cpu",
        "git_clean": False,
        "harness_sha256": "a" * 64,
        "os": "TestOS 1 x86_64",
        "rustc": "rustc 1.88.0",
        "source_state_sha256": "c" * 64,
    }
    return {
        "baseline": "v0.4.3",
        "criterion": {
            "baseline_command": ["just", "bench-save-baseline", "v0.4.3"],
            "criterion_version": "0.7.0",
            "current_command": ["just", "bench-latest"],
            "sample": "new",
            "scope": "release-signal",
            "statistic": "median",
            "suite": "all",
        },
        "measurement": {
            "baseline_commit": "baseline-commit",
            "baseline_git_clean": False,
            "baseline_source_state_sha256": "d" * 64,
            "cargo_lock_sha256": "b" * 64,
            "cpu": "test-cpu",
            "current_commit": "current-commit",
            "current_git_clean": False,
            "current_source_state_sha256": "c" * 64,
            "harness_sha256": "a" * 64,
            "os": "TestOS 1 x86_64",
            "rustc": "rustc 1.88.0",
            "status": "recorded",
        },
        "mode": "shared-current-harness",
        "publication": environment,
        "schema": 2,
        "validation": {
            "baseline_api_compatibility": "la_stack_v0_4_3_api",
            "baseline_commit": "baseline-commit",
            "baseline_git_clean": False,
            "baseline_revision": "passed",
            "baseline_source_state_sha256": "d" * 64,
            "command": ["just", "test-bench-inputs"],
            "current_commit": "current-commit",
            "current_git_clean": False,
            "current_revision": "passed",
            "current_source_state_sha256": "c" * 64,
            "harness": "shared-current",
        },
    }


def _read_harness_provenance(  # noqa: PLR0913
    criterion_dir: Path,
    *,
    baseline: str = "v0.4.3",
    suite: bench_compare.BenchmarkSuite = "all",
    scope: bench_compare.ComparisonScope = "release-signal",
    stat: bench_compare.Statistic = "median",
    sample: str = "new",
) -> bench_compare.HarnessProvenance | None:
    return bench_compare._read_harness_provenance(
        criterion_dir,
        expected_baseline=baseline,
        expected=bench_compare.CriterionSelection(
            suite=suite,
            scope=scope,
            statistic=stat,
            sample=sample,
        ),
    )


def _build_criterion_tree(criterion_dir: Path, stat: str = "median") -> None:
    """Create a fake Criterion directory with exact benchmark results."""
    for d, det, det_exact in [(2, 1.0, 4000.0), (3, 5.0, 21000.0)]:
        group = criterion_dir / f"exact_d{d}"
        _write_estimates(group / "det" / "new" / "estimates.json", stat, det)
        _write_estimates(group / "det_exact" / "new" / "estimates.json", stat, det_exact)
        _write_estimates(group / "det_exact_f64_result" / "new" / "estimates.json", stat, det_exact * 1.1)
        _write_estimates(group / "det_exact_rounded_f64" / "new" / "estimates.json", stat, det_exact * 1.2)
        _write_estimates(group / "solve_exact_f64_result" / "new" / "estimates.json", stat, det_exact * 1.3)
        _write_estimates(group / "solve_exact_rounded_f64" / "new" / "estimates.json", stat, det_exact * 1.4)

    ns_group = criterion_dir / "exact_near_singular_3x3"
    _write_estimates(ns_group / "det_sign_exact" / "new" / "estimates.json", stat, 12000.0)
    _write_estimates(ns_group / "solve_exact_f64_result" / "new" / "estimates.json", stat, 51000.0)
    _write_estimates(ns_group / "solve_exact_rounded_f64" / "new" / "estimates.json", stat, 52000.0)

    random_group = criterion_dir / "exact_random_corpus_d3"
    _write_estimates(random_group / "det_exact" / "new" / "estimates.json", stat, 33000.0)
    _write_estimates(random_group / "solve_exact_f64_result" / "new" / "estimates.json", stat, 54000.0)
    _write_estimates(random_group / "solve_exact_rounded_f64" / "new" / "estimates.json", stat, 55000.0)


def _build_vs_linalg_tree(criterion_dir: Path, stat: str = "median") -> None:
    """Create fake Criterion vs_linalg data with latest la-stack and baseline peer rows."""
    group = criterion_dir / "d2"

    _write_estimates(group / "la_stack_lu_solve" / "new" / "estimates.json", stat, 10.0)
    _write_estimates(group / "la_stack_lu_solve" / "last" / "estimates.json", stat, 20.0)
    _write_estimates(group / "nalgebra_lu_solve" / "last" / "estimates.json", stat, 30.0)
    _write_estimates(group / "faer_lu_solve" / "last" / "estimates.json", stat, 40.0)

    _write_estimates(group / "la_stack_ldlt_solve" / "new" / "estimates.json", stat, 12.0)
    _write_estimates(group / "la_stack_ldlt_solve" / "last" / "estimates.json", stat, 24.0)
    _write_estimates(group / "nalgebra_cholesky_solve" / "last" / "estimates.json", stat, 36.0)
    _write_estimates(group / "faer_ldlt_solve" / "last" / "estimates.json", stat, 48.0)

    # Present only in the baseline; release-signal comparison should not require
    # latest third-party measurements.
    _write_estimates(group / "nalgebra_lu_solve" / "new" / "estimates.json", stat, 31.0)


# ---------------------------------------------------------------------------
# Unit tests for formatting helpers
# ---------------------------------------------------------------------------


class TestFormatTime:
    def test_nanoseconds(self) -> None:
        assert bench_compare._format_time(42.5) == "42.5 ns"

    def test_microseconds(self) -> None:
        assert bench_compare._format_time(1_500.0) == "1.50 µs"

    def test_milliseconds(self) -> None:
        assert bench_compare._format_time(2_500_000.0) == "2.50 ms"


class TestFormatPct:
    def test_ci_separated_point_change_is_not_bold(self) -> None:
        assert bench_compare._format_pct(-10.0, "improvement") == "-10.0%"

    def test_inconclusive_point_improvement_is_not_bold(self) -> None:
        assert bench_compare._format_pct(-10.0, "inconclusive") == "-10.0%"

    def test_small_change(self) -> None:
        assert bench_compare._format_pct(0.5, "inconclusive") == "+0.5%"


# ---------------------------------------------------------------------------
# Group heading
# ---------------------------------------------------------------------------


class TestGroupHeading:
    def test_dimension_group(self) -> None:
        assert bench_compare._group_heading("exact_d3") == "D=3"

    def test_near_singular(self) -> None:
        assert bench_compare._group_heading("exact_near_singular_3x3") == "Near-singular 3x3"

    def test_random_corpus(self) -> None:
        assert bench_compare._group_heading("exact_random_corpus_d4") == "Random corpus D=4"

    def test_large_entries(self) -> None:
        assert bench_compare._group_heading("exact_large_entries_3x3") == "Large entries 3x3"

    def test_hilbert_4x4(self) -> None:
        assert bench_compare._group_heading("exact_hilbert_4x4") == "Hilbert 4x4"

    def test_hilbert_5x5(self) -> None:
        assert bench_compare._group_heading("exact_hilbert_5x5") == "Hilbert 5x5"

    def test_unknown_passthrough(self) -> None:
        assert bench_compare._group_heading("something_else") == "something_else"


def test_exact_registry_only_tracks_supported_direct_determinants() -> None:
    for dimension in (2, 3, 4):
        assert "det_direct" in bench_compare.EXACT_GROUPS[f"exact_d{dimension}"]
    assert "det_direct" not in bench_compare.EXACT_GROUPS["exact_d5"]


# ---------------------------------------------------------------------------
# read_estimate
# ---------------------------------------------------------------------------


def test_read_estimate_success(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    _write_estimates(est, "median", 42.0)
    estimate = bench_compare._read_estimate(est, "median")
    assert estimate.point_ns == 42.0
    assert estimate.ci_lo_ns == pytest.approx(42.0 * 0.9)
    assert estimate.ci_hi_ns == pytest.approx(42.0 * 1.1)


def test_read_estimate_no_ci(tmp_path: Path) -> None:
    """Missing confidence bounds remain unavailable rather than implying certainty."""
    est = tmp_path / "estimates.json"
    est.parent.mkdir(parents=True, exist_ok=True)
    est.write_text(
        json.dumps({"median": {"point_estimate": 99.0}}),
        encoding="utf-8",
    )
    estimate = bench_compare._read_estimate(est, "median")
    assert estimate.point_ns == 99.0
    assert estimate.ci_lo_ns is None
    assert estimate.ci_hi_ns is None


def test_read_estimate_missing_stat(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    _write_estimates(est, "median", 1.0)
    with pytest.raises(KeyError, match="stat 'mean' not found"):
        bench_compare._read_estimate(est, "mean")


def test_read_estimate_malformed_json_names_file(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    est.write_text("{not json", encoding="utf-8")

    with pytest.raises(
        ValueError,
        match=re.escape(f"malformed Criterion estimates JSON in {est}"),
    ):
        bench_compare._read_estimate(est, "median")


def test_read_estimate_missing_point_estimate_names_field(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    est.write_text(json.dumps({"median": {}}), encoding="utf-8")

    with pytest.raises(KeyError, match="field 'point_estimate' for stat 'median' not found"):
        bench_compare._read_estimate(est, "median")


def test_read_estimate_non_numeric_ci_bound_names_field(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    est.write_text(
        json.dumps(
            {
                "median": {
                    "point_estimate": 1.0,
                    "confidence_interval": {"lower_bound": "fast", "upper_bound": 2.0},
                }
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"field 'lower_bound' for stat 'median'.*not numeric"):
        bench_compare._read_estimate(est, "median")


def test_read_estimate_rejects_numeric_overflow(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    est.write_text(json.dumps({"median": {"point_estimate": _OVERFLOWING_TIMING}}), encoding="utf-8")

    with pytest.raises(ValueError, match=r"field 'point_estimate'.*not numeric") as exc_info:
        bench_compare._read_estimate(est, "median")

    assert isinstance(exc_info.value.__cause__, OverflowError)


def test_read_estimate_rejects_partial_confidence_interval(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    est.write_text(
        json.dumps(
            {
                "median": {
                    "point_estimate": 1.0,
                    "confidence_interval": {"lower_bound": 0.9},
                }
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(KeyError, match="field 'upper_bound'"):
        bench_compare._read_estimate(est, "median")


def test_read_estimate_rejects_reversed_confidence_interval(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    _write_estimates(est, "median", 10.0, lower=12.0, upper=11.0)

    with pytest.raises(ValueError, match=r"lower_bound 12\.0 exceeds upper_bound 11\.0"):
        bench_compare._read_estimate(est, "median")


@pytest.mark.parametrize("point", [float("nan"), float("inf"), -1.0, 0.0])
def test_read_estimate_rejects_invalid_timing(tmp_path: Path, point: float) -> None:
    est = tmp_path / "estimates.json"
    est.write_text(json.dumps({"median": {"point_estimate": point}}), encoding="utf-8")

    with pytest.raises(ValueError, match="must be finite and positive"):
        bench_compare._read_estimate(est, "median")


def test_criterion_estimate_rejects_partial_interval_even_when_constructed_directly() -> None:
    with pytest.raises(ValueError, match="both bounds or neither"):
        bench_compare.CriterionEstimate(point_ns=1.0, ci_lo_ns=0.9, ci_hi_ns=None)


# ---------------------------------------------------------------------------
# collect_results / collect_comparisons
# ---------------------------------------------------------------------------


def test_collect_results(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    results = bench_compare._collect_results(tmp_path, "new", "median")
    assert len(results) == 18  # 6 benches x 2 dims + 3 near-singular + 3 random corpus
    groups = {r.group for r in results}
    assert "exact_d2" in groups
    assert "exact_d3" in groups
    assert "exact_random_corpus_d3" in groups
    assert "exact_near_singular_3x3" in groups


def test_collect_results_empty_dir(tmp_path: Path) -> None:
    results = bench_compare._collect_results(tmp_path, "new", "median")
    assert results == []


def test_collect_comparisons(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    # Save the same data as a baseline
    for d, det, det_exact in [(2, 2.0, 8000.0), (3, 10.0, 42000.0)]:
        group = tmp_path / f"exact_d{d}"
        _write_estimates(group / "det" / "v0.3.0" / "estimates.json", "median", det)
        _write_estimates(group / "det_exact" / "v0.3.0" / "estimates.json", "median", det_exact)
        _write_estimates(group / "det_exact_f64" / "v0.3.0" / "estimates.json", "median", det_exact * 1.1)
        _write_estimates(group / "solve_exact_f64" / "v0.3.0" / "estimates.json", "median", det_exact * 1.3)
    random_group = tmp_path / "exact_random_corpus_d3"
    _write_estimates(random_group / "det_exact" / "v0.3.0" / "estimates.json", "median", 66000.0)
    _write_estimates(random_group / "solve_exact_f64" / "v0.3.0" / "estimates.json", "median", 108000.0)

    collection = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median")
    comparisons = collection.comparisons
    assert len(comparisons) == 15  # fixed dimensions plus the stable random corpus
    assert {c.group for c in comparisons} == {
        "exact_d2",
        "exact_d3",
        "exact_random_corpus_d3",
    }
    for c in comparisons:
        assert c.speedup == pytest.approx(c.baseline_ns / c.current_ns)
    assert {(c.bench, c.baseline_bench) for c in comparisons if c.baseline_bench is not None} == {
        ("det_exact_f64_result", "det_exact_f64"),
        ("det_exact_rounded_f64", "det_exact_f64"),
        ("solve_exact_f64_result", "solve_exact_f64"),
        ("solve_exact_rounded_f64", "solve_exact_f64"),
    }

    all_collection = bench_compare._collect_comparisons(
        tmp_path,
        "v0.3.0",
        "median",
        policy=bench_compare.ComparisonPolicy(scope="all-benches"),
    )
    all_comparisons = all_collection.comparisons
    assert len(all_comparisons) == 15
    random_comparisons = [c for c in all_comparisons if c.group == "exact_random_corpus_d3"]
    assert {(c.bench, c.baseline_bench) for c in random_comparisons} == {
        ("det_exact", None),
        ("solve_exact_f64_result", "solve_exact_f64"),
        ("solve_exact_rounded_f64", "solve_exact_f64"),
    }


def test_collect_comparisons_missing_baseline(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    collection = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median", suite="exact")
    assert collection.comparisons == []
    det_gap = next(gap for gap in collection.gaps if gap.group == "exact_d2" and gap.bench == "det")
    assert not det_gap.missing_current
    assert det_gap.missing_baseline


def test_collect_comparisons_records_each_missing_side_in_registry_order(tmp_path: Path) -> None:
    group = tmp_path / "exact_d2"
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 10.0)
    _write_estimates(group / "det_direct" / "last" / "estimates.json", "median", 20.0)

    collection = bench_compare._collect_comparisons(tmp_path, "last", "median", suite="exact")
    first_three = collection.gaps[:3]
    assert [gap.bench for gap in first_three] == ["det", "det_direct", "det_exact"]
    assert (first_three[0].missing_current, first_three[0].missing_baseline) == (False, True)
    assert (first_three[1].missing_current, first_three[1].missing_baseline) == (True, False)
    assert (first_three[2].missing_current, first_three[2].missing_baseline) == (True, True)


def test_collect_comparisons_reports_wholly_absent_selected_suite(tmp_path: Path) -> None:
    group = tmp_path / "exact_d2"
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 10.0)
    _write_estimates(group / "det" / "last" / "estimates.json", "median", 20.0)

    collection = bench_compare._collect_comparisons(
        tmp_path,
        "last",
        "median",
        suite="all",
        policy=bench_compare.ComparisonPolicy(scope="all-benches"),
    )

    assert any(gap.suite == "vs_linalg" and gap.group == "(entire suite)" and gap.bench == "all selected rows" for gap in collection.gaps)


def test_release_signal_comparison_requires_every_canonical_vs_linalg_dimension(tmp_path: Path) -> None:
    group = tmp_path / "d2"
    _write_estimates(group / "la_stack_lu" / "new" / "estimates.json", "median", 10.0)
    _write_estimates(group / "la_stack_lu" / "last" / "estimates.json", "median", 20.0)

    collection = bench_compare._collect_comparisons(tmp_path, "last", "median", suite="vs_linalg")

    missing_groups = {gap.group for gap in collection.gaps}
    assert {"d3", "d4", "d5", "d8", "d16", "d32", "d64"} <= missing_groups


def test_collect_comparisons_classifies_from_criterion_intervals(tmp_path: Path) -> None:
    group = tmp_path / "exact_d2"
    fixtures = {
        "det": ((100.0, 95.0, 105.0), (80.0, 75.0, 85.0), "improvement"),
        "det_direct": ((100.0, 95.0, 105.0), (120.0, 115.0, 125.0), "regression"),
        "det_exact": ((100.0, 90.0, 110.0), (95.0, 85.0, 105.0), "inconclusive"),
    }
    for bench, (baseline, current, _assessment) in fixtures.items():
        _write_estimates(
            group / bench / "last" / "estimates.json",
            "median",
            baseline[0],
            lower=baseline[1],
            upper=baseline[2],
        )
        _write_estimates(
            group / bench / "new" / "estimates.json",
            "median",
            current[0],
            lower=current[1],
            upper=current[2],
        )

    no_ci_path = group / "det_exact_f64_result" / "new" / "estimates.json"
    no_ci_path.parent.mkdir(parents=True, exist_ok=True)
    no_ci_path.write_text(json.dumps({"median": {"point_estimate": 80.0}}), encoding="utf-8")
    _write_estimates(group / "det_exact_f64_result" / "last" / "estimates.json", "median", 100.0)

    comparisons = bench_compare._collect_comparisons(tmp_path, "last", "median", suite="exact").comparisons
    assessments = {comparison.bench: comparison.assessment for comparison in comparisons}
    assert assessments == {
        "det": "improvement",
        "det_direct": "regression",
        "det_exact": "inconclusive",
        "det_exact_f64_result": "unknown",
    }


@pytest.mark.parametrize("compatibility", [None, "unknown-adapter"])
def test_d8_release_signal_rows_are_explicit_and_report_missing_baselines(
    tmp_path: Path,
    compatibility: str | None,
) -> None:
    group = tmp_path / "d8"
    for bench in bench_compare.VS_LINALG_D8_RELEASE_SIGNAL_BENCHES:
        _write_estimates(group / bench / "new" / "estimates.json", "median", 10.0)

    collection = bench_compare._collect_comparisons(
        tmp_path,
        "last",
        "median",
        suite="vs_linalg",
        policy=bench_compare.ComparisonPolicy(baseline_api_compatibility=compatibility),
    )
    special_gaps = [gap for gap in collection.gaps if gap.bench in bench_compare.VS_LINALG_D8_RELEASE_SIGNAL_BENCHES]
    assert [gap.bench for gap in special_gaps] == bench_compare.VS_LINALG_D8_RELEASE_SIGNAL_BENCHES
    assert all(not gap.missing_current and gap.missing_baseline for gap in special_gaps)


def test_v043_adapter_allows_only_known_unavailable_d8_baselines(tmp_path: Path) -> None:
    group = tmp_path / "d8"
    for bench in bench_compare.VS_LINALG_D8_RELEASE_SIGNAL_BENCHES:
        _write_estimates(group / bench / "new" / "estimates.json", "median", 10.0)
    _write_estimates(
        group / "la_stack_det_from_lu_balanced_range" / "v0.4.3" / "estimates.json",
        "median",
        1.0,
    )

    collection = bench_compare._collect_comparisons(
        tmp_path,
        "v0.4.3",
        "median",
        suite="vs_linalg",
        policy=bench_compare.ComparisonPolicy(baseline_api_compatibility="la_stack_v0_4_3_api"),
    )

    special_gaps = [gap for gap in collection.gaps if gap.bench in bench_compare.VS_LINALG_D8_RELEASE_SIGNAL_BENCHES]
    assert [gap.bench for gap in special_gaps] == [
        "la_stack_lu_pivoting",
        "la_stack_lu_ill_conditioned",
        "la_stack_ldlt_ill_conditioned",
    ]
    assert all(not gap.missing_current and gap.missing_baseline for gap in special_gaps)
    assert not any("balanced_range" in comparison.bench for comparison in collection.comparisons)


def test_v043_adapter_still_requires_current_balanced_range_sample(tmp_path: Path) -> None:
    group = tmp_path / "d8"
    _write_estimates(
        group / "la_stack_det_from_lu_balanced_range" / "new" / "estimates.json",
        "median",
        10.0,
    )

    collection = bench_compare._collect_comparisons(
        tmp_path,
        "v0.4.3",
        "median",
        suite="vs_linalg",
        policy=bench_compare.ComparisonPolicy(baseline_api_compatibility="la_stack_v0_4_3_api"),
    )

    balanced_gaps = [gap for gap in collection.gaps if "balanced_range" in gap.bench]
    assert [(gap.bench, gap.missing_current, gap.missing_baseline) for gap in balanced_gaps] == [("la_stack_det_from_ldlt_balanced_range", True, False)]


def test_v043_adapter_validates_current_balanced_range_sample(tmp_path: Path) -> None:
    estimate = tmp_path / "d8" / "la_stack_det_from_lu_balanced_range" / "new" / "estimates.json"
    estimate.parent.mkdir(parents=True)
    estimate.write_text("{not json", encoding="utf-8")

    with pytest.raises(ValueError, match="malformed Criterion estimates JSON"):
        bench_compare._collect_comparisons(
            tmp_path,
            "v0.4.3",
            "median",
            suite="vs_linalg",
            policy=bench_compare.ComparisonPolicy(baseline_api_compatibility="la_stack_v0_4_3_api"),
        )


def test_collect_vs_linalg_release_signal_uses_baseline_peer_context(tmp_path: Path) -> None:
    _build_vs_linalg_tree(tmp_path)
    collection = bench_compare._collect_comparisons(tmp_path, "last", "median", suite="vs_linalg")
    comparisons = collection.comparisons

    assert [c.bench for c in comparisons] == ["la_stack_lu_solve", "la_stack_ldlt_solve"]
    lu = comparisons[0]
    assert lu.baseline_ns == 20.0
    assert lu.current_ns == 10.0
    assert lu.baseline_nalgebra_ns == 30.0
    assert lu.baseline_faer_ns == 40.0

    ldlt = comparisons[1]
    assert ldlt.baseline_nalgebra_ns == 36.0
    assert ldlt.baseline_faer_ns == 48.0


def test_collect_vs_linalg_all_benches_includes_latest_peer_rows(tmp_path: Path) -> None:
    _build_vs_linalg_tree(tmp_path)
    collection = bench_compare._collect_comparisons(
        tmp_path,
        "last",
        "median",
        suite="vs_linalg",
        policy=bench_compare.ComparisonPolicy(scope="all-benches"),
    )
    comparisons = collection.comparisons

    assert [c.bench for c in comparisons] == [
        "la_stack_lu_solve",
        "nalgebra_lu_solve",
        "la_stack_ldlt_solve",
    ]
    assert {(gap.bench, gap.missing_current, gap.missing_baseline) for gap in collection.gaps} == {
        ("faer_lu_solve", True, False),
        ("nalgebra_cholesky_solve", True, False),
        ("faer_ldlt_solve", True, False),
    }


# ---------------------------------------------------------------------------
# Table generation
# ---------------------------------------------------------------------------


def test_snapshot_tables_per_dimension(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    results = bench_compare._collect_results(tmp_path, "new", "median")
    tables = bench_compare._snapshot_tables(results, "median")
    assert "### D=2" in tables
    assert "### D=3" in tables
    assert "### Random corpus D=3" in tables
    assert "### Near-singular 3x3" in tables
    assert "| Benchmark | Median | Criterion CI |" in tables


def test_snapshot_tables_uses_stat_label(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path, stat="mean")
    results = bench_compare._collect_results(tmp_path, "new", "mean")
    tables = bench_compare._snapshot_tables(results, "mean")
    assert "| Benchmark | Mean | Criterion CI |" in tables
    assert "Median" not in tables


def test_comparison_tables_per_dimension(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    for d, det, det_exact in [(2, 2.0, 8000.0), (3, 10.0, 42000.0)]:
        group = tmp_path / f"exact_d{d}"
        _write_estimates(group / "det" / "v0.3.0" / "estimates.json", "median", det)
        _write_estimates(group / "det_exact" / "v0.3.0" / "estimates.json", "median", det_exact)
        _write_estimates(group / "det_exact_f64" / "v0.3.0" / "estimates.json", "median", det_exact * 1.1)
        _write_estimates(group / "solve_exact_f64" / "v0.3.0" / "estimates.json", "median", det_exact * 1.3)

    comparisons = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median").comparisons
    tables = bench_compare._comparison_tables(comparisons, "v0.3.0")
    assert "### D=2" in tables
    assert "### D=3" in tables
    assert "| Benchmark | v0.3.0 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |" in tables
    assert "det_exact_rounded_f64 (vs det_exact_f64)" in tables
    assert "solve_exact_f64_result (vs solve_exact_f64)" in tables


def test_comparison_tables_include_vs_linalg_peer_context(tmp_path: Path) -> None:
    _build_vs_linalg_tree(tmp_path)
    comparisons = bench_compare._collect_comparisons(tmp_path, "last", "median", suite="vs_linalg").comparisons
    tables = bench_compare._comparison_tables(comparisons, "last")

    assert (
        "| Benchmark | last (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | "
        "Point-estimate ratio | last nalgebra | last faer |" in tables
    )
    assert (
        "| la_stack_lu_solve | 20.0 ns [18.0 ns, 22.0 ns] | 10.0 ns [9.0 ns, 11.0 ns] | -50.0% | "
        "faster point estimate; marginal CIs separated | 2.00x | "
        "30.0 ns [27.0 ns, 33.0 ns] | 40.0 ns [36.0 ns, 44.0 ns] |" in tables
    )
    assert (
        "| la_stack_ldlt_solve | 24.0 ns [21.6 ns, 26.4 ns] | 12.0 ns [10.8 ns, 13.2 ns] | -50.0% | "
        "faster point estimate; marginal CIs separated | 2.00x | "
        "36.0 ns [32.4 ns, 39.6 ns] | 48.0 ns [43.2 ns, 52.8 ns] |" in tables
    )


def test_coverage_table_makes_missing_samples_explicit(tmp_path: Path) -> None:
    group = tmp_path / "exact_d2"
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 10.0)
    _write_estimates(group / "det_direct" / "last" / "estimates.json", "median", 20.0)
    collection = bench_compare._collect_comparisons(tmp_path, "last", "median", suite="exact")

    table = bench_compare._coverage_table(collection.gaps, "last")
    assert "## Incomplete Comparison Coverage" in table
    assert "| Exact arithmetic | exact_d2 | det | baseline (`last`) |" in table
    assert "| Exact arithmetic | exact_d2 | det_direct | current (`new`) |" in table
    assert "| Exact arithmetic | exact_d2 | det_exact | current (`new`), baseline (`last`) |" in table


# ---------------------------------------------------------------------------
# Harness provenance
# ---------------------------------------------------------------------------


def test_read_harness_provenance_validates_shared_harness_metadata(tmp_path: Path) -> None:
    _write_harness_provenance(tmp_path)

    provenance = _read_harness_provenance(tmp_path)

    assert provenance == bench_compare.HarnessProvenance(
        schema=1,
        mode="shared-current-harness",
        sha256="a" * 64,
        baseline="v0.4.3",
    )


def test_read_harness_provenance_is_optional(tmp_path: Path) -> None:
    assert _read_harness_provenance(tmp_path) is None


def test_read_schema2_provenance_records_versions_dirty_source_and_both_gates(tmp_path: Path) -> None:
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / ".la-stack-benchmark-harness.json").write_text(
        json.dumps(_schema2_provenance_data()),
        encoding="utf-8",
    )

    provenance = _read_harness_provenance(tmp_path)

    assert provenance is not None
    assert provenance.schema == 2
    assert provenance.criterion is not None
    assert provenance.criterion.criterion_version == "0.7.0"
    markdown = bench_compare._provenance_markdown(provenance)
    rendered = "\n".join(markdown)
    assert "Criterion dependency version: `0.7.0`" in rendered
    assert "Current Git clean: `false`" in rendered
    assert "Validated baseline revision: `baseline-commit`" in rendered
    assert "Baseline API compatibility: `la_stack_v0_4_3_api`" in rendered
    assert "d8/la_stack_det_from_lu_balanced_range" in rendered
    assert "d8/la_stack_det_from_ldlt_balanced_range" in rendered
    assert "exact determinant is one" in rendered
    assert bench_compare._comparison_policy("release-signal", provenance) == bench_compare.ComparisonPolicy(
        scope="release-signal",
        baseline_api_compatibility="la_stack_v0_4_3_api",
    )


def test_historical_asset_provenance_uses_mode_appropriate_gate_wording(tmp_path: Path) -> None:
    data = _schema2_provenance_data()
    data["mode"] = "historical-assets"
    data["measurement"] = {
        "status": "unavailable",
        "reason": "historical release assets do not record the timing environment",
    }
    (tmp_path / ".la-stack-benchmark-harness.json").write_text(json.dumps(data), encoding="utf-8")

    provenance = _read_harness_provenance(tmp_path)

    assert provenance is not None
    rendered = "\n".join(bench_compare._provenance_markdown(provenance))
    assert "historical release assets do not record the timing environment" in rendered
    assert "passed for both referenced source revisions during publication" in rendered
    assert "separately from the historical timing measurements" in rendered
    assert "shared current fixture harness" not in rendered
    assert "both samples under one shared current harness" not in rendered


def test_read_schema2_provenance_requires_criterion_version(tmp_path: Path) -> None:
    data = _schema2_provenance_data()
    criterion = data["criterion"]
    assert isinstance(criterion, dict)
    del cast("dict[str, object]", criterion)["criterion_version"]
    (tmp_path / ".la-stack-benchmark-harness.json").write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(ValueError, match="criterion_version"):
        _read_harness_provenance(tmp_path)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("suite", "exact"),
        ("scope", "all-benches"),
        ("statistic", "mean"),
        ("sample", "v0.4.3"),
    ],
)
def test_read_schema2_provenance_binds_criterion_settings_to_request(
    tmp_path: Path,
    field: str,
    value: str,
) -> None:
    data = _schema2_provenance_data()
    criterion = data["criterion"]
    assert isinstance(criterion, dict)
    cast("dict[str, object]", criterion)[field] = value
    (tmp_path / ".la-stack-benchmark-harness.json").write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(ValueError, match=rf"criterion\.{field}.*does not match requested value"):
        _read_harness_provenance(tmp_path)


def test_read_schema2_provenance_rejects_v043_adapter_for_other_baseline(tmp_path: Path) -> None:
    data = _schema2_provenance_data()
    data["baseline"] = "v0.4.4"
    (tmp_path / ".la-stack-benchmark-harness.json").write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(ValueError, match=r"valid only for baseline 'v0\.4\.3'"):
        _read_harness_provenance(tmp_path, baseline="v0.4.4")


def test_read_harness_provenance_rejects_different_requested_baseline(tmp_path: Path) -> None:
    _write_harness_provenance(tmp_path, baseline="v0.4.3")

    with pytest.raises(ValueError, match="does not match requested Criterion baseline 'last'"):
        _read_harness_provenance(tmp_path, baseline="last")


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("schema", 3, "unsupported or missing schema"),
        ("mode", "independent-harnesses", "unsupported or missing mode"),
        ("sha256", "not-a-digest", "invalid or missing sha256"),
        ("baseline", "", "invalid or missing baseline"),
    ],
)
def test_read_harness_provenance_rejects_malformed_fields(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    data: dict[str, object] = {
        "schema": 1,
        "mode": "shared-current-harness",
        "sha256": "a" * 64,
        "baseline": "v0.4.3",
    }
    data[field] = value
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / ".la-stack-benchmark-harness.json").write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        _read_harness_provenance(tmp_path)


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------


def test_main_snapshot_writes_output(tmp_path: Path) -> None:
    criterion_dir = tmp_path / "criterion"
    _build_criterion_tree(criterion_dir)
    output = tmp_path / "PERFORMANCE.md"

    rc = bench_compare.main(
        [
            "--snapshot",
            "--suite",
            "exact",
            "--scope",
            "all-benches",
            "--criterion-dir",
            str(criterion_dir),
            "--output",
            str(output),
        ]
    )
    assert rc == 0
    assert output.exists()

    text = output.read_text(encoding="utf-8")
    assert "### D=2" in text
    assert "### Random corpus D=3" in text
    assert "### Near-singular 3x3" in text
    assert "just performance-local" in text
    assert "just performance-release" in text
    assert "just performance-github-assets" in text
    assert "just performance-release <current-tag> <previous-tag>" in text
    assert "Older curated release-to-release reports are archived in `docs/archive/performance/`." in text
    assert "git checkout" not in text


def test_main_no_criterion_dir(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rc = bench_compare.main(["--criterion-dir", str(tmp_path / "nonexistent"), "--output", str(tmp_path / "out.md")])
    assert rc == 2
    assert "No Criterion results" in capsys.readouterr().err


def test_main_comparison_no_baseline(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    criterion_dir = tmp_path / "criterion"
    _build_criterion_tree(criterion_dir)

    rc = bench_compare.main(["v0.3.0", "--criterion-dir", str(criterion_dir), "--output", str(tmp_path / "out.md")])
    assert rc == 2
    assert "No comparison data" in capsys.readouterr().err


def test_main_comparison_refuses_incomplete_coverage_before_writing(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    criterion_dir = tmp_path / "criterion"
    group = criterion_dir / "exact_d2"
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 10.0)
    _write_estimates(group / "det" / "v0.4.3" / "estimates.json", "median", 20.0)
    _write_harness_provenance(criterion_dir)
    output = tmp_path / "report.md"

    rc = bench_compare.main(["v0.4.3", "--suite", "exact", "--criterion-dir", str(criterion_dir), "--output", str(output)])

    assert rc == 2
    assert not output.exists()
    error = capsys.readouterr().err
    assert "Incomplete benchmark coverage" in error
    assert "## Incomplete Comparison Coverage" in error


def test_main_rejects_invalid_timing_without_writing_or_traceback(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    criterion_dir = tmp_path / "criterion"
    group = criterion_dir / "exact_d2"
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 0.0)
    _write_estimates(group / "det" / "last" / "estimates.json", "median", 10.0)
    output = tmp_path / "report.md"

    rc = bench_compare.main(["last", "--suite", "exact", "--criterion-dir", str(criterion_dir), "--output", str(output)])

    assert rc == 2
    assert "Invalid Criterion estimate data" in capsys.readouterr().err
    assert not output.exists()


def test_main_rejects_overflowing_timing_without_writing_or_traceback(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    criterion_dir = tmp_path / "criterion"
    current = criterion_dir / "exact_d2" / "det" / "new" / "estimates.json"
    current.parent.mkdir(parents=True)
    current.write_text(json.dumps({"median": {"point_estimate": _OVERFLOWING_TIMING}}), encoding="utf-8")
    _write_estimates(criterion_dir / "exact_d2" / "det" / "last" / "estimates.json", "median", 10.0)
    output = tmp_path / "report.md"

    rc = bench_compare.main(["last", "--suite", "exact", "--criterion-dir", str(criterion_dir), "--output", str(output)])

    assert rc == 2
    assert "Invalid Criterion estimate data" in capsys.readouterr().err
    assert not output.exists()


def test_main_v043_comparison_allows_only_unavailable_balanced_baselines(tmp_path: Path) -> None:
    criterion_dir = tmp_path / "criterion"
    unavailable = bench_compare._V0_4_3_UNAVAILABLE_BASELINE_ROWS
    for dimension in bench_compare.VS_LINALG_CANONICAL_DIMS:
        group = f"d{dimension}"
        benches = {
            *bench_compare.VS_LINALG_LA_STACK_BENCHES,
            *bench_compare.VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM.get(dimension, []),
        }
        for bench in benches:
            _write_estimates(criterion_dir / group / bench / "new" / "estimates.json", "median", 10.0)
            if (group, bench) not in unavailable:
                _write_estimates(criterion_dir / group / bench / "v0.4.3" / "estimates.json", "median", 20.0)

    provenance = _schema2_provenance_data()
    criterion = provenance["criterion"]
    assert isinstance(criterion, dict)
    cast("dict[str, object]", criterion)["suite"] = "vs_linalg"
    (criterion_dir / ".la-stack-benchmark-harness.json").write_text(json.dumps(provenance), encoding="utf-8")
    output = tmp_path / "report.md"

    rc = bench_compare.main(["v0.4.3", "--suite", "vs_linalg", "--criterion-dir", str(criterion_dir), "--output", str(output)])

    assert rc == 0
    rendered = output.read_text(encoding="utf-8")
    assert "d8/la_stack_det_from_lu_balanced_range" in rendered
    assert "d8/la_stack_det_from_ldlt_balanced_range" in rendered
    assert "no speedup is claimed" in rendered


def test_generate_markdown_labels_absent_provenance_unavailable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bench_compare, "_get_git_source_date", lambda _root: "2026-06-01 12:34:56 UTC")
    report = bench_compare._generate_markdown(
        tmp_path,
        "tables",
        bench_compare.ReportSettings(
            baseline_name="last",
            stat="median",
            suite="exact",
            scope="release-signal",
        ),
    )

    assert "**Reproducibility provenance**: unavailable" in report
    assert "CPU, OS, rustc, commit, dependency lock" in report
    assert "performance-improvement claim" in report
    assert "**Source revision timestamp**: 2026-06-01 12:34:56 UTC (deterministic report metadata; not the benchmark measurement time)" in report
    assert "**Benchmark measurement timestamp**: not recorded by Criterion" in report


def test_git_source_date_is_normalized_to_unambiguous_utc(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        bench_compare,
        "run_git_command",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 0, stdout="2026-06-01T05:34:56-07:00\n"),
    )

    assert bench_compare._get_git_source_date(tmp_path) == "2026-06-01 12:34:56 UTC"


def test_main_rejects_malformed_harness_provenance(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    criterion_dir = tmp_path / "criterion"
    criterion_dir.mkdir()
    (criterion_dir / ".la-stack-benchmark-harness.json").write_text("{not json", encoding="utf-8")

    rc = bench_compare.main(["last", "--criterion-dir", str(criterion_dir), "--output", str(tmp_path / "report.md")])

    assert rc == 2
    assert "Invalid benchmark harness provenance" in capsys.readouterr().err


def test_main_rejects_harness_provenance_for_a_different_baseline(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    criterion_dir = tmp_path / "criterion"
    group = criterion_dir / "exact_d2"
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 10.0)
    _write_estimates(group / "det" / "last" / "estimates.json", "median", 20.0)
    _write_harness_provenance(criterion_dir, baseline="v0.4.3")
    output = tmp_path / "report.md"

    rc = bench_compare.main(["last", "--suite", "exact", "--criterion-dir", str(criterion_dir), "--output", str(output)])

    assert rc == 2
    assert "does not match requested Criterion baseline 'last'" in capsys.readouterr().err
    assert not output.exists()
