from __future__ import annotations

import json
import re
from typing import TYPE_CHECKING

import pytest

import bench_compare

if TYPE_CHECKING:
    from pathlib import Path


def _write_estimates(path: Path, stat: str, median: float) -> None:
    """Write a minimal Criterion estimates.json."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                stat: {
                    "point_estimate": median,
                    "confidence_interval": {"lower_bound": median * 0.9, "upper_bound": median * 1.1},
                }
            }
        ),
        encoding="utf-8",
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

    random_group = criterion_dir / "exact_random_percentile_d3"
    _write_estimates(random_group / "det_exact_p95" / "new" / "estimates.json", stat, 33000.0)
    _write_estimates(random_group / "solve_exact_f64_result_p95" / "new" / "estimates.json", stat, 54000.0)
    _write_estimates(random_group / "solve_exact_rounded_f64_p95" / "new" / "estimates.json", stat, 55000.0)


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
    def test_improvement_is_bold(self) -> None:
        assert bench_compare._format_pct(-10.0) == "**-10.0%**"

    def test_regression_not_bold(self) -> None:
        assert bench_compare._format_pct(10.0) == "+10.0%"

    def test_small_change(self) -> None:
        assert bench_compare._format_pct(0.5) == "+0.5%"


# ---------------------------------------------------------------------------
# Group heading
# ---------------------------------------------------------------------------


class TestGroupHeading:
    def test_dimension_group(self) -> None:
        assert bench_compare._group_heading("exact_d3") == "D=3"

    def test_near_singular(self) -> None:
        assert bench_compare._group_heading("exact_near_singular_3x3") == "Near-singular 3x3"

    def test_random_percentile(self) -> None:
        assert bench_compare._group_heading("exact_random_percentile_d4") == "Random percentile D=4"

    def test_large_entries(self) -> None:
        assert bench_compare._group_heading("exact_large_entries_3x3") == "Large entries 3x3"

    def test_hilbert_4x4(self) -> None:
        assert bench_compare._group_heading("exact_hilbert_4x4") == "Hilbert 4x4"

    def test_hilbert_5x5(self) -> None:
        assert bench_compare._group_heading("exact_hilbert_5x5") == "Hilbert 5x5"

    def test_unknown_passthrough(self) -> None:
        assert bench_compare._group_heading("something_else") == "something_else"


# ---------------------------------------------------------------------------
# read_estimate
# ---------------------------------------------------------------------------


def test_read_estimate_success(tmp_path: Path) -> None:
    est = tmp_path / "estimates.json"
    _write_estimates(est, "median", 42.0)
    point, lo, hi = bench_compare._read_estimate(est, "median")
    assert point == 42.0
    assert lo == pytest.approx(42.0 * 0.9)
    assert hi == pytest.approx(42.0 * 1.1)


def test_read_estimate_no_ci(tmp_path: Path) -> None:
    """When confidence_interval is missing, all three values equal the point estimate."""
    est = tmp_path / "estimates.json"
    est.parent.mkdir(parents=True, exist_ok=True)
    est.write_text(
        json.dumps({"median": {"point_estimate": 99.0}}),
        encoding="utf-8",
    )
    point, lo, hi = bench_compare._read_estimate(est, "median")
    assert point == 99.0
    assert lo == 99.0
    assert hi == 99.0


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


# ---------------------------------------------------------------------------
# collect_results / collect_comparisons
# ---------------------------------------------------------------------------


def test_collect_results(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    results = bench_compare._collect_results(tmp_path, "new", "median")
    assert len(results) == 18  # 6 benches x 2 dims + 3 near-singular + 3 random percentile
    groups = {r.group for r in results}
    assert "exact_d2" in groups
    assert "exact_d3" in groups
    assert "exact_random_percentile_d3" in groups
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
    random_group = tmp_path / "exact_random_percentile_d3"
    _write_estimates(random_group / "det_exact_p95" / "v0.3.0" / "estimates.json", "median", 66000.0)
    _write_estimates(random_group / "solve_exact_f64_p95" / "v0.3.0" / "estimates.json", "median", 108000.0)

    comparisons = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median")
    assert len(comparisons) == 12  # 6 benches x 2 dims (near-singular has no baseline)
    assert {c.group for c in comparisons} == {"exact_d2", "exact_d3"}
    for c in comparisons:
        assert c.speedup == pytest.approx(c.baseline_ns / c.current_ns)
    assert {(c.bench, c.baseline_bench) for c in comparisons if c.baseline_bench is not None} == {
        ("det_exact_f64_result", "det_exact_f64"),
        ("det_exact_rounded_f64", "det_exact_f64"),
        ("solve_exact_f64_result", "solve_exact_f64"),
        ("solve_exact_rounded_f64", "solve_exact_f64"),
    }

    all_comparisons = bench_compare._collect_comparisons(
        tmp_path,
        "v0.3.0",
        "median",
        scope="all-benches",
    )
    assert len(all_comparisons) == 15
    random_comparisons = [c for c in all_comparisons if c.group == "exact_random_percentile_d3"]
    assert {(c.bench, c.baseline_bench) for c in random_comparisons} == {
        ("det_exact_p95", None),
        ("solve_exact_f64_result_p95", "solve_exact_f64_p95"),
        ("solve_exact_rounded_f64_p95", "solve_exact_f64_p95"),
    }


def test_collect_comparisons_zero_current(tmp_path: Path) -> None:
    """When the current estimate is zero, speedup should be infinity."""
    group = tmp_path / "exact_d2"
    # Current (new) has a zero point estimate.
    _write_estimates(group / "det" / "new" / "estimates.json", "median", 0.0)
    # Baseline has a normal value.
    _write_estimates(group / "det" / "v0.3.0" / "estimates.json", "median", 5.0)

    comparisons = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median")
    assert len(comparisons) == 1
    c = comparisons[0]
    assert c.speedup == float("inf")
    assert c.pct_change == pytest.approx(-100.0)


def test_collect_comparisons_missing_baseline(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    comparisons = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median")
    assert comparisons == []


def test_collect_vs_linalg_release_signal_uses_baseline_peer_context(tmp_path: Path) -> None:
    _build_vs_linalg_tree(tmp_path)
    comparisons = bench_compare._collect_comparisons(tmp_path, "last", "median")

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
    comparisons = bench_compare._collect_comparisons(tmp_path, "last", "median", scope="all-benches")

    assert [c.bench for c in comparisons] == [
        "la_stack_lu_solve",
        "nalgebra_lu_solve",
        "la_stack_ldlt_solve",
    ]


# ---------------------------------------------------------------------------
# Table generation
# ---------------------------------------------------------------------------


def test_snapshot_tables_per_dimension(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    results = bench_compare._collect_results(tmp_path, "new", "median")
    tables = bench_compare._snapshot_tables(results, "median")
    assert "### D=2" in tables
    assert "### D=3" in tables
    assert "### Random percentile D=3" in tables
    assert "### Near-singular 3x3" in tables
    assert "| Benchmark | Median | 95% CI |" in tables


def test_snapshot_tables_uses_stat_label(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path, stat="mean")
    results = bench_compare._collect_results(tmp_path, "new", "mean")
    tables = bench_compare._snapshot_tables(results, "mean")
    assert "| Benchmark | Mean | 95% CI |" in tables
    assert "Median" not in tables


def test_comparison_tables_per_dimension(tmp_path: Path) -> None:
    _build_criterion_tree(tmp_path)
    for d, det, det_exact in [(2, 2.0, 8000.0), (3, 10.0, 42000.0)]:
        group = tmp_path / f"exact_d{d}"
        _write_estimates(group / "det" / "v0.3.0" / "estimates.json", "median", det)
        _write_estimates(group / "det_exact" / "v0.3.0" / "estimates.json", "median", det_exact)
        _write_estimates(group / "det_exact_f64" / "v0.3.0" / "estimates.json", "median", det_exact * 1.1)
        _write_estimates(group / "solve_exact_f64" / "v0.3.0" / "estimates.json", "median", det_exact * 1.3)

    comparisons = bench_compare._collect_comparisons(tmp_path, "v0.3.0", "median")
    tables = bench_compare._comparison_tables(comparisons, "v0.3.0")
    assert "### D=2" in tables
    assert "### D=3" in tables
    assert "| Benchmark | v0.3.0 | Latest | Change | Speedup |" in tables
    assert "det_exact_rounded_f64 (vs det_exact_f64)" in tables
    assert "solve_exact_f64_result (vs solve_exact_f64)" in tables


def test_comparison_tables_include_vs_linalg_peer_context(tmp_path: Path) -> None:
    _build_vs_linalg_tree(tmp_path)
    comparisons = bench_compare._collect_comparisons(tmp_path, "last", "median")
    tables = bench_compare._comparison_tables(comparisons, "last")

    assert "| Benchmark | last | Latest | Change | Speedup | last nalgebra | last faer |" in tables
    assert "| la_stack_lu_solve | 20.0 ns | 10.0 ns | **-50.0%** | 2.00x | 30.0 ns | 40.0 ns |" in tables
    assert "| la_stack_ldlt_solve | 24.0 ns | 12.0 ns | **-50.0%** | 2.00x | 36.0 ns | 48.0 ns |" in tables


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------


def test_main_snapshot_writes_output(tmp_path: Path) -> None:
    criterion_dir = tmp_path / "criterion"
    _build_criterion_tree(criterion_dir)
    output = tmp_path / "PERFORMANCE.md"

    rc = bench_compare.main(["--snapshot", "--criterion-dir", str(criterion_dir), "--output", str(output)])
    assert rc == 0
    assert output.exists()

    text = output.read_text(encoding="utf-8")
    assert "### D=2" in text
    assert "### Random percentile D=3" in text
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
