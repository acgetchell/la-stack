"""Tests for Criterion dimension-report generation and README updates."""

import argparse
import hashlib
import json
import re
import subprocess
import tomllib
from dataclasses import replace
from types import SimpleNamespace
from typing import TYPE_CHECKING, cast

import pytest

import criterion_dim_plot
from performance_artifacts import (
    ArtifactContext,
    ArtifactPaths,
    PerformanceBundle,
    PerformanceRow,
    ReleasePair,
    ReportSource,
    TimingEstimate,
    write_bundle,
)

if TYPE_CHECKING:
    from pathlib import Path

_OVERFLOWING_TIMING = 10**400


def _toml_dependency_version(data: dict[str, object], name: str) -> str | None:
    for section in ("dependencies", "dev-dependencies", "build-dependencies"):
        table = data.get(section)
        if not isinstance(table, dict):
            continue
        table_dict = cast("dict[str, object]", table)
        value = table_dict.get(name)
        if value is None:
            continue
        if isinstance(value, str):
            return value
        if isinstance(value, dict):
            value_dict = cast("dict[str, object]", value)
            version = value_dict.get("version")
            if isinstance(version, str):
                return version
    return None


def test_detect_versions_matches_cargo_toml() -> None:
    root = criterion_dim_plot._repo_root()
    cargo_toml = root / "Cargo.toml"

    data = tomllib.loads(cargo_toml.read_text(encoding="utf-8"))

    package_version: str | None = None
    package = data.get("package")
    if isinstance(package, dict):
        version = package.get("version")
        if isinstance(version, str):
            package_version = version

    expected_la = package_version or "unknown"
    expected_na = _toml_dependency_version(data, "nalgebra") or "unknown"
    expected_fa = _toml_dependency_version(data, "faer") or "unknown"

    versions = criterion_dim_plot._detect_versions(root)

    assert versions["la-stack"] == expected_la
    assert versions["nalgebra"] == expected_na
    assert versions["faer"] == expected_fa


def test_readme_table_markers_are_stable() -> None:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    assert begin == "<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->"
    assert end == "<!-- BENCH_TABLE:lu_solve:median:new:END -->"


def test_markdown_table_formats_values_and_pct() -> None:
    rows = [
        criterion_dim_plot.Row(
            dim=2,
            la_time=50.0,
            la_lo=45.0,
            la_hi=55.0,
            na_time=100.0,
            na_lo=90.0,
            na_hi=110.0,
            fa_time=200.0,
            fa_lo=180.0,
            fa_hi=220.0,
        ),  # +50.0% vs na, +75.0% vs fa
        criterion_dim_plot.Row(
            dim=64,
            la_time=1_000.0,
            la_lo=900.0,
            la_hi=1_100.0,
            na_time=900.0,
            na_lo=800.0,
            na_hi=1_000.0,
            fa_time=800.0,
            fa_lo=700.0,
            fa_hi=900.0,
        ),  # -11.1% vs na, -25.0% vs fa
    ]

    table = criterion_dim_plot._markdown_table(rows, stat="median")

    assert (
        "| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | reduction vs nalgebra (point est.) | reduction vs faer (point est.) |" in table
    )
    assert "| 2 | 50.000 | 100.000 | 200.000 | +50.0% | +75.0% |" in table
    # thousand separator and sign
    assert "| 64 | 1,000.000 | 900.000 | 800.000 | -11.1% | -25.0% |" in table


def test_row_rejects_zero_peer_time_before_markdown_rendering() -> None:
    with pytest.raises(ValueError, match="na_time must be finite and positive"):
        criterion_dim_plot.Row(
            dim=2,
            la_time=10.0,
            la_lo=9.0,
            la_hi=11.0,
            na_time=0.0,
            na_lo=0.0,
            na_hi=0.0,
            fa_time=100.0,
            fa_lo=90.0,
            fa_hi=110.0,
        )


def test_gp_quote_escapes_backslashes_and_quotes() -> None:
    assert criterion_dim_plot._gp_quote("plain") == "'plain'"
    assert criterion_dim_plot._gp_quote("a'b") == "'a\\'b'"
    assert criterion_dim_plot._gp_quote("a\\b") == "'a\\\\b'"
    assert criterion_dim_plot._gp_quote("a\\'b") == "'a\\\\\\'b'"


def test_maybe_render_plot_handles_gnuplot_failure(capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    # Simulate gnuplot existing but failing to run (CalledProcessError).
    def boom(_req: object) -> None:
        raise criterion_dim_plot.subprocess.CalledProcessError(1, ["gnuplot"])

    monkeypatch.setattr(criterion_dim_plot, "_render_svg_with_gnuplot", boom)

    args = argparse.Namespace(no_plot=False)
    req = criterion_dim_plot.PlotRequest(
        csv_path=criterion_dim_plot.Path("out.csv"),
        out_svg=criterion_dim_plot.Path("out.svg"),
        title="t",
        stat="median",
        dims=(2,),
        la_label="la-stack v0.1.2",
        na_label="nalgebra v0.34.1",
        fa_label="faer v0.24.0",
        log_y=False,
    )

    rc = criterion_dim_plot._maybe_render_plot(args, req, skipped=[])
    assert rc == 1

    captured = capsys.readouterr()
    assert "Wrote CSV instead" in captured.err


def test_update_readme_table_replaces_only_between_markers(tmp_path: Path) -> None:
    marker_begin, marker_end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")

    readme = tmp_path / "README.md"
    readme.write_text(
        f"# Title\nbefore\n{marker_begin}\nold line 1\nold line 2\n{marker_end}\nafter\n",
        encoding="utf-8",
    )

    table_md = "| a |\n|---|\n| 1 |"

    changed = criterion_dim_plot._update_readme_table(readme, marker_begin, marker_end, table_md)
    assert changed is True

    text = readme.read_text(encoding="utf-8")
    assert "old line 1" not in text
    assert "old line 2" not in text
    assert marker_begin in text
    assert marker_end in text
    assert "| a |" in text
    assert f"{marker_begin}\n\n| a |\n|---|\n| 1 |\n\n{marker_end}" in text

    # Re-running with the same content should be a no-op.
    changed_again = criterion_dim_plot._update_readme_table(readme, marker_begin, marker_end, table_md)
    assert changed_again is False


def test_update_readme_table_errors_on_missing_markers(tmp_path: Path) -> None:
    readme = tmp_path / "README.md"
    readme.write_text("# Title\n", encoding="utf-8")

    with pytest.raises(criterion_dim_plot.ReadmeMarkerError, match=r"README markers not found"):
        criterion_dim_plot._update_readme_table(
            readme,
            "<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->",
            "<!-- BENCH_TABLE:lu_solve:median:new:END -->",
            "| x |",
        )


def test_update_readme_table_errors_on_out_of_order_markers(tmp_path: Path) -> None:
    marker_begin, marker_end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")

    readme = tmp_path / "README.md"
    readme.write_text(f"{marker_end}\n{marker_begin}\n", encoding="utf-8")

    with pytest.raises(criterion_dim_plot.ReadmeMarkerError, match=r"out of order"):
        criterion_dim_plot._update_readme_table(readme, marker_begin, marker_end, "| x |")


def test_update_readme_table_errors_on_non_unique_markers(tmp_path: Path) -> None:
    marker_begin, marker_end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")

    readme = tmp_path / "README.md"
    readme.write_text(
        f"{marker_begin}\n{marker_begin}\n{marker_end}\n",
        encoding="utf-8",
    )

    with pytest.raises(criterion_dim_plot.ReadmeMarkerError, match=r"not found or not unique"):
        criterion_dim_plot._update_readme_table(readme, marker_begin, marker_end, "| x |")


def _canonical_benchmark_readme(version: str) -> str:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    return (
        "[docs](https://github.com/acgetchell/la-stack/blob/v0.0.9/README.md)\n"
        f"[csv](https://github.com/acgetchell/la-stack/blob/v{version}/docs/assets/bench/vs_linalg_lu_solve_median.csv)\n"
        f"[provenance](https://github.com/acgetchell/la-stack/blob/v{version}/docs/assets/bench/vs_linalg_lu_solve_median.provenance.json)\n"
        f"[svg](https://raw.githubusercontent.com/acgetchell/la-stack/v{version}/docs/assets/bench/vs_linalg_lu_solve_median.svg)\n"
        f"{begin}\nold table\n{end}\n"
    )


def test_replace_readme_benchmark_asset_versions_updates_only_complete_selected_asset_set() -> None:
    updated = criterion_dim_plot._replace_readme_benchmark_asset_versions(
        _canonical_benchmark_readme("1.2.2"),
        metric="lu_solve",
        stat="median",
        version="1.2.3",
    )

    assert updated.count("v1.2.3/docs/assets/bench/") == 3
    assert "blob/v0.0.9/README.md" in updated
    assert "old table" in updated


def test_replace_readme_benchmark_asset_versions_rejects_incomplete_asset_set() -> None:
    readme = _canonical_benchmark_readme("1.2.2").replace(
        "[provenance](https://github.com/acgetchell/la-stack/blob/v1.2.2/docs/assets/bench/vs_linalg_lu_solve_median.provenance.json)\n",
        "",
    )

    with pytest.raises(criterion_dim_plot.ReadmeBenchmarkLinkError, match="exactly one tag-pinned link"):
        criterion_dim_plot._replace_readme_benchmark_asset_versions(
            readme,
            metric="lu_solve",
            stat="median",
            version="1.2.3",
        )


def _publication_args() -> criterion_dim_plot.PlotCliArgs:
    return criterion_dim_plot.PlotCliArgs(
        metric="lu_solve",
        stat="median",
        sample="new",
        criterion_dir="target/criterion",
        performance_csv=criterion_dim_plot._DEFAULT_PERFORMANCE_CSV,
        out=None,
        csv=None,
        log_y=True,
        no_plot=False,
        update_readme=True,
        readme="README.md",
        allow_partial=False,
    )


def test_readme_publication_rejects_noncanonical_data_and_asset_paths_before_rendering(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    (tmp_path / "README.md").write_text(f"{begin}\nold\n{end}\n", encoding="utf-8")

    assert (
        criterion_dim_plot._validate_readme_target(
            tmp_path,
            replace(_publication_args(), performance_csv="stale-results.csv"),
        )
        == 2
    )
    assert "canonical performance-release input" in capsys.readouterr().err

    assert (
        criterion_dim_plot._validate_readme_target(
            tmp_path,
            replace(_publication_args(), csv="custom.csv"),
        )
        == 2
    )
    assert "canonical CSV/SVG destinations" in capsys.readouterr().err


def test_readme_publication_rejects_no_plot_before_rendering(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    (tmp_path / "README.md").write_text(f"{begin}\nold\n{end}\n", encoding="utf-8")

    assert (
        criterion_dim_plot._validate_readme_target(
            tmp_path,
            replace(_publication_args(), no_plot=True),
        )
        == 2
    )
    assert "--no-plot is exploratory-only" in capsys.readouterr().err


def test_fixture_readme_may_use_custom_asset_destinations(tmp_path: Path) -> None:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    fixture = tmp_path / "README.fixture.md"
    fixture.write_text(f"{begin}\nold\n{end}\n", encoding="utf-8")

    assert (
        criterion_dim_plot._validate_readme_target(
            tmp_path,
            replace(
                _publication_args(),
                readme=str(fixture),
                csv="custom.csv",
                out="custom.svg",
            ),
        )
        == 0
    )


def test_readme_publication_rejects_incomplete_benchmark_links_before_rendering(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (tmp_path / "Cargo.toml").write_text('[package]\nname = "fixture"\nversion = "1.2.3"\n', encoding="utf-8")
    readme = _canonical_benchmark_readme("1.2.2").replace(
        "[svg](https://raw.githubusercontent.com/acgetchell/la-stack/v1.2.2/docs/assets/bench/vs_linalg_lu_solve_median.svg)\n",
        "",
    )
    (tmp_path / "README.md").write_text(readme, encoding="utf-8")

    assert criterion_dim_plot._validate_readme_target(tmp_path, _publication_args()) == 2
    assert "exactly one tag-pinned link" in capsys.readouterr().err


_TEST_COMMIT = "0123456789abcdef0123456789abcdef01234567"


def _timing(value: float) -> TimingEstimate:
    return TimingEstimate(median_ns=value, ci_lower_ns=value * 0.9, ci_upper_ns=value * 1.1)


def _write_benchmark_checkout(root: Path, *, version: str = "0.1.0") -> None:
    (root / "Cargo.toml").write_text(
        f'[package]\nname = "fixture"\nversion = "{version}"\n[dev-dependencies]\ncriterion = "0.7"\nnalgebra = "0.34"\nfaer = "0.22"\n',
        encoding="utf-8",
    )
    (root / "Cargo.lock").write_text("version = 4\n", encoding="utf-8")
    (root / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.88.0"\n', encoding="utf-8")
    (root / "justfile").write_text("test-bench-inputs:\n", encoding="utf-8")
    for relative in (
        ".config/nextest.toml",
        "tests/exact_bench_config.rs",
        "tests/vs_linalg_inputs.rs",
        "benches/vs_linalg.rs",
        "src/lib.rs",
    ):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("// fixture\n", encoding="utf-8")


def _write_performance_bundle(
    root: Path,
    *,
    version: str = "0.1.0",
    omit_dim: int | None = None,
    omit_peer_dim: int | None = None,
    include_contract: bool = True,
) -> ArtifactPaths:
    harness_sha256, missing = criterion_dim_plot._provenance_harness_digest(root)
    source_sha256, source_missing = criterion_dim_plot._source_state_digest(root)
    assert missing == []
    assert source_missing is False
    cargo_lock_sha256 = hashlib.sha256((root / "Cargo.lock").read_bytes()).hexdigest()
    baseline_commit = "fedcba9876543210fedcba9876543210fedcba98"
    baseline_source_sha256 = "b" * 64
    environment: dict[str, object] = {
        "cargo_lock_sha256": cargo_lock_sha256,
        "commit": _TEST_COMMIT,
        "correctness_gate": "passed",
        "cpu": "test-cpu",
        "git_clean": False,
        "harness_sha256": harness_sha256,
        "os": "TestOS 1 arm64",
        "rustc": "rustc 1.88.0",
        "source_state_sha256": source_sha256,
    }
    if include_contract:
        environment["benchmark_contract_sha256"] = criterion_dim_plot.benchmark_contract_digest(root)
    current_tag = f"v{version}"
    baseline_tag = "v0.0.9"
    bundle = PerformanceBundle(
        context=ArtifactContext(
            release=ReleasePair(current=current_tag, baseline=baseline_tag),
            statistic="median",
            suite="all",
            scope="release-signal",
            source=ReportSource(
                version=version,
                commit=_TEST_COMMIT,
                ref="HEAD",
                revision_timestamp="2026-08-20 12:00:00 UTC",
            ),
            benchmark_provenance={
                "baseline": baseline_tag,
                "criterion": {
                    "baseline_command": ["cargo", "bench", "--baseline", baseline_tag],
                    "criterion_version": "0.7.0",
                    "current_command": ["cargo", "bench", "--save-baseline", "current"],
                    "sample": "new",
                    "scope": "release-signal",
                    "statistic": "median",
                    "suite": "all",
                },
                "measurement": {
                    "baseline_api_compatibility": "none",
                    "baseline_commit": baseline_commit,
                    "baseline_git_clean": False,
                    "baseline_source_state_sha256": baseline_source_sha256,
                    "current_commit": _TEST_COMMIT,
                    "current_git_clean": False,
                    "current_source_state_sha256": source_sha256,
                    "status": "recorded",
                    **{key: environment[key] for key in ("cargo_lock_sha256", "cpu", "harness_sha256", "os", "rustc")},
                    **({"benchmark_contract_sha256": environment["benchmark_contract_sha256"]} if "benchmark_contract_sha256" in environment else {}),
                },
                "mode": "shared-current-harness",
                "publication": environment,
                "schema": 2,
                "validation": {
                    "baseline_api_compatibility": "none",
                    "baseline_commit": baseline_commit,
                    "baseline_git_clean": False,
                    "baseline_revision": "passed",
                    "baseline_source_state_sha256": baseline_source_sha256,
                    "command": ["just", "test-bench-inputs"],
                    "current_commit": _TEST_COMMIT,
                    "current_git_clean": False,
                    "current_revision": "passed",
                    "current_source_state_sha256": source_sha256,
                    "harness": "shared-current",
                },
            },
        ),
        rows=tuple(
            PerformanceRow(
                suite="vs_linalg",
                scope="release-signal",
                benchmark_id=f"d{dim}/la_stack_lu_solve",
                group=f"d{dim}",
                benchmark="la_stack_lu_solve",
                baseline_benchmark="la_stack_lu_solve",
                coverage_status="comparable",
                coverage_note="",
                baseline=_timing(float(dim * 6)),
                current=_timing(float(dim * 5)),
                baseline_nalgebra=None if dim == omit_peer_dim else _timing(float(dim * 10)),
                baseline_faer=_timing(float(dim * 20)),
            )
            for dim in criterion_dim_plot.CANONICAL_DIMS
            if dim != omit_dim
        ),
    )
    paths = ArtifactPaths(
        csv=root / criterion_dim_plot._DEFAULT_PERFORMANCE_CSV,
        provenance=root / "target/bench-reports/performance.provenance.json",
    )
    write_bundle(paths, bundle)
    return paths


def test_main_update_readme_happy_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _write_benchmark_checkout(tmp_path)
    performance_paths = _write_performance_bundle(tmp_path)
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")

    def fake_run_git(args: list[str], **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(stdout=f"{_TEST_COMMIT}\n")

    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)
    monkeypatch.setattr(criterion_dim_plot, "run_git_command", fake_run_git)

    def fake_render(request: criterion_dim_plot.PlotRequest) -> None:
        request.out_svg.write_text("<svg/>\n", encoding="utf-8")

    monkeypatch.setattr(criterion_dim_plot, "_render_svg_with_gnuplot", fake_render)

    out_csv = tmp_path / "docs" / "assets" / "bench" / "vs_linalg_lu_solve_median.csv"
    out_svg = tmp_path / "docs" / "assets" / "bench" / "vs_linalg_lu_solve_median.svg"

    rc = criterion_dim_plot.main(
        [
            "--metric",
            "lu_solve",
            "--stat",
            "median",
            "--sample",
            "new",
            "--update-readme",
        ]
    )
    assert rc == 0
    assert out_svg.read_text(encoding="utf-8") == "<svg/>\n"

    # CSV written
    csv_text = out_csv.read_text(encoding="utf-8")
    assert csv_text.startswith("D,la_stack,la_lo,la_hi,nalgebra,na_lo,na_hi,faer,fa_lo,fa_hi\n")
    assert "2,10.0" in csv_text
    assert "64,320.0" in csv_text

    # README updated with computed table
    readme_text = readme.read_text(encoding="utf-8")
    assert "old table" not in readme_text
    assert "| 2 | 10.000 | 20.000 | 40.000 | +50.0% | +75.0% |" in readme_text
    assert "| 64 | 320.000 | 640.000 | 1,280.000 | +50.0% | +75.0% |" in readme_text
    assert readme_text.count("v0.1.0/docs/assets/bench/") == 3
    assert "blob/v0.0.9/README.md" in readme_text

    provenance = json.loads(out_csv.with_suffix(".provenance.json").read_text(encoding="utf-8"))
    assert provenance["measurement"]["status"] == "recorded"
    assert provenance["measurement"]["source"] == "retained performance-release artifact"
    assert provenance["measurement"]["la_stack_sample"] == "current"
    assert provenance["measurement"]["peer_release_context"] == "v0.0.9"
    assert provenance["measurement"]["peer_sample"] == "baseline phase under shared current harness"
    assert provenance["publication"]["correctness_gate"] == "validated-performance-release-artifact"
    assert provenance["publication"]["benchmark_contract"] == "matched"
    assert "git_clean" not in provenance["publication"]
    assert provenance["criterion"]["current_command"][:2] == ["cargo", "bench"]
    assert provenance["performance_artifact"]["csv"] == "target/bench-reports/performance.csv"
    assert provenance["performance_artifact"]["csv_sha256"] == hashlib.sha256(performance_paths.csv.read_bytes()).hexdigest()
    assert re.fullmatch(r"[0-9a-f]{64}", provenance["publication"]["source_state_sha256"])

    published = (readme, out_csv, out_svg, out_csv.with_suffix(".provenance.json"))
    first = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in published}
    assert criterion_dim_plot.main(["--metric", "lu_solve", "--stat", "median", "--sample", "new", "--update-readme"]) == 0
    assert {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in published} == first


def test_dim_parsing_and_discovery(tmp_path: Path) -> None:
    assert criterion_dim_plot._dim_from_group_dir("d2") == 2
    assert criterion_dim_plot._dim_from_group_dir("d10") == 10
    assert criterion_dim_plot._dim_from_group_dir("dx") is None
    assert criterion_dim_plot._dim_from_group_dir("2") is None

    (tmp_path / "d2").mkdir()
    (tmp_path / "d10").mkdir()
    (tmp_path / "not_a_dim").mkdir()
    dims = criterion_dim_plot._discover_dims(tmp_path)
    assert dims == [2, 10]


def test_toml_helpers_read_versions(tmp_path: Path) -> None:
    cargo_toml = tmp_path / "Cargo.toml"
    cargo_toml.write_text(
        (
            "# comment line\n"
            "[package]\n"
            'version = "1.2.3" # inline comment\n'
            "\n"
            "[dependencies]\n"
            'nalgebra = "0.34.0"\n'
            'faer = { version = "0.21.4" }\n'
            "\n"
            "[dev-dependencies]\n"
            'serde = "1.0"'
        ),
        encoding="utf-8",
    )

    assert criterion_dim_plot._read_cargo_package_version(cargo_toml) == "1.2.3"
    deps = criterion_dim_plot._read_cargo_dependency_versions(cargo_toml, {"nalgebra", "faer"})
    assert deps["nalgebra"] == "0.34.0"
    assert deps["faer"] == "0.21.4"


def test_format_legend_label() -> None:
    assert criterion_dim_plot._format_legend_label("la-stack", "0.1.0") == "la-stack v0.1.0"
    assert criterion_dim_plot._format_legend_label("faer", "unknown") == "faer"


def test_read_estimate_errors_and_success(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(
        json.dumps(
            {
                "median": {
                    "point_estimate": 5.0,
                    "confidence_interval": {"lower_bound": 4.0, "upper_bound": 6.0},
                }
            }
        ),
        encoding="utf-8",
    )
    point, lo, hi = criterion_dim_plot._read_estimate(estimates, "median")
    assert (point, lo, hi) == (5.0, 4.0, 6.0)

    with pytest.raises(KeyError, match="stat 'mean' not found"):
        criterion_dim_plot._read_estimate(estimates, "mean")


def test_read_estimate_malformed_json_names_file(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text("{not json", encoding="utf-8")

    with pytest.raises(ValueError, match=re.escape(f"malformed Criterion estimates JSON in {estimates}")):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_missing_point_estimate_names_field(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps({"median": {}}), encoding="utf-8")

    with pytest.raises(KeyError, match="field 'point_estimate' for stat 'median' not found"):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_non_numeric_ci_bound_names_field(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(
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
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_rejects_numeric_overflow(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps({"median": {"point_estimate": _OVERFLOWING_TIMING}}), encoding="utf-8")

    with pytest.raises(ValueError, match=r"field 'point_estimate'.*not numeric") as exc_info:
        criterion_dim_plot._read_estimate(estimates, "median")

    assert isinstance(exc_info.value.__cause__, OverflowError)


def test_read_estimate_rejects_missing_or_partial_confidence_interval(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps({"median": {"point_estimate": 1.0}}), encoding="utf-8")
    with pytest.raises(KeyError, match="field 'confidence_interval'"):
        criterion_dim_plot._read_estimate(estimates, "median")

    estimates.write_text(
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
        criterion_dim_plot._read_estimate(estimates, "median")


@pytest.mark.parametrize(
    ("payload", "field"),
    [
        ({"median": {"point_estimate": True}}, "point_estimate"),
        (
            {
                "median": {
                    "point_estimate": 1.0,
                    "confidence_interval": {"lower_bound": False, "upper_bound": 2.0},
                }
            },
            "lower_bound",
        ),
    ],
)
def test_read_estimate_rejects_boolean_numeric_fields(tmp_path: Path, payload: dict[str, object], field: str) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(TypeError, match=rf"field '{field}' for stat 'median'.*not numeric"):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_rejects_nonfinite_time(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps({"median": {"point_estimate": "NaN"}}), encoding="utf-8")

    with pytest.raises(ValueError, match=r"median\.point_estimate.*finite and positive"):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_rejects_negative_time(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps({"median": {"point_estimate": -1.0}}), encoding="utf-8")

    with pytest.raises(ValueError, match=r"median\.point_estimate.*finite and positive"):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_rejects_zero_time(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(json.dumps({"median": {"point_estimate": 0.0}}), encoding="utf-8")

    with pytest.raises(ValueError, match=r"median\.point_estimate.*finite and positive"):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_rejects_inverted_confidence_interval(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(
        json.dumps(
            {
                "median": {
                    "point_estimate": 5.0,
                    "confidence_interval": {"lower_bound": 6.0, "upper_bound": 4.0},
                }
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="lower bound must be <= upper bound"):
        criterion_dim_plot._read_estimate(estimates, "median")


def test_read_estimate_allows_point_outside_percentile_confidence_interval(tmp_path: Path) -> None:
    estimates = tmp_path / "estimates.json"
    estimates.write_text(
        json.dumps(
            {
                "median": {
                    "point_estimate": 5.0,
                    "confidence_interval": {"lower_bound": 1.0, "upper_bound": 4.0},
                }
            }
        ),
        encoding="utf-8",
    )

    assert criterion_dim_plot._read_estimate(estimates, "median") == (5.0, 1.0, 4.0)


def test_row_rejects_invalid_dimension_and_times() -> None:
    with pytest.raises(ValueError, match="dimension must be positive"):
        criterion_dim_plot.Row(
            dim=0,
            la_time=1.0,
            la_lo=1.0,
            la_hi=1.0,
            na_time=1.0,
            na_lo=1.0,
            na_hi=1.0,
            fa_time=1.0,
            fa_lo=1.0,
            fa_hi=1.0,
        )

    with pytest.raises(ValueError, match="la_time must be finite and positive"):
        criterion_dim_plot.Row(
            dim=2,
            la_time=float("inf"),
            la_lo=1.0,
            la_hi=1.0,
            na_time=1.0,
            na_lo=1.0,
            na_hi=1.0,
            fa_time=1.0,
            fa_lo=1.0,
            fa_hi=1.0,
        )

    row = criterion_dim_plot.Row(
        dim=2,
        la_time=5.0,
        la_lo=1.0,
        la_hi=4.0,
        na_time=1.0,
        na_lo=1.0,
        na_hi=1.0,
        fa_time=1.0,
        fa_lo=1.0,
        fa_hi=1.0,
    )
    assert row.la_time == 5.0


def test_write_csv_and_collect_rows(tmp_path: Path) -> None:
    rows = [
        criterion_dim_plot.Row(
            dim=2,
            la_time=1.0,
            la_lo=0.9,
            la_hi=1.1,
            na_time=2.0,
            na_lo=1.9,
            na_hi=2.1,
            fa_time=3.0,
            fa_lo=2.9,
            fa_hi=3.1,
        )
    ]
    out_csv = tmp_path / "out.csv"
    criterion_dim_plot._write_csv(out_csv, rows)
    text = out_csv.read_text(encoding="utf-8")
    assert text.startswith("D,la_stack,la_lo,la_hi,nalgebra,na_lo,na_hi,faer,fa_lo,fa_hi")
    assert "2,1.0,0.9,1.1,2.0,1.9,2.1,3.0,2.9,3.1" in text

    criterion_dir = tmp_path / "criterion"
    metric = criterion_dim_plot.METRICS["lu_solve"]
    d2 = criterion_dir / "d2"
    d2.mkdir(parents=True)
    # Only la_stack exists; should be skipped.
    (d2 / metric.la_bench / "new").mkdir(parents=True)
    (d2 / metric.la_bench / "new" / "estimates.json").write_text(
        json.dumps({"median": {"point_estimate": 1.0}}),
        encoding="utf-8",
    )
    rows2, skipped = criterion_dim_plot._collect_rows(criterion_dir, [2], metric, "median", "new")
    assert rows2 == []
    assert skipped == ["d2 (missing nalgebra_lu_solve, faer_lu_solve)"]


def test_resolve_paths(tmp_path: Path) -> None:
    root = tmp_path
    resolved = criterion_dim_plot._resolve_under_root(root, "foo/bar.csv")
    assert resolved == root / "foo/bar.csv"

    svg, csv = criterion_dim_plot._resolve_output_paths(root, "lu_solve", "median", None, None)
    assert svg == root / "docs/assets/bench/vs_linalg_lu_solve_median.svg"
    assert csv == root / "docs/assets/bench/vs_linalg_lu_solve_median.csv"


def test_maybe_render_plot_no_plot_path(capsys: pytest.CaptureFixture[str]) -> None:
    args = argparse.Namespace(no_plot=True)
    req = criterion_dim_plot.PlotRequest(
        csv_path=criterion_dim_plot.Path("out.csv"),
        out_svg=criterion_dim_plot.Path("out.svg"),
        title="t",
        stat="median",
        dims=(2,),
        la_label="la",
        na_label="na",
        fa_label="fa",
        log_y=False,
    )
    rc = criterion_dim_plot._maybe_render_plot(args, req, skipped=[])
    assert rc == 0
    captured = capsys.readouterr()
    assert "Wrote CSV: out.csv" in captured.out


def test_maybe_render_plot_success_path(capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    def no_op(_req: object) -> None:
        return None

    monkeypatch.setattr(criterion_dim_plot, "_render_svg_with_gnuplot", no_op)

    args = argparse.Namespace(no_plot=False)
    req = criterion_dim_plot.PlotRequest(
        csv_path=criterion_dim_plot.Path("out.csv"),
        out_svg=criterion_dim_plot.Path("out.svg"),
        title="t",
        stat="median",
        dims=(2,),
        la_label="la",
        na_label="na",
        fa_label="fa",
        log_y=False,
    )

    rc = criterion_dim_plot._maybe_render_plot(args, req, skipped=["d2 (missing)"])
    assert rc == 0
    captured = capsys.readouterr()
    assert "Warning: some dimension groups were skipped:" in captured.out
    assert "Wrote CSV: out.csv" in captured.out
    assert "Wrote SVG: out.svg" in captured.out


def test_maybe_update_readme_errors(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    args = argparse.Namespace(
        update_readme=True,
        readme="missing.md",
        metric="lu_solve",
        stat="median",
        sample="new",
    )
    rc = criterion_dim_plot._maybe_update_readme(tmp_path, args, [])
    assert rc == 2
    captured = capsys.readouterr()
    assert "No such file or directory" in captured.err


def test_main_error_paths(tmp_path: Path) -> None:
    # Missing Criterion directory.
    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(tmp_path / "missing"),
            "--no-plot",
        ]
    )
    assert rc == 2

    # Criterion directory exists but has no usable rows.
    criterion_dir = tmp_path / "criterion"
    (criterion_dir / "d2").mkdir(parents=True)
    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(criterion_dir),
            "--metric",
            "lu_solve",
            "--stat",
            "median",
            "--sample",
            "new",
            "--no-plot",
        ]
    )
    assert rc == 2


def test_main_requires_canonical_dimensions_before_writing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    criterion_dir = tmp_path / "criterion"
    metric = criterion_dim_plot.METRICS["lu_solve"]
    for bench in (metric.la_bench, metric.na_bench, metric.fa_bench):
        estimates = criterion_dir / "d2" / bench / "new" / "estimates.json"
        estimates.parent.mkdir(parents=True, exist_ok=True)
        estimates.write_text(
            json.dumps(
                {
                    "median": {
                        "point_estimate": 1.0,
                        "confidence_interval": {"lower_bound": 0.9, "upper_bound": 1.1},
                    }
                }
            ),
            encoding="utf-8",
        )
    output = tmp_path / "out.csv"
    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)

    rc = criterion_dim_plot.main(["--criterion-dir", str(criterion_dir), "--csv", str(output), "--no-plot"])

    assert rc == 2
    assert not output.exists()
    assert not output.with_suffix(".provenance.json").exists()


def test_main_partial_mode_is_explicit_and_labels_measurement_unavailable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    criterion_dir = tmp_path / "criterion"
    metric = criterion_dim_plot.METRICS["lu_solve"]
    for bench in (metric.la_bench, metric.na_bench, metric.fa_bench):
        estimates = criterion_dir / "d2" / bench / "new" / "estimates.json"
        estimates.parent.mkdir(parents=True, exist_ok=True)
        estimates.write_text(
            json.dumps(
                {
                    "median": {
                        "point_estimate": 1.0,
                        "confidence_interval": {"lower_bound": 0.9, "upper_bound": 1.1},
                    }
                }
            ),
            encoding="utf-8",
        )
    output = tmp_path / "out.csv"
    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)

    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(criterion_dir),
            "--csv",
            str(output),
            "--no-plot",
            "--allow-partial",
        ]
    )

    assert rc == 0
    provenance = json.loads(output.with_suffix(".provenance.json").read_text(encoding="utf-8"))
    assert provenance["measurement"]["status"] == "unavailable"
    assert provenance["publication"]["correctness_gate"] == "not-run-exploratory"


def test_main_rejects_missing_confidence_interval_without_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    criterion_dir = tmp_path / "criterion"
    metric = criterion_dim_plot.METRICS["lu_solve"]
    for bench in (metric.la_bench, metric.na_bench, metric.fa_bench):
        estimates = criterion_dir / "d2" / bench / "new" / "estimates.json"
        estimates.parent.mkdir(parents=True, exist_ok=True)
        estimates.write_text(json.dumps({"median": {"point_estimate": 1.0}}), encoding="utf-8")
    output = tmp_path / "out.csv"
    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)

    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(criterion_dir),
            "--csv",
            str(output),
            "--no-plot",
            "--allow-partial",
        ]
    )

    assert rc == 2
    assert "Invalid Criterion estimate data" in capsys.readouterr().err
    assert not output.exists()


def test_main_rejects_overflowing_timing_without_writing_or_traceback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    criterion_dir = tmp_path / "criterion"
    metric = criterion_dim_plot.METRICS["lu_solve"]
    for bench in (metric.la_bench, metric.na_bench, metric.fa_bench):
        estimates = criterion_dir / "d2" / bench / "new" / "estimates.json"
        estimates.parent.mkdir(parents=True, exist_ok=True)
        estimates.write_text(json.dumps({"median": {"point_estimate": _OVERFLOWING_TIMING}}), encoding="utf-8")
    output = tmp_path / "out.csv"
    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)

    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(criterion_dir),
            "--csv",
            str(output),
            "--no-plot",
            "--allow-partial",
        ]
    )

    assert rc == 2
    assert "Invalid Criterion estimate data" in capsys.readouterr().err
    assert not output.exists()


@pytest.mark.parametrize(
    "failure",
    [
        subprocess.TimeoutExpired(["tool"], 17),
        OSError("working directory unavailable"),
    ],
)
def test_provenance_helpers_treat_timeout_and_os_error_as_unavailable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: Exception,
) -> None:
    def fail_command(*_args: object, **_kwargs: object) -> SimpleNamespace:
        raise failure

    monkeypatch.setattr(criterion_dim_plot, "run_git_command", fail_command)
    assert criterion_dim_plot._git_value(tmp_path, ["rev-parse", "HEAD"]) == "unavailable"
    git_clean, _status_digest = criterion_dim_plot._git_status_metadata(tmp_path)
    assert git_clean is None

    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", fail_command)
    assert criterion_dim_plot._rustc_version(tmp_path) == "unavailable"


def _mock_publication_environment(root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run_git(args: list[str], **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(stdout=f"{_TEST_COMMIT}\n")

    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: root)
    monkeypatch.setattr(criterion_dim_plot, "run_git_command", fake_run_git)


def test_main_publication_fails_closed_when_git_commit_is_unavailable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_benchmark_checkout(tmp_path)
    _write_performance_bundle(tmp_path)
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")
    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)

    def unavailable_git(*_args: object, **_kwargs: object) -> SimpleNamespace:
        message = "git unavailable"
        raise OSError(message)

    monkeypatch.setattr(criterion_dim_plot, "run_git_command", unavailable_git)

    assert criterion_dim_plot.main(["--update-readme"]) == 2
    assert "current checkout commit 'unavailable'" in capsys.readouterr().err
    assert "old table" in readme.read_text(encoding="utf-8")
    assert not (tmp_path / "docs/assets/bench/vs_linalg_lu_solve_median.csv").exists()


def test_main_publication_rejects_missing_performance_bundle(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_benchmark_checkout(tmp_path)
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")
    _mock_publication_environment(tmp_path, monkeypatch)

    assert criterion_dim_plot.main(["--update-readme"]) == 2
    assert "artifact pair is incomplete" in capsys.readouterr().err
    assert "old table" in readme.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "case",
    [
        (3, None, "missing d3/la_stack_lu_solve"),
        (None, 5, "missing nalgebra peer timing"),
    ],
)
def test_main_publication_rejects_incomplete_retained_coverage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    case: tuple[int | None, int | None, str],
) -> None:
    omit_dim, omit_peer_dim, expected = case
    _write_benchmark_checkout(tmp_path)
    _write_performance_bundle(tmp_path, omit_dim=omit_dim, omit_peer_dim=omit_peer_dim)
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")
    _mock_publication_environment(tmp_path, monkeypatch)

    assert criterion_dim_plot.main(["--update-readme"]) == 2
    assert expected in capsys.readouterr().err
    assert "old table" in readme.read_text(encoding="utf-8")


def test_main_publication_rejects_stale_harness_without_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_benchmark_checkout(tmp_path)
    _write_performance_bundle(tmp_path)
    (tmp_path / "benches/vs_linalg.rs").write_text("// changed after measurement\n", encoding="utf-8")
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")
    _mock_publication_environment(tmp_path, monkeypatch)

    assert criterion_dim_plot.main(["--update-readme"]) == 2
    assert "benchmark_contract_sha256" in capsys.readouterr().err
    assert "old table" in readme.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("include_contract", "expected_status"),
    [(True, "matched"), (False, "legacy-retained-artifact")],
)
def test_main_publication_ignores_justfile_only_changes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    include_contract: bool,
    expected_status: str,
) -> None:
    _write_benchmark_checkout(tmp_path)
    _write_performance_bundle(tmp_path, include_contract=include_contract)
    (tmp_path / "justfile").write_text("test-bench-inputs:\n\n# publication help changed\n", encoding="utf-8")
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")
    _mock_publication_environment(tmp_path, monkeypatch)
    monkeypatch.setattr(
        criterion_dim_plot,
        "_render_svg_with_gnuplot",
        lambda request: request.out_svg.write_text("<svg/>\n", encoding="utf-8"),
    )

    assert criterion_dim_plot.main(["--update-readme"]) == 0
    provenance = json.loads((tmp_path / "docs/assets/bench/vs_linalg_lu_solve_median.provenance.json").read_text(encoding="utf-8"))
    assert provenance["publication"]["benchmark_contract"] == expected_status


def test_main_publication_rejects_tampered_performance_csv(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_benchmark_checkout(tmp_path)
    paths = _write_performance_bundle(tmp_path)
    paths.csv.write_text(paths.csv.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
    readme = tmp_path / "README.md"
    readme.write_text(_canonical_benchmark_readme("0.0.8"), encoding="utf-8")
    _mock_publication_environment(tmp_path, monkeypatch)

    assert criterion_dim_plot.main(["--update-readme"]) == 2
    assert "CSV digest mismatch" in capsys.readouterr().err
    assert "old table" in readme.read_text(encoding="utf-8")


def test_staged_publication_leaves_existing_assets_unchanged_when_render_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    csv_path = tmp_path / "benchmark.csv"
    svg_path = tmp_path / "benchmark.svg"
    provenance_path = csv_path.with_suffix(".provenance.json")
    readme = tmp_path / "README.md"
    for path, text in (
        (csv_path, "old csv\n"),
        (svg_path, "old svg\n"),
        (provenance_path, "old provenance\n"),
    ):
        path.write_text(text, encoding="utf-8")
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    readme.write_text(f"before\n{begin}\nold table\n{end}\nafter\n", encoding="utf-8")

    def fail_render(_request: criterion_dim_plot.PlotRequest) -> None:
        raise subprocess.CalledProcessError(1, ["gnuplot"], stderr="render failed")

    monkeypatch.setattr(criterion_dim_plot, "_render_svg_with_gnuplot", fail_render)
    args = criterion_dim_plot.PlotCliArgs(
        metric="lu_solve",
        stat="median",
        sample="new",
        criterion_dir="target/criterion",
        performance_csv=criterion_dim_plot._DEFAULT_PERFORMANCE_CSV,
        out=str(svg_path),
        csv=str(csv_path),
        log_y=False,
        no_plot=False,
        update_readme=True,
        readme=str(readme),
        allow_partial=False,
    )
    row = criterion_dim_plot.Row(2, 1.0, 0.9, 1.1, 2.0, 1.9, 2.1, 3.0, 2.9, 3.1)
    request = criterion_dim_plot.PlotRequest(
        csv_path=csv_path,
        out_svg=svg_path,
        title="title",
        stat="median",
        dims=(2,),
        la_label="la-stack",
        na_label="nalgebra",
        fa_label="faer",
        log_y=False,
    )

    rc = criterion_dim_plot._stage_and_publish_outputs(
        root=tmp_path,
        args=args,
        rows=[row],
        req=request,
        provenance={"schema": 1},
        skipped=[],
    )

    assert rc == 1
    assert csv_path.read_text(encoding="utf-8") == "old csv\n"
    assert svg_path.read_text(encoding="utf-8") == "old svg\n"
    assert provenance_path.read_text(encoding="utf-8") == "old provenance\n"
    assert "old table" in readme.read_text(encoding="utf-8")


def test_artifact_rollback_failure_preserves_backups(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination_one = tmp_path / "one.txt"
    destination_two = tmp_path / "two.txt"
    staged_one = tmp_path / "staged-one.txt"
    staged_two = tmp_path / "staged-two.txt"
    backup_dir = tmp_path / "backups"
    backup_dir.mkdir()
    destination_one.write_text("old one\n", encoding="utf-8")
    destination_two.write_text("old two\n", encoding="utf-8")
    staged_one.write_text("new one\n", encoding="utf-8")
    staged_two.write_text("new two\n", encoding="utf-8")
    original_replace = criterion_dim_plot.Path.replace

    def fail_replacement_and_rollback(source: Path, destination: Path) -> Path:
        if source == staged_two and destination == destination_two:
            msg = "simulated publish failure"
            raise OSError(msg)
        if source == backup_dir / "backup-0" and destination == destination_one:
            msg = "simulated rollback failure"
            raise OSError(msg)
        return original_replace(source, destination)

    monkeypatch.setattr(criterion_dim_plot.Path, "replace", fail_replacement_and_rollback)

    with pytest.raises(criterion_dim_plot.PublicationRollbackError, match="backups preserved") as exc_info:
        criterion_dim_plot._replace_staged_files(
            [(staged_one, destination_one), (staged_two, destination_two)],
            backup_dir,
        )

    assert str(backup_dir) in str(exc_info.value)
    assert (backup_dir / "backup-0").read_text(encoding="utf-8") == "old one\n"
    assert destination_one.read_text(encoding="utf-8") == "new one\n"
    assert destination_two.read_text(encoding="utf-8") == "old two\n"


def test_repo_root_resolution_uses_working_checkout_for_installed_entrypoint(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    nested = tmp_path / "nested"
    nested.mkdir()
    (tmp_path / "Cargo.toml").write_text("[package]\n", encoding="utf-8")
    monkeypatch.chdir(nested)

    assert criterion_dim_plot._repo_root() == tmp_path


def test_publication_paths_reject_output_aliases(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    output = tmp_path / "benchmark.csv"
    args = criterion_dim_plot.PlotCliArgs(
        metric="lu_solve",
        stat="median",
        sample="new",
        criterion_dir="target/criterion",
        performance_csv=criterion_dim_plot._DEFAULT_PERFORMANCE_CSV,
        out=str(output),
        csv=str(output),
        log_y=False,
        no_plot=False,
        update_readme=False,
        readme="README.md",
        allow_partial=False,
    )

    rc = criterion_dim_plot._validate_publication_paths(tmp_path, args, out_svg=output, out_csv=output)

    assert rc == 2
    assert "must use distinct paths" in capsys.readouterr().err
