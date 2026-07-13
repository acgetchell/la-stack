"""Tests for Criterion dimension-report generation and README updates."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tomllib
from dataclasses import replace
from types import SimpleNamespace
from typing import TYPE_CHECKING, cast

import pytest

import criterion_dim_plot

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


@pytest.mark.parametrize(
    ("metric_name", "expected_filter"),
    [
        ("lu_solve", "(la_stack_lu_solve|nalgebra_lu_solve|faer_lu_solve)$"),
        ("dot", "(la_stack_dot|nalgebra_dot|faer_dot)$"),
        ("inf_norm", "(la_stack_inf_norm|nalgebra_inf_norm|faer_inf_norm)$"),
    ],
)
def test_publication_benchmark_command_selects_only_requested_metric(
    metric_name: str,
    expected_filter: str,
) -> None:
    command = criterion_dim_plot._publication_benchmark_command(metric_name)

    assert command[:-2] == criterion_dim_plot._PUBLICATION_BENCHMARK_BASE
    assert command[-2:] == ("--", expected_filter)


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


def _publication_args() -> criterion_dim_plot.PlotCliArgs:
    return criterion_dim_plot.PlotCliArgs(
        metric="lu_solve",
        stat="median",
        sample="new",
        criterion_dir="target/criterion",
        out=None,
        csv=None,
        log_y=True,
        no_plot=False,
        update_readme=True,
        readme="README.md",
        allow_partial=False,
    )


def test_readme_publication_rejects_noncanonical_data_and_asset_paths_before_timing(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    (tmp_path / "README.md").write_text(f"{begin}\nold\n{end}\n", encoding="utf-8")

    assert (
        criterion_dim_plot._validate_readme_target(
            tmp_path,
            replace(_publication_args(), criterion_dir="stale-results"),
        )
        == 2
    )
    assert "requires Criterion output" in capsys.readouterr().err

    assert (
        criterion_dim_plot._validate_readme_target(
            tmp_path,
            replace(_publication_args(), csv="custom.csv"),
        )
        == 2
    )
    assert "canonical CSV/SVG destinations" in capsys.readouterr().err


def test_readme_publication_rejects_no_plot_before_timing(
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


def test_main_update_readme_happy_path(  # noqa: PLR0915
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Create a minimal Criterion directory structure for lu_solve.
    criterion_dir = tmp_path / "target" / "criterion"

    def write_estimates(path: Path, median: float) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {
                    "median": {
                        "point_estimate": median,
                        "confidence_interval": {"lower_bound": median * 0.9, "upper_bound": median * 1.1},
                    }
                }
            ),
            encoding="utf-8",
        )

    for d in criterion_dim_plot.CANONICAL_DIMS:
        la, na, fa = (float(d * 5), float(d * 10), float(d * 20))
        base = criterion_dir / f"d{d}"
        write_estimates(base / "la_stack_lu_solve" / "new" / "estimates.json", la)
        write_estimates(base / "nalgebra_lu_solve" / "new" / "estimates.json", na)
        write_estimates(base / "faer_lu_solve" / "new" / "estimates.json", fa)

    readme = tmp_path / "README.fixture.md"
    marker_begin, marker_end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    readme.write_text(f"# Bench\n{marker_begin}\nplaceholder\n{marker_end}\n", encoding="utf-8")

    (tmp_path / "Cargo.toml").write_text(
        '[package]\nname = "fixture"\nversion = "0.1.0"\n[dev-dependencies]\ncriterion = "0.7"\nnalgebra = "0.34"\nfaer = "0.22"\n',
        encoding="utf-8",
    )
    (tmp_path / "Cargo.lock").write_text("version = 4\n", encoding="utf-8")
    (tmp_path / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.88.0"\n', encoding="utf-8")
    (tmp_path / "justfile").write_text("test-bench-inputs:\n", encoding="utf-8")
    for relative in ("tests/exact_bench_config.rs", "tests/vs_linalg_inputs.rs", "benches/vs_linalg.rs", "src/lib.rs"):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("// fixture\n", encoding="utf-8")

    calls: list[tuple[str, tuple[str, ...]]] = []

    def fake_run_safe(command: str, args: list[str], **_kwargs: object) -> SimpleNamespace:
        calls.append((command, tuple(args)))
        if command == "cargo":
            for dimension in criterion_dim_plot.CANONICAL_DIMS:
                base = criterion_dir / f"d{dimension}"
                write_estimates(base / "la_stack_lu_solve" / "new" / "estimates.json", float(dimension * 5))
                write_estimates(base / "nalgebra_lu_solve" / "new" / "estimates.json", float(dimension * 10))
                write_estimates(base / "faer_lu_solve" / "new" / "estimates.json", float(dimension * 20))
        return SimpleNamespace(stdout="rustc 1.88.0\n" if command == "rustc" else "")

    def fake_run_git(args: list[str], **_kwargs: object) -> SimpleNamespace:
        if "status" in args:
            return SimpleNamespace(stdout=" M src/lib.rs\n")
        return SimpleNamespace(stdout="0123456789abcdef\n")

    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)
    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", fake_run_safe)
    monkeypatch.setattr(criterion_dim_plot, "run_git_command", fake_run_git)

    def fake_render(request: criterion_dim_plot.PlotRequest) -> None:
        request.out_svg.write_text("<svg/>\n", encoding="utf-8")

    monkeypatch.setattr(criterion_dim_plot, "_render_svg_with_gnuplot", fake_render)

    out_csv = tmp_path / "out.csv"
    out_svg = tmp_path / "out.svg"

    rc = criterion_dim_plot.main(
        [
            "--metric",
            "lu_solve",
            "--stat",
            "median",
            "--sample",
            "new",
            "--criterion-dir",
            str(criterion_dir),
            "--csv",
            str(out_csv),
            "--out",
            str(out_svg),
            "--update-readme",
            "--readme",
            str(readme),
        ]
    )
    assert rc == 0
    assert out_svg.read_text(encoding="utf-8") == "<svg/>\n"
    assert calls[:2] == [
        ("just", ("test-bench-inputs",)),
        (
            "cargo",
            (
                "bench",
                "--locked",
                "--features",
                "bench",
                "--bench",
                "vs_linalg",
                "--",
                "(la_stack_lu_solve|nalgebra_lu_solve|faer_lu_solve)$",
            ),
        ),
    ]

    # CSV written
    csv_text = out_csv.read_text(encoding="utf-8")
    assert csv_text.startswith("D,la_stack,la_lo,la_hi,nalgebra,na_lo,na_hi,faer,fa_lo,fa_hi\n")
    assert "2,10.0" in csv_text
    assert "64,320.0" in csv_text

    # README updated with computed table
    readme_text = readme.read_text(encoding="utf-8")
    assert "placeholder" not in readme_text
    assert "| 2 | 10.000 | 20.000 | 40.000 | +50.0% | +75.0% |" in readme_text
    assert "| 64 | 320.000 | 640.000 | 1,280.000 | +50.0% | +75.0% |" in readme_text

    provenance = json.loads(out_csv.with_suffix(".provenance.json").read_text(encoding="utf-8"))
    assert provenance["measurement"]["status"] == "recorded"
    assert provenance["publication"]["correctness_gate"] == "passed"
    assert provenance["publication"]["git_clean"] is False
    assert provenance["criterion"]["benchmark_command"][:3] == ["cargo", "bench", "--locked"]
    assert re.fullmatch(r"[0-9a-f]{64}", provenance["publication"]["source_state_sha256"])


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


def test_main_publication_fails_closed_when_provenance_tool_is_unavailable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    criterion_dir = tmp_path / "target" / "criterion"
    criterion_dir.mkdir(parents=True)
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    readme = tmp_path / "README.fixture.md"
    readme.write_text(f"{begin}\nold table\n{end}\n", encoding="utf-8")
    output = tmp_path / "out.csv"
    row = criterion_dim_plot.Row(2, 1.0, 0.9, 1.1, 2.0, 1.9, 2.1, 3.0, 2.9, 3.1)

    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)
    monkeypatch.setattr(criterion_dim_plot, "_run_publication_benchmarks", lambda _root, _metric: None)
    monkeypatch.setattr(criterion_dim_plot, "_detect_versions", lambda _root: {})
    monkeypatch.setattr(criterion_dim_plot, "_discover_dims", lambda _criterion_dir: [2])
    monkeypatch.setattr(criterion_dim_plot, "_collect_rows", lambda *_args: ([row], []))
    monkeypatch.setattr(
        criterion_dim_plot,
        "_capture_provenance",
        lambda *_args, **_kwargs: {
            "publication": {
                "cargo_lock_sha256": "a" * 64,
                "commit": "commit",
                "cpu": "cpu",
                "git_clean": True,
                "missing_harness_files": [],
                "os": "os",
                "rustc": "unavailable",
                "source_missing": False,
            }
        },
    )

    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(criterion_dir),
            "--csv",
            str(output),
            "--out",
            str(tmp_path / "out.svg"),
            "--update-readme",
            "--readme",
            str(readme),
        ]
    )

    assert rc == 2
    assert "required fields are unavailable: rustc" in capsys.readouterr().err
    assert not output.exists()
    assert "old table" in readme.read_text(encoding="utf-8")


def test_publication_gate_failure_stops_before_timing(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    calls: list[tuple[str, tuple[str, ...]]] = []

    def fail_gate(command: str, args: list[str], **_kwargs: object) -> SimpleNamespace:
        calls.append((command, tuple(args)))
        raise subprocess.CalledProcessError(1, [command, *args], stderr="fixture failure")

    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", fail_gate)

    with pytest.raises(RuntimeError, match="just test-bench-inputs"):
        criterion_dim_plot._run_publication_benchmarks(tmp_path, "lu_solve")

    assert calls == [("just", ("test-bench-inputs",))]


@pytest.mark.parametrize(
    ("failure_kind", "expected_details", "cause_type"),
    [
        ("process", ("timing failed",), subprocess.CalledProcessError),
        ("missing", ("Required executable 'cargo' not found in PATH",), criterion_dim_plot.ExecutableNotFoundError),
        ("timeout", ("timed out after 17 seconds", "timing stalled"), subprocess.TimeoutExpired),
        ("os-error", ("could not start", "working directory unavailable"), OSError),
    ],
)
def test_failed_timing_restores_staged_new_samples(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    failure_kind: str,
    expected_details: tuple[str, ...],
    cause_type: type[Exception],
) -> None:
    old_estimate = tmp_path / "target" / "criterion" / "d2" / "la_stack_lu" / "new" / "estimates.json"
    old_estimate.parent.mkdir(parents=True)
    old_estimate.write_text("old\n", encoding="utf-8")

    def fail_timing(command: str, args: list[str], **_kwargs: object) -> SimpleNamespace:
        if command == "cargo":
            old_estimate.parent.mkdir(parents=True, exist_ok=True)
            old_estimate.write_text("partial\n", encoding="utf-8")
            if failure_kind == "process":
                raise subprocess.CalledProcessError(1, [command, *args], stderr="timing failed")
            if failure_kind == "missing":
                msg = "Required executable 'cargo' not found in PATH"
                raise criterion_dim_plot.ExecutableNotFoundError(msg)
            if failure_kind == "os-error":
                msg = "working directory unavailable"
                raise OSError(msg)
            raise subprocess.TimeoutExpired([command, *args], 17, stderr="timing stalled")
        return SimpleNamespace(stdout="")

    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", fail_timing)

    with pytest.raises(RuntimeError, match="cargo bench") as exc_info:
        criterion_dim_plot._run_publication_benchmarks(tmp_path, "lu_solve")

    assert old_estimate.read_text(encoding="utf-8") == "old\n"
    assert all(detail in str(exc_info.value) for detail in expected_details)
    assert isinstance(exc_info.value.__cause__, cause_type)


def test_partial_staging_failure_restores_every_moved_sample(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    criterion_dir = tmp_path / "target" / "criterion"
    first = criterion_dir / "d2" / "a_bench" / "new" / "estimates.json"
    second = criterion_dir / "d2" / "b_bench" / "new" / "estimates.json"
    for path, text in ((first, "first\n"), (second, "second\n")):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    backup_root = tmp_path / "criterion-backup"
    backup_root.mkdir()
    original_replace = criterion_dim_plot.Path.replace

    def fail_second_move(source: Path, destination: Path) -> Path:
        if source == second.parent and backup_root in destination.parents:
            msg = "simulated second staging failure"
            raise OSError(msg)
        return original_replace(source, destination)

    def fake_mkdtemp(*, prefix: str, **kwargs: object) -> str:
        assert prefix == "la-stack-stale-criterion-"
        assert kwargs == {"dir": criterion_dir.parent}
        return str(backup_root)

    monkeypatch.setattr(criterion_dim_plot.tempfile, "mkdtemp", fake_mkdtemp)
    monkeypatch.setattr(criterion_dim_plot.Path, "replace", fail_second_move)
    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", lambda *_args, **_kwargs: SimpleNamespace(stdout=""))

    with pytest.raises(RuntimeError, match="could not stage existing Criterion samples"):
        criterion_dim_plot._run_publication_benchmarks(tmp_path, "lu_solve")

    assert first.read_text(encoding="utf-8") == "first\n"
    assert second.read_text(encoding="utf-8") == "second\n"
    assert not backup_root.exists()


def test_failed_timing_preserves_backup_when_rollback_fails(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    estimate = tmp_path / "target" / "criterion" / "d2" / "la_stack_lu" / "new" / "estimates.json"
    estimate.parent.mkdir(parents=True)
    estimate.write_text("old\n", encoding="utf-8")
    backup_root = tmp_path / "criterion-backup"
    backup_root.mkdir()

    def fail_timing(command: str, args: list[str], **_kwargs: object) -> SimpleNamespace:
        if command == "cargo":
            estimate.parent.mkdir(parents=True, exist_ok=True)
            estimate.write_text("partial\n", encoding="utf-8")
            raise subprocess.CalledProcessError(1, [command, *args], stderr="timing failed")
        return SimpleNamespace(stdout="")

    original_rmtree = criterion_dim_plot.shutil.rmtree

    def fail_fresh_removal(path: Path) -> None:
        if path == estimate.parent:
            msg = "simulated rollback removal failure"
            raise OSError(msg)
        original_rmtree(path)

    monkeypatch.setattr(criterion_dim_plot.tempfile, "mkdtemp", lambda **_kwargs: str(backup_root))
    monkeypatch.setattr(criterion_dim_plot.shutil, "rmtree", fail_fresh_removal)
    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", fail_timing)

    with pytest.raises(RuntimeError, match="backups preserved") as exc_info:
        criterion_dim_plot._run_publication_benchmarks(tmp_path, "lu_solve")

    assert str(backup_root) in str(exc_info.value)
    preserved = backup_root / "d2" / "la_stack_lu" / "new" / "estimates.json"
    assert preserved.read_text(encoding="utf-8") == "old\n"
    assert estimate.read_text(encoding="utf-8") == "partial\n"


def test_readme_publication_cannot_reuse_stale_new_samples(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    criterion_dir = tmp_path / "target" / "criterion"
    metric = criterion_dim_plot.METRICS["lu_solve"]

    def write_dimension(dimension: int) -> None:
        for bench in (metric.la_bench, metric.na_bench, metric.fa_bench):
            estimates = criterion_dir / f"d{dimension}" / bench / "new" / "estimates.json"
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

    for dimension in criterion_dim_plot.CANONICAL_DIMS:
        write_dimension(dimension)
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    readme = tmp_path / "README.stale-fixture.md"
    readme.write_text(f"{begin}\nstale\n{end}\n", encoding="utf-8")
    output = tmp_path / "benchmark.csv"

    def fake_run_safe(command: str, _args: list[str], **_kwargs: object) -> SimpleNamespace:
        if command == "cargo":
            write_dimension(2)
        return SimpleNamespace(stdout="")

    monkeypatch.setattr(criterion_dim_plot, "_repo_root", lambda: tmp_path)
    monkeypatch.setattr(criterion_dim_plot, "run_safe_command", fake_run_safe)

    rc = criterion_dim_plot.main(
        [
            "--criterion-dir",
            str(criterion_dir),
            "--csv",
            str(output),
            "--update-readme",
            "--readme",
            str(readme),
        ]
    )

    assert rc == 2
    assert not output.exists()
    assert "stale" in readme.read_text(encoding="utf-8")
    assert not (criterion_dir / "d3" / metric.la_bench / "new" / "estimates.json").exists()


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
