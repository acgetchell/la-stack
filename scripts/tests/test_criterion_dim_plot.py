from __future__ import annotations

import argparse
import json
import tomllib
from typing import TYPE_CHECKING, cast

import pytest

import criterion_dim_plot

if TYPE_CHECKING:
    from pathlib import Path


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
            la_lo=0.0,
            la_hi=0.0,
            na_time=100.0,
            na_lo=0.0,
            na_hi=0.0,
            fa_time=200.0,
            fa_lo=0.0,
            fa_hi=0.0,
        ),  # +50.0% vs na, +75.0% vs fa
        criterion_dim_plot.Row(
            dim=64,
            la_time=1_000.0,
            la_lo=0.0,
            la_hi=0.0,
            na_time=900.0,
            na_lo=0.0,
            na_hi=0.0,
            fa_time=800.0,
            fa_lo=0.0,
            fa_hi=0.0,
        ),  # -11.1% vs na, -25.0% vs fa
    ]

    table = criterion_dim_plot._markdown_table(rows, stat="median")

    assert "| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | la-stack vs nalgebra | la-stack vs faer |" in table
    assert "| 2 | 50.000 | 100.000 | 200.000 | +50.0% | +75.0% |" in table
    # thousand separator and sign
    assert "| 64 | 1,000.000 | 900.000 | 800.000 | -11.1% | -25.0% |" in table


def test_markdown_table_handles_zero_nalgebra_time() -> None:
    rows = [
        # nalgebra time of 0 indicates missing/corrupt data; ensure we don't crash.
        criterion_dim_plot.Row(
            dim=2,
            la_time=10.0,
            la_lo=0.0,
            la_hi=0.0,
            na_time=0.0,
            na_lo=0.0,
            na_hi=0.0,
            fa_time=100.0,
            fa_lo=0.0,
            fa_hi=0.0,
        ),
    ]

    table = criterion_dim_plot._markdown_table(rows, stat="median")
    assert "| 2 | 10.000 | 0.000 | 100.000 | n/a | +90.0% |" in table


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
        "\n".join(
            [
                "# Title",
                "before",
                marker_begin,
                "old line 1",
                "old line 2",
                marker_end,
                "after",
                "",
            ]
        ),
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
    readme.write_text("\n".join([marker_end, marker_begin, ""]), encoding="utf-8")

    with pytest.raises(criterion_dim_plot.ReadmeMarkerError, match=r"out of order"):
        criterion_dim_plot._update_readme_table(readme, marker_begin, marker_end, "| x |")


def test_update_readme_table_errors_on_non_unique_markers(tmp_path: Path) -> None:
    marker_begin, marker_end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")

    readme = tmp_path / "README.md"
    readme.write_text(
        "\n".join(
            [
                marker_begin,
                marker_begin,
                marker_end,
                "",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(criterion_dim_plot.ReadmeMarkerError, match=r"not found or not unique"):
        criterion_dim_plot._update_readme_table(readme, marker_begin, marker_end, "| x |")


def test_main_update_readme_no_plot_happy_path(tmp_path: Path) -> None:
    # Create a minimal Criterion directory structure for lu_solve.
    criterion_dir = tmp_path / "criterion"

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

    for d, la, na, fa in [(2, 10.0, 20.0, 40.0), (8, 100.0, 50.0, 200.0)]:
        base = criterion_dir / f"d{d}"
        write_estimates(base / "la_stack_lu_solve" / "new" / "estimates.json", la)
        write_estimates(base / "nalgebra_lu_solve" / "new" / "estimates.json", na)
        write_estimates(base / "faer_lu_solve" / "new" / "estimates.json", fa)

    readme = tmp_path / "README.md"
    marker_begin, marker_end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    readme.write_text("\n".join(["# Bench", marker_begin, "placeholder", marker_end, ""]), encoding="utf-8")

    out_csv = tmp_path / "out.csv"

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
            "--no-plot",
            "--update-readme",
            "--readme",
            str(readme),
        ]
    )
    assert rc == 0

    # CSV written
    csv_text = out_csv.read_text(encoding="utf-8")
    assert csv_text.startswith("D,la_stack,la_lo,la_hi,nalgebra,na_lo,na_hi,faer,fa_lo,fa_hi\n")
    assert "2,10.0" in csv_text
    assert "8,100.0" in csv_text

    # README updated with computed table
    readme_text = readme.read_text(encoding="utf-8")
    assert "placeholder" not in readme_text
    assert "| 2 | 10.000 | 20.000 | 40.000 | +50.0% | +75.0% |" in readme_text
    assert "| 8 | 100.000 | 50.000 | 200.000 | -100.0% | +50.0% |" in readme_text


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
        "\n".join(
            [
                "# comment line",
                "[package]",
                'version = "1.2.3" # inline comment',
                "",
                "[dependencies]",
                'nalgebra = "0.34.0"',
                'faer = { version = "0.21.4" }',
                "",
                "[dev-dependencies]",
                'serde = "1.0"',
            ]
        ),
        encoding="utf-8",
    )

    assert criterion_dim_plot._strip_toml_comment('version = "1.0" # c') == 'version = "1.0"'
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
    assert skipped == ["d2 (missing la_stack_lu_solve, nalgebra_lu_solve, or faer_lu_solve)"]


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
