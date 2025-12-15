from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest

import criterion_dim_plot

if TYPE_CHECKING:
    from pathlib import Path


def test_readme_table_markers_are_stable() -> None:
    begin, end = criterion_dim_plot._readme_table_markers("lu_solve", "median", "new")
    assert begin == "<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->"
    assert end == "<!-- BENCH_TABLE:lu_solve:median:new:END -->"


def test_markdown_table_formats_values_and_pct() -> None:
    rows = [
        # (D, la, la_lo, la_hi, na, na_lo, na_hi)
        (2, 50.0, 0.0, 0.0, 100.0, 0.0, 0.0),  # +50.0%
        (64, 1_000.0, 0.0, 0.0, 900.0, 0.0, 0.0),  # -11.1%
    ]

    table = criterion_dim_plot._markdown_table(rows, stat="median")

    assert "| D | la-stack median (ns) | nalgebra median (ns) | la-stack vs nalgebra |" in table
    assert "| 2 | 50.000 | 100.000 | +50.0% |" in table
    # thousand separator and sign
    assert "| 64 | 1,000.000 | 900.000 | -11.1% |" in table


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

    with pytest.raises(ValueError, match=r"README markers not found"):
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

    with pytest.raises(ValueError, match=r"out of order"):
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

    with pytest.raises(ValueError, match=r"not found or not unique"):
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

    for d, la, na in [(2, 10.0, 20.0), (8, 100.0, 50.0)]:
        base = criterion_dir / f"d{d}"
        write_estimates(base / "la_stack_lu_solve" / "new" / "estimates.json", la)
        write_estimates(base / "nalgebra_lu_solve" / "new" / "estimates.json", na)

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
    assert csv_text.startswith("D,la_stack,la_lo,la_hi,nalgebra,na_lo,na_hi\n")
    assert "2,10.0" in csv_text
    assert "8,100.0" in csv_text

    # README updated with computed table
    readme_text = readme.read_text(encoding="utf-8")
    assert "placeholder" not in readme_text
    assert "| 2 | 10.000 | 20.000 | +50.0% |" in readme_text
    assert "| 8 | 100.000 | 50.000 | -100.0% |" in readme_text
