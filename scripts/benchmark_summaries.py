"""Read historical la-stack summary schemas without rewriting their bytes or hashes."""

import csv
import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from research_repo_tools.criterion import Estimate

from performance_artifacts import ArtifactPaths, TimingEstimate, load_bundle

SCHEMA_VERSION = 1
COLUMNS = (
    "phase",
    "benchmark_id",
    "sample_count",
    "mean_ns",
    "mean_ci_lower_ns",
    "mean_ci_upper_ns",
    "median_ns",
    "median_ci_lower_ns",
    "median_ci_upper_ns",
    "confidence_level",
)


@dataclass(frozen=True, slots=True)
class Measurement:
    """One independently identified baseline or current measurement."""

    phase: Literal["baseline", "current"]
    benchmark_id: str
    mean: TimingEstimate
    median: TimingEstimate

    def __post_init__(self) -> None:
        """Require a phase and a portable, nonempty Criterion identity."""
        if self.phase not in {"baseline", "current"} or not self.benchmark_id.strip():
            msg = "measurement requires a phase and benchmark ID"
            raise ValueError(msg)


def read_object(path: Path) -> dict[str, object]:
    """Read the legacy JSON object without assigning new digest semantics."""
    value = json.loads(path.read_bytes())
    if not isinstance(value, dict):
        raise TypeError(f"expected a legacy JSON object: {path}")
    return cast("dict[str, object]", value)


def full_summary_paths(report: ArtifactPaths) -> ArtifactPaths:
    """Keep the complete summary adjacent to the selected report inputs."""
    return ArtifactPaths(
        csv=report.csv.with_suffix(".full.csv"),
        provenance=report.csv.with_suffix(".full.provenance.json"),
    )


def report_input_paths(report: ArtifactPaths) -> dict[str, Path]:
    """Reserve selected and complete inputs against publication-output aliases."""
    full = full_summary_paths(report)
    return {
        "artifact CSV": report.csv,
        "artifact provenance": report.provenance,
        "complete summary CSV": full.csv,
        "complete summary provenance": full.provenance,
    }


def _validate_rows(rows: list[Measurement] | tuple[Measurement, ...]) -> None:
    if {row.phase for row in rows} != {"baseline", "current"}:
        msg = "complete local summaries require both measurement phases"
        raise ValueError(msg)
    keys = {(row.phase, row.benchmark_id) for row in rows}
    if len(keys) != len(rows):
        msg = "duplicate phase/benchmark identity in local summaries"
        raise ValueError(msg)


def _parse_rows(text: str) -> tuple[Measurement, ...]:
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames != list(COLUMNS):
        msg = "unsupported complete-summary CSV columns"
        raise ValueError(msg)
    rows: list[Measurement] = []
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            msg = "malformed complete-summary CSV row"
            raise ValueError(msg)
        phase = row["phase"]
        if phase not in {"baseline", "current"}:
            raise ValueError(f"invalid measurement phase: {phase!r}")
        if row["sample_count"] != "100" or row["confidence_level"] != "0.95":
            msg = "local summaries require full release sampling and 95% intervals"
            raise ValueError(msg)
        estimates = [
            TimingEstimate(*(Estimate(float(row[field])).point for field in fields))
            for fields in (
                ("mean_ns", "mean_ci_lower_ns", "mean_ci_upper_ns"),
                ("median_ns", "median_ci_lower_ns", "median_ci_upper_ns"),
            )
        ]
        rows.append(Measurement(phase, row["benchmark_id"], *estimates))
    _validate_rows(rows)
    return tuple(rows)


def _validate_saved_summary(report: ArtifactPaths) -> None:
    full = full_summary_paths(report)
    payload = full.csv.read_bytes()
    rows = _parse_rows(payload.decode("utf-8"))
    provenance = read_object(full.provenance)
    expected = {
        "schema_version": SCHEMA_VERSION,
        "columns": list(COLUMNS),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "counts": {phase: sum(row.phase == phase for row in rows) for phase in ("baseline", "current")},
    }
    measurements = provenance.pop("measurements", None)
    if not isinstance(measurements, dict) or type(measurements.get("schema_version")) is not int or measurements != expected:
        msg = "complete-summary digest, schema, or phase counts do not match"
        raise ValueError(msg)
    if provenance != read_object(report.provenance):
        msg = "complete summaries and selected report have different provenance"
        raise ValueError(msg)
    load_bundle(report)


def resolve_report_paths(root: Path, report: ArtifactPaths) -> ArtifactPaths:
    """Use committed summaries after clean, without hiding partial scratch artifacts."""
    default = root / "target/bench-reports/performance.csv"
    if (
        report.csv.resolve() != default.resolve()
        or report.provenance.resolve() != default.with_suffix(".provenance.json").resolve()
        or any(path.exists() or path.is_symlink() for path in report_input_paths(report).values())
    ):
        return report
    directory = root / "docs/performance"
    pointer = directory / "latest.json"
    if not pointer.exists():
        return report
    metadata = read_object(pointer)
    relative = metadata.get("run")
    if type(metadata.get("schema_version")) is not int or metadata.get("schema_version") != SCHEMA_VERSION or not isinstance(relative, str):
        msg = "invalid persistent benchmark pointer"
        raise ValueError(msg)
    run = (directory / relative).resolve()
    if not run.is_relative_to(directory.resolve()) or len(Path(relative).parts) != 2:
        msg = "persistent benchmark pointer escapes its run directory"
        raise ValueError(msg)
    retained = ArtifactPaths(csv=run / "performance.csv", provenance=run / "performance.provenance.json")
    _validate_saved_summary(retained)
    return retained
