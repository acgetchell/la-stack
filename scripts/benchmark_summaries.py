"""Preserve every recorded local benchmark summary outside disposable targets."""

import csv
import hashlib
import io
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from criterion_measurements import positive_number, read_object, validate_measurement
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


def _estimate(directory: Path, statistic: str) -> TimingEstimate:
    estimates = read_object(directory / "estimates.json")
    estimate = estimates[statistic]
    if not isinstance(estimate, dict) or not isinstance(interval := estimate.get("confidence_interval"), dict):
        raise TypeError(f"missing {statistic} confidence interval in {directory}")
    return TimingEstimate(
        median_ns=positive_number(estimate.get("point_estimate")),
        ci_lower_ns=positive_number(interval.get("lower_bound")),
        ci_upper_ns=positive_number(interval.get("upper_bound")),
    )


def collect_measurements(criterion: Path, baseline: str) -> tuple[Measurement, ...]:
    """Read all recorded cases, including peers and groups outside the report registry.

    The caller must remove baseline ``new`` directories before current timing.
    Named baseline samples survive that cleanup. Their phase is never inferred
    from the benchmark's library name or substituted for an absent current run.
    """
    rows: list[Measurement] = []
    for phase, sample in (("baseline", baseline), ("current", "new")):
        for directory in sorted(path for path in criterion.rglob(sample) if path.is_dir()):
            metadata = read_object(directory / "benchmark.json")
            benchmark_id = metadata.get("full_id")
            if not isinstance(benchmark_id, str) or not benchmark_id.strip():
                raise ValueError(f"missing Criterion full_id in {directory}")
            validate_measurement(directory)
            rows.append(Measurement(phase, benchmark_id, _estimate(directory, "mean"), _estimate(directory, "median")))
    _validate_rows(rows)
    return tuple(sorted(rows, key=lambda row: (row.phase, row.benchmark_id)))


def _validate_rows(rows: list[Measurement] | tuple[Measurement, ...]) -> None:
    if {row.phase for row in rows} != {"baseline", "current"}:
        msg = "complete local summaries require both measurement phases"
        raise ValueError(msg)
    keys = {(row.phase, row.benchmark_id) for row in rows}
    if len(keys) != len(rows):
        msg = "duplicate phase/benchmark identity in local summaries"
        raise ValueError(msg)


def _serialize_rows(rows: tuple[Measurement, ...]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(COLUMNS)
    for row in rows:
        writer.writerow(
            (
                row.phase,
                row.benchmark_id,
                100,
                row.mean.median_ns,
                row.mean.ci_lower_ns,
                row.mean.ci_upper_ns,
                row.median.median_ns,
                row.median.ci_lower_ns,
                row.median.ci_upper_ns,
                0.95,
            )
        )
    return output.getvalue()


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
            TimingEstimate(*(positive_number(float(row[field])) for field in fields))
            for fields in (
                ("mean_ns", "mean_ci_lower_ns", "mean_ci_upper_ns"),
                ("median_ns", "median_ci_lower_ns", "median_ci_upper_ns"),
            )
        ]
        rows.append(Measurement(phase, row["benchmark_id"], *estimates))
    _validate_rows(rows)
    return tuple(rows)


def summary_outputs(criterion: Path, baseline: str, report: ArtifactPaths) -> dict[Path, str]:
    """Validate every measured case and bind its CSV to the selected report provenance."""
    bundle = load_bundle(report)
    if bundle.context.benchmark_provenance["mode"] != "shared-current-harness":
        msg = "local summaries require locally measured shared-harness results"
        raise ValueError(msg)
    rows = collect_measurements(criterion, baseline)
    payload = _serialize_rows(rows)
    if _parse_rows(payload) != rows:
        msg = "local summary CSV did not preserve timing precision"
        raise ValueError(msg)
    provenance = read_object(report.provenance)
    provenance["measurements"] = {
        "schema_version": SCHEMA_VERSION,
        "columns": COLUMNS,
        "sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "counts": {phase: sum(row.phase == phase for row in rows) for phase in ("baseline", "current")},
    }
    paths = full_summary_paths(report)
    return {paths.csv: payload, paths.provenance: json.dumps(provenance, indent=2, sort_keys=True) + "\n"}


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


def retained_outputs(performance_dir: Path, report: ArtifactPaths) -> dict[Path, str]:
    """Build an immutable run directory and latest pointer for a release promotion.

    Legacy selected-only artifacts remain renderable, but cannot be presented
    as complete histories or replace the latest complete-summary pointer.
    Include unchanged destinations so promotion can protect every retained file
    before deciding which payloads need writing.
    """
    full = full_summary_paths(report)
    if not full.csv.exists() and not full.provenance.exists():
        return {}
    _validate_saved_summary(report)
    bundle = load_bundle(report)
    release = bundle.context.release
    if any(re.fullmatch(r"v[0-9A-Za-z.+-]+", tag) is None for tag in (release.current, release.baseline)):
        msg = "invalid release identifier for persistent summaries"
        raise ValueError(msg)
    payloads = {
        "performance.csv": report.csv.read_text(encoding="utf-8"),
        "performance.provenance.json": report.provenance.read_text(encoding="utf-8"),
        "performance.full.csv": full.csv.read_text(encoding="utf-8"),
        "performance.full.provenance.json": full.provenance.read_text(encoding="utf-8"),
    }
    digest = hashlib.sha256(json.dumps(payloads, sort_keys=True).encode()).hexdigest()
    relative = Path(f"{release.current}-vs-{release.baseline}") / digest
    if not (performance_dir / relative).resolve().is_relative_to(performance_dir.resolve()):
        msg = "benchmark snapshot escapes the performance directory"
        raise ValueError(msg)
    outputs = {performance_dir / relative / name: text for name, text in payloads.items()}
    for path, text in outputs.items():
        if path.exists() and path.read_text(encoding="utf-8") != text:
            raise ValueError(f"existing benchmark snapshot differs: {path}")
    outputs[performance_dir / "latest.json"] = json.dumps({"schema_version": SCHEMA_VERSION, "run": relative.as_posix()}, indent=2) + "\n"
    return outputs


def resolve_report_paths(root: Path, report: ArtifactPaths) -> ArtifactPaths:
    """Use committed summaries after clean, without hiding partial scratch artifacts."""
    default = root / "target/bench-reports/performance.csv"
    if (
        report.csv.resolve() != default.resolve()
        or report.provenance.resolve() != default.with_suffix(".provenance.json").resolve()
        or any(path.exists() for path in report_input_paths(report).values())
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
