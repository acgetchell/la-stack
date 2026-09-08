"""Validate raw Criterion measurements shared by local and hosted retention."""

import json
import math
from typing import TYPE_CHECKING, cast

from performance_artifacts import TimingEstimate

if TYPE_CHECKING:
    from pathlib import Path


def read_object(path: Path) -> dict[str, object]:
    """Require a JSON object at the raw artifact boundary."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError(f"expected a JSON object in {path}")
    return cast("dict[str, object]", data)


def positive_number(value: object) -> float:
    """Reject booleans, nonnumeric values, and invalid timing numbers."""
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"expected a finite positive timing value, got {value!r}")
    return float(value)


def validate_measurement(directory: Path) -> None:
    """Check full default sampling and the estimates used by report consumers."""
    samples = read_object(directory / "sample.json")
    for field in ("iters", "times"):
        values = samples.get(field)
        if not isinstance(values, list) or len(values) != 100:
            raise ValueError(f"{directory}: {field} must contain 100 samples")
        for value in values:
            positive_number(value)
    estimates = read_object(directory / "estimates.json")
    for statistic in ("mean", "median"):
        estimate = estimates.get(statistic)
        if not isinstance(estimate, dict) or not isinstance(interval := estimate.get("confidence_interval"), dict):
            raise TypeError(f"{directory}: missing {statistic} estimate or confidence interval")
        if interval.get("confidence_level") != 0.95:
            raise ValueError(f"{directory}: expected a 95% confidence interval")
        TimingEstimate(
            median_ns=positive_number(estimate.get("point_estimate")),
            ci_lower_ns=positive_number(interval.get("lower_bound")),
            ci_upper_ns=positive_number(interval.get("upper_bound")),
        )
