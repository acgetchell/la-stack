"""Validate raw Criterion measurements shared by local and hosted retention."""

import json
from typing import TYPE_CHECKING, cast

from research_repo_tools.criterion import Estimate, parse_estimate

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
    return Estimate(cast("float", value)).point


def validate_measurement(directory: Path) -> None:
    """Check full default sampling and the estimates used by report consumers."""
    samples = read_object(directory / "sample.json")
    for field in ("iters", "times"):
        values = samples.get(field)
        if not isinstance(values, list) or len(values) != 100:
            raise ValueError(f"{directory}: {field} must contain 100 samples")
        for index, value in enumerate(values):
            try:
                positive_number(value)
            except ValueError as exc:
                raise ValueError(f"{directory / 'sample.json'}: {field}[{index}]: {exc}") from exc
    payload = (directory / "estimates.json").read_bytes()
    for statistic in ("mean", "median"):
        try:
            estimate = parse_estimate(payload, statistic=statistic)
        except ValueError as exc:
            raise ValueError(f"{directory / 'estimates.json'}: {statistic}: {exc}") from exc
        if estimate.lower is None or estimate.upper is None or estimate.confidence_level != 0.95:
            raise ValueError(f"{directory}: expected a complete 95% confidence interval for {statistic}")
