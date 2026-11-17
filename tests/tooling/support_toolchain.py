"""Declared supplier selection shared by installed-runtime consumer fixtures."""

import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from research_repo_tools.dependency_provenance import check_development_dependency

REPO_ROOT = Path(__file__).resolve().parents[2]


def tooling_requirement(root: Path) -> str:
    """Select the declared supplier after RRT verifies any development source."""
    document = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    (raw,) = document["dependency-groups"]["tooling"]
    requirement = Requirement(raw)
    (pin,) = requirement.specifier
    assert requirement.name == "research-repo-tools"
    assert requirement.url is None
    assert requirement.marker is None
    assert pin.operator == "=="
    assert "*" not in pin.version
    development = check_development_dependency(root, document=document)
    if development is not None:
        return f"{requirement.name} @ git+{development.source}@{development.revision}"
    return f"{requirement.name}{requirement.specifier}"


TOOLING_REQUIREMENT = tooling_requirement(REPO_ROOT)
