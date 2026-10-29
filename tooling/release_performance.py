"""Capture release inputs during preparation and reuse them for README publication.

This temporary adapter composes public research-repo-tools APIs. Release
selection and measurement/publication implementations remain upstream.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import cast

from research_repo_tools.cli import main as tooling_main
from research_repo_tools.common_measurement import CompleteMeasurementPlan
from research_repo_tools.files import replace_many
from research_repo_tools.measurement import load_measurement
from research_repo_tools.process import ExecutableNotFoundError, format_exception_diagnostics, run_git_command
from research_repo_tools.release_metadata import read_package_info
from research_repo_tools.release_pairs import ReleasePair, resolve_pair
from research_repo_tools.releases import published_releases

INPUTS = Path("target/bench-reports/release-inputs.json")


def capture(root: Path, current: str, baseline: str) -> tuple[ReleasePair, str, str]:
    """Resolve independent source assertions with the same policy as measurement."""
    if bool(current) != bool(baseline):
        msg = "current and baseline tags must be provided together"
        raise ValueError(msg)
    package_tag = "v" + read_package_info(root).version
    releases = () if current else published_releases(root, repository="acgetchell/la-stack")
    pair = resolve_pair(
        "explicit" if current else "infer-release",
        package_tag=package_tag,
        releases=releases,
        order="published",
        current=current or None,
        baseline=baseline or None,
    )
    configuration = load_measurement(root, "tooling/measurements/full.toml")
    policy = configuration.measurement if isinstance(configuration, CompleteMeasurementPlan) else configuration
    policy.validate_release(pair.current)
    policy.validate_release(pair.baseline)
    tagged = bool(current) or any(item.tag == pair.current for item in releases)

    def commit(reference: str) -> str:
        return run_git_command(["--no-pager", "rev-parse", "--verify", f"{reference}^{{commit}}"], cwd=root, timeout=30).stdout.strip()

    return pair, commit(f"refs/tags/{pair.current}" if tagged else "HEAD"), commit(f"refs/tags/{pair.baseline}")


def prepare(root: Path, current: str, baseline: str) -> None:
    """Save assertions atomically after all release and commit lookups succeed."""
    pair, current_revision, baseline_revision = capture(root, current, baseline)
    assertions = [pair.current, current_revision, pair.baseline, baseline_revision]
    replace_many({root / INPUTS: (json.dumps({"schema": 1, "assertions": assertions}, sort_keys=True) + "\n").encode("utf-8")})
    print(f"Prepared {pair.current} ({current_revision}) vs {pair.baseline} ({baseline_revision})")
    print(f"Saved release tags and commits to {INPUTS.as_posix()} (relative to the repository root; ignored by Git).")
    print("No shell variables or exports are needed. This command does not set PREVIOUS_TAG, CURRENT_COMMIT, or PREVIOUS_COMMIT in your shell.")
    print("The saved file survives opening a new terminal. After measurements finish, run: just performance-readme (no arguments).")
    print("performance-readme reads the release tags and commits from this saved file.")


def saved_assertions(root: Path) -> list[str]:
    """Require a complete preparation record, independent of measured evidence."""
    document = json.loads((root / INPUTS).read_text(encoding="utf-8"))
    if not isinstance(document, dict) or set(document) != {"schema", "assertions"} or type(document["schema"]) is not int or document["schema"] != 1:
        msg = "invalid release preparation record"
        raise ValueError(msg)
    assertions = document["assertions"]
    if not isinstance(assertions, list) or len(assertions) != 4 or any(not isinstance(value, str) or not value for value in assertions):
        msg = "release preparation requires four nonempty tag/commit assertions"
        raise ValueError(msg)
    pair = ReleasePair(assertions[0], assertions[2])
    if pair.current == pair.baseline or any(re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", assertions[index]) is None for index in (1, 3)):
        msg = "release preparation requires distinct tags and full commit identities"
        raise ValueError(msg)
    return cast("list[str]", assertions)


def publish(root: Path, assertions: list[str]) -> int:
    """Publish retained evidence only; never acquire or promote measurements."""
    if any(assertions) and not all(assertions):
        msg = "provide all four tag/commit assertions, or omit all four to use just performance-prepare"
        raise ValueError(msg)
    current, current_revision, baseline, baseline_revision = assertions if all(assertions) else saved_assertions(root)
    return tooling_main(
        [
            "performance",
            "publish",
            "tooling/performance-readme.toml",
            "--current-release",
            current,
            "--current-revision",
            current_revision,
            "--baseline-release",
            baseline,
            "--baseline-revision",
            baseline_revision,
        ]
    )


def main(argv: list[str] | None = None) -> int:
    """Keep preparation and retryable README publication separate from measurements."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    preparation = commands.add_parser("capture")
    preparation.add_argument("current", nargs="?", default="")
    preparation.add_argument("baseline", nargs="?", default="")
    publication = commands.add_parser("publish")
    for name in ("current", "current_revision", "baseline", "baseline_revision"):
        publication.add_argument(name, nargs="?", default="")
    args = parser.parse_args(argv)
    try:
        if args.command == "capture":
            prepare(Path.cwd(), args.current, args.baseline)
            return 0
        return publish(Path.cwd(), [args.current, args.current_revision, args.baseline, args.baseline_revision])
    except (ExecutableNotFoundError, KeyError, OSError, TypeError, ValueError, subprocess.SubprocessError) as error:
        print(f"Release performance rejected: {format_exception_diagnostics(error, single_line=True)}", file=sys.stderr)
        print("Run just performance-prepare to capture inputs, or supply all four assertions to just performance-readme.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
