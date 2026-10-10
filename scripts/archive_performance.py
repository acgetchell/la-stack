#!/usr/bin/env -S uv run --locked
"""Consumer release selection and publication over shared performance workflows."""

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import dataclass, replace
from pathlib import Path

from research_repo_tools.archives import extract_archive
from research_repo_tools.common_measurement import measure_prepared_pair
from research_repo_tools.complete_runs import render_run, run_identity
from research_repo_tools.evidence import Evidence, serialize_evidence
from research_repo_tools.files import replace_many
from research_repo_tools.measurement import resolve_revision
from research_repo_tools.process import ExecutableNotFoundError, format_exception_diagnostics, run_command
from research_repo_tools.publication import plan_outputs, publish_publication
from research_repo_tools.release_assets import download_release_asset
from research_repo_tools.release_discovery import normalize_tag, stable_published_releases
from research_repo_tools.release_pairs import PairMode, ReleasePair, resolve_pair
from research_repo_tools.run_reports import load_run_report_plan
from research_repo_tools.worktrees import apply_snapshot, capture_snapshot, temporary_worktree

from bench_compare import render_historical_assets, render_release_artifacts
from benchmark_contract import find_project_root
from benchmark_summaries import resolve_report_paths
from performance_artifacts import ArtifactPaths, load_bundle, resolve_shared_harness_compatibility
from performance_runs import (
    REPORT_CONFIG,
    archive_directory,
    measurement_plan,
    retained_run,
    scratch_paths,
    validate_scientific_run,
)


@dataclass(frozen=True)
class ReportId:
    """Legacy release labels, retained without changing historical artifact framing."""

    current_tag: str
    baseline_tag: str


def select_pair(root: Path, mode: PairMode, current: str | None, baseline: str | None, *, local: bool = False) -> ReleasePair:
    """Use publication chronology while retaining la-stack release eligibility."""
    package = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    releases = ()
    if mode != "explicit":
        response = run_command(
            "gh",
            ["release", "list", "--repo", "acgetchell/la-stack", "--limit", "1000", "--json", "tagName,isDraft,isPrerelease,publishedAt"],
            cwd=root,
        )
        releases = stable_published_releases(json.loads(response.stdout))
    # Shared explicit selection describes release pairs. Consumer local audits
    # additionally permit two different sources carrying the same package label.
    if local and current is not None and baseline is not None and normalize_tag(current) == normalize_tag(baseline):
        pair = ReleasePair(current, baseline)
    else:
        pair = resolve_pair(mode, package_tag=package, releases=releases, current=current, baseline=baseline, order="published")
    resolve_shared_harness_compatibility(current=pair.current, baseline=pair.baseline, shared_harness_rational_inputs=True)
    return pair


def measure(root: Path, pair: ReleasePair, suite: str, scope: str) -> Evidence:
    """Capture tracked changes and delegate isolated phases to the published engine."""
    version = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    if pair.current != normalize_tag(version):
        msg = "current release label must match the measured checkout package"
        raise ValueError(msg)
    snapshot = replace(capture_snapshot(root), untracked=())
    # The existing consumer contract deliberately excludes untracked files.
    run_command("git", ["--no-pager", "fetch", "origin", f"refs/tags/{pair.baseline}:refs/tags/{pair.baseline}"], cwd=root)
    baseline_revision = resolve_revision(root, pair.baseline)
    parent = Path(tempfile.mkdtemp(prefix="la-stack-performance-"))
    try:
        with (
            temporary_worktree(root, parent / "baseline", baseline_revision, allow_git_mutations=True) as baseline,
            temporary_worktree(root, parent / "current", snapshot.revision, allow_git_mutations=True) as current,
        ):
            apply_snapshot(current, snapshot)
            # Consumer API-layout adapter: the comparison test moved to its own
            # workspace member and the historical root copy must not compile.
            (baseline / "tests/vs_linalg_inputs.rs").unlink(missing_ok=True)
            plan = measurement_plan(current, pair, suite, scope, parent / "listing")
            evidence = measure_prepared_pair(baseline, current, plan, pair, working_tree=True)
            if replace(capture_snapshot(root), untracked=()) != snapshot:
                msg = "source working tree changed while measurement ran"
                raise ValueError(msg)
    finally:
        # A failed shared worktree cleanup leaves its path intact for recovery.
        if not (parent / "baseline").exists() and not (parent / "current").exists():
            shutil.rmtree(parent)
    validate_scientific_run(evidence)
    return evidence


def save_scratch(root: Path, stem: Path, evidence: Evidence) -> None:
    """Save complete evidence and the shared report without touching legacy files."""
    validate_scientific_run(evidence)
    payload, manifest = scratch_paths(stem)
    data, envelope = serialize_evidence(evidence)
    report = render_run(evidence)
    plan = plan_outputs(
        root,
        {path.relative_to(root).as_posix(): value for path, value in ((payload, data), (manifest, envelope), (stem.with_suffix(".md"), report))},
        inputs={REPORT_CONFIG: (root / REPORT_CONFIG).read_bytes()},
    )
    publish_publication(plan)


def promote(root: Path, stem: Path) -> None:
    """Validate consumer policy and publish the shared report and history plan."""
    evidence = retained_run(root, stem)
    validate_scientific_run(evidence, release=True)
    payload, manifest = scratch_paths(stem)
    kwargs = {}
    if payload.exists() or manifest.exists():
        kwargs = {"payload": payload.relative_to(root).as_posix(), "manifest": manifest.relative_to(root).as_posix()}
    history = load_run_report_plan(root, REPORT_CONFIG, **kwargs)
    outputs = dict(history.outputs)
    if json.loads(outputs[f"{archive_directory(root)}/latest.json"])["run"]["id"] != run_identity(evidence):
        msg = "selected evidence changed while planning report publication"
        raise ValueError(msg)
    publish_publication(history)


def parse_report_id(text: str) -> ReportId:
    """Read the established scientific report header without modifying its bytes."""
    current = re.findall(r"^\*\*la-stack\*\* (v[^\s]+)", text, re.MULTILINE)
    baseline = re.findall(r"^Comparison against baseline \*\*([^*]+)\*\*:", text, re.MULTILINE)
    if len(current) != 1 or len(baseline) != 1:
        msg = "performance report must contain exactly one current and baseline release identity"
        raise ValueError(msg)
    return ReportId(normalize_tag(current[0]), normalize_tag(baseline[0]))


def render_and_promote_artifacts(*, artifacts: ArtifactPaths, output: Path, current: Path, archive_dir: Path) -> ReportId:
    """Read historical CSV/JSON through its original schema; never relabel hashes."""
    bundle = load_bundle(artifacts)
    release = bundle.context.release
    if release.current == release.baseline:
        msg = "cannot promote a same-version local comparison as release documentation"
        raise ValueError(msg)
    text = render_release_artifacts(artifacts).encode()
    root = current.parent.parent.resolve()
    outputs = {current.relative_to(root).as_posix(): text, output.relative_to(root).as_posix(): text}
    if current.exists() and current.read_bytes() != text:
        identity = parse_report_id(current.read_text(encoding="utf-8"))
        if identity != ReportId(release.current, release.baseline):
            archived = archive_dir / f"{identity.current_tag}-vs-{identity.baseline_tag}.md"
            outputs[archived.relative_to(root).as_posix()] = current.read_bytes()
    inputs = {path.relative_to(root).as_posix(): path.read_bytes() for path in (artifacts.csv, artifacts.provenance)}
    publish_publication(plan_outputs(root, outputs, inputs=inputs, immutable=tuple(name for name in outputs if name.startswith("docs/archive/"))))
    return ReportId(release.current, release.baseline)


def compare_assets(root: Path, pair: ReleasePair, output: Path, suite: str, scope: str) -> None:
    """Read published raw archives with shared bounded download/extraction APIs."""
    with tempfile.TemporaryDirectory(prefix="la-stack-assets-") as directory:
        parent = Path(directory)
        criterion = parent / "criterion"
        for phase, tag in (("baseline", pair.baseline), ("current", pair.current)):
            asset = parent / f"{tag}.tar.gz"
            download_release_asset(root, "acgetchell/la-stack", tag, f"la-stack-{tag}-criterion-baseline.tar.gz", asset)
            extracted = parent / phase
            extract_archive(asset, extracted)
            # Preserve native measurements, explicitly retaining their original
            # per-release harnesses; these cannot become a complete common run.
            sample = tag
            destination_sample = "new" if phase == "current" else tag
            for source in (extracted / "criterion").glob(f"**/{sample}/*.json"):
                target = criterion / source.relative_to(extracted / "criterion").parent.parent / destination_sample / source.name
                replace_many({target: source.read_bytes()})
        text = render_historical_assets(criterion, pair, suite, scope)
        publish_publication(plan_outputs(root, {output.relative_to(root).as_posix(): text}, inputs={"Cargo.toml": (root / "Cargo.toml").read_bytes()}))


def main(argv: list[str] | None = None) -> int:
    """Keep routine performance commands thin and explicit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("local", "release", "report", "assets"))
    parser.add_argument("current", nargs="?")
    parser.add_argument("baseline", nargs="?")
    parser.add_argument("--suite", choices=("all", "exact", "vs_linalg"), default="all")
    parser.add_argument("--scope", choices=("release-signal", "all-benches"), default="release-signal")
    parser.add_argument("--stem", type=Path, default=Path("target/bench-reports/performance"))
    args = parser.parse_args(argv)
    root = find_project_root(Path.cwd()).resolve()
    stem = root / args.stem
    try:
        if args.command == "report":
            if args.current or args.baseline:
                parser.error("report reads retained evidence and takes no release labels")
            if any(path.exists() or path.is_symlink() for path in (*scratch_paths(stem), root / archive_directory(root))):
                promote(root, stem)
            else:
                paths = resolve_report_paths(root, ArtifactPaths(stem.with_suffix(".csv"), stem.with_suffix(".provenance.json")))
                render_and_promote_artifacts(
                    artifacts=paths, output=stem.with_suffix(".md"), current=root / "docs/performance.md", archive_dir=root / "docs/archive/performance"
                )
            return 0
        modes: dict[str, PairMode] = {"local": "current-vs-latest", "release": "infer-release", "assets": "published-latest"}
        mode: PairMode = "explicit" if args.current or args.baseline else modes[args.command]
        pair = select_pair(root, mode, args.current or None, args.baseline or None, local=args.command == "local")
        if args.command == "assets":
            compare_assets(root, pair, root / "target/bench-reports/github-assets-performance.md", args.suite, args.scope)
        else:
            if args.command == "release" and pair.current == pair.baseline:
                parser.error("release requires distinct release labels")
            evidence = measure(root, pair, args.suite, args.scope)
            save_scratch(root, stem, evidence)
            if args.command == "release":
                promote(root, stem)
        return 0
    except (OSError, ValueError, TypeError, KeyError, ExecutableNotFoundError, subprocess.SubprocessError, ExceptionGroup) as error:
        print(f"performance: {format_exception_diagnostics(error)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
