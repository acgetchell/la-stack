#!/usr/bin/env -S uv run --locked
"""Consumer release selection and publication over shared performance workflows."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import replace
from pathlib import Path
from typing import Literal

from research_repo_tools.archives import extract_archive
from research_repo_tools.common_measurement import measure_prepared_pair
from research_repo_tools.complete_runs import render_run, run_identity
from research_repo_tools.evidence import Evidence, serialize_evidence
from research_repo_tools.measurement import resolve_revision
from research_repo_tools.postprocess_changelog import MarkdownFormatError, format_markdown
from research_repo_tools.process import ExecutableNotFoundError, format_exception_diagnostics, run_command
from research_repo_tools.publication import GitCheck, GlobCheck, MarkerPair, RunSelectionCheck, plan_outputs, publish_publication, replace_section
from research_repo_tools.release_assets import download_release_asset
from research_repo_tools.release_discovery import normalize_tag, stable_published_releases
from research_repo_tools.release_pairs import PairMode, ReleasePair, resolve_pair
from research_repo_tools.run_presentation import (
    CoordinatePlot,
    RunLayout,
    RunRow,
    RunTableStyle,
    SeriesSelection,
    render_coordinate_svg,
    render_run_table,
    serialize_run_csv,
)
from research_repo_tools.run_reports import load_latest_run, load_run_report_plan
from research_repo_tools.worktrees import apply_snapshot, capture_snapshot, temporary_worktree

from performance_runs import (
    REPORT_CONFIG,
    archive_directory,
    compatibility_adapter,
    find_project_root,
    measurement_plan,
    required_report_ids,
    retained_run,
    scientific_policy,
    scratch_paths,
    validate_scientific_run,
)


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
    compatibility_adapter(pair, root)
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
    parent = Path(tempfile.mkdtemp(prefix="la-stack-performance-")).resolve(strict=True)
    try:
        with (
            temporary_worktree(root, parent / "baseline", baseline_revision, allow_git_mutations=True) as baseline,
            temporary_worktree(root, parent / "current", snapshot.revision, allow_git_mutations=True) as current,
        ):
            apply_snapshot(current, snapshot)
            # Consumer API-layout adapter: the comparison test moved to its own
            # workspace member and the historical root copy must not compile.
            (baseline / "tests/vs_linalg_inputs.rs").unlink(missing_ok=True)
            plan = measurement_plan(current, pair, suite, scope)
            evidence = measure_prepared_pair(baseline, current, plan, pair, working_tree=True)
            if replace(capture_snapshot(root), untracked=()) != snapshot:
                msg = "source working tree changed while measurement ran"
                raise ValueError(msg)
    finally:
        # A failed shared worktree cleanup leaves its path intact for recovery.
        if not (parent / "baseline").exists() and not (parent / "current").exists():
            shutil.rmtree(parent)
    validate_scientific_run(evidence, root=root)
    return evidence


def save_scratch(root: Path, stem: Path, evidence: Evidence) -> None:
    """Save complete evidence and the shared report without touching legacy files."""
    validate_scientific_run(evidence, root=root)
    payload, manifest = scratch_paths(stem)
    data, envelope = serialize_evidence(evidence)
    report = render_run(
        evidence,
        report_path=stem.with_suffix(".md").relative_to(root).as_posix(),
        payload_path=payload.relative_to(root).as_posix(),
        manifest_path=manifest.relative_to(root).as_posix(),
    )
    plan = plan_outputs(
        root,
        {path.relative_to(root).as_posix(): value for path, value in ((payload, data), (manifest, envelope), (stem.with_suffix(".md"), report))},
        inputs={REPORT_CONFIG: (root / REPORT_CONFIG).read_bytes()},
    )
    publish_publication(plan)


def promote(root: Path, stem: Path) -> None:
    """Validate consumer policy and publish the shared report and history plan."""
    evidence = retained_run(root, stem)
    validate_scientific_run(evidence, release=True, root=root)
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


def compare_assets(root: Path, pair: ReleasePair, output: Path) -> None:
    """Compare raw historical archives without inventing missing provenance."""
    from research_repo_tools.criterion import collect_sample, compare_samples, render_comparison  # noqa: PLC0415

    samples = {}
    with tempfile.TemporaryDirectory(prefix="la-stack-assets-") as directory:
        parent = Path(directory)
        for phase, tag in (("baseline", pair.baseline), ("current", pair.current)):
            asset = parent / f"{tag}.tar.gz"
            download_release_asset(root, "acgetchell/la-stack", tag, f"la-stack-{tag}-criterion-baseline.tar.gz", asset)
            extracted = parent / phase
            extract_archive(asset, extracted)
            samples[phase] = collect_sample(extracted / "criterion", tag)
        comparison = compare_samples(samples["baseline"], samples["current"])
        if not comparison.comparisons:
            msg = "release archives have no common benchmark rows"
            raise ValueError(msg)
        text = (
            f"# Published release timings\n\nBaseline: {pair.baseline}; current: {pair.current}.\n\n"
            "Each archive used its original release harness. Source, host and gate provenance are not recorded; "
            "ratios are descriptive and require independent hardware/workload review.\n\n" + render_comparison(comparison)
        ).encode()
        publish_publication(plan_outputs(root, {output.relative_to(root).as_posix(): text}, inputs={"Cargo.toml": (root / "Cargo.toml").read_bytes()}))


def check_release_sample(root: Path, name: str) -> int:
    """Require canonical report identities and complete, matching saved bytes."""
    from research_repo_tools.criterion_baselines import load_saved_sample, validate_saved_sample  # noqa: PLC0415

    plan = load_saved_sample(root, "tooling/criterion-saved.toml")
    required = required_report_ids(root)
    if not required <= set(plan.policy.expected):
        msg = "release inventory omits report consumers"
        raise ValueError(msg)
    saved = validate_saved_sample(Path(os.environ.get("CRITERION_HOME", str(root / plan.criterion_dir))), name, plan.policy, raw_files=plan.raw_files)
    for case in saved.cases:
        if case.full_id in required and case.storage_path != case.full_id:
            msg = f"non-canonical report benchmark storage path for {case.full_id}: {case.storage_path}"
            raise ValueError(msg)
    return len(saved.cases)


def main(argv: list[str] | None = None) -> int:
    """Keep routine performance commands thin and explicit."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("local", "release", "assets"):
        command = commands.add_parser(name)
        command.add_argument("current", nargs="?")
        command.add_argument("baseline", nargs="?")
        if name != "assets":
            command.add_argument("--suite", choices=("all", "exact", "vs_linalg"), default="all")
            command.add_argument("--scope", choices=("release-signal", "all-benches"), default="release-signal")
            command.add_argument("--stem", type=Path, default=Path("target/bench-reports/performance"))
    commands.add_parser("report").add_argument("--stem", type=Path, default=Path("target/bench-reports/performance"))
    commands.add_parser("release-check").add_argument("current", help="Saved Criterion sample name")
    args = parser.parse_args(argv)
    try:
        root = find_project_root(Path.cwd()).resolve()
        if args.command == "release-check":
            print(f"Verified {check_release_sample(root, args.current)} complete saved Criterion benchmarks")
            return 0
        if args.command == "report":
            promote(root, root / args.stem)
            return 0
        modes: dict[str, PairMode] = {"local": "current-vs-latest", "release": "infer-release", "assets": "published-latest"}
        mode: PairMode = "explicit" if args.current or args.baseline else modes[args.command]
        pair = select_pair(root, mode, args.current or None, args.baseline or None, local=args.command == "local")
        if args.command == "assets":
            compare_assets(root, pair, root / "target/bench-reports/github-assets-performance.md")
        else:
            stem = root / args.stem
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


PLOT_CONFIG = "tooling/performance-plot.toml"
README_TABLE_STYLE = RunTableStyle(
    row_heading="Dimension",
    intervals="below",
    timing_legend=True,
    comparison_columns=("ratio",),
    ratio_headings=("Speedup vs nalgebra", "Speedup vs faer"),
    ratio_suffix="\N{MULTIPLICATION SIGN}",
)


def _repo_root() -> Path:
    """Resolve installed commands against the caller's checkout."""
    return find_project_root(Path.cwd()).resolve()


def presentation(root: Path, metric_name: str, statistic: Literal["mean", "median"], *, log_y: bool) -> tuple[RunLayout, CoordinatePlot]:
    """Map the scientific catalog onto the shared coordinate presentation plan."""
    config = tomllib.loads((root / PLOT_CONFIG).read_text(encoding="utf-8"))
    if (
        type(config["schema"]) is not int
        or config["schema"] != 1
        or tuple(config["dimensions"]) != scientific_policy(root).dimensions
        or config["unit"] != "ns"
    ):
        msg = "dimension publication requires the canonical D=2-64 selection"
        raise ValueError(msg)
    catalog = scientific_policy(root)
    metric = catalog.metrics[metric_name]
    benches = dict(zip(("la-stack current", "nalgebra", "faer"), metric.benchmarks, strict=True))
    selected = tuple(
        SeriesSelection(item["name"], item["label"], tuple((str(d), f"d{d}/{benches[item['name']]}") for d in catalog.dimensions)) for item in config["series"]
    )
    if tuple(item.name for item in selected) != tuple(benches):
        msg = "dimension publication requires current la-stack and both baseline peer views"
        raise ValueError(msg)
    layout = RunLayout(
        selected,
        tuple(RunRow(str(d), str(d), d) for d in catalog.dimensions),
        statistic,
        config["unit"],
        comparisons=(("nalgebra", "la-stack current"), ("faer", "la-stack current")),
    )
    settings = config["plot"]
    plot = CoordinatePlot(
        f"{metric.title}: {statistic} time vs dimension", settings["x-label"], settings["y-label"], settings["x-scale"], "log" if log_y else "linear"
    )
    return layout, plot


def _asset_links(document: bytes, *, metric: str, statistic: str, release: str) -> bytes:
    """Advance precisely the three existing la-stack artifact pointers."""
    stem = f"/docs/assets/bench/vs_linalg_{metric}_{statistic}"
    pattern = re.compile(
        rb"(?P<prefix>https://(?:github\.com/acgetchell/la-stack/(?:blob|raw|tree)/|raw\.githubusercontent\.com/acgetchell/la-stack/))"
        rb"v[0-9]+\.[0-9]+\.[0-9]+"
        rb"(?P<asset>" + re.escape(stem.encode()) + rb"(?:\.csv|\.provenance\.json|\.svg))"
    )
    found = sorted(match["asset"] for match in pattern.finditer(document))
    if found != sorted(f"{stem}{suffix}".encode() for suffix in (".csv", ".svg", ".provenance.json")):
        msg = "README requires exactly one tagged pointer to each selected artifact"
        raise ValueError(msg)
    return pattern.sub(lambda match: match["prefix"] + release.encode() + match["asset"], document)


def publish_plot(root: Path, metric: str, statistic: Literal["mean", "median"], *, log_y: bool, update_readme: bool) -> None:
    """Retain original phases, enforce scientific eligibility, and publish once."""
    # The report registry imports only the independent benchmark catalog.
    import performance_runs as policy  # noqa: PLC0415

    stem = root / "target/bench-reports/performance"
    evidence = policy.retained_run(root, stem)
    run = policy.validate_scientific_run(evidence, release=update_readme, root=root)
    verified_inputs = policy.verify_current_source(root, evidence) if update_readme else {}
    layout, plot = presentation(root, metric, statistic, log_y=log_y)
    phases = dict(evidence.sources)
    current_context, peer_context = dict(phases["current"].context), dict(phases["baseline"].context)
    versions = {"la-stack current": current_context["release"], "nalgebra": peer_context["dependency.nalgebra"], "faer": peer_context["dependency.faer"]}
    layout = replace(layout, series=tuple(replace(item, label=f"{item.label} ({versions[item.name]})") for item in layout.series))
    svg = render_coordinate_svg(evidence, layout, plot)
    csv = serialize_run_csv(evidence, layout)
    identity = run_identity(evidence)
    payload, manifest = policy.scratch_paths(stem)
    from_archive = not payload.exists()
    inputs = {
        PLOT_CONFIG: (root / PLOT_CONFIG).read_bytes(),
        policy.POLICY_CONFIG: (root / policy.POLICY_CONFIG).read_bytes(),
        policy.REPORT_CONFIG: (root / policy.REPORT_CONFIG).read_bytes(),
        **verified_inputs,
    }
    archive = root / policy.archive_directory(root)
    if not from_archive:
        inputs.update(zip((payload.relative_to(root).as_posix(), manifest.relative_to(root).as_posix()), serialize_evidence(evidence), strict=True))
    # README provenance must remain valid after scratch cleanup. Require the
    # selected measurement to have been promoted before referencing its archive.
    use_archive = from_archive or update_readme
    if use_archive:
        if not from_archive and run_identity(load_latest_run(root, policy.archive_directory(root))) != identity:
            msg = "README evidence differs from the retained latest run; run performance-doc first"
            raise ValueError(msg)
        payload, manifest = archive / "runs" / identity / "run.json", archive / "runs" / identity / "evidence.json"
        inputs.update({path.relative_to(root).as_posix(): path.read_bytes() for path in archive.rglob("*") if path.is_file()})
        inputs.update(zip((payload.relative_to(root).as_posix(), manifest.relative_to(root).as_posix()), serialize_evidence(evidence), strict=True))
    prefix = f"docs/assets/bench/vs_linalg_{metric}_{statistic}"
    provenance = {
        "schema": 4,
        "artifact": "la-stack dimension plot",
        "run_id": identity,
        "payload": payload.relative_to(root).as_posix(),
        "manifest": manifest.relative_to(root).as_posix(),
        "statistic": statistic,
        "unit": layout.unit,
        "series": {item.name: item.phase for item in run.series if item.name in {selected.name for selected in layout.series}},
    }
    outputs = {f"{prefix}.svg": svg, f"{prefix}.csv": csv, f"{prefix}.provenance.json": (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode()}
    git_checks = ()
    document: bytes | None = None
    if update_readme:
        context = dict(dict(evidence.sources)["current"].context)
        markers = MarkerPair(f"<!-- BENCH_TABLE:{metric}:{statistic}:new:BEGIN -->", f"<!-- BENCH_TABLE:{metric}:{statistic}:new:END -->")
        document = (root / "README.md").read_bytes()
        git_checks = (GitCheck(context["release"], tuple(outputs), allow_missing=True),)
        inputs["pyproject.toml"] = (root / "pyproject.toml").read_bytes()
        outputs["README.md"] = format_markdown(
            _asset_links(
                replace_section(document, markers, render_run_table(evidence, layout, style=README_TABLE_STYLE)),
                metric=metric,
                statistic=statistic,
                release=context["release"],
            ).decode("utf-8"),
            root / "README.md",
            root / "pyproject.toml",
        ).encode("utf-8")
    plan = plan_outputs(root, outputs, inputs=inputs)
    if update_readme and dict(plan.originals)["README.md"] != document:
        msg = "README changed during rendering"
        raise ValueError(msg)
    checks = ()
    if update_readme:
        config = policy.configuration(root)
        checks = tuple(
            GlobCheck(patterns, tuple(json.loads(current_context[field])))
            for field, patterns in (("source-inventory", config.sources), ("harness-inventory", config.harness))
        )
    # Preserve the shared planner's snapshots and add its public freshness
    # checks for the independently validated consumer selections.
    plan = replace(
        plan, git_checks=git_checks, glob_checks=checks, run_checks=(RunSelectionCheck(policy.archive_directory(root), identity),) if use_archive else ()
    )
    publish_publication(plan)


def dimension_main(argv: list[str] | None = None) -> int:
    """Plot complete evidence or publish its canonical README section."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metric", default="lu_solve", help="Metric declared in tooling/performance-policy.toml")
    parser.add_argument("--stat", choices=("mean", "median"), default="median")
    parser.add_argument("--sample", choices=("new",), default="new", help="Retained current series; raw Criterion samples are no longer publication inputs.")
    parser.add_argument("--log-y", action="store_true")
    parser.add_argument("--update-readme", action="store_true")
    args = parser.parse_args(argv)
    try:
        publish_plot(_repo_root(), args.metric, args.stat, log_y=args.log_y, update_readme=args.update_readme)
        return 0
    except (OSError, ValueError, TypeError, KeyError, MarkdownFormatError, ExecutableNotFoundError, subprocess.SubprocessError, ExceptionGroup) as error:
        print(f"dimension publication: {format_exception_diagnostics(error)}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(dimension_main(sys.argv[2:]) if sys.argv[1:2] == ["plot"] else main())
