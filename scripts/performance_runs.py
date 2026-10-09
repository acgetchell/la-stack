"""la-stack scientific policy and rendering over shared complete-run contracts."""

import json
import os
import tomllib
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

from research_repo_tools.common_measurement import CommonHarnessPlan, MeasurementPhase
from research_repo_tools.complete_runs import CompletePolicy, CompleteRun, RunSeries, run_identity, validate_run_evidence
from research_repo_tools.evidence import Evidence, fingerprint_files, load_evidence
from research_repo_tools.measurement import MeasurementConfig, load_measurement
from research_repo_tools.process import run_command
from research_repo_tools.release_pairs import ReleasePair
from research_repo_tools.run_reports import load_latest_run

import bench_compare as report
from benchmark_contract import find_project_root
from performance_artifacts import NO_API_COMPATIBILITY, resolve_shared_harness_compatibility
from release_baseline import parse_benchmark_list

if TYPE_CHECKING:
    from research_repo_tools.criterion import Estimate

REPORT_CONFIG = "tooling/performance-report.toml"
MEASUREMENT_CONFIG = "tooling/performance.toml"
UNAVAILABLE_PREFIXES = ("rational_input_", "canonical_conversion_", "det4_diagnostic_")


def configuration(root: Path | None = None) -> MeasurementConfig:
    """Load the shared measurement schema from the consumer checkout."""
    config = load_measurement(root or find_project_root(), MEASUREMENT_CONFIG)
    if not isinstance(config, MeasurementConfig):
        msg = "la-stack requires schema 1 with dynamically discovered Criterion inventories"
        raise ValueError(msg)  # noqa: TRY004 - wrong configuration schema, not a caller type error
    return config


def archive_directory(root: Path) -> str:
    """Read the history path; shared report APIs validate its path and contents."""
    archive = tomllib.loads((root / REPORT_CONFIG).read_text(encoding="utf-8")).get("archive")
    if not isinstance(archive, str) or not archive:
        msg = "performance report configuration requires an archive path"
        raise ValueError(msg)
    return archive


def selected_ids(suite: str, scope: str, expected: tuple[str, ...]) -> tuple[str, ...]:
    """Require the repository report registry, keeping diagnostic cases in the run."""
    if suite not in {"all", "exact", "vs_linalg"} or scope not in {"all-benches", "release-signal"}:
        msg = "unknown la-stack suite or report scope"
        raise ValueError(msg)
    if scope == "all-benches":
        return tuple(name for name in expected if not name.split("/")[0].removeprefix("d").isdigit() or "/la_stack_" in name)
    required = {
        *(f"{group}/{bench}" for group, benches in report.EXACT_GROUPS.items() for bench in benches if suite != "vs_linalg"),
        *(
            f"d{dimension}/{bench}"
            for dimension in report.VS_LINALG_CANONICAL_DIMS
            for bench in (*report.VS_LINALG_LA_STACK_BENCHES, *report.VS_LINALG_RELEASE_SIGNAL_BENCHES_BY_DIM.get(dimension, []))
            if suite != "exact"
        ),
    }
    missing = required - set(expected)
    if missing:
        raise ValueError(f"inventory omits selected report rows: {sorted(missing)}")
    return tuple(sorted(required))


def phase_environment(root: Path, adapter: str = NO_API_COMPATIBILITY) -> tuple[tuple[str, str], ...]:
    """Capture code-generation overrides and isolate each phase's output tree."""
    toolchain = tomllib.loads((root / "rust-toolchain.toml").read_text(encoding="utf-8"))["toolchain"]["channel"]
    env = {key: os.environ[key] for key in ("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "RUSTC", "RUSTC_WRAPPER", "CARGO_BUILD_TARGET") if key in os.environ}
    env.update(RUSTUP_TOOLCHAIN=os.environ.get("RUSTUP_TOOLCHAIN", toolchain), CARGO_TARGET_DIR="target", CRITERION_HOME="target/criterion")
    key, separator = ("CARGO_ENCODED_RUSTFLAGS", "\x1f") if "CARGO_ENCODED_RUSTFLAGS" in env else ("RUSTFLAGS", " ")
    flags = [env.get(key, ""), "--cap-lints=warn"]
    if adapter != NO_API_COMPATIBILITY:
        flags.append(f"--cfg={adapter}")
    env[key] = separator.join(flag for flag in flags if flag)
    return tuple(sorted(env.items()))


def _historical_case_available(full_id: str) -> bool:
    """Mirror the groups and norm cases disabled by the historical harness cfg."""
    bench = full_id.partition("/")[2]
    return (
        not full_id.startswith(UNAVAILABLE_PREFIXES)
        and bench != "la_stack_norm2"
        and not bench.startswith(("la_stack_norm2_scenario_", "iterative_f64_hypot_norm2_scenario_", "delaunay_scaled_norm_norm2_scenario_"))
    )


def measurement_plan(root: Path, pair: ReleasePair, suite: str, scope: str, listing: Path) -> CommonHarnessPlan:
    """Build the shared plan from native inventories and consumer API policy."""
    compatibility = resolve_shared_harness_compatibility(current=pair.current, baseline=pair.baseline, shared_harness_rational_inputs=True)
    baseline_env = phase_environment(root, compatibility.baseline_api_compatibility)
    current_env = phase_environment(root)
    config = configuration(root)
    inventories = {}
    for phase in ("baseline", "current"):
        # Discover against the current API, then apply the exact groups excluded
        # by the historical cfg. Old method names cannot compile on current code.
        # The shared collector later checks this entire inventory on each source.
        env = os.environ | dict(current_env) | {"CRITERION_HOME": str(listing / phase)}
        ids = []
        for selected in ("vs_linalg", "exact") if suite == "all" else (suite,):
            command = (*config.command[:-1], f"{phase}-{selected}", "list")
            result = run_command(command[0], command[1:], cwd=root, env=env, timeout=config.timeout)
            ids.extend(parse_benchmark_list(result.stdout))
        if len(ids) != len(set(ids)):
            msg = "duplicate semantic benchmark IDs"
            raise ValueError(msg)
        old_api = phase == "baseline" and compatibility.baseline_api_compatibility != NO_API_COMPATIBILITY
        inventories[phase] = tuple(sorted(name for name in ids if not old_api or _historical_case_available(name)))
    return declared_plan(pair, suite, scope, inventories, baseline_env=baseline_env, current_env=current_env, config=config)


def declared_plan(  # noqa: PLR0913
    pair: ReleasePair,
    suite: str,
    scope: str,
    inventories: dict[str, tuple[str, ...]],
    *,
    baseline_env: tuple[tuple[str, str], ...],
    current_env: tuple[tuple[str, str], ...],
    config: MeasurementConfig | None = None,
) -> CommonHarnessPlan:
    """Bind exact case inventories, reference phases, commands and completeness."""
    config = config or configuration()
    compatibility = resolve_shared_harness_compatibility(current=pair.current, baseline=pair.baseline, shared_harness_rational_inputs=True)
    selected = selected_ids(suite, scope, inventories["current"])
    unavailable = compatibility.baseline_api_compatibility != NO_API_COMPATIBILITY
    baseline_rows = tuple(name for name in selected if not unavailable or _historical_case_available(name))
    missing = set(baseline_rows) - set(inventories["baseline"])
    if missing:
        raise ValueError(f"baseline inventory omits selected rows: {sorted(missing)}")
    series = [
        RunSeries("la-stack baseline", "baseline", tuple((name, name) for name in baseline_rows)),
        RunSeries("la-stack current", "current", tuple((name, name) for name in selected)),
    ]
    for library, peer_index in (("nalgebra", 0), ("faer", 1)):
        rows = tuple(
            (name, f"{name.split('/')[0]}/{report.VS_LINALG_BASELINE_PEERS[name.split('/')[1]][peer_index]}")
            for name in selected
            if name.startswith("d") and name.split("/")[1] in report.VS_LINALG_BASELINE_PEERS
        )
        if rows:
            series.append(RunSeries(library, "baseline", rows))
    phases = {
        phase: MeasurementPhase(
            (*config.command[:-1], f"{phase}-{suite}"),
            (*config.command[:-1], "gate"),
            CompletePolicy(inventories[phase], sample_count=100, statistics=("mean", "median"), confidence_level=0.95),
            environment,
        )
        for phase, environment in (("baseline", baseline_env), ("current", current_env))
    }
    config = replace(
        config,
        command=phases["current"].command,
        context=(
            *config.context,
            ("suite", suite),
            ("scope", scope),
            ("baseline-api-compatibility", compatibility.baseline_api_compatibility),
        ),
    )
    return CommonHarnessPlan(config, phases["baseline"], phases["current"], tuple(series))


def _legacy_cargo_commands(suite: str, phase: str) -> str:
    """Validate the frozen pre-Just evidence contract without running its driver."""
    suites = {"vs_linalg": ("-p", "la-stack-comparison", "--features", "bench"), "exact": ("--features", "bench,exact")}
    return json.dumps(
        [
            ("cargo", "bench", "--locked", *features, "--bench", name, "--", *(("la_stack",) if phase == "current" and name == "vs_linalg" else ()), "--noplot")
            for name, features in suites.items()
            if suite in {"all", name}
        ]
    )


def validate_scientific_run(evidence: Evidence, *, release: bool = False) -> CompleteRun:  # noqa: C901
    """Require la-stack's sampling, selected coverage and original reference phases."""
    run = validate_run_evidence(evidence)
    sources = dict(evidence.sources)
    context = dict(sources["current"].context)
    pair = ReleasePair(context["release"], dict(sources["baseline"].context)["release"])
    compatibility = resolve_shared_harness_compatibility(current=pair.current, baseline=pair.baseline, shared_harness_rational_inputs=True)
    if release and pair.current == pair.baseline:
        msg = "cannot promote a same-version local comparison as release documentation"
        raise ValueError(msg)
    if pair.current == pair.baseline and sources["current"].source_sha256 == sources["baseline"].source_sha256:
        msg = "same-version local comparisons require distinct source identities"
        raise ValueError(msg)
    samples = dict(run.phases)
    for sample in samples.values():
        if sample.policy != CompletePolicy(sample.policy.expected):
            msg = "la-stack requires 100 samples, mean and median, and 95% intervals in ns"
            raise ValueError(msg)
    expected = declared_plan(
        pair, context["suite"], context["scope"], {phase: sample.policy.expected for phase, sample in samples.items()}, baseline_env=(), current_env=()
    )
    if run.series != expected.series or context["baseline-api-compatibility"] != compatibility.baseline_api_compatibility:
        msg = "retained series or compatibility policy differs from la-stack declarations"
        raise ValueError(msg)
    legacy = "phase-driver" not in context
    compatible = set(expected.measurement.compatible) - ({"context.tool.just"} if legacy else set())
    if set(run.compatible) != compatible:
        msg = "retained run omits required host, toolchain, or dependency compatibility"
        raise ValueError(msg)
    required = dict(expected.measurement.context)
    if legacy:
        del required["phase-driver"]
        required.update((f"{phase}-cargo-commands", _legacy_cargo_commands(context["suite"], phase)) for phase in ("baseline", "current"))
    for phase, source in evidence.sources:
        actual = dict(source.context)
        if actual.get("phase-driver") != context.get("phase-driver") or any(actual.get(key) != value for key, value in required.items()):
            msg = "retained phase differs from the declared consumer policy"
            raise ValueError(msg)
        for key, command in (("command", getattr(expected, phase).command), ("gate-command", getattr(expected, phase).gate)):
            recorded = tuple(json.loads(actual[key]))
            if legacy:
                # Interpreter paths varied across hosts. Old evidence remains
                # readable without retaining or executing the Python driver.
                command = ("scripts/performance_phase.py", "measure" if key == "command" else "gate", context["suite"], phase)
                recorded = recorded[1:]
            if recorded != command:
                raise ValueError(f"retained {phase} {key} differs from the declared Cargo phase")
    return run


def scratch_paths(stem: Path) -> tuple[Path, Path]:
    """Use new paths and framing, independent of historical CSV/JSON evidence."""
    return stem.with_suffix(".run.json"), stem.with_suffix(".evidence.json")


def retained_run(root: Path, stem: Path) -> Evidence:
    """Resolve validated shared history only when every scratch input is absent."""
    paths = scratch_paths(stem)
    legacy = (stem.with_suffix(".csv"), stem.with_suffix(".provenance.json"), stem.with_suffix(".full.csv"), stem.with_suffix(".full.provenance.json"))
    # Partial new or legacy scratch must never silently select older history.
    present = any(path.exists() or path.is_symlink() for path in (*paths, *legacy))
    evidence = load_evidence(*paths) if present else load_latest_run(root, archive_directory(root))
    validate_scientific_run(evidence)
    return evidence


def verify_current_source(root: Path, evidence: Evidence) -> None:
    """Bind README publication to the measured current library and captured harness."""
    source = dict(evidence.sources)["current"]
    context = dict(source.context)
    config = configuration(root)
    for key, expected in (("source-inventory", source.source_sha256), ("harness-inventory", source.harness_sha256)):
        names = json.loads(context[key])
        patterns = config.sources if key == "source-inventory" else config.harness
        actual = {path.relative_to(root).as_posix() for pattern in patterns for path in root.glob(pattern)}
        if actual != set(names):
            raise ValueError(f"current {key} has stale file membership")
        if fingerprint_files(root, tuple(Path(name) for name in names)) != expected:
            raise ValueError(f"current {key} is stale relative to retained measurements")
    version = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    if f"v{version}" != context["release"]:
        msg = "retained release does not match the current package"
        raise ValueError(msg)


def estimates_by_series(run: CompleteRun, statistic: str = "median") -> dict[str, dict[str, Estimate]]:
    """Select each named view from its actual measurement phase."""
    if statistic not in {"mean", "median"}:
        raise ValueError(f"unsupported statistic: {statistic}")
    phases = {phase: {case.full_id: case for case in sample.cases} for phase, sample in run.phases}
    return {series.name: {label: phases[series.phase][full_id].estimate(statistic) for label, full_id in series.rows} for series in run.series}


def render_scientific_report(evidence: Evidence) -> bytes:
    """Retain la-stack's comparison tables and cautious statistical interpretation."""
    run = validate_scientific_run(evidence)
    sources = dict(evidence.sources)
    current, baseline = (dict(sources[phase].context)["release"] for phase in ("current", "baseline"))
    values = estimates_by_series(run)
    comparisons = []
    current_only = []
    for name, estimate in values["la-stack current"].items():
        group, bench = name.split("/", 1)
        prior = values["la-stack baseline"].get(name)
        if prior is None:
            current_only.append(f"- `{name}`: current {report.format_time(estimate.point)}; baseline API unavailable.")
            continue

        comparisons.append(
            report.Comparison(
                suite="exact" if group.startswith(("exact_", *UNAVAILABLE_PREFIXES)) else "vs_linalg",
                group=group,
                bench=bench,
                baseline_bench=bench,
                baseline=prior,
                current=estimate,
                assessment=report.assess_change(prior, estimate),
                baseline_nalgebra=values.get("nalgebra", {}).get(name),
                baseline_faer=values.get("faer", {}).get(name),
            )
        )
    lines = [
        "# Benchmark Performance",
        "",
        f"**la-stack** {current} · `{sources['current'].revision}`",
        "",
        f"Comparison against baseline **{baseline}**:",
        "",
        "**Statistic**: median; complete mean and median measurements are retained.",
        "",
        "Negative point-estimate change means a smaller current estimate; a baseline/current ratio above 1 has the same meaning.",
        "Marginal Criterion interval separation is not a paired confidence interval or a statistical-significance claim.",
        "",
        "nalgebra/faer reference measurements originate in the baseline phase under the captured current harness; they were not rerun in the current phase.",
        "",
        f"**Complete run**: `{run_identity(evidence)}`",
        "",
    ]
    for phase, source in evidence.sources:
        context = dict(source.context)
        lines.extend(
            [
                f"- {phase}:",
                f"  - Revision: `{source.revision}`.",
                f"  - Source: `{source.source_sha256}`.",
                f"  - Harness: `{source.harness_sha256}`.",
                f"  - Gate status: {context['gate-status']}.",
                f"  - Gate: `{context['gate-command']}`.",
                f"  - Command: `{context['command']}`.",
            ]
        )
    lines.extend(
        [
            "",
            report.comparison_tables(comparisons, baseline),
            "",
            *current_only,
            "",
            "Regenerate from retained evidence with `just performance-doc`; update README with `just performance-readme`.",
            "",
        ]
    )
    return "\n".join(lines).encode()
