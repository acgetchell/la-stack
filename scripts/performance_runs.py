"""la-stack scientific policy and series selection over shared complete runs."""

import json
import os
import re
import tomllib
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, cast

from research_repo_tools.common_measurement import CommonHarnessPlan, MeasurementPhase
from research_repo_tools.complete_runs import CompletePolicy, CompleteRun, RunSeries, validate_run_evidence
from research_repo_tools.criterion_inventory import DiscoveryCommand, DiscoveryPlan, discover_inventory, load_discovery
from research_repo_tools.evidence import Evidence, fingerprint_files, load_evidence
from research_repo_tools.measurement import MeasurementConfig, load_measurement
from research_repo_tools.release_discovery import normalize_tag
from research_repo_tools.release_pairs import ReleasePair
from research_repo_tools.run_reports import load_latest_run

if TYPE_CHECKING:
    from collections.abc import Mapping

REPORT_CONFIG = "tooling/performance-report.toml"
MEASUREMENT_CONFIG = "tooling/performance.toml"
POLICY_CONFIG = "tooling/performance-policy.toml"
NO_API_COMPATIBILITY = "none"
type ReleaseVersion = tuple[int, int, int]


@dataclass(frozen=True)
class PlotMetric:
    """One scientific operation and its current/baseline peer identifiers."""

    title: str
    benchmarks: tuple[str, str, str]


@dataclass(frozen=True)
class ScientificPolicy:
    """Parsed, immutable consumer declarations used by every performance path."""

    dimensions: tuple[int, ...]
    minimum_release: ReleaseVersion
    modern_api_release: ReleaseVersion
    legacy_api_cfg: str
    library_benchmarks: tuple[str, ...]
    unavailable_prefixes: tuple[str, ...]
    unavailable_benchmarks: tuple[str, ...]
    unavailable_benchmark_prefixes: tuple[str, ...]
    exact: Mapping[str, tuple[str, ...]]
    additional: Mapping[str, tuple[str, ...]]
    peers: Mapping[str, tuple[str, str]]
    metrics: Mapping[str, PlotMetric]


def _names(value: object, field: str, *, length: int | None = None) -> tuple[str, ...]:
    """Reject malformed or ambiguous Criterion names before constructing rows."""
    if not isinstance(value, list) or any(not isinstance(item, str) or re.fullmatch(r"[A-Za-z0-9_]+", item) is None for item in value):
        raise ValueError(f"{field} requires a list of Criterion names")
    names = cast("tuple[str, ...]", tuple(value))
    if len(set(names)) != len(names) or (length is not None and len(names) != length):
        raise ValueError(f"{field} requires distinct names" + (f" with exactly {length} entries" if length is not None else ""))
    return names


def _table(value: object, field: str) -> dict[str, object]:
    """Admit only named catalog tables, leaving their values to owning parsers."""
    if not isinstance(value, dict) or any(not isinstance(key, str) or re.fullmatch(r"[A-Za-z0-9_]+", key) is None for key in value):
        raise ValueError(f"{field} requires a named catalog table")
    return cast("dict[str, object]", value)


def _dimensions(value: object) -> tuple[int, ...]:
    """Keep dimensions positive, unique and in the declared plotting order."""
    if not isinstance(value, list) or not value or any(type(item) is not int or item <= 0 for item in value):
        msg = "dimensions requires positive integers"
        raise ValueError(msg)
    dimensions = cast("tuple[int, ...]", tuple(value))
    if tuple(sorted(set(dimensions))) != dimensions:
        msg = "dimensions requires increasing, distinct integers"
        raise ValueError(msg)
    return dimensions


def _metrics(value: object, peers: Mapping[str, tuple[str, str]]) -> Mapping[str, PlotMetric]:
    """Bind each plotted operation to its declared pair of baseline peers."""
    metrics = {}
    for name, record in _table(value, "metrics").items():
        fields = _table(record, f"metrics.{name}")
        title = fields.get("title")
        if not isinstance(title, str) or not title.strip():
            raise ValueError(f"metrics.{name}.title requires nonempty text")
        benchmarks = cast("tuple[str, str, str]", _names(fields.get("benchmarks"), f"metrics.{name}.benchmarks", length=3))
        if peers.get(benchmarks[0]) != benchmarks[1:]:
            raise ValueError(f"metrics.{name}.benchmarks differs from the declared peer pair")
        metrics[name] = PlotMetric(title, benchmarks)
    if not metrics:
        msg = "metrics requires at least one operation"
        raise ValueError(msg)
    return MappingProxyType(metrics)


def _release_threshold(value: object, field: str) -> ReleaseVersion:
    """Parse one canonical, nonnegative release threshold from TOML."""
    if not isinstance(value, list) or len(value) != 3 or any(type(item) is not int or item < 0 for item in value):
        raise ValueError(f"{field} requires three nonnegative integers")
    return cast("ReleaseVersion", tuple(value))


def find_project_root(start: Path | None = None) -> Path:
    """Resolve consumer commands from their Cargo checkout."""
    directory = (start or Path.cwd()).resolve()
    for candidate in (directory, *directory.parents):
        if (candidate / "Cargo.toml").is_file():
            return candidate
    msg = "Could not locate Cargo.toml to determine benchmark root"
    raise FileNotFoundError(msg)


def scientific_policy(root: Path | None = None) -> ScientificPolicy:
    """Parse the consumer-owned scientific catalog before discovery or measurement."""
    policy = tomllib.loads(((root or find_project_root()) / POLICY_CONFIG).read_text(encoding="utf-8"))
    if type(policy.get("schema")) is not int or policy["schema"] != 1:
        msg = "scientific policy requires schema 1"
        raise ValueError(msg)
    minimum = _release_threshold(policy.get("minimum-release"), "minimum-release")
    modern = _release_threshold(policy.get("modern-api-release"), "modern-api-release")
    if modern < minimum:
        msg = "modern API release cannot precede the minimum release"
        raise ValueError(msg)
    cfg = policy.get("legacy-api-cfg")
    if not isinstance(cfg, str) or re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", cfg) is None:
        msg = "legacy API cfg requires a Rust identifier"
        raise ValueError(msg)
    dimensions = _dimensions(policy.get("dimensions"))
    library = _names(policy.get("library-benchmarks"), "library-benchmarks")
    exact = {name: _names(names, f"exact.{name}") for name, names in _table(policy.get("exact"), "exact").items()}
    if not library or not exact or any(not names for names in exact.values()):
        msg = "scientific policy requires library and exact benchmark selections"
        raise ValueError(msg)
    additional = {dimension: _names(names, f"additional.{dimension}") for dimension, names in _table(policy.get("additional"), "additional").items()}
    if not set(additional) <= {str(value) for value in dimensions}:
        msg = "additional benchmark dimensions must occur in dimensions"
        raise ValueError(msg)
    peers = {name: cast("tuple[str, str]", _names(names, f"peers.{name}", length=2)) for name, names in _table(policy.get("peers"), "peers").items()}
    if not set(peers) <= {*library, *(name for names in additional.values() for name in names)}:
        msg = "peers must refer to declared library benchmarks"
        raise ValueError(msg)
    return ScientificPolicy(
        dimensions=dimensions,
        minimum_release=minimum,
        modern_api_release=modern,
        legacy_api_cfg=cfg,
        library_benchmarks=library,
        unavailable_prefixes=_names(policy.get("unavailable-prefixes"), "unavailable-prefixes"),
        unavailable_benchmarks=_names(policy.get("unavailable-benchmarks"), "unavailable-benchmarks"),
        unavailable_benchmark_prefixes=_names(policy.get("unavailable-benchmark-prefixes"), "unavailable-benchmark-prefixes"),
        exact=MappingProxyType(exact),
        additional=MappingProxyType(additional),
        peers=MappingProxyType(peers),
        metrics=_metrics(policy.get("metrics"), peers),
    )


def compatibility_adapter(pair: ReleasePair, root: Path | None = None) -> str:
    """Reject excluded releases and select the historical API cfg."""
    policy = scientific_policy(root)
    versions = {phase: tuple(map(int, normalize_tag(tag)[1:].split("."))) for phase, tag in (("baseline", pair.baseline), ("current", pair.current))}
    if any(version < policy.minimum_release for version in versions.values()):
        minimum_release = ".".join(map(str, policy.minimum_release))
        msg = f"unsupported comparison release; requires v{minimum_release} or newer"
        raise ValueError(msg)
    return policy.legacy_api_cfg if versions["baseline"] < policy.modern_api_release else NO_API_COMPATIBILITY


def required_report_ids(root: Path | None = None) -> set[str]:
    """Keep release discovery accountable to every report and plot selection."""
    return set(load_discovery(root or find_project_root(), "tooling/criterion-inventory.toml").required)


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


def selected_ids(suite: str, scope: str, expected: tuple[str, ...], root: Path | None = None) -> tuple[str, ...]:
    """Require the repository report registry, keeping diagnostic cases in the run."""
    if suite not in {"all", "exact", "vs_linalg"} or scope not in {"all-benches", "release-signal"}:
        msg = "unknown la-stack suite or report scope"
        raise ValueError(msg)
    if scope == "all-benches":
        return tuple(name for name in expected if not name.split("/")[0].removeprefix("d").isdigit() or "/la_stack_" in name)
    policy = scientific_policy(root)
    required = {
        *(f"{group}/{bench}" for group, benches in policy.exact.items() for bench in benches if suite != "vs_linalg"),
        *(
            f"d{dimension}/{bench}"
            for dimension in policy.dimensions
            for bench in (*policy.library_benchmarks, *policy.additional.get(str(dimension), ()))
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


def _historical_case_available(full_id: str, policy: ScientificPolicy) -> bool:
    """Mirror the groups and norm cases disabled by the historical harness cfg."""
    bench = full_id.partition("/")[2]
    return (
        not full_id.startswith(policy.unavailable_prefixes)
        and bench not in policy.unavailable_benchmarks
        and not bench.startswith(policy.unavailable_benchmark_prefixes)
    )


def measurement_plan(root: Path, pair: ReleasePair, suite: str, scope: str) -> CommonHarnessPlan:
    """Build the shared plan from native inventories and consumer API policy."""
    adapter = compatibility_adapter(pair, root)
    policy = scientific_policy(root)
    baseline_env = phase_environment(root, adapter)
    current_env = phase_environment(root)
    config = configuration(root)
    inventories = {}
    for phase in ("baseline", "current"):
        # Discover against the current API, then apply the exact groups excluded
        # by the historical cfg. Old method names cannot compile on current code.
        # The shared collector later checks this entire inventory on each source.
        declarations = load_discovery(root, "tooling/criterion-inventory.toml")
        discovery = DiscoveryPlan(
            tuple(
                DiscoveryCommand(item.name, (*config.command[:-1], f"{phase}-{item.name}", "list"), timeout=config.timeout)
                for item in declarations.suites
                if suite in {"all", item.name}
            )
        )
        ids = discover_inventory(root, discovery, env=os.environ | dict(current_env)).expected
        old_api = phase == "baseline" and adapter != NO_API_COMPATIBILITY
        inventories[phase] = tuple(sorted(name for name in ids if not old_api or _historical_case_available(name, policy)))
    return declared_plan(pair, suite, scope, inventories, baseline_env=baseline_env, current_env=current_env, config=config, root=root)


def declared_plan(  # noqa: PLR0913
    pair: ReleasePair,
    suite: str,
    scope: str,
    inventories: dict[str, tuple[str, ...]],
    *,
    baseline_env: tuple[tuple[str, str], ...],
    current_env: tuple[tuple[str, str], ...],
    config: MeasurementConfig | None = None,
    root: Path | None = None,
) -> CommonHarnessPlan:
    """Bind exact case inventories, reference phases, commands and completeness."""
    config = config or configuration(root)
    policy = scientific_policy(root)
    adapter = compatibility_adapter(pair, root)
    selected = selected_ids(suite, scope, inventories["current"], root)
    unavailable = adapter != NO_API_COMPATIBILITY
    baseline_rows = tuple(name for name in selected if not unavailable or _historical_case_available(name, policy))
    missing = set(baseline_rows) - set(inventories["baseline"])
    if missing:
        raise ValueError(f"baseline inventory omits selected rows: {sorted(missing)}")
    series = [
        RunSeries("la-stack baseline", "baseline", tuple((name, name) for name in baseline_rows)),
        RunSeries("la-stack current", "current", tuple((name, name) for name in selected)),
    ]
    peers = policy.peers
    for library, peer_index in (("nalgebra", 0), ("faer", 1)):
        rows = tuple(
            (name, f"{name.split('/')[0]}/{peers[name.split('/')[1]][peer_index]}") for name in selected if name.startswith("d") and name.split("/")[1] in peers
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
            ("baseline-api-compatibility", adapter),
        ),
    )
    return CommonHarnessPlan(config, phases["baseline"], phases["current"], tuple(series))


def validate_scientific_run(evidence: Evidence, *, release: bool = False, root: Path | None = None) -> CompleteRun:
    """Require la-stack's sampling, selected coverage and original reference phases."""
    run = validate_run_evidence(evidence)
    sources = dict(evidence.sources)
    context = dict(sources["current"].context)
    pair = ReleasePair(context["release"], dict(sources["baseline"].context)["release"])
    adapter = compatibility_adapter(pair, root)
    if release and pair.current == pair.baseline:
        msg = "cannot promote a same-version local comparison as release documentation"
        raise ValueError(msg)
    if pair.current == pair.baseline and sources["current"].source_sha256 == sources["baseline"].source_sha256:
        msg = "same-version local comparisons require distinct source identities"
        raise ValueError(msg)
    samples = dict(run.phases)
    if any(sample.policy != CompletePolicy(sample.policy.expected) for sample in samples.values()):
        msg = "la-stack requires 100 samples, mean and median, and 95% intervals in ns"
        raise ValueError(msg)
    expected = declared_plan(
        pair,
        context["suite"],
        context["scope"],
        {phase: sample.policy.expected for phase, sample in samples.items()},
        baseline_env=(),
        current_env=(),
        root=root,
    )
    if run.series != expected.series or context["baseline-api-compatibility"] != adapter:
        msg = "retained series or compatibility policy differs from la-stack declarations"
        raise ValueError(msg)
    if set(run.compatible) != set(expected.measurement.compatible):
        msg = "retained run omits required host, toolchain, or dependency compatibility"
        raise ValueError(msg)
    required = dict(expected.measurement.context)
    for phase, source in evidence.sources:
        actual = dict(source.context)
        if any(actual.get(key) != value for key, value in required.items()):
            msg = "retained phase differs from the declared consumer policy"
            raise ValueError(msg)
        for key, native_command in (("command", getattr(expected, phase).command), ("gate-command", getattr(expected, phase).gate)):
            recorded = tuple(json.loads(actual[key]))
            if recorded != native_command:
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
    validate_scientific_run(evidence, root=root)
    return evidence


def verify_current_source(root: Path, evidence: Evidence) -> dict[str, bytes]:
    """Bind README publication to the measured current library and captured harness."""
    source = dict(evidence.sources)["current"]
    context = dict(source.context)
    config = configuration(root)
    inputs: dict[str, bytes] = {}
    for key, expected in (("source-inventory", source.source_sha256), ("harness-inventory", source.harness_sha256)):
        names = cast("list[str]", json.loads(context[key]))
        patterns = config.sources if key == "source-inventory" else config.harness
        actual = {path.relative_to(root).as_posix() for pattern in patterns for path in root.glob(pattern)}
        if actual != set(names):
            raise ValueError(f"current {key} has stale file membership")
        observed = {name: (root / name).read_bytes() for name in names}
        if fingerprint_files(root, tuple(Path(name) for name in names)) != expected or any(
            (root / name).read_bytes() != data for name, data in observed.items()
        ):
            raise ValueError(f"current {key} is stale relative to retained measurements")
        inputs.update(observed)
    version = tomllib.loads(inputs["Cargo.toml"].decode("utf-8"))["package"]["version"]
    if f"v{version}" != context["release"]:
        msg = "retained release does not match the current package"
        raise ValueError(msg)
    return inputs
