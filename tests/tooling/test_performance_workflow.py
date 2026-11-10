"""Installed workflow acceptance with independent la-stack scientific declarations."""

import csv
import fnmatch
import hashlib
import io
import json
import os
import shlex
import shutil
import subprocess
import sys
import tomllib
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import pytest
import yaml
from research_repo_tools.candidate_reuse import parse_candidate_policy
from research_repo_tools.cli import main
from research_repo_tools.common_measurement import CompleteMeasurementPlan, DiscoveryPolicy, validate_run_policy
from research_repo_tools.evidence import deterministic_json, load_evidence, publish_evidence
from research_repo_tools.just_inspect import dry_run, inspect_justfile
from research_repo_tools.measurement import load_measurement, measure_pair
from research_repo_tools.performance_progress import ProgressReporter
from research_repo_tools.performance_release import Preparation, serialize_preparation
from research_repo_tools.performance_reports import load_report_plan
from research_repo_tools.process import run_command, run_command_live
from research_repo_tools.publication import publish_publication
from research_repo_tools.release_assets import GitHubRelease, ReleaseTarget, serialize_release_target
from research_repo_tools.release_pairs import ReleasePair
from research_repo_tools.scheduled_runs import ScheduledRun, select_trial
from research_repo_tools.suite_measurement import SuiteMeasurementPlan
from research_repo_tools.suite_runs import comparison_incompatibilities, parse_suite_run
from research_repo_tools.worktrees import capture_snapshot

from tooling.release_candidate import ReviewedAcceptance, archive, capture, freeze, receive, remember_owner, selected_owner, verify

if TYPE_CHECKING:
    from collections.abc import Sequence

    from research_repo_tools.evidence import Evidence

ROOT = Path(__file__).resolve().parents[2]
PAIR = ReleasePair("v0.5.2", "v0.5.1")
DRIVER = """import json, os, sys
from pathlib import Path
phase, suite, action = sys.argv[1:]
if phase == "source":
    phase = Path("src/lib.rs").read_text(encoding="utf-8").strip().removeprefix("// ").split()[0]
assert suite in IDS[phase], "unsupported API was compiled/invoked"
if action == "gate":
    Path("target").mkdir(exist_ok=True)
    with Path("target/gates.txt").open("a", encoding="utf-8", newline="\\n") as stream:
        stream.write("gate\\n")
    raise SystemExit(0)
ids = IDS[phase][suite]
if action == "list":
    for name in ids:
        print(name + ": benchmark")
    raise SystemExit(0)
assert Path("target/gates.txt").exists()
for name in ids:
    path = Path(os.environ["CRITERION_HOME"]) / name / "new"
    path.mkdir(parents=True)
    median = TIMINGS[phase][name]["median"]
    mean = TIMINGS[phase][name]["mean"]
    def estimate(value):
        return {"point_estimate": value, "confidence_interval": {"lower_bound": value * .9, "upper_bound": value * 1.1, "confidence_level": .95}}
    data = {
        "benchmark.json": {"full_id": name},
        "estimates.json": {"mean": estimate(mean), "median": estimate(median)},
        "sample.json": {"iters": [1] * 100, "times": [median] * 100},
        "tukey.json": [median * .5, median * .75, median * 1.25, median * 1.5],
    }
    for filename, value in data.items():
        (path / filename).write_text(json.dumps(value), encoding="utf-8", newline="\\n")
"""


def declaration(root: Path, view: str = "full") -> CompleteMeasurementPlan | SuiteMeasurementPlan:
    plan = load_measurement(root, f"tooling/measurements/{view}.toml")
    assert isinstance(plan, (CompleteMeasurementPlan, SuiteMeasurementPlan))
    return plan


def scientific_ids(view: str = "full") -> tuple[set[str], set[str]]:
    """Derive expectations from the independent domain catalog, not a measured run."""
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    if view == "full":
        return set(catalog["release"]["default"] + catalog["release"]["exact"]), {
            f"d{dimension}/{peer}_lu_solve" for dimension in catalog["dimensions"] for peer in ("nalgebra", "faer")
        }
    library = {
        *(
            f"d{dimension}/{name}"
            for dimension in catalog["dimensions"]
            for name in [*catalog["library-benchmarks"], *catalog["additional"].get(str(dimension), [])]
        ),
    }
    peers = {f"{name.split('/')[0]}/{peer}" for name in library for peer in catalog["peers"].get(name.split("/")[1], [])}
    return library, peers


def baseline_ids(view: str = "full", release: str = PAIR.baseline) -> set[str]:
    """Derive historical expectations from the independent scientific API catalog."""
    library, _peers = scientific_ids(view)
    if view != "full":
        return library
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    version = tuple(map(int, release.removeprefix("v").split(".")))
    return {
        name
        for name in library
        if not any(
            fnmatch.fnmatchcase(name, pattern) and version < tuple(map(int, since.removeprefix("v").split(".")))
            for pattern, since in catalog["availability"].items()
        )
    }


def suite_ids() -> dict[str, set[str]]:
    """Partition the scientific catalog independently of the acquisition declaration."""
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    library, peers = scientific_ids()
    families = {
        "column-basis": {name for name in library if name.startswith("column_basis/")},
        "interval-solve": {name for name in library if name.startswith(("interval_division/", "interval_solve_d"))},
        "polynomial": {name for name in library if name.startswith("polynomial/")},
        "rational-dot": {name for name in library if name.startswith("rational_dot/")},
        "lu-solve": {name for name in library if name.endswith("/la_stack_lu_solve")} | peers,
    }
    separate = set().union(*families.values())
    return {"core-default": set(catalog["release"]["default"]) - separate, "core-exact": set(catalog["release"]["exact"]) - separate, **families}


def synthetic_timings() -> dict[str, dict[str, dict[str, float]]]:
    """Distinguish every case, phase, and statistic so selection errors stay visible."""
    library, peers = scientific_ids()
    diagnostic_library, diagnostic_peers = scientific_ids("non-exact")
    return {
        phase: {
            name: {"median": float(offset + index * 4), "mean": float(offset + index * 4 + 1)}
            for index, name in enumerate(sorted(library | peers | diagnostic_library | diagnostic_peers), start=1)
        }
        for phase, offset in (("current", 0), ("baseline", 10_000))
    }


@pytest.mark.parametrize("release", ["v0.5.0", "v0.5.1", "v0.5.2", "v0.5.3"])
def test_declared_scientific_sampling_coverage_and_original_peer_phases(release: str) -> None:
    plan = declaration(ROOT)
    assert isinstance(plan, SuiteMeasurementPlan)
    expected = suite_ids()
    library, peers = scientific_ids()
    assert len(library) == 1403
    assert len(peers) == 16
    assert set().union(*expected.values()) == library | peers
    assert sum(map(len, expected.values())) == len(library | peers)
    available: set[str] = set()
    for suite in plan.suites:
        assert suite.nominal_seconds_per_case == 8
        assert suite.measurement.baseline.command == suite.measurement.current.command
        assert suite.measurement.current.gate == (*suite.measurement.baseline.gate[:-1], suite.measurement.baseline.gate[-1] + "-current")
        assert suite.measurement.baseline.sources == ("src/**/*.rs", "Cargo.toml", "README.md")
        assert "README.md" not in suite.measurement.sources
        assert dict(suite.measurement.measurement.context)["acquisition-owner"] == "local"
        for phase in (suite.measurement.baseline, suite.measurement.current):
            assert isinstance(phase.policy, DiscoveryPolicy)
            assert set(phase.policy.discovery.required) == expected[suite.name]
            assert phase.policy.sample_count == 100
            assert phase.policy.statistics == ("mean", "median")
            assert phase.policy.confidence_level == 0.95
            assert phase.policy.unit == "ns"
            environment = dict(phase.environment)
            assert environment["RUSTFLAGS"] == environment["CARGO_ENCODED_RUSTFLAGS"] == "--cap-lints=warn"
            assert environment["CARGO_TARGET_DIR"] == "target/build"
            assert environment["CARGO_BUILD_BUILD_DIR"] == "target/build/intermediates"
            assert environment["RUSTUP_TOOLCHAIN"] == tomllib.loads((ROOT / "rust-toolchain.toml").read_text(encoding="utf-8"))["toolchain"]["channel"]
        if suite.decision(release, "{}")["state"] == "available":
            available.update(expected[suite.name])
        assert {"context.os", "context.architecture", "context.cpu", "harness_sha256"} <= set(suite.measurement.measurement.compatible)
    assert available == baseline_ids(release=release) | peers
    lu = next(suite for suite in plan.suites if suite.name == "lu-solve")
    assert {item.name: item.phase for item in lu.measurement.series} == {
        "la-stack baseline": "baseline",
        "la-stack current": "current",
        "nalgebra": "current",
        "faer": "current",
    }
    assert {
        "benches/**/*.rs",
        "tests/exact_bench_config.rs",
        "tooling/performance.just",
        "tooling/measurements/*.toml",
        "pyproject.toml",
        "uv.lock",
        "dependency-sources.toml",
        ".python-version",
        "tests/column_basis.rs",
        "tests/rational_dot.rs",
        "tests/interval_solve.rs",
        "tests/common/interval_solve_oracle.rs",
        "tests/common/proptest_config.rs",
        "tests/polynomial_bench_inputs.rs",
    } <= set(plan.harness)


def test_non_exact_diagnostic_keeps_original_peer_phases() -> None:
    plan = declaration(ROOT, "non-exact")
    assert isinstance(plan, CompleteMeasurementPlan)
    library, peers = scientific_ids("non-exact")
    for phase, required in ((plan.baseline, library | peers), (plan.current, library)):
        assert isinstance(phase.policy, DiscoveryPolicy)
        assert set(phase.policy.discovery.required) == required
    assert {item.name: item.phase for item in plan.series} == {
        "la-stack baseline": "baseline",
        "la-stack current": "current",
        "nalgebra": "baseline",
        "faer": "baseline",
    }


def test_routine_ci_validates_benchmarks_without_acquiring_timings() -> None:
    """Preserve the native CI checks and reserve sampling for explicit requests."""
    workflows = ROOT / ".github/workflows"
    ci = yaml.safe_load((workflows / "ci.yml").read_bytes())
    ci_events = ci.get("on", ci.get(True))
    assert ci_events["push"]["branches"] == ci_events["pull_request"]["branches"] == ["main"]
    build = ci["jobs"]["build"]
    platforms = build["strategy"]["matrix"]["os"]
    assert set(platforms) == {"ubuntu-latest", "macos-latest", "windows-latest"}
    assert build["needs"] == ["prepare", "native"]
    native = ci["jobs"]["native"]
    assert "rust" in native["strategy"]["matrix"]["stage"]
    assert {
        "name": "Run native validation",
        "shell": "bash",
        "env": {"CI_STAGE": "${{ matrix.stage }}"},
        "run": 'just "ci-$CI_STAGE"',
    } in native["steps"]
    ruleset = json.loads((ROOT / ".github/settings/main-ruleset.json").read_bytes())
    required = next(rule["parameters"]["required_status_checks"] for rule in ruleset["rules"] if rule["type"] == "required_status_checks")
    contexts = {check["context"] for check in required}
    assert {f"build ({platform})" for platform in platforms} <= contexts

    commands = dry_run(ROOT, "ci").stderr.splitlines()
    assert any("cargo nextest run --workspace --release --profile ci --all-features --lib --tests" in line for line in commands)
    bench_commands = [line.split(" cargo bench ", 1)[1] for line in commands if " cargo bench " in line]
    compile_commands = [line for line in bench_commands if "--no-run" in line.split()]
    assert any("--workspace --no-run --features bench" in line for line in compile_commands)
    assert any("--no-run --features bench,exact --bench exact" in line for line in compile_commands)
    assert [line for line in bench_commands if line not in compile_commands] == ['--locked --features bench,exact --bench exact -- "$@"']
    recipes = inspect_justfile(ROOT).recipes
    assert "_ci-rust" in {dep["recipe"] for dep in recipes["ci-rust"]["dependencies"]}
    assert {"test-bench-exact", "bench-compile", "test-rust-ci"} <= {dep["recipe"] for dep in recipes["_ci-rust"]["dependencies"]}
    smoke = recipes["test-bench-exact"]["dependencies"]
    assert [(dep["recipe"], dep["arguments"]) for dep in smoke] == [("bench-exact", ["--test"])]

    diagnostic = yaml.safe_load((workflows / "benchmarks.yml").read_bytes())
    assert set(diagnostic.get("on", diagnostic.get(True))) == {"workflow_dispatch"}
    assert diagnostic["permissions"] == {"contents": "read"}
    assert {"bench", "Benchmarks / bench"}.isdisjoint(contexts)
    steps = diagnostic["jobs"]["bench"]["steps"]
    assert [step["run"] for step in steps if "run" in step] == ["just test-bench-inputs", "just bench-exact --noplot"]
    assert steps[-1]["with"]["path"] == "target/criterion/"
    assert steps[-1]["with"]["name"].startswith("exploratory-exact-")
    assert not any(step.get("uses", "").startswith("actions/download-artifact@") for step in steps)


@pytest.mark.parametrize(
    ("pattern", "count"),
    [("column_basis/*", 55), ("interval_division/*", 5), ("interval_solve_d*/*", 30), ("polynomial/*", 5), ("rational_dot/*", 54)],
)
def test_new_api_coverage_is_complete_without_inventing_baseline_timings(pattern: str, count: int) -> None:
    library, _peers = scientific_ids()
    current = {name for name in library if fnmatch.fnmatchcase(name, pattern)}
    assert len(current) == count
    assert current.isdisjoint(baseline_ids(release="v0.5.1"))
    assert current <= baseline_ids(release="v0.5.2")


def test_named_views_keep_every_metric_statistic_axis_and_explicit_input_choice() -> None:
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    for source in ("scratch", "retained", "non-exact"):
        files = sorted((ROOT / "tooling/plots" / source).glob("*.toml"))
        assert len(files) == (len(catalog["metrics"]) if source == "non-exact" else 1) * 4
        for path in files:
            metric, statistic, axis = path.stem.split("-")
            view = tomllib.loads(path.read_text(encoding="utf-8"))
            assert view["statistic"] == statistic
            assert statistic in {"mean", "median"}
            assert view["unit"] == "ns"
            assert view["plot"]["y-scale"] == axis
            assert axis in {"linear", "log"}
            assert [row["coordinate"] for row in view["rows"]] == catalog["dimensions"]
            assert [series["name"] for series in view["series"]] == [
                ("" if source == "non-exact" else "lu-solve/") + name for name in ("la-stack current", "nalgebra", "faer")
            ]
            for series, benchmark in zip(view["series"], catalog["metrics"][metric]["benchmarks"], strict=True):
                assert series["rows"] == {str(d): f"d{d}/{benchmark}" for d in catalog["dimensions"]}
            measurement = "non-exact" if source == "non-exact" else "full"
            assert view["measurement"] == f"tooling/measurements/{measurement}.toml"
            assert ("archive" in view) == (source == "retained")
            if source != "retained":
                assert view["payload"] == f"target/bench-reports/{measurement}/run.json"
                assert view["manifest"] == f"target/bench-reports/{measurement}/evidence.json"
            else:
                assert "payload" not in view
                assert "manifest" not in view


def test_canonical_readme_selects_all_lu_solve_medians_and_peer_over_current_ratios() -> None:
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    view = tomllib.loads((ROOT / "tooling/performance-readme.toml").read_text(encoding="utf-8"))
    assert view["measurement"] == "tooling/measurements/full.toml"
    assert view["trial"] == 1
    assert view["statistic"] == "median"
    assert view["unit"] == "ns"
    assert view["rows"] == [{"key": str(d), "label": str(d), "coordinate": d} for d in catalog["dimensions"]]
    assert [series["name"] for series in view["series"]] == ["lu-solve/la-stack current", "lu-solve/nalgebra", "lu-solve/faer"]
    for series, benchmark in zip(view["series"], catalog["metrics"]["lu_solve"]["benchmarks"], strict=True):
        assert series["rows"] == {str(d): f"d{d}/{benchmark}" for d in catalog["dimensions"]}
    assert view["comparisons"] == [
        {"numerator": "lu-solve/nalgebra", "denominator": "lu-solve/la-stack current"},
        {"numerator": "lu-solve/faer", "denominator": "lu-solve/la-stack current"},
    ]


@pytest.mark.parametrize(
    ("failure", "diagnostic"),
    [
        ("none", ""),
        ("readme", "current compilation reads README.md"),
        ("case-alias", "current compilation reads README.md"),
        ("missing-bench", "missing benchmark-discovery or correctness-gate"),
        ("missing-gate", "missing benchmark-discovery or correctness-gate"),
        ("empty", "missing library input"),
    ],
)
def test_current_input_gate_requires_both_native_profiles(tmp_path: Path, failure: str, diagnostic: str) -> None:
    """Exercise the actual native-Just input guard without compiling or measuring."""
    for name in ("pyproject.toml", "uv.lock", ".python-version"):
        shutil.copyfile(ROOT / name, tmp_path / name)
    for profile in ("debug", "release"):
        if (failure, profile) in {("missing-bench", "release"), ("missing-gate", "debug")}:
            continue
        path = tmp_path / "target/build" / profile / "deps/la_stack-abcdef.d"
        path.parent.mkdir(parents=True)
        data = b"target: src/lib.rs\n"
        if profile == "release":
            if failure == "readme":
                data += b"src/../README.md:\n"
            elif failure == "case-alias":
                data += b"C:\\fixture\\src\\..\\ReadMe.MD:\r\n"
            elif failure == "empty":
                data = b""
        path.write_bytes(data)
    result = run_command(
        "just",
        ["--no-deps", "--justfile", str(ROOT / "tooling/performance.just"), "--working-directory", str(tmp_path), "_readme-input-check"],
        cwd=tmp_path,
        env={**os.environ, "UV_PROJECT_ENVIRONMENT": sys.prefix},
        check=False,
    )
    assert result.returncode == (0 if failure == "none" else 1), result.stderr
    if diagnostic:
        assert diagnostic in result.stderr


def test_checkpoint_resume_preserves_completed_scientific_pairs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Recover the real LU/column-basis declarations without repeating LU samples."""
    pair = ReleasePair("v0.5.3", "v0.5.2")
    _baseline, root = prepared(tmp_path, pair)
    full = declaration(root)
    assert isinstance(full, SuiteMeasurementPlan)
    selected = replace(full, suites=tuple(suite for suite in full.suites if suite.name in {"lu-solve", "column-basis"}))
    checkpoint = root / "target/checkpoint"
    checkpoint.parent.mkdir()
    calls: list[tuple[str, str]] = []
    interrupted = False

    def native(command: str, args: Sequence[str] = (), **kwargs: Any) -> Any:
        nonlocal interrupted
        if args and args[-1] in {"lu-solve", "column-basis"}:
            phase = (cast("Path", kwargs["cwd"]) / "src/lib.rs").read_text(encoding="utf-8").strip().removeprefix("// ")
            if (args[-1], phase) == ("column-basis", "current") and not interrupted:
                interrupted = True
                raise KeyboardInterrupt
            calls.append((args[-1], phase))
        return run_command_live(command, args, **kwargs)

    monkeypatch.setattr("research_repo_tools.common_measurement.run_command_live", native)
    with pytest.raises(KeyboardInterrupt):
        measure_pair(root, selected, pair, working_tree=True, checkpoint=checkpoint, progress=ProgressReporter(mode="none"))
    original = (checkpoint / "events.jsonl").read_bytes()
    assert calls == [("lu-solve", "baseline"), ("lu-solve", "current"), ("column-basis", "baseline")]
    evidence = measure_pair(root, selected, pair, working_tree=True, checkpoint=checkpoint, resume=True, progress=ProgressReporter(mode="none"))
    validate_run_policy(evidence, selected, release=True)
    assert (checkpoint / "events.jsonl").read_bytes().startswith(original)
    assert calls == [
        ("lu-solve", "baseline"),
        ("lu-solve", "current"),
        ("column-basis", "baseline"),
        ("column-basis", "baseline"),
        ("column-basis", "current"),
    ]
    trial = select_trial(evidence, 1)
    for dimension in (2, 3, 4, 5, 8, 16, 32, 64):
        name = f"d{dimension}/la_stack_lu_solve"
        assert comparison_incompatibilities(trial, ("baseline", name), ("current", name)) == ()
    # Input drift rejects before another native measurement and preserves the journal.
    completed = (checkpoint / "events.jsonl").read_bytes()
    before_calls = list(calls)
    (root / "src/lib.rs").write_text("// changed current input\n", encoding="utf-8", newline="\n")
    with pytest.raises(ValueError, match="identity changed"):
        measure_pair(root, selected, pair, working_tree=True, checkpoint=checkpoint, resume=True, progress=ProgressReporter(mode="none"))
    assert calls == before_calls
    assert (checkpoint / "events.jsonl").read_bytes() == completed


def recipe(root: Path, name: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    """Invoke the actual recipe with the already-installed locked consumer environment."""
    return run_command(
        "just",
        ["--no-deps", "--justfile", str(ROOT / "justfile"), "--working-directory", str(root), name, *arguments],
        env={**os.environ, "UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
        # Full-inventory publication exceeded the shared 300-second command
        # budget on macOS and Windows CI; allow 3x for these workflow tests.
        timeout=900,
        check=False,
    )


def test_native_recipes_preserve_scientific_gates_and_peer_reuse() -> None:
    justfile = Path("tooling/performance.just")
    plan = declaration(ROOT)
    assert isinstance(plan, SuiteMeasurementPlan)
    measured_targets: set[str] = set()
    for suite in plan.suites:
        command = suite.measurement.current.command[-1]
        measured = dry_run(ROOT, command, justfile=justfile).stderr
        discovered = dry_run(ROOT, command, ["list"], justfile=justfile).stderr
        assert measured.count("cargo bench --locked") == discovered.count("args+=(--list)")
        assert "--target" in measured
        for manifest in (ROOT / "Cargo.toml", ROOT / "benches/comparison/Cargo.toml"):
            for target in tomllib.loads(manifest.read_text(encoding="utf-8"))["bench"]:
                if f"--bench '{target['name']}'" in measured:
                    measured_targets.add(target["name"])
        gate = dry_run(ROOT, suite.measurement.current.gate[-1], justfile=justfile).stderr
        assert "cargo test --locked" in gate
        if suite.since:
            assert suite.name.replace("-", "_") in gate or "polynomial_bench_inputs" in gate
    targets = {
        bench["name"]
        for manifest in (ROOT / "Cargo.toml", ROOT / "benches/comparison/Cargo.toml")
        for bench in tomllib.loads(manifest.read_text(encoding="utf-8"))["bench"]
    }
    assert measured_targets == targets
    core = dry_run(ROOT, "core-default", justfile=justfile).stderr
    assert "|lu|lu_pivoting|lu_ill_conditioned|" in core
    assert "la_stack_lu_solve" not in core
    lu = dry_run(ROOT, "lu-solve", justfile=justfile).stderr
    assert "(la_stack|nalgebra|faer)_lu_solve$" in lu
    assert lu.count("cargo bench --locked") == 1
    for name in ("performance-local", "performance-local-non-exact"):
        command = dry_run(ROOT, name).stderr
        assert "research-repo-tools performance measure" in command
        assert "--mode current-vs-latest" in command
    command = dry_run(ROOT, "performance-release").stderr
    assert command.count("tooling/release_candidate.py acquire") == 1
    assert command.count("--release-intent") == 1


def prepared(tmp_path: Path, pair: ReleasePair = PAIR) -> tuple[Path, Path]:
    """Replace only native timing/gates with a bounded fixture, retaining real policy."""
    _library, peers = scientific_ids()
    diagnostic_library, diagnostic_peers = scientific_ids("non-exact")
    ids = {
        side: {
            **{name: sorted(cases) for name, cases in suite_ids().items() if side == "current" or cases <= baseline_ids(release=pair.baseline) | peers},
            "all": [],
            "vs_linalg": sorted(diagnostic_library | (diagnostic_peers if side == "baseline" else set())),
        }
        for side in ("baseline", "current")
    }
    for phase in ("baseline", "current"):
        root = tmp_path / phase
        root.mkdir(parents=True)
        for name in ("pyproject.toml", "uv.lock", ".python-version", "justfile"):
            shutil.copyfile(ROOT / name, root / name)
        # Retain every runtime/configuration input. Operational guides link to
        # the full repository, which this bounded generated-README fixture lacks.
        shutil.copytree(ROOT / "tooling", root / "tooling", ignore=shutil.ignore_patterns("__pycache__", "*.md"))
        for pattern in declaration(ROOT).harness:
            for source in ROOT.glob(pattern):
                destination = root / source.relative_to(ROOT)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
        (root / "src").mkdir()
        (root / "src/lib.rs").write_text(f"// {phase}\n", encoding="utf-8", newline="\n")
        (root / "README.md").write_text("# Fixture\n\n<!-- BENCHMARKS:BEGIN -->\nOld table\n<!-- BENCHMARKS:END -->\n", encoding="utf-8", newline="\n")
        (root / ".gitignore").write_text("/target/\n", encoding="utf-8", newline="\n")
        manifest = root / "Cargo.toml"
        contents = manifest.read_text(encoding="utf-8")
        version = tomllib.loads(contents)["package"]["version"]
        manifest.write_text(
            contents.replace(f'version = "{version}"', f'version = "{(pair.baseline if phase == "baseline" else pair.current).removeprefix("v")}"', 1),
            encoding="utf-8",
            newline="\n",
        )
        (root / "tooling/fixture.py").write_text(f"IDS = {ids!r}\nTIMINGS = {synthetic_timings()!r}\n" + DRIVER, encoding="utf-8", newline="\n")
        driver = f"{shlex.quote(sys.executable)} -X utf8 tooling/fixture.py"
        native = 'export CRITERION_HOME := absolute_path(env("CRITERION_HOME", "target/criterion"))\n'
        plan = declaration(ROOT)
        assert isinstance(plan, SuiteMeasurementPlan)
        for suite in plan.suites:
            command = suite.measurement.current.command[-1]
            native += f'\n{command} mode="measure":\n    {driver} source {suite.name} {{{{mode}}}}\n'
        for side in ("baseline", "current"):
            native += f'\n{side}-vs_linalg mode="measure":\n    {driver} {side} vs_linalg {{{{mode}}}}\n'
        for gate in sorted({phase.gate[-1] for suite in plan.suites for phase in (suite.measurement.baseline, suite.measurement.current)}):
            native += f"\n{gate}:\n    {driver} source all gate\n"
        (root / "tooling/performance.just").write_text(native, encoding="utf-8", newline="\n")
        for config in (root / "tooling/measurements").glob("*.toml"):
            config.write_text(
                config.read_text(encoding="utf-8").replace('"tooling/measurements/*.toml"', '"tooling/measurements/*.toml", "tooling/fixture.py"'),
                encoding="utf-8",
                newline="\n",
            )
        run_command("git", ["init", "--quiet"], cwd=root)
        run_command("git", ["add", "."], cwd=root)
        run_command(
            "git",
            [
                "-c",
                "user.name=Fixture",
                "-c",
                "user.email=fixture@example.invalid",
                "-c",
                "commit.gpgsign=false",
                "commit",
                "--allow-empty",
                "--quiet",
                "-m",
                phase,
            ],
            cwd=root,
        )
    baseline, current = tmp_path / "baseline", tmp_path / "current"
    run_command("git", ["fetch", str(baseline), "HEAD"], cwd=current)
    run_command("git", ["-c", "tag.gpgsign=false", "tag", pair.baseline, "FETCH_HEAD"], cwd=current)
    (current / "src/untracked.rs").write_bytes(b"// included working-tree input\n")
    return baseline, current


@pytest.fixture(scope="module")
def measured_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Evidence]:
    baseline, current = prepared(tmp_path_factory.mktemp("performance"))
    asserted = [
        PAIR.current,
        run_command("git", ["--no-pager", "rev-parse", "HEAD"], cwd=current).stdout.strip(),
        PAIR.baseline,
        run_command("git", ["--no-pager", "rev-parse", "HEAD"], cwd=baseline).stdout.strip(),
    ]
    (current / "fixture-identities.json").write_text(json.dumps(asserted), encoding="utf-8", newline="\n")
    progress = io.StringIO()
    events: list[dict[str, object]] = []
    evidence = measure_pair(
        current,
        declaration(current),
        PAIR,
        working_tree=True,
        progress=ProgressReporter(stream=progress, observer=events.append),
    )
    (current / "target").mkdir(exist_ok=True)
    (current / "target/progress.json").write_text(json.dumps(events), encoding="utf-8", newline="\n")
    (current / "target/progress.txt").write_text(progress.getvalue(), encoding="utf-8", newline="\n")
    for path in (ROOT / "tooling").glob("performance-*.toml"):
        shutil.copyfile(path, current / "tooling" / path.name)
    shutil.copytree(ROOT / "tooling/plots", current / "tooling/plots", dirs_exist_ok=True)
    (current / "README.md").write_text("# Fixture\n\n<!-- BENCHMARKS:BEGIN -->\nOld table\n<!-- BENCHMARKS:END -->\n", encoding="utf-8", newline="\n")
    publish_evidence(evidence, current / "target/bench-reports/full/run.json", current / "target/bench-reports/full/evidence.json")
    return current, evidence


def test_shared_progress_reports_validated_inventory_and_unknown_eta(measured_run: tuple[Path, Evidence]) -> None:
    """Exercise the installed reporter with the real consumer plan and captured logs."""
    root, evidence = measured_run
    events = json.loads((root / "target/progress.json").read_text(encoding="utf-8"))
    budget = next(event for event in events if event["event"] == "budget")
    library, peers = scientific_ids()
    baseline = baseline_ids() | peers
    total = len(baseline) + len(library)
    total += len(peers)
    assert total == 2689
    assert budget["total"] == budget["remaining"] == budget["unknown_cases"] == total
    assert budget["completed"] == 0
    assert budget["eta_seconds"] is None
    assert budget["nominal_seconds"] == total * 8
    assert sum(row["total"] for row in budget["inventory"]) == total
    accepted = [event for event in events if event["event"] == "validated"]
    assert accepted[-1]["completed"] == total
    assert accepted[-1]["remaining"] == 0
    assert all(event["completed"] + event["remaining"] == total for event in accepted)
    assert accepted[-1]["eta_seconds"] == 0
    assert {event["state"] for event in events} >= {"discovery", "gate", "acquisition"}
    assert all(event["schema"] == "research-repo-tools/acquisition-progress/v1" for event in events)
    log = (root / "target/progress.txt").read_text(encoding="utf-8")
    assert "ETA unknown" in log
    assert f"{total} validated, 0 remaining" in log
    assert "\x1b" not in log
    assert "\r" not in log
    validate_run_policy(evidence, declaration(root))
    trial = select_trial(evidence, 1)
    for dimension in (2, 3, 4, 5, 8, 16, 32, 64):
        library_id = f"d{dimension}/la_stack_lu_solve"
        assert comparison_incompatibilities(trial, ("baseline", library_id), ("current", library_id)) == ()
        for peer in ("nalgebra", "faer"):
            assert comparison_incompatibilities(trial, ("current", f"d{dimension}/{peer}_lu_solve"), ("current", library_id)) == ()


@pytest.mark.parametrize("baseline", ["v0.5.0", "v0.5.2", "v0.5.3"])
def test_api_availability_uses_original_source_receipts(tmp_path: Path, baseline: str) -> None:
    """The common harness must not change which original revision supports an API."""
    pair = ReleasePair("v0.5.4", baseline)
    _baseline, current = prepared(tmp_path, pair)
    plan = declaration(current)
    assert isinstance(plan, SuiteMeasurementPlan)
    # Shared suites have complete-inventory coverage above. Isolate all four
    # introduction boundaries here while retaining every one of their 149 cases.
    selected = replace(plan, suites=tuple(suite for suite in plan.suites if suite.since))
    evidence = measure_pair(current, selected, pair, working_tree=True, progress=ProgressReporter(mode="none"))
    run = validate_run_policy(evidence, selected)
    assert isinstance(run, ScheduledRun)
    phases = dict(parse_suite_run(select_trial(evidence, 1)).phases)
    introduced = scientific_ids()[0] - baseline_ids(release="v0.5.1")
    assert len(introduced) == 149
    assert set(phases["current"].policy.expected) == introduced
    assert set(phases["baseline"].policy.expected) == (set() if baseline == "v0.5.0" else introduced)


@pytest.mark.parametrize("phase", ["baseline", "current"])
def test_full_plan_rejects_unsupported_release_before_native_commands(tmp_path: Path, phase: str) -> None:
    """Schema 5 captures source receipts before enforcing the native release floor."""
    pair = ReleasePair("v0.4.6" if phase == "current" else PAIR.current, "v0.4.6" if phase == "baseline" else PAIR.baseline)
    _baseline, current = prepared(tmp_path, pair)
    events: list[dict[str, object]] = []
    with pytest.raises(ValueError, match=r"requires v0\.5\.0 or newer"):
        measure_pair(current, declaration(current), pair, working_tree=True, progress=ProgressReporter(mode="none", observer=events.append))
    assert not any(event["state"] in {"discovery", "gate", "acquisition"} for event in events)


def test_supported_new_api_discovery_failure_is_not_treated_as_unavailable(tmp_path: Path) -> None:
    """Once column bases exist, failed native compilation cannot remove coverage."""
    pair = ReleasePair("v0.5.3", "v0.5.2")
    _baseline, current = prepared(tmp_path, pair)
    plan = declaration(current)
    assert isinstance(plan, SuiteMeasurementPlan)
    suite = next(item for item in plan.suites if item.name == "column-basis")
    assert suite.since == pair.baseline == "v0.5.2"
    policy = suite.measurement.baseline.policy
    assert isinstance(policy, DiscoveryPolicy)
    assert len(policy.discovery.required) == 55
    failure = (sys.executable, "-c", 'import sys; sys.stderr.write("fixture column-basis compilation failed\\n"); sys.exit(7)')
    discovery = replace(policy.discovery, suites=(replace(policy.discovery.suites[0], command=failure),))
    baseline = replace(suite.measurement.baseline, policy=replace(policy, discovery=discovery))
    selected = replace(plan, suites=(replace(suite, measurement=replace(suite.measurement, baseline=baseline)),))
    progress = io.StringIO()
    with pytest.raises(subprocess.CalledProcessError) as raised:
        measure_pair(current, selected, pair, working_tree=True, progress=ProgressReporter(mode="json", stream=progress))
    assert raised.value.returncode == 7
    assert b"fixture column-basis compilation failed" in raised.value.stderr
    events = [json.loads(line) for line in progress.getvalue().splitlines()]
    # Discovery failed before an observed inventory could establish the budget.
    assert not any(event["event"] == "budget" for event in events)
    assert {key: events[-1][key] for key in ("event", "state", "suite", "phase", "completed")} == {
        "event": "failed",
        "state": "discovery",
        "suite": "column-basis",
        "phase": "baseline",
        "completed": 0,
    }
    assert not any(event["event"] == "validated" for event in events)
    assert not any(event["state"] in {"gate", "acquisition"} for event in events)


@pytest.mark.parametrize("failed_stage", ["gate", "acquisition"])
def test_shared_progress_keeps_early_failure_unaccepted(tmp_path: Path, failed_stage: str) -> None:
    """Supported native failures remain failures, never completed progress."""
    _baseline, current = prepared(tmp_path)
    plan = declaration(current)
    assert isinstance(plan, SuiteMeasurementPlan)
    failure = (sys.executable, "-c", "raise SystemExit(7)")
    suite = plan.suites[0]
    phase = replace(suite.measurement.baseline, **{("gate" if failed_stage == "gate" else "command"): failure})
    plan = replace(plan, suites=(replace(suite, measurement=replace(suite.measurement, baseline=phase)), *plan.suites[1:]))
    progress = io.StringIO()
    if failed_stage == "gate":
        with pytest.raises(ValueError, match=r"^baseline preflight gate failed$"):
            measure_pair(current, plan, PAIR, working_tree=True, progress=ProgressReporter(mode="json", stream=progress))
    else:
        with pytest.raises(subprocess.CalledProcessError) as raised:
            measure_pair(current, plan, PAIR, working_tree=True, progress=ProgressReporter(mode="json", stream=progress))
        assert raised.value.returncode == 7
    events = [json.loads(line) for line in progress.getvalue().splitlines()]
    failed = events[-1]
    assert failed["event"] == "failed"
    assert failed["state"] == failed_stage
    assert failed["phase"] == "baseline"
    assert failed["completed"] == 0
    assert not any(event["event"] == "validated" for event in events)
    assert not any(event["event"] == "finished" and event["state"] == failed_stage and event["suite"] == suite.name for event in events)
    assert not (current / "target/bench-reports/full/run.json").exists()


@pytest.fixture
def complete_run(tmp_path: Path, measured_run: tuple[Path, Evidence]) -> tuple[Path, Evidence]:
    root, evidence = measured_run
    destination = tmp_path / "current"
    shutil.copytree(root, destination)
    return destination, evidence


def identities(root: Path) -> list[str]:
    """Read assertions captured from the fixture's Git sources before acquisition."""
    return cast("list[str]", json.loads((root / "fixture-identities.json").read_text(encoding="utf-8")))


@pytest.mark.release_candidate
def test_local_candidate_capture_and_final_acceptance_reuse_original_bytes(
    complete_run: tuple[Path, Evidence],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Stub only the separately tested native CI provider; preserve real candidate binding."""
    root, _evidence = complete_run
    current_tag, current_revision, baseline_tag, baseline_revision = identities(root)
    record = serialize_preparation(
        Preparation(
            "acgetchell/la-stack",
            "tooling/measurements/full.toml",
            hashlib.sha256((root / "tooling/measurements/full.toml").read_bytes()).hexdigest(),
            current_tag.removeprefix("v"),
            ReleasePair(current_tag, baseline_tag),
            current_revision,
            baseline_revision,
            working_tree=True,
        )
    )
    (root / "target/bench-reports/release-inputs.json").write_bytes(record)
    owner = tmp_path / "independent owner"
    owner.mkdir()
    freeze(owner, {"preparation.json": record})

    def no_acquisition(*_args: object, **_kwargs: object) -> None:
        pytest.fail("retained candidate consumers must never acquire or preflight")

    monkeypatch.setattr("research_repo_tools.common_measurement.run_command_live", no_acquisition)
    with pytest.raises(FileNotFoundError, match=r"run\.json"):
        capture(root, owner)
    assertions = capture(root, owner, retained_scratch=True)
    assert (owner / "run.json").read_bytes() == (root / "target/bench-reports/full/run.json").read_bytes()
    assert (owner / "evidence.json").read_bytes() == (root / "target/bench-reports/full/evidence.json").read_bytes()
    original = {path.name: path.read_bytes() for path in owner.iterdir()}
    assert_owner_capture(root, owner, assertions, original, monkeypatch)
    publish_publication(
        load_report_plan(
            root,
            "tooling/performance-report.toml",
            payload="target/bench-reports/full/run.json",
            manifest="target/bench-reports/full/evidence.json",
        )
    )
    result = recipe(root, "performance-readme", *identities(root))
    assert result.returncode == 0, result.stdout + result.stderr
    for args in (
        ["add", "."],
        ["-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "final outputs"],
    ):
        run_command("git", args, cwd=root)
    commit = run_command("git", ["--no-pager", "rev-parse", "HEAD"], cwd=root).stdout.strip()
    native = deterministic_json(
        {
            "package": {"name": "la-stack", "version": current_tag.removeprefix("v")},
            "policy": {"repository": "acgetchell/la-stack", "guarantee": "signed-source", "event": "push", "ref": "refs/heads/main", "signing_job": "sign"},
            "run": {"head_sha": commit, "id": 123, "run_attempt": 1, "path": ".github/workflows/ci.yml"},
            "artifact": {"id": 456, "digest": "sha256:" + "a" * 64},
        }
    )
    expected_selections: list[bytes | None] = []

    def native_evidence(*_args: object, expected: bytes | None = None, **_kwargs: object) -> bytes:
        expected_selections.append(expected)
        if expected is not None:
            assert expected == native
        return native

    monkeypatch.setattr("research_repo_tools.candidate_reuse.verify_release_evidence", native_evidence)
    monkeypatch.setattr("research_repo_tools.dependency_provenance.check_release_dependencies", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("tooling.release_candidate.check_release_dependencies", lambda *_args, **_kwargs: None)
    accepted = verify(root, owner)
    original_acceptance = accepted.read_bytes()
    assert verify(root, owner).read_bytes() == original_acceptance
    assert expected_selections == [None, native]
    artifact = archive(owner, commit)
    original_archive = artifact.read_bytes()
    assert archive(owner, commit).read_bytes() == original_archive
    assert {name: (owner / name).read_bytes() for name in original} == original
    assert json.loads(original_acceptance)["binding"]["candidate"]["revision"] == current_revision
    assert json.loads(original_acceptance)["binding"]["final"]["revision"] == commit
    assert_hosted_reuse(root, owner, accepted, artifact, monkeypatch)
    assert expected_selections == [None, native, native]
    # Changed receipt bytes fail against the separately frozen digest.
    (owner / "publication.json").write_bytes(original["publication.json"] + b" ")
    with pytest.raises(RuntimeError, match="candidate verification failed"):
        verify(root, owner)
    assert accepted.read_bytes() == original_acceptance
    assert artifact.read_bytes() == original_archive


def assert_owner_capture(root: Path, owner: Path, assertions: dict[str, str], original: dict[str, bytes], monkeypatch: pytest.MonkeyPatch) -> None:
    """Retain the actual captured assertions while registering automatic lookup."""
    assert capture(root, owner) == assertions
    assert {path.name: path.read_bytes() for path in owner.iterdir()} == original
    with monkeypatch.context() as context:
        context.setattr(Path, "home", lambda: owner.parent / "state home")
        remember_owner(root, owner)
        assert selected_owner(root) == owner.resolve()


def assert_hosted_reuse(root: Path, owner: Path, accepted: Path, artifact: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify the real archive and final binding through the hosted consumer."""
    tmp_path = owner.parent
    original_acceptance, original_archive = accepted.read_bytes(), artifact.read_bytes()
    policy = parse_candidate_policy((accepted.parent / "policy.json").read_bytes())
    current_tag, commit = policy.tag, policy.commit
    assert commit is not None
    native = deterministic_json(json.loads(original_acceptance)["native"])
    assertions = {
        "acquisition_sha256": (owner / "acquisition.sha256").read_text(encoding="utf-8").strip(),
        "publication_sha256": (owner / "publication.sha256").read_text(encoding="utf-8").strip(),
    }
    incoming = tmp_path / "transport with spaces"
    incoming.mkdir()
    files = {
        "policy": accepted.parent / "policy.json",
        "publication": owner / "publication.json",
        "acceptance": accepted,
        "archive": artifact,
    }
    assets = {}
    for identifier, (role, source) in enumerate(files.items(), 1):
        name = "archive.tar.gz" if role == "archive" else f"{role}.json"
        data = source.read_bytes()
        (incoming / name).write_bytes(data)
        assets[role] = {"name": name, "id": identifier, "sha256": hashlib.sha256(data).hexdigest()}
    transfer = json.dumps(
        {
            "schema": 1,
            "target": json.loads(serialize_release_target(ReleaseTarget("acgetchell/la-stack", current_tag, 42, commit, current_tag))),
            "acquisition": assertions["acquisition_sha256"],
            "publication": assertions["publication_sha256"],
            "assets": assets,
        }
    ).encode()
    for key, value in {
        "GITHUB_EVENT_NAME": "workflow_dispatch",
        "GITHUB_REPOSITORY": "acgetchell/la-stack",
        "GITHUB_REF": f"refs/tags/{current_tag}",
        "GITHUB_SHA": commit,
        "RELEASE_TAG": current_tag,
        "GITHUB_WORKFLOW_REF": f"acgetchell/la-stack/.github/workflows/release-benchmarks.yml@refs/tags/{current_tag}",
        "GITHUB_ACTOR": "acgetchell",
        "GITHUB_TRIGGERING_ACTOR": "acgetchell",
    }.items():
        monkeypatch.setenv(key, value)
    received = tmp_path / "hosted acceptance"
    assert receive(root, received, transfer, source=incoming).read_bytes() == original_archive
    assert (received / "acceptance.json").read_bytes() == original_acceptance
    (incoming / "acceptance.json").write_bytes(original_acceptance + b" ")
    with pytest.raises(ValueError, match="transported candidate bytes differ"):
        receive(root, received, transfer, source=incoming)
    assert (received / "acceptance.json").read_bytes() == original_acceptance
    (incoming / "acceptance.json").write_bytes(original_acceptance)

    def changed_native(*_args: object, expected: bytes | None = None, **_kwargs: object) -> bytes:
        assert expected == native
        message = "selected run or attempt changed"
        raise ValueError(message)

    with pytest.MonkeyPatch.context() as changed:
        changed.setattr("research_repo_tools.candidate_reuse.verify_release_evidence", changed_native)
        with pytest.raises(ValueError, match="selected run or attempt changed"):
            receive(root, received, transfer, source=incoming)
    assert (received / "acceptance.json").read_bytes() == original_acceptance
    assert_published_renewal(root, files, transfer, monkeypatch)


def assert_published_renewal(root: Path, files: dict[str, Path], transfer: bytes, monkeypatch: pytest.MonkeyPatch) -> None:
    """Renew real candidate acceptance while retaining every original release asset."""
    owner = files["publication"].parent
    original = {role: path.read_bytes() for role, path in files.items()}
    data = json.loads(transfer)
    tag = data["target"]["tag"]
    data["assets"]["attached"] = {"name": f"la-stack-{tag}-criterion-baseline.tar.gz", "id": 9, "sha256": data["assets"]["archive"]["sha256"]}
    original["attached"] = original["archive"]
    native = json.loads(original["acceptance"])["native"]
    native["run"]["run_attempt"] = 2
    native["artifact"] = {"id": 789, "digest": "sha256:" + "b" * 64}
    renewed_native = deterministic_json(native)
    selections: list[bytes | None] = []

    def native_evidence(*_args: object, expected: bytes | None = None, **_kwargs: object) -> bytes:
        selections.append(expected)
        if expected is not None:
            assert expected == renewed_native
        return renewed_native

    def published(*_args: object, **_kwargs: object) -> GitHubRelease:
        return GitHubRelease(
            repository="acgetchell/la-stack",
            tag=tag,
            identifier=42,
            draft=False,
            immutable=True,
            prerelease=False,
            title=tag,
            assets=tuple((item["name"], item["id"], len(original[role]), item["sha256"]) for role, item in data["assets"].items()),
        )

    def download(_root: Path, _repository: str, _tag: str, name: str, destination: Path, *, expected_sha256: str, **_kwargs: object) -> None:
        role = next(role for role, item in data["assets"].items() if item["name"] == name)
        assert hashlib.sha256(original[role]).hexdigest() == expected_sha256
        destination.write_bytes(original[role])

    with monkeypatch.context() as renewed:
        renewed.setattr("research_repo_tools.candidate_reuse.verify_release_evidence", native_evidence)
        renewed.setattr("tooling.release_candidate.check_reviewed_release", published)
        renewed.setattr("tooling.release_candidate.download_release_asset", download)
        renewed.setenv("GITHUB_WORKFLOW_REF", f"acgetchell/la-stack/.github/workflows/publish-crates.yml@refs/tags/{tag}")
        reviewed = verify(root, owner, "renewed-ci").read_bytes()
        assert reviewed != original["acceptance"]
        assert json.loads(reviewed)["binding"] == json.loads(original["acceptance"])["binding"]
        digest = hashlib.sha256(reviewed).hexdigest()
        assertions = json.dumps(data).encode()
        before, after = owner.parent / "before approval", owner.parent / "after approval"
        assert receive(root, before, assertions, source=ReviewedAcceptance(digest)).read_bytes() == original["archive"]
        frozen = (before / "acceptance.json").read_bytes()
        assert frozen == reviewed
        assert receive(root, after, assertions, source=ReviewedAcceptance(digest, frozen)).read_bytes() == original["archive"]
        assert (after / "acceptance.json").read_bytes() == reviewed
    assert selections == [None, None, renewed_native]
    assert {role: path.read_bytes() for role, path in files.items()} == {role: original[role] for role in files}


def test_working_tree_inputs_include_nonignored_untracked_consumer_sources(complete_run: tuple[Path, Evidence]) -> None:
    root, evidence = complete_run
    snapshot = capture_snapshot(root)
    files = {name: data for name, data, _mode in snapshot.untracked}
    assert files["src/untracked.rs"] == b"// included working-tree input\n"
    assert not any(name.startswith("target/") for name in files)
    observations = [source for _suite, phase, _sample, source in parse_suite_run(select_trial(evidence, 1)).observations if phase == "current"]
    assert len(observations) == 7
    for source in observations:
        assert "src/untracked.rs" in json.loads(dict(source.context)["source-inventory"])


def test_actual_non_exact_measurement_recipe_discovers_and_retains_complete_samples(tmp_path: Path) -> None:
    _baseline, root = prepared(tmp_path)
    shutil.copytree(ROOT / "tooling/plots/non-exact", root / "tooling/plots/non-exact", dirs_exist_ok=True)
    manifest = root / "Cargo.toml"
    current = manifest.read_text(encoding="utf-8")
    manifest.write_text(
        current.replace(f'version = "{PAIR.current.removeprefix("v")}"', f'version = "{PAIR.baseline.removeprefix("v")}"', 1), encoding="utf-8", newline="\n"
    )
    run_command("git", ["-c", "tag.gpgsign=false", "tag", "-d", PAIR.baseline], cwd=root)
    for arguments in (
        ["add", "."],
        ["-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "baseline"],
        ["-c", "tag.gpgsign=false", "tag", PAIR.baseline],
    ):
        run_command("git", arguments, cwd=root)
    manifest.write_text(current, encoding="utf-8", newline="\n")
    (root / "src/lib.rs").write_text("// current changed source\n", encoding="utf-8", newline="\n")
    for arguments in (
        ["add", "."],
        ["-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "current"],
        ["-c", "tag.gpgsign=false", "tag", PAIR.current],
    ):
        run_command("git", arguments, cwd=root)
    asserted = [
        PAIR.current,
        run_command("git", ["--no-pager", "rev-parse", f"{PAIR.current}^{{commit}}"], cwd=root).stdout.strip(),
        PAIR.baseline,
        run_command("git", ["--no-pager", "rev-parse", f"{PAIR.baseline}^{{commit}}"], cwd=root).stdout.strip(),
    ]
    result = recipe(root, "performance-measure", PAIR.current, PAIR.baseline, "non-exact")
    assert result.returncode == 0, result.stdout + result.stderr
    evidence = load_evidence(root / "target/bench-reports/non-exact/run.json", root / "target/bench-reports/non-exact/evidence.json")
    plan = declaration(root, "non-exact")
    assert isinstance(plan, CompleteMeasurementPlan)
    run = validate_run_policy(evidence, plan)
    assert all(name.startswith("d") for _phase, sample in run.phases for name in sample.policy.expected)
    assert (root / "target/bench-reports/non-exact/report.md").is_file()
    assert not (root / "target/bench-reports/full").exists()
    result = recipe(root, "plot-vs-linalg", *asserted, "lu_solve-median-log", "non-exact")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (root / "docs/assets/bench/vs_linalg_non_exact_lu_solve_median_log.csv").is_file()
    assert not (root / "docs/assets/bench/vs_linalg_lu_solve_median_log.csv").exists()


@pytest.mark.parametrize(("statistic", "axis", "source"), [("mean", "linear", "scratch"), ("median", "log", "retained")])
def test_actual_retention_offline_plot_and_readme_recipes(complete_run: tuple[Path, Evidence], statistic: str, axis: str, source: str) -> None:
    root, _evidence = complete_run
    result = recipe(root, "performance-promote")
    assert result.returncode == 0, result.stdout + result.stderr
    if source == "retained":
        shutil.rmtree(root / "target")
        result = recipe(root, "performance-doc")
        assert result.returncode == 0, result.stdout + result.stderr
    result = recipe(root, "plot-vs-linalg", *identities(root), f"lu_solve-{statistic}-{axis}", source)
    assert result.returncode == 0, result.stdout + result.stderr
    csv_text = (root / f"docs/assets/bench/vs_linalg_lu_solve_{statistic}_{axis}.csv").read_text(encoding="utf-8")
    csv_rows = list(csv.DictReader(io.StringIO(csv_text)))
    assert len(csv_rows) == 24
    assert {row["series"] for row in csv_rows} == {"lu-solve/la-stack current", "lu-solve/nalgebra", "lu-solve/faer"}
    assert {row["phase"] for row in csv_rows} == {"current"}
    assert {row["statistic"] for row in csv_rows} == {statistic}
    if source == "retained":
        prepared_inputs = root / "target/bench-reports/release-inputs.json"
        prepared_inputs.parent.mkdir(parents=True)
        current_tag, current_revision, baseline_tag, baseline_revision = identities(root)
        prepared_inputs.write_bytes(
            serialize_preparation(
                Preparation(
                    "acgetchell/la-stack",
                    "tooling/measurements/full.toml",
                    hashlib.sha256((root / "tooling/measurements/full.toml").read_bytes()).hexdigest(),
                    current_tag.removeprefix("v"),
                    ReleasePair(current_tag, baseline_tag),
                    current_revision,
                    baseline_revision,
                    working_tree=True,
                )
            )
        )
    result = recipe(root, "performance-readme", *([] if source == "retained" else identities(root)))
    assert result.returncode == 0, result.stdout + result.stderr
    document = (root / "README.md").read_text(encoding="utf-8")
    assert "Speedup vs nalgebra" in document
    assert "Speedup vs faer" in document
    assert "Complete run evidence" in document
    assert "Evidence envelope" in document
    assert "la-stack v0.5.2" in document
    assert "Statistic: median. Unit: ns." in document
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    rows = [
        [cell.strip() for cell in line.strip("|").split("|")] for line in document.splitlines() if line.startswith("|") and line.split("|")[1].strip().isdigit()
    ]
    assert [row[0] for row in rows] == [str(d) for d in catalog["dimensions"]]
    timings = synthetic_timings()
    for dimension, row in zip(catalog["dimensions"], rows, strict=True):
        assert len(row) == 6
        expected = [
            timings[phase][f"d{dimension}/{benchmark}"]["median"]
            for phase, benchmark in zip(("current", "current", "current"), catalog["metrics"]["lu_solve"]["benchmarks"], strict=True)
        ]
        assert [float(cell.split("<br>")[0]) for cell in row[1:4]] == pytest.approx(expected)
        # Published ratios have six significant digits; the oracle uses the
        # intended peer/current medians, independently of the view declaration.
        assert [float(cell.removesuffix("\N{MULTIPLICATION SIGN}")) for cell in row[4:]] == pytest.approx(
            [expected[1] / expected[0], expected[2] / expected[0]], rel=1e-5
        )
    checked = recipe(root, "markdown-check")
    assert checked.returncode == 0, checked.stdout + checked.stderr


@pytest.mark.parametrize("failure", ["partial-scratch", "stale-source", "stale-harness", "formatter"])
def test_publication_boundaries_leave_all_outputs_unchanged(
    complete_run: tuple[Path, Evidence],
    failure: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root, _evidence = complete_run
    result = recipe(root, "performance-promote")
    assert result.returncode == 0, result.stdout + result.stderr
    if failure == "partial-scratch":
        (root / "target/bench-reports/full/evidence.json").unlink()
    elif failure == "stale-source":
        (root / "src/lib.rs").write_text("// changed\n", encoding="utf-8", newline="\n")
    elif failure == "stale-harness":
        (root / "tooling/performance.just").write_text("# changed\n", encoding="utf-8", newline="\n")
    outputs = [root / "README.md", *(root / "docs").rglob("*")]
    before = {path: path.read_bytes() for path in outputs if path.is_file()}
    if failure == "formatter":
        # Frozen publisher configuration must not change after measurement.
        # Simulate a transient formatter process failure at its real boundary.
        def failed_formatter(*_args: object, **_kwargs: object) -> None:
            raise subprocess.CalledProcessError(7, ["rumdl"], stderr="fixture formatter failure")

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr("research_repo_tools.postprocess_changelog.run_safe_command", failed_formatter)
            current, revision, baseline, previous = identities(root)
            status = main(
                [
                    "--root",
                    str(root),
                    "performance",
                    "readme",
                    "tooling/performance-release.toml",
                    "--current-release",
                    current,
                    "--current-revision",
                    revision,
                    "--baseline-release",
                    baseline,
                    "--baseline-revision",
                    previous,
                ]
            )
        diagnostic = capsys.readouterr()
        result = subprocess.CompletedProcess([], status, diagnostic.out, diagnostic.err)
    else:
        result = recipe(root, "performance-promote") if failure == "partial-scratch" else recipe(root, "performance-readme", *identities(root))
    assert result.returncode != 0
    diagnostics = {
        "stale-source": "current core-default source-inventory differs from measured inputs",
        "stale-harness": "current core-default harness-inventory differs from measured inputs",
        "formatter": "rumdl could not format README.md",
    }
    if failure in diagnostics:
        assert diagnostics[failure] in result.stderr
    if failure == "partial-scratch":
        result = recipe(root, "plot-vs-linalg", *identities(root), "lu_solve-median-log", "scratch")
        assert result.returncode != 0
    assert {path: path.read_bytes() for path in before} == before
    assert not (root / "docs/assets").exists()
    if failure == "formatter":
        # Retry publication after repairing its formatter. Archived measurements
        # remain byte-identical; no timing or discovery command is involved.
        retained = {path: path.read_bytes() for path in (root / "docs/performance-runs").rglob("*") if path.is_file()}
        result = recipe(root, "performance-readme", *identities(root))
        assert result.returncode == 0, result.stdout + result.stderr
        assert {path: path.read_bytes() for path in retained} == retained


@pytest.mark.parametrize("phase", ["baseline", "current"])
def test_actual_measurement_recipe_rejects_old_releases_before_git(phase: str, tmp_path: Path) -> None:
    root = tmp_path / "consumer"
    root.mkdir()
    for name in ("pyproject.toml", "uv.lock", ".python-version", "Cargo.toml", "rust-toolchain.toml"):
        shutil.copyfile(ROOT / name, root / name)
    for pattern in declaration(ROOT).harness:
        for source in ROOT.glob(pattern):
            destination = root / source.relative_to(ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
    (root / "src").mkdir()
    (root / "src/lib.rs").write_text("// fixture\n", encoding="utf-8", newline="\n")
    result = recipe(
        root, "performance-measure", "v0.4.6" if phase == "current" else PAIR.current, "v0.4.6" if phase == "baseline" else PAIR.baseline, "non-exact"
    )
    assert result.returncode != 0
    assert "requires v0.5.0 or newer" in result.stderr
    assert not (root / ".git").exists()
    assert not (root / "target").exists()


def test_legacy_bytes_hashes_are_unchanged() -> None:
    directory = ROOT / "docs/performance/v0.4.6-vs-v0.4.5/39719a1eb599d294dbbc1e46dc6fb0de2028ad2e59aab5f94fb48c7bf3d7f841"
    hashes = {
        "performance.csv": "b263bb34bcfd43b8fc8109110db9f2e908516aef6c01ce16f70f5ac18dc7ef2d",
        "performance.full.csv": "4fa5cc13159f5dc91e711ae67554845bcbd1ace2cfc3e31487fc85f75d3dc459",
        "performance.full.provenance.json": "a43f66e52035aecd8d235dfded916b23c3240529d0c07f793d27322521cfb0d0",
        "performance.provenance.json": "89144a2899ba73f55d009da4fa20b4e6477600902b93cde469383f719c497508",
    }
    assert {name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in hashes} == hashes


def test_pre_declaration_complete_run_is_preserved_as_a_frozen_record() -> None:
    directory = ROOT / "docs/performance-runs/runs/634bc9527c2caf0a431a9bdfac082c190e1094f0b29eb6e3f28c4b96961cdec5"
    hashes = {
        "evidence.json": "6565bbb489662b68df8cbb039947a35de9ecdec0b25aa7f4de03dcf76bf8ed9d",
        "report.md": "ac5ca08784750ea04211dc3142e8b33ecf01fe7c3520d6d5de7cab78187235a7",
        "run.json": "9b860d8a56a4313d60031039fd518f142c1aa2e74933fc9e73fb128a00023b2a",
    }
    evidence = load_evidence(directory / "run.json", directory / "evidence.json")
    with pytest.raises(ValueError, match=r"scheduled run fields must be \['declaration', 'schema', 'trials'\]"):
        validate_run_policy(evidence, declaration(ROOT), release=True)
    assert {name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in hashes} == hashes


def test_shared_saved_comparison_recipe_reports_coverage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Exercise the actual recipe and installed CLI, with no consumer package."""
    for name in ("justfile", "pyproject.toml", "uv.lock", ".python-version", "rust-toolchain.toml"):
        shutil.copyfile(ROOT / name, tmp_path / name)
    for name, sample, value in (("workload/shared", "last", 10), ("workload/shared", "new", 8), ("diagnostic/new", "new", 3)):
        path = tmp_path / "target/criterion" / name / sample
        path.mkdir(parents=True)
        (path / "benchmark.json").write_text(json.dumps({"full_id": name}), encoding="utf-8")
        estimate = {"point_estimate": value, "confidence_interval": {"lower_bound": value * 0.9, "upper_bound": value * 1.1}}
        (path / "estimates.json").write_text(json.dumps({"median": estimate}), encoding="utf-8")
    monkeypatch.setenv("UV_PROJECT_ENVIRONMENT", sys.prefix)
    result = run_command("just", ["--no-deps", "bench-compare", "last"], cwd=tmp_path, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    report = (tmp_path / "target/bench-reports/performance.md").read_text(encoding="utf-8")
    assert "workload/shared" in report
    assert "diagnostic/new" in report
    assert not (tmp_path / "target/bench-reports/performance.csv").exists()
    assert not (tmp_path / "target/bench-reports/performance.provenance.json").exists()
