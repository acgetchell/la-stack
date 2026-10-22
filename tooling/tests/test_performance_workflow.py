"""Installed workflow acceptance with independent la-stack scientific declarations."""

import hashlib
import json
import os
import shlex
import shutil
import sys
import tomllib
from pathlib import Path
from typing import TYPE_CHECKING, cast

import pytest
from research_repo_tools.common_measurement import CommonHarnessPlan, DiscoveryPolicy, measure_prepared_pair, validate_run_policy
from research_repo_tools.complete_runs import RunSeries
from research_repo_tools.evidence import load_evidence, publish_evidence
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.measurement import load_measurement
from research_repo_tools.process import run_command
from research_repo_tools.release_pairs import ReleasePair
from research_repo_tools.worktrees import capture_snapshot

if TYPE_CHECKING:
    import subprocess

    from research_repo_tools.evidence import Evidence

ROOT = Path(__file__).resolve().parents[2]
PAIR = ReleasePair("v0.5.1", "v0.5.0")
DRIVER = """import json, os, sys
from pathlib import Path
phase, suite, action = sys.argv[1:]
if action == "gate":
    Path("target").mkdir(exist_ok=True)
    with Path("target/gates.txt").open("a", encoding="utf-8", newline="\\n") as stream:
        stream.write("gate\\n")
    raise SystemExit(0)
ids = [name for name in IDS[phase] if (name.startswith("d") if suite == "vs_linalg" else not name.startswith("d"))]
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
    }
    for filename, value in data.items():
        (path / filename).write_text(json.dumps(value), encoding="utf-8", newline="\\n")
"""


def declaration(root: Path, view: str = "full") -> CommonHarnessPlan:
    plan = load_measurement(root, f"tooling/measurements/{view}.toml")
    assert isinstance(plan, CommonHarnessPlan)
    return plan


def scientific_ids() -> tuple[set[str], set[str]]:
    """Derive expectations from the independent domain catalog, not a measured run."""
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    library = {
        *(f"{group}/{name}" for group, names in catalog["exact"].items() for name in names),
        *(
            f"d{dimension}/{name}"
            for dimension in catalog["dimensions"]
            for name in [*catalog["library-benchmarks"], *catalog["additional"].get(str(dimension), [])]
        ),
    }
    peers = {f"{name.split('/')[0]}/{peer}" for name in library for peer in catalog["peers"].get(name.split("/")[1], [])}
    return library, peers


def synthetic_timings() -> dict[str, dict[str, dict[str, float]]]:
    """Distinguish every case, phase, and statistic so selection errors stay visible."""
    library, peers = scientific_ids()
    return {
        phase: {
            name: {"median": float(offset + index * 4), "mean": float(offset + index * 4 + 1)} for index, name in enumerate(sorted(library | peers), start=1)
        }
        for phase, offset in (("current", 0), ("baseline", 10_000))
    }


@pytest.mark.parametrize("view", ["full", "non-exact"])
def test_declared_scientific_sampling_coverage_and_original_peer_phases(view: str) -> None:
    library, peers = scientific_ids()
    if view == "non-exact":
        library = {name for name in library if name.startswith("d")}
    plan = declaration(ROOT, view)
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    assert plan.measurement.minimum_release == "v" + ".".join(str(part) for part in catalog["minimum-release"])
    for phase, required in ((plan.baseline, library | peers), (plan.current, library)):
        assert isinstance(phase.policy, DiscoveryPolicy)
        assert set(phase.policy.discovery.required) == required
        assert phase.policy.sample_count == 100
        assert phase.policy.statistics == ("mean", "median")
        assert phase.policy.confidence_level == 0.95
        assert phase.policy.unit == "ns"
        assert phase.gate[-1] == "gate"
        environment = dict(phase.environment)
        assert environment["RUSTFLAGS"] == environment["CARGO_ENCODED_RUSTFLAGS"] == "--cap-lints=warn"
        assert environment["RUSTUP_TOOLCHAIN"] == tomllib.loads((ROOT / "rust-toolchain.toml").read_text(encoding="utf-8"))["toolchain"]["channel"]
    assert {item.name: item.phase for item in plan.series} == {
        "la-stack baseline": "baseline",
        "la-stack current": "current",
        "nalgebra": "baseline",
        "faer": "baseline",
    }
    for item in plan.series:
        assert isinstance(item, RunSeries)
        if item.name.startswith("la-stack"):
            assert dict(item.rows) == {name: name for name in library}
        else:
            assert dict(item.rows) == {
                f"d{dimension}/{name}": f"d{dimension}/{peer}"
                for dimension in catalog["dimensions"]
                for name, peer_names in catalog["peers"].items()
                for peer in peer_names
                if peer.startswith(f"{item.name}_")
            }
    assert {"context.os", "context.architecture", "context.cpu", "harness_sha256"} <= set(plan.measurement.compatible)
    assert {"benches/**/*.rs", "tests/exact_bench_config.rs", "tooling/performance.just", "tooling/measurements/*.toml"} <= set(plan.harness)


def test_named_views_keep_every_metric_statistic_axis_and_explicit_input_choice() -> None:
    catalog = tomllib.loads((ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
    for source in ("scratch", "retained", "non-exact"):
        files = sorted((ROOT / "tooling/plots" / source).glob("*.toml"))
        assert len(files) == len(catalog["metrics"]) * 4
        for path in files:
            metric, statistic, axis = path.stem.split("-")
            view = tomllib.loads(path.read_text(encoding="utf-8"))
            assert view["statistic"] == statistic
            assert statistic in {"mean", "median"}
            assert view["unit"] == "ns"
            assert view["plot"]["y-scale"] == axis
            assert axis in {"linear", "log"}
            assert [row["coordinate"] for row in view["rows"]] == catalog["dimensions"]
            assert [series["name"] for series in view["series"]] == ["la-stack current", "nalgebra", "faer"]
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
    assert view["statistic"] == "median"
    assert view["unit"] == "ns"
    assert view["rows"] == [{"key": str(d), "label": str(d), "coordinate": d} for d in catalog["dimensions"]]
    assert [series["name"] for series in view["series"]] == ["la-stack current", "nalgebra", "faer"]
    for series, benchmark in zip(view["series"], catalog["metrics"]["lu_solve"]["benchmarks"], strict=True):
        assert series["rows"] == {str(d): f"d{d}/{benchmark}" for d in catalog["dimensions"]}
    assert view["comparisons"] == [
        {"numerator": "nalgebra", "denominator": "la-stack current"},
        {"numerator": "faer", "denominator": "la-stack current"},
    ]


def recipe(root: Path, name: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    """Invoke the actual recipe with the already-installed locked consumer environment."""
    return run_command(
        "just",
        ["--no-deps", "--justfile", str(ROOT / "justfile"), "--working-directory", str(root), name, *arguments],
        env={**os.environ, "UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
        check=False,
    )


def test_native_recipes_preserve_scientific_gates_and_peer_reuse() -> None:
    justfile = Path("tooling/performance.just")
    baseline = dry_run(ROOT, "baseline-all", justfile=justfile).stderr
    current = dry_run(ROOT, "current-all", justfile=justfile).stderr
    assert baseline.count("cargo bench --locked") == current.count("cargo bench --locked") == 2
    assert "la_stack --noplot" in current
    assert "la_stack --noplot" not in baseline
    assert "--target" in baseline
    assert "--target" in current
    assert dry_run(ROOT, "baseline-all", ["list"], justfile=justfile).stderr.count("--list") == 2
    assert "--test vs_linalg_inputs --test exact_bench_config" in dry_run(ROOT, "gate", justfile=justfile).stderr
    for name in ("performance-local", "performance-local-non-exact"):
        command = dry_run(ROOT, name).stderr
        assert "research-repo-tools performance measure" in command
        assert "--mode current-vs-latest" in command
        assert "--allow-git-mutations" in command
    command = dry_run(ROOT, "performance-release", [PAIR.current, PAIR.baseline]).stderr
    assert command.count("performance measure") == command.count("performance promote") == 1
    assert "--payload target/bench-reports/full/run.json" in command
    assert "scripts/" not in command


def prepared(tmp_path: Path) -> tuple[Path, Path]:
    """Replace only native timing/gates with a bounded fixture, retaining real policy."""
    library, peers = scientific_ids()
    ids = {"baseline": sorted(library | peers), "current": sorted(library)}
    for phase in ("baseline", "current"):
        root = tmp_path / phase
        root.mkdir(parents=True)
        for name in ("pyproject.toml", "uv.lock", ".python-version", "justfile"):
            shutil.copyfile(ROOT / name, root / name)
        for pattern in declaration(ROOT).harness:
            for source in ROOT.glob(pattern):
                destination = root / source.relative_to(ROOT)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
        (root / "src").mkdir()
        (root / "src/lib.rs").write_text(f"// {phase}\n", encoding="utf-8", newline="\n")
        (root / ".gitignore").write_text("/target/\n", encoding="utf-8", newline="\n")
        if phase == "current":
            (root / "src/untracked.rs").write_bytes(b"// included working-tree input\n")
        manifest = root / "Cargo.toml"
        contents = manifest.read_text(encoding="utf-8")
        version = tomllib.loads(contents)["package"]["version"]
        manifest.write_text(contents.replace(f'version = "{version}"', 'version = "0.5.1"', 1), encoding="utf-8", newline="\n")
        (root / "tooling/fixture.py").write_text(f"IDS = {ids!r}\nTIMINGS = {synthetic_timings()!r}\n" + DRIVER, encoding="utf-8", newline="\n")
        driver = f"{shlex.quote(sys.executable)} -X utf8 tooling/fixture.py"
        native = 'export CRITERION_HOME := absolute_path(env("CRITERION_HOME", "target/criterion"))\n'
        for side in ("baseline", "current"):
            native += f'\n{side}-all mode="measure": ({side}-vs_linalg mode) ({side}-exact mode)\n'
            for suite in ("vs_linalg", "exact"):
                native += f'\n{side}-{suite} mode="measure":\n    {driver} {side} {suite} {{{{mode}}}}\n'
        native += f"\ngate:\n    {driver} current all gate\n"
        (root / "tooling/performance.just").write_text(native, encoding="utf-8", newline="\n")
        for config in (root / "tooling/measurements").glob("*.toml"):
            config.write_text(
                config.read_text(encoding="utf-8").replace('"tooling/measurements/*.toml"', '"tooling/measurements/*.toml", "tooling/fixture.py"'),
                encoding="utf-8",
                newline="\n",
            )
        run_command("git", ["init", "--quiet"], cwd=root)
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
    return tmp_path / "baseline", tmp_path / "current"


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
    evidence = measure_prepared_pair(baseline, current, declaration(current), PAIR, working_tree=True)
    for path in (ROOT / "tooling").glob("performance-*.toml"):
        shutil.copyfile(path, current / "tooling" / path.name)
    shutil.copytree(ROOT / "tooling/plots", current / "tooling/plots")
    (current / "README.md").write_text("# Fixture\n\n<!-- BENCHMARKS:BEGIN -->\nOld table\n<!-- BENCHMARKS:END -->\n", encoding="utf-8", newline="\n")
    publish_evidence(evidence, current / "target/bench-reports/full/run.json", current / "target/bench-reports/full/evidence.json")
    return current, evidence


@pytest.fixture
def complete_run(tmp_path: Path, measured_run: tuple[Path, Evidence]) -> tuple[Path, Evidence]:
    root, evidence = measured_run
    destination = tmp_path / "current"
    shutil.copytree(root, destination)
    return destination, evidence


def identities(root: Path) -> list[str]:
    """Read assertions captured from the fixture's Git sources before acquisition."""
    return cast("list[str]", json.loads((root / "fixture-identities.json").read_text(encoding="utf-8")))


def test_working_tree_inputs_include_nonignored_untracked_consumer_sources(complete_run: tuple[Path, Evidence]) -> None:
    root, evidence = complete_run
    snapshot = capture_snapshot(root)
    files = {name: data for name, data, _mode in snapshot.untracked}
    assert files["src/untracked.rs"] == b"// included working-tree input\n"
    assert not any(name.startswith("target/") for name in files)
    context = dict(dict(evidence.sources)["current"].context)
    assert "src/untracked.rs" in json.loads(context["source-inventory"])


def test_actual_non_exact_measurement_recipe_discovers_and_retains_complete_samples(tmp_path: Path) -> None:
    _baseline, root = prepared(tmp_path)
    shutil.copytree(ROOT / "tooling/plots/non-exact", root / "tooling/plots/non-exact")
    manifest = root / "Cargo.toml"
    current = manifest.read_text(encoding="utf-8")
    manifest.write_text(current.replace('version = "0.5.1"', 'version = "0.5.0"', 1), encoding="utf-8", newline="\n")
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
    run = validate_run_policy(evidence, declaration(root, "non-exact"))
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
    csv = (root / f"docs/assets/bench/vs_linalg_lu_solve_{statistic}_{axis}.csv").read_text(encoding="utf-8")
    assert "nalgebra" in csv
    assert "faer" in csv
    assert "baseline" in csv
    assert statistic in csv
    result = recipe(root, "performance-readme", *identities(root))
    assert result.returncode == 0, result.stdout + result.stderr
    document = (root / "README.md").read_text(encoding="utf-8")
    assert "Speedup vs nalgebra" in document
    assert "Speedup vs faer" in document
    assert "Complete run evidence" in document
    assert "Evidence envelope" in document
    assert "la-stack v0.5.1" in document
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
            for phase, benchmark in zip(("current", "baseline", "baseline"), catalog["metrics"]["lu_solve"]["benchmarks"], strict=True)
        ]
        assert [float(cell.split("<br>")[0]) for cell in row[1:4]] == pytest.approx(expected)
        # Published ratios have six significant digits; the oracle uses the
        # intended peer/current medians, independently of the view declaration.
        assert [float(cell.removesuffix("\N{MULTIPLICATION SIGN}")) for cell in row[4:]] == pytest.approx(
            [expected[1] / expected[0], expected[2] / expected[0]], rel=1e-5
        )
    checked = recipe(root, "markdown-check")
    assert checked.returncode == 0, checked.stderr


@pytest.mark.parametrize("failure", ["partial-scratch", "stale-source", "stale-harness", "formatter"])
def test_publication_boundaries_leave_all_outputs_unchanged(complete_run: tuple[Path, Evidence], failure: str) -> None:
    root, _evidence = complete_run
    result = recipe(root, "performance-promote")
    assert result.returncode == 0, result.stdout + result.stderr
    if failure == "partial-scratch":
        (root / "target/bench-reports/full/evidence.json").unlink()
    elif failure == "stale-source":
        (root / "src/lib.rs").write_text("// changed\n", encoding="utf-8", newline="\n")
    elif failure == "stale-harness":
        (root / "tooling/performance.just").write_text("# changed\n", encoding="utf-8", newline="\n")
    else:
        config = root / "tooling/performance-readme.toml"
        config.write_text(
            config.read_text(encoding="utf-8").replace('formatter = "pyproject.toml"', 'formatter = "bad-formatter.toml"'), encoding="utf-8", newline="\n"
        )
        (root / "bad-formatter.toml").write_text("[invalid TOML\n", encoding="utf-8", newline="\n")
    outputs = [root / "README.md", *(root / "docs").rglob("*")]
    before = {path: path.read_bytes() for path in outputs if path.is_file()}
    result = recipe(root, "performance-promote") if failure == "partial-scratch" else recipe(root, "performance-readme", *identities(root))
    assert result.returncode != 0
    if failure == "stale-source":
        assert "current-sources differ" in result.stderr
    elif failure == "stale-harness":
        assert "current-harness differ" in result.stderr
    elif failure == "formatter":
        assert "rumdl could not format README.md" in result.stderr
    if failure == "partial-scratch":
        result = recipe(root, "plot-vs-linalg", *identities(root), "lu_solve-median-log", "scratch")
        assert result.returncode != 0
    assert {path: path.read_bytes() for path in before} == before
    assert not (root / "docs/assets").exists()


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
    result = recipe(root, "performance-measure", "v0.4.6" if phase == "current" else PAIR.current, "v0.4.6" if phase == "baseline" else PAIR.baseline)
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
    with pytest.raises(ValueError, match="acquisition declaration"):
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
