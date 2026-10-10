"""Installed-package acceptance of la-stack's complete performance workflow."""

import hashlib
import io
import json
import os
import shlex
import shutil
import sys
import tarfile
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from research_repo_tools.common_measurement import CommonHarnessPlan, measure_prepared_pair
from research_repo_tools.complete_runs import RunSeries, render_run, run_identity, serialize_run
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_command
from research_repo_tools.release_pairs import ReleasePair

import archive_performance as workflow
import criterion_dim_plot as plot
import performance_runs as policy
from bench_compare import render_release_artifacts
from performance_artifacts import ArtifactPaths
from release_baseline import required_report_ids

if TYPE_CHECKING:
    from research_repo_tools.evidence import Evidence

ROOT = Path(__file__).resolve().parents[2]
PAIR = ReleasePair("v0.4.6", "v0.4.5")
DRIVER = """import json, os, sys
from pathlib import Path
action, phase = sys.argv[1:]
fault = os.environ.get("FIXTURE_FAULT", "")
if action == "gate":
    if fault == "gate":
        raise SystemExit(1)
    if fault == "stale":
        Path("src/lib.rs").write_text("changed during gate")
else:
    ids = INVENTORIES[phase]
    for name in ids[:-1] if fault == "missing-row" else ids:
        path = Path("target/criterion") / name / "new"
        path.mkdir(parents=True)
        value = 10.0 if phase == "baseline" else 8.0
        estimate = {"point_estimate": value, "confidence_interval": {"lower_bound": value * .9, "upper_bound": value * 1.1, "confidence_level": .95}}
        estimates = {"median": estimate} if fault == "missing-statistic" else {"mean": estimate, "median": estimate}
        count = 99 if fault == "missing-sample" else 100
        sample = {"iters": [1] * count, "times": [value] * count}
        for filename, data in (("benchmark.json", {"full_id": name}), ("estimates.json", estimates), ("sample.json", sample)):
            (path / filename).write_text(json.dumps(data))
"""


def prepared(tmp_path: Path, fault: str = "") -> tuple[Path, Path, CommonHarnessPlan]:
    """Use the real consumer inventory and policies with a bounded timing fixture."""
    ids = tuple(sorted(required_report_ids()))
    inventories = {
        "baseline": tuple(name for name in ids if not name.startswith(policy.UNAVAILABLE_PREFIXES)),
        "current": tuple(name for name in ids if not name.startswith("d") or "/la_stack_" in name),
    }
    for phase in ("baseline", "current"):
        root = tmp_path / phase
        root.mkdir(parents=True)
        (root / "src").mkdir()
        (root / "src/lib.rs").write_text(f"// {phase}\n", encoding="utf-8")
        for pattern in policy.configuration(ROOT).harness:
            for source in ROOT.glob(pattern):
                destination = root / source.relative_to(ROOT)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
        (root / "tooling/fixture.py").write_text(f"INVENTORIES = {inventories!r}\n" + DRIVER, encoding="utf-8")
        (root / "tooling/performance.just").write_text(
            "\n".join(
                f"{recipe}:\n    {shlex.quote(sys.executable)} tooling/fixture.py {action} {phase}\n"
                for recipe, action, phase in (("gate", "gate", "current"), ("baseline-all", "measure", "baseline"), ("current-all", "measure", "current"))
            ),
            encoding="utf-8",
        )
        config = (root / policy.MEASUREMENT_CONFIG).read_text(encoding="utf-8")
        (root / policy.MEASUREMENT_CONFIG).write_text(
            config.replace('    "tooling/performance.just",', '    "tooling/performance.just", "tooling/fixture.py",'), encoding="utf-8"
        )
        run_command("git", ["--no-pager", "init", "--quiet"], cwd=root)
        run_command(
            "git",
            [
                "--no-pager",
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
    plan = policy.declared_plan(
        PAIR, "all", "release-signal", inventories, baseline_env=(("FIXTURE_FAULT", fault),), current_env=(), config=policy.configuration(current)
    )
    return baseline, current, plan


def measured(tmp_path: Path) -> tuple[Path, Evidence]:
    baseline, current, plan = prepared(tmp_path)
    evidence = measure_prepared_pair(baseline, current, plan, PAIR, working_tree=True)
    (current / "tooling").mkdir(exist_ok=True)
    shutil.copyfile(ROOT / policy.REPORT_CONFIG, current / policy.REPORT_CONFIG)
    return current, evidence


@pytest.fixture(scope="module")
def measured_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Evidence]:
    """Measure once; publication regressions operate on independent copies."""
    return measured(tmp_path_factory.mktemp("measured-run"))


@pytest.fixture
def complete_run(tmp_path: Path, measured_run: tuple[Path, Evidence]) -> tuple[Path, Evidence]:
    root, evidence = measured_run
    destination = tmp_path / "current"
    shutil.copytree(root, destination)
    return destination, evidence


@pytest.fixture
def legacy_artifacts(tmp_path: Path) -> ArtifactPaths:
    """Keep historical promotion coverage on the last published CSV schema."""
    history = ROOT / "docs/performance"
    run = json.loads((history / "latest.json").read_bytes())["run"]
    for name in ("performance.csv", "performance.provenance.json"):
        shutil.copyfile(history / run / name, tmp_path / name)
    return ArtifactPaths(tmp_path / "performance.csv", tmp_path / "performance.provenance.json")


def test_common_harness_and_reference_phase(tmp_path: Path) -> None:
    baseline, current, plan = prepared(tmp_path)
    (baseline / "benches/exact.rs").write_text("obsolete harness\n", encoding="utf-8")
    evidence = measure_prepared_pair(baseline, current, plan, PAIR, working_tree=True)
    run = policy.validate_scientific_run(evidence)
    sources = dict(evidence.sources)
    assert sources["baseline"].harness_sha256 == sources["current"].harness_sha256
    assert sources["baseline"].source_sha256 != sources["current"].source_sha256
    assert (baseline / "benches/exact.rs").read_bytes() == (current / "benches/exact.rs").read_bytes()
    assert {series.phase for series in run.series if series.name in {"nalgebra", "faer"}} == {"baseline"}
    assert all(sample.policy.statistics == ("mean", "median") for _, sample in run.phases)


@pytest.mark.parametrize("fault", ["gate", "stale", "missing-row", "missing-statistic", "missing-sample"])
def test_declared_gates_and_completeness_fail_closed(tmp_path: Path, fault: str) -> None:
    baseline, current, plan = prepared(tmp_path, fault)
    with pytest.raises(ValueError, match=r"gate|changed|incomplete|sample|statistic|expected|missing|mean"):
        measure_prepared_pair(baseline, current, plan, PAIR, working_tree=True)
    assert not (current / "target/criterion").exists()
    if fault in {"gate", "stale"}:
        assert not (baseline / "target/criterion").exists()


@pytest.mark.parametrize("report_path", ["docs/performance.md", "docs/reports/current.md"])
def test_shared_retention_repeated_runs_and_offline_report(complete_run: tuple[Path, Evidence], report_path: str) -> None:
    root, first = complete_run
    config = root / policy.REPORT_CONFIG
    config.write_text(
        config.read_text(encoding="utf-8")
        .replace("docs/performance-runs", "docs/performance-history")
        .replace("docs/performance.md", report_path)
        .replace("la-stack complete benchmark measurements", "Configured performance report"),
        encoding="utf-8",
        newline="\n",
    )
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, first)
    assert stem.with_suffix(".md").read_bytes() == render_run(first)
    workflow.promote(root, stem)
    first_report = (root / report_path).read_bytes()
    assert first_report == render_run(first, title="Configured performance report")
    assert b"## Mean (ns)" in first_report
    assert b"## Median (ns)" in first_report
    assert b"nalgebra: measured in **baseline** phase" in first_report
    assert b"faer: measured in **baseline** phase" in first_report
    assert not (root / policy.archive_directory(root) / "current.md").exists()
    if report_path != "docs/performance.md":
        assert not (root / "docs/performance.md").exists()
    retained = root / policy.archive_directory(root) / "runs" / run_identity(first)
    original = {path: path.read_bytes() for path in retained.iterdir()}
    # A second valid measurement with the same release pair has new provenance.
    sources = tuple((phase, replace(source, context=(*source.context, ("experiment", "second")))) for phase, source in first.sources)
    second = replace(first, sources=sources)
    workflow.save_scratch(root, stem, second)
    workflow.promote(root, stem)
    assert len(list((root / policy.archive_directory(root) / "runs").iterdir())) == 2
    assert {path: path.read_bytes() for path in original} == original
    assert (retained / "run.json").is_file()
    assert first_report != (root / report_path).read_bytes()
    assert not (root / "docs/archive/performance").exists()
    shutil.rmtree(root / "target")
    assert run_identity(policy.retained_run(root, stem)) == run_identity(second)
    before = (root / report_path).read_bytes()
    workflow.promote(root, stem)
    assert (root / report_path).read_bytes() == before
    (root / policy.archive_directory(root) / "latest.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match=r"latest|pointer|fields|schema"):
        policy.retained_run(root, stem)


@pytest.mark.parametrize("destination", ["missing", "identical", "different", "directory", "symlink", "relative-symlink"])
def test_legacy_promotion_preserves_prior_report_on_archive_collision(tmp_path: Path, legacy_artifacts: ArtifactPaths, destination: str) -> None:
    root = tmp_path
    stem = root / "target/bench-reports/performance"
    current = root / "docs/performance.md"
    current.parent.mkdir(parents=True, exist_ok=True)
    previous = b"**la-stack** v0.4.5\n\nComparison against baseline **v0.4.4**:\n\nPrior curated report.\n"
    current.write_bytes(previous)
    archived = root / "docs/archive/performance/v0.4.5-vs-v0.4.4.md"
    archived.parent.mkdir(parents=True)
    if destination in {"identical", "different"}:
        archived.write_bytes(previous if destination == "identical" else b"Conflicting archive.\n")
    elif destination == "directory":
        archived.mkdir()
    elif destination == "symlink":
        archived.symlink_to(current)
    elif destination == "relative-symlink":
        archived.symlink_to(os.path.relpath(current, archived.parent))
    original_link_target = archived.readlink() if archived.is_symlink() else None
    before = {path: path.read_bytes() for path in root.rglob("*") if path.is_file() and ".git" not in path.parts}
    if destination in {"missing", "identical"}:
        workflow.render_and_promote_artifacts(artifacts=legacy_artifacts, output=stem.with_suffix(".md"), current=current, archive_dir=archived.parent)
        assert archived.read_bytes() == previous
        assert current.read_text(encoding="utf-8") == render_release_artifacts(legacy_artifacts)
    else:
        errors = {"different": "immutable output differs", "directory": "regular file", "symlink": "symlink", "relative-symlink": "symlink"}
        with pytest.raises(ValueError, match=errors[destination]):
            workflow.render_and_promote_artifacts(artifacts=legacy_artifacts, output=stem.with_suffix(".md"), current=current, archive_dir=archived.parent)
        assert {path: path.read_bytes() for path in root.rglob("*") if path.is_file() and ".git" not in path.parts} == before
        assert not (root / "docs/performance-runs").exists()
        if destination == "directory":
            assert archived.is_dir()
        elif destination in {"symlink", "relative-symlink"}:
            assert archived.is_symlink()
            assert archived.readlink() == original_link_target
            assert archived.samefile(current)


@pytest.mark.parametrize("duplicate", ["**la-stack** v0.4.4", "Comparison against baseline **v0.4.3**:"])
def test_duplicate_legacy_report_identity_prevents_publication(tmp_path: Path, legacy_artifacts: ArtifactPaths, duplicate: str) -> None:
    root = tmp_path
    stem = root / "target/bench-reports/performance"
    current = root / "docs/performance.md"
    current.parent.mkdir(parents=True, exist_ok=True)
    text = f"**la-stack** v0.4.5\n\nComparison against baseline **v0.4.4**:\n\n{duplicate}\n"
    current.write_text(text, encoding="utf-8", newline="\n")
    before = {path: path.read_bytes() for path in root.rglob("*") if path.is_file() and ".git" not in path.parts}
    with pytest.raises(ValueError, match="exactly one"):
        workflow.parse_report_id(text)
    with pytest.raises(ValueError, match="exactly one"):
        workflow.render_and_promote_artifacts(
            artifacts=legacy_artifacts, output=stem.with_suffix(".md"), current=current, archive_dir=root / "docs/archive/performance"
        )
    assert {path: path.read_bytes() for path in root.rglob("*") if path.is_file() and ".git" not in path.parts} == before
    assert not (root / "docs/performance-runs").exists()


def test_partial_scratch_and_stale_source_cannot_select_old_history(tmp_path: Path) -> None:
    root, evidence = measured(tmp_path)
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    policy.scratch_paths(stem)[1].unlink()
    with pytest.raises(FileNotFoundError):
        policy.retained_run(root, stem)
    (root / "src/lib.rs").write_text("// stale\n", encoding="utf-8")
    with pytest.raises(ValueError, match="stale"):
        policy.verify_current_source(root, evidence)


@pytest.mark.parametrize("archive", ["docs/performance-runs", "docs/performance-history"])
@pytest.mark.parametrize("csv_path", ["relative", "dot-relative", "absolute"])
def test_readme_uses_retained_reference_estimates_after_cleanup(
    complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch, archive: str, csv_path: str
) -> None:
    root, evidence = complete_run
    config = root / policy.REPORT_CONFIG
    config.write_text(config.read_text(encoding="utf-8").replace("docs/performance-runs", archive), encoding="utf-8")
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    shutil.rmtree(root / "target")
    monkeypatch.setattr(plot, "_repo_root", lambda: root)
    monkeypatch.setattr(plot, "_render_svg_with_gnuplot", lambda request: request.out_svg.write_bytes(b"<svg/>\n"))
    performance_csv = "target/bench-reports/performance.csv"
    if csv_path == "absolute":
        performance_csv = str(root / performance_csv)
    elif csv_path == "dot-relative":
        performance_csv = f"./{performance_csv}"
    arguments = ["--update-readme", "--log-y", "--performance-csv", performance_csv]
    assert plot.main(arguments) == 0
    assets = root / "docs/assets/bench"
    csv = assets / "vs_linalg_lu_solve_median.csv"
    assert "2,8.0,7.2,8.8,10.0,9.0,11.0,10.0,9.0,11.0" in csv.read_text(encoding="utf-8")
    provenance = json.loads(csv.with_suffix(".provenance.json").read_bytes())
    assert provenance["measurement"]["peer_phase"] == "baseline"
    assert provenance["run_id"] == run_identity(evidence)
    before = {path: path.read_bytes() for path in (*assets.iterdir(), root / "README.md")}
    assert plot.main(arguments) == 0
    assert {path: path.read_bytes() for path in before} == before


@pytest.mark.parametrize(
    ("cleanup", "destination"),
    [
        (False, "run"),
        (False, "evidence"),
        (False, "alias"),
        (False, "source"),
        (False, "harness"),
        (False, "configuration"),
        (True, "run"),
        (True, "evidence"),
        (True, "report"),
        (True, "index"),
        (True, "latest"),
    ],
)
def test_readme_outputs_cannot_replace_complete_run_inputs(
    complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], cleanup: bool, destination: str
) -> None:
    root, evidence = complete_run
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    archive = root / policy.archive_directory(root)
    retained = archive / "runs" / run_identity(evidence)
    payload, manifest = policy.scratch_paths(stem)
    if cleanup:
        shutil.rmtree(root / "target")
        payload, manifest = retained / "run.json", retained / "evidence.json"
    alias = root / "docs/alias.csv"
    if destination == "alias":
        alias.hardlink_to(payload)
    targets = {
        "run": payload,
        "evidence": manifest,
        "alias": alias,
        "source": root / "src/lib.rs",
        "harness": root / "tooling/performance.just",
        "configuration": root / policy.REPORT_CONFIG,
        "report": retained / "report.md",
        "index": archive / "index.json",
        "latest": archive / "latest.json",
    }
    shutil.copyfile(ROOT / "README.md", root / "docs/custom.md")
    before = {path: path.read_bytes() for path in root.rglob("*") if path.is_file() and ".git" not in path.parts}
    monkeypatch.setattr(plot, "_repo_root", lambda: root)
    monkeypatch.setattr(plot, "_render_svg_with_gnuplot", lambda request: request.out_svg.write_bytes(b"<svg/>\n"))
    assert plot.main(["--update-readme", "--readme", "docs/custom.md", "--csv", str(targets[destination])]) == 2
    assert "must use distinct paths" in capsys.readouterr().err
    assert {path: path.read_bytes() for path in root.rglob("*") if path.is_file() and ".git" not in path.parts} == before


def test_complete_report_passes_markdown_recipes_without_rewriting_history(complete_run: tuple[Path, Evidence]) -> None:
    root, evidence = complete_run
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    # Exercise the real recipes against generated, nonignored Markdown in a
    # disposable checkout, reusing the already locked test environment.
    for name in ("pyproject.toml", "uv.lock", ".python-version"):
        shutil.copyfile(ROOT / name, root / name)
    history = root / policy.archive_directory(root)
    originals = {path: path.read_bytes() for path in history.rglob("*") if path.is_file()}
    originals[root / "docs/performance.md"] = (root / "docs/performance.md").read_bytes()
    for recipe in ("markdown-check", "markdown-fix", "markdown-check"):
        result = run_command(
            "just",
            ["--no-deps", "--justfile", str(ROOT / "justfile"), "--working-directory", str(root), recipe],
            env={**os.environ, "UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
    assert {path: path.read_bytes() for path in originals} == originals
    shutil.rmtree(root / "target")
    assert run_identity(policy.retained_run(root, stem)) == run_identity(evidence)


def test_legacy_bytes_hashes_and_report_are_unchanged() -> None:
    directory = ROOT / "docs/performance/v0.4.6-vs-v0.4.5/39719a1eb599d294dbbc1e46dc6fb0de2028ad2e59aab5f94fb48c7bf3d7f841"
    paths = ArtifactPaths(directory / "performance.csv", directory / "performance.provenance.json")
    hashes = {
        "performance.csv": "b263bb34bcfd43b8fc8109110db9f2e908516aef6c01ce16f70f5ac18dc7ef2d",
        "performance.full.csv": "4fa5cc13159f5dc91e711ae67554845bcbd1ace2cfc3e31487fc85f75d3dc459",
        "performance.full.provenance.json": "a43f66e52035aecd8d235dfded916b23c3240529d0c07f793d27322521cfb0d0",
        "performance.provenance.json": "89144a2899ba73f55d009da4fa20b4e6477600902b93cde469383f719c497508",
    }
    assert {name: hashlib.sha256((paths.csv.parent / name).read_bytes()).hexdigest() for name in hashes} == hashes
    assert hashlib.sha256(render_release_artifacts(paths).encode()).hexdigest() == "4ea0ff9bca4c7bb92cd57bd9163b3cacabe78f502abdd2663162cf9562ec6582"


def test_native_command_policy_and_compatibility_adapter() -> None:
    justfile = Path("tooling/performance.just")
    baseline = dry_run(ROOT, "baseline-all", justfile=justfile).stderr
    current = dry_run(ROOT, "current-all", justfile=justfile).stderr
    assert baseline.count("cargo bench --locked") == current.count("cargo bench --locked") == 2
    assert baseline.count("--noplot") == current.count("--noplot") == 2
    assert "la_stack" in current
    assert "la_stack" not in baseline
    assert "--list" not in baseline
    assert dry_run(ROOT, "baseline-all", ["list"], justfile=justfile).stderr.count("--list") == 2
    gate = dry_run(ROOT, "gate", justfile=justfile).stderr
    assert "--test vs_linalg_inputs --test exact_bench_config" in gate
    environment = dict(policy.phase_environment(ROOT, "la_stack_pre_rational_input_api"))
    assert "--cfg=la_stack_pre_rational_input_api" in environment.get("CARGO_ENCODED_RUSTFLAGS", environment.get("RUSTFLAGS", ""))
    assert environment["CARGO_TARGET_DIR"] == "target"
    assert environment["CRITERION_HOME"] == "target/criterion"


@pytest.mark.parametrize("baseline_tag", ["v0.4.4", "v0.4.5", "v0.4.6"])
@pytest.mark.parametrize("scope", ["release-signal", "all-benches"])
def test_native_inventory_and_phase_commands_use_the_measured_checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, baseline_tag: str, scope: str) -> None:
    """Exercise the configured Justfile, stubbing only the expensive Cargo boundary."""
    root = tmp_path / "checkout with spaces"
    shutil.copytree(ROOT / "tooling", root / "tooling")
    shutil.copyfile(ROOT / "rust-toolchain.toml", root / "rust-toolchain.toml")
    binary = tmp_path / "bin"
    binary.mkdir()
    cargo = binary / "cargo"
    # These are the registrations guarded by the historical API cfg in
    # benches/vs_linalg.rs, including the two reference scenario kernels.
    norm_rows = {
        f"d{dimension}/{bench}"
        for dimension in (2, 3, 4, 5)
        for bench in (
            "la_stack_norm2",
            *(
                f"{kernel}_norm2_scenario_{scenario}"
                for kernel in ("la_stack", "iterative_f64_hypot", "delaunay_scaled_norm")
                for scenario in ("descending", "repeated_scale", "sparse", "wide_dynamic_range")
            ),
        )
    }
    ordinary_norm_peers = {
        f"d{dimension}/{bench}" for dimension in (2, 3, 4, 5) for bench in ("iterative_f64_hypot", "delaunay_scaled_norm", "nalgebra_norm", "faer_norm_l2")
    }
    ids = required_report_ids() | norm_rows | ordinary_norm_peers
    exact = sorted(name for name in ids if not name.startswith("d"))
    baseline_ids = sorted(name for name in ids if name.startswith("d"))
    current_ids = sorted(name for name in baseline_ids if "/la_stack_" in name)
    cargo.write_text(
        "#!/bin/sh\n"
        "directory=$(pwd -W 2>/dev/null || pwd -P)\n"
        'printf "%s|%s\\n" "$directory" "$*" >> "$CARGO_RECORD"\n'
        'case "$*" in\n'
        f'  *"--bench exact"*) printf "%s: benchmark\\n" {shlex.join(exact)} ;;\n'
        f'  *"la_stack --noplot"*) printf "%s: benchmark\\n" {shlex.join(current_ids)} ;;\n'
        f'  *"--bench vs_linalg"*) printf "%s: benchmark\\n" {shlex.join(baseline_ids)} ;;\n'
        "esac\n",
        encoding="utf-8",
        newline="\n",
    )
    cargo.chmod(0o755)
    record = tmp_path / "commands.txt"
    monkeypatch.setenv("PATH", f"{binary}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("CARGO_RECORD", str(record))
    pair = ReleasePair("v0.4.6", baseline_tag)
    plan = policy.measurement_plan(root, pair, "all", scope, tmp_path / "inventory")
    assert set(plan.current.policy.expected) == set(current_ids) | set(exact)
    unavailable = norm_rows | {name for name in exact if name.startswith(("rational_input_", "canonical_conversion_", "det4_diagnostic_"))}
    assert set(plan.baseline.policy.expected) == (ids if baseline_tag == "v0.4.6" else ids - unavailable)
    assert ordinary_norm_peers <= set(plan.baseline.policy.expected)
    assert "d2/la_stack_norm2_sq" in plan.baseline.policy.expected
    series = {item.name: {name for name, _ in item.rows} for item in plan.series}
    assert series["la-stack baseline"] == (series["la-stack current"] if baseline_tag == "v0.4.6" else series["la-stack current"] - unavailable)
    if scope == "all-benches":
        assert "d2/la_stack_norm2" in series["la-stack current"]
        assert "d2/la_stack_norm2_scenario_sparse" in series["la-stack current"]
    for command in (plan.baseline.gate, plan.baseline.command, plan.current.command):
        run_command(command[0], command[1:], cwd=root)
    lines = record.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 9
    for line in lines:
        directory, separator, arguments = line.partition("|")
        assert separator == "|"
        assert Path(directory).samefile(root)
        assert "--locked" in arguments.split()
    assert all("--list" in line for line in lines[:4])
    assert all("--list" not in line for line in lines[4:])
    assert "test --locked --workspace --features bench,exact --test vs_linalg_inputs --test exact_bench_config" in lines[4]


def test_pre_just_run_remains_readable_without_its_python_driver(tmp_path: Path) -> None:
    _, evidence = measured(tmp_path)
    run = policy.validate_scientific_run(evidence)
    cargo = {
        "baseline": [
            ["cargo", "bench", "--locked", "-p", "la-stack-comparison", "--features", "bench", "--bench", "vs_linalg", "--", "--noplot"],
            ["cargo", "bench", "--locked", "--features", "bench,exact", "--bench", "exact", "--", "--noplot"],
        ],
        "current": [
            ["cargo", "bench", "--locked", "-p", "la-stack-comparison", "--features", "bench", "--bench", "vs_linalg", "--", "la_stack", "--noplot"],
            ["cargo", "bench", "--locked", "--features", "bench,exact", "--bench", "exact", "--", "--noplot"],
        ],
    }
    sources = []
    for phase, source in evidence.sources:
        context = dict(source.context)
        del context["phase-driver"]
        del context["tool.just"]
        context.update((f"{name}-cargo-commands", json.dumps(commands)) for name, commands in cargo.items())
        for key, action in (("command", "measure"), ("gate-command", "gate")):
            context[key] = json.dumps(["/retired/venv/python", "scripts/performance_phase.py", action, "all", phase])
        sources.append((phase, replace(source, context=tuple(context.items()))))
    legacy = replace(
        evidence,
        sources=tuple(sources),
        payload=serialize_run(replace(run, compatible=tuple(field for field in run.compatible if field != "context.tool.just"))),
    )
    policy.validate_scientific_run(legacy)
    bad = replace(sources[0][1], context=tuple((key, "[]" if key == "baseline-cargo-commands" else value) for key, value in sources[0][1].context))
    with pytest.raises(ValueError, match="consumer policy"):
        policy.validate_scientific_run(replace(legacy, sources=(("baseline", bad), sources[1])))


def test_local_explicit_labels_and_release_eligibility() -> None:
    assert workflow.select_pair(ROOT, "explicit", "v0.4.6", "v0.4.6", local=True) == ReleasePair("v0.4.6", "v0.4.6")
    with pytest.raises(ValueError, match="must differ"):
        workflow.select_pair(ROOT, "explicit", "v0.4.6", "v0.4.6")
    with pytest.raises(ValueError, match=r"v0.4.4|unsupported"):
        workflow.select_pair(ROOT, "explicit", "v0.4.6", "v0.4.3", local=True)


def test_current_gate_failure_prevents_both_timing_phases(tmp_path: Path) -> None:
    baseline, current, plan = prepared(tmp_path)
    plan = replace(plan, current=replace(plan.current, environment=(("FIXTURE_FAULT", "gate"),)))
    with pytest.raises(ValueError, match="current preflight gate failed"):
        measure_prepared_pair(baseline, current, plan, PAIR, working_tree=True)
    assert not (baseline / "target/criterion").exists()
    assert not (current / "target/criterion").exists()


def test_reference_phase_and_same_version_policy(tmp_path: Path) -> None:
    _, evidence = measured(tmp_path)
    run = policy.validate_scientific_run(evidence)
    # Same semantic IDs exist in both phases, so only consumer policy detects
    # an otherwise well-formed la-stack series assigned to the wrong phase.
    series = tuple(RunSeries(item.name, "current", item.rows) if item.name == "la-stack baseline" else item for item in run.series)
    with pytest.raises(ValueError, match="retained series"):
        policy.validate_scientific_run(replace(evidence, payload=serialize_run(replace(run, series=series))))
    # Use a current/current release pair with both sources under the same API.
    baseline, current, plan = prepared(tmp_path / "same")
    pair = ReleasePair("v0.4.6", "v0.4.6")
    inventories = {phase: tuple(case.full_id for case in sample.cases) for phase, sample in run.phases}
    # The older-API fixture deliberately lacks rows; removing them from the
    # current selection allows a narrow all-benches same-version comparison.
    inventories["current"] = tuple(name for name in inventories["current"] if name in inventories["baseline"])
    driver = f"INVENTORIES = {inventories!r}\n" + DRIVER
    (current / "tooling/fixture.py").write_text(driver, encoding="utf-8")
    plan = policy.declared_plan(pair, "all", "all-benches", inventories, baseline_env=(), current_env=(), config=policy.configuration(current))
    same = measure_prepared_pair(baseline, current, plan, pair, working_tree=True)
    policy.validate_scientific_run(same)
    with pytest.raises(ValueError, match="same-version"):
        policy.validate_scientific_run(same, release=True)
    sources = dict(same.sources)
    identical = replace(
        same, sources=(("baseline", replace(sources["baseline"], source_sha256=sources["current"].source_sha256)), ("current", sources["current"]))
    )
    with pytest.raises(ValueError, match="distinct source"):
        policy.validate_scientific_run(identical)


def test_publication_io_failure_preserves_all_published_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root, evidence = measured(tmp_path)
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    originals = {path: path.read_bytes() for path in (root / "docs").rglob("*") if path.is_file()}
    second = replace(
        evidence, sources=tuple((phase, replace(source, context=(*source.context, ("experiment", "second")))) for phase, source in evidence.sources)
    )
    workflow.save_scratch(root, stem, second)
    native_replace = Path.replace
    failed = False

    def fail_once(source: Path, destination: Path) -> Path:
        nonlocal failed
        if destination == root / policy.archive_directory(root) / "latest.json" and not failed:
            failed = True
            msg = "injected publication failure"
            raise OSError(msg)
        return native_replace(source, destination)

    monkeypatch.setattr(Path, "replace", fail_once)
    with pytest.raises(OSError, match="injected publication failure"):
        workflow.promote(root, stem)
    assert failed
    assert {path: path.read_bytes() for path in (root / "docs").rglob("*") if path.is_file()} == originals


def test_unsafe_release_asset_cannot_replace_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output = tmp_path / "report.md"
    output.write_bytes(b"previous report\n")

    def unsafe_asset(_root: Path, _repo: str, _tag: str, _name: str, destination: Path) -> None:
        with tarfile.open(destination, "w:gz") as archive:
            member = tarfile.TarInfo("../../report.md")
            member.size = 4
            archive.addfile(member, io.BytesIO(b"evil"))

    monkeypatch.setattr(workflow, "download_release_asset", unsafe_asset)
    with pytest.raises(ValueError, match=r"unsafe|traversal|relative|member"):
        workflow.compare_assets(tmp_path, PAIR, output, "all", "release-signal")
    assert output.read_bytes() == b"previous report\n"
