"""Installed-package acceptance of la-stack's complete performance workflow."""

import hashlib
import io
import json
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from research_repo_tools import publication
from research_repo_tools.common_measurement import CommonHarnessPlan, measure_prepared_pair
from research_repo_tools.complete_runs import CompletePolicy, RunSeries, render_run, run_identity, serialize_run
from research_repo_tools.criterion import parse_comparison
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.legacy_evidence import convert_csv
from research_repo_tools.process import ExecutableNotFoundError, run_command
from research_repo_tools.release_pairs import ReleasePair

import archive_performance as plot
import archive_performance as workflow
import performance_runs as policy
from performance_runs import required_report_ids

if TYPE_CHECKING:
    from research_repo_tools.evidence import Evidence
    from research_repo_tools.publication import PublicationPlan

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
        "baseline": tuple(name for name in ids if not name.startswith(policy.scientific_policy().unavailable_prefixes)),
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
    shutil.copyfile(ROOT / workflow.PLOT_CONFIG, current / workflow.PLOT_CONFIG)
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
    assert stem.with_suffix(".md").read_bytes() == render_run(
        first,
        report_path="target/bench-reports/performance.md",
        payload_path="target/bench-reports/performance.run.json",
        manifest_path="target/bench-reports/performance.evidence.json",
    )
    workflow.promote(root, stem)
    first_report = (root / report_path).read_bytes()
    assert b"Configured performance report" in first_report
    assert b"Recorded measurement provenance" in first_report
    assert b"run.json" in first_report
    assert b"evidence.json" in first_report
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
def test_readme_uses_original_reference_phases_after_cleanup(complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch, archive: str) -> None:
    root, evidence = complete_run
    config = root / policy.REPORT_CONFIG
    config.write_text(config.read_text(encoding="utf-8").replace("docs/performance-runs", archive), encoding="utf-8")
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    monkeypatch.setattr(plot, "_repo_root", lambda: root)
    monkeypatch.chdir(root.parent)
    arguments = ["--update-readme", "--log-y"]
    assert plot.dimension_main(arguments) == 0
    assets = root / "docs/assets/bench"
    csv = assets / "vs_linalg_lu_solve_median.csv"
    import csv as csv_module  # noqa: PLC0415

    with csv.open(encoding="utf-8", newline="") as stream:
        rows = list(csv_module.DictReader(stream))
    assert len(rows) == 24
    assert {row["phase"] for row in rows if row["series"] in {"nalgebra", "faer"}} == {"baseline"}
    assert {row["phase"] for row in rows if row["series"] == "la-stack current"} == {"current"}
    assert {float(row["coordinate"]) for row in rows} == {2, 3, 4, 5, 8, 16, 32, 64}
    assert {row["unit"] for row in rows} == {"ns"}
    dimension_two = [row for row in rows if float(row["coordinate"]) == 2]
    assert [(float(row["point"]), float(row["lower"]), float(row["upper"])) for row in dimension_two] == [(8.0, 7.2, 8.8), (10.0, 9.0, 11.0), (10.0, 9.0, 11.0)]
    provenance = json.loads(csv.with_suffix(".provenance.json").read_bytes())
    assert provenance["series"]["nalgebra"] == "baseline"
    assert provenance["run_id"] == run_identity(evidence)
    assert b"<svg" in csv.with_suffix(".svg").read_bytes()
    assert b"nalgebra" in csv.with_suffix(".svg").read_bytes()
    before = {path: path.read_bytes() for path in (*assets.iterdir(), root / "README.md")}
    run_command("git", ["--no-pager", "add", "docs/assets/bench"], cwd=root)
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
            "--quiet",
            "-m",
            "assets",
        ],
        cwd=root,
    )
    run_command("git", ["--no-pager", "-c", "tag.gpgsign=false", "tag", PAIR.current], cwd=root)
    shutil.rmtree(root / "target")
    assert plot.dimension_main(arguments) == 0
    assert {path: path.read_bytes() for path in before} == before
    assert provenance["payload"] == f"{archive}/runs/{run_identity(evidence)}/run.json"
    assert provenance["manifest"] == f"{archive}/runs/{run_identity(evidence)}/evidence.json"


@pytest.mark.parametrize(
    ("error", "diagnostic"),
    [
        (ExecutableNotFoundError("Git executable not found"), "Git executable not found"),
        (subprocess.CalledProcessError(128, ["git"], stderr=b"fatal: repository unavailable"), "fatal: repository unavailable"),
        (subprocess.TimeoutExpired(["git"], 30), "command timed out after 30 seconds"),
    ],
    ids=["missing-git", "git-failure", "git-timeout"],
)
def test_plot_git_failure_preserves_publication(
    complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], error: Exception, diagnostic: str
) -> None:
    root, evidence = complete_run
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    assets = root / "docs/assets/bench"
    assets.mkdir(parents=True)
    for suffix in (".svg", ".csv", ".provenance.json"):
        (assets / f"vs_linalg_lu_solve_median{suffix}").write_bytes(b"previous artifact\n")
    before = {path: path.read_bytes() for path in (*assets.iterdir(), root / "README.md")}

    def fail_git(*_args: object, **_kwargs: object) -> None:
        raise error

    monkeypatch.setattr(publication, "run_git_bytes", fail_git)
    monkeypatch.setattr(plot, "_repo_root", lambda: root)
    assert plot.dimension_main(["--update-readme"]) == 2
    output = capsys.readouterr()
    assert output.out == ""
    assert output.err.startswith("dimension publication:")
    assert diagnostic in output.err
    assert "Traceback" not in output.err
    assert {path: path.read_bytes() for path in (*assets.iterdir(), root / "README.md")} == before


def test_stale_source_prevents_dimension_publication(complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch) -> None:
    root, evidence = complete_run
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    original = (root / "README.md").read_bytes()
    (root / "src/lib.rs").write_text("// changed library\n", encoding="utf-8")
    monkeypatch.setattr(plot, "_repo_root", lambda: root)
    assert plot.dimension_main(["--update-readme"]) == 2
    assert (root / "README.md").read_bytes() == original
    assert not (root / "docs/assets/bench").exists()


@pytest.mark.parametrize("retained", [False, True])
def test_readme_requires_matching_retained_evidence(
    complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], retained: bool
) -> None:
    root, evidence = complete_run
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    if retained:
        workflow.promote(root, stem)
        sources = tuple((phase, replace(source, context=(*source.context, ("experiment", "unpromoted")))) for phase, source in evidence.sources)
        workflow.save_scratch(root, stem, replace(evidence, sources=sources))
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    before = {path: path.read_bytes() for path in (root / "docs").rglob("*") if path.is_file()} if retained else {}
    before[root / "README.md"] = (root / "README.md").read_bytes()
    monkeypatch.setattr(plot, "_repo_root", lambda: root)

    assert plot.dimension_main(["--update-readme"]) == 2

    output = capsys.readouterr()
    assert output.out == ""
    assert ("README evidence differs from the retained latest run" if retained else "no retained latest run") in output.err
    assert {path: path.read_bytes() for path in before} == before
    assert not (root / "docs/assets/bench").exists()


def test_readme_edit_during_rendering_is_preserved(complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch) -> None:
    root, evidence = complete_run
    workflow.save_scratch(root, root / "target/bench-reports/performance", evidence)
    workflow.promote(root, root / "target/bench-reports/performance")
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    document = root / "README.md"
    edited = document.read_bytes() + b"\nConcurrent maintainer edit.\n"
    planner = workflow.plan_outputs

    def edit_before_snapshot(directory: Path, outputs: dict[str, bytes], *, inputs: dict[str, bytes]) -> PublicationPlan:
        document.write_bytes(edited)
        return planner(directory, outputs, inputs=inputs)

    monkeypatch.setattr(workflow, "plan_outputs", edit_before_snapshot)
    with pytest.raises(ValueError, match="README changed during rendering"):
        workflow.publish_plot(root, "lu_solve", "median", log_y=True, update_readme=True)
    assert document.read_bytes() == edited
    assert not (root / "docs/assets/bench").exists()


def test_existing_tag_cannot_receive_different_measured_assets(complete_run: tuple[Path, Evidence]) -> None:
    root, evidence = complete_run
    workflow.save_scratch(root, root / "target/bench-reports/performance", evidence)
    workflow.promote(root, root / "target/bench-reports/performance")
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    assets = root / "docs/assets/bench"
    assets.mkdir(parents=True)
    for suffix in (".svg", ".csv", ".provenance.json"):
        (assets / f"vs_linalg_lu_solve_median{suffix}").write_bytes(b"previous tagged artifact\n")
    run_command("git", ["--no-pager", "add", "docs/assets/bench"], cwd=root)
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
            "--quiet",
            "-m",
            "previous artifacts",
        ],
        cwd=root,
    )
    run_command("git", ["--no-pager", "-c", "tag.gpgsign=false", "tag", "v0.4.6"], cwd=root)
    originals = {path: path.read_bytes() for path in (*assets.iterdir(), root / "README.md")}
    with pytest.raises(ValueError, match="exact publication artifact"):
        workflow.publish_plot(root, "lu_solve", "median", log_y=True, update_readme=True)
    assert {path: path.read_bytes() for path in originals} == originals


@pytest.mark.parametrize("scratch_present", [False, True])
@pytest.mark.parametrize("changed", ["source-membership", "unindexed-run", "selected-evidence"])
def test_dimension_publication_rechecks_membership_and_retained_selection(
    complete_run: tuple[Path, Evidence], monkeypatch: pytest.MonkeyPatch, changed: str, scratch_present: bool
) -> None:
    root, evidence = complete_run
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    if not scratch_present:
        shutil.rmtree(root / "target")
    original = (root / "README.md").read_bytes()
    publish = workflow.publish_publication

    def change_before_publication(plan: PublicationPlan) -> tuple[Path, ...]:
        if changed == "source-membership":
            (root / "src/added.rs").write_text("// added after planning\n", encoding="utf-8")
        elif changed == "selected-evidence":
            # Scratch is an input only while present; the archive remains the
            # selected source after cleanup.
            if scratch_present:
                policy.scratch_paths(stem)[0].write_bytes(b"{}\n")
            else:
                (root / policy.archive_directory(root) / "runs" / run_identity(evidence) / "run.json").write_bytes(b"{}\n")
        else:
            unexpected = root / policy.archive_directory(root) / "runs" / ("f" * 64)
            unexpected.mkdir()
            (unexpected / "run.json").write_bytes(b"{}\n")
        return publish(plan)

    monkeypatch.setattr(workflow, "_repo_root", lambda: root)
    monkeypatch.setattr(workflow, "publish_publication", change_before_publication)
    assert workflow.dimension_main(["--update-readme"]) == 2
    assert (root / "README.md").read_bytes() == original
    assert not (root / "docs/assets/bench").exists()


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
    hashes = {
        "performance.csv": "b263bb34bcfd43b8fc8109110db9f2e908516aef6c01ce16f70f5ac18dc7ef2d",
        "performance.full.csv": "4fa5cc13159f5dc91e711ae67554845bcbd1ace2cfc3e31487fc85f75d3dc459",
        "performance.full.provenance.json": "a43f66e52035aecd8d235dfded916b23c3240529d0c07f793d27322521cfb0d0",
        "performance.provenance.json": "89144a2899ba73f55d009da4fa20b4e6477600902b93cde469383f719c497508",
    }
    assert {name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in hashes} == hashes

    retained = convert_csv(
        (directory / "performance.csv").read_bytes(),
        (directory / "performance.provenance.json").read_bytes(),
        (ROOT / "tooling/performance-legacy.toml").read_bytes(),
    )
    comparison = parse_comparison(retained.payload)
    assert len(comparison.comparisons) == 190
    assert len(comparison.missing_baseline) == 35
    assert not comparison.missing_current
    assert "rational_input_d8/solve_big_rational_gaussian" in comparison.missing_baseline
    assert dict(retained.sources)["baseline"].source_sha256 is None
    assert dict(dict(retained.sources)["baseline"].context)["legacy.payload-sha256"] == hashes["performance.csv"]


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
    monkeypatch.chdir(tmp_path)
    pair = ReleasePair("v0.4.6", baseline_tag)
    plan = policy.measurement_plan(root, pair, "all", scope)
    assert isinstance(plan.current.policy, CompletePolicy)
    assert isinstance(plan.baseline.policy, CompletePolicy)
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


def test_local_explicit_labels_and_release_eligibility(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    assert workflow.select_pair(ROOT, "explicit", "v0.4.6", "v0.4.6", local=True) == ReleasePair("v0.4.6", "v0.4.6")
    with pytest.raises(ValueError, match="must differ"):
        workflow.select_pair(ROOT, "explicit", "v0.4.6", "v0.4.6")
    with pytest.raises(ValueError, match=r"v0.4.4|unsupported"):
        workflow.select_pair(ROOT, "explicit", "v0.4.6", "v0.4.3", local=True)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("minimum-release", "[]"),
        ("minimum-release", "[0, 4]"),
        ("minimum-release", "[0, 4, true]"),
        ("minimum-release", "[0, 4, -1]"),
        ("modern-api-release", "[]"),
        ("modern-api-release", "[0, 4, 3]"),
        ("legacy-api-cfg", '"--arbitrary-flag"'),
    ],
)
def test_malformed_release_policy_is_rejected(tmp_path: Path, field: str, value: str) -> None:
    directory = tmp_path / "tooling"
    directory.mkdir()
    text = (ROOT / policy.POLICY_CONFIG).read_text(encoding="utf-8")
    import re  # noqa: PLC0415

    text = re.sub(rf"(?m)^{re.escape(field)} = .*$", f"{field} = {value}", text)
    (tmp_path / policy.POLICY_CONFIG).write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=r"release|cfg"):
        policy.compatibility_adapter(PAIR, tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "diagnostic"),
    [
        ("dimensions", "[]", "dimensions"),
        ("dimensions", "[2, true]", "dimensions"),
        ("dimensions", "[2, -3]", "dimensions"),
        ("dimensions", "[3, 2]", "dimensions"),
        ("dimensions", "[2, 2]", "dimensions"),
        ("library-benchmarks", '"la_stack_lu"', "library-benchmarks"),
        ("library-benchmarks", '["la_stack_lu", "la_stack_lu"]', "library-benchmarks"),
        ("unavailable-prefixes", "[false]", "unavailable-prefixes"),
        ("unavailable-benchmarks", '["../escape"]', "unavailable-benchmarks"),
        ("unavailable-benchmark-prefixes", "[1]", "unavailable-benchmark-prefixes"),
        ("exact_d2", '"det"', "exact.exact_d2"),
        ('"8"', "[false]", "additional.8"),
        ("la_stack_lu", '["nalgebra_lu"]', "peers.la_stack_lu"),
        ("la_stack_lu", '["nalgebra_lu", "faer_lu", "extra"]', "peers.la_stack_lu"),
        ("la_stack_lu", '["nalgebra_lu", "nalgebra_lu"]', "peers.la_stack_lu"),
        ("title", '"   "', "metrics.det_via_lu.title"),
        ("benchmarks", '["la_stack_det_via_lu", "nalgebra_det_via_lu"]', "metrics.det_via_lu.benchmarks"),
        ("benchmarks", '["la_stack_det_via_lu", "nalgebra_det_via_lu", "wrong_peer"]', "declared peer pair"),
    ],
)
def test_malformed_scientific_catalog_is_rejected(tmp_path: Path, field: str, value: str, diagnostic: str) -> None:
    """Bad declarations fail at admission, before inventory indexing or timing."""
    import re  # noqa: PLC0415

    (tmp_path / "tooling").mkdir()
    text = (ROOT / policy.POLICY_CONFIG).read_text(encoding="utf-8")
    text, count = re.subn(rf"(?m)^{re.escape(field)} = .*$", f"{field} = {value}", text, count=1)
    assert count == 1
    (tmp_path / policy.POLICY_CONFIG).write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=re.escape(diagnostic)):
        policy.compatibility_adapter(PAIR, tmp_path)


@pytest.mark.parametrize(
    ("original", "replacement", "diagnostic"),
    [
        ("[exact]", "[renamed_exact]", "exact requires a named catalog table"),
        ('"8" =', '"7" =', "additional benchmark dimensions"),
        ("la_stack_lu =", "unknown_benchmark =", "peers must refer"),
        ("[metrics.det_via_lu]", "[metrics]", "metrics.title requires a named catalog table"),
    ],
)
def test_malformed_scientific_tables_are_rejected(tmp_path: Path, original: str, replacement: str, diagnostic: str) -> None:
    """Missing tables and cross-catalog inconsistencies retain field diagnostics."""
    (tmp_path / "tooling").mkdir()
    text = (ROOT / policy.POLICY_CONFIG).read_text(encoding="utf-8")
    assert original in text
    (tmp_path / policy.POLICY_CONFIG).write_text(text.replace(original, replacement, 1), encoding="utf-8")
    with pytest.raises(ValueError, match=diagnostic):
        policy.scientific_policy(tmp_path)


def test_parsed_scientific_catalog_preserves_immutable_peer_and_metric_rows() -> None:
    """Shared planners consume stable, typed catalog rows rather than raw TOML."""
    catalog = policy.scientific_policy(ROOT)
    assert catalog.dimensions == (2, 3, 4, 5, 8, 16, 32, 64)
    assert catalog.peers["la_stack_lu"] == ("nalgebra_lu", "faer_lu")
    assert catalog.metrics["lu"].benchmarks == ("la_stack_lu", "nalgebra_lu", "faer_lu")
    from typing import cast  # noqa: PLC0415

    with pytest.raises(TypeError):
        cast("dict[str, tuple[str, str]]", catalog.peers)["la_stack_lu"] = ("wrong", "peer")


@pytest.mark.parametrize("flag", [("--suite", "exact"), ("--scope", "all-benches"), ("--stem", "ignored")])
def test_assets_rejects_retired_options(flag: tuple[str, str], capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as error:
        workflow.main(["assets", *flag])
    assert error.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err


@pytest.mark.parametrize("command", ["plot", "local"])
def test_cli_help_and_missing_checkout_diagnostics(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], command: str) -> None:
    monkeypatch.chdir(tmp_path)
    cli = workflow.dimension_main if command == "plot" else workflow.main
    arguments = [] if command == "plot" else [command]
    with pytest.raises(SystemExit) as error:
        cli([*arguments, "--help"])
    assert error.value.code == 0
    output = capsys.readouterr()
    assert "usage:" in output.out
    assert not output.err
    assert cli(arguments) == (2 if command == "plot" else 1)
    output = capsys.readouterr()
    assert "Could not locate Cargo.toml" in output.err
    assert "Traceback" not in output.err


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
        workflow.compare_assets(tmp_path, PAIR, output)
    assert output.read_bytes() == b"previous report\n"


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
