"""Installed-package acceptance of la-stack's complete performance workflow."""

import hashlib
import io
import json
import shutil
import tarfile
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from research_repo_tools.common_measurement import CommonHarnessPlan, measure_prepared_pair
from research_repo_tools.complete_runs import RunSeries, run_identity, serialize_run
from research_repo_tools.process import run_git_command
from research_repo_tools.release_pairs import ReleasePair

import archive_performance as workflow
import criterion_dim_plot as plot
import performance_runs as policy
from bench_compare import render_release_artifacts
from performance_artifacts import ArtifactPaths
from performance_phase import commands
from release_baseline import required_report_ids

if TYPE_CHECKING:
    from research_repo_tools.evidence import Evidence

ROOT = Path(__file__).resolve().parents[2]
PAIR = ReleasePair("v0.4.6", "v0.4.5")
DRIVER = """import json, os, sys
from pathlib import Path
action, suite, phase = sys.argv[1:]
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
        for pattern in policy.HARNESS:
            for source in ROOT.glob(pattern):
                destination = root / source.relative_to(ROOT)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
        (root / "scripts/performance_phase.py").write_text(f"INVENTORIES = {inventories!r}\n" + DRIVER, encoding="utf-8")
        run_git_command(["init", "--quiet"], cwd=root)
        run_git_command(
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
    plan = policy.declared_plan(PAIR, "all", "release-signal", inventories, baseline_env=(("FIXTURE_FAULT", fault),), current_env=())
    return baseline, current, plan


def measured(tmp_path: Path) -> tuple[Path, Evidence]:
    baseline, current, plan = prepared(tmp_path)
    evidence = measure_prepared_pair(baseline, current, plan, PAIR, working_tree=True)
    (current / ".config").mkdir(exist_ok=True)
    shutil.copyfile(ROOT / policy.REPORT_CONFIG, current / policy.REPORT_CONFIG)
    return current, evidence


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
    report = policy.render_scientific_report(evidence)
    assert b"not rerun in the current phase" in report
    assert b"baseline API unavailable" in report
    assert b"paired confidence interval" in report


@pytest.mark.parametrize("fault", ["gate", "stale", "missing-row", "missing-statistic", "missing-sample"])
def test_declared_gates_and_completeness_fail_closed(tmp_path: Path, fault: str) -> None:
    baseline, current, plan = prepared(tmp_path, fault)
    with pytest.raises(ValueError, match=r"gate|changed|incomplete|sample|statistic|expected|missing|mean"):
        measure_prepared_pair(baseline, current, plan, PAIR, working_tree=True)
    assert not (current / "target/criterion").exists()
    if fault in {"gate", "stale"}:
        assert not (baseline / "target/criterion").exists()


def test_shared_retention_repeated_runs_and_offline_report(tmp_path: Path) -> None:
    root, first = measured(tmp_path)
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, first)
    workflow.promote(root, stem)
    first_report = (root / "docs/performance.md").read_bytes()
    retained = root / policy.ARCHIVE / "runs" / run_identity(first)
    original = {path: path.read_bytes() for path in retained.iterdir()}
    # A second valid measurement with the same release pair has new provenance.
    sources = tuple((phase, replace(source, context=(*source.context, ("experiment", "second")))) for phase, source in first.sources)
    second = replace(first, sources=sources)
    workflow.save_scratch(root, stem, second)
    workflow.promote(root, stem)
    assert len(list((root / policy.ARCHIVE / "runs").iterdir())) == 2
    assert {path: path.read_bytes() for path in original} == original
    assert (root / "docs/archive/performance/v0.4.6-vs-v0.4.5.md").read_bytes() == first_report
    shutil.rmtree(root / "target")
    assert run_identity(policy.retained_run(root, stem)) == run_identity(second)
    before = (root / "docs/performance.md").read_bytes()
    workflow.promote(root, stem)
    assert (root / "docs/performance.md").read_bytes() == before
    (root / policy.ARCHIVE / "latest.json").write_text("{}", encoding="utf-8")
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


def test_readme_uses_retained_reference_estimates_after_cleanup(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root, evidence = measured(tmp_path)
    stem = root / "target/bench-reports/performance"
    workflow.save_scratch(root, stem, evidence)
    workflow.promote(root, stem)
    shutil.copyfile(ROOT / "README.md", root / "README.md")
    shutil.rmtree(root / "target")
    monkeypatch.setattr(plot, "_repo_root", lambda: root)
    monkeypatch.setattr(plot, "_render_svg_with_gnuplot", lambda request: request.out_svg.write_bytes(b"<svg/>\n"))
    assert plot.main(["--update-readme", "--log-y"]) == 0
    assets = root / "docs/assets/bench"
    csv = assets / "vs_linalg_lu_solve_median.csv"
    assert "2,8.0,7.2,8.8,10.0,9.0,11.0,10.0,9.0,11.0" in csv.read_text(encoding="utf-8")
    provenance = json.loads(csv.with_suffix(".provenance.json").read_bytes())
    assert provenance["measurement"]["peer_phase"] == "baseline"
    assert provenance["run_id"] == run_identity(evidence)
    before = {path: path.read_bytes() for path in (*assets.iterdir(), root / "README.md")}
    assert plot.main(["--update-readme", "--log-y"]) == 0
    assert {path: path.read_bytes() for path in before} == before


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
    assert all("--locked" in command and "--noplot" in command for command in commands("all", "baseline"))
    assert "la_stack" in commands("vs_linalg", "current")[0]
    assert "la_stack" not in commands("vs_linalg", "baseline")[0]
    environment = dict(policy.phase_environment(ROOT, "la_stack_pre_rational_input_api"))
    assert "--cfg=la_stack_pre_rational_input_api" in environment.get("CARGO_ENCODED_RUSTFLAGS", environment.get("RUSTFLAGS", ""))
    assert environment["CARGO_TARGET_DIR"] == "target"
    assert environment["CRITERION_HOME"] == "target/criterion"


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
    (current / "scripts/performance_phase.py").write_text(driver, encoding="utf-8")
    plan = policy.declared_plan(pair, "all", "all-benches", inventories, baseline_env=(), current_env=())
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
        if destination == root / policy.ARCHIVE / "latest.json" and not failed:
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
