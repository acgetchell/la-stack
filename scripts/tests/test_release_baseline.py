"""Exercise complete release archives and the workflow's publication barriers."""

import json
import os
import shutil
import subprocess
import tarfile
import textwrap
from pathlib import Path

import pytest

import release_baseline

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github/workflows/release-benchmarks.yml"
IDS = {"vs_linalg": ["d2/la_stack_dot", "diagnostic/peer"], "exact": ["exact_d2/det_exact"]}
BASELINE = "v1.2.3"


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) + "\n", encoding="utf-8", newline="\n")


@pytest.fixture
def dataset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    monkeypatch.setattr(release_baseline, "required_report_ids", lambda: {"d2/la_stack_dot", "exact_d2/det_exact"})
    criterion = tmp_path / "target/criterion"
    manifest = tmp_path / "target/release-benchmark-inventory.json"
    write_json(manifest, IDS)
    estimate = {"point_estimate": 2.0, "confidence_interval": {"confidence_level": 0.95, "lower_bound": 1.0, "upper_bound": 3.0}}
    for ids in IDS.values():
        for benchmark in ids:
            for sample in ("new", BASELINE):
                directory = criterion / benchmark / sample
                write_json(directory / "benchmark.json", {"full_id": benchmark})
                write_json(directory / "sample.json", {"iters": [1.0] * 100, "times": [2.0] * 100})
                write_json(directory / "estimates.json", {"mean": estimate, "median": estimate})
                write_json(directory / "tukey.json", [1.0, 2.0, 3.0, 4.0])
    return criterion, manifest


def shell_step(name: str) -> str:
    """Extract a literal run block so tests exercise the workflow's actual shell."""
    text = WORKFLOW.read_text(encoding="utf-8")
    step = text.split(f"      - name: {name}\n", 1)[1].split("\n      - ", 1)[0]
    lines: list[str] = []
    for line in step.split("        run: |\n", 1)[1].splitlines():
        if line.strip() and not line.startswith("          "):
            break
        lines.append(line)
    return textwrap.dedent("\n".join(lines)).rstrip()


def run_shell(script: str, cwd: Path, extra_env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    bash = shutil.which("bash")
    assert bash is not None
    return subprocess.run(  # noqa: S603 - run only checked-in workflow shell in an isolated fixture.
        [bash, "--noprofile", "--norc", "-c", script],
        cwd=cwd,
        env={**os.environ, **extra_env},
        capture_output=True,
        encoding="utf-8",
        check=False,
        timeout=30,
    )


def test_complete_dataset_round_trips_through_workflow_archive(dataset: tuple[Path, Path], tmp_path: Path) -> None:
    criterion, manifest = dataset
    assert release_baseline.validate(criterion, manifest, BASELINE) == 3
    output = tmp_path / "step outputs"
    # Git Bash accepts forward-slash drive paths on Windows as well as POSIX paths.
    result = run_shell(shell_step("Package release Criterion baseline"), tmp_path, {"RELEASE_TAG": BASELINE, "GITHUB_OUTPUT": output.as_posix()})
    assert result.returncode == 0, result.stderr
    asset = f"la-stack-{BASELINE}-criterion-baseline.tar.gz"
    assert output.read_text() == f"asset={asset}\n"
    with tarfile.open(tmp_path / asset) as archive:
        for ids in IDS.values():
            for benchmark in ids:
                for sample in ("new", BASELINE):
                    for filename in release_baseline.RAW_FILES:
                        member = archive.extractfile(f"criterion/{benchmark}/{sample}/{filename}")
                        assert member is not None
                        assert member.read() == (criterion / benchmark / sample / filename).read_bytes()
        member = archive.extractfile("criterion/release-benchmark-inventory.json")
        assert member is not None
        assert json.loads(member.read()) == IDS


@pytest.mark.parametrize("filename", release_baseline.RAW_FILES)
@pytest.mark.parametrize("sample", ["new", BASELINE])
def test_missing_raw_file_blocks_publication(dataset: tuple[Path, Path], sample: str, filename: str) -> None:
    criterion, manifest = dataset
    missing = criterion / "diagnostic/peer" / sample / filename
    missing.unlink()
    if sample == "new" and filename == "benchmark.json":
        with pytest.raises(ValueError, match="incomplete Criterion dataset: diagnostic/peer"):
            release_baseline.validate(criterion, manifest, BASELINE)
    else:
        with pytest.raises(FileNotFoundError) as error:
            release_baseline.validate(criterion, manifest, BASELINE)
        assert error.value.filename == str(missing)


def test_missing_entire_suite_blocks_publication(dataset: tuple[Path, Path]) -> None:
    criterion, manifest = dataset
    shutil.rmtree(criterion / "exact_d2")
    with pytest.raises(ValueError, match="incomplete Criterion dataset: exact_d2/det_exact"):
        release_baseline.validate(criterion, manifest, BASELINE)


def test_saved_baseline_must_match_fresh_samples(dataset: tuple[Path, Path]) -> None:
    criterion, manifest = dataset
    write_json(criterion / "d2/la_stack_dot" / BASELINE / "sample.json", {"iters": [1.0] * 100, "times": [3.0] * 100})
    with pytest.raises(ValueError, match="saved baseline differs"):
        release_baseline.validate(criterion, manifest, BASELINE)


def test_report_consumers_require_canonical_paths(dataset: tuple[Path, Path]) -> None:
    criterion, manifest = dataset
    (criterion / "d2/la_stack_dot").rename(criterion / "d2/misplaced")
    with pytest.raises(ValueError, match="outside its consumer path"):
        release_baseline.validate(criterion, manifest, BASELINE)


@pytest.mark.parametrize("corruption", ["short-samples", "nonfinite", "reversed-interval", "wrong-confidence", "malformed"])
def test_invalid_measurement_blocks_publication(dataset: tuple[Path, Path], corruption: str) -> None:
    criterion, manifest = dataset
    for sample in ("new", BASELINE):
        directory = criterion / "d2/la_stack_dot" / sample
        if corruption == "malformed":
            (directory / "estimates.json").write_text("{broken", encoding="utf-8", newline="\n")
        elif corruption in {"short-samples", "nonfinite"}:
            times = [2.0] * 99 if corruption == "short-samples" else [float("inf")] * 100
            write_json(directory / "sample.json", {"iters": [1.0] * 100, "times": times})
        else:
            estimates = json.loads((directory / "estimates.json").read_text())
            interval = estimates["median"]["confidence_interval"]
            interval["lower_bound" if corruption == "reversed-interval" else "confidence_level"] = 4.0
            write_json(directory / "estimates.json", estimates)
    errors = {
        "short-samples": "times must contain 100 samples",
        "nonfinite": "finite positive timing value",
        "reversed-interval": "confidence interval must be ordered",
        "wrong-confidence": "expected a 95% confidence interval",
        "malformed": "Expecting property name",
    }
    with pytest.raises(ValueError, match=errors[corruption]):
        release_baseline.validate(criterion, manifest, BASELINE)


@pytest.mark.parametrize("baseline", ["../escape", "", "new", "base", "change", "report", "/absolute"])
def test_invalid_baseline_is_rejected(dataset: tuple[Path, Path], baseline: str) -> None:
    with pytest.raises(ValueError, match="safe, nonreserved"):
        release_baseline.validate(*dataset, baseline)


def test_inventory_includes_release_consumers_and_peers() -> None:
    expected = release_baseline.required_report_ids()
    assert {"d64/faer_lu_solve", "d64/nalgebra_lu_solve", "exact_hilbert_5x5/solve_exact", "rational_input_d8/solve_big_rational_gaussian"} <= expected
    with pytest.raises(ValueError, match="inventory omits report consumers"):
        release_baseline.inventory_ids(IDS)


@pytest.mark.parametrize("output", ["", "hello", "a/b: benchmark\na/b: benchmark", "a: benchmark"])
def test_invalid_listing_is_rejected(output: str) -> None:
    with pytest.raises(ValueError, match="unique, nonempty"):
        release_baseline.parse_benchmark_list(output)


@pytest.mark.parametrize("host_newline", ["\n", "\r\n"], ids=["posix", "windows"])
def test_discovery_uses_full_suites_and_never_times_inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, host_newline: str) -> None:
    calls: list[list[str]] = []
    write_text = Path.write_text

    def host_write_text(path: Path, data: str, encoding: str | None = None, errors: str | None = None, newline: str | None = None) -> int:
        # Emulate Windows default text translation on every CI host. An explicit
        # newline policy must preserve identical manifest bytes on both platforms.
        return write_text(path, data, encoding=encoding, errors=errors, newline=host_newline if newline is None else newline)

    def run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        suite = args[args.index("--bench") + 1]
        return subprocess.CompletedProcess(args, 0, stdout="\n".join(f"{name}: benchmark" for name in IDS[suite]))

    monkeypatch.setattr(release_baseline, "run_cargo_command", run)
    monkeypatch.setattr(release_baseline, "required_report_ids", lambda: {"d2/la_stack_dot", "exact_d2/det_exact"})
    monkeypatch.setattr(Path, "write_text", host_write_text)
    manifest = tmp_path / "manifest.json"
    release_baseline.discover(tmp_path, manifest, tmp_path / "criterion")
    assert manifest.read_bytes() == (json.dumps(IDS, indent=2) + "\n").encode("utf-8")
    assert calls == [
        ["bench", "--locked", "--features", "bench", "--bench", "vs_linalg", "--", "--list"],
        ["bench", "--locked", "--features", "bench,exact", "--bench", "exact", "--", "--list"],
    ]


def test_discovery_refuses_stale_measurements(dataset: tuple[Path, Path], tmp_path: Path) -> None:
    criterion, manifest = dataset
    before = manifest.read_bytes()
    with pytest.raises(ValueError, match="fresh Criterion directory"):
        release_baseline.discover(tmp_path, manifest, criterion)
    assert manifest.read_bytes() == before


def test_workflow_fails_closed_and_manual_runs_cannot_publish() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")
    steps = [
        "Validate benchmark inputs",
        "Inventory full release suites",
        "Save comparative Criterion baseline",
        "Save exact Criterion baseline",
        "Validate complete release dataset",
        "Package release Criterion baseline",
        "Upload temporary baseline artifact",
    ]
    positions = [workflow.index(f"- name: {step}") for step in steps]
    assert positions == sorted(positions)
    assert "\n        if:" not in workflow[positions[0] : positions[-1]]
    assert "continue-on-error" not in workflow
    assert "workflow_dispatch:" in workflow
    publisher = workflow.split("  publish-baseline:\n", 1)[1]
    assert "if: ${{ github.event_name == 'release' }}" in publisher
    assert "needs: release-baseline" in publisher
    assert "GH_REPO: ${{ github.repository }}" in publisher
    assert "contents: write" not in workflow.split("  publish-baseline:\n", 1)[0]


@pytest.mark.parametrize(("suite", "step"), [("vs_linalg", "Save comparative Criterion baseline"), ("exact", "Save exact Criterion baseline")])
def test_failed_suite_retains_timing_and_does_not_mark_completion(tmp_path: Path, suite: str, step: str) -> None:
    (tmp_path / "target").mkdir()
    # A shell function substitutes only the benchmark command; timestamps and
    # failure handling execute exactly as checked into the workflow.
    script = "just() { return 7; }\n" + shell_step(step)
    result = run_shell(script, tmp_path, {"RELEASE_TAG": BASELINE})
    assert result.returncode == 7
    assert (tmp_path / f"target/{suite}-started").is_file()
    assert not (tmp_path / f"target/{suite}-finished").exists()
    summary = tmp_path / "suite summary.md"
    result = run_shell(
        shell_step("Suite timing summary"),
        tmp_path,
        {
            "GITHUB_STEP_SUMMARY": summary.as_posix(),
            "COMPARATIVE_OUTCOME": "failure",
            "EXACT_OUTCOME": "skipped",
        },
    )
    assert result.returncode == 0, result.stderr
    assert "| vs_linalg | failure |" in summary.read_text()
    assert "| exact | skipped |" in summary.read_text()
