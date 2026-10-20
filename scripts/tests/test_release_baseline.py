"""Exercise complete release archives and the workflow's publication barriers."""

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import textwrap
import tomllib
from pathlib import Path

import pytest
from research_repo_tools.criterion_inventory import CriterionInventory, load_discovery, serialize_inventory

import archive_performance as workflow
from performance_runs import required_report_ids

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLING_REQUIREMENT = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["dependency-groups"]["tooling"][0]
WORKFLOW = REPO_ROOT / ".github/workflows/release-benchmarks.yml"
PREPARATION = REPO_ROOT / ".github/actions/prepare-release-benchmarks/action.yml"
IDS = {"vs_linalg": ["d2/la_stack_dot", "diagnostic/peer"], "exact": ["exact_d2/det_exact"]}
BASELINE = "v1.2.3"
RAW_FILES = ("benchmark.json", "estimates.json", "sample.json", "tukey.json")


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) + "\n", encoding="utf-8", newline="\n")


@pytest.fixture
def dataset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "Cargo.toml").write_text('[package]\nname = "fixture"\nversion = "1.2.3"\n', encoding="utf-8")
    tooling = tmp_path / "tooling"
    tooling.mkdir()
    shutil.copyfile(REPO_ROOT / "tooling/criterion-saved.toml", tooling / "criterion-saved.toml")
    (tooling / "criterion-inventory.toml").write_text(
        'schema = 1\nrequired = ["d2/la_stack_dot", "exact_d2/det_exact"]\n[[suites]]\nname = "release"\ncommand = ["unused"]\n', encoding="utf-8"
    )
    criterion = tmp_path / "target/criterion"
    manifest = tmp_path / "target/release-benchmark-inventory.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_bytes(serialize_inventory(CriterionInventory(tuple((name, tuple(ids)) for name, ids in IDS.items()))))
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


def shell_step(name: str, source: Path = WORKFLOW) -> str:
    """Extract a literal run block so tests exercise the workflow's actual shell."""
    text = source.read_text(encoding="utf-8")
    indent = "    " if source == PREPARATION else "      "
    step = text.split(f"{indent}- name: {name}\n", 1)[1].split(f"\n{indent}- ", 1)[0]
    lines: list[str] = []
    for line in step.split(f"{indent}  run: |\n", 1)[1].splitlines():
        if line.strip() and not line.startswith(f"{indent}    "):
            break
        lines.append(line)
    return textwrap.dedent("\n".join(lines)).rstrip()


def run_shell(script: str, cwd: Path, extra_env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    bash = shutil.which("bash")
    assert bash is not None
    return subprocess.run(  # noqa: S603 - run only checked-in workflow shell in an isolated fixture.
        [bash, "--noprofile", "--norc", "-c", script],
        cwd=cwd,
        env={**{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}}, **extra_env},
        capture_output=True,
        encoding="utf-8",
        check=False,
        timeout=30,
    )


def test_complete_dataset_round_trips_through_workflow_archive(dataset: tuple[Path, Path], tmp_path: Path) -> None:
    criterion, _manifest = dataset
    assert workflow.main(["release-check", BASELINE]) == 0
    output = tmp_path / "step outputs"
    # Git Bash accepts forward-slash drive paths on Windows as well as POSIX paths.
    result = run_shell(
        shell_step("Package release Criterion baseline"),
        tmp_path,
        {"RELEASE_TAG": BASELINE, "GITHUB_OUTPUT": output.as_posix(), "GITHUB_RUN_ID": "123", "GITHUB_RUN_ATTEMPT": "2"},
    )
    assert result.returncode == 0, result.stderr
    asset = f"la-stack-{BASELINE}-criterion-baseline.tar.gz"
    assert output.read_text(encoding="utf-8") == f"asset={asset}\nartifact-name=bench-baseline-{BASELINE}-123-2\n"
    with tarfile.open(tmp_path / asset) as archive:
        for ids in IDS.values():
            for benchmark in ids:
                for sample in ("new", BASELINE):
                    for filename in RAW_FILES:
                        member = archive.extractfile(f"criterion/{benchmark}/{sample}/{filename}")
                        assert member is not None
                        assert member.read() == (criterion / benchmark / sample / filename).read_bytes()
        member = archive.extractfile("criterion/release-benchmark-inventory.json")
        assert member is not None
        assert json.loads(member.read())["suites"] == IDS


def test_saved_baseline_must_match_fresh_samples(dataset: tuple[Path, Path]) -> None:
    criterion, _manifest = dataset
    write_json(criterion / "d2/la_stack_dot" / BASELINE / "sample.json", {"iters": [1.0] * 100, "times": [3.0] * 100})
    assert workflow.main(["release-check", BASELINE]) == 1


@pytest.mark.parametrize("relocation", ["moved", "swapped"])
def test_report_consumers_require_canonical_paths(dataset: tuple[Path, Path], capsys: pytest.CaptureFixture[str], relocation: str) -> None:
    criterion, _manifest = dataset
    original = criterion / "d2/la_stack_dot"
    moved = original.rename(criterion / "d2/misplaced")
    if relocation == "swapped":
        peer = criterion / "exact_d2/det_exact"
        peer.rename(original)
        moved.rename(peer)
    before = {path: path.read_bytes() for path in criterion.rglob("*") if path.is_file()}

    assert workflow.main(["release-check", BASELINE]) == 1

    output = capsys.readouterr()
    assert output.out == ""
    assert "non-canonical report benchmark storage path" in output.err
    assert "d2/la_stack_dot" in output.err
    assert {path: path.read_bytes() for path in criterion.rglob("*") if path.is_file()} == before


def test_configured_inventory_includes_release_consumers_and_peers() -> None:
    expected = required_report_ids(REPO_ROOT)
    assert {"d64/faer_lu_solve", "d64/nalgebra_lu_solve", "exact_hilbert_5x5/solve_exact", "rational_input_d8/solve_big_rational_gaussian"} <= expected
    assert "d64/faer_dot_native" in expected
    assert "d64/faer_dot" not in expected
    discovery = load_discovery(REPO_ROOT, "tooling/criterion-inventory.toml")
    assert {suite.name for suite in discovery.suites} == {"vs_linalg", "exact"}
    assert all(suite.command[-1] == "list" for suite in discovery.suites)


def test_workflow_requires_preflight_and_isolates_publication_permissions() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")
    preparation = PREPARATION.read_text(encoding="utf-8")
    assert preparation.index("- name: Inventory full release suites") < preparation.index("- name: Validate benchmark inputs")
    assert "run: just test-bench-inputs" in preparation
    assert "run: just bench-release-inventory" in preparation
    assert "continue-on-error" not in preparation
    assert preparation.index("- name: Require fresh release preflight") < preparation.index("- name: Set up declared tools without caches")
    assert "uses: ./.github/actions/setup-tools" in preparation
    assert "cache: false" in preparation
    steps = [
        "Prepare release benchmarks",
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
    assert "required: true" in workflow
    assert "\n  release:" not in workflow
    assert "group: release-benchmarks-${{ inputs.tag }}" in workflow
    assert "cancel-in-progress: false" in workflow
    preflight, remainder = workflow.split("  release-baseline:\n", 1)
    producer, publisher = remainder.split("  publish-baseline:\n", 1)
    assert "needs: validate-release" in producer
    assert "validated-attempt: ${{ github.run_attempt }}" in preflight
    assert "validated-attempt: ${{ needs.validate-release.outputs.validated-attempt }}" in producer
    assert "ref: ${{ github.sha }}" in producer
    assert "ref: ${{ needs.validate-release.outputs.commit }}" not in producer
    assert "contents: write" not in producer
    assert "persist-credentials: false" in producer
    assert "contents: write" in preflight
    assert "actions/checkout@" not in preflight + publisher
    assert "target: ${{ steps.target.outputs.target }}" in preflight
    assert "RELEASE_TARGET: ${{ needs.validate-release.outputs.target }}" in publisher
    assert '--expected-title "$RELEASE_TAG"' in preflight
    assert '--target "$RUNNER_TEMP/release-target.json"' in publisher
    assert "--publish" not in publisher
    assert (preflight + publisher).count(f"--from {TOOLING_REQUIREMENT}") == 2
    assert "uses: ./" not in preflight + publisher
    assert "--clobber" not in workflow
    assert "needs: [validate-release, release-baseline]" in publisher
    assert "GH_REPO: ${{ github.repository }}" in publisher
    assert "contents: write" in publisher
    assert "name: ${{ needs.release-baseline.outputs.artifact-name }}" in publisher


def test_setup_limits_preserve_benchmark_and_tail_budgets() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")
    producer = workflow.split("  release-baseline:\n", 1)[1].split("  publish-baseline:\n", 1)[0]
    job_limit = re.search(r"(?m)^    timeout-minutes: (\d+)$", producer)
    assert job_limit is not None
    setup = producer.split("    steps:\n", 1)[1].split("      - name: Save comparative Criterion baseline", 1)[0]
    setup_steps = re.findall(r"(?ms)^      - .*?(?=^      - |\Z)", setup)
    assert setup_steps
    setup_limits: list[int] = []
    for step in setup_steps:
        limit = re.search(r"(?m)^        timeout-minutes: (\d+)$", step)
        assert limit is not None, f"unbounded setup step: {step}"
        setup_limits.append(int(limit[1]))
    assert all(limit > 0 for limit in setup_limits)
    discovery = load_discovery(REPO_ROOT, "tooling/criterion-inventory.toml")
    discovery_minutes = sum(suite.timeout for suite in discovery.suites) / 60
    # The failed hosted run spent 27 minutes on uncached tool installation.
    # Reserve 35 for tools, each declared cold discovery deadline, and 5 for inputs.
    assert setup_limits[-1] >= 35 + discovery_minutes + 5
    assert "uses: ./.github/actions/prepare-release-benchmarks" in setup
    assert "using: composite" in PREPARATION.read_text(encoding="utf-8")
    benchmark_limits = [
        int(limit) for limit in re.findall(r"(?ms)^      - name: Save (?:comparative|exact) Criterion baseline\n.*?^        timeout-minutes: (\d+)$", producer)
    ]
    assert benchmark_limits == [150, 90]
    assert int(job_limit[1]) - sum(setup_limits) - sum(benchmark_limits) >= 15
    assert int(job_limit[1]) <= 360


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
    assert "| vs_linalg | failure |" in summary.read_text(encoding="utf-8")
    assert "| exact | skipped |" in summary.read_text(encoding="utf-8")


@pytest.mark.parametrize("validated_attempt", ["", "1", "2"])
def test_benchmark_retry_requires_current_preflight(tmp_path: Path, validated_attempt: str) -> None:
    result = run_shell(
        shell_step("Require fresh release preflight", PREPARATION),
        tmp_path,
        {"VALIDATED_ATTEMPT": validated_attempt, "GITHUB_RUN_ATTEMPT": "2"},
    )
    assert result.returncode == (0 if validated_attempt == "2" else 1)


@pytest.mark.parametrize("ref", ["refs/heads/main", f"refs/heads/{BASELINE}", "refs/tags/v9.9.9", ""])
def test_preflight_rejects_dispatch_outside_release_tag(tmp_path: Path, ref: str) -> None:
    output = tmp_path / "outputs"
    result = run_shell(shell_step("Validate draft release target"), tmp_path, {"RELEASE_TAG": BASELINE, "GITHUB_REF": ref, "GITHUB_OUTPUT": output.as_posix()})
    assert result.returncode == 1
    assert f"Dispatch with --ref {BASELINE}" in result.stdout
    assert not output.exists()


GH_FIXTURE = r"""import hashlib, json, os, sys
from pathlib import Path
args = sys.argv[1:]
root = Path(os.environ["RELEASE_FIXTURE"])
tag = os.environ["RELEASE_TAG"]
asset = os.environ["RELEASE_ASSET"]
with (root / "api-calls").open("a", encoding="utf-8") as stream:
    stream.write(json.dumps(args) + "\n")
metadata = {"id": 7, "name": asset, "size": len((root / asset).read_bytes()),
            "digest": "sha256:" + hashlib.sha256((root / asset).read_bytes()).hexdigest(), "state": "uploaded"}
assets = [metadata] if (root / "uploaded").exists() else []
release = {"id": 42, "url": "https://api.github.com/repos/owner/repo/releases/42", "tag_name": tag,
           "name": tag, "draft": True, "prerelease": False, "immutable": False, "assets": assets}
if any(item.endswith("/releases?per_page=100") for item in args):
    value = [[release]]
elif args[0] == "release":
    value = {"databaseId": 42, "tagName": tag}
elif any("uploads.github.com" in item for item in args):
    (root / "uploaded").touch()
    value = {}
elif "PATCH" in args:
    (root / "published").touch()
    value = dict(release, draft=False, immutable=True)
elif args[1].endswith("/releases/assets/7"):
    sys.stdout.buffer.write((root / asset).read_bytes())
    raise SystemExit(0)
elif "/git/ref/tags/" in args[1]:
    value = {"ref": "refs/tags/" + tag, "object": {"type": "commit", "sha": os.environ.get("REMOTE_COMMIT", "a" * 40)}}
elif any("/assets?" in item for item in args):
    value = [assets]
else:
    value = release
print(json.dumps(value))
"""


def test_installed_cli_binds_trusted_target_and_leaves_draft_for_review(tmp_path: Path) -> None:
    """Exercise actual preflight and writer shell with the pinned installed CLI."""
    fixture = tmp_path / "provider with spaces"
    fixture.mkdir()
    program = fixture / "gh_fixture.py"
    program.write_text(GH_FIXTURE, encoding="utf-8")
    if os.name == "nt":
        executable = fixture / "gh.cmd"
        executable.write_text(f'@"{sys.executable}" "%~dp0gh_fixture.py" %*\r\n', encoding="utf-8", newline="")
    else:
        executable = fixture / "gh"
        executable.write_text(f'#!/bin/sh\nexec {shlex.quote(sys.executable)} {shlex.quote(str(program))} "$@"\n', encoding="utf-8")
        executable.chmod(0o755)
    asset = f"la-stack-{BASELINE}-criterion-baseline.tar.gz"
    (fixture / asset).write_bytes(b"raw archive\r\n\xff")
    output = fixture / "outputs"
    env = {
        "PATH": str(fixture) + os.pathsep + os.environ["PATH"],
        "RELEASE_FIXTURE": str(fixture),
        "RELEASE_TAG": BASELINE,
        "RELEASE_ASSET": asset,
        "GH_REPO": "owner/repo",
        "GITHUB_REF": f"refs/tags/{BASELINE}",
        "GITHUB_SHA": "a" * 40,
        "GITHUB_OUTPUT": output.as_posix(),
        "RUNNER_TEMP": fixture.as_posix(),
        "GITHUB_WORKSPACE": fixture.as_posix(),
    }
    # Substitute only uvx's environment acquisition; retain its exact package
    # arguments and execute the already locked package in a real child process.
    runner = (
        'uvx() { [[ "$1 $2 $3 $4 $5" == "--python 3.15 --from '
        f'{TOOLING_REQUIREMENT} research-repo-tools" ]] || return 90; '
        f'shift 5; {shlex.quote(sys.executable)} -m research_repo_tools "$@"; }}\n'
    )
    result = run_shell(runner + shell_step("Validate draft release target"), fixture, env)
    assert result.returncode == 0, result.stderr
    target = output.read_text(encoding="utf-8").strip().removeprefix("target=")
    assert json.loads(target)["release_id"] == 42
    assert json.loads(target)["commit"] == "a" * 40
    assert json.loads(target)["expected_title"] == BASELINE
    stale = run_shell(runner + shell_step("Attach baseline for review"), fixture, {**env, "RELEASE_TARGET": target, "REMOTE_COMMIT": "b" * 40})
    assert stale.returncode != 0
    assert not (fixture / "published").exists()
    result = run_shell(runner + shell_step("Attach baseline for review"), fixture, {**env, "RELEASE_TARGET": target})
    assert result.returncode == 0, result.stderr
    assert (fixture / "uploaded").exists()
    assert not (fixture / "published").exists()
    calls = [json.loads(line) for line in (fixture / "api-calls").read_text(encoding="utf-8").splitlines()]
    assert not any("PATCH" in call for call in calls)
    assert any("/releases/assets/7" in item for call in calls for item in call)
