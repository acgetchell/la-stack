"""Exercise complete release archives and the workflow's publication barriers."""

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest
import yaml
from research_repo_tools.criterion_inventory import CriterionInventory, load_discovery, serialize_inventory
from research_repo_tools.process import run_command

from tests.tooling.support_toolchain import TOOLING_REQUIREMENT

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github/workflows/release-benchmarks.yml"
CATALOG = tomllib.loads((REPO_ROOT / "tooling/performance-policy.toml").read_text(encoding="utf-8"))
REQUIRED = set(CATALOG["release"]["default"] + CATALOG["release"]["exact"])
# Retain optional diagnostic additions without letting a declaration omission
# remove independently required scientific cases from the archive fixture.
DECLARED = set(load_discovery(REPO_ROOT, "tooling/criterion-inventory.toml").required) | set(
    tomllib.loads((REPO_ROOT / "tooling/criterion-saved.toml").read_text(encoding="utf-8"))["required"]
)
IDS = {
    "public-default": sorted(set(CATALOG["release"]["default"])),
    "public-exact": sorted(set(CATALOG["release"]["exact"]) | (DECLARED - REQUIRED)),
}
STORAGE = tomllib.loads((REPO_ROOT / "tooling/criterion-saved.toml").read_text(encoding="utf-8"))["storage-paths"]
BASELINE = "v1.2.3"
RAW_FILES = ("benchmark.json", "estimates.json", "sample.json", "tukey.json")


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) + "\n", encoding="utf-8", newline="\n")


@pytest.fixture
def dataset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    for name in ("pyproject.toml", "uv.lock", ".python-version", "rust-toolchain.toml"):
        shutil.copyfile(REPO_ROOT / name, tmp_path / name)
    (tmp_path / "Cargo.toml").write_text('[package]\nname = "fixture"\nversion = "1.2.3"\n', encoding="utf-8")
    tooling = tmp_path / "tooling"
    tooling.mkdir()
    shutil.copyfile(REPO_ROOT / "tooling/criterion-saved.toml", tooling / "criterion-saved.toml")
    shutil.copyfile(REPO_ROOT / "tooling/criterion-inventory.toml", tooling / "criterion-inventory.toml")
    criterion = tmp_path / "target/criterion"
    manifest = tmp_path / "target/release-benchmark-inventory.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_bytes(serialize_inventory(CriterionInventory(tuple((name, tuple(ids)) for name, ids in IDS.items()))))
    estimate = {"point_estimate": 2.0, "confidence_interval": {"confidence_level": 0.95, "lower_bound": 1.0, "upper_bound": 3.0}}
    for ids in IDS.values():
        for benchmark in ids:
            for sample in ("new", BASELINE):
                directory = criterion / STORAGE.get(benchmark, benchmark) / sample
                write_json(directory / "benchmark.json", {"full_id": benchmark})
                write_json(directory / "sample.json", {"iters": [1.0] * 100, "times": [2.0] * 100})
                write_json(directory / "estimates.json", {"mean": estimate, "median": estimate})
                write_json(directory / "tukey.json", [1.0, 2.0, 3.0, 4.0])
    return criterion


def verify_saved(root: Path) -> subprocess.CompletedProcess[str]:
    """Exercise the installed verifier through the local saved-baseline recipe."""
    return run_command(
        "just",
        ["--no-deps", "--justfile", str(REPO_ROOT / "justfile"), "--working-directory", str(root), "bench-release-check", BASELINE, "target/criterion"],
        env={**os.environ, "UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
        check=False,
    )


def run_shell(script: str, cwd: Path, extra_env: dict[str, str], *, timeout: int = 30) -> subprocess.CompletedProcess[str]:
    bash = shutil.which("bash")
    assert bash is not None
    return subprocess.run(  # noqa: S603 - run only checked-in workflow shell in an isolated fixture.
        [bash, "--noprofile", "--norc", "-c", script],
        cwd=cwd,
        env={**{key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "SHELLOPTS", "BASHOPTS"}}, **extra_env},
        capture_output=True,
        encoding="utf-8",
        check=False,
        timeout=timeout,
    )


def test_complete_saved_baseline_passes_verification(dataset: Path) -> None:
    result = verify_saved(dataset.parent.parent)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"Verified {sum(map(len, IDS.values()))} complete saved Criterion benchmarks" in result.stdout


@pytest.mark.parametrize(
    "benchmark",
    [
        "d2/la_stack_dot",
        "column_basis/8x14/near_dependent",
        "interval_solve_d6/uncertain_dense",
        "polynomial/coincident_lifted/order8/determinant",
        "rational_dot/large/d22/slice",
    ],
)
def test_saved_baseline_must_match_fresh_samples(dataset: Path, benchmark: str) -> None:
    criterion = dataset
    write_json(criterion / STORAGE[benchmark] / BASELINE / "sample.json", {"iters": [1.0] * 100, "times": [3.0] * 100})
    result = verify_saved(criterion.parent.parent)
    assert result.returncode != 0
    assert "saved baseline differs from fresh measurement" in result.stderr


@pytest.mark.parametrize("relocation", ["moved", "swapped"])
def test_report_consumers_require_canonical_paths(dataset: Path, relocation: str) -> None:
    criterion = dataset
    original = criterion / "d2/la_stack_dot"
    moved = original.rename(criterion / "d2/misplaced")
    if relocation == "swapped":
        peer = criterion / "exact_d2/det_exact"
        peer.rename(original)
        moved.rename(peer)
    before = {path: path.read_bytes() for path in criterion.rglob("*") if path.is_file()}

    result = verify_saved(criterion.parent.parent)
    assert result.returncode != 0
    assert "required storage path differs" in result.stderr
    assert "d2/la_stack_dot" in result.stderr
    assert {path: path.read_bytes() for path in criterion.rglob("*") if path.is_file()} == before


def assert_scientific_inventory(root: Path) -> None:
    """Require the independent domain catalog while allowing diagnostic extras."""
    discovery = load_discovery(root, "tooling/criterion-inventory.toml")
    missing = REQUIRED - set(discovery.required)
    assert not missing, f"discovery omitted scientific cases: {sorted(missing)}"
    assert {suite.name for suite in discovery.suites} == {"public-default", "public-exact"}
    assert all(suite.command[-1] == "list" for suite in discovery.suites)
    saved = tomllib.loads((root / "tooling/criterion-saved.toml").read_text(encoding="utf-8"))
    missing = REQUIRED - set(saved["required"])
    assert not missing, f"saved inventory omitted scientific cases: {sorted(missing)}"
    assert set(saved["storage-paths"]) == set(saved["required"])
    assert saved["raw-files"] == list(RAW_FILES)


def test_configured_inventory_includes_public_workloads_without_peers() -> None:
    assert_scientific_inventory(REPO_ROOT)
    assert len(CATALOG["release"]["default"]) == 736
    assert len(CATALOG["release"]["exact"]) == 667
    assert len(REQUIRED) == 1403
    assert not any("/nalgebra_" in name or "/faer_" in name for name in REQUIRED)


@pytest.mark.parametrize("declaration", ["criterion-inventory.toml", "criterion-saved.toml", "both"])
@pytest.mark.parametrize(
    "removed",
    [
        "d3/la_stack_norm2_sq",
        "public_d4/interval_det",
        "canonical_conversion_dyadic_d3/strict_array",
        "column_basis/8x14/near_dependent",
        "interval_solve_d6/uncertain_dense",
        "polynomial/coincident_lifted/order8/determinant",
        "rational_dot/large/d22/slice",
    ],
)
def test_scientific_inventory_detects_a_nonsentinel_omission(tmp_path: Path, declaration: str, removed: str) -> None:
    tooling = tmp_path / "tooling"
    tooling.mkdir()
    for filename in ("criterion-inventory.toml", "criterion-saved.toml"):
        contents = (REPO_ROOT / "tooling" / filename).read_text(encoding="utf-8")
        if declaration in {filename, "both"}:
            assert removed in tomllib.loads(contents)["required"]
            contents = contents.replace(f'    "{removed}",\n', "").replace(f'"{removed}" = "{removed}"\n', "")
            assert removed not in tomllib.loads(contents)["required"]
        (tooling / filename).write_text(contents, encoding="utf-8", newline="\n")
    with pytest.raises(AssertionError, match=f"omitted scientific cases:.*{re.escape(removed)}"):
        assert_scientific_inventory(tmp_path)


def test_workflow_reuses_original_acquisition_and_isolates_write_credentials() -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    assert set(workflow[True]) == {"workflow_dispatch"}
    assert set(workflow[True]["workflow_dispatch"]["inputs"]) == {"tag", "candidate"}
    assert workflow["concurrency"]["cancel-in-progress"] is False
    jobs = workflow["jobs"]
    assert set(jobs) == {"fetch", "verify", "attach"}
    assert jobs["verify"]["needs"] == "fetch"
    assert jobs["attach"]["needs"] == "verify"
    assert jobs["verify"]["permissions"] == {"contents": "read", "checks": "read", "actions": "read"}
    assert "github.triggering_actor == 'acgetchell'" in jobs["fetch"]["if"]
    for name in ("fetch", "attach"):
        job = jobs[name]
        assert job["permissions"]["contents"] == "write"
        assert not any("checkout@" in step.get("uses", "") or step.get("uses", "").startswith("./") for step in job["steps"])
        assert job["steps"][0]["with"]["enable-cache"] is False
        assert job["steps"][0]["with"]["cache-python"] is False
    verifier = jobs["verify"]["steps"]
    assert verifier[0]["with"] == {"ref": "${{ github.sha }}", "fetch-depth": 0, "persist-credentials": False}
    assert verifier[1]["with"] == {"cache": False, "tool": "cargo"}
    assert jobs["verify"]["timeout-minutes"] == 60
    assert verifier[3]["with"]["artifact-ids"] == "${{ needs.fetch.outputs.artifact-id }}"
    assert verifier[4]["run"] == "just release-dependency-builds"
    assert verifier[5]["env"]["CANDIDATE_TRANSFER"] == "${{ inputs.candidate }}"
    assert verifier[5]["run"] == 'just release-receive "$RUNNER_TEMP/accepted-candidate" --incoming "$RUNNER_TEMP/incoming-candidate"'
    attach = jobs["attach"]["steps"]
    assert attach[2]["with"]["artifact-ids"] == "${{ needs.verify.outputs.artifact-id }}"
    assert attach[1]["env"]["ARTIFACT_DIGEST"] == "${{ needs.verify.outputs.artifact-digest }}"
    text = WORKFLOW.read_text(encoding="utf-8")
    runtime_requirements = [shlex.split(line.removesuffix("\\").strip()) for line in text.splitlines() if line.lstrip().startswith("--from ")]
    assert runtime_requirements == [["--from", TOOLING_REQUIREMENT]] * 2
    for forbidden in ("bench-save-baseline", "measure-baseline", "cargo bench", "performance measure", "--clobber", "continue-on-error"):
        assert forbidden not in text
    assert "publish_release_asset(root, target.repository, target.tag, artifact, target=target)" in text
    assert "publish=True" not in text


@pytest.mark.parametrize("validated_attempt", ["", "1", "2"])
def test_retry_requires_current_retrieval(tmp_path: Path, validated_attempt: str) -> None:
    workflow = yaml.safe_load(WORKFLOW.read_bytes())
    result = run_shell(
        workflow["jobs"]["verify"]["steps"][2]["run"],
        tmp_path,
        {"VALIDATED_ATTEMPT": validated_attempt, "GITHUB_RUN_ATTEMPT": "2"},
    )
    assert result.returncode == (0 if validated_attempt == "2" else 1)
