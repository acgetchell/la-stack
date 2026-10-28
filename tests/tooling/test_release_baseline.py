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
from urllib.parse import urlsplit

import pytest
from packaging.requirements import Requirement
from research_repo_tools.criterion_inventory import CriterionInventory, load_discovery, serialize_inventory
from research_repo_tools.host_metadata import HostMetadata, serialize_host
from research_repo_tools.process import run_command

REPO_ROOT = Path(__file__).resolve().parents[2]
SHARED_REQUIREMENT = Requirement(tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["dependency-groups"]["tooling"][0])
# Release orchestration needs the core package independently of contributor extras.
TOOLING_REQUIREMENT = SHARED_REQUIREMENT.name + str(SHARED_REQUIREMENT.specifier)
WORKFLOW = REPO_ROOT / ".github/workflows/release-benchmarks.yml"
PREPARATION = REPO_ROOT / ".github/actions/prepare-release-benchmarks/action.yml"
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
SUITE_ENV = {"RELEASE_TAG": BASELINE, "GITHUB_SHA": "a" * 40, "GITHUB_REPOSITORY": "owner/repo", "GITHUB_RUN_ID": "123", "GITHUB_RUN_ATTEMPT": "2"}


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) + "\n", encoding="utf-8", newline="\n")


@pytest.fixture
def dataset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
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
    return criterion, manifest


def verify_saved(root: Path, *, suite: str | None = None, criterion_dir: str = "target/criterion") -> subprocess.CompletedProcess[str]:
    """Exercise the installed verifier through the actual release-check recipe."""
    recipe = ["bench-release-check", BASELINE, criterion_dir] if suite is None else ["bench-release-check-suite", BASELINE, suite, criterion_dir]
    return run_command(
        "just",
        ["--no-deps", "--justfile", str(REPO_ROOT / "justfile"), "--working-directory", str(root), *recipe],
        env={**os.environ, "UV_NO_SYNC": "1", "UV_PROJECT_ENVIRONMENT": sys.prefix},
        check=False,
    )


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


def installed_runner() -> str:
    """Keep workflow package assertions while using the installed locked CLI."""
    command = f"{shlex.quote(sys.executable)} -I -X utf8 -m research_repo_tools"
    return f'uvx() {{ [[ "$1 $2 $3 $4 $5" == "--python 3.15 --from {TOOLING_REQUIREMENT} research-repo-tools" ]] || return 90; shift 5; {command} "$@"; }}\n'


def package_suites(root: Path, manifest: Path, criterion: Path) -> None:
    """Execute the workflow's installed CLI with all independently required cases."""
    destination = root / "target/release-suites"
    destination.mkdir(parents=True)
    for suite, ids in IDS.items():
        shard = root / "shards" / suite
        output = shard / "target/criterion"
        output.mkdir(parents=True)
        (shard / "tooling").mkdir()
        shutil.copyfile(REPO_ROOT / "tooling/criterion-saved.toml", shard / "tooling/criterion-saved.toml")
        shutil.copyfile(manifest, shard / "target/release-benchmark-inventory.json")
        host = HostMetadata("fixture OS", "fixture architecture", f"CPU for {suite}", 2, 4, None)
        (shard / "target/release-host.json").write_bytes(serialize_host(host))
        for event, timestamp in (("started", "100"), ("finished", "200")):
            (shard / f"target/{suite}-{event}").write_text(timestamp, encoding="utf-8", newline="\n")
        for benchmark in ids:
            storage = STORAGE[benchmark]
            shutil.copytree(criterion / storage, output / storage)
        # Substitute only native version/hash observations. The real uvx command
        # installs the exact published suite package independently of root tooling.
        probes = (
            f'rustc() {{ printf "rustc fixture\\nLLVM fixture\\n"; }}\ncargo() {{ echo "cargo fixture"; }}\nsha256sum() {{ echo "{"f" * 64}  Cargo.lock"; }}\n'
        )
        result = run_shell(
            probes + shell_step("Package Criterion suite"),
            shard,
            {**SUITE_ENV, "BENCH_SUITE": suite, "RUNNER_NAME": suite, "ImageOS": "ubuntu24", "ImageVersion": "fixture-image"},
            timeout=180,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        shutil.move(shard / "target/release-suites" / suite, destination / suite)


def test_complete_dataset_round_trips_through_workflow_archive(dataset: tuple[Path, Path], tmp_path: Path) -> None:
    criterion, manifest = dataset
    criterion = criterion.rename(tmp_path / "original-criterion")
    package_suites(tmp_path, manifest, criterion)
    for suite, ids in IDS.items():
        verified = verify_saved(tmp_path, suite=suite, criterion_dir=f"shards/{suite}/target/criterion")
        assert verified.returncode == 0, verified.stdout + verified.stderr
        assert f"Verified {len(ids)} complete saved Criterion benchmarks" in verified.stdout
    # Assembly must use the downloaded inventory, not an accidentally retained
    # fixture copy or an inventory inferred from the supplied archives.
    manifest.unlink()
    output = tmp_path / "step outputs"
    asset = f"la-stack-{BASELINE}-criterion-baseline.tar.gz"
    script = shell_step("Assemble release Criterion baseline")
    env = {**SUITE_ENV, "GITHUB_WORKSPACE": tmp_path.as_posix(), "GITHUB_OUTPUT": output.as_posix()}
    for assertion in ({"GITHUB_SHA": "b" * 40}, {"GITHUB_RUN_ATTEMPT": "3"}):
        rejected = run_shell(script, tmp_path, {**env, **assertion}, timeout=180)
        assert rejected.returncode != 0
        assert "source revision or collection session differs" in rejected.stderr
        assert not (tmp_path / asset).exists()
        assert not output.exists()
    # Git Bash accepts forward-slash drive paths on Windows as well as POSIX paths.
    result = run_shell(script, tmp_path, env, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assert output.read_text(encoding="utf-8") == f"asset={asset}\nartifact-name=bench-baseline-{BASELINE}-123-2\n"
    assert_archive_contents(tmp_path / asset, criterion)
    # Release readers keep the standard tooling environment. Exercise its
    # documented extraction and full-dataset verification against the new archive.
    run_command(
        "uv",
        ["run", "--locked", "--no-sync", "research-repo-tools", "performance", "extract", str(tmp_path / asset), str(tmp_path / "restored-criterion")],
        cwd=REPO_ROOT,
    )
    verified = verify_saved(tmp_path, criterion_dir="restored-criterion")
    assert verified.returncode == 0, verified.stdout + verified.stderr
    assert f"Verified {sum(map(len, IDS.values()))} complete saved Criterion benchmarks" in verified.stdout


def assert_archive_contents(archive_path: Path, criterion: Path) -> None:
    """Require every raw byte and original runner observation in the final asset."""
    with tarfile.open(archive_path) as archive:
        for ids in IDS.values():
            for benchmark in ids:
                for sample in ("new", BASELINE):
                    for filename in RAW_FILES:
                        member = archive.extractfile(f"{STORAGE.get(benchmark, benchmark)}/{sample}/{filename}")
                        assert member is not None
                        assert member.read() == (criterion / STORAGE.get(benchmark, benchmark) / sample / filename).read_bytes()
        member = archive.extractfile("criterion-inventory.json")
        assert member is not None
        assert json.loads(member.read())["suites"] == IDS
        for index, suite in enumerate(sorted(IDS)):
            member = archive.extractfile(f"release-suites/{index}.json")
            assert member is not None
            retained = member.read()
            record = json.loads(retained)
            assert record["suite"] == suite
            run = record["run"]
            assert run["provenance"]["revision"] == SUITE_ENV["GITHUB_SHA"]
            assert run["session"] == "owner/repo:123:2"
            assert (run["started_at"], run["finished_at"]) == (100, 200)
            assert run["host"]["cpu"] == f"CPU for {suite}"
            context = run["provenance"]["context"]
            assert context["command"] == f"just bench-save-baseline {BASELINE} {suite}"
            assert context["features"] == ("bench,exact" if suite == "public-exact" else "bench")
            assert context["runner"] == suite
            assert context["rustc"] == "rustc fixture;LLVM fixture"
            assert context["cargo"] == "cargo fixture"
            assert context["cargo_lock_sha256"] == "f" * 64
            assert (context["image"], context["image_version"]) == ("ubuntu24", "fixture-image")
            with tarfile.open(archive_path.parent / "target/release-suites" / suite / "criterion.tar.gz") as original:
                original_manifest = original.extractfile("suite.json")
                assert original_manifest is not None
                assert retained == original_manifest.read()


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
def test_saved_baseline_must_match_fresh_samples(dataset: tuple[Path, Path], benchmark: str) -> None:
    criterion, _manifest = dataset
    write_json(criterion / STORAGE[benchmark] / BASELINE / "sample.json", {"iters": [1.0] * 100, "times": [3.0] * 100})
    result = verify_saved(criterion.parent.parent)
    assert result.returncode != 0
    assert "saved baseline differs from fresh measurement" in result.stderr


@pytest.mark.parametrize("relocation", ["moved", "swapped"])
def test_report_consumers_require_canonical_paths(dataset: tuple[Path, Path], relocation: str) -> None:
    criterion, _manifest = dataset
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
        "Capture runner host",
        "Save Criterion suite",
        "Package Criterion suite",
        "Upload temporary suite artifact",
        "Assemble release Criterion baseline",
        "Upload temporary baseline artifact",
    ]
    positions = [workflow.index(f"- name: {step}") for step in steps]
    assert positions == sorted(positions)
    for step in steps:
        body = workflow.split(f"      - name: {step}\n", 1)[1].split("\n      - ", 1)[0]
        assert "\n        if:" not in body
    assert "continue-on-error" not in workflow
    assert "workflow_dispatch:" in workflow
    assert "required: true" in workflow
    assert "\n  release:" not in workflow
    assert "group: release-benchmarks-${{ inputs.tag }}" in workflow
    assert "cancel-in-progress: false" in workflow
    preflight, remainder = workflow.split("  measure-baseline:\n", 1)
    producer, remainder = remainder.split("  release-baseline:\n", 1)
    _assembler, publisher = remainder.split("  publish-baseline:\n", 1)
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
    assert "performance host --output target/release-host.json" in producer
    assert "--host target/release-host.json" in producer
    assert "path: target/release-suites/${{ matrix.suite }}" in producer
    assert "tar -" not in workflow
    assert "tooling.release_benchmarks" not in (REPO_ROOT / "justfile").read_text(encoding="utf-8")


def test_parallel_suites_require_both_current_attempts_before_assembly() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")
    producer, remainder = workflow.split("  measure-baseline:\n", 1)[1].split("  release-baseline:\n", 1)
    assembler = remainder.split("  publish-baseline:\n", 1)[0]
    assert "fail-fast: false" in producer
    assert "suite: [public-default, public-exact]" in producer
    assert "max-parallel: 1" not in producer
    assert "BENCH_SUITE: ${{ matrix.suite }}" in producer
    assert 'just bench-save-baseline "$RELEASE_TAG" "$BENCH_SUITE"' in producer
    assert "performance pack-suite" in producer
    assert 'tooling/criterion-saved.toml "$BENCH_SUITE" --saved-sample "$RELEASE_TAG"' in producer
    assert "ref: ${{ github.sha }}" in assembler
    assert "persist-credentials: false" in assembler
    assert "contents: write" not in assembler
    assert "needs: [validate-release, measure-baseline]" in assembler
    assert "always()" not in assembler
    for suite in IDS:
        assert f"name: bench-suite-${{{{ inputs.tag }}}}-${{{{ github.run_id }}}}-${{{{ github.run_attempt }}}}-{suite}" in assembler
        assert f"path: target/release-suites/{suite}" in assembler
        assert f"target/release-suites/{suite}/criterion.tar.gz" in assembler
    assert "performance assemble-suites tooling/criterion-saved.toml" in assembler
    for job in (producer, assembler):
        assert '--session "$GITHUB_REPOSITORY:$GITHUB_RUN_ID:$GITHUB_RUN_ATTEMPT"' in job
        assert '--revision "$GITHUB_SHA"' in job
    for job in (producer, assembler):
        assert f"--from {TOOLING_REQUIREMENT}" in job


def test_setup_limits_preserve_benchmark_and_tail_budgets() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")
    producer = workflow.split("  measure-baseline:\n", 1)[1].split("  release-baseline:\n", 1)[0]
    job_limit = re.search(r"(?m)^    timeout-minutes: (\d+)$", producer)
    assert job_limit is not None
    setup = producer.split("    steps:\n", 1)[1].split("      - name: Save Criterion suite", 1)[0]
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
    preparation = re.search(r"(?ms)^      - name: Prepare release benchmarks\n.*?^        timeout-minutes: (\d+)$", setup)
    assert preparation is not None
    assert int(preparation[1]) >= 35 + discovery_minutes + 5
    assert "uses: ./.github/actions/prepare-release-benchmarks" in setup
    assert "using: composite" in PREPARATION.read_text(encoding="utf-8")
    benchmark_limit = re.search(r"(?ms)^      - name: Save Criterion suite\n.*?^        timeout-minutes: (\d+)$", producer)
    assert benchmark_limit is not None
    assert int(benchmark_limit[1]) == 240
    for suite in ("default", "exact"):
        # Keep at least 50% headroom beyond the documented planning rate.
        assert int(benchmark_limit[1]) * 60 >= len(CATALOG["release"][suite]) * 12 * 1.5
    assert int(job_limit[1]) - sum(setup_limits) - int(benchmark_limit[1]) >= 20
    assert int(job_limit[1]) <= 360


@pytest.mark.parametrize("suite", ["public-default", "public-exact"])
def test_failed_suite_retains_timing_and_does_not_mark_completion(tmp_path: Path, suite: str) -> None:
    (tmp_path / "target").mkdir()
    # A shell function substitutes only the benchmark command; timestamps and
    # failure handling execute exactly as checked into the workflow.
    script = "just() { return 7; }\n" + shell_step("Save Criterion suite")
    result = run_shell(script, tmp_path, {"RELEASE_TAG": BASELINE, "BENCH_SUITE": suite})
    assert result.returncode == 7
    assert (tmp_path / f"target/{suite}-started").is_file()
    assert not (tmp_path / f"target/{suite}-finished").exists()
    summary = tmp_path / "suite summary.md"
    result = run_shell(
        shell_step("Suite timing summary"),
        tmp_path,
        {
            "GITHUB_STEP_SUMMARY": summary.as_posix(),
            "BENCH_SUITE": suite,
            "SUITE_OUTCOME": "failure",
        },
    )
    assert result.returncode == 0, result.stderr
    assert f"| {suite} | failure |" in summary.read_text(encoding="utf-8")
    assert "| 240 |" in summary.read_text(encoding="utf-8")


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
from urllib.parse import urlsplit
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
elif any(urlsplit(item).hostname == "uploads.github.com" for item in args):
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
    runner = installed_runner()
    result = run_shell(runner + shell_step("Validate draft release target"), fixture, env)
    assert result.returncode == 0, result.stderr
    target = output.read_text(encoding="utf-8").strip().removeprefix("target=")
    assert json.loads(target)["release_id"] == 42
    assert json.loads(target)["commit"] == "a" * 40
    assert json.loads(target)["expected_title"] == BASELINE
    stale = run_shell(runner + shell_step("Attach baseline for review"), fixture, {**env, "RELEASE_TARGET": target, "REMOTE_COMMIT": "b" * 40})
    assert stale.returncode != 0
    assert not (fixture / "uploaded").exists()
    assert not (fixture / "published").exists()
    stale_calls = [json.loads(line) for line in (fixture / "api-calls").read_text(encoding="utf-8").splitlines()]
    assert not any(urlsplit(item).hostname == "uploads.github.com" for call in stale_calls for item in call)
    result = run_shell(runner + shell_step("Attach baseline for review"), fixture, {**env, "RELEASE_TARGET": target})
    assert result.returncode == 0, result.stderr
    assert (fixture / "uploaded").exists()
    assert not (fixture / "published").exists()
    calls = [json.loads(line) for line in (fixture / "api-calls").read_text(encoding="utf-8").splitlines()]
    assert not any("PATCH" in call for call in calls)
    assert any("/releases/assets/7" in item for call in calls for item in call)
