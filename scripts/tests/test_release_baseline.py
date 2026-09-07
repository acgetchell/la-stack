"""Exercise complete release archives and the workflow's publication barriers."""

import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import textwrap
from pathlib import Path

import pytest

import release_baseline

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github/workflows/release-benchmarks.yml"
PREPARATION = REPO_ROOT / ".github/actions/prepare-release-benchmarks/action.yml"
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
    criterion, manifest = dataset
    assert release_baseline.validate(criterion, manifest, BASELINE) == 3
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
            estimates = json.loads((directory / "estimates.json").read_text(encoding="utf-8"))
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


def test_workflow_requires_preflight_and_isolates_publication_permissions() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")
    preparation = PREPARATION.read_text(encoding="utf-8")
    assert preparation.index("- name: Validate benchmark inputs") < preparation.index("- name: Inventory full release suites")
    assert "run: just test-bench-inputs" in preparation
    assert "run: just bench-release-inventory" in preparation
    assert "continue-on-error" not in preparation
    assert preparation.index("- name: Require fresh release preflight") < preparation.index("- name: Install Rust toolchain")
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
    assert sum(setup_limits) <= 30
    assert "uses: ./.github/actions/prepare-release-benchmarks" in setup
    assert "using: composite" in PREPARATION.read_text(encoding="utf-8")
    benchmark_limits = [
        int(limit) for limit in re.findall(r"(?ms)^      - name: Save (?:comparative|exact) Criterion baseline\n.*?^        timeout-minutes: (\d+)$", producer)
    ]
    assert benchmark_limits == [150, 90]
    assert int(job_limit[1]) - sum(setup_limits) - sum(benchmark_limits) >= 15


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


# Only the GitHub boundary is substituted. jq predicates, checksums, shell
# failure handling, and publication ordering execute from the workflow itself.
GITHUB_STUB = r"""
gh() {
  printf '%s\n' "$*" >> api-calls
  [[ "$1" == api ]] || return 90
  case "$2" in
    'repos/owner/repo/releases?per_page=100')
      [[ "$*" == *'--paginate --slurp' ]] || return 91
      cat releases.json
      ;;
    "repos/owner/repo/commits/refs/tags/$RELEASE_TAG")
      [[ "${COMMIT_STATUS:-0}" == 0 ]] || return "$COMMIT_STATUS"
      if [[ -f commit-read ]]; then
        printf '%s\n' "${LATER_COMMIT:-$RELEASE_COMMIT}"
      else
        touch commit-read
        printf '%s\n' "$RELEASE_COMMIT"
      fi
      ;;
    repos/owner/repo/releases/42)
      if [[ "$*" == *'--method PATCH -F draft=false' ]]; then
        touch publish-called
        [[ "${PUBLISH_STATUS:-0}" == 0 ]] || return "$PUBLISH_STATUS"
        cat published.json
      elif [[ -f draft-read ]]; then
        cat later-draft.json
      else
        touch draft-read
        cat draft.json
      fi
      ;;
    'repos/owner/repo/releases/42/assets?per_page=100')
      [[ "$*" == *'--paginate --slurp' ]] || return 92
      if [[ -f upload-called ]]; then
        cat uploaded-assets.json
      else
        cat assets.json
      fi
      ;;
    "https://uploads.github.com/repos/owner/repo/releases/42/assets?name=$RELEASE_ASSET")
      [[ "$*" == *'--method POST --header Content-Type: application/gzip --input '* ]] || return 93
      touch upload-called
      return "${UPLOAD_STATUS:-0}"
      ;;
    *) return 94 ;;
  esac
}
"""


def draft_release(**changes: object) -> dict[str, object]:
    return {"id": 42, "tag_name": BASELINE, "name": BASELINE, "draft": True, "prerelease": False, "immutable": False, **changes}


@pytest.fixture
def github_release(tmp_path: Path) -> dict[str, str]:
    asset = f"la-stack-{BASELINE}-criterion-baseline.tar.gz"
    payload = b"test archive bytes\n"
    (tmp_path / asset).write_bytes(payload)
    metadata = {"name": asset, "state": "uploaded", "size": len(payload), "digest": f"sha256:{hashlib.sha256(payload).hexdigest()}"}
    write_json(tmp_path / "releases.json", [[draft_release(tag_name="v0.1.0")], [draft_release()]])
    write_json(tmp_path / "draft.json", draft_release())
    write_json(tmp_path / "later-draft.json", draft_release())
    write_json(tmp_path / "published.json", draft_release(draft=False, immutable=True))
    write_json(tmp_path / "assets.json", [[]])
    write_json(tmp_path / "uploaded-assets.json", [[{"name": "unrelated.txt"}], [metadata]])
    return {
        "GH_REPO": "owner/repo",
        "RELEASE_TAG": BASELINE,
        "RELEASE_ID": "42",
        "RELEASE_COMMIT": "a" * 40,
        "RELEASE_ASSET": asset,
        "GITHUB_REF": f"refs/tags/{BASELINE}",
        "GITHUB_SHA": "a" * 40,
        "GITHUB_OUTPUT": (tmp_path / "outputs").as_posix(),
    }


def test_preflight_captures_draft_identity_and_peeled_tag_commit(tmp_path: Path, github_release: dict[str, str]) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, github_release)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "outputs").read_text(encoding="utf-8") == f"release-id=42\ncommit={'a' * 40}\n"
    calls = (tmp_path / "api-calls").read_text(encoding="utf-8")
    assert "--paginate --slurp" in calls
    assert f"commits/refs/tags/{BASELINE}" in calls
    assert "--method" not in calls


@pytest.mark.parametrize("ref", ["refs/heads/main", f"refs/heads/{BASELINE}", "refs/tags/v9.9.9", ""])
def test_preflight_rejects_dispatch_outside_release_tag(tmp_path: Path, github_release: dict[str, str], ref: str) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, {**github_release, "GITHUB_REF": ref})
    assert result.returncode == 1
    assert f"Dispatch with --ref {BASELINE}" in result.stdout
    assert not (tmp_path / "api-calls").exists()
    assert not (tmp_path / "outputs").exists()


def test_preflight_rejects_tag_moved_since_dispatch(tmp_path: Path, github_release: dict[str, str]) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, {**github_release, "RELEASE_COMMIT": "b" * 40})
    assert result.returncode == 1
    assert "Release tag no longer matches the workflow commit" in result.stdout
    assert not (tmp_path / "outputs").exists()


@pytest.mark.parametrize("validated_attempt", ["", "1", "2"])
def test_benchmark_retry_requires_current_preflight(tmp_path: Path, validated_attempt: str) -> None:
    result = run_shell(
        shell_step("Require fresh release preflight", PREPARATION),
        tmp_path,
        {"VALIDATED_ATTEMPT": validated_attempt, "GITHUB_RUN_ATTEMPT": "2"},
    )
    assert result.returncode == (0 if validated_attempt == "2" else 1)
    if validated_attempt != "2":
        assert "Rerun all jobs or dispatch again" in result.stdout


@pytest.mark.parametrize("tag", ["", "1.2.3", "v01.2.3", "v1.02.3", "v1.2.03", "v1.2.3-rc.1", "v1.2.3+meta", "main", "v1.2.3\nextra", "$(touch injected)"])
def test_preflight_rejects_invalid_tag_before_api_access(tmp_path: Path, github_release: dict[str, str], tag: str) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, {**github_release, "RELEASE_TAG": tag})
    assert result.returncode == 1
    assert "Expected an existing stable vX.Y.Z tag" in result.stdout
    assert not (tmp_path / "api-calls").exists()
    assert not (tmp_path / "outputs").exists()
    assert not (tmp_path / "injected").exists()


@pytest.mark.parametrize(
    "releases",
    [
        [],
        [draft_release(), draft_release(id=43)],
        [draft_release(draft=False)],
        [draft_release(prerelease=True)],
        [draft_release(immutable=True)],
        [draft_release(immutable=None)],
        [draft_release(tag_name="v9.9.9")],
        [draft_release(name="Wrong title")],
        [draft_release(id="42\ncommit=injected")],
    ],
)
def test_preflight_rejects_unsuitable_release(tmp_path: Path, github_release: dict[str, str], releases: list[dict[str, object]]) -> None:
    write_json(tmp_path / "releases.json", [releases])
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, github_release)
    assert result.returncode != 0
    assert "commits/" not in (tmp_path / "api-calls").read_text(encoding="utf-8")
    assert not (tmp_path / "outputs").exists()


@pytest.mark.parametrize("response", ["{invalid", "null"])
def test_preflight_rejects_malformed_api_response(tmp_path: Path, github_release: dict[str, str], response: str) -> None:
    (tmp_path / "releases.json").write_text(response, encoding="utf-8")
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, github_release)
    assert result.returncode != 0
    assert not (tmp_path / "outputs").exists()


@pytest.mark.parametrize("overrides", [{"COMMIT_STATUS": "4"}, {"RELEASE_COMMIT": ""}, {"RELEASE_COMMIT": "a" * 40 + "\ninjected=true"}])
def test_missing_or_invalid_tag_commit_stops_preflight(tmp_path: Path, github_release: dict[str, str], overrides: dict[str, str]) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Validate draft release target"), tmp_path, {**github_release, **overrides})
    assert result.returncode != 0
    assert not (tmp_path / "outputs").exists()


@pytest.mark.parametrize("reuse", [False, True], ids=["fresh-upload", "interrupted-upload-rerun"])
def test_publisher_verifies_durable_asset_before_publication(tmp_path: Path, github_release: dict[str, str], reuse: bool) -> None:
    if reuse:
        shutil.copyfile(tmp_path / "uploaded-assets.json", tmp_path / "assets.json")
    result = run_shell(GITHUB_STUB + shell_step("Attach baseline and publish draft"), tmp_path, github_release)
    assert result.returncode == 0, result.stderr
    calls = (tmp_path / "api-calls").read_text(encoding="utf-8").splitlines()
    assert calls[-1] == "api repos/owner/repo/releases/42 --method PATCH -F draft=false"
    assert calls[-2] == "api repos/owner/repo/releases/42/assets?per_page=100 --paginate --slurp"
    assert calls.count("api repos/owner/repo/releases/42") == 2
    assert (tmp_path / "upload-called").exists() is not reuse
    assert (tmp_path / "publish-called").exists()
    assert "--clobber" not in " ".join(calls)


@pytest.mark.parametrize("snapshot", ["draft.json", "later-draft.json"])
@pytest.mark.parametrize(
    "changes",
    [{"draft": False}, {"prerelease": True}, {"immutable": True}, {"tag_name": "v9.9.9"}, {"name": "renamed"}, {"id": 43}],
)
def test_changed_release_stops_publication(tmp_path: Path, github_release: dict[str, str], snapshot: str, changes: dict[str, object]) -> None:
    write_json(tmp_path / snapshot, draft_release(**changes))
    result = run_shell(GITHUB_STUB + shell_step("Attach baseline and publish draft"), tmp_path, github_release)
    assert result.returncode != 0
    assert not (tmp_path / "publish-called").exists()
    assert (tmp_path / "upload-called").exists() is (snapshot == "later-draft.json")


@pytest.mark.parametrize("after_upload", [False, True])
def test_moved_tag_stops_publication(tmp_path: Path, github_release: dict[str, str], after_upload: bool) -> None:
    overrides = {"LATER_COMMIT": "b" * 40}
    if not after_upload:
        (tmp_path / "commit-read").touch()
    result = run_shell(GITHUB_STUB + shell_step("Attach baseline and publish draft"), tmp_path, {**github_release, **overrides})
    assert result.returncode != 0
    assert not (tmp_path / "publish-called").exists()
    assert (tmp_path / "upload-called").exists() is after_upload


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("problem", ["missing", "duplicate", "wrong-name", "wrong-digest", "wrong-size", "starter", "no-digest"])
def test_bad_durable_asset_stops_publication(tmp_path: Path, github_release: dict[str, str], existing: bool, problem: str) -> None:
    metadata = json.loads((tmp_path / "uploaded-assets.json").read_text(encoding="utf-8"))[1][0]
    assets = [metadata]
    if problem == "missing":
        assets = []
    elif problem == "duplicate":
        assets.append(metadata)
    else:
        field, value = {
            "wrong-name": ("name", "wrong.tar.gz"),
            "wrong-digest": ("digest", "sha256:" + "0" * 64),
            "wrong-size": ("size", 999),
            "starter": ("state", "starter"),
            "no-digest": ("digest", None),
        }[problem]
        metadata[field] = value
    write_json(tmp_path / "uploaded-assets.json", [assets])
    if existing:
        write_json(tmp_path / "assets.json", [assets])
    result = run_shell(GITHUB_STUB + shell_step("Attach baseline and publish draft"), tmp_path, github_release)
    assert result.returncode != 0
    assert not (tmp_path / "publish-called").exists()
    assert (tmp_path / "upload-called").exists() is (not existing or problem in {"missing", "wrong-name"})


def test_upload_failure_leaves_draft_unpublished(tmp_path: Path, github_release: dict[str, str]) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Attach baseline and publish draft"), tmp_path, {**github_release, "UPLOAD_STATUS": "8"})
    assert result.returncode == 8
    assert (tmp_path / "upload-called").exists()
    assert not (tmp_path / "publish-called").exists()


def test_publication_failure_propagates_without_deleting_uploaded_asset(tmp_path: Path, github_release: dict[str, str]) -> None:
    result = run_shell(GITHUB_STUB + shell_step("Attach baseline and publish draft"), tmp_path, {**github_release, "PUBLISH_STATUS": "9"})
    assert result.returncode != 0
    assert (tmp_path / "upload-called").exists()
    assert (tmp_path / "publish-called").exists()
    assert "DELETE" not in (tmp_path / "api-calls").read_text(encoding="utf-8")
