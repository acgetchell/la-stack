"""Tests for archive_performance.py."""

from __future__ import annotations

import io
import json
import subprocess
import tarfile
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import pytest

import archive_performance
from archive_performance import (
    GenerationConfig,
    generate_and_promote_worktree_report,
    generate_worktree_report,
    main,
    normalize_tag,
    parse_report_id,
    promote_report,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

type RunnerCall = tuple[str, tuple[str, ...], Path | None]


def _result(stdout: str = "") -> SimpleNamespace:
    return SimpleNamespace(stdout=stdout)


def _git(repo_root: Path, *args: str) -> str:
    return archive_performance.run_git_command(list(args), cwd=repo_root).stdout


def _report(version: str, baseline: str) -> str:
    return (
        "# Benchmark Performance\n\n"
        f"**la-stack** v{version} · `abc1234` (release/test) · 2026-06-08 12:00:00 UTC\n"
        "**Statistic**: median\n"
        "**Suite**: all\n"
        "**Scope**: release-signal\n\n"
        "## Benchmark Results\n\n"
        f"Comparison against baseline **{baseline}**:\n\n"
        "Negative change = faster. Speedup > 1.00x = improvement.\n\n"
        "## Exact arithmetic\n\n"
        "| Benchmark | Baseline | Latest | Change | Speedup |\n"
        "|-----------|---------:|-------:|-------:|--------:|\n"
        "| det_exact | 1.0 ns | 0.9 ns | -10.0% | 1.11x |\n"
    )


def _normalized_report(version: str, baseline: str) -> str:
    return archive_performance._normalize_how_to_update(_report(version, baseline))


def test_normalized_report_links_archived_performance_reports() -> None:
    text = _normalized_report("0.4.3", "v0.4.2")
    assert "Older curated release-to-release reports are archived in `docs/archive/performance/`." in text


def _legacy_report(version: str, baseline: str) -> str:
    return (
        _report(version, baseline)
        + "\n"
        + "## How to Update\n\n"
        + "```bash\n"
        + "git checkout v0.2.0\n"
        + "just bench-save-baseline v0.2.0\n"
        + "git checkout main\n"
        + "just bench-compare v0.2.0\n"
        + "```\n"
    )


def _write_baseline_archive(path: Path, *, include_harness_metadata: bool = False) -> None:
    tag = path.name.removeprefix("la-stack-").removesuffix("-criterion-baseline.tar.gz")
    fixture_dir = path.parent / f"baseline-fixture-{tag}"
    criterion_dir = fixture_dir / "criterion"
    criterion_dir.mkdir(parents=True)
    (criterion_dir / "placeholder.txt").write_text("baseline\n", encoding="utf-8")
    sample_dir = criterion_dir / "exact_d2" / "det_exact" / tag
    sample_dir.mkdir(parents=True)
    (sample_dir / "estimates.json").write_text('{"median":{"point_estimate":1.0}}\n', encoding="utf-8")
    if include_harness_metadata:
        (criterion_dir / archive_performance._BENCHMARK_HARNESS_METADATA).write_text(
            json.dumps(
                {
                    "baseline": tag,
                    "mode": "shared-current-harness",
                    "schema": 1,
                    "sha256": "a" * 64,
                }
            ),
            encoding="utf-8",
        )
    with tarfile.open(path, "w:gz") as tar:
        tar.add(criterion_dir, arcname="criterion")


def _write_unsafe_baseline_archive(path: Path) -> None:
    payload = b"escape\n"
    info = tarfile.TarInfo("../escape.txt")
    info.size = len(payload)
    with tarfile.open(path, "w:gz") as tar:
        tar.addfile(info, io.BytesIO(payload))


def _write_current_benchmark_tooling(worktree: Path) -> None:
    (worktree / ".config").mkdir(parents=True, exist_ok=True)
    (worktree / "scripts").mkdir(parents=True, exist_ok=True)
    (worktree / "benches" / "common").mkdir(parents=True, exist_ok=True)
    (worktree / "src").mkdir(parents=True, exist_ok=True)
    (worktree / "tests").mkdir(parents=True, exist_ok=True)
    (worktree / "benches" / "vs_linalg.rs").write_text("fn main() {}\n", encoding="utf-8")
    (worktree / "benches" / "common" / "inputs.rs").write_text("pub const INPUT: f64 = 1.0;\n", encoding="utf-8")
    (worktree / "src" / "lib.rs").write_text("pub fn fixture() {}\n", encoding="utf-8")
    (worktree / "tests" / "exact_bench_config.rs").write_text("// exact fixture\n", encoding="utf-8")
    (worktree / "tests" / "vs_linalg_inputs.rs").write_text("// linalg fixture\n", encoding="utf-8")
    (worktree / ".config" / "nextest.toml").write_text("[profile.ci]\nretries = 1\n", encoding="utf-8")
    (worktree / "Cargo.toml").write_text(
        '[package]\nname = "fixture"\nversion = "0.1.0"\n[dev-dependencies]\ncriterion = "0.7.0"\n',
        encoding="utf-8",
    )
    (worktree / "Cargo.lock").write_text("version = 4\n", encoding="utf-8")
    (worktree / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.96.0"\n', encoding="utf-8")
    (worktree / "justfile").write_text(
        'bench-save-baseline tag suite="all":\nbench-latest: bench-vs-linalg-la-stack bench-exact\n',
        encoding="utf-8",
    )
    (worktree / "scripts" / "bench_compare.py").write_text('parser.add_argument("--suite")\nparser.add_argument("--scope")\n', encoding="utf-8")


def _write_legacy_benchmark_tooling(worktree: Path) -> None:
    _write_current_benchmark_tooling(worktree)
    (worktree / "justfile").write_text("bench-exact:\n", encoding="utf-8")
    (worktree / "scripts" / "bench_compare.py").write_text('parser.add_argument("--output")\n', encoding="utf-8")


def test_shared_benchmark_harness_replaces_baseline_content_and_has_stable_digest(tmp_path: Path) -> None:
    current = tmp_path / "current"
    baseline = tmp_path / "baseline"
    current.mkdir()
    baseline.mkdir()
    _write_current_benchmark_tooling(current)
    _write_current_benchmark_tooling(baseline)
    (baseline / "benches" / "vs_linalg.rs").write_text("fn obsolete() {}\n", encoding="utf-8")
    (baseline / "Cargo.lock").write_text("version = 3\n", encoding="utf-8")
    (baseline / "justfile").write_text("obsolete-benchmark-recipe:\n", encoding="utf-8")
    (baseline / ".config" / "nextest.toml").write_text("[profile.default]\nretries = 0\n", encoding="utf-8")

    digest = archive_performance._install_shared_benchmark_harness(
        source=current,
        destination=baseline,
    )

    assert digest == archive_performance._benchmark_harness_digest(current)
    assert digest == archive_performance._benchmark_harness_digest(baseline)
    assert (baseline / "benches" / "vs_linalg.rs").read_text(encoding="utf-8") == "fn main() {}\n"
    assert (baseline / "justfile").read_text(encoding="utf-8") == (current / "justfile").read_text(encoding="utf-8")
    assert (baseline / ".config" / "nextest.toml").read_text(encoding="utf-8") == (current / ".config" / "nextest.toml").read_text(encoding="utf-8")


def test_github_release_assets_discard_embedded_shared_harness_metadata(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for tag in ("v0.4.2", "v0.4.3"):
        _write_baseline_archive(
            tmp_path / f"la-stack-{tag}-criterion-baseline.tar.gz",
            include_harness_metadata=True,
        )

    def fake_download(*, baseline_tag: str, download_dir: Path, repo_root: Path) -> Path:
        del download_dir, repo_root
        return tmp_path / f"la-stack-{baseline_tag}-criterion-baseline.tar.gz"

    target_worktree = tmp_path / "worktree"
    target_worktree.mkdir()
    monkeypatch.setattr(archive_performance, "_download_release_baseline", fake_download)

    archive_performance._prepare_github_release_assets(
        current_tag="v0.4.3",
        baseline_tag="v0.4.2",
        repo_root=tmp_path,
        target_worktree=target_worktree,
        tmp_dir=tmp_path,
    )

    criterion_dir = target_worktree / "target" / "criterion"
    assert not (criterion_dir / archive_performance._BENCHMARK_HARNESS_METADATA).exists()
    assert (criterion_dir / "exact_d2" / "det_exact" / "new" / "estimates.json").is_file()


def test_purge_selected_new_samples_preserves_named_baselines_and_other_suites(tmp_path: Path) -> None:
    criterion_dir = tmp_path / "criterion"
    exact_new = criterion_dir / "exact_d2" / "det_exact" / "new"
    exact_baseline = criterion_dir / "exact_d2" / "det_exact" / "v0.4.2"
    linalg_new = criterion_dir / "d2" / "la_stack_lu" / "new"
    for directory in (exact_new, exact_baseline, linalg_new):
        directory.mkdir(parents=True)
        (directory / "estimates.json").write_text("{}\n", encoding="utf-8")

    removed = archive_performance._purge_criterion_new_samples(
        criterion_dir=criterion_dir,
        suite="exact",
    )

    assert removed == [exact_new]
    assert not exact_new.exists()
    assert exact_baseline.is_dir()
    assert linalg_new.is_dir()


def test_apply_current_diff_includes_complete_current_tree_without_mutating_index(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    worktree = tmp_path / "worktree"
    repo_root.mkdir()
    _git(repo_root, "init", "--quiet")
    _git(repo_root, "config", "user.name", "Test User")
    _git(repo_root, "config", "user.email", "test@example.com")
    _git(repo_root, "config", "commit.gpgsign", "false")

    (repo_root / ".gitignore").write_text("ignored.bin\n", encoding="utf-8")
    tracked = repo_root / "tracked.txt"
    tracked.write_text("committed\n", encoding="utf-8")
    _git(repo_root, "add", "--", ".gitignore", "tracked.txt")
    _git(repo_root, "commit", "--quiet", "-m", "initial")
    _git(repo_root, "worktree", "add", "--quiet", "--detach", str(worktree), "HEAD")

    tracked.write_text("staged\n", encoding="utf-8")
    _git(repo_root, "add", "--", "tracked.txt")
    tracked.write_text("working tree\n", encoding="utf-8")
    binary_payload = bytes(range(256)) * 2
    binary = repo_root / "--untracked.bin"
    binary.write_bytes(binary_payload)
    link = repo_root / "untracked-link"
    try:
        link.symlink_to(binary.name)
    except OSError as exc:
        pytest.skip(f"symlinks unavailable: {exc}")
    (repo_root / "ignored.bin").write_bytes(b"ignored\x00payload")

    status_before = _git(repo_root, "status", "--porcelain=v1", "--untracked-files=all")
    archive_performance._apply_current_diff_to_worktree(repo_root=repo_root, worktree=worktree)

    assert (worktree / "tracked.txt").read_text(encoding="utf-8") == "working tree\n"
    assert (worktree / binary.name).read_bytes() == binary_payload
    applied_link = worktree / link.name
    assert applied_link.is_symlink()
    assert applied_link.readlink() == Path(binary.name)
    assert not (worktree / "ignored.bin").exists()
    assert _git(repo_root, "show", ":tracked.txt") == "staged\n"
    assert _git(repo_root, "status", "--porcelain=v1", "--untracked-files=all") == status_before


def test_apply_current_diff_preserves_crlf_patch_bytes(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    worktree = tmp_path / "worktree"
    repo_root.mkdir()
    _git(repo_root, "init", "--quiet")
    _git(repo_root, "config", "user.name", "Test User")
    _git(repo_root, "config", "user.email", "test@example.com")
    _git(repo_root, "config", "commit.gpgsign", "false")
    _git(repo_root, "config", "core.autocrlf", "false")

    (repo_root / ".gitattributes").write_text("tracked.txt -text\n", encoding="utf-8")
    tracked = repo_root / "tracked.txt"
    tracked.write_bytes(b"committed\r\n")
    _git(repo_root, "add", "--", ".gitattributes", tracked.name)
    _git(repo_root, "commit", "--quiet", "-m", "initial")
    _git(repo_root, "worktree", "add", "--quiet", "--detach", str(worktree), "HEAD")

    tracked.write_bytes(b"staged\r\n")
    _git(repo_root, "add", "--", tracked.name)
    tracked.write_bytes(b"working tree\r\n")
    index_before = _git(repo_root, "rev-parse", f":{tracked.name}")

    archive_performance._apply_current_diff_to_worktree(repo_root=repo_root, worktree=worktree)

    assert (worktree / tracked.name).read_bytes() == b"working tree\r\n"
    assert _git(repo_root, "rev-parse", f":{tracked.name}") == index_before


def test_apply_current_diff_fails_loudly_and_cleans_temporary_index(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = tmp_path / "repo"
    worktree = tmp_path / "worktree"
    repo_root.mkdir()
    worktree.mkdir()
    temporary_index: Path | None = None

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        nonlocal temporary_index
        assert cwd == repo_root
        env = kwargs["env"]
        temporary_index = Path(env["GIT_INDEX_FILE"])
        assert temporary_index.parent.is_dir()
        if args == ["add", "--all", "--", "."]:
            raise subprocess.CalledProcessError(
                128,
                ["git", *args],
                output="snapshot stdout",
                stderr="cannot snapshot current tree",
            )
        return _result()

    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)

    with pytest.raises(RuntimeError) as exc_info:
        archive_performance._apply_current_diff_to_worktree(repo_root=repo_root, worktree=worktree)

    error = str(exc_info.value)
    assert "command failed (128): git add --all -- ." in error
    assert "snapshot stdout" in error
    assert "cannot snapshot current tree" in error
    assert temporary_index is not None
    assert not temporary_index.parent.exists()


@pytest.mark.parametrize(
    ("failure", "message"),
    [
        (archive_performance.ExecutableNotFoundError("missing tool"), "command could not start: tool --flag: missing tool"),
        (subprocess.TimeoutExpired(["tool", "--flag"], 17, stderr="stalled"), "command timed out after 17 seconds: tool --flag"),
        (OSError("working directory unavailable"), "command could not start: tool --flag: working directory unavailable"),
    ],
)
def test_run_tool_normalizes_launch_timeout_and_os_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: Exception,
    message: str,
) -> None:
    def fail_run(*_args: object, **_kwargs: object) -> SimpleNamespace:
        raise failure

    monkeypatch.setattr(archive_performance, "run_safe_command", fail_run)

    with pytest.raises(RuntimeError) as exc_info:
        archive_performance._run_tool("tool", ["--flag"], cwd=tmp_path)

    assert str(exc_info.value).startswith(message)
    assert exc_info.value.__cause__ is failure


def test_temporary_worktree_cleanup_failure_fails_successful_operation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run_git(args: list[str], *, cwd: Path, timeout: int = 600) -> None:
        del cwd, timeout
        if args[:3] == ["worktree", "remove", "--force"]:
            msg = "cleanup failed"
            raise RuntimeError(msg)

    monkeypatch.setattr(archive_performance, "_run_git", fake_run_git)

    def complete_operation() -> None:
        with archive_performance._temporary_detached_worktree(
            repo_root=tmp_path,
            worktree=tmp_path / "worktree",
            revision="HEAD",
            label="test worktree",
        ):
            pass

    with pytest.raises(RuntimeError, match="failed to remove test worktree"):
        complete_operation()


def test_temporary_worktree_cleanup_does_not_mask_primary_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run_git(args: list[str], *, cwd: Path, timeout: int = 600) -> None:
        del cwd, timeout
        if args[:3] == ["worktree", "remove", "--force"]:
            msg = "cleanup failed"
            raise RuntimeError(msg)

    monkeypatch.setattr(archive_performance, "_run_git", fake_run_git)

    def fail_operation() -> None:
        with archive_performance._temporary_detached_worktree(
            repo_root=tmp_path,
            worktree=tmp_path / "worktree",
            revision="HEAD",
            label="test worktree",
        ):
            msg = "primary failed"
            raise ValueError(msg)

    with pytest.raises(RuntimeError, match="operation failed") as exc_info:
        fail_operation()

    assert "additionally failed to remove test worktree: cleanup failed" in str(exc_info.value)
    assert isinstance(exc_info.value.__cause__, ValueError)


def test_normalize_tag_adds_leading_v() -> None:
    assert normalize_tag("0.4.2") == "v0.4.2"
    assert normalize_tag("v0.4.2") == "v0.4.2"
    assert normalize_tag("v1.2.3-rc.1+build.7") == "v1.2.3-rc.1+build.7"


def test_normalize_tag_rejects_non_semver_baseline_names() -> None:
    with pytest.raises(ValueError, match="semver tag"):
        normalize_tag("last")


def test_parse_report_id_reads_current_and_baseline_tags() -> None:
    report_id = parse_report_id(_report("0.4.2", "v0.4.1"))

    assert report_id.current_tag == "v0.4.2"
    assert report_id.baseline_tag == "v0.4.1"
    assert report_id.archive_name == "v0.4.2-vs-v0.4.1.md"


def test_published_release_pair_discovers_latest_stable_semver_pair(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        assert command == "gh"
        assert args == [
            "release",
            "list",
            "--json",
            "tagName,isDraft,isPrerelease,publishedAt",
            "--limit",
            "100",
        ]
        assert cwd == tmp_path
        return _result(
            "["
            '{"tagName":"v0.4.2","isDraft":false,"isPrerelease":false,"publishedAt":"2026-01-01T00:00:00Z"},'
            '{"tagName":"v0.4.10","isDraft":false,"isPrerelease":false,"publishedAt":"2026-04-01T00:00:00Z"},'
            '{"tagName":"v0.4.11-rc.1","isDraft":false,"isPrerelease":true,"publishedAt":"2026-06-01T00:00:00Z"},'
            '{"tagName":"v0.4.11","isDraft":true,"isPrerelease":false,"publishedAt":"2026-06-02T00:00:00Z"},'
            '{"tagName":"not-semver","isDraft":false,"isPrerelease":false,"publishedAt":"2026-06-03T00:00:00Z"},'
            '{"tagName":"v0.4.3","isDraft":false,"isPrerelease":false,"publishedAt":"2026-03-01T00:00:00Z"}'
            "]"
        )

    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    report_id = archive_performance._published_release_pair(tmp_path)

    assert report_id.current_tag == "v0.4.10"
    assert report_id.baseline_tag == "v0.4.3"


def test_published_release_pair_uses_latest_published_release_not_highest_semver(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        assert command == "gh"
        assert cwd == tmp_path
        return _result(
            "["
            '{"tagName":"v0.5.0","isDraft":false,"isPrerelease":false,"publishedAt":"2026-01-01T00:00:00Z"},'
            '{"tagName":"v0.4.9","isDraft":false,"isPrerelease":false,"publishedAt":"2026-02-01T00:00:00Z"},'
            '{"tagName":"v0.4.8","isDraft":false,"isPrerelease":false,"publishedAt":"2025-12-01T00:00:00Z"}'
            "]"
        )

    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    report_id = archive_performance._published_release_pair(tmp_path)

    assert report_id.current_tag == "v0.4.9"
    assert report_id.baseline_tag == "v0.4.8"


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("isDraft", "false"),
        ("isDraft", 0),
        ("isPrerelease", None),
    ],
)
def test_stable_published_releases_requires_boolean_flags(field: str, value: object) -> None:
    release: dict[str, object] = {
        "tagName": "v0.4.3",
        "isDraft": False,
        "isPrerelease": False,
        "publishedAt": "2026-02-01T00:00:00Z",
    }
    release[field] = value

    with pytest.raises(TypeError, match="boolean isDraft and isPrerelease"):
        archive_performance._stable_published_releases([release])


@pytest.mark.parametrize("published_at", ["2026-02-01T00:00:00", "not-a-timestamp", ""])
def test_stable_published_releases_requires_aware_timestamp(published_at: str) -> None:
    release = {
        "tagName": "v0.4.3",
        "isDraft": False,
        "isPrerelease": False,
        "publishedAt": published_at,
    }

    with pytest.raises((TypeError, ValueError), match="publishedAt"):
        archive_performance._stable_published_releases([release])


def test_stable_published_releases_normalizes_timestamp_to_utc() -> None:
    releases = archive_performance._stable_published_releases(
        [
            {
                "tagName": "v0.4.3",
                "isDraft": False,
                "isPrerelease": False,
                "publishedAt": "2026-02-01T01:00:00+01:00",
            }
        ]
    )

    assert releases[0].published_at.isoformat() == "2026-02-01T00:00:00+00:00"


def test_resolve_archive_request_infer_release_uses_package_version_and_previous_release(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "Cargo.toml").write_text('[package]\nversion = "0.4.3"\n', encoding="utf-8")

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        assert command == "gh"
        assert args[:2] == ["release", "list"]
        assert cwd == tmp_path
        return _result(
            "["
            '{"tagName":"v0.4.1","isDraft":false,"isPrerelease":false,"publishedAt":"2026-01-01T00:00:00Z"},'
            '{"tagName":"v0.4.2","isDraft":false,"isPrerelease":false,"publishedAt":"2026-02-01T00:00:00Z"}'
            "]"
        )

    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    request = archive_performance.resolve_archive_request(
        archive_performance.ArchiveRequestOptions(
            current_tag=None,
            baseline_tag=None,
            published_latest=False,
            infer_release=True,
            current_vs_latest=False,
            worktree_ref="HEAD",
            repo_root=tmp_path,
        )
    )

    assert request.current_tag == "v0.4.3"
    assert request.baseline_tag == "v0.4.2"
    assert request.worktree_ref == "HEAD"
    assert request.tags_to_fetch == ("v0.4.2",)


def test_resolve_archive_request_current_vs_latest_uses_package_version_and_latest_release(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "Cargo.toml").write_text('[package]\nversion = "0.4.3"\n', encoding="utf-8")

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        assert command == "gh"
        assert args[:2] == ["release", "list"]
        assert cwd == tmp_path
        return _result(
            "["
            '{"tagName":"v0.4.1","isDraft":false,"isPrerelease":false,"publishedAt":"2026-01-01T00:00:00Z"},'
            '{"tagName":"v0.4.2","isDraft":false,"isPrerelease":false,"publishedAt":"2026-02-01T00:00:00Z"}'
            "]"
        )

    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    request = archive_performance.resolve_archive_request(
        archive_performance.ArchiveRequestOptions(
            current_tag=None,
            baseline_tag=None,
            published_latest=False,
            infer_release=False,
            current_vs_latest=True,
            worktree_ref="HEAD",
            repo_root=tmp_path,
        )
    )

    assert request.current_tag == "v0.4.3"
    assert request.baseline_tag == "v0.4.2"
    assert request.worktree_ref == "HEAD"
    assert request.tags_to_fetch == ("v0.4.2",)


def test_benchmark_env_uses_current_repo_toolchain(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RUSTUP_TOOLCHAIN", raising=False)
    (tmp_path / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.96.0"\n', encoding="utf-8")

    env = archive_performance._benchmark_env(tmp_path)

    assert env is not None
    assert env["RUSTUP_TOOLCHAIN"] == "1.96.0"


def test_benchmark_env_respects_existing_toolchain_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("RUSTUP_TOOLCHAIN", "nightly")
    (tmp_path / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.96.0"\n', encoding="utf-8")

    assert archive_performance._benchmark_env(tmp_path) is None


def test_parser_rejects_unsupported_scope() -> None:
    with pytest.raises(SystemExit):
        archive_performance.build_parser().parse_args(["v0.4.3", "v0.4.2", "--scope", "quick"])


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("suite", "other", "unsupported benchmark suite"),
        ("scope", "quick", "unsupported comparison scope"),
    ],
)
def test_generation_config_rejects_unsupported_benchmark_selection(
    tmp_path: Path,
    field: str,
    value: str,
    message: str,
) -> None:
    kwargs: dict[str, Any] = {
        "repo_root": tmp_path,
        "current_tag": "v0.4.3",
        "baseline_tag": "v0.4.2",
        "worktree_ref": "HEAD",
    }
    kwargs[field] = value

    with pytest.raises(ValueError, match=message):
        GenerationConfig(**kwargs)


def test_comparison_benchmark_env_preserves_flags_and_selects_v043_adapter(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("CARGO_ENCODED_RUSTFLAGS", raising=False)
    monkeypatch.delenv("RUSTUP_TOOLCHAIN", raising=False)
    monkeypatch.setenv("RUSTFLAGS", "-C target-cpu=native")
    (tmp_path / "rust-toolchain.toml").write_text(
        '[toolchain]\nchannel = "1.96.0"\n',
        encoding="utf-8",
    )

    current = archive_performance._comparison_benchmark_env(tmp_path)
    baseline = archive_performance._comparison_benchmark_env(
        tmp_path,
        baseline_tag="v0.4.3",
    )

    assert current["RUSTFLAGS"] == "-C target-cpu=native --cap-lints=warn"
    assert baseline["RUSTFLAGS"] == ("-C target-cpu=native --cap-lints=warn --cfg=la_stack_v0_4_3_api")
    assert current["RUSTUP_TOOLCHAIN"] == "1.96.0"


def test_comparison_benchmark_env_extends_encoded_rustflags(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("RUSTUP_TOOLCHAIN", "nightly")
    monkeypatch.setenv("CARGO_ENCODED_RUSTFLAGS", "-C\x1ftarget-cpu=native")

    env = archive_performance._comparison_benchmark_env(
        tmp_path,
        baseline_tag="0.4.3",
    )

    assert env["CARGO_ENCODED_RUSTFLAGS"] == ("-C\x1ftarget-cpu=native\x1f--cap-lints=warn\x1f--cfg=la_stack_v0_4_3_api")
    assert env["RUSTUP_TOOLCHAIN"] == "nightly"


@pytest.mark.parametrize("suite", ["exact", "vs_linalg"])
def test_fallback_baseline_cargo_commands_enforce_lockfile(
    tmp_path: Path,
    suite: str,
) -> None:
    worktree = tmp_path / "legacy"
    worktree.mkdir()
    _write_legacy_benchmark_tooling(worktree)

    command, args = archive_performance._baseline_tool_args(
        baseline_tag="v0.4.2",
        suite=suite,
        baseline_worktree=worktree,
    )

    assert command == "cargo"
    assert args[:2] == ["bench", "--locked"]


@pytest.mark.parametrize(
    ("suite", "expected"),
    [
        ("all", ("cargo", "bench", "--locked", "--features", "bench,exact")),
        ("exact", ("cargo", "bench", "--locked", "--features", "bench,exact", "--bench", "exact")),
        ("vs_linalg", ("cargo", "bench", "--locked", "--features", "bench", "--bench", "vs_linalg")),
    ],
)
def test_fallback_current_cargo_command_matches_suite(
    suite: str,
    expected: tuple[str, ...],
) -> None:
    assert archive_performance._fallback_current_command(suite=suite) == expected


def test_promote_report_archives_previous_and_updates_sorted_index(tmp_path: Path) -> None:
    source = tmp_path / "target" / "bench-reports" / "performance.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"

    source.parent.mkdir(parents=True)
    current.parent.mkdir(parents=True)
    archive_dir.mkdir(parents=True)
    source.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
    current.write_text(_report("0.4.1", "v0.4.0"), encoding="utf-8")
    (archive_dir / "v0.3.1-vs-v0.3.0.md").write_text(_report("0.3.1", "v0.3.0"), encoding="utf-8")

    promoted = promote_report(
        source=source,
        current=current,
        archive_dir=archive_dir,
        expected_current_tag="v0.4.2",
        expected_baseline_tag="v0.4.1",
    )

    assert promoted.archive_name == "v0.4.2-vs-v0.4.1.md"
    assert current.read_text(encoding="utf-8") == _normalized_report("0.4.2", "v0.4.1")
    assert (archive_dir / "v0.4.1-vs-v0.4.0.md").read_text(encoding="utf-8") == _normalized_report("0.4.1", "v0.4.0")
    assert (archive_dir / "README.md").read_text(encoding="utf-8") == (
        "# Archived Performance Reports\n\n"
        "Older release-to-release benchmark comparisons are archived here.\n"
        "`docs/PERFORMANCE.md` contains the latest curated comparison.\n\n"
        "- [v0.3.1-vs-v0.3.0](v0.3.1-vs-v0.3.0.md)\n"
        "- [v0.4.1-vs-v0.4.0](v0.4.1-vs-v0.4.0.md)\n"
    )


def test_promote_report_is_idempotent_for_same_release_pair(tmp_path: Path) -> None:
    source = tmp_path / "performance-new.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"

    source.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")

    promote_report(
        source=source,
        current=current,
        archive_dir=archive_dir,
        expected_current_tag="v0.4.2",
        expected_baseline_tag="v0.4.1",
    )

    assert not (archive_dir / "v0.4.2-vs-v0.4.1.md").exists()
    assert "- No archived performance reports yet." in (archive_dir / "README.md").read_text(encoding="utf-8")


def test_promote_report_does_not_overwrite_existing_archive(tmp_path: Path) -> None:
    source = tmp_path / "performance-new.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    archived = archive_dir / "v0.4.1-vs-v0.4.0.md"

    source.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(_report("0.4.1", "v0.4.0"), encoding="utf-8")
    archive_dir.mkdir(parents=True)
    archived.write_text("already archived\n", encoding="utf-8")

    promote_report(
        source=source,
        current=current,
        archive_dir=archive_dir,
        expected_current_tag="v0.4.2",
        expected_baseline_tag="v0.4.1",
    )

    assert archived.read_text(encoding="utf-8") == "already archived\n"


def test_promote_report_rejects_unexpected_release_pair(tmp_path: Path) -> None:
    source = tmp_path / "performance-new.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    source.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")

    with pytest.raises(ValueError, match="does not match requested release pair"):
        promote_report(
            source=source,
            current=current,
            archive_dir=archive_dir,
            expected_current_tag="v0.4.3",
            expected_baseline_tag="v0.4.2",
        )


def test_promote_report_rewrites_legacy_update_instructions(tmp_path: Path) -> None:
    source = tmp_path / "performance-new.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    source.write_text(_legacy_report("0.4.3", "v0.4.2"), encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(_legacy_report("0.4.2", "v0.4.1"), encoding="utf-8")

    promote_report(
        source=source,
        current=current,
        archive_dir=archive_dir,
        expected_current_tag="v0.4.3",
        expected_baseline_tag="v0.4.2",
    )

    current_text = current.read_text(encoding="utf-8")
    archived_text = (archive_dir / "v0.4.2-vs-v0.4.1.md").read_text(encoding="utf-8")
    assert "just performance-local" in current_text
    assert "just performance-release" in current_text
    assert "just performance-github-assets" in current_text
    assert "just performance-release <current-tag> <previous-tag>" in current_text
    assert "git checkout" not in current_text
    assert "just performance-local" in archived_text
    assert "just performance-github-assets" in archived_text
    assert "git checkout" not in archived_text


def test_main_promotes_generated_report_to_docs_performance(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "target" / "bench-reports" / "performance.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    generated = _report("0.4.3", "v0.4.2")

    source.parent.mkdir(parents=True)
    source.write_text(generated, encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
    gates: list[Path] = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        archive_performance,
        "_run_benchmark_input_gate",
        lambda checkout, **_kwargs: gates.append(checkout),
    )

    rc = main(
        [
            "v0.4.3",
            "v0.4.2",
            "--source",
            str(source),
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
        ]
    )

    assert rc == 0
    assert current.read_text(encoding="utf-8") == archive_performance._normalize_how_to_update(generated)
    assert (archive_dir / "v0.4.2-vs-v0.4.1.md").exists()
    assert gates == [tmp_path]
    assert "Current performance report: v0.4.3 vs v0.4.2" in capsys.readouterr().out


def test_main_reports_release_pair_mismatch_to_stderr(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "target" / "bench-reports" / "performance.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    source.parent.mkdir(parents=True)
    source.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
    monkeypatch.setattr(archive_performance, "_run_benchmark_input_gate", lambda *_args, **_kwargs: None)

    rc = main(
        [
            "v0.4.4",
            "v0.4.3",
            "--source",
            str(source),
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
        ]
    )

    captured = capsys.readouterr()
    assert rc == 1
    assert "does not match requested release pair" in captured.err
    assert not current.exists()


def test_main_reraises_unexpected_errors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_unexpected(*, args: object, paths: object, request: object, repo_root: Path) -> archive_performance.ArchiveResult:
        msg = "unexpected test failure"
        raise AssertionError(msg)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "_run_archive_request", fail_unexpected)

    with pytest.raises(AssertionError, match="unexpected test failure"):
        main(["v0.4.3", "v0.4.2"])


def test_main_generates_report_in_temp_worktree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.2", "exact"]:
            assert cwd is not None
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "v0.4.3",
            "v0.4.2",
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--generate-in-temp-worktree",
            "--worktree-ref",
            "v0.4.3",
            "--no-apply-current-diff",
            "--suite",
            "exact",
            "--scope",
            "release-signal",
        ]
    )

    captured = capsys.readouterr()
    assert rc == 0
    assert current.read_text(encoding="utf-8") == _normalized_report("0.4.3", "v0.4.2")
    assert "Generated benchmark report in a temporary worktree" in captured.out
    assert "target/bench-reports/performance.md" not in captured.out
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "v0.4.3" for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-exact",) for kind, args, _ in calls)
    assert any(kind == "uv" and "--suite" in args and args[args.index("--suite") + 1] == "exact" for kind, args, _ in calls)
    assert any(kind == "uv" and args[:2] == ("run", "--locked") for kind, args, _ in calls)
    assert not any(kind == "git" and args[:1] == ("read-tree",) for kind, args, _ in calls)
    assert not any(kind == "git-stdin" for kind, _, _ in calls)


def test_temp_worktree_is_removed_when_benchmark_command_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:2] == ["fetch", "origin"]:
            return _result()
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            if worktree.name == "baseline-worktree":
                _write_legacy_benchmark_tooling(worktree)
            else:
                _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.2"]:
            assert cwd is not None
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "just" and args == ["bench-latest"]:
            raise subprocess.CalledProcessError(42, [command, *args], output="bench stdout", stderr="bench stderr")
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "v0.4.3",
            "v0.4.2",
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--generate-in-temp-worktree",
            "--worktree-ref",
            "HEAD",
            "--no-apply-current-diff",
        ]
    )

    captured = capsys.readouterr()
    assert rc == 1
    assert "command failed (42): just bench-latest" in captured.err
    assert "bench stderr" in captured.err
    assert not current.exists()
    assert any(kind == "git" and args[:3] == ("worktree", "remove", "--force") for kind, args, _ in calls)


def test_generate_report_rejects_unsafe_baseline_archive(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            tag = args[2]
            _write_unsafe_baseline_archive(download_dir / f"la-stack-{tag}-criterion-baseline.tar.gz")
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "v0.4.3",
            "v0.4.2",
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--generate-in-temp-worktree",
            "--worktree-ref",
            "HEAD",
            "--no-apply-current-diff",
            "--github-assets",
        ]
    )

    captured = capsys.readouterr()
    assert rc == 1
    assert "refusing to extract unsafe archive member '../escape.txt'" in captured.err
    assert not (tmp_path / "escape.txt").exists()
    assert not current.exists()
    assert not any(kind in {"just", "uv"} for kind, _, _ in calls)
    assert any(kind == "git" and args[:3] == ("worktree", "remove", "--force") for kind, args, _ in calls)


def test_generate_report_generates_release_baseline_locally(  # noqa: PLR0915
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.delenv("RUSTUP_TOOLCHAIN", raising=False)
    (tmp_path / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.96.0"\n', encoding="utf-8")
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            if worktree.name == "baseline-worktree":
                _write_legacy_benchmark_tooling(worktree)
            else:
                _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.2"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
            assert cwd is not None
            assert "bench-latest" in (cwd / "justfile").read_text(encoding="utf-8")
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
            for sample in ("new", "v0.4.2"):
                estimates = criterion_dir / "exact_d2" / "det_exact" / sample / "estimates.json"
                estimates.parent.mkdir(parents=True, exist_ok=True)
                estimates.write_text("{}\n", encoding="utf-8")
        if command == "just" and args == ["bench-latest"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
            assert cwd is not None
            criterion_dir = cwd / "target" / "criterion" / "exact_d2" / "det_exact"
            assert not (criterion_dir / "new").exists()
            assert (criterion_dir / "v0.4.2" / "estimates.json").is_file()
        if command == "uv":
            assert cwd is not None
            metadata = json.loads((cwd / "target" / "criterion" / archive_performance._BENCHMARK_HARNESS_METADATA).read_text(encoding="utf-8"))
            assert metadata["baseline"] == "v0.4.2"
            assert metadata["mode"] == "shared-current-harness"
            assert metadata["schema"] == 2
            assert metadata["measurement"]["harness_sha256"] == archive_performance._benchmark_harness_digest(cwd)
            assert metadata["measurement"]["current_source_state_sha256"] == archive_performance._source_state_digest(cwd)
            assert metadata["criterion"]["criterion_version"] == "manifest requirement 0.7.0"
            assert metadata["validation"]["baseline_revision"] == "passed"
            assert metadata["validation"]["current_revision"] == "passed"
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "v0.4.3",
            "v0.4.2",
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--generate-in-temp-worktree",
            "--worktree-ref",
            "HEAD",
            "--no-apply-current-diff",
        ]
    )

    captured = capsys.readouterr()
    assert rc == 0
    assert captured.err == ""
    assert current.read_text(encoding="utf-8") == _normalized_report("0.4.3", "v0.4.2")
    assert not any(kind == "gh" for kind, _, _ in calls)
    assert any(kind == "just" and args == ("bench-save-baseline", "v0.4.2") for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)
    assert any(kind == "uv" and "--suite" in args for kind, args, _ in calls)
    baseline_timing_index = next(index for index, (kind, args, _) in enumerate(calls) if kind == "just" and args == ("bench-save-baseline", "v0.4.2"))
    current_timing_index = next(index for index, (kind, args, _) in enumerate(calls) if kind == "just" and args == ("bench-latest",))
    baseline_gate_index = next(
        index
        for index, (kind, args, cwd) in enumerate(calls)
        if kind == "just" and args == ("test-bench-inputs",) and cwd is not None and cwd.name == "baseline-worktree"
    )
    current_gate_index = next(
        index
        for index, (kind, args, cwd) in enumerate(calls)
        if kind == "just" and args == ("test-bench-inputs",) and cwd is not None and cwd.name == "worktree"
    )
    assert baseline_gate_index < baseline_timing_index
    assert current_gate_index < current_timing_index
    assert sum(1 for kind, args, _ in calls if kind == "git" and args[:3] == ("worktree", "remove", "--force")) == 2


def test_generate_report_vs_linalg_suite_uses_copied_current_baseline_recipe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RUSTUP_TOOLCHAIN", raising=False)
    (tmp_path / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.96.0"\n', encoding="utf-8")
    output = tmp_path / "target" / "bench-reports" / "performance.md"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            if worktree.name == "baseline-worktree":
                _write_legacy_benchmark_tooling(worktree)
            else:
                _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.2", "vs_linalg"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
            assert cwd is not None
            assert "bench-latest" in (cwd / "justfile").read_text(encoding="utf-8")
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "just" and args == ["bench-vs-linalg-la-stack"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
        if command == "uv":
            report = Path(args[args.index("--output") + 1])
            report.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    report_id = generate_worktree_report(
        output=output,
        config=GenerationConfig(
            repo_root=tmp_path,
            current_tag="v0.4.3",
            baseline_tag="v0.4.2",
            worktree_ref="HEAD",
            suite="vs_linalg",
            apply_current_diff=False,
        ),
    )

    assert report_id.archive_name == "v0.4.3-vs-v0.4.2.md"
    assert output.read_text(encoding="utf-8") == _normalized_report("0.4.3", "v0.4.2")
    assert any(kind == "just" and args == ("bench-save-baseline", "v0.4.2", "vs_linalg") for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-vs-linalg-la-stack",) for kind, args, _ in calls)
    assert not any(kind == "cargo" for kind, _, _ in calls)
    assert not any(kind == "just" and args == ("bench-exact",) for kind, args, _ in calls)
    assert not any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)


def test_generate_report_vs_linalg_suite_uses_suite_aware_baseline_recipe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output = tmp_path / "target" / "bench-reports" / "performance.md"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.2", "vs_linalg"]:
            assert cwd is not None
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "uv":
            report = Path(args[args.index("--output") + 1])
            report.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    report_id = generate_worktree_report(
        output=output,
        config=GenerationConfig(
            repo_root=tmp_path,
            current_tag="v0.4.3",
            baseline_tag="v0.4.2",
            worktree_ref="HEAD",
            suite="vs_linalg",
            apply_current_diff=False,
        ),
    )

    assert report_id.archive_name == "v0.4.3-vs-v0.4.2.md"
    assert any(kind == "just" and args == ("bench-save-baseline", "v0.4.2", "vs_linalg") for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-vs-linalg-la-stack",) for kind, args, _ in calls)
    assert not any(kind == "cargo" for kind, _, _ in calls)
    assert not any(kind == "just" and args == ("bench-exact",) for kind, args, _ in calls)


def test_main_generates_latest_published_report_from_github_releases(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "gh" and args[:2] == ["release", "list"]:
            return _result(
                "["
                '{"tagName":"v0.4.2","isDraft":false,"isPrerelease":false,"publishedAt":"2026-01-01T00:00:00Z"},'
                '{"tagName":"v0.4.3","isDraft":false,"isPrerelease":false,"publishedAt":"2026-02-01T00:00:00Z"}'
                "]"
            )
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            tag = args[2]
            _write_baseline_archive(download_dir / f"la-stack-{tag}-criterion-baseline.tar.gz")
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--published-latest",
            "--github-assets",
            "--generate-in-temp-worktree",
            "--no-apply-current-diff",
        ]
    )

    assert rc == 0
    assert current.read_text(encoding="utf-8") == _normalized_report("0.4.3", "v0.4.2")
    assert any(
        kind == "git"
        and args
        == (
            "fetch",
            "origin",
            "refs/tags/v0.4.3:refs/tags/v0.4.3",
            "refs/tags/v0.4.2:refs/tags/v0.4.2",
        )
        for kind, args, _ in calls
    )
    fetch_index = next(index for index, (kind, args, _) in enumerate(calls) if kind == "git" and args[:2] == ("fetch", "origin"))
    worktree_index = next(index for index, (kind, args, _) in enumerate(calls) if kind == "git" and args[:3] == ("worktree", "add", "--detach"))
    assert fetch_index < worktree_index
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "v0.4.3" for kind, args, _ in calls)


def test_main_normalizes_explicit_bare_tags_before_fetching_and_checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output = tmp_path / "target" / "bench-reports" / "github-assets-performance.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            tag = args[2]
            _write_baseline_archive(download_dir / f"la-stack-{tag}-criterion-baseline.tar.gz")
        if command == "uv":
            report = Path(args[args.index("--output") + 1])
            report.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "0.4.3",
            "0.4.2",
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--github-assets",
            "--generate-in-temp-worktree",
            "--worktree-ref",
            "0.4.3",
            "--output-only",
            "--output",
            str(output),
        ]
    )

    assert rc == 0
    assert output.read_text(encoding="utf-8") == _normalized_report("0.4.3", "v0.4.2")
    assert not current.exists()
    assert any(
        kind == "git"
        and args
        == (
            "fetch",
            "origin",
            "refs/tags/v0.4.2:refs/tags/v0.4.2",
            "refs/tags/v0.4.3:refs/tags/v0.4.3",
        )
        for kind, args, _ in calls
    )
    fetch_index = next(index for index, (kind, args, _) in enumerate(calls) if kind == "git" and args[:2] == ("fetch", "origin"))
    worktree_index = next(index for index, (kind, args, _) in enumerate(calls) if kind == "git" and args[:3] == ("worktree", "add", "--detach"))
    assert fetch_index < worktree_index
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "v0.4.3" for kind, args, _ in calls)
    assert any(kind == "gh" and args[:3] == ("release", "download", "v0.4.2") for kind, args, _ in calls)
    assert any(kind == "gh" and args[:3] == ("release", "download", "v0.4.3") for kind, args, _ in calls)


def test_main_published_latest_fetch_failure_stops_before_worktree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:2] == ["fetch", "origin"]:
            raise subprocess.CalledProcessError(128, ["git", *args], output="fetch stdout", stderr="missing tag")
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "gh" and args[:2] == ["release", "list"]:
            return _result(
                "["
                '{"tagName":"v0.4.2","isDraft":false,"isPrerelease":false,"publishedAt":"2026-01-01T00:00:00Z"},'
                '{"tagName":"v0.4.3","isDraft":false,"isPrerelease":false,"publishedAt":"2026-02-01T00:00:00Z"}'
                "]"
            )
        return _result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    rc = main(
        [
            "--current",
            str(current),
            "--archive-dir",
            str(archive_dir),
            "--published-latest",
            "--generate-in-temp-worktree",
            "--no-apply-current-diff",
        ]
    )

    captured = capsys.readouterr()
    assert rc == 1
    assert "command failed (128): git fetch origin" in captured.err
    assert "missing tag" in captured.err
    assert not current.exists()
    assert not any(kind == "git" and args[:3] == ("worktree", "add", "--detach") for kind, args, _ in calls)
    assert not any(kind in {"just", "uv"} for kind, _, _ in calls)
    assert not any(kind == "git-stdin" for kind, _, _ in calls)


def test_failed_atomic_replace_preserves_existing_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "performance-new.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    original = _report("0.4.2", "v0.4.1")

    source.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(original, encoding="utf-8")

    def fail_replace(src: Path, dst: Path) -> None:
        msg = f"simulated replace failure for {dst}"
        raise OSError(msg)

    monkeypatch.setattr(archive_performance, "_replace_file", fail_replace)

    with pytest.raises(OSError, match="simulated replace failure"):
        promote_report(
            source=source,
            current=current,
            archive_dir=archive_dir,
            expected_current_tag="v0.4.3",
            expected_baseline_tag="v0.4.2",
        )

    assert current.read_text(encoding="utf-8") == original
    assert not list(current.parent.glob(".PERFORMANCE.md.*.tmp"))


def test_generate_and_promote_uses_temp_worktree_and_current_diff(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RUSTUP_TOOLCHAIN", raising=False)
    (tmp_path / "rust-toolchain.toml").write_text('[toolchain]\nchannel = "1.96.0"\n', encoding="utf-8")
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    current.parent.mkdir(parents=True)
    current.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        if args[:3] == ["diff", "--cached", "--binary"]:
            output_arg = next(arg for arg in args if arg.startswith("--output="))
            Path(output_arg.removeprefix("--output=")).write_bytes(b"diff --git a/README.md b/README.md\n")
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str | bytes, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        assert b"diff --git" in input_data if isinstance(input_data, bytes) else "diff --git" in input_data
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.2"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
            assert cwd is not None
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "just" and args == ["bench-latest"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return _result()

    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    report_id = generate_and_promote_worktree_report(
        current=current,
        archive_dir=archive_dir,
        config=GenerationConfig(
            repo_root=tmp_path,
            current_tag="v0.4.3",
            baseline_tag="v0.4.2",
            worktree_ref="HEAD",
            apply_current_diff=True,
        ),
    )

    assert report_id.archive_name == "v0.4.3-vs-v0.4.2.md"
    assert current.read_text(encoding="utf-8") == _normalized_report("0.4.3", "v0.4.2")
    assert (archive_dir / "v0.4.2-vs-v0.4.1.md").exists()
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "HEAD" for kind, args, _ in calls)
    assert any(kind == "git-stdin" and args == ("apply", "--binary") for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)
    assert any(kind == "git" and args[:3] == ("worktree", "remove", "--force") for kind, args, _ in calls)


def test_generate_and_promote_legacy_published_tag_uses_legacy_commands(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[RunnerCall] = []

    def fake_run_git(args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_legacy_benchmark_tooling(worktree)
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        return _result()

    def fake_run_safe(command: str, args: Sequence[str], cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append((command, tuple(args), cwd))
        if command == "just" and args == ["bench-save-baseline", "v0.4.1"]:
            assert cwd is not None
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
        return _result()

    monkeypatch.setattr(archive_performance, "run_git_command", fake_run_git)
    monkeypatch.setattr(archive_performance, "run_git_command_with_input", fake_run_git_with_input)
    monkeypatch.setattr(archive_performance, "run_safe_command", fake_run_safe)

    report_id = generate_and_promote_worktree_report(
        current=current,
        archive_dir=archive_dir,
        config=GenerationConfig(
            repo_root=tmp_path,
            current_tag="v0.4.2",
            baseline_tag="v0.4.1",
            worktree_ref="v0.4.2",
            apply_current_diff=False,
        ),
    )

    assert report_id.archive_name == "v0.4.2-vs-v0.4.1.md"
    assert current.read_text(encoding="utf-8") == _normalized_report("0.4.2", "v0.4.1")
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "v0.4.2" for kind, args, _ in calls)
    assert any(kind == "cargo" and args == ("bench", "--locked", "--features", "bench,exact") for kind, args, _ in calls)
    assert not any(kind == "just" and args == ("bench-exact",) for kind, args, _ in calls)
    assert not any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)
    assert not any(kind == "uv" and "--suite" in args for kind, args, _ in calls)
    assert not any(kind == "uv" and "--scope" in args for kind, args, _ in calls)
    assert not any(kind == "git" and args[:1] == ("read-tree",) for kind, args, _ in calls)
    assert not any(kind == "git-stdin" for kind, _, _ in calls)
