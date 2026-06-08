"""Tests for archive_performance.py."""

from __future__ import annotations

import io
import subprocess
import tarfile
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import pytest

import archive_performance
from archive_performance import GenerationConfig, generate_and_promote_worktree_report, main, normalize_tag, parse_report_id, promote_report

if TYPE_CHECKING:
    from collections.abc import Sequence

type RunnerCall = tuple[str, tuple[str, ...], Path | None]


def _result(stdout: str = "") -> SimpleNamespace:
    return SimpleNamespace(stdout=stdout)


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


def _write_baseline_archive(path: Path) -> None:
    tag = path.name.removeprefix("la-stack-").removesuffix("-criterion-baseline.tar.gz")
    fixture_dir = path.parent / f"baseline-fixture-{tag}"
    criterion_dir = fixture_dir / "criterion"
    criterion_dir.mkdir(parents=True)
    (criterion_dir / "placeholder.txt").write_text("baseline\n", encoding="utf-8")
    sample_dir = criterion_dir / "exact_d2" / "det_exact" / tag
    sample_dir.mkdir(parents=True)
    (sample_dir / "estimates.json").write_text('{"median":{"point_estimate":1.0}}\n', encoding="utf-8")
    with tarfile.open(path, "w:gz") as tar:
        tar.add(criterion_dir, arcname="criterion")


def _write_unsafe_baseline_archive(path: Path) -> None:
    payload = b"escape\n"
    info = tarfile.TarInfo("../escape.txt")
    info.size = len(payload)
    with tarfile.open(path, "w:gz") as tar:
        tar.addfile(info, io.BytesIO(payload))


def _write_current_benchmark_tooling(worktree: Path) -> None:
    (worktree / "scripts").mkdir(parents=True, exist_ok=True)
    (worktree / "justfile").write_text("bench-latest: bench-vs-linalg-la-stack bench-exact\n", encoding="utf-8")
    (worktree / "scripts" / "bench_compare.py").write_text('parser.add_argument("--suite")\nparser.add_argument("--scope")\n', encoding="utf-8")


def _write_legacy_benchmark_tooling(worktree: Path) -> None:
    (worktree / "scripts").mkdir(parents=True, exist_ok=True)
    (worktree / "justfile").write_text("bench-exact:\n", encoding="utf-8")
    (worktree / "scripts" / "bench_compare.py").write_text('parser.add_argument("--output")\n', encoding="utf-8")


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


def test_main_promotes_generated_report_to_docs_performance(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    source = tmp_path / "target" / "bench-reports" / "performance.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    generated = _report("0.4.3", "v0.4.2")

    source.parent.mkdir(parents=True)
    source.write_text(generated, encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")

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
    assert "Current performance report: v0.4.3 vs v0.4.2" in capsys.readouterr().out


def test_main_reports_release_pair_mismatch_to_stderr(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    source = tmp_path / "target" / "bench-reports" / "performance.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    source.parent.mkdir(parents=True)
    source.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")

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
        if command == "just" and args == ["bench-save-baseline", "v0.4.2"]:
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
    assert any(kind == "uv" and "--suite" in args and args[args.index("--suite") + 1] == "exact" for kind, args, _ in calls)
    assert not any(kind == "git" and args == ("diff", "--binary", "HEAD") for kind, args, _ in calls)
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


def test_generate_report_generates_release_baseline_locally(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
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
            criterion_dir = cwd / "target" / "criterion"
            criterion_dir.mkdir(parents=True)
            (criterion_dir / "baseline.txt").write_text("baseline\n", encoding="utf-8")
        if command == "just" and args == ["bench-latest"]:
            assert kwargs["env"]["RUSTUP_TOOLCHAIN"] == "1.96.0"
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
    assert sum(1 for kind, args, _ in calls if kind == "git" and args[:3] == ("worktree", "remove", "--force")) == 2


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
        if args == ["diff", "--binary", "HEAD"]:
            return _result("diff --git a/README.md b/README.md\n")
        return _result()

    def fake_run_git_with_input(args: Sequence[str], input_data: str, cwd: Path | None = None, **kwargs: Any) -> SimpleNamespace:
        calls.append(("git-stdin", tuple(args), cwd))
        assert "diff --git" in input_data
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
    assert any(kind == "just" and args == ("bench-exact",) for kind, args, _ in calls)
    assert not any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)
    assert not any(kind == "uv" and "--suite" in args for kind, args, _ in calls)
    assert not any(kind == "uv" and "--scope" in args for kind, args, _ in calls)
    assert not any(kind == "git" and args == ("diff", "--binary", "HEAD") for kind, args, _ in calls)
    assert not any(kind == "git-stdin" for kind, _, _ in calls)
