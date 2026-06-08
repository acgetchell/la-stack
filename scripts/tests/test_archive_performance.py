"""Tests for archive_performance.py."""

from __future__ import annotations

import io
import subprocess
import tarfile
from pathlib import Path

import pytest

import archive_performance
from archive_performance import GenerationConfig, generate_and_promote_worktree_report, main, normalize_tag, parse_report_id, promote_report


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


def _write_baseline_archive(path: Path) -> None:
    fixture_dir = path.parent / "baseline-fixture"
    criterion_dir = fixture_dir / "criterion"
    criterion_dir.mkdir(parents=True)
    (criterion_dir / "placeholder.txt").write_text("baseline\n", encoding="utf-8")
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


def test_promote_report_archives_previous_and_updates_sorted_index(tmp_path) -> None:
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
    assert current.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")
    assert (archive_dir / "v0.4.1-vs-v0.4.0.md").read_text(encoding="utf-8") == _report("0.4.1", "v0.4.0")
    assert (archive_dir / "README.md").read_text(encoding="utf-8") == (
        "# Archived Performance Reports\n\n"
        "Older release-to-release benchmark comparisons are archived here.\n"
        "`docs/PERFORMANCE.md` contains the latest curated comparison.\n\n"
        "- [v0.3.1-vs-v0.3.0](v0.3.1-vs-v0.3.0.md)\n"
        "- [v0.4.1-vs-v0.4.0](v0.4.1-vs-v0.4.0.md)\n"
    )


def test_promote_report_is_idempotent_for_same_release_pair(tmp_path) -> None:
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


def test_promote_report_rejects_unexpected_release_pair(tmp_path) -> None:
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


def test_main_promotes_generated_report_to_docs_performance(tmp_path, capsys) -> None:
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
    assert current.read_text(encoding="utf-8") == generated
    assert (archive_dir / "v0.4.2-vs-v0.4.1.md").exists()
    assert "Current performance report: v0.4.3 vs v0.4.2" in capsys.readouterr().out


def test_main_reports_release_pair_mismatch_to_stderr(tmp_path, capsys) -> None:
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


def test_main_generates_report_in_temp_worktree(tmp_path, monkeypatch, capsys) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[tuple[str, tuple[str, ...], Path | None]] = []

    def fake_run_git(args, cwd=None, **kwargs):
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return type("Result", (), {"stdout": ""})()

    def fake_run_git_with_input(args, input_data, cwd=None, **kwargs):
        calls.append(("git-stdin", tuple(args), cwd))
        return type("Result", (), {"stdout": ""})()

    def fake_run_safe(command, args, cwd=None, **kwargs):
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            _write_baseline_archive(download_dir / "la-stack-v0.4.2-criterion-baseline.tar.gz")
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return type("Result", (), {"stdout": ""})()

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
    assert current.read_text(encoding="utf-8") == _report("0.4.3", "v0.4.2")
    assert "Generated benchmark report in a temporary worktree" in captured.out
    assert "target/bench-reports/performance.md" not in captured.out
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "v0.4.3" for kind, args, _ in calls)
    assert any(kind == "uv" and "--suite" in args and args[args.index("--suite") + 1] == "exact" for kind, args, _ in calls)
    assert not any(kind == "git" and args == ("diff", "--binary", "HEAD") for kind, args, _ in calls)
    assert not any(kind == "git-stdin" for kind, _, _ in calls)


def test_temp_worktree_is_removed_when_benchmark_command_fails(tmp_path, monkeypatch, capsys) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[tuple[str, tuple[str, ...], Path | None]] = []

    def fake_run_git(args, cwd=None, **kwargs):
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return type("Result", (), {"stdout": ""})()

    def fake_run_git_with_input(args, input_data, cwd=None, **kwargs):
        calls.append(("git-stdin", tuple(args), cwd))
        return type("Result", (), {"stdout": ""})()

    def fake_run_safe(command, args, cwd=None, **kwargs):
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            _write_baseline_archive(download_dir / "la-stack-v0.4.2-criterion-baseline.tar.gz")
        if command == "just" and args == ["bench-latest"]:
            raise subprocess.CalledProcessError(42, [command, *args], output="bench stdout", stderr="bench stderr")
        return type("Result", (), {"stdout": ""})()

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


def test_generate_report_rejects_unsafe_baseline_archive(tmp_path, monkeypatch, capsys) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[tuple[str, tuple[str, ...], Path | None]] = []

    def fake_run_git(args, cwd=None, **kwargs):
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return type("Result", (), {"stdout": ""})()

    def fake_run_git_with_input(args, input_data, cwd=None, **kwargs):
        calls.append(("git-stdin", tuple(args), cwd))
        return type("Result", (), {"stdout": ""})()

    def fake_run_safe(command, args, cwd=None, **kwargs):
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            _write_unsafe_baseline_archive(download_dir / "la-stack-v0.4.2-criterion-baseline.tar.gz")
        return type("Result", (), {"stdout": ""})()

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
    assert "refusing to extract unsafe archive member '../escape.txt'" in captured.err
    assert not (tmp_path / "escape.txt").exists()
    assert not current.exists()
    assert not any(kind in {"just", "uv"} for kind, _, _ in calls)
    assert any(kind == "git" and args[:3] == ("worktree", "remove", "--force") for kind, args, _ in calls)


def test_generate_report_fails_when_release_baseline_asset_missing(tmp_path, monkeypatch, capsys) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[tuple[str, tuple[str, ...], Path | None]] = []

    def fake_run_git(args, cwd=None, **kwargs):
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        return type("Result", (), {"stdout": ""})()

    def fake_run_git_with_input(args, input_data, cwd=None, **kwargs):
        calls.append(("git-stdin", tuple(args), cwd))
        return type("Result", (), {"stdout": ""})()

    def fake_run_safe(command, args, cwd=None, **kwargs):
        calls.append((command, tuple(args), cwd))
        return type("Result", (), {"stdout": ""})()

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
    assert "release baseline asset was not downloaded" in captured.err
    assert not current.exists()
    assert not any(kind in {"just", "uv"} for kind, _, _ in calls)
    assert any(kind == "git" and args[:3] == ("worktree", "remove", "--force") for kind, args, _ in calls)


def test_failed_atomic_replace_preserves_existing_report(tmp_path, monkeypatch) -> None:
    source = tmp_path / "performance-new.md"
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    original = _report("0.4.2", "v0.4.1")

    source.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
    current.parent.mkdir(parents=True)
    current.write_text(original, encoding="utf-8")

    def fail_replace(src, dst) -> None:
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


def test_generate_and_promote_uses_temp_worktree_and_current_diff(tmp_path, monkeypatch) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    current.parent.mkdir(parents=True)
    current.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
    calls: list[tuple[str, tuple[str, ...], Path | None]] = []

    def fake_run_git(args, cwd=None, **kwargs):
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_current_benchmark_tooling(worktree)
        if args == ["diff", "--binary", "HEAD"]:
            return type("Result", (), {"stdout": "diff --git a/README.md b/README.md\n"})()
        return type("Result", (), {"stdout": ""})()

    def fake_run_git_with_input(args, input_data, cwd=None, **kwargs):
        calls.append(("git-stdin", tuple(args), cwd))
        assert "diff --git" in input_data
        return type("Result", (), {"stdout": ""})()

    def fake_run_safe(command, args, cwd=None, **kwargs):
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            _write_baseline_archive(download_dir / "la-stack-v0.4.2-criterion-baseline.tar.gz")
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.3", "v0.4.2"), encoding="utf-8")
        return type("Result", (), {"stdout": ""})()

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
    assert current.read_text(encoding="utf-8") == _report("0.4.3", "v0.4.2")
    assert (archive_dir / "v0.4.2-vs-v0.4.1.md").exists()
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "HEAD" for kind, args, _ in calls)
    assert any(kind == "git-stdin" and args == ("apply", "--binary") for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)
    assert any(kind == "git" and args[:3] == ("worktree", "remove", "--force") for kind, args, _ in calls)


def test_generate_and_promote_legacy_published_tag_uses_legacy_commands(tmp_path, monkeypatch) -> None:
    current = tmp_path / "docs" / "PERFORMANCE.md"
    archive_dir = tmp_path / "docs" / "archive" / "performance"
    calls: list[tuple[str, tuple[str, ...], Path | None]] = []

    def fake_run_git(args, cwd=None, **kwargs):
        calls.append(("git", tuple(args), cwd))
        if args[:3] == ["worktree", "add", "--detach"]:
            worktree = Path(args[3])
            worktree.mkdir(parents=True)
            _write_legacy_benchmark_tooling(worktree)
        return type("Result", (), {"stdout": ""})()

    def fake_run_git_with_input(args, input_data, cwd=None, **kwargs):
        calls.append(("git-stdin", tuple(args), cwd))
        return type("Result", (), {"stdout": ""})()

    def fake_run_safe(command, args, cwd=None, **kwargs):
        calls.append((command, tuple(args), cwd))
        if command == "gh":
            download_dir = Path(args[args.index("--dir") + 1])
            _write_baseline_archive(download_dir / "la-stack-v0.4.1-criterion-baseline.tar.gz")
        if command == "uv":
            output = Path(args[args.index("--output") + 1])
            output.write_text(_report("0.4.2", "v0.4.1"), encoding="utf-8")
        return type("Result", (), {"stdout": ""})()

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
    assert current.read_text(encoding="utf-8") == _report("0.4.2", "v0.4.1")
    assert any(kind == "git" and args[:3] == ("worktree", "add", "--detach") and args[4] == "v0.4.2" for kind, args, _ in calls)
    assert any(kind == "just" and args == ("bench-exact",) for kind, args, _ in calls)
    assert not any(kind == "just" and args == ("bench-latest",) for kind, args, _ in calls)
    assert not any(kind == "uv" and "--suite" in args for kind, args, _ in calls)
    assert not any(kind == "uv" and "--scope" in args for kind, args, _ in calls)
    assert not any(kind == "git" and args == ("diff", "--binary", "HEAD") for kind, args, _ in calls)
    assert not any(kind == "git-stdin" for kind, _, _ in calls)
