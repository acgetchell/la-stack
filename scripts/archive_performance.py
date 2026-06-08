#!/usr/bin/env -S uv run
"""Promote a benchmark report into docs/PERFORMANCE.md and archive the old one.

Release performance docs have two different lifetimes:

  - ``target/bench-reports/performance.md`` is local scratch output for the
    current machine and branch.
  - ``docs/PERFORMANCE.md`` is the latest curated release-to-release comparison.
  - ``docs/archive/performance/*.md`` stores older curated comparisons.

This script copies a freshly generated local report into ``docs/PERFORMANCE.md``
and archives the previous committed report under a filename derived from the
report metadata, such as ``v0.4.2-vs-v0.4.1.md``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

from subprocess_utils import run_git_command, run_git_command_with_input, run_safe_command

_VERSION_RE = re.compile(r"^\*\*la-stack\*\* v(?P<version>[^\s`]+)", re.MULTILINE)
_BASELINE_RE = re.compile(r"^Comparison against baseline \*\*(?P<baseline>[^*]+)\*\*:", re.MULTILINE)
_SEMVER_IDENTIFIER_RE = r"(?:0|[1-9][0-9]*|[0-9A-Za-z-]*[A-Za-z-][0-9A-Za-z-]*)"
_TAG_RE = re.compile(
    rf"^v?(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
    rf"(?:-{_SEMVER_IDENTIFIER_RE}(?:\.{_SEMVER_IDENTIFIER_RE})*)?"
    r"(?:\+[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?$"
)
_SEMVER_PARTS_RE = re.compile(r"^v?(?P<major>0|[1-9][0-9]*)\.(?P<minor>0|[1-9][0-9]*)\.(?P<patch>0|[1-9][0-9]*)$")

_DEFAULT_SOURCE = "target/bench-reports/performance.md"
_DEFAULT_CURRENT = "docs/PERFORMANCE.md"
_DEFAULT_ARCHIVE_DIR = "docs/archive/performance"
_DEFAULT_SUITE = "all"
_DEFAULT_SCOPE = "release-signal"
_BENCH_TIMEOUT_SECONDS = 7200
_COMMAND_TIMEOUT_SECONDS = 600
_HOW_TO_UPDATE_RE = re.compile(r"(?ms)^## How to Update\n.*\Z")


@dataclass(frozen=True)
class ReportId:
    """Release-pair identity parsed from a benchmark report."""

    current_tag: str
    baseline_tag: str

    @property
    def archive_name(self) -> str:
        """Return the canonical archive filename for this report."""
        return f"{self.current_tag}-vs-{self.baseline_tag}.md"


@dataclass(frozen=True)
class GenerationConfig:
    """Configuration for benchmark report generation in a temp worktree."""

    repo_root: Path
    current_tag: str
    baseline_tag: str
    worktree_ref: str
    suite: str = _DEFAULT_SUITE
    scope: str = _DEFAULT_SCOPE
    apply_current_diff: bool = True


@dataclass(frozen=True)
class ResolvedArchiveRequest:
    """Release pair and worktree ref resolved from CLI arguments."""

    current_tag: str
    baseline_tag: str
    worktree_ref: str
    fetch_tags: bool = False


@dataclass(frozen=True)
class PublishedRelease:
    """Stable GitHub release metadata used to infer release pairs."""

    tag: str
    published_at: str


def normalize_tag(tag: str) -> str:
    """Return *tag* with a leading ``v`` and no surrounding whitespace."""
    normalized = tag.strip()
    if not normalized:
        msg = "tag must not be empty"
        raise ValueError(msg)
    if not normalized.startswith("v"):
        normalized = f"v{normalized}"
    if not _TAG_RE.fullmatch(normalized):
        msg = f"expected a semver tag like v0.4.2, got {tag!r}"
        raise ValueError(msg)
    return normalized


def parse_report_id(text: str) -> ReportId:
    """Parse the current version and baseline tag from a benchmark report."""
    version_match = _VERSION_RE.search(text)
    if version_match is None:
        msg = "could not find la-stack version line in benchmark report"
        raise ValueError(msg)

    baseline_match = _BASELINE_RE.search(text)
    if baseline_match is None:
        msg = "could not find comparison baseline line in benchmark report"
        raise ValueError(msg)

    return ReportId(
        current_tag=normalize_tag(version_match.group("version")),
        baseline_tag=normalize_tag(baseline_match.group("baseline")),
    )


def _semver_sort_key(tag: str) -> tuple[int, int, int]:
    match = _SEMVER_PARTS_RE.fullmatch(normalize_tag(tag))
    if match is None:
        msg = f"expected a stable semver tag like v0.4.2, got {tag!r}"
        raise ValueError(msg)
    return (int(match.group("major")), int(match.group("minor")), int(match.group("patch")))


def _stable_published_releases(releases: object) -> list[PublishedRelease]:
    if not isinstance(releases, list):
        msg = "expected GitHub release list to be a JSON array"
        raise TypeError(msg)

    stable_releases: dict[str, PublishedRelease] = {}
    for release in releases:
        if not isinstance(release, Mapping):
            continue
        release = cast("Mapping[str, Any]", release)
        if release.get("isDraft") or release.get("isPrerelease"):
            continue
        tag_name = release.get("tagName")
        published_at = release.get("publishedAt")
        if not isinstance(tag_name, str) or not isinstance(published_at, str) or not published_at:
            continue
        try:
            normalized = normalize_tag(tag_name)
            _semver_sort_key(normalized)
        except ValueError:
            continue
        stable_releases[normalized] = PublishedRelease(tag=normalized, published_at=published_at)

    return list(stable_releases.values())


def _published_release_pair(repo_root: Path) -> ReportId:
    command = [
        "release",
        "list",
        "--json",
        "tagName,isDraft,isPrerelease,publishedAt",
        "--limit",
        "100",
    ]
    try:
        result = run_safe_command(
            "gh",
            command,
            cwd=repo_root,
            timeout=_COMMAND_TIMEOUT_SECONDS,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_command_failure(["gh", *command], exc)) from exc
    try:
        releases = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        msg = "could not parse GitHub release list JSON"
        raise RuntimeError(msg) from exc
    stable_releases = _stable_published_releases(releases)
    if len(stable_releases) < 2:
        msg = "expected at least two published stable semver releases"
        raise RuntimeError(msg)

    current = max(stable_releases, key=lambda release: release.published_at)
    current_key = _semver_sort_key(current.tag)
    previous_tags = sorted(
        (release.tag for release in stable_releases if _semver_sort_key(release.tag) < current_key),
        key=_semver_sort_key,
    )
    if not previous_tags:
        msg = f"could not find a previous stable semver release before {current.tag}"
        raise RuntimeError(msg)
    return ReportId(current_tag=current.tag, baseline_tag=previous_tags[-1])


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _how_to_update_section() -> str:
    lines = [
        "## How to Update",
        "",
        "Release performance docs are generated in isolated temporary worktrees:",
        "",
        "```bash",
        "# Release PR: update docs/PERFORMANCE.md and archive the previous report",
        "just performance-release <current-tag> <previous-tag>",
        "",
        "# Historical published comparison",
        "just performance-archive-published",
        "",
        "# Explicit historical repair",
        "just performance-archive-published <current-tag> <previous-tag>",
        "```",
        "",
        "For local scratch comparisons, use `just bench-latest` and `just bench-compare`.",
        "Those write `target/bench-reports/performance.md`.",
        "",
        "See `docs/BENCHMARKING.md` for the full comparison workflow.",
        "",
    ]
    return "\n".join(lines)


def _normalize_how_to_update(text: str) -> str:
    section = _how_to_update_section()
    if _HOW_TO_UPDATE_RE.search(text):
        return _HOW_TO_UPDATE_RE.sub(section, text)
    return f"{text.rstrip()}\n\n{section}"


def _replace_file(src: Path, dst: Path) -> None:
    src.replace(dst)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as tmp:
            tmp_path = Path(tmp.name)
            tmp.write(text)
            tmp.flush()
            os.fsync(tmp.fileno())
        _replace_file(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def _archive_readme(archive_dir: Path) -> str:
    reports = sorted(path.name for path in archive_dir.glob("*.md") if path.name != "README.md")
    lines = [
        "# Archived Performance Reports",
        "",
        "Older release-to-release benchmark comparisons are archived here.",
        "`docs/PERFORMANCE.md` contains the latest curated comparison.",
        "",
    ]
    if reports:
        lines.extend(f"- [{name.removesuffix('.md')}]({name})" for name in reports)
    else:
        lines.append("- No archived performance reports yet.")
    return "\n".join(lines) + "\n"


def update_archive_index(archive_dir: Path) -> None:
    """Write a lexicographically sorted archive index."""
    _write_text(archive_dir / "README.md", _archive_readme(archive_dir))


def _format_command_failure(command: list[str], exc: subprocess.CalledProcessError) -> str:
    parts = [f"command failed ({exc.returncode}): {' '.join(command)}"]
    if exc.stdout:
        parts.append(f"stdout:\n{exc.stdout.strip()}")
    if exc.stderr:
        parts.append(f"stderr:\n{exc.stderr.strip()}")
    return "\n".join(parts)


def _run_git(args: list[str], *, cwd: Path, timeout: int = _COMMAND_TIMEOUT_SECONDS) -> None:
    try:
        run_git_command(args, cwd=cwd, timeout=timeout)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_command_failure(["git", *args], exc)) from exc


def _fetch_release_tags(*, repo_root: Path, tags: list[str]) -> None:
    refspecs = [f"refs/tags/{tag}:refs/tags/{tag}" for tag in tags]
    _run_git(["fetch", "origin", *refspecs], cwd=repo_root)


def _run_tool(command: str, args: list[str], *, cwd: Path, timeout: int = _COMMAND_TIMEOUT_SECONDS, env: dict[str, str] | None = None) -> None:
    try:
        run_safe_command(command, args, cwd=cwd, timeout=timeout, env=env)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_command_failure([command, *args], exc)) from exc


def _current_rust_toolchain(repo_root: Path) -> str | None:
    rust_toolchain = repo_root / "rust-toolchain.toml"
    if not rust_toolchain.exists():
        return None
    data = tomllib.loads(_read_text(rust_toolchain))
    toolchain = data.get("toolchain")
    if not isinstance(toolchain, dict):
        return None
    channel = toolchain.get("channel")
    return channel if isinstance(channel, str) else None


def _benchmark_env(repo_root: Path) -> dict[str, str] | None:
    if "RUSTUP_TOOLCHAIN" in os.environ:
        return None
    toolchain = _current_rust_toolchain(repo_root)
    if toolchain is None:
        return None
    env = os.environ.copy()
    env["RUSTUP_TOOLCHAIN"] = toolchain
    return env


def _safe_extract_tar(archive: Path, target_dir: Path) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    target_root = target_dir.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            member_path = (target_dir / member.name).resolve()
            if not member_path.is_relative_to(target_root):
                msg = f"refusing to extract unsafe archive member {member.name!r}"
                raise ValueError(msg)
        tar.extractall(target_dir, filter="data")


def _download_release_baseline(*, baseline_tag: str, download_dir: Path, repo_root: Path) -> Path:
    artifact = download_dir / f"la-stack-{baseline_tag}-criterion-baseline.tar.gz"
    _run_tool(
        "gh",
        [
            "release",
            "download",
            baseline_tag,
            "--pattern",
            artifact.name,
            "--dir",
            str(download_dir),
        ],
        cwd=repo_root,
    )
    if not artifact.exists():
        msg = f"release baseline asset was not downloaded: {artifact}"
        raise FileNotFoundError(msg)
    return artifact


def _generate_release_baseline(*, baseline_tag: str, repo_root: Path, target_worktree: Path, tmp_dir: Path) -> None:
    baseline_worktree = tmp_dir / "baseline-worktree"
    _run_git(["worktree", "add", "--detach", str(baseline_worktree), baseline_tag], cwd=repo_root)
    try:
        _run_tool("just", ["bench-save-baseline", baseline_tag], cwd=baseline_worktree, timeout=_BENCH_TIMEOUT_SECONDS, env=_benchmark_env(repo_root))
        baseline_criterion = baseline_worktree / "target" / "criterion"
        if not baseline_criterion.is_dir():
            msg = f"generated baseline Criterion results were not found: {baseline_criterion}"
            raise FileNotFoundError(msg)
        target_criterion = target_worktree / "target" / "criterion"
        target_criterion.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(baseline_criterion, target_criterion, dirs_exist_ok=True)
    finally:
        try:
            _run_git(["worktree", "remove", "--force", str(baseline_worktree)], cwd=repo_root)
        except RuntimeError as exc:
            print(f"archive-performance: failed to remove baseline worktree: {exc}", file=sys.stderr)


def _prepare_release_baseline(*, baseline_tag: str, repo_root: Path, target_worktree: Path, tmp_dir: Path) -> None:
    try:
        baseline_archive = _download_release_baseline(
            baseline_tag=baseline_tag,
            download_dir=tmp_dir,
            repo_root=repo_root,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"archive-performance: release baseline asset unavailable; generating {baseline_tag} locally ({exc})", file=sys.stderr)
        _generate_release_baseline(
            baseline_tag=baseline_tag,
            repo_root=repo_root,
            target_worktree=target_worktree,
            tmp_dir=tmp_dir,
        )
    else:
        _safe_extract_tar(baseline_archive, target_worktree / "target")


def _apply_current_diff_to_worktree(*, repo_root: Path, worktree: Path) -> None:
    diff = run_git_command(["diff", "--binary", "HEAD"], cwd=repo_root).stdout
    if diff.strip():
        try:
            run_git_command_with_input(["apply", "--binary"], diff, cwd=worktree)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(_format_command_failure(["git", "apply", "--binary"], exc)) from exc


def _has_current_release_signal_tooling(worktree: Path) -> bool:
    justfile = worktree / "justfile"
    bench_compare = worktree / "scripts" / "bench_compare.py"
    if not justfile.exists() or not bench_compare.exists():
        return False

    justfile_text = _read_text(justfile)
    bench_compare_text = _read_text(bench_compare)
    return re.search(r"(?m)^bench-latest(?:[ :]|$)", justfile_text) is not None and '"--suite"' in bench_compare_text and '"--scope"' in bench_compare_text


def _run_benchmarks_and_render_report(*, worktree: Path, report: Path, config: GenerationConfig) -> None:
    benchmark_env = _benchmark_env(config.repo_root)
    if _has_current_release_signal_tooling(worktree):
        _run_tool("just", ["bench-latest"], cwd=worktree, timeout=_BENCH_TIMEOUT_SECONDS, env=benchmark_env)
        _run_tool(
            "uv",
            [
                "run",
                "bench-compare",
                config.baseline_tag,
                "--suite",
                config.suite,
                "--scope",
                config.scope,
                "--output",
                str(report),
            ],
            cwd=worktree,
            timeout=_COMMAND_TIMEOUT_SECONDS,
        )
    else:
        _run_tool("just", ["bench-exact"], cwd=worktree, timeout=_BENCH_TIMEOUT_SECONDS, env=benchmark_env)
        _run_tool(
            "uv",
            [
                "run",
                "bench-compare",
                config.baseline_tag,
                "--output",
                str(report),
            ],
            cwd=worktree,
            timeout=_COMMAND_TIMEOUT_SECONDS,
        )


def _generate_report_in_temp_worktree(
    *,
    config: GenerationConfig,
) -> str:
    with tempfile.TemporaryDirectory(prefix="la-stack-performance-") as tmp:
        tmp_dir = Path(tmp)
        worktree = tmp_dir / "worktree"
        report = tmp_dir / f"{config.current_tag}-vs-{config.baseline_tag}.md"

        _run_git(["worktree", "add", "--detach", str(worktree), config.worktree_ref], cwd=config.repo_root)
        try:
            if config.apply_current_diff:
                _apply_current_diff_to_worktree(repo_root=config.repo_root, worktree=worktree)
            _prepare_release_baseline(
                baseline_tag=config.baseline_tag,
                repo_root=config.repo_root,
                target_worktree=worktree,
                tmp_dir=tmp_dir,
            )
            _run_benchmarks_and_render_report(worktree=worktree, report=report, config=config)
            return _read_text(report)
        finally:
            try:
                _run_git(["worktree", "remove", "--force", str(worktree)], cwd=config.repo_root)
            except RuntimeError as exc:
                print(f"archive-performance: failed to remove temporary worktree: {exc}", file=sys.stderr)


def promote_report(
    *,
    source: Path,
    current: Path,
    archive_dir: Path,
    expected_current_tag: str,
    expected_baseline_tag: str,
) -> ReportId:
    """Archive the old committed report and promote *source* as the current one."""
    source_text = _normalize_how_to_update(_read_text(source))
    source_id = parse_report_id(source_text)
    expected_source_id = ReportId(
        current_tag=normalize_tag(expected_current_tag),
        baseline_tag=normalize_tag(expected_baseline_tag),
    )
    if source_id != expected_source_id:
        msg = (
            "benchmark report does not match requested release pair: "
            f"found {source_id.current_tag} vs {source_id.baseline_tag}, "
            f"expected {expected_source_id.current_tag} vs {expected_source_id.baseline_tag}"
        )
        raise ValueError(msg)

    if current.exists():
        current_text = _normalize_how_to_update(_read_text(current))
        current_id = parse_report_id(current_text)
        if current_id != source_id:
            archive_path = archive_dir / current_id.archive_name
            if not archive_path.exists():
                _write_text(archive_path, current_text)

    _write_text(current, source_text)
    update_archive_index(archive_dir)
    return source_id


def generate_and_promote_worktree_report(
    *,
    current: Path,
    archive_dir: Path,
    config: GenerationConfig,
) -> ReportId:
    """Generate a comparison in a temp worktree, then promote it."""
    current_tag = normalize_tag(config.current_tag)
    baseline_tag = normalize_tag(config.baseline_tag)
    config = GenerationConfig(
        repo_root=config.repo_root,
        current_tag=current_tag,
        baseline_tag=baseline_tag,
        worktree_ref=config.worktree_ref,
        suite=config.suite,
        scope=config.scope,
        apply_current_diff=config.apply_current_diff,
    )
    report_text = _generate_report_in_temp_worktree(
        config=config,
    )
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".md", delete=False) as tmp:
        source = Path(tmp.name)
        tmp.write(report_text)
    try:
        return promote_report(
            source=source,
            current=current,
            archive_dir=archive_dir,
            expected_current_tag=current_tag,
            expected_baseline_tag=baseline_tag,
        )
    finally:
        if source.exists():
            source.unlink()


def resolve_archive_request(
    *,
    current_tag: str | None,
    baseline_tag: str | None,
    published_latest: bool,
    worktree_ref: str,
    repo_root: Path,
) -> ResolvedArchiveRequest:
    """Resolve explicit or latest-published release arguments."""
    if published_latest:
        if current_tag is not None or baseline_tag is not None:
            msg = "do not pass current_tag or baseline_tag with --published-latest"
            raise ValueError(msg)
        published_pair = _published_release_pair(repo_root)
        resolved_worktree_ref = published_pair.current_tag if worktree_ref == "HEAD" else worktree_ref
        return ResolvedArchiveRequest(
            current_tag=published_pair.current_tag,
            baseline_tag=published_pair.baseline_tag,
            worktree_ref=resolved_worktree_ref,
            fetch_tags=True,
        )

    if current_tag is None or baseline_tag is None:
        msg = "current_tag and baseline_tag are required unless --published-latest is used"
        raise ValueError(msg)
    return ResolvedArchiveRequest(
        current_tag=current_tag,
        baseline_tag=baseline_tag,
        worktree_ref=worktree_ref,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Promote a benchmark comparison into docs/PERFORMANCE.md and archive the previous report.",
    )
    parser.add_argument("current_tag", nargs="?", help="Release tag for the new report, e.g. v0.4.3")
    parser.add_argument("baseline_tag", nargs="?", help="Previous release tag used as the comparison baseline, e.g. v0.4.2")
    parser.add_argument(
        "--source",
        default=_DEFAULT_SOURCE,
        help=f"Generated benchmark report to promote (default: {_DEFAULT_SOURCE})",
    )
    parser.add_argument(
        "--current",
        default=_DEFAULT_CURRENT,
        help=f"Committed performance report path (default: {_DEFAULT_CURRENT})",
    )
    parser.add_argument(
        "--archive-dir",
        default=_DEFAULT_ARCHIVE_DIR,
        help=f"Archive directory for older reports (default: {_DEFAULT_ARCHIVE_DIR})",
    )
    parser.add_argument(
        "--generate-in-temp-worktree",
        action="store_true",
        help="Generate the comparison in a temporary detached worktree before promoting it.",
    )
    parser.add_argument(
        "--published-latest",
        action="store_true",
        help="Infer the latest stable published GitHub release and its previous stable release.",
    )
    parser.add_argument(
        "--worktree-ref",
        default="HEAD",
        help="Git ref to check out in the temporary worktree (default: HEAD).",
    )
    parser.add_argument(
        "--no-apply-current-diff",
        action="store_true",
        help="Do not apply the current checkout's tracked diff to the temporary worktree.",
    )
    parser.add_argument(
        "--suite",
        default=_DEFAULT_SUITE,
        help=f"Benchmark suite for --generate-in-temp-worktree (default: {_DEFAULT_SUITE})",
    )
    parser.add_argument(
        "--scope",
        default=_DEFAULT_SCOPE,
        help=f"Comparison scope for --generate-in-temp-worktree (default: {_DEFAULT_SCOPE})",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    args = build_parser().parse_args(argv)
    root = Path.cwd()
    source = Path(args.source)
    current = Path(args.current)
    archive_dir = Path(args.archive_dir)
    if not source.is_absolute():
        source = root / source
    if not current.is_absolute():
        current = root / current
    if not archive_dir.is_absolute():
        archive_dir = root / archive_dir

    try:
        request = resolve_archive_request(
            current_tag=args.current_tag,
            baseline_tag=args.baseline_tag,
            published_latest=args.published_latest,
            worktree_ref=args.worktree_ref,
            repo_root=root,
        )
        if request.fetch_tags:
            _fetch_release_tags(repo_root=root, tags=[request.current_tag, request.baseline_tag])

        if args.generate_in_temp_worktree:
            report_id = generate_and_promote_worktree_report(
                current=current,
                archive_dir=archive_dir,
                config=GenerationConfig(
                    repo_root=root,
                    current_tag=request.current_tag,
                    baseline_tag=request.baseline_tag,
                    worktree_ref=request.worktree_ref,
                    suite=args.suite,
                    scope=args.scope,
                    apply_current_diff=not args.no_apply_current_diff,
                ),
            )
        else:
            report_id = promote_report(
                source=source,
                current=current,
                archive_dir=archive_dir,
                expected_current_tag=request.current_tag,
                expected_baseline_tag=request.baseline_tag,
            )
    except Exception as exc:
        print(f"archive-performance: {exc}", file=sys.stderr)
        return 1

    if args.generate_in_temp_worktree:
        print(f"Generated benchmark report in a temporary worktree and promoted it to {current}")
    else:
        print(f"Promoted {source} to {current}")
    print(f"Current performance report: {report_id.current_tag} vs {report_id.baseline_tag}")
    print(f"Archive directory: {archive_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
