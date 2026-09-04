#!/usr/bin/env -S uv run --locked
"""Generate, retain, render, and promote benchmark comparison reports.

Release performance docs have two different lifetimes:

  - ``target/bench-reports/performance.md`` is local scratch output for the
    current machine and branch.
  - ``target/bench-reports/performance.csv`` and the adjacent provenance JSON
    are validated, reproducible performance-comparison inputs.
  - ``docs/PERFORMANCE.md`` is the latest curated release-to-release comparison.
  - ``docs/archive/performance/*.md`` stores older curated comparisons.

This script renders from a validated artifact reload, copies the result into
``docs/PERFORMANCE.md``, and archives the previous committed report under a
filename derived from the report metadata, such as ``v0.4.2-vs-v0.4.1.md``.
"""

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal, cast

from bench_compare import HOW_TO_UPDATE_SECTION, render_release_artifacts
from benchmark_contract import benchmark_contract_digest
from performance_artifacts import (
    NO_API_COMPATIBILITY,
    ArtifactPaths,
    ensure_distinct_paths,
    load_bundle,
    publish_bundle,
    resolve_shared_harness_compatibility,
)
from subprocess_utils import (
    ExecutableNotFoundError,
    cpu_description,
    format_exception_diagnostics,
    run_git_command,
    run_git_command_with_input,
    run_safe_command,
)

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
_DEFAULT_ARTIFACT_CSV = "target/bench-reports/performance.csv"
_DEFAULT_ARTIFACT_PROVENANCE = "target/bench-reports/performance.provenance.json"
_DEFAULT_CURRENT = "docs/PERFORMANCE.md"
_DEFAULT_ARCHIVE_DIR = "docs/archive/performance"
_DEFAULT_SUITE = "all"
_DEFAULT_SCOPE = "release-signal"
_SUPPORTED_SUITES = ("all", "exact", "vs_linalg")
_SUPPORTED_SCOPES = ("release-signal", "all-benches")
_BENCH_TIMEOUT_SECONDS = 7200
_COMMAND_TIMEOUT_SECONDS = 600
_HOW_TO_UPDATE_RE = re.compile(r"(?ms)^## How to Update\n.*\Z")
_BENCHMARK_HARNESS_DIRS = ("benches",)
_BENCHMARK_HARNESS_FILES = (
    ".config/nextest.toml",
    "Cargo.toml",
    "Cargo.lock",
    "rust-toolchain.toml",
    "justfile",
    "tests/exact_bench_config.rs",
    "tests/vs_linalg_inputs.rs",
)
_BENCHMARK_HARNESS_METADATA = ".la-stack-benchmark-harness.json"
_BENCHMARK_INPUT_GATE = ("just", "test-bench-inputs")
_COMPARISON_LINT_CAP = "--cap-lints=warn"
type BaselineSource = Literal["local", "github-assets"]
type BenchmarkSuite = Literal["all", "exact", "vs_linalg"]
type ComparisonScope = Literal["release-signal", "all-benches"]


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
class PromotionRequest:
    """Destinations and expected identity for one report promotion."""

    current: Path
    archive_dir: Path
    expected: ReportId
    source_path: Path | None = None
    output: Path | None = None
    reserved_paths: Mapping[str, Path] | None = None


@dataclass(frozen=True)
class GenerationConfig:
    """Configuration for benchmark report generation in a temp worktree."""

    repo_root: Path
    current_tag: str
    baseline_tag: str
    worktree_ref: str
    suite: BenchmarkSuite = "all"
    scope: ComparisonScope = "release-signal"
    apply_current_diff: bool = True
    baseline_source: BaselineSource = "local"

    def __post_init__(self) -> None:
        """Reject unsupported benchmark selections at construction time."""
        if self.suite not in _SUPPORTED_SUITES:
            msg = f"unsupported benchmark suite: {self.suite}"
            raise ValueError(msg)
        if self.scope not in _SUPPORTED_SCOPES:
            msg = f"unsupported comparison scope: {self.scope}"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ToolRunOptions:
    """Execution controls for one repository support command."""

    timeout: int = _COMMAND_TIMEOUT_SECONDS
    env: dict[str, str] | None = None
    stream_output: bool = False


@dataclass(frozen=True)
class ResolvedArchiveRequest:
    """Release pair and worktree ref resolved from CLI arguments."""

    current_tag: str
    baseline_tag: str
    worktree_ref: str
    tags_to_fetch: tuple[str, ...] = ()


@dataclass(frozen=True)
class ArchiveRequestOptions:
    """CLI options used to resolve release tags."""

    current_tag: str | None
    baseline_tag: str | None
    published_latest: bool
    infer_release: bool
    current_vs_latest: bool
    worktree_ref: str
    repo_root: Path


@dataclass(frozen=True)
class ArchivePaths:
    """Filesystem paths used by the archive CLI."""

    source: Path
    current: Path
    output: Path
    archive_dir: Path
    artifacts: ArtifactPaths


@dataclass(frozen=True)
class GeneratedReport:
    """A generated benchmark report that remains valid after worktree cleanup."""

    text: str


@dataclass(frozen=True)
class ArchiveResult:
    """Result and destination metadata for a completed archive operation."""

    report_id: ReportId
    action: Literal["output", "promote-generated", "promote-source", "promote-artifacts"]


@dataclass(frozen=True)
class PublishedRelease:
    """Stable GitHub release metadata used to infer release pairs."""

    tag: str
    published_at: datetime


@dataclass(frozen=True)
class BaselineRun:
    """Validated baseline run details needed for final provenance."""

    commit: str
    command: tuple[str, ...]
    harness_sha256: str
    git_clean: bool
    source_state_sha256: str
    api_compatibility: str | None
    shared_harness_rational_inputs: bool = False


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
    version_matches = list(_VERSION_RE.finditer(text))
    if not version_matches:
        msg = "could not find la-stack version line in benchmark report"
        raise ValueError(msg)
    if len(version_matches) != 1:
        msg = f"expected exactly one la-stack version line in benchmark report, found {len(version_matches)}"
        raise ValueError(msg)

    baseline_matches = list(_BASELINE_RE.finditer(text))
    if not baseline_matches:
        msg = "could not find comparison baseline line in benchmark report"
        raise ValueError(msg)
    if len(baseline_matches) != 1:
        msg = f"expected exactly one comparison baseline line in benchmark report, found {len(baseline_matches)}"
        raise ValueError(msg)

    version_match = version_matches[0]
    baseline_match = baseline_matches[0]

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


def _parse_published_at(value: object, *, release_index: int) -> datetime:
    """Parse one GitHub publication timestamp as an aware UTC datetime."""
    if not isinstance(value, str) or not value:
        msg = f"GitHub release entry {release_index} has invalid publishedAt: {value!r}"
        raise TypeError(msg)
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as exc:
        msg = f"GitHub release entry {release_index} has invalid publishedAt: {value!r}"
        raise ValueError(msg) from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        msg = f"GitHub release entry {release_index} publishedAt must include a UTC offset: {value!r}"
        raise ValueError(msg)
    return parsed.astimezone(UTC)


def _stable_published_releases(releases: object) -> list[PublishedRelease]:
    if not isinstance(releases, list):
        msg = "expected GitHub release list to be a JSON array"
        raise TypeError(msg)

    stable_releases: dict[str, PublishedRelease] = {}
    for index, release in enumerate(releases):
        if not isinstance(release, Mapping):
            msg = f"GitHub release entry {index} is not a JSON object"
            raise TypeError(msg)
        release = cast("Mapping[str, Any]", release)
        is_draft = release.get("isDraft")
        is_prerelease = release.get("isPrerelease")
        if not isinstance(is_draft, bool) or not isinstance(is_prerelease, bool):
            msg = f"GitHub release entry {index} must contain boolean isDraft and isPrerelease fields"
            raise TypeError(msg)
        if is_draft or is_prerelease:
            continue
        tag_name = release.get("tagName")
        if not isinstance(tag_name, str) or not tag_name:
            msg = f"GitHub release entry {index} has invalid tagName: {tag_name!r}"
            raise TypeError(msg)
        published_at = _parse_published_at(release.get("publishedAt"), release_index=index)
        try:
            normalized = normalize_tag(tag_name)
            _semver_sort_key(normalized)
        except ValueError:
            continue
        stable_releases[normalized] = PublishedRelease(tag=normalized, published_at=published_at)

    return list(stable_releases.values())


def _github_release_list(repo_root: Path) -> object:
    command = [
        "release",
        "list",
        "--json",
        "tagName,isDraft,isPrerelease,publishedAt",
        "--limit",
        "100",
    ]
    result = _run_tool_output(
        "gh",
        command,
        cwd=repo_root,
    )
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        msg = "could not parse GitHub release list JSON"
        raise RuntimeError(msg) from exc


def _published_stable_releases(repo_root: Path) -> list[PublishedRelease]:
    return _stable_published_releases(_github_release_list(repo_root))


def published_stable_release_tags(repo_root: Path) -> list[str]:
    """Return stable published GitHub release tags for release tooling."""
    return [release.tag for release in _published_stable_releases(repo_root)]


def _latest_published_release(repo_root: Path) -> PublishedRelease:
    stable_releases = _published_stable_releases(repo_root)
    if not stable_releases:
        msg = "expected at least one published stable semver release"
        raise RuntimeError(msg)
    return max(stable_releases, key=lambda release: release.published_at)


def _previous_release_from_list(stable_releases: list[PublishedRelease], current_tag: str) -> PublishedRelease:
    current_key = _semver_sort_key(current_tag)
    previous_releases = sorted(
        (release for release in stable_releases if _semver_sort_key(release.tag) < current_key),
        key=lambda release: _semver_sort_key(release.tag),
    )
    if not previous_releases:
        msg = f"could not find a previous stable semver release before {current_tag}"
        raise RuntimeError(msg)
    return previous_releases[-1]


def _normalize_worktree_ref_for_tag(worktree_ref: str, current_tag: str) -> str:
    try:
        normalized_ref = normalize_tag(worktree_ref)
    except ValueError:
        return worktree_ref
    return current_tag if normalized_ref == current_tag else worktree_ref


def _current_package_tag(repo_root: Path) -> str:
    cargo_toml = repo_root / "Cargo.toml"
    data = tomllib.loads(_read_text(cargo_toml))
    package = data.get("package")
    if not isinstance(package, dict):
        msg = f"could not find [package] in {cargo_toml}"
        raise TypeError(msg)
    version = package.get("version")
    if not isinstance(version, str):
        msg = f"could not find package.version in {cargo_toml}"
        raise TypeError(msg)
    return normalize_tag(version)


def _published_release_pair(repo_root: Path) -> ReportId:
    stable_releases = _published_stable_releases(repo_root)
    if len(stable_releases) < 2:
        msg = "expected at least two published stable semver releases"
        raise RuntimeError(msg)

    current = max(stable_releases, key=lambda release: release.published_at)
    previous = _previous_release_from_list(stable_releases, current.tag)
    return ReportId(current_tag=current.tag, baseline_tag=previous.tag)


def _unpublished_package_release_pair(repo_root: Path) -> ReportId:
    """Pair an unpublished package version with its newest prior release."""
    current = _current_package_tag(repo_root)
    stable_releases = _published_stable_releases(repo_root)
    if not stable_releases:
        msg = "expected at least one published stable semver release"
        raise RuntimeError(msg)
    latest_version = max(stable_releases, key=lambda release: _semver_sort_key(release.tag))
    if _semver_sort_key(current) <= _semver_sort_key(latest_version.tag):
        msg = (
            f"current package tag {current} must be newer than latest published stable release "
            f"{latest_version.tag}; update the package version before running --infer-release, "
            "or pass explicit tags for a repair"
        )
        raise ValueError(msg)
    previous = _previous_release_from_list(stable_releases, current)
    return ReportId(current_tag=current, baseline_tag=previous.tag)


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _normalize_how_to_update(text: str) -> str:
    if _HOW_TO_UPDATE_RE.search(text):
        return _HOW_TO_UPDATE_RE.sub(HOW_TO_UPDATE_SECTION, text)
    return f"{text.rstrip()}\n\n{HOW_TO_UPDATE_SECTION}"


def _replace_file(src: Path, dst: Path) -> None:
    src.replace(dst)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
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


def _restore_file(path: Path, payload: bytes | None) -> None:
    """Restore one file snapshot without relying on the publication writer."""
    if payload is None:
        path.unlink(missing_ok=True)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    restore_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".restore",
            delete=False,
        ) as tmp:
            restore_path = Path(tmp.name)
            tmp.write(payload)
            tmp.flush()
            os.fsync(tmp.fileno())
        restore_path.replace(path)
    finally:
        if restore_path is not None:
            restore_path.unlink(missing_ok=True)


def _snapshot_regular_file(path: Path, *, label: str) -> bytes | None:
    """Snapshot a mutation target while rejecting symlinks and non-files."""
    if path.is_symlink():
        msg = f"{label} must not be a symlink: {path}"
        raise ValueError(msg)
    if not path.exists():
        return None
    if not path.is_file():
        msg = f"{label} must be a regular file: {path}"
        raise ValueError(msg)
    return path.read_bytes()


def _decode_text_snapshot(payload: bytes) -> str:
    """Decode UTF-8 bytes with the universal-newline behavior of ``Path.read_text``."""
    return payload.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")


def _restore_snapshots(snapshots: tuple[tuple[Path, bytes | None], ...]) -> tuple[BaseException, ...]:
    """Attempt every file restoration and return all rollback failures."""
    errors: list[BaseException] = []
    for path, payload in reversed(snapshots):
        try:
            _restore_file(path, payload)
        except BaseException as exc:  # noqa: BLE001
            errors.append(exc)
    return tuple(errors)


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


def _format_command_timeout(command: list[str], exc: subprocess.TimeoutExpired) -> str:
    parts = [f"command timed out after {exc.timeout} seconds: {' '.join(command)}"]
    if exc.stdout:
        parts.append(f"stdout:\n{str(exc.stdout).strip()}")
    if exc.stderr:
        parts.append(f"stderr:\n{str(exc.stderr).strip()}")
    return "\n".join(parts)


def _format_command_start_failure(command: list[str], exc: BaseException) -> str:
    return f"command could not start: {' '.join(command)}: {exc}"


def _run_git(args: list[str], *, cwd: Path, timeout: int = _COMMAND_TIMEOUT_SECONDS) -> None:
    _run_git_output(args, cwd=cwd, timeout=timeout)


@contextmanager
def _temporary_detached_worktree(
    *,
    repo_root: Path,
    worktree: Path,
    revision: str,
    label: str,
) -> Iterator[Path]:
    """Create and remove a detached worktree without masking primary failures."""
    _run_git(["worktree", "add", "--detach", str(worktree), revision], cwd=repo_root)
    primary_error: BaseException | None = None
    try:
        yield worktree
    except BaseException as exc:
        primary_error = exc
        raise
    finally:
        try:
            _run_git(["worktree", "remove", "--force", str(worktree)], cwd=repo_root)
        except RuntimeError as cleanup_error:
            if primary_error is None:
                msg = f"failed to remove {label}: {cleanup_error}"
                raise RuntimeError(msg) from cleanup_error
            msg = f"operation failed ({primary_error}); additionally failed to remove {label}: {cleanup_error}"
            raise RuntimeError(msg) from primary_error


def _run_git_output(
    args: list[str],
    *,
    cwd: Path,
    timeout: int = _COMMAND_TIMEOUT_SECONDS,
    env: dict[str, str] | None = None,
) -> str:
    try:
        return run_git_command(args, cwd=cwd, timeout=timeout, env=env).stdout
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_command_failure(["git", *args], exc)) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(_format_command_timeout(["git", *args], exc)) from exc
    except (ExecutableNotFoundError, OSError) as exc:
        raise RuntimeError(_format_command_start_failure(["git", *args], exc)) from exc


def _fetch_release_tags(*, repo_root: Path, tags: list[str]) -> None:
    refspecs = [f"refs/tags/{tag}:refs/tags/{tag}" for tag in tags]
    _run_git(["fetch", "origin", *refspecs], cwd=repo_root)


def _run_tool_output(
    command: str,
    args: list[str],
    *,
    cwd: Path,
    options: ToolRunOptions | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a support command and normalize all expected launch failures."""
    resolved_options = options or ToolRunOptions()
    try:
        return run_safe_command(
            command,
            args,
            cwd=cwd,
            timeout=resolved_options.timeout,
            env=resolved_options.env,
            capture_output=not resolved_options.stream_output,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_command_failure([command, *args], exc)) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(_format_command_timeout([command, *args], exc)) from exc
    except (ExecutableNotFoundError, OSError) as exc:
        raise RuntimeError(_format_command_start_failure([command, *args], exc)) from exc


def _run_tool(
    command: str,
    args: list[str],
    *,
    cwd: Path,
    options: ToolRunOptions | None = None,
) -> None:
    _run_tool_output(
        command,
        args,
        cwd=cwd,
        options=options,
    )


def _progress(message: str) -> None:
    """Write one immediately visible workflow progress message."""
    print(f"[performance] {message}", file=sys.stderr, flush=True)


def _run_benchmark_input_gate(checkout: Path, *, env: dict[str, str] | None = None) -> None:
    """Run the shared deterministic benchmark-fixture correctness gate."""
    _progress(f"validating benchmark inputs in {checkout.name}")
    _run_tool(
        _BENCHMARK_INPUT_GATE[0],
        list(_BENCHMARK_INPUT_GATE[1:]),
        cwd=checkout,
        options=ToolRunOptions(env=env, stream_output=True),
    )


def _sha256_file(path: Path) -> str:
    """Return a lowercase SHA-256 digest for a required file."""
    if not path.is_file():
        msg = f"required provenance file is missing: {path}"
        raise FileNotFoundError(msg)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _checkout_commit(checkout: Path) -> str:
    """Return the full commit for a checkout, or an explicit unavailable label."""
    commit = _run_git_output(["--no-pager", "rev-parse", "HEAD"], cwd=checkout).strip()
    return commit or "unavailable"


def _git_clean(checkout: Path) -> bool:
    """Return whether Git reports a completely clean checkout."""
    status = _run_git_output(
        ["--no-pager", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=checkout,
    )
    return not status.strip()


def _source_state_digest(checkout: Path) -> str:
    """Hash measured library source content independent of commit cleanliness."""
    source_dir = checkout / "src"
    if not source_dir.is_dir():
        msg = f"measured library source directory is missing: {source_dir}"
        raise FileNotFoundError(msg)
    files = sorted(
        (path for path in source_dir.rglob("*") if path.is_file()),
        key=lambda path: path.relative_to(checkout).as_posix(),
    )
    if not files:
        msg = f"measured library source directory contains no files: {source_dir}"
        raise FileNotFoundError(msg)
    digest = hashlib.sha256()
    for path in files:
        relative = path.relative_to(checkout).as_posix().encode()
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def _rustc_version(checkout: Path) -> str:
    """Return one-line rustc version provenance for the active benchmark toolchain."""
    result = _run_tool_output(
        "rustc",
        ["--version"],
        cwd=checkout,
        options=ToolRunOptions(env=_benchmark_env(checkout)),
    )
    version = result.stdout.strip()
    return version or "unavailable"


def _environment_metadata(checkout: Path, *, harness_sha256: str) -> dict[str, object]:
    """Capture deterministic machine, toolchain, revision, and lock provenance."""
    cpu = cpu_description()
    os_description = " ".join(part for part in (platform.system(), platform.release(), platform.machine()) if part).strip()
    return {
        "benchmark_contract_sha256": benchmark_contract_digest(checkout),
        "cargo_lock_sha256": _sha256_file(checkout / "Cargo.lock"),
        "commit": _checkout_commit(checkout),
        "correctness_gate": "passed",
        "cpu": cpu,
        "git_clean": _git_clean(checkout),
        "harness_sha256": harness_sha256,
        "os": os_description or "unavailable",
        "rustc": _rustc_version(checkout),
        "source_state_sha256": _source_state_digest(checkout),
    }


def _require_recorded_measurement_cpu() -> str:
    """Return the CPU model required for reproducible release measurements."""
    cpu = cpu_description()
    if cpu.casefold() == "unavailable":
        msg = (
            "cannot generate release performance measurements because the CPU model is unavailable; "
            "raw local benchmarks may still be run with the ordinary Criterion recipes"
        )
        raise RuntimeError(msg)
    return cpu


def _criterion_dependency_version(checkout: Path) -> str:
    """Return the resolved Criterion version, falling back to its manifest requirement."""
    lock_data = tomllib.loads(_read_text(checkout / "Cargo.lock"))
    packages = lock_data.get("package")
    if isinstance(packages, list):
        for package in packages:
            if isinstance(package, dict) and package.get("name") == "criterion":
                version = package.get("version")
                if isinstance(version, str) and version:
                    return version

    manifest = tomllib.loads(_read_text(checkout / "Cargo.toml"))
    for section in ("dev-dependencies", "dependencies", "build-dependencies"):
        dependencies = manifest.get(section)
        if not isinstance(dependencies, dict):
            continue
        criterion = dependencies.get("criterion")
        if isinstance(criterion, str) and criterion:
            return f"manifest requirement {criterion}"
        if isinstance(criterion, dict):
            version = criterion.get("version")
            if isinstance(version, str) and version:
                return f"manifest requirement {version}"
    msg = f"Criterion dependency version is unavailable in {checkout / 'Cargo.lock'} and {checkout / 'Cargo.toml'}"
    raise ValueError(msg)


def _criterion_metadata(
    *,
    worktree: Path,
    config: GenerationConfig,
    baseline_command: tuple[str, ...],
    current_command: tuple[str, ...],
) -> dict[str, object]:
    """Record the exact Criterion selection and timing commands."""
    return {
        "baseline_command": list(baseline_command),
        "current_command": list(current_command),
        "criterion_version": _criterion_dependency_version(worktree),
        "sample": "new",
        "scope": config.scope,
        "statistic": "median",
        "suite": config.suite,
    }


def _write_local_run_provenance(
    *,
    worktree: Path,
    config: GenerationConfig,
    baseline_run: BaselineRun,
    current_command: tuple[str, ...],
) -> None:
    """Tie locally generated samples to their shared harness and environment."""
    publication = _environment_metadata(worktree, harness_sha256=baseline_run.harness_sha256)
    measurement = {
        "baseline_api_compatibility": baseline_run.api_compatibility or "none",
        "baseline_commit": baseline_run.commit,
        "benchmark_contract_sha256": publication["benchmark_contract_sha256"],
        "cargo_lock_sha256": publication["cargo_lock_sha256"],
        "cpu": publication["cpu"],
        "current_commit": publication["commit"],
        "current_git_clean": publication["git_clean"],
        "current_source_state_sha256": publication["source_state_sha256"],
        "harness_sha256": baseline_run.harness_sha256,
        "os": publication["os"],
        "rustc": publication["rustc"],
        "baseline_git_clean": baseline_run.git_clean,
        "baseline_source_state_sha256": baseline_run.source_state_sha256,
        "status": "recorded",
    }
    metadata = {
        "baseline": config.baseline_tag,
        "current": config.current_tag,
        "criterion": _criterion_metadata(
            worktree=worktree,
            config=config,
            baseline_command=baseline_run.command,
            current_command=current_command,
        ),
        "measurement": measurement,
        "mode": "shared-current-harness",
        "publication": publication,
        "schema": 2,
        "validation": {
            "baseline_api_compatibility": baseline_run.api_compatibility or "none",
            "baseline_commit": baseline_run.commit,
            "baseline_git_clean": baseline_run.git_clean,
            "baseline_revision": "passed",
            "baseline_source_state_sha256": baseline_run.source_state_sha256,
            "command": list(_BENCHMARK_INPUT_GATE),
            "current_commit": publication["commit"],
            "current_git_clean": publication["git_clean"],
            "current_revision": "passed",
            "current_source_state_sha256": publication["source_state_sha256"],
            "harness": "shared-current",
            "shared_harness_rational_inputs": baseline_run.shared_harness_rational_inputs,
        },
    }
    _write_text(
        worktree / "target" / "criterion" / _BENCHMARK_HARNESS_METADATA,
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
    )


def _write_historical_asset_provenance(
    *,
    worktree: Path,
    config: GenerationConfig,
    baseline_run: BaselineRun,
) -> None:
    """Record validation while explicitly declining to invent historical timing metadata."""
    publication = _environment_metadata(worktree, harness_sha256=baseline_run.harness_sha256)
    metadata = {
        "baseline": config.baseline_tag,
        "current": config.current_tag,
        "criterion": _criterion_metadata(
            worktree=worktree,
            config=config,
            baseline_command=("historical-release-asset", config.baseline_tag),
            current_command=("historical-release-asset", config.current_tag),
        ),
        "measurement": {
            "reason": "the downloaded release assets do not contain schema-2 measurement-environment provenance",
            "status": "unavailable",
        },
        "mode": "historical-assets",
        "publication": publication,
        "schema": 2,
        "validation": {
            "baseline_api_compatibility": baseline_run.api_compatibility or "none",
            "baseline_commit": baseline_run.commit,
            "baseline_git_clean": baseline_run.git_clean,
            "baseline_revision": "passed",
            "baseline_source_state_sha256": baseline_run.source_state_sha256,
            "command": list(_BENCHMARK_INPUT_GATE),
            "current_commit": publication["commit"],
            "current_git_clean": publication["git_clean"],
            "current_revision": "passed",
            "current_source_state_sha256": publication["source_state_sha256"],
            "harness": "shared-current",
            "shared_harness_rational_inputs": baseline_run.shared_harness_rational_inputs,
        },
    }
    _write_text(
        worktree / "target" / "criterion" / _BENCHMARK_HARNESS_METADATA,
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
    )


def _current_rust_toolchain(checkout: Path) -> str | None:
    rust_toolchain = checkout / "rust-toolchain.toml"
    if not rust_toolchain.exists():
        return None
    data = tomllib.loads(_read_text(rust_toolchain))
    toolchain = data.get("toolchain")
    if not isinstance(toolchain, dict):
        return None
    channel = toolchain.get("channel")
    return channel if isinstance(channel, str) else None


def _benchmark_env(checkout: Path) -> dict[str, str] | None:
    if "RUSTUP_TOOLCHAIN" in os.environ:
        return None
    toolchain = _current_rust_toolchain(checkout)
    if toolchain is None:
        return None
    env = os.environ.copy()
    env["RUSTUP_TOOLCHAIN"] = toolchain
    return env


def _append_rustflag(env: dict[str, str], flag: str) -> None:
    """Append one rustc flag without discarding caller-selected codegen flags."""
    encoded = env.get("CARGO_ENCODED_RUSTFLAGS")
    if encoded is not None:
        env["CARGO_ENCODED_RUSTFLAGS"] = "\x1f".join(part for part in (encoded, flag) if part)
        return

    rustflags = env.get("RUSTFLAGS", "").strip()
    env["RUSTFLAGS"] = f"{rustflags} {flag}".strip()


def _shared_harness_rational_inputs(checkout: Path) -> bool:
    """Return whether the installed exact benchmark harness emits rational-input groups."""
    exact_bench = checkout / "benches" / "exact.rs"
    if not exact_bench.is_file():
        return False
    source = _read_text(exact_bench)
    return "fn bench_rational_input<" in source and "rational_input_d{D}" in source


def _baseline_api_compatibility(
    *,
    current_tag: str,
    baseline_tag: str,
    shared_harness_rational_inputs: bool,
) -> str | None:
    """Return the baseline adapter selected by the shared compatibility resolver."""
    resolved = resolve_shared_harness_compatibility(
        current=normalize_tag(current_tag),
        baseline=normalize_tag(baseline_tag),
        shared_harness_rational_inputs=shared_harness_rational_inputs,
    )
    if resolved.baseline_api_compatibility == NO_API_COMPATIBILITY:
        return None
    return resolved.baseline_api_compatibility


def _comparison_benchmark_env(checkout: Path, *, api_compatibility: str | None = None) -> dict[str, str]:
    """Build a comparable benchmark environment for current or historical code."""
    env = _benchmark_env(checkout)
    if env is None:
        env = os.environ.copy()

    # The shared current manifest can enable lints unknown to historical source.
    # Cap diagnostics for both revisions so lint-policy drift cannot prevent a
    # performance comparison; the cap changes diagnostics, not code generation.
    _append_rustflag(env, _COMPARISON_LINT_CAP)

    if api_compatibility is not None:
        _append_rustflag(env, f"--cfg={api_compatibility}")
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


def _copy_criterion_sample(*, criterion_dir: Path, source_sample: str, target_sample: str) -> None:
    copied = 0
    for source in list(criterion_dir.rglob(source_sample)):
        if not source.is_dir() or not (source / "estimates.json").exists():
            continue
        target = source.parent / target_sample
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
        copied += 1
    if copied == 0:
        msg = f"could not find Criterion sample {source_sample!r} under {criterion_dir}"
        raise FileNotFoundError(msg)


def _selected_criterion_groups(criterion_dir: Path, *, suite: str) -> list[Path]:
    """Return selected Criterion group directories in deterministic order."""
    groups: list[Path] = []
    if not criterion_dir.is_dir():
        return groups
    for child in criterion_dir.iterdir():
        if not child.is_dir():
            continue
        is_exact = child.name.startswith("exact_") or re.fullmatch(r"rational_input_d[0-9]+", child.name) is not None
        is_vs_linalg = re.fullmatch(r"d[0-9]+", child.name) is not None
        if (suite in {"all", "exact"} and is_exact) or (suite in {"all", "vs_linalg"} and is_vs_linalg):
            groups.append(child)
    return sorted(groups, key=lambda path: path.name)


def _purge_criterion_new_samples(*, criterion_dir: Path, suite: str) -> list[Path]:
    """Remove stale selected-suite `new` samples while preserving named baselines."""
    removed: list[Path] = []
    for group in _selected_criterion_groups(criterion_dir, suite=suite):
        for sample in sorted(group.glob("*/new")):
            if sample.is_dir():
                shutil.rmtree(sample)
                removed.append(sample)
    return removed


def _benchmark_harness_files(checkout: Path) -> list[Path]:
    """Return every file that defines the comparable benchmark harness."""
    files: list[Path] = []
    for relative in _BENCHMARK_HARNESS_FILES:
        path = checkout / relative
        if not path.is_file():
            msg = f"benchmark harness file is missing: {path}"
            raise FileNotFoundError(msg)
        files.append(path)
    for relative in _BENCHMARK_HARNESS_DIRS:
        directory = checkout / relative
        if not directory.is_dir():
            msg = f"benchmark harness directory is missing: {directory}"
            raise FileNotFoundError(msg)
        files.extend(path for path in directory.rglob("*") if path.is_file())
    return sorted(files, key=lambda path: path.relative_to(checkout).as_posix())


def _benchmark_harness_digest(checkout: Path) -> str:
    """Hash benchmark sources, recipes, dependency resolution, and toolchain."""
    digest = hashlib.sha256()
    for path in _benchmark_harness_files(checkout):
        relative = path.relative_to(checkout).as_posix().encode()
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def _install_shared_benchmark_harness(*, source: Path, destination: Path) -> str:
    """Replace a baseline checkout's harness with the current harness."""
    source_files = _benchmark_harness_files(source)
    for relative in _BENCHMARK_HARNESS_DIRS:
        source_dir = source / relative
        destination_dir = destination / relative
        if destination_dir.exists():
            shutil.rmtree(destination_dir)
        shutil.copytree(source_dir, destination_dir)
    for relative in _BENCHMARK_HARNESS_FILES:
        destination_file = destination / relative
        destination_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / relative, destination_file)

    expected = _benchmark_harness_digest(source)
    actual = _benchmark_harness_digest(destination)
    if actual != expected:
        msg = f"shared benchmark harness copy changed content: expected {expected}, found {actual}"
        raise RuntimeError(msg)
    if not source_files:
        msg = "benchmark harness unexpectedly contained no files"
        raise RuntimeError(msg)
    return expected


def _has_suite_aware_baseline_recipe(worktree: Path) -> bool:
    justfile = worktree / "justfile"
    return justfile.exists() and re.search(r'(?m)^bench-save-baseline\s+tag\s+suite(?:=|"|:|\s)', _read_text(justfile)) is not None


def _baseline_tool_args(*, baseline_tag: str, suite: str, baseline_worktree: Path) -> tuple[str, list[str]]:
    if suite == _DEFAULT_SUITE:
        return ("just", ["bench-save-baseline", baseline_tag])
    if _has_suite_aware_baseline_recipe(baseline_worktree):
        return ("just", ["bench-save-baseline", baseline_tag, suite])
    match suite:
        case "exact":
            return ("cargo", ["bench", "--locked", "--features", "bench,exact", "--bench", "exact", "--", "--save-baseline", baseline_tag])
        case "vs_linalg":
            return ("cargo", ["bench", "--locked", "--features", "bench", "--bench", "vs_linalg", "--", "--save-baseline", baseline_tag])
        case _:
            msg = f"unsupported benchmark suite: {suite}"
            raise ValueError(msg)


def _latest_recipe_args(*, suite: str) -> list[str]:
    match suite:
        case "all":
            return ["bench-latest"]
        case "exact":
            return ["bench-exact"]
        case "vs_linalg":
            return ["bench-vs-linalg-la-stack"]
        case _:
            msg = f"unsupported benchmark suite: {suite}"
            raise ValueError(msg)


def _fallback_current_command(*, suite: str) -> tuple[str, ...]:
    """Return the Cargo command used when current benchmark recipes are unavailable."""
    match suite:
        case "all":
            return ("cargo", "bench", "--locked", "--features", "bench,exact")
        case "exact":
            return ("cargo", "bench", "--locked", "--features", "bench,exact", "--bench", "exact")
        case "vs_linalg":
            return ("cargo", "bench", "--locked", "--features", "bench", "--bench", "vs_linalg")
        case _:
            msg = f"unsupported benchmark suite: {suite}"
            raise ValueError(msg)


def _generate_release_baseline(  # noqa: PLR0913
    *,
    current_tag: str,
    baseline_tag: str,
    suite: str,
    repo_root: Path,
    target_worktree: Path,
    tmp_dir: Path,
) -> BaselineRun:
    baseline_worktree = tmp_dir / "baseline-worktree"
    with _temporary_detached_worktree(
        repo_root=repo_root,
        worktree=baseline_worktree,
        revision=baseline_tag,
        label="baseline worktree",
    ):
        harness_sha256 = _install_shared_benchmark_harness(
            source=target_worktree,
            destination=baseline_worktree,
        )
        baseline_command, baseline_args = _baseline_tool_args(
            baseline_tag=baseline_tag,
            suite=suite,
            baseline_worktree=baseline_worktree,
        )
        shared_harness_rational_inputs = _shared_harness_rational_inputs(target_worktree)
        api_compatibility = _baseline_api_compatibility(
            current_tag=current_tag,
            baseline_tag=baseline_tag,
            shared_harness_rational_inputs=shared_harness_rational_inputs,
        )
        benchmark_env = _comparison_benchmark_env(repo_root, api_compatibility=api_compatibility)
        _run_benchmark_input_gate(baseline_worktree, env=benchmark_env)
        _progress(f"running {suite} baseline benchmarks for {baseline_tag}")
        _run_tool(
            baseline_command,
            baseline_args,
            cwd=baseline_worktree,
            options=ToolRunOptions(
                timeout=_BENCH_TIMEOUT_SECONDS,
                env=benchmark_env,
                stream_output=True,
            ),
        )
        _progress(f"completed {suite} baseline benchmarks for {baseline_tag}")
        baseline_criterion = baseline_worktree / "target" / "criterion"
        if not baseline_criterion.is_dir():
            msg = f"generated baseline Criterion results were not found: {baseline_criterion}"
            raise FileNotFoundError(msg)
        target_criterion = target_worktree / "target" / "criterion"
        target_criterion.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(baseline_criterion, target_criterion, dirs_exist_ok=True)
        return BaselineRun(
            commit=_checkout_commit(baseline_worktree),
            command=(baseline_command, *baseline_args),
            harness_sha256=harness_sha256,
            git_clean=_git_clean(baseline_worktree),
            source_state_sha256=_source_state_digest(baseline_worktree),
            api_compatibility=api_compatibility,
            shared_harness_rational_inputs=shared_harness_rational_inputs,
        )


def _prepare_local_release_baseline(  # noqa: PLR0913
    *,
    current_tag: str,
    baseline_tag: str,
    suite: str,
    repo_root: Path,
    target_worktree: Path,
    tmp_dir: Path,
) -> BaselineRun:
    return _generate_release_baseline(
        current_tag=current_tag,
        baseline_tag=baseline_tag,
        suite=suite,
        repo_root=repo_root,
        target_worktree=target_worktree,
        tmp_dir=tmp_dir,
    )


def _validate_release_revision(
    *,
    current_tag: str,
    revision: str,
    repo_root: Path,
    harness_source: Path,
    tmp_dir: Path,
) -> BaselineRun:
    """Validate one historical revision with the shared current fixture harness."""
    validation_worktree = tmp_dir / "baseline-validation-worktree"
    with _temporary_detached_worktree(
        repo_root=repo_root,
        worktree=validation_worktree,
        revision=revision,
        label="baseline validation worktree",
    ):
        harness_sha256 = _install_shared_benchmark_harness(
            source=harness_source,
            destination=validation_worktree,
        )
        shared_harness_rational_inputs = _shared_harness_rational_inputs(harness_source)
        api_compatibility = _baseline_api_compatibility(
            current_tag=current_tag,
            baseline_tag=revision,
            shared_harness_rational_inputs=shared_harness_rational_inputs,
        )
        _run_benchmark_input_gate(
            validation_worktree,
            env=_comparison_benchmark_env(repo_root, api_compatibility=api_compatibility),
        )
        return BaselineRun(
            commit=_checkout_commit(validation_worktree),
            command=("historical-release-asset", revision),
            harness_sha256=harness_sha256,
            git_clean=_git_clean(validation_worktree),
            source_state_sha256=_source_state_digest(validation_worktree),
            api_compatibility=api_compatibility,
            shared_harness_rational_inputs=shared_harness_rational_inputs,
        )


def _prepare_github_release_assets(*, current_tag: str, baseline_tag: str, repo_root: Path, target_worktree: Path, tmp_dir: Path) -> None:
    baseline_archive = _download_release_baseline(
        baseline_tag=baseline_tag,
        download_dir=tmp_dir,
        repo_root=repo_root,
    )
    current_archive = _download_release_baseline(
        baseline_tag=current_tag,
        download_dir=tmp_dir,
        repo_root=repo_root,
    )
    target_dir = target_worktree / "target"
    _safe_extract_tar(baseline_archive, target_dir)
    _safe_extract_tar(current_archive, target_dir)
    # Published artifacts retain their release-specific harnesses. Never let an
    # embedded or stale local manifest claim that these samples shared one.
    metadata_path = target_dir / "criterion" / _BENCHMARK_HARNESS_METADATA
    if metadata_path.is_symlink() or metadata_path.is_file():
        metadata_path.unlink()
    elif metadata_path.exists():
        msg = f"historical benchmark harness provenance path is not a file: {metadata_path}"
        raise ValueError(msg)
    _copy_criterion_sample(criterion_dir=target_dir / "criterion", source_sample=current_tag, target_sample="new")


def _apply_current_diff_to_worktree(*, repo_root: Path, worktree: Path) -> None:
    # Diff HEAD directly so staged and unstaged tracked changes participate,
    # while unrelated untracked files stay outside the benchmark worktree.
    with tempfile.TemporaryDirectory(prefix="la-stack-current-tree-index-") as tmp:
        temporary_dir = Path(tmp)
        patch_path = temporary_dir / "current-tree.patch"
        _run_git_output(
            ["diff", "--binary", f"--output={patch_path}", "HEAD", "--", "."],
            cwd=repo_root,
        )
        diff = patch_path.read_bytes()
    if diff.strip():
        try:
            run_git_command_with_input(["apply", "--binary"], diff, cwd=worktree)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(_format_command_failure(["git", "apply", "--binary"], exc)) from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(_format_command_timeout(["git", "apply", "--binary"], exc)) from exc
        except (ExecutableNotFoundError, OSError) as exc:
            raise RuntimeError(_format_command_start_failure(["git", "apply", "--binary"], exc)) from exc


def _has_current_release_signal_tooling(worktree: Path) -> bool:
    justfile = worktree / "justfile"
    bench_compare = worktree / "scripts" / "bench_compare.py"
    if not justfile.exists() or not bench_compare.exists():
        return False

    justfile_text = _read_text(justfile)
    bench_compare_text = _read_text(bench_compare)
    return re.search(r"(?m)^bench-latest(?:[ :]|$)", justfile_text) is not None and '"--suite"' in bench_compare_text and '"--scope"' in bench_compare_text


def _render_report(
    *,
    worktree: Path,
    report: Path,
    artifacts: ArtifactPaths,
    config: GenerationConfig,
) -> None:
    """Render Markdown and export validated comparison artifacts."""
    command = [
        "run",
        "--locked",
        "--project",
        str(config.repo_root),
        "bench-compare",
        config.baseline_tag,
        "--repo-root",
        str(worktree),
        "--criterion-dir",
        str(worktree / "target" / "criterion"),
        "--suite",
        config.suite,
        "--scope",
        config.scope,
    ]
    command.extend(
        [
            "--csv-output",
            str(artifacts.csv),
            "--provenance-output",
            str(artifacts.provenance),
        ]
    )
    command.extend(["--output", str(report)])
    _run_tool(
        "uv",
        command,
        cwd=config.repo_root,
    )


def _run_benchmarks_and_render_report(
    *,
    worktree: Path,
    report: Path,
    artifacts: ArtifactPaths,
    config: GenerationConfig,
    baseline_run: BaselineRun,
) -> None:
    benchmark_env = _comparison_benchmark_env(config.repo_root)
    _run_benchmark_input_gate(worktree, env=benchmark_env)
    _purge_criterion_new_samples(
        criterion_dir=worktree / "target" / "criterion",
        suite=config.suite,
    )
    if _has_current_release_signal_tooling(worktree):
        current_command = ("just", *_latest_recipe_args(suite=config.suite))
    else:
        current_command = _fallback_current_command(suite=config.suite)
    _progress(f"running current {config.suite} benchmarks")
    _run_tool(
        current_command[0],
        list(current_command[1:]),
        cwd=worktree,
        options=ToolRunOptions(
            timeout=_BENCH_TIMEOUT_SECONDS,
            env=benchmark_env,
            stream_output=True,
        ),
    )
    _progress(f"completed current {config.suite} benchmarks")
    _write_local_run_provenance(
        worktree=worktree,
        config=config,
        baseline_run=baseline_run,
        current_command=current_command,
    )
    _render_report(worktree=worktree, report=report, artifacts=artifacts, config=config)


def _default_artifact_paths(repo_root: Path) -> ArtifactPaths:
    """Return the canonical retained comparison-artifact pair."""
    return ArtifactPaths(
        csv=repo_root / _DEFAULT_ARTIFACT_CSV,
        provenance=repo_root / _DEFAULT_ARTIFACT_PROVENANCE,
    )


@contextmanager
def _generated_report_in_temp_worktree(
    *,
    config: GenerationConfig,
    published_artifacts: ArtifactPaths,
) -> Iterator[GeneratedReport]:
    """Generate and expose a report before its temporary worktree is removed."""
    if config.baseline_source == "local":
        _require_recorded_measurement_cpu()
    with tempfile.TemporaryDirectory(prefix="la-stack-performance-") as tmp:
        tmp_dir = Path(tmp)
        worktree = tmp_dir / "worktree"
        report = tmp_dir / f"{config.current_tag}-vs-{config.baseline_tag}.md"
        temporary_artifacts = ArtifactPaths(
            csv=tmp_dir / "performance.csv",
            provenance=tmp_dir / "performance.provenance.json",
        )

        with _temporary_detached_worktree(
            repo_root=config.repo_root,
            worktree=worktree,
            revision=config.worktree_ref,
            label="temporary worktree",
        ):
            if config.apply_current_diff:
                _progress("applying staged and unstaged tracked changes; untracked files are excluded")
                _apply_current_diff_to_worktree(repo_root=config.repo_root, worktree=worktree)
            if config.baseline_source == "github-assets":
                _prepare_github_release_assets(
                    current_tag=config.current_tag,
                    baseline_tag=config.baseline_tag,
                    repo_root=config.repo_root,
                    target_worktree=worktree,
                    tmp_dir=tmp_dir,
                )
                baseline_run = _validate_release_revision(
                    current_tag=config.current_tag,
                    revision=config.baseline_tag,
                    repo_root=config.repo_root,
                    harness_source=worktree,
                    tmp_dir=tmp_dir,
                )
                _run_benchmark_input_gate(
                    worktree,
                    env=_comparison_benchmark_env(config.repo_root),
                )
                _write_historical_asset_provenance(
                    worktree=worktree,
                    config=config,
                    baseline_run=baseline_run,
                )
                _render_report(
                    worktree=worktree,
                    report=report,
                    artifacts=temporary_artifacts,
                    config=config,
                )
            else:
                baseline_run = _prepare_local_release_baseline(
                    current_tag=config.current_tag,
                    baseline_tag=config.baseline_tag,
                    suite=config.suite,
                    repo_root=config.repo_root,
                    target_worktree=worktree,
                    tmp_dir=tmp_dir,
                )
                _run_benchmarks_and_render_report(
                    worktree=worktree,
                    report=report,
                    artifacts=temporary_artifacts,
                    config=config,
                    baseline_run=baseline_run,
                )
            report_text = _read_text(report)
            bundle = load_bundle(temporary_artifacts)
            with publish_bundle(published_artifacts, bundle):
                durable_text = render_release_artifacts(published_artifacts)
                if durable_text != report_text:
                    msg = "durable performance-comparison artifacts did not reproduce the generated Markdown"
                    raise ValueError(msg)
                yield GeneratedReport(text=durable_text)


def _current_archive_state(
    *,
    current: Path,
    source_id: ReportId,
    archive_dir: Path,
) -> tuple[str | None, ReportId | None, Path | None]:
    """Return normalized current-report state and its prospective archive path."""
    if not current.exists():
        return None, None, None
    current_text = _normalize_how_to_update(_read_text(current))
    current_id = parse_report_id(current_text)
    archive_path = archive_dir / current_id.archive_name if current_id != source_id else None
    return current_text, current_id, archive_path


def _validate_promotion_paths(request: PromotionRequest, *, index_path: Path, archive_path: Path | None) -> None:
    """Reject aliases and archive-contained scratch output before mutation."""
    named_paths: dict[str, Path] = {
        "current report": request.current,
        "archive index": index_path,
    }
    if request.source_path is not None:
        named_paths["source report"] = request.source_path
    if request.output is not None:
        named_paths["rendered output"] = request.output
        if request.output.resolve(strict=False).is_relative_to(request.archive_dir.resolve(strict=False)):
            msg = f"rendered output must not be written inside the archive directory: {request.output}"
            raise ValueError(msg)
    if archive_path is not None:
        named_paths["archived report"] = archive_path
    if request.reserved_paths is not None:
        named_paths.update(request.reserved_paths)
    ensure_distinct_paths(named_paths)


def _existing_archive_matches(*, archive_path: Path | None, current_id: ReportId | None, current_text: str | None) -> bool:
    """Validate a pre-existing canonical archive before trusting it."""
    if archive_path is None or (not archive_path.exists() and not archive_path.is_symlink()):
        return False
    if current_id is None or current_text is None:
        msg = "existing archive validation invariant violated"
        raise AssertionError(msg)
    archived_payload = _snapshot_regular_file(archive_path, label="existing archived report")
    if archived_payload is None:
        msg = "existing archive snapshot invariant violated"
        raise AssertionError(msg)
    archived_text = _normalize_how_to_update(_decode_text_snapshot(archived_payload))
    try:
        archived_id = parse_report_id(archived_text)
    except (TypeError, ValueError) as exc:
        msg = f"existing archived report does not match the current report being replaced: {archive_path}"
        raise ValueError(msg) from exc
    if archived_id != current_id or archived_text != current_text:
        msg = f"existing archived report does not match the current report being replaced: {archive_path}"
        raise ValueError(msg)
    return True


def _promotion_snapshots(
    request: PromotionRequest,
    *,
    index_path: Path,
    archive_path: Path | None,
    archive_exists: bool,
) -> tuple[tuple[Path, bytes | None], ...]:
    """Snapshot every path mutated by a report promotion."""
    snapshots: list[tuple[Path, bytes | None]] = []
    if archive_path is not None and not archive_exists:
        snapshots.append((archive_path, None))
    snapshots.extend(
        (
            (request.current, _snapshot_regular_file(request.current, label="current report")),
            (index_path, _snapshot_regular_file(index_path, label="archive index")),
        )
    )
    if request.output is not None:
        snapshots.append((request.output, _snapshot_regular_file(request.output, label="rendered output")))
    return tuple(snapshots)


def _promote_report_text(*, source_text: str, request: PromotionRequest) -> ReportId:
    """Archive the old report and atomically promote validated Markdown text."""
    if request.expected.current_tag == request.expected.baseline_tag:
        msg = "cannot promote a same-version local performance comparison as release documentation"
        raise ValueError(msg)
    source_text = _normalize_how_to_update(source_text)
    source_id = parse_report_id(source_text)
    if source_id != request.expected:
        msg = (
            "benchmark report does not match requested release pair: "
            f"found {source_id.current_tag} vs {source_id.baseline_tag}, "
            f"expected {request.expected.current_tag} vs {request.expected.baseline_tag}"
        )
        raise ValueError(msg)

    current_text, current_id, archive_path = _current_archive_state(
        current=request.current,
        source_id=source_id,
        archive_dir=request.archive_dir,
    )
    index_path = request.archive_dir / "README.md"
    _validate_promotion_paths(request, index_path=index_path, archive_path=archive_path)
    archive_exists = _existing_archive_matches(
        archive_path=archive_path,
        current_id=current_id,
        current_text=current_text,
    )
    snapshots = _promotion_snapshots(
        request,
        index_path=index_path,
        archive_path=archive_path,
        archive_exists=archive_exists,
    )
    try:
        if archive_path is not None and not archive_exists:
            if current_text is None:
                msg = "archive promotion invariant violated"
                raise AssertionError(msg)
            _write_text(archive_path, current_text)
        _write_text(request.current, source_text)
        update_archive_index(request.archive_dir)
        if request.output is not None:
            _write_text(request.output, source_text)
    except BaseException as promotion_error:
        rollback_errors = _restore_snapshots(snapshots)
        if rollback_errors:
            group_message = "performance report promotion and rollback failed"
            raise BaseExceptionGroup(group_message, [promotion_error, *rollback_errors]) from None
        raise
    return source_id


def promote_report(
    *,
    source: Path,
    current: Path,
    archive_dir: Path,
    expected_current_tag: str,
    expected_baseline_tag: str,
) -> ReportId:
    """Archive the old committed report and promote *source* as the current one."""
    return _promote_report_text(
        source_text=_read_text(source),
        request=PromotionRequest(
            current=current,
            archive_dir=archive_dir,
            expected=ReportId(
                current_tag=normalize_tag(expected_current_tag),
                baseline_tag=normalize_tag(expected_baseline_tag),
            ),
            source_path=source,
        ),
    )


def generate_and_promote_worktree_report(
    *,
    output: Path,
    current: Path,
    archive_dir: Path,
    config: GenerationConfig,
    artifacts: ArtifactPaths | None = None,
) -> ReportId:
    """Generate a comparison in a temp worktree, then promote it."""
    current_tag = normalize_tag(config.current_tag)
    baseline_tag = normalize_tag(config.baseline_tag)
    if current_tag == baseline_tag:
        msg = "release performance promotion requires distinct current and baseline releases"
        raise ValueError(msg)
    config = GenerationConfig(
        repo_root=config.repo_root,
        current_tag=current_tag,
        baseline_tag=baseline_tag,
        worktree_ref=config.worktree_ref,
        suite=config.suite,
        scope=config.scope,
        apply_current_diff=config.apply_current_diff,
        baseline_source=config.baseline_source,
    )
    published_artifacts = artifacts or _default_artifact_paths(config.repo_root)
    expected = ReportId(current_tag=current_tag, baseline_tag=baseline_tag)
    reserved_paths = {
        "artifact CSV": published_artifacts.csv,
        "artifact provenance": published_artifacts.provenance,
    }
    request = PromotionRequest(
        current=current,
        archive_dir=archive_dir,
        expected=expected,
        output=output,
        reserved_paths=reserved_paths,
    )
    _current_text, _current_id, archive_path = _current_archive_state(
        current=current,
        source_id=expected,
        archive_dir=archive_dir,
    )
    _validate_promotion_paths(
        request,
        index_path=archive_dir / "README.md",
        archive_path=archive_path,
    )
    with _generated_report_in_temp_worktree(
        config=config,
        published_artifacts=published_artifacts,
    ) as generated:
        return _promote_report_text(
            source_text=generated.text,
            request=request,
        )


def generate_worktree_report(
    *,
    output: Path,
    config: GenerationConfig,
    artifacts: ArtifactPaths | None = None,
) -> ReportId:
    """Generate a comparison in a temp worktree and write it to *output*."""
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
        baseline_source=config.baseline_source,
    )
    published_artifacts = artifacts or ArtifactPaths(
        csv=output.with_suffix(".csv"),
        provenance=output.with_suffix(".provenance.json"),
    )
    ensure_distinct_paths(
        {
            "rendered output": output,
            "artifact CSV": published_artifacts.csv,
            "artifact provenance": published_artifacts.provenance,
        }
    )
    with _generated_report_in_temp_worktree(
        config=config,
        published_artifacts=published_artifacts,
    ) as generated:
        report_text = _normalize_how_to_update(generated.text)
        report_id = parse_report_id(report_text)
        expected = ReportId(current_tag=current_tag, baseline_tag=baseline_tag)
        if report_id != expected:
            msg = (
                "benchmark report does not match requested release pair: "
                f"found {report_id.current_tag} vs {report_id.baseline_tag}, "
                f"expected {expected.current_tag} vs {expected.baseline_tag}"
            )
            raise ValueError(msg)
        _write_text(output, report_text)
        return report_id


def render_and_promote_artifacts(
    *,
    artifacts: ArtifactPaths,
    output: Path,
    current: Path,
    archive_dir: Path,
) -> ReportId:
    """Reproduce and promote Markdown using only retained CSV/JSON inputs."""
    bundle = load_bundle(artifacts)
    expected = ReportId(
        current_tag=normalize_tag(bundle.context.release.current),
        baseline_tag=normalize_tag(bundle.context.release.baseline),
    )
    if expected.current_tag == expected.baseline_tag:
        msg = "cannot promote a same-version local performance comparison as release documentation"
        raise ValueError(msg)
    report_text = _normalize_how_to_update(render_release_artifacts(artifacts))
    observed = parse_report_id(report_text)
    if observed != expected:
        msg = (
            "rendered benchmark report does not match retained release pair: "
            f"found {observed.current_tag} vs {observed.baseline_tag}, "
            f"expected {expected.current_tag} vs {expected.baseline_tag}"
        )
        raise ValueError(msg)
    return _promote_report_text(
        source_text=report_text,
        request=PromotionRequest(
            current=current,
            archive_dir=archive_dir,
            expected=expected,
            output=output,
            reserved_paths={
                "artifact CSV": artifacts.csv,
                "artifact provenance": artifacts.provenance,
            },
        ),
    )


def resolve_archive_request(options: ArchiveRequestOptions) -> ResolvedArchiveRequest:
    """Resolve explicit, package-inferred, or latest-published release arguments."""
    current_tag = options.current_tag
    baseline_tag = options.baseline_tag
    worktree_ref = options.worktree_ref
    repo_root = options.repo_root
    published_latest = options.published_latest
    infer_release = options.infer_release
    current_vs_latest = options.current_vs_latest
    requested_modes = sum((published_latest, infer_release, current_vs_latest))
    if requested_modes > 1:
        msg = "choose only one of --published-latest, --infer-release, or --current-vs-latest"
        raise ValueError(msg)

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
            tags_to_fetch=(published_pair.current_tag, published_pair.baseline_tag),
        )

    if infer_release:
        if current_tag is not None or baseline_tag is not None:
            msg = "do not pass current_tag or baseline_tag with --infer-release"
            raise ValueError(msg)
        inferred_pair = _unpublished_package_release_pair(repo_root)
        return ResolvedArchiveRequest(
            current_tag=inferred_pair.current_tag,
            baseline_tag=inferred_pair.baseline_tag,
            worktree_ref=worktree_ref,
            tags_to_fetch=(inferred_pair.baseline_tag,),
        )

    if current_vs_latest:
        if current_tag is not None or baseline_tag is not None:
            msg = "do not pass current_tag or baseline_tag with --current-vs-latest"
            raise ValueError(msg)
        inferred_current = _current_package_tag(repo_root)
        latest = _latest_published_release(repo_root).tag
        return ResolvedArchiveRequest(
            current_tag=inferred_current,
            baseline_tag=latest,
            worktree_ref=worktree_ref,
            tags_to_fetch=(latest,),
        )

    if current_tag is None or baseline_tag is None:
        msg = "current_tag and baseline_tag are required unless an inference mode is used"
        raise ValueError(msg)
    normalized_current = normalize_tag(current_tag)
    normalized_baseline = normalize_tag(baseline_tag)
    return ResolvedArchiveRequest(
        current_tag=normalized_current,
        baseline_tag=normalized_baseline,
        worktree_ref=_normalize_worktree_ref_for_tag(worktree_ref, normalized_current),
        tags_to_fetch=(normalized_baseline,),
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate local benchmark comparisons or promote validated comparison artifacts into release documentation.",
    )
    parser.add_argument("current_tag", nargs="?", help="Current package-version label, e.g. v0.4.3")
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
        "--output",
        default=_DEFAULT_SOURCE,
        help=f"Rendered local comparison report path (default: {_DEFAULT_SOURCE})",
    )
    parser.add_argument(
        "--artifact-csv",
        default=_DEFAULT_ARTIFACT_CSV,
        help=f"Retained comparison CSV path (default: {_DEFAULT_ARTIFACT_CSV})",
    )
    parser.add_argument(
        "--artifact-provenance",
        default=_DEFAULT_ARTIFACT_PROVENANCE,
        help=f"Retained comparison JSON provenance path (default: {_DEFAULT_ARTIFACT_PROVENANCE})",
    )
    parser.add_argument(
        "--archive-dir",
        default=_DEFAULT_ARCHIVE_DIR,
        help=f"Archive directory for older reports (default: {_DEFAULT_ARCHIVE_DIR})",
    )
    parser.add_argument(
        "--generate-in-temp-worktree",
        action="store_true",
        help="Generate the comparison in a temporary detached worktree.",
    )
    parser.add_argument(
        "--published-latest",
        action="store_true",
        help="Infer the latest stable published GitHub release and its previous stable release.",
    )
    parser.add_argument(
        "--infer-release",
        action="store_true",
        help="Require an unpublished Cargo.toml version and compare it with the previous stable published release.",
    )
    parser.add_argument(
        "--current-vs-latest",
        action="store_true",
        help="Infer current_tag from Cargo.toml and baseline_tag from the latest stable published release.",
    )
    parser.add_argument(
        "--github-assets",
        action="store_true",
        help="Compare stored GitHub Release benchmark assets instead of generating the baseline locally.",
    )
    parser.add_argument(
        "--output-only",
        action="store_true",
        help="Write the generated report to --output without promoting docs/PERFORMANCE.md.",
    )
    parser.add_argument(
        "--local-report",
        action="store_true",
        help=(
            "Generate a local report and retained comparison artifacts without promoting release documentation; "
            "the current package version may match the baseline release."
        ),
    )
    parser.add_argument(
        "--promote-artifacts",
        action="store_true",
        help="Render and promote from retained CSV/JSON artifacts without benchmarks or worktrees.",
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
        default=None,
        choices=_SUPPORTED_SUITES,
        help=f"Benchmark suite for --generate-in-temp-worktree (default: {_DEFAULT_SUITE})",
    )
    parser.add_argument(
        "--scope",
        default=None,
        choices=_SUPPORTED_SCOPES,
        help=f"Comparison scope for --generate-in-temp-worktree (default: {_DEFAULT_SCOPE})",
    )
    return parser


def _resolve_cli_paths(root: Path, args: argparse.Namespace) -> ArchivePaths:
    source = Path(args.source)
    current = Path(args.current)
    output = Path(args.output)
    archive_dir = Path(args.archive_dir)
    artifact_csv = Path(args.artifact_csv)
    artifact_provenance = Path(args.artifact_provenance)
    if not source.is_absolute():
        source = root / source
    if not current.is_absolute():
        current = root / current
    if not output.is_absolute():
        output = root / output
    if not archive_dir.is_absolute():
        archive_dir = root / archive_dir
    if not artifact_csv.is_absolute():
        artifact_csv = root / artifact_csv
    if not artifact_provenance.is_absolute():
        artifact_provenance = root / artifact_provenance
    return ArchivePaths(
        source=source,
        current=current,
        output=output,
        archive_dir=archive_dir,
        artifacts=ArtifactPaths(csv=artifact_csv, provenance=artifact_provenance),
    )


def _fetch_required_tags(*, request: ResolvedArchiveRequest, repo_root: Path, include_current: bool) -> None:
    tags_to_fetch = request.tags_to_fetch
    if include_current and request.current_tag not in tags_to_fetch:
        tags_to_fetch = (*tags_to_fetch, request.current_tag)
    if tags_to_fetch:
        _fetch_release_tags(repo_root=repo_root, tags=list(dict.fromkeys(tags_to_fetch)))


def _generation_config(*, args: argparse.Namespace, request: ResolvedArchiveRequest, repo_root: Path) -> GenerationConfig:
    return GenerationConfig(
        repo_root=repo_root,
        current_tag=request.current_tag,
        baseline_tag=request.baseline_tag,
        worktree_ref=request.worktree_ref,
        suite=cast("BenchmarkSuite", args.suite or _DEFAULT_SUITE),
        scope=cast("ComparisonScope", args.scope or _DEFAULT_SCOPE),
        apply_current_diff=not args.no_apply_current_diff and not args.github_assets,
        baseline_source="github-assets" if args.github_assets else "local",
    )


def _validate_generation_request(
    *,
    args: argparse.Namespace,
    request: ResolvedArchiveRequest,
    repo_root: Path,
) -> None:
    """Reject invalid release identities before fetching tags or benchmarking."""
    if not args.local_report and request.current_tag == request.baseline_tag:
        msg = "release performance artifacts require distinct current and baseline releases"
        raise ValueError(msg)

    if args.generate_in_temp_worktree and not args.github_assets and args.current_tag is not None and request.worktree_ref == "HEAD":
        package_tag = _current_package_tag(repo_root)
        if request.current_tag != package_tag:
            msg = (
                f"explicit current tag {request.current_tag} does not match the current Cargo package tag "
                f"{package_tag}; update the package version or pass the matching current tag"
            )
            raise ValueError(msg)


def _run_archive_request(*, args: argparse.Namespace, paths: ArchivePaths, request: ResolvedArchiveRequest, repo_root: Path) -> ArchiveResult:
    _validate_generation_request(args=args, request=request, repo_root=repo_root)
    if args.generate_in_temp_worktree:
        _fetch_required_tags(request=request, repo_root=repo_root, include_current=args.github_assets)
        config = _generation_config(args=args, request=request, repo_root=repo_root)
        if args.output_only:
            return ArchiveResult(
                report_id=generate_worktree_report(
                    output=paths.output,
                    config=config,
                    artifacts=paths.artifacts,
                ),
                action="output",
            )
        return ArchiveResult(
            report_id=generate_and_promote_worktree_report(
                output=paths.output,
                current=paths.current,
                archive_dir=paths.archive_dir,
                config=config,
                artifacts=paths.artifacts,
            ),
            action="promote-generated",
        )

    if args.output_only:
        msg = "--output-only requires --generate-in-temp-worktree"
        raise ValueError(msg)
    if args.github_assets:
        msg = "--github-assets requires --generate-in-temp-worktree"
        raise ValueError(msg)
    _run_benchmark_input_gate(repo_root, env=_benchmark_env(repo_root))
    return ArchiveResult(
        report_id=promote_report(
            source=paths.source,
            current=paths.current,
            archive_dir=paths.archive_dir,
            expected_current_tag=request.current_tag,
            expected_baseline_tag=request.baseline_tag,
        ),
        action="promote-source",
    )


def _validate_cli_preflight(args: argparse.Namespace) -> None:
    """Reject locally invalid options before release discovery or tag fetching."""
    if args.promote_artifacts and any(
        (
            args.current_tag,
            args.baseline_tag,
            args.published_latest,
            args.infer_release,
            args.current_vs_latest,
            args.github_assets,
            args.generate_in_temp_worktree,
            args.output_only,
            args.local_report,
            args.suite,
            args.scope,
        )
    ):
        msg = "--promote-artifacts cannot be combined with release selection, generation, GitHub-asset, output-only, suite, or scope options"
        raise ValueError(msg)
    if args.output_only and not args.generate_in_temp_worktree:
        msg = "--output-only requires --generate-in-temp-worktree"
        raise ValueError(msg)
    if args.local_report and not args.output_only:
        msg = "--local-report requires --output-only"
        raise ValueError(msg)
    if args.local_report and args.github_assets:
        msg = "--local-report cannot be combined with --github-assets"
        raise ValueError(msg)
    if args.github_assets and not args.generate_in_temp_worktree:
        msg = "--github-assets requires --generate-in-temp-worktree"
        raise ValueError(msg)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    args = build_parser().parse_args(argv)
    root = Path.cwd()

    try:
        _validate_cli_preflight(args)
        paths = _resolve_cli_paths(root, args)
        if args.promote_artifacts:
            result = ArchiveResult(
                report_id=render_and_promote_artifacts(
                    artifacts=paths.artifacts,
                    output=paths.output,
                    current=paths.current,
                    archive_dir=paths.archive_dir,
                ),
                action="promote-artifacts",
            )
        else:
            request = resolve_archive_request(
                ArchiveRequestOptions(
                    current_tag=args.current_tag,
                    baseline_tag=args.baseline_tag,
                    published_latest=args.published_latest,
                    infer_release=args.infer_release,
                    current_vs_latest=args.current_vs_latest,
                    worktree_ref=args.worktree_ref,
                    repo_root=root,
                )
            )
            result = _run_archive_request(args=args, paths=paths, request=request, repo_root=root)
    except (
        ExceptionGroup,
        ExecutableNotFoundError,
        OSError,
        TypeError,
        ValueError,
        RuntimeError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ) as exc:
        print(f"archive-performance: {format_exception_diagnostics(exc)}", file=sys.stderr)
        return 1

    if result.action == "output":
        print(f"Generated benchmark report in a temporary worktree and wrote it to {paths.output}")
    elif result.action == "promote-generated":
        print(f"Generated benchmark report in a temporary worktree, wrote it to {paths.output}, and promoted it to {paths.current}")
    elif result.action == "promote-artifacts":
        print(f"Rendered {paths.output} from {paths.artifacts.csv} and promoted it to {paths.current}")
    else:
        print(f"Promoted {paths.source} to {paths.current}")
    print(f"Current performance report: {result.report_id.current_tag} vs {result.report_id.baseline_tag}")
    print(f"Archive directory: {paths.archive_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
