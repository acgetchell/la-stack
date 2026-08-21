"""Update deterministic release-version references from one target Git tag."""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path
from typing import TYPE_CHECKING

import check_docs_version_sync as version_sync
from archive_performance import published_stable_release_tags
from subprocess_utils import ExecutableNotFoundError

if TYPE_CHECKING:
    from collections.abc import Callable

_STABLE_TAG_RE = re.compile(r"^v(?P<major>0|[1-9][0-9]*)\.(?P<minor>0|[1-9][0-9]*)\.(?P<patch>0|[1-9][0-9]*)$")
_TOML_VERSION_RE = re.compile(r'^(?P<prefix>\s*version\s*=\s*")(?P<version>[^"]+)(?P<suffix>"\s*(?:#.*)?)$')
_CITATION_VERSION_RE = re.compile(
    r"^(?P<prefix>version:\s*(?P<quote>['\"]?))"
    r"(?P<version>[0-9A-Za-z][0-9A-Za-z.+-]*)"
    r"(?P<suffix>(?P=quote)\s*(?:#.*)?)$"
)
_CITATION_DATE_RE = re.compile(
    r"^(?P<prefix>date-released:\s*(?P<quote>['\"]?))"
    r"(?P<date>\d{4}-\d{2}-\d{2})"
    r"(?P<suffix>(?P=quote)\s*(?:#.*)?)$"
)
_BENCHMARK_TAG_PAIR_RE = re.compile(
    r"(?P<prefix>just performance-(?:github-assets|local-non-exact|release)[ \t]+)"
    r"v(?P<current>[0-9]+\.[0-9]+\.[0-9]+)"
    r"(?P<separator>[ \t]+)"
    r"v(?P<baseline>[0-9]+\.[0-9]+\.[0-9]+)"
    r"(?=[ \t]|`|$)",
    re.MULTILINE,
)
_BENCHMARK_BASELINE_PROSE_RE = re.compile(
    r"(?P<prefix>This generates a local `v)"
    r"(?P<baseline>[0-9]+\.[0-9]+\.[0-9]+)"
    r"(?P<suffix>` `vs_linalg` baseline)"
)


@dataclass(frozen=True, order=True, slots=True)
class ReleaseTag:
    """A stable release tag with SemVer ordering."""

    major: int
    minor: int
    patch: int
    tag: str = field(compare=False)

    @property
    def version(self) -> str:
        """Return the package version without the leading ``v``."""
        return self.tag.removeprefix("v")


@dataclass(frozen=True, slots=True)
class UpdateSummary:
    """Files and release identities produced by an update."""

    target: ReleaseTag
    previous: ReleaseTag
    changed_paths: tuple[Path, ...]
    release_date: str


@dataclass(frozen=True, slots=True)
class LineReplacement:
    """One fail-closed scalar replacement on a known source line."""

    line_number: int
    pattern: re.Pattern[str]
    group: str
    replacement: str
    allowed: frozenset[str]
    context: str


def parse_release_tag(value: str, *, label: str = "release tag") -> ReleaseTag:
    """Parse one stable ``vX.Y.Z`` release tag."""
    match = _STABLE_TAG_RE.fullmatch(value)
    if match is None:
        msg = f"{label} must be a stable tag in vX.Y.Z form, got {value!r}"
        raise ValueError(msg)
    return ReleaseTag(
        major=int(match.group("major")),
        minor=int(match.group("minor")),
        patch=int(match.group("patch")),
        tag=value,
    )


def select_previous_release_tag(tag_names: list[str], target: ReleaseTag) -> ReleaseTag:
    """Select the newest published stable release before *target*."""
    stable_tags = [parse_release_tag(tag) for tag in tag_names if _STABLE_TAG_RE.fullmatch(tag) is not None]
    if not stable_tags:
        msg = "repository has no published stable vX.Y.Z GitHub releases"
        raise ValueError(msg)
    newer = [tag for tag in stable_tags if tag > target]
    if newer:
        latest = max(newer)
        msg = f"target {target.tag} is older than published stable GitHub release {latest.tag}"
        raise ValueError(msg)
    previous = [tag for tag in stable_tags if tag < target]
    if not previous:
        msg = f"could not find a published stable GitHub release before {target.tag}"
        raise ValueError(msg)
    return max(previous)


def infer_previous_release_tag(root: Path, target: ReleaseTag) -> ReleaseTag:
    """Infer the previous release from published stable GitHub releases."""
    return select_previous_release_tag(published_stable_release_tags(root), target)


def _current_utc_date() -> str:
    """Return today's UTC calendar date for release preparation."""
    return datetime.now(UTC).date().isoformat()


def _validated_date(value: str) -> str:
    """Require one real ISO calendar date and return it unchanged."""
    try:
        parsed = date.fromisoformat(value)
    except ValueError as error:
        msg = f"release date must use YYYY-MM-DD, got {value!r}"
        raise ValueError(msg) from error
    if parsed.isoformat() != value:
        msg = f"release date must use canonical YYYY-MM-DD form, got {value!r}"
        raise ValueError(msg)
    return value


def _replace_line_group(text: str, edit: LineReplacement) -> str:
    lines = text.splitlines(keepends=True)
    if not 1 <= edit.line_number <= len(lines):
        msg = f"{edit.context} has no line {edit.line_number}"
        raise ValueError(msg)
    original_line = lines[edit.line_number - 1]
    body = original_line.rstrip("\r\n")
    ending = original_line[len(body) :]
    match = edit.pattern.fullmatch(body)
    if match is None:
        msg = f"{edit.context}:{edit.line_number} has an unsupported version assignment: {body!r}"
        raise ValueError(msg)
    current = match.group(edit.group)
    if current not in edit.allowed:
        msg = f"{edit.context}:{edit.line_number} has unexpected version {current!r}; expected one of {sorted(edit.allowed)}"
        raise ValueError(msg)
    start, end = match.span(edit.group)
    lines[edit.line_number - 1] = f"{body[:start]}{edit.replacement}{body[end:]}{ending}"
    return "".join(lines)


def _replace_match_groups(match: re.Match[str], replacements: dict[str, str]) -> str:
    updated = match.group(0)
    spans = sorted(((match.start(group) - match.start(), match.end(group) - match.start(), value) for group, value in replacements.items()), reverse=True)
    for start, end, value in spans:
        updated = f"{updated[:start]}{value}{updated[end:]}"
    return updated


def _replace_dependency_versions(text: str, package_name: str, target: ReleaseTag, previous: ReleaseTag, path: Path) -> str:
    allowed = frozenset({target.version, previous.version})
    pattern = version_sync.dependency_regex(package_name)

    def replace(match: re.Match[str]) -> str:
        group = "plain" if match.group("plain") is not None else "table"
        current = match.group(group)
        if current not in allowed:
            msg = f"{path} has unexpected {package_name} dependency version {current!r}; expected one of {sorted(allowed)}"
            raise ValueError(msg)
        return _replace_match_groups(match, {group: target.version})

    return pattern.sub(replace, text)


def _replace_readme_links(text: str, target: ReleaseTag, previous: ReleaseTag, path: Path) -> str:
    allowed = frozenset({target.version, previous.version})

    def replace(match: re.Match[str]) -> str:
        if version_sync.readme_tag_link_is_benchmark_asset(match):
            return match.group(0)
        version = match.group("version")
        if version is not None and version not in allowed:
            msg = f"{path} has unexpected release-pinned link version {version!r}; expected one of {sorted(allowed)}"
            raise ValueError(msg)
        group = "version" if version is not None else "revision"
        replacement = target.version if version is not None else target.tag
        return _replace_match_groups(match, {group: replacement})

    return version_sync.README_TAG_LINK_RE.sub(replace, text)


def _replace_benchmark_tag_pairs(text: str, target: ReleaseTag, previous: ReleaseTag, path: Path) -> str:
    allowed_current = frozenset({target.version, previous.version})

    def replace(match: re.Match[str]) -> str:
        current = match.group("current")
        if current not in allowed_current:
            msg = f"{path} has unexpected benchmark current tag v{current}; expected {target.tag} or {previous.tag}"
            raise ValueError(msg)
        return _replace_match_groups(match, {"current": target.version, "baseline": previous.version})

    return _BENCHMARK_TAG_PAIR_RE.sub(replace, text)


def _replace_benchmark_baseline_prose(text: str, previous: ReleaseTag) -> str:
    """Keep the active specific-release explanation aligned with its command."""
    return _BENCHMARK_BASELINE_PROSE_RE.sub(
        lambda match: _replace_match_groups(match, {"baseline": previous.version}),
        text,
    )


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _metadata_updates(root: Path, target: ReleaseTag, previous: ReleaseTag, release_date: str) -> dict[Path, str]:
    allowed = frozenset({target.version, previous.version})
    cargo_toml = root / "Cargo.toml"
    cargo_lock = root / "Cargo.lock"
    pyproject = root / "pyproject.toml"
    uv_lock = root / "uv.lock"
    citation = root / "CITATION.cff"

    package = version_sync.read_cargo_package_info(cargo_toml)
    project = version_sync.read_python_project_info(pyproject)
    cargo_toml_line = version_sync.toml_table_key_line(cargo_toml, "package", "version")
    cargo_lock_reference = version_sync.cargo_lock_reference(cargo_lock, package)
    pyproject_reference = version_sync.pyproject_reference(pyproject, project)
    uv_lock_reference = version_sync.uv_lock_reference(uv_lock, project)
    citation_reference = version_sync.citation_reference(citation)

    updates = {
        cargo_toml: _replace_line_group(
            _read_text(cargo_toml),
            LineReplacement(
                line_number=cargo_toml_line,
                pattern=_TOML_VERSION_RE,
                group="version",
                replacement=target.version,
                allowed=allowed,
                context=str(cargo_toml),
            ),
        ),
        cargo_lock: _replace_line_group(
            _read_text(cargo_lock),
            LineReplacement(
                line_number=cargo_lock_reference.line,
                pattern=_TOML_VERSION_RE,
                group="version",
                replacement=target.version,
                allowed=allowed,
                context=str(cargo_lock),
            ),
        ),
        pyproject: _replace_line_group(
            _read_text(pyproject),
            LineReplacement(
                line_number=pyproject_reference.line,
                pattern=_TOML_VERSION_RE,
                group="version",
                replacement=target.version,
                allowed=allowed,
                context=str(pyproject),
            ),
        ),
        uv_lock: _replace_line_group(
            _read_text(uv_lock),
            LineReplacement(
                line_number=uv_lock_reference.line,
                pattern=_TOML_VERSION_RE,
                group="version",
                replacement=target.version,
                allowed=allowed,
                context=str(uv_lock),
            ),
        ),
        citation: _replace_line_group(
            _read_text(citation),
            LineReplacement(
                line_number=citation_reference.line,
                pattern=_CITATION_VERSION_RE,
                group="version",
                replacement=target.version,
                allowed=allowed,
                context=str(citation),
            ),
        ),
    }
    citation_line, current_date = version_sync.citation_release_date(citation)
    updates[citation] = _replace_line_group(
        updates[citation],
        LineReplacement(
            line_number=citation_line,
            pattern=_CITATION_DATE_RE,
            group="date",
            replacement=release_date,
            allowed=frozenset({current_date, release_date}),
            context=str(citation),
        ),
    )
    return updates


def _prepare_updates(root: Path, target: ReleaseTag, previous: ReleaseTag, release_date: str) -> dict[Path, str]:
    updates = _metadata_updates(root, target, previous, release_date)
    changelog = root / "CHANGELOG.md"
    changelog_match = version_sync.changelog_release_date(changelog, target.version)
    if changelog_match is not None:
        changelog_line, changelog_date = changelog_match
        updates[changelog] = _replace_changelog_release_date(
            changelog,
            target,
            line=changelog_line,
            current_date=changelog_date,
            release_date=release_date,
        )
    package = version_sync.read_cargo_package_info(root / "Cargo.toml")
    for path in version_sync.iter_active_markdown_files(root):
        original = _read_text(path)
        updated = _replace_dependency_versions(original, package.name, target, previous, path)
        updated = _replace_benchmark_tag_pairs(updated, target, previous, path)
        updated = _replace_benchmark_baseline_prose(updated, previous)
        if path == root / "README.md":
            updated = _replace_readme_links(updated, target, previous, path)
        updates[path] = updated
    return updates


def _replace_changelog_release_date(
    changelog: Path,
    target: ReleaseTag,
    *,
    line: int,
    current_date: str,
    release_date: str,
) -> str:
    """Return a changelog with one target release heading date synchronized."""
    heading_re = re.compile(
        rf"^(?P<prefix>## \[v?{re.escape(target.version)}\] - )"
        r"(?P<date>\d{4}-\d{2}-\d{2})$"
    )
    return _replace_line_group(
        _read_text(changelog),
        LineReplacement(
            line_number=line,
            pattern=heading_re,
            group="date",
            replacement=release_date,
            allowed=frozenset({current_date, release_date}),
            context=str(changelog),
        ),
    )


def _write_text_atomic(path: Path, text: str) -> None:
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
        temporary.chmod(path.stat().st_mode)
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _validate_updated_root(root: Path, target: ReleaseTag, previous: ReleaseTag) -> None:
    mismatches = version_sync.find_version_mismatches(root)
    if mismatches:
        details = "; ".join(
            f"{mismatch.reference.path.relative_to(root)}:{mismatch.reference.line} has {mismatch.reference.version}" for mismatch in mismatches
        )
        msg = f"release-version validation failed after updating to {target.tag}: {details}"
        raise ValueError(msg)
    for path in version_sync.iter_active_markdown_files(root):
        for match in _BENCHMARK_TAG_PAIR_RE.finditer(_read_text(path)):
            if match.group("current") != target.version or match.group("baseline") != previous.version:
                msg = f"{path} contains a benchmark tag pair that does not match {target.tag} against {previous.tag}"
                raise ValueError(msg)
        for match in _BENCHMARK_BASELINE_PROSE_RE.finditer(_read_text(path)):
            if match.group("baseline") != previous.version:
                msg = f"{path} contains active baseline prose that does not match {previous.tag}"
                raise ValueError(msg)


def _publish_transaction(updates: dict[Path, str], validate: Callable[[], None]) -> tuple[Path, ...]:
    originals = {path: _read_text(path) for path in updates}
    changed = tuple(sorted((path for path, text in updates.items() if text != originals[path]), key=str))
    replaced: list[Path] = []
    try:
        for path in changed:
            _write_text_atomic(path, updates[path])
            replaced.append(path)
        validate()
    except BaseException as primary:
        rollback_errors: list[str] = []
        for path in reversed(replaced):
            try:
                _write_text_atomic(path, originals[path])
            except OSError as error:
                rollback_errors.append(f"{path}: {error}")
        if rollback_errors:
            msg = f"release-version update failed ({primary}); rollback also failed: {'; '.join(rollback_errors)}"
            raise RuntimeError(msg) from primary
        raise
    return changed


def update_release_version(
    root: Path,
    tag: str,
    *,
    previous: ReleaseTag | None = None,
    release_date: str | None = None,
) -> UpdateSummary:
    """Update release references transactionally and return a summary."""
    resolved_root = root.resolve()
    target = parse_release_tag(tag, label="target tag")
    previous_release = previous or infer_previous_release_tag(resolved_root, target)
    if previous_release >= target:
        msg = f"previous release {previous_release.tag} must be older than target {target.tag}"
        raise ValueError(msg)
    prepared_date = _validated_date(release_date or _current_utc_date())
    updates = _prepare_updates(resolved_root, target, previous_release, prepared_date)
    changed = _publish_transaction(updates, lambda: _validate_updated_root(resolved_root, target, previous_release))
    return UpdateSummary(target=target, previous=previous_release, changed_paths=changed, release_date=prepared_date)


def sync_changelog_release_date(
    root: Path,
    tag: str,
    *,
    previous: ReleaseTag | None = None,
) -> tuple[tuple[Path, ...], str]:
    """Synchronize a generated changelog heading from ``CITATION.cff``."""
    resolved_root = root.resolve()
    target = parse_release_tag(tag, label="target tag")
    previous_release = previous or infer_previous_release_tag(resolved_root, target)
    package = version_sync.read_cargo_package_info(resolved_root / "Cargo.toml")
    if package.version != target.version:
        msg = f"Cargo.toml version {package.version} does not match target {target.tag}"
        raise ValueError(msg)
    citation = resolved_root / "CITATION.cff"
    _, citation_date = version_sync.citation_release_date(citation)
    changelog = resolved_root / "CHANGELOG.md"
    changelog_match = version_sync.changelog_release_date(changelog, target.version)
    if changelog_match is None:
        msg = f"{changelog} has no generated release heading for {target.tag}"
        raise ValueError(msg)
    changelog_line, changelog_date = changelog_match
    updated = _replace_changelog_release_date(
        changelog,
        target,
        line=changelog_line,
        current_date=changelog_date,
        release_date=citation_date,
    )
    changed = _publish_transaction(
        {changelog: updated},
        lambda: _validate_updated_root(
            resolved_root,
            target,
            previous_release,
        ),
    )
    return changed, citation_date


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="Target stable release tag in vX.Y.Z form")
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Repository root to update (default: current directory)")
    parser.add_argument(
        "--sync-changelog-date",
        action="store_true",
        help="Synchronize the generated changelog heading from CITATION.cff instead of updating release metadata",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Update deterministic release metadata with fail-closed diagnostics."""
    args = parse_args(argv)
    try:
        if args.sync_changelog_date:
            changed_paths, release_date = sync_changelog_release_date(args.root, args.tag)
            if changed_paths:
                print(f"Synchronized CHANGELOG.md release date to {release_date}.")
            else:
                print(f"CHANGELOG.md release date already matches {release_date}.")
            return 0
        summary = update_release_version(args.root, args.tag)
    except (ExecutableNotFoundError, OSError, RuntimeError, subprocess.SubprocessError, TypeError, ValueError, tomllib.TOMLDecodeError) as error:
        print(f"failed to update release version: {error}", file=sys.stderr)
        return 1

    if summary.changed_paths:
        for path in summary.changed_paths:
            print(f"Updated {path.relative_to(args.root.resolve())}")
    else:
        print(f"Release-version references already match {summary.target.tag}.")
    print(f"Previous release: {summary.previous.tag}")
    print(f"CITATION.cff release date: {summary.release_date} (UTC update date)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
