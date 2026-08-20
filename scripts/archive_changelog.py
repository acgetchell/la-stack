#!/usr/bin/env -S uv run
"""Archive completed minor series from CHANGELOG.md into per-minor files.

Parses the full CHANGELOG.md (produced by git-cliff + postprocess-changelog)
into version blocks, groups them by minor series (X.Y), and writes:

  - ``docs/archive/changelog/X.Y.md`` for each completed minor series
  - A trimmed ``CHANGELOG.md`` containing only the preamble, Unreleased,
    the active minor series, and an Archives link section

The active minor is detected from the first tagged release heading after
Unreleased.  All other minors are archived.

Usage:
    archive-changelog                      # default: CHANGELOG.md
    archive-changelog path/to/CHANGELOG.md
    archive-changelog --archive-dir docs/archive/changelog
"""

import argparse
import logging
import os
import re
import sys
import tempfile
from pathlib import Path

from postprocess_changelog import normalize_entry_headings_text, postprocess_text

# Matches ``## [X.Y.Z]`` or ``## [Unreleased]``
_VERSION_HEADING_RE = re.compile(r"^## \[")
_UNRELEASED_HEADING_RE = re.compile(r"^## \[Unreleased\](?:\s|$)")

# Extracts a strict SemVer 2.0.0 version from a ``## [X.Y.Z]`` heading.
_SEMVER_ALNUM_ID = r"(?:(?=[0-9A-Za-z-]*[A-Za-z-])[0-9A-Za-z-]+)"
_SEMVER_PATTERN = (
    r"(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)"
    rf"(?:-(?:(?:0|[1-9]\d*)|{_SEMVER_ALNUM_ID})(?:\.(?:(?:0|[1-9]\d*)|{_SEMVER_ALNUM_ID}))*)?"
    r"(?:\+[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?"
)
_VERSION_RE = re.compile(rf"^## \[({_SEMVER_PATTERN})\](?:\s|$)")

# Matches a reference-style link definition: ``[label]: URL``
_LINK_DEF_RE = re.compile(r"^\[([^\]]+)\]:\s+\S+")

# Archive directory relative to the repository root.
_DEFAULT_ARCHIVE_DIR = "docs/archive/changelog"

LOGGER = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def _minor_key(version: str) -> str:
    """Return the ``X.Y`` minor key for a semver version string.

    Parameters:
        version: A version string like ``0.7.2`` or ``1.2.3-rc.1``.

    Returns:
        The first two numeric components joined by a dot (e.g. ``0.7``).

    Raises:
        ValueError: If *version* does not contain at least two dot-separated components.
    """
    parts = version.split(".")
    if len(parts) < 2:
        msg = f"Expected a version with at least two components (X.Y), got: {version!r}"
        raise ValueError(msg)
    return f"{parts[0]}.{parts[1]}"


def _version_sort_key(label: str) -> tuple[bool, tuple[int, ...], tuple[tuple[int, int | str], ...]]:
    """Return a sort key for a version label that orders by semantic version.

    Non-numeric labels (e.g. ``unreleased``) sort after all numeric versions.
    Numeric parts are compared as integers so that ``0.10`` sorts after ``0.9``.

    Parameters:
        label: A version label like ``0.7.2``, ``0.10``, or ``unreleased``.

    Returns:
        A tuple suitable for use as a sort key.
    """
    label_without_build = label.split("+", 1)[0]
    core, separator, prerelease = label_without_build.partition("-")
    parts = core.split(".")
    try:
        nums = tuple(int(p) for p in parts)
    except ValueError:
        # Non-numeric labels ("unreleased") sort last (True > False).
        return (True, (), ())

    if not separator:
        prerelease_key: tuple[tuple[int, int | str], ...] = ((2, ""),)
    else:
        prerelease_key = tuple((0, int(part)) if part.isdecimal() else (1, part) for part in prerelease.split("."))

    return (False, nums, prerelease_key)


def _extract_link_defs(text: str) -> tuple[str, dict[str, str]]:
    """Separate trailing reference-style link definitions from changelog text.

    git-cliff appends reference-style link definitions at the bottom of
    CHANGELOG.md for every version heading.  When the changelog is split
    into per-version blocks these definitions must be distributed to the
    correct output files so that headings like ``## [0.7.2]`` resolve and
    no unused definitions trigger rumdl MD053.

    Parameters:
        text: The full changelog text.

    Returns:
        A 2-tuple of (*cleaned_text*, *link_defs*) where *link_defs* maps
        lowercase labels to their full definition lines.
    """
    lines = text.rstrip("\n").split("\n")
    link_defs: dict[str, str] = {}

    # Walk backwards from the end, collecting link-def and blank lines.
    i = len(lines) - 1
    while i >= 0:
        line = lines[i]
        m = _LINK_DEF_RE.match(line)
        if m:
            link_defs[m.group(1).lower()] = line
            i -= 1
        elif line.strip() == "":
            i -= 1
        else:
            break

    cleaned = "\n".join(lines[: i + 1])
    return cleaned.rstrip("\n") + "\n", link_defs


def parse_changelog(text: str) -> tuple[str, str, list[tuple[str, str]]]:
    """Split a full changelog into preamble, unreleased block, and version blocks.

    Parameters:
        text: The full contents of CHANGELOG.md.

    Returns:
        A 3-tuple of (preamble, unreleased_block, version_blocks). The
        ``unreleased_block`` is the complete ``## [Unreleased]`` block,
        including its heading. Each item in ``version_blocks`` is a
        ``(semver_label, full_heading_block)`` pair in the order it appears
        (newest first), where ``semver_label`` is only the parsed version text
        (for example ``"0.7.2"``) and ``full_heading_block`` still includes
        the raw ``## [...]`` heading and body.
    """
    lines = text.split("\n")

    # Locate all ``## [`` headings.
    headings: list[int] = []
    for i, line in enumerate(lines):
        if _VERSION_HEADING_RE.match(line):
            headings.append(i)

    if not headings:
        return text, "", []

    preamble = "\n".join(lines[: headings[0]])

    unreleased = ""
    unreleased_line: int | None = None
    version_blocks: list[tuple[str, str]] = []
    version_lines: dict[str, int] = {}

    for idx, start in enumerate(headings):
        end = headings[idx + 1] if idx + 1 < len(headings) else len(lines)
        block = "\n".join(lines[start:end])

        heading_line = lines[start]
        if _UNRELEASED_HEADING_RE.match(heading_line):
            if unreleased_line is not None:
                msg = f"Duplicate Unreleased changelog heading at lines {unreleased_line} and {start + 1}"
                raise ValueError(msg)
            unreleased = block
            unreleased_line = start + 1
        else:
            m = _VERSION_RE.match(heading_line)
            if not m:
                msg = f"Unrecognized changelog version heading at line {start + 1}: {heading_line!r}; expected '## [Unreleased]' or a semantic version"
                raise ValueError(msg)
            version = m.group(1)
            if version in version_lines:
                msg = f"Duplicate changelog version {version!r} at lines {version_lines[version]} and {start + 1}"
                raise ValueError(msg)
            version_lines[version] = start + 1
            version_blocks.append((version, block))

    return preamble, unreleased, version_blocks


def group_by_minor(
    version_blocks: list[tuple[str, str]],
) -> dict[str, list[tuple[str, str]]]:
    """Group version blocks by their ``X.Y`` minor key.

    Preserves insertion order (newest first within each minor).

    Parameters:
        version_blocks: List of ``(version, block_text)`` pairs.

    Returns:
        An ordered dict mapping minor keys to their version blocks.
    """
    groups: dict[str, list[tuple[str, str]]] = {}
    for ver, block in version_blocks:
        key = _minor_key(ver)
        groups.setdefault(key, []).append((ver, block))
    return groups


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def _format_link_defs(link_defs: dict[str, str], labels: set[str]) -> str:
    """Return the subset of *link_defs* whose labels are in *labels*.

    The definitions are returned in reverse-sorted order (matching the
    convention that git-cliff uses: ``[unreleased]`` first, then newest
    version to oldest).
    """
    relevant = [link_defs[label] for label in sorted(link_defs, key=_version_sort_key, reverse=True) if label in labels]
    return "\n".join(relevant) if relevant else ""


def write_archive(
    archive_dir: Path,
    minor: str,
    blocks: list[tuple[str, str]],
    link_defs: dict[str, str] | None = None,
) -> Path:
    """Write an archive file for a single minor series.

    Parameters:
        archive_dir: Directory for archive files.
        minor: The ``X.Y`` minor key.
        blocks: Version blocks belonging to this minor, newest first, using
            the ``(semver_label, full_heading_block)`` shape returned by
            ``parse_changelog``. The archive writer preserves each provided
            block verbatim after the generated archive title.
        link_defs: Optional mapping of lowercase labels to reference-style
            link definition lines.  Only definitions matching versions in
            *blocks* are included.

    Returns:
        The path of the written archive file.
    """
    path = archive_dir / f"{minor}.md"
    text = _build_archive_text(minor, blocks, link_defs)
    _publish_texts({path: text})
    return path


def _build_archive_text(
    minor: str,
    blocks: list[tuple[str, str]],
    link_defs: dict[str, str] | None = None,
) -> str:
    """Build one normalized archive payload without changing the filesystem."""
    parts = [f"# Changelog - {minor}.x\n"]
    for _ver, block in blocks:
        parts.append(block)

    text = "\n".join(parts)

    # Append only the reference-style link definitions for this archive.
    if link_defs:
        versions = {ver.lower() for ver, _ in blocks}
        defs_text = _format_link_defs(link_defs, versions)
        if defs_text:
            text = text.rstrip("\n") + "\n\n" + defs_text

    # Normalize archive output too; archived blocks can preserve historical
    # commit-body indentation that no longer appears in the trimmed root file.
    return postprocess_text(text)


def _stage_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            staged = Path(handle.name)
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        if staged is not None:
            staged.unlink(missing_ok=True)
        raise
    if staged is None:
        msg = "temporary changelog staging completed without a path"
        raise AssertionError(msg)
    return staged


def _stage_bytes(path: Path, payload: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            staged = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        if staged is not None:
            staged.unlink(missing_ok=True)
        raise
    if staged is None:
        msg = "temporary changelog restoration completed without a path"
        raise AssertionError(msg)
    return staged


def _replace_path(source: Path, destination: Path) -> None:
    source.replace(destination)


def _restore_text(path: Path, previous: bytes | None) -> None:
    if previous is None:
        path.unlink(missing_ok=True)
        return
    staged = _stage_bytes(path, previous)
    try:
        _replace_path(staged, path)
    finally:
        staged.unlink(missing_ok=True)


def _publish_texts(payloads: dict[Path, str]) -> None:
    """Publish a set of text files together and roll all of them back on failure."""
    if not payloads:
        return
    previous = {path: path.read_bytes() if path.is_file() else None for path in payloads}
    staged: dict[Path, Path] = {}
    replaced: list[Path] = []
    try:
        for path, text in payloads.items():
            staged[path] = _stage_text(path, text)
        for path, staged_path in staged.items():
            _replace_path(staged_path, path)
            replaced.append(path)
    except BaseException as publication_error:
        rollback_errors: list[BaseException] = []
        for path in reversed(replaced):
            try:
                _restore_text(path, previous[path])
            except BaseException as rollback_error:  # noqa: BLE001
                rollback_errors.append(rollback_error)
        if rollback_errors:
            group_message = "changelog publication and rollback failed"
            raise BaseExceptionGroup(
                group_message,
                [publication_error, *rollback_errors],
            ) from None
        raise
    finally:
        for staged_path in staged.values():
            staged_path.unlink(missing_ok=True)


def _postprocess_existing_archives(archive_dir: Path) -> None:
    """Normalize historical archive files that are not regenerated this run."""
    if not archive_dir.is_dir():
        return

    payloads: dict[Path, str] = {}
    for path in archive_dir.glob("*.md"):
        text = path.read_text(encoding="utf-8")
        normalized = normalize_entry_headings_text(text)
        if normalized != text:
            payloads[path] = normalized
    _publish_texts(payloads)


def _normalized_existing_archive_payloads(archive_dir: Path, *, excluded: set[Path]) -> dict[Path, str]:
    """Return changed historical archive payloads not regenerated this run."""
    if not archive_dir.is_dir():
        return {}
    payloads: dict[Path, str] = {}
    for path in archive_dir.glob("*.md"):
        if path in excluded:
            continue
        text = path.read_text(encoding="utf-8")
        normalized = normalize_entry_headings_text(text)
        if normalized != text:
            payloads[path] = normalized
    return payloads


def build_root(
    preamble: str,
    unreleased: str,
    active_blocks: list[tuple[str, str]],
    archived_minors: list[str],
    archive_dir_rel: str,
) -> str:
    """Assemble the trimmed root CHANGELOG.md content.

    Parameters:
        preamble: Text before the first ``## `` heading.
        unreleased: The full Unreleased block (empty string if absent).
        active_blocks: Version blocks for the active minor series.
        archived_minors: Sorted list of archived ``X.Y`` minor keys.
        archive_dir_rel: Relative path to the archive directory from the changelog file.

    Returns:
        The full text for the trimmed CHANGELOG.md.
    """
    parts: list[str] = [preamble]

    if unreleased:
        parts.append(unreleased)

    for _ver, block in active_blocks:
        parts.append(block)

    if archived_minors:
        # Build the Archives section.
        archive_lines = ["## Archives\n"]
        archive_lines.append("Older releases are archived by minor series:\n")
        archive_lines.extend(f"- [{minor}.x]({archive_dir_rel}/{minor}.md)" for minor in archived_minors)
        archive_lines.append("")
        parts.append("\n".join(archive_lines))

    return postprocess_text("\n".join(parts))


def _archive_dir_link_prefix(archive_dir: Path, changelog_parent: Path) -> str:
    """Return the Markdown link prefix from a changelog to its archive directory."""
    try:
        return archive_dir.relative_to(changelog_parent).as_posix()
    except ValueError:
        try:
            archive_dir_rel = Path(os.path.relpath(archive_dir, changelog_parent)).as_posix()
        except ValueError as err:
            msg = "cannot compute relative archive links because the archive and changelog directories are on different filesystem roots"
            raise ValueError(msg) from err
        if archive_dir_rel == ".." or archive_dir_rel.startswith("../") or Path(archive_dir_rel).is_absolute():
            LOGGER.warning(
                "Archive directory %s is outside changelog directory %s; generated Markdown links use %s",
                archive_dir,
                changelog_parent,
                archive_dir_rel,
            )
        return archive_dir_rel


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def archive_changelog(
    changelog_path: Path,
    archive_dir: Path | None = None,
) -> None:
    """Split a changelog into root + per-minor archive files.

    Parameters:
        changelog_path: Path to the full CHANGELOG.md.
        archive_dir: Directory for archive files.  Defaults to
            ``docs/archive/changelog`` relative to *changelog_path*'s parent.
    """
    if archive_dir is None:
        archive_dir = changelog_path.parent / _DEFAULT_ARCHIVE_DIR

    text = changelog_path.read_text(encoding="utf-8")

    # Separate trailing reference-style link definitions before parsing
    # so they can be distributed to the correct output files.
    text, link_defs = _extract_link_defs(text)

    preamble, unreleased, version_blocks = parse_changelog(text)

    if not version_blocks:
        _postprocess_existing_archives(archive_dir)
        return  # nothing to archive

    groups = group_by_minor(version_blocks)
    minor_keys = list(groups.keys())

    # Active minor = first minor that appears (newest release).
    active_minor = minor_keys[0]

    # Archive every minor except the active one.
    archived_minors = minor_keys[1:]

    if not archived_minors:
        _postprocess_existing_archives(archive_dir)
        return  # only one minor series — nothing to archive yet

    # Validate link portability before creating or modifying archive files.
    # In particular, os.path.relpath() cannot cross Windows volumes.
    archive_dir_rel = _archive_dir_link_prefix(archive_dir, changelog_path.parent)

    payloads = {archive_dir / f"{minor}.md": _build_archive_text(minor, groups[minor], link_defs) for minor in archived_minors}

    root_text = build_root(
        preamble,
        unreleased,
        groups[active_minor],
        sorted(archived_minors, key=_version_sort_key, reverse=True),
        archive_dir_rel,
    )

    # Append reference-style link definitions for active versions.
    if link_defs:
        labels: set[str] = {ver.lower() for ver, _ in groups[active_minor]}
        if unreleased:
            labels.add("unreleased")
        defs_text = _format_link_defs(link_defs, labels)
        if defs_text:
            root_text = root_text.rstrip("\n") + "\n\n" + defs_text + "\n"

    payloads[changelog_path] = root_text
    payloads.update(_normalized_existing_archive_payloads(archive_dir, excluded=set(payloads)))
    _publish_texts(payloads)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point for ``archive-changelog``."""
    parser = argparse.ArgumentParser(
        prog="archive-changelog",
        description="Archive completed minor series from CHANGELOG.md.",
    )
    parser.add_argument(
        "path",
        nargs="?",
        default="CHANGELOG.md",
        help="Path to CHANGELOG.md (default: CHANGELOG.md)",
    )
    parser.add_argument(
        "--archive-dir",
        default=None,
        help=f"Archive output directory (default: {_DEFAULT_ARCHIVE_DIR})",
    )
    args = parser.parse_args()

    changelog = Path(args.path)
    if not changelog.is_file():
        print(f"Error: {changelog} not found", file=sys.stderr)
        sys.exit(1)

    archive_dir = Path(args.archive_dir) if args.archive_dir else None
    archive_changelog(changelog, archive_dir)


if __name__ == "__main__":
    main()
