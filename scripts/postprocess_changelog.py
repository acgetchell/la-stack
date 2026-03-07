#!/usr/bin/env python3
"""Post-process a git-cliff generated CHANGELOG.md.

Applies lightweight markdown hygiene that is difficult to express in
Tera templates:

  1. Strip trailing blank lines (git-cliff Tera templates emit an extra
     trailing newline that triggers markdownlint MD012).

Usage:
    postprocess-changelog                     # default: CHANGELOG.md
    postprocess-changelog path/to/CHANGELOG.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def postprocess(path: Path) -> None:
    """Read *path*, apply hygiene fixes, and write it back."""
    text = path.read_text(encoding="utf-8")

    # 1. Strip trailing blank lines — keep exactly one trailing newline.
    text = text.rstrip("\n") + "\n"

    path.write_text(text, encoding="utf-8")


def main() -> None:
    """CLI entry point for ``postprocess-changelog``."""
    parser = argparse.ArgumentParser(
        prog="postprocess-changelog",
        description="Apply markdown hygiene to a git-cliff generated CHANGELOG.md.",
    )
    parser.add_argument(
        "path",
        nargs="?",
        default="CHANGELOG.md",
        help="Path to CHANGELOG.md (default: CHANGELOG.md)",
    )
    args = parser.parse_args()

    changelog = Path(args.path)
    if not changelog.is_file():
        print(f"Error: {changelog} not found", file=sys.stderr)
        sys.exit(1)

    postprocess(changelog)


if __name__ == "__main__":
    main()
