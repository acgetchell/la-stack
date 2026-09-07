"""Check raw Markdown line lengths using Unicode characters in every locale."""

import argparse
import sys
from pathlib import Path

MAX_LINE_LENGTH = 160


def main() -> int:
    """Check UTF-8 files, preserving the Markdown recipe's table exemption."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()
    failed = False
    for path in args.files:
        try:
            with path.open(encoding="utf-8") as source:
                for line_number, line in enumerate(source, start=1):
                    # Text mode normalizes CRLF; whitespace still counts toward the limit.
                    line = line.removesuffix("\n")
                    if not line.startswith("|") and len(line) > MAX_LINE_LENGTH:
                        print(f"{path}:{line_number}: line length {len(line)} exceeds {MAX_LINE_LENGTH}", file=sys.stderr)
                        failed = True
        except (OSError, UnicodeError) as error:
            print(f"{path}: {error}", file=sys.stderr)
            failed = True
    if failed:
        print("Markdown raw line-length check failed.", file=sys.stderr)
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
