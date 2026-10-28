"""Keep v0.2.4 Semgrep setup/scans in the standard environment until #317.

Only scanner selection is overridden; the installed shared CLI owns all fixture
validation. Remove this adapter after research-repo-tools#165 ships in v0.2.5.
"""

import os
import sys
import sysconfig
from importlib.metadata import version
from pathlib import Path

from research_repo_tools import semgrep_tool
from research_repo_tools.cli import main as tooling_main
from research_repo_tools.process import resolve_executable


def main(argv: list[str] | None = None) -> int:
    """Run the shared CLI with the scanner installed beside this interpreter."""
    if version("research-repo-tools") != "0.2.4":
        msg = "Review/remove the v0.2.4 Semgrep override before upgrading (la-stack#317)."
        raise RuntimeError(msg)
    scripts = Path(sysconfig.get_path("scripts"))
    previous_setup = semgrep_tool.setup
    previous_version = semgrep_tool.VERSION
    previous_path = os.environ.get("PATH")

    def scanner(root: Path) -> Path:
        # Shared setup synchronizes dev before selecting its scanner.
        # The shared CLI prints this module constant after selecting the binary.
        vars(semgrep_tool)["VERSION"] = version("semgrep")
        return resolve_executable(scripts / ("semgrep.exe" if os.name == "nt" else "semgrep"))

    try:
        # The native launcher also resolves pysemgrep through PATH.
        os.environ["PATH"] = os.pathsep.join([str(scripts), *([previous_path] if previous_path else [])])
        vars(semgrep_tool)["setup"] = scanner
        return tooling_main(argv)
    finally:
        vars(semgrep_tool)["setup"] = previous_setup
        vars(semgrep_tool)["VERSION"] = previous_version
        if previous_path is None:
            os.environ.pop("PATH", None)
        else:
            os.environ["PATH"] = previous_path


if __name__ == "__main__":
    sys.exit(main())
