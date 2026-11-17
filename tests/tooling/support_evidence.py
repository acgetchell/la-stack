"""Independent release evidence and provider fixtures shared by consumer tests."""

import hashlib
import io
import json
import os
import runpy
import shlex
import subprocess
import sys
from contextlib import redirect_stderr
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
from research_repo_tools import process, release_assets, release_evidence, release_interfaces, release_validation
from research_repo_tools.cli import main, parser
from research_repo_tools.config import load
from research_repo_tools.release_interfaces import load_evidence_policy

REPO_ROOT = Path(__file__).resolve().parents[2]
REPOSITORY = "acgetchell/la-stack"
COMMIT = "a" * 40
DIGEST = "sha256:" + "b" * 64
RUN_ENDPOINT = "actions/runs/101"
JOBS_ENDPOINT = "actions/runs/101/attempts/2/jobs?per_page=100&page=1"
ARTIFACTS_ENDPOINT = "actions/runs/101/artifacts?per_page=100&page=1"
ARTIFACT_ENDPOINT = "actions/artifacts/301"
RUNS_ENDPOINT = f"actions/workflows/9/runs?event=push&branch=main&head_sha={COMMIT}&per_page=100&page=1"
INVOCATION = f"https://github.com/{REPOSITORY}/actions/runs/101/attempts/2"
POLICY = load_evidence_policy(REPO_ROOT / "tooling/release-evidence.toml")
NATIVE_GATES = (
    ("prepare (ubuntu-latest)", ("Set up declared tools", "Run CI policy checks"), 201),
    ("prepare (macos-latest)", ("Set up declared tools", "Run CI policy checks"), 202),
    ("prepare (windows-latest)", ("Set up declared tools", "Run CI policy checks"), 203),
    ("release-package", ("Verify default and exact Cargo packages",), 204),
    *(
        (f"native-{stage} ({platform})", ("Set up declared tools", "Run native validation"), 206 + platform_index * 4 + stage_index)
        for platform_index, platform in enumerate(("ubuntu-latest", "macos-latest", "windows-latest"))
        for stage_index, stage in enumerate(("rust", "python", "candidate", "dependencies"))
    ),
)
NATIVE_RUNNERS = {name: [name.rsplit(" (", 1)[1].removesuffix(")")] for name, _, _ in NATIVE_GATES if name != "release-package"}
SIGNING_STEPS = (
    "Capture successful native validation identities",
    "Sign original main validation evidence",
    "Upload signed validation evidence",
)
REQUIRED_STEPS = (
    *((index, step) for index, (_, steps, _) in enumerate(NATIVE_GATES) for step in range(len(steps))),
    *((len(NATIVE_GATES), index) for index in range(3)),
)


def install_clean_checkout(patch: pytest.MonkeyPatch, checkout: Path, commit: str) -> None:
    """Supply read-only Git responses to the installed checkout verifier."""
    options = ("--no-pager", "--no-optional-locks", "--no-replace-objects", "-c", "core.fsmonitor=false", "-c", "diff.autoRefreshIndex=false")

    def git(executable: str | Path, args: list[str], *, cwd: Path, timeout: int, **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        assert cwd.samefile(checkout)
        assert timeout == 30
        assert Path(executable).is_absolute()
        assert tuple(args[: len(options)]) == options
        command = tuple(args[len(options) :])
        responses = {
            ("rev-parse", "HEAD"): f"{commit}\n".encode("ascii"),
            ("ls-files", "-v", "-z"): b"H Cargo.toml\0",
            ("submodule", "foreach", "--quiet", "--recursive", shlex.join((Path(executable).as_posix(), *options, "ls-files", "-v", "-z"))): b"",
            ("status", "--porcelain", "-z", "--untracked-files=all", "--ignore-submodules=none"): b"",
        }
        assert command in responses, f"unexpected read-only Git command: {command}"
        return subprocess.CompletedProcess(args, 0, stdout=responses[command], stderr=b"")

    patch.setattr(release_assets, "run_command_bytes", git)


class EvidenceProvider:
    """Keep signed inputs independent from live provider responses."""

    def __init__(self, root: Path) -> None:
        """Create independent signed and provider identities for one CI attempt."""
        self.root = root
        self.run: dict[str, Any] = {
            "id": 101,
            "run_attempt": 2,
            "workflow_id": 9,
            "path": ".github/workflows/ci.yml",
            "name": "CI",
            "event": "push",
            "head_branch": "main",
            "head_sha": COMMIT,
            "pull_requests": [],
            "repository": {"id": 11, "full_name": REPOSITORY},
            "head_repository": {"id": 11, "full_name": REPOSITORY},
            "status": "completed",
            "conclusion": "success",
        }
        gates = [*NATIVE_GATES, ("release-evidence", SIGNING_STEPS, 205)]
        self.jobs: list[dict[str, Any]] = [
            {
                "id": identifier,
                "name": name,
                "run_id": 101,
                "run_attempt": 2,
                "head_branch": "main",
                "head_sha": COMMIT,
                "status": "completed",
                "conclusion": "success",
                "labels": NATIVE_RUNNERS.get(name, ["ubuntu-latest"]).copy(),
                "steps": [{"number": number, "name": step, "status": "completed", "conclusion": "success"} for number, step in enumerate(steps, 1)],
            }
            for name, steps, identifier in gates
        ]
        self.artifact: dict[str, Any] = {
            "id": 301,
            "name": "release-evidence-101-2",
            "expired": False,
            "digest": DIGEST,
            "size_in_bytes": 512,
            "workflow_run": {"id": 101, "head_sha": COMMIT, "head_branch": "main", "repository_id": 11, "head_repository_id": 11},
        }
        self.document: dict[str, Any] = {
            "schema": "research-repo-tools/release-evidence/v1",
            "policy": {
                "repository": REPOSITORY,
                "workflow": ".github/workflows/ci.yml",
                "workflow_name": "CI",
                "ref": "refs/heads/main",
                "event": "push",
                "required_steps": {name: list(steps) for name, steps, _ in NATIVE_GATES},
                "runner_labels": deepcopy(NATIVE_RUNNERS),
                "artifact_prefix": "release-evidence",
                "guarantee": "signed-source",
                "signing_job": {"release-evidence": list(SIGNING_STEPS)},
                "path_in_vcs": "",
            },
            "run": {
                "id": 101,
                "run_attempt": 2,
                "head_sha": COMMIT,
                "head_branch": "main",
                "workflow_id": 9,
                "path": ".github/workflows/ci.yml",
                "name": "CI",
                "event": "push",
                "repository_id": 11,
                "head_repository_id": 11,
            },
            "jobs": {name: identifier for name, _, identifier in NATIVE_GATES},
            "package": {"name": "la-stack", "version": "0.5.42"},
        }
        self.signature = [{"verificationResult": {"signature": {"certificate": {"runInvocationURI": INVOCATION, "buildTrigger": "push"}}}}]
        self.responses: dict[str, Any] = {
            "actions/workflows/ci.yml": {"id": 9, "path": ".github/workflows/ci.yml", "name": "CI", "state": "active"},
            RUNS_ENDPOINT: {"workflow_runs": [self.run]},
            RUN_ENDPOINT: self.run,
            JOBS_ENDPOINT: {"jobs": self.jobs},
            ARTIFACTS_ENDPOINT: {"artifacts": [self.artifact]},
            ARTIFACT_ENDPOINT: self.artifact,
        }
        self.scripted: dict[str, list[Any]] = {}
        self.api_calls: list[str] = []
        self.gh_calls: list[tuple[str, ...]] = []
        self.downloads: list[tuple[str, dict[str, Any], Path, tuple[str, ...]]] = []
        self.document_text: str | None = None
        self.signature_text: str | None = None

    def api(self, endpoint: str, *, repository: str) -> Any:
        """Return an isolated response while recording authenticated requests."""
        assert repository == REPOSITORY
        self.api_calls.append(endpoint)
        if endpoint in self.scripted:
            return deepcopy(self.scripted[endpoint].pop(0))
        assert endpoint in self.responses, f"unexpected authenticated request: {endpoint}"
        return deepcopy(self.responses[endpoint])

    def git(self, args: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        """Supply the exact clean checkout identity to the subprocess verifier."""
        assert "rev-parse" in args or "status" in args
        assert kwargs["cwd"] == self.root
        assert kwargs["timeout"] == 30
        return subprocess.CompletedProcess(args, 0, stdout=COMMIT + "\n" if "rev-parse" in args else "", stderr="")

    def gh(self, *args: str) -> str:
        """Return the selected attestation response without a provider process."""
        assert args[:2] == ("attestation", "verify")
        self.gh_calls.append(args)
        return self.signature_text if self.signature_text is not None else json.dumps(self.signature)

    def download(self, repository: str, artifact: dict[str, Any], directory: Path, names: tuple[str, ...], **bounds: int) -> None:
        """Write the independently authored evidence and signature fixture."""
        # ZIP authentication and extraction are covered by the shared package.
        # This boundary supplies only the adapter's two extracted input files.
        self.downloads.append((repository, deepcopy(artifact), directory, names))
        assert names == ("validation.json", "validation.sigstore.json")
        assert bounds == {"max_bytes": 32 * 1024 * 1024, "max_expanded_bytes": 32 * 1024 * 1024}
        directory.mkdir()
        payload = self.document_text if self.document_text is not None else json.dumps(self.document)
        (directory / "validation.json").write_text(payload, encoding="utf-8")
        (directory / "validation.sigstore.json").write_text('{"fixture":"signed bundle"}', encoding="utf-8")

    def verify(self, identity: str | None = None) -> dict[str, Any]:
        """Exercise the public verifier with independently authored signed data."""
        return json.loads(
            release_evidence.verify_release_evidence(
                POLICY,
                "0.5.42",
                COMMIT,
                lambda endpoint: self.api(endpoint, repository=REPOSITORY),
                expected=identity.encode() if identity is not None else None,
            )
        )

    def seal(self, output: Path) -> None:
        """Exercise the installed producer including its environment/output guards."""
        args = parser().parse_args(
            [
                "release",
                "evidence",
                "create",
                "tooling/release-evidence.toml",
                "--checkout",
                str(self.root),
                "--commit",
                COMMIT,
                "--version",
                "0.5.42",
                "--run-id",
                "101",
                "--attempt",
                "2",
                "--output",
                str(output),
            ]
        )
        release_interfaces.run_evidence(args, load(root=self.root))

    def install(self, patch: pytest.MonkeyPatch) -> None:
        """Fake authenticated services and Git observations, leaving policy/CLI real."""
        patch.setattr(release_interfaces, "api", self.api)
        patch.setattr(release_validation, "gh", self.gh)
        patch.setattr(release_evidence, "download_artifact", self.download)
        install_clean_checkout(patch, self.root, COMMIT)

    def cli(self, *args: str) -> int:
        """Run the installed verifier with the fixture's exact source selection."""
        return main(
            [
                "--root",
                str(self.root),
                "release",
                "evidence",
                "verify",
                "tooling/release-evidence.toml",
                "--checkout",
                str(self.root),
                "--commit",
                COMMIT,
                "--version",
                "0.5.42",
                "--output",
                str(self.root / "selection.json"),
                *args,
            ]
        )


def expected_identity() -> str:
    """Independent approval record for the fixture's literal signed bytes."""
    provider = EvidenceProvider(Path())
    selected = {
        **provider.document,
        "schema": "research-repo-tools/release-selection/v1",
        "jobs": {**provider.document["jobs"], "release-evidence": 205},
        "artifact": {key: provider.artifact[key] for key in ("id", "name", "digest", "size_in_bytes")},
        "manifest_sha256": hashlib.sha256(json.dumps(provider.document).encode()).hexdigest(),
    }
    return json.dumps(selected, sort_keys=True, separators=(",", ":"))


PIN_JSON = expected_identity()


def run_evidence_fixture() -> int:
    """Run the real verifier in a child process with simulated provider responses."""
    root = Path(os.environ["FIXTURE_ROOT"])
    provider = EvidenceProvider(root)
    if os.environ.get("FIXTURE_FAILURE") == "failed-job":
        provider.jobs[3]["steps"][0]["conclusion"] = "failure"
    newline = os.environ.get("FIXTURE_DIAGNOSTIC_NEWLINE", os.linesep)
    sys.stdout.buffer.write(("fixture stdout diagnostic" + newline).encode("utf-8"))
    sys.stdout.flush()
    print("fixture stderr diagnostic", file=sys.stderr)
    diagnostic = io.StringIO()
    script, *arguments = sys.argv[1:]
    with pytest.MonkeyPatch.context() as patch:
        provider.install(patch)
        patch.setattr(process, "run_git_command", provider.git)
        # The uv launcher needs the real test environment; provenance observes
        # only the dependency-free release fixture after process startup.
        patch.delenv("UV_PROJECT", raising=False)
        patch.delenv("UV_PROJECT_ENVIRONMENT", raising=False)
        with redirect_stderr(diagnostic):
            try:
                sys.argv[:] = [script, *arguments]
                runpy.run_path(script, run_name="__main__")
                status = 0
            except SystemExit as error:
                status = error.code
    assert isinstance(status, int)
    print(diagnostic.getvalue(), end="", file=sys.stderr)
    # Just's quiet mode hides stderr. Record the actual verifier outcome so an
    # import/shell failure cannot masquerade as the rejection under test.
    result = {"status": status, "diagnostic": diagnostic.getvalue(), "import_root": sys.path[0], "arguments": arguments}
    (root / "verifier-result.json").write_text(json.dumps(result), encoding="utf-8")
    return status
