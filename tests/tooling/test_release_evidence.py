"""Verify consumer evidence policy through the installed release interfaces."""

import hashlib
import json
import subprocess
from argparse import Namespace
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
from research_repo_tools import release_assets, release_evidence, release_interfaces, release_validation
from research_repo_tools.cli import main
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
    ("build (ubuntu-latest)", "Run CI checks", 201),
    ("build (macos-latest)", "Run CI checks", 202),
    ("build (windows-latest)", "Run CI checks", 203),
    ("release-package", "Verify default and exact Cargo packages", 204),
)
SIGNING_STEPS = (
    "Capture successful native validation identities",
    "Sign original main validation evidence",
    "Upload signed validation evidence",
)
REQUIRED_STEPS = (*((index, 0) for index in range(4)), *((4, index) for index in range(3)))


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
        gates = [*((name, (step,), identifier) for name, step, identifier in NATIVE_GATES), ("release-evidence", SIGNING_STEPS, 205)]
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
                "labels": [name.removeprefix("build (").removesuffix(")")] if name.startswith("build (") else ["ubuntu-latest"],
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
                "required_steps": {name: [step] for name, step, _ in NATIVE_GATES},
                "runner_labels": {name: [name.removeprefix("build (").removesuffix(")")] for name, _, _ in NATIVE_GATES[:3]},
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
        assert repository == REPOSITORY
        self.api_calls.append(endpoint)
        if endpoint in self.scripted:
            return deepcopy(self.scripted[endpoint].pop(0))
        assert endpoint in self.responses, f"unexpected authenticated request: {endpoint}"
        return deepcopy(self.responses[endpoint])

    def git(self, args: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        assert "rev-parse" in args or "status" in args
        assert kwargs["cwd"] == self.root
        assert kwargs["timeout"] == 30
        return subprocess.CompletedProcess(args, 0, stdout=COMMIT + "\n" if "rev-parse" in args else "", stderr="")

    def gh(self, *args: str) -> str:
        assert args[:2] == ("attestation", "verify")
        self.gh_calls.append(args)
        return self.signature_text if self.signature_text is not None else json.dumps(self.signature)

    def download(self, repository: str, artifact: dict[str, Any], directory: Path, names: tuple[str, ...]) -> None:
        # ZIP authentication and extraction are covered by the shared package.
        # This boundary supplies only the adapter's two extracted input files.
        self.downloads.append((repository, deepcopy(artifact), directory, names))
        assert names == ("validation.json", "validation.sigstore.json")
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
        args = Namespace(
            action="evidence",
            operation="create",
            configuration="tooling/release-evidence.toml",
            checkout=self.root,
            commit=COMMIT,
            version="0.5.42",
            run_id=101,
            attempt=2,
            output=output,
            github_output=None,
            tools=None,
            archive=None,
        )
        release_interfaces.run_evidence(args, load(root=self.root))

    def install(self, patch: pytest.MonkeyPatch) -> None:
        """Fake authenticated services and Git observations, leaving policy/CLI real."""
        patch.setattr(release_interfaces, "api", self.api)
        patch.setattr(release_validation, "gh", self.gh)
        patch.setattr(release_evidence, "download_artifact", self.download)
        patch.setattr(release_assets, "run_git_command", self.git)
        patch.setattr(release_assets, "run_git_bytes", lambda *args, **_kwargs: subprocess.CompletedProcess(args, 0, stdout=b"", stderr=b""))

    def cli(self, *args: str) -> int:
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


@pytest.fixture
def provider(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> EvidenceProvider:
    root = tmp_path / "checkout with spaces"
    root.mkdir()
    (root / "Cargo.toml").write_text('[package]\nname = "la-stack"\nversion = "0.5.42"\n', encoding="utf-8")
    (root / "tooling").mkdir()
    (root / "tooling/release-evidence.toml").write_bytes((REPO_ROOT / "tooling/release-evidence.toml").read_bytes())
    result = EvidenceProvider(root)
    result.install(monkeypatch)
    for name, value in {
        "GITHUB_SHA": COMMIT,
        "GITHUB_RUN_ID": "101",
        "GITHUB_RUN_ATTEMPT": "2",
        "GITHUB_REPOSITORY": REPOSITORY,
        "GITHUB_EVENT_NAME": "push",
        "GITHUB_REF": "refs/heads/main",
    }.items():
        monkeypatch.setenv(name, value)
    return result


def test_verify_consumes_exact_main_attempt_with_bound_attestation_flags(provider: EvidenceProvider) -> None:
    selection = provider.verify()
    assert selection["run"]["id"] == 101
    assert selection["artifact"]["digest"] == DIGEST
    assert provider.api_calls.count(JOBS_ENDPOINT) == 2
    assert len(provider.downloads) == 1
    repository, artifact, directory, names = provider.downloads[0]
    assert repository == REPOSITORY
    assert artifact == {"id": 301, "name": "release-evidence-101-2", "digest": DIGEST, "size_in_bytes": 512}
    assert names == ("validation.json", "validation.sigstore.json")
    assert provider.gh_calls == [
        (
            "attestation",
            "verify",
            str(directory / "validation.json"),
            "--bundle",
            str(directory / "validation.sigstore.json"),
            "--repo",
            REPOSITORY,
            "--signer-workflow",
            "acgetchell/la-stack/.github/workflows/ci.yml",
            "--signer-digest",
            COMMIT,
            "--source-digest",
            COMMIT,
            "--source-ref",
            "refs/heads/main",
            "--deny-self-hosted-runners",
            "--predicate-type",
            "https://slsa.dev/provenance/v1",
            "--format",
            "json",
        )
    ]
    assert not directory.exists()


def test_seal_captures_successful_native_jobs_while_its_signing_job_is_pending(provider: EvidenceProvider) -> None:
    provider.run.update(status="in_progress", conclusion=None)
    provider.jobs.pop()  # The sealing job cannot have completed its own upload.
    output = provider.root / "evidence" / "validation.json"
    provider.seal(output)
    assert json.loads(output.read_text(encoding="utf-8")) == provider.document
    assert provider.api_calls == [RUN_ENDPOINT, "actions/workflows/ci.yml", JOBS_ENDPOINT, RUN_ENDPOINT, "actions/workflows/ci.yml"]
    assert provider.downloads == []
    assert provider.gh_calls == []
    original = output.read_bytes()
    provider.seal(output)  # Identical immutable writes are safe retries.
    assert output.read_bytes() == original


def test_seal_rejects_environment_checkout_mismatch_before_creating_output(provider: EvidenceProvider, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GITHUB_SHA", "c" * 40)
    output = provider.root / "evidence" / "validation.json"
    with pytest.raises(ValueError, match="exact executing trusted push/run/attempt"):
        provider.seal(output)
    assert not output.parent.exists()
    assert provider.api_calls == []


def test_seal_rechecks_the_attempt_before_publishing_the_document(provider: EvidenceProvider) -> None:
    provider.scripted[RUN_ENDPOINT] = [deepcopy(provider.run), {**provider.run, "run_attempt": 3}]
    output = provider.root / "evidence" / "validation.json"
    with pytest.raises(ValueError, match="selected run or attempt changed"):
        provider.seal(output)
    assert not output.parent.exists()
    assert JOBS_ENDPOINT in provider.api_calls


@pytest.mark.parametrize("job_index", [0, 3], ids=["native", "package"])
def test_seal_never_writes_a_manifest_for_failed_required_validation(provider: EvidenceProvider, job_index: int) -> None:
    provider.jobs[job_index]["steps"][0]["conclusion"] = "failure"
    output = provider.root / "evidence" / "validation.json"
    with pytest.raises(ValueError, match="did not validate required step"):
        provider.seal(output)
    assert not output.parent.exists()


def test_pinned_recheck_uses_only_the_approved_run_attempt_and_artifact(provider: EvidenceProvider) -> None:
    identity = json.dumps(provider.verify())
    provider.api_calls.clear()
    del provider.responses[RUNS_ENDPOINT]
    assert provider.verify(identity) == json.loads(identity)
    assert provider.api_calls.count(RUN_ENDPOINT) == 2
    assert provider.api_calls.count(ARTIFACT_ENDPOINT) == 2
    assert len(provider.gh_calls) == 2


def test_pinned_recheck_rejects_an_attempt_replaced_by_a_rerun(provider: EvidenceProvider) -> None:
    selection = provider.verify()
    selection["run"]["run_attempt"] = 1
    provider.downloads.clear()
    identity = json.dumps(selection)
    with pytest.raises(ValueError, match="selected run or attempt changed"):
        provider.verify(identity)
    assert provider.downloads == []


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("event", "pull_request"),
        ("event", "workflow_dispatch"),
        ("head_branch", "release"),
        ("head_sha", "c" * 40),
        ("workflow_id", 10),
        ("path", ".github/workflows/other.yml"),
    ],
)
def test_verify_rejects_another_source_event_ref_commit_or_workflow(provider: EvidenceProvider, field: str, value: Any) -> None:
    provider.run[field] = value
    with pytest.raises(ValueError, match="exact commit's trusted main-push"):
        provider.verify()
    assert provider.downloads == []


def test_verify_rejects_a_fork_run_with_matching_commit(provider: EvidenceProvider) -> None:
    provider.run["head_repository"] = {"id": 12, "full_name": "fork/la-stack"}
    with pytest.raises(ValueError, match="repository or head repository differs"):
        provider.verify()
    assert provider.downloads == []


def test_verify_rejects_a_failed_run_without_using_successful_job_fragments(provider: EvidenceProvider) -> None:
    provider.run["conclusion"] = "failure"
    with pytest.raises(ValueError, match="missing or failed main-push validation"):
        provider.verify()
    assert JOBS_ENDPOINT not in provider.api_calls
    assert provider.downloads == []


@pytest.mark.parametrize("job_index", [0, 4], ids=["native", "signing"])
@pytest.mark.parametrize(("field", "value"), [("run_id", 102), ("run_attempt", 1), ("head_sha", "c" * 40), ("head_branch", "release")])
def test_job_identity_must_belong_to_the_exact_run_and_attempt(provider: EvidenceProvider, job_index: int, field: str, value: Any) -> None:
    provider.jobs[job_index][field] = value
    with pytest.raises(ValueError, match="belongs to another run, attempt, ref or commit"):
        provider.verify()
    assert provider.downloads == []


@pytest.mark.parametrize("job_index", [0, 4], ids=["native", "signing"])
def test_failed_jobs_cannot_be_replaced_by_successful_step_fragments(provider: EvidenceProvider, job_index: int) -> None:
    provider.jobs[job_index]["conclusion"] = "failure"
    with pytest.raises(ValueError, match="did not succeed"):
        provider.verify()
    assert provider.downloads == []


def test_native_job_identities_cannot_alias_one_another(provider: EvidenceProvider) -> None:
    provider.jobs[1]["id"] = provider.jobs[0]["id"]
    with pytest.raises(ValueError, match="validation job IDs are not unique"):
        provider.verify()
    assert provider.downloads == []


@pytest.mark.parametrize("job_index", range(5))
@pytest.mark.parametrize("duplicate", [False, True], ids=["missing", "ambiguous"])
def test_every_native_package_and_signing_job_is_unambiguous(provider: EvidenceProvider, job_index: int, duplicate: bool) -> None:
    name = provider.jobs[job_index]["name"]
    if duplicate:
        provider.jobs.append(deepcopy(provider.jobs[job_index]))
    else:
        provider.jobs.pop(job_index)
    with pytest.raises(ValueError, match="missing or ambiguous") as failure:
        provider.verify()
    assert name in str(failure.value)
    assert provider.downloads == []


@pytest.mark.parametrize(("job_index", "step_index"), REQUIRED_STEPS)
@pytest.mark.parametrize("conclusion", ["failure", "skipped"])
def test_required_native_package_and_signing_steps_must_succeed(provider: EvidenceProvider, job_index: int, step_index: int, conclusion: str) -> None:
    step = provider.jobs[job_index]["steps"][step_index]
    step["conclusion"] = conclusion
    with pytest.raises(ValueError, match="did not validate required step") as failure:
        provider.verify()
    assert step["name"] in str(failure.value)
    assert provider.downloads == []


@pytest.mark.parametrize("duplicate", [False, True], ids=["missing", "ambiguous"])
def test_required_step_names_cannot_be_substituted_or_duplicated(provider: EvidenceProvider, duplicate: bool) -> None:
    steps = provider.jobs[3]["steps"]
    if duplicate:
        steps.append(deepcopy(steps[0]))
    else:
        steps[0]["name"] = "Package check skipped"
    with pytest.raises(ValueError, match="did not validate required step Verify default and exact Cargo packages"):
        provider.verify()


@pytest.mark.parametrize(
    ("field", "value"),
    [("id", 102), ("head_sha", "c" * 40), ("head_branch", "release"), ("repository_id", 12), ("head_repository_id", 12)],
)
def test_artifact_association_must_match_the_validated_run(provider: EvidenceProvider, field: str, value: Any) -> None:
    provider.artifact["workflow_run"][field] = value
    with pytest.raises(ValueError, match="artifact is associated with another repository, ref, run or commit"):
        provider.verify()
    assert provider.downloads == []


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("expired", True, "missing, expired or mismatched"),
        ("digest", None, "authoritative SHA-256 digest"),
        ("digest", "sha256:not-a-digest", "authoritative SHA-256 digest"),
        ("size_in_bytes", 0, "artifact size exceeds the release bound"),
        ("size_in_bytes", 32 * 1024 * 1024 + 1, "artifact size exceeds the release bound"),
        ("id", True, "artifact ID must be a positive integer"),
    ],
)
def test_artifact_requires_unexpired_bounded_authoritative_identity(provider: EvidenceProvider, field: str, value: Any, message: str) -> None:
    provider.artifact[field] = value
    with pytest.raises(ValueError, match=message):
        provider.verify()
    assert provider.downloads == []


@pytest.mark.parametrize("inventory", ["missing", "ambiguous", "other-attempt"])
def test_artifact_selection_does_not_fall_back_or_choose_ambiguously(provider: EvidenceProvider, inventory: str) -> None:
    variants = {
        "missing": [],
        "ambiguous": [provider.artifact, deepcopy(provider.artifact)],
        "other-attempt": [{**provider.artifact, "name": "release-evidence-101-1"}],
    }
    provider.responses[ARTIFACTS_ENDPOINT] = {"artifacts": variants[inventory]}
    with pytest.raises(ValueError, match="missing or ambiguous evidence artifact"):
        provider.verify()
    assert provider.downloads == []


def test_pinned_digest_cannot_be_replaced_by_another_valid_digest(provider: EvidenceProvider) -> None:
    identity = json.dumps(provider.verify())
    provider.downloads.clear()
    provider.artifact["digest"] = "sha256:" + "c" * 64
    with pytest.raises(ValueError, match="approved artifact identity changed"):
        provider.verify(identity)
    assert provider.downloads == []


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("schema",), 2),
        (("run", "head_sha"), "c" * 40),
        (("run", "id"), 102),
        (("run", "run_attempt"), 1),
        (("run", "repository_id"), 12),
        (("package", "version"), "0.5.41"),
        (("jobs", "build (ubuntu-latest)"), 999),
    ],
)
def test_signed_manifest_must_equal_current_package_and_native_identities(provider: EvidenceProvider, path: tuple[str, ...], value: Any) -> None:
    target = provider.document
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError, match="retained evidence differs from current exact source"):
        provider.verify()
    assert provider.gh_calls == []
    assert not provider.downloads[0][2].exists()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("runInvocationURI", INVOCATION.replace("/101/", "/102/")),
        ("runInvocationURI", INVOCATION.replace("/attempts/2", "/attempts/1")),
        ("buildTrigger", "pull_request"),
        ("buildTrigger", "workflow_dispatch"),
    ],
)
def test_verified_attestation_must_name_this_push_invocation(provider: EvidenceProvider, field: str, value: str) -> None:
    provider.signature[0]["verificationResult"]["signature"]["certificate"][field] = value
    with pytest.raises(ValueError, match="attestation certificate run, attempt or expected event differs"):
        provider.verify()
    assert not provider.downloads[0][2].exists()


@pytest.mark.parametrize("source", ["pin", "manifest", "attestation"])
def test_duplicate_keys_cannot_replace_pinned_signed_or_verified_claims(provider: EvidenceProvider, source: str) -> None:
    identity = None
    if source == "pin":
        identity = '{"run": {"id": 101, "id": 101}}'
    elif source == "manifest":
        provider.document_text = json.dumps(provider.document).replace('"schema":', '"schema": 1, "schema":', 1)
    else:
        provider.signature_text = json.dumps(provider.signature).replace('"buildTrigger": "push"', '"buildTrigger": "push", "buildTrigger": "push"', 1)
    with pytest.raises(ValueError, match="duplicate JSON key"):
        provider.verify(identity)
    if source == "pin":
        assert provider.api_calls == []
    else:
        assert not provider.downloads[0][2].exists()


def test_successful_verification_rejects_a_rerun_started_during_consumption(provider: EvidenceProvider) -> None:
    provider.responses[RUN_ENDPOINT] = {**provider.run, "run_attempt": 3}
    with pytest.raises(ValueError, match="selected run or attempt changed"):
        provider.verify()
    assert len(provider.gh_calls) == 1
    assert not provider.downloads[0][2].exists()


def test_successful_verification_rechecks_immutable_artifact_identity(provider: EvidenceProvider) -> None:
    provider.responses[ARTIFACT_ENDPOINT] = {**provider.artifact, "digest": "sha256:" + "c" * 64}
    with pytest.raises(ValueError, match="artifact identity changed during verification"):
        provider.verify()
    assert provider.gh_calls == []
    assert provider.downloads == []


@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("timeout", [False, True])
def test_cli_preserves_provider_failure_diagnostics(
    provider: EvidenceProvider, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], binary: bool, timeout: bool
) -> None:
    command = ["gh", "attestation", "verify", "validation.json"]
    output = b"provider response" if binary else "provider response"
    diagnostic = b"artifact 301 expired" if binary else "artifact 301 expired"
    error = (
        subprocess.TimeoutExpired(command, 30, output=output, stderr=diagnostic)
        if timeout
        else subprocess.CalledProcessError(1, command, output=output, stderr=diagnostic)
    )

    def fail(*args: str) -> str:
        raise error

    monkeypatch.setattr(release_validation, "gh", fail)
    assert provider.cli() == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "artifact 301 expired" in captured.err
    assert "provider response" in captured.err
    assert not (provider.root / "selection.json").exists()


@pytest.mark.parametrize("failure", ["invalid-pin", "changed-pin", "failed-job", "signature", "late-artifact-change"])
def test_cli_verification_failure_never_exports_identity(provider: EvidenceProvider, capsys: pytest.CaptureFixture[str], failure: str) -> None:
    selection = provider.verify()
    output = provider.root / "github output"
    original = b"prior=value\n"
    output.write_bytes(original)
    expected = provider.root / "approved.json"
    args = ["--github-output", str(output)]
    if failure in {"invalid-pin", "changed-pin"}:
        selection["artifact"]["digest"] = "sha256:" + "c" * 64
        expected.write_text("{}" if failure == "invalid-pin" else json.dumps(selection), encoding="utf-8")
        args += ["--expected", str(expected)]
    elif failure == "failed-job":
        provider.jobs[3]["steps"][0]["conclusion"] = "failure"
    elif failure == "signature":
        provider.signature[0]["verificationResult"]["signature"]["certificate"]["buildTrigger"] = "workflow_dispatch"
    else:
        provider.responses[ARTIFACT_ENDPOINT] = {**provider.artifact, "digest": "sha256:" + "c" * 64}
    assert provider.cli(*args) == 1
    assert output.read_bytes() == original
    assert not (provider.root / "selection.json").exists()
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("index", range(3))
def test_native_job_requires_its_configured_host(provider: EvidenceProvider, index: int) -> None:
    provider.jobs[index]["labels"] = ["self-hosted"]
    with pytest.raises(ValueError, match="required hosted platform"):
        provider.verify()
    assert provider.downloads == []


def test_signed_manifest_retains_one_mib_limit(provider: EvidenceProvider) -> None:
    provider.document_text = " " * (1024 * 1024 + 1)
    with pytest.raises(ValueError, match=r"size|bound"):
        provider.verify()
    assert provider.gh_calls == []
