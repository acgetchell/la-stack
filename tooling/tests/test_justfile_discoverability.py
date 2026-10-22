"""Regression tests for the Just recipe surface."""

import json
import re
import shlex
import tomllib
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml
from research_repo_tools.cargo_examples import discover_examples
from research_repo_tools.just_inspect import dry_run, inspect_justfile
from research_repo_tools.process import run_command
from research_repo_tools.selection import select_files

REPO_ROOT = Path(__file__).resolve().parents[2]

if TYPE_CHECKING:
    import subprocess


def run_just(
    *args: str,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run the repository's installed Just executable without a shell."""
    return run_command("just", args, cwd=REPO_ROOT, check=check, env=env)


def just_recipes() -> dict[str, dict[str, Any]]:
    """Return parsed recipe metadata from the pinned Just executable."""
    return inspect_justfile(REPO_ROOT).recipes


def test_exact_benchmark_package_excludes_peer_libraries() -> None:
    # No dependency resolution: this also runs before Cargo caches exist in CI.
    metadata = json.loads(
        run_command(
            "cargo",
            ["metadata", "--locked", "--offline", "--no-deps", "--format-version", "1"],
            cwd=REPO_ROOT,
        ).stdout
    )
    packages = {package["name"]: package for package in metadata["packages"]}
    library = packages["la-stack"]
    comparison = packages["la-stack-comparison"]
    dependencies = {dependency["name"] for dependency in library["dependencies"]}

    assert metadata["workspace_default_members"] == [library["id"]]
    assert {"criterion", "num-rational"} <= dependencies
    assert not {"nalgebra", "faer", "la-stack-comparison"} & dependencies
    assert "exact" in {target["name"] for target in library["targets"]}
    assert "vs_linalg" not in {target["name"] for target in library["targets"]}
    assert {"nalgebra", "faer"} <= {dependency["name"] for dependency in comparison["dependencies"]}
    assert "vs_linalg" in {target["name"] for target in comparison["targets"]}


def test_release_runs_each_suite_once_and_reuses_peer_measurements() -> None:
    baseline = dry_run(REPO_ROOT, "bench-save-baseline", ["v0.4.5", "all"])
    # Execute the real recipe's shell branching while replacing only Cargo.
    runner = run_just("--evaluate", "_run").stdout.strip()
    script = 'cargo() { printf "%s\\n" "$*"; }\n' + baseline.stderr.replace(runner + " ", "")
    baseline_run = run_command("bash", ["--noprofile", "--norc", "-euc", script], cwd=REPO_ROOT)
    baseline_commands = [shlex.split(line) for line in baseline_run.stdout.splitlines()]
    recipes = just_recipes()
    current = recipes["bench-latest"]
    assert not current["body"]
    current_commands = []
    for dependency in current["dependencies"]:
        recipe = recipes[dependency["recipe"]]
        assert not recipe["dependencies"]
        commands = dry_run(REPO_ROOT, dependency["recipe"]).stderr
        current_commands.extend(shlex.split(line.removeprefix(runner + " ")) for line in commands.splitlines())

    assert len(baseline_commands) == len(current_commands) == 2
    baseline_comparison, baseline_exact = baseline_commands
    current_comparison, current_exact = current_commands
    for command in (baseline_comparison, current_comparison):
        assert command[command.index("-p") + 1] == "la-stack-comparison"
        assert command[command.index("--bench") + 1] == "vs_linalg"
    for command in (baseline_exact, current_exact):
        assert "-p" not in command
        assert command[command.index("--bench") + 1] == "exact"
    assert baseline_comparison[baseline_comparison.index("--") + 1 :] == ["--noplot", "--save-baseline", "v0.4.5"]
    assert current_comparison[current_comparison.index("--") + 1 :] == ["la_stack", "--noplot"]


def test_ci_enforces_full_python_fixture_lint_policy() -> None:
    """Canonical CI should lint fixtures without narrowing the Ruff configuration."""
    recipes = just_recipes()
    ci_dependencies = {dependency["recipe"] for dependency in recipes["ci"]["dependencies"]}
    python_body = json.dumps(recipes["python-check"]["body"])

    assert "python-check" in ci_dependencies
    assert "research-repo-tools python check" in python_body
    assert "--select" not in python_body
    inventory = set(select_files(REPO_ROOT, include=("*.py", "*.pyi")))
    assert {"tooling/release_evidence.py", "tooling/tests/test_release_evidence.py", "tests/semgrep/scripts/python_portability.py"} <= inventory
    assert not any(".venv" in Path(name).parts for name in inventory)


def test_example_policy_covers_cargo_targets_and_required_features() -> None:
    """Future Cargo examples join automatically; exact targets need explicit policy."""
    policy = tomllib.loads((REPO_ROOT / "tooling/examples.toml").read_text(encoding="utf-8"))
    assert policy["package"] == "la-stack"
    assert "include" not in policy
    assert "exclude" not in policy
    discovered = discover_examples(REPO_ROOT, package="la-stack")
    assert {example.name for example in discovered} >= {"vector_angles", "exact_solve_3x3", "rational_input_5x5"}
    overrides = policy["examples"]
    for example in discovered:
        assert set(example.required_features) <= set(overrides.get(example.name, {}).get("features", []))
    assert set(overrides) <= {example.name for example in discovered}
    assert "validation cargo-examples tooling/examples.toml" in dry_run(REPO_ROOT, "examples").stderr


def test_bare_just_lists_the_complete_documented_surface_after_source_sorting() -> None:
    """Discover new recipes from Just metadata without another maintained help list."""
    recipes = just_recipes()
    listed = run_just().stdout
    assert listed == run_just("--list").stdout
    assert recipes["_default"]["private"]
    assert "default" in recipes["_default"]["attributes"]
    assert "help-workflows" not in recipes
    public = {name: recipe for name, recipe in recipes.items() if not recipe["private"]}
    assert set(re.findall(r"^    ([a-z][a-z0-9-]*)", listed, re.MULTILINE)) == set(public)
    for name, recipe in public.items():
        assert recipe["doc"], name
        line = next(line for line in listed.splitlines() if re.match(rf"^    {re.escape(name)}(?:\s|$)", line))
        for parameter in recipe["parameters"]:
            assert parameter["name"] in line
    source = (REPO_ROOT / "justfile").read_text(encoding="utf-8")
    definitions = re.findall(r"^(?![^\n]*:=)([a-z_][a-z0-9_-]*)(?: [^:\n]*)?:", source, re.MULTILINE)
    assert definitions == sorted(definitions)


def test_dependabot_shared_caller_covers_dependencies_with_explicit_manual_review() -> None:
    """Every dependency file has an explicit automatic or manual review policy."""
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/dependabot-auto-merge.yml").read_bytes())
    events = workflow.get("on", workflow.get(True))
    assert set(events) == {"pull_request_target"}
    assert events["pull_request_target"]["branches"] == ["main"]
    assert workflow["permissions"] == {}
    assert len(workflow["jobs"]) == 1
    job = workflow["jobs"]["approve-and-enable-auto-merge"]
    path, revision = job["uses"].split("@")
    assert path == "acgetchell/research-repo-tools/.github/workflows/dependabot-approve.yml"
    assert re.fullmatch(r"[0-9a-f]{40}", revision)
    assert "steps" not in job
    assert "secrets" not in job
    assert job["permissions"] == {"contents": "write", "pull-requests": "write"}
    assert job["with"]["repository"] == "acgetchell/la-stack"
    policy = json.loads(job["with"]["policy"])
    assert set(policy) == {"cargo", "uv", "github_actions"}
    assert set(policy["cargo"]["files"]) == {"Cargo.toml", "Cargo.lock", "benches/comparison/Cargo.toml"}
    assert set(policy["uv"]["files"]) == {"pyproject.toml", "uv.lock"}
    actions = {
        path.relative_to(REPO_ROOT).as_posix()
        for directory in (REPO_ROOT / ".github/workflows", REPO_ROOT / ".github/actions")
        for path in directory.rglob("*")
        if path.suffix in {".yml", ".yaml"}
    }
    # Publication and signed evidence production retain maintainer review.
    manual_review = {".github/workflows/publish-crates.yml", ".github/workflows/ci.yml"}
    assert manual_review <= actions
    assert set(policy["github_actions"]["files"]) == actions - manual_review


def test_zizmor_workflow_separates_findings_from_sarif_export() -> None:
    """SARIF's successful exit cannot mask findings in the blocking audit."""
    recipes = just_recipes()
    assert "research-repo-tools zizmor check" in json.dumps(recipes["zizmor"]["body"])
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/zizmor.yml").read_bytes())
    steps = workflow["jobs"]["analyze"]["steps"]
    audit = next(step for step in steps if step.get("id") == "audit")
    sarif = next(step for step in steps if step.get("id") == "sarif")
    assert audit["run"] == "just zizmor --require-online"
    assert "continue-on-error" not in audit
    assert "just zizmor --require-online --format sarif" in sarif["run"]
    assert "!cancelled()" in sarif["if"]
    upload = next(step for step in steps if step.get("uses", "").startswith("github/codeql-action/upload-sarif@"))
    assert "steps.sarif.outcome == 'success'" in upload["if"]
    assert "github.actor != 'dependabot[bot]'" in upload["if"]


def test_composite_allowlist_agrees_with_repository_policy() -> None:
    """Pinned references need approval in both workflow and composite checks."""
    policy = json.loads((REPO_ROOT / ".github/actions-allowlist.json").read_bytes())
    approved = {pattern.removesuffix("@*") for pattern in policy["patterns_allowed"]}
    rules = yaml.safe_load((REPO_ROOT / "semgrep.yaml").read_bytes())["rules"]
    rule = next(rule for rule in rules if rule["id"] == "la-stack.github-actions.external-action-approved-allowlist")
    pattern = re.compile(rule["patterns"][0]["pattern-regex"])
    unapproved = {"actions/setup-python", "swatinem/rust-cache", "zizmorcore/zizmor-action", "unapproved/setup-uv"}
    candidates = approved | unapproved | {f"{action}/unapproved" for action in approved}
    for action in candidates:
        for prefix in ("        uses: ", "    - uses: "):
            reference = f"{prefix}{action}@{'a' * 40} # v1.0.0"
            assert bool(pattern.search(reference)) == (action not in approved), reference
    commands = dry_run(REPO_ROOT, "workflow-allowlist-check").stderr
    assert "--policy .github/actions-allowlist.json .github/workflows" in commands
    assert "--config semgrep.yaml .github/actions" in commands


def test_semgrep_sarif_upload_honors_only_in_source_suppressions() -> None:
    """GitHub receives actionable findings with their original metadata across every run."""
    workflow = yaml.safe_load((REPO_ROOT / ".github/workflows/semgrep-sarif.yml").read_bytes())
    steps = workflow["jobs"]["semgrep-sarif"]["steps"]
    scan = next(step for step in steps if step.get("id") == "semgrep")
    command = next(line.strip() for line in scan["run"].splitlines() if line.strip().startswith("jq "))
    arguments = shlex.split(command)
    actionable = {"ruleId": "actionable", "partialFingerprints": {"primaryLocationLineHash": "original"}}
    reviewed = {"ruleId": "reviewed", "suppressions": [{"kind": "inSource"}]}
    external = {"ruleId": "external", "suppressions": [{"kind": "external"}]}
    other_run = {"ruleId": "other-run", "locations": [{"physicalLocation": {"artifactLocation": {"uri": "other.yml"}}}]}
    report = {
        "version": "2.1.0",
        "runs": [
            {"tool": {"driver": {"name": "Semgrep OSS"}}, "results": [reviewed, actionable, external]},
            {"tool": {"driver": {"name": "Semgrep OSS"}}, "results": [other_run, reviewed]},
        ],
    }
    filtered = json.loads(run_command("jq", arguments[1:2], input=json.dumps(report)).stdout)
    assert filtered["version"] == report["version"]
    assert filtered["runs"][0]["results"] == [actionable, external]
    assert filtered["runs"][1]["results"] == [other_run]
    assert filtered["runs"][0]["tool"] == report["runs"][0]["tool"]
    upload = next(step for step in steps if step.get("uses", "").startswith("github/codeql-action/upload-sarif@"))
    assert arguments[2:] == ["semgrep-raw.sarif", ">", upload["with"]["sarif_file"]]
    gate = next(step for step in steps if step["name"] == "Fail on repository rule findings")
    assert gate["if"] == "steps.semgrep.outputs.exit_code != '0'"
    assert gate["run"] == "exit 1"
