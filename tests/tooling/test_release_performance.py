"""Installed preparation saves independent assertions and rejects stale inputs."""

import os
import shlex
import stat
import subprocess
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

import pytest
from research_repo_tools import performance_release
from research_repo_tools.cli import main
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.performance_progress import ProgressReporter
from research_repo_tools.process import run_command, run_git_command
from research_repo_tools.releases import PublishedRelease

from tooling import release_candidate as candidate

ROOT = Path(__file__).resolve().parents[2]
INPUTS = "target/bench-reports/release-inputs.json"


@pytest.fixture
def consumer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "consumer with spaces"
    root.mkdir()
    (root / "Cargo.toml").write_text('[package]\nname = "fixture"\nversion = "0.5.2"\n', encoding="utf-8")
    (root / ".gitignore").write_bytes(b"/target/\n")
    for name in ("tooling/measurements/full.toml", "tooling/performance-release.toml"):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((ROOT / name).read_bytes())
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(tmp_path / "absent-gitconfig"))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    run_git_command(["--no-pager", "init", "--quiet"], cwd=root)
    for tag in ("v0.4.6", "v0.5.0", "v0.5.1", "head"):
        release = "0.5.2" if tag == "head" else tag.removeprefix("v")
        (root / "Cargo.toml").write_text(f'[package]\nname = "fixture"\nversion = "{release}"\n', encoding="utf-8", newline="\n")
        run_git_command(["--no-pager", "add", "Cargo.toml"], cwd=root)
        run_git_command(
            [
                "--no-pager",
                "-c",
                "user.name=Fixture",
                "-c",
                "user.email=fixture@example.invalid",
                "-c",
                "commit.gpgsign=false",
                "commit",
                "--allow-empty",
                "--quiet",
                "-m",
                tag,
            ],
            cwd=root,
        )
        if tag != "head":
            run_git_command(["--no-pager", "-c", "tag.gpgsign=false", "tag", tag], cwd=root)
    monkeypatch.chdir(root)
    return root


def settings(root: Path) -> performance_release.PerformanceRelease:
    return performance_release.load_performance_release(root, "tooling/performance-release.toml")


def releases(monkeypatch: pytest.MonkeyPatch, *tags: str) -> None:
    history = tuple(PublishedRelease(tag, datetime(2026, 10, index, tzinfo=UTC)) for index, tag in enumerate(tags, start=1))
    monkeypatch.setattr(performance_release, "published_releases", lambda _root, **_kwargs: history)


@pytest.mark.parametrize("published", [False, True])
def test_preparation_resolves_prospective_head_or_publication_order(
    consumer: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], published: bool
) -> None:
    current = "v0.5.1" if published else "v0.5.2"
    (consumer / "Cargo.toml").write_text(f'[package]\nname = "fixture"\nversion = "{current[1:]}"\n', encoding="utf-8")
    releases(monkeypatch, *(("v0.5.0", "v0.5.1") if published else ("v0.5.1", "v0.5.0")))
    before = run_git_command(["--no-pager", "show-ref"], cwd=consumer).stdout
    assert main(["performance", "prepare", "tooling/performance-release.toml"]) == 0
    assert f"Saved preparation: {INPUTS}" in capsys.readouterr().out
    record = performance_release.load_preparation(consumer, settings(consumer))
    assert record.identities == {
        "current": (current, run_git_command(["--no-pager", "rev-parse", "v0.5.1" if published else "HEAD"], cwd=consumer).stdout.strip()),
        "baseline": ("v0.5.0", run_git_command(["--no-pager", "rev-parse", "v0.5.0"], cwd=consumer).stdout.strip()),
    }
    assert run_git_command(["--no-pager", "show-ref"], cwd=consumer).stdout == before
    assert not (consumer / "target/bench-reports/full").exists()


def test_explicit_pair_skips_provider_and_failed_preparation_preserves_record(consumer: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(performance_release, "published_releases", lambda *_args, **_kwargs: pytest.fail("explicit pairs need no provider"))
    plan = replace(settings(consumer), mode="explicit", current="v0.5.1", baseline="v0.5.0")
    performance_release.prepare_performance(consumer, plan)
    original = (consumer / INPUTS).read_bytes()
    with pytest.raises(subprocess.CalledProcessError):
        performance_release.prepare_performance(consumer, replace(plan, current="v0.5.2"))
    assert (consumer / INPUTS).read_bytes() == original


@pytest.mark.parametrize("failure", ["empty-history", "provider", "missing-tag", "old-release"])
def test_failed_preparation_never_writes_partial_assertions(
    consumer: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], failure: str
) -> None:
    releases(monkeypatch, *(() if failure == "empty-history" else ("v0.4.6" if failure == "old-release" else "v0.5.1",)))
    if failure == "missing-tag":
        run_git_command(["--no-pager", "tag", "-d", "v0.5.1"], cwd=consumer)
    if failure == "provider":

        def unavailable(*_args: object, **_kwargs: object) -> tuple[PublishedRelease, ...]:
            raise subprocess.CalledProcessError(1, ["gh", "release", "list"], stderr="provider unavailable")

        monkeypatch.setattr(performance_release, "published_releases", unavailable)
    assert main(["performance", "prepare", "tooling/performance-release.toml"]) == 1
    assert not (consumer / INPUTS).exists()
    expected = {
        "empty-history": "a published stable release is required",
        "provider": "provider unavailable",
        "missing-tag": "refs/tags/v0.5.1^{commit}",
        "old-release": "requires v0.5.0 or newer",
    }
    assert expected[failure] in capsys.readouterr().err


@pytest.mark.parametrize("changed", ["head", "baseline", "declaration", "package"])
def test_readme_rejects_stale_preparation_before_publication(consumer: Path, monkeypatch: pytest.MonkeyPatch, changed: str) -> None:
    releases(monkeypatch, "v0.5.1")
    performance_release.prepare_performance(consumer, settings(consumer))
    if changed == "head":
        run_git_command(
            [
                "--no-pager",
                "-c",
                "user.name=Fixture",
                "-c",
                "user.email=fixture@example.invalid",
                "-c",
                "commit.gpgsign=false",
                "commit",
                "--allow-empty",
                "--quiet",
                "-m",
                "later",
            ],
            cwd=consumer,
        )
    elif changed == "baseline":
        run_git_command(["--no-pager", "tag", "--force", "v0.5.1", "HEAD"], cwd=consumer)
    elif changed == "declaration":
        path = consumer / "tooling/measurements/full.toml"
        path.write_bytes(path.read_bytes() + b"\n# changed\n")
    else:
        (consumer / "Cargo.toml").write_bytes(b'[package]\nname="fixture"\nversion="0.5.3"\n')
    monkeypatch.setattr(performance_release, "published_releases", lambda *_args, **_kwargs: pytest.fail("README must not rediscover"))
    with pytest.raises(ValueError, match="stale preparation"):
        performance_release.load_readme_plan(consumer, settings(consumer))


@pytest.mark.parametrize("record", [None, "{", '{"schema": 2}', '{"schema": 1, "assertions": []}'])
def test_missing_or_legacy_preparation_stops_before_readme_and_measurement(consumer: Path, record: str | None) -> None:
    if record is not None:
        path = consumer / INPUTS
        path.parent.mkdir(parents=True)
        path.write_text(record, encoding="utf-8")
    for command in ("readme", "release"):
        assert main(["performance", command, "tooling/performance-release.toml"]) == 1
    assert not (consumer / "target/bench-reports/full").exists()


@pytest.mark.parametrize("owner_kind", ["directory", "file"])
def test_fresh_acquisition_preserves_existing_owner_with_recovery_guidance(consumer: Path, monkeypatch: pytest.MonkeyPatch, owner_kind: str) -> None:
    releases(monkeypatch, "v0.5.1")
    performance_release.prepare_performance(consumer, settings(consumer))
    owner = consumer.parent / "independent owner"
    if owner_kind == "directory":
        owner.mkdir()
    retained = owner / "preparation.json" if owner_kind == "directory" else owner
    original = b"original owner bytes\r\n\xff"
    retained.write_bytes(original)
    monkeypatch.setattr(candidate, "check_release_dependencies", lambda _root: None)
    monkeypatch.setattr(candidate, "freeze", lambda *_args, **_kwargs: pytest.fail("existing owner must be rejected before freezing"))
    monkeypatch.setattr(candidate, "measure_release", lambda *_args, **_kwargs: pytest.fail("existing owner must not start acquisition"))

    with pytest.raises(ValueError, match="acquisition owner record already exists") as stopped:
        candidate.acquire(consumer, owner, "target/checkpoint", resume=False, progress=ProgressReporter(mode="none"))

    assert str(owner) in str(stopped.value)
    assert "new checkpoint/owner directory" in str(stopped.value)
    assert "performance-resume if the original checkpoint still exists" in str(stopped.value)
    assert retained.read_bytes() == original
    if owner_kind == "directory":
        assert list(owner.iterdir()) == [retained]
    assert not (consumer / "target/checkpoint").exists()


def test_fresh_acquisition_preserves_unrelated_directory_creation_error(consumer: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    releases(monkeypatch, "v0.5.1")
    performance_release.prepare_performance(consumer, settings(consumer))
    owner = consumer.parent / "independent owner"
    failure = FileExistsError("unrelated directory creation failure")

    def fail_mkdir(_self: Path, *_args: object, **_kwargs: object) -> None:
        raise failure

    monkeypatch.setattr(candidate, "check_release_dependencies", lambda _root: None)
    monkeypatch.setattr(Path, "mkdir", fail_mkdir)
    with pytest.raises(FileExistsError) as stopped:
        candidate.acquire(consumer, owner, "target/checkpoint", resume=False, progress=ProgressReporter(mode="none"))

    assert stopped.value is failure
    assert not owner.exists()


def test_release_prepares_before_measurement_and_readme_remains_separate() -> None:
    preparation = dry_run(ROOT, "performance-prepare").stderr
    assert 'performance prepare tooling/performance-release.toml "$@"' in preparation
    assert "git fetch" not in preparation
    measurement = dry_run(ROOT, "performance-release").stderr
    assert measurement.count("performance prepare tooling/performance-release.toml --fetch-tags --allow-git-mutations") == 1
    assert measurement.count("--release-intent") == 1
    assert measurement.count("tooling/release_candidate.py acquire") == 1
    assert measurement.index("performance prepare") < measurement.index("tooling/release_candidate.py acquire")
    assert "performance readme" not in measurement
    assert "--progress 'human' --progress-cadence '300'" in measurement
    captured = dry_run(ROOT, "performance-release", ["json", "60"]).stderr
    assert "--progress 'json' --progress-cadence '60'" in captured
    assert captured.count("tooling/release_candidate.py acquire") == 1
    publication = dry_run(ROOT, "performance-readme").stderr
    assert "performance readme tooling/performance-release.toml" in publication
    assert "performance measure" not in publication
    assert "performance prepare" not in publication
    assert "performance release" not in publication


@pytest.mark.parametrize("mode", ["fresh", "resume"])
@pytest.mark.parametrize("failed_phase", ["none", "prepare", "release"])
def test_release_recipe_stops_on_failure_without_publishing_readme(tmp_path: Path, failed_phase: str, mode: str) -> None:
    """Execute the real Bash recipe with a bounded CLI substitute, never fetching or timing."""
    root = tmp_path / "release with spaces"
    root.mkdir()
    directory = tmp_path / "fixture bin"
    directory.mkdir()
    launcher = directory / "uv"
    launcher.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'printf "%s\\n" "$*" >> "$PERFORMANCE_COMMAND_LOG"\n'
        'if [[ " $* " == *" performance $PERFORMANCE_FAIL_PHASE "* || ( "$PERFORMANCE_FAIL_PHASE" == release && " $* " == *" acquire "* ) ]]; then\n'
        '    printf "fixture %s failure\\n" "$PERFORMANCE_FAIL_PHASE" >&2\n'
        "    exit 7\n"
        "fi\n",
        encoding="utf-8",
        newline="\n",
    )
    launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR)
    log = root / "commands.log"
    result = run_command(
        "just",
        ["--no-deps", "--justfile", str(ROOT / "justfile"), "--working-directory", str(root), "performance-release", "human", "300", "target/checkpoint", mode],
        cwd=root,
        env=os.environ
        | {
            "PATH": f"{directory}{os.pathsep}{os.environ['PATH']}",
            "PERFORMANCE_COMMAND_LOG": str(log),
            "PERFORMANCE_FAIL_PHASE": failed_phase,
        },
        check=False,
    )
    commands = [shlex.split(line) for line in log.read_text(encoding="utf-8").splitlines()]
    phases = [
        command[command.index("performance") + 1 :] if "performance" in command else command[command.index("tooling/release_candidate.py") + 1 :]
        for command in commands
    ]
    expected = [["prepare", "tooling/performance-release.toml", "--fetch-tags", "--allow-git-mutations", "--release-intent"]] if mode == "fresh" else []
    if failed_phase != "prepare" or mode == "resume":
        expected.append(
            [
                "acquire",
                "--progress",
                "human",
                "--progress-cadence",
                "300",
                "--checkpoint",
                "target/checkpoint",
            ]
            + (["--resume"] if mode == "resume" else [])
        )
    assert phases == expected
    status = 0 if failed_phase == "none" or (mode == "resume" and failed_phase == "prepare") else 7
    assert result.returncode == status, result.stderr
    if status:
        assert f"fixture {failed_phase} failure" in result.stderr
    assert f"exit {status}" in (root / "target/bench-reports/full/runtime.txt").read_text(encoding="utf-8")
