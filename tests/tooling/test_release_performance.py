"""Release preparation saves assertions; README retries never repeat measurement."""

import subprocess
from datetime import UTC, datetime
from pathlib import Path

import pytest
from research_repo_tools.just_inspect import dry_run
from research_repo_tools.process import run_git_command
from research_repo_tools.releases import PublishedRelease

from tooling import release_performance as adapter

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def consumer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Give preparation real local Git commits and the consumer's declared policy."""
    root = tmp_path / "consumer with spaces"
    root.mkdir()
    (root / "Cargo.toml").write_text('[package]\nname = "fixture"\nversion = "0.5.2"\n', encoding="utf-8")
    config = root / "tooling/measurements/full.toml"
    config.parent.mkdir(parents=True)
    config.write_bytes((ROOT / "tooling/measurements/full.toml").read_bytes())
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(tmp_path / "absent-gitconfig"))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    run_git_command(["--no-pager", "init", "--quiet"], cwd=root)
    for tag in ("v0.5.0", "v0.5.1", "head"):
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


def releases(monkeypatch: pytest.MonkeyPatch, *tags: str) -> None:
    """Keep provider discovery offline while exercising the real shared selector."""
    history = tuple(PublishedRelease(tag, datetime(2026, 10, index, tzinfo=UTC)) for index, tag in enumerate(tags, start=1))
    monkeypatch.setattr(adapter, "published_releases", lambda _root, **_kwargs: history)


@pytest.mark.parametrize("published", [False, True])
def test_preparation_resolves_prospective_head_or_published_tag(
    consumer: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    published: bool,
) -> None:
    current = "v0.5.1" if published else "v0.5.2"
    (consumer / "Cargo.toml").write_text(f'[package]\nname = "fixture"\nversion = "{current[1:]}"\n', encoding="utf-8")
    releases(monkeypatch, *(("v0.5.0", "v0.5.1") if published else ("v0.5.1", "v0.5.0")))
    monkeypatch.setattr(adapter, "tooling_main", lambda _argv: pytest.fail("preparation must not measure or publish"))
    before = run_git_command(["--no-pager", "show-ref"], cwd=consumer).stdout
    assert adapter.main(["capture"]) == 0
    output = capsys.readouterr().out
    assert "target/bench-reports/release-inputs.json (relative to the repository root; ignored by Git)" in output
    assert "No shell variables or exports are needed" in output
    assert "This command does not set PREVIOUS_TAG, CURRENT_COMMIT, or PREVIOUS_COMMIT in your shell" in output
    assert "saved file survives opening a new terminal" in output
    assert "After measurements finish, run: just performance-readme (no arguments)" in output
    assert "performance-readme reads the release tags and commits from this saved file" in output
    expected = [
        current,
        run_git_command(["--no-pager", "rev-parse", "v0.5.1" if published else "HEAD"], cwd=consumer).stdout.strip(),
        "v0.5.0",
        run_git_command(["--no-pager", "rev-parse", "v0.5.0"], cwd=consumer).stdout.strip(),
    ]
    assert adapter.saved_assertions(consumer) == expected
    assert run_git_command(["--no-pager", "show-ref"], cwd=consumer).stdout == before


def test_explicit_pair_skips_provider_and_failed_capture_preserves_record(consumer: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(adapter, "published_releases", lambda _root, **_kwargs: pytest.fail("explicit pairs need no provider"))
    assert adapter.main(["capture", "v0.5.1", "v0.5.0"]) == 0
    original = (consumer / adapter.INPUTS).read_bytes()
    assert adapter.main(["capture", "v0.5.2", "v0.5.0"]) == 1  # Missing current tag.
    assert (consumer / adapter.INPUTS).read_bytes() == original


@pytest.mark.parametrize("failure", ["empty-history", "provider", "missing-tag", "partial-pair", "old-release"])
def test_failed_preparation_never_writes_partial_assertions(
    consumer: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    failure: str,
) -> None:
    releases(monkeypatch, *(() if failure == "empty-history" else ("v0.4.6" if failure == "old-release" else "v0.5.1",)))
    if failure == "missing-tag":
        run_git_command(["--no-pager", "tag", "-d", "v0.5.1"], cwd=consumer)
    if failure == "provider":

        def fail_provider(root: Path, repository: str) -> tuple[PublishedRelease, ...]:
            raise subprocess.CalledProcessError(1, ["gh", "release", "list"], stderr="provider unavailable")

        monkeypatch.setattr(adapter, "published_releases", fail_provider)
    args = ["capture", "v0.5.2"] if failure == "partial-pair" else ["capture"]
    assert adapter.main(args) == 1
    assert not (consumer / adapter.INPUTS).exists()
    output = capsys.readouterr()
    assert "Release performance rejected:" in output.err
    assert "performance-prepare" in output.err
    expected = {
        "empty-history": "a published stable release is required",
        "provider": "provider unavailable",
        "missing-tag": "refs/tags/v0.5.1^{commit}",
        "partial-pair": "current and baseline tags must be provided together",
        "old-release": "requires v0.5.0 or newer",
    }
    assert expected[failure] in output.err


def test_readme_retry_uses_saved_assertions_and_preserves_measurements(consumer: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    assert adapter.main(["capture", "v0.5.1", "v0.5.0"]) == 0
    expected = adapter.saved_assertions(consumer)
    payload = consumer / "target/bench-reports/full/run.json"
    payload.parent.mkdir()
    payload.write_bytes(b"retained measurements must not be touched")
    before = payload.read_bytes()
    calls: list[list[str]] = []

    def publication(argv: list[str]) -> int:
        calls.append(argv)
        return 7 if len(calls) == 1 else 0

    monkeypatch.setattr(adapter, "tooling_main", publication)
    monkeypatch.setattr(adapter, "published_releases", lambda _root, **_kwargs: pytest.fail("publication must not rediscover inputs"))
    for name in ("CURRENT_COMMIT", "PREVIOUS_TAG", "PREVIOUS_COMMIT"):
        monkeypatch.delenv(name, raising=False)
    assert adapter.main(["publish"]) == 7
    assert adapter.main(["publish"]) == 0
    assert (
        calls
        == [
            [
                "performance",
                "publish",
                "tooling/performance-readme.toml",
                "--current-release",
                expected[0],
                "--current-revision",
                expected[1],
                "--baseline-release",
                expected[2],
                "--baseline-revision",
                expected[3],
            ]
        ]
        * 2
    )
    assert payload.read_bytes() == before


@pytest.mark.parametrize("record", [None, "{", '{"schema": 2, "assertions": []}', '{"schema": 1, "assertions": ["v0.5.1", "", "v0.5.0", ""]}'])
def test_missing_or_invalid_preparation_stops_before_publication(consumer: Path, monkeypatch: pytest.MonkeyPatch, record: str | None) -> None:
    if record is not None:
        path = consumer / adapter.INPUTS
        path.parent.mkdir(parents=True)
        path.write_text(record, encoding="utf-8")
    monkeypatch.setattr(adapter, "tooling_main", lambda _argv: pytest.fail("invalid preparation cannot publish"))
    assert adapter.main(["publish"]) == 1
    assert adapter.main(["publish", "v0.5.1"]) == 1


def test_recipes_keep_preparation_measurement_and_publication_separate() -> None:
    preparation = dry_run(ROOT, "performance-prepare").stderr
    assert preparation.index("git fetch origin --tags") < preparation.index("tooling.release_performance capture")
    measurement = dry_run(ROOT, "performance-release").stderr
    assert measurement.count("performance measure") == measurement.count("performance promote") == 1
    assert "performance publish" not in measurement
    assert "release_performance" not in measurement
    publication = dry_run(ROOT, "performance-readme").stderr
    assert "tooling.release_performance publish" in publication
    assert "performance measure" not in publication
