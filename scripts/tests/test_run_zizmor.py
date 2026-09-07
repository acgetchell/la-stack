"""Exercise authentication selection without network access or real credentials."""

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOKEN_ENV = ("ZIZMOR_GITHUB_TOKEN", "GH_TOKEN", "GITHUB_TOKEN")


def _shim(path: Path, body: str) -> None:
    path.write_text("#!/bin/bash\nset -eu\n" + body, encoding="utf-8", newline="\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _run(
    tmp_path: Path,
    env_overrides: dict[str, str],
    *,
    expected: str = "",
    gh_status: int | None = 0,
    scanner_status: int = 0,
) -> subprocess.CompletedProcess[str]:
    bash = shutil.which("bash")
    assert bash is not None, "Bash is required by the repository's Just recipes"
    _shim(
        tmp_path / "zizmor",
        '[[ "${GH_TOKEN:-}" == "$EXPECTED_TOKEN" ]] || exit 91\n'
        '[[ "${ZIZMOR_GITHUB_TOKEN:-}" == "$EXPECTED_TOKEN" ]] || exit 92\n'
        'printf "scanner args: %s\\n" "$*"\n'
        'exit "$SCANNER_STATUS"\n',
    )
    if gh_status is not None:
        _shim(
            tmp_path / "gh",
            '[[ "$*" == "auth token" ]] || exit 93\n'
            'echo "called" >> "$GH_CALL_LOG"\n'
            'echo "test-gh-secret"\n'
            'echo "test-gh-diagnostic-secret" >&2\n'
            'exit "$GH_STATUS"\n',
        )
    env = {
        key: value
        for key, value in os.environ.items()
        if key not in (*TOKEN_ENV, "ZIZMOR_OFFLINE", "ZIZMOR_NO_ONLINE_AUDITS", "BASH_ENV", "SHELLOPTS", "BASHOPTS")
    }
    env.update(
        EXPECTED_TOKEN=expected,
        GH_STATUS=str(gh_status),
        GH_CALL_LOG=(tmp_path / "gh-calls").as_posix(),
        SCANNER_STATUS=str(scanner_status),
    )
    env.update(env_overrides)
    # Bash's PWD supplies a POSIX path on Windows too, without drive-letter colons
    # being interpreted as PATH separators or reintroducing the real gh binary.
    result = subprocess.run(  # noqa: S603 - resolved Bash, fixed script, test-owned environment.
        [
            bash,
            "--noprofile",
            "--norc",
            "-c",
            'cd "$1"; export PATH="$PWD"; cd "$2"; source scripts/run_zizmor.sh',
            "test-zizmor",
            tmp_path.as_posix(),
            REPO_ROOT.as_posix(),
        ],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        encoding="utf-8",
        check=False,
        timeout=10,
    )
    for secret in (*(env_overrides.get(key, "") for key in TOKEN_ENV), "test-gh-secret", "test-gh-diagnostic-secret"):
        if secret:
            assert secret not in result.stdout + result.stderr
    return result


@pytest.mark.parametrize(
    ("tokens", "expected"),
    [
        ({"ZIZMOR_GITHUB_TOKEN": "test-zizmor-secret", "GH_TOKEN": "test-gh-env-secret", "GITHUB_TOKEN": "test-github-secret"}, "test-zizmor-secret"),
        ({"ZIZMOR_GITHUB_TOKEN": "", "GH_TOKEN": "test-gh-env-secret", "GITHUB_TOKEN": "test-github-secret"}, "test-gh-env-secret"),
        ({"GITHUB_TOKEN": "test-github-secret"}, "test-github-secret"),
    ],
)
def test_environment_precedence_and_no_token_logging(tmp_path: Path, tokens: dict[str, str], expected: str) -> None:
    result = _run(tmp_path, tokens, expected=expected)

    assert result.returncode == 0, result.stderr
    assert "authenticated online audits" in result.stdout
    assert "scanner args: --persona regular .github" in result.stdout
    assert not (tmp_path / "gh-calls").exists()


def test_authenticated_gh_fallback(tmp_path: Path) -> None:
    result = _run(tmp_path, {}, expected="test-gh-secret")

    assert result.returncode == 0, result.stderr
    assert "authenticated online audits" in result.stdout
    assert "scanner args: --persona regular .github" in result.stdout
    assert (tmp_path / "gh-calls").read_text(encoding="utf-8") == "called\n"


@pytest.mark.parametrize("gh_status", [1, None])
def test_offline_fallback_without_authentication(tmp_path: Path, gh_status: int | None) -> None:
    result = _run(tmp_path, {}, gh_status=gh_status)

    assert result.returncode == 0, result.stderr
    assert "no GitHub token available; using offline audits" in result.stdout
    assert "scanner args: --offline --persona regular .github" in result.stdout


@pytest.mark.parametrize("offline_env", ["ZIZMOR_OFFLINE", "ZIZMOR_NO_ONLINE_AUDITS"])
def test_explicit_offline_skips_credential_lookup(tmp_path: Path, offline_env: str) -> None:
    result = _run(tmp_path, {offline_env: "true"})

    assert result.returncode == 0, result.stderr
    assert "offline audits requested" in result.stdout
    assert "scanner args: --offline --persona regular .github" in result.stdout
    assert not (tmp_path / "gh-calls").exists()


@pytest.mark.parametrize("gh_status", [0, 1])
def test_scanner_failure_propagates_without_offline_retry(tmp_path: Path, gh_status: int) -> None:
    result = _run(tmp_path, {}, expected="test-gh-secret" if gh_status == 0 else "", gh_status=gh_status, scanner_status=14)

    assert result.returncode == 14
    assert result.stdout.count("scanner args:") == 1
