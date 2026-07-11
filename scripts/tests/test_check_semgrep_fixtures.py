"""Tests for Semgrep fixture-annotation validation."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import check_semgrep_fixtures

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def test_semgrep_results_parses_valid_result_objects(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(
        "SEMGREP_JSON",
        json.dumps({"results": [{"check_id": "rust.foo"}, {"check_id": "rust.bar"}]}),
    )

    results = check_semgrep_fixtures._semgrep_results()

    assert results is not None
    assert [result["check_id"] for result in results.results] == ["rust.foo", "rust.bar"]


def test_semgrep_results_rejects_malformed_result_objects(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    monkeypatch.setenv("SEMGREP_JSON", json.dumps({"results": [{"check_id": "rust.foo"}, "bad"]}))

    results = check_semgrep_fixtures._semgrep_results()

    assert results is None
    assert "result 1 is not an object" in capsys.readouterr().err


def test_main_accepts_matching_annotations(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    fixture = tmp_path / "fixture.rs"
    fixture.write_text(
        "// ruleid: rust.foo, rust.bar\n// ruleid: rust.foo\n",
        encoding="utf-8",
    )
    monkeypatch.setenv(
        "SEMGREP_JSON",
        json.dumps({"results": [{"check_id": "rust.foo"}, {"check_id": "rust.bar"}, {"check_id": "rust.foo"}]}),
    )
    monkeypatch.setattr(check_semgrep_fixtures.sys, "argv", ["check_semgrep_fixtures.py", str(fixture)])

    rc = check_semgrep_fixtures.main()

    assert rc == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_main_reports_missing_check_id(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    fixture = tmp_path / "fixture.rs"
    fixture.write_text("// ruleid: rust.foo\n", encoding="utf-8")
    monkeypatch.setenv("SEMGREP_JSON", json.dumps({"results": [{}]}))
    monkeypatch.setattr(check_semgrep_fixtures.sys, "argv", ["check_semgrep_fixtures.py", str(fixture)])

    rc = check_semgrep_fixtures.main()

    assert rc == 1
    assert "missing string field 'check_id'" in capsys.readouterr().err
