"""Tests for Semgrep fixture-annotation validation."""

import json
from typing import TYPE_CHECKING

import check_semgrep_fixtures

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def _result(check_id: str, line: int, end_line: int | None = None) -> dict[str, object]:
    return {"check_id": check_id, "start": {"line": line}, "end": {"line": line if end_line is None else end_line}}


def test_semgrep_results_parses_valid_result_objects(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(
        "SEMGREP_JSON",
        json.dumps({"results": [_result("rust.foo", 2), _result("rust.bar", 4)]}),
    )

    results = check_semgrep_fixtures._semgrep_results()

    assert results is not None
    assert [result["check_id"] for result in results.results] == ["rust.foo", "rust.bar"]


def test_semgrep_results_rejects_malformed_result_objects(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    monkeypatch.setenv("SEMGREP_JSON", json.dumps({"results": [_result("rust.foo", 2), "bad"]}))

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
        json.dumps({"results": [_result("rust.foo", 2), _result("rust.bar", 2), _result("rust.foo", 3)]}),
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
    captured = capsys.readouterr()
    assert "missing string field 'check_id'" in captured.err
    assert "missing positive integer field 'start.line'" in captured.err
    assert "missing positive integer field 'end.line'" in captured.err


def test_main_rejects_reversed_result_span(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    fixture = tmp_path / "fixture.rs"
    fixture.write_text("// ruleid: rust.foo\nbad();\n", encoding="utf-8")
    monkeypatch.setenv("SEMGREP_JSON", json.dumps({"results": [_result("rust.foo", 4, 2)]}))
    monkeypatch.setattr(check_semgrep_fixtures.sys, "argv", ["check_semgrep_fixtures.py", str(fixture)])

    assert check_semgrep_fixtures.main() == 1
    assert "result 0 has end.line 2 before start.line 4" in capsys.readouterr().err


def test_main_rejects_findings_at_wrong_lines_even_when_rule_counts_match(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    fixture = tmp_path / "fixture.rs"
    fixture.write_text("// ruleid: rust.foo\nbad_one();\n// ruleid: rust.foo\nbad_two();\n", encoding="utf-8")
    monkeypatch.setenv(
        "SEMGREP_JSON",
        json.dumps({"results": [_result("rust.foo", 2), _result("rust.foo", 5)]}),
    )
    monkeypatch.setattr(check_semgrep_fixtures.sys, "argv", ["check_semgrep_fixtures.py", str(fixture)])

    rc = check_semgrep_fixtures.main()

    assert rc == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "rust.foo at line 4: expected finding not reported" in captured.err
    assert "rust.foo at lines 5: unexpected finding" in captured.err


def test_main_matches_overlapping_spans_by_earliest_end_line(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    fixture = tmp_path / "fixture.rs"
    fixture.write_text("// ruleid: rust.foo\nbad_one();\n// ruleid: rust.foo\nbad_two();\n", encoding="utf-8")
    monkeypatch.setenv(
        "SEMGREP_JSON",
        json.dumps({"results": [_result("rust.foo", 2, 4), _result("rust.foo", 2)]}),
    )
    monkeypatch.setattr(check_semgrep_fixtures.sys, "argv", ["check_semgrep_fixtures.py", str(fixture)])

    assert check_semgrep_fixtures.main() == 0


def test_main_matches_markdown_finding_after_blank_line_and_code_fence(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    fixture = tmp_path / "fixture.md"
    fixture.write_text("<!-- ruleid: docs.foo -->\n\n```bash\nbad-command\n```\n", encoding="utf-8")
    monkeypatch.setenv("SEMGREP_JSON", json.dumps({"results": [_result("docs.foo", 4)]}))
    monkeypatch.setattr(check_semgrep_fixtures.sys, "argv", ["check_semgrep_fixtures.py", str(fixture)])

    assert check_semgrep_fixtures.main() == 0
