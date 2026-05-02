from __future__ import annotations

from pathlib import Path

import pytest

from tests.characterization.helpers.text_compare import compare_text_exact


def _write(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


@pytest.mark.local_unit
def test_compare_text_exact_identity_passes(tmp_path):
    """Identical text files return passed=True."""
    ref = _write(tmp_path / "r.txt", "alpha\nbeta\ngamma\n")
    act = _write(tmp_path / "a.txt", "alpha\nbeta\ngamma\n")
    result = compare_text_exact(ref, act)
    assert result.passed
    result.assert_passed()


@pytest.mark.local_unit
def test_compare_text_exact_single_byte_difference_fails(tmp_path):
    """A one-character substitution fails and reports the offending line."""
    ref = _write(tmp_path / "r.txt", "alpha\nbeta\ngamma\n")
    act = _write(tmp_path / "a.txt", "alpha\nBeta\ngamma\n")
    result = compare_text_exact(ref, act)
    assert not result.passed
    assert any("line 2" in d for d in result.differences)


@pytest.mark.local_unit
def test_compare_text_exact_missing_reference_file_fails(tmp_path):
    """Nonexistent reference path fails cleanly with a clear message."""
    ref = tmp_path / "nope.txt"
    act = _write(tmp_path / "a.txt", "alpha\n")
    result = compare_text_exact(ref, act)
    assert not result.passed
    assert "reference" in result.message.lower()


@pytest.mark.local_unit
def test_compare_text_exact_missing_actual_file_fails(tmp_path):
    """Nonexistent actual path fails cleanly with a clear message."""
    ref = _write(tmp_path / "r.txt", "alpha\n")
    act = tmp_path / "nope.txt"
    result = compare_text_exact(ref, act)
    assert not result.passed
    assert "actual" in result.message.lower()


@pytest.mark.local_unit
def test_compare_text_exact_trailing_newline_difference_fails(tmp_path):
    """Files differing only in a trailing newline are NOT byte-equal — must fail."""
    ref = _write(tmp_path / "r.txt", "alpha\nbeta\n")
    act = _write(tmp_path / "a.txt", "alpha\nbeta")
    result = compare_text_exact(ref, act)
    assert not result.passed
