from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.characterization.helpers.json_compare import (
    compare_json_deep,
    compare_json_modulo_paths,
)
from tests.characterization.helpers.path_normalize import canonicalize_workdir_paths


def _write(path: Path, data) -> Path:
    path.write_text(json.dumps(data))
    return path


@pytest.mark.local_unit
def test_compare_json_deep_identity_passes(tmp_path):
    """Identical JSON structures return passed=True."""
    data = {"a": 1, "b": [1, 2, 3], "c": {"d": "e"}}
    ref = _write(tmp_path / "r.json", data)
    act = _write(tmp_path / "a.json", data)
    assert compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_deep_drift_within_tolerance_passes(tmp_path):
    """Numeric drift within abs_tol passes."""
    ref = _write(tmp_path / "r.json", {"v": 1.0})
    act = _write(tmp_path / "a.json", {"v": 1.0 + 1e-9})
    assert compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_deep_drift_outside_tolerance_fails(tmp_path):
    """Drift exceeding tolerance is reported."""
    ref = _write(tmp_path / "r.json", {"v": 1.0})
    act = _write(tmp_path / "a.json", {"v": 2.0})
    assert not compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_deep_nested_path_reported(tmp_path):
    """Mismatch at depth reports the dotted/indexed path (e.g., designs[2].y)."""
    ref = _write(tmp_path / "r.json", {"designs": [{"x": 1}, {"x": 2}, {"x": 3, "y": 4.0}]})
    act = _write(tmp_path / "a.json", {"designs": [{"x": 1}, {"x": 2}, {"x": 3, "y": 9.0}]})
    result = compare_json_deep(ref, act)
    assert not result.passed
    assert any("designs[2].y" in d for d in result.differences)


@pytest.mark.local_unit
def test_compare_json_deep_missing_key_fails(tmp_path):
    """A reference key absent from actual is reported as missing."""
    ref = _write(tmp_path / "r.json", {"a": 1, "b": 2})
    act = _write(tmp_path / "a.json", {"a": 1})
    result = compare_json_deep(ref, act)
    assert not result.passed
    assert any("b" in d and "missing" in d.lower() for d in result.differences)


@pytest.mark.local_unit
def test_compare_json_deep_unexpected_key_fails(tmp_path):
    """An actual key not present in reference is reported as unexpected."""
    ref = _write(tmp_path / "r.json", {"a": 1})
    act = _write(tmp_path / "a.json", {"a": 1, "b": 2})
    result = compare_json_deep(ref, act)
    assert not result.passed
    assert any("b" in d and "unexpected" in d.lower() for d in result.differences)


@pytest.mark.local_unit
def test_compare_json_deep_list_length_mismatch_fails(tmp_path):
    """Lists of different lengths fail."""
    ref = _write(tmp_path / "r.json", {"xs": [1, 2, 3]})
    act = _write(tmp_path / "a.json", {"xs": [1, 2, 3, 4]})
    assert not compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_deep_list_order_significant_by_default(tmp_path):
    """Same elements in different order fail by default."""
    ref = _write(tmp_path / "r.json", {"xs": [1, 2, 3]})
    act = _write(tmp_path / "a.json", {"xs": [3, 2, 1]})
    assert not compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_deep_list_order_insensitive_when_disabled(tmp_path):
    """list_orders_significant=False treats lists as unordered."""
    ref = _write(tmp_path / "r.json", {"xs": [1, 2, 3]})
    act = _write(tmp_path / "a.json", {"xs": [3, 2, 1]})
    assert compare_json_deep(ref, act, list_orders_significant=False).passed


@pytest.mark.local_unit
def test_compare_json_deep_type_mismatch_fails(tmp_path):
    """Different types at the same path (string vs int) fail with a type-mismatch message."""
    ref = _write(tmp_path / "r.json", {"x": "1"})
    act = _write(tmp_path / "a.json", {"x": 1})
    result = compare_json_deep(ref, act)
    assert not result.passed
    assert any("type mismatch" in d.lower() for d in result.differences)


@pytest.mark.local_unit
def test_compare_json_deep_bool_distinct_from_int(tmp_path):
    """JSON true must NOT match JSON 1 (bool is a distinct type)."""
    ref = _write(tmp_path / "r.json", {"x": True})
    act = _write(tmp_path / "a.json", {"x": 1})
    result = compare_json_deep(ref, act)
    assert not result.passed
    assert any("type mismatch" in d.lower() for d in result.differences)


@pytest.mark.local_unit
def test_compare_json_deep_null_handling(tmp_path):
    """null == null passes."""
    ref = _write(tmp_path / "r.json", {"x": None})
    act = _write(tmp_path / "a.json", {"x": None})
    assert compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_deep_null_vs_value_fails(tmp_path):
    """null vs a numeric value fails."""
    ref = _write(tmp_path / "r.json", {"x": None})
    act = _write(tmp_path / "a.json", {"x": 0})
    assert not compare_json_deep(ref, act).passed


@pytest.mark.local_unit
def test_compare_json_modulo_paths_workdir_hash_normalised(tmp_path):
    """Workdir-hash-only differences are erased by the path normalizer."""
    ref = _write(tmp_path / "r.json", {"path": "/work/ab/0123456789abcdef0123456789abcd/sub.txt"})
    act = _write(tmp_path / "a.json", {"path": "/work/cd/aaaaaaaaaaaaaaaaaaaaaaaaaaaa11/sub.txt"})
    assert compare_json_modulo_paths(ref, act, canonicalize_workdir_paths).passed


@pytest.mark.local_unit
def test_compare_json_modulo_paths_real_data_diff_still_caught(tmp_path):
    """Genuine non-path differences are still detected after normalisation."""
    ref = _write(tmp_path / "r.json", {
        "path": "/work/ab/0123456789abcdef0123456789abcd/x", "v": 1,
    })
    act = _write(tmp_path / "a.json", {
        "path": "/work/cd/aaaaaaaaaaaaaaaaaaaaaaaaaaaa11/x", "v": 99,
    })
    assert not compare_json_modulo_paths(ref, act, canonicalize_workdir_paths).passed


@pytest.mark.local_unit
def test_compare_json_deep_missing_reference_file_fails(tmp_path):
    """Nonexistent reference path fails cleanly."""
    act = _write(tmp_path / "a.json", {"x": 1})
    result = compare_json_deep(tmp_path / "nope.json", act)
    assert not result.passed
    assert "reference" in result.message.lower()


@pytest.mark.local_unit
def test_compare_json_deep_missing_actual_file_fails(tmp_path):
    """Nonexistent actual path fails cleanly."""
    ref = _write(tmp_path / "r.json", {"x": 1})
    result = compare_json_deep(ref, tmp_path / "nope.json")
    assert not result.passed
    assert "actual" in result.message.lower()


@pytest.mark.local_unit
def test_compare_json_modulo_paths_missing_reference_file_fails(tmp_path):
    """Nonexistent reference path fails cleanly under JSON-MODULO-PATHS too."""
    act = _write(tmp_path / "a.json", {"x": 1})
    result = compare_json_modulo_paths(tmp_path / "nope.json", act, canonicalize_workdir_paths)
    assert not result.passed
    assert "reference" in result.message.lower()
