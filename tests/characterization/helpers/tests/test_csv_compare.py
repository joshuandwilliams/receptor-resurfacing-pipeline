from __future__ import annotations

from pathlib import Path

import pytest

from tests.characterization.helpers.csv_compare import (
    compare_csv_exact,
    compare_csv_exact_modulo_paths,
    compare_csv_struct,
)
from tests.characterization.helpers.path_normalize import canonicalize_workdir_paths


def _write(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


@pytest.mark.local_unit
def test_compare_csv_exact_identity_passes(tmp_path):
    """Identical CSVs return passed=True."""
    ref = _write(tmp_path / "r.csv", "a,b,c\n1,2.0,x\n3,4.5,y\n")
    act = _write(tmp_path / "a.csv", "a,b,c\n1,2.0,x\n3,4.5,y\n")
    result = compare_csv_exact(ref, act)
    assert result.passed
    result.assert_passed()


@pytest.mark.local_unit
def test_compare_csv_exact_drift_within_tolerance_passes(tmp_path):
    """Numeric drift within abs_tol passes."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,1.0\n2,2.0\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,1.0000000001\n2,2.0\n")
    assert compare_csv_exact(ref, act).passed


@pytest.mark.local_unit
def test_compare_csv_exact_drift_outside_tolerance_fails(tmp_path):
    """Drift exceeding tolerance fails and reports row/column."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,1.0\n2,2.0\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,1.5\n2,2.0\n")
    result = compare_csv_exact(ref, act)
    assert not result.passed
    assert any("row 0" in d and "'b'" in d for d in result.differences)


@pytest.mark.local_unit
def test_compare_csv_exact_missing_column_header_fails(tmp_path):
    """A missing column is reported as a header mismatch, not a cell mismatch."""
    ref = _write(tmp_path / "r.csv", "a,b,c\n1,2,3\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,2\n")
    result = compare_csv_exact(ref, act)
    assert not result.passed
    assert "header" in result.message.lower()


@pytest.mark.local_unit
def test_compare_csv_exact_reordered_header_fails(tmp_path):
    """Same column names in different order is a header mismatch (CSV-EXACT cares about order)."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,2\n")
    act = _write(tmp_path / "a.csv", "b,a\n2,1\n")
    result = compare_csv_exact(ref, act)
    assert not result.passed
    assert "header" in result.message.lower()


@pytest.mark.local_unit
def test_compare_csv_exact_row_count_mismatch_fails(tmp_path):
    """Different row counts fail with a row-count message before any cell check runs."""
    ref = _write(tmp_path / "r.csv", "a\n1\n2\n")
    act = _write(tmp_path / "a.csv", "a\n1\n2\n3\n")
    result = compare_csv_exact(ref, act)
    assert not result.passed
    assert "row count" in result.message.lower()


@pytest.mark.local_unit
def test_compare_csv_exact_row_order_significant(tmp_path):
    """Same rows in different order is a mismatch under CSV-EXACT."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,x\n2,y\n")
    act = _write(tmp_path / "a.csv", "a,b\n2,y\n1,x\n")
    assert not compare_csv_exact(ref, act).passed


@pytest.mark.local_unit
def test_compare_csv_exact_nan_equals_nan(tmp_path):
    """NaN in both reference and actual at the same cell is treated as equal."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,\n2,\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,\n2,\n")
    assert compare_csv_exact(ref, act).passed


@pytest.mark.local_unit
def test_compare_csv_exact_nan_versus_value_fails(tmp_path):
    """NaN in reference vs a finite value in actual is a mismatch."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,3.0\n")
    assert not compare_csv_exact(ref, act).passed


@pytest.mark.local_unit
def test_compare_csv_exact_missing_reference_file_fails(tmp_path):
    """Nonexistent reference path fails cleanly with a clear message."""
    ref = tmp_path / "nope.csv"
    act = _write(tmp_path / "a.csv", "a\n1\n")
    result = compare_csv_exact(ref, act)
    assert not result.passed
    assert "reference" in result.message.lower()


@pytest.mark.local_unit
def test_compare_csv_exact_missing_actual_file_fails(tmp_path):
    """Nonexistent actual path fails cleanly with a clear message."""
    ref = _write(tmp_path / "r.csv", "a\n1\n")
    act = tmp_path / "nope.csv"
    result = compare_csv_exact(ref, act)
    assert not result.passed
    assert "actual" in result.message.lower()


@pytest.mark.local_unit
def test_compare_csv_exact_empty_body_passes(tmp_path):
    """Header-only CSVs with matching headers and zero rows are equal."""
    ref = _write(tmp_path / "r.csv", "a,b,c\n")
    act = _write(tmp_path / "a.csv", "a,b,c\n")
    assert compare_csv_exact(ref, act).passed


@pytest.mark.local_unit
def test_compare_csv_exact_assert_passed_raises_on_failure(tmp_path):
    """assert_passed() on a failed result raises AssertionError."""
    ref = _write(tmp_path / "r.csv", "a\n1\n")
    act = _write(tmp_path / "a.csv", "a\n2\n")
    result = compare_csv_exact(ref, act)
    with pytest.raises(AssertionError):
        result.assert_passed()


@pytest.mark.local_unit
def test_compare_csv_struct_tolerates_row_reordering(tmp_path):
    """CSV-STRUCT sorts both frames by sort_by before comparing — pure row reordering passes."""
    ref = _write(tmp_path / "r.csv", "id,v\n1,10\n2,20\n3,30\n")
    act = _write(tmp_path / "a.csv", "id,v\n3,30\n1,10\n2,20\n")
    assert compare_csv_struct(ref, act, sort_by=["id"]).passed


@pytest.mark.local_unit
def test_compare_csv_struct_still_detects_value_drift(tmp_path):
    """CSV-STRUCT still fails on cell-value mismatches after sorting by sort_by."""
    ref = _write(tmp_path / "r.csv", "id,v\n1,10\n2,20\n")
    act = _write(tmp_path / "a.csv", "id,v\n2,20\n1,99\n")
    result = compare_csv_struct(ref, act, sort_by=["id"])
    assert not result.passed


@pytest.mark.local_unit
def test_compare_csv_struct_invalid_sort_key_fails(tmp_path):
    """sort_by referencing a column not in the CSV fails with a clear message."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,2\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,2\n")
    result = compare_csv_struct(ref, act, sort_by=["does_not_exist"])
    assert not result.passed
    assert "sort_by" in result.message.lower()


# ---------------------------------------------------------------------------
# compare_csv_exact_modulo_paths
# ---------------------------------------------------------------------------

_REF_WORKDIR = "/hpc-home/alice/proj/work/ab/0123456789abcdef0123456789abcd"
_ACT_WORKDIR = "/hpc-home/bob/proj/work/cd/fedcba9876543210fedcba9876543"


@pytest.mark.local_unit
def test_compare_csv_exact_modulo_paths_identity_passes(tmp_path):
    """Two CSVs with no path content compare equal under modulo-paths."""
    ref = _write(tmp_path / "r.csv", "a,b\n1,2.0\n3,4.5\n")
    act = _write(tmp_path / "a.csv", "a,b\n1,2.0\n3,4.5\n")
    result = compare_csv_exact_modulo_paths(ref, act, canonicalize_workdir_paths)
    assert result.passed
    result.assert_passed()


@pytest.mark.local_unit
def test_compare_csv_exact_modulo_paths_normalises_workdir_hash(tmp_path):
    """Different workdir hashes in path cells normalise to the same value."""
    ref = _write(
        tmp_path / "r.csv",
        f"name,pdb_path,score\nseq1,{_REF_WORKDIR}/out.pdb,0.42\n",
    )
    act = _write(
        tmp_path / "a.csv",
        f"name,pdb_path,score\nseq1,{_ACT_WORKDIR}/out.pdb,0.42\n",
    )
    result = compare_csv_exact_modulo_paths(ref, act, canonicalize_workdir_paths)
    assert result.passed, result.differences


@pytest.mark.local_unit
def test_compare_csv_exact_modulo_paths_detects_value_drift(tmp_path):
    """Path cells normalised, but a numeric drift outside tolerance still fails."""
    ref = _write(
        tmp_path / "r.csv",
        f"name,pdb_path,score\nseq1,{_REF_WORKDIR}/out.pdb,0.42\n",
    )
    act = _write(
        tmp_path / "a.csv",
        f"name,pdb_path,score\nseq1,{_ACT_WORKDIR}/out.pdb,0.99\n",
    )
    result = compare_csv_exact_modulo_paths(ref, act, canonicalize_workdir_paths)
    assert not result.passed
    assert any("score" in d for d in result.differences)


@pytest.mark.local_unit
def test_compare_csv_exact_modulo_paths_normalises_hpc_home_prefix(tmp_path):
    """User-home prefixes (/hpc-home/<user>) collapse under the normalizer."""
    ref = _write(
        tmp_path / "r.csv",
        "name,home_path\nseq1,/hpc-home/alice/projects/run1.csv\n",
    )
    act = _write(
        tmp_path / "a.csv",
        "name,home_path\nseq1,/hpc-home/bob/projects/run1.csv\n",
    )
    result = compare_csv_exact_modulo_paths(ref, act, canonicalize_workdir_paths)
    assert result.passed, result.differences


@pytest.mark.local_unit
def test_compare_csv_exact_modulo_paths_missing_reference_fails(tmp_path):
    """Nonexistent reference path fails cleanly with a clear message."""
    ref = tmp_path / "nope.csv"
    act = _write(tmp_path / "a.csv", "a\n1\n")
    result = compare_csv_exact_modulo_paths(ref, act, canonicalize_workdir_paths)
    assert not result.passed
    assert "reference" in result.message.lower()


@pytest.mark.local_unit
def test_compare_csv_exact_modulo_paths_missing_actual_fails(tmp_path):
    """Nonexistent actual path fails cleanly with a clear message."""
    ref = _write(tmp_path / "r.csv", "a\n1\n")
    act = tmp_path / "nope.csv"
    result = compare_csv_exact_modulo_paths(ref, act, canonicalize_workdir_paths)
    assert not result.passed
    assert "actual" in result.message.lower()
