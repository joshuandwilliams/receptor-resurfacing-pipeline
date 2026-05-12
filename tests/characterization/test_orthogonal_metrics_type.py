"""Unit tests for bin/orthogonal_metrics.py (Phase 4 Tier 4.3).

Distinct from the existing tests/characterization/test_orthogonal_metrics.py
(which characterises the existing orthogonal-metrics pipeline outputs).
This file tests the NEW deep-module type only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import orthogonal_metrics as om  # noqa: E402
import pipeline_thresholds as pt  # noqa: E402
import af3_confidence as af3  # noqa: E402
import protein_structure_prediction as psp  # noqa: E402


def _fake_psp(tmp_path):
    p = tmp_path / "f.pdb"; p.write_text("")
    return psp.ProteinStructurePrediction(
        path=p, receptor_chain="A", effector_chain="B")


def _fake_af3(best_ra_eff=2.0):
    a = af3.AF3ConfidenceAggregate(
        af3_output_dir=Path("/dev/null"), ground_truth=Path("/dev/null"))
    a._parsed = True
    if best_ra_eff is not None:
        a._ra_effs = [best_ra_eff]
    a._iptms = [0.5]
    return a


@pytest.mark.local_unit
class TestConstruction:
    def test_basic(self, tmp_path):
        m = om.OrthogonalMetrics(
            mpnn_sequence_id="design_03_seq_01",
            canonical_prediction=_fake_psp(tmp_path), af3=_fake_af3(),
            sc=0.65, bsa=700.0, ddg=-32.0)
        assert m.mpnn_sequence_id == "design_03_seq_01"

    def test_empty_id_raises(self, tmp_path):
        with pytest.raises(ValueError, match="mpnn_sequence_id"):
            om.OrthogonalMetrics(
                mpnn_sequence_id="",
                canonical_prediction=_fake_psp(tmp_path), af3=_fake_af3())

    def test_wrong_psp_type_raises(self):
        with pytest.raises(TypeError, match="canonical_prediction"):
            om.OrthogonalMetrics(
                mpnn_sequence_id="x",
                canonical_prediction="not a psp", af3=_fake_af3())

    def test_wrong_af3_type_raises(self, tmp_path):
        with pytest.raises(TypeError, match="af3"):
            om.OrthogonalMetrics(
                mpnn_sequence_id="x",
                canonical_prediction=_fake_psp(tmp_path), af3="not af3")


@pytest.mark.local_unit
class TestOrthogonalGate:
    """Sc + BSA + ΔΔG only — interface_plddt explicitly NOT in the gate."""

    def _make(self, tmp_path, sc=0.65, bsa=700.0, ddg=-32.0):
        return om.OrthogonalMetrics(
            mpnn_sequence_id="x",
            canonical_prediction=_fake_psp(tmp_path), af3=_fake_af3(),
            sc=sc, bsa=bsa, ddg=ddg)

    def test_all_pass(self, tmp_path):
        m = self._make(tmp_path)
        t = pt.PipelineInternalThresholds.default()
        assert m.passes_orthogonal_filters(t) is True
        assert m.failed_filter_names(t) == []

    def test_sc_too_low(self, tmp_path):
        m = self._make(tmp_path, sc=0.40)
        t = pt.PipelineInternalThresholds.default()
        assert m.passes_orthogonal_filters(t) is False
        assert "sc" in m.failed_filter_names(t)

    def test_bsa_too_low(self, tmp_path):
        m = self._make(tmp_path, bsa=200.0)
        t = pt.PipelineInternalThresholds.default()
        assert m.passes_orthogonal_filters(t) is False
        assert "bsa" in m.failed_filter_names(t)

    def test_ddg_too_high(self, tmp_path):
        m = self._make(tmp_path, ddg=10.0)
        t = pt.PipelineInternalThresholds.default()
        assert m.passes_orthogonal_filters(t) is False
        assert "ddg" in m.failed_filter_names(t)

    def test_missing_metric_fails_gate(self, tmp_path):
        m = self._make(tmp_path, sc=None)
        t = pt.PipelineInternalThresholds.default()
        assert m.passes_orthogonal_filters(t) is False
        assert "missing_sc" in m.failed_filter_names(t)

    def test_multiple_failures_all_reported(self, tmp_path):
        m = self._make(tmp_path, sc=0.40, bsa=200.0, ddg=10.0)
        t = pt.PipelineInternalThresholds.default()
        failed = m.failed_filter_names(t)
        assert "sc" in failed and "bsa" in failed and "ddg" in failed


@pytest.mark.local_unit
class TestAF3Informational:
    def test_disagreement_is_informational_only(self, tmp_path):
        m = om.OrthogonalMetrics(
            mpnn_sequence_id="x",
            canonical_prediction=_fake_psp(tmp_path),
            af3=_fake_af3(best_ra_eff=10.0),
            sc=0.65, bsa=700.0, ddg=-32.0)
        t = pt.PipelineInternalThresholds.default()
        assert m.af3_disagrees(t) is True
        assert m.passes_orthogonal_filters(t) is True


@pytest.mark.local_unit
class TestSummaryRow:
    def test_has_orthogonal_and_af3_keys(self, tmp_path):
        m = om.OrthogonalMetrics(
            mpnn_sequence_id="design_03_seq_01",
            canonical_prediction=_fake_psp(tmp_path),
            af3=_fake_af3(best_ra_eff=2.0),
            sc=0.65, bsa=700.0, hbonds=8, ddg=-32.0)
        row = m.to_summary_row()
        assert row["mpnn_sequence"] == "design_03_seq_01"
        assert row["sc"] == "0.6500"
        assert row["bsa"] == "700.00"
        assert row["hbonds"] == "8"
        assert row["rosetta_ddg"] == "-32.00"
        assert "af3_nomsa_best_ra_eff" in row

    def test_blank_for_missing(self, tmp_path):
        m = om.OrthogonalMetrics(
            mpnn_sequence_id="x",
            canonical_prediction=_fake_psp(tmp_path), af3=_fake_af3(),
            sc=None, bsa=None, ddg=None, hbonds=None)
        row = m.to_summary_row()
        assert row["sc"] == ""
        assert row["bsa"] == ""
        assert row["rosetta_ddg"] == ""
        assert row["hbonds"] == ""


@pytest.mark.local_unit
class TestImmutability:
    def test_frozen(self, tmp_path):
        m = om.OrthogonalMetrics(
            mpnn_sequence_id="x",
            canonical_prediction=_fake_psp(tmp_path), af3=_fake_af3())
        with pytest.raises((AttributeError, Exception)):
            m.sc = 0.99  # type: ignore
