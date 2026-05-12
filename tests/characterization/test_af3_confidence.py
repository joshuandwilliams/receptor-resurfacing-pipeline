"""Unit tests for bin/af3_confidence.py (Phase 4 Tier 0.4).

The AF3 parse logic requires gemmi + a real CIF + a real ground-truth
PDB, which is heavy for a local unit test.  Tier 0 tests cover:
  - aggregate-property arithmetic on hand-injected ra_effs / iptms
  - empty-input / no-data branches
  - to_summary_row schema
  - n_correct_interface threshold check
  - af3_disagrees flag

End-to-end-against-real-AF3-output is tested at hpc-tier (when AF3
runs as part of the orthogonal stage), not here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import af3_confidence as af3  # noqa: E402
import pipeline_thresholds as pt  # noqa: E402


def _make_aggregate_with(ra_effs, iptms, failures=None, cif_paths=None):
    """Build an AF3ConfidenceAggregate with pre-set state, bypassing the
    real-CIF parsing path so we can exercise the aggregate math in
    isolation."""
    a = af3.AF3ConfidenceAggregate(
        af3_output_dir=Path("/dev/null"),
        ground_truth=Path("/dev/null"),
    )
    a._parsed = True
    a._ra_effs = list(ra_effs)
    a._iptms = list(iptms)
    a._failures = list(failures or [])
    a._cif_paths = list(cif_paths or [])
    return a


@pytest.mark.local_unit
class TestAggregateMath:
    def test_best_ra_eff_is_min(self):
        a = _make_aggregate_with([5.0, 3.0, 7.0], [])
        assert a.best_ra_eff == 3.0

    def test_mean_ra_eff(self):
        a = _make_aggregate_with([2.0, 4.0, 6.0], [])
        assert a.mean_ra_eff == 4.0

    def test_best_iptm_is_max(self):
        a = _make_aggregate_with([], [0.5, 0.3, 0.8])
        assert a.best_iptm == 0.8

    def test_mean_iptm(self):
        a = _make_aggregate_with([], [0.2, 0.4, 0.6])
        assert a.mean_iptm == pytest.approx(0.4)

    def test_empty_ra_effs_yields_none(self):
        a = _make_aggregate_with([], [])
        assert a.best_ra_eff is None
        assert a.mean_ra_eff is None

    def test_empty_iptms_yields_none(self):
        a = _make_aggregate_with([3.0], [])
        assert a.best_iptm is None
        assert a.mean_iptm is None

    def test_total_predictions_is_count_of_usable_ra_effs(self):
        a = _make_aggregate_with([1.0, 2.0, 3.0], [0.5])
        assert a.total_predictions == 3


@pytest.mark.local_unit
class TestThresholdMethods:
    def test_n_correct_interface_counts_below_threshold(self):
        a = _make_aggregate_with([2.0, 4.0, 4.999, 5.0, 6.0], [])
        t = pt.PipelineInternalThresholds.default()  # af3_ra_max=5.0
        # Strict < 5.0 → values 2.0, 4.0, 4.999 count
        assert a.n_correct_interface(t) == 3

    def test_n_correct_interface_respects_threshold_override(self):
        a = _make_aggregate_with([2.0, 4.0, 6.0], [])
        t = pt.PipelineInternalThresholds.default().with_overrides(
            {"orthogonal_af3_ra_max": 3.0}
        )
        assert a.n_correct_interface(t) == 1

    def test_n_correct_interface_zero_when_no_data(self):
        a = _make_aggregate_with([], [])
        t = pt.PipelineInternalThresholds.default()
        assert a.n_correct_interface(t) == 0

    def test_af3_disagrees_true_when_best_is_above_threshold(self):
        a = _make_aggregate_with([6.0, 7.0], [])
        t = pt.PipelineInternalThresholds.default()  # af3_ra_max=5.0
        assert a.af3_disagrees(t) is True

    def test_af3_disagrees_false_when_best_is_below_threshold(self):
        a = _make_aggregate_with([4.5, 6.0], [])
        t = pt.PipelineInternalThresholds.default()
        assert a.af3_disagrees(t) is False

    def test_af3_disagrees_false_when_no_data(self):
        """No data → no signal, treated as 'no disagreement'."""
        a = _make_aggregate_with([], [])
        t = pt.PipelineInternalThresholds.default()
        assert a.af3_disagrees(t) is False


@pytest.mark.local_unit
class TestSummaryRow:
    def test_summary_row_schema_matches_existing_csv(self):
        a = _make_aggregate_with([3.5, 4.5], [0.6, 0.7])
        row = a.to_summary_row("test_seq")
        # Must match the schema from bin/parse_af3_output.py:189-198
        expected_keys = {
            "seq_name",
            "af3_nomsa_best_ra_eff",
            "af3_nomsa_mean_ra_eff",
            "af3_nomsa_best_iptm",
            "af3_nomsa_mean_iptm",
            "af3_nomsa_total_predictions",
            "af3_nomsa_failures",
        }
        assert set(row.keys()) == expected_keys
        assert row["seq_name"] == "test_seq"
        assert row["af3_nomsa_best_ra_eff"] == "3.500"  # 3 dp
        assert row["af3_nomsa_best_iptm"] == "0.7000"   # 4 dp
        assert row["af3_nomsa_total_predictions"] == "2"

    def test_summary_row_blank_for_missing(self):
        a = _make_aggregate_with([], [], failures=["no_cif_files_found"])
        row = a.to_summary_row("test_seq")
        assert row["af3_nomsa_best_ra_eff"] == ""
        assert row["af3_nomsa_mean_iptm"] == ""
        assert row["af3_nomsa_total_predictions"] == "0"
        assert "no_cif_files_found" in row["af3_nomsa_failures"]


@pytest.mark.local_unit
class TestFailureHandling:
    def test_failures_round_trip(self):
        a = _make_aggregate_with(
            [], [],
            failures=["no_cif_files_found", "per_prediction_rmsd_errors:2"],
        )
        f = a.failures
        assert "no_cif_files_found" in f
        assert "per_prediction_rmsd_errors:2" in f

    def test_failures_is_a_copy(self):
        """Mutating the returned list does not affect internal state."""
        a = _make_aggregate_with([], [], failures=["a", "b"])
        out = a.failures
        out.append("c")
        assert "c" not in a.failures


@pytest.mark.local_unit
class TestParsingFromEmptyDir:
    def test_missing_dir_emits_failure_tag(self, tmp_path):
        a = af3.AF3ConfidenceAggregate(
            af3_output_dir=tmp_path / "does_not_exist",
            ground_truth=tmp_path / "ref.pdb",
        )
        # Force parse
        assert a.total_predictions == 0
        assert any("af3_output_dir_missing" in f for f in a.failures)

    def test_dir_with_no_cifs_emits_failure_tag(self, tmp_path):
        # Create an empty directory; no CIFs inside.
        empty = tmp_path / "af3_out"
        empty.mkdir()
        a = af3.AF3ConfidenceAggregate(
            af3_output_dir=empty,
            ground_truth=tmp_path / "ref.pdb",
        )
        assert a.total_predictions == 0
        assert any("no_cif_files_found" in f for f in a.failures)
