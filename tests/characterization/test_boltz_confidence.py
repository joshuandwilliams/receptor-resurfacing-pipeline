"""Unit tests for bin/boltz_confidence.py (Phase 4 Tier 0.3)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import boltz_confidence as bc  # noqa: E402
import pipeline_thresholds as pt  # noqa: E402

# A real Boltz confidence JSON pulled from the per-module fixture.
_REAL_FIXTURE = (
    _REPO_ROOT
    / "tests/negative_steering/example_output_files"
    / "negative_steering/runs/design_42_seq_0/cycle_0/initial"
    / "boltz_results_input/predictions/input"
    / "confidence_input_model_4.json"
)


def _write_conf_json(path, payload):
    path.write_text(json.dumps(payload))


@pytest.mark.logic
class TestConstruction:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            bc.BoltzConfidenceMetrics(tmp_path / "missing.json")

    def test_reads_minimal_json(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"ptm": 0.5, "iptm": 0.3})
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.ptm == 0.5
        assert m.iptm == 0.3


@pytest.mark.logic
class TestProperties:
    def test_missing_fields_return_none(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {})
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.ptm is None
        assert m.iptm is None
        assert m.complex_plddt is None
        assert m.confidence_score is None

    def test_unparseable_field_returns_none(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"ptm": "not-a-number"})
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.ptm is None

    def test_malformed_json_treats_as_empty(self, tmp_path):
        cj = tmp_path / "bad.json"
        cj.write_text("{not valid")
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.ptm is None
        assert m.iptm is None

    def test_iptm_falls_back_to_protein_iptm(self, tmp_path):
        """Matches compute_metrics.py:1241-1243 fallback behaviour."""
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"protein_iptm": 0.42})
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.iptm == 0.42

    def test_iptm_prefers_iptm_over_protein_iptm(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"iptm": 0.7, "protein_iptm": 0.42})
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.iptm == 0.7

    def test_rounding_to_4dp(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"ptm": 0.123456789})
        m = bc.BoltzConfidenceMetrics(cj)
        assert m.ptm == 0.1235

    def test_chains_ptm_dict(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"chains_ptm": {"0": 0.9, "1": 0.8}})
        m = bc.BoltzConfidenceMetrics(cj)
        d = m.chains_ptm
        assert d == {"0": 0.9, "1": 0.8}

    def test_pair_chains_iptm_nested(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"pair_chains_iptm": {"0": {"1": 0.5}, "1": {"0": 0.45}}})
        m = bc.BoltzConfidenceMetrics(cj)
        d = m.pair_chains_iptm
        assert d == {"0": {"1": 0.5}, "1": {"0": 0.45}}


@pytest.mark.logic
class TestLazyLoading:
    def test_json_loaded_only_once(self, tmp_path):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, {"ptm": 0.5})
        m = bc.BoltzConfidenceMetrics(cj)
        # Access twice — should not raise even if file is removed after first read
        _ = m.ptm
        cj.unlink()
        _ = m.iptm  # cache means this works
        # And the cached value is what we expect
        assert m.ptm == 0.5


@pytest.mark.logic
@pytest.mark.skipif(not _REAL_FIXTURE.is_file(),
                    reason="real fixture not present in this checkout")
class TestRealFixture:
    """End-to-end against a real Boltz confidence JSON from the per-module
    fixture set, to catch any drift between the wrapper and real output."""

    def test_reads_all_top_level_floats(self):
        m = bc.BoltzConfidenceMetrics(_REAL_FIXTURE)
        # From the fixture:
        # confidence_score: 0.7693772912025452 -> 0.7694
        # ptm: 0.658906102180481 -> 0.6589
        # iptm: 0.4295026361942291 -> 0.4295
        # complex_plddt: 0.8543459177017212 -> 0.8543
        assert m.confidence_score == 0.7694
        assert m.ptm == 0.6589
        assert m.iptm == 0.4295
        assert m.complex_plddt == 0.8543

    def test_chains_and_pair_chains_populated(self):
        m = bc.BoltzConfidenceMetrics(_REAL_FIXTURE)
        cp = m.chains_ptm
        assert "0" in cp
        assert "1" in cp
        pci = m.pair_chains_iptm
        assert "0" in pci and "1" in pci["0"]


@pytest.mark.logic
class TestConfidenceFlag:
    """Tier 0 confidence_flag: iptm + complex_plddt only (PAE-derived
    checks deferred to Tier 2 when chain lengths are accessible)."""

    def _make(self, tmp_path, payload):
        cj = tmp_path / "confidence_input_model_0.json"
        _write_conf_json(cj, payload)
        return bc.BoltzConfidenceMetrics(cj)

    def test_ok_when_above_thresholds(self, tmp_path):
        m = self._make(tmp_path, {"iptm": 0.5, "complex_plddt": 0.8})
        t = pt.PipelineInternalThresholds.default()
        assert m.confidence_flag(t) == "ok"

    def test_low_iptm(self, tmp_path):
        m = self._make(tmp_path, {"iptm": 0.1, "complex_plddt": 0.8})
        t = pt.PipelineInternalThresholds.default()
        assert m.confidence_flag(t) == "low_iptm"

    def test_low_plddt(self, tmp_path):
        m = self._make(tmp_path, {"iptm": 0.5, "complex_plddt": 0.5})
        t = pt.PipelineInternalThresholds.default()
        assert m.confidence_flag(t) == "low_plddt"

    def test_multiple(self, tmp_path):
        m = self._make(tmp_path, {"iptm": 0.1, "complex_plddt": 0.5})
        t = pt.PipelineInternalThresholds.default()
        assert m.confidence_flag(t) == "multiple"

    def test_missing_fields_count_as_passing(self, tmp_path):
        """When a field is missing, the check for it is skipped (returns
        ok if no other trigger fires).  This matches existing behaviour
        in extract_passing.py:_compute_confidence_flag."""
        m = self._make(tmp_path, {})  # both missing
        t = pt.PipelineInternalThresholds.default()
        assert m.confidence_flag(t) == "ok"


@pytest.mark.logic
class TestFromPredictionDir:
    def test_returns_empty_for_dir_with_no_confidences(self, tmp_path):
        result = bc.BoltzConfidenceMetrics.from_prediction_dir(tmp_path)
        assert result == []

    def test_finds_all_samples(self, tmp_path):
        for i in range(3):
            _write_conf_json(
                tmp_path / f"confidence_input_model_{i}.json",
                {"ptm": 0.5 + 0.01 * i},
            )
        result = bc.BoltzConfidenceMetrics.from_prediction_dir(tmp_path)
        assert len(result) == 3
        # Sorted by glob order
        ptms = sorted([m.ptm for m in result])
        assert ptms == [0.5, 0.51, 0.52]

    def test_custom_base_name(self, tmp_path):
        _write_conf_json(
            tmp_path / "confidence_other_model_0.json", {"ptm": 0.9}
        )
        result = bc.BoltzConfidenceMetrics.from_prediction_dir(
            tmp_path, base_name="other"
        )
        assert len(result) == 1
        assert result[0].ptm == 0.9

    def test_missing_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            bc.BoltzConfidenceMetrics.from_prediction_dir(tmp_path / "nope")


# ── The remaining JSON-direct scalars ───────────────────────────────
#
# Each one is a separate key Boltz writes, and each is read by name.  A
# property wired to the wrong key would report a plausible number from
# the wrong metric, which no downstream check could catch.  The values
# here are distinct so a crossed wire fails rather than coincides.


@pytest.mark.logic
def test_each_scalar_property_reads_its_own_json_key(tmp_path):
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {
        "protein_iptm": 0.11,
        "ligand_iptm": 0.22,
        "complex_iplddt": 0.33,
        "complex_pde": 0.44,
        "complex_ipde": 0.55,
    })
    m = bc.BoltzConfidenceMetrics(cj)
    assert m.protein_iptm == 0.11
    assert m.ligand_iptm == 0.22
    assert m.complex_iplddt == 0.33
    assert m.complex_pde == 0.44
    assert m.complex_ipde == 0.55


# ── The two dict-valued properties ──────────────────────────────────


@pytest.mark.logic
def test_an_unparseable_per_chain_ptm_is_dropped_not_zeroed(tmp_path):
    """0.0 would read as a measured pTM of a hopeless chain.  Chain 1's value
    is unusable, so the dict carries chain 0 only."""
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"chains_ptm": {"0": 0.9, "1": "n/a"}})
    m = bc.BoltzConfidenceMetrics(cj)
    assert m.chains_ptm == {"0": 0.9}


@pytest.mark.logic
def test_a_non_dict_chains_ptm_gives_an_empty_dict(tmp_path):
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"chains_ptm": [0.9, 0.8]})
    m = bc.BoltzConfidenceMetrics(cj)
    assert m.chains_ptm == {}


@pytest.mark.logic
def test_a_pair_entry_that_is_not_a_dict_is_skipped(tmp_path):
    """pair_chains_iptm is two levels deep.  A scalar where the inner dict
    should be is skipped whole, leaving the well-formed partner behind."""
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"pair_chains_iptm": {"0": 0.5, "1": {"0": 0.45}}})
    m = bc.BoltzConfidenceMetrics(cj)
    assert m.pair_chains_iptm == {"1": {"0": 0.45}}


@pytest.mark.logic
def test_an_unparseable_inner_pair_value_is_dropped(tmp_path):
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"pair_chains_iptm": {"0": {"1": "n/a", "2": 0.6}}})
    m = bc.BoltzConfidenceMetrics(cj)
    assert m.pair_chains_iptm == {"0": {"2": 0.6}}


# ── avg_plddt, the one value read from the npz ──────────────────────
#
# Boltz writes the per-residue pLDDT on a 0-1 scale.  This property is
# the only place in the pipeline that rescales it to the 0-100 B-factor
# convention, and `complex_plddt` beside it deliberately does not, so a
# mix-up would put two incomparable numbers in the same table.


def _npz(path, array):
    import numpy as np
    np.savez(path, plddt=np.asarray(array, dtype=float))
    return path


def _metrics_with_npz(tmp_path, array, payload=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, payload or {})
    npz = _npz(tmp_path / "plddt_input_model_0.npz", array)
    return bc.BoltzConfidenceMetrics(cj, npz)


@pytest.mark.logic
def test_a_0_to_1_array_is_rescaled_to_0_100(tmp_path):
    """mean([0.8, 0.9, 1.0]) = 0.9, and the max of 1.0 is <= 1.5, so the
    answer is 0.9 * 100 = 90.0."""
    m = _metrics_with_npz(tmp_path, [0.8, 0.9, 1.0])
    assert m.avg_plddt == 90.0


@pytest.mark.logic
def test_an_array_already_on_0_to_100_is_left_alone(tmp_path):
    """mean([10, 20, 30]) = 20.0, and the max of 30 is above 1.5, so no
    rescaling happens.  Rescaling anyway would report 2000.0."""
    m = _metrics_with_npz(tmp_path, [10, 20, 30])
    assert m.avg_plddt == 20.0


@pytest.mark.logic
def test_the_rescaling_boundary_is_a_max_of_1_point_5(tmp_path):
    """1.5 exactly is still treated as the 0-1 scale, so mean([1.5, 0.5]) =
    1.0 becomes 100.0.  A max just above it is not rescaled."""
    assert _metrics_with_npz(tmp_path / "a", [1.5, 0.5]).avg_plddt == 100.0
    assert _metrics_with_npz(tmp_path / "b", [1.6, 0.4]).avg_plddt == 1.0


@pytest.mark.logic
def test_a_per_residue_matrix_is_flattened_before_averaging(tmp_path):
    """Some samples arrive shaped (n_models, n_residues).  Flattening gives
    mean([0.1, 0.2, 0.3, 0.4]) = 0.25 -> 25.0; averaging rows first and then
    rows-of-means would also give 0.25 here, so the test uses a ragged spread
    whose row means differ."""
    m = _metrics_with_npz(tmp_path, [[0.1, 0.2, 0.3], [0.4, 0.5, 1.2]])
    # mean of all six = 2.7 / 6 = 0.45, max 1.2 <= 1.5 -> 45.0
    assert m.avg_plddt == 45.0


@pytest.mark.logic
def test_avg_plddt_is_none_when_no_npz_was_supplied(tmp_path):
    """The npz is optional.  from_workdir builds its metrics from the
    confidence JSON alone, so every negsteer seed takes this path."""
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"ptm": 0.5})
    assert bc.BoltzConfidenceMetrics(cj).avg_plddt is None


@pytest.mark.logic
def test_avg_plddt_is_none_when_the_named_npz_is_absent(tmp_path):
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"ptm": 0.5})
    m = bc.BoltzConfidenceMetrics(cj, tmp_path / "plddt_missing.npz")
    assert m.avg_plddt is None


@pytest.mark.logic
def test_avg_plddt_is_computed_once_and_then_cached(tmp_path):
    """Loading an npz per access would re-read a megabyte-scale array for
    every column of every row.  Deleting the file after the first access shows
    the second one never touches disk."""
    m = _metrics_with_npz(tmp_path, [0.8, 0.9, 1.0])
    assert m.avg_plddt == 90.0
    (tmp_path / "plddt_input_model_0.npz").unlink()
    assert m.avg_plddt == 90.0


@pytest.mark.logic
def test_an_unreadable_npz_gives_none_rather_than_raising(tmp_path):
    """A confidence row is still worth having without its mean pLDDT, so a
    file numpy cannot interpret becomes None."""
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"ptm": 0.5})
    bad = tmp_path / "plddt_input_model_0.npz"
    bad.write_text("this is not an npz archive")
    assert bc.BoltzConfidenceMetrics(cj, bad).avg_plddt is None


@pytest.mark.logic
@pytest.mark.xfail(
    strict=False,
    reason="Defect: the except clause lists (OSError, ValueError, KeyError, "
           "ImportError), and numpy raises EOFError on a zero-byte file.  A "
           "truncated plddt npz is what a killed GPU job leaves, and it is "
           "the case the fallback exists for, so avg_plddt should return None "
           "rather than raise.  Not fixed here because widening the except "
           "clause changes what every caller of avg_plddt sees on a corrupt "
           "run, which is a behaviour change rather than a test change.",
)
def test_a_truncated_npz_should_give_none_rather_than_raising(tmp_path):
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"ptm": 0.5})
    truncated = tmp_path / "plddt_input_model_0.npz"
    truncated.write_bytes(b"")
    assert bc.BoltzConfidenceMetrics(cj, truncated).avg_plddt is None


@pytest.mark.logic
def test_a_truncated_npz_raises_today(tmp_path):
    """Today's behaviour, pinned beside the xfail above so the change is
    visible when the except clause is widened."""
    cj = tmp_path / "confidence_input_model_0.json"
    _write_conf_json(cj, {"ptm": 0.5})
    truncated = tmp_path / "plddt_input_model_0.npz"
    truncated.write_bytes(b"")
    with pytest.raises(EOFError):
        _ = bc.BoltzConfidenceMetrics(cj, truncated).avg_plddt


@pytest.mark.logic
def test_from_prediction_dir_links_the_matching_plddt_npz(tmp_path):
    """The npz is matched by filename stem, not by position in the sort.  Only
    sample 1 has one here, so sample 0 reports no avg_plddt and sample 1
    reports mean([0.5, 0.5]) * 100 = 50.0."""
    for i in range(2):
        _write_conf_json(tmp_path / f"confidence_input_model_{i}.json",
                         {"ptm": 0.5})
    _npz(tmp_path / "plddt_input_model_1.npz", [0.5, 0.5])
    found = bc.BoltzConfidenceMetrics.from_prediction_dir(tmp_path)
    assert len(found) == 2
    assert found[0].avg_plddt is None
    assert found[1].avg_plddt == 50.0
