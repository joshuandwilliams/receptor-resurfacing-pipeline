"""Unit tests for bin/negative_steering_run.py (Phase 4 Tier 5).

Uses MagicMock PSPs (same pattern as test_stage_result.py) so the cross-
stage aggregation logic is testable without real PDBs.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import negative_steering_run as nsr  # noqa: E402
import pipeline_thresholds as pt  # noqa: E402
import stage_result as sr  # noqa: E402
from position_set import PositionSet  # noqa: E402


def _fake_pred(ra_eff, intact_rmsd, contact_positions=()):
    m = MagicMock()
    m.ra_eff_vs.return_value = ra_eff
    m.independent_receptor_rmsd_vs.return_value = intact_rmsd
    cs = PositionSet(list(contact_positions), chain="A", frame="prediction")
    m.receptor_contact_residues.return_value = cs
    return m


def _fake_conf():
    return MagicMock()


def _make_cold_start(predictions, config_id="design_03_seq_01"):
    return sr.StageResult(
        stage_type="cold_start", config_id=config_id,
        predictions=predictions,
        confidences=[_fake_conf() for _ in predictions],
        seed_indices=list(range(len(predictions))))


def _make_steered(predictions, design_id="design_03",
                  mutated_positions=None):
    if mutated_positions is None:
        mutated_positions = PositionSet([5, 7, 12], chain="A", frame="prediction")
    return sr.StageResult(
        stage_type="steered", config_id=design_id,
        predictions=predictions,
        confidences=[_fake_conf() for _ in predictions],
        seed_indices=list(range(len(predictions))),
        mutated_positions=mutated_positions)


def _make_reversion(predictions, rev_id="rev_design_03_abc",
                    applies_to=("design_03",),
                    mutated_positions=None):
    if mutated_positions is None:
        mutated_positions = PositionSet([5, 7, 12], chain="A", frame="prediction")
    return sr.StageResult(
        stage_type="reversion", config_id=rev_id,
        predictions=predictions,
        confidences=[_fake_conf() for _ in predictions],
        seed_indices=list(range(len(predictions))),
        mutated_positions=mutated_positions,
        applies_to_designs=list(applies_to))


@pytest.fixture
def truth():
    return MagicMock()


@pytest.fixture
def contam():
    return PositionSet([5, 7, 12], chain="A", frame="prediction")


# ── Construction ───────────────────────────────────────────────────


@pytest.mark.logic
class TestConstruction:
    def test_basic(self, tmp_path, truth):
        cs = _make_cold_start([_fake_pred(2.0, 1.0)] * 3)
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="design_03_seq_01",
            workdir=tmp_path,
            row_type="steered",
            cold_start=cs,
            truth=truth,
        )
        assert run.mpnn_sequence_id == "design_03_seq_01"
        assert run.num_seeds == 3

    def test_empty_id_raises(self, tmp_path, truth):
        with pytest.raises(ValueError, match="mpnn_sequence_id"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="", workdir=tmp_path,
                row_type="steered",
                cold_start=_make_cold_start([_fake_pred(2.0, 1.0)]),
                truth=truth)

    def test_invalid_row_type_raises(self, tmp_path, truth):
        with pytest.raises(ValueError, match="row_type"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path,
                row_type="bogus",
                cold_start=_make_cold_start([_fake_pred(2.0, 1.0)]),
                truth=truth)

    def test_wrong_cold_start_type_raises(self, tmp_path, truth):
        wrong = _make_steered([_fake_pred(2.0, 1.0)])
        with pytest.raises(ValueError, match="cold_start"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path,
                row_type="steered", cold_start=wrong, truth=truth)


# ── Cold-start path (no steering ran) ──────────────────────────────


@pytest.mark.logic
class TestColdStartPath:
    def test_all_clean_cold_start_is_tier_A(self, tmp_path, truth):
        cs = _make_cold_start([_fake_pred(2.0, 1.0)] * 3)
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, truth=truth)
        t = pt.PipelineInternalThresholds.default()
        assert run.outcome(t) == "no_reversion"
        assert run.n_pass(t) == 3
        assert run.n_seeds_total() == 3
        assert run.tier(t) == "A"

    def test_partial_cold_start_pass_is_tier_B(self, tmp_path, truth):
        cs = _make_cold_start(
            [_fake_pred(2.0, 1.0), _fake_pred(10.0, 1.0), _fake_pred(2.0, 1.0)])
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, truth=truth)
        t = pt.PipelineInternalThresholds.default()
        # 2 clean, 1 fail_structural — but outcome is still no_reversion
        # because steering didn't run; n_pass=2 (only 'clean' counts)
        assert run.outcome(t) == "no_reversion"
        assert run.n_pass(t) == 2
        assert run.tier(t) == "B"

    def test_zero_clean_cold_start_is_tier_none(self, tmp_path, truth):
        cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, truth=truth)
        t = pt.PipelineInternalThresholds.default()
        # The "wrong placement no contamination" case from CL-3 / CL-1
        # outcome=no_reversion BUT tier=none — the bug fix
        assert run.outcome(t) == "no_reversion"
        assert run.n_pass(t) == 0
        assert run.tier(t) == "none"


# ── Steered + reversion path ────────────────────────────────────


@pytest.mark.logic
class TestSteeredOnlyPath:
    """Steering ran but no reversion was triggered (CL-3 rule didn't fire)."""

    def test_all_clean_steered_3of3_tier_A(self, tmp_path, truth):
        cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)  # failed cold start
        steered = {
            "design_03": _make_steered(
                [_fake_pred(2.0, 1.0, contact_positions=[100])] * 3),
        }
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, steered=steered, truth=truth,
            contamination_positions=PositionSet(
                [5, 7, 12], chain="A", frame="prediction"))
        t = pt.PipelineInternalThresholds.default()
        assert run.outcome(t) == "no_reversion"  # reversion didn't run
        assert run.n_pass(t) == 3
        assert run.tier(t) == "A"


@pytest.mark.logic
class TestReversionPath:
    """Steering ran + reversion ran."""

    def test_reversion_pose_holds_is_tier_A(self, tmp_path, truth, contam):
        cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
        # 3 contaminated steered → reversion triggers
        steered = {
            "design_03": _make_steered(
                [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3),
        }
        # Reversion all pose_holds
        reversion = {
            "rev_design_03_abc": _make_reversion(
                [_fake_pred(2.0, 1.0, contact_positions=[100])] * 3,
                applies_to=("design_03",)),
        }
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, steered=steered, reversion=reversion,
            truth=truth, contamination_positions=contam)
        t = pt.PipelineInternalThresholds.default()
        verdicts = run.per_seed_final_verdicts(t)
        # All 3 should be 'reverted' stage with pose_holds verdict
        assert all(e["stage_run"] == "reverted" for e in verdicts)
        assert all(e["verdict"] == "pose_holds" for e in verdicts)
        assert run.outcome(t) == "pose_holds"
        assert run.n_pass(t) == 3
        assert run.tier(t) == "A"

    def test_reversion_new_contamination_is_tier_none(
        self, tmp_path, truth, contam
    ):
        cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
        steered = {
            "design_03": _make_steered(
                [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3),
        }
        # Reversion all new_contamination
        reversion = {
            "rev_design_03_abc": _make_reversion(
                [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3,
                applies_to=("design_03",)),
        }
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, steered=steered, reversion=reversion,
            truth=truth, contamination_positions=contam)
        t = pt.PipelineInternalThresholds.default()
        assert run.outcome(t) == "new_contamination"
        assert run.n_pass(t) == 0
        assert run.tier(t) == "none"

    def test_mixed_reversion_outcomes(self, tmp_path, truth, contam):
        """2 pose_holds, 1 new_contamination → mixed; majority is pose_holds."""
        cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
        steered = {
            "design_03": _make_steered(
                [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3),
        }
        reversion = {
            "rev_design_03_abc": _make_reversion(
                [
                    _fake_pred(2.0, 1.0, contact_positions=[100]),  # pose_holds
                    _fake_pred(2.0, 1.0, contact_positions=[100]),  # pose_holds
                    _fake_pred(2.0, 1.0, contact_positions=[5]),    # new_contam
                ],
                applies_to=("design_03",)),
        }
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, steered=steered, reversion=reversion,
            truth=truth, contamination_positions=contam)
        t = pt.PipelineInternalThresholds.default()
        assert run.outcome(t) == "pose_holds"  # 2 pass > 1 new_contam
        assert run.n_pass(t) == 2
        assert run.tier(t) == "B"


@pytest.mark.logic
class TestReversionForDesign:
    def test_finds_matching_reversion(self, tmp_path, truth, contam):
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        rev_a = _make_reversion([_fake_pred(2.0, 1.0)],
                                rev_id="rev_a", applies_to=("design_03",))
        rev_b = _make_reversion([_fake_pred(2.0, 1.0)],
                                rev_id="rev_b", applies_to=("design_04",))
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, reversion={"rev_a": rev_a, "rev_b": rev_b},
            truth=truth, contamination_positions=contam)
        assert run._reversion_for_design("design_03") is rev_a
        assert run._reversion_for_design("design_04") is rev_b
        assert run._reversion_for_design("design_99") is None

    def test_dedup_one_reversion_applies_to_many(self, tmp_path, truth, contam):
        """Same reversion attempt represents multiple steered designs."""
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        rev = _make_reversion([_fake_pred(2.0, 1.0)],
                              applies_to=("design_03", "design_05"))
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, reversion={"rev": rev},
            truth=truth, contamination_positions=contam)
        # Both designs map to the same reversion
        assert run._reversion_for_design("design_03") is rev
        assert run._reversion_for_design("design_05") is rev


@pytest.mark.logic
class TestImmutability:
    def test_frozen(self, tmp_path, truth):
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        run = nsr.NegativeSteeringRun(
            mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
            cold_start=cs, truth=truth)
        with pytest.raises((AttributeError, Exception)):
            run.mpnn_sequence_id = "new"  # type: ignore


# ── What the stage dictionaries are allowed to hold ────────────────
#
# Every verdict, outcome and tier on this type is computed by walking
# `cold_start`, `steered` and `reversion`.  A stage filed under the wrong
# key silently changes which rule is applied to it — a steered stage read
# as a reversion would turn `contaminated` into `new_contamination` — so
# the checks run at construction rather than at first use.


@pytest.mark.logic
class TestStageDictValidation:
    def test_a_cold_start_that_is_not_a_stage_result_is_rejected(self, tmp_path, truth):
        with pytest.raises(TypeError, match="cold_start must be a StageResult"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
                cold_start="cycle_0/initial", truth=truth)

    def test_a_steered_entry_that_is_not_a_stage_result_is_rejected(
        self, tmp_path, truth
    ):
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        with pytest.raises(TypeError, match="steered\\['design_03'\\]"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
                cold_start=cs, steered={"design_03": "design_03_s0"},
                truth=truth)

    def test_a_non_steered_stage_filed_under_steered_is_rejected(
        self, tmp_path, truth
    ):
        """A cold-start stage in the steered dict would be scored with the
        steered vocabulary (clean_steered / contaminated) instead of its
        own."""
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        other = _make_cold_start([_fake_pred(2.0, 1.0)])
        with pytest.raises(ValueError, match="stage_type='steered'"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
                cold_start=cs, steered={"design_03": other}, truth=truth)

    def test_a_reversion_entry_that_is_not_a_stage_result_is_rejected(
        self, tmp_path, truth
    ):
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        with pytest.raises(TypeError, match="reversion\\['rev_design_03_s0'\\]"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
                cold_start=cs, reversion={"rev_design_03_s0": None},
                truth=truth)

    def test_a_steered_stage_filed_under_reversion_is_rejected(self, tmp_path, truth):
        """Mis-filing here is the dangerous direction: a steered stage read as
        a reversion would report pose_holds for a seed that was never
        reverted."""
        cs = _make_cold_start([_fake_pred(2.0, 1.0)])
        steered_stage = _make_steered([_fake_pred(2.0, 1.0)])
        with pytest.raises(ValueError, match="stage_type='reversion'"):
            nsr.NegativeSteeringRun(
                mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
                cold_start=cs, reversion={"rev_design_03_s0": steered_stage},
                truth=truth)


# ── Truth is required before any verdict ───────────────────────────


@pytest.mark.logic
def test_a_run_without_truth_raises_rather_than_scoring_against_nothing(tmp_path):
    """Every verdict is a comparison against the RFDiffusion design that was
    being reproduced.  Without it ra_eff and the intact-receptor RMSD are both
    None, which would read as no_data for every seed and report tier none for
    a run nobody could score."""
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=_make_cold_start([_fake_pred(2.0, 1.0)] * 3))
    t = pt.PipelineInternalThresholds.default()
    with pytest.raises(ValueError, match="require a 'truth'"):
        run.per_seed_final_verdicts(t)
    with pytest.raises(ValueError, match="require a 'truth'"):
        run.n_pass(t)


# ── A reversion that did not cover every steered seed ──────────────


@pytest.mark.logic
def test_a_seed_the_reversion_skipped_keeps_its_steered_verdict(tmp_path, truth, contam):
    """The engine reverts per (design, source seed), so a reversion can come
    back with fewer seeds than the steered stage it answers.  Seeds 0 and 1
    were reverted and hold the pose; seed 2 has no reversion, so its steered
    verdict stands and it stays clean_steered.  3 of 3 pass-equivalent, which
    is tier A.
    """
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {
        "design_03": sr.StageResult(
            stage_type="steered", config_id="design_03",
            predictions=[
                _fake_pred(2.0, 1.0, contact_positions=[5]),    # contaminated
                _fake_pred(2.0, 1.0, contact_positions=[5]),    # contaminated
                _fake_pred(2.0, 1.0, contact_positions=[100]),  # clean_steered
            ],
            confidences=[_fake_conf() for _ in range(3)],
            seed_indices=[0, 1, 2],
            mutated_positions=PositionSet([5, 7, 12], chain="A",
                                          frame="prediction")),
    }
    reversion = {
        "rev_design_03_s0": sr.StageResult(
            stage_type="reversion", config_id="rev_design_03_s0",
            predictions=[_fake_pred(2.0, 1.0, contact_positions=[100])] * 2,
            confidences=[_fake_conf() for _ in range(2)],
            seed_indices=[0, 1],                       # seed 2 is absent
            mutated_positions=PositionSet([5, 7, 12], chain="A",
                                          frame="prediction"),
            applies_to_designs=["design_03"]),
    }
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, steered=steered, reversion=reversion,
        truth=truth, contamination_positions=contam)
    t = pt.PipelineInternalThresholds.default()

    by_seed = {e["seed_index"]: e for e in run.per_seed_final_verdicts(t)}
    assert by_seed[0]["stage_run"] == "reverted"
    assert by_seed[1]["stage_run"] == "reverted"
    assert by_seed[2]["stage_run"] == "steered"
    assert by_seed[2]["verdict"] == "clean_steered"
    assert run.n_pass(t) == 3
    assert run.tier(t) == "A"


# ── outcome: the failure branches ──────────────────────────────────


@pytest.mark.logic
def test_every_seed_collapsing_after_reversion_is_pose_collapses(
    tmp_path, truth, contam
):
    """Nothing held and nothing was newly contaminated, so the only verdict
    left is pose_collapses: 3 collapses, 0 holds, 0 new contamination."""
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {"design_03": _make_steered(
        [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3)}
    reversion = {"rev_design_03_s0": _make_reversion(
        [_fake_pred(10.0, 1.0, contact_positions=[100])] * 3,
        applies_to=("design_03",))}
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, steered=steered, reversion=reversion,
        truth=truth, contamination_positions=contam)
    t = pt.PipelineInternalThresholds.default()
    assert run.outcome(t) == "pose_collapses"
    assert run.n_pass(t) == 0
    assert run.tier(t) == "none"


# ── outcome_reason ─────────────────────────────────────────────────
#
# The reason string is what a reader sees next to an outcome in the cohort
# table.  An outcome of no_reversion has two quite different causes, and
# only this string separates them.


@pytest.mark.logic
def test_the_reason_for_a_cold_start_run_says_steering_never_ran(tmp_path, truth):
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=_make_cold_start([_fake_pred(2.0, 1.0)] * 3), truth=truth)
    t = pt.PipelineInternalThresholds.default()
    assert run.outcome(t) == "no_reversion"
    assert run.outcome_reason(t) == "cold-start path: no steering ran"


@pytest.mark.logic
def test_the_reason_separates_steering_without_reversion_from_a_cold_start(
    tmp_path, truth, contam
):
    """Both report outcome no_reversion, and they are not the same thing.
    Here steering ran and no design had a majority of its correctly-placed
    seeds contaminated, so no reversion was triggered."""
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {"design_03": _make_steered(
        [_fake_pred(2.0, 1.0, contact_positions=[100])] * 3)}
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, steered=steered, truth=truth,
        contamination_positions=contam)
    t = pt.PipelineInternalThresholds.default()
    assert run.outcome(t) == "no_reversion"
    reason = run.outcome_reason(t)
    assert "no reversion" in reason
    assert reason != "cold-start path: no steering ran"


@pytest.mark.logic
def test_the_reason_for_a_reverted_run_is_the_verdict_breakdown(
    tmp_path, truth, contam
):
    """2 seeds hold and 1 is newly contaminated, so the breakdown reads
    'new_contamination=1; pose_holds=2' — sorted by verdict name."""
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {"design_03": _make_steered(
        [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3)}
    reversion = {"rev_design_03_s0": _make_reversion(
        [
            _fake_pred(2.0, 1.0, contact_positions=[100]),
            _fake_pred(2.0, 1.0, contact_positions=[100]),
            _fake_pred(2.0, 1.0, contact_positions=[5]),
        ],
        applies_to=("design_03",))}
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, steered=steered, reversion=reversion,
        truth=truth, contamination_positions=contam)
    t = pt.PipelineInternalThresholds.default()
    assert run.outcome_reason(t) == "new_contamination=1; pose_holds=2"


# ── tier arithmetic ────────────────────────────────────────────────


@pytest.mark.logic
def test_tier_C_is_exactly_one_passing_seed_out_of_more_than_one(
    tmp_path, truth, contam
):
    """1 of 3 is tier C.  2 of 3 would be B and 3 of 3 A, so C is the only
    tier that pins n_pass == 1 rather than a fraction."""
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {"design_03": _make_steered([
        _fake_pred(2.0, 1.0, contact_positions=[100]),  # clean_steered
        _fake_pred(2.0, 1.0, contact_positions=[5]),    # contaminated
        _fake_pred(10.0, 1.0),                          # pose_collapses
    ])}
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, steered=steered, truth=truth,
        contamination_positions=contam)
    t = pt.PipelineInternalThresholds.default()
    assert run.n_pass(t) == 1
    assert run.n_seeds_total() == 3
    assert run.tier(t) == "C"


@pytest.mark.logic
def test_the_steered_path_totals_seeds_across_designs(tmp_path, truth, contam):
    """Two designs of 3 seeds each is 6 seeds, not 3.  Tiering is per MPNN
    sequence over that product, so 3 passing of 6 is tier B."""
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {
        "design_03": _make_steered(
            [_fake_pred(2.0, 1.0, contact_positions=[100])] * 3,
            design_id="design_03"),
        "design_04": _make_steered(
            [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3,
            design_id="design_04"),
    }
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, steered=steered, truth=truth,
        contamination_positions=contam)
    t = pt.PipelineInternalThresholds.default()
    assert run.n_seeds_total() == 6
    assert run.n_pass(t) == 3
    assert run.tier(t) == "B"


class _ZeroSeedRun(nsr.NegativeSteeringRun):
    """A run reporting no seeds at all.

    No workdir can produce one: StageResult requires at least one prediction,
    and n_seeds_total is a sum of those counts, so it is always >= 1.  The
    `n_seeds <= 0` guard in tier() is therefore only reachable through a
    subclass, and it is worth keeping covered because a zero there would
    otherwise make `n_pass == n_seeds` read as tier A.
    """

    def n_seeds_total(self) -> int:
        return 0


@pytest.mark.logic
def test_a_run_with_no_seeds_is_tier_none_not_tier_A(tmp_path, truth):
    cs = _make_cold_start([_fake_pred(2.0, 1.0)] * 3)
    run = _ZeroSeedRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=cs, truth=truth)
    t = pt.PipelineInternalThresholds.default()
    assert run.n_pass(t) == 3          # the seeds exist in the stage
    assert run.n_seeds_total() == 0     # but the total reports none
    assert run.tier(t) == "none"


@pytest.mark.logic
def test_no_real_run_can_report_zero_seeds(tmp_path, truth):
    """The companion to the test above: the guard is unreachable from data,
    because a constructible run always carries at least one prediction."""
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="x", workdir=tmp_path, row_type="steered",
        cold_start=_make_cold_start([_fake_pred(2.0, 1.0)]), truth=truth)
    assert run.n_seeds_total() >= 1


# ── Identification ─────────────────────────────────────────────────


@pytest.mark.logic
def test_the_repr_names_the_sequence_and_counts_the_stages(tmp_path, truth, contam):
    """The repr is what appears in a pytest failure and in a cohort log line,
    so it has to say which sequence and how much of the chain ran."""
    cs = _make_cold_start([_fake_pred(10.0, 1.0)] * 3)
    steered = {
        "design_03": _make_steered([_fake_pred(2.0, 1.0)], design_id="design_03"),
        "design_04": _make_steered([_fake_pred(2.0, 1.0)], design_id="design_04"),
    }
    reversion = {"rev_design_03_s0": _make_reversion([_fake_pred(2.0, 1.0)])}
    run = nsr.NegativeSteeringRun(
        mpnn_sequence_id="design_03_seq_01", workdir=tmp_path,
        row_type="control_polyA", cold_start=cs, steered=steered,
        reversion=reversion, truth=truth, contamination_positions=contam)
    assert repr(run) == (
        "NegativeSteeringRun(mpnn_sequence_id='design_03_seq_01', "
        "row_type='control_polyA', n_steered_designs=2, n_reversions=1)"
    )
