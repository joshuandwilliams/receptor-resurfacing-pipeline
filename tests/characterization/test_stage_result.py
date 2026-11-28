"""Unit tests for bin/stage_result.py (Phase 4 Tier 4.2).

The verdict logic depends on ra_eff/intact/contact behaviour of
ProteinStructurePrediction, which need real PDB files. To keep the
Tier 4 tests local-only and fast, we use a TestablePSP stand-in that
returns hand-injected metric values. The real PSP wiring is
characterised at hpc-tier when negsteer runs end-to-end.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import pipeline_thresholds as pt  # noqa: E402
import stage_result as sr  # noqa: E402
from position_set import PositionSet  # noqa: E402


def _fake_pred(ra_eff, intact_rmsd, contact_positions=()):
    """Return a MagicMock that quacks like a ProteinStructurePrediction
    for the methods StageResult uses."""
    m = MagicMock()
    m.ra_eff_vs.return_value = ra_eff
    m.independent_receptor_rmsd_vs.return_value = intact_rmsd
    cs = PositionSet(
        positions=list(contact_positions),
        chain="A", frame="prediction",
    )
    m.receptor_contact_residues.return_value = cs
    # spec=ProteinStructurePrediction is not enforced — the type check
    # for StageResult is on the LIST, not on each element.
    return m


def _fake_conf():
    return MagicMock()


def _truth():
    return MagicMock()


# ── Construction validation ────────────────────────────────────────


@pytest.mark.logic
class TestConstruction:
    def test_basic_cold_start(self):
        s = sr.StageResult(
            stage_type="cold_start", config_id="design_03_seq_01",
            predictions=[_fake_pred(2.0, 1.0)],
            confidences=[_fake_conf()],
            seed_indices=[0])
        assert s.stage_type == "cold_start"
        assert s.num_seeds == 1

    def test_invalid_stage_type_raises(self):
        with pytest.raises(ValueError, match="stage_type"):
            sr.StageResult(stage_type="bogus", config_id="c",
                           predictions=[_fake_pred(2.0, 1.0)],
                           confidences=[_fake_conf()], seed_indices=[0])

    def test_empty_config_id_raises(self):
        with pytest.raises(ValueError, match="config_id"):
            sr.StageResult(stage_type="cold_start", config_id="",
                           predictions=[_fake_pred(2.0, 1.0)],
                           confidences=[_fake_conf()], seed_indices=[0])

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="confidences length"):
            sr.StageResult(stage_type="cold_start", config_id="c",
                           predictions=[_fake_pred(2.0, 1.0)],
                           confidences=[_fake_conf(), _fake_conf()],
                           seed_indices=[0])

    def test_empty_predictions_raises(self):
        with pytest.raises(ValueError, match="at least 1"):
            sr.StageResult(stage_type="cold_start", config_id="c",
                           predictions=[], confidences=[], seed_indices=[])

    def test_cold_start_with_mutated_positions_raises(self):
        mp = PositionSet([5], chain="A", frame="prediction")
        with pytest.raises(ValueError, match="cold_start"):
            sr.StageResult(stage_type="cold_start", config_id="c",
                           predictions=[_fake_pred(2.0, 1.0)],
                           confidences=[_fake_conf()], seed_indices=[0],
                           mutated_positions=mp)

    def test_steered_without_mutated_positions_raises(self):
        with pytest.raises(ValueError, match="requires mutated_positions"):
            sr.StageResult(stage_type="steered", config_id="c",
                           predictions=[_fake_pred(2.0, 1.0)],
                           confidences=[_fake_conf()], seed_indices=[0])

    def test_applies_to_designs_only_for_reversion(self):
        mp = PositionSet([5], chain="A", frame="prediction")
        with pytest.raises(ValueError, match="reversion-only"):
            sr.StageResult(stage_type="steered", config_id="c",
                           predictions=[_fake_pred(2.0, 1.0)],
                           confidences=[_fake_conf()], seed_indices=[0],
                           mutated_positions=mp,
                           applies_to_designs=["design_03"])


# ── Per-seed verdicts ──────────────────────────────────────────────


def _make_cold_start(predictions):
    return sr.StageResult(
        stage_type="cold_start", config_id="design_03_seq_01",
        predictions=predictions,
        confidences=[_fake_conf() for _ in predictions],
        seed_indices=list(range(len(predictions))))


def _make_steered(predictions, mutated_positions):
    return sr.StageResult(
        stage_type="steered", config_id="design_03",
        predictions=predictions,
        confidences=[_fake_conf() for _ in predictions],
        seed_indices=list(range(len(predictions))),
        mutated_positions=mutated_positions)


def _make_reversion(predictions, mutated_positions, applies_to=("design_03",)):
    return sr.StageResult(
        stage_type="reversion", config_id="rev_design_03_abc",
        predictions=predictions,
        confidences=[_fake_conf() for _ in predictions],
        seed_indices=list(range(len(predictions))),
        mutated_positions=mutated_positions,
        applies_to_designs=list(applies_to))


@pytest.mark.logic
class TestColdStartVerdicts:
    def test_all_clean(self):
        # ra_eff < 5, intact_rmsd < 5 → "clean"
        preds = [_fake_pred(2.0, 1.0)] * 3
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth()) == ["clean"] * 3

    def test_fail_structural_ra_eff(self):
        # ra_eff >= 5 → fail
        preds = [_fake_pred(10.0, 1.0)] * 3
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth()) == ["fail_structural"] * 3

    def test_fail_structural_intact(self):
        # intact_rmsd >= 5 → fail (receptor not intact)
        preds = [_fake_pred(2.0, 10.0)] * 3
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth()) == ["fail_structural"] * 3

    def test_no_data_when_metrics_none(self):
        preds = [_fake_pred(None, None)]
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth()) == ["no_data"]


@pytest.mark.logic
class TestSteeredVerdicts:
    def _setup(self):
        mp = PositionSet([5, 7, 12], chain="A", frame="prediction")
        contam = PositionSet([5, 7, 12], chain="A", frame="prediction")
        return mp, contam

    def test_clean_steered_no_contamination(self):
        mp, contam = self._setup()
        # No contacts on mutated positions → clean_steered
        preds = [_fake_pred(2.0, 1.0, contact_positions=[100, 101])] * 3
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth(), contam) == ["clean_steered"] * 3

    def test_contaminated(self):
        mp, contam = self._setup()
        # contact at position 5 (in both mutated_positions and contam) → contaminated
        preds = [_fake_pred(2.0, 1.0, contact_positions=[5, 100])] * 3
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth(), contam) == ["contaminated"] * 3

    def test_pose_collapses_structural_fail(self):
        mp, contam = self._setup()
        preds = [_fake_pred(10.0, 1.0)]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth(), contam) == ["pose_collapses"]


@pytest.mark.logic
class TestReversionVerdicts:
    def _setup(self):
        mp = PositionSet([5, 7, 12], chain="A", frame="prediction")
        contam = PositionSet([5, 7, 12], chain="A", frame="prediction")
        return mp, contam

    def test_pose_holds(self):
        mp, contam = self._setup()
        preds = [_fake_pred(2.0, 1.0, contact_positions=[100])] * 3
        s = _make_reversion(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth(), contam) == ["pose_holds"] * 3

    def test_new_contamination(self):
        mp, contam = self._setup()
        preds = [_fake_pred(2.0, 1.0, contact_positions=[5])] * 3
        s = _make_reversion(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth(), contam) == ["new_contamination"] * 3

    def test_pose_collapses(self):
        mp, contam = self._setup()
        preds = [_fake_pred(10.0, 1.0, contact_positions=[100])]
        s = _make_reversion(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.per_seed_verdicts(t, _truth(), contam) == ["pose_collapses"]


# ── CL-3 gating rule (the headline test) ──────────────────────────


@pytest.mark.logic
class TestCL3SteeredGating:
    """CL-3 rule: reversion runs iff n_contaminated >= ceil(n_correctly_placed / 2)."""

    def _setup(self):
        mp = PositionSet([5, 7, 12], chain="A", frame="prediction")
        contam = PositionSet([5, 7, 12], chain="A", frame="prediction")
        return mp, contam

    def test_3of3_correct_1_contaminated_no_reversion(self):
        """3/3 correctly placed, 1 contaminated → 1 < ceil(3/2)=2 → NO reversion."""
        mp, contam = self._setup()
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[100]),     # clean
            _fake_pred(2.0, 1.0, contact_positions=[101]),     # clean
            _fake_pred(2.0, 1.0, contact_positions=[5]),       # contaminated
        ]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.n_correctly_placed(t, _truth()) == 3
        assert s.n_contaminated(t, _truth(), contam) == 1
        assert s.triggers_next_stage(t, _truth(), contam) is False

    def test_3of3_correct_2_contaminated_reversion_runs(self):
        """3/3 correctly placed, 2 contaminated → 2 >= 2 → reversion."""
        mp, contam = self._setup()
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[100]),
            _fake_pred(2.0, 1.0, contact_positions=[5]),
            _fake_pred(2.0, 1.0, contact_positions=[7]),
        ]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.triggers_next_stage(t, _truth(), contam) is True

    def test_2of3_correct_1_contaminated_reversion_runs(self):
        """2/3 correct (1 wrong placement), 1 contaminated → 1 >= ceil(2/2)=1 → reversion."""
        mp, contam = self._setup()
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[100]),  # clean
            _fake_pred(2.0, 1.0, contact_positions=[5]),    # contaminated
            _fake_pred(10.0, 1.0),                          # wrong placement
        ]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.n_correctly_placed(t, _truth()) == 2
        assert s.n_contaminated(t, _truth(), contam) == 1
        assert s.triggers_next_stage(t, _truth(), contam) is True

    def test_2of3_correct_0_contaminated_no_reversion(self):
        """2/3 correctly placed, 0 contaminated → no reversion (tier-B holds)."""
        mp, contam = self._setup()
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[100]),
            _fake_pred(2.0, 1.0, contact_positions=[101]),
            _fake_pred(10.0, 1.0),
        ]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.triggers_next_stage(t, _truth(), contam) is False

    def test_0_correctly_placed_no_reversion(self):
        """0 correctly placed → nothing to revert."""
        mp, contam = self._setup()
        preds = [_fake_pred(10.0, 1.0)] * 3
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.n_correctly_placed(t, _truth()) == 0
        assert s.triggers_next_stage(t, _truth(), contam) is False

    def test_3of3_all_contaminated_reversion_runs(self):
        mp, contam = self._setup()
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[5]),
            _fake_pred(2.0, 1.0, contact_positions=[7]),
            _fake_pred(2.0, 1.0, contact_positions=[12]),
        ]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.triggers_next_stage(t, _truth(), contam) is True


# ── Other gating ─────────────────────────────────────────────────


@pytest.mark.logic
class TestColdStartGating:
    def test_all_clean_skips_steering(self):
        preds = [_fake_pred(2.0, 1.0)] * 3
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        # 3/3 clean → fewer than majority NOT clean → False (skip steering)
        assert s.triggers_next_stage(t, _truth()) is False

    def test_majority_fail_runs_steering(self):
        preds = [_fake_pred(2.0, 1.0), _fake_pred(10.0, 1.0), _fake_pred(10.0, 1.0)]
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        assert s.triggers_next_stage(t, _truth()) is True


@pytest.mark.logic
class TestReversionGating:
    def test_reversion_never_triggers_next_stage(self):
        mp = PositionSet([5], chain="A", frame="prediction")
        contam = PositionSet([5], chain="A", frame="prediction")
        preds = [_fake_pred(2.0, 1.0, contact_positions=[100])]
        s = _make_reversion(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.triggers_next_stage(t, _truth(), contam) is False


# ── n_pass ──────────────────────────────────────────────────────


@pytest.mark.logic
class TestNPass:
    def test_cold_start_pass_counts_clean(self):
        preds = [_fake_pred(2.0, 1.0), _fake_pred(10.0, 1.0), _fake_pred(2.0, 1.0)]
        s = _make_cold_start(preds)
        t = pt.PipelineInternalThresholds.default()
        assert s.n_pass(t, _truth()) == 2

    def test_steered_pass_counts_clean_steered(self):
        mp = PositionSet([5], chain="A", frame="prediction")
        contam = PositionSet([5], chain="A", frame="prediction")
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[100]),  # clean_steered
            _fake_pred(2.0, 1.0, contact_positions=[5]),    # contaminated
            _fake_pred(10.0, 1.0),                          # pose_collapses
        ]
        s = _make_steered(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.n_pass(t, _truth(), contam) == 1

    def test_reversion_pass_counts_pose_holds(self):
        mp = PositionSet([5], chain="A", frame="prediction")
        contam = PositionSet([5], chain="A", frame="prediction")
        preds = [
            _fake_pred(2.0, 1.0, contact_positions=[100]),  # pose_holds
            _fake_pred(2.0, 1.0, contact_positions=[5]),    # new_contamination
        ]
        s = _make_reversion(preds, mp)
        t = pt.PipelineInternalThresholds.default()
        assert s.n_pass(t, _truth(), contam) == 1


# ── The three parallel tuples ──────────────────────────────────────
#
# predictions, confidences and seed_indices are read positionally
# everywhere downstream.  A length mismatch would pair a prediction with
# another seed's confidence, and nothing later could detect it.


@pytest.mark.logic
def test_a_seed_index_per_prediction_is_required():
    with pytest.raises(ValueError, match="seed_indices length"):
        sr.StageResult(stage_type="cold_start", config_id="c",
                       predictions=[_fake_pred(2.0, 1.0)] * 3,
                       confidences=[_fake_conf() for _ in range(3)],
                       seed_indices=[0, 1])


# ── contamination_positions is required where it is used ───────────
#
# Both steered and reversion verdicts turn on the triple intersection
# with the contamination set.  Defaulting it to empty would make every
# seed read clean, so a caller who forgot it gets an error instead of a
# clean sheet.


@pytest.mark.logic
@pytest.mark.parametrize("stage_type", ["steered", "reversion"])
def test_verdicts_without_contamination_positions_raise(stage_type):
    mp = PositionSet([5], chain="A", frame="prediction")
    kwargs = {"applies_to_designs": ["design_03"]} if stage_type == "reversion" else {}
    s = sr.StageResult(stage_type=stage_type, config_id="c",
                       predictions=[_fake_pred(2.0, 1.0)],
                       confidences=[_fake_conf()], seed_indices=[0],
                       mutated_positions=mp, **kwargs)
    t = pt.PipelineInternalThresholds.default()
    with pytest.raises(ValueError, match="requires contamination_positions"):
        s.per_seed_verdicts(t, _truth())


@pytest.mark.logic
def test_steered_gating_without_contamination_positions_raises():
    """triggers_next_stage decides whether reversion runs at all.  Without the
    contamination set it would find nothing contaminated and never revert."""
    mp = PositionSet([5], chain="A", frame="prediction")
    s = _make_steered([_fake_pred(2.0, 1.0, contact_positions=[5])], mp)
    t = pt.PipelineInternalThresholds.default()
    with pytest.raises(ValueError, match="requires contamination_positions"):
        s.triggers_next_stage(t, _truth())


# ── Frames must agree before the sets are intersected ──────────────


@pytest.mark.logic
def test_a_frame_mismatch_is_no_data_for_that_seed_only():
    """The contacts come back in the prediction frame.  Mutated positions read
    from mutations.tsv are in the designed frame, and PositionSet refuses to
    intersect across frames.  One seed's frame problem must not take down the
    stage, so it reads no_data while its neighbour is still scored.
    """
    designed_mp = PositionSet([5], chain="A", frame="designed")
    contam = PositionSet([5], chain="A", frame="prediction")
    s = _make_steered(
        [_fake_pred(2.0, 1.0, contact_positions=[5]),
         _fake_pred(10.0, 1.0, contact_positions=[5])],
        designed_mp)
    t = pt.PipelineInternalThresholds.default()
    # The second seed never reaches the intersection: it fails the structural
    # check first and so is pose_collapses, not no_data.
    assert s.per_seed_verdicts(t, _truth(), contam) == ["no_data", "pose_collapses"]


# ── n_correctly_placed and n_contaminated ──────────────────────────


@pytest.mark.logic
def test_a_seed_with_no_metrics_is_not_counted_as_correctly_placed():
    """2 of 3 seeds have metrics and both are well placed.  The third returns
    None, which must not count either way — silently counting it would inflate
    the denominator of the CL-3 majority rule."""
    preds = [
        _fake_pred(2.0, 1.0),
        _fake_pred(2.0, 1.0),
        _fake_pred(None, None),
    ]
    s = _make_cold_start(preds)
    t = pt.PipelineInternalThresholds.default()
    assert s.n_correctly_placed(t, _truth()) == 2


@pytest.mark.logic
def test_a_cold_start_stage_reports_no_contamination():
    """Cold start runs the native sequence, so it has no mutated positions for
    a contaminating contact to land on.  The count is 0 by definition rather
    than by measurement."""
    contam = PositionSet([5], chain="A", frame="prediction")
    s = _make_cold_start([_fake_pred(2.0, 1.0, contact_positions=[5])] * 3)
    t = pt.PipelineInternalThresholds.default()
    assert s.n_contaminated(t, _truth(), contam) == 0


@pytest.mark.logic
def test_a_reversion_counts_new_contamination_as_its_contaminated_seeds():
    """The verdict word differs by stage — `contaminated` when steered,
    `new_contamination` after reversion — and n_contaminated has to count the
    right one.  2 of 3 seeds here contact a mutated, contaminating position."""
    mp = PositionSet([5, 7], chain="A", frame="prediction")
    contam = PositionSet([5, 7], chain="A", frame="prediction")
    preds = [
        _fake_pred(2.0, 1.0, contact_positions=[5]),
        _fake_pred(2.0, 1.0, contact_positions=[7]),
        _fake_pred(2.0, 1.0, contact_positions=[100]),
    ]
    s = _make_reversion(preds, mp)
    t = pt.PipelineInternalThresholds.default()
    assert s.n_contaminated(t, _truth(), contam) == 2


# ── canonical_prediction ───────────────────────────────────────────
#
# The per-stage representative.  No caller in this repository uses it
# yet, so these tests pin the selection rule before one does.


@pytest.mark.logic
def test_the_canonical_prediction_is_the_lowest_ra_eff_passing_seed():
    """ra_eff 4.0, 2.0, 3.0 and all three clean, so seed 1 wins on 2.0."""
    preds = [_fake_pred(4.0, 1.0), _fake_pred(2.0, 1.0), _fake_pred(3.0, 1.0)]
    s = _make_cold_start(preds)
    t = pt.PipelineInternalThresholds.default()
    assert s.canonical_prediction(t, _truth()) is preds[1]


@pytest.mark.logic
def test_a_passing_seed_beats_a_better_placed_failing_one():
    """Seed 0 is the better pose at ra_eff 1.0 but is contaminated.  Seed 1 is
    clean_steered at 3.0.  The representative has to be the clean one, or the
    stage would be represented by the result it was steered away from."""
    mp = PositionSet([5], chain="A", frame="prediction")
    contam = PositionSet([5], chain="A", frame="prediction")
    preds = [
        _fake_pred(1.0, 1.0, contact_positions=[5]),    # contaminated
        _fake_pred(3.0, 1.0, contact_positions=[100]),  # clean_steered
    ]
    s = _make_steered(preds, mp)
    t = pt.PipelineInternalThresholds.default()
    assert s.canonical_prediction(t, _truth(), contam) is preds[1]


@pytest.mark.logic
def test_with_nothing_passing_the_lowest_ra_eff_still_represents_the_stage():
    """All three collapse, at ra_eff 12.0, 8.0 and 20.0.  The stage still
    needs a representative, and the least-bad pose is seed 1."""
    mp = PositionSet([5], chain="A", frame="prediction")
    contam = PositionSet([5], chain="A", frame="prediction")
    preds = [_fake_pred(12.0, 1.0), _fake_pred(8.0, 1.0), _fake_pred(20.0, 1.0)]
    s = _make_steered(preds, mp)
    t = pt.PipelineInternalThresholds.default()
    assert s.canonical_prediction(t, _truth(), contam) is preds[1]


@pytest.mark.logic
def test_with_no_ra_eff_anywhere_the_first_seed_represents_the_stage():
    """Nothing is measurable, so there is nothing to rank by and the first
    prediction stands in.  Returning None instead would make the caller handle
    an absent representative for a stage that does have predictions."""
    preds = [_fake_pred(None, None), _fake_pred(None, None)]
    s = _make_cold_start(preds)
    t = pt.PipelineInternalThresholds.default()
    assert s.canonical_prediction(t, _truth()) is preds[0]


# ── Deployment: the module finds its own siblings ──────────────────


@pytest.mark.logic
def test_the_module_adds_its_own_directory_to_sys_path():
    """A Nextflow process invokes each bin/ script by absolute path with no
    PYTHONPATH, so stage_result.py has to put bin/ on sys.path itself before
    importing position_set and the two metric types.  Re-executing the file
    with bin/ removed proves the guard works rather than relying on a test
    runner having inserted bin/ already.
    """
    import importlib.util

    bin_dir = _REPO_ROOT / "bin"
    saved = list(sys.path)
    try:
        sys.path[:] = [p for p in sys.path if p != str(bin_dir)]
        spec = importlib.util.spec_from_file_location(
            "stage_result__fresh", bin_dir / "stage_result.py")
        module = importlib.util.module_from_spec(spec)
        # @dataclass resolves its annotations through sys.modules, so the
        # module has to be registered before it is executed.
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        assert str(bin_dir) in sys.path
        assert module.VALID_STAGE_TYPES == frozenset(
            ("cold_start", "steered", "reversion"))
    finally:
        sys.path[:] = saved
        sys.modules.pop("stage_result__fresh", None)
