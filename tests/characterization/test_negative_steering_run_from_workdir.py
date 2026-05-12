"""
tests/characterization/test_negative_steering_run_from_workdir.py
------------------------------------------------------------------
Tests for ``NegativeSteeringRun.from_workdir`` — the deep-form factory
that walks a per-MPNN workdir, discovers cold_start / steered / reversion
prediction directories, and constructs the typed Phase 4 scaffolding.

Most tests synthesize a minimal workdir layout in tmp_path with placeholder
PDB + confidence files.  Tests requiring gemmi parsing of structure data
are skipped if gemmi isn't installed (matching the convention in
test_protein_structure_prediction.py).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_BIN = Path(__file__).resolve().parents[2] / "bin"
if str(_BIN) not in sys.path:
    sys.path.insert(0, str(_BIN))


_MINIMAL_PDB = """ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.560  13.207  10.000  1.00  0.00           C
ATOM      3  C   ALA A   1      13.000  14.661  10.000  1.00  0.00           C
ATOM      4  N   GLY B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM      5  CA  GLY B   1      21.000  20.000  20.000  1.00  0.00           C
ATOM      6  C   GLY B   1      21.500  21.000  20.000  1.00  0.00           C
END
"""


def _make_pred_dir(parent: Path, mutations=None) -> Path:
    """Create a minimal Boltz-output-shaped prediction directory.

    Writes:
      <parent>/boltz_results_input/predictions/input/pdb/input_model_0.pdb
      <parent>/boltz_results_input/predictions/input/confidence_input_model_0.json
      <parent>/mutations.tsv  (if mutations provided)
    """
    pred_root = parent / "boltz_results_input" / "predictions" / "input"
    (pred_root / "pdb").mkdir(parents=True, exist_ok=True)
    (pred_root / "pdb" / "input_model_0.pdb").write_text(_MINIMAL_PDB)
    conf = {
        "confidence_score": 0.6,
        "ptm": 0.5,
        "iptm": 0.4,
        "complex_plddt": 0.7,
    }
    with open(pred_root / "confidence_input_model_0.json", "w") as f:
        json.dump(conf, f)
    if mutations is not None:
        lines = ["pos1\twt\tmut\n"]
        for pos, wt, mut in mutations:
            lines.append(f"{pos}\t{wt}\t{mut}\n")
        (parent / "mutations.tsv").write_text("".join(lines))
    return parent


@pytest.fixture
def synthetic_workdir(tmp_path: Path) -> Path:
    """Build a minimal workdir with 3 cold-start seeds, 2 steered designs
    (2 seeds each), and 1 reversion (2 seeds)."""
    wd = tmp_path / "design_42_seq_0"
    cycle = wd / "cycle_0"
    cycle.mkdir(parents=True)

    # Cold-start: initial / initial_s1 / initial_s2
    _make_pred_dir(cycle / "initial")
    _make_pred_dir(cycle / "initial_s1")
    _make_pred_dir(cycle / "initial_s2")

    # Steered: design_00_s0 + design_00_s1 ; design_01_s0 + design_01_s1
    steered_root = cycle / "steered"
    _make_pred_dir(steered_root / "design_00_s0", mutations=[(5, "K", "R")])
    _make_pred_dir(steered_root / "design_00_s1", mutations=[(5, "K", "R")])
    _make_pred_dir(steered_root / "design_01_s0", mutations=[(7, "I", "L")])
    _make_pred_dir(steered_root / "design_01_s1", mutations=[(7, "I", "L")])

    # Reversion: rev_design_00_s0_s0 + rev_design_00_s0_s1
    rev_root = cycle / "reversions"
    _make_pred_dir(rev_root / "rev_design_00_s0_s0", mutations=[(5, "R", "K")])
    _make_pred_dir(rev_root / "rev_design_00_s0_s1", mutations=[(5, "R", "K")])

    (wd / "row_type.txt").write_text("steered")
    (wd / "run_one_runtime_sec.txt").write_text("1234.5")
    return wd


@pytest.mark.local_unit
def test_from_workdir_missing_returns_none(tmp_path: Path):
    """Workdir without cycle_0 should return None."""
    from negative_steering_run import NegativeSteeringRun
    empty = tmp_path / "empty"
    empty.mkdir()
    run = NegativeSteeringRun.from_workdir(
        workdir=empty, mpnn_sequence_id="x", num_seeds=3,
    )
    assert run is None


@pytest.mark.local_unit
def test_from_workdir_no_pdbs_returns_none(tmp_path: Path):
    """cycle_0 dir present but no Boltz outputs → no cold-start data → None."""
    from negative_steering_run import NegativeSteeringRun
    wd = tmp_path / "design_0_seq_0"
    (wd / "cycle_0" / "initial").mkdir(parents=True)
    run = NegativeSteeringRun.from_workdir(
        workdir=wd, mpnn_sequence_id="design_0_seq_0", num_seeds=3,
    )
    assert run is None


@pytest.mark.local_unit
def test_from_workdir_hydrates_cold_start(synthetic_workdir: Path):
    from negative_steering_run import NegativeSteeringRun
    run = NegativeSteeringRun.from_workdir(
        workdir=synthetic_workdir,
        mpnn_sequence_id="design_42_seq_0",
        num_seeds=3,
    )
    assert run is not None
    assert run.mpnn_sequence_id == "design_42_seq_0"
    assert run.row_type == "steered"
    assert run.run_one_runtime_sec == pytest.approx(1234.5)
    # Cold-start has 3 seeds in our synthetic layout
    assert run.cold_start.stage_type == "cold_start"
    assert len(run.cold_start.predictions) == 3
    assert run.cold_start.seed_indices == (0, 1, 2)


@pytest.mark.local_unit
def test_from_workdir_hydrates_steered(synthetic_workdir: Path):
    from negative_steering_run import NegativeSteeringRun
    run = NegativeSteeringRun.from_workdir(
        workdir=synthetic_workdir,
        mpnn_sequence_id="design_42_seq_0",
        num_seeds=3,
    )
    assert run is not None
    assert set(run.steered.keys()) == {"design_00", "design_01"}
    d00 = run.steered["design_00"]
    assert d00.stage_type == "steered"
    assert d00.seed_indices == (0, 1)
    assert d00.mutated_positions is not None
    assert set(d00.mutated_positions) == {5}
    d01 = run.steered["design_01"]
    assert set(d01.mutated_positions) == {7}


@pytest.mark.local_unit
def test_from_workdir_hydrates_reversion(synthetic_workdir: Path):
    from negative_steering_run import NegativeSteeringRun
    run = NegativeSteeringRun.from_workdir(
        workdir=synthetic_workdir,
        mpnn_sequence_id="design_42_seq_0",
        num_seeds=3,
    )
    assert run is not None
    assert len(run.reversion) == 1
    rev_id, rev_stage = next(iter(run.reversion.items()))
    assert rev_stage.stage_type == "reversion"
    assert rev_stage.applies_to_designs == ("design_00",)
    assert rev_stage.seed_indices == (0, 1)


@pytest.mark.local_unit
def test_from_runs_directory_collects_workdirs(tmp_path: Path):
    """DesignCohort.from_runs_directory walks subdirs and skips ones
    that aren't recognisable as workdirs."""
    from design_cohort import DesignCohort

    runs = tmp_path / "runs"
    # Two valid workdirs:
    for name in ("design_0_seq_0", "design_1_seq_0"):
        wd = runs / name
        cycle = wd / "cycle_0"
        cycle.mkdir(parents=True)
        for sub in ("initial", "initial_s1", "initial_s2"):
            _make_pred_dir(cycle / sub)
    # One bogus dir without cycle_0:
    (runs / "bogus_dir").mkdir(parents=True)

    cohort = DesignCohort.from_runs_directory(runs, num_seeds=3)
    assert len(cohort) == 2
    ids = {r.mpnn_sequence_id for r in cohort}
    assert ids == {"design_0_seq_0", "design_1_seq_0"}
