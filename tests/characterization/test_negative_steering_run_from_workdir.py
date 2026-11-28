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

The engine itself (`negsteer run`, a Singularity container) is not part of
this repository and is never invoked here.  ``from_workdir`` only reads what
the engine left on disk, so the stand-in is the output tree, built by
``fixtures_runtree.negsteer_workdir``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_BIN = Path(__file__).resolve().parents[2] / "bin"
_HERE = Path(__file__).resolve().parent
for _p in (str(_BIN), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from fixtures_runtree import negsteer_workdir  # noqa: E402


@pytest.fixture
def synthetic_workdir(tmp_path: Path) -> Path:
    """Build a minimal workdir with 3 cold-start seeds, 2 steered designs
    (2 seeds each), and 1 reversion (2 seeds)."""
    return negsteer_workdir(
        tmp_path / "design_42_seq_0",
        cold_start=(0, 1, 2),
        steered={
            "design_00_s0": {"mutations": [(5, "K", "R")]},
            "design_00_s1": {"mutations": [(5, "K", "R")]},
            "design_01_s0": {"mutations": [(7, "I", "L")]},
            "design_01_s1": {"mutations": [(7, "I", "L")]},
        },
        reversions={
            "rev_design_00_s0_s0": {"mutations": [(5, "R", "K")]},
            "rev_design_00_s0_s1": {"mutations": [(5, "R", "K")]},
        },
        row_type="steered",
        runtime_sec=1234.5,
    )


@pytest.mark.logic
def test_from_workdir_missing_returns_none(tmp_path: Path):
    """Workdir without cycle_0 should return None."""
    from negative_steering_run import NegativeSteeringRun
    empty = tmp_path / "empty"
    empty.mkdir()
    run = NegativeSteeringRun.from_workdir(
        workdir=empty, mpnn_sequence_id="x", num_seeds=3,
    )
    assert run is None


@pytest.mark.logic
def test_from_workdir_no_pdbs_returns_none(tmp_path: Path):
    """cycle_0 dir present but no Boltz outputs → no cold-start data → None."""
    from negative_steering_run import NegativeSteeringRun
    wd = tmp_path / "design_0_seq_0"
    (wd / "cycle_0" / "initial").mkdir(parents=True)
    run = NegativeSteeringRun.from_workdir(
        workdir=wd, mpnn_sequence_id="design_0_seq_0", num_seeds=3,
    )
    assert run is None


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
def test_from_runs_directory_collects_workdirs(tmp_path: Path):
    """DesignCohort.from_runs_directory walks subdirs and skips ones
    that aren't recognisable as workdirs."""
    from design_cohort import DesignCohort

    runs = tmp_path / "runs"
    # Two valid workdirs:
    for name in ("design_0_seq_0", "design_1_seq_0"):
        negsteer_workdir(runs / name, cold_start=(0, 1, 2))
    # One bogus dir without cycle_0:
    (runs / "bogus_dir").mkdir(parents=True)

    cohort = DesignCohort.from_runs_directory(runs, num_seeds=3)
    assert len(cohort) == 2
    ids = {r.mpnn_sequence_id for r in cohort}
    assert ids == {"design_0_seq_0", "design_1_seq_0"}


# ── Finding the Boltz outputs inside one prediction directory ────────
#
# A seed counts only when its PDB and its confidence sidecar are both
# present.  A half-written directory is what a GPU job killed mid-write
# leaves behind, and taking it as a seed would put a prediction in the
# stage with no confidence to pair it with.


def _run(workdir: Path, **kwargs):
    from negative_steering_run import NegativeSteeringRun
    kwargs.setdefault("mpnn_sequence_id", "design_42_seq_0")
    kwargs.setdefault("num_seeds", 3)
    return NegativeSteeringRun.from_workdir(workdir=workdir, **kwargs)


@pytest.mark.logic
def test_a_pdb_written_without_the_pdb_subdir_is_still_found(tmp_path: Path):
    """Some engine versions write input_model_0.pdb beside the sidecar rather
    than under pdb/.  Both layouts must yield the same three seeds, or a run
    from an older engine version reads as having no cold start at all."""
    wd = negsteer_workdir(
        tmp_path / "flat",
        cold_start={0: {"flat_pdb": True}, 1: {"flat_pdb": True},
                    2: {"flat_pdb": True}},
    )
    run = _run(wd)
    assert run is not None
    assert run.cold_start.seed_indices == (0, 1, 2)


@pytest.mark.logic
def test_a_seed_with_no_pdb_in_either_layout_is_dropped(tmp_path: Path):
    """Seed 1 has its confidence sidecar and no PDB under pdb/ or beside it,
    so the two surviving seeds are 0 and 2."""
    wd = negsteer_workdir(
        tmp_path / "no_pdb",
        cold_start={0: {}, 1: {"write_pdb": False}, 2: {}},
    )
    run = _run(wd)
    assert run is not None
    assert run.cold_start.seed_indices == (0, 2)


@pytest.mark.logic
def test_a_seed_with_no_confidence_sidecar_is_dropped(tmp_path: Path):
    """Seed 0 has a PDB and no confidence JSON, so only seeds 1 and 2 remain
    and the stage reports 2 seeds rather than 3."""
    wd = negsteer_workdir(
        tmp_path / "no_conf",
        cold_start={0: {"write_confidence": False}, 1: {}, 2: {}},
    )
    run = _run(wd)
    assert run is not None
    assert run.cold_start.seed_indices == (1, 2)
    assert run.num_seeds == 2


@pytest.mark.logic
def test_require_pdb_false_drops_the_same_seeds_as_require_pdb_true(tmp_path: Path):
    """``require_pdb`` is inert, which from_workdir's own docstring says: no
    ProteinStructurePrediction can be built without both files, so the seed is
    skipped either way.  The flag is kept for callers that pass it."""
    layout = {0: {}, 1: {"write_pdb": False}, 2: {}}
    strict = _run(negsteer_workdir(tmp_path / "strict", cold_start=layout),
                  require_pdb=True)
    lax = _run(negsteer_workdir(tmp_path / "lax", cold_start=layout),
               require_pdb=False)
    assert strict is not None and lax is not None
    assert strict.cold_start.seed_indices == lax.cold_start.seed_indices == (0, 2)


# ── mutations.tsv ────────────────────────────────────────────────────


@pytest.mark.logic
def test_mutation_positions_are_read_as_the_1_based_designed_frame(tmp_path: Path):
    """mutations.tsv's pos1 column is already 1-based on the designed-frame
    receptor, so it is carried across unchanged.  It is not the engine
    plan.json's `mutated_positions`, which are 0-based and would need +1."""
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={"design_00_s0": {"mutations": [(5, "K", "A"), (12, "D", "A")]}},
    )
    run = _run(wd)
    assert run is not None
    mp = run.steered["design_00"].mutated_positions
    assert mp.frame == "designed"
    assert mp.as_sorted_list() == [5, 12]


@pytest.mark.logic
def test_a_steered_design_without_mutations_tsv_keeps_an_empty_position_set(
    tmp_path: Path,
):
    """StageResult requires mutated_positions on a steered stage.  Dropping
    the stage instead would lose the predictions, so the reader synthesises an
    empty designed-frame set."""
    wd = negsteer_workdir(tmp_path / "wd", steered={"design_00_s0": {}})
    run = _run(wd)
    assert run is not None
    mp = run.steered["design_00"].mutated_positions
    assert mp.is_empty()
    assert mp.frame == "designed"


@pytest.mark.logic
def test_a_reversion_without_mutations_tsv_keeps_an_empty_position_set(
    tmp_path: Path,
):
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={"design_00_s0": {"mutations": [(5, "K", "A")]}},
        reversions={"rev_design_00_s0_s0": {}},
    )
    run = _run(wd)
    assert run is not None
    assert run.reversion["rev_design_00_s0"].mutated_positions.is_empty()


@pytest.mark.logic
def test_an_unparseable_mutation_row_is_skipped_and_the_rest_are_kept(
    tmp_path: Path,
):
    """One bad row must not cost the whole design its mutation set.  Three
    rows go in, two of them parse, so the answer is [7, 9]."""
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={"design_00_s0": {
            "mutations": "pos1\twt\tmut\n7\tK\tA\nnot-a-number\tD\tA\n9\tE\tA\n",
        }},
    )
    run = _run(wd)
    assert run is not None
    assert run.steered["design_00"].mutated_positions.as_sorted_list() == [7, 9]


# ── Subdirectory names the reader must ignore ────────────────────────


@pytest.mark.logic
def test_files_and_unrecognised_names_under_steered_are_ignored(tmp_path: Path):
    """Only `design_NN_sS` directories are stage data.  A log file or a
    scratch directory alongside them must not become a steered design."""
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={"design_00_s0": {"mutations": [(5, "K", "A")]}},
    )
    steered_root = wd / "cycle_0" / "steered"
    (steered_root / "run.log").write_text("engine chatter\n")
    (steered_root / "tmp_scratch").mkdir()
    run = _run(wd)
    assert run is not None
    assert set(run.steered) == {"design_00"}


@pytest.mark.logic
def test_files_and_unrecognised_names_under_reversions_are_ignored(tmp_path: Path):
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={"design_00_s0": {"mutations": [(5, "K", "A")]}},
        reversions={"rev_design_00_s0_s0": {"mutations": [(5, "A", "K")]}},
    )
    rev_root = wd / "cycle_0" / "reversions"
    (rev_root / "reversion.log").write_text("engine chatter\n")
    (rev_root / "rev_design_00").mkdir()        # missing both seed suffixes
    run = _run(wd)
    assert run is not None
    assert set(run.reversion) == {"rev_design_00_s0"}


# ── A stage directory with no usable seed at all ─────────────────────


@pytest.mark.logic
def test_a_steered_design_whose_seeds_are_all_half_written_is_dropped(
    tmp_path: Path,
):
    """design_01's two seeds each lost a file, so it contributes no stage.
    design_00 survives, which is what distinguishes "dropped one design" from
    "failed to read the workdir"."""
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={
            "design_00_s0": {"mutations": [(5, "K", "A")]},
            "design_01_s0": {"write_pdb": False},
            "design_01_s1": {"write_confidence": False},
        },
    )
    run = _run(wd)
    assert run is not None
    assert set(run.steered) == {"design_00"}


@pytest.mark.logic
def test_a_reversion_whose_seeds_are_all_half_written_is_dropped(tmp_path: Path):
    wd = negsteer_workdir(
        tmp_path / "wd",
        steered={"design_00_s0": {"mutations": [(5, "K", "A")]}},
        reversions={
            "rev_design_00_s0_s0": {"write_pdb": False},
            "rev_design_00_s0_s1": {"write_confidence": False},
        },
    )
    run = _run(wd)
    assert run is not None
    assert run.reversion == {}


# ── The sidecars ─────────────────────────────────────────────────────


@pytest.mark.logic
def test_an_unparseable_runtime_sidecar_reads_as_no_runtime(tmp_path: Path):
    """run_one_runtime_sec.txt is diagnostic only.  A truncated or garbled
    value must not fail the whole read, so it becomes None."""
    wd = negsteer_workdir(tmp_path / "wd", runtime_sec="runtime: unknown")
    run = _run(wd)
    assert run is not None
    assert run.run_one_runtime_sec is None


@pytest.mark.logic
def test_an_absent_row_type_sidecar_defaults_to_steered(tmp_path: Path):
    """Only the negative-control rows carry a row_type.txt saying otherwise,
    so the absent case is an ordinary steered row."""
    wd = negsteer_workdir(tmp_path / "wd", row_type=None)
    run = _run(wd)
    assert run is not None
    assert run.row_type == "steered"


@pytest.mark.logic
def test_an_explicit_row_type_argument_overrides_the_sidecar(tmp_path: Path):
    wd = negsteer_workdir(tmp_path / "wd", row_type="steered")
    run = _run(wd, row_type="control_polyA")
    assert run is not None
    assert run.row_type == "control_polyA"


# ── Deployment: the module finds its own siblings ────────────────────


@pytest.mark.logic
def test_the_module_adds_its_own_directory_to_sys_path(tmp_path: Path):
    """A Nextflow process invokes each bin/ script by absolute path with no
    PYTHONPATH, so negative_steering_run.py has to put bin/ on sys.path
    itself before importing stage_result and friends.  Re-executing the file
    with bin/ removed proves the guard works, rather than relying on a test
    runner having inserted bin/ already.
    """
    import importlib.util

    saved = list(sys.path)
    try:
        sys.path[:] = [p for p in sys.path if p != str(_BIN)]
        spec = importlib.util.spec_from_file_location(
            "negative_steering_run__fresh", _BIN / "negative_steering_run.py")
        module = importlib.util.module_from_spec(spec)
        # @dataclass resolves its annotations through sys.modules, so the
        # module has to be registered before it is executed.
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        assert str(_BIN) in sys.path
        assert module.VALID_ROW_TYPES == frozenset(
            ("steered", "control_scrambled", "control_polyA"))
    finally:
        sys.path[:] = saved
        sys.modules.pop("negative_steering_run__fresh", None)
