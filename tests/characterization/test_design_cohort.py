"""Unit tests for bin/design_cohort.py (Phase 4 Tier 6).

DesignCohort is the top of the type hierarchy: every MPNN sequence in one
pipeline run, each carrying its NegativeSteeringRun. It answers the three
questions the campaign is reported in terms of — which tier each sequence
reached, how they rank against each other, and which survive into the
orthogonal stage.

It is also the cohort CSV emitter. `cross_summary_v2.py`, the script the
NEGSTEER_CROSS_SEQUENCE Nextflow process runs, is a thin CLI over
`emit_cross_summary_from_dirs`, so the directory walk here decides which
sequences appear in cross_sequence_summary.csv at all. A sequence silently
dropped by that walk is absent from every figure and every count in the
report, and nothing downstream can tell it was ever meant to be there. The
last sections drive the walk and assert the CONTENTS of the CSV it writes.
"""

from __future__ import annotations

import csv
import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_BIN = _REPO_ROOT / "bin"
_HERE = Path(__file__).resolve().parent
for _p in (str(_BIN), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import design_cohort as dc  # noqa: E402
import fixtures_cohort as fc  # noqa: E402
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


def _truth():
    return MagicMock()


def _make_run(tmp_path, mpnn_id, row_type="steered", n_clean=3, n_seeds=3,
              workdir=None):
    """Build a NegativeSteeringRun on the cold-start path with `n_clean`
    of `n_seeds` clean predictions.

    `workdir` defaults to `tmp_path`. The CSV emitter reads
    `workdir / 'passing_summary.csv'`, so a test about the emitter gives each
    run its own directory."""
    preds = []
    for i in range(n_seeds):
        if i < n_clean:
            preds.append(_fake_pred(2.0, 1.0))  # clean
        else:
            preds.append(_fake_pred(10.0, 1.0))  # fail_structural
    cs = sr.StageResult(
        stage_type="cold_start", config_id=mpnn_id,
        predictions=preds,
        confidences=[_fake_conf() for _ in preds],
        seed_indices=list(range(n_seeds)))
    return nsr.NegativeSteeringRun(
        mpnn_sequence_id=mpnn_id,
        workdir=tmp_path if workdir is None else workdir,
        row_type=row_type, cold_start=cs, truth=_truth())


# ── Construction ───────────────────────────────────────────────────


@pytest.mark.logic
class TestConstruction:
    def test_basic(self, tmp_path):
        cohort = dc.DesignCohort(runs=[
            _make_run(tmp_path, "design_03_seq_01"),
            _make_run(tmp_path, "design_03_seq_02"),
        ])
        assert len(cohort) == 2

    def test_duplicate_ids_raise(self, tmp_path):
        with pytest.raises(ValueError, match="duplicate"):
            dc.DesignCohort(runs=[
                _make_run(tmp_path, "x"),
                _make_run(tmp_path, "x"),
            ])

    def test_wrong_type_in_runs_raises(self):
        with pytest.raises(TypeError):
            dc.DesignCohort(runs=["not a NSR"])

    def test_iterator(self, tmp_path):
        cohort = dc.DesignCohort(runs=[
            _make_run(tmp_path, "a"), _make_run(tmp_path, "b"),
        ])
        ids = [r.mpnn_sequence_id for r in cohort]
        assert ids == ["a", "b"]


# ── Lookup ────────────────────────────────────────────────────────


@pytest.mark.logic
class TestLookup:
    def test_get_existing(self, tmp_path):
        a = _make_run(tmp_path, "a")
        cohort = dc.DesignCohort(runs=[a, _make_run(tmp_path, "b")])
        assert cohort.get("a") is a

    def test_get_missing_returns_none(self, tmp_path):
        cohort = dc.DesignCohort(runs=[_make_run(tmp_path, "a")])
        assert cohort.get("zzz") is None


# ── Row-type filtering ────────────────────────────────────────────


@pytest.mark.logic
class TestRowTypeFiltering:
    def test_steered_and_controls(self, tmp_path):
        cohort = dc.DesignCohort(runs=[
            _make_run(tmp_path, "s1", row_type="steered"),
            _make_run(tmp_path, "s2", row_type="steered"),
            _make_run(tmp_path, "cs", row_type="control_scrambled"),
            _make_run(tmp_path, "cp", row_type="control_polyA"),
        ])
        assert len(cohort.steered_runs()) == 2
        assert len(cohort.control_runs()) == 2
        assert cohort.n_steered() == 2
        assert cohort.n_controls() == 2


# ── Tier breakdown ────────────────────────────────────────────────


@pytest.mark.logic
class TestTierBreakdown:
    def test_full_breakdown(self, tmp_path):
        cohort = dc.DesignCohort(runs=[
            _make_run(tmp_path, "a", n_clean=3, n_seeds=3),  # tier A
            _make_run(tmp_path, "b", n_clean=2, n_seeds=3),  # tier B
            _make_run(tmp_path, "c", n_clean=1, n_seeds=3),  # tier C
            _make_run(tmp_path, "d", n_clean=0, n_seeds=3),  # tier none
        ])
        t = pt.PipelineInternalThresholds.default()
        breakdown = cohort.tier_breakdown(t)
        assert breakdown == {"A": 1, "B": 1, "C": 1, "none": 1}

    def test_excludes_controls(self, tmp_path):
        cohort = dc.DesignCohort(runs=[
            _make_run(tmp_path, "a", n_clean=3, row_type="steered"),
            _make_run(tmp_path, "cs", n_clean=3, row_type="control_scrambled"),
        ])
        t = pt.PipelineInternalThresholds.default()
        assert cohort.tier_breakdown(t) == {"A": 1, "B": 0, "C": 0, "none": 0}


# ── by_tier and survivors ─────────────────────────────────────────


@pytest.mark.logic
class TestByTierAndSurvivors:
    @pytest.fixture
    def cohort(self, tmp_path):
        return dc.DesignCohort(runs=[
            _make_run(tmp_path, "a", n_clean=3),  # A
            _make_run(tmp_path, "b", n_clean=2),  # B
            _make_run(tmp_path, "c", n_clean=1),  # C
            _make_run(tmp_path, "d", n_clean=0),  # none
        ])

    def test_by_tier(self, cohort):
        t = pt.PipelineInternalThresholds.default()
        a_runs = cohort.by_tier("A", t)
        assert len(a_runs) == 1
        assert a_runs[0].mpnn_sequence_id == "a"

    def test_survivors_are_abc(self, cohort):
        t = pt.PipelineInternalThresholds.default()
        survivors = cohort.survivors(t)
        ids = {r.mpnn_sequence_id for r in survivors}
        assert ids == {"a", "b", "c"}  # "d" (tier none) excluded


# ── Ranking ───────────────────────────────────────────────────────


@pytest.mark.logic
class TestRanking:
    def test_ranked_by_composite(self, tmp_path):
        cohort = dc.DesignCohort(runs=[
            _make_run(tmp_path, "tier_c", n_clean=1),  # C
            _make_run(tmp_path, "tier_a", n_clean=3),  # A
            _make_run(tmp_path, "tier_none", n_clean=0),  # none
            _make_run(tmp_path, "tier_b", n_clean=2),  # B
        ])
        t = pt.PipelineInternalThresholds.default()
        ranked = cohort.ranked_by_composite(t)
        # Order should be A, B, C, none
        assert [r.mpnn_sequence_id for r in ranked] == [
            "tier_a", "tier_b", "tier_c", "tier_none"
        ]


@pytest.mark.logic
class TestImmutability:
    def test_frozen(self, tmp_path):
        cohort = dc.DesignCohort(runs=[_make_run(tmp_path, "x")])
        with pytest.raises((AttributeError, Exception)):
            cohort.runs = ()  # type: ignore


# ── Cohort CSV from the runs in hand ─────────────────────────────────
#
# Every passing_summary row below is written by fixtures_cohort.passing_row,
# whose defaults put ra_eff at 2.0 Å and the true Jaccard at 0.80. The
# composite the emitter scores is jaccard - 0.05 x ra_eff, so each row's
# score is 0.80 - 0.05 x 2.0 = 0.70 and the ranking is decided by tier
# alone. n_pass against n_seeds = 3 is what sets the tier: 3/3 is A, 2/3 is
# B, 1/3 is C.

def _sequence_dir(runs: Path, name: str, *, n_pass: int, **kw) -> Path:
    """One runs/<mpnn_sequence>/ holding a one-row passing_summary.csv."""
    d = runs / name
    d.mkdir(parents=True, exist_ok=True)
    fc.write_passing_summary(d, [fc.passing_row(n_pass=n_pass, **kw)])
    return d


def _rows(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


@pytest.mark.logic
def test_the_cohort_csv_carries_one_row_per_run_with_its_tier(tmp_path):
    """The tier comes from each run's own passing_summary.csv, not from the
    cohort: 3 of 3 passing seeds is tier A and 2 of 3 is tier B, so the two
    rank 1 and 2 in that order. A row missing from this file is a sequence
    missing from every figure in the report."""
    runs = tmp_path / "runs"
    cohort = dc.DesignCohort(runs=[
        _make_run(tmp_path, "design_0_seq_0",
                  workdir=_sequence_dir(runs, "design_0_seq_0", n_pass=3)),
        _make_run(tmp_path, "design_0_seq_1",
                  workdir=_sequence_dir(runs, "design_0_seq_1", n_pass=2)),
    ])
    out = tmp_path / "cross_sequence_summary.csv"
    assert cohort.to_cross_summary_csv(out) == 0

    rows = _rows(out)
    assert [r["mpnn_sequence"] for r in rows] == ["design_0_seq_0",
                                                  "design_0_seq_1"]
    assert [r["cross_tier"] for r in rows] == ["A", "B"]
    assert [r["cross_rank_by_composite"] for r in rows] == ["1", "2"]
    # 0.80 - 0.05 x 2.0 = 0.70, the same for both, so tier decided the order.
    assert {r["cross_composite_score"] for r in rows} == {"0.700000"}


@pytest.mark.logic
def test_the_rows_are_named_by_sequence_id_not_by_directory(tmp_path):
    """The emitter passes each run's mpnn_sequence_id, not its workdir name.
    They differ in a published tree, where the workdir is a Nextflow hash
    directory, and the CSV is joined to the MPNN FASTA by sequence id."""
    runs = tmp_path / "runs"
    cohort = dc.DesignCohort(runs=[
        _make_run(tmp_path, "design_7_seq_2",
                  workdir=_sequence_dir(runs, "ab12cd34", n_pass=3)),
    ])
    out = tmp_path / "cross.csv"
    assert cohort.to_cross_summary_csv(out) == 0
    assert [r["mpnn_sequence"] for r in _rows(out)] == ["design_7_seq_2"]


# ── Cohort CSV from a runs directory ─────────────────────────────────
#
# This is the path cross_summary_v2.py takes, and so the path the pipeline
# takes. Every test here is about which sequences the walk admits.

@pytest.mark.logic
def test_every_subdirectory_of_a_runs_directory_becomes_a_row(tmp_path):
    runs = tmp_path / "runs"
    _sequence_dir(runs, "design_0_seq_0", n_pass=3)
    _sequence_dir(runs, "design_1_seq_0", n_pass=1)
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[runs], output_path=out) == 0
    rows = _rows(out)
    assert [r["mpnn_sequence"] for r in rows] == ["design_0_seq_0",
                                                  "design_1_seq_0"]
    assert [r["cross_tier"] for r in rows] == ["A", "C"]


@pytest.mark.logic
def test_a_file_beside_the_sequence_directories_is_not_a_sequence(tmp_path):
    """A published runs/ tree carries sibling files. Admitting one would add a
    row named after it whose passing_summary.csv can never be found."""
    runs = tmp_path / "runs"
    _sequence_dir(runs, "design_0_seq_0", n_pass=3)
    (runs / "README.txt").write_text("published by the negsteer stage\n")
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[runs], output_path=out) == 0
    assert [r["mpnn_sequence"] for r in _rows(out)] == ["design_0_seq_0"]


@pytest.mark.logic
def test_a_directory_name_that_is_not_a_sequence_name_is_skipped(tmp_path):
    """The name becomes a column value and a filename stem downstream, so the
    walk admits only [A-Za-z0-9_.+-]. A directory with a space in it is a
    staging artefact, not a sequence."""
    runs = tmp_path / "runs"
    _sequence_dir(runs, "design_0_seq_0", n_pass=3)
    _sequence_dir(runs, "copy of design_0_seq_0", n_pass=3)
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[runs], output_path=out) == 0
    assert [r["mpnn_sequence"] for r in _rows(out)] == ["design_0_seq_0"]


@pytest.mark.logic
def test_a_runs_directory_that_is_not_there_is_skipped_not_fatal(tmp_path):
    """Nextflow stages the per-cycle runs directories as a list and a cycle
    that produced nothing leaves no directory. Aborting would lose the
    cycles that did run."""
    runs = tmp_path / "runs"
    _sequence_dir(runs, "design_0_seq_0", n_pass=3)
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[tmp_path / "cycle_1_runs", runs], output_path=out) == 0
    assert [r["mpnn_sequence"] for r in _rows(out)] == ["design_0_seq_0"]


@pytest.mark.logic
def test_one_sequence_staged_twice_is_emitted_once(tmp_path):
    """--passing-summary-dir is repeated once per cycle, and a sequence that
    ran in two cycles appears under both. Two rows for one sequence would
    double it in every tier count and give it two ranks."""
    first = tmp_path / "a" / "runs"
    second = tmp_path / "b" / "runs"
    _sequence_dir(first, "design_0_seq_0", n_pass=3)
    _sequence_dir(second, "design_0_seq_0", n_pass=1)
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[first, second], output_path=out) == 0
    rows = _rows(out)
    assert [r["mpnn_sequence"] for r in rows] == ["design_0_seq_0"]
    # The first directory staged wins, so the tier is A rather than C.
    assert rows[0]["cross_tier"] == "A"


@pytest.mark.logic
def test_an_explicitly_named_summary_wins_over_the_walk(tmp_path):
    """--passing-summary NAME=PATH exists to point one sequence at a file
    outside the walked tree. The walk must not then overwrite it with the
    directory of the same name."""
    runs = tmp_path / "runs"
    _sequence_dir(runs, "design_0_seq_0", n_pass=1)      # tier C in the walk
    explicit = _sequence_dir(tmp_path / "elsewhere", "design_0_seq_0",
                             n_pass=3)                   # tier A explicitly
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[runs], output_path=out,
        explicit_pairs=[("design_0_seq_0",
                         explicit / "passing_summary.csv")]) == 0
    rows = _rows(out)
    assert len(rows) == 1
    assert rows[0]["cross_tier"] == "A"


@pytest.mark.logic
def test_the_same_sequence_named_twice_explicitly_is_emitted_once(tmp_path):
    """The CLI appends --passing-summary, so a Nextflow channel that emits a
    sequence twice names it twice. The first naming wins."""
    a = _sequence_dir(tmp_path / "a", "design_0_seq_0", n_pass=3)
    b = _sequence_dir(tmp_path / "b", "design_0_seq_0", n_pass=1)
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[], output_path=out,
        explicit_pairs=[("design_0_seq_0", a / "passing_summary.csv"),
                        ("design_0_seq_0", b / "passing_summary.csv")]) == 0
    rows = _rows(out)
    assert len(rows) == 1
    assert rows[0]["cross_tier"] == "A"


@pytest.mark.logic
def test_nothing_to_emit_exits_two_and_writes_no_csv(tmp_path):
    """An empty cohort CSV would read as a campaign in which every sequence
    failed, which is indistinguishable from a campaign whose runs were never
    staged. Exit 2 says the difference."""
    empty = tmp_path / "runs"
    empty.mkdir()
    out = tmp_path / "cross.csv"
    assert dc.DesignCohort.emit_cross_summary_from_dirs(
        runs_dirs=[empty], output_path=out) == 2
    assert not out.exists()


# ── Walking a runs directory into typed runs ─────────────────────────

@pytest.mark.logic
def test_a_file_in_the_runs_directory_is_not_walked_as_a_workdir(tmp_path):
    """from_runs_directory walks the same published tree, which holds
    cross_sequence_summary.csv and the cohort figures alongside the
    per-sequence directories."""
    runs = tmp_path / "runs"
    runs.mkdir()
    (runs / "cross_sequence_summary.csv").write_text("mpnn_sequence\n")
    assert len(dc.DesignCohort.from_runs_directory(runs)) == 0


# ── Reporting ────────────────────────────────────────────────────────

@pytest.mark.logic
def test_the_repr_states_the_steered_and_control_counts(tmp_path):
    """The repr is what a stage log shows, and the steered/control split is
    the one number that says whether the controls were staged at all."""
    cohort = dc.DesignCohort(runs=[
        _make_run(tmp_path, "design_0_seq_0"),
        _make_run(tmp_path, "design_0_seq_1"),
        _make_run(tmp_path, "input_control_polyA", row_type="control_polyA"),
    ])
    assert repr(cohort) == "DesignCohort(n_steered=2, n_controls=1)"


# ── Importability from a bare interpreter ────────────────────────────

@pytest.mark.logic
def test_the_module_puts_its_own_directory_on_the_path():
    """Every Nextflow process invokes bin/ scripts by absolute path with no
    PYTHONPATH set, so design_cohort has to put its own directory on sys.path
    before importing its siblings. The test suite already has bin/ there, so
    the guard is only observable by executing the file again with bin/ taken
    away, which is as close to the Nextflow invocation as an in-process test
    gets."""
    bin_dir = str(_BIN)
    spec = importlib.util.spec_from_file_location(
        "design_cohort_without_bin_on_the_path", _BIN / "design_cohort.py")
    module = importlib.util.module_from_spec(spec)
    saved = list(sys.path)
    # @dataclass resolves its annotations through sys.modules, so the module
    # has to be registered there before it is executed.
    sys.modules[spec.name] = module
    try:
        sys.path[:] = [p for p in saved if p != bin_dir]
        spec.loader.exec_module(module)
        assert bin_dir in sys.path
        assert module.NegativeSteeringRun is nsr.NegativeSteeringRun
    finally:
        sys.path[:] = saved
        sys.modules.pop(spec.name, None)
