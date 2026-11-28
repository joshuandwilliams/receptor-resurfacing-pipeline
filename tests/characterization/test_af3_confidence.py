"""Cover the AlphaFold3 confidence reader and its pose metric.

The receptor-aligned effector RMSD is the number that decides whether a pose
is correct. It aligns on the receptor and then measures how far the effector
sits from where the reference puts it, so a rigid move of the whole complex
must score zero and only a change in the relative placement may score.

The rotation convention is where this fails silently. `R @ v` and `v @ R` are
both valid-looking products, and the wrong one throws the effector tens of
angstroms off while still returning a float.

Residues are paired by POSITION, not by residue number, because the sources
number differently: AF3 emits per-chain 1-based numbering while RFdiffusion,
Boltz and the cleaned reference PDBs number continuously across chains.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import af3_confidence as af3  # noqa: E402
from fixtures_runtree import ca_line  # noqa: E402

gemmi = pytest.importorskip("gemmi")


def _write(path: Path, rows) -> Path:
    """rows: (chain, resnum, x, y, z)."""
    lines = [ca_line(i + 1, c, rn, x, y, z) for i, (c, rn, x, y, z) in enumerate(rows)]
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _complex_rows(n_rec=8, n_eff=5, dz=4.0, shift=(0.0, 0.0, 0.0),
                  eff_shift=(0.0, 0.0, 0.0), rec_start=1, eff_start=1):
    sx, sy, sz = shift
    ex, ey, ez = eff_shift
    rows = [("A", rec_start + i, i * 3.8 + sx, sy, sz) for i in range(n_rec)]
    rows += [("B", eff_start + j, j * 3.8 + sx + ex, sy + ey, dz + sz + ez)
             for j in range(n_eff)]
    return rows


# ── Chain coordinate reading ─────────────────────────────────────────

@pytest.mark.logic
def test_ca_coords_come_back_in_file_order(tmp_path):
    """Pairing is positional, so the order the reader returns is the frame
    every later comparison happens in."""
    p = _write(tmp_path / "s.pdb", _complex_rows(n_rec=4, n_eff=2))
    st = gemmi.read_structure(str(p))
    coords = af3._chain_ca_coords(st, "A")
    assert len(coords) == 4
    assert coords[2][0] == pytest.approx(2 * 3.8)


@pytest.mark.logic
def test_only_the_requested_chain_is_returned(tmp_path):
    p = _write(tmp_path / "s.pdb", _complex_rows(n_rec=8, n_eff=5))
    st = gemmi.read_structure(str(p))
    assert len(af3._chain_ca_coords(st, "A")) == 8
    assert len(af3._chain_ca_coords(st, "B")) == 5


@pytest.mark.logic
def test_an_absent_chain_gives_no_coordinates(tmp_path):
    p = _write(tmp_path / "s.pdb", _complex_rows())
    st = gemmi.read_structure(str(p))
    assert af3._chain_ca_coords(st, "Z") == []


@pytest.mark.logic
def test_one_coordinate_per_residue_even_with_several_atoms(tmp_path):
    """The reader breaks after the first CA in a residue. Without that an
    altloc would contribute twice and shift every later pairing."""
    p = tmp_path / "s.pdb"
    p.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CB  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  CA  ALA A   2       4.800   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    st = gemmi.read_structure(str(p))
    assert len(af3._chain_ca_coords(st, "A")) == 2


# ── The pose metric ──────────────────────────────────────────────────

@pytest.mark.logic
def test_an_identical_complex_scores_zero(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb", _complex_rows())
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_moving_the_whole_complex_still_scores_zero(tmp_path):
    """This is the test that catches a transposed rotation. Aligning the
    receptors undoes a rigid move of everything exactly, so the effector RMSD
    must be zero. The wrong matrix product throws it tens of angstroms off
    while still returning a float.
    """
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb",
                  _complex_rows(shift=(25.0, -13.0, 7.5)))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
@pytest.mark.parametrize("displacement", [1.0, 3.0, 7.5])
def test_moving_only_the_effector_scores_that_displacement(tmp_path, displacement):
    """Every effector atom moves by the same vector, so the RMSD is the
    length of that vector exactly."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb",
                  _complex_rows(eff_shift=(0.0, displacement, 0.0)))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(displacement, abs=1e-5)


@pytest.mark.logic
def test_displacement_is_the_vector_length_not_a_single_axis(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb", _complex_rows(eff_shift=(3.0, 4.0, 0.0)))
    rmsd, _err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert rmsd == pytest.approx(math.hypot(3.0, 4.0), abs=1e-5)


@pytest.mark.logic
def test_residue_numbering_does_not_affect_the_result(tmp_path):
    """AF3 numbers each chain from 1 and the reference numbers continuously.
    Pairing by residue number instead of position would pair nothing here, or
    pair the wrong residues."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows(rec_start=1, eff_start=9))
    pred = _write(tmp_path / "pred.pdb", _complex_rows(rec_start=1, eff_start=1))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_a_missing_chain_returns_an_error_tag_not_an_exception(tmp_path):
    """A garbage prediction must produce a record carrying an error rather
    than taking the stage down, because one bad design must not discard the
    rest of the batch."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb",
                  [("A", i + 1, i * 3.8, 0.0, 0.0) for i in range(8)])
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert rmsd is None
    assert isinstance(err, str) and err


@pytest.mark.logic
def test_an_unreadable_file_returns_an_error_tag(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    bad = tmp_path / "bad.cif"
    bad.write_text("this is not a structure\n")
    rmsd, err = af3._receptor_aligned_effector_rmsd(bad, ref, "A", "B")
    assert rmsd is None
    assert isinstance(err, str) and err


@pytest.mark.logic
def test_chains_of_different_lengths_do_not_pair_silently(tmp_path):
    """Truncating to the shorter chain would compare residue i of one against
    residue i of the other and report a number for two different molecules."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows(n_rec=8, n_eff=5))
    pred = _write(tmp_path / "pred.pdb", _complex_rows(n_rec=8, n_eff=3))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert (rmsd is None and err) or rmsd == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
def test_the_real_reference_structure_scores_zero_against_itself(tmp_path):
    """One pass over a committed structure rather than a built one, so the
    reader is exercised on a real file with real B-factors and headers."""
    pdbs = sorted(
        p for p in (REPO / "tests").rglob("*.pdb")
        if "data" in p.parts and p.stat().st_size > 2000)
    if not pdbs:
        pytest.skip("no committed multi-chain PDB available")
    src = pdbs[0]
    st = gemmi.read_structure(str(src))
    chains = [ch.name for model in st for ch in model]
    if len(set(chains)) < 2:
        pytest.skip("committed PDB is single-chain")
    a, b = chains[0], chains[1]
    rmsd, err = af3._receptor_aligned_effector_rmsd(src, src, a, b)
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-6)


# ── The aggregate's walk over an AF3 output tree ─────────────────────
#
# Everything above pins the pose metric on one pair of structures. AF3 emits
# one prediction per seed × diffusion sample, and the aggregate walks that set
# and reduces it to the six columns the orthogonal-metrics CSV carries. The
# block below drives the walk on trees where every number is known.


def _cif(path: Path, rows) -> Path:
    """Write `rows` as a real mmCIF at `path`.

    The walk globs for `.cif`, and gemmi chooses its parser from the
    extension, so a `.cif` file holding PDB text raises instead of loading.
    Going through gemmi rather than hand-writing the loop keeps the fixture
    honest about that.
    """
    lines = [ca_line(i + 1, c, rn, x, y, z)
             for i, (c, rn, x, y, z) in enumerate(rows)]
    st = gemmi.read_pdb_string("\n".join(lines) + "\nEND\n")
    st.setup_entities()
    path.parent.mkdir(parents=True, exist_ok=True)
    st.make_mmcif_document().write_file(str(path))
    return path


def _confidence(path: Path, iptm) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"iptm": iptm}))
    return path


class _Thresholds:
    """Stands in for pipeline_thresholds.PipelineThresholds, which the two
    threshold methods take only for this one attribute."""

    def __init__(self, ra_max: float):
        self.orthogonal_af3_ra_max = ra_max


def _tree(tmp_path):
    """One AF3 output directory in the canonical seed-N_sample-M layout.

    The reference is the unshifted complex. Sample 0 is identical to it, so
    its ra_eff is 0.0. Sample 1 has the effector moved 3 A in y with the
    receptor untouched, so its ra_eff is exactly 3.0. iPTMs are 0.80 and
    0.60, giving best 0.80 and mean 0.70.
    """
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    out = tmp_path / "af3_out"
    _cif(out / "seed-1_sample-0" / "model.cif", _complex_rows())
    _confidence(out / "seed-1_sample-0" / "summary_confidences.json", 0.80)
    _cif(out / "seed-1_sample-1" / "model.cif",
         _complex_rows(eff_shift=(0.0, 3.0, 0.0)))
    _confidence(out / "seed-1_sample-1" / "summary_confidences.json", 0.60)
    return out, ref


@pytest.mark.logic
def test_the_aggregate_reduces_the_sample_set_to_best_and_mean(tmp_path):
    """Two samples at ra_eff 0.0 and 3.0, iPTM 0.80 and 0.60. Best ra_eff is
    the minimum because lower is better, and best iPTM is the maximum because
    higher is. Reducing both the same way would invert one of them."""
    out, ref = _tree(tmp_path)
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.total_predictions == 2
    assert agg.best_ra_eff == pytest.approx(0.0, abs=1e-5)
    assert agg.mean_ra_eff == pytest.approx(1.5, abs=1e-5)
    assert agg.best_iptm == pytest.approx(0.80)
    assert agg.mean_iptm == pytest.approx(0.70)
    assert agg.failures == []


@pytest.mark.logic
def test_the_summary_row_is_the_six_columns_at_their_written_precision(
        tmp_path):
    """These strings land in a CSV that is read back by the plots and by the
    cross-sequence summary. ra_eff carries three decimals and iPTM four, and
    the failures column is comma-joined, so the whole dict is pinned rather
    than one field."""
    out, ref = _tree(tmp_path)
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.to_summary_row("design_1_seq_0") == {
        "seq_name": "design_1_seq_0",
        "af3_nomsa_best_ra_eff": "0.000",
        "af3_nomsa_mean_ra_eff": "1.500",
        "af3_nomsa_best_iptm": "0.8000",
        "af3_nomsa_mean_iptm": "0.7000",
        "af3_nomsa_total_predictions": "2",
        "af3_nomsa_failures": "",
    }


@pytest.mark.logic
def test_cif_paths_are_parallel_to_the_ra_eff_values(tmp_path):
    """The docstring calls cif_paths the only route from the aggregate back
    to the individual predictions, which only holds if the two lists line up.
    Sample 0 scored 0.0 and sample 1 scored 3.0."""
    out, ref = _tree(tmp_path)
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    paths = agg.cif_paths
    assert [p.parent.name for p in paths] == ["seed-1_sample-0",
                                             "seed-1_sample-1"]
    assert agg.best_ra_eff == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
def test_the_walk_happens_once_and_is_cached(tmp_path):
    """Documented behaviour: the walk runs on first access to any aggregate
    property. A second sample added afterwards must not appear, because the
    orthogonal-metrics stage reads several properties off one instance and
    they have to describe the same set."""
    out, ref = _tree(tmp_path)
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.total_predictions == 2
    _cif(out / "seed-1_sample-2" / "model.cif", _complex_rows())
    assert agg.total_predictions == 2
    assert len(agg.cif_paths) == 2


@pytest.mark.logic
def test_failures_is_a_copy_the_caller_cannot_corrupt(tmp_path):
    out, ref = _tree(tmp_path)
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    agg.failures.append("invented")
    assert agg.failures == []


# ── What a missing or unusable tree reports ──────────────────────────

@pytest.mark.logic
def test_a_missing_output_directory_is_reported_not_raised(tmp_path):
    """An AF3 job that never produced a directory must come back as a row
    carrying a failure tag, because AF3 is diagnostic and is not in the
    orthogonal gate. Raising here would discard the gate's own metrics."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    absent = tmp_path / "never_ran"
    agg = af3.AF3ConfidenceAggregate.from_output_dir(absent, ref)
    assert agg.failures == [f"af3_output_dir_missing:{absent}"]
    assert agg.best_ra_eff is None
    assert agg.mean_ra_eff is None
    assert agg.best_iptm is None
    assert agg.mean_iptm is None
    assert agg.total_predictions == 0
    assert agg.cif_paths == []


@pytest.mark.logic
def test_an_empty_output_directory_reports_no_cif_files(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    out = tmp_path / "af3_out"
    out.mkdir()
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.failures == ["no_cif_files_found"]


@pytest.mark.logic
def test_the_empty_row_writes_blanks_rather_than_the_word_none(tmp_path):
    """An empty cell is what the CSV readers downstream treat as missing.
    The string 'None' would parse as a value."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    agg = af3.AF3ConfidenceAggregate.from_output_dir(tmp_path / "gone", ref)
    row = agg.to_summary_row("design_0_seq_0")
    assert row["af3_nomsa_best_ra_eff"] == ""
    assert row["af3_nomsa_mean_iptm"] == ""
    assert row["af3_nomsa_total_predictions"] == "0"
    assert row["af3_nomsa_failures"].startswith("af3_output_dir_missing:")


@pytest.mark.logic
def test_an_older_layout_falls_back_to_the_broad_glob(tmp_path):
    """Without a seed-N_sample-M directory the walk globs every .cif, and
    excludes the `<name>_model.cif` duplicate AF3 writes beside the
    per-sample files. Counting that duplicate would double every prediction
    and halve the apparent mean."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    out = tmp_path / "af3_out"
    _cif(out / "pred" / "thing.cif", _complex_rows())
    _cif(out / "thing_model.cif", _complex_rows())
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.total_predictions == 1
    assert [p.name for p in agg.cif_paths] == ["thing.cif"]


@pytest.mark.logic
def test_one_unusable_prediction_is_counted_and_the_rest_survive(tmp_path,
                                                                 capsys):
    """A prediction whose effector chain is absent cannot be scored. It is
    counted in a failure tag and the usable sample still reports, because one
    bad diffusion sample must not discard the seed."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    out = tmp_path / "af3_out"
    _cif(out / "seed-1_sample-0" / "model.cif", _complex_rows())
    _cif(out / "seed-1_sample-1" / "model.cif",
         [("A", i + 1, i * 3.8, 0.0, 0.0) for i in range(8)])
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.total_predictions == 1
    assert agg.failures == ["per_prediction_rmsd_errors:1"]
    assert "RMSD failed" in capsys.readouterr().err


@pytest.mark.logic
def test_every_prediction_failing_reports_both_tags(tmp_path):
    """`no_successful_rmsd` is the one that says the aggregate is empty.
    Without it a reader would see two per-prediction errors and no statement
    that nothing came through."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    out = tmp_path / "af3_out"
    rec_only = [("A", i + 1, i * 3.8, 0.0, 0.0) for i in range(8)]
    _cif(out / "seed-1_sample-0" / "model.cif", rec_only)
    _cif(out / "seed-1_sample-1" / "model.cif", rec_only)
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.failures == ["per_prediction_rmsd_errors:2",
                            "no_successful_rmsd"]
    assert agg.best_ra_eff is None


@pytest.mark.logic
def test_a_missing_confidence_file_does_not_discard_the_pose(tmp_path,
                                                             capsys):
    """ra_eff and iPTM come from different files. A sample with a usable pose
    and no confidence JSON still contributes its ra_eff, and the iPTM gap is
    reported on stderr rather than as a failure tag, because the pose is the
    number that matters."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    out = tmp_path / "af3_out"
    _cif(out / "seed-1_sample-0" / "model.cif", _complex_rows())
    agg = af3.AF3ConfidenceAggregate.from_output_dir(out, ref)
    assert agg.total_predictions == 1
    assert agg.best_ra_eff == pytest.approx(0.0, abs=1e-5)
    assert agg.best_iptm is None
    assert agg.failures == []
    assert "iPTM read failed: confidence_json_not_found" in capsys.readouterr().err


# ── Reading the confidence JSON ──────────────────────────────────────

@pytest.mark.logic
def test_confidence_is_read_from_the_sibling_json(tmp_path):
    cif = tmp_path / "seed-1_sample-0" / "model.cif"
    cif.parent.mkdir(parents=True)
    cif.write_text("")
    _confidence(cif.parent / "summary_confidences.json", 0.7654)
    assert af3._read_confidence(cif) == (pytest.approx(0.7654), None)


@pytest.mark.logic
def test_confidence_is_also_read_from_the_per_model_name(tmp_path):
    """The second candidate drops a `_model` suffix from the stem, so
    `abc_model.cif` is paired with `abc_summary_confidences.json`. That is
    the flat layout AF3 writes at the top of an output directory."""
    cif = tmp_path / "abc_model.cif"
    cif.write_text("")
    _confidence(tmp_path / "abc_summary_confidences.json", 0.5)
    assert af3._read_confidence(cif) == (pytest.approx(0.5), None)


@pytest.mark.logic
def test_no_confidence_file_is_a_named_tag(tmp_path):
    cif = tmp_path / "model.cif"
    cif.write_text("")
    assert af3._read_confidence(cif) == (None, "confidence_json_not_found")


@pytest.mark.logic
def test_an_unparseable_confidence_file_names_the_exception(tmp_path):
    cif = tmp_path / "model.cif"
    cif.write_text("")
    (tmp_path / "summary_confidences.json").write_text("{not json")
    val, err = af3._read_confidence(cif)
    assert val is None
    assert err == "confidence_load_failed:JSONDecodeError"


@pytest.mark.logic
def test_a_null_iptm_reads_as_absent_rather_than_zero(tmp_path):
    """AF3 writes the key with a null value when it has no interface to
    score. Reading that as 0.0 would put a real-looking worst-case iPTM into
    the mean."""
    cif = tmp_path / "model.cif"
    cif.write_text("")
    _confidence(tmp_path / "summary_confidences.json", None)
    assert af3._read_confidence(cif) == (None, "confidence_json_not_found")


# ── The two thresholded readings, on their boundaries ────────────────

def _populated(ra_effs):
    agg = af3.AF3ConfidenceAggregate(Path("/dev/null"), Path("/dev/null"))
    agg._parsed = True
    agg._ra_effs = list(ra_effs)
    return agg


@pytest.mark.logic
def test_n_correct_interface_is_strictly_below_the_cutoff(tmp_path):
    """The docstring says below. A pose at exactly the cutoff is not counted,
    so 1.9 passes, 2.0 does not."""
    agg = _populated([1.9, 2.0, 2.1])
    assert agg.n_correct_interface(_Thresholds(2.0)) == 1


@pytest.mark.logic
def test_n_correct_interface_is_zero_when_nothing_parsed():
    assert _populated([]).n_correct_interface(_Thresholds(5.0)) == 0


@pytest.mark.logic
@pytest.mark.parametrize("best,disagrees", [(1.999, False), (2.0, True),
                                            (2.001, True)])
def test_af3_disagrees_at_or_above_the_cutoff(best, disagrees):
    """The complement of n_correct_interface: that one is strictly below, this
    one is at or above, so a pose exactly on the cutoff counts as a
    disagreement and is not counted as correct. Both readings of the boundary
    are live at once and they do not overlap."""
    assert _populated([best]).af3_disagrees(_Thresholds(2.0)) is disagrees


@pytest.mark.logic
def test_no_data_is_not_a_disagreement():
    """An AF3 job that produced nothing is a gap in the diagnostics, not
    evidence against the design. Reading it as disagreement would turn a
    failed side job into a verdict on the binder."""
    assert _populated([]).af3_disagrees(_Thresholds(2.0)) is False


# ── Receptor-side pairing guards ─────────────────────────────────────

@pytest.mark.logic
def test_too_few_receptor_residues_to_align_is_a_named_tag(tmp_path):
    """Kabsch needs three points. A two-residue receptor would give a
    degenerate rotation, so the tag names the count it got on each side."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb", _complex_rows(n_rec=2))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert rmsd is None
    assert err.startswith("too_few_receptor_residues:2")
    assert "pred=2" in err and "ref=8" in err


@pytest.mark.logic
def test_a_receptor_length_mismatch_warns_and_pairs_on_the_overlap(tmp_path,
                                                                   capsys):
    """Pairing is positional, so unequal receptors pair on the first N. That
    is a silent reinterpretation of the frame, which is why it prints. The
    first six residues of this pred lie exactly on the ref's, so the
    effector still scores zero."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows(n_rec=8))
    pred = _write(tmp_path / "pred.pdb", _complex_rows(n_rec=6))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-5)
    msg = capsys.readouterr().err
    assert "receptor length mismatch" in msg
    assert "pairing on first 6" in msg
