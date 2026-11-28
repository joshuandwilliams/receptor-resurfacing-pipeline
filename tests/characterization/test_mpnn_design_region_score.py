"""Drive bin/mpnn_design_region_score.py against a constructed MPNN run tree.

This module produces the one number the MPNN design round is ranked on. The
global MPNN score averages over the whole complex, which for a 79-residue
receptor plus an 82-residue effector is dominated by the ~150 positions nobody
designed, so two designs whose de novo regions differ completely come out a
hundredth apart. The design-region score averages over the free positions
only, and `mpnn_plots.plot_score_distribution` ranks designs by its median.

So the thing worth testing is which positions get averaged, and the test has
to state the answer before the code runs. It does that through a stub
ProteinMPNN written to disk in tests/characterization/stubs/, picked up by the
module's own top-level import (see `fixtures_mpnn_stage.mpnn_module_with_stub`
for why the stub reaches further than ProteinMPNN itself). The stub returns
log P(aa) = -(index(aa) + 1) in ProteinMPNN's own alphabet, so each position
contributes exactly `index(designed_aa) + 1` to the score: 1 for alanine, 2
for cysteine, 3 for aspartate, 21 for the X that padding inserts. A mean over
a stated set of free positions is then arithmetic, and the test says which
positions those are by listing the fixed ones in the JSONL.

What the stub does NOT cover is ProteinMPNN's network. Nothing here asserts
that a real design-region score has any particular value; it asserts that the
module averages the right positions, maps residues to the right indices, and
writes the result into the right row.

Layouts come from tests/proteinmpnn/example_output_files/mpnn/design_0/ and
the metadata header from that run's sequences/qc_metadata.csv.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Imported with no stub on sys.path, so HAS_MPNN is False here: this is the
# module as it exists outside the RFDiffusion container, and it is the one the
# pure-Python helpers and the HAS_MPNN guard are tested through.
import mpnn_design_region_score as drs  # noqa: E402
from fixtures_mpnn_stage import (  # noqa: E402
    SEQUENCE_METADATA_COLUMNS,
    fixed_positions_jsonl,
    mpnn_checkpoint,
    mpnn_module_with_stub,
    mpnn_seqs_fasta,
    qc_row,
    sequence_metadata_csv,
)
from fixtures_runtree import complex_pdb  # noqa: E402

# Chain A is six residues and chain B is two, so the featurised complex is
# eight positions long. Fixing A1-A3, A6 and all of B leaves exactly two free
# positions, A4 and A5, at global indices 3 and 4.
FIXED_ALL_BUT_TWO = {"A": [1, 2, 3, 6], "B": [1, 2]}
FIXED_EVERYTHING = {"A": [1, 2, 3, 4, 5, 6], "B": [1, 2]}

# "AAACDAAA" over those eight positions: alanine everywhere except cysteine at
# index 3 and aspartate at index 4, which are the two free ones. Per-position
# scores are therefore 1,1,1,2,3,1,1,1.
#   global  = 11/8 = 1.375   (every position)
#   region  = (2+3)/2 = 2.5  (the two free positions)
DESIGNED = "AAACDA/AA"
NATIVE = "AAAAAA/AA"


@pytest.fixture
def mod():
    """bin/mpnn_design_region_score.py with the on-disk stubs importable."""
    return mpnn_module_with_stub()


def _tree(tmp_path, design="0", *, fixed=None, designed=(DESIGNED,),
          native=NATIVE, n_receptor=6, n_effector=2, omit=()):
    """The three input directories MPNN_DESIGN_REGION_SCORE is handed.

    `omit` drops one or more of "pdb", "jsonl" and "fasta" for this design,
    which is how the stage's collected directories look when an upstream
    process failed for one backbone.
    """
    pdb_dir = tmp_path / "pdb_collected"
    mpnn_dir = tmp_path / "mpnn_collected"
    fixed_dir = tmp_path / "fixed_collected"
    for directory in (pdb_dir, mpnn_dir, fixed_dir):
        directory.mkdir(exist_ok=True)
    if "pdb" not in omit:
        complex_pdb(pdb_dir / f"design_{design}.pdb",
                    n_receptor=n_receptor, n_effector=n_effector)
    if "fasta" not in omit:
        mpnn_seqs_fasta(mpnn_dir, design, native=native,
                        designed=list(designed))
    if "jsonl" not in omit:
        fixed_positions_jsonl(fixed_dir, design,
                              fixed if fixed is not None else FIXED_ALL_BUT_TWO)
    return pdb_dir, mpnn_dir, fixed_dir


def _run_main(mod, monkeypatch, tmp_path, rows, *, checkpoint=True, extra=()):
    """Write the metadata CSV, run main(), return the parsed output rows."""
    metadata = sequence_metadata_csv(tmp_path / "qc_metadata.csv", rows)
    out = tmp_path / "scored_metadata.csv"
    argv = ["mpnn_design_region_score.py",
            "--metadata-csv", str(metadata),
            "--pdb-dir", str(tmp_path / "pdb_collected"),
            "--mpnn-fasta-dir", str(tmp_path / "mpnn_collected"),
            "--fixed-positions-dir", str(tmp_path / "fixed_collected"),
            "--output-csv", str(out), *extra]
    if checkpoint:
        argv += ["--checkpoint", str(mpnn_checkpoint(tmp_path / "weights.pt"))]
    monkeypatch.setattr(sys, "argv", argv)
    mod.main()
    with open(out) as handle:
        return list(csv.DictReader(handle)), out


def _score(mod, tmp_path, *, fixed=None, designed=(DESIGNED,), **kwargs):
    """score_sequences_for_design against a one-design tree."""
    pdb_dir, mpnn_dir, fixed_dir = _tree(
        tmp_path, fixed=fixed, designed=designed, **kwargs)
    device = mod.torch.device("cpu")
    model = mod.load_model(str(mpnn_checkpoint(tmp_path / "w.pt")), device)
    fixed_dict = mod.load_fixed_positions(
        str(fixed_dir / "fixed_positions_design_0.jsonl"))
    return mod.score_sequences_for_design(
        str(pdb_dir / "design_0.pdb"),
        str(mpnn_dir / "design_0_mpnn" / "seqs" / "design_0.fa"),
        fixed_dict, model, device)


# ── parse_mpnn_fasta ─────────────────────────────────────────────────

@pytest.mark.logic
def test_parse_mpnn_fasta_keeps_the_native_record_first(tmp_path):
    """Record 0 is ProteinMPNN's native sequence and records 1.. are the
    samples. The scorer drops record 0 by index, so the parser has to preserve
    order rather than, say, sort by header."""
    path = mpnn_seqs_fasta(tmp_path, "0", native="AAA",
                           designed=["CCC", "DDD"])
    entries = drs.parse_mpnn_fasta(str(path))
    assert [seq for _header, seq in entries] == ["AAA", "CCC", "DDD"]
    assert entries[0][0].startswith(">design_0,")
    assert entries[1][0].startswith(">T=0.25, sample=1")


@pytest.mark.logic
def test_parse_mpnn_fasta_joins_a_wrapped_sequence(tmp_path):
    """A FASTA written at 80 columns is a legal FASTA, and splitting one
    sequence into two entries would shift every sample's index by one."""
    path = tmp_path / "wrapped.fa"
    path.write_text(">native\nAAA\nCCC\n>sample\nDD\nEE\nFF\n")
    assert drs.parse_mpnn_fasta(str(path)) == [
        (">native", "AAACCC"), (">sample", "DDEEFF")]


@pytest.mark.logic
def test_parse_mpnn_fasta_of_an_empty_file_is_empty(tmp_path):
    path = tmp_path / "empty.fa"
    path.write_text("")
    assert drs.parse_mpnn_fasta(str(path)) == []


@pytest.mark.logic
def test_parse_mpnn_fasta_drops_sequence_lines_before_any_header(tmp_path):
    """A file whose first line is a sequence has no header to attach it to, so
    it is discarded rather than guessed at. Characterising the behaviour, not
    endorsing it: the result is an empty list and a silent loss."""
    path = tmp_path / "headerless.fa"
    path.write_text("AAACCC\n")
    assert drs.parse_mpnn_fasta(str(path)) == []


@pytest.mark.logic
def test_parse_mpnn_fasta_ignores_blank_lines(tmp_path):
    """Trailing newlines are stripped and the empty string joins harmlessly,
    which is what makes a file with a trailing blank line safe."""
    path = tmp_path / "blanks.fa"
    path.write_text(">a\nAAA\n\n>b\nCCC\n\n")
    assert drs.parse_mpnn_fasta(str(path)) == [(">a", "AAA"), (">b", "CCC")]


# ── load_fixed_positions ─────────────────────────────────────────────

@pytest.mark.logic
def test_load_fixed_positions_returns_the_inner_chain_map(tmp_path):
    """The JSONL is {design_name: {chain: [resnums]}}, and the scorer needs the
    inner map. The design name is dropped because the scorer re-keys it under
    the name parse_PDB derived from the PDB, which is where the two have to
    agree.
    """
    path = fixed_positions_jsonl(tmp_path, "28", {"A": [1, 2], "B": [1]})
    assert drs.load_fixed_positions(str(path)) == {"A": [1, 2], "B": [1]}


@pytest.mark.logic
def test_load_fixed_positions_reads_only_the_first_line(tmp_path):
    """One object per design per file in this pipeline, so a second line is a
    sign of concatenated files. It is ignored rather than merged."""
    path = tmp_path / "two.jsonl"
    path.write_text('{"design_0": {"A": [1]}}\n{"design_1": {"A": [9]}}\n')
    assert drs.load_fixed_positions(str(path)) == {"A": [1]}


@pytest.mark.logic
def test_load_fixed_positions_with_no_fixed_chains(tmp_path):
    """A design with nothing fixed: every position is designable. The inner
    map is empty, which later makes every position free."""
    path = fixed_positions_jsonl(tmp_path, "0", {})
    assert drs.load_fixed_positions(str(path)) == {}


# ── find_checkpoint ──────────────────────────────────────────────────

@pytest.mark.logic
def test_find_checkpoint_prefers_the_020_epoch_weights(tmp_path, monkeypatch):
    """ProteinMPNN ships v_48_002, v_48_010, v_48_020 and v_48_030. The
    pipeline's committed FASTA headers record model_name=v_48_020, so that is
    the one that has to be picked — alphabetical order alone would give
    v_48_002 and silently change every score in the run.
    """
    weights = tmp_path / "vanilla_model_weights"
    weights.mkdir()
    for name in ("v_48_002.pt", "v_48_010.pt", "v_48_020.pt", "v_48_030.pt"):
        (weights / name).write_bytes(b"")
    monkeypatch.setattr(drs, "MPNN_DIR", str(tmp_path))
    assert drs.find_checkpoint() == str(weights / "v_48_020.pt")


@pytest.mark.logic
def test_find_checkpoint_falls_back_to_the_first_pt_in_sorted_order(
        tmp_path, monkeypatch):
    """No file carrying 020. The fallback is sorted order, so b.pt loses to
    a.pt whatever order the directory lists them in."""
    weights = tmp_path / "vanilla_model_weights"
    weights.mkdir()
    for name in ("b.pt", "a.pt"):
        (weights / name).write_bytes(b"")
    monkeypatch.setattr(drs, "MPNN_DIR", str(tmp_path))
    assert drs.find_checkpoint() == str(weights / "a.pt")


@pytest.mark.logic
def test_find_checkpoint_ignores_non_pt_files(tmp_path, monkeypatch):
    weights = tmp_path / "vanilla_model_weights"
    weights.mkdir()
    (weights / "README.md").write_text("weights go here")
    monkeypatch.setattr(drs, "MPNN_DIR", str(tmp_path))
    assert drs.find_checkpoint() is None


@pytest.mark.logic
def test_find_checkpoint_with_no_weights_directory(tmp_path, monkeypatch):
    """Running outside the container. None is the signal main() turns into an
    exit 1, so it must not be an empty string or a path that does not exist."""
    monkeypatch.setattr(drs, "MPNN_DIR", str(tmp_path / "absent"))
    assert drs.find_checkpoint() is None


# ── load_model ───────────────────────────────────────────────────────

@pytest.mark.logic
def test_load_model_takes_its_hyperparameters_from_the_checkpoint(
        mod, tmp_path):
    """noise_level and num_edges were saved with the weights, so reading them
    off the checkpoint rather than hardcoding them is what lets a differently
    trained checkpoint load at all. The model is also put in eval mode, which
    matters because augment_eps adds noise in train mode.
    """
    path = mpnn_checkpoint(tmp_path / "w.pt", noise_level=0.1, num_edges=32)
    model = mod.load_model(str(path), mod.torch.device("cpu"))
    assert model.init_kwargs["augment_eps"] == 0.1
    assert model.init_kwargs["k_neighbors"] == 32
    assert model.init_kwargs["num_letters"] == 21
    assert model.is_eval is True
    assert model.loaded_state_dict == {"stub_layer.weight": [0.0]}


@pytest.mark.logic
def test_load_model_defaults_when_the_checkpoint_omits_them(mod, tmp_path):
    """ProteinMPNN's published vanilla weights carry both keys; a
    locally-retrained checkpoint may not. 0.2 and 48 are the defaults, and 48
    is the k in v_48_020's own name."""
    path = mpnn_checkpoint(tmp_path / "w.pt", noise_level=None, num_edges=None)
    model = mod.load_model(str(path), mod.torch.device("cpu"))
    assert model.init_kwargs["augment_eps"] == 0.2
    assert model.init_kwargs["k_neighbors"] == 48


# ── score_sequences_for_design ───────────────────────────────────────

@pytest.mark.logic
def test_the_region_score_averages_only_the_free_positions(mod, tmp_path):
    """The whole point of the module. Eight positions, two of them free, and
    the free two are the only ones that count: (2+3)/2 = 2.5, against a global
    mean over all eight of 11/8 = 1.375. If the two numbers came out equal the
    module would be computing nothing the global score did not already give.
    """
    results = _score(mod, tmp_path)
    assert results == [(0, pytest.approx(1.375), pytest.approx(2.5), 2)]


@pytest.mark.logic
def test_widening_the_free_region_changes_only_the_region_score(
        mod, tmp_path):
    """Freeing A6 as well adds its alanine, worth 1, to the region mean:
    (2+3+1)/3 = 2.0. The global mean is unchanged at 1.375, because the
    sequence did not change — only which positions were declared designable.
    """
    results = _score(mod, tmp_path, fixed={"A": [1, 2, 3], "B": [1, 2]})
    assert results == [(0, pytest.approx(1.375), pytest.approx(2.0), 3)]


@pytest.mark.logic
def test_every_position_free_makes_the_two_scores_equal(mod, tmp_path):
    """The degenerate case, and a check that the free mask is the only
    difference between the two means."""
    results = _score(mod, tmp_path, fixed={})
    assert results == [(0, pytest.approx(1.375), pytest.approx(1.375), 8)]


@pytest.mark.logic
def test_a_design_with_no_free_positions_scores_none(mod, tmp_path):
    """Everything fixed means there is no design region to score, and None is
    not 0.0: `mpnn_plots` drops a None from its ranking and would rank a 0.0
    first, as the best score in the run."""
    results = _score(mod, tmp_path, fixed=FIXED_EVERYTHING)
    assert results == [(0, pytest.approx(1.375), None, 0)]


@pytest.mark.logic
def test_sequences_are_numbered_from_zero_after_the_native_record(
        mod, tmp_path):
    """main() matches these indices against the metadata CSV's seq column, so
    an off-by-one here would attach every score to the wrong sequence. The
    second sample has tryptophan (index 18, worth 19) at one free position, so
    the two samples cannot be confused: 2.5 against (19+3)/2 = 11.0.
    """
    results = _score(mod, tmp_path,
                     designed=(DESIGNED, "AAAWDA/AA"))
    assert [index for index, *_ in results] == [0, 1]
    assert results[0][2] == pytest.approx(2.5)
    assert results[1][2] == pytest.approx(11.0)


@pytest.mark.logic
@pytest.mark.parametrize("separator", ["/", ":"])
def test_both_chain_separators_are_stripped(mod, tmp_path, separator):
    """ProteinMPNN writes "/" and some downstream tools write ":". A
    separator left in would be mapped to the unknown index 20, worth 21, and
    would shift every position after it by one — so the score would be wrong
    for the right-looking reason.
    """
    results = _score(mod, tmp_path,
                     designed=(DESIGNED.replace("/", separator),))
    assert results[0][2] == pytest.approx(2.5)


@pytest.mark.logic
def test_an_unrecognised_residue_letter_becomes_the_unknown_index(
        mod, tmp_path):
    """Anything outside ProteinMPNN's 21-letter alphabet maps to index 20, so
    a Z at a free position contributes 21 and the region mean is
    (21+3)/2 = 12.0 rather than being skipped. A skipped unknown would hide a
    broken sequence behind a plausible score.
    """
    results = _score(mod, tmp_path, designed=("AAAZDA/AA",))
    assert results[0][2] == pytest.approx(12.0)


@pytest.mark.logic
def test_a_sequence_shorter_than_the_structure_is_padded_with_unknowns(
        mod, tmp_path):
    """Three residues against eight structural positions. The five missing
    ones are padded with index 20, worth 21 each, so both free positions are
    padding and the region mean is 21.0. The score is deliberately terrible:
    a truncated sequence should rank last, not be quietly scored on its first
    three residues.
    """
    results = _score(mod, tmp_path, designed=("AAA",))
    assert results[0][2] == pytest.approx(21.0)
    assert results[0][1] == pytest.approx((3 * 1 + 5 * 21) / 8)


@pytest.mark.logic
def test_a_sequence_longer_than_the_structure_is_truncated(mod, tmp_path):
    """Two extra tryptophans past the end of the structure are dropped rather
    than wrapped around, so the scores match the eight-residue case exactly."""
    results = _score(mod, tmp_path, designed=(DESIGNED + "WW",))
    assert results == [(0, pytest.approx(1.375), pytest.approx(2.5), 2)]


@pytest.mark.logic
def test_a_fasta_with_only_the_native_record_scores_nothing(mod, tmp_path):
    """ProteinMPNN wrote the native and then produced no samples. The result
    is an empty list, which main() then writes as no scores rather than as a
    zero."""
    assert _score(mod, tmp_path, designed=()) == []


# ── main() ───────────────────────────────────────────────────────────

@pytest.mark.logic
def test_main_without_proteinmpnn_exits_one(monkeypatch, tmp_path):
    """The guard that makes every other failure mode moot. Outside the
    RFDiffusion container the import failed, and exiting 1 is right: a
    scored_metadata.csv with an empty score column would be published as a
    successful run and silently un-rank the whole design round.
    """
    monkeypatch.setattr(drs, "HAS_MPNN", False)
    monkeypatch.setattr(sys, "argv", [
        "mpnn_design_region_score.py", "--metadata-csv", "m.csv",
        "--pdb-dir", "p", "--mpnn-fasta-dir", "f", "--fixed-positions-dir", "j",
        "--output-csv", "o.csv"])
    with pytest.raises(SystemExit) as exc:
        drs.main()
    assert exc.value.code == 1


@pytest.mark.logic
def test_main_with_no_checkpoint_anywhere_exits_one(mod, monkeypatch, tmp_path):
    """--checkpoint empty and no weights directory to fall back on."""
    monkeypatch.setattr(mod, "MPNN_DIR", str(tmp_path / "absent"))
    with pytest.raises(SystemExit) as exc:
        _run_main(mod, monkeypatch, tmp_path, [qc_row()], checkpoint=False)
    assert exc.value.code == 1


@pytest.mark.logic
def test_main_falls_back_to_find_checkpoint_when_the_flag_is_empty(
        mod, monkeypatch, tmp_path, capsys):
    """MPNN_DESIGN_REGION_SCORE passes no --checkpoint at all, so the
    discovered path is the one every production run uses. It is logged,
    because which weights produced a score is not recoverable afterwards.
    """
    weights = tmp_path / "mpnn" / "vanilla_model_weights"
    weights.mkdir(parents=True)
    discovered = mpnn_checkpoint(weights / "v_48_020.pt")
    monkeypatch.setattr(mod, "MPNN_DIR", str(tmp_path / "mpnn"))
    _tree(tmp_path)
    rows, _out = _run_main(mod, monkeypatch, tmp_path,
                           [qc_row(design="0", seq="0")], checkpoint=False)
    assert rows[0]["design_region_score"] == "2.5000"
    assert f"Checkpoint: {discovered}" in capsys.readouterr().out


@pytest.mark.logic
def test_main_writes_the_region_score_into_the_matching_row(
        mod, monkeypatch, tmp_path, capsys):
    """The output is the input CSV with one column appended, to four decimal
    places. Row identity is (design, seq), and the second sample's score has
    to land on the second row: 2.5 for sample 0, 11.0 for sample 1.
    """
    _tree(tmp_path, designed=(DESIGNED, "AAAWDA/AA"))
    rows, out = _run_main(mod, monkeypatch, tmp_path, [
        qc_row(design="0", seq="0"), qc_row(design="0", seq="1")])
    assert [(r["design"], r["seq"], r["design_region_score"]) for r in rows] == [
        ("0", "0", "2.5000"), ("0", "1", "11.0000")]
    with open(out) as handle:
        assert next(csv.reader(handle)) == [
            *SEQUENCE_METADATA_COLUMNS, "design_region_score"]
    assert "Design region scores: 2/2 sequences scored" in capsys.readouterr().out


@pytest.mark.logic
def test_main_leaves_the_cell_empty_when_there_is_no_design_region(
        mod, monkeypatch, tmp_path, capsys):
    """A None score is written as the empty string and does not count towards
    the scored tally. `mpnn_plots` reads an empty cell as missing and appends
    that design after the ranked block, which is the behaviour this empty
    cell is what produces.
    """
    _tree(tmp_path, fixed=FIXED_EVERYTHING)
    rows, _out = _run_main(mod, monkeypatch, tmp_path,
                           [qc_row(design="0", seq="0")])
    assert rows[0]["design_region_score"] == ""
    output = capsys.readouterr().out
    assert "region=N/A (0 free)" in output
    assert "Design region scores: 0/1 sequences scored" in output


@pytest.mark.logic
def test_main_scores_each_design_once_however_many_sequences_it_has(
        mod, monkeypatch, tmp_path, capsys):
    """The model runs per design, and the featurisation is the expensive part,
    so a design seen a second time is skipped. With 64 designs times 128
    sequences, scoring per row instead would be 128 times the work."""
    _tree(tmp_path, designed=(DESIGNED, "AAAWDA/AA"))
    rows, _out = _run_main(mod, monkeypatch, tmp_path, [
        qc_row(design="0", seq="0"), qc_row(design="0", seq="1")])
    assert len(rows) == 2
    assert capsys.readouterr().out.count("design_0: scoring...") == 1


@pytest.mark.logic
@pytest.mark.parametrize("omit,expected", [
    (("pdb",), "missing PDB"),
    (("jsonl",), "missing JSONL"),
    (("fasta",), "missing FASTA"),
    (("pdb", "jsonl", "fasta"), "missing PDB, JSONL, FASTA"),
])
def test_main_names_which_inputs_a_design_is_missing(
        mod, monkeypatch, tmp_path, capsys, omit, expected):
    """The three collected directories are built by shell loops over Nextflow
    channels, so one of them can arrive incomplete. Naming which file is
    absent is the difference between a five-minute fix and reading the whole
    stage script. The design is skipped and the row is still written, with an
    empty score.
    """
    _tree(tmp_path, omit=omit)
    rows, _out = _run_main(mod, monkeypatch, tmp_path,
                           [qc_row(design="0", seq="0")])
    assert rows[0]["design_region_score"] == ""
    assert f"design_0: {expected}" in capsys.readouterr().out


@pytest.mark.logic
def test_main_survives_a_design_that_fails_to_score(
        mod, monkeypatch, tmp_path, capsys):
    """One design throwing must not cost the other 63. The exception is
    reported with its traceback and the loop continues, so design_1 is scored
    even though design_0 blew up.
    """
    _tree(tmp_path, "0")
    _tree(tmp_path, "1")
    real = mod.score_sequences_for_design

    def explode(pdb_path, *args, **kwargs):
        if pdb_path.endswith("design_0.pdb"):
            raise RuntimeError("featurisation blew up")
        return real(pdb_path, *args, **kwargs)

    monkeypatch.setattr(mod, "score_sequences_for_design", explode)
    rows, _out = _run_main(mod, monkeypatch, tmp_path, [
        qc_row(design="0", seq="0"), qc_row(design="1", seq="0")])
    assert [r["design_region_score"] for r in rows] == ["", "2.5000"]
    captured = capsys.readouterr()
    assert "design_0: ERROR: featurisation blew up" in captured.out
    assert "RuntimeError" in captured.err


@pytest.mark.logic
def test_main_with_zero_rows_writes_a_four_column_header(
        mod, monkeypatch, tmp_path):
    """Header-only input: nothing passed QC. The output header is hardcoded
    to four columns rather than echoing the ten it was given, so a reader that
    indexes qc_metadata's columns by name finds them missing. Characterising
    the wart: `mpnn_plots` takes its own no-data path on a header-only file and
    never reads those columns, which is why it has never bitten.
    """
    _tree(tmp_path)
    rows, out = _run_main(mod, monkeypatch, tmp_path, [])
    assert rows == []
    assert out.read_text() == "design,seq,mpnn_score,design_region_score\n"


@pytest.mark.logic
def test_main_with_an_absent_metadata_file_raises(mod, monkeypatch, tmp_path):
    """No guard on the input path, which is right for a missing upstream
    file: the stage fails rather than publishing an unranked round."""
    monkeypatch.setattr(sys, "argv", [
        "mpnn_design_region_score.py",
        "--metadata-csv", str(tmp_path / "absent.csv"),
        "--pdb-dir", str(tmp_path), "--mpnn-fasta-dir", str(tmp_path),
        "--fixed-positions-dir", str(tmp_path),
        "--output-csv", str(tmp_path / "o.csv"),
        "--checkpoint", str(mpnn_checkpoint(tmp_path / "w.pt"))])
    with pytest.raises(FileNotFoundError):
        mod.main()


@pytest.mark.logic
def test_main_reports_the_device_it_scored_on(mod, monkeypatch, tmp_path,
                                              capsys):
    """CPU and GPU scoring agree to within floating point, but which one ran
    decides whether a 64-design round takes minutes or hours, so it is logged.
    """
    _tree(tmp_path)
    _run_main(mod, monkeypatch, tmp_path, [qc_row(design="0", seq="0")])
    assert "Device: cpu" in capsys.readouterr().out


# ── CLI ──────────────────────────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("omit", [
    "--metadata-csv", "--pdb-dir", "--mpnn-fasta-dir",
    "--fixed-positions-dir", "--output-csv"])
def test_parse_args_requires_all_five_paths(monkeypatch, omit):
    argv = ["mpnn_design_region_score.py", "--metadata-csv", "m.csv",
            "--pdb-dir", "p", "--mpnn-fasta-dir", "f",
            "--fixed-positions-dir", "j", "--output-csv", "o.csv"]
    index = argv.index(omit)
    del argv[index:index + 2]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):
        drs.parse_args()


@pytest.mark.logic
def test_parse_args_leaves_the_checkpoint_empty_by_default(monkeypatch):
    """The empty string is what sends main() to find_checkpoint, and
    MPNN_DESIGN_REGION_SCORE relies on it."""
    monkeypatch.setattr(sys, "argv", [
        "mpnn_design_region_score.py", "--metadata-csv", "m.csv",
        "--pdb-dir", "p", "--mpnn-fasta-dir", "f",
        "--fixed-positions-dir", "j", "--output-csv", "o.csv"])
    assert drs.parse_args().checkpoint == ""


@pytest.mark.logic
def test_the_alphabet_is_proteinmpnns_own_twenty_one_letters():
    """The index of each residue in this string is what the module looks up in
    the model's output, so it has to be ProteinMPNN's ordering and not
    alphabetical or Kyte-Doolittle. X at index 20 is both the unknown residue
    and the padding value.
    """
    assert drs.ALPHABET == "ACDEFGHIKLMNPQRSTVWYX"
    assert drs.AA_TO_IDX["A"] == 0
    assert drs.AA_TO_IDX["X"] == 20
    assert len(drs.AA_TO_IDX) == 21
