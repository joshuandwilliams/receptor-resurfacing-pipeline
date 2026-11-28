"""Cover contig renumbering against a real complex PDB.

This stage translates a contig written against one numbering into the
numbering of the complex actually on disk. Two things shift the numbers: the
receptor may be renumbered when the complex is assembled, and pLDDT trimming
may remove residues from the N terminus.

Getting the offset wrong fixes the wrong stretch of the receptor. RFdiffusion
accepts it and the run completes, so there is no error to notice.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import build_contigs as bc  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402


def _run(tmp_path: Path, complex_path: Path, contigs: str, **extra):
    """Drive the stage the way Nextflow does: argv in, three files out."""
    argv = ["build_contigs.py", "--complex", str(complex_path),
            "--contigs", contigs]
    for k, v in extra.items():
        argv += [f"--{k.replace('_', '-')}", str(v)]
    cwd = Path.cwd()
    old = sys.argv
    os.chdir(tmp_path)
    try:
        sys.argv = argv
        bc.main()
    finally:
        sys.argv = old
        os.chdir(cwd)


def _outputs(tmp_path: Path):
    contigs = (tmp_path / "updated_contigs.txt").read_text().strip()
    params = json.loads((tmp_path / "updated_params.json").read_text())
    return contigs, params

# ── Receptor references in a contig ──────────────────────────────────

@pytest.mark.logic
def test_every_fixed_boundary_is_collected_in_order():
    """The last reference is what the offset is computed from, so the sort
    matters as much as the collection."""
    assert bc.parse_contig_receptor_refs("A1-400/20-30/A426-435 B", "A") == \
        [1, 400, 426, 435]


@pytest.mark.logic
def test_the_effector_block_contributes_nothing():
    """Only the receptor's own numbering is being translated. Including the
    effector's would compute the offset against the wrong chain."""
    assert bc.parse_contig_receptor_refs("A1-10/5-5/A16-30 B1-90", "A") == \
        [1, 10, 16, 30]


@pytest.mark.logic
def test_a_de_novo_segment_contributes_no_reference():
    """A bare range is a length, not a residue number. Treating 20-30 as
    residues 20 and 30 would corrupt the offset."""
    refs = bc.parse_contig_receptor_refs("A1-400/20-30/A426-435 B", "A")
    assert 20 not in refs and 30 not in refs


@pytest.mark.logic
def test_a_single_residue_reference_is_collected():
    assert bc.parse_contig_receptor_refs("A5/10-10/A20-30 B", "A") == [5, 20, 30]


@pytest.mark.logic
def test_chain_matching_is_case_insensitive():
    assert bc.parse_contig_receptor_refs("a1-10/5-5/a16-30 B", "A") == \
        [1, 10, 16, 30]


@pytest.mark.logic
def test_commas_are_tolerated_as_block_separators():
    assert bc.parse_contig_receptor_refs("A1-10/5-5/A16-30,B", "A") == \
        bc.parse_contig_receptor_refs("A1-10/5-5/A16-30 B", "A")


@pytest.mark.logic
def test_a_contig_with_no_receptor_block_gives_no_references():
    assert bc.parse_contig_receptor_refs("B1-90", "A") == []


# ── Chain extraction ─────────────────────────────────────────────────

@pytest.mark.logic
def test_chain_extraction_returns_that_chain_only(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=12, n_effector=5,
                    receptor_start=1, effector_start=1)
    seq, resnums = bc.extract_chain(p, "A")
    assert len(resnums) == 12
    assert len(seq) == 12, "one letter per residue"


@pytest.mark.logic
def test_chain_extraction_preserves_the_real_residue_numbers(tmp_path):
    """Not 1..N. The offset is computed from the last real number, so a
    reader that renumbered would make every offset zero."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=2,
                    receptor_start=301)
    _seq, resnums = bc.extract_chain(p, "A")
    assert resnums == [301, 302, 303, 304, 305]


@pytest.mark.logic
def test_an_absent_chain_extracts_nothing(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb")
    seq, resnums = bc.extract_chain(p, "Z")
    assert (seq, resnums) == ("", [])


# ── The offset ───────────────────────────────────────────────────────

@pytest.mark.logic
def test_no_offset_when_the_contig_already_matches_the_pdb():
    """Contig ends at 435, PDB ends at 435. Nothing to translate."""
    assert bc.compute_offset("A1-400/20-30/A426-435 B", "A",
                             list(range(1, 436)), {}) == 0


@pytest.mark.logic
def test_the_offset_is_the_gap_between_the_contig_and_the_pdb():
    """Contig written against a 435-residue receptor, PDB holds 400. The
    35-residue difference is what has to be subtracted."""
    assert bc.compute_offset("A1-400/20-30/A426-435 B", "A",
                             list(range(1, 401)), {}) == 35


@pytest.mark.logic
def test_the_offset_never_goes_negative():
    """A PDB longer than the contig means no renumbering happened. A negative
    offset would push references off the start of the chain."""
    assert bc.compute_offset("A1-10 B", "A", list(range(1, 501)), {}) == 0


@pytest.mark.logic
def test_no_references_means_no_offset():
    assert bc.compute_offset("B1-90", "A", list(range(1, 100)), {}) == 0


@pytest.mark.logic
def test_an_empty_pdb_chain_means_no_offset():
    """Nothing to translate against. Guessing here would be inventing a
    numbering."""
    assert bc.compute_offset("A1-400 B", "A", [], {}) == 0


@pytest.mark.logic
def test_trimmed_n_terminal_residues_add_to_the_offset():
    """pLDDT trimming removes residues from the N terminus, which shifts
    everything after it by that many. The contig ends at 435 against a
    400-residue PDB, so 35 comes from renumbering, and 7 more residues were
    trimmed: 35 + 7 = 42.
    """
    assert bc.compute_offset("A1-400/20-30/A426-435 B", "A",
                             list(range(1, 401)),
                             {"n_trimmed_nterm": 7}) == 42


@pytest.mark.logic
def test_a_mapping_without_the_trim_key_adds_nothing():
    """A mapping JSON that records something else must not be read as a
    trim of zero-or-anything-else."""
    assert bc.compute_offset("A1-400/20-30/A426-435 B", "A",
                             list(range(1, 401)),
                             {"something_else": 9}) == 35


# ── Rewriting the contig ─────────────────────────────────────────────

@pytest.mark.logic
def test_the_effector_block_is_rewritten_to_its_real_length():
    """The effector is carried through whole, so its block is replaced with
    the length actually present rather than whatever the user wrote."""
    out = bc.update_contigs("A1-400/20-30/A426-435 B", "A", "B",
                            eff_len=90, rec_len=435, rec_offset=0)
    assert out.endswith("B1-90")


@pytest.mark.logic
def test_receptor_references_are_shifted_by_the_offset():
    """Every fixed boundary moves by the same amount. A de novo range must
    not move, because it is a length."""
    out = bc.update_contigs("A101-200/20-30/A226-235 B", "A", "B",
                            eff_len=50, rec_len=235, rec_offset=100)
    rec_block = out.split()[0]
    assert rec_block.startswith("A1-100")
    assert "/20-30/" in rec_block, "the de novo length must not be shifted"
    assert "A126-135" in rec_block


@pytest.mark.logic
def test_a_zero_offset_leaves_the_receptor_block_alone():
    out = bc.update_contigs("A1-400/20-30/A426-435 B", "A", "B",
                            eff_len=90, rec_len=435, rec_offset=0)
    assert out.split()[0] == "A1-400/20-30/A426-435"


@pytest.mark.logic
def test_the_output_is_canonical_slash_and_space_form():
    out = bc.update_contigs("A1-400/20-30/A426-435,B", "A", "B",
                            eff_len=90, rec_len=435, rec_offset=0)
    assert "," not in out
    assert len(out.split()) == 2


@pytest.mark.logic
def test_a_chain_break_marker_is_carried_through_unshifted():
    """`0` is RFDiffusion's chain break, not a residue reference. Shifting it
    by the offset would turn it into a de novo segment of some length."""
    out = bc.update_contigs("A101-200/0/A226-235 B", "A", "B",
                            eff_len=50, rec_len=235, rec_offset=100)
    assert out.split()[0] == "A1-100/0/A126-135"


@pytest.mark.logic
def test_a_single_residue_reference_is_shifted_and_becomes_a_range():
    """A5 with an offset of 2 is residue 3 of the complex. The rewriter emits
    the canonical two-ended form A3-3, which both parsers read back as one
    residue."""
    out = bc.update_contigs("A5/10-10/A20-30 B", "A", "B",
                            eff_len=50, rec_len=235, rec_offset=2)
    assert out.split()[0] == "A3-3/10-10/A18-28"


@pytest.mark.logic
def test_the_rewritten_contig_parses_back():
    """Round trip through the repository's own parser, so this module and
    contig_utils cannot drift apart."""
    import contig_utils as cu
    out = bc.update_contigs("A101-200/20-30/A226-235 B", "A", "B",
                            eff_len=50, rec_len=235, rec_offset=100)
    descs, chain = cu.parse_block_segments(out.split()[0])
    assert chain == "A"
    assert [d[0] for d in descs] == ["fixed", "denovo", "fixed"]
    assert descs[0] == ("fixed", "A", 1, 100)


# ── Chain resolution ─────────────────────────────────────────────────

@pytest.mark.logic
def test_chains_resolve_to_the_two_present(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=4)
    out = bc.resolve_chains(p, "A", "B")
    # (rec_seq, rec_resnums, eff_seq, eff_resnums) or similar; assert the
    # lengths rather than guessing the arity.
    lens = [len(x) for x in out if isinstance(x, (str, list))]
    assert 10 in lens and 4 in lens


@pytest.mark.logic
def test_resolve_falls_back_to_the_first_two_chains_present(tmp_path):
    """A complex whose chains are not called A and B must still resolve,
    rather than returning two empty sequences that read as empty chains."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=7, n_effector=3,
                    receptor_chain="X", effector_chain="Y")
    out = bc.resolve_chains(p, "A", "B")
    lens = [len(x) for x in out if isinstance(x, (str, list))]
    assert 7 in lens and 3 in lens


@pytest.mark.logic
def test_a_mapping_file_round_trips(tmp_path):
    m = tmp_path / "map.json"
    m.write_text(json.dumps({"A": 1, "B": 2}))
    assert bc.load_mapping(m) == {"A": 1, "B": 2}


@pytest.mark.logic
@pytest.mark.parametrize("name,text", [
    ("placeholder.json", "{}"),       # Nextflow's stand-in for 'no trimming'
    ("empty.json", ""),
])
def test_a_placeholder_mapping_reads_as_no_mapping(tmp_path, name, text):
    """The channel always carries a file, so 'no trim mapping' arrives as an
    empty object. Reading it as a mapping would make `rec_map` truthy and put
    a `.get` default of 0 into the offset instead of skipping the trim
    entirely, which happens to agree today and would not if the default
    changed."""
    p = tmp_path / name
    p.write_text(text)
    assert bc.load_mapping(p) is None


@pytest.mark.logic
def test_an_absent_mapping_path_reads_as_no_mapping(tmp_path):
    assert bc.load_mapping(tmp_path / "never_written.json") is None
    assert bc.load_mapping(None) is None


# ── The stage, end to end ────────────────────────────────────────────

@pytest.mark.logic
def test_the_stage_writes_the_three_files_the_next_stage_reads(tmp_path):
    """RFDiffusion is launched from updated_contigs.txt against
    rfdiffusion_input.pdb, and updated_params.json is what the sequence
    stages read the receptor length and sequence from. A missing one fails
    the run with nothing in the log."""
    cx = complex_pdb(tmp_path / "complex.pdb", n_receptor=30, n_effector=12)
    _run(tmp_path, cx, "A1-10/5-5/A21-30 B")

    contigs, params = _outputs(tmp_path)
    assert contigs == "A1-10/5-5/A21-30 B1-12"
    assert params["receptor_len"] == 30
    assert params["effector_len"] == 12
    assert params["contigs"] == contigs
    assert params["receptor_start_pdb"] == 1
    assert params["offset_applied"] == 0
    assert (tmp_path / "rfdiffusion_input.pdb").read_text() == cx.read_text()


@pytest.mark.logic
def test_the_stage_writes_the_sequences_it_extracted(tmp_path):
    """One letter per residue, in residue order. The fixture builds every
    residue as alanine, so the receptor sequence is 30 As."""
    cx = complex_pdb(tmp_path / "complex.pdb", n_receptor=30, n_effector=12)
    _run(tmp_path, cx, "A1-10/5-5/A21-30 B")
    _contigs, params = _outputs(tmp_path)
    assert params["receptor_seq"] == "A" * 30
    assert params["effector_seq"] == "A" * 12


@pytest.mark.logic
def test_the_stage_translates_a_contig_written_against_the_original_numbering(tmp_path):
    """The contig names residues 101-130 of the uncropped receptor. The
    complex holds 1-30, so every reference moves down by 130 - 30 = 100."""
    cx = complex_pdb(tmp_path / "complex.pdb", n_receptor=30, n_effector=12)
    _run(tmp_path, cx, "A101-110/5-5/A121-130 B")
    contigs, params = _outputs(tmp_path)
    assert params["offset_applied"] == 100
    assert contigs == "A1-10/5-5/A21-30 B1-12"


@pytest.mark.logic
def test_the_stage_adds_the_trim_mapping_to_the_offset(tmp_path):
    """Renumbering accounts for 35 - 30 = 5, and the mapping says 3 more
    residues were cut from the N terminus, so the total shift is 8. A1-20
    lands on A1-12 because the start clamps at residue 1."""
    cx = complex_pdb(tmp_path / "complex.pdb", n_receptor=30, n_effector=12)
    mapping = tmp_path / "rec_trim.json"
    mapping.write_text(json.dumps({"n_trimmed_nterm": 3}))
    _run(tmp_path, cx, "A1-20/5-5/A26-35 B", rec_trim_mapping=mapping)
    contigs, params = _outputs(tmp_path)
    assert params["offset_applied"] == 8
    assert contigs == "A1-12/5-5/A18-27 B1-12"


@pytest.mark.logic
def test_the_stage_falls_back_to_the_chains_actually_present(tmp_path):
    """A complex whose chains are X and Y must still produce a contig. Two
    empty sequences would be written as a zero-length receptor and RFDiffusion
    would fail on a contig naming residues that are not there."""
    cx = complex_pdb(tmp_path / "complex.pdb", n_receptor=18, n_effector=9,
                     receptor_chain="X", effector_chain="Y")
    _run(tmp_path, cx, "A1-6/4-4/A13-18 B")
    _contigs, params = _outputs(tmp_path)
    assert params["receptor_len"] == 18
    assert params["effector_len"] == 9
