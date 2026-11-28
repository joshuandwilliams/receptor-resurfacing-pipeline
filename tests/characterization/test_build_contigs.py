"""Cover contig renumbering against a real complex PDB.

This stage translates a contig written against one numbering into the
numbering of the complex actually on disk. Two things shift the numbers: the
receptor may be renumbered when the complex is assembled, and pLDDT trimming
may remove residues from the N terminus.

Getting the offset wrong fixes the wrong stretch of the receptor. RFdiffusion
accepts it and the run completes, so there is no error to notice.
"""
from __future__ import annotations

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
    import json
    m = tmp_path / "map.json"
    m.write_text(json.dumps({"A": 1, "B": 2}))
    assert bc.load_mapping(m) == {"A": 1, "B": 2}
