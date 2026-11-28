"""Cover the negative-control sequence builders.

These produce the scrambled and polyA controls the campaign compares real
designs against. A control that is accidentally not a control makes the
comparison meaningless, and nothing in the pipeline would report it: the
controls run through exactly the same stages and produce exactly the same
columns.

Two properties carry the whole thing. The scramble must preserve composition,
so length, charge and hydrophobicity are held constant and only the
arrangement moves. And both controls must leave every residue outside the
design region untouched, or they stop being controls for the design region
and become controls for the whole receptor.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import build_control_sequences as bcs  # noqa: E402

RECEPTOR = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"      # 33 residues
DESIGN_1B = list(range(11, 21))                      # positions 11..20


# ── Design-region index file ─────────────────────────────────────────

@pytest.mark.logic
def test_indices_are_read_from_a_comma_list(tmp_path):
    f = tmp_path / "dr.txt"
    f.write_text("# derived from x.pdb\n11,12,13,14\n")
    assert bcs._load_design_region_indices_1b(f) == [11, 12, 13, 14]


@pytest.mark.logic
def test_comment_lines_are_skipped(tmp_path):
    f = tmp_path / "dr.txt"
    f.write_text("# header\n# coordinate system: 1-based\n5,6\n")
    assert bcs._load_design_region_indices_1b(f) == [5, 6]


@pytest.mark.logic
def test_a_non_integer_entry_names_the_file(tmp_path):
    """Returning a short list instead would scramble fewer positions than
    intended and the control would quietly become weaker."""
    f = tmp_path / "dr.txt"
    f.write_text("11,notanumber,13\n")
    with pytest.raises(ValueError) as exc:
        bcs._load_design_region_indices_1b(f)
    assert "dr.txt" in str(exc.value)


# ── Effector FASTA ───────────────────────────────────────────────────

@pytest.mark.logic
def test_a_wrapped_fasta_is_joined(tmp_path):
    f = tmp_path / "eff.fasta"
    f.write_text(">effector\nAAAA\nCCCC\n")
    assert bcs._load_effector_sequence(f) == "AAAACCCC"


@pytest.mark.logic
def test_an_empty_fasta_raises_rather_than_returning_nothing(tmp_path):
    """An empty effector would be written into the control's input and
    predicted as a single chain, which is not the comparison intended."""
    f = tmp_path / "eff.fasta"
    f.write_text(">effector\n\n")
    with pytest.raises(ValueError):
        bcs._load_effector_sequence(f)


# ── The scrambled control ────────────────────────────────────────────

@pytest.mark.logic
def test_the_scramble_preserves_overall_composition():
    """Composition is preserved rather than substituted, so length, charge
    and hydrophobicity are constant and only the arrangement moves."""
    out = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=1)
    assert Counter(out) == Counter(RECEPTOR)
    assert len(out) == len(RECEPTOR)


@pytest.mark.logic
def test_the_scramble_touches_nothing_outside_the_design_region():
    """Otherwise it stops being a control for the design region."""
    out = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=1)
    for i, (a, b) in enumerate(zip(RECEPTOR, out), start=1):
        if i not in DESIGN_1B:
            assert a == b, f"position {i} changed outside the design region"


@pytest.mark.logic
def test_the_scrambled_region_is_a_permutation_of_itself():
    out = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=1)
    before = [RECEPTOR[i - 1] for i in DESIGN_1B]
    after = [out[i - 1] for i in DESIGN_1B]
    assert Counter(before) == Counter(after)


@pytest.mark.logic
def test_the_same_seed_gives_the_same_sequence():
    """Re-running one design twice must give identical output, or a rerun
    silently becomes a different control."""
    a = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=42)
    b = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=42)
    assert a == b


@pytest.mark.logic
def test_different_seeds_generally_give_different_arrangements():
    """Two designs in one cohort must not receive the same shuffle, or the
    controls stop being independent."""
    outs = {bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=s)
            for s in range(12)}
    assert len(outs) > 1


@pytest.mark.logic
def test_an_empty_design_region_returns_the_sequence_unchanged():
    """Documented: the caller decides whether that is an error."""
    assert bcs._build_scrambled_receptor(RECEPTOR, [], seed=1) == RECEPTOR


@pytest.mark.logic
def test_indices_outside_the_receptor_are_dropped_not_fatal():
    """A guard against an upstream off-by-one. Indexing past the end would
    raise, and indexing from 0 would shift the whole region by one."""
    out = bcs._build_scrambled_receptor(RECEPTOR, [1, 2, 999, 0, -5], seed=1)
    assert len(out) == len(RECEPTOR)
    assert Counter(out) == Counter(RECEPTOR)
    assert out[2:] == RECEPTOR[2:]


# ── The polyA control ────────────────────────────────────────────────

@pytest.mark.logic
def test_polya_substitutes_the_design_region_with_alanine():
    out = bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)
    for i in DESIGN_1B:
        assert out[i - 1] in ("A", "G")


@pytest.mark.logic
def test_polya_preserves_native_glycine():
    """Alanine where the backbone needs glycine forces phi/psi out of the
    only Ramachandran well that fits, which would confound the geometric
    comparison the control exists to make.
    """
    seq = "MKTGGGAAAKRT"
    design = list(range(1, len(seq) + 1))
    out = bcs._build_polyA_receptor(seq, design)
    for i, residue in enumerate(seq):
        if residue == "G":
            assert out[i] == "G", f"glycine at {i + 1} must survive"
        else:
            assert out[i] == "A"


@pytest.mark.logic
def test_polya_touches_nothing_outside_the_design_region():
    out = bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)
    for i, (a, b) in enumerate(zip(RECEPTOR, out), start=1):
        if i not in DESIGN_1B:
            assert a == b


@pytest.mark.logic
def test_polya_preserves_length():
    assert len(bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)) == len(RECEPTOR)


@pytest.mark.logic
def test_polya_with_an_empty_region_changes_nothing():
    assert bcs._build_polyA_receptor(RECEPTOR, []) == RECEPTOR


@pytest.mark.logic
def test_polya_is_idempotent():
    once = bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)
    assert bcs._build_polyA_receptor(once, DESIGN_1B) == once


# ── Seed derivation ──────────────────────────────────────────────────

@pytest.mark.logic
def test_the_seed_is_stable_for_one_design():
    a = bcs._default_seed_from_pdb(Path("/x/design_3.pdb"))
    b = bcs._default_seed_from_pdb(Path("/somewhere/else/design_3.pdb"))
    assert a == b, "the seed comes from the stem, not the full path"


@pytest.mark.logic
def test_different_designs_get_different_seeds():
    seeds = {bcs._default_seed_from_pdb(Path(f"/x/design_{i}.pdb"))
             for i in range(20)}
    assert len(seeds) == 20


# ── Output ───────────────────────────────────────────────────────────

@pytest.mark.logic
def test_the_fasta_is_one_header_and_one_body_line(tmp_path):
    out = tmp_path / "nested" / "ctrl.fasta"
    bcs._write_single_record_fasta(out, "control_scrambled", "MKTAY")
    lines = out.read_text().splitlines()
    assert lines == [">control_scrambled", "MKTAY"]


@pytest.mark.logic
def test_the_output_directory_is_created(tmp_path):
    out = tmp_path / "a" / "b" / "c.fasta"
    bcs._write_single_record_fasta(out, "x", "MKT")
    assert out.is_file()
