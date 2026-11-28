"""Cover the receptor residue map and the motif metrics, including a known bug.

`build_receptor_resnum_map` assigns each design atom an input residue number,
and everything structural downstream is keyed off it: motif_rmsd, the
per-region displacements, and the design-versus-fixed contact split.

Its docstring records a bug found on 2026-08-29. Only the TOTAL de novo budget
is known from the design PDB, so the total has to be split across segments,
and the split is a guess. It is exact when at most one segment has
min_len != max_len. With two or more variable segments it is usually wrong and
every residue after the first mis-sized segment is shifted.

The measured consequence, from the same docstring: on the AVR-PikF campaign
the auto-derived contigs had wrong lengths for 62 of 64 and 54 of 64 designs,
inflating motif RMSD to 4.0 and 2.6 angstroms where RFdiffusion's own mapping
gives 0.22 and 0.11. The retained scaffold had not moved at all.

The test below demonstrates it rather than describing it. It is xfail because
fixing it means publishing the .trb file RFdiffusion already writes, which is
a pipeline change rather than a test change.
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

import rfdiffusion_filter as rf  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

# ── Receptor block parsing ───────────────────────────────────────────

@pytest.mark.logic
def test_fixed_and_denovo_segments_are_distinguished():
    segs = rf._parse_receptor_segments("A1-10/5-15/A20-30 B", "A")
    assert segs == [("fixed", 1, 10), ("denovo", 5, 15), ("fixed", 20, 30)]


@pytest.mark.logic
def test_only_the_receptor_block_is_parsed():
    """Substring matching the whole block was rejected because it gives false
    positives for a chain letter appearing elsewhere. The effector block must
    contribute nothing."""
    segs = rf._parse_receptor_segments("A1-10/5-5/A16-25 B1-90", "A")
    assert all(s[0] != "fixed" or s[1] < 90 for s in segs)
    assert ("fixed", 1, 90) not in segs


@pytest.mark.logic
def test_a_single_residue_segment_becomes_a_one_long_range():
    assert ("fixed", 7, 7) in rf._parse_receptor_segments("A7/5-5/A20-30 B", "A")


@pytest.mark.logic
def test_the_chain_break_marker_is_skipped():
    """`0` is RFdiffusion's chain break, not a segment."""
    with_break = rf._parse_receptor_segments("A1-10/0/A20-30 B", "A")
    without = rf._parse_receptor_segments("A1-10/A20-30 B", "A")
    assert with_break == without


@pytest.mark.logic
def test_a_contig_with_no_receptor_block_gives_nothing():
    assert rf._parse_receptor_segments("B1-90", "A") == []


@pytest.mark.logic
def test_chain_matching_is_case_insensitive():
    assert rf._parse_receptor_segments("a1-10/5-5/a20-30 B", "A")


# ── The residue map ──────────────────────────────────────────────────

def _atoms(n, start=1, chain="A"):
    return [{"chain": chain, "resnum": start + i, "x": i * 3.8, "y": 0.0,
             "z": 0.0} for i in range(n)]


@pytest.mark.logic
def test_an_all_fixed_contig_maps_each_atom_to_its_own_residue():
    """No de novo budget to split, so nothing can be guessed wrong."""
    rec = _atoms(10)
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-10 B", "A")
    assert [resnum_map[i] for i in range(10)] == list(range(1, 11))
    assert meta["per_region_lengths"] == []


@pytest.mark.logic
def test_one_variable_segment_absorbs_the_whole_budget_exactly():
    """The case the docstring says is exact: at most one segment with
    min_len != max_len. 25 atoms, 20 fixed, so the de novo region is 5
    whatever its declared range said."""
    rec = _atoms(25)
    _map, meta = rf.build_receptor_resnum_map(rec, "A1-10/3-9/A11-20 B", "A")
    assert meta["per_region_lengths"] == [5]


@pytest.mark.logic
def test_design_region_residues_are_numbered_above_every_fixed_residue():
    """De novo residues take a synthetic base above the highest fixed resnum.
    Downstream calls a residue 'design' when it is absent from the fixed set,
    so a collision would classify a rebuilt residue as retained."""
    rec = _atoms(25)
    resnum_map, _meta = rf.build_receptor_resnum_map(rec, "A1-10/3-9/A11-20 B", "A")
    fixed = set(range(1, 21))
    design_resnums = [v for v in resnum_map.values() if v not in fixed]
    assert design_resnums, "some residues must be design"
    assert min(design_resnums) > max(fixed)


@pytest.mark.logic
def test_the_position_order_covers_every_atom():
    rec = _atoms(25)
    _map, meta = rf.build_receptor_resnum_map(rec, "A1-10/3-9/A11-20 B", "A")
    assert len(meta["receptor_position_order"]) == 25


@pytest.mark.logic
def test_an_unparseable_contig_returns_empty_metadata():
    """Returning a partial map would key every later metric off a guess."""
    resnum_map, meta = rf.build_receptor_resnum_map(_atoms(10), "B1-90", "A")
    assert resnum_map == {}
    assert meta["per_region_lengths"] == []


@pytest.mark.logic
def test_a_fixed_length_denovo_segment_gets_exactly_that_length():
    """6-6 must be exactly 6. The earlier purely-proportional split was
    rejected because it could violate this."""
    rec = _atoms(26)
    _map, meta = rf.build_receptor_resnum_map(rec, "A1-10/6-6/A11-20 B", "A")
    assert meta["per_region_lengths"] == [6]


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Known bug recorded in build_receptor_resnum_map's own docstring, found "
    "2026-08-29. With two or more variable-length de novo segments the "
    "per-segment split is a guess, because only the total budget is "
    "recoverable from the design PDB. Fixing it means publishing the .trb "
    "RFdiffusion already writes, which carries sampled_mask and the exact "
    "input-to-output correspondence. That is a pipeline change, not a test "
    "change, so this demonstrates the bug rather than asserting a fix."))
def test_two_variable_segments_split_the_budget_correctly():
    """The design really has regions of 4 and 8, total 12.

    Both segments declare 2-10, so the code cannot tell 4+8 from 6+6 and
    distributes proportionally by max_len. Every residue after the first
    mis-sized segment is then shifted, which is what inflated motif RMSD to
    4.0 A on a scaffold that had not moved.
    """
    rec = _atoms(32)            # 20 fixed + 12 de novo
    _map, meta = rf.build_receptor_resnum_map(
        rec, "A1-10/2-10/A11-20/2-10 B", "A")
    assert meta["per_region_lengths"] == [4, 8]


# ── Motif metrics ────────────────────────────────────────────────────

@pytest.mark.logic
def test_an_unchanged_motif_has_zero_rmsd():
    """The motif is the fixed segments. If RFdiffusion held them, aligning on
    them gives zero. A non-zero value here means the scaffold moved, which is
    exactly what the metric is for."""
    rec = _atoms(20)
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-20 B", "A")
    fixed = set(range(1, 21))
    motif_rmsd, _ep_in, _ep_des, _com = rf.calc_motif_and_region_metrics(
        rec, rec, resnum_map, fixed, rf._parse_receptor_segments("A1-20 B", "A"),
        meta["per_region_atom_ranges"])
    assert motif_rmsd == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_a_rigidly_moved_motif_still_has_zero_rmsd():
    """Kabsch superposes before measuring, so translating the whole design
    must not register. A non-zero value here would mean the superposition is
    not being applied."""
    rec = _atoms(20)
    moved = [{**a, "x": a["x"] + 25.0, "y": a["y"] - 9.0, "z": a["z"] + 3.0}
             for a in rec]
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-20 B", "A")
    motif_rmsd, _ep_in, _ep_des, _com = rf.calc_motif_and_region_metrics(
        rec, moved, resnum_map, set(range(1, 21)),
        rf._parse_receptor_segments("A1-20 B", "A"),
        meta["per_region_atom_ranges"])
    assert motif_rmsd == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
def test_a_distorted_motif_registers_a_non_zero_rmsd():
    rec = _atoms(20)
    bent = [dict(a) for a in rec]
    for a in bent[10:]:
        a["y"] += 5.0
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-20 B", "A")
    motif_rmsd, _ep_in, _ep_des, _com = rf.calc_motif_and_region_metrics(
        rec, bent, resnum_map, set(range(1, 21)),
        rf._parse_receptor_segments("A1-20 B", "A"),
        meta["per_region_atom_ranges"])
    assert motif_rmsd > 1.0


# ── Reading a real design ────────────────────────────────────────────

@pytest.mark.logic
def test_the_map_lines_up_with_atoms_read_from_a_real_file(tmp_path):
    """One pass through the reader rather than hand-built dicts, so the two
    agree on ordering."""
    pdb = complex_pdb(tmp_path / "d.pdb", n_receptor=20, n_effector=5,
                      receptor_start=1)
    rec = rf.read_ca_atoms(pdb, chain="A")
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-20 B", "A")
    assert len(resnum_map) == len(rec)
    assert len(meta["receptor_position_order"]) == len(rec)
