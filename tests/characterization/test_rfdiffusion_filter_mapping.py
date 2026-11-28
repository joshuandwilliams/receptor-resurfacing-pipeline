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

The per-region section covers the four anchoring cases
`calc_motif_and_region_metrics` documents. Each one decides which input
residues a built region is compared against, and three of them produce a NaN
that must not reach the metrics JSON as a zero. A zero reads as 'the region
was built and did not move', which is the opposite of 'there was nothing to
compare it with'.
"""
from __future__ import annotations

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
def test_an_already_sampled_denovo_length_becomes_a_one_long_range():
    """A bare `5` is a length RFdiffusion has already chosen, so it arrives
    as a de novo segment whose minimum and maximum are both 5. Dropping it
    would leave the whole region to the leftover path."""
    assert rf._parse_receptor_segments("A1-10/5/A20-30 B", "A") == [
        ("fixed", 1, 10), ("denovo", 5, 5), ("fixed", 20, 30)]


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
def test_every_fixed_length_denovo_segment_gets_its_declared_length():
    """Two de novo segments, both already sampled: 6-6 and 4-4. 30 atoms
    minus 20 fixed leaves exactly 10, so there is nothing left to distribute
    and each segment takes what it declared."""
    rec = _atoms(30)
    _map, meta = rf.build_receptor_resnum_map(
        rec, "A1-10/6-6/A11-20/4-4 B", "A")
    assert meta["per_region_lengths"] == [6, 4]


@pytest.mark.logic
def test_a_variable_segment_is_clamped_to_its_declared_range(capsys):
    """36 atoms, 20 fixed, so 16 are de novo. The 6-6 segment takes 6, which
    leaves 10 for a segment declared 2-4. The declared maximum wins, because
    honouring the contig matters more than absorbing the budget, and the
    mismatch is warned about rather than hidden.
    """
    rec = _atoms(36)
    _map, meta = rf.build_receptor_resnum_map(
        rec, "A1-10/6-6/A11-20/2-4 B", "A")
    warning = capsys.readouterr().err
    assert "after clamping to [2, 4]" in warning
    assert "receptor length suggests 10" in warning
    # The 6 atoms the clamp refused to allocate are the leftover, and the
    # leftover is appended to the last region: 4 + 6 = 10.
    assert meta["per_region_lengths"] == [6, 10]


@pytest.mark.logic
def test_a_declared_minimum_larger_than_the_budget_is_truncated_at_the_atoms(capsys):
    """32 atoms, 20 fixed and 6 taken by the 6-6 segment, so 6 are left for a
    segment declared 8-12. The clamp hands it 8, more than exist, and the
    walk stops at the last atom rather than mapping past the end."""
    rec = _atoms(32)
    _map, meta = rf.build_receptor_resnum_map(
        rec, "A1-10/6-6/A11-20/8-12 B", "A")
    assert "receptor length suggests 6" in capsys.readouterr().err
    assert meta["per_region_lengths"] == [6, 6]
    assert sum(meta["per_region_lengths"]) + 20 == 32


@pytest.mark.logic
def test_atoms_the_contig_does_not_account_for_become_a_design_region():
    """An all-fixed contig against a design PDB holding 15 receptor atoms.
    The contig names 10, so the 5 extras cannot be motif and are recorded as
    a de novo region with synthetic numbers above every fixed residue:
    10 + 1000 = 1010 upwards.
    """
    rec = _atoms(15)
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-10 B", "A")
    assert meta["per_region_lengths"] == [5]
    assert meta["per_region_atom_ranges"] == [(10, 15)]
    assert [resnum_map[i] for i in range(10, 15)] == [1010, 1011, 1012, 1013, 1014]


@pytest.mark.logic
def test_an_inverted_denovo_range_is_accepted_as_written():
    """`1-0` is a typo: contig_spec refuses it outright. This parser has its
    own grammar and reads it as min 1, max 0, which makes both segments
    variable with a total maximum of zero. Today the split then falls back to
    an even division, clamped up to each minimum, and the last segment takes
    the rest: 1 and 4 out of a 5-atom budget.
    """
    rec = _atoms(25)
    _map, meta = rf.build_receptor_resnum_map(
        rec, "A1-10/1-0/A11-20/1-0 B", "A")
    assert meta["per_region_lengths"] == [1, 4]


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect, pinned rather than fixed: _parse_receptor_segments accepts an "
    "inverted de novo range such as '1-0', which contig_spec.DeNovoSegment "
    "rejects with 'max_len < min_len'. The two parsers read the same grammar "
    "and only one validates it, so a typo in a hand-written contig reaches "
    "the budget split and produces a mapping for a contig that cannot be "
    "satisfied. Fixing it means raising from a parser the stage calls once "
    "per design, which changes rfdiffusion_filter's failure behaviour, so it "
    "is reported here instead."))
def test_an_inverted_denovo_range_should_be_refused():
    with pytest.raises(ValueError):
        rf._parse_receptor_segments("A1-10/1-0/A11-20 B", "A")


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


@pytest.mark.logic
@pytest.mark.parametrize("input_atoms,design_atoms", [
    ([], "design"),      # no reference to align against
    ("design", []),      # nothing to measure
])
def test_nothing_to_align_gives_nan_rather_than_zero(input_atoms, design_atoms):
    """Zero would read as a perfectly held motif, which is the opposite of
    'the measurement could not be made'. The JSON writer turns NaN into
    null, so the distinction survives into the metrics file."""
    rec = _atoms(20)
    resnum_map, meta = rf.build_receptor_resnum_map(rec, "A1-20 B", "A")
    got = rf.calc_motif_and_region_metrics(
        rec if input_atoms == "design" else [],
        rec if design_atoms == "design" else [],
        resnum_map, set(range(1, 21)),
        rf._parse_receptor_segments("A1-20 B", "A"),
        meta["per_region_atom_ranges"])
    assert math.isnan(got[0])
    assert got[1:] == ([], [], [])


@pytest.mark.logic
def test_a_fixed_residue_missing_from_the_input_is_dropped_from_the_fit():
    """The contig fixes A1-20 but the input PDB only holds 1-10, so the
    contig and the structure disagree. The ten residues that do exist still
    align, and nothing warns: the module's own comment says resolve_contigs
    only warns about a missing chain, not a missing residue.
    """
    design = _atoms(20)
    resnum_map, meta = rf.build_receptor_resnum_map(design, "A1-20 B", "A")
    motif_rmsd, _ep_in, _ep_des, _com = rf.calc_motif_and_region_metrics(
        _atoms(10), design, resnum_map, set(range(1, 21)),
        rf._parse_receptor_segments("A1-20 B", "A"),
        meta["per_region_atom_ranges"])
    assert motif_rmsd == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_a_fixed_set_that_matches_nothing_gives_nan():
    """Not one design atom maps to a residue the caller called fixed, so
    there is no motif to superpose on."""
    design = _atoms(20)
    resnum_map, meta = rf.build_receptor_resnum_map(design, "A1-20 B", "A")
    motif_rmsd, ep_in, ep_des, com = rf.calc_motif_and_region_metrics(
        _atoms(20), design, resnum_map, {999},
        rf._parse_receptor_segments("A1-20 B", "A"),
        meta["per_region_atom_ranges"])
    assert math.isnan(motif_rmsd)
    assert (ep_in, ep_des, com) == ([], [], [])


# ── Per-region spans and displacement ────────────────────────────────

def _region_metrics(contig, n_design, n_input, fixed, seg_descs=None):
    """Both metrics for one contig, with the map and the ranges that the
    stage itself would build from that contig."""
    design = _atoms(n_design)
    resnum_map, meta = rf.build_receptor_resnum_map(design, contig, "A")
    return rf.calc_motif_and_region_metrics(
        _atoms(n_input), design, resnum_map, fixed,
        seg_descs if seg_descs is not None
        else rf._parse_receptor_segments(contig, "A"),
        meta["per_region_atom_ranges"])


@pytest.mark.logic
def test_a_region_that_produced_no_atoms_reports_nan_for_everything():
    """`0-0` between two anchors, and a design with exactly the fixed atom
    count, so the region is empty. An empty region has no span and no centre
    of mass, and reporting 0 for either would read as a region that was
    built and did not move."""
    _rmsd, ep_in, ep_des, com = _region_metrics(
        "A1-10/0-0/A11-20 B", 20, 20, set(range(1, 21)))
    assert [math.isnan(v) for v in ep_in + ep_des + com] == [True] * 3


@pytest.mark.logic
def test_a_one_residue_region_has_no_span_to_measure():
    """A span needs two endpoints. The design side of a one-atom region is
    therefore NaN, and the input side is NaN too because the anchors A1-10
    and A11-20 are adjacent, leaving no input residues in the gap.

    The motif RMSD is the interesting number: inserting one residue pushes
    the second anchor 3.8 A downstream of where the input has it, and the
    best fit over 10 held and 10 shifted atoms splits that evenly, giving
    3.8 / 2 = 1.9 A.
    """
    rmsd, ep_in, ep_des, com = _region_metrics(
        "A1-10/1-1/A11-20 B", 21, 20, set(range(1, 21)))
    assert rmsd == pytest.approx(1.9)
    assert [math.isnan(v) for v in ep_in + ep_des + com] == [True] * 3


@pytest.mark.logic
def test_a_c_terminal_extension_is_compared_with_the_residues_after_the_anchor():
    """The contig ends in a de novo segment, so the region has a leading
    anchor only. The docstring's rule is to compare against the N input
    residues immediately after that anchor, where N is the region's own atom
    count: 5 atoms against input residues 11 to 15.

    Those are exactly where the extension sits, so the displacement is 0,
    while the input span stays NaN because one anchor cannot define a span.
    The design span is the 5-atom region's own reach, 4 * 3.8 = 15.2 A.
    """
    _rmsd, ep_in, ep_des, com = _region_metrics(
        "A1-10/5-5 B", 15, 20, set(range(1, 11)))
    assert math.isnan(ep_in[0])
    assert ep_des[0] == pytest.approx(4 * 3.8)
    assert com[0] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_an_n_terminal_extension_is_compared_with_the_residues_before_the_anchor():
    """The mirror case: a de novo segment before the first fixed one, so the
    region has a trailing anchor only and is compared against the 5 input
    residues immediately before A11. After superposing on A11-20 the
    extension occupies exactly those positions, so the displacement is 0.
    """
    _rmsd, ep_in, ep_des, com = _region_metrics(
        "5-5/A11-20 B", 15, 20, set(range(11, 21)))
    assert math.isnan(ep_in[0])
    assert ep_des[0] == pytest.approx(4 * 3.8)
    assert com[0] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_a_region_with_neither_anchor_has_nothing_to_compare_against():
    """The fourth case in the docstring, 'free-floating chain, very rare':
    a de novo segment with no fixed segment on either side. There are no
    corresponding input residues at all, so the displacement is NaN rather
    than being measured against an arbitrary stretch.

    The contig the stage resolves always supplies at least one anchor, so
    this pins the function's documented contract rather than a state main()
    can reach: the segment list says free-floating while the residue map
    still carries the motif to superpose on.
    """
    _rmsd, ep_in, ep_des, com = _region_metrics(
        "A1-10/5-5 B", 15, 20, set(range(1, 11)),
        seg_descs=[("denovo", 5, 5)])
    assert math.isnan(ep_in[0])
    assert ep_des[0] == pytest.approx(4 * 3.8)
    assert math.isnan(com[0])


@pytest.mark.logic
def test_an_unaccounted_region_beyond_the_contig_gets_a_span_but_no_comparison():
    """The leftover region from an all-fixed contig has no de novo segment in
    the contig, so there are no anchors to look up for it at all. Its own
    span is still real and is reported; the input span and the displacement
    are NaN.
    """
    _rmsd, ep_in, ep_des, com = _region_metrics(
        "A1-10 B", 15, 20, set(range(1, 11)))
    assert math.isnan(ep_in[0])
    assert ep_des[0] == pytest.approx(4 * 3.8)
    assert math.isnan(com[0])


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
