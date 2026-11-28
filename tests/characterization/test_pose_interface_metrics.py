"""
tests/characterization/test_pose_interface_metrics.py
-----------------------------------------------------
Logic-tier tests for bin/pose_interface_metrics.py.

Every metric is asserted on a geometry whose answer is arithmetic.

`fixtures_runtree.complex_pdb` writes two chains of CA atoms on parallel
lines `spacing` apart along x and `separation` apart in z, one CA per
residue.  Receptor residue i sits at (spacing x i, 0, 0) and effector
residue j at (spacing x j, 0, separation), so

    d(i, j) = sqrt(separation^2 + (spacing x (i - j))^2)

and every count at every cutoff follows from counting lattice offsets.
With spacing = 3.8 Å the offsets are:

    |i - j|   0            1        2        3
    d         separation   contributions of 3.8, 7.6, 11.4 in x

The hydrogen-bond and composition tests need real atom names, so those use
a hand-written two-residue PDB instead.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT / "bin") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "bin"))

pytest.importorskip("Bio", reason="biopython not installed")

import pose_interface_metrics as pim  # noqa: E402

from tests.characterization.fixtures_runtree import complex_pdb  # noqa: E402

_SPACING = 3.8


def _two_chains(tmp_path: Path, separation: float, n: int = 3) -> Path:
    return complex_pdb(tmp_path / f"sep_{separation}.pdb",
                       n_receptor=n, n_effector=n,
                       spacing=_SPACING, separation=separation)


def _offset_distance(k: int, separation: float) -> float:
    """The CA-CA distance between residues k places apart along the chains."""
    return math.hypot(separation, _SPACING * k)


# ── Pure helpers ─────────────────────────────────────────────────────


@pytest.mark.logic
def test_dist_stats_are_the_textbook_order_statistics():
    """n, min, max, mean, median and the quartiles of [1, 2, 3, 4]."""
    got = pim._dist_stats(np.array([1.0, 2.0, 3.0, 4.0]))
    assert got == {"n": 4, "min": 1.0, "max": 4.0, "mean": 2.5,
                   "median": 2.5, "p25": 1.75, "p75": 3.25}


@pytest.mark.logic
def test_dist_stats_of_an_empty_array_reports_only_the_count():
    """No contacts means no distribution, so only `n` is emitted.

    `_pretty` then falls back to nan for mean and median rather than
    raising.
    """
    assert pim._dist_stats(np.array([])) == {"n": 0}


@pytest.mark.logic
def test_bin_distribution_bins_are_half_open_on_the_right():
    """Each bin is [lo, hi).  A value on an edge lands in the upper bin.

    4.0 goes to "4.0-5.0" rather than "0.0-4.0", and 8.0 falls off the top
    of the last bin entirely because the final edge is exclusive.
    """
    got = pim._bin_distribution(
        np.array([0.5, 4.0, 4.9, 5.0, 8.0]), [0, 4, 5, 6, 7, 8])
    assert got == {"0.0-4.0": 1, "4.0-5.0": 2, "5.0-6.0": 1,
                   "6.0-7.0": 0, "7.0-8.0": 0}
    assert sum(got.values()) == 4          # the 8.0 is not counted


@pytest.mark.logic
def test_residue_composition_classes_overlap_so_fractions_need_not_sum_to_one():
    """His is polar and aromatic; Phe is hydrophobic and aromatic.

    Five residues: ALA (hydrophobic), SER (polar), ASP (charged),
    PHE (hydrophobic + aromatic), HIS (polar + aromatic).  So
    hydrophobic 2/5, polar 2/5, charged 1/5, aromatic 2/5, which sums to
    1.4.  The module docstring states the overlap and this pins it.
    """
    got = pim._residue_composition([("ALA", 1), ("SER", 2), ("ASP", 3),
                                    ("PHE", 4), ("HIS", 5)])
    assert got == {"hydrophobic": 0.4, "polar": 0.4,
                   "charged": 0.2, "aromatic": 0.4}
    assert sum(got.values()) == pytest.approx(1.4)


@pytest.mark.logic
def test_residue_composition_of_no_interface_residues_is_all_zero():
    assert pim._residue_composition([]) == {
        "hydrophobic": 0.0, "polar": 0.0, "charged": 0.0, "aromatic": 0.0}


@pytest.mark.logic
def test_pairwise_is_the_euclidean_distance_matrix():
    a = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    b = np.array([[0.0, 4.0, 0.0]])
    got = pim._pairwise(a, b)
    assert got.shape == (2, 1)
    assert got.ravel() == pytest.approx([4.0, 5.0])


# ── Chain loading ────────────────────────────────────────────────────


@pytest.mark.logic
def test_load_chains_rejects_a_chain_the_file_does_not_have(tmp_path):
    """The error names the file and both the requested and the found
    chains, so a chain-id mix-up is diagnosable from the message."""
    pdb = _two_chains(tmp_path, 4.0)
    with pytest.raises(ValueError) as exc:
        pim.compute_metrics(pdb, "A", "Z")
    msg = str(exc.value)
    assert "A/Z" in msg
    assert "['A', 'B']" in msg


# ── CA-CA contacts ───────────────────────────────────────────────────


@pytest.mark.logic
def test_ca_contact_count_at_8A_counts_the_offsets_inside_the_cutoff(
        tmp_path):
    """3 + 3 residues at 4 Å separation: the |i - j| <= 1 band, 7 pairs.

        |i - j| = 0 -> 4.000 Å  (3 pairs)  inside
        |i - j| = 1 -> 5.517 Å  (4 pairs)  inside
        |i - j| = 2 -> 8.589 Å  (2 pairs)  outside the 8 Å cutoff

    So 3 + 4 = 7, the minimum is 4.0 and the maximum 5.517.
    """
    assert _offset_distance(2, 4.0) > 8.0 > _offset_distance(1, 4.0)
    m = pim.compute_metrics(_two_chains(tmp_path, 4.0), "A", "B")
    ca = m["ca_contacts_8A"]
    assert ca["count"] == 7
    assert ca["distance_stats"]["n"] == 7
    assert ca["distance_stats"]["min"] == pytest.approx(4.0)
    assert ca["distance_stats"]["max"] == pytest.approx(
        _offset_distance(1, 4.0), abs=1e-3)
    assert ca["distance_histogram"] == {"0.0-4.0": 0, "4.0-5.0": 3,
                                        "5.0-6.0": 4, "6.0-7.0": 0,
                                        "7.0-8.0": 0}


@pytest.mark.logic
def test_ca_contact_count_grows_as_the_chains_close(tmp_path):
    """At 1.5 Å separation the |i - j| = 2 offset also comes inside 8 Å, so
    all 9 pairs count."""
    assert _offset_distance(2, 1.5) < 8.0
    m = pim.compute_metrics(_two_chains(tmp_path, 1.5), "A", "B")
    assert m["ca_contacts_8A"]["count"] == 9


@pytest.mark.logic
def test_no_contacts_at_all_when_the_chains_are_far_apart(tmp_path):
    m = pim.compute_metrics(_two_chains(tmp_path, 30.0), "A", "B")
    assert m["ca_contacts_8A"]["count"] == 0
    assert m["ca_contacts_8A"]["distance_stats"] == {"n": 0}
    assert m["heavy_contacts_5A"]["count"] == 0
    assert m["interface_residues"]["count_a"] == 0
    assert m["interface_residues"]["count_b"] == 0
    assert math.isnan(
        m["shape_complementarity_proxy"]["gap_index_mean_A"])


# ── Heavy-atom contacts and clashes ──────────────────────────────────


@pytest.mark.logic
def test_heavy_contact_count_at_5A_uses_the_tighter_cutoff(tmp_path):
    """The fixture is CA-only, so heavy atoms and CA atoms coincide.

    At 4 Å separation only the three stacked pairs are within 5 Å; the
    |i - j| = 1 offset is 5.517 Å and falls outside.
    """
    m = pim.compute_metrics(_two_chains(tmp_path, 4.0), "A", "B")
    hv = m["heavy_contacts_5A"]
    assert hv["count"] == 3
    assert hv["distance_stats"]["mean"] == pytest.approx(4.0)
    assert hv["distance_histogram"]["4.0-4.5"] == 3
    assert m["size"] == {"residues_a": 3, "residues_b": 3,
                         "heavy_atoms_a": 3, "heavy_atoms_b": 3}


@pytest.mark.logic
def test_clash_counts_at_each_threshold(tmp_path):
    """At 1.5 Å separation exactly the three stacked pairs clash.

    The next offset is 4.085 Å, which is above every clash threshold, so
    all four thresholds report the same three atom pairs.  Thresholds are
    strict: the count is `(d < thr).sum()`.
    """
    assert _offset_distance(1, 1.5) > 3.5
    m = pim.compute_metrics(_two_chains(tmp_path, 1.5), "A", "B")
    assert m["clash_counts"] == {"<2.0A": 3, "<2.5A": 3,
                                 "<3.0A": 3, "<3.5A": 3}


@pytest.mark.logic
def test_clash_thresholds_are_strictly_less_than(tmp_path):
    """A pair exactly on a threshold is not a clash.

    Separation 2.0 puts the three stacked pairs at exactly 2.0 Å, which the
    `< 2.0` threshold excludes and `< 2.5` includes.
    """
    m = pim.compute_metrics(_two_chains(tmp_path, 2.0), "A", "B")
    assert m["clash_counts"]["<2.0A"] == 0
    assert m["clash_counts"]["<2.5A"] == 3


@pytest.mark.logic
def test_no_clashes_once_the_chains_separate(tmp_path):
    m = pim.compute_metrics(_two_chains(tmp_path, 8.0), "A", "B")
    assert set(m["clash_counts"].values()) == {0}


# ── Interface residue identification ─────────────────────────────────


@pytest.mark.logic
def test_interface_residue_count_is_per_chain_and_uses_the_5A_cutoff(
        tmp_path):
    """Every residue of both 3-residue chains has an atom within 5 Å at
    4 Å separation, so both counts are 3.  CA-only ALA residues make the
    composition purely hydrophobic, since Ala is in HYDROPHOBIC alone.
    """
    m = pim.compute_metrics(_two_chains(tmp_path, 4.0), "A", "B")
    ir = m["interface_residues"]
    assert ir["count_a"] == 3
    assert ir["count_b"] == 3
    assert ir["composition_a"] == {"hydrophobic": 1.0, "polar": 0.0,
                                   "charged": 0.0, "aromatic": 0.0}
    assert ir["composition_b"] == ir["composition_a"]


@pytest.mark.logic
def test_interface_residue_count_shrinks_when_only_part_of_the_chain_touches(
        tmp_path):
    """`contact_residues` parks the non-contacting effector residues 500 Å
    away in y, so only the stated receptor positions have a partner.

    The effector's first two residues sit above receptor residues 0 and 1;
    the remaining three are moved out of range.  So chain B contributes two
    interface residues, and chain A the two beneath them.
    """
    pdb = complex_pdb(tmp_path / "partial.pdb", n_receptor=5, n_effector=5,
                      spacing=_SPACING, separation=4.0,
                      contact_residues=(0, 1))
    m = pim.compute_metrics(pdb, "A", "B")
    assert m["interface_residues"]["count_b"] == 2
    assert m["interface_residues"]["count_a"] == 2


# ── Buried surface area ──────────────────────────────────────────────


@pytest.mark.logic
def test_bsa_is_zero_when_nothing_is_buried(tmp_path):
    """BSA is free SASA minus complex SASA per chain.  Chains 30 Å apart
    occlude nothing, so both sides are exactly zero."""
    m = pim.compute_metrics(_two_chains(tmp_path, 30.0), "A", "B")
    assert m["bsa_A2"]["chain_a"] == pytest.approx(0.0, abs=1e-6)
    assert m["bsa_A2"]["chain_b"] == pytest.approx(0.0, abs=1e-6)
    assert m["bsa_A2"]["total"] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_bsa_grows_as_the_chains_approach(tmp_path):
    """Total BSA is the sum of the two per-chain figures and rises
    monotonically as the interface closes."""
    totals = []
    for separation in (8.0, 4.0, 2.0):
        m = pim.compute_metrics(_two_chains(tmp_path, separation), "A", "B")
        bsa = m["bsa_A2"]
        assert bsa["total"] == pytest.approx(bsa["chain_a"] + bsa["chain_b"])
        totals.append(bsa["total"])
    assert totals == sorted(totals)
    assert totals[0] < totals[-1]


# ── Gap index ────────────────────────────────────────────────────────


@pytest.mark.logic
def test_gap_index_is_the_mean_nearest_other_chain_distance(tmp_path):
    """Every interface atom's nearest partner is the atom directly across.

    At 4 Å separation each of the three receptor CAs is 4.0 Å from the
    effector CA above it, and vice versa, so the mean over all six
    nearest-neighbour distances is exactly 4.0.
    """
    m = pim.compute_metrics(_two_chains(tmp_path, 4.0), "A", "B")
    sc = m["shape_complementarity_proxy"]
    assert sc["gap_index_mean_A"] == pytest.approx(4.0)
    assert sc["interface_heavy_atoms_a"] == 3
    assert sc["interface_heavy_atoms_b"] == 3


@pytest.mark.logic
def test_gap_index_falls_as_packing_tightens(tmp_path):
    """Lower is tighter.  Halving the separation halves the gap index."""
    tight = pim.compute_metrics(_two_chains(tmp_path, 2.0), "A", "B")
    loose = pim.compute_metrics(_two_chains(tmp_path, 4.0), "A", "B")
    assert tight["shape_complementarity_proxy"]["gap_index_mean_A"] == (
        pytest.approx(2.0))
    assert (tight["shape_complementarity_proxy"]["gap_index_mean_A"]
            < loose["shape_complementarity_proxy"]["gap_index_mean_A"])


# ── Hydrogen-bond candidates ─────────────────────────────────────────


def _donor_acceptor_pdb(path: Path) -> Path:
    """One Ser on chain A facing one Asp on chain B.

    The only donor/acceptor pair brought inside 3.5 Å is Ser OG at the
    origin against Asp OD1 at (0, 0, 3.0).  Everything else is parked far
    away on purpose:

      * Ser N at z = -50 and Asp N at z = -60, so neither backbone amide
        reaches an acceptor.
      * Ser O at y = 100 and Asp O at y = 200, so neither backbone
        carbonyl is an acceptor in range.
      * Asp OD2 at z = 10, which is 10 Å from OG.

    So donor_a_to_acceptor_b is 1 (OG -> OD1) and donor_b_to_acceptor_a is
    0.  Note Ser OG is in both the donor and the acceptor set, which is why
    the two directions are counted separately rather than deduplicated.
    """
    path.write_text("\n".join([
        "ATOM      1  N   SER A   1       0.000   0.000 -50.000  1.00  0.00           N",
        "ATOM      2  CA  SER A   1       0.000   1.000   0.000  1.00  0.00           C",
        "ATOM      3  C   SER A   1       0.000   2.000   0.000  1.00  0.00           C",
        "ATOM      4  O   SER A   1       0.000 100.000   0.000  1.00  0.00           O",
        "ATOM      5  CB  SER A   1       0.000   0.500   0.000  1.00  0.00           C",
        "ATOM      6  OG  SER A   1       0.000   0.000   0.000  1.00  0.00           O",
        "TER",
        "ATOM      7  N   ASP B   1       0.000   0.000 -60.000  1.00  0.00           N",
        "ATOM      8  CA  ASP B   1       0.000   1.000   3.000  1.00  0.00           C",
        "ATOM      9  C   ASP B   1       0.000   2.000   3.000  1.00  0.00           C",
        "ATOM     10  O   ASP B   1       0.000 200.000   0.000  1.00  0.00           O",
        "ATOM     11  CB  ASP B   1       0.000   0.500   3.000  1.00  0.00           C",
        "ATOM     12  CG  ASP B   1       0.000   0.200   3.000  1.00  0.00           C",
        "ATOM     13  OD1 ASP B   1       0.000   0.000   3.000  1.00  0.00           O",
        "ATOM     14  OD2 ASP B   1       0.000   0.000  10.000  1.00  0.00           O",
        "END",
    ]) + "\n")
    return path


@pytest.mark.logic
def test_hetero_records_are_excluded_from_every_atom_set(tmp_path):
    """Waters and ligands do not count as residues or heavy atoms.

    Both `_ca_atoms` and `_heavy_atoms` skip any residue whose hetero flag
    is set, so adding a water next to the interface changes nothing.  A
    water in range would otherwise inflate the heavy-atom contact count and
    the interface residue count.
    """
    clean = _donor_acceptor_pdb(tmp_path / "clean.pdb")
    with_water = tmp_path / "water.pdb"
    with_water.write_text(clean.read_text().replace(
        "END",
        "HETATM   15  O   HOH A   2       0.000   0.000   1.500  1.00  0.00           O\n"
        "END"))
    before = pim.compute_metrics(clean, "A", "B")
    after = pim.compute_metrics(with_water, "A", "B")
    assert after["size"] == before["size"]
    assert after["heavy_contacts_5A"]["count"] == (
        before["heavy_contacts_5A"]["count"])
    assert after["interface_residues"]["count_a"] == (
        before["interface_residues"]["count_a"])


@pytest.mark.logic
def test_hbond_candidates_count_one_directed_pair(tmp_path):
    pdb = _donor_acceptor_pdb(tmp_path / "hb.pdb")
    m = pim.compute_metrics(pdb, "A", "B")
    hb = m["hbond_candidates_3p5A"]
    assert hb == {"donor_a_to_acceptor_b": 1,
                  "donor_b_to_acceptor_a": 0,
                  "total": 1}


@pytest.mark.logic
def test_hbond_candidates_are_zero_without_a_donor_in_range(tmp_path):
    """The CA-only fixture has no donor or acceptor atom names at all."""
    m = pim.compute_metrics(_two_chains(tmp_path, 2.0), "A", "B")
    assert m["hbond_candidates_3p5A"]["total"] == 0


@pytest.mark.logic
def test_composition_reflects_the_real_residue_names(tmp_path):
    """Ser is polar only; Asp is charged only.

    Both residues are at the interface (OG-OD1 is 3.0 Å), so each chain's
    composition is one residue's class vector.
    """
    pdb = _donor_acceptor_pdb(tmp_path / "hb.pdb")
    m = pim.compute_metrics(pdb, "A", "B")
    assert m["interface_residues"]["count_a"] == 1
    assert m["interface_residues"]["count_b"] == 1
    assert m["interface_residues"]["composition_a"] == {
        "hydrophobic": 0.0, "polar": 1.0, "charged": 0.0, "aromatic": 0.0}
    assert m["interface_residues"]["composition_b"] == {
        "hydrophobic": 0.0, "polar": 0.0, "charged": 1.0, "aromatic": 0.0}
    assert m["size"]["heavy_atoms_a"] == 6
    assert m["size"]["heavy_atoms_b"] == 8


# ── Record shape and pretty-printing ─────────────────────────────────


@pytest.mark.logic
def test_metrics_record_names_its_input_and_chains(tmp_path):
    pdb = _two_chains(tmp_path, 4.0)
    m = pim.compute_metrics(pdb, "A", "B")
    assert m["pdb"] == str(pdb)
    assert m["chains"] == {"a": "A", "b": "B"}
    assert set(m) == {"pdb", "chains", "size", "ca_contacts_8A",
                      "heavy_contacts_5A", "clash_counts",
                      "interface_residues", "bsa_A2",
                      "shape_complementarity_proxy",
                      "hbond_candidates_3p5A"}


@pytest.mark.logic
def test_pretty_reports_every_block(tmp_path):
    m = pim.compute_metrics(_two_chains(tmp_path, 4.0), "A", "B")
    text = pim._pretty(m)
    for fragment in ("chains A/B", "CA contacts <=8A: 7",
                     "heavy contacts <=5A: 3", "clashes:",
                     "interface residues:", "composition A:",
                     "BSA: total", "gap-index", "H-bond candidates"):
        assert fragment in text, fragment


@pytest.mark.logic
def test_pretty_prints_nan_rather_than_failing_on_an_empty_distribution(
        tmp_path):
    """With no contacts there is no mean to print, and `_pretty` reads it
    through `.get(..., nan)`."""
    m = pim.compute_metrics(_two_chains(tmp_path, 30.0), "A", "B")
    text = pim._pretty(m)
    assert "CA contacts <=8A: 0  (mean nan A, median nan A)" in text


# ── CLI ──────────────────────────────────────────────────────────────


@pytest.mark.logic
def test_main_defaults_to_chains_a_and_b_and_writes_json(
        tmp_path, monkeypatch, capsys):
    """Two PDBs in one invocation produce a two-element JSON list.

    `--chains` is omitted, so both default to A:B.
    """
    first = _two_chains(tmp_path, 4.0)
    second = _two_chains(tmp_path, 2.0)
    out = tmp_path / "metrics.json"
    monkeypatch.setattr(sys, "argv", [
        "pose_interface_metrics.py", str(first), str(second),
        "--json-out", str(out)])
    pim.main()

    payload = json.loads(out.read_text())
    assert [p["pdb"] for p in payload] == [str(first), str(second)]
    assert [p["chains"] for p in payload] == [{"a": "A", "b": "B"}] * 2
    assert payload[0]["ca_contacts_8A"]["count"] == 7
    assert payload[1]["ca_contacts_8A"]["count"] == 9
    stdout = capsys.readouterr().out
    assert f"Wrote {out}" in stdout


@pytest.mark.logic
def test_main_accepts_an_explicit_chain_pair_per_pdb(
        tmp_path, monkeypatch, capsys):
    """`--chains B:A` swaps which chain is `a`, mirroring the record."""
    pdb = _two_chains(tmp_path, 4.0)
    monkeypatch.setattr(sys, "argv", [
        "pose_interface_metrics.py", str(pdb), "--chains", "B:A"])
    pim.main()
    assert "chains B/A" in capsys.readouterr().out


@pytest.mark.logic
def test_main_rejects_a_chains_list_of_the_wrong_length(
        tmp_path, monkeypatch):
    """One `--chains` entry is required per PDB, so a mismatch exits."""
    first = _two_chains(tmp_path, 4.0)
    second = _two_chains(tmp_path, 2.0)
    monkeypatch.setattr(sys, "argv", [
        "pose_interface_metrics.py", str(first), str(second),
        "--chains", "A:B"])
    with pytest.raises(SystemExit) as exc:
        pim.main()
    assert "same length" in str(exc.value)
