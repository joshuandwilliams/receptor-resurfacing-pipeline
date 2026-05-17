"""Unit tests for bin/structure_metrics.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import structure_metrics as sm  # noqa: E402


# ── Minimal PDB helpers ─────────────────────────────────────────────


def _atom_line(serial, name, resname, chain, resnum, x, y, z, element=" C"):
    return (
        f"ATOM  {serial:>5d} {name:<4s} {resname:>3s} {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _write_pdb(path: Path, lines):
    path.write_text("".join(lines))


# ── Numpy primitives ────────────────────────────────────────────────


@pytest.mark.local_unit
class TestCentroid:
    def test_basic(self):
        c = sm.centroid(np.array([[0, 0, 0], [2, 4, 6]], dtype=float))
        assert np.allclose(c, [1.0, 2.0, 3.0])

    def test_single_point(self):
        c = sm.centroid(np.array([[5.0, 7.0, 9.0]]))
        assert np.allclose(c, [5.0, 7.0, 9.0])

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            sm.centroid(np.zeros((0, 3)))

    def test_nan_raises(self):
        with pytest.raises(ValueError, match="NaN"):
            sm.centroid(np.array([[1.0, 2.0, float("nan")]]))


@pytest.mark.local_unit
class TestCentroidDistance:
    def test_simple(self):
        d = sm.centroid_distance(
            np.array([0.0, 0.0, 0.0]),
            np.array([3.0, 4.0, 0.0]),
        )
        assert d == pytest.approx(5.0)


# ── PDB readers ─────────────────────────────────────────────────────


@pytest.mark.local_unit
class TestReadChainCaCoords:
    def test_single_chain(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "N",  "ALA", "A", 1, 0, 0, 0),
            _atom_line(2, "CA", "ALA", "A", 1, 1.0, 2.0, 3.0),
            _atom_line(3, "C",  "ALA", "A", 1, 0, 0, 0),
            _atom_line(4, "CA", "GLY", "A", 2, 4.0, 5.0, 6.0),
            _atom_line(5, "CA", "VAL", "A", 3, 7.0, 8.0, 9.0),
        ])
        coords = sm.read_chain_ca_coords(pdb, "A")
        assert coords.shape == (3, 3)
        assert np.allclose(coords[0], [1.0, 2.0, 3.0])
        assert np.allclose(coords[2], [7.0, 8.0, 9.0])

    def test_filters_by_chain(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "B", 1, 5.0, 0.0, 0.0),
            _atom_line(3, "CA", "VAL", "B", 2, 6.0, 0.0, 0.0),
        ])
        a = sm.read_chain_ca_coords(pdb, "A")
        b = sm.read_chain_ca_coords(pdb, "B")
        assert a.shape == (1, 3)
        assert b.shape == (2, 3)

    def test_missing_chain_returns_empty(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [_atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0)])
        coords = sm.read_chain_ca_coords(pdb, "Z")
        assert coords.shape == (0, 3)

    def test_endmdl_stops(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            _atom_line(1, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0)
            + "ENDMDL\n"
            + _atom_line(2, "CA", "GLY", "A", 2, 5.0, 0.0, 0.0)
        )
        coords = sm.read_chain_ca_coords(pdb, "A")
        assert coords.shape == (1, 3)


@pytest.mark.local_unit
class TestChainCentreOfMass:
    def test_basic(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "A", 2, 2.0, 4.0, 6.0),
        ])
        com = sm.chain_centre_of_mass(pdb, "A")
        assert np.allclose(com, [1.0, 2.0, 3.0])

    def test_empty_chain_returns_none(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [_atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0)])
        assert sm.chain_centre_of_mass(pdb, "B") is None


@pytest.mark.local_unit
class TestChainComDistance:
    def test_two_chains(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "A", 2, 2.0, 0.0, 0.0),  # A COM = (1,0,0)
            _atom_line(3, "CA", "GLY", "B", 1, 1.0, 4.0, 0.0),
            _atom_line(4, "CA", "VAL", "B", 2, 1.0, 0.0, 0.0),  # B COM = (1,2,0)
        ])
        assert sm.chain_com_distance(pdb, "A", "B") == pytest.approx(2.0)

    def test_missing_chain_returns_none(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [_atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0)])
        assert sm.chain_com_distance(pdb, "A", "Z") is None


# ── Clash counting ──────────────────────────────────────────────────


@pytest.mark.local_unit
class TestClashCount:
    def test_no_clashes(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "B", 1, 10.0, 0.0, 0.0),
        ])
        assert sm.clash_count(pdb, "A", "B", set()) == (0, 0)

    def test_clash_outside_design(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 5, 0.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "B", 1, 1.5, 0.0, 0.0),
        ])
        # design_region = {10} — residue 5 is NOT in it.
        assert sm.clash_count(pdb, "A", "B", {10}) == (0, 1)

    def test_clash_inside_design(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 5, 0.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "B", 1, 1.5, 0.0, 0.0),
        ])
        # design_region = {5} — residue 5 IS in it.
        assert sm.clash_count(pdb, "A", "B", {5}) == (1, 0)

    def test_skips_hydrogens(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "H",  "ALA", "A", 1, 0.0, 0.0, 0.0, element=" H"),
            _atom_line(2, "CA", "GLY", "B", 1, 1.0, 0.0, 0.0),
        ])
        assert sm.clash_count(pdb, "A", "B", set()) == (0, 0)

    def test_multiple_pairs_per_residue(self, tmp_path):
        # Both A heavy atoms clash with both B heavy atoms (all four
        # cross-chain pairs under 2 Å) — verifies the inner loop counts
        # every pair, not just one per A residue.
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA",  "ALA", "A", 5, 0.0, 0.0, 0.0),
            _atom_line(2, "CB",  "ALA", "A", 5, 0.5, 0.5, 0.0),
            _atom_line(3, "CA",  "GLY", "B", 1, 1.0, 0.0, 0.0),
            _atom_line(4, "CA",  "VAL", "B", 2, 1.0, 1.0, 0.0),
        ])
        assert sm.clash_count(pdb, "A", "B", set()) == (0, 4)

    def test_design_split_mixes(self, tmp_path):
        # Two A residues clashing with one B atom; one A residue is in
        # the design region, the other isn't.
        pdb = tmp_path / "test.pdb"
        _write_pdb(pdb, [
            _atom_line(1, "CA", "ALA", "A", 5, 0.0, 0.0, 0.0),  # in design
            _atom_line(2, "CA", "GLY", "A", 7, 3.0, 0.0, 0.0),  # not
            _atom_line(3, "CA", "VAL", "B", 1, 1.5, 0.0, 0.0),  # clashes with 5
            _atom_line(4, "CA", "VAL", "B", 2, 4.5, 0.0, 0.0),  # clashes with 7
        ])
        # B1 at (1.5,0,0) clashes with BOTH A5 (in design) AND A7 (not).
        # B2 at (4.5,0,0) clashes with A7 only.  Total: 1 in, 2 outside.
        assert sm.clash_count(pdb, "A", "B", {5}) == (1, 2)
