"""Pin the two helpers ported out of the negative-steering engine.

`get_chain_sequence` and `find_contact_residues_heavy` were this repository's
only reason to import from `bin/boltz2_negative_steering.py` and
`bin/compute_metrics.py`. Those belong to the engine, which is installed
separately, so both were reimplemented on `structure_metrics`' own PDB reader.

Reimplemented rather than copied, which means the results have to be checked
against something. These tests assert the properties the pipeline depends on,
on the tracked reference structures. The engine's own copies were compared
directly before they were deleted: 80 of 80 chain sequences and 120 of 120
contact sets agreed, the latter over three cutoffs on 40 structures.

The index these functions return is positional within the chain, not a residue
number. Every design-region and true-interface file is positional too, so a
switch to residue numbers would silently shift the protected set on any
structure that does not start at residue 1.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import structure_metrics as sm  # noqa: E402


def _tracked_complexes(limit: int = 12) -> list[Path]:
    """Tracked PDBs holding at least two chains, so they have an interface."""
    out = subprocess.run(["git", "ls-files", "*.pdb"], cwd=REPO,
                         capture_output=True, text=True).stdout.split()
    picked = []
    for rel in out:
        p = REPO / rel
        if not p.exists():
            continue
        chains = []
        for line in p.read_text(errors="ignore").splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                if line[21] not in chains:
                    chains.append(line[21])
            if len(chains) >= 2:
                break
        if len(chains) >= 2:
            picked.append(p)
        if len(picked) >= limit:
            break
    return picked


COMPLEXES = _tracked_complexes()


def _chains(p: Path) -> list[str]:
    seen: list[str] = []
    for line in p.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
            if line[21] not in seen:
                seen.append(line[21])
    return seen


@pytest.mark.logic
def test_fixture_set_is_not_empty():
    """A silently empty parametrisation would make every test below vacuous."""
    assert COMPLEXES, "no tracked multi-chain PDB found"


@pytest.mark.logic
@pytest.mark.parametrize("pdb", COMPLEXES, ids=lambda p: p.name)
def test_sequence_length_matches_residue_count(pdb: Path):
    """One letter per (resnum, icode), so length equals the residue count.

    An altloc writes the same residue twice. If those were not collapsed the
    sequence would run long and every positional index after the duplicate
    would point at the wrong residue.
    """
    for ch in _chains(pdb)[:2]:
        seq = sm.get_chain_sequence(pdb, ch)
        residues = {(ln[22:26].strip(), ln[26])
                    for ln in pdb.read_text(errors="ignore").splitlines()
                    if ln.startswith("ATOM") and len(ln) > 26 and ln[21] == ch}
        assert len(seq) == len(residues)


@pytest.mark.logic
@pytest.mark.parametrize("pdb", COMPLEXES, ids=lambda p: p.name)
def test_contacts_are_positional_indices_in_range(pdb: Path):
    """Every contact indexes a real position in the receptor sequence."""
    chains = _chains(pdb)
    rec, eff = chains[0], chains[1]
    seq = sm.get_chain_sequence(pdb, rec)
    idx = sm.find_contact_residues_heavy(pdb, rec, eff, 5.0)
    assert idx == sorted(set(idx)), "indices must be sorted and unique"
    assert all(0 <= i < len(seq) for i in idx), (
        f"index outside 0..{len(seq) - 1}")


@pytest.mark.logic
@pytest.mark.parametrize("pdb", COMPLEXES, ids=lambda p: p.name)
def test_contacts_grow_monotonically_with_cutoff(pdb: Path):
    """A wider cutoff can only add contacts, never drop one.

    This is the property that catches an indexing or distance bug, because a
    wrong frame of reference does not usually preserve the ordering.
    """
    chains = _chains(pdb)
    rec, eff = chains[0], chains[1]
    tight = set(sm.find_contact_residues_heavy(pdb, rec, eff, 4.0))
    mid = set(sm.find_contact_residues_heavy(pdb, rec, eff, 5.0))
    wide = set(sm.find_contact_residues_heavy(pdb, rec, eff, 8.0))
    assert tight <= mid <= wide


@pytest.mark.logic
def test_missing_chain_raises_and_names_what_is_there():
    """A silent empty string would read as a chain with no residues."""
    pdb = COMPLEXES[0]
    with pytest.raises(ValueError) as exc:
        sm.get_chain_sequence(pdb, "9")
    assert "Available chains" in str(exc.value)


@pytest.mark.logic
def test_no_contacts_when_the_effector_chain_is_absent():
    """An absent partner yields no interface rather than an exception."""
    pdb = COMPLEXES[0]
    rec = _chains(pdb)[0]
    assert sm.find_contact_residues_heavy(pdb, rec, "9", 5.0) == []
