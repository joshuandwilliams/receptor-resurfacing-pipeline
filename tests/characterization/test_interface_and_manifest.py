"""Cover the weighted-interface maths and the survivor manifest's path remapping.

Two unrelated stages, both with failures that produce a plausible number or a
plausible file rather than an error.

The weighted Jaccard scores two interfaces against each other with a Gaussian
on the Cb-Cb distance. Its peak sits at 4 angstroms, so a pair at exactly that
distance weighs 1.0 and the weight falls away symmetrically either side. A
wrong sigma would still return a number between 0 and 1.

The manifest remaps absolute paths stamped into a cohort CSV on one machine
and read on another. A remap that silently returns the wrong file hands the
next stage a structure from a different design.
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

import compute_interface_metrics as cim  # noqa: E402
import extract_survivor_manifest as esm  # noqa: E402

gemmi = pytest.importorskip("gemmi")


# ── The Gaussian contact weight ──────────────────────────────────────

@pytest.mark.logic
def test_the_weight_peaks_at_four_angstroms():
    """4 A is the mean of the Gaussian, so a pair at exactly that distance
    carries full weight. A shifted mean would quietly re-rank every
    interface."""
    assert cim._gaussian_contact_weight(4.0) == pytest.approx(1.0)


@pytest.mark.logic
@pytest.mark.parametrize("offset", [0.5, 1.0, 2.0, 3.0])
def test_the_weight_is_symmetric_about_the_peak(offset):
    """A pair 1 A closer than the peak weighs the same as one 1 A further.
    Asymmetry would favour either compression or separation."""
    assert cim._gaussian_contact_weight(4.0 - offset) == pytest.approx(
        cim._gaussian_contact_weight(4.0 + offset))


@pytest.mark.logic
def test_the_weight_falls_monotonically_away_from_the_peak():
    further = [cim._gaussian_contact_weight(4.0 + d) for d in (0, 1, 2, 4, 8)]
    assert further == sorted(further, reverse=True)


@pytest.mark.logic
def test_the_weight_matches_the_formula_it_documents():
    """exp(-(d - 4)^2 / 2.25). Recomputed here rather than taken from the
    module, so a change to either constant fails."""
    for d in (2.0, 3.5, 4.0, 6.0, 9.0):
        expected = math.exp(-((d - 4.0) ** 2) / 2.25)
        assert cim._gaussian_contact_weight(d) == pytest.approx(expected)


@pytest.mark.logic
def test_a_distant_pair_weighs_essentially_nothing():
    """Two chains far apart must not accumulate weight. At 20 A the Gaussian
    underflows towards zero, which is why a separated complex scores 0.0
    rather than a small positive number."""
    assert cim._gaussian_contact_weight(20.0) < 1e-30


# ── Cb collection ────────────────────────────────────────────────────

def _res(serial, chain, resnum, resname, atoms):
    out = []
    for name, (x, y, z) in atoms.items():
        el = name[0]
        out.append(
            f"ATOM  {serial:5d}  {name:<3s} {resname} {chain}{resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {el}")
        serial += 1
    return out, serial


@pytest.mark.logic
def test_cb_is_used_when_present_and_ca_when_not(tmp_path):
    """Glycine has no Cb, and an RFdiffusion backbone has none anywhere. The
    Ca fallback is what stops those residues dropping out of the interface
    entirely."""
    lines, s = [], 1
    blk, s = _res(s, "A", 1, "ALA", {"CA": (0, 0, 0), "CB": (1, 0, 0)})
    lines += blk
    blk, s = _res(s, "A", 2, "GLY", {"CA": (3.8, 0, 0)})
    lines += blk
    p = tmp_path / "s.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")

    got = cim._collect_chain_cb_positions(gemmi.read_structure(str(p)), "A")
    assert len(got) == 2
    assert got[0][1][0] == pytest.approx(1.0), "ALA must use its CB"
    assert got[1][1][0] == pytest.approx(3.8), "GLY must fall back to CA"
    assert got[0][2] == "ALA" and got[1][2] == "GLY"


@pytest.mark.logic
def test_the_index_is_positional_not_the_residue_number(tmp_path):
    """Positional indexing is what makes the (i, j) keys comparable between a
    prediction and a ground truth that number differently. Returning seqid
    would make the two label different residues with the same key."""
    lines, s = [], 1
    for k, rn in enumerate((101, 102, 103)):
        blk, s = _res(s, "A", rn, "ALA", {"CA": (k * 3.8, 0, 0), "CB": (k * 3.8, 1, 0)})
        lines += blk
    p = tmp_path / "s.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    got = cim._collect_chain_cb_positions(gemmi.read_structure(str(p)), "A")
    assert [g[0] for g in got] == [0, 1, 2]


@pytest.mark.logic
def test_an_absent_chain_yields_nothing(tmp_path):
    lines, _s = _res(1, "A", 1, "ALA", {"CA": (0, 0, 0)})
    p = tmp_path / "s.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    assert cim._collect_chain_cb_positions(gemmi.read_structure(str(p)), "Z") == []


# ── Manifest path remapping ──────────────────────────────────────────

@pytest.mark.logic
def test_a_workdir_is_found_by_name(tmp_path):
    (tmp_path / "runs" / "design_3_seq_1").mkdir(parents=True)
    (tmp_path / "runs" / "design_4_seq_1").mkdir(parents=True)
    got = esm._find_workdir("design_4_seq_1", str(tmp_path / "runs" / "*"))
    assert got is not None and got.name == "design_4_seq_1"


@pytest.mark.logic
def test_a_workdir_that_is_not_there_is_none(tmp_path):
    (tmp_path / "runs").mkdir()
    assert esm._find_workdir("absent", str(tmp_path / "runs" / "*")) is None


@pytest.mark.logic
@pytest.mark.parametrize("where", ["plan.json", "cycle_0/plan.json"])
def test_the_plan_is_found_in_either_layout(tmp_path, where):
    """The engine has written it at both. Checking only one would make a real
    run look like a missing plan."""
    p = tmp_path / where
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("{}")
    assert esm._resolve_plan_json(tmp_path) == p


@pytest.mark.logic
def test_no_plan_anywhere_is_none(tmp_path):
    assert esm._resolve_plan_json(tmp_path) is None


@pytest.mark.logic
def test_a_stale_absolute_path_is_remapped_onto_the_local_workdir(tmp_path):
    """The cohort CSV stamps absolute paths at aggregation time. A fixture
    built on a Mac and read on the cluster carries /Users/... paths that do
    not resolve there."""
    wd = tmp_path / "runs" / "design_7_seq_2"
    (wd / "cycle_0").mkdir(parents=True)
    target = wd / "cycle_0" / "prediction.pdb"
    target.write_text("ATOM\n")
    stale = "/Users/someone/elsewhere/runs/design_7_seq_2/cycle_0/prediction.pdb"
    assert esm._remap_canonical_pdb(stale, "design_7_seq_2", wd) == target


@pytest.mark.logic
def test_a_remap_whose_target_is_absent_returns_none(tmp_path):
    """None sends the caller down its missing-file path. Returning the
    non-existent Path would hand the next stage a file that is not there."""
    wd = tmp_path / "runs" / "design_7_seq_2"
    wd.mkdir(parents=True)
    stale = "/Users/x/runs/design_7_seq_2/cycle_0/prediction.pdb"
    assert esm._remap_canonical_pdb(stale, "design_7_seq_2", wd) is None


@pytest.mark.logic
def test_a_path_without_the_run_marker_is_not_remapped(tmp_path):
    """Without /runs/<seq_name>/ there is nothing to anchor on, so guessing
    would be inventing a path."""
    assert esm._remap_canonical_pdb(
        "/some/other/place/prediction.pdb", "design_1_seq_1", tmp_path) is None


@pytest.mark.logic
@pytest.mark.parametrize("original", ["", None])
def test_an_empty_original_remaps_to_none(original, tmp_path):
    assert esm._remap_canonical_pdb(original, "d", tmp_path) is None


@pytest.mark.logic
def test_the_runtime_repo_root_is_the_parent_of_bin():
    """extract_survivor_manifest.py lives at <repo>/bin/, so the root is two
    levels up. A wrong root silently remaps every repo path to nothing."""
    assert esm._runtime_repo_root() == REPO


@pytest.mark.logic
def test_a_stale_repo_path_remaps_onto_this_checkout():
    """Anchors on the first canonical top-level directory in the path. The
    target has to exist, so this uses a file that does."""
    here = Path(__file__).resolve()
    stale = f"/somewhere/else/entirely/tests/{here.relative_to(REPO / 'tests')}"
    assert esm._remap_repo_path(stale) == here


@pytest.mark.logic
def test_a_repo_path_with_no_known_marker_is_none():
    assert esm._remap_repo_path("/opt/data/thing.pdb") is None


@pytest.mark.logic
@pytest.mark.parametrize("original", ["", None])
def test_an_empty_repo_path_is_none_before_any_marker_search(original):
    """plan.json can carry a null for either remappable key. Reaching the
    marker loop with it would search the empty string and return the repo
    root itself."""
    assert esm._remap_repo_path(original) is None
