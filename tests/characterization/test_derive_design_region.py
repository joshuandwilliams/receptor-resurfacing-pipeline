"""Cover the design-region derivation, and the 0/1-based seam inside it.

The design region is the set of receptor positions RFDiffusion rebuilt. It is
the complement of the protected set, so an index one out here hands the
steering engine a position it should have left alone.

The seam this module carries is worth stating plainly, because the two halves
disagree on purpose and nothing in the type signatures says so:

    derive_positional_indices()  returns 0-based indices
    format_output()              documents and writes 1-based positions

so `main()` has to add one between them. The sibling module
`derive_true_interface.py` writes 0-based instead, for the same receptor. Both
conventions are live in the repository at once, which is why every test below
names the frame it is asserting in.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import derive_design_region as ddr  # noqa: E402
import derive_true_interface as dti  # noqa: E402

# ── Design id normalisation ──────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("raw", ["design_5", "design_5.pdb", "5"])
def test_every_spelling_of_a_design_id_normalises_the_same(raw):
    assert ddr.normalise_design_id(raw) == {"bare": "design_5",
                                            "with_pdb": "design_5.pdb"}


@pytest.mark.logic
def test_the_two_sibling_modules_normalise_identically():
    """They read the same metrics file. If they disagreed on the key, one
    would silently pick a different design than the other."""
    for raw in ("design_8", "design_8.pdb", "8"):
        assert ddr.normalise_design_id(raw) == dti.normalise_design_id(raw)


# ── Finding the entry ────────────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("key", ["design_7", "design_7.pdb"])
def test_entry_is_found_under_either_key_spelling(key):
    m = {"designs": [{"design": key, "per_design_design_residues": [1]}]}
    assert ddr.find_design_entry(m, "7")["design"] == key


@pytest.mark.logic
def test_a_missing_design_names_the_ones_that_exist():
    m = {"designs": [{"design": "design_1"}, {"design": "design_2"}]}
    with pytest.raises(KeyError) as exc:
        ddr.find_design_entry(m, "design_9")
    assert "design_1" in str(exc.value)


# ── The derivation ───────────────────────────────────────────────────

@pytest.mark.logic
def test_indices_are_zero_based_positions_not_residue_numbers():
    """The receptor starts at 101. Returning residue numbers would hand the
    next stage 101-based values for a 0-based list."""
    entry = {
        "receptor_position_order": [101, 102, 103, 104, 105],
        "per_design_design_residues": [102, 104],
    }
    assert ddr.derive_positional_indices(entry) == [1, 3]


@pytest.mark.logic
def test_the_transform_matches_its_sibling_on_the_same_walk():
    """The two modules differ only in which field they read. Given the same
    position order and the same residue list, they must agree, because the
    design region and the protected set partition the same receptor."""
    order = list(range(20, 30))
    residues = [21, 25, 29]
    a = ddr.derive_positional_indices(
        {"receptor_position_order": order,
         "per_design_design_residues": residues})
    b = dti.derive_positional_indices(
        {"receptor_position_order": order,
         "receptor_contact_residues": residues})
    assert a == b == [1, 5, 9]


@pytest.mark.logic
def test_indices_come_back_sorted():
    entry = {
        "receptor_position_order": [5, 6, 7, 8],
        "per_design_design_residues": [8, 5, 7],
    }
    assert ddr.derive_positional_indices(entry) == [0, 2, 3]


@pytest.mark.logic
@pytest.mark.parametrize("entry,missing", [
    ({"receptor_position_order": [1]}, "per_design_design_residues"),
    ({"per_design_design_residues": [1]}, "receptor_position_order"),
])
def test_a_missing_field_names_the_field(entry, missing):
    with pytest.raises(ValueError) as exc:
        ddr.derive_positional_indices(entry)
    assert missing in str(exc.value)


@pytest.mark.logic
def test_an_empty_design_region_is_an_error():
    """An empty design region means nothing was rebuilt. Writing it as an
    empty file would let a steering run proceed with no region at all."""
    entry = {"receptor_position_order": [1, 2],
             "per_design_design_residues": []}
    with pytest.raises(ValueError):
        ddr.derive_positional_indices(entry)


# ── Output format, which is 1-based ──────────────────────────────────

@pytest.mark.logic
def test_output_declares_the_one_based_convention_in_its_header():
    """The file outlives the run that wrote it. Without the convention in the
    header a reader cannot tell it from the 0-based file the sibling writes."""
    text = ddr.format_output([1, 2, 3], Path("m.json"), "design_1", ranges=True)
    assert "1-based" in text


@pytest.mark.logic
def test_ranges_collapse_contiguous_runs():
    text = ddr.format_output([30, 31, 32, 35, 36, 37], Path("m.json"),
                             "design_1", ranges=True)
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert "30-32" in body[0] and "35-37" in body[0]


@pytest.mark.logic
def test_plain_list_is_a_comma_separated_run():
    text = ddr.format_output([4, 5, 9], Path("m.json"), "design_1", ranges=False)
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert body[0] == "4,5,9"


@pytest.mark.logic
def test_a_single_position_is_not_written_as_a_range():
    text = ddr.format_output([7], Path("m.json"), "design_1", ranges=True)
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert body[0] == "7"


@pytest.mark.logic
def test_output_carries_its_provenance():
    text = ddr.format_output([1], Path("/tmp/rfd.json"), "design_42", ranges=True)
    assert "design_42" in text and "rfd.json" in text


# ── The seam, stated as a test ───────────────────────────────────────

@pytest.mark.logic
def test_derivation_is_zero_based_and_output_is_one_based():
    """The conversion between them is a single +1 in main(), and it is the
    thing most likely to be dropped in a refactor. Index 0 is position 1.
    """
    entry = {"receptor_position_order": [10, 11, 12],
             "per_design_design_residues": [10]}
    zero_based = ddr.derive_positional_indices(entry)
    assert zero_based == [0]

    one_based = [i + 1 for i in zero_based]
    text = ddr.format_output(one_based, Path("m.json"), "d", ranges=False)
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert body[0] == "1"


@pytest.mark.logic
def test_end_to_end_on_the_metrics_shape_the_stage_reads(tmp_path):
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text(json.dumps({"designs": [{
        "design": "design_3.pdb",
        "receptor_position_order": list(range(1, 11)),
        "per_design_design_residues": [4, 5, 6],
    }]}))
    entry = ddr.find_design_entry(json.loads(m.read_text()), "3")
    assert ddr.derive_positional_indices(entry) == [3, 4, 5]
