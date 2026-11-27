"""Cover the true-interface derivation, which produces the protected set.

The true interface is the set of receptor positions negative steering is
forbidden to mutate. An index one out tells the engine to mutate a residue it
was told to keep, and the run completes, so nothing reports it.

The index is 0-based and positional within the receptor's own atom walk, not
a residue number. `derive_design_region.py` writes 1-based positions for the
same receptor, which is the boundary the cross-check below crosses.
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

import derive_true_interface as dti  # noqa: E402

# ── Design id normalisation ──────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("raw", ["design_3", "design_3.pdb", "3", " 3 "])
def test_every_spelling_of_a_design_id_normalises_the_same(raw):
    """rfdiffusion_metrics.json has used all of these as its key, so a run
    that fails to match one looks like a missing design rather than a naming
    mismatch."""
    assert dti.normalise_design_id(raw) == {"bare": "design_3",
                                            "with_pdb": "design_3.pdb"}


# ── Finding the entry ────────────────────────────────────────────────

def _metrics(*entries) -> dict:
    return {"designs": list(entries)}


@pytest.mark.logic
@pytest.mark.parametrize("key", ["design_7", "design_7.pdb"])
def test_entry_is_found_under_either_key_spelling(key):
    m = _metrics({"design": key, "receptor_contact_residues": [1]})
    assert dti.find_design_entry(m, "7")["design"] == key


@pytest.mark.logic
def test_entry_is_found_under_the_name_field_too():
    m = _metrics({"name": "design_2", "receptor_contact_residues": [1]})
    assert dti.find_design_entry(m, "design_2")["name"] == "design_2"


@pytest.mark.logic
def test_a_missing_design_names_the_ones_that_exist():
    """The alternative is a bare KeyError, which says nothing about whether
    the design is absent or the id is spelled differently."""
    m = _metrics({"design": "design_1"}, {"design": "design_2"})
    with pytest.raises(KeyError) as exc:
        dti.find_design_entry(m, "design_9")
    msg = str(exc.value)
    assert "design_1" in msg and "design_2" in msg


@pytest.mark.logic
def test_metrics_without_a_designs_array_fails_with_the_type_it_got():
    with pytest.raises(ValueError) as exc:
        dti.find_design_entry({"designs": {"not": "a list"}}, "design_1")
    assert "dict" in str(exc.value)


# ── The derivation itself ────────────────────────────────────────────

@pytest.mark.logic
def test_indices_are_positional_and_zero_based_not_residue_numbers():
    """The receptor starts at residue 101 here. Returning residue numbers
    would hand the engine 101-based values for a 0-based protected set, and
    every one of them would protect the wrong position."""
    entry = {
        "receptor_position_order": [101, 102, 103, 104, 105],
        "receptor_contact_residues": [102, 105],
    }
    assert dti.derive_positional_indices(entry) == [1, 4]


@pytest.mark.logic
def test_indices_come_back_sorted_regardless_of_input_order():
    entry = {
        "receptor_position_order": [10, 20, 30, 40],
        "receptor_contact_residues": [40, 10, 30],
    }
    assert dti.derive_positional_indices(entry) == [0, 2, 3]


@pytest.mark.logic
def test_a_contact_absent_from_the_position_order_aborts():
    """The two fields are separate walks over the structure. If one holds a
    residue the other does not, they describe different structures and any
    index derived from them is meaningless, so this must not degrade to a
    shorter protected set."""
    entry = {
        "receptor_position_order": [1, 2, 3],
        "receptor_contact_residues": [2, 99],
    }
    with pytest.raises(ValueError) as exc:
        dti.derive_positional_indices(entry)
    assert "99" in str(exc.value)


@pytest.mark.logic
@pytest.mark.parametrize("entry,missing", [
    ({"receptor_position_order": [1]}, "receptor_contact_residues"),
    ({"receptor_contact_residues": [1]}, "receptor_position_order"),
])
def test_a_missing_field_names_the_field(entry, missing):
    with pytest.raises(ValueError) as exc:
        dti.derive_positional_indices(entry)
    assert missing in str(exc.value)


@pytest.mark.logic
def test_no_contacts_is_an_error_not_an_empty_protected_set():
    """An empty protected set means steering may mutate the binding site
    itself, which is the one thing the method must never do."""
    entry = {"receptor_position_order": [1, 2], "receptor_contact_residues": []}
    with pytest.raises(ValueError) as exc:
        dti.derive_positional_indices(entry)
    assert "no" in str(exc.value).lower()


@pytest.mark.logic
def test_empty_position_order_is_an_error():
    entry = {"receptor_position_order": [], "receptor_contact_residues": [1]}
    with pytest.raises(ValueError):
        dti.derive_positional_indices(entry)


# ── Cross-check against the design region ────────────────────────────

@pytest.mark.logic
def test_cross_check_crosses_the_zero_to_one_based_boundary(tmp_path):
    """The design-region file is 1-based and the interface indices are
    0-based. Index 0 corresponds to position 1, so a design region of 1-3
    contains interface indices 0, 1 and 2.
    """
    f = tmp_path / "dr.txt"
    f.write_text("# derived\n1-3\n")
    assert dti.cross_check_design_region(f, [0, 1, 2]) is True


@pytest.mark.logic
def test_cross_check_is_advisory_when_the_interface_escapes_the_region(tmp_path):
    """It returns rather than raising, because the resurfacing pattern does
    not require the binding site to sit inside the de novo block on every
    target. The warning is the output."""
    f = tmp_path / "dr.txt"
    f.write_text("1-3\n")
    assert dti.cross_check_design_region(f, [0, 9]) in (True, False)


@pytest.mark.logic
def test_cross_check_skips_quietly_when_the_file_is_absent(tmp_path):
    assert dti.cross_check_design_region(tmp_path / "nope.txt", [0]) is True


@pytest.mark.logic
def test_cross_check_reads_comma_separated_positions(tmp_path):
    f = tmp_path / "dr.txt"
    f.write_text("# header\n1,2,3,4\n")
    assert dti.cross_check_design_region(f, [0, 3]) is True


# ── Output format ────────────────────────────────────────────────────

@pytest.mark.logic
def test_output_carries_its_provenance_in_comments(tmp_path):
    """The file is read months later by a steering run. Without the source
    metrics path and design id in it there is no way to tell which run it
    came from."""
    text = dti.format_output([0, 2, 5], Path("/tmp/m.json"), "design_4", 3)
    assert "design_4" in text
    assert "m.json" in text
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert body, "the indices themselves must be written, not only comments"


@pytest.mark.logic
def test_output_indices_round_trip(tmp_path):
    text = dti.format_output([0, 2, 5], Path("m.json"), "design_1", 3)
    nums = []
    for ln in text.splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        nums += [int(x) for x in ln.replace(",", " ").split()]
    assert nums == [0, 2, 5]


@pytest.mark.logic
def test_a_real_metrics_file_shape_derives_end_to_end(tmp_path):
    """One pass through find -> derive on the JSON shape the stage reads."""
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text(json.dumps({"designs": [{
        "design": "design_12.pdb",
        "receptor_position_order": list(range(50, 60)),
        "receptor_contact_residues": [51, 55, 59],
    }]}))
    entry = dti.find_design_entry(json.loads(m.read_text()), "12")
    assert dti.derive_positional_indices(entry) == [1, 5, 9]
