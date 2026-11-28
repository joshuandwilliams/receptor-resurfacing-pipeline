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


# ── main(), and the file it leaves behind ────────────────────────────
#
# The helpers above are pinned. What the steering engine reads is the file
# main() writes, so the block below drives main() on input where every value
# is known and asserts the bytes, including the frame the header declares.

def _metrics_file(tmp_path, *, design="design_12.pdb", order=None,
                  contacts=None, name="rfdiffusion_metrics.json"):
    """Receptor residues 50..59 with contacts at 51, 55 and 59, unless
    overridden. The positional indices are therefore 1, 5 and 9."""
    order = list(range(50, 60)) if order is None else order
    contacts = [51, 55, 59] if contacts is None else contacts
    p = tmp_path / name
    p.write_text(json.dumps({"designs": [{
        "design": design,
        "receptor_position_order": order,
        "receptor_contact_residues": contacts,
    }]}))
    return p


def _run(monkeypatch, *argv) -> int:
    monkeypatch.setattr(sys, "argv", ["derive_true_interface.py", *argv])
    return dti.main()


@pytest.mark.logic
def test_main_writes_the_whole_file_byte_for_byte(tmp_path, monkeypatch):
    """The receptor is numbered 50..59 and the contacting residues are 51, 55
    and 59, so the positional indices are 1, 5 and 9. One per line, 0-based,
    with the convention stated in the header. Asserting the whole text keeps
    the provenance, the frame declaration and the body pinned together: this
    file is the protected set, and a reader has nothing else to go on.
    """
    m = _metrics_file(tmp_path)
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 0
    assert out.read_text() == (
        f"# derived from {m}\n"
        "# design: design_12\n"
        "# source field: receptor_contact_residues (3 residues)\n"
        "# coordinate system: 0-based positional indices on the receptor chain,\n"
        "#                    matching prediction PDB sequence order\n"
        "1\n5\n9\n"
    )


@pytest.mark.logic
def test_the_written_indices_are_positional_not_residue_numbers(
        tmp_path, monkeypatch):
    """The receptor starts at 50. Writing residue numbers would hand the
    engine 50-based values for a 0-based protected set, and every one of them
    would protect the wrong position. The file must hold 1, 5, 9 and never
    51, 55, 59."""
    m = _metrics_file(tmp_path)
    out = tmp_path / "iface.txt"
    _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
         "--output", str(out))
    body = [ln for ln in out.read_text().splitlines()
            if ln and not ln.startswith("#")]
    assert body == ["1", "5", "9"]


@pytest.mark.logic
def test_the_first_receptor_position_is_written_as_zero(tmp_path, monkeypatch):
    """The boundary of the 0-based frame. A contact on the receptor's first
    residue is index 0, not 1."""
    m = _metrics_file(tmp_path, contacts=[50])
    out = tmp_path / "iface.txt"
    _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
         "--output", str(out))
    assert out.read_text().splitlines()[-1] == "0"


@pytest.mark.logic
def test_the_header_names_the_source_and_the_design(tmp_path, monkeypatch):
    m = _metrics_file(tmp_path, name="campaign_metrics.json")
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "12",
                "--output", str(out)) == 0
    text = out.read_text()
    assert str(m) in text and "campaign_metrics.json" in text
    # The id is recorded as typed, not normalised to 'design_12'.
    assert "# design: 12\n" in text


@pytest.mark.logic
def test_the_header_count_is_the_source_field_length(tmp_path, monkeypatch):
    """Unlike derive_design_region, this module passes the length of
    receptor_contact_residues explicitly, so the count names the source
    field rather than the number of lines written."""
    m = _metrics_file(tmp_path, contacts=[51, 55])
    out = tmp_path / "iface.txt"
    _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
         "--output", str(out))
    assert "receptor_contact_residues (2 residues)" in out.read_text()


@pytest.mark.logic
def test_main_reports_the_derived_span_on_stderr(tmp_path, monkeypatch,
                                                 capsys):
    m = _metrics_file(tmp_path)
    out = tmp_path / "iface.txt"
    _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
         "--output", str(out))
    err = capsys.readouterr().err
    assert "derived true interface: 3 positions (0-based range 1..9)" in err
    assert str(out) in err


@pytest.mark.logic
def test_an_existing_output_file_is_overwritten(tmp_path, monkeypatch):
    m = _metrics_file(tmp_path)
    out = tmp_path / "iface.txt"
    out.write_text("777\n")
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 0
    assert "777" not in out.read_text()


# ── Failure paths. Each must name its cause and write nothing ────────

@pytest.mark.logic
def test_a_missing_metrics_file_fails_naming_the_path(tmp_path, monkeypatch,
                                                      capsys):
    missing = tmp_path / "absent.json"
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(missing),
                "--design", "design_12", "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "not found" in err and str(missing) in err
    assert not out.exists(), "a failed derivation must leave no protected set"


@pytest.mark.logic
def test_unparseable_json_fails_naming_the_file(tmp_path, monkeypatch, capsys):
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text("]]not json[[")
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 2
    assert "failed to parse" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.logic
def test_an_absent_design_fails_listing_the_designs_that_exist(
        tmp_path, monkeypatch, capsys):
    m = _metrics_file(tmp_path)
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_99",
                "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "design_99" in err and "design_12.pdb" in err
    assert not out.exists()


@pytest.mark.logic
def test_metrics_without_a_designs_array_fails_naming_the_type(
        tmp_path, monkeypatch, capsys):
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text(json.dumps({"designs": "design_12"}))
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 2
    assert "got type str" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.logic
def test_no_contacts_fails_rather_than_writing_an_empty_protected_set(
        tmp_path, monkeypatch, capsys):
    """An empty protected-set file would read as a valid instruction that
    nothing is protected, and negative steering would be free to mutate the
    binding site itself. That is the one thing the method must never do, so
    the stage fails instead of writing the file."""
    m = _metrics_file(tmp_path, contacts=[])
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "receptor_contact_residues' is empty" in err
    assert "not a binder" in err
    assert not out.exists()


@pytest.mark.logic
def test_an_empty_position_order_fails(tmp_path, monkeypatch, capsys):
    m = _metrics_file(tmp_path, order=[], contacts=[51])
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 2
    assert "receptor_position_order' is empty" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.logic
def test_a_contact_absent_from_the_position_order_fails(tmp_path, monkeypatch,
                                                        capsys):
    m = _metrics_file(tmp_path, contacts=[51, 404])
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "404" in err and "atom walks" in err
    assert not out.exists()


# ── The two modules meet at the cross-check ──────────────────────────

@pytest.mark.logic
def test_the_cross_check_reads_a_file_the_sibling_module_wrote(
        tmp_path, monkeypatch, capsys):
    """One pass through both stages on one metrics file, which is the only
    place the 1-based design region and the 0-based interface meet.

    The receptor is numbered 50..59. The design region is residues 53..58,
    which derive_design_region writes as the 1-based range 4-9. The contacts
    are 54 and 57, which this module writes as 0-based 4 and 7. Shifted up
    for the comparison those are positions 5 and 8, both inside 4-9.
    """
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text(json.dumps({"designs": [{
        "design": "design_12.pdb",
        "receptor_position_order": list(range(50, 60)),
        "per_design_design_residues": [53, 54, 55, 56, 57, 58],
        "receptor_contact_residues": [54, 57],
    }]}))
    region = tmp_path / "region.txt"
    monkeypatch.setattr(sys, "argv", [
        "derive_design_region.py", "--metrics-json", str(m),
        "--design", "design_12", "--output", str(region)])
    import derive_design_region as _ddr
    assert _ddr.main() == 0
    assert region.read_text().splitlines()[-1] == "4-9"

    iface = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(iface),
                "--cross-check-design-region", str(region)) == 0
    assert [ln for ln in iface.read_text().splitlines()
            if not ln.startswith("#")] == ["4", "7"]
    err = capsys.readouterr().err
    assert "OK: true interface is a subset" in err


@pytest.mark.logic
@pytest.mark.parametrize("idx,inside", [(2, False), (3, True), (5, True),
                                        (6, False)])
def test_cross_check_boundary_is_the_region_shifted_by_one(tmp_path, idx,
                                                           inside):
    """A design-region file reading 4-6 is 1-based, so it covers 0-based
    indices 3, 4 and 5. Index 2 and index 6 are the positions either side,
    and each must fall outside."""
    f = tmp_path / "dr.txt"
    f.write_text("# derived from somewhere\n4-6\n")
    assert dti.cross_check_design_region(f, [idx]) is inside


@pytest.mark.logic
def test_cross_check_warns_but_leaves_the_output_unchanged(tmp_path,
                                                           monkeypatch,
                                                           capsys):
    """Advisory means advisory. The exit code stays 0 and the file still
    holds the derived indices, because a contact reaching into a fixed
    residue is a finding to report rather than a reason to lose the run."""
    m = _metrics_file(tmp_path)
    region = tmp_path / "dr.txt"
    region.write_text("2-2\n")
    out = tmp_path / "iface.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_12",
                "--output", str(out),
                "--cross-check-design-region", str(region)) == 0
    err = capsys.readouterr().err
    assert "WARNING" in err and "OUTSIDE" in err
    assert [ln for ln in out.read_text().splitlines()
            if not ln.startswith("#")] == ["1", "5", "9"]


@pytest.mark.logic
def test_cross_check_skips_a_file_with_no_parseable_positions(tmp_path,
                                                              capsys):
    """A header-only file is not a region of zero positions. Treating it as
    one would report every interface residue as outside the region, which
    reads as a coordinate-system bug."""
    f = tmp_path / "dr.txt"
    f.write_text("# derived from nowhere\n# design: d\n")
    assert dti.cross_check_design_region(f, [0, 1]) is True
    assert "no parseable positions" in capsys.readouterr().err


@pytest.mark.logic
def test_cross_check_ignores_tokens_that_are_not_positions(tmp_path):
    """A malformed token is skipped rather than raising, so one stray word
    does not turn an advisory check into a lost run. '4-6' still parses, so
    index 3 is still inside."""
    f = tmp_path / "dr.txt"
    f.write_text("4-6,junk,a-b,\n")
    assert dti.cross_check_design_region(f, [3]) is True


@pytest.mark.logic
def test_cross_check_skips_a_file_it_cannot_read(tmp_path, capsys):
    """A directory where a file was expected exists but cannot be read. The
    check reports and returns rather than taking the stage down."""
    d = tmp_path / "dr.txt"
    d.mkdir()
    assert dti.cross_check_design_region(d, [0]) is True
    assert "failed to read" in capsys.readouterr().err
