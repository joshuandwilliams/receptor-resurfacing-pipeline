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


# ── main(), and the file it leaves behind ────────────────────────────
#
# Everything above pins a helper. What a steering run actually consumes is the
# file main() writes, so the block below drives main() on input where every
# value is known and asserts the bytes.

def _metrics_file(tmp_path, *, design="design_3.pdb", order=None,
                  design_residues=None, name="rfdiffusion_metrics.json"):
    """Receptor residues 101..110 with the design region at 104..106, unless
    overridden. Positional indices 3, 4, 5; 1-based positions 4, 5, 6."""
    order = list(range(101, 111)) if order is None else order
    design_residues = ([104, 105, 106] if design_residues is None
                       else design_residues)
    p = tmp_path / name
    p.write_text(json.dumps({"designs": [{
        "design": design,
        "receptor_position_order": order,
        "per_design_design_residues": design_residues,
    }]}))
    return p


def _run(monkeypatch, *argv) -> int:
    monkeypatch.setattr(sys, "argv", ["derive_design_region.py", *argv])
    return ddr.main()


@pytest.mark.logic
def test_main_writes_the_whole_file_byte_for_byte(tmp_path, monkeypatch):
    """The receptor is numbered 101..110 and the rebuilt residues are
    104..106, so the positional indices are 3, 4, 5 and the file must carry
    them as the 1-based range 4-6. Asserting the whole text rather than one
    line keeps the header, the convention line and the body pinned together,
    because a reader months later has only this file to go on.
    """
    m = _metrics_file(tmp_path)
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m),
                "--design", "design_3", "--output", str(out)) == 0
    assert out.read_text() == (
        f"# derived from {m}\n"
        "# design: design_3\n"
        "# source field: per_design_design_residues (3 residues)\n"
        "# coordinate system: 1-based positional indices on the receptor chain,\n"
        "#                    matching prediction PDB sequence order\n"
        "4-6\n"
    )


@pytest.mark.logic
def test_main_writes_a_plain_list_under_no_ranges(tmp_path, monkeypatch):
    m = _metrics_file(tmp_path)
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "3",
                "--output", str(out), "--no-ranges") == 0
    assert out.read_text().splitlines()[-1] == "4,5,6"


@pytest.mark.logic
def test_main_reports_a_gapped_region_as_separate_runs(tmp_path, monkeypatch,
                                                       capsys):
    """Design residues 102, 104 and 105 give positional indices 1, 3, 4, so
    1-based 2, 4, 5. The region is not contiguous and the summary must say
    so, because a gapped design region usually means the contig and the PDB
    disagree."""
    m = _metrics_file(tmp_path, design_residues=[102, 104, 105])
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 0
    assert out.read_text().splitlines()[-1] == "2,4-5"
    err = capsys.readouterr().err
    assert "with gaps" in err
    assert "1-based range 2..5" in err


@pytest.mark.logic
def test_the_header_records_the_design_id_as_it_was_typed(tmp_path,
                                                          monkeypatch):
    """The lookup normalises 'design_3', 'design_3.pdb' and '3' to one key,
    but the provenance line does not: it records the argument. That is the
    behaviour, and it matters because the header is the only link back to the
    command that produced the file."""
    m = _metrics_file(tmp_path)
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "3",
                "--output", str(out)) == 0
    assert "# design: 3\n" in out.read_text()


@pytest.mark.logic
def test_the_header_names_the_source_metrics_path(tmp_path, monkeypatch):
    """The file outlives the run. Without the source path there is no way to
    tell which rfdiffusion_metrics.json it came from."""
    m = _metrics_file(tmp_path, name="run_7_metrics.json")
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 0
    assert str(m) in out.read_text().splitlines()[0]
    assert "run_7_metrics.json" in out.read_text()


@pytest.mark.logic
def test_the_count_in_the_header_is_the_derived_count_not_the_source_count(
        tmp_path, monkeypatch):
    """`# source field: per_design_design_residues (N residues)` reports the
    number of positions written, not the length of the field it names. The
    two are equal whenever the derivation succeeds, because a residue absent
    from the position order aborts the run, so the wording is loose rather
    than wrong. The sibling module passes the source length explicitly
    instead; the two headers are not built the same way.
    """
    m = _metrics_file(tmp_path, design_residues=[104, 105])
    out = tmp_path / "region.txt"
    _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
         "--output", str(out))
    assert "(2 residues)" in out.read_text()


@pytest.mark.logic
def test_an_existing_output_file_is_overwritten(tmp_path, monkeypatch):
    """Documented in --output's help. A stale region silently left in place
    would steer the next cycle against the previous design."""
    m = _metrics_file(tmp_path)
    out = tmp_path / "region.txt"
    out.write_text("99-200\n")
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 0
    assert "99-200" not in out.read_text()


# ── Failure paths. Each must name its cause and write nothing ────────

@pytest.mark.logic
def test_a_missing_metrics_file_fails_naming_the_path(tmp_path, monkeypatch,
                                                      capsys):
    out = tmp_path / "region.txt"
    missing = tmp_path / "absent.json"
    assert _run(monkeypatch, "--metrics-json", str(missing),
                "--design", "design_3", "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "not found" in err and str(missing) in err
    assert not out.exists(), "a failed derivation must leave no output file"


@pytest.mark.logic
def test_unparseable_json_fails_naming_the_file(tmp_path, monkeypatch, capsys):
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text("{not json")
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 2
    assert "failed to parse" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.logic
def test_an_absent_design_fails_listing_the_designs_that_exist(
        tmp_path, monkeypatch, capsys):
    m = _metrics_file(tmp_path, design="design_3.pdb")
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_9",
                "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "design_9" in err and "design_3.pdb" in err
    assert not out.exists()


@pytest.mark.logic
def test_metrics_without_a_designs_array_fails_naming_the_type(
        tmp_path, monkeypatch, capsys):
    m = tmp_path / "rfdiffusion_metrics.json"
    m.write_text(json.dumps({"designs": {"design_3": {}}}))
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 2
    assert "dict" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.logic
def test_an_empty_design_region_fails_rather_than_writing_an_empty_file(
        tmp_path, monkeypatch, capsys):
    """An empty region file would read as a valid instruction to steer
    nothing, and the run would complete reporting no reversion. The failure
    has to be loud."""
    m = _metrics_file(tmp_path, design_residues=[])
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 2
    assert "per_design_design_residues' is empty" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.logic
def test_a_design_residue_absent_from_the_position_order_fails(
        tmp_path, monkeypatch, capsys):
    """The two fields are separate walks over the same structure. If one
    holds a residue the other does not they describe different structures,
    and the derivation must abort rather than emit a shorter region."""
    m = _metrics_file(tmp_path, design_residues=[104, 999])
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out)) == 2
    err = capsys.readouterr().err
    assert "999" in err and "atom walks" in err
    assert not out.exists()


# ── The advisory cross-check against plan.json ───────────────────────
#
# The cross-check compares in the 0-based frame, while the file main() writes
# is 1-based. Both appear below on the same input so the shift is explicit.

@pytest.mark.logic
def test_cross_check_compares_in_the_zero_based_frame(tmp_path, monkeypatch,
                                                      capsys):
    """Design residues 104..106 are positional indices 3, 4, 5. plan.json's
    true_interface_idx is 0-based, so [3, 5] is inside, while the file written
    for the same design says 4-6."""
    m = _metrics_file(tmp_path)
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"true_interface_idx": [3, 5]}))
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out), "--cross-check-plan", str(plan)) == 0
    err = capsys.readouterr().err
    assert "OK: true interface is a subset" in err
    assert "2 inside derived design region, 0 outside" in err
    assert out.read_text().splitlines()[-1] == "4-6"


@pytest.mark.logic
@pytest.mark.parametrize("idx,inside", [(3, True), (2, False), (5, True),
                                        (6, False)])
def test_cross_check_boundary_is_the_first_and_last_derived_index(
        tmp_path, idx, inside):
    """Indices 3, 4 and 5 are the region. 2 and 6 are the positions either
    side of it, and each must fall outside."""
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"true_interface_idx": [idx]}))
    assert ddr.cross_check_with_plan(plan, [3, 4, 5]) is inside


@pytest.mark.logic
def test_cross_check_warns_but_still_writes_the_output(tmp_path, monkeypatch,
                                                       capsys):
    """Advisory means advisory. A mismatch is surfaced and the exit code and
    the file are both unchanged, because a native binding site reaching
    outside the de novo block is legal if unusual."""
    m = _metrics_file(tmp_path)
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"true_interface_idx": [3, 9]}))
    out = tmp_path / "region.txt"
    assert _run(monkeypatch, "--metrics-json", str(m), "--design", "design_3",
                "--output", str(out), "--cross-check-plan", str(plan)) == 0
    err = capsys.readouterr().err
    assert "WARNING" in err and "[9]" in err
    assert out.read_text().splitlines()[-1] == "4-6"


@pytest.mark.logic
def test_cross_check_skips_an_absent_plan(tmp_path, capsys):
    assert ddr.cross_check_with_plan(tmp_path / "nope.json", [0, 1]) is True
    assert "skipping" in capsys.readouterr().err


@pytest.mark.logic
def test_cross_check_skips_an_unreadable_plan(tmp_path, capsys):
    bad = tmp_path / "plan.json"
    bad.write_text("{not json")
    assert ddr.cross_check_with_plan(bad, [0, 1]) is True
    assert "failed to read plan.json" in capsys.readouterr().err


@pytest.mark.logic
def test_cross_check_skips_a_plan_without_a_true_interface(tmp_path, capsys):
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"true_interface_idx": []}))
    assert ddr.cross_check_with_plan(plan, [0, 1]) is True
    assert "no true_interface_idx" in capsys.readouterr().err
