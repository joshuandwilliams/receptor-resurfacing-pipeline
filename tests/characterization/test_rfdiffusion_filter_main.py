"""Drive the RFdiffusion filter stage end to end.

`main()` is the stage. It reads the design PDBs, resolves the contig against
the input structure, scores every design, and writes the three files the rest
of the pipeline reads:

    rfdiffusion_metrics.json   per-design metrics, and the position order
                               every later positional index is counted against
    passing_designs.txt        which designs go forward
    filter_summary.json        the counts the report shows

The tests here assert the CONTENTS of those files, not that they exist. A
stage that writes three well-formed files describing the wrong designs is the
failure mode worth catching, because nothing downstream can tell.
"""
from __future__ import annotations

import json
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


def _stage(tmp_path: Path, n_designs: int = 3, *, separation: float = 4.0,
           n_receptor: int = 20, n_effector: int = 6):
    """An input complex plus `n_designs` design PDBs, in the stage's layout."""
    inp = complex_pdb(tmp_path / "input.pdb", n_receptor=n_receptor,
                      n_effector=n_effector, separation=separation,
                      receptor_start=1)
    designs = tmp_path / "designs"
    designs.mkdir(exist_ok=True)
    for i in range(n_designs):
        complex_pdb(designs / f"design_{i}.pdb", n_receptor=n_receptor,
                    n_effector=n_effector, separation=separation,
                    single_chain=True)
    return inp, designs


def _run(tmp_path, inp, designs, contigs="A1-8/4-4/A13-20 B", **extra):
    argv = ["rfdiffusion_filter.py",
            "--input-pdb", str(inp),
            "--design-dir", str(designs),
            "--contigs", contigs]
    for k, v in extra.items():
        argv += [f"--{k.replace('_', '-')}", str(v)]
    cwd = Path.cwd()
    import os
    os.chdir(tmp_path)
    old = sys.argv
    try:
        sys.argv = argv
        try:
            rf.main()
        except SystemExit as exc:
            return int(exc.code or 0)
        return 0
    finally:
        sys.argv = old
        os.chdir(cwd)


def _outputs(tmp_path):
    metrics = json.loads((tmp_path / "rfdiffusion_metrics.json").read_text())
    summary = json.loads((tmp_path / "filter_summary.json").read_text())
    passing = [ln for ln in (tmp_path / "passing_designs.txt")
               .read_text().splitlines() if ln.strip()]
    return metrics, summary, passing


# ── The empty case ───────────────────────────────────────────────────

@pytest.mark.logic
def test_no_designs_still_writes_all_three_files(tmp_path):
    """Nextflow publishes a fixed file list. Exiting without writing would
    fail the stage with nothing in the log explaining why, so the empty case
    writes well-formed empty files and exits 0.
    """
    inp = complex_pdb(tmp_path / "input.pdb", n_receptor=20, n_effector=6)
    empty = tmp_path / "designs"
    empty.mkdir()
    assert _run(tmp_path, inp, empty) == 0

    metrics, summary, passing = _outputs(tmp_path)
    assert metrics == {"designs": [], "fixed_residues": []}
    assert summary["n_total"] == 0 and summary["n_passing"] == 0
    assert passing == []


# ── The normal case ──────────────────────────────────────────────────

@pytest.mark.logic
def test_every_design_appears_in_the_metrics(tmp_path):
    inp, designs = _stage(tmp_path, n_designs=3)
    _run(tmp_path, inp, designs)
    metrics, summary, _passing = _outputs(tmp_path)
    assert len(metrics["designs"]) == 3
    assert summary["n_total"] == 3


@pytest.mark.logic
def test_each_entry_carries_the_position_order_later_stages_index_against(tmp_path):
    """receptor_position_order is what derive_design_region and
    derive_true_interface count their positional indices against. An entry
    without it makes every later index meaningless."""
    inp, designs = _stage(tmp_path, n_designs=2)
    _run(tmp_path, inp, designs)
    metrics, _summary, _passing = _outputs(tmp_path)
    for entry in metrics["designs"]:
        assert entry["receptor_position_order"], entry.get("design")
        assert len(entry["receptor_position_order"]) == 20


@pytest.mark.logic
def test_the_counts_in_the_summary_agree_with_the_passing_file(tmp_path):
    """Two files reporting the same thing must not disagree. A mismatch would
    make the report say one number and the next stage receive another."""
    inp, designs = _stage(tmp_path, n_designs=4)
    _run(tmp_path, inp, designs)
    _metrics, summary, passing = _outputs(tmp_path)
    assert summary["n_passing"] == len(passing)
    assert summary["n_total"] == summary["n_passing"] + summary["n_filtered"]


@pytest.mark.logic
def test_design_names_are_stems_not_paths(tmp_path):
    """passing_designs.txt is read by name downstream. A full path there would
    not match."""
    inp, designs = _stage(tmp_path, n_designs=2)
    _run(tmp_path, inp, designs)
    _m, _s, passing = _outputs(tmp_path)
    for name in passing:
        assert "/" not in name
        assert name.startswith("design_")


# ── The gate ─────────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_zero_minimum_passes_everything(tmp_path):
    """min_hotspot_frac defaults to 0.0, which is why the gate never fires on
    the committed campaign. The name is also wrong: the fraction is of all
    receptor-effector contacts falling inside the design region, and the
    hotspot list is parsed but not used."""
    inp, designs = _stage(tmp_path, n_designs=3)
    _run(tmp_path, inp, designs, min_hotspot_frac=0.0)
    _m, summary, passing = _outputs(tmp_path)
    assert summary["n_passing"] == 3
    assert len(passing) == 3


@pytest.mark.logic
def test_an_unreachable_minimum_filters_everything(tmp_path):
    """A fraction above 1.0 cannot be met, so every design must be filtered
    and the files must still be written."""
    inp, designs = _stage(tmp_path, n_designs=3)
    _run(tmp_path, inp, designs, min_hotspot_frac=1.5)
    _m, summary, passing = _outputs(tmp_path)
    assert summary["n_passing"] == 0
    assert passing == []
    assert summary["n_filtered"] == 3


@pytest.mark.logic
def test_the_minimum_is_recorded_in_the_summary(tmp_path):
    """The threshold that produced the counts has to travel with them, or a
    rerun at a different setting is indistinguishable."""
    inp, designs = _stage(tmp_path, n_designs=2)
    _run(tmp_path, inp, designs, min_hotspot_frac=0.25)
    _m, summary, _p = _outputs(tmp_path)
    assert summary["min_hotspot_frac"] == pytest.approx(0.25)


@pytest.mark.logic
def test_a_separated_complex_has_no_contacts_and_is_filtered(tmp_path):
    """No interface means the fraction is undefined. The design must not pass
    a gate it cannot be measured against."""
    inp, designs = _stage(tmp_path, n_designs=2, separation=60.0)
    _run(tmp_path, inp, designs, min_hotspot_frac=0.1)
    _m, summary, _p = _outputs(tmp_path)
    assert summary["n_passing"] == 0


# ── Options reaching the computation ─────────────────────────────────

@pytest.mark.logic
def test_the_contact_cutoff_reaches_the_contact_count(tmp_path):
    """A wider cutoff can only find more contacts. If the flag did not reach
    the computation both runs would report the same number."""
    inp, designs = _stage(tmp_path, n_designs=1, separation=4.0)
    _run(tmp_path, inp, designs, rfdiff_contact_cutoff=5.0)
    tight = _outputs(tmp_path)[0]["designs"][0]
    _run(tmp_path, inp, designs, rfdiff_contact_cutoff=20.0)
    wide = _outputs(tmp_path)[0]["designs"][0]
    assert wide["n_contact_pairs"] >= tight["n_contact_pairs"]
    assert wide["n_contact_pairs"] > 0


@pytest.mark.logic
def test_a_hotspot_list_is_accepted_and_does_not_crash(tmp_path):
    """hotspot_residues is parsed and written into the metrics but does not
    gate, which notes/things_to_fix.md records under 'min_hotspot_frac does
    not read the hotspot list'."""
    inp, designs = _stage(tmp_path, n_designs=2)
    assert _run(tmp_path, inp, designs, hotspot="B2,B3") == 0
    _m, summary, _p = _outputs(tmp_path)
    assert summary["n_total"] == 2


@pytest.mark.logic
def test_the_fixed_residue_set_is_written(tmp_path):
    """Downstream calls a residue 'design' when it is absent from this set,
    so an empty one would classify everything as rebuilt."""
    inp, designs = _stage(tmp_path, n_designs=1)
    _run(tmp_path, inp, designs, contigs="A1-8/4-4/A13-20 B")
    metrics, _s, _p = _outputs(tmp_path)
    assert metrics["fixed_residues"], "the retained residues must be recorded"
