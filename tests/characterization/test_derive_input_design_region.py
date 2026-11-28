"""Cover derivation of the design region and true interface from the input.

These two files are the contract between the pipeline and the steering engine.
The design region says which receptor positions were rebuilt. The true
interface says which must never be mutated. An index one out in either tells
the engine to mutate a residue it was told to keep, and the run completes.

The two use different frames on purpose, which is the thing most likely to be
collapsed by a refactor:

    design region    1-based positions on the receptor
    true interface   0-based positional indices

Both are positional within the receptor chain rather than PDB residue numbers,
so a structure that does not start at residue 1 still lines up.

`main()` is the stage the Task 6 negative-control path runs, and the last
section drives it end to end. The two files are the stage's only output, so
the tests there assert what is written into them rather than that they appeared.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import derive_input_design_region as didr  # noqa: E402
from fixtures_runtree import complex_pdb, sequence_pdb  # noqa: E402

# ── Contig parsing to the flat form ──────────────────────────────────

@pytest.mark.logic
def test_anchors_and_denovo_segments_come_back_distinguishable():
    """A de novo segment is marked with the sentinel chain id "_denovo" and
    carries its length bounds where an anchor carries residue numbers.
    Confusing the two would treat a length as a residue range."""
    segs = didr._parse_contigs("A1-10/5-15/A20-30 B", "A", "B")
    assert ("A", 1, 10, False) in segs
    assert ("A", 20, 30, False) in segs
    denovo = [s for s in segs if s[0] == "_denovo"]
    assert denovo == [("_denovo", 5, 15, True)]


@pytest.mark.logic
def test_a_bare_effector_chain_parses_to_a_zero_range():
    """`B` alone means carry the chain through. It has no residue range, so
    it is recorded as 0-0 rather than being dropped."""
    segs = didr._parse_contigs("A1-10/5-5/A16-25 B", "A", "B")
    assert ("B", 0, 0, False) in segs


@pytest.mark.logic
def test_commas_are_tolerated_but_are_not_the_canonical_form():
    """The canonical grammar is '/' within a chain and a space between
    chains. Commas are accepted for back-compatibility only, and this pins
    that they parse to the same thing rather than silently to something
    else."""
    canonical = didr._parse_contigs("A1-10/5-5/A16-25 B", "A", "B")
    with_commas = didr._parse_contigs("A1-10/5-5/A16-25,B", "A", "B")
    assert canonical == with_commas


# ── Receptor anchor walk ─────────────────────────────────────────────

@pytest.mark.logic
def test_only_receptor_anchors_are_walked():
    """The effector's own block must not contribute anchors, or the design
    region would be computed against the wrong chain."""
    segs = didr._parse_contigs("A1-10/5-5/A16-25 B", "A", "B")
    anchors = didr._receptor_positions_by_contigs(segs, "A")
    assert anchors == [(1, 10), (16, 25)]


@pytest.mark.logic
def test_an_all_fixed_receptor_gives_one_anchor():
    segs = didr._parse_contigs("A1-30 B", "A", "B")
    assert didr._receptor_positions_by_contigs(segs, "A") == [(1, 30)]


# ── Design region positions ──────────────────────────────────────────

@pytest.mark.logic
def test_the_documented_example_reproduces():
    """The docstring states: anchors [(1,50), (60,100)] on a 100-residue
    receptor give positions 51..59."""
    got = didr._design_region_positions_1b([(1, 50), (60, 100)], 100)
    assert got == list(range(51, 60))


@pytest.mark.logic
def test_positions_are_one_based():
    """Position 1 is the first residue. Returning 0 would shift the whole
    region by one against the file the engine reads."""
    got = didr._design_region_positions_1b([(4, 10)], 10)
    assert got == [1, 2, 3]
    assert min(got) == 1


@pytest.mark.logic
def test_a_gap_at_each_end_is_found():
    got = didr._design_region_positions_1b([(4, 7)], 10)
    assert got == [1, 2, 3, 8, 9, 10]


@pytest.mark.logic
def test_two_gaps_both_appear():
    got = didr._design_region_positions_1b([(1, 3), (7, 8)], 10)
    assert got == [4, 5, 6, 9, 10]


@pytest.mark.logic
def test_an_all_anchor_receptor_has_no_design_region():
    assert didr._design_region_positions_1b([(1, 10)], 10) == []


@pytest.mark.logic
def test_no_anchors_means_the_whole_receptor_is_rebuilt():
    """Documented edge case. An empty list here would mean nothing is
    rebuilt, which is the opposite of what no anchors means."""
    assert didr._design_region_positions_1b([], 5) == [1, 2, 3, 4, 5]


@pytest.mark.logic
def test_overlapping_anchors_do_not_double_count():
    got = didr._design_region_positions_1b([(1, 5), (3, 7)], 10)
    assert got == [8, 9, 10]


# ── True interface ───────────────────────────────────────────────────

@pytest.mark.logic
def test_the_interface_is_zero_based_positional(tmp_path):
    """Effector residue j sits 4 A above receptor residue j, so a 5 A cutoff
    gives the directly-opposite residues. The receptor starts at 101 here,
    and the answer must still be 0, 1, 2 rather than 101, 102, 103.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=3,
                      spacing=3.8, separation=4.0, receptor_start=101)
    seq = didr._extract_receptor_sequence(pdb, "A")
    idx = didr._compute_true_interface_0b(pdb, "A", "B", 5.0, seq)
    assert idx == [0, 1, 2]


@pytest.mark.logic
def test_the_interface_widens_with_the_cutoff(tmp_path):
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=12, n_effector=4,
                      separation=4.0)
    seq = didr._extract_receptor_sequence(pdb, "A")
    tight = set(didr._compute_true_interface_0b(pdb, "A", "B", 4.5, seq))
    wide = set(didr._compute_true_interface_0b(pdb, "A", "B", 9.0, seq))
    assert tight <= wide


@pytest.mark.logic
def test_a_separated_complex_has_an_empty_interface(tmp_path):
    pdb = complex_pdb(tmp_path / "far.pdb", separation=50.0)
    seq = didr._extract_receptor_sequence(pdb, "A")
    assert didr._compute_true_interface_0b(pdb, "A", "B", 5.0, seq) == []


@pytest.mark.logic
def test_indices_are_sorted_and_unique(tmp_path):
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=12, n_effector=5,
                      separation=4.0)
    seq = didr._extract_receptor_sequence(pdb, "A")
    idx = didr._compute_true_interface_0b(pdb, "A", "B", 8.0, seq)
    assert idx == sorted(set(idx))


# ── The two files, and the frame they disagree on ────────────────────

@pytest.mark.logic
def test_the_design_region_file_carries_its_provenance(tmp_path):
    out = tmp_path / "design_region.txt"
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=2)
    didr._write_design_region(out, [5, 6, 7], "A1-4/3-3/A8-10 B", pdb, "A")
    text = out.read_text()
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert body, "the positions themselves must be written"
    nums = [int(x) for ln in body for x in ln.replace(",", " ").split()]
    assert nums == [5, 6, 7]


@pytest.mark.logic
def test_the_true_interface_file_carries_its_provenance(tmp_path):
    out = tmp_path / "true_interface.txt"
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=2)
    didr._write_true_interface(out, [0, 2, 4], pdb, "A", "B", 5.0)
    body = [ln for ln in out.read_text().splitlines()
            if ln and not ln.startswith("#")]
    nums = [int(x) for ln in body for x in ln.replace(",", " ").split()]
    assert nums == [0, 2, 4]


@pytest.mark.logic
def test_the_two_frames_differ_by_one_on_the_same_residue(tmp_path):
    """The seam, stated as a test. For a receptor where residue 1 is both an
    interface residue and outside every anchor, the true interface records 0
    and the design region records 1.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=1,
                      spacing=3.8, separation=4.0, receptor_start=1)
    seq = didr._extract_receptor_sequence(pdb, "A")
    interface = didr._compute_true_interface_0b(pdb, "A", "B", 5.0, seq)
    assert interface == [0]
    design = didr._design_region_positions_1b([(2, 10)], 10)
    assert design == [1]


# ── The stage, end to end ────────────────────────────────────────────

def _run(tmp_path: Path, pdb: Path, contigs: str, **extra) -> int:
    """Drive `main()` with argv, returning its exit code.

    The two output paths always land in `tmp_path`, so `_files` below can find
    them without each test naming them again.
    """
    argv = ["derive_input_design_region.py",
            "--input-pdb", str(pdb),
            "--contigs", contigs,
            "--design-region-output", str(tmp_path / "input_design_region.txt"),
            "--true-interface-output",
            str(tmp_path / "input_true_interface.txt")]
    for k, v in extra.items():
        argv += [f"--{k.replace('_', '-')}", str(v)]
    old = sys.argv
    try:
        sys.argv = argv
        try:
            return int(didr.main() or 0)
        except SystemExit as exc:
            return int(exc.code or 0)
    finally:
        sys.argv = old


def _files(tmp_path: Path) -> tuple[Path, Path]:
    return (tmp_path / "input_design_region.txt",
            tmp_path / "input_true_interface.txt")


def _numbers(path: Path) -> list[int]:
    """Every index in a written file, whichever of the two formats it is in.

    The design region is one comma-separated line and the true interface is
    one index per line, so reading both the same way keeps a test about the
    indices from also being a test about the separator.
    """
    body = [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")]
    return [int(x) for ln in body for x in ln.replace(",", " ").split()]


def _two_chain_complex(tmp_path: Path, *, separation: float = 4.0) -> Path:
    """A 10-residue receptor with a 3-residue effector lying over its first
    three residues, `separation` angstroms up in z.

    Effector residue j sits directly above receptor residue j, so at a 5 Å
    cutoff the interface is exactly positions 0, 1 and 2: the next receptor
    residue along is sqrt(3.8^2 + 4^2) = 5.52 Å away and misses.
    """
    return complex_pdb(tmp_path / "input.pdb", n_receptor=10, n_effector=3,
                       spacing=3.8, separation=separation)


@pytest.mark.logic
def test_the_two_files_hold_the_indices_the_contigs_and_the_structure_imply(tmp_path):
    """Anchor A1-6 on a 10-residue receptor leaves positions 7..10 to be
    rebuilt, and the effector sits over positions 0..2. Both answers are set
    by the fixture, so a shift in either frame fails here."""
    pdb = _two_chain_complex(tmp_path)
    assert _run(tmp_path, pdb, "A1-6/4-4 B") == 0
    design, interface = _files(tmp_path)
    assert _numbers(design) == [7, 8, 9, 10]
    assert _numbers(interface) == [0, 1, 2]


@pytest.mark.logic
def test_the_headers_record_what_produced_the_numbers(tmp_path):
    """Both files are read by hand when a campaign is debugged, and a size or
    a cutoff that disagrees with the body is how a stale file goes unnoticed.
    """
    pdb = _two_chain_complex(tmp_path)
    _run(tmp_path, pdb, "A1-6/4-4 B", contact_cutoff=5.0)
    design, interface = _files(tmp_path)
    design_header = design.read_text()
    assert "# contigs: A1-6/4-4 B" in design_header
    assert "# design region size: 4 residues" in design_header
    assert "# receptor chain: A" in design_header
    interface_header = interface.read_text()
    assert "# heavy-atom contacts <= 5.0 Angstroms" in interface_header
    assert "# n contacts: 3" in interface_header
    assert "# effector chain: B" in interface_header


@pytest.mark.logic
def test_an_all_anchor_contig_writes_an_empty_design_region(tmp_path):
    """Nothing is rebuilt, so nothing may be scrambled. The file still has to
    be written, because the Nextflow process publishes a fixed file list and
    the control stage downstream reads it for a count of zero rather than
    failing to find it."""
    pdb = _two_chain_complex(tmp_path)
    assert _run(tmp_path, pdb, "A1-10 B") == 0
    design, _interface = _files(tmp_path)
    assert _numbers(design) == []
    assert "# design region size: 0 residues" in design.read_text()


@pytest.mark.logic
def test_the_contact_cutoff_flag_reaches_the_contact_search(tmp_path):
    """At 4.5 Å only the residue directly under each effector residue is in
    contact, at 4.0 Å exactly. At 7.0 Å the next one along, at
    sqrt(3.8^2 + 4^2) = 5.52 Å, joins; the one after that is at
    sqrt(7.6^2 + 4^2) = 8.59 Å and still does not."""
    pdb = _two_chain_complex(tmp_path)
    _run(tmp_path, pdb, "A1-6/4-4 B", contact_cutoff=4.5)
    assert _numbers(_files(tmp_path)[1]) == [0, 1, 2]
    _run(tmp_path, pdb, "A1-6/4-4 B", contact_cutoff=7.0)
    assert _numbers(_files(tmp_path)[1]) == [0, 1, 2, 3]


@pytest.mark.logic
def test_a_separated_complex_writes_an_empty_interface(tmp_path):
    """No interface is a legitimate answer for an input that has not been
    posed yet, and it must be written as zero contacts rather than as a
    crash."""
    pdb = _two_chain_complex(tmp_path, separation=50.0)
    assert _run(tmp_path, pdb, "A1-6/4-4 B") == 0
    _design, interface = _files(tmp_path)
    assert _numbers(interface) == []
    assert "# n contacts: 0" in interface.read_text()


# ── Refusals ─────────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_receptor_chain_that_is_not_in_the_structure_exits_two(tmp_path):
    """Exit 2 and no files. Writing a design region derived from a chain that
    is not there would hand the control stage a set of positions against a
    structure nobody checked."""
    pdb = _two_chain_complex(tmp_path)
    assert _run(tmp_path, pdb, "A1-6/4-4 B", receptor_chain="Z") == 2
    design, interface = _files(tmp_path)
    assert not design.exists() and not interface.exists()


@pytest.mark.logic
def test_an_empty_receptor_chain_exits_two(tmp_path, monkeypatch):
    """The guard is unreachable through the file system today, because
    structure_metrics.get_chain_sequence raises when the chain is absent
    rather than returning an empty string, so the extraction is stubbed. It
    is pinned because a receptor length of zero would make every anchor fall
    outside the chain and the design region come out empty, which is a valid
    answer for a different input."""
    pdb = _two_chain_complex(tmp_path)
    monkeypatch.setattr(didr, "_extract_receptor_sequence", lambda p, c: "")
    assert _run(tmp_path, pdb, "A1-6/4-4 B") == 2
    assert not _files(tmp_path)[0].exists()


@pytest.mark.logic
def test_an_unparseable_contig_exits_two(tmp_path, capsys):
    pdb = _two_chain_complex(tmp_path)
    assert _run(tmp_path, pdb, "A1-6/??/A8-10 B") == 2
    assert "contigs parse failed" in capsys.readouterr().err


@pytest.mark.logic
def test_an_insertion_code_in_the_receptor_exits_two_rather_than_shifting_indices(
        tmp_path, capsys):
    """The sequence walk keys on (resnum, insertion code) and the contact walk
    keys on resnum alone, so a receptor carrying 5 and 5A is one residue
    longer to the first than to the second. Every index from the contact walk
    would then be against a different structure, and the length check in
    find_contact_residues_heavy is what stops that becoming a silent
    off-by-one in the file the steering engine reads.
    """
    pdb = sequence_pdb(tmp_path / "icode.pdb", sequence="MKTAYIAKQR",
                       chain="A", icodes={5: "A"})
    sequence_pdb(pdb, sequence="WVC", chain="B", origin=(0.0, 0.0, 4.0),
                 append=True)
    assert _run(tmp_path, pdb, "A1-6/4-4 B") == 2
    assert "true-interface computation failed" in capsys.readouterr().err


# ── The advisory overlap warning ─────────────────────────────────────

@pytest.mark.logic
def test_rebuilding_over_the_whole_interface_warns(tmp_path, capsys):
    """Anchoring A5-10 leaves positions 1..4 to be rebuilt, which in 0-based
    terms is {0, 1, 2, 3} and covers all three interface positions. The
    controls would then destroy the entire binding face, which the stage
    warns about rather than refusing, because a campaign may intend it."""
    pdb = _two_chain_complex(tmp_path)
    assert _run(tmp_path, pdb, "4-4/A5-10 B") == 0
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "100%" in err


@pytest.mark.logic
def test_an_interface_inside_the_anchors_does_not_warn(tmp_path, capsys):
    """Anchoring A1-6 keeps the whole interface, so the overlap is zero and
    the log must stay quiet. A warning on every run is a warning nobody
    reads."""
    pdb = _two_chain_complex(tmp_path)
    assert _run(tmp_path, pdb, "A1-6/4-4 B") == 0
    assert "WARNING" not in capsys.readouterr().err
