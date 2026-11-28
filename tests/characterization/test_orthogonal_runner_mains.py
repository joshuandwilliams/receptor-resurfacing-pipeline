"""Drive the two per-survivor metric runners end to end, without the tools.

`tests/characterization/test_orthogonal_metric_runners.py` covers the score-file
parser, the binary lookup and the interface residue set. What it leaves is
`main()` for both stages, which is the part a Nextflow process invokes and the
part that decides what reaches `merge_orthogonal_metrics`.

Neither Rosetta nor FreeSASA is installed here, so the Rosetta stage is driven
by writing executable stubs named like the real binaries into tmp_path and
prepending that directory to PATH. The stub decides its own exit code and
whether it writes an output file, which is what makes each failure row
reachable: a relax that exits non-zero, a relax that exits clean with nothing
to show, an InterfaceAnalyzer that writes no score file, and one that writes a
score file while still failing.

Both stages return 0 and write a row whatever happens. That is deliberate. One
survivor failing must leave a row behind for the merge step, because a missing
row would silently shrink the cohort instead of recording a blank metric and a
reason.
"""
from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_biophysical_metrics as bio  # noqa: E402
import run_rosetta_metrics as rose  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

gemmi = pytest.importorskip("gemmi")

_HAS_FREESASA = importlib.util.find_spec("freesasa") is not None

ROSETTA_COLUMNS = ["seq_name", "sc", "rosetta_ddg", "rosetta_failures"]
BIO_COLUMNS = ["seq_name", "bsa", "interface_hbonds", "biophysical_failures"]


def _one_row(path: Path) -> tuple[list[str], dict]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) == 1, "each runner emits exactly one row per survivor"
    return list(reader.fieldnames or []), rows[0]


# ── Rosetta binary stubs ─────────────────────────────────────────────

# A clean FastRelax: read -out:path:all and -s out of argv and write the
# deterministic output name main() expects.
_RELAX_OK = """
import sys, pathlib
a = sys.argv
out = pathlib.Path(a[a.index("-out:path:all") + 1])
src = pathlib.Path(a[a.index("-s") + 1])
(out / (src.stem + "_0001.pdb")).write_text("ATOM      1  CA  ALA A   1\\nEND\\n")
"""

# Exits clean having written nothing at all.
_RELAX_SILENT = "import sys\n"

# Exits clean having written a name main() does not predict, which is what the
# -out:suffix fallback is for.
_RELAX_ODD_NAME = """
import sys, pathlib
a = sys.argv
out = pathlib.Path(a[a.index("-out:path:all") + 1])
(out / "unexpected_name.pdb").write_text("ATOM      1  CA  ALA A   1\\nEND\\n")
"""

# Rosetta logs its [ ERROR ] block to stdout, not stderr, hence both streams.
_RELAX_FAILS = """
import sys
print("[ ERROR ] Mover FastRelax failed")
print("to stderr", file=sys.stderr)
sys.exit(3)
"""


def _ia_body(*, sc: str = "0.66512", ddg: str = "-59.7394",
             exit_code: int = 0, write: bool = True, noisy: bool = False,
             columns: str = "total_score dG_separated sc_value description") -> str:
    """An InterfaceAnalyzer stub writing the score file main() will parse."""
    values = {"total_score": "-412.5", "dG_separated": ddg,
              "sc_value": sc, "description": "relaxed_0001"}
    row = " ".join(values.get(c, "0.0") for c in columns.split())
    return f"""
import sys, pathlib
a = sys.argv
if {noisy!r}:
    print("[ ERROR ] IA could not pack the separated state")
    print("ia to stderr", file=sys.stderr)
if {write!r}:
    p = pathlib.Path(a[a.index("-out:file:score_only") + 1])
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("SEQUENCE: AAAA\\nSCORE: {columns}\\nSCORE: {row}\\n")
sys.exit({exit_code})
"""


def _install_rosetta_stubs(tmp_path, monkeypatch, *, relax: str, ia: str) -> Path:
    """Write both binaries into a fake bin dir and make it the whole PATH.

    The shebang names this interpreter by absolute path, so nothing else has
    to be reachable. The stub filenames use the mpi flavour because that is
    the first one `_find_rosetta_binary` tries, which also pins the flavour
    priority order.
    """
    fake_bin = tmp_path / "fake_bin"
    fake_bin.mkdir(exist_ok=True)
    for name, body in (("rosetta_scripts.mpi.linuxgccrelease", relax),
                       ("InterfaceAnalyzer.mpi.linuxgccrelease", ia)):
        p = fake_bin / name
        p.write_text(f"#!{sys.executable}\n{body}")
        p.chmod(0o755)
    monkeypatch.setenv("PATH", str(fake_bin))
    return fake_bin


def _run_rose(monkeypatch, tmp_path, *, xml: Path | None = None,
              seq_name: str = "design_1_seq_1") -> tuple[int, Path]:
    pdb = tmp_path / "canonical.pdb"
    if not pdb.exists():
        complex_pdb(pdb)
    out_csv = tmp_path / "out" / f"{seq_name}_rosetta.csv"
    argv = ["run_rosetta_metrics.py", "--seq-name", seq_name,
            "--canonical-pdb", str(pdb), "--output-csv", str(out_csv)]
    if xml is not None:
        argv += ["--fast-relax-xml", str(xml)]
    monkeypatch.setattr(sys, "argv", argv)
    return rose.main(), out_csv


@pytest.fixture
def relax_xml(tmp_path) -> Path:
    p = tmp_path / "fastrelax.xml"
    p.write_text("<ROSETTASCRIPTS/>\n")
    return p


# ── Rosetta: the happy row ───────────────────────────────────────────

@pytest.mark.logic
def test_a_complete_rosetta_run_writes_one_fully_known_row(
        monkeypatch, tmp_path, relax_xml):
    """Both values come from the stub's score file, so the row is known
    before the stage runs. sc is written to 4 dp and ddg to 3, which is the
    precision the committed per-survivor CSVs carry."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_OK, ia=_ia_body())
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0

    header, row = _one_row(out_csv)
    assert header == ROSETTA_COLUMNS
    assert row == {"seq_name": "design_1_seq_1", "sc": "0.6651",
                   "rosetta_ddg": "-59.739", "rosetta_failures": ""}


@pytest.mark.logic
def test_the_relaxed_pdb_is_what_interface_analyzer_is_given(
        monkeypatch, tmp_path, relax_xml):
    """Relax is the pre-treatment that makes the separated dG meaningful, so
    scoring the unrelaxed input would inflate it by 10-20 REU. The stub
    records which file it was handed."""
    seen = tmp_path / "ia_input.txt"
    ia = f"""
import sys, pathlib
a = sys.argv
pathlib.Path({str(seen)!r}).write_text(a[a.index("-s") + 1])
p = pathlib.Path(a[a.index("-out:file:score_only") + 1])
p.write_text("SCORE: total_score sc_value dG_separated\\nSCORE: -1.0 0.5 -2.0\\n")
"""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_OK, ia=ia)
    rc, _out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    assert seen.read_text().endswith("canonical_0001.pdb")


@pytest.mark.logic
def test_the_xml_defaults_to_the_one_pinned_beside_the_script(monkeypatch, tmp_path):
    """Nextflow does not pass --fast-relax-xml, so the default has to resolve
    to the committed bin/fastrelax_for_ia.xml. A broken default would turn
    every survivor into a fastrelax_xml_not_found row."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_OK, ia=_ia_body())
    rc, out_csv = _run_rose(monkeypatch, tmp_path)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["rosetta_failures"] == ""
    assert (BIN / "fastrelax_for_ia.xml").is_file()


# ── Rosetta: the failure rows ────────────────────────────────────────

@pytest.mark.logic
def test_a_missing_xml_still_emits_a_row_and_exits_non_zero(monkeypatch, tmp_path):
    """Exit 1 tells Nextflow the stage failed; the row is written anyway so
    the merge step keeps the survivor instead of dropping it."""
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=tmp_path / "absent.xml")
    assert rc == 1
    header, row = _one_row(out_csv)
    assert header == ROSETTA_COLUMNS
    assert row == {"seq_name": "design_1_seq_1", "sc": "",
                   "rosetta_ddg": "", "rosetta_failures": "fastrelax_xml_not_found"}


@pytest.mark.logic
def test_neither_binary_being_installed_names_both(monkeypatch, tmp_path, relax_xml):
    """The two stems are looked up separately, so the row says which halves
    of the Rosetta image are missing rather than just 'rosetta failed'."""
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0, "a missing binary is a recorded row, not a stage crash"
    _header, row = _one_row(out_csv)
    assert row["sc"] == "" and row["rosetta_ddg"] == ""
    assert "rosetta_scripts_binary_not_found" in row["rosetta_failures"]
    assert "interface_analyzer_binary_not_found" in row["rosetta_failures"]


@pytest.mark.logic
def test_a_relax_that_exits_non_zero_carries_its_exit_code(
        monkeypatch, tmp_path, relax_xml, capsys):
    """The code is in the token so a segfault and a bad XML are told apart
    from the CSV. Rosetta logs to stdout, so both streams are echoed to
    stderr where .command.err will pick them up."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_FAILS, ia=_ia_body())
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["rosetta_failures"] == "relax_nonzero_exit:3"
    assert row["sc"] == "" and row["rosetta_ddg"] == ""
    err = capsys.readouterr().err
    assert "[relax-stdout] [ ERROR ] Mover FastRelax failed" in err
    assert "[relax-stderr] to stderr" in err


@pytest.mark.logic
def test_a_relax_that_exits_clean_with_no_pdb_is_its_own_failure(
        monkeypatch, tmp_path, relax_xml):
    """Distinct from a non-zero exit. Rosetta returning 0 having written
    nothing is the shape a silently skipped mover takes."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_SILENT, ia=_ia_body())
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["rosetta_failures"] == "relax_no_output_pdb"


@pytest.mark.logic
def test_an_unexpected_output_name_is_picked_up_rather_than_failed(
        monkeypatch, tmp_path, relax_xml):
    """rosetta_scripts' suffix handling varies between builds. The fallback
    takes whatever PDB appeared in the relax directory, so a renamed output
    still reaches InterfaceAnalyzer."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_ODD_NAME, ia=_ia_body())
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["rosetta_failures"] == ""
    assert row["sc"] == "0.6651"


@pytest.mark.logic
def test_interface_analyzer_writing_no_score_file_is_recorded(
        monkeypatch, tmp_path, relax_xml):
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_OK,
                           ia=_ia_body(write=False))
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["rosetta_failures"] == "ia_no_scorefile_produced"
    assert row["sc"] == "" and row["rosetta_ddg"] == ""


@pytest.mark.logic
def test_a_non_zero_interface_analyzer_that_still_wrote_a_score_file_is_both(
        monkeypatch, tmp_path, relax_xml, capsys):
    """The score file is parsed even though the binary complained, so the
    numbers survive and the exit code is carried alongside them. Dropping
    either half would hide one of the two facts."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_OK,
                           ia=_ia_body(exit_code=1, noisy=True))
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["rosetta_failures"] == "ia_nonzero_exit:1"
    assert row["sc"] == "0.6651" and row["rosetta_ddg"] == "-59.739"
    err = capsys.readouterr().err
    assert "[ia-stdout] [ ERROR ] IA could not pack the separated state" in err
    assert "[ia-stderr] ia to stderr" in err


@pytest.mark.logic
@pytest.mark.parametrize("ia_kwargs,expect_sc,expect_ddg,token", [
    ({"columns": "total_score dG_separated description"},
     "", "-59.739", "sc_column_missing"),
    ({"columns": "total_score sc_value description"},
     "0.6651", "", "ddg_column_missing"),
    ({"sc": "nan_or_junk"}, "", "-59.739", "sc_not_numeric"),
    ({"ddg": "not_a_number"}, "0.6651", "", "ddg_not_numeric"),
])
def test_one_unusable_column_leaves_the_other_intact_and_names_itself(
        monkeypatch, tmp_path, relax_xml, ia_kwargs, expect_sc, expect_ddg, token):
    """A changed Rosetta option set drops or mangles one column rather than
    failing outright. The row then carries whichever metric it did get plus a
    token naming the one it did not, which is strictly more than a bare
    'parse failed'. Coercing the bad one to 0.0 would read as a real
    measurement of a terrible interface."""
    _install_rosetta_stubs(tmp_path, monkeypatch, relax=_RELAX_OK,
                           ia=_ia_body(**ia_kwargs))
    rc, out_csv = _run_rose(monkeypatch, tmp_path, xml=relax_xml)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["sc"] == expect_sc
    assert row["rosetta_ddg"] == expect_ddg
    assert row["rosetta_failures"] == token


@pytest.mark.logic
@pytest.mark.parametrize("stage,token", [
    ("relax", "relax_timed_out"),
    ("ia", "ia_timed_out"),
])
def test_a_tool_that_hangs_becomes_a_timeout_token(monkeypatch, tmp_path, stage, token):
    """Relax gets 900 s and InterfaceAnalyzer 300 s. A hang has to end as a
    row rather than as a stage that never returns, so the cap is raised here
    by making subprocess.run raise rather than by waiting it out."""
    def boom(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd, kwargs.get("timeout", 0))
    monkeypatch.setattr(rose.subprocess, "run", boom)

    if stage == "relax":
        got, failures = rose._run_fast_relax(
            "fake_bin", tmp_path / "x.xml", tmp_path / "in.pdb", tmp_path)
        assert got is None
    else:
        sc, ddg, failures = rose._run_interface_analyzer(
            "fake_bin", tmp_path / "in.pdb", tmp_path)
        assert (sc, ddg) == (None, None)
    assert failures == [token]


# ── Biophysical: hydrogen bonds ──────────────────────────────────────

def _np_pdb(path: Path, receptor: list, effector: list,
            others: dict | None = None) -> Path:
    """A PDB whose atom names and elements are stated explicitly.

    `fixtures_runtree.complex_pdb` writes carbon alphas only, which the
    H-bond counter is right to ignore, so this builder is needed to place
    nitrogens and oxygens at distances chosen in advance. `others` adds
    chains beyond the receptor and effector, for the case where the file
    carries a chain neither metric should see.
    """
    line = ("ATOM  {serial:5d}  {name:<3s} ALA {chain}{resnum:4d}    "
            "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
    lines, serial = [], 1
    chains = [("A", receptor), ("B", effector)]
    chains += sorted((others or {}).items())
    for chain, atoms in chains:
        for resnum, (name, (x, y, z)) in enumerate(atoms, start=1):
            lines.append(line.format(serial=serial, name=name, chain=chain,
                                     resnum=resnum, x=x, y=y, z=z, el=name[0]))
            serial += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


@pytest.mark.logic
def test_the_hbond_count_is_the_number_of_cross_chain_n_o_pairs_in_range(tmp_path):
    """Two facing pairs 3.0 A apart and two diagonal pairs 10 A away, so the
    answer is 2 by construction rather than by measurement."""
    p = _np_pdb(tmp_path / "c.pdb",
                [("N", (0, 0, 0)), ("N", (10, 0, 0))],
                [("O", (0, 0, 3.0)), ("O", (10, 0, 3.0))])
    n, err = bio._compute_hbonds(p, "A", "B")
    assert err is None
    assert n == 2


@pytest.mark.logic
def test_carbon_pairs_do_not_count_however_close_they_are(tmp_path):
    """The criterion is N/O to N/O. A carbon 1 A away is a clash, not an
    H-bond, and counting it would inflate every interface."""
    p = _np_pdb(tmp_path / "c.pdb",
                [("N", (0, 0, 0))],
                [("O", (0, 0, 3.0)), ("CA", (0, 0, 1.0))])
    n, err = bio._compute_hbonds(p, "A", "B")
    assert err is None
    assert n == 1


@pytest.mark.logic
@pytest.mark.parametrize("z,expected", [(3.4, 1), (3.5, 1), (3.6, 0)])
def test_the_distance_criterion_includes_its_own_boundary(tmp_path, z, expected):
    """The comparison is `d <= 3.5`, not strict, so a pair sitting exactly on
    the documented 3.5 A cutoff is counted. That makes the published
    criterion attainable rather than one value short of it."""
    p = _np_pdb(tmp_path / "c.pdb", [("N", (0, 0, 0))], [("O", (0, 0, z))])
    n, err = bio._compute_hbonds(p, "A", "B")
    assert err is None
    assert n == expected


@pytest.mark.logic
def test_an_absent_chain_is_named_rather_than_counted_as_zero_hbonds(tmp_path):
    """Zero would read as a real interface with no polar contacts, so a typo
    in the chain id would look like a measurement."""
    p = _np_pdb(tmp_path / "c.pdb", [("N", (0, 0, 0))], [("O", (0, 0, 3.0))])
    n, err = bio._compute_hbonds(p, "A", "Z")
    assert n is None
    assert err == "chain_not_found"


@pytest.mark.logic
def test_a_structure_that_will_not_load_is_reported_not_raised(tmp_path):
    n, err = bio._compute_hbonds(tmp_path / "absent.pdb", "A", "B")
    assert n is None
    assert err is not None and err.startswith("structure_load_failed:")


@pytest.mark.logic
def test_a_third_chain_contributes_to_neither_metric(tmp_path):
    """A chain that is neither the receptor nor the effector is skipped, not
    folded into whichever list was last assigned. Chain C sits 1 A from both
    polar atoms here, so including it would double the count and move the
    interface set."""
    p = _np_pdb(tmp_path / "c.pdb",
                [("N", (0, 0, 0))],
                [("O", (0, 0, 3.0))],
                others={"C": [("N", (0, 0, 1.0)), ("O", (0, 0, 2.0))]})
    n, err = bio._compute_hbonds(p, "A", "B")
    assert err is None
    assert n == 1

    residues, err = bio._interface_residues(p, "A", "B", 5.0)
    assert err is None
    assert {c for c, _n in residues} == {"A", "B"}


@pytest.mark.logic
def test_hydrogens_do_not_make_an_interface(tmp_path):
    """The interface is defined on heavy atoms, matching the rest of the
    pipeline. The only chain-B atom within the cutoff here is a hydrogen, so
    counting it would invent an interface where the heavy atoms have none."""
    p = _np_pdb(tmp_path / "c.pdb",
                [("N", (0, 0, 0))],
                [("H", (0, 0, 1.0)), ("O", (0, 0, 40.0))])
    residues, err = bio._interface_residues(p, "A", "B", 5.0)
    assert err is None
    assert residues == set()


@pytest.mark.logic
def test_an_interface_on_a_structure_that_is_not_there_is_reported(tmp_path):
    """A prediction that never landed must come back as a token. gemmi is
    lenient about malformed content but not about an absent file, which is
    the shape a failed Boltz run actually takes."""
    residues, err = bio._interface_residues(tmp_path / "absent.pdb", "A", "B", 5.0)
    assert residues is None
    assert err is not None and err.startswith("structure_load_failed:")


# ── Biophysical: buried surface area ─────────────────────────────────

@pytest.mark.logic
@pytest.mark.skipif(_HAS_FREESASA, reason="freesasa is installed, so the "
                                          "import-failure branch is unreachable")
def test_freesasa_not_being_installed_is_a_token_not_a_crash(tmp_path):
    """BSA is one of the two gating orthogonal metrics. Without FreeSASA the
    row must come back with a blank bsa and a named reason, because the merge
    step turns a blank into bsa_missing and gates the design out."""
    p = _np_pdb(tmp_path / "c.pdb", [("N", (0, 0, 0))], [("O", (0, 0, 3.0))])
    bsa, err = bio._compute_bsa(p, "A", "B")
    assert bsa is None
    assert err is not None and err.startswith("freesasa_import_failed:")


@pytest.mark.logic
@pytest.mark.skipif(not _HAS_FREESASA, reason="needs freesasa")
def test_a_separated_complex_buries_nothing(tmp_path):
    """BSA is (SASA_A + SASA_B) - SASA_complex, halved. Two chains far apart
    bury no surface, so the only non-trivial part of the arithmetic cancels
    and the answer is 0."""
    p = complex_pdb(tmp_path / "far.pdb", separation=60.0)
    bsa, err = bio._compute_bsa(p, "A", "B")
    assert err is None
    assert bsa == pytest.approx(0.0, abs=1.0)


def _install_fake_freesasa(monkeypatch, *, structure_raises=None):
    """A FreeSASA whose total area is an arithmetic function of the atoms.

    Real FreeSASA is not installed here, and even where it is its numbers are
    not known in advance, so the three-calculation arithmetic in
    `_compute_bsa` cannot be checked against it.

    This stub makes the area

        area(file) = 10 * n_A + 100 * n_B - 5 * n_A * n_B

    over the chain-A and chain-B atom counts in whatever file it is handed.
    The interaction term is what survives the subtraction: an area that were
    merely additive in the two chains would cancel to exactly 0 for every
    input, so only the cross term tests the formula at all.

        BSA = (area(A alone) + area(B alone) - area(complex)) / 2
            = 5 * n_A * n_B / 2

    That also pins the halving, which turns FreeSASA's both-sides number into
    the per-interface convention the filter threshold is set against.
    """
    mod = types.ModuleType("freesasa")

    class Structure:
        def __init__(self, path):
            if structure_raises is not None:
                raise structure_raises
            atoms = [ln for ln in Path(path).read_text().splitlines()
                     if ln.startswith(("ATOM", "HETATM"))]
            self.n_a = sum(1 for ln in atoms if ln[21] == "A")
            self.n_b = sum(1 for ln in atoms if ln[21] == "B")

    class _Result:
        def __init__(self, st):
            self._st = st

        def totalArea(self):  # noqa: N802 - FreeSASA's own spelling
            return (10.0 * self._st.n_a + 100.0 * self._st.n_b
                    - 5.0 * self._st.n_a * self._st.n_b)

    mod.Structure = Structure
    mod.calc = _Result
    monkeypatch.setitem(sys.modules, "freesasa", mod)


@pytest.mark.logic
def test_bsa_is_half_the_surface_the_two_chains_hide_from_each_other(
        monkeypatch, tmp_path):
    """Two chain-A atoms and three chain-B atoms, so the stub's cross term
    gives 5 * 2 * 3 / 2 = 15.00 angstrom squared. The chain isolation is part
    of the arithmetic: a filter that kept the wrong chain, or dropped the
    halving, lands on a different number."""
    _install_fake_freesasa(monkeypatch)
    p = _np_pdb(tmp_path / "c.pdb",
                [("N", (0, 0, 0)), ("O", (1, 0, 0))],
                [("O", (0, 0, 3.0)), ("N", (1, 0, 3.0)), ("O", (2, 0, 3.0))])
    bsa, err = bio._compute_bsa(p, "A", "B")
    assert err is None
    assert bsa == pytest.approx(15.0)


@pytest.mark.logic
def test_the_bsa_value_is_rounded_to_two_decimals(monkeypatch, tmp_path):
    """One chain-A atom against one chain-B atom gives 5 * 1 * 1 / 2 = 2.5,
    and the stage writes 2 dp. The rounding is in the committed column."""
    _install_fake_freesasa(monkeypatch)
    p = _np_pdb(tmp_path / "c.pdb", [("N", (0, 0, 0))], [("O", (0, 0, 3.0))])
    bsa, err = bio._compute_bsa(p, "A", "B")
    assert err is None
    assert bsa == pytest.approx(2.5)
    rc, out_csv = _run_bio(monkeypatch, tmp_path, p)
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["bsa"] == "2.50"
    assert row["biophysical_failures"] == ""


@pytest.mark.logic
def test_a_freesasa_error_is_a_tagged_row_rather_than_an_exception(
        monkeypatch, tmp_path, capsys):
    """FreeSASA rejects some predicted structures outright. The token names
    the exception type so the cause is diagnosable from the CSV, and the
    traceback still reaches stderr."""
    _install_fake_freesasa(monkeypatch,
                           structure_raises=ValueError("unknown radius"))
    p = _np_pdb(tmp_path / "c.pdb", [("N", (0, 0, 0))], [("O", (0, 0, 3.0))])
    bsa, err = bio._compute_bsa(p, "A", "B")
    assert bsa is None
    assert err is not None and err.startswith("freesasa_failed:ValueError:")
    assert "freesasa-traceback" in capsys.readouterr().err


# ── Biophysical main() ───────────────────────────────────────────────

def _run_bio(monkeypatch, tmp_path, pdb: Path,
             *, seq_name: str = "design_1_seq_1",
             receptor: str = "A", effector: str = "B") -> tuple[int, Path]:
    out_csv = tmp_path / "out" / f"{seq_name}_biophysical.csv"
    monkeypatch.setattr(sys, "argv", [
        "run_biophysical_metrics.py",
        "--seq-name", seq_name,
        "--canonical-pdb", str(pdb),
        "--ground-truth", str(pdb),
        "--receptor-chain", receptor, "--effector-chain", effector,
        "--output-csv", str(out_csv),
    ])
    return bio.main(), out_csv


@pytest.mark.logic
def test_the_biophysical_row_carries_the_hbond_count_and_its_own_columns(
        monkeypatch, tmp_path):
    """The hbond column is known by construction. bsa and the failures
    column move together: a blank bsa is always accompanied by a reason, and
    a populated one never is, whether or not FreeSASA is installed here."""
    p = _np_pdb(tmp_path / "c.pdb",
                [("N", (0, 0, 0)), ("N", (10, 0, 0))],
                [("O", (0, 0, 3.0)), ("O", (10, 0, 3.0))])
    rc, out_csv = _run_bio(monkeypatch, tmp_path, p)
    assert rc == 0

    header, row = _one_row(out_csv)
    assert header == BIO_COLUMNS
    assert row["seq_name"] == "design_1_seq_1"
    assert row["interface_hbonds"] == "2"
    assert (row["bsa"] == "") is bool(row["biophysical_failures"])


@pytest.mark.logic
def test_a_canonical_pdb_that_is_not_there_still_writes_a_row(monkeypatch, tmp_path):
    """Both metrics fail and the stage still exits 0 with a row. Exiting
    non-zero would take down the whole fan-out over one absent prediction,
    and writing no row would quietly shrink the cohort."""
    rc, out_csv = _run_bio(monkeypatch, tmp_path, tmp_path / "absent.pdb")
    assert rc == 0
    header, row = _one_row(out_csv)
    assert header == BIO_COLUMNS
    assert row["seq_name"] == "design_1_seq_1"
    assert row["bsa"] == "" and row["interface_hbonds"] == ""
    assert "structure_load_failed" in row["biophysical_failures"]


@pytest.mark.logic
def test_a_wrong_chain_id_is_reported_in_the_row_rather_than_scored(
        monkeypatch, tmp_path):
    p = _np_pdb(tmp_path / "c.pdb", [("N", (0, 0, 0))], [("O", (0, 0, 3.0))])
    rc, out_csv = _run_bio(monkeypatch, tmp_path, p, effector="Z")
    assert rc == 0
    _header, row = _one_row(out_csv)
    assert row["interface_hbonds"] == ""
    assert "chain_not_found" in row["biophysical_failures"]
