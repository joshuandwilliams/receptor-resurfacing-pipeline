"""Cover the chain-ID normalisation that precedes the pose solver.

The stage exists for one reason. `build_contigs.py` falls back to the first
two alphabetical chains when the chain IDs the user named are absent, but the
contig string it emits still names the user's chains. RFDiffusion would then
receive a PDB whose chain IDs do not match its own contig string, and the run
would either crash much later or quietly diffuse against the wrong chain.

So this module forces the receptor PDB's chain ID to `--receptor-chain` and
the effector PDB's to `--effector-chain`, and the property every test here
asserts is that the chain the contig string will name is the chain actually
present in the output.

The failure mode the column arithmetic protects against is subtle: chain ID
lives in PDB column 22 and the fields either side of it are residue number and
insertion code. Writing one character too early or too late corrupts the
residue numbering rather than the chain, which gemmi will happily read back as
a differently numbered structure.
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

import pose_solver_prepare as psp  # noqa: E402
from fixtures_runtree import ca_line, complex_pdb  # noqa: E402


def _single_chain(path: Path, chain: str, n: int = 6, start: int = 1) -> Path:
    """One chain of `n` CA records under `chain`, which is what each of the
    two inputs to this stage is: a separate receptor and effector PDB."""
    return complex_pdb(path, n_receptor=n, n_effector=0,
                       receptor_chain=chain, receptor_start=start)


def _chains(path: Path) -> set[str]:
    """Every chain ID appearing in column 22 of an ATOM or HETATM record."""
    return {line[21] for line in path.read_text().splitlines()
            if line[:6] in ("ATOM  ", "HETATM") and len(line) >= 22}


def _resnums(path: Path) -> list[str]:
    return [line[22:26] for line in path.read_text().splitlines()
            if line[:6] == "ATOM  " and len(line) >= 26]


def _run(monkeypatch, tmp_path: Path, receptor: Path, effector: Path, **extra):
    """Run main() with its two outputs inside tmp_path. Returns
    (exit_code, receptor_out, effector_out)."""
    rec_out = tmp_path / "receptor_pose.pdb"
    eff_out = tmp_path / "effector_pose.pdb"
    argv = ["pose_solver_prepare.py",
            "--receptor", str(receptor),
            "--effector", str(effector),
            "--receptor-out", str(rec_out),
            "--effector-out", str(eff_out)]
    for key, value in extra.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    monkeypatch.setattr(sys, "argv", argv)
    try:
        code = psp.main()
    except SystemExit as exc:
        code = exc.code
    return code, rec_out, eff_out


# ── relabel_chain: the column arithmetic ─────────────────────────────

@pytest.mark.logic
def test_the_chain_id_is_replaced_and_nothing_either_side_of_it_moves(tmp_path):
    """Column 22 is the chain, columns 23-26 are the residue number. The
    output must differ from the input in exactly the chain column."""
    src = _single_chain(tmp_path / "in.pdb", "X", n=4, start=17)
    out = tmp_path / "out.pdb"
    psp.relabel_chain(src, out, "Q")
    assert _chains(out) == {"Q"}
    assert _resnums(out) == _resnums(src)
    for before, after in zip(src.read_text().splitlines(),
                             out.read_text().splitlines()):
        assert len(before) == len(after)
        assert before[:21] == after[:21]
        assert before[22:] == after[22:]


@pytest.mark.logic
def test_the_atom_count_returned_is_atom_and_hetatm_only(tmp_path):
    """TER and ANISOU get the new chain ID too, because a stale chain in a
    TER record splits the chain when a reader honours it. Neither is an
    atom, so neither may be counted: the count is what main() gates the
    "no atom records" error on."""
    src = tmp_path / "mixed.pdb"
    src.write_text("\n".join([
        "REMARK   1 provenance line",
        ca_line(1, "X", 1, 0.0, 0.0, 0.0),
        "ANISOU    1  CA  ALA X   1     1000   1000   1000      0      0      0",
        ca_line(2, "X", 2, 3.8, 0.0, 0.0),
        "HETATM    3  O   HOH X   3       1.000   0.000   0.000  1.00  0.00      O",
        "TER       4      ALA X   3",
        "END",
    ]) + "\n")
    out = tmp_path / "out.pdb"
    # Two ATOM plus one HETATM; the ANISOU and the TER are relabelled but
    # not counted.
    assert psp.relabel_chain(src, out, "B") == 3
    lines = out.read_text().splitlines()
    assert lines[2][21] == "B", "ANISOU must carry the new chain"
    assert lines[5][21] == "B", "TER must carry the new chain"


@pytest.mark.logic
def test_a_non_atom_record_passes_through_untouched(tmp_path):
    """REMARK, HEADER and CRYST1 all have other things in column 22. The
    receptor PDB a user supplies carries them, and rewriting that column
    would corrupt them."""
    src = tmp_path / "in.pdb"
    src.write_text("REMARK   1 QQQQQQ\nCRYST1    1.000    1.000    1.000\n"
                   + ca_line(1, "X", 1, 0.0, 0.0, 0.0) + "\nEND\n")
    out = tmp_path / "out.pdb"
    assert psp.relabel_chain(src, out, "A") == 1
    lines = out.read_text().splitlines()
    assert lines[0] == "REMARK   1 QQQQQQ"
    assert lines[1] == "CRYST1    1.000    1.000    1.000"
    assert lines[3] == "END"


@pytest.mark.logic
def test_an_empty_file_writes_an_empty_file_and_counts_nothing(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text("")
    out = tmp_path / "out.pdb"
    assert psp.relabel_chain(src, out, "A") == 0
    assert out.read_text() == ""


@pytest.mark.logic
def test_an_atom_line_too_short_to_hold_a_chain_passes_through_uncounted(tmp_path):
    """Today's behaviour, pinned. A record shorter than 22 characters cannot
    be indexed at column 22, so it is copied verbatim and does not increment
    the atom count. It therefore keeps whatever chain it had, or none."""
    src = tmp_path / "in.pdb"
    src.write_text("ATOM  truncated\n" + ca_line(1, "X", 1, 0.0, 0.0, 0.0) + "\n")
    out = tmp_path / "out.pdb"
    assert psp.relabel_chain(src, out, "A") == 1
    assert out.read_text().splitlines()[0] == "ATOM  truncated"


# ── main: the normalisation contract ─────────────────────────────────

@pytest.mark.logic
def test_the_two_outputs_carry_exactly_the_requested_chain_ids(monkeypatch, tmp_path):
    """The whole point of the stage. The inputs name chains X and Y, the
    contig string downstream will name A and B, and after this stage A and B
    are what the PDBs hold."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, rec_out, eff_out = _run(monkeypatch, tmp_path, rec, eff)
    assert code == 0
    assert _chains(rec_out) == {"A"}
    assert _chains(eff_out) == {"B"}


@pytest.mark.logic
def test_a_chain_already_correct_is_left_correct(monkeypatch, tmp_path):
    """The stage must be a no-op on an input that already matches, because it
    runs on every pose whether or not the chain IDs needed fixing."""
    rec = _single_chain(tmp_path / "rec.pdb", "A")
    eff = _single_chain(tmp_path / "eff.pdb", "B")
    code, rec_out, eff_out = _run(monkeypatch, tmp_path, rec, eff)
    assert code == 0
    assert rec_out.read_text() == rec.read_text()
    assert eff_out.read_text() == eff.read_text()


@pytest.mark.logic
def test_non_default_chain_ids_reach_the_output(monkeypatch, tmp_path):
    """The campaign's own params name the receptor and effector chains, so
    the flags have to arrive rather than the A/B defaults always winning."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, rec_out, eff_out = _run(monkeypatch, tmp_path, rec, eff,
                                  receptor_chain="E", effector_chain="C")
    assert code == 0
    assert _chains(rec_out) == {"E"}
    assert _chains(eff_out) == {"C"}


@pytest.mark.logic
def test_a_lowercase_chain_id_is_upcased(monkeypatch, tmp_path):
    """PDB chain IDs are case-sensitive, so 'a' and 'A' are different chains.
    Upcasing means a params file written with a lowercase chain produces the
    same PDB as one written with an uppercase one, rather than a chain the
    contig string cannot name."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, rec_out, eff_out = _run(monkeypatch, tmp_path, rec, eff,
                                  receptor_chain="a", effector_chain="b")
    assert code == 0
    assert _chains(rec_out) == {"A"}
    assert _chains(eff_out) == {"B"}


@pytest.mark.logic
def test_residue_numbering_is_preserved_across_the_rewrite(monkeypatch, tmp_path):
    """The contig string names residues by their PDB number, so renumbering
    here would make every fixed segment point somewhere else."""
    rec = _single_chain(tmp_path / "rec.pdb", "X", n=5, start=403)
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    _code, rec_out, _eff_out = _run(monkeypatch, tmp_path, rec, eff)
    assert [int(n) for n in _resnums(rec_out)] == [403, 404, 405, 406, 407]


@pytest.mark.logic
def test_the_atom_counts_are_reported_per_file(monkeypatch, tmp_path, capsys):
    """The two counts on stdout are the only record in the Nextflow log that
    the stage saw a real structure rather than an empty file."""
    rec = _single_chain(tmp_path / "rec.pdb", "X", n=7)
    eff = _single_chain(tmp_path / "eff.pdb", "Y", n=3)
    _run(monkeypatch, tmp_path, rec, eff)
    printed = capsys.readouterr().out
    assert "(7 atoms, chain A)" in printed
    assert "(3 atoms, chain B)" in printed


# ── main: the guards ─────────────────────────────────────────────────

@pytest.mark.logic
def test_one_chain_id_for_both_inputs_is_refused(monkeypatch, tmp_path):
    """Both chains relabelled to A would make the two-chain complex a single
    chain, and RFDiffusion would diffuse the effector as part of the
    receptor. Refusing is the only safe answer because the stage cannot know
    which of the two the user meant."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, _r, _e = _run(monkeypatch, tmp_path, rec, eff,
                        receptor_chain="A", effector_chain="A")
    assert isinstance(code, str)
    assert "both 'A'" in code


@pytest.mark.logic
def test_a_case_difference_is_not_a_difference(monkeypatch, tmp_path):
    """'A' and 'a' both upcase to A, so they collide. Comparing before
    upcasing would let the collision through."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, _r, _e = _run(monkeypatch, tmp_path, rec, eff,
                        receptor_chain="A", effector_chain="a")
    assert isinstance(code, str) and "both 'A'" in code


@pytest.mark.logic
def test_a_multi_character_chain_id_is_refused(monkeypatch, tmp_path):
    """Column 22 holds one character. A two-character ID would write the
    first into the chain column and shift nothing, silently dropping the
    second, so the PDB would carry a chain the contig string never names."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, _r, _e = _run(monkeypatch, tmp_path, rec, eff,
                        receptor_chain="AB", effector_chain="C")
    assert isinstance(code, str)
    assert "single characters" in code
    assert "receptor='AB'" in code


@pytest.mark.logic
def test_an_empty_chain_id_is_refused_as_well(monkeypatch, tmp_path):
    """Zero characters fails the same single-character check as two."""
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, _r, _e = _run(monkeypatch, tmp_path, rec, eff,
                        receptor_chain="", effector_chain="B")
    assert isinstance(code, str) and "single characters" in code


@pytest.mark.logic
def test_a_receptor_with_no_atoms_fails_the_stage(monkeypatch, tmp_path):
    """An empty receptor would hand the pose solver nothing to dock against
    and the contig string would reference residues that do not exist."""
    rec = tmp_path / "rec.pdb"
    rec.write_text("REMARK   1 header only\nEND\n")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, _r, _e = _run(monkeypatch, tmp_path, rec, eff)
    assert isinstance(code, str)
    assert "receptor PDB" in code and "no atom records" in code


@pytest.mark.logic
def test_an_effector_with_no_atoms_fails_the_stage(monkeypatch, tmp_path):
    rec = _single_chain(tmp_path / "rec.pdb", "X")
    eff = tmp_path / "eff.pdb"
    eff.write_text("END\n")
    code, _r, _e = _run(monkeypatch, tmp_path, rec, eff)
    assert isinstance(code, str)
    assert "effector PDB" in code and "no atom records" in code


@pytest.mark.logic
def test_the_empty_input_check_happens_after_both_files_are_written(
        monkeypatch, tmp_path, capsys):
    """Today's behaviour, pinned. Both outputs exist and both success lines
    are already on stdout by the time the error fires, so a Nextflow log read
    top-down shows 'Wrote receptor_pose.pdb (0 atoms, chain A)' as if it had
    succeeded."""
    rec = tmp_path / "rec.pdb"
    rec.write_text("REMARK   1 header only\nEND\n")
    eff = _single_chain(tmp_path / "eff.pdb", "Y")
    code, rec_out, eff_out = _run(monkeypatch, tmp_path, rec, eff)
    assert isinstance(code, str)
    assert rec_out.is_file() and eff_out.is_file()
    assert "(0 atoms, chain A)" in capsys.readouterr().out
