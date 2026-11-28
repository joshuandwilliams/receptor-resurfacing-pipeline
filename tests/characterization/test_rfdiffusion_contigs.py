"""Drive the contig preprocessing CLI the RESOLVE_CONTIGS process runs.

`contig_utils.resolve_contigs` does the parsing and is covered by
`test_contig_utils.py`. What this file covers is the CLI around it: the one
file it writes, which is the only thing the Nextflow process passes on.

The failure mode worth catching is the file rather than the string. The
resolved contig string is read back by the next process as
`it.text.trim()`, so a stray newline, a dropped block or the printed string
and the written string diverging would all reach RFDiffusion as a different
contig than the one the log reports.

The resolution itself is asserted on a PDB whose chain extents are stated at
the call site, so "A" resolving to "A1-20" is arithmetic over the fixture
rather than agreement with a previous run.
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

import rfdiffusion_contigs as rc  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402


def _pdb(tmp_path: Path, *, receptor_start: int = 1, n_receptor: int = 20,
         n_effector: int = 10) -> Path:
    """A two-chain complex: chain A residues receptor_start..+n_receptor-1,
    chain B residues 1..n_effector."""
    return complex_pdb(tmp_path / "input.pdb", n_receptor=n_receptor,
                       n_effector=n_effector, receptor_start=receptor_start,
                       effector_start=1)


def _run(monkeypatch, pdb: Path, contigs: str, out: Path | None = None):
    """Run main() and return what landed in the output file."""
    argv = ["rfdiffusion_contigs.py", "--contigs", contigs, "--pdb", str(pdb)]
    if out is not None:
        argv += ["--output", str(out)]
    monkeypatch.setattr(sys, "argv", argv)
    rc.main()
    return out


# ── Bare chain letters ───────────────────────────────────────────────

@pytest.mark.logic
def test_a_bare_chain_letter_becomes_that_chain_s_full_extent(monkeypatch, tmp_path):
    """Chain A holds residues 1..20 and chain B holds 1..10 in the fixture,
    so "A B" must resolve to "A1-20 B1-10". RFDiffusion cannot read a bare
    letter, so a passthrough here would fail the run."""
    pdb = _pdb(tmp_path, n_receptor=20, n_effector=10)
    out = _run(monkeypatch, pdb, "A B", tmp_path / "processed.txt")
    assert out.read_text() == "A1-20 B1-10"


@pytest.mark.logic
def test_the_extent_comes_from_the_pdb_not_from_one(monkeypatch, tmp_path):
    """Chain A numbered from 403 must resolve to A403-422, not A1-20. Reading
    the length and assuming 1-based numbering is the mistake this catches."""
    pdb = _pdb(tmp_path, receptor_start=403, n_receptor=20)
    out = _run(monkeypatch, pdb, "A", tmp_path / "processed.txt")
    assert out.read_text() == "A403-422"


@pytest.mark.logic
def test_a_chain_absent_from_the_pdb_warns_and_passes_through(monkeypatch, tmp_path,
                                                              capsys):
    """Today's behaviour, pinned. "Z" is not in the PDB, so it is written out
    unresolved with a WARNING on stderr. RFDiffusion then rejects it, which is
    a louder failure than silently resolving it to something else."""
    pdb = _pdb(tmp_path)
    out = _run(monkeypatch, pdb, "A Z", tmp_path / "processed.txt")
    assert out.read_text() == "A1-20 Z"
    assert "Chain 'Z' not found" in capsys.readouterr().err


# ── Fixed-segment remapping ──────────────────────────────────────────

@pytest.mark.logic
def test_fixed_segments_are_remapped_onto_the_pdb_s_own_numbering(monkeypatch,
                                                                  tmp_path):
    """The user writes positions 1-5 and 16-20 meaning the first five and the
    last five residues of chain A. The PDB numbers that chain from 101, so
    the emitted contig must name 101-105 and 116-120. Passing the user's
    numbers straight through would make RFDiffusion keep residues that do not
    exist."""
    pdb = _pdb(tmp_path, receptor_start=101, n_receptor=20)
    out = _run(monkeypatch, pdb, "A1-5/10-10/A16-20", tmp_path / "processed.txt")
    assert out.read_text() == "A101-105/10-10/A116-120"


@pytest.mark.logic
def test_a_contig_already_in_the_pdb_s_numbering_is_unchanged(monkeypatch, tmp_path):
    """The remap must be idempotent, because the same contig string is
    resolved on every branch and a second pass must not shift it again."""
    pdb = _pdb(tmp_path, receptor_start=101, n_receptor=20)
    out = _run(monkeypatch, pdb, "A101-105/10-10/A116-120",
               tmp_path / "processed.txt")
    assert out.read_text() == "A101-105/10-10/A116-120"


@pytest.mark.logic
def test_a_de_novo_length_range_survives_untouched(monkeypatch, tmp_path):
    """"10-20" is a length range for a segment RFDiffusion will build, not a
    residue range. Remapping it would turn the whole point of the run into a
    fixed segment."""
    pdb = _pdb(tmp_path, n_receptor=20)
    out = _run(monkeypatch, pdb, "A1-5/10-20/A16-20 B", tmp_path / "processed.txt")
    assert out.read_text() == "A1-5/10-20/A16-20 B1-10"


# ── The output file ──────────────────────────────────────────────────

@pytest.mark.logic
def test_the_file_holds_the_resolved_string_with_no_trailing_newline(monkeypatch,
                                                                     tmp_path):
    """RESOLVE_CONTIGS reads this as `it.text.trim()`, so a trailing newline
    would be tolerated. Asserting its absence pins what the file is, because
    anything else reading the file raw gets the contig and nothing else."""
    pdb = _pdb(tmp_path)
    out = _run(monkeypatch, pdb, "A B", tmp_path / "processed.txt")
    text = out.read_text()
    assert text == text.strip()
    assert "\n" not in text


@pytest.mark.logic
def test_the_printed_string_and_the_written_string_are_the_same(monkeypatch,
                                                               tmp_path, capsys):
    """The Nextflow log shows the printed one and RFDiffusion receives the
    written one. If they could diverge, the log would describe a run that did
    not happen."""
    pdb = _pdb(tmp_path, receptor_start=101)
    out = _run(monkeypatch, pdb, "A1-5/10-10/A16-20 B", tmp_path / "processed.txt")
    printed = capsys.readouterr().out
    assert f"Processed contigs: {out.read_text()}" in printed
    assert "Raw contigs:       A1-5/10-10/A16-20 B" in printed


@pytest.mark.logic
def test_the_default_output_name_is_the_one_the_process_publishes(monkeypatch,
                                                                 tmp_path):
    """RESOLVE_CONTIGS declares `path "processed_contigs.txt"` as its only
    output, so the default has to match or the process fails to publish."""
    pdb = _pdb(tmp_path)
    monkeypatch.chdir(tmp_path)
    _run(monkeypatch, pdb, "A")
    assert (tmp_path / "processed_contigs.txt").read_text() == "A1-20"


@pytest.mark.logic
def test_a_nested_output_path_is_honoured(monkeypatch, tmp_path):
    """--output is passed explicitly by some callers, so it must reach the
    write rather than the default winning."""
    pdb = _pdb(tmp_path)
    (tmp_path / "sub").mkdir()
    out = _run(monkeypatch, pdb, "A", tmp_path / "sub" / "named.txt")
    assert out.read_text() == "A1-20"
    assert not (tmp_path / "processed_contigs.txt").exists()


# ── Defects pinned ───────────────────────────────────────────────────

@pytest.mark.logic
def test_a_bare_de_novo_number_raises_rather_than_expanding(monkeypatch, tmp_path):
    """Today's behaviour, pinned. The module docstring advertises `"5" ->
    "5-5"`, but a block with no chain-prefixed fixed segment has no chain to
    resolve against and the parse raises ValueError instead."""
    pdb = _pdb(tmp_path)
    with pytest.raises(ValueError) as exc:
        _run(monkeypatch, pdb, "5", tmp_path / "processed.txt")
    assert "no chain-prefixed fixed segment" in str(exc.value)


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: the module docstring advertises the shorthand '5' -> '5-5' for "
    "a bare de novo length, but a block holding only a bare number raises "
    "ValueError because it carries no chain to resolve against. Not fixed "
    "here because bin/ is unchanged under a characterization pass; the fix "
    "is either to expand a chainless numeric block or to drop the claim "
    "from the docstring."))
def test_a_bare_de_novo_number_should_expand_to_a_single_length(monkeypatch,
                                                                tmp_path):
    pdb = _pdb(tmp_path)
    out = _run(monkeypatch, pdb, "5", tmp_path / "processed.txt")
    assert out.read_text() == "5-5"


@pytest.mark.logic
def test_a_single_residue_shorthand_outside_the_pdb_is_remapped_to_the_first(
        monkeypatch, tmp_path):
    """Today's behaviour, pinned. The docstring advertises `"A403" ->
    "A403-403"`, and the expansion does happen, but the remap then runs: 403
    is not in a chain numbered 1..20, so it collapses onto the chain's first
    residue and the contig silently names a different residue from the one
    asked for."""
    pdb = _pdb(tmp_path, receptor_start=1, n_receptor=20)
    out = _run(monkeypatch, pdb, "A403", tmp_path / "processed.txt")
    assert out.read_text() == "A1-1"


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: a fixed residue absent from the PDB is remapped onto the "
    "chain's first residue instead of being reported, so a typo in a contig "
    "silently keeps the wrong residue. Not fixed here because bin/ is "
    "unchanged under a characterization pass; the fix belongs in "
    "contig_utils.remap_segments_to_pdb, which should refuse a fixed "
    "residue number it cannot find."))
def test_a_fixed_residue_absent_from_the_pdb_should_be_refused(monkeypatch,
                                                               tmp_path):
    pdb = _pdb(tmp_path, receptor_start=1, n_receptor=20)
    with pytest.raises(ValueError):
        _run(monkeypatch, pdb, "A403", tmp_path / "processed.txt")
