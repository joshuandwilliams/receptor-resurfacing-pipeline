"""Cover the negative-control sequence builders.

These produce the scrambled and polyA controls the campaign compares real
designs against. A control that is accidentally not a control makes the
comparison meaningless, and nothing in the pipeline would report it: the
controls run through exactly the same stages and produce exactly the same
columns.

Two properties carry the whole thing. The scramble must preserve composition,
so length, charge and hydrophobicity are held constant and only the
arrangement moves. And both controls must leave every residue outside the
design region untouched, or they stop being controls for the design region
and become controls for the whole receptor.

The last section drives `main()`, which is the stage. Its three FASTAs are the
only thing the negative-steering orchestrator downstream sees, so the tests
there read the sequences back out of the files.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import build_control_sequences as bcs  # noqa: E402
from fixtures_runtree import sequence_pdb  # noqa: E402

RECEPTOR = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"      # 33 residues
DESIGN_1B = list(range(11, 21))                      # positions 11..20


# ── Design-region index file ─────────────────────────────────────────

@pytest.mark.logic
def test_indices_are_read_from_a_comma_list(tmp_path):
    f = tmp_path / "dr.txt"
    f.write_text("# derived from x.pdb\n11,12,13,14\n")
    assert bcs._load_design_region_indices_1b(f) == [11, 12, 13, 14]


@pytest.mark.logic
def test_comment_lines_are_skipped(tmp_path):
    f = tmp_path / "dr.txt"
    f.write_text("# header\n# coordinate system: 1-based\n5,6\n")
    assert bcs._load_design_region_indices_1b(f) == [5, 6]


@pytest.mark.logic
def test_a_file_of_nothing_but_comments_is_rejected(tmp_path):
    """The derive stages write a block of comment lines and then the indices,
    so a file that is all comments is one whose indices never arrived.
    Returning an empty list for it would build two controls identical to the
    design and report success."""
    f = tmp_path / "dr.txt"
    f.write_text("# derived from design_3.pdb\n# design region size: 0\n")
    with pytest.raises(ValueError) as exc:
        bcs._load_design_region_indices_1b(f)
    assert "no non-comment content" in str(exc.value)


@pytest.mark.logic
def test_a_non_integer_entry_names_the_file(tmp_path):
    """Returning a short list instead would scramble fewer positions than
    intended and the control would quietly become weaker."""
    f = tmp_path / "dr.txt"
    f.write_text("11,notanumber,13\n")
    with pytest.raises(ValueError) as exc:
        bcs._load_design_region_indices_1b(f)
    assert "dr.txt" in str(exc.value)


# ── Effector FASTA ───────────────────────────────────────────────────

@pytest.mark.logic
def test_a_wrapped_fasta_is_joined(tmp_path):
    f = tmp_path / "eff.fasta"
    f.write_text(">effector\nAAAA\nCCCC\n")
    assert bcs._load_effector_sequence(f) == "AAAACCCC"


@pytest.mark.logic
def test_an_empty_fasta_raises_rather_than_returning_nothing(tmp_path):
    """An empty effector would be written into the control's input and
    predicted as a single chain, which is not the comparison intended."""
    f = tmp_path / "eff.fasta"
    f.write_text(">effector\n\n")
    with pytest.raises(ValueError):
        bcs._load_effector_sequence(f)


# ── The scrambled control ────────────────────────────────────────────

@pytest.mark.logic
def test_the_scramble_preserves_overall_composition():
    """Composition is preserved rather than substituted, so length, charge
    and hydrophobicity are constant and only the arrangement moves."""
    out = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=1)
    assert Counter(out) == Counter(RECEPTOR)
    assert len(out) == len(RECEPTOR)


@pytest.mark.logic
def test_the_scramble_touches_nothing_outside_the_design_region():
    """Otherwise it stops being a control for the design region."""
    out = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=1)
    for i, (a, b) in enumerate(zip(RECEPTOR, out), start=1):
        if i not in DESIGN_1B:
            assert a == b, f"position {i} changed outside the design region"


@pytest.mark.logic
def test_the_scrambled_region_is_a_permutation_of_itself():
    out = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=1)
    before = [RECEPTOR[i - 1] for i in DESIGN_1B]
    after = [out[i - 1] for i in DESIGN_1B]
    assert Counter(before) == Counter(after)


@pytest.mark.logic
def test_the_same_seed_gives_the_same_sequence():
    """Re-running one design twice must give identical output, or a rerun
    silently becomes a different control."""
    a = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=42)
    b = bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=42)
    assert a == b


@pytest.mark.logic
def test_different_seeds_generally_give_different_arrangements():
    """Two designs in one cohort must not receive the same shuffle, or the
    controls stop being independent."""
    outs = {bcs._build_scrambled_receptor(RECEPTOR, DESIGN_1B, seed=s)
            for s in range(12)}
    assert len(outs) > 1


@pytest.mark.logic
def test_an_empty_design_region_returns_the_sequence_unchanged():
    """Documented: the caller decides whether that is an error."""
    assert bcs._build_scrambled_receptor(RECEPTOR, [], seed=1) == RECEPTOR


@pytest.mark.logic
def test_indices_outside_the_receptor_are_dropped_not_fatal():
    """A guard against an upstream off-by-one. Indexing past the end would
    raise, and indexing from 0 would shift the whole region by one."""
    out = bcs._build_scrambled_receptor(RECEPTOR, [1, 2, 999, 0, -5], seed=1)
    assert len(out) == len(RECEPTOR)
    assert Counter(out) == Counter(RECEPTOR)
    assert out[2:] == RECEPTOR[2:]


# ── The polyA control ────────────────────────────────────────────────

@pytest.mark.logic
def test_polya_substitutes_the_design_region_with_alanine():
    out = bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)
    for i in DESIGN_1B:
        assert out[i - 1] in ("A", "G")


@pytest.mark.logic
def test_polya_preserves_native_glycine():
    """Alanine where the backbone needs glycine forces phi/psi out of the
    only Ramachandran well that fits, which would confound the geometric
    comparison the control exists to make.
    """
    seq = "MKTGGGAAAKRT"
    design = list(range(1, len(seq) + 1))
    out = bcs._build_polyA_receptor(seq, design)
    for i, residue in enumerate(seq):
        if residue == "G":
            assert out[i] == "G", f"glycine at {i + 1} must survive"
        else:
            assert out[i] == "A"


@pytest.mark.logic
def test_polya_touches_nothing_outside_the_design_region():
    out = bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)
    for i, (a, b) in enumerate(zip(RECEPTOR, out), start=1):
        if i not in DESIGN_1B:
            assert a == b


@pytest.mark.logic
def test_polya_preserves_length():
    assert len(bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)) == len(RECEPTOR)


@pytest.mark.logic
def test_polya_with_an_empty_region_changes_nothing():
    assert bcs._build_polyA_receptor(RECEPTOR, []) == RECEPTOR


@pytest.mark.logic
def test_polya_is_idempotent():
    once = bcs._build_polyA_receptor(RECEPTOR, DESIGN_1B)
    assert bcs._build_polyA_receptor(once, DESIGN_1B) == once


# ── Seed derivation ──────────────────────────────────────────────────

@pytest.mark.logic
def test_the_seed_is_stable_for_one_design():
    a = bcs._default_seed_from_pdb(Path("/x/design_3.pdb"))
    b = bcs._default_seed_from_pdb(Path("/somewhere/else/design_3.pdb"))
    assert a == b, "the seed comes from the stem, not the full path"


@pytest.mark.logic
def test_different_designs_get_different_seeds():
    seeds = {bcs._default_seed_from_pdb(Path(f"/x/design_{i}.pdb"))
             for i in range(20)}
    assert len(seeds) == 20


# ── Output ───────────────────────────────────────────────────────────

@pytest.mark.logic
def test_the_fasta_is_one_header_and_one_body_line(tmp_path):
    out = tmp_path / "nested" / "ctrl.fasta"
    bcs._write_single_record_fasta(out, "control_scrambled", "MKTAY")
    lines = out.read_text().splitlines()
    assert lines == [">control_scrambled", "MKTAY"]


@pytest.mark.logic
def test_the_output_directory_is_created(tmp_path):
    out = tmp_path / "a" / "b" / "c.fasta"
    bcs._write_single_record_fasta(out, "x", "MKT")
    assert out.is_file()


# ── The stage, end to end ────────────────────────────────────────────
#
# One receptor serves the whole section. Positions 4..7 of MKTGYIAKQR are
# G, Y, I and A, so the polyA control is MKT + GAAA + KQR with the native
# glycine at position 4 surviving, and the scramble is a permutation of
# {G, Y, I, A} over those four positions and nothing else.
STAGE_RECEPTOR = "MKTGYIAKQR"          # 10 residues
STAGE_DESIGN = "4,5,6,7"
STAGE_POLYA = "MKTGAAAKQR"
STAGE_EFFECTOR = "WVCD"


def _stage(tmp_path: Path, *, receptor: str = STAGE_RECEPTOR,
           design: str = STAGE_DESIGN, effector: str = STAGE_EFFECTOR,
           stem: str = "design_3", chain: str = "A"):
    """The three inputs the stage takes, in the layout Nextflow stages them.

    `stem` is the design PDB's name, which is also what the default scramble
    seed is hashed from, so two calls with the same stem are two runs of the
    same design.
    """
    pdb = sequence_pdb(tmp_path / stem / f"{stem}.pdb", sequence=receptor,
                       chain=chain)
    region = tmp_path / stem / "design_region.txt"
    region.write_text(
        f"# derived from {stem}.pdb\n"
        "# coordinate system: 1-based positional indices on the receptor\n"
        f"{design}\n"
    )
    eff = tmp_path / stem / "effector.fasta"
    eff.write_text(f">effector\n{effector}\n")
    return pdb, region, eff


def _run(pdb: Path, region: Path, eff: Path, out_dir: Path, **extra) -> int:
    argv = ["build_control_sequences.py",
            "--source-pdb", str(pdb),
            "--design-region-file", str(region),
            "--effector-fasta", str(eff),
            "--output-dir", str(out_dir)]
    for k, v in extra.items():
        argv += [f"--{k.replace('_', '-')}", str(v)]
    old = sys.argv
    try:
        sys.argv = argv
        try:
            return int(bcs.main() or 0)
        except SystemExit as exc:
            return int(exc.code or 0)
    finally:
        sys.argv = old


def _fasta(path: Path) -> tuple[str, str]:
    """The (header, sequence) of a single-record FASTA."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    return lines[0], "".join(lines[1:])


def _controls(out_dir: Path) -> tuple[Path, Path, Path]:
    """The three files the negative-steering orchestrator is handed."""
    return (out_dir / "control_scrambled_receptor.fasta",
            out_dir / "control_polyA_receptor.fasta",
            out_dir / "effector.fasta")


@pytest.mark.logic
def test_the_polya_fasta_holds_the_substituted_sequence(tmp_path):
    """Stated at the call site: MKTGYIAKQR with 4..7 polyA-ed is MKTGAAAKQR.
    The glycine at 4 survives and the three residues outside the region at
    either end are untouched."""
    pdb, region, eff = _stage(tmp_path)
    assert _run(pdb, region, eff, tmp_path / "out") == 0
    _scrambled, polya, _effector = _controls(tmp_path / "out")
    header, seq = _fasta(polya)
    assert seq == STAGE_POLYA
    assert header == ">receptor"


@pytest.mark.logic
def test_the_scrambled_fasta_holds_a_permutation_of_the_design_region(tmp_path):
    """Composition over the whole chain, and the flanks character for
    character. A scramble that reached outside positions 4..7 would still
    pass a composition check on its own."""
    pdb, region, eff = _stage(tmp_path)
    _run(pdb, region, eff, tmp_path / "out")
    scrambled, _polya, _effector = _controls(tmp_path / "out")
    _header, seq = _fasta(scrambled)
    assert Counter(seq) == Counter(STAGE_RECEPTOR)
    assert seq[:3] == "MKT" and seq[7:] == "KQR"
    assert Counter(seq[3:7]) == Counter("GYIA")


@pytest.mark.logic
def test_the_effector_is_passed_through_unchanged(tmp_path):
    """The controls test the receptor binding face. An effector that differed
    from the one the real designs were predicted against would make the
    comparison a comparison of two different complexes."""
    pdb, region, eff = _stage(tmp_path)
    _run(pdb, region, eff, tmp_path / "out")
    _s, _p, effector = _controls(tmp_path / "out")
    assert _fasta(effector) == (">effector", STAGE_EFFECTOR)


@pytest.mark.logic
def test_the_record_name_is_receptor_in_both_controls(tmp_path):
    """The orchestrator pairs --receptor-fasta with --effector-fasta by
    record name, so a control named after its control type would not pair."""
    pdb, region, eff = _stage(tmp_path)
    _run(pdb, region, eff, tmp_path / "out")
    scrambled, polya, _e = _controls(tmp_path / "out")
    assert _fasta(scrambled)[0] == ">receptor"
    assert _fasta(polya)[0] == ">receptor"


@pytest.mark.logic
def test_the_default_seed_is_the_hash_of_the_design_pdb_stem(tmp_path):
    """Running with no --seed must give what running with the stem's hash
    gives, or the documented reproducibility guarantee is not the one in
    force."""
    pdb, region, eff = _stage(tmp_path)
    assert _run(pdb, region, eff, tmp_path / "default") == 0
    assert _run(pdb, region, eff, tmp_path / "explicit",
                seed=bcs._default_seed_from_pdb(pdb)) == 0
    assert (_fasta(_controls(tmp_path / "default")[0])
            == _fasta(_controls(tmp_path / "explicit")[0]))


@pytest.mark.logic
def test_rerunning_one_design_writes_the_same_control(tmp_path):
    """Two runs of the same design, same stem, no seed given. A rerun that
    produced a different scramble would silently replace the control a
    campaign's numbers were measured against."""
    pdb, region, eff = _stage(tmp_path)
    _run(pdb, region, eff, tmp_path / "first")
    _run(pdb, region, eff, tmp_path / "second")
    assert (_controls(tmp_path / "first")[0].read_text()
            == _controls(tmp_path / "second")[0].read_text())


@pytest.mark.logic
def test_the_receptor_chain_flag_reaches_the_extraction(tmp_path):
    """The input complex does not always put the receptor on A. Reading the
    default chain regardless would raise, or worse, read the effector."""
    pdb, region, eff = _stage(tmp_path, chain="C")
    assert _run(pdb, region, eff, tmp_path / "out", receptor_chain="C") == 0
    assert _fasta(_controls(tmp_path / "out")[1])[1] == STAGE_POLYA


@pytest.mark.logic
def test_the_output_directory_is_created_by_the_stage(tmp_path):
    pdb, region, eff = _stage(tmp_path)
    assert _run(pdb, region, eff, tmp_path / "deep" / "nested" / "out") == 0
    for path in _controls(tmp_path / "deep" / "nested" / "out"):
        assert path.is_file()


# ── The degenerate design region ─────────────────────────────────────

@pytest.mark.logic
def test_an_all_glycine_design_region_warns_that_the_control_has_no_signal(
        tmp_path, capsys):
    """Every position in the region is already glycine, so the glycine
    carve-out makes the substitution a no-op and the polyA control is the
    design itself. The run still succeeds, because the cohort's other designs
    are unaffected, but it has to say so."""
    pdb, region, eff = _stage(tmp_path, receptor="MKTGGGGKQR")
    assert _run(pdb, region, eff, tmp_path / "out") == 0
    assert "carries no signal" in capsys.readouterr().err
    assert _fasta(_controls(tmp_path / "out")[1])[1] == "MKTGGGGKQR"


@pytest.mark.logic
def test_an_all_alanine_design_region_also_warns(tmp_path, capsys):
    """Today's behaviour, pinned. The condition is `polyA == source`, which an
    all-alanine region meets just as an all-glycine one does."""
    pdb, region, eff = _stage(tmp_path, receptor="MKTAAAAKQR")
    assert _run(pdb, region, eff, tmp_path / "out") == 0
    assert "carries no signal" in capsys.readouterr().err


@pytest.mark.xfail(strict=False, reason=(
    "The no-signal warning always says 'every design-region position is "
    "already glycine', but it fires whenever polyA leaves the sequence "
    "unchanged, which an all-alanine region does too. Not fixed here: "
    "bin/ is not modified by the characterization tier, and the message is "
    "advisory — the exit code and the written FASTAs are correct either way."))
@pytest.mark.logic
def test_the_no_signal_warning_names_the_residue_it_actually_found(
        tmp_path, capsys):
    pdb, region, eff = _stage(tmp_path, receptor="MKTAAAAKQR")
    _run(pdb, region, eff, tmp_path / "out")
    err = capsys.readouterr().err
    assert "glycine" not in err


# ── The invariants the stage refuses to break ────────────────────────

@pytest.mark.logic
def test_a_scramble_that_changed_the_length_is_a_hard_error(tmp_path, monkeypatch):
    """The guard rather than the builder. A control of the wrong length would
    be predicted as a different complex, and every metric the comparison rests
    on would be measured against it, so this is a crash and not a warning."""
    pdb, region, eff = _stage(tmp_path)
    monkeypatch.setattr(bcs, "_build_scrambled_receptor",
                        lambda seq, idx, seed: seq[:-1])
    with pytest.raises(RuntimeError, match="length"):
        _run(pdb, region, eff, tmp_path / "out")


@pytest.mark.logic
def test_a_scramble_that_changed_the_composition_is_a_hard_error(
        tmp_path, monkeypatch):
    """Same length, different residues. The scrambled control exists to hold
    composition constant, so a substitution dressed as a permutation makes it
    a different control than the one reported."""
    pdb, region, eff = _stage(tmp_path)
    monkeypatch.setattr(bcs, "_build_scrambled_receptor",
                        lambda seq, idx, seed: "A" * len(seq))
    with pytest.raises(RuntimeError, match="composition"):
        _run(pdb, region, eff, tmp_path / "out")


@pytest.mark.logic
def test_an_empty_receptor_chain_exits_one(tmp_path, monkeypatch):
    """The guard is unreachable through the file system today, because
    structure_metrics.get_chain_sequence raises when the chain is absent
    rather than returning an empty string, so the extraction is stubbed. It
    is pinned because the guard is what stands between a chain that read as
    empty and three FASTAs built from nothing.
    """
    pdb, region, eff = _stage(tmp_path)
    monkeypatch.setattr(bcs, "_extract_chain_sequence", lambda p, c: "")
    assert _run(pdb, region, eff, tmp_path / "out") == 1
    assert not (tmp_path / "out").exists()


@pytest.mark.logic
def test_an_index_file_with_no_indices_stops_the_stage(tmp_path):
    """A design region of nothing would write two controls identical to the
    design, which reads downstream as a control that passed every metric."""
    pdb, region, eff = _stage(tmp_path)
    region.write_text("# design region size: 0 residues\n")
    with pytest.raises(ValueError, match="no non-comment content"):
        _run(pdb, region, eff, tmp_path / "out")
