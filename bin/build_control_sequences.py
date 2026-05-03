#!/usr/bin/env python3
"""
build_control_sequences.py
--------------------------
Task 6 (P0-6) · Synthesise scrambled and polyA negative-control receptor
sequences for a given RFDiffusion design.

For each parent design we generate two control sequences:

  control_scrambled
      Deterministic seeded permutation of residues at design-region
      positions, preserving length and amino-acid composition (counts of
      each AA stay the same; only their order is shuffled).  Residues
      OUTSIDE the design region are held in place verbatim.  This isolates
      "what does Boltz do when the design region's chemistry is the same
      bag of amino acids in a random order?" — controls the predictor's
      bias toward sequence-pattern matching independent of geometric
      design intent.

  control_polyA
      Substitute every design-region position with alanine (A), EXCEPT
      where the source receptor already carries glycine (G) at that
      position — those stay G.  This isolates "what does Boltz do with
      a backbone whose design region has minimal side-chain chemistry?"
      — controls the predictor's reliance on the designable surface for
      its binding-site call.  The glycine carve-out is deliberate:
      replacing native G with A would force phi/psi into non-G regions
      and confound the geometric signal we want to test, whereas leaving
      native G alone keeps the backbone freedom intact.

Residues OUTSIDE the design region are held verbatim in BOTH controls.
The effector chain is passed through unchanged in both — controls test
the receptor binding face, not the effector.

Inputs
------
--source-pdb         RFDiffusion design PDB.  Receptor sequence is
                     extracted from the chain matching --receptor-chain.
--design-region-file Text file emitted by derive_design_region.py.  One
                     line of comma-separated 1-based positional indices
                     on the receptor chain (comments starting with '#'
                     are skipped).
--effector-fasta     Effector FASTA (single record).  Concatenated into
                     the output as the >effector chain unchanged.
--receptor-chain     Chain ID for receptor in source-pdb (default 'A').
--effector-chain     Chain ID for effector in source-pdb (default 'B';
                     accepted only for symmetry — not actually read.)
--output-dir         Directory to write control_scrambled.fasta and
                     control_polyA.fasta into.
--seed               Seed for the scrambled permutation.  Default is a
                     deterministic hash of the source PDB stem so re-
                     runs of the same design produce identical scrambled
                     sequences (essential for cohort reproducibility).

Output format
-------------
Three single-record FASTAs are written to ``--output-dir``, matching
the receptor.fasta / effector.fasta layout that
``negative_steering_run_one.sh`` expects via its
``--receptor-fasta`` / ``--effector-fasta`` flags:

    control_scrambled_receptor.fasta   — >receptor + scrambled sequence
    control_polyA_receptor.fasta       — >receptor + polyA sequence
    effector.fasta                     — >effector + original effector seq

The control_type ("scrambled" or "polyA") is intentionally NOT encoded
inside the FASTA — it is carried alongside the FASTA at the NextFlow
channel level and stamped into the per-control workdir as
``row_type.txt`` so cross_sequence_summary.py can mark each row.
"""
from __future__ import annotations

import argparse
import hashlib
import random
import sys
from pathlib import Path
from typing import List


# ── Sequence extraction (re-uses the project convention) ────────────

def _extract_chain_sequence(pdb_path: Path, chain_id: str) -> str:
    """
    Return the one-letter sequence for `chain_id` in `pdb_path`.

    Uses the same get_chain_sequence helper from boltz2_negative_steering
    that the rest of the pipeline uses, so chain-letter conventions and
    residue-name handling stay consistent across the codebase.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from boltz2_negative_steering import get_chain_sequence  # noqa: E402
    return get_chain_sequence(pdb_path, chain_id)


# ── Design-region parsing (matches derive_design_region.py output) ──

def _load_design_region_indices_1b(design_region_file: Path) -> List[int]:
    """
    Read the design-region file and return the 1-based positional
    indices on the receptor chain.  Format: comment lines starting
    with '#' are skipped, then a single line of comma-separated ints.
    """
    raw = design_region_file.read_text()
    payload = [
        ln.strip() for ln in raw.splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    if not payload:
        raise ValueError(
            f"{design_region_file}: no non-comment content found"
        )
    # Tolerate the indices being split across multiple lines, although
    # derive_design_region.py writes them on one.
    parts: List[str] = []
    for ln in payload:
        parts.extend(p.strip() for p in ln.split(",") if p.strip())
    try:
        return [int(p) for p in parts]
    except ValueError as e:
        raise ValueError(
            f"{design_region_file}: non-integer entry in indices list: {e}"
        ) from e


# ── Effector FASTA parsing (single-record, header optional) ────────

def _load_effector_sequence(effector_fasta: Path) -> str:
    """Concatenate every non-header line of a single-record FASTA."""
    seq_lines: List[str] = []
    with effector_fasta.open() as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">") or not line:
                continue
            seq_lines.append(line.strip())
    seq = "".join(seq_lines)
    if not seq:
        raise ValueError(f"{effector_fasta}: no sequence found")
    return seq


# ── Control synthesis ──────────────────────────────────────────────

def _build_scrambled_receptor(
    receptor_seq: str,
    design_region_idx_1b: List[int],
    seed: int,
) -> str:
    """
    Permute the residues at design-region positions, leave the rest
    untouched, preserve composition.
    """
    n = len(receptor_seq)
    chars = list(receptor_seq)

    # Convert 1-based positional indices to 0-based.  Filter any that
    # fall outside the receptor chain — defensive against any upstream
    # off-by-one (shouldn't happen but cheap to guard).
    design_idx_0b = sorted(
        i - 1 for i in design_region_idx_1b if 1 <= i <= n
    )
    if not design_idx_0b:
        # Nothing to scramble — return the receptor unchanged.  Caller
        # is responsible for deciding whether this is an error.
        return receptor_seq

    # Pull out the residues at the designated positions, shuffle them
    # under a deterministic Random(seed), put them back in the same
    # positions in the new order.
    rng = random.Random(seed)
    pool = [chars[i] for i in design_idx_0b]
    shuffled = pool[:]
    rng.shuffle(shuffled)
    for src_pos, new_residue in zip(design_idx_0b, shuffled):
        chars[src_pos] = new_residue
    return "".join(chars)


def _build_polyA_receptor(
    receptor_seq: str,
    design_region_idx_1b: List[int],
) -> str:
    """
    Substitute every design-region position with 'A', keeping native
    glycines as 'G' (per the v6 spec — alanine in a glycine-required
    backbone region would force phi/psi out of the only ramachandran
    well that fits, confounding the geometric signal we want to test).
    Residues outside the design region are unchanged.
    """
    n = len(receptor_seq)
    chars = list(receptor_seq)
    design_idx_0b = sorted(
        i - 1 for i in design_region_idx_1b if 1 <= i <= n
    )
    for i in design_idx_0b:
        if chars[i] == "G":
            continue   # preserve native glycine
        chars[i] = "A"
    return "".join(chars)


# ── Output writer ──────────────────────────────────────────────────

def _write_single_record_fasta(
    output_path: Path,
    record_name: str,
    sequence: str,
) -> None:
    """Write a single-record FASTA with one header line and one body line."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.write(f">{record_name}\n")
        f.write(sequence + "\n")


# ── Main ───────────────────────────────────────────────────────────

def _default_seed_from_pdb(source_pdb: Path) -> int:
    """
    Deterministic seed derived from the design PDB stem (e.g. 'design_3').
    Hashing ensures cohort-wide reproducibility — re-running the same
    design twice gives identical scrambled output, but two different
    designs in the same cohort get different shuffles.
    """
    stem = source_pdb.stem.encode("utf-8")
    h = hashlib.sha256(stem).digest()
    return int.from_bytes(h[:8], "big")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Synthesise scrambled and polyA negative-control "
                    "receptor sequences for a given design.",
    )
    ap.add_argument("--source-pdb", required=True, type=Path,
                    help="RFDiffusion design PDB to extract the receptor "
                         "sequence from.")
    ap.add_argument("--design-region-file", required=True, type=Path,
                    help="Text file (1-based comma-separated positional "
                         "indices) emitted by derive_design_region.py.")
    ap.add_argument("--effector-fasta", required=True, type=Path,
                    help="FASTA carrying the effector sequence to pass "
                         "through verbatim into both control FASTAs.")
    ap.add_argument("--receptor-chain", default="A",
                    help="Receptor chain ID in source-pdb (default A).")
    ap.add_argument("--effector-chain", default="B",
                    help="Effector chain ID (accepted for symmetry, not "
                         "actually read — effector comes from --effector-"
                         "fasta).")
    ap.add_argument("--output-dir", required=True, type=Path,
                    help="Output directory for control FASTAs.")
    ap.add_argument("--seed", type=int, default=None,
                    help="Seed for the scrambled permutation.  Defaults "
                         "to a SHA-256 hash of the source-pdb stem so "
                         "re-runs of the same design are reproducible.")
    args = ap.parse_args()

    # ── Load inputs ────────────────────────────────────────────────
    receptor_seq = _extract_chain_sequence(args.source_pdb, args.receptor_chain)
    if not receptor_seq:
        print(f"ERROR: empty receptor chain in {args.source_pdb}",
              file=sys.stderr)
        return 1
    effector_seq = _load_effector_sequence(args.effector_fasta)
    design_region_idx_1b = _load_design_region_indices_1b(
        args.design_region_file
    )

    seed = args.seed if args.seed is not None \
        else _default_seed_from_pdb(args.source_pdb)

    print(f"Source PDB:        {args.source_pdb}")
    print(f"Receptor:          {len(receptor_seq)} residues, "
          f"chain {args.receptor_chain}")
    print(f"Effector:          {len(effector_seq)} residues "
          f"(from {args.effector_fasta})")
    print(f"Design region:     {len(design_region_idx_1b)} positions")
    print(f"Scramble seed:     {seed}")

    # ── Build & write the two controls ────────────────────────────
    scrambled_seq = _build_scrambled_receptor(
        receptor_seq, design_region_idx_1b, seed,
    )
    polyA_seq = _build_polyA_receptor(
        receptor_seq, design_region_idx_1b,
    )

    # Sanity: scrambled must preserve length and composition exactly.
    if len(scrambled_seq) != len(receptor_seq):
        raise RuntimeError(
            f"scrambled length {len(scrambled_seq)} != "
            f"original {len(receptor_seq)}"
        )
    from collections import Counter
    if Counter(scrambled_seq) != Counter(receptor_seq):
        raise RuntimeError("scrambled composition diverged from original")

    # Sanity: polyA must change at least one position when the design
    # region is non-empty.  If every design-region position was already
    # G, the substitution is a no-op — flag that loudly because it
    # means the design region is degenerate and the polyA control
    # carries no signal.
    if polyA_seq == receptor_seq and design_region_idx_1b:
        print(
            "WARNING: polyA control is identical to source receptor — "
            "every design-region position is already glycine.  Control "
            "carries no signal for this design.",
            file=sys.stderr,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scrambled_path = args.output_dir / "control_scrambled_receptor.fasta"
    polyA_path     = args.output_dir / "control_polyA_receptor.fasta"
    effector_path  = args.output_dir / "effector.fasta"
    _write_single_record_fasta(scrambled_path, "receptor", scrambled_seq)
    _write_single_record_fasta(polyA_path,     "receptor", polyA_seq)
    _write_single_record_fasta(effector_path,  "effector", effector_seq)

    # Quick diff summary so the log shows what changed.
    n_scrambled_diff = sum(
        1 for a, b in zip(receptor_seq, scrambled_seq) if a != b
    )
    n_polyA_diff = sum(
        1 for a, b in zip(receptor_seq, polyA_seq) if a != b
    )
    print(f"Wrote {scrambled_path.name}: "
          f"{n_scrambled_diff} positions changed vs source "
          f"(of {len(design_region_idx_1b)} design-region positions)")
    print(f"Wrote {polyA_path.name}: "
          f"{n_polyA_diff} positions changed vs source "
          f"(of {len(design_region_idx_1b)} design-region positions; "
          f"native glycines preserved)")
    print(f"Wrote {effector_path.name}: "
          f"effector passed through unchanged ({len(effector_seq)} residues)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())