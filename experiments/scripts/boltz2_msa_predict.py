"""Boltz-2 with MSA: re-predict one negsteer design for diagnostic comparison.

Runs Boltz-2 with a real MMseqs2 / ColabFold MSA on the steered receptor
sequence + native effector sequence from a single negsteer per-sequence
workdir, then computes the same confidence metrics the production
pipeline computes (compute_metrics.py) and the receptor-aligned
effector RMSD (binding_rmsds()).

This is a diagnostic prediction whose ONLY difference from the
pipeline's no-MSA negsteer prediction is the receptor MSA. Everything
else — the effector template CIF anchoring the effector fold, the
chain layout, the Boltz invocation flags — matches the pipeline.
The point is to isolate the effect of the MSA on Boltz-2's confidence
output for a sequence we already have a no-MSA prediction for.

Compared to the pipeline:
  * Receptor chain gets a real MMseqs2 / ColabFold A3M (mirroring
    experiments/scripts/msa.nf: --use-env 1 --use-templates 0
    --prefilter-mode 1, MMSEQS_IGNORE_INDEX=1).
  * Effector chain still gets a single-sequence A3M AND the effector
    template CIF (extract_effector_template_cif), exactly as
    bin/boltz2_negative_steering.py composes its YAML.
  * Boltz invocation mirrors experiments/scripts/boltz2.nf BOLTZ2_MSA
    (--diffusion_samples 5 --sampling_steps 20 --output_format pdb
    --write_full_pae --use_potentials --no_kernels --override
    --num_workers 0).  --recycling_steps defaults to 3 here (the
    negsteer pipeline default; see boltz2_negative_steering.py
    parser at line 3296), not 20 as in the BOLTZ2_MSA benchmark — we
    want the Boltz settings as close to the negsteer no-MSA run as
    possible so the only delta is the MSA.
  * Metrics use bin/compute_metrics.py via subprocess, exactly as
    boltz2_iterate_steering.py does.
  * ra_eff uses binding_rmsds() from bin/boltz2_negative_steering.py.

Container and reference-data paths are HPC environment constants
sourced from nextflow.config (the canonical declarations) — see the
HPC CONSTANTS block below.  --colabfold-db is the only such path
that has no canonical definition in any committed config (the
benchmark pipeline takes it as a Nextflow param too) and so is
required on the CLI.

Output layout under --outdir:
    msa/<design_id>_receptor.a3m
    inputs/effector_template.cif
    inputs/<design_id>_seed<N>/input.yaml          (+ msa/ siblings)
    predictions/<design_id>/seed<N>/...
    metrics/<design_id>_seed<N>_metrics.json
    <design_id>_summary.json
"""

from __future__ import annotations

import experiments._path_setup  # noqa: F401  - adds bin/ to sys.path

import argparse
import csv
import json
import logging
import os
import shutil
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from boltz2_negative_steering import (
    PRED_EFF_CHAIN,
    PRED_REC_CHAIN,
    binding_rmsds,
    extract_effector_template_cif,
    write_boltz_yaml,
)

LOG = logging.getLogger("boltz2_msa_predict")


# ═══════════════════════════════════════════════════════════════════════
# HPC CONSTANTS
# ───────────────────────────────────────────────────────────────────────
# Container and Boltz CLI defaults sourced from nextflow.config (the
# canonical declarations) and experiments/scripts/boltz2.nf (the
# canonical BOLTZ2_MSA invocation).  These are HPC environment
# constants, not per-run choices, so they live here rather than on the
# CLI.  Override with the corresponding env var if you ever need to
# point at a different container without editing the script.
# ═══════════════════════════════════════════════════════════════════════
# Containers are referenced via symlinks in the repo's containers/ dir
# (see containers/README.md).  The symlinks resolve to the real .img under
# the user's $HOME, which Singularity auto-mounts — so they also resolve
# from inside the boltz2 container when this script launches colabfold
# nested.  parents[2] == repo root (experiments/scripts/<this file>).
_REPO_ROOT = Path(__file__).resolve().parents[2]
BOLTZ2_CONTAINER = Path(os.environ.get(
    "BOLTZ2_CONTAINER",
    str(_REPO_ROOT / "containers" / "boltz2_negsteer.img"),
))
COLABFOLD_CONTAINER = Path(os.environ.get(
    "COLABFOLD_CONTAINER",
    str(_REPO_ROOT / "containers" / "colabfold.img"),
))

# Boltz invocation defaults — match BOLTZ2_MSA in
# experiments/scripts/boltz2.nf so this script's predictions are
# directly comparable to the benchmark MSA run.
BOLTZ_DIFFUSION_SAMPLES = 5     # boltz2.nf:87
BOLTZ_SAMPLING_STEPS    = 20    # boltz2.nf:88
# Recycling matches the negsteer no-MSA run (boltz2_negative_steering.py
# parser default, line 3296), not BOLTZ2_MSA's 20 — we want only the
# MSA to vary between the two runs we are comparing.
BOLTZ_RECYCLING_STEPS_DEFAULT = 3


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Re-predict one negsteer design with Boltz-2 + real MMseqs2 MSA "
            "and emit the same confidence metrics the pipeline computes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Container and Boltz invocation defaults are HPC constants; "
            "see the HPC CONSTANTS block in this script. Set "
            "BOLTZ2_CONTAINER / COLABFOLD_CONTAINER env vars to override.\n"
            "\n"
            "Stages:\n"
            "  --stage msa      MSA search only (CPU partition).  Writes\n"
            "                   <outdir>/msa/<design>_receptor.a3m and\n"
            "                   <outdir>/msa_done.txt.\n"
            "  --stage predict  Boltz + metrics only (GPU partition).\n"
            "                   Reads the A3M path from --msa-a3m or from\n"
            "                   <outdir>/msa_done.txt.\n"
            "  --stage both     Run everything sequentially (default;\n"
            "                   convenient for local testing only)."
        ),
    )
    p.add_argument("--stage", choices=("msa", "predict", "both"),
                   default="both",
                   help="Which stage to run (default: both).  See epilog.")
    p.add_argument("--workdir", required=True, type=Path,
                   help="Path to one published negsteer per-sequence workdir, "
                        "e.g. <results>/negative_steering/runs/design_3_seq_0/. "
                        "Must contain cycle_0/steered_results.csv and "
                        "cycle_0/steered/<design>/receptor.fasta.")
    p.add_argument("--reference-pdb", type=Path,
                   help="Ground-truth complex PDB. Required for --stage "
                        "predict / both: used for ra_eff and for "
                        "extracting the effector template CIF.  Not "
                        "needed for --stage msa.")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Directory to write outputs into. Created if absent.")
    p.add_argument("--colabfold-db", type=Path,
                   help="Path to the ColabFold MMseqs2 sequence database "
                        "directory.  Required for --stage msa / both.")
    p.add_argument("--msa-a3m", type=Path, default=None,
                   help="Path to a pre-computed receptor A3M.  Only used "
                        "with --stage predict.  Defaults to the A3M path "
                        "recorded in <outdir>/msa_done.txt.")
    p.add_argument("--seeds", nargs="+", type=int, default=None,
                   help="Boltz seeds to predict. Default: read the "
                        "boltz_seed values for the chosen steered design's "
                        "sequence_group from cycle_0/plan.json.")
    p.add_argument("--n-recycling", type=int,
                   default=BOLTZ_RECYCLING_STEPS_DEFAULT,
                   help=f"Boltz --recycling_steps "
                        f"(default: {BOLTZ_RECYCLING_STEPS_DEFAULT}, the "
                        f"negsteer no-MSA default).")
    p.add_argument("--truth-rec-chain", default="A",
                   help="Receptor chain ID in --reference-pdb (default: A).")
    p.add_argument("--truth-eff-chain", default="B",
                   help="Effector chain ID in --reference-pdb (default: B).")
    p.add_argument("--threads", type=int, default=10,
                   help="Threads for colabfold_search (default: 10). "
                        "Match this to the SLURM allocation.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose (DEBUG) logging.")
    args = p.parse_args(argv)

    if args.stage in ("msa", "both") and args.colabfold_db is None:
        p.error("--colabfold-db is required for --stage msa / both")
    if args.stage in ("predict", "both") and args.reference_pdb is None:
        p.error("--reference-pdb is required for --stage predict / both")
    return args


# ═══════════════════════════════════════════════════════════════════════
# Steered sequence resolution
# ═══════════════════════════════════════════════════════════════════════
def resolve_steered_inputs(
    workdir: Path,
) -> Tuple[str, str, str, List[int]]:
    """Resolve the steered receptor sequence and seeds for one published
    negsteer per-sequence workdir using the pipeline's own representative
    pick — NOT just the lowest-ra_eff row in steered_results.csv.

    Why this matters: a per-sequence workdir contains MULTIPLE distinct
    steered sequences (one per sequence_group, each a different mutation
    pattern).  Each sequence is predicted with multiple seeds.
    steered_results.csv ranks individual (sequence_group, seed_index)
    rows by raw ra_eff, so its rank=1 row picks both the best sequence
    AND its luckiest seed.  The pipeline does NOT take that row as the
    representative — it applies a tier-then-composite policy across
    sequence_groups (see boltz2_iterate_steering.cmd_aggregate_per_sequence
    and bin/cross_sequence_summary._tier_for_row +
    _pick_tier_none_fallback).  Reading the wrong row here would compare
    the MSA prediction of an unrelated sequence to the no-MSA medians of
    a different sequence — a category-error invalidating the whole
    diagnostic.

    Algorithm — mirrors cross_sequence_summary.py exactly:

      1. If passing_summary.csv has rows, the workdir is tier A/B/C.
         Pick the row with rank_by_composite_score == 1 (the
         cross-sequence summary's representative pick).

      2. Otherwise it is tier-none.  Read aggregated_results.csv and
         pick the row with the maximum composite score
            composite = steered_true_jaccard_median
                        - 0.05 × (reverted_ra_eff_vs_truth_median if
                                  aggregated_verdict == "pose_holds"
                                  else steered_ra_eff_vs_truth_median)
         (matches _pick_tier_none_fallback in cross_sequence_summary.py
         lines 350-440).

    The chosen row carries a base ``design`` name (e.g. ``design_00``,
    with the ``_sX`` seed suffix already stripped by
    cmd_aggregate_per_sequence) and a ``canonical_seed_index``.  All
    seeds in a sequence_group share the same steered receptor sequence
    (only the Boltz seed differs), so we glob
    cycle_0/steered/<design>* / receptor.fasta and read the first match.

    Layout (verified against
    tests/negative_steering/example_output_files/negative_steering/runs/
    design_3_seq_1/):

        <workdir>/                                 ← published per-seq workdir
          aggregated_results.csv                   ← one row per sequence_group
          passing_summary.csv                      ← passing rows only
          inputs/
            receptor.fasta                         ← MPNN receptor (pre-steering)
            effector.fasta                         ← native effector sequence
          cycle_0/
            plan.json                              ← seeds, design metadata
            steered_results.csv                    ← per (sg, seed) raw ra_eff rank
            steered/
              <design_NN[_sX]>/                    ← all seeds in a group share
                receptor.fasta                     ←   the same receptor.fasta
                ...

    PATH_UNCERTAIN: this assumes the published per-sequence workdir is
    intact and matches the layout in the test fixture.  Errors are
    raised with the exact missing path so the user can fix their input.
    """
    if not workdir.is_dir():
        raise FileNotFoundError(f"negsteer workdir not found: {workdir}")

    cycle0 = workdir / "cycle_0"
    if not cycle0.is_dir():
        raise FileNotFoundError(
            f"cycle_0 directory not found under {workdir} — this does not "
            "look like a negsteer per-sequence workdir."
        )

    plan_path = cycle0 / "plan.json"
    if not plan_path.is_file():
        raise FileNotFoundError(
            f"plan.json absent: {plan_path} — cannot recover seeds list."
        )
    plan = json.loads(plan_path.read_text())

    chosen_row = _pick_representative_row(workdir)
    design_name = (chosen_row.get("design") or "").strip()
    sg_str = str(chosen_row.get("sequence_group") or "").strip()
    canonical_seed_index = str(
        chosen_row.get("canonical_seed_index") or ""
    ).strip()
    if not design_name:
        raise ValueError(
            f"Representative row from {workdir} has empty 'design' field; "
            "cannot locate receptor.fasta."
        )

    # PATH_UNCERTAIN: effector.fasta is in <workdir>/inputs/ in the test
    # fixture (split out of the upstream MPNN combined FASTA by the
    # NEGSTEER_RUN_ONE wrapper, modules/negative_steering.nf:166).
    eff_fasta = workdir / "inputs" / "effector.fasta"
    if not eff_fasta.is_file():
        raise FileNotFoundError(
            f"Effector FASTA absent: {eff_fasta} — expected the negsteer "
            "process to have written this when splitting the MPNN combined "
            "FASTA."
        )
    eff_seq = _read_single_fasta_seq(eff_fasta)

    # Special case: cold_start_all_clean / skip_steering=True
    # ------------------------------------------------------
    # When the wild-type prediction was already correct, negsteer
    # writes plan.json with skip_steering=True, no steered/ subdir is
    # produced, and aggregated_results.csv carries one row whose
    # 'design' is literally "initial".  cross_sequence_summary picks
    # this row as the representative (verdict="no_reversion"), so a
    # diagnostic comparison still has to use it.  In that case the
    # "steered" sequence is the wild-type MPNN sequence (the design
    # never needed steering), and the seeds come from cold_start_seeds
    # in plan.json rather than the empty designs[] list.
    if design_name == "initial":
        rec_seq = (plan.get("wild_type_receptor_seq") or "").strip()
        if not rec_seq:
            raise FileNotFoundError(
                f"plan.json at {plan_path} lacks wild_type_receptor_seq, "
                "but the representative is the cold-start 'initial' row. "
                "Cannot recover the wild-type receptor sequence."
            )
        cs_seeds = plan.get("cold_start_seeds") or []
        seeds = sorted({
            int(s["boltz_seed"]) for s in cs_seeds
            if s.get("boltz_seed") is not None
        })
        if not seeds:
            raise FileNotFoundError(
                f"plan.json at {plan_path} has skip_steering=True but no "
                "cold_start_seeds entries — no seeds available to predict."
            )
        return "initial", rec_seq, eff_seq, seeds

    # Regular steered case: cycle_0/steered/<design>[_s<seed>]/receptor.fasta
    # — try the seed-specific dir first (multi-seed runs), fall back to
    # the bare base name (single-seed runs), then to any match across
    # seeds (all seeds in a group share the same FASTA).
    candidates = []
    if canonical_seed_index:
        candidates.append(
            cycle0 / "steered" / f"{design_name}_s{canonical_seed_index}"
            / "receptor.fasta"
        )
    candidates.append(cycle0 / "steered" / design_name / "receptor.fasta")
    rec_fasta: Optional[Path] = next(
        (c for c in candidates if c.is_file()),
        None,
    )
    if rec_fasta is None:
        # Glob for any seed of this base design.  If num_seeds > 1 the
        # earlier _s<seed_index> path catches the canonical seed; this
        # fallback handles odd numbering schemes.
        steered_root = cycle0 / "steered"
        if steered_root.is_dir():
            glob_matches = sorted(
                steered_root.glob(f"{design_name}*/receptor.fasta")
            )
            if glob_matches:
                rec_fasta = glob_matches[0]
    if rec_fasta is None or not rec_fasta.is_file():
        searched = "; ".join(str(c) for c in candidates)
        raise FileNotFoundError(
            f"Steered receptor FASTA absent for representative design "
            f"'{design_name}' (sequence_group={sg_str}, "
            f"canonical_seed_index={canonical_seed_index!r}). "
            f"Looked at: {searched}"
        )

    rec_seq = _read_single_fasta_seq(rec_fasta)

    # Recover the seeds the negsteer run used for this design's
    # sequence_group.  plan.json's "designs" list is the authoritative
    # record of (sequence_group, seed_index, boltz_seed).
    seeds: List[int] = []
    for d in plan.get("designs", []):
        d_sg = str(d.get("sequence_group", "")).strip()
        if d_sg == sg_str:
            try:
                seeds.append(int(d.get("boltz_seed")))
            except (TypeError, ValueError):
                continue
    seeds = sorted(set(seeds))
    if not seeds:
        LOG.warning(
            "Could not match representative sequence_group=%r against "
            "any plan.json design entry; using all plan seeds as fallback.",
            sg_str,
        )
        for d in plan.get("designs", []):
            try:
                seeds.append(int(d.get("boltz_seed")))
            except (TypeError, ValueError):
                continue
        seeds = sorted(set(seeds))

    # Identifier returned for output naming includes the canonical seed
    # suffix when present, matching the on-disk steered/ directory name.
    label = (
        f"{design_name}_s{canonical_seed_index}"
        if canonical_seed_index else design_name
    )
    return label, rec_seq, eff_seq, seeds


def _pick_representative_row(workdir: Path) -> Dict[str, str]:
    """Apply the pipeline's tier-then-composite representative policy
    to one per-sequence workdir.  Mirrors cross_sequence_summary.py:

      * passing_summary.csv non-empty → row with rank_by_composite_score == 1
      * passing_summary.csv empty (tier-none) → aggregated_results.csv
        row with max composite (true_jaccard - 0.05 * ra_eff,
        ra_eff sourced from reverted_* if aggregated_verdict == 'pose_holds',
        else steered_*)
    """
    passing_csv = workdir / "passing_summary.csv"
    if passing_csv.is_file():
        with passing_csv.open() as f:
            rows = list(csv.DictReader(f))
        if rows:
            for r in rows:
                try:
                    if int(float(r.get("rank_by_composite_score") or "")) == 1:
                        return r
                except (TypeError, ValueError):
                    continue
            # No rank=1 found but rows exist → take first row (defensive)
            LOG.warning(
                "passing_summary.csv at %s has rows but none have "
                "rank_by_composite_score == 1; falling back to first row.",
                passing_csv,
            )
            return rows[0]

    # Tier-none fallback — pick max composite from aggregated_results.csv.
    agg_csv = workdir / "aggregated_results.csv"
    if not agg_csv.is_file():
        raise FileNotFoundError(
            f"Neither passing_summary.csv nor aggregated_results.csv "
            f"under {workdir} — cannot pick a representative steered "
            "design."
        )
    with agg_csv.open() as f:
        agg_rows = list(csv.DictReader(f))
    if not agg_rows:
        raise FileNotFoundError(
            f"aggregated_results.csv under {workdir} is empty (header "
            "only); workdir produced no aggregated rows."
        )

    best_row: Optional[Dict[str, str]] = None
    best_composite: Optional[float] = None
    for r in agg_rows:
        # Skip the cold-start "initial" singleton row — it has no
        # sequence_group and no steered metrics.
        if (r.get("aggregated_verdict") or "").strip() == "singleton":
            continue
        if not (r.get("sequence_group") or "").strip():
            continue
        is_pose_holds = (r.get("aggregated_verdict") or "").strip() == "pose_holds"
        ra_key = (
            "reverted_ra_eff_vs_truth_median"
            if is_pose_holds
            else "steered_ra_eff_vs_truth_median"
        )
        ra = _try_float_local(r.get(ra_key))
        tj = _try_float_local(r.get("steered_true_jaccard_median"))
        if ra is None or tj is None:
            continue
        composite = tj - 0.05 * ra
        if best_composite is None or composite > best_composite:
            best_composite = composite
            best_row = r

    if best_row is None:
        raise FileNotFoundError(
            f"No aggregated_results.csv row under {workdir} has both "
            "steered_true_jaccard_median and a finite ra_eff median — "
            "cannot pick a tier-none representative."
        )
    return best_row


def _try_float_local(v) -> Optional[float]:
    if v is None or v == "":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if f != f:
        return None
    return f


def _read_single_fasta_seq(path: Path) -> str:
    lines = path.read_text().splitlines()
    body = [ln.strip() for ln in lines if ln and not ln.startswith(">")]
    seq = "".join(body)
    if not seq:
        raise ValueError(f"No sequence content in FASTA: {path}")
    return seq


# ═══════════════════════════════════════════════════════════════════════
# MSA generation (mirrors experiments/scripts/msa.nf)
# ═══════════════════════════════════════════════════════════════════════
def generate_receptor_msa(
    rec_seq: str,
    out_a3m: Path,
    colabfold_db: Path,
    threads: int,
) -> Path:
    """Run colabfold_search to produce a single chain_A.a3m for the
    receptor sequence.  Mirrors the COLABFOLD_SEARCH process in
    experiments/scripts/msa.nf — same flags, same env var, same use of
    the FASTA header to drive the output filename (>chain_A → chain_A.a3m).
    """
    out_a3m.parent.mkdir(parents=True, exist_ok=True)

    work_dir = out_a3m.parent / f".{out_a3m.stem}_workdir"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    fasta_path = work_dir / "input.fasta"
    fasta_path.write_text(f">chain_A\n{rec_seq}\n")

    msa_out_dir = work_dir / "msa"
    msa_out_dir.mkdir(exist_ok=True)

    if not colabfold_db.is_dir():
        raise FileNotFoundError(
            f"ColabFold MMseqs2 database not found: {colabfold_db}"
        )

    cmd = [
        "singularity", "exec",
        "--bind", f"{work_dir.resolve()}:{work_dir.resolve()}",
        "--bind", f"{colabfold_db.resolve()}:{colabfold_db.resolve()}",
        str(COLABFOLD_CONTAINER),
        "env", "MMSEQS_IGNORE_INDEX=1",
        "colabfold_search",
        "--mmseqs", "mmseqs",
        "--use-env", "1",
        "--use-templates", "0",
        "--threads", str(threads),
        "--prefilter-mode", "1",
        str(fasta_path),
        str(colabfold_db),
        str(msa_out_dir),
    ]
    LOG.info("Running colabfold_search: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(
            f"colabfold_search failed (exit {proc.returncode}). "
            "See stderr above."
        )

    produced = msa_out_dir / "chain_A.a3m"
    if not produced.is_file():
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(
            f"colabfold_search reported success but {produced} is missing. "
            f"Tree of {msa_out_dir}: "
            f"{[p.name for p in msa_out_dir.iterdir()]}"
        )

    shutil.copyfile(produced, out_a3m)
    LOG.info("Wrote receptor MSA: %s (%d lines)",
             out_a3m, sum(1 for _ in out_a3m.open()))
    return out_a3m


# ═══════════════════════════════════════════════════════════════════════
# Boltz invocation (mirrors experiments/scripts/boltz2.nf BOLTZ2_MSA)
# ═══════════════════════════════════════════════════════════════════════
def run_boltz_msa(
    yaml_path: Path,
    out_dir: Path,
    seed: int,
    n_recycling: int,
) -> Path:
    """Invoke boltz predict (Boltz-2) on a single YAML, mirroring the
    BOLTZ2_MSA process flags in experiments/scripts/boltz2.nf.

    Returns the path to a chosen prediction PDB (the model_0 / rank_0
    favoured one) so callers can locate sidecars relative to its parent.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    boltz_args = [
        "boltz", "predict",
        str(yaml_path),
        "--out_dir", str(out_dir),
        "--recycling_steps", str(n_recycling),
        "--diffusion_samples", str(BOLTZ_DIFFUSION_SAMPLES),
        "--sampling_steps", str(BOLTZ_SAMPLING_STEPS),
        "--seed", str(seed),
        "--num_workers", "0",
        "--output_format", "pdb",
        "--write_full_pae",
        "--use_potentials",
        "--no_kernels",
        "--override",
    ]

    cmd = [
        "singularity", "exec", "--nv",
        "--bind", f"{out_dir.resolve()}:{out_dir.resolve()}",
        "--bind", f"{yaml_path.parent.resolve()}:{yaml_path.parent.resolve()}",
        str(BOLTZ2_CONTAINER),
    ] + boltz_args
    LOG.info("Running boltz predict (seed=%d): %s", seed, " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(
            f"boltz predict (seed={seed}) failed exit={proc.returncode}"
        )

    pdbs = sorted(p for p in out_dir.rglob("*.pdb") if "msa" not in p.parts)
    if not pdbs:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(
            f"boltz reported success but produced no PDBs under {out_dir}"
        )
    pdbs.sort(key=lambda p: (
        "model_0" not in p.name and "rank_0" not in p.name,
        p.name,
    ))
    return pdbs[0]


# ═══════════════════════════════════════════════════════════════════════
# compute_metrics.py invocation (mirrors boltz2_iterate_steering)
# ═══════════════════════════════════════════════════════════════════════
def run_compute_metrics(
    pred_pdb: Path,
    chain_lengths: Tuple[int, int],
    out_csv: Path,
) -> Optional[Dict[str, str]]:
    """Invoke bin/compute_metrics.py on the prediction directory of a
    single Boltz output.  Returns the first parsed CSV row as a dict,
    or None on failure.

    Mirrors the subprocess.run pattern in
    boltz2_iterate_steering.cmd_compute_final_metrics: --model boltz2,
    --prediction-dir = the parent of the model PDB (so sidecars are
    picked up), --chain-lengths from the model PDB, default
    --pae-cutoff 10, default --contact-cutoff 5.0, no
    --mutated-positions (this is an unconstrained re-prediction).
    """
    compute_metrics_script = (
        Path(__file__).resolve().parents[2] / "bin" / "compute_metrics.py"
    )
    if not compute_metrics_script.is_file():
        raise FileNotFoundError(
            f"bin/compute_metrics.py not found at {compute_metrics_script}"
        )

    pred_dir = pred_pdb.parent
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rec_len, eff_len = chain_lengths

    inner = [
        "python", str(compute_metrics_script),
        "--model", "boltz2",
        "--prediction-dir", str(pred_dir),
        "--chain-lengths", str(rec_len), str(eff_len),
        "--output-csv", str(out_csv),
        "--receptor-chain", PRED_REC_CHAIN,
        "--effector-chain", PRED_EFF_CHAIN,
    ]
    cmd = [
        "singularity", "exec",
        "--bind", f"{pred_dir.resolve()}:{pred_dir.resolve()}",
        "--bind", f"{compute_metrics_script.parent.resolve()}:"
                  f"{compute_metrics_script.parent.resolve()}",
        "--bind", f"{out_csv.parent.resolve()}:{out_csv.parent.resolve()}",
        str(BOLTZ2_CONTAINER),
    ] + inner
    LOG.info("Running compute_metrics: %s", " ".join(cmd))
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        LOG.warning("compute_metrics.py timed out for %s", pred_pdb)
        return None
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        LOG.warning(
            "compute_metrics.py rc=%d for %s: %s",
            proc.returncode, pred_pdb, (proc.stderr or "").strip()[:300],
        )
        return None

    if not out_csv.is_file():
        LOG.warning("compute_metrics.py produced no CSV at %s", out_csv)
        return None

    with out_csv.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        LOG.warning("compute_metrics.py CSV empty: %s", out_csv)
        return None
    return rows[0]


def chain_lengths_from_pdb(pdb_path: Path) -> Tuple[int, int]:
    """Count CA atoms per chain to derive --chain-lengths for compute_metrics."""
    counts: Dict[str, int] = {}
    with pdb_path.open() as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ch = line[21:22].strip() or " "
                counts[ch] = counts.get(ch, 0) + 1
    rec_len = counts.get(PRED_REC_CHAIN, 0)
    eff_len = counts.get(PRED_EFF_CHAIN, 0)
    if rec_len == 0 or eff_len == 0:
        ordered = sorted(counts.items(), key=lambda kv: -kv[1])
        if len(ordered) >= 2:
            rec_len = rec_len or ordered[0][1]
            eff_len = eff_len or ordered[1][1]
    if rec_len == 0 or eff_len == 0:
        raise RuntimeError(
            f"Could not determine receptor/effector chain lengths from {pdb_path}; "
            f"chain CA counts: {counts}"
        )
    return rec_len, eff_len


# ═══════════════════════════════════════════════════════════════════════
# Per-seed driver
# ═══════════════════════════════════════════════════════════════════════
def _try_float(v) -> Optional[float]:
    if v is None or v == "":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if f != f:
        return None
    return f


def _median(vals: List[float]) -> Optional[float]:
    finite = [v for v in vals if v is not None]
    if not finite:
        return None
    return float(statistics.median(finite))


def predict_one_seed(
    *,
    design_id: str,
    rec_seq: str,
    eff_seq: str,
    receptor_msa_a3m: Path,
    effector_template_cif: Path,
    seed: int,
    outdir: Path,
    reference_pdb: Path,
    truth_rec_chain: str,
    truth_eff_chain: str,
    n_recycling: int,
) -> Dict[str, object]:
    """Run one (design, seed) prediction and return a per-seed result dict."""
    seed_input_dir = outdir / "inputs" / f"{design_id}_seed{seed}"
    seed_input_dir.mkdir(parents=True, exist_ok=True)

    yaml_path = write_boltz_yaml(
        seed_input_dir,
        rec_seq,
        eff_seq,
        PRED_REC_CHAIN,
        PRED_EFF_CHAIN,
        effector_template_cif=effector_template_cif,
    )
    # write_boltz_yaml writes single-seq A3Ms for both chains; replace
    # the receptor's with the real MMseqs2 MSA.  The effector A3M stays
    # single-sequence — the effector chain is anchored by the template
    # CIF, which is the pipeline's convention (the effector is fixed
    # across the experiment so its MSA adds nothing).
    receptor_a3m = seed_input_dir / "msa" / f"chain_{PRED_REC_CHAIN}.a3m"
    shutil.copyfile(receptor_msa_a3m, receptor_a3m)

    pred_dir = outdir / "predictions" / design_id / f"seed{seed}"
    pred_dir.mkdir(parents=True, exist_ok=True)

    seed_result: Dict[str, object] = {
        "seed": seed,
        "yaml": str(yaml_path),
        "prediction_dir": str(pred_dir),
        "status": "ok",
        "error": None,
    }

    try:
        pred_pdb = run_boltz_msa(
            yaml_path=yaml_path,
            out_dir=pred_dir,
            seed=seed,
            n_recycling=n_recycling,
        )
    except Exception as e:
        LOG.error("Boltz prediction failed for %s seed=%d: %s",
                  design_id, seed, e)
        seed_result["status"] = "failed"
        seed_result["error"] = f"boltz: {e}"
        return seed_result

    seed_result["pred_pdb"] = str(pred_pdb)

    try:
        rec_len, eff_len = chain_lengths_from_pdb(pred_pdb)
    except Exception as e:
        LOG.error("Chain-length detection failed for %s seed=%d: %s",
                  design_id, seed, e)
        seed_result["status"] = "failed"
        seed_result["error"] = f"chain_lengths: {e}"
        return seed_result

    metrics_csv = outdir / "metrics" / f"{design_id}_seed{seed}_metrics.csv"
    metrics_row = run_compute_metrics(
        pred_pdb=pred_pdb,
        chain_lengths=(rec_len, eff_len),
        out_csv=metrics_csv,
    )
    if metrics_row is None:
        seed_result["status"] = "failed"
        seed_result["error"] = "compute_metrics: no row produced"
        return seed_result

    metrics_pickup = {
        "complex_plddt":   _try_float(metrics_row.get("complex_plddt")),
        "interface_plddt": _try_float(metrics_row.get("interface_plddt")),
        "ipae":            _try_float(metrics_row.get("ipae")),
        "pae_pass_frac":   _try_float(metrics_row.get("pae_pass_frac")),
        "iptm":            _try_float(metrics_row.get("iptm")),
        "ptm":             _try_float(metrics_row.get("ptm")),
        "actifptm":        _try_float(metrics_row.get("actifptm")),
        "avg_plddt":       _try_float(metrics_row.get("avg_plddt")),
    }

    try:
        ra_eff, ind_rec, ind_eff = binding_rmsds(
            pred_pdb=pred_pdb,
            truth_pdb=reference_pdb,
            pred_rec_chain=PRED_REC_CHAIN,
            pred_eff_chain=PRED_EFF_CHAIN,
            truth_rec_chain=truth_rec_chain,
            truth_eff_chain=truth_eff_chain,
        )
    except Exception as e:
        LOG.warning("ra_eff failed for %s seed=%d: %s",
                    design_id, seed, e)
        ra_eff = ind_rec = ind_eff = None

    seed_result["metrics"] = metrics_pickup
    seed_result["metrics_csv"] = str(metrics_csv)
    seed_result["ra_eff"] = (
        None if ra_eff is None or ra_eff != ra_eff else float(ra_eff)
    )
    seed_result["independent_receptor_rmsd"] = (
        None if ind_rec is None or ind_rec != ind_rec else float(ind_rec)
    )
    seed_result["independent_effector_rmsd"] = (
        None if ind_eff is None or ind_eff != ind_eff else float(ind_eff)
    )

    metrics_json = outdir / "metrics" / f"{design_id}_seed{seed}_metrics.json"
    metrics_json.write_text(json.dumps({
        "design_id": design_id,
        "seed": seed,
        "metrics": metrics_pickup,
        "ra_eff": seed_result["ra_eff"],
        "independent_receptor_rmsd": seed_result["independent_receptor_rmsd"],
        "independent_effector_rmsd": seed_result["independent_effector_rmsd"],
        "raw_compute_metrics_row": metrics_row,
    }, indent=2))
    seed_result["metrics_json"] = str(metrics_json)
    return seed_result


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════
def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    workdir = args.workdir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    design_id = workdir.name
    LOG.info("Stage=%s for %s", args.stage, design_id)

    steered_design_name, rec_seq, eff_seq, plan_seeds = (
        resolve_steered_inputs(workdir)
    )
    LOG.info("Steered design name in negsteer plan: %s", steered_design_name)
    LOG.info("Receptor sequence length: %d", len(rec_seq))
    LOG.info("Effector sequence length: %d", len(eff_seq))

    receptor_msa_a3m = outdir / "msa" / f"{design_id}_receptor.a3m"
    msa_done_path = outdir / "msa_done.txt"

    # ── Stage: msa ────────────────────────────────────────────────────
    if args.stage in ("msa", "both"):
        if not receptor_msa_a3m.is_file():
            generate_receptor_msa(
                rec_seq=rec_seq,
                out_a3m=receptor_msa_a3m,
                colabfold_db=args.colabfold_db.resolve(),
                threads=args.threads,
            )
        else:
            LOG.info("Reusing existing MSA at %s", receptor_msa_a3m)
        msa_done_path.write_text(f"{receptor_msa_a3m.resolve()}\n")
        LOG.info("Wrote %s -> %s", msa_done_path, receptor_msa_a3m)

    if args.stage == "msa":
        LOG.info("Stage=msa complete for %s; exiting before Boltz step.",
                 design_id)
        return 0

    # ── Stage: predict (or both, continued) ───────────────────────────
    reference_pdb = args.reference_pdb.resolve()
    if not reference_pdb.is_file():
        raise FileNotFoundError(f"reference PDB not found: {reference_pdb}")

    seeds = args.seeds if args.seeds else plan_seeds
    if not seeds:
        raise RuntimeError(
            "No seeds resolved (CLI did not specify --seeds and plan.json "
            "had no boltz_seed values)."
        )
    LOG.info("Predicting with seeds: %s", seeds)

    # Resolve the receptor A3M for the predict stage.  Precedence:
    #   --msa-a3m PATH  →  msa_done.txt under --outdir  →  default path
    if args.stage == "predict":
        if args.msa_a3m is not None:
            receptor_msa_a3m = args.msa_a3m.resolve()
        elif msa_done_path.is_file():
            recorded = msa_done_path.read_text().strip().splitlines()
            if recorded:
                receptor_msa_a3m = Path(recorded[0]).resolve()
        if not receptor_msa_a3m.is_file():
            raise FileNotFoundError(
                f"Receptor A3M not found at {receptor_msa_a3m}.  "
                f"Run --stage msa first, or pass --msa-a3m PATH explicitly. "
                f"(msa_done.txt {'exists' if msa_done_path.is_file() else 'absent'} "
                f"at {msa_done_path}.)"
            )
        LOG.info("Using pre-computed receptor MSA: %s", receptor_msa_a3m)

    # Effector template CIF — extract once from --reference-pdb, reused
    # across seeds.  Mirrors what boltz2_negative_steering.py cmd_plan
    # does at cycle_0/effector_template.cif (boltz2_negative_steering.py:
    # extract_effector_template_cif).
    effector_template_cif = outdir / "inputs" / "effector_template.cif"
    if not effector_template_cif.is_file():
        extract_effector_template_cif(
            ground_truth_pdb=reference_pdb,
            truth_eff_chain=args.truth_eff_chain,
            pred_eff_chain=PRED_EFF_CHAIN,
            out_path=effector_template_cif,
        )
        LOG.info("Wrote effector template CIF: %s", effector_template_cif)
    else:
        LOG.info("Reusing existing effector template CIF: %s",
                 effector_template_cif)

    per_seed: List[Dict[str, object]] = []
    for seed in seeds:
        result = predict_one_seed(
            design_id=design_id,
            rec_seq=rec_seq,
            eff_seq=eff_seq,
            receptor_msa_a3m=receptor_msa_a3m,
            effector_template_cif=effector_template_cif,
            seed=seed,
            outdir=outdir,
            reference_pdb=reference_pdb,
            truth_rec_chain=args.truth_rec_chain,
            truth_eff_chain=args.truth_eff_chain,
            n_recycling=args.n_recycling,
        )
        per_seed.append(result)

    ok_seeds = [s for s in per_seed if s.get("status") == "ok"]

    def collect(metric: str) -> List[float]:
        out = []
        for s in ok_seeds:
            m = (s.get("metrics") or {}).get(metric)
            if m is not None:
                out.append(float(m))
        return out

    def collect_top(metric: str) -> List[float]:
        out = []
        for s in ok_seeds:
            v = s.get(metric)
            if v is not None:
                out.append(float(v))
        return out

    summary = {
        "design_id": design_id,
        "steered_design_name_in_plan": steered_design_name,
        "negsteer_workdir": str(workdir),
        "reference_pdb": str(reference_pdb),
        "steered_sequence": rec_seq,
        "effector_sequence": eff_seq,
        "n_seeds": len(per_seed),
        "n_seeds_ok": len(ok_seeds),
        "seeds": list(seeds),
        "complex_plddt_median":   _median(collect("complex_plddt")),
        "interface_plddt_median": _median(collect("interface_plddt")),
        "ipae_median":            _median(collect("ipae")),
        "pae_pass_frac_median":   _median(collect("pae_pass_frac")),
        "iptm_median":            _median(collect("iptm")),
        "ra_eff_median":          _median(collect_top("ra_eff")),
        "per_seed":               per_seed,
    }

    summary_path = outdir / f"{design_id}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))
    LOG.info("Wrote summary: %s", summary_path)
    LOG.info(
        "Done %s: %d/%d seeds OK, complex_plddt_median=%s, ra_eff_median=%s",
        design_id, len(ok_seeds), len(per_seed),
        summary["complex_plddt_median"], summary["ra_eff_median"],
    )

    if not ok_seeds:
        LOG.error("All seeds failed for %s", design_id)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
