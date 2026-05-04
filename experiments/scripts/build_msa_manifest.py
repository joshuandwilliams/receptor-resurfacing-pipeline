"""Build the SLURM-array manifest for the Boltz-2-with-MSA diagnostic.

Walks a 4a-style negsteer output tree and writes a plain-text file
listing one absolute per-sequence workdir path per array task.  Sampling
balances designs that passed the pipeline's own gating (tier A/B/C in
cross_sequence_summary terms) with designs that failed it, stratified
across the three dominant failure modes (low complex_plddt, high ipae,
high ra_eff vs truth).

Inputs assumed under --negsteer-dir:

    runs/<seq_name>/
      aggregated_results.csv
      passing_summary.csv

passing_summary.csv contains zero or more rows per workdir; an empty
passing_summary (header-only) means the sequence cleared no tier in the
pipeline (cross-sequence summary calls this "tier none") — i.e. failing.

aggregated_results.csv carries the per-sequence_group medians the
pipeline uses for ranking.  When picking a representative row for a
multi-sequence_group workdir, prefer rank_by_composite_score == 1
(matches the cross-sequence summary's representative pick); fall back
to the first row.

Failure-mode thresholds (from the user specification):

    pLDDT failure:  steered_complex_plddt_median   <  0.7    (lower = worse)
    ipAE failure:   steered_ipae_median            > 15      (higher = worse)
    ra_eff failure: steered_ra_eff_vs_truth_median >  5      (higher = worse)

For a design failing on more than one metric, the dominant mode is the
one whose value is furthest past threshold (margin in the "worse"
direction, normalised by the threshold so the modes are commensurable).
"""

from __future__ import annotations

import experiments._path_setup  # noqa: F401  - adds bin/ to sys.path

import argparse
import csv
import logging
import random
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

LOG = logging.getLogger("build_msa_manifest")


# ── Failure-mode thresholds ─────────────────────────────────────────────
PLDDT_THRESHOLD = 0.7    # complex_plddt < 0.7 is failing
IPAE_THRESHOLD = 15.0    # ipae > 15 is failing
RA_EFF_THRESHOLD = 5.0   # ra_eff_vs_truth > 5 is failing


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build the SLURM-array manifest for the Boltz-2-with-MSA "
            "diagnostic by sampling passing + failing designs from a "
            "4a negsteer output tree."
        ),
    )
    p.add_argument("--negsteer-dir", required=True, type=Path,
                   help="Path to the negsteer output dir from the 4a HPC "
                        "run.  Expects <negsteer-dir>/runs/<seq_name>/ "
                        "subdirectories with aggregated_results.csv and "
                        "passing_summary.csv.")
    p.add_argument("--outfile", required=True, type=Path,
                   help="Where to write the manifest (plain text).")
    p.add_argument("--n-passing", type=int, default=None,
                   help="Number of passing designs (tier != none) to "
                        "include.  Default: all.")
    p.add_argument("--n-failing", type=int, default=20,
                   help="Number of failing designs (tier none) to include, "
                        "stratified across the three dominant failure "
                        "modes (default: 20).")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for sampling (default: 42).")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose (DEBUG) logging.")
    return p.parse_args(argv)


# ───────────────────────────────────────────────────────────────────────
# Per-workdir loaders
# ───────────────────────────────────────────────────────────────────────
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


def _try_int(v) -> Optional[int]:
    f = _try_float(v)
    return None if f is None else int(f)


def is_passing(workdir: Path) -> bool:
    """A workdir is passing if its passing_summary.csv has at least one
    data row.  Empty (header-only) means tier-none.
    """
    ps = workdir / "passing_summary.csv"
    if not ps.is_file():
        return False
    with ps.open() as f:
        rdr = csv.DictReader(f)
        for _ in rdr:
            return True
    return False


def representative_row(workdir: Path) -> Optional[Dict[str, str]]:
    """Pick the representative aggregated_results.csv row for a workdir.

    Prefers rank_by_composite_score == 1 (the cross-sequence summary
    representative pick); falls back to the first row.  Returns None if
    aggregated_results.csv is missing or empty.
    """
    agg = workdir / "aggregated_results.csv"
    if not agg.is_file():
        return None
    rows: List[Dict[str, str]] = []
    with agg.open() as f:
        rdr = csv.DictReader(f)
        rows = list(rdr)
    if not rows:
        return None
    for r in rows:
        if _try_int(r.get("rank_by_composite_score")) == 1:
            return r
    return rows[0]


def metrics_from_row(row: Dict[str, str]) -> Dict[str, Optional[float]]:
    """Pull the three failure-mode-relevant medians from an
    aggregated_results.csv row.
    """
    return {
        "complex_plddt": _try_float(row.get("steered_complex_plddt_median")),
        "ipae":          _try_float(row.get("steered_ipae_median")),
        "ra_eff":        _try_float(row.get("steered_ra_eff_vs_truth_median")),
    }


def dominant_failure_mode(
    metrics: Dict[str, Optional[float]],
) -> Tuple[Optional[str], Dict[str, Optional[float]]]:
    """Return (mode, normalised_margins).

    Margin = how far the metric is past its threshold in the "worse"
    direction, divided by the threshold.  Negative or None margins mean
    the metric is not failing.  The dominant mode is the one with the
    largest positive normalised margin; if no metric is past threshold,
    the mode is None (meaning: by these three metrics this design is
    not actually failing — surface it for review).
    """
    margins: Dict[str, Optional[float]] = {}

    plddt = metrics.get("complex_plddt")
    margins["plddt"] = (
        (PLDDT_THRESHOLD - plddt) / PLDDT_THRESHOLD if plddt is not None else None
    )

    ipae = metrics.get("ipae")
    margins["ipae"] = (
        (ipae - IPAE_THRESHOLD) / IPAE_THRESHOLD if ipae is not None else None
    )

    ra_eff = metrics.get("ra_eff")
    margins["ra_eff"] = (
        (ra_eff - RA_EFF_THRESHOLD) / RA_EFF_THRESHOLD
        if ra_eff is not None else None
    )

    positive = {k: v for k, v in margins.items() if v is not None and v > 0}
    if not positive:
        return None, margins
    mode = max(positive, key=positive.__getitem__)
    return mode, margins


# ───────────────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────────────
def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    runs_dir = args.negsteer_dir.resolve() / "runs"
    if not runs_dir.is_dir():
        raise FileNotFoundError(
            f"Expected runs/ subdirectory under {args.negsteer_dir}: "
            f"{runs_dir} not found."
        )

    workdirs = sorted(p for p in runs_dir.iterdir() if p.is_dir())
    LOG.info("Found %d candidate per-sequence workdirs under %s",
             len(workdirs), runs_dir)

    passing: List[Path] = []
    failing_by_mode: Dict[str, List[Path]] = {
        "plddt": [],
        "ipae": [],
        "ra_eff": [],
    }
    failing_unclassified: List[Path] = []
    skipped: List[Tuple[Path, str]] = []

    for wd in workdirs:
        # Skip the controls (input_control_polyA, input_control_scrambled)
        # and any helper subdir that is not a real sequence workdir.
        if wd.name.startswith("input_control"):
            skipped.append((wd, "control"))
            continue
        if not (wd / "aggregated_results.csv").is_file():
            skipped.append((wd, "no aggregated_results.csv"))
            continue

        if is_passing(wd):
            passing.append(wd)
            continue

        rep = representative_row(wd)
        if rep is None:
            skipped.append((wd, "no representative row"))
            continue
        metrics = metrics_from_row(rep)
        mode, _margins = dominant_failure_mode(metrics)
        if mode is None:
            failing_unclassified.append(wd)
        else:
            failing_by_mode[mode].append(wd)

    LOG.info("Passing (tier != none):     %d", len(passing))
    for mode, lst in failing_by_mode.items():
        LOG.info("Failing dominant=%s:%s%d",
                 mode, " " * (8 - len(mode)), len(lst))
    LOG.info("Failing (no metric past threshold): %d",
             len(failing_unclassified))
    LOG.info("Skipped:                    %d", len(skipped))

    rng = random.Random(args.seed)

    n_passing_target = (
        len(passing) if args.n_passing is None
        else min(args.n_passing, len(passing))
    )
    chosen_passing = (
        list(passing) if n_passing_target == len(passing)
        else rng.sample(passing, n_passing_target)
    )

    # Stratify failing across the three modes evenly, with overflow
    # spread by remainder.  If a mode is short of its quota, the
    # leftover budget cascades to the next mode in deterministic order
    # so the manifest is reproducible.
    target_total_failing = max(0, args.n_failing)
    mode_order = ["plddt", "ipae", "ra_eff"]
    base_quota = target_total_failing // len(mode_order)
    remainder = target_total_failing % len(mode_order)
    quotas = {m: base_quota + (1 if i < remainder else 0)
              for i, m in enumerate(mode_order)}

    chosen_failing: List[Path] = []
    leftover_budget = 0
    pool_remaining: Dict[str, List[Path]] = {}
    for mode in mode_order:
        pool = list(failing_by_mode[mode])
        rng.shuffle(pool)
        quota = quotas[mode]
        take = min(quota, len(pool))
        chosen_failing.extend(pool[:take])
        leftover_budget += quota - take
        pool_remaining[mode] = pool[take:]

    # Spend leftover budget across whatever pools still have entries,
    # in mode_order, then on failing_unclassified as a last resort.
    if leftover_budget > 0:
        spillover = []
        for mode in mode_order:
            spillover.extend(pool_remaining[mode])
        rng.shuffle(failing_unclassified)
        spillover.extend(failing_unclassified)
        for wd in spillover:
            if leftover_budget == 0:
                break
            if wd in chosen_failing:
                continue
            chosen_failing.append(wd)
            leftover_budget -= 1

    chosen_passing.sort(key=lambda p: p.name)
    chosen_failing.sort(key=lambda p: p.name)

    final_breakdown: Dict[str, int] = {m: 0 for m in mode_order}
    final_breakdown["unclassified"] = 0
    for wd in chosen_failing:
        rep = representative_row(wd)
        mode = (
            dominant_failure_mode(metrics_from_row(rep))[0]
            if rep is not None else None
        )
        final_breakdown[mode if mode in final_breakdown else "unclassified"] += 1

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with args.outfile.open("w") as f:
        f.write(f"# Boltz-2-with-MSA diagnostic manifest\n")
        f.write(f"# Generated:           {now}\n")
        f.write(f"# negsteer-dir:        {args.negsteer_dir.resolve()}\n")
        f.write(f"# Total designs:       {len(chosen_passing)} passing"
                f" + {len(chosen_failing)} failing"
                f" = {len(chosen_passing) + len(chosen_failing)}\n")
        f.write(f"# Passing pool size:   {len(passing)}"
                f" (took {len(chosen_passing)})\n")
        f.write(f"# Failing pool sizes:  "
                f"plddt={len(failing_by_mode['plddt'])} "
                f"ipae={len(failing_by_mode['ipae'])} "
                f"ra_eff={len(failing_by_mode['ra_eff'])} "
                f"unclassified={len(failing_unclassified)}\n")
        f.write(f"# Failing breakdown:   "
                f"plddt={final_breakdown['plddt']} "
                f"ipae={final_breakdown['ipae']} "
                f"ra_eff={final_breakdown['ra_eff']} "
                f"unclassified={final_breakdown['unclassified']}\n")
        f.write(f"# Sampling seed:       {args.seed}\n")
        f.write(f"# Failure thresholds:  "
                f"plddt<{PLDDT_THRESHOLD} ipae>{IPAE_THRESHOLD} "
                f"ra_eff>{RA_EFF_THRESHOLD}\n")
        f.write(f"#\n")
        f.write(f"# --- passing designs ({len(chosen_passing)}) ---\n")
        for wd in chosen_passing:
            f.write(f"{wd.resolve()}\n")
        f.write(f"# --- failing designs ({len(chosen_failing)}) ---\n")
        for wd in chosen_failing:
            f.write(f"{wd.resolve()}\n")

    LOG.info(
        "Wrote manifest with %d entries to %s "
        "(passing=%d, failing=%d; failing breakdown plddt=%d ipae=%d "
        "ra_eff=%d unclassified=%d)",
        len(chosen_passing) + len(chosen_failing),
        args.outfile,
        len(chosen_passing), len(chosen_failing),
        final_breakdown["plddt"], final_breakdown["ipae"],
        final_breakdown["ra_eff"], final_breakdown["unclassified"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
