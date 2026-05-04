"""Build the SLURM-array manifest for the Boltz-2-with-MSA diagnostic.

Walks a 4a-style negsteer output tree and writes a plain-text file
listing one absolute per-sequence workdir path per array task.  Sampling
balances designs that passed the pipeline's own gating (tier A/B/C in
cross_sequence_summary terms) with designs that failed it.

Failing designs are picked by greedy set cover across failure modes:
each design has a SET of modes it fails on (a typical failing design
fails on multiple), and we iteratively pick the design that covers the
most still-under-quota modes until every mode has at least
--n-per-mode covered or the failing pool is exhausted.  This avoids
the under-coverage produced by single-dominant-mode stratification.

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

Failure-mode thresholds:

    pLDDT failure:        steered_complex_plddt_median   <  0.7
    ipAE failure:         steered_ipae_median            > 15
    pae_pass_frac fail:   steered_pae_pass_frac_median   <  0.1
    iPTM failure:         steered_iptm_median            <  0.3
    ra_eff failure:       steered_ra_eff_vs_truth_median >  5
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
from typing import Dict, List, Optional, Set, Tuple

LOG = logging.getLogger("build_msa_manifest")


# ── Failure-mode thresholds ─────────────────────────────────────────────
PLDDT_THRESHOLD = 0.7           # complex_plddt < 0.7 is failing
IPAE_THRESHOLD = 15.0           # ipae > 15 is failing
PAE_PASS_FRAC_THRESHOLD = 0.1   # pae_pass_frac < 0.1 is failing
IPTM_THRESHOLD = 0.3            # iptm < 0.3 is failing
RA_EFF_THRESHOLD = 5.0          # ra_eff_vs_truth > 5 is failing

# Canonical mode order — drives column layout in logs and the manifest
# header.  The names match the metric stems used elsewhere in the repo.
FAILURE_MODES: Tuple[str, ...] = (
    "plddt", "ipae", "pae_pass_frac", "iptm", "ra_eff",
)


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
    p.add_argument("--n-per-mode", type=int, default=5,
                   help="Per-failure-mode coverage target for the failing "
                        "sample.  Greedy set-cover keeps picking failing "
                        "designs until every failure mode has at least "
                        "this many examples (or the failing pool is "
                        "exhausted).  Default: 5.")
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
    """Pull the failure-mode-relevant medians from an
    aggregated_results.csv row.
    """
    return {
        "complex_plddt": _try_float(row.get("steered_complex_plddt_median")),
        "ipae":          _try_float(row.get("steered_ipae_median")),
        "pae_pass_frac": _try_float(row.get("steered_pae_pass_frac_median")),
        "iptm":          _try_float(row.get("steered_iptm_median")),
        "ra_eff":        _try_float(row.get("steered_ra_eff_vs_truth_median")),
    }


def failure_modes_for(metrics: Dict[str, Optional[float]]) -> Set[str]:
    """Return the SET of failure-mode labels this design is failing on.

    A metric that did not parse (None) does NOT count as failing — it
    contributes nothing to coverage and is silently absent from the set.
    """
    modes: Set[str] = set()
    plddt = metrics.get("complex_plddt")
    if plddt is not None and plddt < PLDDT_THRESHOLD:
        modes.add("plddt")
    ipae = metrics.get("ipae")
    if ipae is not None and ipae > IPAE_THRESHOLD:
        modes.add("ipae")
    pae_pass_frac = metrics.get("pae_pass_frac")
    if pae_pass_frac is not None and pae_pass_frac < PAE_PASS_FRAC_THRESHOLD:
        modes.add("pae_pass_frac")
    iptm = metrics.get("iptm")
    if iptm is not None and iptm < IPTM_THRESHOLD:
        modes.add("iptm")
    ra_eff = metrics.get("ra_eff")
    if ra_eff is not None and ra_eff > RA_EFF_THRESHOLD:
        modes.add("ra_eff")
    return modes


def greedy_set_cover(
    candidates: List[Tuple[Path, Set[str]]],
    n_per_mode: int,
    rng: random.Random,
) -> List[Path]:
    """Greedy set-cover sampling: at each step pick the design whose
    failure-mode set covers the largest number of modes that are still
    under quota.  Modes already at >= n_per_mode are dropped from the
    "needed" set when scoring candidates, so a design that fails on a
    saturated mode but also on an under-quota one still gets chosen for
    the under-quota mode it adds.

    Ties are broken at random using rng (so the output is reproducible
    with --seed).  Stops when every mode is covered to n_per_mode OR the
    candidate pool is exhausted.

    Designs that fail on zero modes (their representative row had every
    metric finite and on the right side of every threshold — surfaces
    only when the upstream tier classification disagrees with our
    threshold panel) are skipped: they cannot contribute to coverage.
    """
    coverage: Dict[str, int] = {m: 0 for m in FAILURE_MODES}
    chosen: List[Path] = []
    pool: List[Tuple[Path, Set[str]]] = [
        (wd, modes) for wd, modes in candidates if modes
    ]

    while pool:
        needed = {m for m, c in coverage.items() if c < n_per_mode}
        if not needed:
            break

        best_score = -1
        best_bucket: List[Tuple[Path, Set[str]]] = []
        for wd, modes in pool:
            score = len(modes & needed)
            if score == 0:
                continue
            if score > best_score:
                best_score = score
                best_bucket = [(wd, modes)]
            elif score == best_score:
                best_bucket.append((wd, modes))

        if best_score <= 0:
            # No remaining candidate covers any still-needed mode →
            # the under-quota modes have no examples in the pool.
            break

        pick_wd, pick_modes = rng.choice(best_bucket)
        chosen.append(pick_wd)
        for m in pick_modes:
            coverage[m] += 1
        pool = [(wd, modes) for wd, modes in pool if wd != pick_wd]

    return chosen


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
    failing_with_modes: List[Tuple[Path, Set[str]]] = []
    failing_no_modes: List[Path] = []
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
        modes = failure_modes_for(metrics_from_row(rep))
        if modes:
            failing_with_modes.append((wd, modes))
        else:
            failing_no_modes.append(wd)

    pool_per_mode: Dict[str, int] = {m: 0 for m in FAILURE_MODES}
    for _wd, modes in failing_with_modes:
        for m in modes:
            pool_per_mode[m] += 1

    LOG.info("Passing (tier != none):                    %d", len(passing))
    LOG.info("Failing with at least one mode flagged:    %d",
             len(failing_with_modes))
    LOG.info("Failing pool coverage per mode:")
    for m in FAILURE_MODES:
        LOG.info("  %-14s %d", m, pool_per_mode[m])
    LOG.info("Failing with no metric past threshold:     %d",
             len(failing_no_modes))
    LOG.info("Skipped:                                   %d", len(skipped))

    rng = random.Random(args.seed)

    n_passing_target = (
        len(passing) if args.n_passing is None
        else min(args.n_passing, len(passing))
    )
    chosen_passing = (
        list(passing) if n_passing_target == len(passing)
        else rng.sample(passing, n_passing_target)
    )

    # Shuffle the failing candidate list before set-cover so any
    # deterministic ordering the filesystem produced does not bias
    # tie-breaks before rng.choice() ever runs.
    failing_shuffled = list(failing_with_modes)
    rng.shuffle(failing_shuffled)
    chosen_failing = greedy_set_cover(
        failing_shuffled, args.n_per_mode, rng,
    )

    chosen_passing.sort(key=lambda p: p.name)
    chosen_failing.sort(key=lambda p: p.name)

    # Coverage achieved by the chosen failing set — count how many
    # selected designs fail on each mode.  A single design can
    # contribute to multiple modes, so the column sum can exceed
    # len(chosen_failing).
    chosen_modes = {wd: modes for wd, modes in failing_with_modes}
    final_coverage: Dict[str, int] = {m: 0 for m in FAILURE_MODES}
    for wd in chosen_failing:
        for m in chosen_modes.get(wd, set()):
            final_coverage[m] += 1

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    threshold_summary = (
        f"plddt<{PLDDT_THRESHOLD} "
        f"ipae>{IPAE_THRESHOLD} "
        f"pae_pass_frac<{PAE_PASS_FRAC_THRESHOLD} "
        f"iptm<{IPTM_THRESHOLD} "
        f"ra_eff>{RA_EFF_THRESHOLD}"
    )
    pool_coverage_summary = " ".join(
        f"{m}={pool_per_mode[m]}" for m in FAILURE_MODES
    )
    selected_coverage_summary = " ".join(
        f"{m}={final_coverage[m]}" for m in FAILURE_MODES
    )
    with args.outfile.open("w") as f:
        f.write("# Boltz-2-with-MSA diagnostic manifest\n")
        f.write(f"# Generated:               {now}\n")
        f.write(f"# negsteer-dir:            {args.negsteer_dir.resolve()}\n")
        f.write(f"# Total designs selected:  {len(chosen_passing)} passing"
                f" + {len(chosen_failing)} failing"
                f" = {len(chosen_passing) + len(chosen_failing)}\n")
        f.write(f"# Passing pool:            {len(passing)} "
                f"(took {len(chosen_passing)})\n")
        f.write(f"# Failing pool coverage:   {pool_coverage_summary}\n")
        f.write(f"# Selected failing cover:  {selected_coverage_summary}\n")
        f.write(f"# Per-mode target:         --n-per-mode {args.n_per_mode}\n")
        f.write(f"# Sampling seed:           {args.seed}\n")
        f.write(f"# Failure thresholds:      {threshold_summary}\n")
        f.write("#\n")
        f.write(f"# --- passing designs ({len(chosen_passing)}) ---\n")
        for wd in chosen_passing:
            f.write(f"{wd.resolve()}\n")
        f.write(f"# --- failing designs ({len(chosen_failing)}) ---\n")
        for wd in chosen_failing:
            f.write(f"{wd.resolve()}\n")

    LOG.info(
        "Wrote manifest with %d entries to %s "
        "(passing=%d, failing=%d; selected failing coverage %s; "
        "per-mode target=%d)",
        len(chosen_passing) + len(chosen_failing),
        args.outfile,
        len(chosen_passing), len(chosen_failing),
        selected_coverage_summary, args.n_per_mode,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
