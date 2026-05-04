"""Compare Boltz-2 with-MSA vs no-MSA confidence metrics for the 4a campaign.

Reads the per-design MSA summaries written by boltz2_msa_predict.py
and pairs them with the corresponding no-MSA medians from each
sequence's negsteer aggregated_results.csv.  Emits a comparison CSV
plus a stdout summary table that quantifies the H1 / H2 / H3 signals
the diagnostic was designed to test.

H1 — "no-MSA is the wrong control": both runs produce poor designs.
H2 — "MSA rescues some designs": MSA pLDDT crosses 0.7 from below.
H3 — "MSA improves rankings": deltas correlate with pipeline tier.

The no-MSA values come from the same medians the pipeline used for
its tier classification, so any per-design delta below is "what the
pipeline saw vs what the same design would have looked like with an
MSA".
"""

from __future__ import annotations

import experiments._path_setup  # noqa: F401  - adds bin/ to sys.path

import argparse
import csv
import json
import logging
import sys
from pathlib import Path
from statistics import mean
from typing import Dict, List, Optional, Tuple

LOG = logging.getLogger("compare_msa_vs_nomsa")


# ── Failure-mode thresholds — must match build_msa_manifest.py ─────────
PLDDT_THRESHOLD = 0.7
IPAE_THRESHOLD = 15.0
RA_EFF_THRESHOLD = 5.0


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Compare Boltz-2 with-MSA vs no-MSA confidence metrics for "
            "designs from a 4a negsteer run."
        )
    )
    p.add_argument("--msa-outdir", required=True, type=Path,
                   help="The --outdir from the SLURM array run "
                        "(boltz2_msa_predict.py wrote one "
                        "<design_id>/<design_id>_summary.json under here "
                        "per design).")
    p.add_argument("--negsteer-dir", required=True, type=Path,
                   help="The 4a negsteer output dir (same path you "
                        "passed to build_msa_manifest.py).")
    p.add_argument("--outfile", required=True, type=Path,
                   help="Where to write the comparison CSV.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose (DEBUG) logging.")
    return p.parse_args(argv)


# ───────────────────────────────────────────────────────────────────────
# Loaders
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


def load_msa_summary(summary_path: Path) -> Dict[str, Optional[float]]:
    """Read one boltz2_msa_predict.py per-design summary JSON."""
    data = json.loads(summary_path.read_text())
    return {
        "complex_plddt":  _try_float(data.get("complex_plddt_median")),
        "ipae":           _try_float(data.get("ipae_median")),
        "pae_pass_frac":  _try_float(data.get("pae_pass_frac_median")),
        "iptm":           _try_float(data.get("iptm_median")),
        "ra_eff":         _try_float(data.get("ra_eff_median")),
        "interface_plddt": _try_float(data.get("interface_plddt_median")),
    }


def is_passing(workdir: Path) -> bool:
    """passing_summary.csv has at least one data row → tier != none."""
    ps = workdir / "passing_summary.csv"
    if not ps.is_file():
        return False
    with ps.open() as f:
        rdr = csv.DictReader(f)
        return any(True for _ in rdr)


def representative_row(workdir: Path) -> Optional[Dict[str, str]]:
    """Same representative pick as build_msa_manifest.py: prefer
    rank_by_composite_score == 1, fall back to first row.
    """
    agg = workdir / "aggregated_results.csv"
    if not agg.is_file():
        return None
    with agg.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    for r in rows:
        if _try_int(r.get("rank_by_composite_score")) == 1:
            return r
    return rows[0]


def nomsa_metrics_from_row(row: Dict[str, str]) -> Dict[str, Optional[float]]:
    return {
        "complex_plddt":  _try_float(row.get("steered_complex_plddt_median")),
        "ipae":           _try_float(row.get("steered_ipae_median")),
        "pae_pass_frac":  _try_float(row.get("steered_pae_pass_frac_median")),
        "iptm":           _try_float(row.get("steered_iptm_median")),
        "ra_eff":         _try_float(row.get("steered_ra_eff_vs_truth_median")),
        "interface_plddt": _try_float(
            row.get("steered_interface_plddt_median")
        ),
    }


def dominant_failure_mode(metrics: Dict[str, Optional[float]]) -> str:
    """Same logic as build_msa_manifest.dominant_failure_mode but
    returning a string label suitable for the comparison CSV.
    Returns 'pass' when no metric is past threshold.
    """
    margins: Dict[str, float] = {}
    plddt = metrics.get("complex_plddt")
    if plddt is not None and plddt < PLDDT_THRESHOLD:
        margins["pLDDT"] = (PLDDT_THRESHOLD - plddt) / PLDDT_THRESHOLD
    ipae = metrics.get("ipae")
    if ipae is not None and ipae > IPAE_THRESHOLD:
        margins["ipAE"] = (ipae - IPAE_THRESHOLD) / IPAE_THRESHOLD
    ra_eff = metrics.get("ra_eff")
    if ra_eff is not None and ra_eff > RA_EFF_THRESHOLD:
        margins["ra_eff"] = (ra_eff - RA_EFF_THRESHOLD) / RA_EFF_THRESHOLD
    if not margins:
        return "pass"
    return max(margins, key=margins.__getitem__)


# ───────────────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────────────
COMPARE_FIELDS = [
    "design_id",
    "tier_nomsa",
    "dominant_failure_nomsa",
    "complex_plddt_nomsa",
    "ipae_nomsa",
    "pae_pass_frac_nomsa",
    "iptm_nomsa",
    "ra_eff_nomsa",
    "complex_plddt_msa",
    "ipae_msa",
    "pae_pass_frac_msa",
    "iptm_msa",
    "ra_eff_msa",
    "plddt_delta",
    "ipae_delta",
    "pae_pass_frac_delta",
    "iptm_delta",
    "ra_eff_delta",
]


def _delta(msa: Optional[float], nomsa: Optional[float],
           positive_when_msa_higher: bool) -> Optional[float]:
    if msa is None or nomsa is None:
        return None
    return (msa - nomsa) if positive_when_msa_higher else (nomsa - msa)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    msa_root = args.msa_outdir.resolve()
    negsteer_root = args.negsteer_dir.resolve()
    runs_dir = negsteer_root / "runs"

    # Discover MSA summaries: each design dir under msa-outdir contains
    # <design_id>/<design_id>_summary.json (the SLURM array writes
    # OUTDIR/<design_id>/ — see run_msa_array.slurm.sh).
    summaries: Dict[str, Path] = {}
    for design_dir in sorted(p for p in msa_root.iterdir() if p.is_dir()):
        candidate = design_dir / f"{design_dir.name}_summary.json"
        if candidate.is_file():
            summaries[design_dir.name] = candidate

    if not summaries:
        raise RuntimeError(
            f"No <design_id>/<design_id>_summary.json files found under "
            f"{msa_root}. Did the SLURM array run?"
        )
    LOG.info("Found %d MSA summaries under %s", len(summaries), msa_root)

    rows: List[Dict[str, object]] = []
    missing_workdirs: List[str] = []
    missing_repr: List[str] = []

    for design_id, summary_path in summaries.items():
        workdir = runs_dir / design_id
        if not workdir.is_dir():
            missing_workdirs.append(design_id)
            continue
        rep = representative_row(workdir)
        if rep is None:
            missing_repr.append(design_id)
            continue

        nomsa = nomsa_metrics_from_row(rep)
        msa = load_msa_summary(summary_path)
        passing = is_passing(workdir)
        tier_label = "passing" if passing else "none"
        dominant = dominant_failure_mode(nomsa)

        row = {
            "design_id": design_id,
            "tier_nomsa": tier_label,
            "dominant_failure_nomsa": dominant,
            "complex_plddt_nomsa": nomsa["complex_plddt"],
            "ipae_nomsa":          nomsa["ipae"],
            "pae_pass_frac_nomsa": nomsa["pae_pass_frac"],
            "iptm_nomsa":          nomsa["iptm"],
            "ra_eff_nomsa":        nomsa["ra_eff"],
            "complex_plddt_msa":   msa["complex_plddt"],
            "ipae_msa":            msa["ipae"],
            "pae_pass_frac_msa":   msa["pae_pass_frac"],
            "iptm_msa":            msa["iptm"],
            "ra_eff_msa":          msa["ra_eff"],
            "plddt_delta": _delta(
                msa["complex_plddt"], nomsa["complex_plddt"],
                positive_when_msa_higher=True,
            ),
            "ipae_delta": _delta(
                msa["ipae"], nomsa["ipae"],
                positive_when_msa_higher=False,
            ),
            "pae_pass_frac_delta": _delta(
                msa["pae_pass_frac"], nomsa["pae_pass_frac"],
                positive_when_msa_higher=True,
            ),
            "iptm_delta": _delta(
                msa["iptm"], nomsa["iptm"],
                positive_when_msa_higher=True,
            ),
            "ra_eff_delta": _delta(
                msa["ra_eff"], nomsa["ra_eff"],
                positive_when_msa_higher=False,
            ),
        }
        rows.append(row)

    if missing_workdirs:
        LOG.warning("Missing negsteer workdir for: %s",
                    ", ".join(missing_workdirs))
    if missing_repr:
        LOG.warning("No representative aggregated_results row for: %s",
                    ", ".join(missing_repr))

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    with args.outfile.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COMPARE_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    LOG.info("Wrote %d rows to %s", len(rows), args.outfile)

    _print_summary(rows)
    return 0


# ───────────────────────────────────────────────────────────────────────
# stdout summary
# ───────────────────────────────────────────────────────────────────────
def _mean_or_na(vals: List[Optional[float]]) -> str:
    finite = [v for v in vals if v is not None]
    if not finite:
        return "  n/a"
    return f"{mean(finite):+0.3f}"


def _print_summary(rows: List[Dict[str, object]]) -> None:
    if not rows:
        print("No rows in comparison; nothing to summarise.")
        return

    passing = [r for r in rows if r["tier_nomsa"] == "passing"]
    failing = [r for r in rows if r["tier_nomsa"] != "passing"]

    delta_cols = [
        ("plddt_delta",         "ΔpLDDT  (msa - nomsa, positive=better)"),
        ("ipae_delta",          "ΔipAE   (nomsa - msa, positive=better)"),
        ("pae_pass_frac_delta", "Δpae_pf (msa - nomsa, positive=better)"),
        ("iptm_delta",          "ΔiPTM   (msa - nomsa, positive=better)"),
        ("ra_eff_delta",        "Δra_eff (nomsa - msa, positive=better)"),
    ]

    print()
    print("=" * 70)
    print("Mean Δ per metric: passing vs failing")
    print("=" * 70)
    print(f"{'metric':<42s} {'passing':>10s} {'failing':>10s}")
    for key, label in delta_cols:
        p = _mean_or_na([r[key] for r in passing])
        f = _mean_or_na([r[key] for r in failing])
        print(f"{label:<42s} {p:>10s} {f:>10s}")

    print()
    print("=" * 70)
    print("Failing group: mean Δ per metric, by dominant failure mode")
    print("=" * 70)
    failure_modes = sorted({r["dominant_failure_nomsa"] for r in failing})
    if not failure_modes:
        print("  (no failing designs)")
    else:
        header_modes = "  ".join(f"{m:>10s}" for m in failure_modes)
        print(f"{'metric':<42s} {header_modes}")
        for key, label in delta_cols:
            cells = []
            for mode in failure_modes:
                vals = [r[key] for r in failing
                        if r["dominant_failure_nomsa"] == mode]
                cells.append(f"{_mean_or_na(vals):>10s}")
            print(f"{label:<42s} " + "  ".join(cells))
        for mode in failure_modes:
            n = sum(1 for r in failing
                    if r["dominant_failure_nomsa"] == mode)
            print(f"  n[{mode}] = {n}")

    n_h2 = sum(
        1 for r in rows
        if r["complex_plddt_nomsa"] is not None
        and r["complex_plddt_msa"] is not None
        and r["complex_plddt_nomsa"] < PLDDT_THRESHOLD
        and r["complex_plddt_msa"] >= PLDDT_THRESHOLD
    )
    n_h1 = sum(
        1 for r in rows
        if r["complex_plddt_nomsa"] is not None
        and r["complex_plddt_msa"] is not None
        and r["complex_plddt_nomsa"] < PLDDT_THRESHOLD
        and r["complex_plddt_msa"] < PLDDT_THRESHOLD
    )

    print()
    print("=" * 70)
    print("Hypothesis tallies (complex_plddt threshold = "
          f"{PLDDT_THRESHOLD})")
    print("=" * 70)
    print(f"  H2 (rescued by MSA: nomsa<{PLDDT_THRESHOLD} & msa>="
          f"{PLDDT_THRESHOLD}):   {n_h2:4d}")
    print(f"  H1 (uniformly poor: both < {PLDDT_THRESHOLD}):              "
          f"      {n_h1:4d}")
    print(f"  Total designs compared:                                "
          f"   {len(rows):4d}")
    print()


if __name__ == "__main__":
    sys.exit(main())
