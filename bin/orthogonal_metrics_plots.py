#!/usr/bin/env python3
"""
orthogonal_metrics_plots.py
---------------------------
Production cohort plot script for orthogonal metrics (Task 44).

Reads survivors_with_orthogonal_metrics.csv produced by
NEGSTEER_ORTHOGONAL_METRICS, joined to cross_sequence_summary.csv
from the upstream NEGSTEER_CROSS_SEQUENCE.  Both are required: the
metrics-vs-composite scatter and the combined cohort+orthogonal
summary both join on mpnn_sequence to source upstream composite
scores and cohort-wide context.  Invoked from
modules/negsteer_orthogonal_metrics_plots.nf::ORTHOG_PLOTS.

Lifted verbatim from tests/orthogonal_metrics/test_orthogonal_metrics_plots.py;
the only differences are (a) the test-path fallback in
_resolve_csv_path is removed (production always passes both CSVs
explicitly), (b) both CSV flags are required, and (c) this header.
The two scripts must stay in sync; when iterating on plots, edit
the test script first, validate against real cluster data, then
mirror the changes here.

Plots produced
--------------

orthogonal_af3_vs_boltz.png
    Scatter of af3_nomsa_best_ra_eff (x) vs Boltz
    representative_ra_eff_vs_truth_median (y).  AF3 region failing
    the af3_ra_max threshold shaded pale red — but AF3 is
    INFORMATIONAL only (no longer gates passes_orthogonal_filters).

orthogonal_filter_cascade.png
    Waterfall of the GATING cascade only (sc → bsa → interface_plddt).
    AF3 disagreement is shown as a separate informational annotation
    in the corner — disagreement is reported but doesn't drop
    survivors (negsteer optimises Boltz; AF3 is a sanity check).

orthogonal_metrics_vs_composite.png
    Multi-panel scatter: composite_score (Y) vs each orthogonal
    metric (X).  Six panels: Sc, BSA, interface_plddt, ΔΔG,
    AF3 ra_eff, AF3 iPTM.  Replaces the earlier strip-plot
    distributions and standalone ΔΔG plot — both carried no
    information at small n.  Composite is the upstream ranker, so
    showing every orthogonal metric against it shows whether each
    orthogonal metric agrees with or fights the composite ranking.

orthogonal_combined_cohort_summary.png
    The big combined view: one row per MPNN sequence in cross_summary
    (steered only), sorted by cross_composite_score descending.
    Columns span both the negsteer cascade (Boltz pose, confidence
    suite, n_mutations) and the orthogonal cascade (AF3, Sc, BSA,
    ΔΔG).  Sequences that didn't reach the orthogonal stage have
    grey-hatched cells in those columns.  Tier rendered as left-edge
    coloured stripe.  Designed for wet-lab triage of the full cohort.

Usage
-----
    Invoked by modules/negsteer_orthogonal_metrics_plots.nf::ORTHOG_PLOTS.
    Direct invocation:

        python orthogonal_metrics_plots.py \\
            --survivors-csv     <path to survivors_with_orthogonal_metrics.csv> \\
            --cross-summary-csv <path to cross_sequence_summary.csv> \\
            --outdir            <output directory>

Both --survivors-csv and --cross-summary-csv are required.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator


# ── Threshold defaults (must match params_example.yml) ───────────────────
ORTHOG_AF3_RA_MAX  = 5.0       # mandatory hard drop
ORTHOG_SC_MIN      = 0.55      # Lawrence-Colman shape complementarity
ORTHOG_BSA_MIN     = 600.0     # Å² (Overath et al. 2025)
ORTHOG_PLDDT_MIN   = 0.75      # mean interface pLDDT

# Bennett 2023 reference for ΔΔG_binding on de novo binders.
DDG_BENNETT_REFERENCE = -30.0

COMPOSITE_RA_EFF_WEIGHT = 0.05    # negsteer composite for cross-context


# ── Styling ──────────────────────────────────────────────────────────────
COLOUR_PASS         = "#2CA02C"   # green
COLOUR_FAIL         = "#D62728"   # red
COLOUR_THRESHOLD    = "#D62728"   # red dashed (threshold reference)
COLOUR_SURVIVOR     = "#4C72B0"   # blue
COLOUR_DIAG         = "#888888"   # grey (y=x reference)
COLOUR_BENNETT_REF  = "#1F77B4"   # blue dashed (ΔΔG reference)

COLOUR_TIER = {
    "A":    "#2CA02C",
    "B":    "#FFB000",
    "C":    "#FF7F0E",
    "none": "#B0B0B0",
}

# Sort order for tier-aware row sorts.  Lower number = appears first.
# Mirrors _TIER_ORDER in cross_sequence_summary.py so the cohort
# summary plot's row order matches what `cross_rank_by_composite`
# produces in the CSV.
_TIER_SORT_ORDER = {"A": 0, "B": 1, "C": 2, "none": 3}


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _try_float(v) -> Optional[float]:
    if v is None:
        return None
    if isinstance(v, (int, float)):
        f = float(v)
    else:
        s = str(v).strip()
        if not s:
            return None
        try:
            f = float(s)
        except ValueError:
            return None
    if f != f or f in (float("inf"), float("-inf")):
        return None
    return f


def _try_int(v) -> Optional[int]:
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    try:
        return int(s)
    except ValueError:
        try:
            return int(float(s))
        except ValueError:
            return None


def _parse_flags(s: str) -> List[str]:
    """Parse the comma-separated orthogonal_flags column.  Each flag
    can have a colon-suffix (e.g. 'sc_too_low:0.512' → 'sc_too_low')."""
    if not s:
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        # Strip the value suffix to get just the flag name.
        out.append(tok.split(":", 1)[0])
    return out


def _passes(row: Dict) -> bool:
    """True if passes_orthogonal_filters == '1'."""
    return (row.get("passes_orthogonal_filters") or "").strip() == "1"


def _make_empty_plot(message: str, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center",
            transform=ax.transAxes, fontsize=12, color="grey")
    ax.set_axis_off()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _short_name(name: str) -> str:
    """Compact label: 'design_0_seq_1' → 'd0_s1'."""
    n = name
    if n.startswith("design_"):
        n = n[len("design_"):]
        n = n.replace("_seq_", "_s")
        return f"d{n}"
    return n


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1: AF3 vs Boltz ra_eff scatter (the headline orthogonal plot)
# ═══════════════════════════════════════════════════════════════════════════

def plot_af3_vs_boltz(rows: List[Dict], out_path: str) -> bool:
    """Scatter: AF3-no-MSA best ra_eff (x) vs Boltz median ra_eff (y).

    Points are coloured by negsteer tier (A/B/C/none) — see COLOUR_TIER.
    The threshold line + shaded fail-region convey orthogonal pass/fail
    visually for AF3 ra_eff; readers can see at a glance which tier-A
    candidates pass the orthogonal cross-check vs which fail it (e.g.
    a tier-A point on the right side of the threshold = a candidate
    that passed negsteer's gating but disagreed with AF3).

    Pre-2026-04-29 history: this plot used to bucket points by
    passes_orthogonal_filters (global pass-all). That made every
    point appear "Fail" if it failed any one filter, even in panels
    where its own metric was fine — confusing and misleading.  Tier-
    colouring sidesteps the issue entirely: tier is a discrete
    cohort-level signal, threshold lines speak for per-metric pass/
    fail.
    """
    # Bucket points by tier.
    pts_by_tier: Dict[str, List[Tuple[float, float, str]]] = {
        "A": [], "B": [], "C": [], "none": [],
    }
    pts_missing_af3: List[Tuple[float, str]] = []   # only Boltz available

    for row in rows:
        boltz_ra = _try_float(row.get("representative_ra_eff_vs_truth_median"))
        af3_ra   = _try_float(row.get("af3_nomsa_best_ra_eff"))
        if boltz_ra is None and af3_ra is None:
            continue
        name = row.get("mpnn_sequence", "?")
        if af3_ra is None and boltz_ra is not None:
            pts_missing_af3.append((boltz_ra, name))
            continue
        if boltz_ra is None:
            continue
        tier = (row.get("cross_tier") or "none").strip() or "none"
        if tier not in pts_by_tier:
            tier = "none"
        pts_by_tier[tier].append((af3_ra, boltz_ra, name))

    total_plotted = sum(len(v) for v in pts_by_tier.values())
    if not (total_plotted or pts_missing_af3):
        _make_empty_plot("No ra_eff data available", out_path)
        return False

    fig, ax = plt.subplots(figsize=(7, 6))

    # Compute axis range.  AF3 routinely returns 30–40 Å on bad
    # candidates so the x-axis must accommodate that, not just the
    # threshold * 1.4.
    all_x = [p[0] for v in pts_by_tier.values() for p in v]
    all_y = ([p[1] for v in pts_by_tier.values() for p in v]
             + [p[0] for p in pts_missing_af3])
    if not all_x:
        x_max = ORTHOG_AF3_RA_MAX * 2
    else:
        x_max = max(all_x) * 1.10
    if not all_y:
        y_max = ORTHOG_AF3_RA_MAX * 2
    else:
        y_max = max(all_y) * 1.10
    axis_max = max(x_max, y_max, ORTHOG_AF3_RA_MAX * 1.5)

    # Pale-red shading for the AF3-fail region (x ≥ ORTHOG_AF3_RA_MAX).
    # This is the orthogonal cascade's mandatory hard drop.
    ax.axvspan(ORTHOG_AF3_RA_MAX, axis_max,
               color=COLOUR_FAIL, alpha=0.08, zorder=0)

    # y = x diagonal — perfect cross-model agreement reference.
    ax.plot([0, axis_max], [0, axis_max], color=COLOUR_DIAG,
            linewidth=0.8, linestyle="--", alpha=0.6, zorder=1,
            label="y = x (perfect agreement)")

    # AF3 threshold as a thin reference line (the shading is the
    # primary visual; the line is just the boundary).
    ax.axvline(ORTHOG_AF3_RA_MAX, color=COLOUR_THRESHOLD,
               linestyle="--", linewidth=1.0, alpha=0.7, zorder=1,
               label=f"af3_ra_max = {ORTHOG_AF3_RA_MAX} Å")

    # Plot one scatter per tier with the tier colour.  Order matters
    # for z-stacking — A on top of B on top of C on top of none, so
    # higher-quality tiers aren't hidden by lower-quality ones in
    # clustered regions.
    for tier in ("none", "C", "B", "A"):
        tier_pts = pts_by_tier[tier]
        if not tier_pts:
            continue
        ax.scatter(
            [p[0] for p in tier_pts],
            [p[1] for p in tier_pts],
            s=80, c=COLOUR_TIER[tier], edgecolor="black",
            linewidth=0.5, alpha=0.9,
            label=f"Tier {tier} (n={len(tier_pts)})"
                  if tier != "none" else f"No tier (n={len(tier_pts)})",
            zorder={"none": 3, "C": 4, "B": 5, "A": 6}[tier],
        )

    # Missing AF3 — drawn at the right edge as an open right-arrow.
    if pts_missing_af3:
        for boltz_ra, name in pts_missing_af3:
            ax.scatter(axis_max * 0.97, boltz_ra,
                       marker=">", s=100, c="white",
                       edgecolor="black", linewidth=0.6, alpha=0.85,
                       zorder=4)
        ax.text(axis_max * 0.97, axis_max * 0.02,
                f"AF3 missing\n(n={len(pts_missing_af3)})",
                ha="right", va="bottom", fontsize=8, color="black",
                style="italic")

    # (Survivor-name annotations removed 2026-04-28.  Even with greedy
    # 4-corner anti-overlap placement, labels in clustered regions read
    # as visual noise and don't add information beyond the pass/fail
    # colour.  The combined cohort summary plot carries the per-survivor
    # detail; this scatter is for cohort-level cross-model agreement.)

    ax.set_xlim(0, axis_max)
    ax.set_ylim(0, axis_max)
    ax.set_xlabel("AF3-no-MSA best ra_eff (Å)", fontsize=10)
    ax.set_ylabel("Boltz median ra_eff (Å)", fontsize=10)
    ax.grid(True, linestyle="-", linewidth=0.3, alpha=0.4, color="grey")
    ax.set_axisbelow(True)
    # Legend above the axes so it never overlaps data points.  ncol
    # spreads the entries horizontally so they fit in one row above
    # the plot.
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=4, fontsize=8, framealpha=0.95)
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2: Filter cascade
# ═══════════════════════════════════════════════════════════════════════════

# Order matters: the cascade applies filters in this order so each
# bar represents "would survive if only this and earlier filters were
# applied".  AF3 was demoted from gating to informational on
# 2026-04-28; it appears as a TRAILING bar after the gating cascade,
# rendered in yellow rather than the cascade-red, so the
# informational-vs-gating distinction is visually obvious.
FILTER_ORDER = [
    ("sc_too_low",                 f"sc\n≥ {ORTHOG_SC_MIN}",
     "sc_missing"),
    ("bsa_too_low",                f"+ bsa\n≥ {ORTHOG_BSA_MIN:.0f}",
     "bsa_missing"),
    ("interface_plddt_too_low",    f"+ interface_plddt\n≥ {ORTHOG_PLDDT_MIN}",
     "interface_plddt_missing"),
]

# Informational AF3 flags (trailing bar, not in the gating cascade).
AF3_FLAGS = ("af3_nomsa_ra_eff_too_high", "af3_nomsa_missing")

# Yellow for the AF3-trailing-bar drop.  Distinguishes "lost to AF3
# disagreement" (informational) from cascade-red drops (gating).
COLOUR_AF3_DROP = "#FFB000"


# Columns whose presence indicates a survivor actually entered the
# orthogonal stage (i.e. AF3 / biophysical / Rosetta were attempted).
# `survivors_with_orthogonal_metrics.csv` is produced by passing the
# upstream extended CSV through merge_orthogonal_metrics.py, which
# carries forward EVERY upstream row even if no orthogonal metrics
# were computed for it (those rows get blank orthogonal columns +
# `*_missing` flags).  The cascade should only reflect what was
# actually filtered, so we pre-filter on any-of-these-non-empty.
ORTHOG_PRESENCE_COLS = (
    "af3_nomsa_total_predictions",
    "bsa",
    "sc",
    "rosetta_ddg",
    "representative_interface_plddt_median",
)


def _entered_orthogonal_stage(row: Dict) -> bool:
    """True if any orthogonal-stage column has a non-empty value."""
    for col in ORTHOG_PRESENCE_COLS:
        v = row.get(col, "")
        if v is not None and str(v).strip() != "":
            return True
    return False


def plot_filter_cascade(rows: List[Dict], out_path: str) -> bool:
    """Waterfall: how many survivors pass each gating filter in
    succession.  Reads orthogonal_flags directly to determine which
    filter trips each survivor.

    Pre-filters to rows that actually entered the orthogonal stage —
    the merge step copies through every upstream row regardless of
    whether AF3/biophysical/Rosetta ran, so the raw row count
    overstates what the cascade is actually filtering.

    Layout: standard cumulative-AND cascade for the gating filters
    (sc → bsa → interface_plddt), then a trailing INFORMATIONAL bar
    for AF3 disagreement coloured yellow rather than the cascade-red,
    to make the gating/informational distinction visually obvious."""
    if not rows:
        _make_empty_plot("No survivors in input CSV", out_path)
        return False

    rows = [r for r in rows if _entered_orthogonal_stage(r)]
    if not rows:
        _make_empty_plot(
            "No rows entered the orthogonal stage", out_path,
        )
        return False

    n_total = len(rows)

    # Per-survivor: the set of failing-flag names (any prefix match).
    flag_sets: List[set] = []
    for row in rows:
        flags = set()
        for tok in _parse_flags(row.get("orthogonal_flags", "")):
            # Strip any ":<value>" suffix on the flag name so the
            # in-FILTER_ORDER lookups match.
            flags.add(tok.split(":", 1)[0])
        flag_sets.append(flags)

    # Cumulative-AND survival counts across the GATING cascade.
    survivors = [n_total]
    drops = [0]
    cumulative_excluded: set = set()
    for fail_flag, _label, missing_flag in FILTER_ORDER:
        newly_dropped = 0
        for i, flags in enumerate(flag_sets):
            if i in cumulative_excluded:
                continue
            if fail_flag in flags or missing_flag in flags:
                cumulative_excluded.add(i)
                newly_dropped += 1
        survivors.append(survivors[-1] - newly_dropped)
        drops.append(newly_dropped)

    # AF3 trailing bar: of the rows that survived the gating cascade,
    # how many would also have been dropped by AF3 (informational).
    af3_drop_from_survivors = 0
    for i, flags in enumerate(flag_sets):
        if i in cumulative_excluded:
            continue
        if any(f in flags for f in AF3_FLAGS):
            af3_drop_from_survivors += 1
    af3_remaining = survivors[-1] - af3_drop_from_survivors

    # Append the AF3 step to the bar arrays.
    survivors.append(af3_remaining)
    drops.append(af3_drop_from_survivors)

    labels = (
        ["All survivors"]
        + [f[1] for f in FILTER_ORDER]
        + ["+ AF3 ra_eff\n(informational)"]
    )

    fig, ax = plt.subplots(figsize=(9, 5))
    xs = np.arange(len(labels))

    # Pass bars (always green).
    ax.bar(xs, survivors, color=COLOUR_PASS, edgecolor="black",
           linewidth=0.4, label="Passing", zorder=2)

    # Drop bars: red for the gating cascade steps, yellow for the
    # trailing AF3 step.  Index 0 is the "All survivors" baseline
    # (always 0 drops), 1..len-2 are gating, last is AF3.
    n_steps = len(labels)
    drop_colours = (
        [COLOUR_FAIL] * (n_steps - 1)   # cascade-red for gating drops
        + [COLOUR_AF3_DROP]              # yellow for AF3 trailing drop
    )
    ax.bar(xs, drops, bottom=survivors, color=drop_colours,
           edgecolor="black", linewidth=0.4, alpha=0.85, zorder=2)

    # Annotate counts.
    for i, (s, d) in enumerate(zip(survivors, drops)):
        if s > 0:
            ax.text(i, s / 2, f"{s}", ha="center", va="center",
                    fontsize=10, color="white", fontweight="bold")
        if d > 0:
            ax.text(i, s + d / 2, f"−{d}", ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")

    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Number of survivors", fontsize=10)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.grid(True, linestyle="-", linewidth=0.3, alpha=0.4, color="grey")
    ax.set_axisbelow(True)

    # Custom 3-entry legend so the yellow AF3 bar is explained.
    legend_handles = [
        mpatches.Patch(color=COLOUR_PASS,
                        label="Passing"),
        mpatches.Patch(color=COLOUR_FAIL, alpha=0.85,
                        label="Dropped (gating)"),
        mpatches.Patch(color=COLOUR_AF3_DROP, alpha=0.85,
                        label="Dropped (AF3, informational)"),
    ]
    ax.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, 1.02), ncol=3,
               fontsize=9, framealpha=0.95)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Plot 3: Composite score vs each orthogonal metric (replaces the
#         old standalone Sc / BSA / interface_plddt distributions
#         AND the standalone ΔΔG distribution — both of those were
#         strip-plots that carried no information at small n).
# ═══════════════════════════════════════════════════════════════════════════

# (column_name, display_label, threshold, comparison ('>=' or '<='),
#  units, x-axis style: 'auto' or 'fixed01' for [0,1]-bounded metrics)
METRIC_VS_COMPOSITE_PANELS = [
    ("sc",                                    "Sc (shape complementarity)",
     ORTHOG_SC_MIN,    ">=", "",       "fixed01"),
    ("bsa",                                   "BSA (Å²)",
     ORTHOG_BSA_MIN,   ">=", "Å²",     "auto"),
    ("representative_interface_plddt_median", "Interface pLDDT",
     ORTHOG_PLDDT_MIN, ">=", "",       "fixed01"),
    ("rosetta_ddg",                           "Rosetta ΔΔG (REU)",
     DDG_BENNETT_REFERENCE, "<=", "REU", "auto"),
    ("af3_nomsa_best_ra_eff",                 "AF3-no-MSA best ra_eff",
     ORTHOG_AF3_RA_MAX, "<=", "Å",     "auto"),
    # (AF3-no-MSA best iPTM panel removed 2026-04-28 — out of scope
    # for "metric vs composite"; AF3 iPTM is included as a column in
    # the combined cohort summary instead.)
]


def plot_metrics_vs_composite(
    rows: List[Dict],
    cross_rows_by_seq: Dict[str, Dict],
    out_path: str,
) -> bool:
    """Multi-panel scatter: composite_score (Y) vs each orthogonal
    metric (X), one panel per metric.  Mirrors the negsteer
    composite-vs-confidence plot — composite is the primary ranker,
    so showing every orthogonal metric against it shows whether each
    orthogonal metric agrees with the upstream composite ranking
    (high composite + good metric = upper-corner cluster) or fights
    against it (high composite + bad metric).

    Coloured by passes_orthogonal_filters (green pass / dark grey
    fail).  Threshold drawn as vertical dashed line; the failing
    region shaded pale red.  Survivor names annotated when n is
    small enough.

    Composite score is read from cross_summary (joined on
    mpnn_sequence) since survivors_with_orthogonal_metrics.csv
    doesn't carry it directly."""
    if not rows:
        _make_empty_plot("No survivors to plot", out_path)
        return False

    n_panels = len(METRIC_VS_COMPOSITE_PANELS)
    n_cols = 3
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(4.2 * n_cols, 3.6 * n_rows))
    axes_flat = axes.flatten() if n_panels > 1 else [axes]

    # Collect composite scores once (joined from cross_summary).
    def composite_for(row: Dict) -> Optional[float]:
        seq = row.get("mpnn_sequence", "")
        cross = cross_rows_by_seq.get(seq)
        if not cross:
            return None
        return _try_float(cross.get("cross_composite_score", ""))

    any_data = False
    for ax, (col, label, threshold, op, units, x_style) in zip(
        axes_flat, METRIC_VS_COMPOSITE_PANELS,
    ):
        # Bucket points by tier rather than by passes_orthogonal_filters.
        # The threshold line + shaded region convey per-panel pass/fail
        # visually; tier carries the orthogonal cohort signal.
        # Pre-2026-04-29 history: this loop used to bucket by
        # _passes(row) (global pass-all-filters), which made every
        # point in every panel show "Fail" if any filter failed
        # globally — even in panels where the row's own metric was
        # fine.  Tier-colouring fixes that.
        pts_by_tier: Dict[str, List[Tuple[float, float, str]]] = {
            "A": [], "B": [], "C": [], "none": [],
        }
        for row in rows:
            v = _try_float(row.get(col, ""))
            comp = composite_for(row)
            if v is None or comp is None:
                continue
            tier = (row.get("cross_tier") or "none").strip() or "none"
            if tier not in pts_by_tier:
                tier = "none"
            pts_by_tier[tier].append(
                (v, comp, row.get("mpnn_sequence", "?"))
            )

        total_pts = sum(len(v) for v in pts_by_tier.values())
        if total_pts == 0:
            ax.text(0.5, 0.5, f"No {col}\ndata", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="grey")
            ax.set_axis_off()
            continue
        any_data = True

        all_x = [p[0] for v in pts_by_tier.values() for p in v]
        all_y = [p[1] for v in pts_by_tier.values() for p in v]

        # X-axis range: respect bounded vs unbounded metrics.
        if x_style == "fixed01":
            x_lo, x_hi = 0.0, 1.0
        else:
            x_min = min(all_x)
            x_max = max(all_x)
            if threshold is not None:
                x_min = min(x_min, threshold)
                x_max = max(x_max, threshold)
            pad = max((x_max - x_min) * 0.08, 0.01)
            x_lo, x_hi = x_min - pad, x_max + pad

        # Pale-red failing-region shade (when threshold defined).
        # Threshold info is conveyed via the x-axis label below; no
        # longer added to the per-panel legend (per-panel legends
        # were removed 2026-04-29 in favour of one shared figure-
        # level legend at the top of the grid, freeing the panel
        # interiors of overlapping legend boxes).
        if threshold is not None:
            if op == ">=":
                ax.axvspan(x_lo, threshold, color=COLOUR_FAIL,
                           alpha=0.07, zorder=0)
            else:  # op == "<="
                ax.axvspan(threshold, x_hi, color=COLOUR_FAIL,
                           alpha=0.07, zorder=0)
            ax.axvline(threshold, color=COLOUR_THRESHOLD,
                       linestyle="--", linewidth=1.2, alpha=0.85, zorder=1)

        # Plot one scatter per tier.  Stack A on top so high-quality
        # candidates aren't hidden in clustered regions.  Tier labels
        # are conveyed by the shared figure-level legend.
        for tier in ("none", "C", "B", "A"):
            tier_pts = pts_by_tier[tier]
            if not tier_pts:
                continue
            ax.scatter(
                [p[0] for p in tier_pts], [p[1] for p in tier_pts],
                c=COLOUR_TIER[tier], edgecolor="black",
                linewidth=0.4, s=60, alpha=0.9,
                zorder={"none": 3, "C": 4, "B": 5, "A": 6}[tier],
            )

        # Annotate names when n is small.
        if total_pts <= 12:
            for tier_pts in pts_by_tier.values():
                for x, y, name in tier_pts:
                    ax.annotate(_short_name(name), xy=(x, y),
                                xytext=(4, 3), textcoords="offset points",
                                fontsize=6.5, alpha=0.85)

        ax.set_xlim(x_lo, x_hi)
        # Modest y-pad — the per-panel legend that used to need clearance
        # at the top is gone, so the data can fill more of the panel.
        y_top = max(all_y)
        y_bot = min(all_y)
        y_pad = max((y_top - y_bot) * 0.06, 0.03)
        ax.set_ylim(y_bot - 0.05, y_top + y_pad)
        # Bake the threshold into the x-axis label (was previously a
        # legend entry).  e.g. "AF3 ra_eff (Å)  [thr ≤ 5]".  The unicode
        # ≤/≥ make this read more naturally than ascii <= / >=.
        op_glyph = {"<=": "≤", ">=": "≥"}.get(op, op or "")
        if threshold is not None and op_glyph:
            unit_str = f" {units}" if units else ""
            ax.set_xlabel(
                f"{label}  [thr {op_glyph} {threshold:g}{unit_str}]",
                fontsize=9,
            )
        else:
            ax.set_xlabel(label, fontsize=9)
        ax.set_ylabel("Composite score", fontsize=9)
        ax.grid(True, linestyle="-", linewidth=0.3, alpha=0.35, color="grey")
        ax.set_axisbelow(True)

    # Hide any unused panels (when n_panels not a multiple of n_cols).
    for i in range(n_panels, len(axes_flat)):
        axes_flat[i].set_axis_off()

    if not any_data:
        plt.close(fig)
        _make_empty_plot("No metric / composite data to plot", out_path)
        return False

    # Shared figure-level legend at the top.  Tier colour map +
    # threshold-line indicator.  Replaces the per-panel legends that
    # used to sit in the upper-right of each panel and overlap data.
    # Note: per-panel threshold values now live in the x-axis label
    # of each panel (see set_xlabel above).
    shared_handles = [
        mpatches.Patch(color=COLOUR_TIER["A"], label="Tier A"),
        mpatches.Patch(color=COLOUR_TIER["B"], label="Tier B"),
        mpatches.Patch(color=COLOUR_TIER["C"], label="Tier C"),
        mpatches.Patch(color=COLOUR_TIER["none"], label="No tier"),
        mlines.Line2D([0], [0], color=COLOUR_THRESHOLD, linestyle="--",
                      linewidth=1.5, label="Threshold"),
        mpatches.Patch(color=COLOUR_FAIL, alpha=0.20,
                       label="Failing region"),
    ]
    fig.legend(handles=shared_handles, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), ncol=6,
               fontsize=9, framealpha=0.95,
               bbox_transform=fig.transFigure)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# (Plot 4 — standalone ΔΔG distribution — REMOVED 2026-04-28.
#  ΔΔG is now one panel of plot_metrics_vs_composite above.  The
#  standalone strip-plot/histogram carried no information at small n
#  and duplicated content the new combined plot already shows.)
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# (Plot 5 — standalone per-survivor profile heatmap — REMOVED 2026-04-28.
#  All of its information is now carried by the combined cohort+orthogonal
#  summary plot below, which shows the same per-survivor profile but
#  joined with the upstream negsteer cohort context — strictly a superset.)
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# Plot 6: Combined cohort + orthogonal summary
# ═══════════════════════════════════════════════════════════════════════════
#
# One row per MPNN sequence in cross_summary (steered rows only —
# controls are excluded from the ranking by cross_sequence_summary.py).
# Sorted by cross_composite_score descending.  Each row shows the full
# profile across negsteer + orthogonal metrics; sequences that didn't
# reach the orthogonal stage have grey-hatched orthogonal cells.
#
# Tier (A/B/C/none) is rendered as a left-edge coloured stripe rather
# than a numeric column so the colour vocabulary stays consistent
# (red/green for pass/fail thresholds; light blue for no-threshold
# numerics; grey hatch for missing).

# Column spec (display order, left to right after the tier stripe):
#   (csv_column, source ('cross' / 'survivor' / 'special'),
#    display_label, threshold, comparison, formatter)
# - source 'cross'    → read from cross_summary, joined on mpnn_sequence
# - source 'survivor' → read from survivors_with_orthogonal_metrics.csv
# - source 'special'  → handled in code (e.g. n_mutations is stage-aware)
# - threshold=None means no pass/fail colouring; cell uses light-blue fill
# - formatter is one of 'auto' (3-sig-fig-ish), 'int', '0p2', '0p3', '0f', '1f'

COMBINED_SUMMARY_COLUMNS = [
    # Composite is the sort key — first numerical column.
    ("cross_composite_score",                    "cross",
     "composite",       None,               None,  "0p3"),
    # Boltz pose — receptor-aligned effector RMSD vs truth.
    ("representative_ra_eff_vs_truth_median",    "cross",
     "Boltz ra_eff",    5.0,                "<=",  "0p2"),
    # Receptor's own fold quality.  Sits next to Boltz ra_eff because
    # the two are conceptually paired: ra_eff measures relative
    # positioning (is the binding mode right?) while this measures
    # whether the receptor's own structure survived the steering.
    # Threshold matches the pipeline's existing intact-check at
    # boltz2_iterate_steering.py:1046 (≤5 Å).  See
    # boltz2_negative_steering.py:_compute_pose_metrics for the
    # underlying definition (whole-receptor Kabsch fit).
    ("representative_independent_receptor_rmsd_median",  "cross",
     "rec RMSD",        5.0,                "<=",  "0p2"),
    ("representative_true_jaccard_median",       "cross",
     "true_jaccard",    None,               None,  "0p2"),
    # Boltz confidence (the six selected)
    ("representative_complex_plddt_median",      "cross",
     "complex_plddt",   0.70,               ">=",  "0p2"),
    ("representative_iptm_median",               "cross",
     "iptm",            0.30,               ">=",  "0p2"),
    ("representative_ipae_median",               "cross",
     "ipae",            15.0,               "<=",  "0p2"),
    ("representative_pae_pass_frac_median",      "cross",
     "pae_pass_frac",   0.10,               ">=",  "0p2"),
    ("representative_interface_plddt_median",    "cross",
     "iface_plddt",     0.75,               ">=",  "0p2"),
    ("representative_ipsae_min_15_median",       "cross",
     "ipsae_min_15",    None,               None,  "0p2"),
    # Mutations carried by the representative's final prediction.
    # Stage-aware: cold_start → 0; steering → steered count;
    # reversion → post-reversion surviving count.  Handled in code.
    ("__rep_n_mutations__",                       "special",
     "n_mut",           None,               None,  "int"),
    # Orthogonal block — these are missing for sequences that didn't
    # reach the orthogonal stage; rendered as grey-hatched cells.
    # AF3 iptm column dropped 2026-04-29 (informational-only and
    # redundant with AF3 ra_eff for triage).
    ("af3_nomsa_best_ra_eff",                    "survivor",
     "AF3 ra_eff",      ORTHOG_AF3_RA_MAX,  "<=",  "0p2"),
    ("sc",                                        "survivor",
     "Sc",              ORTHOG_SC_MIN,      ">=",  "0p3"),
    ("bsa",                                       "survivor",
     "BSA",             ORTHOG_BSA_MIN,     ">=",  "0f"),
    ("rosetta_ddg",                               "survivor",
     "ΔΔG",             DDG_BENNETT_REFERENCE, "<=", "1f"),
    ("representative_weighted_jaccard_median",   "cross",
     "weighted_jacc",   None,               None,  "0p2"),
]

# Light-blue neutral fill for no-threshold numeric cells.
COLOUR_NEUTRAL = "#CFE3F3"


def _format_cell(v: float, fmt: str) -> str:
    if fmt == "int":
        return f"{int(round(v))}"
    if fmt == "0p2":
        return f"{v:.2f}"
    if fmt == "0p3":
        return f"{v:.3f}"
    if fmt == "0f":
        return f"{v:.0f}"
    if fmt == "1f":
        return f"{v:.1f}"
    # auto
    if abs(v) >= 100:
        return f"{v:.0f}"
    if abs(v) >= 10:
        return f"{v:.1f}"
    return f"{v:.3f}"


def _representative_n_mutations(cross_row: Dict) -> Optional[float]:
    """Stage-aware count of mutations carried by the representative
    sequence's FINAL prediction.

      cold_start → no reversion attempted AND steered_total_mutations==0
                 → 0
      steering   → no reversion attempted AND steered_total_mutations>0
                 → steered count (representative_total_mutations_median)
      reversion  → reversion ran (any reverted_* populated)
                 → post-reversion count (representative_reverted_total_mutations)

    Mirrors the classify_seed stage detection used elsewhere in the
    plot scripts so n_mut here lines up with the stage-coloured
    dots in the negsteer per-seed plots."""
    rev_verdict = (cross_row.get("representative_reversion_verdict") or "").strip()
    rev_total = _try_int(cross_row.get("representative_reverted_total_mutations", ""))
    steered_total = _try_int(cross_row.get("representative_total_mutations_median", ""))

    # Reversion ran iff verdict is populated OR reverted_total_mutations
    # has a value.  If so, the representative's final prediction is
    # the reverted one and n_mut = post-reversion count.
    if rev_verdict or rev_total is not None:
        return float(rev_total) if rev_total is not None else None

    # No reversion: steered prediction is the final.  Cold-start path
    # has steered_total_mutations==0; both cases are read from the same
    # column.
    if steered_total is not None:
        return float(steered_total)
    return None


def plot_combined_cohort_orthogonal_summary(
    cross_rows: List[Dict],
    survivor_rows_by_seq: Dict[str, Dict],
    out_path: str,
) -> bool:
    """The big one: one row per cross_summary sequence (steered only,
    sort by cross_composite_score desc), one column per metric across
    the negsteer + orthogonal cascades.  Sequences that never reached
    the orthogonal stage have grey-hatched cells in those columns.

    Cell colouring:
      - threshold pass → green
      - threshold fail → red
      - no threshold   → light-blue neutral fill (still shows the value)
      - missing data   → grey hatched
    Tier rendered as left-edge coloured stripe (A/B/C/none)."""
    # Filter to steered rows only (controls have row_type != 'steered').
    steered = [
        r for r in cross_rows
        if (r.get("row_type", "steered") or "steered") == "steered"
    ]
    if not steered:
        _make_empty_plot("No steered rows in cross_summary", out_path)
        return False

    # Sort by (tier, -composite) — strict tier-first, composite within.
    # Tier A → B → C → none.  Within each tier, descending composite.
    # Rows without composite go to the bottom of their tier.  Matches
    # cross_rank_by_composite in the CSV (which itself sorts tier-first
    # at cross_sequence_summary.py:_composite_rank_key).
    def _sort_key(r: Dict):
        tier = (r.get("cross_tier") or "none").strip() or "none"
        tier_order = _TIER_SORT_ORDER.get(tier, _TIER_SORT_ORDER["none"])
        comp = _try_float(r.get("cross_composite_score", ""))
        # None composite → sort to bottom within tier
        comp_key = -(comp if comp is not None else float("-inf"))
        return (tier_order, comp_key)
    steered.sort(key=_sort_key)

    n_seq = len(steered)
    n_metr = len(COMBINED_SUMMARY_COLUMNS)

    # Layout: tier stripe (narrow), then metric columns.
    cell_w = 1.4
    cell_h = 1.0
    stripe_w = 0.35

    fig_w = max(10.0, stripe_w + n_metr * cell_w + 4.0)
    fig_h = max(3.0, 0.42 * n_seq + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for r_idx, cross_row in enumerate(steered):
        seq = cross_row.get("mpnn_sequence", "?")
        survivor = survivor_rows_by_seq.get(seq)
        y = n_seq - 1 - r_idx  # top-down ordering

        # Tier stripe.
        tier = (cross_row.get("cross_tier") or "none").strip() or "none"
        stripe_colour = COLOUR_TIER.get(tier, COLOUR_TIER["none"])
        ax.add_patch(mpatches.Rectangle(
            (0, y), stripe_w, cell_h,
            facecolor=stripe_colour, edgecolor="black", linewidth=0.4,
            alpha=0.85,
        ))

        for c_idx, (col, src, _label, threshold, op, fmt) in enumerate(
            COMBINED_SUMMARY_COLUMNS,
        ):
            x = stripe_w + c_idx * cell_w

            # Source the value.
            if src == "special" and col == "__rep_n_mutations__":
                v = _representative_n_mutations(cross_row)
            elif src == "survivor":
                v = (_try_float(survivor.get(col, "")) if survivor else None)
            else:  # 'cross'
                v = _try_float(cross_row.get(col, ""))

            if v is None:
                # Grey hatched: missing.
                ax.add_patch(mpatches.Rectangle(
                    (x, y), cell_w, cell_h,
                    facecolor="white", edgecolor="grey",
                    linewidth=0.5, hatch="///", alpha=0.45,
                ))
                ax.text(x + cell_w / 2, y + cell_h / 2, "—",
                        ha="center", va="center", fontsize=8, color="grey")
                continue

            # Determine fill: red/green if threshold, light blue otherwise.
            if threshold is None:
                fill = COLOUR_NEUTRAL
                alpha = 0.85
            else:
                if op == ">=":
                    p = v >= threshold
                else:
                    p = v <= threshold
                fill = COLOUR_PASS if p else COLOUR_FAIL
                alpha = 0.55

            ax.add_patch(mpatches.Rectangle(
                (x, y), cell_w, cell_h,
                facecolor=fill, edgecolor="black",
                linewidth=0.4, alpha=alpha,
            ))
            ax.text(x + cell_w / 2, y + cell_h / 2,
                    _format_cell(v, fmt),
                    ha="center", va="center", fontsize=8,
                    color="black", fontweight="bold")

    # Y-axis: sequence labels (with tier prefix marker for redundancy).
    ax.set_yticks([n_seq - 1 - i + 0.5 for i in range(n_seq)])
    yticklabels = []
    for cross_row in steered:
        tier = (cross_row.get("cross_tier") or "none").strip() or "none"
        seq = cross_row.get("mpnn_sequence", "?")
        # Tier letter + short name.  Tier letter is redundant with
        # the stripe but useful for B/W printing or copy-paste.
        marker = "·" if tier == "none" else tier
        yticklabels.append(f"{marker}  {_short_name(seq)}")
    ax.set_yticklabels(yticklabels, fontsize=8.5)

    # X-axis: metric labels at top, with threshold sub-label.
    ax.set_xticks([
        stripe_w + i * cell_w + cell_w / 2
        for i in range(n_metr)
    ])
    xtick_labels = []
    for _col, _src, label, threshold, op, _fmt in COMBINED_SUMMARY_COLUMNS:
        if threshold is None:
            xtick_labels.append(label)
        else:
            xtick_labels.append(f"{label}\n({op}{threshold:g})")
    ax.set_xticklabels(xtick_labels, fontsize=7.5, rotation=30, ha="left")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")

    # Tier-stripe header label.
    ax.text(stripe_w / 2, n_seq + 0.15, "tier",
            ha="center", va="bottom", fontsize=7.5, style="italic")

    ax.set_xlim(0, stripe_w + n_metr * cell_w)
    ax.set_ylim(0, n_seq + 0.6)

    # Legend below.
    legend_handles = [
        mpatches.Patch(color=COLOUR_PASS, alpha=0.55, label="Passes threshold"),
        mpatches.Patch(color=COLOUR_FAIL, alpha=0.55, label="Fails threshold"),
        mpatches.Patch(color=COLOUR_NEUTRAL, alpha=0.85,
                       label="No threshold (informational)"),
        mpatches.Patch(facecolor="white", edgecolor="grey", hatch="///",
                       label="Metric missing"),
        mpatches.Patch(color=COLOUR_TIER["A"], alpha=0.85, label="Tier A"),
        mpatches.Patch(color=COLOUR_TIER["B"], alpha=0.85, label="Tier B"),
        mpatches.Patch(color=COLOUR_TIER["C"], alpha=0.85, label="Tier C"),
        mpatches.Patch(color=COLOUR_TIER["none"], alpha=0.85, label="No tier"),
    ]
    ax.legend(handles=legend_handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.04), ncol=4,
              fontsize=8, framealpha=0.95)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Top-level orchestration
# ═══════════════════════════════════════════════════════════════════════════

def _load_rows(path: Path, label: str) -> List[Dict]:
    print(f"Reading {label}: {path}")
    rows: List[Dict] = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    print(f"  rows: {len(rows)}")
    return rows


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--survivors-csv", required=True, type=Path,
        help="survivors_with_orthogonal_metrics.csv from "
             "NEGSTEER_ORTHOGONAL_METRICS.")
    ap.add_argument("--cross-summary-csv", required=True, type=Path,
        help="cross_sequence_summary.csv from upstream "
             "NEGSTEER_CROSS_SEQUENCE.  Required: provides composite "
             "scores and the full upstream cohort needed by the "
             "metrics-vs-composite and combined-summary plots.")
    ap.add_argument("--outdir", required=True,
        help="Directory to write plots into; created if missing.")
    return ap.parse_args()


def main() -> int:
    args = _parse_args()
    if not args.survivors_csv.exists():
        raise SystemExit(
            f"--survivors-csv path does not exist: {args.survivors_csv}"
        )
    if not args.cross_summary_csv.exists():
        raise SystemExit(
            f"--cross-summary-csv path does not exist: {args.cross_summary_csv}"
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    survivor_rows = _load_rows(args.survivors_csv, "survivors")
    if not survivor_rows:
        print("No survivor rows — nothing to plot.")
        return 0

    n_passing = sum(1 for r in survivor_rows if _passes(r))
    print(f"  passing all orthogonal filters: {n_passing}/{len(survivor_rows)}")

    cross_rows = _load_rows(args.cross_summary_csv, "cross-summary")
    print(f"  upstream cohort: {len(cross_rows)} rows "
          f"(includes controls + tier-none stubs)")

    # Lookup tables for join-by-mpnn_sequence.
    cross_by_seq = {
        r.get("mpnn_sequence", ""): r for r in cross_rows if r.get("mpnn_sequence")
    }
    survivor_by_seq = {
        r.get("mpnn_sequence", ""): r for r in survivor_rows
        if r.get("mpnn_sequence")
    }

    plots = [
        ("AF3 vs Boltz ra_eff scatter (cross-model agreement)",
         "orthogonal_af3_vs_boltz.png",
         lambda p: plot_af3_vs_boltz(survivor_rows, p)),
        ("Filter-cascade waterfall (AF3 informational)",
         "orthogonal_filter_cascade.png",
         lambda p: plot_filter_cascade(survivor_rows, p)),
        ("Composite vs each orthogonal metric "
         "(Sc / BSA / iface_plddt / ΔΔG / AF3 ra_eff)",
         "orthogonal_metrics_vs_composite.png",
         lambda p: plot_metrics_vs_composite(
             survivor_rows, cross_by_seq, p,
         )),
        ("Combined cohort + orthogonal summary (every sequence, "
         "ranked by composite)",
         "orthogonal_combined_cohort_summary.png",
         lambda p: plot_combined_cohort_orthogonal_summary(
             cross_rows, survivor_by_seq, p,
         )),
    ]

    for desc, fname, fn in plots:
        out = outdir / fname
        print(f"Building {desc}...")
        try:
            ok = fn(str(out))
            if ok:
                print(f"  → {out}")
            else:
                print("  (skipped — no data)")
        except Exception as e:  # noqa: BLE001
            print(f"  FAILED: {type(e).__name__}: {e}")

    print(f"\nDone.  Inspect: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())