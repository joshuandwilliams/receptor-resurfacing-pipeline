#!/usr/bin/env python3
"""
test_rosetta_filtering_plots.py
-------------------------------
Standalone iteration harness for ``rosetta_filter_plots.py``.

Reads the metrics JSON produced by a previous ``test_rosetta_filtering``
run and re-emits the diagnostic plot in an output directory of choice.
No GPU, no Nextflow — designed for fast iteration of the plot code
itself against cached test outputs.

Background
----------
The original plot was a 2D scatter of shape complementarity (Sc, x)
vs. dG_separated (y).  But at this point in the pipeline the input
PDBs are RFDiffusion polyvaline backbones with no real sidechains, so
``dG_separated`` is dominated by valine-clash artefacts and is *not*
physically meaningful (see ``rosetta_filter_collect.py`` docstring).
The pipeline correctly uses Sc — purely geometric — as the filter
gate.  Plotting dG_separated alongside it implies a meaning that
isn't there.

So this iteration replaces the scatter with a 1D Sc histogram, with
the filtered region (Sc < threshold) shaded red and extending fully
to the left axis spine.  The "best model" star marker is dropped:
with no second dimension to plot against, "best model" is just the
rightmost bar of the histogram — redundant with the histogram itself,
and removing it kills the gnarliest code in the original file (the
post-render annotation bbox measurement and x-limit extension).

Usage
-----
    python test_rosetta_filtering_plots.py \\
        --metrics tests/rosetta_filtering/receptor_resurfacing_results/rosetta_filtering/rosetta_filter_metrics.json \\
        --outdir  tests/rosetta_filtering/receptor_resurfacing_results/plots_iter

If ``--metrics`` is omitted the script falls back to the canonical
test output path.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── Styling (matches production rosetta_filter_plots.py) ───────────────────
COLOUR_ALL  = "#4C72B0"   # blue (histogram bars — passing region)
COLOUR_WARN = "#FFCCCC"   # pink (filtered region shading)
COLOUR_FAIL = "#D62728"   # red  (filtered-region edge / threshold line)


# ═══════════════════════════════════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════════════════════════════════

def plot_sc_histogram(designs, threshold, out_path):
    """Render the Sc histogram with the filtered region shaded red.

    Strategy:
    - Compute Sc values, drop missing.
    - Set explicit xlim to the data range with a small margin, BEFORE
      drawing the axvspan, so the red region's left edge is the axis
      spine and matplotlib's autoscale doesn't shift it.
    - Histogram with thin black borders on bars.
    - axvspan from (xlim_left, threshold) — guaranteed to touch the
      left spine because we set xlim first.
    - Vertical line at threshold for clarity.
    - Legend with three entries: histogram, filtered region, threshold.
    """
    sc_vals = [d["sc_value"] for d in designs if d.get("sc_value") is not None]
    if not sc_vals:
        _make_empty_plot("No Sc data available", out_path)
        return False

    sc_vals = np.array(sc_vals)
    n_total = len(sc_vals)
    n_pass  = int(np.sum(sc_vals >= threshold))
    n_fail  = n_total - n_pass

    fig, ax = plt.subplots(figsize=(8, 5))

    # ── X-axis limits set FIRST so axvspan can use them deterministically.
    sc_min = float(sc_vals.min())
    sc_max = float(sc_vals.max())
    margin = max((sc_max - sc_min) * 0.10, 0.02)
    xlim_left  = min(sc_min - margin, threshold - margin)
    xlim_right = max(sc_max + margin, threshold + margin)
    ax.set_xlim(xlim_left, xlim_right)

    # ── Filtered region shading.  Use the explicit xlim_left so the
    # span is flush against the left axis spine — no white gap.
    ax.axvspan(xlim_left, threshold, color=COLOUR_WARN, alpha=0.5,
               zorder=1,
               label=f"Sc < {threshold:.2f} (filtered, n={n_fail} of {n_total})")

    # ── Histogram.  Bin edges chosen by Freedman–Diaconis when there's
    # enough data, otherwise a sensible fixed count.
    if n_total >= 12:
        # Freedman–Diaconis: bin width = 2 * IQR / n^(1/3)
        iqr = float(np.percentile(sc_vals, 75) - np.percentile(sc_vals, 25))
        if iqr > 0:
            bin_w = 2.0 * iqr / (n_total ** (1.0 / 3.0))
            n_bins = max(8, min(40, int(np.ceil((sc_max - sc_min) / bin_w))))
        else:
            n_bins = 12
    else:
        n_bins = max(5, min(n_total, 10))

    ax.hist(sc_vals, bins=n_bins, range=(xlim_left, xlim_right),
            color=COLOUR_ALL, edgecolor="black", linewidth=0.5,
            alpha=0.85, zorder=2)

    # ── Threshold marker line (sits on top of histogram + shading).
    ax.axvline(threshold, color=COLOUR_FAIL, linestyle="--", linewidth=1.2,
               zorder=3, label=f"Threshold = {threshold:.2f}")

    # ── Axis cosmetics.
    ax.set_xlabel("Shape Complementarity (Sc)", fontsize=12)
    ax.set_ylabel("Number of designs", fontsize=12)
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    ax.set_axisbelow(True)

    # Y-axis only ever shows integer counts.
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # ── Legend in the upper right.  Pad y-axis up so the legend never
    # sits on top of the tallest histogram bar.
    y_top = ax.get_ylim()[1]
    ax.set_ylim(top=y_top * 1.20)
    ax.legend(fontsize=10, loc="upper right", framealpha=0.95)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def _make_empty_plot(message, out_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center",
            transform=ax.transAxes, fontsize=14, color="grey")
    ax.set_axis_off()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Top-level orchestration
# ═══════════════════════════════════════════════════════════════════════════

def _resolve_metrics_path(supplied):
    if supplied:
        p = Path(supplied)
        if not p.exists():
            raise SystemExit(f"--metrics path does not exist: {p}")
        return p
    fallback = Path(
        "tests/rosetta_filtering/receptor_resurfacing_results/"
        "rosetta_filtering/rosetta_filter_metrics.json"
    )
    if fallback.exists():
        return fallback
    raise SystemExit(
        f"No --metrics passed and fallback {fallback} does not exist. "
        f"Run test_rosetta_filtering first or pass --metrics explicitly."
    )


def _parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--metrics", default=None,
        help="Path to rosetta_filter_metrics.json from a test_rosetta_filtering run.")
    parser.add_argument("--outdir", required=True,
        help="Directory to write plots into; created if missing.")
    return parser.parse_args()


def main():
    args = _parse_args()
    metrics_path = _resolve_metrics_path(args.metrics)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading metrics: {metrics_path}")
    data = json.loads(metrics_path.read_text())
    designs = data.get("designs", [])
    threshold = data.get("sc_threshold", 0.62)
    print(f"  designs:   {len(designs)}")
    print(f"  threshold: {threshold}")

    if not designs:
        print("No designs in metrics — nothing to plot.")
        return 0

    out = outdir / "rosetta_sc_histogram.png"
    print("Building Sc histogram...")
    try:
        ok = plot_sc_histogram(designs, threshold, str(out))
        if ok:
            print(f"  → {out}")
        else:
            print(f"  (no Sc data)")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

    print(f"\nDone.  Inspect: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
