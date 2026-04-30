#!/usr/bin/env python3
"""
rosetta_filter_plots.py
-----------------------
Diagnostic plots from the Rosetta pre-validation physics filter.

Reads rosetta_filter_metrics.json and generates:

    rosetta_sc_histogram.png
        Histogram of shape complementarity (sc_value) across all
        designs.  The Sc < threshold region is shaded red and extends
        flush to the left axis spine.  A dashed vertical line marks
        the threshold.  The legend reports the filtered count as
        "n=<filtered> of <total>".

Why a histogram and not a Sc-vs-dG_separated scatter?
-----------------------------------------------------
At this stage the input PDBs are RFDiffusion polyvaline backbones
with no real sidechains, so Rosetta's energy and packing terms
(``dG_separated``, ``packstat``, ``delta_unsatHbonds``) are
dominated by valine-clash artefacts and are NOT physically
meaningful.  See the docstring caveat in
``rosetta_filter_collect.py``.

The pipeline correctly uses ``sc_value`` (purely geometric) as the
filter gate at this stage; the energy terms only become meaningful
after ProteinMPNN sequence design and downstream relaxation.
Plotting them alongside Sc here would imply meaning that isn't
there, so the diagnostic plot shows just the Sc distribution.

Outputs:
    rosetta_sc_histogram.png
"""

import argparse
import json
import os
import sys

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# ── Styling ──────────────────────────────────────────────────────────────────
COLOUR_ALL  = "#4C72B0"   # blue (histogram bars)
COLOUR_WARN = "#FFCCCC"   # pink (filtered region shading)
COLOUR_FAIL = "#D62728"   # red  (threshold line)


# ── Output filenames (single source of truth) ────────────────────────────────
ALL_PLOT_FILES = [
    "rosetta_sc_histogram.png",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metrics",
        required=True,
        help="Path to rosetta_filter_metrics.json from ROSETTA_FILTER",
    )
    return parser.parse_args()


def make_empty_plot(message, path):
    """Create a placeholder plot with a centred message."""
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(
        0.5, 0.5, message,
        ha="center", va="center", transform=ax.transAxes,
        fontsize=14, color="grey",
    )
    ax.set_axis_off()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def save_fallback_plots(message):
    """Create placeholder plots when no data is available."""
    for name in ALL_PLOT_FILES:
        make_empty_plot(message, name)


# ── Plot: Shape complementarity histogram ────────────────────────────────────

def plot_sc_histogram(designs, threshold):
    """
    Histogram of shape complementarity (Sc) across all designs.

    The Sc < threshold region is shaded red and extends flush to the
    left axis spine — achieved by setting ``ax.set_xlim`` BEFORE the
    ``axvspan`` call so the span's left edge isn't subsequently
    shifted by matplotlib's autoscaling.

    A dashed vertical line marks the threshold.  The legend reports
    ``Sc < <threshold> (filtered, n=<filtered> of <total>)`` plus the
    threshold marker — no separate "all designs" entry, which would
    just be a redundant blue patch.
    """
    out_path = "rosetta_sc_histogram.png"

    sc_vals = [d["sc_value"] for d in designs if d.get("sc_value") is not None]
    if not sc_vals:
        make_empty_plot("No Sc data available", out_path)
        return

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
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # ── Legend in the upper right.  Pad y-axis up so the legend never
    # sits on top of the tallest histogram bar.
    y_top = ax.get_ylim()[1]
    ax.set_ylim(top=y_top * 1.20)
    ax.legend(fontsize=10, loc="upper right", framealpha=0.95)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not HAS_MPL:
        for name in ALL_PLOT_FILES:
            with open(name, "w") as f:
                f.write("matplotlib not available")
        sys.exit(0)

    if not os.path.exists(args.metrics) or os.path.getsize(args.metrics) == 0:
        save_fallback_plots("No metrics data")
        sys.exit(0)

    with open(args.metrics) as f:
        data = json.load(f)

    designs = data.get("designs", [])
    threshold = data.get("sc_threshold", 0.62)

    if not designs:
        save_fallback_plots("No designs in metrics")
        sys.exit(0)

    plot_sc_histogram(designs, threshold)


if __name__ == "__main__":
    main()