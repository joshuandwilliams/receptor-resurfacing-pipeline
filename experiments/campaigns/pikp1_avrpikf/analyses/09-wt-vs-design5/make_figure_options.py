#!/usr/bin/env python3
"""Compose the wild-type versus d26_s1 panels rendered by the ChimeraX script.

render_wt_vs_design5.cxc writes three renders at one camera and one window
size, so every panel here is size-matched by construction. That was the defect
in the previous version of this figure, which was two screen captures taken at
different scales.

    option1_side_by_side   the two complexes side by side, with a key
    option2_with_overlay   the same pair plus a superposition panel
    option3_overlay_only   the superposition alone, most compact

Usage:
    make_figure_options.py [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
RENDERS = os.path.join(HERE, "renders")
INK, MUTED = "#222222", "#666666"

KEY = [("#6495ED", "Fixed Pikp-1 HMA backbone"),
       ("#FF0000", "Design region 1"),
       ("#940094", "Design region 2"),
       ("#DAA520", "Negative-steering mutations"),
       ("#90EE90", "AVR-PikF effector")]


def apply_style():
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 300,
                         "savefig.bbox": "tight", "font.size": 9,
                         "text.color": INK})


def joint_bbox(names, pad=14):
    """One crop box covering the ink in every render.

    The renders are saved on a transparent background at one camera, so the
    alpha channel gives the drawn extent directly. Cropping every panel to the
    SAME box keeps them size-matched, which is the whole point of rendering
    them at one camera, while still removing the dead margin.
    """
    y0 = x0 = 10 ** 9
    y1 = x1 = -1
    for n in names:
        a = mpimg.imread(os.path.join(RENDERS, n))
        ink = a[:, :, 3] > 0.02 if a.shape[2] == 4 else a[:, :, :3].min(2) < 0.98
        ys, xs = np.where(ink)
        y0, x0 = min(y0, ys.min()), min(x0, xs.min())
        y1, x1 = max(y1, ys.max()), max(x1, xs.max())
    h, w = a.shape[:2]
    return (max(0, y0 - pad), min(h, y1 + pad + 1),
            max(0, x0 - pad), min(w, x1 + pad + 1))


def panel(ax, name, title, box=None):
    a = mpimg.imread(os.path.join(RENDERS, name))
    if box:
        y0, y1, x0, x1 = box
        a = a[y0:y1, x0:x1]
    ax.imshow(a)
    ax.set_axis_off()
    ax.set_title(title, fontsize=10, pad=4)


def add_key(fig, ncol=3, y=0.045):
    """Key centred on the FIGURE, not on one axes.

    Anchoring it to a single axes put it under that panel rather than under the
    pair, which is why it read as off-centre.
    """
    handles = [Patch(facecolor=c, edgecolor="none", label=l) for c, l in KEY]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               ncol=ncol, frameon=False, fontsize=8.5)


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print("wrote", path)


def option1_side_by_side(outdir):
    names = ["panel_a_wt.png", "panel_b_d26_s1.png"]
    box = joint_bbox(names)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4))
    panel(axes[0], names[0], "Wild-type Pikp-1 HMA / AVR-PikF", box)
    panel(axes[1], names[1], "d26_s1 / AVR-PikF, negative steered", box)
    for ax, letter in zip(axes, "AB"):
        ax.annotate(letter, xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(-2, 10), textcoords="offset points",
                    fontsize=13, weight="bold", va="top")
    # Tight side by side, with room only for the key beneath.
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.11,
                        wspace=0.02)
    add_key(fig, ncol=3)
    save(fig, outdir, "option1_side_by_side")


def option2_with_overlay(outdir):
    names = ["panel_a_wt.png", "panel_b_d26_s1.png", "panel_c_overlay.png"]
    box = joint_bbox(names)
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.1))
    for ax, n, t in zip(axes, names, ["Wild type", "d26_s1", "Superposed"]):
        panel(ax, n, t, box)
    for ax, letter in zip(axes, "ABC"):
        ax.annotate(letter, xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(-2, 10), textcoords="offset points",
                    fontsize=13, weight="bold", va="top")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.12,
                        wspace=0.02)
    add_key(fig, ncol=5)
    save(fig, outdir, "option2_with_overlay")


def option3_overlay_only(outdir):
    box = joint_bbox(["panel_c_overlay.png"])
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    panel(ax, "panel_c_overlay.png",
          "d26_s1 superposed on the wild-type complex", box)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.14)
    add_key(fig, ncol=2)
    save(fig, outdir, "option3_overlay_only")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    option1_side_by_side(args.outdir)
    option2_with_overlay(args.outdir)
    option3_overlay_only(args.outdir)


if __name__ == "__main__":
    main()
