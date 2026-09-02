#!/usr/bin/env python3
"""Build the thesis panels for the d39 figure from the hand-positioned renders.

The camera and every label position were set by hand in ChimeraX and saved to
positioned.cxs; renders/ holds the three panels written from that session. This
script only crops and lays them out. It draws no A/B letters, because those
belong to the thesis SVG template and are added by assemble_figure.py.

Panel A is one structure and panel B is two. To keep every structure the SAME
size on the page, both panels are written at the same pixel width, TWO
structures wide, with the wild type centred in its half-empty panel. The
assembler scales each panel to one shared display width, so equal panel widths
mean equal structure sizes. Writing panel A one structure wide instead would
scale it up to the full page width and make the wild type twice the size of
either design.

Outputs, into thesis-figures/:
    panel_a.png    wild type, centred, two structures wide
    panel_b.png    d39_s0 and d39_s1 side by side
    key.png        shared colour key, no letter

Usage:
    build_thesis_panels.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RENDERS = os.path.join(HERE, "renders")
OUT = os.path.join(HERE, "thesis-figures")

INK = "#222222"
WT = "panel_a_wt.png"
S0 = "panel_b_d39_s0.png"
S1 = "panel_c_d39_s1.png"

# Design region 2 is outside the framing, so it is not in the key. Arg18 shares
# the helix colour deliberately and is named in the caption rather than here.
KEY = [("#6495ED", "Fixed Pikp-1 HMA backbone"),
       ("#FF0000", "Design region 1"),
       ("#D2B48C", "Helix A, residues 14-26"),
       ("#2F4F4F", "Design region 1 side chains")]


def joint_bbox(names, pad=14):
    """One crop box covering the ink in every render.

    The renders are transparent-background and share one camera, so cropping
    them all to the same box keeps them size-matched while removing the dead
    margin.
    """
    y0 = x0 = 10 ** 9
    y1 = x1 = -1
    for n in names:
        a = mpimg.imread(os.path.join(RENDERS, n))
        ink = a[:, :, 3] > 0.02
        ys, xs = np.where(ink)
        y0, x0 = min(y0, ys.min()), min(x0, xs.min())
        y1, x1 = max(y1, ys.max()), max(x1, xs.max())
    h, w = a.shape[:2]
    return (max(0, y0 - pad), min(h, y1 + pad + 1),
            max(0, x0 - pad), min(w, x1 + pad + 1))


def crop(name, box):
    y0, y1, x0, x1 = box
    return mpimg.imread(os.path.join(RENDERS, name))[y0:y1, x0:x1]


def panel(ax, img, title):
    ax.imshow(img)
    ax.set_axis_off()
    ax.set_title(title, fontsize=10, pad=4, color=INK)


def main():
    os.makedirs(OUT, exist_ok=True)
    # bbox is NOT tight: panel A centres the wild type using two empty grid
    # columns, and a tight box would crop exactly that padding away, leaving A
    # narrower than B and so scaled larger by the assembler.
    plt.rcParams.update({"savefig.dpi": 300, "savefig.bbox": None,
                         "font.size": 9, "text.color": INK})
    box = joint_bbox([WT, S0, S1])
    aspect = (box[3] - box[2]) / (box[1] - box[0])
    panel_w, panel_h = 5.5, 5.5 / aspect

    # Panel A, two structures wide with the wild type in the middle.
    fig = plt.figure(figsize=(panel_w * 2, panel_h))
    gs = fig.add_gridspec(1, 4)
    ax = fig.add_subplot(gs[0, 1:3])
    panel(ax, crop(WT, box), "Wild-type Pikp-1 HMA")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.93, bottom=0.005)
    fig.savefig(os.path.join(OUT, "panel_a.png"), transparent=True, bbox_inches=None)
    plt.close(fig)

    # Panel B, the two designs side by side.
    fig, axes = plt.subplots(1, 2, figsize=(panel_w * 2, panel_h))
    panel(axes[0], crop(S0, box), "d39_s0, not autoactive")
    panel(axes[1], crop(S1, box), "d39_s1, autoactive")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.93, bottom=0.005,
                        wspace=0.01)
    fig.savefig(os.path.join(OUT, "panel_b.png"), transparent=True, bbox_inches=None)
    plt.close(fig)

    # Shared key, at the panels' width so the assembler does not upscale it.
    fig = plt.figure(figsize=(panel_w * 2, 0.5))
    fig.legend(handles=[Patch(facecolor=c, edgecolor="none", label=l)
                        for c, l in KEY],
               loc="center", ncol=2, frameon=False, fontsize=9)
    fig.savefig(os.path.join(OUT, "key.png"), transparent=True, bbox_inches=None)
    plt.close(fig)

    for n in ("panel_a.png", "panel_b.png", "key.png"):
        print("wrote", os.path.join(OUT, n))


if __name__ == "__main__":
    main()
