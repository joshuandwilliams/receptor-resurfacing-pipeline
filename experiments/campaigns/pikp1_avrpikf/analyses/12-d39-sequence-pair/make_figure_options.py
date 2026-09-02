#!/usr/bin/env python3
"""Compose the d39_s0 / d39_s1 panels rendered by the ChimeraX script.

render_d39_pair.cxc writes four renders at one camera and one window size, so
every panel here is size-matched by construction. Same approach as
09-wt-vs-design5/make_figure_options.py.

    option1_three_panel   wild type, d39_s0, d39_s1 side by side, with a key
    option2_designs_only  the two designs only, dropping the wild type
    option3_overlay       all three superposed in one panel

Usage:
    make_figure_options.py [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RENDERS = os.path.join(HERE, "renders")
INK, MUTED = "#222222", "#666666"

# Same colours as 09-wt-vs-design5, plus helix A, which is what this figure is
# about. The effector is absent here, so it is not in the key.
KEY = [("#6495ED", "Fixed Pikp-1 HMA backbone"),
       ("#FF0000", "Design region 1"),
       ("#D2B48C", "Helix A, residues 14-26"),
       ("#2F4F4F", "Side chains differing between the two sequences")]


def apply_style():
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 300,
                         "savefig.bbox": "tight", "font.size": 9,
                         "text.color": INK})


def joint_bbox(names, pad=14):
    """One crop box covering the ink in every render.

    The renders are saved on a transparent background at one camera, so the
    alpha channel gives the drawn extent directly. Cropping every panel to the
    SAME box keeps them size-matched, which is the point of rendering them at
    one camera, while still removing the dead margin.
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


def keep_bottom(box, fraction=0.5):
    """Trim the crop box to its lower part.

    The top of the domain is the same beta sheet in all three panels and
    carries nothing the figure is about, so it is cropped away here rather
    than in ChimeraX. Cropping after the render keeps one camera and one
    joint box, so the panels stay size-matched.
    """
    y0, y1, x0, x1 = box
    return (int(y1 - (y1 - y0) * fraction), y1, x0, x1)


def fig_size(box, n_cols=1, n_rows=1, panel_w=6.6, key_rows=2):
    """Figure size that matches the crop aspect, so no panel carries dead space.

    imshow holds the image aspect fixed, so an axes wider or taller than the
    crop pads the difference with whitespace. Sizing the figure from the crop
    instead means the panels sit directly against each other.
    """
    y0, y1, x0, x1 = box
    aspect = (x1 - x0) / (y1 - y0)
    panel_h = panel_w / aspect
    return (panel_w * n_cols, panel_h * n_rows + 0.34 * key_rows)


def panel(ax, name, title, box=None):
    a = mpimg.imread(os.path.join(RENDERS, name))
    if box:
        y0, y1, x0, x1 = box
        a = a[y0:y1, x0:x1]
    ax.imshow(a)
    ax.set_axis_off()
    ax.set_title(title, fontsize=10, pad=4)


def letters(axes):
    for ax, letter in zip(axes, "ABCD"):
        ax.annotate(letter, xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(-2, 10), textcoords="offset points",
                    fontsize=13, weight="bold", va="top")


def add_key(fig, ncol=3, y=0.0):
    handles = [Patch(facecolor=c, edgecolor="none", label=l) for c, l in KEY]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, y),
               ncol=ncol, frameon=False, fontsize=8.5)


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print("wrote", path)


def option1_three_panel(outdir):
    """Wild type on top, the two designs side by side beneath as one panel B.

    Cropped to the lower half, where design region 1 and helix A sit. Stacking
    all three made the figure too tall for the page, and the two designs are
    the comparison, so they share a letter and sit next to each other.
    """
    names = ["panel_a_wt.png", "panel_b_d39_s0.png", "panel_c_d39_s1.png"]
    box = joint_bbox(names)
    fig = plt.figure(figsize=fig_size(box, 2, 2))
    # 2 x 4 so the single top panel can be centred over the pair beneath while
    # every panel keeps the same width.
    gs = fig.add_gridspec(2, 4)
    ax_a = fig.add_subplot(gs[0, 1:3])
    ax_b = fig.add_subplot(gs[1, 0:2])
    ax_c = fig.add_subplot(gs[1, 2:4])
    panel(ax_a, names[0], "Wild-type Pikp-1 HMA", box)
    panel(ax_b, names[1], "d39_s0, not autoactive", box)
    panel(ax_c, names[2], "d39_s1, autoactive", box)
    letters([ax_a, ax_b])
    fig.subplots_adjust(left=0.005, right=0.995, top=0.95, bottom=0.10,
                        wspace=0.02, hspace=0.14)
    add_key(fig, ncol=2)
    save(fig, outdir, "option1_three_panel")


def option2_designs_only(outdir):
    names = ["panel_b_d39_s0.png", "panel_c_d39_s1.png"]
    box = joint_bbox(names)
    fig, axes = plt.subplots(2, 1, figsize=fig_size(box, 1, 2))
    panel(axes[0], names[0], "d39_s0, not autoactive", box)
    panel(axes[1], names[1], "d39_s1, autoactive", box)
    letters(axes)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.95, bottom=0.11,
                        hspace=0.12)
    add_key(fig, ncol=2)
    save(fig, outdir, "option2_designs_only")


def option3_overlay(outdir):
    names = ["panel_d_overlay.png"]
    box = joint_bbox(names)
    fig, ax = plt.subplots(1, 1, figsize=fig_size(box, 1, 1, panel_w=6.0))
    panel(ax, names[0], "All three superposed on the receptor", box)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.16)
    add_key(fig, ncol=2)
    save(fig, outdir, "option3_overlay")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    option1_three_panel(args.outdir)
    option2_designs_only(args.outdir)


if __name__ == "__main__":
    main()
