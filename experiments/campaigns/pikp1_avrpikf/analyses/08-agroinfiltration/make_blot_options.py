#!/usr/bin/env python3
"""Reformat the AVR-PikF cohort western blot into the Banfield lab house style.

The house style, taken from the blot panels in Maidment (2019) figures 5.3-5.5
and De la Concepcion (2020), is:

    * the coloured ladder lane cropped out, replaced by numeric kDa marks with
      short ticks down the left-hand side
    * one tight strip per molecular-weight region rather than one tall panel
      spanning the whole gel
    * the antibody on the right of each strip, with the protein it detects
      beneath it in brackets
    * a Ponceau strip at the bottom, labelled on the right, as the loading
      control
    * a construct table above the strips, with lane identities as rows

The source image already carries rotated lane labels and a full-height ladder,
which is what this script replaces.

Every panel is drawn with aspect="equal" and its axes box is sized from the
crop's true pixel height, so no strip is ever stretched. An earlier version
used aspect="equal", which fills the box regardless of the crop's shape and
therefore changed the proportions of the blot. Rescaling one axis of a gel
image is not a presentational choice, it is a misrepresentation, and the
height_ratios below exist to make it impossible.

Geometry is hand-set below rather than detected. Ladder positions are read off
the printed colour standard by eye, which is how the lab does it. Lane centres
ARE measured, from the evenly spaced sample lanes in the Ponceau strip, because
the gel has dead space to the right of the last lane and assuming an even split
of the cropped width pushed every label to the right of its lane.

On Ponceau brightness. Rossner and Yamada (J Cell Biol 2004) allow adjustments
of brightness and contrast provided they are applied to the whole image and do
not obscure or eliminate information present in the original, and require that
non-linear adjustments such as gamma are disclosed in the legend. PONCEAU_GAIN
below is a single linear scale factor applied to the entire Ponceau panel and
nothing else, so it is within that rule. It is set to 1.0 by default. Anything
above 1.0 should be stated in the figure legend.

Usage:
    make_blot_options.py [--raw PATH] [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED = "#222222", "#666666"

# --- geometry, in the NATIVE source images --------------------------------
# western_blot_raw.png is a composite that was itself upsampled from these two
# originals, so reading it meant resampling twice. blot.jpg is 1248x1248 and
# ponceau.png is 1500x2000, and everything below is measured on them directly.
# That is the most real detail available. The lane region is only ~670 px wide
# in the original, so at thesis width the blot prints near 100 dpi however it
# is handled, and the only way past that is a rescan of the membrane.
BLOT_SRC = "blot.jpg"
PONCEAU_SRC = "ponceau.png"

# Lane grid, from make_blot_figure.py, which derived it from the slide the
# crops came from. Pitch checked against the receptor bands: 57.3 px matches
# the measured 53-59 px spacings.
# Pitch transferred from the composite, where the eight crisp receptor bands
# gave a least-squares fit with residuals under 2.6 px, then mapped back into
# blot.jpg coordinates. The membrane was photographed at a slight angle, so the
# true pitch is a little wider than the 57.3 the earlier grid used and the
# error accumulates towards the right-hand lanes.
LANE0, LANE_PITCH, N_LANES = 303.1, 57.57, 12
# Sample lanes with a margin either side, matching the framing the composite
# had. The coloured ladder runs to x=274, so the left edge sits just clear of
# it, and the right edge matches BLOT_CROP in make_blot_figure.py. Cropping
# tight to the outermost lane centres, as an earlier version did, cut the blot
# off flush against lanes 1 and 12.
LANE_X0, LANE_X1 = 278, 1000

# Ladder band rows, detected from the coloured marker strip of blot.jpg.
LADDER_Y = {180: 369, 130: 396, 100: 448, 70: 491, 55: 554,
            40: 619, 35: 683, 25: 735, 15: 824, 10: 885}

# The two originals are framed differently. make_blot_figure.py butts
# blot.jpg[195:1000] against ponceau.png[276:1084], so that pair fixes the
# mapping between them and the Ponceau crop follows from the blot one.
BLOT_REF, PONCEAU_REF = (195, 1000), (276, 1084)
# Cropped to the loading band itself, which is darkest at row 927. The full
# 891-1027 window is mostly blank membrane below it, and at true aspect that
# made the Ponceau taller than both blots together. Maidment's figures show a
# thin strip, which is what this now matches.
PONCEAU_Y0, PONCEAU_Y1 = 850, 965

# Ponceau.jpeg in the run tree and data/ponceau.png are the same photograph
# (mean pixel difference 0.05), so these are read off data/ponceau.png. Same
# PageRuler standard as the blot, detected by hue in the ladder strip at
# x 295-365: red 70 kDa at 870, with 100 at 830, 55 at 937, 40 at 1002,
# 35 at 1070 and 15 at 1222. The loading band row sits at about 927, so the
# blue 55 bar lies immediately beside it and the red 70 bar is the next one up.
# Those are the two marks used. 40 at 1002 was tried first and sits well below
# the bands, far enough to read as unrelated to them.
PONCEAU_LADDER_Y = {100: 830, 70: 870, 55: 937, 40: 1002, 35: 1070, 15: 1222}
PONCEAU_MARKS = [70, 55]

RECEPTOR_CROP = (352, 412)           # spans the 180 and 130 marks
# Centred on the control band, whose visible mass sits near row 815, with the
# 15 kDa mark at 824 falling just below centre. Extended up to row 722 so the
# 25 kDa mark at 735 is inside the crop as well: one tick alone gives the strip
# no scale.
EFFECTOR_CROP = (722, 855)
FULL_CROP = (352, 900)

# Linear brightness on the Ponceau panel only. 1.0 leaves it untouched.
PONCEAU_GAIN = 1.0

LANES = ["d39_s0", "d17_s0", "d47_s0", "d39_s1", "d26_s1", "d18_s0",
         "d49_s0", "d40_s0", "d6_s1", "d14_s0", "Tie1 (-)", "OsHIPP43 (+)"]

# TO CONFIRM: the antibody is not in the lab record yet, so it is a placeholder
# drawn in grey. The two control lanes are taken from the 3 July expression
# sheet's mix order, which puts Tie1 at lane 11 and OsHIPP43 at lane 12. Tie1
# came from another group and is assumed to carry a different tag, which is why
# it gives no signal; OsHIPP43 is the ~18 kDa band. Provisional pending the
# construct records.
ANTIBODY = "$\\alpha$-FLAG"
DETECTS = "(resurfaced Pikp-1)"
UNCONFIRMED = True


def ponceau_panel_height():
    """Ponceau crop height expressed in blot pixels, so the stack is faithful."""
    w = ponceau_x(LANE_X1) - ponceau_x(LANE_X0)
    return (PONCEAU_Y1 - PONCEAU_Y0) * (LANE_X1 - LANE_X0) / w


def crop_height(crop_range):
    """Pixel height of a crop, used to size its axes box faithfully."""
    return crop_range[1] - crop_range[0]


def ponceau_x(x):
    """Map an x in blot.jpg onto the matching x in ponceau.png."""
    b0, b1 = BLOT_REF
    p0, p1 = PONCEAU_REF
    return p0 + (x - b0) / (b1 - b0) * (p1 - p0)


# Inches reserved around the stack for the labels that sit outside it.
PAD_LEFT = 0.40        # kDa numbers and their ticks
PAD_RIGHT = 1.55       # antibody and Ponceau labels
PAD_TOP = 0.85         # rotated lane names
PAD_BOTTOM = 0.10
GAP = 0.06             # inches between strips


def stack_axes(fig_w, heights_px, gap=GAP, pad_top=PAD_TOP):
    """Create a figure and a column of axes whose boxes match their crops.

    subplots + height_ratios cannot do this. The axes width follows the figure
    width, so with aspect="equal" matplotlib letterboxes each image inside a
    box of the wrong shape and the leftovers read as gaps between the strips.
    Placing the axes explicitly makes every box exactly as tall as its crop is
    wide, so the image fills it and the only space between strips is `gap`.
    """
    panel_w = fig_w - PAD_LEFT - PAD_RIGHT
    scale = panel_w / (LANE_X1 - LANE_X0)          # inches per source pixel
    heights_in = [h * scale for h in heights_px]
    fig_h = pad_top + sum(heights_in) + gap * (len(heights_in) - 1) + PAD_BOTTOM
    fig = plt.figure(figsize=(fig_w, fig_h))
    axes, y = [], fig_h - pad_top
    for h in heights_in:
        y -= h
        axes.append(fig.add_axes(
            [PAD_LEFT / fig_w, y / fig_h, panel_w / fig_w, h / fig_h]))
        y -= gap
    return fig, axes


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "text.color": INK, "axes.grid": False,
    })


def crop(img, y0, y1):
    return img.crop((LANE_X0, y0, LANE_X1, y1))


def ponceau(_img=None):
    """The Ponceau panel from its own original, optionally brightened."""
    im = Image.open(os.path.join(HERE, "data", PONCEAU_SRC)).convert("RGB")
    box = (int(round(ponceau_x(LANE_X0))), PONCEAU_Y0,
           int(round(ponceau_x(LANE_X1))), PONCEAU_Y1)
    a = np.asarray(im.crop(box)).astype(float)
    if PONCEAU_GAIN != 1.0:
        a = np.clip(a * PONCEAU_GAIN, 0, 255)
    return a.astype(np.uint8)


def draw_ponceau_marks(ax):
    """kDa marks on the Ponceau strip, in the same style as draw_strip.

    The Ponceau has its own ladder coordinates because it is a different
    photograph from blot.jpg, so LADDER_Y does not apply to it.
    """
    h = PONCEAU_Y1 - PONCEAU_Y0
    for kda in PONCEAU_MARKS:
        y = PONCEAU_LADDER_Y[kda] - PONCEAU_Y0
        if not (0 <= y <= h):
            continue
        ax.annotate(f"{kda}", xy=(0, y), xycoords=("axes fraction", "data"),
                    xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=8.5, color=INK)
        ax.plot([-0.004, 0], [y, y], transform=ax.get_yaxis_transform(),
                color=INK, lw=1.0, clip_on=False)


def lane_centres(n_lanes=None):
    """Lane centres in axes fraction of the cropped panel.

    Measured, not assumed. The gel carries dead space to the right of the last
    lane, so splitting the cropped width evenly puts every label right of its
    lane and the error grows across the gel.
    """
    n = n_lanes or len(LANES)
    xs = LANE0 + LANE_PITCH * np.arange(n)
    return (xs - LANE_X0) / (LANE_X1 - LANE_X0)


def draw_strip(ax, img, y0, y1, marks, label, sublabel=None, grey=False):
    ax.imshow(np.asarray(crop(img, y0, y1)), aspect="equal", interpolation="lanczos")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor(INK)
        spine.set_linewidth(0.8)
    h = y1 - y0
    for kda in marks:
        y = LADDER_Y[kda] - y0
        if not (0 <= y <= h):
            continue
        ax.annotate(f"{kda}", xy=(0, y), xycoords=("axes fraction", "data"),
                    xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=8.5, color=INK)
        ax.plot([-0.004, 0], [y, y], transform=ax.get_yaxis_transform(),
                color=INK, lw=1.0, clip_on=False)
    colour = MUTED if (grey and UNCONFIRMED) else INK
    ax.annotate(label, xy=(1, 0.62), xycoords="axes fraction",
                xytext=(8, 0), textcoords="offset points",
                ha="left", va="center", fontsize=9, color=colour)
    if sublabel:
        ax.annotate(sublabel, xy=(1, 0.32), xycoords="axes fraction",
                    xytext=(8, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8, color=colour)


def draw_lane_labels(ax, rotation=90, fontsize=8.5):
    for x, name in zip(lane_centres(), LANES):
        ax.annotate(name, xy=(x, 1), xycoords="axes fraction",
                    xytext=(0, 5), textcoords="offset points",
                    ha="center" if rotation else "center",
                    va="bottom", rotation=rotation, fontsize=fontsize,
                    color=INK)


def draw_construct_table(fig, ax, rows):
    """A +/- table above the strips, as in Maidment figure 5.3C."""
    xs = lane_centres()
    for i, (name, marks) in enumerate(rows):
        y = 1.06 + 0.075 * (len(rows) - 1 - i)
        ax.annotate(name, xy=(0, y), xycoords=("axes fraction", "axes fraction"),
                    xytext=(-8, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=8.5, color=INK)
        for x, m in zip(xs, marks):
            ax.annotate(m, xy=(x, y), xycoords="axes fraction",
                        ha="center", va="center", fontsize=9, color=INK)


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    # bbox_inches="tight" would re-crop the explicit geometry, so it is off.
    fig.savefig(path, facecolor="white", bbox_inches=None)
    plt.close(fig)
    print("wrote", path)


def option1_two_strips(img, outdir):
    """Receptor strip, effector strip, Ponceau. The house layout."""
    fig, axes = stack_axes(8.4, [crop_height(RECEPTOR_CROP),
                                 crop_height(EFFECTOR_CROP),
                                 ponceau_panel_height()])
    draw_strip(axes[0], img, *RECEPTOR_CROP, marks=[180, 130],
               label=ANTIBODY, sublabel=DETECTS, grey=True)
    draw_strip(axes[1], img, *EFFECTOR_CROP, marks=[25, 15],
               label=ANTIBODY, sublabel="(control lanes)", grey=True)
    axes[2].imshow(ponceau(img), aspect="equal", interpolation="lanczos")
    axes[2].set_xticks([]); axes[2].set_yticks([])
    for s in axes[2].spines.values():
        s.set_edgecolor(INK); s.set_linewidth(0.8)
    draw_ponceau_marks(axes[2])
    axes[2].annotate("Ponceau", xy=(1, 0.5), xycoords="axes fraction",
                     xytext=(8, 0), textcoords="offset points",
                     ha="left", va="center", fontsize=9, color=INK)
    draw_lane_labels(axes[0])
    axes[0].annotate("kDa", xy=(0, 1), xycoords="axes fraction",
                     xytext=(-6, 4), textcoords="offset points",
                     ha="right", va="bottom", fontsize=8.5, weight="bold")
    save(fig, outdir, "option1_two_strips")


def option2_construct_table(img, outdir):
    """Receptor strip and Ponceau, with a +/- construct table above."""
    # Extra head room for the construct table above the lane names.
    fig, axes = stack_axes(8.4, [crop_height(RECEPTOR_CROP),
                                 ponceau_panel_height()], pad_top=1.55)
    draw_strip(axes[0], img, *RECEPTOR_CROP, marks=[180, 130],
               label=ANTIBODY, sublabel=DETECTS, grey=True)
    axes[1].imshow(ponceau(img), aspect="equal", interpolation="lanczos")
    axes[1].set_xticks([]); axes[1].set_yticks([])
    for s in axes[1].spines.values():
        s.set_edgecolor(INK); s.set_linewidth(0.8)
    axes[1].annotate("Ponceau", xy=(1, 0.5), xycoords="axes fraction",
                     xytext=(8, 0), textcoords="offset points",
                     ha="left", va="center", fontsize=9, color=INK)
    resurfaced = ["+"] * 10 + ["-", "-"]
    effector = ["+"] * 12
    draw_construct_table(fig, axes[0], [
        ("resurfaced Pikp-1:HF", resurfaced),
        ("myc:AVR-PikF", effector)])
    for x, name in zip(lane_centres(), LANES):
        axes[0].annotate(name, xy=(x, 1.22), xycoords="axes fraction",
                         ha="center", va="bottom", rotation=90, fontsize=8.5)
    axes[0].annotate("kDa", xy=(0, 1), xycoords="axes fraction",
                     xytext=(-6, 4), textcoords="offset points",
                     ha="right", va="bottom", fontsize=8.5, weight="bold")
    save(fig, outdir, "option2_construct_table")


def option3_receptor_only(img, outdir):
    """The most compact reading: only the band that carries the result."""
    fig, axes = stack_axes(8.4, [crop_height(RECEPTOR_CROP),
                                 ponceau_panel_height()])
    draw_strip(axes[0], img, *RECEPTOR_CROP, marks=[180, 130],
               label=ANTIBODY, sublabel=DETECTS, grey=True)
    axes[1].imshow(ponceau(img), aspect="equal", interpolation="lanczos")
    axes[1].set_xticks([]); axes[1].set_yticks([])
    for s in axes[1].spines.values():
        s.set_edgecolor(INK); s.set_linewidth(0.8)
    axes[1].annotate("Ponceau", xy=(1, 0.5), xycoords="axes fraction",
                     xytext=(8, 0), textcoords="offset points",
                     ha="left", va="center", fontsize=9, color=INK)
    draw_lane_labels(axes[0])
    axes[0].annotate("kDa", xy=(0, 1), xycoords="axes fraction",
                     xytext=(-6, 4), textcoords="offset points",
                     ha="right", va="bottom", fontsize=8.5, weight="bold")
    # d18_s0 is the design that did not express, so it is called out.
    x = lane_centres()[LANES.index("d18_s0")]
    axes[0].annotate("no band", xy=(x, 0.5), xycoords="axes fraction",
                     xytext=(0, -34), textcoords="offset points",
                     ha="center", fontsize=8, color="#D55E00",
                     arrowprops=dict(arrowstyle="->", color="#D55E00", lw=1.1))
    save(fig, outdir, "option3_receptor_only")


def option4_full_relabelled(img, outdir):
    """Minimal intervention: keep the full panel, replace the ladder lane."""
    fig, axes = stack_axes(8.4, [crop_height(FULL_CROP),
                                 ponceau_panel_height()], gap=0.03)
    draw_strip(axes[0], img, *FULL_CROP,
               marks=[180, 130, 100, 70, 55, 40, 35, 25, 15, 10],
               label=ANTIBODY, sublabel=DETECTS, grey=True)
    axes[1].imshow(ponceau(img), aspect="equal", interpolation="lanczos")
    axes[1].set_xticks([]); axes[1].set_yticks([])
    for s in axes[1].spines.values():
        s.set_edgecolor(INK); s.set_linewidth(0.8)
    axes[1].annotate("Ponceau", xy=(1, 0.5), xycoords="axes fraction",
                     xytext=(8, 0), textcoords="offset points",
                     ha="left", va="center", fontsize=9, color=INK)
    draw_lane_labels(axes[0])
    axes[0].annotate("kDa", xy=(0, 1), xycoords="axes fraction",
                     xytext=(-6, 4), textcoords="offset points",
                     ha="right", va="bottom", fontsize=8.5, weight="bold")
    save(fig, outdir, "option4_full_relabelled")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(HERE, "data", BLOT_SRC))
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    img = Image.open(args.raw).convert("RGB")
    option1_two_strips(img, args.outdir)
    option2_construct_table(img, args.outdir)
    option3_receptor_only(img, args.outdir)
    option4_full_relabelled(img, args.outdir)


if __name__ == "__main__":
    main()
