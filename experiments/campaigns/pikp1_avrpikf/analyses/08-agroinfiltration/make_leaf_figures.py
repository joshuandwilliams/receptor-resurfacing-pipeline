#!/usr/bin/env python3
"""Candidate leaf figures for the AVR-PikF agroinfiltration result.

Four options, following the two conventions used in the group's theses:

  * Maidment (2019, fig 5.1): one representative leaf, background removed,
    construct labels placed around the leaf next to their spots.
  * Stone (2025, fig 4-1): one representative leaf on the dark UV background,
    dashed dividers separating the infiltration sectors, labels outside.

Spot positions come from the annotated slide, mapped back through the slide's
crop rectangle onto full-resolution pixels, so the labels sit where the
infiltrations actually were rather than being placed by eye.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage  # noqa: E402

HERE = Path(__file__).parent
DATA = HERE / "data"
# The UV leaf photographs are lab record, not analysis output, so they are
# read from the run tree. leaf_uv.tif here was a renamed copy of DSC_9957.
AGRO = HERE / ".." / ".." / "runs" / "crystal_full_test_contig" / \
    "agroinfiltrations" / "7dpi_5conditions"
LEAF_UV = AGRO / "DSC_9957.tif"
PLOTS = HERE / "plots"

INK = "#222222"
LABEL_SIZE = 8

# The slide's five annotation groups, renamed for the thesis. "Design N" on the
# slide is the AGROINFILTRATION id, not the RFDiffusion design number: see
# data/design_id_map.csv. So Design 5 is d26_s1 and Design 6 is d18_s0 (d6_s1 is
# agroinfiltration 009 and is not on this leaf).
RENAME = {
    "ControlsPik-HIPP43 : AVR-PikFPik-HIPP43 : PWL2": "control",
    "Design 5 : AVR-PikF": "d26_s1\nAVR-PikF",
    "Design 5 : PWL2": "d26_s1\nPWL2",
    "Design 6 : AVR-PikF": "d18_s0\nAVR-PikF",
    "Design 6 : PWL2": "d18_s0\nPWL2",
}
CONTROL_NAMES = ["Pik-HIPP43\nPWL2", "Pik-HIPP43\nAVR-PikF"]


def load_spots():
    return json.loads((DATA / "spots.json").read_text())


def leaf_clusters(spots, k=6):
    """Group the 36 spots into leaves by simple 2x3 gridding of their centres."""
    xs = np.array([s["x"] for s in spots], float)
    ys = np.array([s["y"] for s in spots], float)
    col = np.digitize(xs, [xs.min() + (xs.max() - xs.min()) * f for f in (1 / 3, 2 / 3)])
    row = (ys > (ys.min() + ys.max()) / 2).astype(int)
    groups = {}
    for i, s in enumerate(spots):
        groups.setdefault((row[i], col[i]), []).append(s)
    return groups


def leaf_mask_full(arr, spots, thresh=55, open_px=7, close_px=15):
    """Mask of ONE leaf, as a single filled connected component.

    A plain luminance cut fails twice over. It keeps every speck of debris on the
    tray, which is what produced the green grain around the earlier cut-outs, and
    it drops the leaf's darkest tissue. The fix is ordered: threshold above the
    background grain (the image's median luminance is 17 and the leaves are well
    over 150), open to delete isolated specks, keep only the component the
    infiltration spots sit in, then close and hole-fill so dark veins and edge
    tissue inside that component come back.

    Closing BEFORE labelling is what broke the first attempt: a large structuring
    element bridged the debris field into one blob spanning the whole tray, so
    the "leaf" component became the entire background.
    """
    lum = arr.astype(np.float32).max(axis=2)
    m = ndimage.binary_opening(lum > thresh, structure=np.ones((open_px, open_px)))
    lab, n = ndimage.label(m)
    if not n:
        return np.zeros(lum.shape, bool)
    ids = [lab[int(s["y"]), int(s["x"])] for s in spots]
    ids = [i for i in ids if i > 0]
    if not ids:
        return np.zeros(lum.shape, bool)
    keep = max(set(ids), key=ids.count)
    leaf = ndimage.binary_closing(lab == keep, structure=np.ones((close_px, close_px)))
    return ndimage.binary_fill_holes(leaf)


def mask_bbox(mask, pad=40):
    ys, xs = np.nonzero(mask)
    return [max(int(xs.min()) - pad, 0), max(int(ys.min()) - pad, 0),
            int(xs.max()) + pad, int(ys.max()) + pad]


def crop_leaf(arr, spots, pad=40):
    """Bounding box of the leaf itself, not of its infiltration spots.

    Padding the spots was what clipped the leaf tips, since the topmost spot sits
    well below the top of the leaf.
    """
    return mask_bbox(leaf_mask_full(arr, spots), pad)


def to_white_bg(arr, box, mask):
    """Crop, then replace everything outside the leaf with white."""
    sub = arr[box[1]:box[3], box[0]:box[2]].copy()
    m = mask[box[1]:box[3], box[0]:box[2]]
    out = np.full_like(sub, 255)
    out[m] = sub[m]
    return out


def spot_labels(spots):
    """(x, y, text) for each spot, numbering the two controls per leaf."""
    out, n_ctrl = [], 0
    for s in sorted(spots, key=lambda s: (s["y"], s["x"])):
        name = RENAME[s["label"]]
        if name == "control":
            name = CONTROL_NAMES[n_ctrl % 2]
            n_ctrl += 1
        out.append((s["x"], s["y"], name))
    return out


def draw_labelled(ax, img, box, spots, white_bg):
    ax.imshow(img)
    ax.set_xlim(0, box[2] - box[0])
    ax.set_ylim(box[3] - box[1], 0)
    cx = (box[2] - box[0]) / 2
    for x, y, text in spot_labels(spots):
        lx, ly = x - box[0], y - box[1]
        side = -1 if lx < cx else 1
        tx = lx + side * (0.30 * (box[2] - box[0]))
        colour = INK if white_bg else "white"
        ax.annotate(text, xy=(lx, ly), xytext=(tx, ly),
                    ha="right" if side < 0 else "left", va="center",
                    fontsize=LABEL_SIZE, color=colour,
                    arrowprops=dict(arrowstyle="-", color=colour, lw=0.7,
                                    shrinkA=2, shrinkB=8))
    ax.axis("off")


def main():
    PLOTS.mkdir(exist_ok=True)
    spots = load_spots()
    arr = np.asarray(Image.open(LEAF_UV).convert("RGB"))
    groups = leaf_clusters(spots)

    strongest = groups[(1, 1)]      # bottom middle
    weakest = groups[(1, 0)]        # bottom left

    masks = {k: leaf_mask_full(arr, g) for k, g in groups.items()}

    # ── Option 1: Maidment style, background removed ────────────────────────
    m1 = masks[(1, 1)]
    box = mask_bbox(m1)
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    draw_labelled(ax, to_white_bg(arr, box, m1), box, strongest, True)
    fig.savefig(PLOTS / "option1_single_white.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # ── Option 2: Stone style, dark UV background kept ──────────────────────
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    draw_labelled(ax, arr[box[1]:box[3], box[0]:box[2]], box, strongest, False)
    fig.savefig(PLOTS / "option2_single_dark.png", dpi=300, bbox_inches="tight",
                facecolor="black")
    plt.close(fig)

    # ── Option 3: strongest against weakest, labelled once ──────────────────
    m2 = masks[(1, 0)]
    b2 = mask_bbox(m2)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.2))
    draw_labelled(axes[0], to_white_bg(arr, box, m1), box, strongest, True)
    axes[1].imshow(to_white_bg(arr, b2, m2))
    axes[1].axis("off")
    axes[0].set_title("Strongest response", fontsize=10, color=INK)
    axes[1].set_title("Weakest response", fontsize=10, color=INK)
    fig.savefig(PLOTS / "option3_strong_vs_weak.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    # ── Option 4: all six leaves, one labelled ──────────────────────────────
    order = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 7.0))
    for ax, key in zip(axes.ravel(), order):
        g = groups[key]
        mm = masks[key]
        bb = mask_bbox(mm)
        ax.imshow(to_white_bg(arr, bb, mm))
        ax.axis("off")
        for x, y, _ in spot_labels(g):
            ax.plot(x - bb[0], y - bb[1], "o", ms=7, mfc="none", mec=INK, mew=0.8)
    fig.tight_layout()
    fig.savefig(PLOTS / "option4_all_leaves.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)

    for p in sorted(PLOTS.glob("option*.png")):
        print(f"  plots/{p.name}")


if __name__ == "__main__":
    main()
