#!/usr/bin/env python3
"""PROVISIONAL western blot figure, in the group's thesis house style.

Layout follows Maidment (2019, fig 5.1C) and Stone (2025, fig 4-5): lane
identities above the membrane, kDa ladder values on the left aligned to the
detected marker bands, antibody on the right, Ponceau loading control beneath.

THREE THINGS ARE UNCONFIRMED and are drawn as placeholders:
  * the ladder product. Ten bands with a single red band fourth from the top and
    a green terminal band match Thermo PageRuler Prestained (26616), but the
    group's methods list three ladders in use, so this is inference from the band
    pattern rather than a record.
  * the antibody this membrane was probed with.
  * what the two control lanes contained.

Lane identities come from data/design_id_map.csv, on the assumption that lanes
1-10 run in agroinfiltration order.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage  # noqa: E402

HERE = Path(__file__).parent
DATA = HERE / "data"
PLOTS = HERE / "plots"

INK = "#222222"
# Thermo PageRuler Prestained (26616). PROVISIONAL — confirm against the lab record.
LADDER_KDA = [180, 130, 100, 70, 55, 40, 35, 25, 15, 10]

# Both crops are the ones used on slide 7 of AI_Pikp1_HMA_Prelim_2_JoshW.pptx,
# converted from its srcRect fractions. Reusing them keeps the blot and the
# Ponceau strip in the same lane-to-lane registration the slide already had.
BLOT_CROP = (195, 323, 1000, 965)       # image7.jpg, 1248x1248
PONCEAU_CROP = (276, 891, 1084, 1027)   # image6.png, 1500x2000; band row only
DPI_OUT = 300
LADDER_X = (200, 285)
LANE0, LANE_PITCH, N_LANES = 302, 57.3, 12


def ladder_rows(arr):
    """y of each coloured marker band, from saturation in the ladder strip."""
    strip = arr[:, LADDER_X[0]:LADDER_X[1]].astype(int)
    sat = strip.max(2) - strip.min(2)
    lab, n = ndimage.label((sat > 40).sum(1) > 8)
    return [float(np.nonzero(lab == i)[0].mean()) for i in range(1, n + 1)]


def lane_names():
    rows = list(csv.DictReader((DATA / "design_id_map.csv").open()))
    names = [r["design"] for r in rows]
    # Lanes 11 and 12 from the 3 July expression sheet's mix order. Tie1 came
    # from another group and is assumed to carry a different tag, which is why
    # it gives no alpha-FLAG signal; OsHIPP43 is the ~18 kDa band. Both
    # identities are provisional pending the construct records.
    return names + ["Tie1 (-)", "OsHIPP43 (+)"]


def main():
    PLOTS.mkdir(exist_ok=True)
    blot = np.asarray(Image.open(DATA / "blot.jpg").convert("RGB"))
    ponceau = np.asarray(Image.open(DATA / "ponceau.png").convert("RGB"))

    rows = ladder_rows(blot)
    x0, y0, x1, y1 = BLOT_CROP
    sub = blot[y0:y1, x0:x1]

    px0, py0, px1, py1 = PONCEAU_CROP
    pon = ponceau[py0:py1, px0:px1]

    # Both strips are drawn at exactly the same width and butted together. The
    # axes heights are derived from each image's own aspect, because imshow
    # defaults to aspect="equal" and would otherwise shrink the wide, short
    # Ponceau strip to fit its box instead of filling it.
    fig_w = 7.2
    left, right = 0.16, 0.03
    strip_w = fig_w * (1 - left - right)
    blot_h = strip_w * sub.shape[0] / sub.shape[1]
    pon_h = strip_w * pon.shape[0] / pon.shape[1]
    top_pad, bot_pad = 0.62, 0.10
    fig_h = top_pad + blot_h + pon_h + bot_pad
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI_OUT)

    ax = fig.add_axes([left, (bot_pad + pon_h) / fig_h,
                       1 - left - right, blot_h / fig_h])
    ax.imshow(sub, aspect="auto")
    ax.axis("off")

    names = lane_names()
    for i, nm in enumerate(names[:N_LANES]):
        lx = LANE0 + i * LANE_PITCH - x0
        if 0 <= lx <= x1 - x0:
            ax.text(lx, -12, nm, rotation=90, ha="center", va="bottom",
                    fontsize=7, color=INK)

    for kda, yy in zip(LADDER_KDA, rows):
        ly = yy - y0
        if 0 <= ly <= y1 - y0:
            ax.text(-14, ly, f"{kda}", ha="right", va="center", fontsize=7, color=INK)
    ax.text(-14, -12, "kDa", ha="right", va="bottom", fontsize=7, color=INK)

    pax = fig.add_axes([left, bot_pad / fig_h, 1 - left - right, pon_h / fig_h])
    pax.imshow(pon, aspect="auto")
    pax.axis("off")
    pax.text(-0.01, 0.5, "Ponceau", transform=pax.transAxes, ha="right",
             va="center", fontsize=8, color=INK)

    out = PLOTS / "western_blot_draft.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"{out}  ladder bands at y={[round(r) for r in rows]}")


if __name__ == "__main__":
    main()
