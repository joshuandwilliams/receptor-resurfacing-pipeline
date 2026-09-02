#!/usr/bin/env python3
"""Compose the RFdiffusion design-region figure from four structure renders.

Layout is a large input panel with the colour key beside it, above a row of
three design panels.  Panels keep their source aspect and relative scale rather
than being normalised to a common height, because the point of the figure is
that the rebuilt design region differs in length between designs.  Panels are
top-aligned so that difference reads downwards from a common start.

Source images are the slide renders, which are only a few hundred pixels wide,
so the composite is capped at their native resolution rather than upscaled.

Usage:
    build_design_region_figure.py --input in.png --design 17 d17.png \\
        --design 39 d39.png --design 26 d26.png --out fig.png
"""

import argparse

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

INK = "#222222"          # matches the analysis figures' text colour
MUTED = "#666666"        # MUTED there; used for the qualifying sublabel
LABEL_SIZE = 10          # axes.titlesize in confidence_common.apply_style
LEGEND_SIZE = 9          # font.size there
SUBLABEL_SIZE = 8
SUBLABEL_RISE = 0.15     # inches; how far the bold label sits above the sublabel

# Sampled from the renders themselves so the key cannot drift from the cartoons.
KEY = [("#FF0000", "Design region 1"),
       ("#940094", "Design region 2"),
       ("#8FBDFF", "Fixed Pikp-1 HMA backbone"),
       ("#A7FFA7", "AVR-PikF effector")]

MARGIN = 0.06            # inches
LABEL_PAD = 0.30         # inches reserved above each panel row for its label
ROW_GAP = 0.20
COL_GAP = 0.18
INPUT_SCALE = 1.45       # input panel drawn larger than the design panels
LEGEND_MIN_W = 2.15      # inches; enough for the longest key entry

# Save at print resolution rather than at the sources' own scale. The renders
# are low-resolution whatever we do, but the text and the key are vector-drawn
# here and there is no reason to penalise them for it.
DPI_OUT = 300


def load(path, w_in, h_in):
    """Read a render, crop its transparent border, and keep its physical size.

    The renders carry a lot of empty margin, which would otherwise become
    whitespace inside the figure. Cropping shrinks the drawn box by exactly the
    fraction cropped, so the panels stay on a common scale and the design region
    lengths remain comparable between them.
    """
    im = mpimg.imread(path)
    h, w = im.shape[:2]
    if im.shape[2] == 4:
        ys, xs = np.nonzero(im[:, :, 3] > 0.02)
        y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    else:
        y0, y1, x0, x1 = 0, h, 0, w
    im = im[y0:y1, x0:x1]
    # Flatten onto white before it reaches imshow. A transparent render carries
    # arbitrary RGB in its fully-transparent pixels, and lanczos resampling
    # interpolates the colour channels without reference to alpha, so that junk
    # is dragged into the visible edges as white and dark speckle. Compositing
    # here makes the alpha channel irrelevant to the interpolation. Opaque
    # sources are unaffected.
    if im.shape[2] == 4:
        a = im[:, :, 3:4]
        im = im[:, :, :3] * a + (1.0 - a)
    return im, w_in * (x1 - x0) / w, h_in * (y1 - y0) / h


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", nargs=3, metavar=("PNG", "W_IN", "H_IN"), required=True)
    ap.add_argument("--design", nargs=5, action="append",
                    metavar=("LABEL", "SUBLABEL", "PNG", "W_IN", "H_IN"), required=True)
    ap.add_argument("--out-a", required=True, help="panel A png (input + key)")
    ap.add_argument("--out-b", required=True, help="panel B png (design row)")
    args = ap.parse_args()

    inp, inp_w, inp_h = load(args.input[0], float(args.input[1]), float(args.input[2]))
    inp_w, inp_h = inp_w * INPUT_SCALE, inp_h * INPUT_SCALE
    # Sizes come from the source layout, NOT from a common height. Forcing the
    # designs to one height would normalise away the very difference being shown.
    designs = [(lab, sub, *load(png, float(w), float(h)))
               for lab, sub, png, w, h in args.design]

    d_w = [w for _, _, _, w, _ in designs]
    d_h = max(h for _, _, _, _, h in designs)
    row_w = sum(d_w) + COL_GAP * (len(designs) - 1)
    top_w = inp_w + COL_GAP + max(LEGEND_MIN_W, 0.0)
    legend_w = LEGEND_MIN_W

    content_w = max(row_w, top_w)

    # One figure per panel. The A/B letters are added by assemble_figure.py in the
    # thesis SVG template, at the template's own font and position, rather than
    # being drawn into the raster here.
    def new_fig(width, height):
        return plt.figure(figsize=(width, height), dpi=DPI_OUT)

    def place(fig, fw, fh, im, x, y_top, w, h, label=None, sublabel=None):
        """x and y_top in inches from the figure's top-left.

        The label sits in the band reserved by LABEL_PAD above the panel. It is
        bold, since it names the design, and the sublabel beneath it is smaller
        and lighter because it qualifies the label rather than competing with it.
        Panel A passes neither, so its band collapses to nothing.
        """
        ax = fig.add_axes([x / fw, 1.0 - (y_top + h) / fh, w / fw, h / fh])
        ax.imshow(im, interpolation="lanczos")
        ax.axis("off")
        if sublabel:
            fig.text((x + w / 2) / fw, 1.0 - (y_top - 0.06) / fh, sublabel,
                     ha="center", va="bottom", fontsize=SUBLABEL_SIZE, color=MUTED)
        if label:
            y = (y_top - 0.06) - (SUBLABEL_RISE if sublabel else 0.0)
            fig.text((x + w / 2) / fw, 1.0 - y / fh, label, ha="center",
                     va="bottom", fontsize=LABEL_SIZE, color=INK, fontweight="bold")

    # ── Panel A: input structure and colour key ─────────────────────────────
    # No label here: the caption names the panel, and the colour key already
    # says what the regions are.
    fa_w = content_w + 2 * MARGIN
    fa_h = inp_h + 2 * MARGIN
    figa = new_fig(fa_w, fa_h)
    top_x = MARGIN + (content_w - top_w) / 2.0
    place(figa, fa_w, fa_h, inp, top_x, MARGIN, inp_w, inp_h)
    lax = figa.add_axes([(top_x + inp_w + COL_GAP) / fa_w,
                         1.0 - (MARGIN + inp_h) / fa_h,
                         legend_w / fa_w, inp_h / fa_h])
    lax.axis("off")
    lax.legend(handles=[Patch(facecolor=c, edgecolor="#333333", linewidth=0.6, label=t)
                        for c, t in KEY],
               loc="center left", frameon=False, fontsize=LEGEND_SIZE,
               labelcolor=INK, handlelength=1.4, handleheight=1.0,
               borderpad=0.0, labelspacing=0.9)
    figa.savefig(args.out_a, dpi=DPI_OUT, facecolor="white")

    # ── Panel B: the three designs, top-aligned at a common scale ───────────
    fb_w = content_w + 2 * MARGIN
    fb_h = d_h + LABEL_PAD + 2 * MARGIN
    figb = new_fig(fb_w, fb_h)
    x = MARGIN + (content_w - row_w) / 2.0
    for lab, sub, im, w, h in designs:
        place(figb, fb_w, fb_h, im, x, MARGIN + LABEL_PAD, w, h,
              label=lab, sublabel=sub)
        x += w + COL_GAP
    figb.savefig(args.out_b, dpi=DPI_OUT, facecolor="white")

    native = min([inp.shape[1] / inp_w] + [im.shape[1] / w for _, _, im, w, _ in designs])
    print(f"{args.out_a}  {fa_w * DPI_OUT:.0f}x{fa_h * DPI_OUT:.0f}px\n"
          f"{args.out_b}  {fb_w * DPI_OUT:.0f}x{fb_h * DPI_OUT:.0f}px\n"
          f"(sources are only {native:.0f} px/in, so the renders are interpolated)")


if __name__ == "__main__":
    main()
