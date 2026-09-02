#!/usr/bin/env python3
"""Crop the contact close-up out of the full-view render.

recapture_contact.cxc renders panel_d_contact_full.png at the SAME camera as
the overview panels, with nothing hidden. The close-up is then just a crop of
that image, which is what keeps the rest of the complex visible behind the
contact.

Earlier attempts re-aimed the ChimeraX camera and hid the distant cartoon
instead. Both looked wrong: the framing pulled the complex out of the frame and
the hiding left white space where the structure should have been. A plain 2D
crop of the full view has neither problem.

The box is a fraction of the render so it survives a change of output size.
Tuned by eye on the 4800 x 4560 render: the Asp42-Arg8 contact and both labels
sit around x 2300-3400, y 1800-3600, and the box below keeps the effector sheet
on the left and the design region 1 loop on the right for context.

Usage:
    crop_contact.py [--frac x0 y0 x1 y1]
"""
import argparse
import os

from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "renders", "panel_d_contact_full.png")
OUT = os.path.join(HERE, "renders", "panel_d_contact.png")

# x0, y0, x1, y1 as fractions of the full render.
FRAC = (0.32, 0.29, 0.86, 0.86)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frac", nargs=4, type=float, default=None,
                    help="x0 y0 x1 y1 as fractions of the full render")
    args = ap.parse_args()
    frac = tuple(args.frac) if args.frac else FRAC

    im = Image.open(SRC)
    im.load()
    w, h = im.size
    box = (int(frac[0] * w), int(frac[1] * h),
           int(frac[2] * w), int(frac[3] * h))
    im.crop(box).save(OUT)
    print(f"{SRC} {w}x{h}")
    print(f"  crop {box} -> {OUT} "
          f"{box[2] - box[0]}x{box[3] - box[1]}")


if __name__ == "__main__":
    main()
