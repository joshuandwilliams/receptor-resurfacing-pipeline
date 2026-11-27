#!/usr/bin/env python3
"""Draw the ReceptorResurfacing hex logo and favicon.

Matches the hex sticker frame of the author's three other tool repos, so the
four sit together: a pointy-top hexagon with a corner at the top and the
bottom, rounded corners, a green rim, the package name arched over the top in
white, and one motif below it.

The frame is not re-derived here. It is the geometry worked out for
structure-negative-steering's logo, measured off CDALabeller rather than
eyeballed: a hexagon 208 x 240 inside a 240 px square, which is a pointy-top
hexagon of circumradius 120 (width = 120 * sqrt(3) = 207.8), with corner arcs
of r = 20 px at R = 120, an 8 px rim, and #4c6433 as the rim colour, which is
the most common rim pixel in both reference stickers.

Three things about that frame are easy to get wrong and are handled below.
A stroke is centred on its path, so the polygon is inset by half the rim width
or the sticker comes out 226 px wide. Stroke widths are in points, which are
absolute, so they go through _pt() rather than being scaled by BASE_PX/240.
And rounding pulls the top and bottom vertices in but not the left and right
extremes, because those are flat edges, so the result is normalised to its
bounding box.

The motif is what this pipeline does, and it is deliberately the inverse of
NegativeSteering's. There, an effector moves OFF a surface it should not be on.
Here a receptor's surface is rebuilt so an effector can sit ON it: the receptor
holds a pocket carved into its upper face, the rim of that pocket is drawn in
white because it is the part the pipeline redesigns, white ticks mark the
designed positions, and the effector is seated in it.

    python3 docs/make_logo.py

Writes docs/logo.png (240 px, for the README) and docs/favicon.ico (a 16/32/
48/64/128/256 px stack, so a browser picks its own frame rather than scaling
the 240 down).
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402

HERE = Path(__file__).resolve().parent

# The shared palette. Taken from structure-negative-steering, which sampled it
# from CDALabeller and AutoCDAScorer rather than choosing it again.
INK = "#000000"       # hexagon fill, and what carves the pocket
RIM = "#4c6433"       # hexagon border
RECEPTOR = "#d8d677"  # the receptor being resurfaced
EFFECTOR = "#e28060"  # the effector it is being built to bind
WHITE = "#fbfbfb"     # the redesigned surface, the ticks, and the name

# Axis units: the canvas spans (-1, 1) on both axes and maps to 240 px, so one
# unit is 120 px. Every size below is written in px-at-240 and converted, so
# the numbers can be compared directly with a measurement off a reference PNG.
PX = 1.0 / 120.0                  # one px at 240, in axis units

RIM_W = 8.0                       # rim width in px, matching the references
TARGET_CORNER = 20.0              # outer corner radius in px, fitted to them

# Matplotlib centres a stroke on its path, so the drawn shape runs half the rim
# width wider than the polygon on every side. The polygon is therefore inset by
# that half-width, and its corner radius reduced by it, so the OUTER edge of
# the rim lands on 208 x 240 with 20 px corners.
_HALF = RIM_W / 2.0 * PX
HEX_CORNER = TARGET_CORNER * PX - _HALF

# Rounding a corner pulls the outline in along the bisector by
# r * (1/sin(60) - 1), so the hexagon is scaled up by exactly that before being
# drawn and the result still fills its bounding box.
_PULL_IN = 1.0 / math.sin(math.radians(60)) - 1.0
HEX_R = (1.0 - _HALF) + _PULL_IN * HEX_CORNER

OUT_PX = 240
HEX_W_PX = 208


def rounded_hexagon(R: float, r: float, per_arc: int = 48) -> np.ndarray:
    """A pointy-top hexagon whose corners are circular arcs of radius r.

    Pointy-top means vertices at 90 and 270 degrees, so a corner sits at the
    top and the bottom and the flat edges run down the left and right sides.
    """
    verts = [(R * math.cos(math.radians(a)), R * math.sin(math.radians(a)))
             for a in (90, 30, -30, -90, -150, 150)]
    pts: list[tuple[float, float]] = []
    n = len(verts)
    for i in range(n):
        p0, p1, p2 = verts[i - 1], verts[i], verts[(i + 1) % n]
        v1 = np.array(p0) - np.array(p1)
        v1 /= np.linalg.norm(v1)
        v2 = np.array(p2) - np.array(p1)
        v2 /= np.linalg.norm(v2)
        ang = math.acos(float(np.clip(v1 @ v2, -1.0, 1.0)))
        tangent = r / math.tan(ang / 2.0)
        t1 = np.array(p1) + v1 * tangent
        t2 = np.array(p1) + v2 * tangent
        bisector = v1 + v2
        bisector /= np.linalg.norm(bisector)
        centre = np.array(p1) + bisector * (r / math.sin(ang / 2.0))
        a1 = math.atan2(t1[1] - centre[1], t1[0] - centre[0])
        a2 = math.atan2(t2[1] - centre[1], t2[0] - centre[0])
        sweep = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
        for s in np.linspace(0.0, 1.0, per_arc):
            th = a1 + sweep * s
            pts.append((centre[0] + r * math.cos(th),
                        centre[1] + r * math.sin(th)))
    return np.asarray(pts)


def arched_text(ax, text: str, radius: float, y0: float, size: float,
                spread_deg: float) -> None:
    """Set `text` along a circular arc, each glyph rotated to the tangent.

    Matplotlib has no text-on-a-path, so each glyph is placed and rotated by
    hand. The angular step is proportional to each glyph's own advance width,
    measured from its outline, rather than one step per character. Even steps
    leave visible gaps around narrow letters, and this name has five of them.
    """
    fp = FontProperties(size=size, weight="bold")

    def advance(s: str) -> float:
        if s == " ":
            return TextPath((0, 0), "n", prop=fp).get_extents().width * 0.6
        return TextPath((0, 0), s, prop=fp).get_extents().width

    # Kerning is ignored, so a pair is measured together and the single-glyph
    # widths are scaled to match the whole string. Without this the arc runs
    # long.
    widths = [advance(c) for c in text]
    gap = advance("nn") - 2 * advance("n")
    total = sum(widths) + gap * (len(text) - 1)
    steps = [(w + gap) / total * spread_deg for w in widths]

    ang = 90.0 + spread_deg / 2.0
    for ch, step in zip(text, steps):
        centre = ang - step / 2.0
        rad = math.radians(centre)
        ax.text(radius * math.cos(rad), y0 + radius * math.sin(rad), ch,
                ha="center", va="center", color=WHITE, fontsize=size,
                fontweight="bold", rotation=centre - 90.0,
                rotation_mode="anchor", zorder=7)
        ang -= step


def blob(cx: float, cy: float, rx: float, ry: float, rot: float = 0.0,
         lobes: int = 5, wobble: float = 0.05, phase: float = 0.0) -> np.ndarray:
    """A closed rounded outline, as a stand-in for a folded domain.

    A plain ellipse reads as a ball, and a cartoon backbone trace turns to mush
    at 16 px. This is an ellipse with a few very shallow lobes, which is enough
    to look like a protein surface and still survives the favicon.
    """
    t = np.linspace(0, 2 * np.pi, 240, endpoint=False)
    r = 1.0 + wobble * np.sin(lobes * t + phase)
    x, y = rx * r * np.cos(t), ry * r * np.sin(t)
    c, s = math.cos(rot), math.sin(rot)
    return np.column_stack([cx + x * c - y * s, cy + x * s + y * c])


# The resurfaced patch. Rather than carving a dish, which reads as a mouth
# once a round effector sits above it, the rebuilt region is drawn along the
# receptor's own top edge: the outline it already has, restroked in white over
# the span the pipeline redesigns.
PATCH_X = 0.46          # half-width of the redesigned span, in axis units


def top_edge(outline: np.ndarray, cy: float, x_half: float) -> np.ndarray:
    """The upper side of a closed outline, within +/- x_half of centre.

    Keeping every point in the x window and sorting by x is not enough: the
    outline is a closed loop, so that returns the underside as well and the
    path zigzags between the two, which strokes as a filled block. Selecting
    on y above the outline's own centre takes the upper side alone.
    """
    sel = outline[(np.abs(outline[:, 0]) <= x_half) & (outline[:, 1] > cy)]
    return sel[np.argsort(sel[:, 0])]


def draw(ax, *, with_text: bool, ribs: bool = True) -> None:
    ax.add_patch(Polygon(rounded_hexagon(HEX_R, HEX_CORNER), closed=True,
                         facecolor=INK, edgecolor=RIM,
                         linewidth=_pt(RIM_W), joinstyle="round",
                         zorder=1))

    # The receptor fills the wide middle band, so the name can arc through the
    # upper taper where the hexagon narrows.
    ax.add_patch(Polygon(blob(0.0, -0.42, 0.64, 0.26, phase=1.1),
                         closed=True, facecolor=RECEPTOR, edgecolor="none",
                         zorder=3))

    # Ribs, so the receptor reads as a folded chain rather than a puddle. They
    # are the first thing to close up when the canvas shrinks, so the small
    # favicon frames leave them out rather than carry a smear.
    if ribs:
        for dx in (-0.34, -0.17, 0.17, 0.34):
            ax.plot([dx, dx], [-0.58, -0.38], color=INK,
                    solid_capstyle="round", linewidth=_pt(2.4), zorder=4)

    # The redesigned surface: the receptor's own top edge, restroked in white
    # across the span the pipeline rebuilds. The effector is drawn over its
    # middle, so the white shows on both flanks and reads as a patch the
    # partner is sitting on rather than as a line under it.
    edge = top_edge(blob(0.0, -0.42, 0.64, 0.26, phase=1.1), -0.42, PATCH_X)
    ax.plot(edge[:, 0], edge[:, 1], color=WHITE, solid_capstyle="round",
            linewidth=_pt(4.6), zorder=5)

    # The effector, seated on the patch. Its centre sits a little below the
    # edge so it reads as resting in the surface rather than hovering over it.
    ax.add_patch(Polygon(blob(0.0, -0.055, 0.235, 0.205, phase=0.4),
                         closed=True, facecolor=EFFECTOR, edgecolor="none",
                         zorder=6))

    if with_text:
        # Font size is in points too, so it goes through the same conversion.
        # ReceptorResurfacing is 19 glyphs against NegativeSteering's 16, so
        # it takes a wider spread and a smaller size to clear the rim.
        arched_text(ax, "ReceptorResurfacing", radius=0.670, y0=-0.015,
                    size=_pt(12.5), spread_deg=112.0)

    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_aspect("equal")
    ax.axis("off")


# Everything is drawn on one canvas this size and then resized down. Line
# widths and arrow heads are in points, which are absolute, so drawing a 16 px
# figure directly gives a logo whose strokes are twenty times too heavy for it.
# Draw once, big, and let the resampler do the shrinking.
BASE_PX = 512
DPI = 100


def _pt(px_at_240: float) -> float:
    """A width in px-at-240, in points on the drawing canvas.

    The canvas is BASE_PX/DPI inches across and spans 2 axis units, so one
    unit is (BASE_PX / 2 / DPI) * 72 points. Scaling stroke widths by
    BASE_PX/240 instead is not the point conversion and makes every stroke
    39% too heavy.
    """
    return px_at_240 * PX * (BASE_PX / 2.0 / DPI) * 72.0


def render(px: int, *, with_text: bool, ribs: bool = True):
    from PIL import Image
    fig = plt.figure(figsize=(BASE_PX / DPI, BASE_PX / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    fig.patch.set_alpha(0.0)
    draw(ax, with_text=with_text, ribs=ribs)
    tmp = HERE / ".logo_tmp.png"
    fig.savefig(tmp, dpi=DPI, transparent=True)
    plt.close(fig)
    im = Image.open(tmp).convert("RGBA")
    tmp.unlink()

    # Normalise to the reference's bounding box: crop to the drawn sticker,
    # scale it to (208, 240) at this size, and centre it in a square canvas.
    a = np.asarray(im)
    ys, xs = np.nonzero(a[..., 3] > 8)
    box = im.crop((int(xs.min()), int(ys.min()), int(xs.max()) + 1,
                   int(ys.max()) + 1))
    w = max(1, round(px * HEX_W_PX / OUT_PX))
    box = box.resize((w, px), Image.LANCZOS)
    out = Image.new("RGBA", (px, px), (0, 0, 0, 0))
    out.paste(box, ((px - w) // 2, 0))
    return out


def main() -> None:
    render(240, with_text=True).save(HERE / "logo.png")

    # The name is unreadable below 256 px and the ribs close up below 64, so
    # the small frames drop each in turn. A browser picks the frame nearest the
    # size it wants, so each one is drawn for that size rather than scaled down
    # from the 240.
    frames = []
    for px in (16, 32, 48, 64, 128, 256):
        frames.append(render(px, with_text=px >= 256, ribs=px >= 64))
    frames[-1].save(HERE / "favicon.ico", format="ICO",
                    sizes=[(f.width, f.height) for f in frames],
                    append_images=frames[:-1])
    print(f"wrote {HERE / 'logo.png'} and {HERE / 'favicon.ico'}")


if __name__ == "__main__":
    main()
