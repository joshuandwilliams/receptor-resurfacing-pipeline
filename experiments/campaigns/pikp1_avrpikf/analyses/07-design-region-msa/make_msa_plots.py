#!/usr/bin/env python3
"""Two-panel MSA figure for the two design regions of the ten selected designs.

Rebuilt from the ClustalW alignments the pipeline's run_alignment.py produced
(cached in data/), so ClustalW is not re-run and the alignments cannot drift
from the ones already reported.

Three deliberate departures from run_alignment.py:

  * Both panels use its PALE colour scheme.  The original applied the pale
    colours to region 1 and the saturated ones to region 2, which made the two
    look like different kinds of data.
  * Rows are ordered by the selection table rather than by ClustalW's output
    order, which differs between the two regions.  A shared row order is what
    lets a reader read across the panels.
  * Region 2 carries an extra reference row, SNK-EKE, a known engineered Pikp-1
    variant at these positions.  It has no structure in this campaign, so its
    sequence is supplied directly rather than extracted from a PDB.
  * Region 1 is read from region1_anchored.aln rather than region1.aln.
    ClustalW put every gap at the right-hand end, which left the second beta
    strand scattered across three columns; see data/anchor_alignment.py for
    why no gap penalty fixes that and what is done instead. Region 2 is a
    fixed six residues and still comes straight from ClustalW.

A colour key is drawn beneath the panels, since the palette is Clustal-X by
physicochemical group and nothing on the figure said so.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

HERE = Path(__file__).parent
DATA = HERE / "data"
PLOTS = HERE / "plots"

WILD_TYPE_ID = "crystal_input"
WILD_TYPE_LABEL = "Wild Type"

# Extra reference row for region 2 only, drawn directly under the wild type.
SNK_EKE = ("SNK-EKE", "EQAKED")

# Selection-table order, so the two panels share a row order and both match
# the order the designs are reported in.
DESIGN_ORDER = ["d39_s0", "d17_s0", "d47_s0", "d39_s1", "d26_s1",
                "d18_s0", "d49_s0", "d40_s0", "d6_s1", "d14_s0"]

# run_alignment.py's pale Clustal-X palette, applied to both panels here.
# Same six side-chain classes as the composition figure in
# 03-mpnn-sequence-design, so one key serves both. Colours are pale tints of
# the saturated Okabe-Ito hues used there, since letters sit on top here.
CLASS_OF = {
    **{a: "Aliphatic" for a in "AVLIM"},
    **{a: "Aromatic" for a in "FWY"},
    **{a: "Polar" for a in "STNQC"},
    **{a: "Positive" for a in "KRH"},
    **{a: "Negative" for a in "DE"},
    **{a: "Special" for a in "GP"},
}
# Exactly the fills used by the composition figure, not tints of them, so the
# two figures are the same colour. Letters sit on top here, so the text colour
# is chosen per cell from the fill's luminance.
CLASS_COLOUR = {
    "Aliphatic": "#0072B2",
    "Aromatic": "#56B4E9",
    "Polar": "#009E73",
    "Positive": "#CC79A7",
    "Negative": "#E69F00",
    "Special": "#999999",
}


def text_on(hex_colour):
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return "#222222" if 0.299 * r + 0.587 * g + 0.114 * b > 0.6 else "#ffffff"
RESIDUE_COLOUR_PALE = {aa: CLASS_COLOUR[c] for aa, c in CLASS_OF.items()}
RESIDUE_COLOUR_PALE["-"] = "#ffffff"


def read_aln(path: Path) -> dict[str, str]:
    """Minimal CLUSTAL reader. These alignments are single-block."""
    seqs: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("CLUSTAL") or line.startswith(" "):
            continue
        parts = line.split()
        if len(parts) == 2:
            seqs[parts[0]] = seqs.get(parts[0], "") + parts[1]
    return seqs


def label_for(mid: str) -> str:
    if mid == WILD_TYPE_ID:
        return WILD_TYPE_LABEL
    m = re.match(r"d(\d+)_s(\d+)$", mid)
    return f"Design {m.group(1)} Sequence {m.group(2)}" if m else mid


def ordered_rows(seqs: dict[str, str], extra: tuple[str, str] | None,
                 by_length: bool = False, order: list[str] | None = None) -> list:
    """Wild type first, then any extra reference, then the designs.

    by_length sorts the designs longest first rather than by the selection
    table. For region 1 that groups rows whose loops are the same size, so a
    column compares turn against turn instead of turn against strand, which is
    what made the alignment look ragged. Region 2 is a fixed six residues, so it
    cannot derive that order itself and takes region 1's via order, which keeps
    a design on the same row in both panels.
    """
    rows = [(WILD_TYPE_LABEL, seqs[WILD_TYPE_ID])]
    if extra is not None:
        width = len(seqs[WILD_TYPE_ID])
        name, seq = extra
        rows.append((name, seq.ljust(width, "-")[:width]))
    if order is None:
        order = DESIGN_ORDER
        if by_length:
            order = sorted((m for m in DESIGN_ORDER if m in seqs),
                           key=lambda m: (-len(seqs[m].replace("-", "")), m))
    for mid in order:
        if mid in seqs:
            rows.append((label_for(mid), seqs[mid]))
    return rows



def draw_panel(ax, rows, title):
    ncol = max(len(s) for _, s in rows)
    nrow = len(rows)
    for i, (_, seq) in enumerate(rows):
        y = nrow - 1 - i
        for j, aa in enumerate(seq):
            ax.add_patch(Rectangle((j, y), 1, 1,
                                   facecolor=RESIDUE_COLOUR_PALE.get(aa, "#dddddd"),
                                   edgecolor="white", linewidth=0.5))
            if aa != "-":
                ax.text(j + 0.5, y + 0.5, aa, ha="center", va="center",
                        fontsize=8, family="monospace",
                        color=text_on(RESIDUE_COLOUR_PALE.get(aa, "#dddddd")))
    ax.set_xlim(0, ncol)
    ax.set_ylim(0, nrow)
    ax.set_yticks([nrow - 0.5 - i for i in range(nrow)])
    ax.set_yticklabels([n for n, _ in rows], fontsize=8, family="monospace")
    ax.set_xticks([j + 0.5 for j in range(ncol)])
    ax.set_xticklabels(range(1, ncol + 1), fontsize=7)
    ax.set_xlabel("Alignment column")
    # Left-aligned, because the two panels share an x-limit but not a column
    # count, so a centred title would float away from the narrower panel's data.
    ax.set_title(title, fontsize=10, loc="left")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


# Groups as Clustal-X colours them, for the key.
COLOUR_KEY = [(CLASS_COLOUR[c], c) for c in
              ["Aliphatic", "Aromatic", "Polar", "Positive", "Negative",
               "Special"]]


def draw_key(path, ncol=6, pad_to=None):
    """The colour key as its own image, so the assembler can stack it.

    pad_to widens the saved PNG to the panels' pixel width. The assembler
    scales every image to one shared display width, so a narrow key would be
    upscaled and its text would come out larger than the panel labels.
    """
    fig = plt.figure(figsize=(7.4, 0.5))
    handles = [Patch(facecolor=c, edgecolor="#bbbbbb", linewidth=0.5, label=l)
               for c, l in COLOUR_KEY]
    fig.legend(handles=handles, loc="center", ncol=ncol, frameon=False,
               fontsize=8, handlelength=1.4, columnspacing=1.6)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    if pad_to:
        im = Image.open(path)
        if im.width < pad_to:
            canvas = Image.new("RGB", (pad_to, im.height), "white")
            canvas.paste(im, ((pad_to - im.width) // 2, 0))
            canvas.save(path)
            im = canvas
        print(f"{path}  colour key, padded to {im.width}px")
        return
    print(f"{path}  colour key")


def main():
    PLOTS.mkdir(exist_ok=True)
    a1 = read_aln(DATA / "region1_anchored.aln")
    r1 = ordered_rows(a1, None, by_length=True)
    r1_order = sorted((m for m in DESIGN_ORDER if m in a1),
                      key=lambda m: (-len(a1[m].replace("-", "")), m))
    r2 = ordered_rows(read_aln(DATA / "region2.aln"), SNK_EKE, order=r1_order)

    n1 = max(len(s) for _, s in r1)
    n2 = max(len(s) for _, s in r2)

    # One PNG per panel. The A/B letters belong in the thesis SVG template, added
    # by assemble_figure.py, not drawn into the plot. Both panels are rendered on
    # region 1's column range so that when the assembler scales them to a shared
    # display width the cells end up the same size in each.
    for rows, ncols_used, title, name in [
            (r1, n1, "Design region 1", "region1_msa.png"),
            (r2, n2, "Design region 2", "region2_msa.png")]:
        fig, ax = plt.subplots(figsize=(n1 * 0.42 + 3.2, len(rows) * 0.42 + 0.9))
        draw_panel(ax, rows, title)
        ax.set_xlim(0, n1)
        if ncols_used != n1:
            # Only the x position, so matplotlib still places the label just
            # under the tick labels. set_label_coords sets y in axes fractions
            # too, which drops the label much further on a taller panel.
            ax.xaxis.label.set_x(ncols_used / (2.0 * n1))
        out = PLOTS / name
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"{out}  {len(rows)} rows x {ncols_used} cols")
    panel_w = Image.open(PLOTS / "region1_msa.png").width
    draw_key(PLOTS / "colour_key.png", pad_to=panel_w)


if __name__ == "__main__":
    main()
