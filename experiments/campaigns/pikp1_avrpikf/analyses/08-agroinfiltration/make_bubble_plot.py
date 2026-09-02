#!/usr/bin/env python3
"""Cell death bubble plot for the AVR-PikF agroinfiltration, 7 dpi.

Layout follows the dot plots in Maidment (2019, figures 5.3B to 5.5B), turned
on its side: one row per infiltration condition, cell death score along the x
axis, and a dot at each score sized by how many infiltration sites received it.
Twelve conditions on the x axis made the figure far too wide for the page.

The three files join like this. `treatment` in the CDALabeller output is the
condition number in column A of the infiltration spreadsheet; the condition
names its receptor as Pik_R_NN; and NN is the agro_id in design_id_map.csv,
which gives the RFdiffusion design and ProteinMPNN sequence.

Five of the ten selected designs are here. R01, R05, R06 and R10 gave almost no
cell death at 5 dpi, so they were repeated at 7 dpi together with R04, which
was autoactive. The other five are not in this experiment.

Row order is the two Pik-HIPP43 conditions first, then each design beside its
own PWL2 control, so autoactivity reads as a pair rather than having to be
looked up elsewhere. Everything carrying PWL2 is a control and is coloured as
one, since PWL2 was not expected to be recognised by anything here; only a
design with AVR-PikF is a test condition.

Usage:
    make_bubble_plot.py [--outdir plots]
"""

import argparse
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
# The infiltration sheet is lab record, not analysis output, so it is read
# from the run tree rather than duplicated here.
AGRO = os.path.join(HERE, "..", "..", "runs", "crystal_full_test_contig",
                    "agroinfiltrations", "7dpi_5conditions")
CONDITIONS = os.path.join(
    AGRO, "20260629_coinfiltration_calculation.xlsx")

INK, MUTED = "#222222", "#666666"
CONTROL, DESIGN = "#D62728", "#1F6FB4"
MAX_SCORE = 6


def load():
    cond = (pd.read_excel(CONDITIONS,
                          header=None).iloc[15:27, [0, 1]])
    cond.columns = ["treatment", "condition"]
    cond["treatment"] = cond.treatment.astype(int)

    idmap = pd.read_csv(os.path.join(DATA, "design_id_map.csv"))
    idmap["agro_id"] = idmap.agro_id.astype(str).str.zfill(3)
    to_design = dict(zip(idmap.agro_id, idmap.design))

    def receptor(c):
        m = re.search(r"Pik_R_0?(\d+)", c)
        if not m:
            return "Pik-HIPP43", None
        agro = (m.group(1).lstrip("0") or "0").zfill(3)
        return f"Pik_R_{agro}", to_design.get(agro)

    cond[["receptor", "design"]] = cond.condition.apply(
        lambda c: pd.Series(receptor(c)))
    cond["effector"] = ["PWL2" if "PWL2" in c else "AVR-PikF"
                        for c in cond.condition]

    scores = pd.read_csv(os.path.join(DATA, "cell_death_scores.csv"))
    scores["treatment"] = scores.treatment.astype(int)
    return cond, scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "plots"))
    ap.add_argument("--counts", action="store_true",
                    help="write the site count inside each bubble")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.labelsize": 9, "axes.edgecolor": MUTED,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.grid": False, "legend.frameon": False})

    cond, scores = load()
    controls = cond[cond.design.isna()].sort_values("treatment")
    designs = cond[cond.design.notna()]

    # Controls first, then each design next to its own PWL2 pair.
    order = list(controls.treatment)
    for d in designs[designs.effector == "AVR-PikF"].sort_values("treatment").itertuples():
        order.append(d.treatment)
        pair = designs[(designs.design == d.design) & (designs.effector == "PWL2")]
        order += list(pair.treatment)

    info = cond.set_index("treatment")
    AREA = 24                       # dot area per infiltration site
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    labels = []
    for y, t in enumerate(order):
        r = info.loc[t]
        # A control is anything not testing a design against AVR-PikF, so the
        # per-design PWL2 rows are controls too and are coloured as such.
        is_test = isinstance(r.design, str) and r.effector == "AVR-PikF"
        colour = DESIGN if is_test else CONTROL
        counts = scores[scores.treatment == t].score.value_counts()
        ax.scatter(counts.index, [y] * len(counts), s=AREA * counts.values,
                   color=colour, alpha=0.75, edgecolor=INK, linewidth=0.5)
        # The house style carries a size legend only, but their n runs to 30-60
        # leaves per condition. At 6 per design a reader will try to count the
        # dots, so the count is written above each. Inside the bubble does not
        # work: at one site the digit is as wide as the circle.
        if args.counts:
            for v, c in counts.items():
                # Offset clears the bubble, whose radius in points is
                # sqrt(area / pi) and so grows with the count.
                pad = math.sqrt(AREA * c / math.pi) + 3.5
                ax.annotate(str(c), (v, y), xytext=(0, pad),
                            textcoords="offset points", ha="center",
                            va="bottom", fontsize=7, color=INK, zorder=5)
        receptor = r.design if isinstance(r.design, str) else "Pik-HIPP43"
        labels.append(f"{receptor} + {r.effector}")

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()                # first condition at the top
    ax.set_ylim(len(order) - 0.4, -0.85)   # head room for the top row's labels
    ax.set_xlim(-0.5, MAX_SCORE + 0.5)
    ax.set_xticks(range(MAX_SCORE + 1))
    ax.set_xlabel("Cell death severity")
    # A rule after the controls and after each design/control pair, so the
    # rows group visually the way they group logically.
    for y in range(1, len(order) - 1, 2):
        ax.axhline(y + 0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)))

    sizes = sorted(scores.groupby(["treatment", "score"]).size().unique())
    sizes = [sizes[0], sizes[len(sizes) // 2], sizes[-1]]
    handles = [plt.scatter([], [], s=AREA * n, color=MUTED, alpha=0.75,
                           edgecolor=INK, linewidth=0.5, label=str(n))
               for n in sizes]
    size_leg = ax.legend(handles=handles, loc="upper left",
                         bbox_to_anchor=(1.02, 1.0), labelspacing=1.1,
                         fontsize=8.5, title="Infiltration sites",
                         title_fontsize=8.5)
    size_leg._legend_box.align = "left"
    ax.add_artist(size_leg)

    # Colour legend as its own artist, dropped well below the size one rather
    # than sharing its spacing, which crowded the two together.
    colour_handles = [
        plt.scatter([], [], s=AREA * 4, color=DESIGN, alpha=0.75,
                    edgecolor=INK, linewidth=0.5, label="Resurfaced HMAs"),
        plt.scatter([], [], s=AREA * 4, color=CONTROL, alpha=0.75,
                    edgecolor=INK, linewidth=0.5, label="Control")]
    colour_leg = ax.legend(handles=colour_handles, loc="upper left",
                           bbox_to_anchor=(1.02, 0.52), labelspacing=0.9,
                           fontsize=8.5)
    colour_leg._legend_box.align = "left"
    fig.tight_layout()
    name = "cell_death_bubble_counts" if args.counts else "cell_death_bubble"
    out = os.path.join(args.outdir, f"{name}.png")
    fig.savefig(out)
    print("wrote", out)

    print("\nmedian score per condition:")
    for t in order:
        r = info.loc[t]
        name = r.design if isinstance(r.design, str) else "Pik-HIPP43"
        print(f"  {t:2d}  {name:8s} + {r.effector:9s} "
              f"median {scores[scores.treatment == t].score.median():.1f}")


if __name__ == "__main__":
    main()
