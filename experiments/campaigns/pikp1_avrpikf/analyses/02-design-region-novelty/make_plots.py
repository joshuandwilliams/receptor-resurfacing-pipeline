#!/usr/bin/env python3
"""Candidate plots for the design-region novelty comparison.

Six candidates, deliberately overlapping, so the useful one can be picked rather
than assumed. Run with no arguments; everything lands in plots/.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import novelty_common as nc  # noqa: E402

nc.apply_style()
nc.PLOTS.mkdir(exist_ok=True)

regions = nc.load_regions()
designs = nc.load_designs()
spread = nc.region_com_spread()

RUN_HANDLES = [Line2D([], [], color=c, linewidth=8, alpha=0.55, label=r)
               for r, c in nc.RUN_COLOR.items()]


def save(fig, name):
    fig.savefig(nc.PLOTS / name)
    plt.close(fig)
    print(f"  plots/{name}")


# ── 1. The plot as asked for: COM displacement per region, per run ───────────
fig, ax = plt.subplots(figsize=(8.0, 4.0))
nc._box_by_run(ax, regions, "com_displacement")
ax.set_xlabel("design region (numbered within each run)")
ax.set_ylabel("centre-of-mass displacement (Å)")
ax.set_title("1. COM displacement per design region")
ax.legend(handles=RUN_HANDLES, loc="upper left")
save(fig, "01_com_displacement_by_region.png")

# ── 2. Pooled across regions, one box per run ────────────────────────────────
fig, ax = plt.subplots(figsize=(5.4, 4.0))
nc._box_by_run(ax, regions, "com_displacement", by_region=False)
ax.set_ylabel("centre-of-mass displacement (Å)")
ax.set_title("2. COM displacement pooled across regions")
save(fig, "02_com_displacement_pooled.png")

# ── 3. The confound: displacement tracks region length ───────────────────────
fig, ax = plt.subplots(figsize=(5.8, 4.4))
for run, c in nc.RUN_COLOR.items():
    sub = regions[regions.run == run].dropna(subset=["length", "com_displacement"])
    ax.scatter(sub.length, sub.com_displacement, s=22, color=c, alpha=0.6,
               edgecolor="white", linewidth=0.4, label=run, zorder=3)
ok = regions.dropna(subset=["length", "com_displacement"])
rho = np.corrcoef(ok.length.rank(), ok.com_displacement.rank())[0, 1]
ax.set_xlabel("design region length (residues)")
ax.set_ylabel("centre-of-mass displacement (Å)")
ax.set_title(f"3. Displacement against region length  (ρ = {rho:+.2f})")
ax.legend(loc="upper left")
save(fig, "03_com_vs_length.png")

# ── 4. Length-normalised, so short regions are not penalised ─────────────────
fig, ax = plt.subplots(figsize=(5.4, 4.0))
nc._box_by_run(ax, regions, "com_per_residue", by_region=False)
ax.set_ylabel("COM displacement per designed residue (Å)")
ax.set_title("4. Length-normalised displacement")
save(fig, "04_com_per_residue.png")

# ── 5. Span ratio: stretched or pulled tight relative to the input gap ───────
fig, ax = plt.subplots(figsize=(5.4, 4.0))
nc._box_by_run(ax, regions.dropna(subset=["span_ratio"]), "span_ratio", by_region=False)
ax.axhline(1.0, color=nc.MUTED, linestyle="--", linewidth=1.2, zorder=1)
ax.set_ylabel("design span ÷ input span")
ax.set_title("5. Does the built region span the same reach?")
save(fig, "05_span_ratio.png")

# ── 6. Diversity: how different are designs from EACH OTHER? ─────────────────
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
for run, c in nc.RUN_COLOR.items():
    sub = spread[spread.run == run]
    axes[0].scatter(sub.region, sub.mean_pairwise_com_dist, s=44, color=c,
                    alpha=0.8, edgecolor="white", linewidth=0.5, label=run, zorder=3)
axes[0].set_xlabel("design region")
axes[0].set_ylabel("mean pairwise COM distance (Å)")
axes[0].set_title("6a. Diversity between designs, per region")
axes[0].legend(loc="upper right")

nc._box_by_run(axes[1], designs.assign(region=1), "total_designed", by_region=False)
axes[1].set_ylabel("total designed residues per design")
axes[1].set_title("6b. How much was rebuilt at all")
fig.tight_layout()
save(fig, "06_design_diversity.png")

# ── Run-level summary table ──────────────────────────────────────────────────
summary = (designs.groupby("run", sort=False)
           .agg(designs=("design", "size"),
                passing=("passes_filter", "sum"),
                median_total_designed=("total_designed", "median"),
                median_motif_rmsd=("motif_rmsd", "median"),
                median_coverage=("design_region_coverage", "median"))
           .assign(pass_rate=lambda d: (d.passing / d.designs).round(3)))
summary.to_csv(nc.PLOTS.parent / "run_summary.csv")
print(summary.to_string())

by_region = (regions.groupby(["run", "region"], sort=False)
             .agg(n=("com_displacement", "size"),
                  median_length=("length", "median"),
                  median_com=("com_displacement", "median"),
                  median_com_per_res=("com_per_residue", "median")).round(3))
by_region.to_csv(nc.PLOTS.parent / "region_summary.csv")
print()
print(by_region.to_string())


# ── 7. The design-level view: does the contig give RFDiffusion room? ─────────
# Region-level plots split the very thing the contig controls. Aggregating to
# the design keeps "how much was rebuilt" and "how far it moved" together.
agg = (regions.assign(w=lambda d: d.com_displacement * d.length)
       .groupby(["run", "design"], sort=False)
       .agg(total_designed=("length", "sum"),
            weighted_com=("w", "sum"),
            max_com=("com_displacement", "max"))
       .reset_index())
agg["mean_com_per_design"] = agg.weighted_com / agg.total_designed

fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.0))
for ax, col, title, ylab in [
        (axes[0], "total_designed", "7a. Room to work",
         "designed residues per design"),
        (axes[1], "weighted_com", "7b. Total geometric change",
         "Σ (COM displacement × region length) (Å·res)"),
        (axes[2], "max_com", "7c. Largest single-region move",
         "max COM displacement (Å)")]:
    nc._box_by_run(ax, agg.assign(region=1), col, by_region=False)
    ax.set_title(title)
    ax.set_ylabel(ylab)
    ax.set_xticklabels([t.get_text().replace(" (", "\n(") for t in ax.get_xticklabels()])
fig.tight_layout()
save(fig, "07_design_level_room.png")

# ── 8. Diversity per run: are the designs different from each other? ─────────
fig, ax = plt.subplots(figsize=(6.0, 4.0))
nc._box_by_run(ax, spread.assign(region=1), "mean_pairwise_com_dist", by_region=False)
ax.set_ylabel("mean pairwise COM distance between designs (Å)")
ax.set_title("8. Design-to-design diversity, one point per region")
ax.set_xticklabels([t.get_text().replace(" (", "\n(") for t in ax.get_xticklabels()])
save(fig, "08_diversity_by_run.png")

print()
print(agg.groupby("run", sort=False)[
    ["total_designed", "weighted_com", "max_com"]].median().round(2).to_string())
