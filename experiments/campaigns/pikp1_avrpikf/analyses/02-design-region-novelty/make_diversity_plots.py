#!/usr/bin/env python3
"""One figure for how region size and placement relate to design diversity.

The model in diversity_model.py finds a positive size effect that varies in
strength between regions, and no placement effect that survives a region-level
test. Two plain panels say that without asking the reader to interpret
coefficients. Statistics go in the caption, not the axes.

Reads the .trb-corrected metrics, so run trb_metrics.py first. On the
uncorrected pipeline output the size effect disappears and scaffold drift looks
like the dominant term, which is the residue-mapping bug rather than a result.

Run with no arguments. The figure lands in plots/ and the fitted numbers in
diversity_model_stats.txt.
"""
from __future__ import annotations

import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import diversity_model as dm  # noqa: E402
import novelty_common as nc  # noqa: E402

warnings.filterwarnings("ignore")
nc.apply_style()
nc.PLOTS.mkdir(exist_ok=True)

df = dm.build_table()

fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)

# A. Size. Bigger regions do give more varied backbones, though the slope
# differs between regions, so the line is the average and not a law.
ax = axes[0]
for run, c in nc.RUN_COLOR.items():
    g = df[df.run == run]
    ax.scatter(g.length, g.chamfer_div, s=13, color=c, alpha=0.45,
               edgecolor="white", linewidth=0.25, zorder=3, label=run)
b = np.polyfit(df.length, df.chamfer_div, 1)
xs = np.linspace(df.length.min(), df.length.max(), 20)
ax.plot(xs, np.polyval(b, xs), color=nc.INK, linewidth=1.6, zorder=4)
ax.set_xlabel("design region size (residues)")
ax.set_ylabel("diversity between designs (Å)")
ax.set_title("A. Bigger regions, more varied backbones")
ax.legend(loc="upper left", markerscale=1.8, handletextpad=0.2)

# B. Placement. Terminal regions are anchored on one side only, and look more
# diverse, but there are three of them so the difference does not hold up.
ax = axes[1]
groups = [("interior\n(9 regions)", df[df.terminal == 0]),
          ("terminal\n(3 regions)", df[df.terminal == 1])]
bp = ax.boxplot([g.chamfer_div.values for _, g in groups], patch_artist=True,
                widths=0.5, showfliers=False,
                medianprops=dict(color=nc.INK, linewidth=1.6))
for patch in bp["boxes"]:
    patch.set_facecolor(nc.MUTED)
    patch.set_alpha(0.25)
    patch.set_edgecolor(nc.MUTED)
for i, (_, g) in enumerate(groups, start=1):
    jit = np.random.default_rng(i).uniform(-0.15, 0.15, len(g))
    ax.scatter(i + jit, g.chamfer_div, s=10,
               color=[nc.RUN_COLOR[r] for r in g.run], alpha=0.45, zorder=3,
               edgecolor="white", linewidth=0.25)
ax.set_xticks([1, 2])
ax.set_xticklabels([lbl for lbl, _ in groups])
ax.set_xlabel("how the region is anchored")
ax.set_title("B. Placement, not resolved by 3 regions")
ax.grid(axis="x", visible=False)

fig.tight_layout()
fig.savefig(nc.PLOTS / "09_what_drives_diversity.png")
plt.close(fig)
print("  plots/09_what_drives_diversity.png")


# ── Numbers for the caption ──────────────────────────────────────────────────
lines = []
common = dm.fit(dm.varying_length(df), "log_chamfer")
varying = dm.fit_random_slope(df, "log_chamfer")
lrt = dm.slope_lrt(df, "log_chamfer")
slopes = dm.per_region_slopes(df, "log_chamfer")
perm_p, n_perm = dm.permutation_p(df, "log_chamfer")

ci = varying.conf_int().loc["length_within"]
lines.append("SIZE")
lines.append(f"  per-region slopes positive: {slopes.attrs['n_positive']}/{len(slopes)}"
             f"  (sign test p = {slopes.attrs['sign_test_p']:.2f})")
lines.append(f"  common slope   {lrt['beta_common']:+.4f}")
lines.append(f"  varying slope  {lrt['beta_varying']:+.4f} "
             f"[{ci[0]:+.4f}, {ci[1]:+.4f}]  p = {varying.pvalues['length_within']:.2f}")
lines.append(f"  LRT for slope variation chi2 = {lrt['chi2']:.2f}, p = {lrt['p']:.4f}")

lines.append("PLACEMENT")
ri = dm.fit(df, "log_chamfer")
lines.append(f"  terminal {ri.params['terminal']:+.3f}, Wald p = "
             f"{ri.pvalues['terminal']:.3f}, exact permutation p = {perm_p:.3f} "
             f"({n_perm} assignments, 3 terminal regions of 12)")

lines.append("SCAFFOLD DRIFT (now a trivial range after the .trb correction)")
d = ri.conf_int().loc["motif_rmsd_within"]
lines.append(f"  motif_rmsd_within {ri.params['motif_rmsd_within']:+.3f} "
             f"[{d[0]:+.3f}, {d[1]:+.3f}]  p = {ri.pvalues['motif_rmsd_within']:.1e}")

text = "\n".join(lines)
(nc.PLOTS.parent / "diversity_model_stats.txt").write_text(text + "\n")
slopes.to_csv(nc.PLOTS.parent / "diversity_size_slopes_by_region.csv", index=False)
print()
print(text)


# ── 10. All twelve regions at once ───────────────────────────────────────────
# Size, anchoring and what the region was cut out of are three separate claims,
# and with twelve regions they cannot be separated by a model. Drawing every
# region with all three properties shown lets the confounding be read directly.
import ss_context as sc  # noqa: E402

ss = sc.annotate()
reg = (df.groupby(["run", "region", "region_id"], as_index=False)
         .agg(size=("length", "median"), diversity=("chamfer_div", "median"))
         .merge(ss[["region_id", "terminal", "ss_class"]], on="region_id"))
reg.to_csv(nc.PLOTS.parent / "region_context_summary.csv", index=False)

SS_COLOR = {"strand": "#0072B2", "loop": "#E69F00", "mixed": "#999999",
            "helix": "#009E73"}
SHORT = {"full (2 regions)": "full", "auto 5 A (4 regions)": "5Å",
         "auto 3 A (6 regions)": "3Å"}

fig, ax = plt.subplots(figsize=(7.4, 5.0))
# The three single-residue regions sit almost on top of each other, so labels
# are staggered when a point is close to one already drawn.
placed = []
for row in reg.sort_values(["size", "diversity"]).itertuples():
    ax.scatter(row.size, row.diversity, s=150 if row.terminal else 90,
               marker="D" if row.terminal else "o",
               color=SS_COLOR[row.ss_class], alpha=0.85, zorder=3,
               edgecolor="white", linewidth=1.0)
    # Shift sideways, never vertically, so a label always sits on the same
    # row as the point it names.
    crowded = sum(1 for x, y in placed
                  if abs(x - row.size) < 1.0 and abs(y - row.diversity) < 0.22)
    ax.annotate(f"{SHORT[row.run]}·{row.region}",
                (row.size, row.diversity), textcoords="offset points",
                xytext=(11 + 36 * crowded, -3), fontsize=8, color=nc.MUTED)
    placed.append((row.size, row.diversity))

ax.set_xlabel("design region size (median residues)")
ax.set_ylabel("diversity between designs (median Å)")
ax.set_title("10. Every design region, with all three properties\n"
             "diamonds are terminal, colour is what the region replaced")

handles = [plt.Line2D([], [], marker="o", linestyle="", color=SS_COLOR["strand"],
                      label="replaces strand", markersize=9),
           plt.Line2D([], [], marker="o", linestyle="", color=SS_COLOR["loop"],
                      label="replaces loop", markersize=9),
           plt.Line2D([], [], marker="o", linestyle="", color=nc.MUTED,
                      label="interior (2 anchors)", markersize=9),
           plt.Line2D([], [], marker="D", linestyle="", color=nc.MUTED,
                      label="terminal (1 anchor)", markersize=9)]
ax.legend(handles=handles, loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(nc.PLOTS / "10_region_context.png")
plt.close(fig)
print("\n  plots/10_region_context.png")
print(reg.sort_values("diversity").to_string(index=False))
