#!/usr/bin/env python3
"""Region length against design diversity, for the `full` contig run only.

Design region 1 has a 10-20 residue allowance, so RFdiffusion picked a different
length in different designs. That is the one comparison here free of any
between-run or between-region confound.

Diversity is measured within a length. Each design is compared only with the
other designs that came out the same length, and the mean of those distances is
its diversity. The obvious alternative, distance to every other design, cannot
be used here: RFdiffusion sampled lengths unevenly, centred near 15, so designs
at either end of the allowance score as distant simply for being a minority.
That produces a spurious U with high values at both 10 and 20. The within-length
measure removes it and leaves a monotonic rise.

Run with no arguments. Needs the .trb-corrected metrics, so run trb_metrics.py
first with PDB_ROOT and TRB_ROOT set.
"""
from __future__ import annotations

import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import diversity_model as dm  # noqa: E402
import novelty_common as nc  # noqa: E402
import trb_metrics as tm  # noqa: E402

warnings.filterwarnings("ignore")
nc.apply_style()

RUN = "full (2 regions)"
INTERIOR_C = "#0072B2"
TERMINAL_C = "#D55E00"

def within_length_diversity(region: int) -> pd.DataFrame:
    """Per design, mean chamfer distance to same-length designs of one region."""
    d = tm_json[region - 1]
    pts, lens = {}, {}
    for des in d:
        c = des["per_region_coords"][region - 1]
        if c:
            pts[des["design"]] = np.asarray(c, float)
            lens[des["design"]] = des["per_region_lengths"][region - 1]
    names = list(pts)
    n = len(names)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = dm._chamfer(pts[names[i]], pts[names[j]])
    L = np.array([lens[k] for k in names])
    rows = []
    for i in range(n):
        o = [j for j in range(n) if j != i and L[j] == L[i]]
        if o:
            rows.append({"design": names[i], "length": L[i],
                         "within": D[i, o].mean()})
    return pd.DataFrame(rows)


designs = json.loads((nc.DATA / f"{nc.RUNS[RUN]}{tm.OUT_SUFFIX}").read_text())["designs"]
tm_json = [designs, designs]

r1 = within_length_diversity(1)
r2 = within_length_diversity(2)
r1["log_within"] = np.log(r1.within)


def summarise(g):
    s = g.groupby("length").within.agg(n="size", mean="mean", sd="std").reset_index()
    s["ci95"] = 1.96 * s["sd"].fillna(0) / np.sqrt(s.n)
    return s


s1, s2 = summarise(r1), summarise(r2)

# A quadratic was fitted first, because the to-all measure looked U shaped. On
# the within-length measure the quadratic term is not supported, so the reported
# fit is a straight line on the raw scale, drawn straight rather than as an
# exponentiated log fit.
r1["Lc"] = r1.length - r1.length.mean()
r1["L2"] = r1.Lc ** 2
lin = smf.ols("within ~ Lc", r1).fit()
quad = smf.ols("within ~ Lc + L2", r1).fit()
ftest = quad.compare_f_test(lin)
per_res = lin.params["Lc"]
rho, rho_p = spearmanr(s1.length, s1["mean"])

fig, ax = plt.subplots(figsize=(7.0, 4.6))

xs = np.linspace(r1.length.min(), r1.length.max(), 120)
pred = lin.predict(pd.DataFrame({"Lc": xs - r1.length.mean()}))
ax.plot(xs, pred, color=INTERIOR_C, linewidth=1.8, zorder=3)

ax.errorbar(s1.length, s1["mean"], yerr=s1.ci95, fmt="o", color=INTERIOR_C,
            markersize=6, capsize=3, linewidth=1.2, zorder=4,
            markeredgecolor="white", markeredgewidth=0.6,
            label="region 1, interior (10-20 allowance)")
ax.errorbar(s2.length, s2["mean"], yerr=s2.ci95, fmt="D", color=TERMINAL_C,
            markersize=8, capsize=3, linewidth=1.2, zorder=4,
            markeredgecolor="white", markeredgewidth=0.6,
            label="region 2, C-terminal (6 residues, fixed)")

ax.set_xlabel("design region length RFdiffusion built (residues)")
ax.set_ylabel("diversity between same-length designs (Å)")
ax.set_title(f"Longer rebuilt regions give more varied backbones\n"
             f"{per_res:+.3f} A per added residue")
ax.legend(loc="upper left", fontsize=8.5)
fig.tight_layout()
fig.savefig(nc.PLOTS / "11_length_vs_diversity.png")
plt.close(fig)

stats = "\n".join([
    f"full region 1, n={len(r1)} designs, lengths {int(r1.length.min())}-{int(r1.length.max())}",
    "diversity = mean chamfer distance to same-length designs only",
    f"  linear     slope {per_res:+.4f} A/res  p={lin.pvalues['Lc']:.2e}  "
    f"R2={lin.rsquared:.3f}",
    f"  quadratic  quad  {quad.params['L2']:+.5f}  p={quad.pvalues['L2']:.3f}  "
    f"(F-test vs linear p={ftest[1]:.3f}, not supported)",
    f"  group-level Spearman over the {len(s1)} lengths: rho={rho:+.2f}, p={rho_p:.4f}",
    f"  diversity at 10 res {s1[s1.length == 10]['mean'].iloc[0]:.2f} A, "
    f"at 20 res {s1[s1.length == 20]['mean'].iloc[0]:.2f} A",
    f"full region 2 (C-terminal, 6 res, n={len(r2)}): "
    f"{s2['mean'].iloc[0]:.2f} +- {s2.ci95.iloc[0]:.2f} A",
    "",
    "NOTE: distance-to-all-designs instead of within-length produces a spurious",
    "U shape, because RFdiffusion sampled lengths unevenly around 15 and the",
    "extremes score as distant for being a minority rather than for varying.",
])
(nc.PLOTS.parent / "length_diversity_stats.txt").write_text(stats + "\n")
print("  plots/11_length_vs_diversity.png")
print(stats)
