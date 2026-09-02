"""Mixed model for how design region size and placement drive diversity.

Plot 08 collapses diversity to one number per region, which leaves 12 points to
carry both predictors. Here diversity is computed per (design, region) instead,
so region identity becomes a random effect and the two predictors separate.

Size varies within a region, because the contigs allow length ranges, so its
effect is estimable within region and is free of the run confound. Placement is
fixed per region, so it rests on 12 regions and is reported as descriptive.

Two response variables. `centroid_dev` matches the existing centroid-based
measure. `chamfer_div` is a symmetric nearest-neighbour distance between two
designs' Ca coordinate sets, which is defined for unequal lengths and responds
to shape rather than only to where the centre landed.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import novelty_common as nc
import trb_metrics as tm

# The last contig segment is a design block in all three runs, so exactly one
# region per run is anchored on a single side. Everything else has fixed
# scaffold on both flanks.
N_DESIGN_REGIONS = {label: sum(1 for k, *_ in nc.parse_contig(c) if k == "design")
                    for label, c in nc.CONTIGS.items()}


def _chamfer(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric mean nearest-neighbour distance between two Ca point sets."""
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return 0.5 * (d.min(axis=1).mean() + d.min(axis=0).mean())


def _load(run: str) -> dict:
    """Prefer the .trb-corrected metrics, and say so if they are missing.

    The pipeline's own per_region_lengths are a proportional guess and are wrong
    for most designs of both auto runs, which shifts every downstream residue.
    See trb_metrics.py. Fitting this model on the uncorrected file measures the
    bug, so it is refused rather than silently allowed.
    """
    corrected = nc.DATA / f"{run}{tm.OUT_SUFFIX}"
    if not corrected.exists():
        raise FileNotFoundError(
            f"{corrected.name} not found. Run trb_metrics.py first, with "
            f"PDB_ROOT and TRB_ROOT set.")
    return json.loads(corrected.read_text())


def build_table() -> pd.DataFrame:
    """One row per (run, design, region) with both diversity responses."""
    rows = []
    for label, run in nc.RUNS.items():
        d = _load(run)
        n_regions = N_DESIGN_REGIONS[label]

        coords: dict[int, list] = {}
        meta: dict[int, list] = {}
        for des in d["designs"]:
            lens = des.get("per_region_lengths") or []
            for i, c in enumerate(des.get("per_region_coords") or []):
                if not c:
                    continue
                coords.setdefault(i, []).append(np.asarray(c, float))
                meta.setdefault(i, []).append(
                    (des["design"], lens[i] if i < len(lens) else np.nan,
                     nc._num(des.get("motif_rmsd"))))

        for i, pts in coords.items():
            cents = np.asarray([p.mean(axis=0) for p in pts])
            mean_cent = cents.mean(axis=0)
            dev = np.linalg.norm(cents - mean_cent, axis=1)

            n = len(pts)
            cham = np.zeros((n, n))
            for a in range(n):
                for b in range(a + 1, n):
                    cham[a, b] = cham[b, a] = _chamfer(pts[a], pts[b])
            cham_mean = cham.sum(axis=1) / (n - 1) if n > 1 else np.full(n, np.nan)

            for j, (design, length, rmsd) in enumerate(meta[i]):
                rows.append({
                    "run": label,
                    "region": i + 1,
                    "region_id": f"{label}|{i + 1}",
                    "design": design,
                    "length": float(length),
                    "motif_rmsd": rmsd,
                    "terminal": int(i + 1 == n_regions),
                    "centroid_dev": float(dev[j]),
                    "chamfer_div": float(cham_mean[j]),
                })

    df = pd.DataFrame(rows)
    # Within/between decomposition. The `_within` term varies across designs
    # inside one region and is what the 64 designs per region actually inform.
    # The `_between` term is the region mean and rests on 12 regions. Without
    # this split the motif_rmsd coefficient would mix the per-design effect with
    # the between-run difference in scaffold drift, which are different claims.
    for col in ("length", "motif_rmsd"):
        mean = df.groupby("region_id")[col].transform("mean")
        df[f"{col}_between"] = mean
        df[f"{col}_within"] = df[col] - mean

    # Both responses are non-negative distances and right skewed. On the raw
    # scale the chamfer residuals reach skew 2.0, so the log scale is the primary
    # one and the raw scale is kept as a sensitivity check.
    df["log_centroid"] = np.log(df.centroid_dev)
    df["log_chamfer"] = np.log(df.chamfer_div)
    return df


TERMS = ["length_within", "length_between", "terminal",
         "motif_rmsd_within", "motif_rmsd_between"]


def fit(df: pd.DataFrame, response: str, terms: list[str] | None = None):
    """Mixed model with a random intercept per design region."""
    import statsmodels.formula.api as smf

    terms = TERMS if terms is None else terms
    sub = df.dropna(subset=[response, *terms])
    model = smf.mixedlm(f"{response} ~ " + " + ".join(terms), sub,
                        groups=sub["region_id"].to_numpy())
    # lbfgs hits a singular Hessian here. powell, bfgs, nm and cg all converge
    # and agree on the variance components to three decimals.
    return model.fit(reml=True, method="powell")


def varying_length(df: pd.DataFrame) -> pd.DataFrame:
    """Rows from the regions whose contig actually allows a length range.

    Four of the twelve regions are pinned to one length, so they carry no
    information about size and are dropped before anything is asked about it.
    """
    return df[df.groupby("region_id").length.transform("nunique") > 1]


def fit_random_slope(df: pd.DataFrame, response: str):
    """Same model but with the size slope free to vary between regions.

    A common slope assumes every region responds to extra residues the same way.
    Testing that against a random slope is what says whether a pooled size
    coefficient means anything.
    """
    import statsmodels.formula.api as smf

    sub = varying_length(df).dropna(subset=[response, *TERMS])
    return smf.mixedlm(f"{response} ~ " + " + ".join(TERMS), sub,
                       groups=sub["region_id"].to_numpy(),
                       re_formula="~length_within").fit(reml=False, method="powell")


def slope_lrt(df: pd.DataFrame, response: str) -> dict:
    """Likelihood ratio test of the random slope against a common slope."""
    import statsmodels.formula.api as smf
    from scipy.stats import chi2

    sub = varying_length(df).dropna(subset=[response, *TERMS])
    formula = f"{response} ~ " + " + ".join(TERMS)
    g = sub["region_id"].to_numpy()
    m0 = smf.mixedlm(formula, sub, groups=g).fit(reml=False, method="powell")
    m1 = smf.mixedlm(formula, sub, groups=g,
                     re_formula="~length_within").fit(reml=False, method="powell")
    stat = 2 * (m1.llf - m0.llf)
    return {"chi2": float(stat), "df": 2, "p": float(chi2.sf(stat, 2)),
            "aic_common": float(m0.aic), "aic_varying": float(m1.aic),
            "beta_common": float(m0.params["length_within"]),
            "beta_varying": float(m1.params["length_within"])}


def per_region_slopes(df: pd.DataFrame, response: str) -> pd.DataFrame:
    """Size slope fitted separately inside each region, with a sign test.

    The pooled coefficient can be positive while most regions are flat or
    negative, so the individual slopes are what show whether the effect is a
    general one or a couple of regions carrying the average.
    """
    from scipy.stats import binomtest, linregress

    rows = []
    for rid, g in varying_length(df).groupby("region_id"):
        lr = linregress(g.length_within, g[response])
        rows.append({"region_id": rid, "run": g.run.iloc[0],
                     "region": int(g.region.iloc[0]), "n": len(g),
                     "n_lengths": int(g.length.nunique()), "slope": lr.slope,
                     "se": lr.stderr, "p": lr.pvalue})
    out = pd.DataFrame(rows)
    n_pos = int((out.slope > 0).sum())
    out.attrs["sign_test_p"] = float(binomtest(n_pos, len(out)).pvalue)
    out.attrs["n_positive"] = n_pos
    return out


def permutation_p(df: pd.DataFrame, response: str) -> tuple[float, int]:
    """Exact test for the `terminal` effect, permuting at the region level.

    `terminal` is fixed within a region, so 768 design rows carry no more
    information about it than the 12 regions do and the Wald p-value is
    anticonservative. Each run has exactly one terminal region by construction,
    so the null permutes which region in each run is labelled terminal. That is
    2 x 4 x 6 = 48 assignments, small enough to enumerate exactly.
    """
    from itertools import product

    obs = fit(df, response).params["terminal"]
    per_run = {r: sorted(g.region.unique()) for r, g in df.groupby("run")}
    runs = list(per_run)

    null = []
    for combo in product(*(per_run[r] for r in runs)):
        pick = dict(zip(runs, combo))
        d = df.assign(terminal=[int(pick[r] == reg)
                                for r, reg in zip(df.run, df.region)])
        null.append(fit(d, response).params["terminal"])

    null = np.asarray(null)
    p = float((np.abs(null) >= abs(obs)).mean())
    return p, len(null)


def coef_table(res, response: str) -> pd.DataFrame:
    ci = res.conf_int()
    out = pd.DataFrame({
        "response": response,
        "term": res.params.index,
        "estimate": res.params.values,
        "lo": ci[0].values,
        "hi": ci[1].values,
        "p": res.pvalues.values,
    })
    return out[~out.term.str.contains("Group Var")].reset_index(drop=True)
