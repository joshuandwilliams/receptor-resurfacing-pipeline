"""
Shared helpers for the correlates-of-success analyses.

Imported by each analysis QMD via:
    import sys; sys.path.insert(0, ".")
    from helpers import *
"""

import json
import random
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from lightgbm import LGBMClassifier, LGBMRegressor
from scipy.stats import pearsonr, pointbiserialr
from sklearn.model_selection import GroupShuffleSplit, cross_val_score
from IPython.display import display

warnings.filterwarnings("ignore")

# ── Reproducibility ───────────────────────────────────────────────────────────
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

# ── Paths ─────────────────────────────────────────────────────────────────────
HPC_BASE = Path(
    "/Volumes/HPC-Home/receptor_design/receptor-resurfacing-pipeline"
    "/experiments/campaigns/pikp1_avrpikf/runs"
)
RUNS = {
    "crystal": "crystal_full_test_contig",
    # posed run available but excluded from primary analysis:
    # only 13 hits gives near-chance CV AUC and unreliable SHAP values.
}
OUT_DIR = Path("plots")
OUT_DIR.mkdir(exist_ok=True)

# ── Amino-acid sets ───────────────────────────────────────────────────────────
HYDROPHOBIC = set("VILMFYW")
POSITIVE    = set("KR")
NEGATIVE    = set("DE")
AROMATIC    = set("FYW")


def aa_features(seq: str) -> dict:
    """Scalar composition features from a single amino-acid string."""
    if not seq or pd.isna(seq):
        return dict(hydrophob_frac=np.nan, net_charge=np.nan,
                    has_proline=np.nan, aromatic_frac=np.nan)
    n = len(seq)
    return dict(
        hydrophob_frac = sum(c in HYDROPHOBIC for c in seq) / n,
        net_charge     = sum(c in POSITIVE for c in seq) - sum(c in NEGATIVE for c in seq),
        has_proline    = int("P" in seq),
        aromatic_frac  = sum(c in AROMATIC for c in seq) / n,
    )


# ── Data loading ──────────────────────────────────────────────────────────────

def load_run(label: str, run_dir_name: str) -> pd.DataFrame:
    """
    Load and merge all per-design data for one run.
    Returns one row per (design_idx, seq_idx) — 128 steered rows.

    Columns added beyond the raw source files
    -----------------------------------------
    dG_separated_win   : dG_separated winsorised at 1st/99th percentile.
                         Raw max is +1088 REU (severe Rosetta clashes).
    loop1_len          : Length of de novo segment 1, from
                         design_region_length_observed. This is what
                         RFDiffusion chose within the 10-20 constraint.
    mpnn_dr_ratio      : design_region_score / mpnn_score. >1 means the
                         de novo positions are harder to design than the
                         fixed scaffold.
    rep_stage          : Pipeline stage of the representative design.
                         "initial"       — initial Boltz prediction selected
                                           as representative; 100% hit rate
                                           in this dataset by construction.
                         "steered"       — a steering-cycle design is the
                                           representative; outcome=no_reversion.
                         "pose_holds"    — steering ran, reversion triggered,
                                           mutations confirmed genuine.
                         "pose_collapses"— reversion triggered, complex
                                           collapsed; non-specific binding.
    hit                : 1 if cross_tier in {A, B, C}, else 0.
    run / run_bin      : run label and binary (crystal=1, posed=0).
    """
    base = HPC_BASE / run_dir_name / "results"

    # Negsteer summary
    neg = pd.read_csv(base / "negative_steering" / "cross_sequence_summary.csv")
    neg = neg[neg["row_type"] == "steered"].copy()
    neg["design_idx"] = (
        neg["source_passing_summary"].str.extract(r"design_(\d+)_seq")[0].astype(float)
    )
    neg["seq_idx"] = (
        neg["source_passing_summary"].str.extract(r"seq_(\d+)/")[0].astype(float)
    )

    # MPNN / sequence metadata
    meta = pd.read_csv(base / "sequences" / "scored_metadata.csv")

    # Rosetta filter metrics
    rf_raw = json.load(open(base / "rosetta_filtering" / "rosetta_filter_metrics.json"))
    rf = pd.DataFrame(rf_raw["designs"])
    rf["design_idx"] = rf["design_stem"].str.extract(r"design_(\d+)")[0].astype(float)
    rf = rf[["design_idx", "sc_value", "dG_separated", "dSASA_int",
             "dG_dSASA_density", "packstat", "delta_unsatHbonds", "nres_int"]]

    # Merge — negsteer and meta share column names; keep meta versions
    df = (
        neg.merge(meta, left_on=["design_idx", "seq_idx"],
                  right_on=["design", "seq"], how="left", suffixes=("_neg", ""))
           .merge(rf, on="design_idx", how="left")
    )
    df.drop(columns=[c for c in df.columns if c.endswith("_neg")], inplace=True)

    # ── Derived features ───────────────────────────────────────────────────

    lo, hi = df["dG_separated"].quantile([0.01, 0.99])
    df["dG_separated_win"] = df["dG_separated"].clip(lo, hi)

    df["loop1_len"]        = df["design_region_length_observed"].str.split("|").str[0].astype(float)
    df["num_changes_loop"] = df["num_changes"].astype(str).str.split("|").str[0].astype(float)
    df["num_changes_tail"] = df["num_changes"].astype(str).str.split("|").str[1].astype(float)

    loop_seqs = df["designed_residues"].str.split("|").str[0]
    loop_feats = loop_seqs.apply(aa_features).apply(pd.Series)
    loop_feats.columns = [f"loop_{c}" for c in loop_feats.columns]
    df = pd.concat([df, loop_feats], axis=1)

    df["mpnn_dr_ratio"] = df["design_region_score"] / df["mpnn_score"].replace(0, np.nan)

    _reversion_outcomes = {"pose_holds", "pose_collapses", "new_contamination"}

    def _stage(row):
        if row["rep_design"] == "initial":
            return "initial"
        if row["rep_outcome"] in _reversion_outcomes:
            return "reversion"
        return "steered"

    df["rep_stage"] = df.apply(_stage, axis=1)

    df["hit"] = df["cross_tier"].isin(["A", "B", "C"]).astype(int)
    df["run"] = label
    return df


# ── Feature lists ─────────────────────────────────────────────────────────────
#
# Stage colours for SHAP plots
# red    = RFD / Rosetta  (computed on RFD backbone with Rosetta sidechains)
# blue   = MPNN sequence  (computed on the designed sequence)
# green  = Boltz structural (computed from Boltz complex prediction)
# yellow = contact geometry (also from Boltz, but contact-specific)
# purple = steering dynamics

FEATS_RFD = [
    "sc_value",           # Rosetta shape complementarity (backbone geometry)
    "dG_separated_win",   # Rosetta binding energy, winsorised
    "dSASA_int",          # Buried interface area (Å²)
    "dG_dSASA_density",   # Binding energy per unit buried area
    "packstat",           # Interface packing quality
    "delta_unsatHbonds",  # Unsatisfied H-bonds at interface
    "nres_int",           # Number of interface residues
    "loop1_len",          # RFDiffusion-chosen loop length (10-20 range)
]

FEATS_MPNN = [
    "mpnn_score",          # Global MPNN sequence log-probability
    "design_region_score", # MPNN log-probability on de novo positions only
    "mpnn_dr_ratio",       # design_region_score / mpnn_score
    # num_changes_loop is excluded: r ≈ 1 with loop1_len (12% recovery rate)
    "loop_hydrophob_frac", # Fraction hydrophobic in designed loop
    "loop_net_charge",     # Net charge of designed loop
    "loop_has_proline",    # Binary: proline present in designed loop
    "loop_aromatic_frac",  # Fraction aromatic in designed loop
]

FEATS_BOLTZ_STRUCTURAL = [
    "rep_independent_receptor_rmsd_median",  # Receptor RMSD alone vs in-complex
    "rep_independent_effector_rmsd_median",  # Binder RMSD alone vs in-complex
    "rep_iptm_median",            # Interface pTM score
    "rep_ptm_median",             # Overall pTM score
    "rep_interface_plddt_median", # Per-residue pLDDT at interface
    "rep_ipae_median",            # Interface PAE
    "rep_pae_pass_frac_median",   # Fraction of PAE pairs below threshold
]

FEATS_CONTACT = [
    "rep_n_contact_residues_median", # Total receptor residues contacted
    "rep_wrong_jaccard_median",      # Contact overlap with scrambled receptor
    "rep_n_shared_wrong_median",     # Raw count of shared-wrong contacts
    # rep_true_jaccard_median  EXCLUDED: ingredient of cross_composite_score
    # rep_n_shared_true_median EXCLUDED: numerator of true_jaccard
    # jaccard_specificity      EXCLUDED: derived from true_jaccard
    # rep_weighted_jaccard     EXCLUDED: composite of true/wrong jaccard
]

FEATS_STEERING = [
    "rep_total_mutations_median", # Mutations introduced by steering (0 if initial)
    "rep_cycle",                  # Steering cycle of representative
    # rep_stage            NOT a feature — explained in glossary
    # rep_n_pass           EXCLUDED: definitional (must be > 0 to hit)
    # rep_ra_eff_vs_truth  EXCLUDED: ingredient of cross_composite_score
    # rep_outcome_encoded  EXCLUDED: rep_stage captures this more precisely
]

STAGE_COLOURS = (
    {f: "#E57373" for f in FEATS_RFD}
    | {f: "#64B5F6" for f in FEATS_MPNN}
    | {f: "#81C784" for f in FEATS_BOLTZ_STRUCTURAL}
    | {f: "#FFD54F" for f in FEATS_CONTACT}
    | {f: "#BA68C8" for f in FEATS_STEERING}
)


# ── Modelling helpers ─────────────────────────────────────────────────────────

def make_X(df: pd.DataFrame, feature_list: list) -> pd.DataFrame:
    return df[feature_list].apply(pd.to_numeric, errors="coerce")


def backbone_cv(model, X: pd.DataFrame, y: pd.Series,
                groups: pd.Series, scoring: str, n_splits: int = 5):
    """
    Cross-validate with GroupShuffleSplit so that sibling sequences
    (same backbone, different MPNN sample) always land in the same fold.
    """
    gss = GroupShuffleSplit(n_splits=n_splits, test_size=0.2, random_state=SEED)
    scores = cross_val_score(model, X, y, cv=gss, groups=groups, scoring=scoring)
    return scores.mean(), scores.std()


def fit_and_shap(model, X: pd.DataFrame, y: pd.Series,
                 background_frac: float = 0.4):
    """Fit model, return (explainer, shap_values) using interventional SHAP."""
    model.fit(X, y)
    bg = X.sample(frac=background_frac, random_state=SEED)
    explainer = shap.TreeExplainer(
        model, data=bg, feature_perturbation="interventional"
    )
    sv = explainer.shap_values(X)
    if isinstance(sv, list):
        sv = sv[1]
    return explainer, sv


def lgbm_classifier(**kwargs):
    defaults = dict(n_estimators=200, max_depth=4, learning_rate=0.05,
                    subsample=0.8, colsample_bytree=0.8,
                    random_state=SEED, n_jobs=1, verbose=-1)
    return LGBMClassifier(**{**defaults, **kwargs})


def lgbm_regressor(**kwargs):
    defaults = dict(n_estimators=200, max_depth=4, learning_rate=0.05,
                    subsample=0.8, colsample_bytree=0.8,
                    random_state=SEED, n_jobs=1, verbose=-1)
    return LGBMRegressor(**{**defaults, **kwargs})


# ── Plotting helpers ──────────────────────────────────────────────────────────

def beeswarm(sv, X: pd.DataFrame, title: str,
             save_name: str = None, max_display: int = 20):
    expl = shap.Explanation(
        values=sv, base_values=np.zeros(len(sv)),
        data=X.values, feature_names=list(X.columns),
    )
    fig, ax = plt.subplots(figsize=(10, 7))
    shap.plots.beeswarm(expl, max_display=max_display, show=False, plot_size=None)
    ax = plt.gca()
    ax.set_title(title, fontsize=11, pad=10)
    for lbl in ax.get_yticklabels():
        lbl.set_color(STAGE_COLOURS.get(lbl.get_text(), "#888888"))
    plt.tight_layout()
    if save_name:
        plt.savefig(OUT_DIR / f"{save_name}.png", bbox_inches="tight")
    plt.show()


def dependence(sv, X: pd.DataFrame, feature: str,
               colour_feature: str, title: str, save_name: str = None):
    expl = shap.Explanation(
        values=sv, base_values=np.zeros(len(sv)),
        data=X.values, feature_names=list(X.columns),
    )
    fidx = list(X.columns).index(feature)
    fig, ax = plt.subplots(figsize=(7, 5))
    if colour_feature and colour_feature in X.columns:
        cidx = list(X.columns).index(colour_feature)
        shap.plots.scatter(expl[:, fidx], color=expl[:, cidx], show=False, ax=ax)
    else:
        shap.plots.scatter(expl[:, fidx], show=False, ax=ax)
    ax.set_title(title, fontsize=11)
    plt.tight_layout()
    if save_name:
        plt.savefig(OUT_DIR / f"{save_name}.png", bbox_inches="tight")
    plt.show()


def top_features(sv, feature_names: list, n: int = 3) -> list:
    return [feature_names[i]
            for i in np.argsort(np.abs(sv).mean(axis=0))[::-1][:n]]
