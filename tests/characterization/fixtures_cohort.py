"""Build negative-steering cohort CSVs whose figures' content is known before
the plot code runs.

Driving the plot scripts with a copy of a real run only proves they still
agree with their own last output. These builders go the other way: each one
takes the answer — this seed passed, that one collapsed on reversion, this
sequence is tier A — and writes the CSV row that must produce it. A figure
that then disagrees is a failure rather than a new blessed reference.

Three row builders, one per file the plot scripts read:

`cross_row` builds a cross_sequence_summary.csv row (the cohort's tier and
composite rank).

`agg_row` builds an aggregated_results.csv row (one sequence_group's medians,
with a steered and a reverted side so "which prediction is final" is testable).

`seed_row` builds a raw_per_seed_results.csv row from the (stage, outcome) pair
that `classify_seed` must return for it. The metric values are chosen to sit
on the correct side of every threshold in extract_passing, so the expected
classification is a property of the call site, not of the numbers.

`passing_row` builds a passing_summary.csv row, the per-sequence file the
cohort aggregation reads. Its column set comes from
`extract_passing.OUTPUT_FIELDS`, which is the code that writes the file.

`cohort_tree` assembles the first three into the runs/ layout the two scripts
walk, and `write_passing_summary` writes one sequence's passing_summary.csv
for the aggregation to read.

Column names are taken from the committed reference set under
tests/negative_steering/example_output_files/, never from docstrings. Each
builder rejects an override naming a column it does not know, so a typo
cannot quietly invent a column that the production code would never read.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List

_BIN = Path(__file__).resolve().parents[2] / "bin"
if str(_BIN) not in sys.path:
    sys.path.insert(0, str(_BIN))

from extract_passing import OUTPUT_FIELDS as PASSING_SUMMARY_COLUMNS  # noqa: E402

# ── cross_sequence_summary.csv ───────────────────────────────────────────
# Subset of the 63-column header in
# example_output_files/negative_steering/cross_sequence_summary.csv, limited
# to the columns bin/negsteer_plots.py and bin/negsteer_within_sequence_plots.py
# actually read.
_CROSS_DEFAULTS: Dict[str, str] = {
    "mpnn_sequence": "",
    "row_type": "steered",
    "cross_rank_by_composite": "",
    "cross_tier": "none",
    "cross_composite_score": "",
    "rep_cycle": "0",
    "rep_pathway": "cycle_0",
    "rep_sequence_group": "0",
    "rep_n_seeds": "3",
    "rep_n_pass": "0",
    "rep_outcome": "no_reversion",
    "rep_total_mutations_median": "0",
    "rep_steered_mutations_chimerax": "",
    "rep_steered_mutations_aa": "",
    "rep_reverted_mutations_majority_chimerax": "",
    "rep_ra_eff_vs_truth_median": "",
    "rep_true_jaccard_median": "",
    "rep_weighted_jaccard_median": "",
    "rep_complex_plddt_median": "",
    "rep_avg_plddt_median": "",
    "rep_interface_plddt_median": "",
    "rep_ipae_median": "",
    "rep_pae_pass_frac_median": "",
    "rep_iptm_median": "",
    "rep_ptm_median": "",
    "rep_actifptm_median": "",
    "rep_ipsae_min_median": "",
    "rep_receptor_intact_majority": "",
}

# ── aggregated_results.csv ───────────────────────────────────────────────
# Subset of the 280-column header in
# example_output_files/.../runs/<seq>/aggregated_results.csv.
_AGG_DEFAULTS: Dict[str, str] = {
    "cycle": "0",
    "pathway": "cycle_0",
    "sequence_group": "0",
    "design": "design_00",
    "n_seeds": "3",
    "n_pass": "0",
    "outcome": "no_reversion",
    "steered_total_mutations_median": "0",
    "steered_mutations_chimerax": "",
    "steered_mutations_aa": "",
    "steered_ra_eff_vs_truth_median": "",
    "steered_true_jaccard_median": "",
    "steered_weighted_jaccard_median": "",
    "steered_complex_plddt_median": "",
    "steered_interface_plddt_median": "",
    "steered_ipae_median": "",
    "steered_pae_pass_frac_median": "",
    "steered_iptm_median": "",
    "steered_receptor_intact_majority": "",
    "reverted_ra_eff_vs_truth_median": "",
    "reverted_true_jaccard_median": "",
    "reverted_weighted_jaccard_median": "",
    "reverted_complex_plddt_median": "",
    "reverted_interface_plddt_median": "",
    "reverted_ipae_median": "",
    "reverted_pae_pass_frac_median": "",
    "reverted_iptm_median": "",
    "reverted_receptor_intact_majority": "",
}

# ── raw_per_seed_results.csv ─────────────────────────────────────────────
# Subset of the 87-column header in
# example_output_files/.../runs/<seq>/raw_per_seed_results.csv.
_SEED_DEFAULTS: Dict[str, str] = {
    "cycle": "0",
    "pathway": "cycle_0",
    "design": "design_00",
    "sequence_group": "0",
    "seed_index": "0",
    "status": "ok",
    "steered_total_mutations": "",
    "steered_mutations_chimerax": "",
    "steered_mutations_aa": "",
    "steered_ra_eff_vs_truth": "",
    "steered_receptor_intact": "",
    "steered_true_jaccard": "",
    "steered_weighted_jaccard": "",
    "steered_avg_plddt": "",
    "steered_complex_plddt": "",
    "steered_interface_plddt": "",
    "steered_ptm": "",
    "steered_iptm": "",
    "steered_actifptm": "",
    "steered_pae_mean": "",
    "steered_ipae": "",
    "steered_pae_pass_frac": "",
    "steered_ipsae_ab": "",
    "steered_ipsae_ba": "",
    "steered_ipsae_min": "",
    "steered_ipsae_ab_15": "",
    "steered_ipsae_ba_15": "",
    "steered_ipsae_min_15": "",
    "reversion_verdict": "",
    "reverted_total_mutations": "",
    "reverted_mutations_chimerax": "",
    "reverted_ra_eff_vs_truth": "",
    "reverted_receptor_intact": "",
    "reverted_true_jaccard": "",
    "reverted_weighted_jaccard": "",
    "reverted_avg_plddt": "",
    "reverted_complex_plddt": "",
    "reverted_interface_plddt": "",
    "reverted_ptm": "",
    "reverted_iptm": "",
    "reverted_actifptm": "",
    "reverted_pae_mean": "",
    "reverted_ipae": "",
    "reverted_pae_pass_frac": "",
    "reverted_ipsae_ab": "",
    "reverted_ipsae_ba": "",
    "reverted_ipsae_min": "",
    "reverted_ipsae_ab_15": "",
    "reverted_ipsae_ba_15": "",
    "reverted_ipsae_min_15": "",
}

# ── passing_summary.csv ──────────────────────────────────────────────────
# The full header, taken from the module that writes the file rather than
# restated here. The populated subset below is what decides a sequence's tier
# and composite rank in the cross-sequence aggregation:
#
#   n_pass == n_seeds       tier A, so 3 of 3 is the default
#   ra_eff / true_jaccard   the composite score, true_jaccard - 0.05 * ra_eff
#   outcome                 which prediction counts as the sequence's final
#
# Every other column is blank, which is what the writer emits for a metric a
# run did not produce.
_PASSING_DEFAULTS: Dict[str, str] = {col: "" for col in PASSING_SUMMARY_COLUMNS}
_PASSING_DEFAULTS.update({
    "rank_by_composite_score": "1",
    "rank_by_ra_eff": "1",
    "cycle": "0",
    "pathway": "cycle_0",
    "sequence_group": "0",
    "design": "initial",
    "n_seeds": "3",
    "n_pass": "3",
    "outcome": "no_reversion",
    "ra_eff_vs_truth_median": "2.0",
    "true_jaccard_median": "0.80",
    "confidence_flag": "ok",
})

# Metric blocks, one per outcome `classify_seed` can return. Each sits on the
# stated side of every threshold in extract_passing (ra_eff < 5 Å, receptor
# intact, complex_plddt >= 0.70, ipae <= 15, pae_pass_frac >= 0.10,
# iptm >= 0.30), so the outcome named at the call site is what the classifier
# must produce.
#
#   pass               structure and confidence both clear
#   off_target         ra_eff 30 Å, confidence clears
#   poor_prediction    ra_eff 2 Å, every confidence metric fails
#   multiple_failures  both sides fail
#   no_data            ra_eff is the literal "nan" a prediction that never
#                      landed writes, which must read as absent, not as a
#                      failure
#   new_contamination  the engine's own label; metrics are never consulted
_OUTCOME_METRICS: Dict[str, Dict[str, object]] = {
    "pass":              {"ra": 2.0, "intact": 1, "plddt": 0.92, "ipae": 4.0,
                          "paepf": 0.80, "iptm": 0.72, "tj": 0.70, "wj": 0.66},
    "off_target":        {"ra": 30.0, "intact": 1, "plddt": 0.92, "ipae": 4.0,
                          "paepf": 0.80, "iptm": 0.72, "tj": 0.08, "wj": 0.05},
    "poor_prediction":   {"ra": 2.0, "intact": 1, "plddt": 0.20, "ipae": 40.0,
                          "paepf": 0.01, "iptm": 0.10, "tj": 0.60, "wj": 0.55},
    "multiple_failures": {"ra": 30.0, "intact": 0, "plddt": 0.20, "ipae": 40.0,
                          "paepf": 0.01, "iptm": 0.10, "tj": 0.04, "wj": 0.02},
    "no_data":           {"ra": "nan", "intact": 1, "plddt": 0.92, "ipae": 4.0,
                          "paepf": 0.80, "iptm": 0.72, "tj": 0.70, "wj": 0.66},
    "new_contamination": {"ra": 30.0, "intact": 0, "plddt": 0.20, "ipae": 40.0,
                          "paepf": 0.01, "iptm": 0.10, "tj": 0.03, "wj": 0.01},
}

# What the steered side of a reversion-stage seed says. Deliberately a healthy
# prediction, and deliberately distinct from every `ra` in _OUTCOME_METRICS, so
# a test can tell whether a reader took the final (reverted) value or the
# intermediate (steered) one.
STEERED_BEFORE_REVERSION = {"ra": 2.5, "tj": 0.75, "wj": 0.70}

ALL_OUTCOMES = tuple(_OUTCOME_METRICS)


def _checked(defaults: Dict[str, str], overrides: Dict[str, object],
             what: str) -> Dict[str, str]:
    unknown = sorted(set(overrides) - set(defaults))
    if unknown:
        raise KeyError(
            f"{what}: unknown column(s) {unknown}. Known columns come from the "
            f"committed header under tests/negative_steering/"
            f"example_output_files/; add the column there first."
        )
    row = dict(defaults)
    row.update({k: "" if v is None else str(v) for k, v in overrides.items()})
    return row


def write_csv(path: Path, rows: Iterable[Dict[str, str]],
              columns: Iterable[str]) -> Path:
    """Write `rows` with `columns` as the header, creating parent dirs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(columns)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in cols})
    return path


# ── Row builders ─────────────────────────────────────────────────────────

def cross_row(name: str, **overrides) -> Dict[str, str]:
    """One cross_sequence_summary.csv row for MPNN sequence `name`."""
    overrides.setdefault("mpnn_sequence", name)
    return _checked(_CROSS_DEFAULTS, overrides, "cross_row")


def passing_row(**overrides) -> Dict[str, str]:
    """One passing_summary.csv row, every OUTPUT_FIELDS column present.

    The default is a tier-A sequence: 3 of 3 seeds passed, ra_eff 2.0 A and
    true Jaccard 0.80, so a test that wants a different tier says so by
    naming `n_pass`.
    """
    return _checked(_PASSING_DEFAULTS, overrides, "passing_row")


def write_passing_summary(path: Path, rows: List[Dict[str, str]]) -> Path:
    """One sequence's passing_summary.csv at `path`, parents created."""
    return write_csv(path, rows, PASSING_SUMMARY_COLUMNS)


def agg_row(**overrides) -> Dict[str, str]:
    """One aggregated_results.csv row."""
    return _checked(_AGG_DEFAULTS, overrides, "agg_row")


def agg_row_for(outcome: str, *, steered_ra: float, steered_tj: float,
                reverted_ra: float | None = None,
                reverted_tj: float | None = None,
                **overrides) -> Dict[str, str]:
    """An aggregated_results.csv row with both metric sides populated.

    `outcome` is the engine's verdict, which is what decides whether the
    reverted or the steered side is the sequence's final prediction. The two
    sides are given different numbers so a reader taking the wrong one is
    visible rather than coincidentally right.
    """
    row = {
        "outcome": outcome,
        "steered_ra_eff_vs_truth_median": steered_ra,
        "steered_true_jaccard_median": steered_tj,
        "steered_complex_plddt_median": 0.92,
        "steered_ipae_median": 4.0,
        "steered_pae_pass_frac_median": 0.80,
        "steered_iptm_median": 0.72,
        "steered_interface_plddt_median": 0.80,
        "steered_weighted_jaccard_median": steered_tj,
        "steered_receptor_intact_majority": 1,
    }
    if reverted_ra is not None:
        row["reverted_ra_eff_vs_truth_median"] = reverted_ra
        row["reverted_complex_plddt_median"] = 0.55
        row["reverted_ipae_median"] = 22.0
        row["reverted_pae_pass_frac_median"] = 0.02
        row["reverted_iptm_median"] = 0.12
        row["reverted_interface_plddt_median"] = 0.40
        row["reverted_receptor_intact_majority"] = 0
    if reverted_tj is not None:
        row["reverted_true_jaccard_median"] = reverted_tj
        row["reverted_weighted_jaccard_median"] = reverted_tj
    row.update(overrides)
    return agg_row(**row)


def seed_row(*, stage: str, outcome: str, seed_index: int = 0,
             sequence_group: str = "0", **overrides) -> Dict[str, str]:
    """One raw_per_seed_results.csv row that `classify_seed` must sort into
    exactly `(stage, outcome)`.

    stage == "reversion" fills the reverted_* columns from `outcome` and the
    steered_* columns from STEERED_BEFORE_REVERSION, because a reverted seed
    has both predictions and the reverted one is the last made.

    stage == "no_data" leaves steered_total_mutations unparseable, which is the
    only way the classifier cannot tell which stage a seed ended at.
    """
    if outcome not in _OUTCOME_METRICS:
        raise KeyError(f"seed_row: unknown outcome {outcome!r}")
    if stage not in ("cold_start", "steering", "reversion", "no_data"):
        raise KeyError(f"seed_row: unknown stage {stage!r}")
    m = _OUTCOME_METRICS[outcome]

    row: Dict[str, object] = {
        "sequence_group": sequence_group,
        "seed_index": seed_index,
    }
    if stage == "no_data":
        row["steered_total_mutations"] = ""
        row.update(_metric_block("steered", m))
    elif stage == "reversion":
        row["steered_total_mutations"] = 3
        row["steered_mutations_chimerax"] = "/A:5,7,26"
        row["steered_mutations_aa"] = "A5T,K7E,L26M"
        row.update(_metric_block("steered", {**_OUTCOME_METRICS["pass"],
                                             **STEERED_BEFORE_REVERSION}))
        row["reversion_verdict"] = (
            "new_contamination" if outcome == "new_contamination" else "pose_holds"
        )
        row["reverted_total_mutations"] = 2
        row["reverted_mutations_chimerax"] = "/A:7,26"
        row.update(_metric_block("reverted", m))
    else:
        row["steered_total_mutations"] = 0 if stage == "cold_start" else 3
        if stage == "steering":
            row["steered_mutations_chimerax"] = "/A:5,7,26"
            row["steered_mutations_aa"] = "A5T,K7E,L26M"
        row.update(_metric_block("steered", m))

    row.update(overrides)
    return _checked(_SEED_DEFAULTS, row, "seed_row")


def _metric_block(prefix: str, m: Dict[str, object]) -> Dict[str, object]:
    """The per-seed metric columns of one stage's column family."""
    return {
        f"{prefix}_ra_eff_vs_truth": m["ra"],
        f"{prefix}_receptor_intact": m["intact"],
        f"{prefix}_true_jaccard": m["tj"],
        f"{prefix}_weighted_jaccard": m["wj"],
        f"{prefix}_complex_plddt": m["plddt"],
        f"{prefix}_avg_plddt": float(m["plddt"]) * 100 if m["plddt"] != "" else "",
        f"{prefix}_interface_plddt": m["plddt"],
        f"{prefix}_ptm": 0.80,
        f"{prefix}_iptm": m["iptm"],
        f"{prefix}_actifptm": 0.65,
        f"{prefix}_pae_mean": 8.5,
        f"{prefix}_ipae": m["ipae"],
        f"{prefix}_pae_pass_frac": m["paepf"],
        f"{prefix}_ipsae_min": 0.44,
        f"{prefix}_ipsae_ab": 0.46,
        f"{prefix}_ipsae_ba": 0.42,
        f"{prefix}_ipsae_min_15": 0.48,
        f"{prefix}_ipsae_ab_15": 0.50,
        f"{prefix}_ipsae_ba_15": 0.46,
    }


def cold_start_baseline_row(ra_eff: float) -> Dict[str, str]:
    """The cycle-0 initial row: no sequence_group, no seed_index.

    Both scripts treat a blank sequence_group as the single pre-steering
    prediction. The cohort script drops it; the within-sequence script keeps
    its ra_eff as the trajectory plot's cold-start point.
    """
    return _checked(_SEED_DEFAULTS, {
        "sequence_group": "",
        "seed_index": "",
        "status": "initial",
        "steered_total_mutations": 0,
        "steered_ra_eff_vs_truth": ra_eff,
        "steered_receptor_intact": 1,
    }, "cold_start_baseline_row")


# ── Tree builders ────────────────────────────────────────────────────────

def write_cross_summary(path: Path, rows: List[Dict[str, str]]) -> Path:
    return write_csv(path, rows, _CROSS_DEFAULTS)


def write_sequence_dir(runs_dir: Path, name: str, *,
                       agg_rows: List[Dict[str, str]] | None = None,
                       seed_rows: List[Dict[str, str]] | None = None) -> Path:
    """One runs/<mpnn_sequence>/ directory.

    Omitting agg_rows or seed_rows omits that file, which is how a partial run
    looks on disk and is a path both scripts have to survive.
    """
    d = runs_dir / name
    d.mkdir(parents=True, exist_ok=True)
    if agg_rows is not None:
        write_csv(d / "aggregated_results.csv", agg_rows, _AGG_DEFAULTS)
    if seed_rows is not None:
        write_csv(d / "raw_per_seed_results.csv", seed_rows, _SEED_DEFAULTS)
    return d


def three_seeds(*, stage: str, outcome: str, sequence_group: str = "0",
                **overrides) -> List[Dict[str, str]]:
    """The three Boltz seeds of one sequence_group, all with the same fate."""
    return [seed_row(stage=stage, outcome=outcome, seed_index=i,
                     sequence_group=sequence_group, **overrides)
            for i in range(3)]


def cohort_tree(root: Path) -> Dict[str, Path]:
    """The cohort both plot scripts read, with every answer stated here.

    Four sequences and two controls:

      design_0_seq_0  tier A, rank 1, no_reversion. Final = the steered
                      prediction: ra_eff 2.0 Å, Jaccard 0.70, composite 0.60.
                      Three passing seeds in sg 0 plus three in sg 1.
      design_0_seq_1  tier B, rank 2, pose_collapses. Reversion ran, so final
                      = the reverted prediction: ra_eff 30.5 Å, Jaccard 0.05.
                      Its steered side still reads 2.0 Å / 0.70, which is the
                      value a reader must NOT take.
      design_1_seq_0  tier none, no rank. Present in runs/ and in the summary,
                      so its composite has to be recomputed from the medians.
      design_2_seq_0  in runs/ only, absent from the summary: the synthesised
                      path.
      input_control_polyA / input_control_scrambled
                      controls, ra_eff 25 Å and 22 Å, far above the 5 Å line
                      where the diagnostic plot says controls belong.

    Returns the paths the two CLIs take as arguments.
    """
    runs = root / "runs"
    cross = root / "cross_sequence_summary.csv"

    write_cross_summary(cross, [
        cross_row("design_0_seq_0", cross_tier="A", cross_rank_by_composite=1,
                  cross_composite_score=0.60, rep_sequence_group="0",
                  rep_outcome="no_reversion", rep_n_pass=3,
                  rep_ra_eff_vs_truth_median=2.0, rep_true_jaccard_median=0.70,
                  rep_weighted_jaccard_median=0.66,
                  rep_complex_plddt_median=0.92, rep_ipae_median=4.0,
                  rep_pae_pass_frac_median=0.80, rep_iptm_median=0.72,
                  rep_steered_mutations_chimerax="/A:5,7,26"),
        cross_row("design_0_seq_1", cross_tier="B", cross_rank_by_composite=2,
                  cross_composite_score=0.40, rep_sequence_group="0",
                  rep_outcome="pose_collapses", rep_n_pass=1,
                  rep_ra_eff_vs_truth_median=2.0, rep_true_jaccard_median=0.70,
                  rep_complex_plddt_median=0.92, rep_ipae_median=4.0,
                  rep_pae_pass_frac_median=0.80, rep_iptm_median=0.72,
                  rep_steered_mutations_chimerax="/A:5,7,26",
                  rep_reverted_mutations_majority_chimerax="/A:7,26"),
        cross_row("design_1_seq_0", cross_tier="none",
                  rep_sequence_group="", rep_outcome="no_reversion",
                  rep_ra_eff_vs_truth_median=12.0,
                  rep_true_jaccard_median=0.10,
                  rep_complex_plddt_median=0.40, rep_ipae_median=30.0,
                  rep_pae_pass_frac_median=0.01, rep_iptm_median=0.10,
                  rep_steered_mutations_chimerax="/A:31"),
        cross_row("input_control_polyA", row_type="control_polyA",
                  rep_ra_eff_vs_truth_median=25.0,
                  rep_true_jaccard_median=0.01),
        cross_row("input_control_scrambled", row_type="control_scrambled",
                  rep_ra_eff_vs_truth_median=22.0,
                  rep_true_jaccard_median=0.02),
    ])

    write_sequence_dir(
        runs, "design_0_seq_0",
        agg_rows=[agg_row_for("no_reversion", steered_ra=2.0, steered_tj=0.70,
                              sequence_group="0", n_pass=3,
                              steered_mutations_chimerax="/A:5,7,26"),
                  agg_row_for("no_reversion", steered_ra=8.0, steered_tj=0.20,
                              sequence_group="1", n_pass=0)],
        seed_rows=[cold_start_baseline_row(30.1)]
        + three_seeds(stage="steering", outcome="pass", sequence_group="0")
        + three_seeds(stage="steering", outcome="off_target", sequence_group="1"),
    )
    write_sequence_dir(
        runs, "design_0_seq_1",
        agg_rows=[agg_row_for("pose_collapses", steered_ra=2.0, steered_tj=0.70,
                              reverted_ra=30.5, reverted_tj=0.05,
                              sequence_group="0", n_pass=1,
                              steered_mutations_chimerax="/A:5,7,26")],
        seed_rows=[cold_start_baseline_row(29.4)]
        + three_seeds(stage="reversion", outcome="off_target", sequence_group="0"),
    )
    write_sequence_dir(
        runs, "design_1_seq_0",
        agg_rows=[agg_row_for("no_reversion", steered_ra=12.0, steered_tj=0.10,
                              sequence_group="0"),
                  agg_row(outcome="singleton", sequence_group="2")],
        seed_rows=[cold_start_baseline_row(31.0)]
        + three_seeds(stage="steering", outcome="multiple_failures",
                      sequence_group="0"),
    )
    write_sequence_dir(
        runs, "design_2_seq_0",
        agg_rows=[agg_row_for("no_reversion", steered_ra=4.0, steered_tj=0.50,
                              sequence_group="0", n_pass=2)],
        seed_rows=three_seeds(stage="cold_start", outcome="pass",
                              sequence_group="0"),
    )
    for ctrl, ra in (("input_control_polyA", 25.0),
                     ("input_control_scrambled", 22.0)):
        write_sequence_dir(
            runs, ctrl,
            agg_rows=[agg_row_for("no_reversion", steered_ra=ra, steered_tj=0.01,
                                  sequence_group="0")],
            seed_rows=three_seeds(stage="steering", outcome="off_target"),
        )

    design_region = root / "input_design_region.txt"
    design_region.write_text(
        "# derived from input structure input_complex.pdb + contigs\n"
        "# contigs: A1-32/10-20/A46-72/6-6 C\n"
        "33,34,35,36,37\n"
    )
    return {"cross_csv": cross, "runs_dir": runs,
            "design_region": design_region}
