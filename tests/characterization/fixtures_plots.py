"""Inputs for the five plot CLIs, with the expected figure content at the call site.

Each of `bin/rfdiffusion_plots.py`, `bin/rosetta_filter_plots.py`,
`bin/mpnn_plots.py`, `bin/orthogonal_metrics_plots.py` and
`bin/pose_solver_plots.py` is a standalone CLI that reads one or more CSV/JSON
files and writes PNGs. A figure test is only worth anything if the right answer
is known before the code runs, so every builder here takes the values that go
into the rows rather than hiding them in a default. A test that wants design 7
to rank first says so in the call.

Column and key names are taken from the producers and the committed reference
sets, never invented:

    rfdiffusion   bin/rfdiffusion_filter.py, and
                  tests/rfdiffusion/example_output_files/rfdiffusion/
                  rfdiffusion_metrics.json
    rosetta       bin/rosetta_filter_collect.py, and
                  tests/rosetta_filtering/example_output_files/
                  rosetta_filtering/rosetta_filter_metrics.json
    mpnn          bin/mpnn_sequence_qc.py / bin/mpnn_design_region_score.py, and
                  tests/proteinmpnn/example_output_files/sequences/
                  scored_metadata.csv + mpnn_cluster_counts.csv
    orthogonal    bin/merge_orthogonal_metrics.py (which columns it appends) and
                  bin/cross_sequence_summary.py, with the committed headers in
                  tests/orthogonal_metrics/
    pose solver   bin/pose_solver.py and bin/pose_interface_metrics.py, with the
                  committed tests/pose_solver/*/solved_pose_results.json and
                  tests/pose_solver/interface_metrics_out/comparison.json

Two helpers go with them. `assert_png_not_blank` is the floor every figure test
needs, because a crashed plot and a blank canvas both leave a file behind.
`figures_kept` suspends the modules' `plt.close` so a test can read back the
tick labels and bar heights a plot function set, which is the only way to assert
WHICH row a plot ranked first when the ordering exists nowhere but inside the
figure.
"""
from __future__ import annotations

import contextlib
import csv
import json
from pathlib import Path
from typing import Iterable, Sequence

from tests.characterization.helpers.png_compare import count_distinct_colours

# ── Figure assertions ────────────────────────────────────────────────

def assert_png_not_blank(path: Path) -> None:
    """The file exists and has more than one pixel colour.

    Nextflow publishes a fixed figure list, so a stage that writes nothing fails
    with nothing in the log. A stage that writes a blank canvas is worse: it
    passes. Counting distinct colours separates "rendered something" from "wrote
    an empty axes", which is what a bare `exists()` misses. Shares
    `count_distinct_colours` with the negsteer figure tests so both tiers apply
    the same definition of blank.
    """
    assert path.exists(), f"no figure written to {path}"
    assert path.stat().st_size > 0, f"empty file at {path}"
    n = count_distinct_colours(path)
    assert n > 1, f"{path.name} is a single flat colour ({n} distinct)"


@contextlib.contextmanager
def figures_kept():
    """Suspend `plt.close` so the figure a plot function built stays readable.

    Every plot function here calls `savefig` then `plt.close`, which discards the
    artists. Row order, bar heights and tick labels then exist only inside the
    PNG. Disabling `close` leaves the figure registered with pyplot, so a test
    can read what the module set. Yields a callable returning the live figures in
    creation order; the one a plot function saved is the last, since the
    throwaway figures some helpers build (label measurement, dendrogram leaf
    order) are created first.
    """
    import matplotlib.pyplot as plt

    real_close = plt.close
    real_close("all")
    plt.close = lambda *a, **k: None          # noqa: E731 - deliberate stub
    try:
        yield lambda: [plt.figure(n) for n in plt.get_fignums()]
    finally:
        plt.close = real_close
        plt.close("all")


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


# ══════════════════════════════════════════════════════════════════════
# bin/rfdiffusion_plots.py
# ══════════════════════════════════════════════════════════════════════

# rfdiffusion_filter numbers de novo residues from a synthetic base well above
# any real receptor resnum, so they sort after every fixed residue. That is why
# per_region_lengths exists at all: sorting the design resnums cannot recover
# the region boundaries, and the sort-based fallback in
# _split_design_into_regions collapses every region into one.
DE_NOVO_BASE = 1000


def receptor_layout(
    *,
    fixed_runs: Sequence[int],
    region_lengths: Sequence[int],
    de_novo_base: int = DE_NOVO_BASE,
) -> tuple[list[int], list[int], list[int]]:
    """Return (receptor_position_order, design_resnums, fixed_resnums).

    `fixed_runs` and `region_lengths` interleave in file order: fixed_runs[0]
    residues, then region_lengths[0] de novo residues, then fixed_runs[1], and so
    on. Fixed residues are numbered 1.. in order. De novo residues are numbered
    from `de_novo_base`, which is the arrangement in the committed
    rfdiffusion_metrics.json (receptor_position_order there runs
    1..32, 1072..1090, 46..72, 1091..1096).

    So the position index of any residue is its index in the returned order, and
    a contact's in-region/out-of-region classification is decided purely by
    whether its resnum is in the returned design list.
    """
    order: list[int] = []
    design: list[int] = []
    fixed: list[int] = []
    next_fixed = 1
    next_dn = de_novo_base
    for i in range(max(len(fixed_runs), len(region_lengths))):
        if i < len(fixed_runs):
            for _ in range(fixed_runs[i]):
                order.append(next_fixed)
                fixed.append(next_fixed)
                next_fixed += 1
        if i < len(region_lengths):
            for _ in range(region_lengths[i]):
                order.append(next_dn)
                design.append(next_dn)
                next_dn += 1
    return order, design, fixed


def region_coords(region_lengths: Sequence[int], offset: float = 0.0) -> list[list[list[float]]]:
    """Per-region CA coordinates, translated `offset` angstroms along x.

    Residue i of region r sits at (i * 3.8 + offset, r * 10, 0). Two designs
    built from the same `region_lengths` therefore differ by a pure translation,
    and their design-region RMSD is exactly the difference of their offsets. That
    makes the clustering plot's distance matrix and dendrogram order arithmetic.
    """
    return [
        [[i * 3.8 + offset, r * 10.0, 0.0] for i in range(length)]
        for r, length in enumerate(region_lengths)
    ]


def rfdiff_design(
    index: int,
    *,
    fixed_runs: Sequence[int] = (5, 5),
    region_lengths: Sequence[int] = (4,),
    contact_resnums: Sequence[int] = (),
    effector_contact_residues: Sequence[int] = (),
    com_displacement: Sequence[float | None] | None = None,
    passes_filter: bool = True,
    frac_contacts_in_design: float | None = None,
    n_contacts_in_design: int | None = None,
    n_contact_pairs: int | None = None,
    coord_offset: float = 0.0,
    include_per_region_lengths: bool = True,
    include_coords: bool = True,
) -> dict:
    """One entry of `designs` in rfdiffusion_metrics.json.

    `contact_resnums` are receptor resnums as `receptor_layout` numbered them, so
    a contact inside the design region is one of the de novo numbers and a
    contact outside it is one of the fixed numbers. The specificity figures
    default to the counts implied by that split, so a call that states its
    contacts does not also have to restate its fractions.

    `include_per_region_lengths=False` and `include_coords=False` build the
    legacy-metrics shapes the module still has fallbacks for.
    """
    order, design_res, _fixed = receptor_layout(
        fixed_runs=fixed_runs, region_lengths=region_lengths)
    design_set = set(design_res)
    contacts = list(contact_resnums)
    n_in = sum(1 for r in contacts if r in design_set)
    n_pairs = n_contact_pairs if n_contact_pairs is not None else len(contacts)
    if n_contacts_in_design is None:
        n_contacts_in_design = n_in
    if frac_contacts_in_design is None:
        frac_contacts_in_design = (n_contacts_in_design / n_pairs) if n_pairs else 0.0

    entry: dict = {
        "design": f"design_{index}.pdb",
        "receptor_position_order": order,
        "per_design_design_residues": design_res,
        "receptor_contact_residues": contacts,
        "effector_contact_residues": list(effector_contact_residues),
        "n_contact_pairs": n_pairs,
        "n_contacts_in_design": n_contacts_in_design,
        "n_contacts_outside_design": max(n_pairs - n_contacts_in_design, 0),
        "frac_contacts_in_design": frac_contacts_in_design,
        "passes_filter": passes_filter,
    }
    if include_per_region_lengths:
        entry["per_region_lengths"] = list(region_lengths)
    if include_coords:
        prc = region_coords(region_lengths, coord_offset)
        entry["per_region_coords"] = prc
        entry["design_region_coords"] = [c for region in prc for c in region]
    if com_displacement is not None:
        entry["design_region_com_displacement"] = list(com_displacement)
    return entry


def rfdiffusion_metrics_json(
    path: Path,
    designs: Sequence[dict],
    *,
    fixed_residues: Sequence[int] = (),
    design_residues: Sequence[int] = (),
    min_hotspot_frac: float = 0.0,
) -> Path:
    """The metrics JSON `rfdiffusion_plots.py --metrics` reads.

    `min_hotspot_frac` is the specificity threshold the plot draws as a dashed
    line, and a value of 0 suppresses the line, which is what the committed
    reference set carries.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "designs": list(designs),
        "fixed_residues": list(fixed_residues),
        "design_residues": list(design_residues),
        "hotspot_residues": [],
        "min_hotspot_frac": min_hotspot_frac,
    }, indent=2))
    return path


# ══════════════════════════════════════════════════════════════════════
# bin/rosetta_filter_plots.py
# ══════════════════════════════════════════════════════════════════════

# The threshold rosetta_filter_plots.py falls back to when the metrics JSON
# carries no sc_threshold key. The committed reference set carries 0.5, so the
# default is only reached by a metrics file written before the key existed.
ROSETTA_DEFAULT_SC_THRESHOLD = 0.62


def rosetta_design(index: int, *, sc_value: float | None = 0.60,
                   passes_filter: bool | None = None, **physics) -> dict:
    """One entry of `designs` in rosetta_filter_metrics.json.

    The plot reads `sc_value` and nothing else, which is the module's own
    documented choice: at this stage the inputs are RFDiffusion polyvaline
    backbones, so `dG_separated`, `packstat` and `delta_unsatHbonds` are
    valine-clash artefacts. They are still written here, because the committed
    metrics file carries them and a fixture that dropped them would be a
    statement about the plot rather than about the file.

    `sc_value=None` is the shape Rosetta's InterfaceAnalyzer leaves behind when
    it fails on a design, and the plot's `d.get("sc_value") is not None` guard
    exists for exactly that. `passes_filter` defaults to nothing being claimed
    beyond what sc_value implies against the committed 0.5 threshold.
    """
    entry = {
        "design": f"design_{index}.pdb",
        "design_stem": f"design_{index}",
        "sc_value": sc_value,
        "dG_separated": 58.046,
        "dSASA_int": 1313.2,
        "dG_dSASA_density": 4.42,
        "packstat": 0.0,
        "delta_unsatHbonds": 1.0,
        "nres_int": 57,
        "passes_filter": (passes_filter if passes_filter is not None
                          else (sc_value is not None and sc_value >= 0.5)),
    }
    entry.update(physics)
    return entry


def rosetta_filter_metrics_json(path: Path, designs: Sequence[dict], *,
                                sc_threshold: float | None = 0.5) -> Path:
    """The metrics JSON `rosetta_filter_plots.py --metrics` reads.

    `sc_threshold=None` omits the key, which is the input that exercises the
    module's 0.62 fallback.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict = {"designs": list(designs)}
    if sc_threshold is not None:
        payload["sc_threshold"] = sc_threshold
    path.write_text(json.dumps(payload, indent=2))
    return path


# ══════════════════════════════════════════════════════════════════════
# bin/mpnn_plots.py
# ══════════════════════════════════════════════════════════════════════

MPNN_METADATA_COLUMNS = [
    "design", "seq", "mpnn_score", "native_residues", "designed_residues",
    "num_changes", "design_region_length_observed", "design_region_spec",
    "corrected_receptor", "receptor_total_length", "design_region_score",
]


def mpnn_row(
    *,
    design: str,
    seq: str = "0",
    mpnn_score: str | float = 1.0,
    design_region_score: str | float = 1.0,
    designed_residues: str = "AAAA",
    native_residues: str = "AAAA",
) -> dict:
    """One row of scored_metadata.csv.

    `designed_residues` and `native_residues` are the pipe-joined per-region
    strings `pipeline_correct_sequences.py` writes, e.g. "AVSVTEKGGWSDIT|GDIVVP".
    The AA-composition plot counts the letters and ignores the separator; the
    physicochem plot strips the separator before reducing each string to one
    number. Both scores are written as strings because the module reads the CSV
    with csv.DictReader and parses them itself, so "nan" and "" are reachable
    inputs.
    """
    designed_flat = designed_residues.replace("|", "")
    native_flat = native_residues.replace("|", "")
    return {
        "design": design,
        "seq": seq,
        "mpnn_score": str(mpnn_score),
        "native_residues": native_residues,
        "designed_residues": designed_residues,
        "num_changes": str(sum(1 for a, b in zip(designed_flat, native_flat) if a != b)),
        "design_region_length_observed": str(len(designed_flat)),
        "design_region_spec": f"{len(designed_flat)}-{len(designed_flat)}",
        "corrected_receptor": designed_flat,
        "receptor_total_length": str(len(designed_flat)),
        "design_region_score": str(design_region_score),
    }


def mpnn_metadata_csv(path: Path, rows: Sequence[dict]) -> Path:
    return _write_csv(path, MPNN_METADATA_COLUMNS, rows)


def mpnn_cluster_csv(path: Path, points: Sequence[tuple[float, int]],
                     n_sequences: int = 2) -> Path:
    """mpnn_cluster_counts.csv: one (identity threshold, cluster count) per row.

    The threshold column is a fraction in the file and the plot multiplies it by
    100 for the x-axis, so `points=[(0.3, 1), (1.0, 4)]` means an x-range of
    30-100%.
    """
    return _write_csv(
        path, ["threshold", "n_clusters", "n_sequences"],
        [{"threshold": f"{t:.2f}", "n_clusters": str(n), "n_sequences": str(n_sequences)}
         for t, n in points],
    )


# ══════════════════════════════════════════════════════════════════════
# bin/pose_solver_plots.py
# ══════════════════════════════════════════════════════════════════════

# The weight constants pose_solver.py records alongside the loss breakdown.
# Values are its module-level defaults; the plot annotates each bar with them.
POSE_LOSS_WEIGHTS = {
    "W_VALIDITY": 1e6,
    "W_DIST": 1e4,
    "W_LOWER": 1e4,
    "W_EXCL": 1e4,
    "W_INTERP": 10.0,
    "W_PAIR_SC_CLASH": 1e4,
}


def pose_pair(*, binder_res: int, target_res: int, achieved_d: float,
              min_d: float = 4.0, max_d: float = 6.0,
              satisfied: bool | None = None) -> dict:
    """One entry of `pair_results`.

    Omitting `satisfied` leaves the key out of the dict, which is the case the
    plot falls back to computing `min_d <= achieved_d <= max_d` itself.
    """
    pair = {
        "binder_res": binder_res, "target_res": target_res,
        "achieved_d": achieved_d, "min_d": min_d, "max_d": max_d,
    }
    if satisfied is not None:
        pair["satisfied"] = satisfied
    return pair


def pose_contig(*, contig: str, design_region: str, n_residues: int,
                join_gap: int) -> dict:
    """One value of the `contigs` map. `design_region` is the comma-separated
    resnum-range string pose_solver.py writes, e.g. "71-78" or "33-49, 69-78",
    with "(none)" for an empty region."""
    return {"contig": contig, "design_region": design_region,
            "n_residues": n_residues, "join_gap": join_gap}


def pose_results(
    *,
    loss_breakdown: dict | None = None,
    loss_weights: dict | None = None,
    pair_results: Sequence[dict] = (),
    contigs: dict | None = None,
    final_loss: float = 0.0,
    all_pairs_satisfied: bool = True,
) -> dict:
    """solved_pose_results.json as pose_solver.py writes it.

    Note what the committed fixtures show: the three results JSONs under
    tests/pose_solver/ carry no `loss_breakdown` or `loss_weights` at all, so a
    results file without them is a real input and the loss plot's placeholder
    path is a real path, not a contrivance.
    """
    results: dict = {
        "pair_results": list(pair_results),
        "all_pairs_satisfied": all_pairs_satisfied,
        "final_loss": final_loss,
        "contigs": dict(contigs) if contigs else {},
    }
    if loss_breakdown is not None:
        results["loss_breakdown"] = dict(loss_breakdown)
    if loss_weights is not None:
        results["loss_weights"] = dict(loss_weights)
    return results


def pose_results_json(path: Path, **kwargs) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(pose_results(**kwargs), indent=2))
    return path


def pose_interface_metrics(
    *,
    bsa_chain_a: float = 1100.0,
    bsa_chain_b: float = 1600.0,
    bsa_total: float | None = None,
    distance_histogram: dict | None = None,
    gap_index: float | None = 3.0,
    composition_a: dict | None = None,
    composition_b: dict | None = None,
) -> dict:
    """One entry of the list bin/pose_interface_metrics.py emits.

    Histogram keys are "<lo>-<hi>" strings; the dashboard colours each bar by
    the bin's lower edge, so the key, not the value, decides the colour.
    `gap_index=None` leaves gap_index_mean_A null, which is the "no gap-index"
    path.
    """
    if bsa_total is None:
        bsa_total = bsa_chain_a + bsa_chain_b
    if distance_histogram is None:
        distance_histogram = {"3.0-3.5": 10, "3.5-4.0": 20, "4.0-4.5": 30}
    if composition_a is None:
        composition_a = {"hydrophobic": 0.35, "polar": 0.26,
                         "charged": 0.39, "aromatic": 0.04}
    if composition_b is None:
        composition_b = {"hydrophobic": 0.45, "polar": 0.19,
                         "charged": 0.35, "aromatic": 0.19}
    return {
        "bsa_A2": {"chain_a": bsa_chain_a, "chain_b": bsa_chain_b,
                   "total": bsa_total},
        "heavy_contacts_5A": {"distance_histogram": dict(distance_histogram)},
        "shape_complementarity_proxy": {"gap_index_mean_A": gap_index},
        "interface_residues": {"composition_a": dict(composition_a),
                               "composition_b": dict(composition_b)},
    }


def pose_metrics_json(path: Path, entries: Sequence[dict] | dict) -> Path:
    """Write the interface-metrics JSON. A bare dict is written as a dict, which
    main() normalises to a one-entry list, so both shapes are exercisable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = entries if isinstance(entries, dict) else list(entries)
    path.write_text(json.dumps(payload, indent=2))
    return path


def pose_restart_losses_csv(path: Path,
                            rows: Sequence[tuple[int, float, float]]) -> Path:
    """solved_pose_restart_losses.csv: (restart_idx, restart_loss,
    best_loss_so_far) per row, header-only when `rows` is empty. The DE path
    writes exactly one row and the L-BFGS-B path writes one per restart."""
    return _write_csv(
        path, ["restart_idx", "restart_loss", "best_loss_so_far"],
        [{"restart_idx": str(i), "restart_loss": f"{loss:.6f}",
          "best_loss_so_far": f"{best:.6f}"} for i, loss, best in rows],
    )


# ══════════════════════════════════════════════════════════════════════
# bin/orthogonal_metrics_plots.py
# ══════════════════════════════════════════════════════════════════════

# survivors_with_orthogonal_metrics.csv is cross_sequence_summary_with_
# interface_metrics.csv plus the columns merge_orthogonal_metrics.py appends
# (its ALL_NEW).  Only the columns the plots read are listed, so the fixture
# stays a statement about the plots rather than a copy of a 67-column header.
ORTHOG_SURVIVOR_COLUMNS = [
    "mpnn_sequence", "row_type", "cross_tier",
    "rep_ra_eff_vs_truth_median", "rep_interface_plddt_median",
    "af3_nomsa_best_ra_eff", "af3_nomsa_total_predictions",
    "bsa", "sc", "rosetta_ddg",
    "orthogonal_flags", "passes_orthogonal_filters",
]

# Every cross_summary column a plot reads.  rep_reversion_verdict and
# rep_reverted_total_mutations are read by _rep_n_mutations but are NOT emitted
# by cross_sequence_summary.py (its rep_ columns come from extract_passing.py,
# which writes total_mutations_median already switched to the post-reversion
# count).  They are kept here, defaulting to empty, so the documented branch is
# reachable in a test while the default fixture matches the real header.
ORTHOG_CROSS_COLUMNS = [
    "mpnn_sequence", "row_type", "cross_tier", "cross_composite_score",
    "rep_total_mutations_median",
    "rep_ra_eff_vs_truth_median", "rep_independent_receptor_rmsd_median",
    "rep_true_jaccard_median", "rep_complex_plddt_median", "rep_iptm_median",
    "rep_ipae_median", "rep_pae_pass_frac_median",
    "rep_interface_plddt_median", "rep_ipsae_min_15_median",
    "rep_weighted_jaccard_median",
    "rep_reversion_verdict", "rep_reverted_total_mutations",
]

# Column values that pass every thresholded Boltz-2 confidence metric in
# COMBINED_SUMMARY_COLUMNS.  The combined-summary row sort promotes rows that
# pass all of them above rows that fail any, so a test about tier ordering needs
# a baseline that does not accidentally trip the tie-breaker.
BOLTZ_CONF_ALL_PASS = {
    "rep_complex_plddt_median": "0.90",
    "rep_iptm_median": "0.80",
    "rep_ipae_median": "5.00",
    "rep_pae_pass_frac_median": "0.90",
    "rep_interface_plddt_median": "0.90",
}


def orthogonal_survivor_row(
    *,
    mpnn_sequence: str,
    cross_tier: str = "none",
    boltz_ra_eff: str | float = "",
    af3_best_ra_eff: str | float = "",
    sc: str | float = "",
    bsa: str | float = "",
    rosetta_ddg: str | float = "",
    interface_plddt: str | float = "",
    af3_total_predictions: str | int = "",
    flags: Sequence[str] = (),
    passes: str | None = None,
    row_type: str = "steered",
) -> dict:
    """One row of survivors_with_orthogonal_metrics.csv.

    `flags` are the orthogonal_flags tokens merge_orthogonal_metrics.py writes:
    sc_missing, sc_too_low:<v>, bsa_missing, bsa_too_low:<v>,
    af3_nomsa_missing, af3_nomsa_ra_eff_too_high:<v>. `passes` defaults to the
    value the merge step would derive, which is 1 exactly when no non-AF3 flag
    is present, so a fixture cannot state a pass that contradicts its flags.
    """
    flag_list = list(flags)
    if passes is None:
        gating = [f for f in flag_list if not f.startswith("af3_nomsa")]
        passes = "0" if gating else "1"
    return {
        "mpnn_sequence": mpnn_sequence,
        "row_type": row_type,
        "cross_tier": cross_tier,
        "rep_ra_eff_vs_truth_median": str(boltz_ra_eff),
        "rep_interface_plddt_median": str(interface_plddt),
        "af3_nomsa_best_ra_eff": str(af3_best_ra_eff),
        "af3_nomsa_total_predictions": str(af3_total_predictions),
        "bsa": str(bsa),
        "sc": str(sc),
        "rosetta_ddg": str(rosetta_ddg),
        "orthogonal_flags": ",".join(flag_list),
        "passes_orthogonal_filters": passes,
    }


def orthogonal_survivors_csv(path: Path, rows: Sequence[dict]) -> Path:
    return _write_csv(path, ORTHOG_SURVIVOR_COLUMNS, rows)


def orthogonal_cross_row(
    *,
    mpnn_sequence: str,
    row_type: str = "steered",
    cross_tier: str = "none",
    composite: str | float = "",
    boltz_conf: dict | None = None,
    total_mutations: str | int = "",
    reversion_verdict: str = "",
    reverted_total_mutations: str | int = "",
    **extra: str | float,
) -> dict:
    """One row of cross_sequence_summary.csv.

    `boltz_conf` overrides the Boltz-2 confidence columns; the default is
    BOLTZ_CONF_ALL_PASS, so the combined-summary tie-breaker stays out of the way
    unless a test sets out to trip it. `extra` sets any other rep_ column by
    name, e.g. rep_true_jaccard_median.
    """
    row = {col: "" for col in ORTHOG_CROSS_COLUMNS}
    row.update({
        "mpnn_sequence": mpnn_sequence,
        "row_type": row_type,
        "cross_tier": cross_tier,
        "cross_composite_score": str(composite),
        "rep_total_mutations_median": str(total_mutations),
        "rep_reversion_verdict": reversion_verdict,
        "rep_reverted_total_mutations": str(reverted_total_mutations),
    })
    row.update({k: str(v) for k, v in (boltz_conf or BOLTZ_CONF_ALL_PASS).items()})
    row.update({k: str(v) for k, v in extra.items()})
    return row


def orthogonal_cross_summary_csv(path: Path, rows: Sequence[dict]) -> Path:
    return _write_csv(path, ORTHOG_CROSS_COLUMNS, rows)
