#!/usr/bin/env python3
"""
validate_params.py
------------------
Fail-fast parameter validator for the Nextflow pipeline.  Invoked at
the very top of `workflow {}` in main.nf, BEFORE any process is
dispatched.  Any constraint violation prints the full list of errors
to stderr and exits non-zero, blocking submission of SLURM jobs that
would otherwise burn GPU time before failing.

Usage
=====
    python3 bin/validate_params.py --params-json <path>

The JSON file is expected to be a flat object of `params.*` key:value
pairs as serialised by Nextflow's `groovy.json.JsonOutput.toJson(params)`.
Keys not in PARAM_SPECS are reported as warnings (unknown params are
not fatal — they may be legitimate experimental flags).

Spec model
==========
Each entry in PARAM_SPECS is a `ParamSpec(name, kind, **kwargs)` where
`kind` is one of:

  - `int_range(min, max)`     — integer with inclusive bounds
  - `float_range(min, max)`   — float with inclusive bounds; either
                                 bound may be None for one-sided
  - `choice(values)`          — string must be in `values`
  - `bool`                    — JSON true/false
  - `non_empty_str`           — non-empty string
  - `optional_path`           — None OR a non-empty string (existence
                                 checked at Nextflow level, not here)
  - `regex(pattern)`          — string matches pattern
  - `any`                     — no constraint (escape hatch for params
                                 still needing a spec — flagged in
                                 coverage report)
  - `custom(fn)`              — callable accepting the value, returning
                                 None on success or an error string on
                                 failure

Coverage report
===============
Pass `--report` to print which params have specs and which don't.  Use
this to track coverage as new params are added.

Tests
=====
See tests/characterization/test_validate_params.py for unit coverage of
each validator kind.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, List, Optional, Sequence, Union


# ── Validator kinds ──────────────────────────────────────────────────


@dataclass
class ParamSpec:
    """Single-parameter validation spec.

    `kind` is the validator name; `kwargs` carries kind-specific
    arguments (bounds, allowed values, etc.).  `note` is a one-line
    rationale surfaced in the coverage report.
    """
    name: str
    kind: str
    kwargs: dict = field(default_factory=dict)
    note: str = ""


# Each validator returns None on success or an error string on failure.
# Validators take (value, kwargs) and may need to coerce.

def _validate_int_range(v: Any, kw: dict) -> Optional[str]:
    try:
        iv = int(v)
    except (TypeError, ValueError):
        return f"expected integer, got {v!r} ({type(v).__name__})"
    lo, hi = kw.get("min"), kw.get("max")
    if lo is not None and iv < lo:
        return f"value {iv} < min {lo}"
    if hi is not None and iv > hi:
        return f"value {iv} > max {hi}"
    return None


def _validate_float_range(v: Any, kw: dict) -> Optional[str]:
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return f"expected number, got {v!r} ({type(v).__name__})"
    if fv != fv:  # NaN
        return "value is NaN"
    lo, hi = kw.get("min"), kw.get("max")
    if lo is not None and fv < lo:
        return f"value {fv} < min {lo}"
    if hi is not None and fv > hi:
        return f"value {fv} > max {hi}"
    return None


def _validate_choice(v: Any, kw: dict) -> Optional[str]:
    values: Sequence = kw["values"]
    if v not in values:
        return f"value {v!r} not in allowed set {sorted(values)}"
    return None


def _validate_bool(v: Any, kw: dict) -> Optional[str]:
    if not isinstance(v, bool):
        return f"expected bool, got {v!r} ({type(v).__name__})"
    return None


def _validate_non_empty_str(v: Any, kw: dict) -> Optional[str]:
    if not isinstance(v, str):
        return f"expected string, got {type(v).__name__}"
    if not v.strip():
        return "empty string"
    return None


def _validate_optional_path(v: Any, kw: dict) -> Optional[str]:
    if v is None:
        return None
    if not isinstance(v, str):
        return f"expected None or string path, got {type(v).__name__}"
    if not v.strip():
        return "empty path string"
    return None


def _validate_regex(v: Any, kw: dict) -> Optional[str]:
    if not isinstance(v, str):
        return f"expected string, got {type(v).__name__}"
    pat = kw["pattern"]
    if not re.fullmatch(pat, v):
        return f"value {v!r} does not match pattern {pat!r}"
    return None


def _validate_any(v: Any, kw: dict) -> Optional[str]:
    return None


def _validate_custom(v: Any, kw: dict) -> Optional[str]:
    fn: Callable[[Any], Optional[str]] = kw["fn"]
    return fn(v)


VALIDATORS = {
    "int_range":      _validate_int_range,
    "float_range":    _validate_float_range,
    "choice":         _validate_choice,
    "bool":           _validate_bool,
    "non_empty_str":  _validate_non_empty_str,
    "optional_path":  _validate_optional_path,
    "regex":          _validate_regex,
    "any":            _validate_any,
    "custom":         _validate_custom,
}


def _check_num_seeds_odd(v: Any) -> Optional[str]:
    """negsteer_num_seeds must be odd.  An even seed count lets per-seed
    verdict majority votes tie, breaking the unambiguous-winner contract
    in _classify_outcome / _per_seed_verdict_breakdown.  See the comment
    in main.nf workflow body and the audit Q121 close."""
    err = _validate_int_range(v, {"min": 1, "max": 99})
    if err is not None:
        return err
    if int(v) % 2 == 0:
        return f"value {v} must be ODD (3 recommended); even values let " \
               f"per-seed verdict majority ties bypass downgrade-on-tie."
    return None


# ── Param specs ──────────────────────────────────────────────────────
#
# Coverage: most pipeline params have meaningful specs.  A handful
# remain marked `any` with a TODO note — these are coverage gaps to
# close as the threshold audit (Task 47 in design_audit.md) decides
# each value's valid range.  An unspecced param does NOT block the
# pipeline; the coverage report (--report) lists them so reviewers
# can prioritise.

PARAM_SPECS: List[ParamSpec] = [
    # ── Identification ──────────────────────────────────────────────
    ParamSpec("project_name",  "regex",
              {"pattern": r"[A-Za-z0-9_.\-]+"},
              "alphanumeric + ._- only; used as a directory name"),
    ParamSpec("outdir",        "non_empty_str", {},
              "output directory; existence not checked here (workflow creates it)"),

    # ── Branch A / B inputs ─────────────────────────────────────────
    ParamSpec("pdb_file",      "optional_path", {},
              "pre-docked complex; Branch B input"),
    ParamSpec("receptor_input","optional_path", {},
              "Branch A input — receptor PDB"),
    ParamSpec("effector_input","optional_path", {},
              "Branch A input — effector PDB"),
    ParamSpec("receptor_chain","regex", {"pattern": r"[A-Za-z]"},
              "single chain letter"),
    ParamSpec("effector_chain","regex", {"pattern": r"[A-Za-z]"},
              "single chain letter"),
    ParamSpec("rfdiff_output_receptor_chain","regex",
              {"pattern": r"[A-Za-z]"}, "single chain letter"),
    ParamSpec("rfdiff_output_effector_chain","regex",
              {"pattern": r"[A-Za-z]"}, "single chain letter"),

    # ── RFDiffusion ─────────────────────────────────────────────────
    ParamSpec("num_designs",   "int_range", {"min": 1, "max": 10000},
              "RFDiffusion candidate count"),
    ParamSpec("rfdiff_iterations","int_range", {"min": 1, "max": 200},
              "RFDiffusion denoising steps"),
    ParamSpec("rfdiff_checkpoint","non_empty_str", {},
              "checkpoint filename inside /opt/RFdiffusion/models/"),
    ParamSpec("min_hotspot_frac","float_range", {"min": 0.0, "max": 1.0},
              "min fraction of contacts inside design region"),
    ParamSpec("symmetry",      "choice", {"values": ["none", "C2", "C3", "C4", "C5", "C6"]},
              "RFDiffusion symmetry mode"),
    ParamSpec("order",         "int_range", {"min": 1, "max": 12},
              "symmetry order"),
    ParamSpec("add_potential", "bool", {},
              "enable guiding potentials during diffusion"),
    ParamSpec("rfdiff_guide_scale","float_range", {"min": 0.0, "max": 100.0},
              "global potential multiplier"),
    ParamSpec("rfdiff_guide_decay","choice",
              {"values": ["constant", "linear", "quadratic", "cubic"]},
              "weight decay schedule"),
    ParamSpec("rfdiff_interface_weight","float_range", {"min": 0.0, "max": 100.0}),
    ParamSpec("rfdiff_rog_weight","float_range", {"min": 0.0, "max": 100.0}),
    ParamSpec("rfdiff_rog_min_dist","float_range", {"min": 0.0, "max": 50.0}),

    # ── Rosetta pre-MPNN filter ─────────────────────────────────────
    ParamSpec("sc_threshold",  "float_range", {"min": 0.0, "max": 1.0},
              "shape complementarity threshold"),
    ParamSpec("stop_after_rosetta","bool", {},
              "halt cleanly after Sc filtering"),

    # ── ProteinMPNN ─────────────────────────────────────────────────
    ParamSpec("num_seqs",      "int_range", {"min": 1, "max": 1000},
              "MPNN sequences per RFDiffusion design"),
    ParamSpec("mpnn_sampling_temp","float_range", {"min": 0.0, "max": 5.0},
              "MPNN sampling temperature"),
    ParamSpec("rm_aa",         "regex", {"pattern": r"[ACDEFGHIKLMNPQRSTVWY]*"},
              "amino acids to exclude (one-letter codes)"),
    ParamSpec("mpnn_top_n",    "int_range", {"min": 0, "max": 100000},
              "select top N by MPNN score (0 = all)"),

    # ── Negative steering — strategy ────────────────────────────────
    ParamSpec("negsteer_mode", "choice",
              {"values": ["strong", "mild", "conservative", "alanine"]},
              "steering mutation strategy"),
    ParamSpec("negsteer_max_mutations","int_range", {"min": 1, "max": 50},
              "receptor residues mutated per cycle"),
    ParamSpec("negsteer_candidate_pool_size","int_range", {"min": 1, "max": 200},
              "candidate residue pool size"),

    # ── Negative steering — counts ──────────────────────────────────
    ParamSpec("negsteer_n_designs","int_range", {"min": 1, "max": 200},
              "steered designs per MPNN sequence"),
    ParamSpec("negsteer_num_seeds","custom", {"fn": _check_num_seeds_odd},
              "Boltz seeds per sequence — MUST be odd (audit Q121)"),
    ParamSpec("negsteer_n_cycles","int_range", {"min": 1, "max": 10},
              "steering cycles (single-cycle = 1)"),

    # ── Negative steering — Boltz hyperparameters ───────────────────
    ParamSpec("negsteer_diffusion_samples","int_range", {"min": 1, "max": 50},
              "diffusion samples per seed"),
    ParamSpec("negsteer_recycling_steps","int_range", {"min": 1, "max": 20},
              "Boltz recycling steps"),
    ParamSpec("negsteer_rmsd_threshold","float_range", {"min": 0.0, "max": 20.0},
              "initial-RMSD threshold for skip_steering"),
    ParamSpec("negsteer_contact_cutoff","float_range", {"min": 0.0, "max": 20.0},
              "Heavy-atom contact cutoff (Å)"),
    ParamSpec("negsteer_no_kernels","bool", {}),

    # ── Negative steering — post-processing ─────────────────────────
    ParamSpec("negsteer_postprocess_rmsd_threshold","float_range",
              {"min": 0.0, "max": 20.0}),
    ParamSpec("negsteer_postprocess_metric_column","choice",
              {"values": ["steered_ra_eff_vs_truth",
                          "reverted_ra_eff_vs_truth"]}),
    ParamSpec("negsteer_postprocess_contact_cutoff","float_range",
              {"min": 0.0, "max": 20.0}),

    # ── Negative controls ───────────────────────────────────────────
    ParamSpec("run_negative_controls","bool", {}),
    ParamSpec("negsteer_controls_n_designs","int_range", {"min": 1, "max": 50}),
    ParamSpec("controls_warning_ipsae_max","float_range", {"min": 0.0, "max": 1.0}),
    ParamSpec("controls_warning_ra_eff_min","float_range", {"min": 0.0, "max": 20.0}),

    # ── Orthogonal filter thresholds ────────────────────────────────
    ParamSpec("orthogonal_filter_sc_min","float_range", {"min": 0.0, "max": 1.0},
              "Lawrence-Colman Sc threshold"),
    ParamSpec("orthogonal_filter_bsa_min","float_range", {"min": 0.0, "max": 5000.0},
              "BSA threshold (Å²)"),
    ParamSpec("orthogonal_filter_plddt_min","float_range", {"min": 0.0, "max": 1.0},
              "interface pLDDT threshold (informational only after audit)"),
    ParamSpec("orthogonal_filter_af3_ra_max","float_range", {"min": 0.0, "max": 20.0},
              "AF3 ra_eff threshold"),

    # ── Interface metrics ───────────────────────────────────────────
    ParamSpec("interface_plddt_trim_threshold","float_range", {"min": 0.0, "max": 100.0}),
    ParamSpec("interface_intact_threshold","float_range", {"min": 0.0, "max": 20.0}),
    ParamSpec("weighted_jaccard_pair_cutoff","float_range", {"min": 0.0, "max": 20.0}),

    # ── Sequence QC ─────────────────────────────────────────────────
    ParamSpec("max_poly_x",    "int_range", {"min": 1, "max": 100}),
    ParamSpec("min_pct_identity","float_range", {"min": 0.0, "max": 100.0}),
    ParamSpec("max_pct_identity","float_range", {"min": 0.0, "max": 100.0}),

    # ── Pose Solver (Branch A docking) ──────────────────────────────
    ParamSpec("rfdiff_contact_cutoff","float_range", {"min": 0.0, "max": 30.0}),

    # Hard-pin CA-CA pair restraints.  Empty string allowed everywhere;
    # cross-param "non-empty in Branch A" is enforced in
    # _branch_a_pose_solver_pairs_present below.
    #
    # Format: space-separated 'CHAIN+RES-CHAIN+RES' tokens with optional
    # '@DIST' suffix overriding the per-pair max distance, e.g.
    #   "A73-B31 A71-B33 A8-B22"
    #   "A73-B48@4.5 A71-B50 A8-B39"
    ParamSpec("pose_solver_pairs", "regex",
              {"pattern": r"^(\s*[A-Za-z]\d+-[A-Za-z]\d+(@\d+(\.\d+)?)?(?:\s+|$))*\s*$"},
              "space-separated 'CHAIN+RES-CHAIN+RES[@DIST]' pair tokens"),
    ParamSpec("pose_solver_exclusions", "regex",
              {"pattern": r"^(\s*[A-Za-z]\d+-[A-Za-z]\d+@\d+(\.\d+)?(?:\s+|$))*\s*$"},
              "space-separated 'CHAIN+RES-CHAIN+RES@DIST' exclusion tokens "
              "(@DIST is mandatory for exclusions)"),
    ParamSpec("pose_solver_min_pair_distance","float_range", {"min": 0.0, "max": 30.0},
              "minimum CA-CA distance for pair constraints (Å)"),
    ParamSpec("pose_solver_max_pair_distance","float_range", {"min": 0.0, "max": 30.0},
              "maximum CA-CA distance for pair constraints (Å); "
              "default 6.0 matches the natural-interface CA-CA range "
              "documented in notes/pipeline_notes/pipeline_notes16.md §4.1"),
    ParamSpec("pose_solver_pair_sc_clash_cutoff","float_range", {"min": 0.0, "max": 10.0},
              "heavy-atom clash cutoff between paired-residue sidechains (Å); "
              "set 0 to disable"),
    ParamSpec("pose_solver_clash_cutoff","float_range", {"min": 0.0, "max": 10.0},
              "heavy-atom distance counted as a clash in the post-hoc report (Å)"),
    ParamSpec("pose_solver_contact_cutoff","float_range", {"min": 0.0, "max": 30.0},
              "CA-CA cutoff for the contact heatmap (Å)"),
    ParamSpec("pose_solver_n_restarts","int_range", {"min": 1, "max": 100000},
              "number of L-BFGS-B random restarts (or DE evaluations / 90)"),
    ParamSpec("pose_solver_use_de", "bool", {},
              "use differential evolution instead of L-BFGS-B random restarts"),
    ParamSpec("pose_solver_global_interp", "bool", {},
              "penalise ALL CA atoms inside the opposite hull (vs pair-only)"),
    ParamSpec("pose_solver_interp_weight","float_range", {"min": 0.0, "max": 1e6},
              "weight on the interpenetration penalty"),
    ParamSpec("pose_solver_contig_design_region", "regex",
              {"pattern": r"^(\d+(-\d+)?(,\s*\d+(-\d+)?)*)?$"},
              "cosmetic-only annotation for the contact heatmap: comma-"
              "separated residue numbers + ranges"),
    ParamSpec("stop_after_pose_solve", "bool", {},
              "halt after POSE_SOLVER_PLOTS for manual inspection"),

    # ── Orthogonal cascade ──────────────────────────────────────────
    # Which cross_tier values get the AF3 + biophysical + Rosetta
    # cascade.  Default 'all' includes tier-none failed designs
    # (diagnostic).  'abc' restricts to survivors (cheaper for GPU).
    ParamSpec("orthogonal_tier_filter","choice",
              {"values": ["all", "abc"]},
              "which cross_tier values get the orthogonal cascade"),

    # ── Infrastructure ──────────────────────────────────────────────
    # Container / model-store paths are non_empty_str — existence is
    # checked by the workflow when the container is invoked, not here.
    ParamSpec("max_boltz2_parallel","int_range", {"min": 1, "max": 100}),
    ParamSpec("max_af3_parallel","int_range", {"min": 1, "max": 100},
              "maxForks cap for AF3_NOMSA_ON_SURVIVORS"),
    ParamSpec("rfdiff_container",  "non_empty_str", {},
              "Singularity image path for RFDiffusion + MPNN + scipy "
              "(used by RFDIFFUSION, POSE_SOLVE, POSE_SOLVER_PREPARE, "
              "POSE_SOLVER_PLOTS, BUILD_CONTIGS)"),
    ParamSpec("rosetta_container", "non_empty_str", {},
              "Singularity image path for Rosetta"),
    ParamSpec("boltz2_container",  "non_empty_str", {},
              "Singularity image path for Boltz2 negsteer"),
    ParamSpec("colabfold_container","non_empty_str", {},
              "Singularity image path for ColabFold"),
    ParamSpec("af3_package_id",    "regex",
              {"pattern": r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}"
                          r"-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"},
              "NBI source-package UUID for AlphaFold 3"),
    ParamSpec("af3_model_dir",     "non_empty_str", {},
              "directory containing af3.bin"),
    ParamSpec("af3_db_v3",         "non_empty_str", {},
              "AlphaFold 3 reference database root"),
    ParamSpec("af2_data_dir",      "non_empty_str", {},
              "AlphaFold 2 BFD + MGnify database root"),

    # ── Coverage gaps (specs TBD per threshold audit) ──────────────
    # Listed explicitly so the coverage report surfaces them, not just
    # treated as "unknown".  Add a real spec as each is decided.
    ParamSpec("contigs",                "any", {}, "free-form RFDiffusion contig string"),
    ParamSpec("hotspot",                "any", {}, "RFDiffusion hotspot residue list (string)"),
    ParamSpec("receptor_seq",           "any", {}, "auto-derived if blank"),
    ParamSpec("effector_seq",           "any", {}, "auto-derived if blank"),
    ParamSpec("negsteer_protected_set_source", "any", {},
              "protected-set source mode — TODO spec"),
    ParamSpec("af3_nomsa_seeds",        "any", {}, "list of 3 ints"),
]


# ── Validation engine ────────────────────────────────────────────────


def _branch_a_pose_solver_pairs_present(params: dict) -> Optional[str]:
    """Cross-param rule: when Branch A is active (params.receptor_input
    set), the user MUST supply at least one pose_solver_pairs token.
    The pose solver is constraint-driven; blind docking is not
    supported.  Returns an error message if violated, else None.
    """
    if not params.get("receptor_input"):
        return None  # Not Branch A; pose solver is skipped.
    pairs = (params.get("pose_solver_pairs") or "").strip()
    if not pairs:
        return (
            "Branch A (params.receptor_input set) requires "
            "params.pose_solver_pairs to contain at least one "
            "CA-CA pair restraint.  Example:\n"
            "    pose_solver_pairs: \"A73-B31 A71-B33 A8-B22\"\n"
            "Use 2-4 pairs covering the intended interface.  See "
            "notes/pipeline_notes/pipeline_notes16.md for guidance on "
            "choosing pair residues and distance bounds."
        )
    return None


def _pose_solver_pair_distance_bounds_ordered(params: dict) -> Optional[str]:
    """params.pose_solver_min_pair_distance must be < ...max_pair_distance."""
    lo = params.get("pose_solver_min_pair_distance")
    hi = params.get("pose_solver_max_pair_distance")
    if lo is None or hi is None:
        return None
    try:
        if float(lo) >= float(hi):
            return (
                f"params.pose_solver_min_pair_distance ({lo}) must be strictly "
                f"less than params.pose_solver_max_pair_distance ({hi})."
            )
    except (TypeError, ValueError):
        return None  # type errors are caught by the per-param spec
    return None


def validate_params(params: dict) -> List[str]:
    """Return a list of error strings (empty list = all valid)."""
    errors: List[str] = []
    spec_names = {s.name for s in PARAM_SPECS}

    for spec in PARAM_SPECS:
        if spec.name not in params:
            # Missing keys are allowed — they'll fall through to
            # Nextflow defaults.  Only validate keys that are set.
            continue
        validator = VALIDATORS[spec.kind]
        err = validator(params[spec.name], spec.kwargs)
        if err is not None:
            errors.append(f"  - params.{spec.name}: {err}")

    # Cross-param rules — only meaningful once per-param shapes are valid.
    for rule in (_branch_a_pose_solver_pairs_present,
                 _pose_solver_pair_distance_bounds_ordered):
        msg = rule(params)
        if msg is not None:
            errors.append(f"  - {msg}")

    # Unknown params: warn but don't fail.
    unknown = sorted(k for k in params.keys() if k not in spec_names)
    if unknown:
        print(f"NOTE: {len(unknown)} param(s) have no spec entry (not validated): "
              f"{', '.join(unknown[:10])}{'...' if len(unknown) > 10 else ''}",
              file=sys.stderr)

    return errors


def coverage_report() -> str:
    lines = []
    by_kind: dict = {}
    for s in PARAM_SPECS:
        by_kind.setdefault(s.kind, []).append(s.name)
    lines.append(f"Total spec entries: {len(PARAM_SPECS)}")
    for kind in sorted(by_kind):
        lines.append(f"  {kind:<16s} : {len(by_kind[kind])}")
    gap = sorted(by_kind.get("any", []))
    if gap:
        lines.append("")
        lines.append(f"Coverage gaps ({len(gap)} param(s) marked 'any'):")
        for name in gap:
            lines.append(f"  - {name}")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Validate Nextflow params against PARAM_SPECS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--params-json", type=Path,
                   help="Path to a JSON file containing the flat params object.")
    p.add_argument("--report", action="store_true",
                   help="Print coverage report and exit.")
    args = p.parse_args(argv)

    if args.report:
        print(coverage_report())
        return 0

    if args.params_json is None:
        p.error("--params-json is required unless --report is set.")

    if not args.params_json.is_file():
        print(f"ERROR: params JSON not found: {args.params_json}", file=sys.stderr)
        return 2
    try:
        params = json.loads(args.params_json.read_text())
    except json.JSONDecodeError as e:
        print(f"ERROR: malformed params JSON: {e}", file=sys.stderr)
        return 2
    if not isinstance(params, dict):
        print(f"ERROR: params JSON must be an object, got {type(params).__name__}",
              file=sys.stderr)
        return 2

    errors = validate_params(params)
    if errors:
        print("Parameter validation FAILED:", file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)
        return 1
    print(f"Parameter validation OK ({len(PARAM_SPECS)} specs checked, "
          f"{len(params)} params provided).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
