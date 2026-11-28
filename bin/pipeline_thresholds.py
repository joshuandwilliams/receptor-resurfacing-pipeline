"""
pipeline_thresholds.py
----------------------
PipelineInternalThresholds — the catalogue of every hard-coded research
threshold used across `bin/`.  Defaults live here, in one Python
module.  Callers can override them programmatically but normally don't.
These are research-methodology defaults chosen for this pipeline.

This module is where the pending threshold audit lands.  See
the design audit in https://github.com/joshuandwilliams/resurfacing-experiments.  Changing a threshold here is a research-
methodology change.  Changing a `PipelineParams` field is a user-
contract change with a different audit trail (see
bin/pipeline_params.py).

Design:
- Immutable.  Every override returns a new instance.
- Defaults are Python literals on the dataclass, version-controlled
  with the code.  Chosen over YAML so the methodology stays bundled
  with the code that uses it.
- Consumed by passing an instance through every method that does a
  threshold check, e.g. `StageResult.triggers_next_stage(thresholds)`
  and `OrthogonalMetrics.passes_orthogonal_filters(thresholds)`.

Imports nothing from the other typed modules in `bin/`, so any of them
can import it without a cycle.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
from typing import Any, Dict


@dataclass(frozen=True)
class PipelineInternalThresholds:
    """The catalogue of internal research-methodology thresholds.

    Every default is the value currently hard-coded somewhere in the
    `bin/` tree, or used as a default argument in a current function.
    The initial population mirrors what the codebase does today.  The
    pending threshold audit will revise individual values.

    All fields are immutable.  To produce a modified instance, use
    `with_overrides({...})`.  Never assign in place.
    """

    # ── Structural thresholds (negsteer gating) ───────────────────────
    #
    # `steered_ra_eff` and `reverted_ra_eff` set the "correctly placed"
    # bar for steered and reverted predictions respectively.  Today both
    # are 5.0 Å.  The comparison is strict: a prediction is correctly
    # placed when ra_eff < cutoff.  See
    # stage_result.StageResult.n_correctly_placed.
    # `cold_start_rmsd_threshold` mirrors
    # `params.negsteer_rmsd_threshold` (default 6.0 — main.nf:165),
    # which the negative-steering engine takes as --rmsd-threshold.  No
    # module in `bin/` reads this field yet.
    # `intact_threshold` is the receptor-fold-quality bar used by the
    # intact check (independent_receptor_rmsd < this).
    steered_ra_eff: float = 5.0
    reverted_ra_eff: float = 5.0
    cold_start_rmsd_threshold: float = 6.0
    intact_threshold: float = 5.0

    # ── Confidence thresholds (consumed by extract_passing) ───────────
    #
    # `extract_passing._compute_confidence_flag` lazily imports these
    # four defaults from here, and repeats the same literals as a
    # fallback so it stays importable without this module.
    # pae_pass_frac_min is the floor on the fraction of inter-chain
    # residue pairs whose PAE is below the cutoff.
    iptm_min: float = 0.30
    complex_plddt_min: float = 0.70
    ipae_max: float = 15.0
    pae_pass_frac_min: float = 0.10

    # ── Orthogonal gate (Sc + BSA + ΔΔG only) ─────────────────────────
    #
    # interface_plddt is Boltz-derived, so it is informational only and
    # NOT part of the gate.  bin/orthogonal_metrics_plots.py reads these
    # defaults from here and repeats the same literals as a standalone
    # fallback.  The ΔΔG default is the Bennett et al. 2023 calibration
    # reference.  af3_ra_max is informational only (NOT gated) but kept
    # here for the cohort-summary plot.
    orthogonal_sc_min: float = 0.55
    orthogonal_bsa_min: float = 600.0
    orthogonal_ddg_max: float = -30.0
    orthogonal_af3_ra_max: float = 5.0
    # Read only by the cohort-summary plot.  Explicitly NOT part of the
    # orthogonal_filters gate (see
    # the glossary in resurfacing-experiments,
    # passes_orthogonal_filters).
    orthogonal_plddt_min_informational: float = 0.75

    # ── Contact and interface geometry ────────────────────────────────
    #
    # contact_cutoff is the heavy-atom distance threshold for declaring
    # two residues "in contact".  Used by
    # structure_metrics.find_contact_residues_heavy and the negsteer
    # contamination check.  The weighted_jaccard_* constants match
    # compute_interface_metrics.py, where _WJ_MU=4.0 and
    # _WJ_TWO_SIGMA_SQ=2.25, so sigma = sqrt(2.25/2) ≈ 1.06.
    # weighted_jaccard_pair_cutoff mirrors
    # params.weighted_jaccard_pair_cutoff (default 8.0 — main.nf:253).
    contact_cutoff: float = 4.5
    weighted_jaccard_mu: float = 4.0
    weighted_jaccard_two_sigma_sq: float = 2.25
    weighted_jaccard_pair_cutoff: float = 8.0

    # ── Composite score ───────────────────────────────────────────────
    #
    # composite = true_jaccard − composite_ra_eff_weight × ra_eff.
    # Default 0.05, matching negsteer_plots.COMPOSITE_RA_EFF_WEIGHT.
    composite_ra_eff_weight: float = 0.05

    # ── Interface and control diagnostics ─────────────────────────────
    #
    # interface_plddt_trim_threshold: per-residue pLDDT below which an
    # interface residue's plddt is excluded from the interface_plddt
    # average.  Default 50.0 from main.nf:242.
    # controls_warning_*: thresholds that flag a negative control as
    # SUSPICIOUSLY GOOD, a sign the control isn't behaving as a control.
    # Defaults from main.nf:234-236.
    interface_plddt_trim_threshold: float = 50.0
    controls_warning_ipsae_max: float = 0.5
    controls_warning_ra_eff_min: float = 5.0

    # ── Sequence QC ───────────────────────────────────────────────────
    #
    # MPNN sequence quality control limits.  Defaults from the main.nf
    # parameters max_poly_x / min_pct_identity / max_pct_identity
    # (main.nf:275-277).
    max_poly_x: int = 5
    min_pct_identity: float = 0.0
    max_pct_identity: float = 100.0

    # ── Sc filter (Rosetta pre-MPNN) ──────────────────────────────────
    #
    # Shape complementarity threshold for the Rosetta pre-MPNN filter
    # pass.  Default 0.62 from main.nf:99 (Overath et al., 2025).
    sc_threshold: float = 0.62

    # ── Reversion gating ──────────────────────────────────────────────
    #
    # Reversion is gated per design, not per seed.  It fires when
    # n_correctly_placed > 0 and n_contaminated is at least
    # ceil(n_correctly_placed / 2).  The old per-seed rule let a single
    # contaminated seed trigger reversion for the whole design.  The
    # rule lives in StageResult.triggers_next_stage and reuses the
    # thresholds above, specifically steered_ra_eff for "correctly
    # placed" and contact_cutoff for contamination detection.  It needs
    # no threshold field of its own.  ceil(n / 2) is arithmetic, not a
    # configurable value.

    def as_dict(self) -> Dict[str, Any]:
        """Flat snapshot of all threshold values."""
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def with_overrides(
        self, overrides: Dict[str, Any]
    ) -> "PipelineInternalThresholds":
        """Return a NEW instance with the given overrides applied.

        The original instance is unchanged (frozen dataclass).  Raises
        ``TypeError`` if any override key is not a known threshold.
        That catches typos.
        """
        valid = {f.name for f in fields(self)}
        unknown = set(overrides) - valid
        if unknown:
            raise TypeError(
                "Unknown threshold override(s): "
                f"{sorted(unknown)}.  Valid: {sorted(valid)}."
            )
        return replace(self, **overrides)

    @classmethod
    def default(cls) -> "PipelineInternalThresholds":
        """Return an instance with all defaults, the research-
        methodology values currently in use.  Equivalent to ``cls()``
        but more discoverable from a caller."""
        return cls()
