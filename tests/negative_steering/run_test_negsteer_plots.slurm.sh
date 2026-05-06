#!/bin/bash
#SBATCH --job-name="negsteer_plots_test"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=test_plots_%j.out
#SBATCH --error=test_plots_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/negative_steering

# =============================================================================
# run_test_negsteer_plots_slurm.sh
# -----------------------------------------------------------------------------
# Re-renders the negative-steering plots from a previous test_negative_steering
# run's cross_sequence_summary.csv, plus per-sequence aggregated_results.csv
# files under runs/, plus optionally the upstream rfdiffusion metrics for
# design-region context.  Designed for fast iteration of the plot code itself
# against cached test outputs.
#
# Architecture note
# -----------------
# cross_sequence_summary.csv contains only sequences that produced at least
# one passing row.  Sequences with zero passing rows (tier-none cohort) are
# absent from the cross-summary but DO have aggregated_results.csv in their
# runs/<name>/ directory.  This wrapper passes both --csv (for cross-summary)
# and --runs-dir (for the full cohort) so the plots show every sequence in
# the cohort, not just the passing subset.
#
# Outputs (under <outdir>):
#   negsteer_tier_landscape.png       primary-ranker landscape
#                                     (composite per sequence, all tiers)
#   negsteer_seed_verdicts.png        stacked seed verdicts per sequence
#                                     (pose_holds / clean_steered / collapses
#                                     / contamination / no_data)
#   negsteer_ra_eff_vs_jaccard.png    scatter, "better" arrow toward upper-
#                                     left, single composite=0 reference,
#                                     two-panel (steered + off-scale strip)
#   negsteer_filter_cascade.png       confidence-filter waterfall — labels
#                                     read from --*-min flags so the plot
#                                     matches the values pinned in this run
#   negsteer_controls_diagnostic.png  steered-vs-control boxplots; ra_eff
#                                     panel shaded pale-red above threshold
#   negsteer_mutation_impact.png      per-input-PDB-position median composite
#                                     across the full PDB length, with the
#                                     design region shaded for context
#
# Usage:
#   sbatch tests/negative_steering/run_test_negsteer_plots_slurm.sh
#
# Optional CLI flags forwarded to the python script (override defaults):
#   --csv               PATH   override cross_sequence_summary.csv path
#   --runs-dir          PATH   override the per-sequence runs directory
#   --rfdiff-metrics    PATH   override rfdiffusion_metrics.json path
#   --complex-plddt-min FLOAT  override complex_plddt threshold
#   --ipae-max          FLOAT  override ipae threshold
#   --pae-pass-frac-min FLOAT  override pae_pass_frac threshold
#   --iptm-min          FLOAT  override iptm threshold
#   --outdir            PATH   override output dir
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
TEST_DIR="${PIPELINE_DIR}/tests/negative_steering"

# nextflow.config overrides params.outdir, so the real results dir is
# ${projectDir}/${params.project_name}_results = receptor_resurfacing_results.
RESULTS_DIR="${TEST_DIR}/receptor_resurfacing_results"
CROSS_CSV="${RESULTS_DIR}/negative_steering/cross_sequence_summary.csv"
RUNS_DIR="${RESULTS_DIR}/negative_steering/runs"

# Input-PDB design-region file produced by DERIVE_INPUT_INDICES — used
# by the mutation-impact plot to shade design-region positions.  These
# are 1-based input-PDB residue indices, NOT the 1000+ contig-internal
# numbers that show up in rfdiffusion_metrics.json.  Optional.
INPUT_DESIGN_REGION="${RESULTS_DIR}/negative_steering/controls_inputs/input_design_region.txt"

OUTDIR="${RESULTS_DIR}/plots_iter"

TEST_SCRIPT="${TEST_DIR}/test_negsteer_plots.py"

# ── Threshold defaults pinned to extract_passing.py values ───────────────
# Override at the SBATCH command line if you want to visualise a different
# cascade.  These exact values are baked into _compute_confidence_flag in
# extract_passing.py so the cascade matches the pipeline's own filter.
COMPLEX_PLDDT_MIN=0.70
IPAE_MAX=15.0
PAE_PASS_FRAC_MIN=0.10
IPTM_MIN=0.30

# ── Container ─────────────────────────────────────────────────────────────
# Same image as production NEGSTEER_PLOTS process.  Plotting depends only
# on matplotlib + numpy + stdlib (csv / json) — no GPU, no boltz2 / gemmi
# / MDAnalysis required.
RFDIFF_CONTAINER="/hpc-home/jowillia/singularity/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [ ! -f "${TEST_SCRIPT}" ]; then
    echo "ERROR: test script not found: ${TEST_SCRIPT}" >&2
    exit 1
fi

if [ ! -f "${CROSS_CSV}" ]; then
    echo "ERROR: cross_sequence_summary.csv not found: ${CROSS_CSV}" >&2
    echo "       Run test_negative_steering first to generate it." >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# Build optional --runs-dir arg (cohort plots include tier-none seqs when
# present; degrades gracefully to passing-only when absent).
RUNS_DIR_ARG=()
if [ -d "${RUNS_DIR}" ]; then
    RUNS_DIR_ARG=(--runs-dir "${RUNS_DIR}")
else
    echo "NOTE: runs/ not found at ${RUNS_DIR}"
    echo "      Cohort plots will only show sequences with passing rows;"
    echo "      tier-none sequences will be invisible."
fi

# Build optional --input-design-region arg.
INPUT_DR_ARG=()
if [ -f "${INPUT_DESIGN_REGION}" ]; then
    INPUT_DR_ARG=(--input-design-region "${INPUT_DESIGN_REGION}")
else
    echo "NOTE: input_design_region.txt not found at ${INPUT_DESIGN_REGION}"
    echo "      Mutation-impact plot will render without design-region shading."
fi

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Negative-Steering Plots Iteration — variant test"
echo "============================================================"
echo "Pipeline dir:      ${PIPELINE_DIR}"
echo "Test dir:          ${TEST_DIR}"
echo "Cross CSV:         ${CROSS_CSV}"
echo "Runs dir:          ${RUNS_DIR}"
echo "Input design region: ${INPUT_DESIGN_REGION}"
echo "Output dir:        ${OUTDIR}"
echo "Thresholds:        plddt≥${COMPLEX_PLDDT_MIN} ipae≤${IPAE_MAX}"
echo "                   paepf≥${PAE_PASS_FRAC_MIN} iptm≥${IPTM_MIN}"
echo "Date:              $(date)"
echo "Node:              $(hostname)"
echo "============================================================"

singularity exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --env MPLCONFIGDIR=/tmp \
    "${RFDIFF_CONTAINER}" \
    python "${TEST_SCRIPT}" \
        --csv               "${CROSS_CSV}" \
        "${RUNS_DIR_ARG[@]}" \
        "${INPUT_DR_ARG[@]}" \
        --complex-plddt-min "${COMPLEX_PLDDT_MIN}" \
        --ipae-max          "${IPAE_MAX}" \
        --pae-pass-frac-min "${PAE_PASS_FRAC_MIN}" \
        --iptm-min          "${IPTM_MIN}" \
        --outdir            "${OUTDIR}" \
        "$@"

echo ""
echo "============================================================"
echo "Plots written to: ${OUTDIR}"
echo "Finished: $(date)"
echo "============================================================"
