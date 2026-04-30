#!/bin/bash
#SBATCH --job-name="nf_test_orthogonal"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk

# -----------------------------------------------------------------------------
# Launch wrapper for the orthogonal-metrics test (P0-31).
#
# Runs the full post-negative-steering cascade on a completed
# test_negative_steering output dir.  Expects:
#
#   tests/negative_steering/receptor_resurfacing_results/negative_steering/
#       cross_sequence_summary.csv
#       runs/design_*_seq_*/{plan.json, effector_template.cif, cycle_0/...}
#
# Five stages downstream of that:
#   1. NEGSTEER_INTERFACE_METRICS     — P0-29
#   2. EXTRACT_SURVIVOR_MANIFEST
#   3. AF3-no-MSA         \\
#      biophysical        |  — three streams, each fan-out per survivor
#      Rosetta            /
#   4. NEGSTEER_ORTHOGONAL_METRICS    — merge + filter
#
# The launcher itself is lightweight (jic-medium, 2 cpu, 4 GB).  Heavy
# work happens in the spawned process jobs — AF3_NOMSA_ON_SURVIVORS hits
# jic-gpu with up to max_af3_parallel=30 concurrent jobs.  Wall-clock
# budget of 24h is generous for a 20-survivor test; real wall will be
# dominated by AF3 (~5-10 min per survivor) plus queue waits.
# -----------------------------------------------------------------------------

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
NEXTFLOW_IMG="/hpc-home/jowillia/singularity/NextFlow/NextFlow.img"

TEST_DIR="${PIPELINE_DIR}/tests/orthogonal_metrics"
NXF_HOME="${PIPELINE_DIR}/nxf_home"
NXF_WORK="${TEST_DIR}/work"
NXF_TEMP="${TEST_DIR}/tmp"

# ── Java ──────────────────────────────────────────────────────────────────
export JAVA_HOME="${PIPELINE_DIR}/jdk-17.0.2"
export PATH="${JAVA_HOME}/bin:${PATH}"

# ── Nextflow environment ──────────────────────────────────────────────────
export NXF_OFFLINE=true
export NXF_PLUGINS_DEFAULT=false
export NXF_HOME="${NXF_HOME}"
export NXF_WORK="${NXF_WORK}"
export NXF_TEMP="${NXF_TEMP}"

mkdir -p "${NXF_HOME}" "${NXF_WORK}" "${NXF_TEMP}"

# ── Nextflow binary ───────────────────────────────────────────────────────
NEXTFLOW_BIN="${NXF_HOME}/nextflow"
if [ ! -x "${NEXTFLOW_BIN}" ]; then
    echo "Extracting nextflow binary from container..."
    singularity exec "${NEXTFLOW_IMG}" cat /usr/local/bin/nextflow > "${NEXTFLOW_BIN}"
    chmod +x "${NEXTFLOW_BIN}"
fi

# ── Seed plugin cache if needed ───────────────────────────────────────────
if [ ! -d "${NXF_HOME}/plugins" ] || [ -z "$(ls -A "${NXF_HOME}/plugins" 2>/dev/null)" ]; then
    echo "Seeding NXF_HOME from container..."
    singularity exec --bind "${NXF_HOME}:/mnt/out" "${NEXTFLOW_IMG}" \
        bash -c "cp -r /opt/nextflow/* /mnt/out" 2>/dev/null || true
fi

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Orthogonal-Metrics Test — Nextflow Launcher (P0-31)"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Test dir:      ${TEST_DIR}"
echo "NXF_WORK:      ${NXF_WORK}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

"${NEXTFLOW_BIN}" run "${TEST_DIR}/test_orthogonal_metrics.nf" \
    -c "${PIPELINE_DIR}/nextflow.config" \
    -resume \
    "$@"

echo ""
echo "============================================================"
echo "Test finished: $(date)"
echo "============================================================"
