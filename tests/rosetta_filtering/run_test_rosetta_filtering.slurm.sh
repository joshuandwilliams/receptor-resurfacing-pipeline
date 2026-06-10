#!/bin/bash
#SBATCH --job-name="nf_test_rosetta"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/rosetta_filtering

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
NEXTFLOW_IMG="${PIPELINE_DIR}/containers/NextFlow.img"

TEST_DIR="${PIPELINE_DIR}/tests/rosetta_filtering"
NXF_HOME="${PIPELINE_DIR}/nxf_home"
NXF_WORK="${TEST_DIR}/work"
NXF_TEMP="${TEST_DIR}/tmp"

# ── Java ──────────────────────────────────────────────────────────────────
export JAVA_HOME="${PIPELINE_DIR}/containers/jdk-17.0.2"
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
echo "Rosetta Filtering Module Test — Nextflow Launcher"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Test dir:      ${TEST_DIR}"
echo "NXF_WORK:      ${NXF_WORK}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

"${NEXTFLOW_BIN}" run "${TEST_DIR}/test_rosetta_filtering.nf" \
    -c "${PIPELINE_DIR}/nextflow.config" \
    -resume \
    "$@"

echo ""
echo "============================================================"
echo "Test finished: $(date)"
echo "============================================================"