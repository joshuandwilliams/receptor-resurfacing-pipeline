#!/bin/bash
#SBATCH --job-name="nf_receptor_resurfacing"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=nextflow_pipeline_%j.out
#SBATCH --error=nextflow_pipeline_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk

set -euo pipefail

# ── Usage check ───────────────────────────────────────────────────────────
if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch run_pipeline.slurm.sh <path/to/params.yml> [extra nextflow args]"
    echo ""
    echo "Params file can specify either:"
    echo "  pdb_file:          Pre-docked complex PDB (skips the pose solver)"
    echo "  receptor_input:    Receptor monomer PDB (+ effector_input → pose solver docks them)"
    echo "  effector_input:    Effector monomer PDB (Branch A also requires pose_solver_pairs)"
    echo ""
    echo "Both inputs must be PDB files."
    exit 1
fi

PARAMS_FILE="$(realpath "$1")"
shift

if [ ! -f "${PARAMS_FILE}" ]; then
    echo "ERROR: params file not found: ${PARAMS_FILE}"
    exit 1
fi

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
NEXTFLOW_IMG="${PIPELINE_DIR}/containers/NextFlow.img"
# Combined container providing RFDiffusion, ProteinMPNN and MMseqs2.
# Its filename/recipe still bundle HADDOCK3, but the haddock3 CLI is
# unused — Branch A docking is now the pose solver (modules/pose_solver.nf).
# Variable name kept as RFDIFF_CONTAINER for backwards compatibility with
# the .nf modules that reference it as params.rfdiff_container.
RFDIFF_CONTAINER="${PIPELINE_DIR}/containers/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img"
AF2_DATA_DIR="/nbi/Reference-Data/AlphaFold/db-v2.3.2"

EXPERIMENT_DIR="$(dirname "${PARAMS_FILE}")"

# ── Java ──────────────────────────────────────────────────────────────────
export JAVA_HOME="${PIPELINE_DIR}/containers/jdk-17.0.2"
export PATH="${JAVA_HOME}/bin:${PATH}"

# ── Nextflow environment ──────────────────────────────────────────────────
export NXF_OFFLINE=true
export NXF_PLUGINS_DEFAULT=false
export NXF_HOME="${PIPELINE_DIR}/nxf_home"
export NXF_WORK="${EXPERIMENT_DIR}/work"
export NXF_TEMP="${EXPERIMENT_DIR}/tmp"

mkdir -p "${NXF_HOME}" "${NXF_WORK}" "${NXF_TEMP}"

# ── Extract nextflow binary if needed ─────────────────────────────────────
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
echo "Receptor Resurfacing Pipeline v0.4.0 — Nextflow Launcher"
echo "============================================================"
echo "Params file:     ${PARAMS_FILE}"
echo "Experiment dir:  ${EXPERIMENT_DIR}"
echo "Pipeline dir:    ${PIPELINE_DIR}"
echo "NXF_HOME:        ${NXF_HOME}"
echo "NXF_WORK:        ${NXF_WORK}"
echo "NXF_TEMP:        ${NXF_TEMP}"
echo "JAVA_HOME:       ${JAVA_HOME}"
echo "Java:            $(java -version 2>&1 | head -1)"
echo "sbatch:          $(which sbatch)"
echo "Nextflow:        ${NEXTFLOW_BIN}"
echo "Date:            $(date)"
echo "Node:            $(hostname)"
echo "============================================================"

"${NEXTFLOW_BIN}" run "${PIPELINE_DIR}/main.nf" \
    -params-file "${PARAMS_FILE}" \
    --rfdiff_container "${RFDIFF_CONTAINER}" \
    --af2_data_dir "${AF2_DATA_DIR}" \
    -resume \
    "$@"

echo ""
echo "============================================================"
echo "Pipeline finished: $(date)"
echo "============================================================"
