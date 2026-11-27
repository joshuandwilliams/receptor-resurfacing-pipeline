#!/bin/bash
#SBATCH --job-name="nf_resample_mpnn"
# No -p: the partition is site-specific.  tests/run_tests.sh passes
# RECEPTOR_PARTITION_CPU from site.conf as `sbatch -p`, where it beats both
# the environment and a header line.  Submitting this script directly takes
# the cluster's default partition.
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=nextflow_pipeline_%j.out
#SBATCH --error=nextflow_pipeline_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user: sbatch(1) mails the submitting user when only --mail-type
# is given, so an address here would be site-specific and redundant.

set -euo pipefail

# ── Usage check ───────────────────────────────────────────────────────────
if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch run_resample_mpnn.slurm.sh <path/to/params.yml> [extra nextflow args]"
    echo ""
    echo "Runs resample_mpnn.nf — regenerates ProteinMPNN sequences for a"
    echo "single already-generated RFDiffusion design and carries them"
    echo "through QC / negative-steering / orthogonal validation."
    echo ""
    echo "Params file must specify:"
    echo "  source_run_dir:    the ORIGINAL run's results/ dir that produced"
    echo "                     the design being resampled"
    echo "  source_design_id:  RFDiffusion design stem, e.g. 'design_26'"
    echo "  pdb_file:          the SAME pre-docked complex PDB the source run used"
    echo "  contigs:           the SAME contig spec string the source run used"
    echo ""
    echo "See resample_mpnn.nf for the full parameter list."
    exit 1
fi

PARAMS_FILE="$(realpath "$1")"
shift

if [ ! -f "${PARAMS_FILE}" ]; then
    echo "ERROR: params file not found: ${PARAMS_FILE}"
    exit 1
fi

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="${RECEPTOR_PIPELINE_DIR:-$(cd "." && pwd)}"
NEXTFLOW_IMG="${PIPELINE_DIR}/containers/NextFlow.img"
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
echo "Receptor Resurfacing Pipeline — Resample-MPNN Launcher"
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

"${NEXTFLOW_BIN}" run "${PIPELINE_DIR}/resample_mpnn.nf" \
    -params-file "${PARAMS_FILE}" \
    --rfdiff_container "${RFDIFF_CONTAINER}" \
    --af2_data_dir "${AF2_DATA_DIR}" \
    -resume \
    "$@"

echo ""
echo "============================================================"
echo "Pipeline finished: $(date)"
echo "============================================================"
