#!/bin/bash
#SBATCH --job-name="nf_test_pose_solver"
# No -p: the partition is site-specific.  tests/run_tests.sh passes
# RECEPTOR_PARTITION_CPU from site.conf as `sbatch -p`, where it beats both
# the environment and a header line.  Submitting this script directly takes
# the cluster's default partition.
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user: sbatch(1) mails the submitting user when only --mail-type
# is given, so an address here would be site-specific and redundant.
# No --chdir: the path is site-specific.  tests/run_tests.sh passes it as
# `sbatch --chdir`, which is what sets the $PWD this script resolves its
# inputs from.  Submit from this directory if you run it by hand.

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="${RECEPTOR_PIPELINE_DIR:-$(cd "../.." && pwd)}"
NEXTFLOW_IMG="${PIPELINE_DIR}/containers/NextFlow.img"

TEST_DIR="${PIPELINE_DIR}/tests/pose_solver"
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
echo "Pose Solver Module Test — Nextflow Launcher"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Test dir:      ${TEST_DIR}"
echo "NXF_WORK:      ${NXF_WORK}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

"${NEXTFLOW_BIN}" run "${TEST_DIR}/test_pose_solver.nf" \
    -c "${PIPELINE_DIR}/nextflow.config" \
    -resume \
    "$@"

echo ""
echo "============================================================"
echo "Test finished: $(date)"
echo "============================================================"
