#!/bin/bash
#SBATCH --job-name=af3_unsteered
#SBATCH --partition=jic-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=1-22%22
#SBATCH --output=logs/af3_unsteered_%A_%a.out
#SBATCH --error=logs/af3_unsteered_%A_%a.err
#
# AlphaFold3 without an MSA on the ORIGINAL ProteinMPNN sequences of the 22
# designs the pipeline never steered.
#
# The pipeline's own AF3 pass ran on each survivor's steered sequence, so it
# cannot say whether AlphaFold3 and Boltz-2 agree on the designs Boltz-2 posed
# correctly with no steering at all. This fills that gap.
#
# The pipeline requested 20 cpus and 64 GB for this process, but its own trace
# shows 1m25s realtime at 7.8 GB peak RSS, so 8 cpus and 16 GB is ample and
# lets these tasks share partly-used nodes rather than waiting for empty ones.
#
# Usage, from this directory:
#   mkdir -p logs
#   sbatch af3_unsteered.slurm.sh

set -euo pipefail

# SLURM copies the script into /var/spool, so BASH_SOURCE points at the spool
# copy rather than the submission directory. SLURM_SUBMIT_DIR is the one that
# survives, with a fallback for running this by hand.
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
AF3_PACKAGE=e8edb411-7374-4342-b9f1-408da41fc197
AF3_MODEL_DIR="${HOME}/singularity/AlphaFold3"
AF3_DB_DIR="${HOME}/singularity/AlphaFold3/af3_db"

UNIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${HERE}/units.txt")
if [[ -z "${UNIT}" ]]; then
    echo "No unit at line ${SLURM_ARRAY_TASK_ID} of units.txt" >&2
    exit 1
fi

WORKDIR="${HERE}/output/${UNIT}"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# Same JAX memory settings the pipeline exports before every AF3 call.
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2

source package "${AF3_PACKAGE}"

echo "AF3-no-MSA, unsteered sequence, ${UNIT}"
run_alphafold.py \
    --json_path="${HERE}/inputs/${UNIT}.json" \
    --model_dir="${AF3_MODEL_DIR}" \
    --db_dir="${AF3_DB_DIR}" \
    --output_dir="${WORKDIR}/output" \
    --norun_data_pipeline

echo "done ${UNIT}"
find output -name "*.cif" | head -20
