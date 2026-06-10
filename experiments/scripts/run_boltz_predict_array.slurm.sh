#!/bin/bash
#SBATCH --job-name="boltz_predict_array"
#SBATCH -p jic-gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
# Per-array-task stdout/stderr land in OUTDIR/logs/ — see redirect
# below; --output / --error in the SBATCH header are deliberately
# omitted so we control the path after CLI parsing.
#
# Resources sourced verbatim from
# experiments/scripts/nextflow.config:110-117 (BOLTZ2_MSA):
#   queue          = 'jic-gpu'         → -p jic-gpu
#   cpus           = 10                → -c 10
#   memory         = '32 GB'           → --mem=32G
#   time           = '12h'             → --time=12:00:00
#   clusterOptions = '--gres=gpu:1'    → --gres=gpu:1
#
# Submit AFTER the MSA stage:
#   #SBATCH --dependency=afterok:${MSA_JOB_ID}
# (added on the command line — see the chained-submission usage block
# at the top of run_msa_search_array.slurm.sh).

set -euo pipefail

# ── Usage ─────────────────────────────────────────────────────────────
if [ "$#" -lt 4 ]; then
    cat >&2 <<'EOF'
Usage:
  sbatch --array=0-N%K --dependency=afterok:${MSA_JOB} \
      run_boltz_predict_array.slurm.sh \
      MANIFEST_PATH REFERENCE_PDB OUTDIR MSA_STAGEDIR

Where:
  MANIFEST_PATH   Same manifest used for the MSA stage (one absolute
                  negsteer workdir path per line; lines starting with
                  # are skipped).  $SLURM_ARRAY_TASK_ID picks the line.
  REFERENCE_PDB   Ground-truth complex PDB (used by binding_rmsds and
                  by extract_effector_template_cif).
  OUTDIR          Root output directory for the predict stage.  Each
                  task writes into OUTDIR/<workdir_basename>/.
  MSA_STAGEDIR    Same staging dir given to run_msa_search_array.slurm.sh.
                  This task reads MSA_STAGEDIR/<design_id>/msa_done.txt.

Chained submission (full pipeline):
  N=$(grep -cv '^#' MANIFEST_PATH)
  MSA_JOB=$(sbatch --parsable --array=0-$((N-1))%8 \
      run_msa_search_array.slurm.sh \
      MANIFEST_PATH MSA_STAGEDIR COLABFOLD_DB)
  sbatch --array=0-$((N-1))%4 --dependency=afterok:${MSA_JOB} \
      run_boltz_predict_array.slurm.sh \
      MANIFEST_PATH REFERENCE_PDB OUTDIR MSA_STAGEDIR
EOF
    exit 2
fi

MANIFEST_PATH="$(realpath "$1")"
REFERENCE_PDB="$(realpath "$2")"
OUTDIR="$(realpath "$3")"
MSA_STAGEDIR="$(realpath "$4")"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "ERROR: manifest not found: ${MANIFEST_PATH}" >&2
    exit 1
fi
if [ ! -f "${REFERENCE_PDB}" ]; then
    echo "ERROR: reference PDB not found: ${REFERENCE_PDB}" >&2
    exit 1
fi
if [ ! -d "${MSA_STAGEDIR}" ]; then
    echo "ERROR: MSA stagedir not found: ${MSA_STAGEDIR}" >&2
    exit 1
fi
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: \$SLURM_ARRAY_TASK_ID is unset — submit with --array=…" >&2
    exit 1
fi

mkdir -p "${OUTDIR}/logs"
LOG_OUT="${OUTDIR}/logs/predict_task_${SLURM_ARRAY_TASK_ID}.out"
LOG_ERR="${OUTDIR}/logs/predict_task_${SLURM_ARRAY_TASK_ID}.err"
exec > >(tee -a "${LOG_OUT}") 2> >(tee -a "${LOG_ERR}" >&2)

# ── Pick this task's manifest line ────────────────────────────────────
WORKDIR="$(awk -v i="${SLURM_ARRAY_TASK_ID}" '
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    {if (n++ == i) {print; exit}}
' "${MANIFEST_PATH}")"

if [ -z "${WORKDIR}" ]; then
    echo "ERROR: manifest ${MANIFEST_PATH} has no entry at index "\
"${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

DESIGN_ID="$(basename "${WORKDIR}")"
TASK_OUTDIR="${OUTDIR}/${DESIGN_ID}"
TASK_MSA_DIR="${MSA_STAGEDIR}/${DESIGN_ID}"
MSA_DONE="${TASK_MSA_DIR}/msa_done.txt"

mkdir -p "${TASK_OUTDIR}"

# ── Verify the MSA stage produced an A3M for this design ─────────────
if [ ! -f "${MSA_DONE}" ]; then
    echo "ERROR: msa_done.txt absent for ${DESIGN_ID} at ${MSA_DONE}." >&2
    echo "       The MSA search has not completed for this design." >&2
    exit 1
fi
MSA_A3M="$(head -n 1 "${MSA_DONE}" | tr -d '[:space:]')"
if [ ! -f "${MSA_A3M}" ]; then
    echo "ERROR: msa_done.txt at ${MSA_DONE} points at ${MSA_A3M},"\
" but that file does not exist." >&2
    exit 1
fi

# ── Repo + script paths ───────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

# Boltz2 container symlink, mirroring params.boltz2_container in nextflow.config.
BOLTZ2_CONTAINER="${REPO_ROOT}/containers/boltz2_negsteer.img"

PY_SCRIPT="${SCRIPT_DIR}/boltz2_msa_predict.py"
if [ ! -f "${PY_SCRIPT}" ]; then
    echo "ERROR: ${PY_SCRIPT} not found" >&2
    exit 1
fi

# ── Banner ────────────────────────────────────────────────────────────
echo "──────────────────────────────────────────────────────────────"
echo "Boltz predict array task ${SLURM_ARRAY_TASK_ID}"
echo "──────────────────────────────────────────────────────────────"
echo "  Date:           $(date)"
echo "  Node:           $(hostname)"
echo "  Manifest:       ${MANIFEST_PATH}"
echo "  Workdir:        ${WORKDIR}"
echo "  Design id:      ${DESIGN_ID}"
echo "  Reference PDB:  ${REFERENCE_PDB}"
echo "  Outdir:         ${TASK_OUTDIR}"
echo "  MSA A3M:        ${MSA_A3M}"
echo "  Repo root:      ${REPO_ROOT}"
echo "  Container:      ${BOLTZ2_CONTAINER}"
echo "──────────────────────────────────────────────────────────────"

singularity exec --nv \
    --bind "${REPO_ROOT}:${REPO_ROOT}" \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "$(dirname "${REFERENCE_PDB}"):$(dirname "${REFERENCE_PDB}")" \
    --bind "${OUTDIR}:${OUTDIR}" \
    --bind "${MSA_STAGEDIR}:${MSA_STAGEDIR}" \
    "${BOLTZ2_CONTAINER}" \
    python "${PY_SCRIPT}" \
        --stage         predict \
        --workdir       "${WORKDIR}" \
        --reference-pdb "${REFERENCE_PDB}" \
        --outdir        "${TASK_OUTDIR}" \
        --msa-a3m       "${MSA_A3M}" \
        --threads       "${SLURM_CPUS_PER_TASK:-10}"
RC=$?

echo ""
echo "──────────────────────────────────────────────────────────────"
echo "Boltz predict array task ${SLURM_ARRAY_TASK_ID} done"
echo "  exit code:  ${RC}"
echo "  date:       $(date)"
echo "──────────────────────────────────────────────────────────────"

exit "${RC}"
