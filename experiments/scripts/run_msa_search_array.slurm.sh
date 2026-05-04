#!/bin/bash
#SBATCH --job-name="msa_search_array"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
# Per-array-task stdout/stderr land in MSA_STAGEDIR/logs/ — see the
# redirect below; the SBATCH --output / --error directives are
# intentionally omitted so we can derive the path after parsing CLI args.
#
# Resources sourced verbatim from
# experiments/scripts/nextflow.config:65-70 (COLABFOLD_SEARCH):
#   queue  = 'jic-medium'      → -p jic-medium
#   cpus   = 20                → -c 20
#   memory = '128 GB'          → --mem=128G
#   time   = '4h'              → --time=4:00:00
# No GPU.

set -euo pipefail

# ── Usage ─────────────────────────────────────────────────────────────
if [ "$#" -lt 3 ]; then
    cat >&2 <<'EOF'
Usage:
  sbatch --array=0-N%K run_msa_search_array.slurm.sh \
      MANIFEST_PATH MSA_STAGEDIR COLABFOLD_DB

Where:
  MANIFEST_PATH   Plain text manifest, one absolute negsteer workdir
                  path per line (lines starting with # are skipped).
                  $SLURM_ARRAY_TASK_ID selects the line (0-indexed).
                  Build with build_msa_manifest.py.
  MSA_STAGEDIR    Root staging dir for MSA outputs.  Each task writes
                  MSA_STAGEDIR/<workdir_basename>/msa/<...>_receptor.a3m
                  and MSA_STAGEDIR/<workdir_basename>/msa_done.txt.
  COLABFOLD_DB    ColabFold MMseqs2 sequence database directory.

Chain with the predict stage:
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
MSA_STAGEDIR="$(realpath "$2")"
COLABFOLD_DB="$(realpath "$3")"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "ERROR: manifest not found: ${MANIFEST_PATH}" >&2
    exit 1
fi
if [ ! -d "${COLABFOLD_DB}" ]; then
    echo "ERROR: ColabFold DB directory not found: ${COLABFOLD_DB}" >&2
    exit 1
fi
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: \$SLURM_ARRAY_TASK_ID is unset — submit with --array=…" >&2
    exit 1
fi

mkdir -p "${MSA_STAGEDIR}/logs"
LOG_OUT="${MSA_STAGEDIR}/logs/msa_task_${SLURM_ARRAY_TASK_ID}.out"
LOG_ERR="${MSA_STAGEDIR}/logs/msa_task_${SLURM_ARRAY_TASK_ID}.err"
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
TASK_OUTDIR="${MSA_STAGEDIR}/${DESIGN_ID}"
mkdir -p "${TASK_OUTDIR}"

# ── Repo + script paths ───────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

# Boltz2 container is the canonical runtime for resolve_steered_inputs
# (it imports bin/boltz2_negative_steering, which needs BioPython /
# numpy / gemmi from the boltz2_negsteer image).  The MSA itself is
# run via a nested singularity exec by the Python script, against the
# colabfold container — boltz2_msa_predict.py composes that call from
# its HPC CONSTANTS block.  Path mirrors nextflow.config line 21.
BOLTZ2_CONTAINER="/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/boltz2_negsteer.img"

PY_SCRIPT="${SCRIPT_DIR}/boltz2_msa_predict.py"
if [ ! -f "${PY_SCRIPT}" ]; then
    echo "ERROR: ${PY_SCRIPT} not found" >&2
    exit 1
fi

# ── Banner ────────────────────────────────────────────────────────────
echo "──────────────────────────────────────────────────────────────"
echo "MSA search array task ${SLURM_ARRAY_TASK_ID}"
echo "──────────────────────────────────────────────────────────────"
echo "  Date:           $(date)"
echo "  Node:           $(hostname)"
echo "  Manifest:       ${MANIFEST_PATH}"
echo "  Workdir:        ${WORKDIR}"
echo "  Design id:      ${DESIGN_ID}"
echo "  MSA stagedir:   ${TASK_OUTDIR}"
echo "  ColabFold DB:   ${COLABFOLD_DB}"
echo "  Repo root:      ${REPO_ROOT}"
echo "  Container:      ${BOLTZ2_CONTAINER}"
echo "──────────────────────────────────────────────────────────────"

singularity exec \
    --bind "${REPO_ROOT}:${REPO_ROOT}" \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "${MSA_STAGEDIR}:${MSA_STAGEDIR}" \
    --bind "${COLABFOLD_DB}:${COLABFOLD_DB}" \
    "${BOLTZ2_CONTAINER}" \
    python "${PY_SCRIPT}" \
        --stage         msa \
        --workdir       "${WORKDIR}" \
        --outdir        "${TASK_OUTDIR}" \
        --colabfold-db  "${COLABFOLD_DB}" \
        --threads       "${SLURM_CPUS_PER_TASK:-20}"
RC=$?

# Final invariant: msa_done.txt must exist on success.  If the Python
# script reported success but the file isn't there, treat it as a
# failure so the predict-stage dependency does NOT advance.
if [ "${RC}" -eq 0 ] && [ ! -f "${TASK_OUTDIR}/msa_done.txt" ]; then
    echo "ERROR: ${TASK_OUTDIR}/msa_done.txt missing after a clean exit." >&2
    RC=1
fi

echo ""
echo "──────────────────────────────────────────────────────────────"
echo "MSA search array task ${SLURM_ARRAY_TASK_ID} done"
echo "  exit code:  ${RC}"
echo "  date:       $(date)"
echo "──────────────────────────────────────────────────────────────"

exit "${RC}"
