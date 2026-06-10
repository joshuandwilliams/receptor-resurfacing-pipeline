#!/bin/bash
#SBATCH --job-name="rr_env_check"
#SBATCH -p jic-gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --gres=gpu:1
#SBATCH --output=env_check_%j.out
#SBATCH --error=env_check_%j.err
#
# =============================================================================
# Environment / software pre-flight check for the receptor-resurfacing pipeline
# =============================================================================
#
# Runs on the GPU queue so it can ALSO verify GPU access — both on the host
# (nvidia-smi) and from inside a container (torch.cuda.is_available() with
# `singularity exec --nv`).
#
# It verifies that every piece of software the pipeline relies on is installed
# and reachable through the containers/ symlinks (see containers/README.md):
#   - Singularity / SLURM tooling on the host
#   - GPU driver + GPU passthrough into a container
#   - each container symlink resolves to a runnable image
#   - representative tools inside each image (torch, java, pytest, ...)
#   - the OpenJDK symlink used as JAVA_HOME for Nextflow
#   - the Nextflow runner image
#   - reference-data directories (AlphaFold DBs, AF3 model dir)  [advisory]
#   - the native AlphaFold3 `source package`                     [advisory]
#
# Usage:
#   sbatch scripts/check_environment.slurm.sh           # submit from repo root
#   PIPELINE_DIR=/path/to/repo sbatch scripts/check_environment.slurm.sh
#
# Exit status: 0 if all REQUIRED checks pass; 1 if any required check fails.
# Advisory [WARN] checks never change the exit status.
# =============================================================================

set -uo pipefail   # NOT -e: we want to run every check and summarise at the end

# ── Locate the pipeline directory ─────────────────────────────────────────
# Prefer an explicit PIPELINE_DIR, then the dir sbatch was launched from,
# then the current working dir.  Validate by looking for nextflow.config.
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
if [ ! -f "${PIPELINE_DIR}/nextflow.config" ]; then
    echo "ERROR: '${PIPELINE_DIR}' does not look like the pipeline root" >&2
    echo "       (no nextflow.config there)." >&2
    echo "       Submit from the repo root, or set PIPELINE_DIR explicitly:" >&2
    echo "         PIPELINE_DIR=/path/to/repo sbatch scripts/check_environment.slurm.sh" >&2
    exit 2
fi
PIPELINE_DIR="$(cd "${PIPELINE_DIR}" && pwd)"
CONTAINERS_DIR="${PIPELINE_DIR}/containers"

# Reference data + AF3 native package (mirror nextflow.config).
AF2_DATA_DIR="/nbi/Reference-Data/AlphaFold/db-v2.3.2"
AF3_DB_V3="/nbi/Reference-Data/AlphaFold/db-v3.0.0"
AF3_MODEL_DIR="${HOME}/singularity/AlphaFold3"
AF3_PACKAGE_ID="e8edb411-7374-4342-b9f1-408da41fc197"

# ── Result tallies + helpers ──────────────────────────────────────────────
PASS=0; FAIL=0; WARN=0
ok()   { printf '  [ OK ] %s\n' "$*"; PASS=$((PASS+1)); }
bad()  { printf '  [FAIL] %s\n' "$*"; FAIL=$((FAIL+1)); }
warn() { printf '  [WARN] %s\n' "$*"; WARN=$((WARN+1)); }
hdr()  { printf '\n── %s  [+%ss] ──────────────────────────────\n' "$*" "${SECONDS}"; }

# ── Banner ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Receptor Resurfacing Pipeline — environment check"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Containers:    ${CONTAINERS_DIR}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "Job ID:        ${SLURM_JOB_ID:-<none>}"
echo "============================================================"

# ── Host tooling ──────────────────────────────────────────────────────────
hdr "Host tooling"
if command -v singularity >/dev/null 2>&1; then
    ok "singularity present: $(singularity --version 2>&1 | head -1)"
    HAVE_SINGULARITY=1
else
    bad "singularity not found on PATH — no container will run"
    HAVE_SINGULARITY=0
fi
for tool in sbatch squeue scontrol; do
    if command -v "${tool}" >/dev/null 2>&1; then
        ok "SLURM '${tool}' present"
    else
        warn "SLURM '${tool}' not on PATH"
    fi
done

# ── GPU on the host ───────────────────────────────────────────────────────
hdr "GPU (host)"
if command -v nvidia-smi >/dev/null 2>&1; then
    if nvidia-smi -L >/dev/null 2>&1; then
        ok "nvidia-smi sees: $(nvidia-smi -L 2>/dev/null | head -1)"
    else
        bad "nvidia-smi present but no GPU visible (allocation issue?)"
    fi
else
    bad "nvidia-smi not found — is this running on the GPU queue?"
fi
echo "  (CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}, " \
     "SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-<unset>})"

# ── Container symlinks resolve + images load ──────────────────────────────
# filename in containers/  →  human label
CONTAINERS=(
    "NextFlow.img|Nextflow runner"
    "HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img|RFDiffusion / ProteinMPNN / MMseqs2"
    "Rosetta.img|Rosetta"
    "boltz2_negsteer.img|Boltz-2 / negative steering"
    "colabfold.img|ColabFold MSA"
    "pytest_runner.img|pytest characterisation suite"
)
hdr "Container symlinks resolve"
for entry in "${CONTAINERS[@]}"; do
    f="${entry%%|*}"; label="${entry##*|}"
    path="${CONTAINERS_DIR}/${f}"
    if [ -e "${path}" ]; then
        if [ -L "${path}" ]; then
            ok "${f} → $(readlink "${path}")  (${label})"
        else
            ok "${f} present (regular file)  (${label})"
        fi
    elif [ -L "${path}" ]; then
        bad "${f}: dangling symlink → $(readlink "${path}")  (${label})"
    else
        bad "${f}: missing from containers/  (${label})"
    fi
done

# ── Container images: one parallel exec per image ─────────────────────────
# Cold-starting a multi-GB .img over the networked /hpc-home mount is the
# slow part, so each image is loaded EXACTLY ONCE (all of its checks folded
# into a single exec) and all images run CONCURRENTLY — wall-time is the
# slowest single image, not the sum.  Every exec is wrapped in `timeout` so
# a hung mount surfaces as a FAIL instead of stalling the whole job.
# The boltz2 image is run with --nv and doubles as the in-container GPU
# check (torch.cuda.is_available()).
hdr "Container images (parallel, one exec each)"
EXEC_TIMEOUT=420
if [ "${HAVE_SINGULARITY}" -ne 1 ]; then
    warn "skipping all in-container checks (no singularity)"
else
    CRES="$(mktemp -d)"
    # spawn <order> <required|optional> <label> <nv:0|1> <image-file> <cmd...>
    spawn() {
        local order="$1" sev="$2" label="$3" nv="$4" f="$5"; shift 5
        local p="${CONTAINERS_DIR}/${f}" nvflag=""
        [ "${nv}" = "1" ] && nvflag="--nv"
        (
            if [ ! -e "${p}" ]; then
                printf 'WARN|%s|%s\n' "${label}" "image not available"
            elif detail="$(timeout "${EXEC_TIMEOUT}" singularity exec ${nvflag} "${p}" "$@" 2>/dev/null)"; then
                printf 'OK|%s|%s\n' "${label}" "${detail}"
            elif [ $? -eq 124 ]; then
                printf '%s|%s|%s\n' "$([ "${sev}" = optional ] && echo WARN || echo FAIL)" \
                    "${label}" "timed out after ${EXEC_TIMEOUT}s"
            elif [ "${sev}" = optional ]; then
                printf 'WARN|%s|%s\n' "${label}" "not found"
            else
                printf 'FAIL|%s|%s\n' "${label}" "check failed"
            fi
        ) > "${CRES}/${order}" &
    }

    spawn 1 required "NextFlow.img (java + nextflow dist)" 0 NextFlow.img \
        bash -lc 'v=$(java -version 2>&1 | head -1); test -f /usr/local/bin/nextflow && printf "%s + nextflow dist" "$v"'
    spawn 2 required "HADDOCK img (torch)" 0 HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img \
        python -c 'import torch; print("torch", torch.__version__)'
    spawn 3 required "boltz2 img (torch + GPU via --nv)" 1 boltz2_negsteer.img \
        bash -lc 'command -v boltz >/dev/null || { echo "boltz CLI missing"; exit 1; }; python -c "import torch,sys; ok=torch.cuda.is_available(); print(\"torch\", torch.__version__, \"| cuda\", ok, \"|\", torch.cuda.get_device_name(0) if ok else \"no GPU\"); sys.exit(0 if ok else 1)"'
    spawn 4 required "pytest img (pytest)" 0 pytest_runner.img \
        python -c 'import pytest; print("pytest", pytest.__version__)'
    spawn 5 optional "Rosetta img (rosetta app on PATH)" 0 Rosetta.img \
        bash -lc 'command -v relax.default.linuxgccrelease || command -v score_jd2.default.linuxgccrelease'
    spawn 6 optional "colabfold img (colabfold_search/mmseqs)" 0 colabfold.img \
        bash -lc 'command -v colabfold_search || command -v mmseqs'

    wait
    for o in 1 2 3 4 5 6; do
        IFS='|' read -r st label detail < "${CRES}/${o}"
        case "${st}" in
            OK)   ok   "${label} — ${detail}" ;;
            WARN) warn "${label}: ${detail}" ;;
            *)    bad  "${label}: ${detail}" ;;
        esac
    done
    rm -rf "${CRES}"
fi

# ── OpenJDK symlink (JAVA_HOME for Nextflow) ──────────────────────────────
hdr "OpenJDK (JAVA_HOME)"
JDK_DIR="${CONTAINERS_DIR}/jdk-17.0.2"
if [ -e "${JDK_DIR}/bin/java" ]; then
    if "${JDK_DIR}/bin/java" -version >/dev/null 2>&1; then
        ok "jdk-17.0.2 → $("${JDK_DIR}/bin/java" -version 2>&1 | head -1)"
    else
        bad "jdk-17.0.2/bin/java present but '-version' failed"
    fi
elif [ -L "${JDK_DIR}" ]; then
    bad "jdk-17.0.2: dangling symlink → $(readlink "${JDK_DIR}")"
else
    bad "jdk-17.0.2: missing from containers/"
fi

# ── Reference data + native AF3 (advisory) ────────────────────────────────
hdr "Reference data + AlphaFold3 (advisory)"
for d in "${AF2_DATA_DIR}" "${AF3_DB_V3}" "${AF3_MODEL_DIR}"; do
    if [ -d "${d}" ]; then ok "present: ${d}"; else warn "missing: ${d}"; fi
done
# `source package` is a shell builtin so it can't be `timeout`-wrapped
# directly — run it inside a time-bounded login subshell instead.
if timeout 120 bash -lc "source package ${AF3_PACKAGE_ID} >/dev/null 2>&1 && command -v run_alphafold >/dev/null 2>&1"; then
    ok "AF3 native package loads (run_alphafold on PATH)"
else
    warn "AF3 native package '${AF3_PACKAGE_ID}' not loadable here (only needed for AF3 orthogonal validation)"
fi

# ── Summary ───────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "Summary:  ${PASS} OK   ${FAIL} FAIL   ${WARN} WARN"
echo "============================================================"
if [ "${FAIL}" -gt 0 ]; then
    echo "RESULT: FAIL — ${FAIL} required check(s) failed (see [FAIL] above)."
    exit 1
fi
echo "RESULT: PASS — all required software is installed and reachable."
[ "${WARN}" -gt 0 ] && echo "(${WARN} advisory warning(s) — review if you need those features.)"
exit 0
