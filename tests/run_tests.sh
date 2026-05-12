#!/bin/bash
# -----------------------------------------------------------------------------
# tests/run_tests.sh
# -----------------------------------------------------------------------------
# Fan-out dispatcher: submit per-module Nextflow test SLURM jobs in parallel.
# Runs on the HPC login node (NOT itself a SLURM job — just calls sbatch
# repeatedly and exits).
#
# Usage:
#   ./tests/run_tests.sh --modules MOD [MOD ...] [--with-plots] [--dry-run]
#
# Examples:
#   ./tests/run_tests.sh --modules rfdiffusion proteinmpnn
#   ./tests/run_tests.sh --modules negative_steering --with-plots
#   ./tests/run_tests.sh --modules rfdiffusion --dry-run        # print, don't submit
#
# Recognised modules: rfdiffusion, proteinmpnn, rosetta_filtering,
# negative_steering, orthogonal_metrics, haddock.
#
# --with-plots adds the per-module plot SLURM script as a dependency
# (afterok:<workflow-job-id>) so plots only run if the workflow test
# succeeds.  Modules without a plot script (e.g. haddock) silently skip.
#
# Each module is its own independent sbatch — they run in parallel
# subject to cluster capacity.  Output / error logs land in each
# module's tests/<module>/ directory (slurm_<jobid>.out / .err).
# -----------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VALID_MODULES=(
    rfdiffusion
    proteinmpnn
    rosetta_filtering
    negative_steering
    orthogonal_metrics
    haddock
)

MODULES=()
WITH_PLOTS=0
DRY_RUN=0

usage() {
    cat <<EOF
Usage: $(basename "$0") --modules MOD [MOD ...] [--with-plots] [--dry-run]

  --modules     One or more of: ${VALID_MODULES[*]}
  --with-plots  Also submit the per-module plot SLURM script with an
                afterok:<workflow-job> dependency.
  --dry-run     Print what would be submitted, but do not call sbatch.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --modules)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                MODULES+=("$1"); shift
            done
            ;;
        --with-plots)  WITH_PLOTS=1; shift ;;
        --dry-run)     DRY_RUN=1;    shift ;;
        -h|--help)     usage; exit 0 ;;
        *)             echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [ ${#MODULES[@]} -eq 0 ]; then
    echo "ERROR: --modules requires at least one module" >&2
    usage
    exit 2
fi

# Validate every requested module name.
for m in "${MODULES[@]}"; do
    if [[ ! " ${VALID_MODULES[*]} " == *" $m "* ]]; then
        echo "ERROR: unknown module: $m" >&2
        echo "       valid: ${VALID_MODULES[*]}" >&2
        exit 2
    fi
done

submit_one() {
    local module="$1"
    local mod_dir="${SCRIPT_DIR}/${module}"
    local workflow_script="${mod_dir}/run_test_${module}.slurm.sh"

    if [ ! -f "${workflow_script}" ]; then
        echo "[${module}] SKIP: workflow script not found: ${workflow_script}" >&2
        return
    fi

    echo "[${module}] sbatch ${workflow_script}"
    if [ "${DRY_RUN}" = "1" ]; then
        local workflow_jid="DRY_RUN_JOBID"
    else
        # sbatch with --parsable emits just the job id on stdout, which
        # makes the dependency wiring below trivial.
        local workflow_jid
        workflow_jid=$(sbatch --parsable "${workflow_script}")
        echo "[${module}]   workflow job id: ${workflow_jid}"
    fi

    if [ "${WITH_PLOTS}" = "1" ]; then
        # Plot test SLURM scripts vary by module name:
        #   tests/<module>/run_test_<module>_plots.slurm.sh
        # except a couple where the naming differs slightly — list them
        # all here.  Modules without a plot script (haddock) silently skip.
        local plot_scripts=()
        case "${module}" in
            rfdiffusion|proteinmpnn|rosetta_filtering)
                # Singular short name: run_test_<short>_plots.slurm.sh
                local short="${module/proteinmpnn/mpnn}"
                short="${short/rosetta_filtering/rosetta_filtering}"
                plot_scripts+=("${mod_dir}/run_test_${short}_plots.slurm.sh")
                ;;
            negative_steering)
                plot_scripts+=("${mod_dir}/run_test_negsteer_plots.slurm.sh")
                plot_scripts+=("${mod_dir}/run_test_negsteer_within_sequence_plots.slurm.sh")
                ;;
            orthogonal_metrics)
                plot_scripts+=("${mod_dir}/run_test_orthogonal_metrics_plot.slurm.sh")
                ;;
            haddock)
                # No plot tests yet.
                ;;
        esac

        for plot_script in "${plot_scripts[@]:-}"; do
            # Skip the empty-array sentinel emitted by `${arr[@]:-}` under set -u.
            [ -z "${plot_script}" ] && continue
            if [ ! -f "${plot_script}" ]; then
                echo "[${module}]   SKIP plot test: ${plot_script} (not found)" >&2
                continue
            fi
            echo "[${module}]   --with-plots: sbatch ${plot_script} (afterok:${workflow_jid})"
            if [ "${DRY_RUN}" != "1" ]; then
                local plot_jid
                plot_jid=$(sbatch --parsable \
                    --dependency=afterok:"${workflow_jid}" \
                    "${plot_script}")
                echo "[${module}]   plot job id:     ${plot_jid}"
            fi
        done
    fi
}

echo "============================================================"
echo "run_tests.sh"
echo "  Modules:    ${MODULES[*]}"
echo "  With plots: ${WITH_PLOTS}"
echo "  Dry-run:    ${DRY_RUN}"
echo "============================================================"

for m in "${MODULES[@]}"; do
    submit_one "$m"
done

echo "============================================================"
if [ "${DRY_RUN}" = "1" ]; then
    echo "Dry-run complete — nothing was submitted."
else
    echo "All requested jobs submitted.  Monitor with: squeue -u \$USER"
    echo "When a module finishes, refresh its fixture with:"
    echo "  sbatch tests/update_example_dataset.slurm.sh \\"
    echo "        --module <module> \\"
    echo "        --updated-output-folder tests/<module>/receptor_resurfacing_results"
fi
echo "============================================================"
