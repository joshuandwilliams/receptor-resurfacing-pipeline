#!/bin/bash
#
# capture_reference_set.sh — replace example_output_files/ with the latest
# receptor_resurfacing_results/ from a test_rfdiffusion.nf run.
#
# Run this ON THE HPC after run_test_rfdiffusion.slurm.sh completes:
#
#   bash tests/rfdiffusion/capture_reference_set.sh
#
# Then on Mac, copy back from the mounted volume:
#
#   cp -r /Volumes/HPC-Home/receptor_design/receptor-resurfacing-pipeline/tests/rfdiffusion/example_output_files/ \
#        tests/rfdiffusion/example_output_files/
#
# Then commit, push, and sync:
#
#   git add tests/rfdiffusion/example_output_files/
#   git commit -m "tests/rfdiffusion: regenerate reference set — Complex_beta_ckpt.pt + pikp1_avrpikf_complex.pdb"
#   git push
#   ./scripts/sync_to_hpc.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="${SCRIPT_DIR}/receptor_resurfacing_results"
REF="${SCRIPT_DIR}/example_output_files"

if [ ! -d "${RESULTS}" ]; then
    echo "ERROR: receptor_resurfacing_results/ not found at ${RESULTS}" >&2
    echo "       Run run_test_rfdiffusion.slurm.sh first." >&2
    exit 1
fi

# Sanity-check that the key output files are present before overwriting ref.
for f in rfdiffusion/rfdiffusion_metrics.json rfdiffusion/filter_summary.json rfdiffusion/passing_designs.txt; do
    if [ ! -f "${RESULTS}/${f}" ]; then
        echo "ERROR: expected output file not found: ${RESULTS}/${f}" >&2
        echo "       The test run may not have completed successfully." >&2
        exit 1
    fi
done

echo "Replacing example_output_files/ with new reference set..."
rm -rf "${REF}"
cp -r "${RESULTS}" "${REF}"

echo ""
echo "Done. New reference set is at: ${REF}"
echo ""
echo "Next steps on Mac:"
echo "  cp -r /Volumes/HPC-Home/receptor_design/receptor-resurfacing-pipeline/tests/rfdiffusion/example_output_files/ \\"
echo "       tests/rfdiffusion/example_output_files/"
echo "  git add tests/rfdiffusion/example_output_files/"
echo "  git commit -m \"tests/rfdiffusion: regenerate reference set — Complex_beta_ckpt.pt + pikp1_avrpikf_complex.pdb\""
echo "  git push"
echo "  ./scripts/sync_to_hpc.sh"
