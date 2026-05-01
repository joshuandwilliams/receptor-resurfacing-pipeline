/*
 * =============================================================================
 * modules/negsteer_biophysical_metrics.nf
 *
 * P0-31 · Biophysical orthogonal metrics on survivors.
 *
 * One CPU process per survivor, fanned out in parallel with the AF3 and
 * Rosetta streams.  Three metrics:
 *
 *   bsa                 Buried surface area (Å²), via FreeSASA.
 *                       Target: > 600 Å² (Overath et al. 2025).
 *   interface_plddt     Mean per-residue pLDDT over interface residues.
 *                       Reads the cached Boltz prediction's B-factor
 *                       column — no new prediction run.  Target > 0.75.
 *   interface_hbonds    Interface H-bond count via MDAnalysis
 *                       HydrogenBondAnalysis.
 *
 * Runs inside boltz2_negsteer.img (FreeSASA, MDAnalysis, gemmi, numpy
 * all pre-installed and pinned).  Input: one per-survivor tuple.
 * Output: single-row CSV with the three metrics plus a failures column.
 * =============================================================================
 */

process NEGSTEER_BIOPHYSICAL_METRICS {
    tag "${seq_name}"
    label 'cpu'

    publishDir "${params.outdir}/orthogonal_metrics/biophysical/${seq_name}",
        mode: 'copy', pattern: "biophysical_summary.csv"

    input:
    tuple val(seq_name),
          path(canonical_pdb),
          path(ground_truth_pdb)
    path metrics_script

    output:
    path "biophysical_summary.csv", emit: summary_csv

    script:
    """
    set -euo pipefail

    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python ${metrics_script} \\
            --seq-name ${seq_name} \\
            --canonical-pdb ${canonical_pdb} \\
            --ground-truth ${ground_truth_pdb} \\
            --receptor-chain ${params.rfdiff_output_receptor_chain} \\
            --effector-chain ${params.rfdiff_output_effector_chain} \\
            --contact-cutoff ${params.negsteer_postprocess_contact_cutoff} \\
            --output-csv biophysical_summary.csv
    """
}
