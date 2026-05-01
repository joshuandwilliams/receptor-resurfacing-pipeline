/*
 * =============================================================================
 * modules/negsteer_manifest.nf
 *
 * P0-31 · Per-survivor manifest extraction.
 *
 * Single CPU process EXTRACT_SURVIVOR_MANIFEST that reads the P0-29
 * extended cross-sequence CSV (cross_sequence_summary_with_interface_metrics.csv)
 * plus every per-sequence workdir's plan.json, and emits a manifest CSV
 * keyed by seq_name with the per-survivor inputs the three orthogonal-
 * metrics streams need:
 *
 *   seq_name
 *   canonical_pdb_abs              — Boltz canonical prediction PDB
 *   ground_truth_abs               — RFDiffusion ground-truth PDB
 *   effector_template_cif_abs      — effector template CIF (reused from Boltz)
 *   receptor_seq, effector_seq     — chain sequences extracted from the PDB
 *
 * The downstream workflow consumes this CSV via splitCsv(header: true)
 * to fan out one tuple per survivor into AF3, biophysical, and Rosetta
 * streams.
 *
 * Promoted out of test_orthogonal_metrics.nf into this module by P0-34
 * (todo_list6) so the production main.nf and the test harness import the
 * same code.  The bin/extract_survivor_manifest.py script is unchanged.
 *
 * Runs inside boltz2_negsteer.img (uses the same numpy/gemmi/biopython
 * stack as the other orthogonal-metrics processes).  ~1 minute wall on
 * cohorts up to ~200 rows.
 * =============================================================================
 */

process EXTRACT_SURVIVOR_MANIFEST {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}/orthogonal_metrics", mode: 'copy',
        pattern: "survivor_manifest.csv"

    input:
    path extended_csv
    path workdirs, stageAs: 'workdirs/*'
    path manifest_script

    output:
    path "survivor_manifest.csv", emit: manifest

    script:
    """
    set -euo pipefail

    singularity exec \\
            --bind \${PWD}:\${PWD} \\
            --bind ${projectDir}:${projectDir} \\
            ${params.boltz2_container} \\
        python ${manifest_script} \\
            --input-csv ${extended_csv} \\
            --workdirs-glob 'workdirs/*' \\
            --receptor-chain ${params.rfdiff_output_receptor_chain} \\
            --effector-chain ${params.rfdiff_output_effector_chain} \\
            --output-manifest survivor_manifest.csv

    echo "Manifest built:"
    head -1 survivor_manifest.csv
    echo "(row count: \$(tail -n +2 survivor_manifest.csv | wc -l))"
    """
}
