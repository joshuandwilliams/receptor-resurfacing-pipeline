/*
 * =============================================================================
 * modules/negsteer_interface_metrics.nf
 *
 * P0-29 · Interface-restricted RMSD (iRMSD), fnat, DockQ, and 15 Å iPSAE.
 * P0-38 · Smooth distance-weighted Jaccard contact overlap.
 *
 * Single process NEGSTEER_INTERFACE_METRICS runs after NEGSTEER_CROSS_SEQUENCE.
 * Reads cross_sequence_summary.csv plus every per-sequence workdir (for
 * plan.json → ground-truth PDB resolution), and emits a new CSV with
 * eight additional metric columns appended per row:
 *
 *   irmsd, fnat, dockq                          (DockQ Python API)
 *   ipsae_ab_15, ipsae_ba_15, ipsae_min_15      (15 Å PAE cutoff)
 *   intact_core                                 (pLDDT-trimmed RMSD filter)
 *   weighted_jaccard                            (P0-38 — pair-level
 *                                                Cβ–Cβ Gaussian-weighted
 *                                                Jaccard, generalises
 *                                                the existing residue-set
 *                                                true_jaccard for Task 8
 *                                                AUROC comparison)
 *   interface_metrics_failures                  (comma-separated reasons)
 *
 * Per todo_list3 P0-29 acceptance criterion: rows where a metric cannot
 * be computed are flagged in the failures column, not silently dropped.
 *
 * Runs inside boltz2_negsteer.img (DockQ + numpy + gemmi already installed).
 * CPU-only, ~1-3 s per row, so runs as a single job rather than fanning out.
 * =============================================================================
 */

process NEGSTEER_INTERFACE_METRICS {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}/negative_steering", mode: 'copy',
        pattern: "cross_sequence_summary_with_interface_metrics.csv"

    input:
    path cross_sequence_summary_csv
    path per_sequence_workdirs, stageAs: 'workdirs/*'

    output:
    path "cross_sequence_summary_with_interface_metrics.csv",
         emit: extended_csv

    script:
    """
    set -euo pipefail

    singularity exec \\
            --bind \${PWD}:\${PWD} \\
            --bind ${projectDir}:${projectDir} \\
            ${params.boltz2_container} \\
        python ${projectDir}/bin/compute_interface_metrics.py \\
            --input-csv ${cross_sequence_summary_csv} \\
            --output-csv cross_sequence_summary_with_interface_metrics.csv \\
            --workdirs-glob 'workdirs/*' \\
            --plddt-threshold ${params.interface_plddt_trim_threshold} \\
            --intact-threshold ${params.interface_intact_threshold} \\
            --receptor-chain ${params.receptor_chain} \\
            --effector-chain ${params.effector_chain} \\
            --weighted-jaccard-pair-cutoff ${params.weighted_jaccard_pair_cutoff}

    echo "Interface metrics written."
    head -1 cross_sequence_summary_with_interface_metrics.csv
    echo "(row count: \$(tail -n +2 cross_sequence_summary_with_interface_metrics.csv | wc -l))"
    """
}
