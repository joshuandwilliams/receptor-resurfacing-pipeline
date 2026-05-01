/*
 * =============================================================================
 * modules/negsteer_orthogonal_metrics.nf
 *
 * P0-31 · Aggregate the three parallel orthogonal-metrics streams into
 * a single survivors_with_orthogonal_metrics.csv.
 *
 * Inputs:
 *   cross_sequence_summary_with_interface_metrics.csv  (from P0-29)
 *   af3_nomsa_summaries.csv      (concat of AF3_PARSE_OUTPUT rows)
 *   biophysical_summaries.csv    (concat of NEGSTEER_BIOPHYSICAL_METRICS rows)
 *   rosetta_summaries.csv        (concat of NEGSTEER_ROSETTA_METRICS rows)
 *
 * Output columns added to each row:
 *   af3_nomsa_best_ra_eff, af3_nomsa_mean_ra_eff, af3_nomsa_best_iptm,
 *   af3_nomsa_mean_iptm, af3_nomsa_n_correct_interface, af3_nomsa_failures,
 *   bsa, interface_plddt, interface_hbonds, biophysical_failures,
 *   sc, rosetta_ddg, rosetta_failures,
 *   orthogonal_flags, passes_orthogonal_filters
 *
 * Filter cascade:
 *   - af3_nomsa_best_ra_eff < orthogonal_filter_af3_ra_max → hard drop on fail
 *   - sc              >= orthogonal_filter_sc_min     → flag on fail
 *   - bsa             >= orthogonal_filter_bsa_min    → flag on fail
 *   - interface_plddt >= orthogonal_filter_plddt_min  → flag on fail
 *
 * passes_orthogonal_filters is 1 only when AF3 passes AND no flags present.
 * =============================================================================
 */

process NEGSTEER_ORTHOGONAL_METRICS {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}/orthogonal_metrics", mode: 'copy',
        pattern: "survivors_with_orthogonal_metrics.csv"

    input:
    path cross_sequence_summary_with_interface_metrics_csv
    // Per-survivor summary CSVs have fixed filenames
    // (af3_nomsa_summary.csv, biophysical_summary.csv, rosetta_summary.csv)
    // because each upstream process emits one file with that name.  When
    // the cohort has ≥2 survivors, staging them all into the same
    // directory causes Nextflow to throw an "input file name collision"
    // error.  The 'inputs/?/*' pattern tells Nextflow to put each file
    // in its own auto-numbered subdirectory (inputs/1/file.csv,
    // inputs/2/file.csv, ...) so they don't clobber each other.  The
    // merge script joins by the seq_name column inside each CSV, not
    // by filename, so the disambiguation has no semantic effect.
    path af3_summaries,         stageAs: 'af3_inputs/?/*'
    path biophysical_summaries, stageAs: 'biophysical_inputs/?/*'
    path rosetta_summaries,     stageAs: 'rosetta_inputs/?/*'
    path merge_script

    output:
    path "survivors_with_orthogonal_metrics.csv",
         emit: final_csv

    script:
    """
    set -euo pipefail

    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python ${merge_script} \\
            --input-csv ${cross_sequence_summary_with_interface_metrics_csv} \\
            --af3-summaries-glob 'af3_inputs/*/*.csv' \\
            --biophysical-summaries-glob 'biophysical_inputs/*/*.csv' \\
            --rosetta-summaries-glob 'rosetta_inputs/*/*.csv' \\
            --filter-af3-ra-max ${params.orthogonal_filter_af3_ra_max} \\
            --filter-sc-min ${params.orthogonal_filter_sc_min} \\
            --filter-bsa-min ${params.orthogonal_filter_bsa_min} \\
            --filter-plddt-min ${params.orthogonal_filter_plddt_min} \\
            --output-csv survivors_with_orthogonal_metrics.csv

    echo "Orthogonal metrics merged."
    head -1 survivors_with_orthogonal_metrics.csv
    echo "(row count: \$(tail -n +2 survivors_with_orthogonal_metrics.csv | wc -l))"
    """
}


/*
 * ORTHOG_PLOTS
 * ────────────
 * Diagnostic plot suite for orthogonal metrics.  Companion to
 * NEGSTEER_PLOTS in modules/negative_steering.nf — that one plots
 * negsteer cohort summaries; this one plots the orthogonal-metrics
 * cascade and joins them back to the upstream cohort context.
 *
 * Mirrors tests/orthogonal_metrics/test_orthogonal_metrics_plots.py
 * exactly — bin/orthogonal_metrics_plots.py is a verbatim lift with
 * only the production docstring header changed.
 *
 * Four PNGs:
 *   1. orthogonal_af3_vs_boltz.png             cross-model ra_eff agreement
 *   2. orthogonal_filter_cascade.png           gating cascade waterfall
 *                                              with trailing AF3 informational bar
 *   3. orthogonal_metrics_vs_composite.png     2x3 scatter: composite vs each metric
 *   4. orthogonal_combined_cohort_summary.png  the big combined view (every
 *                                              upstream cross_summary sequence,
 *                                              negsteer + orthogonal columns)
 *
 * Inputs are both CSVs (no per-sequence workdirs needed — orthogonal
 * plot inputs are aggregated tables, no glob discovery required).
 *
 * Container: rfdiff_container — same as NEGSTEER_PLOTS; only depends
 * on stdlib csv + numpy + matplotlib, no GPU.
 */
process ORTHOG_PLOTS {
    tag "orthog_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path survivors_csv
    path cross_summary_csv
    // Stage the plot script as a path input so Nextflow's cache key
    // includes its CONTENT.  Without this, edits to
    // bin/orthogonal_metrics_plots.py do not invalidate the cache
    // (the script command line below is byte-identical regardless of
    // the file's contents at ${projectDir}/bin/), and `-resume` will
    // silently skip re-rendering with the new script.  Same class of
    // gotcha as AF3_PARSE_OUTPUT (negsteer_af3_nomsa.nf).  Audit the
    // rest of the pipeline for the same pattern as part of Task 56.
    path plot_script

    output:
    path "orthogonal_*.png", emit: plots

    script:
    """
    set -euo pipefail

    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rfdiff_container} \\
        python ${plot_script} \\
            --survivors-csv     ${survivors_csv} \\
            --cross-summary-csv ${cross_summary_csv} \\
            --outdir            .
    """
}
