/*
 * =============================================================================
 * Rosetta pre-validation physics filtering module
 * =============================================================================
 * Runs Rosetta InterfaceAnalyzer on RFDiffusion backbone designs to compute
 * interface shape complementarity (Sc), binding energy (dG_separated), and
 * buried surface area (dSASA_int).  Designs below the Sc threshold are
 * filtered out before ProteinMPNN sequence design.
 *
 * Processes:
 *   ROSETTA_SC              Per-design InterfaceAnalyzer (parallelised)
 *   ROSETTA_FILTER          Collect scores, apply Sc threshold, emit passing PDBs
 *   ROSETTA_FILTER_PLOTS    Diagnostic plots (Sc distribution, dG distribution)
 *
 * Container notes:
 *   Rosetta is installed inside a separate Singularity container (Rosetta.img).
 *   The binary suffix is auto-detected at runtime because it differs between
 *   source and binary builds (e.g. .static.linuxgccrelease vs
 *   .default.linuxgccrelease).
 *
 *   --bind $PWD:$PWD is required for writable CWD access on NBI HPC.
 *
 * Input PDBs:
 *   Split two-chain PDBs from RFDIFFUSION_FILTER (chain A = receptor,
 *   chain B = effector, renumbered from 1).  These are polyvaline backbones,
 *   so no FastRelax or score_jd2 pre-scoring is needed — InterfaceAnalyzer
 *   runs directly on the backbone structures.
 * =============================================================================
 */


/*
 * ROSETTA_SC
 * ----------
 * Run InterfaceAnalyzer on a single RFDiffusion split PDB to compute
 * interface metrics: sc_value, dG_separated, dSASA_int, packstat, etc.
 *
 * One process per design PDB — Nextflow schedules these as parallel SLURM
 * jobs automatically.
 *
 * The binary suffix is auto-detected inside the container so the process
 * works with both source and binary Rosetta builds.
 */
process ROSETTA_SC {
    tag "${design_pdb.baseName}"
    label 'cpu'

    input:
    path design_pdb

    output:
    tuple path(design_pdb), path("interface_scores_${design_pdb.baseName}.sc"), emit: pdb_and_scores

    script:
    """
    # ── Detect Rosetta binary suffix ──────────────────────────────────────
    ROSETTA_BIN=/opt/rosetta/main/source/bin
    ROSETTA_DB=/opt/rosetta/main/database

    SUFFIX=\$(singularity exec ${params.rosetta_container} bash -c \\
        "ls \${ROSETTA_BIN}/InterfaceAnalyzer.*.linuxgcc*release 2>/dev/null | head -1 | sed 's|.*/InterfaceAnalyzer||'" \\
    )

    if [ -z "\$SUFFIX" ]; then
        echo "ERROR: Could not detect InterfaceAnalyzer binary suffix"
        exit 1
    fi
    echo "Detected binary suffix: \$SUFFIX"

    # ── Run InterfaceAnalyzer ─────────────────────────────────────────────
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        ${params.rosetta_container} \\
        \${ROSETTA_BIN}/InterfaceAnalyzer\${SUFFIX} \\
            -database \$ROSETTA_DB \\
            -s ${design_pdb} \\
            -no_optH false \\
            -ignore_unrecognized_res \\
            -pack_separated \\
            -pack_input \\
            -compute_interface_sc true \\
            -add_regular_scores_to_scorefile \\
            -use_input_sc \\
            -out:file:score_only interface_scores_${design_pdb.baseName}.sc \\
            -out:no_nstruct_label

    echo "InterfaceAnalyzer completed for ${design_pdb.baseName}"
    """
}


/*
 * ROSETTA_FILTER
 * ──────────────
 * Collect InterfaceAnalyzer score files from all designs, apply the Sc
 * threshold, and emit:
 *   - rosetta_filter_metrics.json  (all designs, for plotting)
 *   - rosetta_passing_designs.txt  (filenames of designs that pass)
 *   - rosetta_filter_summary.json  (counts)
 *   - passing/*.pdb                (PDB files for downstream ProteinMPNN)
 *
 * Designs with sc_value >= sc_threshold pass.  dG_separated and dSASA_int
 * are recorded for all designs but only used for reference/plotting.
 */
process ROSETTA_FILTER {
    tag "rosetta_filter"
    label 'cpu'

    publishDir "${params.outdir}/rosetta_filtering", mode: 'copy',
        pattern: '{rosetta_filter_metrics.json,rosetta_filter_summary.json,rosetta_passing_designs.txt}'
    publishDir "${params.outdir}/rosetta_filtering", mode: 'copy',
        pattern: 'passing/*.pdb'

    input:
    path pdb_and_score_files   // staged into work dir, discovered by glob in collect script
    val  sc_threshold

    output:
    path "rosetta_filter_metrics.json",    emit: metrics
    path "rosetta_filter_summary.json",    emit: filter_summary
    path "rosetta_passing_designs.txt",    emit: passing_list
    path "passing/*.pdb",                  emit: passing_pdbs, optional: true

    script:
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        ${params.rosetta_container} \\
        python ${projectDir}/bin/rosetta_filter_collect.py \\
            --score-dir . \\
            --pdb-dir . \\
            --sc-threshold ${sc_threshold}
    """
}


/*
 * ROSETTA_FILTER_PLOTS
 * ────────────────────
 * Two diagnostic plots from the Rosetta pre-validation filter:
 *   1. Shape complementarity (Sc) distribution with threshold line
 *   2. dG_separated distribution (pre-MPNN backbone binding energy landscape)
 */
process ROSETTA_FILTER_PLOTS {
    tag "rosetta_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path metrics

    output:
    path "rosetta_*.png", emit: plots

    script:
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rosetta_container} \\
        python ${projectDir}/bin/rosetta_filter_plots.py \\
            --metrics ${metrics}
    """
}
