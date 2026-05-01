/*
 * =============================================================================
 * HADDOCK3 module
 * =============================================================================
 * Docking, cluster validation, best model selection, interface analysis,
 * hotspot extraction, and diagnostic plots.
 *
 * Restraints are auto-generated from the de novo regions in the RFDiffusion
 * contig string.
 *
 * All Python logic lives in bin/ scripts. Scripts are called via
 * singularity exec so that all pipeline dependencies (including matplotlib)
 * are available consistently across all processes.
 * =============================================================================
 */


/*
 * HADDOCK3_PREPARE
 * ----------------
 * Merge receptor and effector PDBs into a pair for HADDOCK.
 * Auto-generate ambiguous restraints from the contig string.
 * Parse contigs to find the de novo gaps on the receptor and use
 * them as active residues.
 */
process HADDOCK3_PREPARE {
    tag "haddock_prep"
    label 'cpu'

    publishDir "${params.outdir}/haddock/data", mode: 'copy'

    input:
    path receptor_pdb
    path effector_pdb
    val  receptor_chain
    val  effector_chain
    val  contigs
    val  effector_active_residues
    path prepare_script

    output:
    path "receptor_haddock.pdb",    emit: receptor_pdb_out
    path "effector_haddock.pdb",    emit: effector_pdb_out
    path "ambig_restraints.tbl",    emit: restraints

    script:
    def eff_res_arg = effector_active_residues ? "--effector-active-residues '${effector_active_residues}'" : ""
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${prepare_script} \\
            --receptor ${receptor_pdb} \\
            --effector ${effector_pdb} \\
            --contigs "${contigs}" \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain} \\
            ${eff_res_arg}
    """
}


/*
 * HADDOCK3_DOCK
 * -------------
 * Run HADDOCK3 docking.
 * collect_haddock3_dock.py post-processes the run directory to extract the
 * best model, CAPRI scores, and cluster summary.
 */
process HADDOCK3_DOCK {
    tag "haddock3"
    label 'cpu_haddock'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path receptor_pdb
    path effector_pdb
    path restraints
    val  haddock_sampling
    val  haddock_seletop
    path collect_script

    output:
    path "run/run-haddock/",      emit: run_dir
    path "best_model.pdb",        emit: best_model
    path "best_cluster*.pdb",     emit: cluster_models,    optional: true
    path "capri_scores.tsv",      emit: capri_scores,      optional: true
    path "cluster_summary.txt",   emit: cluster_summary,   optional: true
    path "haddock_report.json",   emit: haddock_report

    script:
    """
    mkdir -p data run

    cp ${receptor_pdb} data/receptor.pdb
    cp ${effector_pdb} data/effector.pdb
    cp ${restraints}   data/ambig_restraints.tbl

    # ── HADDOCK3 config ──────────────────────────────────────────────────
    cat > run/docking.cfg << HADDOCK_CFG
run_dir = "run-haddock"
ncores = ${task.cpus}
mode = "local"

molecules = [
    "../data/receptor.pdb",
    "../data/effector.pdb"
]

[topoaa]

[rigidbody]
ambig_fname = "../data/ambig_restraints.tbl"
sampling = ${haddock_sampling}
concat = 20

[seletop]
select = ${haddock_seletop}

[flexref]
ambig_fname = "../data/ambig_restraints.tbl"
concat = 20

[emref]
ambig_fname = "../data/ambig_restraints.tbl"
concat = 20

[clustfcc]
clust_cutoff = 0.6
min_population = ${params.haddock_min_cluster_size}
# Single source of truth: params.haddock_min_cluster_size in nextflow.config.
# The same value is plumbed via --min-cluster-size to bin/collect_haddock3_dock.py
# and bin/haddock3_plots.py below so all three sites stay in lockstep.

[seletopclusts]
top_models = 4

[caprieval]
HADDOCK_CFG

    # ── Run HADDOCK3 ─────────────────────────────────────────────────────
    # --pid: run in own PID namespace so Singularity init reaps zombies
    WORKDIR=\${PWD}
    cd run
    singularity exec --pid --bind \${WORKDIR}:\${WORKDIR} ${params.rfdiff_container} \\
        haddock3 docking.cfg
    cd ..

    # ── Collect results ───────────────────────────────────────────────────
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${collect_script} \\
            --run-dir run/run-haddock \\
            --min-cluster-size ${params.haddock_min_cluster_size}
    """
}


/*
 * HADDOCK3_PLOTS
 * ──────────────
 * Score vs BSA scatter, cluster size bar chart, and per-cluster interface
 * contact heatmaps (receptor + effector).
 * cluster_summary.txt is read directly from clustfcc.tsv so cluster sizes
 * reflect true membership rather than just caprieval representatives.
 */
process HADDOCK3_PLOTS {
    tag "haddock_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path capri_scores
    path cluster_summary
    path run_dir
    val  contigs
    val  receptor_chain
    val  effector_active_residues
    path plots_script

    output:
    path "haddock_*.png", emit: plots

    script:
    def eff_res_arg = effector_active_residues ? "--effector-active-residues '${effector_active_residues}'" : ""
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rfdiff_container} \\
        python ${plots_script} \\
            --capri-scores ${capri_scores} \\
            --cluster-summary ${cluster_summary} \\
            --run-dir ${run_dir} \\
            --contigs '${contigs}' \\
            --receptor-chain ${receptor_chain} \\
            --min-cluster-size ${params.haddock_min_cluster_size} \\
            ${eff_res_arg}
    """
}


/*
 * EXTRACT_HOTSPOTS
 * ----------------
 * Analyse the HADDOCK best model to identify effector residues at the
 * interface (within contact distance of receptor).
 * Outputs the hotspot string in RFDiffusion format: "B24,B25,..."
 */
process EXTRACT_HOTSPOTS {
    tag "extract_hotspots"
    label 'cpu'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path best_model
    val  receptor_chain
    val  effector_chain
    val  contact_cutoff
    val  receptor_seq
    val  effector_seq
    path extract_script

    output:
    path "hotspot_string.txt",      emit: hotspot_string
    path "interface_analysis.json", emit: interface_json

    script:
    def rec_seq_arg = receptor_seq ? "--receptor-seq '${receptor_seq}'" : ""
    def eff_seq_arg = effector_seq ? "--effector-seq '${effector_seq}'" : ""
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${extract_script} \\
            --complex ${best_model} \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain} \\
            --cutoff ${contact_cutoff} \\
            ${rec_seq_arg} \\
            ${eff_seq_arg}
    """
}


/*
 * BUILD_CONTIGS
 * ─────────────
 * Given the user's original contig string, the trim mappings, and the
 * docked complex, auto-update the contig string so that:
 *   - Effector length is derived from the PDB (no manual "B1-123")
 *   - Receptor residue numbers are adjusted for any renumbering offset
 *     between the contig string and the docked PDB
 *   - The contig uses correct numbering for the docked PDB
 * Also extracts final receptor/effector sequences from the complex.
 *
 * Note: rec_trim_mapping / eff_trim_mapping inputs are currently always
 * empty placeholder JSONs ({}) produced by WRITE_DUMMY_MAPPING.  They
 * existed originally to feed pLDDT-trimming offsets from AF2 monomer
 * (since removed) and are kept in the API surface so the same plumbing
 * can be reused if disordered-region trimming is reinstated later.
 */
process BUILD_CONTIGS {
    tag "build_contigs"
    label 'cpu'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path complex_pdb
    val  receptor_chain
    val  effector_chain
    val  user_contigs
    path rec_trim_mapping
    path eff_trim_mapping
    path build_script

    output:
    path "updated_contigs.txt",    emit: contigs_txt
    path "updated_params.json",    emit: updated_params
    path "rfdiffusion_input.pdb",  emit: rfdiffusion_pdb

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${build_script} \\
            --complex ${complex_pdb} \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain} \\
            --contigs "${user_contigs}" \\
            --rec-trim-mapping ${rec_trim_mapping} \\
            --eff-trim-mapping ${eff_trim_mapping}
    """
}
