/*
 * =============================================================================
 * HADDOCK3 module
 * =============================================================================
 * Docking, cluster qualification, per-cluster metrics, BSA + pair-contact
 * driven cluster selection (with optional manual override).  Restructured
 * in Session 7 per notes/design_audit.md Q129-Q157.
 *
 * Restraints come from user params, NOT from the contig string:
 *   - haddock_contact_pairs           — hard pin (unambig_restraints.tbl)
 *   - haddock_receptor_active_residues / haddock_effector_active_residues
 *     — soft preference (ambig_restraints.tbl, AND-of-OR with 50% nrest)
 *
 * Auto-pick at cluster selection uses (-pair_contact_fraction, -bsa)
 * lexicographic; user can override with params.haddock_chosen_cluster.
 *
 * EXTRACT_HOTSPOTS is gone (deleted in commit 3); hotspots for RFDiffusion
 * come ONLY from params.hotspot (user biological knowledge).
 */


/*
 * HADDOCK3_PREPARE
 * ----------------
 * Relabel input PDBs to chain A (receptor) / B (effector) and write the
 * AIR + unambig restraint files from the user's params.
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
    val  contact_pairs
    val  receptor_active_residues
    val  effector_active_residues
    val  pair_distance
    val  contigs
    path prepare_script

    output:
    path "receptor_haddock.pdb",     emit: receptor_pdb_out
    path "effector_haddock.pdb",     emit: effector_pdb_out
    path "ambig_restraints.tbl",     emit: ambig_restraints
    path "unambig_restraints.tbl",   emit: unambig_restraints
    path "restraints_summary.json",  emit: restraints_summary

    script:
    def pairs_arg = contact_pairs        ? "--contact-pairs '${contact_pairs}'" : ""
    def rec_arg   = receptor_active_residues
                                          ? "--receptor-active-residues '${receptor_active_residues}'" : ""
    def eff_arg   = effector_active_residues
                                          ? "--effector-active-residues '${effector_active_residues}'" : ""
    def contigs_arg = contigs            ? "--contigs '${contigs}'" : ""
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${prepare_script} \\
            --receptor ${receptor_pdb} \\
            --effector ${effector_pdb} \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain} \\
            ${pairs_arg} \\
            ${rec_arg} \\
            ${eff_arg} \\
            --pair-distance '${pair_distance}' \\
            ${contigs_arg}
    """
}


/*
 * HADDOCK3_DOCK
 * -------------
 * Run HADDOCK3 docking with the two restraint files; post-process via
 * collect_haddock3_dock.py to extract one best-model PDB per qualifying
 * cluster and write haddock_report.json.
 *
 * The CNS engine reads ambig_fname for AIRs and unambig_fname for hard
 * pin pairs; either file may be empty.
 */
process HADDOCK3_DOCK {
    tag "haddock3"
    label 'cpu_haddock'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path receptor_pdb
    path effector_pdb
    path ambig_restraints
    path unambig_restraints
    val  haddock_sampling
    val  haddock_seletop
    path collect_script

    output:
    path "run/run-haddock/",         emit: run_dir
    path "best_cluster*.pdb",        emit: cluster_models
    path "capri_scores.tsv",         emit: capri_scores,      optional: true
    path "cluster_summary.txt",      emit: cluster_summary,   optional: true
    path "haddock_report.json",      emit: haddock_report

    script:
    """
    mkdir -p data run

    cp ${receptor_pdb}        data/receptor.pdb
    cp ${effector_pdb}        data/effector.pdb
    cp ${ambig_restraints}    data/ambig_restraints.tbl
    cp ${unambig_restraints}  data/unambig_restraints.tbl

    # ── HADDOCK3 config ──────────────────────────────────────────────────
    # ambig_fname is always wired; the AIRs file may be empty (no-op).
    # unambig_fname is wired only when there are contact pairs.
    HAS_UNAMBIG=\$([ -s data/unambig_restraints.tbl ] && echo "yes" || echo "no")

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
\$( [ "\$HAS_UNAMBIG" = "yes" ] && echo 'unambig_fname = "../data/unambig_restraints.tbl"' )
sampling = ${haddock_sampling}
concat = 20

[seletop]
select = ${haddock_seletop}

[flexref]
ambig_fname = "../data/ambig_restraints.tbl"
\$( [ "\$HAS_UNAMBIG" = "yes" ] && echo 'unambig_fname = "../data/unambig_restraints.tbl"' )
concat = 20

[emref]
ambig_fname = "../data/ambig_restraints.tbl"
\$( [ "\$HAS_UNAMBIG" = "yes" ] && echo 'unambig_fname = "../data/unambig_restraints.tbl"' )
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

    # ── Collect cluster best-models ──────────────────────────────────────
    # (ambig_restraints.tbl and unambig_restraints.tbl are already staged
    # into the workdir by Nextflow as `path` inputs; they flow into
    # HADDOCK_CLUSTER_METRICS via the HADDOCK3_PREPARE output channels.)
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${collect_script} \\
            --run-dir run/run-haddock \\
            --min-cluster-size ${params.haddock_min_cluster_size}
    """
}


/*
 * HADDOCK_CLUSTER_METRICS
 * -----------------------
 * Compute BSA / hbonds / COM / AIR satisfaction / pair contact fraction /
 * clash counts (in/out design region) for every qualifying cluster.
 * Runs in boltz2_container (needs numpy + freesasa + MDAnalysis).
 *
 * Sc is NOT computed here — see HADDOCK_CLUSTER_SC sibling process which
 * runs in rosetta_container and emits cluster_sc.json.  HaddockRun
 * merges the two files at construction time.
 */
process HADDOCK_CLUSTER_METRICS {
    tag "haddock_metrics"
    label 'cpu_haddock_metrics'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path haddock_report
    path cluster_models
    path restraints_summary
    path ambig_restraints
    path unambig_restraints
    path metrics_script
    path bin_dir

    output:
    path "cluster_metrics.json", emit: cluster_metrics
    path "best_cluster*.pdb",    emit: cluster_models_pass

    script:
    """
    # haddock_cluster_metrics expects its scripts side-by-side in bin/
    # (it shells out to run_biophysical_metrics.py).  Symlink everything
    # into a single workdir so the relative paths line up.
    for f in ${bin_dir}/*.py ${bin_dir}/*.xml; do
        ln -sf "\${f}" .
    done

    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python haddock_cluster_metrics.py \\
            --workdir . \\
            --receptor-chain A \\
            --effector-chain B
    """
}


/*
 * HADDOCK_CLUSTER_SC
 * ------------------
 * Compute Rosetta shape complementarity (Sc) for every qualifying
 * cluster's best model.  Runs in rosetta_container; output is merged
 * onto HaddockCluster fields downstream via cluster_sc.json.
 */
process HADDOCK_CLUSTER_SC {
    tag "haddock_sc"
    label 'cpu_haddock_metrics'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path haddock_report
    path cluster_models
    path sc_script
    path bin_dir

    output:
    path "cluster_sc.json", emit: cluster_sc

    script:
    """
    for f in ${bin_dir}/*.py ${bin_dir}/*.xml; do
        ln -sf "\${f}" .
    done

    singularity exec --bind \${PWD}:\${PWD} ${params.rosetta_container} \\
        python haddock_cluster_sc.py \\
            --workdir . \\
            --receptor-chain A \\
            --effector-chain B
    """
}


/*
 * SELECT_HADDOCK_CLUSTER
 * ----------------------
 * Pick the cluster downstream stages consume — auto-pick by
 * (-pair_contact_fraction, -bsa) or honour params.haddock_chosen_cluster
 * when set.  Emits selected_complex.pdb + selected_cluster_id.txt + a
 * sorted cluster_metrics_table.csv for inspection.
 */
process SELECT_HADDOCK_CLUSTER {
    tag "haddock_select"
    label 'cpu'

    publishDir "${params.outdir}/haddock", mode: 'copy'

    input:
    path haddock_report
    path cluster_metrics
    path cluster_sc
    path restraints_summary
    path cluster_models
    path receptor_pdb
    path effector_pdb
    val  contact_pairs
    val  receptor_active_residues
    val  effector_active_residues
    val  pair_distance
    val  chosen_cluster_id
    path select_script

    output:
    path "selected_complex.pdb",         emit: selected_pdb
    path "selected_cluster_id.txt",      emit: selected_id
    path "cluster_metrics_table.csv",    emit: cluster_table

    script:
    def chosen_arg = chosen_cluster_id != null && "${chosen_cluster_id}" != "null"
                     ? "--chosen-cluster-id ${chosen_cluster_id}" : ""
    def pairs_arg = contact_pairs ? "--contact-pairs '${contact_pairs}'" : ""
    def rec_arg   = receptor_active_residues
                    ? "--receptor-active-residues '${receptor_active_residues}'" : ""
    def eff_arg   = effector_active_residues
                    ? "--effector-active-residues '${effector_active_residues}'" : ""
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${select_script} \\
            --workdir . \\
            --receptor-pdb ${receptor_pdb} \\
            --effector-pdb ${effector_pdb} \\
            ${pairs_arg} \\
            ${rec_arg} \\
            ${eff_arg} \\
            --pair-distance '${pair_distance}' \\
            ${chosen_arg}
    """
}


/*
 * HADDOCK3_PLOTS
 * ──────────────
 * Score vs BSA scatter, cluster size bar chart, and per-cluster interface
 * contact heatmaps (receptor + effector).
 *
 * The active-residue inputs are renamed (haddock_ prefix) and the contig
 * input is gone — plotting code no longer derives active residues from
 * the contig string.
 */
process HADDOCK3_PLOTS {
    tag "haddock_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path capri_scores
    path cluster_summary
    path run_dir
    path restraints_summary
    path cluster_metrics
    path cluster_sc
    val  contact_pairs
    val  effector_active_residues
    path plots_script

    output:
    path "haddock_*.png", emit: plots

    script:
    def pairs_arg = contact_pairs
                    ? "--contact-pairs '${contact_pairs}'" : ""
    def eff_arg = effector_active_residues
                  ? "--effector-active-residues '${effector_active_residues}'" : ""
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rfdiff_container} \\
        python ${plots_script} \\
            --capri-scores ${capri_scores} \\
            --cluster-summary ${cluster_summary} \\
            --run-dir ${run_dir} \\
            --restraints-summary ${restraints_summary} \\
            --cluster-metrics ${cluster_metrics} \\
            --cluster-sc ${cluster_sc} \\
            --min-cluster-size ${params.haddock_min_cluster_size} \\
            ${pairs_arg} \\
            ${eff_arg}
    """
}


/*
 * BUILD_CONTIGS
 * ─────────────
 * Given the user's original contig string, the trim mappings, and the
 * SELECT-chosen docked complex, auto-update the contig string so that:
 *   - Effector length is derived from the PDB (no manual "B1-123")
 *   - Receptor residue numbers are adjusted for any renumbering offset
 *     between the contig string and the docked PDB
 *   - The contig uses correct numbering for the docked PDB
 * Also extracts final receptor/effector sequences from the complex.
 *
 * Now operates on selected_complex.pdb (chosen by SELECT_HADDOCK_CLUSTER)
 * instead of HADDOCK_DOCK's old single best_model.pdb.
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
