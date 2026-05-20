/*
 * =============================================================================
 * Pose Solver module
 * =============================================================================
 * Constraint-driven rigid-body docking via bin/pose_solver.py.  Returns
 * one deterministic pose (no clusters to rank), driven by 2-4
 * user-supplied CA-CA pair constraints, in minutes.  See
 * notes/pipeline_notes/pipeline_notes16.md for the geometric rationale
 * behind the defaults.
 *
 * Process chain (Branch A: separate receptor + effector PDBs):
 *
 *   POSE_SOLVER_PREPARE     — normalise input chain IDs
 *           │
 *           ▼
 *   POSE_SOLVE              — run pose_solver.py → solved_pose_posed.pdb
 *           │
 *           ├─── POSE_INTERFACE_METRICS   (BSA, clashes, gap-index, H-bonds)
 *           │
 *           ├─── POSE_SOLVER_PLOTS        (loss breakdown, pair distances,
 *           │                              interface dashboard, restart-loss
 *           │                              curve, contig comparison)
 *           │
 *           ▼
 *   BUILD_CONTIGS           — derive updated contig string for RFDiffusion
 *
 * Chain convention: the user's params.receptor_chain / params.effector_chain
 * are authoritative end-to-end.  POSE_SOLVER_PREPARE rewrites both input
 * PDBs to use those chain IDs; pose_solver.py preserves them in the
 * solved complex; BUILD_CONTIGS reads them unchanged.
 */


/*
 * POSE_SOLVER_PREPARE
 * -------------------
 * Normalise the two input PDBs to use the user-declared chain IDs.
 * Guards against the two-monomers-both-named-A case that would
 * otherwise produce an output complex with duplicate chain IDs.
 */
process POSE_SOLVER_PREPARE {
    tag "pose_prep"
    label 'cpu'

    publishDir "${params.outdir}/pose_solver/data", mode: 'copy'

    input:
    path receptor_pdb
    path effector_pdb
    val  receptor_chain
    val  effector_chain
    path prepare_script

    output:
    path "receptor_pose.pdb", emit: receptor_pdb_out
    path "effector_pose.pdb", emit: effector_pdb_out

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${prepare_script} \\
            --receptor ${receptor_pdb} \\
            --effector ${effector_pdb} \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain}
    """
}


/*
 * POSE_SOLVE
 * ----------
 * Run bin/pose_solver.py.  Outputs:
 *   solved_pose_posed.pdb     - docked complex (receptor + effector,
 *                               original chain IDs preserved)
 *   solved_pose_results.json  - pair distances, loss, clashes, contigs
 *   solved_pose_heatmap.png   - contact / clash heatmap
 *
 * Per-restart loss history is dumped to solved_pose_restart_losses.csv
 * for the POSE_SOLVER_PLOTS restart-curve plot.
 */
process POSE_SOLVE {
    tag "pose_solve"
    label 'cpu'

    publishDir "${params.outdir}/pose_solver", mode: 'copy'

    input:
    path receptor_pdb
    path effector_pdb
    val  receptor_chain
    val  effector_chain
    val  pairs
    val  exclusions
    val  min_pair_distance
    val  max_pair_distance
    val  pair_sc_clash_cutoff
    val  clash_cutoff
    val  contact_cutoff
    val  n_restarts
    val  use_de
    val  global_interp
    val  interp_weight
    val  contig_design_region
    path solver_script

    output:
    path "solved_pose_posed.pdb",            emit: posed_pdb
    path "solved_pose_results.json",         emit: results_json
    path "solved_pose_heatmap.png",          emit: heatmap
    path "solved_pose_restart_losses.csv",   emit: restart_losses

    script:
    def excl_arg     = exclusions          ? "--exclusions '${exclusions}'" : ""
    def de_arg       = use_de              ? "--use-de" : ""
    def global_arg   = global_interp       ? "--global-interp" : ""
    def contig_arg   = contig_design_region
                       ? "--contig-design-region '${contig_design_region}'" : ""
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${solver_script} \\
            --binder ${receptor_pdb} \\
            --target ${effector_pdb} \\
            --binder-chain ${receptor_chain} \\
            --target-chain ${effector_chain} \\
            --pairs '${pairs}' \\
            ${excl_arg} \\
            --min-pair-distance ${min_pair_distance} \\
            --max-pair-distance ${max_pair_distance} \\
            --pair-sc-clash-cutoff ${pair_sc_clash_cutoff} \\
            --clash-cutoff ${clash_cutoff} \\
            --contact-cutoff ${contact_cutoff} \\
            --n-restarts ${n_restarts} \\
            --interp-weight ${interp_weight} \\
            ${global_arg} \\
            ${de_arg} \\
            ${contig_arg} \\
            --output-prefix solved_pose
    """
}


/*
 * POSE_INTERFACE_METRICS
 * ----------------------
 * Run bin/pose_interface_metrics.py on the solved complex.  Computes
 * BSA (Shrake-Rupley), gap-index, CA / heavy contact counts, clash
 * counts at 2/2.5/3/3.5 Å, interface residue counts + composition,
 * and H-bond donor->acceptor candidates.  Output is a single JSON
 * consumed by POSE_SOLVER_PLOTS.
 *
 * Runs in boltz2_container (provides biopython).
 */
process POSE_INTERFACE_METRICS {
    tag "pose_metrics"
    label 'cpu'

    publishDir "${params.outdir}/pose_solver", mode: 'copy'

    input:
    path posed_pdb
    val  receptor_chain
    val  effector_chain
    path metrics_script

    output:
    path "interface_metrics.json", emit: metrics_json

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python ${metrics_script} \\
            ${posed_pdb} \\
            --chains ${receptor_chain}:${effector_chain} \\
            --json-out interface_metrics.json
    """
}


/*
 * POSE_SOLVER_PLOTS
 * -----------------
 * Render the five diagnostic plots:
 *   1. pose_loss_breakdown.png       - per-term loss contributions
 *   2. pose_pair_distances.png       - pair-constraint distance whiskers
 *   3. pose_interface_dashboard.png  - 4-panel interface metrics
 *   4. pose_restart_loss_curve.png   - loss vs restart index
 *   5. pose_contig_comparison.png    - three contig options stacked
 *
 * The heatmap is already produced by POSE_SOLVE; not regenerated here.
 */
process POSE_SOLVER_PLOTS {
    tag "pose_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path results_json
    path metrics_json
    path restart_losses
    path posed_pdb
    val  receptor_chain
    val  effector_chain
    path plots_script

    output:
    path "pose_*.png", emit: plots

    script:
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rfdiff_container} \\
        python ${plots_script} \\
            --results-json ${results_json} \\
            --metrics-json ${metrics_json} \\
            --restart-losses ${restart_losses} \\
            --posed-pdb ${posed_pdb} \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain}
    """
}


/*
 * BUILD_CONTIGS
 * -------------
 * Given the user's original contig string and the solved complex PDB:
 *   - Detects renumbering + pLDDT trimming offsets
 *   - Replaces bare effector chain references with actual PDB length
 *   - Adjusts receptor residue numbers for any renumbering offset
 *   - Extracts receptor/effector sequences from the complex
 *
 * Note: rec_trim_mapping / eff_trim_mapping inputs are currently
 * always empty placeholder JSONs ({}) produced by WRITE_DUMMY_MAPPING.
 * They existed originally to feed pLDDT-trimming offsets from AF2
 * monomer (since removed) and are kept in the API surface so the same
 * plumbing can be reused if disordered-region trimming is reinstated
 * later.
 */
process BUILD_CONTIGS {
    tag "build_contigs"
    label 'cpu'

    publishDir "${params.outdir}/pose_solver", mode: 'copy'

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
