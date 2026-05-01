/*
 * =============================================================================
 * RFDiffusion module
 * =============================================================================
 * Runs backbone generation, filters designs by interface quality, and
 * produces diagnostic plots.
 *
 * All Python logic lives in bin/ scripts. Scripts are called via
 * singularity exec so that all pipeline dependencies (including matplotlib)
 * are available consistently across all processes.
 *
 * Container notes (Singularity 3.8 + .img format on NBI HPC):
 *   - --writable-tmpfs does NOT work (Permission denied)
 *   - --no-home breaks CWD access under /hpc-home
 *   - --env HOME is blocked by admin policy
 *   - --bind $PWD:$PWD is REQUIRED to make work dir writable
 *   - --bind local:/container/path works for specific writable paths
 *   - MPLCONFIGDIR=/tmp for matplotlib cache (always writable)
 */


/*
 * RFDIFFUSION
 * -----------
 * Generate backbone designs using RFDiffusion.
 *
 * The user contig string is preprocessed by rfdiffusion_contigs.py to
 * convert user-friendly notation (e.g. bare chain letters, single-residue
 * shorthands) into the format expected by RFDiffusion's native parser.
 */
process RFDIFFUSION {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/rfdiffusion", mode: 'copy'

    input:
    path pdb_file
    val  contigs
    val  hotspot
    val  num_designs
    val  iterations
    // Stage bin scripts as path inputs so Nextflow content-hashes them
    // for the task-cache key.  Without this, edits to bin/*.py do not
    // invalidate the cache (see comment on AF3_PARSE_OUTPUT in
    // modules/negsteer_af3_nomsa.nf).
    path contigs_script

    output:
    path "design_*.pdb",  emit: design_pdbs
    path "traj/",         emit: trajectories, optional: true

    script:
    def raw_contigs = contigs.replaceAll("'", "")
    def hotspot_arg = hotspot ? "\"ppi.hotspot_res=[${hotspot}]\"" : ""
    """
    cp ${pdb_file} input.pdb
    mkdir -p traj

    # ─── Preprocess contigs ───────────────────────────────────────────────
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${contigs_script} \\
            --contigs "${raw_contigs}" \\
            --pdb input.pdb \\
            --output processed_contigs.txt

    PROCESSED_CONTIGS=\$(cat processed_contigs.txt)

    echo "Raw contigs:       ${raw_contigs}"
    echo "Processed contigs: \${PROCESSED_CONTIGS}"

    # ─── Schedules cache ──────────────────────────────────────────────────
    mkdir -p schedules_cache
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        bash -c "cp /opt/RFdiffusion/schedules/*.pkl \${PWD}/schedules_cache/ 2>/dev/null; true"

    # ─── Run RFDiffusion ──────────────────────────────────────────────────
    singularity exec --nv \\
        --bind \${PWD}:\${PWD} \\
        --bind \${PWD}/schedules_cache:/opt/RFdiffusion/schedules \\
        ${params.rfdiff_container} \\
        python /opt/RFdiffusion/run_inference.py \\
            inference.output_prefix="./design" \\
            inference.input_pdb="input.pdb" \\
            inference.num_designs=${num_designs} \\
            "contigmap.contigs=[\${PROCESSED_CONTIGS}]" \\
            ${hotspot_arg} \\
            diffuser.T=${iterations}
    """
}


/*
 * RFDIFFUSION_FILTER
 * ──────────────────
 * Analyse each RFDiffusion design and filter by the fraction of
 * effector (hotspot) contacts that land inside the design region
 * vs outside it.
 *
 * Also computes design-region RMSD and per-residue contact maps for
 * downstream diagnostic plots.
 *
 * Designs that fail the filter are still included in the metrics JSON
 * (for plotting) but are excluded from passing_designs.txt and from the
 * filtered PDB output channel.
 */
process RFDIFFUSION_FILTER {
    tag "rfdiff_filter"
    label 'cpu'

    publishDir "${params.outdir}/rfdiffusion", mode: 'copy',
        pattern: '{rfdiffusion_metrics.json,filter_summary.json,passing_designs.txt}'
    publishDir "${params.outdir}/rfdiffusion", mode: 'copy',
        pattern: 'split/*.pdb'
    publishDir "${params.outdir}/rfdiffusion", mode: 'copy',
        pattern: 'passing/*.pdb'

    input:
    path design_pdbs
    path input_pdb
    val  contigs
    val  hotspot
    val  receptor_chain
    val  effector_chain
    val  contact_cutoff
    val  min_hotspot_frac
    path filter_script

    output:
    path "rfdiffusion_metrics.json",   emit: metrics
    path "filter_summary.json",        emit: filter_summary
    path "passing_designs.txt",        emit: passing_list
    path "passing/*.pdb",              emit: passing_pdbs, optional: true
    path "split/*.pdb",                emit: split_pdbs,   optional: true

    script:
    def hotspot_arg = hotspot ? "--hotspot \"${hotspot}\"" : ""
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${filter_script} \\
            --input-pdb ${input_pdb} \\
            --design-dir . \\
            --contigs "${contigs}" \\
            ${hotspot_arg} \\
            --receptor-chain ${receptor_chain} \\
            --effector-chain ${effector_chain} \\
            --contact-cutoff ${contact_cutoff} \\
            --min-hotspot-frac ${min_hotspot_frac}

    # ─── Copy passing designs to a subdirectory for the output channel ────
    mkdir -p passing
    if [ -s passing_designs.txt ]; then
        while IFS= read -r pdb; do
            cp "\${pdb}" passing/
        done < passing_designs.txt
    else
        echo "WARNING: No designs passed the filter"
    fi
    """
}


/*
 * RFDIFFUSION_PLOTS
 * ─────────────────
 * Three diagnostic plots from the RFDiffusion filter analysis:
 *   1. Design-region RMSD vs interface contacts (scatter, pass/fail coloured)
 *   2. Receptor contact position heatmap (2D, sorted by interface position)
 *   3. Fraction of contacts inside design region (bar chart with threshold)
 */
process RFDIFFUSION_PLOTS {
    tag "rfdiff_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path metrics
    path plots_script

    output:
    path "rfdiff_*.png", emit: plots

    script:
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rfdiff_container} \\
        python ${plots_script} \\
            --metrics ${metrics}
    """
}
