/*
 * =============================================================================
 * ProteinMPNN module
 * =============================================================================
 * Fixed position generation, sequence design, correction, QC, design-region
 * scoring, sequence clustering, top-N selection, and diagnostic plots.
 *
 * All Python logic lives in bin/ scripts.  Scripts are called via
 * singularity exec so that all pipeline dependencies are available.
 *
 * Container notes (Singularity 3.8 + .img format on NBI HPC):
 *   - --bind $PWD:$PWD is REQUIRED to make work dir writable
 *   - MPLCONFIGDIR=/tmp for matplotlib cache (always writable)
 */


process MPNN_FIXED_POSITIONS {
    tag "${design_pdb.baseName}"
    label 'cpu'

    input:
    path design_pdb
    val  receptor_seq
    val  effector_seq
    val  contigs
    val  receptor_start_pdb
    path correct_script

    output:
    tuple path(design_pdb), path("fixed_positions_*.jsonl"), emit: pdb_and_jsonl

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${correct_script} \\
            --mode gen_fixed_positions \\
            --pdb_file ${design_pdb} \\
            --output_dir . \\
            --receptor_seq "${receptor_seq}" \\
            --effector_seq "${effector_seq}" \\
            --contigs "${contigs}" \\
            --receptor_chain "${params.receptor_chain}" \\
            --receptor_start_pdb ${receptor_start_pdb} \\
            --num_designs 1 \\
            --num_seqs ${params.num_seqs}
    """
}


/*
 * PROTEINMPNN
 * -----------
 * Run ProteinMPNN sequence design.  --save_score 1 saves per-position
 * log-probabilities to scores/*.npz for design-region scoring downstream.
 */
process PROTEINMPNN {
    tag "${pdb.baseName}"
    label 'cpu'

    publishDir "${params.outdir}/mpnn/${pdb.baseName}", mode: 'copy'

    input:
    tuple path(pdb), path(jsonl)
    val   num_seqs
    val   sampling_temp
    val   rm_aa

    output:
    tuple path(pdb), path(jsonl), path("${pdb.baseName}_mpnn/"), emit: pdb_jsonl_mpnn
    path "${pdb.baseName}_mpnn/", emit: mpnn_results

    script:
    def pdb_stem = pdb.baseName
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        ${params.rfdiff_container} \\
        python /opt/ProteinMPNN/protein_mpnn_run.py \\
            --pdb_path ${pdb} \\
            --out_folder ./${pdb_stem}_mpnn \\
            --num_seq_per_target ${num_seqs} \\
            --sampling_temp ${sampling_temp} \\
            --omit_AAs "${rm_aa}" \\
            --fixed_positions_jsonl ${jsonl} \\
            --save_score 1 \\
            --seed 42
    """
}


process SEQUENCE_CORRECTION {
    tag "seq_correction"
    label 'cpu'

    publishDir "${params.outdir}/sequences", mode: 'copy',
        pattern: '{af2_fastas/**,sequence_metadata.csv,mpnn_corrected.fasta}'

    input:
    path mpnn_dirs
    val  receptor_seq
    val  effector_seq
    val  contigs
    val  num_designs
    val  num_seqs
    val  receptor_start_pdb
    path correct_script

    output:
    path "af2_fastas/",              emit: af2_fastas
    path "sequence_metadata.csv",    emit: metadata_csv
    path "mpnn_corrected.fasta",     emit: corrected_fasta

    script:
    """
    mkdir -p mpnn_combined
    for d in ${mpnn_dirs}; do
        if [ -d "\$d" ]; then
            design_name=\$(basename "\$d" | sed 's/_mpnn\$//')
            mkdir -p "mpnn_combined/\${design_name}/seqs"
            if [ -d "\$d/seqs" ]; then
                cp "\$d/seqs/"*.fa "mpnn_combined/\${design_name}/seqs/" 2>/dev/null || true
            fi
        fi
    done

    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${correct_script} \\
            --mode correct \\
            --mpnn_dir mpnn_combined \\
            --output_dir . \\
            --receptor_seq "${receptor_seq}" \\
            --effector_seq "${effector_seq}" \\
            --contigs "${contigs}" \\
            --receptor_chain "${params.receptor_chain}" \\
            --receptor_start_pdb ${receptor_start_pdb} \\
            --num_designs ${num_designs} \\
            --num_seqs ${num_seqs}
    """
}


process SEQUENCE_QC {
    tag "seq_qc"
    label 'cpu'

    publishDir "${params.outdir}/sequences", mode: 'copy', pattern: 'qc_*'

    input:
    path af2_fastas_dir
    path metadata_csv
    val  receptor_seq
    val  max_poly_x
    val  min_pct_identity
    val  max_pct_identity
    path qc_script

    output:
    path "qc_fastas/",           emit: qc_fastas
    path "qc_metadata.csv",      emit: qc_metadata
    path "qc_report.txt",        emit: qc_report

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${qc_script} \\
            --af2-fastas-dir ${af2_fastas_dir} \\
            --metadata-csv ${metadata_csv} \\
            --receptor-seq "${receptor_seq}" \\
            --max-poly-x ${max_poly_x} \\
            --min-pct-identity ${min_pct_identity} \\
            --max-pct-identity ${max_pct_identity} \\
            --output-dir .
    """
}


/*
 * MPNN_DESIGN_REGION_SCORE
 * ------------------------
 * Run ProteinMPNN's model directly to score each designed sequence,
 * extracting per-position log-probabilities and computing the mean
 * over only the free (designable) positions.
 *
 * Requires: design PDBs, MPNN FASTAs (with designed sequences),
 * and fixed_positions JSONLs.
 */
process MPNN_DESIGN_REGION_SCORE {
    tag "design_region_score"
    label 'cpu'

    publishDir "${params.outdir}/sequences", mode: 'copy',
        pattern: 'scored_metadata.csv'

    input:
    path qc_metadata
    path design_pdbs
    path mpnn_dirs
    path fixed_jsonls
    path score_script

    output:
    path "scored_metadata.csv", emit: scored_metadata

    script:
    """
    # Collect design PDBs into a single directory
    mkdir -p pdb_collected
    for f in ${design_pdbs}; do
        if [ -f "\$f" ]; then
            ln -s "\$(readlink -f \$f)" "pdb_collected/\$(basename \$f)"
        fi
    done

    # Collect MPNN output dirs (contain seqs/ subdirs with FASTAs)
    mkdir -p mpnn_collected
    for d in ${mpnn_dirs}; do
        if [ -d "\$d" ]; then
            ln -s "\$(readlink -f \$d)" "mpnn_collected/\$(basename \$d)"
        fi
    done

    # Collect fixed position JSONLs
    mkdir -p fixed_collected
    for j in ${fixed_jsonls}; do
        if [ -f "\$j" ]; then
            ln -s "\$(readlink -f \$j)" "fixed_collected/\$(basename \$j)"
        fi
    done

    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${score_script} \\
            --metadata-csv ${qc_metadata} \\
            --pdb-dir pdb_collected \\
            --mpnn-fasta-dir mpnn_collected \\
            --fixed-positions-dir fixed_collected \\
            --output-csv scored_metadata.csv
    """
}


/*
 * MPNN_CLUSTER
 * ------------
 * Cluster design-region sequences at multiple identity thresholds
 * using MMseqs2.  Outputs cluster counts for diversity plotting.
 */
process MPNN_CLUSTER {
    tag "mmseqs_cluster"
    label 'cpu'

    publishDir "${params.outdir}/sequences", mode: 'copy',
        pattern: 'mpnn_cluster_counts.csv'

    input:
    path metadata
    path cluster_script

    output:
    path "mpnn_cluster_counts.csv", emit: cluster_counts

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${cluster_script} \\
            --metadata-csv ${metadata} \\
            --output-csv mpnn_cluster_counts.csv
    """
}


/*
 * MPNN_PLOTS
 * ----------
 * Diagnostic plots:
 *   1. Per-design global MPNN score boxplots (ranked by median)
 *   2. Per-design design-region MPNN score boxplots (same rank order)
 *   3. Sequence diversity: MMseqs2 clusters vs identity threshold
 *   4. AA composition: per-design heatmap + ranked boxplot (design region)
 */
process MPNN_PLOTS {
    tag "mpnn_plots"
    label 'cpu'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path metadata
    path cluster_csv
    val  receptor_seq
    val  contigs
    path plots_script

    output:
    path "mpnn_*.png", emit: plots

    script:
    """
    singularity exec \\
        --bind \${PWD}:\${PWD} \\
        --env MPLCONFIGDIR=/tmp \\
        ${params.rfdiff_container} \\
        python ${plots_script} \\
            --metadata ${metadata} \\
            --cluster-csv ${cluster_csv} \\
            --receptor-seq "${receptor_seq}" \\
            --contigs "${contigs}"
    """
}


process MPNN_SELECT_TOP {
    tag "select_top_${top_n}"
    label 'cpu'

    publishDir "${params.outdir}/sequences", mode: 'copy', pattern: 'top_*'

    input:
    path qc_fastas_dir
    path qc_metadata
    val  top_n
    path select_script

    output:
    path "top_fastas/",          emit: top_fastas
    path "top_metadata.csv",     emit: top_metadata
    path "selection_report.txt", emit: report

    script:
    """
    singularity exec --bind \${PWD}:\${PWD} ${params.rfdiff_container} \\
        python ${select_script} \\
            --qc-fastas-dir ${qc_fastas_dir} \\
            --qc-metadata ${qc_metadata} \\
            --top-n ${top_n} \\
            --output-dir .
    """
}
