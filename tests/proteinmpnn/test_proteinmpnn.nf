#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_proteinmpnn.nf — Isolated test for the ProteinMPNN module
 * =============================================================================
 * Mirrors main.nf Branch B end-to-end up to the MPNN block:
 *
 *   1. Load the input complex PDB from this test's data/ directory.
 *   2. Call the *production* RESOLVE_CONTIGS to convert the user-supplied
 *      contig string into PDB-resolved coordinates.
 *   3. Call the *production* EXTRACT_SEQUENCES to extract receptor and
 *      effector sequences plus receptor_start_pdb from the input PDB.
 *   4. Feed everything into the MPNN block exactly as main.nf would.
 *
 * Both EXTRACT_SEQUENCES and RESOLVE_CONTIGS are imported from
 * modules/preprocessing — no inline duplicates.  This guarantees the test
 * exercises the same code paths the production pipeline does.
 *
 * Inputs (canonical, no upstream chaining)
 * ----------------------------------------
 *   data/input_complex.pdb
 *       Receptor+effector complex PDB.  Reference structure for
 *       RESOLVE_CONTIGS and EXTRACT_SEQUENCES.
 *   data/rosetta_passing/design_*.pdb
 *       Curated split-PDB fixtures captured from a discovery run's
 *       ROSETTA_FILTER ``passing/`` output.  See data/README.md.
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults
// ---------------------------------------------------------------------------

// Canonical input paths: this test's own data/ directory.  No fallback to
// upstream test outputs — per-module tests run from committed fixtures.
params.design_pdbs = "${projectDir}/data/rosetta_passing/design_*.pdb"
params.input_pdb   = "${projectDir}/data/input_complex.pdb"

params.receptor_chain    = "A"
params.effector_chain    = "C"
params.contigs           = "A1-32/10-20/A46-72/6-6 C"

params.num_seqs           = 2
params.mpnn_sampling_temp = 0.25
params.rm_aa              = "C"
params.mpnn_top_n         = 8

params.max_poly_x         = 5
params.min_pct_identity   = 0.0
params.max_pct_identity   = 100.0

params.project_name       = "test_proteinmpnn"
params.outdir             = "${projectDir}/results"

// ---------------------------------------------------------------------------
// Includes — production modules only, no inline duplicates
// ---------------------------------------------------------------------------

include { EXTRACT_SEQUENCES        } from '../../modules/preprocessing'
include { RESOLVE_CONTIGS          } from '../../modules/preprocessing'

include { MPNN_FIXED_POSITIONS     } from '../../modules/proteinmpnn'
include { PROTEINMPNN              } from '../../modules/proteinmpnn'
include { SEQUENCE_CORRECTION      } from '../../modules/proteinmpnn'
include { SEQUENCE_QC              } from '../../modules/proteinmpnn'
include { MPNN_DESIGN_REGION_SCORE } from '../../modules/proteinmpnn'
include { MPNN_CLUSTER             } from '../../modules/proteinmpnn'
include { MPNN_PLOTS               } from '../../modules/proteinmpnn'
include { MPNN_SELECT_TOP          } from '../../modules/proteinmpnn'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    log.info "ProteinMPNN test — input_pdb : ${params.input_pdb}"
    log.info "ProteinMPNN test — design_pdbs glob: ${params.design_pdbs}"

    input_pdb_ch = Channel.fromPath(params.input_pdb, checkIfExists: true)

    // ── Resolve the user contig against the input PDB ────────────────
    // Mirrors main.nf Branch B.
    RESOLVE_CONTIGS(
        input_pdb_ch,
        params.contigs,
        Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
    )
    contigs_ch = RESOLVE_CONTIGS.out.resolved_contigs
        .map { it.text.trim() }
        .first()

    // ── Extract sequences (incl. receptor_start_pdb) ─────────────────
    // Production EXTRACT_SEQUENCES from modules/preprocessing.  Reads
    // the user's input chain letters from params.receptor_chain /
    // params.effector_chain — these refer to the chain IDs in the
    // input PDB, not the A/B convention RFDiffusion uses for split
    // design PDBs.  The MPNN processes downstream handle the split-PDB
    // chain layout internally.
    EXTRACT_SEQUENCES(
        input_pdb_ch,
        params.receptor_chain,
        params.effector_chain
    )

    // Parse receptor_seq, effector_seq, receptor_start_pdb in one map
    // and split into value channels (.first() converts queue → value
    // so MPNN_FIXED_POSITIONS, which fans out across all designs, can
    // read each value repeatedly).
    seqs_ch = EXTRACT_SEQUENCES.out.sequences_json.map { json_file ->
        def slurp = new groovy.json.JsonSlurper()
        def data  = slurp.parseText(json_file.text)
        return [data.receptor_seq, data.effector_seq, data.receptor_start_pdb ?: 1]
    }.first()
    receptor_seq_ch       = seqs_ch.map { it[0] }
    effector_seq_ch       = seqs_ch.map { it[1] }
    receptor_start_pdb_ch = seqs_ch.map { it[2] }

    // ── Load design PDBs ─────────────────────────────────────────────
    design_pdbs_ch = Channel.fromPath(params.design_pdbs, checkIfExists: true)
    num_designs_ch = design_pdbs_ch.count()

    // ── Step 1: Fixed positions ──────────────────────────────────────
    MPNN_FIXED_POSITIONS(
        design_pdbs_ch,
        receptor_seq_ch,
        effector_seq_ch,
        contigs_ch,
        receptor_start_pdb_ch,
        Channel.value(file("${projectDir}/bin/pipeline_correct_sequences.py"))
    )

    // ── Step 2: ProteinMPNN ──────────────────────────────────────────
    PROTEINMPNN(
        MPNN_FIXED_POSITIONS.out.pdb_and_jsonl,
        params.num_seqs,
        params.mpnn_sampling_temp,
        params.rm_aa
    )

    // ── Step 3: Sequence correction ──────────────────────────────────
    SEQUENCE_CORRECTION(
        PROTEINMPNN.out.mpnn_results.collect(),
        receptor_seq_ch,
        effector_seq_ch,
        contigs_ch,
        num_designs_ch,
        params.num_seqs,
        receptor_start_pdb_ch,
        Channel.value(file("${projectDir}/bin/pipeline_correct_sequences.py"))
    )

    // ── Step 4: Quality control ──────────────────────────────────────
    SEQUENCE_QC(
        SEQUENCE_CORRECTION.out.af2_fastas,
        SEQUENCE_CORRECTION.out.metadata_csv,
        receptor_seq_ch,
        params.max_poly_x,
        params.min_pct_identity,
        params.max_pct_identity,
        Channel.value(file("${projectDir}/bin/mpnn_sequence_qc.py"))
    )

    // ── Step 5: Design-region scoring ────────────────────────────────
    design_pdbs_for_score = PROTEINMPNN.out.pdb_jsonl_mpnn.map { it[0] }.collect()
    fixed_jsonls_ch       = PROTEINMPNN.out.pdb_jsonl_mpnn.map { it[1] }.collect()
    mpnn_dirs_ch          = PROTEINMPNN.out.pdb_jsonl_mpnn.map { it[2] }.collect()

    MPNN_DESIGN_REGION_SCORE(
        SEQUENCE_QC.out.qc_metadata,
        design_pdbs_for_score,
        mpnn_dirs_ch,
        fixed_jsonls_ch,
        Channel.value(file("${projectDir}/bin/mpnn_design_region_score.py"))
    )

    // ── Step 6: Sequence clustering ──────────────────────────────────
    MPNN_CLUSTER(
        MPNN_DESIGN_REGION_SCORE.out.scored_metadata,
        Channel.value(file("${projectDir}/bin/mpnn_cluster_sequences.py"))
    )

    // ── Step 7: Diagnostic plots ─────────────────────────────────────
    MPNN_PLOTS(
        MPNN_DESIGN_REGION_SCORE.out.scored_metadata,
        MPNN_CLUSTER.out.cluster_counts,
        receptor_seq_ch,
        contigs_ch,
        Channel.value(file("${projectDir}/bin/mpnn_plots.py"))
    )

    // ── Step 8 (optional): Top-N selection ───────────────────────────
    if (params.mpnn_top_n > 0) {
        MPNN_SELECT_TOP(
            SEQUENCE_QC.out.qc_fastas,
            MPNN_DESIGN_REGION_SCORE.out.scored_metadata,
            params.mpnn_top_n,
            Channel.value(file("${projectDir}/bin/mpnn_select_top.py"))
        )
    }
}

workflow.onComplete {
    log.info """
    =============================================================
    ProteinMPNN Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "ProteinMPNN test failed: ${workflow.errorMessage}"
}
