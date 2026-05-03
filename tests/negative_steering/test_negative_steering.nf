#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_negative_steering.nf — Isolated test for the Negative-Steering module
 * =============================================================================
 * Runs the three processes in modules/negative_steering.nf in isolation:
 *
 *   NEGSTEER_DERIVE_INDICES  →  NEGSTEER_RUN_ONE  →  NEGSTEER_CROSS_SEQUENCE
 *
 * Uses pre-generated MPNN FASTAs, the parent RFDiffusion design PDBs, and
 * a real rfdiffusion_metrics.json as input — the test runs the full
 * single-cycle chain end-to-end on one GPU node per MPNN sequence.
 *
 * Inputs (canonical, no upstream chaining)
 * ----------------------------------------
 *   data/fastas/<design_stem>_seq_<N>.fasta
 *       MPNN-corrected FASTAs (>receptor + >effector).  One file per
 *       (design, sequence) pair.  Filename must match
 *       /^design_\\d+_seq_\\d+\\.fasta$/ for the parser regex below.
 *   data/design_pdbs/design_*.pdb
 *       Parent RFDiffusion design PDBs (chain A = Cα-only designed
 *       receptor, chain B = native effector).
 *   data/rfdiffusion_metrics.json
 *       RFDIFFUSION_FILTER metrics blob with one entry per design_*
 *       containing receptor_contact_residues, receptor_position_order,
 *       and per_design_design_residues.
 *   data/input_complex.pdb
 *       Receptor+effector complex PDB.  Used only when controls are
 *       enabled (DERIVE_INPUT_INDICES).  Must match the contigs string.
 *
 * The bin/ directory should be symlinked into the test dir so the
 * module can find it:
 *
 *   ln -s ../../bin tests/negative_steering/bin
 *
 * Usage:
 *   sbatch tests/negative_steering/run_test_negative_steering_slurm.sh
 *
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults — override via params.yml or --param on the command line
// ---------------------------------------------------------------------------

params.project_name   = "test_negative_steering"
params.outdir         = "${projectDir}/results"

// Canonical input paths: this test's own data/ directory.  No fallback to
// upstream test outputs — per-module tests run from committed fixtures.
params.fastas_dir          = "${projectDir}/data/fastas"
params.design_pdb_dir      = "${projectDir}/data/design_pdbs"
params.rfdiffusion_metrics = "${projectDir}/data/rfdiffusion_metrics.json"
params.input_pdb           = "${projectDir}/data/input_complex.pdb"

// ── Chain identifiers ─────────────────────────────────────────────────
// receptor_chain / effector_chain refer to the INPUT PDB convention
// (the receptor+effector complex going into RFDiffusion).  For this
// test target, the input PDB has effector on chain C, so we set
// effector_chain="C" — matching the contigs string and the documented
// convention in params_example.yml.
//
// rfdiff_output_* refer to the post-RFDiffusion split convention,
// which is HARDCODED to A/B by rfdiffusion_filter.py.  Every process
// that consumes split design PDBs (NEGSTEER_RUN_ONE etc.) must use
// these labels, NOT receptor_chain / effector_chain.
params.receptor_chain               = "A"
params.effector_chain               = "C"
params.rfdiff_output_receptor_chain = "A"
params.rfdiff_output_effector_chain = "B"

// ── Steering strategy ──────────────────────────────────────────────────
// Default to the same flags as the production main.nf so the test
// exercises the production path.  Reduce n_designs for faster testing
// if needed via --negsteer_n_designs <N>.
params.negsteer_mode                 = "mild"
params.negsteer_max_mutations        = 6
params.negsteer_candidate_pool_size  = 10
params.negsteer_protected_set_source = "design_region_union"
params.negsteer_n_designs            = 4
params.negsteer_num_seeds            = 3
params.negsteer_n_cycles             = 1

// ── Boltz prediction hyperparameters ──────────────────────────────────
// Test defaults match production; lower diffusion_samples/recycling if
// the test box is slow.
params.negsteer_diffusion_samples    = 5
params.negsteer_recycling_steps      = 3
params.negsteer_rmsd_threshold       = 6.0
params.negsteer_contact_cutoff       = 4.5

// Match production — P0-35 closed with kernels DISABLED because the
// cuEquivariance path is ~5% slower on the ~155-residue complexes.
// See pipeline_notes7 for the measurement.
params.negsteer_no_kernels           = true

// ── Post-processing thresholds ────────────────────────────────────────
params.negsteer_postprocess_rmsd_threshold = 5.0
params.negsteer_postprocess_metric_column  = "steered_ra_eff_vs_truth"
params.negsteer_postprocess_contact_cutoff = 5.0

// ── Negative controls (Task 6) ────────────────────────────────────────
// Default ON to validate the controls path end-to-end alongside the
// steered runs.  The 2 control GPU jobs run in parallel with the main
// steered fan-out (max_boltz2_parallel covers them), so they don't add
// to wall time on a non-contended queue.  Per-sequence runtimes are
// measured separately so this doesn't pollute the inner-loop timing.
//
// Controls require:
//   - input_pdb  : receptor+effector complex PDB (the same one fed to
//                  RFDiffusion to define the design problem)
//   - contigs    : the contigs string passed to RFDiffusion, so the
//                  controls' design region matches the steered runs'.
params.run_negative_controls         = true
params.negsteer_controls_n_designs   = 1
// Contigs string for DERIVE_INPUT_INDICES.  Trailing chain-letter
// token "C" matches the input PDB's effector chain (params.effector_chain
// above), per the documented convention.
params.contigs                       = "A1-32/10-20/A46-72/6-6 C"

// ── Infrastructure ────────────────────────────────────────────────────
// params.boltz2_container is inherited from nextflow.config.
params.max_boltz2_parallel     = 24

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { NEGSTEER_DERIVE_INDICES         } from '../../modules/negative_steering'
include { NEGSTEER_RUN_ONE                } from '../../modules/negative_steering'
include { NEGSTEER_CROSS_SEQUENCE         } from '../../modules/negative_steering'
include { NEGSTEER_PLOTS                  } from '../../modules/negative_steering'
include { NEGSTEER_WITHIN_SEQUENCE_PLOTS  } from '../../modules/negative_steering'
include { DERIVE_INPUT_INDICES            } from '../../modules/negsteer_controls'
include { NEGSTEER_CONTROLS               } from '../../modules/negsteer_controls'


/*
 * EMPTY_DESIGN_REGION_PLACEHOLDER
 * ───────────────────────────────
 * Emits a 0-byte file named input_design_region.txt to stand in for
 * DERIVE_INPUT_INDICES.out.design_region when controls are disabled
 * (DERIVE_INPUT_INDICES is only invoked under `if (run_negative_controls)`
 * in this test workflow).
 *
 * NEGSTEER_PLOTS expects a `path input_design_region` input
 * unconditionally — its mutation-impact plot uses the file to shade
 * design-region positions on the input PDB.  bin/negsteer_plots.py
 * handles an empty file gracefully (parser returns an empty position
 * set; plot renders without the shading), so the placeholder is the
 * minimal way to keep the plot process always-on without changing
 * the production process signature.
 */
process EMPTY_DESIGN_REGION_PLACEHOLDER {
    tag "empty_design_region"
    label 'cpu'

    output:
    path "input_design_region.txt", emit: design_region

    script:
    """
    : > input_design_region.txt
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
//
// Same shape as the Step-4 block in main.nf, but with the upstream
// RFDIFFUSION/MPNN stages replaced by static file inputs.

workflow {

    log.info "Negative-steering test — input sources:"
    log.info "  fastas_dir           : ${params.fastas_dir}"
    log.info "  design_pdb_dir       : ${params.design_pdb_dir}"
    log.info "  rfdiffusion_metrics  : ${params.rfdiffusion_metrics}"
    log.info "  input_pdb            : ${params.input_pdb}"
    log.info "  run_negative_controls: ${params.run_negative_controls}"
    log.info "  negsteer_no_kernels  : ${params.negsteer_no_kernels}"

    // ── Resolve inputs ────────────────────────────────────────────────
    fastas_dir_ch   = Channel.fromPath(params.fastas_dir, checkIfExists: true, type: 'dir')
    metrics_json_ch = Channel.fromPath(params.rfdiffusion_metrics, checkIfExists: true)

    // ── FASTA stream [design_stem, seq_name, mpnn_fasta] ──────────────
    // Same parsing logic as in main.nf so the test exercises the exact
    // production regex.
    fasta_stream_ch = fastas_dir_ch
        .flatMap { fdir ->
            fdir.listFiles()
                .findAll { it.name.endsWith(".fasta") }
                .collect { f ->
                    def stem = f.name.replaceAll(/\.fasta$/, '')
                    def m = stem =~ /^(design_\d+)_seq_\d+$/
                    def design_stem = m ? m[0][1] : null
                    [design_stem, stem, f]
                }
        }
        .filter { trip -> trip[0] != null }

    // ── Design-PDB map [design_stem, design_pdb] ──────────────────────
    // Glob every design_*.pdb under the data dir and key each by its
    // basename ("design_2", "design_3", etc.).  Used both to feed the
    // index-derivation step (one task per design) and to join into
    // the per-sequence input tuple downstream.
    design_pdb_map_ch = Channel.fromPath("${params.design_pdb_dir}/design_*.pdb", checkIfExists: true)
        .map { pdb -> [pdb.baseName, pdb] }

    // ── Derive index files for every PDB unconditionally ──────────────
    // Earlier versions filtered down to "designs that have FASTAs"
    // before deriving indices, which required reading fasta_stream_ch
    // twice (once for .unique() of stems, once for the three-way join
    // below).  Queue channels are single-consumer in Nextflow, which
    // caused half the per-sequence runs to silently disappear.  The
    // fix is to skip the filter entirely — index derivation is cheap
    // (~5s per design, stdlib + JSON) and removing the multi-consumer
    // dependency makes the topology bulletproof.  Designs without a
    // matching FASTA simply have their indices computed and discarded
    // by the downstream .combine(by:0) join.
    NEGSTEER_DERIVE_INDICES(
        design_pdb_map_ch.map { stem, pdb -> pdb },
        metrics_json_ch.first(),
        Channel.value(file("${projectDir}/bin/derive_design_region.py")),
        Channel.value(file("${projectDir}/bin/derive_true_interface.py"))
    )

    indices_map_ch = NEGSTEER_DERIVE_INDICES.out.indices

    // ── Three-way join → per-sequence input tuple ─────────────────────
    run_one_input_ch = fasta_stream_ch
        .combine(design_pdb_map_ch, by: 0)
        .combine(indices_map_ch,    by: 0)
        .map { stem, seq_name, fasta, pdb, region, iface ->
            tuple(seq_name, stem, fasta, pdb, region, iface)
        }

    // ── Compose plan_extra_args from params ───────────────────────────
    plan_extra_args = [
        "--mode ${params.negsteer_mode}",
        "--max-mutations ${params.negsteer_max_mutations}",
        "--candidate-pool-size ${params.negsteer_candidate_pool_size}",
        "--protected-set-source ${params.negsteer_protected_set_source}",
        "--n-designs ${params.negsteer_n_designs}",
        "--num-seeds ${params.negsteer_num_seeds}",
        "--diffusion-samples ${params.negsteer_diffusion_samples}",
        "--recycling-steps ${params.negsteer_recycling_steps}",
        "--rmsd-threshold ${params.negsteer_rmsd_threshold}",
        "--contact-cutoff ${params.negsteer_contact_cutoff}",
        params.negsteer_no_kernels ? "--no-kernels" : "",
    ].findAll { it }.join(" ")

    // ── Fan out — one GPU job per MPNN sequence ──────────────────────
    // NEGSTEER_RUN_ONE consumes split design PDBs from RFDiffusion,
    // which always have receptor=A, effector=B (hardcoded by
    // rfdiffusion_filter.py:write_split_pdb).  Use the rfdiff_output_*
    // params here, not params.receptor_chain / params.effector_chain
    // (which describe the input PDB convention, A+C for this target).
    NEGSTEER_RUN_ONE(
        run_one_input_ch,
        params.rfdiff_output_receptor_chain,
        params.rfdiff_output_effector_chain,
        plan_extra_args,
        params.negsteer_postprocess_rmsd_threshold,
        params.negsteer_postprocess_metric_column,
        params.negsteer_postprocess_contact_cutoff,
        Channel.value(file("${projectDir}/bin/negative_steering_run_one.sh"))
    )

    // ── Negative controls (Task 6) ────────────────────────────────────
    // Mirrors the controls block in main.nf.  Two control GPU jobs:
    // scrambled receptor + polyA receptor.  Both should fail to
    // produce passing rows under negative-steering — that's the
    // diagnostic point.  If a control accidentally tier-A's, that's
    // a strong signal of pipeline rubber-stamping.
    if (params.run_negative_controls) {

        controls_plan_extra_args = [
            "--mode ${params.negsteer_mode}",
            "--max-mutations ${params.negsteer_max_mutations}",
            "--candidate-pool-size ${params.negsteer_candidate_pool_size}",
            "--protected-set-source ${params.negsteer_protected_set_source}",
            "--n-designs ${params.negsteer_controls_n_designs}",
            "--num-seeds ${params.negsteer_num_seeds}",
            "--diffusion-samples ${params.negsteer_diffusion_samples}",
            "--recycling-steps ${params.negsteer_recycling_steps}",
            "--rmsd-threshold ${params.negsteer_rmsd_threshold}",
            "--contact-cutoff ${params.negsteer_contact_cutoff}",
            params.negsteer_no_kernels ? "--no-kernels" : "",
        ].findAll { it }.join(" ")

        // Derive the input-structure design region + true interface once.
        // Channels for the input PDB and contigs string — single-emission
        // value channels so DERIVE_INPUT_INDICES sees one task.
        input_pdb_ch = Channel.fromPath(params.input_pdb, checkIfExists: true)
        contigs_ch   = Channel.value(params.contigs)

        DERIVE_INPUT_INDICES(
            input_pdb_ch,
            contigs_ch,
            params.receptor_chain,
            params.effector_chain,
            params.negsteer_contact_cutoff,
            Channel.value(file("${projectDir}/bin/derive_input_design_region.py"))
        )

        // Build a 2-element input channel for NEGSTEER_CONTROLS — same
        // Cartesian-product pattern as main.nf:
        //   2 control types × 1 input bundle = 2 process invocations.
        controls_input_ch = Channel.of(
                ["input_control_scrambled", "scrambled"],
                ["input_control_polyA",     "polyA"]
            )
            .combine(DERIVE_INPUT_INDICES.out.input_pdb.first())
            .combine(DERIVE_INPUT_INDICES.out.design_region.first())
            .combine(DERIVE_INPUT_INDICES.out.true_interface.first())
            .map { name, type, pdb, region, iface ->
                tuple(name, type, pdb, region, iface)
            }

        NEGSTEER_CONTROLS(
            controls_input_ch,
            params.receptor_chain,
            params.effector_chain,
            controls_plan_extra_args,
            params.negsteer_postprocess_rmsd_threshold,
            params.negsteer_postprocess_metric_column,
            params.negsteer_postprocess_contact_cutoff,
            Channel.value(file("${projectDir}/bin/build_control_sequences.py")),
            Channel.value(file("${projectDir}/bin/negative_steering_run_one.sh"))
        )

        // Mix steered + control workdirs into one channel feeding the
        // cross-sequence aggregator.  cross_sequence_summary.py picks
        // up the row_type.txt sidecar from each workdir and tags rows
        // appropriately; control rows are excluded from cross-rank
        // scoring but appear in the output for diagnosis.
        all_workdirs_ch = NEGSTEER_RUN_ONE.out.per_sequence_workdir
            .mix(NEGSTEER_CONTROLS.out.per_control_workdir)

    } else {
        all_workdirs_ch = NEGSTEER_RUN_ONE.out.per_sequence_workdir
    }

    // ── Cross-sequence aggregation ───────────────────────────────────
    per_sequence_workdirs_ch = all_workdirs_ch.collect()
    NEGSTEER_CROSS_SEQUENCE(
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/cross_sequence_summary.py"))
    )

    // ── Diagnostic plots ─────────────────────────────────────────────
    // Mirrors main.nf "Step 4b": cohort + within-sequence plots in
    // parallel.  DERIVE_INPUT_INDICES only runs when controls are
    // enabled, so when it doesn't, an empty placeholder file feeds
    // the cohort plot's --input-design-region slot.  bin/negsteer_plots.py
    // handles the empty file by rendering the mutation-impact plot
    // without design-region shading.
    if (params.run_negative_controls) {
        design_region_ch = DERIVE_INPUT_INDICES.out.design_region.first()
    } else {
        EMPTY_DESIGN_REGION_PLACEHOLDER()
        design_region_ch = EMPTY_DESIGN_REGION_PLACEHOLDER.out.design_region.first()
    }
    NEGSTEER_PLOTS(
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        per_sequence_workdirs_ch,
        design_region_ch,
        Channel.value(file("${projectDir}/bin/negsteer_plots.py"))
    )
    NEGSTEER_WITHIN_SEQUENCE_PLOTS(
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/negsteer_within_sequence_plots.py"))
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    =============================================================
    Negative-Steering Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Negative-steering test failed: ${workflow.errorMessage}"
}
