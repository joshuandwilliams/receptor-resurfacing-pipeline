#!/usr/bin/env nextflow

/*
 * =============================================================================
 * Receptor Resurfacing Pipeline — MPNN Resample Sibling Workflow
 * =============================================================================
 *
 * Regenerates ProteinMPNN sequences for a SINGLE, already-generated
 * RFDiffusion design and carries them through the rest of the pipeline
 * (correction, QC, Boltz2 negative-steering, orthogonal validation)
 * exactly as main.nf does from Step 3 onward.
 *
 * Use case: a design has already been wet-lab tested and you want more
 * sequence variants of the SAME backbone (same contig, same realised
 * design region) rather than a fresh RFDiffusion batch.
 *
 * Skips RFDiffusion + Rosetta pre-validation entirely by pointing
 * straight at the already-published artefacts from the source run:
 *   <source_run_dir>/rosetta_filtering/passing/<design_id>.pdb
 *     — the Rosetta-passing split PDB ProteinMPNN designs against.
 *   <source_run_dir>/rfdiffusion/<design_id>.pdb
 *     — the raw (un-split) RFDiffusion output used as negsteer's
 *       ground truth.
 *   <source_run_dir>/rfdiffusion/rfdiffusion_metrics.json
 *     — batch metrics NEGSTEER_DERIVE_INDICES reads to recover the
 *       design region / true interface for <design_id> (the file
 *       still contains every design from the source run; only the
 *       one record matching design_id is used).
 *
 * params.pdb_file and params.contigs must be the SAME complex PDB and
 * contig spec the source run used — they regenerate the fixed-position
 * mask (which depends on contigs + the realised PDB, not on anything
 * RFDiffusion-specific) and drive the negative-control derivation,
 * exactly as they do in main.nf's Branch B.
 *
 * IMPORTANT: Step 3 (ProteinMPNN) → Step 5 (orthogonal metrics) below
 * is a straight copy of main.nf's wiring — same modules, same channel
 * shapes.  If that wiring changes in main.nf, mirror the change here.
 *
 * Author: Josh Williams
 * Institute: John Innes Centre / The Sainsbury Laboratory
 * =============================================================================
 */

nextflow.enable.dsl = 2

import groovy.json.JsonOutput

// ── Project ─────────────────────────────────────────────────────────────
params.project_name = "receptor_resurfacing_resample"
params.outdir        = "${launchDir}/${params.project_name}_results"
if (!new File(params.outdir.toString()).isAbsolute()) {
    params.outdir = file(params.outdir).toAbsolutePath().normalize().toString()
}

// ── Source design to resample ─────────────────────────────────────────
// source_run_dir   : the ORIGINAL run's results/ dir that produced the
//                     design being resampled, e.g.
//                     ".../runs/crystal_full_test_contig/results"
// source_design_id : RFDiffusion design stem, e.g. "design_26"
params.source_run_dir   = null
params.source_design_id = null

// ── Original docking inputs (MUST match the source run exactly) ───────
params.pdb_file        = null   // Pre-docked complex PDB the source run used
params.contigs          = null   // Same contig spec string the source run used
params.receptor_chain    = "A"
params.effector_chain    = "B"

// ── Post-RFDiffusion chain convention (see main.nf) ────────────────────
// RFDiffusion's split function always writes A=receptor, B=effector
// regardless of the input PDB's chain labels.  The source run's split/
// passing PDBs follow this convention, so NEGSTEER_RUN_ONE needs it too.
params.rfdiff_output_receptor_chain = "A"
params.rfdiff_output_effector_chain = "B"

// ── ProteinMPNN ─────────────────────────────────────────────────────────
params.num_seqs           = 8
params.mpnn_sampling_temp = 0.25
params.rm_aa               = "C"
params.mpnn_top_n          = 0

// ── ProteinMPNN sequence quality control ──────────────────────────────
params.max_poly_x       = 5
params.min_pct_identity = 0.0
params.max_pct_identity = 100.0

// ═════════════════════════════════════════════════════════════════════════
// Boltz2 negative-steering (identical defaults to main.nf)
// ═════════════════════════════════════════════════════════════════════════
params.negsteer_mode                 = "mild"
params.negsteer_max_mutations        = 6
params.negsteer_candidate_pool_size  = 10
params.negsteer_protected_set_source = "design_region_union"
params.negsteer_n_designs            = 20
params.negsteer_num_seeds            = 3
params.negsteer_n_cycles             = 1
params.negsteer_diffusion_samples    = 5
params.negsteer_recycling_steps      = 3
params.negsteer_rmsd_threshold       = 6.0
params.negsteer_contact_cutoff       = 4.5
params.negsteer_no_kernels           = false

params.negsteer_postprocess_rmsd_threshold = 5.0
params.negsteer_postprocess_metric_column  = "steered_ra_eff_vs_truth"
params.negsteer_postprocess_contact_cutoff = 5.0

// ── Negative controls ─────────────────────────────────────────────────
params.run_negative_controls       = true
params.negsteer_controls_n_designs = 1
params.controls_warning_ipsae_max  = 0.5
params.controls_warning_ra_eff_min = 5.0

// ── Interface-restricted metrics ────────────────────────────────────
params.interface_plddt_trim_threshold = 50.0
params.interface_intact_threshold     = 5.0

// ── Weighted Jaccard contact overlap ────────────────────────────────
params.weighted_jaccard_pair_cutoff = 8.0

// ── Orthogonal validation metrics ───────────────────────────────────
params.af3_nomsa_seeds = [42, 123, 456]

params.orthogonal_filter_sc_min     = 0.55
params.orthogonal_filter_bsa_min    = 600
params.orthogonal_filter_plddt_min  = 0.75
params.orthogonal_filter_af3_ra_max = 5.0

// ── Execution ─────────────────────────────────────────────────────────
params.max_boltz2_parallel = 5

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------
if (!params.source_run_dir) {
    error "Provide --source_run_dir (the ORIGINAL run's results/ dir that produced the design being resampled)"
}
if (!params.source_design_id) {
    error "Provide --source_design_id (RFDiffusion design stem, e.g. 'design_26')"
}
if (!params.pdb_file) {
    error "Provide --pdb_file (the pre-docked complex PDB the source run used)"
}
if (!params.contigs) {
    error "Provide --contigs (the SAME contig spec string the source run used)"
}

// bin/pipeline_correct_sequences.py's `correct` mode (called by
// SEQUENCE_CORRECTION below) assumes designs are named design_0.. design_
// {num_designs-1} sequentially and loops `for m in range(num_designs)` —
// it does not scan mpnn_dir for whatever design names actually exist.
// source_design_id is an arbitrary pre-existing index (e.g. "design_26"),
// so num_designs must be at least source_design_index + 1 for the loop to
// ever reach it; the intervening m values (0..25) just produce harmless
// "not found, skipping" warnings.
def _design_idx_match = (params.source_design_id =~ /^design_(\d+)$/)
if (!_design_idx_match) {
    error "params.source_design_id must look like 'design_<N>' (got '${params.source_design_id}')"
}
def source_design_index = _design_idx_match[0][1].toInteger()

// ---------------------------------------------------------------------------
// Include modules — identical set to main.nf's Step 3 → Step 5
// ---------------------------------------------------------------------------
include { EXTRACT_SEQUENCES } from './modules/preprocessing'
include { RESOLVE_CONTIGS   } from './modules/preprocessing'

include { MPNN_FIXED_POSITIONS     } from './modules/proteinmpnn'
include { PROTEINMPNN               } from './modules/proteinmpnn'
include { SEQUENCE_CORRECTION       } from './modules/proteinmpnn'
include { SEQUENCE_QC               } from './modules/proteinmpnn'
include { MPNN_DESIGN_REGION_SCORE  } from './modules/proteinmpnn'
include { MPNN_CLUSTER              } from './modules/proteinmpnn'
include { MPNN_SELECT_TOP           } from './modules/proteinmpnn'
include { MPNN_PLOTS                } from './modules/proteinmpnn'

include { NEGSTEER_DERIVE_INDICES        } from './modules/negative_steering'
include { NEGSTEER_RUN_ONE                } from './modules/negative_steering'
include { NEGSTEER_CROSS_SEQUENCE         } from './modules/negative_steering'
include { NEGSTEER_PLOTS                  } from './modules/negative_steering'
include { NEGSTEER_WITHIN_SEQUENCE_PLOTS  } from './modules/negative_steering'
include { NEGSTEER_CONTROLS               } from './modules/negsteer_controls'
include { DERIVE_INPUT_INDICES            } from './modules/negsteer_controls'
include { NEGSTEER_INTERFACE_METRICS      } from './modules/negsteer_interface_metrics'
include { EXTRACT_SURVIVOR_MANIFEST       } from './modules/negsteer_manifest'
include { AF3_SETUP_DB                    } from './modules/negsteer_af3_nomsa'
include { AF3_NOMSA_ON_SURVIVORS          } from './modules/negsteer_af3_nomsa'
include { AF3_PARSE_OUTPUT                } from './modules/negsteer_af3_nomsa'
include { NEGSTEER_BIOPHYSICAL_METRICS    } from './modules/negsteer_biophysical_metrics'
include { NEGSTEER_ROSETTA_METRICS        } from './modules/negsteer_rosetta_metrics'
include { NEGSTEER_ORTHOGONAL_METRICS     } from './modules/negsteer_orthogonal_metrics'
include { ORTHOG_PLOTS                    } from './modules/negsteer_orthogonal_metrics'


workflow {

    // ── Parameter validation ─────────────────────────────────────────────
    // Same validator main.nf uses.  Keys with no spec here (RFDiffusion /
    // pose-solver params this workflow never sets) are simply absent from
    // the JSON dump and skipped — see bin/validate_params.py:validate_params.
    def _params_dump = file("${workflow.workDir}/.params_validation.json")
    _params_dump.parent.mkdirs()
    _params_dump.text = JsonOutput.toJson(params)
    def _validator = "${projectDir}/bin/validate_params.py"
    def _proc = ["python3", _validator, "--params-json", _params_dump.toString()]
        .execute()
    _proc.waitFor()
    def _stderr = _proc.err.text
    def _stdout = _proc.in.text
    if (_proc.exitValue() != 0) {
        error("Parameter validation failed (see bin/validate_params.py PARAM_SPECS):\n" +
              _stderr + _stdout)
    }
    log.info(_stderr.trim())

    // =====================================================================
    // Resolve the ORIGINAL docking inputs — mirrors main.nf Branch B, so
    // the fixed-position mask and negsteer's input-structure controls are
    // byte-identical to what the source run used.
    // =====================================================================
    rfdiff_pdb_ch = Channel.fromPath(params.pdb_file, checkIfExists: true)

    RESOLVE_CONTIGS(
        rfdiff_pdb_ch,
        params.contigs,
        Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
    )
    contigs_ch = RESOLVE_CONTIGS.out.resolved_contigs
        .map { it.text.trim() }
        .first()

    EXTRACT_SEQUENCES(
        rfdiff_pdb_ch,
        params.receptor_chain,
        params.effector_chain
    )

    seqs_ch = EXTRACT_SEQUENCES.out.sequences_json.map { json_file ->
        def slurp = new groovy.json.JsonSlurper()
        def data  = slurp.parseText(json_file.text)
        return [data.receptor_seq, data.effector_seq, data.receptor_start_pdb ?: 1]
    }.first()
    receptor_seq_ch       = seqs_ch.map { it[0] }
    effector_seq_ch       = seqs_ch.map { it[1] }
    receptor_start_pdb_ch = seqs_ch.map { it[2] }

    // =====================================================================
    // Step 0 (replaces RFDiffusion + Rosetta): point straight at the
    // already-published artefacts for the single design being resampled.
    // =====================================================================
    rosetta_passing_pdbs_ch = Channel.fromPath(
        "${params.source_run_dir}/rosetta_filtering/passing/${params.source_design_id}.pdb",
        checkIfExists: true
    )

    // "RFDIFFUSION.out.design_pdbs" equivalent — raw design PDB used as
    // negsteer's ground truth.
    design_pdbs_raw_ch = Channel.fromPath(
        "${params.source_run_dir}/rfdiffusion/${params.source_design_id}.pdb",
        checkIfExists: true
    )

    // "RFDIFFUSION_FILTER.out.metrics" equivalent — batch metrics JSON;
    // NEGSTEER_DERIVE_INDICES only reads the record for source_design_id.
    rfdiff_metrics_ch = Channel.fromPath(
        "${params.source_run_dir}/rfdiffusion/rfdiffusion_metrics.json",
        checkIfExists: true
    )

    // =====================================================================
    // Step 3: ProteinMPNN            (verbatim main.nf wiring from here)
    // =====================================================================
    MPNN_FIXED_POSITIONS(
        rosetta_passing_pdbs_ch,
        receptor_seq_ch,
        effector_seq_ch,
        contigs_ch,
        receptor_start_pdb_ch,
        Channel.value(file("${projectDir}/bin/pipeline_correct_sequences.py"))
    )

    PROTEINMPNN(
        MPNN_FIXED_POSITIONS.out.pdb_and_jsonl,
        params.num_seqs,
        params.mpnn_sampling_temp,
        params.rm_aa
    )

    SEQUENCE_CORRECTION(
        PROTEINMPNN.out.mpnn_results.collect(),
        receptor_seq_ch,
        effector_seq_ch,
        contigs_ch,
        source_design_index + 1,   // see note above: must reach source_design_id's index
        params.num_seqs,
        receptor_start_pdb_ch,
        Channel.value(file("${projectDir}/bin/pipeline_correct_sequences.py"))
    )

    SEQUENCE_QC(
        SEQUENCE_CORRECTION.out.af2_fastas,
        SEQUENCE_CORRECTION.out.metadata_csv,
        receptor_seq_ch,
        params.max_poly_x,
        params.min_pct_identity,
        params.max_pct_identity,
        Channel.value(file("${projectDir}/bin/mpnn_sequence_qc.py"))
    )

    // ── Design-region scoring ────────────────────────────────────────
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

    // ── Sequence clustering ──────────────────────────────────────────
    MPNN_CLUSTER(
        MPNN_DESIGN_REGION_SCORE.out.scored_metadata,
        Channel.value(file("${projectDir}/bin/mpnn_cluster_sequences.py"))
    )

    MPNN_PLOTS(
        MPNN_DESIGN_REGION_SCORE.out.scored_metadata,
        MPNN_CLUSTER.out.cluster_counts,
        receptor_seq_ch,
        contigs_ch,
        Channel.value(file("${projectDir}/bin/mpnn_plots.py"))
    )

    // =====================================================================
    // Step 3b: Select top MPNN sequences (if mpnn_top_n > 0)
    // =====================================================================
    if (params.mpnn_top_n > 0) {
        MPNN_SELECT_TOP(
            SEQUENCE_QC.out.qc_fastas,
            MPNN_DESIGN_REGION_SCORE.out.scored_metadata,
            params.mpnn_top_n,
            Channel.value(file("${projectDir}/bin/mpnn_select_top.py"))
        )
        mpnn_fastas_dir_ch = MPNN_SELECT_TOP.out.top_fastas
    } else {
        mpnn_fastas_dir_ch = SEQUENCE_QC.out.qc_fastas
    }

    // =====================================================================
    // Step 4: Boltz2 Negative-Steering Validation (per MPNN sequence)
    // =====================================================================

    // -- 4a: FASTA stream [design_stem, seq_name, mpnn_fasta] --------
    fasta_stream_ch = mpnn_fastas_dir_ch
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

    // -- 4b: Design-PDB map [design_stem, design_pdb] ----------------
    design_pdb_map_ch = design_pdbs_raw_ch
        .map { pdb -> [pdb.baseName, pdb] }   // "design_26" -> file

    // -- 4c: Derive index files for the resampled design -------------
    NEGSTEER_DERIVE_INDICES(
        design_pdb_map_ch.map { stem, pdb -> pdb },
        rfdiff_metrics_ch.first(),   // value-channel broadcast
        Channel.value(file("${projectDir}/bin/derive_design_region.py")),
        Channel.value(file("${projectDir}/bin/derive_true_interface.py"))
    )

    indices_map_ch = NEGSTEER_DERIVE_INDICES.out.indices

    // -- 4d: Three-way join into the per-sequence input tuple --------
    run_one_input_ch = fasta_stream_ch
        .combine(design_pdb_map_ch, by: 0)
        .combine(indices_map_ch,    by: 0)
        .map { stem, seq_name, fasta, pdb, region, iface ->
            tuple(seq_name, stem, fasta, pdb, region, iface)
        }

    // -- 4e: Compose plan_extra_args string from params --------------
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

    if (params.negsteer_n_cycles != 1) {
        log.warn "params.negsteer_n_cycles = ${params.negsteer_n_cycles}; " +
                 "the NextFlow integration currently runs single-cycle only."
    }

    // -- 4f: Fan out - one GPU SLURM job per MPNN sequence -----------
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

    // -- 4f-bis: Negative controls ------------------------------------
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

        DERIVE_INPUT_INDICES(
            rfdiff_pdb_ch,
            contigs_ch,
            params.receptor_chain,
            params.effector_chain,
            params.negsteer_contact_cutoff,
            Channel.value(file("${projectDir}/bin/derive_input_design_region.py"))
        )

        controls_input_ch = Channel.of(
                ["input_control_scrambled", "scrambled"],
                ["input_control_polyA",     "polyA"]
            )
            .combine(rfdiff_pdb_ch.first())
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

        all_workdirs_ch = NEGSTEER_RUN_ONE.out.per_sequence_workdir
            .mix(NEGSTEER_CONTROLS.out.per_control_workdir)

    } else {
        all_workdirs_ch = NEGSTEER_RUN_ONE.out.per_sequence_workdir
    }

    // -- 4g: Cross-sequence aggregation ------------------------------
    per_sequence_workdirs_ch = all_workdirs_ch.collect()

    NEGSTEER_CROSS_SEQUENCE(
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin")),
        MPNN_DESIGN_REGION_SCORE.out.scored_metadata
    )

    NEGSTEER_PLOTS(
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        per_sequence_workdirs_ch,
        DERIVE_INPUT_INDICES.out.design_region.first(),
        Channel.value(file("${projectDir}/bin/negsteer_plots.py"))
    )
    NEGSTEER_WITHIN_SEQUENCE_PLOTS(
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/negsteer_within_sequence_plots.py"))
    )

    // =====================================================================
    // Step 5: Post-negative-steering orthogonal metrics
    // =====================================================================
    NEGSTEER_INTERFACE_METRICS(
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/compute_interface_metrics.py"))
    )
    extended_csv_ch = NEGSTEER_INTERFACE_METRICS.out.extended_csv

    EXTRACT_SURVIVOR_MANIFEST(
        extended_csv_ch,
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/extract_survivor_manifest.py"))
    )

    manifest_records_ch = EXTRACT_SURVIVOR_MANIFEST.out.manifest
        .splitCsv(header: true)
        .map { r ->
            [
                r.seq_name,
                r.canonical_pdb_abs,
                r.ground_truth_abs,
                r.effector_template_cif_abs,
                r.receptor_seq,
                r.effector_seq,
            ]
        }

    // ── AF3-no-MSA stream ─────────────────────────────────────────────
    AF3_SETUP_DB()

    af3_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],
            rec[4],
            rec[5],
            file(rec[3]),
            file(rec[2]),
        )
    }

    AF3_NOMSA_ON_SURVIVORS(
        af3_input_ch,
        AF3_SETUP_DB.out.db_dir,
        AF3_SETUP_DB.out.ready_flag,
    )
    AF3_PARSE_OUTPUT(
        AF3_NOMSA_ON_SURVIVORS.out.prediction,
        Channel.value(file("${projectDir}/bin/parse_af3_output.py")),
    )

    // ── Biophysical stream ────────────────────────────────────────────
    biophysical_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],
            file(rec[1]),
            file(rec[2]),
        )
    }
    NEGSTEER_BIOPHYSICAL_METRICS(
        biophysical_input_ch,
        Channel.value(file("${projectDir}/bin/run_biophysical_metrics.py"))
    )

    // ── Rosetta stream ────────────────────────────────────────────────
    rosetta_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],
            file(rec[1]),
        )
    }
    fastrelax_xml_ch = Channel.value(
        file("${projectDir}/bin/fastrelax_for_ia.xml")
    )
    NEGSTEER_ROSETTA_METRICS(
        rosetta_input_ch,
        fastrelax_xml_ch,
        Channel.value(file("${projectDir}/bin/run_rosetta_metrics.py"))
    )

    // ── Merge the three streams into the final survivors CSV ──────────
    NEGSTEER_ORTHOGONAL_METRICS(
        extended_csv_ch,
        AF3_PARSE_OUTPUT.out.summary_csv.collect().ifEmpty([]),
        NEGSTEER_BIOPHYSICAL_METRICS.out.summary_csv.collect(),
        NEGSTEER_ROSETTA_METRICS.out.summary_csv.collect(),
        Channel.value(file("${projectDir}/bin/merge_orthogonal_metrics.py"))
    )

    ORTHOG_PLOTS(
        NEGSTEER_ORTHOGONAL_METRICS.out.final_csv,
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        Channel.value(file("${projectDir}/bin/orthogonal_metrics_plots.py")),
    )

    // Terminal outputs:
    //   ${params.outdir}/negative_steering/cross_sequence_summary.csv
    //   ${params.outdir}/negative_steering/cross_sequence_summary_with_interface_metrics.csv
    //   ${params.outdir}/orthogonal_metrics/survivors_with_orthogonal_metrics.csv
    //   ${params.outdir}/plots/orthogonal_*.png
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    def control_warnings = []
    if (params.run_negative_controls) {
        def summary_csv = new File(
            "${params.outdir}/negative_steering/cross_sequence_summary.csv"
        )
        if (summary_csv.exists()) {
            def lines = summary_csv.readLines()
            if (lines.size() >= 2) {
                def header = lines[0].split(',') as List
                def i_row_type = header.indexOf('row_type')
                def i_ipsae    = header.indexOf('rep_ipsae_min_median')
                def i_ra_eff   = header.indexOf('rep_ra_eff_vs_truth_median')
                def i_seq      = header.indexOf('mpnn_sequence')
                if (i_row_type < 0) {
                    log.warn "control miscalibration check skipped: no row_type column in ${summary_csv}"
                } else {
                    lines.drop(1).each { ln ->
                        def cells = ln.split(',', -1) as List
                        if (cells.size() <= i_row_type) return
                        def rt = cells[i_row_type]
                        if (!rt?.startsWith('control_')) return
                        def seq = (i_seq >= 0 && cells.size() > i_seq) ? cells[i_seq] : '?'
                        def ipsae_str  = (i_ipsae  >= 0 && cells.size() > i_ipsae)  ? cells[i_ipsae]  : ''
                        def ra_eff_str = (i_ra_eff >= 0 && cells.size() > i_ra_eff) ? cells[i_ra_eff] : ''
                        def hits = []
                        if (ipsae_str) {
                            try {
                                def v = ipsae_str as Double
                                if (v > params.controls_warning_ipsae_max) {
                                    hits << "ipSAE_min=${ipsae_str} > ${params.controls_warning_ipsae_max}"
                                }
                            } catch (NumberFormatException e) { /* ignore */ }
                        }
                        if (ra_eff_str) {
                            try {
                                def v = ra_eff_str as Double
                                if (v < params.controls_warning_ra_eff_min) {
                                    hits << "ra_eff=${ra_eff_str} < ${params.controls_warning_ra_eff_min}"
                                }
                            } catch (NumberFormatException e) { /* ignore */ }
                        }
                        if (hits) {
                            control_warnings << "  ${seq} (${rt}): ${hits.join('; ')}"
                        }
                    }
                }
            }
        } else {
            log.warn "control miscalibration check skipped: ${summary_csv} not found"
        }
    }

    def banner = """
    =============================================================
    Resample-MPNN Pipeline Complete
    =============================================================
    Project:           ${params.project_name}
    Output:            ${params.outdir}
    Source design:     ${params.source_design_id} (${params.source_run_dir})
    Duration:          ${workflow.duration}
    Success:           ${workflow.success}

    Terminal outputs:
      negative_steering/cross_sequence_summary.csv
      negative_steering/cross_sequence_summary_with_interface_metrics.csv
      orthogonal_metrics/survivors_with_orthogonal_metrics.csv
    =============================================================
    """.stripIndent()
    log.info banner

    if (control_warnings) {
        log.warn """
        =============================================================
        ⚠  RANKER MISCALIBRATION WARNING (Task 6 controls)
        =============================================================
        ${control_warnings.size()} control row(s) scored as good binders.
${control_warnings.join('\n')}
        =============================================================
        """.stripIndent()
    } else if (params.run_negative_controls) {
        log.info "Negative controls: no miscalibration warnings (all controls correctly scored as non-binders)."
    }
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
