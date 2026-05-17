#!/usr/bin/env nextflow

/*
 * =============================================================================
 * Receptor Resurfacing Pipeline (v0.3.0)
 * =============================================================================
 *
 * [HADDOCK3 Docking] → RFDiffusion →
 *   Rosetta Physics Filter → ProteinMPNN → Boltz2
 *
 * Inputs are PDB only:
 *   - Two separate PDBs (receptor + effector): docked with HADDOCK3
 *   - One pre-docked complex PDB: skips HADDOCK and proceeds directly
 *
 * Sequences, effector length, and hotspots are derived automatically where
 * possible, minimising manual parameter entry.
 *
 * Author: Josh Williams
 * Institute: John Innes Centre / The Sainsbury Laboratory
 * =============================================================================
 */

nextflow.enable.dsl = 2

import groovy.json.JsonOutput

// ---------------------------------------------------------------------------
// Parameter defaults
// ---------------------------------------------------------------------------

// ── Project ─────────────────────────────────────────────────────────────
params.project_name       = "receptor_resurfacing"
params.outdir             = "${launchDir}/${params.project_name}_results"

// Normalise params.outdir to an absolute path.  When the user supplies
// a relative outdir via -params-file (e.g. "./results"), Nextflow stores
// the literal string and downstream "${params.outdir}/foo" interpolation
// produces a relative path.  Relative paths stamped into shared CSVs
// (e.g. rep_canonical_pdb via cross_sequence_summary.py's
// --published-runs-dir) then resolve against the per-task container CWD
// in downstream processes, not the launch dir — so files at the right
// host path are reported missing.  Forcing absolute resolution here
// once means every "${params.outdir}/..." use is safe regardless of
// what the user supplied.
if (!new File(params.outdir.toString()).isAbsolute()) {
    params.outdir = file(params.outdir).toAbsolutePath().normalize().toString()
}

// ── Inputs ──────────────────────────────────────────────────────────────
// Provide EITHER a complex PDB (receptor+effector already docked)
// OR separate receptor/effector files (each can be .pdb or .fasta).
// Both inputs must be PDB.  Provide either a single pre-docked complex
// (params.pdb_file) or two separate PDBs (params.receptor_input +
// params.effector_input), in which case HADDOCK3 will dock them.

params.pdb_file           = null   // Pre-docked complex PDB
params.receptor_input     = null   // Receptor PDB
params.effector_input     = null   // Effector PDB

params.receptor_chain     = "A"
params.effector_chain     = "B"

// ── Post-RFDiffusion chain convention ────────────────────────────────
// RFDiffusion's split function (rfdiffusion_filter.py:write_split_pdb)
// HARDCODES the output chain labels to A=receptor, B=effector,
// regardless of what the input PDB or contigs called them.  Every
// process that consumes the split design PDBs (post-RFDiffusion onward)
// must therefore use these fixed labels, NOT params.receptor_chain /
// params.effector_chain (which describe the INPUT PDB convention).
// Exposed as parameters for visibility, but in practice always A/B.
params.rfdiff_output_receptor_chain = "A"
params.rfdiff_output_effector_chain = "B"

// ── Sequences (auto-derived if not provided) ────────────────────────────
params.receptor_seq       = null
params.effector_seq       = null

// ── RFDiffusion ─────────────────────────────────────────────────────────
params.contigs            = null   // RFDiffusion contig string
params.hotspot            = ""     // Auto-derived from HADDOCK if blank
params.num_designs        = 10
params.rfdiff_iterations  = 50
params.rfdiff_checkpoint  = "Complex_beta_ckpt.pt" // Filename inside /opt/RFdiffusion/models/
params.min_hotspot_frac   = 0.0     // Min fraction of contacts inside design region (0.0 = no filtering)
params.symmetry                = "none"
params.order                   = 1
params.add_potential           = false       // Enable guiding potentials during diffusion
params.rfdiff_guide_scale      = 2           // Global potential multiplier
params.rfdiff_guide_decay      = "quadratic" // Weight decay schedule: constant | linear | quadratic | cubic
params.rfdiff_interface_weight = 1.0         // interface_ncontacts weight (0 to disable)
params.rfdiff_rog_weight       = 0.5         // monomer_ROG weight (0 to disable)
params.rfdiff_rog_min_dist     = 5           // monomer_ROG minimum Rg floor (Å)

// ── Rosetta pre-validation ──────────────────────────────────────────────
params.sc_threshold       = 0.62   // Min shape complementarity to pass (Overath et al., 2025)
params.stop_after_rosetta = false  // Stop pipeline after Rosetta Sc filtering; resume with -resume

// ── ProteinMPNN ─────────────────────────────────────────────────────────
params.num_seqs           = 8
params.mpnn_sampling_temp = 0.25
params.rm_aa              = "C"
params.mpnn_top_n         = 0      // Select top N sequences by MPNN score for Boltz2 (0 = use all)

// ── Boltz2 Negative Steering ────────────────────────────────────────────
// The Boltz2 validation stage runs the negative-steering pipeline per
// MPNN sequence (plan → predict-one×N → collect → reversion pass →
// finalize → postprocess).  Each MPNN sequence becomes one GPU SLURM
// job; concurrency across sequences is capped by params.max_boltz2_parallel.
// See modules/negative_steering.nf for the process shape.
//
// Each sequence's parent RFDiffusion design PDB serves as its ground
// truth.  `plan` is driven sequence-only (receptor/effector FASTAs
// split from the MPNN output) against that ground truth, using
// pre-computed protected-set files (true interface + design region)
// derived from rfdiffusion_metrics.json by NEGSTEER_DERIVE_INDICES.
//
// All knobs below are forwarded to boltz2_negative_steering.py's
// `plan` stage.  The ones you are most likely to change are
// negsteer_mode, negsteer_max_mutations, negsteer_n_designs and
// negsteer_num_seeds.

// ── Steering strategy ──────────────────────────────────────────────────
params.negsteer_mode                 = "mild"     // strong | mild | conservative | alanine
                                                  // mild is the default for resurfacing designs:
                                                  // D/E/K/R substitutions only, which are
                                                  // aggressive enough to force binding-site
                                                  // migration without introducing volumetrically
                                                  // incompatible side chains.
params.negsteer_max_mutations        = 6          // Number of receptor residues to mutate per
                                                  // design.  If strict protection reduces the
                                                  // candidate pool below this, the pipeline
                                                  // transparently drops to the pool size (see
                                                  // max_mutations_effective in plan.json).
params.negsteer_candidate_pool_size  = 10         // Size of the pool of residues sampled from.
                                                  // Larger than max_mutations => designs sample
                                                  // different subsets across runs, increasing
                                                  // positional diversity.
params.negsteer_protected_set_source = "design_region_union"
                                                  // true_interface | design_region_union
                                                  // design_region_union (default) excludes
                                                  // both the true interface and the RFDiffusion
                                                  // design region from the candidate pool.
                                                  // Required for resurfacing designs, where
                                                  // mutations must never land inside the
                                                  // designed block (see notes5.md / notes6.md).

// ── Design / seed counts ───────────────────────────────────────────────
params.negsteer_n_designs            = 20         // Steered designs generated per MPNN sequence.
params.negsteer_num_seeds            = 3          // Independent Boltz seeds per unique sequence.
                                                  // Total predictions per MPNN sequence =
                                                  // n_designs × num_seeds (≤ 60 by default).
params.negsteer_n_cycles             = 1          // Number of steering cycles.  The NextFlow
                                                  // integration currently runs single-cycle
                                                  // (N=1) only; multi-cycle will be added when
                                                  // the negative-steering codebase is
                                                  // restructured (see notes6 open item 4).

// ── Boltz prediction hyperparameters ───────────────────────────────────
params.negsteer_diffusion_samples    = 5          // Diffusion samples per seed (Boltz-side parallel).
params.negsteer_recycling_steps      = 3          // Recycling steps.
params.negsteer_rmsd_threshold       = 6.0        // Initial-RMSD threshold under which steering
                                                  // is skipped (plan writes a one-line summary).
params.negsteer_contact_cutoff       = 4.5        // Heavy-atom contact cutoff (Å) for interface
                                                  // detection in plan.
params.negsteer_no_kernels           = false      // Pass --no-kernels to boltz.  Default false:
                                                  // the boltz2_negsteer.img container includes
                                                  // the pinned cuEquivariance cu12 stack
                                                  // (cuequivariance==0.9.1), giving Boltz the
                                                  // CUDA triangle-attention kernels for ~20-30%
                                                  // faster inference vs pure PyTorch.  Set true
                                                  // as a fast debugging fallback.

// ── Post-processing gate ───────────────────────────────────────────────
// Forwarded to compute-final-metrics.  Rows with
// steered_ra_eff_vs_truth >= this threshold don't get confidence
// metrics populated (saves time on designs that are already
// structurally wrong).
params.negsteer_postprocess_rmsd_threshold = 5.0
params.negsteer_postprocess_metric_column  = "steered_ra_eff_vs_truth"
params.negsteer_postprocess_contact_cutoff = 5.0

// ── Negative controls (Task 6, P0-6) ───────────────────────────────
// Per-design negative controls — one scrambled + one polyA receptor
// per parent RFDiffusion design — run through the same negative-
// steering inner loop as a real MPNN sequence.  See
// modules/negsteer_controls.nf for what each control does and
// bin/build_control_sequences.py for the construction logic.
params.run_negative_controls         = true   // Default ON for production;
                                              // tests flip to false to keep
                                              // wall-clock low.
params.negsteer_controls_n_designs   = 1      // Steered designs per control.
                                              // Production runs 20 per real
                                              // MPNN sequence; controls just
                                              // need a diagnostic reading,
                                              // so 1 is enough and 20× cheaper.
// Miscalibration warning thresholds, applied to control rows in
// cross_sequence_summary.csv at workflow.onComplete time.  A control
// that the ranker thinks is a good binder (low ra_eff OR high ipSAE)
// indicates the ranker is sequence-blind on this target — emit a
// warning, do NOT fail the run.  Per todo_list3 P0-6 acceptance.
params.controls_warning_ipsae_max    = 0.5    // ipSAE > this on a control
                                              // → suspicious
params.controls_warning_ra_eff_min   = 5.0    // ra_eff < this on a control
                                              // → suspicious

// ── Interface-restricted metrics (P0-29) ────────────────────────────
// Parameters for NEGSTEER_INTERFACE_METRICS, which adds iRMSD / fnat /
// DockQ / ipsae_15 / intact_core columns to cross_sequence_summary.csv.
params.interface_plddt_trim_threshold = 50.0
params.interface_intact_threshold     = 5.0

// ── Weighted Jaccard contact overlap (P0-38) ────────────────────────
// Forwarded to compute_interface_metrics.py.  Hard distance ceiling
// (Å) for receptor–effector residue pairs entering the weighted_jaccard
// sum.  At d=8 the Gaussian weight exp(-(d-4)^2/2.25) is ≈ 8e-4, so the
// cutoff is numerically inert — only there to bound the all-vs-all loop.
// The Gaussian's μ=4 and σ²=1.125 are deliberately NOT exposed as flags;
// re-tuning them would silently invalidate cross-cohort weighted_jaccard
// comparisons in Task 8 AUROC.  See compute_interface_metrics.py.
params.weighted_jaccard_pair_cutoff   = 8.0

// ── Orthogonal validation metrics (P0-31) ───────────────────────────
// AF3-no-MSA, FreeSASA/MDAnalysis biophysical, Rosetta Sc/ΔΔG.  All
// run on survivors of the filter cascade that emerges from P0-29.
//
// AF3 seed list: matches Boltz's num_seeds=3 (the negative-steering
// default).  AF3 does not expose a CLI for seed count on the NBI
// build — seeds come via the modelSeeds JSON field.
params.af3_nomsa_seeds = [42, 123, 456]

// Filter cascade thresholds — first-pass, to be retuned after Task 8
// AUROC (three-target cohort).  Per todo_list3 P0-31: below-threshold
// survivors are FLAGGED in orthogonal_flags, not dropped — except
// af3_nomsa agreement failure, which is a hard drop.
params.orthogonal_filter_sc_min     = 0.55   // Lawrence-Colman Sc
params.orthogonal_filter_bsa_min    = 600    // Å² (Overath et al.)
params.orthogonal_filter_plddt_min  = 0.75   // mean interface pLDDT
params.orthogonal_filter_af3_ra_max = 5.0    // ra_eff threshold, hard drop

// ── Quality control ─────────────────────────────────────────────────────
params.max_poly_x         = 5
params.min_pct_identity   = 0.0
params.max_pct_identity   = 100.0

// ── HADDOCK3 ────────────────────────────────────────────────────────────
params.haddock_sampling   = 10000  // Rigid-body sampling (10000=semi-blind)
params.haddock_seletop    = 400    // Top N rigid-body models passed to flexref
params.rfdiff_contact_cutoff = 8.0 // Cα–Cα cutoff (Å) for rfdiffusion_filter contact
                                   // detection on RFDiffusion-designed complexes.
// HADDOCK restraints (per Session 7 grill-me, notes/design_audit.md Q137).
// At least one of haddock_contact_pairs or haddock_receptor_active_residues
// must be non-empty in Branch A; validate_params.py enforces this.
params.haddock_contact_pairs            = ""    // Hard CA-CA pin pairs, e.g. "A25-C42 A13-C94"
params.haddock_receptor_active_residues = ""    // Soft AIR receptor side, e.g. "25,35,40-44"
params.haddock_effector_active_residues = ""    // Soft AIR effector side; empty = entire chain
params.haddock_pair_distance            = "2,2,4"  // Global pair distance "target,lo_dev,hi_dev"
params.haddock_chosen_cluster           = null  // Set to a cluster_id to override auto-pick
params.stop_after_haddock               = false // Halt after HADDOCK_PLOTS for manual inspection

// ── Infrastructure ──────────────────────────────────────────────────────
// All container paths (rfdiff_container, rosetta_container, boltz2_container,
// colabfold_container) are declared in nextflow.config — single source of
// truth across the main pipeline and per-module test scripts.
//
// Note: the boltz2_container image (boltz2_negsteer.img) provides the
// Boltz-2 CLI plus pinned cuEquivariance CUDA kernels, numpy, BioPython,
// and the orthogonal-metrics analysis packages (DockQ, freesasa,
// MDAnalysis).  The negative-steering pipeline runs every phase inside
// this container.

params.max_boltz2_parallel = 5

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------

if (!params.pdb_file && !params.receptor_input) {
    error "Provide either --pdb_file (complex) or --receptor_input + --effector_input"
}
if (params.receptor_input && !params.effector_input) {
    error "If --receptor_input is provided, --effector_input is also required"
}
if (!params.contigs) {
    error "Please provide --contigs (RFDiffusion contig string)"
}

// ---------------------------------------------------------------------------
// Include modules
// ---------------------------------------------------------------------------

include { EXTRACT_SEQUENCES                    } from './modules/preprocessing'
include { RESOLVE_CONTIGS                      } from './modules/preprocessing'
include { WRITE_DUMMY_MAPPING                  } from './modules/preprocessing'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_REC } from './modules/preprocessing'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_EFF } from './modules/preprocessing'

include { HADDOCK3_PREPARE                     } from './modules/haddock'
include { HADDOCK3_DOCK                        } from './modules/haddock'
include { HADDOCK3_PLOTS                       } from './modules/haddock'
include { HADDOCK_CLUSTER_METRICS              } from './modules/haddock'
include { SELECT_HADDOCK_CLUSTER               } from './modules/haddock'
include { BUILD_CONTIGS                        } from './modules/haddock'

include { RFDIFFUSION                          } from './modules/rfdiffusion'
include { RFDIFFUSION_FILTER                   } from './modules/rfdiffusion'
include { RFDIFFUSION_PLOTS                    } from './modules/rfdiffusion'
include { ROSETTA_SC                           } from './modules/rosetta_filtering'
include { ROSETTA_FILTER                       } from './modules/rosetta_filtering'
include { ROSETTA_FILTER_PLOTS                 } from './modules/rosetta_filtering'
include { MPNN_FIXED_POSITIONS                 } from './modules/proteinmpnn'
include { PROTEINMPNN                          } from './modules/proteinmpnn'
include { SEQUENCE_CORRECTION                  } from './modules/proteinmpnn'
include { SEQUENCE_QC                          } from './modules/proteinmpnn'
include { MPNN_DESIGN_REGION_SCORE             } from './modules/proteinmpnn'
include { MPNN_CLUSTER                         } from './modules/proteinmpnn'
include { MPNN_SELECT_TOP                      } from './modules/proteinmpnn'
include { MPNN_PLOTS                           } from './modules/proteinmpnn'
include { NEGSTEER_DERIVE_INDICES               } from './modules/negative_steering'
include { NEGSTEER_RUN_ONE                      } from './modules/negative_steering'
include { NEGSTEER_CROSS_SEQUENCE               } from './modules/negative_steering'
include { NEGSTEER_PLOTS                        } from './modules/negative_steering'
include { NEGSTEER_WITHIN_SEQUENCE_PLOTS        } from './modules/negative_steering'
include { NEGSTEER_CONTROLS                     } from './modules/negsteer_controls'
include { DERIVE_INPUT_INDICES                  } from './modules/negsteer_controls'
include { NEGSTEER_INTERFACE_METRICS            } from './modules/negsteer_interface_metrics'
include { EXTRACT_SURVIVOR_MANIFEST             } from './modules/negsteer_manifest'
include { AF3_SETUP_DB                          } from './modules/negsteer_af3_nomsa'
include { AF3_NOMSA_ON_SURVIVORS                } from './modules/negsteer_af3_nomsa'
include { AF3_PARSE_OUTPUT                      } from './modules/negsteer_af3_nomsa'
include { NEGSTEER_BIOPHYSICAL_METRICS          } from './modules/negsteer_biophysical_metrics'
include { NEGSTEER_ROSETTA_METRICS              } from './modules/negsteer_rosetta_metrics'
include { NEGSTEER_ORTHOGONAL_METRICS           } from './modules/negsteer_orthogonal_metrics'
include { ORTHOG_PLOTS                          } from './modules/negsteer_orthogonal_metrics'


// ---------------------------------------------------------------------------
// Helper: detect if file is FASTA or PDB
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    // ── Parameter validation ─────────────────────────────────────────────
    // Fail fast before any compute is dispatched.  bin/validate_params.py
    // applies all PARAM_SPECS (range/choice/regex/custom) and prints every
    // error before exiting non-zero; if it returns non-zero we abort the
    // entire workflow.  Coverage report: `python3 bin/validate_params.py --report`.
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
    // BRANCH A: Separate receptor + effector PDBs → HADDOCK3 docking
    // =====================================================================
    if (params.receptor_input) {

        rec_file = Channel.fromPath(params.receptor_input, checkIfExists: true)
        eff_file = Channel.fromPath(params.effector_input, checkIfExists: true)

        // ── Receptor: PDB used directly, dummy mapping for BUILD_CONTIGS ──
        receptor_pdb_ch = rec_file
        // Write an empty dummy mapping (exactly 2 bytes: "{}") so
        // BUILD_CONTIGS gets a real file.  load_mapping() in build_contigs.py
        // treats files of <=2 bytes as "no mapping" and returns None.
        WRITE_DUMMY_MAPPING_REC(Channel.value("receptor"))
        rec_trim_mapping_ch = WRITE_DUMMY_MAPPING_REC.out.mapping

        // ── Effector: PDB used directly, dummy mapping for BUILD_CONTIGS ──
        effector_pdb_ch = eff_file
        WRITE_DUMMY_MAPPING_EFF(Channel.value("effector"))
        eff_trim_mapping_ch = WRITE_DUMMY_MAPPING_EFF.out.mapping

        // ── HADDOCK3 docking ────────────────────────────────────────────
        // Per Session 7 (notes/design_audit.md A131): the contig string
        // is NOT a HADDOCK restraint source.  Restraints come from the
        // user's contact_pairs + receptor/effector active residue params.
        HADDOCK3_PREPARE(
            receptor_pdb_ch,
            effector_pdb_ch,
            params.receptor_chain,
            params.effector_chain,
            params.haddock_contact_pairs,
            params.haddock_receptor_active_residues,
            params.haddock_effector_active_residues,
            params.haddock_pair_distance,
            params.contigs,
            Channel.value(file("${projectDir}/bin/haddock3_prepare.py"))
        )

        HADDOCK3_DOCK(
            HADDOCK3_PREPARE.out.receptor_pdb_out,
            HADDOCK3_PREPARE.out.effector_pdb_out,
            HADDOCK3_PREPARE.out.ambig_restraints,
            HADDOCK3_PREPARE.out.unambig_restraints,
            params.haddock_sampling,
            params.haddock_seletop,
            Channel.value(file("${projectDir}/bin/collect_haddock3_dock.py"))
        )

        HADDOCK3_PLOTS(
            HADDOCK3_DOCK.out.capri_scores,
            HADDOCK3_DOCK.out.cluster_summary,
            HADDOCK3_DOCK.out.run_dir,
            params.haddock_receptor_active_residues,
            params.haddock_effector_active_residues,
            Channel.value(file("${projectDir}/bin/haddock3_plots.py"))
        )

        // ── Per-cluster metrics (BSA, Sc, COM, AIR/pair satisfaction,
        //    clash counts) — computed once and consumed by SELECT below.
        HADDOCK_CLUSTER_METRICS(
            HADDOCK3_DOCK.out.haddock_report,
            HADDOCK3_DOCK.out.cluster_models,
            HADDOCK3_PREPARE.out.restraints_summary,
            HADDOCK3_PREPARE.out.ambig_restraints,
            HADDOCK3_PREPARE.out.unambig_restraints,
            Channel.value(file("${projectDir}/bin/haddock_cluster_metrics.py")),
            Channel.value(file("${projectDir}/bin"))
        )

        // ── Stop-and-resume gate (per Session 7 A134) ───────────────────
        // When params.stop_after_haddock is true, halt before SELECT so
        // the user can inspect plots + per-cluster best_cluster*.pdb files
        // (already published to ${outdir}/haddock/) and pick a cluster.
        // To resume after picking: re-run with
        //     params.haddock_chosen_cluster: N
        // (Nextflow -resume reuses cached HADDOCK_DOCK + HADDOCK_CLUSTER_METRICS
        //  outputs; only SELECT and downstream re-execute.)
        if (params.stop_after_haddock) {
            log.warn(
                "stop_after_haddock=true → pipeline will halt after HADDOCK_PLOTS.\n" +
                "  Inspect ${params.outdir}/haddock/best_cluster*.pdb and\n" +
                "  ${params.outdir}/haddock/cluster_metrics.json, then resume with\n" +
                "  params.haddock_chosen_cluster: <id> in your params file and\n" +
                "  re-run with --resume.\n" +
                "  (Cleared by setting stop_after_haddock=false on resume.)"
            )
            return
        }

        // ── Select chosen cluster (auto or user-specified) ───────────────
        SELECT_HADDOCK_CLUSTER(
            HADDOCK3_DOCK.out.haddock_report,
            HADDOCK_CLUSTER_METRICS.out.cluster_metrics,
            HADDOCK3_PREPARE.out.restraints_summary,
            HADDOCK3_DOCK.out.cluster_models,
            receptor_pdb_ch,
            effector_pdb_ch,
            params.haddock_contact_pairs,
            params.haddock_receptor_active_residues,
            params.haddock_effector_active_residues,
            params.haddock_pair_distance,
            params.haddock_chosen_cluster,
            Channel.value(file("${projectDir}/bin/select_haddock_cluster.py"))
        )

        // ── Build updated contigs & extract sequences ───────────────────
        BUILD_CONTIGS(
            SELECT_HADDOCK_CLUSTER.out.selected_pdb,
            params.receptor_chain,
            params.effector_chain,
            params.contigs,
            rec_trim_mapping_ch,
            eff_trim_mapping_ch,
            Channel.value(file("${projectDir}/bin/build_contigs.py"))
        )

        // Hotspot for RFDiffusion: user-only.  Auto-derivation from the
        // docked complex was removed in Session 7 (A143/A145) because the
        // contact pattern HADDOCK chose is already visible to RFDiffusion
        // from the input PDB, so feeding it back as a hotspot bias just
        // reinforces that pattern and reduces design diversity.  Set
        // params.hotspot only when you have biological knowledge of
        // target residues you want contacted (e.g. homologous binding sites).
        hotspot_ch = Channel.value(params.hotspot ?: "")

        // Channels for downstream steps
        rfdiff_pdb_ch   = BUILD_CONTIGS.out.rfdiffusion_pdb
        contigs_ch      = BUILD_CONTIGS.out.contigs_txt.map { it.text.trim() }.first()
        updated_params_ch = BUILD_CONTIGS.out.updated_params

    }
    // =====================================================================
    // BRANCH B: Pre-docked complex PDB (direct mode)
    // =====================================================================
    else {
        rfdiff_pdb_ch = Channel.fromPath(params.pdb_file, checkIfExists: true)
        hotspot_ch    = Channel.value(params.hotspot ?: "")

        // Resolve the raw contig string to PDB-resolved coordinates.
        // Branch A does this inside BUILD_CONTIGS (which also handles
        // HADDOCK-driven coordinate remapping); Branch B has no HADDOCK
        // step, so we just call contig_utils.resolve_contigs() via the
        // RESOLVE_CONTIGS process.  Downstream consumers (including
        // pipeline_correct_sequences.py) can then assume contigs_ch is
        // always PDB-resolved regardless of which branch produced it.
        RESOLVE_CONTIGS(
            rfdiff_pdb_ch,
            params.contigs,
            Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
        )
        contigs_ch = RESOLVE_CONTIGS.out.resolved_contigs
            .map { it.text.trim() }
            .first()

        // Extract sequences (incl. receptor_start_pdb) from the input PDB.
        EXTRACT_SEQUENCES(
            rfdiff_pdb_ch,
            params.receptor_chain,
            params.effector_chain
        )
        updated_params_ch = EXTRACT_SEQUENCES.out.sequences_json
    }

    // =====================================================================
    // Resolve sequences and receptor_start_pdb from the JSON producer
    // =====================================================================

    // Parse receptor_seq, effector_seq, and receptor_start_pdb from the
    // JSON written by BUILD_CONTIGS (Branch A) or EXTRACT_SEQUENCES
    // (Branch B).  Both producers write the same key set.
    //
    // We need these as VALUE channels because MPNN_FIXED_POSITIONS fans
    // out across all RFDiffusion designs and must be able to read each
    // value repeatedly.  .map() on a queue channel produces a queue
    // channel that is consumed after the first emission, so we wrap the
    // shared parse in .first() to convert queue → value before splitting.
    seqs_ch = updated_params_ch.map { json_file ->
        def slurp = new groovy.json.JsonSlurper()
        def data  = slurp.parseText(json_file.text)
        return [data.receptor_seq, data.effector_seq, data.receptor_start_pdb ?: 1]
    }.first()
    receptor_seq_ch       = seqs_ch.map { it[0] }
    effector_seq_ch       = seqs_ch.map { it[1] }
    receptor_start_pdb_ch = seqs_ch.map { it[2] }

    // =====================================================================
    // Step 1: RFDiffusion
    // =====================================================================
    RFDIFFUSION(
        rfdiff_pdb_ch,
        contigs_ch,
        hotspot_ch,
        params.num_designs,
        params.rfdiff_iterations,
        params.rfdiff_checkpoint,
        params.add_potential,
        params.rfdiff_guide_scale,
        params.rfdiff_guide_decay,
        params.rfdiff_interface_weight,
        params.rfdiff_rog_weight,
        params.rfdiff_rog_min_dist,
        Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
    )

    RFDIFFUSION_FILTER(
        RFDIFFUSION.out.design_pdbs,
        rfdiff_pdb_ch,
        contigs_ch,
        hotspot_ch,
        params.receptor_chain,
        params.effector_chain,
        params.rfdiff_contact_cutoff,
        params.min_hotspot_frac,
        Channel.value(file("${projectDir}/bin/rfdiffusion_filter.py"))
    )

    design_pdbs_ch = RFDIFFUSION_FILTER.out.passing_pdbs.flatten()

    RFDIFFUSION_PLOTS(
        RFDIFFUSION_FILTER.out.metrics,
        Channel.value(file("${projectDir}/bin/rfdiffusion_plots.py"))
    )

    // =====================================================================
    // Step 2: Rosetta Pre-Validation Physics Filter
    // =====================================================================
    // Run InterfaceAnalyzer on each split two-chain PDB from RFDiffusion
    // to compute shape complementarity (Sc), dG_separated, and dSASA_int.
    // Filter out designs with Sc below threshold before ProteinMPNN.
    //
    // Uses split PDBs (chain A = receptor, chain B = effector) so that
    // InterfaceAnalyzer can identify the interface correctly.

    split_pdbs_ch = RFDIFFUSION_FILTER.out.split_pdbs.flatten()

    ROSETTA_SC(
        split_pdbs_ch
    )

    ROSETTA_FILTER(
        ROSETTA_SC.out.pdb_and_scores.collect(),
        params.sc_threshold,
        Channel.value(file("${projectDir}/bin/rosetta_filter_collect.py"))
    )

    rosetta_passing_pdbs_ch = ROSETTA_FILTER.out.passing_pdbs.flatten()

    ROSETTA_FILTER_PLOTS(
        ROSETTA_FILTER.out.metrics,
        Channel.value(file("${projectDir}/bin/rosetta_filter_plots.py"))
    )

    if (!params.stop_after_rosetta) {

    // =====================================================================
    // Step 3: ProteinMPNN
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
        params.num_designs,
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
    // NOTE: the old `boltz2_metadata_ch` (MPNN_*.out.*_metadata) is not
    // wired into the negative-steering stage. The metadata is still
    // published to ${params.outdir}/sequences/ by the upstream MPNN
    // processes; a future aggregation step will re-join it with the
    // cross-sequence triage table by (design, seq) key.

    // =====================================================================
    // Step 4: Boltz2 Negative-Steering Validation (per MPNN sequence)
    // =====================================================================
    //
    // Shape of the per-sequence input tuple (one per MPNN FASTA):
    //   [seq_name, design_stem, mpnn_fasta, design_pdb,
    //    design_region_file, true_interface_file]
    //
    // Built via three keyed channels joined on design_stem:
    //   (a) FASTA stream:   [design_stem, seq_name, mpnn_fasta]
    //   (b) Design-PDB map: [design_stem, design_pdb]
    //   (c) Indices map:    [design_stem, design_region, true_interface]
    //                       (produced by NEGSTEER_DERIVE_INDICES)
    //
    // Only designs that produced MPNN FASTAs get their indices
    // derived - if Rosetta filtered out a design, there is nothing
    // downstream to negative-steer against, so we skip the derivation
    // to save CPU.

    // -- 4a: FASTA stream [design_stem, seq_name, mpnn_fasta] --------
    fasta_stream_ch = mpnn_fastas_dir_ch
        .flatMap { fdir ->
            fdir.listFiles()
                .findAll { it.name.endsWith(".fasta") }
                .collect { f ->
                    // Expect MPNN's "design_<N>_seq_<M>.fasta" naming.
                    // Filenames that do not match are filtered out a
                    // step below so stray files do not poison the join.
                    def stem = f.name.replaceAll(/\.fasta$/, '')
                    def m = stem =~ /^(design_\d+)_seq_\d+$/
                    def design_stem = m ? m[0][1] : null
                    [design_stem, stem, f]
                }
        }
        .filter { trip -> trip[0] != null }

    // -- 4b: Design-PDB map [design_stem, design_pdb] ----------------
    // Keyed map across all RFDiffusion outputs.  Used both to feed
    // index derivation and to join into the per-sequence tuple.
    design_pdb_map_ch = RFDIFFUSION.out.design_pdbs
        .flatten()
        .map { pdb -> [pdb.baseName, pdb] }   // "design_3" -> file

    // -- 4c: Derive index files for every RFDiffusion design ---------
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
        RFDIFFUSION_FILTER.out.metrics.first(),   // value-channel broadcast
        Channel.value(file("${projectDir}/bin/derive_design_region.py")),
        Channel.value(file("${projectDir}/bin/derive_true_interface.py"))
    )

    // NEGSTEER_DERIVE_INDICES emits tuple(stem, design_region, true_iface)
    indices_map_ch = NEGSTEER_DERIVE_INDICES.out.indices

    // -- 4d: Three-way join into the per-sequence input tuple --------
    // Use combine(by: 0) instead of join() because join() is one-to-one
    // (consumes one entry from each side per match), whereas combine
    // (by: 0) is one-to-many keyed on the join key — so a single PDB
    // can pair with multiple FASTAs that share its design stem.
    run_one_input_ch = fasta_stream_ch
        .combine(design_pdb_map_ch, by: 0)
        .combine(indices_map_ch,    by: 0)
        .map { stem, seq_name, fasta, pdb, region, iface ->
            tuple(seq_name, stem, fasta, pdb, region, iface)
        }

    // -- 4e: Compose plan_extra_args string from params --------------
    // Passed verbatim to boltz2_negative_steering.py plan via the
    // bin/negative_steering_run_one.sh orchestrator.  Keep the flag
    // values surface-level here so users can see what knobs they are
    // paying for in every run.
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
                 "the NextFlow integration currently runs single-cycle only. " +
                 "Multi-cycle will be added when the negative-steering " +
                 "codebase is restructured (notes6 open item 4)."
    }

    // -- 4f: Fan out - one GPU SLURM job per MPNN sequence -----------
    // NEGSTEER_RUN_ONE consumes split design PDBs from RFDiffusion,
    // which always have receptor=A, effector=B (hardcoded by
    // rfdiffusion_filter.py:write_split_pdb).  Use the rfdiff_output_*
    // params here, not params.receptor_chain / params.effector_chain
    // (which describe the input PDB convention).
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

    // -- 4f-bis: Negative controls (Task 6 revision — v7) ------------
    // TWO negative controls per cohort run against the INPUT STRUCTURE
    // (the receptor+effector complex going into RFDiffusion), not
    // against RFDiffusion design outputs.  They answer "if RFDiffusion
    // did a terrible job, what metrics would those sequences get?" —
    // a property of the target + its contigs, not of any particular
    // design.  Per-design pathology is addressed later by Task 8 AUROC.
    //
    // Gated on params.run_negative_controls — flip to false in test
    // harnesses to keep wall-clock low.
    if (params.run_negative_controls) {

        // Compose the controls plan_extra_args — identical to the main
        // args except --n-designs is overridden to params.negsteer_
        // controls_n_designs (default 1; controls only need one
        // steered design to produce a diagnostic reading).
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
        // Both branches of the preprocessing workflow (A: HADDOCK-driven;
        // B: pre-complex) produce rfdiff_pdb_ch and contigs_ch by this
        // point in the graph, so DERIVE_INPUT_INDICES is branch-agnostic.
        DERIVE_INPUT_INDICES(
            rfdiff_pdb_ch,
            contigs_ch,
            params.receptor_chain,
            params.effector_chain,
            params.negsteer_contact_cutoff,
            Channel.value(file("${projectDir}/bin/derive_input_design_region.py"))
        )

        // Build a 2-element input channel for NEGSTEER_CONTROLS.  Each
        // element is (control_name, control_type, input_pdb,
        //             design_region, true_interface).
        // The .combine() here is a Cartesian product: 2 control types ×
        // 1 input-structure bundle = 2 process invocations.
        // .first() on each combine target converts the upstream single-
        // emission queue to a value channel — safe for operator re-use
        // and semantically identical (we're combining with exactly one
        // upstream element either way).
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

        // Mix steered + control workdirs into one channel feeding the
        // cross-sequence aggregator.  cross_sequence_summary.py picks
        // up the row_type.txt sidecar from each workdir and tags the
        // resulting row appropriately; control rows are excluded from
        // cross-rank scoring but appear in the output for diagnosis.
        all_workdirs_ch = NEGSTEER_RUN_ONE.out.per_sequence_workdir
            .mix(NEGSTEER_CONTROLS.out.per_control_workdir)

    } else {
        all_workdirs_ch = NEGSTEER_RUN_ONE.out.per_sequence_workdir
    }

    // -- 4g: Cross-sequence aggregation ------------------------------
    // Collect every per-sequence workdir (steered + controls if
    // enabled) and hand them to the tier-then-composite aggregator.
    // Single global CPU task.
    //
    // The collected workdir list is bound to a local channel because
    // it feeds two downstream consumers (NEGSTEER_CROSS_SEQUENCE here
    // and NEGSTEER_INTERFACE_METRICS / EXTRACT_SURVIVOR_MANIFEST in
    // step 5 below).  Nextflow caches `.collect()` results, but binding
    // it to a name also makes the dataflow easier to read.
    per_sequence_workdirs_ch = all_workdirs_ch.collect()

    NEGSTEER_CROSS_SEQUENCE(
        per_sequence_workdirs_ch,
        // Phase 4: route through the typed CLI entry point
        // (DesignCohort.emit_cross_summary_from_dirs).  Bit-identical
        // CSV output relative to cross_sequence_summary.py — the typed
        // wrapper delegates to the same aggregate() for column-level
        // construction.  See bin/cross_summary_v2.py for the routing.
        //
        // Pass the entire bin/ directory (not just cross_summary_v2.py)
        // so Nextflow content-hashes every transitively-imported file.
        // Closes the indirect-import cache-busting gap.
        Channel.value(file("${projectDir}/bin")),
        MPNN_DESIGN_REGION_SCORE.out.scored_metadata
    )

    // ── Step 4b: Diagnostic plots from cross_summary + per-sequence
    //            workdirs.  Two parallel processes — cohort-level and
    //            within-sequence (per-seed dispersion).  Both consume
    //            the same upstream channels; both run in parallel
    //            with NEGSTEER_INTERFACE_METRICS below (no downstream
    //            consumer of plot artifacts).
    //
    //            Cohort plots also need input_design_region.txt (for
    //            mutation-impact shading) — emitted by
    //            DERIVE_INPUT_INDICES in step 3.
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
    // Step 5: Post-negative-steering orthogonal metrics (P0-29 + P0-31)
    // =====================================================================
    // Channel topology mirrors tests/orthogonal_metrics/test_orthogonal_metrics.nf
    // exactly — see that file for the cached-data smoke test that
    // validates this graph end-to-end on design_0..design_3 inputs.
    //
    //   NEGSTEER_CROSS_SEQUENCE
    //      │
    //      ▼
    //   NEGSTEER_INTERFACE_METRICS               (P0-29)
    //      │
    //      ▼
    //   EXTRACT_SURVIVOR_MANIFEST                (per-survivor tuples)
    //      │
    //      ├──→  AF3_SETUP_DB ──→ AF3_NOMSA_ON_SURVIVORS ──→ AF3_PARSE_OUTPUT
    //      ├──→  NEGSTEER_BIOPHYSICAL_METRICS                            (fan-out)
    //      └──→  NEGSTEER_ROSETTA_METRICS                                (fan-out)
    //      │
    //      ▼
    //   NEGSTEER_ORTHOGONAL_METRICS              (merge + filter cascade)

    // ── P0-29: interface-restricted metrics ───────────────────────────
    NEGSTEER_INTERFACE_METRICS(
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/compute_interface_metrics.py"))
    )
    extended_csv_ch = NEGSTEER_INTERFACE_METRICS.out.extended_csv

    // ── Build per-survivor manifest ───────────────────────────────────
    EXTRACT_SURVIVOR_MANIFEST(
        extended_csv_ch,
        per_sequence_workdirs_ch,
        Channel.value(file("${projectDir}/bin/extract_survivor_manifest.py"))
    )

    // ── Fan-out: one record per survivor ──────────────────────────────
    // The manifest CSV is consumed via splitCsv.  Each record becomes
    // three downstream tuples — one per orthogonal-metrics stream.
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
            rec[0],              // seq_name
            rec[4],              // receptor_seq
            rec[5],              // effector_seq
            file(rec[3]),        // effector_template_cif path
            file(rec[2]),        // ground_truth_pdb path
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
            rec[0],              // seq_name
            file(rec[1]),        // canonical_pdb
            file(rec[2]),        // ground_truth_pdb (unused but kept symmetric)
        )
    }
    NEGSTEER_BIOPHYSICAL_METRICS(
        biophysical_input_ch,
        Channel.value(file("${projectDir}/bin/run_biophysical_metrics.py"))
    )

    // ── Rosetta stream ────────────────────────────────────────────────
    rosetta_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],              // seq_name
            file(rec[1]),        // canonical_pdb
        )
    }
    // FastRelax XML — pinned single file shared across all rosetta
    // tasks.  Value channel re-emits per consumer so each fan-out
    // task gets its own staged copy in its work dir.
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
        AF3_PARSE_OUTPUT.out.summary_csv.collect(),
        NEGSTEER_BIOPHYSICAL_METRICS.out.summary_csv.collect(),
        NEGSTEER_ROSETTA_METRICS.out.summary_csv.collect(),
        Channel.value(file("${projectDir}/bin/merge_orthogonal_metrics.py"))
    )

    // ── Orthogonal-metrics diagnostic plots ───────────────────────────
    // Mirrors the Step-4b NEGSTEER_PLOTS/NEGSTEER_WITHIN_SEQUENCE_PLOTS
    // pattern: render the diagnostic suite as a parallel terminal stage
    // alongside the metric CSVs, so the wet-lab triage view is
    // immediately available without a separate plotting step.
    ORTHOG_PLOTS(
        NEGSTEER_ORTHOGONAL_METRICS.out.final_csv,
        NEGSTEER_CROSS_SEQUENCE.out.cross_summary,
        Channel.value(file("${projectDir}/bin/orthogonal_metrics_plots.py")),
    )

    // Terminal outputs:
    //   ${params.outdir}/negative_steering/cross_sequence_summary.csv
    //   ${params.outdir}/negative_steering/cross_sequence_summary_with_interface_metrics.csv
    //   ${params.outdir}/orthogonal_metrics/survivors_with_orthogonal_metrics.csv
    //   ${params.outdir}/plots/orthogonal_*.png  (4 diagnostic plots)
    //   plus per-sequence workdirs under ${params.outdir}/negative_steering/runs/
    //   and per-survivor AF3/biophysical/Rosetta artifacts under
    //   ${params.outdir}/orthogonal_metrics/{af3_nomsa,biophysical,rosetta}/<seq_name>/

    } // end if (!params.stop_after_rosetta)
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    // ── Task 6 / P0-6 control miscalibration check ────────────────
    // Scan the final cross_sequence_summary.csv for control rows
    // (row_type starting with 'control_').  Any control where the
    // ranker thinks the sequence is a good binder — low ra_eff or
    // high ipSAE — is a sign the pipeline is sequence-blind on this
    // target and the ranker thresholds need retightening.  Per the
    // P0-6 acceptance: WARN, do not fail the run.
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
    Pipeline Complete
    =============================================================
    Project:    ${params.project_name}
    Output:     ${params.outdir}
    Duration:   ${workflow.duration}
    Success:    ${workflow.success}

    Terminal outputs:
      negative_steering/cross_sequence_summary.csv
        — tier-then-composite ranked summary, one row per MPNN sequence
      negative_steering/cross_sequence_summary_with_interface_metrics.csv
        — same rows, plus P0-29 columns (iRMSD, fnat, DockQ, ipSAE_15, intact_core)
      orthogonal_metrics/survivors_with_orthogonal_metrics.csv
        — survivors only, plus P0-31 columns (AF3-no-MSA, BSA, Sc, ΔΔG, ...)
          and the orthogonal_flags / passes_orthogonal_filters cascade columns
    =============================================================
    """.stripIndent()
    log.info banner

    if (control_warnings) {
        log.warn """
        =============================================================
        ⚠  RANKER MISCALIBRATION WARNING (Task 6 controls)
        =============================================================
        ${control_warnings.size()} control row(s) scored as good binders.
        This suggests the ranker is sequence-blind on this target — a
        scrambled or polyA receptor should NOT pass binding criteria.
        Investigate before promoting any steered design to wet lab:
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
