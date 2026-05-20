#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_pose_solver.nf — Isolated test for the Pose Solver module
 * =============================================================================
 * Exercises the full Branch A docking chain on real biological inputs
 * (Pikp-1_HMA + avr-Pia) using the pair set validated in
 * tests/pose_solver/run_example.sh:
 *
 *   POSE_SOLVER_PREPARE → POSE_SOLVE → POSE_INTERFACE_METRICS
 *                                    → POSE_SOLVER_PLOTS → BUILD_CONTIGS
 *
 * Usage:
 *   sbatch tests/pose_solver/run_test_pose_solver.slurm.sh
 *
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults
// ---------------------------------------------------------------------------

// Real biological inputs: Pikp-1_HMA (78-residue HMA domain) and
// avr-Pia (68-residue MAX-fold effector) as separate monomer PDBs.
// Both renumbered to start at residue 1 in ChimeraX before saving.
params.receptor_input    = "${projectDir}/data/Pikp-1_HMA.pdb"
params.effector_input    = "${projectDir}/data/avr-pia.pdb"
params.receptor_chain    = "A"
params.effector_chain    = "B"

// Contig: A1-32 + A50-68 fixed, two de novo regions (33-49 → 10-30,
// 69-78 → 10 residues).
params.contigs           = "A1-32/10-30/A50-68/10-10 B"

// Pose solver constraints — same pair set as run_example.sh.
//   A73-B31, A71-B33: anti-parallel β-strand H-bond contacts
//   A8-B22:           rotation anchor that builds shape complementarity
params.pose_solver_pairs               = "A73-B31 A71-B33 A8-B22"
params.pose_solver_exclusions          = "A8-B33@4.5"
params.pose_solver_min_pair_distance   = 3.5
params.pose_solver_max_pair_distance   = 6.0
params.pose_solver_pair_sc_clash_cutoff = 2.0
params.pose_solver_clash_cutoff        = 2.0
params.pose_solver_contact_cutoff      = 8.0
// Short test run — 200 restarts converges in <1 min on the validated
// avr-Pia case; production default is 1000.
params.pose_solver_n_restarts          = 200
params.pose_solver_use_de              = false
params.pose_solver_global_interp       = true
params.pose_solver_interp_weight       = 10.0
params.pose_solver_contig_design_region = "33-49,69-78"

params.rfdiff_contact_cutoff = 8.0
params.project_name          = "test_pose_solver"
params.outdir                = "${projectDir}/results"

// Infrastructure — container paths inherited from nextflow.config.

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { POSE_SOLVER_PREPARE                       } from '../../modules/pose_solver'
include { POSE_SOLVE                                } from '../../modules/pose_solver'
include { POSE_INTERFACE_METRICS                    } from '../../modules/pose_solver'
include { POSE_SOLVER_PLOTS                         } from '../../modules/pose_solver'
include { BUILD_CONTIGS                             } from '../../modules/pose_solver'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_REC } from '../../modules/preprocessing'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_EFF } from '../../modules/preprocessing'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    rec_ch = Channel.fromPath(params.receptor_input, checkIfExists: true)
    eff_ch = Channel.fromPath(params.effector_input, checkIfExists: true)

    WRITE_DUMMY_MAPPING_REC(Channel.value("receptor"))
    WRITE_DUMMY_MAPPING_EFF(Channel.value("effector"))

    POSE_SOLVER_PREPARE(
        rec_ch,
        eff_ch,
        params.receptor_chain,
        params.effector_chain,
        Channel.value(file("${projectDir}/../../bin/pose_solver_prepare.py"))
    )

    POSE_SOLVE(
        POSE_SOLVER_PREPARE.out.receptor_pdb_out,
        POSE_SOLVER_PREPARE.out.effector_pdb_out,
        params.receptor_chain,
        params.effector_chain,
        params.pose_solver_pairs,
        params.pose_solver_exclusions,
        params.pose_solver_min_pair_distance,
        params.pose_solver_max_pair_distance,
        params.pose_solver_pair_sc_clash_cutoff,
        params.pose_solver_clash_cutoff,
        params.pose_solver_contact_cutoff,
        params.pose_solver_n_restarts,
        params.pose_solver_use_de,
        params.pose_solver_global_interp,
        params.pose_solver_interp_weight,
        params.pose_solver_contig_design_region,
        Channel.value(file("${projectDir}/../../bin/pose_solver.py"))
    )

    POSE_INTERFACE_METRICS(
        POSE_SOLVE.out.posed_pdb,
        params.receptor_chain,
        params.effector_chain,
        Channel.value(file("${projectDir}/../../bin/pose_interface_metrics.py"))
    )

    POSE_SOLVER_PLOTS(
        POSE_SOLVE.out.results_json,
        POSE_INTERFACE_METRICS.out.metrics_json,
        POSE_SOLVE.out.restart_losses,
        POSE_SOLVE.out.posed_pdb,
        params.receptor_chain,
        params.effector_chain,
        Channel.value(file("${projectDir}/../../bin/pose_solver_plots.py"))
    )

    BUILD_CONTIGS(
        POSE_SOLVE.out.posed_pdb,
        params.receptor_chain,
        params.effector_chain,
        params.contigs,
        WRITE_DUMMY_MAPPING_REC.out.mapping,
        WRITE_DUMMY_MAPPING_EFF.out.mapping,
        Channel.value(file("${projectDir}/../../bin/build_contigs.py"))
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    =============================================================
    Pose Solver Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pose solver test failed: ${workflow.errorMessage}"
}
