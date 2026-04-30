/*
 * =============================================================================
 * modules/negsteer_rosetta_metrics.nf
 *
 * P0-31 · Rosetta orthogonal metrics on survivors.
 *
 * Two-stage Rosetta call per survivor:
 *
 *   1. FastRelax (1 repeat, Cartesian, ref2015_cart) on the bound
 *      complex to remove clashes and bring the predicted structure
 *      onto Rosetta's energy surface.  XML pinned at
 *      bin/fastrelax_for_ia.xml; Bennett 2023's calibration setup.
 *   2. InterfaceAnalyzer on the relaxed PDB to compute Lawrence-
 *      Colman shape complementarity (Sc) and the separated ΔG
 *      (rosetta_ddg).
 *
 * Without the relax pre-treatment, ΔΔG values are inflated 10-20 REU
 * by clashes and idealised-bond-geometry artefacts in the predicted
 * input — so the resulting numbers are not comparable to published
 * thresholds (Bennett 2023's −30 REU for de novo binders) and not
 * particularly meaningful in absolute terms.  With FastRelax, they
 * are.  Both stages live inside run_rosetta_metrics.py; this module
 * just wires the inputs.
 *
 * Metrics emitted:
 *   sc                Lawrence-Colman shape complementarity (> 0.55 typical).
 *   rosetta_ddg       InterfaceAnalyzer separated-ΔG on the relaxed pose.
 *
 * Runs inside Rosetta.img.  ~3-6 min per survivor (relax dominates;
 * IA call is ~15 s).  Survivor fan-out is CPU — the queue is
 * jic-medium, not jic-gpu.
 * =============================================================================
 */

process NEGSTEER_ROSETTA_METRICS {
    tag "${seq_name}"
    label 'cpu'

    publishDir "${params.outdir}/orthogonal_metrics/rosetta/${seq_name}",
        mode: 'copy', pattern: "rosetta_summary.csv"

    input:
    tuple val(seq_name),
          path(canonical_pdb)
    // FastRelax XML staged as a process input so Nextflow puts it
    // inside the per-task work dir (already covered by --bind ${PWD}).
    // Without staging, the script's default path resolution would
    // point at ${projectDir}/bin/fastrelax_for_ia.xml — outside the
    // container's bound paths.
    path fastrelax_xml

    output:
    path "rosetta_summary.csv", emit: summary_csv

    script:
    """
    set -euo pipefail

    singularity exec --bind \${PWD}:\${PWD} ${params.rosetta_container} \\
        python ${projectDir}/bin/run_rosetta_metrics.py \\
            --seq-name ${seq_name} \\
            --canonical-pdb ${canonical_pdb} \\
            --receptor-chain ${params.rfdiff_output_receptor_chain} \\
            --effector-chain ${params.rfdiff_output_effector_chain} \\
            --fast-relax-xml ${fastrelax_xml} \\
            --output-csv rosetta_summary.csv
    """
}