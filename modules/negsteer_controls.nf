/*
 * =============================================================================
 * modules/negsteer_controls.nf
 *
 * Task 6 (P0-6, revised in v7) · Scrambled and polyA negative controls.
 *
 * TWO negative-control variants are generated once per cohort and run
 * through the SAME negative-steering inner loop as a real MPNN sequence,
 * BOTH against the input structure (the receptor+effector complex going
 * into RFDiffusion — "rfdiffusion_input.pdb"):
 *
 *   control_scrambled  — the design region of the input receptor
 *                        sequence is permuted under a deterministic
 *                        seed, preserving length and amino-acid
 *                        composition.  Tests the ranker's bias toward
 *                        sequence-pattern matching independent of
 *                        geometric design intent.
 *
 *   control_polyA      — design-region positions are substituted with
 *                        alanine, EXCEPT where the source receptor
 *                        already carries glycine (those stay G to avoid
 *                        forcing alanine into glycine-only Φ/ψ wells).
 *                        Tests the ranker's reliance on the designable
 *                        surface for its binding-site call.
 *
 * Outside the design region the receptor is held verbatim in BOTH
 * controls.  The effector chain is passed through unchanged (controls
 * test the receptor binding face, not the effector).
 *
 * The "design region on the input structure" is the gap between anchor
 * segments in the contigs string — i.e. positions on the native
 * receptor that RFDiffusion is being asked to REPLACE with de novo
 * chemistry.  Derived by DERIVE_INPUT_INDICES (process below) via the
 * bin/derive_input_design_region.py script.
 *
 * TOTAL controls per cohort = 2 (not 2N).  The headline question they
 * answer — "if RFDiffusion did a terrible job, what metrics would those
 * sequences get?" — is a property of the target + its contigs, not of
 * any particular RFDiffusion design.  Per-design pathology is covered
 * by Task 8 AUROC instead.
 *
 * Re-uses bin/negative_steering_run_one.sh as the orchestrator so
 * controls go through identical Boltz settings, identical post-
 * processing, and produce an identical passing_summary.csv shape.
 * The only difference is the row_type.txt sidecar this process writes
 * alongside the passing_summary, which cross_sequence_summary.py
 * picks up to mark these rows as 'control_scrambled' or 'control_polyA'.
 * =============================================================================
 */


/*
 * DERIVE_INPUT_INDICES
 * ────────────────────
 * Produce the two index files the NEGSTEER_CONTROLS process consumes:
 *   input_design_region.txt   — 1-based positional indices on the receptor
 *                               chain (matches derive_design_region.py format)
 *   input_true_interface.txt  — 0-based positional indices (matches
 *                               derive_true_interface.py format)
 *
 * Runs once per pipeline invocation (not fan-out).  Consumes the
 * rfdiffusion_input.pdb and the resolved contigs string — both produced
 * by either Branch A (BUILD_CONTIGS) or Branch B (RESOLVE_CONTIGS) of
 * the preprocessing workflow.
 */
process DERIVE_INPUT_INDICES {
    tag "derive_input_indices"
    label 'cpu'

    publishDir "${params.outdir}/negative_steering/controls_inputs", mode: 'copy'

    input:
    path  input_pdb
    val   contigs
    val   receptor_chain
    val   effector_chain
    val   contact_cutoff

    output:
    path "input_design_region.txt",  emit: design_region
    path "input_true_interface.txt", emit: true_interface
    // Re-emit the input PDB so downstream consumers can pick it up
    // from this process's output channels without sharing the
    // upstream rfdiff_pdb_ch with a second operator consumer.
    path "${input_pdb.name}",        emit: input_pdb

    script:
    """
    set -euo pipefail

    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python ${projectDir}/bin/derive_input_design_region.py \\
            --input-pdb           ${input_pdb} \\
            --contigs             "${contigs}" \\
            --receptor-chain      ${receptor_chain} \\
            --effector-chain      ${effector_chain} \\
            --contact-cutoff      ${contact_cutoff} \\
            --design-region-output input_design_region.txt \\
            --true-interface-output input_true_interface.txt

    echo "Input-structure design region + true interface derived:"
    ls -la input_design_region.txt input_true_interface.txt
    """
}


process NEGSTEER_CONTROLS {
    tag "${control_name}"
    label 'gpu'

    publishDir "${params.outdir}/negative_steering/runs", mode: 'copy'

    input:
    // control_name = "input_control_<scrambled|polyA>"
    // control_type = "scrambled" | "polyA"
    // input_pdb    = rfdiffusion_input.pdb (receptor + effector complex)
    // design_region = input_design_region.txt from DERIVE_INPUT_INDICES
    // true_interface = input_true_interface.txt from DERIVE_INPUT_INDICES
    tuple val(control_name), val(control_type),
          path(input_pdb), path(design_region), path(true_interface)
    val  receptor_chain
    val  effector_chain
    val  plan_extra_args
    val  postprocess_rmsd_threshold
    val  postprocess_metric_column
    val  postprocess_contact_cutoff

    output:
    path "${control_name}/", emit: per_control_workdir

    script:
    """
    set -euo pipefail

    mkdir -p ${control_name}/inputs

    # ─── Extract effector sequence from the input PDB ─────────────────
    # The input PDB carries both chains; we need the effector as a
    # single-record FASTA to hand to build_control_sequences.py.  Use
    # the same get_chain_sequence helper the rest of the pipeline
    # uses, so chain-letter conventions stay consistent.
    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python - <<PYEOF
import sys
from pathlib import Path
sys.path.insert(0, "${projectDir}/bin")
from boltz2_negative_steering import get_chain_sequence
seq = get_chain_sequence(Path("${input_pdb}"), "${effector_chain}")
if not seq:
    raise SystemExit("ERROR: effector chain empty in ${input_pdb}")
out = Path("${control_name}/inputs/source_effector.fasta")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(f">effector\\n{seq}\\n")
print(f"wrote {out} ({len(seq)} effector residues)")
PYEOF

    # ─── Build the control receptor + companion FASTAs ───────────────
    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python ${projectDir}/bin/build_control_sequences.py \\
            --source-pdb         "${input_pdb}" \\
            --design-region-file "${design_region}" \\
            --effector-fasta     "${control_name}/inputs/source_effector.fasta" \\
            --receptor-chain     "${receptor_chain}" \\
            --effector-chain     "${effector_chain}" \\
            --output-dir         "${control_name}/inputs"

    # The synthesised files the orchestrator will consume.  build_control_
    # sequences.py writes three single-record FASTAs:
    #   control_scrambled_receptor.fasta
    #   control_polyA_receptor.fasta
    #   effector.fasta
    # We pick the receptor matching this control_type and use the shared
    # effector.fasta for both flavours.
    CONTROL_RECEPTOR="${control_name}/inputs/control_${control_type}_receptor.fasta"
    EFFECTOR_FASTA="${control_name}/inputs/effector.fasta"
    if [[ ! -f "\$CONTROL_RECEPTOR" || ! -f "\$EFFECTOR_FASTA" ]]; then
        echo "ERROR: build_control_sequences did not produce expected files" >&2
        ls -la "${control_name}/inputs/" >&2
        exit 1
    fi

    # ─── Stamp the row_type sidecar BEFORE running the orchestrator ──
    # cross_sequence_summary.py reads this file from the same directory
    # as passing_summary.csv.  Writing it up-front (rather than after
    # the orchestrator finishes) means even a crashed inner loop leaves
    # the workdir in a state where the partial output is correctly
    # identified as a control rather than mistaken for a steered run.
    echo "control_${control_type}" > ${control_name}/row_type.txt

    # ─── Invoke the inline orchestrator ──────────────────────────────
    # Same flags as NEGSTEER_RUN_ONE.  --ground-truth is the input PDB
    # because, for the controls, the "true" complex is the input
    # structure itself — we're asking "what does Boltz say about this
    # input structure if we wreck the native sequence at the would-be
    # design region?" and the true-interface indices were derived
    # against this exact PDB by DERIVE_INPUT_INDICES.
    bash ${projectDir}/bin/negative_steering_run_one.sh \\
        --seq-name                     "${control_name}" \\
        --ground-truth                 "${input_pdb}" \\
        --receptor-chain               "${receptor_chain}" \\
        --effector-chain               "${effector_chain}" \\
        --receptor-fasta               "\$CONTROL_RECEPTOR" \\
        --effector-fasta               "\$EFFECTOR_FASTA" \\
        --true-interface-indices-file  "${true_interface}" \\
        --design-region-indices-file   "${design_region}" \\
        --workdir                      "${control_name}" \\
        --bin-dir                      "${projectDir}/bin" \\
        --boltz-container              "${params.boltz2_container}" \\
        --plan-extra-args              "${plan_extra_args}" \\
        --postprocess-rmsd-threshold   "${postprocess_rmsd_threshold}" \\
        --postprocess-metric-column    "${postprocess_metric_column}" \\
        --postprocess-contact-cutoff   "${postprocess_contact_cutoff}"

    # Same invariant check as NEGSTEER_RUN_ONE — fail loudly if the
    # final deliverable is missing.
    if [[ ! -f "${control_name}/passing_summary.csv" ]]; then
        echo "ERROR: ${control_name}/passing_summary.csv not produced" >&2
        ls -la ${control_name}/ || true
        exit 1
    fi
    """
}
