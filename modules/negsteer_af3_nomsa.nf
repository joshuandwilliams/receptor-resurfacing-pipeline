/*
 * =============================================================================
 * modules/negsteer_af3_nomsa.nf
 *
 * P0-31 · AF3-no-MSA orthogonal validation on survivors.
 *
 * Three processes:
 *
 *   AF3_SETUP_DB               — idempotently build $HOME/af3_db symlink farm
 *                                (combined v3.0.0 + v2.3.2 BFD/MGnify).
 *                                Runs once, feeds all survivor fan-outs.
 *
 *   AF3_NOMSA_ON_SURVIVORS     — per survivor, write an AF3 JSON with:
 *                                  modelSeeds    = params.af3_nomsa_seeds (3 by default)
 *                                  unpairedMsa   = ""     (no MSA)
 *                                  pairedMsa     = ""
 *                                  templates     = effector template CIF
 *                                                   (receptor gets no template,
 *                                                    matching Boltz)
 *                                then run `run_alphafold.py --norun_data_pipeline`
 *                                and emit the entire output/ directory plus the
 *                                JSON that was used.
 *
 *   AF3_PARSE_OUTPUT           — per survivor, parse the mmCIF predictions and
 *                                compute ra_eff vs the RFDiffusion reference.
 *                                Emits one CSV row per survivor keyed by
 *                                mpnn_sequence.
 *
 * Parallelisation: AF3 runs natively (source package) — no container.
 * maxForks is capped at params.max_af3_parallel (default 30), much wider
 * than Boltz because AF3 prediction without MSA is GPU-only and brief
 * (no long data-pipeline phase).
 *
 * Parameter alignment with Boltz where possible:
 *   - seeds: 3 (matches Boltz num_seeds=3)
 *   - diffusion_samples: 5 per seed (AF3 default; Boltz matches at 5)
 *   - recycles: 10 (AF3 default; Boltz uses 3 — irreducible mismatch on
 *     the NBI AF3 build which doesn't expose a CLI flag for this)
 *   - effector template: YES (forced, same CIF Boltz used)
 *   - receptor template: NO (matches Boltz)
 * =============================================================================
 */


/*
 * AF3_SETUP_DB — build the combined AF3 database directory.
 *
 * Same pattern as the standalone af3.nf module: on NBI the databases
 * are split across two reference-data dirs (v3.0.0 and v2.3.2), and
 * AF3 requires a single --db_dir.  This process builds that flat union
 * via symlinks under $HOME/af3_db.
 *
 * Emits a sentinel path so downstream processes can depend on it
 * without passing the full directory through the channel.
 */
process AF3_SETUP_DB {
    tag "${params.project_name}"
    label 'cpu'

    output:
    path "af3_db_ready.flag", emit: ready_flag
    val  "${System.getenv('HOME')}/af3_db", emit: db_dir

    script:
    """
    set -euo pipefail

    AF3_DATA_DIR="\${HOME}/af3_db"
    mkdir -p "\${AF3_DATA_DIR}"

    # ─── Link v3.0.0 databases ────────────────────────────────────────
    for f in ${params.af3_db_v3}/*; do
        ln -sfn "\$f" "\${AF3_DATA_DIR}/\$(basename "\$f")"
    done

    # ─── Link BFD/MGnify from v2.3.2 (AF3 still needs these) ─────────
    ln -sfn ${params.af2_data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta \\
        "\${AF3_DATA_DIR}/bfd-first_non_consensus_sequences.fasta"
    ln -sfn ${params.af2_data_dir}/mgnify/mgy_clusters_2022_05.fa \\
        "\${AF3_DATA_DIR}/mgy_clusters_2022_05.fa"

    echo "AF3 database dir ready: \${AF3_DATA_DIR}"
    # Trailing `|| true` swallows SIGPIPE (141) propagated by pipefail
    # when head closes early.  Same pattern as line 224 below.
    ls "\${AF3_DATA_DIR}" | head -20 || true

    touch af3_db_ready.flag
    """
}


/*
 * AF3_NOMSA_ON_SURVIVORS — one GPU job per survivor.
 *
 * Input tuple per survivor:
 *   seq_name              (e.g. "design_3_seq_1")
 *   receptor_seq          MPNN-designed receptor (post-reversion)
 *   effector_seq          native effector from ground truth
 *   effector_template_cif the same CIF Boltz used (reused from workdir)
 *   ground_truth_pdb      RFDiffusion reference, for downstream parse step
 *
 * The AF3 seed list is passed as a Groovy list param (params.af3_nomsa_seeds)
 * and rendered into the modelSeeds JSON field.  AF3 does not take the seed
 * count via CLI, so this is the only knob.
 */
process AF3_NOMSA_ON_SURVIVORS {
    tag "${seq_name}"
    label 'gpu'

    publishDir "${params.outdir}/orthogonal_metrics/af3_nomsa/${seq_name}",
        mode: 'copy', pattern: '{input.json,output/**}'

    input:
    tuple val(seq_name),
          val(receptor_seq),
          val(effector_seq),
          path(effector_template_cif),
          path(ground_truth_pdb)
    val  af3_db_dir
    path db_ready_flag

    output:
    tuple val(seq_name),
          path("output"),
          path("input.json"),
          path(ground_truth_pdb),
          emit: prediction

    script:
    // Render the seed list as a JSON array.  Nextflow can stringify a
    // Groovy list, but we want the raw JSON form inside the heredoc —
    // join explicitly to avoid relying on toString() formatting quirks.
    def seeds_json = "[" + params.af3_nomsa_seeds.collect { it.toString() }.join(", ") + "]"
    """
    set -euo pipefail

    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export TF_FORCE_UNIFIED_MEMORY=true
    export XLA_CLIENT_MEM_FRACTION=3.2

    mkdir -p output

    # ─── Load AF3 environment (native, not via container) ────────────
    source package ${params.af3_package_id}

    # ─── Write AF3 JSON input ────────────────────────────────────────
    # Effector template passed via the `templates` field on chain B,
    # matching Boltz's force-template convention.  Receptor (chain A)
    # has no template, matching Boltz.
    #
    # unpairedMsa and pairedMsa empty + --norun_data_pipeline below
    # together implement the no-MSA mode.
    #
    # AF3's input parser expects each template entry to provide the
    # mmCIF inline under the "mmcif" key (the older "mmcifPath" key
    # was removed).  Inlining requires JSON-escaping the CIF text
    # (newlines, quotes, etc.), so build the JSON in Python with
    # stdlib json instead of a bash heredoc.  Nextflow has already
    # interpolated \${seq_name}, \${seeds_json}, \${receptor_seq},
    # \${effector_seq}, and \${effector_template_cif} into the script
    # body before bash sees it; the PYEOF heredoc is single-quoted so
    # bash performs no further expansion inside.
    python3 - <<'PYEOF' > input.json
import json
with open("${effector_template_cif}") as f:
    effector_cif_text = f.read()
payload = {
    "name": "${seq_name}",
    "modelSeeds": ${seeds_json},
    "dialect": "alphafold3",
    "version": 1,
    "sequences": [
        {
            "protein": {
                "id": "A",
                "sequence": "${receptor_seq}",
                "unpairedMsa": "",
                "pairedMsa": "",
                "templates": [],
            }
        },
        {
            "protein": {
                "id": "B",
                "sequence": "${effector_seq}",
                "unpairedMsa": "",
                "pairedMsa": "",
                "templates": [
                    {
                        "mmcif": effector_cif_text,
                        "queryIndices": [],
                        "templateIndices": [],
                    }
                ],
            }
        },
    ],
}
print(json.dumps(payload, indent=2))
PYEOF

    echo "AF3-no-MSA prediction for ${seq_name}"
    echo "Seeds: ${seeds_json}"
    echo "Effector template: ${effector_template_cif}"
    # Show JSON structure with the (potentially huge) inlined mmcif
    # text replaced by a length stub for readability.
    python3 - <<'PYEOF'
import json
with open("input.json") as f:
    d = json.load(f)
for chain in d["sequences"]:
    for tpl in chain["protein"].get("templates", []):
        if "mmcif" in tpl:
            tpl["mmcif"] = f"<{len(tpl['mmcif'])} chars of mmcif text>"
print(json.dumps(d, indent=2))
PYEOF

    run_alphafold.py \\
        --json_path=\${PWD}/input.json \\
        --model_dir="${params.af3_model_dir}" \\
        --db_dir="${af3_db_dir}" \\
        --output_dir=\${PWD}/output \\
        --norun_data_pipeline

    echo "AF3 output files for ${seq_name}:"
    find output -type f | head -30 || true
    """
}


/*
 * AF3_PARSE_OUTPUT — post-process a single AF3 prediction.
 *
 * For each seed×diffusion-sample pair AF3 wrote under output/, load the
 * mmCIF and compute receptor-aligned effector RMSD vs the RFDiffusion
 * reference.  Emit a single-row CSV per survivor with:
 *   seq_name, af3_nomsa_best_ra_eff, af3_nomsa_mean_ra_eff,
 *   af3_nomsa_n_correct_interface, af3_nomsa_best_iptm,
 *   af3_nomsa_mean_iptm, af3_nomsa_failures
 *
 * Runs in boltz2_negsteer.img for its numpy/gemmi/biopython stack.
 * CPU-only, fast (< 10 s per survivor).
 */
process AF3_PARSE_OUTPUT {
    tag "${seq_name}"
    label 'cpu'

    publishDir "${params.outdir}/orthogonal_metrics/af3_nomsa/${seq_name}",
        mode: 'copy', pattern: "af3_nomsa_summary.csv"

    input:
    tuple val(seq_name),
          path(output_dir),
          path(input_json),
          path(ground_truth_pdb)
    // Stage the parser script as a path input so Nextflow's task-cache
    // hash includes its CONTENT.  Without this, edits to
    // bin/parse_af3_output.py do not invalidate the cache (the script
    // command line interpolated below is byte-identical regardless of
    // the file's contents at ${projectDir}/bin/), and `-resume` will
    // silently run on stale parses.  This was the cause of the
    // 2026-04-29 silent-cache-hit incident where a parser bugfix
    // appeared to deploy successfully but didn't actually re-execute.
    path parse_script

    output:
    path "af3_nomsa_summary.csv", emit: summary_csv

    script:
    """
    set -euo pipefail

    singularity exec --bind \${PWD}:\${PWD} ${params.boltz2_container} \\
        python ${parse_script} \\
            --seq-name ${seq_name} \\
            --af3-output-dir ${output_dir} \\
            --ground-truth ${ground_truth_pdb} \\
            --receptor-chain ${params.rfdiff_output_receptor_chain} \\
            --effector-chain ${params.rfdiff_output_effector_chain} \\
            --ra-eff-threshold ${params.orthogonal_filter_af3_ra_max} \\
            --output-csv af3_nomsa_summary.csv
    """
}