# 17 — Storage bloat audit (negative steering focus)

Audit of the storage footprint of the `receptor-resurfacing-pipeline` runs
under `experiments/campaigns/*/runs/`, motivated by two observations: (a)
finished runs occupy several GB despite being conceptually small, and (b)
deleting a single finished run takes ~15 minutes due to high inode count.
The audit is scoped to the negative-steering subsystem (which the human
identified as the largest consumer) and uses one representative run,
`pikp1_avrpia/runs/v1_4b`, with a confirmation pass on `v1_4a` of the same
campaign.

This report does **not** implement code-level remediation; it describes the
problem, lists the cleanup that should be applied to already-finished runs
(see §6), and sketches a forward-looking code change that should be made
before further runs are started (§7).

## Summary

Five remediable causes of bloat are identified. Together they account for
the majority of bytes and ~85% of the inodes in a finished run.

| # | Cause                                                          | Bytes recoverable per run | Inodes recoverable | Ease |
|---|----------------------------------------------------------------|---------------------------|--------------------|------|
| 1 | `NEGSTEER_RUN_ONE` `publishDir mode: 'copy'` duplicates work→results | ~1.5 GB                  | ~150,000           | 5/5 (delete `work/`) |
| 2 | Boltz keeps 5 diffusion samples; only the top-ranked (`model_0`) is read | ~1.0 GB                  | ~135,000           | 4/5 (post-hoc prune) |
| 3 | Boltz scratch dirs (`processed/`, `lightning_logs/`, `msa/`) are kept    | ~50 MB                   | ~40,000            | 4/5 (post-hoc prune) |
| 4 | AF3 15-sample ensemble retained in `work/` only                          | ~1.0 GB                  | ~7,500             | 5/5 (delete `work/`) |
| 5 | `work/`, `.nextflow/cache/`, `.nextflow/plr/` retained after completion  | ~3 GB depending on stage | ~250,000           | 5/5 (delete) |

Realistic combined recovery on `v1_4b`: roughly **3 GB → 0.5 GB** (≈ 80%
fewer bytes), and **300k → 50k inodes** (≈ 85% fewer files). The 15-minute
delete on a finished run should drop to ~2–3 minutes.

The largest single issue is #5 (delete `work/` after completion). Issues #1,
#2, #3 cascade: each Boltz call's bloat is multiplied by ~18 calls per
sequence, ~100 sequences per run, and ×2 by `mode: 'copy'`.

## 1 — Methodology

The HPC-Home filesystem is mounted via SSHFS, so `find`, `du`, and recursive
listings are slow (sequential per-file metadata roundtrips). Sizing was
therefore inferred without traversal:

- **`results/trace.txt`** — Nextflow's per-task trace records `wchar`
  (write-bytes-by-process) and `rchar` per task. This gives the cumulative
  output volume per process across a run without traversing the tree.
- **`ls -la <single dir>`** — file sizes per direct child, summed by
  inspection. Used to spot-check individual Boltz output dirs and confirm
  the trace.txt picture.
- **Source code reading** — `bin/boltz2_negative_steering.py:799-863`
  (the `_invoke_boltz` function), `bin/boltz2_iterate_steering.py:1412-1437`
  (metrics computation), `modules/negative_steering.nf:122-229`
  (`NEGSTEER_RUN_ONE` process declaration). Determined which files in the
  Boltz output are actually consumed by downstream code.

For one task only (a representative `NEGSTEER_RUN_ONE` work directory),
`ls -R | wc -l` was run to confirm file count: **1,605 entries** for one
sequence's workdir.

## 2 — The negsteer output structure

Per pipeline run, `NEGSTEER_RUN_ONE` is fired once per surviving MPNN
sequence (62–118 sequences depending on the campaign). Each task creates
the following layout under `results/negative_steering/runs/design_X_seq_Y/`
(also duplicated in `work/<hash>/design_X_seq_Y/`):

```
design_X_seq_Y/
├── inputs/                       # tiny — receptor.fasta, effector.fasta
├── logs/
├── pathways.json, passing_summary.csv, …  # tiny aggregator outputs
└── cycle_0/
    ├── plan.json, prefilter.json, contaminated.json, reversion_plan.json,
    │   reversion_results*.json, summary.txt, *.csv         # tiny
    ├── effector_template.cif                               # ~20 KB
    ├── initial_prediction.pdb, initial_prediction_s1.pdb,
    │   initial_prediction_s2.pdb                           # ~88 KB × 3
    ├── true_interface_residues.txt, wrong_interface_residues.txt
    ├── initial/                  # ← Boltz raw output, seed 0
    ├── initial_s1/               # ← Boltz raw output, seed 1
    ├── initial_s2/               # ← Boltz raw output, seed 2
    ├── steered/                  # one subfolder per (candidate × seed)
    │   ├── design_00_s0/
    │   │   ├── prediction.pdb, mutations.tsv, result.json   # tiny
    │   │   └── boltz_results_input/   # ← Boltz raw output
    │   ├── design_00_s1/  …
    │   └── design_02_s2/         # 3 candidates × 3 seeds = 9 subfolders
    ├── reversions/               # 0–N reversion candidates × 3 seeds
    │   └── rev_design_NN_sX_sY/
    │       └── boltz_results_input/   # ← Boltz raw output
    └── contamination_scratch/    # tiny CSVs
```

Each `boltz_results_input/` directory is the raw `--out_dir` Boltz wrote
to. Its internal layout is fixed by Boltz, not by this pipeline:

```
boltz_results_input/
├── lightning_logs/version_0/hparams.yaml      # 4 KB — not read after Boltz
├── msa/                                       # empty dir
├── processed/
│   ├── manifest.json                          # 0.6 KB — not read
│   ├── structures/input.npz                   # 4.6 KB — not read
│   ├── msa/input_*.npz                        # ~1.5 KB — not read
│   ├── mols/input.pkl                         # 5 B — not read
│   └── templates/                             # empty
└── predictions/input/                         # the only consumed dir
    ├── input_model_0.pdb                      # 88 KB — KEPT
    ├── input_model_1.pdb                      # 88 KB — UNUSED
    ├── input_model_2.pdb                      # 88 KB — UNUSED
    ├── input_model_3.pdb                      # 88 KB — UNUSED
    ├── input_model_4.pdb                      # 88 KB — UNUSED
    ├── confidence_input_model_0.json          # 0.65 KB — KEPT
    ├── confidence_input_model_[1-4].json      # 4 × 0.65 KB — UNUSED
    ├── pae_input_model_0.npz                  # 73 KB — KEPT
    ├── pae_input_model_[1-4].npz              # 4 × 73 KB — UNUSED
    ├── pde_input_model_0.npz                  # 65 KB — KEPT
    ├── pde_input_model_[1-4].npz              # 4 × 65 KB — UNUSED
    ├── plddt_input_model_0.npz                # 0.74 KB — KEPT
    └── plddt_input_model_[1-4].npz            # 4 × 0.74 KB — UNUSED
```

Total per Boltz call: ~1.15 MB; of that ~910 KB (≈ 80%) is models 1–4.

A fully-running negsteer sequence (i.e. one that reaches reversion) makes
**3 + 9 + 6 = 18 Boltz calls**, so ~21 MB of Boltz raw output per sequence
of which ~16 MB is models 1–4. Multiply by 50 fully-running sequences →
~0.8 GB of unused Boltz samples per run. Doubled by issue #1 to ~1.6 GB.

## 3 — Issue #1: `mode: 'copy'` duplicates the workdir

**Where.** `modules/negative_steering.nf:126`:

```groovy
publishDir "${params.outdir}/negative_steering/runs", mode: 'copy'
…
output:
path "${seq_name}/", emit: per_sequence_workdir
```

`mode: 'copy'` makes Nextflow byte-copy every file inside `${seq_name}/`
from the task's work directory into `${params.outdir}/negative_steering/runs/`.
Because `output: path "${seq_name}/"` declares the entire folder (rather
than a curated `path "${seq_name}/passing_summary.csv"` etc.), every
intermediate file — including all of §2's Boltz raw output — is copied.

**Confirmed by direct comparison.** For `pikp1_avrpia/runs/v1_4b`:

- `work/14/ad505…/design_44_seq_1/cycle_0/initial_prediction.pdb` is 88,613 bytes.
- `results/negative_steering/runs/design_44_seq_1/cycle_0/initial_prediction.pdb` is 88,613 bytes.
- Same for every other file in the tree. They are independent copies, not symlinks.

Other modules with `publishDir mode: 'copy'` are smaller producers and
present the same duplication-with-`work/` pattern but at lower volume:
`rfdiffusion.nf` (~10 MB doubled), `proteinmpnn.nf`, `negsteer_*.nf`.

**Recovery.** `rm -rf <run>/work` (after the run is `OK`). This is the
single largest win and is the foundation of the cleanup script in §6.

**Why not switch to `mode: 'symlink'`?** It would save the bytes
immediately but ties `results/` lifetime to `work/` lifetime — deleting
`work/` would silently break `results/`. Given that the goal is to delete
`work/`, `mode: 'copy'` is actually the correct choice; the bug is in
*not* deleting `work/` afterwards.

## 4 — Issue #2: only `model_0` of N diffusion samples is consumed

**Where.** `bin/boltz2_negative_steering.py:799-863`:

```python
boltz_args = [
    "boltz", "predict",
    str(yaml_path),
    "--out_dir", str(out_dir),
    "--recycling_steps", str(recycling_steps),
    "--diffusion_samples", str(diffusion_samples),   # 5 in current params
    …
    "--write_full_pae",
    "--use_potentials",
    "--override",
]
…
all_pdbs = sorted(out_dir.rglob("*.pdb"))
…
pdbs.sort(key=lambda p: (
    "model_0" not in p.name and "rank_0" not in p.name,
    p.name,
))
return pdbs[0]
```

`negsteer_diffusion_samples: 5` causes Boltz to run 5 diffusion trajectories,
rank them internally, and emit them as `input_model_0.pdb` … `input_model_4.pdb`
where **`model_0` is the highest-confidence sample** (Boltz's internal
ranking, not just the sample index). The picker preferentially keeps
`model_0`/`rank_0`, which is correct: this is how the multi-sample
ensembling improves prediction quality. Reducing `negsteer_diffusion_samples`
to 1 would skip the ranking and likely produce a lower-quality top pick;
it is **not** a recommended fix.

The samples that are *not* picked (`model_1` through `model_4`) and their
PAE / PDE / plddt / confidence sidecars are written but never read. The
downstream metrics computation at `bin/boltz2_iterate_steering.py:1416`:

```python
sidecar = _locate_boltz_sidecar_pdb(pdb)
pred_dir = sidecar.parent if sidecar is not None else pdb.parent
…
cmd = [
    sys.executable, str(compute_metrics_script),
    "--model", "boltz2",
    "--prediction-dir", str(pred_dir),
    …
]
```

points `compute_metrics.py` at the directory containing the *single* picked
PDB and reads only the corresponding sidecar files for that one model.

**Recovery.** Once Boltz has returned and `pdbs[0]` has been chosen, the
files matching `*_model_[1-9]*` and `*_rank_[1-9]*` inside every
`boltz_results_input/predictions/*/` directory are dead weight and can be
deleted. Two ways to do this:

- **Retroactively** (existing finished runs): cleanup script in §6.
- **Forward-looking** (future runs): inline cleanup in `_invoke_boltz`
  immediately after `pdbs[0]` is selected. See §7.

## 5 — Issue #3: Boltz scratch dirs are never read after Boltz returns

`boltz_results_input/processed/`, `boltz_results_input/lightning_logs/`,
and `boltz_results_input/msa/` are Boltz's internal staging output — none
of them are read by `boltz2_negative_steering.py`,
`boltz2_iterate_steering.py`, or any of the negsteer post-processing
modules. Confirmed by `grep -rn "processed\|lightning_logs" bin/ modules/`.

Bytes-wise this is a small contributor (~12 KB per Boltz call ≈ 50 MB
across a run). **Inode-wise it is significant**: each Boltz call adds
~9 small files plus 6 directory entries, all <1 KB on a network FS.

**Recovery.** Delete these subdirs alongside the model_[1-9] cleanup in §4.

## 6 — Issue #4: AF3 ensemble samples retained in `work/`

`AF3_NOMSA_ON_SURVIVORS` runs AlphaFold-3 with three seeds (`af3_nomsa_seeds:
[42, 123, 456]`) × 5 samples each. Per task, `work/<hash>/output/design_X_seq_Y/`
contains:

- top-level: `model.cif` (97 KB), `confidences.json` (192 KB), `data.json`
  (21 KB), `summary_confidences.json` (0.3 KB), `ranking_scores.csv`
  (0.4 KB), `TERMS_OF_USE.md` (13 KB).
- 15 subdirs (`seed-{42,123,456}_sample-{0..4}/`) each with their own
  per-sample model + confidence files.

The `publishDir` for AF3 (`modules/negsteer_af3_nomsa.nf`) is correctly
selective — `results/orthogonal_metrics/af3_nomsa/design_X_seq_Y/`
contains only `af3_nomsa_summary.csv` + `input.json`. So AF3 does **not**
have the duplication problem of #1. But the full ensemble persists in
`work/` until `work/` is deleted. ~9.8 MB × 99 tasks ≈ 970 MB per run.

**Recovery.** Solved by deleting `work/` (issue #1's fix). No code change
needed.

## 7 — Issue #5: `work/`, `.nextflow/cache/`, `.nextflow/plr/` are kept indefinitely

Nextflow retains `work/` to support `--resume`. The cache directory
(`.nextflow/cache/<run-uuid>/`) maps task hashes to metadata, and
`.nextflow/plr/` is the per-launch resume index. None of this is needed
once the run is complete and a re-run is not anticipated.

**Important relationship:** `.nextflow/cache/` is only useful when `work/`
exists. The cache stores task-hash → output-metadata mappings; on resume,
Nextflow consults the cache to identify cached tasks, then verifies the
matching files in `work/` are still present. If `work/` is deleted, every
cache entry is invalid and Nextflow re-runs every task on resume regardless
of cache state. So once `work/` is deleted, `.nextflow/cache/` and
`.nextflow/plr/` are orphan state with no use to either resume or audit.

`.nextflow/history` (single line per run, status + timing + run-name + the
exact command) is small and useful as audit trail. Keep it.

`.nextflow.log` (run log, can be ~1 MB) is useful for post-mortem of
failed runs. Keep it for now; can be revisited later.

`tmp/` (typically <100 KB, contains `libjansi-64.so` and similar) — keep
or delete, doesn't matter.

## 8 — Run fleet inventory

The 12 run directories under `experiments/campaigns/*/runs/` divide into
three categories based on `.nextflow/history` and `.out` log inspection.

### Detection rules

| Signal                                                      | Means…                                       |
|-------------------------------------------------------------|----------------------------------------------|
| `.nextflow/history` is missing                              | Run was never executed (only `params.yml`)   |
| `.nextflow/history` last column = `OK`                      | Pipeline completed (every task ended OK)     |
| `.nextflow/history` last column = `ERR`                     | One or more tasks failed; stage may be partial|
| `.out` log contains `Pipeline finished:` line               | Wrapper saw `nextflow run` exit 0 (only on `OK`) |
| `results/negative_steering/runs/design_*/` non-empty        | Negsteer produced output worth cleaning      |
| `results/orthogonal_metrics/af3_nomsa/` exists & non-empty  | AF3 stage completed                          |

Empirically, in the current fleet:

- **All `OK` runs** also reached AF3 (`af3_nomsa/` present and non-empty).
  None of them exited early due to filter-failure-with-zero-survivors.
  The user's hypothetical "OK with nothing past Rosetta" did not occur in
  this batch — but the cleanup script should still tolerate it (use
  `find … -delete`, which is a no-op on missing paths, never `rm` of a
  hard-coded path).
- **The one `ERR` run** (`pikp1_pwt7/v1_4e`) reached negsteer (with 0
  designs surviving the filter) and then died at the
  `NEGSTEER_WITHIN_SEQUENCE_PLOTS` step because there was nothing to plot.
  Its `work/` is large and has nothing to recover by re-running, so it
  should be cleaned.

### Categorised list

**Eligible for cleanup (8):**

| Run                              | Status | Negsteer designs | AF3 ran | Notes                |
|----------------------------------|--------|-----------------:|---------|----------------------|
| `pikp1_avrpia/runs/v1_4a`        | OK     | 62               | yes     |                      |
| `pikp1_avrpia/runs/v1_4b`        | OK     | 102              | yes     | representative run   |
| `pikp1_avrpikf/runs/v1_4a`       | OK     | 118              | yes     |                      |
| `pikp1_pby2/runs/v1_4a`          | OK     | 114              | yes     |                      |
| `pikp1_pby2/runs/v1_4e`          | OK     | 116              | yes     |                      |
| `pikp1_pwt3/runs/v1_4e`          | OK     | 26               | yes     |                      |
| `pikp1_pwt7/runs/v1_4a`          | OK     | 30               | yes     |                      |
| `pikp1_pwt7/runs/v1_4e`          | ERR    | 0                | no      | 0 designs passed; clean for `work/` only |

**Not eligible — never executed (4), preserved for future runs:**

- `pikp1_avrpikf/runs/v2_4b`
- `pikp1_pby2/runs/v1_4b`
- `pikp1_pwt3/runs/v1_4b`
- `pikp1_pwt7/runs/v1_4b`

The cleanup driver hard-codes the 8 eligible paths above. The 4 unexecuted
runs are skipped — they have no `.nextflow/history`, so the per-run script
will refuse to act on them defensively in case the driver list ever drifts.

## 9 — Cleanup actions (retroactive)

For each eligible run, perform in order:

1. **Delete `work/` and Nextflow caches** — recovers the bulk of bytes and
   inodes (issues #1, #4, #5).

   ```
   rm -rf <run>/work
   rm -rf <run>/.nextflow/cache
   rm -rf <run>/.nextflow/plr
   ```

2. **Prune unused diffusion samples in `results/`** — recovers ~80% of the
   bytes inside `boltz_results_input/predictions/*/` (issue #2).

   ```
   find <run>/results/negative_steering/runs \
     -path '*/boltz_results_input/predictions/*' \
     \( -name '*_model_[1-9]*' -o -name '*_rank_[1-9]*' \) \
     -delete
   ```

   The `*_model_[1-9]*` glob matches `input_model_1.pdb`,
   `pae_input_model_1.npz`, `pde_input_model_1.npz`,
   `plddt_input_model_1.npz`, and `confidence_input_model_1.json` — i.e.
   the PDB and all four sidecars per sample. Same for samples 2–4
   (and 5–9, which the current params don't generate but the pattern
   handles defensively if the param ever increases).

3. **Prune Boltz scratch dirs in `results/`** — recovers inodes (issue #3).

   ```
   find <run>/results/negative_steering/runs \
     -type d \( -name processed -o -name lightning_logs -o -name msa \) \
     -path '*/boltz_results_input/*' \
     -prune -exec rm -rf {} +
   ```

4. **Verify the user-facing PDBs survived.** Sanity check that the four
   ChimeraX-relevant PDB types are still present:

   - `results/rfdiffusion/design_*.pdb` (RFDiffusion outputs)
   - `results/negative_steering/runs/*/cycle_0/initial_prediction*.pdb` (cold-start)
   - `results/negative_steering/runs/*/cycle_0/steered/*/prediction.pdb` (steered)
   - `results/negative_steering/runs/*/cycle_0/reversions/*/prediction.pdb` (reverted)

   None of these are touched by steps 1–3.

The cleanup script is at `scripts/clean_finished_run.sh` and the driver is
at `scripts/clean_all_finished_runs.sh`.

### Where to run the script

The HPC-Home filesystem is exposed to the local Mac via SSHFS, which adds
a network roundtrip per file metadata operation. `find … -delete` over
SSHFS is *catastrophically* slow for 100k-file trees — minutes per run
becomes hours. **Run the cleanup script directly on the HPC login node**
(`ssh hpc-home`), where the filesystem is local. The script itself is
filesystem-agnostic; only the runtime cost differs.

## 10 — Forward-looking code change (recommended, not yet implemented)

Once the retroactive cleanup is done, future runs should not re-create the
issue. The single most valuable code change is to inline the §9.2/§9.3
prune at the moment the picker chooses `model_0`. Suggested patch site:
`bin/boltz2_negative_steering.py`, immediately after the `pdbs.sort(...)
return pdbs[0]` block at line 858–863.

Sketch (uncommitted):

```python
def _invoke_boltz(...):
    ...
    pdbs.sort(key=lambda p: (
        "model_0" not in p.name and "rank_0" not in p.name,
        p.name,
    ))
    chosen = pdbs[0]

    # Prune unused diffusion samples and Boltz internal scratch
    # before returning. The picker has already committed to `chosen`;
    # nothing downstream reads the other samples or the scratch dirs.
    chosen_stem = chosen.stem  # e.g. "input_model_0"
    for p in out_dir.rglob("*"):
        if p.is_file():
            name = p.name
            # keep only files belonging to chosen_stem
            if "_model_" in name or "_rank_" in name:
                if chosen_stem.split("_model_")[-1] not in name and \
                   chosen_stem.split("_rank_")[-1] not in name:
                    p.unlink(missing_ok=True)
    for d_name in ("processed", "lightning_logs", "msa"):
        d = out_dir / d_name
        if d.is_dir():
            shutil.rmtree(d, ignore_errors=True)

    return chosen
```

This **does not change any pipeline-visible behaviour** — `compute_metrics.py`
already only reads the sidecars belonging to the chosen model. It only
deletes files that the existing code was going to write to results/ but
never read.

A complementary change would add `workflow.onComplete` cleanup in
`main.nf`:

```groovy
workflow.onComplete {
    if (workflow.success) {
        // Best-effort: remove work/ and orphan caches once all stages
        // have completed and publishDir has finished copying.
        ['work', '.nextflow/cache', '.nextflow/plr'].each { p ->
            new File("${workflow.launchDir}/${p}").deleteDir()
        }
    }
}
```

Both are deferred to a separate task and not part of this audit's
deliverables.

## 11 — Validation after cleanup

For each cleaned run, confirm:

- `work/` does not exist; `results/` is unchanged in non-Boltz directories
  (`mpnn/`, `rfdiffusion/`, `orthogonal_metrics/`, `plots/`, etc.).
- One representative `cycle_0/initial/boltz_results_input/predictions/input/`
  has only `*_model_0.*` files left (4 PDBs and 4×4 sidecars per Boltz
  call should be gone).
- Top-level `cycle_0/initial_prediction*.pdb` and per-candidate
  `prediction.pdb` files are intact (the user-facing ChimeraX inputs).
- `cross_sequence_summary.csv` and
  `cross_sequence_summary_with_interface_metrics.csv` open and have the
  same row count as before cleanup.

## 12 — Open questions

1. **Should the forward-looking `_invoke_boltz` cleanup hook be applied
   before or after the next set of runs?** Doing it before would make
   future runs lighter from the start; doing it after means one more
   round of retroactive cleanup but no risk of changing behaviour
   mid-investigation.

2. **Are `.nextflow.log` files worth preserving long-term?** They're
   ~1 MB per run and useful for post-mortem of failures. Could be
   compressed to `.log.gz` (~100 KB) at cleanup time — small win, decide
   if worth the complexity.

3. **`workflow.onComplete` automation** — does the user want this added
   to `main.nf` so runs self-clean, or prefer the explicit external
   script? External script is auditable and easier to dry-run; automation
   is one less thing to remember.

4. **HADDOCK module** — `modules/haddock.nf` has 5 `publishDir mode: 'copy'`
   declarations. None of the current 8 eligible runs used the HADDOCK
   path (they all used `INPUT MODE 2: pre-docked complex PDB`). If
   future campaigns use HADDOCK, that subsystem may have its own
   bloat profile to audit separately.
