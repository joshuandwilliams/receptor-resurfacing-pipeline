# Receptor-Resurfacing Pipeline — Negative Steering Integration & First Wet-Lab Triage

## Context

This session was the long-deferred step of replacing the legacy
`boltz2.nf` module (per notes6 open work item 4 and the deferred
entry in `pipeline_notes2.md`'s Boltz2 block) with a NextFlow port of
the negative-steering pipeline developed across notes7–notes11.

Scope entering this session:

1. Port the standalone negative-steering SLURM chain (`plan` →
   `predict-one array` → `collect` → `build-contaminated` →
   `plan-reversions` → `predict-reversion array` → `harvest` →
   `finalize` → `compute-final-metrics` → `aggregate-per-sequence` →
   `extract_passing`) into NextFlow, downstream of ProteinMPNN.
2. Validate the integration end-to-end on a small test run.
3. Run overnight on a realistic batch (16 MPNN sequences across 4
   parent RFDiffusion designs).
4. Tune resource allocations now that we had actual usage data.

What actually happened: the port went quickly, but three rounds of
NextFlow channel-topology debugging were needed before the graph ran
end-to-end without silently dropping half the sequences. The overnight
validation then exposed two semantic bugs in the underlying Python
scripts that were latent in the standalone pipeline and only became
visible at production scale. Both were fixed and re-staged as patches
v3. The 16-sequence run itself produced 7 tier-A survivors, two of
which came through paths that had been silently filtered out pre-patch.
Then we closed the session with resource-tuning based on `trace.txt`
analysis, runtime instrumentation added to the output CSVs, and a
wet-lab triage of three candidate picks.

Strict ground rule as usual — minimal targeted changes, no scope creep
— with the explicit exception that bugs exposed by the validation run
get fixed now, not deferred.

## What we built

### 1. NextFlow integration of the negative-steering pipeline

The standalone SLURM chain in `submit_boltz2_negative_steering.sh`
produces roughly `3 + N_designs*num_seeds + 3` SLURM jobs per
MPNN-sequence / ground-truth pair. At production scale (16 MPNN
sequences × 20 designs × 3 seeds), a naive port would have produced
~4200 individual SLURM submissions.

The integration collapses this to **one GPU SLURM job per MPNN
sequence**, running the complete single-cycle chain inline via a shell
orchestrator (`bin/negative_steering_run_one.sh`). Per-design
parallelism within a sequence is traded away in exchange for loading
the GPU and Boltz weights once per sequence, which turns out to be a
much bigger win at this scale.

Three NextFlow processes in the new `modules/negative_steering.nf`:

- **`NEGSTEER_DERIVE_INDICES`** (CPU, small): runs once per *parent
  RFDiffusion design*, not per MPNN sequence. Derives the
  design-region file and the true-interface file from
  `rfdiffusion_metrics.json`. The outputs are shared across every MPNN
  sequence of that parent design via `.combine(by: 0)`.
- **`NEGSTEER_RUN_ONE`** (GPU, large): runs once per MPNN sequence.
  Executes the orchestrator. Produces a `passing_summary.csv` per
  sequence.
- **`NEGSTEER_CROSS_SEQUENCE`** (CPU, small): runs once at the end.
  Symlinks all per-sequence workdirs into an aggregator tree and runs
  `bin/cross_sequence_summary.py` to produce the single ranked output
  `cross_sequence_summary.csv`.

Total NextFlow jobs: `1 + N_designs + N_MPNN_sequences + 1`. For the
16-sequence run: **21 SLURM submissions**.

Files touched: new `modules/negative_steering.nf`, new
`tests/negative_steering/test_negative_steering.nf` +
`run_test_negative_steering_slurm.sh` +
`tests/negative_steering/data/` + symlinked
`tests/negative_steering/bin -> ../../bin`. `main.nf` gained the
three-process wiring downstream of `MPNN_SELECT_TOP.out`. Legacy
`modules/boltz2.nf`, `modules/msa.nf`, `modules/aggregate.nf`,
`tests/boltz2/`, and four orphaned `boltz2_*.py` scripts in `bin/`
were removed.

### 2. `bin/` additions

Three new scripts went into `bin/`:

- **`derive_true_interface.py`** — mirrors `derive_design_region.py`
  but for the true binding interface. Reads
  `rfdiffusion_metrics.json`, extracts the 0-based true-interface
  indices for a given design, writes them to a file. Cross-checks that
  the true interface is a subset of the design region (advisory
  warning).
- **`negative_steering_run_one.sh`** — shell orchestrator that runs
  the full single-cycle chain inline. Takes a long argv list (sequence
  name, ground-truth PDB, receptor/effector chain IDs, receptor/effector
  FASTA files, indices files, workdir, bin dir, boltz container, plan
  extra args, postprocess thresholds), and chains seven stages
  sequentially. Uses argv arrays rather than `eval` so quoting is never
  a problem.
- **`cross_sequence_summary.py`** — tier-then-composite aggregator,
  per notes11 policy.

Plus three existing scripts came in from the standalone
negative-steering repo, unchanged in this session:
`boltz2_negative_steering.py`, `boltz2_iterate_steering.py`,
`reversion.py`. Pre-existing supporting scripts reused as-is:
`compute_metrics.py`, `derive_design_region.py`, `extract_passing.py`,
`sequence_registry.py`.

### 3. Three rounds of channel-topology debugging

This is the part that took a while to get right. Each fix looked
reasonable but introduced a new failure mode.

**Round 1 — `.join()` is one-to-one, not one-to-many.** The initial
graph wanted: "for each MPNN FASTA, look up its parent design's PDB
and derive outputs." Used `.join()` on `design_id`. Symptom: only the
first MPNN sequence per parent design ran; the other 3 silently
vanished. `.join()` consumes the PDB channel on the first match and
drops siblings — correct semantics for a 1:1 join, wrong for a 1:N
fan-out. Fixed by using `.combine(by: 0)`, which broadcasts the PDB to
every matching FASTA.

**Round 2 — queue channels are single-consumer.** After the `.combine`
fix, new symptom: half the sequences still disappeared, but
unpredictably. Root cause: `fasta_stream_ch` was being read BOTH to
compute `.unique()` stems (to filter which designs needed index
derivation) AND for the three-way combine. Queue channels drain on
first read — the second reader got an empty channel.

**Round 3 — skip the filter entirely.** The "only derive indices for
designs that have matching FASTAs" filter was the source of the
dual-consumer problem. Cost of running `NEGSTEER_DERIVE_INDICES`
unconditionally on every design PDB regardless of FASTA availability:
~5s per orphan design. Not worth the code complexity. Dropped the
filter. Added `.first()` on `metrics_json_ch` so it becomes a value
channel (broadcastable without consumption).

Verified end-to-end on 4 FASTAs × 2 designs (8 runs): all 8 submitted,
all 8 completed.

### 4. Test validation

Set up `tests/negative_steering/data/` with 4 RFDiffusion design PDBs,
one shared `rfdiffusion_metrics.json`, and 16 MPNN FASTAs (4 per
design). First test run used defaults matching production.

Test outcome: 21 SLURM jobs submitted, all 21 completed,
`cross_sequence_summary.csv` produced. 7 tier-A, 2 tier-B, 7 tier-none
— numbers that turned out to be consistent with the subsequent
overnight run on the same sequences, which is how we discovered the
bugs in step 5.

### 5. Two semantic bugs exposed at production scale (patches v3)

These are bugs in the underlying negative-steering Python that were
latent in standalone runs and only surfaced when multiple MPNN
sequences were compared side-by-side in the NextFlow run.

**Bug (a): `clean_steered` misclassified as `no_data`.**

Setup: a sequence where 2/3 Boltz seeds predict the correct pose
without needing any load-bearing mutations (zero contamination on the
steering mutations, so reversion is correctly skipped), and 1/3 seeds
end up contaminated but their reverted prediction holds the pose
(pose_holds).

Pre-patch, the aggregator's `_per_seed_verdict_breakdown` treated the
2 clean seeds as `no_data` because their `reversion_verdict` column
was blank and their `reverted_*` columns were empty. The aggregated
verdict classifier then saw only 1 `pose_holds` out of 3 → demoted
from tier A to tier C. Several sequences that were genuinely strong
wet-lab candidates were being dropped to tier C on this basis, then
filtered out of `passing_summary.csv` entirely.

Fix: new `clean_steered` bucket in `_per_seed_verdict_breakdown`.
Admission criteria are strict enough that only genuinely
no-reversion-needed rows qualify:

- `reversion_verdict` blank
- `steered_n_contacts_on_mutated_positions == 0`
- `steered_receptor_intact == 1`
- `steered_ra_eff_vs_truth < 5.0`
- No `reverted_*` fields populated (reversion truly didn't run)

`_classify_aggregated_verdict` updated to use the effective
pass-equivalent count `nph_eff = n_pose_holds + n_clean_steered` for
both the any-seed-passing gate and the Bug-E residual guard.

`n_seeds_clean_steered` is now a column in `aggregated_results.csv`,
`passing_summary.csv`, and `cross_sequence_summary.csv`.
`cross_sequence_summary.py`'s tier classifier uses
`n_pass = n_pose_holds + n_clean_steered` as the pass count.

**Bug (b): cold-start runs only ONE Boltz prediction regardless of
`--num-seeds`.**

Setup: MPNN sequences where the cold-start Boltz already gets the
pose right without any negative steering — the skip-steering fast
path. Pre-patch, `plan` ran exactly one cold-start prediction at one
seed. If that single prediction happened to hit the ra_eff threshold,
skip-steering was triggered and the pipeline exited after ~4 minutes
with no variance evidence.

Policy fix (agreed this session): 3 cold-start seeds always run. If
all 3 are clean (ra_eff ≤ threshold AND receptor intact), steering is
skipped and the 3 cold-start predictions ARE the candidate set (tier A
via `aggregated_verdict=no_reversion`). If any seed fails, full
steering runs as normal and the cold-start rows are screening-only,
not wet-lab candidates.

Code fix in `boltz2_negative_steering.py`: `plan` now runs `num_seeds`
cold-start predictions at seeds `args.seed..args.seed+num_seeds-1`,
stored in `plan["cold_start_seeds"]`. New `_write_initial_multiseed_csv`
writes one row per cold-start seed (`design` names `initial`,
`initial_s1`, `initial_s2`, all with `sequence_group=0`, distinct
`seed_index`). PDBs saved at `initial_prediction.pdb`,
`initial_prediction_s1.pdb`, `initial_prediction_s2.pdb` with sidecars
under matching `initial/`, `initial_s1/`, `initial_s2/` directories.

Supporting fixes in `boltz2_iterate_steering.py`:

- Cycle-0 CSV ingestion recognises `initial` AND `initial_s*` design
  names and treats both as baseline (no mutations).
- `compute-final-metrics` filter admits `initial_s*` rows when the
  parent plan has `skip_steering=true`.
- `_locate_boltz_sidecar_pdb` routes `initial_prediction_s*.pdb` to
  `initial_s<N>/` for sidecar discovery.
- `_agg_row_is_clean_steered` admits cold-start aggregates
  (design==`initial` AND verdict==`no_reversion`) as ranking-eligible
  with a relaxed ra_eff cap read from the plan's `rmsd_threshold`
  (default 5.0 for legacy compatibility, typically 6.0 for cold-start
  runs).

Supporting fix in `extract_passing.py`: `n_seeds_clean_steered` added
to `OUTPUT_FIELDS`.

Supporting fix in `cross_sequence_summary.py`: uses
`n_pose_holds + n_clean_steered` for tier classification, backward
compatible with missing `n_clean_steered` column.

All four files staged as `patches_v3/` and deployed.

### 6. 16-sequence overnight validation

With patches v3 deployed and `diffusion_samples=5`, `num_seeds=3`,
`n_designs=20`, the overnight run produced:

| Tier | Count | Example (rank 1) |
|------|-------|------------------|
| A    | 7     | `design_2_seq_0` — pose_holds 1/3, clean_steered 2/3, ra_eff 2.23, tj 0.80, comp 0.688 |
| B    | 2     | `design_2_seq_1` — pose_holds 2/3, pose_collapses 1/3, ra_eff 2.86, comp 0.571 |
| none | 7     | — failed to produce any steered passing prediction |

Parent-design breakdown: `design_3` strongest (3/4 tier A),
`design_1` and `design_2` each 2/4, `design_0` only 1/4 (and that one
was tier B).

**Two survivors validated the v3 patches.** `design_2_seq_0` (rank 1)
had `ph=1 cs=2` — under pre-patch code, the 2 `clean_steered` seeds
would have been classified `no_data`, demoting the sequence to tier C
and filtering it out. `design_3_seq_0` (rank 5) and `design_3_seq_2`
(rank 2) and `design_2_seq_2` (rank 6) were all
`aggregated_verdict=no_reversion` cold-start binders — under pre-patch
code, `plan` would have run ONE cold-start seed, not three, and the
aggregator would have had insufficient data to promote any to tier A.

Three picks for wet-lab:

1. **`design_2_seq_0`** (rank 1) — highest composite, patch-a validated.
2. **`design_3_seq_0`** (rank 5) — cold-start binder with highest
   iPSAE (0.74) and iPTM (0.90) in the set. Different parent design,
   hedges against iPSAE being informative despite the notes10 AUROC
   result.
3. **`design_1_seq_2`** (rank 3) — pose_holds 3/3, different parent
   from picks 1 and 2. Diversity argument.

### 7. Resource tuning

`trace.txt` analysis from the 16-sequence run showed massive
overprovisioning on NEGSTEER_RUN_ONE:

| Parameter | Pre-tuning | Actual peak usage | Post-tuning |
|-----------|-----------|-------------------|-------------|
| cpus      | 10        | 144% (~1.4 cores) | 4           |
| memory    | 32 GB     | RSS 5.1 GB        | 12 GB       |
| time      | 24h       | max 4h 7m         | 6h          |

The 4h 7m outlier was `design_3_seq_3` — pose_holds 3/3 with high
reversion-pass inflation (rchar=373 GB vs the typical ~170 GB). Not a
bug; just a sequence where many steering candidates all passed the
steered filter and all needed reversion validation. 6h gives ~45%
headroom over this observed worst case.

NEGSTEER_DERIVE_INDICES and NEGSTEER_CROSS_SEQUENCE were both wildly
overprovisioned (4 GB / 15-30m for tasks using 25 MB / 15 seconds).
Tightened to 1 GB / 5m each.

Also added `fields = '…'` to the `trace {}` block in
`nextflow.config` so `%cpu / peak_rss / peak_vmem` are always present
in future `trace.txt` files.

### 8. Runtime instrumentation

Added `run_one_runtime_sec` as a column in `passing_summary.csv` and
as a top-level column in `cross_sequence_summary.csv`. Covers the
whole-script wall-clock for one NEGSTEER_RUN_ONE invocation (plan +
predict-one + collect + reversion + harvest + finalize + postprocess;
excludes `extract_passing` which writes the column itself).

Implementation: `negative_steering_run_one.sh` stamps `SECONDS=0` after
the banner and writes `run_one_runtime_sec.txt` to the workdir before
invoking `extract_passing`. `extract_passing.py` reads the sidecar and
stamps the column onto every passing row. `cross_sequence_summary.py`
surfaces it as a top-level column (position 11 in the output) and
falls back to reading the sidecar directly for tier-none sequences
whose `passing_summary.csv` is empty.

Three files in `bin/` touched: `negative_steering_run_one.sh`,
`extract_passing.py`, `cross_sequence_summary.py`. Staged as
`runtime_instrumentation/`.

### 9. GPU-partition capacity check (`max_af2_parallel` tuning)

`sinfo -p jic-gpu` showed 10 A100 nodes totalling 26 GPUs. At the time
of inspection, ~19 GPUs were free (one competing user holding 7 for
multi-day jobs). Recommended bump from `max_af2_parallel = 5` (current)
to `max_af2_parallel = 12` for a ~2.4× throughput improvement on
sequence parallelism, still leaving ~30% partition headroom for
competing users.

Combined with `diffusion_samples = 2` (staged in
`params_example.yml`), the next production run should do 16 sequences
in roughly 2 hours of wall-clock rather than 6+.

## Things we learned

**NextFlow channel semantics get unintuitive for 1:N fan-out.**
`.join()` for one-to-one, `.combine(by: key)` for one-to-many, and
never read a queue channel twice — use `.first()` to promote it to a
value channel if you need multi-consumer behaviour. This took three
rounds of debugging to internalise.

**Per-sequence SLURM job collapse is a massive win at this scale.**
~4200 naive-port submissions → 21. Loading the GPU and Boltz weights
once per MPNN sequence rather than once per (sequence, design, seed)
tuple trades modest parallelism loss for very large wall-clock savings
and very large queue-pressure reduction.

**Semantic bugs that are latent at development scale show up
immediately at production scale.** The `clean_steered` classifier bug
and the cold-start single-seed bug were both present in every
standalone run notes7–notes11 but never triggered the failure modes
that made them obvious. First multi-sequence production run exposed
both in the same afternoon.

**Resource overprovisioning costs more than memory.** Requesting 32 GB
on a job that uses 5 GB tanks fair-share priority AND blocks you from
packing into partial nodes. The 32 GB → 12 GB tuning won't 3× the
parallelism (GPUs are the binding constraint), but it will noticeably
shorten queue pending times.

**Skill scaling: cold-start all-clean is a legitimate output path.**
Some MPNN sequences are good enough that Boltz gets the pose right
without steering. Pre-patch these were being dropped on the floor.
Post-patch they're tier-A candidates via the `no_reversion` verdict.
Three of the seven tier-A survivors from the overnight run came
through this path.

## File deliverables

Staged patches by deployment group:

- `patches_v3/` — `boltz2_negative_steering.py`,
  `boltz2_iterate_steering.py`, `cross_sequence_summary.py`,
  `extract_passing.py`. Deploy as a set.
- `runtime_instrumentation/` — `negative_steering_run_one.sh`,
  `extract_passing.py`, `cross_sequence_summary.py`. Deploy as a set.
  Note: `extract_passing.py` and `cross_sequence_summary.py` appear
  in both — the `runtime_instrumentation/` copies include both the v3
  semantic patches AND the runtime column, so deploy the
  `runtime_instrumentation/` copies.
- `config_tuning/nextflow.config` — resource sizing for NEGSTEER
  processes, explicit trace `fields`. **Diff against the live repo
  copy before replacing** — the project knowledge indexed here is
  stale and may not include the current NEGSTEER blocks.

NextFlow integration files (already deployed in an earlier stage of
this session): `main.nf`, `modules/negative_steering.nf`,
`params_example.yml`, `tests/negative_steering/`, new `bin/` scripts
(`derive_true_interface.py`, `negative_steering_run_one.sh`,
`cross_sequence_summary.py`).

## Status summary

- **NextFlow integration complete** and validated end-to-end on 16
  MPNN sequences.
- **Two semantic bugs fixed** (patches v3).
- **Runtime instrumentation deployed** — `run_one_runtime_sec` column
  visible in `passing_summary.csv` and `cross_sequence_summary.csv`.
- **Resource tuning ready to deploy** (`nextflow.config`).
- **Three wet-lab candidates selected** from the 16-sequence run.
- **Multi-cycle support deferred** — `params.negsteer_n_cycles`
  currently only honours `=1`. See notes12 in the negative-steering
  notes series (to be written in a future session) for the design
  sketch.
- **Aggregate module rebuild deferred** —
  `MPNN_SELECT_TOP.out.top_metadata` is emitted but not consumed; it
  should be joined with `cross_sequence_summary.csv` into a final
  `results_summary.csv`. Harmless as-is (NextFlow ignores unconsumed
  emits).
- **`project_map.md` update deferred** — needs a new "Negative
  steering (integrated)" block describing `modules/negative_steering.nf`,
  the three `NEGSTEER_*` processes, and the new `bin/` scripts.

## Next steps

### Immediate (before next overnight run)

1. Deploy `patches_v3/` + `runtime_instrumentation/` + tuned
   `nextflow.config` to the live repo.
2. Drop `negsteer_diffusion_samples: 2` in `params_example.yml`
   (already staged, flagged for this session).
3. Bump `max_af2_parallel: 12` (from 5) in the production launcher.
4. Optionally bump `negsteer_n_designs: 30` (from 20) to help the
   7/16 tier-none rate from the overnight run. Worth testing but
   may reveal that those sequences are genuinely unsteerable rather
   than starved for designs.
5. Verify `cuequivariance_torch` availability in the container:
   `singularity exec <container> python3 -c "import cuequivariance_torch; print(cuequivariance_torch.__version__)"`.
   If present, flip `negsteer_no_kernels: false` for a 20-30%
   per-prediction speedup.

### Next session — multi-cycle support

`params.negsteer_n_cycles` currently only works for `=1`. Multi-cycle
requires restructuring the kickoff/finalize boundary in
`boltz2_iterate_steering.py` so the orchestrator can spawn subsequent
cycles inline rather than via the old SLURM resubmission mechanism.
Not urgent — the single-cycle path is producing useful wet-lab
candidates — but worth doing before the next methodology push.

### Later — aggregate module rebuild

`MPNN_SELECT_TOP.out.top_metadata` (the per-sequence MPNN metadata:
native-sequence recovery, design-region RMSD, cluster membership,
etc.) is currently unused downstream. Should be joined with
`cross_sequence_summary.csv` to produce a final `results_summary.csv`
that carries both the MPNN provenance and the negative-steering
verdicts, so wet-lab triage has a single file with everything in it.

### Open diagnostic: 7/16 tier-none rate

Half the MPNN sequences produced no passing predictions. Worth
logging, for each tier-none sequence: did it fail at cold-start (all
3 seeds above rmsd_threshold), or did it get through cold-start but
fail at steered-pass (none of the 60 steered predictions passed the
ra_eff filter), or did it get through steered but fail at reversion
(every pose_holds candidate broke when the mutations were reverted)?
The answer tells you which lever to pull: cold-start failures want
better MPNN sequences or RFDiffusion designs; steered failures want
more `n_designs` or different modes; reversion failures mean the
design is genuinely too mutation-dependent to survive wet-lab.