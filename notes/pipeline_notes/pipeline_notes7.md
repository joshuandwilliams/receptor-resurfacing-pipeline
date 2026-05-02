# pipeline_notes6.md — Kernel A/B result, per-prediction cost breakdown

## Context

Container `boltz2_negsteer.img` was rebuilt in the previous session (see pipeline_notes5.md) with NGC-style CUDA installation and all pip-nvidia-cu12 wheels purged. Rebuild validated: `cuequivariance_ops_torch` imports, matmul runs on GPU, smoke tests pass. Before committing to `negsteer_no_kernels = false` in production, wanted to measure whether the kernels actually deliver the hoped-for 20–30% speedup on this workload.

## The A/B experiment

Three runs of `test_negative_steering.nf` with `n_designs=3`, `num_seeds=3`, `recycling_steps=3`, `mild` mode:

| Run | diffusion_samples | kernels | purpose |
|---|---|---|---|
| A | 5 | ON | production-equivalent, kernels enabled |
| B | 5 | OFF | production-equivalent, kernels disabled |
| C | 2 | OFF | quantify diffusion vs fixed cost breakdown |

Timings read from each sequence's `run_one_runtime_sec.txt`, filtered to the 0-reversion path (the 3-cold-start + 9-steered = 12-prediction path) for comparability. Cold-start-only sequences (skip_steering path, 3 predictions) also recorded.

### Raw per-sequence wall times (seconds)

| sequence | A (d=5 ON) | B (d=5 OFF) | C (d=2 OFF) |
|---|---|---|---|
| design_0_seq_0 | 1358 (3 rev) | 1034 (0 rev) | 852 |
| design_0_seq_1 | 1084 | 1047 | 863 |
| design_0_seq_2 | 1080 | 1029 | 860 |
| design_1_seq_0 | 262 (cold-only) | 227 (cold-only) | 191 (cold-only) |
| design_1_seq_1 | 1088 | 1018 | 860 |

design_0_seq_0 in run A ran 3 reversions (9 extra predictions) and is excluded from the 0-reversion comparison.

### Per-prediction time, 0-reversion path (12 predictions/sequence)

| | A (d=5 ON) | B (d=5 OFF) | C (d=2 OFF) |
|---|---|---|---|
| mean sequence wall time | 1085 s | 1031 s | 861 s |
| per prediction | 90.4 s | 85.9 s | 71.8 s |

## Result 1: kernels are net SLOWER on this workload

Kernels ON vs kernels OFF at d=5: **+54 s per sequence, +4.5 s per prediction, ~5% slower with kernels on.**

Consistent across all four matched 0-reversion sequences (not noise). Code path confirmed clean — `boltz/model/layers/triangular_mult.py` and `triangular_attention/primitives.py` call `cuequivariance_torch.primitives.triangle.{triangle_multiplicative_update, triangle_attention}` directly under `use_kernels=True`, no try/except fallback. Triton cache populated on disk during runs, confirming kernels actually JIT-compile and execute.

Why they're slower isn't definitively traced, but most likely: cuEquivariance kernels are tuned for AlphaFold-scale tensor shapes (long MSA × large pair representation). Our complexes are small (~75-residue receptor + ~80-residue effector). At that size, kernel dispatch and memory-layout overhead likely exceed the arithmetic win.

**Decision: keep `negsteer_no_kernels = true` in production.** The kernel path works, it just isn't faster for our inputs. Rebuilt container is still valuable for unrelated reasons (modernized CUDA stack, cleaner dep management, the rebuild bug lessons). Keeping the fallback costs nothing at runtime.

## Result 2: where prediction time actually goes

Comparing B (d=5 OFF, 85.9 s/prediction) vs C (d=2 OFF, 71.8 s/prediction): dropping 3 diffusion samples saves 14.1 s per prediction, or **~4.7 s per diffusion sample**.

Extrapolating: at d=5, diffusion cost is 5 × 4.7 = 23.5 s per prediction. Total is 85.9 s. So diffusion is **~27% of per-prediction wall time**. The remaining ~73% (~62 s) is fixed per Boltz call:

- Model weight load (GPU cold-start)
- Input featurization (tokenization, template embedding, MSA processing if any)
- Pairformer trunk forward pass (runs ONCE per call, independent of diffusion_samples)
- Confidence head
- Output writing (PDB, CIF, JSON, full PAE)
- Singularity exec overhead, Python import chain

This explains why `--diffusion_samples 5 → 2` only speeds things up ~16% end-to-end, not 60%: only the diffusion denoiser scales with sample count. Everything else is sunk cost per `boltz predict` invocation.

## Implications for 10k-prediction experiments

Current per-prediction time at production config (d=5, 3 recycling, kernels OFF): ~86 s. At 10k predictions: ~240 GPU-hours per experiment.

The ~60 s fixed overhead per call is where the biggest theoretical speedup lies, NOT the diffusion loop. Possible paths forward (filed as future task):

1. **Batch multiple inputs through one `boltz predict` invocation.** If it accepts a directory of YAMLs, weight-load cost amortizes across all of them. Could save 5-15 s per prediction depending on cold-start cost. At 10k predictions × 10 s saved = ~28 GPU-hours per experiment.

2. **Persistent Python process holding the model.** Invoke Boltz's inference API directly in a loop instead of spawning `boltz predict` per prediction. Larger engineering lift, larger potential win — model load only once per 10k calls.

3. **Lower `recycling_steps` from 3 to 1.** Pairformer trunk runs 3× under current config. Dropping to 1 would meaningfully speed up the trunk, at the cost of prediction quality.

4. **`torch.compile()` on the trunk.** If Boltz supports it, 1.3–1.5× speedup on the trunk after a one-time warmup.

Options 1 and 4 are worth investigating without changing prediction quality. Option 3 is a quality trade-off that should be benchmarked against ground truth before adoption. Option 2 is the highest-lift / highest-reward.

## Misleading baseline from earlier notes

The todo_list4 acceptance criteria quoted `notes13 baseline (median 105 min with kernels off)` and targeted `~75 min with kernels on`. That target was based on a hope/estimate from cuEquivariance documentation, not a measured result on this pipeline. The measured result is the opposite direction: kernels are ~5% slower, not 20–30% faster. The 105 min baseline was also at different parameters (overnight 16-sequence production run, not 3-design test). Neither the target nor the baseline is directly comparable to the clean A/B above.

## Operational notes

- `run_one_runtime_sec.txt` is written by `negative_steering_run_one.sh` BEFORE the final `extract_passing` step, so it's readable from the work dir as soon as the Boltz + metrics work is done — no need to wait for the full nextflow pipeline to complete. Handy for mid-run diagnosis.
- Boltz predictions with different reversion counts aren't directly comparable; always filter to matched reversion counts when doing wall-time comparisons.
- Each NEGSTEER_RUN_ONE task gets one A100 exclusively via `SLURM_JOB_GPUS`; `CUDA_VISIBLE_DEVICES=0` is always set within the cgroup but refers to a different physical GPU per job. When debugging apparent contention, always check `SLURM_JOB_GPUS`, not `CUDA_VISIBLE_DEVICES`.
- `nvidia-smi` from a plain SSH session to a GPU node may not show processes or utilization because of cgroup isolation. Use `srun --jobid=<id> --overlap --pty nvidia-smi` to see what a running job sees.

## Status on the two fixes baked into the rebuilt container

1. **Stale cuBLAS preload fix** (pip nvidia-*-cu12 wheel purge): validated. Kernels load, GPU matmuls succeed. Not actually exercised in production since we're leaving `no_kernels=true`, but the container is clean regardless.

2. **pynvml FutureWarning fix** (uninstall deprecated pynvml PyPI package): validated. No more warning at torch import.

## Files affected this session

- `/mnt/project/pipeline_notes6.md` — this document
- No code changes; production config keeps `negsteer_no_kernels = true`

## Open questions → future tasks

- How much of the 62 s fixed overhead per `boltz predict` is model load vs. CPU-side setup vs. I/O? Diagnose with `strace -c` or a timed wrapper around `boltz predict` before deciding which batching strategy (1 vs 2 above) is worth pursuing.
- Does `boltz predict` accept a directory of YAMLs? Check `boltz predict --help` and `main.py` around line 1000+. If yes, this is the cheapest win available.
- Does cuequivariance kernel slowdown persist on larger complexes (e.g. multi-chain / longer receptors)? Only worth re-testing if the pipeline starts processing substantially larger substrates.
