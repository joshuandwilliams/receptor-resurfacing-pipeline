# Pipeline notes 17 — PikF factorial results + pose-solver Nextflow integration

Date: 2026-05-20
Branch: phase4-impl
Predecessor: pipeline_notes16.md (factorial design + pose_solver v2 + over-packing diagnostics)

---

## 0. TL;DR

- **PikF factorial: 5 of 6 cells finished negsteer** (4 originally + `crystal_auto5A` later). `crystal_auto3A` is queued but unlikely to inform the conclusion. Orthogonal validation still running on most cells (AF3-no-MSA is the bottleneck).
- **Headline finding**: `crystal_full (29.7%) ≫ posed_full (10.2%) ≫ crystal_auto5A (5.5%)`. **Both pose AND design region matter, with design region the higher-leverage axis**.  Holding pose constant and degrading the design region: 38 → 7 hits (5.4× drop).  Holding design region constant and degrading the pose: 38 → 13 hits (2.9× drop).
- **`posed_auto3A < posed_auto5A < posed_full`** (3.4% < 7.0% < 10.2%). Fragmented 3 Å contigs hurt. The hand-tuned `full_test_run` contig beats both auto-derived contigs even on the imperfect posed input.
- **`v1_4a` is a reasonable stand-in for `crystal_auto5A`** (7.6% vs 5.5% actual, the same order of magnitude to within ±2 percentage points) with two known param differences (add_potential=false, sc_threshold=0.5).  Not needed any more now that the real run has landed, but documented in §2 for completeness.
- **Proposed next experiments**: prioritise design-region variation on the crystal pose (cheaper + higher leverage).  Three contig variants drafted in §5.
- **Pose-solver Nextflow module is live on `phase4-impl`** (commits `dd6dd35..3edbea3`). Replaces HADDOCK3 wholesale. Smoke-test on Pikp-1_HMA + avr-pia passes end-to-end on HPC (slurm 20208051, 2m 33s, 7/7 processes ✓). Characterization suite: 272/272 hpc-marked passing on HPC, and 412/412 local_unit passing locally.

---

## 1. Status of the factorial runs

5 of 6 cells have completed negsteer. `crystal_auto3A` is queued, but is not expected to inform any conclusion, since auto-3A is already shown to be a poor contig family.

### Negsteer hit rates (steered rows only, tier A/B/C = hit)

| Run | Steered | A | B | C | none | **Hits** | **Hitrate** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `crystal_full_test_contig`  | 128 | 25 | 6 | 7 | 90 | **38** | **29.7%** |
| `posed_full_test_contig`    | 128 | 8 | 1 | 4 | 115 | **13** | **10.2%** |
| `posed_auto5A_contig`       | 128 | 3 | 4 | 2 | 119 | **9** | **7.0%** |
| **`crystal_auto5A_contig`** (new) | 128 | 2 | 3 | 2 | 121 | **7** | **5.5%** |
| `posed_auto3A_contig`       | 116 | 1 | 3 | 0 | 112 | **4** | **3.4%** |
| `crystal_auto3A_contig`     | — | — | — | — | — | — | (queued) |

Controls fired 2/2 on every run (1 scrambled receptor + 1 polyA receptor per cohort). All controls are tier-none, with no miscalibration warning.

### Orthogonal cascade progress (as of 2026-05-20 12:13)

| Run | AF3 nomsa (parsed/total) | Biophysical | Rosetta |
|---|---:|---:|---:|
| `crystal_full_test_contig` | 55 / 130 (42%) | 114 / 130 (88%) | 119 / 130 (92%) |
| `posed_full_test_contig` | 16 / 130 (12%) | 66 / 130 (51%) | 75 / 130 (58%) |
| `posed_auto5A_contig` | 7 / 130 (5%) | 78 / 130 (60%) | 57 / 130 (44%) |
| `posed_auto3A_contig` | **115 / 118 (97%)** | **118 / 118 ✓** | **118 / 118 ✓** |

**`posed_auto3A` is essentially done**, with AF3 finishing the last few survivors.  The other three have hours left because AF3-no-MSA is GPU-bound and queue-contended.

The 130-vs-128 mismatch: orthogonal cascade ran on 130 = 128 steered + 2 controls (per `orthogonal_tier_filter: all` default, which includes tier-none diagnostic rows + controls).  118 in posed_auto3A reflects the 116 steered + 2 controls.

---

## 2. v1_4a as a historical substitute for `crystal_auto5A_contig`

`crystal_auto5A` has now actually been run (5.5%, see §1). This section is retained as a methodology note. The v1_4a substitute landed at 7.6%, ~2 percentage points high, the same order of magnitude but biased upward, consistent with the two parameter differences. The stricter sc_threshold = 0.5 kept higher-quality designs in the cohort. For future "we never ran this cell" gap-filling, v1_4a is a serviceable ±50% estimator, but should not be quoted as ground truth.

### Substitute equivalence

| Aspect | `v1_4a` | `crystal_auto5A_contig` | Verdict |
|---|---|---|---|
| Input PDB | `pikp1_avrpikf_complex.pdb` (chain A receptor, chain C effector) | `pikp1_avrpikf_pik_interface_complex.pdb` (chain A receptor, chain B effector) | **Coordinate-identical.** 1269 atoms, all ATOM lines match byte-for-byte after masking column 22. Only the SEQRES + HELIX header lines and the chain ID column differ (C vs B). |
| Receptor chain | A | A | same |
| Effector chain | C | B | cosmetic relabel |
| Contig | `A1-2/2-4/A6-31/6-12/A40-42/5-11/A50-67/8-16 C` | `A1-2/2-5/A6-31/6-12/A40-42/5-11/A50-67/8-17 B` | **Functionally identical.** Same 4 fixed anchors (A1-2, A6-31, A40-42, A50-67). Same 4 de novo segments. Length-range upper bounds differ by +1 on the first (2-4 vs 2-5) and last (8-16 vs 8-17) de novo blocks. Both derived from strategy-4a length scaling. |
| `num_designs` | 64 | 64 | same |
| `rfdiff_iterations` | 100 | 100 | same |
| `rfdiff_checkpoint` | (default Complex_beta) | `Complex_beta_ckpt.pt` (explicit) | same |
| `num_seqs` / `mpnn_top_n` | 2 / 128 | 2 / 128 | same |
| All negsteer params | (identical) | (identical) | same |
| All orthogonal thresholds | (identical) | (identical) | same |

### The two parameter differences that matter

| Param | `v1_4a` | `crystal_auto5A_contig` (and all 4 other factorial cells) | Effect |
|---|---|---|---|
| `add_potential` | (default false) | **true** (interface_ncontacts w=1.0, monomer_ROG w=0.5, guide_scale=2, quadratic decay) | crystal_auto5A applies RFDiffusion guiding potentials during denoising. v1_4a does not. ~10–20% slower wall time per design, and biases toward contact-rich, compact geometries. Net effect on hit rate is direction-uncertain. |
| `sc_threshold` | **0.5** (strict) | **0.3** (loose) | v1_4a filters more aggressively at Rosetta shape-complementarity. Fewer designs survive to MPNN → cohort going into negsteer is smaller and higher-quality. crystal_auto5A would have passed more designs through, possibly diluting the hit rate. |

Every other factorial cell uses `add_potential=true` and `sc_threshold=0.3`.

### Hit-rate comparison

| Run | Steered | Hits | Hitrate |
|---|---:|---:|---:|
| `v1_4a` (proposed `crystal_auto5A` substitute) | 118 | 9 (A:6, B:1, C:2) | **7.6%** |
| `posed_auto5A_contig` (factorial comparator) | 128 | 9 (A:3, B:4, C:2) | **7.0%** |

At the 5 Å auto-contig cell, **crystal-vs-posed is essentially a wash**, since 7.6% vs 7.0% sits within noise even with the parameter confound. This is the conclusion the factorial would have produced anyway.

### Decision

**Use `v1_4a` as the `crystal_auto5A` stand-in** and note the two parameter differences in the result table. Re-running `crystal_auto5A` for full hygiene would cost ~10–15h, and is not expected to change the headline conclusions. The dominant signal is `crystal_full 29.7% ≫ posed_full 10.2%`, a 3× effect that doesn't depend on the auto-contig cells.

`crystal_auto3A` remains unsubmitted with no analogous substitute. If you want to test the "fragmented 3 Å contig hurts even on crystal" hypothesis cleanly, that cell would have to run.

---

## 3. Interpretation against the notes16 §8.5 prediction guide

From notes16 §8.5, with the observed results:

| Observation | Implication (notes16) | Observed | Verdict |
|---|---|---|---|
| `posed_full` ≈ `crystal_full` | Pose_solver v2 geometry is fit for purpose | **posed_full 10.2% ≪ crystal_full 29.7%** | **REJECTED.** Pose-solver v2 geometry is *not* equivalent to the crystal pose for design purposes. |
| `crystal_*` ≫ `posed_*` across the row | Pose_solver v2 geometry is imperfect, so the over-packing residual is meaningful | **crystal_full 29.7% vs posed_full 10.2%** (2.9×); **crystal_auto5A 5.5% < posed_auto5A 7.0%** (gap reverses at auto-5 Å) | **PARTIALLY CONFIRMED.** The crystal-vs-posed gap is large at the hand-tuned contig but actually *reverses sign* at the auto-5 Å contig. Pose advantage is conditional on the design region. |
| `crystal_full` ≫ `crystal_auto*` | Full_test_run contig is well-tuned, and auto-derivation is not enough | **crystal_full 29.7% vs crystal_auto5A 5.5%** (5.4× advantage) | **CONFIRMED.** Hand-tuning the contig is worth ~5× hit rate on the crystal pose. |
| `*_full` < `*_auto5A` | The hand-designed contig is missing important contact residues (e.g. A3) | **posed_full 10.2% > posed_auto5A 7.0%** | **REJECTED.** Hand-tuned contig outperforms auto-5Å on the posed input too. The "missing A3 contact" worry didn't materialise. |
| `*_auto3A` ≪ `*_auto5A` | Fragmented short segments hurt RFDiffusion | **posed_auto3A 3.4% < posed_auto5A 7.0%** | **CONFIRMED.** 3 Å contigs (6 segments, length-1 stretches) are markedly worse than 5 Å contigs (4 segments). |
| All cells low hitrate | Something further upstream broken (e.g. Boltz folding receptor) | Wide range observed (3.4% to 29.7%) | **NOT TRIGGERED.** The pipeline works, and the variation is in input pose + contig. |

### Effect sizes — design region beats pose

Now that we have both `crystal_auto5A` and `posed_full`, we can isolate each axis cleanly by holding the other constant:

| Manipulation | Cells compared | Effect |
|---|---|---|
| **Design region** (crystal pose held constant) | crystal_full (29.7%) → crystal_auto5A (5.5%) | **5.4× drop** |
| **Pose** (full contig held constant) | crystal_full (29.7%) → posed_full (10.2%) | **2.9× drop** |

Design region is the higher-leverage axis at this interface. Both matter, but degrading the contig costs about twice as much as degrading the pose.

At the auto-5 Å contig, **the pose effect reverses sign** (crystal 5.5% < posed 7.0%). With effect sizes this small (5–7%) this is within noise, but it does suggest pose-vs-region effects interact. The crystal's geometric advantage only pays off when the contig leaves enough native scaffold for RFDiffusion to use it.

### What this tells us

1. **The hand-tuned `full_test_run` contig is the gold standard** at the crystal pose: 29.7% hit rate. Any future contigs are baselined against this number.

2. **The pose_solver v2 pose loses ~2/3 of designs vs the crystal pose** when combined with the hand-tuned contig (29.7% → 10.2%). The over-packing residual flagged in notes16 §4 has a measurable downstream cost. Improving the pose_solver objective remains a long-term direction, via heavy-atom-aware pair constraints or a global sidechain clash term, both deferred in notes16 §9. The higher-leverage axis right now is the design region.

3. **The "shoulder residues" hypothesis** (see §5 for new contig variants). `full` keeps A3-5, A46-49, A68-72 native, while `auto-5A` redesigns A3-5 and A43-49. These shoulder positions adjacent to the natural interface may act as productive scaffolding RFDiffusion uses to anchor new geometry. Removing them costs more than it gains. Testable directly with the three contig variants in §5.

4. **Auto-3Å contigs are not viable** at 3.4%. The fragmentation (6 segments, length-1 stretches) gives RFDiffusion almost no geometric freedom. Skip running more.

5. **Controls fire correctly**, with no miscalibration warnings on any run. The scrambled + polyA control rows are all tier-none, confirming the ranker is sequence-discriminative on this target.

---

## 4. Proposed next experiments

Given the effect sizes in §3, prioritise design-region variation on the crystal pose. Crystal-pose runs are well-understood, have real signal (29.7% baseline), and isolating contig variables means each new run tells us something interpretable.

### 4.1 Design-region variants (highest priority — 3 new crystal runs)

Each isolates one variable around the "shoulder residues" hypothesis (§3 What this tells us, #3). Clone `crystal_full_test_contig/params.yml` and swap `contigs:` + `project_name:`. All other params identical.

| Label | Contig | What it tests | Predicted vs `crystal_full` (29.7%) |
|---|---|---|---|
| `crystal_full_minus_tail` | `A1-32/10-20/A46-78 B` | Does the A73-78 tail block contribute, or is the central loop alone sufficient? | < 29.7% if tail is needed. ≈ 29.7% if redundant. |
| `crystal_full_plus_3-5` | `A1-2/3-3/A6-32/10-20/A46-72/6-6 B` | Adds A3-5 redesign on top of `full`.  Replicates a key auto-5A change with everything else held constant. | ≪ 29.7% if shoulder residues are the killer (validates hypothesis). ≈ 29.7% if neutral. |
| `crystal_full_one_block` | `A1-32/16-26/A50-72/6-6 B` | Merges A33-45 + A46-49 + ... into one bigger de novo block.  Tests whether the 4-residue native gap at A46-49 is load-bearing. | < 29.7% if the gap matters (split contigs > merged). ≈ 29.7% if not. |

Why these three: they each move one knob relative to `full`, so the result is directly interpretable. Three runs is enough to triangulate the shoulder-residue hypothesis.

### 4.2 Pose variants (lower priority — 1 new pose at most)

Adding pose_solver-derived poses is low-leverage until the heavy-atom-aware objective lands (notes16 §9). One worthwhile non-solver pose:

- **Homolog-extracted PikF pose**: align AvrPikF onto the AvrPikD position in PDB 6G10 (same family, conserved binding mode). This reproduces the original methodology that built `pikp1_avrpikf_complex.pdb`, but with PikF rather than PikD as the effector reference. Gives a third pose that's *not* a CA-only pose_solver output. Then run the `full` contig on it to see whether a natively-derived pose closes the crystal-vs-posed gap (29.7% → 10.2% with v2 → ?% with homolog).

Skip more pose_solver variants until the objective is improved.

### 4.3 Stop running

- **`crystal_auto3A_contig`**: let the queued run finish if already submitted, but don't submit more 3-Å auto-contigs. Confirmed unviable at 3.4% on posed. It will not be informative at crystal.
- **More pose_solver variants** at the current CA-only objective. Pose v2 already matches the crystal on gap-index and clash counts. The remaining BSA gap is structural and won't shrink without heavy-atom-aware constraints.

### 4.4 Analysis infrastructure

Once 5–8 (pose, region) cells have run, build a summary table joining:
- `hits`, `hitrate`
- input-pose metrics: `BSA_total`, `gap_index`, `heavy_contacts_5A`, `clashes_lt_2A`
- contig metrics: `n_design_residues`, `n_design_segments`, `length_range_min`, `length_range_max`, `largest_segment_length`

Even with N=8 first-order trends should be visible (scatter of `n_design_residues vs hitrate` coloured by pose, etc). A two-line Python script can pull the inputs from `solved_pose_results.json` + `interface_metrics.json` + the contig string itself.

### 4.5 Pose-solver code improvements (notes16 §9, still deferred)

- **Heavy-atom-aware pair constraints** (LogSumExp soft-min instead of CA-CA distance). Closes the over-packing loophole structurally.
- **Global sidechain-clash penalty** vs only on pair residues. Pair-only is necessary but not sufficient (notes16 §6 v3).
- **Re-solve PikF** at the new defaults, then re-run `posed_full_test_contig` to see if the crystal/posed gap closes.

### 4.6 When orthogonal cascade finishes

- **Refresh §1 with the orthogonal-validated hit counts** (`survivors_with_orthogonal_metrics.csv`). Negsteer tier is the first gate. The orthogonal cascade further filters for AF3 agreement, BSA, Sc, and interface pLDDT. Expect the ranking to be preserved but absolute numbers to drop.
- **Cross-reference designs across `crystal_full` and `posed_full`**. Are the 13 posed_full hits a subset of the 38 crystal_full hits, meaning the same sequences from a worse pose? Or are they different sequences exploring different design space?

---

## 5. Pose-solver Nextflow integration (this session)

Replaced the HADDOCK3 Nextflow module wholesale with a pose-solver module on `phase4-impl`. Four atomic commits:

```
3edbea3 tests/pose_solver: flesh out reference_manifest.txt with real artefacts
608c74b pose_solver: add 5-plot diagnostic suite
9f9e29a tests: move tests/haddock/ to tests/pose_solver/
e1f23d4 pose_solver: replace HADDOCK Nextflow module with pose-solver pipeline
dd6dd35 pose_solver: move scripts to bin/, preserve target chain, emit loss breakdown + restart CSV
```

### What changed

- **`bin/`**: deleted 9 `haddock*.py` scripts. Added `pose_solver.py`, `pose_interface_metrics.py`, `derive_pose_contigs.py` (moved from `tests/haddock/pose_solver/`), `pose_solver_prepare.py` (chain-normalisation step), and `pose_solver_plots.py` (5-plot diagnostic suite).
- **`modules/`**: deleted `modules/haddock.nf`. Added `modules/pose_solver.nf` with 4 new processes (`POSE_SOLVER_PREPARE`, `POSE_SOLVE`, `POSE_INTERFACE_METRICS`, `POSE_SOLVER_PLOTS`) plus the migrated `BUILD_CONTIGS`.
- **`main.nf`**: replaced 9 `haddock_*` params with 12 `pose_solver_*` params (defaults per notes16, where `max_pair_distance=6.0` matches natural CA-CA range). Branch A rewired through the pose-solver chain.
- **`bin/validate_params.py`**: dropped 8 haddock specs + 2 cross-param rules. Added 12 pose_solver specs + 2 new cross-param rules (Branch A requires `pose_solver_pairs` non-empty, and `min < max` pair distance).
- **`nextflow.config`**: dropped 6 HADDOCK withName blocks (including the 64-CPU exclusive-node sizing on `HADDOCK3_DOCK`). pose_solver runs on 8 CPUs in 2h. Manifest description updated.
- **`tests/`**: `tests/haddock/` → `tests/pose_solver/`. This preserves all campaign result subdirs from notes16, namely `Pikp-1_HMA_PikF_v2/`, `PikF_2pair_min4_max6_scclash/`, `interface_metrics_out/`, etc. New `test_pose_solver.nf` + slurm launcher. New characterization tests in `test_validate_params.py`. Dispatcher (`run_tests.sh`) updated.
- **Chain convention fix**: `pose_solver.py` now preserves the `--target-chain` ID end-to-end (drop hardcoded `tc_out='B'`). New `POSE_SOLVER_PREPARE` step normalises both input PDBs to `params.receptor_chain` / `params.effector_chain`. Closes a latent bug in the old HADDOCK pipeline. The docked complex was always written as A/B regardless of the user's claimed chain IDs. `build_contigs.py` then silently fell back to the first two alphabetical chains.
- **Loss-breakdown emission**: `pose_solver.py` now writes per-term loss components (`validity`, `dist_upper`, `dist_lower`, `excl`, `interp`, `pair_sc_clash`) to `solved_pose_results.json`. It prints them unprompted at end of every run. Per-restart loss history is dumped to `solved_pose_restart_losses.csv` for the convergence plot.

### Five new diagnostic plots (consumed by `POSE_SOLVER_PLOTS`)

1. `pose_loss_breakdown.png`: horizontal bars per penalty term with weight constants annotated.
2. `pose_pair_distances.png`: whisker rows showing allowed window + achieved distance per pair + β-strand reference.
3. `pose_interface_dashboard.png`, with 4 panels: BSA bars (with native ~2700 Å² reference), clash-coloured heavy-atom distance histogram, gap-index vs native ~3.41 Å reference, interface residue composition.
4. `pose_restart_loss_curve.png`: per-restart scatter + best-so-far step + "converged at restart N" annotation.
5. `pose_contig_comparison.png`: three contig options (gap-0, gap-1, gap-2) as residue bars over the binder.

### Verification

- **HPC characterization tests**: 272 passed, 26 skipped (all pre-existing absences in negsteer/orthogonal fixtures), 429 deselected in 34.7s. `slurm_20208055.out`.
- **Pose-solver Nextflow smoke test**: `Success: true`, duration 2m 33s, all 7 processes ✓ on Pikp-1_HMA + avr-pia. `slurm_20208051.out`. Final loss = 0.0 (all 3 pair constraints satisfied), `BUILD_CONTIGS` correctly derived `B1-68` from the preserved chain B (chain-preservation working as designed).
- **Fixture curated**: `tests/pose_solver/example_output_files/` populated with 19 artefacts (4 Nextflow run-artefacts + 5 plots + 2 chain-normalised inputs + 5 solved-pose outputs + 3 BUILD_CONTIGS outputs). Manifest in `tests/pose_solver/reference_manifest.txt`.

### One legacy reference retained

`nextflow.config:31` still points at `/hpc-home/.../HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img`, the literal filename of the Singularity image on HPC. The container itself provides RFDiffusion + MPNN + scipy + matplotlib, and the haddock3 CLI inside it is unused. Rename the `.img` on HPC and update the path when convenient. The same path appears in several `tests/*/run_test_*_plots.slurm.sh` SLURM scripts, which will need the same one-line update after the rename.

---

## 6. File index

New / modified this session (high-signal only, see commit log for full):

| Path | Status |
|---|---|
| `bin/pose_solver.py` | moved + modified (chain preservation, loss breakdown, restart CSV) |
| `bin/pose_solver_prepare.py` | new (chain normalisation) |
| `bin/pose_interface_metrics.py` | moved (renamed from `interface_metrics.py`) |
| `bin/pose_solver_plots.py` | new (5-plot suite) |
| `bin/derive_pose_contigs.py` | moved (renamed from `derive_contigs.py`) |
| `bin/haddock*.py` | **deleted** (9 files) |
| `modules/pose_solver.nf` | new (4 processes + BUILD_CONTIGS) |
| `modules/haddock.nf` | **deleted** |
| `main.nf` | rewritten Branch A, haddock params → pose_solver params |
| `bin/validate_params.py` | rewritten haddock specs → pose_solver specs |
| `nextflow.config` | dropped haddock withName blocks, added pose_solver blocks |
| `params_example.yml` | rewrote HADDOCK section as pose_solver |
| `tests/full_test_run/params_full_test.yml` | dropped HADDOCK params block |
| `tests/pose_solver/` | new home for moved tests/haddock/pose_solver/ campaign results |
| `tests/pose_solver/test_pose_solver.nf` | new (smoke test) |
| `tests/pose_solver/example_output_files/` | new (19-file curated fixture) |
| `tests/characterization/test_validate_params.py` | swapped haddock tests for pose_solver tests |
| `tests/characterization/test_haddock_*.py` | **deleted** (2 files) |
| `tests/run_tests.sh` | haddock → pose_solver in dispatcher |
| `scripts/sync_from_hpc.sh` | haddock → pose_solver in module list |
| `notes/pipeline_notes/pipeline_notes17.md` | new (this file) |
