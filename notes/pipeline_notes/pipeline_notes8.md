# pipeline_notes8.md

Session handover — succeeds todo_list6 / pipeline_notes7. Documents what
landed in the v7 planning session, what's staged for deployment, and
what the next chat needs to pick up. Written 2026-04-24.

---

## Session outcome in one paragraph

Closed out the v6 P0 block by staging code for Tasks 6 / 29 / 31 / 34
/ 38; they flip to DONE once the test cascade runs green. Produced
`todo_list7.docx` restructuring priorities around the user's 1–15
execution plan: Task 32 promoted to P0 (blocks Task 7 throughput),
Tasks 43–49 + 52–53 added (plots, multi-cycle, codebase/threshold
audits, target selection), Tasks 15 + 19 demoted, Tasks 17 + 39
re-scoped as standalone-by-hand analyses, Task 40 discarded. Revised
Task 6 from 2N per-design controls (wrong semantics) to 2 total
controls against the input structure — matches the user's "if
RFDiffusion did a terrible job, what metrics would it get" framing.
Updated test NF files with an upstream-first fallback-to-data/ chaining
pattern so the tests form a natural cascade. User submitted the
`test_rfdiffusion` job (8 designs × 50 iterations) at the end of
the session; that's the first step of the v7 execution plan.

---

## Reference documents

- **`todo_list7.docx`** — authoritative task list going forward.
- **`todo_list6.docx`** — predecessor; still useful for anything v7
  doesn't explicitly revise.
- **`pipeline_notes7.md`** — predecessor; kernel A/B closure, 86 s per
  Boltz prediction, 62 s fixed-overhead decomposition.

---

## What flipped to DONE vs. will-close-on-tests

v6 P0 items now at "will-close-when-tests-pass" state:

| Task | Status | Flips to DONE when |
|------|--------|---------------------|
| 6 | REVISED | `test_negative_steering` → full pipeline run shows 2 control rows (scrambled + polyA) in `cross_sequence_summary.csv` |
| 29 | STAGED | `test_orthogonal_metrics` emits the 8 interface-metric columns |
| 31 | STAGED | `test_orthogonal_metrics` emits `survivors_with_orthogonal_metrics.csv` |
| 34 | STAGED | full `nextflow run main.nf` on small input produces the orthogonal-metrics CSV end-to-end |
| 38 | STAGED | `test_orthogonal_metrics` shows `weighted_jaccard` column populated |

v6 items already DONE (no work here): 30 (kernels rebuilt), 35 (kernel
perf audit — negative-but-decisive), 2 and every other task marked DONE
in v6.

---

## Code staged for deployment

All files live at `/mnt/user-data/outputs/` in the previous chat turn's
output. Directory layout (drop into the live repo preserving paths):

```
bin/
  build_control_sequences.py          NEW (scrambled / polyA synth, ~340 LOC, 8 unit tests)
  compute_interface_metrics.py        MODIFIED (P0-29 + P0-38 weighted_jaccard, 766 LOC)
  cross_sequence_summary.py           MODIFIED (row_type sidecar pickup, +43 LOC)
  derive_input_design_region.py       NEW (contigs → input design region, ~320 LOC, 7 unit tests)

modules/
  negsteer_controls.nf                NEW (DERIVE_INPUT_INDICES + NEGSTEER_CONTROLS)
  negsteer_interface_metrics.nf       MODIFIED (forwards --weighted-jaccard-pair-cutoff)
  negsteer_manifest.nf                NEW (EXTRACT_SURVIVOR_MANIFEST, shared between prod and test)

tests/
  rosetta_filtering/test_rosetta_filtering.nf   MODIFIED (upstream chaining)
  proteinmpnn/test_proteinmpnn.nf               MODIFIED (upstream chaining)
  negative_steering/test_negative_steering.nf   MODIFIED (upstream chaining, controls OFF,
                                                          negsteer_no_kernels flipped true)
  orthogonal_metrics/test_orthogonal_metrics.nf MODIFIED (upstream chaining, path fix,
                                                          controls OFF)

main.nf                               MODIFIED (Steps 2+3+4 of session + Task 6 revision:
                                                 2-total input-structure controls, onComplete
                                                 miscalibration scanner; 987 lines, 132/132
                                                 braces balanced)
params_example.yml                    MODIFIED (P0-29, P0-38, Task 6 params documented)
todo_list7.docx                       NEW
```

### Notes on specific files

- **`main.nf` Task 6 implementation**: `DERIVE_INPUT_INDICES` runs once
  per cohort after whichever preprocessing branch produced
  `rfdiff_pdb_ch` and `contigs_ch`. A 2-element
  `Channel.of("scrambled", "polyA")` is combined with the input-structure
  bundle to give exactly 2 controls per cohort. Gated on
  `params.run_negative_controls = true`. The onComplete block scans the
  final `cross_sequence_summary.csv` for control rows breaching
  `controls_warning_ipsae_max=0.5` or `controls_warning_ra_eff_min=5.0`
  and emits a ⚠ RANKER MISCALIBRATION banner. Warns only, does not fail.

- **`derive_input_design_region.py`**: added, NOT replacing the existing
  `derive_design_region.py` or `derive_true_interface.py`. The per-
  design path (reading `rfdiffusion_metrics.json`) still uses those
  scripts. The new script parses the contigs string + input PDB to
  compute: (a) the design region as the receptor positions between
  contig anchor segments (positions RFDiffusion was told to replace),
  (b) the true interface via heavy-atom contacts on the input complex.
  Emits files in the exact formats `derive_design_region.py` /
  `derive_true_interface.py` write, so downstream code consumes them
  identically.

- **`negsteer_controls.nf` NEGSTEER_CONTROLS process input tuple**:
  `(control_name, control_type, input_pdb, design_region, true_interface)`
  — 5 elements, NOT the v6 7-element per-design shape. Effector is
  extracted inline from the input_pdb via `get_chain_sequence`. The
  control uses `input_pdb` as `--ground-truth` since the controls run
  against the input structure.

- **`compute_interface_metrics.py` weighted_jaccard implementation**:
  per-pair Cβ–Cβ distance machinery (Cα fallback for glycine / Cα-only
  RFDiffusion designs). Hard-coded Gaussian parameters (μ=4.0,
  2σ²=2.25) — deliberately NOT CLI-tunable to keep cohort
  comparability. Only the pair inclusion cutoff (default 8 Å,
  numerically inert since weight at 8 Å is ≈8e-4) is a CLI flag.
  Returns NaN for empty+empty (matches existing jaccard() convention),
  0.0 for empty+nonempty, value otherwise.

- **Test NF chaining pattern** (four files, consistent idiom):
  ```groovy
  params.upstream_<stage>_outdir = "${projectDir}/../<prev_test>/results"
  def _upstream_glob = "${params.upstream_<stage>_outdir}/<path>/*.pdb"
  def _fallback_glob = "${projectDir}/data/*.pdb"
  params.<input> = files(_upstream_glob).size() > 0 \
      ? _upstream_glob \
      : _fallback_glob
  ```
  Uses `files()` (returns a list) not `file()`. Nextflow's CLI param
  priority means `--<input> /direct/path` on the CLI still wins —
  confirmed by the docs. Each test logs which source it selected at
  workflow start.

- **`test_negative_steering.nf` extras**:
  - `negsteer_no_kernels = true` (was `false`) — matches production per
    P0-35's finding that cuEquivariance kernels are ~5% slower on the
    ~155-residue complexes.
  - `run_negative_controls = false` — the negsteer test times the inner
    loop, not the loop + 2 extra control GPU jobs.

- **`test_orthogonal_metrics.nf` path fix**: the v6 version had a stale
  default `../negative_steering/receptor_resurfacing_results/negative_steering`.
  Now derives cleanly from `upstream_negsteer_outdir` → `results/negative_steering/`.

---

## Task 6 revision — what changed from the earlier implementation

The earlier (v6) implementation generated 2N controls, one per design,
with receptor source = the design PDB. This session replaced it with
2 total controls per cohort, receptor source = the input structure
(rfdiffusion_input.pdb).

**Why**: the wet-lab-decision-relevant question ("is the ranker fooled
by sequences the predictor has every reason to reject?") is answered
by 2 total controls across a cohort. Per-design pathology is a
narrower question that Task 8 AUROC answers better. User framing:
"if RFDiffusion did a terrible job, what metrics would those
sequences get?"

**Design region for controls**: on the input receptor, the "design
region" is the **gap between anchor segments** in the contigs string
— the positions that RFDiffusion would rip out and replace with de
novo chemistry. For contigs `A1-50 5-15 A60-100 B` on a 100-residue
receptor, the design region is positions 51–59 on the input receptor.

**polyA semantics**: design-region positions → A, EXCEPT native
glycines are preserved as G (avoids forcing alanine into
glycine-only Φ/ψ wells that glycine specifically enables).

---

## Priority changes in v7 vs. v6

### Promoted

- **Task 32 (RFDiffusion parallelisation)** P1 → P0. Blocks Task 7
  throughput. Begins with a `trace.txt` diagnostic step confirming
  the serial-execution hypothesis; refactor only if confirmed.

### Demoted

- **Task 15 (delta-to-cold-start columns)** P2 → P3. User rationale:
  "wouldn't be used for ranking. Diagnostic only. Metric overkill
  for now." Revisit if Task 8 AUROC suggests value.
- **Task 19 (candidate-pool distance bias)** P1 → P3. User rationale:
  "something to be comparing rather than front-loading." Run as a
  comparison experiment after Task 7 on one of the three targets.

### Re-scoped as standalone (NOT pipeline integration)

- **Task 17 (submit_final_validation)** — standalone wrapper script,
  run by hand on survivors or experimental structures. Not every
  analysis needs it.
- **Task 39 (DockQ vs experimental structure)** — standalone analysis
  script on survivors vs. ground-truth binding structure. Feeds Task 8.

### New tasks

- **Task 43** P1 — Negative-steering plots (per-sequence + per-cohort).
- **Task 44** P1 — Orthogonal-metrics plots (distributions + scatter pairs).
- **Task 45** P1 — Multi-cycle negative steering (full implementation push).
- **Task 46** P1 — Codebase audit + cleanup (structural pass).
- **Task 47** P1 — Threshold audit (policy pass).
- **Task 48** P1 — Regression-test gate after 46/47.
- **Task 49** P1 — 3-target selection + contigs.
- **Task 52** P3 — Negsteer 10-candidate blocking (early termination).
- **Task 53** P3 — Shell-radius contig mode.

### Discarded

- **Task 40** — Full-length receptor autoactivity filter. User
  rationale: "not convinced the whole NLR will be well predicted."

### Still parked, unchanged

- **Task 7** P1 — cohort run (gated on 47 + 48 + 49).
- **Task 8** P2 — cross-target per-metric AUROC (after Task 7).
- **Task 18** P2 — noise-halo (gated on Task 8).
- **Tasks 20 / 21 / 22** P2 — post-wet-lab validation (blocked on hits).
- **Task 42** P2 — `boltz predict` batching (absorbs "runtime
  optimisation" umbrella).

---

## Execution order (the 1–15 plan from the user)

This is the shape of the next few weeks / months. Numbers are steps,
not task IDs — maps to task IDs in the right column.

| # | Step | Task |
|---|------|------|
| 1 | Individual section tests | closes 6 / 29 / 31 / 34 / 38 |
| 2 | Full pipeline run (small) | 34 |
| 3 | Plots refresh | 43 + 44 |
| 4 | RFDiffusion parallelisation | 32 |
| 5 | Multi-cycle negative steering | 45 |
| 6 | Codebase audit / cleanup | 46 |
| 7 | Threshold audit | 47 |
| 8 | Regression-test gate | 48 |
| 9 | 3-target selection + contigs | 49 |
| 10 | Cohort run on 3 targets | 7 |
| 11 | submit_final_validation standalone | 17 (re-scoped) |
| 12 | DockQ vs experimental standalone | 39 (re-scoped) |
| 13 | Cross-target AUROC | 8 |
| 14 | Noise-halo thresholds | 18 |
| 15 | Delta-to-cold-start columns (maybe) | 15 |

No formal "codebase freeze" gate — user opted for implicit discipline.
The implicit rule: once Task 7 starts, no pipeline-behaviour changes
until all three cohorts complete.

---

## Open at session end

### In flight

- **`test_rfdiffusion` submitted** — user's change: 8 designs, 50
  iterations (up from default 4 and 25). Expected 30–40 minutes if the
  jic-gpu queue is light. This is step 1 of the execution plan.

### Session one-liners not yet done

None — everything discussed in the session that was agreed has been
staged.

### Pending in v7 but not started

The next chat should pick up with whatever step 1 produces:

1. **If `test_rfdiffusion` passes**: proceed to
   `test_rosetta_filtering` (it should automatically pick up the
   8 designs from `../rfdiffusion/results/rfdiffusion/split/` via
   the new chaining pattern).
2. **If anything fails**: diagnose against the staged code; likely
   suspects would be the contigs-parser in
   `derive_input_design_region.py` on an unusual contig form, or
   the `files()` glob pattern in the chained test NFs misbehaving
   on an empty directory.

### Design decisions flagged for review once cohort data exists

- `controls_warning_ipsae_max = 0.5` and `controls_warning_ra_eff_min = 5.0`
  were set on gut feel. Task 47 (threshold audit) will confirm
  against the design_24 variance data.
- Weighted-Jaccard formulation: implemented as per-pair Cβ–Cβ
  (strictly more discriminating) rather than residue-set weighted.
  This was a design call made during Step 3 — both the existing
  `true_jaccard` and the new `weighted_jaccard` flow through cross-
  sequence summary, so Task 8 AUROC will decide which wins.

---

## Things the next chat will likely need to know

1. **Brace balance conventions**: `main.nf` at 132 / 132,
   `negsteer_controls.nf` at 61 / 61. Every edit that changes control
   flow should re-check with the same regex-strip-then-count approach.

2. **`test_orthogonal_metrics` smoke-test path**: if a full run is
   too expensive, the test runs against cached negsteer output. Path:
   `tests/negative_steering/results/negative_steering/`. The chaining
   default points there.

3. **HADDOCK is shelved**: the user confirmed HADDOCK is not part of
   the current pipeline. It stays in the codebase but is NOT part of
   the test cascade. `test_haddock` was deliberately left un-chained.

4. **cuEquivariance kernels are OFF in production**: per P0-35's
   closure. `negsteer_no_kernels = true` in both production params
   and the negsteer test defaults. The rebuilt container
   `boltz2_negsteer.img` at `/hpc-home/jowillia/singularity/...`
   still has the kernel path available — we just don't take it.

5. **`contact_cutoff` default of 8 Å is the right ballpark** for
   RFDiffusion filtering. Cα–Cα because designs are polyvaline at
   this stage. Do NOT tune this for test variety — it's a calibrated
   threshold the downstream pipeline assumes. Task 47 revisits it
   with real data. Knobs safe to vary: `num_designs`, `rfdiff_iterations`,
   `hotspot`.

6. **MPNN chaining default** uses `qc_fastas/` (always produced),
   not `top_fastas/` (only when `mpnn_top_n > 0`). If you set
   `--mpnn_top_n N` in the MPNN test, switch the chaining default
   or pass `--fastas_dir` directly to the negsteer test.

7. **Task 45 (multi-cycle) has two unresolved decisions** the v7
   spec flagged: (a) convergence criterion (no-new-candidates vs
   metric-stalls vs fixed-count); (b) whether controls get multiple
   cycles. Spec recommends "no-new-candidates with fixed-count
   fallback" and "controls stay single-cycle" — but these are to be
   finalised during implementation.

8. **Nextflow DSL2 `params` priority** behaves as
   expected: CLI `--foo bar` overrides script-level `params.foo =`
   defaults, and closures that read `params.foo` in the script see
   the post-override value. The `files(glob).size() > 0 ? a : b`
   pattern used for test chaining relies on this and it works.

---

## If the session runs out of context in the new chat

The minimum context to carry forward:

- This file (pipeline_notes8.md).
- `todo_list7.docx` — task priorities.
- The 10 files under `/mnt/user-data/outputs/` from the previous
  session (the staged code).
- The user's 1–15 execution plan (§Execution order above).

Everything else in the v7 session can be re-derived from those.
