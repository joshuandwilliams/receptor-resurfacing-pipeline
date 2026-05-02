# 14 — Phase 2 Revision: Per-Module Tests as the Primary Safety Net

A revision to `notes/inventory/11_phase_2_plan.md`. The original plan stands as the record of what we set out to do; this document supersedes the parts that have changed and explains why.

The pivot in one sentence: **the per-module tests, fed from curated fixture inputs, become the primary characterization safety net; the full pipeline run is rare and milestone-only.**

---

## 1. Why this revision

The Phase 2 plan as written assumed a single canonical reference run — the supervisor-demo, 32 designs × 4 MPNN sequences, ~6 hours — would be both the source of reference outputs *and* the verification target during Phase 3+ iteration.

That assumption was wrong about iteration cost. Each Phase 3 verification cycle would have meant a 6-hour HPC run. Across the dozens of Phase 3 commits ahead, the cumulative cost is unworkable.

Two further problems surfaced as we worked through the discovery process:

- **Path coverage was implicit, not designed.** We had to retroactively go searching for which sequences in the supervisor-demo cohort exercised which paths through negative steering (cold-start outcomes, contamination, reversion succeed/fail). The trio we ended up with (`input_control_polyA`, `design_0_seq_0`, `design_13_seq_2`) was a best-effort selection, not a coverage-guaranteed fixture.

- **Per-module tests have been silently chained.** Tests like `test_negative_steering.nf` were written to consume upstream tests' outputs as a convenience during one-shot manual testing. This means each per-module test isn't actually isolated — running a downstream test depends on having recently run the upstream tests. For *characterization* purposes, where reproducibility is everything, this coupling is a defect.

The pivot below addresses all three.

---

## 2. The new strategy

### 2.1 Two-tier verification

| Tier | What | Cost | Cadence |
|---|---|---|---|
| Per-module tests (each in isolation) | Stage-by-stage characterization on fixed input fixtures. Reference outputs from each `tests/<module>/results/`. | 5–30 minutes per module; full sweep ~1–2 hours. | Every Phase 3 commit (or batch of commits). |
| Full pipeline run | End-to-end smoke test that the modules still wire together correctly. | 2–4 hours (with reduced `negsteer_n_designs`). | Only at major milestones (e.g., end of Phase 3, end of Phase 4). |

### 2.2 Each per-module test runs from a fixed fixture

Every per-module test gets a `data/` directory under `tests/<module>/` containing the inputs that test consumes. Inputs are committed to the repo. No upstream chaining for characterization purposes — if a test currently auto-discovers upstream outputs, it must fall back to `data/` whenever those upstream outputs are absent (the existing fallback logic in `test_negative_steering.nf` does this; other tests may need adjustment).

The fixture inputs for each per-module test are **curated** — chosen to exercise specific paths through that module's logic. They are not a random sampling.

### 2.3 The reference set under the new model

`tests/full_test_run/example_output_files/` is replaced by a per-module structure:

```
tests/<module>/example_output_files/      <- new location
```

Each per-module example_output_files/ mirrors the structure that module's test produces. The `tests/full_test_run/example_output_files/` directory remains only for the rare full-pipeline-run characterization (and may eventually be deleted in favour of relying on the per-module sets).

The characterization tests in `tests/characterization/` are restructured accordingly:

- Each `test_<stage>.py` resolves its `reference_root` to `tests/<module>/example_output_files/` rather than `tests/full_test_run/example_output_files/`.
- Each `test_<stage>.py` resolves its `output_root` to `tests/<module>/results/` (the directory the per-module test writes to) rather than a fresh full-pipeline run.
- The parametrize lists for per-sequence tests update to whatever sequence names the new fixtures produce.

### 2.4 The discovery run

Before any of the above can happen, we need to *find* which inputs to put into each per-module test's fixture. This is the discovery run — one expensive pass whose purpose is to surface candidates.

**Params for the discovery run:**

- `rfdiffusion_n_designs = 64` (vs. supervisor-demo 32, for more design diversity)
- `mpnn_seqs_per_design = 2` (vs. supervisor-demo 4)
- All 128 MPNN sequences forward to negative steering (no top-N filtering, or N=128)
- `negsteer_n_designs = 4` (vs. production 20 — the cost-driving parameter)
- All other params match production
- HADDOCK branch off (Branch B only)

Estimated runtime: 2–4 hours, dominated by RFDiffusion generation and 128 × 4 negsteer mutation jobs.

**What the discovery run answers:**

For Stage 2 (RFDiffusion):
- Which design indices produce contigs that pass / fail downstream? (Need at least one of each for `tests/rfdiffusion/data/` to cover the filter path.)

For Stage 3 (Rosetta filtering):
- Which RFDiffusion designs pass / fail Rosetta filtering?

For Stage 4 (MPNN):
- Which designs produce sequences that go on to clean cold-start vs. failed cold-start (i.e., produce viable input for negative steering vs. not)?

For Stage 5 (negative steering):
- For each (MPNN sequence, design index) pair: what path did it take?
  - Cold-start fail (effector placed incorrectly)
  - Cold-start pass, steering produces no contamination
  - Cold-start pass, steering contamination, reversion succeeds (correct placement persists after reverting contaminating mutations)
  - Cold-start pass, steering contamination, reversion fails

For Stage 6 (orthogonal metrics):
- Which sequences make it through to orthogonal metrics? (Need at least one survivor for the orthogonal-metrics test to have meaningful inputs.)

### 2.5 Fixture curation procedure

After the discovery run completes:

1. **Inventory per-module path coverage.** For each pipeline stage, identify which inputs from the discovery run exercise which paths.
2. **Choose minimum-coverage fixtures.** For each per-module test, select the smallest set of inputs that hits every path the module is meant to exercise:
   - `tests/rfdiffusion/data/` — the input PDB and contigs the discovery run used. Reference outputs cover passing and failing designs.
   - `tests/rosetta_filtering/data/` — RFDiffusion outputs from the discovery run, copied in. Reference outputs cover the filter pass/fail split.
   - `tests/proteinmpnn/data/` — selected RFDiffusion designs from the discovery run.
   - `tests/negative_steering/data/` — selected (MPNN sequence, design PDB, metrics JSON) bundles known to exercise each negsteer path. Aim for 4–6 sequences total covering: cold-start fail, cold-start pass + clean steering, cold-start pass + contaminated + reversion success, cold-start pass + contaminated + reversion fail. The two controls are generated by the test itself from `params.input_pdb`.
   - `tests/orthogonal_metrics/data/` — survivors from the discovery run's negsteer output.
3. **Commit fixtures to repo.** These are the canonical test inputs going forward. They are static; they only change when the test fixture deliberately needs updating.
4. **Run each per-module test on HPC** with its committed fixture. Capture the output as the per-module reference set.
5. **Subtractive rebuild** of each `tests/<module>/example_output_files/` from its per-module test's `results/`, using the same approach as the original supervisor-demo rebuild (mirror structure, exclude PDBs/NPZs/MSA/predictions/ subtrees). Commit.

### 2.6 The 184 existing characterization tests

The tests written in Prompt 3 are tied to:

- The `tests/full_test_run/example_output_files/` reference path
- The trio `input_control_polyA, design_0_seq_0, design_13_seq_2`

Both will change under the new strategy. The tests themselves are still useful — the *comparators* are correct, the *strategies* are correct — but the resolution paths and parametrize lists need updating.

This is a substantive but mechanical rewrite. Proposed at the appropriate point in §4 below.

---

## 3. What stays from the original Phase 2 plan

Most of the framework. To be explicit:

- `pyproject.toml`, marker registry, `--strict-markers` — unchanged.
- The four comparator helpers (`csv_compare`, `json_compare`, `png_compare`, `path_normalize`) plus `text_compare` and `compare_csv_exact_modulo_paths` — unchanged.
- The 67 `local_unit` tests — unchanged.
- `ComparisonResult` dataclass and all fixture infrastructure — unchanged.
- The two-wave structure (Wave 1 = output-file pinning; Wave 2 = Python-internal state) — unchanged.
- The manual reference-update procedure (plan §8.1) — applies to per-module reference sets too, just per-module instead of cohort-wide.
- The Verification Queue entries from §7.2 of the original plan — unchanged.

---

## 4. Revised execution sequence

Replaces plan §6.

### 4.1 Pre-conditions (already met)

- ✅ Cache-busting audit complete (commit `phase-2.4-cache-busting-complete`).
- ✅ Comparators + unit tests written.
- ✅ Pre-HPC-roundtrip audit complete (`audit_pre_hpc_roundtrip.md`); JAVA_HOME fixed.
- ✅ Repo synced to HPC.

### 4.2 New work to do

1. **Decouple per-module tests from upstream chaining.** Where a per-module test currently reads from `tests/<previous_module>/results/`, change it to read from `tests/<this_module>/data/` as the primary path, falling back to upstream-chained outputs only as a convenience for ad-hoc local testing. For characterization purposes, the `data/` path is canonical.

2. **Run the discovery run on HPC.** Params per §2.4 above. Capture the full `results/` tree.

3. **Inventory the discovery run's outputs** to identify path coverage per stage.

4. **Curate per-module fixtures** per §2.5. Commit to the repo.

5. **Run each per-module test on HPC** with its committed fixture. Verify each test produces a consistent output tree.

6. **Subtractive rebuild** of each `tests/<module>/example_output_files/`. Commit.

7. **Restructure the existing 184 characterization tests** to use per-module reference paths and updated parametrize lists. The framework code (`conftest.py`, fixtures) updates so `reference_root` and `output_root` are *per-stage* fixtures, parameterized over stages. Each test file then resolves its own paths via those fixtures.

8. **Round-trip verify.** Run each per-module test on HPC; run `pytest -m hpc` against the per-module outputs; expect green.

9. **Safety-net validation commit.** Delete `bin/sequence_registry.py` (and the `main` empty file if it still exists). Re-run all per-module tests + characterization. Expect zero diffs.

### 4.3 Tag sequence under the revision

The existing tags (`phase-2.1-plan-complete` through `phase-2.4-cache-busting-complete` plus `phase-2.5-reference-rebuild-complete`) stand. Going forward:

- `phase-2.6-revision-plan-complete` — this document.
- `phase-2.7-discovery-run-complete` — after the discovery run executes and is curated.
- `phase-2.8-fixtures-and-references-complete` — after per-module fixtures and reference sets are in place.
- `phase-2.9-tests-restructured-complete` — after the 184 characterization tests are updated to the per-module structure.
- `phase-2-complete` — after the safety-net validation deletion commit produces zero diffs.

The previously-anticipated `phase-2.6-hpc-tier-complete` (tests green on full pipeline run) is superseded — that scenario now happens at `phase-2.9` with per-module verification.

---

## 5. What to retire from the existing reference set

`tests/full_test_run/example_output_files/` will be deleted once the per-module reference sets are in place. Retiring it earlier would leave the existing 184 characterization tests broken (they currently resolve against it). The order of operations is:

1. Build per-module references (§4.2 step 6).
2. Restructure tests to point at per-module references (§4.2 step 7).
3. Verify tests pass (§4.2 step 8).
4. Then delete `tests/full_test_run/example_output_files/`.

Until step 4, the existing reference set stays — it's just no longer the canonical iteration target.

---

## 6. Open questions

These don't block writing this document, but they will need answering during execution.

1. **HADDOCK branch.** Branch A is out of Wave 1 scope. Should `tests/haddock/` get a per-module reference set under the new strategy at all, or stay deferred until Branch A is redesigned? *Suggested answer:* defer.

2. **Disk-space cost of per-module reference sets.** Each `tests/<module>/example_output_files/` mirrors that module's run outputs. Five or six modules × moderate output size could add up. Worth measuring after the first per-module rebuild to see whether further trimming (similar to the supervisor-demo subtractive rebuild) is needed.

3. **Stochasticity in module outputs.** Even with fixed inputs, some Boltz-2 prediction outputs may vary slightly run-to-run if seeds aren't fully pinned. The first per-module round-trip will reveal which outputs need the comparator strategy demoting (e.g., JSON-DEEP → JSON-STRUCT).

4. **The discovery run's params file.** Needs writing as a new `params_discovery.yml`. Suggest creating it as part of §4.2 step 2 — it's a small YAML file but worth keeping in the repo for reproducibility.

---

## 7. What this document supersedes in the original Phase 2 plan

| Plan §  | Status |
|---|---|
| §1 (framing, what these tests are/aren't) | Unchanged. |
| §2 (file-by-file scope) | **Restructured** — file-by-file scope now applies per-module rather than cohort-wide. |
| §3 (test framework architecture) | Mostly unchanged; conftest fixture resolution updates per §4.2 step 7. |
| §4 (comparator specifications) | Unchanged. |
| §5 (Wave 2) | Unchanged in intent; deferred timing per the new sequence. |
| §6 (entry conditions, execution sequence) | **Replaced** by §4 of this document. |
| §7 (risks, Verification Queue) | Unchanged. |
| §8 (deliverables) | Adjusted — Wave 1 deliverables are now per-module reference sets, not a unified one. |
| §9 (decisions made) | Unchanged historical record. |

The original plan stays in the repo as a record of the original design and the reasoning behind each choice. This document is the active spec.
