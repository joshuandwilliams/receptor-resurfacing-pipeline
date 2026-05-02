# pipeline_notes9.md

Session handover — succeeds todo_list8 / pipeline_notes8. Documents
what landed in the long debugging session that ran the v7 test
cascade end-to-end for the first time, performed a major architectural
refactor of the metric-computation path, and submitted both an
overnight test_negative_steering rerun (with controls) and a parallel
production main.nf run. Written 2026-04-25 evening.

---

## Session outcome in one paragraph

Ran the v7 test cascade (rfdiffusion → rosetta → mpnn → negsteer)
end-to-end for the first time, hitting and fixing ~13 bugs along the
way (cascading chaining defaults, AF3 mmCIF format issues, missing
function `find_contact_residues_heavy`, malformed input PDB with
chains in separate MODEL blocks, chain-convention confusion between
input-PDB and RFDiffusion-output conventions). Performed a major
architectural refactor: four metric families (15 Å iPSAE,
interface_plddt, intact_core, weighted_jaccard) moved from post-hoc
orthogonal-metrics into per-prediction computation inside negative-
steering. This avoids the sidecar-file path-resolution issues that
had been silently leaving columns blank, and surfaces all
per-prediction metrics at the same aggregation level
(`representative_*_median`). Submitted overnight: a
test_negative_steering rerun with controls ON + the metric refactor
patches deployed, plus a parallel production main.nf run on PikP1 +
AvrPikF. User asleep at write time; both runs expected complete by
morning.

---

## Reference documents

- **`todo_list8.docx`** — authoritative task list going forward.
- **`todo_list7.docx`** — predecessor (created this session — was
  described in pipeline_notes8 but never actually saved on disk,
  reconstructed from the spec).
- **`pipeline_notes8.md`** — predecessor; documented the v7 staging.
- **All 6 patched files at `/mnt/user-data/outputs/`** from this
  session — these are the metric-refactor patches plus the
  `find_contact_residues_heavy` fix. Should already be live in `bin/`
  by morning.

---

## What flipped status this session

| Task | v7 → v8 | Trigger |
|------|---------|---------|
| 6 | STAGED → IN PROGRESS | Controls wired into test rerun, running tonight |
| 29 | STAGED → REFACTORED | Metrics moved per-prediction |
| 31 | STAGED → PARTIALLY REFACTORED | DockQ stays post-hoc; rest moved per-prediction |
| 34 | STAGED → IN PROGRESS | Production main.nf submitted in parallel |
| 38 | STAGED → REFACTORED | weighted_jaccard now in compute_metrics.py |
| 32 | OPEN → ✖ CLOSED | Per-design ~35s makes parallelisation overhead not worth it |

---

## Bugs fixed this session

Listed in roughly the order they were hit:

1. **Chaining default outdir mismatch** in 4 test files (rosetta,
   mpnn, negsteer, orthogonal_metrics) — pointing at the wrong
   upstream directory.

2. **`rfdiffusion.nf` publishDir double-`split/`** path bug.

3. **`rosetta_filtering.nf` missing publishDir for `passing/*.pdb`**
   — added.

4. **`test_proteinmpnn.nf` input_pdb hardcoded** to a stale path —
   repointed to `${projectDir}/../rfdiffusion/data/`.

5. **`test_negative_steering.nf` chaining preference**: prefer
   `top_fastas/` when present, fall back to `qc_fastas/` otherwise.

6. **AF3 model-dir path mismatch** (`$HOME/af3_models` in params,
   `/hpc-home/jowillia/singularity/AlphaFold3` in container) — fixed
   in nextflow.config.

7. **AF3 input format `mmcifPath` → `mmcif` inline**: rewrote
   `negsteer_af3_nomsa.nf` to inline the mmCIF text in the JSON via a
   Python heredoc rather than the file-path form (which AF3 expected
   in a different schema version).

8. **gemmi minimal mmCIF missing `_entity_poly_seq` block**:
   `boltz2_negative_steering.py:extract_effector_template_cif` now
   uses `gemmi.MmcifOutputGroups(True)` to emit a complete mmCIF
   document.

9. **Chain convention architectural fix**: introduced
   `params.rfdiff_output_receptor_chain="A"` and
   `params.rfdiff_output_effector_chain="B"` alongside the existing
   `params.effector_chain` (which describes the input PDB
   convention). RFDiffusion's `write_split_pdb` hardcodes A/B
   regardless of input — that's the convention used by everything
   downstream of `split/`. Updated `main.nf` (NEGSTEER_RUN_ONE call)
   and `test_negative_steering.nf`.

10. **`run_test_negative_steering.slurm.sh` had been overwritten** —
    user restored from backup.

11. **Negative controls weren't wired into `test_negative_steering.nf`**
    — added the controls workflow block.

12. **Missing function `find_contact_residues_heavy`** in
    `compute_metrics.py`. Pre-existing latent bug. Function was
    imported by `derive_input_design_region.py` but had never been
    written. Only fired when controls actually ran for the first
    time (controls are the only path that exercises
    DERIVE_INPUT_INDICES). Wrote ~60 lines, unit-tested on synthetic
    PDB.

13. **Input PDB had non-standard MODEL/ENDMDL structure**:
    `af3_pikp1_native_avrpikf_complex.pdb` had chains A and C in
    separate MODEL blocks. AF3 emits PDBs this way by default.
    `extract_sequences_gemmi` does `break  # first model only`, so
    only saw chain A; my plain-text `_parse_heavy_atoms_by_chain`
    ignored MODEL boundaries, so saw both. User fixed in place
    (stripped MODEL/ENDMDL/END, kept ATOM and TER, re-added END).
    Backup saved as `.original_2model`.

---

## Major architectural refactor — metric relocation

This is the structural change worth highlighting. Four metric
families moved from post-hoc orthogonal-metrics into per-prediction
computation:

**Moved to per-prediction (inside `compute_metrics.py`)**:
- `ipsae_ab_15`, `ipsae_ba_15`, `ipsae_min_15` — 15 Å iPSAE
- `interface_plddt` — mean Cα pLDDT over interface residues
- `intact_core` — binary Kabsch-aligned receptor-core RMSD < threshold
- `weighted_jaccard` — Gaussian-weighted Cβ-Cβ contact overlap vs native

**Kept in orthogonal-metrics (genuinely post-hoc)**:
- DockQ family (`irmsd`, `fnat`, `dockq`) — needs DockQ binary
- `bsa`, `interface_hbonds` — needs FreeSASA, MDAnalysis
- Rosetta `sc`, `rosetta_ddg` — needs Rosetta InterfaceAnalyzer

**Why the refactor**: the post-hoc path was failing to find sidecar
`pae_*.npz` files for the non-DockQ metrics, leaving columns blank in
`cross_sequence_summary_with_interface_metrics.csv`. Per-prediction
computation runs at the moment the PAE matrix is already loaded — the
metric is essentially free. The new metrics flow through the existing
aggregator → passing_summary → cross_sequence_summary chain as
`representative_*_median` columns.

**Files patched** (all in `/mnt/user-data/outputs/`):
- `compute_metrics.py` — added the 4 metric families + 6 new fields +
  `find_contact_residues_heavy`
- `boltz2_iterate_steering.py` — added 6 new keys to `METRIC_KEYS_RAW`
  + `_AGG_CONTINUOUS_METRICS` (5) + `_AGG_BINARY_METRICS` (intact_core)
  + passes `--native-pdb` to compute_metrics.py subprocess
- `compute_interface_metrics.py` — trimmed to DockQ family only;
  failures column renamed `dockq_failures`; old CLI flags kept as
  deprecated for backward compat
- `run_biophysical_metrics.py` — removed `interface_plddt` from output
- `extract_passing.py` — added 6 median + 1 majority columns
- `merge_orthogonal_metrics.py` — sources `interface_plddt` from
  `representative_interface_plddt_median` instead of
  `biophysical_summary.csv`

**Cache impact**: every cached `NEGSTEER_RUN_ONE` job invalidated
because the script signature changed. Full re-run required (~12
GPU-hours for 8 sequences + 2 controls).

---

## Currently in flight (overnight 25 → 26 April)

1. **test_negative_steering rerun** — submitted as
   `sbatch tests/negative_steering/run_test_negative_steering.slurm.sh`.
   Has controls ON, all metric-refactor patches deployed. ~12
   GPU-hours expected. Will produce: 8 normal-sequence rows + 2
   control rows in cross_sequence_summary.csv. onComplete banner
   should warn if any control row scores ipSAE > 0.5 or ra_eff < 5 Å.

2. **Production main.nf small-scale run** — submitted as
   `sbatch /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/run_pipeline.slurm.sh ./params.yml`.
   Same target PDB. params.yml composed from test values + overnight-
   sized parameters (8 designs × 16 MPNN seqs → 8 negsteer sequences).
   `-resume` capable for morning recovery.

If both runs complete cleanly, Tasks 6 / 29 / 31 / 34 / 38 all flip
to DONE and the v9 revision can move on to step 3 of the 1–15 plan
(plots refresh — Tasks 43 + 44).

---

## Open at session end

### Action items if both runs succeed

1. Verify `cross_sequence_summary.csv` has the 6 new
   `representative_*_median` columns populated for tier-A and
   tier-B rows.
2. Verify the 2 control rows are present and have `row_type` populated
   correctly.
3. Verify `survivors_with_orthogonal_metrics.csv` from the production
   run has DockQ family + BSA + Rosetta columns populated.
4. Move on to step 3 of the 1–15 plan: write/update Task 43 (negsteer
   plots) and Task 44 (orthogonal-metrics plots).

### Action items if either run fails

1. Diagnose against the staged code — the metric refactor changed a
   lot of column names, so most failure modes will be DictReader/
   DictWriter mismatches in the aggregator chain.
2. If DERIVE_INPUT_INDICES failed: check the
   `find_contact_residues_heavy` function in compute_metrics.py is
   what actually shipped (the function was new this session).
3. If a control ran but didn't produce row_type column: check that
   `cross_sequence_summary.py` correctly reads `row_type.txt`
   sidecars from each per-sequence workdir.
4. AF3 fix re-validation: if an orthogonal-metrics run hit AF3, check
   that the `mmcif` inline format works on the cluster's AF3 build
   (the fix worked in isolation; cluster build version not verified
   end-to-end).

### Items the next chat will likely need to know

1. **Patched files at `/mnt/user-data/outputs/`** are the source of
   truth — these should be live in `bin/` by morning. There are 6 of
   them.

2. **Chain conventions**: `params.effector_chain` describes the
   INPUT PDB. `params.rfdiff_output_effector_chain` describes
   RFDiffusion split-output PDBs (always A/B because
   `write_split_pdb` hardcodes them). Most processes consume the
   input convention; only NEGSTEER_RUN_ONE consumes the rfdiff_output
   convention. Don't mix them up — this caused real confusion and
   required the architectural fix in this session.

3. **AF3 PDB output format**: AF3 emits PDBs with chains in separate
   MODEL blocks. gemmi parsers see only the first model; plain-text
   parsers see everything. This is the bug class flagged for Task 46
   audit. Standardising to one convention would prevent a recurrence.

4. **`find_contact_residues_heavy` is positional-index based**, NOT
   seqid-based. It walks receptor residues in PDB order and returns
   `[(positional_index_0b, closest_distance), ...]`. The
   `expected_rec_seq` argument is for sanity-checking that the
   walked-residue count matches the FASTA sequence length.

5. **`derive_input_design_region.py`** is a different script from
   `derive_design_region.py`. The former is for control rows
   (computing the design region by parsing contigs + input PDB).
   The latter is for steered designs (reading
   `rfdiffusion_metrics.json`).

6. **Cache state**: cached NEGSTEER_RUN_ONE jobs from before this
   session are invalidated. If the test rerun is incomplete in the
   morning, `-resume` will only resume from completed-this-session
   jobs, not pre-this-session ones.

7. **Dual-run setup**: the test cascade and the production main.nf
   share the negsteer modules. They run in different output
   directories (`tests/negative_steering/results/` vs `./results/`)
   so their `.nextflow/` caches don't conflict. With
   `max_boltz2_parallel: 16` per instance, the two together can
   queue up to 32 NEGSTEER_RUN_ONE jobs.

---

## Things the next chat will need to be aware of

- User's working pattern: "do not give me hacky fixes — fully reason
  any proposed fix and assess it and retry again if not 9/10
  confident." Reading this and following it is non-negotiable.
- Files the user has confirmed paths for verbatim:
  - `/Volumes/HPC-Home/receptor_design/receptor-resurfacing-pipeline/tests/rfdiffusion/receptor_resurfacing_results/rfdiffusion/split/`
  - `/hpc-home/jowillia/singularity/AlphaFold3/af3.bin`
  - `$HOME/af3_db` exists (AF3 database)
  - `.slurm.sh` extension (NOT `_slurm.sh`)
- The previous-session journal lives at `journal.txt` in
  `/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/`.
- Long previous transcripts are at `/mnt/transcripts/` and can be
  read incrementally if a specific historical decision needs to be
  re-derived.

---

## Closing note

v8 / pipeline_notes9 is a stabilisation revision — no new methodology,
no new metrics, just running what was staged in v7 and fixing what
broke. The metric refactor is the substantive architectural change:
it cleans up a sidecar-path-resolution issue that would have re-
surfaced repeatedly, and it surfaces all per-prediction metrics at
the same aggregation level. Once the overnight runs land cleanly,
v9 can move on to plots.
