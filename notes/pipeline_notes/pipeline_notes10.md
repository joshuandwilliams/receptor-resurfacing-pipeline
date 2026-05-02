# Pipeline Notes 10 — Session 26 Apr 2026

Orthogonal-metrics test closed; first end-to-end production run staged for resume.

## Headline

Two big things landed this session.

- Orthogonal-metrics test (Step 2 of the 1–15 plan) is now closed. AF3, Rosetta, biophysical, DockQ all populate cleanly for the two test survivors.
- First production run on the PikP1 + AvrPikF target ran end-to-end through negsteer, but produced an empty survivor manifest because of a chain-param wiring bug. Patches written and staged; resume blocked on GPU queue contention.

Net effect: pipeline is structurally complete from input PDB through orthogonal metrics. Step 3 (plot scripts, Tasks 43+44) is unblocked and starts in the next session.

## Bug cascade resolved (in order)

### 1. AF3 model_dir path

`nextflow.config:28` had `af3_model_dir = "${HOME}/af3_models"` — wrong path. Real location is `/hpc-home/jowillia/singularity/AlphaFold3/` (contains `af3.bin`). Patched config at `/mnt/user-data/outputs/nextflow.config`.

### 2. Effector-template mmCIF generation (boltz2_negative_steering.py)

Three sub-bugs in `extract_effector_template_cif` surfaced sequentially in AF3:

- **(a)** Missing `_entity_poly_seq` block + `?` for `pdbx_seq_one_letter_code` + `.` for every `label_seq_id`. Root cause: `setup_entities()` identifies polymer chains but does NOT populate `Entity.full_sequence` or `Residue.label_seq` from PDB input. Fix: prune-and-rename approach (start from parsed structure, remove unwanted chains rather than building from scratch); explicitly populate `entity.full_sequence = [r.name for r in polymer]` and `r.label_seq = i`.
- **(b)** HETATM water rows in effector chain (controls' input PDB has waters in chain C) caused validation to fail. Fix: use gemmi's built-in `st.remove_waters()`.
- **(c)** AF3's `templates.py:get_polymer_features` requires `release_date`. AF3 reads ONLY `_pdbx_audit_revision_history.revision_date` (confirmed by reading source: `alphafold3/structure/mmcif.py:225` — function `get_release_date`). Initial fix appended raw text `loop_` block after `gemmi.write` but parser didn't see it. Final fix: use `block.init_mmcif_loop("_pdbx_audit_revision_history.", [...])` to add loop INSIDE the CIF document, then `loop.add_row(["1", "'Structure model'", "1", "0", "2024-01-01"])`.

Validation block re-parses written CIF with `gemmi.cif.read()` and confirms `find_values()` returns the date — same lookup AF3 uses. Patched file at `/mnt/user-data/outputs/boltz2_negative_steering.py`.

### 3. Regenerated CIFs in 22 negsteer work directories

Avoided full negsteer rerun (~12 GPU-hr) by regenerating all `effector_template.cif` files in place. Used `singularity exec /hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/boltz2_negsteer.img /usr/local/bin/python3.11` (NOT `/usr/bin/python3` — that lacks numpy/gemmi). 18 design CIFs regenerated using rfdiffusion split outputs as source with chain B; 4 control CIFs regenerated using `tests/negative_steering/receptor_resurfacing_results/negative_steering/controls_inputs/af3_pikp1_native_avrpikf_complex.pdb` with chain C→B remap.

Verification: all CIFs have 82 `_entity_poly_seq` rows, full one-letter sequence, zero `.` in `label_seq_id`.

### 4. Nextflow input file collision in NEGSTEER_ORTHOGONAL_METRICS

Each upstream module emitted fixed-name CSVs (`af3_nomsa_summary.csv`, `biophysical_summary.csv`, `rosetta_summary.csv`); collecting across 2 survivors caused name collision. Fix: tag outputs with `${seq_name}`. Three modules patched; `merge_orthogonal_metrics.py` uses globs (`*.csv`) so no script change needed.

### 5. Rosetta InterfaceAnalyzer crash (exit 255)

`ROSETTA_CRASH.log`: "Option matching `-compute_separated_dG` not found in command line top-level context". Rosetta 2025.37 dropped that flag (now default behaviour). Fix in `bin/run_rosetta_metrics.py`: mirror working invocation from `modules/rosetta_filtering.nf:ROSETTA_SC`. Added `-database /opt/rosetta/main/database`; dropped `-compute_separated_dG`, `-compute_packstat`, `-tracer_data_print`, `-out:file:silent_struct_type binary`, `-unmute`, `-interface A_B`; used `-pack_separated -pack_input -compute_interface_sc true -add_regular_scores_to_scorefile -use_input_sc -no_optH false -ignore_unrecognized_res -out:file:score_only`.

### 6. Biophysical script (FreeSASA + H-bond)

Initial errors swallowed by broad `except` clauses: `freesasa_failed:AttributeError` and `hba_run_failed:NoDataError`. Added traceback printing to surface real errors:

- **FreeSASA**: `module 'freesasa' has no attribute 'Calc'` — freesasa>=2 dropped the `Calc` class. Fix: `freesasa.calc(structure)` (lowercase, top-level function).
- **H-bonds**: `MDAnalysis NoDataError: This Universe does not contain charge information` — Boltz PDBs are heavy-atom only, no hydrogens. Fix: replaced MDAnalysis HBA with gemmi-based heavy-atom geometric detector (N/O atoms within 3.5Å between chains). Standard fallback used by FreeContact/PISA/BioPython NeighborSearch.

## Orthogonal-metrics test result

Both real survivors produce clean numbers. Pipeline plumbing closed.

| Metric | Survivor 1 | Survivor 2 |
|---|---|---|
| `bsa` (Å²) | 889.76 | 801.66 |
| `interface_hbonds` | 11 | 9 |
| `sc` (Lawrence-Colman) | 0.628 | 0.693 |
| `rosetta_ddg` | -15.85 | +24.12 |
| `af3_nomsa_best_ra_eff` (Å) | 37.4 | 36.4 |

All `*_failures` columns empty. `passes_orthogonal_filters = 0` for both survivors — design-quality result, not pipeline bug. AF3 ra_eff ~37Å is far above the 5Å pass threshold, expected since negsteer steered against Boltz biases specifically.

## Production run (PikP1 + AvrPikF)

### Setup

- Path: `/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/full_test_run/`
- Input: `af3_pikp1_native_avrpikf_complex.pdb` (mode 2, no HADDOCK)
- Contigs: `A1-32/10-20/A46-72/6-6 C32-113`
- `num_designs=8`, `num_seqs=16` → 128 MPNN sequences
- `mpnn_top_n=8` → 8 sequences progressed to negative steering
- SLURM job 19566777

### Outcome of first run attempt

Ran end-to-end through negsteer. Cross-sequence summary produced 10 rows:

- 1 Tier A row: `design_2_seq_13` (composite 0.599, ra_eff 2.29Å, intact_core=1)
- 1 Tier B row
- 8 "tier=none" rows (no passing predictions)

But: orthogonal-metrics step produced an empty `survivor_manifest.csv`. Root cause: chain-param wiring bug.

### The chain-param wiring bug

`extract_survivor_manifest.py` was being called with `--receptor-chain A --effector-chain C` (`params.receptor_chain` / `params.effector_chain` from production `params.yml`). Per the documented convention in `params_example.yml`, those describe the INPUT PDB. But `extract_survivor_manifest` reads from the prediction PDB, which uses A/B regardless (Boltz output convention).

Result: `_extract_chain_seq()` couldn't find chain C in the prediction PDB → both real survivors skipped with `seq_extraction_failed` → empty manifest → empty orthogonal_metrics output.

Why this slipped through earlier: the orthogonal-metrics test had `params.effector_chain="B"` set directly because that matched the test's prediction PDBs. Test passed even with the latent wiring bug. Production correctly sets `"C"` per the documented convention, which exposed the bug.

### Fix: switched four module files to use `rfdiff_output_*` params

Per `main.nf:680-695`, the `rfdiff_output_receptor_chain` (A) and `rfdiff_output_effector_chain` (B) params describe the prediction-PDB convention. All processes operating on prediction PDBs should use these, not `params.receptor_chain`/`effector_chain`. Patched four files:

- `modules/negsteer_manifest.nf` — line 55-56
- `modules/negsteer_af3_nomsa.nf` — line 267-268
- `modules/negsteer_biophysical_metrics.nf` — line 48-49
- `modules/negsteer_rosetta_metrics.nf` — line 43-44

All four staged at `/mnt/user-data/outputs/`. Drop into `modules/`, resume — content-keyed cache invalidates only those four processes (plus the merge step). All upstream stages cached.

### Resume command

```bash
cd /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/full_test_run
sbatch /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/run_pipeline.slurm.sh ./params.yml --resume
```

Currently blocked on GPU queue contention (someone else submitted a large job Sunday evening). When unblocked: resume should re-run EXTRACT_SURVIVOR_MANIFEST → AF3_PARSE_OUTPUT → BIOPHYSICAL → ROSETTA → ORTHOGONAL_METRICS only. Negsteer cache survives.

## Files staged at `/mnt/user-data/outputs/`

Latest version of each, ready to deploy:

- `nextflow.config` — AF3 model_dir fix (already deployed pre-production)
- `boltz2_negative_steering.py` — `extract_effector_template_cif` rewrite + validation
- `run_rosetta_metrics.py` — mirror working ROSETTA_SC invocation
- `run_biophysical_metrics.py` — FreeSASA API fix + gemmi H-bond detector
- `negsteer_manifest.nf` — chain-param fix (NEEDS DEPLOY before resume)
- `negsteer_af3_nomsa.nf` — chain-param fix + `seq_name` in filename (NEEDS DEPLOY)
- `negsteer_biophysical_metrics.nf` — same (NEEDS DEPLOY)
- `negsteer_rosetta_metrics.nf` — same (NEEDS DEPLOY)

## Next session — Step 3: plot scripts

Per the 1–15 execution plan, Step 3 is plot scripts (Tasks 43 + 44). Stated approach: develop Python scripts standalone against existing test output CSVs outside Nextflow; iterate freely; integrate into pipeline only when happy.

### Inputs available for plot development

- `cross_sequence_summary_with_interface_metrics.csv` (68 columns) — use for Task 43 (negsteer plots)
- `survivors_with_orthogonal_metrics.csv` (83 columns) — use for Task 44 (orthogonal-metrics plots)

### Task 43 — negsteer plots

- tier distribution (A/B/C/none counts)
- composite-score histogram
- ra_eff vs ipSAE scatter
- contamination rates
- control-row diagnostics (polyA, scrambled — flag controls passing filters)

### Task 44 — orthogonal-metrics plots

- AF3 ra_eff vs Boltz ra_eff scatter
- Sc histogram
- BSA histogram
- interface pLDDT
- ΔΔG distribution
- weighted_jaccard

## Environmental gotchas (re-statement for future Claude)

- Login node has stdlib only. Anything importing numpy, gemmi, freesasa, MDAnalysis MUST run inside the boltz container.
- Container has TWO pythons: `/usr/bin/python3` (system, no scientific stack) and `/usr/local/bin/python3.11` (has numpy/gemmi/etc — use this).
- AF3 image: `/software/e8edb411-7374-4342-b9f1-408da41fc197/e8edb411-7374-4342-b9f1-408da41fc197.img`. Internal package path: `/alphafold3_venv/lib/python3.11/site-packages/alphafold3/`. Source package via `source package e8edb411-...` does NOT make `alphafold3` importable from login-node python — must read source via singularity exec.
- Chain conventions: `params.effector_chain` describes input PDB; `params.rfdiff_output_effector_chain` describes split / prediction output (always B); receptor=A throughout.
- Production outdir: `${projectDir}/${params.project_name}_results` — NOT `${params.outdir}` from `params.yml` (config overrides).
- Production run launches from `tests/full_test_run/` — `work/` and `.nextflow.log` live there, not in project root.
