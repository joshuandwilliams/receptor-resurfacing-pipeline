# 03 — Dependency Graph

Two sections: (a) which Python files in `bin/` import which other Python files in `bin/`, and (b) which `.nf` modules invoke which `bin/` scripts. A third section calls out subprocess-mediated dependencies (Python files that exec each other via `subprocess.run`), which the import graph alone does not capture.

---

## A. Python imports inside `bin/`

Direct `import X` and `from X import …` edges, restricted to other files in `bin/`. Standard-library and third-party imports omitted.

```
boltz2_negative_steering.py        ← compute_metrics.py    (compute_effector_interface_residues; lazy, inside cmd_predict_one)
boltz2_iterate_steering.py         ← boltz2_negative_steering.py   (lazy import of whole module + named: get_chain_sequence, run_boltz)
build_control_sequences.py         ← boltz2_negative_steering.py   (get_chain_sequence)
collect_haddock3_dock.py           ← haddock_utils.py
compute_interface_metrics.py       ← compute_metrics.py    (compute_ipsae)
cross_sequence_summary.py          ← extract_passing.py    (extract_row)
derive_input_design_region.py      ← compute_metrics.py    (find_contact_residues_heavy; lazy, inside _compute_true_interface_0b)
                                   ← boltz2_negative_steering.py   (get_chain_sequence)
extract_hotspots.py                ← haddock_utils.py
haddock3_plots.py                  ← haddock_utils.py
reversion.py                       ← boltz2_negative_steering.py   (write_boltz_yaml; binding_rmsds, find_contact_residues_heavy, jaccard)
rfdiffusion_contigs.py             ← contig_utils.py
rfdiffusion_filter.py              ← contig_utils.py
```

### "Hub" / "leaf" summary

```
                     ┌────────────────────────────────────────┐
                     │   boltz2_negative_steering.py  (HUB)    │
                     │   imported by 4 files; itself imports   │
                     │   compute_metrics                       │
                     └────────────────────────────────────────┘
                              ▲           ▲           ▲
                              │           │           │
            boltz2_iterate_steering   reversion   build_control_sequences
                                                       │
                                                       └─ derive_input_design_region

                ┌─────────────────────────────────┐
                │   compute_metrics.py            │
                │   imported by 3 files; also     │
                │   invoked via subprocess (×2)   │
                └─────────────────────────────────┘
                  ▲                ▲
                  │                │
   boltz2_negative_steering   compute_interface_metrics, derive_input_design_region

                ┌─────────────────────────────────┐
                │   haddock_utils.py              │
                │   imported by 3 files           │
                └─────────────────────────────────┘
                  ▲           ▲              ▲
            collect_haddock3_dock   haddock3_plots   extract_hotspots

                ┌─────────────────────────────────┐
                │   contig_utils.py               │
                │   imported by 2 files           │
                └─────────────────────────────────┘
                  ▲                          ▲
        rfdiffusion_contigs          rfdiffusion_filter

                ┌─────────────────────────────────┐
                │   extract_passing.py            │
                │   imported by 1 file            │
                └─────────────────────────────────┘
                  ▲
        cross_sequence_summary
```

### Files in `bin/` that nothing else in `bin/` imports

These are either Nextflow-direct entry points (one-shot CLIs) or genuine orphans:

| File | Status |
|---|---|
| `boltz2_iterate_steering.py` | Top-level; invoked by negative_steering_run_one.sh + as subprocess from itself across kickoff steps. Not imported. |
| `boltz2_negative_steering.py` | Top-level CLI; **also a hub library** (4 importers). |
| `build_contigs.py` | Top-level CLI for `BUILD_CONTIGS`. |
| `compute_interface_metrics.py` | Top-level CLI for `NEGSTEER_INTERFACE_METRICS`. |
| `compute_metrics.py` | Top-level CLI; **also a hub library** (3 importers + 2 subprocess callers). |
| `cross_sequence_summary.py` | Top-level CLI for `NEGSTEER_CROSS_SEQUENCE`. |
| `derive_design_region.py` | Top-level CLI for `NEGSTEER_DERIVE_INDICES`. |
| `derive_input_design_region.py` | Top-level CLI for `DERIVE_INPUT_INDICES`. |
| `derive_true_interface.py` | Top-level CLI for `NEGSTEER_DERIVE_INDICES`. |
| `extract_hotspots.py` | Top-level CLI for `EXTRACT_HOTSPOTS`. |
| `extract_passing.py` | CLI; **also imported once** by `cross_sequence_summary.py`. |
| `extract_survivor_manifest.py` | Top-level CLI for `EXTRACT_SURVIVOR_MANIFEST`. |
| `merge_orthogonal_metrics.py` | Top-level CLI for `NEGSTEER_ORTHOGONAL_METRICS`. |
| `mpnn_*.py` (6 files) | All top-level CLIs for the MPNN module's processes. |
| `parse_af3_output.py` | Top-level CLI for `AF3_PARSE_OUTPUT`. |
| `pipeline_correct_sequences.py` | Top-level CLI for `SEQUENCE_CORRECTION` + `MPNN_FIXED_POSITIONS`. |
| `reversion.py` | Top-level (called via subcommands inside boltz2_iterate_steering's harvest-reversions); also imports from boltz2_negative_steering. |
| `rfdiffusion_contigs.py` | Top-level CLI for `RESOLVE_CONTIGS`. |
| `rfdiffusion_filter.py` | Top-level CLI for `RFDIFFUSION_FILTER`. |
| `rfdiffusion_plots.py` | Top-level CLI for `RFDIFFUSION_PLOTS`. |
| `rosetta_filter_collect.py` | Top-level CLI for `ROSETTA_FILTER`. |
| `rosetta_filter_plots.py` | Top-level CLI for `ROSETTA_FILTER_PLOTS`. |
| `run_biophysical_metrics.py` | Top-level CLI for `NEGSTEER_BIOPHYSICAL_METRICS`. |
| `run_rosetta_metrics.py` | Top-level CLI for `NEGSTEER_ROSETTA_METRICS`. |
| `negsteer_plots.py` / `negsteer_within_sequence_plots.py` / `orthogonal_metrics_plots.py` | Top-level CLIs for the three plot processes. |
| **`sequence_registry.py`** | **Genuine orphan: zero importers, zero invocations from .nf or .sh, last referenced in pipeline_notes9 as "confirmed dead code".** |
| `haddock3_prepare.py` | Top-level CLI; not imported. |

---

## B. Nextflow → bin invocations

Each `.nf` process either runs Python via `python ${projectDir}/bin/<script>.py` or shells out via `bash ${projectDir}/bin/<script>.sh`. Listed by source file.

### main.nf
- Pure orchestration (no direct script calls of its own — all script invocations happen inside the included process bodies). Passes `Channel.value(file("${projectDir}/bin/parse_af3_output.py"))`, `Channel.value(file("${projectDir}/bin/orthogonal_metrics_plots.py"))`, `Channel.value(file("${projectDir}/bin/fastrelax_for_ia.xml"))` into downstream processes (Bug 7 cache-invalidation pattern from pipeline_notes14).

### modules/preprocessing.nf
- `EXTRACT_SEQUENCES` — inline Python (no bin/ script).
- `RESOLVE_CONTIGS` → `bin/rfdiffusion_contigs.py`.
- `WRITE_DUMMY_MAPPING` — inline `printf '{}'`.

### modules/haddock.nf
- `HADDOCK3_PREPARE` → `bin/haddock3_prepare.py`.
- `HADDOCK3_DOCK` → `bin/collect_haddock3_dock.py`.
- `HADDOCK3_PLOTS` → `bin/haddock3_plots.py`.
- `EXTRACT_HOTSPOTS` → `bin/extract_hotspots.py`.
- `BUILD_CONTIGS` → `bin/build_contigs.py`.

### modules/rfdiffusion.nf
- `RFDIFFUSION` — inline RFDiffusion invocation (no bin/ script).
- `RFDIFFUSION_FILTER` → `bin/rfdiffusion_filter.py`.
- `RFDIFFUSION_PLOTS` → `bin/rfdiffusion_plots.py`.

### modules/rosetta_filtering.nf
- `ROSETTA_SC` — inline Rosetta InterfaceAnalyzer (no bin/ script).
- `ROSETTA_FILTER` → `bin/rosetta_filter_collect.py`.
- `ROSETTA_FILTER_PLOTS` → `bin/rosetta_filter_plots.py`.

### modules/proteinmpnn.nf
- `MPNN_FIXED_POSITIONS` → `bin/pipeline_correct_sequences.py`.
- `PROTEINMPNN` — inline ProteinMPNN invocation.
- `SEQUENCE_CORRECTION` → `bin/pipeline_correct_sequences.py`.
- `SEQUENCE_QC` → `bin/mpnn_sequence_qc.py`.
- `MPNN_DESIGN_REGION_SCORE` → `bin/mpnn_design_region_score.py`.
- `MPNN_CLUSTER` → `bin/mpnn_cluster_sequences.py`.
- `MPNN_SELECT_TOP` → `bin/mpnn_select_top.py`.
- `MPNN_PLOTS` → `bin/mpnn_plots.py`.

### modules/negative_steering.nf
- `NEGSTEER_DERIVE_INDICES` → `bin/derive_design_region.py` AND `bin/derive_true_interface.py`.
- `NEGSTEER_RUN_ONE` → `bash bin/negative_steering_run_one.sh` (which in turn runs `bin/boltz2_negative_steering.py`, `bin/boltz2_iterate_steering.py`, `bin/extract_passing.py` — see Section C).
- `NEGSTEER_CROSS_SEQUENCE` → `bin/cross_sequence_summary.py` (which transitively imports `bin/extract_passing.py`).
- `NEGSTEER_PLOTS` → `bin/negsteer_plots.py`.
- `NEGSTEER_WITHIN_SEQUENCE_PLOTS` → `bin/negsteer_within_sequence_plots.py`.

### modules/negsteer_controls.nf
- `DERIVE_INPUT_INDICES` → `bin/derive_input_design_region.py`.
- `NEGSTEER_CONTROLS` → `bin/build_control_sequences.py` AND `bash bin/negative_steering_run_one.sh` (per control variant).

### modules/negsteer_manifest.nf
- `EXTRACT_SURVIVOR_MANIFEST` → `bin/extract_survivor_manifest.py`.

### modules/negsteer_interface_metrics.nf
- `NEGSTEER_INTERFACE_METRICS` → `bin/compute_interface_metrics.py`.

### modules/negsteer_af3_nomsa.nf
- `AF3_SETUP_DB` — inline shell (no bin/ script).
- `AF3_NOMSA_ON_SURVIVORS` — inline AF3 invocation.
- `AF3_PARSE_OUTPUT` → `bin/parse_af3_output.py`.

### modules/negsteer_biophysical_metrics.nf
- `NEGSTEER_BIOPHYSICAL_METRICS` → `bin/run_biophysical_metrics.py`.

### modules/negsteer_rosetta_metrics.nf
- `NEGSTEER_ROSETTA_METRICS` → `bin/run_rosetta_metrics.py` (with `bin/fastrelax_for_ia.xml` as data input).

### modules/negsteer_orthogonal_metrics.nf
- `NEGSTEER_ORTHOGONAL_METRICS` → `bin/merge_orthogonal_metrics.py`.
- `ORTHOG_PLOTS` → `bin/orthogonal_metrics_plots.py`.

---

## C. Subprocess-mediated edges (NOT visible in import graph)

Several Python scripts invoke other `bin/` scripts via `subprocess.run`. These are real runtime dependencies but invisible to static import analysis:

| Caller | Subprocess target | Where |
|---|---|---|
| `bin/boltz2_iterate_steering.py` | `bin/compute_metrics.py` | inside `cmd_iterate_collect_finalize`, `cmd_harvest_reversions`, `cmd_compute_final_metrics` (see lines ~1413, ~2723, ~5331). Path resolved via `--compute-metrics-script` arg or `SCRIPT_DIR / "compute_metrics.py"`. |
| `bin/reversion.py` | `bin/compute_metrics.py` | inside `harvest_reversion_results` (line ~751). Same arg pattern. |
| `bin/boltz2_iterate_steering.py` | `boltz` (Boltz binary) | via the imported `run_boltz` from `boltz2_negative_steering`. |
| `bin/negative_steering_run_one.sh` | `bin/boltz2_negative_steering.py`, `bin/boltz2_iterate_steering.py`, `bin/extract_passing.py` | the orchestrator's stage chain. |
| `bin/run_rosetta_metrics.py` | rosetta_scripts / InterfaceAnalyzer binaries | via `_run_fast_relax` / `_run_interface_analyzer`. |
| `bin/mpnn_cluster_sequences.py` | `mmseqs` binary | via `run_mmseqs_cluster`. |

---

## D. Test-tree → bin/ relationships

Test workflows include the production process bodies directly (so they exercise the same code paths) — no test-only Python is wired into production. But two test Python files duplicate production scripts:

- `tests/orthogonal_metrics/merge_orthogonal_metrics.py` is a divergent copy of `bin/merge_orthogonal_metrics.py` (see `05_findings.md`).
- `tests/proteinmpnn/test_mpnn_plots.py` ↔ `bin/mpnn_plots.py`, `tests/rfdiffusion/test_rfdiffusion_plots.py` ↔ `bin/rfdiffusion_plots.py`, `tests/rosetta_filtering/test_rosetta_filtering_plots.py` ↔ `bin/rosetta_filter_plots.py`, `tests/negative_steering/test_negsteer_plots.py` ↔ `bin/negsteer_plots.py`, `tests/negative_steering/test_negsteer_within_sequence_plots.py` ↔ `bin/negsteer_within_sequence_plots.py`, `tests/orthogonal_metrics/test_orthogonal_metrics_plots.py` ↔ `bin/orthogonal_metrics_plots.py`. These pairs are the documented "test-iterates-fast, production-mirrors-after-signoff" workflow established in pipeline_notes11.

---

## E. Cycle / coupling notes

**No import cycles**: the lazy-import pattern in `boltz2_iterate_steering.py` (`_lazy_bns()` and per-function `from boltz2_negative_steering import …`) is in place explicitly to avoid an import-time round-trip with `boltz2_negative_steering` (which pulls heavy third-party packages and is slow to import). With those guards, the static graph is a DAG.

**The `boltz2_negative_steering.py` hub problem**: 4 other bin/ files import named symbols from it. That makes the file's exported API surface implicit (no `__all__`), so any rename or signature change ripples — and because `boltz2_iterate_steering.py` does `from boltz2_negative_steering import get_chain_sequence` and `import run_boltz` from inside specific function bodies (lines 2560 + 2644), the rippling is invisible at the top of the file.

**The `compute_metrics.py` dual-mode problem**: it's both a library (3 importers via `from compute_metrics import …`) and a subprocess target (2 callers via `subprocess.run`, with the path passed as a CLI flag `--compute-metrics-script`). This means rename or move of the file requires touching not just import sites but also CLI default paths in two callers' argparse setups.
