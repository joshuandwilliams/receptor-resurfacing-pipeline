# tests/orthogonal_metrics/data/

Canonical inputs for `test_orthogonal_metrics.nf`. Two per-sequence
workdirs from the discovery run plus a trimmed
`cross_sequence_summary.csv`, with all absolute HPC paths rewritten to
the local staging location by `scripts/rewrite_fixture_paths.py`.

## Layout

```
data/
└── negsteer_run/
    ├── cross_sequence_summary.csv        <- 2 rows: design_28_seq_1, design_62_seq_0
    └── runs/
        ├── design_28_seq_1/               <- Class 1 (cold_start_all_clean), tier A
        │   ├── aggregated_results.csv
        │   ├── passing_summary.csv
        │   ├── pathways.json
        │   ├── all_results_multicycle{,_with_metrics}.csv
        │   ├── raw_per_seed_results.csv
        │   ├── run_one_runtime_sec.txt
        │   ├── inputs/{receptor,effector}.fasta
        │   └── cycle_0/
        │       ├── plan.json              (paths rewritten to local)
        │       ├── effector_template.cif
        │       ├── initial_prediction[_s1, _s2].pdb
        │       ├── contaminated.json / reversion_*.json
        │       ├── kickoff_distances.json / prefilter.json
        │       ├── steered_results.csv
        │       └── initial[, _s1, _s2]/input/input.yaml
        └── design_62_seq_0/               <- Class 2 (steered_clean tiered), tier A
            └── ...                        same shape, plus cycle_0/steered/<N>/prediction.pdb tree
```

## Sequence → outcome class mapping

| Sequence | Class | Cross-rank | Why curated |
|---|---|---|---|
| `design_28_seq_1` | 1 — `cold_start_all_clean` | 1 | Top representative. Exercises the survivor-manifest extraction on a representative whose `canonical_pdb` is `cycle_0/initial_prediction.pdb`. |
| `design_62_seq_0` | 2 — `steered_clean` tiered | n/a | Closest-to-passing on orthogonal metrics (`af3_nomsa_best_ra_eff=14.53`, only one orthogonal flag in the discovery run). Exercises the survivor-manifest extraction on a representative whose `canonical_pdb` is `cycle_0/steered/design_02_s0/prediction.pdb` (different code path from Class 1's `initial_prediction.pdb`). |

The pair covers both representative-PDB shapes (cold-start vs steered)
that `bin/extract_survivor_manifest.py` resolves.

## File-by-file

### `negsteer_run/cross_sequence_summary.csv`

The cohort table emitted by `NEGSTEER_CROSS_SEQUENCE`, **trimmed to
the two curated sequences**. Drives `EXTRACT_SURVIVOR_MANIFEST`.

The `rep_canonical_pdb` column has been **rewritten** to
point at the local staging:
- HPC `/hpc-home/.../tests/full_test_run/results/negative_steering/runs/<seq>/...`
- Local `<repo>/tests/orthogonal_metrics/data/negsteer_run/runs/<seq>/...`

### `negsteer_run/runs/<seq>/`

Per-sequence workdirs published by `NEGSTEER_RUN_ONE`. Each carries
the data the orthogonal-metrics cascade needs:

- `cycle_0/plan.json` — `ground_truth`, `effector_template_cif`, and
  `cold_start_seeds[*].pdb_path` rewritten to local paths.
  `boltz_container`, `receptor_fasta_override`, `effector_fasta_override`,
  and every `designs[*].dir` / `designs[*].yaml` blanked out (these
  pointed at HPC scratch / a Singularity image; the orthogonal-metrics
  test cascade does not read them).
- `cycle_0/effector_template.cif` — used by `AF3_NOMSA_ON_SURVIVORS`.
- `cycle_0/initial_prediction[_s<N>].pdb` — the cold-start prediction
  PDB(s).
- `cycle_0/steered/<label>/prediction.pdb` (for `design_62_seq_0` only) —
  the steered prediction the cohort representative points at.
- `aggregated_results.csv` — read by `compute_interface_metrics.py`.

### Tarball excludes (HPC-side)

These exclusion patterns were applied when building the source tarball:

- `*.npz` — Boltz internal weight checkpoints; not read by any
  orthogonal_metrics module (verified).
- `*.a3m`, `msa/` — MSA inputs; not used by AF3-no-MSA stream.
- `predictions/` — raw Boltz outputs; the cascade only reads PDBs via
  the `canonical_pdb` path.
- `contamination_scratch/` — intermediate files from the
  build-contaminated stage.

## Path-rewriting

The discovery-run outputs contain absolute HPC paths in
`cycle_0/plan.json` and the cohort CSV. `scripts/rewrite_fixture_paths.py`
(stdlib-only, idempotent) was applied at staging time to rewrite those
to the local Mac/repo paths. If you re-pull the fixture from HPC, re-run:

```bash
python3 scripts/rewrite_fixture_paths.py \
    --staged-runs-dir tests/orthogonal_metrics/data/negsteer_run/runs \
    --cross-csv       tests/orthogonal_metrics/data/negsteer_run/cross_sequence_summary.csv \
    --design-pdbs-dir tests/negative_steering/data/design_pdbs \
    --hpc-prefix      /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/full_test_run/results
```

The rewriter points `plan["ground_truth"]` at
`tests/negative_steering/data/design_pdbs/design_<N>.pdb`. The two
fixtures share that directory rather than duplicating the parent
design PDBs into both `data/` trees — the per-module-test isolation
contract is preserved because the path is absolute (not a relative
import) and either fixture can be re-staged independently with the
rewriter pointing at a different `--design-pdbs-dir`.

## Path coverage

See `notes/inventory/15_discovery_run_path_coverage.md` §Orthogonal
metrics for the full path-coverage analysis. Note Coverage gap 3:
production thresholds yield 0 passes on the discovery cohort because
`af3_nomsa_ra_eff` is a warning-style gate that AF3-no-MSA's accuracy
floor regularly trips. The pass branch of `merge_orthogonal_metrics.py`
is exercisable via param override (e.g. `--orthogonal_filter_af3_ra_max=50.0`)
without changing fixtures.

## Source

Discovery run `tests/full_test_run/results/`, discovery date
2026-05-02, params per `params_full_test.yml`. Path-rewriter applied
on staging.
