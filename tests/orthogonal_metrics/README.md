# tests/orthogonal_metrics/

Per-module test for the orthogonal-metrics validation stage.

## What it exercises

The full post-negsteer cascade:

`NEGSTEER_INTERFACE_METRICS` → `EXTRACT_SURVIVOR_MANIFEST` →
{ `AF3_SETUP_DB` → `AF3_NOMSA_ON_SURVIVORS` → `AF3_PARSE_OUTPUT`,
  `NEGSTEER_BIOPHYSICAL_METRICS`,
  `NEGSTEER_ROSETTA_METRICS` } →
`NEGSTEER_ORTHOGONAL_METRICS`

Three independent-method validators (AF3-no-MSA, biophysical,
Rosetta) run in parallel on cohort survivors with a representative
PDB on disk. The merge step combines the three streams and computes
`passes_orthogonal_filters`.

## Inputs (from `data/`)

| Path | What it is |
|---|---|
| `data/negsteer_run/cross_sequence_summary.csv` | Cohort table from a prior negsteer fixture run. |
| `data/negsteer_run/runs/<seq_name>/` | Per-sequence workdirs: `plan.json`, `effector_template.cif`, `cycle_0/...`. |

See `data/README.md` for the per-workdir contract. The shape mirrors
`tests/negative_steering/results/negative_steering/`.

## Outputs (under `results/`)

- `negsteer_interface_metrics/cross_sequence_summary_extended.csv` —
  cohort table with iRMSD/fnat/DockQ/iPSAE_15/intact_core/weighted_jaccard
  appended.
- `negsteer_manifest/survivor_manifest.csv` — fan-out manifest.
- `af3_nomsa/<seq_name>/...` — AF3 predictions per survivor.
- `negsteer_biophysical_metrics/<seq_name>/summary.csv` — BSA /
  interface_plddt / interface_hbonds.
- `negsteer_rosetta_metrics/<seq_name>/summary.csv` — Sc / dG_separated
  / dSASA_int / ΔΔG.
- `negsteer_orthogonal_metrics/cross_sequence_summary_orthogonal.csv` —
  final merged + filtered cohort table.

## Running on HPC

```bash
sbatch tests/orthogonal_metrics/run_test_orthogonal_metrics.slurm.sh
```
