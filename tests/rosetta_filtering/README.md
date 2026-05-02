# tests/rosetta_filtering/

Per-module test for the Rosetta filtering stage.

## What it exercises

`ROSETTA_SC` → `ROSETTA_FILTER` → `ROSETTA_FILTER_PLOTS` from
`modules/rosetta_filtering.nf`, in isolation. Computes the Rosetta
shape-complementarity score (Sc) for each input split-PDB, partitions
designs into `passing/` vs `failing/` against `params.sc_threshold`,
and emits diagnostic plots.

## Inputs (from `data/`)

| Path | What it is |
|---|---|
| `data/rfdiffusion_split/design_*.pdb` | Split two-chain PDBs (chain A receptor, chain B effector) — the same shape as `RFDIFFUSION_FILTER`'s `split/` output. |

See `data/README.md` for the per-file contract.

## Outputs (under `results/`)

- `rosetta_filtering/scored/design_*.pdb` — every input PDB after
  Rosetta scoring (with the score sidecar).
- `rosetta_filtering/passing/design_*.pdb` — designs with Sc above
  threshold; canonical input feed for the MPNN test.
- `rosetta_filtering/failing/design_*.pdb` — below-threshold designs.
- `rosetta_filtering/rosetta_filter_metrics.json` — aggregate metrics.
- `rosetta_filter_plots/*.png` — diagnostic plots.

## Running on HPC

```bash
sbatch tests/rosetta_filtering/run_test_rosetta_filtering.slurm.sh
```
