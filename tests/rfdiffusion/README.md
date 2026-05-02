# tests/rfdiffusion/

Per-module test for the RFDiffusion stage.

## What it exercises

`RFDIFFUSION` → `RFDIFFUSION_FILTER` → `RFDIFFUSION_PLOTS` from
`modules/rfdiffusion.nf`, in isolation. Generates N candidate
backbones around the fixed effector pose, filters by hotspot
occupancy, splits each surviving design into chain-A receptor +
chain-B effector PDBs, and emits the per-design metrics blob plus
diagnostic plots.

## Inputs (from `data/`)

| Path | What it is |
|---|---|
| `data/af3_pikp1_native_avrpikf_complex.pdb` | Pre-docked receptor+effector complex PDB (chain A receptor, chain C effector). |

The contigs string and hotspot are wired into `test_rfdiffusion.nf`
directly. See `data/README.md` for details.

## Outputs (under `results/`)

- `rfdiffusion/raw/design_*.pdb` — every backbone generated.
- `rfdiffusion/split/design_*.pdb` — surviving designs split into
  receptor + effector chains.
- `rfdiffusion/passing/design_*.pdb` — designs above the
  hotspot-occupancy threshold.
- `rfdiffusion/rfdiffusion_metrics.json` — per-design metrics
  (`receptor_contact_residues`, `receptor_position_order`,
  `per_design_design_residues`, …).
- `rfdiffusion_plots/*.png` — diagnostic plots.

The `split/` directory is the canonical input feed for downstream
per-module tests (Rosetta filtering, MPNN, negative steering) once
fixtures are curated.

## Running on HPC

```bash
sbatch tests/rfdiffusion/run_test_rfdiffusion.slurm.sh
```
