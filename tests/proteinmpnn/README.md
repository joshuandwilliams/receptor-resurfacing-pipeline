# tests/proteinmpnn/

Per-module test for the ProteinMPNN sequence-design stage.

## What it exercises

The full MPNN block from `modules/proteinmpnn.nf` plus the upstream
preprocessing helpers (`EXTRACT_SEQUENCES`, `RESOLVE_CONTIGS` from
`modules/preprocessing.nf`):

`MPNN_FIXED_POSITIONS` → `PROTEINMPNN` → `SEQUENCE_CORRECTION` →
`SEQUENCE_QC` → `MPNN_DESIGN_REGION_SCORE` → `MPNN_CLUSTER` →
`MPNN_PLOTS` → (optional) `MPNN_SELECT_TOP`.

Mirrors the Branch B path in `main.nf` end-to-end up through MPNN.

## Inputs (from `data/`)

| Path | What it is |
|---|---|
| `data/input_complex.pdb` | Receptor+effector complex PDB. Used by `RESOLVE_CONTIGS` and `EXTRACT_SEQUENCES` to derive sequences and `receptor_start_pdb`. |
| `data/rosetta_passing/design_*.pdb` | Rosetta-survivor split-PDB backbones, ready for MPNN. |

See `data/README.md` for the per-file contract.

## Outputs (under `results/`)

- `proteinmpnn/<design>/` — raw MPNN outputs per design.
- `sequences/qc_fastas/design_<N>_seq_<M>.fasta` — QC-passing
  MPNN-corrected FASTAs (always produced).
- `sequences/top_fastas/design_<N>_seq_<M>.fasta` — top-N selected
  FASTAs (only when `params.mpnn_top_n > 0`).
- `mpnn_metadata/*.csv` — sequence-correction metadata, design-region
  scores, cluster counts.
- `mpnn_plots/*.png` — diagnostic plots.

The `qc_fastas/` and `top_fastas/` directories are the canonical
upstream feed for the negative-steering test once fixtures are curated.

## Running on HPC

```bash
sbatch tests/proteinmpnn/run_test_proteinmpnn.slurm.sh
```
