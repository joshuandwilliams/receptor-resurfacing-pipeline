# tests/negative_steering/

Per-module test for the negative-steering validation stage.

## What it exercises

The three processes in `modules/negative_steering.nf`:

`NEGSTEER_DERIVE_INDICES` → `NEGSTEER_RUN_ONE` → `NEGSTEER_CROSS_SEQUENCE`

plus the diagnostic plots (`NEGSTEER_PLOTS`,
`NEGSTEER_WITHIN_SEQUENCE_PLOTS`) and, when
`params.run_negative_controls = true` (default), the controls path
from `modules/negsteer_controls.nf` (`DERIVE_INPUT_INDICES`,
`NEGSTEER_CONTROLS`).

For each MPNN sequence: one GPU job runs the full single-cycle chain
end-to-end (cold-start prediction → identify wrong-interface →
mutate → re-predict → revert protected-set mutations → re-predict
→ harvest verdicts). Outputs are aggregated into the cohort
`cross_sequence_summary.csv`.

## Inputs (from `data/`)

| Path | What it is |
|---|---|
| `data/input_complex.pdb` | Receptor+effector complex (controls path only). |
| `data/rfdiffusion_metrics.json` | Per-design metrics blob from `RFDIFFUSION_FILTER`. |
| `data/design_pdbs/design_*.pdb` | Parent RFDiffusion design PDBs. |
| `data/fastas/design_<N>_seq_<M>.fasta` | MPNN-corrected FASTAs (one per cohort row). |

See `data/README.md` for the per-file contract and the regex the
fasta filenames must satisfy.

## Outputs (under `results/`)

- `negative_steering/runs/<seq_name>/` — per-sequence workdirs
  containing `plan.json`, `effector_template.cif`, `cycle_0/...`.
- `negative_steering/cross_sequence_summary.csv` — cohort table
  (drives the orthogonal-metrics test).
- `negative_steering/passing_summary.csv` — tier A/B/C survivors.
- `negsteer_plots/*.png` — cohort-level plots.
- `negsteer_within_sequence_plots/*.png` — per-sequence plots.

## Running on HPC

The `bin/` directory must be symlinked into this test dir so the
module can find it (the symlink is part of the committed repo
state):

```bash
ln -s ../../bin tests/negative_steering/bin   # one-time, if missing
sbatch tests/negative_steering/run_test_negative_steering.slurm.sh
```
