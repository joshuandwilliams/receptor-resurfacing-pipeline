# experiments/

Exploratory work for designing and analysing binder-design campaigns. Sits
alongside the production Nextflow pipeline (`main.nf`, `bin/`, `modules/`)
but does not feed back into it — nothing under `experiments/` is invoked by
the pipeline, and the pipeline does not depend on anything here.

## Layout

- `scripts/` — reusable helpers shared across campaigns. Importable Python
  package (has `__init__.py`).
- `param_derivation/` — code for deriving pipeline input parameters,
  particularly RFDiffusion contigs. Kept separate from `scripts/` because
  it is a substantial, distinct concern.
- `inputs/` — shared input assets used across campaigns (PDBs, UniProt
  FASTAs, reference structures). See `inputs/README.md` for naming and
  provenance conventions.
- `campaigns/` — one subdirectory per campaign. Each campaign holds both
  its driving script and its outputs (`results/`, `analyses/`). See
  `campaigns/README.md`.

## Mac-authoritative workflow

The repo follows a Mac-authoritative, HPC-runs workflow (see
`notes/remediation_state.md`). Edits happen on the Mac and propagate to
HPC via `scripts/sync_to_hpc.sh`. Campaign scripts that submit SLURM jobs
should assume the same: author on Mac, sync, run on HPC.

## Imports

Campaign scripts will need to import from two places:

1. `experiments/scripts/` — campaign-shared helpers.
2. `bin/` — production helpers like `compute_metrics`, contact-residue
   logic, etc.

`bin/` is not a Python package and the repo root is not on `sys.path` by
default. A future prompt will introduce a small bootstrap module under
`experiments/scripts/` that puts both locations on `sys.path` so campaign
scripts can do `from experiments.scripts import …` and
`from compute_metrics import …` uniformly. Until that lands, campaign
scripts should not be authored — defer to the bootstrap prompt.

## Relationship to `tests/`

`tests/` holds per-module characterisation tests for the production
pipeline. `experiments/` holds exploratory design work. They do not share
fixtures or helpers; do not cross-import.
