# Receptor Resurfacing Pipeline

A Nextflow DSL2 pipeline for *de novo* binder design: takes a docked
receptor–effector complex, redesigns the receptor surface around the
effector pose, and outputs candidate binder sequences with predicted
quality metrics from multiple orthogonal validators.

## Status

This codebase is undergoing structured remediation following
[`notes/codebase_remediation_plan.md`](notes/codebase_remediation_plan.md).
Active work happens on the `remediation` branch; `main` is the pristine
baseline. The single source of truth for *where the work currently
stands* is [`notes/remediation_state.md`](notes/remediation_state.md) —
read that at the start of every session.

## What the pipeline does

The input is a complex PDB containing a receptor (the protein to
redesign) and an effector (the target it must bind). Branch B takes a
pre-docked complex directly; Branch A docks a receptor + effector pair
with HADDOCK3 first.

From the complex, RFDiffusion generates candidate receptor backbones
around the fixed effector pose, Rosetta filters them by interface
shape complementarity, and ProteinMPNN designs amino-acid sequences for
the survivors. Each designed sequence is then validated through a
negative-steering Boltz2 stage (cold-start prediction → identify
wrong-interface contacts → mutate → re-predict → revert
contamination → re-predict) and graded into tiers A/B/C/none.

Cohort survivors are passed through three orthogonal validators
(AlphaFold3 no-MSA, biophysical metrics, Rosetta InterfaceAnalyzer)
plus a per-row interface-metrics step (iRMSD, fnat, DockQ, weighted
Jaccard). The output is a ranked cohort table of candidate binder
sequences with their quality metrics.

For the full stage-by-stage breakdown, see
[`notes/inventory/04_functional_categorization.md`](notes/inventory/04_functional_categorization.md).
For domain vocabulary, see
[`notes/inventory/06_ubiquitous_language.md`](notes/inventory/06_ubiquitous_language.md).

## Repo layout

```
bin/         Pipeline scripts (Python, called from Nextflow processes)
modules/     Nextflow modules — one per pipeline stage
tests/       Per-module tests + characterization tests + reference data
containers/  Singularity definition files for HPC execution
notes/       Inventory documents, remediation plan, decisions, glossaries
scripts/     Repo tooling (sync_to_hpc.sh, etc.)
main.nf      Top-level workflow
nextflow.config   Container paths and per-process SLURM resources
params_example.yml   Annotated parameter starter
run_pipeline.slurm.sh   SLURM submission wrapper
```

## Running the pipeline

The pipeline runs on HPC; the local Mac repo is for editing only.
Submit a run with:

```bash
sbatch run_pipeline.slurm.sh ./params_full_test.yml
```

Copy [`params_example.yml`](params_example.yml) as a starting point
for your own parameter file — it documents every parameter inline.
A reference small-scale params file lives at
[`tests/full_test_run/params_full_test.yml`](tests/full_test_run/params_full_test.yml).

## Running the tests

Two tiers:

```bash
# Local — comparator framework + helper unit tests (Mac, fast)
pytest -m local_unit

# HPC — characterization tests against a fresh pipeline run
RECEPTOR_OUTPUT_ROOT=<run-output-dir> pytest -m hpc
```

See [`tests/characterization/README.md`](tests/characterization/README.md)
for the full framing, marker tiers, and how to update reference
outputs after intentional behavior changes.

## Where to look for what

| Question | Document |
|---|---|
| How does the pipeline work? | [`notes/inventory/04_functional_categorization.md`](notes/inventory/04_functional_categorization.md) |
| What does *\<term\>* mean? | [`notes/inventory/06_ubiquitous_language.md`](notes/inventory/06_ubiquitous_language.md) |
| Where is the remediation at? | [`notes/remediation_state.md`](notes/remediation_state.md) |
| What's the testing strategy? | [`notes/inventory/14_phase_2_revision_per_module_tests.md`](notes/inventory/14_phase_2_revision_per_module_tests.md) |
| How do I sync to HPC? | [`scripts/sync_to_hpc.sh`](scripts/sync_to_hpc.sh) |
| How are containers built? | [`containers/README.md`](containers/README.md) |
| What are the known issues? | [`notes/inventory/05_findings.md`](notes/inventory/05_findings.md) |

## Development workflow

Develop on Mac with Claude Code, sync to HPC with
`./scripts/sync_to_hpc.sh`, run the pipeline or tests there, iterate.
The HPC has no git; the Mac repo is authoritative and the only place
commits happen. Round-trips are slow — plan thoroughly before
transferring, batch related changes, and lean on `pytest -m local_unit`
plus static analysis (`ruff`, `vulture`, `radon`) for fast feedback
during development.

## License / contact

TODO — license not yet specified.

Maintainer: Joshua Williams (joshuandwilliams@outlook.com).
