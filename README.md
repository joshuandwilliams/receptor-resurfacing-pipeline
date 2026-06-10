# Receptor Resurfacing Pipeline

A Nextflow DSL2 pipeline for *de novo* binder design: takes a docked
receptor–effector complex, redesigns the receptor surface around the
effector pose, and outputs candidate binder sequences with predicted
quality metrics from multiple orthogonal validators.

## Status

This codebase is undergoing structured remediation following
[`notes/codebase_remediation_plan.md`](notes/codebase_remediation_plan.md).
Active work happens on the `phase4-impl` branch; `main` is the
safe-fallback baseline (Phase 4 spec, no implementation code). Start
every session by reading
[`notes/SESSION_HANDOFF.md`](notes/SESSION_HANDOFF.md) for orientation;
[`notes/remediation_state.md`](notes/remediation_state.md) is the
authoritative session log and branch state.

## What the pipeline does

The input is PDB only. Branch B takes a single pre-docked complex PDB
(receptor + effector already positioned). Branch A takes two separate
monomer PDBs and docks them with the constraint-driven pose solver
before redesign.

### Branch A — pose solver input contract

When you supply `params.receptor_input` + `params.effector_input`,
each must be a **monomer PDB** (not a pre-aligned complex split into
chains). The constraint-driven pose solver
([`modules/pose_solver.nf`](modules/pose_solver.nf)) places the two
monomers into a compact, clash-free arrangement that satisfies
user-supplied CA–CA pair restraints, giving RFDiffusion a good
starting point. Its job is **geometric placement** only — energetics
are redesigned downstream, so the pose needs to be geometrically
sensible, not energetically optimal.

Branch A **requires** `params.pose_solver_pairs`: one or more
hard-pinned CA–CA contacts, written as space-separated,
chain-prefixed residue pairs — e.g. `"A73-B31 A71-B33 A8-B22"`, with
an optional `@distance` target per pair (e.g. `"A73-B48@4.5"`). `A`
is the receptor, `B` the effector. See the "Pose Solver" block in
[`params_example.yml`](params_example.yml) for the rest of the knobs
(exclusions, min/max pair distance, clash and contact cutoffs,
restart count, the differential-evolution toggle, and interpolation
weight).

An optional manual-pick checkpoint lets you inspect the solved pose
before the expensive design stages: run with
`params.stop_after_pose_solve: true`, open
`${outdir}/pose_solver/solved_pose_posed.pdb` in ChimeraX, then resume
with `stop_after_pose_solve: false` (`-resume` reuses the solved pose).

From this complex, RFDiffusion generates candidate receptor backbones
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
tests/       Per-module tests, characterization tests, reference data,
             and the curated fixtures each per-module test runs against
containers/  Singularity .def recipes; image + JDK symlinks (git-ignored) go here
notes/       Inventory documents, remediation plan, decisions, glossaries
scripts/     Repo tooling (sync_to_hpc.sh, check_environment.slurm.sh, …)
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

Three layers, each with a different cost and a different question it answers:

- **`pytest -m local_unit`** (Mac, fast). Tests the comparator framework itself — the helpers in `tests/characterization/helpers/` that pin reference outputs against fresh runs. Run on every commit during development.
- **Per-module Nextflow tests** (HPC, minutes per stage). Each `tests/<module>/test_<module>.nf` runs one pipeline stage end-to-end against a curated fixture in `tests/<module>/data/`. The fixtures are hand-picked from a discovery pipeline run to exercise every distinct path through that stage. Run with the per-module wrappers, e.g. `sbatch tests/negative_steering/run_test_negative_steering.slurm.sh`. This is the primary characterization safety net during refactoring.
- **`pytest -m hpc`** (HPC, slow). Cohort-level characterization tests that pin output files from a fresh pipeline run against committed reference outputs. Currently parametrized over a legacy sequence trio; restructuring against the per-module reference sets is in progress.

```bash
# Local
pytest -m local_unit

# Per-module on HPC
sbatch tests/<module>/run_test_<module>.slurm.sh

# Cohort-level characterization on HPC
RECEPTOR_OUTPUT_ROOT=<run-output-dir> pytest -m hpc
```

See [`tests/characterization/README.md`](tests/characterization/README.md)
for the comparator framework, marker tiers, and how to update reference
outputs after intentional behaviour changes.
[`notes/inventory/14_phase_2_revision_per_module_tests.md`](notes/inventory/14_phase_2_revision_per_module_tests.md)
covers the per-module testing strategy in full.

## Where to look for what

| Question | Document |
|---|---|
| How does the pipeline work? | [`notes/inventory/04_functional_categorization.md`](notes/inventory/04_functional_categorization.md) |
| What does *\<term\>* mean? | [`notes/inventory/06_ubiquitous_language.md`](notes/inventory/06_ubiquitous_language.md) |
| Where is the remediation at? | [`notes/remediation_state.md`](notes/remediation_state.md) |
| What's the testing strategy? | [`notes/inventory/14_phase_2_revision_per_module_tests.md`](notes/inventory/14_phase_2_revision_per_module_tests.md) |
| Why these test fixtures? | [`notes/inventory/15_discovery_run_path_coverage.md`](notes/inventory/15_discovery_run_path_coverage.md) |
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