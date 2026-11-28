# scripts/

Operational helpers for running the pipeline on a cluster and for maintaining
its test fixtures. Nothing here is invoked by `main.nf`. Everything in `bin/`
is. Everything here is run by hand.

Two scripts that used to live here have moved to
[`resurfacing-experiments`](https://github.com/joshuandwilliams/resurfacing-experiments).
Both are downstream of the pipeline rather than part of it.
`make_gblocks.py` reverse-translates chosen sequences for gBlock ordering.
`RFDiffusion_JoshW_v4_incremental_save_2.ipynb` is the Colab notebook this
pipeline replaced.

## Site configuration

| | |
|---|---|
| `site_conf.sh` | Sourced by the others. Loads `site.conf` without clobbering anything already in the environment. |

Precedence, highest first: a value already exported (including one `sbatch`
carried from the submitting shell), then `site.conf`, then each consumer's own
default. That ordering is what lets a one-off run override a setting on the
command line without editing `site.conf`.

Copy `site.conf.example` to `site.conf` and fill in the sixteen `RECEPTOR_*`
values it lists:

- three SLURM partitions, plus the account and the QOS
- the GPU request string
- four AlphaFold database and model paths
- the negative-steering install directory and container image
- the output root
- three values the sync scripts need

`site.conf` is gitignored. A checkout therefore carries no cluster's queue
names, account or filesystem layout.

No script in this repository hardcodes a partition.
`check_environment.slurm.sh` needs a GPU one, so submit it as
`sbatch -p "$RECEPTOR_PARTITION_GPU" scripts/check_environment.slurm.sh`. A
`sbatch -p` on the command line beats both the environment and a header
line.

## Cluster operations

| | |
|---|---|
| `check_environment.slurm.sh` | Pre-flight check on a GPU node. Confirms every container, model weight and CUDA runtime the pipeline needs is present and loadable, before a 72-hour run discovers it is not. |
| `clean_finished_run.sh` | Reclaims disk and inodes from one finished run. Refuses a run that was never executed, warns loudly on a failed one, then deletes `work/`, `.nextflow/cache` and `.nextflow/plr`. Those three are orphaned once `work/` is gone. |

`clean_finished_run.sh` takes a run directory, which lives in
`resurfacing-experiments` under `experiments/campaigns/<campaign>/runs/<run>/`.

## Fixture maintenance

The characterization suite reads committed example outputs under
`tests/<module>/example_output_files/`. These three regenerate them. All three
touch the cluster, so read what they do before running one.

| | |
|---|---|
| `sync_to_hpc.sh` | rsync the checkout to the cluster over the `slurm` SSH alias. `--dry` first. |
| `sync_from_hpc.sh` | The reverse. Pulls regenerated fixtures back so they can be committed. `--module <name>` narrows it to one. |
| `refresh_orthog_fixture.sh` | Rebuilds the `orthogonal_metrics` input fixture from the newest negative-steering test outputs, so that test runs against a full cohort rather than the staged two-survivor mini-fixture. |
| `rewrite_fixture_paths.py` | A fixture staged from the cluster carries absolute cluster paths inside `plan.json` and `cross_sequence_summary.csv`. This rewrites them to the local checkout, which is what makes the staged tree usable off the cluster. |

Order for a full refresh: `sync_to_hpc.sh`, then
`tests/update_example_dataset.slurm.sh` on the cluster, then
`sync_from_hpc.sh`, then `rewrite_fixture_paths.py`.

## Teaching

| | |
|---|---|
| `pipeline_intro_for_students.ipynb` | Walks through what each stage does and what its output looks like, on the committed example outputs. Needs no cluster. |
