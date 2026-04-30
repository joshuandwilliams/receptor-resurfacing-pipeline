# Cluster submission guide

Short reference for submitting commands on the HPC cluster for this project.

## Environment

- **HPC cluster** with SLURM.  All heavy work dispatched as sbatch jobs.
- **Singularity container** with Boltz2, Boltz1, Chai1, ColabFold at
  `/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/benchmark_models.img`
- **Working directory**: `/hpc-home/jowillia/receptor_design/negative_steering/`
- **Runs directory**: `runs/<experiment_name>/` under the working dir.
- **Pipeline scripts** are at the root of the working dir.  User edits
  them locally, copies (scps or rsyncs) to the cluster, and runs from
  there.

## Submit scripts and what they do

Five submit scripts, each wrapping a different phase of work:

### `submit_boltz2_negative_steering.sh`
Main steering pipeline.  Generates K unique substituted-receptor designs,
predicts each × `num_seeds` with Boltz, flags contaminated designs, runs
reversion predictions, classifies verdicts.  Produces
`<workdir>/steered_results.csv` and `<workdir>/summary.txt`.

**Typical usage:**
```bash
./submit_boltz2_negative_steering.sh \
    --ground-truth resurface_pipeline_test/design_3.pdb \
    --receptor A --effector B \
    --receptor-fasta resurface_pipeline_test/design_3_seq_3_receptor.fasta \
    --workdir runs/design_3_dedup_test_vN \
    --mode mild --max-mutations 4 --candidate-pool-size 8 \
    --n-designs 50 --num-seeds 3 --n-cycles 1 \
    --diffusion-samples 5 \
    --protected-set-source design_region_union \
    --true-interface-indices-file resurface_pipeline_test/design_3_true_interface.txt \
    --design-region-indices-file resurface_pipeline_test/design_3_design_region.txt
```

Submits the whole chain (plan → steering array → collect → kickoff
distances → reversion → finalize) as linked sbatch jobs.  Returns
immediately after queueing.

### `submit_postprocess.sh`
Post-processing chain to run AFTER the main queue has drained.  Four
phases: aggregate → compute-final-metrics → aggregate-per-sequence →
extract_passing.  Single SLURM CPU job, ~3 minutes.

**Typical usage:**
```bash
./submit_postprocess.sh --experiment-root runs/design_3_dedup_test_vN
```

Optional flags:
- `--force-recompute` — re-run compute-final-metrics on every row,
  ignoring any cached output.  Use when the postprocess code has
  changed and cached rows are stale.  (Default is `--skip-existing`.)
- `--rmsd-threshold`, `--metric-column`, `--contact-cutoff`,
  `--populate-all` / `--no-populate-all`.

Produces:
- `<experiment-root>/raw_per_seed_results.csv`
- `<experiment-root>/aggregated_results.csv`
- `<experiment-root>/passing_summary.csv`

### `submit_multiseed_variance_test_slurm.sh`
Standalone diagnostic: run one receptor-effector sequence through
Boltz N times with different seeds.  Three stages: plan (GPU),
predict-one (GPU array), harvest (CPU in container).

**Typical usage:**
```bash
./submit_multiseed_variance_test_slurm.sh \
    --ground-truth  resurface_pipeline_test/design_3.pdb \
    --receptor A --effector B \
    --receptor-fasta resurface_pipeline_test/design_3_seq_3_receptor.fasta \
    --workdir runs/variance_test_design_3_baseline \
    --n-seeds 100
```

Produces `<workdir>/per_seed_metrics.csv` and
`<workdir>/variance_summary.csv`.

### `submit_boltz2_iterate_steering.sh`
Internal — used by cycle-N+1 iteration chains.  Not usually called
directly by the user.

### `submit_final_validation.sh`
Post-pipeline tool to re-validate survivors from a completed run
with multiple seeds.  Currently unused since the main pipeline runs
with `num_seeds >= 3` by default.

## Typical workflow

1. **Edit pipeline scripts locally** (in this chat).
2. **Copy changed files to the cluster** — scp / rsync
   the edited files into the working dir.
3. **Submit the experiment** — one of the submit scripts above.
4. **Wait for the queue to drain** — check with `squeue -u $USER`.
   Wait until no `b2ns_*` or `b2is_*` jobs remain.
5. **Run postprocess** — `./submit_postprocess.sh --experiment-root ...`.
6. **Download the output CSVs** locally for inspection / upload to chat.

## File deployment reminder

When the pipeline code (`boltz2_iterate_steering.py`, `reversion.py`,
`extract_postprocess.sh`, etc.) has changed, BOTH the Python file AND
any affected shell scripts need to ship together.  A previous
deployment confusion happened where `submit_postprocess.sh` was updated
on the cluster but the Python file wasn't — the old classifier kept
running.  Diagnostic: the `aggregated_verdict_reason` string in
`aggregated_results.csv` changes when the classifier changes; if the
new phrasing doesn't appear, the new code didn't run.

MD5-checksum deployed files against local stages if uncertain.

## Queue monitoring

```bash
squeue -u $USER                              # all my jobs
squeue -u $USER -o '%.10i %.20j %.10T %R'    # with status and reason
sacct -j <JOBID>                             # post-mortem on a completed job
tail -f runs/<exp>/logs/<jobname>_<jobid>.out   # follow a running job
```

## Sanity check a completed experiment

After postprocess:

```bash
# How many rows at each stage?
wc -l runs/<exp>/raw_per_seed_results.csv          # N_designs × num_seeds + 1
wc -l runs/<exp>/aggregated_results.csv            # unique_sequences + 1 singleton
wc -l runs/<exp>/passing_summary.csv               # passing rows

# Verdict tally on aggregate
python3 -c "
import csv, collections
rows = list(csv.DictReader(open('runs/<exp>/aggregated_results.csv')))
print(collections.Counter(r.get('aggregated_verdict','') for r in rows))
"
```

## Non-cluster testing

Any pure-Python phase (`aggregate`, `aggregate-per-sequence`,
`extract_passing.py`) can be run directly with `python3` outside the
container — no GPU or singularity required.  Only the Boltz
prediction steps and `compute-final-metrics` need the container.
