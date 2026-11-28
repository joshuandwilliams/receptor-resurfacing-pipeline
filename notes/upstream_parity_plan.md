# Upstream parity plan

What to implement here so this repo reaches the state
`structure-negative-steering` reached across its 43 merged pull requests. That
is #1–#33 and #35–#44, since #34 does not exist, plus the split that produced
`negative-steering-analyses`.

Written 2026-11-27. Every claim below was checked against the file it names.
Where a number appears, the command that produced it is given.

---

## Decisions taken, 2026-11-27

**The engine does not live here.** `structure-negative-steering` is installed
per site and called as a tool. The installation is at
`/jic/scratch/groups/Mark-Banfield/JoshW_Software/negative_steering`, locally
visible as `/Volumes/Group-Scratch/Mark-Banfield/...`. It holds a checkout at
`4f70608`, the merge of PR #38, with a built `image/negsteer.img` and its own
`run_NegativeSteering.slurm.sh`.

That supersedes §4 of this plan. Every one of the six defects §4 lists was
fixed upstream in PR #16, #18, #22, #28, #30 or #31, all of which are below
#38, so the installed engine already has them. There is nothing to port, so
the vendored copies get deleted instead.

The engine exposes a versioned tool interface: `negsteer run`, writing
`negsteer_result.json` at output contract version 1, with five named outputs
(`negsteer/__init__.py`). `--effector-template` is already opt-in and off by
default there, which is the §7.3 question answered. Upstream also left a
62-line replacement Nextflow module in its history, deleted in `5807225` with
the stated reason that it belongs in this repository. Recover it with
`git show 5807225^:contrib/negative_steering.nf`.

The installed checkout is **12 commits behind** `origin/main`, which is PRs #39
to #44. Those are the Docker build fixes, keeping `mols.tar` in the image, the
container resolver and the single-run assertions. Updating it is a change on
the cluster and has not been made.

**The shared library modules decouple fully.** `pipeline_thresholds.py` is
rewritten for this repository's own stages, `contig_utils.py` stays as
RFDiffusion support, and `compute_metrics.py` folds into
`compute_interface_metrics.py`. No engine code and no engine imports. The
coupling to undo is real. `pipeline_thresholds` is imported by 11 files here,
`contig_utils` by 6 and `compute_metrics` by 5.

**The analyses split into their own repository**, with this one as a submodule.
This has not been executed, for a reason found after the decision was taken.
See the blocker below.

## The blocker on the analyses split

The twelve numbered analyses read their inputs from
`experiments/campaigns/pikp1_avrpikf/runs/crystal_full_test_contig/results/`
and `.../agroinfiltrations/`. Both are matched by `.gitignore:101` and
`.gitignore:153` and carry **zero tracked files**, and they exist on local
disk only.

So `git mv` offers them no protection, a new repository would not receive
them, and nothing in version control would notice if they were lost. They are
the primary sources behind every figure and number in Chapter 5, and the
analyses resolve them by relative path, so moving the analyses without them
breaks every extractor.

What is safe, and what is not:

- The thesis renders **without touching this repository**. Chapter 5's
  figures resolve to `chapters/05-.../figures/thesis-figures/`, 28 files that
  physically exist there, and the only mention of this repo in that chapter is
  a GitHub URL in prose at line 188. Moving tracked files here cannot break a
  render.
- 16 of the 118 tracked images under the campaign's analyses are byte-identical
  to a file in the thesis figure tree. Those are the preserved provenance.
- Four are not, so the pipeline copy has already drifted from what the thesis
  shows: `03-mpnn-sequence-design/thesis-figures/loop_depth_composition.png`,
  `07-design-region-msa/thesis-figures/design_region_msa.png`,
  `06-.../thesis-figures/constraint_level_and_rank.png` and
  `09-.../thesis-figures/wt_vs_design5_panel.png`. For the first of those the
  thesis-matching bytes do survive, under a different name, as
  `03-mpnn-sequence-design/figure-options/composition_stack_fixed.png`.

## Status, end of 2026-11-27

Done and committed: §1.3, §2, §3, §5, §6.1, §6.2, §6.4, and §8.2's guards.
Not done: the engine decouple and rewiring, the analyses split, §6.3, §6.5,
§8.1, §8.3 to §8.6, §9.

---

## 0. What is already done, so nobody redoes it

Two upstream fixes are already present here. Both were ported by hand at some
earlier point, independently of the rest.

- **The reversion majority gate.** `bin/boltz2_iterate_steering.py:1562-1563`
  computes `majority_threshold = math.ceil(n_correctly_placed / 2)` and skips a
  design below it. That is upstream PR17's rule. Upstream's copy is at
  `bin/reversion.py:1430-1444` and uses the same formula.
- **The plot-module logging and placeholder writes.** Neither
  `bin/negsteer_plots.py` nor `bin/negsteer_within_sequence_plots.py` contains
  an `or True` lambda wrapper (`grep -n "or True"` returns nothing in both).
  Every no-data branch calls `_make_empty_plot`, which calls `savefig`, then
  returns False. That is upstream PR18 and PR21.

One upstream defect class cannot occur here. Upstream PRs 40 and 43 fixed shell
scripts that scraped `boltz_container:` out of a YAML with grep and sed and got
back the literal string `null`. No `.sh` or `.nf` file in this repo extracts
that key (`grep -rn boltz_container --include='*.sh' --include='*.nf'` finds no
extraction). Containers here come from Nextflow params, set at
`nextflow.config:36-40`. The defect is absent because the architecture differs,
not because it was fixed.

---

## 1. Prerequisites, before any of the rest

**1.1 Deal with the dirty working tree.** `git status --porcelain` reports 29
modified tracked files and 26 untracked paths. The modified set includes
`main.nf` (+402 lines, adding the `NEGSTEER_ONLY` and `MAIN_PIPELINE` entry
workflows and `params.negsteer_rehydrate_from`), nine `bin/` modules, and
`containers/pytest_runner.def`. None of the work below should start on top of
uncommitted changes of that size, because a failed test then has two candidate
causes.

**1.2 Decide what happens to `main`.** `phase4-impl` is 71 commits ahead of
`main` and 0 behind, and `git merge-base --is-ancestor main phase4-impl`
succeeds, so it fast-forwards cleanly. Upstream did all 43 PRs against its
`main`. Doing the same here means merging `phase4-impl` to `main` first, or
accepting that this repo's default branch stays behind.

**1.3 Install the declared test dependency.** `scikit-image` is declared at
`pyproject.toml:14` and is not installed in the active interpreter
(`/Users/jowillia/miniforge3`, Python 3.12). Two collection targets fail on
`ModuleNotFoundError: No module named 'skimage'`:
`tests/characterization/helpers/tests/test_png_compare.py` and
`tests/characterization/test_plots.py`.

---

## 2. The test suite is red, and the cause is a pandas 3 change

This is first because nothing else can be measured until the suite is green.

**2.1 Fix the path-normalisation comparator under pandas 3.**
`python -m pytest -q -m local_unit` (ignoring the two skimage files) reports
**2 failed, 410 passed, 7 skipped, 272 deselected**. Both failures are in
`tests/characterization/helpers/tests/test_csv_compare.py`. The failing cases
are `test_compare_csv_exact_modulo_paths_normalises_workdir_hash` and
`test_compare_csv_exact_modulo_paths_normalises_hpc_home_prefix`.

The mechanism: `tests/characterization/helpers/csv_compare.py:125` guards the
normalisation with `if df[col].dtype == object`. The installed pandas is 3.0.1,
where a text column read by `read_csv` has dtype `str`, not `object`. I
confirmed it directly, and `pd.read_csv(...)['home_path'].dtype` prints `str`
while `dtype == object` is False. So the normaliser at
`tests/characterization/helpers/path_normalize.py:18` is never applied, and two
paths that should collapse to `<HOME>` and `<WORKDIR>` are compared raw.

The normaliser itself is correct. The regex at `path_normalize.py:13-15` does
match `/hpc-home/<user>`, so only the dtype guard is incorrect.

The fix is to select string columns by `not is_numeric_dtype(...)` rather than
by identity against `object`. `csv_compare.py:153` already uses that predicate
for a different purpose in `_string_column_set`.
`pyproject.toml:14` declares `pandas>=2.0` with no upper bound, so this also
needs a decision on whether to pin below 3 or support both. Supporting both is
one predicate.

**2.2 Stop `pytest tests/` crashing.** `python -m pytest --collect-only -q
tests/` ends in `INTERNALERROR`, not a collection error. The cause is
`tests/negative_steering/test_negsteer_plots.py:37`, which raises `SystemExit`
at module import. It is a command-line harness rather than a test module, and
pytest's assertion-rewriting importer does not handle `SystemExit` at import.
Five files share this shape and contain zero `def test_` functions between
them:

| file | `def test_` count |
|:---|---:|
| `tests/negative_steering/test_negsteer_plots.py` | 0 |
| `tests/negative_steering/test_negsteer_within_sequence_plots.py` | 0 |
| `tests/orthogonal_metrics/test_orthogonal_metrics_plots.py` | 0 |
| `tests/proteinmpnn/test_mpnn_plots.py` | 0 |
| `tests/rfdiffusion/test_rfdiffusion_plots.py` | 0 |
| `tests/rosetta_filtering/test_rosetta_filter_plots.py` | 0 |

Renaming them off the `test_*.py` pattern is the fix, because they are
harnesses and the name is what makes pytest import them.

**2.3 Rename the markers to what they gate on, as upstream PR1 did.**
`pyproject.toml:18-27` declares eight markers. Five carry no test:

| marker | uses |
|:---|---:|
| `local_unit` | 184 |
| `hpc` | 48 |
| `wave1` | 48 |
| `local_integration` | 0 |
| `wave2` | 0 |
| `requires_populated_reversion` | 0 |
| `requires_controls` | 0 |
| `slow` | 0 |

Upstream collapsed its two laptop tiers into one marker `logic` and renamed
`hpc` to `check_test_run`, because `hpc` names a machine while the real
requirement is a finished run. The same rename applies here. `wave1` and
`wave2` are remediation-phase labels rather than requirements, so they go with
it.

**2.4 Widen `testpaths`.** `pyproject.toml:23` sets
`testpaths = ["tests/characterization"]`, so a bare `pytest` run collects 691
tests and never looks at the other seven directories. That is only safe to
widen after 2.2.

---

## 3. Coverage, lint and CI, none of which exist here

**3.1 Add coverage with a ratcheted gate.** There is no coverage configuration
anywhere. `grep -ni "cov" pyproject.toml` returns nothing, there is no
`.coveragerc` and no `codecov.yml`, and neither `coverage` nor `pytest-cov` is
in the `test` extra at `pyproject.toml:10-16`. Upstream runs coverage over
`bin/` and `negsteer/` on every invocation through `addopts` and gates it with
`fail_under`, raising the number as tests land. Its final value is 100.0.

Starting here means measuring first and setting the gate to whatever the first
honest run reports. Upstream's own sequence was 69.5 → 100 over PRs 1, 7,
13–21, so the gate is a ratchet rather than a target to hit at once.

**3.2 Add `.github/workflows/test.yml` and `lint.yml`.** There is no `.github`
directory. Upstream's test workflow installs the `test` extra, runs
`pytest -m logic --cov-report=xml`, and uploads to Codecov only when
`CODECOV_TOKEN` is set. Its second job stub-runs the whole DAG on Nextflow
pinned to 25.10.4 with `-stub-run -profile local`, and asserts the cluster
profile keeps its queue, `clusterOptions` and `maxForks` through profile
resolution.

The Nextflow pin matters and this repo is where it is defined.
`containers/NextFlow.def:23-24` downloads `nextflow-25.10.4-dist`. Upstream's
workflow comment says the pin comes from that file in this repo, which is
accurate.

**3.3 Adopt ruff with upstream's rule set.** No lint configuration exists
here: no `[tool.ruff]`, no `.flake8`, no `.pre-commit-config.yaml`. Under
upstream's exact rules (`--select E,F,W,I --ignore E501 --line-length 100`)
this repo reports **498 errors**, 225 of them auto-fixable:

| rule | count | what it is |
|:---|---:|:---|
| I001 | 125 | import block unsorted |
| E701 | 112 | multiple statements on one line (colon) |
| E702 | 72 | multiple statements on one line (semicolon) |
| E402 | 44 | module import not at top of file |
| F541 | 39 | f-string with no placeholder |
| W292 | 26 | no newline at end of file |
| F401 | 19 | imported and unused |
| E741 | 16 | ambiguous variable name |
| F841 | 15 | assigned and never used |
| F811 | 12 | redefinition of unused name |
| F821 | 6 | undefined name |
| other | 12 | E401, E731, W291, W293 |

By directory: `scripts/` 183, `bin/` 153, `experiments/` 132, `tests/` 28.

Two of these are worth separating from the mechanical bulk. The six F821
undefined names are all in analysis helpers and notebooks, not in `bin/`.
Three are `shap` in
`experiments/campaigns/pikp1_avrpikf/analyses/05-steerability-correlates/helpers.py:297,305,307`
and three are `design_residues_pdb` in
`scripts/RFDiffusion_JoshW_v4_incremental_save_2.ipynb`. The F811
redefinitions are in notebooks and in
`analyses/_archive/design_region_conservation/make_figure_options.py`. Upstream
PR7 found exactly this class in its own tests, where a redefined test class
silently discarded three earlier tests.

**3.4 Fix the one real syntax incompatibility.**
`experiments/campaigns/pikp1_avrpikf/analyses/08-agroinfiltration/data/map_spots.py:55`
uses a backslash escape inside an f-string, which is Python 3.12 syntax.
`pyproject.toml:8` declares `requires-python = ">=3.10"`, so this file cannot
be parsed on the floor of the declared range. Ruff reports it as
`invalid-syntax` under both 3.10 targets.

---

## 4. The six engine defects to port

This repo vendors an older copy of part of the upstream engine. Six basenames
exist in both `bin/` directories, and all six have diverged:

| file | upstream | here | differing lines |
|:---|---:|---:|---:|
| `boltz2_negative_steering.py` | 3348 | 3596 | 1664 |
| `compute_metrics.py` | 1510 | 1695 | 826 |
| `extract_passing.py` | 548 | 694 | 695 |
| `reversion.py` | 1849 | 1307 | 1360 |
| `pipeline_thresholds.py` | 164 | 168 | 154 |
| `negative_steering_run_one.sh` | 488 | 458 | 70 |

The divergence is two-way, so none of this is a file copy. It is six targeted
patches.

**One structural fact governs where each patch lands.** Upstream split its
monolithic `bin/boltz2_iterate_steering.py` into `negsteer_coldstart.py`,
`negsteer_aggregate.py`, `pathway_labels.py` and the CLI half of
`reversion.py`. This repo never took that split. The monolith is still here at
`bin/boltz2_iterate_steering.py`. This repo's `bin/reversion.py` is the
pre-split library half, with no `main`, no `build_parser` and no `cmd_*`
functions. So an upstream patch to `reversion.py` or `negsteer_aggregate.py`
has to be re-aimed at the matching function inside the monolith.

**4.1 Emit `templates:` at the YAML root (upstream PR22).**
`bin/boltz2_negative_steering.py:750-761` appends the templates block at
six-space indent, inside the effector's `- protein:` entry. The comment on
lines 751-753 asserts that placement is correct. It is not. Boltz 2.2.1 reads
`templates` only from the schema root. A chain entry is inspected only for the
key naming its entity type, so the nested block is dropped with no warning.
Upstream verified this from Boltz's own run records, where both committed
fixtures carry `"templates": []`.

This is the same defect `notes/things_to_fix.md` already records under
"Effector template is silently ignored". Fixing it changes predictions, so it
is listed again in §7 as a decision rather than a free fix.

**4.2 Resolve thresholds once, not per call (upstream PR30).**
`bin/extract_passing.py:202-205` does `import sys` then
`sys.path.insert(0, str(Path(__file__).resolve().parent))` inside
`_compute_confidence_flag`, with no `not in sys.path` guard. The function is
called once per row from `extract_row`. Every call appends another duplicate
entry to `sys.path` and re-resolves `__file__` through `realpath`. Upstream
measured the per-call cost at 66.89 µs against 0.81 µs once resolved at import.

**4.3 Delete the dead majority threshold (upstream PR18).**
`bin/extract_passing.py:297-300` computes `threshold` from a four-line
conditional on `threshold_majority_frac`, and line 303 overwrites it
unconditionally with `(n // 2) + 1`.

I checked whether the two agree. At the default `frac=0.5` they agree for
n=1..5, so there is no live consequence today. At any other fraction they
diverge. At `frac=0.7, n=3` the conditional gives 3 and the overwrite gives 2.
At `frac=0.34, n=2` the conditional gives 1 and the overwrite gives 2. So the
parameter is dead rather than wrong. All six occurrences of
`threshold_majority_frac` in the repo are inside `extract_passing.py` itself,
so nothing passes a non-default value. The fix is deletion, which is what
upstream did.

**4.4 Index the Boltz sidecars once (upstream PR28).**
`bin/compute_metrics.py` runs three chained `glob.glob` calls per sidecar kind
per PDB, inside the `for pdb_path in pdb_files:` loop. The call sites are
confidence at 1270-1272, pLDDT at 1296-1298 and PAE at 1325-1327. The third
call in each triple is `recursive=True` over `pred_dir`, which is
loop-invariant. At the benchmark's 25 models that is 75 recursive walks of one
directory. This repo has no `_sidecar_index` and no `_sidecar` function
(`grep -c` returns 0). Upstream measured 4.2× on `parse_boltz2` end to end.

The preference order has to be preserved exactly, because the caller takes the
first hit: beside the PDB, then the sibling category directory, then shallowest
anywhere.

**4.5 Remove the `id(rec)` sort tie-break (upstream PR31).**
`bin/cross_sequence_summary.py:617` and `:627` end their rank keys with
`id(rec)`. `list.sort` is already stable and preserves input order, so a third
key of a memory address replaces that with something that varies between
processes. Upstream documented this rather than fixing it, because deleting it
changes committed output. Here the stakes are higher than upstream's. Its
benchmark has one sequence per target, while a campaign cohort here is 128
sequences, so ties are reachable. The ranks are columns in
`cross_sequence_summary.csv`.

**4.6 Write a header-only CSV on zero rows (upstream PR16).**
In `bin/boltz2_iterate_steering.py` the guard `if rows:` opens at line 4335 at
four-space indent, and its body writes the CSV through line 4381. The
`print(f"Wrote {out_csv} ({len(rows)} rows)")` at line 4382 is back at
four-space indent. So a zero-row aggregate prints that it wrote a file it never
created, and the next reader gets `FileNotFoundError` having been told
otherwise. Upstream's fix also needed a companion change to
`_aggregate_csv_fieldnames`, which has no fallback when there are no rows to
intersect against.

**4.7 Stop the runner script aborting on its own closing banner.**
`bin/negative_steering_run_one.sh:37` sets `-uo pipefail` without `-e`. Lines
286 and 360 then do `set +e`, and lines 291 and 365 do `set -e`. That turns
errexit **on** for the remainder of the script rather than restoring the
previous state. Line 291 sits inside the per-design predict loop, which runs
whenever there is at least one design, so errexit is normally on by the end.

The last functional line is a cosmetic
`ls -1 "$WORKDIR_ABS"/*.csv 2>/dev/null | sed 's/^/    /'` with no `|| true`.
With errexit on and `pipefail` set, a workdir holding no CSV makes that
pipeline return 1 and aborts the script before its closing banner and its
`exit 0`. I reproduced the mechanism in isolation. The demo script printed
`reached tail`, never printed `AFTER ls`, and exited 1.

The exposure is a run that completes with every design failing. Upstream's
fixed copy appends `|| true` to that line and keeps errexit off throughout.

**4.8 Audit, not patch.** `bin/pipeline_thresholds.py` differs on 154 lines and
every one is prose. No numeric threshold or field name differs between the two
copies. It needs no change beyond whatever §6 does to comments.

---

## 5. Portability: one gitignored site file

Nothing here is parameterised by site. Upstream PR23 moved every
site-specific setting into one gitignored `site.conf`, with a tracked
`site.conf.example`, a `scripts/site_conf.sh` loader and a `scripts/submit.sh`.
This repo has none of those four files, and no `.env`.

**5.1 Parameterise the partitions.** `nextflow.config` has 30 `queue =` lines
across 29 `withName:` blocks, 27 occurrences of `jic-medium` and 5 of
`jic-gpu`. Upstream's approach reads `NEGSTEER_PARTITION_GPU`,
`NEGSTEER_PARTITION_CPU` and `NEGSTEER_GPU_CLUSTER_OPTIONS` from the
environment, and its CI asserts both the set case and the bare case. A null
queue is safe. Nextflow 25.10.4's `SlurmExecutor.groovy:98` guards the
directive with `if( task.config.queue )`, so null omits `-p` rather than
submitting `-p null`.

**5.2 Remove the hardcoded paths, email and ticket id from infrastructure.**
The counts split cleanly between infrastructure, which must be parameterised,
and run provenance, which is a record of what executed and should stay:

| pattern | tracked files total | infrastructure | data and provenance |
|:---|---:|---:|---:|
| `jowillia` | 719 | 4 | 715 |
| `/hpc-home` | 682 | 9 | 673 |
| `jic-medium` | 24 | 4 | 20 |
| `@nbi.ac.uk` | 20 | 2 | 18 |
| `jic-gpu` | 14 | 3 | 11 |
| AF3 package UUID | 7 | 2 | 5 |
| `/nbi/` | 5 | 3 | 2 |
| `jic-short` | 2 | 0 | 2 |

The infrastructure sites in `run_pipeline.slurm.sh` are `:12` (`--mail-user`),
`:38` (`PIPELINE_DIR`) and `:46` (`AF2_DATA_DIR`). `nextflow.config` adds
`:46-47` (the two AlphaFold reference-data paths) and `:44` (the AF3 package
UUID). The remaining sites are `params_example.yml:5`,
`scripts/clean_all_finished_runs.slurm.sh:9,14`,
`scripts/check_environment.slurm.sh:57-58`,
`tests/characterization/run_pytest.slurm.sh:3,12,13,42` and
`tests/update_example_dataset.slurm.sh:3,11`.

Upstream deleted the mail address rather than moving it, because `sbatch(1)`
says `--mail-user` defaults to the submitting user, so `--mail-type` alone is
correct everywhere.

**5.3 A third site-path form exists here that upstream never had.**
`experiments/campaigns/pikp1_avrpikf/analyses/06-confidence-and-biophysical-metrics/sequence_match_extract.py:36-38`
hardcodes `/Volumes/HPC-Home/receptor_design/...`, the macOS mount of the
cluster home. Any site handling has to cover it alongside `/hpc-home` and
`/nbi/`.

---

## 6. Repository setup and documentation

**6.1 Licence and citation.** There is no `LICENSE` and no `CITATION.cff`.
Upstream PR5 added MIT and a CFF file. `notes/things_to_fix.md` already notes
that MIT is a precondition for nf-core membership.

**6.2 Collapse the four version strings.** They disagree:

| file:line | value |
|:---|:---|
| `pyproject.toml:7` | `0.0.0` |
| `nextflow.config:3` | `v0.2.0` (header comment) |
| `nextflow.config:371` | `0.2.0` (`manifest.version`) |
| `params_example.yml:2` | `v0.3.0` |
| `main.nf:5` | `v0.4.0` |

`notes/things_to_fix.md` records three of these, but there are four distinct
values across five sites. Upstream settled on one development version,
`0.0.0.9000`, written once.

**6.3 Logo and favicon.** Upstream PRs 5 and 9 generated a hex sticker from a
committed script, `docs/make_logo.py`, with the frame measured off CDALabeller
rather than eyeballed. There is no `docs/` directory here at all.

**6.4 Cut the README and add the badges.** `README.md` is 154 lines and 885
words. Upstream cut its own from 476 lines to 219, moved the detail that
records a decision into `docs/implementation-notes.md`, and replaced
hand-written test-count and coverage badges with Codecov. This README's Status
section at `:7-17` points at `notes/codebase_remediation_plan.md`,
`notes/SESSION_HANDOFF.md` and `notes/remediation_state.md` as the entry point.
`SESSION_HANDOFF.md:1` is dated 2026-05-14 and describes an in-flight full
pipeline test, so the README currently directs a reader to a six-month-old
session brief as orientation.

**6.5 Separate the living notes from the durable ones.** `notes/` holds 72
tracked files. Roughly 47 are per-session state: `SESSION_HANDOFF.md`,
`remediation_state.md`, `pipeline_notes/pipeline_notes1-16.md`,
`negsteer_notes/notes1-11.md`, nine `todo_list*.docx`. Roughly 18 are durable
reference, including `phase4_architecture_spec.md`, `design_audit.md`,
`negsteer_notes/cluster_submission_guide.md` and
`inventory/01_module_map.md` through `06_ubiquitous_language.md`. Seven are
point-in-time tool dumps, including `inventory/07_ruff_report.txt`, which §3.3
supersedes.

README references notes at `:11`, `:15`, `:16`, `:69`, `:71`, `:81`, `:125` and
`:132-139`. Several `tests/*/README.md` files reference `notes/inventory/`
directly, so this is a repointing job rather than a deletion.

Upstream hit the same problem from the other side. Its PR20 removed about sixty
citations of internal notes from code comments, because
`docs/development_notes/` was gitignored and the files existed only on one
machine. The equivalent check here is whether any tracked file cites a
`notes/` path that is itself untracked.

---

## 7. Decisions that change the scope, and that I am not taking alone

**7.1 Does the engine become a submodule?** `negative-steering-analyses` pins
`structure-negative-steering` as a submodule at `engine/`, per its
`.gitmodules`. That is the model that makes §4 a one-time migration rather than
a recurring port. Against it: this repo's copies have diverged in both
directions, its `reversion.py` is pre-split, and its monolith holds the stages
upstream has already broken up. Keeping the vendored copy means §4's six
patches now and the same exercise after the next upstream change.

**7.2 Do the campaign analyses split out?** Only `pikp1_avrpikf` has an
`analyses/` tree, and it holds 579 of the 598 tracked files in that campaign,
across 12 numbered analyses plus an `_archive/` of 3 retired ones. The other
four campaigns have `inputs/` and `runs/` only. Splitting mirrors what
`negative-steering-analyses` did. Keeping them means this repo carries both the
pipeline and the thesis analyses.

**7.3 Is the effector template fixed, or recorded and left?** §4.1 is a
two-line change whose consequence is that every prediction made after it
differs from every prediction made before it. Upstream made the template
opt-in and defaulted it off, set all 167 of its configs to `false`, and
published no rerun. Its justification was that every committed number was
produced template-free, so the results stand. The same route is available
here. The AF3 arm has the same net outcome by a different mechanism.
`modules/negsteer_af3_nomsa.nf:192-193` emits `"queryIndices": []` and
`"templateIndices": []`. `notes/things_to_fix.md` records that AF3 3.0.0
builds template features by iterating that mapping, so an empty one leaves the
positions and masks zero. Since `things_to_fix.md` also records a decision
already taken to remove the AF3 stage, fixing its template may be moot.

---

## 8. Analyses, if they stay here

**8.1 Add a render harness.** There is no `_quarto.yml`, no `Makefile` and no
render script. Each of the 12 numbered analyses is an independently rendered
Quarto document with its `.html` and `_files/libs/` committed beside it, driven
by ad-hoc `make_*.py` and `*_extract.py` helpers. Upstream's reproduction check
renders all nine of its analyses and compares every committed CSV and figure
byte for byte. Nothing here can make that claim, because nothing renders them
together.

**8.2 Stop three extractors silently shrinking committed CSVs (upstream PR11).**
Each globs an untracked pipeline output directory and overwrites a tracked CSV
with whatever it finds, exiting 0:

| script | glob | writes (tracked) |
|:---|:---|:---|
| `01-design-region-backbones/dssp_extract.py:31` | `runs/.../results/rfdiffusion/design_*.pdb` | `design_secondary_structure.csv` |
| `01-design-region-backbones/structural_correspondence_extract.py:88` | same | `structural_correspondence.csv`, `..._per_design.csv` |
| `03-mpnn-sequence-design/depth_extract.py:97` | same | `loop_depth.csv` |

I confirmed the first one end to end. `git ls-files` on
`experiments/campaigns/pikp1_avrpikf/runs/crystal_full_test_contig/results/rfdiffusion/`
returns 0 files, and `design_secondary_structure.csv` is tracked. So on any
checkout without that local output, a rerun replaces a 5189-line committed CSV
with a two-row one and reports success. Upstream's guard compares against what
is committed and refuses to shrink it, with `--allow-partial` to override.

**8.3 Repoint one stale analysis script.**
`06-confidence-and-biophysical-metrics/sequence_match_extract.py:34-35` sets
`CONF` to `../04-confidence-cost` and `CAMPAIGN` to `../05-campaign-confidence`.
Neither exists. The sibling directories are `04-negative-steering-outcomes` and
`05-steerability-correlates`. Those two names are upstream's analysis names, so
this file was copied from there and its sibling references were never
repointed. It does not crash on import. I ran `--help` and it printed its
usage, so the failure arrives when it reads those paths.

**8.4 Deduplicate the reference structures (upstream PR2).** There are 207
tracked `.pdb` files holding 127 distinct contents, so 80 are extra copies. One
hash appears eight times: the AVR-PikF reference complex, stored under
`inputs/` twice, under four different `analyses/*/data/` directories, under
`_archive/`, and as `tests/pose_solver/data/PikF.pdb`. Upstream's reasoning was
that four copies of the ground truth is four chances for two arms to disagree
about what truth is. It consolidated into `data/reference_complexes/` and
repointed 166 configs.

**8.5 Rewrite the prose (upstream PR6).** The three largest analysis documents
all sit in `05-steerability-correlates/`:

| file | prose words | mean words/sentence | semicolons | em dashes |
|:---|---:|---:|---:|---:|
| `variables-glossary.qmd` | 1159 | 16.3 | 4 | 22 |
| `boltz-confidence-by-stage.qmd` | 1092 | 25.4 | 2 | 18 |
| `05-steerability-correlates.qmd` | 1105 | 22.6 | 0 | 16 |

Upstream's constraint was no colon, semicolon or dash joining two independent
clauses, and it checked that every code chunk stayed byte-identical and every
figure regenerated unchanged.

**8.6 Retire the superseded steerability document.**
`notes/things_to_fix.md` records that
`05-steerability-correlates.qmd` is the LightGBM and SHAP version, which
defines a hit as `cross_tier` in A/B/C across all 128 designs. It also records
that `steerability_features.py` replaced it with a per-feature Mann-Whitney
over the 106 steered designs. The three `shap` F821 errors in §3.3 are in that
document's `helpers.py`, which is consistent with it being dead code.

---

## 9. Test-fixture weight, recorded rather than proposed

`tests/` holds 2910 tracked files. `tests/negative_steering/` is 2537 of them,
and `tests/negative_steering/example_output_files/` alone is about 2511, a
committed copy of eleven per-design run output trees. 624 tracked files under
`tests/` contain a literal `/hpc-home` or `/nbi/` path, mostly inside captured
`input.yaml` and `result.json` fixtures.

Upstream's PR16 took the opposite approach, building fixture run trees whose
answers are stated up front, with the layouts taken from a real run rather than
from docstrings. Its stated reason was that testing against the real tree ties
the suite to a 15 GB cluster directory. Such a test only proves the code still
agrees with its own last output.

Whether to follow that here is a larger question than parity.
`notes/inventory/17_storage_bloat_audit.md` already covers the storage side, so
this section is a record rather than a proposal.

---

## Sequencing

1. §1 prerequisites, then §2 to get a green suite with an honest number.
2. §3 coverage, lint and CI, so everything after is gated.
3. §4 the six engine patches, each with a test that fails before it.
4. §5 site configuration, then §6 repository setup.
5. §7 answered before §8 starts, because 7.2 decides whether §8 happens here.
