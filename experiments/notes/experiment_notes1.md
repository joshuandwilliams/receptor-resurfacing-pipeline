# Experiment Notes 01 — Scaffolding the experiments/ workflow and first overnight runs

First session establishing the `experiments/` directory as a home for
binder design campaigns, building the supporting tooling, and
submitting the first two overnight pipeline runs (Strategy 4a and
Strategy 4b for `pikp1_avrpikf`). Session ended with the 4a run
completed but mostly tier-none, and the 4b run discovered to have
been overwriting the same output directory due to a copy-paste error
in its params file.

## Headline

By session close:

- A new `experiments/` directory is scaffolded, branched, and
  documented. Layout is `scripts/`, `param_derivation/`, `inputs/`,
  `campaigns/`, with each campaign holding `inputs/`, `runs/`, and
  `analyses/`. Every campaign run lives in its own `runs/v<N>_<slug>/`
  subdirectory with its own `params.yml`.
- A full `experiments/README.md` codifies framing, the six-stage
  campaign lifecycle, and naming/placement conventions. Future-self
  audience.
- The first contig derivation script, `experiments/param_derivation/
  contigs_4a.py`, is built and produces a real contig for the
  pikp1_avrpikf campaign.
- A `_path_setup.py` helper makes `bin/` importable from any campaign
  script so production pipeline logic doesn't get duplicated.
- Two pipeline runs submitted overnight (4a and 4b for pikp1_avrpikf).
  4a completed cleanly. 4b was unintentionally overwriting 4a's
  output directory because the params file's `outdir` and
  `project_name` weren't actually updated — paste failure on the
  user's side. Both need re-running.
- Several infrastructure issues uncovered and fixed along the way:
  the bare `experiments/` line in `.gitignore` that was silently
  hiding everything; a missing `--chdir` directive interaction in
  `run_pipeline.slurm.sh`; missing `experiments/` excludes in
  `sync_to_hpc.sh`.

## Architecture decisions

Captured here because they shape every campaign that follows.

### One worktree per branch, not branch-switching

A `git worktree` was set up at `~/Documents/GitHub/receptor-pipeline-experiments/`
running the `experiments` branch alongside the existing
`receptor-resurfacing-pipeline/` worktree on `remediation`. This lets
the user run two Claude Code sessions concurrently without uncommitted
work on one branch blocking the other. The user is solo on this repo,
so the slight overhead of two folders pays for itself in cognitive
isolation between the two streams.

Implication: `sync_to_hpc.sh` will sync from whichever worktree it's
launched from, and `HPC_DEST` is hard-coded — so syncing from
remediation while experiments is in flight would wipe `experiments/`
from HPC. The discipline is: only sync from the experiments worktree
during campaign work. Documented but not enforced.

### Campaign = receptor–effector pair, runs nested inside

Each campaign is one binder–target pair (`pikp1_avrpikf`,
`pikp1_avrpia`, etc.). Multiple parameter explorations against the
same pair go in `runs/v<N>_<slug>/` subdirectories rather than as
sibling campaigns. The campaign-level `analyses/` directory is for
cross-run comparison; each run has its own `analyses/` inside
`runs/<run>/` for run-specific work.

The user originally proposed campaign = single run, but switched to
this nested model when the conversation surfaced the realistic
likelihood of multiple parameter explorations per pair (different
contig strategies in particular).

### Naming convention

- Campaigns: `<receptor>_<effector>` lowercase (`pikp1_avrpikf`).
- Runs: `v<N>_<slug>`, monotonic version + descriptor
  (`v1_4a`, `v2_4b`).
- Inputs: `<campaign>_<purpose>` (`pikp1_avrpikf_complex.pdb`,
  `pikp1_avrpikf_4a_contig.txt`).
- Provenance: `<filename>.notes` plain text alongside any input
  whose filename can't carry full provenance.

User caught a mid-session inconsistency: the campaign-naming
examples showed `avrpikF` (capital F) while the actual directory
was `avrpikf` (lowercase). Fixed before commit.

### Selective gitignore, not whole-folder ignore

The original `.gitignore` had a bare `experiments/` line that was
silently hiding the entire scaffolding tree from git. Removed that
and replaced with narrower rules that ignore only run outputs
(`experiments/campaigns/*/runs/*/results/`, `analyses/`, `work/`,
`.nextflow*`, `*.out`, `*.err`). This means READMEs, scripts,
`__init__.py`, params files, and contig artefacts ARE tracked —
which is the right shape, since those are the bits that benefit
from version control.

### Inputs go in the campaign, not the run

Curated inputs (assembled complex PDB, contig files, residue lists
for 4b/4c) live in `<campaign>/inputs/` and are reused across runs
of the campaign. The single exception is `params.yml`, which lives
at the run root because it's run-specific by definition.

Shared inputs reusable across campaigns (e.g. an experimental
structure of a receptor used in multiple campaigns) go in
top-level `experiments/inputs/`. None of those exist yet.

### Empty directories pinned with .gitkeep

Empty directories don't survive a fresh worktree or clone (git
doesn't track them). The campaign skeleton uses `.gitkeep` files in
empty subdirectories to ensure the structure is durable. The user
discovered this the painful way when creating the experiments
worktree — the originally-created campaign directories had vanished
because they were empty.

## The contig derivation script

`experiments/param_derivation/contigs_4a.py` implements Strategy 4a
(interface-facing side chains, automatic geometric selection). The
spec went through a critical inversion mid-build — see "Bug 1" below.

### Final semantics

Reads a complex PDB, identifies binder residues whose heavy
side-chain atoms come within a distance cutoff (default 5.0 Å) of
the target chain, and emits an RFDiffusion contig string with those
residues as **design regions** and the rest as **anchor blocks**.

The intent: redesign the binding interface itself, keep the rest of
the binder native. NOT the opposite, even though this was the user's
description in passing and got encoded wrong on the first build.

### Reused helpers from bin/

- `is_surface_exposed` from `bin/boltz2_negative_steering.py` for
  the outward-facing test (Cβ-vector heuristic against
  neighbour-centroid). Glycine excluded.
- `find_contact_residues_heavy` from same file for the
  heavy-side-chain-distance test. Note: this helper compares all
  heavy atoms (backbone + side chain) on both sides, while strict
  4a wants side-chain-heavy on the binder. The pre-filter by
  outward-facing residues mitigates the difference in practice.
  Documented in the script's `.notes` output.

### CLI flags

```
--complex-pdb PATH       Required.
--binder-chain CHAR      Required.
--target-chain CHAR      Required.
--output PATH            Required.
--cutoff FLOAT           Default 5.0 Å.
--min-multiplier FLOAT   Default 0.7.
--max-multiplier FLOAT   Default 1.5.
--merge-threshold INT    Default 3.
```

### Length-range formula for design regions

For each design region, the token emitted is `<min>-<max>` where:

- `native_length == 1` → `1-1`
- `native_length > 1`:
  - `min = max(1, round(native × min_multiplier))`
  - `max = round(native × max_multiplier)`

Asymmetric (more upward stretch than downward compression) because
a too-short region disrupts fold; a slightly-too-long region just
adds slack.

### Bridging

After identifying interface residues, internal anchor blocks
(bordered by design regions on both sides) of length
≤ `--merge-threshold` get absorbed into the surrounding design.
The merged region's native length becomes the sum of (design +
absorbed anchor + design), and the formula is applied to that sum.

This handles the helical periodic pattern where every 2nd or 3rd
residue is interface-facing — without bridging, a single helix
becomes a string of `1-1` tokens with no design freedom.

Terminal anchor blocks (one side of the binder with no design region
beyond it) are NEVER bridged regardless of length. The bridging rule
only applies to internal anchors.

### Future enhancement noted in module docstring

Density-based merging — sliding-window analysis of design-residue
density, marking high-density windows as design regardless of
internal anchor lengths — is a noted alternative. Not implemented
because the fixed-threshold approach handles the dominant case
(helical periodicity) cleanly, and the density approach has more
parameters to tune.

### Output for pikp1_avrpikf

Final contig (with `--merge-threshold 2`):

```
A1-2/2-4/A6-31/6-12/A40-42/5-11/A50-67/8-16 C
```

The user inspected the design regions in ChimeraX and decided to
override the default threshold from 3 to 2 specifically for this
campaign. The single anchor that changed at threshold 2 was
`A40-42` — a 3-residue loop between two β-strands at the interface.
Keeping it as a fixed anchor gives RFDiffusion a structural pivot
to design around rather than letting the entire strand-loop-strand
region float.

This is a real campaign-design insight worth carrying forward: the
all-β-strand interface here may benefit from threshold 2 as a
default; α-helical interfaces likely want threshold 3. Worth
codifying in campaign READMEs as they get populated.

## Six issues fixed this session

### Issue 1 — Strategy 4a contig built with inverted semantics

The user's original description of 4a included the phrase "select
all binder residues with side-chains that face into the interface,"
and the assistant interpreted "select" as "select to anchor" because
earlier conversation framed the four strategies as "anchor sets."
The first-build contig had interface-facing residues as anchors and
the gaps as design — the opposite of what 4a is for.

Caught when the user ran the ChimeraX visualisation: the cyan design
regions were on the back face of the binder, not the interface. Real
mental-model mismatch, not a code bug — the geometry was right; the
assignment was inverted.

**Fix**: rewrote the contig assembly logic to treat interface-facing
residues as design and gaps as anchors. The "gap-filling" concept
from the README (which had also been encoded backwards) was replaced
with bridging — absorbing short anchor blocks INTO surrounding design
regions. N-terminal trimming was removed (under the inverted logic
the N-terminus is now an anchor by default, nothing to trim). The
README's Stage 4 framing also got things backwards and is on a list
of follow-up edits.

The lesson for future strategy descriptions: "select" / "designate" /
"choose" are ambiguous between "fix" and "redesign." Always pin down
which.

### Issue 2 — Bare `experiments/` rule in `.gitignore` hiding the entire scaffolding

After Claude Code finished the initial scaffolding pass, `git status`
reported "nothing to commit." Investigation showed the root
`.gitignore` had a bare `experiments/` line under "Archive / scratch
directories" — leftover from when `experiments/` was used as a scratch
dump. Every README, every `__init__.py`, the entire tree was being
ignored.

**Fix**: removed the bare line, added narrower rules for run outputs
only:

```
experiments/campaigns/*/runs/*/results/
experiments/campaigns/*/runs/*/analyses/
experiments/campaigns/*/runs/*/work/
experiments/campaigns/*/runs/*/.nextflow*
experiments/campaigns/*/runs/*/*.out
experiments/campaigns/*/runs/*/*.err
```

Source files now tracked correctly.

### Issue 3 — Empty directories don't survive worktree creation

When the user created the experiments worktree, the five originally-
created campaign directories vanished because they were empty. Git
materialises only what's tracked; empty directories aren't tracked.
Briefly looked like branch-switching had deleted them — it hadn't,
they just hadn't existed in git.

**Fix**: added `.gitkeep` to every empty intentional subdirectory in
`experiments/`. Documented the failure mode in the experiments README
and in the discussion of branch hygiene with the user.

### Issue 4 — `run_pipeline.slurm.sh` had a leftover `--chdir`

First overnight submission failed within seconds with
`ERROR: params file not found: /hpc-home/.../receptor-resurfacing-pipeline/params.yml`.
The user had submitted from inside `runs/v1_4a/` with `./params.yml`
as a relative argument. Investigation showed
`run_pipeline.slurm.sh` had:

```
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline
```

SLURM did the chdir BEFORE running the script, so `./params.yml`
resolved against the repo root, not the launch directory. The same
directive also caused slurm logs to land at the repo root instead
of inside the run directory.

**Fix**: removed the `--chdir` line. The script self-anchors
internally for finding `main.nf`, so the directive was redundant.
Likely propagated from the per-module test slurm scripts via
copy-paste.

After fix: launch from inside `runs/v1_4a/` works, slurm logs land
inside the run directory, params relative path resolves correctly.

### Issue 5 — `sync_to_hpc.sh` missing `experiments/` excludes

The script had excludes for `tests/*/results/`, `tests/*/work/`, etc.
but no parallel rules for `experiments/campaigns/*/runs/*/`. With
`--delete` enabled, syncing from Mac while a campaign was producing
output on HPC would wipe the run's `results/` and `work/` from HPC
mid-run.

Caught before the second sync. The user spotted it and asked whether
the sync would clobber the running 4a outputs. It would have.

**Fix**: added the parallel exclude family:

```
experiments/campaigns/*/runs/*/results/
experiments/campaigns/*/runs/*/work/
experiments/campaigns/*/runs/*/tmp/
experiments/campaigns/*/runs/*/.nextflow*
experiments/campaigns/*/runs/*/*.out
experiments/campaigns/*/runs/*/*.err
```

Header comment block updated to document.

### Issue 6 — User's v2_4b params.yml never actually updated

The user reported in the morning that 4a's results showed mostly
tier-none, but worse: 4a's outputs were visible in BOTH `v1_4a/` and
`v2_4b/` directories. Diagnosis: the user had copied `v1_4a/params.yml`
to `v2_4b/params.yml` and then thought they were pasting an updated
version from chat into the file. The clipboard "copy" action from
the chat UI sometimes silently fails. The result was that v2_4b
inherited v1_4a's `outdir`, so both runs wrote to the same place.

This is not a code bug. It's a UX failure mode worth flagging in
operating constraints.

**Fix path** (next session, not this one):

1. Delete the corrupted v2_4b output directories on HPC.
2. Verify the params files now actually differ:
   `diff v1_4a/params.yml v2_4b/params.yml`
3. Resubmit both runs.

Before re-running, the v1_4a results should be considered tainted —
they may include 4b's writes if any landed before 4a finished. Worth
checking the timestamps on the output files to see if there's
contamination.

## The 4a results — what 4a is telling us

The 4a run completed before the v2_4b discovery, and the user
inspected it. Pattern is striking:

- 3/128 surviving designs after the full pipeline. All tier-none.
- 58/118 negsteer failures from `complex_pLDDT < 0.7` (49%).
- 26 more from `ipAE > 15` (22%).
- 22 more from `ra_eff < 5` (19%).
- 3 more from `pae_pass_frac < 0.1`.
- 9 sequences passing negsteering. Almost nothing reaches reversion.
- pLDDT distribution clusters in the 0.6–0.8 band — sitting right
  around the threshold, not separating into clear pass/fail.
- `pae_pass_frac` mostly ~0.

The user proposed a hypothesis: Boltz-2 without MSA may simply not
be reliable enough for a binder this aggressively redesigned, which
would invalidate the whole negative-steering approach.

The session ended before this was diagnosed. The discussion-worthy
question is: which of the following is happening?

1. The designs are genuinely bad. Contig too aggressive, or contig
   in the wrong place. Boltz is correctly rejecting them.
2. The designs are OK but Boltz-no-MSA can't tell. Confidence is
   degraded enough that the threshold catches everything.
3. The designs are mixed and Boltz-no-MSA is too noisy to separate
   them. Threshold tuning could rescue genuinely good designs.
4. The thresholds are calibrated for a different distribution
   (production binders against simpler targets) and don't transfer
   to this binder–effector pair.

The discriminating experiment proposed: take a small sample of
designs spanning the failure modes plus the 9 survivors, and
re-predict with AF3-WITH-MSA on the receptor. If AF3-with-MSA gives
high confidence on designs Boltz-no-MSA flagged as bad, it's a
calibration issue. If AF3-with-MSA also flags them as bad, the
designs really are bad.

Worth doing before deciding what to change for the next round.

A second informative thing: look at the 9 survivors specifically.
What metric profiles or sequence features distinguish them? If they
cluster in a specific region of the contig, that's a signal about
which residues to handle differently next time.

User's note on 4b: "4b always produces lots of successes, because
the design regions are regions that are variable naturally and the
design regions are smaller so the whole interface doesn't necessarily
change." This means the eventual 4b vs 4a comparison will be
informative — if 4b succeeds while 4a fails, the diagnosis is
"4a's contig is too aggressive for this binder," not "Boltz-no-MSA
is broken."

## What's running where

At session close, two runs were submitted but the v2_4b file had
the bug above:

- `experiments/campaigns/pikp1_avrpikf/runs/v1_4a/` — completed,
  results inspected. May be partially contaminated by v2_4b's
  writes; needs verification.
- `experiments/campaigns/pikp1_avrpikf/runs/v2_4b/` — was a
  duplicate of v1_4a due to params.yml not actually updating.
  Outputs unusable. Need to clean and resubmit with corrected
  params.

Strategies 4c and 4d weren't reached. The user planned to do them
in the morning before deciding the post-mortem on 4a/4b was the
priority instead.

## Files created or modified this session

In the experiments worktree (committed and pushed):

- `experiments/README.md` — full docs for the experiments workflow
  (framing, layout, six-stage lifecycle, conventions). ~195 lines.
- `experiments/_path_setup.py` — repo-root discovery helper that
  adds `bin/` to `sys.path`. Walks upward from `__file__` looking
  for `main.nf` AND `bin/` together. Caches `REPO_ROOT` at import.
  Idempotent.
- `experiments/scripts/__init__.py` — empty, makes `scripts/`
  importable.
- `experiments/scripts/README.md`, `experiments/inputs/README.md`,
  `experiments/param_derivation/README.md`,
  `experiments/campaigns/README.md` — placeholder READMEs.
- `experiments/param_derivation/contigs_4a.py` — Strategy 4a contig
  derivation script.
- `experiments/campaigns/pikp1_<five-effector>/` — five campaign
  skeletons, each with `inputs/`, `runs/`, `analyses/`, README, and
  `.gitkeep` files in empty directories.
- `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_complex.pdb`
  — assembled complex PDB, generated in ChimeraX from a two-MODEL
  source PDB via `combine #1.1 #1.2`.
- `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_complex.notes`
  — provenance for the complex.
- `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_4a_contig.txt`
  and `.notes` — contig string and derivation provenance.
- `experiments/campaigns/pikp1_avrpikf/runs/v1_4a/params.yml` — full
  pipeline params for the 4a run.
- `experiments/campaigns/pikp1_avrpikf/runs/v2_4b/params.yml` — was
  meant to be the 4b run's params but didn't actually update from
  v1_4a. Needs fixing.

In `pyproject.toml`:

- Added empty `experiments` optional-dependencies extra.

In `.gitignore`:

- Removed bare `experiments/` rule.
- Added narrower `experiments/campaigns/*/runs/*/...` excludes.

In `scripts/sync_to_hpc.sh` (committed on remediation, merged into
experiments):

- Added `experiments/campaigns/*/runs/*/...` excludes mirroring
  `tests/*/...`.

In `run_pipeline.slurm.sh` (committed on remediation, merged into
experiments):

- Removed leftover `#SBATCH --chdir=...` directive.

## Branching workflow

Two worktrees, two branches:

- `~/Documents/GitHub/receptor-resurfacing-pipeline/` — `remediation`
  branch. The original codebase work continues here.
- `~/Documents/GitHub/receptor-pipeline-experiments/` — `experiments`
  branch. All campaign work happens here.

When a fix is needed in shared infrastructure (`scripts/sync_to_hpc.sh`,
`run_pipeline.slurm.sh`, `.gitignore`, `pyproject.toml`), the pattern is:

1. Edit on the remediation worktree.
2. Commit, push.
3. Switch to experiments worktree.
4. `git fetch origin && git merge origin/remediation`.
5. Push.
6. Sync.

Two such fixes happened this session (`--chdir` removal and the
`sync_to_hpc.sh` excludes). Both merged cleanly with no conflicts.

The eventual integration of experiments back into remediation is a
future-session question, deferred until campaigns are stable.

## Tasks for the next session

1. **Diagnose the 4a result.** AF3-with-MSA spot-check on a sample
   of designs across the failure modes plus the 9 survivors. Pinpoint
   whether the issue is the designs, the predictor, or the
   thresholds.
2. **Re-run 4b cleanly.** Verify the params actually differ from
   4a's, clean the output directory on HPC, resubmit. Use the diff
   command before resubmission to catch any further paste failures.
3. **Verify 4a outputs aren't contaminated** by v2_4b's parallel
   writes. Check timestamps; if contamination is suspected,
   re-run 4a too.
4. **Fix the README's Stage 4 description.** Currently frames
   strategies as "anchor sets" when they're actually about which
   residues are designable. Same inversion that bit the script.
5. **Add 4d and 4c contig derivation scripts.** 4d is mechanically
   close to 4a. 4c is more involved (needs a curated input).
6. **Codify the threshold-2-for-β-strand-interface insight** in
   campaign READMEs once those start being populated.

## Operating constraints — note for next session

A few things the user reasonably called out this session:

- **Filler phrases like "Wait — actually..." waste tokens and feel
  uncommitted.** When the user has already done a step, the response
  should pick up from there cleanly, not re-narrate the decision
  process.
- **Long commit messages should be provided as a separate paste
  block** so the user can copy them into nano without manually
  parsing them out of prose. Adopted mid-session.
- **The "copy" button in the chat UI sometimes silently fails.** The
  user has to verify the destination file actually contains what
  was supposedly pasted. Caused the v2_4b duplicate-output incident.
  Diff-against-source before any submission that depends on a
  recent paste.
- **Never sync from the wrong worktree.** `sync_to_hpc.sh` is
  branch-blind and uses `--delete`. Always sync from the experiments
  worktree during campaign work.
- **Verify the actual params.yml on disk** before submission, not
  the one Claude generated. The diff between two run params files
  should show only the expected differences (project_name, outdir,
  contigs).

Hard rules carried forward:

- Never guess. Read the source.
- Verify with a real command, not a claim.
- Treat user observations as truth.
- Never say "must be" without proof in the same message.
- For shared-infrastructure fixes, edit on the appropriate branch
  (remediation), then merge. Don't paper over with worktree-local
  edits.
