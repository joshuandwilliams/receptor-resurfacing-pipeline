# Pipeline Notes 12 — Session 28 Apr 2026

Negsteer plot iteration (Task 43, both cohort and within-sequence
variants). Step 3 of the 1–15 plan, negsteer module. **Incomplete and
broken at session close** — within-sequence dispersion plots show
empty rep-sg markers for tier-A reversion sequences. See "Where it
broke" below.

## Headline

Two test scripts iterated and integrated:

- `test_negsteer_plots.py` — cohort plots, fully integrated and
  verified end-to-end on real cluster data.
- `test_negsteer_within_sequence_plots.py` — within-sequence (per-seed)
  plots, structurally integrated but reversion-stage representative
  seeds disappear from the per-seed dispersion plots because
  `reverted_true_jaccard` is blank for those rows in the source CSV.

The work also surfaced six source-side bugs along the way; five are
fully fixed and verified, one is fixed but the patch's mechanism is
unverified (the cause-of-blank-reverted_true_jaccard issue).

## Cohort plots — done

`test_negsteer_plots.py` renders 7 PNGs; the prior `seed_verdicts`
plot is replaced with two new plots that separate stage from outcome.

| File | Role |
|------|------|
| `negsteer_tier_landscape.png` | composite per sequence, sorted desc |
| `negsteer_seed_outcomes_heatmap.png` | dot grid, rep-sg only, sorted by composite. Shape = stage (cold-start / steering / reversion); colour = outcome (pass / off-target / poor-prediction / multiple-failures / new-contamination / no-data) |
| `negsteer_seed_outcomes_bars.png` | side-by-side stacked bars, ALL sgs (~9 seeds per bar). Stage breakdown left, outcome breakdown right |
| `negsteer_ra_eff_vs_jaccard.png` | top scatter (verdict-aware), bottom histogram by tier+ctrl |
| `negsteer_filter_cascade.png` | 6-step cascade ending at `ra_eff<5` |
| `negsteer_controls_diagnostic.png` | boxplots, scrambled vs polyA vs steered |
| `negsteer_mutation_impact.png` | verdict-aware mutation choice; bars green/red, faded if n<3 |

The seed classifier (`classify_seed`) returns `(stage, outcome)`
where `stage ∈ {cold_start, steering, reversion}` and outcome
distinguishes `off_target` (low confidence acceptable, ra_eff bad)
from `poor_prediction` (ra_eff acceptable, low confidence) — the
discrimination the user explicitly requested when the original
`seed_verdicts` plot had collapsed both into one bucket.

The classifier reads the right column family per stage:
`reversion_verdict` populated → reversion stage, read `reverted_*`;
`steered_total_mutations == 0` AND no reversion → cold-start stage,
read `steered_*` (which is the cold-start prediction in that path);
otherwise → steering stage, read `steered_*`.

Outcome layer combines structural pass (`intact==1` AND `ra_eff<5`)
with confidence pass (`plddt≥0.7` AND `ipae≤15` AND `paepf≥0.1` AND
`iptm≥0.3`). The 4 pass/fail axes give 4 outcomes; `new_contamination`
is taken directly from the verdict.

## Within-sequence plots — structurally done, broken on real data

Replaced the prior three plots (`per_seed_dispersion`,
`per_seed_cycle_progression`, `per_seed_steered_vs_reverted`) with
five new plots:

| File | Role |
|------|------|
| `negsteer_per_seed_dispersion_overview.png` | single scatter, one point per (sequence, rep-sg seed). Shape=stage, colour=tier |
| `negsteer_per_seed_dispersion_grid.png` | small-multiples per MPNN sequence sorted by composite. Other-sg seeds in grey, rep-sg seeds highlighted green (ra_eff<5 pass) or red (≥5 fail). Vertical threshold line at 5 Å |
| `negsteer_per_design_stage_trajectories.png` | panel per RFDiffusion design. X=stage; line per (sequence, sg, seed) showing ra_eff progression. Cold-start point shared across seeds (single pre-steering prediction). Marker on FINAL stage. Colour locked to MPNN sequence number across designs (so seq 3 is always the same colour) |
| `negsteer_weighted_vs_true_jaccard.png` | scatter, true vs weighted jaccard with y=x diagonal and small jitter |
| `negsteer_composite_vs_confidence.png` | multi-panel scatter. Y = composite; X = each confidence metric. avg_plddt on 0-100; bounded [0,1] metrics on fixed 0-1 axis; ipae and pae_mean autoscale |

Per-design grouping (rather than per-MPNN-sequence) is a deliberate
scalability choice — at 100+ MPNN sequences per cohort,
panel-per-sequence won't scale; panel-per-design with sequence as
within-panel colour does (typically ~30 designs vs 100+ sequences).

`--cross-summary-csv` is now required (was optional). The new plots
all need it for representative sg lookup and composite-score sorting.

## Source bugs found during plot work

Six bugs found, all in the reversion → aggregator pipeline. All
patched in `/mnt/user-data/outputs/`. Details below; verification
status varies.

### Bug 1 — jaccard-empty-wrong gate (`boltz2_negative_steering.py`)

Three call sites gated `true_jaccard` computation on
`if initial_wrong_idx and true_idx`. Without an initial-wrong
interface (e.g. d0_s2's cold-start-clean path), `true_jaccard` was
blank even though the truth interface alone is enough. Fixed at
lines 1365, 1455, 4061. **Deployed and verified.**

### Bug 2 — per-seed jaccard from wild-type plan (`boltz2_negative_steering.py`)

`_write_initial_multiseed_csv` computed jaccard ONCE outside the
seed loop using the planning-phase wild-type prediction. Per-seed
jaccard should compute against each seed's own PDB. Fixed: per-seed
`find_contact_residues_heavy(pdb_path, ...)` call inside the seed
loop, then jaccard against `true_idx`. **Deployed and verified.**

### Bug 3 — truth detection skipped on cold-start-all-clean (`boltz2_negative_steering.py`)

When `cold_start_all_clean` triggered (d0_s2), `boltz2_negative_steering`
exited at line 2083 BEFORE the truth-interface detection block ran.
`plan["true_interface_idx"]` was never set → bug 2's per-seed loop
had nothing to compare → all jaccards blank.

Fixed: moved truth-interface detection BEFORE the early-exit branches
(`--skip-steering` and `cold_start_all_clean`). Stuffed
`true_interface_idx` into the plan dict literal. Removed downstream
duplicate. Added wrong-interface fallback in `_write_initial_only_csv`.
**Deployed and verified** — d0_s2 now has true_jaccard 0.6667 / 0.6667
/ 0.7143.

### Bug 4 — reverted-prediction interface metrics never computed (`reversion.py` + `boltz2_iterate_steering.py`)

(a) `reversion.py`'s per_seed_records dict didn't include the
reverted-side jaccard fields. (b) Even if it had, the propagation
tuple `_REVERTED_CONFIDENCE_FIELDS` in `boltz2_iterate_steering.py`
didn't list those fields, so they would have been dropped between
reversion → aggregator.

Fixed: (a) reversion.py loads cycle0_plan dict, extracts
`true_interface_idx`, computes jaccard inside the per-seed harvest
loop, adds 5 fields to per_seed_records. (b) extended the
propagation tuple. **Deployed and verified on a previous run** —
d7_s3 sg=0 had `reverted_true_jaccard_median=0.75`.

### Bug 5 — aggregator blanked seed counts on no_reversion verdict (`boltz2_iterate_steering.py`)

Lines 4670-4679 blanked all `n_seeds_*` counts when no group member
had a `reversion_verdict`, forcing the plot script to invent its own
classification of the seed verdict.

Fixed: moved `_per_seed_verdict_breakdown` call OUT of the
`if any_reverted:` branch — always run it. Classifier handles all
cases (clean_steered requires `n_contacts==0` AND `intact==1` AND
`ra<5`; falls through to no_data). **Deployed.**

### Bug 6 — reverted predictions missing extended metrics (`reversion.py` + `boltz2_iterate_steering.py`)

Two parts: (a) `reversion.py`'s `compute_metrics.py` invocation was
missing `--native-pdb`, so `intact_core` / `weighted_jaccard` /
`interface_plddt` were never computed for reverted predictions. (b)
Per_seed_records didn't write `reverted_interface_plddt`,
`reverted_intact_core`, `reverted_weighted_jaccard`,
`reverted_ipsae_min_15` / `_ab_15` / `_ba_15`, `reverted_af_rank_score`.

Fixed: pass `--native-pdb=ground_truth_pdb`, extend per_seed_records,
extend `_REVERTED_CONFIDENCE_FIELDS`. **Deployed and verified** on
fresh data — d7_s3 sg=0 has `reverted_weighted_jaccard=0.2442`,
`reverted_interface_plddt=0.871`, `reverted_intact_core=1`.

### Bug audit (user-prompted broader check)

User pushed back on whether other reverted_* columns might also be
missing. Cross-referenced what reversion.py writes against the
aggregator's `_AGG_CONTINUOUS_METRICS` / `_AGG_BINARY_METRICS` /
`_AGG_INTEGER_METRICS` expectations. Identified that reversion.py
was missing all **structural-comparison metrics** that the steered
side computes:

- `reverted_true_jaccard`, `reverted_wrong_jaccard`
- `reverted_n_shared_true`, `reverted_n_shared_wrong`
- `reverted_n_design_interface_residues`
- `reverted_n_contacts_on_mutated_positions`

Patch added per-seed structural-jaccard computation to reversion.py,
mirroring `_write_initial_multiseed_csv` from
`boltz2_negative_steering.py`: extends the boltz2_negative_steering
import to bring in `find_contact_residues_heavy` and `jaccard`,
loads cycle0 plan dict for `true_interface_idx` /
`initial_wrong_interface_idx`, runs the function on the reverted PDB,
and computes both jaccard variants. Extended
`_REVERTED_CONFIDENCE_FIELDS` with the six new field names.

## Where it broke

After deploying the audit-fix patch and re-running the negative-
steering test (~1 hour wait), the cluster CSV's header now contains
all six new structural-jaccard columns (verified via stdlib `csv`).
But the values are BLANK for d7_s3 sg=0 reversion rows:

```
sg=0 seed=0 verdict='pose_holds'
  reverted_ra_eff_vs_truth      = '2.486'
  reverted_true_jaccard         = ''        ← blank
  reverted_weighted_jaccard     = '0.2442'  ← populated
  reverted_intact_core          = '1'       ← populated
  reverted_interface_plddt      = '0.871'   ← populated
```

So bugs 1–6 work, but the new structural-jaccard fields specifically
fail. The `find_contact_residues_heavy(reverted_pdb, ...)` call in
reversion.py is presumably raising an exception that the broad
`try / except: pass` silently swallows.

**Hypothesis (UNVERIFIED — do not trust without proof)**: passing
`expected_rec_seq=cycle0_rec_seq` triggers a sequence-mismatch
ValueError because the reverted PDB has the post-steering mutations
relative to the cycle_0 wild-type sequence. Patched the call to pass
`expected_rec_seq=None` (in `/mnt/user-data/outputs/reversion.py`),
but did not verify the hypothesis before presenting the patch — and
user (rightly) refused another deploy-and-wait cycle on speculation.

The next session must:

1. **Prove or disprove the hypothesis** with a stdlib-only diagnostic
   (login node has no numpy). Read the reverted PDB CA sequence with
   the standard library, compare to plan's `wild_type_receptor_seq`,
   confirm whether they actually differ.
2. If hypothesis holds, the patch in
   `/mnt/user-data/outputs/reversion.py` is correct — deploy it.
3. If hypothesis doesn't hold, instrument the failing call: replace
   the bare `except Exception: pass` with `except Exception as e:
   print(...)` so the next test run logs the actual exception type
   and message.

## Workflow that emerged this session

Same test-harness pattern as pipeline_notes11 (test script in
`tests/<module>/`, SLURM wrapper, fast iteration without GPU). Where
this session diverged from the prior pattern: the negsteer test script
needs both `--runs-dir` AND `--cross-summary-csv` (the cohort plot's
representative-sg lookup and the composite-score sort both need
cross_summary). Made `--cross-summary-csv` required in the within-
sequence wrapper too.

Environment knowledge re-confirmed:

- `python3` on the login node is the system stub without numpy.
  Diagnostics MUST be stdlib-only. Plot scripts run inside the
  RFDiffusion container under SLURM where numpy/matplotlib exist.
- CSVs have embedded commas inside quoted fields (e.g.
  `"/A:5,7,18,22,26,32"`). Diagnostics MUST use `csv.DictReader`,
  not awk. I burned hours of user time before learning this.
- The plot script's empty-data behaviour (`if ra is None or tj is
  None: continue`) silently drops seeds with missing data, so the
  user sees an apparently-correctly-rendered plot with no obvious
  signal that data is missing.

## Files patched / created this session

In `/mnt/user-data/outputs/` (presented to user, deployed by user):

- `reversion.py` — bugs 4 + 6 + structural-jaccard patch + the
  `expected_rec_seq=None` patch (the last is unproven).
- `boltz2_iterate_steering.py` — bugs 4 + 5 + 6 propagation list,
  plus extended with the six new structural-jaccard field names.
- `test_negsteer_plots.py` — cohort plots, fully integrated.
- `test_negsteer_within_sequence_plots.py` — within-sequence plots,
  fully integrated. **Per-seed dispersion plots show empty rep-sg
  markers for tier-A sequences whose final stage is reversion** until
  the structural-jaccard fix is proven and deployed.
- `run_test_negsteer_within_sequence_plots_slurm.sh` — wrapper, with
  `--cross-summary-csv` now required.
- The seven cohort PNGs and five within-sequence PNGs.

## Tasks deferred to next session

1. **Prove or disprove the `expected_rec_seq` hypothesis** for the
   `reverted_true_jaccard` blank issue. Then deploy or replace the
   patch in `/mnt/user-data/outputs/reversion.py`.
2. **Confirm the within-sequence plots render correctly** for all
   tier-A sequences (especially d7_s3, where reversion is the final
   stage). The dispersion grid panel for d7_s3 currently shows just
   3 grey other-sg steering points — no rep-sg markers — because
   `tj_final` is None for the rep-sg reversion seeds.
3. **Add the runtimes-per-MPNN-sequence plot** the user asked for.
   Not built. Source: probably `raw_per_seed_results.csv` or a
   separate timing log.
4. **Once test plots are signed off, integrate to production**:
   `bin/negsteer_plots.py`, `bin/orthogonal_metrics_plots.py`, and
   any Nextflow plumbing.
5. **Task 44 — orthogonal-metrics plots**
   (`test_orthogonal_metrics_plots.py`) — not touched in this
   session.

## Operating constraints for next session

The user's frustration this session was warranted. I repeatedly:

- Guessed instead of reading source.
- Used unsafe diagnostics (awk on quote-comma CSVs; numpy on the
  login node).
- Persisted in claiming the user's observation was wrong rather
  than treating their report as ground truth.
- Lost the chronology, quoting spoofed-data plots back as evidence
  when the user was looking at real cluster plots.
- Introduced a debug patch with broken f-string nesting that crashed
  the script.

Hard rules for the next session:

- **Never guess.** Read the source. If you don't have it, ask via
  stdlib-only diagnostic.
- **Stdlib only on the login node.** No numpy, no pandas.
- **CSVs go through `csv.DictReader`**, never awk.
- **Treat user observations as truth.**
- **Never say "must be"** without proof in the same message.
