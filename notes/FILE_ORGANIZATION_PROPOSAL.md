---
status: DRAFT
scope: cleanup-proposal
last_verified: 2026-05-19
---

# File Organization Proposal

A concrete proposal for re-shaping `notes/` and the surrounding
files so that future Claude sessions (and future humans) can find
authoritative information without spelunking through 25+ session
logs.

The trigger for this document: the repo has built up enough notes,
one-off scripts, and superseded design docs that there is no longer
a fast way to answer "what is the current truth about X."  The
`notes/` directory mixes timeless reference material with
session-by-session handover narratives.  Frontmatter and a directory
shape will fix this without losing history.

**No deletions.**  Every move proposed below is `git mv` into an
archive subtree so history and `git log --follow` survive intact.

The companion document `notes/PIPELINE.md` is the new authoritative
how-it-works reference (added in the same session as this proposal).

---

## 1. The problem in one paragraph

Authoritative content (`phase4_architecture_spec.md`, the glossary
`inventory/06_ubiquitous_language.md`, the HADDOCK protocol log) and
session-snapshot content (`pipeline_notes15.md`,
`SESSION_HANDOFF.md`, 14 numbered predecessor `pipeline_notesN.md`
files, 11 numbered `notesN.md` files under `negsteer_notes/`,
multiple `.docx` todo lists) sit at the same nesting level in
`notes/`.  Nothing in filenames or directory structure tells the
reader (Claude or otherwise) which to treat as current.  The
inventory documents `01-17` under `notes/inventory/` are a
point-in-time snapshot from `2026-05-12` but read as if they were
current.  Old planning prompts (`12_phase_2_prompts.md`,
`13_prompt_3_hpc_tier_tests.md`, `16_prompt_phase_2_9.md`) and
static-analysis dumps (`07-09`) clutter alongside the actually
useful `06_ubiquitous_language.md`.

## 2. The four-layer fix

### Layer 1 — Frontmatter labels on every `.md` under `notes/`

Add a YAML header to every markdown file in `notes/`:

```yaml
---
status: AUTHORITATIVE | TRANSIENT | HISTORICAL | DRAFT
scope: pipeline-overview | architecture | glossary | subsystem:<name> | strategy | session-log | campaign-log | cleanup-proposal
last_verified: YYYY-MM-DD
supersedes: <path>      # optional
superseded_by: <path>   # optional
---
```

Four labels exhaustively partition the file set:

- **AUTHORITATIVE** — current truth.  Cite this.
- **TRANSIENT** — current state-of-the-work that is expected to be
  overwritten each session (`SESSION_HANDOFF.md`,
  `remediation_state.md`).
- **HISTORICAL** — point-in-time record kept for traceability; do
  not treat as current.  Includes archived session logs, old
  planning docs, grill-me transcripts after their decisions have
  been folded into authoritative docs.
- **DRAFT** — work-in-progress, not yet authoritative.  Becomes
  AUTHORITATIVE on next promotion.

Frontmatter is grep-able and makes "where is the canonical version
of X" answerable by `rg "scope: subsystem:negsteer" notes/`.

### Layer 2 — `CLAUDE.md` at the repo root

Claude Code auto-loads `CLAUDE.md` from the repo root into every
session's context.  This is the cheapest possible mechanism to
give every future session the right starting set without paying
that cost in tokens manually.

Proposed `CLAUDE.md` content (sketch — exact wording TBD on
approval):

```markdown
# Project: Receptor Resurfacing Pipeline

## Read first (in this order)

1. `notes/SESSION_HANDOFF.md` — current session state (TRANSIENT)
2. `notes/remediation_state.md` — branch + phase state (TRANSIENT)
3. `notes/PIPELINE.md` — how the pipeline works, end-to-end (AUTHORITATIVE)

## Authoritative references (current truth)

- `notes/PIPELINE.md`                          — pipeline overview
- `notes/phase4_architecture_spec.md`          — deep-module architecture
- `notes/inventory/06_ubiquitous_language.md`  — glossary
- `notes/codebase_remediation_plan.md`         — strategy
- `notes/subsystems/haddock_protocol_experiments.md` — HADDOCK subsystem
- `notes/subsystems/rfdiffusion_parameters.md`       — RFD parameters
- `experiments/README.md`                            — campaigns layout

## Historical (do NOT treat as current truth)

- `notes/design_audit.md`             — grill-me transcripts; decisions are in PIPELINE.md and the spec
- `notes/archive/**`                  — superseded session logs and planning docs
- `notes/inventory/01-05, 07-17`      — Phase-1 inventory snapshot (2026-05-12)

## Special memory files

The user maintains a Claude-Code auto-memory directory at
`~/.claude/projects/-Users-jowillia-Documents-GitHub-receptor-resurfacing-pipeline/memory/`.
Read these specifically before:

- touching negsteer outcomes → `project_no_reversion_semantics.md`
- writing any contig example → `project_contig_string_format.md`

## When in doubt

Code wins over docs.  Update the relevant authoritative file in the
same commit as any behavior change.
```

### Layer 3 — Directory shape

Reorganise `notes/` so the path itself tells you what kind of
document it is:

```
notes/
  PIPELINE.md                       AUTHORITATIVE
  SESSION_HANDOFF.md                TRANSIENT
  remediation_state.md              TRANSIENT
  codebase_remediation_plan.md      AUTHORITATIVE (strategy)
  phase4_architecture_spec.md       AUTHORITATIVE (architecture)
  glossary.md                       AUTHORITATIVE (rename / symlink of inventory/06)
  FILE_ORGANIZATION_PROPOSAL.md     DRAFT (this file)

  subsystems/                       AUTHORITATIVE per subsystem
    haddock_protocol_experiments.md
    rfdiffusion_parameters.md
    negsteer_pipeline.md            (NEW — currently no canonical doc)
    pose_solver.md                  (NEW — currently only in pipeline_notes15)

  operations/                       AUTHORITATIVE (how-to)
    cluster_submission.md           (verify currency; or archive)
    storage_audit.md                (from inventory/17)

  research_memos/                   reference research
    metric_reliability_memo1.docx
    metric_reliability_prompt1.md

  reference/                        external reference material
    lab_handbooks/
      agroinfiltration_coip_handbook.docx
      protoplast_assays_handbook.docx

  campaigns/                        per-campaign running notes
    pikp1_avrpia_pose_solved.md     (when the current campaign closes; promotes pipeline_notes15)

  inventory_snapshots/              HISTORICAL point-in-time scans
    2026-05-12/                     (the existing 01-17 + 07-09)

  archive/                          HISTORICAL — kept for traceability
    pipeline_session_logs/          pipeline_notes1.md … pipeline_notes14.md
    negsteer_session_logs/          notes1.md … notes11.md (notes11 if folded), merge_map.md
    phase1/                         10_phase_1_synthesis.md, audit_pre_hpc_roundtrip.md
    phase2_planning/                11, 12, 13, 14, 15, 16
    todo_lists/                     todo_list2-10.docx
    weekly_reports/                 weekly_progress_apr15_apr22_1.docx
    grill_me/                       design_audit.md (annotate as HISTORICAL append-only)
```

### Layer 4 — Naming convention for ephemera

Files matching these patterns are by convention **session logs** —
HISTORICAL the moment a successor exists:

- `pipeline_notes\d+\.md`
- `notes\d+\.md`
- `todo_list\d+\.docx`
- `experiment_notes\d+\.md`
- `weekly_progress_.*\.docx`

The convention should be documented in `CLAUDE.md` so future logs go
straight into `notes/archive/<category>/` the moment they are
superseded.

### Special case: one-off analysis scripts

The pattern that prompted "scripts built for a specific quick check
that aren't reproducible" lives under
`experiments/campaigns/<campaign>/`.  Currently the two untracked
files there are:

- `experiments/campaigns/pikp1_avrpia/analysis1_sidechain_clashes.py`
- `experiments/campaigns/pikp1_avrpia/analysis2_misfold_localisation.py`

Both hardcode threshold values, contig-derived fixed-segment
lengths, and run-specific paths.  They run, they answer a question,
and then their context is forgotten.

Proposed convention: `experiments/campaigns/<campaign>/analyses/`
with a sibling `analysis_<N>.yml` per script:

```yaml
# experiments/campaigns/pikp1_avrpia/analyses/analysis1.yml
question: "How often are sidechain clashes in the pose_solved 3a vs 5a runs?"
script: analysis1_sidechain_clashes.py
inputs:
  runs:
    - runs/pose_solved_3a/
    - runs/pose_solved_5a/
fixed_segment_lengths_3a: [6, 27, 2, 24, 6]
fixed_segment_lengths_5a: [6, 27, 2, 24, 6]
rmsd_threshold: 6.0
clash_distance: 2.5
outputs:
  - clash_summary.csv
  - clash_histogram.png
ran_on: 2026-05-19
ran_by: joshua
```

Even a 10-line YAML next to each script restores reproducibility
without forcing the script into the main pipeline.  Convention to be
documented in `experiments/README.md`.

---

## 3. Inventory of every file in `notes/` with proposed action

The table covers everything currently in `notes/`.  "Action" is the
single move recommended; "rationale" is the one-line reason.

### Top level

| File | Status | Action | Rationale |
| --- | --- | --- | --- |
| `notes/SESSION_HANDOFF.md` | TRANSIENT | keep; add frontmatter | per-session brief; users read first |
| `notes/remediation_state.md` | TRANSIENT | keep; add frontmatter | branch + phase state |
| `notes/codebase_remediation_plan.md` | AUTHORITATIVE | keep; add frontmatter | remediation strategy |
| `notes/phase4_architecture_spec.md` | AUTHORITATIVE | keep; add frontmatter | architectural contract |
| `notes/haddock_protocol_experiments.md` | AUTHORITATIVE | move → `subsystems/`; add frontmatter | HADDOCK subsystem reference |
| `notes/design_audit.md` | HISTORICAL | annotate with frontmatter `status: HISTORICAL`; add a note at top: "decisions folded into PIPELINE.md and phase4_architecture_spec.md; do not edit retroactively"; later move → `archive/grill_me/` | 274 KB grill-me record; decisions live elsewhere now |
| `notes/agroinfiltration_coip_handbook.docx` | reference | move → `reference/lab_handbooks/` | wet-lab reference |
| `notes/protoplast_assays_handbook.docx` | reference | move → `reference/lab_handbooks/` | wet-lab reference |

### `notes/inventory/` (Phase-1 snapshot, 2026-05-12)

| File | Status | Action | Rationale |
| --- | --- | --- | --- |
| `01_module_map.md` | HISTORICAL snapshot | move → `inventory_snapshots/2026-05-12/`; refresh post-Phase-4 | per-file LOC scan |
| `02_function_inventory.md` | HISTORICAL snapshot | same | per-file function signatures |
| `03_dependency_graph.md` | HISTORICAL snapshot | same | import + subprocess edges |
| `04_functional_categorization.md` | HISTORICAL (superseded by `PIPELINE.md` §3) | move → `inventory_snapshots/2026-05-12/`; add `superseded_by: notes/PIPELINE.md` | 7-category decomposition |
| `05_findings.md` | HISTORICAL snapshot | move; do entry-by-entry RESOLVED/OPEN/OBSOLETE pass | findings list |
| `06_ubiquitous_language.md` | **AUTHORITATIVE** | **rename/symlink to `notes/glossary.md` at top level** | the single most reusable file in `notes/`; should not be buried |
| `07_ruff_*`, `08_vulture_*`, `09_radon_*` | HISTORICAL | move → `inventory_snapshots/2026-05-12/` | static-analysis dumps |
| `10_phase_1_synthesis.md` | HISTORICAL | move → `archive/phase1/` | Phase-1 close briefing |
| `11_phase_2_plan.md` | HISTORICAL | move → `archive/phase2_planning/` | superseded by `14` |
| `12_phase_2_prompts.md` | HISTORICAL | move → `archive/phase2_planning/` | one-off prompts |
| `13_prompt_3_hpc_tier_tests.md` | HISTORICAL | move → `archive/phase2_planning/` | one-off prompt |
| `14_phase_2_revision_per_module_tests.md` | mixed | extract testing-strategy sections → `operations/testing_strategy.md`; move rest → `archive/phase2_planning/` | testing strategy still authoritative |
| `15_discovery_run_path_coverage.md` | mixed | extract post-rename trail → glossary; move rest → `archive/phase2_planning/` | fixture coverage notes |
| `16_prompt_phase_2_9.md` | HISTORICAL | move → `archive/phase2_planning/` | one-off prompt |
| `17_storage_bloat_audit.md` | AUTHORITATIVE (operational) | move → `operations/storage_audit.md` | open recommendations |
| `audit_pre_hpc_roundtrip.md` | HISTORICAL | move → `archive/phase1/` | pre-HPC audit (Wave 1) |
| `rfdiffusion_parameter_audit.md` | AUTHORITATIVE (subsystem) | move → `subsystems/rfdiffusion_parameters.md`; add frontmatter | RFD param analysis (2026-05-06) |

### `notes/negsteer_notes/`

| File | Status | Action | Rationale |
| --- | --- | --- | --- |
| `negsteer_pipeline_merge_map.md` | HISTORICAL | move → `archive/negsteer_session_logs/` | "two-codebases" merge plan, integration done |
| `notes1.md`–`notes10.md` | HISTORICAL | move → `archive/negsteer_session_logs/` | pre-integration session logs |
| `notes11.md` | mixed | review: if any unique content not yet in `phase4_architecture_spec.md` / `06_ubiquitous_language.md` / `PIPELINE.md`, extract; then move → `archive/negsteer_session_logs/` | latest notes; may have residual canonical content |
| `todo_list2.docx`–`todo_list9.docx` | HISTORICAL | move → `archive/todo_lists/` | superseded by `todo_list10` |
| `todo_list10.docx` | TRANSIENT | keep accessible (top level or `archive/todo_lists/` tagged as latest); audit each item RESOLVED/OPEN | last-known todo list |
| `weekly_progress_apr15_apr22_1.docx` | HISTORICAL | move → `archive/weekly_reports/` | one-off |
| `cluster_submission_guide.md` | unclear | verify currency post-Nextflow integration; if usable, move → `operations/cluster_submission.md`; if obsolete, archive | HPC submission cheat-sheet for standalone negsteer chain |
| `metric_reliability_prompt1.md` | reference | move → `research_memos/` | literature-search prompt |
| `metric_reliability_memo1.docx` | reference | move → `research_memos/` | the produced memo |

### `notes/pipeline_notes/`

| File | Status | Action | Rationale |
| --- | --- | --- | --- |
| `pipeline_notes1.md`–`pipeline_notes14.md` | HISTORICAL | move → `archive/pipeline_session_logs/` | each succeeded by its successor |
| `pipeline_notes15.md` | TRANSIENT | keep at top level until the avrpia campaign closes; then roll into `subsystems/pose_solver.md` and/or `campaigns/pikp1_avrpia_pose_solved.md` | current campaign running log |

---

## 4. One-off / unreproducible script audit

Quick scan of scripts outside `bin/` and `scripts/` that smell ad-hoc.

| Path | Smell | Recommended action |
| --- | --- | --- |
| `experiments/campaigns/pikp1_avrpia/analysis1_sidechain_clashes.py` | hardcoded paths, thresholds, contig fixed-lengths; no CLI; no YAML | move into `experiments/campaigns/pikp1_avrpia/analyses/`; add `analysis1.yml` recording inputs/outputs/thresholds; commit |
| `experiments/campaigns/pikp1_avrpia/analysis2_misfold_localisation.py` | same as above | same as above |
| `experiments/scripts/boltz2_msa_predict.py` | "diagnostic" docstring; campaign-specific design | wrap in a `boltz2_msa_predict.yml` with explicit inputs; or move to `experiments/campaigns/<campaign>/analyses/` |
| `experiments/scripts/build_msa_manifest.py` | same | same |
| `experiments/scripts/compare_msa_vs_nomsa.py` | same | same |
| `experiments/param_derivation/contigs_4a.py` | named after one strategy / one campaign | move to `experiments/campaigns/<campaign>/` and document |
| `experiments/scripts/render_prompt.py` + template files | actually reproducible (config + template) | leave; mark as reproducible in `experiments/README.md` |
| `tests/haddock/pose_solver/pose_solver.py` + `run_example.sh` | CLI + example, reproducible | leave; flag location oddity (pose_solver lives under `tests/` despite being a solver, not a test) |
| `scripts/refresh_orthog_fixture.sh` | committed, parameterised | leave; document in `experiments/README.md` |
| `scripts/rewrite_fixture_paths.py` | committed, reusable | leave |

The most important single change here is the
`experiments/campaigns/<campaign>/analyses/` convention.  Without
it, this category will grow without bound.

---

## 5. Ranked action list

In order of impact-per-effort:

1. **Add `CLAUDE.md` at the repo root** (~30 min).  Highest leverage:
   every future Claude session reads it for free.  No risk; can be
   updated freely.
2. **Add frontmatter to the seven top-level `notes/*.md`** that
   exist today (~15 min).  Makes role machine-detectable.
3. **Add a header to `notes/design_audit.md`** explicitly flagging
   it as a historical record whose decisions now live in
   `PIPELINE.md` and `phase4_architecture_spec.md` (~5 min).  Stops
   future Claude treating it as current.
4. **Promote `inventory/06_ubiquitous_language.md` to `notes/glossary.md`**
   (`git mv` + update cross-references in 3-4 docs, ~15 min).  The
   most cross-referenced file in `notes/` should not be 2 levels deep.
5. **Move `pipeline_notes1.md`-`pipeline_notes14.md`** into
   `notes/archive/pipeline_session_logs/` (~5 min `git mv`).  Keep
   `pipeline_notes15.md` at the current level.
6. **Move `notes/negsteer_notes/notes1.md`-`notes10.md`** into
   `notes/archive/negsteer_session_logs/` (~5 min `git mv`).
7. **Move all `todo_list*.docx`** into `notes/archive/todo_lists/`
   except `todo_list10.docx` (audit pending) (~5 min).
8. **Move static-analysis artefacts** (`07_ruff_*`, `08_vulture_*`,
   `09_radon_*`) into `notes/inventory_snapshots/2026-05-12/`
   (~5 min).
9. **Move `inventory/rfdiffusion_parameter_audit.md`** to
   `notes/subsystems/rfdiffusion_parameters.md` and move
   `haddock_protocol_experiments.md` to `notes/subsystems/`
   (~5 min).
10. **Archive Phase-1/2 planning files** (`10`-`16`,
    `audit_pre_hpc_roundtrip.md`, `negsteer_pipeline_merge_map.md`)
    (~5 min).
11. **Move both `.docx` lab handbooks** to
    `notes/reference/lab_handbooks/` (~3 min).
12. **Set up `experiments/campaigns/pikp1_avrpia/analyses/`** with
    YAML sidecars for the two existing analysis scripts; document
    the convention in `experiments/README.md` (~30 min).
13. **Refresh `inventory/01_module_map.md` + `02_function_inventory.md`
    post-Phase-4** (~60 min: re-run tooling, audit outputs).  Lower
    priority but will diverge increasingly from reality.
14. **Audit `inventory/05_findings.md`** entry-by-entry against
    current code: mark RESOLVED / OPEN / OBSOLETE (~30 min).
15. **Defer:** writing `subsystems/pose_solver.md`.  Wait until the
    avrpia campaign closes; the pose_solver is still iterating.

Items 1-12 together: ~2.5 h of `git mv` + small writes; eliminates
~30 stale files from the top-level `notes/` view; gives every future
Claude session an unambiguous "read this first" path.

---

## 6. What this proposal explicitly does NOT do

- It does not delete anything.  Every move is a `git mv` into
  `notes/archive/`.
- It does not rewrite `notes/design_audit.md`.  The grill-me
  transcript is preserved verbatim as a historical record.
- It does not propose changes to code organisation or naming.  The
  inventory's §10 flagged ambiguities (F1-F11) and the Findings
  document are out of scope for this proposal.
- It does not change `bin/`, `modules/`, `main.nf`, or test fixtures.
- It does not propose any operational change to the pipeline.

Approval needed for: items 1-12 in §5.  Items 13-15 are independent
follow-ups.
