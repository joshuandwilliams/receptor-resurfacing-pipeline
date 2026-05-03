# experiments/

Tooling and working space for designing binder campaigns.

## What this is for

The Banfield lab studies plant immune receptors with the goal of developing wheat varieties that recognise novel and changing pathogens. The production pipeline (top-level `main.nf` and the modules under `modules/`) takes a receptor–effector pair and designs receptor variants predicted to bind the target. This `experiments/` directory is where those design efforts are planned, run, and analysed.

A campaign produces a list of designed receptor variants for one receptor–effector pair. Promising candidates from a campaign are taken forward to wet-lab validation by transient expression in *N. benthamiana* via agroinfiltration. The wet-lab work happens outside this repo; `experiments/` ends at the candidate-selection step.

This directory exists separately from the production pipeline because campaign design is iterative and exploratory in a way the pipeline itself is not. The pipeline is a fixed transformation from inputs to outputs; campaigns are the work of choosing those inputs well, running the pipeline against them, and making sense of what comes back. Keeping that work in its own directory keeps the production code clean and gives campaigns a durable home.

Nothing in `experiments/` feeds back into the production pipeline as code, fixtures, or defaults. The relationship is one-way: `experiments/` consumes the pipeline.

## Directory layout

```
experiments/
  scripts/           Reusable analysis and plotting helpers, importable
                     from any campaign script. Add helpers here when the
                     same code is needed across two or more campaigns.

  param_derivation/  Code for deriving pipeline input parameters from
                     structural inputs — most substantially the
                     RFDiffusion contig strings, which encode which
                     binder regions are anchored and which are
                     designable. Separated from scripts/ because contig
                     derivation is a substantial concern in its own
                     right with multiple strategies (see lifecycle
                     Stage 4).

  inputs/            Shared structural and sequence assets used across
                     campaigns: PDB structures (experimental and
                     AF3-predicted), UniProt sequences, reference
                     materials. Anything reusable across campaigns
                     targeting the same protein lives here, not inside
                     a campaign directory.

  campaigns/         One subdirectory per campaign. Each campaign is
                     one receptor–effector pair and may contain
                     multiple runs against that pair (different
                     parameter choices, different contig strategies).

    <campaign>/
      README.md      What this campaign is asking, and notes on its
                     overall trajectory.
      inputs/        Campaign-specific inputs reusable across this
                     campaign's runs: the assembled complex, curated
                     residue lists for the contig strategies, etc.
      runs/          One subdirectory per pipeline run. Each run is
                     one invocation of the pipeline with one parameter
                     set. See conventions section below for run naming.
      analyses/      Cross-run analyses for this campaign — comparisons
                     across runs, hypothesis-testing across the
                     campaign's full trajectory.
```

Each campaign run produces both `results/` (raw pipeline outputs) and `analyses/` (derived plots and tables specific to that run) inside its own subdirectory under `runs/`. Cross-run work — comparing one run's behaviour to another's — lives in the campaign-level `analyses/` directory, not in any individual run.

## Relationship to the rest of the repo

Campaign scripts can import from two places:

- `experiments/scripts/` — helpers shared across campaigns. Imported as `from experiments.scripts.<module> import …`.
- `bin/` — the production pipeline's Python modules, which contain domain logic worth reusing (contact-residue calculations, cohort loaders, tier classifiers, etc.). The import pattern needed to make `bin/` reachable from campaign scripts is documented separately — see the conventions section below.

Reusing `bin/` rather than copying logic into `experiments/` matters. The production pipeline already has a documented duplication problem (see `notes/inventory/05_findings.md`); adding a fifth contig parser or a fifth contact-residue function inside `experiments/` would make that worse.

## Mac and HPC

Some campaign work runs locally on Mac (analysis of existing pipeline outputs, plotting, parameter derivation from structures); some runs on HPC (the pipeline itself, AF3 structure prediction). The Mac side is authoritative; HPC is treated as a compute target. The repo-wide `scripts/sync_to_hpc.sh` handles the round-trip — see its source for which paths are included and excluded.

Outputs of pipeline runs (`runs/<name>/results/`) are produced on HPC and pulled back to Mac for analysis. They are gitignored (`experiments/campaigns/*/runs/*/results/` and the corresponding `analyses/` and `work/` patterns) — large prediction artefacts don't belong in git.

## Campaign lifecycle

A campaign goes through six stages, from choosing a target to selecting candidates for wet-lab validation. The first two produce assets that are reusable across any campaign targeting the same proteins (and so live in `inputs/`); the rest are campaign-specific.

### 1. Resource gathering

Identify what structural and sequence data already exists for the chosen receptor and effector. PDB entries, UniProt records, published structures of related complexes. The output of this stage is a sense of what the campaign has to work with and what it will need to generate.

### 2. Structure acquisition

Secure 3D structures for both the receptor and the effector. Where an experimental structure exists, use it. Where it doesn't, predict it with AlphaFold3 on HPC. Predicted structures live in `inputs/` alongside experimental ones, with provenance recorded so future campaigns know which is which and where each came from.

### 3. Complex assembly

Position the receptor and effector together so the intended design interface is geometrically correct. Currently a manual ChimeraX step: align the two structures using homologous solved complexes as templates. HADDOCK is a possible future automation but is not the current workflow.

The assembled complex is a hand-curated artefact, not a pipeline output. Capture provenance in the campaign README — which homologous structures were used as templates, which ChimeraX session produced the alignment — otherwise the campaign isn't reproducible.

### 4. Contig generation

Choose which binder regions are anchored (kept fixed) and which are designable (handed to RFDiffusion). The contig string passed to the pipeline is the encoding of that choice.

Four strategies of increasing flexibility, ordered from most constrained to most:

- **4a. Interface-facing side chains.** Automatic, geometric. Anchor any binder residue whose side chain points into the interface within a distance cutoff. On regular secondary structure this naturally produces a periodic pattern (e.g. every second or third residue on an α-helix).
- **4b. Biologically variable contact residues.** Curated. Anchor residues known from the literature to be important for the canonical binding interaction *and* known to be variable across natural sequence diversity. Identifying them requires reading the field's existing characterisation of the receptor–effector pair.
- **4c. Binder residues near specified target residues.** Mixed. Specify a set of effector residues of interest, then anchor any binder residue within a distance cutoff of those target residues. The anchor set is defined by the *target* side of the interface rather than the binder side.
- **4d. Whole-interface redesign.** Automatic, geometric. Anchor everything *not* on the target-facing face of the binder; the whole interfacial surface is in scope for redesign.

Which strategy works best is itself an open empirical question. The first rounds of agroinfiltration are intended to investigate this — running campaigns under different strategies and seeing which produces designs that behave well in the wet lab. Until that's resolved, no strategy is the default; the choice is part of what each campaign is testing.

**Gap-filling.** Strategies 4a, 4b, and 4c can produce anchor sets so dense that RFDiffusion has no contiguous designable stretches to work with. The fix is to expand the designable region around each selected anchor, converting a sparse pinned set into longer designable segments.

**Constitutive activity — interpretive caveat, not a hard exclusion.** The binders worked with here are often constitutively active (firing without the effector present) when inward-pointing or conserved residues are changed — those changes destabilise the autoinhibited resting state. Strategy 4d is most exposed to this risk because its broad selection sweeps in conserved residues by default; 4a sidesteps inward-facing residues by construction.

This is not encoded as an exclusion mask handed to RFDiffusion. It's a consideration applied at candidate selection (Stage 6) and during downstream interpretation: designs with many changes to internally-facing residues are flagged as more likely to be autoactive. Specific residues with known autoactivity effects exist in the literature, but their relevance to a designed variant is uncertain — RFDiffusion may have changed the surrounding region enough that the known effect doesn't reproduce in context. So the caveat is about awareness during interpretation, not a constraint on the design space.

### 5. Pipeline execution

Choose pipeline parameters and run on HPC. The principal heuristic for parameter choice is **breadth beats depth**: across the pipeline's stochastic stages, success rate per design is improved more by exploring more starting points than by exhausting any individual one.

Concretely:

- More RFDiffusion designs.
- Fewer MPNN sequences per design.
- More sequences entering negative steering.
- Smaller negative-steering search per sequence.

The same compute budget gets distributed across more independent attempts rather than concentrated on fewer thorough ones. The pipeline's bottleneck is finding *any* viable design, not refining a marginal one — expected yield is higher when bets are spread.

A campaign may have multiple runs in `runs/` because parameter choice is a knob worth turning across attempts. A second run with a different contig strategy, or a different breadth/depth ratio, is a continuation of the same campaign rather than a new one.

### 6. Candidate selection

Pick a list of designed receptor variants to take forward to agroinfiltration. **This is a judgement step, not a ranking-cutoff step.** The pipeline produces a ranked output, but the top N by composite score is not the right candidate list.

The pipeline's output metrics (Boltz confidence, jaccard scores, Rosetta ΔΔG, AF3 cross-check agreement, and others) measure different aspects of design quality. Which of those aspects actually predicts agroinfiltration success is unknown. The right strategy is to deliberately pick a *diverse* set of candidates spanning different metric profiles, treating each batch of agroinfiltrations as an opportunity to learn which metrics matter.

Selection at this stage is also where the constitutive-activity caveat from Stage 4 is applied: designs with many changes to internally-facing residues are flagged as more likely to be autoactive, and that flag is one input to whether a candidate is taken forward.

The wet lab is part of the experimental loop here: each round tests both the candidates *and*, implicitly, the hypotheses about what makes a candidate good. Over time this should converge on better selection criteria — but only if selection is varied enough each round to discriminate between competing hypotheses. The first such hypothesis on the queue is which contig strategy from Stage 4 produces the best wet-lab outcomes.

The selection rationale (which candidates were picked and why, in terms of metric profile) is itself a campaign artefact worth preserving. When a candidate succeeds or fails in agroinfiltration, you want to be able to look back and check what its profile was. This is the kind of analysis that lives in the campaign-level `analyses/` directory.

Tooling for this stage is mostly visualisation and exploration, not automation. The point of Stage 6 is human judgement informed by hypotheses, so the right helpers are ones that *display* candidate metrics in ways that make trade-offs visible (parallel coordinates plots, metric-vs-metric scatter plots, tier breakdowns), rather than ones that pick candidates automatically.

## Conventions

This section codifies the naming and placement decisions that have settled. They aren't load-bearing — nothing breaks if a campaign drifts from them — but consistency makes campaigns easier to navigate and compare, especially when looking back at older work.

### Campaign naming

Campaign directories are named `<receptor>_<effector>`, lowercase, underscore-separated. Examples: `pikp1_avrpia`, `pikp1_avrpikf`, `pikp1_pby2`. The receptor comes first because campaigns are usually grouped mentally by receptor (one receptor against several effectors).

Once a receptor–effector pair has a campaign directory, all subsequent work against that pair lives inside it as additional runs. Starting a new campaign directory for an already-targeted pair would only happen after a substantial rethink that makes the old campaign's runs no longer comparable to the new ones — and even then, it's a strong signal worth pausing on.

### Run naming

Each subdirectory under `<campaign>/runs/` is named `v<N>_<slug>`. The version number `v<N>` is monotonic within the campaign — `v1`, `v2`, `v3`, never reused. The slug is a short descriptor of what makes this run distinct from the others, typically the contig strategy and any other notable parameter choice.

Examples:

```
runs/v1_4a/
runs/v2_4b/
runs/v3_4a_more_designs/
runs/v4_4c_tightened_cutoff/
```

The version number is the canonical reference; the slug is a hint for humans skimming the directory. Dates aren't in the run name because each run's outputs are timestamped internally — looking at the run is enough to recover when it was produced.

### Input placement

All inputs to a campaign — the assembled complex from Stage 3, any curated residue lists for Stage 4 strategies, the binder's exclusion-relevant residue annotations, anything else — live in `<campaign>/inputs/`. They are placed there once and reused across runs of that campaign. Input filenames should reflect the campaign name and the strategy or purpose they serve, e.g. `pikp1_avrpia_complex.pdb`, `pikp1_avrpia_4b_residues.txt`, `pikp1_avrpia_4c_target_residues.txt`.

The single exception is `params.yml`. Each run's `params.yml` is the parameter file actually passed to the pipeline for that run, and it lives at the run root (`<campaign>/runs/<run>/params.yml`) rather than in `<campaign>/inputs/`. It's run-specific by definition — different runs exist precisely because they have different `params.yml` files — so it travels with the run, not with the campaign.

Shared assets that are reusable across multiple campaigns (e.g. an experimental or AF3-predicted structure of a receptor used in several campaigns) live in `experiments/inputs/` at the top level, not inside any campaign directory.

### File naming inside `experiments/inputs/`

Top-level shared inputs need enough provenance encoded in their filenames that future-you can tell what they are without opening them. The conventions:

- **Structures:** `<protein>_<source>.pdb`, where `<source>` is `experimental` for PDB-derived structures or `af3` for AlphaFold3-predicted ones. Examples: `pikp1_experimental.pdb`, `avrpia_af3.pdb`.
- **Sequences:** `<protein>_uniprot.fasta` for UniProt-sourced sequences, with the UniProt accession recorded inside the file or in a sibling `.notes` file.
- **Provenance notes:** any input where the filename can't carry full provenance (which PDB entry, which UniProt version, which AF3 run produced it) gets a sibling `.notes` plain-text file with the same stem — e.g. `pikp1_experimental.pdb` and `pikp1_experimental.notes`. The `.notes` file is short: where it came from, when, and any caveats.

These conventions are deliberately loose. The point is provenance recoverability, not strict schema compliance.

### Mac and HPC split

The repo lives on both the Mac (authoritative) and the HPC (compute target). Different parts of the campaign workflow happen on each:

- **Mac:** resource gathering (Stage 1), complex assembly in ChimeraX (Stage 3), contig generation (Stage 4), candidate selection and analysis (Stage 6), all README and provenance writing.
- **HPC:** AlphaFold3 structure prediction (Stage 2), the pipeline itself (Stage 5), and any other GPU-bound or large-compute work.

Everything in `experiments/` other than pipeline outputs is written on the Mac and synced to the HPC. Pipeline outputs (`runs/<run>/results/`) are produced on the HPC and pulled back to the Mac. The repo-wide `scripts/sync_to_hpc.sh` handles both directions; consult its source for which paths it includes and excludes.

When a campaign is in flight, the working assumption is that the Mac copy is the source of truth for everything except the most recent run's `results/` directory, which may exist on the HPC before it has been pulled back.

### Importing from `bin/`

> **TODO:** the helper module that makes `bin/` importable from campaign scripts has not been written yet. When it lands, this section will document the import recipe (likely a small `experiments/_path_setup.py` module imported at the top of each campaign script that adds the repo's `bin/` directory to `sys.path`). Until then, campaigns that need to reuse `bin/` logic will need to handle path setup themselves.

Once the helper exists, the rule will be: never copy logic from `bin/` into `experiments/`. The production pipeline already has a documented duplication problem; adding more copies inside `experiments/` would compound it.
