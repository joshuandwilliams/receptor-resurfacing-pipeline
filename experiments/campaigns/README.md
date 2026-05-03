# experiments/campaigns/

One subdirectory per binder-design campaign. A campaign couples a driving
script (parameter generation, pipeline submission, post-hoc analysis) to
its outputs in a single place, so that re-reading a campaign months later
does not require cross-referencing scattered locations.

## Per-campaign layout

```
campaigns/<campaign_name>/
    README.md         # what this campaign is investigating, key params,
                      # status, links to related issues / notes
    <driver>.py       # the campaign script(s)
    results/          # pipeline / analysis outputs (large, often gitignored)
    analyses/         # secondary analyses on top of results/
                      # (notebooks, plot scripts, summary CSVs)
```

`results/` holds raw outputs: cohort CSVs, per-design PDBs, prediction
artefacts, plots. `analyses/` holds secondary work that consumes
`results/` — notebooks comparing tier distributions, custom plots beyond
what the pipeline emits, manual triage CSVs.

## Naming

Campaign directory names should encode receptor and effector:
`<receptor>_<effector>` (e.g. `pikp1_avrpia`). Lowercase,
underscore-separated. Where multiple campaigns target the same pair with
different parameter regimes, append a short discriminator
(`pikp1_avrpia_long_contigs`).

## Outputs and `.gitignore`

Campaign outputs (large PDBs, prediction CIFs, work-tree artefacts) are
governed by the root `.gitignore`. Do NOT add a `.gitignore` inside
`experiments/`. If campaign outputs leak into git, the root `.gitignore`
is the place to fix it.

## Imports

Campaign drivers may import from `experiments/scripts/` and from `bin/`.
The bootstrap that wires `sys.path` for the latter is documented in
`experiments/README.md` (and will land in a future prompt).

## What does not live here

- Generic helpers — those go under `experiments/scripts/`.
- Shared input assets (the receptor PDB, the effector FASTA) — those go
  under `experiments/inputs/`. Campaigns should reference shared inputs
  by relative path, not copy them.
