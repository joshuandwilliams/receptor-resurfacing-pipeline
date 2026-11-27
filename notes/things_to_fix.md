# Things to fix

Surfaced while writing the Chapter 5 methods, 2026-09-04. Nothing here is
actioned. Ordered roughly by how much it changes results.

---

## Changes that would alter results

**AF3 stage to be removed.** `modules/negsteer_af3_nomsa.nf` runs AlphaFold3
natively via `source package`, not containerised, against hardcoded
`/nbi/Reference-Data/AlphaFold/db-v3.0.0` and `db-v2.3.2` paths. It is the only
stage that cannot run off JIC infrastructure. Decision taken to remove it.

**Effector template is silently ignored.** `write_boltz_yaml` emits the
`templates:` block nested inside the chain entry, but Boltz-2 2.2.1 reads
`templates` only from the YAML root (`data/parse/schema.py:1603`). Every
prediction in both this repo and `structure-negative-steering` therefore ran
template-free. Verified in that repo by byte-identical outputs across paired
template and template-free runs. Fixing it means re-running anything being
compared.

**The AlphaFold3 effector template is also inert.** `negsteer_af3_nomsa.nf`
writes the effector template into every `input.json`, but with empty
`queryIndices` and `templateIndices`. AF3 3.0.0 builds template features by
iterating `query_to_hit_mapping.items()` in
`alphafold3/data/templates.py:get_polymer_features`. An empty mapping therefore
leaves the atom positions and masks all zero and the template sequence all
gaps. The effector was untemplated in every AF3 run, the same net outcome as
the Boltz-2 nested-block bug but for a different reason. Fix by emitting the
index mapping, or drop the template and say so.

**Controls run one steering set, real sequences run three.**
`negsteer_controls_n_designs: 1` against `negsteer_n_designs: 3`. The controls
get a third of the attempts, so "the controls failed" is a weaker statement
than it looks. Either match them or state the asymmetry wherever controls are
reported.

**Shape complementarity filter scores bare backbone.** RFdiffusion writes
N, CA, C and O only, with de novo positions labelled glycine. InterfaceAnalyzer
therefore builds side chains for the retained scaffold and scores the design
region on backbone alone. The threshold had to drop from the 0.62
default to 0.3, and across the AVR-PikF campaign Sc ran 0.444 to 0.611, so the
gate never fired. Either drop the gate and keep the value as a diagnostic, or
replace it with something that works on backbones.

**Two RMSD thresholds.** The skip check uses `negsteer_rmsd_threshold` (6.0)
while a correct pose is `negsteer_postprocess_rmsd_threshold` (5.0), so a
sequence can be skipped without qualifying as correctly posed. Worth
collapsing to one parameter unless the split is deliberate.

**Top-N selection ignores the design-region score.** `mpnn_select_top.py:38`
sorts on `mpnn_score`, the whole-receptor value dominated by the fixed
scaffold. `mpnn_design_region_score.py` computes a design-region-confined score
and writes it into the metadata that the selection step receives, but nothing
reads it there. Decide which should rank.

**MMseqs2 clustering coverage.** `-c 0.8 --cov-mode 1` leaves up to a fifth of
the shorter sequence unaligned, so pipeline cluster counts merge sequences that
are not identical. Chapter 5 uses a separate exact-identity analysis instead.

**Interface pLDDT never reaches the output table, and its threshold is dead.**
`run_biophysical_metrics.py` computes `interface_plddt`, but
`merge_orthogonal_metrics.py:132` excludes it from `BIO_COLS` on purpose, and
line 136 notes that `--filter-plddt-min` is "accepted for backward compat". So
`orthogonal_filter_plddt_min: 0.75` raises no flag and gates nothing. The
`rep_interface_plddt_median` column in the table is a different quantity, coming
from the Boltz-2 confidence metrics upstream. Either wire the flag up or drop
the parameter.

**Metrics computed but never reported.** The pipeline pays for several metrics
that no analysis or figure uses. Decide whether to keep, expose or remove them.

| Metric | Where computed | Status |
|:---|:---|:---|
| DockQ, iRMSD, fnat | `compute_interface_metrics.py` | in the table, unused |
| weighted Jaccard, `intact_core` | `compute_interface_metrics.py` | in the table, unused |
| ipSAE at 15 Å (`ipsae_min_15`) | `compute_interface_metrics.py` | in the table, unused. The reported ipSAE is the 10 Å one from `compute_metrics.py` |
| interface hydrogen bonds | `run_biophysical_metrics.py` (MDAnalysis) | in the table, unused |
| pTM | Boltz-2 confidence JSON | in the table, unused |
| interface pLDDT | `run_biophysical_metrics.py` | never reaches the table, see above |

Only shape complementarity, buried surface area and Rosetta ΔG are reported
from the orthogonal stage, alongside pLDDT, ipTM and ipSAE from Boltz-2. The
MDAnalysis dependency exists solely for the unused hydrogen bond count.

## Naming and metadata

**`min_hotspot_frac` does not read the hotspot list.** The gate is
`frac_in >= min_hotspot_frac`, the fraction of all receptor-effector Cα-Cα
contacts within 8 Å falling inside the design region. `hotspot_residues` is
parsed and written to the metrics JSON but never used. Rename.

**Three conflicting version strings.** `manifest.version = '0.2.0'` and the
header comment v0.2.0 in `nextflow.config`, v0.4.0 in the `main.nf` header,
v0.3.0 in `params_example.yml`.

**Stale "polyvaline" comments.** `modules/rosetta_filtering.nf` describes its
inputs as polyvaline backbones. They are backbone-only with glycine at de novo
positions. This wording propagated into a thesis draft before being caught.

## Portability and packaging

**Cluster settings split across two files.** `nextflow.config` holds the
executor, queues, `clusterOptions` and container paths.
`run_pipeline.slurm.sh` separately hardcodes `#SBATCH -p jic-medium`,
`PIPELINE_DIR=/hpc-home/jowillia/...` and a mail address.

**No nf-core scaffolding.** No `nextflow_schema.json`, `.nf-core.yml` or
`CITATIONS.md`. `nf-core pipelines create` gives the structure, lint and CI
without any obligation to submit. Full membership is a much higher bar, and not
worth pursuing on current evidence. It needs MIT, a `TEMPLATE` branch, two
reviewers including one from core or maintainers, and a test profile that runs
end to end on trivial data. That last requirement is hard for GPU stages with
large weights. Rosetta
is also free for academics only, which sits awkwardly against nf-core's
open-to-all stance.

**No `environment.yml`.** Analyses run across two environments, the
`receptor-analysis` Jupyter kernel (Python 3.11.15) and the miniforge base
install (Python 3.12.10) that the standalone scripts use. Neither is pinned in
the repo.

## Analyses

**`analyses/05-steerability-correlates/05-steerability-correlates.qmd` is
superseded.** It is the LightGBM and SHAP version, defining a hit as
`cross_tier` in A/B/C across all 128 designs (38 hits). `steerability_features.py`
replaced it with a per-feature Mann-Whitney over the 106 steered designs
(18 corrected). Retire or annotate the old document.

**`analyses/10-gblock-selection` covers the d26 variant campaign, not the
original ten.** The selection of the ten agroinfiltrated designs is described
in the thesis but is not scripted anywhere in the repo.
