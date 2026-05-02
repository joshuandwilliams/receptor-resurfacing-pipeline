# Prompt: Metric reliability analysis for protein structure prediction ranking

## Your role

You are an expert computational structural biologist advising on a protein complex structure prediction pipeline. I need you to conduct a thorough literature search and reasoning-based analysis to determine which metrics I should use to rank predicted protein–protein binding poses, specifically in the context of de novo designed binder interfaces predicted by Boltz2.

The bottom line: I want a ranking metric whose values correspond to real structural accuracy. I don't care if the metric values are low — I'd rather have an honest metric that says "this is uncertain" than one that confidently says "this is great" when the underlying prediction is wrong.

## Context: what the pipeline does

I have an **RFDiffusion-designed protein binder** with a **ProteinMPNN-designed sequence**. When I predict the complex structure with Boltz2, the binding partner (effector) docks to the **wrong surface** of the receptor (~30 Å from the designed binding site). My pipeline introduces 3 mild point mutations on the wrong-interface surface to disrupt it, re-predicts with Boltz2, and checks whether the effector relocates to the designed binding site. This approach is called "negative steering" — it's a computational technique for improving structure prediction accuracy on designed protein complexes.

The pipeline is effective — 57/100 steered designs find the correct interface (receptor-aligned effector RMSD < 5 Å vs the design model). But most of those have contaminating mutations (the steering mutations themselves contact the effector in the new pose). Since the final protein will be expressed from the original ProteinMPNN sequence without the steering mutations, only designs where the mutations DON'T contact the effector are trustworthy structural predictions.

## The specific problem I need you to solve

I currently rank surviving designs by **ipSAE** (interface-predicted Symmetrised Aligned Error), which is derived from Boltz2's PAE matrix. My data shows this metric is unreliable for distinguishing correct from incorrect interface placement:

| Sequence | ra_eff (Å) | ipSAE | iPTM | Interface placement |
|---|---|---|---|---|
| WT sequence (100 seeds, ds=5) | 29.8 ± 0.7 | 0.787 ± 0.055 | 0.893 ± 0.020 | WRONG (native-like surface) |
| design_24 (10 seeds, validated) | 2.51 ± 0.53 | 0.755 ± 0.055 | 0.862 ± 0.019 | CORRECT (designed surface) |
| design_12_s0 (partially reverted) | 2.62 | 0.522 | 0.754 | CORRECT |
| design_10_s2 (partially reverted) | 3.56 | 0.245 | 0.649 | CORRECT |

**Boltz2 assigns higher confidence to the wrong interface than the correct one.** The wild-type sequence confidently places the effector on the wrong surface (ipSAE 0.787), while design_24 correctly places it on the designed surface with lower confidence (ipSAE 0.755). Designs with a single remaining steering mutation that correctly find the designed interface have even lower confidence (ipSAE 0.25–0.52).

This is a fundamental calibration problem: the wrong interface is a naturally occurring protein surface that Boltz2 has seen in its training data, while the designed interface is novel.

## Data I'm providing

I'm providing the following datasets from my pipeline runs. Please study them carefully.

1. **per_seed_metrics.csv** — 10-seed validation of design_24 (the best clean steered design, carrying mutations K5E M22E L36D). Every seed finds the correct interface. Shows within-sequence variance of all metrics.

2. **variance_summary.csv** — Aggregate statistics from the 10-seed validation.

3. **steered_results.csv** — Full steered results from the num_seeds=3 run (v4). 30 unique sequences × 3 seeds = 90 predictions. Contains structural RMSDs, Jaccard overlaps with true/wrong interfaces, sequence_group, seed_index.

4. **steered_results_aggregate.csv** — Per-sequence-group aggregates from v4.

5. **all_results_multicycle_with_metrics.csv** — The full pipeline output from v4 including confidence metrics (ipSAE, iPTM, pLDDT, PAE), contamination analysis, and reversion verdicts.

6. **passing_summary.csv** — The ranked survivors from v4.

## What I need you to do

### Part 1: Literature review

Search the literature for:

1. **ipSAE / iPTM / pTM reliability on de novo designed interfaces.** Specifically:
   - Do the original Boltz2, AlphaFold2, AlphaFold-Multimer, and RFDiffusion papers discuss confidence metric reliability for designed (non-native) interfaces?
   - Has anyone shown that ipSAE or iPTM correlates with experimentally determined binding for designed protein complexes?
   - Is there literature on confidence metrics being systematically higher for wrong-but-familiar interfaces vs correct-but-novel interfaces? This is a known concern in the protein design field.

2. **DockQ, iRMSD, lDDT-PPI, fnat as structural quality metrics.** These are reference-dependent structural metrics:
   - Which is most predictive of actual complex formation based on published benchmarks?
   - Which handles disordered termini and flexible loops best?
   - Which is standard in the protein design and docking communities?

3. **AlphaFold2/Boltz2 confidence calibration.** Is it documented that structure prediction confidence metrics are calibrated on native complexes and may be miscalibrated on computationally designed ones?

4. **What do successful de novo binder design campaigns actually use for computational ranking?** What metrics predicted experimental success in published design campaigns (e.g. from the Baker lab, Rosetta-based campaigns, or recent diffusion-model-based design work)?

### Part 2: Analysis of my data

Using the datasets I've provided:

1. **Within-sequence metric stability.** For designs predicted 3 times (v4), how stable is each metric across seeds? Compute coefficient of variation for ra_eff, ipSAE, iPTM, pLDDT, PAE, true_jaccard, wrong_jaccard. Which metrics are reproducible enough to be useful for ranking?

2. **Metric correlation with correct interface placement.** Across all designs, does ipSAE / iPTM / pLDDT correlate with ra_eff? Or are they orthogonal? Describe the relationship in detail.

3. **Can confidence metrics distinguish correct from wrong interface?** Using the v4 data where we have both correct-interface and wrong-interface predictions from the same experimental conditions, test whether any confidence metric has discriminative power. Consider the separation between distributions.

4. **What about within the correct-interface regime only?** For designs that DO find the correct interface (ra_eff < 5 Å), does higher ipSAE correlate with lower ra_eff or higher true_jaccard? In other words, once I've already identified correct-interface designs, does ipSAE help rank them by quality?

5. **Composite metrics.** Would a combination of metrics (e.g. ra_eff × true_jaccard, or a weighted sum) be more robust than any single metric? What about metrics that don't require a reference structure (for future deployment where the design model may be the only reference)?

### Part 3: Recommendation

Based on the literature and data analysis, provide:

1. **A specific ranking metric or formula** that I should use to rank surviving designs. Justify it with both literature precedent and my data.

2. **Which metrics should be used as pass/fail filters** (not for ranking, just for screening out failed predictions). Specify thresholds.

3. **What to do when a high-quality reference structure is not available.** My pipeline requires a design-model reference for ra_eff and true_jaccard. For future deployment on designs where only the RFDiffusion backbone is available, what metric(s) should I use?

4. **Whether ipSAE should be retained, replaced, or repurposed.** Be specific about when (if ever) ipSAE is informative in a protein design context.

5. **Any metrics I'm not currently computing that I should add.** For example, interface buried surface area, shape complementarity, per-residue confidence at interface positions, predicted ΔΔG, hydrogen bond counts, etc.

## Important technical constraints

- The pipeline predicts with **Boltz2**, not AlphaFold2. Boltz2's native confidence outputs are PAE, pTM, iPTM, pLDDT. ipSAE is derived from PAE post-hoc by our pipeline code.
- The **reference structure is an RFDiffusion Cα-only backbone**, not a solved crystal structure. There are no side chains in the reference, only Cα positions. This limits which reference-dependent metrics can be computed.
- The **receptor sequence is ProteinMPNN-designed**, so it is a computationally designed sequence, not a natural one. This matters because confidence metrics trained on natural proteins may not transfer to designed sequences.
- The biological system is a **plant immune receptor (NLR HMA domain)** binding a recognition target. The designed binder is an RFDiffusion-designed variant of this HMA domain. This is standard plant molecular biology research on receptor–ligand recognition.
- I'm ultimately interested in **whether the designed complex will form as predicted** — i.e. whether the computational structure prediction is accurate enough to guide experimental construct selection.

## Format of your response

Please structure your response as:
1. Literature findings (with citations and DOIs where possible)
2. Data analysis results (with specific numbers from my datasets)
3. Concrete recommendations with justification
4. A summary table of "metric → use as ranking / filter / don't use / compute but don't act on"

Be thorough and detailed. I'd rather have a long, careful analysis than a quick answer. This decision determines which protein constructs get selected for experimental characterisation — incorrect ranking wastes significant experimental resources.
