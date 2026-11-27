# Pipeline Notes 6 — RFDiffusion Design Strategies, Classifier Flagging, and Scoring Methodology

## Context

This session defined design region strategies for RFDiffusion-based receptor resurfacing. The goal is to avoid constitutive signalling (autoactivity) in plant NLR immune receptors. The discussion covered the technical methodology. It also covered a practical problem, namely content classifiers flagging conversations about computational protein design in plant pathology.

## Classifier Flagging Analysis

### Observable Patterns

**Individual phrases that don't flag:**
- "Nextflow pipeline for resurfacing" — clean
- "NLR receptor" — clean  
- "I want to avoid autoactivity in the plant" — clean
- "Nextflow pipeline for resurfacing an NLR plant immune receptor to recognise fungal effectors" — clean (with plant/fungal context)

**Phrases that trigger flagging when combined:**
- "Nextflow pipeline for resurfacing" + "NLR receptor" — flagged when paired
- "I want to avoid autoactivity in the plant" + "I want to make sure RFDiffusion has enough flexibility to work with" — flagged when combined

**Project context effect:**
Flagging appears more frequent once the user is added to a project holding extensive technical notes. Those notes cover negative steering, contamination detection, effectors, receptors, and computational design workflows. The classifier likely responds to the denser domain-specific vocabulary that project context loads.

### Mitigation Strategies

**Per-message disambiguation:**
- Include "plant immune" with "receptor" 
- Include "fungal" or "oomycete" with "effector"
- Use "avoid constitutive signalling" instead of "autoactivity"
- Use "conformational search space" instead of "flexibility to work with"
- Frontload context: "I'm a plant pathology researcher at the John Innes Centre working on computational design of disease-resistance receptors"

**Consistent feedback:**
Use the thumbs-down button on flagged responses. It signals the classifier tuning team directly, with conversation context attached.

## RFDiffusion Design Region Strategies

### Current Approach
User is currently avoiding:
- Residues whose side chains point inwards (based on native structure)
- Residues not on the natural interface

This pairs with the negative steering validation workflow. It exploits the bias of structure prediction models toward canonical binding sites.

### Four Strategies Assessed

#### Strategy 1: Alternating Outward-Facing Residues Only
**Description:** Only diffuse every other residue across interface regions where sidechains point outward.

**Subject specialist perspective:** Most biologically conservative for preserving autoinhibition. It risks surfaces that are geometrically plausible but biologically blind to specificity-determining features. Sidechain orientation filters built on a single conformational state miss the allosteric rearrangements of activation.

**Computational perspective:** Under-utilizes RFDiffusion's strengths. Diffusing alternating residues with frozen intervening positions gives little backbone flexibility. It amounts to sidechain-level redesign with extra steps. ProteinMPNN alone suits this scope better.

**Reviewer perspective:** Would question why RFDiffusion is used at all under such constraints. Would expect benchmarking against sidechain-only redesign to show what backbone diffusion contributes.

#### Strategy 2: Known-Important + Variable Residues + Neighbors  
**Description:** Target residues known to be important for binding in homologues, highly variable positions, plus surrounding residues.

**Subject specialist perspective:** Most defensible for Pik-style scaffolds with mapped allelic variation. Uses evolutionary signal to find recognition plasticity. Limited to scaffolds with resolved allelic data, but the highest-signal option within that scope.

**Computational perspective:** Technically tractable, but RFDiffusion works best on contiguous diffusable segments (4-6 residues minimum). May need to widen "surrounding residues" to get continuity. Strongest case for negative steering validation.

**Reviewer perspective:** Would probe whether "important/variable" residue selection is independent of the validation test set. Information leaks if the same allelic series both guides design and tests recognition.

#### Strategy 3: Whole-Interface Redesign
**Description:** Redesign entire binding interface while maintaining intact backbone scaffold for topology constraint.

**Subject specialist perspective:** Highest autoactivity risk, though the risk varies by scaffold. Relatively safe for integrated HMA domains (Pik-1), where the recognition module is decoupled from NBD regulatory machinery. Dangerous for canonical sensor NLRs, where the LRR does double duty (recognition + NBD restraint).

**Computational perspective:** Optimal use of RFDiffusion's capabilities. Canonical binder design recipe with fixed scaffold + effector conditioning. Requires high throughput because success rates are low (single-digit percent). Needs a downstream autoactivity filter. That filter predicts the full-length receptor and compares its NBD-LRR interface to the native autoinhibited state.

**Reviewer perspective:** Would demand wet-lab autoactivity data. The computational story must end with falsifiable predictions tested in planta.

#### Strategy 4: Hotspot-Driven Shell Design (Proposed)
**Description:** Specify hotspot residues on effector, let RFDiffusion choose receptor residues to redesign within bounded shell (e.g., 8 Å radius from effector hotspots).

**Advantages:**
- Doesn't pre-commit to native binding mode, so it can discover non-native recognition geometries
- Geometrically bounded redesign region keeps the diffusable area contiguous
- Autoactivity risk tunable via shell radius (6 Å conservative, 10 Å aggressive)  
- Integrates naturally with negative steering, since redesigned residues that predictors push toward canonical sites become negative examples

**Requirements:** Initial effector placement (from RFDiffusion, docking, or AF3 prediction of homologous complex)

## Current Task Context: Reproducing Rational Engineering

### Problem Setup
The current validation target is binding demonstrated by rational engineering, a single residue change. This creates specific challenges:

**Structure predictor limitations:** AF3/Boltz-2 cannot resolve the functional effect of a single residue change, especially with limited plant-specific training data. The models predict near-identical structures for wild-type and functionally different point mutants.

**Validation paradox:** Telling RFDiffusion exactly which residue to change hands it the answer. Successful "recovery" is then uninformative as a pipeline test.

### Reframed Evaluation Question
**Not:** "Can we recover the rational design's exact solution?"  
**But:** "Do our designs predict as well as the validated positive control?"

This is legitimate for RFDiffusion inpainting, which requires:
- Correctly oriented receptor-effector complexes
- Binding face specification  
- Relative orientation fixed

These are standard input requirements for the method, not "giving away the answer."

## Scoring Methodology

### Primary Geometric Metrics
**Correct placement score:** Does Boltz-2 place the effector where RFDiffusion intended?  
**Weighted Jaccard contact overlap:** Predicted interface contacts vs. intended contacts. Uses smooth distance weighting (e.g., `exp(-(d-4)²/2.25)`) instead of hard cutoffs, for robustness to predictor noise.

**DockQ vs. experimental structure:** Most biologically grounded metric. It compares the predicted complex against the experimentally solved rational design structure. That is ground truth validation, not self-consistency.

### Secondary Orthogonal Metrics
**Rosetta interface energy (ddG):** Most informative energy metric, with caveats about minimization sensitivity and hydrogen bond bias.  
**Buried surface area (BSA):** Sanity check filter (too small = probably not binder) rather than ranker.  
**Interface packing density:** Complements ddG by capturing shape complementarity independent of polar contacts.

### Demoted Confidence Metrics
**ipSAE/ipTM:** Relegated to "plus factor" rather than primary ranker due to:
- Wild fluctuation between samples of identical sequences
- Dramatic changes when off-interface residues move during negative steering
- No consistency as affinity predictors

### Composite Ranking Strategy
**Tier-then-composite method:** DockQ vs. experimental structure is the dominant term, giving ground truth validation. Geometric consistency metrics support it (placement + weighted Jaccard). Orthogonal energy metrics act as secondary filters.

## Key Insights

### On Predictor Capabilities
Structure predictors excel at geometry prediction (what the PDB teaches) but fail at affinity discrimination (what the PDB doesn't teach). Leverage the strength, acknowledge the limitation.

### On Validation Design
For single-residue rational designs, geometric primary metrics with energetic checks validate more defensibly than confidence-based ranking. Comparison against the experimental structure gives genuine ground truth, not a self-consistency check.

### On Autoactivity Risk Management
Risk varies dramatically by scaffold architecture. Integrated HMA domains (Pik-1) allow aggressive interface redesign. Canonical sensor NLRs (Sr35/ZAR1-like), where the LRR restrains the NBD, require conservative approaches. A computational autoactivity filter can predict the full-length receptor and compare NBD-LRR interface geometry to the native autoinhibited state.

### On Negative Steering Integration
Predictor bias toward canonical sites becomes a feature, not a bug. Redesigned positions that predictors revert to canonical interfaces supply the negative steering examples the validation workflow needs.

## Next Steps

### Recommended Strategy Implementation
1. **Primary:** Strategy 4 (hotspot-driven shell) with 8 Å radius from effector hotspots
2. **Secondary:** Strategy 2 (known-important positions) if sufficient allelic data exists
3. **Validation:** Tier-then-composite ranking with DockQ vs. experimental structure as dominant metric
4. **Autoactivity filter:** Full-length receptor NBD-LRR interface RMSD vs. native autoinhibited state

### Pipeline Development Priorities
1. Implement weighted Jaccard contact overlap (smooth distance weighting)
2. Integrate experimental structure DockQ calculation
3. Validate composite scoring against rational design positive control
4. Develop full-length receptor autoactivity filter for high-priority designs

### Wet-lab Validation
Progress from geometric+energetic computational filters to plant cell pack (PCP) throughput triage. Agroinfiltration of top candidates then gives final validation.
