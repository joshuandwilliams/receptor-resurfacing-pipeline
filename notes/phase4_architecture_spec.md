# Phase 4 Architecture Spec

The written specification for the Phase 4 codebase rebuild.  Strategy B (deep modules) — the whole codebase is organised around a small set of domain types, each carrying both its data and the behaviour around it.  Naming convention: domain nouns; descriptive multi-word names are fine when they reduce ambiguity.

This document grows as the Phase 4 grill-me works through it.  Final form is the architectural contract: implementation in the migration phases follows this spec exactly.  Anything not in the spec is internal detail of the type it lives in.

Live document; updated as decisions land.  Cross-references to `notes/design_audit.md` Session 6+ for the reasoning behind each decision.

---

## Step 1 — Candidate types (locked)

A scan of the codebase for nouns that recur across modules.  The criterion: a concept becomes a deep module if it carries multiple related fields AND has behaviour around it AND appears in more than one pipeline stage.

The list below is the agreed candidate set after Session 6 grill-me.  Interfaces for each type are defined in Step 2.

### The types

1. **`ProteinStructurePrediction`** — one prediction's structural payload.  Wraps a PDB (Boltz output) or a CIF (AF3 output).  Knows: chain sequences, atom coordinates, contact residues, RMSDs against another prediction.  Does NOT know predictor-specific confidence metrics.

2. **`BoltzConfidenceMetrics`** — the confidence metrics produced by one Boltz prediction (pLDDT, PAE, ipTM, iPSAE, pae_pass_frac, interface_plddt, …).  Attaches to a `ProteinStructurePrediction`.

3. **`AF3ConfidenceAggregate`** — the aggregated confidence metrics produced by ONE AF3 call (which itself runs seed × diffusion-sample = many predictions internally).  Exposes `best_ra_eff`, `mean_ra_eff`, `best_iptm`, `mean_iptm`, `n_correct_interface`, etc.  Attaches to a `ProteinStructurePrediction`.

4. **`DesignedBackbone`** — one RFDiffusion output: the receptor backbone PDB (no sequence yet), its ContigSpec, its design_id (e.g. `design_03`), its Rosetta-filter metrics.  Knows: design region positions, whether it passes the Sc filter.

5. **`DesignedSequence`** — one MPNN-designed sequence overlaid on a `DesignedBackbone`.  Carries the corrected receptor amino-acid string, the slice at design positions, the MPNN score.  Identified by `design_<N>_seq_<M>`.

6. **`PositionSet`** — a set of residue positions on the receptor, *with frame*.  Frame is one of `native`, `designed`, `prediction`.  Operations between sets in different frames raise; conversion is explicit.  Used for: design region, true interface, wrong interface, protected set, contaminated positions, mutated positions.

7. **`ContigSpec`** — the RFDiffusion contig string parsed into structured segments (fixed / denovo, start, end, length).  Knows: total receptor length, denovo segments, fixed segments, which native positions are anchored.  Used by RFDiffusion, MPNN sequence correction, HADDOCK, indices derivation.

8. **`StageResult`** — the result of one negsteer stage for one configuration: a list of `ProteinStructurePrediction`s (one per seed, count = `num_seeds`, NOT fixed at 3 — `num_seeds` is a user parameter, validated odd-only) plus the collective verdict for that stage (e.g. clean / contaminated / pose_collapses).  The stage's verdict drives whether the next stage runs.

9. **`NegativeSteeringRun`** — one full negsteer experiment for one MPNN sequence.  Owns its sequence of `StageResult`s (cold-start → steered-per-design → reversion-per-contaminated-design), plus the per-sequence outcome label and `n_pass` count.  Externally exposes: outcome, n_pass, tier, canonical PDB, passing-summary rows.

10. **`DesignCohort`** — the set of all MPNN sequences in one pipeline run, each with its `NegativeSteeringRun`, aggregated and ranked.  Knows tier A/B/C/none breakdowns, the composite-score ranking, who the survivors are, how to emit cross_sequence_summary.csv.  (Renamed from `Cohort` for clarity.)

11. **`OrthogonalMetrics`** — the AF3 + biophysical (BSA, hbonds) + Rosetta (Sc, ΔΔG) measurements for one `DesignCohort` survivor.  Knows the `passes_orthogonal_filters` gate (Sc + BSA + ΔΔG only — interface_plddt is NOT orthogonal, removed in commit `47cb9f2`).

12. **`PipelineParams`** — the user-facing parameter set.  Everything `bin/validate_params.py` checks.  Sourced from `params.X` in Nextflow / `params_example.yml`.  Knows how to validate itself.

13. **`PipelineInternalThresholds`** — the internal-threshold catalogue.  Every hard-coded threshold currently scattered through `bin/` (ra_eff < 5 Å, intact threshold, ipae cutoffs, contact cutoff, composite weights, etc.).  Default values live here; the user CAN override but normally doesn't.  This type is the home of the future "Task 47 threshold audit."

### Confirmed scoping decisions

- **AF3 and Boltz outputs share the structure type, not the metric types.**  ProteinStructurePrediction is predictor-agnostic; `BoltzConfidenceMetrics` and `AF3ConfidenceAggregate` are predictor-specific.
- **PositionSet carries its frame.**  Single type, `frame` field, mixing-frames operations raise.  Conversion is explicit (`set.in_frame("native")`).
- **`StageResult` is in the hierarchy.**  The three-tier-deep hierarchy is: `ProteinStructurePrediction` (one prediction) → `StageResult` (one stage's predictions for one configuration) → `NegativeSteeringRun` (one MPNN sequence's full chain of stages) → `DesignCohort` (all MPNN sequences in one run).
- **Pipeline parameters split.**  Two types: `PipelineParams` (user contract) and `PipelineInternalThresholds` (research methodology).

### Things considered and dropped from the type list

- `NativeComplex` (the input PDB) — used to derive things, but doesn't carry active behaviour.  Stays as a `Path` + a `ProteinStructurePrediction` view.
- `ReversionPlan` / `ReversionAttempt` — internal detail of `NegativeSteeringRun`.
- `Contamination` — just a `PositionSet` with a particular meaning.
- Plot constants (`COLOUR_TIER`, `_TIER_SORT_ORDER`, etc.) — constants, not types.  Go in a `plot_lib` module.
- `PerSeedResult` — discussed and rejected.  The right intermediate concept is `StageResult` (a collection of predictions for one configuration) rather than "one seed's mini-workflow."  See Session 6 below.

## User clarifications carried forward from Session 6

These are domain-truths the user has stated during the grill-me that the architecture spec MUST honour.  Going against any of these is a bug, not a stylistic choice.

### CL-1.  `no_reversion` is NOT cold-start

The aggregate outcome `no_reversion` fires whenever no seed in a sequence-group needed reversion.  This can happen in two distinct paths: (a) cold-start `skip_steering` (all cold-start seeds were structurally fine, steering never ran), OR (b) steering ran but every seed's steered prediction had zero contamination on mutated positions.  `no_reversion` is NOT a synonym for cold-start.  See `~/.claude/.../memory/project_no_reversion_semantics.md`.

### CL-2.  Stages run all `num_seeds` predictions; they are not "per-seed mini-workflows"

Each pipeline stage (cold-start, steered-per-design, reversion-per-unique-reverted-sequence) runs ALL `num_seeds` predictions.  Seed indices are RNG controls for individual predictions, NOT threads that flow through stages.  The decision about whether the next stage runs is made on the *collective* result of the current stage's predictions, not seed-by-seed.

### CL-3.  Reversion gating: majority-of-correctly-placed rule (**TO BE IMPLEMENTED — current code is wrong**)

Reversion should run for a design's reverted sequence if and only if the contaminated count is at least a majority of the correctly-placed seeds for that design.

For each steered design:
1. `n_correctly_placed` = number of seeds whose steered prediction is intact AND has `ra_eff < threshold`.
2. `n_contaminated` = number of those correctly-placed seeds that ALSO have mutated positions contacting the effector.
3. Reversion runs IFF `n_correctly_placed > 0` AND `n_contaminated >= ceil(n_correctly_placed / 2)`.

When reversion runs, it always runs all `num_seeds` predictions for the reverted sequence (no per-seed gating at the reversion stage).

Examples (with `num_seeds=3`):
- 3 correctly placed, 1 contaminated → 1 < ceil(3/2)=2 → **no reversion** (Tier-B-shape result stands; the 1 contaminated seed counts as a failure, NOT as part of `n_pass`).
- 3 correctly placed, 2 contaminated → 2 ≥ 2 → **reversion runs** (3 seeds of reversion).
- 2 correctly placed (1 wrong placement), 1 contaminated → 1 ≥ ceil(2/2)=1 → **reversion runs**.
- 2 correctly placed, 0 contaminated → **no reversion** (Tier-B holds).
- 0 correctly placed → **no reversion** (nothing to revert).

**Code today does it per-seed**: any single contaminated (design, seed) entry queues reversion for that design's reverted sequence.  This drives down Tier B / Tier A outcomes by re-running reversion when it shouldn't.  Phase 4 must change `cmd_build_contaminated` to apply the majority rule on grouping by design.

A new per-seed verdict bucket is needed for "correctly placed but contaminated and didn't trigger reversion under the majority rule" — call it `contaminated_skipped` for now.  It contributes neither to `n_pass` (so the design doesn't claim the contaminated seed as a pass) nor to the reversion-derived verdict set.  Naming to be confirmed in Step 2.

### CL-4.  Reversion ALWAYS runs `num_seeds` predictions (when it runs)

When reversion runs for a unique reverted sequence, it always runs all `num_seeds` predictions.  Never 1, never 2.  Verified in `reversion.py:449` (`for seed_idx in range(num_seeds)`).  This is the one piece of current code that already matches the intended behaviour.

### CL-5.  `num_seeds` is a parameter, not the constant 3

`num_seeds` is validated odd-only (current odd-validation in `bin/validate_params.py`, custom validator `_check_num_seeds_odd`).  Default is 3 for ease of explanation, but the spec uses `num_seeds` throughout, not the literal `3`.  No type or method should hard-code 3.

---

## Step 2 — Per-type interface specs

For each of the 13 types we work through:
- **Domain meaning** — one paragraph in non-CS language.  User confirms.
- **Fields** — the data the type carries.
- **Methods** — the behaviour the type exposes.
- **Boundary** — what the type explicitly does NOT do.
- **Code today** — where this concept lives in the current codebase (so the migration knows what to consolidate).

Order: ProteinStructurePrediction first (it's the spine), then build outward.

### 2.1 `ProteinStructurePrediction`

**Domain meaning.**  One prediction's structural payload — a Boltz PDB or an AF3 mmCIF — together with its chain composition.  Holds the sequence and atom coordinates and exposes everything derivable from them: per-chain sequences, contact residues between chains, RMSDs against another prediction, jaccard overlap of interfaces.  Does NOT carry predictor-specific confidence metrics (those are `BoltzConfidenceMetrics` / `AF3ConfidenceAggregate`) and does NOT know about its position in any larger experiment (no awareness of "I'm cold-start seed 0").

**Construction.**

```
ProteinStructurePrediction(
    path,                      # Path to .pdb (Boltz) or .cif (AF3).
    receptor_chain,            # str, e.g. "A".
    effector_chain,            # str, e.g. "B".
    predictor = "boltz" | "af3" | None,  # diagnostic only; no behaviour
)
```

Validates `path` exists at construction (fails fast with a clear error).  Atom coordinates and chain sequences are loaded lazily on first access and cached for the instance's lifetime.

**Methods — sequence access.**
- `receptor_sequence() -> str` — receptor chain's amino-acid sequence (1-letter).
- `effector_sequence() -> str` — effector chain's amino-acid sequence.
- `chain_sequence(chain: str) -> str` — escape hatch for the rare case where a caller really does need a non-receptor/non-effector chain.

**Methods — composition.**
- `receptor_length() -> int`
- `effector_length() -> int`

**Methods — coordinate access.**
- `receptor_ca_coords() -> np.ndarray` (shape `(N, 3)`)
- `effector_ca_coords() -> np.ndarray`

**Methods — contact analysis.**
- `receptor_contact_residues(cutoff: float) -> PositionSet` — receptor residues whose heavy atoms come within `cutoff` Å of any effector heavy atom.  Returns a `PositionSet` in frame `prediction`.

**Methods — pair operations (asymmetric, `_vs_` reads as "self's relation to other").**

All take `other: ProteinStructurePrediction` as the argument.  Both predictions need to be the same type, so either can be called on — values would differ where the operation is genuinely asymmetric (e.g. RMSD with `self` as the mobile structure aligned onto `other`).

- `ra_eff_vs(other) -> float` — receptor-aligned effector RMSD.  Aligns receptor of `self` onto receptor of `other`, then measures effector atom displacement.  The headline structural pose metric.
- `independent_receptor_rmsd_vs(other) -> float` — receptor-only Cα RMSD (`self` aligned to `other` on receptor atoms).
- `independent_effector_rmsd_vs(other) -> float` — effector-only Cα RMSD.
- `jaccard_interface_vs(other, cutoff: float) -> float` — jaccard overlap between this prediction's receptor contact residues and `other`'s receptor contact residues, at the given cutoff.
- `weighted_jaccard_interface_vs(other, sigma: float, mu: float) -> float` — Gaussian-weighted variant.

**Boundary — what this type explicitly does NOT do.**
- Confidence metrics (pLDDT, PAE, ipTM, etc.) — those live on the predictor-specific confidence types.
- Threshold checks (`is_intact()`, `passes_ra_eff_threshold()`, `is_contaminated()`) — those are policy decisions, not properties of the prediction.  Callers do `pred.ra_eff_vs(truth) < thresholds.steered_ra_eff` themselves.  See open question below.
- Knowledge of the prediction's role (cold-start, steered, reverted).  That's `StageResult`'s job.
- Knowledge of mutations applied — that's `DesignedSequence`.
- Loading or interpreting sidecar JSON files (confidence/PAE/plddt) — those are the confidence types' responsibility.
- Writing files — read-only.

**Code today (consolidation targets — every one of these collapses into this type).**

| Function | File | Notes |
|---|---|---|
| `get_chain_sequence` | `bin/boltz2_negative_steering.py` | The canonical version; used by many. |
| `read_ca_atoms` | `bin/boltz2_negative_steering.py` | Returns CAEntry dataclass list. |
| `find_contact_residues_heavy` | `bin/boltz2_negative_steering.py` | The negsteer version is canonical (per Session 5 Q113). |
| `jaccard`, `weighted_jaccard` primitives | `bin/boltz2_negative_steering.py` | |
| `binding_rmsds` | `bin/boltz2_negative_steering.py` | |
| `read_ca_atoms` | `bin/rfdiffusion_filter.py` | Independent reimplementation. |
| `_chain_ca_coords` | `bin/parse_af3_output.py` | AF3 CIF reader. |
| `_read_ca_chain` | `bin/reversion.py` | Length-only guard helper. |
| `find_contact_residues_heavy` | `bin/compute_metrics.py` | Second copy. |
| `_extract_chain_seq` | `bin/extract_survivor_manifest.py` | Third sequence extractor. |
| `_extract_chain_sequence` | `bin/build_control_sequences.py` | Fourth. |
| `get_pdb_sequence` | `bin/pipeline_correct_sequences.py` | Fifth. |
| `_WJ_MU=4.0, _WJ_TWO_SIGMA_SQ=2.25, _gaussian_contact_weight, _collect_chain_cb_positions, _build_weighted_pair_map` | `bin/compute_metrics.py` AND `bin/compute_interface_metrics.py` | Duplicated weighted-jaccard implementation. |

### 2.2 `ContigSpec`

**Domain meaning.**  An RFDiffusion contig string parsed into structured per-chain segments.  Canonical input form: `A1-10/5-7/A15-20 B` (see memory `project_contig_string_format.md`).

- Within a chain, segments are slash-separated.
- Fixed segments are chain-prefixed native residue ranges: `A1-10` = chain A, native residues 1 through 10.
- De novo segments have no chain prefix and are written as length ranges: `5-7` = "fill in a de novo region of length 5 to 7 residues."
- Bare chain identifier (e.g. ` B` after a space) = entire chain included, no design regions.
- Different chains are space-separated.  Never commas.

ContigSpec represents the *resolved* form: after RFDiffusion has sampled a specific length for each de novo segment, the contig becomes a concrete layout.  The constraint form (with length ranges) is an input string to RFDiffusion and is NOT a runtime type — it's just plumbed through as text.  Each DesignedBackbone has exactly one ContigSpec describing its resolved layout.

**Construction.**

```
ContigSpec(
    chains,          # ordered list of ContigChain (one per chain in the contig)
    contig_string,   # original resolved string, kept for round-trip
)

ContigChain(
    chain_id,        # str, e.g. "A"
    segments,        # ordered list of FixedSegment | DeNovoSegment
)

FixedSegment(start: int, end: int)            # native residue numbers; inclusive
DeNovoSegment(length: int)                    # specific resolved length (not a range)
```

Static constructor: `ContigSpec.from_resolved_string(s: str) -> ContigSpec` — the canonical parser.

**Methods — frame conversion (per Q7 decision).**

Per-position conversion lives here, not on PositionSet.  PositionSet's `in_frame()` calls these.

- `designed_to_native(pos: int, chain: str) -> int | None` — returns the native residue number, or `None` if `pos` falls in a de novo segment (no native counterpart).
- `native_to_designed(pos: int, chain: str) -> int | None` — returns the designed residue number, or `None` if `pos` is not anchored by any fixed segment.
- `prediction_to_designed(pos: int, chain: str) -> int` — Boltz/AF3 always emit per-chain 1-based; for the receptor chain this is identity (Boltz numbering matches designed numbering when the receptor PDB went straight into Boltz).  Documented for clarity.

**Methods — structural queries.**

- `chain_length(chain: str) -> int` — total residues in the resolved chain (sum of fixed-segment lengths + de novo lengths).
- `fixed_segments(chain: str) -> list[FixedSegment]`
- `denovo_segments(chain: str) -> list[DeNovoSegment]`
- `design_region(chain: str) -> PositionSet` — the de novo positions in the *designed* frame.  Convenience for the common ask.
- `fixed_anchor_positions(chain: str) -> PositionSet` — the fixed positions in the *designed* frame (i.e. positions that DO have native counterparts).
- `is_design_region_position(pos: int, chain: str) -> bool` — predicate.

**Boundary — what this type explicitly does NOT do.**

- Carry the constraint form (length ranges).  ContigSpec is always resolved.  The constraint is an input string used by RFDiffusion's launcher only.
- Know about residue identities or coordinates — it's pure positional metadata.
- Validate against any specific PDB (no "does the contig match what's in this file" check).  Validation is the caller's job at the boundary between contig and PDB.

**Code today (consolidation targets).**

Four contig parsers today.  ContigSpec absorbs them:

| Parser | File | Notes |
|---|---|---|
| `parse_block_segments` | `bin/contig_utils.py` | Closest to the canonical parser; used by `rfdiffusion_contigs.py` and `rfdiffusion_filter.py`. |
| `_parse_contigs` | `bin/derive_input_design_region.py` | Independent reimplementation. |
| `parse_contig_segments` | `bin/haddock3_prepare.py` | HADDOCK-specific subset. |
| `parse_contig_segments` | `bin/pipeline_correct_sequences.py` | MPNN-side reimplementation; subtly different. |

---

### 2.3 `PositionSet`

**Domain meaning.**  A set of receptor residue positions, tagged with which numbering frame the integers refer to.  Frames in use today: `native` (input PDB's residue numbers; may have gaps), `designed` (RFDiffusion's 1-based contiguous renumbering), `prediction` (Boltz / AF3's 1-based per-chain).  Used everywhere the codebase talks about "which residues" — design region, true interface, wrong interface, protected set, mutated positions, contaminated positions, contact residues.

**Construction.**

```
PositionSet(
    positions,      # iterable of ints; deduplicated and sorted on construction
    chain,          # str, e.g. "A"
    frame,          # one of "native" | "designed" | "prediction"
    contig = None,  # optional ContigSpec; required for in_frame() conversion
)
```

Immutable.  Every operation returns a new `PositionSet`.

**Methods — set operations.**

Both operands must be in the same frame and on the same chain.  Mismatched frame or chain raises a clear error.

- `union(other) -> PositionSet`
- `intersection(other) -> PositionSet`
- `difference(other) -> PositionSet`
- `__contains__(pos: int) -> bool` (so `pos in myset` works)
- `__or__`, `__and__`, `__sub__` operator aliases.

**Methods — frame conversion.**

- `in_frame(target_frame: str) -> PositionSet` — returns a new PositionSet in `target_frame`.  Requires `self.contig is not None`; raises a clear error otherwise.  Default behaviour for positions with no counterpart in the target frame: silently dropped (returned set may be smaller).  Pass `strict=True` to raise on any unmappable position instead.

**Methods — predicates and access.**

- `is_empty() -> bool`
- `__len__() -> int`
- `__iter__()` yields positions in sorted order.
- `as_sorted_list() -> list[int]`

**Methods — I/O.**

- `to_chimerax_string() -> str` — `"/A:5,7,22"` format used widely in the codebase.
- `PositionSet.from_chimerax_string(s, frame, contig=None) -> PositionSet`
- `to_text_file(path)` — writes integers one per line, sorted.  Matches current `*_design_region.txt` / `*_true_interface.txt` files.  Chain and frame are NOT serialised (callers know them); for round-trip discipline, the reader must specify both.
- `PositionSet.from_text_file(path, chain, frame, contig=None) -> PositionSet`

**Boundary — what this type explicitly does NOT do.**

- Carry per-residue properties (amino acid identity, coordinates, hotspot flags).  Positions are just integers; semantic meaning is conveyed by the variable name (`design_region`, `true_interface`, `mutated_positions`).
- Bounds-check positions against any specific receptor length.  Callers validate against their specific receptor when needed.
- Silently convert between frames.  Cross-frame operations raise.
- Mutate after construction.

**Code today (consolidation targets).**

Position-set logic is scattered across these files today; PositionSet absorbs all of it:

| Where | What |
|---|---|
| `bin/derive_design_region.py` | Builds the design-region position set; writes `_design_region.txt`. |
| `bin/derive_true_interface.py` | Builds the true-interface position set; writes `_true_interface.txt`. |
| `bin/derive_input_design_region.py` | Same for the input PDB (control workflow). |
| `bin/boltz2_negative_steering.py` | `find_contact_residues_heavy`, mutated-position arithmetic, protected-set union. |
| `bin/boltz2_iterate_steering.py` | `mutated_positions_set & contact_set` for contamination detection. |
| `bin/extract_passing.py` | `_summarize_mutations_across_seeds` parses and aggregates `/A:5,7,22` strings. |
| `bin/compute_metrics.py` | Second copy of `find_contact_residues_heavy`. |

---

### 2.4 `BoltzConfidenceMetrics`

**Domain meaning.**  The confidence-flavoured measurements Boltz emits alongside each prediction — pLDDT variants, PAE summaries, ipTM, ipSAE.  Pairs one-to-one with a `ProteinStructurePrediction` from Boltz.  Does NOT carry the structural payload (that's `ProteinStructurePrediction`).

**Construction.**

```
BoltzConfidenceMetrics.from_prediction_dir(path: Path) -> BoltzConfidenceMetrics
```

Reads `confidence_<name>_model_<sample>.json`, `pae_<name>_model_<sample>.json`, `plddt_<name>_model_<sample>.json` from the Boltz output directory.  Lazy reads (matches ProteinStructurePrediction's behaviour); validates the directory exists at construction.

**Methods — exposed metrics (read-only properties).**

Selected-sample metrics — the ones currently surfaced through the pipeline:
- `avg_plddt`, `complex_plddt`, `interface_plddt`
- `pae_mean`, `ipae`, `pae_pass_frac`
- `ptm`, `iptm`, `actifptm`
- `ipsae_ab`, `ipsae_ba`, `ipsae_min`, `ipsae_ab_15`, `ipsae_ba_15`, `ipsae_min_15`

Each is a `float | None`; `None` if the sidecar JSON was unparseable or the field was missing.

**Methods — derived.**
- `confidence_flag(thresholds) -> str` — returns one of `ok | low_pass_frac | low_iptm | low_plddt | high_ipae | multiple`.  Threshold values come from `PipelineInternalThresholds`.

**Boundary.**
- Does not own the structure (no atom coords, no contacts).  Composes naturally with `ProteinStructurePrediction`.
- Does not aggregate across multiple predictions — that aggregation happens at the `StageResult` level.
- Does not interpret per-residue arrays (per-position pLDDT, full PAE matrices) — only the scalar summaries the pipeline currently uses.  If a future stage wants per-residue, extend then.

**Code today.**  `bin/compute_metrics.py` reads these JSONs and emits the scalars used downstream.  Several plot scripts and analysis steps re-read the same JSONs.  Consolidates here.

---

### 2.5 `AF3ConfidenceAggregate`

**Domain meaning.**  The aggregated confidence + structural quality from ONE AF3-no-MSA call.  Unlike Boltz (where one call → one kept prediction), AF3 emits multiple per seed × diffusion sample, and `parse_af3_output.py` aggregates across them into best / mean summaries.  This type wraps that aggregate.  Pairs with the *best-ra_eff* `ProteinStructurePrediction` from the AF3 set — for the cases (rare) where the single best prediction needs structural follow-up.

**Construction.**

```
AF3ConfidenceAggregate.from_output_dir(path: Path, ground_truth: Path) -> AF3ConfidenceAggregate
```

Walks the AF3 output tree, parses every per-sample mmCIF + its `summary_confidences.json`, computes ra_eff against `ground_truth` for each, then aggregates.

**Methods — exposed metrics (read-only properties).**
- `best_ra_eff`, `mean_ra_eff` — min and mean across all (seed × sample) predictions.
- `best_iptm`, `mean_iptm` — max and mean across all.
- `n_correct_interface(threshold) -> int` — count of predictions where `ra_eff < threshold`.  Threshold from `PipelineInternalThresholds`.
- `total_predictions` — total mmCIFs parsed.
- `failures` — comma-separated failure reasons.

**Methods — access to underlying predictions.**
- `best_prediction() -> ProteinStructurePrediction` — the single mmCIF with the lowest ra_eff against ground truth.
- `all_predictions() -> list[ProteinStructurePrediction]` — full set.

**Boundary.**
- Does not gate (the AF3 ra_eff is informational only in `passes_orthogonal_filters` per Session 5 Q118).
- Does not compare against Boltz — that's done outside, in the cohort summary.

**Code today.**  `bin/parse_af3_output.py` is the canonical existing implementation; `bin/orthogonal_metrics_plots.py` reads the aggregate fields.  Consolidates here.

---

### 2.6 `DesignedBackbone`

**Domain meaning.**  One RFDiffusion output — a designed receptor backbone that has STRUCTURE (atom coordinates) but no SEQUENCE yet (the residue identities for the de novo region come later from MPNN).  Carries the contig that produced it and the Rosetta-filter score that gated it into the next stage.

**Construction.**

```
DesignedBackbone(
    pdb_path,          # Path to the RFDiffusion-output PDB
    design_id,         # str, e.g. "design_03"
    contig,            # ContigSpec (resolved form, this design's actual layout)
    structure,         # ProteinStructurePrediction (composed; reads the PDB)
    sc_score,          # float | None — Rosetta InterfaceAnalyzer Sc; None if filter pass didn't run
    rfdiff_metrics,    # dict | None — raw rfdiffusion_metrics.json entry for this design
)
```

Composes `ProteinStructurePrediction` (per Q8 — composition, not inheritance).  All structural access goes through `backbone.structure.X`.

**Methods.**
- `design_region() -> PositionSet` — delegates to `contig.design_region(chain="A")`, returns the de novo positions on the receptor in `designed` frame.
- `fixed_anchor_positions() -> PositionSet` — delegates to `contig.fixed_anchor_positions(chain="A")`.
- `receptor_length() -> int` — delegates to `contig.chain_length("A")`.
- `receptor_chain` / `effector_chain` — properties on the composed structure.

**Boundary.**
- Carries no sequence.  Adding sequence produces a `DesignedSequence`.
- Carries no MPNN-side data (scores, QC).
- Owns the structural file, but does NOT own any Boltz prediction OF itself (those come later in the negsteer stage and live as `ProteinStructurePrediction`s inside `StageResult`s).
- Filter checks (`if sc_score >= sc_threshold`) live in caller code using `PipelineInternalThresholds`.

**Code today.**  `bin/rfdiffusion_filter.py` and `bin/rosetta_filter_collect.py` produce the data that populates this type; consumers in `bin/proteinmpnn/*` read it.  Replaces the implicit "design_id + path + sc" tuple passed through dicts.

---

### 2.7 `DesignedSequence`

**Domain meaning.**  One MPNN-designed sequence overlaid on one `DesignedBackbone`.  The receptor amino-acid string with native residues at fixed positions and MPNN-designed residues at de novo positions, plus MPNN's quality metadata for it.  Identified by `design_<N>_seq_<M>` — N is the backbone, M is the MPNN sequence index within that backbone.

**Construction.**

```
DesignedSequence(
    backbone,             # DesignedBackbone
    sequence_id,          # str, e.g. "design_03_seq_01"
    corrected_receptor,   # str — full receptor sequence; native at fixed, MPNN at design positions
    designed_residues,    # str — pipe-joined MPNN slice at design positions only (e.g. "RGHIK")
    native_residues,      # str — pipe-joined native slice at the same positions (for comparison)
    mpnn_score,           # float
    qc_metadata,          # dict — poly-X count, % identity vs native, etc.
)
```

**Methods.**
- `to_fasta() -> str` — receptor record in FASTA form.
- `design_region() -> PositionSet` — delegates to `backbone.design_region()`.
- `mutations_vs_native() -> list[tuple[int, str, str]]` — list of `(designed_position, native_aa, designed_aa)` differences at design positions.
- `n_changes() -> int` — count of positions where designed != native.

**Boundary.**
- Carries no Boltz prediction OF this sequence.  Those come from the negsteer stage and live in `StageResult.predictions`.
- Does not run QC itself — QC happens upstream; this type just carries the metadata.
- Filter checks (`passes_qc`) live in caller code using `PipelineInternalThresholds`.

**Code today.**  `bin/pipeline_correct_sequences.py` produces `mpnn_corrected.fasta` + the metadata rows that become `scored_metadata.csv`; `bin/mpnn_sequence_qc.py` adds the QC metadata.  Today this is plumbed through CSVs and dicts; DesignedSequence absorbs that.

---

### 2.8 `StageResult`

**Domain meaning.**  One stage's worth of predictions for one configuration.  A stage is `cold_start`, `steered`, or `reversion`.  A configuration is "the receptor sequence being predicted in this stage" — for cold-start it's the wild-type baseline for one MPNN sequence; for steered it's one steered-design slot; for reversion it's one unique reverted sequence.  Each stage runs `num_seeds` predictions for the configuration, with the user-configured `num_seeds` parameter (validated odd-only).

The collective verdict of the stage drives whether the *next* stage runs.  Each stage has its own gating rule (see `triggers_next_stage` below).

**Construction.**

```
StageResult(
    stage_type,         # "cold_start" | "steered" | "reversion"
    config_id,          # str — meaning is stage-type-specific:
                        #   cold_start: the MPNN sequence ID, e.g. "design_03_seq_01"
                        #   steered:    the steered-design slot, e.g. "design_03"
                        #   reversion:  the reverted-sequence label, e.g. "rev_design_03_<hash>"
    predictions,        # list[ProteinStructurePrediction] of length num_seeds
    confidences,        # list[BoltzConfidenceMetrics] parallel to predictions
    seed_indices,       # list[int] — the RNG seeds used (parallel to predictions)
    mutated_positions,  # PositionSet | None — only for steered/reversion;
                        #   the receptor positions that differ from cold-start
    applies_to_designs, # list[str] | None — REVERSION ONLY (per Step 4.3):
                        #   the steered design_ids whose contamination
                        #   produced this unique reverted sequence.  Used by
                        #   NegativeSteeringRun._reversion_for_design() to
                        #   map a steered design to its reversion attempt,
                        #   accounting for the dedup-by-reverted-sequence
                        #   that write_reversion_plan does today.
)
```

`num_seeds` is the parameter, not a constant.  `len(predictions) == len(confidences) == len(seed_indices) == num_seeds`.

**Methods — per-seed.**

- `per_seed_verdicts(thresholds, contamination_positions=None) -> list[str]` — applies the per-seed verdict rule for the stage type:
  - `cold_start`: each seed → `clean` (intact + ra_eff < cutoff) or `fail_structural`.
  - `steered`: each seed → `clean_steered` (intact + ra_eff < cutoff + no mutated contacts) or `contaminated` (correct placement BUT mutated positions contact effector) or `pose_collapses` (intact failed or ra_eff ≥ cutoff) or `no_data`.
  - `reversion`: each seed → `pose_holds` (reverted prediction intact + ra_eff < cutoff + no new contamination) or `new_contamination` (reverted prediction has mutated positions contacting effector) or `pose_collapses` (intact failed or ra_eff ≥ cutoff) or `no_data`.

  `contamination_positions` is the gated position set (design_region ∪ true_interface) — required for steered and reversion stages.

**Methods — stage-level.**

- `triggers_next_stage(thresholds, contamination_positions=None) -> bool` — applies the stage-specific gating rule:
  - `cold_start`: returns `True` (steering should run) iff fewer than majority of cold-start seeds are `clean`.  (Existing rule.)
  - `steered`: returns `True` (reversion should run for this design) iff `n_correctly_placed > 0` AND `n_contaminated >= ceil(n_correctly_placed / 2)`.  (NEW rule from CL-3 — replaces the per-seed gate currently in `cmd_build_contaminated`.)
  - `reversion`: returns `False` always (no next stage).

- `n_correctly_placed(thresholds) -> int` — count of seeds with intact + ra_eff < cutoff.  Used by the steered-stage gating rule.
- `n_contaminated(thresholds, contamination_positions) -> int` — count of correctly-placed seeds with mutated positions in `contamination_positions`.  Used by the steered-stage gating rule.
- `n_pass(thresholds, contamination_positions=None) -> int` — count of seeds whose final verdict is pass-equivalent.  Pass-equivalent depends on stage and on whether reversion ran (see `NegativeSteeringRun.aggregate_n_pass` for cross-stage logic).
- `canonical_prediction() -> ProteinStructurePrediction` — the representative one for this stage (best by composite, or by `ra_eff`, depending on stage policy).

**Boundary.**

- Does not know about other stages.  Cross-stage logic (mapping reversion verdicts back to the steered seeds they replaced; deciding final per-seed verdicts after reversion runs) lives in `NegativeSteeringRun`.
- Does not orchestrate Boltz runs.  StageResult is constructed from already-completed predictions; it is read-only.
- Does not own thresholds.  All threshold checks take `PipelineInternalThresholds` as an argument.

**Code today.**

The per-seed and per-stage logic currently lives scattered across `bin/boltz2_negative_steering.py` (`cmd_plan`, `cmd_finalize`), `bin/boltz2_iterate_steering.py` (`cmd_build_contaminated`, `_per_seed_verdict_breakdown`, `_classify_outcome`), and the orchestration shell `bin/negative_steering_run_one.sh`.  StageResult absorbs the per-stage classification and gating logic; the orchestration stays in shell + Nextflow (see `NegativeSteeringRun` §2.9 below).

**Implementation note — the CL-3 fix lands here.**  The Phase 4 migration of `cmd_build_contaminated` becomes: instead of iterating per (design, seed) and emitting a contaminated entry per match, group by design, build a `StageResult` for each design's steered seeds, call `stage_result.triggers_next_stage(thresholds, contamination_positions)` to decide whether to queue reversion for that design.  Per-seed verdicts that don't trigger reversion (contaminated seeds in a design that doesn't pass the majority rule) carry the `contaminated` label into the aggregate, where they count as failures, not as `n_pass`.

---

### 2.9 `NegativeSteeringRun`

**Domain meaning.**  One MPNN sequence's complete negsteer experiment: the cold-start `StageResult`, every steered-design `StageResult`, every reversion `StageResult`, plus the aggregation logic that produces the single per-sequence outcome.  Read-only — constructed from a completed workdir, does NOT execute anything (orchestration stays in `negative_steering_run_one.sh` + the Nextflow process body per Q10 decision).

One per MPNN sequence in a pipeline run, plus one per negative control (`control_scrambled`, `control_polyA`).

**Construction.**

```
NegativeSteeringRun.from_workdir(workdir: Path) -> NegativeSteeringRun
```

Parses `cycle_0/plan.json`, `aggregated_results.csv`, `passing_summary.csv`, the per-stage Boltz outputs under `cycle_0/initial*`, `cycle_0/steered/*`, `cycle_0/reversions/*`, and the sidecars (`row_type.txt`, `run_one_runtime_sec.txt`).

```
NegativeSteeringRun(
    mpnn_sequence_id,    # str — e.g. "design_03_seq_01" or "control_scrambled"
    workdir,             # Path — the per-sequence workdir under runs/
    row_type,            # "steered" | "control_scrambled" | "control_polyA"
    cold_start,          # StageResult (stage_type="cold_start", num_seeds predictions)
    steered,             # dict[design_id -> StageResult] — one per steered-design slot
    reversion,           # dict[reverted_seq_id -> StageResult] — one per unique reverted sequence
    run_one_runtime_sec, # float — wall-clock of the whole negsteer chain
)
```

**Methods — per-seed cross-stage aggregation.**

- `per_seed_final_verdicts(thresholds) -> list[dict]` — produces the final verdict for each (design × seed).  Rule: if reversion ran for that design, use the reversion verdict; otherwise use the steered verdict.  Each entry: `{design_id, seed_index, stage_run: "steered" | "reverted", verdict}`.

- `n_pass(thresholds) -> int` — count of seeds whose final verdict is pass-equivalent (`clean_steered` or `pose_holds`).
- `n_seeds_total() -> int` — same as `cold_start.num_seeds` (= `num_seeds` parameter).

**Methods — sequence-level outcome.**

- `outcome(thresholds) -> str` — aggregate label.  Returns one of `no_reversion` | `pose_holds` | `pose_collapses` | `new_contamination` | `singleton`.  Computed from per-seed final verdicts (see `_classify_outcome` in current code, which the migration extracts here).
- `outcome_reason(thresholds) -> str` — diagnostic explanation of the outcome decision.
- `tier(thresholds) -> str` — `A` | `B` | `C` | `none`.  Derived from `n_pass / n_seeds_total` per the post-rename rule (no `no_reversion → A` shortcut; see CL-3).

**Methods — representative selection (for cohort plots).**

- `canonical_prediction(thresholds) -> ProteinStructurePrediction | None` — the "best" single prediction across all stages, used as the representative in `DesignCohort` and orthogonal-stage processing.  Selected per current tier-then-composite policy.

**Methods — workdir-cached state (added per Step 4.1).**

- `contamination_positions() -> PositionSet` — cached load of `design_region ∪ true_interface` from `cycle_0/plan.json`.  Required as input to per-seed verdict and CL-3 majority-rule checks.
- `_reversion_for_design(design_id: str) -> StageResult | None` — private helper (added per Step 4.3); walks `self.reversion.values()` looking for the stage whose `applies_to_designs` includes `design_id`.  Returns the matching reversion StageResult, or `None` if this design didn't trigger reversion.

**Boundary.**

- Read-only.  Does NOT run Boltz, does NOT submit SLURM jobs, does NOT manage subprocesses.
- Does not know about other NegativeSteeringRuns (no cohort-level methods here — those live in `DesignCohort`).
- Threshold checks (`tier`, `outcome`) take `PipelineInternalThresholds` as an argument; no hidden state.
- Does not own orthogonal metrics for this sequence; those are a separate type computed from `canonical_prediction()` downstream.

**Code today.**  `bin/cross_sequence_summary.py`'s representative-picking + per-sequence aggregation + tier logic; `bin/boltz2_iterate_steering.py`'s `_classify_outcome` + `_per_seed_verdict_breakdown` + aggregator chain.  All of that read-side logic absorbs into this type.

---

### 2.10 `DesignCohort`

**Domain meaning.**  Every MPNN sequence in one pipeline run, plus the negative controls — each carrying its `NegativeSteeringRun`.  The thing `cross_sequence_summary.csv` represents.  Knows tier breakdowns, the composite-score ranking across sequences, and which sequences are "survivors" eligible for the orthogonal validation stage.

**Construction.**

```
DesignCohort.from_runs_directory(runs_dir: Path) -> DesignCohort
DesignCohort(runs: list[NegativeSteeringRun])
```

`from_runs_directory` walks `<runs_dir>/<seq_name>/` and constructs one `NegativeSteeringRun` per subdirectory.

**Methods — sequence selection.**

- `steered_runs() -> list[NegativeSteeringRun]` — only the steered rows (excludes controls).
- `control_runs() -> list[NegativeSteeringRun]` — only the controls.
- `by_tier(tier, thresholds) -> list[NegativeSteeringRun]` — sequences at the given tier (A/B/C/none).
- `survivors(thresholds) -> list[NegativeSteeringRun]` — tier A/B/C steered rows.  The set fed to the orthogonal stage.

**Methods — ranking.**

- `ranked_by_composite(thresholds) -> list[NegativeSteeringRun]` — sorted by composite score descending.  Used by `cross_rank_by_composite`.
- `ranked_by_ra_eff(thresholds) -> list[NegativeSteeringRun]` — parallel ranking by raw `ra_eff`.
- `tier_breakdown(thresholds) -> dict[str, int]` — `{"A": 3, "B": 5, "C": 1, "none": 63}` for one pipeline run.

**Methods — output.**

- `to_cross_summary_csv(path, thresholds, scored_metadata=None)` — writes `cross_sequence_summary.csv`.  Optional `scored_metadata` join (the MPNN-sequence join from the rename commit `05f0c07`).

**Methods — MPNN-sequence join.**

- `attach_designed_sequences(scored_metadata: dict)` — populates each NegativeSteeringRun's `DesignedSequence` reference from `scored_metadata.csv`.  Optional; without it, `corrected_receptor`/`designed_residues`/`native_residues` columns come out blank in the CSV.

**Boundary.**

- Does not know about orthogonal metrics.  Survivor selection produces a list of `NegativeSteeringRun`s that the orthogonal stage processes externally.
- Does not generate plots; plot scripts call methods on this type to get the data they need.
- Does not own thresholds; takes `PipelineInternalThresholds` as a method argument.

**Code today.**  `bin/cross_sequence_summary.py`'s `aggregate` function + CLI; the cohort-summary plot scripts' iteration over CSV rows.  Replaces dict-of-fieldnames flowing through CSVs with method calls on this type.

---

### 2.11 `OrthogonalMetrics`

**Domain meaning.**  The AF3 + biophysical + Rosetta measurements for ONE cohort survivor.  Independent of Boltz (the original prediction tool) so disagreement is a meaningful flag.  Knows the `passes_orthogonal_filters` gate — strictly Sc + BSA + ΔΔG (interface_plddt is Boltz-derived and NOT in the orthogonal gate per Session 5 Q118 and commit `47cb9f2`).

One per cohort survivor (tier A/B/C steered run).

**Construction.**

```
OrthogonalMetrics.from_orthogonal_outputs(workdir, mpnn_sequence_id) -> OrthogonalMetrics
```

Reads the per-survivor outputs of `NEGSTEER_BIOPHYSICAL_METRICS`, `NEGSTEER_ROSETTA_METRICS`, `AF3_PARSE_OUTPUT` from the standard locations.

```
OrthogonalMetrics(
    mpnn_sequence_id,        # str
    canonical_prediction,    # ProteinStructurePrediction — the one being validated
    af3,                     # AF3ConfidenceAggregate
    bsa,                     # float | None — Rosetta dSASA_int (Å²)
    hbonds,                  # int | None — interface H-bond count
    sc,                      # float | None — Rosetta InterfaceAnalyzer shape complementarity
    ddg,                     # float | None — Rosetta dG_separated (ΔΔG of binding)
)
```

**Methods — the gate.**

- `passes_orthogonal_filters(thresholds) -> bool` — True iff `sc >= thresholds.orthogonal_sc_min` AND `bsa >= thresholds.orthogonal_bsa_min` AND `ddg <= thresholds.orthogonal_ddg_max`.  AF3 is NOT in this gate; it's surfaced separately as informational disagreement.
- `failed_filter_names(thresholds) -> list[str]` — returns the names of any orthogonal metrics that failed the gate (e.g. `["sc", "ddg"]`).  Empty if `passes_orthogonal_filters` is True.

**Methods — informational.**

- `af3_disagrees(thresholds) -> bool` — True iff `af3.best_ra_eff >= thresholds.orthogonal_af3_ra_max`.  Diagnostic only; never used in the gate.
- `to_summary_row() -> dict` — flat dict for CSV output, matches existing `survivors_with_orthogonal_metrics.csv` schema.

**Boundary.**

- Does not own the prediction it validates (just a reference to a `ProteinStructurePrediction`).
- Does not own thresholds; takes `PipelineInternalThresholds` as a method argument.
- Does not depend on `DesignCohort` — produced one-per-survivor, consumed by cohort-summary plot or further analysis.

**Code today.**  `bin/run_biophysical_metrics.py`, `bin/run_rosetta_metrics.py`, `bin/merge_orthogonal_metrics.py`, `bin/parse_af3_output.py`.  Per Session 5 the test-copy divergence (`tests/orthogonal_metrics/merge_orthogonal_metrics.py`) is eliminated as part of this migration.

---

### 2.12 `PipelineParams`

**Domain meaning.**  The validated, user-facing parameter set for one pipeline run.  Sourced from Nextflow `params.X` (resolved from `params_example.yml` + user overrides).  Every parameter has a type, a valid range, and a default; the type validates itself at construction and refuses to operate with bad inputs.  This is the user contract — changing a field touches user-visible behaviour.

**Construction.**

```
PipelineParams.from_nextflow_json(path: Path) -> PipelineParams
PipelineParams.from_dict(d: dict) -> PipelineParams
```

Both routes validate via the existing `bin/validate_params.py` `PARAM_SPECS` catalogue (66 specs, 59 with concrete constraints).  Construction raises `ParamValidationError` (carrying all error messages, not just the first) when any spec fails.

**Fields (sample — the full set mirrors `PARAM_SPECS`).**

- Branch/input paths: `pdb_file`, `receptor_input`, `effector_input`, `receptor_chain`, `effector_chain`, …
- RFDiffusion: `num_designs`, `rfdiff_iterations`, `rfdiff_checkpoint`, `symmetry`, guiding-potential weights, …
- ProteinMPNN: `num_seqs`, `mpnn_sampling_temp`, `rm_aa`, `mpnn_top_n`, …
- Negsteer: `negsteer_mode`, `negsteer_n_designs`, `negsteer_num_seeds` (odd-validated), `negsteer_n_cycles`, …
- Controls + filters: `run_negative_controls`, `controls_warning_*`, `orthogonal_filter_*`, …
- Infrastructure: `max_boltz2_parallel`, `outdir`, `project_name`.

**Methods.**

- `as_dict() -> dict` — flat snapshot of all values.
- `validate() -> list[str]` — re-runs every spec; returns the empty list on success.

**Boundary.**

- Carries USER-facing params only.  Internal research thresholds (ra_eff cutoffs, intact thresholds, composite weights) live in `PipelineInternalThresholds`.  The split is by audit trail: changing a `PipelineParams` field affects the user contract; changing a `PipelineInternalThresholds` field is a research-methodology change.
- Read-only after construction.

**Code today.**  `bin/validate_params.py` (already exists from commit `05f0c07`).  This type is a thin runtime wrapper around the existing validator's `PARAM_SPECS`.

---

### 2.13 `PipelineInternalThresholds`

**Domain meaning.**  The catalogue of every hard-coded threshold currently scattered through `bin/`.  Defaults live here, in one Python module.  The user CAN override via parameter file or environment, but normally doesn't — these are research-methodology defaults the user has chosen as good values for their pipeline.

This is the home of the upcoming Task 47 threshold audit.  Changing a threshold here is a research-methodology change with an audit trail.

**Construction.**

```
PipelineInternalThresholds.default() -> PipelineInternalThresholds
PipelineInternalThresholds.from_overrides(d: dict) -> PipelineInternalThresholds  # default ∪ overrides
```

Defaults are Python literals defined in the type's module (chosen because they're part of the methodology and live in version control with the code, not in a YAML the user might edit per-environment).

**Fields (initial enumeration; populated/extended during the Task 47 audit).**

Structural thresholds:
- `steered_ra_eff = 5.0` — ra_eff threshold for "correctly placed" steered prediction.
- `reverted_ra_eff = 5.0` — same for reverted prediction.
- `cold_start_rmsd_threshold = 6.0` — currently `negsteer_rmsd_threshold` param default.
- `intact_threshold = 5.0` — receptor-fold-quality threshold (independent_receptor_rmsd).

Confidence thresholds (used by `confidence_flag` and friends):
- `iptm_min = 0.30`
- `complex_plddt_min = 0.70`
- `ipae_max = 15.0`
- `pae_pass_frac_min = 0.10`

Orthogonal thresholds (the gate):
- `orthogonal_sc_min = 0.55`
- `orthogonal_bsa_min = 600.0`
- `orthogonal_ddg_max = ?` (Bennett reference; populate during Task 47)
- `orthogonal_af3_ra_max = 5.0` (informational)

Contact / interface:
- `contact_cutoff = 4.5` — heavy-atom contact cutoff (Å).
- `weighted_jaccard_pair_cutoff = 8.0`
- `weighted_jaccard_sigma = 1.5` / `weighted_jaccard_mu = 4.0` (existing constants).

Composite score:
- `composite_ra_eff_weight = 0.05` — currently the hard-coded `COMPOSITE_RA_EFF_WEIGHT` constant.

Filter values:
- `interface_plddt_trim_threshold = 50.0`
- `controls_warning_ipsae_max = 0.5`
- `controls_warning_ra_eff_min = 5.0`

(Complete enumeration done during the Task 47 audit; this list is the starting set.)

**Methods.**

- `as_dict() -> dict`
- `with_overrides(d: dict) -> PipelineInternalThresholds` — returns a new instance with overrides applied; original is unchanged (immutable).

**Boundary.**

- Carries internal research thresholds only.  User-facing params live in `PipelineParams`.
- Read-only after construction (override produces a new instance).
- Does not enforce thresholds itself; consumed by methods on other types (`StageResult.triggers_next_stage(thresholds)`, `OrthogonalMetrics.passes_orthogonal_filters(thresholds)`, etc.).

**Code today.**  Currently scattered: hard-coded literal `5.0` in `cmd_build_contaminated`'s `rmsd_threshold`, `0.05` in `_ranking_composite`, etc.  Phase 4 collects every such literal into this type as defaults, then refactors each usage site to read from a passed `PipelineInternalThresholds` instance.


---

## Step 3 — Type relationships

Three views: composition (who contains what), dependency (who needs to be implemented before what), and data flow (how a pipeline run threads through the types).

### Composition tree

```
DesignCohort
└── NegativeSteeringRun  [one per MPNN sequence + per control]
    ├── DesignedSequence  [the MPNN sequence being negsteered]
    │   └── DesignedBackbone
    │       ├── ProteinStructurePrediction  [composed — backbone-only PDB]
    │       └── ContigSpec
    ├── StageResult  [cold_start]
    │   ├── list[ProteinStructurePrediction]
    │   └── list[BoltzConfidenceMetrics]
    ├── dict[design_id → StageResult]  [steered]
    │   └── list[ProteinStructurePrediction] + list[BoltzConfidenceMetrics]
    └── dict[reverted_seq_id → StageResult]  [reversion — only when CL-3 majority rule fires]
        └── list[ProteinStructurePrediction] + list[BoltzConfidenceMetrics]

OrthogonalMetrics  [one per DesignCohort survivor; sibling of NegativeSteeringRun, not contained]
├── ProteinStructurePrediction  [the canonical prediction being validated]
└── AF3ConfidenceAggregate
    └── list[ProteinStructurePrediction]  [internal — the seed×sample mmCIFs]

PositionSet  [embedded wherever a position-set is returned/passed; carries
              optional ContigSpec reference for frame conversion]

PipelineParams              [singleton per pipeline run; passed top-down]
PipelineInternalThresholds  [singleton per pipeline run; passed to every method
                             that does a threshold check]
```

### Dependency-order table (drives migration order)

Each row lists what each type DEPENDS ON (must exist first).  Tier 0 = no deps; migrate first.

| Type | Depends on | Tier |
|---|---|---|
| `ContigSpec` | — | 0 |
| `BoltzConfidenceMetrics` | — | 0 |
| `AF3ConfidenceAggregate` | (uses `ProteinStructurePrediction` internally — forward reference; OK to implement in parallel with PSP) | 0 |
| `PipelineParams` | — | 0 |
| `PipelineInternalThresholds` | — | 0 |
| `PositionSet` | `ContigSpec` (optional) | 1 |
| `ProteinStructurePrediction` | `PositionSet` (return type for contact methods) | 2 |
| `DesignedBackbone` | `ProteinStructurePrediction`, `ContigSpec`, `PositionSet` | 3 |
| `DesignedSequence` | `DesignedBackbone`, `PositionSet` | 4 |
| `StageResult` | `ProteinStructurePrediction`, `BoltzConfidenceMetrics`, `PositionSet`, `PipelineInternalThresholds` | 4 |
| `OrthogonalMetrics` | `ProteinStructurePrediction`, `AF3ConfidenceAggregate`, `PipelineInternalThresholds` | 4 |
| `NegativeSteeringRun` | `StageResult`, `DesignedSequence`, `PipelineInternalThresholds` | 5 |
| `DesignCohort` | `NegativeSteeringRun`, `OrthogonalMetrics`, `PipelineInternalThresholds` | 6 |

**Migration tier order:** 0 → 1 → 2 → 3 → 4 → 5 → 6.  Within a tier, types are independent and can land in any order.

### Data flow through one pipeline run

1. **Entry.**  `PipelineParams.from_nextflow_json(...)` validates user input; pipeline aborts on error before any compute.  `PipelineInternalThresholds.default()` loaded alongside.
2. **Branch A (HADDOCK) or Branch B (pre-docked).**  Input PDB(s) become a `ProteinStructurePrediction` (the input complex).
3. **RFDiffusion.**  Produces N `DesignedBackbone`s (each carrying its `ContigSpec`).  Rosetta filter populates `sc_score` and culls those below threshold.
4. **ProteinMPNN.**  For each surviving `DesignedBackbone`, produces M `DesignedSequence`s.  QC + cluster + top-N pick.
5. **Negative steering — per `DesignedSequence`.**  Driven by `bin/negative_steering_run_one.sh` (orchestration unchanged per Q10; the read-side types just consume its workdir at the end).
   - Cold-start: predict the wild-type baseline `num_seeds` times → `StageResult` (cold_start).
   - Check `cold_start.triggers_next_stage(thresholds)` (existing majority-clean rule).  If False (= majority clean): skip steering; the run is `outcome="no_reversion"` with `n_pass = n_seeds`.
   - Else: for each steered design slot, predict `num_seeds` times → `StageResult` (steered).
   - For each steered `StageResult`, check `triggers_next_stage(thresholds, contamination_positions)` (the NEW CL-3 majority rule).  If True: queue reversion for this design.
   - Reversion: for each unique reverted sequence, predict `num_seeds` times → `StageResult` (reversion).
   - Construct `NegativeSteeringRun.from_workdir(workdir)` once all stages are done.
6. **Cohort aggregation.**  `DesignCohort.from_runs_directory(runs_dir)` builds the cohort.  Tier-ranking, composite-score ranking, survivor selection.  `cross_sequence_summary.csv` emitted.
7. **Orthogonal validation — per survivor.**  For each `survivor` in `cohort.survivors(thresholds)`: run AF3, biophysical, Rosetta against `survivor.canonical_prediction(thresholds)`.  Build `OrthogonalMetrics.from_orthogonal_outputs(workdir, survivor.mpnn_sequence_id)`.
8. **Cohort summary plot.**  Reads `DesignCohort` + per-survivor `OrthogonalMetrics`; renders the cohort-summary table with the three-axis column groups.

The whole pipeline is linear — no loops, no shared mutable state between stages.  Each type is constructed once from already-produced outputs and is read-only thereafter.

## Step 4 — Code-fit validation

Four representative pieces of current code rewritten on paper against the new types.  The point is to surface gaps in the spec BEFORE migration starts, not to ship the rewrites here.

### 4.1 Contamination check + reversion gating (the CL-3 fix)

**Current code.**  `bin/boltz2_iterate_steering.py:cmd_build_contaminated` walks each prefilter-passing (design, seed) candidate; for each candidate that's intact AND ra_eff < threshold AND has at least one mutated position contacting effector, it appends a contaminated entry.  Result: per-(design, seed) gating — even one contaminated seed triggers reversion for the design.  This is the bug CL-3 fixes.

**Rewrite using the new types.**

```python
# After Phase 4 migration, this becomes a method on NegativeSteeringRun
# during construction, OR a CLI subcommand that takes a workdir and
# produces contaminated.json the same way as today.

def designs_needing_reversion(
    self: NegativeSteeringRun,
    thresholds: PipelineInternalThresholds,
    contamination_positions: PositionSet,  # design_region ∪ true_interface
) -> list[str]:
    """Returns the design IDs (steered slots) that should have reversion run.
    Applies CL-3: per-design majority rule, not per-seed."""
    needing = []
    for design_id, steered_stage in self.steered.items():
        if steered_stage.triggers_next_stage(thresholds, contamination_positions):
            needing.append(design_id)
    return needing

# StageResult.triggers_next_stage for stage_type == "steered" is:
def triggers_next_stage(self, thresholds, contamination_positions):
    if self.stage_type == "steered":
        n_correctly_placed = self.n_correctly_placed(thresholds)
        n_contaminated = self.n_contaminated(thresholds, contamination_positions)
        if n_correctly_placed == 0:
            return False  # nothing to revert
        return n_contaminated >= math.ceil(n_correctly_placed / 2)
    ...
```

**Discoveries (gaps surfaced).**
- ✅ `contamination_positions` is a `PositionSet` — needs to come from somewhere.  In current code it's built from `cycle_0/plan.json` (`design_region ∪ true_interface`).  Means: `NegativeSteeringRun` needs a `contamination_positions(self)` method that loads + caches this from the workdir.  Add to §2.9.
- ✅ `StageResult.n_correctly_placed` and `n_contaminated` need to access the per-seed metrics — those live in the `confidences` list and the `predictions` list (heavy-atom-contacts come from `ProteinStructurePrediction.receptor_contact_residues(cutoff)`).  Cutoff comes from `thresholds.contact_cutoff`.  Spec already supports this.
- ⚠ The "mutated positions" needed by `n_contaminated` come from `DesignedSequence.mutations_vs_native()` — but the StageResult doesn't directly know its DesignedSequence.  Fix: `NegativeSteeringRun.designs_needing_reversion(...)` knows both, and either passes mutated_positions into `triggers_next_stage`, or `StageResult` carries a `mutated_positions: PositionSet | None` field at construction (per current §2.8 spec).

**Conclusion.**  Spec covers this with one tiny addition (loading `contamination_positions` once from the workdir into `NegativeSteeringRun`).  CL-3 rule implementable cleanly.

### 4.2 Cross-sequence aggregation (`cross_sequence_summary.aggregate` rewrite)

**Current code.**  `bin/cross_sequence_summary.py:aggregate` reads every per-sequence `passing_summary.csv`, picks a representative row per sequence (tier-then-composite), cross-ranks across sequences, writes `cross_sequence_summary.csv` with `rep_*` columns.  ~200 LOC.

**Rewrite using the new types.**

```python
def write_cross_summary(
    runs_dir: Path,
    output: Path,
    thresholds: PipelineInternalThresholds,
    scored_metadata: dict | None = None,
):
    cohort = DesignCohort.from_runs_directory(runs_dir)
    if scored_metadata:
        cohort.attach_designed_sequences(scored_metadata)
    cohort.to_cross_summary_csv(output, thresholds)

# DesignCohort.to_cross_summary_csv is:
def to_cross_summary_csv(self, path, thresholds):
    ranked = self.ranked_by_composite(thresholds)
    rows = []
    for cross_rank, run in enumerate(ranked, start=1):
        row = {
            "mpnn_sequence": run.mpnn_sequence_id,
            "row_type": run.row_type,
            "cross_rank_by_composite": cross_rank,
            "cross_tier": run.tier(thresholds),
            "cross_composite_score": run.composite_score(thresholds),
            ...
            "rep_canonical_pdb": str(run.canonical_prediction(thresholds).path),
            "rep_ra_eff_vs_truth_median": run.median_ra_eff(thresholds),
            ...
        }
        if run.designed_sequence is not None:
            row["corrected_receptor"] = run.designed_sequence.corrected_receptor
            row["designed_residues"] = run.designed_sequence.designed_residues
            row["native_residues"] = run.designed_sequence.native_residues
        rows.append(row)
    csv.DictWriter(...).writerows(rows)
```

**Discoveries.**
- ✅ The `rep_*` column names map cleanly to method calls on `NegativeSteeringRun` (`run.canonical_prediction(...)`, `run.median_ra_eff(...)`, etc.).  Spec covers it.
- ⚠ `NegativeSteeringRun` needs convenience methods like `median_ra_eff`, `median_iptm`, etc. for every column that ends up in the rep_* group.  Either: (a) a wide explicit set of methods, (b) a single `aggregated_metric(name, thresholds)` method that takes a metric name.  Option (a) is more discoverable but verbose; (b) is dynamic but less safe.  Defer to first migration commit.
- ✅ The tier-none fallback (`_pick_rep_from_aggregated` in current code) becomes part of `NegativeSteeringRun.canonical_prediction(thresholds)` — it knows what to do when passing_summary.csv is empty.
- ✅ Sort by `(tier_order, -composite)` is `DesignCohort.ranked_by_composite(thresholds)`.

**Conclusion.**  Spec covers this.  One implementation choice (wide methods vs metric-name method) to make at migration time.

### 4.3 Per-seed cross-stage verdict aggregation

**Current code.**  `bin/boltz2_iterate_steering.py:_per_seed_verdict_breakdown` + `_classify_outcome` together produce the final per-seed verdict tally for one MPNN sequence (used to compute `n_pass` and the aggregate `outcome` label).  Branches on whether the reversion verdict column is populated for each seed.

**Rewrite using the new types.**

```python
# Method on NegativeSteeringRun:
def per_seed_final_verdicts(
    self, thresholds, contamination_positions
) -> list[dict]:
    """Returns one entry per (design, seed). If reversion ran for that
    design, use the reversion verdict; otherwise use the steered verdict."""
    entries = []

    # Special case: cold-start path (no steering)
    if not self.steered:  # cold-start was sufficient
        cold_verdicts = self.cold_start.per_seed_verdicts(thresholds)
        for seed_idx, verdict in zip(self.cold_start.seed_indices, cold_verdicts):
            entries.append({
                "design_id": "initial",
                "seed_index": seed_idx,
                "stage_run": "cold_start",
                "verdict": verdict,  # "clean" or "fail_structural"
            })
        return entries

    # Normal case: steered + maybe reversion
    for design_id, steered_stage in self.steered.items():
        steered_verdicts = steered_stage.per_seed_verdicts(
            thresholds, contamination_positions
        )
        # Does this design's reversion exist?
        rev_stage = self._reversion_for_design(design_id)  # may be None
        if rev_stage is None:
            # No reversion for this design — steered verdicts stand
            for seed_idx, sv in zip(steered_stage.seed_indices, steered_verdicts):
                entries.append({"design_id": design_id, "seed_index": seed_idx,
                                "stage_run": "steered", "verdict": sv})
        else:
            rev_verdicts = rev_stage.per_seed_verdicts(
                thresholds, contamination_positions
            )
            # Reverted verdicts override steered for matching seed_index
            for seed_idx, rv in zip(rev_stage.seed_indices, rev_verdicts):
                entries.append({"design_id": design_id, "seed_index": seed_idx,
                                "stage_run": "reverted", "verdict": rv})
    return entries

def n_pass(self, thresholds, contamination_positions) -> int:
    return sum(1 for e in self.per_seed_final_verdicts(thresholds, contamination_positions)
               if e["verdict"] in ("clean_steered", "pose_holds", "clean"))
```

**Discoveries.**
- ⚠ The mapping from a steered design to its reversion is not obvious.  Current code uses a label hash to dedupe by reverted sequence.  Need a `NegativeSteeringRun._reversion_for_design(design_id)` helper that walks the reversion stages and finds the one whose `mutated_positions` corresponds to this design's reversions.  Tracked: add this private method to §2.9.
- ⚠ When a design's reversion deduplicates with another design (same reverted sequence), the same reversion stage is referenced by multiple designs.  The current code maps via `all_label_seeds` in reversion_metadata.json — needs to carry through to the rewritten type.  Spec needs: `StageResult` for reversion carries `applies_to_designs: list[str]`.  Add to §2.8.
- ✅ The "verdict overrides steered" logic is clean once the dedup mapping is right.

**Conclusion.**  Two small spec additions: `_reversion_for_design` on NegativeSteeringRun and `applies_to_designs` field on reversion-typed StageResult.

### 4.4 Orthogonal filter check (`merge_orthogonal_metrics._apply_filters` rewrite)

**Current code.**  `bin/merge_orthogonal_metrics.py:_apply_filters` reads a survivor's row, checks `sc >= 0.55 AND bsa >= 600 AND ddg <= reference`, returns flags list.  AF3 demoted to flag-only (post Phase 5 Q118 fix).

**Rewrite using the new types.**

```python
def write_survivor_orthogonal_metrics(
    survivors: list[NegativeSteeringRun],
    runs_dir: Path,
    output_csv: Path,
    thresholds: PipelineInternalThresholds,
):
    rows = []
    for survivor in survivors:
        om = OrthogonalMetrics.from_orthogonal_outputs(
            runs_dir / survivor.mpnn_sequence_id, survivor.mpnn_sequence_id
        )
        passes = om.passes_orthogonal_filters(thresholds)
        failed = om.failed_filter_names(thresholds)
        af3_disagree = om.af3_disagrees(thresholds)
        rows.append({
            "mpnn_sequence": survivor.mpnn_sequence_id,
            "passes_orthogonal_filters": passes,
            "failed_filters": ",".join(failed),
            "af3_disagrees": af3_disagree,
            **om.to_summary_row(),
        })
    csv.DictWriter(...).writerows(rows)
```

**Discoveries.**
- ✅ Clean.  `OrthogonalMetrics` already owns the gate; the caller just iterates over survivors.
- ✅ The test/production divergence (the long-standing `merge_orthogonal_metrics.py` problem) disappears because there's only one home for the gate logic — it's a method on `OrthogonalMetrics`, not a script that two copies of can drift.

**Conclusion.**  Cleanest rewrite of the four.  Spec covers this without modification.

### Step 4 summary — spec additions needed

The four code-fit validations surfaced three small spec additions:

1. **§2.8 `StageResult`** — for stage_type=`reversion`, carry an `applies_to_designs: list[str]` field that records which steered design_ids' contamination produced this unique reverted sequence.
2. **§2.9 `NegativeSteeringRun`** — add `contamination_positions(self) -> PositionSet` (cached load of design_region ∪ true_interface from the workdir's plan.json).
3. **§2.9 `NegativeSteeringRun`** — add private `_reversion_for_design(design_id) -> StageResult | None` helper that uses `applies_to_designs` to find the right reversion stage for a given design.

No major structural changes needed.  Spec survives contact with reality.

---

## Phase 4 implementation notes (post-spec)

Once Steps 1–4 are locked, the migration follows the cadence agreed in Session 6 (Option 3: design-first, then module-by-module migration).

Implementation tasks (to be created as TaskCreate entries when migration begins):

1. Foundation library: `ContigSpec`, `PositionSet`, `BoltzConfidenceMetrics`, `AF3ConfidenceAggregate` — small types, no orchestration coupling.
2. Structural type: `ProteinStructurePrediction` — consolidates all chain/Cα/contact code.
3. Design types: `DesignedBackbone`, `DesignedSequence` — RFDiffusion + MPNN stage consolidation.
4. Negsteer types: `StageResult`, `NegativeSteeringRun` — the read-side rewrite, includes CL-3 majority-rule fix to `cmd_build_contaminated`.
5. Cohort + orthogonal: `DesignCohort`, `OrthogonalMetrics` — replaces the CSV-and-dict aggregation in `cross_sequence_summary.py` and `merge_orthogonal_metrics.py`.  Test/production divergence eliminated.
6. Params + thresholds: `PipelineParams`, `PipelineInternalThresholds` — wraps `validate_params.py`, collects scattered literals into the threshold catalogue.
7. Test/prod plot-script consolidation — separate but related cleanup.

Each step gets its own commit, with characterization tests re-run between commits (per the workflow established in `tests/run_tests.sh` and `tests/update_example_dataset.slurm.sh`).
