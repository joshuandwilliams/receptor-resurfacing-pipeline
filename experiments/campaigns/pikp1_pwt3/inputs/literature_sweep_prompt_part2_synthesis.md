# Literature sweep and structure analysis — Part 2: Synthesis

You are conducting the synthesis phase of a structured literature, PDB, and UniProt sweep to inform the design of a binder-resurfacing campaign. Part 1 produced a fetch manifest at `experiments/campaigns/pikp1_pwt3/inputs/pikp1_pwt3_fetch_manifest.md` listing papers that could not be accessed online. The user has now placed the requested PDFs into `experiments/inputs/context/`. Your task here is to read those new PDFs, combine them with everything you already fetched online in Part 1, and write the full literature-sweep report.

This report is the curated input to a downstream pipeline. Its quality determines whether the campaign succeeds. Be thorough, be specific, and surface uncertainty rather than hide it.

---

## Inputs

**Receptor**: Pikp-1 HMA
**Receptor UniProt**: E9KPB5
**Receptor domain boundaries**: 186-263

**Effector**: PWT3

**Prior knowledge to take into account**:
- All deposited Pikp-1 HMA crystal structures truncate or disorder C-terminal residues critical for binding. Use rung 5 of the structural-input fallback ladder — an AlphaFold3 monomer prediction of native Pikp-1 HMA, residues 186-263 — rather than dropping to an engineered crystal. This receptor-side decision is shared across every Pikp-1-HMA-binder campaign, including this one.
- The AF3 prediction is renumbered from 1 (position 1 = Pikp-1 residue 186; position 78 = residue 263). The contig must be written in this 1-78 numbering.
- There is no deposited structure for the PWT3 effector. Chain C of the assembled complex must therefore come from an AlphaFold3 monomer prediction of PWT3 (or, alternatively, an AF3-multimer prediction of the Pikp-HMA / PWT3 complex).
- PWT3 is expected to be a MAX-fold effector based on the project context. If the AF3 monomer prediction is reasonably confident in the MAX fold, complex assembly can proceed by ChimeraX alignment of the predicted PWT3 onto a published Pik-HMA / MAX-effector complex template (e.g. 6G10). If AF3 confidence is poor or the fold prediction is ambiguous, fall back to AF3-multimer.
- The Pikp-1 HMA binder is shared across this campaign and the other pikp1_* campaigns; the receptor-side anchor / flex-region analysis from those campaigns is reusable here.

---

## What to read first

1. The fetch manifest at `experiments/campaigns/pikp1_pwt3/inputs/pikp1_pwt3_fetch_manifest.md` — header records the recommended source-PDB and numbering convention plus the AF3 prediction plan carried forward from Part 1.
2. `experiments/inputs/context/` — list its contents. For each gap-listed paper from the manifest, find the matching PDF (filename convention `<firstauthor>_<year>_<shortjournal>.pdf`) and read it. If a gap-listed paper is missing from the context folder, note it and proceed without it; do not block.
3. Any sibling pikp1_* campaign's literature-sweep report (in `experiments/campaigns/pikp1_*/inputs/*_literature_sweep.md`) — the receptor-side analysis is shared, and citing it explicitly avoids re-deriving the same conclusions.

---

## Required output structure

Produce a Markdown document with exactly the following six top-level sections, in this order. Use the section names verbatim. Do not add extra top-level sections.

### 1. Summary

A short paragraph (3–6 sentences) stating: what receptor and effector this report covers, the goal of the campaign, and the headline findings — which structures are usable, what the structural-input plan is (acknowledging that the effector requires AF3 prediction), and which receptor regions are the primary redesign targets.

### 2. Available structures

Catalogue every relevant structure you find. Use Markdown tables. Separate tables for:

- Native receptor structures (no engineering mutations).
- Engineered or variant receptor structures (with mutations) — note the mutations explicitly.
- Effector structures (alone or in complex with anything) — for this campaign, this section will primarily list related MAX-fold effectors used as templates rather than PWT3 itself.
- Any complexes between the receptor and any effector (target or related) — used as alignment templates.

For each entry, record at minimum: PDB ID, complex contents, domain boundaries, resolved ATOM range if it differs from SEQRES, key mutations if any, and the publication.

### 3. Structural-input decision

Lay out the decision about what structural input the downstream pipeline will use, treating it as a **fallback ladder**. The rungs, in preference order, are listed above.

State explicitly **which rung this campaign lands on, and why**. Since PWT3 has no deposited structure, the chain C source must be either rung 6 (AF3 monomer of effector) or rung 7 (AF3-multimer of complex). Specify:

- The PWT3 sequence to use (UniProt accession + boundary or full precursor).
- Whether to use AF3 monomer (with subsequent ChimeraX alignment to a Pik-HMA / MAX-effector template) or AF3 multimer.
- Validation steps for the AF3 prediction (e.g. compare predicted fold to related MAX effectors; check disulfide bond placement; verify pLDDT in the binding-interface residues).
- Which Pik-HMA / MAX-effector complex template to use for the ChimeraX alignment, and why that template is appropriate (does PWT3 bind the same face as AVR-PikD? — if uncertain, flag the assumption).

### 4. Binding interfaces

Decompose the predicted binding interface into named sub-interfaces based on the literature on related MAX-fold effectors and the AF3 prediction's contact pattern. For each sub-interface:

- A short structural description (which secondary-structure elements, approximate residue range).
- **Receptor-side residues** that participate, listed individually with their identities, citations, and the specific role each plays.
- **Effector-side residues** that participate (predicted from the AF3 model and the homology to known MAX-fold complexes).
- **Specificity role**: is this interface conserved across related effector variants or specificity-determining?
- **Design implication**: should this region be left fixed, redesigned as a primary target, or redesigned as a secondary target?

If the literature gives a specific reason why related Pik-HMA receptors fail to recognise PWT3 (or its closest characterised relative), state it explicitly. Mark inferred claims clearly — for an effector with no deposited structure, much of the residue-level analysis will be inferred from homology and from the AF3 prediction rather than direct mutagenesis evidence.

### 5. Numbering translation

Translate between numbering conventions explicitly. Provide:

- The numbering convention used in the literature for the receptor and any related effector.
- The numbering convention the pipeline will use.
- A translation table mapping each key residue from literature numbering to pipeline numbering, with the offset stated clearly.
- For PWT3 specifically, since the AF3 prediction will renumber from 1, state explicitly the offset between the AF3 chain numbering and the original UniProt / precursor numbering.

### 6. Contig string design

Produce an RFDiffusion contig string for this campaign, using **pipeline numbering** (not literature numbering). The contig syntax is:

```
A<start>-<end>/<min>-<max>/A<start>-<end>/... B
```

**Design principle: the flexible regions must contain the residues to be redesigned.** Do not anchor on or copy any contig string from a sibling pikp1_* campaign — the contact pattern in *this* campaign's complex (especially given chain C comes from an AF3 prediction with its own uncertainty) determines the flex / anchor regions.

Provide:

- A **primary contig** representing the recommended design.
- One or more **alternative contigs** representing more conservative or more aggressive variants.
- For each contig, a breakdown of what each segment corresponds to structurally.

Flag any boundary you propose that you suspect may fall mid-secondary-structure. Acknowledge that an AF3-derived chain C means the contact pattern that drives flex-region selection has uncertainty itself; recommend re-running contact derivation if the assembled complex's binding pose differs materially from the assumed Pik-HMA / MAX-effector template.

---

## Quality criteria

- **Be specific**. Cite PDB IDs, residue numbers, paper authors and years.
- **Surface uncertainty**. For PWT3 specifically, much of the analysis is inferred from homology — be clear about which claims are direct (from a PWT3-specific source) and which are inferred (from related MAX effectors).
- **Distinguish fact from recommendation**. Sections 2, 4, and 5 are factual; sections 3 and 6 are recommendations.
- **Numbering must be consistent**.
- **Flag campaign-blockers prominently**. The lack of a deposited PWT3 structure is itself a soft blocker — the campaign can proceed but with greater uncertainty than for effectors with crystal structures.

## Output format

Deliver the report as a single Markdown document at `experiments/campaigns/pikp1_pwt3/inputs/pikp1_pwt3_literature_sweep.md`. Use the section structure prescribed above. Use Markdown tables for structure inventories and renumbering translations. Use code blocks for contig strings. Cite references inline by author + year on first mention, with a final References section listing each citation in full at the end of the document.

Also write a sibling provenance file `pikp1_pwt3_literature_sweep.notes` recording: the date the sweep was completed, the source-PDB and numbering-convention determinations, the names of the primary papers cited, and which papers came from the context folder vs online.

When the report is written, delete the fetch manifest. Summarise back in chat: which fallback-ladder rung the campaign lands on, the numbering convention, the primary contig, any campaign-blockers, and any gap-listed papers that were missing from the context folder.

## What not to do

- Do not modify any input complex PDB, its `.notes` file, an existing contig, or any pipeline scripts.
- Do not modify or move any PDFs in `experiments/inputs/context/`.
- Do not run the pipeline, derive new contigs from scripts, or submit anything to HPC.
- Do not commit or push.
- Do not produce a v2 alongside v1 or append a "changes from context/" section.
- Do not derive contig strings from any pre-existing contig file. Build them from the §4 binding-interface analysis.
