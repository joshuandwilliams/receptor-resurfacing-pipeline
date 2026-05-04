# Literature sweep and structure analysis — Part 2: Synthesis

You are conducting the synthesis phase of a structured literature, PDB, and UniProt sweep to inform the design of a binder-resurfacing campaign. Part 1 produced a fetch manifest at `experiments/campaigns/pikp1_pwt7/inputs/pikp1_pwt7_fetch_manifest.md` listing papers that could not be accessed online. The user has now placed the requested PDFs into `experiments/inputs/context/`. Your task here is to read those new PDFs, combine them with everything you already fetched online in Part 1, and write the full literature-sweep report.

This report is the curated input to a downstream pipeline. Its quality determines whether the campaign succeeds. Be thorough, be specific, and surface uncertainty rather than hide it.

---

## Inputs

**Receptor**: Pikp-1 HMA
**Receptor UniProt**: E9KPB5
**Receptor domain boundaries**: 186-263

**Effector**: PWT7

**Prior knowledge to take into account**:
- All deposited Pikp-1 HMA crystal structures truncate or disorder C-terminal residues critical for binding. Use rung 5 of the structural-input fallback ladder — an AlphaFold3 monomer prediction of native Pikp-1 HMA, residues 186-263 — rather than dropping to an engineered crystal. This receptor-side decision is shared across every Pikp-1-HMA-binder campaign, including this one.
- The AF3 prediction is renumbered from 1 (position 1 = Pikp-1 residue 186; position 78 = residue 263). The contig must be written in this 1-78 numbering.
- The user has indicated PWT7 binds Pikp-HMA at a different face than the canonical Pik-HMA / AVR-Pik face. Verify this in the literature; if confirmed, the structural-input decision and contig design must follow a different geometry from the pikp1_avrpikf campaign — specifically, the alignment template will not be 6G10 (which positions the effector on the canonical AVR-Pik face), and the flex / anchor regions will need to track whichever residues actually contact PWT7.
- PWT7 is expected to be a MAX-fold effector based on the project context. Complex assembly should follow the ChimeraX-alignment recipe used for other pikp1_* campaigns, but with the alignment template chosen to reflect the actual PWT7 binding face.
- The Pikp-1 HMA binder is shared across this campaign and the other pikp1_* campaigns; the receptor-side anchor / flex-region analysis from those campaigns is reusable here, but the binding-face determination is campaign-specific and may move which Pikp-HMA residues are flex vs anchor.

---

## What to read first

1. The fetch manifest at `experiments/campaigns/pikp1_pwt7/inputs/pikp1_pwt7_fetch_manifest.md` — header records the recommended source-PDB and numbering convention plus the alignment-template choice for the PWT7 binding face, carried forward from Part 1.
2. `experiments/inputs/context/` — list its contents. For each gap-listed paper from the manifest, find the matching PDF (filename convention `<firstauthor>_<year>_<shortjournal>.pdf`) and read it. If a gap-listed paper is missing from the context folder, note it and proceed without it; do not block.
3. Any sibling pikp1_* campaign's literature-sweep report (in `experiments/campaigns/pikp1_*/inputs/*_literature_sweep.md`) — the receptor-side analysis is shared, and citing it explicitly avoids re-deriving the same conclusions.

---

## Required output structure

Produce a Markdown document with exactly the following six top-level sections, in this order. Use the section names verbatim. Do not add extra top-level sections.

### 1. Summary

A short paragraph (3–6 sentences) stating: what receptor and effector this report covers, the goal of the campaign, and the headline findings — which structures are usable, what the structural-input plan is (acknowledging the non-canonical binding face), and which receptor regions are the primary redesign targets.

### 2. Available structures

Catalogue every relevant structure you find. Use Markdown tables. Separate tables for:

- Native receptor structures (no engineering mutations).
- Engineered or variant receptor structures (with mutations) — note the mutations explicitly.
- Effector structures (alone or in complex with anything).
- Any complexes between the receptor and any effector (target or related) — used as alignment templates. **Be especially careful to identify which complex template represents the PWT7 binding face, since this determines the contig.**

For each entry, record at minimum: PDB ID, complex contents, domain boundaries, resolved ATOM range, key mutations, and the publication.

### 3. Structural-input decision

Lay out the decision about what structural input the downstream pipeline will use, treating it as a **fallback ladder**. State explicitly **which rung this campaign lands on, and why**.

Specify the alignment-template recipe explicitly:

- Which complex template will be used to position the AF3-predicted Pikp-1 HMA against PWT7? Justify the choice — the template must reflect the actual PWT7 binding face, not the canonical AVR-Pik face.
- If no published Pik-HMA / PWT7 complex exists, state which closest-related complex provides the alignment frame.
- If the binding face cannot be confidently determined from the available literature, fall back to AF3-multimer of the complex and flag this as a soft blocker.

### 4. Binding interfaces

Decompose the binding interface into named sub-interfaces based on the structural and mutagenesis literature for PWT7 and any closest-related Pik-HMA / MAX-effector complex. For each sub-interface:

- A short structural description (which secondary-structure elements, approximate residue range).
- **Receptor-side residues** that participate, listed individually with their identities, citations, and the specific role each plays. **Note explicitly whether these are the same residues as in the canonical AVR-Pik face or different ones.**
- **Effector-side residues** that participate.
- **Specificity role**: is this interface conserved across related effector variants or specificity-determining?
- **Design implication**: should this region be left fixed, redesigned as a primary target, or redesigned as a secondary target?

If the literature gives a specific reason why the native Pikp-1 fails to recognise PWT7 (a polymorphism, a steric clash, a missing contact at the alternative face), state it explicitly.

### 5. Numbering translation

Translate between numbering conventions explicitly. Provide:

- The numbering convention used in the literature for the receptor and effector.
- The numbering convention the pipeline will use.
- A translation table mapping each key residue from literature numbering to pipeline numbering, with the offset stated clearly.

### 6. Contig string design

Produce an RFDiffusion contig string for this campaign, using **pipeline numbering** (not literature numbering). The contig syntax is:

```
A<start>-<end>/<min>-<max>/A<start>-<end>/... B
```

**Design principle: the flexible regions must contain the residues to be redesigned.** Do not anchor on or copy any contig string from a sibling pikp1_* campaign — the contact pattern in *this* campaign's complex (which uses a different binding face) determines the flex / anchor regions, and is expected to be substantially different from the canonical AVR-Pik face used in pikp1_avrpikf.

Provide:

- A **primary contig** representing the recommended design, derived from the contact pattern at the PWT7 binding face.
- One or more **alternative contigs** representing more conservative or more aggressive variants.
- For each contig, a breakdown of what each segment corresponds to structurally.

Flag any boundary you propose that you suspect may fall mid-secondary-structure.

---

## Quality criteria

- **Be specific**. Cite PDB IDs, residue numbers, paper authors and years.
- **Surface uncertainty**. The non-canonical binding face means much of the analysis depends on whether the Part-1 binding-face determination was confirmed. Mark inferred claims clearly.
- **Distinguish fact from recommendation**.
- **Numbering must be consistent**.
- **Flag campaign-blockers prominently**. If the binding face cannot be determined, the campaign cannot proceed without further structural input.

## Output format

Deliver the report as a single Markdown document at `experiments/campaigns/pikp1_pwt7/inputs/pikp1_pwt7_literature_sweep.md`. Use the section structure prescribed above. Use Markdown tables for structure inventories and renumbering translations. Use code blocks for contig strings. Cite references inline by author + year on first mention, with a final References section listing each citation in full at the end of the document.

Also write a sibling provenance file `pikp1_pwt7_literature_sweep.notes` recording: the date the sweep was completed, the source-PDB and numbering-convention determinations, the names of the primary papers cited, and which papers came from the context folder vs online.

When the report is written, delete the fetch manifest. Summarise back in chat: which fallback-ladder rung the campaign lands on, the numbering convention, the primary contig, any campaign-blockers, and any gap-listed papers that were missing from the context folder.

## What not to do

- Do not modify any input complex PDB, its `.notes` file, an existing contig, or any pipeline scripts.
- Do not modify or move any PDFs in `experiments/inputs/context/`.
- Do not run the pipeline, derive new contigs from scripts, or submit anything to HPC.
- Do not commit or push.
- Do not produce a v2 alongside v1 or append a "changes from context/" section.
- Do not derive contig strings from any pre-existing contig file. Build them from the §4 binding-interface analysis.
