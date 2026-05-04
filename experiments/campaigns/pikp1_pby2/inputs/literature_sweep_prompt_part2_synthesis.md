# Literature sweep and structure analysis — Part 2: Synthesis

You are conducting the synthesis phase of a structured literature, PDB, and UniProt sweep to inform the design of a binder-resurfacing campaign. Part 1 produced a fetch manifest at `experiments/campaigns/pikp1_pby2/inputs/pikp1_pby2_fetch_manifest.md` listing papers that could not be accessed online. The user has now placed the requested PDFs into `experiments/inputs/context/`. Your task here is to read those new PDFs, combine them with everything you already fetched online in Part 1, and write the full literature-sweep report.

This report is the curated input to a downstream pipeline. Its quality determines whether the campaign succeeds. Be thorough, be specific, and surface uncertainty rather than hide it.

---

## Inputs

**Receptor**: Pikp-1 HMA
**Receptor UniProt**: E9KPB5
**Receptor domain boundaries**: 186-263

**Effector**: PBY2

**Prior knowledge to take into account**:
- All deposited Pikp-1 HMA crystal structures truncate or disorder C-terminal residues critical for binding. Use rung 5 of the structural-input fallback ladder — an AlphaFold3 monomer prediction of native Pikp-1 HMA, residues 186-263 — rather than dropping to an engineered crystal. This receptor-side decision is shared across every Pikp-1-HMA-binder campaign, including this one.
- The AF3 prediction is renumbered from 1 (position 1 = Pikp-1 residue 186; position 78 = residue 263). The contig must be written in this 1-78 numbering.
- The user has only limited information about PBY2; identifying the effector (full name, organism, MAX-fold family membership, any deposited PDB / UniProt records) is part of the discovery work in Part 1. PBY2 is expected to be a MAX-fold effector based on the project context.
- If a Pik-HMA / PBY2 (or a related Pik-HMA / MAX-effector) complex template exists, complex assembly should follow the same ChimeraX-alignment recipe used for other pikp1_* campaigns. If no such template exists, fall back to AF3-multimer for the assembly.
- The Pikp-1 HMA binder is shared across this campaign and the other pikp1_* campaigns; the receptor-side anchor / flex-region analysis from those campaigns is reusable here.

---

## What to read first

1. The fetch manifest at `experiments/campaigns/pikp1_pby2/inputs/pikp1_pby2_fetch_manifest.md` — header records source-PDB and numbering-convention determinations (or recommendations) carried forward from Part 1.
2. `experiments/inputs/context/` — list its contents. For each gap-listed paper from the manifest, find the matching PDF (filename convention `<firstauthor>_<year>_<shortjournal>.pdf`) and read it. If a gap-listed paper is missing from the context folder, note it and proceed without it; do not block.
3. Any sibling pikp1_* campaign's literature-sweep report (in `experiments/campaigns/pikp1_*/inputs/*_literature_sweep.md`) — the receptor-side analysis is shared, and citing it explicitly avoids re-deriving the same conclusions.

---

## Required output structure

Produce a Markdown document with exactly the following six top-level sections, in this order. Use the section names verbatim. Do not add extra top-level sections; if you have additional content, fit it as a sub-section under one of these.

### 1. Summary

A short paragraph (3–6 sentences) stating: what receptor and effector this report covers, the goal of the campaign (redesign the receptor's binding surface to gain or improve recognition of the target effector), and the headline findings — which structures are usable, what the structural-input plan is, and which receptor regions are the primary redesign targets.

### 2. Available structures

Catalogue every relevant structure you find. Use Markdown tables. Separate tables for:

- Native receptor structures (no engineering mutations relative to the target protein).
- Engineered or variant receptor structures (with mutations) — note the mutations explicitly.
- Effector structures (alone or in complex with anything).
- Any complexes between the receptor and any effector (target or related).

For each entry, record at minimum: PDB ID, complex contents (what chains are what), domain boundaries, resolved ATOM range if it differs from SEQRES, key mutations if any, and the publication. If a structure has known issues (disordered terminal residues, missing loops, non-physiological constructs), say so in a notes column.

### 3. Structural-input decision

This is the most important section. Lay out the decision about what structural input the downstream pipeline will use, treating it as a **fallback ladder**. The rungs, in preference order:

1. A native receptor / target-effector complex crystal — directly usable as the input PDB.
2. A native receptor / related-effector complex crystal — usable if the related effector can be substituted for the target by structural alignment.
3. A native receptor structure (alone) plus a target-effector structure (alone) — usable if the two can be aligned onto a related complex template to generate a positioned input PDB.
4. An engineered receptor crystal — usable only if the mutations are not in the binding interface; flag explicitly which mutations are present and whether they affect the design region.
5. An AlphaFold3 monomer prediction of the native receptor — used when no native crystal exists or all native crystals have material problems (disorder in the design region, wrong domain boundaries, etc.).
6. An AlphaFold3 monomer prediction of the effector — used when no effector crystal exists.
7. An AlphaFold3-multimer prediction of the complex — used as a last resort when no complex template exists at all.

State explicitly **which rung this campaign lands on, and why**. If a higher rung is theoretically available but unusable, say what disqualifies it.

**If no usable complex template exists at all**, flag this prominently. The downstream workflow assembles the input PDB by aligning structures onto a complex template; if no complex template exists, the campaign needs a different assembly strategy (likely AF3-multimer of the complex).

If the chosen rung requires AF3 prediction, specify the input: which sequence, which boundaries, and what validation should be done on the predicted structure.

### 4. Binding interfaces

Decompose the binding interface into named sub-interfaces based on the structural and mutagenesis literature. For each sub-interface:

- A short structural description (which secondary-structure elements, approximate residue range).
- **Receptor-side residues** that participate, listed individually with their identities, citations, and the specific role each plays.
- **Effector-side residues** that participate, with the same level of detail, and noting any polymorphisms across effector variants if relevant.
- **Specificity role**: is this interface conserved across effector variants (structural/scaffold role) or specificity-determining (varies between variants)?
- **Design implication**: should this region be left fixed, redesigned as a primary target, or redesigned as a secondary target? Justify with reference to the literature.

If the literature gives a specific reason why the target effector evades the native receptor (a polymorphism, a steric clash, a missing contact), state it explicitly. This is the central insight that drives which interface is the primary redesign target.

### 5. Numbering translation

Translate between numbering conventions explicitly. Provide:

- The numbering convention used in the literature (full-length receptor positions, full-length effector positions).
- The numbering convention the pipeline will use (typically: receptor domain renumbered from 1; effector chain renumbered from 1).
- A translation table mapping each key residue from literature numbering to pipeline numbering, with the offset stated clearly.

### 6. Contig string design

Produce an RFDiffusion contig string for this campaign, using **pipeline numbering** (not literature numbering). The contig syntax is:

```
A<start>-<end>/<min>-<max>/A<start>-<end>/... B
```

Where `A<start>-<end>` is a fixed (anchored) region of the receptor in pipeline numbering, `<min>-<max>` is a flexible (de novo) region with a length range for RFDiffusion sampling, and `B` is the effector chain (entirely fixed).

**Design principle: the flexible regions must contain the residues to be redesigned.** Do not anchor on or copy any contig string from a sibling pikp1_* campaign — the contact pattern in *this* campaign's complex determines the flex / anchor regions.

Provide:

- A **primary contig** representing the recommended design.
- One or more **alternative contigs** representing more conservative or more aggressive variants.
- For each contig, a breakdown of what each segment corresponds to structurally.

Flag any boundary you propose that you suspect may fall mid-secondary-structure.

---

## Quality criteria

- **Be specific**. Cite PDB IDs, residue numbers, paper authors and years.
- **Surface uncertainty**. If the literature is conflicting, sparse, or absent, say so.
- **Distinguish fact from recommendation**. Sections 2, 4, and 5 are factual; sections 3 and 6 are recommendations.
- **Numbering must be consistent**. Every residue number must be unambiguous.
- **Flag campaign-blockers prominently**.

## Output format

Deliver the report as a single Markdown document at `experiments/campaigns/pikp1_pby2/inputs/pikp1_pby2_literature_sweep.md`. Use the section structure prescribed above. Use Markdown tables for structure inventories and renumbering translations. Use code blocks for contig strings. Cite references inline by author + year on first mention, with a final References section listing each citation in full at the end of the document.

Also write a sibling provenance file `pikp1_pby2_literature_sweep.notes` recording: the date the sweep was completed, the source-PDB and numbering-convention determinations, the names of the primary papers cited, and which papers came from the context folder vs online.

When the report is written, delete the fetch manifest. Summarise back in chat: which fallback-ladder rung the campaign lands on, the numbering convention, the primary contig, any campaign-blockers, and any gap-listed papers that were missing from the context folder.

## What not to do

- Do not modify any input complex PDB, its `.notes` file, an existing contig, or any pipeline scripts.
- Do not modify or move any PDFs in `experiments/inputs/context/`.
- Do not run the pipeline, derive new contigs from scripts, or submit anything to HPC.
- Do not commit or push.
- Do not produce a v2 alongside v1 or append a "changes from context/" section.
- Do not derive contig strings from any pre-existing contig file. Build them from the §4 binding-interface analysis.
