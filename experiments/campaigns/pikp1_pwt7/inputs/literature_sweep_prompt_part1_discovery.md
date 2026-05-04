# Literature sweep and structure analysis — Part 1: Discovery

You are conducting the discovery phase of a structured literature, PDB, and UniProt sweep to inform the design of a binder-resurfacing campaign. The eventual report (produced in Part 2) will be used to (a) decide which structures to use as input to a structural alignment, (b) curate which residues should be redesigned vs preserved, and (c) define an RFDiffusion contig string that encodes those decisions.

This part is **discovery only**. You will identify every source you need, attempt to fetch each one, and produce a short gap-list of papers you could not access. **Do not write the report yet** — that is Part 2.

---

## Inputs

**Receptor**: Pikp-1 HMA
**Receptor UniProt**: E9KPB5
**Receptor domain boundaries**: 186-263

**Effector**: PWT7

**Prior knowledge to take into account**:
- All deposited Pikp-1 HMA crystal structures truncate or disorder C-terminal residues critical for binding. Use rung 5 of the structural-input fallback ladder — an AlphaFold3 monomer prediction of native Pikp-1 HMA, residues 186-263 — rather than dropping to an engineered crystal. This receptor-side decision is shared across every Pikp-1-HMA-binder campaign, including this one.
- The AF3 prediction is renumbered from 1 (position 1 = Pikp-1 residue 186; position 78 = residue 263). The contig must be written in this 1-78 numbering.
- The user has indicated PWT7 binds Pikp-HMA at a different face than the canonical Pik-HMA / AVR-Pik face. Verify this in the literature during Part 1; if confirmed, the structural-input decision and contig design must follow a different geometry from the pikp1_avrpikf campaign — specifically, the alignment template will not be 6G10 (which positions the effector on the canonical AVR-Pik face), and the flex / anchor regions will need to track whichever residues actually contact PWT7.
- PWT7 is expected to be a MAX-fold effector based on the project context. Complex assembly should follow the ChimeraX-alignment recipe used for other pikp1_* campaigns, but with the alignment template chosen to reflect the actual PWT7 binding face.
- Identifying the effector (full name, organism, family membership, any deposited PDB / UniProt records) is part of the discovery work in Part 1; the user has limited information.
- The Pikp-1 HMA binder is shared across this campaign and the other pikp1_* campaigns; the receptor-side anchor / flex-region analysis from those campaigns is reusable here, but the binding-face determination is campaign-specific and may move which Pikp-HMA residues are flex vs anchor.

---

## Sources to consult

Search across all of the following. Cite specific papers, PDB entries, and UniProt records throughout the manifest you produce.

1. **Plant protein structures CSV** at `experiments/inputs/plant_protein_structures.csv`. Read this first as a starting reference for what structures already exist for this receptor and effector. Search PDB independently for any PWT7-related entries the CSV may have missed.

2. **Scientific literature**. Identify the primary structural and functional papers for both the receptor and the effector. Pay particular attention to:
   - Structural studies of PWT7 (or any closely-related effector) bound to a Pik-HMA or any other rice HMA — needed to confirm the "different binding face" prior knowledge.
   - Mutagenesis studies that map specific PWT7 residues to binding affinity, specificity, or virulence.
   - Survey / review papers cataloguing MAX-fold effectors and the rice-blast effector arsenal.
   - Engineering studies that have already attempted to redirect Pik-HMA recognition toward a non-canonical binding face.

3. **Protein Data Bank**. For every relevant structure, you will eventually need: PDB ID, the complex contents, domain boundaries, the resolved ATOM range, resolution, and the publication.

4. **UniProt**. Canonical sequences for the receptor and effector. Domain annotations, signal peptides, post-translational modifications, feature annotations relevant to binding.

---

## Campaign context to read first

Before scoping the sweep, read:

- `experiments/README.md` — campaign-lifecycle and naming conventions.
- The campaign's own README at `experiments/campaigns/pikp1_pwt7/README.md` if present.
- Any existing complex PDB at `experiments/campaigns/pikp1_pwt7/inputs/` (if not present, note that complex assembly is an outcome of this sweep, not an input).
- Any existing contig file in the campaign's `inputs/` directory (if present).
- Outputs of the parallel pikp1_* campaigns' literature sweeps in their `inputs/` directories — the Pikp-1 HMA binder-side analysis is shared.
- `experiments/inputs/context/` — the shared paper PDF folder. List its contents now. Filename convention is `<firstauthor>_<year>_<shortjournal>.pdf`. You will check this folder for any paper you can't fetch online.

---

## What to do

1. Read all the campaign context above.
2. Plan the literature/PDB/UniProt sweep. Enumerate every paper, PDB entry, and UniProt record you intend to consult. **Be thorough.**
3. For each source, attempt to fetch it. Use web search to find it, web fetch to retrieve it. Prefer open-access mirrors over publisher sites.
4. For each source you cannot fetch online:
   - Check `experiments/inputs/context/` for a matching PDF using the filename convention.
   - If no match exists in the context folder, add it to the gap-list.
5. Write a fetch manifest at:
   `experiments/campaigns/pikp1_pwt7/inputs/pikp1_pwt7_fetch_manifest.md`

   It should contain:

   - **A header** recording the *recommended* source PDB and numbering convention, including the alignment-template choice for the PWT7 binding face.
   - **A "Papers I need but cannot access" section.** Each entry: citation, best URL, one- or two-sentence justification, suggested context-folder filename.
   - **A "Papers I read from context/" section** if any.
   - **An "Effector identification" section** specifically for this campaign, since PWT7 needs disambiguating: state what PWT7 is, the source organism, the MAX-fold family it belongs to (if known), the canonical UniProt accession, and the deposited PDB(s) if any. If the identification is uncertain, say so.
   - **A "Binding-face determination" section** stating whether the literature confirms PWT7 binds Pikp-HMA at a non-canonical face, with the supporting citation. If the face cannot be determined from the literature, note this as a critical gap that needs the user's input.

6. **Stop.** Do not write the report. End your turn with a brief chat summary:
   - The effector identification (resolved or unresolved).
   - The binding-face determination.
   - Number of sources successfully fetched online.
   - Number of sources read from the context folder.
   - Number of sources in the gap-list.
   - The recommended source-PDB and numbering convention.

---

## Quality criteria

- The gap-list determines what the user must do between Part 1 and Part 2. Be honest about which papers genuinely contain campaign-relevant detail.
- Justifications must be specific.
- If you find a paper open-access, use that and do not gap-list the publisher version.
- If PWT7 cannot be identified or the binding face cannot be determined, surface this clearly.

## What not to do

- Do not write the literature-sweep report.
- Do not modify any input complex PDB, its `.notes` file, an existing contig, or any pipeline scripts.
- Do not run the pipeline, derive new contigs, or submit anything to HPC.
- Do not commit or push.
- Do not include in the gap-list papers that you successfully fetched online or that are present in the context folder.
