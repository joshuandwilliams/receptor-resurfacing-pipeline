# Literature sweep and structure analysis — Part 1: Discovery

You are conducting the discovery phase of a structured literature, PDB, and UniProt sweep to inform the design of a binder-resurfacing campaign. The eventual report (produced in Part 2) will be used to (a) decide which structures to use as input to a structural alignment, (b) curate which residues should be redesigned vs preserved, and (c) define an RFDiffusion contig string that encodes those decisions.

This part is **discovery only**. You will identify every source you need, attempt to fetch each one, and produce a short gap-list of papers you could not access. **Do not write the report yet** — that is Part 2.

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
- Identifying the effector (full name, organism, family membership, any preprint or unpublished data) is part of the discovery work in Part 1; the user has limited information.
- The Pikp-1 HMA binder is shared across this campaign and the other pikp1_* campaigns; the receptor-side anchor / flex-region analysis from those campaigns is reusable here.

---

## Sources to consult

Search across all of the following. Cite specific papers, PDB entries, and UniProt records throughout the manifest you produce.

1. **Plant protein structures CSV** at `experiments/inputs/plant_protein_structures.csv`. Read this first as a starting reference. PWT3 is unlikely to appear since no structure exists, but the CSV may list related effectors.

2. **Scientific literature**. Identify the primary structural and functional papers for the receptor, and for PWT3 specifically. Pay particular attention to:
   - Survey / review papers cataloguing MAX-fold effectors and rice-blast effector families (PWL family in particular if PWT3 is from that family — verify in Part 1).
   - Any preprints, accepted manuscripts, or unpublished data describing PWT3 sequence, function, or structural prediction.
   - Mutagenesis or pathogenicity studies that map PWT3 sequence variation to virulence.
   - Engineering studies that have already attempted to redirect Pik-HMA recognition toward MAX-fold effectors that lack a deposited structure.

3. **Protein Data Bank**. For every relevant structure, you will eventually need: PDB ID, the complex contents, domain boundaries, the resolved ATOM range, resolution, and the publication. PWT3 itself has no PDB; concentrate on related effector / Pik-HMA complexes that can serve as alignment templates.

4. **UniProt**. Canonical sequences for the receptor and the effector. Domain annotations, signal peptides, post-translational modifications, feature annotations relevant to binding.

---

## Campaign context to read first

Before scoping the sweep, read:

- `experiments/README.md` — campaign-lifecycle and naming conventions.
- The campaign's own README at `experiments/campaigns/pikp1_pwt3/README.md` if present.
- Any existing complex PDB at `experiments/campaigns/pikp1_pwt3/inputs/` (if not present, note that complex assembly is an outcome of this sweep, not an input).
- Any existing contig file in the campaign's `inputs/` directory (if present).
- Outputs of the parallel pikp1_* campaigns' literature sweeps in their `inputs/` directories — the Pikp-1 HMA binder-side analysis is shared.
- `experiments/inputs/context/` — the shared paper PDF folder. List its contents now. Filename convention is `<firstauthor>_<year>_<shortjournal>.pdf`. You will check this folder for any paper you can't fetch online.

---

## What to do

1. Read all the campaign context above.
2. Plan the literature/PDB/UniProt sweep. Enumerate every paper, PDB entry, and UniProt record you intend to consult. **Be thorough.**
3. For each source, attempt to fetch it. Use web search to find it, web fetch to retrieve it. Prefer open-access mirrors (PMC, eLife, PLoS, bioRxiv) over publisher sites.
4. For each source you cannot fetch online:
   - Check `experiments/inputs/context/` for a matching PDF using the filename convention.
   - If no match exists in the context folder, add it to the gap-list.
5. Write a fetch manifest at:
   `experiments/campaigns/pikp1_pwt3/inputs/pikp1_pwt3_fetch_manifest.md`

   It should contain:

   - **A header** recording the *recommended* source PDB and numbering convention (no complex exists yet).
   - **A "Papers I need but cannot access" section.** Each entry: citation, best URL, one- or two-sentence justification, suggested context-folder filename.
   - **A "Papers I read from context/" section** if any.
   - **An "Effector identification" section** specifically for this campaign, since PWT3 needs disambiguating: state what PWT3 is, the source organism, the MAX-fold family it belongs to (if known), the canonical UniProt accession (if any), and confirm whether any structure exists. If the identification is uncertain, say so.
   - **An "AF3 prediction plan" section** specifying: the input sequence to use for the PWT3 AF3 monomer prediction (UniProt accession + boundary or full preprocessor sequence), the validation steps (e.g. compare the AF3-predicted MAX fold to a related deposited MAX effector), and whether AF3-multimer should be run as a backup.

6. **Stop.** Do not write the report. End your turn with a brief chat summary:
   - The effector identification (resolved or unresolved).
   - Number of sources successfully fetched online.
   - Number of sources read from the context folder.
   - Number of sources in the gap-list.
   - The recommended source-PDB / numbering convention and AF3 prediction plan.

---

## Quality criteria

- The gap-list determines what the user must do between Part 1 and Part 2. Be honest about which papers genuinely contain campaign-relevant detail.
- Justifications must be specific.
- If you find a paper open-access, use that and do not gap-list the publisher version.
- If PWT3 cannot be identified from the open literature, surface this clearly — the user may need to provide a private reference or unpublished sequence.

## What not to do

- Do not write the literature-sweep report.
- Do not modify any input complex PDB, its `.notes` file, an existing contig, or any pipeline scripts.
- Do not run the pipeline, derive new contigs, or submit anything to HPC.
- Do not commit or push.
- Do not include in the gap-list papers that you successfully fetched online or that are present in the context folder.
