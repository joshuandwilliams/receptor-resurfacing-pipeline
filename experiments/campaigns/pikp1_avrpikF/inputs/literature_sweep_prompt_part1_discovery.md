# Literature sweep and structure analysis — Part 1: Discovery

You are conducting the discovery phase of a structured literature, PDB, and UniProt sweep to inform the design of a binder-resurfacing campaign. The eventual report (produced in Part 2) will be used to (a) decide which structures to use as input to a structural alignment, (b) curate which residues should be redesigned vs preserved, and (c) define an RFDiffusion contig string that encodes those decisions.

This part is **discovery only**. You will identify every source you need, attempt to fetch each one, and produce a short gap-list of papers you could not access. **Do not write the report yet** — that is Part 2.

---

## Inputs

**Receptor**: Pikp-1 HMA
**Receptor UniProt**: E9KPB5
**Receptor domain boundaries**: 186-263

**Effector**: AVR-PikF

**Prior knowledge to take into account**:
- All deposited Pikp-1 HMA crystal structures are missing C-terminal residues that are critical for binding. The native crystal 6G10 truncates before this region. The engineered structures (SNK-EKE: 7QPX, 7QZD) appear from their construct boundaries to include the full C-terminus, but ATOM records show several C-terminal residues are unresolved (disordered) in the deposited coordinates. None of the deposited structures provides a complete C-terminal interface.
- Therefore the structural-input decision for THIS campaign overrides the default fallback ladder: use rung 5 (AlphaFold3 monomer prediction of native Pikp-1 HMA, residues 186-263) rather than rung 4 (engineered crystal). The user has already verified the AF3 prediction is reliable for this domain — there are multiple high-resolution crystal structures of closely-related sequences to anchor it, and the C-terminal residues that are disordered in the crystals can be modelled by AF3 with reasonable confidence.
- The AF3 prediction is renumbered from 1 (position 1 = Pikp-1 residue 186; position 78 = residue 263). The contig must be written in this 1-78 numbering.
- The contig should redesign both Interface 2 (β2/β3 region) AND Interface 3 (C-terminal region). Do NOT treat the engineered SNK-EKE residues at Interface 3 as a fixed anchor — the campaign explicitly aims to find alternatives to that engineered solution rather than build on top of it. The C-terminal region must be inside a flexible block.

---

## Sources to consult

Search across all of the following. Cite specific papers, PDB entries, and UniProt records throughout the manifest you produce.

1. **Plant protein structures CSV** at `experiments/inputs/plant_protein_structures.csv`. Read this first as a starting reference for what structures already exist for this receptor and effector. **Treat it as a starting point, not exhaustive** — verify each relevant entry by retrieving the actual PDB record, and search PDB independently for any structures the CSV may have missed.

2. **Scientific literature**. Identify the primary structural and functional papers for both the receptor and the effector. Pay particular attention to:
   - Structural studies (crystal, cryo-EM, NMR) of the receptor, the effector, or any complex between them.
   - Mutagenesis studies that map specific residues to binding affinity, specificity, or function.
   - Allelic-series or variant studies if the receptor or effector is part of a family with multiple known forms.
   - Engineering studies that have already attempted to redirect binding specificity.

3. **Protein Data Bank**. For every relevant structure, you will eventually need: PDB ID, the complex contents, domain boundaries, the resolved ATOM range, resolution, and the publication.

4. **UniProt**. Canonical sequences for the receptor and effector. Domain annotations, signal peptides, post-translational modifications, feature annotations relevant to binding.

---

## Campaign context to read first

Before scoping the sweep, read:

- `experiments/README.md` — campaign-lifecycle and naming conventions.
- The campaign's own README at `experiments/campaigns/pikp1_avrpikf/README.md` if present.
- The assembled input complex PDB at `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_complex.pdb` and its `.notes` file. Determine (a) the source PDB it was assembled from, and (b) the residue-numbering convention used. Both feed the eventual report's structural-input decision and numbering-translation sections.
- Any existing contig file in the campaign's `inputs/` directory — it tells you which numbering convention the campaign is currently working in.
- `experiments/inputs/context/` — the shared paper PDF folder. List its contents now. Filename convention is `<firstauthor>_<year>_<shortjournal>.pdf` (e.g. `maqbool_2015_elife.pdf`). You will check this folder for any paper you can't fetch online.

---

## What to do

1. Read all the campaign context above.
2. Plan the literature/PDB/UniProt sweep. Enumerate every paper, PDB entry, and UniProt record you intend to consult to write the eventual report. **Be thorough.** Include sources you have alternative-source coverage for; the gap-list is the only signal of what the user needs to fetch manually.
3. For each source, attempt to fetch it. Use web search to find it, web fetch to retrieve it. Prefer open-access mirrors (PMC, eLife, PLoS, bioRxiv) over publisher sites (ScienceDirect, Wiley, Nature subscription content). Do not retry a publisher URL after a 403; if an open-access mirror exists, use it.
4. For each source you cannot fetch online:
   - Check `experiments/inputs/context/` for a matching PDF using the filename convention. Match on first-author surname, year, and short journal token. If you find one, read it.
   - If no match exists in the context folder, add it to the gap-list.
5. Write a fetch manifest at:
   `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_fetch_manifest.md`

   It should contain:

   - **A header** recording the source-PDB and numbering-convention determinations from step 1.
   - **A "Papers I need but cannot access" section.** Each entry:
     - Citation (author, year, journal).
     - Best URL you found (canonical, e.g. DOI or PMC link).
     - One- or two-sentence justification: what the report needs from this paper. Be specific — name the residue, polymorphism, or mechanism, not "background".
     - Suggested context-folder filename, following the convention `<firstauthor>_<year>_<shortjournal>.pdf`.
   - **A "Papers I read from context/" section** if any, so the user can see at a glance that the local cache was used. Just citation + filename; no justification needed.

6. **Stop.** Do not write the report. End your turn with a brief chat summary:
   - Number of sources successfully fetched online.
   - Number of sources read from the context folder.
   - Number of sources in the gap-list.
   - The source-PDB and numbering-convention determinations.

---

## Quality criteria

- The gap-list determines what the user must do between Part 1 and Part 2. A bloated gap-list wastes their time; a sparse one produces a half-blind report. Be honest about which papers genuinely contain campaign-relevant detail you cannot get from open-access sources.
- Justifications must be specific enough that the user can decide per-paper whether it is worth fetching. "May contain useful information" is not specific enough.
- If you find a paper open-access on PMC, eLife, PLoS, bioRxiv, or similar, use that and do not gap-list the publisher version.

## What not to do

- Do not write the literature-sweep report. Part 2 writes the report.
- Do not modify the input complex PDB, its `.notes` file, the existing contig, or any pipeline scripts.
- Do not run the pipeline, derive new contigs, or submit anything to HPC.
- Do not commit or push.
- Do not include in the gap-list papers that you successfully fetched online or that are present in the context folder.
