# Literature sweep — pikp1_avrpia binder-resurfacing campaign

## 1. Summary

This campaign aims to redesign the binding face of the Pikp-1 integrated HMA domain
(UniProt **E9KPB5**, residues 186–263) so that it recognises the rice-blast effector **AVR-Pia**
(UniProt **A0A6P8B9Y5**, mature chain residues 20–85) with greater affinity and specificity
than wild-type Pikp-1 currently achieves. Pikp-1 already cross-reacts weakly with AVR-Pia in
vitro and confers partial disease resistance in rice (Varden et al. 2019 *JBC*); the campaign
goal is to discover redesigned receptor variants with stronger, more robust recognition. Unlike
the canonical Pik-HMA / AVR-Pik interactions, AVR-Pia binds Pikp-1 HMA at the **α1/β2
face** — the same surface RGA5-HMA uses for AVR1-CO39 recognition — which is on the
opposite side of the HMA fold from the AVR-Pik binding face (Guo et al. 2018 *PNAS*; Varden
2019). The deposited native Pikp-HMA / AVR-Pia complex (PDB **6Q76**, 1.90 Å; Varden 2019)
is directly usable as the complex template for chain-B placement; chain A uses the shared
AF3 monomer prediction of native Pikp-1 HMA(186–263) (rung 5), consistent with all
pikp1\_\* campaigns. The **primary redesign targets** are the second half of α1 and the start of
β2 (pipeline positions ~14–42), which carry the five structurally characterised Pikp-HMA
contacts to AVR-Pia; the remainder of the HMA domain is anchored. There are **no
campaign-blockers**: a usable complex template, a complete effector structure, and a
well-characterised mutagenesis literature are all present.

---

## 2. Available structures

Pipeline numbering (chain A 1–78, position 1 = Pikp-1 residue 186, position 78 = Pikp-1
residue 263) is *only* used internally by this campaign; residue numbers in the tables below
use the **literature numbering** as deposited in each PDB entry. Where chain A deposited
numbering coincides with pipeline numbering this is noted explicitly.

### 2.1 Native Pikp-1 HMA structures (no engineering mutations)

| PDB | Complex contents | Construct (Pikp-1 full-length residues) | Resolved ATOM range | Resolution | Reference | Notes for this campaign |
|-----|------------------|-----------------------------------------|---------------------|------------|-----------|-------------------------|
| 5A6P | Pikp-HMA apo (dimer) | 186–258 | 186–258 | 2.10 Å | Maqbool 2015 *eLife* | Truncates before C-terminal tail; no Interface 3 and no AVR-Pia binding-face context. |
| 5A6W | Pikp-HMA / AVR-PikD | 186–258 | 186–258 | 1.60 Å | Maqbool 2015 *eLife* | Same truncation as 5A6P. AVR-PikD binds the opposite face from AVR-Pia. |
| 6G10 | Pikp-HMA / AVR-PikD | 186–263 | C-terminal residues disordered in ATOM records | 1.35 Å | De la Concepcion 2018 *Nat Plants* | **Not the alignment template for this campaign** — AVR-PikD uses the opposite HMA face. The α1/β2 face (the AVR-Pia face) is fully resolved here and can serve as a cross-validation reference for AF3 chain A quality in the redesign region. |
| 6G11 | Pikp-HMA / AVR-PikE | 186–263 | as 6G10 | 1.90 Å | De la Concepcion 2018 *Nat Plants* | Same as 6G10; not used as campaign template. |
| **6Q76** | **Pikp-HMA / AVR-Pia** | **186–263** | **Chain A: ATOM residues 186–260 in Pikp-1 full-length numbering (75 residues; gap at 199–200; C-terminal 261–263 not resolved). Chain B: ATOM residues 18–85 (68 residues, no gaps).** | **1.90 Å** | **Varden 2019 *JBC*** | **Campaign complex template.** The only deposited structure of Pikp-1 HMA bound to this campaign's target effector. Binding face confirmed as α1/β2. **6Q76 chain A is numbered 186–260 (Pikp-1 full-length); pipeline numbering is therefore literature − 185 (e.g., pipeline 1 = 6Q76 chain A residue 186).** The missing positions 76–78 (residues 261–263) are NOT part of the AVR-Pia interface and do not affect chain-B placement or the redesign region. |

**Disqualifying constraint on chain A (shared across all pikp1\_\* campaigns).** All five entries
above either truncate or disorder C-terminal residues. Although positions 76–78 are not part
of the AVR-Pia interface, the campaign adopts the shared pikp1\_\* policy of using an AF3
monomer prediction for chain A so that the full 186–263 construct is available in the pipeline
input PDB.

### 2.2 Engineered or variant Pikp-1 / Pik-allelic HMA structures

| PDB | Complex contents | Mutations vs Pikp-1 WT | Construct | Resolution | Reference | Notes for this campaign |
|-----|------------------|------------------------|-----------|------------|-----------|-------------------------|
| 6R8K | Pikp-HMA-NK-KE / AVR-PikD | Asn261Lys + Lys262Glu | 186–263 | 1.60 Å | De la Concepcion 2019 *eLife* | Mutations in Interface 3 (pipeline 76, 77); not relevant to the AVR-Pia face. |
| 6R8M | Pikp-HMA-NK-KE / AVR-PikE | same | 186–263 | 1.85 Å | De la Concepcion 2019 *eLife* | As above. |
| 7A8W | Pikp-HMA-NK-KE / AVR-PikC | same | 186–263 | 2.15 Å | Maidment 2021 *PLoS Pathog* | As above. |
| 7A8X | Pikh-HMA / AVR-PikC | Pikh natural allele (Asn261Lys) | 186–263 | 2.30 Å | Maidment 2021 *PLoS Pathog* | Pikh/Pikp differ at Interface 3 only; not relevant to AVR-Pia binding face. |
| 7QPX | Pikp-HMA-SNK-EKE / AVR-PikC | Ser258Glu + Asn261Lys + Lys262Glu | 186–264 | 2.05 Å | Maidment 2023 *eLife* | Engineering in Interface 3; not relevant to AVR-Pia face. |
| 7QZD | Pikp-HMA-SNK-EKE / AVR-PikF | same | 186–264 | 2.20 Å | Maidment 2023 *eLife* | As above. |
| 6FU9 | Pikm-HMA / AVR-PikD | Pikm natural allele | 186–264 | 1.20 Å | De la Concepcion 2018 *Nat Plants* | Highest-resolution HMA fold reference for AF3 validation. **Important note:** Pikm-HMA does NOT bind AVR-Pia — the residue equivalent to Pikp-1 Asp-217 (pipeline 32) is histidine in Pikm (Maidment 2025 *bioRxiv*), explaining loss of the conserved salt bridge to AVR-Pia Arg-43. |
| 6FUB | Pikm-HMA / AVR-PikE | Pikm | 186–264 | 1.30 Å | De la Concepcion 2018 *Nat Plants* | Reference for HMA fold. |
| 6FUD | Pikm-HMA / AVR-PikA | Pikm | 186–264 | 1.30 Å | De la Concepcion 2018 *Nat Plants* | Reference for HMA fold. |
| 7BNT | Ancestral Pik-1 HMA / AVR-PikD | Reconstructed ancestral sequence | n/a | 1.32 Å | Białas 2021 *eLife* | Background on HMA evolution. Not used directly. |
| 8B2R | Engineered RGA5-HMA mutant / AVR-PikF | Six engineered RGA5-HMA mutations | n/a | 1.22 Å | Bentham 2023 *Plant Cell* | Cross-system engineering precedent; not a direct template. |

### 2.3 Effector structures — AVR-Pia alone or in complex

| PDB | Effector | Complex partner | Construct (effector precursor residues) | Method | Resolution | Reference | Notes for this campaign |
|-----|----------|------------------|-----------------------------------------|--------|------------|-----------|-------------------------|
| 2MYW | AVR-Pia (WT) | None (free) | ~20–85 (mature) | NMR | — | de Guillen 2015 *PLoS Pathog* | AVR-Pia free-solution structure; MAX effector β-sandwich fold established. |
| 2N37 | AVR-Pia (WT) | None (free) | 20–85 (mature, 66 residues) | NMR | — | Ose 2015 *J Biomol NMR* | Alternate NMR deposit; assigns backbone resonances used in Ortiz 2017 and Maidment 2025 NMR titrations. |
| 5JHJ | AVR-Pia mutant H3 (F24S/T46N) | None (free) | ~20–85 | NMR | — | Ortiz 2017 *Plant Cell* | Structure of the naturally occurring inactive variant; backbone similar to WT — mutations alter surface properties not fold. Used to validate that binding surface is altered by F24S. |
| **6Q76** | **AVR-Pia (WT)** | **Pikp-1 HMA** | **18–85 (68 residues, all modelled, no gaps)** | **X-ray** | **1.90 Å** | **Varden 2019 *JBC*** | **Source of chain B in the assembled complex.** The binding face (β2 of AVR-Pia antiparallel to HMA β2) is directly visible. **Chain B uses AVR-Pia precursor numbering directly (verified: residue 24 = Phe, residue 41 = Tyr, residue 43 = Arg, residue 46 = Thr, residue 85 = Tyr — all consistent with the precursor sequence at those positions).** |
| **9RSV** | **AVR-Pia (WT)** | **OsHPP09-HMA** | **~20–85 (70 residues)** | **X-ray** | **1.65 Å** | **Maidment 2025 *bioRxiv*** | Higher-resolution analog complex on the same binding face as 6Q76. RMSD 0.357 Å vs 6Q76 AVR-Pia chain; interface area 521.7 Å² vs 460.7 Å² in 6Q76 due to more side-chain H-bonds. Can serve as secondary cross-validation template for chain-B conformation. Chain B numbering uses precursor numbering (residues R43, T46, Y41, R23 labeled in Fig 6). |

**Natural variation: AVR-Pia-H3 (F24S/T46N).** This naturally occurring allele abolishes binding to
Pikp-HMA, OsRGA5-HMA, OsHPP09-HMA, OsHPP10-HMA, OsHPP11-HMA, and OsHIPP21-HMA
(Cesari 2013, Ortiz 2017, Maidment 2025). F24S alone is sufficient to abolish RGA5-HMA binding
(Cesari 2013). These mutations define which face of AVR-Pia is the functionally relevant
binding surface.

### 2.4 RGA5-HMA structures (structural comparators for the αA/β2 binding face)

| PDB | Complex contents | Resolution | Reference | Notes for this campaign |
|-----|------------------|------------|-----------|-------------------------|
| 5ZNE | RGA5-HMA apo | 1.78 Å | Guo 2018 *PNAS* | RGA5-HMA homodimer via αA/β2 surface — the same face that binds AVR-Pia. Illustrates the dimerization geometry competed by effector binding. |
| 5ZNG | RGA5-HMA / AVR1-CO39 | 2.19 Å | Guo 2018 *PNAS* | **Structural comparator.** AVR1-CO39 binds the αA/β2 face of RGA5-HMA — the same orientation as AVR-Pia binding to Pikp-HMA (Guo 2018, Maidment 2025). Key RGA5-HMA contacts: V1025, D1026 (area 1); V1028, E1029 (area 2); R1012, I1030 (area 3). KD = 5.4 µM (ITC; Guo 2018). KD of RGA5-HMA / AVR-Pia = 1.8 µM (Guo 2018). |
| 7DV8 | RGA5-HMA2 | 2.45 Å | Liu 2021 *PNAS* | RGA5 variant; background only. |
| 9RSV | AVR-Pia / OsHPP09-HMA | 1.65 Å | Maidment 2025 *bioRxiv* | See §2.3 above. OsHPP09-HMA is an sHMA, not a NLR-integrated domain, but uses the identical binding face. |

---

## 3. Structural-input decision

### 3.1 Fallback ladder analysis

**Rung 1: native receptor / target-effector complex crystal.** PDB **6Q76** (Pikp-1 HMA /
AVR-Pia; Varden 2019; 1.90 Å) is available. Chain A (Pikp-HMA) is resolved for residues
**186–260 in Pikp-1 full-length numbering (= pipeline 1–75; gap at residues 199–200)**
and chain B (AVR-Pia) for **residues 18–85 in precursor numbering (68 residues, no gaps)**.
For the AVR-Pia campaign the AVR-Pia binding face lies entirely within the resolved chain A
range (the binding contacts are at pipeline positions 19, 27, 32, 34, 41, all ≤75), so the
missing C-terminal residues (pipeline 76–78) do not compromise the complex geometry.

However, the shared pikp1\_\* campaign policy mandates rung 5 (AF3 monomer prediction)
for **chain A** because all deposited Pikp-1 HMA structures truncate or disorder the C-terminal
tail. Even though positions 76–78 are non-contacts for AVR-Pia, the decision is retained for
consistency: it avoids creating a different class of input PDB for this campaign and ensures
the full 186–263 sequence is present in the pipeline input. Chain B is not subject to this
constraint and sources directly from 6Q76.

**Rung 2–4.** Not required; rung-1 effector coordinates are available.

**Rung 5 (AF3 monomer of native receptor — chain A).** An AF3 monomer prediction of
native Pikp-1 HMA (UniProt **E9KPB5**, residues 186–263) is the **chosen source for chain A**,
reused from the sibling pikp1\_\* campaigns. Renumbered 1–78, position 1 = Pikp-1 residue
186. Validation procedure: superpose onto 6Q76 chain A (positions 1–73 overlap, RMSD on
Cα should be <1.0 Å) and onto 6FU9 chain A (Pikm-HMA / AVR-PikD; 1.20 Å; the
highest-resolution Pik-HMA reference).

### 3.2 Chosen structural input

| Chain | Source | Rung | Residue range in complex | Numbering convention |
|-------|--------|------|--------------------------|----------------------|
| **A** (Pikp-1 HMA) | AF3 monomer prediction of UniProt E9KPB5 (186–263) | 5 | 1–78 | Pipeline: position 1 = Pikp-1 186; pipeline = literature − 185. (The AF3 PDB at `experiments/inputs/structures/AF3/Pikp-1_HMA.pdb` is already deposited with chain A renumbered 1–78.) |
| **B** (AVR-Pia) | 6Q76 chain B | 1 | 18–85 as deposited (68 residues, all modelled, no gaps) | Precursor numbering (verified — see §5.3) |

### 3.3 Complex assembly recipe

1. Run AF3 monomer prediction for Pikp-1 HMA (186–263 of E9KPB5). Renumber output 1–78.
2. Open 6Q76 in ChimeraX. Extract chain B (AVR-Pia) coordinates.
3. Superpose the AF3 chain A model (residues 1–78) onto 6Q76 chain A (residues 186–260
   in Pikp-1 full-length numbering) by Cα alignment over the structurally equivalent residues
   (AF3 pipeline 1–75 ↔ 6Q76 chain A 186–260, with the offset pipeline = literature − 185).
   Note that the two chains do **NOT** share residue numbers — superposition must be done
   on Cα geometry, not on residue-number identity. This places the AF3 model at the correct
   binding pose for chain B placement.
4. Combine the AF3 chain A (positions 1–78) with 6Q76 chain B coordinates to produce
   `pikp1_avrpia_complex.pdb`. Chain A contributes positions 1–78 (including the AF3-modelled
   C-terminal tail 76–78 absent from 6Q76). Chain B is taken directly from 6Q76 without modification.
5. Verify the binding pose: inspect that Asp-217 (pipeline 32) is in salt-bridge geometry with
   AVR-Pia Arg-43 (the conserved contact present in 6Q76, 9RSV, and predicted by structural
   alignment with 5ZNG; Varden 2019, Maidment 2025, Guo 2018).

**Note on 9RSV as cross-validation.** 9RSV (Maidment 2025) at 1.65 Å provides higher-resolution
coordinates of AVR-Pia and independently confirms the binding geometry. If the assembled complex
from step 4 shows any ambiguity in chain-B conformation, aligning 9RSV chain B (AVR-Pia) onto the
pose is a valid cross-check. 9RSV uses the same β2-antiparallel binding mode and AVR-Pia is
essentially identical (RMSD 0.357 Å vs 6Q76 chain B).

### 3.4 Why 6G10 is not used

The user's prior knowledge explicitly flags that **6G10 is not the right template** for AVR-Pia,
and this is confirmed by the structural literature. In 6G10 (Pikp-HMA / AVR-PikD), AVR-PikD
docks via β4 of the HMA against β3 of the effector — the opposite face of the HMA from
where AVR-Pia binds. Aligning using 6G10 would place an AVR-Pia model at the wrong
interface. 6Q76 is the only correct template for this campaign.

---

## 4. Binding interfaces

The Pikp-1 HMA / AVR-Pia binding interface is a **single surface** on the α1/β2 face of the
HMA fold. This is the same face used by RGA5-HMA to bind AVR1-CO39 (Guo 2018) and by
OsHPP09-HMA to bind AVR-Pia (Maidment 2025). It is the **opposite face** from the canonical
AVR-Pik binding surface (Interface 1/2/3 in De la Concepcion 2018). Receptor-side residue
numbers use Pikp-1 full-length numbering; pipeline positions are given in parentheses.
Effector residue numbers use AVR-Pia **precursor numbering** (signal peptide 1–19; mature
20–85) throughout.

### 4.1 The α1/β2 face — the AVR-Pia binding surface

**Structural description.** The binding interface involves the second half of α1 (Pikp-1
approximately 199–215, pipeline positions approximately 14–30) and the first part of β2
(Pikp-1 approximately 216–227, pipeline positions approximately 31–42). The binding
geometry is dominated by **antiparallel β-strand alignment**: β2 of AVR-Pia (MAX effector
β-sheet 2) runs antiparallel to β2 of Pikp-HMA, forming a continuous antiparallel β-sheet
incorporating β6, β1, and β2 of the MAX fold together with the four-stranded β-sheet of
the HMA domain (Guo 2018; Maidment 2025). The α1 helix of Pikp-HMA provides additional
side-chain contacts from above. Buried interface area: **460.7 Å²** in 6Q76 (Varden 2019;
Maidment 2025), compared to 521.7 Å² in 9RSV (OsHPP09-HMA / AVR-Pia) and 492.8 Å²
in 5ZNG (RGA5-HMA / AVR1-CO39).

This surface also mediates **Pikp-1 HMA homodimerization** in solution and in crystal
(Guo 2018: "in the Pikp-1HMA/AVR-PikD complex, the β2/αA surface is free for
self-interaction and, consequently, Pikp-1HMA occurs as a dimer in the Pikp-1HMA/AVR-PikD
crystal structure"). AVR-Pia and AVR1-CO39 each compete with RGA5-HMA self-interaction
by binding the same surface (Guo 2018). Redesigning this face may therefore alter Pikp-1
HMA dimerization behavior, though the significance of HMA homodimerization in the context
of the full-length NLR receptor is unknown (Guo 2018).

**Receptor-side residues (Pikp-1 full-length / pipeline position):**

| Pikp-1 residue | Identity | Pipeline pos. | Secondary element | Contact type | AVR-Pia residue | Notes |
|---|---|---|---|---|---|---|
| **Ser-204** | Ser | **19** | α1 (mid-helix) | Water-mediated H-bond | Tyr-41 (side chain hydroxyl) | Varden 2019. OsHPP09-HMA has Gln20 at the equivalent structural position, which makes a **direct** H-bond to Tyr-41 (Maidment 2025). Ser-204 makes a weaker, water-mediated equivalent contact, contributing to the lower affinity of Pikp-HMA vs OsHPP09-HMA. |
| **Ser-212** | Ser | **27** | α1 (C-terminal end) | H-bond | Tyr-85 (side chain hydroxyl) | Varden 2019. OsHPP09-HMA has Asp28 at the equivalent position, which instead forms a salt bridge to AVR-Pia Arg-23 (Maidment 2025) — a stronger, different interaction. |
| **Asp-217** | Asp | **32** | α1/β2 junction | Salt bridge / H-bond | Arg-43 (side chain) | Varden 2019; Maidment 2025. **The most conserved contact.** Asp-217 is structurally equivalent to OsHPP09-HMA Asp33 (salt bridge to Arg-43) and OsRGA5-HMA Asp1026 (water-mediated contact to AVR1-CO39 Thr-41). **Pikm-HMA has His at this position and does NOT bind AVR-Pia** (Maidment 2025). Loss of this Asp is likely to abolish binding. |
| **Val-219** | Val | **34** | β2 (N-terminal end) | Main-chain / backbone contact | Tyr-41 (backbone) | Varden 2019. Part of the β-strand antiparallel pairing. |
| **Arg-226** | Arg | **41** | β2 (mid-strand) or β2-β3 loop | H-bond | Leu-38 (backbone) | Varden 2019; Maidment 2025. OsHPP09-HMA Lys19 is in the equivalent structural position to OsPikp-1 Arg203 (pipeline 18), but it is Arg-226 (pipeline 41) in Pikp-1 — from the β2-β3 loop region — that makes the analogous H-bond to AVR-Pia Leu-38 (Maidment 2025). |

**Effector-side residues (AVR-Pia precursor numbering):**

| AVR-Pia residue | Identity | MAX fold element | Contact type | HMA contact | Notes |
|---|---|---|---|---|---|
| Arg-23 | Arg | N-terminal loop / pre-β1 | Salt bridge | OsHPP09 Asp28 (Maidment 2025); not directly assigned in Varden 2019 for Pikp | Mutation R43G abolishes RGA5-HMA binding (Ortiz 2017). |
| **Phe-24** | Phe | β1 | Hydrophobic (+ surface patch) | Contact with HMA β2/α1 face | **F24S abolishes binding to RGA5-HMA and all H(I)PPs tested** (Cesari 2013; Maidment 2025). Present in the binding-surface hydrophobic patch (de Guillen 2015). |
| **Leu-38** | Leu | β2 | Backbone H-bond donor | Arg-226 / pipeline 41 (Varden 2019) | Backbone pairing in the antiparallel β2/β2 alignment. |
| **Tyr-41** | Tyr | β2 (C-terminal end) | Side-chain H-bond + backbone | Ser-204 (water-mediated, Varden 2019); OsHPP09 Gln20 (direct H-bond, Maidment 2025); Val-219 backbone | Central to the antiparallel β-sheet interface. |
| **Arg-43** | Arg | β2-β3 loop | Salt bridge | **Asp-217 (pipeline 32)** | Key conserved contact. In 9RSV, Arg-43 forms salt bridge to OsHPP09 Asp33 (= Pikp-1 Asp-217). |
| **Thr-46** | Thr | β2-β3 loop | — | — | T46N in H3 variant. No direct HMA contact described, but T46N contributes to loss of binding alongside F24S (Cesari 2013; Maidment 2025). May affect local β2 conformation. |
| Tyr-85 | Tyr | β6 | H-bond | Ser-212 (Varden 2019) | Part of the β6/β1/β2 MAX effector surface that contacts the HMA β-sheet. |

**Specificity role.** AVR-Pia is natively recognised by RGA4/RGA5 (Cesari 2013), not by
Pik-1/Pik-2. Pikp-HMA's cross-reactivity with AVR-Pia at ~7.8 µM (comparable to
RGA5-HMA; Guo 2018) is weak relative to AVR-Pik effectors. The campaign exploits the fact
that OsHPP09-HMA achieves KD = 115 nM for the same effector through additional side-chain
contacts at the same face. The key bottleneck is that Pikp-HMA Ser-204 (pipeline 19) makes
only a water-mediated contact where OsHPP09 Gln-20 makes a direct H-bond, and the Asp28
(OsHPP09) / Ser-212 (Pikp-1) difference at pipeline 27 eliminates a salt bridge with
AVR-Pia Arg-23.

**Design implication. Primary redesign target: pipeline positions 14–42 (second half of α1 and
all of β2).** This is the region that contains all five structurally characterised contacts and that
could, in principle, adopt higher-affinity geometries analogous to OsHPP09-HMA. Asp-217
(pipeline 32) is the most conserved contact across all HMA domains that bind AVR-Pia
(Pikp-1, OsHPP09, OsRGA5); its identity should be treated with care — RFDiffusion is free
to change it in the flex block but designs that discard the Asp here are likely to have reduced
AVR-Pia binding.

### 4.2 Non-contact regions — anchor

The following regions have **no documented contact** with AVR-Pia and should be anchored:

- **β1 (pipeline 1–13, Pikp-1 186–198):** The short N-terminal β-strand. No AVR-Pia contacts.
  Structurally equivalent to the RGA5-HMA β1 region which also lacks AVR-Pia / AVR1-CO39
  contacts in both 5ZNG and 6Q76.
- **β2-tail / β3 / α2 / β4 / C-terminal tail (pipeline 43–78, Pikp-1 228–263):** These
  elements form the AVR-Pik binding surface (Interfaces 2 residue Glu-230, Interface 3) and
  the structural scaffold (α2 helix, β4 main-chain). None contact AVR-Pia. Anchoring this block
  preserves the overall HMA fold and the non-AVR-Pia-facing surfaces.

**Constitutive-activity caveat.** α1 (the primary redesign region) is a surface helix and
its redesign is unlikely to destabilise the HMA fold or trigger autoactivity. The risk is lower
than for Interface 3 redesign (which is C-terminal and adjacent to known autoactivity-relevant
positions). Nevertheless, RFDiffusion designs with many changes to positions that pack
against β3 or α2 internally should be flagged for inspection.

---

## 5. Numbering translation

### 5.1 Conventions

| Convention | What it means | Offset |
|---|---|---|
| **Pikp-1 full-length (literature)** | Pikp-1 residues 1–1142; HMA domain 186–263 | — |
| **Pipeline numbering (chain A)** | 1–78, position 1 = Pikp-1 residue 186 | pipeline = literature − 185 |
| **6Q76 chain A** | **Residues 186–260 as deposited (75 residues; gap at 199–200; uses Pikp-1 full-length numbering)** | pipeline = 6Q76 chain A − 185 (i.e. they do NOT coincide; superposition by Cα geometry) |
| **AVR-Pia precursor (literature)** | Residues 1–85; signal peptide 1–19; mature chain 20–85 | — |
| **6Q76 chain B** | **Residues 18–85 as deposited (68 residues, all modelled, no gaps). Uses precursor numbering directly (verified — see §5.3)** | chain B = precursor |
| **9RSV chain B** | Uses precursor numbering directly (R23, R43, T46, Y41 labeled in Maidment 2025 Fig 6) | chain B = precursor |

### 5.2 Receptor key-residue translation table

| Pikp-1 residue (full) | Identity | Pipeline position | Secondary element | Role |
|---|---|---|---|---|
| 186 | — | 1 | β1 (N-term) | Scaffold |
| 198 | — | 13 | β1/α1 transition | Scaffold; proposed anchor boundary (verify) |
| 204 | **Ser** | **19** | α1 | AVR-Pia contact: water-mediated to Tyr-41 |
| 212 | **Ser** | **27** | α1 (C-terminal region) | AVR-Pia contact: H-bond to Tyr-85 |
| 217 | **Asp** | **32** | α1/β2 junction | **Conserved AVR-Pia contact: salt bridge to Arg-43** |
| 219 | **Val** | **34** | β2 (N-end) | AVR-Pia contact: backbone/Tyr-41 |
| 226 | **Arg** | **41** | β2-β3 loop | AVR-Pia contact: H-bond to Leu-38 backbone |
| 228 | — | 43 | β2-β3 loop or β3 start | Proposed anchor boundary (verify) |
| 245–251 | — | 60–66 | α2 | Scaffold helix; not AVR-Pia contact |
| 252–258 | — | 67–73 | β4 | AVR-Pik Interface 3 region; not AVR-Pia contact |
| 263 | Asp | 78 | C-terminal | End of construct |

### 5.3 Effector chain B numbering — confirmed

**Confirmed by direct inspection of `experiments/inputs/structures/6Q76.pdb`:** chain B is
deposited with residues 18–85 in **AVR-Pia precursor numbering** (68 residues, all modelled,
no gaps). The construct begins two residues upstream of the predicted signal-peptide cleavage
(precursor 20) and runs to the C-terminus (precursor 85). Verified residue identities:

| Chain B residue | Identity | Role |
|---|---|---|
| 24 | Phe | F24S abolishes binding (Cesari 2013) |
| 38 | Leu | β2 backbone H-bond donor |
| 41 | Tyr | β2 H-bond / backbone partner |
| 43 | Arg | Salt bridge to Pikp-1 Asp-217 (pipeline 32) |
| 46 | Thr | T46N in H3 variant |
| 85 | Tyr | β6 H-bond to Pikp-1 Ser-212 |

This matches the precursor labels used by Maidment 2025 Fig 6 exactly. No renumbering of
chain B is required for complex assembly — 6Q76 chain B coordinates can be lifted directly
into `pikp1_avrpia_complex.pdb`.

### 5.4 Effector key-residue table

| AVR-Pia precursor residue | Identity | MAX fold element | Binding role | Notes |
|---|---|---|---|---|
| 23 | Arg | Loop/pre-β1 | Salt bridge partner (OsHPP09 D28; not confirmed for Pikp) | |
| **24** | **Phe** | **β1** | **Hydrophobic surface** | F24S abolishes binding (Cesari 2013, Maidment 2025) |
| 38 | Leu | β2 | Backbone H-bond acceptor (Arg-226/pipeline 41) | |
| 41 | Tyr | β2 | H-bond to HMA α1 (Ser-204/pipeline 19; Glu-204 of OsHPP09) | Central to β2/β2 pairing |
| **43** | **Arg** | β2-β3 loop | **Salt bridge to Asp-217 / pipeline 32** | Most conserved effector-side contact |
| **46** | **Thr** | β2-β3 loop | Indirect | T46N in H3; loss of binding (Cesari 2013) |
| 85 | Tyr | β6 | H-bond to Ser-212 (Varden 2019) | On β6 of MAX effector |

---

## 6. Contig string design

All contigs are written in **pipeline numbering** (chain A 1–78). Chain B in the contig syntax
corresponds to the effector chain in the pipeline (6Q76 chain B in the assembled complex).

### 6.1 Design principle from §4

The literature analysis pins down two distinct regions on Pikp-1 HMA:

- **N-terminal scaffold (positions 1–13)** — β1 and β1-α1 loop. No AVR-Pia contacts. *Anchor.*
- **AVR-Pia binding face (positions ~14–42)** — second half of α1 and β2. Contains all five
  characterised contacts: Ser-204(19), Ser-212(27), Asp-217(32), Val-219(34), Arg-226(41).
  *Flex* — this is the entire redesign target for this campaign.
- **C-terminal scaffold (positions 43–78)** — β3, α2, β4, and C-terminal tail. No AVR-Pia
  contacts. Includes the AVR-Pik Interface 2/3 residues that are irrelevant here. *Anchor.*

This produces a **single contiguous flex block** covering the entire AVR-Pia binding face. Unlike
the pikp1\_avrpikf campaign (which required two separate flex blocks for Interface 2 and
Interface 3), the AVR-Pia interface is one unified surface and benefits from a single block that
allows RFDiffusion to optimise the complementarity of α1 and β2 together.

### 6.2 Primary contig

```
A1-13/22-36/A43-78 B
```

**Pikp-1 full-length (literature) cross-reference**: `A186-198/22-36/A228-263 B` —
i.e. anchor residues 186–198 + flex 22–36 (covering native residues 199–227) + anchor
residues 228–263. The flex-block length range (22–36) is unchanged between numberings.
Use the pipeline-numbered form for RFDiffusion; use the literature form when checking
positions against figures in 6Q76, 6G10, or any Banfield-lab paper.

| Segment | Pipeline positions | Pikp-1 residues | Native length | Sampled length | Structural element | Role |
|---|---|---|---|---|---|---|
| `A1-13` | 1–13 | 186–198 | 13 | 13 (anchor) | β1 + β1-α1 loop | **Anchored.** Preserves β1 scaffold native. |
| `22-36` | ~14–42 (flex) | ~199–227 | ~29 | 22–36 | α1 (second half) + β2 | **Flexible.** Encompasses all AVR-Pia contacts: Ser-204(19), Ser-212(27), Asp-217(32), Val-219(34), Arg-226(41). |
| `A43-78` | 43–78 | 228–263 | 36 | 36 (anchor) | β3 + β3-α2-β4 transition + α2 + β4 + C-terminal tail | **Anchored.** Preserves the C-terminal half of the HMA domain native. |
| `B` | full chain B | AVR-Pia | — | — | MAX effector | **Fully fixed.** RFDiffusion does not modify the target. |

**Why this contig follows from §4.** The flex block covers the entire α1/β2 face confirmed by
6Q76 and validated by the OsHPP09/AVR-Pia comparison in 9RSV and the AVR1-CO39/RGA5
comparison in 5ZNG. The anchor blocks preserve the structural scaffold on both sides.

**Boundary check.** Both boundaries warrant ChimeraX verification before committing:
- **Position 13 / 14 boundary (anchor → flex).** β1 of Pikp-HMA spans approximately
  Pikp-1 189–195 (pipeline positions ~4–10). The β1-α1 loop and the start of α1 follow.
  Confirm position 13 falls in the β1-α1 loop or at the very end of β1 — it should not be
  mid-helix. If α1 begins at position 11 or 12, the boundary may need to shift to position 10 or 11.
- **Position 42 / 43 boundary (flex → anchor).** β2 of Pikp-HMA spans approximately
  Pikp-1 216–226 (pipeline ~31–41). The β2-β3 loop follows. Confirm position 43 is in or
  at the start of the β2-β3 loop — it should not be mid-strand. If the β2-β3 loop starts at
  position 41 or 42, the boundary is appropriate; if β2 extends to position 44 or 45, shift
  the boundary accordingly.

### 6.3 Alternative contig — conservative, native-length flex

```
A1-13/29-29/A43-78 B
```

**Literature cross-reference**: `A186-198/29-29/A228-263 B`.

Same anchor blocks as the primary contig but the flex block is fixed at the native segment
length (29 residues). Removes RFDiffusion's freedom to insert or delete residues in the
redesigned region; only sequence changes are sampled. Useful as a comparator to the primary
contig to disentangle the contribution of length variation from sequence variation alone.

### 6.4 Alternative contig — aggressive, broader flex range and extended N-terminal freedom

```
A1-10/20-40/A43-78 B
```

**Literature cross-reference**: `A186-195/20-40/A228-263 B`.

| Segment | Pipeline positions | Sampled length | Role |
|---|---|---|---|
| `A1-10` | 1–10 | anchor | β1 only; drops the β1-α1 loop from the anchor |
| `20-40` | ~11–42 (flex) | 20–40 | Entire α1 + β2; native length ~32, sampling 20–40 gives ±10 residues |
| `A43-78` | 43–78 | anchor | As primary |

Extends the flex block to include the full α1 helix from the start (position 11), giving
RFDiffusion freedom to reshape the entire α1-β2 surface. Wider length range accommodates
more structural exploration. Use if the primary contig yields designs that are too constrained
around the N-terminal part of the binding surface.

### 6.5 Alternative contig — split, α1 and β2 flexed separately

```
A1-13/12-20/A32-41/A42-78 B
```

**Literature cross-reference**: `A186-198/12-20/A217-226/A227-263 B`.

This anchors Asp-217 (32) and its immediate neighbourhood in the β2 anchor block, while
flexing only α1 (positions 14–31) and separately the β2 region including Arg-226 (41).
However, because Asp-217 (32) is the most conserved contact across all AVR-Pia binding HMA
domains (and its loss in Pikm-HMA explains non-binding), anchoring it in the native identity
ensures designs retain this contact. Use only if other contigs produce too many designs that
lack the Asp-217–Arg-43 salt bridge.

**Note:** The approach of using two separate flex blocks with a zero-anchor between them (one for
α1, one for β2) is structurally less sensible for AVR-Pia than for AVR-PikF, because the α1/β2
binding face is a single continuous surface rather than two geometrically distinct sub-interfaces.

### 6.6 Alternative contig — tight flex (conservative search)

```
A1-13/25-33/A43-78 B
```

**Literature cross-reference**: `A186-198/25-33/A228-263 B`.

Same as primary but narrower flex range (22–36 → 25–33; ±4 from native 29). Useful as a
tighter search if broad sampling in the primary contig converges on too wide a design distribution.

---

## References

- Białas A, Langner T, Harant A, Contreras MP, Stevenson CEM, Banfield MJ, Kourelis J,
  Kamoun S, Wu C-H, Krasileva KV. *Two NLR immune receptors acquired high-affinity binding
  to a fungal effector through convergent evolution of their integrated domain*. eLife 10,
  e66961 (2021). doi:10.7554/eLife.66961.
- Bentham AR, De la Concepcion JC, Benjumea JV, Kourelis J, Jones S, Mendel M, Stubbs J,
  Stevenson CEM, Maidment JHR, Youles M, Erickson FL, Sornay C, Adamiak G, Brunkard JO,
  Banfield MJ, Kamoun S. *Allelic compatibility in plant immune receptors facilitates engineering
  of new effector recognition specificities*. The Plant Cell 35, 3809–3827 (2023).
  doi:10.1093/plcell/koad204.
- Cesari S, Thilliez G, Ribot C, Chalvon V, Michel C, Jauneau A, Rivas S, Alaux L, Kanzaki H,
  Okuyama Y, Morel J-B, Fournier E, Tharreau D, Terauchi R, Kroj T. *The rice resistance
  protein pair RGA4/RGA5 recognizes the Magnaporthe oryzae effectors AVR-Pia and
  AVR1-CO39 by direct binding*. The Plant Cell 25, 1463–1481 (2013).
  doi:10.1105/tpc.113.112714. **[PMC3663280]**
- Cesari S, Bernoux M, Moncuquet P, Kroj T, Dodds PN. *A novel conserved mechanism for
  plant NLR protein pairs: the 'integrated decoy' hypothesis*. Frontiers in Plant Science 5,
  606 (2014). doi:10.3389/fpls.2014.00606.
- de Guillen K, Ortiz-Vallejo D, Gracy J, Fournier E, Kroj T, Padilla A. *Structure analysis
  uncovers a highly diverse but structurally conserved effector family in phytopathogenic
  fungi*. PLoS Pathogens 11, e1005228 (2015). doi:10.1371/journal.ppat.1005228.
- De la Concepcion JC, Franceschetti M, Maqbool A, Saitoh H, Terauchi R, Kamoun S,
  Banfield MJ. *Polymorphic residues in rice NLRs expand binding and response to effectors of
  the blast pathogen*. Nature Plants 4, 576–585 (2018). doi:10.1038/s41477-018-0194-x.
  **[from experiments/inputs/context/ as delaconcepcion_2019_natplants.pdf]**
- De la Concepcion JC, Franceschetti M, MacLean D, Terauchi R, Kamoun S, Banfield MJ.
  *Protein engineering expands the effector recognition profile of a rice NLR immune receptor*.
  eLife 8, e47713 (2019). doi:10.7554/eLife.47713.
- Guo L, Cesari S, de Guillen K, Chalvon V, Mammri L, Ma M, Meusnier I, Bonnot F, Padilla A,
  Peng Y-L, Liu J, Kroj T. *Specific recognition of two MAX effectors by integrated HMA domains
  in plant immune receptors involves distinct binding surfaces*. Proceedings of the National
  Academy of Sciences 115, 11637–11642 (2018). doi:10.1073/pnas.1810705115.
  **[from experiments/inputs/context/ as guo_2018_pnas.pdf]**
- Maidment JHR, Franceschetti M, Maqbool A, Saitoh H, Jantasuriyarat C, Kamoun S,
  Terauchi R, Banfield MJ. *Multiple variants of the fungal effector AVR-Pik bind the HMA
  domain of the rice protein OsHIPP19, providing a foundation to engineer plant defense*.
  Journal of Biological Chemistry 296, 100371 (2021). doi:10.1016/j.jbc.2021.100371.
- Maidment JHR, Shimizu M, Bentham AR, Vera S, Franceschetti M, Longya A, Stevenson CEM,
  De la Concepcion JC, Białas A, Kamoun S, Terauchi R, Banfield MJ. *Effector
  target-guided engineering of an integrated domain expands the disease resistance profile of
  a rice NLR immune receptor*. eLife 12, e81123 (2023). doi:10.7554/eLife.81123.
- Maidment JHR, Vera S, Franceschetti M, Bentham AR, De la Concepcion JC, Banfield MJ.
  *The allelic rice immune receptor Pikh confers extended resistance to strains of the blast
  fungus through a single polymorphism in the effector binding interface*. PLoS Pathogens 17,
  e1009368 (2021). doi:10.1371/journal.ppat.1009368.
- Maidment JHR, Saile SC, Bocquet A, Thivolle C, Bourcet L, Planel L-F, Gelin M, Kroj T,
  Padilla A, de Guillen K, Cesari S. *The Magnaporthe oryzae MAX effector AVR-Pia binds a
  novel group of rice HMA domain-containing proteins*. bioRxiv (2025).
  doi:10.1101/2025.07.11.664054. **[from experiments/inputs/context/ as
  maidment_2025_bioarxiv.pdf]**
- Maqbool A, Saitoh H, Franceschetti M, Stevenson CEM, Uemura A, Kanzaki H, Kamoun S,
  Terauchi R, Banfield MJ. *Structural basis of pathogen recognition by an integrated HMA
  domain in a plant NLR immune receptor*. eLife 4, e08709 (2015). doi:10.7554/eLife.08709.
- Ortiz D, de Guillen K, Cesari S, Chalvon V, Gracy J, Padilla A, Kroj T. *Recognition of the
  Magnaporthe oryzae effector AVR-Pia by the decoy domain of the rice NLR immune receptor
  RGA5*. The Plant Cell 29, 156–168 (2017). doi:10.1105/tpc.16.00435. **[PMC5304345]**
- Ose T, Oikawa A, Nakamura Y, Maenaka K, Higuchi Y, Satoh Y, Fujiwara S, Demura M,
  Sone T, Kamiya M. *Solution structure of an avirulence protein, AVR-Pia, from Magnaporthe
  oryzae*. Journal of Biomolecular NMR 63, 229–235 (2015). doi:10.1007/s10858-015-9979-7.
- Varden FA, Saitoh H, Yoshino K, Franceschetti M, Kamoun S, Terauchi R, Banfield MJ.
  *Cross-reactivity of a rice NLR immune receptor to distinct effectors from the rice blast
  pathogen Magnaporthe oryzae provides partial disease resistance*. Journal of Biological
  Chemistry 294, 13006–13016 (2019). doi:10.1074/jbc.RA119.007730. **[PMC6721932]**
- Zdrzałek R, Stevenson CEM, Pham Z, Schimmel L, Kourelis J, Kamoun S, Banfield MJ.
  *Bioengineering a plant NLR immune receptor with a robust binding interface toward a
  conserved fungal pathogen effector*. PNAS 121, e2402872121 (2024).
  doi:10.1073/pnas.2402872121.
- UniProt entries: **E9KPB5** (Pikp-1, *Oryza sativa* subsp. japonica, 1142 aa; HMA annotated
  188–257); **A0A6P8B9Y5** (AVR-Pia, *Pyricularia grisea*, 85 aa; signal peptide 1–19).
