# Literature sweep — pikp1_pwt7 binder-resurfacing campaign

## 1. Summary

This campaign aims to redesign the binding face of the Pikp-1 integrated HMA domain
(UniProt **E9KPB5**, residues 186–263) to gain recognition of **PWT7**, an avirulence
effector from *Pyricularia oryzae* (Avena/Lolium blast pathotypes, isolate Br58; 98 aa
precursor, signal peptide 1–20, mature chain 21–98 = 78 residues; Asuke et al. 2023
*MPMI*). Like the pikp1\_pby2 campaign, this is a de novo gain-of-recognition design:
PWT7 is naturally recognised by wheat TKP **Rwt7** through its N-terminal HMA integrated
domain, and Pikp-1 HMA shows no natural binding to PWT7. The prior knowledge claim
— that PWT7 binds Pikp-HMA at a **different face from the canonical Pik-HMA/AVR-Pik
face** — is **confirmed**: the Rwt7-HMA/PWT7 crystal structure (PDB **9TFP**, 1.6 Å;
Yu et al. 2025 *bioRxiv*) reveals that PWT7 contacts the **N-terminal β1/M1 face** of
the HMA domain, involving M1 of Rwt7-HMA and "two faces" of Rwt7-HMA, with PWT7
residues N39 and D40 (precursor numbering) as the central side-chain contacts. This is
structurally distinct from all three binding modes seen in sibling campaigns (AVR-Pik β4
face, AVR-Pia β2/αA face, PBY2 β3 face). Chain A uses the shared rung-5 AF3 monomer
prediction of Pikp-1 HMA(186–263); chain B sources from **9TFP chain B** (PWT7; Banfield
lab HPUB data). The **primary redesign target is the N-terminal β1 region of Pikp-1 HMA
(pipeline positions ~1–30)**, which corresponds to the M1-involving "two faces" of
Rwt7-HMA that PWT7 contacts — a region that was anchored as Interface 1 scaffold in all
sibling campaigns. **No campaign-blockers** exist; the structural template (9TFP) and
PWT7 coordinates are Banfield lab data accessible to the user.

---

## 2. Available structures

Pipeline numbering (chain A 1–78, position 1 = Pikp-1 residue 186) is used internally
only. All table entries use literature numbering.

### 2.1 Native Pikp-1 HMA structures (no engineering mutations)

Full catalogue identical to sibling pikp1\_\* campaigns; see `pikp1_avrpikf_literature_sweep.md`.
Key summary:

| PDB | Complex contents | Resolved ATOM range | Resolution | Reference | Notes for this campaign |
|-----|------------------|--------------------|------------|-----------|-------------------------|
| 5A6P | Pikp-HMA apo | 186–258 | 2.10 Å | Maqbool 2015 *eLife* | C-term truncated; shared constraint. |
| 5A6W | Pikp-HMA / AVR-PikD | 186–258 | 1.60 Å | Maqbool 2015 *eLife* | Canonical AVR-Pik face; **NOT the PWT7 face**. |
| 6G10 | Pikp-HMA / AVR-PikD | 186–263 | 1.35 Å | De la Concepcion 2018 *Nat Plants* | **NOT the alignment template** for this campaign; canonical face only. |
| 6G11 | Pikp-HMA / AVR-PikE | 186–263 | 1.90 Å | De la Concepcion 2018 *Nat Plants* | Same comment as 6G10. |
| 6Q76 | Pikp-HMA / AVR-Pia | 186–263 | 1.90 Å | Varden 2019 *JBC* | AVR-Pia face; not applicable here. |

### 2.2 Engineered or variant Pikp-1 / Pik-allelic HMA structures

Full table in `pikp1_avrpikf_literature_sweep.md`. None directly relevant; all target the
canonical AVR-Pik face.

### 2.3 Effector structures — PWT7, Rwt7-HMA/PWT7 complexes, and related MAX templates

| PDB | Effector / complex | Method | Resolution | Reference | Campaign relevance |
|-----|--------------------|--------|------------|-----------|-------------------|
| **9TFP** | **Rwt7-HMA / PWT7 (WT complex)** | X-ray | **1.6 Å** | **Yu et al. 2025 *bioRxiv*** | **Primary complex template and effector source.** Confirms PWT7 MAX-fold. M1/N-terminal binding face confirmed. HPUB; Banfield lab depositor. **Chain assignment (verified from `experiments/inputs/structures/9TFP_refmac_final.pdb`):** chains A and C are PWT7 (residues 19–98 in PWT7 precursor numbering; N39 confirmed at residue 39 = ASN, D40 at residue 40 = ASP); chains B and D are Rwt7-HMA (residues 1–74 in 1-based domain numbering; M1 confirmed at residue 1 = MET). |
| **9TFT** | **Rwt7-HMA+ / PWT7 (engineered)** | X-ray | **0.99 Å** | **Yu et al. 2025 *bioRxiv*** | Highest-resolution PWT7 structure; alternative source for chain B. HPUB; Banfield lab. |
| 9TFQ | Rmo2-HMA+ / PBY2 / PWT7 tripartite | X-ray | 1.4 Å | Yu et al. 2025 *bioRxiv* | PWT7 in tripartite complex with Rmo2-HMA+ and PBY2; useful for understanding cross-specificity. HPUB. |
| 5A6W | AVR-PikD in Pikp-HMA complex | X-ray | 1.60 Å | Maqbool 2015 | MAX-effector fold reference for cross-comparison with PWT7 fold. |
| 2MYW | AVR-Pia free structure (NMR) | NMR | — | de Guillen 2015 *PLoS Pathog* | Canonical MAX fold reference for AF3 validation comparison. |

**No deposited structure of PWT7 alone.** The only PWT7 coordinates are in 9TFP and 9TFT (both HPUB).

### 2.4 Complexes used as alignment templates

| PDB | Receptor | Effector | Resolution | Reference | Status for this campaign |
|-----|----------|----------|------------|-----------|--------------------------|
| **9TFP** | **Rwt7-HMA (wheat TKP; chains B and D, ATOM residues 1–74)** | **PWT7 (chains A and C, ATOM residues 19–98)** | **1.6 Å** | **Yu et al. 2025 *bioRxiv*** | **Campaign complex template.** M1 of Rwt7-HMA contacts PWT7 at the N-terminal face — different from all canonical Pikp-HMA complexes. Cross-system alignment: superpose AF3 Pikp-1 HMA onto Rwt7-HMA from 9TFP **chain B (or D)**. HPUB. |
| 6G10 | Pikp-HMA | AVR-PikD | 1.35 Å | De la Concepcion 2018 | **NOT used.** Positions effector at canonical C-terminal AVR-Pik face, opposite from PWT7 binding face. |
| 6FU9 | Pikm-HMA | AVR-PikD | 1.20 Å | De la Concepcion 2018 | NOT used. Same rationale as 6G10. |

---

## 3. Structural-input decision

### 3.1 Fallback ladder analysis

**Rung 1: native receptor / target-effector complex crystal.** No Pikp-1 HMA / PWT7
complex exists. → **Not available.**

**Rung 2: native receptor / related-effector complex.** No HMA-domain protein that
naturally binds PWT7 in the same way as Pikp-1 HMA. Rwt7-HMA/PWT7 (9TFP) is a
**cross-system rung-2 analog** (the receptor partner is Rwt7-HMA from barley TKP, not
Pikp-1 HMA from rice NLR, but both are HMA domains with the same βαββαβ fold). Same
approach used in pikp1\_pby2 (which used 9TFO as the Rmo2-HMA/PBY2 cross-system
template). **Chosen for this campaign.**

**Rung 5 (AF3 monomer of native receptor — chain A).** As for all pikp1\_\* campaigns.
AF3 monomer of Pikp-1 HMA (E9KPB5, 186–263), pipeline 1–78.

**Rung 6 (AF3 monomer of effector — optional fallback).** If 9TFP chain B is not
accessible, run AF3 monomer of PWT7 (mature chain residues 21–98 of precursor; 78
residues). Validate fold against 2MYW (AVR-Pia NMR, RMSD < 3 Å over ≥40 residues
= MAX-like confirmed). If MAX-like, align onto 9TFP via the Rwt7-HMA core alignment
as above. If not MAX-like, fall to rung 7.

**Rung 7 (AF3-multimer).** Fallback if 9TFP and AF3-monomer both fail or give
non-MAX fold. Run AF3-multimer of Pikp-1 HMA (186–263) + PWT7 (mature 21–98).

### 3.2 Chosen structural input

| Chain | Source | Rung | Residue range | Numbering convention |
|-------|--------|------|---------------|----------------------|
| **A** (Pikp-1 HMA) | AF3 monomer of UniProt E9KPB5 (186–263) | 5 | 1–78 | Pipeline: 1 = Pikp-1 186; pipeline = literature − 185 |
| **B** (PWT7) | **9TFP chain A (or C)** — both are PWT7. The literature label "chain B" is the Rwt7-HMA receptor, not the effector. | 1 cross-system | 19–98 as deposited in 9TFP chain A | **PWT7 precursor numbering, verified directly: chain A residue 39 = ASN (N39, binding-critical; N39R abolishes binding), residue 40 = ASP (D40). Signal peptide 1–20, mature 21–98.** |

### 3.3 Alignment template and complex assembly recipe

The binding face confirmed by 9TFP is the **N-terminal β1 face** of Rwt7-HMA
(involving M1), which maps to the N-terminus of Pikp-1 HMA (pipeline position 1 =
Lys-186). **6G10 is not used.**

Assembly:
1. Run AF3 monomer of Pikp-1 HMA (E9KPB5, 186–263). Renumber 1–78.
2. Access 9TFP. Extract **Rwt7-HMA (chain B, residues 1–74; first residue Met-1)** and
   **PWT7 (chain A, residues 19–98)** coordinates. (Chains C and D are equivalent
   copies in the asymmetric-unit dimer if needed.)
3. Superpose the AF3 Pikp-1 HMA model onto **9TFP chain B (Rwt7-HMA)** using Cα
   alignment over the conserved βαββαβ HMA core. Pay particular attention to matching
   β1 (pipeline ~3–10 of Pikp-1 HMA) with the equivalent β1 of Rwt7-HMA, since the
   M1/β1 N-terminal region is the key contact face.
4. Adopt PWT7 (**9TFP chain A**) at its crystallographic position relative to the
   now-placed AF3 chain A. Write `pikp1_pwt7_complex.pdb`. PWT7 keeps its 9TFP
   precursor numbering (19–98); the chain label may be relabelled to "B" for pipeline
   convention.
5. **Verify in ChimeraX**: confirm that the N-terminal region of the AF3 Pikp-1 HMA
   model (pipeline positions 1–15) faces PWT7 in the assembled model, and that
   N39/D40 of PWT7 (binding-critical contacts) are proximal to pipeline positions
   ~1–15. Measure inter-atomic distances between PWT7 N39/D40 and the nearest
   Pikp-1 HMA residues in this region.

### 3.4 Cross-validation

Compare the assembled 9TFP-derived model to the AF3-multimer of Pikp-1 HMA + PWT7
(run as a validation step). If both independently support the N-terminal face as the
PWT7 contact region, the contig flex region is confirmed.

---

## 4. Binding interfaces

Numbering: receptor side uses both Pikp-1 full-length (literature) and pipeline (in
parentheses). Effector side uses PWT7 **precursor numbering** (signal peptide 1–20;
mature 21–98) throughout, matching Yu et al. 2025 convention.

### 4.1 The N-terminal β1/M1 face — the PWT7 binding surface

**Structural description (CONFIRMED from Yu et al. 2025, PDB 9TFP).** The
Rwt7-HMA/PWT7 interface is unlike all previously characterised HMA/MAX-effector
complexes in this project. Yu et al. 2025 state: "the Rwt7HMA/PWT7 interface displayed
many polar side-chain interactions across **two faces** of Rwt7HMA. Central to the
interaction interface were the sidechains of **N39 and D40 from PWT7**, while the
**first residue of Rwt7HMA (M1) formed extensive interactions with the effector**."
This binding mode involves the N-terminal face of the HMA domain, centred on M1 (the
first residue of the βαββαβ fold). The "two faces" indicates a more dispersed interface
than single-β-strand pairing seen in PBY2/Rmo2-HMA (β3/β2) or AVR-Pia/Pikp-HMA
(β2/αA). Full atomic details of both faces require the 9TFP supplementary data
(fig. S4B/C of Yu et al. 2025, HPUB) — confirm contact residues in ChimeraX.

**Comparison with sibling campaigns (CONFIRMED):**

| Campaign effector | HMA face used | HMA elements | Effector elements | Source |
|---|---|---|---|---|
| AVR-PikD/F (avrpikf) | C-terminal | β4, pipeline ~67–73 | β3 of MAX effector | De la Concepcion 2018 |
| AVR-Pia (avrpia) | αA/β2 (dimerization face) | α1+β2, pipeline ~19–42 | β2/β1/β6 of MAX | Guo 2018 *PNAS* |
| PBY2 (pby2) | β3 backbone face | β3, pipeline ~43–49 | β2 of PBY2 MAX | Yu et al. 2025 |
| **PWT7 (this campaign)** | **N-terminal β1/M1 face** | **β1 + "two faces", pipeline ~1–30** | **β2 of PWT7 MAX (N19/D20)** | **Yu et al. 2025** |

**The N-terminal face (pipeline ~1–30) was anchored in all sibling campaigns.** This
campaign uniquely targets it.

**Receptor-side residues (Pikp-1 HMA / pipeline position):**

Pikp-1 HMA residues equivalent to M1 of Rwt7-HMA and the "two faces" contacts cannot
be assigned individually without the 9TFP contact map. The following are **INFERRED**
by structural equivalence of M1 in the Rwt7-HMA fold to the corresponding N-terminal
residues in Pikp-1 HMA pipeline numbering:

| Pikp-1 residue (full) | Identity | Pipeline pos. | HMA element | Role (INFERRED) |
|---|---|---|---|---|
| **Lys-186** | Lys | **1** | N-terminal/pre-β1 | **Equivalent of M1 of Rwt7-HMA; predicted central contact** |
| Lys-188 | Lys | 3 | β1 start | β1 N-terminal face |
| Leu-189 | Leu | 4 | β1 | β1 face |
| Gln-190 | Gln | 5 | β1 | β1 face |
| Lys-191 | Lys | 6 | β1 | β1 face (Pikm has Lys here also) |
| Val-193 | Val | 8 | β1 | β1 hydrophobic face |
| Val-195 | Val | 10 | β1-α1 transition | β1/α1 junction |
| α1 residues (~196–215) | various | ~11–30 | α1 | Second face candidate (INFERRED) |

*Caveat: All assignments above are INFERRED by fold alignment of Rwt7-HMA onto Pikp-1
HMA. Atomic contact residues require ChimeraX analysis of 9TFP once accessible.*

**This is a fundamentally different target region from all sibling campaigns.** In the
AVR-Pik system, the β1 region (pipeline ~3–10) is part of Interface 1 — a scaffold
region that mediates conserved main-chain contacts and was explicitly left fixed in
pikp1\_avrpikf. For PWT7, β1 must be redesigned.

**Effector-side residues (PWT7 precursor numbering):**

| PWT7 residue (precursor) | Identity | Mature position | MAX element | Role |
|---|---|---|---|---|
| **N39** | Asn | 19 (mature) | β2 region | **Central side-chain contact to Rwt7-HMA** |
| **D40** | Asp | 20 (mature) | β2 region | **Central side-chain contact to Rwt7-HMA** |
| (other contacts TBD from 9TFP) | — | — | — | Require 9TFP contact analysis |

**Binding-critical mutation confirmed (CONFIRMED from Yu et al. 2025):** PWT7 **N39R**
abolishes binding to Rwt7-HMA in vitro (aSEC: no peak shift; ITC: no binding) and in
planta (barley spray assay: full susceptibility on Rwt7 wheat; wheat protoplast
luciferase assay: no immune activation). KD Rwt7-HMA/PWT7 = **2.92 nM** by ITC.

**Specificity role.** PWT7 variants (PWT7Ele with N19/D20 equivalent; PWT7Lom,
PWT7Set with variable positions elsewhere) all retain avirulence on wheat with Rwt3
(Asuke 2023 *MPMI*, Fig. 3A). The variable positions (highlighted yellow in Asuke 2023
Fig. 3B) are outside the N39/D40 region, confirming that N39/D40 are the
specificity-determining contacts. The N39/D40 binding-critical contacts are conserved
across all functional PWT7 alleles.

**Design implication.** The β1/N-terminal face of Pikp-1 HMA (pipeline ~1–30) is the
**primary redesign target.** Pikp-1 HMA currently does not bind PWT7, presumably
because its β1/N-terminal residues do not complement PWT7's N39/D40 contacts. The goal
is to redesign this face to provide productive contacts with PWT7. The rest of the HMA
(β2-α2-β4-C-tail, pipeline ~31–78) is anchored.

### 4.2 Non-contact regions — anchor

- **α1 (pipeline ~11–30):** The α1 helix is the "second face" candidate. If 9TFP
  analysis shows α1 contacts, part of it enters the flex block (see §6). Under the
  narrow primary contig assumption, α1 is anchored.
- **β2/β3/α2/β4/C-tail (pipeline ~31–78):** Canonical AVR-Pik Interfaces 2 and 3 lie
  here. No PWT7 contacts expected based on the N-terminal binding face.

---

## 5. Numbering translation

### 5.1 Conventions

| Convention | What it means | Offset |
|---|---|---|
| **Pikp-1 full-length (literature)** | Residues 1–1142; HMA 186–263 | — |
| **Pipeline numbering (chain A)** | 1–78; pipeline = full-length − 185 | pipeline = literature − 185 |
| **PWT7 precursor numbering** | 1–98; signal peptide 1–20; mature 21–98 | — |
| **PWT7 mature numbering** | 1–78; mature residue = precursor − 20 | N39 precursor = N19 mature |
| **9TFP chain assignments (verified)** | Chains A and C = PWT7 (residues 19–98 in PWT7 precursor numbering; N39 confirmed at residue 39, D40 at residue 40). Chains B and D = Rwt7-HMA (residues 1–74 in 1-based domain numbering; M1 at residue 1). | PWT7 chains use precursor numbering directly |

### 5.2 Receptor key-residue translation table

| Pikp-1 residue (full) | Identity | Pipeline pos. | HMA element | PWT7-campaign role |
|---|---|---|---|---|
| **186** | Lys | **1** | Pre-β1 / N-terminus | **M1-equivalent; PRIMARY target** |
| 187 | — | 2 | N-terminal loop | β1 approach |
| 188 | Lys | 3 | β1 start | β1 N-terminal face |
| 189 | Leu | 4 | β1 | β1 face |
| 190 | Gln | 5 | β1 | β1 face |
| 191 | Lys | 6 | β1 | β1 face |
| 193 | Val | 8 | β1 | β1 face |
| 195 | Val | 10 | β1-α1 | β1/α1 junction |
| ~196–215 | various | ~11–30 | α1 | Second face candidate (confirm from 9TFP) |
| 218 | Ser | 33 | β2 | **NOT redesign target** (Interface 2 in AVR-Pik, anchored here) |
| 262 | Lys | 77 | C-tail | **NOT redesign target** (Interface 3 in AVR-Pik, anchored here) |

### 5.3 Effector key-residue table (PWT7)

| PWT7 precursor | Identity | Mature equiv. | MAX fold element | Role |
|---|---|---|---|---|
| **N39** | Asn | **N19** (mature) | β2 | **Central contact; N39R abolishes binding** |
| **D40** | Asp | **D20** (mature) | β2 | Central contact |
| R23/C24 | Arg/Cys | R3/C4 (mature) | β1? | Near mature N-terminus; structural |
| C44 | Cys | C24 (mature) | β2/β3 junction | Disulfide bond candidate |
| C76 | Cys | C56 (mature) | β5/β6 region | Disulfide bond candidate |

---

## 6. Contig string design

All contigs are written in **pipeline numbering** (chain A 1–78). Chain B is PWT7 (fixed).

The flex region covers the **N-terminal β1/M1 face** of Pikp-1 HMA — the face that
contacts PWT7 in the Rwt7-HMA/PWT7 complex. In Pikp-1 HMA pipeline numbering:
- β1: ~positions 3–10 (Interface 1 in AVR-Pik; here the PRIMARY redesign target)
- Pre-β1/M1 equivalent: positions 1–2
- α1 N-terminal portion: positions 11–30 (potential second face)
- β2 through C-tail: positions 31–78 (anchor)

**⚠ Boundary caveat (both contigs):** The flex/anchor boundary at positions ~12–13
(β1-α1 junction) and ~30–31 (α1-β2 junction) must be verified in ChimeraX on the
assembled complex before committing. If either boundary falls mid-helix, shift to the
nearest loop residue.

**⚠ Contact caveat:** The precise contact residues on Pikp-1 HMA require qtPISA
analysis of the assembled 9TFP-derived complex in ChimeraX. The contigs below are
derived from the confirmed binding face (N-terminal β1/M1 region) but should be
revised if the 9TFP contact map identifies residues outside pipeline 1–30.

### 6.1 Design principle

PWT7 binds Rwt7-HMA via M1 (N-terminal) + two faces. On Pikp-1 HMA, M1 corresponds
to pipeline position 1 (Lys-186). The β1 strand (pipeline ~3–10) + the beginning of
α1 (~11–20) constitute the "two faces" hypothesis. The rest of the HMA (β2 through
C-tail, pipeline ~31–78) is not expected to contact PWT7.

- **N-terminal face (pipeline 1–30):** β1 + α1 (N-terminal). *Flex — PRIMARY redesign
  target.* Includes the M1-equivalent and the most likely second face.
- **C-terminal scaffold (pipeline 31–78):** β2, β3, α2, β4, C-tail. *Anchor.* This is
  the canonical AVR-Pik Interface 2/3 region — not the PWT7 contact face.

### 6.2 Primary contig

```
15-28/A31-78 B
```

**Pikp-1 full-length (literature) cross-reference**: `15-28/A216-263 B` — i.e. flex
15–28 (over native residues 186–215, which is Pikp-1 Lys-186 through residue 215),
followed by anchor 216–263. The flex-block length range (15–28) is unchanged between
numberings. Note that the flex block is the **N-terminal** segment of chain A (no
preceding anchor), so the flex spans the chain N-terminus.

| Segment | Pipeline positions | Pikp-1 residues | Native length | Sampled length | Structural element | Role |
|---|---|---|---|---|---|---|
| `15-28` | 1–30 (flex) | 186–215 | 30 | 15–28 | Pre-β1 + β1 + α1 | **Flexible.** Covers M1-equivalent (pipeline 1), full β1 (~3–10), and α1 N-terminal portion (~11–30). Native length 30 residues, sampling 15–28 (±6 from native) gives RFDiffusion freedom to reshape the N-terminal face. |
| `A31-78` | 31–78 | 216–263 | 48 | anchor | β2 + β3 + α2 + β4 + C-tail | **Anchored.** Preserves the C-terminal scaffold and canonical AVR-Pik contact regions native. |
| `B` | full chain B | PWT7 | — | — | MAX effector | **Fully fixed.** |

**Boundary checks:**
- Position 30/31 (flex → anchor): α1 ends around Pikp-1 residue ~215 (pipeline ~30)
  and β2 begins at ~216–218 (pipeline ~31–33). This should be a loop/turn region.
  Verify in ChimeraX that position 31 is in the α1-β2 loop, not mid-helix.
- Beginning of flex at pipeline 1 is the chain N-terminus — no preceding anchor needed.

### 6.3 Alternative contig — conservative, β1-only flex

```
6-12/A13-78 B
```

**Literature cross-reference**: `6-12/A198-263 B`.

| Segment | Pipeline positions | Sampled length | Role |
|---|---|---|---|
| `6-12` | 1–12 (flex) | 6–12 | **Flexible.** Covers pre-β1 + β1 only (~12 residues, native). Conservative target; focuses on the M1-equivalent face. |
| `A13-78` | 13–78 | anchor | α1 through C-tail anchored. |

Use if the primary contig's flex region is too broad and produces structurally incoherent designs. **Boundary caveat:** position 12/13 may fall at the β1-α1 junction or into early α1 — verify in ChimeraX.

### 6.4 Alternative contig — aggressive, two-block design (IF second face confirmed from 9TFP)

```
6-12/A13-25/6-12/A31-78 B
```

**Literature cross-reference**: `6-12/A198-210/6-12/A216-263 B`.

| Segment | Pipeline positions | Sampled length | Role |
|---|---|---|---|
| `6-12` | 1–12 (flex) | 6–12 | **Flexible.** β1 face. |
| `A13-25` | 13–25 | anchor | β1-α1 loop + early α1 (non-contact if second face is not here). |
| `6-12` | 26–30 (flex) | 6–12 | **Flexible.** Second face, if at α1 C-terminal region (~pipeline 26–30). |
| `A31-78` | 31–78 | anchor | β2 through C-tail. |

**Note:** This contig is valid only if ChimeraX analysis of 9TFP shows contacts at a
second, non-adjacent HMA region (pipeline ~26–30). Do not use unless 9TFP contact
map confirms a second contact cluster in this region. The two-block design places a
zero-residue anchor between the flex blocks, which may not be structurally sensible;
adjust boundaries based on actual contact residue positions.

### 6.5 Alternative contig — generous flex with wide range

```
20-38/A41-78 B
```

**Literature cross-reference**: `20-38/A226-263 B`.

Flex: pipeline 1–40 (pre-β1 + β1 + full α1 + β2 beginning), native ~40, sampling
20–38. Use if the primary contig underperforms; gives RFDiffusion maximum freedom on
the N-terminal half of the HMA.

**Boundary caveat:** position 40/41 may fall in the β2 strand — verify in ChimeraX.

### 6.6 Relationship to sibling campaign contigs

This campaign is structurally unique: all other pikp1\_\* campaigns anchor pipeline
positions 1–30 or most of that range, while this campaign **flexes them**. The
asymmetry relative to siblings is a direct consequence of the different PWT7 binding
face and is expected.

---

## References

- Asuke S, Horie A, Komatsu K, Mori R, Vy TTP, Inoue Y, Jiang Y, Tatematsu Y,
  Shimizu M, Tosa Y. *Loss of PWT7, Located on a Supernumerary Chromosome, Is
  Associated with Parasitic Specialization of Pyricularia oryzae on Wheat.* Molecular
  Plant–Microbe Interactions 36, 716–725 (2023). doi:10.1094/MPMI-06-23-0078-R.
  PMID: 37432132. **[from context/ as asuke_2023_mpmi.pdf]**
- Asuke S, Tagle AG, Hyon G-S, Koizumi S, Murakami T, Horie A, Niwamoto D, Katayama E,
  Shibata M, Takahashi Y, Islam MT, Matsuoka Y, Yamaji N, Shimizu M, Terauchi R,
  Hisano H, Sato K, Tosa Y. *Evolution of HMA-integrated tandem kinases accompanied
  by expansion of target pathogens.* *bioRxiv* (2025). doi:10.64898/2025.12.15.692859.
  **[from context/ as asuke_2025_bioarxiv.pdf]**
- De la Concepcion JC, Franceschetti M, Maqbool A, Saitoh H, Terauchi R, Kamoun S,
  Banfield MJ. *Polymorphic residues in rice NLRs expand binding and response to
  effectors of the blast pathogen.* Nature Plants 4, 576–585 (2018).
  doi:10.1038/s41477-018-0194-x. **[from context/ as delaconcepcion_2019_natplants.pdf]**
- de Guillen K, Ortiz-Vallejo D, Gracy J, Fournier E, Kroj T, Padilla A. *Structure
  analysis uncovers a highly diverse but structurally conserved effector family in
  phytopathogenic fungi.* PLoS Pathogens 11, e1005228 (2015).
  doi:10.1371/journal.ppat.1005228.
- Guo L, Cesari S, de Guillen K, Chalvon V, Mammri L, Ma M, Meusnier I, Bonnot F,
  Padilla A, Peng Y-L, Liu J, Kroj T. *Specific recognition of two MAX effectors by
  integrated HMA domains in plant immune receptors involves distinct binding surfaces.*
  Proceedings of the National Academy of Sciences 115, 11637–11642 (2018).
  doi:10.1073/pnas.1810705115. **[from context/ as guo_2018_pnas.pdf]**
- Maqbool A, Saitoh H, Franceschetti M, Stevenson CEM, Uemura A, Kanzaki H, Kamoun S,
  Terauchi R, Banfield MJ. *Structural basis of pathogen recognition by an integrated
  HMA domain in a plant NLR immune receptor.* eLife 4, e08709 (2015).
  doi:10.7554/eLife.08709.
- Varden FA, Saitoh H, Yoshino K, Franceschetti M, Kamoun S, Terauchi R, Banfield MJ.
  *Cross-reactivity of a rice NLR immune receptor to distinct effectors from the rice
  blast pathogen Magnaporthe oryzae provides partial disease resistance.* Journal of
  Biological Chemistry 294, 13006–13016 (2019). doi:10.1074/jbc.RA119.007730.
- **Yu DS, Zdrzałek R, Katayama E, Akiyama H, Daykin L, Williams NJ, Goodridge I,
  Asuke S, Banfield MJ.** *Engineering plant tandem kinase immune receptors expands
  effector recognition profiles.* *bioRxiv* (2025). doi:10.64898/2025.12.15.694194.
  **[from context/ as yu_2025_bioarxiv.pdf]**
- UniProt entry: **E9KPB5** (Pikp-1, *Oryza sativa* subsp. japonica, 1142 aa).
- PWT7 sequence: 98 aa precursor (signal peptide 1–20; mature 21–98 = 78 residues).
  No public UniProt accession. DDBJ project accessions **PRJDB14584 and PRJDB14749**
  (Asuke et al. 2023).
- PDB entries: **9TFP** (Rwt7-HMA/PWT7, 1.6 Å, HPUB; Banfield lab depositor),
  **9TFT** (Rwt7-HMA+/PWT7, 0.99 Å, HPUB), **9TFQ** (Rmo2-HMA+/PBY2/PWT7,
  1.4 Å, HPUB).
