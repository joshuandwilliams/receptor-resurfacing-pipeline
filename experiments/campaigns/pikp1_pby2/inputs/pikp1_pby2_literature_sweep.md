# Literature sweep — pikp1_pby2 binder-resurfacing campaign

## 1. Summary

This campaign aims to redesign the binding face of the Pikp-1 integrated HMA domain
(UniProt **E9KPB5**, residues 186–263) to gain recognition of **PBY2**, a MAX-fold
avirulence effector from *Pyricularia oryzae* (barley blast pathotype, isolate Br48;
GenBank **LC905591–LC905597**). PBY2 is natively recognised by the barley tandem
kinase protein (TKP) **Rmo2** via its integrated N-terminal HMA domain; Pikp-1 HMA
currently shows no binding to PBY2. Unlike the other pikp1\_\* campaigns — which
redesign Interfaces 2 or 3 (AVR-Pik face) or the α1/β2 face (AVR-Pia face) — PBY2
binds HMA through **β3 of the HMA domain** (predominantly main-chain contacts), a
third spatially distinct interface (Yu et al. 2025 *bioRxiv*). The campaign therefore
targets a different region of Pikp-1 HMA from all sibling campaigns. The receptor
chain (chain A) uses the shared rung-5 AF3 monomer prediction of Pikp-1 HMA(186–263).
The effector chain (chain B) is lifted from **PDB 9TFO chain B** (Rmo2-HMA/PBY2
complex at 2.2 Å; Yu et al. 2025). Complex assembly requires a **cross-system alignment**:
superposing the AF3 Pikp-1 HMA model onto Rmo2-HMA from 9TFO to position PBY2.
**No standard Pikp-HMA/AVR-Pik template (e.g. 6G10) is appropriate.** The primary
redesign target is the **β3 strand of Pikp-1 HMA (pipeline positions ~43–49)**, with
the adjacent β2-β3 and β3-α2 loops. There are **no campaign-blockers**: PBY2
coordinates, Rmo2-HMA/PBY2 structure, and mutagenesis data are all present.

---

## 2. Available structures

Pipeline numbering (chain A 1–78, position 1 = Pikp-1 residue 186, position 78 = Pikp-1
residue 263) is *only* used internally; residue numbers in the tables below use literature
numbering as deposited in each PDB entry.

### 2.1 Native Pikp-1 HMA structures (no engineering mutations)

Same set as all pikp1\_\* campaigns. All deposits truncate or disorder C-terminal
residues; see sibling pikp1\_avrpikf sweep for full table. Key entries:

| PDB | Complex contents | Resolved ATOM range | Resolution | Reference | Notes for this campaign |
|-----|------------------|--------------------|------------|-----------|-------------------------|
| 5A6P | Pikp-HMA apo (dimer) | 186–258 | 2.10 Å | Maqbool 2015 *eLife* | C-terminus truncated; same constraint as all pikp1_* campaigns. |
| 5A6W | Pikp-HMA / AVR-PikD | 186–258 | 1.60 Å | Maqbool 2015 *eLife* | AVR-Pik binding face (β4/β3 antiparallel) — NOT the PBY2 face. |
| 6G10 | Pikp-HMA / AVR-PikD | 186–263 | 1.35 Å | De la Concepcion 2018 *Nat Plants* | **Not the alignment template for this campaign.** AVR-Pik face; PBY2 binds via a spatially distinct β3 face. |
| 6Q76 | Pikp-HMA / AVR-Pia | 186–263 | 1.90 Å | Varden 2019 *JBC* | AVR-Pia face (β2/αA); also not applicable here. |

### 2.2 Rmo2-HMA and Rwt7-HMA structures — TKP integrated HMA domains

| PDB | Complex contents | Receptor | Resolution | Reference | Notes for this campaign |
|-----|------------------|----------|------------|-----------|-------------------------|
| **9TFO** | **Rmo2-HMA / PBY2** | **Rmo2-HMA in chains A and B (ATOM residues 2–76, 1-based domain numbering); PBY2 in chains C and D (ATOM residues 24–84 and 22–84 respectively, AVR-Pia-style precursor numbering)** | **2.2 Å** | **Yu et al. 2025 *bioRxiv*** | **Campaign complex template.** The only available structure of PBY2 bound to an HMA domain. **HPUB (on hold until publication).** Binding interface: β3 of Rmo2-HMA (main-chain) + α1 side chains **D15, K18** (verified: chain A residue 15 = ASP, residue 18 = LYS) / β2 of PBY2 + side chains **E42, S48** (verified: chain D residue 42 = GLU, residue 48 = SER). |
| 9TFP | Rwt7-HMA / PWT7 | Rwt7 (wheat TKP), residues 1–77 | 1.6 Å | Yu et al. 2025 *bioRxiv* | Companion structure for the sibling pikp1_pwt7 campaign. Different (spatially distinct) binding interface from 9TFO. HPUB. |
| 9TFQ | Rmo2-HMA+ / PBY2 / PWT7 (tripartite) | Rmo2-HMA+ (engineered) | 1.4 Å | Yu et al. 2025 *bioRxiv* | Dual-specificity engineered receptor. HPUB. |
| 9TFS | Rwt7-HMA+ / PBY2 | Rwt7-HMA+ (engineered) | 1.8 Å | Yu et al. 2025 *bioRxiv* | Gain-of-binding to PBY2 by resurfaced Rwt7. Provides additional PBY2 binding geometry. HPUB. |
| 9TFT | Rwt7-HMA+ / PWT7 | Rwt7-HMA+ (engineered) | 0.99 Å | Yu et al. 2025 *bioRxiv* | High-resolution PWT7 complex; relevant to sibling pwt7 campaign. HPUB. |

### 2.3 Effector structures — PBY2 alone or in complex

| PDB | Effector | Complex partner | Method | Resolution | Reference | Notes |
|-----|----------|-----------------|--------|------------|-----------|-------|
| **9TFO** (**chain D** — see note) | **PBY2 (WT)** | **Rmo2-HMA** | X-ray | **2.2 Å** | **Yu et al. 2025 *bioRxiv*** | **Source of effector chain.** MAX-fold effector; signal peptide 1–20, mature chain 21–85. Binding-critical residue: **E42** (verified at chain D residue 42 = GLU; E42R abolishes binding to Rmo2-HMA and resistance in planta). **Note on chain assignment:** in 9TFO chains A and B are both Rmo2-HMA (asymmetric-unit dimer); PBY2 is deposited as chains C (residues 24–84) and D (residues 22–84). Chain D is more complete and is the recommended source. |

No standalone PBY2 NMR or crystal structure exists outside the 9TFO complex.

**Key effector natural variant.** PBY2 (Br48 isolate) and PBY1 (MZ5-1-6 isolate)
are identical at the nucleotide level (Asuke et al. 2025 *bioRxiv*) and are both recognised
by the same Rmo2 alleles, confirming the binding surface is conserved between them.

### 2.4 Engineered or variant Pikp-1 / Pik-allelic HMA structures

Not directly relevant to this campaign (all existing engineering targets the AVR-Pik
face, Interface 2 or 3 — not the β3 face used by PBY2). See pikp1\_avrpikf sweep for
the full table. The engineered Rmo2-HMA+/PBY2 structure (9TFQ and 9TFS) demonstrates
that the PBY2 binding interface can be successfully resurfaced; those experiments serve
as proof-of-concept for this campaign.

---

## 3. Structural-input decision

### 3.1 Fallback ladder analysis

**Rung 1: native receptor / target-effector complex crystal.** No Pikp-1 HMA / PBY2
complex exists. Pikp-1 HMA does not bind PBY2 (no interaction detected; PBY2 is
recognised by barley Rmo2-HMA, not by rice Pikp-1 HMA). → **Not available.**

**Rung 2: native receptor / related-effector complex crystal.** PBY2 has been
crystallised bound to the related HMA domain of Rmo2 (PDB **9TFO**, 2.2 Å; Yu et al.
2025). Rmo2-HMA (barley TKP, residues 1–77) and Pikp-1 HMA (rice NLR, pipeline
residues 1–78) share the canonical HMA fold (βαββαβ topology) and sufficient sequence
similarity that the two domains superpose with sub-ångström RMSD on their shared
secondary-structure elements (comparable to the Pikp-1 / RGA5 RMSD of 0.90 Å
reported by Guo et al. 2018 *PNAS* for two other structurally related HMA domains).

The situation here is a **cross-system rung 2 analog**: the complex template uses a
different receptor (Rmo2-HMA, not Pikp-1 HMA), but the target effector (PBY2) is the
same and the HMA fold is conserved. Complex assembly proceeds by superposing the
AF3 Pikp-1 HMA chain A onto Rmo2-HMA from 9TFO chain A, then adopting the PBY2
coordinates from 9TFO chain B at their crystallographically determined position. **This
landing is rung 2 cross-system (with the caveat described below).**

**Caveat.** Because PBY2 does not currently bind Pikp-1 HMA, the assembled complex
from the Rmo2 superposition is a **starting-geometry model**, not a validated binding
pose. The actual β3 surface of Pikp-1 HMA differs in sequence from Rmo2-HMA, so the
precise PBY2 docking pose may shift. The contig must therefore give RFDiffusion
sufficient flexibility around the β3 region to find alternative geometries, not lock in
the 9TFO-derived pose as authoritative.

**Chain assignment in 9TFO (verified from `experiments/inputs/structures/9TFO_refmac_final.pdb`):**
chains A and B are both Rmo2-HMA (residues 2–76, 1-based domain numbering; an
asymmetric-unit dimer); chains C and D are both PBY2 (chain C residues 24–84, chain D
residues 22–84; precursor numbering with E42 confirmed at residue 42). Earlier text in
this sweep that refers to "9TFO chain B" as PBY2 is **incorrect** — PBY2 is chain C or
chain D, and chain D is the more complete copy.

**Rungs 3–4.** Rung 3 (native receptor alone + effector alone aligned onto a complex
template) is effectively what is being done: AF3 Pikp-1 HMA (receptor alone) + PBY2
from 9TFO (effector from a complex with a different receptor), assembled by structural
superposition. The distinction between rung 2 and rung 3 here is semantic; the approach
is the same.

**Rung 5 (AF3 monomer of native receptor — chain A).** Used, as for all pikp1\_\*
campaigns. An AF3 monomer prediction of native Pikp-1 HMA (UniProt **E9KPB5**,
residues 186–263) is the chosen source for chain A.

**Rung 7 (AF3-multimer of the complex).** Fall-back if 9TFO chain B (PBY2 coordinates)
cannot be accessed. AF3-multimer of Pikp-1 HMA + PBY2 would then provide the only
available binding-pose estimate. PBY2 sequence can be retrieved from GenBank
**LC905591–LC905597** (Asuke et al. 2025).

### 3.2 Chosen structural input

| Chain | Source | Rung | Residue range | Numbering convention |
|-------|--------|------|---------------|----------------------|
| **A** (Pikp-1 HMA) | AF3 monomer of UniProt E9KPB5 (186–263) | 5 | 1–78 | Pipeline: 1 = Pikp-1 186; pipeline = literature − 185 |
| **B** (PBY2) | **9TFO chain D** (PBY2 in the Rmo2-HMA / PBY2 complex; chain D preferred over chain C as the more complete copy, residues 22–84 vs 24–84) | 1 (cross-system) | 22–84 as deposited in 9TFO chain D | **PBY2 precursor numbering, verified directly: chain D residue 42 = GLU (E42, binding-critical), residue 48 = SER (S48)** |

### 3.3 Complex assembly recipe

1. Run AF3 monomer prediction for Pikp-1 HMA (E9KPB5 residues 186–263). Renumber 1–78.
2. Access PDB 9TFO. Extract Rmo2-HMA (**chain A**, residues 2–76) and PBY2 (**chain D**,
   residues 22–84 — chain D is the more complete copy of the two PBY2 chains in the
   asymmetric unit) coordinates.
3. Superpose the AF3 Pikp-1 HMA model onto 9TFO chain A (Rmo2-HMA) using Cα
   alignment over the conserved βαββαβ core. The β3 strand (Rmo2-HMA ~residues
   29–38) and its flanking secondary-structure elements are the most important
   structural anchors for this alignment.
4. Extract PBY2 (9TFO chain D) at its crystallographic position relative to the now-placed
   AF3 chain A. Write `pikp1_pby2_complex.pdb` with AF3 chain A (positions 1–78) and
   PBY2 (chain D residues 22–84). PBY2 keeps its 9TFO precursor numbering; the chain
   label may be relabelled to "B" for pipeline convention but the residue numbers are
   preserved.
5. **Verify in ChimeraX**: confirm that β3 of the AF3 Pikp-1 HMA model faces PBY2 β2
   in the assembled complex, and that the D15/K18-equivalent residues of Pikp-1 HMA
   α1 (~pipeline 18–22) are positioned near PBY2 E42 (~precursor position 42). The
   inter-atomic distances will not match the 9TFO crystal contacts exactly (Pikp-1 HMA
   and Rmo2-HMA differ in sequence), but the overall topology should be preserved.

### 3.4 Why no standard AVR-Pik or AVR-Pia template is used

- **6G10 (canonical Pik-HMA / AVR-Pik template)** positions AVR-Pik so that effector
  β3 faces HMA β4 — the opposite end of the HMA β-sheet from where PBY2 docks.
  Using 6G10 would place PBY2 at the wrong face of Pikp-1 HMA.
- **6Q76 (AVR-Pia template)** positions AVR-Pia at the αA/β2 face — also wrong for
  PBY2.
- The only valid structural basis is **9TFO**, which is the Rmo2-HMA/PBY2 complex.

---

## 4. Binding interfaces

PBY2 engages a **single, structurally defined interface** on Rmo2-HMA that is
**spatially distinct** from both the canonical AVR-Pik interface (β4 of HMA / β3 of
effector) and the AVR-Pia / AVR1-CO39 interface (β2 of HMA / β2 of effector; αA
dimerization face) (Yu et al. 2025 *bioRxiv*; Guo et al. 2018 *PNAS*). Receptor-side
residue numbers use Rmo2-HMA's own 1-based domain numbering (1–77), with the
**inferred equivalent Pikp-1 HMA pipeline positions** given in parentheses (derived by
structural alignment of the shared HMA βαββαβ fold; exact equivalences require ChimeraX
confirmation). Effector residue numbers use PBY2 **precursor numbering** (signal
peptide 1–20; mature chain 21–85; GenBank LC905591).

### 4.1 The β3/α1 face — the PBY2 binding surface

**Structural description.** The PBY2 binding interface on Rmo2-HMA is centred on the
**β3 strand of the HMA fold** (~residues 29–38 in Rmo2-HMA domain numbering).
Binding involves **antiparallel β-strand pairing**: β2 of PBY2 (the MAX effector) aligns
antiparallel with β3 of Rmo2-HMA, forming a cross-family hydrogen-bond ladder mediated
predominantly by **main-chain contacts**. Secondary side-chain interactions involve
residues **D15 and K18** from the Rmo2-HMA α1 helix (~α1 C-terminal portion) with
**E42 and S48** from PBY2. Buried interface area: not explicitly reported in the main
text; based on the structure cartoon in Fig. 1C of Yu et al. 2025 the interface appears
smaller than the canonical AVR-Pik interfaces (~1000 Å²).

**Comparison with other binding modes.** Three HMA binding surfaces are now known:
(i) β4/β3 pairing: AVR-Pik effectors + Pikp-HMA (Maqbool 2015 *eLife*; De la
Concepcion 2018 *Nat Plants*); (ii) β2/β2 (αA face) pairing: AVR-Pia and AVR1-CO39
+ RGA5-HMA and Pikp-HMA (Varden 2019 *JBC*; Guo 2018 *PNAS*); (iii) **β3/β2
pairing**: PBY2 + Rmo2-HMA (Yu et al. 2025). PWT7 binds Rwt7-HMA via a fourth,
distinct mode involving extensive side-chain contacts across two faces of Rwt7-HMA
(Yu et al. 2025). The β3 face used by PBY2 is flanked by the AVR-Pik face (β4) on
one side and the AVR-Pia face (β2/αA) on the other; it is spatially accessible in the
native Pikp-1 HMA structure.

**Receptor-side residues (Rmo2-HMA domain numbering / inferred Pikp-1 HMA pipeline position):**

| Rmo2-HMA residue | Identity | Inferred Pikp-1 pipeline equiv. | Secondary element | Contact type | PBY2 partner |
|---|---|---|---|---|---|
| β3 backbone (~29–38) | (strand, multiple residues) | ~pipeline 43–49 | β3 strand | **Main-chain H-bonds** (dominant) | β2 of PBY2 (backbone) |
| **D15** | Asp | ~pipeline 18–20 | α1 (C-terminal portion) | **Side-chain H-bond / salt bridge** | PBY2 **E42** |
| **K18** | Lys | ~pipeline 21–22 | α1 / α1-β2 junction | **Side-chain H-bond** | PBY2 **S48** |

*Caveat: inferred equivalences between Rmo2-HMA and Pikp-1 HMA positions are based
on fold alignment of the shared βαββαβ topology (1–2 residue uncertainty; confirm in
ChimeraX when 9TFO is open).*

**Effector-side residues (PBY2 precursor numbering):**

| PBY2 residue | Identity | MAX fold element | Contact type | HMA partner | Notes |
|---|---|---|---|---|---|
| β2 backbone (~35–48) | (strand) | β2 | Main-chain H-bonds | β3 of Rmo2-HMA | Dominant contacts |
| **E42** | Glu | β2 (within or near β2) | Side-chain salt bridge | D15 of Rmo2-HMA | **E42R abolishes binding to Rmo2-HMA in vitro and recognition in planta** (Yu et al. 2025) |
| **S48** | Ser | β2-β3 loop region | Side-chain H-bond | K18 of Rmo2-HMA | Contributing but not individually tested for essentiality |

**Specificity role.** Rmo2-HMA has nanomolar affinity for PBY2 (K_D = 53.5 nM by ITC;
Yu et al. 2025) but shows **no binding to PWT7** (no peak shift in aSEC; flat ITC
curve). Pikp-1 HMA presumably shows no binding to PBY2 (untested directly, but
consistent with the lack of cross-recognition between Pikp-1 and barley blast in planta).
The β3 face of Pikp-1 HMA is currently inert to PBY2 — this campaign's goal is to
introduce recognition on this face through resurfacing.

**Design implication.** **Primary redesign target: β3 of Pikp-1 HMA (pipeline ~43–49)
plus the flanking β2-β3 and β3-α2 loops.** The backbone geometry of β3 determines
the main-chain H-bond ladder to PBY2 β2. The adjacent α1 positions equivalent to
Rmo2-HMA D15/K18 (~pipeline 18–22) could additionally be flexed to allow side-chain
redesign, but because the Yu et al. 2025 interface is "predominantly main-chain," the
β3 strand is the higher-priority target.

This redesign region (β3, ~pipeline 43–49) **has not been flexed in any sibling
pikp1\_\* campaign**: the avrpikf campaign anchors it (positions 43–78 are in the C-terminal
anchor), and the avrpia campaign also anchors it (positions 43–78 anchored). This
campaign requires a distinct contig from all siblings.

### 4.2 Non-contact regions — anchor

- **N-terminal scaffold (pipeline 1–~41):** β1, α1, β2 — not the PBY2 contact face.
  Includes the α1 positions equivalent to D15/K18 of Rmo2-HMA; see §6 for whether to
  include these in a flex block.
- **C-terminal scaffold (pipeline ~50–78):** β3-α2 loop, α2, β4, C-terminal tail — not
  part of the PBY2/Rmo2-HMA interface. Also includes the canonical AVR-Pik Interface 3
  (positions 69–78), which is irrelevant here.

**Constitutive-activity caveat.** β3 of Pikp-1 HMA (pipeline ~43–49) lies between the
β2/β3 region (Interface 2 of AVR-Pik) and α2 (structural scaffold). Changing β3
backbone/side-chain identity should not directly destabilise the HMA fold (β3 is a
surface-exposed strand) nor should it affect conserved inward-facing residues that are
known autoactivity determinants. This risk level is low; flag any designs with unexpectedly
many changes to β2 or α2-adjacent residues.

---

## 5. Numbering translation

### 5.1 Conventions

| Convention | What it means | Offset |
|---|---|---|
| **Pikp-1 full-length (literature)** | Pikp-1 residues 1–1142; HMA 186–263 | — |
| **Pipeline numbering (chain A)** | 1–78; position = full-length − 185 | pipeline = literature − 185 |
| **Rmo2-HMA domain numbering** | 1–77; the N-terminal HMA domain of barley TKP Rmo2 | — |
| **PBY2 precursor numbering** | 1–85; signal peptide 1–20; mature 21–85 | — |
| **9TFO chain assignments (verified)** | Chains A and B = Rmo2-HMA (residues 2–76, 1-based domain numbering); chains C and D = PBY2 (chain C 24–84, chain D 22–84, both in PBY2 precursor numbering with signal peptide 1–20 and mature 21–85) | E42 verified at chain D residue 42 = GLU |

### 5.2 Receptor key-residue translation table

Pipeline positions are given in the context of the **PBY2 binding face**. The β3 strand is
the primary contact; its exact boundaries in pipeline numbering require secondary-structure
assignment on the AF3 model in ChimeraX.

| Pikp-1 residue (full) | Pipeline pos. | Rmo2-HMA structural equivalent | Role |
|---|---|---|---|
| ~Lys-203 / Thr-204 | ~18–19 | Rmo2-HMA D15 (α1) | Side-chain contact to PBY2 E42 (inferred by fold alignment; verify) |
| ~Arg-206 / Val-207 | ~21–22 | Rmo2-HMA K18 (α1/β2 junction) | Side-chain contact to PBY2 S48 (inferred; verify) |
| ~228 (β3 start) | ~43 | Rmo2-HMA β3 start (~29) | β3/β2 main-chain H-bond ladder begins |
| ~234 (β3 end) | ~49 | Rmo2-HMA β3 end (~37) | β3/β2 main-chain H-bond ladder ends |

*All Pikp-1 HMA equivalences to Rmo2-HMA are inferred from fold alignment of the shared
βαββαβ topology. Offsets are ±1–2 residues until confirmed by ChimeraX superposition of
AF3 chain A onto 9TFO chain A.*

### 5.3 Effector chain assignment and numbering — confirmed

**Confirmed by direct inspection of `experiments/inputs/structures/9TFO_refmac_final.pdb`:**

- 9TFO chains **A and B** are both Rmo2-HMA (residues 2–76 in 1-based HMA domain
  numbering; identical 75-residue sequence). Asymmetric-unit dimer.
- 9TFO chains **C and D** are both PBY2 (chain C residues 24–84, chain D residues 22–84).
  Both use **PBY2 precursor numbering** directly: residue 42 = GLU (E42, the binding-critical
  glutamate); residue 48 = SER (S48). Chain D is the more complete copy.

The pipeline-input source for chain B (effector) is **9TFO chain D**. No renumbering of
PBY2 is required for complex assembly.

### 5.4 Effector key-residue table

| PBY2 residue (precursor) | Identity | MAX fold element | Binding role | Notes |
|---|---|---|---|---|
| β2 backbone (~35–48 precursor) | (strand) | β2 | Main-chain H-bond ladder to HMA β3 | Dominant interaction mode |
| **E42** | Glu | β2 | Salt bridge to Rmo2-HMA D15 | **E42R fully abolishes binding and resistance** (Yu et al. 2025) |
| **S48** | Ser | β2/β3 region | H-bond to Rmo2-HMA K18 | Secondary contact |

---

## 6. Contig string design

All contigs are written in **pipeline numbering** (chain A 1–78). Chain B in the contig
syntax is the fixed PBY2 effector chain.

### 6.1 Design principle from §4

The β3 strand of Pikp-1 HMA (~pipeline 43–49) is the primary redesign target. The
adjacent β2-β3 and β3-α2 loops provide geometry for main-chain H-bond pairing with
PBY2 β2. The rest of the HMA (β1, α1, β2, α2, β4, C-tail) is anchored.

Unlike the AVR-Pik campaigns (Interface 2 at ~33–49 and Interface 3 at ~69–78) and
the AVR-Pia campaign (α1/β2 face at ~14–42), the PBY2 campaign focuses on β3
(~43–49). The anchor blocks are therefore:
- **N-terminal anchor (1–~41):** β1 + α1 + β2 (non-contact face)
- **C-terminal anchor (~50–78):** α2 + β4 + C-terminal tail (non-contact face)

Note that the α1 positions equivalent to Rmo2-HMA D15/K18 (~pipeline 18–22) are in
the N-terminal anchor in the primary contig. The contacts they make with PBY2 E42/S48
are "limited" (Yu et al. 2025, "only limited polar side-chain interactions"), so anchoring
them is a reasonable primary choice; an alternative aggressive contig flexes them.

### 6.2 Primary contig

```
A1-41/10-18/A52-78 B
```

**Pikp-1 full-length (literature) cross-reference**: `A186-226/10-18/A237-263 B` —
i.e. anchor 186–226, flex 10–18 (over native residues 227–236), anchor 237–263. Use the
pipeline-numbered form for RFDiffusion; use the literature form when checking positions
against published Pikp-HMA figures.

| Segment | Pipeline positions | Pikp-1 residues | Native length | Sampled length | Structural element | Role |
|---|---|---|---|---|---|---|
| `A1-41` | 1–41 | 186–226 | 41 | anchor | β1 + α1 + β2 | **Anchored.** Non-PBY2-contact face. Includes α1 equivalents of Rmo2-HMA D15/K18 as native side chains. |
| `10-18` | ~42–51 (flex) | ~227–236 | ~10 | 10–18 | β2-β3 loop + β3 + β3-α2 loop | **Flexible.** Encompasses the β3 strand (~43–49) and immediately flanking loops. Native length ~10 residues; sampling 10–18 (±4 from native) gives RFDiffusion freedom to reshape the strand geometry. |
| `A52-78` | 52–78 | 237–263 | 27 | anchor | α2 + β4 + C-terminal tail | **Anchored.** Structural scaffold and canonical AVR-Pik Interface 3. |
| `B` | full chain B | PBY2 | — | — | MAX effector | **Fully fixed.** |

**Why this contig follows from §4.** The flex block captures the entire β3 strand and its
flanking loops, which are the primary PBY2-contact elements on Rmo2-HMA. The anchor
blocks preserve the HMA fold and the non-PBY2-contact surfaces. The length range
(10–18) gives RFDiffusion room to vary both sequence and local strand geometry, which
is important because PBY2 does not currently bind Pikp-1 HMA and the redesign starts
from a non-binding baseline.

**Boundary check (verify in ChimeraX before committing):**
- **41/42 boundary (anchor → flex).** The β2 strand ends and the β2-β3 loop begins
  around Pikp-1 residue ~228 (pipeline ~43). Position 41 should fall at the end of β2 or
  in the β2-β3 loop — not mid-strand. If β2 extends to position 43, shift the boundary to
  42 or 43.
- **51/52 boundary (flex → anchor).** The β3-α2 loop and start of α2 begin around
  pipeline ~49–54. Position 52 should fall at the end of the β3-α2 loop, not mid-helix.
  Confirm in ChimeraX.

### 6.3 Alternative contig — conservative, native-length flex

```
A1-41/10-10/A52-78 B
```

**Literature cross-reference**: `A186-226/10-10/A237-263 B`.

Fixed at native segment length (10 residues). Only sequence changes sampled; no
length variation. Useful as a comparator to the primary contig to assess whether length
flexibility or sequence changes drive any phenotype difference.

### 6.4 Alternative contig — aggressive, also flex α1 side-chain contacts

```
A1-17/5-7/A23-41/10-18/A52-78 B
```

**Literature cross-reference**: `A186-202/5-7/A208-226/10-18/A237-263 B`.

| Segment | Pipeline positions | Sampled length | Role |
|---|---|---|---|
| `A1-17` | 1–17 | anchor | β1 + early α1 |
| `5-7` | ~18–22 (flex) | 5–7 | **Flexible α1 side chains.** Covers the region equivalent to Rmo2-HMA D15 and K18 (~pipeline 18–22, native ~5 residues). |
| `A23-41` | 23–41 | anchor | Late α1 + β2 |
| `10-18` | ~42–51 (flex) | 10–18 | **Flexible β3** (as primary contig) |
| `A52-78` | 52–78 | anchor | α2 + β4 + C-tail |

Flexes both the side-chain contact region of α1 (D15/K18 equivalents, ~pipeline 18–22)
and β3 (~42–51). Adds complexity but allows RFDiffusion to co-optimize the side-chain
and backbone contacts simultaneously. Use if the primary contig produces designs that
retain the β3 geometry but fail to recreate the D15/K18-equivalent contacts.

**Boundary caveat:** the 17/18 and 22/23 boundaries both fall within α1 (mid-helix).
This is a deliberate trade-off for accessing the contact region; the user must verify
in ChimeraX that these boundaries do not fall on a kinked or critical structural position.

### 6.5 Alternative contig — tight flex with narrow range

```
A1-41/8-12/A50-78 B
```

**Literature cross-reference**: `A186-226/8-12/A235-263 B`.

Narrow range (8–12, ±2 from native). Conservative search: reduces sampling breadth.
Useful if broad sampling in the primary contig produces structurally incoherent designs.

### 6.6 Alternative contig — generous flex with extended β3 and flanking regions

```
A1-40/12-24/A53-78 B
```

**Literature cross-reference**: `A186-225/12-24/A238-263 B`.

Wider flex range (12–24, -3/+11 from native ~10) and slightly extended at both
boundaries. Useful if the primary contig underperforms on the β3/PBY2 geometry due to
insufficient length freedom.

---

## References

- Asuke S, Tagle AG, Hyon G-S, Koizumi S, Murakami T, Horie A, Niwamoto D, Katayama E,
  Shibata M, Takahashi Y, Islam MT, Matsuoka Y, Yamaji N, Shimizu M, Terauchi R, Hisano H,
  Sato K, Tosa Y. *Evolution of HMA-integrated tandem kinases accompanied by expansion of
  target pathogens.* *bioRxiv* (2025). doi:10.64898/2025.12.15.692859. **[from
  experiments/inputs/context/ as asuke_2025_bioarxiv.pdf]**
- De la Concepcion JC, Franceschetti M, Maqbool A, Saitoh H, Terauchi R, Kamoun S,
  Banfield MJ. *Polymorphic residues in rice NLRs expand binding and response to effectors of
  the blast pathogen.* Nature Plants 4, 576–585 (2018). doi:10.1038/s41477-018-0194-x.
  **[from experiments/inputs/context/ as delaconcepcion_2019_natplants.pdf]**
- Guo L, Cesari S, de Guillen K, Chalvon V, Mammri L, Ma M, Meusnier I, Bonnot F, Padilla A,
  Peng Y-L, Liu J, Kroj T. *Specific recognition of two MAX effectors by integrated HMA
  domains in plant immune receptors involves distinct binding surfaces.* Proceedings of the
  National Academy of Sciences 115, 11637–11642 (2018). doi:10.1073/pnas.1810705115.
  **[from experiments/inputs/context/ as guo_2018_pnas.pdf]**
- Maqbool A, Saitoh H, Franceschetti M, Stevenson CEM, Uemura A, Kanzaki H, Kamoun S,
  Terauchi R, Banfield MJ. *Structural basis of pathogen recognition by an integrated HMA
  domain in a plant NLR immune receptor.* eLife 4, e08709 (2015).
  doi:10.7554/eLife.08709.
- Nga NTT, Inoue Y, Chuma I, Hyon G-S, Okada K, Vy TTP, Kusaba M, Tosa Y. *Identification
  of a novel locus Rmo2 conditioning resistance in barley to host-specific subgroups of
  Magnaporthe oryzae.* Phytopathology 102, 674–682 (2012). doi:10.1094/PHYTO-11-11-0297-R.
- Reveguk T, Fatiukha A, Potapenko E, Reveguk I, Sela H, Klymiuk V, Li Y, Pozniak C,
  Wicker T, Coaker G, Fahima T. *Tandem kinase proteins across the plant kingdom.* Nature
  Genetics 57, 254–262 (2025). doi:10.1038/s41588-024-02032-x.
- Varden FA, Saitoh H, Yoshino K, Franceschetti M, Kamoun S, Terauchi R, Banfield MJ.
  *Cross-reactivity of a rice NLR immune receptor to distinct effectors from the rice blast
  pathogen Magnaporthe oryzae provides partial disease resistance.* Journal of Biological
  Chemistry 294, 13006–13016 (2019). doi:10.1074/jbc.RA119.007730.
- **Yu DS, Zdrzałek R, Katayama E, Akiyama H, Daykin L, Williams NJ, Goodridge I, Asuke S,
  Banfield MJ.** *Engineering plant tandem kinase immune receptors expands effector recognition
  profiles.* *bioRxiv* (2025). doi:10.64898/2025.12.15.694194. **[from
  experiments/inputs/context/ as yu_2025_bioarxiv.pdf]**
- UniProt entry: **E9KPB5** (Pikp-1, *Oryza sativa* subsp. japonica, 1142 aa; HMA 188–257).
- GenBank accessions: **LC905591–LC905597** (Asuke et al. 2025; includes PBY2 gene sequence).
- PDB entries: **9TFO** (Rmo2-HMA/PBY2, 2.2 Å, HPUB), **9TFP** (Rwt7-HMA/PWT7, 1.6 Å,
  HPUB), **9TFQ** (Rmo2-HMA+/PBY2/PWT7 tripartite, 1.4 Å, HPUB), **9TFS**
  (Rwt7-HMA+/PBY2, 1.8 Å, HPUB), **9TFT** (Rwt7-HMA+/PWT7, 0.99 Å, HPUB).
