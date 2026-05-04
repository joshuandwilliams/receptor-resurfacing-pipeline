# Literature sweep — pikp1_pwt3 binder-resurfacing campaign

## 1. Summary

This campaign aims to redesign the binding face of the Pikp-1 integrated HMA domain
(UniProt **E9KPB5**, residues 186–263) to gain recognition of **PWT3**, an avirulence
effector from *Pyricularia oryzae* (barley blast pathotype, isolate Br58; DDBJ accession
range LC202650–LC202657 for various allele types; UniProt **A0A223ZP76** for the B-type
from Br48). The campaign goal is de novo gain of recognition: Pikp-1 HMA currently shows
no known binding to PWT3, and PWT3 is not naturally recognised by any HMA-containing
immune receptor — its native recognising gene in wheat, **Rwt3** (syn. Rmg6), is a
**canonical CC-NLR without an integrated HMA domain** (Arora et al. 2023 *Nature Plants*).
**⚠ Campaign-wide flag:** The prior knowledge assumption that PWT3 is a MAX-fold effector
is **contradicted by the published literature**: Inoue et al. 2017 *Science* explicitly state
that PWT3 "lacked similarity to known proteins or protein domains", and Le Naour-Vernet
et al. 2023 *PLoS Pathogens* classify PWT3 as a "non-MAX effector." No deposited
structure for PWT3 exists. **The structural input for this campaign cannot be finalised
until an AF3 monomer prediction of PWT3 is run and its fold validated.** If AF3 predicts
a MAX-like β-sandwich fold, complex assembly can proceed by ChimeraX alignment onto
a Pik-HMA/MAX-effector template (6G10); if not, AF3-multimer is the only route (rung 7).
The receptor chain uses the shared rung-5 AF3 monomer prediction of Pikp-1 HMA(186–263).
Provisional contig strings are provided in §6 conditional on the MAX-fold assumption;
they must be revised once the AF3 prediction is validated.

---

## 2. Available structures

Pipeline numbering (chain A 1–78, position 1 = Pikp-1 residue 186, position 78 =
Pikp-1 residue 263) is used internally only; all table entries use literature numbering.

### 2.1 Native Pikp-1 HMA structures (no engineering mutations)

The full catalogue is identical to sibling pikp1\_\* campaigns and is documented in
`pikp1_avrpikf_literature_sweep.md`. Summary of key entries:

| PDB | Complex contents | Resolved ATOM range | Resolution | Reference | Notes for this campaign |
|-----|------------------|--------------------|------------|-----------|-------------------------|
| 5A6P | Pikp-HMA apo | 186–258 | 2.10 Å | Maqbool 2015 *eLife* | C-term truncated. |
| 5A6W | Pikp-HMA / AVR-PikD | 186–258 | 1.60 Å | Maqbool 2015 *eLife* | First Pik-HMA/MAX complex; canonical AVR-Pik binding face (β4/β3). |
| **6G10** | **Pikp-HMA / AVR-PikD** | **186–263** | **1.35 Å** | **De la Concepcion 2018 *Nat Plants*** | **Preferred alignment template if PWT3 proves MAX-like**. Canonical AVR-Pik binding face. C-terminal residues disordered per prior inspection. |
| 6G11 | Pikp-HMA / AVR-PikE | 186–263 | 1.90 Å | De la Concepcion 2018 *Nat Plants* | Related complex; same binding face as 6G10. |
| 6Q76 | Pikp-HMA / AVR-Pia | 186–263 | 1.90 Å | Varden 2019 *JBC* | AVR-Pia face (different from AVR-Pik); for context only. |

**Disqualifying constraint (shared with all pikp1\_\* campaigns):** All deposits truncate
or disorder C-terminal residues (Pikp-1 residues ~261–263, pipeline positions 76–78).
Chain A uses AF3 rung 5.

### 2.2 Engineered or variant Pikp-1 / Pik-allelic HMA structures

Full table in `pikp1_avrpikf_literature_sweep.md`. Not directly relevant to this
campaign. Key note: 7QPX/7QZD (Pikp-SNK-EKE/AVR-PikF) demonstrate that Interface
3 resurfacing can expand recognition of a non-cognate MAX effector — the closest
engineering precedent for the concept underlying this campaign.

### 2.3 Effector structures — PWT3 and related MAX-fold templates

| PDB | Effector | Complex contents | Method | Resolution | Reference | Campaign relevance |
|-----|----------|-----------------|--------|------------|-----------|-------------------|
| **None** | **PWT3 (Ao-type / B-type)** | — | — | — | Inoue et al. 2017 *Science* | **No deposited structure. Only computationally predicted.** Fold unknown; must be determined by AF3. |
| 5A6W | AVR-PikD | Pikp-HMA | X-ray | 1.60 Å | Maqbool 2015 | Canonical MAX-effector structure with AVR-Pik binding face; reference if PWT3 is MAX-like. |
| 6G10 | AVR-PikD | Pikp-HMA | X-ray | 1.35 Å | De la Concepcion 2018 | Best Pikp-HMA / MAX complex template for alignment (if PWT3 is MAX-like). |
| 2MYW | AVR-Pia | None (free) | NMR | — | de Guillen 2015 *PLoS Pathog* | Free MAX-effector structure for fold comparison during AF3 validation. |
| 7B1I | AVR-PikF | OsHIPP19-HMA | X-ray | 1.90 Å | Maidment 2021 *JBC* | Alternative MAX-effector template (canonical binding face). |

**No PWT3 crystal structure or NMR structure exists.** The campaign is entirely dependent
on AF3 structure prediction.

### 2.4 Complexes between receptor and effectors (alignment template candidates)

| PDB | Receptor | Effector | Resolution | Reference | Status for this campaign |
|-----|----------|----------|------------|-----------|--------------------------|
| **6G10** | Pikp-HMA | AVR-PikD | 1.35 Å | De la Concepcion 2018 | **Preferred template if PWT3 is MAX-like** and binds the canonical AVR-Pik face (β4 of HMA / β3 of MAX) |
| 6FU9 | Pikm-HMA | AVR-PikD | 1.20 Å | De la Concepcion 2018 | Highest-resolution Pik-HMA/MAX-effector complex; secondary validation template |
| 7B1I | OsHIPP19-HMA | AVR-PikF | 1.90 Å | Maidment 2021 *JBC* | AVR-PikF binding face; alternative template |

---

## 3. Structural-input decision

### 3.1 Fallback ladder analysis

**Rung 1: native receptor / target-effector complex crystal.** No Pikp-1 HMA / PWT3
complex. → **Not available.**

**Rung 2: native receptor / related-effector complex crystal.** PWT3 has no deposited
structure and no complex with any HMA domain. No related effector complex that could
substitute exists because PWT3 is not known to interact with any HMA protein. →
**Not available.**

**Rung 3 and 4.** Not applicable. No PWT3 structure of any kind exists.

**Rung 5 (AF3 monomer of native receptor — chain A).** Used, as for all pikp1\_\*
campaigns. AF3 monomer prediction of Pikp-1 HMA (UniProt E9KPB5, residues 186–263),
pipeline numbering 1–78.

**Rung 6 (AF3 monomer of effector — chain B).** **Required.** PWT3 has no deposited
structure; chain B must come from AF3 monomer prediction of PWT3.

**Rung 7 (AF3-multimer of the complex).** **Fallback, and recommended to run
regardless.** If the AF3 monomer prediction of PWT3 does not produce a recognisable
MAX-like fold, or if the resulting ChimeraX alignment generates an implausible binding
pose, AF3-multimer of Pikp-1 HMA + PWT3 is the only route.

### 3.2 The non-MAX fold problem

The published literature explicitly states that PWT3 **lacked similarity to known
proteins or protein domains** (Inoue et al. 2017 *Science*; reiterated in Inoue et al.
2021 *New Phytologist*), and Le Naour-Vernet et al. 2023 *PLoS Pathogens* classify
PWT3 as "non-MAX." This contradicts the prior knowledge for this campaign.

Two scenarios must be considered:

**Scenario A: AF3 predicts a MAX-like fold for PWT3.** Despite the absence of sequence
similarity to MAX effectors, AF3 might predict a β-sandwich fold structurally similar to
the MAX family. This is plausible because MAX effectors are defined by structure, not
sequence, and the sequence space of β-sandwich folds is vast. If pLDDT is high (>70
over the folded core) and RMSD to 2MYW (AVR-Pia, NMR) or AVR-PikD in 5A6W is
<3 Å over ≥40 equivalent Cα atoms, treat PWT3 as functionally MAX-like and proceed
with the ChimeraX alignment route.

**Scenario B: AF3 predicts a non-MAX fold.** If the predicted fold is clearly different
from the β-sandwich (e.g. helical bundle, TIM barrel, or disordered), or if pLDDT is
low across the mature chain, the ChimeraX-alignment route is invalid. AF3-multimer
of Pikp-1 HMA + PWT3 is then the only structural basis for complex assembly.

**Scenario C (worst case): AF3 confidence is too low to determine either.** Mark as
campaign-blocked pending manual analysis of the AF3 multimer prediction.

### 3.3 Chosen structural input

| Chain | Source | Rung | Sequence input | Notes |
|-------|--------|------|----------------|-------|
| **A** (Pikp-1 HMA) | AF3 monomer of UniProt E9KPB5 (186–263) | 5 | E9KPB5 residues 186–263 | Shared with all pikp1_* campaigns; pipeline 1–78 |
| **B** (PWT3) | AF3 monomer prediction | 6 | **Ao-type PWT3** mature chain (DDBJ LC202650 or equivalent) — see §3.4 | Conditional on fold validation; fallback is AF3-multimer (rung 7) |

### 3.4 PWT3 sequence for AF3 input

Two PWT3 variants are relevant:

- **Ao-type (functional, avirulent):** from Br58 and Lolium isolates; triggers Rwt3-mediated
  resistance in wheat. Sequence in DDBJ under LC202650–LC202657 (Inoue et al. 2017
  *Science*). **This is the design target** — the form Pikp-1 HMA needs to recognise.
  Mature chain: residues after signal peptide cleavage (cleavage at approximately
  residue 18–20 by SignalP prediction). Use 123-residue mature chain.
- **B-type (virulent):** from Br48; 12 amino acid substitutions vs Ao-type. UniProt
  A0A223ZP76 (Br48 strain, 141 aa, signal peptide 1–18, mature 19–141 = 123 aa).
  Suitable as a proxy for fold prediction (the 12-residue divergence will not change
  the overall fold), but the Ao-type is preferred for design.

**Preferred input:** Ao-type PWT3 mature chain from DDBJ (LC202650 or one of the
LC202650–LC202657 accessions encoding the Ao-type). If the Ao-type sequence is not
conveniently retrievable, UniProt A0A223ZP76 mature chain (residues 19–141) is an
acceptable substitute for fold prediction.

### 3.5 Complex assembly recipe (CONDITIONAL on AF3 validation)

**If Scenario A (MAX-like fold, ≥40 Cα RMSD < 3 Å to 2MYW or 5A6W):**

1. Superpose AF3 PWT3 monomer onto the effector chain of **6G10** (AVR-PikD in
   Pikp-HMA/AVR-PikD; chain B of 6G10) using Cα alignment over the MAX core.
2. Adopt the Pikp-HMA coordinates from 6G10 to establish the HMA/effector geometry.
3. Superpose the AF3 Pikp-1 HMA (chain A, pipeline 1–78) onto 6G10 chain A.
4. Combine AF3 chain A + repositioned AF3 PWT3 chain B to produce
   `pikp1_pwt3_complex.pdb`.
5. Verify: are the canonical Pik-HMA Interface 2 and Interface 3 residues
   (from De la Concepcion 2018; Maqbool 2015) geometrically positioned to contact
   the equivalent region of PWT3 β3 in the aligned model?

**If Scenario B (non-MAX fold):**

1. Run AF3-multimer of Pikp-1 HMA (E9KPB5 186–263) + PWT3 mature chain
   (Ao-type, ~123 aa) on HPC.
2. Inspect the predicted interface in ChimeraX. Determine which HMA face contacts PWT3.
3. Extract the best-scoring complex model as the structural input.
4. Return to Part 2 for interface analysis — §4 and §6 will need to be revised based
   on the multimer result.

**Validation step (both scenarios):** after assembly, inspect whether Pikp-1 HMA
α-helices and β-strands in the proposed binding interface show good backbone geometry.
Check pLDDT >70 in the predicted binding region of the AF3 PWT3 chain.

---

## 4. Binding interfaces

**⚠ This entire section is INFERRED and CONDITIONAL.** PWT3 has no deposited structure
and no characterised interaction with any HMA domain. The natural resistance gene
Rwt3 is a CC-NLR without HMA; PWT3 is not naturally detected via an HMA bait. All
claims in this section derive from:
(a) homology reasoning from related MAX-fold / Pik-HMA complexes (**labelled: INFERRED
from MAX homology**), or
(b) direct literature evidence for PWT3 biology (**labelled: CONFIRMED from PWT3 sources**).

No mutagenesis data mapping PWT3 residues to any HMA-binding interaction exists.
No structural data for any PWT3/HMA complex exists.

### 4.1 Canonical AVR-Pik binding face (Interfaces 1–3) — CONDITIONAL on MAX fold

**Structural description (INFERRED from MAX homology, Scenario A only).** If the AF3
monomer of PWT3 reveals a MAX-fold β-sandwich and the ChimeraX alignment onto 6G10
positions PWT3 in the canonical binding pose (β3 of PWT3 antiparallel to β4 of
Pikp-HMA), the interface would decompose into the same three sub-interfaces defined
by De la Concepcion et al. 2018 for AVR-PikD/Pikp-HMA. A summary of those
interfaces is reproduced below for contig-design purposes (full analysis is in
`pikp1_avrpikf_literature_sweep.md`).

**Interface 1 — N-terminus / β1 area (pipeline 3–10):**
- Structurally conserved across all AVR-Pik variants; not specificity-determining.
- **Design implication:** Anchor. Disrupting this interface would risk destablising the
  HMA fold.

**Interface 2 — β2/β3 region (pipeline ~33–49, Pikp-1 residues ~218–234):**
- Receptor-side (INFERRED from AVR-Pik analogy): Ser-218 (33), Asp-224 (39),
  Glu-230 (45), Val-232 (47) — the canonical Pikp-1 Interface 2 contacts for His-46
  of AVR-Pik effectors.
- Effector-side (UNKNOWN for PWT3): no data. If PWT3 is MAX-like, the equivalent
  of AVR-Pik position 46 in PWT3 would need to be identified from the AF3 model.
  The Ao-type PWT3 and AVR-PikD share no obvious sequence identity at this region.
- **Specificity role (INFERRED):** For AVR-Pik effectors, Interface 2 determines
  allele-specific recognition (Kanzaki 2012 *Plant Journal*; De la Concepcion 2018).
  Whether Interface 2 residues of PWT3 would contact Pikp-1 HMA at all is unknown.
- **Design implication (CONDITIONAL):** If PWT3 binds the AVR-Pik canonical face,
  Interface 2 is a **primary redesign target** because the Pikp-1 Interface 2 residues
  will not match PWT3 without redesign.

**Interface 3 — β4 to C-terminal tail (pipeline ~69–78, Pikp-1 residues ~254–263):**
- Receptor-side (INFERRED): Met-254 (69), Ser-258 (73), Ala-260 (75), Asn-261 (76),
  Lys-262 (77), Asp-263 (78) — the Lys262 pocket.
- Effector-side (UNKNOWN for PWT3): no data. For AVR-PikD, Interface 3 involves
  Glu-53, Tyr-71, Ser-72, Trp-74 forming the "Lys262 pocket." Whether the equivalent
  pocket exists in PWT3's β-strand 4 / C-terminal region is unknown.
- **Specificity role (INFERRED):** For AVR-PikF, the unique Met78Lys substitution in
  the effector disrupts the Lys262 pocket (Longya 2019 *MPMI*). The equivalent
  polymorphism status in PWT3 is unknown.
- **Design implication (CONDITIONAL):** If PWT3 binds the canonical face, Interface 3
  is a **primary redesign target** because all available Pikp-HMA crystals disorder
  exactly this region, and no functional Pik-HMA/PWT3 interaction is known.

### 4.2 Non-canonical binding mode — SCENARIO B

If AF3 predicts a non-MAX fold for PWT3, the canonical Interface 1/2/3 framework
does not apply. The binding face and the relevant Pikp-1 HMA residues would need to
be determined from the AF3-multimer contacts. **No prediction is possible here without
the multimer output.**

### 4.3 Why PWT3 is a de novo engineering target (CONFIRMED from PWT3 sources)

- PWT3 is recognised by Rwt3, a CC-NLR with no HMA integrated domain (Arora et al.
  2023 *Nature Plants*; Inoue et al. 2017 *Science*). Therefore PWT3 has **no natural
  history of HMA binding** to anchor the design on.
- The Ao-type PWT3 (functional avirulence gene, triggers Rwt3 resistance) and the
  B-type (virulent, 12 amino acid substitutions) differ at 12 positions in the mature
  chain (Inoue et al. 2017 *Science*). The substituted positions in the B-type are
  likely surface-exposed on the functional protein, possibly at a site relevant to Rwt3
  recognition — but their relevance to a potential HMA-binding surface on PWT3 is
  unknown.
- No mutagenesis data maps any PWT3 residue to recognition or binding function at the
  molecular level. Inoue et al. 2021 *New Phytologist* focuses entirely on PWT4
  biology and provides no new PWT3 structural or binding data.

---

## 5. Numbering translation

### 5.1 Conventions

| Convention | What it means | Offset |
|---|---|---|
| **Pikp-1 full-length (literature)** | Pikp-1 residues 1–1142; HMA 186–263 | — |
| **Pipeline numbering (chain A)** | 1–78; pipeline = full-length − 185 | pipeline = literature − 185 |
| **PWT3 precursor numbering** | 1–141 (UniProt A0A223ZP76 B-type; Ao-type same length); signal peptide 1–18; mature 19–141 | — |
| **AF3 output numbering (chain B)** | If full precursor (1–141) submitted: AF3 output = precursor numbering directly. If mature chain only (19–141) submitted: output is renumbered 1–123, and precursor residue = AF3 position + 18. | Offset: AF3 mature position + 18 = precursor position |

### 5.2 Receptor key-residue translation (CONDITIONAL on Interface 2/3 being the target)

These positions are inherited from the pikp1\_avrpikf analysis and are relevant only
if PWT3 binds the canonical AVR-Pik face.

| Pikp-1 residue (full) | Identity | Pipeline pos. | Interface | Role (from AVR-Pik analogy — INFERRED) |
|---|---|---|---|---|
| 218 | Ser | 33 | 2 | H-bond to effector His-46 (AVR-PikD); equivalent in PWT3 unknown |
| 224 | Asp | 39 | 2 | Salt bridge to effector Arg-64 (AVR-PikD) |
| 230 | Glu | 45 | 2 | Canonical Pikp/Pikm polymorphism site |
| 232 | Val | 47 | 2 | Hydrophobic contact to effector His-46 |
| 254 | Met | 69 | 3 | β4 main-chain pairing |
| 258 | Ser | 73 | 3 | SNK-EKE site; C-terminal geometry |
| 261 | Asn | 76 | 3 | NK-KE/SNK-EKE engineering site |
| 262 | Lys | 77 | 3 | Central Interface 3 contact (Lys262 pocket) |
| 263 | Asp | 78 | 3 | C-terminal anchor |

### 5.3 Effector numbering (UNKNOWN until AF3 prediction)

Because PWT3 has no deposited structure and no HMA-binding data, the effector-side
residue table cannot be completed. If the AF3 monomer or multimer reveals an interface,
the following placeholder positions can be populated:

| PWT3 residue (precursor) | Identity (Ao-type) | AF3 mature position | Role | Status |
|---|---|---|---|---|
| TBD | TBD | TBD | Predicted from AF3 multimer / Scenario A alignment | **UNKNOWN — requires AF3** |

---

## 6. Contig string design

**⚠ ALL CONTIGS IN THIS SECTION ARE CONDITIONAL on the MAX-fold assumption
(Scenario A). If the AF3 monomer of PWT3 does not produce a MAX-like fold, or
if the AF3-multimer contacts indicate a different binding face, the contigs must be
re-derived from the multimer contact map before use.**

**⚠ Because chain B comes from an AF3 prediction (not a crystal structure), the contact
pattern that drives flex-region selection has inherent uncertainty. The provisional
contigs below are derived from the assumption that PWT3, if MAX-like, docks at the
canonical Pik-HMA AVR-Pik face (Interface 2 + Interface 3) — the same face as
AVR-PikF. Recommend re-running contact derivation from the assembled complex in
ChimeraX once the AF3 prediction is available.**

### 6.1 Design principle (CONDITIONAL — Scenario A, MAX-fold, canonical AVR-Pik face)

If PWT3 is MAX-like and binds the canonical AVR-Pik face, the relevant Pikp-1 HMA
regions to redesign are the same as established for AVR-PikF (sibling pikp1\_avrpikf
campaign): Interface 2 (β2/β3 region, pipeline ~33–49) and Interface 3 (β4/C-tail,
pipeline ~69–78). Interface 1 (β1, pipeline ~1–10) is conserved and should be anchored.
The α2 scaffold helix (pipeline ~60–66) is structural and should also be anchored.

Under this assumption:
- **N-terminal scaffold (positions 1–32):** β1 + α1 + end of α1. *Anchor.*
- **Interface 2 (positions 33–49):** β2 + β2-β3 loop + β3 N-terminus. *Flex.*
- **β3-tail / α2 scaffold (positions 50–68):** β3 remainder + α2. *Anchor.*
- **Interface 3 (positions 69–78):** β4 + C-terminal tail. *Flex.*

### 6.2 Primary contig (CONDITIONAL — Scenario A, canonical face)

```
A1-32/12-20/A50-68/8-14 B
```

**Pikp-1 full-length (literature) cross-reference**: `A186-217/12-20/A235-253/8-14 B` —
i.e. anchor 186–217, flex 12–20 (over native residues 218–234), anchor 235–253, flex 8–14
(over native residues 254–263). Identical anchor/flex layout to the sibling pikp1_avrpikf
primary contig, as expected under the canonical-AVR-Pik-face Scenario A assumption.

| Segment | Pipeline positions | Pikp-1 residues | Native length | Sampled length | Structural element | Role |
|---|---|---|---|---|---|---|
| `A1-32` | 1–32 | 186–217 | 32 | anchor | β1 + α1 + start of β2 | **Anchored.** Preserves Interface 1 and scaffold native. |
| `12-20` | ~33–49 (flex) | ~218–234 | 17 | 12–20 | β2 + β2-β3 loop + β3 N-terminus | **Flexible.** Interface 2 region. INFERRED as redesign target (canonical face assumption). |
| `A50-68` | 50–68 | 235–253 | 19 | anchor | β3 tail + α2 + β3-α2-β4 transition | **Anchored.** α2 scaffold. |
| `8-14` | ~69–78 (flex) | ~254–263 | 10 | 8–14 | β4 + C-terminal tail | **Flexible.** Interface 3 region. INFERRED as redesign target (canonical face assumption). |
| `B` | full chain B | PWT3 (AF3 predicted) | — | — | MAX effector (conditional) | **Fully fixed.** |

**Boundary check (verify after AF3 complex assembly in ChimeraX):**
- Position 32/33: should fall in the β1-α1-β2 loop or at the start of β2. Mid-strand is a concern — verify in ChimeraX.
- Position 49/50: should be in the β3-α2 loop, not mid-strand. Verify in ChimeraX.
- Position 68/69: should be one residue before β4 starts; verify position of Met-254 (pipeline 69) in the AF3 model.

### 6.3 Alternative contig — conservative, Interface 3 only (CONDITIONAL)

```
A1-67/8-14 B
```

**Literature cross-reference**: `A186-252/8-14 B`.

Redesigns only Interface 3 (pipeline 69–78); leaves Interface 2 native. Rationale: if
PWT3's recognition barrier is dominated by the C-terminal tail geometry (analogous to
AVR-PikF Lys78 disrupting the Lys262 pocket), Interface 3 alone may be sufficient.
Under-powered, but included as the minimum intervention.

### 6.4 Alternative contig — more aggressive, extended flex ranges (CONDITIONAL)

```
A1-32/14-25/A50-68/8-18 B
```

**Literature cross-reference**: `A186-217/14-25/A235-253/8-18 B`.

Wider length ranges for both flex blocks; gives RFDiffusion more freedom to search.
Use if the primary contig fails to converge on viable designs.

### 6.5 Alternative contig — fixed-length flex (CONDITIONAL)

```
A1-32/17-17/A50-68/10-10 B
```

**Literature cross-reference**: `A186-217/17-17/A235-253/10-10 B`.

Fixed at native lengths. Only sequence changes sampled; length variation removed.
Useful as a comparator to the primary contig.

### 6.6 Non-MAX alternative — contig from AF3-multimer (Scenario B)

**If Scenario B applies (non-MAX fold confirmed or low pLDDT), this section
replaces §6.2–6.5 entirely.**

After running AF3-multimer of Pikp-1 HMA + PWT3 on HPC:
1. Identify the predicted contact residues on Pikp-1 HMA using qtPISA or equivalent.
2. Map contact residues to pipeline numbering (pipeline = literature − 185).
3. Define flex blocks to encompass the contact region, with anchor blocks preserving
   non-contact faces.
4. Return here and replace the conditional contigs with contact-derived contigs before
   committing to the pipeline run.

The contig format will be the same `A<start>-<end>/<min>-<max>/... B` syntax, but
derived from the actual multimer contact pattern rather than from MAX-effector analogy.

---

## References

- Arora S, Steed A, Goddard R, Gaurav K, O'Hara T, Schoen A, Rawat N, Elkot AF,
  Korolev AV, Chinoy C, Nicholson MH, Asuke S, Antoniou-Kourounioti R, Steuernagel B,
  Yu G, Awal R, Forner-Martínez L, Wingen L, Baggs E, Clarke J, Saunders DGO,
  Krasileva KV, Tosa Y, Jones JDG, Tiwari VK, Wulff BBH, Nicholson P.
  *A wheat kinase and immune receptor form host-specificity barriers against the blast
  fungus.* Nature Plants 9, 385–392 (2023). doi:10.1038/s41477-023-01357-5.
  **[PMC10027608; open access]**
- De la Concepcion JC, Franceschetti M, Maqbool A, Saitoh H, Terauchi R, Kamoun S,
  Banfield MJ. *Polymorphic residues in rice NLRs expand binding and response to
  effectors of the blast pathogen.* Nature Plants 4, 576–585 (2018).
  doi:10.1038/s41477-018-0194-x. **[from context/ as delaconcepcion_2019_natplants.pdf]**
- de Guillen K, Ortiz-Vallejo D, Gracy J, Fournier E, Kroj T, Padilla A. *Structure
  analysis uncovers a highly diverse but structurally conserved effector family in
  phytopathogenic fungi.* PLoS Pathogens 11, e1005228 (2015).
  doi:10.1371/journal.ppat.1005228.
- **Inoue Y, Vy TTP, Yoshida K, Asano H, Mitsuoka C, Asuke S, Anh VL, Cumagun CJR,
  Chuma I, Terauchi R, Kato K, Mitchell T, Valent B, Farman M, Tosa Y.**
  *Evolution of the wheat blast fungus through functional losses in a host specificity
  determinant.* Science 357, 80–83 (2017). doi:10.1126/science.aam9654.
  **[from context/ as inoue_2017_science.pdf]**
- **Inoue Y, Vy TTP, Tani D, Tosa Y.** *Suppression of wheat blast resistance by an
  effector of Pyricularia oryzae is counteracted by a host specificity resistance gene
  in wheat.* New Phytologist 229, 488–500 (2021). doi:10.1111/nph.16894.
  **[from context/ as inoue_2021_newphytol.pdf]** *(Note: this paper primarily
  characterises PWT4, not PWT3. The relevant PWT3 content is the
  confirmation that "both genes encoded small secreted proteins with no functional
  domains.")*
- Kanzaki H, Yoshida K, Saitoh H, Fujisaki K, Hirabuchi A, Alaux L, Fournier E,
  Tharreau D, Terauchi R. *Arms race co-evolution of Magnaporthe oryzae AVR-Pik and
  rice Pik genes driven by their physical interactions.* The Plant Journal 72, 894–907
  (2012). doi:10.1111/j.1365-313X.2012.05110.x.
  **[from context/ as kanzaki_2012_plantj.pdf]**
- Le Naour-Vernet M, Were VM, Foster AJ, Langner T, Malmgren A, Harant A, Asuke S,
  Reyes-Avila S, Gupta DR, Jensen C, Ma W, Mahmud NU, Mehebub MS, Mulenga RM,
  Muzahid ANM, Paul SK, Rabby SMF, Rahat AAM, Ryder L, Shrestha R-K, Sichilima S,
  Soanes DM, Singh PK, Bentley AR, Saunders DGO, Tosa Y, Croll D, Lamour KH,
  Islam T, Tembo B, Win J, Talbot NJ, Burbano HA, Kamoun S.
  *Adaptive evolution in virulence effectors of the rice blast fungus Pyricularia oryzae.*
  PLoS Pathogens 19, e1011294 (2023). doi:10.1371/journal.ppat.1011294. **[open access]**
  *(Classifies PWT3 as "non-MAX effector.")*
- Longya A, Chaipanya C, Franceschetti M, Maidment JHR, Banfield MJ, Jantasuriyarat C.
  *Gene Duplication and Mutation in the Emergence of a Novel Aggressive Allele of the
  AVR-Pik Effector in the Rice Blast Fungus.* Molecular Plant–Microbe Interactions 32,
  740–749 (2019). doi:10.1094/MPMI-09-18-0245-R.
  **[from context/ as longya_2019_mpmi.pdf]**
- Maqbool A, Saitoh H, Franceschetti M, Stevenson CEM, Uemura A, Kanzaki H, Kamoun S,
  Terauchi R, Banfield MJ. *Structural basis of pathogen recognition by an integrated
  HMA domain in a plant NLR immune receptor.* eLife 4, e08709 (2015).
  doi:10.7554/eLife.08709.
- Maidment JHR, Franceschetti M, Maqbool A, Saitoh H, Jantasuriyarat C, Kamoun S,
  Terauchi R, Banfield MJ. *Multiple variants of the fungal effector AVR-Pik bind the
  HMA domain of the rice protein OsHIPP19, providing a foundation to engineer plant
  defense.* Journal of Biological Chemistry 296, 100371 (2021).
  doi:10.1016/j.jbc.2021.100371.
- UniProt entry: **E9KPB5** (Pikp-1, *Oryza sativa* subsp. japonica, 1142 aa).
- UniProt entry: **A0A223ZP76** (PWT3, *Pyricularia oryzae* strain Br48 B-type, 141 aa;
  signal peptide 1–18; mature 19–141).
- DDBJ/GenBank accessions **LC202650–LC202657, LC215053, LC215054, LC229726**
  (Inoue et al. 2017; PWT3 and PWT4 allele sequences from various P. oryzae isolates).
