# Literature sweep — pikp1_avrpikf binder-resurfacing campaign

## 1. Summary

This campaign aims to redesign the binding face of the Pikp-1 integrated HMA
domain (UniProt E9KPB5, residues 186–263) so that it recognises the
"stealthy" rice-blast effector AVR-PikF (Longya et al. 2019), which evades
all naturally-occurring Pik alleles tested to date by combining
Interface-2-disrupting substitutions at AVR-Pik positions 46/47/48 with a
unique Met78Lys substitution that disrupts Interface 3. Many high-resolution
crystal structures of the Pik-HMA / AVR-Pik system exist (Maqbool et al.
2015, De la Concepcion et al. 2018, Maidment et al. 2021 *JBC*, Maidment
et al. 2023), but every deposited Pikp-1-HMA structure either truncates or
disorders one or more residues in the C-terminal Interface 3 region, so the
campaign drops to **fallback rung 5** — an AlphaFold3 monomer prediction of
native Pikp-1 HMA(186–263) — and pairs it with AVR-PikF coordinates lifted
from PDB **7B1I** (OsHIPP19-HMA / AVR-PikF; Maidment et al. 2021 *JBC*).
The **primary redesign targets are Interface 2** (β2/β3 region; Pikp-1
residues 218–232 = pipeline positions 33–47) **and Interface 3** (β4 to
C-terminus; Pikp-1 residues 254–263 = pipeline positions 69–78). There are
**no campaign-blockers**: a viable structural input, a usable complex
template, a high-quality effector structure, and an established mutagenesis
literature are all present.

## 2. Available structures

Pipeline numbering (chain A 1–78, where 1 = Pikp-1 residue 186 and
78 = residue 263) is *only* used internally by this campaign; the tables in
this section use literature numbering as deposited in each PDB entry.

### 2.1 Native Pikp-1 HMA structures (no engineering mutations)

| PDB | Complex contents | Construct (full-length residues) | Resolved ATOM range | Resolution | Reference | Notes for this campaign |
|-----|------------------|-----------------------------------|---------------------|------------|-----------|-------------------------|
| 5A6P | Pikp-HMA apo (dimer)              | 186–258                                           | 186–258                       | 2.10 Å | Maqbool 2015 *eLife* | Construct truncates before Lys262, so no Interface 3 packing visible. |
| 5A6W | Pikp-HMA / AVR-PikD               | 186–258                                           | 186–258                       | 1.60 Å | Maqbool 2015 *eLife* | First Pik-HMA / AVR-Pik complex; same construct truncation as 5A6P — Interface 3 not represented. |
| 6G10 | Pikp-HMA / AVR-PikD               | 186–263 (5-residue C-terminal extension)          | per user inspection of ATOM records, C-terminal residues critical to the binding interface are missing | 1.35 Å | De la Concepcion 2018 *Nature Plants* | The "C-terminal extension" construct that revealed Interface 3, but the user has verified that the deposited ATOM records still do not resolve all residues required for the campaign's design region. |
| 6G11 | Pikp-HMA / AVR-PikE               | 186–263                                           | as 6G10                       | 1.90 Å | De la Concepcion 2018 *Nature Plants* | Pikp does not respond to AVR-PikE in planta; complex captured only in vitro and shows interface rearrangements (Asn46 of effector flipped out of binding pocket, Gln259/Ala260 of receptor "looped out"). |
| 6Q76 | Pikp-HMA / AVR-Pia                | 186–263                                           | through 263 per Varden 2019    | 1.90 Å | Varden 2019 *JBC* | AVR-Pia binds a *different face* of Pikp-HMA (adjacent to α1 / β2, ~460 Å² interface) — not the AVR-Pik face used here. Useful only as a control. |

**Disqualifying constraint (carried forward from the user's prior inspection
of the ATOM records).** All five entries either truncate or disorder
C-terminal residues that lie inside Interface 3 (Lys262 and the surrounding
β4 / C-terminal tail). 5A6P and 5A6W truncate before residue 259. 6G10,
6G11 and 6Q76 have the long construct but the deposited ATOM coordinates
do not resolve every C-terminal residue the campaign needs to anchor or
redesign. **No native Pikp-1 HMA crystal provides a complete Interface 3.**

### 2.2 Engineered or variant Pikp-1 / Pik-allelic HMA structures

| PDB | Complex contents | Mutations vs Pikp-1 wild-type | Construct | Resolution | Reference | Notes for this campaign |
|-----|------------------|-------------------------------|-----------|------------|-----------|-------------------------|
| 6R8K | Pikp-HMA-NK-KE / AVR-PikD         | Asn261Lys + Lys262Glu (Pikp full-length numbering) | 186–263   | 1.60 Å | De la Concepcion 2019 *eLife* | "NK-KE" engineering. Mutations sit at Interface 3 (positions 76, 77 in pipeline numbering). |
| 6R8M | Pikp-HMA-NK-KE / AVR-PikE         | Asn261Lys + Lys262Glu        | 186–263   | 1.85 Å | De la Concepcion 2019 *eLife* | NK-KE expands recognition to AVR-PikE. |
| 7A8W | Pikp-HMA-NK-KE / AVR-PikC         | Asn261Lys + Lys262Glu        | 186–263   | 2.15 Å | Maidment 2021 *PLoS Pathog* | NK-KE engineering reused as a molecular-replacement template. |
| 7A8X | Pikh-HMA / AVR-PikC               | Pikh natural allele (Asn261Lys vs Pikp at Interface 3) | 186–263   | 2.30 Å | Maidment 2021 *PLoS Pathog* | Pikh-HMA differs from Pikp-HMA by one amino acid (Asn261Lys); no other Pik-allele polymorphism mapped here. |
| 7QPX | Pikp-HMA-SNK-EKE / AVR-PikC       | Ser258Glu + Asn261Lys + Lys262Glu                  | 186–264   | 2.05 Å | Maidment 2023 *eLife* | "SNK-EKE" — the published engineered solution that *does* recognise AVR-PikC and AVR-PikF. The campaign explicitly aims to find alternatives to this solution rather than build on it. |
| 7QZD | Pikp-HMA-SNK-EKE / AVR-PikF       | Ser258Glu + Asn261Lys + Lys262Glu                  | 186–264   | 2.20 Å | Maidment 2023 *eLife* | SNK-EKE bound to the campaign's target effector. Per the user's verified prior knowledge, ATOM records show several C-terminal residues are unresolved (disordered) in the deposited coordinates, so even this construct does not give a complete Interface 3 in usable form. |
| 6FU9 | Pikm-HMA / AVR-PikD               | Pikm natural allele (multiple polymorphisms vs Pikp across all three interfaces) | 186–264 | 1.20 Å | De la Concepcion 2018 *Nature Plants* | Highest-resolution Pik-HMA / AVR-Pik complex known. Pikm-HMA is the template used by De la Concepcion 2018 for naming Interface 1 / 2 / 3 residues; comparison to Pikp-HMA reveals the polymorphic determinants. |
| 6FUB | Pikm-HMA / AVR-PikE               | Pikm natural allele                               | 186–264   | 1.30 Å | De la Concepcion 2018 *Nature Plants* | Captures the conformational rearrangement at Interface 2 when Asn46 (effector) flips out. |
| 6FUD | Pikm-HMA / AVR-PikA               | Pikm natural allele                               | 186–264   | 1.30 Å | De la Concepcion 2018 *Nature Plants* | AVR-PikA shares positions 46/47/48 with AVR-PikF; informative as an Interface-2 analogue for the F variant. |
| 7BNT | Ancestral Pik-1 HMA / AVR-PikD    | Reconstructed ancestral sequence                  | n/a       | 1.32 Å | Białas 2021 *eLife* | Useful for understanding which residues were under positive selection across the Pik family. |
| 8B2R | RGA5-HMA-mut / AVR-PikF           | Six engineered RGA5-HMA mutations                 | n/a       | 1.22 Å | Bentham 2023 *Plant Cell* | The other published precedent for engineering AVR-PikF recognition into a non-cognate HMA scaffold. RGA5 is a distinct rice NLR; not directly usable as a Pikp template but indicates which contact pattern works for AVR-PikF in a Pikm-1-like geometry. |

### 2.3 Effector structures (AVR-Pik family, alone or in complex)

| PDB | Effector | Complex partner | Construct (effector residues, precursor numbering) | Resolution | Reference | Notes for this campaign |
|-----|----------|------------------|-----------------------------------------------------|------------|-----------|-------------------------|
| 7B1I | AVR-PikF | OsHIPP19-HMA (chain B, ATOM residues 1–76; gap at 14–15) | **Chain C ATOM residues 33–113 in AVR-Pik precursor numbering (81 residues, no gaps)** | 1.90 Å | Maidment 2021 *JBC* | **Source of chain C in the assembled complex.** Structure shows AVR-PikF with **K78 in place of M78** (verified: 7B1I chain C residue 78 = LYS). Demonstrates that the AVR-PikF Met78Lys substitution is compatible with Pik-HMA-fold binding *if* the receptor can accommodate Lys78 — OsHIPP19 does, Pikp-1 does not. |
| 7QZD | AVR-PikF | Pikp-HMA-SNK-EKE                                    | resolved range as deposited                         | 2.20 Å | Maidment 2023 *eLife* | The only AVR-PikF / Pik-HMA complex in the PDB; uses the engineered SNK-EKE Pikp scaffold. |

No structure of AVR-PikF alone is deposited.

### 2.4 Complex templates available for chain assembly

The pipeline complex `pikp1_avrpikf_complex.pdb` was assembled with chain A
as an AF3 prediction of native Pikp-1 HMA(186–263) and chain C as AVR-PikF
lifted from 7B1I. The user's complex notes record this directly. The
positioning of chain A relative to chain C was anchored on the Pik-HMA /
AVR-Pik binding pose established in the De la Concepcion 2018 structures
(6G10 / 6FU9), which are highly homologous to the OsHIPP19 / AVR-PikF
binding mode in 7B1I (Maidment 2021 *JBC* explicitly compares the two).
**A usable complex template therefore exists**; this is *not* a
campaign-blocking constraint.

## 3. Structural-input decision

The fallback ladder lands on **rung 5: an AlphaFold3 monomer prediction of
the native Pikp-1 HMA, residues 186–263**, with chain C built from rung 1
data (the AVR-PikF crystal structure 7B1I). This decision is carried over
from the user's verified prior knowledge and is justified by the
literature inventory above:

- **Rung 1 (native receptor / target-effector complex crystal)** — does not
  exist for Pikp-1 / AVR-PikF. AVR-PikF has only ever been crystallised
  with OsHIPP19-HMA (7B1I) and with the engineered Pikp-SNK-EKE scaffold
  (7QZD).
- **Rung 2 (native receptor / related-effector complex crystal)** — would
  be 6G10 (Pikp-HMA / AVR-PikD) or 6G11 (Pikp-HMA / AVR-PikE). Both have
  the right binding pose but neither resolves the full Interface 3 in the
  ATOM records, per the user's inspection. Disqualifying for a campaign
  that explicitly targets Interface 3.
- **Rung 3 (native receptor alone + target-effector alone)** — there is no
  Pikp-1 HMA apo structure that resolves the full C-terminus, and no
  AVR-PikF apo structure exists.
- **Rung 4 (engineered receptor crystal)** — would be 7QZD (Pikp-SNK-EKE /
  AVR-PikF). Two disqualifications: (i) the SNK-EKE mutations Ser258Glu /
  Asn261Lys / Lys262Glu sit *inside* the campaign's primary redesign target
  (Interface 3, positions 73, 76, 77), so anchoring on this structure
  would lock in the very residues the campaign is trying to redesign;
  (ii) per the user's verified inspection, 7QZD's ATOM records also leave
  several C-terminal residues unresolved.
- **Rung 5 (AF3 monomer prediction of native receptor)** — usable, and is
  the chosen rung. Validation against the high-resolution Pikm-HMA /
  AVR-PikD crystal (6FU9 at 1.20 Å) and the closely-related Pikp-HMA /
  AVR-PikD crystal (6G10 at 1.35 Å) is straightforward — most of the
  domain has crystallographic ground truth, and AF3 only needs to do
  competent modelling on the C-terminal tail residues (~258–263) that
  are missing from the deposits. The user has already verified this
  prediction is reliable for this domain.

**AF3 prediction specification.**

- Input sequence: residues 186–263 of UniProt **E9KPB5** (Pikp-1, *Oryza
  sativa subsp. japonica*).
- Boundary justification: residue 186 begins the construct used by the
  Banfield-lab Pik-HMA crystallography line (Maqbool 2015, all subsequent
  De la Concepcion / Maidment papers). Residue 263 is the natural
  C-terminus of the HMA fold and the last residue Lys262 plus the
  flanking Asp263 are the load-bearing Interface 3 contacts (De la
  Concepcion 2018).
- AF3 numbering: output is renumbered 1–78 with position 1 = Pikp-1
  residue 186 and position 78 = Pikp-1 residue 263. The pipeline contig
  is written in this 1–78 numbering.
- Validation step on AF3 model: superpose onto chain A of 6G10 over
  residues 186–258 (positions 1–73) and over Pikm-HMA chain in 6FU9
  for full-length comparison. RMSD on Cα over the resolved range should
  be <1.0 Å for the model to be considered fit for this campaign;
  C-terminal residues 258–263 (positions 73–78) cannot be backbone-
  validated against any deposit and must be accepted on the strength of
  AF3's pLDDT in that region.

## 4. Binding interfaces

The Pik-HMA / AVR-Pik binding interface decomposes into three sub-interfaces
named **Interface 1**, **Interface 2**, and **Interface 3** in the
nomenclature of De la Concepcion et al. 2018. Both Interface 2 and
Interface 3 carry effector polymorphisms; Interface 1 does not.
Receptor-side residues are given in **Pikp-1 full-length numbering**
(186–263) with the **pipeline (1–78) position** in parentheses.
Where a residue is described in the literature using Pikm-HMA numbering
(De la Concepcion 2018), I give the structurally equivalent Pikp-HMA
residue (which may be a different amino acid identity, since Pikp and
Pikm differ at multiple positions across all three interfaces). Effector
residues are given in AVR-Pik **precursor (1–113) numbering** — i.e. the
numbering chain C inherits from PDB 7B1I.

### 4.1 Interface 1 — N-terminus / β1 area (Pikp-1 ~188–195, positions ~3–10)

- **Structural description.** The N-terminal turn and β1 of the HMA fold
  pack against the C-terminal tail of the effector's MAX fold. Per De la
  Concepcion 2018, this interface contributes a single weak side-chain
  H-bond and one hydrophobic side-chain contact in the Pikm-HMA / AVR-PikD
  structure (6FU9); the rest of the contacts are main-chain.
- **Pikp/Pikm sequence divergence at Interface 1.** Pikm-HMA Glu188,
  Met189 and Lys191 (the residues called out by name in De la Concepcion
  2018 Fig. 3a) correspond to Pikp-1 **Lys188, Gln189, Ile191** — none
  of the Pikm Interface-1 side chains is conserved in Pikp-1. Lys195
  (Pikm) corresponds to Pikp-1 **Val195**. The interface in Pikp-1 is
  therefore largely main-chain–mediated; the side-chain contacts are
  Pikm-specific. This is consistent with both alleles binding AVR-PikD
  in vitro (Maqbool 2015, De la Concepcion 2018) — Interface 1 is too
  divergent in side-chain identity between alleles to be the
  specificity-determining surface.
- **Effector-side residues.** Ile49, Thr69, Asp66 (Pikm-HMA contacts);
  none are polymorphic across the AVR-Pik allelic series (Longya et al.
  2019, Kanzaki et al. 2012).
- **Specificity role.** Conserved across all AVR-Pik variants. Provides
  scaffold contacts; not specificity-determining.
- **Design implication.** **Leave fixed.** The β1 backbone and the
  surrounding scaffold should be preserved native — disrupting the
  geometry here would compromise the HMA fold itself.

### 4.2 Interface 2 — β2/β3 region (Pikp-1 218–232, positions 33–47)

- **Structural description.** β2 (Pikp-1 ~218–227, positions 33–42),
  the β2–β3 loop, and the start of β3 (Pikp-1 ~228–234, positions 43–49).
  The most extensive of the three sub-interfaces by buried area
  (De la Concepcion 2018, supplementary interface analysis).
- **Receptor-side residues (Pikp / pipeline / role).**
  - **Ser218** (33, H-bond to effector His46 in AVR-PikD complexes —
    Maqbool 2015; this is the residue contacting the polymorphic position
    46 of the effector).
  - **Asp224** (39, β2/β3 loop, equivalent to Pikm Asp225 which
    salt-bridges effector Arg64 — De la Concepcion 2018).
  - **Glu230** (45, H-bonds effector His46 in Pikp-HMA / AVR-PikD —
    Maqbool 2015; **the canonical Pikp/Pikm polymorphism site** — Pikm
    has Val231 here, lacking H-bond capacity to His46).
  - **Val232** (47, hydrophobic contact with effector His46 — Maqbool 2015).
- **Effector-side residues.** **Position 46** (His in AVR-PikD/E; **Asn in
  AVR-PikA/F/C** — Longya 2019); position 47 (Pro in D/E/C; **Ala in
  A/F** — Longya 2019); position 48 (Gly in D/E/C; **Asp in A/F** —
  Longya 2019). Position 67 (Ala in D/E/A/F; **Asp in C only** — defines
  AVR-PikC by steric clash, De la Concepcion 2018). Plus the conserved
  scaffold contacts Arg64 and Asp66.
- **Specificity role.** **Specificity-determining for the Pikp/Pikm
  allelic split and for the AVR-PikD/E/A series.** The Pikp Glu230 →
  Pikm Val231 polymorphism shifts which AVR-Pik variants give a
  productive H-bond network; Pikm-HMA's bulk binding affinity to
  AVR-PikE and AVR-PikA depends on additional contacts outside this
  sub-interface (Interface 3 — see §4.3) compensating for the loss of
  the Glu230 H-bond.
- **Specificity role for AVR-PikF specifically.** AVR-PikF inherits the
  AVR-PikA-like Asn46 / Ala47 / Asp48 / Ala67 set (Longya 2019). Pikp-HMA
  *cannot* productively bind AVR-PikA or AVR-PikF at Interface 2 — the
  H46→N46 substitution alone is sufficient to abolish the canonical
  Ser218 / Glu230 H-bond network (De la Concepcion 2018 Fig. 5a–c).
- **Design implication.** **Primary redesign target.** Pikp-HMA Interface 2
  is intrinsically incompatible with AVR-PikF's polymorphic loop; the
  contig must allow RFDiffusion to rebuild this region to find an
  alternative residue set that accommodates Asn46 / Ala47 / Asp48.
  Pikm-HMA already provides one natural solution (the Glu230→Val231
  substitution plus Interface 3 compensation), and SNK-EKE provides
  another via Interface-3 changes alone (Maidment 2023) — the campaign
  aims for de novo alternatives to both.

### 4.3 Interface 3 — β4 to C-terminal tail (Pikp-1 254–263, positions 69–78)

- **Structural description.** β4 (Pikp-1 ~252–258, positions 67–73) plus
  the C-terminal extension residues 259–263 (positions 74–78). β4 forms
  main-chain hydrogen bonds with β3 of the effector MAX fold, and the
  C-terminal tail places the receptor's Lys262 side-chain into a pocket
  on the effector lined by Glu53, Tyr71, Ser72 and Trp74 (De la
  Concepcion 2018, Fig. 4).
- **Receptor-side residues (Pikp / pipeline / role).**
  - β4 main-chain (Met254–Ser258, positions 69–73): main-chain H-bond
    pairing with effector β3.
  - **Ala260** (75, packing residue; Pikp-specific conformation that has
    to "loop out" to maintain Lys262 in the effector pocket — De la
    Concepcion 2018, Fig. 4b).
  - **Asn261** (76, Pikp-specific; loops out to maintain the C-terminal
    geometry — De la Concepcion 2018. **Pikh-HMA differs from Pikp-HMA
    by Asn261Lys at this single position**, which is sufficient to
    expand Pikh's recognition profile to AVR-PikC binding in vitro —
    Maidment 2021 *PLoS Pathog*.)
  - **Lys262** (77, the central Interface 3 contact; salt bridge / H-bonds
    to effector Glu53 and Ser72 — De la Concepcion 2018, Maqbool 2015).
  - **Asp263** (78, C-terminal anchor of the construct).
- **Effector-side residues.** Glu53, Tyr71, Ser72, Trp74 (the "Lys262
  pocket"; conserved across AVR-PikD/E/A/B). **Position 78** is **Met in
  AVR-PikA/D/E/C and Lys in AVR-PikF** (Longya 2019) — this is the
  defining AVR-PikF substitution. Met78 packs with the β4 / Lys262 side
  of the receptor in cognate complexes; Lys78 introduces both bulk and a
  positive charge directly opposite the receptor's Lys262, electrostatically
  clashing and (per Longya 2019, gel filtration) preventing complex
  formation entirely with both Pikp-HMA and Pikm-HMA in vitro.
- **Specificity role.** **Specificity-determining for AVR-PikF
  recognition.** This is the *unique* lesion that distinguishes AVR-PikF
  from AVR-PikA. AVR-PikA differs from AVR-PikF only at this position
  (Met78 vs Lys78) and is recognised by Pikm; AVR-PikF is not.
- **Design implication.** **Primary redesign target.** A successful
  AVR-PikF binder needs a residue set at positions 69–78 that
  accommodates Lys78. The published SNK-EKE solution achieves this with
  S258E + N261K + K262E (positions 73, 76, 77), which the campaign
  explicitly does *not* anchor on. The C-terminal block of the contig
  must therefore be flexible — including Lys262 — so RFDiffusion can
  search for alternatives.

### 4.4 The α2 scaffold helix is not a redesign target

- **Structural description.** α2 of the HMA fold spans approximately
  Pikp-1 residues 245–251 (positions 60–66), nestled between β3 and β4 on
  the *opposite* face of the domain from the AVR-Pik binding surface.
- **Specificity role.** None of the AVR-Pik / Pik-HMA contacts described
  by Maqbool 2015 or De la Concepcion 2018 involve α2 side chains. Its
  role is structural — positioning β4 (and therefore Lys262) correctly
  relative to the effector pocket.
- **Design implication.** **Leave fixed.** Disrupting α2 would risk
  destabilising the HMA fold; redesigning α2 buys nothing for
  binding-surface specificity.

### 4.5 Why AVR-PikF evades Pikp-1 — synthesis

Pikp-HMA fails to recognise AVR-PikF for **two compounded reasons**, one
at each polymorphic interface:

1. **Interface 2.** AVR-PikF has Asn46 (not His46), so the Pikp-specific
   Ser218 / Glu230 H-bond network to His46 cannot form (De la Concepcion
   2018, Maqbool 2015). This alone is enough to abolish recognition —
   AVR-PikA, which shares Asn46 with AVR-PikF, is also not recognised by
   Pikp (Kanzaki 2012, Longya 2019).
2. **Interface 3.** AVR-PikF additionally has Lys78 (not Met78), placing
   a positive charge directly across from Pikp's Lys262. This disrupts
   the Lys262 pocket interactions and is the unique AVR-PikF lesion that
   prevents binding even by Pikm-HMA, which can otherwise compensate for
   Interface 2 disruption (Longya 2019: Pikm-HMA does not bind AVR-PikF
   in gel-filtration assays).

Both interfaces must therefore be redesigned for the campaign to recover
AVR-PikF recognition.

## 5. Numbering translation

| Convention                                      | What it means                                                                                                              |
|--------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| **Receptor literature numbering**                | Pikp-1 full-length (1–1142) numbering, used by every Pik-HMA paper since Maqbool 2015. Pikp-1 HMA spans residues 186–263. |
| **Receptor pipeline numbering (chain A)**        | 1–78, with position 1 = Pikp-1 residue 186 and position 78 = Pikp-1 residue 263. **Offset: pipeline = literature − 185.**  |
| **Effector literature numbering**                | AVR-Pik precursor (1–113), with signal peptide 1–21 and mature protein 22–113. Polymorphic positions 46, 47, 48, 67, 78 are in this numbering. |
| **Effector pipeline numbering (chain C)**        | Inherits 7B1I chain C numbering, which is AVR-Pik precursor numbering — so **identical to literature**. **Resolved range in 7B1I is 33–113 (verified by inspection of `experiments/inputs/structures/7B1I.pdb`); the assembled `pikp1_avrpikf_complex.pdb` carries chain C with residues 32–113.**|

### 5.1 Receptor key-residue translation table

| Pikp-1 residue (full) | Identity (Pikp) | Identity (Pikm equivalent) | Pipeline position | Interface | Specificity-determining? |
|-----------------------|------------------|------------------------------|--------------------|-----------|---------------------------|
| 188                   | Lys              | Glu                          | 3                 | 1         | No (Pikm side-chain only) |
| 189                   | Gln              | Met                          | 4                 | 1         | No (Pikm side-chain only) |
| 191                   | Ile              | Lys                          | 6                 | 1         | No (Pikm side-chain only) |
| 195                   | Val              | Lys                          | 10                | 1/2       | No (Pikm side-chain only) |
| 218                   | Ser              | Ser (conserved)             | 33                | 2         | Yes — H-bond to effector His46 |
| 224                   | Asp              | Asp (conserved at 225)      | 39                | 2         | Yes — salt bridge to effector Arg64 |
| 230                   | **Glu**          | **Val (at 231)**            | 45                | 2         | **Yes — canonical Pikp/Pikm polymorphism** |
| 232                   | Val              | Val (conserved)             | 47                | 2         | Yes — hydrophobic contact to His46 |
| 254                   | Met              | Met (conserved)             | 69                | 3         | No (β4 main-chain pairing, conserved) |
| 258                   | **Ser**          | Gln (at 259)                | 73                | 3         | **Yes — SNK-EKE engineering site** |
| 260                   | Ala              | Val (at 261)                | 75                | 3         | Yes — Pikp-specific packing |
| 261                   | **Asn**          | Lys (at 262)                | 76                | 3         | **Yes — NK-KE/SNK-EKE engineering site; Pikh polymorphism** |
| 262                   | **Lys**          | (no equivalent; Pikm Lys is at 262 directly) | 77 | 3 | **Yes — central Interface 3 contact (in *both* alleles)** |
| 263                   | Asp              | (Pikm extends to 264)        | 78                | 3         | Yes — C-terminal anchor |

### 5.2 Effector key-residue translation table

Chain C numbering matches AVR-Pik precursor numbering directly, so no
offset is needed.

| Residue (precursor) | AVR-PikF identity | AVR-PikD / AVR-PikA reference | Interface | Notes |
|---------------------|-------------------|-------------------------------|-----------|-------|
| 46                  | Asn (N)           | His (D, E) / Asn (A, F)       | 2         | F shares with A; Pikp-HMA cannot form the canonical His46 H-bond network with Asn46. |
| 47                  | Ala (A)           | Pro (D, E, C) / Ala (A, F)    | 2         | F shares with A. |
| 48                  | Asp (D)           | Gly (D, E, C) / Asp (A, F)    | 2         | F shares with A. |
| 67                  | Ala (A)           | Ala (A, D, E, F) / Asp (C only) | 2 (steric, AVR-PikC-only) | F is wild-type at this position. |
| 78                  | **Lys (K)**       | Met (D, E, A, C)              | 3         | **AVR-PikF defining substitution.** Introduces a positive charge at Pikp-HMA Lys262 contact. |

## 6. Contig string design

All contigs are written in pipeline numbering (chain A 1–78). Chain B in
the contig syntax corresponds to the effector chain in the pipeline; in
the assembled complex this is chain C, but RFDiffusion's contig syntax
reads it as the second chain regardless of label.

### 6.1 Design principle from §4

The literature analysis pins down four regions on Pikp-1 HMA in pipeline
numbering:

- **N-terminal scaffold (positions 1–32)** — covers Interface 1 (positions
  3–10) and α1, neither of which is specificity-determining for AVR-Pik
  recognition. *Anchor.*
- **Interface 2 (positions 33–49)** — covers β2, the β2/β3 loop, and the
  β3 N-terminus. Contains the canonical Pikp-1 Interface-2 residues
  Ser218(33), Asp224(39), Glu230(45) and Val232(47) that contact effector
  position 46. *Flex* — this is the region that must be redesigned to
  accommodate AVR-PikF's Asn46 / Ala47 / Asp48 polymorphisms.
- **β3-tail / α2 scaffold (positions 50–68)** — covers the rest of β3,
  α2, and the residue immediately before β4. None of these contact the
  effector; they are the structural cradle that positions Interface 3.
  *Anchor.*
- **Interface 3 (positions 69–78)** — covers β4 and the C-terminal tail.
  Contains Met254(69), Ser258(73), Ala260(75), Asn261(76), Lys262(77),
  and Asp263(78) — i.e. the entire SNK-EKE engineering window plus the
  flanking residues. *Flex* — this is the region that must be redesigned
  to accommodate AVR-PikF's Lys78 polymorphism, with the user's prior
  knowledge mandating that the C-terminal residues sit inside a flexible
  block (not anchored to the SNK-EKE engineered solution).

### 6.2 Primary contig

```
A1-32/12-20/A50-68/8-14 C
```

**Pikp-1 full-length (literature) cross-reference**: `A186-217/12-20/A235-253/8-14 C` —
i.e. anchor 186–217, flex 12–20 (over native residues 218–234), anchor 235–253, flex 8–14
(over native residues 254–263). The flex-block length ranges (12–20 and 8–14) are
unchanged between numberings. Use the pipeline-numbered form for RFDiffusion; use the
literature form when checking positions against figures in 6G10/6FU9/7B1I or
Banfield-lab papers.

| Segment      | Pipeline positions | Pikp-1 residues  | Native length | Sampled length | Structural element                                | Role                                                                                                        |
|--------------|--------------------|--------------------|---------------|-----------------|---------------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| `A1-32`      | 1–32              | 186–217          | 32            | 32 (anchor)    | N-terminus, β1, α1, end of α1                    | **Anchored.** Holds Interface 1 and the entire α1 scaffold helix native.                                    |
| `12-20`      | 33–49 (flex)      | 218–234          | 17            | 12–20          | β2, β2/β3 loop, β3 N-terminus                    | **Flexible.** Encompasses Ser218, Asp224, Glu230, Val232 — all canonical Interface 2 contacts.              |
| `A50-68`     | 50–68             | 235–253          | 19            | 19 (anchor)    | β3 tail, α2, β3-α2-β4 transition                 | **Anchored.** Holds the α2 scaffold native.                                                                 |
| `8-14`       | 69–78 (flex)      | 254–263          | 10            | 8–14           | β4, C-terminal tail                               | **Flexible.** Encompasses Met254, Ser258, Ala260, Asn261, **Lys262**, Asp263 — all of Interface 3.          |
| `C`          | full chain C      | AVR-PikF         | n/a           | n/a            | full effector                                     | **Fully fixed** — RFDiffusion does not modify the target.                                                   |

**Why this contig follows from §4.** The flex blocks are exactly the two
specificity-determining regions identified in §4.2 and §4.3; the anchor
blocks cover the N-terminal Interface 1 (non-specificity-determining,
§4.1) and the α2 scaffold (structural, §4.4). The C-terminal tail
including Lys262 sits inside the second flex block, satisfying the
user's prior-knowledge requirement.

**Boundary check (mid-secondary-structure caveat).** The user must verify
in ChimeraX before committing. Two boundaries warrant attention:

- Position 32 / 33 boundary (anchor → flex). Per De la Concepcion 2018
  Fig. 1c, β2 of the HMA fold begins at Pikp-1 residue ~217–218 (position
  32–33). The boundary at 32/33 may therefore fall at the very start of
  β2. In ChimeraX, confirm position 32 is the end of the β1-α1-β2 loop
  (or the first residue of β2 with no upstream H-bond partners).
- Position 68 / 69 boundary (anchor → flex). Per De la Concepcion 2018
  the start of β4 is Met254 (position 69). The boundary at 68/69 should
  therefore be one residue *before* β4 begins, in the loop linking α2
  and β4. Confirm in ChimeraX.

### 6.3 Alternative contig — conservative, Interface 3 only

```
A1-67/8-14 C
```

**Literature cross-reference**: `A186-252/8-14 C`.

| Segment   | Pipeline positions | Pikp-1 residues | Sampled length | Role                                                                                                       |
|-----------|--------------------|------------------|----------------|------------------------------------------------------------------------------------------------------------|
| `A1-67`   | 1–67              | 186–252         | anchor          | **Anchored.** Holds Interface 1, α1, all of Interface 2 native, β3 tail and α2 scaffold native.            |
| `8-14`    | 68–78 (flex)      | 253–263         | 8–14           | **Flexible.** Interface 3 only.                                                                            |

Redesigns *only* Interface 3. Leaves Interface 2 native, which means the
designs will retain Pikp's intrinsic incompatibility with the AVR-PikF
Asn46 / Ala47 / Asp48 set — so binding will only be recovered if
Interface 3 redesign alone is sufficient to make up the Interface 2
deficit. On the De la Concepcion 2018 / Longya 2019 evidence this is
unlikely to suffice (Pikp-HMA loses AVR-PikA recognition because of
Interface 2 alone, before Interface 3 enters the picture), but the
SNK-EKE result (Maidment 2023) is a published proof-of-principle that
Interface-3-only changes can in principle expand recognition to AVR-PikF.
Included as a deliberately under-powered comparator: if Interface 3
redesign alone is sufficient, that would be a strong and surprising
result; if not, nothing has been lost.

### 6.4 Alternative contig — aggressive, also redesign Interface 2's β3 portion

```
A1-32/12-20/A50-67/9-15 C
```

**Literature cross-reference**: `A186-217/12-20/A235-252/9-15 C`.

| Segment    | Pipeline positions | Pikp-1 residues | Sampled length | Role                                                                                                          |
|------------|--------------------|------------------|----------------|---------------------------------------------------------------------------------------------------------------|
| `A1-32`    | 1–32              | 186–217         | anchor          | **Anchored.** As primary.                                                                                     |
| `12-20`    | 33–49 (flex)      | 218–234         | 12–20          | **Flexible.** Interface 2.                                                                                    |
| `A50-67`   | 50–67             | 235–252         | anchor          | **Anchored.** α2 scaffold only — boundary moved one residue earlier than primary, releasing position 68.       |
| `9-15`     | 68–78 (flex)      | 253–263         | 9–15           | **Flexible.** Interface 3 plus the residue immediately before β4 (β3-α2-β4 transition).                       |

Same as the primary contig but extends the Interface 3 flex block by one
residue at its N-terminal end (position 68 is included). This gives
RFDiffusion freedom to adjust the geometry of the β3-α2-β4 transition
that positions β4 against the effector. Use only if the primary contig
underperforms on β4 placement.

### 6.5 Alternative contig — tighter flex ranges (more conservative search)

```
A1-32/17-17/A50-68/10-10 C
```

**Literature cross-reference**: `A186-217/17-17/A235-253/10-10 C`.

Same anchor blocks as the primary contig but each flex block is fixed at
the native segment length (17 and 10 residues respectively). Removes
RFDiffusion's freedom to insert or delete residues in the redesigned
regions. Useful as a comparator to the primary contig if the campaign
wants to disentangle the contribution of length variation from sequence
variation.

### 6.6 Alternative contig — generous flex ranges (broader search)

```
A1-32/14-25/A50-68/8-18 C
```

**Literature cross-reference**: `A186-217/14-25/A235-253/8-18 C`.

Same anchor blocks as the primary contig but each flex block is given a
wider length range (12–20 → 14–25 for Interface 2; 8–14 → 8–18 for
Interface 3). Useful if the primary contig converges on too narrow a
design distribution.

## References

- Białas A, Langner T, Harant A, Contreras MP, Stevenson CEM, Banfield MJ,
  Kourelis J, Kamoun S, Wu C-H, Krasileva KV. *Two NLR immune receptors
  acquired high-affinity binding to a fungal effector through convergent
  evolution of their integrated domain*. eLife 10, e66961 (2021).
  doi:10.7554/eLife.66961.
- Bentham AR, De la Concepcion JC, Benjumea JV, Kourelis J, Jones S,
  Mendel M, Stubbs J, Stevenson CEM, Maidment JHR, Youles M, Erickson FL,
  Sornay C, Adamiak G, Brunkard JO, Banfield MJ, Kamoun S. *Allelic
  compatibility in plant immune receptors facilitates engineering of new
  effector recognition specificities*. The Plant Cell 35, 3809–3827
  (2023). doi:10.1093/plcell/koad204.
- De la Concepcion JC, Franceschetti M, Maqbool A, Saitoh H, Terauchi R,
  Kamoun S, Banfield MJ. *Polymorphic residues in rice NLRs expand binding
  and response to effectors of the blast pathogen*. Nature Plants 4,
  576–585 (2018). doi:10.1038/s41477-018-0194-x. **[from
  experiments/inputs/context/]**
- De la Concepcion JC, Franceschetti M, MacLean D, Terauchi R, Kamoun S,
  Banfield MJ. *Protein engineering expands the effector recognition
  profile of a rice NLR immune receptor*. eLife 8, e47713 (2019).
  doi:10.7554/eLife.47713.
- Kanzaki H, Yoshida K, Saitoh H, Fujisaki K, Hirabuchi A, Alaux L,
  Fournier E, Tharreau D, Terauchi R. *Arms race co-evolution of
  Magnaporthe oryzae AVR-Pik and rice Pik genes driven by their physical
  interactions*. The Plant Journal 72, 894–907 (2012).
  doi:10.1111/j.1365-313X.2012.05110.x. **[from
  experiments/inputs/context/]**
- Longya A, Chaipanya C, Franceschetti M, Maidment JHR, Banfield MJ,
  Jantasuriyarat C. *Gene Duplication and Mutation in the Emergence of a
  Novel Aggressive Allele of the AVR-Pik Effector in the Rice Blast
  Fungus*. Molecular Plant–Microbe Interactions 32, 740–749 (2019).
  doi:10.1094/MPMI-09-18-0245-R. **[from
  experiments/inputs/context/]**
- Maidment JHR, Franceschetti M, Maqbool A, Saitoh H, Jantasuriyarat C,
  Kamoun S, Terauchi R, Banfield MJ. *Multiple variants of the fungal
  effector AVR-Pik bind the HMA domain of the rice protein OsHIPP19,
  providing a foundation to engineer plant defense*. Journal of Biological
  Chemistry 296, 100371 (2021). doi:10.1016/j.jbc.2021.100371.
- Maidment JHR, Shimizu M, Bentham AR, Vera S, Franceschetti M, Longya A,
  Stevenson CEM, De la Concepcion JC, Białas A, Kamoun S, Terauchi R,
  Banfield MJ. *Effector target-guided engineering of an integrated
  domain expands the disease resistance profile of a rice NLR immune
  receptor*. eLife 12, e81123 (2023). doi:10.7554/eLife.81123.
- Maidment JHR, Vera S, Franceschetti M, Bentham AR, De la Concepcion JC,
  Banfield MJ. *The allelic rice immune receptor Pikh confers extended
  resistance to strains of the blast fungus through a single polymorphism
  in the effector binding interface*. PLoS Pathogens 17, e1009368 (2021).
  doi:10.1371/journal.ppat.1009368.
- Maqbool A, Saitoh H, Franceschetti M, Stevenson CEM, Uemura A, Kanzaki
  H, Kamoun S, Terauchi R, Banfield MJ. *Structural basis of pathogen
  recognition by an integrated HMA domain in a plant NLR immune
  receptor*. eLife 4, e08709 (2015). doi:10.7554/eLife.08709.
- Varden FA, De la Concepcion JC, Maidment JHR, Franceschetti M, Banfield
  MJ. *Cross-reactivity of a rice NLR immune receptor to distinct
  effectors from the rice blast pathogen Magnaporthe oryzae provides
  partial disease resistance*. Journal of Biological Chemistry 294,
  13006–13016 (2019). doi:10.1074/jbc.RA119.007730.
- Yoshida K, Saitoh H, Fujisawa S, Kanzaki H, Matsumura H, Yoshida K,
  Tosa Y, Chuma I, Takano Y, Win J, Kamoun S, Terauchi R. *Association
  genetics reveals three novel avirulence genes from the rice blast
  fungal pathogen Magnaporthe oryzae*. The Plant Cell 21, 1573–1591
  (2009). doi:10.1105/tpc.109.066324.
- Zdrzałek R, Stevenson CEM, Pham Z, Schimmel L, Kourelis J, Kamoun S,
  Banfield MJ. *Bioengineering a plant NLR immune receptor with a robust
  binding interface toward a conserved fungal pathogen effector*. PNAS
  121, e2402872121 (2024). doi:10.1073/pnas.2402872121.
- UniProt entries: **E9KPB5** (Pikp-1, *Oryza sativa subsp. japonica*),
  **C4B8C0** (canonical AVR-Pik, *Magnaporthe oryzae*).
