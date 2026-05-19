# Pipeline notes 16 — Pose-solver diagnostics, fixes, and PikF factorial campaign

Date: 2026-05-19
Branch: phase4-impl
Predecessor: pipeline_notes15.md (cold-start RMSD analysis, proposed Pia sidechain/misfold analyses — those analyses are deferred)

---

## 0. TL;DR

- `pose_solved_5a_2` (pikp1_avrpia, hotspot removed, sc_threshold 0.3) **gave no improvement** in hitrate over 5a or 3a. Hotspot and SC filter stringency are not the Pia bottleneck.
- Investigation pivoted from "RFDiffusion params" to "input pose geometry".
- Diagnosed the pose_solver as **over-packing the interface**: its objective is CA-only (no sidechain awareness), and `--max-pair-distance 4.0` is unrealistically tight relative to natural interface CA-CAs (which cluster at 4.5–7 Å).
- Fixed: re-solved PikF at `--max-pair-distance 6.0`, replacing the user's third pair (A35-B44, 6.58 Å on the AF3 ground truth) with A38-B42. Resulting pose's gap-index = 3.45 Å (vs AF3 3.41 — within 1 %). Stored at `tests/haddock/pose_solver/Pikp-1_HMA_PikF_v2/`.
- Added pair-residue sidechain-clash penalty (CLI: `--pair-sc-clash-cutoff`); fixed a pre-existing 2-pair crash bug. The pair-only sidechain penalty is necessary but not sufficient to prevent global over-packing.
- Built **2 × 3 factorial campaign for PikF**: {CRYSTAL, POSED v2} × {full_test_run contig, auto-5 Å contig, auto-3 Å contig} = 6 runs. All cloned from `pose_solved_5a/params.yml` with `hotspot=""`, `sc_threshold=0.3`. Ready for submission on the HPC.

---

## 1. 5a_2 outcome — RFDiffusion knobs are not the bottleneck

`pose_solved_5a_2` (pikp1_avrpia, hotspot="" instead of B22,B24,B31,B33; sc_threshold=0.3 instead of 0.5; otherwise identical to 5a) produced **no improvement in hitrate** over `pose_solved_5a` or `pose_solved_3a`.

This rules out:
- hotspot mis-specification
- over-strict SC filter

… as the dominant failure mode for Pia. The cold-start receptor RMSD failures noted in `pipeline_notes15.md` (83 % in 3a, 51 % in 5a) cannot be remediated at the RFDiffusion / Rosetta layer. The bottleneck has to be **upstream**: either the input pose, the design region, or the MPNN/Boltz refolding behaviour.

Decision: pivot to investigating input pose geometry vs. the natural complex.

---

## 2. New diagnostic tool: `interface_metrics.py`

Added `tests/haddock/pose_solver/interface_metrics.py` for systematic comparison of two-chain complex interfaces. Metrics:

- CA-CA contacts ≤ 8 Å (count + distance distribution)
- Heavy-atom contacts ≤ 5 Å (count + histogram by 0.5-Å bins)
- Clash counts at 2, 2.5, 3, 3.5 Å
- Interface residue count per chain
- Composition fractions (hydrophobic / polar / charged / aromatic) per side
- BSA via Shrake-Rupley (biopython SASA)
- **Gap-index**: mean nearest-other-chain heavy-atom distance over interface atoms — captures global packing density
- H-bond candidates (donor-acceptor heavy-atom pairs ≤ 3.5 Å, geometric only)

CLI:
```
python interface_metrics.py PDB1 [PDB2 ...] --chains A:B [A:C ...] [--json-out FILE]
```

Outputs at `tests/haddock/pose_solver/interface_metrics_out/`.

---

## 3. Comparing pose_solver outputs to the natural complex

### 3.1 Benchmark choice

The "natural" reference is `tests/full_test_run/af3_pikp1_native_avrpikf_complex.pdb` — the AF3-predicted PikP1+AvrPikF complex used as input to the original successful full pipeline test run (`PikP1_AvrPikF_full_test`). User confirmation:

> The PikF complex used as input was just created by aligning with PikD from PDB 6G10 (and also the AF3 Pikp-1 with the 6G10 Pikp-1 because it has a few extra residues). The pose exactly matches a solved structure.

Confirmed identical (atom-coordinate diff: one trivial `-0.000`/`0.000` sign) to `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_pik_interface_complex.pdb` (effector chain renamed C→B). All comparisons use this PDB as ground truth.

### 3.2 First-pass comparison (3-way)

`solved_pose_posed.pdb` (pose_solved Pia, v1 4-pair max 4 Å) vs AF3 PikF:

| Metric | AF3 PikF | pose_solved Pia v1 |
|---|---:|---:|
| BSA total | 2717 Å² | 1420 Å² |
| Gap-index mean | 3.41 Å | 3.12 Å (−0.29 Å) |
| Contacts ≤ 5 Å | 857 | 1137 (+33 %) |
| Aromatic on chain A | 4 % | 0 % |

Initial reading: BSA much smaller, packing tighter, aromatics absent.

**User correction:** "Pia lacks an important loop region that forms a large part of the binding interface with PikF, so it's natural that the BSA might be smaller." The BSA gap is partly explained by Pia geometry; gap-index and aromatic absence are still informative signals.

### 3.3 Apples-to-apples PikF solve

To remove the Pia-vs-PikF confound, ran the pose_solver on PikF directly (same receptor, same effector geometry as the AF3 input). Pair set provided by user: `A73-B76 A71-B78 A35-B44`. Output at `tests/haddock/pose_solver/Pikp-1_HMA_PikF/solved_pose_posed.pdb` (v1).

Metrics vs AF3 PikF:

| Metric | AF3 PikF (truth) | pose_solved PikF v1 (4 Å, A35-B44) |
|---|---:|---:|
| Pair distances achieved | 5.02 / 4.69 / 6.58 (natural) | 3.98 / 3.73 / 3.96 |
| CA contacts ≤ 8 Å | 52 | 89 (+71 %) |
| Heavy contacts ≤ 5 Å | 857 | 1327 (+55 %) |
| Gap-index | 3.41 | **3.18 (−0.23)** |
| Contacts per iface heavy atom | 3.06 | **3.75 (+22 %)** |
| BSA total | 2717 | 2373 (−13 %) |
| Aromatic on chain A | 4 % | 0 % |

The over-packing was real and consistent across Pia and PikF — not a Pia-specific artefact.

---

## 4. Root cause: solver is CA-only

Direct read of `tests/haddock/pose_solver/pose_solver.py` `_obj()` (lines 199–264):

| Penalty term | Operates on |
|---|---|
| `validity_loss` (W_VALIDITY = 1e6) | Signed depth of pair-residue **CA** inside CA-derived convex hull |
| `dist_loss` (W_DIST + W_LOWER = 1e4 each) | Pair-residue **CA-CA** distance vs [min, max] |
| `excl_loss` (W_EXCL = 1e4) | Exclusion pair **CA-CA** distance |
| `interp_loss` (W_INTERP = 10) | **CA** points inside opposite CA hull |

**No sidechain or heavy-atom term in the objective.** Heavy-atom clashes are computed post-hoc in `analyse_contacts()` for reporting only.

### 4.1 Natural CA-CA range is 4.5–7 Å, not ≤ 4 Å

Quick analysis on the AF3 ground truth (`A73-B76`, `A71-B78`, `A35-B44`):

| Pair | CA-CA on AF3 ground truth | vs `--max-pair-distance 4.0` |
|---|---:|---|
| A73 – B76 | 5.02 Å | 1.02 Å too tight |
| A71 – B78 | 4.69 Å | 0.69 Å too tight |
| A35 – B44 | **6.58 Å** | **2.58 Å too tight (A35 not in top-20 closest pairs)** |

Distribution of A-vs-B CA-CA pairs on the AF3 ground truth:

| Cutoff | Pairs |
|---|---:|
| ≤ 4.0 Å | 4 |
| ≤ 5.0 Å | 8 |
| ≤ 6.0 Å | 17 |
| ≤ 7.0 Å | 33 |
| ≤ 8.0 Å | 52 |

Of the 4 pairs ≤ 4 Å on the AF3 input, three (A77/78–B49/50/51 at 1.83–3.15 Å) are AF3 placement artefacts at the disordered C-terminal region. **Effectively the natural interface has zero CA-CA pairs < 4 Å outside the disordered region.** Forcing pair atoms into a ≤ 4 Å window therefore guarantees an unnatural conformation.

### 4.2 AF3 sub-2 Å clashes are localised to the disordered C-terminus

The AF3 ground truth has 55 atomic clashes < 2 Å — looks like the natural complex is also clashy. Closer look:

- A76 contacts B71: 17 atom-pairs < 2 Å
- A77 contacts B50: 16
- A78 contacts B51/B50/B52: 19
- All other A-residue clashes: 3 atom-pairs (A43-B66 ×2, A3-B49 ×1)

**52 of 55 sub-2 Å AF3 clashes are at the A76-78 disordered C-terminus.** The natural complex has effectively no atomic clashes elsewhere. So a "good" pose should have very few clashes at non-disordered residues — the pose_solver v1 (44 such clashes, none in disorder) is way over.

### 4.3 Closest realistic AF3 pairs for choosing anchors

| Pair | CA-CA |
|---|---:|
| A76 – B50 | 4.52 |
| A71 – B78 | 4.69 ✓ user chose |
| A70 – B79 | 4.74 |
| A78 – B49 | 4.81 |
| A73 – B76 | 5.02 ✓ user chose |
| A76 – B49 | 5.47 |
| A77 – B51 | 5.50 |
| **A38 – B42** | **5.52** |
| A37 – B42 | 5.62 |
| A70 – B78 | 5.63 |

A38-B42 (5.52 Å) is the closest realistic anchor in the N-terminal half of the interface — replacement candidate for the bad A35-B44.

---

## 5. Pose_solver code changes

### 5.1 Bug fix: 2-pair crash

`solve()` unconditionally unpacked `p1, p2, p3 = b_arr[pb[0]], b_arr[pb[1]], b_arr[pb[2]]` before the `if len(pb) >= 3:` check, so 2-pair runs crashed with `IndexError`. Moved unpacking inside the branch.

### 5.2 New: pair sidechain-clash penalty

Added in `pose_solver.py`:

- New constants:
  - `W_PAIR_SC_CLASH = 1e4`
  - `PAIR_SC_CLASH_TOL = 2.0`  (Å)
  - `_BACKBONE_ATOMS = frozenset({'N', 'CA', 'C', 'O', 'OXT'})`
- New function `read_sidechain_heavy(pdb, chain)` — returns `{resnum: np.ndarray}` of sidechain heavy atoms only.
- New CLI arg `--pair-sc-clash-cutoff` (default 2.0; set 0 to disable).
- `_obj()` extended: for each pair, computes min heavy-atom distance between the binder pair-residue sidechain and the target pair-residue sidechain; quadratic penalty if below cutoff.
- `solve()` extended: builds flat sidechain arrays + per-pair offsets, passes through to `_obj()`. GLY pairs are silently skipped (with a printed warning when active).
- Progress printout now reports `Pair sidechain-clash: W=1e+04 cutoff=2.00Å (active pairs: N/N)` when active.

### 5.3 Limitation noted

The pair-only sidechain check protects only the constraint residues — surrounding residues can still over-pack. Demonstrated in v3 (see § 6) where the local check was satisfied but global metrics were worse than v1.

A future change would either:
- extend the penalty to all interface heavy-atom pairs (O(N×M) per iteration; ~10–30 s extra at 1000 restarts), or
- switch pair constraints to heavy-atom min-distance (LogSumExp soft-min) instead of CA-CA.

Not implemented this session.

---

## 6. PikF pose iterations

All under `tests/haddock/pose_solver/`. Same binder (`../data/Pikp-1_HMA.pdb`), target (`../data/PikF.pdb`, chains A + B, derived from `pikp1_avrpikf_pik_interface_complex.pdb` via ChimeraX `save .../PikF.pdb #2/B`), 1000 restarts, `--global-interp`.

| Dir | Pairs | min / max (Å) | sidechain-clash | Outcome |
|---|---|---|---|---|
| `Pikp-1_HMA_PikF/` | A73-B76, A71-B78, A35-B44 | 3.5 / 4.0 | off (CA-only) | v1: over-packed; A35-B44 pulled effector around |
| `Pikp-1_HMA_PikF_v2/` | A73-B76, A71-B78, A38-B42 | 3.5 / 6.0 | off | **Best pose** — gap-index 3.45 Å (vs AF3 3.41); used as POSED input in factorial |
| `PikF_2pair_min4_max6_scclash/` | A73-B76, A71-B78 | 4.0 / 6.0 | on, 2.0 Å | v3 — worse over-packing than v1 (gap-index 3.08, 67 sub-2 Å clashes); confirms 2 pairs underconstrained |

Full metrics:

| Metric | AF3 truth | v1 (3pair, 4 Å) | **v2 (3pair, 6 Å, A38-B42)** | v3 (2pair, 4–6 Å, scclash) |
|---|---:|---:|---:|---:|
| CA contacts ≤ 8 Å | 52 | 89 | **56** ✓ | 61 |
| Heavy contacts ≤ 5 Å | 857 | 1327 | **764** ✓ | 1074 |
| Clashes < 2 Å | 55 (disorder) | 44 | **27** ✓ | 67 |
| Clashes < 2.5 Å | 95 | 107 | **53** ✓ | 127 |
| Clashes < 3.0 Å | 158 | 229 | **101** ✓ | 221 |
| Gap-index (Å) | 3.41 | 3.18 | **3.45** ✓ | 3.08 |
| Contacts / iface atom | 3.06 | 3.75 | **3.07** ✓ | 4.78 |
| BSA (Å²) | 2717 | 2373 | 1793 | 1262 |
| Iface residues A / eff | 23 / 31 | 29 / 37 | 24 / 25 | 14 / 21 |
| H-bond candidates | 31 | 45 | 15 | 30 |

**v2 wins.** Selected as POSED input for the factorial.

Copied to `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_pose_solved_v2.pdb` for HPC use.

---

## 7. New contig deriver

Added `tests/haddock/pose_solver/derive_contigs.py` — re-runs the auto-contig logic on an already-saved pose at any heavy-atom contact cutoff. Replicates the strategy-4a length-range scaling used in `pose_solved_5a/params.yml` (`int(N*0.7 + 0.5)` to `int(N*1.5 + 0.5)`, with N=1 → 1-1). Join-gap = 1 (segments separated by a single fixed residue merge).

Usage:
```
python derive_contigs.py PDB CHAIN_A CHAIN_B CUTOFF [CUTOFF ...]
```

---

## 8. The PikF factorial campaign (6 runs)

### 8.1 Design

| | **CRYSTAL** (`pikp1_avrpikf_pik_interface_complex.pdb`) | **POSED v2** (`pikp1_avrpikf_pose_solved_v2.pdb`) |
|---|---|---|
| **FULL** contig | `crystal_full_test_contig` | `posed_full_test_contig` |
| **5 Å auto** | `crystal_auto5A_contig` | `posed_auto5A_contig` |
| **3 Å auto** | `crystal_auto3A_contig` | `posed_auto3A_contig` |

Auto-derived contigs are per-input (not shared) so each input's design region matches its own clash pattern.

### 8.2 Contig strings (all effector chain B)

| Label | Contig | Residues redesigned (orig) | De novo length range |
|---|---|---|---|
| FULL | `A1-32/10-20/A46-72/6-6 B` | 19 (33-45, 73-78) | 16-26 |
| CRYSTAL 5 Å | `A1-2/2-5/A6-31/6-12/A40-42/5-11/A50-67/8-17 B` | 29 (3-5, 32-39, 43-49, 68-78) | 21-45 |
| CRYSTAL 3 Å | `A1-2/1-1/A4-31/1-3/A34-38/1-1/A40-42/1-1/A44-67/2-5/A71-73/4-8 B` | 13 (3, 32-33, 39, 43, 68-70, 74-78) | 10-19 |
| POSED 5 Å | `A1-2/2-5/A6-31/6-14/A41-42/4-8/A48-68/7-15 B` | 27 (3-5, 32-40, 43-47, 69-78) | 19-42 |
| POSED 3 Å | `A1-2/1-1/A4-32/3-6/A37-38/1-1/A40-42/1-1/A44-70/1-1/A72-73/4-8 B` | 13 (3, 33-36, 39, 43, 71, 74-78) | 11-18 |

Observations worth flagging when interpreting:

- **FULL is the only contig that keeps A3 fixed.** All four auto-derived contigs flag A3 (LYS — long sidechain) as a contact at 3 Å and 5 Å. FULL was hand-designed without it; that may be under-designing in retrospect.
- **3 Å contigs are fragmented** (6 segments, several length-1 stretches). Length-1 segments give RFDiffusion essentially zero geometric flexibility — just a residue substitution.
- **5 Å contigs are large** (27–29 residues redesigned, de novo up to 45). More freedom to fix clashes, but slower and more failure modes (ROG, contact loss).
- **CRYSTAL and POSED diverge most at 5 Å** in the middle of the receptor — CRYSTAL uses A50-67 / A40-42 fixed; POSED uses A48-68 / A41-42 fixed. Reflects the slightly different effector positioning in the two poses.
- **CRYSTAL and POSED converge at 3 Å** — same 13 residues, slightly different placements.

### 8.3 Shared parameter overrides

All six clone `experiments/campaigns/pikp1_avrpia/runs/pose_solved_5a/params.yml` with:

- `pdb_file`, `project_name`, `outdir`, `contigs` — per cell
- `hotspot: ""` (was `B22,B24,B31,B33`)
- `sc_threshold: 0.3` (was 0.5)
- `effector_chain: "B"`

Everything else identical to 5a: 64 designs, 2 MPNN seqs/design, 128 forward to negsteer, `rfdiff_iterations: 100`, `rfdiff_contact_cutoff: 8.0`, etc.

### 8.4 Submission

```
sbatch run_pipeline.slurm.sh ./experiments/campaigns/pikp1_avrpikf/runs/crystal_full_test_contig/params.yml
sbatch run_pipeline.slurm.sh ./experiments/campaigns/pikp1_avrpikf/runs/crystal_auto5A_contig/params.yml
sbatch run_pipeline.slurm.sh ./experiments/campaigns/pikp1_avrpikf/runs/crystal_auto3A_contig/params.yml
sbatch run_pipeline.slurm.sh ./experiments/campaigns/pikp1_avrpikf/runs/posed_full_test_contig/params.yml
sbatch run_pipeline.slurm.sh ./experiments/campaigns/pikp1_avrpikf/runs/posed_auto5A_contig/params.yml
sbatch run_pipeline.slurm.sh ./experiments/campaigns/pikp1_avrpikf/runs/posed_auto3A_contig/params.yml
```

### 8.5 Interpretation guide

| Observation | Implies |
|---|---|
| `posed_full` hitrate ≈ `crystal_full` hitrate | Pose_solver v2 geometry is fit for purpose; pose itself is not the bottleneck |
| `crystal_*` ≫ `posed_*` across the row | Pose_solver v2 geometry still imperfect — over-packing residual is meaningful |
| `crystal_full` ≫ `crystal_auto*` | Full_test_run contig is well-tuned; auto-derivation is not enough |
| `*_full` < `*_auto5A` | The hand-designed contig is missing important contact residues (e.g. A3) |
| `*_auto3A` ≪ `*_auto5A` | Fragmented short segments hurt RFDiffusion; broader design region needed |
| All cells low hitrate | Something further upstream (e.g. Boltz folding the receptor) is broken |

---

## 9. What's deferred

- **Pia analyses 1 (sidechain clashes) and 2 (misfold localisation)** from `pipeline_notes15.md` — files `experiments/campaigns/pikp1_avrpia/analysis1_sidechain_clashes.py` and `analysis2_misfold_localisation.py` exist (created in a prior session) but not yet run / not integrated here.
- **Re-solving Pia at max-pair-distance 6.0** to give Pia its v2 pose equivalent. The current pose_solved Pia is still over-packed. Worth doing if the PikF factorial implicates pose quality.
- **Heavy-atom-aware pair constraints** (LogSumExp min-distance instead of CA-CA). Would close the over-packing loophole structurally. Larger code change.
- **Global sidechain-clash penalty** (vs only on pair residues). The current pair-only check is necessary but not sufficient — see v3 result.

---

## 10. File index

New / modified this session:

| Path | Status |
|---|---|
| `tests/haddock/pose_solver/pose_solver.py` | modified (sidechain clash penalty + bug fix) |
| `tests/haddock/pose_solver/interface_metrics.py` | new |
| `tests/haddock/pose_solver/derive_contigs.py` | new |
| `tests/haddock/pose_solver/Pikp-1_HMA_PikF/` | new (v1 PikF pose, 3-pair max 4 Å, over-packed) |
| `tests/haddock/pose_solver/Pikp-1_HMA_PikF_v2/` | new (v2 PikF pose, 3-pair max 6 Å, A38-B42 — used in factorial) |
| `tests/haddock/pose_solver/PikF_2pair_min4_max6_scclash/` | new (v3 PikF pose, 2-pair with sc-clash; worse) |
| `tests/haddock/pose_solver/interface_metrics_out/` | new (comparison json/txt outputs) |
| `tests/haddock/data/PikF.pdb` | new (effector for solver, chain B of pik_interface_complex) |
| `experiments/campaigns/pikp1_avrpikf/inputs/pikp1_avrpikf_pose_solved_v2.pdb` | new (POSED input for factorial) |
| `experiments/campaigns/pikp1_avrpikf/runs/{crystal,posed}_{full_test,auto5A,auto3A}_contig/params.yml` | new (6 campaign runs) |
| `notes/pipeline_notes/pipeline_notes16.md` | new (this file) |
