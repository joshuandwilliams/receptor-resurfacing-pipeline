# Constraint-driven pose solver

Finds the rigid-body placement of the effector relative to the receptor
that satisfies user-specified CA–CA contact constraints, **without any
clash penalty**.  The resulting clashes at the solved pose are the
biologically meaningful design region — the receptor residues that need
to be redesigned by RFDiffusion to accommodate the effector.

## Quick start

```bash
cd tests/haddock/pose_solver
bash run_example.sh
```

Opens `solved_pose_posed.pdb` in ChimeraX:
```
open solved_pose_posed.pdb
```

## How it works

1. You specify 2–4 pair constraints `ARESNUM-BRESNUM@TARGET_DIST` (CA–CA Å).
2. The solver minimises `Σ (|dist(r_i, T(e_i))| − target_i)²` over the 6
   rigid-body DOFs (rotation + translation) of the effector.
3. No clash penalty — the effector lands wherever the pair constraints say.
4. Heavy-atom clashes (<2 Å) at the solved pose are reported as the
   **design region**.

## Usage

```bash
python3 pose_solver.py \
    --receptor  ../data/receptor.pdb \
    --effector  ../data/effector.pdb  \
    --pairs "A71-B33 A72-B32 A73-B31" \
    --target-distance 3.8             \
    --contig-design-region "33-49,69-78" \
    --n-restarts 200
```

### Key flags

| Flag | Default | Meaning |
|---|---|---|
| `--pairs` | required | Space-separated `ARESNUM-BRESNUM` or `ARESNUM-BRESNUM@DIST` |
| `--target-distance` | 3.8 | Default CA–CA target (Å) when `@DIST` not specified |
| `--clash-cutoff` | 2.0 | Heavy-atom distance counted as a clash |
| `--contact-cutoff` | 8.0 | CA–CA distance shown on heatmap |
| `--contig-design-region` | "" | Residue ranges already planned for redesign, e.g. `"33-49,69-78"` |
| `--n-restarts` | 200 | More restarts → better global minimum (slower) |
| `--output-prefix` | solved_pose | Prefix for output files |

## Typical target distances

| Contact type | CA–CA (Å) |
|---|---|
| β-strand H-bond (anti-parallel) | 3.8–4.8 |
| Close sidechain contact | 4.5–6.0 |
| Electrostatic / polar interaction | 5.0–8.0 |
| Hotspot "within reach" | 6.0–9.0 |

## Outputs

| File | Contents |
|---|---|
| `*_posed.pdb` | Complex with effector transformed to constraint-optimal position |
| `*_results.json` | Pair distances (target vs achieved), clash set, design region |
| `*_heatmap.png` | Contact heatmap; stars=pair pins; red borders=clashes; orange bar=design region |

## Requirements

```bash
pip install numpy scipy matplotlib
```
