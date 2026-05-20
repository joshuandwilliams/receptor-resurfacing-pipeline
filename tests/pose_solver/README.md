# Constraint-driven pose solver

Finds the rigid-body placement of the effector relative to the receptor
that satisfies user-specified CA–CA contact constraints, **without any
clash penalty**.  The resulting clashes at the solved pose are the
biologically meaningful design region — the receptor residues that need
to be redesigned by RFDiffusion to accommodate the effector.

## Quick start

```bash
cd tests/pose_solver
bash run_example.sh
```

The Nextflow integration lives in `modules/pose_solver.nf` and is
exercised by `tests/pose_solver/test_pose_solver.nf`.  Use the shell
script below to iterate on pair choices outside the pipeline.

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
python3 ../../bin/pose_solver.py \
    --binder    ./data/Pikp-1_HMA.pdb \
    --target    ./data/avr-pia.pdb \
    --pairs "A73-B31 A71-B33 A8-B22" \
    --min-pair-distance 3.5 \
    --max-pair-distance 6.0 \
    --contig-design-region "33-49,69-78" \
    --n-restarts 1000 \
    --global-interp
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
