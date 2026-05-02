# tests/rosetta_filtering/data/

Canonical inputs for `test_rosetta_filtering.nf` — split two-chain
PDBs from a discovery-run `RFDIFFUSION_FILTER` step's `split/`
directory, hand-picked to span Sc above and below
`params.sc_threshold = 0.5`.

## Layout

```
data/
└── rfdiffusion_split/
    ├── design_0.pdb     <- strong pass  (sc=0.569)
    ├── design_1.pdb     <- marginal pass (sc=0.504, just above 0.5)
    ├── design_14.pdb    <- boundary fail (sc=0.488, just below 0.5)
    └── design_59.pdb    <- clear fail   (sc=0.405)
```

Per-file shape:

- Two chains: A = receptor (polyvaline backbone from RFDiffusion),
  B = effector. Renumbered from 1.
- Filename pattern `design_<N>.pdb` (matched by the `design_*.pdb`
  glob in `test_rosetta_filtering.nf:34`).

## Path coverage

- **Pass branch** of `bin/rosetta_filter_collect.py`: exercised by
  `design_0` (clearly above threshold) and `design_1` (marginal).
- **Fail branch**: exercised by `design_14` (marginal) and `design_59`
  (clearly below threshold).
- The marginal cases either side of `sc_threshold=0.5` exercise the
  threshold gate exactly; the clear cases exercise the typical
  pass/fail metric ranges.

See `notes/inventory/15_discovery_run_path_coverage.md` §Rosetta
filtering for the candidate-selection rationale and verification
commands.

## Source

Discovery run `tests/full_test_run/results/rfdiffusion/split/`,
discovery date 2026-05-02, params per `params_full_test.yml`.
