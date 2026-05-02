# tests/proteinmpnn/data/

Canonical inputs for `test_proteinmpnn.nf`. Two Rosetta-passing split
PDBs hand-picked because their downstream MPNN sequences span the
negsteer outcome partition (one feeds the top-ranked
`cold_start_all_clean` survivor, one feeds a tier-none / poorly-placed
sequence — the per-module test stops at MPNN_SELECT_TOP, so the
downstream split surfaces only via the per-MPNN-sequence metadata).

## Layout

```
data/
├── input_complex.pdb              <- receptor+effector complex
└── rosetta_passing/
    ├── design_0.pdb               <- parent for design_0_seq_0 / _seq_1
    └── design_28.pdb              <- parent for design_28_seq_0 / _seq_1
```

### `input_complex.pdb`

The same receptor+effector complex PDB that fed RFDiffusion. Used by
`RESOLVE_CONTIGS` and `EXTRACT_SEQUENCES` to derive the receptor
sequence, effector sequence, and `receptor_start_pdb` offset.

Identical to `tests/rfdiffusion/data/af3_pikp1_native_avrpikf_complex.pdb`;
the file is duplicated rather than symlinked to keep test fixtures
self-contained per `notes/inventory/14_phase_2_revision_per_module_tests.md` §2.4.

### `rosetta_passing/design_<N>.pdb`

Split two-chain PDBs from `tests/full_test_run/results/rosetta_filtering/passing/`:

| File | sc_value | Rationale |
|---|---:|---|
| `design_0.pdb`  | 0.569 | Parent design for `design_0_seq_0` (Class 3a — `steered_clean` not tier-eligible) and `design_0_seq_1` (Class 1 — `cold_start_all_clean`). |
| `design_28.pdb` | 0.567 | Parent design for `design_28_seq_1` (Class 1, cross-rank 1 in the discovery-run cohort — top representative). |

Per-file shape:
- Two chains: A = receptor backbone, B = effector. Renumbered from 1.
- Filename pattern `design_<N>.pdb` (matched by the `design_*.pdb`
  glob in `test_proteinmpnn.nf:39`).

## Path coverage

- The MPNN per-module test exercises MPNN_FIXED_POSITIONS →
  PROTEINMPNN → SEQUENCE_CORRECTION → SEQUENCE_QC →
  MPNN_DESIGN_REGION_SCORE → MPNN_CLUSTER → MPNN_PLOTS →
  MPNN_SELECT_TOP. The "selected for top-N" vs "QC-dropped" partition
  is exercised by running the workflow on these two parent designs
  and observing which of the resulting `num_seqs * 2 = 8` sequences
  pass `mpnn_top_n=8` selection.

See `notes/inventory/15_discovery_run_path_coverage.md` §ProteinMPNN
for candidate-selection rationale and verification commands. Note the
path-correction in that section: MPNN's `top_metadata.csv` is
published to `${params.outdir}/sequences/`, not `${params.outdir}/mpnn/`.

## Source

Discovery run `tests/full_test_run/results/`, discovery date
2026-05-02, params per `params_full_test.yml`.
