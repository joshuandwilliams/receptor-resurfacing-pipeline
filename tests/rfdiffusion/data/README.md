# tests/rfdiffusion/data/

Canonical inputs for `test_rfdiffusion.nf`. RFDiffusion is the first
stage in the pipeline, so this directory holds *primary* test inputs
(no upstream module produces them).

## Contents

| File | Shape | Source |
|---|---|---|
| `af3_pikp1_native_avrpikf_complex.pdb` | Receptor+effector complex PDB. Chain A = receptor (Pikp-1), chain C = effector (AvrPikF). | Hand-curated from an AF3 prediction; committed to repo. |
| `af3_pikp1_native_avrpikf_complex.pdb.original_2model` | Two-model PDB this single-model file was derived from. Reference only; the test reads only the `.pdb`. | Same source. |

## Curation status

- **No discovery-run input is needed for this test.** The required
  fixtures (input PDB, contigs string baked into `test_rfdiffusion.nf`,
  chain labels) are independent of the discovery run's output.
- The reference output set for `tests/rfdiffusion/example_output_files/`
  will be built later, in Phase 2.8, by running this test on its own
  `data/` inputs and capturing the result.

See `notes/inventory/15_discovery_run_path_coverage.md` §RFDiffusion
for the path-coverage analysis. Note the report's Coverage Gap 1 —
`passes_filter == false` is not exercisable from this fixture (no
discovery run has produced a failed RFDiffusion design); that branch
will be covered by a Phase 2.9 unit-style test, not by this fixture.

## Notes

- The contigs string and chain labels are wired into
  `test_rfdiffusion.nf` directly; only the PDB file lives here.
- Replace the PDB by overriding `--pdb_file` on the command line if
  you need to test against a different complex.
