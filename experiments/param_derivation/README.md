# experiments/param_derivation/

Code for deriving pipeline input parameters from upstream data, with an
emphasis on RFDiffusion contig strings. Kept separate from `scripts/`
because contig derivation is a substantial concern with its own logic
(hotspot selection, chain mapping, motif segmentation) rather than a
small shared utility.

What will live here: routines that take a receptor PDB / hotspot list /
HADDOCK output and emit the parameter values a campaign feeds to the
production pipeline (chiefly `rfdiff_contigs`, but also adjacent
parameters like designed-region indices and chain assignments).

Read alongside the production contig code in `bin/contig_utils.py` and
`bin/build_contigs.py` — the goal here is exploratory derivation logic
that may eventually feed back into those, not a parallel implementation.
