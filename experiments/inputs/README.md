# experiments/inputs/

Shared library of structural and sequence inputs used across campaigns.
PDBs, UniProt FASTAs, reference structures, and any other input asset
that more than one campaign might consume.

## Naming convention

Files should be named `<system>_<source>_<descriptor>.<ext>`, lowercase,
underscore-separated. Examples:

- `pikp1_alphafold_monomer.pdb`
- `pikp1_uniprot_q6k4q3.fasta`
- `avrpia_pdb_6r8m_chainB.pdb`

The `<source>` component should make the provenance scannable from the
filename alone (`alphafold`, `pdb`, `uniprot`, `colabfold`, `manual`,
etc.). When in doubt, err on the side of a longer, more specific name.

## Provenance metadata

Every file should be accompanied by a sibling `<basename>.meta.yml` with
at minimum:

```yaml
source_url: https://...        # download URL or DOI
source_id: 6R8M                # PDB ID, UniProt accession, etc.
acquired: 2026-05-03           # date downloaded (YYYY-MM-DD)
acquired_by: jowillia
notes: |                       # any post-download processing
    Stripped waters and HETATMs; renumbered chain B to start at 1.
```

If a file is derived from another (e.g. a chain extraction or a
renumbering), record the parent file in `notes:` and put the derivation
script under `scripts/` rather than running it ad-hoc.

## Subdirectories

Group by asset type as the library grows:

- `structures/` — PDB / mmCIF files.
- `sequences/` — FASTA files.
- `references/` — curated reference structures used for RMSD / metric
  calculations.

Subdirectories beyond these should only be added when there is a clear
category that does not fit existing buckets.
