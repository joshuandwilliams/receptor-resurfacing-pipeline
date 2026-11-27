# ⛔ STOP — THIS ANALYSIS IS **NOT** THE THESIS CHAPTER 5 CAMPAIGN ⛔

## READ THIS BEFORE USING ANY FILE IN THIS DIRECTORY

Everything here belongs to the run

    crystal_full_contig_d26_seqvar

This is a **sequence-variation run on the single backbone design_26**.
It is a **DIFFERENT RUN** from the main AVR-PikF campaign written up in
thesis Chapter 5.

## How to recognise that you are in the wrong place

Every row in the tables here looks like this:

- `mpnn_sequence` values are all `design_26_seq_<N>` — one backbone, 127 sequences
- `corrected_receptor` is **84 aa for every row, with no exceptions**
- the design region 1 field of `designed_residues` is **19 residues for every row**
- the gBlocks are **280 bp for every row**, core 252 bp

Those constants are **correct for this run**. They are constant because this
run varied the sequence on one fixed backbone.

## Why that matters

The thesis Chapter 5 ten selected designs are named `d39_s0`, `d17_s0`,
`d47_s0`, `d39_s1`, `d26_s1`, `d18_s0`, `d49_s0`, `d40_s0`, `d6_s1`,
`d14_s0`. They come from **nine different RFdiffusion backbones** and have
**design region 1 lengths from 10 to 19 residues**.

None of that is in this directory. The `design_26_seq_<N>` naming here is
**not** the `d<backbone>_s<seq>` naming used in the thesis. Do not map one
onto the other.

If you compare the constant 19-residue design region 1 in these files against
the 10-to-19 range in thesis `@tbl-selected-designs`, you will conclude there
is a contradiction in the ordered constructs. **There is not.** You are simply
reading the wrong run.

## What this directory IS good for

`../../../../../scripts/make_gblocks.py` is the shared reverse-translation and
Golden Gate tailing script, and it is the same code used for the real campaign.
The codon table, the tails and the domestication logic can be read from here.
The **sequences and the lengths** cannot.
