# The analyses moved

The twelve Chapter 5 analyses that were here now live in
[`receptor-resurfacing-analyses`](https://github.com/joshuandwilliams/receptor-resurfacing-analyses),
with their history, so this repository holds the pipeline and that one holds
what was concluded from running it.

Same split as
[`negative-steering-analyses`](https://github.com/joshuandwilliams/negative-steering-analyses).

## The run data stayed here

`../runs/crystal_full_test_contig/results/` and `.../agroinfiltrations/` did
not move. They are gitignored, they have never been tracked, and they are
roughly 4.5 GB. They are also the primary source behind every Chapter 5 figure
and number, so they are the one thing in this campaign that cannot be
regenerated.

The analyses find them through `analyses/_run_root.py` in that repository,
which checks `$RECEPTOR_RUN_ROOT`, then a `pipeline/` submodule, then a
sibling checkout of this repository. The sibling case is what works on the
machine they were written on, so a side-by-side clone needs no configuration.
