# experiments/scripts/

Reusable Python helpers shared across campaigns. Importable as a package
(`from experiments.scripts import …`) — the empty `__init__.py` marks it
as such.

What lives here: small, campaign-agnostic utilities. Examples of things
that will accumulate here over time:

- A `sys.path` bootstrap that exposes `bin/` to campaign scripts (future
  prompt).
- I/O helpers for reading campaign result trees.
- Plotting primitives shared across campaign analyses.

What does NOT live here: anything campaign-specific (that goes under
`campaigns/<name>/`), and anything substantial enough to warrant its own
subsystem — `param_derivation/` is the model for the latter.
