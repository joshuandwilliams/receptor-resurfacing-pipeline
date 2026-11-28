# Known defects

Every entry here was found while writing the characterization suite, verified
by reading the code, and left unfixed on purpose. Each one is pinned by a test,
so it cannot change without a test changing too.

Two kinds of pin are used.

**An `xfail`** states the intended behaviour and is expected to fail. It is
paired with a passing test of today's behaviour, so the defect cannot be
half-fixed silently: fixing it flips the `xfail` to a pass and breaks its
partner, which forces both to be looked at.

**A passing test** asserting today's behaviour is used where there is no clean
statement of intent, only a design decision worth recording so it is not
mistaken for a bug by the next reader.

`xfail`s are `strict=False`, so a fix does not fail the suite. Running
`pytest -m logic -ra` prints every reason.

## Why these are not fixed

Six of the eight change what a stage writes. The AVR-PikF campaign's
`passing_summary.csv`, `cross_sequence_summary.csv` and
`rfdiffusion_metrics.json` are committed in
[`resurfacing-experiments`](https://github.com/joshuandwilliams/resurfacing-experiments)
and are the primary source behind the Chapter 5 figures. Fixing a defect that
changes one of those columns means the figures and the tables stop agreeing,
and there is no way to tell from the file which version produced it.

So the order is: pin the defect, then fix it in a change that also re-runs the
campaign and replaces the committed outputs together. Not before.

Defects 4 and 8 are the exceptions, and both can be fixed on their own.

MMseqs2 was present when the campaign ran, so defect 4's missing-binary path
was never taken and no committed output depends on it.

Defect 8 only fires on a contig the grammar does not allow. Every campaign
contig is valid, and `derive_pose_contigs.strat4a_range` cannot emit an
inverted range, so refusing one changes no committed output either.

The evidence for that is the campaign's own committed cluster counts. The
fallback reports `n_clusters == n_sequences`, and
`crystal_full_test_contig.mpnn_cluster_counts.csv` reports 40, 55, 64 and 78
clusters out of 128 sequences at the four thresholds, which only real
clustering produces. `HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.def:98-102`
installs it and puts `/opt/mmseqs/bin` on `PATH`.

---

## 1. Tied cross-ranks are ordered by heap address

`bin/cross_sequence_summary.py:623`, `:633`

Both rank keys end in `id(rec)`:

```python
return (tier_order, -comp, id(rec))
```

`id()` is the object's memory address. For two records in the same tier with
an exactly equal composite score, which one gets the lower
`cross_rank_by_composite` depends on where CPython happened to allocate them.
Re-running the same stage on the same input can swap them.

The module's own comment says so. Input order would be deterministic and is
what a reader assumes a stable rank means.

Pinned by `tests/characterization/test_cross_sequence_summary.py::test_tied_cross_ranks_must_not_depend_on_object_identity`.

## 2. A `nan` median is reported as a confident prediction

`bin/extract_passing.py:219-239`

Every confidence check has the shape:

```python
if float(out.get("iptm_median", 1)) < iptm_min:
    triggers.append("low_iptm")
```

`float("nan") < anything` is `False` by IEEE 754, so a `nan` median triggers
nothing. With no triggers the function returns `"ok"`, which is the same value
a genuinely confident sequence gets. A sequence whose confidence could not be
measured is indistinguishable in `confidence_flag` from one that passed.

The `try`/`except (ValueError, TypeError)` around each check catches an
unparseable string but not a parseable `nan`.

Pinned by `tests/characterization/test_extract_passing.py::test_a_nan_metric_must_not_be_reported_as_confident`.

## 3. A `nan` DockQ metric leaves the row looking complete

`bin/compute_interface_metrics.py:390-398`

```python
if dq_err:
    failures.append(dq_err)
else:
    out["irmsd"] = f"{irmsd:.3f}" if not math.isnan(irmsd) else ""
```

When DockQ returns without an error but with `nan` for a metric, the column is
written blank and nothing is appended to `failures`, so `dockq_failures` ends
up empty. The row then carries blank metrics and no stated reason, and
`main()` counts it as fully populated.

The module's own docstring at line 58 says "Failures are flagged, never
silently dropped". This is the one path where that does not hold.

Pinned by `tests/characterization/test_interface_and_manifest_mains.py::test_a_nan_metric_must_not_be_reported_as_a_success`.

## 4. A missing MMseqs2 binary fails the stage instead of falling back

`bin/mpnn_cluster_sequences.py:87`

`run_mmseqs_cluster` catches `subprocess.CalledProcessError` and nothing else.
A binary that is not on `PATH` raises `FileNotFoundError`, which propagates out
of `main()` and fails `MPNN_CLUSTER`.

The module's stated design is to warn and fall back to reporting `n_seqs`, and
a container built without MMseqs2 is the likelier of the two failures.

Pinned by `tests/characterization/test_mpnn_cluster_sequences.py::test_an_absent_mmseqs_binary_should_fall_back_rather_than_raise`.

## 5. An empty MPNN score cell fails the whole QC stage

`bin/mpnn_sequence_qc.py:134`

```python
score = float(row.get("mpnn_score", "999"))
```

`dict.get` returns the default only when the key is absent. A column that is
present but empty gives `""`, which reaches `float("")` and raises
`ValueError` out of `main()`. One unparsed score fails `SEQUENCE_QC` for the
whole cohort.

The `999` default exists to treat an unknown score as worst-possible. Applying
it to an empty cell too would match that intent.

Pinned by `tests/characterization/test_mpnn_sequence_qc.py::test_an_unparseable_mpnn_score_should_not_fail_the_stage`.

## 6. A contig that fixes nothing frees the entire receptor, and exits 0

`bin/pipeline_correct_sequences.py:391`, writer at `:473`

`parse_contig_segments` returning an empty list is indistinguishable from a
contig that genuinely fixes nothing. `generate_fixed_positions` carries on,
`fixed_resnums` stays empty, and it writes a `fixed_positions_*.jsonl` that
frees every receptor residue. ProteinMPNN then redesigns the whole receptor
instead of the design region, and the stage exits 0.

A contig that resolves to no segments for the named receptor chain is always a
specification error. Aborting, or at minimum warning, would catch a mistyped
chain letter at the stage that made it rather than in the metrics forty
minutes later.

Pinned by `tests/characterization/test_pipeline_correct_sequences.py::test_an_empty_segment_list_should_be_reported_as_an_error`.

## 7. The de novo length split is a guess when more than one segment varies

`bin/rfdiffusion_filter.py`, `build_receptor_resnum_map`
Also `bin/pipeline_correct_sequences.py:455-467`, same arithmetic

The design PDB records only the TOTAL number of de novo residues. With one
variable-length segment the total IS that segment's length, so the mapping is
exact. With two or more, the total has to be divided, and the code divides it
proportionally by `max_len`:

```python
actual_len = round(total_denovo_actual * seg['max_len'] / total_max)
```

A design whose real regions are 4 and 8 cannot be told from 6 and 6. Every
residue after the first mis-sized segment is shifted, and everything keyed off
the residue map shifts with it: `motif_rmsd`, the per-region displacements,
and the design-versus-fixed contact split.

**This one has been measured.** On the AVR-PikF campaign the auto-derived
contigs gave wrong lengths for 62 of 64 and 54 of 64 designs, inflating motif
RMSD to 4.0 and 2.6 Å where RFdiffusion's own mapping gives 0.22 and 0.11 Å.
The retained scaffold had not moved at all. Recorded in
`build_receptor_resnum_map`'s docstring, found 2026-08-29.

The fix is not a better guess. RFdiffusion already writes a `.trb` file
carrying `sampled_mask` and the exact input-to-output correspondence, and the
pipeline does not publish it. Publishing it removes the guess entirely. That
is a change to `modules/`, not to `bin/`, which is why it is pinned here.

Pinned by `tests/characterization/test_rfdiffusion_filter_mapping.py::test_two_variable_segments_split_the_budget_correctly`.

## 8. Only one of the two contig parsers validates the grammar

`bin/rfdiffusion_filter.py`, `_parse_receptor_segments`

`contig_spec.DeNovoSegment` refuses a segment whose maximum is below its
minimum:

```python
raise ValueError(f"DeNovoSegment max_len={self.max_len} < min_len={self.min_len}")
```

`_parse_receptor_segments` reads the same grammar with its own code and
accepts it. A mistyped `1-0` becomes `("denovo", 1, 0)`, which the budget
split then treats as a variable-length segment whose maximum is zero. The
proportional division divides by a total maximum of zero, falls back to an
even split, clamps each share up to the minimum, and hands the remainder to
the last segment. So a contig that cannot be satisfied still produces a
residue map, and everything keyed off that map is computed from it.

Refusing it means raising from a parser the filter calls once per design, so
a typo that reaches a cluster run would fail the stage instead of writing
wrong numbers. That is the right behaviour and it is a change to what the
stage does, not to a test, which is why it is pinned here.

Pinned by `tests/characterization/test_rfdiffusion_filter_mapping.py::test_an_inverted_denovo_range_should_be_refused`, with
`test_an_inverted_denovo_range_is_accepted_as_written` asserting today's
behaviour.

---

# Design decisions that look like defects

These are deliberate. They are listed so the next reader does not "fix" one.

## `--filter-plddt-min` is required and has no effect

`bin/merge_orthogonal_metrics.py:133-135`, `:154`

`interface_pLDDT` is a Boltz-2 confidence metric, and negative steering
optimises against Boltz-2. Gating on it would be scoring a prediction against
itself, which is the opposite of an orthogonal check, so it is deliberately
not checked. The flag stays `required=True` so the Nextflow caller does not
have to change.

The consequence is real: no `plddt` flag is ever emitted, and raising
`--filter-plddt-min` to 0.99 changes nothing. Pinned by
`test_orthogonal_merge_and_collect.py::test_interface_plddt_is_never_checked_and_its_threshold_is_dead`.

## `--min-hotspot-frac` does not read the hotspot list

`bin/rfdiffusion_filter.py:811`, `:836`

The name is wrong rather than the code. `frac_in` is the fraction of all
receptor-effector contacts whose receptor residue lies OUTSIDE the fixed set,
which is to say inside the design region:

```python
n_in = sum(1 for c in contacts
           if resnum_map.get(c[0], -1) >= 0
           and resnum_map[c[0]] not in fixed_residues)
```

`hotspot_residues` is parsed from `--hotspot`, written into
`rfdiffusion_metrics.json`, and never consulted by the gate. The default is
`0.0`, so the gate never fires on the committed campaign either way.

## AF3 disagreement is flagged but does not gate

`bin/merge_orthogonal_metrics.py`

It gated until 2026-04-28. Negative steering optimises against Boltz, so an
AF3 disagreement is a sanity-check warning rather than a kill criterion. The
flag `af3_nomsa_ra_eff_too_high` is emitted and `passes_orthogonal_filters`
stays `1`.

## Rosetta ddG is carried but never thresholded

`bin/merge_orthogonal_metrics.py`

A wildly favourable and a wildly unfavourable ddG both pass. The column is
reported for inspection, not used as a filter.

## `mpnn_sequence_qc` requires `--receptor-seq` and never reads it

`bin/mpnn_sequence_qc.py:29`

`args.receptor_seq` appears nowhere in the module. `SEQUENCE_QC` interpolates
the entire native receptor sequence onto the command line to satisfy it. Not a
correctness problem, but it is a long argument on every invocation for
nothing.

## An empty result writes a short, fixed header

`bin/mpnn_sequence_qc.py:177`, `bin/mpnn_design_region_score.py:253`

On the normal path both modules write whatever columns the input carried. When
nothing passes, both write a hardcoded header instead, seven columns in one
and four in the other. So an empty output file has a different schema from a
populated one.

Nothing downstream reads columns from an empty file today, which is why it has
not bitten.
