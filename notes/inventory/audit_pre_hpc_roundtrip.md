---
title: Pre-flight audit — Wave 1 hpc-tier round-trip
date: 2026-05-02
scope: analytical only (no code changes)
---

# Pre-flight audit before HPC round-trip

Verifies the codebase no longer depends on three pieces of pre-git environment
state that were cleaned up by hand:

- the in-tree `jdk-17.0.2/` directory (moved to `~/singularity/jdk-17.0.2/`)
- the project-local `nxf_home/` directory (deleted; Nextflow now defaults to `$HOME/.nextflow`)
- repo-root symlinks for `slurm` and `java` tools (removed)

Each section grep'd `**/*.nf`, `**/*.config`, `**/*.sh`, `**/*.py`, and
`containers/*.def`.

---

## PART 1 — JDK references

### Findings

On HPC, the macOS mount `/Volumes/HPC-Home/` corresponds to `/hpc-home/jowillia/`.
`PIPELINE_DIR` in every slurm launcher resolves to
`/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline`.
The old in-tree JDK (`${PIPELINE_DIR}/jdk-17.0.2`) is gone; the new
canonical location is `/hpc-home/jowillia/singularity/jdk-17.0.2/`.

1. **`run_pipeline.slurm.sh:50`** — `export JAVA_HOME="${PIPELINE_DIR}/jdk-17.0.2"`
   - Resolves to `/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/jdk-17.0.2` — **BROKEN** (directory no longer exists at that path).
   - Fix: `export JAVA_HOME="/hpc-home/jowillia/singularity/jdk-17.0.2"` (or, more portably, derive it from `$HOME/singularity/jdk-17.0.2`).
2. **`run_pipeline.slurm.sh:51`** — `export PATH="${JAVA_HOME}/bin:${PATH}"`
   - Inherits the broken `JAVA_HOME` from the line above. **BROKEN** (transitively).
   - Fix: same as #1; this line is correct once `JAVA_HOME` is fixed.
3. **`run_pipeline.slurm.sh:87`** — `echo "JAVA_HOME:       ${JAVA_HOME}"`
   - Diagnostic echo only. Not broken on its own, but will print the wrong path until #1 is fixed.
4. **`tests/rfdiffusion/run_test_rfdiffusion.slurm.sh:26`** — `export JAVA_HOME="${PIPELINE_DIR}/jdk-17.0.2"` — **BROKEN**.
5. **`tests/rfdiffusion/run_test_rfdiffusion.slurm.sh:27`** — `export PATH="${JAVA_HOME}/bin:${PATH}"` — **BROKEN** (transitively).
6. **`tests/haddock/run_test_haddock.slurm.sh:26`** — same pattern — **BROKEN**.
7. **`tests/haddock/run_test_haddock.slurm.sh:27`** — **BROKEN** (transitively).
8. **`tests/orthogonal_metrics/run_test_orthogonal_metrics.slurm.sh:51`** — same pattern — **BROKEN**.
9. **`tests/orthogonal_metrics/run_test_orthogonal_metrics.slurm.sh:52`** — **BROKEN** (transitively).
10. **`tests/negative_steering/run_test_negative_steering.slurm.sh:39`** — same pattern — **BROKEN**.
11. **`tests/negative_steering/run_test_negative_steering.slurm.sh:40`** — **BROKEN** (transitively).
12. **`tests/rosetta_filtering/run_test_rosetta_filtering.slurm.sh:26`** — same pattern — **BROKEN**.
13. **`tests/rosetta_filtering/run_test_rosetta_filtering.slurm.sh:27`** — **BROKEN** (transitively).
14. **`tests/proteinmpnn/run_test_proteinmpnn.slurm.sh:24`** — same pattern — **BROKEN**.
15. **`tests/proteinmpnn/run_test_proteinmpnn.slurm.sh:25`** — **BROKEN** (transitively).
16. **`containers/NextFlow.def:7`** — `export JAVA_HOME=/usr/lib/jvm/java-17-openjdk`
   - Path inside the Singularity image, populated by `yum install -y java-17-openjdk-headless` on line 13. **OK** — unrelated to the host-side JDK move.
17. **`containers/NextFlow.def:8`** — `export PATH=$JAVA_HOME/bin:/usr/local/bin:$PATH` — **OK** (in-container).
18. **`containers/NextFlow.def:13`** — `yum install -y java-17-openjdk-headless ...` — **OK** (in-container provisioning).
19. **`containers/NextFlow.def:16`** — `export JAVA_HOME=/usr/lib/jvm/java-17-openjdk` (in `%post`) — **OK** (in-container).
20. **`containers/NextFlow.def:17`** — `export PATH=$JAVA_HOME/bin:$PATH` (in `%post`) — **OK** (in-container).

`bin/`, `*.py`, and `nextflow.config` had no JDK references.
A second case-insensitive sweep on `jdk` / `JDK` confirmed there are no other matches.

### Action items

- Update **all 7** slurm launchers (1× root, 6× tests) to point `JAVA_HOME` at the new singularity-tree JDK before the round-trip. Recommend a single shared variable (e.g. `JDK_DIR="${HOME}/singularity/jdk-17.0.2"`) so future moves only require one edit.
- Container provisioning (`containers/NextFlow.def`) is unaffected.

---

## PART 2 — `nxf_home` references

### Findings

The project-local `nxf_home/` directory was deleted, but every slurm
launcher still computes `NXF_HOME="${PIPELINE_DIR}/nxf_home"`, then runs
`mkdir -p "${NXF_HOME}"` and seeds the directory from the Nextflow
container. Net effect on first HPC run: the deleted directory is
**recreated in-tree** rather than Nextflow defaulting to `$HOME/.nextflow`.
Nothing breaks per se, but the cleanup is silently undone.

1. **`run_pipeline.slurm.sh:56`** — `export NXF_HOME="${PIPELINE_DIR}/nxf_home"` — **NEEDS-CHECK** (will recreate the deleted dir).
2. **`run_pipeline.slurm.sh:60`** — `mkdir -p "${NXF_HOME}" "${NXF_WORK}" "${NXF_TEMP}"` — recreates `nxf_home/`.
3. **`run_pipeline.slurm.sh:63`** — `NEXTFLOW_BIN="${NXF_HOME}/nextflow"` — depends on the recreated dir.
4. **`run_pipeline.slurm.sh:71-74`** — plugin-cache seeding into `${NXF_HOME}/plugins`.
5. **`run_pipeline.slurm.sh:84`** — diagnostic echo.
6. **`tests/rfdiffusion/run_test_rfdiffusion.slurm.sh:21,32,36,39,47-49`** — same `nxf_home`-recreation pattern. **NEEDS-CHECK**.
7. **`tests/haddock/run_test_haddock.slurm.sh:21,32,36,39,47-49`** — same. **NEEDS-CHECK**.
8. **`tests/orthogonal_metrics/run_test_orthogonal_metrics.slurm.sh:46,57,61,64,72-74`** — same. **NEEDS-CHECK**.
9. **`tests/negative_steering/run_test_negative_steering.slurm.sh:17` (comment),`34,45,49,52,60-62`** — same. The comment on line 17 ("using the cached binary under nxf_home/") explicitly documents the in-tree convention. **NEEDS-CHECK**.
10. **`tests/rosetta_filtering/run_test_rosetta_filtering.slurm.sh:21,32,36,39,47-49`** — same. **NEEDS-CHECK**.
11. **`tests/proteinmpnn/run_test_proteinmpnn.slurm.sh:20,29,33,35,42-44`** — same. **NEEDS-CHECK**.
12. **`containers/NextFlow.def:18`** — `export NXF_HOME=/opt/nextflow` (in `%post`, inside the image) — **OK** (unrelated to the host directory).
13. **`containers/NextFlow.def:20`** — `mkdir -p $NXF_HOME/plugins` (in-container) — **OK**.
14. **`containers/NextFlow.def:32`** — `export NXF_PLUGINS_DIR=$NXF_HOME/plugins` (in-container) — **OK**.
15. **`containers/NextFlow.def:46-47`** — `%test` block, `/tmp/.nxf_test` — **OK**.

No code references `NXF_HOME` outside slurm shell scripts and the container def.

### Action items

- **Decision required, not a blocker for the round-trip.** The slurm launchers will recreate the project-local `nxf_home/` on first run. If the cleanup was intentional (use Nextflow's default `$HOME/.nextflow`), redirect `NXF_HOME` accordingly in the launchers (and update the binary-extraction + plugin-seed blocks to match). If the cleanup was just to keep the repo tidy in git but you're fine with it being recreated at runtime, no change needed — just ensure `nxf_home/` is in `.gitignore`.
- The Wave 1 round-trip will work either way; flag for follow-up.

---

## PART 3 — Slurm tool references

### Findings

1. **`run_pipeline.slurm.sh:18`** — `echo "Usage: sbatch run_pipeline_slurm.sh ..."` — usage string. Bare on `$PATH`. **OK**.
2. **`run_pipeline.slurm.sh:89`** — `echo "sbatch:          $(which sbatch)"` — diagnostic. Bare on `$PATH`. **OK**.
3. **`bin/boltz2_iterate_steering.py:30`** — docstring mentioning `sbatch`. **OK**.
4. **`bin/boltz2_iterate_steering.py:798`** — comment. **OK**.
5. **`bin/boltz2_iterate_steering.py:988`** — docstring (`sbatch`). **OK**.
6. **`bin/boltz2_iterate_steering.py:991`** — comment. **OK**.
7. **`bin/boltz2_iterate_steering.py:1092`** — comment ("sbatch is on PATH"). Confirms bare-on-PATH invocation downstream. **OK**.
8. **`bin/boltz2_iterate_steering.py:1910`** — docstring. **OK**.
9. **`bin/boltz2_iterate_steering.py:2002`** — docstring. **OK**.
10. **`bin/boltz2_iterate_steering.py:2004`** — comment. **OK**.
11. **`bin/boltz2_iterate_steering.py:5570`** — argparse help string. **OK**.
12. **`tests/rfdiffusion/test_rfdiffusion.nf:15`** — usage comment (`sbatch tests/...`). **OK**.
13. **`tests/rfdiffusion/run_test_rfdiffusion_plots.slurm.sh:31`** — usage comment. **OK**.
14. **`tests/negative_steering/run_test_negsteer_within_sequence_plots.slurm.sh:49`** — usage comment. **OK**.
15. **`tests/negative_steering/test_negative_steering.nf:41`** — usage comment. **OK**.
16. **`tests/negative_steering/run_test_negsteer_plots.slurm.sh:51`** — usage comment. **OK**.
17. **`tests/haddock/test_haddock.nf:11`** — usage comment. **OK**.
18. **`tests/rosetta_filtering/run_test_rosetta_filtering_plots.slurm.sh:25`** — usage comment. **OK**.
19. **`tests/rosetta_filtering/test_rosetta_filtering.nf:31`** — usage comment. **OK**.
20. **`tests/proteinmpnn/run_test_mpnn_plots.slurm.sh:28`** — usage comment. **OK**.
21. **`tests/orthogonal_metrics/run_test_orthogonal_metrics_plot.slurm.sh:52`** — usage comment. **OK**.
22. **`tests/orthogonal_metrics/test_orthogonal_metrics.nf:35`** — usage comment. **OK**.

The actual `sbatch` invocations from inside `boltz2_iterate_steering.py`
(e.g. lines 1092 onwards) call the bare `sbatch` binary on `$PATH`,
which is the documented convention ("we're OUTSIDE the container here,
so sbatch is on PATH"). No absolute or project-relative paths to slurm
tools anywhere in the tree. No matches for `sacct`, `srun`, `scancel`,
`squeue`, `sinfo`, `sstat`, or `salloc`.

### Action items

No issues found.

---

## PART 4 — Java tool references

### Findings

1. **`run_pipeline.slurm.sh:88`** — `echo "Java:            $(java -version 2>&1 | head -1)"`
   - Bare `java` on `$PATH`. After the JDK fix in Part 1, this resolves via `${JAVA_HOME}/bin/java`. **OK after Part 1 is fixed; currently BROKEN-ADJACENT** (it'll print whatever the system Java is, or fail silently if there is none).
2. **`containers/NextFlow.def:7`** — `export JAVA_HOME=/usr/lib/jvm/java-17-openjdk` — in-container env. **OK**.
3. **`containers/NextFlow.def:13`** — `yum install -y java-17-openjdk-headless ...` — in-container provisioning. **OK**.
4. **`containers/NextFlow.def:16`** — `export JAVA_HOME=/usr/lib/jvm/java-17-openjdk` (`%post`) — **OK**.

No invocations of `javac` or `jar` anywhere. No project-relative `./java` etc. No host-side absolute paths to `java` outside what's covered by the JDK fix.

### Action items

- None beyond Part 1 — fixing `JAVA_HOME` automatically resolves the bare `java` lookup on line 88.

---

## PART 5 — Sanity check on path resolution after the cache-busting commit

### Findings

1. **`tests/*/bin` symlinks (all 6)** — verified pointing to `../../bin`:
   - `tests/haddock/bin -> ../../bin`
   - `tests/negative_steering/bin -> ../../bin`
   - `tests/orthogonal_metrics/bin -> ../../bin`
   - `tests/proteinmpnn/bin -> ../../bin`
   - `tests/rfdiffusion/bin -> ../../bin`
   - `tests/rosetta_filtering/bin -> ../../bin`

2. **`modules/preprocessing.nf` (cache-busted module check)** — `RESOLVE_CONTIGS` declares `path contigs_script` (line 125) and references `${contigs_script}` in its script body (line 133).
   - **Call site in `main.nf:432-436`**:
     ```groovy
     RESOLVE_CONTIGS(
         rfdiff_pdb_ch,
         params.contigs,
         Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
     )
     ```
   - **Call site in `tests/proteinmpnn/test_proteinmpnn.nf:110-114`**:
     ```groovy
     RESOLVE_CONTIGS(
         input_pdb_ch,
         params.contigs,
         Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
     )
     ```
   - Pattern is identical. `${projectDir}` resolves to the repo root in both contexts (the production `projectDir` for `main.nf`, the test workflow's own `projectDir` for the test — but `bin/` exists under both because of the symlink, so `file("${projectDir}/bin/rfdiffusion_contigs.py")` resolves either way).

3. No inconsistencies observed in the spot check.

### Action items

No issues found.

---

## Overall verdict

**Blocker present.** The 7 slurm launchers all set `JAVA_HOME` to the
deleted in-tree path (`${PIPELINE_DIR}/jdk-17.0.2`); fix that before the
HPC round-trip. The `nxf_home` recreation behaviour is non-blocking but
warrants a follow-up decision; the `sbatch` and Part-5 path-resolution
checks are clean.
