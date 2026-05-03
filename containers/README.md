# Containers

This directory holds the Singularity definition (`.def`) files that build the
container images used by the pipeline. The corresponding image files
(`.img` / `.sif`) live on the HPC software node under
`/hpc-home/jowillia/singularity/` (the exact paths are wired up in
`nextflow.config` via `params.rfdiff_container`, `params.rosetta_container`,
`params.boltz2_container`, and in `run_pipeline.slurm.sh` via `NEXTFLOW_IMG`).
The `.def` files are version-controlled in this repository; the built images
themselves are not — they are large binary build artefacts (multi-GB) and are
rebuilt from the `.def` files when needed.

> Note on file extension: NBI HPC uses Singularity 3.8 with the older `.img`
> (squashfs) format; modern Singularity / Apptainer produces `.sif`. The
> definition files are format-agnostic — choose the extension to match your
> Singularity version when building.

## Image inventory

| `.def` file              | Pipeline stage(s) / module(s) consuming it                                                                                                                                                                                                                                                                                                                                                                                       | Python | Contents                                                                                                                                                                                                                          |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NextFlow.def`           | Workflow runner (`run_pipeline.slurm.sh` and every `tests/*/run_test_*.slurm.sh` invokes `nextflow` from this image; modules themselves do not reference it).                                                                                                                                                                                                                                                                   | n/a    | Java 17 (OpenJDK headless) + Nextflow 25.10.4 dist binary, with `nf-tower` plugin pre-cached and `NXF_OFFLINE=true`.                                                                                                              |
| `LRR_Pipeline.def`       | `params.rfdiff_container` — used by `rfdiffusion.nf`, `proteinmpnn.nf`, `preprocessing.nf`, `haddock.nf`, the plotting steps in `negative_steering.nf` and `negsteer_orthogonal_metrics.nf`.                                                                                                                                                                                                                                       | 3.10   | PyTorch 2.4.0 + CUDA 12.4, RFDiffusion (with model weights + schedules), ProteinMPNN, ColabDesign v1.1.1, DGL, e3nn 0.5.5, HADDOCK3 (bundled CNS binary) + OpenMM 8.2 + pdbfixer, MMseqs2; numpy / scipy / pandas / biopython / py3Dmol. |
| `Rosetta.def`            | `params.rosetta_container` — used by `rosetta_filtering.nf` and `negsteer_rosetta_metrics.nf`.                                                                                                                                                                                                                                                                                                                                  | 3.12   | Rosetta 3.15 source build, release mode, four apps only: `rosetta_scripts`, `score_jd2`, `InterfaceAnalyzer`, `relax`. Plus numpy / scipy / pandas / biopython / matplotlib for analysis.                                         |
| `boltz2_negsteer.def`    | `params.boltz2_container` — used by `negative_steering.nf`, `negsteer_biophysical_metrics.nf`, `negsteer_orthogonal_metrics.nf` (DockQ step), `negsteer_interface_metrics.nf`, `negsteer_controls.nf`, `negsteer_manifest.nf`, and the post-processing step in `negsteer_af3_nomsa.nf`.                                                                                                                                            | 3.11   | Boltz-2 (with pre-downloaded weights: `boltz2_conf.ckpt`, `boltz2_aff.ckpt`, `mols/`), PyTorch 2.6.0+cu124 against system-level CUDA 12.8 runtime, cuEquivariance 0.9.1 (cu12) kernels, trifast 0.1.13, DockQ, freesasa, MDAnalysis, gemmi, scikit-learn. |
| `pytest_runner.def`      | Characterisation test suite under `tests/characterization/` (every `local_unit`, `local_integration`, and `hpc`-tier test invoked via `pytest`); not referenced from `nextflow.config`.                                                                                                                                                                                                                                          | 3.11   | pytest, pandas, numpy, Pillow (`PIL`), scikit-image (`skimage`). No PyTorch, no CUDA, no pipeline runtime — the comparators only diff files produced by upstream stages. Package set matches the `[project.optional-dependencies]` "test" extra in `pyproject.toml`. |

## Building

Build each image on the HPC software node (where Singularity is installed and
fakeroot / `sudo` is available — compute nodes typically lack the privileges
for `singularity build`). The canonical command is:

```bash
sudo singularity build <name>.img <name>.def
```

For example:

```bash
sudo singularity build LRR_Pipeline.img       LRR_Pipeline.def
sudo singularity build Rosetta.img            Rosetta.def
sudo singularity build boltz2_negsteer.img    boltz2_negsteer.def
sudo singularity build NextFlow.img           NextFlow.def
sudo singularity build pytest_runner.img      pytest_runner.def
```

Substitute `.sif` for `.img` if your Singularity / Apptainer version produces
SIF format. After building, copy or symlink the image into the path expected
by `nextflow.config` (or override the corresponding `params.*_container`
variable on the command line).

Build-time notes worth knowing:

- `boltz2_negsteer.def` performs pre-flight URL reachability checks before the
  ~30-minute install, downloads Boltz-2 weights into `/opt/boltz_cache` so
  compute nodes can run air-gapped, and aborts the build on any version-pin
  drift (torch 2.6.0+cu124, numpy <2, pandas 2.x). GPU-gated import checks
  (`trifast`, `cuequivariance_ops_torch`) are reported as `INFO` on a CPU-only
  build node and only resolve at runtime on a GPU node.
- `Rosetta.def` compiles Rosetta 3.15 from source (~5 GB download, then a
  multi-hour `scons -j8` build). The `bin/` entries are symlinks into
  `build/`, so do **not** delete `build/` post-compile.
- `LRR_Pipeline.def` downloads RFDiffusion weights and schedules into
  `/opt/RFdiffusion/`; HADDOCK3 ships its own CNS binary via the PyPI wheel.
