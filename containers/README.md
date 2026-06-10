# Containers

This directory holds the Singularity definition (`.def`) files that build the
container images used by the pipeline. The `.def` files are version-controlled;
the built images themselves are **not** — they are large binary build artefacts
(multi-GB) and are rebuilt from the `.def` files when needed.

## How the pipeline finds its images

Every component that runs a container references it **by a fixed filename
inside this `containers/` directory** — never by an absolute path on any one
machine. `nextflow.config` uses `${projectDir}/containers/<name>.img`, the
SLURM launchers use `${PIPELINE_DIR}/containers/<name>.img`, and the JDK is
`${PIPELINE_DIR}/containers/jdk-17.0.2`.

The actual `.img` / `.sif` files (and the JDK) are **not** committed. On a given
machine, each expected filename is a **symlink** into wherever the real images
live (on the NBI HPC software node that is
`/hpc-home/<user>/singularity/...`). These symlinks are git-ignored, so they
never reach GitHub and never break another user's checkout.

**To install the pipeline:** drop your own built images (or symlinks to them)
into this directory using the exact filenames below. For example, on the NBI
HPC:

```bash
cd containers/
ln -s /hpc-home/$USER/singularity/NextFlow/NextFlow.img                                       NextFlow.img
ln -s /hpc-home/$USER/singularity/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img  HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img
ln -s /hpc-home/$USER/singularity/Rosetta/Rosetta.img                                         Rosetta.img
ln -s /hpc-home/$USER/singularity/Boltz1_Boltz2_Chai1_ColabFold/boltz2_negsteer.img           boltz2_negsteer.img
ln -s /hpc-home/$USER/singularity/ColabFold/colabfold.img                                     colabfold.img
ln -s /hpc-home/$USER/singularity/pytest/pytest_runner.img                                    pytest_runner.img
ln -s /hpc-home/$USER/singularity/jdk-17.0.2                                                  jdk-17.0.2
```

Expected filenames (must match exactly):

| Symlink / file in `containers/`              | Provides                          |
| -------------------------------------------- | --------------------------------- |
| `NextFlow.img`                               | `params` — Nextflow runner image  |
| `HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img`| `params.rfdiff_container`          |
| `Rosetta.img`                                | `params.rosetta_container`        |
| `boltz2_negsteer.img`                        | `params.boltz2_container`         |
| `colabfold.img`                              | `params.colabfold_container`      |
| `pytest_runner.img`                          | characterisation test suite       |
| `jdk-17.0.2`                                 | OpenJDK 17 (`JAVA_HOME` for Nextflow) |

Run `scripts/check_environment.slurm.sh` (an sbatch job on the GPU queue) to
confirm every symlink resolves and every tool inside the images is reachable.

> Note on file extension: NBI HPC uses Singularity 3.8 with the older `.img`
> (squashfs) format; modern Singularity / Apptainer produces `.sif`. The
> definition files are format-agnostic — choose the extension to match your
> Singularity version when building. If you build `.sif` images, either name
> the symlinks `*.img` anyway or update the `params.*_container` filenames.

## Image inventory

| `.def` file              | Pipeline stage(s) / module(s) consuming it                                                                                                                                                                                                                                                                                                                                                                                       | Python | Contents                                                                                                                                                                                                                          |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NextFlow.def`           | Workflow runner (`run_pipeline.slurm.sh` and every `tests/*/run_test_*.slurm.sh` invokes `nextflow` from this image; modules themselves do not reference it).                                                                                                                                                                                                                                                                   | n/a    | Java 17 (OpenJDK headless) + Nextflow 25.10.4 dist binary, with `nf-tower` plugin pre-cached and `NXF_OFFLINE=true`.                                                                                                              |
| `HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.def` | `params.rfdiff_container` — used by `rfdiffusion.nf`, `proteinmpnn.nf`, `preprocessing.nf`, `pose_solver.nf`, the plotting steps in `negative_steering.nf` and `negsteer_orthogonal_metrics.nf`.                                                                                                                                                                                                                                       | 3.10   | PyTorch 2.4.0 + CUDA 12.4, RFDiffusion (with model weights + schedules), ProteinMPNN, ColabDesign v1.1.1, DGL, e3nn 0.5.5, HADDOCK3 (bundled CNS binary) + OpenMM 8.2 + pdbfixer, MMseqs2; numpy / scipy / pandas / biopython / py3Dmol. |
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
sudo singularity build HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img  HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.def
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
- `HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.def` downloads RFDiffusion weights
  and schedules into `/opt/RFdiffusion/`; HADDOCK3 ships its own CNS binary via
  the PyPI wheel.
