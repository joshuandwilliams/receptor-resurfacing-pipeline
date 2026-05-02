# pipeline_notes5.md — cuEquivariance kernel fix (boltz2_negsteer container)

## Goal

Enable cuEquivariance CUDA kernels inside the `boltz2_negsteer` Singularity/Apptainer container so Boltz-2 inference runs 20–30% faster than the CPU-fallback path. Relevant because production experiments run ~10,000 Boltz predictions; compounded, the speedup is meaningful.

## Symptom

With the pre-existing container, any attempt to `import cuequivariance_ops_torch` inside Python on a GPU node failed with:

```
ImportError: libcue_ops.so: cannot open shared object file: No such file or directory
```

A closer look at the preceding line:

```
Error while loading libcue_ops.so: /usr/local/lib/python3.11/site-packages/cuequivariance_ops/lib/libcue_ops.so: undefined symbol: cublasGemmGroupedBatchedEx, version libcublas.so.12
```

The workaround was to invoke the pipeline with `--negsteer_no_kernels true`, which skips kernel-based ops and routes through a slower Python/PyTorch path. Worked, but left 20–30% on the table.

## Root cause (established during debugging)

Three independent facts combined:

1. **cuEquivariance 0.9.1's `libcue_ops.so` requires `cublasGemmGroupedBatchedEx`**, a symbol introduced in cuBLAS 12.5. Confirmed via `nm -D` on the `.so`.

2. **torch 2.6.0+cu124 ships cuBLAS 12.4.5.8** bundled as the `nvidia-cublas-cu12` pip wheel. This wheel does NOT contain the required symbol.

3. **`torch/__init__.py`'s `_preload_cuda_deps()` dlopens the bundled cuBLAS at import time**, populating process memory with the stale 12.4 library. When `libcue_ops.so` is subsequently loaded and asks the dynamic linker for `libcublas.so.12`, Linux's policy of one SONAME per process means the already-loaded 12.4 version is reused — LD_LIBRARY_PATH changes at that point are inert.

So even if we installed cuBLAS 12.8 somewhere on disk, torch would always win at runtime by loading 12.4 first.

## Approaches that did NOT work (and why)

These were tried, confirmed broken, and abandoned:

- **`pip install --upgrade nvidia-cublas-cu12==12.6.*`** — succeeded in isolation but produced `ResolutionImpossible` errors on every subsequent pip install because torch 2.6 declares a hard `==12.4.5.8` pin on `nvidia-cublas-cu12`.
- **Pinning upgraded NVIDIA wheels in the constraints file** — same `ResolutionImpossible` cascade.
- **`LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH` in `%environment`** — inert. Once torch's `_preload_cuda_deps` has loaded the pip-bundled cuBLAS into the process, the dynamic linker won't look elsewhere. Proven with the `ctypes.CDLL` vs `import torch; ctypes.CDLL` comparison during sandbox testing: the raw dlopen succeeded against system 12.8, but the same dlopen after `import torch` failed with the undefined-symbol error.
- **Overriding `LD_LIBRARY_PATH` at `apptainer exec` invocation** — same reason. Too late; torch has already loaded stale cuBLAS by the time the user's code runs.

## The architectural fix: NGC-style build

Pattern lifted from NVIDIA's NGC PyTorch containers (e.g. 25.01): install CUDA at the **system level** via the NVIDIA rhel9 yum repo, install torch with `--no-deps` so it doesn't pull the pip wheels, rely on `ldconfig` + system `/usr/local/cuda/lib64` to satisfy torch's runtime dlopens.

Concretely in `%post`:

1. Add NVIDIA rhel9 yum repo (`cuda-rhel9.repo`).
2. `dnf install` the minimal CUDA 12.8 runtime packages (package names verified against NVIDIA's CUDA-Q build script, cuSPARSELt download page, and NCCL download page — naming is inconsistent across libraries, not guessable).
3. `pip install --no-deps torch==2.6.0` to prevent torch pulling its declared nvidia-*-cu12 wheels.
4. Install torch's non-CUDA python deps explicitly (filelock, sympy, etc.).
5. Install the rest of the stack (cuequivariance, trifast, boltz, DockQ, freesasa, MDAnalysis) with constraints pinning versions.

**But this alone was insufficient.** Subsequent pip installs re-resolved torch's declared deps and pulled in all 14 nvidia-*-cu12 wheels anyway — the `--no-deps` only applied to the literal `pip install torch` command. Once those wheels exist on disk, torch's `_preload_cuda_deps` finds them via entry-point metadata and loads them in preference to the system libs.

**So the final step that actually fixes it:** after all installs complete, purge the 12 stale nvidia-*-cu12 wheels from the image.

### What the purge does

```bash
pip uninstall -y \
    nvidia-cublas-cu12 nvidia-cuda-cupti-cu12 nvidia-cuda-nvrtc-cu12 \
    nvidia-cuda-runtime-cu12 nvidia-cudnn-cu12 nvidia-cufft-cu12 \
    nvidia-curand-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12 \
    nvidia-nccl-cu12 nvidia-nvjitlink-cu12 nvidia-nvtx-cu12

rm -rf /usr/local/lib/python3.11/site-packages/nvidia/{cublas,cuda_cupti,cuda_nvrtc,cuda_runtime,cudnn,cufft,curand,cusolver,cusparse,nccl,nvjitlink,nvtx}
```

Two wheels kept intentionally:
- `nvidia-cusparselt-cu12==0.6.2` — torch wants this exact version for `libcusparseLt.so.0` at runtime. The system repo's `libcusparselt0` is a newer version. Keeping the pip wheel satisfies torch's version-sensitive dlopen.
- `nvidia-ml-py` — not a `-cu12` wheel. Provides the `pynvml` import namespace.

### Verification step baked in

After the purge, `%post` enumerates the 12 forbidden wheels explicitly (not via regex — that bit us once) and fails the build if any are still present:

```bash
STALE_WHEELS=$(pip list 2>/dev/null | awk '{print $1}' \
    | grep -E "^nvidia-(cublas|cuda-cupti|cuda-nvrtc|cuda-runtime|cudnn|cufft|curand|cusolver|cusparse|nccl|nvjitlink|nvtx)-cu12$" || true)
if [ -n "$STALE_WHEELS" ]; then
    echo "FATAL: stale nvidia-*-cu12 wheels still present:"; echo "$STALE_WHEELS"
    exit 1
fi
```

Then re-verifies `import torch` still works (it does, because torch was installed `--no-deps` so pip doesn't record these as required).

## The pynvml side-fix

torch 2.6's `torch/cuda/__init__.py` does `import pynvml`. Two PyPI packages provide that module name:
- `pynvml` — original third-party wrapper, deprecated.
- `nvidia-ml-py` — NVIDIA's own maintained package, also exporting `pynvml`.

With both installed, torch finds the deprecated one first and emits `FutureWarning` on every `import torch`. Uninstalling the `pynvml` package leaves `nvidia-ml-py` as the sole provider; `import pynvml` resolves to its shim (`/usr/local/lib/python3.11/site-packages/pynvml.py`), and the warning disappears.

Added to the same purge block in `%post`.

## Validation path (writable sandbox before bake-in)

Container rebuilds take ~hour. To avoid another failed rebuild, the fix was validated on a writable sandbox derived from the existing broken image:

1. `apptainer build --sandbox boltz2_sandbox <broken.img>` — unpacks the squashfs into a mutable directory tree.
2. `apptainer exec --writable --fakeroot boltz2_sandbox pip uninstall -y <the 12 wheels>`.
3. `rm -rf` the leftover `site-packages/nvidia/*` dirs.
4. `apptainer exec --writable --fakeroot boltz2_sandbox pip uninstall -y pynvml`.
5. Submit GPU srun and verify:
   - `import torch` OK, `torch.cuda.is_available()` True, no pynvml warning
   - `import cuequivariance_ops_torch` OK (the critical test — this was the failing import)
   - `import trifast` OK
   - A real cuBLAS matmul on GPU returns a valid result

All six criteria passed in the sandbox. Then and only then was the .def updated and a full rebuild kicked off. The rebuilt image passed the same six checks.

Notable sandbox logistics:
- Built on NBI-HPC (software node) where `/tmp/jowillia` exists and sandbox builds work. GPU nodes can't see `/tmp` across nodes, so sandbox had to be copied to `/hpc-home/jowillia/singularity/`.
- Per-file `cp` across NFS filesystems with 512 KB block size is extremely slow (135k files → ~45 min copy, ~50 GB destination due to block-rounding on a 21 GB source).
- Lesson for next time: build the sandbox directly at its final destination via `apptainer build --sandbox /final/path <img>`. Squashfs decompression is sequential and ~10× faster than cp'ing the unpacked tree.

## Current state (post-fix)

Container location: `/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/boltz2_negsteer.img`

Smoke test on rebuilt image passes end-to-end:
```
torch: 2.6.0+cu124 cuda: True
pynvml from: /usr/local/lib/python3.11/site-packages/pynvml.py
ceot OK
trifast import OK, attrs: ['autotune', 'autotune_helpers', 'torch', 'triangle_attention', 'triton']
cuBLAS matmul OK, shape: torch.Size([8, 64, 64]) sum: 31220.7265625
```

Next: real negative-steering test with kernels enabled (drop `--negsteer_no_kernels true`), compare wall-clock to the v3 run that used the CPU-fallback path. **Result pending as of writing**; will be appended when available.

## Key design principles learned

- **Don't trust pip's `--no-deps`** when the downstream installs will re-resolve the dep tree. Either uninstall after, or use `uv pip install --no-deps --system ...` and build a lockfile.
- **Python dynamic linking is one-shot per SONAME per process.** Once a `.so` is loaded, LD_LIBRARY_PATH changes can't redirect subsequent dlopens of the same SONAME. This is why the NGC approach requires physically removing the competing libs, not just re-ordering paths.
- **Always enumerate package names explicitly in verification regexes.** Lazy wildcards like `^nvidia-.*-cu12` match things you didn't mean to match. The verification step is specifically the one place where laziness is most expensive — when it fires, the build fails at the very end, after hours.
- **Sandbox-test architectural changes before rebuilding from scratch.** Even a writable sandbox copy that takes 45 min to set up is cheap compared to a 1-hour rebuild that fails at the end.
- **Verify package names against authoritative sources, not pattern-matching.** CUDA package names on NVIDIA's rhel9 repo are inconsistent: `cuda-cudart-12-8`, `libcublas-12-8`, `libcudnn9-cuda-12`, `libcusparselt0`, `libnccl`. Each library was verified from NVIDIA's own docs (CUDA-Q, cuSPARSELt, NCCL download pages) before adding to the dnf list.

## Files touched

- `/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/boltz2_negsteer.def` — NGC-style rebuild with system CUDA 12.8 + pip wheel purge + pynvml uninstall
- `/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/boltz2_negsteer.img` — rebuilt from above .def
