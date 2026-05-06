# RFDiffusion Parameter Audit

*Written 2026-05-06. Covers modules/rfdiffusion.nf and the bin/ scripts as of
the experiments branch.*

---

## What the pipeline is doing

The pipeline uses RFDiffusion's **motif-scaffolding / partial diffusion** mode,
not classical binder hallucination. The contig string specifies fixed receptor
segments (e.g. `A1-32`, `A50-68`) that anchor the known scaffold, and design
segments between them (e.g. `10-30`, `8-10`) that are fully diffused from
noise. The effector chain (`B`) is included as a fixed whole-chain entry at the
end of the contig, so RFDiffusion conditions the design on the existing
receptor-effector complex geometry.

This is the correct approach for receptor loop redesign (changing which effector
a receptor recognises). The contig preprocessing pipeline (`rfdiffusion_contigs.py`
→ `contig_utils.resolve_contigs()`) correctly handles user-friendly notation and
converts it to RFDiffusion-native format.

---

## What is missing or unused

### 1. No checkpoint specified — base model used for PPI work

**Status: should be fixed.**

`inference.ckpt_override_path` is never set in `modules/rfdiffusion.nf`. The
container default is almost certainly `Base_ckpt.pt`, which is optimised for
monomer scaffolding. For receptor-effector interface design (a fixed partner
chain is present in the contig), the **Complex_beta checkpoint** is more
appropriate — it was fine-tuned on protein-protein complexes and produces
better interface geometries.

```
inference.ckpt_override_path=/opt/RFdiffusion/models/Complex_base_ckpt.pt
```

Note: the container (`LRR_Pipeline.def`) downloads `Base_ckpt.pt` and
`Complex_base_ckpt.pt` only — `Complex_beta_ckpt.pt` is not present.

### 2. `params.add_potential = true` is dead code

**Status: should be wired up or removed.**

Declared at `main.nf:84` but never passed to the `RFDIFFUSION` process. The
intent was presumably to add guiding potentials (see section below for what
these do). The parameter is ignored at runtime. Either remove it or wire it
into the RFDIFFUSION call with the appropriate potential string.

### 3. Hotspot residues often empty in mode 2

**Status: worth addressing per-campaign.**

In input mode 2 (pre-docked complex), `hotspot: ""` is the default and no
hotspot residues are auto-derived before RFDiffusion runs. The `ppi.hotspot_res`
argument is therefore never passed. This means RFDiffusion gets no explicit
guidance about which effector residues to target with the design region.

For campaigns where the interface is well-characterised (e.g. from a crystal
structure or from HADDOCK), manually specifying key effector residues as
hotspots is likely to improve interface-focused outputs:

```yaml
hotspot: "B45,B67,B89"
```

### 4. No guiding potentials

**Status: easy win once add_potential is wired up.**

RFDiffusion supports energy-like potentials that bias the trajectory during
diffusion without requiring a separate checkpoint. See the detailed section
below for usage guidance.

### 5. `min_hotspot_frac: 0.0` — no enrichment at the RFDiffusion filter stage

**Status: noted; not a priority.**

All designs pass the hotspot-fraction filter regardless of interface quality.
The downstream Boltz2 validation provides the meaningful quality gate. If
compute becomes a constraint, setting this to 0.3–0.5 would enrich the
ProteinMPNN input pool at negligible cost.

### 6. `contigmap.provide_seq` — unused capability

**Status: worth knowing about; apply case-by-case.**

If a position in the design region is known to require a specific amino acid
(a conserved contact residue, a motif anchor), the sequence can be fixed there
while still designing the backbone:

```
contigmap.provide_seq=[A33-35]
```

No current campaign needs this, but it is a useful escape hatch.

---

## Guiding potentials — full description

RFDiffusion's guiding potentials are soft energy terms added to the score at
each denoising step. They do not change the model weights — they steer the
trajectory by nudging the gradient toward configurations with lower potential
energy, analogous to a weak force applied throughout diffusion. The key
parameters are:

```
potentials.guiding_potentials=["type:<name>,weight:<w>[,<extra_args>]"]
potentials.guide_scale=<scale>      # global multiplier (default 1)
potentials.guide_decay=<decay>      # how weight decays over timesteps
                                    # options: "constant", "linear",
                                    #          "quadratic", "cubic"
```

`guide_decay` controls when the potential is strongest. `"quadratic"` (the
recommended default) concentrates the guidance early in diffusion when global
topology is being set, then fades as fine structure is locked in.

### Available potential types

**`interface_ncontacts`** — promotes contacts between the design chain and
a fixed partner chain. This is the most directly useful potential for the
receptor-resurfacing task. Weight 1–2 is a reasonable starting point.

```
potentials.guiding_potentials=["type:interface_ncontacts,weight:1"]
```

**`monomer_ROG`** — penalises high radius of gyration, encouraging compact
designs. Prevents diffuse, extended loops. Useful when design segments are
long (>15 residues) and there is no strong steric constraint from the
fixed scaffold to keep them compact.

```
potentials.guiding_potentials=["type:monomer_ROG,weight:1,min_dist:5"]
```

`min_dist` sets a lower bound on the allowed Rg (in Å); values below this are
not penalised (prevents over-collapsing short loops).

**`olig_contacts`** — for oligomers; not relevant here.

### Combining potentials

Multiple potentials are passed as a list:

```
potentials.guiding_potentials=["type:interface_ncontacts,weight:1","type:monomer_ROG,weight:0.5,min_dist:5"]
```

### Practical notes

- Potentials add compute per step (extra gradient pass) — expect ~10–20 %
  slower wall time.
- Too-high weights can destabilise diffusion and produce clashing or
  geometrically odd structures. Start with weight 1 and increase cautiously.
- `guide_decay="quadratic"` is almost always better than `"constant"` —
  constant-weight potentials tend to over-constrain the late fine-structure
  steps.
- The most useful combination for this pipeline is likely
  `interface_ncontacts` (to drive interface contacts) plus optionally
  `monomer_ROG` for longer design regions.

---

## Changes made (2026-05-06)

### `params.rfdiff_checkpoint` wired through the pipeline

The following files were changed to expose checkpoint selection as a
first-class parameter:

- **`main.nf`**: Added `params.rfdiff_checkpoint = "Complex_base_ckpt.pt"`
  in the RFDiffusion params block. Threaded it as a new positional argument
  to the `RFDIFFUSION(...)` process call.

- **`modules/rfdiffusion.nf`**: Added `val checkpoint` to the `RFDIFFUSION`
  process input block. Built `def ckpt_arg` in the script section:
  ```
  def ckpt_arg = checkpoint ? "inference.ckpt_override_path=/opt/RFdiffusion/models/${checkpoint}" : ""
  ```
  Appended `${ckpt_arg}` to the `run_inference.py` argument list (between
  the hotspot arg and `diffuser.T`). An empty string is passed through cleanly
  as a no-op, preserving backward compatibility.

- **`params_example.yml`**: Added `rfdiff_checkpoint: "Complex_base_ckpt.pt"`
  with a comment listing the three available options
  (`Base_ckpt.pt`, `Complex_base_ckpt.pt`, `Complex_beta_ckpt.pt`).

- **`tests/full_test_run/params_full_test.yml`**: Added
  `rfdiff_checkpoint: "Complex_base_ckpt.pt"`.

The default is `Complex_beta_ckpt.pt` — the PPI-optimised checkpoint now
present in the rebuilt container (`LRR_Pipeline.img`). `Complex_base_ckpt.pt`
remains available as a fallback.

Campaign `params.yml` files under `experiments/` are intentionally NOT updated
here — they will pick up the `main.nf` default until explicitly overridden.

### Guiding potentials wired through the pipeline

The following files were changed to expose `interface_ncontacts` and
`monomer_ROG` potential guidance as first-class parameters:

- **`main.nf`**: Replaced dead `params.add_potential = true` with six
  explicit params (default `add_potential = false` so existing campaigns
  are unaffected):
  ```
  params.add_potential           = false
  params.rfdiff_guide_scale      = 2
  params.rfdiff_guide_decay      = "quadratic"
  params.rfdiff_interface_weight = 1.0
  params.rfdiff_rog_weight       = 0.5
  params.rfdiff_rog_min_dist     = 5
  ```
  All six threaded into the `RFDIFFUSION(...)` process call.

- **`modules/rfdiffusion.nf`**: Added six `val` inputs. In the script block,
  builds `potential_arg` from the inputs — each potential is included only
  if its weight > 0, so individual potentials can be disabled by setting
  their weight to 0:
  ```groovy
  def potential_arg = ""
  if (add_potential) {
      def pot_list = []
      if ((interface_weight as Double) > 0) pot_list << "type:interface_ncontacts,weight:${interface_weight}"
      if ((rog_weight as Double) > 0)       pot_list << "type:monomer_ROG,weight:${rog_weight},min_dist:${rog_min_dist}"
      if (pot_list) {
          def pot_str = pot_list.collect { "'${it}'" }.join(",")
          potential_arg = "\"potentials.guiding_potentials=[${pot_str}]\" potentials.guide_scale=${guide_scale} potentials.guide_decay=${guide_decay}"
      }
  }
  ```

- **`params_example.yml`**: All six params added with full inline documentation.
  Default `add_potential: true` in this file — new campaigns copying the
  example will have potentials enabled.

- **`tests/full_test_run/params_full_test.yml`**: All six params added,
  `add_potential: true`.

---

## RFDiffusion version landscape (as of May 2026)

### Version summary

| Version | Released | Primary use case | Relevant here? |
|---|---|---|---|
| **RFDiffusion v1** | 2023 | Motif scaffolding, binder design, partial diffusion | Yes — what the pipeline uses |
| **RFDiffusion2** | Apr 2025 | Enzyme active site scaffolding (atom-level motifs, rotamer inference) | No — different problem |
| **RFDiffusion3** | Dec 2025 | All-atom design (protein-protein, protein-DNA, small molecules, enzymes) | Not yet — too new |

### RFDiffusion v1 (what we use)

Still the most experimentally validated and practically deployable tool for
protein-protein interface design and receptor loop redesign. Has not been
superseded for this task. The canonical repo is `RosettaCommons/RFdiffusion`
(the authoritative source); `sokrypton/RFdiffusion` is a Colab-friendly
mirror with identical model weights, not an independent development line.
For HPC/Nextflow use, the RosettaCommons repo is the right reference.

Our container clones from `sokrypton/RFdiffusion` — the weights are the same
but it is worth knowing the RosettaCommons repo is where issues and upstream
fixes land.

### RFDiffusion2 (not relevant to this pipeline)

Released April 2025 (*Nature Methods* 2025). Specialised for **enzyme active
site scaffolding** — it accepts unindexed atomic motifs (functional groups
without pre-specified residue indices) and infers both rotamers and the
sequence position of the motif. Uses flow matching rather than DDPM.
Repo: `RosettaCommons/RFdiffusion2`. Not aimed at protein-protein interface
redesign.

### RFDiffusion3 (watch for late 2026)

Released December 2025 (bioRxiv Sep 2025; code in `RosettaCommons/foundry`).
Fundamental shift: **individual atoms (backbone + side chain) are the
diffusion units**, not residue frames. Benchmarks show matched or superior
performance to v1 on protein-protein binding tasks, and ~10× faster batch
throughput than v2. Theoretically the right long-term successor for this
pipeline, but as of May 2026 it is only ~5 months old. Documentation and
community workflows for receptor loop redesign are not yet mature. Revisit
in late 2026.

### β-pairing targeted RFDiffusion fine-tune (potentially relevant)

Published *Nature Communications* 2025. A fine-tuned RFDiffusion checkpoint
that biases diffusion toward β-strand pairing, improving binder design
success rates at β-strand-mediated interfaces. Given that the HMA/MAX effector
system involves a β-sandwich fold, this fine-tune may be worth testing once
the main pipeline is validated. The checkpoint and weights should be available
via the Baker lab / IPD.

---

## Remaining recommended changes

### 1. Guiding potentials — **done (2026-05-06)**

`params.add_potential` and associated weight/scale/decay params wired through
the full pipeline. See "Changes made" section for details.

### 2. Hotspot derivation in mode 2 — **higher effort**

The lightweight fix is to document that `hotspot:` should always be set
manually in mode 2 campaign params when the interface is known. No code
change required.

The proper fix — auto-deriving hotspots from the input complex before
RFDiffusion runs — requires a new pipeline process (contact extraction from
the pre-docked PDB, analogous to the existing `EXTRACT_HOTSPOTS` that runs
post-HADDOCK in mode 1). That is a meaningful addition: new bin script, new
Nextflow process, branching logic in `main.nf`.

### 3. `contigmap.provide_seq` — **no code change needed**

Already supported by RFDiffusion. If a specific amino acid is required at a
position in the design region (conserved contact residue, motif anchor), pass
it directly in the `run_inference.py` call:

```
contigmap.provide_seq=[A33-35]
```

Apply per-campaign as needed; no pipeline changes required.
