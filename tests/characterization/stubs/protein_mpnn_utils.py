"""A stand-in for the three ProteinMPNN symbols bin/mpnn_design_region_score.py imports.

`mpnn_design_region_score.py` runs inside the RFDiffusion container, where
/opt/ProteinMPNN is on sys.path and `parse_PDB`, `tied_featurize` and
`ProteinMPNN` come from the real `protein_mpnn_utils`. Off the cluster that
import fails, the module sets HAS_MPNN = False, and `main()` exits 1 before
reaching any of the scoring arithmetic — which is the only part of the module
worth testing.

This file supplies those three names so the module can be imported with
HAS_MPNN = True. It is not on sys.path by default: a test puts its directory
there for the duration of one `importlib` load (see
`fixtures_mpnn_stage.mpnn_module_with_stub`), so nothing else in the suite can
pick it up by accident.

The point is not to imitate the network. It is to make the scorer's arithmetic
predictable, so a test can state the design-region score before the code runs:

`parse_PDB` really reads the PDB, so the chain ids and per-chain sequence
lengths come from the fixture structure rather than from a hardcoded shape.

`tied_featurize` derives the free/fixed mask the way the real one does, from
1-based positions within each chain, so a test states which residues are fixed
and the free-position count follows.

`ProteinMPNN.__call__` returns log-probabilities that depend only on the amino
acid, not on the position or the coordinates: log P(aa) = -(index(aa) + 1) in
ProteinMPNN's own alphabet. The scorer's per-position score is
`-log_probs[i, S[i]]`, so each position contributes exactly
`index(designed_aa) + 1` — 1 for alanine, 3 for cysteine, 21 for the X that
padding inserts. A mean over a stated set of free positions is then
arithmetic, and a bug in the sequence-to-index conversion changes the answer
rather than being averaged away.
"""
from __future__ import annotations

import os

import torch

# ProteinMPNN's own residue alphabet, index 20 being the unknown/padding X.
# mpnn_design_region_score.py carries the same string; the stub repeats it
# because the real protein_mpnn_utils does too.
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"

# Backbone atoms per residue in the featurised coordinate tensor (N, CA, C, O).
N_BACKBONE_ATOMS = 4

# The real tied_featurize returns a 20-tuple. The scorer indexes 0, 1, 2, 4,
# 5, 10 and 12 out of it and ignores the rest, so the rest stay None: a test
# that starts depending on one of them fails loudly instead of reading a
# plausible-looking zero.
N_FEATURIZE_OUTPUTS = 20

_THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def parse_PDB(pdb_path, **kwargs):  # noqa: N802 - name fixed by the real module
    """One entry per PDB, carrying `name` and a seq_chain_X per chain.

    Only CA records are read, which is all the fixture structures from
    `fixtures_runtree.complex_pdb` contain. Chains keep PDB file order, and
    `name` is the basename without its extension, because that is the key
    `mpnn_design_region_score.py` looks the fixed-position dict up under.
    """
    chains: dict[str, list[str]] = {}
    order: list[str] = []
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21]
            if chain not in chains:
                chains[chain] = []
                order.append(chain)
            chains[chain].append(_THREE_TO_ONE.get(line[17:20].strip(), "X"))

    entry: dict = {"name": os.path.basename(pdb_path).rsplit(".", 1)[0],
                   "num_of_chains": len(order)}
    for chain in order:
        entry[f"seq_chain_{chain}"] = "".join(chains[chain])
        entry[f"coords_chain_{chain}"] = {
            f"{atom}_chain_{chain}": [[0.0, 0.0, 0.0]] * len(chains[chain])
            for atom in ("N", "CA", "C", "O")
        }
    entry["seq"] = "".join(entry[f"seq_chain_{c}"] for c in order)
    return [entry]


def tied_featurize(pdb_dict_list, device, chain_dict, fixed_position_dict=None,
                   **kwargs):
    """The 20-tuple, with the seven entries the scorer reads filled in.

    Chains are concatenated in sorted id order, matching the `chain_ids.sort()`
    the scorer does before building its chain_dict, so position i of the
    returned tensors is position i of the scorer's view.

    `chain_M_pos` (index 10) is the free/fixed mask: 0 for a position listed in
    `fixed_position_dict`, 1 otherwise. Those listed positions are 1-based
    indices within their own chain, which is how ProteinMPNN's own driver
    writes the fixed_positions JSONL that `load_fixed_positions` reads back.
    """
    entry = pdb_dict_list[0]
    chain_ids = sorted(k[len("seq_chain_"):] for k in entry
                       if k.startswith("seq_chain_"))
    seq = "".join(entry[f"seq_chain_{c}"] for c in chain_ids)
    length = len(seq)

    fixed_for_name = (fixed_position_dict or {}).get(entry["name"]) or {}
    free: list[float] = []
    chain_index: list[float] = []
    for number, chain in enumerate(chain_ids, start=1):
        fixed_here = set(fixed_for_name.get(chain) or [])
        for position in range(1, len(entry[f"seq_chain_{chain}"]) + 1):
            free.append(0.0 if position in fixed_here else 1.0)
            chain_index.append(float(number))

    out: list = [None] * N_FEATURIZE_OUTPUTS
    out[0] = torch.zeros(1, length, N_BACKBONE_ATOMS, 3, device=device)
    out[1] = torch.tensor(
        [[ALPHABET.index(aa) if aa in ALPHABET else 20 for aa in seq]],
        dtype=torch.long, device=device)
    out[2] = torch.ones(1, length, device=device)
    out[4] = torch.ones(1, length, device=device)
    out[5] = torch.tensor([chain_index], dtype=torch.float32, device=device)
    out[10] = torch.tensor([free], dtype=torch.float32, device=device)
    out[12] = torch.arange(length, device=device).unsqueeze(0)
    return tuple(out)


class ProteinMPNN:
    """log P(aa) = -(index(aa) + 1), independent of position and structure.

    The constructor keyword arguments, the loaded state dict and the `eval()`
    call are all recorded rather than used, so a test can assert that
    `load_model` read `noise_level` and `num_edges` out of the checkpoint and
    put the model in eval mode.
    """

    def __init__(self, **kwargs):
        self.init_kwargs = dict(kwargs)
        self.loaded_state_dict = None
        self.device = None
        self.is_eval = False

    def load_state_dict(self, state_dict):
        self.loaded_state_dict = state_dict

    def to(self, device):
        self.device = device
        return self

    def eval(self):
        self.is_eval = True
        return self

    def __call__(self, X, S, mask, chain_M, residue_idx, chain_encoding, randn):
        length = S.shape[1]
        row = [-(index + 1.0) for index in range(len(ALPHABET))]
        return torch.tensor([[row] * length], dtype=torch.float32)
