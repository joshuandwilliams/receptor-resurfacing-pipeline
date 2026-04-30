#!/usr/bin/env python3
"""
mpnn_design_region_score.py
---------------------------
Compute design-region-specific MPNN scores by running ProteinMPNN's
model on each design PDB with each designed sequence.

Uses ProteinMPNN's tied_featurize() and model forward pass to get
per-position log-probabilities, then computes the mean over free
(designable) positions only.

Must run inside the ProteinMPNN container (/opt/ProteinMPNN/).

Outputs:
    Updated metadata CSV with 'design_region_score' column.
"""

import argparse
import csv
import json
import os
import sys
import traceback

MPNN_DIR = "/opt/ProteinMPNN"
if os.path.isdir(MPNN_DIR):
    sys.path.insert(0, MPNN_DIR)

import numpy as np

try:
    import torch
    from protein_mpnn_utils import ProteinMPNN, parse_PDB, tied_featurize
    HAS_MPNN = True
except ImportError as e:
    print(f"WARNING: Cannot import ProteinMPNN: {e}", file=sys.stderr)
    HAS_MPNN = False


ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
AA_TO_IDX = {aa: i for i, aa in enumerate(ALPHABET)}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-csv", required=True)
    parser.add_argument("--pdb-dir", required=True)
    parser.add_argument("--mpnn-fasta-dir", required=True)
    parser.add_argument("--fixed-positions-dir", required=True)
    parser.add_argument("--checkpoint", default="")
    parser.add_argument("--output-csv", required=True)
    return parser.parse_args()


def find_checkpoint():
    weight_dir = os.path.join(MPNN_DIR, "vanilla_model_weights")
    if os.path.isdir(weight_dir):
        for f in sorted(os.listdir(weight_dir)):
            if f.endswith(".pt") and "020" in f:
                return os.path.join(weight_dir, f)
        for f in sorted(os.listdir(weight_dir)):
            if f.endswith(".pt"):
                return os.path.join(weight_dir, f)
    return None


def load_model(checkpoint_path, device):
    checkpoint = torch.load(checkpoint_path, map_location=device,
                            weights_only=False)
    model = ProteinMPNN(
        num_letters=21, node_features=128, edge_features=128,
        hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
        augment_eps=checkpoint.get("noise_level", 0.2),
        k_neighbors=checkpoint.get("num_edges", 48),
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    model.eval()
    return model


def parse_mpnn_fasta(fasta_path):
    entries = []
    header, seq_lines = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def load_fixed_positions(jsonl_path):
    """Load fixed positions JSONL.
    Returns dict: {chain_id: list_of_fixed_resnums}
    """
    with open(jsonl_path) as f:
        data = json.loads(f.readline().strip())
    return list(data.values())[0]


def score_sequences_for_design(pdb_path, fasta_path, fixed_dict,
                               model, device):
    """Score all designed sequences for one design PDB.

    Args:
        pdb_path: Path to design PDB (two-chain: A=receptor, B=effector)
        fasta_path: MPNN output FASTA with designed sequences
        fixed_dict: {chain_id: [list of fixed resnums]} from JSONL
        model: Loaded ProteinMPNN model
        device: torch device

    Returns:
        list of (seq_idx, global_score, design_region_score, n_free)
    """
    # Parse PDB — returns list of dicts with seq_chain_X, coords_chain_X
    pdb_dict_list = parse_PDB(pdb_path)
    entry = pdb_dict_list[0]
    pdb_name = entry["name"]

    # Identify chains from the parsed dict
    chain_ids = []
    for key in entry:
        if key.startswith("seq_chain_"):
            chain_ids.append(key.replace("seq_chain_", ""))
    chain_ids.sort()

    # chain_dict: {name: (masked_chains, visible_chains)}
    # All chains are "masked" (passed through decoder), fixed positions
    # are handled separately via fixed_position_dict
    chain_dict = {pdb_name: [chain_ids, []]}

    # fixed_position_dict: {name: {chain: [resnums]}}
    fixed_position_dict = {pdb_name: fixed_dict}

    # Featurize
    result = tied_featurize(
        pdb_dict_list, device, chain_dict,
        fixed_position_dict=fixed_position_dict,
    )
    # tied_featurize returns 20 values. The key ones (verified empirically):
    #   [0]  X              (1, L, 4, 3)  coordinates
    #   [1]  S              (1, L)         native sequence indices
    #   [2]  mask           (1, L)         valid position mask
    #   [4]  chain_M        (1, L)         chain-level mask (1=designable chain)
    #   [5]  chain_encoding (1, L)         chain ID encoding
    #   [10] chain_M_pos    (1, L)         position-level mask (0=fixed, 1=free)
    #   [12] residue_idx    (1, L)         residue index encoding
    X             = result[0]
    S_native      = result[1]
    mask          = result[2]
    chain_M       = result[4]
    chain_encoding = result[5]
    chain_M_pos   = result[10]
    residue_idx   = result[12]

    seq_len = S_native.shape[1]

    # chain_M_pos: 0.0 = fixed, 1.0 = free (designable)
    free_mask = chain_M_pos.squeeze(0).cpu().numpy().astype(bool)
    mask_np = mask.squeeze(0).cpu().numpy().astype(bool)

    # Parse designed sequences
    entries = parse_mpnn_fasta(fasta_path)
    designed_entries = entries[1:]  # skip native (index 0)

    results = []

    for seq_idx, (header, full_seq) in enumerate(designed_entries):
        # Convert sequence to indices — handle chain separators
        clean_seq = full_seq.replace("/", "").replace(":", "")
        S = torch.tensor(
            [[AA_TO_IDX.get(aa, 20) for aa in clean_seq]],
            dtype=torch.long, device=device,
        )

        # Pad/truncate to match structure
        if S.shape[1] < seq_len:
            pad = torch.full((1, seq_len - S.shape[1]), 20,
                             dtype=torch.long, device=device)
            S = torch.cat([S, pad], dim=1)
        elif S.shape[1] > seq_len:
            S = S[:, :seq_len]

        # Forward pass — model.forward() requires a randn tensor for
        # decoder noise (set to zeros for deterministic scoring)
        randn = torch.zeros(X.shape[0], X.shape[1], device=device)

        with torch.no_grad():
            log_probs = model(X, S, mask, chain_M, residue_idx,
                              chain_encoding, randn)
            # log_probs: (1, L, 21)

        lp = log_probs.squeeze(0).cpu().numpy()  # (L, 21)
        s = S.squeeze(0).cpu().numpy()            # (L,)

        # Per-position score: negative log-prob of the designed AA
        per_pos = np.zeros(seq_len)
        for i in range(seq_len):
            if s[i] < 21 and mask_np[i]:
                per_pos[i] = -lp[i, s[i]]

        # Global score: mean over all valid positions
        valid = mask_np[:seq_len]
        global_score = float(per_pos[valid].mean()) if valid.any() else 0.0

        # Design region score: mean over free positions only
        fm = free_mask[:seq_len] & valid
        n_free = int(fm.sum())
        design_score = float(per_pos[fm].mean()) if n_free > 0 else None

        results.append((seq_idx, global_score, design_score, n_free))

    return results


def main():
    args = parse_args()

    if not HAS_MPNN:
        print("ERROR: ProteinMPNN not available", file=sys.stderr)
        sys.exit(1)

    checkpoint_path = args.checkpoint or find_checkpoint()
    if not checkpoint_path:
        print("ERROR: No checkpoint found", file=sys.stderr)
        sys.exit(1)
    print(f"Checkpoint: {checkpoint_path}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    model = load_model(checkpoint_path, device)
    print("Model loaded")

    # Read metadata
    rows = []
    with open(args.metadata_csv) as f:
        for row in csv.DictReader(f):
            rows.append(row)

    if not rows:
        with open(args.output_csv, "w") as f:
            f.write("design,seq,mpnn_score,design_region_score\n")
        return

    for row in rows:
        row["design_region_score"] = ""

    designs_seen = set()
    n_scored = 0

    for row in rows:
        design_idx = row["design"]
        if design_idx in designs_seen:
            continue
        designs_seen.add(design_idx)
        design_name = f"design_{design_idx}"

        pdb_path = os.path.join(args.pdb_dir, f"{design_name}.pdb")
        jsonl_path = os.path.join(args.fixed_positions_dir,
                                  f"fixed_positions_{design_name}.jsonl")
        fasta_path = os.path.join(args.mpnn_fasta_dir,
                                  f"{design_name}_mpnn", "seqs",
                                  f"{design_name}.fa")

        missing = []
        if not os.path.exists(pdb_path): missing.append("PDB")
        if not os.path.exists(jsonl_path): missing.append("JSONL")
        if not os.path.exists(fasta_path): missing.append("FASTA")
        if missing:
            print(f"  {design_name}: missing {', '.join(missing)}")
            continue

        fixed_dict = load_fixed_positions(jsonl_path)

        n_free_total = 0
        for chain, resnums in fixed_dict.items():
            n_free_total += 0  # we'll get this from the model

        print(f"  {design_name}: scoring...")

        try:
            results = score_sequences_for_design(
                pdb_path, fasta_path, fixed_dict, model, device)

            for seq_idx, global_s, region_s, n_free in results:
                for r in rows:
                    if r["design"] == design_idx and r["seq"] == str(seq_idx):
                        if region_s is not None:
                            r["design_region_score"] = f"{region_s:.4f}"
                            n_scored += 1
                        region_str = f"{region_s:.4f}" if region_s is not None else "N/A"
                        print(f"    seq_{seq_idx}: global={global_s:.4f}, "
                              f"region={region_str} ({n_free} free)")
                        break

        except Exception as e:
            print(f"  {design_name}: ERROR: {e}")
            traceback.print_exc()

    # Write output
    out_fieldnames = list(rows[0].keys())
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nDesign region scores: {n_scored}/{len(rows)} sequences scored")


if __name__ == "__main__":
    main()