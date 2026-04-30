#!/usr/bin/env python3
"""
sequence_registry.py
--------------------
Lightweight registry for tracking which receptor sequences have been
predicted by Boltz, how many seeds each has, and where the results live.

Used by both the steered-prediction and reversion passes to avoid
duplicate Boltz predictions.  A sequence that was predicted as a cold
start doesn't need to be predicted again during reversion, and 40
designs that revert to the same wild-type sequence only need one set
of predictions.

The registry is a JSON file at the experiment root keyed by a hash of
the full receptor sequence.  Each entry stores:
    - sequence: the full receptor sequence (for collision checking)
    - prediction_dirs: list of paths where Boltz outputs live
    - seeds_used: list of Boltz seeds that have been run
    - source: how this sequence entered the registry (e.g. "cold_start",
              "steered", "reversion")
    - labels: list of design labels that share this sequence

Thread safety: the registry is read at plan time (single process) and
updated at plan time (single process), so no locking is needed.  The
predict-one array tasks only READ plan.json, they never touch the
registry.

Pure stdlib — no numpy, no Bio.  Safe to import from host-side code.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional


def _seq_hash(sequence: str) -> str:
    """Stable hash of a receptor sequence for registry keying."""
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()[:16]


def load_registry(experiment_root: Path) -> Dict:
    """Load the sequence registry from disk, or return an empty one."""
    path = experiment_root / "sequence_registry.json"
    if path.exists():
        try:
            return json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def save_registry(experiment_root: Path, registry: Dict) -> Path:
    """Write the registry to disk.  Returns the path written."""
    path = experiment_root / "sequence_registry.json"
    path.write_text(json.dumps(registry, indent=2))
    return path


def lookup(registry: Dict, sequence: str) -> Optional[Dict]:
    """Look up a sequence in the registry.
    
    Returns the entry dict if found, None otherwise.
    Checks for hash collisions by comparing the full sequence.
    """
    h = _seq_hash(sequence)
    entry = registry.get(h)
    if entry is not None and entry.get("sequence") == sequence:
        return entry
    return None


def register(
    registry: Dict,
    sequence: str,
    prediction_dirs: List[str],
    seeds_used: List[int],
    source: str,
    labels: List[str],
) -> Dict:
    """Register a sequence as predicted.
    
    If the sequence is already registered, appends the new prediction
    dirs, seeds, and labels to the existing entry.  Returns the
    (possibly updated) entry.
    """
    h = _seq_hash(sequence)
    existing = registry.get(h)
    
    if existing is not None and existing.get("sequence") == sequence:
        # Merge into existing entry
        existing["prediction_dirs"].extend(prediction_dirs)
        existing["seeds_used"].extend(seeds_used)
        existing["labels"].extend(labels)
        # Deduplicate
        existing["prediction_dirs"] = list(dict.fromkeys(
            existing["prediction_dirs"]))
        existing["seeds_used"] = sorted(set(existing["seeds_used"]))
        existing["labels"] = list(dict.fromkeys(existing["labels"]))
        existing["sources"] = list(set(
            existing.get("sources", [existing.get("source", "")]) + [source]
        ))
        return existing
    
    # New entry
    entry = {
        "sequence": sequence,
        "hash": h,
        "prediction_dirs": list(prediction_dirs),
        "seeds_used": sorted(seeds_used),
        "source": source,
        "sources": [source],
        "labels": list(labels),
    }
    registry[h] = entry
    return entry


def count_existing_seeds(registry: Dict, sequence: str) -> int:
    """How many seeds have already been predicted for this sequence?"""
    entry = lookup(registry, sequence)
    if entry is None:
        return 0
    return len(entry.get("seeds_used", []))


def get_existing_dirs(registry: Dict, sequence: str) -> List[str]:
    """Return the prediction directories for an already-predicted sequence."""
    entry = lookup(registry, sequence)
    if entry is None:
        return []
    return list(entry.get("prediction_dirs", []))


def deduplicate_sequences(
    sequences: List[str],
    labels: List[str],
    registry: Dict,
    num_seeds: int,
) -> Dict[str, Dict]:
    """Given a list of sequences (with corresponding labels), determine
    which need new predictions and how many seeds each needs.
    
    Returns a dict keyed by sequence, with values:
        {
            "labels": [labels sharing this sequence],
            "existing_seeds": int,
            "new_seeds_needed": int,
            "total_seeds_target": num_seeds,
            "already_predicted": bool,
        }
    
    Sequences that already have >= num_seeds in the registry need
    zero new predictions.  Those with fewer need (num_seeds - existing)
    new seeds.
    """
    # Group labels by sequence
    by_seq: Dict[str, List[str]] = {}
    for seq, label in zip(sequences, labels):
        by_seq.setdefault(seq, []).append(label)
    
    result = {}
    for seq, seq_labels in by_seq.items():
        existing = count_existing_seeds(registry, seq)
        needed = max(0, num_seeds - existing)
        result[seq] = {
            "labels": seq_labels,
            "existing_seeds": existing,
            "new_seeds_needed": needed,
            "total_seeds_target": num_seeds,
            "already_predicted": existing >= num_seeds,
        }
    return result
