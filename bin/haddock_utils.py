#!/usr/bin/env python3
"""
haddock_utils.py
----------------
Shared utilities for HADDOCK3 post-processing scripts.

Provides TSV parsing, column lookup, PDB I/O, model-stem matching,
and NumPy-accelerated interface contact detection used by both
collect_haddock3_dock.py and haddock3_plots.py.
"""

import glob
import gzip
import os
import shutil

import numpy as np


# ── TSV Parsing ──────────────────────────────────────────────────────────────

def parse_capri_tsv(capri_file):
    """Parse capri_ss.tsv into (rows, header). Each row is a dict keyed by column name."""
    lines = open(capri_file).readlines()
    if len(lines) < 2:
        return [], []
    header = lines[0].strip().split("\t")
    rows = []
    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) >= len(header):
            rows.append(dict(zip(header, parts)))
    return rows, header


def parse_clustfcc_tsv(clust_file, min_cluster_size=4):
    """
    Parse clustfcc.tsv into cluster membership data.

    The default min_cluster_size matches HADDOCK3's `min_population` default
    and the haddock.nf clustfcc block.  Callers should pass an explicit value
    if `min_population` in the HADDOCK3 config has been overridden.

    Returns:
        clusters          : {cid: [(model_name, score), ...]}
        n_clusters        : number of named clusters (excluding '-')
        best_qual_cluster : cid with lowest mean score among clusters
                            with >= min_cluster_size members, or None
    """
    lines = open(clust_file).readlines()
    if len(lines) < 2:
        return {}, 0, None
    header = lines[0].strip().split("\t")
    clusters = {}

    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) < len(header):
            continue
        row = dict(zip(header, parts))

        cid = _first_col(row, ("cluster_id", "cluster-id", "clust_id", "cluster"))
        if cid is None:
            continue

        model = _first_col(row, ("model_name", "model", "structure", "pdb"))
        if model is None:
            model = f"model_{sum(len(v) for v in clusters.values())}"

        score = _first_float(row, ("score", "haddock-score", "total"))
        clusters.setdefault(cid, []).append((model, score))

    named = {c: m for c, m in clusters.items() if c != "-"}
    n_clusters = len(named)

    # Best qualifying cluster: lowest mean score among large-enough clusters
    qualifying = {c: m for c, m in named.items() if len(m) >= min_cluster_size}
    best_qual_cluster = None
    if qualifying:
        best_qual_cluster = min(qualifying, key=lambda c: cluster_mean_score(qualifying[c]))

    return clusters, n_clusters, best_qual_cluster


# ── Column Helpers ───────────────────────────────────────────────────────────

def _first_col(row, candidates):
    """Return the value of the first matching column name, or None."""
    for col in candidates:
        if col in row:
            return row[col]
    return None


def _first_float(row, candidates):
    """Return the float value of the first matching column, or None."""
    for col in candidates:
        if col in row:
            try:
                return float(row[col])
            except (ValueError, TypeError):
                pass
    return None


def get_score(row):
    """Extract the HADDOCK score from a capri/cluster row dict."""
    return _first_float(row, ("score", "haddock-score", "total"))


def get_model_name(row):
    """Extract the model filename from a capri row dict."""
    name = _first_col(row, ("model", "model_name", "structure", "pdb"))
    return name if name is not None else next(iter(row.values()), None)


def get_numeric_col(data, name_options):
    """Extract a numeric column from parsed TSV data. Returns list of floats or None."""
    for name in name_options:
        if data and name in data[0]:
            try:
                return [float(row[name]) for row in data]
            except (ValueError, TypeError):
                continue
    return None


def cluster_mean_score(members):
    """Mean HADDOCK score for a list of (model, score) tuples, ignoring None."""
    scored = [s for (_, s) in members if s is not None]
    return sum(scored) / len(scored) if scored else float("inf")


# ── Model Stem Matching ─────────────────────────────────────────────────────

def model_stem(model_path):
    """
    Bare filename stem for matching against PDB files on disk.
    Strips directory components and .pdb/.pdb.gz extensions.
    e.g. '../emref/emref_42.pdb' → 'emref_42'
    """
    base = os.path.basename(model_path)
    for ext in (".pdb.gz", ".pdb"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


# ── PDB I/O ──────────────────────────────────────────────────────────────────

def iter_pdb_lines(pdb_path):
    """Yield lines from a plain or gzip-compressed PDB file."""
    if pdb_path.endswith(".gz"):
        with gzip.open(pdb_path, "rt", errors="replace") as fh:
            yield from fh
    else:
        with open(pdb_path, errors="replace") as fh:
            yield from fh


def extract_heavy_atoms(pdb_path, chain_ids):
    """
    Extract heavy-atom coordinates for one or more chains in a single pass.

    Parameters
    ----------
    pdb_path  : path to PDB or .pdb.gz
    chain_ids : str or iterable of chain ID characters

    Returns
    -------
    dict : {chain_id: {resnum: [(x, y, z), ...]}}
    """
    if isinstance(chain_ids, str):
        chain_ids = (chain_ids,)
    chain_set = set(chain_ids)
    atoms = {c: {} for c in chain_set}

    for line in iter_pdb_lines(pdb_path):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        if chain not in chain_set:
            continue
        atom_name = line[12:16].strip()
        if atom_name.startswith("H") or (len(atom_name) > 1 and atom_name[1] == "H"):
            continue
        try:
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atoms[chain].setdefault(resnum, []).append((x, y, z))
    return atoms


def copy_pdb(src, dst):
    """Copy a PDB or .pdb.gz to dst, decompressing .gz if needed."""
    if src.endswith(".gz"):
        with gzip.open(src, "rb") as gz:
            with open(dst, "wb") as f:
                f.write(gz.read())
    else:
        shutil.copy(src, dst)


# ── Contact Detection (NumPy-accelerated) ────────────────────────────────────

def contacted_residues(query_atoms, target_atoms, cutoff):
    """
    Return the set of query residue numbers that have at least one heavy atom
    within *cutoff* Å of any target heavy atom.

    Uses NumPy broadcasting for vectorised distance computation.

    Parameters
    ----------
    query_atoms  : {resnum: [(x,y,z), ...]} — residues to test
    target_atoms : {resnum: [(x,y,z), ...]} — partner chain atoms
    cutoff       : distance threshold in Angstroms

    Returns
    -------
    set of query residue numbers at the interface
    """
    if not query_atoms or not target_atoms:
        return set()

    target_coords = np.array([c for coords in target_atoms.values() for c in coords])
    cutoff_sq = cutoff ** 2
    contacted = set()

    for resnum, coords in query_atoms.items():
        arr = np.array(coords)  # (M, 3)
        diff = arr[:, np.newaxis, :] - target_coords[np.newaxis, :, :]  # (M, K, 3)
        if (diff ** 2).sum(axis=2).min() < cutoff_sq:
            contacted.add(resnum)

    return contacted


def find_interface_residues(rec_atoms, eff_atoms, cutoff):
    """
    Return (interface_eff_residues, interface_rec_residues) as sorted lists.
    Uses NumPy-accelerated distance calculation.
    """
    rec_contacts = contacted_residues(rec_atoms, eff_atoms, cutoff)
    eff_contacts = contacted_residues(eff_atoms, rec_atoms, cutoff)
    return sorted(eff_contacts), sorted(rec_contacts)


# ── PDB Index Builder ────────────────────────────────────────────────────────

def collect_pdb_index(complex_dir=None, run_dir=None):
    """
    Build a {stem: filepath} index for docked complex PDB files.

    Search order: run_dir (seletopclusts → emref → full), then complex_dir.
    """
    index = {}

    def _index_dir(d):
        if not d or not os.path.isdir(d):
            return
        for pattern in ("**/*.pdb.gz", "**/*.pdb"):
            for fpath in glob.glob(os.path.join(d, pattern), recursive=True):
                stem = model_stem(fpath)
                if stem not in index or fpath.endswith(".pdb"):
                    index[stem] = fpath

    if run_dir and os.path.isdir(run_dir):
        for step_glob in ("*seletopclusts*", "*emref*"):
            for step_dir in sorted(glob.glob(os.path.join(run_dir, step_glob))):
                _index_dir(step_dir)
        if index:
            return index
        _index_dir(run_dir)
        if index:
            return index

    _index_dir(complex_dir)
    return index