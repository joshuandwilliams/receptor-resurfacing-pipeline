"""
haddock_run.py
--------------
HaddockRun — one HADDOCK3 docking call's complete state.

Per Session 7 grill-me (notes/design_audit.md Q153, Q151).  Tier-N type
analogous to NegativeSteeringRun: a typed wrapper over the on-disk state
of one HADDOCK run, with a from_workdir factory and the auto-pick logic.

Composition:
- Receptor + effector input PDBs (always monomers per Q136).
- Restraints (flat fields per Q151): contact_pairs, receptor_active_residues,
  effector_active_residues, pair_distance.
- Tuple[HaddockCluster, ...] for the qualifying clusters.
- Optional chosen_cluster_id (None ⇒ auto-pick by HaddockCluster.auto_pick_key).

The "design region at HADDOCK time" for clash bookkeeping comes from
``receptor_active_residues ∪ {p.rec_resnum for p in contact_pairs}``.
Per A139 the contig string is intentionally NOT a HaddockRun field —
clash bookkeeping is fully expressible from the restraint params alone.

Restraint parsing
=================
The string formats are:
- Contact pair: ``"A25-C42 A13-C94"`` — chain-prefixed receptor resnum, dash,
  chain-prefixed effector resnum; space-separated pairs.  Case-insensitive
  on chain letters.
- Active residues: ``"25,35,40-44"`` — comma-separated residue numbers
  and ``N-M`` ranges; no chain prefix (chain is given by params.receptor_chain
  or params.effector_chain).
- Pair distance: ``"2,2,4"`` — three positive floats (target, lo_dev, hi_dev).

Per A148 the validator owns input-format checking via regex; the parsers
below are defensive but assume the input has passed validation.
"""

from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Set, Tuple

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from haddock_cluster import HaddockCluster  # noqa: E402


# ── Restraint-string parsers (importable; used by haddock3_prepare too) ──


_PAIR_RE = re.compile(r"^([A-Za-z])(\d+)-([A-Za-z])(\d+)$")


def parse_contact_pairs(s: str) -> Tuple[Tuple[str, int, str, int], ...]:
    """Parse ``"A25-C42 A13-C94"`` into a tuple of
    ``(rec_chain, rec_resnum, eff_chain, eff_resnum)`` tuples.

    Empty / whitespace-only input returns an empty tuple.  Chain letters
    are uppercased.
    """
    if not s or not s.strip():
        return ()
    pairs: List[Tuple[str, int, str, int]] = []
    for token in s.split():
        m = _PAIR_RE.match(token)
        if not m:
            raise ValueError(
                f"parse_contact_pairs: bad pair {token!r}; "
                f"expected format like 'A25-C42'"
            )
        rec_chain, rec_resnum, eff_chain, eff_resnum = m.groups()
        pairs.append(
            (rec_chain.upper(), int(rec_resnum),
             eff_chain.upper(), int(eff_resnum))
        )
    return tuple(pairs)


def parse_active_residues(s: str) -> Tuple[int, ...]:
    """Parse ``"25,35,40-44"`` into a sorted tuple of unique residue
    numbers.  Empty input returns an empty tuple.
    """
    if not s or not s.strip():
        return ()
    out: Set[int] = set()
    for token in s.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            lo, hi = token.split("-", 1)
            lo_i, hi_i = int(lo), int(hi)
            if lo_i > hi_i:
                raise ValueError(
                    f"parse_active_residues: bad range {token!r} (lo > hi)"
                )
            out.update(range(lo_i, hi_i + 1))
        else:
            out.add(int(token))
    return tuple(sorted(out))


def parse_pair_distance(s: str) -> Tuple[float, float, float]:
    """Parse ``"2,2,4"`` into ``(target, lo_dev, hi_dev)``.  All values
    must be non-negative.
    """
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 3:
        raise ValueError(
            f"parse_pair_distance: expected 3 comma-separated floats, "
            f"got {s!r}"
        )
    try:
        target, lo_dev, hi_dev = (float(p) for p in parts)
    except ValueError as e:
        raise ValueError(f"parse_pair_distance: {e}") from None
    if min(target, lo_dev, hi_dev) < 0.0:
        raise ValueError(
            f"parse_pair_distance: all values must be >= 0, got {s!r}"
        )
    return (target, lo_dev, hi_dev)


# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class HaddockRun:
    """One HADDOCK3 docking run wrapped as a typed read-only view.

    Construct directly with parsed restraint fields + a list of
    HaddockCluster instances; or use ``from_workdir`` to hydrate from a
    HADDOCK3_DOCK output directory (haddock_report.json +
    cluster_metrics.json + best_cluster*.pdb files).
    """
    receptor_pdb: Path
    effector_pdb: Path
    # ── Restraints (flat) ───────────────────────────────────────────
    contact_pairs: Tuple[Tuple[str, int, str, int], ...] = ()
    receptor_active_residues: Tuple[int, ...] = ()
    effector_active_residues: Tuple[int, ...] = ()
    pair_distance: Tuple[float, float, float] = (2.0, 2.0, 4.0)
    # Contig-derived design region (post-commit-3 amendment): the
    # canonical source of "which receptor residues are in the design
    # region" for clash bookkeeping.  Empty tuple means "no contig was
    # supplied to haddock3_prepare.py" and design_region falls back to
    # pair receptor halves + receptor_active_residues.
    contig_design_region: Tuple[int, ...] = ()
    # ── Outputs ─────────────────────────────────────────────────────
    qualifying_clusters: Tuple[HaddockCluster, ...] = field(default_factory=tuple)
    chosen_cluster_id: Optional[int] = None

    def __post_init__(self) -> None:
        # Coerce to tuples in case the caller passed lists.
        if not isinstance(self.qualifying_clusters, tuple):
            object.__setattr__(
                self, "qualifying_clusters",
                tuple(self.qualifying_clusters),
            )
        if not isinstance(self.contact_pairs, tuple):
            object.__setattr__(self, "contact_pairs", tuple(self.contact_pairs))
        if not isinstance(self.receptor_active_residues, tuple):
            object.__setattr__(
                self, "receptor_active_residues",
                tuple(self.receptor_active_residues),
            )
        if not isinstance(self.effector_active_residues, tuple):
            object.__setattr__(
                self, "effector_active_residues",
                tuple(self.effector_active_residues),
            )
        if not isinstance(self.contig_design_region, tuple):
            object.__setattr__(
                self, "contig_design_region",
                tuple(self.contig_design_region),
            )

        for c in self.qualifying_clusters:
            if not isinstance(c, HaddockCluster):
                raise TypeError(
                    "qualifying_clusters must contain HaddockCluster instances"
                )
        if len(self.pair_distance) != 3:
            raise ValueError(
                f"pair_distance must be a 3-tuple, got {self.pair_distance!r}"
            )
        if self.chosen_cluster_id is not None:
            ids = {c.cluster_id for c in self.qualifying_clusters}
            if self.chosen_cluster_id not in ids:
                raise ValueError(
                    f"chosen_cluster_id {self.chosen_cluster_id} not in "
                    f"qualifying_clusters {sorted(ids)}"
                )

    # ── Derived sets ────────────────────────────────────────────────

    @property
    def design_region(self) -> Set[int]:
        """Receptor residues that count as "in the design region" for
        clash bookkeeping.

        Preferred source (post-commit-3 amendment): contig_design_region,
        derived by haddock3_prepare.py from the contig string + receptor
        PDB.  Falls back to ``receptor halves of contact_pairs ∪
        receptor_active_residues`` when no contig was passed in.
        """
        if self.contig_design_region:
            return set(self.contig_design_region)
        rec_from_pairs = {p[1] for p in self.contact_pairs}
        return rec_from_pairs | set(self.receptor_active_residues)

    @property
    def has_restraints(self) -> bool:
        return bool(self.contact_pairs) or bool(self.receptor_active_residues)

    # ── Cluster selection ───────────────────────────────────────────

    @property
    def selected(self) -> HaddockCluster:
        """The HaddockCluster the pipeline picks downstream.

        Per Session 7 A134/A154:
        - If ``chosen_cluster_id`` is set, return that cluster.
        - Otherwise auto-pick by lexicographic ``(-pair_contact_fraction,
          -bsa)`` — highest pair satisfaction first, BSA as tie-breaker.

        Raises if there are no qualifying clusters.
        """
        if not self.qualifying_clusters:
            raise ValueError("HaddockRun.selected: no qualifying clusters")
        if self.chosen_cluster_id is not None:
            for c in self.qualifying_clusters:
                if c.cluster_id == self.chosen_cluster_id:
                    return c
            raise AssertionError("post_init should have caught this")
        return min(self.qualifying_clusters, key=lambda c: c.auto_pick_key)

    def to_rfdiffusion_input(self) -> Path:
        """The PDB file path that goes to BUILD_CONTIGS → RFDiffusion."""
        return self.selected.best_model_pdb

    # ── Factory ─────────────────────────────────────────────────────

    @classmethod
    def from_workdir(
        cls,
        workdir: Path,
        receptor_pdb: Path,
        effector_pdb: Path,
        contact_pairs: str = "",
        receptor_active_residues: str = "",
        effector_active_residues: str = "",
        pair_distance: str = "2,2,4",
        chosen_cluster_id: Optional[int] = None,
    ) -> "HaddockRun":
        """Hydrate from a HADDOCK3_DOCK output directory.

        Expects:
        - ``<workdir>/haddock_report.json`` — written by
          collect_haddock3_dock.py; carries the cluster_models map.
        - ``<workdir>/cluster_metrics.json`` — written by
          haddock_cluster_metrics.py; carries the per-cluster metric values
          keyed by cluster_id.
        - ``<workdir>/best_cluster*.pdb`` — per-cluster best model PDBs.

        Restraint params come in as strings (the same format as the YAML
        params) so the factory is the single conversion point.

        Returns a HaddockRun with the parsed restraints and the
        constructed HaddockCluster tuple.  Bridge distances are NOT
        filled here — they're a post-construction step (Q152(i)).
        """
        workdir = Path(workdir)
        report_path = workdir / "haddock_report.json"
        metrics_path = workdir / "cluster_metrics.json"
        sc_path = workdir / "cluster_sc.json"
        restraints_path = workdir / "restraints_summary.json"
        if not report_path.is_file():
            raise FileNotFoundError(
                f"HaddockRun.from_workdir: {report_path} not found"
            )
        if not metrics_path.is_file():
            raise FileNotFoundError(
                f"HaddockRun.from_workdir: {metrics_path} not found"
            )

        with open(report_path) as f:
            report = json.load(f)
        with open(metrics_path) as f:
            metrics = json.load(f)
        # cluster_sc.json (sibling output from HADDOCK_CLUSTER_SC in
        # rosetta_container) is optional — if missing, sc stays None on
        # every HaddockCluster.
        sc_by_cluster: dict = {}
        if sc_path.is_file():
            with open(sc_path) as f:
                sc_by_cluster = json.load(f)
        # restraints_summary.json may be absent for legacy / synthetic fixtures.
        contig_design_region: Tuple[int, ...] = ()
        if restraints_path.is_file():
            with open(restraints_path) as f:
                restraints = json.load(f)
            contig_design_region = tuple(restraints.get("contig_design_region", ()))

        cluster_models = report.get("cluster_models", {})
        clusters: List[HaddockCluster] = []
        for cid_str, info in cluster_models.items():
            cid = int(cid_str)
            m = metrics.get(str(cid)) or metrics.get(cid)
            if m is None:
                # No metrics for this cluster — skip.  collect_haddock3_dock
                # should ensure metrics are present for every reported
                # cluster, so this is a defensive guard.
                continue
            best_pdb = workdir / info["filename"]
            # Sc may come from either cluster_metrics.json (legacy single-
            # process flow) or cluster_sc.json (post-commit-4 split).
            sc_value = m.get("sc")
            if sc_value is None and str(cid) in sc_by_cluster:
                sc_value = sc_by_cluster[str(cid)].get("sc")
            cluster = HaddockCluster(
                cluster_id=cid,
                size=int(info["size"]),
                mean_haddock_score=float(info["mean_score"]),
                best_model_pdb=best_pdb,
                bsa=float(m["bsa"]),
                com_distance=float(m["com_distance"]),
                air_satisfaction_count=int(m["air_satisfaction_count"]),
                air_total_count=int(m["air_total_count"]),
                pair_contact_fraction=float(m["pair_contact_fraction"]),
                pair_total_count=int(m["pair_total_count"]),
                clashes_in_design_region=int(m["clashes_in_design_region"]),
                clashes_outside_design_region=int(m["clashes_outside_design_region"]),
                sc=None if sc_value is None else float(sc_value),
            )
            clusters.append(cluster)

        return cls(
            receptor_pdb=Path(receptor_pdb),
            effector_pdb=Path(effector_pdb),
            contact_pairs=parse_contact_pairs(contact_pairs),
            receptor_active_residues=parse_active_residues(receptor_active_residues),
            effector_active_residues=parse_active_residues(effector_active_residues),
            pair_distance=parse_pair_distance(pair_distance),
            contig_design_region=contig_design_region,
            qualifying_clusters=tuple(clusters),
            chosen_cluster_id=chosen_cluster_id,
        )
