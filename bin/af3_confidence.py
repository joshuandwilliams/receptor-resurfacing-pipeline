"""
af3_confidence.py
-----------------
AF3ConfidenceAggregate — the aggregated structural + confidence summary
from ONE AF3-no-MSA call.

Boltz gives one prediction per PDB.  AF3 instead emits one per
seed × diffusion sample, so this type aggregates across that set.  It
reports best and mean ra_eff against ground truth, best and mean iptm,
and the count of predictions whose interface placement passes the
threshold.

`bin/parse_af3_output.py` is the CLI over this type.  It imports
AF3ConfidenceAggregate and adds the one column this class does not
produce, `af3_nomsa_n_correct_interface`, which needs a thresholds
object the class takes only as a method argument.

Scope limitation: this class exposes aggregate numbers, not the
structures behind them.  `ProteinStructurePrediction` exists in
`bin/protein_structure_prediction.py`, but no `best_prediction()` or
`all_predictions()` accessor was ever added here.  Use the `cif_paths`
property to reach the underlying files.

Imports numpy and gemmi, and no other module in bin/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ═════════════════════════════════════════════════════════════════════════
# Private helpers
# ═════════════════════════════════════════════════════════════════════════


def _chain_ca_coords(structure, chain_id: str) -> List[np.ndarray]:
    """Return Cα (x, y, z) coords for the given chain, in the order they
    appear in the structure.

    Callers pair these by position, not by residue number, because the
    sources number residues differently.  AF3 emits per-chain 1-based
    numbering.  RFDiffusion, Boltz and the cleaned ground-truth PDBs
    number continuously across chains.
    """
    coords: List[np.ndarray] = []
    for model in structure:
        for chain in model:
            if chain.name != chain_id:
                continue
            for res in chain:
                for atom in res:
                    if atom.name == "CA":
                        coords.append(np.array(
                            [atom.pos.x, atom.pos.y, atom.pos.z]
                        ))
                        break
    return coords


def _receptor_aligned_effector_rmsd(
    pred_cif_path: Path,
    ref_pdb_path: Path,
    receptor_chain: str,
    effector_chain: str,
) -> Tuple[Optional[float], Optional[str]]:
    """Return receptor-aligned effector Cα RMSD, or (None, error_tag)."""
    import gemmi
    try:
        pred = gemmi.read_structure(str(pred_cif_path))
        ref = gemmi.read_structure(str(ref_pdb_path))
    except Exception as e:  # noqa: BLE001
        return None, f"structure_load_failed:{type(e).__name__}"

    try:
        pred_rec = _chain_ca_coords(pred, receptor_chain)
        pred_eff = _chain_ca_coords(pred, effector_chain)
        ref_rec = _chain_ca_coords(ref, receptor_chain)
        ref_eff = _chain_ca_coords(ref, effector_chain)
    except Exception as e:  # noqa: BLE001
        return None, f"chain_parse_failed:{type(e).__name__}"

    rec_n = min(len(pred_rec), len(ref_rec))
    eff_n = min(len(pred_eff), len(ref_eff))
    if rec_n < 3:
        return None, (f"too_few_receptor_residues:{rec_n} "
                      f"(pred={len(pred_rec)}, ref={len(ref_rec)})")
    if len(pred_rec) != len(ref_rec):
        print(f"[{pred_cif_path.parent.name}] receptor length mismatch: "
              f"pred={len(pred_rec)} ref={len(ref_rec)}; "
              f"pairing on first {rec_n}", file=sys.stderr)

    P = np.array(pred_rec[:rec_n])
    Q = np.array(ref_rec[:rec_n])
    P_c, P_m = P - P.mean(axis=0), P.mean(axis=0)
    Q_c, Q_m = Q - Q.mean(axis=0), Q.mean(axis=0)
    H = P_c.T @ Q_c
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T

    if eff_n < 3:
        return None, (f"too_few_effector_residues:{eff_n} "
                      f"(pred={len(pred_eff)}, ref={len(ref_eff)})")
    if len(pred_eff) != len(ref_eff):
        print(f"[{pred_cif_path.parent.name}] effector length mismatch: "
              f"pred={len(pred_eff)} ref={len(ref_eff)}; "
              f"pairing on first {eff_n}", file=sys.stderr)

    P_eff = np.array(pred_eff[:eff_n])
    Q_eff = np.array(ref_eff[:eff_n])
    P_eff_aligned = (P_eff - P_m) @ R.T + Q_m
    diff = P_eff_aligned - Q_eff
    return float(np.sqrt((diff * diff).sum(axis=1).mean())), None


def _read_confidence(cif_path: Path) -> Tuple[Optional[float], Optional[str]]:
    """Read iPTM from the AF3 summary_confidences.json sibling of `cif_path`."""
    for candidate in [
        cif_path.parent / "summary_confidences.json",
        cif_path.parent
            / f"{cif_path.stem.replace('_model', '')}_summary_confidences.json",
    ]:
        if candidate.is_file():
            try:
                with open(candidate) as f:
                    data = json.load(f)
                iptm = data.get("iptm")
                if iptm is not None:
                    return float(iptm), None
            except (json.JSONDecodeError, OSError) as e:
                return None, f"confidence_load_failed:{type(e).__name__}"
    return None, "confidence_json_not_found"


# ═════════════════════════════════════════════════════════════════════════
# AF3ConfidenceAggregate
# ═════════════════════════════════════════════════════════════════════════


class AF3ConfidenceAggregate:
    """The aggregated AF3-no-MSA summary for one survivor.

    Walks the AF3 output tree and computes ra_eff against the supplied
    ground-truth PDB for every (seed × diffusion-sample) prediction,
    then aggregates to best, mean and count-passing.

    The walk happens on first access to any aggregate property, and
    the result is cached for later accesses.

    Construction:
        AF3ConfidenceAggregate.from_output_dir(
            af3_output_dir,
            ground_truth,
            receptor_chain="A",
            effector_chain="B",
        )

    `__init__` takes the same arguments.  Tests populate the private
    caches directly to supply parsed values.
    """

    def __init__(
        self,
        af3_output_dir: Path,
        ground_truth: Path,
        receptor_chain: str = "A",
        effector_chain: str = "B",
    ):
        self._af3_output_dir = Path(af3_output_dir)
        self._ground_truth = Path(ground_truth)
        self._receptor_chain = receptor_chain
        self._effector_chain = effector_chain
        # Caches populated by _ensure_parsed.
        self._parsed = False
        self._ra_effs: List[float] = []
        self._iptms: List[float] = []
        self._cif_paths: List[Path] = []
        self._failures: List[str] = []

    # ── Aggregate properties ─────────────────────────────────────────

    @property
    def best_ra_eff(self) -> Optional[float]:
        self._ensure_parsed()
        return min(self._ra_effs) if self._ra_effs else None

    @property
    def mean_ra_eff(self) -> Optional[float]:
        self._ensure_parsed()
        return (sum(self._ra_effs) / len(self._ra_effs)
                if self._ra_effs else None)

    @property
    def best_iptm(self) -> Optional[float]:
        self._ensure_parsed()
        return max(self._iptms) if self._iptms else None

    @property
    def mean_iptm(self) -> Optional[float]:
        self._ensure_parsed()
        return (sum(self._iptms) / len(self._iptms)
                if self._iptms else None)

    @property
    def total_predictions(self) -> int:
        """Count of predictions that yielded a usable ra_eff."""
        self._ensure_parsed()
        return len(self._ra_effs)

    @property
    def failures(self) -> List[str]:
        """Failure tags from parsing, e.g.
        ['no_cif_files_found', 'per_prediction_rmsd_errors:2'].
        `to_summary_row` joins them with commas for the CSV."""
        self._ensure_parsed()
        return list(self._failures)

    @property
    def cif_paths(self) -> List[Path]:
        """All CIF paths that yielded a usable ra_eff, in the same order
        as the ra_eff values.  This is the only route from the
        aggregate back to the individual predictions."""
        self._ensure_parsed()
        return list(self._cif_paths)

    # ── Methods that take thresholds ─────────────────────────────────

    def n_correct_interface(self, thresholds) -> int:
        """Count of predictions whose ra_eff is below
        thresholds.orthogonal_af3_ra_max."""
        self._ensure_parsed()
        cutoff = thresholds.orthogonal_af3_ra_max
        return sum(1 for r in self._ra_effs if r < cutoff)

    def af3_disagrees(self, thresholds) -> bool:
        """True iff best_ra_eff is at or above
        thresholds.orthogonal_af3_ra_max.  Returns False when there is
        no data.  Diagnostic only.  AF3 is NOT in the orthogonal gate,
        which is Sc + BSA + ΔΔG."""
        br = self.best_ra_eff
        if br is None:
            return False  # no data → no disagreement signal
        return br >= thresholds.orthogonal_af3_ra_max

    # ── Output shape ─────────────────────────────────────────────────

    def to_summary_row(self, seq_name: str) -> Dict[str, str]:
        """The six `af3_nomsa_*` columns this class can fill on its own.

        bin/parse_af3_output.py writes these plus
        `af3_nomsa_n_correct_interface`, which it adds itself because
        that count needs a thresholds object."""
        return {
            "seq_name": seq_name,
            "af3_nomsa_best_ra_eff":
                f"{self.best_ra_eff:.3f}" if self.best_ra_eff is not None else "",
            "af3_nomsa_mean_ra_eff":
                f"{self.mean_ra_eff:.3f}" if self.mean_ra_eff is not None else "",
            "af3_nomsa_best_iptm":
                f"{self.best_iptm:.4f}" if self.best_iptm is not None else "",
            "af3_nomsa_mean_iptm":
                f"{self.mean_iptm:.4f}" if self.mean_iptm is not None else "",
            "af3_nomsa_total_predictions": str(self.total_predictions),
            "af3_nomsa_failures": ",".join(self.failures),
        }

    # ── Constructor that delays parsing ──────────────────────────────

    @classmethod
    def from_output_dir(
        cls,
        af3_output_dir: Path,
        ground_truth: Path,
        receptor_chain: str = "A",
        effector_chain: str = "B",
    ) -> "AF3ConfidenceAggregate":
        return cls(af3_output_dir, ground_truth, receptor_chain, effector_chain)

    # ── Parsing ──────────────────────────────────────────────────────

    def _ensure_parsed(self) -> None:
        if self._parsed:
            return
        self._parsed = True

        if not self._af3_output_dir.is_dir():
            self._failures.append(f"af3_output_dir_missing:{self._af3_output_dir}")
            return

        # Canonical AF3-v2 layout is seed-N_sample-M/model.cif.  Older
        # layouts fall back to a broader glob, which excludes the
        # top-level <name>_model.cif duplicate.
        cifs = sorted(
            self._af3_output_dir.rglob("seed-*_sample-*/model.cif")
        )
        if not cifs:
            cifs = [
                p for p in sorted(self._af3_output_dir.rglob("*.cif"))
                if p.name == "model.cif"
                or (p.name.endswith(".cif")
                    and not p.name.endswith("_model.cif"))
            ]
        if not cifs:
            self._failures.append("no_cif_files_found")
            return

        per_prediction_errors = 0
        for cif in cifs:
            ra_eff, rmsd_err = _receptor_aligned_effector_rmsd(
                cif, self._ground_truth,
                self._receptor_chain, self._effector_chain,
            )
            iptm, iptm_err = _read_confidence(cif)
            if ra_eff is None:
                per_prediction_errors += 1
                err_msg = rmsd_err or "(no error message returned)"
                print(
                    f"[{cif.parent.name}] RMSD failed: {err_msg}",
                    file=sys.stderr,
                )
                continue
            self._ra_effs.append(ra_eff)
            self._cif_paths.append(cif)
            if iptm is not None:
                self._iptms.append(iptm)
            elif iptm_err:
                print(
                    f"[{cif.parent.name}] iPTM read failed: {iptm_err}",
                    file=sys.stderr,
                )

        if per_prediction_errors:
            self._failures.append(
                f"per_prediction_rmsd_errors:{per_prediction_errors}"
            )
        if not self._ra_effs:
            self._failures.append("no_successful_rmsd")
