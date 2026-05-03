"""Make the repo's bin/ directory importable from campaign scripts.

The production pipeline's Python modules live under bin/ as standalone
scripts that Nextflow invokes. They are not packaged, so bin/ is not on
sys.path by default and modules like compute_metrics or
boltz2_negative_steering cannot be imported directly from elsewhere in
the repo.

Campaign scripts under experiments/ should reuse those modules rather
than copying their logic. Duplication of production logic into ad-hoc
analysis scripts is a known pain point in this repo (see
notes/inventory/05_findings.md, "Surprising duplication") and is the
specific problem this helper exists to avoid.

Import this module once at the top of a campaign script:

    import experiments._path_setup  # noqa: F401  - adds bin/ to sys.path

After that, `from compute_metrics import find_contact_residues_heavy`
and similar imports from bin/ work normally.

The module also exposes REPO_ROOT (a pathlib.Path) for callers that
need the repo root for other purposes (locating notes/, params files,
etc.) without re-walking the filesystem.
"""

from __future__ import annotations

import sys
from pathlib import Path


def _discover_repo_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "main.nf").is_file() and (candidate / "bin").is_dir():
            return candidate
    raise RuntimeError(
        "Could not locate repo root: no ancestor of "
        f"{start} contains both 'main.nf' and 'bin/'. "
        "experiments/_path_setup.py expects to live inside the "
        "receptor-pipeline repo."
    )


REPO_ROOT: Path = _discover_repo_root(Path(__file__).resolve().parent)

_BIN_DIR = str(REPO_ROOT / "bin")
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)
