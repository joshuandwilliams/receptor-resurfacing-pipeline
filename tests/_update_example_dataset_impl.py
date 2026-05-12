#!/usr/bin/env python3
"""
_update_example_dataset_impl.py
--------------------------------
Diff-and-replace tool for tests/<module>/example_output_files/.

Reads tests/<module>/reference_manifest.txt (a positive include list of
glob patterns) and replaces example_output_files/ with the subset of
the freshly-run --updated-output-folder that matches the manifest.

Prints a diff (added / removed / changed paths) to stdout, then performs
the replace.  Designed to run non-interactively under SLURM — there is
no confirmation prompt.  See tests/update_example_dataset.slurm.sh for
the wrapper that submits this.

Manifest format
===============
- Plain text, one glob per line.
- Lines starting with `#` are comments; blank lines are ignored.
- Globs match paths RELATIVE to the updated output folder.  Patterns
  support `*`, `?`, `**` (recursive across path segments — pathlib
  semantics).
- Order does not matter.  Duplicate matches are deduplicated.
- A manifest line that matches ZERO files in the updated output folder
  is a hard error: the manifest claims a file the run did not produce
  (or did not produce in the expected location).  Fail fast.

Behaviour
=========
1. Read manifest → list of glob patterns.
2. For each pattern, glob it against updated_output_folder/.
3. Union all matched paths.  Verify every pattern matched ≥1 file
   (otherwise fail with a clear error).
4. Compare against current example_output_files/ (also recursively
   walked); compute added / removed / changed sets.  "Changed" means
   the file appears on both sides but content (or just size+mtime, by
   default) differs.
5. Print a one-line summary then per-set listings.
6. Replace example_output_files/ with the matched subset:
   - Write to example_output_files.new/
   - rm -rf example_output_files/
   - mv example_output_files.new/ → example_output_files/
7. Exit 0 on success.

Local --dry-run
===============
Pass --dry-run to print the diff and skip the replace.  Useful for
sanity-checking a new manifest against the current example_output_files/
without writing anything.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
from pathlib import Path
from typing import List, Set, Tuple


def _read_manifest(path: Path) -> List[str]:
    if not path.is_file():
        raise SystemExit(f"ERROR: manifest not found: {path}")
    patterns: List[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        patterns.append(line)
    if not patterns:
        raise SystemExit(f"ERROR: manifest has no patterns: {path}")
    return patterns


def _match_one(root: Path, pattern: str) -> Set[Path]:
    """Glob `pattern` (relative) against `root`, returning relative-Path matches."""
    # pathlib.Path.glob handles ** recursive; '*' matches a single segment.
    matched = set()
    for p in root.glob(pattern):
        if p.is_file():
            matched.add(p.relative_to(root))
    return matched


def _match_all(root: Path, patterns: List[str]) -> Tuple[Set[Path], List[str]]:
    """Return (union of all matches, list of patterns that matched ZERO files)."""
    union: Set[Path] = set()
    empty: List[str] = []
    for pat in patterns:
        hits = _match_one(root, pat)
        if not hits:
            empty.append(pat)
            continue
        union.update(hits)
    return union, empty


def _walk_files(root: Path) -> Set[Path]:
    """All files under `root`, as relative Paths.  Empty set if root absent."""
    if not root.is_dir():
        return set()
    out: Set[Path] = set()
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            out.add(Path(dirpath, name).relative_to(root))
    return out


def _hash_file(p: Path, chunk: int = 1 << 16) -> str:
    h = hashlib.sha1()
    with open(p, "rb") as f:
        while True:
            data = f.read(chunk)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def _classify(matched: Set[Path], current: Set[Path],
              src: Path, dst: Path,
              ) -> Tuple[Set[Path], Set[Path], Set[Path], Set[Path]]:
    """Return (added, removed, changed, unchanged) sets of relative paths."""
    added = matched - current
    removed = current - matched
    intersect = matched & current
    changed: Set[Path] = set()
    unchanged: Set[Path] = set()
    for rel in intersect:
        s = src / rel
        d = dst / rel
        try:
            if s.stat().st_size != d.stat().st_size:
                changed.add(rel)
                continue
        except OSError:
            changed.add(rel)
            continue
        # Same size — hash to confirm.  Hashing every file in the
        # 2500-file negsteer fixture is acceptable on HPC (single-digit
        # seconds for CSV-sized payloads).
        try:
            if _hash_file(s) != _hash_file(d):
                changed.add(rel)
            else:
                unchanged.add(rel)
        except OSError:
            changed.add(rel)
    return added, removed, changed, unchanged


def _print_set(label: str, paths: Set[Path], limit: int = 200) -> None:
    print(f"\n  {label}: {len(paths)}")
    if not paths:
        return
    for p in sorted(paths)[:limit]:
        print(f"    {p}")
    if len(paths) > limit:
        print(f"    ... ({len(paths) - limit} more — full list suppressed)")


def _replace_atomically(matched: Set[Path], src: Path, dst: Path) -> None:
    """Stage `matched` files under <dst>.new/, then atomic swap with <dst>."""
    staging = dst.parent / (dst.name + ".new")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True, exist_ok=True)
    for rel in sorted(matched):
        s = src / rel
        d = staging / rel
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(s, d)
    # Atomic swap: rm old, mv new
    if dst.exists():
        shutil.rmtree(dst)
    staging.rename(dst)


def main(argv: List[str] = None) -> int:
    p = argparse.ArgumentParser(
        description="Diff + replace tests/<module>/example_output_files/ "
                    "from a freshly-run output folder, gated by the "
                    "per-module reference_manifest.txt.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--module", required=True,
                   help="Module name (e.g. rfdiffusion, negative_steering).  "
                        "Manifest is read from tests/<module>/reference_manifest.txt.")
    p.add_argument("--updated-output-folder", required=True, type=Path,
                   help="Path to the freshly-run output folder, typically "
                        "tests/<module>/receptor_resurfacing_results/.")
    p.add_argument("--tests-root", type=Path, default=None,
                   help="Override tests/ root (default: parent of this script).")
    p.add_argument("--dry-run", action="store_true",
                   help="Print diff and exit; do not modify example_output_files/.")
    args = p.parse_args(argv)

    tests_root = (args.tests_root or Path(__file__).resolve().parent).resolve()
    module_dir = tests_root / args.module
    if not module_dir.is_dir():
        print(f"ERROR: module dir not found: {module_dir}", file=sys.stderr)
        return 2

    manifest_path = module_dir / "reference_manifest.txt"
    example_dir = module_dir / "example_output_files"
    src_dir = args.updated_output_folder.resolve()

    if not src_dir.is_dir():
        print(f"ERROR: --updated-output-folder not found: {src_dir}",
              file=sys.stderr)
        return 2

    print("=" * 60)
    print(f"update_example_dataset — module={args.module}")
    print(f"  manifest:           {manifest_path}")
    print(f"  source (new):       {src_dir}")
    print(f"  target (existing):  {example_dir}")
    print(f"  dry-run:            {args.dry_run}")
    print("=" * 60)

    patterns = _read_manifest(manifest_path)
    print(f"\nManifest patterns: {len(patterns)}")

    matched, empty_patterns = _match_all(src_dir, patterns)
    if empty_patterns:
        print("\nERROR: the following manifest patterns matched ZERO files "
              "in the updated output folder.  This typically means the test "
              "run did not produce the expected outputs, OR the manifest "
              "is out of date.  Fix one of those before re-running.",
              file=sys.stderr)
        for pat in empty_patterns:
            print(f"  - {pat}", file=sys.stderr)
        return 3

    print(f"Files matched by manifest: {len(matched)}")

    current = _walk_files(example_dir)
    print(f"Files in current example_output_files/: {len(current)}")

    added, removed, changed, unchanged = _classify(
        matched, current, src_dir, example_dir,
    )

    print("\nDiff summary:")
    print(f"  added (in new, not in existing):   {len(added)}")
    print(f"  removed (in existing, not in new): {len(removed)}")
    print(f"  changed (in both, content differs):{len(changed)}")
    print(f"  unchanged:                         {len(unchanged)}")

    _print_set("ADDED",   added)
    _print_set("REMOVED", removed)
    _print_set("CHANGED", changed)

    if args.dry_run:
        print("\n[--dry-run] no files modified.")
        return 0

    print(f"\nReplacing {example_dir} ...")
    _replace_atomically(matched, src_dir, example_dir)
    print(f"Done.  example_output_files/ now contains {len(matched)} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
