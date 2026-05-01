from __future__ import annotations

import re

# Nextflow workdir: …/work/AB/HASH where AB is two hex chars and HASH is the
# per-task hash (30 hex chars in current Nextflow versions; quantifier widened
# to {28,40} to absorb minor version variation).
_WORKDIR_RE = re.compile(
    r"(?:/[^/\s\"']+)*?/work/[0-9a-f]{2}/[0-9a-f]{28,40}"
)

# User-home prefixes on the platforms this pipeline runs on.
_HOME_RE = re.compile(
    r"/(?:hpc-home|Users|home)/[^/\s\"']+"
)


def canonicalize_workdir_paths(s: str) -> str:
    """Rewrite Nextflow workdir prefixes and user-home prefixes to canonical placeholders.

    Workdir-hash substitution must run before the home substitution so the
    longer, more specific match wins (see plan §4.4).
    """
    s = _WORKDIR_RE.sub("<WORKDIR>", s)
    s = _HOME_RE.sub("<HOME>", s)
    return s
