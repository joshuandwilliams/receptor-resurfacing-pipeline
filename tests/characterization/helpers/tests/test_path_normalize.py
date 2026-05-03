from __future__ import annotations

import pytest

from tests.characterization.helpers.path_normalize import canonicalize_workdir_paths


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_workdir_hash_replaced():
    """Standalone /work/AB/HASH/... is rewritten to <WORKDIR>/..."""
    s = "/work/ab/0123456789abcdef0123456789abcd/file.txt"
    assert canonicalize_workdir_paths(s) == "<WORKDIR>/file.txt"


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_workdir_with_hpc_home_prefix():
    """A workdir path nested under /hpc-home/<user>/... has the hpc-home prefix absorbed."""
    s = "/hpc-home/jowillia/receptor_design/x/work/50/97ffcec2913113c9106d53dfe9d00b/sub"
    assert canonicalize_workdir_paths(s) == "<WORKDIR>/sub"


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_hpc_home_replaced():
    """/hpc-home/<user>/... outside a workdir becomes <HOME>/..."""
    s = "/hpc-home/jowillia/notes/foo.md"
    assert canonicalize_workdir_paths(s) == "<HOME>/notes/foo.md"


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_macos_home_replaced():
    """/Users/<user>/... becomes <HOME>/... (Mac home prefix)."""
    s = "/Users/jowillia/dev/receptor"
    assert canonicalize_workdir_paths(s) == "<HOME>/dev/receptor"


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_linux_home_replaced():
    """/home/<user>/... becomes <HOME>/... (typical Linux home prefix)."""
    s = "/home/jowillia/dev/receptor"
    assert canonicalize_workdir_paths(s) == "<HOME>/dev/receptor"


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_idempotent():
    """Applying the function twice equals applying it once."""
    s = "/hpc-home/jowillia/x/work/50/97ffcec2913113c9106d53dfe9d00b/sub"
    once = canonicalize_workdir_paths(s)
    twice = canonicalize_workdir_paths(once)
    assert once == twice


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_non_path_strings_pass_through():
    """Strings that don't contain workdir or home prefixes must come back unchanged."""
    for s in [
        "A50/0-50",
        "design_42_seq_3",
        "1234567890",
        "hello world",
        "RFDIFFUSION_FILTER",
        "motif_rmsd",
    ]:
        assert canonicalize_workdir_paths(s) == s


@pytest.mark.local_unit
def test_canonicalize_workdir_paths_short_hash_not_matched():
    """A workdir-shaped path with a hash shorter than the regex floor must NOT match."""
    s = "/work/ab/0123456789abcdef0123/sub"  # 20 hex chars; below the {28,40} floor
    assert canonicalize_workdir_paths(s) == s
