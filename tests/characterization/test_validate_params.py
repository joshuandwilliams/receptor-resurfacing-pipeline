"""Unit tests for bin/validate_params.py.

These exercise each validator kind and the top-level validate_params()
entrypoint.  Marked local_unit so they run in the cheap Mac sweep.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

# Make bin/ importable so we can pull in the validator module directly.
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import validate_params as vp  # noqa: E402


# ── Per-validator kind ──────────────────────────────────────────────


@pytest.mark.local_unit
class TestIntRange:
    def test_in_range(self):
        assert vp._validate_int_range(5, {"min": 1, "max": 10}) is None

    def test_below_min(self):
        err = vp._validate_int_range(0, {"min": 1, "max": 10})
        assert err and "< min 1" in err

    def test_above_max(self):
        err = vp._validate_int_range(11, {"min": 1, "max": 10})
        assert err and "> max 10" in err

    def test_non_integer(self):
        err = vp._validate_int_range("abc", {"min": 1, "max": 10})
        assert err and "expected integer" in err

    def test_one_sided_min(self):
        assert vp._validate_int_range(99, {"min": 1, "max": None}) is None


@pytest.mark.local_unit
class TestFloatRange:
    def test_in_range(self):
        assert vp._validate_float_range(2.5, {"min": 0.0, "max": 5.0}) is None

    def test_nan_rejected(self):
        err = vp._validate_float_range(float("nan"), {"min": 0.0, "max": 5.0})
        assert err and "NaN" in err

    def test_string_coerced(self):
        assert vp._validate_float_range("2.5", {"min": 0.0, "max": 5.0}) is None

    def test_below_min(self):
        err = vp._validate_float_range(-0.1, {"min": 0.0, "max": 5.0})
        assert err and "< min 0.0" in err


@pytest.mark.local_unit
class TestChoice:
    def test_in_set(self):
        assert vp._validate_choice("mild", {"values": ["strong", "mild"]}) is None

    def test_not_in_set(self):
        err = vp._validate_choice("loud", {"values": ["strong", "mild"]})
        assert err and "loud" in err


@pytest.mark.local_unit
class TestBool:
    def test_true(self):
        assert vp._validate_bool(True, {}) is None

    def test_false(self):
        assert vp._validate_bool(False, {}) is None

    def test_int_rejected(self):
        # Python's bool is an int subclass, so isinstance(1, bool) is False
        # but bool(1) is True.  We want strict bool: int 1 should fail.
        err = vp._validate_bool(1, {})
        assert err and "expected bool" in err


@pytest.mark.local_unit
class TestNonEmptyStr:
    def test_normal(self):
        assert vp._validate_non_empty_str("hello", {}) is None

    def test_empty(self):
        err = vp._validate_non_empty_str("", {})
        assert err and "empty" in err

    def test_whitespace_only(self):
        err = vp._validate_non_empty_str("   ", {})
        assert err and "empty" in err


@pytest.mark.local_unit
class TestOptionalPath:
    def test_none(self):
        assert vp._validate_optional_path(None, {}) is None

    def test_path(self):
        assert vp._validate_optional_path("/tmp/x.pdb", {}) is None

    def test_empty(self):
        err = vp._validate_optional_path("", {})
        assert err and "empty" in err


@pytest.mark.local_unit
class TestRegex:
    def test_matches(self):
        assert vp._validate_regex("A", {"pattern": r"[A-Za-z]"}) is None

    def test_no_match(self):
        err = vp._validate_regex("AB", {"pattern": r"[A-Za-z]"})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestNumSeedsOdd:
    """Validator for the negsteer_num_seeds-must-be-odd rule (audit Q121)."""

    def test_three_ok(self):
        assert vp._check_num_seeds_odd(3) is None

    def test_one_ok(self):
        assert vp._check_num_seeds_odd(1) is None

    def test_two_rejected(self):
        err = vp._check_num_seeds_odd(2)
        assert err and "must be ODD" in err

    def test_four_rejected(self):
        err = vp._check_num_seeds_odd(4)
        assert err and "must be ODD" in err

    def test_zero_rejected_by_range(self):
        err = vp._check_num_seeds_odd(0)
        assert err and "< min 1" in err

    def test_huge_rejected_by_range(self):
        err = vp._check_num_seeds_odd(101)
        assert err and "> max 99" in err


# ── End-to-end ──────────────────────────────────────────────────────


@pytest.mark.local_unit
def test_validate_params_happy_path():
    params = {
        "project_name":            "test_project",
        "outdir":                  "/tmp/out",
        "num_designs":             10,
        "negsteer_num_seeds":      3,
        "negsteer_mode":           "mild",
        "sc_threshold":            0.62,
    }
    errors = vp.validate_params(params)
    assert errors == [], f"expected no errors, got {errors}"


@pytest.mark.local_unit
def test_validate_params_collects_all_errors():
    """All errors are collected — first failure doesn't short-circuit."""
    params = {
        "num_designs":        -1,            # below min
        "negsteer_num_seeds": 2,             # not odd
        "negsteer_mode":      "screaming",   # not in choices
        "sc_threshold":       1.5,           # above max
    }
    errors = vp.validate_params(params)
    assert len(errors) == 4
    assert any("num_designs" in e for e in errors)
    assert any("negsteer_num_seeds" in e for e in errors)
    assert any("negsteer_mode" in e for e in errors)
    assert any("sc_threshold" in e for e in errors)


@pytest.mark.local_unit
def test_validate_params_missing_keys_pass_through():
    """Missing keys are not errors — Nextflow defaults will fill in."""
    errors = vp.validate_params({})
    assert errors == []


@pytest.mark.local_unit
def test_validate_params_unknown_keys_warn_only(capsys):
    """Keys without a spec produce a NOTE on stderr but no error."""
    errors = vp.validate_params({"some_experimental_flag": 42})
    assert errors == []
    captured = capsys.readouterr()
    assert "some_experimental_flag" in captured.err
    assert "no spec entry" in captured.err


@pytest.mark.local_unit
def test_main_exit_codes(tmp_path):
    """End-to-end via main(): success returns 0, failure returns 1."""
    good = tmp_path / "good.json"
    good.write_text(json.dumps({"num_designs": 10, "negsteer_num_seeds": 3}))
    assert vp.main(["--params-json", str(good)]) == 0

    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"negsteer_num_seeds": 4}))
    assert vp.main(["--params-json", str(bad)]) == 1

    missing = tmp_path / "missing.json"
    assert vp.main(["--params-json", str(missing)]) == 2


@pytest.mark.local_unit
def test_coverage_report_runs():
    """--report path doesn't crash and surfaces coverage gaps."""
    text = vp.coverage_report()
    assert "Total spec entries" in text
    assert "any" in text  # the coverage-gap section is present


# ── Pose Solver restraint params ────────────────────────────────────


@pytest.mark.local_unit
class TestPoseSolverPairsRegex:
    spec = next(s for s in vp.PARAM_SPECS if s.name == "pose_solver_pairs")
    pat = spec.kwargs["pattern"]

    def test_empty_ok(self):
        assert vp._validate_regex("", {"pattern": self.pat}) is None

    def test_single_pair(self):
        assert vp._validate_regex("A25-C42", {"pattern": self.pat}) is None

    def test_multiple_pairs(self):
        assert vp._validate_regex("A25-C42 A13-C94", {"pattern": self.pat}) is None

    def test_with_dist_override(self):
        assert vp._validate_regex("A73-B48@4.5 A71-B50",
                                  {"pattern": self.pat}) is None

    def test_bad_missing_dash(self):
        err = vp._validate_regex("A25 C42", {"pattern": self.pat})
        assert err and "does not match" in err

    def test_bad_missing_chain(self):
        err = vp._validate_regex("25-C42", {"pattern": self.pat})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestPoseSolverExclusionsRegex:
    spec = next(s for s in vp.PARAM_SPECS if s.name == "pose_solver_exclusions")
    pat = spec.kwargs["pattern"]

    def test_empty_ok(self):
        assert vp._validate_regex("", {"pattern": self.pat}) is None

    def test_with_dist(self):
        assert vp._validate_regex("A8-B33@4.5", {"pattern": self.pat}) is None

    def test_multiple_with_dist(self):
        assert vp._validate_regex("A8-B33@4.5 A12-B17@6",
                                  {"pattern": self.pat}) is None

    def test_missing_dist_rejected(self):
        # @DIST is mandatory for exclusions.
        err = vp._validate_regex("A8-B33", {"pattern": self.pat})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestPoseSolverContigDesignRegionRegex:
    spec = next(s for s in vp.PARAM_SPECS
                if s.name == "pose_solver_contig_design_region")
    pat = spec.kwargs["pattern"]

    def test_empty_ok(self):
        assert vp._validate_regex("", {"pattern": self.pat}) is None

    def test_singles(self):
        assert vp._validate_regex("25,35,40", {"pattern": self.pat}) is None

    def test_with_range(self):
        assert vp._validate_regex("33-49,69-78", {"pattern": self.pat}) is None

    def test_bad_letters(self):
        err = vp._validate_regex("A25,A30", {"pattern": self.pat})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestBranchAPoseSolverPairsRequired:
    def test_branch_b_is_ok(self):
        # Branch B (receptor_input absent) — no rule applies.
        errors = vp.validate_params({"pdb_file": "/tmp/x.pdb"})
        assert errors == []

    def test_branch_a_with_no_pairs_fails(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
        })
        assert any("pose_solver_pairs" in e for e in errors)

    def test_branch_a_with_empty_pairs_fails(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "   ",  # whitespace-only also rejected
        })
        assert any("pose_solver_pairs" in e for e in errors)

    def test_branch_a_with_pairs_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A25-C42",
        })
        assert errors == []

    def test_branch_a_with_multiple_pairs_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A73-B31 A72-B32 A71-B33",
        })
        assert errors == []

    def test_branch_a_with_dist_overrides_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A73-B48@4.5 A71-B50 A8-B39",
        })
        assert errors == []


@pytest.mark.local_unit
class TestPoseSolverPairDistanceBoundsOrdered:
    def test_default_ordering_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A25-C42",
            "pose_solver_min_pair_distance": 3.5,
            "pose_solver_max_pair_distance": 6.0,
        })
        assert errors == []

    def test_min_equal_max_fails(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A25-C42",
            "pose_solver_min_pair_distance": 4.0,
            "pose_solver_max_pair_distance": 4.0,
        })
        assert any("min_pair_distance" in e and "max_pair_distance" in e
                   for e in errors)

    def test_min_greater_than_max_fails(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A25-C42",
            "pose_solver_min_pair_distance": 8.0,
            "pose_solver_max_pair_distance": 4.0,
        })
        assert any("min_pair_distance" in e and "max_pair_distance" in e
                   for e in errors)
