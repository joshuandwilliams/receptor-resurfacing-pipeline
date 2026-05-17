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


# ── HADDOCK Session 7 restraint params ──────────────────────────────


@pytest.mark.local_unit
class TestHaddockContactPairsRegex:
    spec = next(s for s in vp.PARAM_SPECS if s.name == "haddock_contact_pairs")
    pat = spec.kwargs["pattern"]

    def test_empty_ok(self):
        assert vp._validate_regex("", {"pattern": self.pat}) is None

    def test_single_pair(self):
        assert vp._validate_regex("A25-C42", {"pattern": self.pat}) is None

    def test_multiple_pairs(self):
        assert vp._validate_regex("A25-C42 A13-C94", {"pattern": self.pat}) is None

    def test_bad_missing_dash(self):
        err = vp._validate_regex("A25 C42", {"pattern": self.pat})
        assert err and "does not match" in err

    def test_bad_missing_chain(self):
        err = vp._validate_regex("25-C42", {"pattern": self.pat})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestHaddockActiveResiduesRegex:
    spec = next(s for s in vp.PARAM_SPECS
                if s.name == "haddock_receptor_active_residues")
    pat = spec.kwargs["pattern"]

    def test_empty_ok(self):
        assert vp._validate_regex("", {"pattern": self.pat}) is None

    def test_singles(self):
        assert vp._validate_regex("25,35,40", {"pattern": self.pat}) is None

    def test_with_range(self):
        assert vp._validate_regex("25,40-44", {"pattern": self.pat}) is None

    def test_bad_letters(self):
        err = vp._validate_regex("A25,A30", {"pattern": self.pat})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestHaddockPairDistanceRegex:
    spec = next(s for s in vp.PARAM_SPECS if s.name == "haddock_pair_distance")
    pat = spec.kwargs["pattern"]

    def test_default(self):
        assert vp._validate_regex("2,2,4", {"pattern": self.pat}) is None

    def test_floats(self):
        assert vp._validate_regex("2.5,1.5,4.0", {"pattern": self.pat}) is None

    def test_wrong_count(self):
        err = vp._validate_regex("2,4", {"pattern": self.pat})
        assert err and "does not match" in err

    def test_no_plusminus(self):
        err = vp._validate_regex("2±2±4", {"pattern": self.pat})
        assert err and "does not match" in err


@pytest.mark.local_unit
class TestHaddockChosenClusterCustom:
    # Pull the lambda dynamically per call so Python doesn't bind it as
    # a method via the descriptor protocol when accessed on the class.
    @staticmethod
    def _fn():
        spec = next(s for s in vp.PARAM_SPECS
                    if s.name == "haddock_chosen_cluster")
        return spec.kwargs["fn"]

    def test_null_ok(self):
        fn = self._fn()
        assert fn(None) is None
        assert fn("null") is None

    def test_int_ok(self):
        assert self._fn()(3) is None

    def test_zero_rejected(self):
        assert self._fn()(0) is not None

    def test_string_rejected(self):
        assert self._fn()("3") is not None


@pytest.mark.local_unit
class TestBranchACrossParamRule:
    def test_no_branch_a_is_ok(self):
        # Branch B (receptor_input absent) — no rule applies.
        errors = vp.validate_params({"pdb_file": "/tmp/x.pdb"})
        assert errors == []

    def test_branch_a_with_no_restraints_fails(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
        })
        assert any("at least one HADDOCK restraint" in e for e in errors)

    def test_branch_a_with_pairs_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "haddock_contact_pairs": "A25-C42",
        })
        assert errors == []

    def test_branch_a_with_receptor_active_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "haddock_receptor_active_residues": "25,40-44",
        })
        assert errors == []

    def test_branch_a_with_only_effector_active_ok(self):
        # Post-commit-3 amendment: effector-only is valid (the "I want
        # this target face involved" mode).  Receptor design region for
        # clash bookkeeping comes from the contig at HADDOCK time.
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "haddock_effector_active_residues": "42,94",
        })
        assert errors == []

    def test_branch_a_with_pairs_and_effector_active_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "haddock_contact_pairs": "A73-B31 A72-B32 A71-B33",
            "haddock_effector_active_residues": "20-26",
        })
        assert errors == []


@pytest.mark.local_unit
class TestChosenClusterBranchAOnly:
    def test_branch_b_with_chosen_fails(self):
        errors = vp.validate_params({
            "pdb_file": "/tmp/x.pdb",
            "haddock_chosen_cluster": 2,
        })
        assert any("only applies in Branch A" in e for e in errors)

    def test_branch_a_with_chosen_ok(self):
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "haddock_contact_pairs": "A25-C42",
            "haddock_chosen_cluster": 2,
        })
        assert errors == []
