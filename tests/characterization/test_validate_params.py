"""Unit tests for bin/validate_params.py.

These exercise each validator kind and the top-level validate_params()
entrypoint.  Marked logic so they run in the cheap Mac sweep.
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


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
class TestChoice:
    def test_in_set(self):
        assert vp._validate_choice("mild", {"values": ["strong", "mild"]}) is None

    def test_not_in_set(self):
        err = vp._validate_choice("loud", {"values": ["strong", "mild"]})
        assert err and "loud" in err


@pytest.mark.logic
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


@pytest.mark.logic
class TestNonEmptyStr:
    def test_normal(self):
        assert vp._validate_non_empty_str("hello", {}) is None

    def test_empty(self):
        err = vp._validate_non_empty_str("", {})
        assert err and "empty" in err

    def test_whitespace_only(self):
        err = vp._validate_non_empty_str("   ", {})
        assert err and "empty" in err


@pytest.mark.logic
class TestOptionalPath:
    def test_none(self):
        assert vp._validate_optional_path(None, {}) is None

    def test_path(self):
        assert vp._validate_optional_path("/tmp/x.pdb", {}) is None

    def test_empty(self):
        err = vp._validate_optional_path("", {})
        assert err and "empty" in err


@pytest.mark.logic
class TestRegex:
    def test_matches(self):
        assert vp._validate_regex("A", {"pattern": r"[A-Za-z]"}) is None

    def test_no_match(self):
        err = vp._validate_regex("AB", {"pattern": r"[A-Za-z]"})
        assert err and "does not match" in err


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
def test_validate_params_missing_keys_pass_through():
    """Missing keys are not errors — Nextflow defaults will fill in."""
    errors = vp.validate_params({})
    assert errors == []


@pytest.mark.logic
def test_validate_params_unknown_keys_warn_only(capsys):
    """Keys without a spec produce a NOTE on stderr but no error."""
    errors = vp.validate_params({"some_experimental_flag": 42})
    assert errors == []
    captured = capsys.readouterr()
    assert "some_experimental_flag" in captured.err
    assert "no spec entry" in captured.err


@pytest.mark.logic
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


@pytest.mark.logic
def test_coverage_report_runs():
    """--report path doesn't crash and surfaces coverage gaps."""
    text = vp.coverage_report()
    assert "Total spec entries" in text
    assert "any" in text  # the coverage-gap section is present


# ── Pose Solver restraint params ────────────────────────────────────


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
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


@pytest.mark.logic
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

    def test_unusable_bounds_are_left_to_the_per_param_spec(self):
        """The ordering rule cannot compare two values it cannot read.
        Reporting "min must be less than max" on top of the two type errors
        would be a third message about the same mistake, so the rule stands
        down and only the float_range errors are raised.
        """
        errors = vp.validate_params({
            "receptor_input": "/tmp/r.pdb",
            "effector_input": "/tmp/e.pdb",
            "pose_solver_pairs": "A25-C42",
            "pose_solver_min_pair_distance": "close",
            "pose_solver_max_pair_distance": "far",
        })
        assert len(errors) == 2
        assert all("expected number" in e for e in errors)


# ── A value of the wrong type, kind by kind ─────────────────────────
#
# Nextflow serialises params with JsonOutput.toJson, so a param can
# arrive as null, as a number where a string was meant, or as a list.
# Every validator has to name the type it got rather than raise, because
# validate_params reports the whole list at once and a raise would hide
# every error after the first.


@pytest.mark.logic
@pytest.mark.parametrize("value", [None, [1.0], {"min": 1}])
def test_a_float_range_reports_an_uncoercible_value(value):
    err = vp._validate_float_range(value, {"min": 0.0, "max": 1.0})
    assert err and "expected number" in err
    assert type(value).__name__ in err


@pytest.mark.logic
def test_a_non_empty_str_reports_a_non_string():
    err = vp._validate_non_empty_str(42, {})
    assert err == "expected string, got int"


@pytest.mark.logic
def test_an_optional_path_reports_a_value_that_is_neither_none_nor_a_string():
    """None is the "not supplied" case and passes.  A number is a mistake, and
    it must not be coerced into a path."""
    err = vp._validate_optional_path(42, {})
    assert err == "expected None or string path, got int"


@pytest.mark.logic
def test_a_regex_reports_a_non_string():
    """re.fullmatch would raise TypeError on a non-string, which would abort
    the whole validation run partway through its error list."""
    err = vp._validate_regex(42, {"pattern": r"[A-Za-z]"})
    assert err == "expected string, got int"


# ── The `any` escape hatch ──────────────────────────────────────────


@pytest.mark.logic
def test_the_any_validator_accepts_whatever_it_is_given():
    """`any` is the marker for a param whose valid range the threshold audit
    has not decided.  It must accept every shape, because an unspecced param
    does not block the pipeline."""
    for value in ("A1-5/10-10/A16-20 B", 7, None, [1, 2, 3], {"a": 1}):
        assert vp._validate_any(value, {}) is None


@pytest.mark.logic
def test_a_param_specced_as_any_produces_no_error():
    """contigs carries a free-form RFDiffusion contig string, listed as a
    coverage gap rather than left unknown, so the coverage report surfaces
    it."""
    assert vp.validate_params({"contigs": "A1-5/10-10/A16-20 B"}) == []
    assert "contigs" in vp.coverage_report()


# ── main(): the argument and file-shape paths ───────────────────────
#
# main() is invoked from the top of workflow MAIN_PIPELINE before any
# process is dispatched, and its exit code is the whole gate.  Each
# non-zero code has to be reachable for the right reason.


@pytest.mark.logic
def test_report_prints_the_coverage_report_and_exits_zero(capsys):
    assert vp.main(["--report"]) == 0
    out = capsys.readouterr().out
    assert "Total spec entries" in out
    assert "Coverage gaps" in out


@pytest.mark.logic
def test_neither_report_nor_params_json_is_an_argparse_error():
    """argparse exits 2 for a usage error, which is distinct from the 2 that
    a missing or malformed params file returns.  Both are non-zero, so the
    pipeline is blocked either way."""
    with pytest.raises(SystemExit) as exc:
        vp.main([])
    assert exc.value.code == 2


@pytest.mark.logic
def test_malformed_params_json_returns_two(tmp_path, capsys):
    """Nextflow writes this file itself, so malformed JSON means the harness
    is broken rather than a param being wrong.  It must not read as "no
    params set", which would pass with an empty error list."""
    bad = tmp_path / "params.json"
    bad.write_text('{"num_designs": 10,')
    assert vp.main(["--params-json", str(bad)]) == 2
    assert "malformed params JSON" in capsys.readouterr().err


@pytest.mark.logic
@pytest.mark.parametrize("payload", ["[1, 2, 3]", '"a string"', "null"])
def test_params_json_that_is_not_an_object_returns_two(tmp_path, payload, capsys):
    """JsonOutput.toJson(params) always yields a flat object.  Anything else
    would make `k in params` a membership test on the wrong container, and
    every spec would silently be skipped."""
    f = tmp_path / "params.json"
    f.write_text(payload)
    assert vp.main(["--params-json", str(f)]) == 2
    assert "must be an object" in capsys.readouterr().err


@pytest.mark.logic
def test_a_valid_run_reports_the_spec_and_param_counts(tmp_path, capsys):
    """The success line is the only record in the Nextflow log that validation
    happened at all, so it carries both counts."""
    f = tmp_path / "params.json"
    f.write_text(json.dumps({"num_designs": 10, "negsteer_num_seeds": 3}))
    assert vp.main(["--params-json", str(f)]) == 0
    err = capsys.readouterr().err
    assert f"{len(vp.PARAM_SPECS)} specs checked" in err
    assert "2 params provided" in err
