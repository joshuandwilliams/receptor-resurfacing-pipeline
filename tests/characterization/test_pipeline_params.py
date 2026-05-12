"""Unit tests for bin/pipeline_params.py (Phase 4 Tier 0.2)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import pipeline_params as pp  # noqa: E402


# Minimal valid params set that satisfies the specs we test against.
_BASE = {
    "project_name": "test",
    "num_designs": 10,
    "negsteer_num_seeds": 3,
    "negsteer_mode": "mild",
    "sc_threshold": 0.62,
}


@pytest.mark.local_unit
class TestFromDictHappyPath:
    def test_constructs_with_valid_dict(self):
        params = pp.PipelineParams.from_dict(_BASE)
        assert isinstance(params, pp.PipelineParams)

    def test_values_accessible_via_brackets(self):
        params = pp.PipelineParams.from_dict(_BASE)
        assert params["num_designs"] == 10

    def test_values_accessible_via_get(self):
        params = pp.PipelineParams.from_dict(_BASE)
        assert params.get("num_designs") == 10
        assert params.get("absent_key", "default") == "default"

    def test_contains(self):
        params = pp.PipelineParams.from_dict(_BASE)
        assert "num_designs" in params
        assert "absent_key" not in params

    def test_as_dict_round_trip(self):
        params = pp.PipelineParams.from_dict(_BASE)
        out = params.as_dict()
        assert out == _BASE
        # And it's a copy, not the internal mapping
        out["num_designs"] = 999
        assert params["num_designs"] == 10


@pytest.mark.local_unit
class TestFromDictValidation:
    def test_single_invalid_raises(self):
        bad = {**_BASE, "negsteer_num_seeds": 4}  # not odd
        with pytest.raises(pp.ParamValidationError) as exc_info:
            pp.PipelineParams.from_dict(bad)
        assert "negsteer_num_seeds" in str(exc_info.value)

    def test_multiple_errors_all_collected(self):
        bad = {**_BASE,
               "negsteer_num_seeds": 4,         # not odd
               "num_designs": -1,               # below min
               "negsteer_mode": "screaming"}    # not in choices
        with pytest.raises(pp.ParamValidationError) as exc_info:
            pp.PipelineParams.from_dict(bad)
        msg = str(exc_info.value)
        assert "negsteer_num_seeds" in msg
        assert "num_designs" in msg
        assert "negsteer_mode" in msg
        assert exc_info.value.errors  # list is populated

    def test_validation_error_message_includes_count(self):
        bad = {**_BASE,
               "negsteer_num_seeds": 4,
               "num_designs": -1}
        with pytest.raises(pp.ParamValidationError) as exc_info:
            pp.PipelineParams.from_dict(bad)
        assert "2 error(s)" in str(exc_info.value)


@pytest.mark.local_unit
class TestFromNextflowJson:
    def test_reads_valid_json(self, tmp_path):
        path = tmp_path / "params.json"
        path.write_text(json.dumps(_BASE))
        params = pp.PipelineParams.from_nextflow_json(path)
        assert params["num_designs"] == 10

    def test_missing_file_raises_filenotfound(self, tmp_path):
        path = tmp_path / "does_not_exist.json"
        with pytest.raises(FileNotFoundError):
            pp.PipelineParams.from_nextflow_json(path)

    def test_malformed_json_raises_valueerror(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text("{not valid json")
        with pytest.raises(ValueError):
            pp.PipelineParams.from_nextflow_json(path)

    def test_non_object_json_raises_typeerror(self, tmp_path):
        path = tmp_path / "list.json"
        path.write_text("[1, 2, 3]")
        with pytest.raises(TypeError):
            pp.PipelineParams.from_nextflow_json(path)

    def test_invalid_values_in_json_raise_validation_error(self, tmp_path):
        bad = {**_BASE, "negsteer_num_seeds": 4}
        path = tmp_path / "params.json"
        path.write_text(json.dumps(bad))
        with pytest.raises(pp.ParamValidationError):
            pp.PipelineParams.from_nextflow_json(path)


@pytest.mark.local_unit
class TestImmutability:
    def test_frozen_dataclass(self):
        """Cannot reassign the values field after construction."""
        params = pp.PipelineParams.from_dict(_BASE)
        with pytest.raises((AttributeError, Exception)):
            params.values = {}  # type: ignore


@pytest.mark.local_unit
class TestRevalidate:
    def test_validate_returns_empty_on_good_params(self):
        params = pp.PipelineParams.from_dict(_BASE)
        assert params.validate() == []
