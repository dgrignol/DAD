"""
test_pipeline_config.py

Purpose
Validate the shared configuration infrastructure introduced to replace the
numbered 00/01 scripts.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_pipeline_config.py -q
"""

from __future__ import annotations

from pathlib import Path
import sys

import pytest


PREPROCESSING_DIR = (
    Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
)
sys.path.insert(0, str(PREPROCESSING_DIR))

import pipeline_config  # noqa: E402


def _minimal_config(tmp_path: Path) -> dict:
    """Build a valid core config with separate raw and derivative roots."""
    return {
        "project": {
            "root": str(tmp_path / "project"),
            "data_root": str(tmp_path / "raw"),
            "derivatives_root": str(tmp_path / "outputs"),
            "docs_root": str(tmp_path / "project" / "docs"),
            "params_root": str(tmp_path / "project" / "params"),
        },
        "layout": {
            "subject_bids_id": "sub-01",
            "session_bids_id": "ses-01",
            "raw_subject_folder_name": "sub01_date",
            "derivatives_preprocessing_template": (
                "preprocessing/{subject_bids_id}/{session_bids_id}"
            ),
            "derivatives_analysis_template": (
                "analysis/{subject_bids_id}/{session_bids_id}"
            ),
            "derivatives_reports_template": (
                "reports/{subject_bids_id}/{session_bids_id}"
            ),
        },
        "data_discovery": {
            "meg_fif_glob": "meg/sub01_run*.fif",
            "reported_trigger_csv_glob": "debug*.csv",
        },
        "runtime": {
            "overwrite_existing": False,
            "stop_on_error": True,
        },
    }


def test_derivative_templates_resolve_under_authoritative_root(tmp_path):
    """Derivative templates must be rooted at project.derivatives_root."""
    cfg = _minimal_config(tmp_path)
    resolved = pipeline_config.build_resolved_paths(cfg)

    expected = (
        tmp_path / "outputs" / "preprocessing" / "sub-01" / "ses-01"
    ).resolve()
    assert Path(resolved["derivatives_preprocessing_dir"]) == expected
    assert Path(resolved["subject_data_dir"]) == (
        tmp_path / "raw" / "sub01_date"
    ).resolve()


def test_raw_tree_output_guard_rejects_nested_output(tmp_path):
    """No pipeline output may be written inside the raw dataset tree."""
    raw_root = tmp_path / "raw"
    with pytest.raises(ValueError, match="read-only raw-data tree"):
        pipeline_config.assert_output_path_outside_raw(
            raw_root / "sub01" / "generated.json",
            raw_root,
        )


def test_stage_validation_requires_only_requested_sections(tmp_path):
    """Core config should not require ICA/epoching/export for an inventory run."""
    cfg = _minimal_config(tmp_path)
    errors, warnings = pipeline_config.validate_config(cfg)
    assert errors == []
    assert warnings == []

    errors, _warnings = pipeline_config.validate_config(
        cfg,
        required_sections=("trigger_parser",),
    )
    assert errors == ["`trigger_parser` must be a mapping."]
