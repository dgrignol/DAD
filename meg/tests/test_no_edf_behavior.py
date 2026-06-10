"""
test_no_edf_behavior.py

Purpose
Validate the no-EDF preflight contract implemented in
`code/preprocessing_basic/dataset_inventory.py`.

Contract under test
- When eye model is enabled and configured with `skip_if_missing=true`,
  absence of EDF/eye files must not hard-fail the pipeline.
- The report payload must explicitly mark eye-model skip and include a reason.

Usage examples
1) Run only this test:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_no_edf_behavior.py -q

2) Run all tests:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests -q
"""

from __future__ import annotations

from pathlib import Path
import sys

import mne
import numpy as np


PREPROCESSING_DIR = (
    Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
)
sys.path.insert(0, str(PREPROCESSING_DIR))

import dataset_inventory  # noqa: E402


def test_inventory_marks_eye_model_skip_when_no_edf(tmp_path):
    """
    No-EDF case should be valid and explicitly reported as eye-model skip.

    The setup includes:
    - one valid FIF with expected stim channel,
    - one matching reported-trigger CSV,
    - no eye files in configured eye search roots.
    """
    subject_data_dir = tmp_path / "data" / "sub01"
    reports_csv_dir = tmp_path / "experiment" / "output_files" / "sub01"
    subject_data_dir.mkdir(parents=True, exist_ok=True)
    reports_csv_dir.mkdir(parents=True, exist_ok=True)

    # Minimal FIF that passes step-01 header checks.
    sfreq = 1000.0
    ch_names = ["MEG001", "STI101"]
    ch_types = ["mag", "stim"]
    data = np.zeros((2, 1000), dtype=float)
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
    raw = mne.io.RawArray(data, info, verbose="ERROR")
    fif_path = subject_data_dir / "SUB01_block01_meg.fif"
    raw.save(fif_path, overwrite=True, verbose="ERROR")

    # Minimal reported-trigger CSV so "has_reported_trigger_csv_files" is true.
    csv_path = (
        reports_csv_dir
        / "debug_actual_triggers_occlusion_v18_rescueTraject_sub01_block01_run01.csv"
    )
    csv_path.write_text("trigger,trial,frame,seconds,label\n31,1,1,0.0,trial_condition_onset\n")

    default_cfg_path = tmp_path / "pipeline_default.yaml"
    default_cfg_path.write_text("dummy: true\n")

    cfg = {
        "layout": {
            "subject_bids_id": "sub-01",
            "session_bids_id": "ses-01",
        },
        "resolved_paths": {
            "data_root": str(tmp_path / "data"),
            "subject_data_dir": str(subject_data_dir),
            "derivatives_root": str(tmp_path / "derivatives"),
        },
        "data_discovery": {
            "meg_fif_glob": "SUB01_block*_meg.fif",
            "reported_trigger_csv_glob": "debug_actual_triggers*_sub01_block*_run*.csv",
            "reported_trigger_search_roots": [str(reports_csv_dir)],
            "eye_search_glob_roots": [str(subject_data_dir)],
            "eye_search_extensions": [".edf", ".asc", ".tsv"],
        },
        "trigger_parser": {
            "stim_channel": "STI101",
        },
        "eye_model": {
            "enabled": True,
            "skip_if_missing": True,
            "skip_reason_when_missing": "No EDF/eye-tracker input found under configured eye search roots.",
        },
    }

    payload = dataset_inventory.build_inventory_payload(
        cfg=cfg,
        default_config_path=default_cfg_path,
        subject_config_path=None,
        warnings=[],
    )

    assert payload["status"] == "ok"
    assert payload["checks"]["eye_data_available"] is False
    assert payload["counts"]["eye_files"] == 0
    assert payload["eye_model_status"]["enabled"] is True
    assert payload["eye_model_status"]["skip_if_missing"] is True
    assert payload["eye_model_status"]["will_skip_due_to_missing_eye_data"] is True
    assert (
        payload["eye_model_status"]["skip_reason"]
        == "No EDF/eye-tracker input found under configured eye search roots."
    )
    assert payload["validation"]["is_valid"] is True
