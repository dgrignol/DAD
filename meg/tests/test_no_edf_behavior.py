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
    meg_dir = subject_data_dir / "meg"
    input_dir = subject_data_dir / "input_files"
    output_dir = subject_data_dir / "output_files"
    meg_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Minimal FIF that passes step-01 header checks.
    sfreq = 1000.0
    ch_names = ["MEG001", "STI101"]
    ch_types = ["mag", "stim"]
    data = np.zeros((2, 1000), dtype=float)
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
    raw = mne.io.RawArray(data, info, verbose="ERROR")
    fif_path = meg_dir / "sub01_run01_150426.fif"
    raw.save(fif_path, overwrite=True, verbose="ERROR")

    # Required saved scientific inputs are discovered by explicit category.
    (input_dir / "MovDot_Sub01_V28_rescueTraject.mat").write_bytes(b"MAT")
    (input_dir / "Sub01_TrialStruct_v21_rescueTraject.mat").write_bytes(b"MAT")
    (output_dir / "MoveDot1_occlusion_v18_rescueTraject_SUB01_BLOCK01.mat").write_bytes(
        b"MAT"
    )
    csv_path = (
        output_dir
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
        "subject_profile": {
            "ocular_modalities": {
                "saved_edf": False,
            },
        },
        "data_discovery": {
            "subject_folder_glob_template": "sub{subject_id}_*",
            "meg_fif_glob": "meg/sub01_run*.fif",
            "trajectory_mat_glob": (
                "input_files/MovDot_Sub01_V*_rescueTraject.mat"
            ),
            "trial_struct_mat_glob": (
                "input_files/Sub01_TrialStruct_v*_rescueTraject.mat"
            ),
            "behavioral_mat_glob": (
                "output_files/MoveDot1_occlusion_v*_rescueTraject_SUB01_BLOCK*.mat"
            ),
            "reported_trigger_csv_glob": (
                "output_files/debug_actual_triggers*_sub01_block*_run*.csv"
            ),
            "eyetracker_edf_glob": "eyetracker/*.edf",
            "empty_room_fif_glob": "empty_room_150426.fif",
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
