"""
test_epoch_export_smoke.py

Purpose
Run a small synthetic Chat-4 workflow through epoching and dRSA MAT export.
The test verifies run lineage, explicit trigger exclusions, epoch metadata,
stable MEG channel order, bad-channel NaNs, non-MEG omission, and manifest
updates without reading or writing participant data.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_epoch_export_smoke.py -q
"""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys

import mne
import numpy as np
from scipy.io import loadmat


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))


def _load_script_module(filename: str, module_name: str):
    """Load a numbered script as an importable module for direct smoke tests."""
    spec = importlib.util.spec_from_file_location(module_name, CODE_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


EPOCHING = _load_script_module("04_epoching.py", "dad_epoching_smoke")
EXPORT = _load_script_module(
    "06_export_neural_for_drsa.py",
    "dad_export_smoke",
)


def _event_row(
    *,
    code: int,
    sample: int,
    label: str,
    trial: int,
    source_name: str,
) -> dict[str, str]:
    """Build one complete canonical trigger CSV row."""
    return {
        "trigger_code": str(code),
        "event_label": label,
        "sample_start_raw": str(sample),
        "sample_end_raw": str(sample + 1),
        "raw_sample_frame": str(sample),
        "source_fif": source_name,
        "source_fif_path": f"/raw/{source_name}",
        "acquisition_run": "1",
        "is_transition_warning": "False",
        "transition_partner_code": "",
        "expected_trigger_code": str(code),
        "transition_decision_source": "",
        "transition_resolution_confidence": "",
        "unresolved_transition": "False",
        "behavioral_block": "1",
        "behavioral_run": "1",
        "behavioral_trial": str(trial),
        "behavioral_context_source": "reported_alignment",
        "replay_status": "original",
        "is_replay_trial": "False",
        "is_gaze_break": str(code == 150),
        "is_escape_termination": str(code == 152),
    }


def test_synthetic_epoching_and_drsa_export(tmp_path):
    """Tiny FIF data completes epoching/export with stable metadata contracts."""
    raw_root = tmp_path / "raw"
    subject_dir = raw_root / "sub01_date"
    raw_meg = subject_dir / "meg"
    input_files = subject_dir / "input_files"
    raw_meg.mkdir(parents=True)
    input_files.mkdir(parents=True)
    derivatives = tmp_path / "derivatives"
    preprocessing = derivatives / "preprocessing" / "sub-01" / "ses-01"
    analysis = derivatives / "analysis" / "sub-01" / "ses-01"
    reports = derivatives / "reports" / "sub-01" / "ses-01"
    sss = preprocessing / "sss"
    events_dir = preprocessing / "events"
    sss.mkdir(parents=True)
    events_dir.mkdir(parents=True)
    reports.mkdir(parents=True)

    source_name = "sub01_run01_block01_date.fif"
    source_path = raw_meg / source_name
    cleaned_path = sss / "sub01_run01_block01_date_ica.fif"
    info = mne.create_info(
        ["MEG0111", "MEG0112", "MISC001", "STI101"],
        sfreq=100.0,
        ch_types=["mag", "mag", "misc", "stim"],
    )
    rng = np.random.default_rng(42)
    data = rng.normal(scale=1e-13, size=(4, 900))
    data[2] = rng.normal(scale=1.0, size=900)
    data[3] = 0.0
    raw = mne.io.RawArray(data, info, verbose=False)
    raw.info["bads"] = ["MEG0112"]
    raw.save(source_path, overwrite=True, verbose=False)
    raw.save(cleaned_path, overwrite=True, verbose=False)

    observed_model = input_files / "MovDot_Sub01_Vtest_rescueTraject.mat"
    observed_model.touch()
    canonical_path = tmp_path / "canonical_meg.txt"
    canonical_path.write_text("MEG0111\nMEG0112\n", encoding="utf-8")

    event_rows = [
        _event_row(
            code=31,
            sample=100,
            label="trial_condition_onset_always_visible",
            trial=1,
            source_name=source_name,
        ),
        _event_row(
            code=1,
            sample=110,
            label="trial_sequence_identity",
            trial=1,
            source_name=source_name,
        ),
        _event_row(
            code=31,
            sample=300,
            label="trial_condition_onset_always_visible",
            trial=2,
            source_name=source_name,
        ),
        _event_row(
            code=2,
            sample=310,
            label="trial_sequence_identity",
            trial=2,
            source_name=source_name,
        ),
        _event_row(
            code=150,
            sample=320,
            label="gaze_break",
            trial=2,
            source_name=source_name,
        ),
        _event_row(
            code=31,
            sample=500,
            label="trial_condition_onset_always_visible",
            trial=3,
            source_name=source_name,
        ),
        _event_row(
            code=2,
            sample=510,
            label="trial_sequence_identity",
            trial=3,
            source_name=source_name,
        ),
    ]
    events_csv = events_dir / "02_trigger_events.csv"
    with events_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(event_rows[0]))
        writer.writeheader()
        writer.writerows(event_rows)

    manifest_path = reports / "source_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "derivative_lineage": {
                    "schema_version": "1.0",
                    "stages": {
                        "ica_manual": {
                            "processing_pass": "apply",
                            "status": "completed",
                            "qc": {"manual_exclude": []},
                            "records": [
                                {
                                    "source_path": str(source_path),
                                    "source_filename": source_name,
                                    "source_stem": source_path.stem,
                                    "acquisition_run_label": 1,
                                    "raw_relative_path": (
                                        f"sub01_date/meg/{source_name}"
                                    ),
                                    "input_path": str(source_path),
                                    "output_path": str(cleaned_path),
                                    "status": "completed",
                                }
                            ],
                        }
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    cfg = {
        "subject": {"id": "01"},
        "subject_profile": {
            "identity_range": [1, 2],
            "iti_sec_range": [0.5, 1.0],
            "versions": {"stimuli": "Vtest"},
        },
        "layout": {
            "subject_bids_id": "sub-01",
            "session_bids_id": "ses-01",
        },
        "resolved_paths": {
            "project_root": str(tmp_path),
            "data_root": str(raw_root),
            "subject_data_dir": str(subject_dir),
            "subject_input_files_dir": str(input_files),
            "derivatives_root": str(derivatives),
            "derivatives_preprocessing_dir": str(preprocessing),
            "derivatives_analysis_dir": str(analysis),
            "derivatives_reports_dir": str(reports),
        },
        "data_discovery": {
            "meg_fif_glob": "meg/sub01_run*.fif",
            "trajectory_mat_glob": (
                "input_files/MovDot_Sub01_V*_rescueTraject.mat"
            ),
        },
        "epoching": {
            "event_id_map": {"always_visible": 31},
            "tmin_sec": -0.5,
            "tmax_sec": 1.0,
            "baseline_sec": None,
            "detrend": 1,
            "reject_by_annotation": True,
            "amplitude_reject": {"mag": 1e-10},
            "trial_exclusion": {
                "gaze_break_codes": [150],
                "catch_type_1_codes": [116, 117],
                "catch_type_2_codes": [113, 118, 201, 202],
                "esc_termination_codes": [152],
                "exclude_gaze_break": True,
                "exclude_catch_type_1": True,
                "exclude_catch_type_2": True,
                "exclude_esc_termination": True,
                "exclude_replay_trials": False,
                "exclude_unresolved_transitions": True,
            },
        },
        "export": {
            "canonical_meg_channels_file": str(canonical_path),
            "tmin_sec": 0.0,
            "tmax_sec": 1.0,
            "sampling_rate_hz": 50.0,
            "missing_identity_policy": "error",
            "bad_or_missing_channel_policy": "nan",
        },
    }

    epoch_report = EPOCHING.run_epoching(
        cfg,
        source_manifest=manifest_path,
        overwrite=False,
    )
    assert epoch_report["counts"]["n_final_epochs"] == 2
    assert epoch_report["counts"]["trial_policy_exclusions"] == {
        "gaze_break": 1
    }
    assert epoch_report["baseline_correction_applied"] is False

    export_report = EXPORT.run_export(
        cfg,
        source_manifest=manifest_path,
        overwrite=False,
    )
    condition = export_report["conditions"][0]
    assert condition["tensor_shape"][:2] == [2, 2]
    mat = loadmat(condition["mat_path"], simplify_cells=True)
    assert mat["neural"].shape[:2] == (2, 2)
    assert np.isfinite(mat["neural"][:, 0, :]).all()
    assert np.isnan(mat["neural"][:, 1, :]).all()
    assert condition["ordered_channel_labels"] == ["MEG0111", "MEG0112"]
    assert "MISC001" not in condition["ordered_channel_labels"]
    assert "STI101" not in condition["ordered_channel_labels"]
    assert mat["expected_model_source_condition"] == "occluded_nondeviant"
    assert "predicted_model_path" not in mat
    assert export_report["trajectory_models"]["expected_source_condition"] == (
        "occluded_nondeviant"
    )
    assert "predicted" not in export_report["trajectory_models"]

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["derivative_lineage"]["stages"]["epoching"]["status"] == (
        "completed"
    )
    assert manifest["derivative_lineage"]["stages"]["export_drsa"]["status"] == (
        "completed"
    )
