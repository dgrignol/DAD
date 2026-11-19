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


def test_working_cohort_excludes_p15_without_changing_accepted_accounting_roster():
    """Current analysis defaults and immutable Stage 2/3 coverage remain distinct."""
    import yaml
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] / "params/pipeline_default.yaml").read_text())
    assert cfg["analysis_cohort"]["included_subjects"] == list(pipeline_config.WORKING_SUBJECT_IDS)
    assert len(pipeline_config.WORKING_SUBJECT_IDS) == 35 and "15" not in pipeline_config.WORKING_SUBJECT_IDS
    assert len(pipeline_config.ANALYSIS_SUBJECT_IDS) == 36 and "15" in pipeline_config.ANALYSIS_SUBJECT_IDS
    assert "15" in pipeline_config.SUPPORTED_SUBJECT_IDS and "15" in cfg["subject_profiles"]
    assert cfg["matlab_drsa"]["group_average_subjects"] == list(pipeline_config.WORKING_SUBJECT_IDS)


@pytest.mark.parametrize("cutoff", [None, .01])
def test_final_highpass_is_optional_and_independent_of_fit_configuration(tmp_path, cutoff):
    """Null is a supported final setting, while fit settings remain separate."""
    cfg = _minimal_config(tmp_path)
    cfg["preprocessing"] = {"highpass_hz": cutoff, "lowpass_hz": None, "resample_sfreq_hz": None}
    errors, _ = pipeline_config.validate_config(cfg)
    assert not errors


@pytest.mark.parametrize("cutoff", [0, -1, float("nan"), float("inf"), True, "none"])
def test_invalid_final_highpass_cannot_masquerade_as_omitted(tmp_path, cutoff):
    cfg = _minimal_config(tmp_path)
    cfg["preprocessing"] = {"highpass_hz": cutoff}
    errors, _ = pipeline_config.validate_config(cfg)
    assert any("preprocessing.highpass_hz" in error for error in errors)


def _minimal_config(tmp_path: Path) -> dict:
    """
    Build a valid multisubject core config with isolated raw/output roots.

    The helper includes all supported profile keys because the default
    configuration contract validates the complete supported cohort even when a
    test selects only one participant.
    """
    profiles = {}
    for subject_id in pipeline_config.SUPPORTED_SUBJECT_IDS:
        profiles[subject_id] = {
            "external_folder": {
                "expected_name": f"sub{subject_id}_date",
            },
            "versions": {
                "stimuli": "Vtest",
                "schedule": "vtest",
                "runtime": "vtest",
            },
            "condition_triggers": {
                "always_visible": 31,
                "occluded_nondeviant": 41,
                "occluded_deviant": 51,
            },
            "identity_range": [1, 20],
            "reported_trigger_csv_pattern": (
                "output_files/debug_actual_triggers_occlusion_vtest_"
                "rescueTraject_sub{subject_id}_block*_run*.csv"
            ),
            "trigger_policy": {
                "required_subject_codes": [31, 41, 51],
                "allowed_fixed_codes": [111, 150, 151, 152],
                "optional_event_codes": [150, 151, 152],
                "gaze_break_code": 150,
                "replay_start_code": 151,
                "esc_termination_code": 152,
            },
            "ocular_modalities": {
                "eog_channels_in_fif": [],
                "saved_edf": False,
            },
            "completed_behavioral_blocks": [1],
            "known_exceptions": [],
        }

    return {
        "project": {
            "root": str(tmp_path / "project"),
            "data_root": str(tmp_path / "raw"),
            "derivatives_root": str(tmp_path / "outputs"),
            "docs_root": str(tmp_path / "project" / "docs"),
            "params_root": str(tmp_path / "project" / "params"),
        },
        "subject": {
            "id": "01",
        },
        "subject_profiles": profiles,
        "layout": {
            "subject_bids_id": "sub-01",
            "session_bids_id": "ses-01",
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
            "subject_folder_glob_template": "sub{subject_id}_*",
            "task_meg_fif_glob_template": "meg/sub{subject_id}_run*.fif",
            "trajectory_mat_glob_template": (
                "input_files/MovDot_Sub{subject_id}_V*_rescueTraject.mat"
            ),
            "trial_struct_mat_glob_template": (
                "input_files/Sub{subject_id}_TrialStruct_v*_rescueTraject.mat"
            ),
            "behavioral_mat_glob_template": (
                "output_files/MoveDot*_SUB{subject_id}_BLOCK*.mat"
            ),
            "reported_trigger_csv_glob_template": (
                "output_files/debug*_sub{subject_id}_block*_run*.csv"
            ),
            "eyetracker_edf_glob_template": "eyetracker/*.edf",
            "empty_room_fif_glob_template": "empty_room_{subject_folder_date}.fif",
            "meg_fif_glob": "meg/sub01_run*.fif",
            "trajectory_mat_glob": (
                "input_files/MovDot_Sub01_V*_rescueTraject.mat"
            ),
            "trial_struct_mat_glob": (
                "input_files/Sub01_TrialStruct_v*_rescueTraject.mat"
            ),
            "behavioral_mat_glob": "output_files/MoveDot*_SUB01_BLOCK*.mat",
            "reported_trigger_csv_glob": "debug*.csv",
            "eyetracker_edf_glob": "eyetracker/*.edf",
        },
        "runtime": {
            "overwrite_existing": False,
            "stop_on_error": True,
        },
    }


def test_derivative_templates_resolve_under_authoritative_root(tmp_path):
    """Derivative templates must be rooted at project.derivatives_root."""
    cfg = _minimal_config(tmp_path)
    (tmp_path / "raw" / "sub01_date").mkdir(parents=True)
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


def test_matlab_drsa_validation_accepts_explicit_pe_model_sources(tmp_path):
    """Shared config validation accepts registered PE model-source names."""
    cfg = _minimal_config(tmp_path)
    cfg["matlab_drsa"] = {
        "enabled": True,
        "matlab_bin": "/Applications/MATLAB_R2020a.app/bin/matlab",
        "use_arch_x86_64": False,
        "analyses": ["alt"],
        "drsa_type": "corr",
        "conditions_to_run": ["occluded_deviant"],
        "sensor_groups_to_run": ["occipital"],
        "model_source": "pe_prior_constrained_50",
        "model_sources": [
            "prediction_error",
            "pe_prior_constrained_50",
            "pe_belief_update",
        ],
        "pcr_strategy": "baseline_pcr_border",
        "smooth_neural": True,
        "smooth_window_samples": 15,
        "smoothing_plots": False,
        "plot_smoothing_window": 9,
        "group_average_subjects": ["01"],
        "occlusion_timing": {
            "frame_frequency_hz": 120.0,
            "first_complete_occlusion_frame": 130,
            "reappearance_frame": 191,
        },
    }

    errors, _warnings = pipeline_config.validate_config(cfg)

    assert errors == []


def test_site_config_precedes_subject_and_tuning_layers(tmp_path):
    """Cluster paths merge separately while tuning retains final precedence."""
    cfg = _minimal_config(tmp_path)
    default_raw = tmp_path / "default_raw"
    site_raw = tmp_path / "cluster_raw"
    site_derivatives = tmp_path / "cluster_derivatives"
    (site_raw / "sub01_date").mkdir(parents=True)
    cfg["project"]["data_root"] = str(default_raw)

    import yaml

    default_path = tmp_path / "default.yaml"
    site_path = tmp_path / "site.yaml"
    subject_path = tmp_path / "subject.yaml"
    tuning_path = tmp_path / "tuning.yaml"
    default_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    site_path.write_text(
        yaml.safe_dump(
            {
                "project": {
                    "data_root": str(site_raw),
                    "derivatives_root": str(site_derivatives),
                },
                "runtime": {"overwrite_existing": True},
            }
        ),
        encoding="utf-8",
    )
    subject_path.write_text(
        yaml.safe_dump({"subject": {"id": "01"}}),
        encoding="utf-8",
    )
    tuning_path.write_text(
        yaml.safe_dump({"runtime": {"overwrite_existing": False}}),
        encoding="utf-8",
    )

    resolved, _warnings, _default, _subject = (
        pipeline_config.load_and_resolve_config(
            default_path,
            subject_path,
            tuning_path,
            site_config_path=site_path,
        )
    )
    assert resolved["resolved_paths"]["data_root"] == str(site_raw.resolve())
    assert resolved["resolved_paths"]["derivatives_root"] == str(
        site_derivatives.resolve()
    )
    assert resolved["runtime"]["overwrite_existing"] is False
    assert resolved["_config_sources"]["site"] == str(site_path.resolve())


@pytest.mark.parametrize(
    ("selector", "expected"),
    [
        (1, "01"),
        ("01", "01"),
        ("sub01", "01"),
        ("sub-07", "07"),
    ],
)
def test_subject_selector_normalization(selector, expected):
    """Supported selector spellings resolve to the canonical two-digit ID."""
    assert pipeline_config.normalize_subject_id(selector) == expected


def test_subject_selector_rejects_out_of_range_values():
    """The external dataset currently supports participants 01 through 39."""
    with pytest.raises(ValueError, match="01-39"):
        pipeline_config.normalize_subject_id("40")


@pytest.mark.parametrize("field,value", [("n_components", 11), ("n_components", True),
    ("threshold_std", 0), ("threshold_std", float("nan")), ("min_peak_distance_sec", -1),
    ("tmin_sec", .25), ("tmax_sec", 0), ("long_peak_width_sec", float("inf"))])
def test_blink_aid_settings_reject_invalid_values_without_changing_final_ica(tmp_path, field, value):
    """The supplementary cap/threshold block is independent of final ICA tuning."""
    cfg = _minimal_config(tmp_path)
    cfg["artifact_annotation"] = {"enabled": True, "muscle_threshold": 10., "muscle_threshold_by_run": {},
        "muscle_channel_type": "mag", "muscle_min_length_good_sec": .2, "muscle_filter_freq_hz": [110., 140.],
        "eog_mode": "disabled", "eog_channel_names": [], "blink_channel": None,
        "blink_auto": {"n_components": 10, "threshold_std": 5., "min_peak_distance_sec": .25,
            "tmin_sec": -.25, "tmax_sec": .25, "long_peak_width_sec": .5}}
    cfg["ica"] = {"n_components": None}
    errors, _ = pipeline_config.validate_config(cfg, required_sections=("artifact_annotation",))
    assert errors == []
    cfg["artifact_annotation"]["blink_auto"][field] = value
    errors, _ = pipeline_config.validate_config(cfg, required_sections=("artifact_annotation",))
    assert any("blink_auto." + field in error for error in errors)
    assert cfg["ica"]["n_components"] is None


# -----------------------------------------------------------------------------
# Optional independent eye-aid settings must leave older configs usable
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("field,value", [("median_window_sec", 0), ("min_duration_sec", -1),
    ("threshold_std", float("nan")), ("threshold_std", True), ("merge_gap_sec", -1),
    ("padding_sec", float("inf")), ("unexpected", 1)])
def test_eye_settings_validate_without_affecting_blink_contract(tmp_path, field, value):
    cfg = _minimal_config(tmp_path)
    cfg["artifact_annotation"] = {"enabled": True, "muscle_threshold": 10., "muscle_threshold_by_run": {},
        "muscle_channel_type": "mag", "muscle_min_length_good_sec": .2, "muscle_filter_freq_hz": [110., 140.],
        "eog_mode": "disabled", "eog_channel_names": [], "blink_channel": None}
    assert pipeline_config.validate_config(cfg, required_sections=("artifact_annotation",))[0] == []
    cfg["artifact_annotation"]["eye_movement_auto"] = {"merge_gap_sec": 0., "padding_sec": 0.}
    assert pipeline_config.validate_config(cfg, required_sections=("artifact_annotation",))[0] == []
    cfg["artifact_annotation"]["eye_movement_auto"][field] = value
    errors, _ = pipeline_config.validate_config(cfg, required_sections=("artifact_annotation",))
    assert any("eye_movement_auto" in error for error in errors)


def test_tuning_config_is_merged_after_subject_config(tmp_path):
    """Accepted review choices override defaults without replacing the selector."""
    cfg = _minimal_config(tmp_path)
    cfg["maxwell_filter"] = {
        "enabled": True,
        "crosstalk_file": str(tmp_path / "ct.fif"),
        "calibration_file": str(tmp_path / "cal.dat"),
        "bad_channel_policy": "reference_run",
        "reference_acquisition_run": None,
        "manual_bad_channels": {"add": [], "remove": [], "by_run": {}},
        "destination": {
            "policy": "central_run",
            "acquisition_run": None,
        },
        "find_bad_channels": {},
        "filter": {},
    }
    cfg["artifact_annotation"] = {
        "enabled": True,
        "muscle_threshold": 10.0,
        "muscle_threshold_by_run": {},
        "muscle_channel_type": "mag",
        "muscle_min_length_good_sec": 0.2,
        "muscle_filter_freq_hz": [110.0, 140.0],
        "eog_mode": "auto",
        "eog_channel_names": [],
        "blink_channel": None,
    }
    cfg["ica"] = {
        "enabled": True,
        "method": "fastica",
        "n_components": 30,
        "random_state": 97,
        "max_iter": "auto",
        "fit_resample_sfreq_hz": 200.0,
        "fit_l_freq_hz": 1.0,
        "fit_h_freq_hz": 40.0,
        "reuse_prefit": True,
        "manual_exclude": None,
        "eog_channel_names": [],
        "eog_diagnostics": True,
    }
    raw_subject = tmp_path / "raw" / "sub01_date"
    raw_subject.mkdir(parents=True)
    default_path = tmp_path / "default.yaml"
    subject_path = tmp_path / "subject.yaml"
    tuning_path = tmp_path / "tuning.yaml"
    import yaml

    default_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    subject_path.write_text(yaml.safe_dump({"subject": {"id": "01"}}), encoding="utf-8")
    tuning_path.write_text(
        yaml.safe_dump(
            {
                "subject": {"id": "01"},
                "artifact_annotation": {
                    "muscle_threshold_by_run": {"run01": 8.0}
                },
                "ica": {"manual_exclude": [0, 3]},
            }
        ),
        encoding="utf-8",
    )

    resolved, _warnings, _default, _subject = (
        pipeline_config.load_and_resolve_config(
            default_path,
            subject_path,
            tuning_path,
            required_sections=("maxwell_filter", "artifact_annotation", "ica"),
        )
    )
    assert resolved["subject"]["id"] == "01"
    assert resolved["artifact_annotation"]["muscle_threshold_by_run"] == {
        "run01": 8.0
    }
    assert resolved["ica"]["manual_exclude"] == [0, 3]
    assert resolved["_config_sources"]["tuning"] == str(tuning_path.resolve())
