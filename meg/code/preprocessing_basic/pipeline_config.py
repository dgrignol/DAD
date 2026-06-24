#!/usr/bin/env python3
"""
pipeline_config.py

Purpose
Load, merge, validate, and resolve configuration for the DAD MEG pipeline.
This module is infrastructure, not a numbered preprocessing stage. Pipeline
scripts import its public functions, while the command-line interface is useful
for validation and for writing one resolved-configuration provenance snapshot.

Usage examples
1) Validate default and subject configuration:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/pipeline_config.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --validate-only --print-summary

2) Write the resolved configuration used for a pipeline run:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/pipeline_config.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/03_tuning_choices.yaml \
     --out-json /tmp/resolved_config.json

3) Validate cluster paths without changing scientific defaults:
   python .../pipeline_config.py --config .../pipeline_default.yaml \
     --site-config .../pipeline_cluster.yaml \
     --subject-config .../sub-02.yaml --validate-only --print-summary

Data flow
default YAML + optional site YAML + optional subject YAML + optional tuning YAML
-> recursive merge -> subject-profile selection -> validation -> external
subject-folder discovery -> absolute raw/derivative paths -> in-memory
configuration or JSON snapshot.

Configuration layers are merged in this order: default, machine/site paths,
subject selector, then reviewed tuning decisions.

Key safety contract
The configured raw-data root is read-only. All output paths must pass
`assert_output_path_outside_raw()` before a script creates directories or
writes files.
"""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable, Mapping

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment dependency.
    raise SystemExit(
        "PyYAML is required for pipeline configuration. Install it with "
        "`python -m pip install pyyaml`."
    ) from exc


SUPPORTED_SUBJECT_IDS = tuple(f"{number:02d}" for number in range(1, 8))
CORE_REQUIRED_SECTIONS = (
    "project",
    "subject",
    "subject_profiles",
    "layout",
    "data_discovery",
    "runtime",
)

PROFILE_REQUIRED_KEYS = (
    "external_folder",
    "versions",
    "condition_triggers",
    "identity_range",
    "reported_trigger_csv_pattern",
    "trigger_policy",
    "ocular_modalities",
    "completed_behavioral_blocks",
    "known_exceptions",
)

DISCOVERY_TEMPLATE_TO_GLOB = {
    "task_meg_fif_glob_template": "meg_fif_glob",
    "trajectory_mat_glob_template": "trajectory_mat_glob",
    "trial_struct_mat_glob_template": "trial_struct_mat_glob",
    "behavioral_mat_glob_template": "behavioral_mat_glob",
    "reported_trigger_csv_glob_template": "reported_trigger_csv_glob",
    "eyetracker_edf_glob_template": "eyetracker_edf_glob",
}


# -----------------------------------------------------------------------------
# YAML loading and merge
# -----------------------------------------------------------------------------


def read_yaml_mapping(path: str | Path) -> dict[str, Any]:
    """
    Read one YAML file and require a mapping at the document root.

    Parameters
    ----------
    path
        YAML file to read.

    Returns
    -------
    dict
        Parsed mapping. An empty YAML document becomes an empty mapping.
    """
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Config file not found: {resolved}")

    with resolved.open("r", encoding="utf-8") as handle:
        parsed = yaml.safe_load(handle)

    if parsed is None:
        return {}
    if not isinstance(parsed, dict):
        raise ValueError(
            f"Top-level YAML object must be a mapping in {resolved}; "
            f"got {type(parsed).__name__}."
        )
    return parsed


def deep_merge_dicts(
    base: dict[str, Any],
    override: Mapping[str, Any],
) -> dict[str, Any]:
    """
    Recursively merge mappings without mutating either input.

    Nested mappings are merged. Lists and scalar values in `override` replace
    the corresponding default value in full.
    """
    merged = deepcopy(base)
    for key, override_value in override.items():
        base_value = merged.get(key)
        if isinstance(base_value, dict) and isinstance(override_value, Mapping):
            merged[key] = deep_merge_dicts(base_value, override_value)
        else:
            merged[key] = deepcopy(override_value)
    return merged


# -----------------------------------------------------------------------------
# Subject selection and profile-derived stage configuration
# -----------------------------------------------------------------------------


def normalize_subject_id(value: Any) -> str:
    """
    Normalize a supported subject selector to a two-digit identifier.

    Accepted examples are ``1``, ``"01"``, ``"sub01"``, and ``"sub-01"``.
    Only participants 01 through 07 are valid for the current external
    multisubject dataset.
    """
    if isinstance(value, bool):
        raise ValueError("Subject ID must identify one of 01-07, not bool.")
    text = str(value).strip().lower()
    if text.startswith("sub-"):
        text = text[4:]
    elif text.startswith("sub"):
        text = text[3:]
    if not text.isdigit():
        raise ValueError(f"Subject ID must identify one of 01-07; got {value!r}.")
    normalized = f"{int(text):02d}"
    if normalized not in SUPPORTED_SUBJECT_IDS:
        raise ValueError(
            "Subject ID must identify one of "
            f"{SUPPORTED_SUBJECT_IDS[0]}-{SUPPORTED_SUBJECT_IDS[-1]}; "
            f"got {value!r}."
        )
    return normalized


def apply_subject_profile(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """
    Select one verified subject profile and derive stage-facing configuration.

    Scientific facts such as condition codes, identity range, ocular
    modalities, and version labels live only in ``subject_profiles``. This
    function copies the selected profile to ``subject_profile`` for provenance
    and derives the generic stage fields consumed by existing scripts.
    """
    resolved_cfg = deepcopy(dict(cfg))
    subject = resolved_cfg.get("subject")
    profiles = resolved_cfg.get("subject_profiles")
    if not isinstance(subject, Mapping):
        raise ValueError("`subject` must be a mapping with an `id` selector.")
    if "id" not in subject:
        raise ValueError("Missing required configuration key: `subject.id`.")
    if not isinstance(profiles, Mapping):
        raise ValueError("`subject_profiles` must be a mapping.")

    subject_id = normalize_subject_id(subject["id"])
    profile = profiles.get(subject_id)
    if not isinstance(profile, Mapping):
        raise ValueError(f"No subject profile is configured for subject {subject_id}.")

    resolved_cfg["subject"] = deepcopy(dict(subject))
    resolved_cfg["subject"]["id"] = subject_id
    resolved_cfg["subject_profile"] = deepcopy(dict(profile))

    # Derive labels and discovery globs while keeping filesystem resolution
    # separate until the configuration has passed schema validation.
    layout = resolved_cfg.setdefault("layout", {})
    layout["subject_bids_id"] = f"sub-{subject_id}"
    layout.setdefault("session_bids_id", "ses-01")

    discovery = resolved_cfg.setdefault("data_discovery", {})
    format_vars = {
        "subject_id": subject_id,
        "subject_bids_id": layout["subject_bids_id"],
        "runtime_version": profile.get("versions", {}).get("runtime", ""),
    }
    for template_key, target_key in DISCOVERY_TEMPLATE_TO_GLOB.items():
        template = discovery.get(template_key)
        if isinstance(template, str):
            discovery[target_key] = template.format(**format_vars)
    reported_pattern = profile.get("reported_trigger_csv_pattern")
    if isinstance(reported_pattern, str):
        discovery["reported_trigger_csv_glob"] = reported_pattern.format(
            **format_vars
        )

    condition_triggers = profile.get("condition_triggers", {})
    identity_range = profile.get("identity_range", [])
    trigger_policy = profile.get("trigger_policy", {})
    trigger = resolved_cfg.setdefault("trigger_parser", {})
    if isinstance(identity_range, list) and len(identity_range) == 2:
        trigger["expected_sequence_code_min"] = identity_range[0]
        trigger["expected_sequence_code_max"] = identity_range[1]
    if isinstance(profile.get("versions"), Mapping):
        trigger["runtime_version"] = profile["versions"].get("runtime")
    if isinstance(trigger_policy, Mapping):
        trigger["required_trigger_codes"] = list(
            trigger_policy.get("required_subject_codes", [])
        )
        trigger["allowed_fixed_trigger_codes"] = list(
            trigger_policy.get("allowed_fixed_codes", [])
        )
        trigger["optional_trigger_codes"] = list(
            trigger_policy.get("optional_event_codes", [])
        )
        for key in (
            "gaze_break_code",
            "replay_start_code",
            "esc_termination_code",
        ):
            if key in trigger_policy:
                trigger[key] = trigger_policy[key]

    epoching = resolved_cfg.setdefault("epoching", {})
    if isinstance(condition_triggers, Mapping):
        epoching["event_id_map"] = deepcopy(dict(condition_triggers))
    if isinstance(identity_range, list) and len(identity_range) == 2:
        epoching["identity_range"] = list(identity_range)

    ocular = profile.get("ocular_modalities", {})
    if isinstance(ocular, Mapping):
        ica = resolved_cfg.setdefault("ica", {})
        ica["eog_channel_names"] = list(ocular.get("eog_channels_in_fif", []))
        artifact = resolved_cfg.setdefault("artifact_annotation", {})
        artifact["eog_channel_names"] = list(
            ocular.get("eog_channels_in_fif", [])
        )
        eye_model = resolved_cfg.setdefault("eye_model", {})
        eye_model["saved_edf_expected"] = bool(ocular.get("saved_edf", False))

    return resolved_cfg


# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------


def _require_mapping(
    cfg: Mapping[str, Any],
    section: str,
    errors: list[str],
) -> Mapping[str, Any] | None:
    """Return a required mapping section and append a readable error if absent."""
    value = cfg.get(section)
    if not isinstance(value, Mapping):
        errors.append(f"`{section}` must be a mapping.")
        return None
    return value


def _require_keys(
    mapping: Mapping[str, Any],
    section: str,
    keys: Iterable[str],
    errors: list[str],
) -> None:
    """Append one validation error for each absent required key."""
    for key in keys:
        if key not in mapping:
            errors.append(f"Missing required configuration key: `{section}.{key}`.")


def _expect_type(
    value: Any,
    expected: type | tuple[type, ...],
    field: str,
    errors: list[str],
) -> None:
    """Append a type error while treating bool as distinct from numeric values."""
    numeric_expected = expected in (int, float) or (
        isinstance(expected, tuple) and any(item in (int, float) for item in expected)
    )
    if numeric_expected and isinstance(value, bool):
        errors.append(f"`{field}` has invalid type bool.")
        return
    if not isinstance(value, expected):
        expected_name = (
            "/".join(item.__name__ for item in expected)
            if isinstance(expected, tuple)
            else expected.__name__
        )
        errors.append(
            f"`{field}` must be {expected_name}; got {type(value).__name__}."
        )


def validate_config(
    cfg: Mapping[str, Any],
    *,
    required_sections: Iterable[str] = (),
) -> tuple[list[str], list[str]]:
    """
    Validate shared fields plus explicitly requested stage sections.

    The configuration library deliberately does not require ICA, epoching, or
    export settings for a trigger-only command. Each caller declares the
    sections it consumes through `required_sections`.

    Returns
    -------
    errors, warnings
        Fatal errors and non-fatal scientific/operational warnings.
    """
    errors: list[str] = []
    warnings: list[str] = []
    requested = tuple(dict.fromkeys((*CORE_REQUIRED_SECTIONS, *required_sections)))

    sections: dict[str, Mapping[str, Any]] = {}
    for section in requested:
        value = _require_mapping(cfg, section, errors)
        if value is not None:
            sections[section] = value
    if errors:
        return errors, warnings

    project = sections["project"]
    subject = sections["subject"]
    subject_profiles = sections["subject_profiles"]
    layout = sections["layout"]
    data_discovery = sections["data_discovery"]
    runtime = sections["runtime"]

    _require_keys(subject, "subject", ("id",), errors)
    _require_keys(
        project,
        "project",
        ("root", "data_root", "derivatives_root", "docs_root", "params_root"),
        errors,
    )
    _require_keys(
        layout,
        "layout",
        (
            "subject_bids_id",
            "session_bids_id",
            "derivatives_preprocessing_template",
            "derivatives_analysis_template",
            "derivatives_reports_template",
        ),
        errors,
    )
    _require_keys(
        data_discovery,
        "data_discovery",
        (
            "subject_folder_glob_template",
            *DISCOVERY_TEMPLATE_TO_GLOB.keys(),
            "empty_room_fif_glob_template",
            "meg_fif_glob",
            "trajectory_mat_glob",
            "trial_struct_mat_glob",
            "behavioral_mat_glob",
            "reported_trigger_csv_glob",
            "eyetracker_edf_glob",
        ),
        errors,
    )
    _require_keys(runtime, "runtime", ("overwrite_existing", "stop_on_error"), errors)
    if errors:
        return errors, warnings

    try:
        selected_subject_id = normalize_subject_id(subject["id"])
    except ValueError as exc:
        errors.append(str(exc))
        selected_subject_id = None

    missing_profiles = [
        subject_id
        for subject_id in SUPPORTED_SUBJECT_IDS
        if not isinstance(subject_profiles.get(subject_id), Mapping)
    ]
    if missing_profiles:
        errors.append(
            "`subject_profiles` must contain mappings for all supported "
            f"subjects; missing: {missing_profiles}."
        )
    for subject_id in SUPPORTED_SUBJECT_IDS:
        profile = subject_profiles.get(subject_id)
        if not isinstance(profile, Mapping):
            continue
        _require_keys(
            profile,
            f"subject_profiles.{subject_id}",
            PROFILE_REQUIRED_KEYS,
            errors,
        )
        identity_range = profile.get("identity_range")
        if (
            not isinstance(identity_range, list)
            or len(identity_range) != 2
            or any(not isinstance(value, int) for value in identity_range)
            or identity_range[0] > identity_range[1]
        ):
            errors.append(
                f"`subject_profiles.{subject_id}.identity_range` must be "
                "[minimum, maximum] integers."
            )
        condition_triggers = profile.get("condition_triggers")
        if (
            not isinstance(condition_triggers, Mapping)
            or not condition_triggers
            or any(not isinstance(value, int) for value in condition_triggers.values())
        ):
            errors.append(
                f"`subject_profiles.{subject_id}.condition_triggers` must be "
                "a non-empty mapping of condition names to integer codes."
            )
        completed_blocks = profile.get("completed_behavioral_blocks")
        if not isinstance(completed_blocks, list) or any(
            not isinstance(value, int) for value in completed_blocks
        ):
            errors.append(
                f"`subject_profiles.{subject_id}.completed_behavioral_blocks` "
                "must contain integer block identifiers."
            )
        reported_pattern = profile.get("reported_trigger_csv_pattern")
        if not isinstance(reported_pattern, str) or not reported_pattern:
            errors.append(
                f"`subject_profiles.{subject_id}.reported_trigger_csv_pattern` "
                "must be a non-empty string."
            )
        trigger_policy = profile.get("trigger_policy")
        if not isinstance(trigger_policy, Mapping):
            errors.append(
                f"`subject_profiles.{subject_id}.trigger_policy` must be a mapping."
            )
        else:
            _require_keys(
                trigger_policy,
                f"subject_profiles.{subject_id}.trigger_policy",
                (
                    "required_subject_codes",
                    "allowed_fixed_codes",
                    "optional_event_codes",
                    "gaze_break_code",
                    "replay_start_code",
                    "esc_termination_code",
                ),
                errors,
            )
            for list_key in (
                "required_subject_codes",
                "allowed_fixed_codes",
                "optional_event_codes",
            ):
                values = trigger_policy.get(list_key)
                if not isinstance(values, list) or any(
                    not isinstance(value, int) for value in values
                ):
                    errors.append(
                        f"`subject_profiles.{subject_id}.trigger_policy."
                        f"{list_key}` must contain integer codes."
                    )
            for code_key in (
                "gaze_break_code",
                "replay_start_code",
                "esc_termination_code",
            ):
                if not isinstance(trigger_policy.get(code_key), int):
                    errors.append(
                        f"`subject_profiles.{subject_id}.trigger_policy."
                        f"{code_key}` must be an integer."
                    )
    if (
        selected_subject_id is not None
        and selected_subject_id not in subject_profiles
    ):
        errors.append(
            f"No subject profile is configured for subject {selected_subject_id}."
        )

    for key in ("root", "data_root", "derivatives_root", "docs_root", "params_root"):
        _expect_type(project[key], str, f"project.{key}", errors)
    for key in (
        "subject_bids_id",
        "session_bids_id",
        "derivatives_preprocessing_template",
        "derivatives_analysis_template",
        "derivatives_reports_template",
    ):
        _expect_type(layout[key], str, f"layout.{key}", errors)
    for key in (
        "subject_folder_glob_template",
        *DISCOVERY_TEMPLATE_TO_GLOB.keys(),
        "empty_room_fif_glob_template",
        "meg_fif_glob",
        "trajectory_mat_glob",
        "trial_struct_mat_glob",
        "behavioral_mat_glob",
        "reported_trigger_csv_glob",
        "eyetracker_edf_glob",
    ):
        _expect_type(data_discovery[key], str, f"data_discovery.{key}", errors)
    _expect_type(
        runtime["overwrite_existing"],
        bool,
        "runtime.overwrite_existing",
        errors,
    )
    _expect_type(runtime["stop_on_error"], bool, "runtime.stop_on_error", errors)

    if "trigger_parser" in requested:
        trigger = sections["trigger_parser"]
        _require_keys(
            trigger,
            "trigger_parser",
            (
                "stim_channel",
                "runtime_version",
                "expected_sequence_code_min",
                "expected_sequence_code_max",
                "required_trigger_codes",
                "allowed_fixed_trigger_codes",
                "optional_trigger_codes",
                "gaze_break_code",
                "replay_start_code",
                "esc_termination_code",
                "enable_deep_confrontation_report",
                "deep_report_output_subdir",
                "deep_report_anomaly_window_ms",
                "deep_report_transition_window_ms",
            ),
            errors,
        )
        if not errors:
            _expect_type(
                trigger["stim_channel"],
                str,
                "trigger_parser.stim_channel",
                errors,
            )
            _expect_type(
                trigger["runtime_version"],
                str,
                "trigger_parser.runtime_version",
                errors,
            )
            seq_min = trigger["expected_sequence_code_min"]
            seq_max = trigger["expected_sequence_code_max"]
            _expect_type(seq_min, int, "trigger_parser.expected_sequence_code_min", errors)
            _expect_type(seq_max, int, "trigger_parser.expected_sequence_code_max", errors)
            if isinstance(seq_min, int) and isinstance(seq_max, int) and seq_min > seq_max:
                errors.append(
                    "`trigger_parser.expected_sequence_code_min` must be <= "
                    "`trigger_parser.expected_sequence_code_max`."
                )
            required_codes = trigger["required_trigger_codes"]
            if not isinstance(required_codes, list) or not required_codes:
                errors.append(
                    "`trigger_parser.required_trigger_codes` must be a non-empty list."
                )
            elif any(not isinstance(code, int) for code in required_codes):
                errors.append(
                    "`trigger_parser.required_trigger_codes` must contain only integers."
                )
            for list_key in (
                "allowed_fixed_trigger_codes",
                "optional_trigger_codes",
            ):
                values = trigger[list_key]
                if not isinstance(values, list) or any(
                    not isinstance(code, int) for code in values
                ):
                    errors.append(
                        f"`trigger_parser.{list_key}` must contain only integers."
                    )
            for code_key in (
                "gaze_break_code",
                "replay_start_code",
                "esc_termination_code",
            ):
                _expect_type(
                    trigger[code_key],
                    int,
                    f"trigger_parser.{code_key}",
                    errors,
                )
            _expect_type(
                trigger["enable_deep_confrontation_report"],
                bool,
                "trigger_parser.enable_deep_confrontation_report",
                errors,
            )
            _expect_type(
                trigger["deep_report_output_subdir"],
                str,
                "trigger_parser.deep_report_output_subdir",
                errors,
            )
            for window_key in (
                "deep_report_anomaly_window_ms",
                "deep_report_transition_window_ms",
            ):
                _expect_type(
                    trigger[window_key],
                    (int, float),
                    f"trigger_parser.{window_key}",
                    errors,
                )
                if (
                    isinstance(trigger[window_key], (int, float))
                    and not isinstance(trigger[window_key], bool)
                    and trigger[window_key] <= 0
                ):
                    errors.append(
                        f"`trigger_parser.{window_key}` must be greater than zero."
                    )

    if "epoching" in requested:
        epoching = sections["epoching"]
        _require_keys(
            epoching,
            "epoching",
            (
                "enabled",
                "tmin_sec",
                "tmax_sec",
                "baseline_sec",
                "detrend",
                "reject_by_annotation",
                "event_id_map",
                "identity_range",
                "amplitude_reject",
                "trial_exclusion",
            ),
            errors,
        )
        if not errors:
            tmin = epoching["tmin_sec"]
            tmax = epoching["tmax_sec"]
            _expect_type(epoching["enabled"], bool, "epoching.enabled", errors)
            _expect_type(tmin, (int, float), "epoching.tmin_sec", errors)
            _expect_type(tmax, (int, float), "epoching.tmax_sec", errors)
            if (
                isinstance(tmin, (int, float))
                and not isinstance(tmin, bool)
                and isinstance(tmax, (int, float))
                and not isinstance(tmax, bool)
                and tmin >= tmax
            ):
                errors.append("`epoching.tmin_sec` must be smaller than `tmax_sec`.")
            baseline = epoching["baseline_sec"]
            if baseline is not None and (
                not isinstance(baseline, list)
                or len(baseline) != 2
                or any(
                    value is not None
                    and (
                        not isinstance(value, (int, float))
                        or isinstance(value, bool)
                    )
                    for value in baseline
                )
            ):
                errors.append(
                    "`epoching.baseline_sec` must be null or a two-item "
                    "numeric/null list."
                )
            _expect_type(epoching["detrend"], int, "epoching.detrend", errors)
            _expect_type(
                epoching["reject_by_annotation"],
                bool,
                "epoching.reject_by_annotation",
                errors,
            )
            event_map = epoching["event_id_map"]
            if not isinstance(event_map, Mapping) or any(
                not isinstance(name, str) or not isinstance(code, int)
                for name, code in event_map.items()
            ):
                errors.append(
                    "`epoching.event_id_map` must map condition names to "
                    "integer trigger codes."
                )
            identity_range = epoching["identity_range"]
            if (
                not isinstance(identity_range, list)
                or len(identity_range) != 2
                or any(not isinstance(value, int) for value in identity_range)
                or identity_range[0] > identity_range[1]
            ):
                errors.append(
                    "`epoching.identity_range` must be [minimum, maximum] integers."
                )
            amplitude_reject = epoching["amplitude_reject"]
            if not isinstance(amplitude_reject, Mapping) or any(
                not isinstance(name, str)
                or not isinstance(value, (int, float))
                or isinstance(value, bool)
                or value <= 0
                for name, value in amplitude_reject.items()
            ):
                errors.append(
                    "`epoching.amplitude_reject` must map channel types to "
                    "positive numeric thresholds."
                )
            trial_exclusion = epoching["trial_exclusion"]
            if not isinstance(trial_exclusion, Mapping):
                errors.append("`epoching.trial_exclusion` must be a mapping.")
            else:
                required_trial_keys = (
                    "gaze_break_codes",
                    "catch_type_1_codes",
                    "catch_type_2_codes",
                    "esc_termination_codes",
                    "exclude_gaze_break",
                    "exclude_catch_type_1",
                    "exclude_catch_type_2",
                    "exclude_esc_termination",
                    "exclude_replay_trials",
                    "exclude_unresolved_transitions",
                )
                _require_keys(
                    trial_exclusion,
                    "epoching.trial_exclusion",
                    required_trial_keys,
                    errors,
                )
                for key in (
                    "gaze_break_codes",
                    "catch_type_1_codes",
                    "catch_type_2_codes",
                    "esc_termination_codes",
                ):
                    values = trial_exclusion.get(key)
                    if not isinstance(values, list) or any(
                        not isinstance(value, int) for value in values
                    ):
                        errors.append(
                            f"`epoching.trial_exclusion.{key}` must contain "
                            "integer trigger codes."
                        )
                for key in (
                    "exclude_gaze_break",
                    "exclude_catch_type_1",
                    "exclude_catch_type_2",
                    "exclude_esc_termination",
                    "exclude_replay_trials",
                    "exclude_unresolved_transitions",
                ):
                    if key in trial_exclusion:
                        _expect_type(
                            trial_exclusion[key],
                            bool,
                            f"epoching.trial_exclusion.{key}",
                            errors,
                        )

    if "export" in requested:
        export = sections["export"]
        _require_keys(
            export,
            "export",
            (
                "enabled",
                "output_format",
                "tmin_sec",
                "tmax_sec",
                "sampling_rate_hz",
                "canonical_meg_channels_file",
                "missing_identity_policy",
                "bad_or_missing_channel_policy",
                "include_epochs_metadata",
            ),
            errors,
        )
        if not errors:
            _expect_type(export["enabled"], bool, "export.enabled", errors)
            if export["output_format"] != "mat":
                errors.append("`export.output_format` must be 'mat'.")
            for key in ("tmin_sec", "tmax_sec", "sampling_rate_hz"):
                _expect_type(export[key], (int, float), f"export.{key}", errors)
            if (
                isinstance(export["tmin_sec"], (int, float))
                and isinstance(export["tmax_sec"], (int, float))
                and export["tmin_sec"] >= export["tmax_sec"]
            ):
                errors.append("`export.tmin_sec` must be smaller than `export.tmax_sec`.")
            if (
                isinstance(export["sampling_rate_hz"], (int, float))
                and export["sampling_rate_hz"] <= 0
            ):
                errors.append("`export.sampling_rate_hz` must be positive.")
            _expect_type(
                export["canonical_meg_channels_file"],
                str,
                "export.canonical_meg_channels_file",
                errors,
            )
            if export["missing_identity_policy"] not in ("error", "nan"):
                errors.append(
                    "`export.missing_identity_policy` must be 'error' or 'nan'."
                )
            if export["bad_or_missing_channel_policy"] != "nan":
                errors.append(
                    "`export.bad_or_missing_channel_policy` must be 'nan'."
                )
            _expect_type(
                export["include_epochs_metadata"],
                bool,
                "export.include_epochs_metadata",
                errors,
            )

    if "matlab_drsa" in requested:
        matlab_drsa = sections["matlab_drsa"]
        _require_keys(
            matlab_drsa,
            "matlab_drsa",
            (
                "enabled",
                "matlab_bin",
                "use_arch_x86_64",
                "analyses",
                "drsa_type",
                "conditions_to_run",
                "sensor_groups_to_run",
                "model_source",
                "pcr_strategy",
                "smooth_neural",
                "smooth_window_samples",
            ),
            errors,
        )
        if not errors:
            _expect_type(matlab_drsa["enabled"], bool, "matlab_drsa.enabled", errors)
            _expect_type(matlab_drsa["matlab_bin"], str, "matlab_drsa.matlab_bin", errors)
            _expect_type(
                matlab_drsa["use_arch_x86_64"],
                bool,
                "matlab_drsa.use_arch_x86_64",
                errors,
            )
            analyses = matlab_drsa["analyses"]
            if not isinstance(analyses, list) or not analyses or any(
                value not in ("real", "alt", "simulations") for value in analyses
            ):
                errors.append(
                    "`matlab_drsa.analyses` must contain one or more of "
                    "'real', 'alt', or 'simulations'."
                )
            if matlab_drsa["drsa_type"] not in ("corr", "PCR"):
                errors.append("`matlab_drsa.drsa_type` must be 'corr' or 'PCR'.")
            valid_conditions = {
                "always_visible",
                "occluded_nondeviant",
                "occluded_deviant",
            }
            conditions = matlab_drsa["conditions_to_run"]
            if not isinstance(conditions, list) or not conditions or any(
                value not in valid_conditions for value in conditions
            ):
                errors.append(
                    "`matlab_drsa.conditions_to_run` must contain valid condition names."
                )
            groups = matlab_drsa["sensor_groups_to_run"]
            if not isinstance(groups, list) or not groups or any(
                not isinstance(value, str) or not value for value in groups
            ):
                errors.append(
                    "`matlab_drsa.sensor_groups_to_run` must contain non-empty strings."
                )
            valid_model_sources = ("observed", "predicted", "prediction_error")
            if matlab_drsa["model_source"] not in valid_model_sources:
                errors.append(
                    "`matlab_drsa.model_source` must be 'observed', 'predicted', "
                    "or 'prediction_error'."
                )
            model_sources = matlab_drsa.get("model_sources")
            if model_sources is not None:
                if isinstance(model_sources, str):
                    model_sources_to_check = [model_sources]
                else:
                    model_sources_to_check = model_sources
                if (
                    not isinstance(model_sources_to_check, list)
                    or not model_sources_to_check
                    or any(value not in valid_model_sources for value in model_sources_to_check)
                ):
                    errors.append(
                        "`matlab_drsa.model_sources` must contain one or more of "
                        "'observed', 'predicted', or 'prediction_error'."
                    )
            _expect_type(
                matlab_drsa["pcr_strategy"],
                str,
                "matlab_drsa.pcr_strategy",
                errors,
            )
            _expect_type(
                matlab_drsa["smooth_neural"],
                bool,
                "matlab_drsa.smooth_neural",
                errors,
            )
            _expect_type(
                matlab_drsa["smooth_window_samples"],
                int,
                "matlab_drsa.smooth_window_samples",
                errors,
            )
            if (
                isinstance(matlab_drsa["smooth_window_samples"], int)
                and matlab_drsa["smooth_window_samples"] <= 0
            ):
                errors.append(
                    "`matlab_drsa.smooth_window_samples` must be positive."
                )

    if "eye_model" in requested:
        eye_model = sections["eye_model"]
        _require_keys(eye_model, "eye_model", ("enabled", "skip_if_missing"), errors)
        if not errors:
            _expect_type(eye_model["enabled"], bool, "eye_model.enabled", errors)
            _expect_type(
                eye_model["skip_if_missing"],
                bool,
                "eye_model.skip_if_missing",
                errors,
            )

    if "maxwell_filter" in requested:
        maxwell = sections["maxwell_filter"]
        _require_keys(
            maxwell,
            "maxwell_filter",
            (
                "enabled",
                "crosstalk_file",
                "calibration_file",
                "bad_channel_policy",
                "reference_acquisition_run",
                "manual_bad_channels",
                "destination",
                "find_bad_channels",
                "filter",
            ),
            errors,
        )
        if not errors:
            _expect_type(maxwell["enabled"], bool, "maxwell_filter.enabled", errors)
            for key in ("crosstalk_file", "calibration_file"):
                _expect_type(maxwell[key], str, f"maxwell_filter.{key}", errors)
            policy = maxwell["bad_channel_policy"]
            if policy not in ("reference_run", "union", "per_run"):
                errors.append(
                    "`maxwell_filter.bad_channel_policy` must be one of "
                    "'reference_run', 'union', or 'per_run'."
                )
            reference_run = maxwell["reference_acquisition_run"]
            if reference_run is not None:
                _expect_type(
                    reference_run,
                    int,
                    "maxwell_filter.reference_acquisition_run",
                    errors,
                )
            manual_bads = maxwell["manual_bad_channels"]
            if not isinstance(manual_bads, Mapping):
                errors.append("`maxwell_filter.manual_bad_channels` must be a mapping.")
            else:
                _require_keys(
                    manual_bads,
                    "maxwell_filter.manual_bad_channels",
                    ("add", "remove", "by_run"),
                    errors,
                )
                for key in ("add", "remove"):
                    values = manual_bads.get(key)
                    if not isinstance(values, list) or any(
                        not isinstance(item, str) for item in values
                    ):
                        errors.append(
                            f"`maxwell_filter.manual_bad_channels.{key}` must "
                            "contain channel-name strings."
                        )
                if not isinstance(manual_bads.get("by_run"), Mapping):
                    errors.append(
                        "`maxwell_filter.manual_bad_channels.by_run` must be a mapping."
                    )
            destination = maxwell["destination"]
            if not isinstance(destination, Mapping):
                errors.append("`maxwell_filter.destination` must be a mapping.")
            else:
                _require_keys(
                    destination,
                    "maxwell_filter.destination",
                    (
                        "policy",
                        "acquisition_run",
                    ),
                    errors,
                )
                destination_policy = destination.get("policy")
                if destination_policy != "central_run":
                    errors.append(
                        "`maxwell_filter.destination.policy` must be "
                        "'central_run'."
                    )
                destination_run = destination.get("acquisition_run")
                if destination_run is not None:
                    errors.append(
                        "`maxwell_filter.destination.acquisition_run` must be "
                        "null because the central run is selected from data."
                    )
            find_cfg = maxwell["find_bad_channels"]
            if not isinstance(find_cfg, Mapping):
                errors.append("`maxwell_filter.find_bad_channels` must be a mapping.")
            filter_cfg = maxwell["filter"]
            if not isinstance(filter_cfg, Mapping):
                errors.append("`maxwell_filter.filter` must be a mapping.")

    if "artifact_annotation" in requested:
        artifact = sections["artifact_annotation"]
        _require_keys(
            artifact,
            "artifact_annotation",
            (
                "enabled",
                "muscle_threshold",
                "muscle_threshold_by_run",
                "muscle_channel_type",
                "muscle_min_length_good_sec",
                "muscle_filter_freq_hz",
                "eog_mode",
                "eog_channel_names",
                "blink_channel",
            ),
            errors,
        )
        if not errors:
            _expect_type(
                artifact["enabled"],
                bool,
                "artifact_annotation.enabled",
                errors,
            )
            _expect_type(
                artifact["muscle_threshold"],
                (int, float),
                "artifact_annotation.muscle_threshold",
                errors,
            )
            if not isinstance(artifact["muscle_threshold_by_run"], Mapping):
                errors.append(
                    "`artifact_annotation.muscle_threshold_by_run` must be a mapping."
                )
            _expect_type(
                artifact["muscle_channel_type"],
                str,
                "artifact_annotation.muscle_channel_type",
                errors,
            )
            _expect_type(
                artifact["muscle_min_length_good_sec"],
                (int, float),
                "artifact_annotation.muscle_min_length_good_sec",
                errors,
            )
            frequencies = artifact["muscle_filter_freq_hz"]
            if (
                not isinstance(frequencies, list)
                or len(frequencies) != 2
                or any(
                    not isinstance(value, (int, float)) or isinstance(value, bool)
                    for value in frequencies
                )
                or frequencies[0] >= frequencies[1]
            ):
                errors.append(
                    "`artifact_annotation.muscle_filter_freq_hz` must be an "
                    "increasing two-number list."
                )
            if artifact["eog_mode"] not in (
                "auto",
                "disabled",
                "diagnostics_only",
                "annotate",
            ):
                errors.append(
                    "`artifact_annotation.eog_mode` must be one of 'auto', "
                    "'disabled', 'diagnostics_only', or 'annotate'."
                )
            eog_names = artifact["eog_channel_names"]
            if not isinstance(eog_names, list) or any(
                not isinstance(item, str) for item in eog_names
            ):
                errors.append(
                    "`artifact_annotation.eog_channel_names` must contain "
                    "channel-name strings."
                )
            blink_channel = artifact["blink_channel"]
            if blink_channel is not None and not isinstance(blink_channel, str):
                errors.append(
                    "`artifact_annotation.blink_channel` must be null or a string."
                )

    if "ica" in requested:
        ica = sections["ica"]
        _require_keys(
            ica,
            "ica",
            (
                "enabled",
                "method",
                "n_components",
                "random_state",
                "max_iter",
                "fit_resample_sfreq_hz",
                "fit_l_freq_hz",
                "fit_h_freq_hz",
                "reuse_prefit",
                "manual_exclude",
                "eog_channel_names",
                "eog_diagnostics",
            ),
            errors,
        )
        if not errors:
            _expect_type(ica["enabled"], bool, "ica.enabled", errors)
            _expect_type(ica["method"], str, "ica.method", errors)
            _expect_type(
                ica["n_components"],
                (int, float),
                "ica.n_components",
                errors,
            )
            _expect_type(ica["random_state"], int, "ica.random_state", errors)
            if not isinstance(ica["max_iter"], (int, str)):
                errors.append("`ica.max_iter` must be an integer or string.")
            for key in (
                "fit_resample_sfreq_hz",
                "fit_l_freq_hz",
                "fit_h_freq_hz",
            ):
                _expect_type(ica[key], (int, float), f"ica.{key}", errors)
            _expect_type(ica["reuse_prefit"], bool, "ica.reuse_prefit", errors)
            manual_exclude = ica["manual_exclude"]
            if manual_exclude is not None and (
                not isinstance(manual_exclude, list)
                or any(not isinstance(item, int) for item in manual_exclude)
            ):
                errors.append(
                    "`ica.manual_exclude` must be null or a list of component indices."
                )
            eog_names = ica["eog_channel_names"]
            if not isinstance(eog_names, list) or any(
                not isinstance(item, str) for item in eog_names
            ):
                errors.append(
                    "`ica.eog_channel_names` must contain channel-name strings."
                )
            _expect_type(
                ica["eog_diagnostics"],
                bool,
                "ica.eog_diagnostics",
                errors,
            )

    return errors, warnings


# -----------------------------------------------------------------------------
# Path resolution and raw-data write protection
# -----------------------------------------------------------------------------


def _format_derivative_path(
    derivatives_root: Path,
    template: str,
    format_vars: Mapping[str, str],
) -> Path:
    """
    Resolve one derivative template beneath the authoritative derivative root.

    Templates must be relative. This prevents a subject override from silently
    redirecting outputs to an arbitrary absolute location.
    """
    formatted = Path(template.format(**format_vars))
    if formatted.is_absolute():
        raise ValueError(f"Derivative template must be relative: {template}")
    return (derivatives_root / formatted).resolve()


def is_relative_to(path: str | Path, parent: str | Path) -> bool:
    """Return whether `path` is equal to or located below `parent`."""
    resolved_path = Path(path).expanduser().resolve()
    resolved_parent = Path(parent).expanduser().resolve()
    try:
        resolved_path.relative_to(resolved_parent)
    except ValueError:
        return False
    return True


def assert_output_path_outside_raw(
    output_path: str | Path,
    raw_root: str | Path,
) -> Path:
    """
    Reject any output located in or below the configured raw-data root.

    Returns the resolved output path so callers can use the validated value.
    This guard must run before directory creation or file writing.
    """
    resolved_output = Path(output_path).expanduser().resolve()
    resolved_raw = Path(raw_root).expanduser().resolve()
    if is_relative_to(resolved_output, resolved_raw):
        raise ValueError(
            "Refusing to write inside the read-only raw-data tree: "
            f"output={resolved_output}, raw_root={resolved_raw}"
        )
    return resolved_output


def discover_external_subject_folder(
    data_root: str | Path,
    subject_id: Any,
    *,
    subject_folder_glob_template: str = "sub{subject_id}_*",
) -> Path:
    """
    Resolve exactly one external raw folder for a supported subject.

    Discovery intentionally uses only the stable ``subNN_`` prefix. Dates in
    filenames and migration notes are provenance facts, not safe selectors;
    this is important for the known sub07 date-label discrepancy.
    """
    normalized_id = normalize_subject_id(subject_id)
    root = Path(data_root).expanduser().resolve()
    pattern = subject_folder_glob_template.format(subject_id=normalized_id)
    matches = sorted(path.resolve() for path in root.glob(pattern) if path.is_dir())
    if not matches:
        raise FileNotFoundError(
            f"No external subject folder matched {pattern!r} under {root}."
        )
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous external subject folders for subject {normalized_id}: "
            + ", ".join(str(path) for path in matches)
        )
    return matches[0]


def resolve_subject_discovery_paths(
    cfg: dict[str, Any],
    resolved_paths: Mapping[str, str],
) -> None:
    """
    Attach exact raw search roots and the folder-date empty-room glob in place.

    Existing stage scripts consume ``reported_trigger_search_roots`` and
    ``eye_search_glob_roots``. Deriving them here preserves that API while
    ensuring all searches remain inside the uniquely discovered subject tree.
    """
    discovery = cfg["data_discovery"]
    subject_dir = Path(resolved_paths["subject_data_dir"])
    folder_date = str(resolved_paths["subject_folder_date"])

    # The configured globs include their output_files/eyetracker prefixes, so
    # legacy stage helpers must search from the subject root.
    discovery["reported_trigger_search_roots"] = [str(subject_dir)]
    discovery["eye_search_glob_roots"] = [str(subject_dir)]
    discovery["empty_room_search_roots"] = [resolved_paths["data_root"]]
    discovery["empty_room_fif_glob"] = str(
        discovery["empty_room_fif_glob_template"]
    ).format(subject_folder_date=folder_date)
    discovery["subject_data_dir"] = str(subject_dir)


def build_resolved_paths(cfg: Mapping[str, Any]) -> dict[str, str]:
    """
    Discover the external subject folder and resolve canonical output paths.

    `project.derivatives_root` is authoritative. Layout templates are relative
    to that root and must not include an additional leading `derivatives/`.
    """
    project = cfg["project"]
    layout = cfg["layout"]
    format_vars = {
        "subject_bids_id": str(layout["subject_bids_id"]),
        "session_bids_id": str(layout["session_bids_id"]),
    }

    project_root = Path(project["root"]).expanduser().resolve()
    data_root = Path(project["data_root"]).expanduser().resolve()
    derivatives_root = Path(project["derivatives_root"]).expanduser().resolve()
    subject_data_dir = discover_external_subject_folder(
        data_root,
        cfg["subject"]["id"],
        subject_folder_glob_template=str(
            cfg["data_discovery"]["subject_folder_glob_template"]
        ),
    )
    if "_" not in subject_data_dir.name:
        raise ValueError(
            "External subject folder must preserve its date suffix: "
            f"{subject_data_dir}"
        )
    subject_folder_date = subject_data_dir.name.split("_", 1)[1]

    resolved = {
        "project_root": str(project_root),
        "data_root": str(data_root),
        "subject_data_dir": str(subject_data_dir),
        "subject_folder_name": subject_data_dir.name,
        "subject_folder_date": subject_folder_date,
        "subject_meg_dir": str((subject_data_dir / "meg").resolve()),
        "subject_input_files_dir": str((subject_data_dir / "input_files").resolve()),
        "subject_output_files_dir": str((subject_data_dir / "output_files").resolve()),
        "subject_eyetracker_dir": str((subject_data_dir / "eyetracker").resolve()),
        "derivatives_root": str(derivatives_root),
        "derivatives_preprocessing_dir": str(
            _format_derivative_path(
                derivatives_root,
                str(layout["derivatives_preprocessing_template"]),
                format_vars,
            )
        ),
        "derivatives_analysis_dir": str(
            _format_derivative_path(
                derivatives_root,
                str(layout["derivatives_analysis_template"]),
                format_vars,
            )
        ),
        "derivatives_reports_dir": str(
            _format_derivative_path(
                derivatives_root,
                str(layout["derivatives_reports_template"]),
                format_vars,
            )
        ),
    }

    for key in (
        "derivatives_root",
        "derivatives_preprocessing_dir",
        "derivatives_analysis_dir",
        "derivatives_reports_dir",
    ):
        assert_output_path_outside_raw(resolved[key], data_root)
    return resolved


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------


def load_and_resolve_config(
    config_path: str | Path,
    subject_config_path: str | Path | None = None,
    tuning_config_path: str | Path | None = None,
    *,
    site_config_path: str | Path | None = None,
    required_sections: Iterable[str] = (),
) -> tuple[dict[str, Any], list[str], Path, Path | None]:
    """
    Load layered YAML, validate it, and attach resolved paths.

    Merge order is default -> site -> subject -> tuning. Site overrides are
    intended for machine-specific paths such as the CIMeC cluster mount and
    must not contain participant decisions. Tuning remains the final layer so
    accepted scientific choices take precedence in the apply pass.

    Parameters
    ----------
    config_path
        Required default YAML.
    subject_config_path
        Optional subject override YAML.
    tuning_config_path
        Optional review/apply tuning YAML. It is merged after the subject
        selector so accepted thresholds and ICA exclusions take precedence.
    site_config_path
        Optional machine/site override YAML, merged immediately after the
        default config and before subject/tuning settings.
    required_sections
        Stage-specific sections required by the caller.
    """
    default_path = Path(config_path).expanduser().resolve()
    merged = read_yaml_mapping(default_path)

    site_path: Path | None = None
    if site_config_path is not None:
        site_path = Path(site_config_path).expanduser().resolve()
        merged = deep_merge_dicts(merged, read_yaml_mapping(site_path))

    subject_path: Path | None = None
    if subject_config_path is not None:
        subject_path = Path(subject_config_path).expanduser().resolve()
        merged = deep_merge_dicts(merged, read_yaml_mapping(subject_path))

    tuning_path: Path | None = None
    if tuning_config_path is not None:
        tuning_path = Path(tuning_config_path).expanduser().resolve()
        merged = deep_merge_dicts(merged, read_yaml_mapping(tuning_path))

    merged = apply_subject_profile(merged)
    errors, warnings = validate_config(
        merged,
        required_sections=required_sections,
    )
    if errors:
        raise ValueError(
            "Invalid configuration:\n" + "\n".join(f"- {item}" for item in errors)
        )

    merged["resolved_paths"] = build_resolved_paths(merged)
    resolve_subject_discovery_paths(merged, merged["resolved_paths"])

    expected_folder = merged["subject_profile"].get("external_folder", {}).get(
        "expected_name"
    )
    actual_folder = merged["resolved_paths"]["subject_folder_name"]
    if expected_folder and expected_folder != actual_folder:
        warnings.append(
            "Discovered subject folder differs from the profile's documented "
            f"name: expected={expected_folder}, discovered={actual_folder}. "
            "The unique subNN_* match is authoritative for raw discovery."
        )
    merged["_config_sources"] = {
        "default": str(default_path),
        "site": str(site_path) if site_path else None,
        "subject": str(subject_path) if subject_path else None,
        "tuning": str(tuning_path) if tuning_path else None,
    }
    return merged, warnings, default_path, subject_path


def write_resolved_config(path: str | Path, cfg: Mapping[str, Any]) -> Path:
    """Write one stable JSON snapshot after enforcing raw-data write safety."""
    raw_root = cfg["resolved_paths"]["data_root"]
    resolved_path = assert_output_path_outside_raw(path, raw_root)
    resolved_path.parent.mkdir(parents=True, exist_ok=True)
    with resolved_path.open("w", encoding="utf-8") as handle:
        json.dump(cfg, handle, indent=2, sort_keys=True)
    return resolved_path


# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the configuration validation/snapshot CLI."""
    parser = argparse.ArgumentParser(
        description="Validate and resolve DAD MEG pipeline configuration."
    )
    parser.add_argument("--config", required=True, help="Default YAML config.")
    parser.add_argument(
        "--site-config",
        default=None,
        help="Optional machine/site path override YAML.",
    )
    parser.add_argument("--subject-config", default=None, help="Subject override YAML.")
    parser.add_argument(
        "--tuning-config",
        default=None,
        help="Optional review/apply parameter override YAML.",
    )
    parser.add_argument(
        "--require-section",
        action="append",
        default=[],
        help="Additional required section; may be supplied more than once.",
    )
    parser.add_argument("--out-json", default=None, help="Resolved config JSON output.")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate without writing --out-json.",
    )
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--print-resolved", action="store_true")
    return parser


def main() -> int:
    """Validate configuration and optionally write/print its resolved form."""
    args = build_arg_parser().parse_args()
    try:
        cfg, warnings, _default_path, _subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            args.tuning_config,
            site_config_path=args.site_config,
            required_sections=args.require_section,
        )
        if args.out_json and not args.validate_only:
            write_resolved_config(args.out_json, cfg)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2

    if args.print_summary:
        resolved = cfg["resolved_paths"]
        print("Configuration: valid")
        print(
            f"Subject/session: {cfg['layout']['subject_bids_id']} / "
            f"{cfg['layout']['session_bids_id']}"
        )
        print(f"Raw subject directory: {resolved['subject_data_dir']}")
        print(f"Derivative root: {resolved['derivatives_root']}")
        for warning in warnings:
            print(f"Warning: {warning}")
    if args.print_resolved:
        print(json.dumps(cfg, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
