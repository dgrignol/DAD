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
     --out-json /tmp/resolved_config.json

Data flow
default YAML + optional subject YAML -> recursive merge -> validation ->
absolute raw/derivative paths -> in-memory configuration or JSON snapshot.

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


CORE_REQUIRED_SECTIONS = ("project", "layout", "data_discovery", "runtime")


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
    layout = sections["layout"]
    data_discovery = sections["data_discovery"]
    runtime = sections["runtime"]

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
            "raw_subject_folder_name",
            "derivatives_preprocessing_template",
            "derivatives_analysis_template",
            "derivatives_reports_template",
        ),
        errors,
    )
    _require_keys(
        data_discovery,
        "data_discovery",
        ("meg_fif_glob", "reported_trigger_csv_glob"),
        errors,
    )
    _require_keys(runtime, "runtime", ("overwrite_existing", "stop_on_error"), errors)
    if errors:
        return errors, warnings

    for key in ("root", "data_root", "derivatives_root", "docs_root", "params_root"):
        _expect_type(project[key], str, f"project.{key}", errors)
    for key in (
        "subject_bids_id",
        "session_bids_id",
        "raw_subject_folder_name",
        "derivatives_preprocessing_template",
        "derivatives_analysis_template",
        "derivatives_reports_template",
    ):
        _expect_type(layout[key], str, f"layout.{key}", errors)
    _expect_type(
        data_discovery["meg_fif_glob"],
        str,
        "data_discovery.meg_fif_glob",
        errors,
    )
    _expect_type(
        data_discovery["reported_trigger_csv_glob"],
        str,
        "data_discovery.reported_trigger_csv_glob",
        errors,
    )
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
                "expected_sequence_code_min",
                "expected_sequence_code_max",
                "required_trigger_codes",
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

    if "epoching" in requested:
        epoching = sections["epoching"]
        _require_keys(
            epoching,
            "epoching",
            ("tmin_sec", "tmax_sec", "baseline_sec", "event_id_map"),
            errors,
        )
        if not errors:
            tmin = epoching["tmin_sec"]
            tmax = epoching["tmax_sec"]
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


def build_resolved_paths(cfg: Mapping[str, Any]) -> dict[str, str]:
    """
    Resolve canonical subject and derivative paths.

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
    subject_data_dir = (data_root / str(layout["raw_subject_folder_name"])).resolve()

    resolved = {
        "project_root": str(project_root),
        "data_root": str(data_root),
        "subject_data_dir": str(subject_data_dir),
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
    *,
    required_sections: Iterable[str] = (),
) -> tuple[dict[str, Any], list[str], Path, Path | None]:
    """
    Load default/subject YAML, merge, validate, and attach resolved paths.

    Parameters
    ----------
    config_path
        Required default YAML.
    subject_config_path
        Optional subject override YAML.
    required_sections
        Stage-specific sections required by the caller.
    """
    default_path = Path(config_path).expanduser().resolve()
    merged = read_yaml_mapping(default_path)

    subject_path: Path | None = None
    if subject_config_path is not None:
        subject_path = Path(subject_config_path).expanduser().resolve()
        merged = deep_merge_dicts(merged, read_yaml_mapping(subject_path))

    errors, warnings = validate_config(
        merged,
        required_sections=required_sections,
    )
    if errors:
        raise ValueError(
            "Invalid configuration:\n" + "\n".join(f"- {item}" for item in errors)
        )

    merged["resolved_paths"] = build_resolved_paths(merged)
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
    parser.add_argument("--subject-config", default=None, help="Subject override YAML.")
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
