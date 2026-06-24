#!/usr/bin/env python3
"""
dataset_inventory.py

Purpose
Discover and validate external multisubject inputs before expensive MEG
processing. This module is reusable infrastructure, not a numbered scientific
stage. It also builds the durable source manifest written by run_pipeline.py.

Usage examples
1) Run read-only preflight for one subject:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/dataset_inventory.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-07.yaml \
     --print-summary

2) Write optional QC reports and a durable provenance manifest:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/dataset_inventory.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-07.yaml \
     --out-json /tmp/dataset_inventory.json \
     --out-md /tmp/dataset_inventory.md \
     --source-manifest /tmp/source_manifest.json

Data flow
resolved configuration -> explicit source-category discovery -> lazy FIF-header
inspection -> optional inventory QC and durable source manifest.

Scientific boundary
MEG files are acquisition runs; behavioral MAT, trigger CSV, and EDF files are
organized by experiment block. Discovery keeps those lists separate and never
zips them positionally. The manifest reserves an empty association table for a
later trigger-derived mapping.

Safety
The external raw tree is opened read-only. Every output passes
``assert_output_path_outside_raw()`` before directories are created or files
are written.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import mne

from pipeline_config import (
    assert_output_path_outside_raw,
    load_and_resolve_config,
)


INVENTORY_SCHEMA_VERSION = "4.0"
SOURCE_MANIFEST_SCHEMA_VERSION = "2.0"


# -----------------------------------------------------------------------------
# Explicit external-source discovery
# -----------------------------------------------------------------------------


def _iso_utc_now() -> str:
    """Return a UTC ISO8601 timestamp with explicit Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _glob_files(root: Path, pattern: str) -> list[Path]:
    """Return sorted unique files matching one configured relative glob."""
    return sorted({path.resolve() for path in root.glob(pattern) if path.is_file()})


def discover_meg_fif_files(
    subject_data_dir: str | Path,
    meg_fif_glob: str,
) -> list[Path]:
    """
    Return task acquisition FIFs matching the subject-specific task pattern.

    The configured pattern is rooted at the subject directory and starts with
    ``meg/subNN_run``. Root-level empty-room files therefore cannot enter the
    task-run list.
    """
    return _glob_files(Path(subject_data_dir).expanduser().resolve(), meg_fif_glob)


def discover_subject_sources(cfg: Mapping[str, Any]) -> dict[str, list[Path]]:
    """
    Reproduce all source discovery directly from resolved configuration.

    Required scientific source categories stay separate. In particular, no
    positional relationship is inferred among task FIFs, behavioral blocks,
    debug-trigger chunks, or EDF files.
    """
    resolved = cfg["resolved_paths"]
    discovery = cfg["data_discovery"]
    subject_root = Path(resolved["subject_data_dir"]).resolve()
    data_root = Path(resolved["data_root"]).resolve()

    return {
        "task_meg_fif": discover_meg_fif_files(
            subject_root,
            str(discovery["meg_fif_glob"]),
        ),
        "empty_room_fif": _glob_files(
            data_root,
            str(discovery["empty_room_fif_glob"]),
        ),
        "trajectory_mat": _glob_files(
            subject_root,
            str(discovery["trajectory_mat_glob"]),
        ),
        "trial_struct_mat": _glob_files(
            subject_root,
            str(discovery["trial_struct_mat_glob"]),
        ),
        "behavioral_mat": _glob_files(
            subject_root,
            str(discovery["behavioral_mat_glob"]),
        ),
        "reported_trigger_csv": _glob_files(
            subject_root,
            str(discovery["reported_trigger_csv_glob"]),
        ),
        "eyetracker_edf": _glob_files(
            subject_root,
            str(discovery["eyetracker_edf_glob"]),
        ),
    }


def summarize_fif_file(fif_path: str | Path, stim_channel: str) -> dict[str, Any]:
    """
    Read only a FIF header and summarize acquisition/channel metadata.

    The raw object is not preloaded. MISC channels are reported explicitly so
    the left-eye analog traces remain visible for future EDF-to-MEG alignment.
    """
    path = Path(fif_path).expanduser().resolve()
    raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
    channel_types = raw.get_channel_types()
    return {
        "path": str(path),
        "name": path.name,
        "size_bytes": path.stat().st_size,
        "sfreq_hz": float(raw.info["sfreq"]),
        "n_channels": int(raw.info["nchan"]),
        "n_samples": int(raw.n_times),
        "duration_sec": float(raw.n_times / raw.info["sfreq"]),
        "stim_channel_expected": stim_channel,
        "stim_channel_present": stim_channel in raw.ch_names,
        "channel_type_counts": dict(Counter(channel_types)),
        "eog_channel_names": [
            name
            for name, kind in zip(raw.ch_names, channel_types, strict=True)
            if kind == "eog"
        ],
        "misc_channel_names": [
            name
            for name, kind in zip(raw.ch_names, channel_types, strict=True)
            if kind == "misc"
        ],
    }


# -----------------------------------------------------------------------------
# Inventory validation and optional QC rendering
# -----------------------------------------------------------------------------


def _validate_discovered_sources(
    cfg: Mapping[str, Any],
    sources: Mapping[str, list[Path]],
) -> tuple[list[str], list[str]]:
    """
    Validate required source categories without using prepared trial counts.

    Presence and uniqueness are checked for the subject-specific trajectory and
    TrialStruct files. Behavioral and trigger files only need to be non-empty;
    their saved contents determine presented trials, replay, and ESC handling.
    """
    errors: list[str] = []
    warnings: list[str] = []
    subject_dir = Path(cfg["resolved_paths"]["subject_data_dir"])

    if len(sources["task_meg_fif"]) == 0:
        errors.append(f"No task MEG FIF files found under {subject_dir / 'meg'}.")
    if len(sources["trajectory_mat"]) != 1:
        errors.append(
            "Expected exactly one subject-specific presented trajectory MAT; "
            f"found {len(sources['trajectory_mat'])}."
        )
    if len(sources["trial_struct_mat"]) != 1:
        errors.append(
            "Expected exactly one subject-specific TrialStruct MAT; "
            f"found {len(sources['trial_struct_mat'])}."
        )
    if not sources["behavioral_mat"]:
        errors.append("No saved behavioral MAT files were discovered.")
    if not sources["reported_trigger_csv"]:
        errors.append("No saved debug trigger CSV files were discovered.")

    ocular = cfg.get("subject_profile", {}).get("ocular_modalities", {})
    saved_edf_expected = bool(ocular.get("saved_edf", False))
    if saved_edf_expected and not sources["eyetracker_edf"]:
        errors.append("The subject profile expects saved EDF files, but none were found.")
    if not saved_edf_expected and sources["eyetracker_edf"]:
        warnings.append(
            "EDF files were found although the selected subject profile marks "
            "saved EDF data as unavailable."
        )
    if len(sources["empty_room_fif"]) > 1:
        errors.append(
            "More than one date-matched empty-room recording was discovered: "
            f"{[str(path) for path in sources['empty_room_fif']]}"
        )
    return errors, warnings


def build_inventory_payload(
    cfg: Mapping[str, Any],
    default_config_path: Path,
    subject_config_path: Path | None,
    warnings: list[str],
) -> dict[str, Any]:
    """
    Discover configured inputs and return a complete read-only preflight.

    The inventory remains optional QC data. Later stages can reproduce source
    discovery from configuration and do not consume this payload as a
    scientific input.
    """
    resolved = cfg["resolved_paths"]
    trigger_cfg = cfg.get("trigger_parser", {})
    eye_cfg = cfg.get("eye_model", {})
    sources = discover_subject_sources(cfg)
    errors, discovery_warnings = _validate_discovered_sources(cfg, sources)
    inventory_warnings = [*warnings, *discovery_warnings]

    stim_channel = str(trigger_cfg.get("stim_channel", "STI101"))
    fif_summaries: list[dict[str, Any]] = []
    for fif_path in sources["task_meg_fif"]:
        try:
            fif_summaries.append(summarize_fif_file(fif_path, stim_channel))
        except Exception as exc:  # pylint: disable=broad-except
            errors.append(f"Could not read FIF header {fif_path}: {exc}")

    stim_ok = bool(fif_summaries) and all(
        item["stim_channel_present"] for item in fif_summaries
    )
    if fif_summaries and not stim_ok:
        missing = [
            item["path"]
            for item in fif_summaries
            if not item["stim_channel_present"]
        ]
        errors.append(f"Stim channel {stim_channel!r} missing from: {missing}")

    eye_available = bool(sources["eyetracker_edf"])
    eye_enabled = bool(eye_cfg.get("enabled", False))
    eye_skip_if_missing = bool(eye_cfg.get("skip_if_missing", True))
    eye_will_skip = eye_enabled and eye_skip_if_missing and not eye_available
    status = "error" if errors else ("warning" if inventory_warnings else "ok")

    counts = {key: len(paths) for key, paths in sources.items()}
    # Compatibility aliases keep existing report readers stable.
    counts["meg_fif_files"] = counts["task_meg_fif"]
    counts["reported_trigger_csv_files"] = counts["reported_trigger_csv"]
    counts["eye_files"] = counts["eyetracker_edf"]

    files = {key: [str(path) for path in paths] for key, paths in sources.items()}
    files["meg_fif"] = fif_summaries
    files["eye"] = files["eyetracker_edf"]

    return {
        "inventory_schema_version": INVENTORY_SCHEMA_VERSION,
        "generated_utc": _iso_utc_now(),
        "status": status,
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "subject_profile": cfg.get("subject_profile", {}),
        "config_sources": {
            "default": str(default_config_path),
            "subject": str(subject_config_path) if subject_config_path else None,
        },
        "paths": {
            "data_root": resolved["data_root"],
            "subject_data_dir": resolved["subject_data_dir"],
            "derivatives_root": resolved["derivatives_root"],
        },
        "discovery": {
            "subject_folder_pattern": cfg["data_discovery"][
                "subject_folder_glob_template"
            ],
            "meg_fif_glob": cfg["data_discovery"]["meg_fif_glob"],
            "reported_trigger_csv_glob": cfg["data_discovery"][
                "reported_trigger_csv_glob"
            ],
            "stim_channel": stim_channel,
        },
        "checks": {
            "subject_data_dir_exists": Path(resolved["subject_data_dir"]).is_dir(),
            "has_meg_fif_files": bool(sources["task_meg_fif"]),
            "stim_channel_checks_passed": stim_ok,
            "subject_trajectory_unique": len(sources["trajectory_mat"]) == 1,
            "trial_struct_unique": len(sources["trial_struct_mat"]) == 1,
            "behavioral_mat_found": bool(sources["behavioral_mat"]),
            "reported_trigger_csv_found": bool(sources["reported_trigger_csv"]),
            "eye_data_available": eye_available,
            "empty_room_kept_separate": not bool(
                set(sources["task_meg_fif"]) & set(sources["empty_room_fif"])
            ),
        },
        "counts": counts,
        "eye_model_status": {
            "enabled": eye_enabled,
            "saved_edf_expected": bool(
                cfg.get("subject_profile", {})
                .get("ocular_modalities", {})
                .get("saved_edf", False)
            ),
            "skip_if_missing": eye_skip_if_missing,
            "will_skip_due_to_missing_eye_data": eye_will_skip,
            "skip_reason": (
                eye_cfg.get(
                    "skip_reason_when_missing",
                    "No configured eye-tracker input was found.",
                )
                if eye_will_skip
                else None
            ),
        },
        "files": files,
        "validation": {
            "is_valid": not errors,
            "warnings": inventory_warnings,
            "errors": errors,
        },
        "runtime": {
            "python_version": sys.version,
            "platform": platform.platform(),
        },
    }


def inventory_markdown(payload: Mapping[str, Any]) -> str:
    """Render a compact human-readable view of an inventory payload."""
    lines = [
        "# Dataset Inventory",
        "",
        f"- status: `{payload['status']}`",
        f"- subject: `{payload['subject_bids_id']}`",
        f"- session: `{payload['session_bids_id']}`",
        f"- generated_utc: `{payload['generated_utc']}`",
        f"- raw subject directory: `{payload['paths']['subject_data_dir']}`",
        f"- derivatives root: `{payload['paths']['derivatives_root']}`",
        "",
        "## Source Counts",
    ]
    for key in (
        "task_meg_fif",
        "empty_room_fif",
        "trajectory_mat",
        "trial_struct_mat",
        "behavioral_mat",
        "reported_trigger_csv",
        "eyetracker_edf",
    ):
        lines.append(f"- {key}: **{payload['counts'][key]}**")

    lines.extend(["", "## Checks"])
    lines.extend(f"- {key}: `{value}`" for key, value in payload["checks"].items())
    lines.extend(["", "## Task MEG Acquisition Runs"])
    if payload["files"]["meg_fif"]:
        for item in payload["files"]["meg_fif"]:
            lines.append(
                f"- `{item['path']}`: {item['sfreq_hz']} Hz, "
                f"{item['n_channels']} channels, EOG={item['eog_channel_names']}, "
                f"MISC={item['misc_channel_names']}"
            )
    else:
        lines.append("- none")

    lines.extend(["", "## Validation"])
    lines.append(
        "- warnings: "
        + (
            "; ".join(payload["validation"]["warnings"])
            if payload["validation"]["warnings"]
            else "none"
        )
    )
    lines.append(
        "- errors: "
        + (
            "; ".join(payload["validation"]["errors"])
            if payload["validation"]["errors"]
            else "none"
        )
    )
    lines.append("")
    return "\n".join(lines)


def write_inventory_reports(
    payload: Mapping[str, Any],
    *,
    out_json: str | Path | None,
    out_md: str | Path | None,
) -> dict[str, str]:
    """Write requested optional inventory formats outside the raw tree."""
    raw_root = payload["paths"]["data_root"]
    outputs: dict[str, str] = {}
    if out_json is not None:
        json_path = assert_output_path_outside_raw(out_json, raw_root)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        outputs["json"] = str(json_path)
    if out_md is not None:
        md_path = assert_output_path_outside_raw(out_md, raw_root)
        md_path.parent.mkdir(parents=True, exist_ok=True)
        md_path.write_text(inventory_markdown(payload), encoding="utf-8")
        outputs["markdown"] = str(md_path)
    return outputs


# -----------------------------------------------------------------------------
# Durable source manifest
# -----------------------------------------------------------------------------


def _source_record(path: Path, data_root: Path) -> dict[str, Any]:
    """Build a path-preserving source record with a raw-root-relative path."""
    return {
        "source_filename": path.name,
        "source_path": str(path),
        "raw_relative_path": str(path.relative_to(data_root)),
        "size_bytes": path.stat().st_size,
    }


def _extract_filename_number(name: str, label: str) -> int | None:
    """Extract a case-insensitive numeric filename label such as run or block."""
    match = re.search(rf"{re.escape(label)}(\d+)", name, flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def build_source_manifest_payload(
    cfg: Mapping[str, Any],
    default_config_path: Path,
    subject_config_path: Path | None,
    warnings: list[str],
) -> dict[str, Any]:
    """
    Build durable source lineage without asserting cross-modality alignment.

    Acquisition run labels are parsed only from MEG filenames. Behavioral
    block/run labels are parsed only from their own filenames. The association
    table remains empty until trigger-derived evidence is available.
    """
    sources = discover_subject_sources(cfg)
    resolved = cfg["resolved_paths"]
    data_root = Path(resolved["data_root"]).resolve()

    task_runs = []
    for path in sources["task_meg_fif"]:
        record = _source_record(path, data_root)
        record["acquisition_run_label"] = _extract_filename_number(path.name, "run")
        task_runs.append(record)

    behavior_files = []
    for path in sources["behavioral_mat"]:
        record = _source_record(path, data_root)
        record["behavioral_block_label"] = _extract_filename_number(
            path.name,
            "block",
        )
        behavior_files.append(record)

    trigger_files = []
    for path in sources["reported_trigger_csv"]:
        record = _source_record(path, data_root)
        record["behavioral_block_label"] = _extract_filename_number(
            path.name,
            "block",
        )
        record["within_block_run_label"] = _extract_filename_number(
            path.name,
            "run",
        )
        trigger_files.append(record)

    return {
        "source_manifest_schema_version": SOURCE_MANIFEST_SCHEMA_VERSION,
        "generated_utc": _iso_utc_now(),
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "subject_profile": cfg["subject_profile"],
        "config_sources": {
            "default": str(default_config_path),
            "subject": str(subject_config_path) if subject_config_path else None,
        },
        "raw_dataset": {
            "read_only": True,
            "data_root": str(data_root),
            "subject_folder_name": resolved["subject_folder_name"],
            "subject_folder_path": resolved["subject_data_dir"],
            "subject_folder_date": resolved["subject_folder_date"],
        },
        "sources": {
            "task_meg_acquisition_runs": task_runs,
            "empty_room_recordings": [
                _source_record(path, data_root)
                for path in sources["empty_room_fif"]
            ],
            "presented_trajectory_files": [
                _source_record(path, data_root)
                for path in sources["trajectory_mat"]
            ],
            "trial_struct_files": [
                _source_record(path, data_root)
                for path in sources["trial_struct_mat"]
            ],
            "behavioral_block_outputs": behavior_files,
            "debug_trigger_files": trigger_files,
            "eyetracker_edf_files": [
                _source_record(path, data_root)
                for path in sources["eyetracker_edf"]
            ],
        },
        "trigger_derived_associations": {
            "status": "pending",
            "method": (
                "Populate from observed MEG triggers and saved behavioral/debug "
                "records; never infer by positional list order."
            ),
            "meg_acquisition_run_to_behavioral_segments": [],
        },
        "warnings": list(warnings),
    }


def write_source_manifest(
    payload: Mapping[str, Any],
    path: str | Path,
) -> Path:
    """Write one durable source manifest after applying the shared raw guard."""
    raw_root = payload["raw_dataset"]["data_root"]
    resolved = assert_output_path_outside_raw(path, raw_root)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    resolved.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return resolved


# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the standalone read-only discovery and manifest CLI."""
    parser = argparse.ArgumentParser(description="Inventory DAD MEG subject inputs.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--out-json", default=None)
    parser.add_argument("--out-md", default=None)
    parser.add_argument("--source-manifest", default=None)
    parser.add_argument("--print-summary", action="store_true")
    return parser


def main() -> int:
    """Run discovery, optional FIF preflight, and requested serialization."""
    args = build_arg_parser().parse_args()
    try:
        cfg, warnings, default_path, subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            site_config_path=args.site_config,
        )
        payload = build_inventory_payload(
            cfg,
            default_path,
            subject_path,
            warnings,
        )
        write_inventory_reports(
            payload,
            out_json=args.out_json,
            out_md=args.out_md,
        )
        if args.source_manifest:
            manifest = build_source_manifest_payload(
                cfg,
                default_path,
                subject_path,
                warnings,
            )
            write_source_manifest(manifest, args.source_manifest)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Dataset inventory error: {exc}", file=sys.stderr)
        return 2

    if args.print_summary:
        print(f"Inventory status: {payload['status']}")
        print(f"Task MEG FIF files: {payload['counts']['task_meg_fif']}")
        print(f"Behavioral MAT files: {payload['counts']['behavioral_mat']}")
        print(f"Trigger CSV files: {payload['counts']['reported_trigger_csv']}")
        print(f"Eye EDF files: {payload['counts']['eyetracker_edf']}")
        print(f"Empty-room FIF files: {payload['counts']['empty_room_fif']}")
    return 0 if payload["validation"]["is_valid"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
