#!/usr/bin/env python3
"""
dataset_inventory.py

Purpose
Discover and validate subject inputs before expensive MEG processing. This is
a reusable preflight module, not a mandatory numbered transformation stage.
`run_pipeline.py` calls it automatically, and it can also be run directly when
an inventory report is useful.

Usage examples
1) Print a preflight summary without writing a report:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/dataset_inventory.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --print-summary

2) Write optional JSON and Markdown inventories:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/dataset_inventory.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --out-json /tmp/dataset_inventory.json \
     --out-md /tmp/dataset_inventory.md

Data flow
resolved configuration -> file discovery -> lazy FIF-header inspection ->
machine-readable inventory and preflight status.

Safety
All configured raw files are opened read-only. Report outputs are checked with
`assert_output_path_outside_raw()` before any directory is created.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import mne

from pipeline_config import (
    assert_output_path_outside_raw,
    load_and_resolve_config,
)


INVENTORY_SCHEMA_VERSION = "2.0"


# -----------------------------------------------------------------------------
# Discovery helpers
# -----------------------------------------------------------------------------


def _iso_utc_now() -> str:
    """Return a UTC ISO8601 timestamp with explicit Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def discover_meg_fif_files(
    subject_data_dir: str | Path,
    meg_fif_glob: str,
) -> list[Path]:
    """Return sorted unique FIF paths matching the configured subject glob."""
    root = Path(subject_data_dir).expanduser().resolve()
    return sorted(
        {path.resolve() for path in root.glob(meg_fif_glob) if path.is_file()}
    )


def discover_files(
    search_roots: Iterable[str | Path],
    *,
    glob_pattern: str | None = None,
    extensions: Iterable[str] = (),
) -> list[Path]:
    """
    Discover files across configured roots using a glob or extension filter.

    Trigger CSV discovery is normally non-recursive and glob-based. Eye files
    can be nested, so extension-based discovery is recursive.
    """
    discovered: set[Path] = set()
    normalized_exts = {str(ext).lower() for ext in extensions}
    for root_value in search_roots:
        root = Path(root_value).expanduser().resolve()
        if not root.is_dir():
            continue
        if glob_pattern is not None:
            discovered.update(
                path.resolve()
                for path in root.glob(glob_pattern)
                if path.is_file()
            )
        if normalized_exts:
            discovered.update(
                path.resolve()
                for path in root.rglob("*")
                if path.is_file() and path.suffix.lower() in normalized_exts
            )
    return sorted(discovered)


def summarize_fif_file(fif_path: str | Path, stim_channel: str) -> dict[str, Any]:
    """
    Read only the FIF header and summarize acquisition/channel metadata.

    The raw object is not preloaded, so this preflight does not read the full
    signal into memory or alter the source file.
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
    }


# -----------------------------------------------------------------------------
# Inventory construction
# -----------------------------------------------------------------------------


def build_inventory_payload(
    cfg: Mapping[str, Any],
    default_config_path: Path,
    subject_config_path: Path | None,
    warnings: list[str],
) -> dict[str, Any]:
    """
    Discover configured inputs and return a complete preflight payload.

    This signature intentionally remains straightforward for tests and other
    pipeline modules. The returned inventory is data, not a required derivative
    consumed by later stages.
    """
    resolved = cfg["resolved_paths"]
    discovery = cfg["data_discovery"]
    trigger_cfg = cfg.get("trigger_parser", {})
    eye_cfg = cfg.get("eye_model", {})

    subject_data_dir = Path(resolved["subject_data_dir"]).resolve()
    stim_channel = str(trigger_cfg.get("stim_channel", "STI101"))
    fif_files = discover_meg_fif_files(
        subject_data_dir,
        str(discovery["meg_fif_glob"]),
    )
    trigger_csvs = discover_files(
        discovery.get("reported_trigger_search_roots", []),
        glob_pattern=str(discovery["reported_trigger_csv_glob"]),
    )
    eye_files = discover_files(
        discovery.get("eye_search_glob_roots", []),
        extensions=discovery.get("eye_search_extensions", []),
    )

    fif_summaries: list[dict[str, Any]] = []
    errors: list[str] = []
    inventory_warnings = list(warnings)
    for fif_path in fif_files:
        try:
            fif_summaries.append(summarize_fif_file(fif_path, stim_channel))
        except Exception as exc:  # pylint: disable=broad-except
            errors.append(f"Could not read FIF header {fif_path}: {exc}")

    subject_dir_exists = subject_data_dir.is_dir()
    has_fifs = bool(fif_files)
    stim_ok = bool(fif_summaries) and all(
        item["stim_channel_present"] for item in fif_summaries
    )
    if not subject_dir_exists:
        errors.append(f"Subject data directory does not exist: {subject_data_dir}")
    if subject_dir_exists and not has_fifs:
        errors.append(
            f"No FIF files matched {discovery['meg_fif_glob']!r} in "
            f"{subject_data_dir}."
        )
    if fif_summaries and not stim_ok:
        missing = [
            item["path"]
            for item in fif_summaries
            if not item["stim_channel_present"]
        ]
        errors.append(f"Stim channel {stim_channel!r} missing from: {missing}")
    if not trigger_csvs:
        inventory_warnings.append(
            "No reported trigger CSV files found; trigger extraction may still "
            "use observed STI data, but reported-vs-observed QA is unavailable."
        )

    eye_available = bool(eye_files)
    eye_enabled = bool(eye_cfg.get("enabled", False))
    eye_skip_if_missing = bool(eye_cfg.get("skip_if_missing", True))
    eye_will_skip = eye_enabled and eye_skip_if_missing and not eye_available
    status = "error" if errors else ("warning" if inventory_warnings else "ok")

    return {
        "inventory_schema_version": INVENTORY_SCHEMA_VERSION,
        "generated_utc": _iso_utc_now(),
        "status": status,
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "config_sources": {
            "default": str(default_config_path),
            "subject": str(subject_config_path) if subject_config_path else None,
        },
        "paths": {
            "data_root": resolved["data_root"],
            "subject_data_dir": str(subject_data_dir),
            "derivatives_root": resolved["derivatives_root"],
        },
        "discovery": {
            "meg_fif_glob": discovery["meg_fif_glob"],
            "reported_trigger_csv_glob": discovery["reported_trigger_csv_glob"],
            "stim_channel": stim_channel,
        },
        "checks": {
            "subject_data_dir_exists": subject_dir_exists,
            "has_meg_fif_files": has_fifs,
            "stim_channel_checks_passed": stim_ok,
            "reported_trigger_csv_found": bool(trigger_csvs),
            "eye_data_available": eye_available,
        },
        "counts": {
            "meg_fif_files": len(fif_files),
            "reported_trigger_csv_files": len(trigger_csvs),
            "eye_files": len(eye_files),
        },
        "eye_model_status": {
            "enabled": eye_enabled,
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
        "files": {
            "meg_fif": fif_summaries,
            "reported_trigger_csv": [str(path) for path in trigger_csvs],
            "eye": [str(path) for path in eye_files],
        },
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
        "## Counts",
        f"- MEG FIF files: **{payload['counts']['meg_fif_files']}**",
        f"- reported trigger CSV files: **{payload['counts']['reported_trigger_csv_files']}**",
        f"- eye files: **{payload['counts']['eye_files']}**",
        "",
        "## Checks",
    ]
    lines.extend(
        f"- {key}: `{value}`" for key, value in payload["checks"].items()
    )
    lines.extend(["", "## MEG FIF Files"])
    if payload["files"]["meg_fif"]:
        for item in payload["files"]["meg_fif"]:
            lines.append(
                f"- `{item['path']}`: {item['sfreq_hz']} Hz, "
                f"{item['n_channels']} channels, EOG={item['eog_channel_names']}"
            )
    else:
        lines.append("- none")
    lines.extend(["", "## Validation"])
    if payload["validation"]["warnings"]:
        lines.append("- warnings:")
        lines.extend(f"  - {item}" for item in payload["validation"]["warnings"])
    else:
        lines.append("- warnings: none")
    if payload["validation"]["errors"]:
        lines.append("- errors:")
        lines.extend(f"  - {item}" for item in payload["validation"]["errors"])
    else:
        lines.append("- errors: none")
    lines.append("")
    return "\n".join(lines)


def write_inventory_reports(
    payload: Mapping[str, Any],
    *,
    out_json: str | Path | None,
    out_md: str | Path | None,
) -> dict[str, str]:
    """Write requested inventory formats after enforcing the raw-root guard."""
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
# Command-line interface
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the optional standalone preflight CLI."""
    parser = argparse.ArgumentParser(description="Inventory DAD MEG subject inputs.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--out-json", default=None)
    parser.add_argument("--out-md", default=None)
    parser.add_argument("--print-summary", action="store_true")
    return parser


def main() -> int:
    """Run preflight discovery and optionally serialize the inventory."""
    args = build_arg_parser().parse_args()
    try:
        cfg, warnings, default_path, subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
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
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Dataset inventory error: {exc}", file=sys.stderr)
        return 2

    if args.print_summary:
        print(f"Inventory status: {payload['status']}")
        print(f"MEG FIF files: {payload['counts']['meg_fif_files']}")
        print(f"Trigger CSV files: {payload['counts']['reported_trigger_csv_files']}")
        print(f"Eye files: {payload['counts']['eye_files']}")
    return 0 if payload["validation"]["is_valid"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
