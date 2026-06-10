#!/usr/bin/env python3
"""
02_trigger_extraction_custom.py

Purpose
Extract trigger events from MEG FIF files using custom logic aligned to this
experiment's trigger design, including transition-trigger resolution behavior.

This step reads the trigger channel from each FIF file, reconstructs non-zero
trigger windows, resolves contiguous transition pairs using expected-trigger-
first logic when expected trigger CSVs are available, and writes both machine-
readable outputs and human-readable companion reports.

The step also runs the deep subject-level reported-vs-observed confrontation
report from:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers_subject.py
so the same execution produces trigger derivatives for the pipeline and the
full trigger-QA report bundle.

Usage examples
1) Run with default report/output locations:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --print-summary

2) Run with explicit companion report paths:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --out-report-json /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/02_trigger_extraction_custom_report.json \
     --out-report-md /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/02_trigger_extraction_custom_report.md

3) Run with explicit deep-report output directory:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --deep-report-out-dir /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/02_triggers_extraction_custum_report_deeper

4) Run without deep confrontation report:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --disable-deep-trigger-report

Inputs
- Default YAML config and optional subject override YAML.
- FIF files discovered in subject data directory.
- Optional expected trigger CSV files from configured search roots.

Outputs
- Machine-readable event table:
  /derivatives/preprocessing/sub-01/ses-01/events/02_trigger_events.csv
  Includes both global-frame fields (`sample_start`, `sample_end`,
  `time_start_sec`, `time_end_sec`) and source-FIF-frame fields
  (`sample_start_raw`, `sample_end_raw`, `time_start_raw_sec`,
  `time_end_raw_sec`).
- Machine-readable event JSON payload:
  /derivatives/preprocessing/sub-01/ses-01/events/02_trigger_events.json
- Transition warning table:
  /derivatives/preprocessing/sub-01/ses-01/events/02_transition_warnings.csv
- Step companion reports (JSON + Markdown).
- Updated subject-level summary report (JSON + Markdown).
- Deep confrontation report bundle (alignment/anomalies/plots/Markdown):
  /derivatives/reports/sub-01/ses-01/02_triggers_extraction_custum_report_deeper/

Assumptions
- pipeline_config.py is available in the same folder and exposes the shared
  configuration and raw-data write-safety APIs.
- MNE-Python is available.
- triggers_check subject-level comparator script is available at the path above
  when deep report generation is enabled.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import subprocess
import sys
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import mne
import numpy as np

from pipeline_config import (
    assert_output_path_outside_raw,
    load_and_resolve_config,
)


STEP_ID = "02_trigger_extraction_custom"
STEP_SCRIPT_NAME = "02_trigger_extraction_custom.py"


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def _iso_utc_now() -> str:
    """Return current UTC timestamp in ISO8601 format with Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _sha256_for_file(path: Path) -> str:
    """Compute SHA256 checksum for provenance metadata."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _build_file_fingerprint(path: Path) -> dict[str, Any]:
    """Build file metadata payload with path/size/mtime/sha256."""
    if not path.exists():
        return {
            "path": str(path),
            "exists": False,
            "size_bytes": None,
            "mtime_utc": None,
            "sha256": None,
        }

    stat_info = path.stat()
    mtime_utc = datetime.fromtimestamp(stat_info.st_mtime, tz=timezone.utc)
    return {
        "path": str(path),
        "exists": True,
        "size_bytes": stat_info.st_size,
        "mtime_utc": mtime_utc.isoformat().replace("+00:00", "Z"),
        "sha256": _sha256_for_file(path),
    }


def _write_json_file(path: Path, payload: Any) -> None:
    """Write JSON payload to disk with stable formatting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _write_text_file(path: Path, text: str) -> None:
    """Write UTF-8 text content to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(text)


def _write_csv_file(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    """Write row dictionaries to CSV using provided field ordering."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _load_json_mapping_if_exists(path: Path) -> dict[str, Any]:
    """Load JSON mapping if path exists; otherwise return empty mapping."""
    if not path.exists():
        return {}
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # pylint: disable=broad-except
        return {}
    if not isinstance(parsed, dict):
        return {}
    return parsed


def _fmt_optional(value: Any) -> str:
    """Format optional value for Markdown display."""
    if value is None:
        return "-"
    return str(value)


def _extract_subject_number(subject_bids_id: str) -> Optional[int]:
    """Extract numeric subject identifier from BIDS-like subject label."""
    match = re.search(r"(\d+)", subject_bids_id)
    if not match:
        return None
    return int(match.group(1))


def run_deep_trigger_confrontation_report(
    *,
    cfg: dict[str, Any],
    subject_data_dir: Path,
    reported_csv_files: list[Path],
    derivatives_reports_dir: Path,
    output_dir_override: Optional[Path],
    variant: Optional[str],
) -> dict[str, Any]:
    """
    Run the subject-level deep reported-vs-observed trigger confrontation script.

    This keeps step-02 as the single entry point for both pipeline-ready trigger
    derivatives and the deeper trigger QA report bundle.
    """
    subject_bids_id = str(cfg["layout"]["subject_bids_id"])
    subject_num = _extract_subject_number(subject_bids_id)
    if subject_num is None:
        return {
            "status": "skipped",
            "reason": f"Could not parse numeric id from subject label: {subject_bids_id}",
            "enabled": True,
        }

    if not reported_csv_files:
        return {
            "status": "skipped",
            "reason": "No reported trigger CSV files found; deep confrontation cannot run.",
            "enabled": True,
        }

    project_root = Path(cfg["project"]["root"]).expanduser().resolve()
    deep_script_path = (
        project_root.parent
        / "triggers_check"
        / "compare_reported_vs_observed_triggers_subject.py"
    )
    if not deep_script_path.exists():
        return {
            "status": "skipped",
            "reason": f"Deep comparison script not found: {deep_script_path}",
            "enabled": True,
        }

    csv_parent_counts = Counter(path.parent.resolve() for path in reported_csv_files)
    csv_dir = sorted(
        csv_parent_counts.items(),
        key=lambda item: (-item[1], str(item[0])),
    )[0][0]

    deep_out_dir = (
        output_dir_override.resolve()
        if output_dir_override is not None
        else derivatives_reports_dir / "02_triggers_extraction_custum_report_deeper"
    )
    deep_out_dir.mkdir(parents=True, exist_ok=True)

    trigger_cfg = cfg["trigger_parser"]
    cmd = [
        sys.executable,
        str(deep_script_path),
        "--subject",
        str(subject_num),
        "--csv-dir",
        str(csv_dir),
        "--fif-dir",
        str(subject_data_dir),
        "--out-dir",
        str(deep_out_dir),
        "--trigger-channel",
        str(trigger_cfg["stim_channel"]),
        "--min-duration",
        str(float(trigger_cfg.get("min_duration_sec", 0.0))),
        "--transition-max-samples",
        str(int(trigger_cfg.get("transition_max_samples", 1))),
        "--expected-lookahead",
        str(int(trigger_cfg.get("expected_lookahead", 6))),
    ]
    if variant:
        cmd.extend(["--variant", variant])

    completed = subprocess.run(  # noqa: S603
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )

    stem = f"sub{subject_num:02d}_all_blocks_reported_vs_observed"
    outputs = {
        "alignment_csv": str((deep_out_dir / f"{stem}.csv").resolve()),
        "anomalies_csv": str((deep_out_dir / f"{stem}_anomalies.csv").resolve()),
        "report_md": str((deep_out_dir / f"{stem}.md").resolve()),
        "plots_dir": str((deep_out_dir / f"{stem}_plots").resolve()),
    }
    output_exists = {key: Path(path).exists() for key, path in outputs.items()}

    status = "ok" if completed.returncode == 0 else "error"
    reason = (
        "Deep confrontation report generated."
        if status == "ok"
        else f"Deep confrontation script failed with return code {completed.returncode}."
    )

    return {
        "enabled": True,
        "status": status,
        "reason": reason,
        "script_path": str(deep_script_path),
        "command": cmd,
        "return_code": completed.returncode,
        "csv_dir": str(csv_dir),
        "fif_dir": str(subject_data_dir),
        "out_dir": str(deep_out_dir),
        "outputs": outputs,
        "output_exists": output_exists,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-4000:],
    }


# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------


@dataclass
class ObservedWindow:
    """
    One observed non-zero trigger window reconstructed from FIF step events.

    Frame semantics
    - `start_sample` / `end_sample` and `start_sec` / `end_sec` are in the
      global concatenated frame used by step-02 alignment logic.
    - `start_sample_raw` / `end_sample_raw` and `start_sec_raw` /
      `end_sec_raw` are in the source FIF frame (includes `raw.first_samp`).
    """

    idx_global: int
    trigger: int
    start_sample: int
    end_sample: int
    start_sec: float
    end_sec: float
    start_sample_raw: int
    end_sample_raw: int
    start_sec_raw: float
    end_sec_raw: float
    width_samples: int
    source_fif: str
    source_fif_path: str


@dataclass
class ExpectedTrigger:
    """One expected trigger entry parsed from reported trigger CSV sequence."""

    idx_global: int
    trigger: int
    source_csv: str
    block: Optional[int]
    run: Optional[int]
    row_in_csv: int


# -----------------------------------------------------------------------------
# Trigger extraction with overlap support
# -----------------------------------------------------------------------------


def _merge_samples(min_duration: float, sfreq: float) -> int:
    """Convert minimum duration (seconds) into merge threshold samples."""
    min_samples = min_duration * sfreq
    if min_samples > 0:
        merge = int(min_samples // 1)
        if merge == min_samples:
            merge -= 1
        return merge
    return 0


def _find_stim_steps_1d(
    data: np.ndarray,
    first_samp: int,
    *,
    pad_stop: int | None = 0,
    merge: int = 0,
) -> np.ndarray:
    """Detect step changes in one-dimensional trigger data while keeping sample indexes."""
    data = np.asarray(data)
    if data.size == 0:
        return np.empty((0, 3), dtype=np.int64)

    changed = np.diff(data) != 0
    idx = np.where(changed)[0]
    if idx.size == 0:
        steps = np.empty((0, 3), dtype=np.int64)
    else:
        pre_step = data[idx]
        idx = idx + 1
        post_step = data[idx]
        idx = idx + first_samp
        steps = np.c_[idx, pre_step, post_step]

    if pad_stop is not None:
        last_val = data[-1]
        if last_val != pad_stop:
            last_idx = len(data) + first_samp
            steps = np.append(steps, [[last_idx, last_val, pad_stop]], axis=0)

    if merge != 0 and steps.size:
        diff = np.diff(steps[:, 0])
        idx_merge = diff <= abs(merge)
        if np.any(idx_merge):
            where = np.where(idx_merge)[0]
            keep = np.logical_not(idx_merge)
            if merge > 0:
                steps[where + 1, 1] = steps[where, 1]
                keep = np.append(keep, True)
            else:
                steps[where, 2] = steps[where + 1, 2]
                keep = np.insert(keep, 0, True)
            is_step = steps[:, 1] != steps[:, 2]
            keep = np.logical_and(keep, is_step)
            steps = steps[keep]

    return steps


def _mask_steps(steps: np.ndarray, mask: int | None, mask_type: str) -> np.ndarray:
    """Apply bitmask filtering to step values."""
    if steps.size == 0 or mask is None:
        return steps

    if mask_type == "not_and":
        mask = np.bitwise_not(mask)
    elif mask_type != "and":
        raise ValueError(f"mask_type must be 'and' or 'not_and', got {mask_type!r}")

    steps = steps.copy()
    steps[:, 1:] = np.bitwise_and(steps[:, 1:], mask)
    steps = steps[steps[:, 1] != steps[:, 2]]
    return steps


def _as_structured(arr: np.ndarray) -> np.ndarray:
    """Convert 2D numeric array to 1D structured view for row-wise set operations."""
    if arr.size == 0:
        return arr
    return np.ascontiguousarray(arr).view(
        np.dtype((np.void, arr.dtype.itemsize * arr.shape[1]))
    ).ravel()


def find_events_with_overlaps(
    raw: mne.io.BaseRaw,
    *,
    stim_channel: str,
    min_duration: float = 0.0,
    mask: int | None = None,
    mask_type: str = "and",
) -> np.ndarray:
    """
    Extract step events and explicitly preserve non-zero->non-zero overlap transitions.

    This mirrors the trigger-check utilities used in this project so transition
    pulses are represented before custom resolution logic is applied.
    """
    events = mne.find_events(
        raw,
        stim_channel=stim_channel,
        output="step",
        consecutive="increasing",
        min_duration=min_duration,
        shortest_event=1,
        mask=mask,
        mask_type=mask_type,
        verbose=False,
    )

    picks = mne.pick_channels(raw.info["ch_names"], include=[stim_channel], ordered=False)
    if len(picks) == 0:
        return events

    data = raw.get_data(picks=picks)
    merge = _merge_samples(min_duration, float(raw.info["sfreq"]))
    extra_steps = []

    for channel_data in data:
        if channel_data.size == 0:
            continue

        data_int = channel_data.astype(np.int64, copy=False)
        if data_int.min() < 0:
            data_int = np.abs(data_int)

        steps = _find_stim_steps_1d(data_int, raw.first_samp, pad_stop=0, merge=merge)
        steps = _mask_steps(steps, mask, mask_type)
        if steps.size == 0:
            continue

        overlap_mask = np.logical_and(steps[:, 1] > 0, steps[:, 2] > 0)
        overlap_steps = steps[overlap_mask]
        if overlap_steps.size:
            extra_steps.append(overlap_steps)

    if extra_steps:
        extra_arr = np.vstack(extra_steps).astype(np.int64)
        if events.size == 0:
            events = extra_arr
        else:
            merged = np.vstack([events.astype(np.int64), extra_arr])
            merged_view = _as_structured(merged)
            _, unique_idx = np.unique(merged_view, return_index=True)
            events = merged[np.sort(unique_idx)]

    if events.size:
        order = np.argsort(events[:, 0], kind="mergesort")
        events = events[order]

    return events


def reconstruct_observed_windows_for_fif(
    fif_path: Path,
    stim_channel: str,
    min_duration_sec: float,
    global_time_offset_sec: float,
    global_sample_offset: int,
) -> tuple[list[ObservedWindow], dict[str, Any]]:
    """
    Reconstruct observed non-zero trigger windows from one FIF file.

    Data flow
    - FIF trigger channel -> step events (with overlaps) -> non-zero windows.
    - Local sample/time -> global sample/time via offsets.

    Outputs
    - windows: list of observed windows with global coordinates.
    - file_summary: metadata/check payload for this FIF.
    """
    raw = mne.io.read_raw_fif(str(fif_path), preload=False, verbose="ERROR")
    ch_types = raw.get_channel_types()
    eog_names = [
        ch_name for ch_name, ch_type in zip(raw.ch_names, ch_types, strict=True) if ch_type == "eog"
    ]

    file_summary = {
        "path": str(fif_path),
        "sfreq_hz": float(raw.info["sfreq"]),
        "n_channels": int(raw.info["nchan"]),
        "n_samples": int(raw.n_times),
        "duration_sec": float(raw.n_times) / float(raw.info["sfreq"]),
        "stim_channel_expected": stim_channel,
        "stim_channel_present": stim_channel in raw.ch_names,
        "eog_channel_names": eog_names,
    }

    if stim_channel not in raw.ch_names:
        return [], file_summary

    events = find_events_with_overlaps(
        raw,
        stim_channel=stim_channel,
        min_duration=min_duration_sec,
    )

    sfreq = float(raw.info["sfreq"])
    first_samp = int(raw.first_samp)
    local_t0 = first_samp / sfreq

    windows: list[ObservedWindow] = []
    active_trigger: Optional[int] = None
    active_start: Optional[int] = None

    for sample, prev_val, new_val in events:
        sample_i = int(sample)
        prev_i = int(prev_val)
        new_i = int(new_val)

        if prev_i == 0 and new_i > 0:
            if active_trigger is not None and active_start is not None:
                global_start = global_sample_offset + (active_start - first_samp)
                global_end = global_sample_offset + (sample_i - first_samp)
                start_sec = global_time_offset_sec + ((active_start / sfreq) - local_t0)
                end_sec = global_time_offset_sec + ((sample_i / sfreq) - local_t0)
                windows.append(
                    ObservedWindow(
                        idx_global=0,
                        trigger=active_trigger,
                        start_sample=global_start,
                        end_sample=global_end,
                        start_sec=start_sec,
                        end_sec=end_sec,
                        start_sample_raw=active_start,
                        end_sample_raw=sample_i,
                        start_sec_raw=(active_start / sfreq),
                        end_sec_raw=(sample_i / sfreq),
                        width_samples=max(0, global_end - global_start),
                        source_fif=fif_path.name,
                        source_fif_path=str(fif_path),
                    )
                )
            active_trigger = new_i
            active_start = sample_i
        elif prev_i > 0 and new_i == 0:
            if active_trigger == prev_i and active_start is not None:
                global_start = global_sample_offset + (active_start - first_samp)
                global_end = global_sample_offset + (sample_i - first_samp)
                start_sec = global_time_offset_sec + ((active_start / sfreq) - local_t0)
                end_sec = global_time_offset_sec + ((sample_i / sfreq) - local_t0)
                windows.append(
                    ObservedWindow(
                        idx_global=0,
                        trigger=active_trigger,
                        start_sample=global_start,
                        end_sample=global_end,
                        start_sec=start_sec,
                        end_sec=end_sec,
                        start_sample_raw=active_start,
                        end_sample_raw=sample_i,
                        start_sec_raw=(active_start / sfreq),
                        end_sec_raw=(sample_i / sfreq),
                        width_samples=max(0, global_end - global_start),
                        source_fif=fif_path.name,
                        source_fif_path=str(fif_path),
                    )
                )
            active_trigger = None
            active_start = None
        elif prev_i > 0 and new_i > 0:
            if active_trigger == prev_i and active_start is not None:
                global_start = global_sample_offset + (active_start - first_samp)
                global_end = global_sample_offset + (sample_i - first_samp)
                start_sec = global_time_offset_sec + ((active_start / sfreq) - local_t0)
                end_sec = global_time_offset_sec + ((sample_i / sfreq) - local_t0)
                windows.append(
                    ObservedWindow(
                        idx_global=0,
                        trigger=active_trigger,
                        start_sample=global_start,
                        end_sample=global_end,
                        start_sec=start_sec,
                        end_sec=end_sec,
                        start_sample_raw=active_start,
                        end_sample_raw=sample_i,
                        start_sec_raw=(active_start / sfreq),
                        end_sec_raw=(sample_i / sfreq),
                        width_samples=max(0, global_end - global_start),
                        source_fif=fif_path.name,
                        source_fif_path=str(fif_path),
                    )
                )
            active_trigger = new_i
            active_start = sample_i

    return windows, file_summary


# -----------------------------------------------------------------------------
# Discovery helpers and expected-sequence parsing
# -----------------------------------------------------------------------------


def _parse_fif_block_bounds(filename: str) -> tuple[int, int]:
    """Parse block bounds from names like SUB01_block08_09_meg.fif or SUB01_block10_meg.fif."""
    match = re.search(r"_block(\d+)(?:_(\d+))?(?:_meg)?\.fif(?:\.gz)?$", filename, flags=re.IGNORECASE)
    if not match:
        return (10**9, 10**9)
    start = int(match.group(1))
    end = int(match.group(2)) if match.group(2) else start
    return (start, end)


def _parse_block_run_from_csv_name(filename: str) -> tuple[int, int]:
    """Parse block/run numbers from filenames containing blockXX_runYY tokens."""
    match = re.search(r"block(\d+)_run(\d+)", filename, flags=re.IGNORECASE)
    if not match:
        return (10**9, 10**9)
    return (int(match.group(1)), int(match.group(2)))


def discover_fif_files(subject_data_dir: Path, meg_fif_glob: str) -> list[Path]:
    """Discover and deterministically sort FIF files for trigger extraction."""
    files = [p.resolve() for p in subject_data_dir.glob(meg_fif_glob) if p.is_file()]
    return sorted(files, key=lambda p: (*_parse_fif_block_bounds(p.name), p.name))


def discover_reported_trigger_csvs(search_roots: list[str], csv_glob: str) -> list[Path]:
    """Discover and deterministically sort reported trigger CSV files."""
    found: set[Path] = set()
    for root_text in search_roots:
        root = Path(root_text).expanduser().resolve()
        if not root.exists() or not root.is_dir():
            continue
        for path in root.glob(csv_glob):
            if path.is_file():
                found.add(path.resolve())
    return sorted(found, key=lambda p: (*_parse_block_run_from_csv_name(p.name), p.name))


def load_expected_trigger_sequence(csv_paths: list[Path]) -> list[ExpectedTrigger]:
    """
    Load expected trigger stream by concatenating reported trigger CSVs in block/run order.

    Expected CSV schema requires at least a `trigger` column.
    """
    expected: list[ExpectedTrigger] = []
    next_idx = 1
    for csv_path in csv_paths:
        block, run = _parse_block_run_from_csv_name(csv_path.name)
        with csv_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "trigger" not in reader.fieldnames:
                continue
            for row_idx, row in enumerate(reader, start=1):
                trigger_text = (row.get("trigger") or "").strip()
                if not trigger_text:
                    continue
                try:
                    trigger_code = int(trigger_text)
                except ValueError:
                    continue
                expected.append(
                    ExpectedTrigger(
                        idx_global=next_idx,
                        trigger=trigger_code,
                        source_csv=str(csv_path),
                        block=(None if block >= 10**9 else block),
                        run=(None if run >= 10**9 else run),
                        row_in_csv=row_idx,
                    )
                )
                next_idx += 1
    return expected


# -----------------------------------------------------------------------------
# Transition resolution and event assembly
# -----------------------------------------------------------------------------


def trigger_label(code: int, seq_min: int, seq_max: int) -> str:
    """Map trigger code to canonical label used in machine-readable outputs."""
    fixed_map = {
        31: "trial_condition_onset_always_visible",
        41: "trial_condition_onset_occluded_nondeviant",
        51: "trial_condition_onset_occluded_deviant",
        111: "occlusion_start",
        112: "occlusion_complete",
        114: "occlusion_end_start",
        115: "occlusion_end_complete",
        116: "catch1_disappear",
        117: "catch1_reappear",
        118: "catch_question_start",
        201: "catch_question_correct",
        202: "catch_question_incorrect",
        113: "catch_question_timeout",
        150: "gaze_break",
        151: "replay_start",
        152: "esc_termination",
    }
    if code in fixed_map:
        return fixed_map[code]
    if seq_min <= code <= seq_max:
        return "trial_sequence_identity"
    return f"unknown_trigger_{code}"


def _find_expected_match_index(
    expected: list[ExpectedTrigger],
    start_idx: int,
    candidate_codes: set[int],
    max_lookahead: int,
) -> Optional[int]:
    """Find first expected index in lookahead window matching any candidate code."""
    if start_idx >= len(expected):
        return None
    end = min(len(expected), start_idx + max_lookahead + 1)
    for idx in range(start_idx, end):
        if expected[idx].trigger in candidate_codes:
            return idx
    return None


def resolve_observed_windows_to_events(
    observed_windows: list[ObservedWindow],
    expected_sequence: list[ExpectedTrigger],
    transition_max_samples: int,
    seq_min: int,
    seq_max: int,
    expected_lookahead: int = 6,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    """
    Resolve contiguous transition pairs and emit final event rows.

    Resolution policy
    - Primary: expected-trigger-first (drop non-intended, keep intended), symmetric by order.
    - Fallback: width heuristic when expected trigger is unavailable/ambiguous.
    - Final fallback: deterministic keep-first with unresolved warning.
    """
    events: list[dict[str, Any]] = []
    transition_warnings: list[dict[str, Any]] = []
    decision_counts: Counter[str] = Counter()

    i = 0
    expected_idx = 0
    next_event_idx = 1
    next_warning_idx = 1

    while i < len(observed_windows):
        cur = observed_windows[i]
        nxt = None
        if i + 1 < len(observed_windows):
            candidate = observed_windows[i + 1]
            if (
                candidate.source_fif_path == cur.source_fif_path
                and cur.end_sample == candidate.start_sample
            ):
                nxt = candidate

        if nxt is not None:
            pair_codes = {cur.trigger, nxt.trigger}
            match_idx = _find_expected_match_index(
                expected=expected_sequence,
                start_idx=expected_idx,
                candidate_codes=pair_codes,
                max_lookahead=expected_lookahead,
            )

            expected_code = None
            if match_idx is not None:
                expected_code = expected_sequence[match_idx].trigger

            keep = cur
            drop = nxt
            decision_source = "deterministic_fallback"
            unresolved = True
            confidence = 0.0

            if expected_code is not None and expected_code in pair_codes:
                if cur.trigger == expected_code:
                    keep, drop = cur, nxt
                else:
                    keep, drop = nxt, cur
                decision_source = "expected_trigger"
                unresolved = False
                confidence = 1.0
                expected_idx = match_idx + 1
            else:
                # Width fallback when expected identity is unavailable/ambiguous.
                cur_short = cur.width_samples <= transition_max_samples
                nxt_short = nxt.width_samples <= transition_max_samples
                if cur_short and not nxt_short:
                    keep, drop = nxt, cur
                    decision_source = "width_fallback"
                    unresolved = False
                    confidence = 0.7
                elif nxt_short and not cur_short:
                    keep, drop = cur, nxt
                    decision_source = "width_fallback"
                    unresolved = False
                    confidence = 0.7
                else:
                    # Deterministic fallback: keep first, drop second.
                    keep, drop = cur, nxt
                    decision_source = "deterministic_fallback"
                    unresolved = True
                    confidence = 0.0

            decision_counts[decision_source] += 1

            event_row = {
                "event_index": next_event_idx,
                "trigger_code": keep.trigger,
                "event_label": trigger_label(keep.trigger, seq_min, seq_max),
                "sample_start": keep.start_sample,
                "sample_end": keep.end_sample,
                "time_start_sec": keep.start_sec,
                "time_end_sec": keep.end_sec,
                "sample_start_raw": keep.start_sample_raw,
                "sample_end_raw": keep.end_sample_raw,
                "time_start_raw_sec": keep.start_sec_raw,
                "time_end_raw_sec": keep.end_sec_raw,
                "width_samples": keep.width_samples,
                "source_fif": keep.source_fif,
                "source_fif_path": keep.source_fif_path,
                "is_transition_warning": True,
                "transition_partner_code": drop.trigger,
                "expected_trigger_code": expected_code,
                "transition_decision_source": decision_source,
                "transition_resolution_confidence": confidence,
                "unresolved_transition": unresolved,
            }
            events.append(event_row)
            next_event_idx += 1

            transition_warnings.append(
                {
                    "warning_index": next_warning_idx,
                    "source_fif": cur.source_fif,
                    "source_fif_path": cur.source_fif_path,
                    "order_pattern": f"{cur.trigger}_then_{nxt.trigger}",
                    "dropped_code": drop.trigger,
                    "kept_code": keep.trigger,
                    "expected_trigger_code": expected_code,
                    "decision_source": decision_source,
                    "unresolved_transition": unresolved,
                    "transition_resolution_confidence": confidence,
                    "dropped_start_sample": drop.start_sample,
                    "dropped_end_sample": drop.end_sample,
                    "dropped_start_sec": drop.start_sec,
                    "dropped_end_sec": drop.end_sec,
                    "dropped_width_samples": drop.width_samples,
                    "kept_start_sample": keep.start_sample,
                    "kept_end_sample": keep.end_sample,
                    "kept_start_sec": keep.start_sec,
                    "kept_end_sec": keep.end_sec,
                }
            )
            next_warning_idx += 1

            i += 2
            continue

        # Non-transition window path.
        matched_idx = _find_expected_match_index(
            expected=expected_sequence,
            start_idx=expected_idx,
            candidate_codes={cur.trigger},
            max_lookahead=expected_lookahead,
        )
        expected_code = None
        if matched_idx is not None and expected_sequence[matched_idx].trigger == cur.trigger:
            expected_code = cur.trigger
            expected_idx = matched_idx + 1
            decision_counts["direct_expected_match"] += 1
        else:
            decision_counts["no_expected_match"] += 1

        events.append(
            {
                "event_index": next_event_idx,
                "trigger_code": cur.trigger,
                "event_label": trigger_label(cur.trigger, seq_min, seq_max),
                "sample_start": cur.start_sample,
                "sample_end": cur.end_sample,
                "time_start_sec": cur.start_sec,
                "time_end_sec": cur.end_sec,
                "sample_start_raw": cur.start_sample_raw,
                "sample_end_raw": cur.end_sample_raw,
                "time_start_raw_sec": cur.start_sec_raw,
                "time_end_raw_sec": cur.end_sec_raw,
                "width_samples": cur.width_samples,
                "source_fif": cur.source_fif,
                "source_fif_path": cur.source_fif_path,
                "is_transition_warning": False,
                "transition_partner_code": None,
                "expected_trigger_code": expected_code,
                "transition_decision_source": None,
                "transition_resolution_confidence": None,
                "unresolved_transition": False,
            }
        )
        next_event_idx += 1
        i += 1

    return events, transition_warnings, dict(decision_counts)


# -----------------------------------------------------------------------------
# Report builders
# -----------------------------------------------------------------------------


def build_step_report_payload(
    cfg: dict[str, Any],
    default_config_path: Path,
    subject_config_path: Optional[Path],
    fif_summaries: list[dict[str, Any]],
    discovered_csv_paths: list[Path],
    expected_sequence: list[ExpectedTrigger],
    events: list[dict[str, Any]],
    transition_warnings: list[dict[str, Any]],
    decision_counts: dict[str, int],
    deep_trigger_report: dict[str, Any],
) -> dict[str, Any]:
    """
    Build companion report payload for step 02 trigger extraction.

    Includes discovery, checks, event statistics, and transition-resolution
    provenance fields required by TRIGGERS_EVENTS.md.
    """
    trigger_cfg = cfg["trigger_parser"]
    required_codes = [int(code) for code in trigger_cfg["required_trigger_codes"]]
    observed_codes = [int(row["trigger_code"]) for row in events]
    observed_code_set = set(observed_codes)

    missing_required = sorted(code for code in required_codes if code not in observed_code_set)
    required_codes_passed = len(missing_required) == 0

    seq_min = int(trigger_cfg["expected_sequence_code_min"])
    seq_max = int(trigger_cfg["expected_sequence_code_max"])
    dynamic_codes_present = sorted(code for code in observed_code_set if seq_min <= code <= seq_max)

    stim_channel = trigger_cfg["stim_channel"]
    stim_checks_passed = bool(fif_summaries) and all(
        summary["stim_channel_present"] for summary in fif_summaries
    )

    warnings: list[str] = []
    errors: list[str] = []
    status = "ok"

    if not fif_summaries:
        status = "error"
        errors.append("No FIF files were discovered for trigger extraction.")

    if fif_summaries and not stim_checks_passed:
        status = "error"
        missing_files = [
            summary["path"] for summary in fif_summaries if not summary["stim_channel_present"]
        ]
        errors.append(
            f"Stim channel `{stim_channel}` missing in one or more FIF files: {missing_files}"
        )

    if missing_required:
        if status == "ok":
            status = "warning"
        warnings.append(
            "Some configured required trigger codes are absent in extracted events: "
            f"{missing_required}."
        )

    unresolved_count = sum(1 for row in events if bool(row["unresolved_transition"]))
    if unresolved_count > 0:
        if status == "ok":
            status = "warning"
        warnings.append(
            f"Unresolved transition fallback used for {unresolved_count} event(s)."
        )

    if not discovered_csv_paths:
        if status == "ok":
            status = "warning"
        warnings.append(
            "No reported trigger CSV files found; expected-trigger-first decisions may rely on fallback rules."
        )

    deep_status = str(deep_trigger_report.get("status", "skipped"))
    deep_reason = str(deep_trigger_report.get("reason", ""))
    if deep_status == "error":
        if status == "ok":
            status = "warning"
        warnings.append(f"Deep confrontation report failed: {deep_reason}")
    elif deep_status == "skipped":
        if status == "ok":
            status = "warning"
        warnings.append(f"Deep confrontation report skipped: {deep_reason}")

    decision_counter = Counter(decision_counts)

    return {
        "report_schema_version": "1.0",
        "step": {
            "step_id": STEP_ID,
            "script_name": STEP_SCRIPT_NAME,
            "status": status,
            "generated_utc": _iso_utc_now(),
        },
        "inputs": {
            "default_config": _build_file_fingerprint(default_config_path),
            "subject_config": (
                _build_file_fingerprint(subject_config_path) if subject_config_path else None
            ),
            "reported_trigger_csv_files": [str(path) for path in discovered_csv_paths],
        },
        "inventory": {
            "fif_files_count": len(fif_summaries),
            "reported_trigger_csv_count": len(discovered_csv_paths),
            "expected_sequence_length": len(expected_sequence),
            "stim_channel_expected": stim_channel,
        },
        "checks": {
            "stim_channel_present_all_fifs": stim_checks_passed,
            "required_trigger_codes_passed": required_codes_passed,
            "required_trigger_codes_missing": missing_required,
            "dynamic_sequence_range_min": seq_min,
            "dynamic_sequence_range_max": seq_max,
            "dynamic_sequence_codes_present": dynamic_codes_present,
            "transition_resolution_warnings_count": len(transition_warnings),
            "unresolved_transition_count": unresolved_count,
            "deep_trigger_report_status": deep_status,
        },
        "stats": {
            "raw_window_count": sum(summary.get("raw_window_count", 0) for summary in fif_summaries),
            "resolved_event_count": len(events),
            "decision_counts": dict(decision_counter),
            "trigger_code_counts": dict(Counter(observed_codes)),
        },
        "files": {
            "fif_summaries": fif_summaries,
            "deep_trigger_report": deep_trigger_report,
        },
        "validation": {
            "is_valid": status != "error",
            "warnings": warnings,
            "warning_count": len(warnings),
            "errors": errors,
            "error_count": len(errors),
        },
        "runtime": {
            "python_version": sys.version,
            "platform": platform.platform(),
            "argv": sys.argv,
        },
    }


def build_step_report_markdown(report_payload: dict[str, Any], output_paths: dict[str, str]) -> str:
    """Render step-02 companion report Markdown."""
    step = report_payload["step"]
    checks = report_payload["checks"]
    stats = report_payload["stats"]
    inv = report_payload["inventory"]
    validation = report_payload["validation"]
    deep_report = report_payload["files"].get("deep_trigger_report", {})

    lines: list[str] = []
    lines.append("# 02_trigger_extraction_custom Companion Report")
    lines.append("")
    lines.append("## Step")
    lines.append(f"- step_id: `{step['step_id']}`")
    lines.append(f"- status: `{step['status']}`")
    lines.append(f"- generated_utc: `{step['generated_utc']}`")
    lines.append("")

    lines.append("## Inventory")
    lines.append(f"- fif_files_count: **{inv['fif_files_count']}**")
    lines.append(f"- reported_trigger_csv_count: **{inv['reported_trigger_csv_count']}**")
    lines.append(f"- expected_sequence_length: **{inv['expected_sequence_length']}**")
    lines.append(f"- stim_channel_expected: `{inv['stim_channel_expected']}`")
    lines.append("")

    lines.append("## Checks")
    lines.append(f"- stim_channel_present_all_fifs: `{checks['stim_channel_present_all_fifs']}`")
    lines.append(f"- required_trigger_codes_passed: `{checks['required_trigger_codes_passed']}`")
    lines.append(
        f"- required_trigger_codes_missing: `{checks['required_trigger_codes_missing']}`"
    )
    lines.append(
        f"- dynamic_sequence_codes_present: `{checks['dynamic_sequence_codes_present']}`"
    )
    lines.append(
        f"- transition_resolution_warnings_count: `{checks['transition_resolution_warnings_count']}`"
    )
    lines.append(f"- unresolved_transition_count: `{checks['unresolved_transition_count']}`")
    lines.append("")

    lines.append("## Stats")
    lines.append(f"- raw_window_count: `{stats['raw_window_count']}`")
    lines.append(f"- resolved_event_count: `{stats['resolved_event_count']}`")
    lines.append(f"- decision_counts: `{stats['decision_counts']}`")
    lines.append("")

    lines.append("## Validation")
    lines.append(f"- is_valid: `{validation['is_valid']}`")
    lines.append(f"- warning_count: `{validation['warning_count']}`")
    lines.append(f"- error_count: `{validation['error_count']}`")
    if validation["warnings"]:
        lines.append("- warnings:")
        for warning in validation["warnings"]:
            lines.append(f"  - {warning}")
    else:
        lines.append("- warnings: none")
    if validation["errors"]:
        lines.append("- errors:")
        for err in validation["errors"]:
            lines.append(f"  - {err}")
    else:
        lines.append("- errors: none")
    lines.append("")

    lines.append("## Deep Trigger Confrontation Report")
    lines.append(f"- status: `{deep_report.get('status', '-')}`")
    lines.append(f"- reason: `{deep_report.get('reason', '-')}`")
    lines.append(f"- script_path: `{_fmt_optional(deep_report.get('script_path'))}`")
    lines.append(f"- csv_dir: `{_fmt_optional(deep_report.get('csv_dir'))}`")
    lines.append(f"- fif_dir: `{_fmt_optional(deep_report.get('fif_dir'))}`")
    lines.append(f"- out_dir: `{_fmt_optional(deep_report.get('out_dir'))}`")
    deep_outputs = deep_report.get("outputs", {})
    if isinstance(deep_outputs, dict) and deep_outputs:
        lines.append(f"- deep_alignment_csv: `{_fmt_optional(deep_outputs.get('alignment_csv'))}`")
        lines.append(f"- deep_anomalies_csv: `{_fmt_optional(deep_outputs.get('anomalies_csv'))}`")
        lines.append(f"- deep_report_md: `{_fmt_optional(deep_outputs.get('report_md'))}`")
        lines.append(f"- deep_plots_dir: `{_fmt_optional(deep_outputs.get('plots_dir'))}`")
    lines.append("")

    lines.append("## Output Files")
    lines.append(f"- events_csv: `{output_paths['events_csv']}`")
    lines.append(f"- events_json: `{output_paths['events_json']}`")
    lines.append(f"- transition_warnings_csv: `{output_paths['warnings_csv']}`")
    lines.append(f"- companion_report_json: `{output_paths['report_json']}`")
    lines.append(f"- companion_report_md: `{output_paths['report_md']}`")
    lines.append("")

    return "\n".join(lines)


def update_subject_summary_payload(
    *,
    cfg: dict[str, Any],
    report_payload: dict[str, Any],
    report_json_path: Path,
    report_md_path: Path,
    events_csv_path: Path,
    existing_summary_payload: dict[str, Any],
) -> dict[str, Any]:
    """Upsert step-02 compact status block inside subject-level summary payload."""
    updated = deepcopy(existing_summary_payload) if existing_summary_payload else {}
    steps = updated.get("steps")
    if not isinstance(steps, dict):
        steps = {}

    checks = report_payload["checks"]
    stats = report_payload["stats"]
    validation = report_payload["validation"]
    deep_report = report_payload.get("files", {}).get("deep_trigger_report", {})
    deep_outputs = deep_report.get("outputs", {}) if isinstance(deep_report, dict) else {}
    steps[STEP_ID] = {
        "status": report_payload["step"]["status"],
        "generated_utc": report_payload["step"]["generated_utc"],
        "required_trigger_codes_passed": checks["required_trigger_codes_passed"],
        "missing_required_codes": checks["required_trigger_codes_missing"],
        "stim_channel_checks_passed": checks["stim_channel_present_all_fifs"],
        "transition_resolution_warnings_count": checks["transition_resolution_warnings_count"],
        "unresolved_transition_count": checks["unresolved_transition_count"],
        "resolved_event_count": stats["resolved_event_count"],
        "warnings_count": validation["warning_count"],
        "errors_count": validation["error_count"],
        "companion_report_json": str(report_json_path),
        "companion_report_md": str(report_md_path),
        "events_csv": str(events_csv_path),
        "deep_trigger_report_status": deep_report.get("status"),
        "deep_trigger_report_reason": deep_report.get("reason"),
        "deep_trigger_report_md": deep_outputs.get("report_md"),
        "deep_trigger_alignment_csv": deep_outputs.get("alignment_csv"),
    }

    layout = cfg["layout"]
    updated["report_schema_version"] = "1.0"
    updated["subject_bids_id"] = layout["subject_bids_id"]
    updated["session_bids_id"] = layout["session_bids_id"]
    updated["updated_utc"] = _iso_utc_now()
    updated["steps"] = steps
    return updated


def build_subject_summary_markdown(summary_payload: dict[str, Any]) -> str:
    """Render subject-level summary Markdown including available step blocks."""
    steps = summary_payload.get("steps", {})

    lines: list[str] = []
    lines.append("# Subject Summary Report")
    lines.append("")
    lines.append(f"- subject: `{summary_payload.get('subject_bids_id', '-')}`")
    lines.append(f"- session: `{summary_payload.get('session_bids_id', '-')}`")
    lines.append(f"- updated_utc: `{summary_payload.get('updated_utc', '-')}`")
    lines.append("")

    lines.append("## Step Overview")
    lines.append("| step | status | key_check | warnings_count | errors_count | generated_utc |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for step_name in sorted(steps.keys()):
        info = steps[step_name]
        key_check = "-"
        if step_name == STEP_ID and "required_trigger_codes_passed" in info:
            key_check = (
                "required_trigger_codes_passed=true"
                if info.get("required_trigger_codes_passed")
                else "required_trigger_codes_passed=false"
            )
        lines.append(
            "| "
            + " | ".join(
                [
                    step_name,
                    str(info.get("status", "-")),
                    key_check,
                    str(info.get("warnings_count", "-")),
                    str(info.get("errors_count", "-")),
                    str(info.get("generated_utc", "-")),
                ]
            )
            + " |"
        )
    lines.append("")

    if STEP_ID in steps:
        info = steps[STEP_ID]
        lines.append(f"## {STEP_ID}")
        lines.append(
            f"- required_trigger_codes_passed: `{info.get('required_trigger_codes_passed', False)}`"
        )
        lines.append(f"- missing_required_codes: `{info.get('missing_required_codes', [])}`")
        lines.append(f"- stim_channel_checks_passed: `{info.get('stim_channel_checks_passed', False)}`")
        lines.append(
            f"- transition_resolution_warnings_count: `{info.get('transition_resolution_warnings_count', '-')}`"
        )
        lines.append(f"- unresolved_transition_count: `{info.get('unresolved_transition_count', '-')}`")
        lines.append(f"- resolved_event_count: `{info.get('resolved_event_count', '-')}`")
        lines.append(f"- companion_report_json: `{_fmt_optional(info.get('companion_report_json'))}`")
        lines.append(f"- companion_report_md: `{_fmt_optional(info.get('companion_report_md'))}`")
        lines.append(f"- events_csv: `{_fmt_optional(info.get('events_csv'))}`")
        lines.append(
            f"- deep_trigger_report_status: `{_fmt_optional(info.get('deep_trigger_report_status'))}`"
        )
        lines.append(
            f"- deep_trigger_report_reason: `{_fmt_optional(info.get('deep_trigger_report_reason'))}`"
        )
        lines.append(
            f"- deep_trigger_report_md: `{_fmt_optional(info.get('deep_trigger_report_md'))}`"
        )
        lines.append(
            f"- deep_trigger_alignment_csv: `{_fmt_optional(info.get('deep_trigger_alignment_csv'))}`"
        )
        lines.append("")

    return "\n".join(lines)


# -----------------------------------------------------------------------------
# CLI and main
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Build command-line parser for step-02 trigger extraction."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract custom trigger events from FIF files with expected-trigger-first "
            "transition resolution and report generation."
        )
    )
    parser.add_argument("--config", required=True, help="Path to default config YAML.")
    parser.add_argument(
        "--subject-config",
        default=None,
        help="Optional path to subject override config YAML.",
    )
    parser.add_argument(
        "--out-events-csv",
        default=None,
        help=(
            "Optional events CSV output path. Defaults to "
            "<derivatives_preprocessing_dir>/events/02_trigger_events.csv."
        ),
    )
    parser.add_argument(
        "--out-events-json",
        default=None,
        help=(
            "Optional events JSON output path. Defaults to "
            "<derivatives_preprocessing_dir>/events/02_trigger_events.json."
        ),
    )
    parser.add_argument(
        "--out-warnings-csv",
        default=None,
        help=(
            "Optional transition warnings CSV path. Defaults to "
            "<derivatives_preprocessing_dir>/events/02_transition_warnings.csv."
        ),
    )
    parser.add_argument(
        "--out-report-json",
        default=None,
        help=(
            "Optional companion report JSON path. Defaults to "
            "<derivatives_reports_dir>/02_trigger_extraction_custom_report.json."
        ),
    )
    parser.add_argument(
        "--out-report-md",
        default=None,
        help=(
            "Optional companion report Markdown path. Defaults to report JSON path "
            "with .md suffix."
        ),
    )
    parser.add_argument(
        "--out-subject-summary-json",
        default=None,
        help=(
            "Optional subject summary JSON path. Defaults to "
            "<derivatives_reports_dir>/<subject_bids_id>_summary_report.json."
        ),
    )
    parser.add_argument(
        "--out-subject-summary-md",
        default=None,
        help=(
            "Optional subject summary Markdown path. Defaults to summary JSON path "
            "with .md suffix."
        ),
    )
    parser.add_argument(
        "--disable-deep-trigger-report",
        action="store_true",
        help=(
            "Skip deep reported-vs-observed confrontation report generation "
            "from step 02."
        ),
    )
    parser.add_argument(
        "--deep-report-out-dir",
        default=None,
        help=(
            "Optional output directory for deep trigger confrontation report files. "
            "Defaults to <derivatives_reports_dir>/02_triggers_extraction_custum_report_deeper."
        ),
    )
    parser.add_argument(
        "--print-summary",
        action="store_true",
        help="Print compact step summary to stdout.",
    )
    return parser


def _print_compact_summary(report_payload: dict[str, Any], output_paths: dict[str, str]) -> None:
    """Print concise trigger extraction outcome."""
    step = report_payload["step"]
    checks = report_payload["checks"]
    stats = report_payload["stats"]
    deep = report_payload.get("files", {}).get("deep_trigger_report", {})
    print(f"Step {STEP_ID} status: {step['status']}")
    print(f"Resolved events: {stats['resolved_event_count']}")
    print(f"Transition warnings: {checks['transition_resolution_warnings_count']}")
    print(f"Required trigger coverage passed: {checks['required_trigger_codes_passed']}")
    print(f"Deep trigger report status: {deep.get('status', '-')}")
    print(f"Deep trigger report dir: {_fmt_optional(deep.get('out_dir'))}")
    print(f"Events CSV: {output_paths['events_csv']}")


def main() -> int:
    """Run trigger extraction, write outputs, and update reports."""
    args = build_arg_parser().parse_args()

    try:
        cfg, _, default_path, subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            required_sections=("trigger_parser",),
        )
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Trigger extraction error: invalid config: {exc}", file=sys.stderr)
        return 2

    resolved = cfg["resolved_paths"]
    trigger_cfg = cfg["trigger_parser"]

    subject_data_dir = Path(resolved["subject_data_dir"]).resolve()
    derivatives_preproc_dir = Path(resolved["derivatives_preprocessing_dir"]).resolve()
    derivatives_reports_dir = Path(resolved["derivatives_reports_dir"]).resolve()
    subject_bids_id = cfg["layout"]["subject_bids_id"]

    events_dir = derivatives_preproc_dir / "events"
    events_csv_path = (
        Path(args.out_events_csv).expanduser().resolve()
        if args.out_events_csv
        else events_dir / "02_trigger_events.csv"
    )
    events_json_path = (
        Path(args.out_events_json).expanduser().resolve()
        if args.out_events_json
        else events_dir / "02_trigger_events.json"
    )
    warnings_csv_path = (
        Path(args.out_warnings_csv).expanduser().resolve()
        if args.out_warnings_csv
        else events_dir / "02_transition_warnings.csv"
    )

    report_json_path = (
        Path(args.out_report_json).expanduser().resolve()
        if args.out_report_json
        else derivatives_reports_dir / f"{STEP_ID}_report.json"
    )
    report_md_path = (
        Path(args.out_report_md).expanduser().resolve()
        if args.out_report_md
        else report_json_path.with_suffix(".md")
    )

    summary_json_path = (
        Path(args.out_subject_summary_json).expanduser().resolve()
        if args.out_subject_summary_json
        else derivatives_reports_dir / f"{subject_bids_id}_summary_report.json"
    )
    summary_md_path = (
        Path(args.out_subject_summary_md).expanduser().resolve()
        if args.out_subject_summary_md
        else summary_json_path.with_suffix(".md")
    )
    deep_report_out_dir = (
        Path(args.deep_report_out_dir).expanduser().resolve()
        if args.deep_report_out_dir
        else (
            derivatives_reports_dir
            / str(trigger_cfg.get("deep_report_output_subdir", "02_triggers_extraction_custum_report_deeper"))
        ).resolve()
    )

    # Enforce the repository-wide contract that the source dataset is read-only.
    raw_root = Path(resolved["data_root"]).resolve()
    for output_path in (
        events_csv_path,
        events_json_path,
        warnings_csv_path,
        report_json_path,
        report_md_path,
        summary_json_path,
        summary_md_path,
        deep_report_out_dir,
    ):
        assert_output_path_outside_raw(output_path, raw_root)

    fif_files = discover_fif_files(subject_data_dir, cfg["data_discovery"]["meg_fif_glob"])
    reported_csv_files = discover_reported_trigger_csvs(
        cfg["data_discovery"].get("reported_trigger_search_roots", []),
        cfg["data_discovery"]["reported_trigger_csv_glob"],
    )
    expected_sequence = load_expected_trigger_sequence(reported_csv_files)

    all_windows: list[ObservedWindow] = []
    fif_summaries: list[dict[str, Any]] = []
    global_time_offset = 0.0
    global_sample_offset = 0
    next_window_idx = 1

    for fif_path in fif_files:
        windows, file_summary = reconstruct_observed_windows_for_fif(
            fif_path=fif_path,
            stim_channel=trigger_cfg["stim_channel"],
            min_duration_sec=float(trigger_cfg.get("min_duration_sec", 0.0)),
            global_time_offset_sec=global_time_offset,
            global_sample_offset=global_sample_offset,
        )

        for window in windows:
            window.idx_global = next_window_idx
            all_windows.append(window)
            next_window_idx += 1

        file_summary["raw_window_count"] = len(windows)
        fif_summaries.append(file_summary)

        # Advance offsets using file length even when stim channel is missing.
        global_time_offset += float(file_summary["duration_sec"])
        global_sample_offset += int(file_summary["n_samples"])

    seq_min = int(trigger_cfg["expected_sequence_code_min"])
    seq_max = int(trigger_cfg["expected_sequence_code_max"])
    transition_max_samples = int(trigger_cfg.get("transition_max_samples", 1))

    events, transition_warnings, decision_counts = resolve_observed_windows_to_events(
        observed_windows=all_windows,
        expected_sequence=expected_sequence,
        transition_max_samples=transition_max_samples,
        seq_min=seq_min,
        seq_max=seq_max,
    )

    events_csv_rows = [
        {
            "event_index": row["event_index"],
            "trigger_code": row["trigger_code"],
            "event_label": row["event_label"],
            "sample_start": row["sample_start"],
            "sample_end": row["sample_end"],
            "time_start_sec": f"{row['time_start_sec']:.6f}",
            "time_end_sec": f"{row['time_end_sec']:.6f}",
            # Source-FIF frame (includes raw.first_samp), used by step-04 epoching.
            "sample_start_raw": row["sample_start_raw"],
            "sample_end_raw": row["sample_end_raw"],
            "time_start_raw_sec": f"{row['time_start_raw_sec']:.6f}",
            "time_end_raw_sec": f"{row['time_end_raw_sec']:.6f}",
            "width_samples": row["width_samples"],
            "source_fif": row["source_fif"],
            "source_fif_path": row["source_fif_path"],
            "is_transition_warning": row["is_transition_warning"],
            "transition_partner_code": row["transition_partner_code"],
            "expected_trigger_code": row["expected_trigger_code"],
            "transition_decision_source": row["transition_decision_source"],
            "transition_resolution_confidence": row["transition_resolution_confidence"],
            "unresolved_transition": row["unresolved_transition"],
        }
        for row in events
    ]

    warnings_csv_rows = [
        {
            "warning_index": row["warning_index"],
            "source_fif": row["source_fif"],
            "source_fif_path": row["source_fif_path"],
            "order_pattern": row["order_pattern"],
            "dropped_code": row["dropped_code"],
            "kept_code": row["kept_code"],
            "expected_trigger_code": row["expected_trigger_code"],
            "decision_source": row["decision_source"],
            "unresolved_transition": row["unresolved_transition"],
            "transition_resolution_confidence": row["transition_resolution_confidence"],
            "dropped_start_sample": row["dropped_start_sample"],
            "dropped_end_sample": row["dropped_end_sample"],
            "dropped_start_sec": f"{row['dropped_start_sec']:.6f}",
            "dropped_end_sec": f"{row['dropped_end_sec']:.6f}",
            "dropped_width_samples": row["dropped_width_samples"],
            "kept_start_sample": row["kept_start_sample"],
            "kept_end_sample": row["kept_end_sample"],
            "kept_start_sec": f"{row['kept_start_sec']:.6f}",
            "kept_end_sec": f"{row['kept_end_sec']:.6f}",
        }
        for row in transition_warnings
    ]

    _write_csv_file(
        events_csv_path,
        fieldnames=list(events_csv_rows[0].keys()) if events_csv_rows else [
            "event_index",
            "trigger_code",
            "event_label",
            "sample_start",
            "sample_end",
            "time_start_sec",
            "time_end_sec",
            "sample_start_raw",
            "sample_end_raw",
            "time_start_raw_sec",
            "time_end_raw_sec",
            "width_samples",
            "source_fif",
            "source_fif_path",
            "is_transition_warning",
            "transition_partner_code",
            "expected_trigger_code",
            "transition_decision_source",
            "transition_resolution_confidence",
            "unresolved_transition",
        ],
        rows=events_csv_rows,
    )

    _write_csv_file(
        warnings_csv_path,
        fieldnames=list(warnings_csv_rows[0].keys()) if warnings_csv_rows else [
            "warning_index",
            "source_fif",
            "source_fif_path",
            "order_pattern",
            "dropped_code",
            "kept_code",
            "expected_trigger_code",
            "decision_source",
            "unresolved_transition",
            "transition_resolution_confidence",
            "dropped_start_sample",
            "dropped_end_sample",
            "dropped_start_sec",
            "dropped_end_sec",
            "dropped_width_samples",
            "kept_start_sample",
            "kept_end_sample",
            "kept_start_sec",
            "kept_end_sec",
        ],
        rows=warnings_csv_rows,
    )

    events_json_payload = {
        "generated_utc": _iso_utc_now(),
        "step_id": STEP_ID,
        "events": events,
        "transition_warnings": transition_warnings,
    }
    _write_json_file(events_json_path, events_json_payload)

    deep_enabled_cfg = bool(trigger_cfg.get("enable_deep_confrontation_report", True))
    deep_variant_cfg = trigger_cfg.get("deep_report_variant")
    deep_variant = str(deep_variant_cfg) if deep_variant_cfg else None

    if args.disable_deep_trigger_report:
        deep_trigger_report = {
            "enabled": False,
            "status": "skipped",
            "reason": "Disabled by --disable-deep-trigger-report.",
        }
    elif not deep_enabled_cfg:
        deep_trigger_report = {
            "enabled": False,
            "status": "skipped",
            "reason": "Disabled by trigger_parser.enable_deep_confrontation_report=false.",
        }
    else:
        deep_trigger_report = run_deep_trigger_confrontation_report(
            cfg=cfg,
            subject_data_dir=subject_data_dir,
            reported_csv_files=reported_csv_files,
            derivatives_reports_dir=derivatives_reports_dir,
            output_dir_override=deep_report_out_dir,
            variant=deep_variant,
        )

    step_report = build_step_report_payload(
        cfg=cfg,
        default_config_path=default_path,
        subject_config_path=subject_path,
        fif_summaries=fif_summaries,
        discovered_csv_paths=reported_csv_files,
        expected_sequence=expected_sequence,
        events=events,
        transition_warnings=transition_warnings,
        decision_counts=decision_counts,
        deep_trigger_report=deep_trigger_report,
    )

    output_paths = {
        "events_csv": str(events_csv_path),
        "events_json": str(events_json_path),
        "warnings_csv": str(warnings_csv_path),
        "report_json": str(report_json_path),
        "report_md": str(report_md_path),
    }

    _write_json_file(report_json_path, step_report)
    _write_text_file(report_md_path, build_step_report_markdown(step_report, output_paths))

    existing_summary = _load_json_mapping_if_exists(summary_json_path)
    updated_summary = update_subject_summary_payload(
        cfg=cfg,
        report_payload=step_report,
        report_json_path=report_json_path,
        report_md_path=report_md_path,
        events_csv_path=events_csv_path,
        existing_summary_payload=existing_summary,
    )
    _write_json_file(summary_json_path, updated_summary)
    _write_text_file(summary_md_path, build_subject_summary_markdown(updated_summary))

    if args.print_summary:
        _print_compact_summary(step_report, output_paths)

    if step_report["step"]["status"] == "error":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
