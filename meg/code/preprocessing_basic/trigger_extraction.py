#!/usr/bin/env python3
"""
trigger_extraction.py

Purpose
Implement reusable trigger reconstruction, transition handling, behavioral
segment parsing, and reported-vs-observed alignment for the DAD MEG pipeline.
The numbered scientific entry point is ``02_trigger_extraction_custom.py``.

Data flow
1. Resolve one subject profile through ``pipeline_config.py``.
2. Discover acquisition FIFs and exact runtime-specific debug CSV files.
3. Reconstruct every non-zero STI101 window in its source FIF/raw frame.
4. Resolve one-sample non-zero transition pulses using nearby reported codes
   as QC evidence, with width and deterministic fallbacks.
5. Align each saved behavioral block/run CSV segment in chronological order to
   the remaining observed event stream. Explicitly reviewed acquisition/block
   boundaries constrain both transition resolution and alignment when present.
6. Attach block/run/trial/replay metadata where alignment supports it; infer
   block/run for observed-only events only from unambiguous same-FIF context.
7. Render the deep Stage 1 report from the canonical alignment, keeping raw
   alignment statuses while separately marking next-frame matched gaze-break
   trials that the configured exclusion policy removes.
8. Write observed events, full alignment QC, source associations, reports, and
   an optional source-manifest update outside the read-only raw tree.

Usage examples
1) Run through the numbered stage:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-04.yaml

2) Run parser tests:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_trigger_parser.py -q

Key assumptions
- FIF filenames identify acquisition runs but do not define behavioral blocks.
- CSV filenames identify behavioral block and within-block run segments.
- CSV rows describe events actually emitted/saved, including replay and ESC;
  no prepared schedule count is used.
- Trigger 151 starts the replay phase for the remainder of that saved run
  segment, matching the runtime implementation that appends replay trials.
- Alignment metadata is QC/provenance. Observed STI101 events remain present
  even when no CSV row aligns to them.
- P37 has a user-confirmed truncated first acquisition. Its profile explicitly
  maps recordings to blocks so the unrecorded tail cannot consume later FIFs.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import sys
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

import mne
import numpy as np

from deep_trigger_report import generate_legacy_deep_report
from pipeline_config import assert_output_path_outside_raw, load_and_resolve_config


STEP_ID = "02_trigger_extraction_custom"
STEP_SCRIPT_NAME = "02_trigger_extraction_custom.py"
EVENTS_CSV_NAME = "02_trigger_events.csv"
EVENTS_JSON_NAME = "02_trigger_events.json"
ALIGNMENT_CSV_NAME = "02_trigger_alignment.csv"
WARNINGS_CSV_NAME = "02_transition_warnings.csv"
ASSOCIATIONS_JSON_NAME = "02_trigger_source_associations.json"


# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------


@dataclass
class ObservedWindow:
    """
    One non-zero STI101 window reconstructed from a source FIF.

    ``start_sample_raw`` and ``end_sample_raw`` include ``raw.first_samp`` and
    are the authoritative sample frame for later source-FIF processing.
    ``start_sample`` and ``end_sample`` are a deterministic subject-order frame
    retained for compatibility and cross-file reporting only.
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
    acquisition_run: Optional[int] = None


@dataclass
class ExpectedTrigger:
    """
    One saved debug-trigger row used as alignment/QC evidence.

    The first six fields preserve the legacy test/API contract. Additional
    fields retain actual presentation metadata without becoming event sources.
    """

    idx_global: int
    trigger: int
    source_csv: str
    block: Optional[int]
    run: Optional[int]
    row_in_csv: int
    trial: Optional[int] = None
    frame: Optional[int] = None
    seconds: Optional[float] = None
    label: str = ""
    replay_status: str = "original"
    is_replay_trial: bool = False
    segment_index: Optional[int] = None


@dataclass
class ReportedSegment:
    """One validated behavioral block/run CSV and its actual saved rows."""

    segment_index: int
    block: int
    run: int
    source_csv: Path
    events: list[ExpectedTrigger] = field(default_factory=list)
    termination_status: str = "saved_end"


@dataclass
class AlignmentOp:
    """One sequence-alignment operation using local segment/event indexes."""

    op_type: str
    reported_idx: Optional[int]
    observed_idx: Optional[int]


# -----------------------------------------------------------------------------
# Generic parsing and file helpers
# -----------------------------------------------------------------------------


def _iso_utc_now() -> str:
    """Return a UTC ISO8601 timestamp with a Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _optional_int(value: Any) -> Optional[int]:
    """Parse a possibly empty integer cell."""
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    return int(float(text))


def _optional_float(value: Any) -> Optional[float]:
    """Parse a possibly empty floating-point cell."""
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    return float(text)


def _parse_acquisition_run(filename: str) -> Optional[int]:
    """Parse the acquisition run label from an external task FIF filename."""
    match = re.search(r"(?:^|_)run(\d+)", filename, flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def _parse_block_run_from_csv_name(filename: str) -> tuple[int, int]:
    """Parse behavioral block and within-block run labels from a debug CSV."""
    match = re.search(r"block(\d+)_run(\d+)", filename, flags=re.IGNORECASE)
    if not match:
        raise ValueError(
            "Debug trigger CSV filename must contain blockNN_runNN: "
            f"{filename}"
        )
    return int(match.group(1)), int(match.group(2))


def _write_json(path: Path, payload: Any) -> None:
    """Write stable UTF-8 JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    """Write a machine-readable CSV even when the row list is empty."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def discover_fif_files(subject_data_dir: Path, meg_fif_glob: str) -> list[Path]:
    """Discover task FIFs and sort by acquisition run, never by block."""
    files = {
        path.resolve()
        for path in subject_data_dir.glob(meg_fif_glob)
        if path.is_file()
    }
    return sorted(
        files,
        key=lambda path: (
            _parse_acquisition_run(path.name) is None,
            _parse_acquisition_run(path.name) or 10**9,
            path.name,
        ),
    )


def discover_reported_trigger_csvs(
    search_roots: list[str],
    csv_glob: str,
) -> list[Path]:
    """Discover exact profile-derived debug CSV candidates."""
    found: set[Path] = set()
    for root_text in search_roots:
        root = Path(root_text).expanduser().resolve()
        if not root.is_dir():
            continue
        found.update(
            path.resolve() for path in root.glob(csv_glob) if path.is_file()
        )
    return sorted(
        found,
        key=lambda path: (*_parse_block_run_from_csv_name(path.name), path.name),
    )


def load_reported_trigger_segments(
    csv_paths: list[Path],
    *,
    subject_id: str,
    runtime_version: str,
    runtime_by_block: Optional[Mapping[int, str]] = None,
    replay_start_code: int = 151,
    esc_termination_code: int = 152,
) -> list[ReportedSegment]:
    """
    Parse debug CSVs as validated behavioral segments.

    Files are not treated as interchangeable chunks. Each unique block/run
    pair becomes one segment, filename runtime/subject tokens are verified, and
    every actual saved row is retained. Replay state begins after trigger 151
    within a segment. Trigger 152 marks an ESC-terminated segment. A block
    runtime mapping supports sessions that changed experiment version mid-run.
    """
    filename_re = re.compile(
        rf"^debug_actual_triggers_occlusion_(?P<runtime>v\d+_[A-Za-z0-9]+)"
        rf"_sub{re.escape(subject_id)}_block"
        r"(?P<block>\d+)_run(?P<run>\d+)\.csv$",
        flags=re.IGNORECASE,
    )
    seen: set[tuple[int, int]] = set()
    segments: list[ReportedSegment] = []
    next_expected_index = 1

    for csv_path in csv_paths:
        match = filename_re.match(csv_path.name)
        if not match:
            raise ValueError(
                "Debug trigger CSV does not match the selected subject profile "
                f"(runtime={runtime_version}, subject={subject_id}): {csv_path.name}"
            )
        block = int(match.group("block"))
        run = int(match.group("run"))
        selected_runtime = (
            runtime_by_block.get(block, runtime_version)
            if runtime_by_block is not None
            else runtime_version
        )
        expected_runtime = (
            f"{selected_runtime}_rescueTraject"
            if re.fullmatch(r"v\d+", selected_runtime, flags=re.IGNORECASE)
            else selected_runtime
        )
        if match.group("runtime").lower() != expected_runtime.lower():
            raise ValueError(
                "Debug trigger CSV runtime differs from the selected subject "
                f"profile for block {block}: expected {expected_runtime}, "
                f"found {match.group('runtime')} in {csv_path.name}."
            )
        key = (block, run)
        if key in seen:
            raise ValueError(
                f"Duplicate debug trigger segment for block {block}, run {run}."
            )
        seen.add(key)

        segment = ReportedSegment(
            segment_index=len(segments) + 1,
            block=block,
            run=run,
            source_csv=csv_path,
        )
        replay_phase = False
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            required = {"trigger", "trial", "frame", "seconds", "label"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Debug trigger CSV missing columns {sorted(missing)}: {csv_path}"
                )

            for row_in_csv, row in enumerate(reader, start=1):
                trigger = _optional_int(row.get("trigger"))
                if trigger is None:
                    continue
                if trigger == replay_start_code:
                    replay_status = "replay_start"
                    is_replay_trial = False
                    replay_phase = True
                else:
                    replay_status = "replay" if replay_phase else "original"
                    is_replay_trial = replay_phase

                segment.events.append(
                    ExpectedTrigger(
                        idx_global=next_expected_index,
                        trigger=trigger,
                        source_csv=str(csv_path),
                        block=block,
                        run=run,
                        row_in_csv=row_in_csv,
                        trial=_optional_int(row.get("trial")),
                        frame=_optional_int(row.get("frame")),
                        seconds=_optional_float(row.get("seconds")),
                        label=(row.get("label") or "").strip(),
                        replay_status=replay_status,
                        is_replay_trial=is_replay_trial,
                        segment_index=segment.segment_index,
                    )
                )
                next_expected_index += 1

        if not segment.events:
            raise ValueError(f"Debug trigger CSV has no trigger rows: {csv_path}")
        if any(event.trigger == esc_termination_code for event in segment.events):
            segment.termination_status = "esc_terminated"
        segments.append(segment)

    return sorted(segments, key=lambda item: (item.block, item.run))


def flatten_reported_segments(
    segments: Iterable[ReportedSegment],
) -> list[ExpectedTrigger]:
    """Return a boundary-preserving chronological view for transition lookup."""
    return [event for segment in segments for event in segment.events]


def load_expected_trigger_sequence(csv_paths: list[Path]) -> list[ExpectedTrigger]:
    """
    Legacy parser retained for focused transition tests.

    Production code uses ``load_reported_trigger_segments`` because runtime and
    subject validation are mandatory there.
    """
    expected: list[ExpectedTrigger] = []
    for csv_path in csv_paths:
        block, run = _parse_block_run_from_csv_name(csv_path.name)
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            for row_in_csv, row in enumerate(reader, start=1):
                trigger = _optional_int(row.get("trigger"))
                if trigger is None:
                    continue
                expected.append(
                    ExpectedTrigger(
                        idx_global=len(expected) + 1,
                        trigger=trigger,
                        source_csv=str(csv_path),
                        block=block,
                        run=run,
                        row_in_csv=row_in_csv,
                    )
                )
    return expected


# -----------------------------------------------------------------------------
# Trigger-channel reconstruction
# -----------------------------------------------------------------------------


def _merge_samples(min_duration: float, sfreq: float) -> int:
    """Convert minimum event duration to MNE's step-merge threshold."""
    min_samples = min_duration * sfreq
    if min_samples > 0:
        merge = int(min_samples // 1)
        return merge - 1 if merge == min_samples else merge
    return 0


def _find_stim_steps_1d(
    data: np.ndarray,
    first_samp: int,
    *,
    pad_stop: int | None = 0,
    merge: int = 0,
) -> np.ndarray:
    """Detect all value changes in one trigger channel."""
    data = np.asarray(data)
    if data.size == 0:
        return np.empty((0, 3), dtype=np.int64)
    changed = np.diff(data) != 0
    indexes = np.where(changed)[0]
    if indexes.size:
        previous = data[indexes]
        indexes = indexes + 1
        current = data[indexes]
        steps = np.c_[indexes + first_samp, previous, current]
    else:
        steps = np.empty((0, 3), dtype=np.int64)

    if pad_stop is not None and data[-1] != pad_stop:
        steps = np.append(
            steps,
            [[len(data) + first_samp, data[-1], pad_stop]],
            axis=0,
        )

    if merge != 0 and steps.size:
        differences = np.diff(steps[:, 0])
        merge_indexes = differences <= abs(merge)
        if np.any(merge_indexes):
            where = np.where(merge_indexes)[0]
            keep = np.logical_not(merge_indexes)
            if merge > 0:
                steps[where + 1, 1] = steps[where, 1]
                keep = np.append(keep, True)
            else:
                steps[where, 2] = steps[where + 1, 2]
                keep = np.insert(keep, 0, True)
            steps = steps[np.logical_and(keep, steps[:, 1] != steps[:, 2])]
    return steps.astype(np.int64, copy=False)


def find_events_with_overlaps(
    raw: mne.io.BaseRaw,
    *,
    stim_channel: str,
    min_duration: float = 0.0,
    mask: int | None = None,
    mask_type: str = "and",
) -> np.ndarray:
    """Return MNE step events plus explicit non-zero to non-zero transitions."""
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
    picks = mne.pick_channels(
        raw.info["ch_names"],
        include=[stim_channel],
        ordered=False,
    )
    if len(picks) == 0:
        return events

    data = raw.get_data(picks=picks)
    merge = _merge_samples(min_duration, float(raw.info["sfreq"]))
    extra_steps: list[np.ndarray] = []
    for channel_data in data:
        values = channel_data.astype(np.int64, copy=False)
        if values.size and values.min() < 0:
            values = np.abs(values)
        steps = _find_stim_steps_1d(
            values,
            raw.first_samp,
            pad_stop=0,
            merge=merge,
        )
        if mask is not None and steps.size:
            effective_mask = np.bitwise_not(mask) if mask_type == "not_and" else mask
            if mask_type not in ("and", "not_and"):
                raise ValueError(f"Unsupported mask_type: {mask_type}")
            steps[:, 1:] = np.bitwise_and(steps[:, 1:], effective_mask)
            steps = steps[steps[:, 1] != steps[:, 2]]
        overlaps = steps[(steps[:, 1] > 0) & (steps[:, 2] > 0)]
        if overlaps.size:
            extra_steps.append(overlaps)

    if extra_steps:
        combined = np.vstack([events.astype(np.int64), *extra_steps])
        events = np.unique(combined, axis=0)
        events = events[np.argsort(events[:, 0], kind="mergesort")]
    return events.astype(np.int64, copy=False)


def reconstruct_observed_windows_for_fif(
    fif_path: Path,
    stim_channel: str,
    min_duration_sec: float,
    global_time_offset_sec: float,
    global_sample_offset: int,
) -> tuple[list[ObservedWindow], dict[str, Any]]:
    """Read one FIF lazily and reconstruct all observed non-zero windows."""
    raw = mne.io.read_raw_fif(str(fif_path), preload=False, verbose="ERROR")
    channel_types = raw.get_channel_types()
    summary = {
        "path": str(fif_path),
        "source_fif": fif_path.name,
        "acquisition_run": _parse_acquisition_run(fif_path.name),
        "sfreq_hz": float(raw.info["sfreq"]),
        "first_samp": int(raw.first_samp),
        "n_samples": int(raw.n_times),
        "duration_sec": float(raw.n_times / raw.info["sfreq"]),
        "stim_channel_expected": stim_channel,
        "stim_channel_present": stim_channel in raw.ch_names,
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
    if stim_channel not in raw.ch_names:
        summary["raw_window_count"] = 0
        return [], summary

    steps = find_events_with_overlaps(
        raw,
        stim_channel=stim_channel,
        min_duration=min_duration_sec,
    )
    sfreq = float(raw.info["sfreq"])
    first_samp = int(raw.first_samp)
    windows: list[ObservedWindow] = []
    active_trigger: Optional[int] = None
    active_start: Optional[int] = None

    def append_window(end_sample_raw: int) -> None:
        """Close the current active pulse in both source and subject frames."""
        if active_trigger is None or active_start is None:
            return
        local_start = active_start - first_samp
        local_end = end_sample_raw - first_samp
        windows.append(
            ObservedWindow(
                idx_global=0,
                trigger=active_trigger,
                start_sample=global_sample_offset + local_start,
                end_sample=global_sample_offset + local_end,
                start_sec=global_time_offset_sec + (local_start / sfreq),
                end_sec=global_time_offset_sec + (local_end / sfreq),
                start_sample_raw=active_start,
                end_sample_raw=end_sample_raw,
                start_sec_raw=active_start / sfreq,
                end_sec_raw=end_sample_raw / sfreq,
                width_samples=max(0, end_sample_raw - active_start),
                source_fif=fif_path.name,
                source_fif_path=str(fif_path),
                acquisition_run=_parse_acquisition_run(fif_path.name),
            )
        )

    for sample, previous, current in steps:
        sample_i = int(sample)
        previous_i = int(previous)
        current_i = int(current)
        if previous_i == 0 and current_i > 0:
            append_window(sample_i)
            active_trigger = current_i
            active_start = sample_i
        elif previous_i > 0 and current_i == 0:
            append_window(sample_i)
            active_trigger = None
            active_start = None
        elif previous_i > 0 and current_i > 0:
            append_window(sample_i)
            active_trigger = current_i
            active_start = sample_i

    summary["raw_window_count"] = len(windows)
    return windows, summary


# -----------------------------------------------------------------------------
# Transition resolution
# -----------------------------------------------------------------------------


def trigger_label(
    code: int,
    seq_min: int,
    seq_max: int,
    condition_triggers: Optional[Mapping[str, int]] = None,
) -> str:
    """
    Map profile-derived condition, identity, and fixed codes to labels.

    When a subject profile supplies condition codes, those codes are
    authoritative. This avoids labeling v20 identity 31 as the legacy v18
    always-visible condition. Calls without a profile retain the historical
    cross-generation labels used by standalone parser diagnostics.
    """
    condition_labels = {
        "always_visible": "trial_condition_onset_always_visible",
        "occluded_nondeviant": "trial_condition_onset_occluded_nondeviant",
        "occluded_deviant": "trial_condition_onset_occluded_deviant",
    }
    if condition_triggers:
        for name, trigger_code in condition_triggers.items():
            if code == int(trigger_code):
                return condition_labels.get(name, f"trial_condition_onset_{name}")

    # Keep both generations recognizable in standalone parser tests.
    legacy_conditions = {
        31: "trial_condition_onset_always_visible",
        41: "trial_condition_onset_occluded_nondeviant",
        51: "trial_condition_onset_occluded_deviant",
        71: "trial_condition_onset_always_visible",
        81: "trial_condition_onset_occluded_nondeviant",
        91: "trial_condition_onset_occluded_deviant",
    }
    fixed = {
        111: "occlusion_start",
        112: "occlusion_complete",
        113: "catch_question_timeout",
        114: "occlusion_end_start",
        115: "occlusion_end_complete",
        116: "catch1_disappear",
        117: "catch1_reappear",
        118: "catch_question_start",
        150: "gaze_break",
        151: "replay_start",
        152: "esc_termination",
        201: "catch_question_correct",
        202: "catch_question_incorrect",
    }
    if code in fixed:
        return fixed[code]
    if seq_min <= code <= seq_max:
        return "trial_sequence_identity"
    if not condition_triggers and code in legacy_conditions:
        return legacy_conditions[code]
    return f"unknown_trigger_{code}"


def _find_expected_match_index(
    expected: list[ExpectedTrigger],
    start_idx: int,
    candidate_codes: set[int],
    max_lookahead: int,
) -> Optional[int]:
    """Find the nearest reported code matching one observed candidate."""
    stop = min(len(expected), start_idx + max_lookahead + 1)
    for index in range(start_idx, stop):
        if expected[index].trigger in candidate_codes:
            return index
    return None


def _event_row(
    window: ObservedWindow,
    *,
    event_index: int,
    seq_min: int,
    seq_max: int,
    condition_triggers: Optional[Mapping[str, int]],
    transition_partner_code: Optional[int] = None,
    expected_trigger_code: Optional[int] = None,
    transition_decision_source: Optional[str] = None,
    transition_resolution_confidence: Optional[float] = None,
    unresolved_transition: bool = False,
    gaze_break_code: int = 150,
    esc_termination_code: int = 152,
) -> dict[str, Any]:
    """Convert one authoritative observed window to a mutable event record."""
    label = trigger_label(
        window.trigger,
        seq_min,
        seq_max,
        condition_triggers,
    )
    return {
        "event_index": event_index,
        "trigger_code": window.trigger,
        "event_label": label,
        "sample_start": window.start_sample,
        "sample_end": window.end_sample,
        "time_start_sec": window.start_sec,
        "time_end_sec": window.end_sec,
        "sample_start_raw": window.start_sample_raw,
        "sample_end_raw": window.end_sample_raw,
        "raw_sample_frame": window.start_sample_raw,
        "time_start_raw_sec": window.start_sec_raw,
        "time_end_raw_sec": window.end_sec_raw,
        "width_samples": window.width_samples,
        "source_fif": window.source_fif,
        "source_fif_path": window.source_fif_path,
        "acquisition_run": window.acquisition_run,
        "is_transition_warning": transition_partner_code is not None,
        "transition_partner_code": transition_partner_code,
        "expected_trigger_code": expected_trigger_code,
        "transition_decision_source": transition_decision_source,
        "transition_resolution_confidence": transition_resolution_confidence,
        "unresolved_transition": unresolved_transition,
        "behavioral_block": None,
        "behavioral_run": None,
        "behavioral_trial": None,
        "behavioral_frame": None,
        "reported_seconds": None,
        "reported_label": None,
        "reported_source_csv": None,
        "reported_row_in_csv": None,
        "reported_alignment_status": "UNALIGNED_OBSERVED",
        "reported_alignment_confidence": 0.0,
        "behavioral_context_source": None,
        "replay_status": "unknown",
        "is_replay_trial": None,
        "is_gaze_break": window.trigger == gaze_break_code,
        "is_escape_termination": window.trigger == esc_termination_code,
    }


def resolve_observed_windows_to_events(
    observed_windows: list[ObservedWindow],
    expected_sequence: list[ExpectedTrigger],
    transition_max_samples: int,
    seq_min: int,
    seq_max: int,
    expected_lookahead: int = 6,
    condition_triggers: Optional[Mapping[str, int]] = None,
    gaze_break_code: int = 150,
    esc_termination_code: int = 152,
    reviewed_keep_both_transitions: Optional[list[Mapping[str, Any]]] = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    """
    Collapse contiguous non-zero transition pairs without synthesizing events.

    Reported codes are used first to choose which observed pulse to retain.
    Width evidence is the fallback. Ambiguous pairs keep the first pulse and
    remain explicitly unresolved. Exact, manually reviewed source/sample
    overrides can retain both members of a transition pair. All warnings use
    the same CSV schema; ``retained_both`` distinguishes these reviewed pairs
    from an actual drop despite the legacy kept/dropped field names.
    """
    events: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    decisions: Counter[str] = Counter()
    observed_index = 0
    expected_index = 0
    reviewed_keep_both_transitions = reviewed_keep_both_transitions or []

    while observed_index < len(observed_windows):
        current = observed_windows[observed_index]
        following: Optional[ObservedWindow] = None
        if observed_index + 1 < len(observed_windows):
            candidate = observed_windows[observed_index + 1]
            if (
                candidate.source_fif_path == current.source_fif_path
                and current.end_sample_raw == candidate.start_sample_raw
            ):
                following = candidate

        if following is None:
            matched = _find_expected_match_index(
                expected_sequence,
                expected_index,
                {current.trigger},
                expected_lookahead,
            )
            expected_code = None
            if matched is not None:
                expected_code = expected_sequence[matched].trigger
                expected_index = matched + 1
                decisions["direct_expected_match"] += 1
            else:
                decisions["no_expected_match"] += 1
            events.append(
                _event_row(
                    current,
                    event_index=len(events) + 1,
                    seq_min=seq_min,
                    seq_max=seq_max,
                    condition_triggers=condition_triggers,
                    expected_trigger_code=expected_code,
                    gaze_break_code=gaze_break_code,
                    esc_termination_code=esc_termination_code,
                )
            )
            observed_index += 1
            continue

        pair_codes = {current.trigger, following.trigger}
        matched = _find_expected_match_index(
            expected_sequence,
            expected_index,
            pair_codes,
            expected_lookahead,
        )
        expected_code = (
            expected_sequence[matched].trigger if matched is not None else None
        )

        # Keep both pulses only for an exact human-reviewed source/sample
        # override, and require the behavior log to contain them consecutively.
        keep_both_rule = next(
            (
                rule
                for rule in reviewed_keep_both_transitions
                if str(rule.get("source_fif")) == current.source_fif
                and int(rule.get("first_code", -1)) == current.trigger
                and int(rule.get("second_code", -1)) == following.trigger
                and int(rule.get("first_sample_raw", -1)) == current.start_sample_raw
                and int(rule.get("second_sample_raw", -1)) == following.start_sample_raw
            ),
            None,
        )
        if keep_both_rule is not None:
            if (
                matched is None
                or matched + 1 >= len(expected_sequence)
                or expected_sequence[matched].trigger != current.trigger
                or expected_sequence[matched + 1].trigger != following.trigger
            ):
                raise ValueError(
                    "Reviewed keep-both transition no longer matches two "
                    "consecutive reported triggers: "
                    f"{current.source_fif} samples "
                    f"{current.start_sample_raw}/{following.start_sample_raw}."
                )
            for window, partner, expected in (
                (current, following.trigger, expected_sequence[matched]),
                (following, current.trigger, expected_sequence[matched + 1]),
            ):
                events.append(
                    _event_row(
                        window,
                        event_index=len(events) + 1,
                        seq_min=seq_min,
                        seq_max=seq_max,
                        condition_triggers=condition_triggers,
                        transition_partner_code=partner,
                        expected_trigger_code=expected.trigger,
                        transition_decision_source="reviewed_keep_both",
                        transition_resolution_confidence=1.0,
                        unresolved_transition=False,
                        gaze_break_code=gaze_break_code,
                        esc_termination_code=esc_termination_code,
                    )
                )
            decisions["reviewed_keep_both"] += 1
            warnings.append(
                {
                    "warning_index": len(warnings) + 1,
                    "source_fif": current.source_fif,
                    "source_fif_path": current.source_fif_path,
                    "acquisition_run": current.acquisition_run,
                    "order_pattern": f"{current.trigger}_then_{following.trigger}",
                    "dropped_code": following.trigger,
                    "kept_code": current.trigger,
                    "retained_both": True,
                    "retained_codes": [current.trigger, following.trigger],
                    "expected_trigger_code": current.trigger,
                    "decision_source": "reviewed_keep_both",
                    "unresolved_transition": False,
                    "transition_resolution_confidence": 1.0,
                    "dropped_start_sample_raw": following.start_sample_raw,
                    "dropped_end_sample_raw": following.end_sample_raw,
                    "dropped_width_samples": following.width_samples,
                    "kept_start_sample_raw": current.start_sample_raw,
                    "kept_end_sample_raw": current.end_sample_raw,
                    "kept_width_samples": current.width_samples,
                    "review_reason": str(keep_both_rule["reason"]),
                }
            )
            expected_index = matched + 2
            observed_index += 2
            continue

        keep, drop = current, following
        source = "deterministic_fallback"
        confidence = 0.0
        unresolved = True

        if expected_code in pair_codes:
            keep, drop = (
                (current, following)
                if current.trigger == expected_code
                else (following, current)
            )
            source = "expected_trigger"
            confidence = 1.0
            unresolved = False
            expected_index = int(matched) + 1
        else:
            current_short = current.width_samples <= transition_max_samples
            following_short = following.width_samples <= transition_max_samples
            if current_short != following_short:
                keep, drop = (
                    (following, current)
                    if current_short
                    else (current, following)
                )
                source = "width_fallback"
                confidence = 0.7
                unresolved = False

        decisions[source] += 1
        events.append(
            _event_row(
                keep,
                event_index=len(events) + 1,
                seq_min=seq_min,
                seq_max=seq_max,
                condition_triggers=condition_triggers,
                transition_partner_code=drop.trigger,
                expected_trigger_code=expected_code,
                transition_decision_source=source,
                transition_resolution_confidence=confidence,
                unresolved_transition=unresolved,
                gaze_break_code=gaze_break_code,
                esc_termination_code=esc_termination_code,
            )
        )
        warnings.append(
            {
                "warning_index": len(warnings) + 1,
                "source_fif": current.source_fif,
                "source_fif_path": current.source_fif_path,
                "acquisition_run": current.acquisition_run,
                "order_pattern": f"{current.trigger}_then_{following.trigger}",
                "dropped_code": drop.trigger,
                "kept_code": keep.trigger,
                # Keep the export schema identical to reviewed-pair warnings,
                # including when the first reviewed pair appears later in a FIF.
                "retained_both": False,
                "retained_codes": [keep.trigger],
                "expected_trigger_code": expected_code,
                "decision_source": source,
                "unresolved_transition": unresolved,
                "transition_resolution_confidence": confidence,
                "dropped_start_sample_raw": drop.start_sample_raw,
                "dropped_end_sample_raw": drop.end_sample_raw,
                "dropped_width_samples": drop.width_samples,
                "kept_start_sample_raw": keep.start_sample_raw,
                "kept_end_sample_raw": keep.end_sample_raw,
                "kept_width_samples": keep.width_samples,
                "review_reason": "",
            }
        )
        observed_index += 2

    return events, warnings, dict(decisions)


# -----------------------------------------------------------------------------
# Segment-aware reported-vs-observed alignment
# -----------------------------------------------------------------------------


def semiglobal_align_codes(
    reported_codes: list[int],
    observed_codes: list[int],
    *,
    score_match: int = 3,
    score_mismatch: int = -3,
    score_gap: int = -2,
    allow_trailing_reported_gap: bool = False,
) -> tuple[list[AlignmentOp], int, int]:
    """
    Align one complete reported segment to any observed subsequence.

    Leading observed events are free because acquisition files can begin before
    a behavioral segment. Reported gaps are retained as QC rows and never
    create observed events. At a confirmed truncated recording's end, a free
    reported suffix allows the recorded prefix to end at the actual cutoff.
    This prevents repeated later codes from pulling a missing tail into the
    middle of an otherwise matching prefix.
    """
    n_reported = len(reported_codes)
    n_observed = len(observed_codes)
    scores = np.zeros((n_reported + 1, n_observed + 1), dtype=np.int32)
    moves = np.zeros((n_reported + 1, n_observed + 1), dtype=np.uint8)
    for row in range(1, n_reported + 1):
        scores[row, 0] = scores[row - 1, 0] + score_gap
        moves[row, 0] = 2
    moves[0, 1:] = 4

    for row in range(1, n_reported + 1):
        reported_code = reported_codes[row - 1]
        for column in range(1, n_observed + 1):
            diagonal = scores[row - 1, column - 1] + (
                score_match
                if reported_code == observed_codes[column - 1]
                else score_mismatch
            )
            up = scores[row - 1, column] + score_gap
            left = scores[row, column - 1] + score_gap
            best = max(diagonal, up, left)
            scores[row, column] = best
            moves[row, column] = 1 if best == diagonal else (2 if best == up else 3)

    if allow_trailing_reported_gap:
        end_row = int(np.argmax(scores[:, n_observed]))
        end_column = n_observed
    else:
        end_row = n_reported
        end_column = int(np.argmax(scores[n_reported, :]))
    best_score = int(scores[end_row, end_column])
    row, column = end_row, end_column
    reversed_ops: list[AlignmentOp] = []
    while row > 0:
        move = int(moves[row, column])
        if move == 1:
            reversed_ops.append(AlignmentOp("PAIR", row - 1, column - 1))
            row -= 1
            column -= 1
        elif move == 2:
            reversed_ops.append(AlignmentOp("REPORTED_ONLY", row - 1, None))
            row -= 1
        elif move == 3:
            reversed_ops.append(AlignmentOp("OBSERVED_ONLY", None, column - 1))
            column -= 1
        elif move == 4:
            column -= 1
        else:
            raise RuntimeError(
                f"Unexpected alignment backtrack at row={row}, column={column}."
            )
    reversed_ops.reverse()
    reversed_ops.extend(AlignmentOp("REPORTED_ONLY", index, None)
                        for index in range(end_row, n_reported))
    return reversed_ops, best_score, end_column


def align_reported_segments_to_events(
    events: list[dict[str, Any]],
    segments: list[ReportedSegment],
    *,
    alignment_slack_events: int = 250,
    recording_truncated: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """
    Align each actual saved CSV segment to the remaining observed event stream.

    Behavioral metadata is attached to paired observed rows. Multi-block FIFs
    naturally receive multiple segment associations. A partial final segment
    consumes only its saved rows. An observed ESC 152 missing from the saved
    CSV remains OBSERVED_ONLY but can inherit an unambiguous preceding segment
    label from the same FIF. ``recording_truncated`` permits a free reported
    suffix only when fewer observed events remain than reported segment rows.
    """
    alignment_rows: list[dict[str, Any]] = []
    observed_cursor = 0
    paired_observed: set[int] = set()
    accounted_observed: set[int] = set()
    segment_summaries: list[dict[str, Any]] = []

    for segment in segments:
        reported = segment.events
        max_window = min(
            len(events) - observed_cursor,
            max(len(reported) * 3, len(reported) + alignment_slack_events),
        )
        observed_slice = events[observed_cursor : observed_cursor + max_window]
        operations, score, end_index = semiglobal_align_codes(
            [item.trigger for item in reported],
            [int(item["trigger_code"]) for item in observed_slice],
            allow_trailing_reported_gap=(recording_truncated and len(observed_slice) < len(reported)),
        )

        matched_count = 0
        mismatch_count = 0
        reported_only_count = 0
        observed_only_count = 0
        associated_fifs: Counter[tuple[str, Optional[int]]] = Counter()

        for operation in operations:
            reported_item = (
                reported[operation.reported_idx]
                if operation.reported_idx is not None
                else None
            )
            observed_global_index = (
                observed_cursor + operation.observed_idx
                if operation.observed_idx is not None
                else None
            )
            observed_item = (
                events[observed_global_index]
                if observed_global_index is not None
                else None
            )
            if observed_global_index is not None:
                accounted_observed.add(observed_global_index)

            status = operation.op_type
            if reported_item is not None and observed_item is not None:
                is_match = reported_item.trigger == int(observed_item["trigger_code"])
                if is_match:
                    matched_count += 1
                    status = (
                        "MATCH_WARNING"
                        if observed_item["is_transition_warning"]
                        else "MATCH"
                    )
                else:
                    mismatch_count += 1
                    status = "MISMATCH"

                observed_item.update(
                    {
                        "behavioral_block": reported_item.block,
                        "behavioral_run": reported_item.run,
                        "behavioral_trial": reported_item.trial,
                        "behavioral_frame": reported_item.frame,
                        "reported_seconds": reported_item.seconds,
                        "reported_label": reported_item.label,
                        "reported_source_csv": reported_item.source_csv,
                        "reported_row_in_csv": reported_item.row_in_csv,
                        "reported_alignment_status": status,
                        "reported_alignment_confidence": 1.0 if is_match else 0.0,
                        "behavioral_context_source": "reported_alignment",
                        "replay_status": reported_item.replay_status,
                        "is_replay_trial": reported_item.is_replay_trial,
                    }
                )
                paired_observed.add(observed_global_index)
                if is_match:
                    associated_fifs[
                        (
                            str(observed_item["source_fif"]),
                            observed_item["acquisition_run"],
                        )
                    ] += 1
            elif reported_item is not None:
                reported_only_count += 1
            elif observed_item is not None:
                observed_only_count += 1

            alignment_rows.append(
                {
                    "alignment_index": len(alignment_rows) + 1,
                    "status": status,
                    "behavioral_block": (
                        reported_item.block if reported_item is not None else None
                    ),
                    "behavioral_run": (
                        reported_item.run if reported_item is not None else None
                    ),
                    "behavioral_trial": (
                        reported_item.trial if reported_item is not None else None
                    ),
                    "reported_trigger_code": (
                        reported_item.trigger if reported_item is not None else None
                    ),
                    "reported_label": (
                        reported_item.label if reported_item is not None else None
                    ),
                    "reported_frame": (
                        reported_item.frame if reported_item is not None else None
                    ),
                    "reported_seconds": (
                        reported_item.seconds if reported_item is not None else None
                    ),
                    "reported_source_csv": (
                        reported_item.source_csv if reported_item is not None else None
                    ),
                    "reported_row_in_csv": (
                        reported_item.row_in_csv if reported_item is not None else None
                    ),
                    "reported_replay_status": (
                        reported_item.replay_status
                        if reported_item is not None
                        else None
                    ),
                    "observed_event_index": (
                        observed_item["event_index"]
                        if observed_item is not None
                        else None
                    ),
                    "observed_trigger_code": (
                        observed_item["trigger_code"]
                        if observed_item is not None
                        else None
                    ),
                    "source_fif": (
                        observed_item["source_fif"]
                        if observed_item is not None
                        else None
                    ),
                    "acquisition_run": (
                        observed_item["acquisition_run"]
                        if observed_item is not None
                        else None
                    ),
                    "raw_sample_frame": (
                        observed_item["raw_sample_frame"]
                        if observed_item is not None
                        else None
                    ),
                    "behavioral_context_source": (
                        observed_item["behavioral_context_source"]
                        if observed_item is not None
                        else None
                    ),
                }
            )

        denominator = max(1, len(reported))
        segment_summaries.append(
            {
                "segment_index": segment.segment_index,
                "behavioral_block": segment.block,
                "behavioral_run": segment.run,
                "source_csv": str(segment.source_csv),
                "termination_status": segment.termination_status,
                "reported_event_count": len(reported),
                "matched_event_count": matched_count,
                "mismatch_count": mismatch_count,
                "reported_only_count": reported_only_count,
                "observed_only_count": observed_only_count,
                "match_fraction": matched_count / denominator,
                "alignment_score": score,
                "associated_sources": [
                    {
                        "source_fif": source_fif,
                        "acquisition_run": acquisition_run,
                        "matched_event_count": count,
                    }
                    for (source_fif, acquisition_run), count in sorted(
                        associated_fifs.items(),
                        key=lambda item: (
                            item[0][1] is None,
                            item[0][1] or 10**9,
                            item[0][0],
                        ),
                    )
                ],
            }
        )
        observed_cursor += end_index

    # Observed events without a paired CSV row remain scientific events. Infer
    # behavioral context only when an exact same-FIF neighbor makes it
    # unambiguous; the event remains OBSERVED_ONLY QC evidence.
    for observed_index, observed_item in enumerate(events):
        if observed_index in accounted_observed:
            continue
        preceding = next(
            (
                candidate
                for candidate in reversed(events[:observed_index])
                if candidate["source_fif_path"] == observed_item["source_fif_path"]
                and candidate["reported_alignment_status"]
                in {"MATCH", "MATCH_WARNING"}
            ),
            None,
        )
        following = next(
            (
                candidate
                for candidate in events[observed_index + 1 :]
                if candidate["source_fif_path"] == observed_item["source_fif_path"]
                and candidate["reported_alignment_status"]
                in {"MATCH", "MATCH_WARNING"}
            ),
            None,
        )
        context = None
        context_source = None
        preceding_segment = (
            (preceding["behavioral_block"], preceding["behavioral_run"])
            if preceding is not None
            else None
        )
        following_segment = (
            (following["behavioral_block"], following["behavioral_run"])
            if following is not None
            else None
        )
        if (
            preceding is not None
            and following is not None
            and preceding_segment == following_segment
        ):
            context = preceding
            context_source = "inferred_same_segment_neighbors"
        elif observed_item["is_escape_termination"] and preceding is not None:
            context = preceding
            context_source = "inferred_esc_from_preceding_segment"
        elif preceding is not None and following is None:
            context = preceding
            context_source = "inferred_trailing_same_fif"
        elif following is not None and preceding is None:
            context = following
            context_source = "inferred_leading_same_fif"

        if context is not None:
            observed_item.update(
                {
                    "behavioral_block": context["behavioral_block"],
                    "behavioral_run": context["behavioral_run"],
                    "replay_status": context["replay_status"],
                    "is_replay_trial": context["is_replay_trial"],
                    "behavioral_context_source": context_source,
                }
            )
            for segment_summary in segment_summaries:
                if (
                    segment_summary["behavioral_block"]
                    == context["behavioral_block"]
                    and segment_summary["behavioral_run"]
                    == context["behavioral_run"]
                ):
                    segment_summary["observed_only_count"] += 1
                    if observed_item["is_escape_termination"]:
                        segment_summary[
                            "termination_status"
                        ] = "observed_esc_after_saved_rows"
                    break

        alignment_rows.append(
            {
                "alignment_index": len(alignment_rows) + 1,
                "status": "OBSERVED_ONLY",
                "behavioral_block": observed_item["behavioral_block"],
                "behavioral_run": observed_item["behavioral_run"],
                "behavioral_trial": None,
                "reported_trigger_code": None,
                "reported_label": None,
                "reported_frame": None,
                "reported_seconds": None,
                "reported_source_csv": None,
                "reported_row_in_csv": None,
                "reported_replay_status": None,
                "observed_event_index": observed_item["event_index"],
                "observed_trigger_code": observed_item["trigger_code"],
                "source_fif": observed_item["source_fif"],
                "acquisition_run": observed_item["acquisition_run"],
                "raw_sample_frame": observed_item["raw_sample_frame"],
                "behavioral_context_source": observed_item[
                    "behavioral_context_source"
                ],
            }
        )

    status_counts = Counter(row["status"] for row in alignment_rows)
    return events, alignment_rows, {
        "segments": segment_summaries,
        "status_counts": dict(status_counts),
        "paired_observed_event_count": len(paired_observed),
        "unpaired_observed_event_count": len(events) - len(paired_observed),
    }


def process_trigger_windows(
    windows: list[ObservedWindow],
    segments: list[ReportedSegment],
    *,
    trigger_config: Mapping[str, Any],
    subject_profile: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int],
           list[dict[str, Any]], dict[str, Any]]:
    """Resolve transitions and align within explicit acquisition boundaries.

    With no ``acquisition_block_map`` the existing subject-wide behavior is
    retained, including multi-block FIFs. A reviewed map associates each
    acquisition run with a disjoint list of behavioral blocks. Both the
    expected-code cursor and the alignment cursor restart at each boundary.
    Missing reported rows never create observed events. A known truncated
    acquisition labels only the reported suffix after its last observed event
    as unavailable recording tail; internal missing rows remain unexplained.
    Raw sample frames and source paths pass through unchanged.
    """
    block_map = subject_profile.get("acquisition_block_map")
    truncated = set(subject_profile.get("truncated_acquisition_runs", []))
    if block_map:
        observed_runs = {window.acquisition_run for window in windows}
        if observed_runs != set(block_map):
            raise ValueError(
                "Acquisition/block map must exactly cover recorded runs; "
                f"configured={sorted(block_map)}, observed={observed_runs}."
            )
        configured_blocks = [block for blocks in block_map.values() for block in blocks]
        if len(set(configured_blocks)) != len(configured_blocks) or set(configured_blocks) != {s.block for s in segments}:
            raise ValueError("Acquisition/block map must cover each saved behavioral block exactly once.")
        groups = [
            (run, [window for window in windows if window.acquisition_run == run],
             [segment for segment in segments if segment.block in blocks])
            for run, blocks in sorted(block_map.items())
        ]
    else:
        groups = [(None, windows, segments)]

    # Each group produces canonical events, not synthetic fillers. Reindex
    # before alignment so references remain unique across all source FIFs.
    events: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    decisions: Counter[str] = Counter()
    alignment_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    boundaries: list[dict[str, Any]] = []
    paired_count = 0
    seq_min, seq_max = [int(value) for value in subject_profile["identity_range"]]
    for run, selected_windows, selected_segments in groups:
        selected_events, selected_warnings, selected_decisions = resolve_observed_windows_to_events(
            selected_windows, flatten_reported_segments(selected_segments),
            int(trigger_config.get("transition_max_samples", 1)), seq_min, seq_max,
            int(trigger_config.get("expected_lookahead", 6)),
            subject_profile["condition_triggers"],
            int(trigger_config["gaze_break_code"]), int(trigger_config["esc_termination_code"]),
            trigger_config.get("reviewed_keep_both_transitions", []),
        )
        for index, row in enumerate(selected_events, start=len(events) + 1):
            row["event_index"] = index
        for index, row in enumerate(selected_warnings, start=len(warnings) + 1):
            row["warning_index"] = index
        selected_events, selected_alignment, summary = align_reported_segments_to_events(
            selected_events, selected_segments,
            alignment_slack_events=int(trigger_config.get("alignment_slack_events", 250)),
            recording_truncated=run in truncated,
        )

        # Classification follows the observed stream's end, independently of
        # expected trial counts. It does not explain an earlier isolated gap.
        tail_count = 0
        if run is not None:
            last_observed = max((index for index, row in enumerate(selected_alignment)
                                 if row["observed_event_index"] is not None), default=-1)
            for index, row in enumerate(selected_alignment):
                row["constrained_acquisition_run"] = run
                is_tail = run in truncated and index > last_observed and row["status"] == "REPORTED_ONLY"
                row["reported_only_reason"] = "recording_tail_unavailable" if is_tail else None
                tail_count += int(is_tail)
            boundaries.append({
                "acquisition_run": run,
                "behavioral_blocks": list(block_map[run]),
                "declared_truncated": run in truncated,
                "recording_tail_unavailable_count": tail_count,
                "alignment_status_counts": summary["status_counts"],
            })
        for index, row in enumerate(selected_alignment, start=len(alignment_rows) + 1):
            row["alignment_index"] = index
        events.extend(selected_events)
        warnings.extend(selected_warnings)
        decisions.update(selected_decisions)
        alignment_rows.extend(selected_alignment)
        summaries.extend(summary["segments"])
        paired_count += summary["paired_observed_event_count"]
    return events, warnings, dict(decisions), alignment_rows, {
        "segments": summaries,
        "status_counts": dict(Counter(row["status"] for row in alignment_rows)),
        "paired_observed_event_count": paired_count,
        "unpaired_observed_event_count": len(events) - paired_count,
        "acquisition_boundaries": boundaries,
    }


def build_source_associations(
    events: list[dict[str, Any]],
    alignment_summary: Mapping[str, Any],
) -> dict[str, Any]:
    """Aggregate evidence-backed behavioral segments for each acquisition FIF."""
    per_source: dict[tuple[str, str, Optional[int]], dict[str, Any]] = {}
    for event in events:
        key = (
            str(event["source_fif_path"]),
            str(event["source_fif"]),
            event["acquisition_run"],
        )
        record = per_source.setdefault(
            key,
            {
                "source_fif_path": key[0],
                "source_fif": key[1],
                "acquisition_run": key[2],
                "observed_event_count": 0,
                "aligned_event_count": 0,
                "behavioral_segments": defaultdict(int),
            },
        )
        record["observed_event_count"] += 1
        block = event["behavioral_block"]
        run = event["behavioral_run"]
        alignment_status = event["reported_alignment_status"]
        if (
            block is not None
            and run is not None
            and alignment_status in {"MATCH", "MATCH_WARNING"}
        ):
            record["aligned_event_count"] += 1
            record["behavioral_segments"][(int(block), int(run))] += 1

    associations: list[dict[str, Any]] = []
    for record in per_source.values():
        segments = [
            {
                "behavioral_block": block,
                "behavioral_run": run,
                "aligned_event_count": count,
            }
            for (block, run), count in sorted(record["behavioral_segments"].items())
        ]
        associations.append(
            {
                **{key: value for key, value in record.items() if key != "behavioral_segments"},
                "behavioral_segments": segments,
                "association_status": "aligned" if segments else "unresolved",
            }
        )

    associations.sort(
        key=lambda item: (
            item["acquisition_run"] is None,
            item["acquisition_run"] or 10**9,
            item["source_fif"],
        )
    )
    return {
        "schema_version": "1.0",
        "generated_utc": _iso_utc_now(),
        "method": (
            "Observed STI101 events aligned to validated runtime-specific saved "
            "debug-trigger segments; no positional FIF/block zipping."
        ),
        "alignment_status_counts": alignment_summary["status_counts"],
        "meg_acquisition_run_to_behavioral_segments": associations,
    }


# -----------------------------------------------------------------------------
# Reports and source-manifest update
# -----------------------------------------------------------------------------


def build_report(
    *,
    cfg: Mapping[str, Any],
    fif_summaries: list[dict[str, Any]],
    csv_paths: list[Path],
    segments: list[ReportedSegment],
    events: list[dict[str, Any]],
    transition_warnings: list[dict[str, Any]],
    decision_counts: Mapping[str, int],
    alignment_summary: Mapping[str, Any],
    manifest_update: Mapping[str, Any],
    deep_report_result: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the compact machine-readable stage report."""
    trigger_cfg = cfg["trigger_parser"]
    profile = cfg["subject_profile"]
    policy = profile["trigger_policy"]
    observed_codes = Counter(int(event["trigger_code"]) for event in events)
    required_codes = [int(code) for code in policy["required_subject_codes"]]
    missing_required = [code for code in required_codes if code not in observed_codes]
    seq_min, seq_max = [int(value) for value in profile["identity_range"]]
    out_of_policy = sorted(
        code
        for code in observed_codes
        if not (
            seq_min <= code <= seq_max
            or code in {int(value) for value in policy["allowed_fixed_codes"]}
            or code in {int(value) for value in profile["condition_triggers"].values()}
        )
    )
    errors: list[str] = []
    warnings: list[str] = []
    if not fif_summaries:
        errors.append("No task FIF files were discovered.")
    if any(not item["stim_channel_present"] for item in fif_summaries):
        errors.append("STI101 is missing from one or more task FIF files.")
    if not csv_paths:
        warnings.append("No runtime-specific debug trigger CSV files were discovered.")
    if missing_required:
        warnings.append(
            f"Profile-required subject codes are absent: {missing_required}."
        )
    if out_of_policy:
        warnings.append(f"Observed trigger codes outside profile policy: {out_of_policy}.")
    if alignment_summary["status_counts"].get("MISMATCH", 0):
        warnings.append("Reported-vs-observed trigger mismatches remain.")
    if alignment_summary["status_counts"].get("REPORTED_ONLY", 0):
        warnings.append("Some saved reported rows have no observed STI101 pulse.")
    if alignment_summary["status_counts"].get("OBSERVED_ONLY", 0):
        warnings.append(
            "Some observed STI101 events have no saved reported CSV row."
        )
    unresolved_transitions = sum(
        bool(item["unresolved_transition"]) for item in transition_warnings
    )
    if unresolved_transitions:
        warnings.append(
            f"{unresolved_transitions} transition pulse(s) remain unresolved."
        )

    status = "error" if errors else ("warning" if warnings else "ok")
    return {
        "report_schema_version": "2.0",
        "step": {
            "step_id": STEP_ID,
            "script_name": STEP_SCRIPT_NAME,
            "status": status,
            "generated_utc": _iso_utc_now(),
        },
        "subject": {
            "subject_bids_id": cfg["layout"]["subject_bids_id"],
            "runtime_version": profile["versions"]["runtime"],
            "condition_triggers": profile["condition_triggers"],
            "identity_range": profile["identity_range"],
            "reported_trigger_csv_pattern": profile[
                "reported_trigger_csv_pattern"
            ],
            "trigger_policy": policy,
            "acquisition_block_map": profile.get("acquisition_block_map"),
            "truncated_acquisition_runs": profile.get("truncated_acquisition_runs", []),
        },
        "inventory": {
            "task_fif_count": len(fif_summaries),
            "debug_csv_count": len(csv_paths),
            "behavioral_segment_count": len(segments),
            "reported_event_count": sum(len(segment.events) for segment in segments),
            "observed_raw_window_count": sum(
                item["raw_window_count"] for item in fif_summaries
            ),
            "resolved_observed_event_count": len(events),
        },
        "checks": {
            "stim_channel": trigger_cfg["stim_channel"],
            "stim_channel_present_all_fifs": bool(fif_summaries)
            and all(item["stim_channel_present"] for item in fif_summaries),
            "missing_required_subject_codes": missing_required,
            "out_of_policy_codes": out_of_policy,
            "transition_warning_count": len(transition_warnings),
            "unresolved_transition_count": unresolved_transitions,
            "alignment_status_counts": alignment_summary["status_counts"],
        },
        "stats": {
            "trigger_code_counts": dict(observed_codes),
            "transition_decision_counts": dict(decision_counts),
            "segment_alignment": alignment_summary["segments"],
            "acquisition_boundaries": alignment_summary.get("acquisition_boundaries", []),
        },
        "files": {
            "fif_summaries": fif_summaries,
            "debug_csv_files": [str(path) for path in csv_paths],
            "source_manifest_update": dict(manifest_update),
            "deep_trigger_report": dict(deep_report_result),
        },
        "validation": {
            "is_valid": not errors,
            "warnings": warnings,
            "errors": errors,
        },
        "runtime": {
            "python_version": sys.version,
            "platform": platform.platform(),
            "argv": sys.argv,
        },
    }


def report_markdown(
    report: Mapping[str, Any],
    outputs: Mapping[str, str],
) -> str:
    """Render a concise companion Markdown report."""
    step = report["step"]
    subject = report["subject"]
    inventory = report["inventory"]
    checks = report["checks"]
    validation = report["validation"]
    lines = [
        "# Trigger Extraction Report",
        "",
        f"- subject: `{subject['subject_bids_id']}`",
        f"- runtime: `{subject['runtime_version']}`",
        f"- status: `{step['status']}`",
        f"- generated_utc: `{step['generated_utc']}`",
        "",
        "## Evidence",
        f"- task FIF files: `{inventory['task_fif_count']}`",
        f"- behavioral CSV segments: `{inventory['behavioral_segment_count']}`",
        f"- reported rows: `{inventory['reported_event_count']}`",
        f"- observed raw windows: `{inventory['observed_raw_window_count']}`",
        f"- resolved observed events: `{inventory['resolved_observed_event_count']}`",
        f"- alignment statuses: `{checks['alignment_status_counts']}`",
        f"- transition warnings: `{checks['transition_warning_count']}`",
        "",
        "## Validation",
        f"- valid: `{validation['is_valid']}`",
    ]
    lines.extend(f"- warning: {item}" for item in validation["warnings"])
    lines.extend(f"- error: {item}" for item in validation["errors"])
    boundaries = report.get("stats", {}).get("acquisition_boundaries", [])
    if boundaries:
        lines.extend(["", "## Reviewed Acquisition Boundaries"])
        for boundary in boundaries:
            lines.append(
                f"- acquisition {boundary['acquisition_run']} -> behavioral blocks "
                f"{boundary['behavioral_blocks']}; alignment "
                f"`{boundary['alignment_status_counts']}`; logged rows beyond "
                f"recording end: {boundary['recording_tail_unavailable_count']}."
            )
    lines.extend(
        [
            "",
            "## Outputs",
            *[f"- {name}: `{path}`" for name, path in outputs.items()],
            "",
        ]
    )
    return "\n".join(lines)


def update_source_manifest(
    manifest_path: Optional[Path],
    associations_payload: Mapping[str, Any],
    *,
    raw_root: Path,
    associations_path: Path,
) -> dict[str, Any]:
    """
    Update a durable source manifest when one exists.

    The manifest is optional provenance, not an input dependency. A standalone
    stage run succeeds and records ``not_updated_missing`` when no manifest was
    supplied or found.
    """
    if manifest_path is None:
        return {"status": "not_requested"}
    if not manifest_path.exists():
        return {
            "status": "not_updated_missing",
            "path": str(manifest_path),
        }
    assert_output_path_outside_raw(manifest_path, raw_root)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Source manifest must contain a JSON object: {manifest_path}")
    updated = deepcopy(payload)
    updated["trigger_derived_associations"] = {
        "status": "resolved",
        "generated_utc": associations_payload["generated_utc"],
        "method": associations_payload["method"],
        "association_artifact": str(associations_path),
        "alignment_status_counts": associations_payload[
            "alignment_status_counts"
        ],
        "meg_acquisition_run_to_behavioral_segments": associations_payload[
            "meg_acquisition_run_to_behavioral_segments"
        ],
    }
    _write_json(manifest_path, updated)
    return {"status": "updated", "path": str(manifest_path)}


# -----------------------------------------------------------------------------
# Command-line orchestration
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the trigger extraction CLI."""
    parser = argparse.ArgumentParser(
        description="Extract and align observed DAD MEG trigger events."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--source-manifest", default=None)
    parser.add_argument("--out-events-csv", default=None)
    parser.add_argument("--out-events-json", default=None)
    parser.add_argument("--out-alignment-csv", default=None)
    parser.add_argument("--out-warnings-csv", default=None)
    parser.add_argument("--out-associations-json", default=None)
    parser.add_argument("--out-report-json", default=None)
    parser.add_argument("--out-report-md", default=None)
    parser.add_argument(
        "--deep-report-out-dir",
        default=None,
        help=(
            "Optional override for the legacy deep report bundle directory. "
            "The default is configured by trigger_parser.deep_report_output_subdir."
        ),
    )
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument(
        "--disable-deep-trigger-report",
        action="store_true",
        help=(
            "Skip generation of the legacy deep Markdown/CSV/plot bundle. "
            "Canonical event and alignment derivatives are still written."
        ),
    )
    return parser


def _output_path(
    value: Optional[str],
    default: Path,
) -> Path:
    """Resolve an optional CLI output override."""
    return Path(value).expanduser().resolve() if value else default.resolve()


def main(argv: Optional[list[str]] = None) -> int:
    """Run the complete read-only-input trigger extraction stage."""
    args = build_arg_parser().parse_args(argv)
    try:
        cfg, _, _, _ = load_and_resolve_config(
            args.config,
            args.subject_config,
            site_config_path=args.site_config,
            required_sections=("trigger_parser",),
        )
        resolved = cfg["resolved_paths"]
        trigger_cfg = cfg["trigger_parser"]
        profile = cfg["subject_profile"]
        subject_id = cfg["subject"]["id"]
        raw_root = Path(resolved["data_root"]).resolve()
        subject_data_dir = Path(resolved["subject_data_dir"]).resolve()
        events_dir = Path(resolved["derivatives_preprocessing_dir"]) / "events"
        reports_dir = Path(resolved["derivatives_reports_dir"])
        deep_report_dir = _output_path(
            args.deep_report_out_dir,
            reports_dir / str(trigger_cfg["deep_report_output_subdir"]),
        )

        paths = {
            "events_csv": _output_path(
                args.out_events_csv,
                events_dir / EVENTS_CSV_NAME,
            ),
            "events_json": _output_path(
                args.out_events_json,
                events_dir / EVENTS_JSON_NAME,
            ),
            "alignment_csv": _output_path(
                args.out_alignment_csv,
                events_dir / ALIGNMENT_CSV_NAME,
            ),
            "transition_warnings_csv": _output_path(
                args.out_warnings_csv,
                events_dir / WARNINGS_CSV_NAME,
            ),
            "associations_json": _output_path(
                args.out_associations_json,
                events_dir / ASSOCIATIONS_JSON_NAME,
            ),
            "report_json": _output_path(
                args.out_report_json,
                reports_dir / f"{STEP_ID}_report.json",
            ),
            "report_md": _output_path(
                args.out_report_md,
                reports_dir / f"{STEP_ID}_report.md",
            ),
        }
        for output_path in paths.values():
            assert_output_path_outside_raw(output_path, raw_root)
        assert_output_path_outside_raw(deep_report_dir, raw_root)

        if args.source_manifest == "":
            manifest_path = None
        elif args.source_manifest:
            manifest_path = Path(args.source_manifest).expanduser().resolve()
        else:
            candidate = reports_dir / "source_manifest.json"
            manifest_path = candidate.resolve() if candidate.exists() else None
        if manifest_path is not None:
            assert_output_path_outside_raw(manifest_path, raw_root)

        fif_files = discover_fif_files(
            subject_data_dir,
            cfg["data_discovery"]["meg_fif_glob"],
        )
        csv_paths = discover_reported_trigger_csvs(
            cfg["data_discovery"]["reported_trigger_search_roots"],
            cfg["data_discovery"]["reported_trigger_csv_glob"],
        )
        segments = load_reported_trigger_segments(
            csv_paths,
            subject_id=subject_id,
            runtime_version=profile["versions"]["runtime"],
            runtime_by_block={
                int(block): str(version)
                for block, version in profile.get("runtime_by_block", {}).items()
            },
            replay_start_code=int(trigger_cfg["replay_start_code"]),
            esc_termination_code=int(trigger_cfg["esc_termination_code"]),
        )
        all_windows: list[ObservedWindow] = []
        fif_summaries: list[dict[str, Any]] = []
        sample_offset = 0
        time_offset = 0.0
        for fif_path in fif_files:
            windows, summary = reconstruct_observed_windows_for_fif(
                fif_path,
                trigger_cfg["stim_channel"],
                float(trigger_cfg.get("min_duration_sec", 0.0)),
                time_offset,
                sample_offset,
            )
            for window in windows:
                window.idx_global = len(all_windows) + 1
                all_windows.append(window)
            fif_summaries.append(summary)
            sample_offset += int(summary["n_samples"])
            time_offset += float(summary["duration_sec"])

        events, transition_warnings, decision_counts, alignment_rows, alignment_summary = (
            process_trigger_windows(
                all_windows, segments, trigger_config=trigger_cfg, subject_profile=profile,
            )
        )
        associations = build_source_associations(events, alignment_summary)

        event_fields = list(events[0].keys()) if events else [
            "event_index",
            "trigger_code",
            "event_label",
            "raw_sample_frame",
            "source_fif",
            "acquisition_run",
            "behavioral_block",
            "behavioral_run",
            "replay_status",
        ]
        csv_events = []
        for event in events:
            row = dict(event)
            for key in (
                "time_start_sec",
                "time_end_sec",
                "time_start_raw_sec",
                "time_end_raw_sec",
            ):
                row[key] = f"{float(row[key]):.6f}"
            csv_events.append(row)

        _write_csv(paths["events_csv"], csv_events, event_fields)
        _write_csv(
            paths["alignment_csv"],
            alignment_rows,
            list(alignment_rows[0].keys()) if alignment_rows else [
                "alignment_index",
                "status",
            ],
        )
        _write_csv(
            paths["transition_warnings_csv"],
            transition_warnings,
            list(transition_warnings[0].keys()) if transition_warnings else [
                "warning_index",
                "source_fif",
                "acquisition_run",
                "order_pattern",
                "dropped_code",
                "kept_code",
                "decision_source",
            ],
        )
        _write_json(
            paths["events_json"],
            {
                "schema_version": "2.0",
                "generated_utc": _iso_utc_now(),
                "step_id": STEP_ID,
                "events": events,
                "alignment": alignment_rows,
                "transition_warnings": transition_warnings,
            },
        )
        _write_json(paths["associations_json"], associations)
        deep_report_enabled = bool(
            trigger_cfg.get("enable_deep_confrontation_report", True)
        ) and not args.disable_deep_trigger_report
        if deep_report_enabled:
            trial_exclusion = cfg.get("epoching", {}).get(
                "trial_exclusion", {}
            )
            deep_report_result = generate_legacy_deep_report(
                subject_id=subject_id,
                raw_root=raw_root,
                output_dir=deep_report_dir,
                stim_channel=str(trigger_cfg["stim_channel"]),
                events=events,
                alignment_rows=alignment_rows,
                transition_warnings=transition_warnings,
                segments=segments,
                segment_summaries=alignment_summary["segments"],
                anomaly_window_ms=float(
                    trigger_cfg.get("deep_report_anomaly_window_ms", 500.0)
                ),
                transition_window_ms=float(
                    trigger_cfg.get("deep_report_transition_window_ms", 3.0)
                ),
                gaze_break_codes=trial_exclusion.get(
                    "gaze_break_codes",
                    [profile["trigger_policy"]["gaze_break_code"]],
                ),
                gaze_break_exclusion_enabled=bool(
                    trial_exclusion.get("exclude_gaze_break", False)
                ),
            )
        else:
            deep_report_result = {
                "status": "disabled",
                "output_dir": str(deep_report_dir),
            }
        manifest_update = update_source_manifest(
            manifest_path,
            associations,
            raw_root=raw_root,
            associations_path=paths["associations_json"],
        )
        report = build_report(
            cfg=cfg,
            fif_summaries=fif_summaries,
            csv_paths=csv_paths,
            segments=segments,
            events=events,
            transition_warnings=transition_warnings,
            decision_counts=decision_counts,
            alignment_summary=alignment_summary,
            manifest_update=manifest_update,
            deep_report_result=deep_report_result,
        )
        _write_json(paths["report_json"], report)
        paths["report_md"].parent.mkdir(parents=True, exist_ok=True)
        paths["report_md"].write_text(
            report_markdown(
                report,
                {
                    **{name: str(path) for name, path in paths.items()},
                    "legacy_deep_report": deep_report_result.get(
                        "report_md",
                        deep_report_result.get("output_dir"),
                    ),
                },
            ),
            encoding="utf-8",
        )

        if args.print_summary:
            print(f"Step {STEP_ID} status: {report['step']['status']}")
            print(f"Runtime: {profile['versions']['runtime']}")
            print(f"Observed events: {len(events)}")
            print(f"Behavioral segments: {len(segments)}")
            print(f"Alignment: {alignment_summary['status_counts']}")
            print(f"Legacy deep report: {deep_report_result['status']}")
            if deep_report_result.get("report_md"):
                print(f"Deep report Markdown: {deep_report_result['report_md']}")
            print(f"Events CSV: {paths['events_csv']}")
        return 0 if report["validation"]["is_valid"] else 2
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Trigger extraction error: {exc}", file=sys.stderr)
        return 2


__all__ = [
    "AlignmentOp",
    "ExpectedTrigger",
    "ObservedWindow",
    "ReportedSegment",
    "align_reported_segments_to_events",
    "build_arg_parser",
    "build_source_associations",
    "discover_fif_files",
    "discover_reported_trigger_csvs",
    "find_events_with_overlaps",
    "flatten_reported_segments",
    "load_expected_trigger_sequence",
    "load_reported_trigger_segments",
    "main",
    "process_trigger_windows",
    "reconstruct_observed_windows_for_fif",
    "resolve_observed_windows_to_events",
    "semiglobal_align_codes",
    "trigger_label",
]
