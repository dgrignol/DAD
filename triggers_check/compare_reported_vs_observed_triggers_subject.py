#!/usr/bin/env python3
"""
Compare reported vs observed triggers across a whole subject in one pass.

Overview
- This subject-level script reads all reported trigger CSV files for one subject,
  ordered by `block` then `run`, concatenates them into one reported sequence,
  reads the corresponding FIF files one at a time (trigger channel only),
  concatenates observed trigger windows, and runs alignment/anomaly checks once
  on the full subject timeline.
- It writes a single report bundle for the subject:
  1) full alignment CSV,
  2) anomaly CSV,
  3) markdown report with block/run start markers.

Important memory behavior
- FIF files are never fully preloaded.
- At any time, only one FIF and only the selected trigger channel are read.
- The script stores only compact trigger-window metadata in memory.

Required CSV columns
- `trigger`, `trial`, `frame`, `seconds`, `label`

Filename expectations
- Reported CSVs: `debug_actual_triggers*_subXX_blockYY_runZZ.csv`
- FIF files: `SUBXX_blockYY.fif` or `SUBXX_blockYY_ZZ.fif`
  (the second form means one FIF contains two blocks).

Usage examples
1) Run subject-level comparison with explicit folders
   python compare_reported_vs_observed_triggers_subject.py \
     --subject 1 \
     --csv-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/output_files/sub01 \
     --fif-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/data/sub01_260415 \
     --out-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/reports/sub01

2) Filter reported CSV candidates by variant token
   python compare_reported_vs_observed_triggers_subject.py \
     --subject 1 \
     --csv-dir /path/to/csv \
     --fif-dir /path/to/fif \
     --out-dir /path/to/out \
     --variant rescueTraject

3) Use non-default trigger channel and transition-collapse threshold
   python compare_reported_vs_observed_triggers_subject.py \
     --subject 1 \
     --csv-dir /path/to/csv \
     --fif-dir /path/to/fif \
     --out-dir /path/to/out \
     --trigger-channel STI101 \
     --transition-max-samples 1

4) Open interactive scroller for manual trigger timeline inspection
   python compare_reported_vs_observed_triggers_subject.py \
     --subject 1 \
     --csv-dir /path/to/csv \
     --fif-dir /path/to/fif \
     --out-dir /path/to/out \
     --open-scroller \
     --scroller-duration 12

Status behavior
- MATCH_WARNING marks reverse transition warnings where a collapsed short pulse
  carried the reported trigger before the kept observed trigger:
  0 -> reported -> random -> 0
  These rows are treated as matched-with-warning (yellow) instead of mismatch.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import mne
import numpy as np


# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------


@dataclass
class FifEntry:
    """One discovered FIF file and the block span encoded in its filename."""

    path: Path
    block_start: int
    block_end: int


@dataclass
class ReportedTrigger:
    """One trigger row from reported debug CSV, extended with block/run context."""

    idx_1based: int
    trigger: int
    trial: Optional[int]
    frame: Optional[int]
    seconds: Optional[float]
    label: str
    block: int
    run: int
    source_csv: str


@dataclass
class RunSegment:
    """One run segment in the concatenated reported stream."""

    block: int
    run: int
    source_csv: Path
    fif_path: Path
    reported_start_idx0: int
    reported_end_idx0_excl: int


@dataclass
class ObservedWindow:
    """One observed trigger window in the concatenated subject timeline."""

    idx_1based: int
    trigger: int
    start_sample: int
    end_sample: int
    start_s: float
    end_s: float
    source_fif: str

    @property
    def width_samples(self) -> int:
        return self.end_sample - self.start_sample

    @property
    def width_ms(self) -> float:
        return 1000.0 * (self.end_s - self.start_s)


@dataclass
class TransitionPulseWarning:
    """One collapsed short transition pulse warning."""

    warning_idx_1based: int
    dropped_trigger: int
    dropped_start_sample: int
    dropped_end_sample: int
    dropped_start_s: float
    dropped_end_s: float
    dropped_width_samples: int
    next_trigger: int
    next_start_sample: int
    next_start_s: float
    source_fif: str
    detail: str
    window_start_s: Optional[float] = None
    window_end_s: Optional[float] = None
    plot_path: str = ""


@dataclass
class AlignmentOp:
    """One alignment operation pairing reported and/or observed entries."""

    op_type: str  # PAIR, REPORTED_ONLY, OBSERVED_ONLY
    reported_idx0: Optional[int]
    observed_idx0: Optional[int]


@dataclass
class Anomaly:
    """One row-level anomaly in the final concatenated alignment."""

    anomaly_idx_1based: int
    anomaly_type: str
    align_row_1based: int
    detail: str
    reported_idx: Optional[int]
    reported_trigger: Optional[int]
    reported_label: str
    reported_block: Optional[int]
    reported_run: Optional[int]
    observed_idx: Optional[int]
    observed_trigger: Optional[int]
    observed_source_fif: str
    observed_width_samples: Optional[int]
    observed_start_s: Optional[float]
    anchor_s: Optional[float] = None
    anchor_source: str = ""
    window_start_s: Optional[float] = None
    window_end_s: Optional[float] = None
    plot_path: str = ""


@dataclass
class ReportedAnchor:
    """One reported trigger mapped to a global-time anchor for scroller overlays."""

    align_row_1based: int
    reported_idx: int
    reported_trigger: int
    reported_label: str
    status: str
    anchor_s: Optional[float]
    anchor_source: str


@dataclass
class RunStartMarker:
    """One run-start marker for the report (block/run + alignment anchor)."""

    block: int
    run: int
    source_csv: str
    source_fif: str
    reported_start_idx: int
    align_row: Optional[int]
    align_status: str
    observed_idx: Optional[int]
    observed_start_s: Optional[float]


@dataclass
class FifTimelineSegment:
    """One FIF segment mapped into the concatenated subject timeline."""

    path: Path
    source_fif: str
    sfreq: float
    first_samp: int
    n_times: int
    local_t0_s: float
    global_time_offset_s: float
    global_sample_offset: int

    @property
    def global_start_s(self) -> float:
        return self.global_time_offset_s

    @property
    def global_end_s(self) -> float:
        return self.global_time_offset_s + (self.n_times / self.sfreq)


# -----------------------------------------------------------------------------
# Parsing helpers
# -----------------------------------------------------------------------------


def _parse_optional_int(value: str | None) -> Optional[int]:
    """Parse optional integer CSV cell."""
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    return int(value)


def _parse_optional_float(value: str | None) -> Optional[float]:
    """Parse optional float CSV cell."""
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    return float(value)


# -----------------------------------------------------------------------------
# Stim-step extraction with overlap support (same core logic as single-run tool)
# -----------------------------------------------------------------------------


def _merge_samples(min_duration: float, sfreq: float) -> int:
    """Convert min duration (seconds) to merge-threshold samples."""
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
    """Detect 1D stim-step changes while preserving sample indices."""
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
        v = data[-1]
        if v != pad_stop:
            last_idx = len(data) + first_samp
            steps = np.append(steps, [[last_idx, v, pad_stop]], axis=0)

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
    """View a 2D numeric array as a 1D structured vector for row-wise isin()."""
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
    Return MNE step events augmented with non-zero -> non-zero overlap steps.

    This preserves overlap transitions explicitly, consistent with prior checks.
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

        non_zero = (steps[:, 1] != 0) & (steps[:, 2] != 0)
        steps = steps[non_zero]
        if steps.size:
            extra_steps.append(steps)

    if not extra_steps:
        return events

    extra_steps = np.vstack(extra_steps)
    if events.size:
        existing = _as_structured(events)
        candidate = _as_structured(extra_steps)
        missing_mask = ~np.isin(candidate, existing)
        extra_steps = extra_steps[missing_mask]

    if extra_steps.size == 0:
        return events

    combined = np.vstack([events, extra_steps]) if events.size else extra_steps.copy()
    combined = np.unique(combined, axis=0)
    combined = combined[np.argsort(combined[:, 0])]
    return combined.astype(events.dtype, copy=False)


# -----------------------------------------------------------------------------
# Discovery and loading
# -----------------------------------------------------------------------------


def discover_reported_csvs(
    csv_dir: Path,
    *,
    subject: int,
    variant: Optional[str],
) -> list[tuple[int, int, Path]]:
    """Discover and sort reported CSV files by block then run."""
    csv_re = re.compile(
        rf"debug_actual_triggers.*_sub{subject:02d}_block(?P<block>\d{{2}})_run(?P<run>\d{{2}})\.csv$",
        re.IGNORECASE,
    )

    candidates = sorted(csv_dir.glob(f"debug_actual_triggers*_sub{subject:02d}_block*_run*.csv"))
    out: list[tuple[int, int, Path]] = []
    for path in candidates:
        name = path.name
        if variant and f"_{variant}_" not in name:
            continue
        match = csv_re.match(name)
        if not match:
            continue
        out.append((int(match.group("block")), int(match.group("run")), path))

    out.sort(key=lambda x: (x[0], x[1], x[2].name))
    if not out:
        raise FileNotFoundError(
            f"No reported CSVs found in {csv_dir} for subject {subject:02d}."
        )
    return out


def discover_fif_entries(fif_dir: Path, *, subject: int) -> list[FifEntry]:
    """Discover FIF files and decode block spans from filenames."""
    fif_re = re.compile(
        rf"sub{subject:02d}_block(?P<b1>\d{{2}})(?:_(?P<b2>\d{{2}}))?\.fif$",
        re.IGNORECASE,
    )

    out: list[FifEntry] = []
    for path in sorted(fif_dir.glob("*.fif")):
        match = fif_re.match(path.name)
        if not match:
            continue
        b1 = int(match.group("b1"))
        b2 = int(match.group("b2")) if match.group("b2") else b1
        if b2 < b1:
            b1, b2 = b2, b1
        out.append(FifEntry(path=path, block_start=b1, block_end=b2))

    out.sort(key=lambda e: (e.block_start, e.block_end, e.path.name))
    if not out:
        raise FileNotFoundError(
            f"No FIF files matching sub{subject:02d}_block*.fif found in {fif_dir}."
        )
    return out


def _select_fif_for_block(block: int, fif_entries: list[FifEntry]) -> FifEntry:
    """Resolve exactly one FIF entry containing the requested block number."""
    matches = [e for e in fif_entries if e.block_start <= block <= e.block_end]
    if not matches:
        raise FileNotFoundError(f"No FIF file covers block {block:02d}.")
    if len(matches) > 1:
        raise ValueError(
            "Ambiguous FIF mapping for block "
            f"{block:02d}: "
            + "; ".join(m.path.name for m in matches)
        )
    return matches[0]


def load_reported_subject_sequence(
    csv_specs: list[tuple[int, int, Path]],
    *,
    fif_entries: list[FifEntry],
) -> tuple[list[ReportedTrigger], list[RunSegment], list[int]]:
    """
    Load and concatenate all reported rows for the subject in block/run order.

    Returns
    - full reported trigger list
    - run segments with source CSV and mapped FIF
    - sorted block list appearing in CSVs
    """
    reported: list[ReportedTrigger] = []
    run_segments: list[RunSegment] = []
    seen_blocks: set[int] = set()

    for block, run, csv_path in csv_specs:
        seen_blocks.add(block)
        fif_entry = _select_fif_for_block(block, fif_entries)

        start_idx0 = len(reported)
        with csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            required = {"trigger", "trial", "frame", "seconds", "label"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"Reported CSV missing columns: {sorted(missing)} in {csv_path}")

            for row in reader:
                reported.append(
                    ReportedTrigger(
                        idx_1based=len(reported) + 1,
                        trigger=int(row["trigger"]),
                        trial=_parse_optional_int(row.get("trial")),
                        frame=_parse_optional_int(row.get("frame")),
                        seconds=_parse_optional_float(row.get("seconds")),
                        label=(row.get("label") or "").strip(),
                        block=block,
                        run=run,
                        source_csv=str(csv_path),
                    )
                )

        end_idx0_excl = len(reported)
        if end_idx0_excl == start_idx0:
            raise ValueError(f"Reported CSV has no rows: {csv_path}")

        run_segments.append(
            RunSegment(
                block=block,
                run=run,
                source_csv=csv_path,
                fif_path=fif_entry.path,
                reported_start_idx0=start_idx0,
                reported_end_idx0_excl=end_idx0_excl,
            )
        )

    return reported, run_segments, sorted(seen_blocks)


# -----------------------------------------------------------------------------
# Observed reconstruction and concatenation
# -----------------------------------------------------------------------------


def reconstruct_observed_windows_one_fif(
    *,
    fif_path: Path,
    trigger_channel: str,
    min_duration: float,
) -> tuple[list[ObservedWindow], float, int, int]:
    """
    Reconstruct observed trigger windows from one FIF file.

    Returns windows with local sample/time coordinates based on FIF indexing,
    plus `(sfreq, first_samp, n_times)` for downstream global concatenation.
    """
    raw = mne.io.read_raw_fif(fif_path, preload=False, verbose="ERROR")
    if trigger_channel not in raw.ch_names:
        raise ValueError(
            f"Trigger channel {trigger_channel!r} not found in {fif_path.name}."
        )

    sfreq = float(raw.info["sfreq"])
    first_samp = int(raw.first_samp)
    n_times = int(raw.n_times)

    events = find_events_with_overlaps(
        raw,
        stim_channel=trigger_channel,
        min_duration=min_duration,
    )

    windows: list[ObservedWindow] = []
    active_trigger: Optional[int] = None
    active_start: Optional[int] = None

    for sample, prev_val, new_val in events:
        sample_i = int(sample)
        prev_i = int(prev_val)
        new_i = int(new_val)

        if prev_i == 0 and new_i > 0:
            if active_trigger is not None and active_start is not None:
                windows.append(
                    ObservedWindow(
                        idx_1based=len(windows) + 1,
                        trigger=active_trigger,
                        start_sample=active_start,
                        end_sample=sample_i,
                        start_s=active_start / sfreq,
                        end_s=sample_i / sfreq,
                        source_fif=str(fif_path),
                    )
                )
            active_trigger = new_i
            active_start = sample_i

        elif prev_i > 0 and new_i == 0:
            if active_trigger == prev_i and active_start is not None:
                windows.append(
                    ObservedWindow(
                        idx_1based=len(windows) + 1,
                        trigger=active_trigger,
                        start_sample=active_start,
                        end_sample=sample_i,
                        start_s=active_start / sfreq,
                        end_s=sample_i / sfreq,
                        source_fif=str(fif_path),
                    )
                )
            active_trigger = None
            active_start = None

        elif prev_i > 0 and new_i > 0:
            if active_trigger == prev_i and active_start is not None:
                windows.append(
                    ObservedWindow(
                        idx_1based=len(windows) + 1,
                        trigger=active_trigger,
                        start_sample=active_start,
                        end_sample=sample_i,
                        start_s=active_start / sfreq,
                        end_s=sample_i / sfreq,
                        source_fif=str(fif_path),
                    )
                )
            active_trigger = new_i
            active_start = sample_i

    return windows, sfreq, first_samp, n_times


def collapse_transition_windows(
    windows: list[ObservedWindow],
    *,
    transient_max_samples: int,
) -> tuple[list[ObservedWindow], list[TransitionPulseWarning]]:
    """Collapse short contiguous transition pulses (0->A->B->0 transients)."""
    if transient_max_samples < 1 or len(windows) < 2:
        return windows, []

    keep = [True] * len(windows)
    warnings: list[TransitionPulseWarning] = []

    for idx in range(len(windows) - 1):
        cur = windows[idx]
        nxt = windows[idx + 1]
        if cur.width_samples <= transient_max_samples and cur.end_sample == nxt.start_sample:
            keep[idx] = False
            warnings.append(
                TransitionPulseWarning(
                    warning_idx_1based=len(warnings) + 1,
                    dropped_trigger=cur.trigger,
                    dropped_start_sample=cur.start_sample,
                    dropped_end_sample=cur.end_sample,
                    dropped_start_s=cur.start_s,
                    dropped_end_s=cur.end_s,
                    dropped_width_samples=cur.width_samples,
                    next_trigger=nxt.trigger,
                    next_start_sample=nxt.start_sample,
                    next_start_s=nxt.start_s,
                    source_fif=cur.source_fif,
                    detail=(
                        "Collapsed short contiguous transition pulse "
                        f"{cur.trigger}->{nxt.trigger} (line-update transient)."
                    ),
                )
            )

    collapsed: list[ObservedWindow] = []
    for old, is_keep in zip(windows, keep):
        if not is_keep:
            continue
        collapsed.append(
            ObservedWindow(
                idx_1based=len(collapsed) + 1,
                trigger=old.trigger,
                start_sample=old.start_sample,
                end_sample=old.end_sample,
                start_s=old.start_s,
                end_s=old.end_s,
                source_fif=old.source_fif,
            )
        )

    return collapsed, warnings


def build_concatenated_observed_sequence(
    *,
    used_fif_entries: list[FifEntry],
    trigger_channel: str,
    min_duration: float,
    transition_max_samples: int,
) -> tuple[list[ObservedWindow], list[TransitionPulseWarning], float, list[FifTimelineSegment]]:
    """
    Build one observed sequence by concatenating FIF files in block-name order.

    Memory note: one FIF is read and processed at a time.
    """
    observed_all: list[ObservedWindow] = []
    warnings_all: list[TransitionPulseWarning] = []
    timeline_segments: list[FifTimelineSegment] = []

    sfreq_ref: Optional[float] = None
    time_offset_s = 0.0
    sample_offset = 0

    for fif_entry in used_fif_entries:
        local_windows, sfreq, first_samp, n_times = reconstruct_observed_windows_one_fif(
            fif_path=fif_entry.path,
            trigger_channel=trigger_channel,
            min_duration=min_duration,
        )
        if sfreq_ref is None:
            sfreq_ref = sfreq
        elif abs(sfreq - sfreq_ref) > 1e-9:
            raise ValueError(
                "All FIF files must have the same sampling frequency for concatenation. "
                f"Got {sfreq_ref} and {sfreq}."
            )

        local_windows, local_warnings = collapse_transition_windows(
            local_windows,
            transient_max_samples=transition_max_samples,
        )

        local_t0_s = first_samp / sfreq
        timeline_segments.append(
            FifTimelineSegment(
                path=fif_entry.path,
                source_fif=fif_entry.path.name,
                sfreq=sfreq,
                first_samp=first_samp,
                n_times=n_times,
                local_t0_s=local_t0_s,
                global_time_offset_s=time_offset_s,
                global_sample_offset=sample_offset,
            )
        )

        for warning in local_warnings:
            warnings_all.append(
                TransitionPulseWarning(
                    warning_idx_1based=len(warnings_all) + 1,
                    dropped_trigger=warning.dropped_trigger,
                    dropped_start_sample=sample_offset + (warning.dropped_start_sample - first_samp),
                    dropped_end_sample=sample_offset + (warning.dropped_end_sample - first_samp),
                    dropped_start_s=time_offset_s + (warning.dropped_start_s - local_t0_s),
                    dropped_end_s=time_offset_s + (warning.dropped_end_s - local_t0_s),
                    dropped_width_samples=warning.dropped_width_samples,
                    next_trigger=warning.next_trigger,
                    next_start_sample=sample_offset + (warning.next_start_sample - first_samp),
                    next_start_s=time_offset_s + (warning.next_start_s - local_t0_s),
                    source_fif=warning.source_fif,
                    detail=warning.detail,
                )
            )

        for win in local_windows:
            observed_all.append(
                ObservedWindow(
                    idx_1based=len(observed_all) + 1,
                    trigger=win.trigger,
                    start_sample=sample_offset + (win.start_sample - first_samp),
                    end_sample=sample_offset + (win.end_sample - first_samp),
                    start_s=time_offset_s + (win.start_s - local_t0_s),
                    end_s=time_offset_s + (win.end_s - local_t0_s),
                    source_fif=win.source_fif,
                )
            )

        sample_offset += n_times
        time_offset_s += n_times / sfreq

    if sfreq_ref is None:
        raise RuntimeError("No FIF files were processed.")

    return observed_all, warnings_all, sfreq_ref, timeline_segments


# -----------------------------------------------------------------------------
# Alignment and anomaly detection
# -----------------------------------------------------------------------------


def semiglobal_align(
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
    *,
    score_match: int = 3,
    score_mismatch: int = -3,
    score_gap: int = -2,
) -> tuple[list[AlignmentOp], int, int]:
    """
    Align full reported sequence to any observed subsequence.

    Returns
    - alignment operations in forward order
    - best score
    - observed index (1-based DP coordinate) where alignment ends
    """
    m = len(reported)
    n = len(observed)

    dp = [[0] * (n + 1) for _ in range(m + 1)]
    bt = [["STOP"] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + score_gap
        bt[i][0] = "U"
    for j in range(1, n + 1):
        dp[0][j] = 0
        bt[0][j] = "L0"

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            is_match = reported[i - 1].trigger == observed[j - 1].trigger
            diag_score = dp[i - 1][j - 1] + (score_match if is_match else score_mismatch)
            up_score = dp[i - 1][j] + score_gap
            left_score = dp[i][j - 1] + score_gap

            best = diag_score
            move = "D"
            if up_score > best:
                best = up_score
                move = "U"
            if left_score > best:
                best = left_score
                move = "L"

            dp[i][j] = best
            bt[i][j] = move

    best_end_j = 0
    best_score = dp[m][0]
    for j in range(1, n + 1):
        if dp[m][j] > best_score:
            best_score = dp[m][j]
            best_end_j = j

    i = m
    j = best_end_j
    rev_ops: list[AlignmentOp] = []
    while i > 0:
        move = bt[i][j]
        if move == "D":
            rev_ops.append(AlignmentOp("PAIR", i - 1, j - 1))
            i -= 1
            j -= 1
        elif move == "U":
            rev_ops.append(AlignmentOp("REPORTED_ONLY", i - 1, None))
            i -= 1
        elif move == "L":
            rev_ops.append(AlignmentOp("OBSERVED_ONLY", None, j - 1))
            j -= 1
        elif move == "L0":
            j -= 1
        else:
            raise RuntimeError(f"Unexpected backtrack move {move!r} at i={i}, j={j}")

    rev_ops.reverse()
    return rev_ops, best_score, best_end_j


def align_subject_by_runs(
    *,
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
    run_segments: list[RunSegment],
) -> tuple[list[AlignmentOp], list[RunStartMarker], int, int]:
    """
    Align concatenated subject sequence run-by-run to keep memory bounded.

    Data flow
    - For each run segment, align only that run's reported rows against the
      remaining observed suffix.
    - Convert per-run op indices to global indices.
    - Advance observed cursor by run alignment end (`end_j`) to preserve order.
    """
    all_ops: list[AlignmentOp] = []
    run_markers: list[RunStartMarker] = []

    obs_cursor = 0
    score_total = 0

    for seg in run_segments:
        rep_start = seg.reported_start_idx0
        rep_end = seg.reported_end_idx0_excl
        rep_run = reported[rep_start:rep_end]
        obs_suffix = observed[obs_cursor:]

        if rep_run and obs_suffix:
            ops_run, score_run, end_j = semiglobal_align(rep_run, obs_suffix)
        elif rep_run and not obs_suffix:
            ops_run = [AlignmentOp("REPORTED_ONLY", i, None) for i in range(len(rep_run))]
            score_run = -2 * len(rep_run)
            end_j = 0
        else:
            ops_run = []
            score_run = 0
            end_j = 0

        score_total += score_run

        run_start_global_rep_idx = rep_run[0].idx_1based if rep_run else None
        run_start_align_row: Optional[int] = None
        run_start_status = "UNMAPPED"
        run_start_obs_idx: Optional[int] = None
        run_start_obs_s: Optional[float] = None

        for op in ops_run:
            rep_idx0 = rep_start + op.reported_idx0 if op.reported_idx0 is not None else None
            obs_idx0 = obs_cursor + op.observed_idx0 if op.observed_idx0 is not None else None
            all_ops.append(AlignmentOp(op.op_type, rep_idx0, obs_idx0))

            if (
                run_start_global_rep_idx is not None
                and rep_idx0 is not None
                and reported[rep_idx0].idx_1based == run_start_global_rep_idx
                and run_start_align_row is None
            ):
                run_start_align_row = len(all_ops)
                if op.op_type == "PAIR" and obs_idx0 is not None:
                    rep = reported[rep_idx0]
                    obs = observed[obs_idx0]
                    run_start_status = "MATCH" if rep.trigger == obs.trigger else "MISMATCH"
                    run_start_obs_idx = obs.idx_1based
                    run_start_obs_s = obs.start_s
                else:
                    run_start_status = op.op_type

        run_markers.append(
            RunStartMarker(
                block=seg.block,
                run=seg.run,
                source_csv=seg.source_csv.name,
                source_fif=seg.fif_path.name,
                reported_start_idx=(run_start_global_rep_idx or -1),
                align_row=run_start_align_row,
                align_status=run_start_status,
                observed_idx=run_start_obs_idx,
                observed_start_s=run_start_obs_s,
            )
        )

        obs_cursor += end_j

    return all_ops, run_markers, score_total, obs_cursor


def _infer_anchor_for_reported_only(
    ops: list[AlignmentOp],
    observed: list[ObservedWindow],
    row_idx0: int,
) -> tuple[Optional[float], str]:
    """Infer global-time anchor for REPORTED_ONLY rows from nearby observed rows."""
    prev_obs: Optional[ObservedWindow] = None
    next_obs: Optional[ObservedWindow] = None

    for j in range(row_idx0 - 1, -1, -1):
        obs_idx0 = ops[j].observed_idx0
        if obs_idx0 is not None:
            prev_obs = observed[obs_idx0]
            break

    for j in range(row_idx0 + 1, len(ops)):
        obs_idx0 = ops[j].observed_idx0
        if obs_idx0 is not None:
            next_obs = observed[obs_idx0]
            break

    if prev_obs is not None and next_obs is not None:
        return ((prev_obs.start_s + next_obs.start_s) / 2.0, "neighbor_midpoint")
    if prev_obs is not None:
        return (prev_obs.start_s, "previous_observed")
    if next_obs is not None:
        return (next_obs.start_s, "next_observed")
    return (None, "unresolved")


def build_reported_anchors(
    *,
    ops: list[AlignmentOp],
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
    row_status_by_align_row: Optional[dict[int, str]] = None,
) -> list[ReportedAnchor]:
    """
    Build per-reported-event global-time anchors for timeline overlays/scroller.

    Why this exists
    - The anomaly list only contains non-match rows, but the scroller needs
      a time anchor for all reported triggers to draw faded vertical reference
      lines consistently across the full timeline.
    """
    anchors: list[ReportedAnchor] = []
    for row_idx0, op in enumerate(ops):
        if op.reported_idx0 is None:
            continue

        row_1based = row_idx0 + 1
        rep = reported[op.reported_idx0]

        if op.observed_idx0 is not None:
            obs = observed[op.observed_idx0]
            anchor_s = obs.start_s
            anchor_source = "observed_start"
        else:
            anchor_s, anchor_source = _infer_anchor_for_reported_only(
                ops, observed, row_idx0
            )

        status = op.op_type
        if row_status_by_align_row is not None:
            status = row_status_by_align_row.get(row_1based, status)

        anchors.append(
            ReportedAnchor(
                align_row_1based=row_1based,
                reported_idx=rep.idx_1based,
                reported_trigger=rep.trigger,
                reported_label=rep.label,
                status=status,
                anchor_s=anchor_s,
                anchor_source=anchor_source,
            )
        )

    return anchors


def detect_anomalies(
    *,
    ops: list[AlignmentOp],
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
    row_status_by_align_row: Optional[dict[int, str]] = None,
    include_width_outliers: bool,
) -> list[Anomaly]:
    """
    Build anomaly rows from final global alignment operations.

    Notes
    - Rows marked as MATCH_WARNING in the alignment table are intentionally
      excluded from MISMATCH anomalies. This covers both transition orders
      where a short line-update transient is collapsed:
      0->reported->unexpected->0 and 0->unexpected->reported->0.
    """
    anomalies: list[Anomaly] = []

    for row_idx0, op in enumerate(ops):
        row_1based = row_idx0 + 1
        rep = reported[op.reported_idx0] if op.reported_idx0 is not None else None
        obs = observed[op.observed_idx0] if op.observed_idx0 is not None else None

        if op.op_type == "PAIR":
            assert rep is not None and obs is not None
            if rep.trigger != obs.trigger:
                if (
                    row_status_by_align_row is not None
                    and row_status_by_align_row.get(row_1based) == "MATCH_WARNING"
                ):
                    continue
                anomalies.append(
                    Anomaly(
                        anomaly_idx_1based=len(anomalies) + 1,
                        anomaly_type="MISMATCH",
                        align_row_1based=row_1based,
                        detail=(
                            f"Reported trigger {rep.trigger} differs from observed "
                            f"trigger {obs.trigger}."
                        ),
                        reported_idx=rep.idx_1based,
                        reported_trigger=rep.trigger,
                        reported_label=rep.label,
                        reported_block=rep.block,
                        reported_run=rep.run,
                        observed_idx=obs.idx_1based,
                        observed_trigger=obs.trigger,
                        observed_source_fif=Path(obs.source_fif).name,
                        observed_width_samples=obs.width_samples,
                        observed_start_s=obs.start_s,
                    )
                )

        elif op.op_type == "REPORTED_ONLY":
            assert rep is not None
            anchor_s, anchor_source = _infer_anchor_for_reported_only(ops, observed, row_idx0)
            anomalies.append(
                Anomaly(
                    anomaly_idx_1based=len(anomalies) + 1,
                    anomaly_type="REPORTED_ONLY",
                    align_row_1based=row_1based,
                    detail="Reported trigger has no aligned observed pulse.",
                    reported_idx=rep.idx_1based,
                    reported_trigger=rep.trigger,
                    reported_label=rep.label,
                    reported_block=rep.block,
                    reported_run=rep.run,
                    observed_idx=None,
                    observed_trigger=None,
                    observed_source_fif="",
                    observed_width_samples=None,
                    observed_start_s=None,
                    anchor_s=anchor_s,
                    anchor_source=anchor_source,
                )
            )

        elif op.op_type == "OBSERVED_ONLY":
            assert obs is not None
            anomalies.append(
                Anomaly(
                    anomaly_idx_1based=len(anomalies) + 1,
                    anomaly_type="OBSERVED_ONLY",
                    align_row_1based=row_1based,
                    detail="Observed trigger has no aligned reported entry.",
                    reported_idx=None,
                    reported_trigger=None,
                    reported_label="",
                    reported_block=None,
                    reported_run=None,
                    observed_idx=obs.idx_1based,
                    observed_trigger=obs.trigger,
                    observed_source_fif=Path(obs.source_fif).name,
                    observed_width_samples=obs.width_samples,
                    observed_start_s=obs.start_s,
                    anchor_s=obs.start_s,
                    anchor_source="observed_start",
                )
            )

    if include_width_outliers:
        aligned_widths = [
            observed[op.observed_idx0].width_samples
            for op in ops
            if op.observed_idx0 is not None
        ]
        if aligned_widths:
            values, counts = np.unique(np.asarray(aligned_widths), return_counts=True)
            modal_width = int(values[np.argmax(counts)])
            for row_idx0, op in enumerate(ops):
                if op.observed_idx0 is None:
                    continue
                obs = observed[op.observed_idx0]
                if obs.width_samples == modal_width:
                    continue
                rep = reported[op.reported_idx0] if op.reported_idx0 is not None else None
                anomalies.append(
                    Anomaly(
                        anomaly_idx_1based=len(anomalies) + 1,
                        anomaly_type="WIDTH_OUTLIER",
                        align_row_1based=row_idx0 + 1,
                        detail=(
                            f"Observed width is {obs.width_samples} samples; "
                            f"modal width is {modal_width} samples."
                        ),
                        reported_idx=(rep.idx_1based if rep else None),
                        reported_trigger=(rep.trigger if rep else None),
                        reported_label=(rep.label if rep else ""),
                        reported_block=(rep.block if rep else None),
                        reported_run=(rep.run if rep else None),
                        observed_idx=obs.idx_1based,
                        observed_trigger=obs.trigger,
                        observed_source_fif=Path(obs.source_fif).name,
                        observed_width_samples=obs.width_samples,
                        observed_start_s=obs.start_s,
                        anchor_s=obs.start_s,
                        anchor_source="observed_start",
                    )
                )

    # Set anchors for mismatch rows after construction (all PAIR rows with mismatch).
    for anomaly in anomalies:
        if anomaly.anchor_s is not None:
            continue
        if anomaly.anomaly_type == "MISMATCH" and anomaly.observed_start_s is not None:
            anomaly.anchor_s = anomaly.observed_start_s
            anomaly.anchor_source = "observed_start"

    return anomalies


# -----------------------------------------------------------------------------
# Plotting and report formatting helpers
# -----------------------------------------------------------------------------


def _safe_slug(text: str) -> str:
    """Return a filesystem-safe lowercase slug."""
    out = []
    for ch in text.lower():
        if ch.isalnum():
            out.append(ch)
        elif ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_")


def _md_escape(text: str) -> str:
    """Escape markdown table separators in free-text cells."""
    return text.replace("|", "\\|")


def _path_for_report(report_path: Path, target: str) -> str:
    """Format absolute path relative to report directory when possible."""
    if not target:
        return ""
    target_path = Path(target).expanduser()
    if not target_path.is_absolute():
        return target_path.as_posix()
    try:
        return Path(os.path.relpath(target_path, start=report_path.parent)).as_posix()
    except Exception:
        return target_path.as_posix()


def _status_bg_color(status: str) -> str:
    """Map status string to requested box background color."""
    st = status.upper()
    if "WARNING" in st:
        return "#f5e7a1"  # yellow
    if st == "MATCH":
        return "#b6e3b6"  # green
    if st in {"MISMATCH", "REPORTED_ONLY", "OBSERVED_ONLY"}:
        return "#f3b0b0"  # red
    return "#dddddd"


def _status_badge(status: str) -> str:
    """Return HTML badge span for markdown tables."""
    color = _status_bg_color(status)
    label = _md_escape(status)
    return (
        f'<span style="background:{color};padding:2px 6px;'
        f'border-radius:4px;border:1px solid #999;">{label}</span>'
    )


def _find_segment_by_source(
    timeline_segments: list[FifTimelineSegment],
    source_fif_name: str,
) -> Optional[FifTimelineSegment]:
    """Find timeline segment by FIF basename."""
    for seg in timeline_segments:
        if seg.source_fif == source_fif_name:
            return seg
    return None


def _find_segment_by_global_time(
    timeline_segments: list[FifTimelineSegment],
    t_global_s: float,
) -> Optional[FifTimelineSegment]:
    """Find timeline segment containing the given global time."""
    for seg in timeline_segments:
        if seg.global_start_s <= t_global_s <= seg.global_end_s:
            return seg
    return None


def _extract_window_from_segment(
    *,
    segment: FifTimelineSegment,
    trigger_channel: str,
    window_start_s: float,
    window_end_s: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Extract trigger-channel samples for a global-time window from one FIF segment.

    Returns
    - times in global concatenated seconds
    - trigger data
    - actual extracted start/end in global seconds
    """
    raw = mne.io.read_raw_fif(segment.path, preload=False, verbose="ERROR")
    sfreq = segment.sfreq

    local_start_s = segment.local_t0_s + (window_start_s - segment.global_time_offset_s)
    local_end_s = segment.local_t0_s + (window_end_s - segment.global_time_offset_s)

    abs_min = segment.first_samp
    abs_max = segment.first_samp + segment.n_times - 1
    start_abs = max(abs_min, int(round(local_start_s * sfreq)))
    end_abs = min(abs_max, int(round(local_end_s * sfreq)))
    if end_abs <= start_abs:
        raise RuntimeError("Window extraction collapsed to empty interval.")

    rel_start = start_abs - segment.first_samp
    rel_stop = end_abs - segment.first_samp + 1
    trig_data = raw.get_data(
        picks=[trigger_channel],
        start=rel_start,
        stop=rel_stop,
        reject_by_annotation="omit",
    )[0]
    local_times = (np.arange(rel_start, rel_stop) + segment.first_samp) / sfreq
    global_times = segment.global_time_offset_s + (local_times - segment.local_t0_s)
    return global_times, trig_data, float(global_times[0]), float(global_times[-1])


def _extract_window_across_segments(
    *,
    timeline_segments: list[FifTimelineSegment],
    trigger_channel: str,
    window_start_s: float,
    window_end_s: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Extract trigger samples for one global window across one or more FIF segments.

    This keeps the memory model bounded: only the overlapping segment slices are
    read, and no full-subject trigger array is materialized.
    """
    if window_end_s <= window_start_s:
        raise ValueError("window_end_s must be greater than window_start_s.")

    times_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []

    for segment in timeline_segments:
        if segment.global_end_s < window_start_s or segment.global_start_s > window_end_s:
            continue

        seg_start = max(window_start_s, segment.global_start_s)
        seg_end = min(window_end_s, segment.global_end_s)
        if seg_end <= seg_start:
            continue

        times_s, trig_data, _actual_start, _actual_end = _extract_window_from_segment(
            segment=segment,
            trigger_channel=trigger_channel,
            window_start_s=seg_start,
            window_end_s=seg_end,
        )

        # Remove boundary duplicates caused by sample rounding at segment joins.
        if times_parts:
            prev_last = float(times_parts[-1][-1])
            keep = times_s > prev_last
            if not np.any(keep):
                continue
            times_s = times_s[keep]
            trig_data = trig_data[keep]

        times_parts.append(times_s)
        data_parts.append(trig_data)

    if not times_parts:
        raise RuntimeError("No data extracted for the requested global window.")

    all_times = np.concatenate(times_parts)
    all_data = np.concatenate(data_parts)
    return all_times, all_data, float(all_times[0]), float(all_times[-1])


def open_trigger_scroller(
    *,
    timeline_segments: list[FifTimelineSegment],
    trigger_channel: str,
    reported_anchors: list[ReportedAnchor],
    duration_s: float,
    start_s: Optional[float],
    step_s: Optional[float],
) -> None:
    """
    Open an interactive trigger-channel scroller over concatenated FIF timeline.

    Controls
    - Slider: move visible-window start time.
    - Left/Right or A/D: move backward/forward by step_s.
    - Home/End: jump to timeline start/end.
    """
    backend = str(plt.get_backend()).lower()
    if "agg" in backend or "inline" in backend:
        print(
            "Scroller not opened: matplotlib backend is non-interactive "
            f"({plt.get_backend()})."
        )
        return

    if not timeline_segments:
        raise RuntimeError("Cannot open scroller: no timeline segments available.")
    if duration_s <= 0:
        raise ValueError(f"scroller duration must be > 0, got {duration_s}")

    min_t = float(timeline_segments[0].global_start_s)
    max_t = float(timeline_segments[-1].global_end_s)
    full_span = max_t - min_t
    if full_span <= 0:
        raise RuntimeError("Cannot open scroller: timeline has non-positive duration.")
    if duration_s > full_span:
        duration_s = full_span
    if step_s is None or step_s <= 0:
        step_s = duration_s / 2.0

    finite_anchors = np.asarray(
        [a.anchor_s for a in reported_anchors if a.anchor_s is not None and np.isfinite(a.anchor_s)],
        dtype=float,
    )

    # Default start follows first anchor when available.
    if start_s is None:
        if finite_anchors.size:
            start_s = float(finite_anchors[0] - duration_s * 0.25)
        else:
            start_s = min_t

    def _clamp_start(v: float) -> float:
        return max(min_t, min(float(v), max_t - duration_s))

    start_s = _clamp_start(start_s)

    fig, ax = plt.subplots(figsize=(12, 4.5))
    plt.subplots_adjust(bottom=0.23)
    (line,) = ax.step([], [], where="post", color="#1f77b4", linewidth=0.9)
    ax.set_xlabel("Global Time (s)")
    ax.set_ylabel(trigger_channel)
    ax.set_title(
        f"Trigger scroller ({trigger_channel}) | "
        "Faded lines = reported-trigger anchors"
    )
    ax.grid(True, alpha=0.25)

    anchor_artists: list = []

    def _draw_window(left: float) -> None:
        nonlocal anchor_artists
        left = _clamp_start(left)
        right = left + duration_s

        try:
            times_s, trig_data, actual_start, actual_end = _extract_window_across_segments(
                timeline_segments=timeline_segments,
                trigger_channel=trigger_channel,
                window_start_s=left,
                window_end_s=right,
            )
        except Exception:
            return

        line.set_data(times_s, trig_data)
        ax.set_xlim(actual_start, actual_end)

        if trig_data.size:
            y_min = float(np.nanmin(trig_data))
            y_max = float(np.nanmax(trig_data))
            pad = max((y_max - y_min) * 0.1, 1.0)
            ax.set_ylim(y_min - pad, y_max + pad)

        for artist in anchor_artists:
            artist.remove()
        anchor_artists = []

        if finite_anchors.size:
            in_window = finite_anchors[
                (finite_anchors >= actual_start) & (finite_anchors <= actual_end)
            ]
            for t in in_window:
                anchor_artists.append(
                    ax.axvline(t, color="#555555", alpha=0.16, linewidth=0.8)
                )

        fig.canvas.draw_idle()

    slider_ax = fig.add_axes([0.11, 0.08, 0.78, 0.04])
    slider = Slider(
        slider_ax,
        "Window start (s)",
        min_t,
        max_t - duration_s,
        valinit=start_s,
    )

    def _on_slider(val: float) -> None:
        _draw_window(float(val))

    def _on_key(event) -> None:
        left = float(slider.val)
        if event.key in ("right", "d"):
            slider.set_val(_clamp_start(left + step_s))
        elif event.key in ("left", "a"):
            slider.set_val(_clamp_start(left - step_s))
        elif event.key == "home":
            slider.set_val(min_t)
        elif event.key == "end":
            slider.set_val(max_t - duration_s)

    slider.on_changed(_on_slider)
    fig.canvas.mpl_connect("key_press_event", _on_key)
    _draw_window(start_s)
    plt.show()


def plot_anomaly_windows(
    *,
    anomalies: list[Anomaly],
    observed: list[ObservedWindow],
    timeline_segments: list[FifTimelineSegment],
    trigger_channel: str,
    window_ms: float,
    out_dir: Path,
) -> None:
    """Save per-anomaly trigger-channel windows around anomaly anchors."""
    out_dir.mkdir(parents=True, exist_ok=True)
    half_window_s = window_ms / 1000.0

    for anomaly in anomalies:
        if anomaly.anchor_s is None or not np.isfinite(anomaly.anchor_s):
            continue

        seg: Optional[FifTimelineSegment] = None
        if anomaly.observed_source_fif:
            seg = _find_segment_by_source(timeline_segments, anomaly.observed_source_fif)
        if seg is None:
            seg = _find_segment_by_global_time(timeline_segments, anomaly.anchor_s)
        if seg is None:
            continue

        req_start = anomaly.anchor_s - half_window_s
        req_end = anomaly.anchor_s + half_window_s
        req_start = max(req_start, seg.global_start_s)
        req_end = min(req_end, seg.global_end_s)
        if req_end <= req_start:
            continue

        try:
            times_s, trig_data, actual_start, actual_end = _extract_window_from_segment(
                segment=seg,
                trigger_channel=trigger_channel,
                window_start_s=req_start,
                window_end_s=req_end,
            )
        except Exception:
            continue

        fig, ax = plt.subplots(figsize=(10, 3.4))
        ax.step(times_s, trig_data, where="post", color="#1f77b4", linewidth=1.0)
        ax.axvline(anomaly.anchor_s, color="#d62728", linestyle="--", linewidth=1.2)

        for obs in observed:
            if Path(obs.source_fif).name != seg.source_fif:
                continue
            if actual_start <= obs.start_s <= actual_end:
                ax.axvline(obs.start_s, color="#bbbbbb", linewidth=0.7, alpha=0.8)

        status_color = _status_bg_color(anomaly.anomaly_type)
        ax.text(
            0.01,
            0.98,
            f"STATUS: {anomaly.anomaly_type}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": status_color,
                "edgecolor": "#777777",
                "alpha": 0.95,
            },
        )
        ax.text(
            0.01,
            0.86,
            f"Source FIF: {seg.source_fif}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=7,
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "edgecolor": "#cccccc",
                "alpha": 0.85,
            },
        )

        ax.set_xlabel("Global Time (s)")
        ax.set_ylabel(trigger_channel)
        ax.set_title(
            f"{anomaly.anomaly_type} | row {anomaly.align_row_1based} | "
            f"rep block/run: {anomaly.reported_block}/{anomaly.reported_run}"
        )
        ax.grid(True, alpha=0.25)
        ax.set_xlim(actual_start, actual_end)
        fig.tight_layout()

        filename = (
            f"{anomaly.anomaly_idx_1based:04d}_"
            f"{_safe_slug(anomaly.anomaly_type)}_"
            f"row{anomaly.align_row_1based:05d}.png"
        )
        plot_path = out_dir / filename
        fig.savefig(plot_path, dpi=150)
        plt.close(fig)

        anomaly.window_start_s = actual_start
        anomaly.window_end_s = actual_end
        anomaly.plot_path = str(plot_path)


def plot_transition_warning_windows(
    *,
    warnings: list[TransitionPulseWarning],
    timeline_segments: list[FifTimelineSegment],
    trigger_channel: str,
    window_ms: float,
    out_dir: Path,
) -> None:
    """Save transition-warning trigger windows with yellow status box."""
    if not warnings:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    half_window_s = window_ms / 1000.0

    for warning in warnings:
        seg = _find_segment_by_source(timeline_segments, Path(warning.source_fif).name)
        if seg is None:
            seg = _find_segment_by_global_time(timeline_segments, warning.dropped_start_s)
        if seg is None:
            continue

        req_start = max(seg.global_start_s, warning.dropped_start_s - half_window_s)
        req_end = min(seg.global_end_s, warning.dropped_start_s + half_window_s)
        if req_end <= req_start:
            continue

        try:
            times_s, trig_data, actual_start, actual_end = _extract_window_from_segment(
                segment=seg,
                trigger_channel=trigger_channel,
                window_start_s=req_start,
                window_end_s=req_end,
            )
        except Exception:
            continue

        fig, ax = plt.subplots(figsize=(10, 3.4))
        ax.step(times_s, trig_data, where="post", color="#1f77b4", linewidth=1.0)
        ax.axvline(warning.dropped_start_s, color="#ff7f0e", linestyle="--", linewidth=1.3)
        ax.axvline(warning.next_start_s, color="#d62728", linestyle=":", linewidth=1.1)
        ax.axvspan(warning.dropped_start_s, warning.dropped_end_s, color="#ff7f0e", alpha=0.2)

        ax.text(
            0.01,
            0.98,
            "STATUS: WARNING",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": _status_bg_color("WARNING"),
                "edgecolor": "#777777",
                "alpha": 0.95,
            },
        )
        ax.text(
            0.01,
            0.86,
            "Orange region: dropped transient pulse",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=7,
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "edgecolor": "#cccccc",
                "alpha": 0.85,
            },
        )

        ax.set_xlabel("Global Time (s)")
        ax.set_ylabel(trigger_channel)
        ax.set_title(
            f"TRANSITION_WARNING {warning.warning_idx_1based} | "
            f"dropped={warning.dropped_trigger} next={warning.next_trigger}"
        )
        ax.grid(True, alpha=0.25)
        ax.set_xlim(actual_start, actual_end)
        fig.tight_layout()

        filename = (
            f"{warning.warning_idx_1based:04d}_transition_warning_"
            f"{warning.dropped_trigger}_to_{warning.next_trigger}.png"
        )
        plot_path = out_dir / filename
        fig.savefig(plot_path, dpi=150)
        plt.close(fig)

        warning.window_start_s = actual_start
        warning.window_end_s = actual_end
        warning.plot_path = str(plot_path)


# -----------------------------------------------------------------------------
# Output writers
# -----------------------------------------------------------------------------


def _fmt_optional_int(value: Optional[int]) -> str:
    return "" if value is None else str(value)


def _fmt_optional_float(value: Optional[float], digits: int = 6) -> str:
    return "" if value is None else f"{value:.{digits}f}"


def build_alignment_rows(
    *,
    ops: list[AlignmentOp],
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
    run_markers: list[RunStartMarker],
    transition_warnings: list[TransitionPulseWarning],
) -> tuple[list[dict[str, str]], dict[str, int], dict[int, list[str]]]:
    """
    Build full row-wise output with run/block marker labels.

    Transition-warning handling
    - If a PAIR row corresponds to a collapsed short transition pulse and the
      reported trigger participates in that transition, mark status as
      MATCH_WARNING.
    - This covers both orders:
      0 -> reported -> unexpected -> 0
      0 -> unexpected -> reported -> 0
    """
    stats = {
        "MATCH": 0,
        "MATCH_WARNING": 0,
        "MISMATCH": 0,
        "REPORTED_ONLY": 0,
        "OBSERVED_ONLY": 0,
    }
    rows: list[dict[str, str]] = []

    # Index warnings by kept-window start sample to support row-level status
    # reclassification for reverse-transition patterns.
    warnings_by_next_start: dict[int, list[TransitionPulseWarning]] = {}
    for warning in transition_warnings:
        warnings_by_next_start.setdefault(warning.next_start_sample, []).append(warning)

    run_start_repidx = {m.reported_start_idx for m in run_markers if m.reported_start_idx > 0}
    block_start_repidx = set()
    first_run_by_block: dict[int, int] = {}
    for m in run_markers:
        if m.reported_start_idx <= 0:
            continue
        if m.block not in first_run_by_block or m.run < first_run_by_block[m.block]:
            first_run_by_block[m.block] = m.run
            block_start_repidx.add(m.reported_start_idx)

    markers_by_row: dict[int, list[str]] = {}

    for out_idx, op in enumerate(ops, start=1):
        rep = reported[op.reported_idx0] if op.reported_idx0 is not None else None
        obs = observed[op.observed_idx0] if op.observed_idx0 is not None else None

        if op.op_type == "PAIR":
            assert rep is not None and obs is not None
            if rep.trigger == obs.trigger:
                status = "MATCH"
            else:
                status = "MISMATCH"
            # Reclassify to MATCH_WARNING whenever this aligned observed pulse
            # is the kept side of a collapsed transition and the reported code
            # appears in either transition leg.
            for warning in warnings_by_next_start.get(obs.start_sample, []):
                if warning.next_trigger != obs.trigger:
                    continue
                if rep.trigger in {warning.dropped_trigger, warning.next_trigger}:
                    status = "MATCH_WARNING"
                    break
        elif op.op_type == "REPORTED_ONLY":
            status = "REPORTED_ONLY"
        else:
            status = "OBSERVED_ONLY"

        stats[status] += 1

        markers: list[str] = []
        if rep is not None:
            if rep.idx_1based in block_start_repidx:
                markers.append("BLOCK_START")
            if rep.idx_1based in run_start_repidx:
                markers.append("RUN_START")
        markers_by_row[out_idx] = markers

        rows.append(
            {
                "align_row": str(out_idx),
                "status": status,
                "segment_marker": ",".join(markers),
                "reported_idx": _fmt_optional_int(rep.idx_1based if rep else None),
                "reported_block": _fmt_optional_int(rep.block if rep else None),
                "reported_run": _fmt_optional_int(rep.run if rep else None),
                "reported_trigger": _fmt_optional_int(rep.trigger if rep else None),
                "reported_trial": _fmt_optional_int(rep.trial if rep else None),
                "reported_frame": _fmt_optional_int(rep.frame if rep else None),
                "reported_seconds": _fmt_optional_float(rep.seconds if rep else None, digits=6),
                "reported_label": rep.label if rep else "",
                "reported_source_csv": Path(rep.source_csv).name if rep else "",
                "observed_idx": _fmt_optional_int(obs.idx_1based if obs else None),
                "observed_trigger": _fmt_optional_int(obs.trigger if obs else None),
                "observed_start_s": _fmt_optional_float(obs.start_s if obs else None, digits=6),
                "observed_end_s": _fmt_optional_float(obs.end_s if obs else None, digits=6),
                "observed_width_ms": _fmt_optional_float(obs.width_ms if obs else None, digits=3),
                "observed_source_fif": Path(obs.source_fif).name if obs else "",
            }
        )

    return rows, stats, markers_by_row


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write generic rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("No rows to write.")
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_anomalies_csv(path: Path, anomalies: list[Anomaly]) -> None:
    """Write anomaly CSV for machine inspection."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "anomaly_idx",
        "anomaly_type",
        "align_row",
        "detail",
        "reported_idx",
        "reported_block",
        "reported_run",
        "reported_trigger",
        "reported_label",
        "observed_idx",
        "observed_trigger",
        "observed_source_fif",
        "observed_width_samples",
        "observed_start_s",
        "anchor_s",
        "anchor_source",
        "window_start_s",
        "window_end_s",
        "plot_path",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for a in anomalies:
            writer.writerow(
                {
                    "anomaly_idx": a.anomaly_idx_1based,
                    "anomaly_type": a.anomaly_type,
                    "align_row": a.align_row_1based,
                    "detail": a.detail,
                    "reported_idx": _fmt_optional_int(a.reported_idx),
                    "reported_block": _fmt_optional_int(a.reported_block),
                    "reported_run": _fmt_optional_int(a.reported_run),
                    "reported_trigger": _fmt_optional_int(a.reported_trigger),
                    "reported_label": a.reported_label,
                    "observed_idx": _fmt_optional_int(a.observed_idx),
                    "observed_trigger": _fmt_optional_int(a.observed_trigger),
                    "observed_source_fif": a.observed_source_fif,
                    "observed_width_samples": _fmt_optional_int(a.observed_width_samples),
                    "observed_start_s": _fmt_optional_float(a.observed_start_s, digits=6),
                    "anchor_s": _fmt_optional_float(a.anchor_s, digits=6),
                    "anchor_source": a.anchor_source,
                    "window_start_s": _fmt_optional_float(a.window_start_s, digits=6),
                    "window_end_s": _fmt_optional_float(a.window_end_s, digits=6),
                    "plot_path": a.plot_path,
                }
            )


def write_markdown_report(
    path: Path,
    *,
    subject: int,
    csv_dir: Path,
    fif_dir: Path,
    out_csv: Path,
    out_anomalies_csv: Path,
    anomaly_plots_dir: Path,
    transition_plots_dir: Path,
    transition_plots_disabled: bool,
    n_csv_files: int,
    n_fif_files: int,
    run_markers: list[RunStartMarker],
    alignment_rows: list[dict[str, str]],
    stats: dict[str, int],
    anomalies: list[Anomaly],
    transition_warnings: list[TransitionPulseWarning],
    score_total: int,
    observed_cursor_end: int,
    anomaly_plot_window_ms: float,
    transition_plot_window_ms: float,
) -> None:
    """Write single subject-level markdown report with start-marker tables."""
    path.parent.mkdir(parents=True, exist_ok=True)

    n_reported = len([r for r in alignment_rows if r["reported_idx"]])
    n_observed = len([r for r in alignment_rows if r["observed_idx"]])
    effective_match = stats.get("MATCH", 0) + stats.get("MATCH_WARNING", 0)
    match_rate = (effective_match / n_reported * 100.0) if n_reported else 0.0

    block_starts: list[RunStartMarker] = []
    seen_blocks: set[int] = set()
    for marker in run_markers:
        if marker.block in seen_blocks:
            continue
        seen_blocks.add(marker.block)
        block_starts.append(marker)

    csv_dir_disp = _path_for_report(path, str(csv_dir))
    fif_dir_disp = _path_for_report(path, str(fif_dir))
    out_csv_disp = _path_for_report(path, str(out_csv))
    out_anom_disp = _path_for_report(path, str(out_anomalies_csv))
    anomaly_plots_disp = _path_for_report(path, str(anomaly_plots_dir))
    transition_plots_disp = _path_for_report(path, str(transition_plots_dir))

    lines: list[str] = []
    lines.append("# Subject Trigger Confrontation Report")
    lines.append("")
    lines.append("## Inputs")
    lines.append(f"- Subject: **sub{subject:02d}**")
    lines.append(f"- Reported CSV directory: `{csv_dir_disp}`")
    lines.append(f"- FIF directory: `{fif_dir_disp}`")
    lines.append(f"- Discovered reported CSV files: **{n_csv_files}**")
    lines.append(f"- Discovered FIF files used: **{n_fif_files}**")
    lines.append("")

    lines.append("## Summary")
    lines.append(f"- Total reported events: **{n_reported}**")
    lines.append(f"- Total observed windows (inside aligned rows): **{n_observed}**")
    lines.append(f"- Aggregate alignment score (sum across runs): **{score_total}**")
    lines.append(f"- MATCH: **{stats['MATCH']}**")
    lines.append(f"- MATCH_WARNING: **{stats.get('MATCH_WARNING', 0)}**")
    lines.append(f"- MISMATCH: **{stats['MISMATCH']}**")
    lines.append(f"- REPORTED_ONLY: **{stats['REPORTED_ONLY']}**")
    lines.append(f"- OBSERVED_ONLY: **{stats['OBSERVED_ONLY']}**")
    lines.append(
        f"- Reported match rate (MATCH + MATCH_WARNING): **{match_rate:.1f}%**"
    )
    lines.append(f"- Transition warnings (collapsed pulses): **{len(transition_warnings)}**")
    lines.append(f"- Observed cursor end index after final run: **{observed_cursor_end}**")
    lines.append("")

    lines.append("## Artifacts")
    lines.append(f"- Alignment CSV: `{out_csv_disp}`")
    lines.append(f"- Anomaly CSV: `{out_anom_disp}`")
    lines.append(f"- Anomaly plots directory: `{anomaly_plots_disp}`")
    lines.append(f"- Transition warning plots directory: `{transition_plots_disp}`")
    lines.append(f"- This report: `{_path_for_report(path, str(path))}`")
    lines.append("")

    lines.append("## Block Starts")
    if not block_starts:
        lines.append("No block-start markers available.")
    else:
        lines.append("| block | first_run | source_csv | source_fif | reported_start_idx | align_row | status | observed_idx | observed_start_s |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for marker in block_starts:
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(marker.block),
                        str(marker.run),
                        marker.source_csv,
                        marker.source_fif,
                        str(marker.reported_start_idx),
                        _fmt_optional_int(marker.align_row),
                        _status_badge(marker.align_status),
                        _fmt_optional_int(marker.observed_idx),
                        _fmt_optional_float(marker.observed_start_s, digits=6),
                    ]
                )
                + " |"
            )
    lines.append("")

    lines.append("## Run Starts")
    if not run_markers:
        lines.append("No run-start markers available.")
    else:
        lines.append("| block | run | source_csv | source_fif | reported_start_idx | align_row | status | observed_idx | observed_start_s |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for marker in run_markers:
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(marker.block),
                        str(marker.run),
                        marker.source_csv,
                        marker.source_fif,
                        str(marker.reported_start_idx),
                        _fmt_optional_int(marker.align_row),
                        _status_badge(marker.align_status),
                        _fmt_optional_int(marker.observed_idx),
                        _fmt_optional_float(marker.observed_start_s, digits=6),
                    ]
                )
                + " |"
            )
    lines.append("")

    lines.append("## Anomaly Plots")
    lines.append(f"- Plot window around each anomaly: ±{anomaly_plot_window_ms:.1f} ms.")
    if not anomalies:
        lines.append("No anomalies detected.")
    else:
        lines.append("")
        lines.append(
            "| anomaly_idx | status | align_row | reported_block | reported_run | "
            "reported_trigger | observed_trigger | anchor_s | anchor_source | plot_path |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for a in anomalies:
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(a.anomaly_idx_1based),
                        _md_escape(a.anomaly_type),
                        str(a.align_row_1based),
                        _fmt_optional_int(a.reported_block),
                        _fmt_optional_int(a.reported_run),
                        _fmt_optional_int(a.reported_trigger),
                        _fmt_optional_int(a.observed_trigger),
                        _fmt_optional_float(a.anchor_s, digits=6),
                        _md_escape(a.anchor_source),
                        _md_escape(_path_for_report(path, a.plot_path)),
                    ]
                )
                + " |"
            )
        lines.append("")
        for a in anomalies:
            if not a.plot_path:
                continue
            lines.append(
                f"### Anomaly {a.anomaly_idx_1based}: {a.anomaly_type} "
                f"(row {a.align_row_1based})"
            )
            lines.append(f"- Status: {_status_badge(a.anomaly_type)}")
            lines.append(f"- Detail: {a.detail}")
            lines.append(
                f"![Anomaly {a.anomaly_idx_1based}]"
                f"({_path_for_report(path, a.plot_path)})"
            )
            lines.append("")

    lines.append("## Non-MATCH Rows")
    non_match_rows = [r for r in alignment_rows if r["status"] != "MATCH"]
    if not non_match_rows:
        lines.append("All aligned rows are MATCH.")
    else:
        lines.append(
            "| align_row | segment_marker | status | reported_block | reported_run | "
            "reported_trigger | reported_label | observed_trigger | observed_start_s | observed_source_fif |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in non_match_rows:
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["align_row"],
                        row["segment_marker"],
                        _status_badge(row["status"]),
                        row["reported_block"],
                        row["reported_run"],
                        row["reported_trigger"],
                        row["reported_label"].replace("|", "\\|"),
                        row["observed_trigger"],
                        row["observed_start_s"],
                        row["observed_source_fif"],
                    ]
                )
                + " |"
            )
    lines.append("")

    lines.append("## Transition Warnings")
    lines.append(f"- Plot window around each warning: ±{transition_plot_window_ms:.1f} ms.")
    if transition_plots_disabled:
        lines.append("- Transition warning plots were disabled.")
    lines.append(
        "- Yellow status/background indicates warning context. "
        "Orange region in plots marks the collapsed transient pulse."
    )
    if not transition_warnings:
        lines.append("No transition warnings.")
    else:
        lines.append("")
        lines.append(
            "| idx | status | dropped | next | dropped_start_s | width_samples | "
            "source_fif | detail | plot_path |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for w in transition_warnings:
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(w.warning_idx_1based),
                        "WARNING",
                        str(w.dropped_trigger),
                        str(w.next_trigger),
                        f"{w.dropped_start_s:.6f}",
                        str(w.dropped_width_samples),
                        Path(w.source_fif).name,
                        w.detail.replace("|", "\\|"),
                        _md_escape(_path_for_report(path, w.plot_path)),
                    ]
                )
                + " |"
            )
        lines.append("")
        if not transition_plots_disabled:
            for w in transition_warnings:
                if not w.plot_path:
                    continue
                lines.append(
                    f"### Transition Warning {w.warning_idx_1based}: "
                    f"{w.dropped_trigger}->{w.next_trigger}"
                )
                lines.append(f"- Status: {_status_badge('WARNING')}")
                lines.append(
                    f"![Transition Warning {w.warning_idx_1based}]"
                    f"({_path_for_report(path, w.plot_path)})"
                )
                lines.append("")
    lines.append("")

    lines.append("## Full Alignment Table")
    lines.append(
        "- Status colors: green = MATCH, red = MISMATCH/REPORTED_ONLY/OBSERVED_ONLY, "
        "yellow = warning statuses (e.g., MATCH_WARNING)."
    )
    lines.append("")
    lines.append(
        "| align_row | segment_marker | status | reported_idx | reported_block | "
        "reported_run | reported_trigger | reported_trial | reported_frame | "
        "reported_seconds | reported_label | reported_source_csv | observed_idx | "
        "observed_trigger | observed_start_s | observed_end_s | observed_width_ms | observed_source_fif |"
    )
    lines.append(
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
    )
    for row in alignment_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["align_row"],
                    row["segment_marker"],
                    _status_badge(row["status"]),
                    row["reported_idx"],
                    row["reported_block"],
                    row["reported_run"],
                    row["reported_trigger"],
                    row["reported_trial"],
                    row["reported_frame"],
                    row["reported_seconds"],
                    row["reported_label"].replace("|", "\\|"),
                    row["reported_source_csv"],
                    row["observed_idx"],
                    row["observed_trigger"],
                    row["observed_start_s"],
                    row["observed_end_s"],
                    row["observed_width_ms"],
                    row["observed_source_fif"],
                ]
            )
            + " |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# -----------------------------------------------------------------------------
# CLI and main
# -----------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser for subject-level concatenated comparison."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare all reported CSVs and FIF trigger windows for one subject "
            "as one concatenated sequence."
        )
    )
    parser.add_argument("--subject", type=int, required=True, help="Subject number, e.g. 1.")
    parser.add_argument(
        "--csv-dir",
        type=Path,
        required=True,
        help="Directory containing reported CSV files for the subject.",
    )
    parser.add_argument(
        "--fif-dir",
        type=Path,
        required=True,
        help="Directory containing subject FIF files (blockYY or blockYY_ZZ naming).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory for single subject-level report artifacts.",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="Optional variant token filter for CSV names (e.g., rescueTraject).",
    )
    parser.add_argument(
        "--trigger-channel",
        type=str,
        default="STI101",
        help="Trigger channel to read from FIF files (default: STI101).",
    )
    parser.add_argument(
        "--min-duration",
        type=float,
        default=0.0,
        help="Minimum event duration in seconds for step merge logic (default: 0).",
    )
    parser.add_argument(
        "--transition-max-samples",
        type=int,
        default=1,
        help="Max width (samples) treated as transition transient (default: 1).",
    )
    parser.add_argument(
        "--include-width-outliers",
        action="store_true",
        help="Include WIDTH_OUTLIER anomalies (disabled by default).",
    )
    parser.add_argument(
        "--plot-window-ms",
        type=float,
        default=500.0,
        help="Half-window size in ms around anomaly anchors (default: 500).",
    )
    parser.add_argument(
        "--transition-plot-window-ms",
        type=float,
        default=3.0,
        help="Half-window size in ms around transition warnings (default: 3).",
    )
    parser.add_argument(
        "--disable-plot-transition-pulses",
        action="store_true",
        help="Disable transition-warning plot generation (warnings still reported).",
    )
    parser.add_argument(
        "--open-scroller",
        action="store_true",
        help=(
            "Open an interactive trigger-channel scroller in a separate window "
            "after outputs are generated."
        ),
    )
    parser.add_argument(
        "--scroller-duration",
        type=float,
        default=10.0,
        help="Visible window duration in seconds for the scroller (default: 10).",
    )
    parser.add_argument(
        "--scroller-start",
        type=float,
        default=None,
        help=(
            "Initial scroller window start time in seconds. "
            "Defaults to the first reported-trigger anchor."
        ),
    )
    parser.add_argument(
        "--scroller-step",
        type=float,
        default=None,
        help=(
            "Keyboard navigation step in seconds (left/right, a/d). "
            "Defaults to half the scroller duration."
        ),
    )
    return parser


def main() -> None:
    """Run subject-level concatenated comparison and write one report bundle."""
    args = build_parser().parse_args()

    subject = int(args.subject)
    csv_dir = args.csv_dir.expanduser().resolve()
    fif_dir = args.fif_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_specs = discover_reported_csvs(csv_dir, subject=subject, variant=args.variant)
    fif_entries = discover_fif_entries(fif_dir, subject=subject)

    reported, run_segments, blocks = load_reported_subject_sequence(
        csv_specs,
        fif_entries=fif_entries,
    )

    used_fif_by_name: dict[str, FifEntry] = {}
    for seg in run_segments:
        used_fif_by_name[seg.fif_path.name] = _select_fif_for_block(seg.block, fif_entries)
    used_fif_entries = sorted(
        used_fif_by_name.values(),
        key=lambda e: (e.block_start, e.block_end, e.path.name),
    )

    observed, transition_warnings, sfreq, timeline_segments = build_concatenated_observed_sequence(
        used_fif_entries=used_fif_entries,
        trigger_channel=args.trigger_channel,
        min_duration=args.min_duration,
        transition_max_samples=args.transition_max_samples,
    )

    if not observed:
        raise RuntimeError("No observed windows reconstructed from subject FIF files.")

    ops, run_markers, score_total, obs_cursor_end = align_subject_by_runs(
        reported=reported,
        observed=observed,
        run_segments=run_segments,
    )

    alignment_rows, stats, _markers = build_alignment_rows(
        ops=ops,
        reported=reported,
        observed=observed,
        run_markers=run_markers,
        transition_warnings=transition_warnings,
    )

    row_status_by_align_row = {
        int(row["align_row"]): row["status"] for row in alignment_rows if row["align_row"]
    }

    reported_anchors: list[ReportedAnchor] = []
    if args.open_scroller:
        reported_anchors = build_reported_anchors(
            ops=ops,
            reported=reported,
            observed=observed,
            row_status_by_align_row=row_status_by_align_row,
        )

    anomalies = detect_anomalies(
        ops=ops,
        reported=reported,
        observed=observed,
        row_status_by_align_row=row_status_by_align_row,
        include_width_outliers=args.include_width_outliers,
    )

    stem = f"sub{subject:02d}_all_blocks_reported_vs_observed"
    out_csv = out_dir / f"{stem}.csv"
    out_anom_csv = out_dir / f"{stem}_anomalies.csv"
    out_report = out_dir / f"{stem}.md"
    anomaly_plots_dir = out_dir / f"{stem}_plots"
    transition_plots_dir = anomaly_plots_dir / "transition_pulses"

    plot_anomaly_windows(
        anomalies=anomalies,
        observed=observed,
        timeline_segments=timeline_segments,
        trigger_channel=args.trigger_channel,
        window_ms=args.plot_window_ms,
        out_dir=anomaly_plots_dir,
    )
    if not args.disable_plot_transition_pulses:
        plot_transition_warning_windows(
            warnings=transition_warnings,
            timeline_segments=timeline_segments,
            trigger_channel=args.trigger_channel,
            window_ms=args.transition_plot_window_ms,
            out_dir=transition_plots_dir,
        )

    write_csv(out_csv, alignment_rows)
    write_anomalies_csv(out_anom_csv, anomalies)
    write_markdown_report(
        out_report,
        subject=subject,
        csv_dir=csv_dir,
        fif_dir=fif_dir,
        out_csv=out_csv,
        out_anomalies_csv=out_anom_csv,
        anomaly_plots_dir=anomaly_plots_dir,
        transition_plots_dir=transition_plots_dir,
        transition_plots_disabled=args.disable_plot_transition_pulses,
        n_csv_files=len(csv_specs),
        n_fif_files=len(used_fif_entries),
        run_markers=run_markers,
        alignment_rows=alignment_rows,
        stats=stats,
        anomalies=anomalies,
        transition_warnings=transition_warnings,
        score_total=score_total,
        observed_cursor_end=obs_cursor_end,
        anomaly_plot_window_ms=args.plot_window_ms,
        transition_plot_window_ms=args.transition_plot_window_ms,
    )

    print(f"Subject: sub{subject:02d}")
    print(f"Discovered CSV files: {len(csv_specs)}")
    print(f"Discovered block numbers: {', '.join(f'{b:02d}' for b in blocks)}")
    print(f"FIF files used: {len(used_fif_entries)}")
    print(f"Sampling frequency: {sfreq:.3f} Hz")
    print(f"Reported events (concatenated): {len(reported)}")
    print(f"Observed windows (concatenated): {len(observed)}")
    print(f"Aggregate alignment score: {score_total}")
    print(
        "Status counts: "
        f"MATCH={stats['MATCH']} "
        f"MATCH_WARNING={stats.get('MATCH_WARNING', 0)} "
        f"MISMATCH={stats['MISMATCH']} "
        f"REPORTED_ONLY={stats['REPORTED_ONLY']} "
        f"OBSERVED_ONLY={stats['OBSERVED_ONLY']}"
    )
    print(f"Transition warnings: {len(transition_warnings)}")
    print(f"Anomalies: {len(anomalies)}")
    print(f"Saved alignment CSV: {out_csv}")
    print(f"Saved anomaly CSV: {out_anom_csv}")
    print(f"Saved anomaly plots: {anomaly_plots_dir}")
    if args.disable_plot_transition_pulses:
        print("Transition warning plots disabled.")
    else:
        print(f"Saved transition warning plots: {transition_plots_dir}")
    print(f"Saved report: {out_report}")

    if args.open_scroller:
        print(
            "Opening interactive scroller. "
            "Use slider or left/right (a/d) keys to navigate."
        )
        open_trigger_scroller(
            timeline_segments=timeline_segments,
            trigger_channel=args.trigger_channel,
            reported_anchors=reported_anchors,
            duration_s=args.scroller_duration,
            start_s=args.scroller_start,
            step_s=args.scroller_step,
        )


if __name__ == "__main__":
    main()
