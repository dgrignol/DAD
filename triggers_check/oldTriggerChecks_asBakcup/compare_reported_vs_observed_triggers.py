#!/usr/bin/env python3
"""
Compare reported trigger logs against observed trigger windows extracted from a FIF.

This script is designed for trigger-confrontation workflows where a MATLAB debug
CSV ("reported" triggers) must be checked against the trigger channel stored in
an MEG FIF file ("observed" triggers). The script builds explicit trigger windows
(start/end samples and times) from the FIF stim channel, aligns the reported and
observed trigger sequences, and writes:

1) A row-wise comparison CSV with alignment status per event.
2) An anomaly CSV (MISMATCH / REPORTED_ONLY / OBSERVED_ONLY, plus optional
   width outliers) with plotting anchors.
3) A Markdown report summarizing metrics and embedding anomaly plots.
4) One trigger-channel PNG per anomaly centered on a local time window.
5) Optional interactive trigger-channel scroller with faded vertical lines at
   reported-trigger anchors.

Key behavior
- FIF extraction uses MNE step events and includes non-zero -> non-zero overlaps.
- Observed triggers are represented as windows with start/end sample and time.
- Sequence alignment is semi-global:
  - full reported sequence must be aligned;
  - observed sequence can have free unmatched prefix/suffix;
  - insertions/deletions are preserved as REPORTED_ONLY / OBSERVED_ONLY rows.

Inputs
- reported CSV columns expected: trigger, trial, frame, label
- FIF file with a stim channel (default: STI101)

Outputs
- comparison CSV with columns including:
  status, reported_trigger, observed_trigger, observed_start_s, observed_end_s
- anomaly CSV with anomaly metadata + plot paths
- Markdown report with summary metrics, anomaly table, and full alignment table
- anomaly plot PNGs (one per anomaly, default window +/-500 ms)

Usage examples
1) Basic comparison
   python compare_reported_vs_observed_triggers.py \
     --reported-csv /path/debug_actual_triggers_occlusion_v16_sub51_block01_run01.csv \
     --fif /path/Subj51_4blocks_fakeGazeBreakReplay.fif \
     --out-csv /tmp/sub51_block01_run01_confrontation.csv \
     --out-report /tmp/sub51_block01_run01_confrontation_report.md

2) Explicit stim channel
   python compare_reported_vs_observed_triggers.py \
     --reported-csv /path/debug.csv \
     --fif /path/data.fif \
     --trigger-channel STI101 \
     --out-csv /tmp/confrontation.csv \
     --out-report /tmp/confrontation.md

3) Custom anomaly window (+/-750 ms) and explicit plot directory
   python compare_reported_vs_observed_triggers.py \
     --reported-csv /path/debug.csv \
     --fif /path/data.fif \
     --out-csv /tmp/confrontation.csv \
     --out-report /tmp/confrontation.md \
     --plot-window-ms 750 \
     --plots-dir /tmp/confrontation_plots

4) Open interactive scroller after writing outputs
   python compare_reported_vs_observed_triggers.py \
     --reported-csv /path/debug.csv \
     --fif /path/data.fif \
     --out-csv /tmp/confrontation.csv \
     --out-report /tmp/confrontation.md \
     --open-scroller \
     --scroller-duration 10

Assumptions and constraints
- Reported CSV rows are ordered by emission time from the experiment script.
- Observed windows are ordered by FIF sample time.
- Alignment is sequence-based; no absolute time synchronization with trial/frame
  metadata is required for this first-pass confrontation.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
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
class ReportedTrigger:
    """One trigger row from the debug CSV produced by the experiment."""

    idx_1based: int
    trigger: int
    trial: Optional[int]
    frame: Optional[int]
    label: str


@dataclass
class ObservedWindow:
    """One observed trigger window reconstructed from the FIF stim channel."""

    idx_1based: int
    trigger: int
    start_sample: int
    end_sample: int
    start_s: float
    end_s: float

    @property
    def width_samples(self) -> int:
        return self.end_sample - self.start_sample

    @property
    def width_ms(self) -> float:
        return 1000.0 * (self.end_s - self.start_s)


@dataclass
class AlignmentOp:
    """One alignment operation pairing reported and/or observed entries."""

    op_type: str  # PAIR, REPORTED_ONLY, OBSERVED_ONLY
    reported_idx0: Optional[int]
    observed_idx0: Optional[int]


@dataclass
class Anomaly:
    """One anomaly used for reporting and trigger-window plotting."""

    anomaly_idx_1based: int
    anomaly_type: str
    align_row_1based: int
    detail: str
    reported_idx: Optional[int]
    reported_trigger: Optional[int]
    reported_label: str
    observed_idx: Optional[int]
    observed_trigger: Optional[int]
    observed_width_samples: Optional[int]
    anchor_s: float
    anchor_source: str
    window_start_s: Optional[float] = None
    window_end_s: Optional[float] = None
    plot_path: str = ""


@dataclass
class ReportedAnchor:
    """
    One reported-trigger anchor mapped into FIF time for plotting.

    For rows that align to an observed pulse, anchor_s is the observed start time.
    For REPORTED_ONLY rows, anchor_s is inferred from neighboring observed rows.
    """

    align_row_1based: int
    reported_idx: int
    reported_trigger: int
    reported_label: str
    status: str
    anchor_s: float
    anchor_source: str


# -----------------------------------------------------------------------------
# Helpers for robust stim-step extraction with overlap support
# -----------------------------------------------------------------------------


def _merge_samples(min_duration: float, sfreq: float) -> int:
    """Convert min duration in seconds to an MNE-compatible merge threshold."""
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
    """Detect 1D stim step changes while preserving sample indices."""
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
    """Apply MNE-like bitmask filtering to step values."""
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
    Return MNE step events, augmented with non-zero -> non-zero overlap steps.

    This mirrors the behavior used in older DAD trigger-inspection scripts,
    allowing overlap transitions to be represented explicitly.
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

        # Keep only transitions where both pre and post values are non-zero.
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
# Input loading
# -----------------------------------------------------------------------------


def _parse_optional_int(value: str | None) -> Optional[int]:
    """Parse an integer-like cell, treating empty values as missing."""
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    return int(value)


def load_reported_csv(path: Path) -> list[ReportedTrigger]:
    """Load reported trigger rows from the debug CSV file."""
    out: list[ReportedTrigger] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"trigger", "trial", "frame", "label"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Reported CSV missing columns: {sorted(missing)}")

        for idx, row in enumerate(reader, start=1):
            out.append(
                ReportedTrigger(
                    idx_1based=idx,
                    trigger=int(row["trigger"]),
                    trial=_parse_optional_int(row.get("trial")),
                    frame=_parse_optional_int(row.get("frame")),
                    label=(row.get("label") or "").strip(),
                )
            )
    if not out:
        raise ValueError(f"Reported CSV has no rows: {path}")
    return out


def reconstruct_observed_windows(
    fif_path: Path,
    trigger_channel: str,
    min_duration: float,
) -> tuple[list[ObservedWindow], float, int]:
    """
    Extract observed trigger windows from FIF step events.

    Data flow
    1) Read FIF metadata + stim channel.
    2) Detect step transitions with overlap support.
    3) Convert transitions into explicit windows (start/end per trigger pulse).
    """
    raw = mne.io.read_raw_fif(fif_path, preload=False, verbose="ERROR")
    if trigger_channel not in raw.ch_names:
        raise ValueError(
            f"Trigger channel {trigger_channel!r} not found in FIF. "
            f"Available channels include: {', '.join(raw.ch_names[-10:])}"
        )

    sfreq = float(raw.info["sfreq"])
    events = find_events_with_overlaps(
        raw,
        stim_channel=trigger_channel,
        min_duration=min_duration,
    )

    windows: list[ObservedWindow] = []
    active_trigger: Optional[int] = None
    active_start: Optional[int] = None

    # Convert step events into trigger windows while preserving overlap boundaries.
    for sample, prev_val, new_val in events:
        sample_i = int(sample)
        prev_i = int(prev_val)
        new_i = int(new_val)

        if prev_i == 0 and new_i > 0:
            # Standard pulse onset.
            if active_trigger is not None and active_start is not None:
                # Defensive close in case an onset arrives before prior close.
                windows.append(
                    ObservedWindow(
                        idx_1based=len(windows) + 1,
                        trigger=active_trigger,
                        start_sample=active_start,
                        end_sample=sample_i,
                        start_s=active_start / sfreq,
                        end_s=sample_i / sfreq,
                    )
                )
            active_trigger = new_i
            active_start = sample_i

        elif prev_i > 0 and new_i == 0:
            # Standard pulse offset.
            if active_trigger == prev_i and active_start is not None:
                windows.append(
                    ObservedWindow(
                        idx_1based=len(windows) + 1,
                        trigger=active_trigger,
                        start_sample=active_start,
                        end_sample=sample_i,
                        start_s=active_start / sfreq,
                        end_s=sample_i / sfreq,
                    )
                )
            active_trigger = None
            active_start = None

        elif prev_i > 0 and new_i > 0:
            # Overlap transition: close previous at this sample and open new.
            if active_trigger == prev_i and active_start is not None:
                windows.append(
                    ObservedWindow(
                        idx_1based=len(windows) + 1,
                        trigger=active_trigger,
                        start_sample=active_start,
                        end_sample=sample_i,
                        start_s=active_start / sfreq,
                        end_s=sample_i / sfreq,
                    )
                )
            active_trigger = new_i
            active_start = sample_i

    return windows, sfreq, raw.first_samp


# -----------------------------------------------------------------------------
# Sequence alignment
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

    Alignment model
    - Reported sequence must be fully consumed.
    - Observed prefix/suffix can be skipped without penalty.
    - Internal gaps are penalized.

    Returns
    - alignment operations in forward order
    - best alignment score
    - observed index (1-based in DP table) where alignment ends
    """
    m = len(reported)
    n = len(observed)

    dp = [[0] * (n + 1) for _ in range(m + 1)]
    bt = [["STOP"] * (n + 1) for _ in range(m + 1)]

    # Initialization: full reported must be aligned; observed prefix is free.
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + score_gap
        bt[i][0] = "U"  # reported-only
    for j in range(1, n + 1):
        dp[0][j] = 0
        bt[0][j] = "L0"  # free observed-prefix skip

    # Fill DP matrix with deterministic tie-breaking: D > U > L.
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

    # Free observed suffix: choose best terminal column in last reported row.
    best_end_j = 0
    best_score = dp[m][0]
    for j in range(1, n + 1):
        if dp[m][j] > best_score:
            best_score = dp[m][j]
            best_end_j = j

    # Backtrack until the full reported sequence is consumed.
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
            # Free prefix skip; move left.
            j -= 1
        else:
            raise RuntimeError(f"Unexpected backtrack move {move!r} at i={i}, j={j}")

    rev_ops.reverse()
    return rev_ops, best_score, best_end_j


# -----------------------------------------------------------------------------
# Anomaly detection and plotting
# -----------------------------------------------------------------------------


def _infer_anchor_for_reported_only(
    ops: list[AlignmentOp],
    observed: list[ObservedWindow],
    row_idx0: int,
) -> tuple[float, str]:
    """
    Infer a plotting anchor for REPORTED_ONLY rows from neighboring observed rows.

    REPORTED_ONLY rows do not have a direct observed timestamp in the FIF stream.
    For visualization, anchor the missing event to the midpoint between the
    closest observed rows on the left/right when possible.
    """
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
    return (float("nan"), "unresolved")


def build_reported_anchors(
    ops: list[AlignmentOp],
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
) -> list[ReportedAnchor]:
    """
    Build per-reported-event anchors in FIF time for plotting overlays.

    Data flow
    - MATCH / MISMATCH rows: anchor to observed window start.
    - REPORTED_ONLY rows: infer anchor from neighbor observed windows.
    """
    anchors: list[ReportedAnchor] = []
    for row_idx0, op in enumerate(ops):
        if op.reported_idx0 is None:
            continue
        rep = reported[op.reported_idx0]
        if op.observed_idx0 is not None:
            obs = observed[op.observed_idx0]
            anchor_s = obs.start_s
            anchor_source = "observed_start"
        else:
            anchor_s, anchor_source = _infer_anchor_for_reported_only(
                ops, observed, row_idx0
            )
        anchors.append(
            ReportedAnchor(
                align_row_1based=row_idx0 + 1,
                reported_idx=rep.idx_1based,
                reported_trigger=rep.trigger,
                reported_label=rep.label,
                status=op.op_type,
                anchor_s=anchor_s,
                anchor_source=anchor_source,
            )
        )
    return anchors


def detect_anomalies(
    ops: list[AlignmentOp],
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
    *,
    include_width_outliers: bool,
) -> list[Anomaly]:
    """
    Detect alignment anomalies and optional width outliers.

    Reported anomalies:
    - MISMATCH
    - REPORTED_ONLY
    - OBSERVED_ONLY
    Optional:
    - WIDTH_OUTLIER (observed pulse width differs from modal width)
    """
    anomalies: list[Anomaly] = []

    # First pass: alignment-based anomalies.
    for row_idx0, op in enumerate(ops):
        row_1based = row_idx0 + 1
        rep = reported[op.reported_idx0] if op.reported_idx0 is not None else None
        obs = observed[op.observed_idx0] if op.observed_idx0 is not None else None

        if op.op_type == "PAIR":
            assert rep is not None and obs is not None
            if rep.trigger != obs.trigger:
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
                        observed_idx=obs.idx_1based,
                        observed_trigger=obs.trigger,
                        observed_width_samples=obs.width_samples,
                        anchor_s=obs.start_s,
                        anchor_source="observed_start",
                    )
                )

        elif op.op_type == "REPORTED_ONLY":
            assert rep is not None
            anchor_s, anchor_source = _infer_anchor_for_reported_only(
                ops, observed, row_idx0
            )
            anomalies.append(
                Anomaly(
                    anomaly_idx_1based=len(anomalies) + 1,
                    anomaly_type="REPORTED_ONLY",
                    align_row_1based=row_1based,
                    detail="Reported trigger has no aligned observed pulse in FIF.",
                    reported_idx=rep.idx_1based,
                    reported_trigger=rep.trigger,
                    reported_label=rep.label,
                    observed_idx=None,
                    observed_trigger=None,
                    observed_width_samples=None,
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
                    detail="Observed trigger has no aligned reported entry in debug CSV.",
                    reported_idx=None,
                    reported_trigger=None,
                    reported_label="",
                    observed_idx=obs.idx_1based,
                    observed_trigger=obs.trigger,
                    observed_width_samples=obs.width_samples,
                    anchor_s=obs.start_s,
                    anchor_source="observed_start",
                )
            )

    # Second pass: pulse-width outliers among aligned observed rows.
    if include_width_outliers:
        aligned_observed_widths = [
            observed[op.observed_idx0].width_samples
            for op in ops
            if op.observed_idx0 is not None
        ]
        if aligned_observed_widths:
            modal_width = Counter(aligned_observed_widths).most_common(1)[0][0]
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
                        reported_idx=(rep.idx_1based if rep is not None else None),
                        reported_trigger=(rep.trigger if rep is not None else None),
                        reported_label=(rep.label if rep is not None else ""),
                        observed_idx=obs.idx_1based,
                        observed_trigger=obs.trigger,
                        observed_width_samples=obs.width_samples,
                        anchor_s=obs.start_s,
                        anchor_source="observed_start",
                    )
                )

    return anomalies


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


def open_trigger_scroller(
    *,
    raw: mne.io.BaseRaw,
    trigger_channel: str,
    reported_anchors: list[ReportedAnchor],
    duration_s: float,
    start_s: Optional[float],
    step_s: Optional[float],
) -> None:
    """
    Open an interactive trigger-channel scroller with reported-trigger markers.

    Controls
    - Slider: move the visible window start time.
    - Left/Right or A/D: move backward/forward by step_s.
    - Home/End: jump to recording start/end.
    """
    backend = str(plt.get_backend()).lower()
    if "agg" in backend or "inline" in backend:
        print(
            "Scroller not opened: matplotlib backend is non-interactive "
            f"({plt.get_backend()})."
        )
        return

    if duration_s <= 0:
        raise ValueError(f"scroller duration must be > 0, got {duration_s}")

    sfreq = float(raw.info["sfreq"])
    signal = raw.get_data(picks=[trigger_channel], reject_by_annotation="omit")[0]
    times = (np.arange(raw.n_times) + raw.first_samp) / sfreq
    min_t = float(times[0])
    max_t = float(times[-1])
    full_span = max_t - min_t
    if full_span <= 0:
        raise RuntimeError("Cannot open scroller: recording has non-positive duration.")
    if duration_s > full_span:
        duration_s = full_span
    if step_s is None or step_s <= 0:
        step_s = duration_s / 2.0

    # Choose a useful default start based on the first finite reported anchor.
    if start_s is None:
        finite_anchors = [
            a.anchor_s for a in reported_anchors if np.isfinite(a.anchor_s)
        ]
        if finite_anchors:
            start_s = finite_anchors[0] - duration_s * 0.25
        else:
            start_s = min_t
    start_s = max(min_t, min(start_s, max_t - duration_s))

    fig, ax = plt.subplots(figsize=(12, 4.5))
    plt.subplots_adjust(bottom=0.23)

    # Plot full signal once; scrolling is done by x-limits for responsiveness.
    ax.step(times, signal, where="post", color="#1f77b4", linewidth=0.8)

    # Faded vertical lines mark where triggers are reported (requested behavior).
    for anchor in reported_anchors:
        if not np.isfinite(anchor.anchor_s):
            continue
        ax.axvline(anchor.anchor_s, color="#555555", alpha=0.16, linewidth=0.8)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(trigger_channel)
    ax.set_title(
        f"Trigger scroller ({trigger_channel}) | "
        "Faded lines = reported-trigger anchors"
    )
    ax.grid(True, alpha=0.25)
    ax.set_xlim(start_s, start_s + duration_s)

    # Slider controls the left edge of the visible time window.
    slider_ax = fig.add_axes([0.11, 0.08, 0.78, 0.04])
    slider = Slider(
        slider_ax,
        "Window start (s)",
        min_t,
        max_t - duration_s,
        valinit=start_s,
    )

    def _clamp_start(v: float) -> float:
        return max(min_t, min(float(v), max_t - duration_s))

    def _update_from_slider(val: float) -> None:
        left = _clamp_start(val)
        ax.set_xlim(left, left + duration_s)
        fig.canvas.draw_idle()

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

    slider.on_changed(_update_from_slider)
    fig.canvas.mpl_connect("key_press_event", _on_key)
    plt.show()


def plot_anomaly_windows(
    *,
    raw: mne.io.BaseRaw,
    trigger_channel: str,
    observed: list[ObservedWindow],
    reported_anchors: list[ReportedAnchor],
    anomalies: list[Anomaly],
    window_ms: float,
    out_dir: Path,
) -> None:
    """
    Save one trigger-channel plot per anomaly around its anchor (±window_ms).

    The plot includes:
    - raw trigger channel signal (step plot),
    - anchor line for the anomaly,
    - nearby observed trigger onsets with trigger-code labels.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    sfreq = float(raw.info["sfreq"])
    win_samples = int(round((window_ms / 1000.0) * sfreq))
    if win_samples < 1:
        raise ValueError(f"window_ms too small: {window_ms}")

    abs_min_sample = raw.first_samp
    abs_max_sample = raw.first_samp + raw.n_times - 1

    for anom in anomalies:
        if not np.isfinite(anom.anchor_s):
            continue

        anchor_sample = int(round(anom.anchor_s * sfreq))
        start_sample = max(abs_min_sample, anchor_sample - win_samples)
        end_sample = min(abs_max_sample, anchor_sample + win_samples)
        if end_sample <= start_sample:
            continue

        rel_start = start_sample - raw.first_samp
        rel_stop = end_sample - raw.first_samp + 1  # inclusive -> exclusive
        trig_data = raw.get_data(
            picks=[trigger_channel],
            start=rel_start,
            stop=rel_stop,
            reject_by_annotation="omit",
        )[0]
        times_s = (np.arange(rel_start, rel_stop) + raw.first_samp) / sfreq

        fig, ax = plt.subplots(figsize=(10, 3.2))
        ax.step(times_s, trig_data, where="post", linewidth=1.0, color="#1f77b4")
        ax.axvline(anom.anchor_s, color="#d62728", linestyle="--", linewidth=1.2)

        # Add faded lines for reported anchors in this local window.
        for rep_anchor in reported_anchors:
            if not np.isfinite(rep_anchor.anchor_s):
                continue
            if (start_sample / sfreq) <= rep_anchor.anchor_s <= (end_sample / sfreq):
                ax.axvline(rep_anchor.anchor_s, color="#555555", alpha=0.16, linewidth=0.8)

        # Annotate nearby observed trigger onsets for context in this window.
        for obs in observed:
            if start_sample <= obs.start_sample <= end_sample:
                ax.axvline(obs.start_s, color="#bbbbbb", linewidth=0.7, alpha=0.8)
                ax.text(
                    obs.start_s,
                    np.nanmax(trig_data) * 0.95 if trig_data.size else 1.0,
                    str(obs.trigger),
                    rotation=90,
                    va="top",
                    ha="center",
                    fontsize=7,
                    color="#444444",
                )

        ax.set_xlabel("Time (s)")
        ax.set_ylabel(trigger_channel)
        ax.set_title(
            f"{anom.anomaly_type} (row {anom.align_row_1based}) | "
            f"anchor={anom.anchor_s:.3f}s | {anom.detail}"
        )
        ax.grid(True, alpha=0.25)
        ax.set_xlim(start_sample / sfreq, end_sample / sfreq)
        fig.tight_layout()

        filename = (
            f"{anom.anomaly_idx_1based:03d}_"
            f"{_safe_slug(anom.anomaly_type)}_"
            f"row{anom.align_row_1based:03d}.png"
        )
        plot_path = out_dir / filename
        fig.savefig(plot_path, dpi=150)
        plt.close(fig)

        anom.window_start_s = start_sample / sfreq
        anom.window_end_s = end_sample / sfreq
        anom.plot_path = str(plot_path)


def write_anomaly_csv(path: Path, anomalies: list[Anomaly]) -> None:
    """Write anomaly summary rows, including plot file paths."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "anomaly_idx",
        "anomaly_type",
        "align_row",
        "detail",
        "reported_idx",
        "reported_trigger",
        "reported_label",
        "observed_idx",
        "observed_trigger",
        "observed_width_samples",
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
                    "reported_trigger": _fmt_optional_int(a.reported_trigger),
                    "reported_label": a.reported_label,
                    "observed_idx": _fmt_optional_int(a.observed_idx),
                    "observed_trigger": _fmt_optional_int(a.observed_trigger),
                    "observed_width_samples": _fmt_optional_int(a.observed_width_samples),
                    "anchor_s": f"{a.anchor_s:.6f}" if np.isfinite(a.anchor_s) else "",
                    "anchor_source": a.anchor_source,
                    "window_start_s": (
                        f"{a.window_start_s:.6f}" if a.window_start_s is not None else ""
                    ),
                    "window_end_s": (
                        f"{a.window_end_s:.6f}" if a.window_end_s is not None else ""
                    ),
                    "plot_path": a.plot_path,
                }
            )


# -----------------------------------------------------------------------------
# Output formatting and writing
# -----------------------------------------------------------------------------


def _fmt_optional_int(value: Optional[int]) -> str:
    """Format optional int as text for CSV/Markdown output."""
    return "" if value is None else str(value)


def build_output_rows(
    ops: list[AlignmentOp],
    reported: list[ReportedTrigger],
    observed: list[ObservedWindow],
) -> tuple[list[dict[str, str]], dict[str, int], Optional[int], Optional[int]]:
    """Build row-wise output records and aggregate alignment stats."""
    rows: list[dict[str, str]] = []

    stats = {
        "MATCH": 0,
        "MISMATCH": 0,
        "REPORTED_ONLY": 0,
        "OBSERVED_ONLY": 0,
    }

    observed_indices_used: list[int] = []

    for out_idx, op in enumerate(ops, start=1):
        rep = reported[op.reported_idx0] if op.reported_idx0 is not None else None
        obs = observed[op.observed_idx0] if op.observed_idx0 is not None else None

        if op.op_type == "PAIR":
            assert rep is not None and obs is not None
            status = "MATCH" if rep.trigger == obs.trigger else "MISMATCH"
        elif op.op_type == "REPORTED_ONLY":
            status = "REPORTED_ONLY"
        elif op.op_type == "OBSERVED_ONLY":
            status = "OBSERVED_ONLY"
        else:
            raise ValueError(f"Unknown alignment op type: {op.op_type}")

        stats[status] += 1
        if obs is not None:
            observed_indices_used.append(obs.idx_1based)

        rows.append(
            {
                "align_row": str(out_idx),
                "status": status,
                "reported_idx": _fmt_optional_int(rep.idx_1based if rep else None),
                "reported_trigger": _fmt_optional_int(rep.trigger if rep else None),
                "reported_trial": _fmt_optional_int(rep.trial if rep else None),
                "reported_frame": _fmt_optional_int(rep.frame if rep else None),
                "reported_label": rep.label if rep else "",
                "observed_idx": _fmt_optional_int(obs.idx_1based if obs else None),
                "observed_trigger": _fmt_optional_int(obs.trigger if obs else None),
                "observed_start_sample": _fmt_optional_int(obs.start_sample if obs else None),
                "observed_end_sample": _fmt_optional_int(obs.end_sample if obs else None),
                "observed_start_s": "" if obs is None else f"{obs.start_s:.6f}",
                "observed_end_s": "" if obs is None else f"{obs.end_s:.6f}",
                "observed_width_ms": "" if obs is None else f"{obs.width_ms:.3f}",
            }
        )

    observed_start = min(observed_indices_used) if observed_indices_used else None
    observed_end = max(observed_indices_used) if observed_indices_used else None
    return rows, stats, observed_start, observed_end


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write alignment rows to CSV."""
    if not rows:
        raise ValueError("No rows available for CSV output.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _md_escape(text: str) -> str:
    """Escape markdown table separators in free-text cells."""
    return text.replace("|", "\\|")


def write_markdown_report(
    path: Path,
    *,
    reported_csv: Path,
    fif_path: Path,
    trigger_channel: str,
    out_csv: Path,
    anomaly_csv: Path,
    n_reported: int,
    n_observed: int,
    alignment_score: int,
    stats: dict[str, int],
    observed_start_idx: Optional[int],
    observed_end_idx: Optional[int],
    rows: list[dict[str, str]],
    anomalies: list[Anomaly],
    plot_window_ms: float,
) -> None:
    """Write a human-readable markdown report for inspection."""
    path.parent.mkdir(parents=True, exist_ok=True)

    match_rate = (stats["MATCH"] / n_reported * 100.0) if n_reported else 0.0

    lines: list[str] = []
    lines.append("# Trigger Confrontation Report")
    lines.append("")
    lines.append("## Inputs")
    lines.append(f"- Reported CSV: `{reported_csv}`")
    lines.append(f"- FIF file: `{fif_path}`")
    lines.append(f"- Trigger channel: `{trigger_channel}`")
    lines.append("")
    lines.append("## Summary")
    lines.append(f"- Reported events: **{n_reported}**")
    lines.append(f"- Observed windows in FIF: **{n_observed}**")
    lines.append(f"- Alignment score: **{alignment_score}**")
    lines.append(f"- MATCH: **{stats['MATCH']}**")
    lines.append(f"- MISMATCH: **{stats['MISMATCH']}**")
    lines.append(f"- REPORTED_ONLY (missing on FIF): **{stats['REPORTED_ONLY']}**")
    lines.append(f"- OBSERVED_ONLY (extra on FIF inside aligned span): **{stats['OBSERVED_ONLY']}**")
    lines.append(f"- Reported match rate: **{match_rate:.1f}%**")
    if observed_start_idx is not None and observed_end_idx is not None:
        lines.append(
            "- Aligned observed window index span: "
            f"**{observed_start_idx}..{observed_end_idx}**"
        )
    lines.append("")
    lines.append("## Artifacts")
    lines.append(f"- Full alignment CSV: `{out_csv}`")
    lines.append(f"- Anomaly CSV: `{anomaly_csv}`")
    if anomalies:
        plot_dir = Path(anomalies[0].plot_path).parent if anomalies[0].plot_path else None
        if plot_dir is not None:
            lines.append(f"- Anomaly plots directory: `{plot_dir}`")
    lines.append("")

    # Highlight non-MATCH rows first for quick triage.
    lines.append("## Non-MATCH Rows")
    non_match_rows = [r for r in rows if r["status"] != "MATCH"]
    if not non_match_rows:
        lines.append("All aligned rows are MATCH.")
    else:
        lines.append("| align_row | status | reported_trigger | reported_label | observed_trigger | observed_start_s | observed_end_s |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for r in non_match_rows:
            lines.append(
                "| "
                + " | ".join(
                    [
                        r["align_row"],
                        r["status"],
                        r["reported_trigger"] or "",
                        _md_escape(r["reported_label"]),
                        r["observed_trigger"] or "",
                        r["observed_start_s"] or "",
                        r["observed_end_s"] or "",
                    ]
                )
                + " |"
            )
    lines.append("")

    # Summarize anomalies by type and include plot references.
    lines.append("## Anomaly Plots")
    lines.append(f"- Trigger window around each anomaly: ±{plot_window_ms:.1f} ms.")
    if not anomalies:
        lines.append("No anomalies were detected; no plots generated.")
    else:
        type_counts = Counter(a.anomaly_type for a in anomalies)
        lines.append(
            "- Anomaly counts: "
            + ", ".join(f"{k}={v}" for k, v in sorted(type_counts.items()))
        )
        lines.append("")
        lines.append(
            "| anomaly_idx | type | align_row | reported_trigger | observed_trigger | "
            "anchor_s | window_start_s | window_end_s | plot_path |"
        )
        lines.append(
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"
        )
        for a in anomalies:
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(a.anomaly_idx_1based),
                        a.anomaly_type,
                        str(a.align_row_1based),
                        _fmt_optional_int(a.reported_trigger),
                        _fmt_optional_int(a.observed_trigger),
                        (f"{a.anchor_s:.6f}" if np.isfinite(a.anchor_s) else ""),
                        (
                            f"{a.window_start_s:.6f}"
                            if a.window_start_s is not None
                            else ""
                        ),
                        (
                            f"{a.window_end_s:.6f}"
                            if a.window_end_s is not None
                            else ""
                        ),
                        _md_escape(a.plot_path),
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
            lines.append(f"- Detail: {a.detail}")
            lines.append(f"- Anchor source: `{a.anchor_source}`")
            lines.append(f"![Anomaly {a.anomaly_idx_1based}]({a.plot_path})")
            lines.append("")

    # Full alignment table kept after anomaly summaries for complete traceability.
    lines.append("")

    # Full table for complete inspection.
    lines.append("## Full Alignment Table")
    lines.append(
        "| align_row | status | reported_idx | reported_trigger | reported_trial | "
        "reported_frame | reported_label | observed_idx | observed_trigger | "
        "observed_start_s | observed_end_s | observed_width_ms |"
    )
    lines.append(
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
    )
    for r in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    r["align_row"],
                    r["status"],
                    r["reported_idx"],
                    r["reported_trigger"],
                    r["reported_trial"],
                    r["reported_frame"],
                    _md_escape(r["reported_label"]),
                    r["observed_idx"],
                    r["observed_trigger"],
                    r["observed_start_s"],
                    r["observed_end_s"],
                    r["observed_width_ms"],
                ]
            )
            + " |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# -----------------------------------------------------------------------------
# CLI entrypoint
# -----------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare reported debug trigger CSV against observed FIF trigger windows "
            "and write CSV + Markdown confrontation outputs."
        )
    )
    parser.add_argument("--reported-csv", type=Path, required=True, help="Path to debug trigger CSV.")
    parser.add_argument("--fif", type=Path, required=True, help="Path to FIF file.")
    parser.add_argument(
        "--trigger-channel",
        type=str,
        default="STI101",
        help="Stim/trigger channel name in FIF (default: STI101).",
    )
    parser.add_argument(
        "--min-duration",
        type=float,
        default=0.0,
        help="Minimum duration in seconds for event detection merge (default: 0).",
    )
    parser.add_argument("--out-csv", type=Path, required=True, help="Output confrontation CSV path.")
    parser.add_argument("--out-report", type=Path, required=True, help="Output markdown report path.")
    parser.add_argument(
        "--anomaly-csv",
        type=Path,
        default=None,
        help=(
            "Output anomaly CSV path. Defaults to <out_csv stem>_anomalies.csv "
            "in the same directory as --out-csv."
        ),
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=None,
        help=(
            "Directory where anomaly plots are written. Defaults to "
            "<out_report stem>_plots next to --out-report."
        ),
    )
    parser.add_argument(
        "--plot-window-ms",
        type=float,
        default=500.0,
        help="Half-window size in ms around each anomaly anchor (default: 500).",
    )
    parser.add_argument(
        "--no-width-outliers",
        dest="include_width_outliers",
        action="store_false",
        help="Disable WIDTH_OUTLIER anomaly detection based on pulse-width mode.",
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
    parser.set_defaults(include_width_outliers=True)
    return parser


def main() -> None:
    """Run reported-vs-observed confrontation and write outputs."""
    args = build_parser().parse_args()

    # Resolve paths once for stable output provenance.
    reported_csv = args.reported_csv.expanduser().resolve()
    fif_path = args.fif.expanduser().resolve()
    out_csv = args.out_csv.expanduser().resolve()
    out_report = args.out_report.expanduser().resolve()
    anomaly_csv = (
        args.anomaly_csv.expanduser().resolve()
        if args.anomaly_csv is not None
        else out_csv.with_name(f"{out_csv.stem}_anomalies.csv")
    )
    plots_dir = (
        args.plots_dir.expanduser().resolve()
        if args.plots_dir is not None
        else out_report.with_name(f"{out_report.stem}_plots")
    )

    # Load reported debug triggers.
    reported = load_reported_csv(reported_csv)

    # Extract observed trigger windows from FIF stim channel.
    observed, sfreq, first_samp = reconstruct_observed_windows(
        fif_path,
        args.trigger_channel,
        args.min_duration,
    )
    if not observed:
        raise RuntimeError("No observed trigger windows were reconstructed from the FIF.")

    # Align the full reported sequence to the observed sequence.
    ops, alignment_score, end_j = semiglobal_align(reported, observed)
    reported_anchors = build_reported_anchors(ops, reported, observed)

    # Build row-level confrontation output + summary metrics.
    rows, stats, observed_start_idx, observed_end_idx = build_output_rows(ops, reported, observed)

    # Detect anomalies and generate trigger-channel windows around each anomaly.
    anomalies = detect_anomalies(
        ops,
        reported,
        observed,
        include_width_outliers=args.include_width_outliers,
    )
    raw_for_plots = mne.io.read_raw_fif(fif_path, preload=False, verbose="ERROR")
    plot_anomaly_windows(
        raw=raw_for_plots,
        trigger_channel=args.trigger_channel,
        observed=observed,
        reported_anchors=reported_anchors,
        anomalies=anomalies,
        window_ms=args.plot_window_ms,
        out_dir=plots_dir,
    )

    # Persist machine-readable + human-readable artifacts.
    write_csv(out_csv, rows)
    write_anomaly_csv(anomaly_csv, anomalies)
    write_markdown_report(
        out_report,
        reported_csv=reported_csv,
        fif_path=fif_path,
        trigger_channel=args.trigger_channel,
        out_csv=out_csv,
        anomaly_csv=anomaly_csv,
        n_reported=len(reported),
        n_observed=len(observed),
        alignment_score=alignment_score,
        stats=stats,
        observed_start_idx=observed_start_idx,
        observed_end_idx=observed_end_idx,
        rows=rows,
        anomalies=anomalies,
        plot_window_ms=args.plot_window_ms,
    )

    # Console summary for quick terminal feedback.
    print(f"Reported events: {len(reported)}")
    print(f"Observed windows: {len(observed)} (sfreq={sfreq:.3f} Hz, first_samp={first_samp})")
    print(f"Alignment score: {alignment_score} (end_j={end_j})")
    print(
        "Status counts: "
        f"MATCH={stats['MATCH']} "
        f"MISMATCH={stats['MISMATCH']} "
        f"REPORTED_ONLY={stats['REPORTED_ONLY']} "
        f"OBSERVED_ONLY={stats['OBSERVED_ONLY']}"
    )
    print(f"Anomalies detected: {len(anomalies)}")
    print(f"Saved confrontation CSV: {out_csv}")
    print(f"Saved anomaly CSV: {anomaly_csv}")
    print(f"Saved anomaly plots in: {plots_dir}")
    print(f"Saved confrontation report: {out_report}")

    # Optionally open an interactive scroller for manual timeline inspection.
    if args.open_scroller:
        print(
            "Opening interactive scroller. "
            "Use slider or left/right (a/d) keys to navigate."
        )
        open_trigger_scroller(
            raw=raw_for_plots,
            trigger_channel=args.trigger_channel,
            reported_anchors=reported_anchors,
            duration_s=args.scroller_duration,
            start_s=args.scroller_start,
            step_s=args.scroller_step,
        )


if __name__ == "__main__":
    main()
