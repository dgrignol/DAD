#!/usr/bin/env python3
"""
deep_trigger_report.py

Purpose
Render the legacy subject-level reported-vs-observed trigger report bundle
from the canonical outputs of ``trigger_extraction.py``. This preserves the
familiar report filenames, block/run tables, anomaly plots, and transition
plots without independently rediscovering or realigning FIF and CSV files.

Data flow
1. Receive canonical in-memory events, alignment rows, transition warnings,
   behavioral segments, and segment summaries from step 02.
2. Convert them to the legacy tabular presentation schema.
   Preserve reviewed acquisition boundaries and unavailable-tail labels.
   Resolve rows already identified by the canonical parser as unavailable
   tails of reviewed truncated acquisitions, without creating observations.
   A reported-only event may receive a separate Stage 1 adjudication when a
   matched gaze-break trigger follows on the next logged frame in that same
   trial, provided the configured policy excludes gaze-break trials.
3. Read only short STI101 windows from the source FIF named by each anomaly or
   transition warning. Repeated rows in one unavailable recording tail share
   one representative cutoff plot. Transition plots show the canonical
   alignment row for the retained pulse and, when the next reported event is
   the dropped code, the corresponding reported-only candidate row. They also
   show each transition code as an 8-bit binary value with its decimal value.
   Reviewed pairs that preserve both pulses use keep1/keep2 labels instead of
   drop.
4. Write the alignment CSV, anomaly CSV, plots, and Markdown report under the
   configured derivatives report directory.

Usage
This module is called by ``02_trigger_extraction_custom.py``. A direct Python
example for tests or diagnostics is:

    from deep_trigger_report import generate_legacy_deep_report
    result = generate_legacy_deep_report(
        subject_id="04",
        raw_root="/Volumes/MORWUR/Projects/DAMIANO/Attention/data",
        output_dir="/tmp/deep-report",
        stim_channel="STI101",
        events=events,
        alignment_rows=alignment_rows,
        transition_warnings=transition_warnings,
        segments=segments,
        segment_summaries=segment_summaries,
        gaze_break_codes=[150],
        gaze_break_exclusion_enabled=True,
    )

Safety and assumptions
- The canonical step-02 alignment is authoritative; this module never zips
  acquisition FIFs and behavioral blocks or performs a second alignment.
- Source FIF files are opened read-only and never preloaded.
- Every output path is checked with the shared raw-tree guard before writing.
- Plot failures are reported per artifact and do not remove the underlying
  machine-readable anomaly or transition row.
- A matched next-frame gaze break can resolve a reported-only row for trial
  inclusion only when the configured exclusion policy excludes that trial.
  The original `REPORTED_ONLY` alignment status is preserved and remains in
  the raw alignment counts.
"""

from __future__ import annotations

import csv
import math
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import mne
import numpy as np

from pipeline_config import assert_output_path_outside_raw


# -----------------------------------------------------------------------------
# Generic formatting and output helpers
# -----------------------------------------------------------------------------


def _safe_text(value: Any) -> str:
    """Return an empty string for null values and text otherwise."""
    return "" if value is None else str(value)


def _markdown_escape(value: Any) -> str:
    """Escape free text for a Markdown table cell."""
    return _safe_text(value).replace("|", "\\|").replace("\n", " ")


def _format_trigger_code(code: Any) -> str:
    """Format a byte-range trigger as zero-padded binary and decimal."""
    value = int(code)
    return f"{value:08b} ({value})"


def _transition_alignment_label(warning: Mapping[str, Any]) -> str:
    """Format retained and dropped-code alignment row references for reports."""
    labels = []
    if warning.get("kept_align_row") is not None:
        first_label = "keep1" if warning.get("retained_both") else "kept"
        labels.append(f"{first_label} {warning['kept_align_row']}")
    if warning.get("dropped_align_row") is not None:
        if warning.get("retained_both"):
            labels.append(f"keep2 {warning['dropped_align_row']}")
        else:
            labels.append(f"dropped-code candidate {warning['dropped_align_row']}")
    return "; ".join(labels) if labels else "unavailable"


def _status_color(status: str) -> str:
    """Map an alignment status to the legacy report color."""
    if "WARNING" in status:
        return "#f5e7a1"
    if status == "MATCH":
        return "#b6e3b6"
    if status in {"MISMATCH", "REPORTED_ONLY", "OBSERVED_ONLY"}:
        return "#f3b0b0"
    return "#dddddd"


def _status_badge(status: str) -> str:
    """Render the legacy colored status badge."""
    return (
        f'<span style="background:{_status_color(status)};padding:2px 6px;'
        f'border-radius:4px;border:1px solid #999;">'
        f"{_markdown_escape(status)}</span>"
    )


def _relative_to_report(report_path: Path, target: Path) -> str:
    """Return a stable path relative to the Markdown report directory."""
    return Path(os.path.relpath(target, start=report_path.parent)).as_posix()


def _write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    """Write a stable CSV, including its header when no rows exist."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# -----------------------------------------------------------------------------
# Canonical-to-legacy table conversion
# -----------------------------------------------------------------------------


def _event_lookup(events: Iterable[Mapping[str, Any]]) -> dict[int, Mapping[str, Any]]:
    """Index canonical observed events by one-based event identifier."""
    return {int(event["event_index"]): event for event in events}


def build_legacy_alignment_rows(
    alignment_rows: list[dict[str, Any]],
    events: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """
    Convert canonical alignment rows to the legacy full-table schema.

    Segment markers are derived only from reported rows. Observed-only rows
    retain inferred behavioral context but do not become synthetic run starts.
    Canonical unavailable-tail reasons tied to reviewed acquisition boundaries
    resolve the report adjudication while retaining the raw REPORTED_ONLY status.
    """
    events_by_index = _event_lookup(events)
    seen_blocks: set[int] = set()
    seen_segments: set[tuple[int, int]] = set()
    reported_index = 0
    legacy_rows: list[dict[str, Any]] = []

    for row in alignment_rows:
        observed_index = row.get("observed_event_index")
        observed = (
            events_by_index.get(int(observed_index))
            if observed_index is not None and _safe_text(observed_index)
            else None
        )
        block = row.get("behavioral_block")
        run = row.get("behavioral_run")
        markers: list[str] = []
        has_reported = row.get("reported_trigger_code") is not None
        if has_reported:
            reported_index += 1
            if block is not None and int(block) not in seen_blocks:
                markers.append("BLOCK_START")
                seen_blocks.add(int(block))
            if block is not None and run is not None:
                segment_key = (int(block), int(run))
                if segment_key not in seen_segments:
                    markers.append("RUN_START")
                    seen_segments.add(segment_key)

        # The canonical reason comes from the reviewed truncated-acquisition
        # profile. Report that accepted explanation without inventing a pulse.
        adjudication = ""
        adjudication_detail = ""
        if row["status"] == "REPORTED_ONLY":
            adjudication = "REVIEW_REQUIRED"
            if (
                row.get("reported_only_reason") == "recording_tail_unavailable"
                and row.get("constrained_acquisition_run") is not None
                and observed_index is None
            ):
                adjudication = "RESOLVED_RECORDING_TAIL_UNAVAILABLE"
                adjudication_detail = (
                    "Accepted truncated-acquisition boundary explains this logged "
                    "event: the recording ended before it occurred. No observed "
                    "MEG event is created; raw REPORTED_ONLY status is retained."
                )

        legacy_rows.append(
            {
                "align_row": int(row["alignment_index"]),
                "status": row["status"],
                "segment_marker": ",".join(markers),
                "reported_idx": reported_index if has_reported else "",
                "reported_block": block,
                "reported_run": run,
                "reported_trigger": row.get("reported_trigger_code"),
                "reported_trial": row.get("behavioral_trial"),
                "reported_frame": row.get("reported_frame"),
                "reported_seconds": row.get("reported_seconds"),
                "reported_label": row.get("reported_label"),
                "reported_source_csv": (
                    Path(row["reported_source_csv"]).name
                    if row.get("reported_source_csv")
                    else ""
                ),
                "reported_replay_status": row.get("reported_replay_status"),
                "observed_idx": row.get("observed_event_index"),
                "observed_trigger": row.get("observed_trigger_code"),
                "observed_start_s": (
                    f"{float(observed['time_start_sec']):.6f}"
                    if observed is not None
                    else ""
                ),
                "observed_end_s": (
                    f"{float(observed['time_end_sec']):.6f}"
                    if observed is not None
                    else ""
                ),
                "observed_width_ms": (
                    f"{1000.0 * (float(observed['time_end_sec']) - float(observed['time_start_sec'])):.3f}"
                    if observed is not None
                    else ""
                ),
                "observed_source_fif": row.get("source_fif"),
                "acquisition_run": row.get("acquisition_run"),
                "raw_sample_frame": row.get("raw_sample_frame"),
                "constrained_acquisition_run": row.get("constrained_acquisition_run"),
                "reported_only_reason": row.get("reported_only_reason"),
                "stage1_adjudication": adjudication,
                "adjudication_partner_align_row": "",
                "adjudication_gap_ms": "",
                "adjudication_detail": adjudication_detail,
                "behavioral_context_source": row.get(
                    "behavioral_context_source"
                ),
            }
        )
    return legacy_rows


def _same_reported_trial(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    """Return whether two report rows identify the same trial instance."""
    required_fields = ("reported_block", "reported_run", "reported_trial")
    if any(left.get(field) in (None, "") for field in required_fields):
        return False
    if any(
        _safe_text(left.get(field)) != _safe_text(right.get(field))
        for field in required_fields
    ):
        return False
    for field in ("reported_replay_status", "reported_source_csv"):
        left_value = _safe_text(left.get(field))
        right_value = _safe_text(right.get(field))
        if left_value and right_value and left_value != right_value:
            return False
    return True


def _annotate_gaze_break_exclusions(
    legacy_rows: list[dict[str, Any]],
    *,
    gaze_break_codes: Iterable[int],
    gaze_break_exclusion_enabled: bool,
) -> list[dict[str, Any]]:
    """Flag a missing trigger when its matched next-frame gaze break excludes the trial.

    The raw alignment status remains ``REPORTED_ONLY``. Automatic adjudication
    requires a matched condition onset in the same trial, the next logged
    reported event to be a matched configured gaze-break code exactly one frame
    later, and the configured trial policy to exclude gaze-break trials. This
    leaves missing gaze-break triggers for review and preserves any separate
    adjudication already attached to a confirmed unavailable recording tail.
    """
    configured_break_codes = {int(code) for code in gaze_break_codes}
    if not gaze_break_exclusion_enabled or not configured_break_codes:
        return legacy_rows

    reported_positions = [
        index
        for index, row in enumerate(legacy_rows)
        if row.get("reported_trigger") not in (None, "")
    ]
    next_reported_position = {
        left: right
        for left, right in zip(reported_positions, reported_positions[1:])
    }

    for position, row in enumerate(legacy_rows):
        if row.get("status") != "REPORTED_ONLY":
            continue
        if row.get("reported_only_reason"):
            continue

        try:
            missing_code = int(row["reported_trigger"])
        except (TypeError, ValueError):
            continue
        # A missing gaze-break event cannot itself prove that a gaze break was
        # recorded in MEG; it remains an acquisition discrepancy for review.
        if missing_code in configured_break_codes:
            continue

        partner_position = next_reported_position.get(position)
        if partner_position is None:
            continue
        partner = legacy_rows[partner_position]
        if not _same_reported_trial(row, partner):
            continue
        if (
            partner.get("status") != "MATCH"
            or partner.get("reported_label") != "gaze_break"
        ):
            continue
        try:
            partner_code = int(partner["reported_trigger"])
            observed_code = int(partner["observed_trigger"])
            current_frame = int(row["reported_frame"])
            partner_frame = int(partner["reported_frame"])
            current_report_index = int(row["reported_idx"])
            partner_report_index = int(partner["reported_idx"])
            gap_ms = (
                float(partner["reported_seconds"])
                - float(row["reported_seconds"])
            ) * 1000.0
        except (KeyError, TypeError, ValueError):
            continue
        if (
            partner_code not in configured_break_codes
            or observed_code not in configured_break_codes
            or partner_frame != current_frame + 1
            or partner_report_index != current_report_index + 1
            or not math.isfinite(gap_ms)
            or gap_ms <= 0
        ):
            continue

        # The interval can only be treated as an excluded trial when its
        # condition onset itself is present in the observed data.
        has_matched_condition_onset = any(
            _same_reported_trial(row, candidate)
            and candidate.get("reported_label") == "trial_condition_onset"
            and candidate.get("status") == "MATCH"
            and candidate.get("observed_trigger") not in (None, "")
            for candidate in legacy_rows[:position]
        )
        if not has_matched_condition_onset:
            continue

        has_intervening_condition_onset = any(
            candidate.get("reported_label") == "trial_condition_onset"
            and _same_reported_trial(row, candidate)
            for candidate in legacy_rows[position + 1 : partner_position]
        )
        if has_intervening_condition_onset:
            continue

        row["stage1_adjudication"] = "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK"
        row["adjudication_partner_align_row"] = partner["align_row"]
        row["adjudication_gap_ms"] = f"{gap_ms:.3f}"
        row["adjudication_detail"] = (
            "Matched gaze-break trigger occurs on the next reported frame "
            "in the same trial before the next condition onset. The configured "
            "trial policy excludes this trial; raw REPORTED_ONLY status is retained."
        )
    return legacy_rows


def build_anomaly_rows(
    legacy_rows: list[dict[str, Any]],
    events: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build anomaly context without anchoring a missing tail in a later FIF.

    A known unavailable tail has no measured pulse to plot. Its contextual
    anchor is the last recorded event within its explicit acquisition boundary.
    The detail text identifies that distinction in the CSV and Markdown.
    """
    events_by_index = _event_lookup(events)
    anomalies: list[dict[str, Any]] = []

    for row_index, row in enumerate(legacy_rows):
        if row["status"] not in {"MISMATCH", "REPORTED_ONLY", "OBSERVED_ONLY"}:
            continue

        observed_index = row.get("observed_idx")
        anchor_event = (
            events_by_index.get(int(observed_index))
            if observed_index is not None and _safe_text(observed_index)
            else None
        )
        anchor_source = "observed_event"
        if anchor_event is None:
            for candidate in legacy_rows[row_index + 1 :]:
                constrained_run = row.get("constrained_acquisition_run")
                if constrained_run is not None and candidate.get("acquisition_run") != constrained_run:
                    continue
                candidate_index = candidate.get("observed_idx")
                if candidate_index is not None and _safe_text(candidate_index):
                    anchor_event = events_by_index[int(candidate_index)]
                    anchor_source = "next_observed"
                    break
        if anchor_event is None:
            for candidate in reversed(legacy_rows[:row_index]):
                constrained_run = row.get("constrained_acquisition_run")
                if constrained_run is not None and candidate.get("acquisition_run") != constrained_run:
                    continue
                candidate_index = candidate.get("observed_idx")
                if candidate_index is not None and _safe_text(candidate_index):
                    anchor_event = events_by_index[int(candidate_index)]
                    anchor_source = "previous_observed"
                    break

        detail = {
            "MISMATCH": "Reported and observed trigger codes differ.",
            "REPORTED_ONLY": "Reported trigger has no aligned observed pulse.",
            "OBSERVED_ONLY": "Observed trigger has no aligned reported entry.",
        }[row["status"]]
        if row.get("reported_only_reason") == "recording_tail_unavailable":
            detail = (
                "Logged event belongs to the unavailable tail of a confirmed "
                "truncated acquisition. Any plot shows the last recorded trigger "
                "as cutoff context; the missing event has no recorded waveform."
            )
            anchor_source = "recorded_cutoff_context"
        if row.get("stage1_adjudication") == "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK":
            detail = _safe_text(row.get("adjudication_detail"))
        anomalies.append(
            {
                "anomaly_idx": len(anomalies) + 1,
                "anomaly_type": row["status"],
                "align_row": row["align_row"],
                "detail": detail,
                "reported_only_reason": row.get("reported_only_reason"),
                "constrained_acquisition_run": row.get("constrained_acquisition_run"),
                "stage1_adjudication": row.get("stage1_adjudication"),
                "adjudication_partner_align_row": row.get(
                    "adjudication_partner_align_row"
                ),
                "adjudication_gap_ms": row.get("adjudication_gap_ms"),
                "adjudication_detail": row.get("adjudication_detail"),
                "reported_idx": row["reported_idx"],
                "reported_trigger": row["reported_trigger"],
                "reported_label": row["reported_label"],
                "reported_block": row["reported_block"],
                "reported_run": row["reported_run"],
                "observed_idx": row["observed_idx"],
                "observed_trigger": row["observed_trigger"],
                "observed_source_fif": row["observed_source_fif"],
                "observed_width_samples": (
                    anchor_event.get("width_samples")
                    if anchor_event is not None
                    else None
                ),
                "observed_start_s": row["observed_start_s"],
                "anchor_source": anchor_source if anchor_event is not None else "",
                "anchor_source_fif_path": (
                    anchor_event.get("source_fif_path")
                    if anchor_event is not None
                    else None
                ),
                "anchor_raw_sample_frame": (
                    anchor_event.get("raw_sample_frame")
                    if anchor_event is not None
                    else None
                ),
                "plot_representative_align_row": "",
                "plot_path": "",
                "plot_error": "",
            }
        )
    return anomalies


# -----------------------------------------------------------------------------
# Read-only short-window plotting
# -----------------------------------------------------------------------------


def _read_stim_window(
    fif_path: Path,
    stim_channel: str,
    center_sample_raw: int,
    half_window_ms: float,
    raw: Any | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Read one bounded trigger window, reusing an open FIF handle when given."""
    owns_raw = raw is None
    if owns_raw:
        raw = mne.io.read_raw_fif(str(fif_path), preload=False, verbose="ERROR")
    try:
        if stim_channel not in raw.ch_names:
            raise ValueError(
                f"Trigger channel {stim_channel!r} missing from {fif_path.name}."
            )
        sfreq = float(raw.info["sfreq"])
        half_samples = max(1, int(round(half_window_ms * sfreq / 1000.0)))
        center_relative = int(center_sample_raw) - int(raw.first_samp)
        start = max(0, center_relative - half_samples)
        stop = min(int(raw.n_times), center_relative + half_samples + 1)
        if stop <= start:
            raise ValueError("Requested plot window is outside the source FIF.")
        data = raw.get_data(picks=[stim_channel], start=start, stop=stop)[0]
        raw_samples = np.arange(start, stop) + int(raw.first_samp)
        times_ms = 1000.0 * (raw_samples - int(center_sample_raw)) / sfreq
        return times_ms, data, sfreq
    finally:
        if owns_raw:
            raw.close()


def _save_stim_plot(
    *,
    output_path: Path,
    fif_path: Path,
    stim_channel: str,
    center_sample_raw: int,
    half_window_ms: float,
    title: str,
    status: str,
    highlighted_span_samples: tuple[int, int] | None = None,
    right_annotation: str | None = None,
    raw: Any | None = None,
) -> None:
    """Render one non-interactive trigger plot and optional right-side key."""
    times_ms, data, sfreq = _read_stim_window(
        fif_path,
        stim_channel,
        center_sample_raw,
        half_window_ms,
        raw=raw,
    )
    figure = Figure(figsize=(13, 3.6) if right_annotation else (10, 3.4))
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    axis.step(times_ms, data, where="post", color="#1f77b4", linewidth=1.0)
    axis.axvline(0.0, color="#d62728", linestyle="--", linewidth=1.2)
    if highlighted_span_samples is not None:
        span_start = (
            1000.0
            * (highlighted_span_samples[0] - center_sample_raw)
            / sfreq
        )
        span_end = (
            1000.0
            * (highlighted_span_samples[1] - center_sample_raw)
            / sfreq
        )
        axis.axvspan(span_start, span_end, color="#ff7f0e", alpha=0.25)
    axis.text(
        0.01,
        0.97,
        f"STATUS: {status}",
        transform=axis.transAxes,
        va="top",
        fontsize=8,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": _status_color(status),
            "edgecolor": "#777777",
            "alpha": 0.95,
        },
    )
    axis.set_title(title)
    axis.set_xlabel("Time from anchor (ms)")
    axis.set_ylabel(stim_channel)
    axis.grid(True, alpha=0.25)
    if right_annotation:
        # Reserve a fixed panel beside the waveform so the two binary codes
        # remain stacked and aligned in every transition plot.
        figure.subplots_adjust(left=0.08, right=0.73, bottom=0.18, top=0.83)
        figure.text(
            0.77,
            0.68,
            right_annotation,
            ha="left",
            va="top",
            fontsize=11,
            family="monospace",
            linespacing=1.6,
        )
    else:
        figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=150)
    figure.clear()


def _render_anomaly_plots(
    anomalies: list[dict[str, Any]],
    *,
    plots_dir: Path,
    stim_channel: str,
    half_window_ms: float,
) -> None:
    """Render anomaly plots, sharing one cutoff plot per unavailable tail.

    Tail rows that refer to the same acquisition boundary, behavioral run, and
    recorded cutoff use the first row's plot. Existing duplicate images from
    earlier report generations are removed as each duplicate is encountered.
    """
    raw_cache: dict[str, Any] = {}
    tail_representatives: dict[tuple[str, ...], dict[str, str]] = {}
    try:
        for anomaly in anomalies:
            anomaly["plot_representative_align_row"] = ""
            fif_text = anomaly.get("anchor_source_fif_path")
            raw_sample = anomaly.get("anchor_raw_sample_frame")
            if not fif_text or raw_sample is None:
                anomaly["plot_error"] = "No observed source anchor is available."
                continue
            filename = (
                f"{int(anomaly['anomaly_idx']):04d}_"
                f"{str(anomaly['anomaly_type']).lower()}_"
                f"row{int(anomaly['align_row']):05d}.png"
            )
            plot_path = plots_dir / filename
            is_tail = (
                anomaly.get("reported_only_reason")
                == "recording_tail_unavailable"
            )
            tail_key: tuple[str, ...] | None = None
            if is_tail:
                tail_key = (
                    _safe_text(anomaly.get("constrained_acquisition_run")),
                    _safe_text(anomaly.get("reported_block")),
                    _safe_text(anomaly.get("reported_run")),
                    str(Path(fif_text).resolve()),
                    str(int(raw_sample)),
                )
                representative = tail_representatives.get(tail_key)
                if representative is not None:
                    anomaly["plot_path"] = representative["plot_path"]
                    anomaly["plot_representative_align_row"] = representative[
                        "align_row"
                    ]
                    anomaly["detail"] += (
                        " The single cutoff plot for this tail is shown for "
                        f"alignment row {representative['align_row']}."
                    )
                    if plot_path.is_file():
                        plot_path.unlink()
                    continue
                anomaly["plot_representative_align_row"] = str(
                    anomaly["align_row"]
                )
            try:
                if fif_text not in raw_cache:
                    raw_cache[fif_text] = mne.io.read_raw_fif(
                        str(fif_text), preload=False, verbose="ERROR"
                    )
                _save_stim_plot(
                    output_path=plot_path,
                    fif_path=Path(fif_text),
                    stim_channel=stim_channel,
                    center_sample_raw=int(raw_sample),
                    half_window_ms=half_window_ms,
                    title=(
                        f"{anomaly['anomaly_type']} | row {anomaly['align_row']} | "
                        f"block/run {anomaly['reported_block']}/{anomaly['reported_run']}"
                    ),
                    status=str(anomaly["anomaly_type"]),
                    raw=raw_cache[fif_text],
                )
                anomaly["plot_path"] = str(plot_path)
                if is_tail and tail_key is not None:
                    tail_representatives[tail_key] = {
                        "align_row": str(anomaly["align_row"]),
                        "plot_path": str(plot_path),
                    }
            except Exception as exc:  # pylint: disable=broad-except
                anomaly["plot_error"] = str(exc)
    finally:
        for raw in raw_cache.values():
            raw.close()


def _render_transition_plots(
    warnings: list[dict[str, Any]],
    *,
    plots_dir: Path,
    stim_channel: str,
    half_window_ms: float,
) -> list[dict[str, Any]]:
    """Plot each raw pair, labeling reviewed retained pulses keep1/keep2.

    Warning dictionaries retain the historical kept/dropped code fields for
    compatibility. ``retained_both`` determines whether a pulse was actually
    discarded and controls the visible labels and output filename.
    """
    rendered: list[dict[str, Any]] = []
    raw_cache: dict[str, Any] = {}
    try:
        for warning in warnings:
            row = dict(warning)
            state = "_both_retained" if row.get("retained_both") else ""
            filename = (
                f"{int(row['warning_index']):04d}_transition_warning_"
                f"{row['dropped_code']}_to_{row['kept_code']}{state}.png"
            )
            plot_path = plots_dir / filename
            source_path = str(row["source_fif_path"])
            try:
                if source_path not in raw_cache:
                    raw_cache[source_path] = mne.io.read_raw_fif(
                        source_path, preload=False, verbose="ERROR"
                    )
                if row.get("retained_both"):
                    codes = ", ".join(str(code) for code in row["retained_codes"])
                    transition_summary = f"reviewed; retained both codes {codes}"
                    first_label, second_label = "keep1:", "keep2:"
                else:
                    transition_summary = (
                        f"kept {row['kept_code']}, dropped {row['dropped_code']}"
                    )
                    first_label, second_label = "kept:", "drop:"
                _save_stim_plot(
                    output_path=plot_path,
                    fif_path=Path(source_path),
                    stim_channel=stim_channel,
                    center_sample_raw=int(row["dropped_start_sample_raw"]),
                    half_window_ms=half_window_ms,
                    title=(
                        f"Transition {row['order_pattern']} | {transition_summary}\n"
                        f"Alignment rows: {_transition_alignment_label(row)}"
                    ),
                    status="WARNING",
                    right_annotation=(
                        f"{first_label} {_format_trigger_code(row['kept_code'])}\n"
                        f"{second_label} {_format_trigger_code(row['dropped_code'])}"
                    ),
                    highlighted_span_samples=(
                        int(row["dropped_start_sample_raw"]),
                        int(row["dropped_end_sample_raw"]),
                    ),
                    raw=raw_cache[source_path],
                )
                row["plot_path"] = str(plot_path)
                row["plot_error"] = ""
            except Exception as exc:  # pylint: disable=broad-except
                row["plot_path"] = ""
                row["plot_error"] = str(exc)
            rendered.append(row)
    finally:
        for raw in raw_cache.values():
            raw.close()
    return rendered


def _add_transition_alignment_references(
    warnings: list[dict[str, Any]],
    events: list[dict[str, Any]],
    alignment_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Attach alignment row numbers for each kept pulse and adjacent dropped code."""
    event_by_source_sample: dict[tuple[str, int, int], int] = {}
    for event in events:
        sample = event.get("sample_start_raw", event.get("raw_sample_frame"))
        if sample is None:
            continue
        event_by_source_sample[
            (
                str(event.get("source_fif_path", "")),
                int(sample),
                int(event["trigger_code"]),
            )
        ] = int(event["event_index"])

    alignment_by_event = {
        int(row["observed_event_index"]): row
        for row in alignment_rows
        if row.get("observed_event_index") is not None
    }
    ordered_alignment = sorted(
        alignment_rows, key=lambda row: int(row["alignment_index"])
    )
    referenced: list[dict[str, Any]] = []
    for warning in warnings:
        row = dict(warning)
        kept_event_index = event_by_source_sample.get(
            (
                str(row.get("source_fif_path", "")),
                int(row["kept_start_sample_raw"]),
                int(row["kept_code"]),
            )
        )
        kept_alignment = alignment_by_event.get(kept_event_index)
        if kept_alignment is not None:
            row["kept_align_row"] = int(kept_alignment["alignment_index"])

            if row.get("retained_both"):
                dropped_event_index = event_by_source_sample.get(
                    (
                        str(row.get("source_fif_path", "")),
                        int(row["dropped_start_sample_raw"]),
                        int(row["dropped_code"]),
                    )
                )
                dropped_alignment = alignment_by_event.get(dropped_event_index)
                if dropped_alignment is not None:
                    row["dropped_align_row"] = int(
                        dropped_alignment["alignment_index"]
                    )
            else:
                # A dropped pulse can be an independently reported trigger. Link it
                # only when it is the next reported row and alignment left it unmatched.
                next_reported = next(
                    (
                        candidate
                        for candidate in ordered_alignment
                        if int(candidate["alignment_index"])
                        > int(kept_alignment["alignment_index"])
                        and candidate.get("reported_trigger_code") is not None
                    ),
                    None,
                )
                if (
                    next_reported is not None
                    and next_reported.get("observed_event_index") is None
                    and int(next_reported["reported_trigger_code"])
                    == int(row["dropped_code"])
                ):
                    row["dropped_align_row"] = int(next_reported["alignment_index"])
        referenced.append(row)
    return referenced


# -----------------------------------------------------------------------------
# Markdown report construction
# -----------------------------------------------------------------------------


def _first_rows_by_segment(
    rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return first reported alignment row for each block and block/run."""
    block_rows: list[dict[str, Any]] = []
    run_rows: list[dict[str, Any]] = []
    seen_blocks: set[int] = set()
    seen_runs: set[tuple[int, int]] = set()
    for row in rows:
        if not row["reported_idx"]:
            continue
        block = int(row["reported_block"])
        run = int(row["reported_run"])
        if block not in seen_blocks:
            block_rows.append(row)
            seen_blocks.add(block)
        if (block, run) not in seen_runs:
            run_rows.append(row)
            seen_runs.add((block, run))
    return block_rows, run_rows


def _append_start_table(
    lines: list[str],
    heading: str,
    rows: list[dict[str, Any]],
    *,
    include_run: bool,
) -> None:
    """Append a legacy block-start or run-start table."""
    lines.extend(["", f"## {heading}"])
    run_column = " | run" if include_run else " | first_run"
    lines.append(
        f"| block{run_column} | source_csv | source_fif | reported_start_idx | "
        "align_row | status | observed_idx | observed_start_s |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    _safe_text(row["reported_block"]),
                    _safe_text(row["reported_run"]),
                    _markdown_escape(row["reported_source_csv"]),
                    _markdown_escape(row["observed_source_fif"]),
                    _safe_text(row["reported_idx"]),
                    _safe_text(row["align_row"]),
                    _status_badge(str(row["status"])),
                    _safe_text(row["observed_idx"]),
                    _safe_text(row["observed_start_s"]),
                ]
            )
            + " |"
        )


def _write_markdown_report(
    report_path: Path,
    *,
    subject_id: str,
    legacy_rows: list[dict[str, Any]],
    anomalies: list[dict[str, Any]],
    transition_warnings: list[dict[str, Any]],
    segments: list[Any],
    segment_summaries: list[dict[str, Any]],
    alignment_csv: Path,
    anomalies_csv: Path,
    plots_dir: Path,
    anomaly_window_ms: float,
    transition_window_ms: float,
) -> None:
    """Write the complete legacy-style deep Markdown report."""
    status_counts = Counter(str(row["status"]) for row in legacy_rows)
    reported_count = sum(bool(row["reported_idx"]) for row in legacy_rows)
    observed_count = sum(bool(row["observed_idx"]) for row in legacy_rows)
    matched = status_counts["MATCH"] + status_counts["MATCH_WARNING"]
    match_rate = 100.0 * matched / max(1, reported_count)
    auto_resolved_gaze_break = sum(
        row.get("stage1_adjudication")
        == "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK"
        for row in legacy_rows
    )
    resolved_recording_tail = sum(
        row.get("stage1_adjudication") == "RESOLVED_RECORDING_TAIL_UNAVAILABLE"
        for row in legacy_rows
    )
    reported_only_review_required = sum(
        row.get("status") == "REPORTED_ONLY"
        and row.get("stage1_adjudication") == "REVIEW_REQUIRED"
        for row in legacy_rows
    )
    score_total = sum(
        int(summary.get("alignment_score", 0))
        for summary in segment_summaries
    )
    block_rows, run_rows = _first_rows_by_segment(legacy_rows)

    lines = [
        "# Subject Trigger Confrontation Report",
        "",
        "## Inputs",
        f"- Subject: **sub{subject_id}**",
        f"- Discovered reported CSV files: **{len(segments)}**",
        f"- Discovered FIF files used: **{len({row['observed_source_fif'] for row in legacy_rows if row['observed_source_fif']})}**",
        "- Alignment source: canonical step-02 segment-aware alignment",
        f"- Logged rows in unavailable recording tails: **{sum(row.get('reported_only_reason') == 'recording_tail_unavailable' for row in legacy_rows)}**",
        "",
        "## Summary",
        f"- Total reported events: **{reported_count}**",
        f"- Total observed windows (inside alignment rows): **{observed_count}**",
        f"- Aggregate alignment score (sum across runs): **{score_total}**",
        f"- MATCH: **{status_counts['MATCH']}**",
        f"- MATCH_WARNING: **{status_counts['MATCH_WARNING']}**",
        f"- MISMATCH: **{status_counts['MISMATCH']}**",
        f"- REPORTED_ONLY: **{status_counts['REPORTED_ONLY']}**",
        "- REPORTED_ONLY resolved for trial exclusion by matched gaze break: "
        f"**{auto_resolved_gaze_break}**",
        "- REPORTED_ONLY resolved as confirmed unavailable recording tail: "
        f"**{resolved_recording_tail}**",
        f"- REPORTED_ONLY still requiring review: **{reported_only_review_required}**",
        f"- OBSERVED_ONLY: **{status_counts['OBSERVED_ONLY']}**",
        f"- Reported match rate (MATCH + MATCH_WARNING): **{match_rate:.1f}%**",
        f"- Transition warnings (collapsed pulses): **{len(transition_warnings)}**",
        "",
        "## Artifacts",
        f"- Alignment CSV: `{_relative_to_report(report_path, alignment_csv)}`",
        f"- Anomaly CSV: `{_relative_to_report(report_path, anomalies_csv)}`",
        f"- Anomaly plots directory: `{_relative_to_report(report_path, plots_dir)}`",
        f"- Transition warning plots directory: `{_relative_to_report(report_path, plots_dir / 'transition_pulses')}`",
        f"- This report: `{report_path.name}`",
    ]
    _append_start_table(lines, "Block Starts", block_rows, include_run=False)
    _append_start_table(lines, "Run Starts", run_rows, include_run=True)

    lines.extend(
        [
            "",
            "## Anomaly Plots",
            f"- Plot window: +/-{anomaly_window_ms:.1f} ms; rows in one unavailable "
            "recording tail share a single cutoff plot.",
        ]
    )
    if not anomalies:
        lines.append("No anomalies detected.")
    else:
        lines.extend(
            [
                "",
                "| anomaly_idx | status | align_row | Stage 1 adjudication | paired 150 row | "
                "gap_ms | reported_block | reported_run | reported_trigger | "
                "observed_trigger | anchor_source | plot_for_align_row | plot_path |",
                "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for anomaly in anomalies:
            plot_path = (
                _relative_to_report(report_path, Path(anomaly["plot_path"]))
                if anomaly["plot_path"]
                else ""
            )
            lines.append(
                "| "
                + " | ".join(
                    [
                        _safe_text(anomaly["anomaly_idx"]),
                        _status_badge(str(anomaly["anomaly_type"])),
                        _safe_text(anomaly["align_row"]),
                        _safe_text(anomaly["stage1_adjudication"]),
                        _safe_text(anomaly["adjudication_partner_align_row"]),
                        _safe_text(anomaly["adjudication_gap_ms"]),
                        _safe_text(anomaly["reported_block"]),
                        _safe_text(anomaly["reported_run"]),
                        _safe_text(anomaly["reported_trigger"]),
                        _safe_text(anomaly["observed_trigger"]),
                        _safe_text(anomaly["anchor_source"]),
                        _safe_text(anomaly["plot_representative_align_row"]),
                        _markdown_escape(plot_path),
                    ]
                )
                + " |"
            )
        embedded_anomaly_plots: set[str] = set()
        for anomaly in anomalies:
            lines.extend(
                [
                    "",
                    f"### Anomaly {anomaly['anomaly_idx']}: "
                    f"{anomaly['anomaly_type']} (row {anomaly['align_row']})",
                    f"- Status: {_status_badge(str(anomaly['anomaly_type']))}",
                    f"- Detail: {_markdown_escape(anomaly['detail'])}",
                ]
            )
            if anomaly.get("stage1_adjudication"):
                lines.append(
                    f"- Stage 1 adjudication: `{anomaly['stage1_adjudication']}`"
                )
                if anomaly.get("adjudication_partner_align_row"):
                    lines.append(
                        "- Matched gaze-break alignment row: "
                        f"{anomaly['adjudication_partner_align_row']} "
                        f"({anomaly['adjudication_gap_ms']} ms later)."
                    )
            if anomaly["plot_path"]:
                relative_plot = _relative_to_report(
                    report_path, Path(anomaly["plot_path"])
                )
                if relative_plot in embedded_anomaly_plots:
                    lines.append(
                        "- Shared cutoff plot: "
                        f"[alignment row {anomaly['plot_representative_align_row']}]"
                        f"({relative_plot})"
                    )
                else:
                    lines.append(
                        f"![Anomaly {anomaly['anomaly_idx']}]({relative_plot})"
                    )
                    embedded_anomaly_plots.add(relative_plot)
            elif anomaly["plot_error"]:
                lines.append(f"- Plot error: `{anomaly['plot_error']}`")

    non_matches = [
        row for row in legacy_rows if row["status"] != "MATCH"
    ]
    lines.extend(
        [
            "",
            "## Non-MATCH Rows",
            "| align_row | segment_marker | status | Stage 1 adjudication | paired 150 row | "
            "gap_ms | reported_block | reported_run | reported_trigger | "
            "reported_label | observed_trigger | observed_start_s | "
            "observed_source_fif |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in non_matches:
        lines.append(
            "| "
            + " | ".join(
                [
                    _safe_text(row["align_row"]),
                    _safe_text(row["segment_marker"]),
                    _status_badge(str(row["status"])),
                    _safe_text(row["stage1_adjudication"]),
                    _safe_text(row["adjudication_partner_align_row"]),
                    _safe_text(row["adjudication_gap_ms"]),
                    _safe_text(row["reported_block"]),
                    _safe_text(row["reported_run"]),
                    _safe_text(row["reported_trigger"]),
                    _markdown_escape(row["reported_label"]),
                    _safe_text(row["observed_trigger"]),
                    _safe_text(row["observed_start_s"]),
                    _markdown_escape(row["observed_source_fif"]),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Transition Warnings",
            f"- Plot window around each warning: +/-{transition_window_ms:.1f} ms.",
        ]
    )
    if not transition_warnings:
        lines.append("No transition warnings detected.")
    else:
        lines.extend(
            [
                "",
                "| idx | align_row reference | status | dropped / other (binary (decimal)) | kept / first (binary (decimal)) | expected | decision_source | unresolved | confidence | raw_sample_frame | source_fif | detail | plot_path |",
                "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for warning in transition_warnings:
            plot_path = (
                _relative_to_report(report_path, Path(warning["plot_path"]))
                if warning.get("plot_path")
                else ""
            )
            lines.append(
                "| "
                + " | ".join(
                    [
                        _safe_text(warning["warning_index"]),
                        _markdown_escape(_transition_alignment_label(warning)),
                        _status_badge("WARNING"),
                        _format_trigger_code(warning["dropped_code"]),
                        _format_trigger_code(warning["kept_code"]),
                        _safe_text(warning["expected_trigger_code"]),
                        _safe_text(warning["decision_source"]),
                        _safe_text(warning["unresolved_transition"]),
                        _safe_text(warning["transition_resolution_confidence"]),
                        _safe_text(warning["dropped_start_sample_raw"]),
                        _markdown_escape(warning["source_fif"]),
                        _markdown_escape(
                            (
                                f"Reviewed transition {warning['order_pattern']}; "
                                f"retained both codes {warning['retained_codes']}. "
                                f"{warning.get('review_reason', '')}"
                                if warning.get("retained_both")
                                else (
                                    f"Collapsed {warning['order_pattern']}; kept "
                                    f"{warning['kept_code']} and dropped "
                                    f"{warning['dropped_code']}."
                                )
                            )
                        ),
                        _markdown_escape(plot_path),
                    ]
                )
                + " |"
            )
        for warning in transition_warnings:
            if not warning.get("plot_path"):
                continue
            lines.extend(
                [
                    "",
                    f"### Transition Warning {warning['warning_index']}: {warning['order_pattern']}",
                    f"- Status: {_status_badge('WARNING')}",
                    f"- Alignment rows: {_transition_alignment_label(warning)}",
                    f"![Transition Warning {warning['warning_index']}]"
                    f"({_relative_to_report(report_path, Path(warning['plot_path']))})",
                ]
            )

    lines.extend(
        [
            "",
            "## Full Alignment Table",
            "- Status colors: green = MATCH, red = MISMATCH/REPORTED_ONLY/OBSERVED_ONLY, yellow = MATCH_WARNING.",
            "",
            "| align_row | segment_marker | status | reported_idx | reported_block | reported_run | reported_trigger | reported_trial | reported_frame | reported_seconds | reported_label | reported_source_csv | observed_idx | observed_trigger | observed_start_s | observed_end_s | observed_width_ms | observed_source_fif | acquisition_run | raw_sample_frame |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    table_keys = [
        "align_row",
        "segment_marker",
        "status",
        "reported_idx",
        "reported_block",
        "reported_run",
        "reported_trigger",
        "reported_trial",
        "reported_frame",
        "reported_seconds",
        "reported_label",
        "reported_source_csv",
        "observed_idx",
        "observed_trigger",
        "observed_start_s",
        "observed_end_s",
        "observed_width_ms",
        "observed_source_fif",
        "acquisition_run",
        "raw_sample_frame",
    ]
    for row in legacy_rows:
        values = [
            _status_badge(str(row[key]))
            if key == "status"
            else _markdown_escape(row[key])
            for key in table_keys
        ]
        lines.append("| " + " | ".join(values) + " |")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# -----------------------------------------------------------------------------
# Public renderer
# -----------------------------------------------------------------------------


def generate_legacy_deep_report(
    *,
    subject_id: str,
    raw_root: str | Path,
    output_dir: str | Path,
    stim_channel: str,
    events: list[dict[str, Any]],
    alignment_rows: list[dict[str, Any]],
    transition_warnings: list[dict[str, Any]],
    segments: list[Any],
    segment_summaries: list[dict[str, Any]],
    anomaly_window_ms: float = 500.0,
    transition_window_ms: float = 3.0,
    gaze_break_codes: Iterable[int] = (150,),
    gaze_break_exclusion_enabled: bool = False,
) -> dict[str, Any]:
    """
    Generate the complete legacy-format deep trigger report bundle.

    ``gaze_break_codes`` and ``gaze_break_exclusion_enabled`` come from the
    resolved trial-exclusion policy. When enabled, the report separately marks
    a reported-only event as resolved for trial inclusion only if a matched
    configured gaze-break event follows on the next reported frame in the same
    trial and a matched condition onset anchors that trial. Raw alignment
    statuses are never changed.

    Returns output paths, raw alignment counts, and counts of auto-resolved and
    review-required reported-only rows for embedding in the step report.
    """
    resolved_output_dir = assert_output_path_outside_raw(output_dir, raw_root)
    stem = f"sub{subject_id}_all_blocks_reported_vs_observed"
    alignment_csv = resolved_output_dir / f"{stem}.csv"
    anomalies_csv = resolved_output_dir / f"{stem}_anomalies.csv"
    report_md = resolved_output_dir / f"{stem}.md"
    plots_dir = resolved_output_dir / f"{stem}_plots"
    transition_plots_dir = plots_dir / "transition_pulses"
    for output_path in (
        alignment_csv,
        anomalies_csv,
        report_md,
        plots_dir,
        transition_plots_dir,
    ):
        assert_output_path_outside_raw(output_path, raw_root)

    legacy_rows = build_legacy_alignment_rows(alignment_rows, events)
    _annotate_gaze_break_exclusions(
        legacy_rows,
        gaze_break_codes=gaze_break_codes,
        gaze_break_exclusion_enabled=gaze_break_exclusion_enabled,
    )
    anomalies = build_anomaly_rows(legacy_rows, events)
    _render_anomaly_plots(
        anomalies,
        plots_dir=plots_dir,
        stim_channel=stim_channel,
        half_window_ms=anomaly_window_ms,
    )
    warnings_with_alignment_rows = _add_transition_alignment_references(
        transition_warnings,
        events,
        alignment_rows,
    )
    rendered_warnings = _render_transition_plots(
        warnings_with_alignment_rows,
        plots_dir=transition_plots_dir,
        stim_channel=stim_channel,
        half_window_ms=transition_window_ms,
    )

    alignment_fields = list(legacy_rows[0].keys()) if legacy_rows else [
        "align_row",
        "status",
    ]
    anomaly_fields = [
        "anomaly_idx",
        "anomaly_type",
        "align_row",
        "detail",
        "reported_only_reason",
        "constrained_acquisition_run",
        "stage1_adjudication",
        "adjudication_partner_align_row",
        "adjudication_gap_ms",
        "adjudication_detail",
        "reported_idx",
        "reported_trigger",
        "reported_label",
        "reported_block",
        "reported_run",
        "observed_idx",
        "observed_trigger",
        "observed_source_fif",
        "observed_width_samples",
        "observed_start_s",
        "anchor_source",
        "plot_representative_align_row",
        "plot_path",
        "plot_error",
    ]
    anomaly_csv_rows = [
        {key: anomaly.get(key) for key in anomaly_fields}
        for anomaly in anomalies
    ]
    _write_csv(alignment_csv, legacy_rows, alignment_fields)
    _write_csv(anomalies_csv, anomaly_csv_rows, anomaly_fields)
    _write_markdown_report(
        report_md,
        subject_id=subject_id,
        legacy_rows=legacy_rows,
        anomalies=anomalies,
        transition_warnings=rendered_warnings,
        segments=segments,
        segment_summaries=segment_summaries,
        alignment_csv=alignment_csv,
        anomalies_csv=anomalies_csv,
        plots_dir=plots_dir,
        anomaly_window_ms=anomaly_window_ms,
        transition_window_ms=transition_window_ms,
    )

    status_counts = dict(Counter(row["status"] for row in legacy_rows))
    auto_resolved_gaze_break_count = sum(
        row.get("stage1_adjudication")
        == "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK"
        for row in legacy_rows
    )
    resolved_recording_tail_count = sum(
        row.get("stage1_adjudication") == "RESOLVED_RECORDING_TAIL_UNAVAILABLE"
        for row in legacy_rows
    )
    reported_only_review_required_count = sum(
        row.get("status") == "REPORTED_ONLY"
        and row.get("stage1_adjudication") == "REVIEW_REQUIRED"
        for row in legacy_rows
    )
    plot_errors = [
        item["plot_error"]
        for item in [*anomalies, *rendered_warnings]
        if item.get("plot_error")
    ]
    return {
        "status": "warning" if plot_errors else "ok",
        "output_dir": str(resolved_output_dir),
        "report_md": str(report_md),
        "alignment_csv": str(alignment_csv),
        "anomalies_csv": str(anomalies_csv),
        "plots_dir": str(plots_dir),
        "alignment_status_counts": status_counts,
        "reported_only_resolved_by_gaze_break": auto_resolved_gaze_break_count,
        "reported_only_resolved_recording_tail": resolved_recording_tail_count,
        "reported_only_review_required": reported_only_review_required_count,
        "anomaly_count": len(anomalies),
        "transition_warning_count": len(rendered_warnings),
        "anomaly_plot_count": len(
            {item["plot_path"] for item in anomalies if item["plot_path"]}
        ),
        "transition_plot_count": sum(
            bool(item.get("plot_path")) for item in rendered_warnings
        ),
        "plot_errors": plot_errors,
    }


__all__ = [
    "build_anomaly_rows",
    "build_legacy_alignment_rows",
    "generate_legacy_deep_report",
]
