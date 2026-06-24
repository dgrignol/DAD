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
3. Read only short STI101 windows from the source FIF named by each anomaly or
   transition warning.
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
    )

Safety and assumptions
- The canonical step-02 alignment is authoritative; this module never zips
  acquisition FIFs and behavioral blocks or performs a second alignment.
- Source FIF files are opened read-only and never preloaded.
- Every output path is checked with the shared raw-tree guard before writing.
- Plot failures are reported per artifact and do not remove the underlying
  machine-readable anomaly or transition row.
"""

from __future__ import annotations

import csv
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
                "behavioral_context_source": row.get(
                    "behavioral_context_source"
                ),
            }
        )
    return legacy_rows


def build_anomaly_rows(
    legacy_rows: list[dict[str, Any]],
    events: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build anomaly rows and raw-data plot anchors for non-match statuses."""
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
                candidate_index = candidate.get("observed_idx")
                if candidate_index is not None and _safe_text(candidate_index):
                    anchor_event = events_by_index[int(candidate_index)]
                    anchor_source = "next_observed"
                    break
        if anchor_event is None:
            for candidate in reversed(legacy_rows[:row_index]):
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
        anomalies.append(
            {
                "anomaly_idx": len(anomalies) + 1,
                "anomaly_type": row["status"],
                "align_row": row["align_row"],
                "detail": detail,
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
) -> tuple[np.ndarray, np.ndarray, float]:
    """Read one bounded trigger-channel window around an absolute raw sample."""
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
) -> None:
    """Render one non-interactive trigger-channel plot from a source FIF."""
    times_ms, data, sfreq = _read_stim_window(
        fif_path,
        stim_channel,
        center_sample_raw,
        half_window_ms,
    )
    figure = Figure(figsize=(10, 3.4))
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
    """Render one plot per anomaly when a source anchor is available."""
    for anomaly in anomalies:
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
        try:
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
            )
            anomaly["plot_path"] = str(plot_path)
        except Exception as exc:  # pylint: disable=broad-except
            anomaly["plot_error"] = str(exc)


def _render_transition_plots(
    warnings: list[dict[str, Any]],
    *,
    plots_dir: Path,
    stim_channel: str,
    half_window_ms: float,
) -> list[dict[str, Any]]:
    """Render one short raw STI101 plot per collapsed transition pair."""
    rendered: list[dict[str, Any]] = []
    for warning in warnings:
        row = dict(warning)
        filename = (
            f"{int(row['warning_index']):04d}_transition_warning_"
            f"{row['dropped_code']}_to_{row['kept_code']}.png"
        )
        plot_path = plots_dir / filename
        try:
            _save_stim_plot(
                output_path=plot_path,
                fif_path=Path(row["source_fif_path"]),
                stim_channel=stim_channel,
                center_sample_raw=int(row["dropped_start_sample_raw"]),
                half_window_ms=half_window_ms,
                title=(
                    f"Transition {row['order_pattern']} | "
                    f"kept {row['kept_code']}, dropped {row['dropped_code']}"
                ),
                status="WARNING",
                highlighted_span_samples=(
                    int(row["dropped_start_sample_raw"]),
                    int(row["dropped_end_sample_raw"]),
                ),
            )
            row["plot_path"] = str(plot_path)
            row["plot_error"] = ""
        except Exception as exc:  # pylint: disable=broad-except
            row["plot_path"] = ""
            row["plot_error"] = str(exc)
        rendered.append(row)
    return rendered


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
        "",
        "## Summary",
        f"- Total reported events: **{reported_count}**",
        f"- Total observed windows (inside alignment rows): **{observed_count}**",
        f"- Aggregate alignment score (sum across runs): **{score_total}**",
        f"- MATCH: **{status_counts['MATCH']}**",
        f"- MATCH_WARNING: **{status_counts['MATCH_WARNING']}**",
        f"- MISMATCH: **{status_counts['MISMATCH']}**",
        f"- REPORTED_ONLY: **{status_counts['REPORTED_ONLY']}**",
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
            f"- Plot window around each anomaly: +/-{anomaly_window_ms:.1f} ms.",
        ]
    )
    if not anomalies:
        lines.append("No anomalies detected.")
    else:
        lines.extend(
            [
                "",
                "| anomaly_idx | status | align_row | reported_block | reported_run | reported_trigger | observed_trigger | anchor_source | plot_path |",
                "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
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
                        _safe_text(anomaly["reported_block"]),
                        _safe_text(anomaly["reported_run"]),
                        _safe_text(anomaly["reported_trigger"]),
                        _safe_text(anomaly["observed_trigger"]),
                        _safe_text(anomaly["anchor_source"]),
                        _markdown_escape(plot_path),
                    ]
                )
                + " |"
            )
        for anomaly in anomalies:
            lines.extend(
                [
                    "",
                    f"### Anomaly {anomaly['anomaly_idx']}: {anomaly['anomaly_type']} (row {anomaly['align_row']})",
                    f"- Status: {_status_badge(str(anomaly['anomaly_type']))}",
                    f"- Detail: {_markdown_escape(anomaly['detail'])}",
                ]
            )
            if anomaly["plot_path"]:
                lines.append(
                    f"![Anomaly {anomaly['anomaly_idx']}]"
                    f"({_relative_to_report(report_path, Path(anomaly['plot_path']))})"
                )
            elif anomaly["plot_error"]:
                lines.append(f"- Plot error: `{anomaly['plot_error']}`")

    non_matches = [
        row for row in legacy_rows if row["status"] != "MATCH"
    ]
    lines.extend(
        [
            "",
            "## Non-MATCH Rows",
            "| align_row | segment_marker | status | reported_block | reported_run | reported_trigger | reported_label | observed_trigger | observed_start_s | observed_source_fif |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
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
                "| idx | status | dropped | kept | expected | decision_source | unresolved | confidence | raw_sample_frame | source_fif | detail | plot_path |",
                "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
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
                        _status_badge("WARNING"),
                        _safe_text(warning["dropped_code"]),
                        _safe_text(warning["kept_code"]),
                        _safe_text(warning["expected_trigger_code"]),
                        _safe_text(warning["decision_source"]),
                        _safe_text(warning["unresolved_transition"]),
                        _safe_text(warning["transition_resolution_confidence"]),
                        _safe_text(warning["dropped_start_sample_raw"]),
                        _markdown_escape(warning["source_fif"]),
                        _markdown_escape(
                            f"Collapsed {warning['order_pattern']}; kept "
                            f"{warning['kept_code']} and dropped "
                            f"{warning['dropped_code']}."
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
) -> dict[str, Any]:
    """
    Generate the complete legacy-format deep trigger report bundle.

    Returns a compact result suitable for embedding in the step report.
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
    anomalies = build_anomaly_rows(legacy_rows, events)
    _render_anomaly_plots(
        anomalies,
        plots_dir=plots_dir,
        stim_channel=stim_channel,
        half_window_ms=anomaly_window_ms,
    )
    rendered_warnings = _render_transition_plots(
        transition_warnings,
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
        "anomaly_count": len(anomalies),
        "transition_warning_count": len(rendered_warnings),
        "anomaly_plot_count": sum(bool(item["plot_path"]) for item in anomalies),
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
