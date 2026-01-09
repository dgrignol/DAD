#!/usr/bin/env python3
"""
Plot the trigger/stim channel in an interactive, scrollable window.

This script loads a FIF file, extracts a trigger channel, and opens a
scrollable/zoomable matplotlib view of the signal. When a detected-trigger CSV
is available, it overlays red markers for missing expected triggers (rows with
NaN start/end), using midpoints between neighboring events. With --subject set,
the default CSV is derivatives/triggers/subXX/actual_triggers_subXX.csv; HTML
report parsing is available as a fallback or via --plot-missing.

Inputs
- FIF file: data/subXX/MEG/subXX_runYY.fif (default when --subject is set)
- Detected-trigger CSV: derivatives/triggers/subXX/actual_triggers_subXX.csv
  (default when --subject is set)
- Report HTML (optional override): reports/subXX/subXX_report.html

Outputs
- Interactive plot window; no files are written.

Usage
  python plot_trigger_chan.py --subject 4 --run 2
  python plot_trigger_chan.py --subject 4 --run 2 --plot-missing derivatives/triggers/sub04/actual_triggers_sub04.csv
  python plot_trigger_chan.py --subject 4 --run 2 --plot-missing reports/sub04/sub04_report.html
  python plot_trigger_chan.py --fif data/sub04/MEG/sub04_run02.fif --trigger-channel STI101
  python plot_trigger_chan.py --subject 4 --run 2 --duration 30 --start 120
"""

import argparse
import csv
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import Slider
import mne


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Open a scrollable, zoomable plot focused on the trigger channel."
    )
    parser.add_argument(
        "--fif",
        type=Path,
        default=Path("data/dRSA_Att_test00.fif"),
        help=(
            "Path to the raw FIF file to inspect "
            "(overridden by --subject/--run when using the default)."
        ),
    )
    parser.add_argument(
        "--subject",
        type=int,
        default=None,
        help="Subject number used to build the default FIF path (e.g., 4).",
    )
    parser.add_argument(
        "--run",
        type=int,
        default=1,
        help="Run number used to build the default FIF path (default: 1).",
    )
    parser.add_argument(
        "--trigger-channel",
        type=str,
        default="STI101",
        help="Channel name carrying triggers to display.",
    )
    parser.add_argument(
        "--preload",
        action="store_true",
        help="Preload the raw data into memory before plotting.",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=10.0,
        help="Window duration in seconds for the plot.",
    )
    parser.add_argument(
        "--start",
        type=float,
        default=0.0,
        help="Start time (s) for the plot.",
    )
    parser.add_argument(
        "--no-block",
        dest="block",
        action="store_false",
        help="Run the plot in a detached child process and return immediately.",
    )
    parser.add_argument(
        "--time-origin",
        choices=("raw", "absolute"),
        default=None,
        help=(
            "Time origin for the x-axis: 'raw' starts at 0 s (first sample), "
            "'absolute' uses raw.first_samp/sfreq. Defaults to 'absolute' when "
            "missing triggers are plotted, otherwise 'raw'."
        ),
    )
    parser.add_argument(
        "--plot-missing",
        type=Path,
        default=None,
        help=(
            "Override the missing-trigger source (CSV or report HTML). When "
            "omitted and --subject is set, defaults to "
            "derivatives/triggers/subXX/actual_triggers_subXX.csv, falling "
            "back to reports/subXX/subXX_report.html if needed."
        ),
    )
    parser.set_defaults(block=True)
    return parser


def _spawn_detached_child() -> None:
    # Relaunch this script as a detached process so the caller returns immediately.
    child_env = os.environ.copy()
    child_env["PLOT_TRIGGER_CHAN_CHILD"] = "1"
    child_args = [sys.executable, str(Path(__file__).resolve())]
    for arg in sys.argv[1:]:
        if arg == "--no-block":
            continue
        child_args.append(arg)
    subprocess.Popen(child_args, env=child_env)


def _parse_float(value: str) -> float:
    # Accept "NaN" values from HTML tables while still parsing normal floats.
    value = value.strip()
    if not value or value.lower() == "nan":
        return float("nan")
    return float(value)


def _parse_int(value: str):
    """Parse an integer cell, treating empty/NaN as missing."""
    value = value.strip()
    if not value or value.lower() == "nan":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _find_trigger_table(html, start_idx=0):
    """Return the first trigger-table HTML block after start_idx, if any."""
    for match in re.finditer(
        r"<table>.*?</table>", html[start_idx:], flags=re.DOTALL | re.IGNORECASE
    ):
        table = match.group(0)
        if re.search(r"<th>\s*Trigger\s*</th>", table, flags=re.IGNORECASE):
            return table
    return None


def _extract_trigger_rows(report_path: Path, run_id=None):
    """
    Parse the report HTML and extract trigger-table rows for a run.

    The report is expected to follow the structure produced by
    inspect_fif_report.py, with a "Detected triggers" table under each run.
    Returns (header, rows), where header may be None if not found.
    """
    # Parse the HTML report and locate the trigger table for the selected run.
    html = report_path.read_text(encoding="utf-8", errors="ignore")
    table_html = None
    if run_id is not None:
        section_id = f"Run_{run_id:02d}-Detected_triggers"
        section_match = re.search(
            rf'id="{re.escape(section_id)}"', html, flags=re.IGNORECASE
        )
        if section_match:
            table_html = _find_trigger_table(html, start_idx=section_match.start())
        else:
            # If the report is multi-run but the requested run is missing, fail fast.
            if re.search(r'id="Run_\d{2}-Detected_triggers"', html, flags=re.IGNORECASE):
                raise ValueError(
                    f"Run {run_id:02d} not found in report: {report_path}"
                )
    if table_html is None:
        table_html = _find_trigger_table(html, start_idx=0)
    if table_html is None:
        raise ValueError(f"No trigger table found in report: {report_path}")

    header = None
    rows = []
    for row_html in re.findall(r"<tr>.*?</tr>", table_html, flags=re.DOTALL | re.IGNORECASE):
        if re.search(r"<th", row_html, flags=re.IGNORECASE):
            # Capture the header row so we can map columns by name.
            header_cells = re.findall(
                r"<th[^>]*>(.*?)</th>", row_html, flags=re.DOTALL | re.IGNORECASE
            )
            if header_cells:
                header = [
                    re.sub(r"<[^>]+>", "", cell).strip() for cell in header_cells
                ]
            continue
        cells = re.findall(r"<td[^>]*>(.*?)</td>", row_html, flags=re.DOTALL | re.IGNORECASE)
        if not cells:
            continue
        cells = [re.sub(r"<[^>]+>", "", cell).strip() for cell in cells]
        rows.append(cells)
    return header, rows


def _extract_trigger_records_from_report(report_path: Path, run_id):
    """
    Extract (trigger, start, end) records from the report trigger table.

    Column indices are resolved from the header when available; otherwise, the
    legacy order (trigger/start/end) is used as a fallback.
    """
    header, rows = _extract_trigger_rows(report_path, run_id)
    trig_idx = start_idx = end_idx = None
    if header:
        # Resolve column indices by name so added columns (e.g., trial counters) do not shift.
        header_map = {
            re.sub(r"\s+", " ", name.strip().lower()): idx
            for idx, name in enumerate(header)
        }

        def pick_index(*keys):
            for key in keys:
                if key in header_map:
                    return header_map[key]
            return None

        trig_idx = pick_index("trigger")
        start_idx = pick_index("start (s)", "start", "start_s")
        end_idx = pick_index("end (s)", "end", "end_s")
    if trig_idx is None or start_idx is None or end_idx is None:
        # Fall back to the legacy column order when no header is available.
        trig_idx, start_idx, end_idx = 0, 1, 2
    records = []
    for cells in rows:
        if len(cells) <= max(trig_idx, start_idx, end_idx):
            continue
        trig_val = _parse_int(cells[trig_idx])
        if trig_val is None:
            continue
        try:
            start = _parse_float(cells[start_idx])
            end = _parse_float(cells[end_idx])
        except ValueError:
            continue
        records.append((trig_val, start, end))
    return records


def _missing_triggers_from_records(records):
    """Compute missing-trigger markers from a sequence of trigger records."""
    missing = []
    for idx, (trig_val, start, end) in enumerate(records):
        if not (math.isnan(start) or math.isnan(end)):
            continue
        # Use the midpoint between the nearest known neighbors as a marker.
        prev_end = None
        for j in range(idx - 1, -1, -1):
            prev_end = records[j][2]
            if not math.isnan(prev_end):
                break
        if prev_end is None or math.isnan(prev_end):
            continue
        next_start = None
        for j in range(idx + 1, len(records)):
            next_start = records[j][1]
            if not math.isnan(next_start):
                break
        if next_start is None or math.isnan(next_start):
            continue
        midpoint = (prev_end + next_start) / 2.0
        missing.append((midpoint, trig_val))
    return missing


def _resolve_report_path(report_path: Path, subject=None) -> Path:
    """Resolve a report path using reports/subXX fallbacks when possible."""
    if report_path.exists():
        return report_path
    candidates = []
    if subject is not None:
        candidates.append(Path("reports") / f"sub{subject:02d}" / report_path.name)
    match = re.match(r"sub(\d{2})", report_path.stem, flags=re.IGNORECASE)
    if match:
        sub_id = match.group(1)
        if report_path.parent == Path("."):
            candidates.append(Path("reports") / f"sub{sub_id}" / report_path.name)
        if report_path.parent == Path("reports"):
            candidates.append(report_path.parent / f"sub{sub_id}" / report_path.name)
    if not candidates:
        return report_path
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return report_path


def _resolve_missing_csv_path(csv_path: Path, subject=None) -> Path:
    """Resolve a CSV path using derivatives/triggers/subXX fallbacks."""
    if csv_path.exists():
        return csv_path
    candidates = []
    if subject is not None:
        candidates.append(
            Path("derivatives/triggers") / f"sub{subject:02d}" / csv_path.name
        )
    match = re.match(r"actual_triggers_sub(\d{2})", csv_path.stem, flags=re.IGNORECASE)
    if match:
        sub_id = match.group(1)
        if csv_path.parent == Path("."):
            candidates.append(
                Path("derivatives/triggers") / f"sub{sub_id}" / csv_path.name
            )
        if csv_path.parent == Path("derivatives/triggers"):
            candidates.append(csv_path.parent / f"sub{sub_id}" / csv_path.name)
    if not candidates:
        return csv_path
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return csv_path


def _extract_trigger_records_from_csv(csv_path: Path, run_id):
    """Extract (trigger, start, end) records from a detected-trigger CSV."""
    records = []
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"No CSV header found in: {csv_path}")
        fieldnames = [name.strip() for name in reader.fieldnames if name]
        field_map = {name.lower(): name for name in fieldnames}
        trig_key = field_map.get("trigger")
        start_key = field_map.get("start_s") or field_map.get("start")
        end_key = field_map.get("end_s") or field_map.get("end")
        run_key = field_map.get("run")
        if trig_key is None or start_key is None or end_key is None:
            raise ValueError(f"CSV missing trigger/start/end columns: {csv_path}")

        for row in reader:
            if run_id is not None and run_key is not None:
                run_val = _parse_int(row.get(run_key, ""))
                if run_val is None or run_val != run_id:
                    continue
            trig_val = _parse_int(row.get(trig_key, ""))
            if trig_val is None:
                continue
            try:
                start = _parse_float(row.get(start_key, ""))
                end = _parse_float(row.get(end_key, ""))
            except ValueError:
                continue
            records.append((trig_val, start, end))
    return records


def _select_report_path(report_arg, subject, run_id):
    """
    Decide which report HTML to use for missing-trigger markers.

    Returns the resolved path (which may not exist yet) or None if no
    subject/report context is available.
    """
    if report_arg is not None:
        return _resolve_report_path(report_arg.expanduser(), subject)
    if subject is None:
        return None
    report_dir = Path("reports") / f"sub{subject:02d}"
    subject_report = report_dir / f"sub{subject:02d}_report.html"
    if subject_report.exists():
        return subject_report
    if run_id is not None:
        run_report = report_dir / f"sub{subject:02d}_run{run_id:02d}_report.html"
        if run_report.exists():
            return run_report
    return subject_report


def _missing_triggers_from_report(report_path: Path, run_id):
    """
    Identify missing expected triggers from the report table.

    Missing triggers are rows with NaN start/end values; run_id selects the
    run-specific table when the report contains multiple runs.
    """
    records = _extract_trigger_records_from_report(report_path, run_id)
    return _missing_triggers_from_records(records)


def _missing_triggers_from_csv(csv_path: Path, run_id):
    """Identify missing expected triggers from a detected-trigger CSV."""
    records = _extract_trigger_records_from_csv(csv_path, run_id)
    return _missing_triggers_from_records(records)


def _guess_missing_source_type(path: Path) -> str:
    """Return 'csv' or 'html' based on the file suffix (defaults to csv)."""
    suffix = path.suffix.lower()
    if suffix == ".html":
        return "html"
    if suffix == ".csv":
        return "csv"
    return "csv"


def _resolve_missing_source_path(path: Path, subject):
    """Resolve a missing-trigger source path and classify it as csv/html."""
    source_type = _guess_missing_source_type(path)
    if source_type == "html":
        resolved = _resolve_report_path(path, subject)
    else:
        resolved = _resolve_missing_csv_path(path, subject)
    return resolved, source_type


def _default_missing_source(subject, run_id):
    """Select the default missing-trigger source for a subject/run."""
    if subject is None:
        return None, None
    csv_path = (
        Path("derivatives/triggers")
        / f"sub{subject:02d}"
        / f"actual_triggers_sub{subject:02d}.csv"
    )
    if csv_path.exists():
        return csv_path, "csv"
    report_path = _select_report_path(None, subject, run_id)
    if report_path is not None and report_path.exists():
        return report_path, "html"
    return csv_path, "csv"


def main() -> None:
    # Parse CLI arguments and handle detached plotting mode if requested.
    args = build_parser().parse_args()
    if not args.block and os.environ.get("PLOT_TRIGGER_CHAN_CHILD") != "1":
        _spawn_detached_child()
        return

    # Resolve the raw data path (subject/run can override the default).
    fif_path = args.fif.expanduser()
    default_fif = Path("data/dRSA_Att_test00.fif")
    if args.subject is not None and fif_path == default_fif:
        fif_path = Path(
            f"data/sub{args.subject:02d}/MEG/sub{args.subject:02d}_run{args.run:02d}.fif"
        )
    if not fif_path.exists():
        raise FileNotFoundError(f"FIF file not found: {fif_path}")

    # Load the raw file and pull out the requested trigger channel.
    raw = mne.io.read_raw_fif(fif_path, preload=args.preload)

    if args.trigger_channel not in raw.ch_names:
        raise ValueError(
            f"Channel {args.trigger_channel!r} not found in data. "
            f"Available channels: {', '.join(raw.ch_names)}"
        )

    sfreq = float(raw.info["sfreq"])
    data = raw.get_data(picks=[args.trigger_channel])[0]
    if data.size == 0:
        raise ValueError("No samples found in the trigger channel.")

    # Resolve the missing-trigger source (CSV by default, HTML as fallback).
    missing_required = args.plot_missing is not None
    if args.plot_missing is not None:
        missing_path, missing_kind = _resolve_missing_source_path(
            args.plot_missing.expanduser(), args.subject
        )
    else:
        missing_path, missing_kind = _default_missing_source(args.subject, args.run)
    if missing_path is not None and not missing_path.exists():
        if missing_required:
            raise FileNotFoundError(f"Missing-trigger source not found: {missing_path}")
        print(
            f"Missing-trigger source not found: {missing_path} "
            "(skipping missing-trigger markers)"
        )
        missing_path = None

    # Optionally extract missing triggers from the selected source (per run).
    missing_points = []
    if missing_path is not None:
        try:
            if missing_kind == "csv":
                missing_points = _missing_triggers_from_csv(missing_path, args.run)
            else:
                missing_points = _missing_triggers_from_report(missing_path, args.run)
        except ValueError as exc:
            if missing_required:
                raise
            print(f"{exc} (skipping missing-trigger markers)")
            missing_path = None
            missing_points = []

    # Determine the x-axis origin. "absolute" aligns with raw.first_samp.
    time_origin = args.time_origin
    if time_origin is None:
        time_origin = "absolute" if missing_path is not None else "raw"
    time_offset = raw.first_samp / sfreq if time_origin == "absolute" else 0.0

    missing_times = np.array([t for t, _ in missing_points], dtype=float)
    missing_values = np.array([v for _, v in missing_points], dtype=float)

    if time_origin == "raw" and missing_times.size:
        missing_times = missing_times - (raw.first_samp / sfreq)

    # Clip missing points to the actual data range.
    data_duration = (data.size - 1) / sfreq
    data_start = time_offset
    data_end = time_offset + data_duration
    if missing_times.size:
        in_range = (missing_times >= data_start) & (missing_times <= data_end)
        missing_times = missing_times[in_range]
        missing_values = missing_values[in_range]

    # Determine the initial window and its bounds.
    window_samples = max(1, int(round(args.duration * sfreq)))
    window_samples = min(window_samples, data.size)
    max_start_idx = max(0, data.size - window_samples)
    max_start = max_start_idx / sfreq + time_offset
    start_idx = int(round((args.start - time_offset) * sfreq))
    start_idx = max(0, min(start_idx, max_start_idx))
    start_time = start_idx / sfreq + time_offset

    times = np.arange(data.size) / sfreq + time_offset
    ymin, ymax = float(np.min(data)), float(np.max(data))
    if missing_values.size:
        ymin = min(ymin, float(np.min(missing_values)))
        ymax = max(ymax, float(np.max(missing_values)))
    y_pad = max(1.0, 0.05 * (ymax - ymin)) if ymax > ymin else 1.0

    # Render the initial window and optional missing-trigger markers.
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    line, = ax.plot(
        times[start_idx : start_idx + window_samples],
        data[start_idx : start_idx + window_samples],
        drawstyle="steps-post",
        linewidth=1.0,
    )
    ax.set(
        xlim=(start_time, start_time + (window_samples / sfreq)),
        ylim=(ymin - y_pad, ymax + y_pad),
        xlabel="Time (s)",
        ylabel="Trigger value",
        title=f"Trigger channel: {args.trigger_channel}",
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, alpha=0.3)
    window_duration = window_samples / sfreq
    missing_scatter = None
    if missing_times.size:
        in_view = (missing_times >= start_time) & (
            missing_times <= start_time + window_duration
        )
        missing_scatter = ax.scatter(
            missing_times[in_view],
            missing_values[in_view],
            color="red",
            marker="x",
            s=36,
            label="Missing triggers",
        )
        ax.legend(loc="upper right")
    ax.set_xlim(start_time, start_time + window_duration)
    ax.set_ylim(ymin - y_pad, ymax + y_pad)

    # Add a slider for horizontal scrolling.
    if max_start_idx > 0:
        slider_ax = fig.add_axes([0.12, 0.07, 0.78, 0.04])
        slider = Slider(
            slider_ax,
            "Start (s)",
            valmin=time_offset,
            valmax=max_start,
            valinit=start_time,
        )

        def update(val):
            # Update the plot window when the slider moves.
            new_start = float(slider.val)
            new_start_idx = int(round((new_start - time_offset) * sfreq))
            new_start_idx = max(0, min(new_start_idx, max_start_idx))
            new_end_idx = new_start_idx + window_samples
            line.set_data(
                times[new_start_idx:new_end_idx],
                data[new_start_idx:new_end_idx],
            )
            ax.set_xlim(
                new_start_idx / sfreq + time_offset,
                (new_start_idx + window_samples) / sfreq + time_offset,
            )
            if missing_scatter is not None:
                view_start = new_start_idx / sfreq + time_offset
                view_end = view_start + window_duration
                in_view = (missing_times >= view_start) & (missing_times <= view_end)
                if np.any(in_view):
                    offsets = np.column_stack(
                        (missing_times[in_view], missing_values[in_view])
                    )
                    missing_scatter.set_offsets(offsets)
                    missing_scatter.set_visible(True)
                else:
                    missing_scatter.set_offsets(np.empty((0, 2)))
                    missing_scatter.set_visible(False)
            fig.canvas.draw_idle()

        def on_scroll(event):
            # Mouse-wheel scroll nudges the window by 10% of its duration.
            if event.inaxes is None:
                return
            step = 0.1 * window_samples / sfreq
            direction = getattr(event, "step", None)
            if direction is None:
                if event.button == "up":
                    direction = 1
                elif event.button == "down":
                    direction = -1
                else:
                    return
            new_start = min(
                max(slider.val + direction * step, slider.valmin), slider.valmax
            )
            if new_start != slider.val:
                slider.set_val(new_start)

        slider.on_changed(update)
        fig.canvas.mpl_connect("scroll_event", on_scroll)
    plt.show(block=args.block)


if __name__ == "__main__":
    main()
