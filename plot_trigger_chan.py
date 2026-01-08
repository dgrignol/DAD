#!/usr/bin/env python3
"""Open an interactive, scrollable plot for the trigger channel."""

import argparse
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
            "--plot-missing is used, otherwise 'raw'."
        ),
    )
    parser.add_argument(
        "--plot-missing",
        type=Path,
        default=None,
        help="Plot missing triggers in red using the trigger table from a report HTML.",
    )
    parser.set_defaults(block=True)
    return parser


def _spawn_detached_child() -> None:
    child_env = os.environ.copy()
    child_env["PLOT_TRIGGER_CHAN_CHILD"] = "1"
    child_args = [sys.executable, str(Path(__file__).resolve())]
    for arg in sys.argv[1:]:
        if arg == "--no-block":
            continue
        child_args.append(arg)
    subprocess.Popen(child_args, env=child_env)


def _parse_float(value: str) -> float:
    value = value.strip()
    if value.lower() == "nan":
        return float("nan")
    return float(value)


def _extract_trigger_rows(report_path: Path):
    html = report_path.read_text(encoding="utf-8", errors="ignore")
    table_html = None
    for match in re.finditer(r"<table>.*?</table>", html, flags=re.DOTALL | re.IGNORECASE):
        table = match.group(0)
        if re.search(r"<th>\s*Trigger\s*</th>", table, flags=re.IGNORECASE):
            table_html = table
            break
    if table_html is None:
        raise ValueError(f"No trigger table found in report: {report_path}")

    rows = []
    for row_html in re.findall(r"<tr>.*?</tr>", table_html, flags=re.DOTALL | re.IGNORECASE):
        if re.search(r"<th", row_html, flags=re.IGNORECASE):
            continue
        cells = re.findall(r"<td[^>]*>(.*?)</td>", row_html, flags=re.DOTALL | re.IGNORECASE)
        if not cells:
            continue
        cells = [re.sub(r"<[^>]+>", "", cell).strip() for cell in cells]
        rows.append(cells)
    return rows


def _missing_triggers_from_report(report_path: Path):
    rows = _extract_trigger_rows(report_path)
    records = []
    for cells in rows:
        if len(cells) < 3:
            continue
        try:
            trig_val = int(cells[0])
        except ValueError:
            continue
        try:
            start = _parse_float(cells[1])
            end = _parse_float(cells[2])
        except ValueError:
            continue
        records.append((trig_val, start, end))

    missing = []
    for idx, (trig_val, start, end) in enumerate(records):
        if not (math.isnan(start) or math.isnan(end)):
            continue
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


def main() -> None:
    args = build_parser().parse_args()
    if not args.block and os.environ.get("PLOT_TRIGGER_CHAN_CHILD") != "1":
        _spawn_detached_child()
        return

    fif_path = args.fif.expanduser()
    default_fif = Path("data/dRSA_Att_test00.fif")
    if args.subject is not None and fif_path == default_fif:
        fif_path = Path(
            f"data/sub{args.subject:02d}/MEG/sub{args.subject:02d}_run{args.run:02d}.fif"
        )
    if not fif_path.exists():
        raise FileNotFoundError(f"FIF file not found: {fif_path}")

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

    time_origin = args.time_origin
    if time_origin is None:
        time_origin = "absolute" if args.plot_missing is not None else "raw"
    time_offset = raw.first_samp / sfreq if time_origin == "absolute" else 0.0

    missing_points = []
    if args.plot_missing is not None:
        report_path = args.plot_missing.expanduser()
        if not report_path.exists():
            raise FileNotFoundError(f"Report not found: {report_path}")
        missing_points = _missing_triggers_from_report(report_path)

    missing_times = np.array([t for t, _ in missing_points], dtype=float)
    missing_values = np.array([v for _, v in missing_points], dtype=float)

    if time_origin == "raw" and missing_times.size:
        missing_times = missing_times - (raw.first_samp / sfreq)

    data_duration = (data.size - 1) / sfreq
    data_start = time_offset
    data_end = time_offset + data_duration
    if missing_times.size:
        in_range = (missing_times >= data_start) & (missing_times <= data_end)
        missing_times = missing_times[in_range]
        missing_values = missing_values[in_range]

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
