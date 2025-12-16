#!/usr/bin/env python3
"""Generate an MNE HTML report to inspect a FIF recording."""

import math
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import mne


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Load a FIF file with MNE and generate an HTML report."
    )
    parser.add_argument(
        "--fif",
        type=Path,
        default=Path("data/dRSA_Att_test00.fif"),
        help="Path to the raw FIF file to inspect (overridden by --subject).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("reports/dRSA_Att_test00_report.html"),
        help="Where to write the HTML report.",
    )
    parser.add_argument(
        "--no-browser",
        dest="open_browser",
        action="store_false",
        help="Only save the report; do not open it in a browser.",
    )
    parser.add_argument(
        "--preload",
        action="store_true",
        help="Preload the raw data into memory before plotting.",
    )
    parser.add_argument(
        "--browse-raw",
        action="store_true",
        help="Open an interactive Raw.plot() browser for the continuous signal.",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=10.0,
        help="Window duration in seconds when browsing interactively.",
    )
    parser.add_argument(
        "--n-channels",
        type=int,
        default=30,
        help="Number of channels to show at once in the interactive browser.",
    )
    parser.add_argument(
        "--start",
        type=float,
        default=0.0,
        help="Start time (s) for the interactive browser.",
    )
    parser.add_argument(
        "--trigger-channel",
        type=str,
        default="STI101",
        help="Channel name carrying triggers to snapshot.",
    )
    parser.add_argument(
        "--trigger-start",
        type=float,
        default=200.0,
        help="Start time (s) of the trigger snapshot window.",
    )
    parser.add_argument(
        "--trigger-stop",
        type=float,
        default=260.0,
        help="End time (s) of the trigger snapshot window.",
    )
    parser.add_argument(
        "--event-min-duration",
        type=float,
        default=0.0,
        help="Minimum event duration (s) to keep when detecting triggers.",
    )
    parser.add_argument(
        "--trigger-csv",
        type=Path,
        nargs="+",
        default=[Path("derivatives/triggers/actual_triggers_sub99.csv")],
        help=(
            "One or more CSV paths. First is the output for detected triggers. "
            "Additional paths (if any) are existing actual trigger CSVs to compare, "
            "each adding a separate match column."
        ),
    )
    parser.add_argument(
        "--subject",
        type=int,
        default=None,
        help="Subject number used to locate expected trigger files (e.g., 99).",
    )
    parser.add_argument(
        "--run",
        type=int,
        default=1,
        help="Run number for expected triggers (default: 1).",
    )
    parser.set_defaults(open_browser=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()

    fif_path = args.fif.expanduser()
    if args.subject is not None and fif_path == Path("data/dRSA_Att_test00.fif"):
        fif_path = Path(f"data/raw_sub{args.subject:02d}.fif")
    if not fif_path.exists():
        raise FileNotFoundError(f"FIF file not found: {fif_path}")

    out_path = args.out.expanduser()

    raw = mne.io.read_raw_fif(fif_path, preload=args.preload)

    report = mne.Report(title=f"Raw overview: {fif_path.name}")
    report.add_raw(raw=raw, title="Raw data", psd=True, butterfly=False)

    fig_sensors = raw.plot_sensors(kind="topomap", show=False)
    report.add_figure(
        fig=fig_sensors,
        title="Sensor layout",
        section="Sensors",
    )

    if args.trigger_start >= args.trigger_stop:
        raise ValueError("trigger-start must be less than trigger-stop.")
    if args.trigger_channel not in raw.ch_names:
        raise ValueError(
            f"Channel {args.trigger_channel!r} not found in data. "
            f"Available channels: {', '.join(raw.ch_names)}"
        )
    # Snapshot of the trigger/stim channel over a specified window.
    start_idx, stop_idx = raw.time_as_index([args.trigger_start, args.trigger_stop])
    trig_data = raw.get_data(
        picks=[args.trigger_channel],
        start=start_idx,
        stop=stop_idx,
        reject_by_annotation="omit",
    )[0]
    times = raw.times[start_idx:stop_idx]
    fig_trig, ax = plt.subplots(figsize=(10, 3))
    ax.step(times, trig_data, where="post", linewidth=1)
    ax.set(
        xlim=(args.trigger_start, args.trigger_stop),
        xlabel="Time (s)",
        ylabel="Amplitude",
        title=f"{args.trigger_channel} from {args.trigger_start:.1f}s to {args.trigger_stop:.1f}s",
    )
    ax.grid(True, alpha=0.3)
    fig_trig.tight_layout()
    report.add_figure(
        fig=fig_trig,
        title=f"Trigger channel ({args.trigger_channel})",
        section="Triggers",
    )

    # Detect trigger events and summarize in the report.
    min_dur = args.event_min_duration if args.event_min_duration > 0 else 0.0
    events = mne.find_events(
        raw,
        stim_channel=args.trigger_channel,
        output="step",
        min_duration=min_dur,
        shortest_event=1,
        verbose=False,
    )
    if events.size > 0:
        sfreq = raw.info["sfreq"]
        allowed_triggers = set(range(1, 82)) | {100, 101, 102, 103, 150, 151, 201}
        match_triggers = allowed_triggers - {150, 151}
        pending = {}
        rows = []
        for sample, prev, new in events:
            if prev == 0 and new > 0:
                pending[new] = sample  # rising edge
            elif prev > 0 and new == 0 and prev in pending:
                start_sample = pending.pop(prev)
                start = start_sample / sfreq
                end = sample / sfreq
                trig_val = int(prev)
                allowed = 1 if trig_val in allowed_triggers else 0
                rows.append((trig_val, start, end, allowed))

        rows.sort(key=lambda x: x[1])
        all_allowed = all(r[3] == 1 for r in rows) if rows else False

        # Load expected trigger order for comparison (if subject provided).
        expected_seq = []
        if args.subject is not None:
            expected_dir = Path("derivatives/triggers")
            for block_idx in (1, 2):
                exp_path = expected_dir / (
                    f"expected_triggers_sub{args.subject:02d}_run{args.run:02d}_block{block_idx}.csv"
                )
                if not exp_path.exists():
                    continue
                with exp_path.open() as f:
                    next(f, None)  # skip header
                    for line in f:
                        try:
                            trig_val = int(line.split(",")[0].strip())
                        except ValueError:
                            continue
                        if trig_val in match_triggers:
                            expected_seq.append(trig_val)

        def align_records(records):
            """Align actual records to expected_seq; insert missing expected as synthetic rows."""
            combined_rows = []
            matches = []
            expected_indices = []
            exp_idx = 0
            act_idx = 0
            while act_idx < len(records):
                trig, start, end, allowed = records[act_idx]
                if trig not in match_triggers:
                    combined_rows.append((trig, start, end, allowed))
                    matches.append(float("nan"))
                    expected_indices.append(None)
                    act_idx += 1
                    continue

                if exp_idx >= len(expected_seq):
                    combined_rows.append((trig, start, end, allowed))
                    matches.append(float("nan"))
                    expected_indices.append(None)
                    act_idx += 1
                    continue

                try:
                    found_idx = expected_seq.index(trig, exp_idx)
                except ValueError:
                    combined_rows.append((trig, start, end, allowed))
                    matches.append(0)
                    expected_indices.append(None)
                    act_idx += 1
                    continue

                # Insert any missing expected triggers before the found match.
                for missing_idx in range(exp_idx, found_idx):
                    missing_trig = expected_seq[missing_idx]
                    combined_rows.append((missing_trig, float("nan"), float("nan"), 1))
                    matches.append("MISSING")
                    expected_indices.append(missing_idx)

                # Add the aligned actual trigger.
                combined_rows.append((trig, start, end, allowed))
                matches.append(1)
                expected_indices.append(found_idx)
                exp_idx = found_idx + 1
                act_idx += 1

            # Append any remaining expected triggers as missing.
            for remaining_idx in range(exp_idx, len(expected_seq)):
                missing_trig = expected_seq[remaining_idx]
                combined_rows.append((missing_trig, float("nan"), float("nan"), 1))
                matches.append("MISSING")
                expected_indices.append(remaining_idx)

            return combined_rows, matches, expected_indices

        combined_rows, matches_primary, expected_idx_list = (
            align_records(rows) if expected_seq else (rows, [float("nan")] * len(rows), [None] * len(rows))
        )

        def load_trigger_list_from_csv(path: Path):
            trigs = []
            try:
                with path.open() as f:
                    next(f, None)
                    for line in f:
                        try:
                            val = int(line.split(",")[0].strip())
                        except ValueError:
                            continue
                        if val in match_triggers:
                            trigs.append(val)
            except FileNotFoundError:
                pass
            return trigs

        def compute_matches_for_rows(actual_trigs):
            matches = []
            act_ptr = 0
            for (_, _, _, _), exp_idx in zip(combined_rows, expected_idx_list):
                if exp_idx is None:
                    matches.append(float("nan"))
                    continue
                if exp_idx >= len(expected_seq):
                    matches.append("MISSING")
                    continue
                try:
                    found_idx = actual_trigs.index(expected_seq[exp_idx], act_ptr)
                except ValueError:
                    matches.append("MISSING")
                    continue
                matches.append(1)
                act_ptr = found_idx + 1
            return matches

        # Prepare match columns: primary detection + any additional files.
        match_columns = [("match", matches_primary)]
        trigger_csv_paths = [p.expanduser() for p in args.trigger_csv]
        extra_paths = trigger_csv_paths[1:]
        for extra_path in extra_paths:
            actual_trigs = load_trigger_list_from_csv(extra_path)
            label = extra_path.stem
            if "_" in label:
                label = label.split("_")[-1]
            match_columns.append((f"match {label}", compute_matches_for_rows(actual_trigs)))

        # Save trigger table to CSV (first path is output).
        csv_path = trigger_csv_paths[0]
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        header_cols = ["trigger", "start_s", "end_s", "allowed"] + [col for col, _ in match_columns]
        with csv_path.open("w", encoding="utf-8") as f:
            f.write(",".join(header_cols) + "\n")
            for row_idx, (trig, start, end, allowed) in enumerate(combined_rows):
                start_str = "NaN" if isinstance(start, float) and math.isnan(start) else f"{start:.6f}"
                end_str = "NaN" if isinstance(end, float) and math.isnan(end) else f"{end:.6f}"
                match_vals = []
                for _, match_list in match_columns:
                    val = match_list[row_idx]
                    if isinstance(val, float) and math.isnan(val):
                        match_vals.append("NaN")
                    else:
                        match_vals.append(str(val))
                f.write(f"{trig},{start_str},{end_str},{allowed}," + ",".join(match_vals) + "\n")
        print(f"Saved detected triggers CSV to {csv_path}")

        def build_html_table(records, match_cols):
            header = (
                "<tr><th>Trigger</th><th>Start (s)</th>"
                "<th>End (s)</th><th>Allowed (0/1)</th>"
            )
            for col_name, _ in match_cols:
                header += f"<th>{col_name}</th>"
            header += "</tr>"

            def match_cell(val):
                if isinstance(val, float) and math.isnan(val):
                    return "<td style='background:#ffd699;text-align:center;'>NaN</td>"
                if val == "MISSING":
                    return "<td style='background:#f8c5c5;text-align:center;'>MISSING</td>"
                color = "#c8f7c5" if val == 1 else "#f8c5c5"
                return f"<td style='background:{color};text-align:center;'>{val}</td>"

            body_parts = []
            for row_idx, (trig, start, end, allowed) in enumerate(records):
                start_str = "NaN" if isinstance(start, float) and math.isnan(start) else f"{start:.3f}"
                end_str = "NaN" if isinstance(end, float) and math.isnan(end) else f"{end:.3f}"
                row_html = (
                    f"<tr><td>{trig}</td><td>{start_str}</td>"
                    f"<td>{end_str}</td>"
                    f"<td style='background:{'#c8f7c5' if allowed else '#f8c5c5'};"
                    f"text-align:center;'>{allowed}</td>"
                )
                for _, match_list in match_cols:
                    row_html += match_cell(match_list[row_idx])
                row_html += "</tr>"
                body_parts.append(row_html)
            body = "".join(body_parts)
            return (
                "<style>table,th,td{border:1px solid #ccc;border-collapse:collapse;}"
                "th,td{padding:4px 8px;font-family:monospace;}</style>"
                f"<table>{header}{body}</table>"
            )

        status_html = (
            "<p style='color:green;font-weight:bold;'>Check passed: all triggers allowed.</p>"
            if all_allowed
            else "<p style='color:red;font-weight:bold;'>Check not passed: some triggers are not allowed.</p>"
            if rows
            else "<p>No complete on/off trigger pairs to check.</p>"
        )

        html = (
            build_html_table(combined_rows, match_columns) + status_html
            if rows
            else f"<p>Events detected on {args.trigger_channel}, but no complete on/off pairs.</p>"
            )

        report.add_html(
            html=html,
            title="Detected triggers",
            section="Triggers",
        )
    else:
        report.add_html(
            html=f"<p>No events detected on {args.trigger_channel}.</p>",
            title="Detected triggers",
            section="Triggers",
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    report.save(out_path, overwrite=True, open_browser=args.open_browser)
    print(f"Saved report to {out_path}")
    if args.open_browser:
        print("Report opened in your default browser.")
    if args.browse_raw:
        print("Launching interactive Raw.plot() browser...")
        raw.plot(
            duration=args.duration,
            n_channels=args.n_channels,
            start=args.start,
            scalings="auto",
            block=True,
            title=f"Raw browser: {fif_path.name}",
        )


if __name__ == "__main__":
    main()
