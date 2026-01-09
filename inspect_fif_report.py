#!/usr/bin/env python3
"""
Generate an MNE HTML report to inspect MoveDot1 FIF recordings.

This report includes:
- Raw-data overviews and sensor layouts.
- Trigger-channel snapshots over a requested time window.
- Trigger tables extracted from the stim channel, optionally aligned against
  expected trigger sequences to highlight missing/extra events.

Multi-run behavior:
- When --subject is provided and --fif is left at its default, the script scans
  all run FIF files in data/subXX/MEG and uses expected trigger blocks from
  derivatives/triggers/subXX. Each run becomes a separate report section and
  contributes rows to a unified detected-trigger CSV with trial/trial_exp,
  run, and block columns.

Outputs:
- HTML report (single run or multi-run).
- Unified detected-trigger CSV (stacked runs) with trial/trial_exp, run, and block
  columns: derivatives/triggers/subXX/actual_triggers_subXX.csv

Trial numbering:
- "trial" increments on detected trial-start triggers (1-80 or 102) and skips
  missing expected rows (NaN start/end).
- "trial_exp" increments on expected trial-start triggers (including missing
  expected rows) when expected triggers are available.

Examples:
    # All runs for subject 04, default output path + unified CSV:
    ./inspect_fif_report.py --subject 4

    # Single run via explicit FIF:
    ./inspect_fif_report.py --subject 4 --run 2 --fif data/sub04/MEG/sub04_run02.fif

    # Save report without opening a browser:
    ./inspect_fif_report.py --subject 4 --no-browser --out reports/sub04_report.html
"""

import csv
import math
import argparse
import re
from bisect import bisect_right
from pathlib import Path

import matplotlib.pyplot as plt
import mne
import numpy as np

from trigger_events import find_events_with_overlaps


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Load one or more FIF files with MNE and generate an HTML report."
    )
    parser.add_argument(
        "--fif",
        type=Path,
        default=Path("data/dRSA_Att_test00.fif"),
        help=(
            "Path to the raw FIF file to inspect. If --subject is set and this stays "
            "at the default, all runs in data/subXX/MEG are processed."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help=(
            "Where to write the HTML report. Defaults to "
            "reports/subXX/subXX_report.html when --subject is set (all runs), "
            "reports/subXX/subXX_runYY_report.html for a single run, "
            "otherwise reports/dRSA_Att_test00_report.html."
        ),
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
        default=None,
        help=(
            "One or more CSV paths. First is the unified output for detected triggers "
            "(stacked runs with trial/trial_exp, run, and block columns). Additional "
            "paths (if any) are "
            "existing actual-trigger CSVs to compare, each adding a separate match "
            "column. When omitted, defaults to derivatives/triggers/subXX/"
            "actual_triggers_subXX.csv when --subject is set, otherwise "
            "derivatives/triggers/actual_triggers_sub99.csv."
        ),
    )
    parser.add_argument(
        "--subject",
        type=int,
        default=None,
        help=(
            "Subject number used to locate expected trigger files under "
            "derivatives/triggers/subXX; when paired with the default --fif, "
            "all runs for the subject are processed."
        ),
    )
    parser.add_argument(
        "--run",
        type=int,
        default=1,
        help=(
            "Run number for expected triggers when inspecting a single run "
            "(default: 1; ignored when scanning all runs)."
        ),
    )
    parser.set_defaults(open_browser=True)
    return parser


def _parse_run_id_from_name(path: Path) -> int | None:
    """Extract a run number from a FIF filename like sub04_run02.fif."""
    match = re.search(r"_run(\d+)", path.stem)
    if not match:
        return None
    return int(match.group(1))


def _discover_subject_runs(subject: int) -> list[tuple[int, Path]]:
    """
    Discover all run FIF files for a subject.

    Returns a sorted list of (run_id, fif_path) tuples.
    """
    meg_dir = Path("data") / f"sub{subject:02d}" / "MEG"
    if not meg_dir.exists():
        raise FileNotFoundError(f"MEG directory not found: {meg_dir}")

    fif_paths = sorted(meg_dir.glob(f"sub{subject:02d}_run*.fif*"))
    if not fif_paths:
        raise FileNotFoundError(f"No run FIF files found in: {meg_dir}")

    run_map = {}
    for path in fif_paths:
        run_id = _parse_run_id_from_name(path)
        if run_id is None:
            raise ValueError(f"Could not parse run number from FIF name: {path}")
        if run_id in run_map:
            raise ValueError(
                f"Multiple FIF files found for run {run_id:02d}: {run_map[run_id]} and {path}"
            )
        run_map[run_id] = path

    return sorted(run_map.items())


def _resolve_run_path(path: Path, run_id: int | None) -> Path:
    """Return a run-specific path by replacing or appending _runXX in the name."""
    if run_id is None:
        return path
    name = path.name
    if re.search(r"_run\d+", name):
        name = re.sub(r"_run\d+", f"_run{run_id:02d}", name)
        return path.with_name(name)
    if path.suffix:
        return path.with_name(f"{path.stem}_run{run_id:02d}{path.suffix}")
    return path.with_name(f"{path.name}_run{run_id:02d}")


def _resolve_run_paths(
    paths: list[Path],
    run_id: int | None,
    *,
    skip_first: bool,
) -> list[Path]:
    """Apply run-specific naming to a list of paths when needed."""
    if run_id is None:
        return list(paths)
    resolved = []
    for idx, path in enumerate(paths):
        if skip_first and idx == 0:
            resolved.append(path)
            continue
        resolved.append(_resolve_run_path(path, run_id))
    return resolved


def _format_match_label(path: Path) -> str:
    """Create a stable match-column label from a CSV path."""
    label = path.stem
    label = re.sub(r"_run\d+$", "", label)
    label = re.sub(r"_block\d+$", "", label)
    if "_" in label:
        label = label.split("_")[-1]
    return label


def _parse_csv_int(value: str) -> int | None:
    """Parse an integer CSV cell, treating empty/NaN as missing."""
    if value is None:
        return None
    value = value.strip()
    if not value or value.lower() == "nan":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_csv_float(value: str) -> float | None:
    """Parse a float CSV cell, treating empty/NaN as missing."""
    if value is None:
        return None
    value = value.strip()
    if not value or value.lower() == "nan":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _assign_blocks(
    expected_idx_list: list[int | None],
    block_ranges: list[tuple[int, int, int]],
) -> list[int | None]:
    """
    Map each combined row to a block index using expected-index ranges.

    Rows without an expected index inherit the most recent block to keep
    non-matchable triggers grouped with their surrounding block.
    """
    if not block_ranges:
        return [None] * len(expected_idx_list)
    current_block = block_ranges[0][0]
    block_by_row: list[int | None] = []
    for exp_idx in expected_idx_list:
        if exp_idx is not None:
            for block_idx, start_idx, end_idx in block_ranges:
                if start_idx <= exp_idx <= end_idx:
                    current_block = block_idx
                    break
        block_by_row.append(current_block)
    return block_by_row


def _write_unified_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Write the unified detected-trigger CSV with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def _validate_expected_triggers(
    subject: int,
    run_entries: list[tuple[int, Path]],
    expected_dir: Path,
) -> None:
    """Ensure every run has expected trigger block files available."""
    missing_runs = []
    for run_id, _ in run_entries:
        pattern = (
            f"expected_triggers_sub{subject:02d}_run{run_id:02d}_block*.csv"
        )
        if not list(expected_dir.glob(pattern)):
            missing_runs.append(run_id)
    if missing_runs:
        missing_str = ", ".join(f"{run_id:02d}" for run_id in missing_runs)
        raise FileNotFoundError(
            "Missing expected trigger files for run(s): "
            f"{missing_str} in {expected_dir}"
        )


def _process_trigger_events(
    *,
    args: argparse.Namespace,
    raw: mne.io.BaseRaw,
    report: mne.Report,
    run_id: int | None,
    trigger_base_dir: Path,
    trigger_subject_dir: Path | None,
    use_subject_runs: bool,
    triggers_section: str,
    csv_rows: list[list[str]],
    extra_csv_paths: list[Path],
    extra_match_labels: list[str],
) -> None:
    """
    Detect triggers for one run, append unified CSV rows, and add HTML/figures.

    Expected triggers follow the MoveDot1 mapping documented in
    experiment/trigger_codes.md. Trial counters are derived from trial-start
    codes (1-80 for non-catch trials, 102 for catch trials).
    """
    # Detect trigger events and summarize them in the report.
    min_dur = args.event_min_duration if args.event_min_duration > 0 else 0.0
    events = find_events_with_overlaps(
        raw,
        stim_channel=args.trigger_channel,
        output="step",
        min_duration=min_dur,
        shortest_event=1,
        verbose=False,
    )
    if events.size > 0:
        sfreq = raw.info["sfreq"]
        # Allowed/expected trigger ranges for this experiment.
        allowed_triggers = set(range(1, 82)) | {100, 101, 102, 103, 150, 151, 201}
        match_triggers = allowed_triggers - {150, 151}
        # Trial-start triggers per experiment/trigger_codes.md.
        trial_start_trigs = set(range(1, 81)) | {102}

        # Load expected trigger order for comparison (if subject provided).
        # This creates a single flattened expected sequence plus block index ranges.
        expected_seq = []
        expected_seq_blocks = []
        block_ranges = []
        if args.subject is not None:
            if run_id is None:
                raise ValueError("Run number required to load expected triggers.")
            pattern = (
                f"expected_triggers_sub{args.subject:02d}_run{run_id:02d}_block*.csv"
            )

            def load_block_seq(path: Path):
                seq = []
                with path.open() as f:
                    next(f, None)  # skip header
                    for line in f:
                        try:
                            trig_val = int(line.split(",")[0].strip())
                        except ValueError:
                            continue
                        if trig_val in match_triggers:
                            seq.append(trig_val)
                return seq

            # Load each block's expected triggers and keep block index ordering.
            block_entries = []
            expected_dirs = []
            if trigger_subject_dir is not None:
                expected_dirs.append(trigger_subject_dir)
            if not use_subject_runs:
                expected_dirs.append(trigger_base_dir)

            for expected_dir in expected_dirs:
                block_entries = []
                for exp_path in expected_dir.glob(pattern):
                    match = re.search(r"_block(\d+)$", exp_path.stem)
                    if not match:
                        continue
                    block_idx = int(match.group(1))
                    seq = load_block_seq(exp_path)
                    if seq:
                        block_entries.append((block_idx, seq))
                if block_entries:
                    break

            block_entries.sort(key=lambda x: x[0])
            expected_seq_blocks = block_entries
            for block_idx, seq in expected_seq_blocks:
                start_idx = len(expected_seq)
                expected_seq.extend(seq)
                end_idx = len(expected_seq) - 1
                block_ranges.append((block_idx, start_idx, end_idx))

        # Decide how to resolve overlapping non-zero steps when two triggers overlap.
        # We use expected order (when available) to choose which trigger "owns" the overlap.
        def resolve_overlap(prev_val, new_val, expected_idx):
            if not expected_seq:
                return "new"
            prev_val = int(prev_val)
            new_val = int(new_val)
            if prev_val not in match_triggers or new_val not in match_triggers:
                return "both"
            if expected_idx >= len(expected_seq):
                return "new"
            next_expected = expected_seq[expected_idx]
            if next_expected == prev_val:
                if (
                    expected_idx + 1 < len(expected_seq)
                    and expected_seq[expected_idx + 1] == new_val
                ):
                    return "both"
                return "prev"
            if next_expected == new_val:
                return "new"
            try:
                prev_idx = expected_seq.index(prev_val, expected_idx)
            except ValueError:
                prev_idx = None
            try:
                new_idx = expected_seq.index(new_val, expected_idx)
            except ValueError:
                new_idx = None
            if prev_idx is None and new_idx is None:
                return "new"
            if prev_idx is None:
                return "new"
            if new_idx is None:
                return "prev"
            return "prev" if prev_idx < new_idx else "new"

        # Advance the expected pointer when a trigger matches the expected sequence.
        def advance_expected_idx(trig_val, expected_idx):
            if trig_val not in match_triggers or expected_idx >= len(expected_seq):
                return expected_idx
            try:
                found_idx = expected_seq.index(trig_val, expected_idx)
            except ValueError:
                return expected_idx
            return found_idx + 1

        # Reconstruct trigger on/off windows from the step events.
        active_val = None
        active_start = None
        expected_idx = 0
        rows = []
        for sample, prev, new in events:
            if prev == 0 and new > 0:
                active_val = int(new)
                active_start = sample
            elif prev > 0 and new > 0:
                choice = resolve_overlap(prev, new, expected_idx)
                if choice in ("prev", "both"):
                    if (
                        active_val is not None
                        and active_start is not None
                        and int(prev) == active_val
                    ):
                        start = active_start / sfreq
                        end = sample / sfreq
                        trig_val = int(prev)
                        allowed = 1 if trig_val in allowed_triggers else 0
                        rows.append((trig_val, start, end, allowed))
                        expected_idx = advance_expected_idx(trig_val, expected_idx)
                if choice in ("new", "both"):
                    active_val = int(new)
                    active_start = sample
                else:
                    active_val = None
                    active_start = None
            elif prev > 0 and new == 0:
                if active_val is None or active_start is None or int(prev) != active_val:
                    continue
                start = active_start / sfreq
                end = sample / sfreq
                trig_val = int(prev)
                allowed = 1 if trig_val in allowed_triggers else 0
                rows.append((trig_val, start, end, allowed))
                expected_idx = advance_expected_idx(trig_val, expected_idx)
                active_val = None
                active_start = None

        # Keep rows time-sorted and compute a simple "all triggers allowed" check.
        rows.sort(key=lambda x: x[1])
        all_allowed = all(r[3] == 1 for r in rows) if rows else False

        record_end_time = (raw.first_samp + raw.n_times - 1) / sfreq

        # Compute catch-trial summary stats (hits, misses, false alarms, etc.).
        def compute_catch_stats(records):
            """
            Compute catch-trial accuracy stats for one run.

            Catch windows are 100 -> next trial start (1-80/102). Gaze breaks (150)
            inside a catch window exclude that window and any responses within it.
            """
            catch_trig = 100
            response_trig = 201
            gaze_break_trig = 150
            replay_trig = 151

            events_by_time = []
            for trig, start, _, _ in records:
                if np.isnan(start):
                    continue
                events_by_time.append((float(start), trig))
            events_by_time.sort(key=lambda x: x[0])

            trial_start_times = [t for t, trig in events_by_time if trig in trial_start_trigs]
            catch_times = [t for t, trig in events_by_time if trig == catch_trig]
            response_times = [t for t, trig in events_by_time if trig == response_trig]
            gaze_break_times = [t for t, trig in events_by_time if trig == gaze_break_trig]
            replay_times = [t for t, trig in events_by_time if trig == replay_trig]

            catch_windows = []
            for catch_time in catch_times:
                end_time = record_end_time
                if trial_start_times:
                    idx = bisect_right(trial_start_times, catch_time)
                    if idx < len(trial_start_times):
                        end_time = trial_start_times[idx]
                if end_time > catch_time:
                    catch_windows.append((catch_time, end_time))

            def in_any_window(timestamp, windows):
                for start, end in windows:
                    if start <= timestamp < end:
                        return True
                return False

            valid_catch_windows = []
            discarded_catch_windows = []
            for start, end in catch_windows:
                if any(start <= t < end for t in gaze_break_times):
                    discarded_catch_windows.append((start, end))
                else:
                    valid_catch_windows.append((start, end))

            hits = sum(
                1
                for start, end in valid_catch_windows
                if any(start <= t < end for t in response_times)
            )
            total_catches = len(valid_catch_windows)
            misses = total_catches - hits
            false_alarms = sum(
                1
                for t in response_times
                if not in_any_window(t, valid_catch_windows)
                and not in_any_window(t, discarded_catch_windows)
            )

            non_catch_intervals = []
            for idx in range(len(trial_start_times) - 1):
                start = trial_start_times[idx]
                end = trial_start_times[idx + 1]
                if any(start <= ct < end for ct in catch_times):
                    continue
                non_catch_intervals.append((start, end))

            correct_rejections = sum(
                1
                for start, end in non_catch_intervals
                if not any(start <= t < end for t in response_times)
            )
            accuracy = hits / total_catches if total_catches else float("nan")

            return {
                "total_catches": total_catches,
                "hits": hits,
                "misses": misses,
                "false_alarms": false_alarms,
                "correct_rejections": correct_rejections,
                "accuracy": accuracy,
                "catch_count": len(catch_times),
                "response_count": len(response_times),
                "trial_start_count": len(trial_start_times),
                "discarded_catch_windows": len(discarded_catch_windows),
                "replay_count": len(replay_times),
            }

        # Align actual rows to the expected sequence, inserting placeholders for missing.
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

        def compute_trial_columns(records, expected_idx_list):
            """
            Compute per-row trial counters for actual and expected sequences.

            "trial" increments on detected trial starts (1-80/102) and skips
            missing expected rows; "trial_exp" increments on expected trial
            starts (including missing expected rows).
            """
            trial = 0
            trial_exp = 0
            trial_vals = []
            trial_exp_vals = []
            has_expected = bool(expected_seq)
            for (trig, start, end, _), exp_idx in zip(records, expected_idx_list):
                # Missing expected rows are represented by NaN start/end timestamps.
                missing = (
                    start is None
                    or end is None
                    or (isinstance(start, float) and math.isnan(start))
                    or (isinstance(end, float) and math.isnan(end))
                )
                if trig in trial_start_trigs and not missing:
                    trial += 1
                trial_vals.append(trial)
                if not has_expected:
                    trial_exp_vals.append(None)
                    continue
                if exp_idx is not None and exp_idx < len(expected_seq):
                    if expected_seq[exp_idx] in trial_start_trigs:
                        # Count expected trial starts even when the row is missing.
                        trial_exp += 1
                trial_exp_vals.append(trial_exp)
            return trial_vals, trial_exp_vals

        combined_rows, matches_primary, expected_idx_list = (
            align_records(rows)
            if expected_seq
            else (rows, [float("nan")] * len(rows), [None] * len(rows))
        )
        # Compute per-row trial counters for the CSV/table outputs.
        trial_numbers, trial_exp_numbers = compute_trial_columns(
            combined_rows, expected_idx_list
        )
        catch_stats = compute_catch_stats(rows)

        # Load triggers from existing CSVs to build extra match columns.
        def load_trigger_list_from_csv(path: Path, run_id: int | None):
            trigs = []
            try:
                with path.open(newline="") as f:
                    reader = csv.DictReader(f)
                    if not reader.fieldnames:
                        return trigs
                    fieldnames = [name.strip() for name in reader.fieldnames if name]
                    field_map = {name.lower(): name for name in fieldnames}
                    trig_key = field_map.get("trigger")
                    run_key = field_map.get("run")
                    start_key = field_map.get("start_s")
                    end_key = field_map.get("end_s")
                    if trig_key is None:
                        return trigs

                    for row in reader:
                        # Skip rows from other runs when the CSV includes run labels.
                        if run_id is not None and run_key is not None:
                            run_val = _parse_csv_int(row.get(run_key, ""))
                            if run_val is None or run_val != run_id:
                                continue
                        # Skip missing rows (NaN start/end) to avoid false matches.
                        if start_key and end_key:
                            start_val = _parse_csv_float(row.get(start_key, ""))
                            end_val = _parse_csv_float(row.get(end_key, ""))
                            if start_val is None or end_val is None:
                                continue
                        trig_val = _parse_csv_int(row.get(trig_key, ""))
                        if trig_val is None:
                            continue
                        if trig_val in match_triggers:
                            trigs.append(trig_val)
            except FileNotFoundError:
                pass
            return trigs

        # Build a per-row match vector for an external "actual trigger" CSV.
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
        extra_paths = _resolve_run_paths(extra_csv_paths, run_id, skip_first=False)
        for extra_path, label in zip(extra_paths, extra_match_labels):
            actual_trigs = load_trigger_list_from_csv(extra_path, run_id)
            match_columns.append((f"match {label}", compute_matches_for_rows(actual_trigs)))

        # Assign block indices and append rows to the unified CSV accumulator.
        block_by_row = _assign_blocks(expected_idx_list, block_ranges)
        run_str = str(run_id) if run_id is not None else "NaN"
        for row_idx, (trig, start, end, allowed) in enumerate(combined_rows):
            block_val = block_by_row[row_idx]
            block_str = str(block_val) if block_val is not None else "NaN"
            trial_val = trial_numbers[row_idx]
            trial_exp_val = trial_exp_numbers[row_idx]
            trial_str = str(trial_val)
            trial_exp_str = "NaN" if trial_exp_val is None else str(trial_exp_val)
            start_str = (
                "NaN"
                if isinstance(start, float) and math.isnan(start)
                else f"{start:.6f}"
            )
            end_str = (
                "NaN"
                if isinstance(end, float) and math.isnan(end)
                else f"{end:.6f}"
            )
            match_vals = []
            for _, match_list in match_columns:
                val = match_list[row_idx]
                if isinstance(val, float) and math.isnan(val):
                    match_vals.append("NaN")
                else:
                    match_vals.append(str(val))
            csv_rows.append(
                [
                    trial_str,
                    trial_exp_str,
                    run_str,
                    block_str,
                    str(trig),
                    start_str,
                    end_str,
                    str(allowed),
                ]
                + match_vals
            )

        # Build an HTML table with per-row status and match indicators.
        def build_html_table(records, match_cols, trial_vals, trial_exp_vals):
            header = (
                "<tr><th>Trial</th><th>Trial_exp</th><th>Trigger</th><th>Start (s)</th>"
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
                trial_str = str(trial_vals[row_idx])
                trial_exp_val = trial_exp_vals[row_idx]
                trial_exp_str = "NaN" if trial_exp_val is None else str(trial_exp_val)
                start_str = (
                    "NaN"
                    if isinstance(start, float) and math.isnan(start)
                    else f"{start:.3f}"
                )
                end_str = (
                    "NaN"
                    if isinstance(end, float) and math.isnan(end)
                    else f"{end:.3f}"
                )
                row_html = (
                    f"<tr><td>{trial_str}</td><td>{trial_exp_str}</td><td>{trig}</td>"
                    f"<td>{start_str}</td>"
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

        # Summarize missing expected triggers by trigger value (if any).
        missing_summary_html = ""
        if expected_seq:
            total_counts = {}
            for trig_val in expected_seq:
                total_counts[trig_val] = total_counts.get(trig_val, 0) + 1
            missing_counts = {}
            for (trig, _, _, _), match_val in zip(combined_rows, matches_primary):
                if match_val == "MISSING":
                    missing_counts[trig] = missing_counts.get(trig, 0) + 1
            if missing_counts:
                parts = []
                for trig in sorted(missing_counts):
                    missing = missing_counts[trig]
                    total = total_counts.get(trig, 0)
                    parts.append(f"{trig} -> {missing} missing out of {total}")
                missing_summary_html = (
                    "<p style='margin-top:4px;'>"
                    "<strong>Missing triggers:</strong> "
                    + "; ".join(parts)
                    + "</p>"
                )

        # Create a compact catch-trial summary panel and plot.
        catch_summary_html = ""
        catch_fig = None
        if rows:
            accuracy = catch_stats["accuracy"]
            if math.isnan(accuracy):
                accuracy_pct = 0.0
                accuracy_label = "n/a"
                accuracy_text = "Accuracy n/a (no catch trials detected)"
            else:
                accuracy_pct = accuracy * 100.0
                accuracy_label = f"{accuracy_pct:.1f}%"
                accuracy_text = f"Accuracy {accuracy_pct:.1f} %"

            catch_summary_html = (
                "<div style='margin-top:8px;'>"
                f"<p><strong>{accuracy_text}</strong></p>"
                f"<p>HIT: {catch_stats['hits']} out of {catch_stats['total_catches']}</p>"
                f"<p>MISS: {catch_stats['misses']}</p>"
                f"<p>FALSE ALARM: {catch_stats['false_alarms']}</p>"
                f"<p>CORRECT REJECTION: {catch_stats['correct_rejections']}</p>"
                f"<p>DISCARDED CATCH WINDOWS (gaze break): {catch_stats['discarded_catch_windows']}</p>"
                f"<p>REPLAYS (151): {catch_stats['replay_count']}</p>"
                "<p><strong>LEGEND:</strong><br>"
                "HIT -> responses (201) in between glitch start (100) and the subsequent trial start (1-80/102);<br>"
                "MISS -> periods in between glitch start (100) and the subsequent trial start (1-80/102) without a response (201).<br>"
                "FALSE ALARM -> responses outside a glitch start (100) and the subsequent trial start (1-80/102) period.<br>"
                "CORRECT REJECTION -> no response outside a glitch start (100) and the subsequent trial start (1-80/102) period.<br>"
                "NOTE -> if a gaze break (150) occurs between glitch start (100) and the subsequent trial start "
                "(1-80/102), that catch window is excluded from HIT/MISS/total glitches, and any responses "
                "within it are excluded from FALSE ALARM and CORRECT REJECTION counts."
                "</p>"
                "</div>"
            )

            fig, axes = plt.subplots(1, 2, figsize=(7, 2.5))
            axes[0].barh([0], [accuracy_pct], color="#4caf50")
            axes[0].set_xlim(0, 100)
            axes[0].set_yticks([0])
            axes[0].set_yticklabels(["Accuracy"])
            axes[0].set_xlabel("Percent")
            axes[0].set_title("Catch accuracy")
            label_x = min(accuracy_pct + 2, 98)
            axes[0].text(label_x, 0, accuracy_label, va="center")

            false_alarms = catch_stats["false_alarms"]
            max_fa = max(1, false_alarms)
            axes[1].barh([0], [false_alarms], color="#f4a261")
            axes[1].set_xlim(0, max_fa)
            axes[1].set_yticks([0])
            axes[1].set_yticklabels(["False alarms"])
            axes[1].set_xlabel("Count")
            axes[1].set_title("False alarms")
            fa_label_x = false_alarms + max_fa * 0.02
            if fa_label_x > max_fa:
                fa_label_x = max_fa * 0.98
            axes[1].text(fa_label_x, 0, str(false_alarms), va="center")
            fig.tight_layout()
            catch_fig = fig

        # Overall status summary for allowed triggers.
        status_html = (
            "<p style='color:green;font-weight:bold;'>Check passed: all triggers allowed.</p>"
            if all_allowed
            else "<p style='color:red;font-weight:bold;'>Check not passed: some triggers are not allowed.</p>"
            if rows
            else "<p>No complete on/off trigger pairs to check.</p>"
        )

        html = (
            build_html_table(
                combined_rows, match_columns, trial_numbers, trial_exp_numbers
            )
            + status_html
            + missing_summary_html
            + catch_summary_html
            if rows
            else f"<p>Events detected on {args.trigger_channel}, but no complete on/off pairs.</p>"
        )

        report.add_html(
            html=html,
            title="Detected triggers",
            section=triggers_section,
        )
        if catch_fig is not None:
            report.add_figure(
                fig=catch_fig,
                title="Catch accuracy summary",
                section=triggers_section,
            )
    else:
        report.add_html(
            html=f"<p>No events detected on {args.trigger_channel}.</p>",
            title="Detected triggers",
            section=triggers_section,
        )


def main() -> None:
    """Generate the HTML report and trigger CSV outputs."""
    args = build_parser().parse_args()
    trigger_base_dir = Path("derivatives/triggers")
    trigger_subject_dir = (
        trigger_base_dir / f"sub{args.subject:02d}" if args.subject is not None else None
    )
    # Resolve unified CSV output and any extra CSVs used for match columns.
    default_trigger_csv = (
        trigger_subject_dir / f"actual_triggers_sub{args.subject:02d}.csv"
        if trigger_subject_dir is not None and args.subject is not None
        else Path("derivatives/triggers/actual_triggers_sub99.csv")
    )
    if args.trigger_csv is None:
        output_csv_path = default_trigger_csv
        extra_csv_paths = []
    else:
        output_csv_path = args.trigger_csv[0].expanduser()
        extra_csv_paths = [p.expanduser() for p in args.trigger_csv[1:]]
    extra_match_labels = [_format_match_label(path) for path in extra_csv_paths]
    # Include trial counters first so CSVs align with the report table order.
    csv_header = (
        [
            "trial",
            "trial_exp",
            "run",
            "block",
            "trigger",
            "start_s",
            "end_s",
            "allowed",
            "match",
        ]
        + [f"match {label}" for label in extra_match_labels]
    )
    csv_rows: list[list[str]] = []

    # Resolve which runs to process and validate expected trigger coverage.
    default_fif = Path("data/dRSA_Att_test00.fif")
    fif_path = args.fif.expanduser()
    use_subject_runs = args.subject is not None and fif_path == default_fif
    if use_subject_runs:
        run_entries = _discover_subject_runs(args.subject)
        if trigger_subject_dir is None or not trigger_subject_dir.exists():
            raise FileNotFoundError(
                f"Expected trigger directory not found: {trigger_subject_dir}"
            )
        _validate_expected_triggers(args.subject, run_entries, trigger_subject_dir)
    else:
        if not fif_path.exists():
            raise FileNotFoundError(f"FIF file not found: {fif_path}")
        run_id = _parse_run_id_from_name(fif_path)
        if run_id is None and args.subject is not None:
            run_id = args.run
        run_entries = [(run_id, fif_path)]

    # Choose a default output path that matches single-run or multi-run usage.
    if args.out is None:
        if args.subject is not None:
            if use_subject_runs:
                out_path = (
                    Path("reports")
                    / f"sub{args.subject:02d}"
                    / f"sub{args.subject:02d}_report.html"
                )
            else:
                run_id = run_entries[0][0] if run_entries else args.run
                out_path = (
                    Path("reports")
                    / f"sub{args.subject:02d}"
                    / f"sub{args.subject:02d}_run{run_id:02d}_report.html"
                )
        else:
            out_path = Path("reports/dRSA_Att_test00_report.html")
    else:
        out_path = args.out.expanduser()

    if args.trigger_start >= args.trigger_stop:
        raise ValueError("trigger-start must be less than trigger-stop.")

    # Initialize the report once, then append per-run sections.
    report_title = (
        f"Raw overview: sub{args.subject:02d} ({len(run_entries)} runs)"
        if use_subject_runs
        else f"Raw overview: {run_entries[0][1].name}"
    )
    report = mne.Report(title=report_title)

    browse_raw_target = None

    for run_id, fif_path in run_entries:
        # Load the raw data for this run; preloading is optional.
        raw = mne.io.read_raw_fif(fif_path, preload=args.preload)
        if browse_raw_target is None:
            browse_raw_target = (raw, fif_path)

        # Build per-run titles/sections so multi-run reports stay readable.
        run_label = f"Run {run_id:02d}" if run_id is not None else "Recording"
        title_suffix = f" ({run_label})" if use_subject_runs else ""
        sensors_section = run_label if use_subject_runs else "Sensors"
        triggers_section = run_label if use_subject_runs else "Triggers"

        # Start with raw data + PSD overview for this run.
        report.add_raw(
            raw=raw,
            title=f"Raw data{title_suffix}",
            psd=True,
            butterfly=False,
        )

        # Add a sensor layout snapshot for quick channel topology checks.
        fig_sensors = raw.plot_sensors(kind="topomap", show=False)
        report.add_figure(
            fig=fig_sensors,
            title=f"Sensor layout{title_suffix}",
            section=sensors_section,
        )

        if args.trigger_channel not in raw.ch_names:
            raise ValueError(
                f"Channel {args.trigger_channel!r} not found in data. "
                f"Available channels: {', '.join(raw.ch_names)}"
            )
        # Snapshot of the trigger/stim channel over a specified window.
        # Clamp the window to the available time span to avoid out-of-range indexing.
        first_time = float(raw.times[0])
        last_time = float(raw.times[-1])
        req_start = max(args.trigger_start, first_time)
        req_stop = min(args.trigger_stop, last_time)
        if req_start >= req_stop:
            raise ValueError(
                f"Requested trigger window ({args.trigger_start}-{args.trigger_stop}s) "
                f"is outside data range ({first_time:.3f}-{last_time:.3f}s)."
            )
        # time_as_index returns absolute samples (incl. first_samp) on older MNE.
        # Convert to relative indices for get_data(start/stop).
        start_idx, stop_idx = raw.time_as_index([req_start, req_stop])
        start_idx -= raw.first_samp
        stop_idx -= raw.first_samp
        if start_idx < 0 or stop_idx <= start_idx:
            raise ValueError(
                f"Computed trigger window indices invalid: start={start_idx}, stop={stop_idx}"
            )
        # Use the raw samples directly (do not drop annotated segments) so times
        # stay aligned with detected events.
        trig_data = raw.get_data(
            picks=[args.trigger_channel],
            start=start_idx,
            stop=stop_idx,
            reject_by_annotation="omit",  # match find_events behaviour
        )[0]
        # Build the time axis directly from sample indices so it matches the
        # event detector (which also operates on sample indices).
        # Event times from mne.find_events are based on absolute sample numbers,
        # i.e., they include raw.first_samp. Mirror that here so the plot aligns.
        times = (np.arange(start_idx, stop_idx) + raw.first_samp) / raw.info["sfreq"]
        fig_trig, ax = plt.subplots(figsize=(10, 3))
        ax.step(times, trig_data, where="post", linewidth=1)
        # Y ticks every 5 units for clearer trigger levels.
        ymin, ymax = ax.get_ylim()
        ax.set_yticks(np.arange(0, max(5, ymax + 1), 5))
        ax.set(
            xlim=(args.trigger_start, args.trigger_stop),
            xlabel="Time (s)",
            ylabel="Amplitude",
            title=(
                f"{args.trigger_channel} from {args.trigger_start:.1f}s "
                f"to {args.trigger_stop:.1f}s{title_suffix}"
            ),
        )
        ax.grid(True, alpha=0.3)
        fig_trig.tight_layout()
        report.add_figure(
            fig=fig_trig,
            title=f"Trigger channel ({args.trigger_channel}){title_suffix}",
            section=triggers_section,
        )

        # Detect triggers, append unified CSV rows, and add summary tables/plots.
        _process_trigger_events(
            args=args,
            raw=raw,
            report=report,
            run_id=run_id,
            trigger_base_dir=trigger_base_dir,
            trigger_subject_dir=trigger_subject_dir,
            use_subject_runs=use_subject_runs,
            triggers_section=triggers_section,
            csv_rows=csv_rows,
            extra_csv_paths=extra_csv_paths,
            extra_match_labels=extra_match_labels,
        )

    # Write the unified CSV, then save the report.
    _write_unified_csv(output_csv_path, csv_header, csv_rows)
    print(f"Saved detected triggers CSV to {output_csv_path}")

    # Write the report and optionally open a browser / interactive raw viewer.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    report.save(out_path, overwrite=True, open_browser=args.open_browser)
    print(f"Saved report to {out_path}")
    if args.open_browser:
        print("Report opened in your default browser.")
    if args.browse_raw and browse_raw_target is not None:
        if use_subject_runs and len(run_entries) > 1:
            print(
                "Launching Raw.plot() for the first run only; use --fif to inspect others."
            )
        else:
            print("Launching interactive Raw.plot() browser...")
        browse_raw, browse_path = browse_raw_target
        browse_raw.plot(
            duration=args.duration,
            n_channels=args.n_channels,
            start=args.start,
            scalings="auto",
            block=True,
            title=f"Raw browser: {browse_path.name}",
        )


if __name__ == "__main__":
    main()
