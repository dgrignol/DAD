# Trigger Pipeline (MoveDot1)

This project uses four Python scripts to generate, detect, and inspect trigger timelines for the MoveDot1 experiment. The logic mirrors the MATLAB implementation in `experiment/MoveDot1_experiment_vX.m` and the trigger definitions in `experiment/trigger_codes.md`.

## Overview

The pipeline answers two questions:
1. What triggers should have occurred, given the subject/run inputs?
2. What triggers actually occurred in the recorded FIF data?

Expected triggers come from the stimulus and trial-structure files. Actual triggers are extracted from the stim channel in the FIF file, then aligned against the expected sequence.

## Script Roles

1) `expected_triggers.py`
- Reconstructs the `condMatrix` and trial ordering to match the MATLAB logic.
- Generates per-block expected triggers (trial onset, catch start/end, trial end).
- Excludes subject-driven events (responses, gaze breaks, replays).
- Output: `derivatives/triggers/subXX/expected_triggers_subXX_runYY_blockN.csv`.

2) `inspect_fif_report.py`
- Loads one or more FIF files with MNE, creates an HTML report, and snapshots the trigger channel.
- Detects trigger steps (including overlaps) and converts them to on/off intervals.
- Aligns actual triggers to expected sequences and flags missing triggers.
- Writes detected triggers to a unified CSV with `trial`, `trial_exp`, `run`, and
  `block` columns. The report
  tables are generated from the in-memory detection results, not from the CSV.
- Output: HTML report + `derivatives/triggers/subXX/actual_triggers_subXX.csv`.

3) `plot_trigger_chan.py`
- Opens an interactive, scrollable plot of the trigger channel.
- Overlays missing triggers from the unified CSV by default; use `--plot-missing`
  to override the CSV/HTML source.
- Useful for quick visual inspection of timing and anomalies.

4) `trigger_events.py`
- Helper used by `inspect_fif_report.py`.
- Extends `mne.find_events(..., output="step")` to include non-zero-to-non-zero overlaps.
- Important when trigger pulses overlap in time.

## Typical Workflow

1) Generate expected triggers (per block):
```
python expected_triggers.py --subject 4 --run 1
```

1b) Generate expected triggers for all runs for a subject:
```
python expected_triggers.py --subject 4
```

2) Detect and report actual triggers (single run):
```
python inspect_fif_report.py --subject 4 --run 1 --trigger-channel STI101
```

3) Detect and report actual triggers (all runs for a subject):
```
python inspect_fif_report.py --subject 4 --trigger-channel STI101
```

4) Plot the trigger channel (missing triggers are shown if the report exists):
```
python plot_trigger_chan.py --subject 4 --run 1
```

## Inputs and Outputs

Inputs:
- FIF data: `data/subXX/MEG/subXX_runYY.fif`
- Stimulus and trial structure:
  - `experiment/input_files/MovDot_SubXX.mat`
  - `experiment/input_files/SubXX_TrialStruct.mat`

Outputs:
- Expected triggers: `derivatives/triggers/subXX/expected_triggers_subXX_runYY_blockN.csv`
- Actual triggers: `derivatives/triggers/subXX/actual_triggers_subXX.csv` (with `trial`,
  `trial_exp`, `run`, and `block` columns)
- Report: `reports/..._report.html` (single run or multi-run per subject)

## Notes

- Trigger mapping (e.g., 1-80, 100-103, 150-151, 201) is defined in
  `experiment/trigger_codes.md` and implemented in `experiment/MoveDot1_experiment_vX.m`.
- Expected triggers include only deterministic events; response/gaze break
  triggers appear only in the actual-trigger extraction.
- Missing triggers are shown in the report table as rows with NaN start/end.
  `plot_trigger_chan.py` places a red marker midway between neighboring events.
- The report tables are built from detected triggers in memory; the unified
  `actual_triggers_subXX.csv` file is an output, not an input. CSVs are only
  read when you pass `--trigger-csv` to add extra match columns.
- When running multi-run mode (`--subject` with default `--fif`), the report
  includes one section per run and the script errors if any run lacks expected
  trigger files in `derivatives/triggers/subXX`.

## Matching Procedure Used by `inspect_fif_report.py`

The trigger table in the report is built by aligning detected trigger windows
against the expected trigger sequence for the same run. The procedure is:

1) Detect step events with overlaps
- Uses `trigger_events.find_events_with_overlaps(..., output="step")` to include
  non-zero-to-non-zero transitions that `mne.find_events` would normally skip.
- Each step event is `(sample, prev_value, new_value)`.

2) Reconstruct trigger on/off windows
- Starts a window when `prev == 0` and `new > 0`.
- Ends a window when `prev > 0` and `new == 0`.
- When `prev > 0` and `new > 0`, an overlap has occurred; the resolver decides
  whether the overlap belongs to the previous trigger, the new trigger, or both.

3) Resolve overlaps using expected order
- Only trigger values in the allowed set are considered for matching; gaze
  breaks (150) and replays (151) are excluded from order matching.
- If either overlapping value is not in the matchable set, both are kept.
- If the next expected trigger matches `prev`, the overlap is assigned to `prev`
  (unless the following expected trigger is `new`, in which case it is split).
- If the next expected trigger matches `new`, the overlap is assigned to `new`.
- Otherwise, the resolver searches the expected sequence ahead and assigns the
  overlap to whichever of `prev`/`new` appears first; if neither appears, it
  defaults to the new trigger.

4) Advance the expected pointer
- After each completed on/off window for a matchable trigger, the pointer moves
  to the next instance of that trigger in the expected sequence.

5) Align actual windows to expected sequence
- Walks the actual windows in order, and for each matchable trigger:
  - Finds the next occurrence in the expected sequence starting at the pointer.
  - Inserts "MISSING" rows for expected triggers skipped between the pointer
    and the found index (these rows have NaN start/end).
  - Marks the actual row as a match and advances the pointer.
- If no expected occurrence is found, the row is marked as an unexpected
  trigger (match value 0).
- After all actual rows, any remaining expected triggers are appended as
  "MISSING" rows.

6) Build additional match columns
- For each extra actual-trigger CSV, the same expected-index mapping is used.
- Rows without an expected index remain NaN.
- Rows with an expected index are marked `1` if the trigger is found in-order
  in the external CSV; otherwise they are marked "MISSING".
