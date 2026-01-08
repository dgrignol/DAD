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
- Output: `derivatives/triggers/expected_triggers_subXX_runYY_blockN.csv`.

2) `inspect_fif_report.py`
- Loads a FIF file with MNE, creates an HTML report, and snapshots the trigger channel.
- Detects trigger steps (including overlaps) and converts them to on/off intervals.
- Optionally aligns actual triggers to expected sequences and flags missing triggers.
- Writes detected triggers to CSV (per block when expected files exist).
- Output: HTML report + `derivatives/triggers/actual_triggers_subXX_runYY[_blockN].csv`.

3) `plot_trigger_chan.py`
- Opens an interactive, scrollable plot of the trigger channel.
- Can overlay missing triggers from the report HTML (`--plot-missing`).
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

2) Detect and report actual triggers:
```
python inspect_fif_report.py --subject 4 --run 1 --trigger-channel STI101
```

3) Plot the trigger channel and optionally mark missing triggers:
```
python plot_trigger_chan.py --subject 4 --run 1 --plot-missing reports/sub04/sub04_run01_report.html
```

## Inputs and Outputs

Inputs:
- FIF data: `data/subXX/MEG/subXX_runYY.fif`
- Stimulus and trial structure:
  - `experiment/input_files/MovDot_SubXX.mat`
  - `experiment/input_files/SubXX_TrialStruct.mat`

Outputs:
- Expected triggers: `derivatives/triggers/expected_triggers_subXX_runYY_blockN.csv`
- Actual triggers: `derivatives/triggers/actual_triggers_subXX_runYY[_blockN].csv`
- Report: `reports/..._report.html`

## Notes

- Trigger mapping (e.g., 1-80, 100-103, 150-151, 201) is defined in
  `experiment/trigger_codes.md` and implemented in `experiment/MoveDot1_experiment_vX.m`.
- Expected triggers include only deterministic events; response/gaze break
  triggers appear only in the actual-trigger extraction.
- Missing triggers are shown in the report table as rows with NaN start/end.
  `plot_trigger_chan.py` places a red marker midway between neighboring events.
