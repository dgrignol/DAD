# compare_reported_vs_observed_triggers_subject.py

## Purpose
This subject-level version:
1. Loads **all reported CSVs** for one subject in filename order (block, then run).
2. Maps each block to the correct FIF (supports `SUBXX_blockYY.fif` and `SUBXX_blockYY_ZZ.fif`).
3. Reads FIF files **one at a time**, using only the trigger channel.
4. Concatenates reported and observed trigger streams.
5. Runs one global confrontation and outputs one report bundle.

## Memory model
- FIFs are never all loaded together.
- One FIF at a time is opened/read.
- Only compact trigger-window metadata is kept for global alignment.

## Inputs
- Reported CSV filename pattern: `debug_actual_triggers*_subXX_blockYY_runZZ.csv`
- FIF filename pattern: `SUBXX_blockYY.fif` or `SUBXX_blockYY_ZZ.fif`

## Outputs
In `--out-dir`:
- `subXX_all_blocks_reported_vs_observed.csv`
- `subXX_all_blocks_reported_vs_observed_anomalies.csv`
- `subXX_all_blocks_reported_vs_observed.md`
- `subXX_all_blocks_reported_vs_observed_plots/` (anomaly windows)
- `subXX_all_blocks_reported_vs_observed_plots/transition_pulses/` (transition warnings, unless disabled)

## Status colors
Status background colors are used in report badges and plot status boxes:
- Green: `MATCH`
- Yellow: `MATCH_WARNING`, `WARNING`
- Red: `MISMATCH`, `REPORTED_ONLY`, `OBSERVED_ONLY`

`MATCH_WARNING` is used for collapsed short transition-pulse cases in either
order:
- `0 -> reported -> unexpected -> 0`
- `0 -> unexpected -> reported -> 0`
These rows are not listed as mismatch anomalies.

## Example command
```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers_subject.py \
  --subject 1 \
  --csv-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/output_files/sub01 \
  --fif-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/data/sub01_260415 \
  --out-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/reports/sub01 \
  --variant rescueTraject
```

## Scroller example
```bash
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/compare_reported_vs_observed_triggers_subject.py \
  --subject 1 \
  --csv-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/output_files/sub01 \
  --fif-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/data/sub01_260415 \
  --out-dir /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/triggers_check/reports/sub01 \
  --open-scroller \
  --scroller-duration 12 \
  --scroller-step 2
```
Scroller controls:
- Slider: move the visible window start time.
- Keyboard: `left/right` or `a/d` for step navigation, `home/end` for start/end.
- Faded vertical lines indicate reported-trigger anchors.

## Full flags
| Flag | Required | Meaning |
| --- | --- | --- |
| `--subject` | yes | Subject ID (e.g., `1`). |
| `--csv-dir` | yes | Folder with reported CSV files. |
| `--fif-dir` | yes | Folder with subject FIF files. |
| `--out-dir` | yes | Output folder for single subject bundle. |
| `--variant` | no | Filter CSV names by token (e.g., `rescueTraject`). |
| `--trigger-channel` | no | Trigger channel in FIF (default `STI101`). |
| `--min-duration` | no | Min duration (seconds) for step merge logic. |
| `--transition-max-samples` | no | Max transient width (samples) for transition collapse. |
| `--include-width-outliers` | no | Add `WIDTH_OUTLIER` anomalies. |
| `--plot-window-ms` | no | Half-window size (ms) for anomaly plots (default `500`). |
| `--transition-plot-window-ms` | no | Half-window size (ms) for transition-warning plots (default `3`). |
| `--disable-plot-transition-pulses` | no | Skip transition-warning plot generation (warnings remain in report). |
| `--open-scroller` | no | Open interactive trigger-channel scroller after outputs are generated. |
| `--scroller-duration` | no | Visible scroller window duration in seconds (default `10`). |
| `--scroller-start` | no | Initial scroller window start time in seconds (default: first reported anchor). |
| `--scroller-step` | no | Keyboard navigation step in seconds (default: half of duration). |
