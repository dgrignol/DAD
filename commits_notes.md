# 16 Jan 2026 22:34

## Simulation pipeline scripts

### `simulations/build_movdot_simulation_inputs.m`
- Resolve stimulus/output paths relative to the repo, load MovDot_SubXX stimuli, and split trials into deviant vs. non-deviant inputs.
- Save dot1/dot2 path arrays with frame-count validation plus metadata (units, dot colors, condition definition) in the output .mat files.

### `simulations/build_simulation_report_assets.m`
- Configure participant inputs, ensure condition-specific dot paths exist, and load cached results or compute dRSA matrices per condition.
- Generate grouped figures and persist report assets/parameters in `simulations/output/subXX` for the report generator.

### `simulations/generate_PIPELINE_simulation_report.m`
- Load precomputed report assets and assemble a PDF with parameter summaries and figure sections via MATLAB Report Generator.
- Validate required assets and write the final report into `simulations/report/subXX`.

### `simulations/messy_PIPELINE_simulation.m`
- Build condition-specific inputs, load dot paths, and run trial-locked triggered subsampling for dRSA computations.
- Plot dot-path and distance diagnostics, then save dRSA matrices/diagonal summaries and PNG outputs.

## Simulation debug plotting

### `simulations/debug/plot_paths.m`
- Render time-graded scatter plots for dot trajectories with optional alpha/axis controls.
- Auto-sanitize the title into a filename and save PNGs into `simulations/debug`.

### `simulations/debug/plot_position_time_distance.m`
- Compute time-by-time mean distance matrices across trials and visualize them as heatmaps.
- Save diagnostics to `simulations/debug` with sanitized titles for repeatability.

## dRSA functions

### `simulations/functions/dRSA_PCR.m`
- Implement PCR-based dRSA with optional autocorrelation and model-regression steps per timepoint.
- Support extra PCA preprocessing and multiple component-selection modes for large models.

### `simulations/functions/dRSA_average.m`
- Average dRSA diagonals across a symmetric time window defined by `AverageTime` and `fs`, handling out-of-range samples as NaNs.

### `simulations/functions/dRSA_border.m`
- Compute autocorrelation regression borders by running dRSA autocorrelations, averaging diagonals, and thresholding by variance.

### `simulations/functions/dRSA_computeRDM.m`
- Build time-resolved RDMs via pdist (or correlation fast-path), with optional normalization for PCR workflows.

### `simulations/functions/dRSA_concatenate.m`
- Concatenate trial/segment data into a features×time matrix from arrays, cell matrices, or .mat paths.
- Generate and append a "Concatenation mask" that flags boundaries, with diagnostic plotting.

### `simulations/functions/dRSA_coreFunction.m`
- Centralize dRSA execution by filling defaults, validating subsamples, and computing neural/model RDMs as needed.
- Dispatch to correlation or PCR engines, including an autocorrelation-only path when requested.

### `simulations/functions/dRSA_corr.m`
- Compute correlation-based dRSA matrices for each requested model against the neural RDM.

### `simulations/functions/dRSA_fastpdist.m`
- Provide a vectorized correlation-distance implementation for fast RDM construction.

### `simulations/functions/dRSA_random_subsampling.m`
- Generate random non-overlapping subsamples from an availability mask with spacing and repetition checks.
- Return indices plus a visualization mask and emit warnings when repetition rates are high.

### `simulations/functions/dRSA_resample_subsamples.m`
- Resample existing subsample sets with replacement to build iteration tensors at a configurable keep percentage.

### `simulations/functions/dRSA_rescaleRDM.m`
- Rescale RDM values to the [0, 1] interval and mean-center each column.

### `simulations/functions/dRSA_standardizeRDM.m`
- Mean-center and z-score RDM values to normalize variance for downstream analysis.

### `simulations/functions/dRSA_subsampling_diagnostics.m`
- Report mask statistics, iteration overlap (Jaccard), coverage, and start-point entropy with summary plots.

### `simulations/functions/dRSA_trialwise_subsampling.m`
- Sample subsamples trial-by-trial, limiting to one window per trial with repetition control and visualization masks.

### `simulations/functions/dRSA_triggered_subsampling.m`
- Time-lock subsamples to trigger points with pre/post windows, spacing constraints, and repetition tracking.

# 16 Jan 2026 22:25

## Experiment documentation

### `experiment/MoveDot1_experiment_vX.m`
- Indent the top banner line and adjust spacing in the output path comment; no runtime behavior changes.

## Trigger analysis tooling

### `inspect_fif_report.py`
- Update replay expectations so each replay-start (151) opens a replay window until the next block's first expected trigger, and tag gaze breaks (150) plus the missing 81 as GazeBreak in the match column.
- Build replay windows for every 151 using expected block order and an end-of-recording fallback so matching pauses per replay block.
- Handle gaze breaks, replay starts, and response triggers (201) explicitly during alignment so they stay visible without advancing expected matching.
- Add a gaze-break count to replay summaries and add HTML match styling for RESPONSE and GazeBreak rows.

## Reports

### `reports/sub98/sub98_report.html`
- Regenerate the report with updated trigger matching so GazeBreak/REPLAY/RESPONSE labels and replay summary counts are included.
- Refresh MNE section IDs and the footer timestamp from the new output.

# 14 Jan 2026 12:50

## Replay trigger semantics and reporting

### `experiment/MoveDot1_experiment_vX.m`
- Update trigger documentation so replay-start (151) is described as firing once per block before the first replay trial, keeping the 100 ms pause note.
- Move the replay-start guard into the per-block loop so each block can emit one 151 before its first replayed trial.
- Comment out the MEG photodiode flicker block in the video playback section.

### `experiment/trigger_codes.md`
- Clarify that trigger 151 is emitted once per block, immediately before the first replayed trial in that block, with the 100 ms pause preserved.

## Trigger analysis tooling

### `expected_triggers.py`
- Fix run-count detection by using the first axis of `TrialStruct` so 2D run-by-trial arrays no longer inflate the run count.

### `inspect_fif_report.py`
- Add replay-rule documentation and a usage line in the script header, aligning expectations with trigger-code definitions.
- Introduce replay validation checks (150-to-81 violations, 151 ordering/counts) and surface them in a new HTML summary panel.
- Mark triggers between the first 151 and the first expected block-2 trigger as REPLAY so expected matching pauses during replay segments.
- Add a block column to the HTML trigger table and color-code REPLAY matches for clarity.
- Centralize missing-row detection and reuse it for event extraction and trial numbering.

## Reports

### `reports/sub98/sub98_report.html` (new)
- New trigger report for subject 98 with replay-aware matching and summary panels.

# 13 Jan 2026 11:57

## Eye-tracker testing modes and gaze-break handling

### `experiment/MoveDot1_experiment_vX.m`
- Extend script-level documentation with gaze-break (150) details, usage steps,
  and testing mode behavior for ignore vs fake tracking, aligned to
  `experiment/trigger_codes.md`.
- Add `Conf.IgnoreEyeTracker`, `Conf.FakeEyeTracker`, and
  `Conf.fakeGazeBreakRate`, then derive `useEyelink`, `useFakeEyeTracker`, and
  `useGazeMonitoring` so real vs simulated gaze is resolved in one place.
- Gate EyeLink setup, recording, and messaging on `useEyelink` to avoid
  hardware I/O during testing while preserving MEG-trigger flow.
- Schedule per-trial fake gaze breaks with a random onset that respects
  `Conf.fixBreakToleranceFrames`, generate out-of-window samples during the
  break, and feed them into the existing fixation-break logic/output.
- Preserve inline trigger-code annotations for gaze-break (150) and replay
  start (151) so the testing path stays consistent with the experiment
  definitions.

# 13 Jan 2026 11:25

## Experiment runtime trigger behavior

### `experiment/MoveDot1_experiment_vX.m`
- Add script-level documentation with usage, inputs/outputs, and trigger logic assumptions aligned to `experiment/trigger_codes.md`.
- Introduce a run-level guard so the replay-start trigger (151) is emitted once before the first replay trial, then pause 100 ms to keep separation from subsequent trial-start pulses.

## Trigger reference updates

### `experiment/trigger_codes.md`
- Clarify that trigger 151 fires once before the first replayed trial begins and includes a 100 ms pause after the pulse.

# 9 Jan 2026 16:23

## Documentation and workflow guidance

### `AGENTS.md`
- Add a requirement to include extensive script-level documentation with usage
  examples, and tighten wording around section-level comments describing data flow.

### `trigger_pipeline.md`
- Update expected/actual trigger output paths to use `derivatives/triggers/subXX`
  and clarify that the actual-trigger CSV is now unified per subject.
- Add multi-run usage examples for expected-trigger generation and report
  creation, including default plotting behavior.
- Document the detailed matching procedure used by `inspect_fif_report.py`,
  covering overlap resolution, expected alignment, missing rows, and extra match
  columns.

## Trigger analysis tooling

### `expected_triggers.py`
- Expand the module docstring with inputs/outputs, run selection behavior, and
  usage examples aligned to `experiment/MoveDot1_experiment_vX.m` and
  `experiment/trigger_codes.md`.
- Make `--run` optional so omission generates all runs found in `TrialStruct`,
  while re-seeding per run to keep outputs identical to standalone invocations.
- Add `resolve_runs` and update CSV output/logging to include run IDs explicitly.

### `inspect_fif_report.py`
- Add multi-run processing when `--subject` is set and `--fif` is default:
  discover all `subXX_runYY.fif` files, require expected block CSVs per run, and
  generate a single report with per-run sections.
- Replace per-run/per-block CSV outputs with one unified detected-trigger CSV
  containing `trial`, `trial_exp`, `run`, and `block` columns populated from
  in-memory alignment results.
- Introduce helpers for run-specific path resolution, CSV parsing, and block
  assignment; external match CSVs now filter by run and skip missing rows to
  avoid false matches.
- Add trial counters based on trial-start triggers (1-80/102) for both detected
  and expected sequences, and surface them in the HTML tables.
- Update default report naming for multi-run output and restrict `--browse-raw`
  to the first run with a clear log message.

### `plot_trigger_chan.py`
- Expand script documentation with inputs/outputs and usage examples, and switch
  the default missing-trigger source to the unified CSV when `--subject` is set.
- Allow `--plot-missing` to accept CSV or report HTML, with subject/run-aware
  path resolution and robust parsing based on table headers.
- Default to absolute time origin when missing markers are plotted and skip
  overlays cleanly when sources are absent.

## Reports

### `reports/sub04/sub04_report.html` (new)
- New multi-run MNE report for subject 04 with per-run raw, sensor layout, and
  trigger-channel sections plus trigger tables that include trial/trial_exp
  columns.
- Footer indicates creation on 2026-01-09 11:33:03 with MNE-Python 1.10.1.

## Python bytecode artifacts

### `__pycache__/inspect_fif_report.cpython-313.pyc`
- Updated Python 3.13 bytecode cache generated by the trigger-report script; no
  textual diff available.

# 8 Jan 2025 15:45

## Changes made to the experiment code and documentation

### `experiment/MoveDot1_experiment_vX.m`
- Delay the `desiredFrameRate` assignment until stimulus FPS is loaded, then set
  `Conf.refrate` from `Dat.Cfg.fps` and mirror it into `desiredFrameRate`.
- Add a safety guard (for non-practice runs) that checks for existing output
  files (`*_SUBxx_RUNyy` and `*_firstBlock`, both with and without `.mat`) and
  aborts the run with a detailed error if any exist, preventing accidental
  overwrites.
- Precompute catch-phase durations in whole frames (fixation and occlusion
  phases) to avoid repeated floating-point comparisons and to keep frame
  boundaries consistent.
- Replace all catch-phase comparisons with the precomputed frame counts,
  including the "last frame" checks and occlusion boundary logic.
- Add a refresh-rate sanity check immediately after opening the Psychtoolbox
  window: compute the measured Hz from `GetFlipInterval` and warn if it differs
  from the expected rate by more than 1 Hz or 5%.
- Update the EyeLink file preamble text to the new project label
  ("Damiano Attention Dot Project").

### `experiment/lib/Config.m`
- Increase `frameFrequency` from 60 to 120 to reflect the higher display
  refresh rate expected by the experiment configuration.

### `experiment/trigger_codes.md`
- Rewrite the trial-onset mapping to make the block/condition logic explicit:
  blockCondition 1 uses 1-20 (condition 0) or 21-40 (condition 45), and
  blockCondition 2 uses 41-60 or 61-80.
- Clarify catch start/end semantics as "glitch" boundaries rather than generic
  start/end wording.
- Expand response trigger notes: 201 is expected; 202 is flagged as unexpected
  and is not in `TriggerValues`.
- Add notes about how `condition` maps to `directionVariance` and how
  `blockCondition` is sourced from `BlockOrder`.

### `experiment/TODO.md`
- Add a new "SAVING ISSUE" checklist with concrete failure modes and mitigation
  ideas (backup save, saving order, catch output indexing risks, teardown
  ordering, and partial save on ESC).
- Add a task to update scripts for the new data folder layout (notably
  `inspect_fif_report.py` for `data/subXX/MEG/subXX_runYY.fif`).
- Add tasks to build a trigger-fix script for subject 99 simulation issues and
  to improve trigger checking for replays (using 150/151 and summary stats).
- Mark two items in the "FROM DAVIDE" section as completed.

## Trigger analysis tooling

### `trigger_events.py` (new)
- Introduce `find_events_with_overlaps`, a helper that wraps
  `mne.find_events(..., output="step")` and adds non-zero-to-non-zero steps so
  overlapping triggers are not lost.
- Implement step merging, mask handling, and duplicate filtering to keep the
  output compatible with MNE event arrays while preserving overlaps.

### `inspect_fif_report.py`
- Switch event detection from `mne.find_events` to
  `find_events_with_overlaps`, enabling detection of overlapping trigger
  transitions.
- Update default FIF selection: when `--subject` is provided and `--fif` is the
  default, the script now targets
  `data/subXX/MEG/subXX_runYY.fif` (with `--run`).
- Expand expected-trigger loading to auto-discover block files via a glob
  pattern, track per-block ranges, and maintain block-order concatenation.
- Add overlap resolution logic that uses the expected sequence to decide which
  overlapping trigger to close or start, and advance the expected index
  accordingly.
- Add automatic CSV output routing: when using the default output path and
  `--subject` is set, write per-block outputs if expected block files exist;
  otherwise write a single per-run output.
- Add an HTML "Missing triggers" summary that counts and reports missing
  expected triggers (e.g., "101 -> X missing out of Y") beneath the table.

### `plot_trigger_chan.py` (new)
- Provide a CLI tool for interactive, scrollable visualization of the trigger
  channel with a slider and mouse-wheel scrolling.
- Support detached plotting (`--no-block`) by spawning a child Python process.
- Add `--plot-missing` to read a report HTML, infer missing-trigger midpoints,
  and overlay them in red on the plot.
- Allow selecting the time origin (`raw` vs `absolute`) and handle clipping of
  missing points to the data range.

## Reports

### `reports/dRSA_Att_test00_report.html` (new)
- New MNE report for the `dRSA_Att_test00` data with a "Triggers" section,
  including an embedded trigger-channel plot and a trigger table.
- Includes a "Check passed: all triggers allowed" status.
- Footer indicates creation on 2026-01-07 16:06:05 with MNE-Python 1.11.0.

### `reports/sub04_run01_report.html` (new)
- New MNE report for subject 04, run 01 with trigger-channel visualization and
  a trigger table that includes a `match` column against expected triggers.
- The table shows multiple 201 entries with match=0 and several 101 entries
  marked "MISSING".
- Includes a missing-trigger summary line (101 missing out of 20) and a
  "Check passed: all triggers allowed" status.
- Footer indicates creation on 2026-01-08 10:55:22 with MNE-Python 1.10.1.

### `reports/sub99_report.html`
- Update the report title to reference `sub99_run01.fif` instead of
  `raw_sub99.fif`.
- Trigger table now includes two match columns: `match` and `match debug1`,
  consistent with comparing against multiple actual-trigger sources.
- Footer timestamp updated to 2026-01-07 16:06:40 (MNE-Python 1.11.0).

## Python bytecode artifacts

### `__pycache__/inspect_fif_report.cpython-313.pyc`
### `__pycache__/trigger_events.cpython-313.pyc`
- Binary `.pyc` files generated by Python 3.13 when running the scripts above.
  These are compiled caches with no textual diff.
# 9 Jan 2026 18:28

## Experiment planning updates

### `experiment/TODO.md`
- Replace the initial "Questions" header with the active checklist and add a
  task to verify non-deviant paths stay within constraints.
- Move the practice-response timing, drift/saccade handling, and false-positive
  feedback items into DONE so the active list reflects remaining work.

## Stimulus generation logic

### `experiment/stimuli_generation_v5.m`
- Add script-level documentation covering purpose, usage, inputs/outputs, and
  assumptions, including the baseline no-deviant path behavior.
- Track a baseline no-deviant direction/position alongside the real path and
  reset it when starting a new path segment.
- Update per-frame angle and position updates to keep baseline angles free of
  deviants for likelihood trials, then compute baseline positions in parallel.
- Expand boundary checks to fail paths if either real or baseline trajectories
  exit bounds, and advance the baseline position after successful frames.

## Commit notes

### `commits_notes.md`
- Add a new entry documenting the staged experiment checklist and
  stimuli-generation changes.

# 12 Jan 2026 19:34

## Catch-trial generation

### `experiment/CreateInputFiles_v10.m`
- Add script-level documentation with usage, inputs/outputs, and catch-jitter
  assumptions for the input/output file workflow.
- Precompute the dot-rect center and max distance to support distance-based
  jitter scaling with an `eps` guard against divide-by-zero.
- Introduce configurable jitter parameters (`Catch.JitterBase`,
  `Catch.JitterSlope`, `Catch.JitterScaling`) to control occlusion catch
  amplitude and scaling mode (`off`, `linear`, `log`, `exp`, with `on` as
  linear).
- Replace the fixed alternating y-offset with per-dot, perpendicular jitter
  driven by `AngleDirection + 90`, alternating sign each frame and scaling by
  distance from the rect center before boundary validation.

## Experiment runtime defaults

### `experiment/MoveDot1_experiment_vX.m`
- Flip the default run-mode flags to a practice configuration: disable MEG,
  eye-tracking, and debug, while enabling practice mode.

### `experiment/MoveDot1_experiment_vX.m~` (new)
- Add a full snapshot of the experiment script (editor backup) mirroring the
  current defaults and logic.

## Stimulus configuration and generation

### `experiment/lib/Config.m`
- Retune curvature settings for 120 Hz by halving `curvFactor` to 1.05 and
  randomizing the curviness valence per path (`isCurvValenceRand = true`).

### `experiment/stimuli_generation_v5.m`
- Expand assumptions to note per-dot curviness valence randomization when
  enabled in `Config`.
- Update per-path curvyness updates to compute signed factors independently for
  each dot, preserving per-dot randomness across paths.

# 13 Jan 2026 10:33

## Fixation monitoring

### `experiment/MoveDot1_experiment_vX.m`
- Reduce the fixation window radius from 10 to 7 degrees by updating
  `Conf.fixWindowDeg`, tightening the eye-tracking break threshold.

## Repository cleanup

### `experiment/MoveDot1_experiment_vX.m~` (deleted)
- Remove the editor backup copy of the experiment script to avoid tracking a
  redundant snapshot.

## Commit notes

### `commits_notes.md`
- Add a new entry covering the fixation window update and backup-file removal.
# 13 Jan 2026 19:45

## MEG video sync gating

### `experiment/MoveDot1_experiment_vX.m`
- Add a per-frame `needsVideoSync` flag that is set only when a trigger pulse is queued, tying `RegWrVideoSync` to actual trigger activity.
- Gate all `Datapixx('RegWrVideoSync')` calls on `Conf.MEG && needsVideoSync` to avoid extra sync writes on frames without scheduled triggers.
- Clarify the trigger scheduling data flow with inline comments so MEG sync behavior stays aligned with trigger emission.
