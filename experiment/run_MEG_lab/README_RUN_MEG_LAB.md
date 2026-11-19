# run_MEG_lab Launcher Notes (V30/v23/v25 dot-offset blank)

This folder is a standalone lab copy built from the current `experiment/` sources and launched with:

- `runExpReadable.m` (main launcher)

## Quick start (interactive MATLAB)

1. Open MATLAB.
2. `cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/run_MEG_lab')`
3. `run('runExpReadable.m')`

## Quick start (non-interactive shell)

Preferred for MATLAB R2020a on Apple Silicon:

```bash
arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/run_MEG_lab'); run('runExpReadable.m');"
```

Alternative direct launch (if Intel wrapper is not required on your machine):

```bash
/Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/run_MEG_lab'); run('runExpReadable.m');"
```

## Active lab defaults

- Stimulus/schedule/runtime stack: `V30 / v23 / v25_dotOffsetBlank`
- Unique path IDs per condition: `40` (`trialsPerCondition = 40`)
- MEG mode: enabled
- Real EyeLink: enabled
- Fake eye tracker: disabled
- Fixation replay threshold: `fixWindowDeg = 3.0` in `lib/Config_runtime_v25_dotOffsetBlank.m`
- Fixation replay policy: `replayLimit = 1`
- Start block default: `targetBlock = 1` in `runExpReadable.m`
- Photodiode patch: enabled by default in MEG mode, physically validated at
  `28 x 25 px` and `[193 0]` relative to the PTB window top-left
  (`rect=[193 0 221 25]`).
- Photodiode pulses: four frames for condition/trial onset; two frames for
  marked occlusion, type-1 catch, and catch-question events.
- Trial onset: condition codes `71`, `81`, and `91` use a continuous,
  frame-scheduled ITI-to-onset path. Catch-question code `118` intentionally
  remains on the unscheduled v21 path.
- All-trial guard: the first trial of every run is rendered once into the back
  buffer and erased without a flip. Before every retained `71`/`81`/`91`, v25
  requires at least `3 ms` before the requested flip deadline. A late attempt
  shows one extra low-patch fixation frame and retries, with a maximum of six
  retries (about `50 ms` at `120 Hz`). Exhaustion is saved and the trial
  continues; it does not abort acquisition.
- Warm-up: once per MATLAB experiment launch, code `153` marks 10 cycles of
  four high-patch frames plus one explicit low frame (about `0.42 s` at
  `120 Hz`). No dot or task event is shown, and the warm-up is not repeated
  for every run or block.
- PTB startup: the lab preset keeps `SkipSyncTests=0` and applies
  `SyncTestSettings(0.0005, 50, 0.1, 5)` before opening the window. Applied
  and previous values are saved in `Conf`, and previous values are restored
  at cleanup. Do not add a manual `Screen('CloseAll')` startup call.
- Dot offset: v25 presents a fixation-only, low-patch frame on the VBL
  immediately following the final motion frame, before EyeLink stop/start or
  next-trial preparation. No new trigger is emitted. The console prints a
  `Dot-offset summary`, and the MAT output saves `dot_offset_*` timing fields.
  The full configured ITI still follows preparation, so preparation can add
  fixation time but can no longer leave the endpoint dot visible.
- Refresh validation: v25 measures the Psychtoolbox inter-frame interval and
  reads DataPixx video status; it does not change display configuration.
  Configure and confirm `120 Hz` separately with
  `display_refresh_rate_utility.m` before acquisition. The default runtime
  tolerance is `+/- 0.5 Hz`.

The valid v24 replacement audit measured 433 condition onsets at
`18.76 +/- 0.47 ms`, with 14 rollover frames across seven trials, no guard
exhaustion, and a successful block-10 stop/resume test. V25 preserves that
onset path and changes only the dot-offset display transition plus its saved
diagnostics.

## Short v25 readiness check

Before the participant run, use an available prepared test subject and execute
at least several ordinary trials with the same real EyeLink and MEG settings.
Do not change the real-acquisition overrides in `runExpReadable.m`.

Confirm all of the following:

- The endpoint dot disappears immediately, without a visible stationary hold.
- The following screen is the ordinary fixation-only baseline.
- Each completed run prints `Dot-offset summary` with the median near one
  measured frame (about `8.34 ms` at `120 Hz`) and `positiveMissed=0`.
- Startup reports `rect=[193 0 221 25]`, real EyeLink, `SkipSyncTests=0`, and
  approximately `120 Hz` for PTB and DataPixx.
- The saved output contains finite `dot_offset_stimulus_onset_time` and
  `dot_endpoint_hold_sec` for every completed trial.

## Resume from a later block

Edit `targetBlock` in `runExpReadable.m` before launch.

Example:

- set `targetBlock = 3` to resume from block 3 onward.

## Regeneration behavior (skip-if-exists)

`runExpReadable.m` currently uses:

- `overwriteStimuliInput = false`
- `overwriteScheduleInput = false`

When existing files are found and overwrite is `false`, `runExpReadable.m` skips:

- Step 1 (`stimuli_generation_V30_rescueTraject.m`) if both stimulus files for `targetSub` exist.
- Step 2 (`CreateInputFiles_v23_rescueTraject.m`) if the schedule file for `targetSub` exists.

If files are missing, the corresponding generation step runs normally.
