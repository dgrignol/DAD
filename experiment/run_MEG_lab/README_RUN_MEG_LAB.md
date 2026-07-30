# run_MEG_lab Launcher Notes (V30/v23/v21 photodiode)

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

- Stimulus/schedule/runtime stack: `V30 / v23 / v21_photodiode`
- Unique path IDs per condition: `40` (`trialsPerCondition = 40`)
- MEG mode: enabled
- Real EyeLink: enabled
- Fake eye tracker: disabled
- Fixation replay threshold: `fixWindowDeg = 3.0` in `lib/Config_runtime_v21_photodiode.m`
- Fixation replay policy: `replayLimit = 1`
- Start block default: `targetBlock = 1` in `runExpReadable.m`
- Photodiode patch: enabled by default in MEG mode, `50 x 50 px`, at
  `[0 0]` relative to the PTB window top-left.
- Photodiode pulses: four frames for condition/trial onset; two frames for
  marked occlusion, type-1 catch, and catch-question events.

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
