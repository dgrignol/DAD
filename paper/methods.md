# Methods

## Current implementation update (V30 generation / v23 schedule / v21 photodiode runtime)

The active experiment stack is now:
- generation: `stimuli_generation_V30_rescueTraject.m` with `Config_stimuli_generation_V30_rescueTraject.m` in /lib
- schedule creation: `CreateInputFiles_v23_rescueTraject.m` with `Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m` in /lib
- runtime: `MoveDot1_experiment_occlusion_v21_photodiode.m` with `Config_runtime_v21_photodiode.m` in /lib
- triggers in trigger_codes_occlusion_v10_rescueTraject.md

Behavioral and UI updates applied in the current runtime:
- Catch-question layout is centered more strongly around screen center:
  - question vertical fraction default: `0.42`
  - YES/NO vertical fraction default: `0.58`
  - NO/YES horizontal fractions default: `0.38` / `0.62`
- ITI visual separation supports four runtime-selectable modes (`itiVisualMode`):
  - `screen_gray`: gray screen + gray arena rectangle
  - `fixation_gray`: black/black screen with gray fixation cross (Thaler et al., 2013 )
  - `fixation_short`: black/black ITI with fixation-only color pulse (default onset `100 ms`, duration `300 ms`)
  - `off`: no visual-separation cue
  with configurable gray RGB (`itiGreyColor`, default `[90 90 90]`).
- `screen_gray`/`fixation_gray` visual-separation cue is time-windowed, not full-ITI:
  - `itiVisualSeparationDur = [0.300, 0.700] s` by default
  - if `ITIRangeSec(2) < itiVisualSeparationDur(2) + 0.300`, runtime warns and forces `itiVisualMode='off'`.
- `fixation_short` uses tunable parameters `itiFixationShortDelaySec` (default `0.100 s`) and `itiFixationShortDurationSec` (default `0.300 s`) and does not use the `300 ms` tail guard.
- ITI jitter runtime preset in `Config_runtime_v21_photodiode.m` is `ITIRangeSec = [1.000, 1.500]` (seconds); this remains fully tunable and validated at runtime.
- Catch rates are now `15%` in schedule config:
  - `catchRateType1Run1 = 0.15`
  - `catchRateType2Runs23 = 0.15`
- Optional timed blink breaks were added:
  - toggle: `miniBlinkBreakEnabled`
  - interval: every `miniBlinkBreakEveryNTrials` (default `10`) executed trials
  - duration: `miniBlinkBreakDurationSec` (default `5 s`)
  - text: `"Blink break"` and `"Restart in: <seconds>"`
  - suppression rule: blink break is skipped when another scheduled interruption already occurs (run1->run2 transition message or end-of-block message).
- A `500 ms` black/black fixation cushion is shown at moving-dot to non-dot transitions
  (e.g., before/after catch screens and message screens) to separate non-dot visual events from dot-related responses.
- To avoid double transitions, the ITI visual-separation cue is suppressed whenever the next event is a catch screen, blink-break screen, or run/block message.
- End-of-block message now appends current block catch accuracy as an integer percent (e.g., `Current block catch accuracy: 83%`; `N/A` if no scored catches).
- The dedicated final message screen is removed from runtime flow; the final block uses the same end-of-block message pathway.
- Central fixation rendering in the current runtime uses an ABC/bullseye-style fixation cross (Thaler et al., 2013 ) with a crosshair, replacing the simpler cross used in earlier runtimes.

Photodiode synchronization added in the v21 runtime:
- Trial-condition triggers `71`, `81`, and `91` are emitted only after the
  first stimulus frame is fully drawn and immediately before the patch draw
  and corresponding `Screen('Flip')`.
- The patch is high for four displayed frames at condition/trial onset and
  for two displayed frames at occlusion boundaries (`111`, `112`, `114`,
  `115`), type-1 catch disappear/reappear events (`116`, `117`), and
  catch-question onset (`118`). At `120 Hz`, these defaults correspond to
  approximately `33.3 ms` and `16.7 ms`, respectively.
- The patch is drawn last on every experiment-controlled display flip. It is
  held at its low color outside marked pulses, including an explicit low
  frame after a pulse and during cleanup.
- The runtime preset enables the patch by default when `MEG=true`. Patch
  size (`photodiodePatchSizePx`, default `[50 50]` pixels), top-left-relative
  location (`photodiodePatchPositionPx`, default `[0 0]`), high/low colors,
  and both pulse widths are configurable and validated against the opened
  display bounds.
- Trigger values and the DataPixx output schedule are unchanged from v20;
  v21 changes the trial-condition emission point and adds the visual
  photodiode marker.

## Methods changes

The main methods sections describe the V30 stimulus generation and v23
scheduling shared by the current v21 stack and the v20 stack used for
sub04-sub07. The v21 photodiode runtime is the current launch target, but it
is not assigned retroactively to those participants. Earlier participants
were run with pilot stacks. The table below reports only data actually
saved/presented to participants; generated or prepared extra blocks are not
interpreted as part of the methods.

| Subject | Presented stack | Saved blocks | Source trials completed | Replay presentations | Total presentations | Scored catch questions | Notes |
|---|---|---:|---:|---:|---:|---:|---|
| sub01 | V28 stimuli / v21 schedule / v18 runtime | 11 | 726 | 14 | 740 | 65 |  |
| sub02 | V29 stimuli / v22 schedule / v19 runtime | 10 | 690 | 48 | 738 | 83 |  |
| sub03 | V29 stimuli / v22 schedule / v19 runtime | 10 | 690 | 11 | 701 | 88 |  |
| sub04 | V30 stimuli / v23 schedule / v20 runtime | 5 | 661 | 17 | 678 | 80 | block 02 run 3 stopped by escape key after 18/47 planned trials |
| sub05 | V30 stimuli / v23 schedule / v20 runtime | 5 | 690 | 10 | 700 | 89 |  |
| sub06 | V30 stimuli / v23 schedule / v20 runtime | 5 | 690 | 16 | 706 | 90 |  |
| sub07 | V30 stimuli / v23 schedule / v20 runtime | 5 | 690 | 23 | 713 | 89 |  |

Source trials completed are non-replay scheduled presentations. Total presentations include fixation-break replays. Scored catch questions are the completed catch-question counts saved in the run summaries.

Stimulus-generation changes:
- The unique trajectory pool increased from `20` sequence IDs per condition in sub01-sub03 to `40` sequence IDs per condition in sub04-sub07.
- Condition trigger codes changed with the runtime: v18/v19 used `31` (`always_visible`), `41` (`occluded_nondeviant`), and `51` (`occluded_deviant`); v20 introduced `71`, `81`, and `91`, which v21 retains. Sequence identity triggers were emitted at frame `12`; the usable sequence range was `1-20` for sub01-sub03 and `1-40` for sub04-sub07.
- Baseline and post-deviance curvature windows were stable across the saved rescueTraject inputs: `[-0.8, -0.3755] U [0.3755, 0.8]`, with `randomizeCurvatureOnDeviant = true` and `flipCurvatureOnDeviant = false`.
- The signed deviant-turn window changed across early inputs. The saved presented files show `[-81,-15] U [15,81] deg` for sub01, `[-81,-10] U [10,81] deg` for sub02, and `[-81,0] U [0,81] deg` for sub03-sub07.
- `Cfg.DirChange = 30` was saved in the input files but should not be reported as a fixed deviant turn. The generator sampled actual deviant turns from `deviantSignedTurnWindows`, saved the resulting `xySeqs`, and runtime consumed the precomputed `xy` trajectories rather than recomputing paths from `Cfg.DirChange`.
- Turn-magnitude rescue was enabled for all saved rescueTraject inputs with scale grid `[1.0, 0.85, 0.70, 0.55, 0.40, 0.25, 0.10, 0.0]`; when needed, it reduced the sampled turn magnitude while preserving sign to keep a trajectory feasible.

Schedule, task, and display changes:
- sub01 used v21 scheduling with `20` base trials per run and `10%` type-1/type-2 catch insertion. sub02-sub03 used v22 scheduling with `20` base trials per run and `15%` catches. sub04-sub07 used v23 scheduling with `40` base trials per run and `15%` catches.
- sub01 used a shorter ITI range, `[0.5, 1.0] s`, with black/black fixation and no saved `itiVisualMode`. sub02 used `[1.0, 1.5] s` with `screen_gray` visual separation from `0.3-0.7 s`. sub03-sub07 used `[1.0, 1.5] s` with `fixation_short`: a fixation-only cue after `0.1 s` for `0.3 s`.
- Catch-question layout changed after sub01. sub01 used a wide layout (`top/bottom y = 0.12/0.88`, left/right x = `0.18/0.82`). sub02-sub07 used the centered layout (`top/bottom y = 0.42/0.58`, left/right x = `0.38/0.62`).
- A `0.5 s` black/black message/catch cushion and automatic blink breaks every `10` completed trials for `5 s` were present from sub02 onward and absent in sub01.
- The dedicated final message screen was enabled for sub01 and disabled from sub02 onward.
- The fixation acceptance window was `7 deg` for all participants except sub02, which used `3 deg`.
- Run-family color cueing was enabled for all participants except sub03, where it was disabled.
- The v20 runtime used for sub04-sub07 and the current v21 runtime render an ABC/bullseye-style fixation cross (Thaler et al., 2013 ) with a crosshair; v18/v19 used a simpler fixation cross.

Eye tracking and ocular recordings:
- For sub01, the eye tracker operated online during the experiment, so fixation-break detection and the associated replay logic were active. The saved eye-tracker data were subsequently lost, however, and offline analysis of ocular activity for this participant relied on EOG recordings instead.
- From sub02 onward, EOG was not collected. Offline gaze/ocular analyses for sub02-sub07 therefore relied on the saved eye-tracker data.
- Consequently, the presence of online fixation-break replays for sub01 does not imply that eye-tracker data were available for that participant's offline analysis.

## Stimulus generation (V30 rescueTraject)

Stimuli were generated with `stimuli_generation_V30_rescueTraject.m` and `Config_stimuli_generation_V30_rescueTraject.m`. Each trial contained a single moving dot in a `10 x 10 deg` arena. Dot diameter was `0.51442 deg`, speed was `3.73 deg/s`, frame rate was `120 Hz`, and trial duration was `2.67 s` (`320` frames). Viewing distance at runtime was `1000 mm` by default.

The dot center was bounded by arena margins to keep the full dot inside the `10 x 10 deg` window (`x,y` in approximately `[0.257, 9.743] deg`). With fixation at arena center (`5,5 deg`), this yields an approximate eccentricity span of:
- minimum about `0.557 deg` (from fixation-exclusion constraint in generation)
- maximum about `6.71 deg` (near arena corners)

For each subject, generation was deterministic (`rng(targetSubjectID)`) and produced `40` unique sequence IDs per presented condition. Presented condition labels and codes were:
- `always_visible` (`-1`)
- `occluded_nondeviant` (`0`)
- `occluded_deviant` (`45`, from `directionVariance = [0, 45]`)

For each sequence, nondeviant and deviant trajectories shared the same initial direction and were forced identical through the deviance frame. Baseline curvature was sampled from `[-0.8, -0.3755] U [0.3755, 0.8]` (signed heading increment in `deg/frame`; negative and positive values bend in opposite directions, and excluding near-zero values avoids nearly straight trajectories). Deviant signed turn at onset was sampled from `[-81, 0] U [0, 81]` deg. Post-onset deviant curvature was re-sampled from the same curvature windows (`randomizeCurvatureOnDeviant = true`, `flipCurvatureOnDeviant = false`). Since zero-turn samples are now allowed, generation applies an effective-deviance guard that rejects candidates with zero turn and unchanged post-onset curvature relative to the nondeviant branch. Deviant suffixes were translated at the splice frame to preserve positional continuity.

Generation used fixation-zone collision handling (`fixationCollisionMode = 'move'`): trajectories entering a central exclusion radius of `0.55721 deg` were minimally translated out of the zone if feasible; otherwise they were re-sampled.

## Occlusion model and timing

The current default occlusion model is `pathband` (not twocircle). Metadata also includes an alpha profile fallback, but runtime defaults to `pathband`.

Path-band settings:
- width = `max(1.10 * dotDiameter, dotDiameter + 0.005)` = `0.565862 deg`
- terminal style = `straight`
- straight start backshift scale = `0.0`

Timing constraints:
- deviance frame fixed at `130` (~`1.083 s`)
- first fully occluded frame fixed at `130`
- full occlusion required through frame `190` inclusive
- nominal reappearance search starts at frame `191`

`occlusion_start_frame`, `occlusion_end_frame`, and `occlusion_end_complete_frame` are computed geometrically from trajectory and band overlap, while preserving the fixed full-occlusion anchor at frame `130`.

## Schedule and block-resume execution

Scheduling is created by `CreateInputFiles_v23_rescueTraject.m` with `Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m`.

Latest-stack per-block schedule:
- `3` runs per block
- Run 1: shuffled `always_visible`
- Run 2: shuffled first half of `occluded_nondeviant` + first half of `occluded_deviant`
- Run 3: shuffled remaining halves of both occluded conditions

With `40` base trials per condition, base run lengths are `40/40/40` before catch insertion. Runtime (`MoveDot1_experiment_occlusion_v21_photodiode.m`) supports block resume by starting at user-selected `iBlock`. The actual number of saved blocks and trials per participant is reported in the Methods changes section.

## Catch trials and task

Two catch types are inserted upstream in `CatchPlan`:

Type 1 (run 1 only):
- count per block: `round(0.15 * 40) = 6`
- trial begins as `always_visible`
- disappearance onset sampled from `0.30-1.00 s`
- invisible duration: `0.50 s` (`60` frames at `120 Hz`)
- on reappearance, trajectory either remains on source path or switches to a plausible deviant branch with probability `0.50` (when matched sequence exists)
- expected response: YES if branch changed, NO otherwise

Type 2 (runs 2+3 pooled):
- count per block: `round(0.15 * (40 + 40)) = 12` total across runs 2 and 3
- no added disappearance manipulation (standard occlusion trial)
- expected response: YES for `occluded_deviant`, NO for `occluded_nondeviant`
- catches are inserted into random slots; run-2/run-3 split varies by block

Per complete v23 block, the planned catch counts are:
- type 1: `6`
- type 2: `12`
- total catches: `18`

Actual completed catch-question counts per participant are reported in the Methods changes section.

## Trial procedure and responses

Each trial displayed one trajectory. A central fixation cross (Thaler et al., 2013 ) was drawn during the motion phase, the ITI blank phase, and catch-question screens. No response was required on non-catch trials.

On catch trials, a post-trial question was shown for `4.0 s`:
- schedule default prompt: `"Has the dot changed its course?"`
- runtime preset in `Config_runtime_v21_photodiode.m` overrides this with: `"Did the trajectory change while occluded?"`

Response mapping:
- YES: code `1` (`RightArrow`, `8`, `y`; MEG button value `8`, blue)
- NO: code `2` (`LeftArrow`, `1`, `n`; MEG button value `1`, red)

Catch feedback (`500 ms`):
- correct: green fixation
- incorrect: red fixation
- timeout: `"Too slow"` text

Post-trial ITI is jittered with a uniform draw in a configurable range:
- runtime default: `[1.000, 1.500] s`
- runtime override: `ITIRangeSec = [min, max]`
- implementation detail: all ITIs are precomputed at block start (before trial execution) with deterministic seeding from participant identity plus block/run offsets.
- replay behavior: replayed trials reuse the ITI of their originating source trial via stored source-trial mapping.
- `itiVisualMode` options are `screen_gray`, `fixation_gray`, `fixation_short`, and `off`.
- for `screen_gray`/`fixation_gray`, ITI visual-separation window defaults to `0.300-0.700 s` (`itiVisualSeparationDur`).
- for `fixation_short`, the fixation color pulse starts after `0.100 s` (`itiFixationShortDelaySec`) and lasts `0.300 s` (`itiFixationShortDurationSec`) by default.
- when `ITIRangeSec(2)` is too short to safely include `itiVisualSeparationDur` plus a `300 ms` tail, runtime warns and automatically switches to `itiVisualMode='off'` (applies to `screen_gray`/`fixation_gray`).
- `screen_gray`/`fixation_gray` ITI visual-separation cue is suppressed if the next event is catch question/feedback, blink break, or run/block message, to avoid gray-black-gray double transitions.
- a `500 ms` black/black fixation cushion is inserted before and after non-dot message/catch screens.
- optional mini blink break:
  - a `5 s` countdown screen (`Blink break`, `Restart in: ...`) can be enabled every 10 completed trials
  - blink break is skipped when another scheduled interruption is already imminent (run1->run2 transition screen or end-of-block screen)

The previous grid-mask ITI path (blank + low-luminance grid + blank) is kept in code as commented legacy settings for possible future reuse.

Current runtime display background settings are all black:
- `Conf.background = [0 0 0]`
- `Conf.rectColor = [0 0 0]`
- `Conf.rectBorderColor = [0 0 0]`
so, outside the dot/fixation/mask elements, the screen is black.

## Trial timing

Per-trial timing is:
- Motion epoch: `2670 ms` (`320` frames at `120 Hz`)
- Post-trial ITI: uniform jitter in `[1000, 1500] ms` by default

For catch trials:
- Pre-question black/black cushion: `500 ms`
- Catch question: `0-4000 ms` (ends on response or timeout)
- Catch feedback: `500 ms`
- Post-feedback black/black cushion: `500 ms`
- Then standard jittered ITI from the same run precomputed ITI vector

So:
- Non-catch trial cycle (motion + ITI): `3670-4170 ms` (default ITI range)
- Catch trial cycle (motion + cushions + question + feedback + ITI): `5170-9170 ms` (default ITI range)
- Estimated per-block trial-loop duration (`138` trials: `120` non-catch + `18` catch; excluding breaks/messages/replays): approximately `533-665 s` (`8.9-11.1 min`).
- Additional blink-break overhead when enabled: every shown blink break adds about `6 s` (`500 ms` pre-cushion + `5 s` countdown + `500 ms` post-cushion); breaks that coincide with already scheduled transition/end-of-block interruption points are skipped.

Within occluded motion trials (frame-index anchored):
- Deviance / first full occlusion at frame `130` (`~1083 ms`)
- Full occlusion required through frame `190`
- Nominal reappearance search start at frame `191` (`~1592 ms`)

## Run messaging and color cue

Run-family color cue was enabled (`runColorCueEnabled = true`) with two colors matched for luminance difference with background:
- green `[0 255 0]`
- yellow `[242 223 0]`

Counterbalancing:
- odd subjects: run 1 green, runs 2/3 yellow
- even subjects: run 1 yellow, runs 2/3 green

Message flow defaults:
- start message enabled, waits for `1/8` response (`duration = 0`, no timeout)
- run1->run2 transition message: `5 s` with countdown enabled
- end-of-block message enabled, waits for `1/8` response
- optional between-block calibration-choice gate enabled (`1=skip`, `8=calibrate`) when eye-tracker session is active
- no separate final message screen in runtime flow; final block goes through the same end-of-block message pathway
- post-message fixation: `3 s`
- response lockout on message gates: `1 s`

## Trigger/event codes (v10 rescueTraject map)

- Condition onset: `71` (`always_visible`), `81` (`occluded_nondeviant`), `91` (`occluded_deviant`)
- Sequence identity: `1..N` at frame `12` (`N` usually `40`)
- Occlusion events: `111` (start), `112` (complete), `114` (end start), `115` (end complete)
- Type-1 catch visibility events: `116` (disappear), `117` (reappear)
- Catch question status: `118` (question start), then exactly one of `201` (correct), `202` (incorrect), `113` (timeout)
- Gaze/replay/termination: `150` (gaze break), `151` (first replay start in run), `152` (ESC termination)

Runtime validates that dynamic sequence codes (`1..N`) do not overlap fixed trigger values.
