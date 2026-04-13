# Methods

## Stimulus generation (V27 block-resume)

Stimuli were generated with `stimuli_generation_V27_blockResume.m` and `Config_stimuli_generation_V27_blockResume.m`. Each trial contained a single moving dot in a `10 x 10 deg` arena. Dot diameter was `0.51442 deg`, speed was `3.73 deg/s`, frame rate was `120 Hz`, and trial duration was `2.67 s` (`320` frames). Viewing distance at runtime was `1000 mm` by default.

The dot center was bounded by arena margins to keep the full dot inside the `10 x 10 deg` window (`x,y` in approximately `[0.257, 9.743] deg`). With fixation at arena center (`5,5 deg`), this yields an approximate eccentricity span of:
- minimum about `0.557 deg` (from fixation-exclusion constraint in generation)
- maximum about `6.71 deg` (near arena corners)

For each subject, generation was deterministic (`rng(targetSubjectID)`) and produced `20` unique sequence IDs per presented condition. Presented condition labels and codes were:
- `always_visible` (`-1`)
- `occluded_nondeviant` (`0`)
- `occluded_deviant` (`45`, from `directionVariance = [0, 45]`)

For each sequence, nondeviant and deviant trajectories shared the same initial direction and were forced identical through the deviance frame. Baseline curvature was sampled from `[-0.8, -0.3755] U [0.3755, 0.8]`. Deviant signed turn at onset was sampled from `[-81, -10] U [10, 81]` deg. Post-onset deviant curvature was re-sampled from the same curvature windows (`randomizeCurvatureOnDeviant = true`, `flipCurvatureOnDeviant = false`). Deviant suffixes were translated at the splice frame to preserve positional continuity.

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

Scheduling was created by `CreateInputFiles_v20_threeRunsPerBlock_catch_blockResume.m` with `Config_schedule_CreateInputV20_MoveDotV17_blockResume.m`.

Default schedule:
- `10` blocks
- `3` runs per block
- Run 1: shuffled `always_visible`
- Run 2: shuffled first half of `occluded_nondeviant` + first half of `occluded_deviant`
- Run 3: shuffled remaining halves of both occluded conditions

With `20` base trials per condition, base run lengths are `20/20/20` before catch insertion. Runtime (`MoveDot1_experiment_occlusion_v17_blockResume.m`) supports block resume by starting at user-selected `iBlock` and continuing through the final scheduled block.

## Catch trials and task

Two catch types are inserted upstream in `CatchPlan`:

Type 1 (run 1 only):
- count per block: `round(0.10 * 20) = 2`
- trial begins as `always_visible`
- disappearance onset sampled from `0.30-1.00 s`
- invisible duration: `0.50 s` (`60` frames at `120 Hz`)
- on reappearance, trajectory either remains on source path or switches to a plausible deviant branch with probability `0.50` (when matched sequence exists)
- expected response: YES if branch changed, NO otherwise

Type 2 (runs 2+3 pooled):
- count per block: `round(0.10 * (20 + 20)) = 4` total across runs 2 and 3
- no added disappearance manipulation (standard occlusion trial)
- expected response: YES for `occluded_deviant`, NO for `occluded_nondeviant`
- catches are inserted into random slots; run-2/run-3 split varies by block

Default totals over 10 blocks:
- type 1: `20`
- type 2: `40`
- total catches: `60`

## Trial procedure and responses

Each trial displayed one trajectory. A central fixation cross was drawn during the motion phase, the ITI blank phase, and catch-question screens. No response was required on non-catch trials.

On catch trials, a post-trial question was shown for `4.0 s`:
- schedule default prompt: `"Has the dot changed its course?"`
- runtime preset in `Config_runtime_v17_blockResume.m` overrides this with: `"Did the trajectory change while occluded?"`

Response mapping:
- YES: code `1` (`RightArrow`, `8`, `y`; MEG button value `8`, blue)
- NO: code `2` (`LeftArrow`, `1`, `n`; MEG button value `1`, red)

Catch feedback (`500 ms`):
- correct: green fixation
- incorrect: red fixation
- timeout: `"Too slow"` text

Post-trial ITI is jittered with a uniform draw in a configurable range:
- runtime default: `[0.500, 1.000] s`
- runtime override: `ITIRangeSec = [min, max]`
- implementation detail: all ITIs are precomputed at block start (before trial execution) with deterministic seeding from participant identity plus block/run offsets.
- replay behavior: replayed trials reuse the ITI of their originating source trial via stored source-trial mapping.

The previous grid-mask ITI path (blank + low-luminance grid + blank) is kept in code as commented legacy settings for possible future reuse.

Current runtime display background settings are all black:
- `Conf.background = [0 0 0]`
- `Conf.rectColor = [0 0 0]`
- `Conf.rectBorderColor = [0 0 0]`
so, outside the dot/fixation/mask elements, the screen is black.

## Trial timing

Per-trial timing is:
- Motion epoch: `2670 ms` (`320` frames at `120 Hz`)
- Post-trial ITI: uniform jitter in `[500, 1000] ms` by default (fixation visible)

For catch trials:
- Catch question: `0-4000 ms` (ends on response or timeout)
- Catch feedback: `500 ms`
- Then standard jittered ITI from the same run precomputed ITI vector

So:
- Non-catch trial cycle (motion + ITI): `3170-3670 ms` (default ITI range)
- Catch trial cycle (motion + question + feedback + ITI): `3670-8170 ms` (default ITI range)
- Estimated per-block trial-loop duration (`66` trials: `60` non-catch + `6` catch; excluding breaks/messages/replays): approximately `212-269 s` (`3.5-4.5 min`).

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
- final message: `5 s`
- post-message fixation: `3 s`
- response lockout on message gates: `1 s`

## Trigger/event codes (v8 block-resume map)

- Condition onset: `31` (`always_visible`), `41` (`occluded_nondeviant`), `51` (`occluded_deviant`)
- Sequence identity: `1..N` at frame `12` (`N` usually `20`)
- Occlusion events: `111` (start), `112` (complete), `114` (end start), `115` (end complete)
- Type-1 catch visibility events: `116` (disappear), `117` (reappear)
- Catch question and response: `118` (question start), `101` (question end), `201` (YES), `202` (NO), `113` (timeout)
- Gaze/replay/termination: `150` (gaze break), `151` (first replay start in run), `152` (ESC termination)

Runtime validates that dynamic sequence codes (`1..N`) do not overlap fixed trigger values.
