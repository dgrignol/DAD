# Stimuli generation versions: path capabilities

## Input file parameter summary

This table summarizes key generation parameters per subject input file. Values
are read from `inferred_repro.config` (and `inferred_repro.stimulusTypeConfig`
for `pathDuration`/`directionVariance`) when available; if missing, the lookup
falls back to `repro` and then to the top-level `Cfg`. When none of these
structures expose a field (currently Subj 04 and 05), the cell is left as `—`.

| Input file | Subj | Source | curvFactor | isCurvFactorRand | trialsPerCondition | trialDuration | pathDuration | directionVariance | dotSpeedDegPerFrame | minDistanceBetweenDots | flipCurvatureOnDeviant | dotWidth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `experiment/input_files/MovDot_Sub98.mat` | 98 | inferred_repro | 1.05 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub97.mat` | 97 | inferred_repro | 1.05 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub96.mat` | 96 | inferred_repro | 1.05 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub95.mat` | 95 | inferred_repro | 1.05 | 0 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub94.mat` | 94 | inferred_repro | 1.088 | 1 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub93.mat` | 93 | inferred_repro | 1.05 | 0 | 10000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub92.mat` | 92 | inferred_repro | 0.9 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub91.mat` | 91 | inferred_repro | 0.4 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub90.mat` | 90 | inferred_repro | 0.9 | 0 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub89.mat` | 89 | inferred_repro | 0.85 | 0 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub88.mat` | 88 | inferred_repro | 0.85 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | 0.5144 |
| `experiment/input_files/MovDot_Sub87.mat` | 87 | repro | 0.45 | 0 | 20 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0.5144 |
| `experiment/input_files/MovDot_Sub05.mat` | 05 | Cfg | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub04.mat` | 04 | Cfg | — | — | — | — | — | — | — | — | — | — |

This document summarizes how each `stimuli_generation_vXX.m` script differs in the
kinds of dot-motion paths it can generate. It focuses on *capabilities* and how
parameter options affect the output paths, not on any specific parameter values.

## Base behavior (all versions)

- **Inputs and outputs**: Each script reads parameters from `Config` (plus user
  dialog inputs), seeds the RNG with the subject ID, and saves `xySeqs` (plus a
  `Cfg` struct) to `Config.inputDirectory` using `Config.stimuliFileName`.
- **Coordinate system**: `xySeqs(...).xy` is `frames x 4` with `[x1 y1 x2 y2]`
  in visual degrees within the dot rectangle. Positions are constrained to the
  rectangular stimulus area.
- **Conditions**: Paths are generated for each combination of
  `directionVariance` × `pathDuration` × `trialsPerCondition` drawn from the
  stimulus-type config (`Config.likelihood`, `Config.path_duration`,
  `Config.path_duration_norm`).
- **Stimulus types (direction change logic)**:
  - `likelihood`: a deviant direction change occurs at a deviant onset (with
    optional onset variance). The deviant magnitude comes from the
    `directionVariance` condition. A *no-deviant baseline* path is computed and
    used to validate boundary constraints.
  - `path_duration`: direction change angles are sampled uniformly; the exact
    onset timing (segment start vs. early-trial) depends on the version.
  - `path_duration_norm`: direction change angles are sampled from a normal
    distribution with width `directionChange`; the onset timing depends on the
    version.
- **Motion update**: Direction evolves frame-by-frame by combining the current
  direction with a per-dot curvature term and any stimulus-type direction change,
  then stepping by `dotSpeedDegPerFrame`.
- **Path metadata**: `pathAll` marks frames that start a new segment or contain
  the direction change, `curvyness` stores per-frame curvature, and
  `AngleDirection` stores per-frame direction angles.
- **Bounds handling**: All versions enforce the dot-rectangle boundaries, but
  *how* they enforce it (rejecting trials vs. placing paths) differs by version.

## Quick pick (choose by path requirement)

- **Segmented subpaths within a trial** (with optional connected segments and
  path-duration jitter): use **v05**.
- **Uniform start positions + simple whole-trial rejection** (no min-distance
  constraint): use **v06**.
- **No rejection bias + geometry-bounded curvature** (path placed to fit bounds):
  use **v07**.
- **Min-distance constraint enforced across frames** (vectorized generation):
  use **v08**.
- **Curvature flips after deviant onset** (S-shaped deviant paths): use **v09**.
- **Export no-deviant baseline paths for deviant conditions**: use **v10**.
- **Random post-deviant curvature in configurable +/- range**: use **v12**.

## Version differences from the base

### v05 — `stimuli_generation_v05.m`

**Adds**
- **Segmented subpaths per trial**: `pathDuration` controls the length of each
  subpath, and `trialDuration / pathDuration` determines how many subpaths are
  generated within a trial. *Example*: shorter `pathDuration` yields more
  direction-change segments per trial.
- **Path-duration jitter**: `Config.pathDurationVariance` jitters each subpath
  length, producing variable segment durations within the same trial.
- **Connected vs. reset segments**: `Config.isConnectedPaths` toggles whether
  each subpath starts where the previous one ended (`true`) or re-seeds near the
  center (`false`). *Example*: disabling connection produces discontinuities and
  more start markers in `pathAll`.
- **Center-biased starts**: `Config.focusAroundCenterFactor` concentrates starts
  around the rectangle center instead of uniform placement.
- **Spread-out enforcement**: `Config.spreadoutToleranceInterval` and
  `Config.spreadoutTolerance` enforce a grid-occupancy balance across trials,
  repeating generation until coverage is within tolerance. This biases paths
  toward uniform spatial coverage.
- **Per-path curvature updates**: curvature can be re-sampled each subpath using
  `Config.isCurvValenceRand` / `Config.isCurvFactorRand`, not just per trial.
- **Min-distance enforcement**: `Config.minDistanceBetweenDots` is checked at
  start and per frame; paths are regenerated if dots get too close.

**Omits (relative to other versions)**
- **Uniform full-rectangle starts** (v06+).
- **Whole-trial single-path generation** (v06+).
- **Geometry-bounded curvature without rejection** (v07).
- **Curvature flip at deviant onset** (v09+).
- **Exported no-deviant baseline file** (v10).

### v06 — `stimuli_generation_v06.m`

**Adds**
- **Uniform start positions**: both dots are initialized uniformly inside the
  rectangle. *Example*: removing center bias yields more edge-started paths.
- **Single continuous path per trial**: `pathDuration` no longer segments the
  path; it is stored as condition metadata only.
- **Whole-trial boundary rejection**: if any frame of the real or baseline path
  leaves bounds, the entire trial is regenerated.
- **Per-trial curvature**: curvature is sampled once per dot per trial via
  `Config.curvFactor` and the `isCurvValenceRand` / `isCurvFactorRand` flags.

**Omits (relative to other versions)**
- **Segmented subpaths, connection toggles, and path-duration jitter** (v05).
- **Spread-out grid enforcement** (v05).
- **Min-distance enforcement** (v05, v08+).
- **Geometry-bounded curvature / fit-in-bounds placement** (v07).

### v07 — `stimuli_generation_v07.m`

**Adds**
- **Fit-to-bounds placement**: paths are generated in relative coordinates, then
  a start position is chosen so *both* the deviant and no-deviant paths fit in
  bounds—no rejection sampling required. This removes boundary-related bias.
- **Geometry-bounded curvature**: curvature magnitude is constrained by the dot
  rectangle and trial length (computed `curvMinDeg`/`curvMaxDeg`). *Example*:
  shorter trials permit tighter curvature; larger rectangles allow gentler
  curvature while staying in bounds.
- **Per-trial, per-dot curvature with independent signs** (not scaled by
  `Config.curvFactor`).

**Omits (relative to other versions)**
- **Rejection-based boundary enforcement** (v06, v08–v10).
- **Min-distance enforcement** (v05, v08+).
- **`Config.curvFactor` scaling of curvature** (v05, v06, v08–v10).

### v08 — `stimuli_generation_v08.m`

**Adds**
- **Vectorized path generation**: per-frame angles and positions are computed
  with vectorized updates (functional change, same path family as v06).
- **Min-distance enforcement across frames**: `Config.minDistanceBetweenDots`
  rejects trials where dots approach too closely at any frame. *Example*: higher
  minimum distance yields more separated paths and higher rejection rates.
- **Existing-output guard**: can plot existing paths without regenerating (no
  path-shape change, but useful for inspection).

**Omits (relative to other versions)**
- **Geometry-bounded curvature / fit-in-bounds placement** (v07).
- **Curvature flip at deviant onset** (v09+).

### v09 — `stimuli_generation_v09.m`

**Adds**
- **Curvature flip on deviant**: `flipCurvatureOnDeviant` (local flag) flips the
  curvature sign from the deviant frame onward *for the deviant path only*.
  *Example*: a path that curves clockwise before the deviant can curve
  counterclockwise after, producing an S-shaped trajectory; the baseline path
  retains the original curvature for boundary checks.
- **Per-frame curvature record**: `curvyness` reflects the flip on a per-frame
  basis when enabled.

**Omits (relative to other versions)**
- **Exported no-deviant baseline file** (v10).

### v10 — `stimuli_generation_v10.m`

**Adds**
- **No-deviant baseline export for deviant conditions**: saves
  `MovDot_SubXX_predicted.mat` containing `xySeqsPredicted`, i.e., the
  no-deviant baseline paths that share the same start points and curvature but
  remove the deviant direction change. *Example*: enables direct deviant vs.
  baseline comparisons in analyses without regenerating paths.

**Omits (relative to other versions)**
- **No additional path-generation omissions beyond v09** (v10 retains v09
  options, including curvature flip).

### v12 — `stimuli_generation_v12.m`

**Adds**
- **Randomized post-deviant curvature option**: `randomizeCurvatureOnDeviant`
  can replace the post-onset curvature with a newly sampled value for each dot
  in the `[-deviantCurvatureRange, +deviantCurvatureRange]` interval, applied
  from deviant onset onward on the deviant path only. *Example*: a path that
  begins with gentle clockwise curvature can switch to a stronger or weaker
  curvature (including opposite sign) after the deviant frame.
- **Explicit precedence between deviant-curvature modes**: when both options are
  enabled, randomization takes precedence over `flipCurvatureOnDeviant`; the
  no-deviant baseline still keeps the original curvature profile.

**Omits (relative to other versions)**
- **No additional path-generation omissions beyond v10** (v12 retains v10
  options, including predicted no-deviant baseline export).

## Notes on parameter effects that change paths

These parameters exist across versions, but *only matter* when the version uses
that capability:

- **`directionVariance` (likelihood)**: larger values produce larger deviant
  turns at the deviant onset frame.
- **`deviantOnset` / `deviantOnsetVariance`**: shift the deviant turn earlier or
  later within the trial (likelihood only).
- **`directionChange` (path_duration_norm)**: wider distributions yield larger
  early-trial direction changes.
- **`curvFactor` + `isCurvValenceRand` / `isCurvFactorRand`**: increase curvature
  magnitude and/or randomize curvature sign (v05, v06, v08-v10, v12).
- **`randomizeCurvatureOnDeviant` + `deviantCurvatureRange`**: in v12, replace
  post-deviant curvature with a random value in a configurable interval.
- **`minDistanceBetweenDots`**: enforces separation when the version checks it
  (v05, v08-v10, and v12); larger values yield more separated trajectories.
