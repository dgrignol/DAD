# Stimuli generation versions: path capabilities

## Input file parameter summary

This table summarizes key generation parameters per subject input file. Values
are read from `inferred_repro.config` (and `inferred_repro.stimulusTypeConfig`
for `pathDuration`/`directionVariance`) when available; if missing, the lookup
falls back to `repro` and then to the top-level `Cfg`. When none of these
structures expose a field (currently Subj 04 and 05), the cell is left as `—`.

| Input file | Subj | Source | curvFactor | isCurvFactorRand | trialsPerCondition | trialDuration | pathDuration | directionVariance | dotSpeedDegPerFrame | minDistanceBetweenDots | flipCurvatureOnDeviant | randomizeCurvatureOnDeviant | deviantCurvatureRange | dotWidth | isCurvValenceRand | directionChange | pathDurationVariance | isConnectedPaths | focusAroundCenterFactor | spreadoutTolerance | spreadoutToleranceInterval | deviantSignedTurnWindows | initialCurvatureWindows | deviantCurvatureWindows | deviantDisplacementRadiusRangeDeg | deviantDisplacementAngleWindowsDeg | deviantDisplacementMode | enableDrsaProxyGate | applyDrsaProxyGateToNondeviantOnly | deviantOnset | deviantOnsetVariance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `experiment/input_files/MovDot_Sub98.mat` | 98 | inferred_repro | 1.05 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub97.mat` | 97 | inferred_repro | 1.05 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub96.mat` | 96 | inferred_repro | 1.05 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub95.mat` | 95 | inferred_repro | 1.05 | 0 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub94.mat` | 94 | inferred_repro | 1.088 | 1 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub93.mat` | 93 | inferred_repro | 1.05 | 0 | 10000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub92.mat` | 92 | inferred_repro | 0.9 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub91.mat` | 91 | inferred_repro | 0.4 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 1 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub90.mat` | 90 | inferred_repro | 0.9 | 0 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 0 | — | — | — | — | 0.5 | [40, 120] | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub89.mat` | 89 | inferred_repro | 0.85 | 0 | 1000 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 0 | — | — | — | — | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub88.mat` | 88 | inferred_repro | 0.85 | 0 | 20 | 2.667 | 2.667 | [0, 45] | 0.03108 | — | — | — | — | 0.5144 | 0 | — | — | — | — | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub87.mat` | 87 | repro | 0.45 | 0 | 20 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | — | — | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub86.mat` | 86 | repro | 0.45 | 0 | 1000 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | — | — | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub85.mat` | 85 | repro | 0.45 | 1 | 100 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub84.mat` | 84 | repro | 0.45 | 1 | 100 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub83.mat` | 83 | repro | 0.45 | 1 | 100 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub82.mat` | 82 | repro | 0.45 | 1 | 100 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub81.mat` | 81 | repro | 0.45 | 1 | 100 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub80.mat` | 80 | repro | 0.45 | 1 | 100 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub79.mat` | 79 | repro | 0.45 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | 1 | 0 | (0.45) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub78.mat` | 78 | repro | 0.8 | 1 | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | (1) | 1 | 0.35 | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | — | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub77.mat` | 77 | repro | 0.8 | 1 | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | (1) | 1 | 0.35 | 0.5144 | 1 | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | — | — | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub76.mat` | 76 | repro | 0.8 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | (1) | 1 | (0.35) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | [-0.8, -0.3755, 0.3755, 0.8] | [-0.8, -0.3755, 0.3755, 0.8] | — | — | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub75.mat` | 75 | repro | 0.8 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | (1) | 1 | (0.35) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | [-0.8, -0.3755, 0.3755, 0.8] | [-0.8, -0.3755, 0.3755, 0.8] | [0.5, 1.8] | [-140, -40, 40, 140] | — | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub74.mat` | 74 | repro | 0.8 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | 1.029 | (1) | 1 | (0.35) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | [-0.8, -0.3755, 0.3755, 0.8] | [-0.8, -0.3755, 0.3755, 0.8] | [0.5, 1.8] | [-140, -40, 40, 140] | freeStart | — | — | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub73.mat` | 73 | repro | 0.8 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | (1.029) | (1) | 1 | (0.35) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | [-0.8, -0.3755, 0.3755, 0.8] | [-0.8, -0.3755, 0.3755, 0.8] | ([0.5, 1.8]) | ([-140, -40, 40, 140]) | off | 0 | (1) | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub72.mat` | 72 | repro | 0.8 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | (1.029) | (1) | 1 | (0.35) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | [-0.8, -0.3755, 0.3755, 0.8] | [-0.8, -0.3755, 0.3755, 0.8] | [0.5, 1.8] | [-140, -40, 40, 140] | freeStart | 1 | 0 | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub71.mat` | 71 | repro | 0.8 | (1) | 50 | 2.67 | 2.667 | [0, 45] | 0.03108 | (1.029) | (1) | 1 | (0.35) | 0.5144 | (1) | — | (0) | (0) | (0.5) | (0.5) | ([40, 120]) | [-80, -20, 20, 80] | [-0.8, -0.3755, 0.3755, 0.8] | [-0.8, -0.3755, 0.3755, 0.8] | [0.5, 1.8] | [-140, -40, 40, 140] | constrained | 1 | 0 | 0.5 | 0 |
| `experiment/input_files/MovDot_Sub05.mat` | 05 | Cfg | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — |
| `experiment/input_files/MovDot_Sub04.mat` | 04 | Cfg | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — |

This document summarizes how each `stimuli_generation_vXX.m` script (including
`stimuli_generation_V20.m`) differs in the
kinds of dot-motion paths it can generate. It focuses on *capabilities* and how
parameter options affect the output paths, not on any specific parameter values.

## Table parameter guide

This section documents each parameter currently included in the input-file table.
Descriptions emphasize path consequences and include examples when useful.

Legend: values in round brackets `(value)` are set in the saved config/repro but not used by that row's generator version or stimulus mode.

### Deviant control quick-toggle map

- **`[TOGGLE TURN]` `deviantSignedTurnWindows` + `directionVariance`**: turn deviance is active when a likelihood condition uses nonzero `directionVariance` and a non-empty allowed turn definition; set `directionVariance` to `[0]` to remove turn deviance.
- **`[TOGGLE CURVATURE]` `flipCurvatureOnDeviant` / `randomizeCurvatureOnDeviant`**: `flipCurvatureOnDeviant = true` enables sign-flip at onset, `randomizeCurvatureOnDeviant = true` enables post-onset resampling and takes precedence; set both to `false` for no curvature deviance.
- **`[TOGGLE DISPLACEMENT]` `deviantDisplacementMode` + displacement geometry fields**: displacement is active when mode is not `off` and geometry is non-empty; set mode to `off` to disable, or keep mode active and set radius to `[0, 0]` / angles to `[]` to suppress sampling.
- **`[DISPLACEMENT GEOMETRY]` `deviantDisplacementRadiusRangeDeg = [rInner rOuter]` and `deviantDisplacementAngleWindowsDeg = [aMin1 aMax1; ...]`**: defines annular-sector sampling for post-deviant translation.
- **`[DISPLACEMENT MODE]` `deviantDisplacementMode`**: `constrained` uses configured annular sectors, `freeStart` enables broad relocation sectors, and `off` disables displacement.
- **`[TOGGLE dRSA GATE]` `enableDrsaProxyGate` + `applyDrsaProxyGateToNondeviantOnly`**: enable gating with `enableDrsaProxyGate = true`; scope to nondeviant trials by setting `applyDrsaProxyGateToNondeviantOnly = true`.

- **`curvFactor`**: Scales curvature magnitude for versions that use factor-based curvature updates. Example: a higher factor yields more strongly bending trajectories frame to frame.
- **`isCurvFactorRand`**: Toggles whether curvature magnitude is resampled stochastically instead of fixed by the configured factor. Example: when enabled, trial-to-trial paths show wider bend-strength variability.
- **`trialsPerCondition`**: Controls how many trajectories are generated per condition combination. Example: increasing it gives denser sampling for analysis while any single path shape is unchanged.
- **`trialDuration`**: Sets the total path length in time, which determines how many motion updates occur per trial. Example: longer duration creates longer trajectories with more accumulated curvature.
- **`pathDuration`**: Defines the segment duration for segmented versions or condition metadata for single-path versions. Example: in segmented scripts, shorter values produce more within-trial direction-change segments.
- **`directionVariance`**: Defines likelihood-condition turn-deviance availability/magnitude when explicit signed windows are not overriding schedule construction. Example: larger nonzero values increase fallback deviant-turn magnitude.
- **`dotSpeedDegPerFrame`**: Sets spatial step size per frame in visual degrees. Example: higher speed stretches the path farther each frame and can increase boundary pressure.
- **`minDistanceBetweenDots`**: Minimum inter-dot separation for versions that enforce distance constraints. Example: increasing it rejects close approaches and yields more separated paired trajectories.
- **`flipCurvatureOnDeviant`**: Toggles post-deviant curvature sign inversion on the deviant path. Example: a clockwise pre-deviant bend can become counterclockwise after deviant onset, creating an S-like turn.
- **`randomizeCurvatureOnDeviant`**: Toggles post-deviant curvature resampling on the deviant path. Example: after onset, a gentle bend can switch to a stronger, weaker, or opposite-sign curvature.
- **`deviantCurvatureRange`**: Sets the symmetric sampling interval for randomized post-deviant curvature (used when randomization mode is active and explicit windows are not used). Example: wider ranges permit larger curvature jumps after the deviant frame.
- **`dotWidth`**: Defines dot size used for rendering and geometric constraints. Example: larger dots reduce usable interior area and can increase boundary or distance-based rejections.
- **`isCurvValenceRand`**: Toggles whether curvature sign is randomly sampled (+/-) or alternated deterministically in legacy/randomized-curvature flows. Example: enabling it increases clockwise/counterclockwise variability across paths.
- **`directionChange`**: Sets the width/scale of direction-change sampling in path_duration_norm mode. Example: larger values produce broader turn-angle variability at change points.
- **`pathDurationVariance`**: Jitters per-segment duration around pathDuration in segmented generators. Example: nonzero variance yields irregular segment lengths within a trial.
- **`isConnectedPaths`**: Toggles whether a new segment starts from the previous segment endpoint or from a newly seeded point. Example: disabling connection produces discontinuous multi-segment trajectories.
- **`focusAroundCenterFactor`**: Controls how strongly start-position sampling is concentrated near the stimulus center in center-biased versions. Example: smaller values keep starts more central.
- **`spreadoutTolerance`**: Sets the allowed fractional deviation from ideal grid occupancy during spread-out balancing. Example: smaller tolerance enforces more uniform spatial coverage across generated paths.
- **`spreadoutToleranceInterval`**: Encodes the absolute accepted occupancy interval derived from spread target and tolerance. Example: frames outside this interval trigger extra regeneration in spread-enforced versions.
- **`deviantSignedTurnWindows`**: Explicit signed turn windows for likelihood deviance. Example: excluding near-zero turns enforces stronger clockwise/counterclockwise deviant deflections.
- **`initialCurvatureWindows`**: Allowed signed intervals for baseline curvature sampling at trial start. Example: excluding near-zero intervals avoids near-straight starts.
- **`deviantCurvatureWindows`**: Allowed signed intervals for post-onset curvature sampling in randomization mode. Example: matching baseline windows keeps curvature magnitude regimes comparable before/after onset.
- **`deviantDisplacementRadiusRangeDeg`**: Inner/outer annulus radii for deviant displacement candidates. Example: increasing outer radius allows larger post-onset spatial jumps.
- **`deviantDisplacementAngleWindowsDeg`**: Allowed signed angle sectors for deviant displacement relative to pre-deviant heading. Example: [-140 -40; 40 140] limits jumps to side/back sectors.
- **`deviantDisplacementMode`**: Selects displacement behavior at deviant onset (`off`, `constrained`, or `freeStart`). Example: `off` disables displacement, `constrained` samples configured annular sectors, and `freeStart` permits broad relocation sectors.
- **`enableDrsaProxyGate`**: Enables the dRSA-proxy acceptance gate used during candidate-trial filtering. Example: when enabled, candidates that increase position-direction coupling can be rejected before acceptance.
- **`applyDrsaProxyGateToNondeviantOnly`**: Scopes dRSA-proxy gating to nondeviant-only trials when enabled. Example: setting this true keeps deviant-condition acceptance unconstrained by the proxy gate.
- **`deviantOnset`**: Base onset frame/time for deviant events in likelihood mode. Example: moving onset later shifts where turn/curvature/displacement manipulations begin along the trajectory.
- **`deviantOnsetVariance`**: Jitter around deviant onset for likelihood mode. Example: nonzero variance spreads deviant timing across trials instead of locking every trial to one onset frame.

## Base behavior (all versions)

- **Inputs and outputs**: Each script reads parameters from `Config` (plus user
  dialog inputs), seeds the RNG with the subject ID, and saves `xySeqs` (plus a
  `Cfg` struct) to `Config.inputDirectory` using `Config.stimuliFileName`.
- **Coordinate system**: `xySeqs(...).xy` stores per-frame dot coordinates in
  visual degrees within the dot rectangle. Column count depends on version:
  two-dot versions store `[x1 y1 x2 y2]`, while single-dot versions store
  `[x y]`.
- **Conditions**: Paths are generated for each combination of
  `directionVariance` × `pathDuration` × `trialsPerCondition` drawn from the
  stimulus-type config (`Config.likelihood`, `Config.path_duration`,
  `Config.path_duration_norm`).
- **Stimulus types (direction change logic)**:
  - `likelihood`: a deviant direction change occurs at a deviant onset (with
    optional onset variance). The deviant turn schedule can be version- and
    config-specific: legacy behavior derives turns from `directionVariance`,
    while v14 also supports explicit signed-turn windows. A *no-deviant
    baseline* path is computed and used to validate boundary constraints.
  - `path_duration`: direction change angles are sampled uniformly; the exact
    onset timing (segment start vs. early-trial) depends on the version.
  - `path_duration_norm`: direction change angles are sampled from a normal
    distribution with width `directionChange`; the onset timing depends on the
    version.
- **Motion update**: Direction evolves frame-by-frame by combining the current
  direction with a curvature term and any stimulus-type direction change, then
  stepping by `dotSpeedDegPerFrame`.
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
- **Reduce boundary-selection bias while keeping v12-style deviant curvature
  options**: use **v13**.
- **Constrain deviant turns with signed threshold-style windows**: use **v14**.
- **Set explicit start and post-deviant curvature windows**: use **v15**.
- **Add deviant onset displacement in an annular-sector region while keeping
  turn/curvature controls available**: use **v16_Displacement**.
- **Generate a one-dot trajectory set with v16-style turn/curvature/displacement
  controls and Config-toggled dRSA gate scope**: use **v17**.
- **Build the three-condition occlusion paradigm (always-visible, occluded
  nondeviant, occluded deviant) with twocircle default and configurable
  deviance-centered timing window `deviance-X` to `deviance+Y`**: use **v18**.
- **Build the occlusion paradigm with guaranteed matched pre-deviance branch
  for occluded deviant (`pre-deviance identical`, `post-deviance translated
  deviant suffix`)**: use **v19**.
- **Build the occlusion paradigm with fixed-frame timing (`deviance=130`,
  `full occlusion 130..190`) and geometric gradual reappearance from
  frame 191 onward**: use **V20**.
- **Generate occlusion stimuli directly from config (no source subject MAT)
  while keeping fixed-frame occlusion geometry (`deviance=130`,
  `full occlusion 130..190`, nominal reappearance search from frame 191)**:
  use **V21**.
- **Reduce residual position-direction coupling while preserving constant
  within-trial curvature**: use **v14**.
- **Experimentation version: half-length paths to reduce
  boundary-conditioned placement coupling**: use **v15_experimentalPathScale**.

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
- **Boundary-safe feasible placement and curvature-floor safeguards** (v13).

### v13 — `stimuli_generation_v13.m`

**Adds**
- **Boundary-safe feasible placement with relative paths**: trajectories are
  generated around the origin, then start positions are sampled from feasible
  ranges so both deviant and no-deviant baseline paths remain in bounds.
  *Example*: boundary validity is achieved by construction instead of selecting
  rare start points that survive full-trial boundary rejection.
- **Geometry-aware curvature floor option**: a local switch can clamp
  near-zero per-trial curvature magnitudes to a rectangle-derived minimum when
  needed. *Example*: very straight trajectories that cannot fit the full trial
  are suppressed before placement, reducing regeneration loops.
- **Optional baseline min-distance enforcement**: the analysis-only no-deviant
  baseline path can be checked against `minDistanceBetweenDots` in addition to
  the observed path. *Example*: predicted paths can be kept within the same
  geometric separation regime as presented paths.
- **Attempt guard for incompatible settings**: generation aborts after a
  configurable maximum attempts per trial, which adds explicit failure behavior
  when constraints are mutually incompatible (no direct path-shape effect).

**Omits (relative to other versions)**
- **Whole-trial boundary-hit rejection as the primary bounds mechanism** (v06,
  v08-v12).
- **dRSA-proxy-aware trial-selection gate** (v14).

### v14 — `stimuli_generation_v14.m`

**Adds**
- **Signed deviant-turn windows for likelihood mode**: optional
  `deviantSignedTurnWindows` lets you define allowed signed turn intervals and
  sample deviant turns only from those intervals; when omitted, v14 falls back
  to the legacy `directionVariance`-derived schedule. *Example*: exclude
  near-zero turns so deviants always exceed a chosen absolute-turn threshold
  while preserving both clockwise and counterclockwise options.
- **dRSA-proxy-aware trial gate**: candidate trials are evaluated with a
  reduced dRSA proxy that matches the target cross-model branch
  (position-RDM euclidean vs direction-RDM cosine), and non-improving
  candidates are regenerated. *Example*: trajectories remain smooth, but the
  accepted trial bank is explicitly shaped to reduce position-dot1 vs
  direction-dot1 coupling in corr-style matrices.
- **Relative-improvement acceptance rule with optional hard cap**: the gate can
  require a minimum score improvement over the current accepted-trial proxy,
  plus an optional absolute score ceiling. *Example*: generation can enforce a
  monotonic decline in the proxy score instead of only filtering out extreme
  outliers.
- **Condition-scoped and frame-strided proxy control**: gating can be limited
  to nondeviant likelihood conditions and evaluated on sampled frames to
  control compute load. *Example*: nondeviant trajectories can be decorrelated
  more aggressively while deviant manipulations remain unchanged.
- **Expanded reproducibility metadata for proxy-gate parameters and failures**:
  saved repro fields include proxy-gate settings and rejection counts.
  *Example*: downstream analyses can separate boundary/min-distance failures
  from dRSA-proxy gate pressure.

**Omits (relative to other versions)**
- **Within-trial piecewise curvature updates** (not used; curvature remains
  constant within each trial except optional deviant-point modulation).

### v15 — `stimuli_generation_v15.m`

**Adds**
- **Explicit baseline curvature windows**: initial per-dot curvature is sampled
  from `Config.initialCurvatureWindows` as signed interval windows.
  *Example*: enforcing `[-0.8, -0.3755] U [0.3755, 0.8]` removes near-zero
  curvatures at trial start while preserving clockwise and counterclockwise
  paths.
- **Explicit post-deviant curvature windows**: when deviant curvature
  randomization is enabled, post-onset curvature is sampled from
  `Config.deviantCurvatureWindows` instead of a single symmetric scalar range.
  *Example*: using the same signed windows before and after deviant onset keeps
  curvature magnitudes comparable across trial phases.
- **Signed deviant-turn windows retained from v14**: v15 keeps explicit
  likelihood deviant-turn control via `deviantSignedTurnWindows` in addition to
  curvature-window control. *Example*: you can jointly constrain turn-angle and
  curvature magnitudes at deviant onset.

**Omits (relative to other versions)**
- **Local path-scale displacement control** (v15_experimentalPathScale).

### v16_Displacement — `stimuli_generation_v16_Displacement.m`

**Adds**
- **Mode-controlled annular-sector deviant displacement**:
  deviant onset can include an instantaneous observed-path translation sampled
  from signed-angle windows relative to local pre-deviant heading, with mode
  controls that either disable displacement, keep radius-constrained sampling,
  or force full-circle displacement sampling at a board-scale radius.
  *Example*: switch from side/back constrained offsets to full-angle,
  screen-spanning offsets.
- **Composability with existing deviant turn and curvature logic**:
  displacement is applied after the trajectory is generated, so legacy turn
  scheduling and curvature modulation can stay active, be disabled, or be used
  jointly with displacement. *Example*: retain turn-angle deviations while also
  shifting the post-onset segment start.
- **Observed-only displacement with unchanged no-deviant baseline path**:
  baseline trajectory generation and export stay on the pre-existing
  no-displacement branch, preserving analysis comparability for predicted-path
  references.
- **Deviant-onset marker support for displacement events**: `pathAll` marks the
  onset frame when a displacement jump occurs, even when turn-angle change is
  disabled.
- **Reproducibility metadata for displacement settings**: saved repro
  parameters include displacement mode, configured/effective annulus bounds,
  signed-angle windows, and displacement enable state.

**Omits (relative to other versions)**
- **Local whole-path step scaling** (v15_experimentalPathScale).

### v17 — `stimuli_generation_v17.m`

**Adds**
- **Single-dot generation mode with v16-equivalent deviant controls**:
  keeps signed-turn scheduling, baseline/deviant curvature-window sampling,
  optional deviant curvature flip/resampling, and annular-sector deviant
  displacement, but generates only one observed dot trajectory per trial.
  *Example*: likelihood deviants can still combine turn and displacement
  manipulations while exporting a one-dot path matrix.
- **Single-dot output structure**: exports `xySeqs(...).xy` as `frames x 2`
  (`[x y]`) with matching one-column `pathAll`, `curvyness`, and
  `AngleDirection` metadata.
  *Example*: downstream analyses can operate on one trajectory stream without
  slicing out dot-2 columns.
- **Single-dot dRSA proxy gate**: preserves v14/v16-style trial acceptance
  gating by scoring position-vs-direction coupling on the generated dot, with
  gate activation and scope controlled from Config (`enableDrsaProxyGate`,
  `applyDrsaProxyGateToNondeviantOnly`).
  *Example*: keep gate `off` to avoid any condition-specific acceptance
  pressure, or enable it and choose nondeviant-only vs all-condition gating
  without editing the script.
- **Predicted baseline export retained in one-dot format**: deviant
  conditions still save no-deviant baseline trajectories for analysis, now as
  single-dot paths.
  *Example*: deviant-vs-baseline comparisons remain available after removing
  dot-2 generation.

**Omits (relative to other versions)**
- **Second-dot trajectory generation and metadata channels** (v05-v16 variants).
- **Inter-dot minimum-distance enforcement** (v05, v08-v16).
- **Two-dot plotting/preview rendering and `[x1 y1 x2 y2]` output layout**
  (v05-v16 variants).

### v18 — `stimuli_generation_v18.m`

**Adds**
- **Three-condition occlusion dataset builder from an existing one-dot
  source dataset**: creates `always_visible`, `occluded_nondeviant`, and
  `occluded_deviant` trial sets with matched trajectory pairing.
- **Configurable deviance-centered occlusion timing window**:
  user inputs `X` and `Y` to define `deviance-X` (occlusion start) and
  `deviance+Y` (reappearance start), defaulting to `0.25 s` and `0.25 s`.
- **Twocircle occlusion metadata for runtime rendering**:
  per trial exports center/radii/contact frames and post-deactivation frame
  with pre-circle active only before deviance and post-circle active from
  deviance onward.
- **Event-frame metadata for MEG triggers**:
  exports start/complete/end-start/end-complete frame fields so runtime
  scripts can emit occlusion event triggers deterministically.
- **Predicted-branch rule for occluded deviant simulations**:
  writes predicted trials where predicted path is forced to equal observed
  until reappearance start, then diverges from that frame onward.
- **Alpha fallback visibility profiles retained**:
  exports `visibility_alpha` and `visibility_geom` alongside twocircle fields
  to support optional alpha-based presentation paths.

**Omits (relative to other versions)**
- **Direct de-novo trajectory synthesis**:
  v18 transforms an existing one-dot source dataset rather than generating
  raw trajectories from scratch (v17 and earlier generation role).

### v19 — `stimuli_generation_v19.m`

**Adds**
- **Matched pre-deviance construction for occluded deviant**:
  occluded deviant trials are explicitly built so frames up to deviance are
  identical to the paired nondeviant trajectory.
- **Position-continuous deviant splice rule**:
  post-deviance branch is taken from deviant source trajectories and
  translated so the splice at deviance has no positional jump.
- **Non-interactive overwrite control**:
  supports optional `overwriteExisting=true` for batch regeneration without
  interactive overwrite prompts.
- **Updated reproducibility metadata**:
  outputs include a v19-specific construction rule label in `repro`.

**Keeps from v18**
- Three-condition occlusion outputs (`always_visible`, `occluded_nondeviant`,
  `occluded_deviant`) and twocircle metadata.
- Configurable `deviance-X` to `deviance+Y` timing and trigger-aligned event
  frame metadata.
- Predicted branch rule where predicted occluded deviant matches observed
  until reappearance start, then diverges.

### v20 — `stimuli_generation_V20.m`

**Adds**
- **Fixed-frame occlusion timeline with exact frame targets**:
  deviance is fixed at frame `130`, first full occlusion is fixed at frame
  `130`, and full occlusion is enforced through frame `190` (inclusive).
  The nominal reappearance search starts at frame `191`.
- **Fixed pre-occluder geometry at deviance**:
  the pre-deviance occluder is centered at the frame-130 dot position and
  uses a radius equal to the moving-dot radius.
- **Post-occluder hold and geometric reappearance**:
  the post-deviance occluder activates at frame `130`, has radius computed
  to fully occlude the dot for frames `130..190`, and remains active until
  trial end so reappearance is gradual and geometry-driven.
- **Single-deviance splice rule at frame 130**:
  the occluded-deviant pre-branch is copied from the paired nondeviant
  trial up to frame `130`; the source deviant post-branch is re-timed to
  start at frame `130` and translated for positional continuity, preventing
  a second inherited deviance point.
- **Twocircle-derived reappearance metadata**:
  `occlusion_end_frame` and `occlusion_end_complete_frame` are derived from
  geometric visibility state instead of a hard post-occluder deactivation
  frame.
- **Batch-safe fixed-frame controls**:
  supports non-interactive `fixedDevianceFrame` and
  `fixedOcclusionEndFrame` inputs (defaults `130` and `190`) with
  `overwriteExisting=true` for batch regeneration.

**Keeps from v19**
- Three-condition occlusion outputs (`always_visible`, `occluded_nondeviant`,
  `occluded_deviant`) and matched pre-deviance construction.
- Twocircle default + alpha fallback metadata and predicted output generation.

**Usage examples**
- Generate a V20 occlusion dataset non-interactively from repo root:
  `/Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); addpath('lib'); sourceSubjectID=73; targetSubjectID=99; fixedDevianceFrame=130; fixedOcclusionEndFrame=190; overwriteExisting=true; stimuli_generation_V20;"`
- If R2020a host-architecture detection fails on Apple Silicon, force Intel mode:
  `arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); addpath('lib'); sourceSubjectID=73; targetSubjectID=99; fixedDevianceFrame=130; fixedOcclusionEndFrame=190; overwriteExisting=true; stimuli_generation_V20;"`
- Create an occlusion preview video (twocircle mode) from the generated dataset:
  `python3 control-center/presentation/make_stimulus_movie_occlusion.py --mat-file /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/input_files/MovDot_Sub99.mat --condition occluded_deviant --mode twocircle --num-trials 1 --trial-selection first --show-event-markers --format both --output-stem /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/control-center/presentation/sub99_occluded_deviant_twocircle_v20_preview`

### v21 — `stimuli_generation_V21.m`

**Adds**
- **Config-driven one-dot trajectory synthesis with no source subject files**:
  V21 generates nondeviant/deviant/predicted trajectory triplets directly from
  `Config_stimuli_generation_V21` controls (curvature windows, deviant-turn windows,
  speed/frame geometry) instead of transforming pre-existing source MAT files.
- **Dedicated occlusion config class (`Config_stimuli_generation_V21`)**:
  keeps only one-dot parameters needed for this paradigm and introduces
  explicit fixed-frame occlusion timing constants.
- **Fixed-frame occlusion timeline preserved from V20**:
  uses `fixedDevianceFrame = 130`, `fixedOcclusionEndFrame = 190`, and a
  nominal geometric reappearance search start at frame `191`.
- **V20-compatible twocircle metadata and trigger-aligned event fields**:
  output trial structs remain compatible with
  `MoveDot1_experiment_occlusion_v1.m` and `trigger_codes.md`.
- **Predicted occlusion output in current tooling style**:
  writes `MovDot_SubXX_predicted.mat` with `xySeqsPredicted` entries labeled
  `occluded_deviant_predicted`, including the same occlusion metadata fields.

**Keeps from v20**
- Three observed occlusion conditions:
  `always_visible`, `occluded_nondeviant`, `occluded_deviant`.
- Pre-occluder at frame-130 position with dot-sized radius and post-occluder
  active from frame `130` through trial end.
- Geometric derivation of `occlusion_end_frame` and
  `occlusion_end_complete_frame`.

**Usage examples (MATLAB code style)**
- Interactive generation from `experiment/`:
  `addpath('lib');`
  `stimuli_generation_V21;`
- Override fixed-frame controls before generation:
  `addpath('lib');`
  `fixedDevianceFrame = 130;`
  `fixedOcclusionEndFrame = 190;`
  `overwriteExisting = true;`
  `stimuli_generation_V21;`

### v15_experimentalPathScale — `stimuli_generation_v15_experimentalPathScale.m`

**Adds**
- **Experimentation version: half-length paths to reduce boundary-conditioned placement coupling**:
  introduces a local step-scale control that shortens both observed and
  no-deviant baseline trajectories while preserving frame count, deviant
  onset logic, and constant-curvature-with-optional-deviant modulation.
  *Example*: trajectories occupy a tighter interior footprint, so feasible
  start placement has less edge pressure while path timing remains unchanged.
- **Single shared scaled-step path pipeline**: the same scaled displacement is
  used for observed-path stepping, no-deviant baseline generation, boundary
  feasibility ranges, and min-distance validation. *Example*: predicted
  baseline trajectories stay geometrically matched to the presented paths
  under the reduced-extent setting.
- **Repro metadata for step scaling**: saved outputs include the local scaling
  setting and effective step so downstream analyses can recover the exact
  spatial scaling regime used during generation.

**Usage notes**
- `pathScale` is a local scalar inside `stimuli_generation_v15_experimentalPathScale.m` (default
  `0.5`); adjust it there to tune spatial extent without editing `Config.m`.
- Keep `framesPerTrial` and deviant timing logic unchanged; only displacement
  per frame is scaled.
- If compared against v14 outputs, interpret differences as spatial-extent
  effects rather than curvature/deviant-rule changes.

**Omits (relative to other versions)**
- **Within-trial piecewise curvature updates** (not used; curvature remains
  constant within each trial except optional deviant-point modulation).
- **Explicit signed-window deviant-turn control in likelihood mode** (v14).
