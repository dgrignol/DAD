# MoveDot1 Experiment Blueprint (Work Guide)

Status: Draft (living document)
Last updated: 2026-02-24

## Why This Document Exists

This project can easily drift into "cool analysis" without a clean experiment that
actually tests a falsifiable scientific hypothesis. This document is meant to:

- Keep the project centered on a *data-collection-ready* experiment.
- Translate high-level theory (predictive coding in dynamic vision) into
  concrete, testable hypotheses.
- Make design tradeoffs explicit (what we choose, why we choose it, and what
  it costs).
- Map each scientific question to: stimulus requirements, trigger timing,
  acquisition constraints, and analysis requirements.

If you change the experiment (stimulus generation, trial structure, triggers,
timing), update this document and link to the exact file/parameter that changed.

## How To Use This Doc

Suggested workflow (example):

1. Start in "Scientific Questions and Hypotheses" and pick the smallest set of
   hypotheses that are worth a full MEG study.
2. Read "Experiment Factors" and "Stimulus Specification" to verify the design
   manipulations and the constraints are compatible.
3. Use "Implementation Map" to locate the canonical code that implements each
   part (stimuli, trial struct, runtime, triggers).
4. Before collecting real data, run the "Quality Gates" checklist.

## Canonical Sources Of Truth In This Repo

- Stimulus generation: `experiment/stimuli_generation_v14.m` + `experiment/lib/Config.m`
- Trial structure (catch timing, block order): `experiment/CreateInputFiles_v11.m`
- Experiment runtime (presentation, responses, eye tracker, MEG triggers):
  `experiment/MoveDot1_experiment_vX.m`
- Trigger definitions: `experiment/trigger_codes.md`
- Trigger verification pipeline: `trigger_pipeline.md` (+ `expected_triggers.py`,
  `inspect_fif_report.py`, `plot_trigger_chan.py`, `trigger_events.py`)
- MEG lab procedure (operational protocol): `experiment/MEG_lab_procedure.md`

The descriptions below should match these files. When in doubt, trust the code.

## Scientific Questions And Hypotheses

### Core Question

When the visual system observes a smooth, dynamic stimulus (moving dots), does
neural activity primarily reflect:

- an *updated sensory representation* of what just happened (lagged encoding), or
- an *anticipatory representation* of what is about to happen (predictive encoding),
  and how do these representations evolve within a trial?

### Manipulation: Predictability Disruption (Deviant)

We introduce a deviant (an abrupt change in direction and optionally curvature)
mid-trial. This allows dissociating:

- representation of the ongoing predicted trajectory,
- representation of the new (deviant) trajectory,
- and a candidate "prediction error" (PE) representation driven by mismatch.

### Hypotheses (draft, refine before preregistration)

H1: Within predictable trials, the timing of the strongest stimulus-brain match
shifts toward smaller positive lags and/or negative lags across trial time
(stronger prior -> more anticipatory representation).

Competing H1a (sharpening): lagged sensory representations strengthen over trial
time (stronger encoding of what is observed as the system becomes tuned).

Competing H1b (no prediction): representations remain primarily lagged and stable
over trial time.

H2: After the deviant onset, there is a transient PE-like representation before
the representation switches to the new (deviant) trajectory.

H3 (attention check / optional scientific target): Attended-dot representations
are stronger and/or earlier than unattended-dot representations (block-level
attention manipulation).

Each hypothesis must be linked to a concrete analysis metric (see "Analysis
Plan") and a pre-specified contrast.

## Experiment Factors (What We Manipulate)

These are the manipulations currently implemented in code:

1. Predictability condition (likelihood stimulus type)
   - Non-deviant / more predictable: `directionVariance = 0`
   - Deviant / less predictable: `directionVariance = 45`
   - Implemented in `experiment/lib/Config.m` and `experiment/stimuli_generation_v14.m`

2. Attention block condition (which dot is task-relevant)
   - Block condition 1: attend DotColor1
   - Block condition 2: attend DotColor2
   - Implemented in `experiment/CreateInputFiles_v11.m` (BlockOrder) and
     `experiment/MoveDot1_experiment_vX.m` (block messages + trigger mapping)

3. Catch trials (behavioral engagement / attention enforcement)
   - Two catch types exist in code: fixation catch and occlusion catch.
   - As of `experiment/CreateInputFiles_v11.m`, fixation catches are disabled
     (`Catch.CatchRatioFix = 0`) and occlusion catches dominate.
   - Catch timing and overlap constraints are in `experiment/CreateInputFiles_v11.m`.

## Trial Timeline (Conceptual)

The precise presentation logic lives in `experiment/MoveDot1_experiment_vX.m`,
but conceptually:

1. Trial starts: first frame of dot motion (MEG onset trigger is sent here).
2. Smooth dot motion for ~2.67 s (120 Hz by default).
3. Deviant onset (deviant condition only) occurs at:
   - `Config.deviantOnset * Config.trialDuration`
   - With current defaults: `0.5 * 2.67 s ~= 1.335 s`
4. Trial ends: last frame of dot motion (MEG end trigger is sent here).
5. Catch trials (if present) introduce a short within-trial manipulation window
   and use separate catch triggers.

## Stimulus Specification

### What A "Stimulus" Is Here

For each trial, the stimulus is a time series of dot positions:

- `xySeqs(...).xy` is `framesPerTrial x 4` with columns `[x1 y1 x2 y2]`
- Units are visual degrees inside the dot rectangle `[0..rectSize]`
- The runtime converts degrees to pixels and draws both dots each frame

This matches `experiment/stimuli_generation_v14.m` (generation) and
`experiment/MoveDot1_experiment_vX.m` (presentation).

### Current Generation Defaults (as-of v14 + Config.m)

Key parameters (verify in `experiment/lib/Config.m`):

- Trial duration: `Config.trialDuration = 2.67` s
- Frame frequency: `Config.frameFrequency = 120` Hz
- Dot rectangle size: `Config.dotRectSize = [10, 10]` deg
- Dot speed: `Config.dotSpeedDegPerSec = 3.73` deg/s
- Deviant onset: `Config.deviantOnset = 0.5` (fraction of trial)
- Predictability levels: `Config.likelihood.directionVariance = [0, 45]`
- Two-dot separation: `Config.minDistanceBetweenDots = Config.dotWidth * 2`

Important v14-specific property:

- v14 uses boundary-safe placement (feasible start-position ranges) and an
  optional dRSA-proxy gate to reduce position-vs-direction cross-correlation
  (see `experiment/stimuli_generation_v14.m` header and local settings).

### Design Choices And Rationale (Editable)

The theoretical goals impose constraints on the stimuli.

1. Smooth motion with simple latent parameters
   - Constant speed per trial: supports stable inference and a strong prior.
   - Mostly constant curvature within trial: supports prediction of trajectory
     evolution after early frames.

2. Deviant mid-trial (predictability break)
   - Provides a controlled "unexpected event" to probe prediction error and
     switching dynamics.
   - Implemented by an instantaneous direction change at deviant onset; v14 can
     also flip or randomize curvature post-onset for the deviant path only.

3. Two dots with block-level attention
   - Supports an attention manipulation without changing stimulus statistics
     across blocks (only task relevance changes).
   - Requires enforcing a minimum inter-dot distance to avoid near-overlaps.

4. Analysis-driven constraint: model separability
   - For dRSA (and for regression/PCR that partials out correlated models),
     position and direction model RDMs must not be too collinear.
   - v14 explicitly shapes the accepted trial bank to reduce cross-model
     coupling, rather than hoping random sampling produces enough separability.

If a stimulus change improves one goal but harms another, record it in the
"Decision Log" with the expected impact on hypotheses.

## Behavioral Task (Participant View)

Current intent:

- Participants fixate centrally.
- Two dots move simultaneously within a rectangle.
- Block instruction indicates which dot is relevant (DotColor1 vs DotColor2).
- Catch events occur occasionally to enforce attention and/or fixation.
  - Current implementation uses occlusion-based catches by default.
- Responses are recorded via button box (MEG) or keyboard (non-MEG testing).

If catch trials are meant to test attention rather than just enforce it, that
must be specified as a hypothesis and analysis plan.

## Triggers And Logging (MEG + Eye Tracking)

Trigger definitions are canonical in `experiment/trigger_codes.md`.

High-level summary:

- Trial onset:
  - 1-80 encode attention block condition + predictability condition + sequence
  - 102 overrides onset when the trial contains at least one catch
- Catch window:
  - 100 catch start (first frame of glitching)
  - 101 catch end (first frame after glitch finishes)
- Trial end:
  - 81 end of no-catch trial
  - 103 end of catch trial
- Eye-tracking abort/replay:
  - 150 gaze break (trial aborted)
  - 151 replay start (once per block, before first replay trial)
- Responses:
  - 201 typical response pulse (response value 1 recoded as 201)

Implementation notes:

- Runtime emission is in `experiment/MoveDot1_experiment_vX.m`.
- The expected-vs-actual trigger validation workflow is documented in
  `trigger_pipeline.md`.

## Implementation Map (What To Edit For Which Change)

If you want to change:

- Trajectory statistics (curvature, speed, deviant, dot spacing):
  - Edit `experiment/lib/Config.m`
  - Regenerate stimuli with `experiment/stimuli_generation_v14.m`
- Catch frequency, durations, overlap rules, block randomization:
  - Edit `experiment/CreateInputFiles_v11.m`
  - Regenerate trial struct `experiment/input_files/SubXX_TrialStruct.mat`
- Presentation logic, timing, triggers, eye-tracker behavior:
  - Edit `experiment/MoveDot1_experiment_vX.m`
- Trigger meaning / numbering:
  - Edit `experiment/trigger_codes.md` and keep MATLAB + Python consistent
- Trigger verification:
  - Follow and extend `trigger_pipeline.md`

## Common Commands (Regenerate Inputs / Run QC)

MATLAB is installed at `/Applications/MATLAB_R2020a.app/bin/matlab`. Prefer
invoking it via the absolute path.

Generate stimuli (creates `MovDot_SubXX.mat` and `MovDot_SubXX_predicted.mat`):

```
/Applications/MATLAB_R2020a.app/bin/matlab -batch "run('experiment/stimuli_generation_v14.m')"
```

Generate trial structure (creates `SubXX_TrialStruct.mat`):

```
/Applications/MATLAB_R2020a.app/bin/matlab -batch "run('experiment/CreateInputFiles_v11.m')"
```

If MATLAB R2020a has issues starting on Apple Silicon, force Intel mode:

```
arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "run('experiment/stimuli_generation_v14.m')"
```

Trigger QC (expected vs actual):

```
python expected_triggers.py --subject 4 --run 1
python inspect_fif_report.py --subject 4 --run 1 --trigger-channel STI101
python plot_trigger_chan.py --subject 4 --run 1
```

## Data Products (What Files Exist)

Inputs (stimulus/trial structure):

- `experiment/input_files/MovDot_SubXX.mat`
  - Generated by `experiment/stimuli_generation_v14.m`
- `experiment/input_files/MovDot_SubXX_predicted.mat`
  - Analysis-only baseline paths for deviant conditions (v10+)
- `experiment/input_files/SubXX_TrialStruct.mat`
  - Generated by `experiment/CreateInputFiles_v11.m`

Runtime outputs:

- `experiment/output_files/MoveDot1_SUBXX_RUNYY.mat`
  - Contains `CatchOutput`, `Catch`, `output`, `Conf` (see runtime script header)
- `experiment/output_files/debug_actual_triggers_subXX_runYY.csv`
  - Only when `Conf.debug = 1` (runtime script)

Trigger QC outputs (Python pipeline):

- `derivatives/triggers/subXX/expected_triggers_subXX_runYY_blockN.csv`
- `derivatives/triggers/subXX/actual_triggers_subXX.csv`
- `reports/..._report.html`

## Analysis Plan (dRSA-First, Draft)

This section is intentionally high-level. Before real data collection, lock the
minimal analysis commitments required to test H1/H2 (and optionally H3).

### Candidate Stimulus Models (time-resolved)

Models should be defined per dot, per frame, and per condition:

- Position model: dot (x,y) over time
- Direction model: unit tangent / angle over time
- Predicted-vs-observed models around deviant onset:
  - "Predicted path" (no-deviant baseline)
  - "Deviant path" (observed)
  - PE model (difference between predicted and observed)

For interpretability, plan for controlling shared variance between models:

- Regress out / partial out correlated model RDMs (e.g., PCR approach)
- Explicitly handle autocorrelation-limited timing precision (model and brain)

### Primary Metrics (examples)

- dRSA peak lag (sign and magnitude): does neural representation lead or lag?
- dRSA peak amplitude: strength of representation
- Evolution of peak lag/amplitude over trial time (pre-deviant)
- Post-deviant time course: PE transient, then switch to deviant encoding

### Minimal Prereg Items To Lock (before real data)

- Exact condition contrasts (predictable vs deviant; attended vs unattended)
- Primary time windows of interest (pre-registered or justified)
- The operational definition of "predictive" vs "lagged" representation
- The regression/partialing strategy across correlated stimulus models
- Handling of gaze breaks/replays and catch trials (exclude vs model separately)

## Quality Gates Before Collecting Real MEG Data

Stimulus sanity:

- Confirm paths stay within bounds by construction (v13/v14 feasible placement).
- Confirm deviant timing and magnitude are as intended.
- Confirm min distance between dots holds across frames.
- Confirm predicted baseline file matches deviant pre-onset exactly.

Trigger sanity (end-to-end):

- Run a short test session (`Conf.testingMode = 1`) and log debug triggers.
- Verify trial onset/end triggers match `experiment/trigger_codes.md`.
- If using eye tracker, verify 150/151 behavior (real or fake modes).
- Run the Python trigger pipeline on a recorded sample and inspect the report.

Behavior sanity:

- Verify catch trials are perceivable and response collection works.
- Decide whether missed catches are acceptable and how they affect analysis.

## Open Decisions / Questions (Keep Short, Resolve Early)

These should be resolved before committing to a full data-collection campaign:

- How many unique trajectories per condition are needed for generalization vs
  stable dRSA estimates? (`Config.trialsPerCondition` and trial struct runs)
- Deviant definition: direction change only vs direction + curvature modulation
  (v14 supports both; pick one and justify).
- Catch trials: purely engagement vs a scientific manipulation (then prereg).
- Whether fixation catches should be re-enabled (currently disabled in v11).
- Whether speed should be fixed (current) or varied (adds complexity and model
  coupling; requires analysis justification).

## Milestones (Suggested)

1. Lock hypotheses + minimal analysis commitments (H1/H2, optionally H3).
2. Lock stimulus parameter set (Config) and generate final input files.
3. Lock trial-struct generation parameters (catch frequency/durations).
4. Pilot with triggers + eye tracker; run QC pipeline; fix issues.
5. Finalize prereg-ready protocol and start real data collection.

## Decision Log (Fill As You Decide)

Template:

- YYYY-MM-DD: Decision summary
  - What changed:
  - Why:
  - Expected impact on hypotheses:
  - Code refs (file + parameter):
