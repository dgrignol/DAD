# MoveDot1 Experiment Blueprint (Work Guide)

Status: Draft (living document)
Last updated: 2026-02-26

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

### Hypotheses (draft)

In the predictable condition (non-deviant):

H1a prediction: after the first few hundreds milliseconds of only lagged representation, the brain starts to figure out the the trajectory parameters and form a stronger prior, consequently we should observe a predictive peak emerging. 
Concrete example (illustrative):
Early in the trial: best brain/stimulus match at lag ≈ +80 ms
Later in the same trial: match emerges at lag ≈ +20 ms or even lag < 0 indexing a predictive representation.
Interpretation: as the trajectory parameters become more certain within a trial, the brain forms an internal model (“stronger prior”), so anticipatory representations will emerge.

Competing H1b (no prediction): representations remain primarily lagged and stable
over trial time.

H2a (dampening): As a stronger prior is formed and anticipatory representations emerges, lagged representations starts to fade, as they are explained away by the prediction. Lagged peak reduces intensity over the course of the trial.

Competing H2b (sharpening): lagged sensory representations strengthen over trial
time (stronger encoding of what is observed as the system becomes tuned to the upcoming stimulation).


In the non-predictable condition (deviant):

H2 (vanilla PC): After the deviant onset, there is a transient PE-like lagged representation before
the predictive representation switches to the new (deviant) trajectory.

Competing H2a (Sharpened representation): After the deviant onset, brain continue representing predicted path for a small period and then swithces to the new (deviant) trajectory.

-TODO: complete the attention manipulation-related hyphotheses:
    Attention manipulation (trial-by-trial with two dots or one dot and central task): 

    H3 (Attention strengthen representation): Attended-dot predictive representations
    are stronger and/or earlier than unattended-dot.

    H3b (vanilla PC): worse lagged representations for predictable unattended, than predictable attended; ...



Each hypothesis must be linked to a concrete analysis metric (see "Analysis
Plan") and a pre-specified contrast.

## Experiment Factors (What We Manipulate)

These are the manipulations currently implemented in code:

1. Deviance condition (likelihood stimulus type)
   - Non-deviant: no deviant turn at the deviant onset (implemented as `directionVariance = 0`)
   - Deviant: deviant turn at the deviant onset (implemented as `directionVariance = 45` in current defaults)
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
- Deviance levels: `Config.likelihood.directionVariance = [0, 45]` (non-deviant, deviant)
- Two-dot separation: `Config.minDistanceBetweenDots = Config.dotWidth * 2`

Important v14-specific property:

- v14 uses boundary-safe placement (feasible start-position ranges) and an
  optional dRSA-proxy gate to reduce position-vs-direction cross-correlation
  (see `experiment/stimuli_generation_v14.m` header and local settings).

### Design Choices And Rationale (Editable)

The theoretical goals impose constraints on the stimuli.

1. We might need to add an intertrial period for baseline correction! https://www.sciencedirect.com/science/article/pii/S0165027021000157

2. Smooth motion with simple latent parameters
   - Constant speed per trial: supports stable inference and a strong prior.
   - Mostly constant curvature within trial: supports prediction of trajectory
     evolution after early frames.

3. Deviant mid-trial (predictability break)
   - Provides a controlled "unexpected event" to probe prediction error and
     switching dynamics.
   - Implemented by an instantaneous direction change at deviant onset; v14 can also flip or randomize curvature post-onset for the deviant path only.

4. Two dots with block-level attention
   - Supports an attention manipulation without changing stimulus statistics
     across blocks (only task relevance changes).
   - Requires enforcing a minimum inter-dot distance to avoid near-overlaps.

5. Analysis-driven constraint: model separability
   - For dRSA (and for regression/PCR that partials out correlated models),
     position and direction model RDMs must not be too collinear.
   - v14 explicitly shapes the accepted trial bank to reduce cross-model
     coupling, rather than hoping random sampling produces enough separability.

6. Deviance schedule across blocks/runs (**design choice under consideration**)
   - Motivation: a stable environment (only non-deviant paths) should let the
     trajectory prior keep strengthening across trials; a volatile environment
     (deviants can occur) should reduce confidence in the non-deviant path.
     This is expected to affect H1, including whether dRSA peak strength
     increases monotonically over within-trial time.
   - Option A (fully mixed, current default): deviant and non-deviant trials are
     interleaved within each attention block.
   - Option B (stable vs volatile contexts):
     - Stable blocks/runs: 100% non-deviant trials.
     - Volatile blocks/runs: deviant and non-deviant trials mixed (tentative:
       50/50).
   - Practical note: because blocks are already used for attention (attend
     DotColor1 vs DotColor2), adding stable/volatile as an additional *block*
     factor likely requires either (a) 4 blocks per run (attention x context),
     or (b) implementing stable vs volatile at the *run* level while keeping
     the existing two attention blocks per run.
   - Analysis implication: compare H1 within non-deviant trials between stable
     vs volatile contexts (e.g., slope/trend of dRSA peak lag/amplitude across
     within-trial time).
   - Tradeoff: a 50/50 deviant probability reduces "oddball surprise" and may
     reduce PE magnitude; if H2 is primary, consider lowering deviant
     probability in volatile blocks.

7. Inter-trial masking (**design choice under consideration**)
   - Consider adding a brief mask between trials to reduce sensory
     carryover from the previous trajectory (aftereffects).
   - Key tradeoff: stronger masking can reduce carryover but increases trial
     length and may alter participant state/arousal.
   - If adopted, lock mask type and duration in runtime code and document the
     expected impact on H1/H2 contrasts.

8. Deviant definition: direction-only vs direction+curvature
   (**design choice under consideration**)
   - Current framework supports either a pure direction change or a combined
     direction+curvature change at deviant onset.
   - Direction-only simplifies the formation of a new prior for participants, but adding curvature should not complicate that much.
   - Consistency with free parameter at the start of trial.
   - Adding curvature may increase ecological richness and decrease autocorrelation and cross-correlation.
   - Currently adopted the direction+curvature.

9. Catch role: enforcement-only vs scientific manipulation
   (**design choice under consideration**)
   - Option A: keep catch trials as compliance/engagement checks only.
   - Option B: treat catch effects as a hypothesis-driven factor (requires
     explicit preregistered contrasts and power planning).
   - This choice affects both participant instructions and the inferential
     scope of behavioral + MEG analyses.

10. Fixation catch policy (currently disabled)
   (**design choice under consideration**)
   - Fixation catches are currently disabled (`Catch.CatchRatioFix = 0`), while
     occlusion catches are used.
   - Re-enabling fixation catches may improve fixation monitoring but can alter
     attentional state and task demands.
   - If enabled, define target frequency/timing constraints and verify trigger
     behavior for mixed catch types.

11. Speed policy: fixed vs variable speed
    (**design choice under consideration**)
    - Fixed speed (current) supports cleaner latent-parameter inference and
      simpler model separability.
    - Variable speed may broaden generalization but can introduce extra
      covariance between direction, position, and time-varying kinematics.
    - If variable speed is introduced, predefine the speed distribution and the
      strategy for controlling added model coupling in analysis.

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
- Evolution of peak lag/amplitude over trial time
- Post-deviant time course only: PE

### Minimal Items To Lock (before real data)

- [ ] Exact condition contrasts (non-deviant vs deviant; attended vs unattended):
  Deviance and volatility context (H1/H2):
  - H1 (within non-deviant trials): contrast the evolution of dRSA peak lag/amplitude across within-trial time in the pre-deviant portion of the trial (early vs late pre-deviant windows, justified), or test a monotonic trend (e.g., peak amplitude increasing over within-trial time).
  - H1 (context modulation; if using stable vs volatile blocks/runs): use only non-deviant trials and contrast stable context (all non-deviant) vs volatile context (mixed deviant/non-deviant, tentative 50/50) on the H1 readout (e.g., slope/trend of dRSA peak lag/amplitude across within-trial time). Working prediction: stable context shows a stronger monotonic increase in peak strength and/or a stronger shift toward predictive lags.
  - Optional H1 between-trial-type check (in matched pre-onset windows): non-deviant vs deviant to test whether the *possibility* of a deviant (and/or the deviant trial type itself) reduces anticipatory shift even before the onset.
  - H2 (deviant response): time-lock to deviant onset and test a post-onset sequence in deviant trials:
    - predicted-path model (no-deviant baseline) vs deviant-path model (observed) vs PE model (predicted minus observed).
    - primary contrasts: early post-onset PE > (predicted, deviant); later post-onset deviant > predicted; and deviant vs non-deviant differences in the same post-onset windows (use a pseudo-onset in non-deviant trials for alignment).

  Attention (H3, optional):
  - attended vs unattended dot: compare dRSA peak lag/amplitude for the attended vs unattended dot as a function of block instruction (attend DotColor1 vs attend DotColor2), collapsing across deviance unless an interaction is explicitly targeted.
  - Optional interaction: (attended minus unattended) differs between deviant and non-deviant trials, especially post-onset.

- [x] Operational definition of "predictive" vs "lagged" representation (via dRSA peak lag, where `lag = t_brain - t_stim`):
  *predictive = peak lag < 0 ms, or peak lag < +100 ms (timing incompatible with a purely feedforward sensory delay);*
  *lagged = peak lag > +100 ms.*
- [ ] The regression/partialing strategy across correlated stimulus models
- [x] Handling of gaze breaks/replays and catch trials: 
  *exclude gaze-break trials; treat replayed trials as normal trials.*

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

- Verify catch trials are perceivable and response collection works. Run small Pilot with lab members.
- Decide which missed catches rate is acceptable.

## Open Decisions / Questions (Keep Short, Resolve Early)

These should be resolved before committing to a full data-collection campaign:

- How many unique trajectories per condition are needed for generalization vs
  stable dRSA estimates? (`Config.trialsPerCondition` and trial struct runs)
- Deviant definition: direction change only vs direction + curvature modulation
  (v14 supports both; pick one and justify).
- Deviance schedule / block context: keep a fully mixed design, or add stable
  contexts (all non-deviant) plus volatile contexts (mixed deviant/non-deviant;
  tentative 50/50). If adopted, decide whether context is implemented as blocks
  (attention x context implies 4 blocks per run) or as run types, and decide how
  many stable/volatile blocks or runs per subject (tentative minimum: 1 stable
  run + 1 volatile run; better: 2+2, counterbalanced order).
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
