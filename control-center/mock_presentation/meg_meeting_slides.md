---
title: "MoveDot1 MEG Experiment"
subtitle: "Dynamic representations, deviance, and predictive coding"
author: "DAD project - MEG meeting"
date: "2026-02-24"
---

## Literature state (brief)

- Predictive coding frameworks propose that perception combines top-down predictions and bottom-up prediction errors.
- Most work has focused on static or short discrete events, not continuously evolving trajectories.
- Dynamic representational analyses (RSA / dRSA) make it possible to test *when* a representation is lagged vs predictive.
- This project targets dynamic vision with MEG by modeling continuously changing dot trajectories and testing representational timing around deviance.

## Scientific gap and why now

- Key unresolved question: does neural coding of smooth motion remain mostly lagged, or does it become anticipatory within-trial?
- A controlled deviant (direction change mid-trial) allows us to dissociate:
  - ongoing predicted trajectory representation,
  - new deviant trajectory representation,
  - transient prediction-error-like representation.
- Terminology note (why \"expected\" nondeviant is still valid at 50/50 deviant events):
  - event-level probability is 50% deviant vs nondeviant,
  - but trajectory-level continuation probability is sharply peaked for the smooth nondeviant continuation relative to any *single* alternative direction.
- The current codebase now supports this test with:
  - explicit deviant and predicted baseline paths,
  - trigger-QC infrastructure,
  - dRSA-ready simulation and analysis scripts.

## Hypotheses and falsifiable contrasts

- H1 (predictive drift): during predictable motion, peak dRSA lag shifts toward smaller positive lags and possibly negative lags.
- H1a (sharpening alternative): lagged sensory representation gets stronger later in trial.
- H1b (null alternative): lag structure remains stable and mostly lagged.
- H2 (post-deviant dynamics): a transient PE-like representation appears first, then representation switches to deviant trajectory.
- Optional H3: attended-dot representations are stronger and/or earlier than unattended-dot representations.

Primary falsifiable contrasts:

- Pre-deviant lag slope: H1 vs H1a/H1b.
- Post-deviant temporal ordering: PE transient before deviant-path dominance.

## Experimental design overview

- Stimulus: two moving dots in a bounded visual rectangle; smooth motion within each trial.
- Predictability manipulation:
  - nondeviant (directionVariance = 0),
  - deviant (directionVariance = 45).
- Attention manipulation: block-wise attended dot (DotColor1 vs DotColor2).
- Catch logic (current v11 trial struct settings):
  - occlusion catches enabled,
  - fixation catches currently disabled (`Catch.CatchRatioFix = 0`).
- Trial organization:
  - 2 condition rows x 100 trajectories each in Sub80 input,
  - 10 runs x 2 attention blocks in trial-struct generation.

## Stimulus specification (Sub80)

Values from `experiment/stimuli_generation_versions.md` and current config/runtime assumptions:

- `trialDuration = 2.67 s`, `fps = 120`, `framesPerTrial = 320`.
- `pathDuration = 2.667 s` (likelihood stimulus type).
- `directionVariance = [0, 45]`.
- `dotSpeedDegPerFrame = 0.03108`.
- `curvFactor = 0.45`, `isCurvFactorRand = 1`.
- `minDistanceBetweenDots = 1.029 deg`.
- Deviant curvature flags: `flipCurvatureOnDeviant = 1`, `randomizeCurvatureOnDeviant = 0`, `deviantCurvatureRange = 0.45`.
- `dotWidth = 0.5144 deg`.

## Trial timeline and trigger map

Conceptual trial timeline:

- Motion onset (`t = 0`): trial onset trigger.
- Deviant onset (~mid-trial): expected around frame ~160 (~1.33 s for 2.67 s trial).
- Trial end: end trigger.

Canonical trigger summary (`experiment/trigger_codes.md`):

- Onset without catch: 1-80 (encodes block, condition, sequence).
- Onset with catch: 102.
- Catch start/end: 100 / 101.
- Trial end: 81 (no catch), 103 (with catch).
- Eye-tracker control: 150 (gaze break), 151 (replay start).
- Response: 201 typical (`response = 1`), 202 possible but flagged unexpected.

## Planned analysis pipeline (dRSA-first)

- Build time-resolved stimulus models per dot:
  - position,
  - direction,
  - predicted path,
  - deviant path,
  - PE (predicted vs observed mismatch).
- Compute dRSA matrices with **x = neural time** and **y = model time** for each model contrast.
- Then extract diagonal/lag summaries over trial time.
- Primary metrics:
  - peak lag sign and magnitude,
  - peak amplitude,
  - pre-deviant lag evolution,
  - post-deviant PE-to-deviant transition timing.
- Model-separability control:
  - partial/regression approach for correlated model RDMs,
  - explicit handling of autocorrelation-limited timing precision.
- Pre-data commitments to lock before full collection:
  - exact contrasts,
  - time windows,
  - exclusion/modeling of replay and catch trials,
  - regression strategy across correlated models.

## Expected result A: pre-deviant lag evolution (synthetic)

Interpretation target:

- H1 predicts a drift from lagged toward less-lagged/anticipatory representation.
- H1a predicts stronger lagged coding instead.
- H1b predicts stable lag profile.
- All panels are dRSA matrices (x = neural time, y = model time) in the pre-deviant epoch.

![](figures/fig03_mock_predrift_lag.png){ width=90% }

## Expected result B: post-deviant transition (synthetic)

Interpretation target:

- Early PE-like transient after deviant onset.
- Subsequent dominance of deviant-path representation over predicted-path representation.
- Matrix view should show an evolving lag profile around the deviant epoch.
- Left: neural(deviant) vs predicted model, middle: neural(deviant) vs deviant model, right: neural(deviant) vs PE model.

![](figures/fig04_mock_postdeviant_switch.png){ width=92% }

## Expected result C: realistic Sub80 trajectory assets

Real data-derived presentation assets:

- Sample trajectories for nondeviant vs deviant conditions.
- Timecourse of deviant-vs-predicted divergence (near-zero pre-onset, larger post-onset).
- Combined view connects geometry-level stimulus differences to model-level divergence.

![](figures/fig01_paths_sub80_realistic.png){ width=57% } ![](figures/fig02_deviant_predicted_divergence.png){ width=39% }

## Risks, open decisions, immediate milestones, and core references

Open decisions / risks:

- Required unique trajectories and repetitions for stable dRSA in real MEG data.
- Final preregistration choices for model partialing and timing windows.
- Handling replayed trials and missed catches in main analysis.
- Potential residual position-direction model coupling despite v14 controls.

Immediate milestones:

- Lock H1/H2 (and optional H3) test definitions and contrasts.
- Lock final stimulus + trial-struct parameter set.
- Pilot with trigger and eye-tracker QC, then freeze protocol.

Core references:

- Friston, K. (2005). A theory of cortical responses.
- Bastos, A. et al. (2012). Canonical microcircuits for predictive coding.
- Kriegeskorte, N. et al. (2008). Representational similarity analysis.
- Alink, A. et al. (2010). Predictive coding in visual cortex.
- Kok, P. et al. (2012). Less is more: expectation sharpens representations.
