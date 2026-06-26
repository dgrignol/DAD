# dRSA Sliding Lag Surfaces

## Purpose

`09_plot_drsa_sliding_lag_surfaces.py` is an exploratory post-processing stage
for Stage-07 dRSA matrices. It turns each selected model matrix into a
lag-by-model-time surface by averaging diagonal samples inside overlapping
model-time windows.

The convention is the same as Stage 08:

```text
lag = neural time - model time
sample = matrix[model_time_index, model_time_index + lag]
```

Stage-07 matrices are interpreted as `matrix[model_time, neural_time]`.

## Initial Scope

The first version targets occipital `A_dRSA_alt` position summaries:

- always_visible | actual visible position
- occluded_nondeviant | actual nondeviant position
- occluded_deviant | actual deviant position
- occluded_deviant | expected/predicted position
- occluded_deviant | prediction error position

Each target summary receives:

- a static lag-by-model-time heatmap
- an interactive Plotly 3-D waterfall

The stage also writes two peak-timecourse comparison figures:

- actual-position comparison: actual visible, actual nondeviant, actual deviant
- occluded-deviant model comparison: actual deviant, expected/predicted,
  prediction error

## Peak Timecourses

Peak timecourses are extracted from each lag-by-model-time surface. For every
model-time window:

```text
peak_magnitude = max over lags of mean dRSA(window, lag)
peak_lag = lag at peak_magnitude
```

The peak-timecourse figure has two panels. The top panel overlays peak
magnitude over model time. The bottom panel overlays peak lag over model time.
By default, peak lag is masked when the corresponding peak magnitude is below
20 percent of the maximum finite peak magnitude in the report. This avoids
showing an apparent latency when the lag curve is weak or flat.

## Interpretation Risks

These surfaces are descriptive and should not be treated as independent
time-resolved statistical tests without additional subject-level inference.

- Overlapping windows share model-time samples, so adjacent model-time rows are
  not independent.
- dRSA matrices, stimulus trajectories, and MEG responses can all be temporally
  autocorrelated; smooth lag ridges may reflect temporal structure rather than a
  sharply timed neural delay.
- Large positive and negative lags can have fewer valid diagonal samples near
  matrix edges.
- Prediction-error matrices can have undefined model-time support before the
  paths diverge, so early windows may be intentionally empty.
- The model-anchored rule restricts model-time rows only; neural-time columns
  can fall outside the displayed model-time window.
- Peak-lag overlays are argmax summaries of each window's averaged lag curve.
  They are useful for orientation but can jump because of noise, plateaus,
  smoothing, unequal valid sample counts, or multiple local maxima.
- Masked peak-lag timecourses are still descriptive. A threshold removes weak
  windows from the latency panel, but it does not make the remaining latency
  estimates inferentially reliable.

## Practical Reading

Use the heatmap as the primary comparison view because all five targets share
one color range within a report. Use the 3-D waterfall as an inspection tool for
rotating the same surface and checking whether apparent ridges are broad,
stable, or driven by a few windows. Use the peak-timecourse comparisons when
the main question is whether deviant and nondeviant conditions differ in peak
magnitude or peak lag after reappearance.
