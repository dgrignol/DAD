# dRSA Autocorrelation-Removal Strategies (2026-03-13)

## Scope

This document explains the four `params.PCR.RegressStrategy` backends implemented in:

- `simulations/functions/dRSA_PCR.m`

and how they are configured and propagated in:

- `simulations/functions/dRSA_coreFunction.m`
- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m`

The focus is on three new alternatives:

1. `ridge_full_autocorr`
2. `ridge_tapered_autocorr`
3. `ar1_prewhite_ridge`

while preserving legacy behavior via:

4. `baseline_pcr_border`

## Where to read first in code

- Strategy entry point and high-level flow:
  - `simulations/functions/dRSA_PCR.m:63-172`
- New strategy switch:
  - `simulations/functions/dRSA_PCR.m:102-165`
- Baseline (legacy) path:
  - `simulations/functions/dRSA_PCR.m:103-121`
  - helper: `simulations/functions/dRSA_PCR.m:284-312`
  - helper: `simulations/functions/dRSA_PCR.m:397-448`
- Ridge paths:
  - `simulations/functions/dRSA_PCR.m:122-162`
  - ridge solver helper: `simulations/functions/dRSA_PCR.m:515-582`
- AR(1) prewhitening helper path:
  - `simulations/functions/dRSA_PCR.m:207-260`

## Strategy catalog

### 1) `baseline_pcr_border` (legacy)

### Concept

This is the historical implementation: build a predictor matrix with:

- target model at current model-time sample,
- nuisance predictors from other models,
- nuisance own-model lag predictors selected via border (`dRSA_border`),

then run PCA and regress neural RDM on PCA scores, and project coefficients back.

### Practical implementation

- Selected in strategy switch:
  - `simulations/functions/dRSA_PCR.m:103-121`
- Border-style autocorr nuisance construction:
  - `simulations/functions/dRSA_PCR.m:284-312`
- Other-model nuisance block:
  - `simulations/functions/dRSA_PCR.m:268-282`
- Optional normalization preserving old behavior:
  - `simulations/functions/dRSA_PCR.m:377-395`
- Original PCR math:
  - `simulations/functions/dRSA_PCR.m:397-448`

### Notes

- This path still requires `Autocorrborder` when `RegressAutocor = 1`.
- In `PE_simulation_diff_1Dot.m`, border is computed only for this strategy:
  - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:667-676`

## 2) `ridge_full_autocorr`

### Concept

Replace PCA+OLS with ridge regression while keeping explicit nuisance regressors. The design matrix contains:

- first column: target predictor (`XTest`),
- other-model regressors over the local time window,
- own-model lag regressors over the same window, excluding lag-0.

Ridge stabilizes coefficients when nuisance regressors are collinear.

### Practical implementation

- Strategy branch:
  - `simulations/functions/dRSA_PCR.m:122-162`
- Full-window autocorr predictor block:
  - `simulations/functions/dRSA_PCR.m:138-140`
  - helper implementation: `simulations/functions/dRSA_PCR.m:314-327`
- Ridge solve (first predictor unpenalized):
  - `simulations/functions/dRSA_PCR.m:515-582`

### Key parameter

- `params.PCR.RidgeLambdaFactor` controls `lambda = factor * nObs`:
  - read at `simulations/functions/dRSA_PCR.m:184-196`

## 3) `ridge_tapered_autocorr`

### Concept

Same regression structure as `ridge_full_autocorr`, but with lag-dependent ridge penalties (not predictor scaling). Near lags get lighter penalty and far lags get stronger penalty. This preserves standardization while still favoring local autocorr control and reducing hard-edge artifacts.

### Practical implementation

- Strategy branch and penalty-profile call:
  - `simulations/functions/dRSA_PCR.m:142-155`
- Tapered penalty helper details:
  - `simulations/functions/dRSA_PCR.m:342-386`
- Penalty-vector assembly helper:
  - `simulations/functions/dRSA_PCR.m:388-416`

### Key parameter

- `params.PCR.TaperSigmaFactor` controls Gaussian width relative to local half-window:
  - read and applied at `simulations/functions/dRSA_PCR.m:362-365`
- `params.PCR.AutoPenaltyStrength` controls how much extra shrinkage far lags receive:
  - read and applied at `simulations/functions/dRSA_PCR.m:367-370`
- `params.PCR.NormalizeAutoPenaltyPerRow` keeps mean autocorr penalty stable across rows:
  - read and applied at `simulations/functions/dRSA_PCR.m:373-383`

## 4) `ar1_prewhite_ridge`

### Concept

First, reduce temporal autocorrelation in model and neural RDM time series via pooled AR(1) prewhitening. Then run ridge using:

- target predictor,
- other-model regressors,
- no own-model lag regressors.

This shifts autocorr handling from explicit own-lag nuisance columns to a preprocessing step.

### Practical implementation

- Strategy-level prewhitening entry:
  - `simulations/functions/dRSA_PCR.m:207-229`
- AR(1) whitening helper:
  - `simulations/functions/dRSA_PCR.m:231-260`
- Strategy branch excluding own-lag nuisance block:
  - `simulations/functions/dRSA_PCR.m:148-151`

### Key parameter

- `params.PCR.AR1Clip` clamps pooled AR estimate to keep whitening stable:
  - `simulations/functions/dRSA_PCR.m:216-220`
  - `simulations/functions/dRSA_PCR.m:254-255`

## Why Ridge (Not Lasso/Elastic Net) In This dRSA Pipeline

### Problem structure in this dataset

At each model-time sample, the nuisance design matrix contains many lag regressors that are highly correlated with each other and with nearby nuisance terms. In this setting, we need a method that:

- handles strong multicollinearity,
- keeps coefficient estimates stable across adjacent lags/timepoints,
- does not create arbitrary sparse lag selection that can look like step-like artifacts,
- keeps runtime manageable because this regression is solved repeatedly for many model/time combinations.

### Ridge vs Lasso in this context

Ridge and Lasso both regularize, but their bias is different:

- Ridge (`L2`) shrinks correlated nuisance predictors together.
  - This is well-aligned with autocorr nuisance structure, where nearby lags are expected to move as groups.
  - It produces smoother, more stable nuisance control and avoids hard on/off feature selection.
- Lasso (`L1`) promotes sparsity and tends to select one/few predictors among correlated groups.
  - With dense, correlated lag regressors, this selection can be unstable across adjacent rows/timepoints.
  - That instability can appear as piecewise transitions or band-like structure, which is the opposite of the target visual behavior.

### Why not Elastic Net (for now)

Elastic Net (`L1 + L2`) is a valid compromise and can reduce the worst Lasso instability, but in this pipeline it adds practical costs:

- an extra hyperparameter (`alpha`) plus lambda tuning,
- heavier optimization (iterative solvers) inside an already large repeated-regression loop,
- more tuning degrees of freedom that increase implementation risk and reproducibility complexity.

Given the goal (artifact reduction with predictable behavior and low implementation risk), ridge-based strategies were the most suitable first implementation choice.

### Practical decision rationale for this project

The implemented strategies were chosen to maximize stability and interpretability under strong nuisance collinearity:

- `ridge_full_autocorr`: stable baseline replacement for PCR with explicit autocorr nuisance terms,
- `ridge_tapered_autocorr`: same ridge stability plus lag-distance-aware shrinkage,
- `ar1_prewhite_ridge`: prewhitening plus ridge to reduce autocorr before regression.

This keeps the target-effect coefficient interpretation consistent while minimizing sparse-selection side effects.

## Shared design matrix flow (all non-baseline strategies)

For each tested model and model-time sample (`iT`):

1. Build `XTest` (target column):
   - `simulations/functions/dRSA_PCR.m:99-100`
2. Build window indices around `iT`:
   - `simulations/functions/dRSA_PCR.m:125`
   - helper: `simulations/functions/dRSA_PCR.m:262-266`
3. Build other-model nuisance block:
   - `simulations/functions/dRSA_PCR.m:127-133`
   - helper: `simulations/functions/dRSA_PCR.m:268-282`
4. Build strategy-specific autocorr block and corresponding penalty profile:
   - `simulations/functions/dRSA_PCR.m:135-167`
5. Concatenate and solve weighted ridge, then keep first beta row:
   - `simulations/functions/dRSA_PCR.m:169-174`
   - helper: `simulations/functions/dRSA_PCR.m:515-582`

## New parameters and defaults

### In `dRSA_coreFunction` defaults

Defined at:

- `simulations/functions/dRSA_coreFunction.m:90-94`
  - plus weighted-penalty defaults at `simulations/functions/dRSA_coreFunction.m:97-98`

Validation of allowed strategy names:

- `simulations/functions/dRSA_coreFunction.m:119-127`

### In one-dot script defaults

Set early so naming and run metadata are deterministic:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:135-170`

## Output naming and overwrite protection

To keep baseline and alternative outputs separate, non-baseline strategies append a strategy tag in `resultsBase`:

- tag injection:
  - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:350-360`
- strategy-to-tag mapping helper:
  - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:1767-1802`

Examples:

- `..._pcrridgeFullAuto_...`
- `..._pcrridgeTaperedAuto_...`
- `..._pcrar1PrewhiteRidge_...`

## Reproducibility metadata additions

Stored in each results MAT (`repro` struct):

- `pcrRegressStrategy`
- `pcrRidgeLambdaFactor`
- `pcrTaperSigmaFactor`
- `pcrStandardizePredictors`
- `pcrAR1Clip`
- `pcrAutoPenaltyStrength`
- `pcrNormalizeAutoPenaltyPerRow`

Code:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:876-882`

## Existing-results matching updates

The existing-results summary now compares PCR strategy too:

- call site:
  - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:387-393`
- comparison logic:
  - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:1387-1480`

## Validation outputs kept (Sub73)

Generated via three nondeviant PCR runs (one per strategy):

- `ridge_full_autocorr`
- `ridge_tapered_autocorr`
- `ar1_prewhite_ridge`

Artifacts live under:

- `simulations/output/Sub73_oneDot/nondeviant/pcr/results/`
- `simulations/output/Sub73_oneDot/nondeviant/pcr/matrices/commCbar/`
- `simulations/output/Sub73_oneDot/nondeviant/pcr/matrices/sepCbar/`
- `simulations/output/Sub73_oneDot/nondeviant/pcr/diagonal/`

with strategy-tagged filenames, for example:

- `dRSA_sub73_nondeviant_pcr_oneDot_obsOnly_pcrridgeFullAuto_results.mat`
- `dRSA_sub73_nondeviant_pcr_oneDot_obsOnly_pcrridgeTaperedAuto_results.mat`
- `dRSA_sub73_nondeviant_pcr_oneDot_obsOnly_pcrar1PrewhiteRidge_results.mat`
