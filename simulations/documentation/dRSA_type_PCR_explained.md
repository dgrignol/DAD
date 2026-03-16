# `dRSAtype = 'PCR'` Explained

## Scope

This note explains how `dRSAtype = 'PCR'` works in this repository, with emphasis on the legacy PCR path:

- `params.PCR.RegressStrategy = 'baseline_pcr_border'`

It is written for someone who wants both:

1. the conceptual reason PCR exists here, and
2. the exact implementation flow in the current code.

Primary code paths covered:

- `simulations/functions/dRSA_coreFunction.m`
- `simulations/functions/dRSA_PCR.m`
- `simulations/functions/dRSA_corr.m`
- `simulations/functions/dRSA_computeRDM.m`
- `simulations/functions/dRSA_border.m`
- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m`

## Executive Summary

In this codebase, both `corr` and `PCR` start from the same RDM representation:

- raw model/neural time series are converted into time-resolved RDM columns,
- each RDM column is a vector of pairwise distances across subsamples/trials,
- the final dRSA result is a matrix of `model time x neural time`.

The difference is what happens after the RDMs exist:

- `dRSAtype = 'corr'` directly correlates model-RDM columns with neural-RDM columns.
- `dRSAtype = 'PCR'` builds a regression problem at each model timepoint, keeps the tested model as predictor 1, adds nuisance predictors for other models and own-model autocorrelation, runs PCA-based regression, and stores only the coefficient for the tested predictor.

So the practical meaning is:

- `corr` asks: "How similar are these two RDM time courses?"
- `PCR` asks: "How much unique variance does this model explain after controlling for nuisance structure?"

## Shared Data Representation First

Before comparing `corr` and `PCR`, it helps to be explicit about what the code calls an RDM.

### Raw input shape

`dRSA_coreFunction` expects:

- `Y`: neural/data signal as `features x time`, or a precomputed RDM
- `model`: cell array of model signals, each also `features x time`, or precomputed RDMs

See:

- `simulations/functions/dRSA_coreFunction.m:5-23`

### What `CurrSubsamples` does

When raw time series are passed in, `CurrSubsamples` tells the code which samples belong to the same relative within-trial position. Each iteration slice is reshaped to `nSubsamples x subSampleDuration` in:

- `simulations/functions/dRSA_coreFunction.m:166-172`

Then `dRSA_computeRDM` loops over each within-trial timepoint and computes pairwise distances across subsamples/trials:

- `simulations/functions/dRSA_computeRDM.m:19-30`

Important detail:

- each column of the RDM is a condensed pairwise-distance vector,
- so an RDM has shape `nPairs x nTimepoints`,
- where `nPairs = nSubsamples * (nSubsamples - 1) / 2`.

This means the dRSA code does not compare raw signals directly. It compares how trial geometry evolves over time.

### Normalization before PCR

`dRSA_computeRDM` normalizes RDMs when either:

- `params.normalization` is explicitly set, or
- `params.dRSAtype == 'PCR'`

See:

- `simulations/functions/dRSA_computeRDM.m:86-107`

If PCR is requested and no normalization is provided, the code defaults to:

- `params.normalization = 'Standardize'`

at:

- `simulations/functions/dRSA_computeRDM.m:88-97`

This is important because the PCR code assumes the RDM columns are already centered/scaled before regression.

## `corr` vs `PCR`

| Aspect | `dRSAtype = 'corr'` | `dRSAtype = 'PCR'` with `baseline_pcr_border` |
| --- | --- | --- |
| Main computation | Direct column-wise correlation | PCA-based regression coefficient for the tested predictor |
| Core file | `simulations/functions/dRSA_corr.m` | `simulations/functions/dRSA_PCR.m` |
| One row of the final dRSA matrix | Correlation of one model-time column with all neural-time columns | Beta of one tested model-time predictor against all neural-time columns |
| Controls overlap with other models | No | Yes, if `params.PCR.RessModel = 1` |
| Controls own-model temporal autocorrelation | No | Yes, if `params.PCR.RegressAutocor = 1` and a border is available |
| Output units | Correlation coefficients | Regression weights, not correlations |
| Collinearity handling | None beyond raw correlation | PCA reduction before regression |

### What `corr` does exactly

`dRSA_corr` is very simple:

- it loops over `params.modelToTest`,
- and for each tested model it computes:
  - `corr(mRDM{iModel}, nRDM)`

See:

- `simulations/functions/dRSA_corr.m:12-22`

Because `corr` in MATLAB correlates columns of `mRDM{iModel}` with columns of `nRDM`, the result is a full `model time x neural time` matrix.

What it does **not** do:

- it does not regress out other models,
- it does not regress out autocorrelation,
- it does not try to isolate unique model variance.

### What `PCR` does instead

At each tested model timepoint `iT`, legacy PCR builds a predictor matrix like this:

1. first column = the model currently being tested at time `iT`
2. additional columns = nuisance predictors from other models in a local time window
3. additional columns = nuisance predictors from the tested model's own autocorrelation, chosen using the border logic

Then it:

1. runs PCA on that predictor matrix,
2. regresses the neural RDM on PCA scores,
3. projects the solution back to the original predictor space,
4. keeps only the coefficient of predictor 1, the tested model.

See:

- `simulations/functions/dRSA_PCR.m:92-181`
- `simulations/functions/dRSA_PCR.m:465-513`

## Why PCR Exists In This Pipeline

PCR exists here to solve two concrete problems that plain `corr` cannot solve well.

### Problem 1: models overlap with each other

In this repository, the tested model is usually not alone. The script builds:

- `params.modelToTest = 1:numel(model)`
- `params.modeltoRegressout{iModel} = setdiff(1:numel(model), iModel)`

in:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:663-685`

So when a model is tested, every other model becomes a nuisance regressor by default.

Why this matters:

- position and direction models are not guaranteed to be independent,
- predicted and observed models can also share structure,
- a plain correlation can look strong simply because the tested model covaries with another model that better matches the neural data.

PCR is the code's attempt to estimate the tested model's **unique** contribution.

### Problem 2: model RDMs are strongly autocorrelated across time

Nearby model timepoints tend to look similar. If that autocorrelation is not handled, a model can appear to match neural timepoints just because both sides have smooth temporal structure.

The legacy PCR strategy addresses this by:

- estimating a model-specific autocorrelation border with `dRSA_border`,
- then adding own-model nuisance predictors from selected lagged timepoints.

See:

- `simulations/functions/dRSA_border.m:63-119`
- `simulations/functions/dRSA_PCR.m:107-125`
- `simulations/functions/dRSA_PCR.m:297-325`

### Why PCA is used at all

Once nuisance predictors are added, the design matrix can become wide and highly collinear. PCA is used to compress or rotate this predictor set into a lower-dimensional score space before regression.

That is the "PCR" part:

- Principal Component Regression

The code implements it in:

- `simulations/functions/dRSA_PCR.m:465-513`

## Step-By-Step: Exact Legacy PCR Flow In This Repository

This section follows the actual execution path used by `PE_simulation_diff_1Dot.m`.

### Step 1: the one-dot script decides which dRSA branch to run

The top-level script allows:

- `dRSAtypeToRun = 'corr'`
- `dRSAtypeToRun = 'PCR'`

See:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:211-231`

For PCR runs, it sets the default strategy to the legacy path:

- `params.PCR.RegressStrategy = 'baseline_pcr_border'`

at:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:140-173`

### Step 2: the script builds the model list and PCR control structure

For the standard one-dot run, the script builds model arrays such as:

- observed position dot1
- observed direction dot1
- optionally predicted position/direction in the deviant condition

See:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:644-668`

Then it fills the core dRSA parameters:

- `params.fs`
- `params.AverageTime`
- `params.modelToTest`
- `params.Var`
- `params.modelDistMeasure`
- `params.neuralDistMeasure`
- `params.dRSAtype`

and the PCR flags:

- `params.PCR.AdditionalPCA`
- `params.PCR.RegressAutocor`
- `params.PCR.RessModel`

in:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:656-680`

Practical defaults in this script:

- `AdditionalPCA = 1`
- `RegressAutocor = 1`
- `RessModel = 1`

That is a key practical point:

- the function default for `AdditionalPCA` is `0` in `dRSA_coreFunction`,
- but this script usually overrides it to `1`.

References:

- `simulations/functions/dRSA_coreFunction.m:86-99`
- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:673-680`

### Step 3: the script defines which other models should be regressed out

For each tested model, the nuisance list is simply every other model:

- `params.modeltoRegressout{iModel} = setdiff(1:numel(model), iModel)`

See:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:683-685`

So if the run contains four models, each tested model regresses out the other three.

### Step 4: if PCR is selected, the script computes `Autocorrborder`

The one-dot script only computes a real border when the strategy is the legacy one:

- `baseline_pcr_border` -> call `dRSA_border(...)`
- otherwise -> store `NaN` borders because the newer strategies do not use border-based autocorr regression

See:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:688-697`

### Step 5: `dRSA_border` estimates the border from model autocorrelation

`dRSA_border` works like this:

1. for each tested model and each subsampling iteration, it temporarily sets:
   - `params.dRSAtype = 'corr'`
2. it calls `dRSA_coreFunction(model{iModel}, {model{iModel}}, params, ...)`
3. `dRSA_coreFunction` detects the special "model against itself" case and computes one RDM only
4. `dRSA_corr` then produces a self-dRSA matrix, which in practice is the model's time-by-time autocorrelation structure
5. those matrices are averaged across iterations
6. `dRSA_average` averages across the main diagonal neighborhood
7. the border is chosen where the averaged autocorrelation drops below a threshold

Code path:

- `simulations/functions/dRSA_border.m:69-90`
- `simulations/functions/dRSA_coreFunction.m:176-188`
- `simulations/functions/dRSA_coreFunction.m:223-232`
- `simulations/functions/dRSA_corr.m:12-22`
- `simulations/functions/dRSA_average.m:12-38`

#### Border threshold detail that matters in practice

The border threshold is not `params.Var` directly.

`dRSA_border` computes:

- `regvalue = sqrt(params.Var)`

then compares autocorrelation against:

- `regvalue * peak`

See:

- `simulations/functions/dRSA_border.m:55-58`
- `simulations/functions/dRSA_border.m:108-119`

So with the one-dot script default:

- `params.Var = 0.1`

the amplitude threshold is about:

- `sqrt(0.1) = 0.316`

times the peak autocorrelation, not `0.1` times the peak.

### Step 6: the main script calls `dRSA_coreFunction`

For each iteration, the script passes:

- `Y`
- `model`
- `params`
- `CurrSubsamples`
- `Autocorrborder`

into:

- `dRSA_coreFunction`

See:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:717-742`
- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:744-779`

### Step 7: `dRSA_coreFunction` fills defaults and validates PCR settings

`dRSA_coreFunction` is the dispatcher. It:

1. parses optional inputs
2. fills missing defaults
3. checks that PCR is valid
4. decides which RDMs must be computed
5. computes neural/model RDMs
6. dispatches to either `dRSA_corr` or `dRSA_PCR`

Relevant lines:

- input parsing: `simulations/functions/dRSA_coreFunction.m:65-77`
- default filling: `simulations/functions/dRSA_coreFunction.m:79-115`
- PCR validation: `simulations/functions/dRSA_coreFunction.m:118-134`
- decide which model RDMs are needed: `simulations/functions/dRSA_coreFunction.m:138-156`
- compute RDMs: `simulations/functions/dRSA_coreFunction.m:159-216`
- final dispatch: `simulations/functions/dRSA_coreFunction.m:220-240`

Important PCR defaults defined here:

- `params.PCR.Method = 'ExplainedVar'`
- `params.PCR.Methodfactor = 0.1`
- `params.PCR.AdditionalPCA = 0`
- `params.PCR.RegressAutocor = 1`
- `params.PCR.RessModel = 1`
- `params.PCR.RegressStrategy = 'baseline_pcr_border'`

See:

- `simulations/functions/dRSA_coreFunction.m:86-98`

### Step 8: `dRSA_coreFunction` converts raw signals into RDMs

For both `corr` and `PCR`, the code first computes:

- `nRDM` from the current neural/data signal `Y`
- `mRDM{iModel}` for all tested models and any models that will be regressed out

See:

- `simulations/functions/dRSA_coreFunction.m:140-156`
- `simulations/functions/dRSA_coreFunction.m:189-213`

This is why PCR can regress out nuisance models: their RDMs are explicitly computed even if they are not the currently tested model.

### Step 9: `dRSA_PCR` starts the regression loop

Inside `dRSA_PCR`, the code loops over:

1. each tested model
2. each model timepoint `iT`

and fills one row of the final dRSA matrix at a time.

See:

- `simulations/functions/dRSA_PCR.m:67-91`
- `simulations/functions/dRSA_PCR.m:92-181`

At each row:

- `XTest = mRDM_work{iModel}(:, iT)`

is the tested predictor.

Reference:

- `simulations/functions/dRSA_PCR.m:101-105`

### Step 10: PCR assembles nuisance regressors

For the legacy branch, nuisance regressors come from two places.

#### 10a. Other-model nuisance regressors

The local time window around `iT` is:

- `idxWindow = local_window_indices(iT, nTimePoints, tRadius)`

See:

- `simulations/functions/dRSA_PCR.m:113-116`
- `simulations/functions/dRSA_PCR.m:275-279`

Then `local_build_other_model_regressors`:

1. takes each nuisance model listed in `models2regressout`
2. slices that model's RDM over `idxWindow`
3. optionally reduces that block with `AdditionalPCA`
4. concatenates all nuisance blocks column-wise

See:

- `simulations/functions/dRSA_PCR.m:113-116`
- `simulations/functions/dRSA_PCR.m:281-295`

#### 10b. Own-model autocorrelation nuisance regressors

For `baseline_pcr_border`, the code does **not** use every lag in the local window.

Instead, it uses the border estimate to take columns from the outer left and outer right parts of the window:

- left side: from `iT - windowSamples` up to `iT - windowSamples + regborder`
- right side: from `iT + windowSamples - regborder` up to `iT + windowSamples`

See:

- `simulations/functions/dRSA_PCR.m:107-111`
- `simulations/functions/dRSA_PCR.m:297-325`

This is a very repository-specific implementation detail. The legacy autocorr control is therefore:

- border-based,
- asymmetric at the edges if indices fall outside the valid range,
- and based on outer-window lag slices, not on "all nonzero lags".

### Step 11: optional Additional PCA is applied to nuisance blocks

Whenever `params.PCR.AdditionalPCA == 1`, each nuisance block is reduced separately before the final PCR step.

That logic lives in:

- `simulations/functions/dRSA_PCR.m:419-443`

What it does:

1. run `pca(blockIn)`
2. keep components whose explained variance is `> 0.1`
3. keep at least one component

This threshold is hard-coded in `local_reduce_predictor_block`; it does **not** use `params.PCR.Methodfactor`.

That is important because there are two different PCA decisions in this pipeline:

- `AdditionalPCA` for nuisance-block reduction, with a fixed `0.1` threshold
- `params.PCR.Method` / `Methodfactor` for the final PCR solve

### Step 12: nuisance blocks may be normalized again

In the legacy path, after nuisance blocks are concatenated, the code normalizes `xRegressout` only when `AdditionalPCA` is enabled.

See:

- `simulations/functions/dRSA_PCR.m:118-124`
- `simulations/functions/dRSA_PCR.m:445-463`

The normalization method is:

- `Standardize`, or
- `Rescale`

from `params.normalization`.

This preserves the old behavior.

### Step 13: PCR builds the full design matrix

The predictor matrix is:

- `Xx = [XTest, xRegressout]`

where:

- column 1 is always the tested predictor,
- all later columns are nuisance regressors.

See:

- `simulations/functions/dRSA_PCR.m:122-124`

The response matrix is the full neural RDM across all neural timepoints:

- `nRDM_work`

passed into:

- `local_run_original_pcr(Xx, nRDM_work, params)`

at:

- `simulations/functions/dRSA_PCR.m:122-125`

### Step 14: the legacy PCR solve happens in `local_run_original_pcr`

This helper is the heart of the historical method:

- `simulations/functions/dRSA_PCR.m:465-513`

Its sub-steps are:

1. choose how many PCA components to keep
2. regress `nRDM` on `PCAScores`
3. project the solution back to predictor space
4. keep only the coefficient of predictor 1

#### Component selection choices

The code supports four PCR methods:

- `FixedComp`
- `MinCompPCR`
- `ExplainedVar`
- `CumulativeVar`

See:

- `simulations/functions/dRSA_coreFunction.m:32-45`
- `simulations/functions/dRSA_PCR.m:475-505`

In this repository's current defaults, if nothing else is specified, the active choice is:

- `params.PCR.Method = 'ExplainedVar'`
- `params.PCR.Methodfactor = 0.1`

from:

- `simulations/functions/dRSA_coreFunction.m:87-88`

#### The regression itself

After PCA, the code does:

- `betaPCR = PCAScores \\ nRDM`

which is a least-squares solve with:

- predictors = PCA scores
- responses = all neural timepoints at once

See:

- `simulations/functions/dRSA_PCR.m:508-509`

Then it back-projects:

- `temporarydRSA = PCALoadings * betaPCR`

and returns only:

- `temporarydRSA(1, :)`

because row 1 corresponds to `XTest`.

See:

- `simulations/functions/dRSA_PCR.m:509-512`

That returned row becomes:

- `dRSAmat(iT, :, iMod)`

at:

- `simulations/functions/dRSA_PCR.m:180-181`

So the final output is still a `model time x neural time x model` tensor, but its entries are regression weights, not correlations.

#### A subtle implementation detail for `ExplainedVar` and `CumulativeVar`

For `ExplainedVar` and `CumulativeVar`, the current code reduces `PCALoadings` to its first row before back-projection:

- `PCALoadings = PCALoadings(1, 1:imax)`

See:

- `simulations/functions/dRSA_PCR.m:483-503`

That is consistent with the fact that the function only keeps the first predictor's coefficient anyway, but it is worth knowing if you are tracing the exact current implementation.

## How Autocorrelation Is Handled

### In `corr`

Normal `corr` analysis does not explicitly remove autocorrelation. It simply correlates model and neural RDM columns.

Reference:

- `simulations/functions/dRSA_corr.m:12-22`

The only special autocorrelation-related use of `corr` is inside `dRSA_border`, where the model is correlated with itself to estimate the border:

- `simulations/functions/dRSA_border.m:69-90`

### In legacy `PCR`

Autocorrelation handling is explicit and model-specific:

1. estimate a border from each model's self-dRSA/autocorrelation profile
2. use that border to pick own-model nuisance regressors
3. include those nuisance columns in the PCR design matrix

Code path:

- border estimation: `simulations/functions/dRSA_border.m:63-119`
- border-only used for legacy strategy: `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:688-697`
- nuisance assembly: `simulations/functions/dRSA_PCR.m:107-125`
- border-based autocorr helper: `simulations/functions/dRSA_PCR.m:297-325`

### Practical meaning of the border

`Autocorrborder(iModel)` controls how much of the outer local window is used as own-model nuisance regressors.

It is not:

- a time in seconds,
- a correlation value,
- or a count of all lags outside zero.

It is a model-specific sample index derived from the averaged autocorrelation profile and then used to cut left/right lag slices in the baseline PCR implementation.

## Where Strategy, PCA, And Ridge Choices Happen

Even though this document focuses on legacy PCR, the current codebase supports multiple backends. The branch point is worth knowing.

### Strategy choice

`dRSA_coreFunction` validates the strategy name here:

- `simulations/functions/dRSA_coreFunction.m:123-134`

`dRSA_PCR` then reads it and switches behavior here:

- `simulations/functions/dRSA_PCR.m:83-90`
- `simulations/functions/dRSA_PCR.m:106-178`

For legacy PCR, the relevant branch is:

- `case 'baseline_pcr_border'`

at:

- `simulations/functions/dRSA_PCR.m:106-125`

### Final PCR component-selection choice

The legacy PCR method itself is chosen in:

- `simulations/functions/dRSA_PCR.m:475-505`

using:

- `params.PCR.Method`
- `params.PCR.Methodfactor`

### Ridge-related choices

The same file also contains the newer ridge-based branches:

- `ridge_full_autocorr`
- `ridge_tapered_autocorr`
- `ar1_prewhite_ridge`

See:

- `simulations/functions/dRSA_PCR.m:126-175`
- `simulations/functions/dRSA_PCR.m:515-586`

Those parameters are present in the current code, but they do not change the behavior of legacy `baseline_pcr_border`.

## Key Parameters In Practice

### Parameters that matter directly for legacy `baseline_pcr_border`

| Parameter | What it does in practice | Where used |
| --- | --- | --- |
| `params.dRSAtype` | Chooses direct correlation vs PCR dispatch | `simulations/functions/dRSA_coreFunction.m:223-240` |
| `params.modelToTest` | Which model indices will get a full dRSA matrix | `simulations/functions/dRSA_coreFunction.m:147-155`, `simulations/functions/dRSA_PCR.m:69-75` |
| `params.modeltoRegressout` | For each tested model, which other models become nuisance regressors | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:683-685`, `simulations/functions/dRSA_PCR.m:98-99` |
| `params.fs` | Converts seconds to samples for windowing | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:658-662`, `simulations/functions/dRSA_PCR.m:71-72` |
| `params.AverageTime` | Sets the half-window radius used for local nuisance selection and for border estimation | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:660-662`, `simulations/functions/dRSA_PCR.m:71-72`, `simulations/functions/dRSA_border.m:15-16` |
| `params.Var` | Sets the autocorr-border threshold indirectly via `sqrt(params.Var)` | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:663-664`, `simulations/functions/dRSA_border.m:55-58`, `simulations/functions/dRSA_border.m:116-119` |
| `params.modelDistMeasure` | Distance used to build each model RDM | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:665-668`, `simulations/functions/dRSA_coreFunction.m:207-212` |
| `params.neuralDistMeasure` | Distance used to build the neural/data RDM | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:669-670`, `simulations/functions/dRSA_coreFunction.m:195-197` |
| `params.normalization` | Chooses `Standardize` vs `Rescale` for PCR RDM normalization and for nuisance-block normalization after `AdditionalPCA` | `simulations/functions/dRSA_computeRDM.m:88-107`, `simulations/functions/dRSA_PCR.m:445-463` |
| `params.PCR.RegressStrategy` | Chooses legacy border-based PCR vs newer alternatives | `simulations/functions/dRSA_coreFunction.m:123-130`, `simulations/functions/dRSA_PCR.m:187-195` |
| `params.PCR.RegressAutocor` | Enables/disables own-model autocorr nuisance regressors | `simulations/functions/dRSA_PCR.m:306-313` |
| `params.PCR.RessModel` | Enables/disables nuisance regressors from other models | `simulations/functions/dRSA_PCR.m:281-295` |
| `params.PCR.AdditionalPCA` | Reduces each nuisance block separately before the final PCR solve | `simulations/functions/dRSA_PCR.m:419-443` |
| `params.PCR.Method` | Chooses how the final PCR PCA dimensionality is selected | `simulations/functions/dRSA_PCR.m:475-505` |
| `params.PCR.Methodfactor` | Threshold/count/fraction used by `Method` | `simulations/functions/dRSA_PCR.m:475-505` |
| `Autocorrborder` | Supplies the model-specific border used by legacy autocorr nuisance selection | `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:688-697`, `simulations/functions/dRSA_PCR.m:297-325` |

### Parameters present in the current code but not used by legacy `baseline_pcr_border`

These matter only for the newer ridge-based strategies:

- `params.PCR.RidgeLambdaFactor`
- `params.PCR.TaperSigmaFactor`
- `params.PCR.StandardizePredictors`
- `params.PCR.AR1Clip`
- `params.PCR.AutoPenaltyStrength`
- `params.PCR.NormalizeAutoPenaltyPerRow`

Defined/validated in:

- `simulations/functions/dRSA_coreFunction.m:92-98`
- `simulations/functions/dRSA_coreFunction.m:123-133`

Used only in ridge/AR(1) branches of:

- `simulations/functions/dRSA_PCR.m:126-175`
- `simulations/functions/dRSA_PCR.m:197-209`
- `simulations/functions/dRSA_PCR.m:342-417`
- `simulations/functions/dRSA_PCR.m:515-586`

### Two naming quirks worth knowing

1. The code uses the field name:
   - `params.PCR.RessModel`
   - not `RegressModel`
2. A comment in `dRSA_coreFunction` says `AdditionalPCR`, but the actual field used everywhere is:
   - `params.PCR.AdditionalPCA`

References:

- `simulations/functions/dRSA_coreFunction.m:27-30`
- `simulations/functions/dRSA_coreFunction.m:87-90`
- `simulations/functions/dRSA_PCR.m:419-443`

## Follow-The-Code Map

If you want to trace one legacy PCR run in order, this is the shortest reliable path.

1. Start at the script entry and defaults:
   - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:140-173`
   - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:211-231`
2. See where models and parameters are assembled:
   - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:644-685`
3. See where `Autocorrborder` is computed for the legacy branch:
   - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:688-697`
4. Read how `dRSA_border` computes that border:
   - `simulations/functions/dRSA_border.m:63-119`
5. See the actual dRSA call into the core function:
   - `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:717-742`
6. In `dRSA_coreFunction`, trace defaults, RDM creation, and dispatch:
   - `simulations/functions/dRSA_coreFunction.m:79-115`
   - `simulations/functions/dRSA_coreFunction.m:138-216`
   - `simulations/functions/dRSA_coreFunction.m:220-240`
7. In `dRSA_computeRDM`, see how time-resolved RDM columns are formed:
   - `simulations/functions/dRSA_computeRDM.m:19-30`
   - `simulations/functions/dRSA_computeRDM.m:86-107`
8. In `dRSA_PCR`, trace the baseline branch:
   - `simulations/functions/dRSA_PCR.m:92-125`
9. Then inspect the two nuisance builders:
   - other models: `simulations/functions/dRSA_PCR.m:281-295`
   - autocorr border terms: `simulations/functions/dRSA_PCR.m:297-325`
10. Finally inspect the legacy PCR math:
   - `simulations/functions/dRSA_PCR.m:465-513`

## Minimal Runnable Examples

These examples follow the style used by the project scripts: define a few workspace variables, then `run(...)` the script. They use the absolute MATLAB binary requested in the repository instructions.

### Example 1: `corr`

```bash
/Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'); clear params; participantNumber=73; dRSAtypeToRun='corr'; suppressDispText=0; run('simulations/scripts/one-dot/PE_simulation_diff_1Dot.m');"
```

If MATLAB R2020a needs Intel mode on Apple Silicon:

```bash
arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'); clear params; participantNumber=73; dRSAtypeToRun='corr'; suppressDispText=0; run('simulations/scripts/one-dot/PE_simulation_diff_1Dot.m');"
```

### Example 2: legacy `PCR`

This example makes the legacy PCR settings explicit instead of relying only on defaults.

```bash
/Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'); clear params; params=struct(); params.PCR=struct(); params.PCR.RegressStrategy='baseline_pcr_border'; params.PCR.AdditionalPCA=1; params.PCR.RegressAutocor=1; params.PCR.RessModel=1; participantNumber=73; dRSAtypeToRun='PCR'; suppressDispText=0; run('simulations/scripts/one-dot/PE_simulation_diff_1Dot.m');"
```

Apple Silicon fallback:

```bash
arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'); clear params; params=struct(); params.PCR=struct(); params.PCR.RegressStrategy='baseline_pcr_border'; params.PCR.AdditionalPCA=1; params.PCR.RegressAutocor=1; params.PCR.RessModel=1; participantNumber=73; dRSAtypeToRun='PCR'; suppressDispText=0; run('simulations/scripts/one-dot/PE_simulation_diff_1Dot.m');"
```

## Troubleshooting

### Common artifacts or misinterpretations

#### 1. "PCR values look smaller than corr, so PCR failed."

Not necessarily. The scales are different:

- `corr` returns correlation coefficients,
- `PCR` returns regression weights for the tested predictor.

So direct magnitude comparison is often misleading.

Relevant code:

- `simulations/functions/dRSA_corr.m:19-21`
- `simulations/functions/dRSA_PCR.m:508-512`

#### 2. "I still see strong near-diagonal structure in PCR."

That can happen because the legacy strategy does not regress all neighboring lags symmetrically. It uses the border-based outer-window slices from `local_build_baseline_autocorr_regressors`.

Check first:

- `Autocorrborder`
- `params.Var`
- `params.AverageTime`
- whether `params.PCR.RegressAutocor` is actually `1`

Relevant code:

- `simulations/functions/dRSA_border.m:55-58`
- `simulations/functions/dRSA_border.m:108-119`
- `simulations/functions/dRSA_PCR.m:297-325`

#### 3. "PCR looks too much like corr."

Usually check these first:

- `params.PCR.RessModel`
- `params.PCR.RegressAutocor`
- `params.modeltoRegressout`
- whether the run really used `dRSAtypeToRun = 'PCR'`
- whether `Autocorrborder` was computed instead of left empty/`NaN`

Relevant code:

- `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m:671-697`
- `simulations/functions/dRSA_coreFunction.m:233-239`

#### 4. "Borders are odd or NaN."

If the averaged autocorrelation never drops below the threshold inside the inspected lag range, border detection can become unstable or return `NaN`.

Check first:

- `params.Var`
- the model autocorr plots from `dRSA_border`
- whether the chosen model has meaningful temporal variation

Relevant code:

- `simulations/functions/dRSA_border.m:116-125`
- `simulations/functions/dRSA_border.m:151-221`

#### 5. "The nuisance regression feels too aggressive or too weak."

The first settings to inspect are:

- `params.AverageTime`
- `params.Var`
- `params.PCR.AdditionalPCA`
- `params.PCR.Method`
- `params.PCR.Methodfactor`

Why:

- `AverageTime` changes the local window size,
- `Var` changes the border threshold,
- `AdditionalPCA` changes nuisance-block compression,
- `Method` and `Methodfactor` change the final PCR dimensionality.

#### 6. "I changed ridge settings and nothing happened."

If the strategy is still:

- `params.PCR.RegressStrategy = 'baseline_pcr_border'`

then ridge-only settings do nothing.

Check first:

- `params.PCR.RegressStrategy`

Relevant code:

- `simulations/functions/dRSA_PCR.m:106-125`
- `simulations/functions/dRSA_PCR.m:126-175`

## Bottom Line

In this repository, legacy `dRSAtype = 'PCR'` is best understood as:

- time-resolved RDM regression,
- with target-model coefficient extraction,
- after nuisance control for other models and own-model temporal autocorrelation,
- using PCA to stabilize a highly collinear predictor set.

If you want to understand the exact implementation, the most important lines are:

- `simulations/functions/dRSA_coreFunction.m:220-240`
- `simulations/functions/dRSA_PCR.m:92-125`
- `simulations/functions/dRSA_PCR.m:297-325`
- `simulations/functions/dRSA_PCR.m:465-513`
- `simulations/functions/dRSA_border.m:63-119`

