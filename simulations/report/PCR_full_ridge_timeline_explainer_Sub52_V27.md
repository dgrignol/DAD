# PCR full ridge (`ridge_full_autocorr`) and timeline selection in the one-dot occlusion simulation

This note explains, in detail, what the `PCR` mode with `params.PCR.RegressStrategy = 'ridge_full_autocorr'` does in the current sandbox simulation code, and exactly which timeline segment is analyzed (and why).

Context used for this explanation:
- Script: `simulations/scripts/one-dot/PE_simulation_diff_1Dot_occlusion_v1.m`
- PCR backend: `simulations/functions/dRSA_PCR.m`
- Latest sandbox run: `Sub52`, V27 block-resume input, `dRSAtypeToRun='PCR'`, strategy `ridge_full_autocorr`
- Result file inspected:
  `simulations/output/Sub52_oneDot_occlusion/occluded_deviant/pcr/results/dRSA_sub52_occluded_deviant_pcr_oneDot_pcrridgeFullAuto_results.mat`

---

## 1) What `ridge_full_autocorr` is doing conceptually

For each tested model and each model-time sample `iT`, the code builds one regression problem where:
- The first predictor is the tested model at time `iT`.
- Additional nuisance predictors include:
1. Other models around `iT` (local window).
2. The same tested model at other nearby lags (same local window, excluding lag 0).
- The response is the neural RDM across all neural times.

Then it solves a **multi-response ridge regression** and keeps only the coefficient of the **first predictor** (the tested effect). That coefficient row becomes one row in the dRSA matrix.

So this strategy tries to estimate the unique contribution of the model at `iT` while controlling for:
- overlap with other models,
- temporal autocorrelation of the tested model itself.

---

## 2) Exact design matrix in `ridge_full_autocorr`

In `dRSA_PCR.m`:
- `XTest = mRDM_work{iModel}(:, iT)` (tested predictor; first column)
- `idxWindow = local_window_indices(iT, nTimePoints, tRadius)`
- `xModelRegressout` is built from nuisance models on that same window
- `xAutocorrelation` is built from the tested model on that same window, excluding `iT`

`local_build_full_autocorr_regressors(...)` does:
- Start from full local window around `iT`
- Remove `iT` itself
- Keep all other lag columns as autocorr nuisance regressors

So, for `ridge_full_autocorr`, the matrix is:
- `X = [XTest, xModelRegressout, xAutocorrelation]`

Important detail:
- `AdditionalPCA` may reduce **other-model** nuisance blocks (`xModelRegressout`) through `local_reduce_predictor_block`.
- The own-lag autocorr block in this strategy is not PCA-reduced; it is included explicitly as lag columns.

---

## 3) How the ridge penalty is applied

`local_build_penalty_weights(...)` and `local_ridge_first_predictor(...)` implement weighted ridge:
- Tested predictor (first column): penalty weight `0` (unpenalized)
- Other-model nuisance columns: penalty weight `1`
- Full-autocorr lag nuisance columns: penalty weight `1` for this strategy (uniform)

Penalty scale:
- `lambda = RidgeLambdaFactor * nObs`
- In this run: `RidgeLambdaFactor = 0.06`

Solve:
- `beta = (X'X + lambda * D)^(-1) X'Y`
- `D = diag(penaltyWeights)` with first diagonal forced to `0`
- If enabled (default), predictors are z-scored before solving, then beta is mapped back to original scale.

Finally:
- The algorithm stores only `beta(firstPredictor, :)`.
- That row is the dRSA row for the current model-time `iT` against all neural times.

---

## 4) Two different "time windows" are involved (critical distinction)

There are two timeline concepts, and it is easy to mix them up.

### A) Global analyzed timeline (matrix axes)
This is the actual data window entering RDMs and shown on the dRSA matrix axes.

It is controlled in `PE_simulation_diff_1Dot_occlusion_v1.m` by:
- whether streams are cut before dRSA,
- and by trial-locked subsampling options (`PreTrigger=0`, `PostTrigger=trialLen-1`).

Result: every dRSA matrix time axis spans the selected trial window from its first included frame to its last included frame.

### B) Local nuisance-regression window inside each matrix row
Inside PCR, for each row `iT`, nuisance columns are collected from:
- `idxWindow = (iT - tRadius) : (iT + tRadius)`, clamped to valid range.

So the global timeline says "what samples are in the analysis", while `tRadius` says "how wide the local nuisance context is when estimating each row".

---

## 5) Which timeline is used in this current sandbox run

From inspected repro metadata (Sub52 V27 run):
- `cutPostDev = 0` (false)
- `trialLen = 320`
- `occlusionEventFrames.fullDisappearanceFrame = 130`
- `occlusionEventFrames.reappearanceStartFrame = 191`
- `peCutFrame = 191`
- `peCutSource = 'peBranchCut'`
- `peTrialLenBeforeCut = 320`
- `peTrialLenAfterCut = 130`
- `fs = 120 Hz`

So, in this run:

### Standard (non-PE) branch
- Because `cutPostDev=false`, standard observed streams remain full-length.
- Global timeline analyzed: frames `1:320`.
- Guide times are therefore absolute-trial times:
  - full disappearance at frame 130 -> `1.075 s`
  - reappearance start at frame 191 -> `1.5833 s`

### PE branch (deviant-only)
- PE branch is forcibly post-reappearance even when `cutPostDev=false`.
- The script applies an internal PE cut at frame `191`.
- Global timeline analyzed for PE dRSA: frames `191:320` only (130 samples).
- PE axis time zero corresponds to absolute frame 191.
- Therefore:
  - reappearance starts at `0 s` in PE matrices,
  - full disappearance frame 130 is out of PE range and appears as `NaN` guide time.

This behavior is intentional and is explicitly enforced in the script.

---

## 6) Why this timeline choice exists

The PE signal is defined as:
- `PE = observed - predicted`

In this experiment logic, PE of interest is the post-occlusion/post-reappearance mismatch phase. If pre-deviance samples were mixed into PE matrices, the PE estimate would include long periods where prediction error is not the target phenomenon, diluting interpretability.

So the code enforces:
- Standard matrices can remain full-trial (for continuity and comparability) when `cutPostDev=false`.
- PE matrices always use post-reappearance data so that PE models and PE neural targets are computed on the specific segment where deviance-driven mismatch is expected.

This is also robust to experiment-version changes because the cut frame is metadata-driven:
1. `reappearance_start_frame`
2. fallback `occlusion_end_frame`
3. only if metadata missing: fractional fallback (`deviantOnset`, default `0.5`)

Given your note that deviance starts earlier in the new experiment version, this metadata-first rule is exactly why the code now tracks event timing correctly instead of relying on a fixed fraction.

---

## 7) How large the local PCR window is in this run

`dRSA_PCR` sets:
- `tRadius = round(params.AverageTime * params.fs)`

### Standard branch
- `params.AverageTime = floor((trialLen-1)/2) / fs`
- with `trialLen=320`, `fs=120`:
  - `AverageTime = 159/120 = 1.325 s`
  - `tRadius = 159`
- So each row can use almost the whole trial as local context (except boundary clipping).

### PE branch
- PE branch sets its own `AverageTime` from `peTrialLenAfterCut`.
- with `peTrialLenAfterCut=130`, `fs=120`:
  - `AverageTime = 64/120 = 0.5333 s`
  - `tRadius = 64`
- Again, this is effectively near-full coverage of the PE timeline for each row.

Interpretation: with `ridge_full_autocorr`, the tested effect at time `iT` is adjusted against many own-lag nuisance terms from almost the entire currently analyzed window.

---

## 8) Practical implication for reading the output matrices

For `occluded_deviant` in this run:
- Standard matrices (`base`/`predicted`) are on full trial time.
- PE matrices are on post-reappearance time only.

So do not compare PE matrix x/y-axis zero directly to standard matrix x/y-axis zero as the same absolute frame. In PE plots, `0 s` means "first frame of reappearance window" (absolute frame 191 in this run), not trial onset.

---

## 9) Short summary

`PCR` + `ridge_full_autocorr` here means:
- Build explicit nuisance-rich design matrix per model-time row.
- Use ridge regularization with unpenalized tested predictor and uniformly penalized nuisance predictors.
- Keep only tested predictor beta across neural times.

Timeline used:
- Standard dRSA: full trial (when `cutPostDev=false`).
- PE dRSA: mandatory post-reappearance segment (metadata-driven cut), because PE is defined and interpreted specifically in that deviance-relevant window.

