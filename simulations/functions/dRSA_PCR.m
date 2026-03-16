function dRSAmat = dRSA_PCR(mRDM, nRDM, Autocorrborder, params)

%% DRSA_PCR Compute dRSA matrices using baseline PCR or ridge-based alternatives.
% Usage examples:
%   % 1) Original baseline behavior (default strategy).
%   params.PCR.RegressStrategy = 'baseline_pcr_border';
%   dRSAmat = dRSA_PCR(mRDM, nRDM, Autocorrborder, params);
%
%   % 2) Ridge replacement with explicit full-window autocorr regressors.
%   params.PCR.RegressStrategy = 'ridge_full_autocorr';
%   params.PCR.RidgeLambdaFactor = 0.06;
%   dRSAmat = dRSA_PCR(mRDM, nRDM, [], params);
%
%   % 3) Ridge with tapered autocorr regressors.
%   params.PCR.RegressStrategy = 'ridge_tapered_autocorr';
%   params.PCR.TaperSigmaFactor = 0.33;
%   dRSAmat = dRSA_PCR(mRDM, nRDM, [], params);
%
%   % 4) AR(1) prewhitening + ridge regression.
%   params.PCR.RegressStrategy = 'ar1_prewhite_ridge';
%   params.PCR.AR1Clip = 0.99;
%   params.PCR.RidgeLambdaFactor = 0.04;
%   dRSAmat = dRSA_PCR(mRDM, nRDM, [], params);
%
% Inputs:
%   mRDM: 1 x nModels cell, each model RDM is nFeatures x nTimepoints.
%   nRDM: neural/data RDM, nFeatures x nTimepoints.
%   Autocorrborder: nModels x 1 vector with regression border samples
%       (required only for baseline_pcr_border when RegressAutocor = 1).
%   params: struct with standard dRSA PCR fields plus strategy-specific fields:
%       .modelToTest
%       .modeltoRegressout
%       .AverageTime
%       .fs
%       .normalization (optional)
%       .PCR.RegressAutocor
%       .PCR.RessModel
%       .PCR.AdditionalPCA
%       .PCR.Method
%       .PCR.Methodfactor
%       .PCR.RegressStrategy
%       .PCR.RidgeLambdaFactor
%       .PCR.TaperSigmaFactor
%       .PCR.StandardizePredictors
%       .PCR.AR1Clip
%       .PCR.AutoPenaltyStrength
%       .PCR.NormalizeAutoPenaltyPerRow
%
% Output:
%   dRSAmat: nTimepoints x nTimepoints x nModels matrix.
%
% Strategy overview:
%   baseline_pcr_border:
%       Original implementation (PCR + autocorr border + least-squares in PCA space).
%   ridge_full_autocorr:
%       Ridge solve with target predictor + other-model regressors + full-window
%       own-model lag regressors (excluding lag 0).
%   ridge_tapered_autocorr:
%       Same as ridge_full_autocorr, but own-model lag regressors receive
%       lag-dependent ridge penalties (near lags penalized less, far lags
%       penalized more) to avoid hard discontinuities while preserving
%       predictor standardization.
%   ar1_prewhite_ridge:
%       AR(1)-prewhiten model and neural RDM time series first, then solve ridge
%       with target predictor + other-model regressors (no own-lag regressors).

%% Preparation
% Data flow: params -> strategy + regularization settings -> reusable loop constants.
iMod = 0;
nModels = length(params.modelToTest);
nTimePoints = size(nRDM, 2);
tRadius = max(0, round(params.AverageTime * params.fs));

% Preallocate final output tensor.
dRSAmat = zeros(nTimePoints, nTimePoints, nModels);

% Keep legacy default for normalization when AdditionalPCA is enabled.
if isfield(params.PCR, 'AdditionalPCA') && isequal(params.PCR.AdditionalPCA, 1) && ...
        ~isfield(params, 'normalization')
    params.normalization = 'Standardize';
end

% Resolve strategy and ridge configuration defaults once at function entry.
strategy = local_get_pcr_strategy(params);
ridgeLambdaFactor = local_get_ridge_lambda(params, strategy);
standardizePredictors = local_get_bool_field(params.PCR, 'StandardizePredictors', 1);

% Prepare strategy-specific RDM views once.
% For AR(1) strategy, this prewhitens all RDMs before any regression steps.
[nRDM_work, mRDM_work] = local_prepare_strategy_rdms(nRDM, mRDM, params, strategy);

%% Main dRSA loop
% Data flow: model index -> model-time predictors -> regression weights per neural time.
for iModel = params.modelToTest
    iMod = iMod + 1;
    fprintf('\n Model: %04d ', iModel);

    % For each tested model, we track which other model indices are nuisance regressors.
    models2regressout = params.modeltoRegressout{iMod};

    % Loop all model-time samples (rows of the dRSA matrix).
    for iT = 1:nTimePoints
        % Target predictor at current model-time sample (first column in design matrix).
        XTest = mRDM_work{iModel}(:, iT);

        switch strategy
            case 'baseline_pcr_border'
                %% Baseline: original PCR + border logic
                % 1) Build own-model autocorr nuisance block outside border.
                xAutocorrelation = local_build_baseline_autocorr_regressors( ...
                    mRDM_work{iModel}, iT, nTimePoints, params, Autocorrborder, iModel);

                % 2) Build other-model nuisance block across local averaging window.
                idxWindow = local_window_indices(iT, nTimePoints, tRadius);
                xModelRegressout = local_build_other_model_regressors( ...
                    mRDM_work, models2regressout, idxWindow, params);

                % 3) Combine nuisance blocks and normalize exactly as in baseline path.
                xRegressout = [xModelRegressout, xAutocorrelation];
                xRegressout = local_normalize_if_requested(xRegressout, params);

                % 4) Build full predictor matrix and run original PCR routine.
                Xx = [XTest, xRegressout];
                betaFirst = local_run_original_pcr(Xx, nRDM_work, params);

            case {'ridge_full_autocorr', 'ridge_tapered_autocorr', 'ar1_prewhite_ridge'}
                %% Ridge alternatives: explicit design matrix + ridge solve
                % 1) Shared local window around iT for nuisance predictors.
                idxWindow = local_window_indices(iT, nTimePoints, tRadius);

                % 2) Other-model nuisance predictors (optional via params.PCR.RessModel).
                if params.PCR.RessModel
                    xModelRegressout = local_build_other_model_regressors( ...
                        mRDM_work, models2regressout, idxWindow, params);
                else
                    xModelRegressout = [];
                end

                % 3) Own-model autocorr predictors differ by strategy.
                autoPenaltyWeights = [];
                if params.PCR.RegressAutocor
                    switch strategy
                        case 'ridge_full_autocorr'
                            [xAutocorrelation, ~] = local_build_full_autocorr_regressors( ...
                                mRDM_work{iModel}, idxWindow, iT);
                            % Full strategy keeps uniform penalties across autocorr lags.
                            autoPenaltyWeights = ones(1, size(xAutocorrelation, 2));

                        case 'ridge_tapered_autocorr'
                            [xAutocorrelation, autoIdx] = local_build_full_autocorr_regressors( ...
                                mRDM_work{iModel}, idxWindow, iT);
                            % Tapered strategy keeps predictors unchanged and shapes only
                            % the ridge penalty profile as a function of lag distance.
                            autoPenaltyWeights = local_build_tapered_autocorr_penalty_weights( ...
                                autoIdx, iT, tRadius, params);

                        case 'ar1_prewhite_ridge'
                            % AR(1) variant intentionally excludes own-lag nuisance terms.
                            xAutocorrelation = [];
                            autoPenaltyWeights = [];
                    end
                else
                    xAutocorrelation = [];
                    autoPenaltyWeights = [];
                end

                % 4) Assemble design matrix with target predictor in the first column.
                Xx = [XTest, xModelRegressout, xAutocorrelation];
                penaltyWeights = local_build_penalty_weights( ...
                    size(xModelRegressout, 2), size(xAutocorrelation, 2), autoPenaltyWeights);

                % 5) Solve ridge and keep only the coefficient of the first predictor.
                betaFirst = local_ridge_first_predictor( ...
                    Xx, nRDM_work, ridgeLambdaFactor, standardizePredictors, penaltyWeights);

            otherwise
                error('Unhandled PCR strategy: %s', strategy);
        end

        % Fill one dRSA row: model-time iT against all neural times.
        dRSAmat(iT, :, iMod) = betaFirst(:)';
    end
end

end

function strategy = local_get_pcr_strategy(params)
% LOCAL_GET_PCR_STRATEGY Read strategy with backward-compatible default.
if isfield(params, 'PCR') && isfield(params.PCR, 'RegressStrategy') && ...
        ~isempty(params.PCR.RegressStrategy)
    strategy = char(params.PCR.RegressStrategy);
else
    strategy = 'baseline_pcr_border';
end
end

function lambdaFactor = local_get_ridge_lambda(params, strategy)
% LOCAL_GET_RIDGE_LAMBDA Choose lambda factor, with strategy-aware defaults.
if isfield(params.PCR, 'RidgeLambdaFactor') && ~isempty(params.PCR.RidgeLambdaFactor)
    lambdaFactor = params.PCR.RidgeLambdaFactor;
    return;
end
switch strategy
    case 'ar1_prewhite_ridge'
        lambdaFactor = 0.04;
    otherwise
        lambdaFactor = 0.06;
end
end

function out = local_get_bool_field(s, fieldName, defaultValue)
% LOCAL_GET_BOOL_FIELD Read scalar logical-like fields with default fallback.
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    out = logical(s.(fieldName));
else
    out = logical(defaultValue);
end
end

function [nRDM_work, mRDM_work] = local_prepare_strategy_rdms(nRDM, mRDM, params, strategy)
% LOCAL_PREPARE_STRATEGY_RDMS Apply strategy-level preprocessing once.
nRDM_work = nRDM;
mRDM_work = mRDM;

if ~strcmp(strategy, 'ar1_prewhite_ridge')
    return;
end

% For AR(1) strategy, prewhiten both neural and model RDM time series.
ar1Clip = 0.99;
if isfield(params.PCR, 'AR1Clip') && ~isempty(params.PCR.AR1Clip)
    ar1Clip = params.PCR.AR1Clip;
end

nRDM_work = local_prewhiten_rdm(nRDM, ar1Clip);
for iModel = 1:numel(mRDM)
    if isempty(mRDM{iModel})
        continue;
    end
    mRDM_work{iModel} = local_prewhiten_rdm(mRDM{iModel}, ar1Clip);
end
end

function Xw = local_prewhiten_rdm(X, ar1Clip)
% LOCAL_PREWHITEN_RDM Apply pooled AR(1) whitening across time.
%
% Sub-steps:
%   1) Estimate per-feature AR(1) slopes using lagged least-squares.
%   2) Pool with median (robust to outliers).
%   3) Clip pooled slope for numerical stability.
%   4) Filter along time: x_t <- x_t - phi * x_{t-1}.
if size(X, 2) < 2
    Xw = X;
    return;
end

Xprev = X(:, 1:end-1);
Xnext = X(:, 2:end);
phiFeat = sum(Xnext .* Xprev, 2) ./ (sum(Xprev .^ 2, 2) + eps);
phiFeat = phiFeat(isfinite(phiFeat));
if isempty(phiFeat)
    phi = 0;
else
    phi = median(phiFeat);
end

% Clip to avoid explosive filters near unit roots.
phi = max(min(phi, abs(ar1Clip)), -abs(ar1Clip));

Xw = zeros(size(X));
Xw(:, 1) = X(:, 1);
Xw(:, 2:end) = X(:, 2:end) - phi .* X(:, 1:end-1);
end

function idxWindow = local_window_indices(iT, nTimePoints, tRadius)
% LOCAL_WINDOW_INDICES Return valid local window indices around iT.
idxWindow = (iT - tRadius):(iT + tRadius);
idxWindow(idxWindow < 1 | idxWindow > nTimePoints) = [];
end

function xModelRegressout = local_build_other_model_regressors(mRDM, models2regressout, idxWindow, params)
% LOCAL_BUILD_OTHER_MODEL_REGRESSORS Concatenate nuisance blocks from other models.
%
% Sub-steps:
%   1) Slice each nuisance model over the local time window.
%   2) Optionally reduce each slice with AdditionalPCA.
%   3) Concatenate slices column-wise.
xModelRegressout = [];

for iReg = models2regressout
    regressBlock = mRDM{iReg}(:, idxWindow);
    regressBlock = local_reduce_predictor_block(regressBlock, params);
    xModelRegressout = [xModelRegressout, regressBlock]; %#ok<AGROW>
end
end

function xAutocorrelation = local_build_baseline_autocorr_regressors( ...
        mRDM_model, iT, nTimePoints, params, Autocorrborder, iModel)
% LOCAL_BUILD_BASELINE_AUTOCORR_REGRESSORS Build baseline border-based autocorr terms.
%
% This function preserves the original border logic:
%   LeftSide  = [iT-window, ..., iT-window+regborder]
%   RightSide = [iT+window-regborder, ..., iT+window]
%
% If RegressAutocor is disabled or border is NaN/empty, returns [].
if ~params.PCR.RegressAutocor
    xAutocorrelation = [];
    return;
end
if isempty(Autocorrborder) || numel(Autocorrborder) < iModel || isnan(Autocorrborder(iModel))
    xAutocorrelation = [];
    return;
end

regborder = Autocorrborder(iModel);
windowSamples = params.AverageTime * params.fs;
leftSide = (iT - windowSamples):(iT - windowSamples + regborder);
rightSide = (iT + windowSamples - regborder):(iT + windowSamples);

autoIdx = [leftSide, rightSide];
autoIdx(autoIdx < 1 | autoIdx > nTimePoints) = [];

xAutocorrelation = mRDM_model(:, autoIdx);
xAutocorrelation = local_reduce_predictor_block(xAutocorrelation, params);
end

function [xAutocorrelation, autoIdx] = local_build_full_autocorr_regressors(mRDM_model, idxWindow, iT)
% LOCAL_BUILD_FULL_AUTOCORR_REGRESSORS Build full-window own-model lag regressors.
%
% Sub-steps:
%   1) Start from the full local window.
%   2) Remove lag-0 column (current iT) so target predictor remains unique.
%   3) Slice model RDM at the remaining lag indices.
autoIdx = idxWindow(idxWindow ~= iT);
if isempty(autoIdx)
    xAutocorrelation = [];
    return;
end
xAutocorrelation = mRDM_model(:, autoIdx);
end

function autoPenaltyWeights = local_build_tapered_autocorr_penalty_weights(autoIdx, iT, tRadius, params)
% LOCAL_BUILD_TAPERED_AUTOCORR_PENALTY_WEIGHTS Build lag-dependent ridge penalties.
%
% Core idea:
%   Keep predictor columns unchanged, but increase ridge penalty for distant
%   lags so far autocorr terms shrink more strongly than near autocorr terms.
%
% Sub-steps:
%   1) Convert lag offsets to Gaussian proximity weights g(lag) in (0,1].
%   2) Convert proximity to penalty multipliers:
%        penalty = 1 + kappa * (1 - g)
%      where kappa is AutoPenaltyStrength.
%   3) Optionally normalize penalties to mean 1 within each row/timepoint so
%      total autocorr penalty energy stays comparable across edge/center rows.
if isempty(autoIdx)
    autoPenaltyWeights = [];
    return;
end

sigmaFactor = 0.33;
if isfield(params.PCR, 'TaperSigmaFactor') && ~isempty(params.PCR.TaperSigmaFactor)
    sigmaFactor = params.PCR.TaperSigmaFactor;
end
sigma = max(1, sigmaFactor * max(1, tRadius));

kappa = 2;
if isfield(params.PCR, 'AutoPenaltyStrength') && ~isempty(params.PCR.AutoPenaltyStrength)
    kappa = params.PCR.AutoPenaltyStrength;
end
kappa = max(0, kappa);

normalizePerRow = local_get_bool_field(params.PCR, 'NormalizeAutoPenaltyPerRow', 1);

lagAbs = abs(autoIdx - iT);
gaussianProximity = exp(-0.5 * (lagAbs ./ sigma) .^ 2);
autoPenaltyWeights = 1 + kappa .* (1 - gaussianProximity);

if normalizePerRow
    % Keep mean penalty at 1 to avoid row-wise scale drift near matrix edges.
    denom = mean(autoPenaltyWeights);
    if isfinite(denom) && denom > 0
        autoPenaltyWeights = autoPenaltyWeights ./ denom;
    end
end
end

function penaltyWeights = local_build_penalty_weights(nModelRegCols, nAutoCols, autoPenaltyWeights)
% LOCAL_BUILD_PENALTY_WEIGHTS Assemble per-column ridge penalties for X.
%
% Column layout in X:
%   1) tested predictor (unpenalized),
%   2) other-model nuisance predictors,
%   3) own-model autocorr nuisance predictors.
%
% Sub-steps:
%   1) Start with target predictor weight = 0.
%   2) Assign uniform weight = 1 for other-model nuisance predictors.
%   3) Use strategy-specific autocorr penalty weights when available; fall
%      back to uniform weights when not provided.
penaltyWeights = 0;

if nModelRegCols > 0
    penaltyWeights = [penaltyWeights, ones(1, nModelRegCols)];
end

if nAutoCols > 0
    if isempty(autoPenaltyWeights)
        autoPenaltyWeights = ones(1, nAutoCols);
    end
    if numel(autoPenaltyWeights) ~= nAutoCols
        error('Autocorr penalty weight count (%d) does not match autocorr columns (%d).', ...
            numel(autoPenaltyWeights), nAutoCols);
    end
    penaltyWeights = [penaltyWeights, reshape(autoPenaltyWeights, 1, [])];
end
end

function blockOut = local_reduce_predictor_block(blockIn, params)
% LOCAL_REDUCE_PREDICTOR_BLOCK Optional AdditionalPCA reduction for predictor blocks.
%
% Sub-steps:
%   1) If AdditionalPCA is off, return the block unchanged.
%   2) If block has <= 1 column, PCA is unnecessary; return unchanged.
%   3) Run PCA and keep components with explained variance > 0.1.
%   4) Guarantee at least one component to avoid empty predictor blocks.
blockOut = blockIn;

if isempty(blockIn)
    return;
end
if ~(isfield(params.PCR, 'AdditionalPCA') && isequal(params.PCR.AdditionalPCA, 1))
    return;
end
if size(blockIn, 2) <= 1
    return;
end

[~, score, ~, ~, explained] = pca(blockIn);
imax = sum(explained > 0.1);
imax = max(1, min(imax, size(score, 2)));
blockOut = score(:, 1:imax);
end

function xRegressout = local_normalize_if_requested(xRegressout, params)
% LOCAL_NORMALIZE_IF_REQUESTED Preserve legacy normalization behavior.
%
% The original code normalized only when AdditionalPCA was enabled.
if isempty(xRegressout)
    return;
end
if ~(isfield(params.PCR, 'AdditionalPCA') && isequal(params.PCR.AdditionalPCA, 1))
    return;
end

if strcmp(params.normalization, 'Standardize')
    xRegressout = dRSA_standardizeRDM(xRegressout);
elseif strcmp(params.normalization, 'Rescale')
    xRegressout = dRSA_rescaleRDM(xRegressout);
else
    error('Unknown Normalization method');
end
end

function betaFirst = local_run_original_pcr(Xx, nRDM, params)
% LOCAL_RUN_ORIGINAL_PCR Execute the legacy PCR algorithm and return beta(1,:).
%
% Sub-steps:
%   1) Run PCA according to params.PCR.Method.
%   2) Regress neural RDM on PCA scores.
%   3) Project coefficients back to original predictor space.
%   4) Return coefficient row for the first predictor only.
[nFeature, nPredictors] = size(Xx);

if strcmp(params.PCR.Method, 'FixedComp')
    [PCALoadings, PCAScores] = pca(Xx, 'NumComponents', params.PCR.Methodfactor);

elseif strcmp(params.PCR.Method, 'MinCompPCR')
    minComp = min(nFeature, nPredictors);
    minComp = floor(params.PCR.Methodfactor * minComp);
    [PCALoadings, PCAScores] = pca(Xx, 'NumComponents', minComp);

elseif strcmp(params.PCR.Method, 'ExplainedVar')
    [PCALoadings, PCAScores, ~, ~, explained] = pca(Xx);
    imax = sum(explained > params.PCR.Methodfactor);
    if imax >= nFeature
        error('More Predictors than observations. This will lead to overfitting.');
    end
    imax = max(1, min(imax, size(PCAScores, 2)));
    PCAScores = PCAScores(:, 1:imax);
    PCALoadings = PCALoadings(1, 1:imax);

elseif strcmp(params.PCR.Method, 'CumulativeVar')
    [PCALoadings, PCAScores, ~, ~, explained] = pca(Xx);
    threshold = params.PCR.Methodfactor * sum(explained);
    imax = find(cumsum(explained) >= threshold, 1);
    if imax >= nFeature
        error('More Predictors than observations. This will lead to overfitting.');
    end
    imax = max(1, min(imax, size(PCAScores, 2)));
    PCAScores = PCAScores(:, 1:imax);
    PCALoadings = PCALoadings(1, 1:imax);

else
    error('Unknown params.PCR.Method: %s', params.PCR.Method);
end

betaPCR = PCAScores \ nRDM;
temporarydRSA = PCALoadings * betaPCR;

% The first row corresponds to the target predictor column in Xx.
betaFirst = temporarydRSA(1, :);
end

function betaFirst = local_ridge_first_predictor(X, Y, lambdaFactor, standardizePredictors, penaltyWeights)
% LOCAL_RIDGE_FIRST_PREDICTOR Solve multiresponse ridge and keep first beta row.
%
% Inputs:
%   X: nObs x nPredictors design matrix (first column = tested effect).
%   Y: nObs x nResponses response matrix (all neural-time points).
%   lambdaFactor: scalar, lambda = lambdaFactor * nObs.
%   standardizePredictors: true/false.
%   penaltyWeights: 1 x nPredictors vector with per-column ridge weights.
%       Convention: target predictor weight should be 0.
%
% Sub-steps:
%   1) Drop non-finite nuisance columns (but require finite target column).
%   2) Align penalty weights to surviving predictor columns.
%   3) Optionally z-score predictors for numerical stability.
%   4) Mean-center responses.
%   5) Solve weighted ridge in standardized predictor space.
%   6) Convert beta back to original predictor scale.
X = double(X);
Y = double(Y);

if isempty(X)
    betaFirst = zeros(1, size(Y, 2));
    return;
end

validCols = all(isfinite(X), 1);
if ~validCols(1)
    betaFirst = zeros(1, size(Y, 2));
    return;
end

if nargin < 5 || isempty(penaltyWeights)
    penaltyWeights = [0, ones(1, size(X, 2) - 1)];
end
penaltyWeights = reshape(penaltyWeights, 1, []);
if numel(penaltyWeights) ~= size(X, 2)
    error('penaltyWeights length (%d) must match predictor count (%d).', ...
        numel(penaltyWeights), size(X, 2));
end

X = X(:, validCols);
penaltyWeights = penaltyWeights(validCols);
firstColIdx = 1;

if standardizePredictors
    muX = mean(X, 1);
    stdX = std(X, 0, 1);
    stdX(stdX < eps) = 1;
    Xproc = bsxfun(@rdivide, bsxfun(@minus, X, muX), stdX);
else
    stdX = ones(1, size(X, 2));
    Xproc = X;
end

Yproc = bsxfun(@minus, Y, mean(Y, 1));

nObs = size(Xproc, 1);
lambda = lambdaFactor * nObs;

% Weighted ridge penalty:
%   lambda * beta' * diag(penaltyWeights) * beta
% Applied in standardized predictor space when standardization is enabled.
penalty = diag(penaltyWeights);
penalty(firstColIdx, firstColIdx) = 0;

lhs = (Xproc' * Xproc) + lambda * penalty;
rhs = Xproc' * Yproc;
betaProc = lhs \ rhs;

betaFirst = betaProc(firstColIdx, :) ./ stdX(firstColIdx);
end
