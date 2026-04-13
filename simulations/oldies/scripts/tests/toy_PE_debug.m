% toy_PE_debug.m
%
% Purpose:
%   Build the smallest reproducible dRSA debug run for the PE position model
%   autocorrelation used in simulations/scripts/PE_simulation_diff.m.
%
% Canonical PE definition (same as the main PE script):
%   PE_signal = observed_signal - predicted_signal
%
% Scope of this toy script:
%   - Use deviant observed and predicted paths for one participant.
%   - Enforce the PE post-deviance window (same rule as PE_simulation_diff).
%   - Focus only on PE position dot1.
%   - Run corr-style dRSA with PE model against itself (autocorrelation).
%   - Plot only the resulting time x time matrix.
%
% Data flow summary:
%   1) Load (or build) observed and predicted deviant center-relative paths.
%   2) Cut both streams to the post-deviance window.
%   3) Compute PE position dot1 = observed dot1 - predicted dot1.
%   4) Concatenate trials, build trial-locked subsamples, run dRSA autocorr.
%   5) Plot the PE autocorrelation matrix.
%
% Usage example (from this folder):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/tests');
%   toy_PE_debug;
%
% Usage example (from anywhere):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/tests');
%   toy_PE_debug;
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat
%   - experiment/input_files/MovDot_SubXX_predicted.mat
%   - simulations/input/MovDot_SubXX_deviant.mat
%   - simulations/input/MovDot_SubXX_predicted_deviant.mat
%
% Optional output:
%   - simulations/scripts/tests/derivatives/subXX/toy_PE_debug_subXX_*.png
%   - simulations/scripts/tests/derivatives/subXX/toy_PE_debug_positions_subXX_*.png
%
% Assumptions:
%   - Center-relative path fields exist in both deviant files.
%   - Path arrays are trials x 2 x time.
%   - Observed and predicted deviant arrays are shape-aligned.

clear all
close all

%% Resolve repository paths and add dRSA dependencies
% Data flow: script path -> repo root -> helper folders on MATLAB path.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
addpath(fullfile(repoRoot, 'simulations'));
addpath(fullfile(repoRoot, 'simulations', 'debug'));
addpath(fullfile(repoRoot, 'simulations', 'functions'));

%% User configuration
% Data flow: editable config -> input selection, PE windowing, and plotting.
participantNumber = 76;
deviantOnset = 0.5; % fraction of trial duration where deviance starts.
enforcePostDevCut = true; % keep aligned with PE_simulation_diff PE branch.
plotSampleRateHz = 120;
suppressDispText = 1;
figureVisibility = 'on'; % set to 'off' for headless batch usage.
useCommonColorLimits = true;
matrixColorLimits = [0 1];
positionWindowSamples = 20; % number of post-deviant samples shown for start/end windows.
recenterWindowsAtDeviant = true; % true to shift each trial so deviant-point sample is at [0, 0].
titleLayoutMode = 'compact'; % 'compact' for short titles, 'legacy' for verbose titles.
saveFigure = true;
outputDpi = 300;
testsDerivativesDir = fullfile(scriptDir, 'derivatives');
local_ensure_dir(testsDerivativesDir);
if ~(isscalar(recenterWindowsAtDeviant) && ...
        (islogical(recenterWindowsAtDeviant) || isnumeric(recenterWindowsAtDeviant)))
    error('recenterWindowsAtDeviant must be a scalar logical/numeric value.');
end
if isstring(titleLayoutMode) && isscalar(titleLayoutMode)
    titleLayoutMode = char(titleLayoutMode);
end
if ~ischar(titleLayoutMode)
    error('titleLayoutMode must be ''compact'' or ''legacy''.');
end
titleLayoutMode = lower(strtrim(titleLayoutMode));
validTitleLayoutModes = {'compact', 'legacy'};
if ~ismember(titleLayoutMode, validTitleLayoutModes)
    error('titleLayoutMode must be one of: %s', strjoin(validTitleLayoutModes, ', '));
end

%% Resolve input paths and build missing deviant simulation files
% Data flow: participant id -> expected files -> optional build step.
simulationInputDir = fullfile(repoRoot, 'simulations', 'input');
observedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
predictedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d_predicted.mat', participantNumber));
deviantSimulationFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_deviant.mat', participantNumber));
predictedSimulationFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_predicted_deviant.mat', participantNumber));

if ~isfile(observedInputFile)
    error('Observed input file not found: %s', observedInputFile);
end
if ~isfile(predictedInputFile)
    error('Predicted input file not found: %s', predictedInputFile);
end

if ~isfile(deviantSimulationFile)
    build_movdot_simulation_inputs(observedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'suppressDispText', suppressDispText);
end
if ~isfile(predictedSimulationFile)
    build_movdot_simulation_inputs(predictedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'AllowMissingConditions', true, ...
        'suppressDispText', suppressDispText);
end

%% Load observed/predicted deviant dot1 paths and validate alignment
% Data flow: deviant files -> center-relative arrays -> size checks.
observedData = load(deviantSimulationFile);
predictedData = load(predictedSimulationFile);
requiredField = 'dot1GreenPathsCenterRelative';
if ~isfield(observedData, requiredField)
    error('Observed deviant file lacks %s: %s', requiredField, deviantSimulationFile);
end
if ~isfield(predictedData, requiredField)
    error('Predicted deviant file lacks %s: %s', requiredField, predictedSimulationFile);
end

observedDot1 = observedData.dot1GreenPathsCenterRelative;
predictedDot1 = predictedData.dot1GreenPathsCenterRelative;
if ~isequal(size(observedDot1), size(predictedDot1))
    error(['Observed/predicted dot1 shapes must match. ', ...
        'Observed=%s Predicted=%s'], ...
        mat2str(size(observedDot1)), mat2str(size(predictedDot1)));
end

%% Enforce PE post-deviance window (same PE branch rule as main script)
% Data flow: full-trial observed+predicted dot1 -> PE window-aligned paths.
cutFrame = [];
if enforcePostDevCut
    [observedDot1, cutFrame] = local_cut_postdeviant(observedDot1, deviantOnset);
    [predictedDot1, ~] = local_cut_postdeviant(predictedDot1, deviantOnset);
end

%% Build PE position dot1 stream
% Data flow: observed - predicted -> PE dot1 -> concatenated dRSA stream.
peDot1 = observedDot1 - predictedDot1;
plotConcat = 0; % keep toy run minimal and avoid concatenation plots.
dataPEPositionDot1 = dRSA_concatenate(peDot1, [], plotConcat, ...
    'suppressDispText', suppressDispText);

%% Plot start/end post-deviance position windows (observed, predicted, PE)
% Data flow: post-deviance streams -> start/end windows -> trial-colored overlays.
figPositionWindows = local_plot_position_windows( ...
    observedDot1, predictedDot1, peDot1, positionWindowSamples, ...
    figureVisibility, participantNumber, cutFrame, enforcePostDevCut, ...
    recenterWindowsAtDeviant, titleLayoutMode);

%% Build trial-locked subsamples for dRSA
% Data flow: PE trial length + concatenated stream length -> trigger windows.
trialLen = size(peDot1, 3);
totalTime = size(dataPEPositionDot1, 2);
if mod(totalTime, trialLen) ~= 0
    error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
end

trialStarts = 1:trialLen:totalTime;
maskSubsampling = true(1, totalTime);
maskTrigger = false(1, totalTime);
maskTrigger(trialStarts) = true;

opt.PreTrigger = 0;
opt.PostTrigger = trialLen - 1;
opt.spacing = 0;
opt.nSubSamples = numel(trialStarts);
opt.nIter = 1;
opt.checkRepetition = 0;

subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt, ...
    'suppressDispText', suppressDispText);

%% Run corr-style dRSA autocorrelation (PE position dot1 with itself)
% Data flow: PE dot1 stream + subsamples -> dRSA autocorrelation matrix.
params = struct();
params.nIter = 1;
params.fs = plotSampleRateHz;
params.AverageTime = floor((trialLen - 1) / 2) / params.fs;
params.modelToTest = 1;
params.Var = 0.1;
params.modelNames = {'PE position dot1'};
% dRSA_coreFunction autocorr mode expects a single distance string here.
params.modelDistMeasure = 'euclidean';
params.neuralDistMeasure = 'euclidean';
params.dRSAtype = 'corr';

modelPE = {dataPEPositionDot1};
dRSA_PE_autocorr = dRSA_coreFunction(dataPEPositionDot1, modelPE, params, ...
    'CurrSubsamples', subsamples(:, :, 1), ...
    'Autocorrborder', [], ...
    'suppressDispText', suppressDispText);
dRSA_PE_autocorr = reshape(dRSA_PE_autocorr, ...
    size(dRSA_PE_autocorr, 1), size(dRSA_PE_autocorr, 2));

%% Plot PE autocorrelation matrix
% Data flow: autocorrelation matrix -> time-axis heatmap for inspection.
nTime = size(dRSA_PE_autocorr, 1);
timeSeconds = (0:nTime - 1) / params.fs;

fig = figure('Name', 'Toy PE debug: position dot1 autocorr', ...
    'NumberTitle', 'off', 'Visible', figureVisibility);
imagesc(timeSeconds, timeSeconds, dRSA_PE_autocorr);
set(gca, 'YDir', 'normal');
axis image;
if useCommonColorLimits
    caxis(matrixColorLimits);
end
colorbar;
hold on;
plot(timeSeconds, timeSeconds, 'w-', 'LineWidth', 1);
hold off;
xlabel('Time in PE neural stream (s)');
ylabel('Time in PE model stream (s)');
if enforcePostDevCut
    title(sprintf('PE position dot1 autocorr (sub%02d, postDev from frame %d)', ...
        participantNumber, cutFrame));
else
    title(sprintf('PE position dot1 autocorr (sub%02d, full trial)', participantNumber));
end

%% Optional figure export
% Data flow: on-screen matrix figure -> deterministic output path -> PNG.
if saveFigure
    subjectLabel = sprintf('sub%02d', participantNumber);
    outputDir = fullfile(testsDerivativesDir, subjectLabel);
    local_ensure_dir(outputDir);
    if enforcePostDevCut
        runTag = sprintf('postDevFrom_%d', cutFrame);
    else
        runTag = 'fullTrial';
    end
    outFile = fullfile(outputDir, ...
        sprintf('toy_PE_debug_%s_%s.png', subjectLabel, runTag));
    outFilePositions = fullfile(outputDir, ...
        sprintf('toy_PE_debug_positions_%s_%s.png', subjectLabel, runTag));
    print(fig, outFile, '-dpng', sprintf('-r%d', outputDpi));
    print(figPositionWindows, outFilePositions, '-dpng', sprintf('-r%d', outputDpi));
    fprintf('Saved toy PE debug matrix to:\n  %s\n', outFile);
    fprintf('Saved toy PE position windows to:\n  %s\n', outFilePositions);
end

%% Local helpers
function fig = local_plot_position_windows( ...
        observedDot1, predictedDot1, peDot1, positionWindowSamples, ...
        figVisibility, participantNumber, cutFrame, enforcePostDevCut, ...
        recenterWindowsAtDeviant, titleLayoutMode)
% LOCAL_PLOT_POSITION_WINDOWS Plot post-deviance start/end windows by trial.
%
% Purpose:
%   Show spatial structure of dot1 paths across trials in two time windows:
%   early post-deviance and late post-deviance. Three datasets are shown:
%   observed, predicted, and PE.
%
% Inputs:
%   observedDot1         : trials x 2 x time observed post-deviance positions.
%   predictedDot1        : trials x 2 x time predicted post-deviance positions.
%   peDot1               : trials x 2 x time PE positions.
%   positionWindowSamples: number of samples used for each window.
%   figVisibility        : 'on' or 'off'.
%   participantNumber    : numeric subject id for title.
%   cutFrame             : 1-based cut frame when post-deviance is enforced.
%   enforcePostDevCut    : true when post-deviance cut is active.
%   recenterWindowsAtDeviant : true to shift each trial so sample 1 is at [0, 0].
%   titleLayoutMode      : 'compact' or 'legacy' title strategy.
%
% Output:
%   fig : figure handle with a 2x3 panel layout.
%
% Data flow:
%   post-deviance paths -> optional per-trial recentering -> start/end windows
%   -> per-trial color overlays.

datasetsRaw = {observedDot1, predictedDot1, peDot1};
if recenterWindowsAtDeviant
    datasets = cellfun(@local_recenter_paths_at_deviant, datasetsRaw, 'UniformOutput', false);
    local_assert_recentered_origins(datasets, 1e-10);
else
    datasets = datasetsRaw;
end

nTime = size(datasets{1}, 3);
windowLen = min(positionWindowSamples, nTime);
startIdx = 1:windowLen;
endIdx = (nTime - windowLen + 1):nTime;
nTrials = size(datasets{1}, 1);
trialColors = local_trial_colors(nTrials);

if strcmp(titleLayoutMode, 'compact')
    colTitles = {'Observed', 'Predicted', 'PE'};
    fig = figure( ...
        'Name', 'Toy PE debug: position windows (dot1)', ...
        'NumberTitle', 'off', ...
        'Visible', figVisibility, ...
        'Position', [80, 60, 1500, 980]);
    tl = tiledlayout(2, 3, 'Padding', 'normal', 'TileSpacing', 'compact');
else
    colTitles = {'Observed dot1', 'Predicted dot1', 'PE dot1 (observed - predicted)'};
    fig = figure( ...
        'Name', 'Toy PE debug: position windows (dot1)', ...
        'NumberTitle', 'off', ...
        'Visible', figVisibility);
    tl = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
end

startLimits = local_window_limits(datasets, startIdx);
endLimits = local_window_limits(datasets, endIdx);

rowAxesTop = gobjects(1, 3);
rowAxesBottom = gobjects(1, 3);
for iCol = 1:3
    ax = nexttile(iCol);
    rowAxesTop(iCol) = ax;
    local_plot_trial_window(ax, datasets{iCol}, startIdx, trialColors);
    if strcmp(titleLayoutMode, 'compact')
        title(ax, colTitles{iCol}, 'FontWeight', 'bold', 'FontSize', 10);
        if iCol == 1
            ylabel(ax, sprintf('Start window\nY (visual degrees)'));
        else
            ylabel(ax, '');
        end
    else
        title(ax, sprintf('%s | Start window (%d samples after dev)', colTitles{iCol}, windowLen));
        ylabel(ax, 'Y (visual degrees)');
    end
    xlabel(ax, 'X (visual degrees)');
    xlim(ax, startLimits(1:2));
    ylim(ax, startLimits(3:4));
end

for iCol = 1:3
    ax = nexttile(3 + iCol);
    rowAxesBottom(iCol) = ax;
    local_plot_trial_window(ax, datasets{iCol}, endIdx, trialColors);
    if strcmp(titleLayoutMode, 'compact')
        title(ax, '', 'FontSize', 10);
        if iCol == 1
            ylabel(ax, sprintf('End window\nY (visual degrees)'));
        else
            ylabel(ax, '');
        end
    else
        title(ax, sprintf('%s | End window (%d samples before trial end)', colTitles{iCol}, windowLen));
        ylabel(ax, 'Y (visual degrees)');
    end
    xlabel(ax, 'X (visual degrees)');
    xlim(ax, endLimits(1:2));
    ylim(ax, endLimits(3:4));
end

if strcmp(titleLayoutMode, 'compact')
    if enforcePostDevCut
        titleText = sprintf('Sub%02d dot1 windows | post-dev frame %d | same color = same trial', ...
            participantNumber, cutFrame);
    else
        titleText = sprintf('Sub%02d dot1 windows | same color = same trial', participantNumber);
    end
    title(tl, titleText, 'FontWeight', 'bold', 'FontSize', 11);
else
    if enforcePostDevCut
        title(tl, sprintf(['Sub%02d trial-colored windows (postDev from frame %d). ', ...
            'Same color = same trial across all panels'], participantNumber, cutFrame));
    else
        title(tl, sprintf(['Sub%02d trial-colored windows (full trial). ', ...
            'Same color = same trial across all panels'], participantNumber));
    end
end
end

function pathsCentered = local_recenter_paths_at_deviant(paths)
% LOCAL_RECENTER_PATHS_AT_DEVIANT Shift each trial to origin at sample 1.
%
% Purpose:
%   Recenter each path in a post-deviance stream so the first sample
%   (deviant point after cutting) is exactly at [0, 0].
%
% Input:
%   paths : trials x 2 x time tensor.
%
% Output:
%   pathsCentered : paths shifted by per-trial sample-1 anchor.
%
% Data flow:
%   trial tensor -> sample-1 anchor per trial -> broadcast subtraction.
if isempty(paths)
    pathsCentered = paths;
    return;
end
anchor = paths(:, :, 1);
pathsCentered = paths - repmat(anchor, 1, 1, size(paths, 3));
end

function local_assert_recentered_origins(datasets, tolerance)
% LOCAL_ASSERT_RECENTERED_ORIGINS Verify recentered sample-1 origins.
%
% Inputs:
%   datasets  : cell array of recentered path tensors.
%   tolerance : absolute tolerance for near-zero origin check.
%
% Data flow:
%   sample-1 coordinates -> absolute max -> strict assertion.
for iData = 1:numel(datasets)
    origins = datasets{iData}(:, :, 1);
    maxAbsOrigin = max(abs(origins(:)));
    assert(maxAbsOrigin <= tolerance, ...
        'Recentering failed for dataset %d (max abs origin %.3g > %.3g).', ...
        iData, maxAbsOrigin, tolerance);
end
end

function local_plot_trial_window(ax, paths, windowIdx, trialColors)
% LOCAL_PLOT_TRIAL_WINDOW Plot one dataset window with one color per trial.
%
% Inputs:
%   ax         : target axes.
%   paths      : trials x 2 x time tensor.
%   windowIdx  : sample indices in the third dimension to draw.
%   trialColors: nTrials x 3 RGB color matrix.
%
% Data flow:
%   selected window samples per trial -> colored line overlays with endpoint marks.

hold(ax, 'on');
nTrials = size(paths, 1);
for iTrial = 1:nTrials
    xVals = squeeze(paths(iTrial, 1, windowIdx));
    yVals = squeeze(paths(iTrial, 2, windowIdx));
    plot(ax, xVals, yVals, '-', 'Color', trialColors(iTrial, :), 'LineWidth', 1.0);
    plot(ax, xVals(1), yVals(1), 'o', ...
        'Color', trialColors(iTrial, :), 'MarkerFaceColor', trialColors(iTrial, :), ...
        'MarkerSize', 3);
    plot(ax, xVals(end), yVals(end), 's', ...
        'Color', trialColors(iTrial, :), 'MarkerFaceColor', trialColors(iTrial, :), ...
        'MarkerSize', 3);
end
axis(ax, 'equal');
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');
end

function limits = local_window_limits(datasets, windowIdx)
% LOCAL_WINDOW_LIMITS Compute common axis limits across datasets for a window.
%
% Inputs:
%   datasets : cell array of path tensors (trials x 2 x time).
%   windowIdx: indices of the selected time window.
%
% Output:
%   limits : [xmin xmax ymin ymax] with padding.
%
% Data flow:
%   all dataset points in the selected window -> min/max -> padded limits.

xAll = [];
yAll = [];
for iData = 1:numel(datasets)
    data = datasets{iData};
    xAll = [xAll; reshape(data(:, 1, windowIdx), [], 1)];
    yAll = [yAll; reshape(data(:, 2, windowIdx), [], 1)];
end

xMin = min(xAll);
xMax = max(xAll);
yMin = min(yAll);
yMax = max(yAll);

xSpan = max(xMax - xMin, eps);
ySpan = max(yMax - yMin, eps);
padX = 0.05 * xSpan;
padY = 0.05 * ySpan;
limits = [xMin - padX, xMax + padX, yMin - padY, yMax + padY];
end

function colors = local_trial_colors(nTrials)
% LOCAL_TRIAL_COLORS Build visually distinct per-trial colors.
%
% Input:
%   nTrials : number of trial colors needed.
%
% Output:
%   colors  : nTrials x 3 RGB matrix.
if nTrials <= 0
    colors = zeros(0, 3);
    return;
end
colors = hsv(nTrials);
end

function [pathsCut, cutFrame] = local_cut_postdeviant(paths, deviantOnset)
% LOCAL_CUT_POSTDEVIANT Keep the post-deviant portion of each trial.
%
% Usage example:
%   [pathsCut, cutFrame] = local_cut_postdeviant(paths, 0.5);
%
% Inputs:
%   paths        : trials x features x time center-relative positions.
%   deviantOnset : fraction of trial duration (0-1) marking deviant onset.
%
% Outputs:
%   pathsCut     : paths sliced from deviant onset to the end.
%   cutFrame     : 1-based frame index where the cut begins.
%
% Data flow:
%   paths -> onset fraction -> frame index -> post-deviant-only paths.
if isempty(paths)
    pathsCut = paths;
    cutFrame = [];
    return;
end
if ~isscalar(deviantOnset) || deviantOnset < 0 || deviantOnset > 1
    error('deviantOnset must be a scalar between 0 and 1.');
end
trialLen = size(paths, 3);
if trialLen < 1
    pathsCut = paths;
    cutFrame = [];
    return;
end
cutFrame = round(deviantOnset * trialLen);
cutFrame = max(1, min(trialLen, cutFrame));
pathsCut = paths(:, :, cutFrame:end);
end

function local_ensure_dir(dirPath)
% LOCAL_ENSURE_DIR Create output directory if it is missing.
%
% Input:
%   dirPath : directory path to create.
if isempty(dirPath)
    return;
end
if ~exist(dirPath, 'dir')
    mkdir(dirPath);
end
end
