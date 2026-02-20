% toy_PE_diff.m
%
% Purpose:
%   Build a compact toy visualization of prediction-error (PE) geometry using
%   a small subset of deviant trials and their predicted counterparts.
%
%   For each selected trial and dot:
%     - deviant path is plotted
%     - predicted path is plotted
%     - PE path is computed as PE = deviant - predicted
%     - PE is then centered on the deviant path and overlaid as:
%         deviant + PE
%     - PE is also plotted from the deviant start point only:
%         deviant_origin + (PE - PE_at_start)
%     - PE vectors (arrows) are drawn from each deviant sample using the PE
%       displacement at that sample.
%
% Why this script exists:
%   The main PE simulation plot uses scatter points and can look visually
%   "straight/dotted" depending on path structure and axis limits. This toy
%   script isolates a few trials and explicitly overlays deviant/predicted/PE
%   components on the same axes with dynamic limits for easier inspection.
%
% Example usage (from this folder):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/hypothesis_testing');
%   toy_PE_diff;
%
% Example usage (from anywhere):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/hypothesis_testing');
%   toy_PE_diff;
%
% Inputs:
%   - simulations/input/MovDot_SubXX_deviant.mat
%   - simulations/input/MovDot_SubXX_predicted_deviant.mat
%   - If missing, files are generated from:
%       experiment/input_files/MovDot_SubXX.mat
%       experiment/input_files/MovDot_SubXX_predicted.mat
%
% Outputs:
%   - Figure saved to:
%       simulations/output/subXX/paths/toy_PE_diff_subXX_trials_<...>.png
%
% Key assumptions:
%   - Center-relative path fields exist:
%       dot1GreenPathsCenterRelative
%       dot2YellowPathsCenterRelative
%   - Path arrays are trials x 2 x time in visual degrees.
%   - Predicted deviant and observed deviant arrays are shape-aligned.

clear all
close all

%% Resolve repository paths and add dependencies
% Data flow: script path -> repo root -> helper functions on MATLAB path.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot, 'simulations'));
addpath(fullfile(repoRoot, 'simulations', 'debug'));
addpath(fullfile(repoRoot, 'simulations', 'functions'));

%% User configuration
% Data flow: user-editable config -> file loading/cropping/plot export behavior.
participantNumber = 87;
trialIndices = [1, 2]; % "a couple of paths only" by default.
enforcePostDev = true; % true keeps only post-deviance samples for toy visualization.
deviantOnset = 0.5; % fraction of trial duration where deviance starts.
arrowStride = 16; % draw one PE vector every N time samples.
lineWidth = 1.6;
figureVisibility = 'on'; % 'on' for interactive inspection, 'off' for headless export.
saveFigure = true;
outputDpi = 300;

%% Resolve input files and build missing simulation inputs if needed
% Data flow: expected input paths -> optional build step -> concrete files for loading.
simulationInputDir = fullfile(repoRoot, 'simulations', 'input');
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
observedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
predictedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d_predicted.mat', participantNumber));

deviantSimulationFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_deviant.mat', participantNumber));
predictedSimulationFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_predicted_deviant.mat', participantNumber));

if ~isfile(deviantSimulationFile)
    build_movdot_simulation_inputs(observedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript);
end
if ~isfile(predictedSimulationFile)
    build_movdot_simulation_inputs(predictedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'AllowMissingConditions', true);
end

%% Load observed/predicted deviant center-relative paths
% Data flow: deviant files -> dot path tensors -> toy subset extraction.
observedData = load(deviantSimulationFile);
predictedData = load(predictedSimulationFile);
requiredFields = {'dot1GreenPathsCenterRelative', 'dot2YellowPathsCenterRelative'};
if ~all(isfield(observedData, requiredFields))
    error('Observed deviant file missing required center-relative fields: %s', deviantSimulationFile);
end
if ~all(isfield(predictedData, requiredFields))
    error('Predicted deviant file missing required center-relative fields: %s', predictedSimulationFile);
end

obsDot1 = observedData.dot1GreenPathsCenterRelative;
obsDot2 = observedData.dot2YellowPathsCenterRelative;
predDot1 = predictedData.dot1GreenPathsCenterRelative;
predDot2 = predictedData.dot2YellowPathsCenterRelative;

if ~isequal(size(obsDot1), size(predDot1)) || ~isequal(size(obsDot2), size(predDot2))
    error(['Observed and predicted deviant paths must have identical shapes. ', ...
        'Observed dot1=%s predicted dot1=%s'], mat2str(size(obsDot1)), mat2str(size(predDot1)));
end

%% Optional post-deviance crop (kept explicit for reproducible toy interpretation)
% Data flow: full-trial paths -> post-deviance paths (optional).
cutFrame = [];
if enforcePostDev
    [obsDot1, cutFrame] = local_cut_postdeviant(obsDot1, deviantOnset);
    [obsDot2, ~] = local_cut_postdeviant(obsDot2, deviantOnset);
    [predDot1, ~] = local_cut_postdeviant(predDot1, deviantOnset);
    [predDot2, ~] = local_cut_postdeviant(predDot2, deviantOnset);
end

%% Keep only the requested small trial subset
% Data flow: full condition tensors -> selected trials -> toy plotting tensors.
trialIndices = local_validate_trial_indices(trialIndices, size(obsDot1, 1));
obsDot1 = obsDot1(trialIndices, :, :);
obsDot2 = obsDot2(trialIndices, :, :);
predDot1 = predDot1(trialIndices, :, :);
predDot2 = predDot2(trialIndices, :, :);

%% Compute PE paths and center PE on deviant points
% Data flow: observed/predicted -> PE displacement -> PE centered on deviant.
peDot1 = obsDot1 - predDot1;
peDot2 = obsDot2 - predDot2;
peCenteredDot1 = obsDot1 + peDot1;
peCenteredDot2 = obsDot2 + peDot2;
peFromOriginDot1 = local_anchor_pe_to_deviant_origin(obsDot1, peDot1);
peFromOriginDot2 = local_anchor_pe_to_deviant_origin(obsDot2, peDot2);

%% Plot overlay for dot1 and dot2 on a single figure
% Data flow: selected paths + PE displacements -> 1x2 overlay figure.
subjectLabel = sprintf('sub%02d', participantNumber);
fig = figure('Name', sprintf('Toy PE diff (%s)', subjectLabel), ...
    'NumberTitle', 'off', 'Visible', figureVisibility);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile(1);
local_plot_toy_overlay(ax1, obsDot1, predDot1, peCenteredDot1, peFromOriginDot1, peDot1, ...
    'dot1', trialIndices, arrowStride, lineWidth);

ax2 = nexttile(2);
local_plot_toy_overlay(ax2, obsDot2, predDot2, peCenteredDot2, peFromOriginDot2, peDot2, ...
    'dot2', trialIndices, arrowStride, lineWidth);

if enforcePostDev
    sgtitle(fig, sprintf('Toy PE diff (%s) | trials %s | postDev from frame %d', ...
        subjectLabel, local_trial_label(trialIndices), cutFrame));
else
    sgtitle(fig, sprintf('Toy PE diff (%s) | trials %s | full trial', ...
        subjectLabel, local_trial_label(trialIndices)));
end

%% Save figure to subject paths output folder
% Data flow: in-memory figure -> deterministic output filename -> PNG.
if saveFigure
    pathsOutputDir = fullfile(repoRoot, 'simulations', 'output', subjectLabel, 'paths');
    local_ensure_dir(pathsOutputDir);
    if enforcePostDev
        runTag = sprintf('postDevFrom_%d', cutFrame);
    else
        runTag = 'fullTrial';
    end
    outBase = sprintf('toy_PE_diff_%s_trials_%s_%s', ...
        subjectLabel, local_trial_label(trialIndices), runTag);
    outFile = fullfile(pathsOutputDir, [outBase '.png']);
    print(fig, outFile, '-dpng', sprintf('-r%d', outputDpi));
    fprintf('Saved toy PE figure to:\n  %s\n', outFile);
end

%% Local helpers
function local_plot_toy_overlay(ax, observedPaths, predictedPaths, centeredPEPaths, ...
        peFromOriginPaths, peDisplacements, ...
        dotLabel, trialIndices, arrowStride, lineWidth)
% LOCAL_PLOT_TOY_OVERLAY Overlay deviant, predicted, and centered-PE paths.
%
% Inputs:
%   ax              : target axes.
%   observedPaths   : selected deviant paths (trials x 2 x time).
%   predictedPaths  : selected predicted paths (trials x 2 x time).
%   centeredPEPaths : observed + PE (trials x 2 x time).
%   peFromOriginPaths : start-anchored PE path using deviant trial origin.
%   peDisplacements : PE vectors (trials x 2 x time).
%   dotLabel        : dot identifier used in title.
%   trialIndices    : original trial indices shown.
%   arrowStride     : spacing for PE vector quiver arrows.
%   lineWidth       : base line width for all curve layers.
%
% Data flow:
%   selected per-trial arrays -> layered line plots + PE quiver vectors.
hold(ax, 'on');

obsColor = [0.00, 0.45, 0.74];
predColor = [0.85, 0.33, 0.10];
peColor = [0.47, 0.67, 0.19];
peOriginColor = [0.20, 0.20, 0.20];
peVectorColor = [0.15, 0.60, 0.15];

hObserved = gobjects(1);
hPredicted = gobjects(1);
hPECentered = gobjects(1);
hPEFromOrigin = gobjects(1);
hPEVector = gobjects(1);

nTrials = size(observedPaths, 1);
for iTrial = 1:nTrials
    xObs = squeeze(observedPaths(iTrial, 1, :));
    yObs = squeeze(observedPaths(iTrial, 2, :));
    xPred = squeeze(predictedPaths(iTrial, 1, :));
    yPred = squeeze(predictedPaths(iTrial, 2, :));
    xPeCentered = squeeze(centeredPEPaths(iTrial, 1, :));
    yPeCentered = squeeze(centeredPEPaths(iTrial, 2, :));
    xPeFromOrigin = squeeze(peFromOriginPaths(iTrial, 1, :));
    yPeFromOrigin = squeeze(peFromOriginPaths(iTrial, 2, :));

    hObsCurr = plot(ax, xObs, yObs, '-', 'Color', obsColor, 'LineWidth', lineWidth);
    hPredCurr = plot(ax, xPred, yPred, '--', 'Color', predColor, 'LineWidth', lineWidth);
    hPeCurr = plot(ax, xPeCentered, yPeCentered, '-.', 'Color', peColor, 'LineWidth', lineWidth);
    hPeOriginCurr = plot(ax, xPeFromOrigin, yPeFromOrigin, ':', ...
        'Color', peOriginColor, 'LineWidth', lineWidth);

    arrowIdx = 1:max(1, arrowStride):numel(xObs);
    q = quiver(ax, xObs(arrowIdx), yObs(arrowIdx), ...
        squeeze(peDisplacements(iTrial, 1, arrowIdx)), ...
        squeeze(peDisplacements(iTrial, 2, arrowIdx)), ...
        0, 'Color', peVectorColor, 'LineWidth', 0.9, 'MaxHeadSize', 0.8);

    % Capture legend handles once to avoid repeated legend entries.
    if iTrial == 1
        hObserved = hObsCurr;
        hPredicted = hPredCurr;
        hPECentered = hPeCurr;
        hPEFromOrigin = hPeOriginCurr;
        hPEVector = q;
    end

    % Mark each selected trial start point for orientation.
    scatter(ax, xObs(1), yObs(1), 24, obsColor, 'filled');
end

local_apply_dynamic_limits(ax, observedPaths, predictedPaths, centeredPEPaths, peFromOriginPaths);
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'X (visual degrees)');
ylabel(ax, 'Y (visual degrees)');
title(ax, sprintf('%s | trials %s', dotLabel, local_trial_label(trialIndices)));
legend(ax, [hObserved, hPredicted, hPECentered, hPEFromOrigin, hPEVector], ...
    {'Deviant', 'Predicted', 'Deviant + PE (pointwise)', ...
    'PE + deviant origin (start-anchored)', 'PE vector at deviant'}, ...
    'Location', 'best');
hold(ax, 'off');
end

function local_apply_dynamic_limits(ax, observedPaths, predictedPaths, centeredPEPaths, peFromOriginPaths)
% LOCAL_APPLY_DYNAMIC_LIMITS Use data-driven limits with margin.
%
% Data flow:
%   all plotted coordinates -> global min/max -> padded axis limits.
xAll = [
    reshape(observedPaths(:, 1, :), [], 1);
    reshape(predictedPaths(:, 1, :), [], 1);
    reshape(centeredPEPaths(:, 1, :), [], 1);
    reshape(peFromOriginPaths(:, 1, :), [], 1)];
yAll = [
    reshape(observedPaths(:, 2, :), [], 1);
    reshape(predictedPaths(:, 2, :), [], 1);
    reshape(centeredPEPaths(:, 2, :), [], 1);
    reshape(peFromOriginPaths(:, 2, :), [], 1)];

xMin = min(xAll);
xMax = max(xAll);
yMin = min(yAll);
yMax = max(yAll);

xRange = xMax - xMin;
yRange = yMax - yMin;
xPad = max(0.5, 0.08 * xRange); % minimum pad avoids near-flat clipping.
yPad = max(0.5, 0.08 * yRange);

xlim(ax, [xMin - xPad, xMax + xPad]);
ylim(ax, [yMin - yPad, yMax + yPad]);
end

function peFromOriginPaths = local_anchor_pe_to_deviant_origin(observedPaths, peDisplacements)
% LOCAL_ANCHOR_PE_TO_DEVIANT_ORIGIN Anchor PE trajectories to deviant start point.
%
% Inputs:
%   observedPaths   : trials x 2 x time deviant trajectories.
%   peDisplacements : trials x 2 x time PE = deviant - predicted.
%
% Outputs:
%   peFromOriginPaths : trials x 2 x time where each PE trial starts exactly
%       at the deviant first sample (start-anchored PE shape).
%
% Data flow:
%   PE displacements -> subtract first PE sample -> add deviant origin.
if ~isequal(size(observedPaths), size(peDisplacements))
    error('local_anchor_pe_to_deviant_origin requires shape-aligned inputs.');
end

peFromOriginPaths = zeros(size(peDisplacements));
nTrials = size(peDisplacements, 1);
for iTrial = 1:nTrials
    devOrigin = observedPaths(iTrial, :, 1); % 1x2 origin for this trial.
    peStart = peDisplacements(iTrial, :, 1);
    peFromOriginPaths(iTrial, :, :) = devOrigin + (peDisplacements(iTrial, :, :) - peStart);
end
end

function trialIndices = local_validate_trial_indices(trialIndicesIn, nTrials)
% LOCAL_VALIDATE_TRIAL_INDICES Keep valid, unique trial indices in range.
if isempty(trialIndicesIn)
    trialIndices = 1:min(2, nTrials);
    return;
end
trialIndices = unique(trialIndicesIn(:)', 'stable');
trialIndices = trialIndices(trialIndices >= 1 & trialIndices <= nTrials);
if isempty(trialIndices)
    error('No valid trial indices remain after range filtering (1..%d).', nTrials);
end
end

function out = local_trial_label(trialIndices)
% LOCAL_TRIAL_LABEL Build a compact trial label for titles/filenames.
parts = cell(1, numel(trialIndices));
for i = 1:numel(trialIndices)
    parts{i} = sprintf('%d', trialIndices(i));
end
out = strjoin(parts, '-');
end

function local_ensure_dir(dirPath)
% LOCAL_ENSURE_DIR Create directory (and parents) when missing.
if ~exist(dirPath, 'dir')
    mkdir(dirPath);
end
end

function [pathsCut, cutFrame] = local_cut_postdeviant(paths, deviantOnset)
% LOCAL_CUT_POSTDEVIANT Keep the post-deviant segment of each trial.
%
% Inputs:
%   paths        : trials x features x time.
%   deviantOnset : scalar fraction in [0, 1].
%
% Outputs:
%   pathsCut     : post-deviance trial windows.
%   cutFrame     : 1-based frame index where cut starts.
if isempty(paths)
    pathsCut = paths;
    cutFrame = [];
    return;
end
if ~isscalar(deviantOnset) || deviantOnset < 0 || deviantOnset > 1
    error('deviantOnset must be a scalar in [0, 1].');
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
