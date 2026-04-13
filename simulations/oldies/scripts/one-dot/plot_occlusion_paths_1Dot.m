% plot_occlusion_paths_1Dot.m
%
% Purpose:
%   Minimal path-inspection script for one-dot occlusion datasets.
%
%   Given a subject number, this script loads:
%     - always_visible
%     - occluded_nondeviant
%     - occluded_deviant
%
%   from experiment/input_files/MovDot_SubXX.mat, then plots the same sampled
%   trial identities across all three conditions using shared x/y limits.
%
% Usage example (interactive):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/one-dot');
%   plot_occlusion_paths_1Dot
%
% Usage example (non-interactive):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/one-dot'); participantNumber=70; numTrialsToPlot=12; randomSeed=7; run('plot_occlusion_paths_1Dot.m');"
%
% Inputs (override before run if needed):
%   participantNumber : subject ID (default 70)
%   numTrialsToPlot   : max shared trials plotted per condition (default 20)
%   randomSeed        : RNG seed for deterministic shared trial selection (default 70)
%   figVisibility     : 'on' or 'off' (default 'on')
%
% Outputs:
%   - On-screen figure with 3 panels (always/nondeviant/deviant), shared limits.
%   - Saved PNG under:
%       simulations/output/SubXX_oneDot_occlusion/paths/
%
% Key assumptions:
%   - MovDot_SubXX.mat contains one-dot trajectories (xy columns [x y]).
%   - condition_label exists and includes the three occlusion labels.
%   - sequence metadata is used when available to align trial identities;
%     otherwise, shared index fallback is used.

%% Resolve caller overrides and paths
% Data flow: caller overrides -> defaults -> input/output paths and addpath setup.
clearvars -except participantNumber numTrialsToPlot randomSeed figVisibility;
close all;

if ~exist('participantNumber', 'var')
    participantNumber = 70;
end
if ~exist('numTrialsToPlot', 'var')
    numTrialsToPlot = 20;
end
if ~exist('randomSeed', 'var')
    randomSeed = 70;
end
if ~exist('figVisibility', 'var') || isempty(figVisibility)
    figVisibility = 'on';
end

participantNumber = round(double(participantNumber));
numTrialsToPlot = round(double(numTrialsToPlot));
randomSeed = round(double(randomSeed));
if participantNumber < 0
    error('participantNumber must be >= 0.');
end
if numTrialsToPlot < 1
    error('numTrialsToPlot must be >= 1.');
end
if ~ismember(lower(char(figVisibility)), {'on', 'off'})
    error('figVisibility must be ''on'' or ''off''.');
end

rng(randomSeed);

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(fullfile(repoRoot, 'simulations', 'debug'));

inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
if ~isfile(inputFile)
    error('Input file not found: %s', inputFile);
end

outputDir = fullfile(repoRoot, 'simulations', 'output', ...
    sprintf('Sub%02d_oneDot_occlusion', participantNumber), 'paths');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load one-dot paths for each condition
% Data flow: input MAT -> center-relative one-dot arrays + sequence metadata.
[alwaysPaths, alwaysMeta] = local_load_condition_paths(inputFile, 'always_visible');
[nondevPaths, nondevMeta] = local_load_condition_paths(inputFile, 'occluded_nondeviant');
[devPaths, devMeta] = local_load_condition_paths(inputFile, 'occluded_deviant');

%% Sample same trial identities across all three conditions
% Data flow: condition arrays + sequence vectors -> shared sample tensors.
[alwaysSample, nondevSample, devSample] = local_sample_shared_trials( ...
    alwaysPaths, alwaysMeta.sequenceValues, ...
    nondevPaths, nondevMeta.sequenceValues, ...
    devPaths, devMeta.sequenceValues, ...
    numTrialsToPlot);

%% Compute shared axis limits and render figure
% Data flow: sampled arrays -> pooled limits -> fixed-axis 3-panel figure.
axisLimits = local_compute_shared_axis_limits({alwaysSample, nondevSample, devSample});

fig = figure('Name', sprintf('Occlusion paths sub%02d', participantNumber), ...
    'NumberTitle', 'off', 'Visible', figVisibility);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile(1);
plot_paths(alwaysSample, 'ParentAxes', ax1, 'SavePlot', false, ...
    'Title', 'Always visible: dot1', 'SampleRateHz', alwaysMeta.fps, ...
    'AxisLimits', axisLimits);

ax2 = nexttile(2);
plot_paths(nondevSample, 'ParentAxes', ax2, 'SavePlot', false, ...
    'Title', 'Occluded nondeviant: dot1', 'SampleRateHz', nondevMeta.fps, ...
    'AxisLimits', axisLimits);

ax3 = nexttile(3);
plot_paths(devSample, 'ParentAxes', ax3, 'SavePlot', false, ...
    'Title', 'Occluded deviant: dot1', 'SampleRateHz', devMeta.fps, ...
    'AxisLimits', axisLimits);

sgtitle(sprintf('Sub%02d shared-trial paths (n=%d)', participantNumber, size(alwaysSample, 1)));

outFile = fullfile(outputDir, sprintf('sub%02d_shared_paths_n%02d.png', ...
    participantNumber, size(alwaysSample, 1)));
print(fig, outFile, '-dpng', '-r300');
fprintf('Saved path comparison figure: %s\n', outFile);

%% Local helpers
function [paths, meta] = local_load_condition_paths(inputFile, conditionLabel)
% LOCAL_LOAD_CONDITION_PATHS Load one-dot center-relative paths for one condition.
S = load(inputFile, 'xySeqs', 'Cfg');
if ~isfield(S, 'xySeqs') || isempty(S.xySeqs)
    error('Input has no xySeqs: %s', inputFile);
end
if ~isfield(S, 'Cfg') || ~isfield(S.Cfg, 'rectSize')
    error('Input missing Cfg.rectSize: %s', inputFile);
end
if ~isfield(S.Cfg, 'fps') || isempty(S.Cfg.fps)
    error('Input missing Cfg.fps: %s', inputFile);
end

trials = S.xySeqs(:);
if ~isfield(trials, 'condition_label')
    error('Trials in %s must include condition_label for occlusion plotting.', inputFile);
end
labels = arrayfun(@(s) string(s.condition_label), trials, 'UniformOutput', true);

mask = strcmp(labels, conditionLabel);
condTrials = trials(mask);
if isempty(condTrials)
    error('No trials found for condition %s in %s.', conditionLabel, inputFile);
end

frameCounts = arrayfun(@(s) size(s.xy, 1), condTrials);
if any(frameCounts ~= frameCounts(1))
    error('Inconsistent frame counts in %s for %s.', conditionLabel, inputFile);
end

nTrials = numel(condTrials);
nFrames = frameCounts(1);
paths = zeros(nTrials, 2, nFrames);
for i = 1:nTrials
    xy = double(condTrials(i).xy);
    if size(xy, 2) ~= 2
        error('One-dot plot helper expects xy columns [x y].');
    end
    paths(i, 1, :) = reshape(xy(:, 1), 1, 1, []);
    paths(i, 2, :) = reshape(xy(:, 2), 1, 1, []);
end

rectSize = double(S.Cfg.rectSize(:)');
centerShift = reshape(rectSize ./ 2, [1, 2, 1]);
paths = bsxfun(@minus, paths, centerShift);

sequenceValues = nan(nTrials, 1);
for i = 1:nTrials
    if isfield(condTrials(i), 'sequence') && ~isempty(condTrials(i).sequence) ...
            && isfinite(double(condTrials(i).sequence))
        sequenceValues(i) = round(double(condTrials(i).sequence));
    end
end

meta = struct( ...
    'fps', double(S.Cfg.fps), ...
    'sequenceValues', sequenceValues);
end

function [alwaysSample, nondevSample, devSample] = local_sample_shared_trials( ...
        alwaysPaths, seqAlways, nondevPaths, seqNondev, devPaths, seqDev, sampleCount)
% LOCAL_SAMPLE_SHARED_TRIALS Sample same trial identities across all conditions.
nAlways = size(alwaysPaths, 1);
nNondev = size(nondevPaths, 1);
nDev = size(devPaths, 1);

useSeq = ...
    numel(seqAlways) == nAlways && numel(seqNondev) == nNondev && numel(seqDev) == nDev && ...
    all(isfinite(seqAlways)) && all(isfinite(seqNondev)) && all(isfinite(seqDev));

if useSeq
    commonSeq = intersect(intersect(unique(seqAlways), unique(seqNondev)), unique(seqDev));
    if ~isempty(commonSeq)
        nPick = min(sampleCount, numel(commonSeq));
        drawSeq = commonSeq(randperm(numel(commonSeq), nPick));
        idxAlways = local_seq_to_indices(seqAlways, drawSeq);
        idxNondev = local_seq_to_indices(seqNondev, drawSeq);
        idxDev = local_seq_to_indices(seqDev, drawSeq);
        if all(idxAlways > 0) && all(idxNondev > 0) && all(idxDev > 0)
            alwaysSample = alwaysPaths(idxAlways, :, :);
            nondevSample = nondevPaths(idxNondev, :, :);
            devSample = devPaths(idxDev, :, :);
            return;
        end
    end
end

nCommon = min([nAlways, nNondev, nDev]);
if nCommon < 1
    error('No shared trials available for plotting.');
end
nPick = min(sampleCount, nCommon);
drawIdx = randperm(nCommon, nPick);
alwaysSample = alwaysPaths(drawIdx, :, :);
nondevSample = nondevPaths(drawIdx, :, :);
devSample = devPaths(drawIdx, :, :);
end

function idx = local_seq_to_indices(sequenceVector, sequenceList)
% LOCAL_SEQ_TO_INDICES Resolve indices for requested sequence IDs.
idx = zeros(1, numel(sequenceList));
for i = 1:numel(sequenceList)
    match = find(sequenceVector == sequenceList(i), 1, 'first');
    if isempty(match)
        idx(i) = 0;
    else
        idx(i) = match;
    end
end
end

function axisLimits = local_compute_shared_axis_limits(pathCell)
% LOCAL_COMPUTE_SHARED_AXIS_LIMITS Compute shared [xmin xmax ymin ymax].
xMin = inf;
xMax = -inf;
yMin = inf;
yMax = -inf;
for i = 1:numel(pathCell)
    p = pathCell{i};
    if isempty(p)
        continue;
    end
    xVals = reshape(p(:, 1, :), [], 1);
    yVals = reshape(p(:, 2, :), [], 1);
    xMin = min(xMin, min(xVals));
    xMax = max(xMax, max(xVals));
    yMin = min(yMin, min(yVals));
    yMax = max(yMax, max(yVals));
end

if ~isfinite(xMin) || ~isfinite(xMax) || ~isfinite(yMin) || ~isfinite(yMax)
    axisLimits = [];
    return;
end

xRange = max(eps, xMax - xMin);
yRange = max(eps, yMax - yMin);
xPad = 0.05 * xRange;
yPad = 0.05 * yRange;
axisLimits = [xMin - xPad, xMax + xPad, yMin - yPad, yMax + yPad];
end
