% PIPELINE_BASE_NONDEVIANT_DOT1_MINIMAL
%
% Purpose:
%   Minimal dRSA script derived from:
%     - simulations/scripts/step_by_step/pipeline_base.m
%     - simulations/scripts/PE_simulation_diff.m
%
%   This version keeps only the essential pieces so you can inspect each
%   step clearly. It uses nondeviant trials only and dot1 observed streams:
%     1) dot1 observed position (x,y)
%     2) dot1 observed direction (cos(theta), sin(theta))
%
%   Input compatibility:
%   - Works with both 2-dot files (xy columns [x1 y1 x2 y2]) and
%     1-dot files from stimuli_generation_v17.m (xy columns [x y]).
%   - The script auto-detects the input format and always extracts dot1
%     (first two columns) from nondeviant trials.
%
%   It runs exactly two PCR analyses with autocorrelation regression:
%     A) Position test:
%        - Neural/data stream: dot1 observed position
%        - Model tested:       dot1 observed position
%        - Regressed out:      dot1 observed direction
%     B) Direction test:
%        - Neural/data stream: dot1 observed direction
%        - Model tested:       dot1 observed direction
%        - Regressed out:      dot1 observed position
%
% Data flow summary:
%   1) Load MovDot_SubXX.mat, detect 1-dot vs 2-dot format, and extract
%      nondeviant dot1 center-relative trials.
%   2) Build dot1 position and direction streams.
%   3) Concatenate trials with dRSA_concatenate.
%   4) Build trial-locked subsamples with dRSA_triggered_subsampling.
%   5) Compute autocorrelation borders with dRSA_border.
%   6) Run the two minimal dRSA tests and plot the two resulting matrices.
%
% Usage examples:
%   % Example 1: run from repo root
%   run('simulations/scripts/step_by_step/pipeline_base_nondeviant_dot1_minimal.m');
%
%   % Example 2: run from anywhere
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/step_by_step');
%   pipeline_base_nondeviant_dot1_minimal;
%
%   % Example 3: override participant before running the script
%   participantNumber = 73;
%   dRSAtypeToRun = 'PCR';
%   run('simulations/scripts/step_by_step/pipeline_base_nondeviant_dot1_minimal.m');
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat
%
% Outputs:
%   - In-memory variables:
%       dRSA_position_test
%       dRSA_direction_test
%   - One simple figure with 2 subplots (no lag/diagonal plot).

% Keep optional caller-provided config overrides (e.g., participantNumber).
clearvars -except participantNumber plotConcat suppressDispText plotSampleRateHz dRSAtypeToRun
close all

%% Resolve repo paths and add required functions
% Data flow: script path -> repo root -> add simulation dependencies.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(fullfile(repoRoot, 'simulations', 'functions'));
addpath(fullfile(repoRoot, 'simulations'));

%% User configuration (minimal)
% Data flow: these settings determine which participant/input/models are used.
if ~exist('participantNumber', 'var') || isempty(participantNumber)
    participantNumber = 73;
end
if ~exist('plotConcat', 'var') || isempty(plotConcat)
    plotConcat = 0;          % 0 = no concatenation diagnostics figure
end
if ~exist('suppressDispText', 'var') || isempty(suppressDispText)
    suppressDispText = 1;    % 1 = quieter console output
end
if ~exist('plotSampleRateHz', 'var') || isempty(plotSampleRateHz)
    plotSampleRateHz = 120;  % same sample rate used in PE_simulation_diff.m
end
if ~exist('dRSAtypeToRun', 'var') || isempty(dRSAtypeToRun)
    dRSAtypeToRun = 'corr';  % can be 'corr' or 'PCR'
end

%% Load nondeviant dot1 trials with auto-detection (1-dot vs 2-dot input)
% Data flow: MovDot_SubXX.mat -> detect xy format -> extract nondeviant dot1 trials.
rawInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
if ~isfile(rawInputFile)
    error('Input file not found: %s', rawInputFile);
end
[dot1Trials, inputDotMode, xyColumnCount] = local_load_dot1_nondeviant_trials(rawInputFile);
if ~suppressDispText
    fprintf('Detected %s input format (%d xy columns).\n', inputDotMode, xyColumnCount);
end

%% Build observed dot1 position and direction streams
% Data flow: dot1 trials -> position stream and direction stream -> model cell.
dataPosition = dRSA_concatenate(dot1Trials, [], plotConcat, ...
    'suppressDispText', suppressDispText);

dot1DirectionTrials = local_positions_to_direction(dot1Trials);
dataDirection = dRSA_concatenate(dot1DirectionTrials, [], plotConcat, ...
    'suppressDispText', suppressDispText);

% Two-model minimal set requested by the user.
model = {dataPosition, dataDirection};
modelNames = {'dot1 observed position', 'dot1 observed direction'};

%% Build trial-locked subsamples (same logic as PE_simulation_diff.m)
% Data flow: trial length + total time -> trigger masks -> subsample indices.
trialLen = size(dot1Trials, 3);
totalTime = size(dataPosition, 2);
[maskSubsampling, maskTrigger, opt, subsamples] = ...
    local_build_trial_subsamples(trialLen, totalTime);

%% Configure shared dRSA parameters
% Data flow: global settings -> params used in border + two dRSA runs.
paramsBase = struct();
paramsBase.modelNames = modelNames;
paramsBase.nIter = size(subsamples, 3);
paramsBase.fs = plotSampleRateHz;
avgHalfWindowSamples = floor((trialLen - 1) / 2);
paramsBase.AverageTime = avgHalfWindowSamples / paramsBase.fs;
paramsBase.Var = 0.1;
paramsBase.modelDistMeasure = {'euclidean', 'cosine'};
paramsBase.dRSAtype = dRSAtypeToRun;
paramsBase.PCR.AdditionalPCA = 1;
paramsBase.PCR.RegressAutocor = 1;
paramsBase.PCR.RessModel = 1;

% Build one shared autocorrelation border from both models.
paramsBorder = paramsBase;
paramsBorder.modelToTest = [1 2];
paramsBorder.modeltoRegressout = {2, 1};
Autocorrborder = [];
if ~strcmp(paramsBorder.dRSAtype, 'corr')
    Autocorrborder = dRSA_border(model, subsamples, paramsBorder, ...
        'suppressDispText', suppressDispText, 'plotAutocorr', 'off');
end

%% Test 1: position model on position data (regress out direction)
% Data flow: position stream + model config -> iteration loop -> mean dRSA matrix.
paramsPosition = paramsBase;
paramsPosition.modelToTest = 1;
paramsPosition.modeltoRegressout = {2};
paramsPosition.neuralDistMeasure = paramsPosition.modelDistMeasure{1};
dRSA_position_test = local_run_drsa_iterations( ...
    dataPosition, model, paramsPosition, subsamples, Autocorrborder, suppressDispText);

%% Test 2: direction model on direction data (regress out position)
% Data flow: direction stream + model config -> iteration loop -> mean dRSA matrix.
paramsDirection = paramsBase;
paramsDirection.modelToTest = 2;
paramsDirection.modeltoRegressout = {1};
paramsDirection.neuralDistMeasure = paramsDirection.modelDistMeasure{2};
dRSA_direction_test = local_run_drsa_iterations( ...
    dataDirection, model, paramsDirection, subsamples, Autocorrborder, suppressDispText);

%% Simple plot of the two resulting dRSA matrices (no lag plot)
% Data flow: dRSA outputs -> two heatmaps for quick visual inspection.
figure('Name', sprintf('Minimal dRSA (sub%02d, nondeviant, dot1)', participantNumber), ...
    'NumberTitle', 'off');

subplot(1, 2, 1);
imagesc(dRSA_position_test(:, :, 1));
set(gca, 'YDir', 'normal');
axis image;
title('Position test (regress out direction)');
xlabel('Model time (samples)');
ylabel('Neural time (samples)');
colorbar;

subplot(1, 2, 2);
imagesc(dRSA_direction_test(:, :, 1));
set(gca, 'YDir', 'normal');
axis image;
title('Direction test (regress out position)');
xlabel('Model time (samples)');
ylabel('Neural time (samples)');
colorbar;

if ~suppressDispText
    fprintf('Finished minimal dRSA run for sub%02d (nondeviant, dot1).\n', participantNumber);
    fprintf('Subsamples: %d x %d x %d\n', size(subsamples, 1), size(subsamples, 2), size(subsamples, 3));
end

%% Local helpers
function [dot1Trials, inputDotMode, xyColumnCount] = ...
        local_load_dot1_nondeviant_trials(rawInputFile)
% LOCAL_LOAD_DOT1_NONDEVIANT_TRIALS Load nondeviant dot1 paths from MovDot_SubXX.mat.
%
% Inputs:
%   rawInputFile : full path to experiment/input_files/MovDot_SubXX.mat.
%
% Outputs:
%   dot1Trials   : trials x 2 x time center-relative positions for dot1.
%   inputDotMode : 'single-dot' when xy has [x y], 'two-dot' when xy has
%                  [x1 y1 x2 y2] (or more).
%   xyColumnCount: number of xy columns detected in source trials.
%
% Data flow:
%   MAT file -> xySeqs -> valid trials -> nondeviant selection -> dot1
%   extraction -> recenter by Cfg.rectSize/2.
sourceData = load(rawInputFile);
if ~isfield(sourceData, 'xySeqs')
    error('Expected variable xySeqs in %s.', rawInputFile);
end
if ~isfield(sourceData, 'Cfg') || ~isfield(sourceData.Cfg, 'rectSize')
    error('Cfg.rectSize missing in %s.', rawInputFile);
end

allTrials = sourceData.xySeqs(:);
if isempty(allTrials) || ~isfield(allTrials, 'xy')
    error('xySeqs is empty or missing field "xy" in %s.', rawInputFile);
end

% Keep only trials with numeric non-empty xy and at least x/y columns.
validMask = arrayfun(@(s) isnumeric(s.xy) && ~isempty(s.xy) && ...
    size(s.xy, 1) > 0 && size(s.xy, 2) >= 2, allTrials);
allTrials = allTrials(validMask);
if isempty(allTrials)
    error('No valid xy trials found in %s.', rawInputFile);
end

% Use canonical nondeviant definition from the PE pipeline: condition == 0.
if isfield(allTrials, 'condition')
    condValues = arrayfun(@(s) s.condition, allTrials);
    allTrials = allTrials(condValues == 0);
    if isempty(allTrials)
        error('No nondeviant trials (condition == 0) found in %s.', rawInputFile);
    end
end

xyColumnCount = size(allTrials(1).xy, 2);
if xyColumnCount >= 4
    inputDotMode = 'two-dot';
elseif xyColumnCount >= 2
    inputDotMode = 'single-dot';
else
    error('Unsupported xy format with %d columns in %s.', xyColumnCount, rawInputFile);
end

% All selected trials must have the same number of frames.
frameCounts = arrayfun(@(s) size(s.xy, 1), allTrials);
if any(frameCounts ~= frameCounts(1))
    error('Inconsistent frame counts across nondeviant trials in %s.', rawInputFile);
end

nTrials = numel(allTrials);
nFrames = frameCounts(1);
dot1Trials = zeros(nTrials, 2, nFrames);
for iTrial = 1:nTrials
    xy = double(allTrials(iTrial).xy);
    dot1Trials(iTrial, 1, :) = reshape(xy(:, 1), 1, 1, []);
    dot1Trials(iTrial, 2, :) = reshape(xy(:, 2), 1, 1, []);
end

% Convert absolute coordinates [0..rectSize] to center-relative coordinates.
rectSize = double(sourceData.Cfg.rectSize(:)');
if numel(rectSize) ~= 2
    error('Cfg.rectSize must contain [width height] in %s.', rawInputFile);
end
centerShift = reshape(rectSize ./ 2, [1, 2, 1]);
dot1Trials = bsxfun(@minus, dot1Trials, centerShift);
end

function direction = local_positions_to_direction(paths)
% LOCAL_POSITIONS_TO_DIRECTION Convert x/y paths to unit-vector direction signals.
%
% Inputs:
%   paths : trials x 2 x time center-relative positions.
%
% Outputs:
%   direction : trials x 2 x time, columns = [cos(theta), sin(theta)].
%
% Assumption:
%   Direction at each sample is estimated from frame-to-frame displacement.
if isempty(paths)
    direction = paths;
    return;
end
if size(paths, 3) < 2
    % Edge case: with one sample there is no displacement direction.
    direction = zeros(size(paths));
    return;
end

dx = diff(paths(:, 1, :), 1, 3);
dy = diff(paths(:, 2, :), 1, 3);

% Pad first sample so output keeps the same time length as input.
dx = cat(3, dx(:, :, 1), dx);
dy = cat(3, dy(:, :, 1), dy);

angle = atan2(dy, dx);
direction = cat(2, cos(angle), sin(angle));
end

function [maskSubsampling, maskTrigger, opt, subsamples] = ...
        local_build_trial_subsamples(trialLen, totalTime)
% LOCAL_BUILD_TRIAL_SUBSAMPLES Build trial-locked subsamples with one trigger/trial.
%
% Inputs:
%   trialLen  : number of samples in each trial.
%   totalTime : concatenated stream length in samples.
%
% Outputs:
%   maskSubsampling : 1 x totalTime logical availability mask.
%   maskTrigger     : 1 x totalTime logical trigger mask (trial starts).
%   opt             : options used by dRSA_triggered_subsampling.
%   subsamples      : nSubSamples x trialLen x 1 index tensor.
if mod(totalTime, trialLen) ~= 0
    error('totalTime (%d) must be a multiple of trialLen (%d).', totalTime, trialLen);
end

maskSubsampling = true(1, totalTime);
trialStarts = 1:trialLen:totalTime;
maskTrigger = false(1, totalTime);
maskTrigger(trialStarts) = true;

opt.PreTrigger = 0;
opt.PostTrigger = trialLen - 1;
opt.spacing = 0;
opt.nSubSamples = numel(trialStarts);
opt.nIter = 1;
opt.checkRepetition = 0;

subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt, ...
    'suppressDispText', 1);
end

function dRSA = local_run_drsa_iterations( ...
        Y, model, params, subsamples, Autocorrborder, suppressDispText)
% LOCAL_RUN_DRSA_ITERATIONS Run dRSA_coreFunction across all subsampling iterations.
%
% Inputs:
%   Y              : neural/data stream (features x time).
%   model          : model cell array.
%   params         : dRSA parameters for this test.
%   subsamples     : nSubSamples x subSampleDuration x nIter indices.
%   Autocorrborder : border vector from dRSA_border (empty for corr mode).
%   suppressDispText : 0/1 pass-through to dRSA_coreFunction.
%
% Output:
%   dRSA : time x time x nModelTests averaged across iterations.
nIter = size(subsamples, 3);
dRSA_Iter = [];

for iIter = 1:nIter
    CurrSubsamples = subsamples(:, :, iIter);
    dRSAma = dRSA_coreFunction(Y, model, params, ...
        'CurrSubsamples', CurrSubsamples, ...
        'Autocorrborder', Autocorrborder, ...
        'suppressDispText', suppressDispText);
    dRSA_Iter(iIter, :, :, :) = dRSAma;
end

dRSA = mean(dRSA_Iter, 1);
dRSA = reshape(dRSA, size(dRSA, 2), size(dRSA, 3), size(dRSA, 4));
end
