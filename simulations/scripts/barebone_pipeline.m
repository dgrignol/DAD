% BAREBONE_PIPELINE
%
% Purpose:
%   Minimal showcase pipeline for running dRSA directly on real dot-path
%   stimuli. The script keeps the original barebone steps (concatenate ->
%   subsample -> dRSA -> diagonal), but replaces synthetic random inputs
%   with stimulus-driven position data.
%
% What this script runs:
%   - Neural data (Y): dot1 position (x/y over time).
%   - Models: dot1 position and dot2 position.
%   - dRSA types: both 'corr' and 'PCR' in one run.
%
% Example usage (from anywhere in MATLAB):
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/barebone_pipeline.m');
%
% Example usage (non-interactive on this machine):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch "run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/barebone_pipeline.m')"
%
% Inputs:
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat
%   - expected fields: dot1GreenPathsCenterRelative, dot2YellowPathsCenterRelative
%
% Outputs in workspace:
%   - dRSA_results.corr / dRSA_results.PCR
%       .dRSA (time x time x model)
%       .dRSA_diagonal (model x time)
%       .params (run parameters)
%
% Key assumptions:
%   - Dot paths are trial x 2 x time (x/y features).
%   - simulations/functions is available on the MATLAB path.
%   - This is a demo pipeline (quick showcase), not final simulation analysis.

%% Resolve paths and demo configuration
% Data flow: script location -> simulations folder -> default input file.
scriptDir = fileparts(mfilename('fullpath'));
simDir = fileparts(scriptDir);
addpath(scriptDir);
addpath(simDir);
addpath(fullfile(simDir, 'functions'));

% Change this file to your preferred stimulus input.
stimulusFile = fullfile(simDir, 'input', 'MovDot_Sub87_nondeviant.mat');

% Run both dRSA variants for the showcase.
dRSATypesToRun = {'corr', 'PCR'};

% Keep barebone random-subsampling step, but with safe defaults for trial data.
subsampleDuration = []; % [] => auto (<= trial length - 2 to avoid boundary mask violations)
subsampleSpacing = 10;
nSubSamples = 10;
nIter = 10;
sampleRateHz = 120;
averageTimeSec = 2;
suppressDispText = 1;

%% Load stimuli and extract dot positions
% Data flow: stimulus .mat -> trial-wise dot1/dot2 position arrays.
simulationData = load(stimulusFile);
if isfield(simulationData, 'dot1GreenPathsCenterRelative') ...
        && isfield(simulationData, 'dot2YellowPathsCenterRelative')
    dot1Trials = simulationData.dot1GreenPathsCenterRelative;
    dot2Trials = simulationData.dot2YellowPathsCenterRelative;
elseif isfield(simulationData, 'dot1GreenPaths') ...
        && isfield(simulationData, 'dot2YellowPaths')
    % Fallback for older/pre-center-relative files.
    dot1Trials = simulationData.dot1GreenPaths;
    dot2Trials = simulationData.dot2YellowPaths;
else
    error(['Stimulus file must contain dot1/dot2 paths. Expected fields: ' ...
        'dot1GreenPathsCenterRelative and dot2YellowPathsCenterRelative.']);
end

%% Prepare data and position models
% Data flow: trial x feature x time -> concatenated feature x time arrays.
% [data mask] = dRSA_concatenate(data,mask) % see documentation with help dRSA_concatenate
[data, mask] = dRSA_concatenate(dot1Trials, [], 0, 'suppressDispText', suppressDispText);
modelDot1 = dRSA_concatenate(dot1Trials, [], 0, 'suppressDispText', suppressDispText);
modelDot2 = dRSA_concatenate(dot2Trials, [], 0, 'suppressDispText', suppressDispText);
model = {modelDot1, modelDot2};

%% Optional module: create random subsamples
% Data flow: concatenation mask -> valid starts -> index windows for dRSA.
trialLen = size(dot1Trials, 3);
maskSubsampling = logical(mask.mask(1, :)); % true = available, false = unavailable
opt.SubSampleDur = subsampleDuration;
if isempty(opt.SubSampleDur)
    % Keep windows inside trial segments (boundary mask zeros at segment edges).
    opt.SubSampleDur = min(200, max(20, trialLen - 2));
end
opt.spacing = subsampleSpacing; % required spacing between subsamples (in points)
opt.nSubSamples = nSubSamples; % number of subsamples per iteration
opt.nIter = nIter; % number of iterations
opt.checkRepetition = 0; % true/false (default in helper is true)

[subsamples, ~, usedIterations] = dRSA_random_subsampling(maskSubsampling, opt);
if usedIterations < 1
    error('No valid subsampling iterations were generated. Reduce SubSampleDur/nSubSamples/spacing.');
end
subsamples = subsamples(:, :, 1:usedIterations);

% include illustration
% mask is maskTypes*time
% output: nSubsamples*subSampleDuration*iterations
% for later: option to provide predefined time points (e.g. for predictable vs. unpredictable time points)

%% dRSA (dot1 as neural; dot1 and dot2 as models)
% Data flow: neural/model inputs + params + subsamples -> dRSA matrices.
Y = data; % neural data: dot1 position

% Shared parameters for both corr and PCR runs.
paramsBase.nIter = size(subsamples, 3);
paramsBase.AverageTime = averageTimeSec; % in s
paramsBase.fs = sampleRateHz; % samples per second
paramsBase.modelToTest = [1 2]; % models to test: dot1, dot2
paramsBase.Var = 0.1; % variance threshold for border estimation
paramsBase.modelDistMeasure = {'euclidean', 'euclidean'};
paramsBase.neuralDistMeasure = 'euclidean';
paramsBase.modelNames = {'position dot1', 'position dot2'};
paramsBase.modeltoRegressout = {[2], [1]};
paramsBase.PCR.AdditionalPCA = 1;
paramsBase.PCR.RegressAutocor = 1;
paramsBase.PCR.RessModel = 1;

dRSA_results = struct();
for iType = 1:numel(dRSATypesToRun)
    params = paramsBase;
    params.dRSAtype = dRSATypesToRun{iType};

    % For PCR it is better to calculate the border outside the iteration loop.
    Autocorrborder = [];
    if ~strcmp(params.dRSAtype, 'corr')
        Autocorrborder = dRSA_border(model, subsamples, params, ...
            'suppressDispText', suppressDispText);
    end

    % Wrapper for many subsamples (most cases, except e.g. Ayman).
    dRSA_Iter = [];
    for iIter = 1:params.nIter
        CurrSubsamples = subsamples(:, :, iIter); % nSubsamples*subSampleDuration*iterations
        dRSAma = dRSA_coreFunction(Y, model, params, ...
            'CurrSubsamples', CurrSubsamples, ...
            'Autocorrborder', Autocorrborder, ...
            'suppressDispText', suppressDispText);
        dRSA_Iter(iIter, :, :, :) = dRSAma;
    end

    % Average dRSA matrices across iterations.
    dRSA = mean(dRSA_Iter, 1); % mean across first dim
    dRSA = reshape(dRSA, size(dRSA, 2), size(dRSA, 3), size(dRSA, 4));

    % Module for averaging across time.
    dRSA_diagonal = dRSA_average(dRSA, params);

    dRSA_results.(params.dRSAtype).dRSA = dRSA;
    dRSA_results.(params.dRSAtype).dRSA_diagonal = dRSA_diagonal;
    dRSA_results.(params.dRSAtype).params = params;

    %% Stats and plots (per dRSA type)
    % Data flow: diagonal outputs -> quick lag overview figure.
    figure('Name', sprintf('dRSA diagonal (%s)', params.dRSAtype), 'NumberTitle', 'off');
    plot(dRSA_diagonal(1, :), 'LineWidth', 1.5);
    hold on;
    plot(dRSA_diagonal(2, :), 'LineWidth', 1.5);
    hold off;
    xlabel('Time (samples)');
    ylabel('dRSA (diagonal)');
    legend(params.modelNames, 'Location', 'best');
    title(sprintf('barebone pipeline (%s)', params.dRSAtype));
end
