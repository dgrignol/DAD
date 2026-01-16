% PIPELINE_simulation.m
%
% Purpose:
%   Run dRSA simulations on dot-motion inputs. This script loads
%   condition-specific dot-path data (dot 1 + dot 2, visual degrees),
%   prepares concatenated inputs, and executes the dRSA workflow.
%
% Example usage (from simulations/ in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   PIPELINE_simulation;
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   PIPELINE_simulation;
%
% Configuration (edit at the start of the script):
%   - participantNumber selects MovDot_SubXX.mat and the derived simulation
%     inputs (MovDot_SubXX_nondeviant/deviant.mat).
%   - inputCondition selects which condition to load for the simulation run.
%   Example configuration:
%     participantNumber = 98;
%     inputCondition = 'deviant';
%
% Inputs:
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat with:
%       dot1GreenPaths (trials × 2 × time) and dot2YellowPaths (same shape)
%
% Outputs:
%   - dRSA matrices and figures (see sections below).
%
% Key assumptions:
%   - dot1GreenPaths and dot2YellowPaths are in visual degrees.
%   - Trials have uniform frame counts within each condition.
%   - The functions in simulations/functions are on the MATLAB path.
%   - Paths are resolved relative to this script, not the MATLAB cwd.
%
% Preparatory steps (done by individual user):
% - smoothing data etc
% - concatenate runs (e.g. for movie)
% - average across runs/repetitions
% - create mask (e.g. on TPs with too many NaN repetitions)
clear all
close all

%% Resolve repo paths and add dependencies
% Data flow: script location -> repo root -> addpath for simulation helpers.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'functions'));
addpath(scriptDir);
addpath(fullfile(scriptDir, 'debug'));
addpath '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'
%% Participant configuration
% Data flow: participant number + condition -> input filenames -> simulation data.
participantNumber = 98;
inputCondition = 'nondeviant'; % 'nondeviant' or 'deviant' to match output files.
inputCondition = lower(inputCondition);
validConditions = {'nondeviant', 'deviant'};
if ~ismember(inputCondition, validConditions)
    error('inputCondition must be one of: %s', strjoin(validConditions, ', '));
end
simulationInputDir = fullfile(scriptDir, 'input');

%% Create input file needed from paths
% Data flow: MovDot_SubXX.mat -> split into condition-specific dot paths.
inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));

build_movdot_simulation_inputs(inputFile, ...
    'OutputDir', simulationInputDir, ...
    'MoveDotScript', fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m'));

%% Load condition-specific dot paths (visual degrees)
% Data flow: condition .mat -> dot1/dot2 arrays -> placeholders for simulation.
simulationInputFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_%s.mat', participantNumber, inputCondition));
% Extract subject and condition from the input filename for output naming.
[~, inputBase, ~] = fileparts(simulationInputFile);
tokens = regexp(inputBase, '^MovDot_Sub([0-9]+)_?(.*)$', 'tokens', 'once');
if isempty(tokens)
    error('Input filename does not match expected pattern: %s', inputBase);
end
subjectLabel = sprintf('sub%02d', str2double(tokens{1}));
if numel(tokens) > 1 && ~isempty(tokens{2})
    conditionLabel = tokens{2};
else
    conditionLabel = 'unknown';
end
simulationData = load(simulationInputFile);
dot1GreenPaths = simulationData.dot1GreenPaths;
dot2YellowPaths = simulationData.dot2YellowPaths;

% Prepare per-dot inputs (trial × features × time).
dot1Trials = dot1GreenPaths;
dot2Trials = dot2YellowPaths;

% % debug
% simulationInputFile
% size(dot1GreenPaths, 3)
% which plot_paths -all

plot_paths(dot1GreenPaths, 'Title', 'Dot 1 paths');
plot_paths(dot2YellowPaths, 'Title', 'Dot 2 paths');

plot_position_time_distance(dot1GreenPaths, 'Title', 'Dot 1 distance');
plot_position_time_distance(dot2YellowPaths, 'Title', 'Dot 2 distance');

%% Prepare data, models, and trial-locked subsamples
% Data flow: trial arrays -> concatenated time series -> trial-start triggers.
%[data mask] = dRSA_concatenate(data,mask) %see documentation with help dRSA_concatenate
data = dRSA_concatenate(dot1Trials); %see documentation with help dRSA_concatenate
model1 = data;
model2 = dRSA_concatenate(dot2Trials);
% Build a trigger mask at the first sample of each trial (trial-locked dRSA).
trialLen = size(dot1GreenPaths, 3);
totalTime = size(data, 2);
maskSubsampling = true(1, totalTime); % allow full-trial windows
if mod(totalTime, trialLen) ~= 0
    error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
end
trialStarts = 1:trialLen:totalTime;
maskTrigger = false(1, totalTime);
maskTrigger(trialStarts) = true;

% Triggered subsampling: each subsample spans the full trial.
opt.PreTrigger = 0;
opt.PostTrigger = trialLen - 1;
opt.spacing = 0;
opt.nSubSamples = numel(trialStarts);
opt.nIter = 1;
opt.checkRepetition = 0;

subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);
% Resample subsamples with replacement to create multiple iterations.
% resampleIterations = 50;
% resamplePercent = 50;
% subsamples = dRSA_resample_subsamples(subsamples, resampleIterations, resamplePercent);

% include illustration
% mask is maskTypes*time
% output: nSubsamples*subSampleDuration*iterations
% for later: option to provide predefined time points (e.g. for predictable vs. unpredictable time points)



%% Simulations dRSA

Y = data; % data;
model = {model1, model2};
% wrapper for many subsamples (most cases, except e.g. Ayman)


% In case we do the PCR, it is better to calculate the border outside, because otherwise we would need to recalculate it
% for each iteration
params.nIter = 1;%resampleIterations;  %how many Iterations?
params.fs = 120; %framerate or how many samples fit into 1 second
% Use an integer sample window for averaging to avoid non-integer sizes.
avgHalfWindowSamples = floor((trialLen - 1) / 2);
params.AverageTime = avgHalfWindowSamples / params.fs; %in s
params.modelToTest = [1 2];  %array of models to test
params.Var = 0.1; % how much variance? 
params.modelDistMeasure = {'euclidean', 'euclidean'};
params.neuralDistMeasure = 'euclidean'; % pdist expects a string/char or function handle.

Autocorrborder = dRSA_border(model, subsamples, params);

%For the PCR
params.dRSAtype = 'PCR';% 'PCR';% 'corr'; %
params.modeltoRegressout = {[] []};%{2 1}; % {[4 5] [1 3] [1 2]};

% Regress out autocorrelation
params.PCR.RegressAutocor = 1 
params.PCR.RessModel = 0

%the other stuff use default values. Can be changed, see documentation of PCR function

dRSA_Iter = [];

for iIter = 1:params.nIter
    CurrSubsamples = subsamples(:,:,iIter); % subsamples is nSubsamples*subSampleDuration*iterations
    dRSAma = dRSA_coreFunction(Y,model, params, ...
        'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);%
    % Y is features*time or a finished RDM of subsamples x subsamples x time
    % models: 1*nModels cell arrray (also used for regress out models, automatically regresses out other models)
    
    
    %In this funciton we have default values
    % we also have varagin: If we want, we can add the autocorrelation and already Subsamples at the end, but we dont need to
    
    
    dRSA_Iter(iIter,:,:,:) = dRSAma;
end  %of nIter

% average dRSAmats (do by summing current with previous summed)  
 % ----------------------------------QUESTION:  include this into averaging function?
dRSA = mean(dRSA_Iter ,1); % 1 = Fmean across first dim
dRSA = reshape(dRSA, size(dRSA,2), size(dRSA,3), size(dRSA,4));

%% Plot dRSA matrices
% Data flow: dRSA (time x time x model) -> per-model heatmaps.
nModels = size(dRSA, 3);
nTime = size(dRSA, 1);
figMatrices = figure('Name', 'dRSA matrices', 'NumberTitle', 'off');
tiledlayout(1, nModels, 'Padding', 'compact', 'TileSpacing', 'compact');
for iModel = 1:nModels
    nexttile;
    imagesc(dRSA(:, :, iModel));
    set(gca, 'YDir', 'normal');
    axis image;
    colorbar;
    title(sprintf('dRSA model %d', iModel));
    xlabel('Time (samples)');
    ylabel('Time (samples)');
    hold on;
    plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
    hold off;
end


% plot diagonal
params.AverageTime = 2; %in s, how much should be left and right of the zerolag middle? 
% already set above: params.fs = 120; %framerate or how many samples fit into 1 second
dRSA_diagonal = dRSA_average(dRSA, params);

figLag = figure('Name', 'dRSA diagonal', 'NumberTitle', 'off');
plot(dRSA_diagonal(1,:));

%% Save matrices and plots (300 dpi)
outputDir = fullfile('output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

dRSAtypeLabel = lower(params.dRSAtype);
resultsBase = sprintf('dRSA_%s_%s_%s', subjectLabel, conditionLabel, dRSAtypeLabel);
save(fullfile(outputDir, [resultsBase '_results.mat']), 'dRSA', 'dRSA_diagonal', 'params');

print(figMatrices, fullfile(outputDir, [resultsBase '_matrices.png']), '-dpng', '-r300');
print(figLag, fullfile(outputDir, [resultsBase '_diagonal.png']), '-dpng', '-r300');

% matrices
% lag plots
%
