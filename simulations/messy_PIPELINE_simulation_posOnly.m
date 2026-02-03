% messy_PIPELINE_simulation_posOnly.m
%
% Purpose:
%   Run position-only dRSA simulations on dot-motion inputs. This script
%   loads dot-path data (dot 1 + dot 2, visual degrees) for both
%   nondeviant and deviant conditions, prepares concatenated inputs, and
%   executes the dRSA workflow using dot 1 position as the neural data and
%   dot 1 + dot 2 positions as the models for each condition. For the
%   deviant condition, it can also add "wannabe deviant" position models
%   from MovDot_SubXX_wannabeDev.mat (deviant-only paths without deviation).
%
% Example usage (from simulations/ in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   messy_PIPELINE_simulation_posOnly;
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   messy_PIPELINE_simulation_posOnly;
%
% Configuration (edit at the start of the script):
%   - participantNumber selects MovDot_SubXX.mat and the derived simulation
%     inputs (MovDot_SubXX_nondeviant/deviant.mat).
%   - inputConditions selects which conditions to run (default: both).
%   - shuffleDot2Trials decouples dot2 trial order from dot1 (diagnostic).
%   - resampleSubsamples toggles trial subsampling across iterations to
%     reduce memory use; resampleIterations, resampleSubsampleCount, and
%     resampleWithReplacement control the resampling behavior.
%   - output filenames include subject/condition, dRSA type, sampling rate,
%     averaging window, and shuffle/resampling flags.
%   Example configuration:
%     participantNumber = 98;
%     inputConditions = {'nondeviant', 'deviant'};
%
% Inputs:
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat with:
%       dot1GreenPaths (trials × 2 × time) and dot2YellowPaths (same shape)
%   - experiment/input_files/MovDot_SubXX_wannabeDev.mat (deviant-only)
%       used to build additional position models for deviant analyses.
%
% Outputs:
%   - simulations/output/dRSA_<subject>_<condition>_<dRSAtype>_fs<fs>_avg<avg>_<shuffle>_<resample>_results.mat
%     containing dRSA matrices, diagonal summaries, parameters, trigger masks,
%     subsamples, and input dot-trial data needed to reproduce the run.
%   - simulations/output/dRSA_<subject>_<condition>_<dRSAtype>_fs<fs>_avg<avg>_<shuffle>_<resample>_matrices.png
%   - simulations/output/dRSA_<subject>_<condition>_<dRSAtype>_fs<fs>_avg<avg>_<shuffle>_<resample>_diagonal.png
%   - Debug plots in simulations/debug (paths, time-distance, center-distance)
%     saved separately per condition.
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
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
addpath(fullfile(scriptDir, 'functions'));
addpath(scriptDir);
addpath(fullfile(scriptDir, 'debug'));
addpath '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'
%% Participant configuration
% Data flow: participant number + condition list -> input filenames -> simulation data.
participantNumber = 88;
inputConditions = {'nondeviant', 'deviant'}; % conditions to run.
inputConditions = cellfun(@lower, inputConditions, 'UniformOutput', false);
shuffleDot2Trials = false; % true to shuffle dot2 trials relative to dot1 before dRSA
resampleSubsamples = false; % true to use per-iteration trial subsets for dRSA
resampleIterations = 100; % number of resampled iterations (only if resampleSubsamples)
resampleSubsampleCount = 100; % number of triggers per iteration (only if resampleSubsamples)
resampleWithReplacement = false; % true = bootstrap-style resampling, false = unique subset per iteration
validConditions = {'nondeviant', 'deviant'};
if ~all(ismember(inputConditions, validConditions))
    error('inputConditions must be drawn from: %s', strjoin(validConditions, ', '));
end
simulationInputDir = fullfile(scriptDir, 'input');

%% Create input file needed from paths
% Data flow: MovDot_SubXX.mat -> split into condition-specific dot paths.
inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));

build_movdot_simulation_inputs(inputFile, ...
    'OutputDir', simulationInputDir, ...
    'MoveDotScript', moveDotScript);

%% Load plotting config for center-distance diagnostics
% Data flow: input file -> Cfg.rectSize -> debug plot parameters.
cfgData = load(inputFile, 'Cfg');
rectSize = [];
if isfield(cfgData, 'Cfg') && isfield(cfgData.Cfg, 'rectSize')
    % Cfg.rectSize is Config.dotRectSize (visual degrees, origin at [0 0]).
    rectSize = cfgData.Cfg.rectSize;
end

%% Run position-only dRSA per condition
% Data flow: condition list -> per-condition load -> position-only dRSA outputs.
for iCond = 1:numel(inputConditions)
    inputCondition = inputConditions{iCond};
    useWannabeModels = strcmp(inputCondition, 'deviant');

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
conditionSuffix = sprintf('(%s)', conditionLabel);
simulationData = load(simulationInputFile);
dot1GreenPaths = simulationData.dot1GreenPaths;
dot2YellowPaths = simulationData.dot2YellowPaths;

% Prepare per-dot inputs (trial × features × time).
dot1Trials = dot1GreenPaths;
dot2Trials = dot2YellowPaths; %size(dot2YellowPaths) -> (20     2   320)

%% Optional trial shuffling (dot2 relative to dot1)
% Data flow: dot2Trials -> permuted trial order -> decoupled dot1/dot2 pairing.
% Purpose: test whether dot1-dot2 structure is driven by trial-locked coupling.
shuffleOrder = []; % store permutation for reproducibility (empty when unshuffled).
if shuffleDot2Trials
    rng('shuffle'); % avoid reusing the subject-seeded RNG for the shuffle
    shuffleOrder = randperm(size(dot2Trials, 1));
    dot2Trials = dot2Trials(shuffleOrder, :, :);
end

% % debug
% simulationInputFile
% size(dot1GreenPaths, 3)
% which plot_paths -all

plot_paths(dot1GreenPaths, 'Title', sprintf('Dot 1 paths %s', conditionSuffix));
plot_paths(dot2YellowPaths, 'Title', sprintf('Dot 2 paths %s', conditionSuffix));

% plot_position_distribution(dot1GreenPaths, ...
%     'Title', sprintf('Dot 1 dot position distribution %s', conditionSuffix));
% plot_position_distribution(dot2YellowPaths, ...
%     'Title', sprintf('Dot 2 dot position distribution %s', conditionSuffix));
% 
% plot_position_time_distance(dot1GreenPaths, ...
%     'Title', sprintf('Dot 1 distance %s', conditionSuffix));
% plot_position_time_distance(dot2YellowPaths, ...
%     'Title', sprintf('Dot 2 distance %s', conditionSuffix));

%% Plot mean distance from center for each timepoint
% % Data flow: Cfg.rectSize -> center definition -> distance profiles.
% plot_position_distribution(dot1GreenPaths, ...
%     'Title', sprintf('Dot 1 distance from center %s', conditionSuffix), ...
%     'RectSize', rectSize);
% plot_position_distribution(dot2YellowPaths, ...
%     'Title', sprintf('Dot 2 distance from center %s', conditionSuffix), ...
%     'RectSize', rectSize);

%% Optional wannabe deviant position diagnostics (deviant only)
% Data flow: wannabe deviant file -> dot paths -> debug plots and model inputs.
dot1TrialsWannabe = [];
dot2TrialsWannabe = [];
wannabeSimulationFile = '';
if useWannabeModels
    wannabeInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
        sprintf('MovDot_Sub%02d_wannabeDev.mat', participantNumber));
    build_movdot_simulation_inputs(wannabeInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'AllowMissingConditions', true);
    wannabeSimulationFile = fullfile(simulationInputDir, ...
        sprintf('MovDot_Sub%02d_wannabeDev_deviant.mat', participantNumber));
    if ~isfile(wannabeSimulationFile)
        error('Wannabe deviant input file not found: %s', wannabeSimulationFile);
    end
    wannabeData = load(wannabeSimulationFile);
    if ~isfield(wannabeData, 'dot1GreenPaths') || ~isfield(wannabeData, 'dot2YellowPaths')
        error('Wannabe deviant file lacks expected dot paths: %s', wannabeSimulationFile);
    end
    dot1TrialsWannabe = wannabeData.dot1GreenPaths;
    dot2TrialsWannabe = wannabeData.dot2YellowPaths;
    if ~isequal(size(dot1TrialsWannabe), size(dot1Trials)) || ...
            ~isequal(size(dot2TrialsWannabe), size(dot2Trials))
        error(['Wannabe deviant trials do not match deviant trials size. ' ...
            'Expected %s, got %s.'], mat2str(size(dot1Trials)), ...
            mat2str(size(dot1TrialsWannabe)));
    end

    plot_paths(dot1TrialsWannabe, ...
        'Title', sprintf('Dot 1 wannabe paths %s', conditionSuffix));
    plot_paths(dot2TrialsWannabe, ...
        'Title', sprintf('Dot 2 wannabe paths %s', conditionSuffix));

    % plot_position_distribution(dot1TrialsWannabe, ...
    %     'Title', sprintf('Dot 1 wannabe position distribution %s', conditionSuffix));
    % plot_position_distribution(dot2TrialsWannabe, ...
    %     'Title', sprintf('Dot 2 wannabe position distribution %s', conditionSuffix));
    %
    % plot_position_time_distance(dot1TrialsWannabe, ...
    %     'Title', sprintf('Dot 1 wannabe distance %s', conditionSuffix));
    % plot_position_time_distance(dot2TrialsWannabe, ...
    %     'Title', sprintf('Dot 2 wannabe distance %s', conditionSuffix));
    %
    % plot_position_distribution(dot1TrialsWannabe, ...
    %     'Title', sprintf('Dot 1 wannabe distance from center %s', conditionSuffix), ...
    %     'RectSize', rectSize);
    % plot_position_distribution(dot2TrialsWannabe, ...
    %     'Title', sprintf('Dot 2 wannabe distance from center %s', conditionSuffix), ...
    %     'RectSize', rectSize);
end

%% Prepare data, models, and trial-locked subsamples
% Data flow: trial arrays -> per-dot position models -> concatenated time series.
% Position models use raw x/y paths for dot 1 and dot 2.
% [data mask] = dRSA_concatenate(data,mask) %see documentation with help dRSA_concatenate
dataPosition = dRSA_concatenate(dot1Trials); %see documentation with help dRSA_concatenate

% Position models (per dot).
modelPositionDot1 = dataPosition;
modelPositionDot2 = dRSA_concatenate(dot2Trials);
modelPositionDot1Wannabe = [];
modelPositionDot2Wannabe = [];
if useWannabeModels
    modelPositionDot1Wannabe = dRSA_concatenate(dot1TrialsWannabe);
    modelPositionDot2Wannabe = dRSA_concatenate(dot2TrialsWannabe);
end

% Build a trigger mask at the first sample of each trial (trial-locked dRSA)
trialLen = size(dot1Trials, 3);
totalTime = size(dataPosition, 2);
maskSubsampling = true(1, totalTime); % allow full-trial windows
if mod(totalTime, trialLen) ~= 0
    error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
end
trialStarts = 1:trialLen:totalTime;
maskTrigger = false(1, totalTime);
maskTrigger(trialStarts) = true;

% Triggered subsampling: each subsample spans the full trial
opt.PreTrigger = 0;
opt.PostTrigger = trialLen - 1;
opt.spacing = 0;
opt.nSubSamples = numel(trialStarts);
opt.nIter = 1;
opt.checkRepetition = 0;

subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);
%% Optional subsample resampling (memory control)
% Data flow: full trial subsamples -> (optional) resampled subsets -> iteration tensor.
if resampleSubsamples
    baseSubsamples = subsamples(:, :, 1);
    baseCount = size(baseSubsamples, 1);
    if resampleSubsampleCount > baseCount
        error('resampleSubsampleCount (%d) exceeds available subsamples (%d).', ...
            resampleSubsampleCount, baseCount);
    end
    if resampleWithReplacement
        resamplePercent = 100 * (resampleSubsampleCount / baseCount);
        subsamples = dRSA_resample_subsamples(baseSubsamples, resampleIterations, resamplePercent);
    else
        subsamples = zeros(resampleSubsampleCount, size(baseSubsamples, 2), resampleIterations);
        for iIter = 1:resampleIterations
            drawIdx = randperm(baseCount, resampleSubsampleCount);
            subsamples(:, :, iIter) = baseSubsamples(drawIdx, :);
        end
    end
end

% include illustration
% mask is maskTypes*time
% output: nSubsamples*subSampleDuration*iterations
% for later: option to provide predefined time points (e.g. for predictable vs. unpredictable time points)



%% Simulations dRSA (position-only data)
% Data flow: dot1 position data -> dRSA against position models.
Y = dataPosition; % position data from dot 1
if useWannabeModels
    model = {modelPositionDot1, modelPositionDot2, ...
        modelPositionDot1Wannabe, modelPositionDot2Wannabe};
    modelNames = {'position dot1', 'position dot2', ...
        'wannabe position dot1', 'wannabe position dot2'};
else
    model = {modelPositionDot1, modelPositionDot2};
    modelNames = {'position dot1', 'position dot2'};
end

% In case we do the PCR, it is better to calculate the border outside, because otherwise we would need to recalculate it
% for each iteration
params.nIter = size(subsamples, 3);  % how many iterations?
params.fs = 120; %framerate or how many samples fit into 1 second
% Use an integer sample window for averaging to avoid non-integer sizes.
avgHalfWindowSamples = floor((trialLen - 1) / 2);
params.AverageTime = avgHalfWindowSamples / params.fs; %in s
if useWannabeModels
    params.modelToTest = [1 2 3 4];  %array of models to test
else
    params.modelToTest = [1 2];  %array of models to test
end
params.Var = 0.1; % how much variance? 
if useWannabeModels
    params.modelDistMeasure = {'euclidean', 'euclidean', 'euclidean', 'euclidean'};
else
    params.modelDistMeasure = {'euclidean', 'euclidean'};
end
params.neuralDistMeasure = 'euclidean'; % pdist expects a string/char or function handle.

%For the PCR
params.dRSAtype = 'corr';% 'PCR';% 'corr'; %
if useWannabeModels
    params.modeltoRegressout = {[2 3 4] [1 3 4] [1 2 4] [1 2 3]};
else
    params.modeltoRegressout = {[2] [1]};
end
params.PCR.AdditionalPCA = 1;
% params.normalizazion = 'Rescale';

% Regress out autocorrelation
params.PCR.RegressAutocor = 0; % 1 to get the autocorr
params.PCR.RessModel = 1;

% Autocorrelation border is only needed for non-corr dRSA.
Autocorrborder = []; % empty signals "no border" to downstream code paths.
if ~strcmp(params.dRSAtype, 'corr') % Autocorrelation not necessary for 'corr' type.
    Autocorrborder = dRSA_border(model, subsamples, params, ...
        'suppressDispText', suppressDispText);
end
%the other stuff use default values. Can be changed, see documentation of PCR function

%% Run dRSA for position data.
dRSA_Iter = [];
for iIter = 1:params.nIter
    CurrSubsamples = subsamples(:,:,iIter); % subsamples is nSubsamples*subSampleDuration*iterations
    dRSAma = dRSA_coreFunction(Y, model, params, ...
        'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);%
    % Y is features*time or a finished RDM of subsamples x subsamples x time
    % models: 1*nModels cell arrray (also used for regress out models, automatically regresses out other models)
    
    %In this funciton we have default values
    % we also have varagin: If we want, we can add the autocorrelation and already Subsamples at the end, but we dont need to
    
    dRSA_Iter(iIter,:,:,:) = dRSAma;
end  %of nIter

% average dRSAmats (do by summing current with previous summed)  
% ----------------------------------QUESTION:  include this into averaging function?
dRSA_position = mean(dRSA_Iter ,1); % 1 = Fmean across first dim
dRSA_position = reshape(dRSA_position, size(dRSA_position,2), size(dRSA_position,3), size(dRSA_position,4));

% Preserve parameters used for dRSA before reusing them for diagonal averaging.
paramsCore = params;

%% Plot dRSA matrices (position data)
% Data flow: dRSA per model (time x time x model) -> grid of heatmaps.
nModels = size(dRSA_position, 3);
nTime = size(dRSA_position, 1);
rowName = 'Data: position dot1';
figMatrices = figure('Name', 'dRSA matrices', 'NumberTitle', 'off');
tiledlayout(1, nModels, 'Padding', 'compact', 'TileSpacing', 'compact');
for iModel = 1:nModels
    nexttile(iModel);
    imagesc(dRSA_position(:, :, iModel));
    set(gca, 'YDir', 'normal');
    axis image;
    colorbar;
    title(sprintf('dRSA %s', modelNames{iModel}));
    xlabel(sprintf('Time in %s (samples)', modelNames{iModel}));
    ylabel(sprintf('Time in %s (samples)', rowName));
    hold on;
    plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
    hold off;
end

%% Plot diagonal
% Data flow: dRSA matrices -> diagonal averages -> summary lines.
paramsDiagonal = paramsCore;
paramsDiagonal.AverageTime = 2; %in s, how much should be left and right of the zerolag middle?
% already set above: params.fs = 120; %framerate or how many samples fit into 1 second
dRSA_diagonal_position = dRSA_average(dRSA_position, paramsDiagonal);

figLag = figure('Name', 'dRSA diagonal', 'NumberTitle', 'off');
hold on;
for iModel = 1:nModels
    plot(dRSA_diagonal_position(iModel, :), 'LineWidth', 1.5);
end
hold off;
xlabel('Time (samples)');
ylabel('dRSA (diagonal)');
legend(modelNames, 'Location', 'best');

%% Save matrices, plots, and reproducibility metadata (300 dpi)
% Data flow: dRSA outputs + run metadata -> descriptive filenames -> .mat/.png outputs.
outputDir = fullfile(scriptDir, 'output');
subjectOutputDir = fullfile(outputDir, subjectLabel);
if ~exist(subjectOutputDir, 'dir')
    mkdir(subjectOutputDir);
end

dRSAtypeLabel = lower(paramsCore.dRSAtype);
shuffleTag = 'noShufDot2';
if shuffleDot2Trials
    shuffleTag = 'shufDot2';
end
resampleTag = 'rsNone';
if resampleSubsamples
    resampleMode = 'sub';
    if resampleWithReplacement
        resampleMode = 'boot';
    end
    resampleTag = sprintf('rs%d_i%d_%s', ...
        resampleSubsampleCount, resampleIterations, resampleMode);
end
runTagParts = {subjectLabel, conditionLabel, dRSAtypeLabel, ...
    sprintf('fs%d', paramsCore.fs), sprintf('avg%d', avgHalfWindowSamples), ...
    shuffleTag, resampleTag};
resultsBase = strjoin([{'dRSA'}, runTagParts], '_');
resultsBase = regexprep(resultsBase, '\s+', '');

repro = struct();
repro.timestamp = datestr(now, 'yyyymmddTHHMMSS');
repro.scriptPath = mfilename('fullpath');
repro.repoRoot = repoRoot;
repro.matlabVersion = version;
repro.participantNumber = participantNumber;
repro.subjectLabel = subjectLabel;
repro.inputCondition = inputCondition;
repro.conditionLabel = conditionLabel;
repro.simulationInputFile = simulationInputFile;
repro.rawInputFile = inputFile;
repro.shuffleDot2Trials = shuffleDot2Trials;
repro.shuffleOrder = shuffleOrder;
repro.resampleSubsamples = resampleSubsamples;
repro.resampleIterations = resampleIterations;
repro.resampleSubsampleCount = resampleSubsampleCount;
repro.resampleWithReplacement = resampleWithReplacement;
repro.trialLen = trialLen;
repro.totalTime = totalTime;
repro.nTrials = size(dot1Trials, 1);
repro.nModels = nModels;
repro.modelNames = modelNames;
repro.paramsCore = paramsCore;
repro.paramsDiagonal = paramsDiagonal;
repro.triggerOptions = opt;
repro.maskSubsampling = maskSubsampling;
repro.maskTrigger = maskTrigger;
repro.subsamples = subsamples;
repro.rectSize = rectSize;
repro.dot1Trials = dot1Trials;
repro.dot2Trials = dot2Trials;
repro.dot1TrialsWannabe = dot1TrialsWannabe;
repro.dot2TrialsWannabe = dot2TrialsWannabe;
repro.dataPosition = dataPosition;
repro.models = model;
repro.autocorrBorder = Autocorrborder;
repro.wannabeSimulationFile = wannabeSimulationFile;

save(fullfile(subjectOutputDir, [resultsBase '_results.mat']), ...
    'dRSA_position', ...
    'dRSA_diagonal_position', ...
    'params', 'paramsDiagonal', 'repro');

print(figMatrices, fullfile(subjectOutputDir, [resultsBase '_matrices.png']), '-dpng', '-r300');
print(figLag, fullfile(subjectOutputDir, [resultsBase '_diagonal.png']), '-dpng', '-r300');

end

% matrices
% lag plots
%
