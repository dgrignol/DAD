% PIPELINE_simulation_posOnly.m
%
% Purpose:
%   Run position-only dRSA simulations on dot-motion inputs. This script
%   loads center-relative dot-path data (dot 1 + dot 2, visual degrees) for
%   both nondeviant and deviant conditions, prepares concatenated inputs,
%   and executes the dRSA workflow using dot 1 position as the neural data
%   and dot 1 + dot 2 positions as the models for each condition. For the
%   deviant condition, it can also add "predicted deviant" position models
%   from MovDot_SubXX_predicted.mat (deviant-only paths without deviation).
%
% Example usage (from simulations/scripts in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   PIPELINE_simulation_posOnly;
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   PIPELINE_simulation_posOnly;
%
% Configuration (edit at the start of the script):
%   - participantNumber selects MovDot_SubXX.mat and the derived simulation
%     inputs (MovDot_SubXX_nondeviant/deviant.mat).
%   - inputConditions selects which conditions to run (default: both).
%   - shuffleDot2Trials decouples dot2 trial order from dot1 (diagnostic).
%   - centerRelativeRectSize overrides Cfg.rectSize when recentering.
%   - resampleSubsamples toggles trial subsampling across iterations to
%     reduce memory use; resampleIterations, resampleSubsampleCount, and
%     resampleWithReplacement control the resampling behavior.
%   - saveOutputs toggles whether .mat/.png outputs are written (true by default).
%   - output filenames strip any "run_default" token from condition labels
%     (e.g., run_default_deviant -> deviant) for cleaner naming.
%   - output filenames include subject/condition, dRSA type, and (only when
%     enabled) shuffle/resampling flags.
%   - output directories are organized by subject/condition/dRSAtype with
%     subfolders for results, matrices, and diagonals; dot-path overviews
%     are stored under subject-level paths/.
%   Example configuration:
%     participantNumber = 98;
%     inputConditions = {'nondeviant', 'deviant'};
%
% Inputs:
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat with:
%       dot1GreenPathsCenterRelative and dot2YellowPathsCenterRelative
%   - experiment/input_files/MovDot_SubXX_predicted.mat (deviant-only)
%       used to build additional position models for deviant analyses.
%
% Outputs (organized, when saveOutputs = true):
%   - simulations/output/subXX/paths/paths_subXX_allConditions_all_conditions.png
%   - simulations/output/subXX/<condition>/<dRSAtype>/results/
%       dRSA_<subject>_<condition>_<dRSAtype>_[<shuffle>][_<resample>]_results.mat
%   - simulations/output/subXX/<condition>/<dRSAtype>/diagonal/
%       dRSA_<subject>_<condition>_<dRSAtype>_[<shuffle>][_<resample>]_diagonal.png
%   - simulations/output/subXX/<condition>/<dRSAtype>/matrices/
%       dRSA_<subject>_<condition>_<dRSAtype>_[<shuffle>][_<resample>]_matrices.png
%   - Debug plots in simulations/debug (paths, time-distance, center-distance).
%
% Key assumptions:
%   - dot1GreenPathsCenterRelative and dot2YellowPathsCenterRelative are in visual degrees.
%   - Paths are centered so fixation is (0,0).
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
% Data flow: script location -> simulations dir -> repo root -> addpath for helpers.
scriptDir = fileparts(mfilename('fullpath'));
simDir = fileparts(scriptDir);
repoRoot = fileparts(simDir);
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
addpath(scriptDir);
addpath(simDir);
addpath(fullfile(simDir, 'functions'));
addpath(fullfile(simDir, 'debug'));
addpath '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD'
%% Participant configuration
% Data flow: participant number + condition list -> input filenames -> simulation data.
participantNumber = 87;
inputConditions = {'nondeviant', 'deviant'}; % conditions to run.
inputConditions = cellfun(@lower, inputConditions, 'UniformOutput', false);
shuffleDot2Trials = false; % true to shuffle dot2 trials relative to dot1 before dRSA
centerRelativeRectSize = []; % optional [width height] in degrees if Cfg.rectSize missing
resampleSubsamples = false; % true to use per-iteration trial subsets for dRSA
resampleIterations = 100; % number of resampled iterations (only if resampleSubsamples)
resampleSubsampleCount = 100; % number of triggers per iteration (only if resampleSubsamples)
resampleWithReplacement = false; % true = bootstrap-style resampling, false = unique subset per iteration
saveOutputs = true; % true to save outputs into simulations/output; false to compute only.
validConditions = {'nondeviant', 'deviant'};
if ~all(ismember(inputConditions, validConditions))
    error('inputConditions must be drawn from: %s', strjoin(validConditions, ', '));
end
simulationInputDir = fullfile(simDir, 'input');

%% Create input file needed from paths
% Data flow: MovDot_SubXX.mat -> split into condition-specific dot paths.
inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));

build_movdot_simulation_inputs(inputFile, ...
    'OutputDir', simulationInputDir, ...
    'MoveDotScript', moveDotScript, ...
    'RectSize', centerRelativeRectSize);

%% Load plotting config for center-distance diagnostics
% Data flow: input file -> Cfg.rectSize -> debug plot parameters.
cfgData = load(inputFile, 'Cfg');
rectSize = [];
if isfield(cfgData, 'Cfg') && isfield(cfgData.Cfg, 'rectSize')
    % Cfg.rectSize is Config.dotRectSize (visual degrees, origin at [0 0]).
    rectSize = cfgData.Cfg.rectSize;
end

%% Plot dot paths for all conditions in a single figure
% Data flow: per-condition .mat files -> dot paths -> plot_paths_wrap figure.
n = 20; % number of paths to plot per condition
plotSampleRateHz = 120; % keep in sync with params.fs for seconds-based colorbars.
allConditions = {'nondeviant', 'deviant'};
pathsByCondition = struct();
for iCond = 1:numel(allConditions)
    condName = allConditions{iCond};
    condFile = fullfile(simulationInputDir, ...
        sprintf('MovDot_Sub%02d_%s.mat', participantNumber, condName));
    condData = load(condFile);
    if ~isfield(condData, 'dot1GreenPathsCenterRelative') ...
            || ~isfield(condData, 'dot2YellowPathsCenterRelative')
        error(['Center-relative paths not found for condition: %s. ' ...
            'Rebuild inputs or provide centerRelativeRectSize.'], condName);
    end
    pathsByCondition.(condName).dot1 = condData.dot1GreenPathsCenterRelative;
    pathsByCondition.(condName).dot2 = condData.dot2YellowPathsCenterRelative;
end

figAllPaths = plot_paths_wrap( ...
    pathsByCondition.nondeviant.dot1, ...
    pathsByCondition.nondeviant.dot2, ...
    pathsByCondition.deviant.dot1, ...
    pathsByCondition.deviant.dot2, ...
    'SampleCount', n, ...
    'SampleRateHz', plotSampleRateHz, ...
    'Title', sprintf('Dot paths (sub%02d, all conditions)', participantNumber));

%% Run position-only dRSA per condition
% Data flow: condition list -> per-condition load -> position-only dRSA outputs.
for iCond = 1:numel(inputConditions)
    inputCondition = inputConditions{iCond};
    usePredictedModels = strcmp(inputCondition, 'deviant');

    %% Load condition-specific dot paths (center-relative visual degrees)
    % Data flow: condition .mat -> center-relative dot arrays -> placeholders for simulation.
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
    if ~isfield(simulationData, 'dot1GreenPathsCenterRelative') ...
            || ~isfield(simulationData, 'dot2YellowPathsCenterRelative')
        error(['Center-relative paths not found. Rebuild inputs or provide ' ...
            'centerRelativeRectSize if Cfg.rectSize is missing.']);
    end
    dot1GreenPaths = simulationData.dot1GreenPathsCenterRelative;
    dot2YellowPaths = simulationData.dot2YellowPathsCenterRelative;

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

    %% Optional predicted deviant position diagnostics (deviant only)
    % Data flow: predicted deviant file -> center-relative paths -> debug plots and model inputs.
dot1TrialsPredicted = [];
dot2TrialsPredicted = [];
predictedSimulationFile = '';
if usePredictedModels
    predictedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
        sprintf('MovDot_Sub%02d_predicted.mat', participantNumber));
    build_movdot_simulation_inputs(predictedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'AllowMissingConditions', true);
    predictedSimulationFile = fullfile(simulationInputDir, ...
        sprintf('MovDot_Sub%02d_predicted_deviant.mat', participantNumber));
    if ~isfile(predictedSimulationFile)
        error('Predicted deviant input file not found: %s', predictedSimulationFile);
    end
        predictedData = load(predictedSimulationFile);
        if ~isfield(predictedData, 'dot1GreenPathsCenterRelative') ...
                || ~isfield(predictedData, 'dot2YellowPathsCenterRelative')
            error('Predicted deviant file lacks center-relative paths: %s', predictedSimulationFile);
        end
        dot1TrialsPredicted = predictedData.dot1GreenPathsCenterRelative;
        dot2TrialsPredicted = predictedData.dot2YellowPathsCenterRelative;
    if ~isequal(size(dot1TrialsPredicted), size(dot1Trials)) || ...
            ~isequal(size(dot2TrialsPredicted), size(dot2Trials))
        error(['Predicted deviant trials do not match deviant trials size. ' ...
            'Expected %s, got %s.'], mat2str(size(dot1Trials)), ...
            mat2str(size(dot1TrialsPredicted)));
    end

    % plot_position_distribution(dot1TrialsPredicted, ...
    %     'Title', sprintf('Dot 1 predicted position distribution %s', conditionSuffix));
    % plot_position_distribution(dot2TrialsPredicted, ...
    %     'Title', sprintf('Dot 2 predicted position distribution %s', conditionSuffix));
    %
    % plot_position_time_distance(dot1TrialsPredicted, ...
    %     'Title', sprintf('Dot 1 predicted distance %s', conditionSuffix));
    % plot_position_time_distance(dot2TrialsPredicted, ...
    %     'Title', sprintf('Dot 2 predicted distance %s', conditionSuffix));
    %
    % plot_position_distribution(dot1TrialsPredicted, ...
    %     'Title', sprintf('Dot 1 predicted distance from center %s', conditionSuffix), ...
    %     'RectSize', rectSize);
    % plot_position_distribution(dot2TrialsPredicted, ...
    %     'Title', sprintf('Dot 2 predicted distance from center %s', conditionSuffix), ...
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
modelPositionDot1Predicted = [];
modelPositionDot2Predicted = [];
if usePredictedModels
    modelPositionDot1Predicted = dRSA_concatenate(dot1TrialsPredicted);
    modelPositionDot2Predicted = dRSA_concatenate(dot2TrialsPredicted);
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
if usePredictedModels
    model = {modelPositionDot1, modelPositionDot2, ...
        modelPositionDot1Predicted, modelPositionDot2Predicted};
    modelNames = {'position dot1', 'position dot2', ...
        'predicted position dot1', 'predicted position dot2'};
else
    model = {modelPositionDot1, modelPositionDot2};
    modelNames = {'position dot1', 'position dot2'};
end
params.modelNames = modelNames; % human-readable labels for diagnostics/plots.
% In case we do the PCR, it is better to calculate the border outside, because otherwise we would need to recalculate it
% for each iteration
params.nIter = size(subsamples, 3);  % how many iterations?
params.fs = plotSampleRateHz; % framerate or how many samples fit into 1 second
% Use an integer sample window for averaging to avoid non-integer sizes.
avgHalfWindowSamples = floor((trialLen - 1) / 2);
params.AverageTime = avgHalfWindowSamples / params.fs; %in s
if usePredictedModels
    params.modelToTest = [1 2 3 4];  %array of models to test
else
    params.modelToTest = [1 2];  %array of models to test
end
params.Var = 0.1; % how much variance? 
if usePredictedModels
    params.modelDistMeasure = {'euclidean', 'euclidean', 'euclidean', 'euclidean'};
else
    params.modelDistMeasure = {'euclidean', 'euclidean'};
end
params.neuralDistMeasure = 'euclidean'; % pdist expects a string/char or function handle.

%For the PCR
params.dRSAtype = 'corr';% 'PCR';% 'corr'; %
if usePredictedModels
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
    % dRSA uses corr(modelRDM, neuralRDM): x = neural/data time, y = model time.
    xlabel(sprintf('Time in %s (samples)', rowName));
    ylabel(sprintf('Time in %s (samples)', modelNames{iModel}));
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
if saveOutputs
    outputDir = fullfile(simDir, 'output');
    subjectOutputDir = fullfile(outputDir, subjectLabel);
    conditionOutputDir = fullfile(subjectOutputDir, conditionLabel);
    dRSAtypeLabel = lower(paramsCore.dRSAtype);
    analysisOutputDir = fullfile(conditionOutputDir, dRSAtypeLabel);
    resultsOutputDir = fullfile(analysisOutputDir, 'results');
    diagonalOutputDir = fullfile(analysisOutputDir, 'diagonal');
    matricesOutputDir = fullfile(analysisOutputDir, 'matrices');
    pathsOutputDir = fullfile(subjectOutputDir, 'paths');
    local_ensure_dir(resultsOutputDir);
    local_ensure_dir(diagonalOutputDir);
    local_ensure_dir(matricesOutputDir);
    local_ensure_dir(pathsOutputDir);

    optionalTags = {};
    if shuffleDot2Trials
        optionalTags{end + 1} = 'shufDot2';
    end
    if resampleSubsamples
        resampleMode = 'sub';
        if resampleWithReplacement
            resampleMode = 'boot';
        end
        optionalTags{end + 1} = sprintf('rs%d_i%d_%s', ...
            resampleSubsampleCount, resampleIterations, resampleMode);
    end
    runTagParts = [{subjectLabel, conditionLabel, dRSAtypeLabel}, optionalTags];
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
    repro.dot1TrialsPredicted = dot1TrialsPredicted;
    repro.dot2TrialsPredicted = dot2TrialsPredicted;
    repro.dataPosition = dataPosition;
    repro.models = model;
    repro.autocorrBorder = Autocorrborder;
    repro.predictedSimulationFile = predictedSimulationFile;

    save(fullfile(resultsOutputDir, [resultsBase '_results.mat']), ...
        'dRSA_position', ...
        'dRSA_diagonal_position', ...
        'params', 'paramsDiagonal', 'repro');

    print(figMatrices, fullfile(matricesOutputDir, [resultsBase '_matrices.png']), '-dpng', '-r300');
    print(figLag, fullfile(diagonalOutputDir, [resultsBase '_diagonal.png']), '-dpng', '-r300');
    if iCond == 1
        pathsBase = strjoin({'paths', subjectLabel, 'allConditions'}, '_');
        pathsBase = regexprep(pathsBase, '\s+', '');
        print(figAllPaths, fullfile(pathsOutputDir, [pathsBase '_all_conditions.png']), ...
            '-dpng', '-r300');
    end
end

end

% matrices
% lag plots
%

%% Local helpers
function local_ensure_dir(dirPath)
% LOCAL_ENSURE_DIR Create the directory if it is missing (including parents).
if isempty(dirPath)
    return;
end
if ~exist(dirPath, 'dir')
    mkdir(dirPath);
end
end
