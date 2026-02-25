% pipeline_recursive.m
%
% Purpose:
%   Select the "best" trial subset by scanning random batches and scoring
%   the max cross-correlation between dot1 position data and dot2 position
%   model (via dRSA_coreFunction). After N iterations, the subset with the
%   minimum max correlation is re-used to run a full dRSA analysis (position
%   + direction models) once, with plots for inspection.
%
% Example usage (from simulations/scripts in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   pipeline_recursive;
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   pipeline_recursive;
%
% Configuration (edit at the start of the script):
%   - participantNumber selects MovDot_SubXX.mat and the derived simulation
%     inputs (MovDot_SubXX_nondeviant/deviant.mat).
%   - inputCondition selects which condition to load for the simulation run.
%   - nRecursiveRuns controls how many subsets to score (N times).
%   - subsetTrialCount sets the number of trials per subset (unique within).
%   - shuffleDot2Trials decouples dot2 trial order from dot1 (diagnostic).
%   - rngSeed controls reproducibility of subset sampling.
%   - kBest controls how many top subsets (lowest scores) are kept in memory.
%   - nPlotBest controls how many top-scoring subsets are plotted at the end.
%   - saveEvery controls how often to checkpoint best subsets to disk.
%   - scoreOption selects the scoring metric used to rank subsets.
%   - bestSubsampleFile stores top-k subset indices for resuming scans.
%   - suppressDispText silences function-level status text (0 = show, 1 = suppress).
%   - plotSwitch controls figure visibility (1 = show, 0 = hidden but saved).
%
% Inputs:
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat with:
%       dot1GreenPaths (trials × 2 × time) and dot2YellowPaths (same shape)
%
% Outputs:
%   - bestResults: dRSA outputs for the selected subset.
%   - Figures: only for the selected subset (paths + dRSA matrices).
%   - best_subsample_<participantNumber>_<inputCondition>_<scoreOption>.mat:
%     checkpoint with top-k scores (mean/max) + subset indices for resuming scans.
%   - Duplicate subsets (order-insensitive) are skipped when updating top-k.
%
% Key assumptions:
%   - dot1GreenPaths and dot2YellowPaths are in visual degrees.
%   - Trials have uniform frame counts within each condition.
%   - The functions in simulations/functions are on the MATLAB path.
%   - Paths are resolved relative to this script, not the MATLAB cwd.
%   - Best subsets are the lowest scoreOption values (smaller is better).
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
addpath(scriptDir);
addpath(simDir);
addpath(fullfile(simDir, 'functions'));
addpath(fullfile(simDir, 'debug'));
addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');

%% Participant configuration
% Data flow: participant number + condition -> input filenames -> simulation data.
participantNumber = 89;
inputCondition = 'nondeviant'; % 'nondeviant' or 'deviant' to match output files.
inputCondition = lower(inputCondition);
shuffleDot2Trials = false; % true to shuffle dot2 trials relative to dot1 before dRSA

% Recursive dRSA configuration (N runs with unique trials per subset).
nRecursiveRuns = 20000; % N times to run dRSA on different trial subsets.
subsetTrialCount = 20; % number of trials per run (unique within each subset).

% Best-subset tracking + checkpointing.
kBest = 1000; % number of best scores to keep in memory.
nPlotBest = 20; % number of top-scoring subsets to plot at the end.
saveEvery = 1000; % save checkpoint every N iterations.
statusEvery = 250; % print status every N iterations.

% Scoring configuration for subset ranking.
% "mean" -> mean(abs(values(:))); "max" -> max(abs(values(:))).
scoreOption = 'max';

% Random seed for subset sampling (set numeric for reproducibility).
% Use 'shuffle' to re-seed from the current time.
rngSeed = 'shuffle';

% Control console output from helper functions (0 = show, 1 = suppress).
suppressDispText = 1;
plotSwitch = 0; % show figures when 1; when 0, create hidden figures for saving.

validConditions = {'nondeviant', 'deviant'};
if ~ismember(inputCondition, validConditions)
    error('inputCondition must be one of: %s', strjoin(validConditions, ', '));
end
validScoreOptions = {'mean', 'max'};
if ~ismember(scoreOption, validScoreOptions)
    error('scoreOption must be one of: %s', strjoin(validScoreOptions, ', '));
end
if ~isscalar(nPlotBest) || nPlotBest < 1 || nPlotBest ~= floor(nPlotBest)
    error('nPlotBest must be a positive integer.');
end
simulationInputDir = fullfile(simDir, 'input');
outputDir = fullfile(simDir, 'output', sprintf('sub%d', participantNumber));
bestSubsampleFile = fullfile(outputDir, ...
    sprintf('best_subsample_%d_%s_%s.mat', participantNumber, inputCondition, scoreOption));
reportPattern = 'settings_report_*.md';

%% Create input file needed from paths
% Data flow: MovDot_SubXX.mat -> split into condition-specific dot paths.
inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
build_movdot_simulation_inputs(inputFile, ...
    'OutputDir', simulationInputDir, ...
    'MoveDotScript', fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m'), ...
    'suppressDispText', suppressDispText);

%% Load condition-specific dot paths (visual degrees)
% Data flow: condition .mat -> dot1/dot2 arrays -> placeholders for simulation.
simulationInputFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_%s.mat', participantNumber, inputCondition));
simulationData = load(simulationInputFile);
dot1GreenPaths = simulationData.dot1GreenPaths;
dot2YellowPaths = simulationData.dot2YellowPaths;

% Prepare per-dot inputs (trial × features × time).
dot1TrialsAll = dot1GreenPaths;
dot2TrialsAll = dot2YellowPaths;

%% Optional trial shuffling (dot2 relative to dot1)
% Data flow: dot2Trials -> permuted trial order -> decoupled dot1/dot2 pairing.
% Purpose: test whether dot1-dot2 structure is driven by trial-locked coupling.
if shuffleDot2Trials
    rng('shuffle'); % avoid reusing the subject-seeded RNG for the shuffle
    shuffleOrder = randperm(size(dot2TrialsAll, 1));
    dot2TrialsAll = dot2TrialsAll(shuffleOrder, :, :);
end

%% Validate trial counts and seed RNG for subset sampling
% Data flow: total trials -> sampling constraints -> rng seeded.
nTrialsTotal = size(dot1TrialsAll, 1);
if subsetTrialCount > nTrialsTotal
    error('subsetTrialCount (%d) exceeds available trials (%d).', ...
        subsetTrialCount, nTrialsTotal);
end
if ischar(rngSeed) || isstring(rngSeed)
    rng(char(rngSeed));
else
    rng(rngSeed);
end

%% Resume or initialize best-subset tracking
% Data flow: checkpoint file -> best score values + subset indices in memory.
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

bestScores = inf(kBest, 1);
bestScoresMean = inf(kBest, 1);
bestScoresMax = inf(kBest, 1);
bestSubsetIdx = nan(kBest, subsetTrialCount);
runsCompleted = 0;

if exist(bestSubsampleFile, 'file')
    resumeData = load(bestSubsampleFile);
    if isfield(resumeData, 'participantNumber') && ...
            resumeData.participantNumber ~= participantNumber
        error('Checkpoint participantNumber (%d) does not match current (%d).', ...
            resumeData.participantNumber, participantNumber);
    end
    if isfield(resumeData, 'inputCondition') && ...
            ~strcmp(resumeData.inputCondition, inputCondition)
        error('Checkpoint inputCondition (%s) does not match current (%s).', ...
            resumeData.inputCondition, inputCondition);
    end
    if isfield(resumeData, 'subsetTrialCount') && ...
            resumeData.subsetTrialCount ~= subsetTrialCount
        error('Checkpoint subsetTrialCount (%d) does not match current (%d).', ...
            resumeData.subsetTrialCount, subsetTrialCount);
    end
    if isfield(resumeData, 'scoreOption') && ...
            ~strcmp(resumeData.scoreOption, scoreOption)
        error('Checkpoint scoreOption (%s) does not match current (%s).', ...
            resumeData.scoreOption, scoreOption);
    end

    if isfield(resumeData, 'bestScores')
        loadedScores = resumeData.bestScores(:);
    else
        loadedScores = [];
    end
    if isfield(resumeData, 'bestScoresMean')
        loadedScoresMean = resumeData.bestScoresMean(:);
    else
        loadedScoresMean = [];
    end
    if isfield(resumeData, 'bestScoresMax')
        loadedScoresMax = resumeData.bestScoresMax(:);
    else
        loadedScoresMax = [];
    end
    if isfield(resumeData, 'bestSubsetIdx')
        loadedSubsetIdx = resumeData.bestSubsetIdx;
    else
        loadedSubsetIdx = [];
    end

    if isempty(loadedScores) && ~isempty(loadedScoresMean) && ~isempty(loadedScoresMax)
        switch scoreOption
            case 'mean'
                loadedScores = loadedScoresMean;
            case 'max'
                loadedScores = loadedScoresMax;
        end
    end

    if ~isempty(loadedScores) && ~isempty(loadedSubsetIdx)
        if size(loadedSubsetIdx, 2) ~= subsetTrialCount
            error('Checkpoint subset size (%d) does not match current (%d).', ...
                size(loadedSubsetIdx, 2), subsetTrialCount);
        end
        if numel(loadedScores) ~= size(loadedSubsetIdx, 1)
            error('Checkpoint score count does not match subset row count.');
        end
        if ~isempty(loadedScoresMean) && numel(loadedScoresMean) ~= size(loadedSubsetIdx, 1)
            error('Checkpoint mean score count does not match subset row count.');
        end
        if ~isempty(loadedScoresMax) && numel(loadedScoresMax) ~= size(loadedSubsetIdx, 1)
            error('Checkpoint max score count does not match subset row count.');
        end

        if numel(loadedScores) > kBest
            [~, sortIdx] = sort(loadedScores, 'ascend');
            sortIdx = sortIdx(1:kBest);
            loadedScores = loadedScores(sortIdx);
            if ~isempty(loadedScoresMean)
                loadedScoresMean = loadedScoresMean(sortIdx);
            end
            if ~isempty(loadedScoresMax)
                loadedScoresMax = loadedScoresMax(sortIdx);
            end
            loadedSubsetIdx = loadedSubsetIdx(sortIdx, :);
        end
        bestScores(1:numel(loadedScores)) = loadedScores;
        if ~isempty(loadedScoresMean)
            bestScoresMean(1:numel(loadedScoresMean)) = loadedScoresMean;
        end
        if ~isempty(loadedScoresMax)
            bestScoresMax(1:numel(loadedScoresMax)) = loadedScoresMax;
        end
        bestSubsetIdx(1:size(loadedSubsetIdx, 1), :) = loadedSubsetIdx;
    end

    if isfield(resumeData, 'runsCompleted')
        runsCompleted = resumeData.runsCompleted;
    end
end

%% Scan random subsets and score dot1 vs dot2 position cross-correlation
% Data flow: subset indices -> position-only dRSA -> max correlation score.
startRun = runsCompleted + 1;
skipScan = false;
if startRun > nRecursiveRuns
    warning('Checkpoint already has %d iterations; no new runs scheduled.', runsCompleted);
    % No new runs: skip the scan and proceed to plot the best saved subset.
    skipScan = true;
end

if ~skipScan
    for runIdx = startRun:nRecursiveRuns
        % Select a random subset of trials (unique within the subset).
        subsetIdx = randperm(nTrialsTotal, subsetTrialCount);
        dot1Trials = dot1TrialsAll(subsetIdx, :, :);
        dot2Trials = dot2TrialsAll(subsetIdx, :, :);

        %% Prepare position data/model for scoring
        % Data flow: trial arrays -> concatenated dot1 data + dot2 model.
        dataPosition = dRSA_concatenate(dot1Trials, [], 0, ...
            'suppressDispText', suppressDispText); % position data from dot 1
        modelPositionDot2 = dRSA_concatenate(dot2Trials, [], 0, ...
            'suppressDispText', suppressDispText); % position model from dot 2

        %% Build trigger mask and subsamples (trial-locked)
        % Data flow: trial length -> trigger mask -> subsample windows.
        trialLen = size(dot1Trials, 3);
        totalTime = size(dataPosition, 2);
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
        subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt, ...
            'suppressDispText', suppressDispText);

        %% Configure position-only dRSA parameters for scoring
        % Data flow: position data + dot2 model -> single-model dRSA matrix.
        scoreParams.nIter = 1; % force single-iteration dRSA per subset
        scoreParams.fs = 120; % framerate (samples per second)
        avgHalfWindowSamples = floor((trialLen - 1) / 2);
        scoreParams.AverageTime = avgHalfWindowSamples / scoreParams.fs; % seconds
        scoreParams.modelToTest = 1; % only dot2 position model
        scoreParams.Var = 0.1; % variance parameter
        scoreParams.modelDistMeasure = {'euclidean'};
        scoreParams.neuralDistMeasure = scoreParams.modelDistMeasure{1};
        scoreParams.dRSAtype = 'corr';

        % Use dot2 as the only model for position cross-correlation scoring.
        scoreModel = {modelPositionDot2};
        dRSA_score = dRSA_coreFunction(dataPosition, scoreModel, scoreParams, ...
            'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', [], ...
            'suppressDispText', suppressDispText);

        %% Store score and subset index (top-k only)
        % Data flow: dRSA score matrix -> mean/max scores -> best arrays (kBest).
        scoreValues = abs(dRSA_score(:));
        scoreMean = mean(scoreValues);
        scoreMax = max(scoreValues);
        switch scoreOption
            case 'mean'
                score = scoreMean;
            case 'max'
                score = scoreMax;
        end

        %% Reject duplicate subsets (order-insensitive) before updating top-k
        % Data flow: candidate subset -> sorted index -> membership against stored sets.
        finiteMask = isfinite(bestScores);
        isDuplicate = false;
        if any(finiteMask)
            sortedSubset = sort(subsetIdx);
            sortedStored = sort(bestSubsetIdx(finiteMask, :), 2);
            % Order-insensitive equality: subset is duplicate if the sorted row matches.
            isDuplicate = ismember(sortedSubset, sortedStored, 'rows');
        end

        if ~isDuplicate
            insertIdx = find(isinf(bestScores), 1, 'first');
            if ~isempty(insertIdx)
                bestScores(insertIdx) = score;
                bestScoresMean(insertIdx) = scoreMean;
                bestScoresMax(insertIdx) = scoreMax;
                bestSubsetIdx(insertIdx, :) = subsetIdx;
            else
                [worstScore, worstIdx] = max(bestScores);
                if score < worstScore
                    bestScores(worstIdx) = score;
                    bestScoresMean(worstIdx) = scoreMean;
                    bestScoresMax(worstIdx) = scoreMax;
                    bestSubsetIdx(worstIdx, :) = subsetIdx;
                end
            end
        end
        runsCompleted = runIdx;

        %% Checkpoint best subsets to disk
        % Data flow: best arrays + metadata -> .mat file for resuming runs.
        if mod(runIdx, saveEvery) == 0 || runIdx == nRecursiveRuns
            save(bestSubsampleFile, 'bestScores', 'bestScoresMean', ...
                'bestScoresMax', 'bestSubsetIdx', ...
                'runsCompleted', 'participantNumber', 'inputCondition', ...
                'subsetTrialCount', 'nRecursiveRuns', 'kBest', 'saveEvery', ...
                'rngSeed', 'scoreOption');
        end

        %% Report progress for long runs
        % Data flow: iteration counters -> console status line.
        if mod(runIdx, statusEvery) == 0 || runIdx == startRun
            fprintf('Status: run %d/%d (scoreOption=%s).\n', ...
                runIdx, nRecursiveRuns, scoreOption);
        end
    end
end

%% Select top-N subsets with minimum scores
% Data flow: best score values -> top-N subset indices -> data slices for analysis.
finiteMask = isfinite(bestScores);
if ~any(finiteMask)
    error('No valid scores available to select a best subset.');
end
finiteIdx = find(finiteMask);
availableScores = bestScores(finiteMask);
[~, sortIdx] = sort(availableScores, 'ascend');
nPlot = min(nPlotBest, numel(sortIdx));
topIdx = sortIdx(1:nPlot);
bestScoresSelected = availableScores(topIdx);
bestSubsetIdxSelected = bestSubsetIdx(finiteIdx(topIdx), :);

% Extract trials for each selected subset.
dot1TrialsSelected = cell(nPlot, 1);
dot2TrialsSelected = cell(nPlot, 1);
for iPlot = 1:nPlot
    dot1TrialsSelected{iPlot} = dot1TrialsAll(bestSubsetIdxSelected(iPlot, :), :, :);
    dot2TrialsSelected{iPlot} = dot2TrialsAll(bestSubsetIdxSelected(iPlot, :), :, :);
end

%% Prepare data, models, and trial-locked subsamples for final analysis
% Data flow: selected trials -> per-subset models -> dRSA outputs per subset.
dRSA_position = cell(nPlot, 1);
dRSA_direction = cell(nPlot, 1);
paramsPerSubset = cell(nPlot, 1);
trialLenPerSubset = nan(nPlot, 1);
pathPlotFiles = cell(nPlot, 1);
matrixPlotFiles = cell(nPlot, 1);
avgMatrixPlotFile = '';

%% Configure figure visibility for plotting
% Data flow: plotSwitch -> figure visibility flag for all plots in this script.
if plotSwitch == 0
    figureVisibility = 'off';
else
    figureVisibility = 'on';
end

for iPlot = 1:nPlot
    %% Prepare data, models, and trial-locked subsamples
    % Data flow: trial arrays -> position + direction models -> concatenated time series.
    % Position models use raw x/y paths; direction models use per-timepoint angles.
    dot1Trials = dot1TrialsSelected{iPlot};
    dot2Trials = dot2TrialsSelected{iPlot};

    dataPosition = dRSA_concatenate(dot1Trials, [], 0, ...
        'suppressDispText', suppressDispText); % position data from dot 1

    % Position models (per dot).
    modelPositionDot1 = dataPosition;
    modelPositionDot2 = dRSA_concatenate(dot2Trials, [], 0, ...
        'suppressDispText', suppressDispText);

    % Direction models (per dot).
    % Data flow: dot positions -> frame-to-frame displacement -> angle -> unit vectors.
    dot1Dx = diff(dot1Trials(:, 1, :), 1, 3); % frame-to-frame x displacement (trials x 1 x time-1)
    dot1Dy = diff(dot1Trials(:, 2, :), 1, 3); % frame-to-frame y displacement (trials x 1 x time-1)
    dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx); % pad so length matches original time
    dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy); % pad so length matches original time
    dot1Angle = atan2(dot1Dy, dot1Dx); % per-timepoint direction angle in radians
    dot1Direction = cat(2, cos(dot1Angle), sin(dot1Angle)); % unit vectors [cos; sin]
    modelDirectionDot1 = dRSA_concatenate(dot1Direction, [], 0, ...
        'suppressDispText', suppressDispText);

    dot2Dx = diff(dot2Trials(:, 1, :), 1, 3); % frame-to-frame x displacement for dot 2
    dot2Dy = diff(dot2Trials(:, 2, :), 1, 3); % frame-to-frame y displacement for dot 2
    dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx); % pad so length matches original time
    dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy); % pad so length matches original time
    dot2Angle = atan2(dot2Dy, dot2Dx); % per-timepoint direction angle in radians
    dot2Direction = cat(2, cos(dot2Angle), sin(dot2Angle)); % unit vectors [cos; sin]
    modelDirectionDot2 = dRSA_concatenate(dot2Direction, [], 0, ...
        'suppressDispText', suppressDispText);

    dataDirection = modelDirectionDot1; % use dot 1 direction as observed data

    %% Build trigger mask and subsamples (trial-locked)
    % Data flow: trial length -> trigger mask -> subsample windows.
    trialLen = size(dot1Trials, 3);
    totalTime = size(dataPosition, 2);
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
    subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt, ...
        'suppressDispText', suppressDispText);

    %% Configure models and core dRSA parameters for final analysis
    % Data flow: models + subsamples -> dRSA params -> autocorrelation border.
    model = {modelPositionDot1, modelPositionDot2, modelDirectionDot1, modelDirectionDot2};
    modelNames = {'position dot1', 'position dot2', 'direction dot1', 'direction dot2'};

    params.nIter = 1; % force single-iteration dRSA for the selected subset
    params.fs = 120; % framerate (samples per second)
    avgHalfWindowSamples = floor((trialLen - 1) / 2);
    params.AverageTime = avgHalfWindowSamples / params.fs; % seconds
    params.modelToTest = [1 2 3 4]; % indices of models to test
    params.Var = 0.1; % variance parameter
    params.modelDistMeasure = {'euclidean', 'euclidean', 'cosine', 'cosine'};
    % Set neural distance to match the data-representing model (position vs direction).
    % Set later

    % PCR configuration.
    params.dRSAtype = 'corr';
    params.modeltoRegressout = {[2 3 4] [1 3 4] [1 2 4] [1 2 3]};
    params.PCR.AdditionalPCA = 1;
    params.PCR.RegressAutocor = 1;
    params.PCR.RessModel = 1;

    % Autocorrelation border is only needed for non-corr dRSA.
    Autocorrborder = []; % empty signals "no border" to downstream code paths.
    if ~strcmp(params.dRSAtype, 'corr') % Autocorrelation not necessary for 'corr' type.
        Autocorrborder = dRSA_border(model, subsamples, params, ...
            'suppressDispText', suppressDispText);
    end

    %% Compute shared axis limits for dot path plots
    % Data flow: dot1/dot2 trials -> pooled x/y extents -> axis limits for both plots.
    dot1X = dot1Trials(:, 1, :);
    dot1Y = dot1Trials(:, 2, :);
    dot2X = dot2Trials(:, 1, :);
    dot2Y = dot2Trials(:, 2, :);
    xMin = min([dot1X(:); dot2X(:)]);
    xMax = max([dot1X(:); dot2X(:)]);
    yMin = min([dot1Y(:); dot2Y(:)]);
    yMax = max([dot1Y(:); dot2Y(:)]);
    xPad = 0.05 * max(xMax - xMin, eps);
    yPad = 0.05 * max(yMax - yMin, eps);
    sharedAxisLimits = [xMin - xPad, xMax + xPad, yMin - yPad, yMax + yPad];

    %% Plot subset paths for the selected subset
    % Data flow: selected trials -> path plots for dot1 and dot2.
    figPaths = figure('Name', sprintf('Top %d/%d subset: dot paths', iPlot, nPlot), ...
        'NumberTitle', 'off', 'Visible', figureVisibility);
    tiledlayout(figPaths, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    plot_paths_no_save(dot1Trials, 'Title', sprintf('Top %d/%d: Dot 1 paths', iPlot, nPlot), ...
        'AxisLimits', sharedAxisLimits, 'suppressDispText', suppressDispText);
    nexttile;
    plot_paths_no_save(dot2Trials, 'Title', sprintf('Top %d/%d: Dot 2 paths', iPlot, nPlot), ...
        'AxisLimits', sharedAxisLimits, 'suppressDispText', suppressDispText);
    % Save path plot to output directory for reporting.
    pathPlotFiles{iPlot} = fullfile(outputDir, sprintf('paths_top_%02d.png', iPlot));
    exportgraphics(figPaths, pathPlotFiles{iPlot}, 'Resolution', 300);

    %% Run dRSA for position data (single iteration)
    % Data flow: position data + models -> dRSA matrices (time x time x model).
    % Match neural distance to the position model (model 1) when data are positions.
    params.neuralDistMeasure = params.modelDistMeasure{1};
    dRSA_position{iPlot} = dRSA_coreFunction(dataPosition, model, params, ...
        'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', Autocorrborder, ...
        'suppressDispText', suppressDispText);

    %% Run dRSA for direction data (single iteration)
    % Data flow: direction data + models -> dRSA matrices (time x time x model).
    % Match neural distance to the direction model (model 3) when data are directions.
    params.neuralDistMeasure = params.modelDistMeasure{3};
    dRSA_direction{iPlot} = dRSA_coreFunction(dataDirection, model, params, ...
        'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', Autocorrborder, ...
        'suppressDispText', suppressDispText);

    %% Plot dRSA matrices for the selected subset
    % Data flow: dRSA per data type -> grid of heatmaps for inspection.
    nModels = size(dRSA_position{iPlot}, 3);
    nTime = size(dRSA_position{iPlot}, 1);
    dRSAAll = {dRSA_position{iPlot}, dRSA_direction{iPlot}};
    rowNames = {'Data: position dot1', 'Data: direction dot1'};

    figMatrices = figure('Name', sprintf('Top %d/%d subset: dRSA matrices', iPlot, nPlot), ...
        'NumberTitle', 'off', 'Visible', figureVisibility);
    tiledlayout(figMatrices, 2, nModels, 'Padding', 'compact', 'TileSpacing', 'compact');
    rowAxes = gobjects(numel(dRSAAll), 1);
    for iRow = 1:numel(dRSAAll)
        for iModel = 1:nModels
            nexttile((iRow - 1) * nModels + iModel);
            imagesc(dRSAAll{iRow}(:, :, iModel));
            set(gca, 'YDir', 'normal');
            axis image;
            colorbar;
            title(sprintf('dRSA %s', modelNames{iModel}));
            % dRSA uses corr(modelRDM, neuralRDM): x = neural/data time, y = model time.
            xlabel(sprintf('Time in %s (samples)', rowNames{iRow}));
            ylabel(sprintf('Time in %s (samples)', modelNames{iModel}));
            hold on;
            plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
            hold off;
            if iModel == 1
                rowAxes(iRow) = gca;
            end
        end
    end
    % Save dRSA matrices plot to output directory for reporting.
    matrixPlotFiles{iPlot} = fullfile(outputDir, sprintf('dRSA_matrices_top_%02d.png', iPlot));
    exportgraphics(figMatrices, matrixPlotFiles{iPlot}, 'Resolution', 300);
    for iRow = 1:numel(rowAxes)
        if isgraphics(rowAxes(iRow))
            text(rowAxes(iRow), -0.35, 0.5, rowNames{iRow}, ...
                'Units', 'normalized', 'Rotation', 90, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end

    %% Capture per-subset metadata
    % Data flow: subset metadata + dRSA outputs -> per-subset containers.
    paramsPerSubset{iPlot} = params;
    trialLenPerSubset(iPlot) = trialLen;
end

%% Compute and plot average dRSA matrices across selected subsets
% Data flow: per-subset dRSA -> mean across subsets -> averaged plot + save.
nModels = size(dRSA_position{1}, 3);
nTime = size(dRSA_position{1}, 1);
dRSA_position_all = cat(4, dRSA_position{:});
dRSA_direction_all = cat(4, dRSA_direction{:});
dRSA_position_avg = mean(dRSA_position_all, 4, 'omitnan');
dRSA_direction_avg = mean(dRSA_direction_all, 4, 'omitnan');
dRSAAllAvg = {dRSA_position_avg, dRSA_direction_avg};
rowNamesAvg = {'Data: position dot1 (avg)', 'Data: direction dot1 (avg)'};

figMatricesAvg = figure('Name', 'Average dRSA matrices (top subsets)', ...
    'NumberTitle', 'off', 'Visible', figureVisibility);
tiledlayout(figMatricesAvg, 2, nModels, 'Padding', 'compact', 'TileSpacing', 'compact');
rowAxesAvg = gobjects(numel(dRSAAllAvg), 1);
for iRow = 1:numel(dRSAAllAvg)
    for iModel = 1:nModels
        nexttile((iRow - 1) * nModels + iModel);
        imagesc(dRSAAllAvg{iRow}(:, :, iModel));
        set(gca, 'YDir', 'normal');
        axis image;
        colorbar;
        title(sprintf('dRSA %s', modelNames{iModel}));
        % dRSA uses corr(modelRDM, neuralRDM): x = neural/data time, y = model time.
        xlabel(sprintf('Time in %s (samples)', rowNamesAvg{iRow}));
        ylabel(sprintf('Time in %s (samples)', modelNames{iModel}));
        hold on;
        plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
        hold off;
        if iModel == 1
            rowAxesAvg(iRow) = gca;
        end
    end
end
for iRow = 1:numel(rowAxesAvg)
    if isgraphics(rowAxesAvg(iRow))
        text(rowAxesAvg(iRow), -0.35, 0.5, rowNamesAvg{iRow}, ...
            'Units', 'normalized', 'Rotation', 90, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
avgMatrixPlotFile = fullfile(outputDir, 'dRSA_matrices_avg.png');
exportgraphics(figMatricesAvg, avgMatrixPlotFile, 'Resolution', 300);

%% Store final outputs in the workspace
% Data flow: best subset metadata + dRSA outputs -> bestResults struct.
bestResults.subsetIdx = bestSubsetIdxSelected;
bestResults.bestScore = bestScoresSelected;
bestResults.scoreOption = scoreOption;
bestResults.bestScoresMean = bestScoresMean;
bestResults.bestScoresMax = bestScoresMax;
bestResults.dRSA_position = dRSA_position;
bestResults.dRSA_direction = dRSA_direction;
bestResults.modelNames = modelNames;
bestResults.paramsCore = paramsPerSubset;
bestResults.trialLen = trialLenPerSubset;
bestResults.nPlotBest = nPlot;
bestResults.bestSubsetIdx = bestSubsetIdxSelected(1, :);
bestResults.bestScoreSingle = bestScoresSelected(1);

%% Update settings report with plots and subset summary
% Data flow: selected subset metadata + saved plot filenames -> Markdown report append.
reportCandidates = dir(fullfile(outputDir, reportPattern));
if isempty(reportCandidates)
    reportFile = fullfile(outputDir, ...
        sprintf('settings_report_%s.md', datestr(now, 'yyyymmddTHHMMSS')));
    warning('No existing settings_report_*.md found; creating %s.', reportFile);
else
    [~, newestIdx] = max([reportCandidates.datenum]);
    reportFile = fullfile(outputDir, reportCandidates(newestIdx).name);
end

reportLines = {};
reportLines{end+1} = '';
reportLines{end+1} = '## Recursive pipeline plots and subset summary';
reportLines{end+1} = '';
reportLines{end+1} = sprintf('- Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
reportLines{end+1} = sprintf('- Selected subsets: %d', nPlot);
reportLines{end+1} = sprintf('- Score option: %s', scoreOption);
reportLines{end+1} = sprintf('- Score min (selected): %.6f', min(bestScoresSelected));
reportLines{end+1} = sprintf('- Score max (selected): %.6f', max(bestScoresSelected));
reportLines{end+1} = '';
reportLines{end+1} = '### Plot gallery';
reportLines{end+1} = '';
for iPlot = 1:nPlot
    if ~isempty(pathPlotFiles{iPlot})
        reportLines{end+1} = sprintf('##### Subsample %d - %s', iPlot, inputCondition);
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('![Dot paths top %d](%s)', ...
            iPlot, pathPlotFiles{iPlot});
        reportLines{end+1} = '';
    end
    if ~isempty(matrixPlotFiles{iPlot})
        reportLines{end+1} = sprintf('![dRSA matrices top %d](%s)', ...
            iPlot, matrixPlotFiles{iPlot});
        reportLines{end+1} = '';
    end
end
if ~isempty(avgMatrixPlotFile)
    reportLines{end+1} = sprintf('##### Average dRSA - %s', inputCondition);
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('![dRSA matrices average](%s)', avgMatrixPlotFile);
    reportLines{end+1} = '';
end
reportLines{end+1} = '### Selected subset indices and scores';
reportLines{end+1} = '';
reportLines{end+1} = '| Rank | Score (selected) | Score mean | Score max | Trial indices |';
reportLines{end+1} = '| --- | --- | --- | --- | --- |';
for iPlot = 1:nPlot
    subsetIdx = bestSubsetIdxSelected(iPlot, :);
    subsetText = sprintf('%d ', subsetIdx);
    subsetText = strtrim(subsetText);
    reportLines{end+1} = sprintf('| %d | %.6f | %.6f | %.6f | `%s` |', ...
        iPlot, bestScoresSelected(iPlot), bestScoresMean(finiteIdx(topIdx(iPlot))), ...
        bestScoresMax(finiteIdx(topIdx(iPlot))), subsetText);
end
reportLines{end+1} = '';

if exist(reportFile, 'file')
    fid = fopen(reportFile, 'a');
    if fid == -1
        error('Could not open report file for appending: %s', reportFile);
    end
    fprintf(fid, '%s\n', reportLines{:});
    fclose(fid);
else
    fid = fopen(reportFile, 'w');
    if fid == -1
        error('Could not create report file: %s', reportFile);
    end
    fprintf(fid, '# DAD Settings Report\n\n');
    fprintf(fid, '%s\n', reportLines{:});
    fclose(fid);
end

%% Local helper: plot paths without saving to disk
% Data flow: paths -> scatter plot with time-graded colors (no file output).
function figHandle = plot_paths_no_save(paths, varargin)
    % plot_paths_no_save
    %
    % Purpose:
    %   Replicate the visualization style of simulations/debug/plot_paths.m
    %   without writing a PNG to disk.
    %
    % Example usage:
    %   fig = plot_paths_no_save(dot1Trials, 'Title', 'Dot 1 paths');
    %
    % Inputs:
    %   paths : numeric array, trials × 2 × time (features are [x, y])
    %
    % Name/value options:
    %   'Title'      : figure title string (default: '')
    %   'Colormap'   : colormap name or n×3 array (default: parula)
    %   'MarkerSize' : scatter marker size (default: 8)
    %   'UseAlpha'   : true/false to enable transparency (default: false)
    %   'AlphaValue' : alpha for marker face/edge when UseAlpha = true (default: 0.1)
    %   'AxisEqual'  : true/false to enforce equal axes (default: true)
    %   'AxisLimits' : 1×4 [xmin xmax ymin ymax] (default: [])
    %   'suppressDispText' : suppress console output (default: false)
    %
    % Outputs:
    %   figHandle : handle to the created figure

    %% Input validation and options
    parser = inputParser;
    parser.addRequired('paths', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addParameter('Title', '', @(x) ischar(x) || isstring(x));
    parser.addParameter('Colormap', parula, @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && size(x, 2) == 3));
    parser.addParameter('MarkerSize', 8, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('UseAlpha', false, @(x) islogical(x) && isscalar(x));
    parser.addParameter('AlphaValue', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    parser.addParameter('AxisEqual', true, @(x) islogical(x) && isscalar(x));
    parser.addParameter('AxisLimits', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 4));
    parser.addParameter('suppressDispText', false, ...
        @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
    parser.parse(paths, varargin{:});

    opts = parser.Results;
    % Note: suppressDispText is accepted for compatibility; this helper is quiet.

    %% Flatten trials into a single scatter set
    nTrials = size(paths, 1);
    nTime = size(paths, 3);

    xVals = reshape(paths(:, 1, :), [], 1);
    yVals = reshape(paths(:, 2, :), [], 1);
    timeIdx = repmat(1:nTime, nTrials, 1);
    timeIdx = timeIdx(:);

    % Build per-point color map (time -> RGB).
    if ischar(opts.Colormap) || isstring(opts.Colormap)
        cmap = feval(char(opts.Colormap), nTime);
    else
        cmap = opts.Colormap;
        if size(cmap, 1) ~= nTime
            cmap = interp1(linspace(1, nTime, size(cmap, 1)), cmap, 1:nTime);
        end
    end
    pointColors = cmap(timeIdx, :);

    %% Plot paths with time-graded color
    figHandle = gcf;
    ax = gca;
    hold(ax, 'on');

    if opts.UseAlpha
        scatter(ax, xVals, yVals, opts.MarkerSize, pointColors, 'filled', ...
            'MarkerFaceAlpha', opts.AlphaValue, 'MarkerEdgeAlpha', opts.AlphaValue);
    else
        scatter(ax, xVals, yVals, opts.MarkerSize, pointColors, 'filled');
    end

    colormap(ax, cmap);
    caxis(ax, [1 nTime]);
    cb = colorbar(ax);
    cb.Label.String = 'Time (samples)';

    xlabel(ax, 'X (visual degrees)');
    ylabel(ax, 'Y (visual degrees)');
    title(ax, char(opts.Title));

    if opts.AxisEqual
        axis(ax, 'equal');
    end
    if ~isempty(opts.AxisLimits)
        axis(ax, opts.AxisLimits);
    end

    hold(ax, 'off');
end
