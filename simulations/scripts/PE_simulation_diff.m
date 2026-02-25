% PE_simulation_diff.m
%
% Purpose:
%   Run the counterfactual dRSA simulation pipeline and add a deviant-only
%   prediction-error (PE) analysis where PE is defined directly at the signal
%   level (not as a difference between precomputed RDMs).
%
% Canonical PE definition used in this script:
%   PE_signal = observed_signal - predicted_signal
%
% Scope and design:
%   - Keep the same core pipeline shape as counterfactual_simulation.m.
%   - Keep standard nondeviant/deviant counterfactual outputs.
%   - Add PE dRSA outputs for deviant using four PE models:
%       1) PE position dot1
%       2) PE position dot2
%       3) PE direction dot1
%       4) PE direction dot2
%   - Evaluate PE models against three PE neural targets:
%       a) neuralPE        = observed neural signal - predicted neural signal
%       b) neuralPredicted = predicted neural signal
%       c) neuralObserved  = observed neural signal
%   - Keep dRSA mode corr-only in this script version.
%   - Do not add any separate "dRSA from precomputed RDM" branch.
%
% Existing-results decision system:
%   Configure existingResultsAction at the top of this script:
%     [] : interactive prompt (calls input with explicit 3 choices)
%      1 : reuse existing dataset and replot
%      2 : full rerun and overwrite outputs (default for batch safety)
%      3 : just error out if results already exist
%
% PE post-deviance enforcement rule:
%   PE analyses always run post-deviance only.
%   - If cutPostDev is true: PE uses globally cut post-deviance streams.
%   - If cutPostDev is false: PE applies an internal PE-only post-deviance cut
%     to both observed and predicted streams before PE_signal is computed.
%
% Data flow summary:
%   1) Build/load center-relative condition inputs from MovDot_SubXX.mat.
%   2) Build observed position/direction streams for each condition.
%   3) Build predicted streams for deviant condition.
%   4) Run standard counterfactual dRSA.
%   5) Build signal-level PE streams (observed - predicted), enforce PE
%      post-deviance windows, then run deviant PE dRSA for the three neural
%      targets above.
%   6) Save standard outputs + PE matrices + PE path plots + reproducibility
%      metadata with PE-specific post-deviance fields.
%
% Usage examples:
%   % Example 1: run from this folder.
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   PE_simulation_diff;
%
%   % Example 2: run from anywhere.
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   PE_simulation_diff;
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat
%   - experiment/input_files/MovDot_SubXX_predicted.mat (deviant branch)
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat (derived)
%
% Outputs:
%   - Standard counterfactual outputs:
%       simulations/output/subXX/<condition>/<dRSAtype>/{results,diagonal,matrices}
%   - PE matrices:
%       simulations/output/subXX/deviant/<dRSAtype>/matrices/PE/commCbar/*.png
%       simulations/output/subXX/deviant/<dRSAtype>/matrices/PE/sepCbar/*.png
%   - PE paths figure:
%       simulations/output/subXX/paths/*PE_paths*.png
%   - Results MAT fields:
%       dRSA_PE.neuralPE
%       dRSA_PE.neuralPredicted
%       dRSA_PE.neuralObserved
%
% Assumptions:
%   - Dot path arrays are trials x 2 x time, center-relative visual degrees.
%   - Trial lengths are uniform within each condition.
%   - This script is intentionally corr-only.
%
% Preparatory steps (user-owned):
%   - smoothing / concatenation / run averaging / masking before simulation
clear all
close all

%% Resolve repo paths and add dependencies
% Data flow: script location -> repo root -> addpath for simulation helpers.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
addpath(fullfile(repoRoot, 'simulations', 'functions'));
addpath(fullfile(repoRoot, 'simulations'));
addpath(fullfile(repoRoot, 'simulations', 'debug'));
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
% Data flow: deviantOnset fraction -> cut frame -> optional post-deviant-only trials.
cutPostDev = false; % true to keep only post-deviant samples in each trial
deviantOnset = 0.5; % fraction of trial duration marking deviant onset (0-1)
textScaleFactor = 2; % multiply axes/title/label text size across generated figures.
dRSAtypeToRun = 'corr'; % corr-only in this script (guard below).
% Existing-results behavior:
%   [] = interactive prompt, 1 = reuse/replot, 2 = rerun/overwrite, 3 = error out.
existingResultsAction = 2; % default non-interactive full rerun for batch execution.
suppressDispText = 1; % 1 to suppress function-level status text, 0 to show.
validConditions = {'nondeviant', 'deviant'};
if ~all(ismember(inputConditions, validConditions))
    error('inputConditions must be drawn from: %s', strjoin(validConditions, ', '));
end
if ~strcmpi(dRSAtypeToRun, 'corr')
    error('PE_simulation_diff only supports dRSAtypeToRun = ''corr''.');
end
if ~(isempty(existingResultsAction) || ...
        (isnumeric(existingResultsAction) && isscalar(existingResultsAction) ...
        && ismember(existingResultsAction, [1, 2, 3])))
    error('existingResultsAction must be [], 1, 2, or 3.');
end
simulationInputDir = fullfile(repoRoot, 'simulations', 'input');
outputDir = fullfile(repoRoot, 'simulations', 'output');

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
    if cutPostDev
        % Data flow: full paths -> deviant onset fraction -> post-deviant-only paths.
        [pathsByCondition.(condName).dot1, ~] = local_cut_postdeviant( ...
            pathsByCondition.(condName).dot1, deviantOnset);
        [pathsByCondition.(condName).dot2, ~] = local_cut_postdeviant( ...
            pathsByCondition.(condName).dot2, deviantOnset);
    end
end

figAllPaths = plot_paths_wrap( ...
    pathsByCondition.nondeviant.dot1, ...
    pathsByCondition.nondeviant.dot2, ...
    pathsByCondition.deviant.dot1, ...
    pathsByCondition.deviant.dot2, ...
    'SampleCount', n, ...
    'SampleRateHz', plotSampleRateHz, ...
    'Title', sprintf('Dot paths (sub%02d, all conditions)', participantNumber));

%% Run standard + PE dRSA per condition
% Data flow: condition list -> standard models -> optional deviant PE outputs.
for iCond = 1:numel(inputConditions)
    inputCondition = inputConditions{iCond};
    usePredictedModels = strcmp(inputCondition, 'deviant');
    reuseExistingResults = false;
    loadedResults = struct();
    existingActionUsed = [];
    cutFrame = [];
    dRSA_PE = struct(); % deviant-only; remains empty for nondeviant.
    peModelNames = {};
    peTargetLabels = {};
    figPEMatricesCommon = gobjects(0, 1);
    figPEMatricesSeparate = gobjects(0, 1);
    figPEPaths = gobjects(1, 1);
    pePathFile = '';
    peDot1 = [];
    peDot2 = [];
    peDirectionDot1 = [];
    peDirectionDot2 = [];
    peCutFrame = [];
    peCutSource = 'notApplicable';
    pePostDevEnforced = usePredictedModels;
    peObservedPredictedSizeMatchBeforeCut = [];
    peObservedPredictedSizeMatchAfterCut = [];
    peTrialLenBeforeCut = [];
    peTrialLenAfterCut = [];
    peMaskSubsampling = [];
    peMaskTrigger = [];
    peSubsamples = [];
    peTriggerOptions = struct();

    %% Resolve condition labels and output paths
    % Data flow: condition -> labels -> deterministic output filenames.
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
    % Normalize output labels by dropping run_default tokens and cleanup underscores.
    conditionLabel = regexprep(conditionLabel, '(^|_)run_default(_|$)', '$1');
    conditionLabel = regexprep(conditionLabel, '^_|_$', '');
    conditionLabel = regexprep(conditionLabel, '_+', '_');
    if isempty(conditionLabel)
        conditionLabel = 'unknown';
    end
    conditionSuffix = sprintf('(%s)', conditionLabel);
    dRSAtypeLabel = lower(dRSAtypeToRun);
    subjectOutputDir = fullfile(outputDir, subjectLabel);
    if cutPostDev
        subjectOutputDir = fullfile(subjectOutputDir, 'postDev');
    end
    conditionOutputDir = fullfile(subjectOutputDir, conditionLabel);
    analysisOutputDir = fullfile(conditionOutputDir, dRSAtypeLabel);
    pathsOutputDir = fullfile(subjectOutputDir, 'paths');
    resultsOutputDir = fullfile(analysisOutputDir, 'results');
    diagonalOutputDir = fullfile(analysisOutputDir, 'diagonal');
    matricesOutputDir = fullfile(analysisOutputDir, 'matrices');
    peMatricesOutputDir = fullfile(matricesOutputDir, 'PE');
    peMatricesCommDir = '';
    peMatricesSepDir = '';
    % Data flow: model-group flag -> concrete matrix folders for export.
    if usePredictedModels
        baseMatricesDir = fullfile(matricesOutputDir, 'base');
        predictedMatricesDir = fullfile(matricesOutputDir, 'predicted');
        matricesBaseCommDir = fullfile(baseMatricesDir, 'commCbar');
        matricesBaseSepDir = fullfile(baseMatricesDir, 'sepCbar');
        matricesPredCommDir = fullfile(predictedMatricesDir, 'commCbar');
        matricesPredSepDir = fullfile(predictedMatricesDir, 'sepCbar');
        peMatricesCommDir = fullfile(peMatricesOutputDir, 'commCbar');
        peMatricesSepDir = fullfile(peMatricesOutputDir, 'sepCbar');
    else
        matricesBaseCommDir = fullfile(matricesOutputDir, 'commCbar');
        matricesBaseSepDir = fullfile(matricesOutputDir, 'sepCbar');
        matricesPredCommDir = '';
        matricesPredSepDir = '';
    end
    local_ensure_dir(resultsOutputDir);
    local_ensure_dir(diagonalOutputDir);
    local_ensure_dir(matricesBaseCommDir);
    local_ensure_dir(matricesBaseSepDir);
    if usePredictedModels
        local_ensure_dir(matricesPredCommDir);
        local_ensure_dir(matricesPredSepDir);
        local_ensure_dir(peMatricesCommDir);
        local_ensure_dir(peMatricesSepDir);
    end
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
    runTagParts = [{subjectLabel, conditionLabel, lower(dRSAtypeToRun)}, optionalTags];
    resultsBase = strjoin([{'dRSA'}, runTagParts], '_');
    resultsBase = regexprep(resultsBase, '\s+', '');
    if usePredictedModels
        pePathFile = fullfile(pathsOutputDir, sprintf('%s_%s_PE_paths.png', resultsBase, conditionLabel));
    end
    resultsFile = fullfile(resultsOutputDir, [resultsBase '_results.mat']);
    legacyResultsFile = fullfile(subjectOutputDir, [resultsBase '_results.mat']);
    resultsFileToCheck = resultsFile;
    if ~isfile(resultsFile) && isfile(legacyResultsFile)
        % Backward-compat: allow reuse of flat outputs but write new plots to folders.
        resultsFileToCheck = legacyResultsFile;
    end

    %% Optional reuse of existing results
    % Data flow: existing results file -> parameter summary -> user-selected run mode.
    if isfile(resultsFileToCheck)
        existingMeta = load(resultsFileToCheck, 'repro', 'params');
        local_print_existing_match_summary( ...
            resultsFileToCheck, participantNumber, inputCondition, dRSAtypeToRun, ...
            shuffleDot2Trials, resampleSubsamples, resampleIterations, ...
            resampleSubsampleCount, resampleWithReplacement, ...
        usePredictedModels, plotSampleRateHz, cutPostDev, deviantOnset, existingMeta);
        existingActionUsed = local_resolve_existing_results_action(existingResultsAction);
        if existingActionUsed == 3
            error(['Existing results found and existingResultsAction = 3. ' ...
                'Stopping without rerun: %s'], resultsFileToCheck);
        elseif existingActionUsed == 1
            loadedResults = load(resultsFileToCheck);
            requiredVars = {'dRSA_position', 'dRSA_direction', ...
                'dRSA_diagonal_position', 'dRSA_diagonal_direction', 'params'};
            if usePredictedModels
                % Deviant PE analyses require dRSA_PE with all three targets.
                requiredVars = [requiredVars, {'dRSA_PE'}];
            end
            missingVars = requiredVars(~isfield(loadedResults, requiredVars));
            if isempty(missingVars) && usePredictedModels
                peFieldNames = {'neuralPE', 'neuralPredicted', 'neuralObserved'};
                if ~isstruct(loadedResults.dRSA_PE) || ...
                        ~all(isfield(loadedResults.dRSA_PE, peFieldNames))
                    missingVars = [missingVars, {'dRSA_PE.(neuralPE/neuralPredicted/neuralObserved)'}];
                end
            end
            if isempty(missingVars)
                reuseExistingResults = true;
            else
                warning(['Existing file is missing required variables (%s). ' ...
                    'Falling back to full re-run.'], strjoin(missingVars, ', '));
            end
        else
            % existingActionUsed == 2 -> full rerun by design.
            fprintf('existingResultsAction=2: full rerun for %s.\n', conditionSuffix);
        end
    end

    if reuseExistingResults
        % Load saved matrices/metadata and skip expensive dRSA recomputation.
        dRSA_position = loadedResults.dRSA_position;
        dRSA_direction = loadedResults.dRSA_direction;
        dRSA_diagonal_position = loadedResults.dRSA_diagonal_position;
        dRSA_diagonal_direction = loadedResults.dRSA_diagonal_direction;
        params = loadedResults.params;
        if isfield(loadedResults, 'paramsDiagonal')
            paramsDiagonal = loadedResults.paramsDiagonal;
        else
            paramsDiagonal = params;
        end
        if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'paramsCore')
            paramsCore = loadedResults.repro.paramsCore;
        else
            paramsCore = params;
        end
        if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'modelNames')
            modelNames = loadedResults.repro.modelNames;
        else
            modelNames = local_default_model_names(size(dRSA_position, 3));
        end
        if usePredictedModels && isfield(loadedResults, 'dRSA_PE')
            dRSA_PE = loadedResults.dRSA_PE;
            peTargetLabels = fieldnames(dRSA_PE);
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peModelNames')
                peModelNames = loadedResults.repro.peModelNames;
            else
                peModelNames = {'PE position dot1', 'PE position dot2', ...
                    'PE direction dot1', 'PE direction dot2'};
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'pePathFile')
                pePathFile = loadedResults.repro.pePathFile;
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peDot1') ...
                    && isfield(loadedResults.repro, 'peDot2')
                peDot1 = loadedResults.repro.peDot1;
                peDot2 = loadedResults.repro.peDot2;
                figPEPaths = local_plot_pe_paths(peDot1, peDot2, plotSampleRateHz, ...
                    sprintf('PE paths (sub%02d, %s)', participantNumber, inputCondition), 'on');
            end
        end
    else
        %% Load condition-specific dot paths (center-relative visual degrees)
        % Data flow: condition .mat -> center-relative dot arrays -> placeholders for simulation.
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
    dot2Trials = dot2YellowPaths;

    %% Optional trial shuffling (dot2 relative to dot1)
    % Data flow: dot2Trials -> permuted trial order -> decoupled dot1/dot2 pairing.
    % Purpose: test whether dot1-dot2 structure is driven by trial-locked coupling.
    shuffleOrder = []; % store permutation for reproducibility (empty when unshuffled).
    if shuffleDot2Trials
        rng('shuffle'); % avoid reusing the subject-seeded RNG for the shuffle
        shuffleOrder = randperm(size(dot2Trials, 1));
        dot2Trials = dot2Trials(shuffleOrder, :, :);
    end

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
    end

    %% Optional post-deviant truncation
    % Data flow: full trial paths -> deviant onset -> post-deviant-only trials.
    if cutPostDev
        [dot1Trials, cutFrame] = local_cut_postdeviant(dot1Trials, deviantOnset);
        [dot2Trials, ~] = local_cut_postdeviant(dot2Trials, deviantOnset);
        if usePredictedModels
            [dot1TrialsPredicted, ~] = local_cut_postdeviant(dot1TrialsPredicted, deviantOnset);
            [dot2TrialsPredicted, ~] = local_cut_postdeviant(dot2TrialsPredicted, deviantOnset);
        end
    end

    %% Prepare data, models, and trial-locked subsamples
    % Data flow: trial arrays -> per-dot position + direction models -> concatenated time series.
    % Position models use raw x/y paths; direction models use per-timepoint angles.
    plotConcat = 0; % suppress concatenation plots for batch runs
    dataPosition = dRSA_concatenate(dot1Trials, [], plotConcat, ...
        'suppressDispText', suppressDispText); % see documentation with help dRSA_concatenate

    % Position models (per dot).
    modelPositionDot1 = dataPosition;
    modelPositionDot2 = dRSA_concatenate(dot2Trials, [], plotConcat, ...
        'suppressDispText', suppressDispText);
    modelPositionDot1Predicted = [];
    modelPositionDot2Predicted = [];
    modelDirectionDot1Predicted = [];
    modelDirectionDot2Predicted = [];
    if usePredictedModels
        modelPositionDot1Predicted = dRSA_concatenate(dot1TrialsPredicted, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        modelPositionDot2Predicted = dRSA_concatenate(dot2TrialsPredicted, [], plotConcat, ...
            'suppressDispText', suppressDispText);
    end

    % Direction models (per dot).
    % Data flow: dot positions -> frame-to-frame displacement -> unit vectors.
    dot1Direction = local_positions_to_direction(dot1Trials);
    dot2Direction = local_positions_to_direction(dot2Trials);
    modelDirectionDot1 = dRSA_concatenate(dot1Direction, [], plotConcat, ...
        'suppressDispText', suppressDispText);
    modelDirectionDot2 = dRSA_concatenate(dot2Direction, [], plotConcat, ...
        'suppressDispText', suppressDispText);
    dataDirection = modelDirectionDot1; % use dot 1 direction as observed neural stream.

    % PE branch placeholders.
    dataPositionPredicted = [];
    dataDirectionPredicted = [];
    dataPositionObservedPEWindow = [];
    dataDirectionObservedPEWindow = [];
    dataPositionPE = [];
    dataDirectionPE = [];
    modelPE = {};
    paramsPE = struct();

    if usePredictedModels
        predictedDot1Direction = local_positions_to_direction(dot1TrialsPredicted);
        predictedDot2Direction = local_positions_to_direction(dot2TrialsPredicted);
        modelDirectionDot1Predicted = dRSA_concatenate(predictedDot1Direction, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        modelDirectionDot2Predicted = dRSA_concatenate(predictedDot2Direction, [], plotConcat, ...
            'suppressDispText', suppressDispText);

        %% Enforce post-deviance windows for PE branch (mandatory)
        % Data flow: observed+predicted streams -> post-deviance cut -> PE models/targets.
        peObservedDot1 = dot1Trials;
        peObservedDot2 = dot2Trials;
        pePredictedDot1 = dot1TrialsPredicted;
        pePredictedDot2 = dot2TrialsPredicted;
        peObservedPredictedSizeMatchBeforeCut = ...
            isequal(size(peObservedDot1), size(pePredictedDot1)) && ...
            isequal(size(peObservedDot2), size(pePredictedDot2));
        if ~peObservedPredictedSizeMatchBeforeCut
            error('Observed and predicted streams must match before PE post-deviance cut.');
        end
        peTrialLenBeforeCut = size(peObservedDot1, 3);

        if cutPostDev
            % Global cut already applied above; keep traceability for PE metadata.
            peCutSource = 'globalCutPostDev';
            peCutFrame = cutFrame;
        else
            % Apply mandatory PE-only post-deviance cut when global cut is disabled.
            [peObservedDot1, peCutFrame] = local_cut_postdeviant(peObservedDot1, deviantOnset);
            [peObservedDot2, ~] = local_cut_postdeviant(peObservedDot2, deviantOnset);
            [pePredictedDot1, ~] = local_cut_postdeviant(pePredictedDot1, deviantOnset);
            [pePredictedDot2, ~] = local_cut_postdeviant(pePredictedDot2, deviantOnset);
            peCutSource = 'peBranchCut';
        end

        peObservedPredictedSizeMatchAfterCut = ...
            isequal(size(peObservedDot1), size(pePredictedDot1)) && ...
            isequal(size(peObservedDot2), size(pePredictedDot2));
        if ~peObservedPredictedSizeMatchAfterCut
            error('Observed and predicted streams diverged after PE post-deviance cut.');
        end
        peTrialLenAfterCut = size(peObservedDot1, 3);

        %% Build signal-level PE models (observed - predicted)
        % Data flow: PE-cut observed/predicted position+direction streams -> PE model set.
        peDot1 = peObservedDot1 - pePredictedDot1;
        peDot2 = peObservedDot2 - pePredictedDot2;
        modelPEPositionDot1 = dRSA_concatenate(peDot1, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        modelPEPositionDot2 = dRSA_concatenate(peDot2, [], plotConcat, ...
            'suppressDispText', suppressDispText);

        peObservedDirectionDot1 = local_positions_to_direction(peObservedDot1);
        peObservedDirectionDot2 = local_positions_to_direction(peObservedDot2);
        pePredictedDirectionDot1 = local_positions_to_direction(pePredictedDot1);
        pePredictedDirectionDot2 = local_positions_to_direction(pePredictedDot2);
        peDirectionDot1 = peObservedDirectionDot1 - pePredictedDirectionDot1;
        peDirectionDot2 = peObservedDirectionDot2 - pePredictedDirectionDot2;
        modelPEDirectionDot1 = dRSA_concatenate(peDirectionDot1, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        modelPEDirectionDot2 = dRSA_concatenate(peDirectionDot2, [], plotConcat, ...
            'suppressDispText', suppressDispText);

        modelPE = {modelPEPositionDot1, modelPEPositionDot2, ...
            modelPEDirectionDot1, modelPEDirectionDot2};
        peModelNames = {'PE position dot1', 'PE position dot2', ...
            'PE direction dot1', 'PE direction dot2'};

        %% Build PE neural targets from the same post-deviance window
        % Data flow: PE-cut observed/predicted streams -> observed/predicted/PE targets.
        dataPositionObservedPEWindow = dRSA_concatenate(peObservedDot1, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        dataDirectionObservedPEWindow = dRSA_concatenate(peObservedDirectionDot1, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        dataPositionPredicted = dRSA_concatenate(pePredictedDot1, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        dataDirectionPredicted = dRSA_concatenate(pePredictedDirectionDot1, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        dataPositionPE = dataPositionObservedPEWindow - dataPositionPredicted;
        dataDirectionPE = dataDirectionObservedPEWindow - dataDirectionPredicted;

        % Additional PE paths figure (deviant only), saved under subject-level paths.
        figPEPaths = local_plot_pe_paths(peDot1, peDot2, plotSampleRateHz, ...
            sprintf('PE paths (sub%02d, %s)', participantNumber, inputCondition), 'on');
    end

    %% Build trial-locked subsamples
    % Data flow: trial length -> trigger mask -> subsamples tensor for dRSA.
    trialLen = size(dot1Trials, 3);
    totalTime = size(dataPosition, 2);
    [maskSubsampling, maskTrigger, opt, subsamples] = local_build_trial_subsamples(trialLen, totalTime);

    %% Optional subsample resampling (memory control)
    % Data flow: full trial subsamples -> (optional) resampled subsets -> iteration tensor.
    if resampleSubsamples
        subsamples = local_resample_subsamples( ...
            subsamples, resampleSubsampleCount, resampleIterations, resampleWithReplacement);
    end

    %% Build PE-specific subsamples (deviant only)
    % Data flow: PE post-deviance trial length -> PE-specific trigger tensor.
    if usePredictedModels
        peTotalTime = size(dataPositionObservedPEWindow, 2);
        [peMaskSubsampling, peMaskTrigger, peTriggerOptions, peSubsamples] = ...
            local_build_trial_subsamples(peTrialLenAfterCut, peTotalTime);
        if resampleSubsamples
            peSubsamples = local_resample_subsamples( ...
                peSubsamples, resampleSubsampleCount, resampleIterations, resampleWithReplacement);
        end
    end

    %% Configure dRSA models and parameters
    % Data flow: position/direction models -> dRSA parameters -> per-iteration core calls.
    Y = dataPosition; % position data from dot 1
    model = {modelPositionDot1, modelPositionDot2, modelDirectionDot1, modelDirectionDot2};
    modelNames = {'position dot1', 'position dot2', 'direction dot1', 'direction dot2'};
    if usePredictedModels
        model = {modelPositionDot1, modelPositionDot2, modelDirectionDot1, modelDirectionDot2, ...
            modelPositionDot1Predicted, modelPositionDot2Predicted, ...
            modelDirectionDot1Predicted, modelDirectionDot2Predicted};
        modelNames = {'position dot1', 'position dot2', 'direction dot1', 'direction dot2', ...
            'predicted position dot1', 'predicted position dot2', ...
            'predicted direction dot1', 'predicted direction dot2'};
    end
    
    params.modelNames = modelNames;
    % Keep corr-only mode explicit in this script version.
    params.nIter = size(subsamples, 3);  % how many iterations?
    params.fs = plotSampleRateHz; % framerate or how many samples fit into 1 second
    % Use an integer sample window for averaging to avoid non-integer sizes.
    avgHalfWindowSamples = floor((trialLen - 1) / 2);
    params.AverageTime = avgHalfWindowSamples / params.fs; % in s
    params.modelToTest = 1:numel(model); % array of models to test
    params.Var = 0.1; % how much variance?
    params.modelDistMeasure = {'euclidean', 'euclidean', 'cosine', 'cosine'};
    if usePredictedModels
        params.modelDistMeasure = [params.modelDistMeasure, ...
            {'euclidean', 'euclidean', 'cosine', 'cosine'}];
    end
    params.neuralDistMeasure = params.modelDistMeasure{1}; % pdist expects a string/char or function handle.

    params.dRSAtype = dRSAtypeToRun; % corr-only guard is enforced at config time.
    params.PCR.AdditionalPCA = 1;
    % Keep PCR fields populated for API compatibility with downstream calls.
    params.PCR.RegressAutocor = 1; % 1 to get the autocorr
    params.PCR.RessModel = 1;

    params.modeltoRegressout = cell(1, numel(model));
    for iModel = 1:numel(model)
        params.modeltoRegressout{iModel} = setdiff(1:numel(model), iModel);
    end

    % Autocorrelation border is only needed for non-corr dRSA.
    Autocorrborder = []; % empty signals "no border" to downstream code paths.
    if ~strcmp(params.dRSAtype, 'corr')
        Autocorrborder = dRSA_border(model, subsamples, params, ...
            'suppressDispText', suppressDispText, 'plotAutocorr', 'on');
    end

    %% Configure PE dRSA parameters (deviant only)
    % Data flow: shared params -> PE-specific model list + fixed model metrics.
    if usePredictedModels
        paramsPE = params;
        paramsPE.nIter = size(peSubsamples, 3);
        avgHalfWindowSamplesPE = floor((peTrialLenAfterCut - 1) / 2);
        paramsPE.AverageTime = avgHalfWindowSamplesPE / paramsPE.fs;
        paramsPE.modelNames = peModelNames;
        paramsPE.modelToTest = 1:numel(modelPE);
        paramsPE.modelDistMeasure = {'euclidean', 'euclidean', 'cosine', 'cosine'};
        paramsPE.modeltoRegressout = cell(1, numel(modelPE));
        for iModel = 1:numel(modelPE)
            paramsPE.modeltoRegressout{iModel} = setdiff(1:numel(modelPE), iModel);
        end
    end

    %% Run dRSA for position data
    % Data flow: position data -> distance metric aligned to position model -> dRSA.
    params.neuralDistMeasure = params.modelDistMeasure{1};
    dRSA_Iter = [];
    for iIter = 1:params.nIter
        CurrSubsamples = subsamples(:, :, iIter); % subsamples is nSubsamples*subSampleDuration*iterations
        dRSAma = dRSA_coreFunction(Y, model, params, ...
            'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);
        dRSA_Iter(iIter, :, :, :) = dRSAma;
    end
    dRSA_position = mean(dRSA_Iter, 1);
    dRSA_position = reshape(dRSA_position, size(dRSA_position, 2), size(dRSA_position, 3), size(dRSA_position, 4));

    %% Run dRSA for direction data
    % Data flow: direction data -> distance metric aligned to direction model -> dRSA.
    params.neuralDistMeasure = params.modelDistMeasure{3};
    Y = dataDirection; % direction data from dot 1
    dRSA_Iter = [];
    for iIter = 1:params.nIter
        CurrSubsamples = subsamples(:, :, iIter); % subsamples is nSubsamples*subSampleDuration*iterations
        dRSAma = dRSA_coreFunction(Y, model, params, ...
            'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);
        dRSA_Iter(iIter, :, :, :) = dRSAma;
    end
    dRSA_direction = mean(dRSA_Iter, 1);
    dRSA_direction = reshape(dRSA_direction, size(dRSA_direction, 2), size(dRSA_direction, 3), size(dRSA_direction, 4));

    %% Run deviant-only PE dRSA using three neural targets
    % Data flow: (PE neural target stream, PE models, PE subsamples) -> dRSA PE matrices.
    if usePredictedModels
        peTargetLabels = {'neuralPE', 'neuralPredicted', 'neuralObserved'};
        peTargetPositionData = {dataPositionPE, dataPositionPredicted, dataPositionObservedPEWindow};
        peTargetDirectionData = {dataDirectionPE, dataDirectionPredicted, dataDirectionObservedPEWindow};
        for iTarget = 1:numel(peTargetLabels)
            currTargetLabel = peTargetLabels{iTarget};

            % Position row for current PE neural target.
            paramsPE.neuralDistMeasure = paramsPE.modelDistMeasure{1};
            Y = peTargetPositionData{iTarget};
            dRSA_Iter = [];
            for iIter = 1:paramsPE.nIter
                CurrSubsamples = peSubsamples(:, :, iIter);
                dRSAma = dRSA_coreFunction(Y, modelPE, paramsPE, ...
                    'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);
                dRSA_Iter(iIter, :, :, :) = dRSAma;
            end
            dRSA_position_pe = mean(dRSA_Iter, 1);
            dRSA_position_pe = reshape(dRSA_position_pe, ...
                size(dRSA_position_pe, 2), size(dRSA_position_pe, 3), size(dRSA_position_pe, 4));

            % Direction row for current PE neural target.
            paramsPE.neuralDistMeasure = paramsPE.modelDistMeasure{3};
            Y = peTargetDirectionData{iTarget};
            dRSA_Iter = [];
            for iIter = 1:paramsPE.nIter
                CurrSubsamples = peSubsamples(:, :, iIter);
                dRSAma = dRSA_coreFunction(Y, modelPE, paramsPE, ...
                    'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);
                dRSA_Iter(iIter, :, :, :) = dRSAma;
            end
            dRSA_direction_pe = mean(dRSA_Iter, 1);
            dRSA_direction_pe = reshape(dRSA_direction_pe, ...
                size(dRSA_direction_pe, 2), size(dRSA_direction_pe, 3), size(dRSA_direction_pe, 4));

            % Save current neural-target matrices into the PE output struct.
            dRSA_PE.(currTargetLabel) = struct( ...
                'position', dRSA_position_pe, ...
                'direction', dRSA_direction_pe, ...
                'modelNames', {peModelNames}, ...
                'rowNames', {{'Data: position PE target', 'Data: direction PE target'}});
        end
    end

    % Preserve parameters used for dRSA before reusing them for diagonal averaging.
    paramsCore = params;
    end

    %% Plot dRSA matrices (position vs direction data)
    % Data flow: dRSA per data type (time x time x model) -> grid of heatmaps.
    % Ticks are fixed every 0.5 s and include midpoint label 'M' for readability.
    nModels = size(dRSA_position, 3);
    nTime = size(dRSA_position, 1);
    timeSeconds = (0:nTime-1) / paramsCore.fs;
    dRSAAll = {dRSA_position, dRSA_direction};
    rowNames = {'Data: position dot1', 'Data: direction dot1'};
    baseModelIdx = 1:4;
    predictedModelIdx = 5:nModels;
    hasPredictedModels = nModels > numel(baseModelIdx);
    figLabels = {};
    if hasPredictedModels
        figLabels = {'base', 'predicted'};
    else
        figLabels = {'base'};
    end
    cbarLimits = [0 1];
    % Always create both matrix variants on screen for direct comparison.
    figMatricesCommon = local_plot_dRSA_matrices( ...
        dRSAAll, modelNames, rowNames, timeSeconds, figLabels, ...
        baseModelIdx, predictedModelIdx, true, cbarLimits, 'on');
    figMatricesSeparate = local_plot_dRSA_matrices( ...
        dRSAAll, modelNames, rowNames, timeSeconds, figLabels, ...
        baseModelIdx, predictedModelIdx, false, cbarLimits, 'on');

    %% Plot deviant-only PE matrices (three neural targets, each 2x4)
    % Data flow: PE dRSA struct -> per-target 2x4 matrix figures -> handle arrays.
    if usePredictedModels
        if isempty(peModelNames)
            peModelNames = {'PE position dot1', 'PE position dot2', ...
                'PE direction dot1', 'PE direction dot2'};
        end
        % PE matrices can use a shorter post-deviance window than standard matrices.
        peNTime = size(dRSA_PE.neuralPE.position, 1);
        peTimeSeconds = (0:peNTime-1) / paramsCore.fs;
        figPEMatricesCommon = local_plot_pe_matrices( ...
            dRSA_PE, peModelNames, peTimeSeconds, true, cbarLimits, 'on');
        figPEMatricesSeparate = local_plot_pe_matrices( ...
            dRSA_PE, peModelNames, peTimeSeconds, false, cbarLimits, 'on');
    end

    %% Plot diagonal
    % Data flow: dRSA matrices -> diagonal averages -> summary lines.
    paramsDiagonal = paramsCore;
    paramsDiagonal.AverageTime = 2; % in s, how much should be left and right of the zerolag middle?
    dRSA_diagonal_position = dRSA_average(dRSA_position, paramsDiagonal);
    dRSA_diagonal_direction = dRSA_average(dRSA_direction, paramsDiagonal);

    figLag = figure('Name', 'dRSA diagonal', 'NumberTitle', 'off');
    hold on;
    plot(dRSA_diagonal_position(1, :), 'LineWidth', 1.5);
    plot(dRSA_diagonal_direction(3, :), 'LineWidth', 1.5);
    hold off;
    xlabel('Time (samples)');
    ylabel('dRSA (diagonal)');
    legend(rowNames, 'Location', 'best');

    %% Save matrices, plots, and reproducibility metadata (300 dpi)
    % Data flow: dRSA outputs + run metadata -> descriptive filenames -> .mat/.png outputs.
    if ~reuseExistingResults
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
        repro.existingResultsActionConfig = existingResultsAction;
        repro.existingResultsActionUsed = existingActionUsed;
        repro.shuffleDot2Trials = shuffleDot2Trials;
        repro.shuffleOrder = shuffleOrder;
        repro.resampleSubsamples = resampleSubsamples;
        repro.resampleIterations = resampleIterations;
        repro.resampleSubsampleCount = resampleSubsampleCount;
        repro.resampleWithReplacement = resampleWithReplacement;
        repro.cutPostDev = cutPostDev;
        repro.deviantOnset = deviantOnset;
        repro.cutFrame = cutFrame;
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
        repro.dataDirection = dataDirection;
        repro.models = model;
        repro.autocorrBorder = Autocorrborder;
        repro.predictedSimulationFile = predictedSimulationFile;
        repro.peModelNames = peModelNames;
        repro.peTargetLabels = peTargetLabels;
        repro.peDot1 = peDot1;
        repro.peDot2 = peDot2;
        repro.peDirectionDot1 = peDirectionDot1;
        repro.peDirectionDot2 = peDirectionDot2;
        repro.pePathFile = pePathFile;
        repro.pePostDevEnforced = pePostDevEnforced;
        repro.peCutSource = peCutSource;
        repro.peCutFrame = peCutFrame;
        repro.peTrialLenBeforeCut = peTrialLenBeforeCut;
        repro.peTrialLenAfterCut = peTrialLenAfterCut;
        repro.peObservedPredictedSizeMatchBeforeCut = peObservedPredictedSizeMatchBeforeCut;
        repro.peObservedPredictedSizeMatchAfterCut = peObservedPredictedSizeMatchAfterCut;
        repro.peMaskSubsampling = peMaskSubsampling;
        repro.peMaskTrigger = peMaskTrigger;
        repro.peSubsamples = peSubsamples;
        repro.peTriggerOptions = peTriggerOptions;

        save(resultsFile, ...
            'dRSA_position', 'dRSA_direction', ...
            'dRSA_diagonal_position', 'dRSA_diagonal_direction', ...
            'dRSA_PE', 'params', 'paramsDiagonal', 'repro');
    else
        fprintf('Using existing dataset for %s; only re-plotting figures.\n', conditionSuffix);
    end

    % Maximize figures before export so matrices occupy more pixels in saved PNGs.
    for iFigPrep = 1:numel(figMatricesCommon)
        local_scale_figure_text(figMatricesCommon(iFigPrep), textScaleFactor);
        local_prepare_figure_for_export(figMatricesCommon(iFigPrep));
    end
    for iFigPrep = 1:numel(figMatricesSeparate)
        local_scale_figure_text(figMatricesSeparate(iFigPrep), textScaleFactor);
        local_prepare_figure_for_export(figMatricesSeparate(iFigPrep));
    end
    for iFigPrep = 1:numel(figPEMatricesCommon)
        local_scale_figure_text(figPEMatricesCommon(iFigPrep), textScaleFactor);
        local_prepare_figure_for_export(figPEMatricesCommon(iFigPrep));
    end
    for iFigPrep = 1:numel(figPEMatricesSeparate)
        local_scale_figure_text(figPEMatricesSeparate(iFigPrep), textScaleFactor);
        local_prepare_figure_for_export(figPEMatricesSeparate(iFigPrep));
    end
    % Keep lag plots at their native size/text for readability consistency.
    if iCond == 1
        local_scale_figure_text(figAllPaths, textScaleFactor);
        local_prepare_figure_for_export(figAllPaths);
    end
    if usePredictedModels && isgraphics(figPEPaths, 'figure')
        local_scale_figure_text(figPEPaths, textScaleFactor);
        local_prepare_figure_for_export(figPEPaths);
    end

    % Save both colorbar variants (guard against invalid figure handles).
    if numel(figLabels) == 1
        local_safe_print(figMatricesCommon(1), ...
            fullfile(matricesBaseCommDir, [resultsBase '_matrices_commCbar.png']));
        local_safe_print(figMatricesSeparate(1), ...
            fullfile(matricesBaseSepDir, [resultsBase '_matrices_sepCbar.png']));
    else
        local_safe_print(figMatricesCommon(1), ...
            fullfile(matricesBaseCommDir, [resultsBase '_matrices_base_commCbar.png']));
        local_safe_print(figMatricesCommon(2), ...
            fullfile(matricesPredCommDir, [resultsBase '_matrices_predicted_commCbar.png']));
        local_safe_print(figMatricesSeparate(1), ...
            fullfile(matricesBaseSepDir, [resultsBase '_matrices_base_sepCbar.png']));
        local_safe_print(figMatricesSeparate(2), ...
            fullfile(matricesPredSepDir, [resultsBase '_matrices_predicted_sepCbar.png']));
    end
    if usePredictedModels
        peTargetOrder = {'neuralPE', 'neuralPredicted', 'neuralObserved'};
        for iTarget = 1:numel(peTargetOrder)
            targetLabel = peTargetOrder{iTarget};
            local_safe_print(figPEMatricesCommon(iTarget), ...
                fullfile(peMatricesCommDir, sprintf('%s_PE_%s_commCbar.png', resultsBase, targetLabel)));
            local_safe_print(figPEMatricesSeparate(iTarget), ...
                fullfile(peMatricesSepDir, sprintf('%s_PE_%s_sepCbar.png', resultsBase, targetLabel)));
        end
        pePathFile = fullfile(pathsOutputDir, sprintf('%s_%s_PE_paths.png', resultsBase, conditionLabel));
        local_safe_print(figPEPaths, pePathFile);
    end
    local_safe_print(figLag, fullfile(diagonalOutputDir, [resultsBase '_diagonal.png']));
    if iCond == 1
        pathsBase = strjoin({'paths', subjectLabel, 'allConditions'}, '_');
        pathsBase = regexprep(pathsBase, '\s+', '');
        print(figAllPaths, fullfile(pathsOutputDir, [pathsBase '_all_conditions.png']), ...
            '-dpng', '-r300');
    end
end

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

function local_safe_print(figHandle, outPath)
% LOCAL_SAFE_PRINT Print a figure if the handle is valid (avoid print errors).
% Data flow: figure handle -> validity check -> print to disk.
if isempty(figHandle) || ~isgraphics(figHandle, 'figure')
    warning('Skipping print; invalid figure handle for %s', outPath);
    return;
end
print(figHandle, outPath, '-dpng', '-r300');
end

function direction = local_positions_to_direction(paths)
% LOCAL_POSITIONS_TO_DIRECTION Convert x/y path samples to unit-vector directions.
%
% Inputs:
%   paths : trials x 2 x time center-relative positions.
%
% Outputs:
%   direction : trials x 2 x time where columns are [cos(theta), sin(theta)].
%
% Data flow:
%   positions -> frame-to-frame displacement -> angle -> unit-vector direction.
if isempty(paths)
    direction = paths;
    return;
end
if size(paths, 3) < 2
    % Edge case: one-sample trials have no displacement; keep neutral vectors.
    direction = zeros(size(paths));
    return;
end
dx = diff(paths(:, 1, :), 1, 3);
dy = diff(paths(:, 2, :), 1, 3);
dx = cat(3, dx(:, :, 1), dx); % keep time length aligned to original paths.
dy = cat(3, dy(:, :, 1), dy);
angle = atan2(dy, dx);
direction = cat(2, cos(angle), sin(angle));
end

function [maskSubsampling, maskTrigger, opt, subsamples] = ...
        local_build_trial_subsamples(trialLen, totalTime)
% LOCAL_BUILD_TRIAL_SUBSAMPLES Build trial-locked subsamples for dRSA.
%
% Inputs:
%   trialLen  : number of samples per trial window.
%   totalTime : total concatenated sample length.
%
% Outputs:
%   maskSubsampling : logical mask of eligible time points.
%   maskTrigger     : logical trigger vector with one trigger per trial start.
%   opt             : dRSA_triggered_subsampling options used.
%   subsamples      : nSubSamples x subSampleDuration x nIter trigger windows.
if mod(totalTime, trialLen) ~= 0
    error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
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

subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);
end

function subsamples = local_resample_subsamples( ...
        subsamples, resampleSubsampleCount, resampleIterations, resampleWithReplacement)
% LOCAL_RESAMPLE_SUBSAMPLES Optional resampling wrapper for subsample tensors.
%
% Inputs:
%   subsamples              : base subsample tensor from triggered sampling.
%   resampleSubsampleCount  : number of subsamples per iteration.
%   resampleIterations      : number of iterations.
%   resampleWithReplacement : true for bootstrap-style resampling.
%
% Outputs:
%   subsamples : resampled subsample tensor (count x duration x iterations).
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

function figMatrices = local_plot_dRSA_matrices( ...
        dRSAAll, modelNames, rowNames, timeSeconds, figLabels, ...
        baseModelIdx, predictedModelIdx, useCommonCbarLimits, cbarLimits, figVisibility)
% LOCAL_PLOT_DRSA_MATRICES Render dRSA matrix figures with standardized ticks.
%
% Usage example:
%   figH = local_plot_dRSA_matrices(dRSAAll, modelNames, rowNames, timeSeconds, ...
%       {'base'}, 1:4, [], true, [0 1], 'on');
%
% Inputs:
%   dRSAAll             : 1x2 cell {dRSA_position, dRSA_direction}.
%   modelNames          : 1xnModels cell of model labels.
%   rowNames            : 1x2 cell of row labels.
%   timeSeconds         : 1xnTime time axis in seconds.
%   figLabels           : figure groups ('base', 'predicted').
%   baseModelIdx        : indices of base models in the third dRSA dimension.
%   predictedModelIdx     : indices of predicted models in the third dRSA dimension.
%   useCommonCbarLimits : true for shared row colorbars + fixed limits.
%   cbarLimits          : [min max] colorbar limits used when sharing.
%   figVisibility       : 'on' to show figures, 'off' to keep hidden.
%
% Output:
%   figMatrices         : array of figure handles, one per entry in figLabels.
%
% Data flow:
%   dRSA tensors -> per-group tiled figures -> standardized ticks every 0.5 s
%   + midpoint marker 'M' -> optional shared or per-subplot colorbars.

figMatrices = gobjects(numel(figLabels), 1);
midSample = ceil(numel(timeSeconds) / 2);
midTimeSeconds = timeSeconds(midSample);
timeTickStepSeconds = 0.5;
timeTicks = 0:timeTickStepSeconds:timeSeconds(end);
if isempty(timeTicks) || abs(timeTicks(end) - timeSeconds(end)) > eps(max(timeSeconds(end), 1))
    timeTicks = [timeTicks, timeSeconds(end)];
end
timeTicks = unique([timeTicks, midTimeSeconds]);
timeTickLabels = cellstr(compose('%.3g', timeTicks));
[~, midTickIdx] = min(abs(timeTicks - midTimeSeconds));
timeTickLabels{midTickIdx} = 'M';

for iFig = 1:numel(figLabels)
    if strcmp(figLabels{iFig}, 'base')
        modelIdx = baseModelIdx;
    else
        modelIdx = predictedModelIdx;
    end

    figMatrices(iFig) = figure( ...
        'Name', sprintf('dRSA matrices (%s)', figLabels{iFig}), ...
        'NumberTitle', 'off', ...
        'Visible', figVisibility);
    tiledlayout(2, numel(modelIdx), 'Padding', 'compact', 'TileSpacing', 'compact');
    rowAxes = gobjects(numel(dRSAAll), numel(modelIdx));

    for iRow = 1:numel(dRSAAll)
        for iModel = 1:numel(modelIdx)
            tileIdx = (iRow - 1) * numel(modelIdx) + iModel;
            modelCol = modelIdx(iModel);
            nexttile(tileIdx);
            imagesc(timeSeconds, timeSeconds, dRSAAll{iRow}(:, :, modelCol));
            set(gca, 'YDir', 'normal');
            axis image;
            if ~useCommonCbarLimits
                colorbar;
            end
            title(sprintf('dRSA %s', modelNames{modelCol}));
            % dRSA uses corr(modelRDM, neuralRDM): x = neural/data time, y = model time.
            xlabel(sprintf('Time in %s (s)', rowNames{iRow}));
            ylabel(sprintf('Time in %s (s)', modelNames{modelCol}));
            xticks(timeTicks);
            yticks(timeTicks);
            xticklabels(timeTickLabels);
            yticklabels(timeTickLabels);
            local_apply_offset_ticks(gca, timeTicks, midTickIdx);
            hold on;
            plot(timeSeconds, timeSeconds, 'w-', 'LineWidth', 1);
            hold off;
            if useCommonCbarLimits
                caxis(cbarLimits);
            end
            rowAxes(iRow, iModel) = gca;
        end
    end

    if useCommonCbarLimits
        for iRow = 1:numel(dRSAAll)
            cbarHandle = colorbar(rowAxes(iRow, end), 'Location', 'eastoutside');
            cbarHandle.Limits = cbarLimits;
        end
    end
end
end

function figPEMatrices = local_plot_pe_matrices( ...
        dRSA_PE, peModelNames, timeSeconds, useCommonCbarLimits, cbarLimits, figVisibility)
% LOCAL_PLOT_PE_MATRICES Plot deviant PE dRSA matrices for each PE neural target.
%
% Purpose:
%   Render one figure per PE neural target. Each figure always uses a
%   2x4 layout: row 1 = position neural stream, row 2 = direction neural
%   stream, and columns = the four fixed PE models.
%
% Inputs:
%   dRSA_PE             : struct with fields neuralPE/neuralPredicted/neuralObserved.
%   peModelNames        : 1x4 model labels for PE model columns.
%   timeSeconds         : 1xnTime time axis in seconds.
%   useCommonCbarLimits : true for shared [min max] per-row colorbars.
%   cbarLimits          : [min max] limits when useCommonCbarLimits=true.
%   figVisibility       : 'on' or 'off'.
%
% Output:
%   figPEMatrices       : 3x1 array of figure handles in fixed target order.
%
% Data flow:
%   PE dRSA struct -> target-specific matrix pairs -> standardized 2x4 figures.

targetOrder = {'neuralPE', 'neuralPredicted', 'neuralObserved'};
rowNames = {'PE neural target (position)', 'PE neural target (direction)'};
figPEMatrices = gobjects(numel(targetOrder), 1);

midSample = ceil(numel(timeSeconds) / 2);
midTimeSeconds = timeSeconds(midSample);
timeTickStepSeconds = 0.5;
timeTicks = 0:timeTickStepSeconds:timeSeconds(end);
if isempty(timeTicks) || abs(timeTicks(end) - timeSeconds(end)) > eps(max(timeSeconds(end), 1))
    timeTicks = [timeTicks, timeSeconds(end)];
end
timeTicks = unique([timeTicks, midTimeSeconds]);
timeTickLabels = cellstr(compose('%.3g', timeTicks));
[~, midTickIdx] = min(abs(timeTicks - midTimeSeconds));
timeTickLabels{midTickIdx} = 'M';

for iTarget = 1:numel(targetOrder)
    targetLabel = targetOrder{iTarget};
    if ~isfield(dRSA_PE, targetLabel)
        continue;
    end
    dRSAAll = {dRSA_PE.(targetLabel).position, dRSA_PE.(targetLabel).direction};
    figPEMatrices(iTarget) = figure( ...
        'Name', sprintf('PE dRSA matrices (%s)', targetLabel), ...
        'NumberTitle', 'off', ...
        'Visible', figVisibility);
    tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    rowAxes = gobjects(numel(dRSAAll), numel(peModelNames));
    for iRow = 1:numel(dRSAAll)
        for iModel = 1:numel(peModelNames)
            tileIdx = (iRow - 1) * numel(peModelNames) + iModel;
            nexttile(tileIdx);
            imagesc(timeSeconds, timeSeconds, dRSAAll{iRow}(:, :, iModel));
            set(gca, 'YDir', 'normal');
            axis image;
            if ~useCommonCbarLimits
                colorbar;
            end
            title(sprintf('%s | %s', local_pe_target_title(targetLabel), peModelNames{iModel}));
            % dRSA uses corr(modelRDM, neuralRDM): x = neural/data time, y = model time.
            xlabel(sprintf('Time in %s (s)', rowNames{iRow}));
            ylabel(sprintf('Time in %s (s)', peModelNames{iModel}));
            xticks(timeTicks);
            yticks(timeTicks);
            xticklabels(timeTickLabels);
            yticklabels(timeTickLabels);
            local_apply_offset_ticks(gca, timeTicks, midTickIdx);
            hold on;
            plot(timeSeconds, timeSeconds, 'w-', 'LineWidth', 1);
            hold off;
            if useCommonCbarLimits
                caxis(cbarLimits);
            end
            rowAxes(iRow, iModel) = gca;
        end
    end
    if useCommonCbarLimits
        for iRow = 1:numel(dRSAAll)
            cbarHandle = colorbar(rowAxes(iRow, end), 'Location', 'eastoutside');
            cbarHandle.Limits = cbarLimits;
        end
    end
end
end

function figHandle = local_plot_pe_paths(peDot1, peDot2, sampleRateHz, titleText, figVisibility)
% LOCAL_PLOT_PE_PATHS Plot PE trajectories (observed minus predicted) for both dots.
%
% Purpose:
%   Visualize trajectory-level PE in a compact 1x2 layout so dot1 and dot2
%   can be inspected side-by-side with a shared design language used in the
%   existing path plotting helpers.
%
% Inputs:
%   peDot1        : trials x 2 x time position PE for dot1.
%   peDot2        : trials x 2 x time position PE for dot2.
%   sampleRateHz  : sampling rate used to label colorbar in seconds.
%   titleText     : figure-level title.
%   figVisibility : 'on' or 'off'.
%
% Output:
%   figHandle     : figure handle for export.
%
% Data flow:
%   PE paths -> two plot_paths subplots -> single PE_paths figure handle.

figHandle = figure('Name', char(titleText), 'NumberTitle', 'off', 'Visible', figVisibility);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax1 = nexttile(1);
plot_paths(peDot1, 'ParentAxes', ax1, 'SavePlot', false, ...
    'Title', 'PE paths: dot1', 'SampleRateHz', sampleRateHz, ...
    'AxisLimits', []);
ax2 = nexttile(2);
plot_paths(peDot2, 'ParentAxes', ax2, 'SavePlot', false, ...
    'Title', 'PE paths: dot2', 'SampleRateHz', sampleRateHz, ...
    'AxisLimits', []);
sgtitle(figHandle, char(titleText));
end

function titleText = local_pe_target_title(targetLabel)
% LOCAL_PE_TARGET_TITLE Map PE target field names to readable panel titles.
switch targetLabel
    case 'neuralPE'
        titleText = 'Neural target: PE';
    case 'neuralPredicted'
        titleText = 'Neural target: predicted';
    case 'neuralObserved'
        titleText = 'Neural target: observed';
    otherwise
        titleText = ['Neural target: ' targetLabel];
end
end

function local_apply_offset_ticks(axHandle, timeTicks, midTickIdx)
% LOCAL_APPLY_OFFSET_TICKS Keep numeric ticks while offsetting crowded labels.
% Data flow: axes ticks -> hide selected labels -> re-add with manual offsets.
baseLabels = cellstr(compose('%.3g', timeTicks));
idxToOffset = unique([midTickIdx, numel(timeTicks)]);
xLabels = baseLabels;
yLabels = baseLabels;
xLabels(idxToOffset) = {''};
yLabels(idxToOffset) = {''};
xticklabels(axHandle, xLabels);
yticklabels(axHandle, yLabels);

xLim = xlim(axHandle);
yLim = ylim(axHandle);
xSpan = xLim(2) - xLim(1);
ySpan = yLim(2) - yLim(1);
xLabelYOffset = 0.07 * ySpan;
yLabelXOffset = 0.09 * xSpan;

for iTick = 1:numel(idxToOffset)
    currIdx = idxToOffset(iTick);
    if currIdx == midTickIdx
        tickText = 'M';
    else
        tickText = baseLabels{currIdx};
    end
    % Shift special x tick labels slightly downward.
    text(axHandle, timeTicks(currIdx), yLim(1) - xLabelYOffset, tickText, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'Clipping', 'off');
    % Shift special y tick labels slightly left.
    text(axHandle, xLim(1) - yLabelXOffset, timeTicks(currIdx), tickText, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'Clipping', 'off');
end
end

function local_prepare_figure_for_export(figHandle)
% LOCAL_PREPARE_FIGURE_FOR_EXPORT Maximize figure and keep on-screen size for print.
% Data flow: figure handle -> full-screen window state/position -> print-ready sizing.
if ~isgraphics(figHandle, 'figure')
    return;
end
set(figHandle, 'PaperPositionMode', 'auto');
try
    set(figHandle, 'WindowState', 'maximized');
catch
    % Fallback for environments where WindowState is unavailable.
    oldUnits = figHandle.Units;
    set(figHandle, 'Units', 'pixels');
    screenRect = get(groot, 'ScreenSize');
    set(figHandle, 'Position', screenRect);
    set(figHandle, 'Units', oldUnits);
end
drawnow;
end

function local_scale_figure_text(figHandle, scaleFactor)
% LOCAL_SCALE_FIGURE_TEXT Scale axes/colorbar/title/label text in a figure.
% Data flow: figure handle -> axes/colorbar handles -> uniformly enlarged text.
if ~isgraphics(figHandle, 'figure') || ~isscalar(scaleFactor) || scaleFactor <= 0
    return;
end
allAxes = findall(figHandle, 'Type', 'axes');
for iAx = 1:numel(allAxes)
    ax = allAxes(iAx);
    if isprop(ax, 'FontSize')
        ax.FontSize = ax.FontSize * scaleFactor;
    end
end
allText = findall(figHandle, 'Type', 'text');
for iText = 1:numel(allText)
    txt = allText(iText);
    if isprop(txt, 'FontSize')
        txt.FontSize = txt.FontSize * scaleFactor;
    end
end
end

function local_print_existing_match_summary( ...
        resultsFile, participantNumber, inputCondition, dRSAtypeToRun, ...
        shuffleDot2Trials, resampleSubsamples, resampleIterations, ...
        resampleSubsampleCount, resampleWithReplacement, ...
        usePredictedModels, plotSampleRateHz, cutPostDev, deviantOnset, existingMeta)
% LOCAL_PRINT_EXISTING_MATCH_SUMMARY Print parameter-by-parameter reuse summary.
% Data flow: current config + existing metadata -> one-line comparison per parameter.
fprintf('\nExisting results file found:\n  %s\n', resultsFile);

existingParticipant = [];
existingCondition = [];
existingType = [];
existingShuffle = [];
existingResample = [];
existingResampleIter = [];
existingResampleCount = [];
existingResampleReplacement = [];
existingUsePredicted = [];
existingFs = [];
existingCutPostDev = [];
existingDeviantOnset = [];

if isfield(existingMeta, 'repro')
    repro = existingMeta.repro;
    if isfield(repro, 'participantNumber')
        existingParticipant = repro.participantNumber;
    end
    if isfield(repro, 'inputCondition')
        existingCondition = repro.inputCondition;
    end
    if isfield(repro, 'shuffleDot2Trials')
        existingShuffle = repro.shuffleDot2Trials;
    end
    if isfield(repro, 'resampleSubsamples')
        existingResample = repro.resampleSubsamples;
    end
    if isfield(repro, 'resampleIterations')
        existingResampleIter = repro.resampleIterations;
    end
    if isfield(repro, 'resampleSubsampleCount')
        existingResampleCount = repro.resampleSubsampleCount;
    end
    if isfield(repro, 'resampleWithReplacement')
        existingResampleReplacement = repro.resampleWithReplacement;
    end
    if isfield(repro, 'nModels')
        existingUsePredicted = repro.nModels > 4;
    end
    if isfield(repro, 'cutPostDev')
        existingCutPostDev = repro.cutPostDev;
    end
    if isfield(repro, 'deviantOnset')
        existingDeviantOnset = repro.deviantOnset;
    end
    if isfield(repro, 'paramsCore')
        if isfield(repro.paramsCore, 'dRSAtype')
            existingType = repro.paramsCore.dRSAtype;
        end
        if isfield(repro.paramsCore, 'fs')
            existingFs = repro.paramsCore.fs;
        end
    end
end
if isempty(existingType) && isfield(existingMeta, 'params') && isfield(existingMeta.params, 'dRSAtype')
    existingType = existingMeta.params.dRSAtype;
end
if isempty(existingFs) && isfield(existingMeta, 'params') && isfield(existingMeta.params, 'fs')
    existingFs = existingMeta.params.fs;
end

local_print_match_line('participantNumber', participantNumber, existingParticipant);
local_print_match_line('inputCondition', inputCondition, existingCondition);
local_print_match_line('dRSAtypeToRun', dRSAtypeToRun, existingType);
local_print_match_line('shuffleDot2Trials', shuffleDot2Trials, existingShuffle);
local_print_match_line('resampleSubsamples', resampleSubsamples, existingResample);
local_print_match_line('resampleIterations', resampleIterations, existingResampleIter);
local_print_match_line('resampleSubsampleCount', resampleSubsampleCount, existingResampleCount);
local_print_match_line('resampleWithReplacement', resampleWithReplacement, existingResampleReplacement);
local_print_match_line('usePredictedModels', usePredictedModels, existingUsePredicted);
local_print_match_line('plotSampleRateHz', plotSampleRateHz, existingFs);
local_print_match_line('cutPostDev', cutPostDev, existingCutPostDev);
local_print_match_line('deviantOnset', deviantOnset, existingDeviantOnset);
end

function local_print_match_line(paramName, currentValue, existingValue)
% LOCAL_PRINT_MATCH_LINE Print one parameter comparison line.
isExistingMissing = isempty(existingValue);
isMatch = ~isExistingMissing && isequaln(currentValue, existingValue);
statusLabel = 'MISMATCH';
if isExistingMissing
    statusLabel = 'MISSING';
elseif isMatch
    statusLabel = 'MATCH';
end
fprintf('  - %-24s current=%-10s existing=%-10s [%s]\n', ...
    paramName, local_value_to_string(currentValue), ...
    local_value_to_string(existingValue), statusLabel);
end

function out = local_value_to_string(in)
% LOCAL_VALUE_TO_STRING Convert scalar/cell/string values to compact text.
if isempty(in)
    out = 'N/A';
elseif islogical(in) && isscalar(in)
    out = mat2str(in);
elseif isnumeric(in) && isscalar(in)
    out = num2str(in);
elseif ischar(in)
    out = in;
elseif isstring(in) && isscalar(in)
    out = char(in);
elseif iscell(in) && isscalar(in)
    out = local_value_to_string(in{1});
else
    out = mat2str(in);
end
end

function action = local_resolve_existing_results_action(existingResultsAction)
% LOCAL_RESOLVE_EXISTING_RESULTS_ACTION Resolve interactive vs configured reuse action.
%
% Inputs:
%   existingResultsAction : [] for interactive prompt, or scalar 1/2/3.
%
% Outputs:
%   action:
%     1 = load existing dataset + replot
%     2 = full rerun + overwrite dataset and plots
%     3 = just error out
if isempty(existingResultsAction)
    action = local_prompt_existing_results_action();
    return;
end
action = existingResultsAction;
end

function userChoice = local_prompt_existing_results_action()
% LOCAL_PROMPT_EXISTING_RESULTS_ACTION Ask whether to reuse, rerun, or stop.
% Output:
%   1 = load existing and replot
%   2 = rerun and overwrite all outputs
%   3 = just error out
userChoice = 0;
while ~ismember(userChoice, [1, 2, 3])
    userInput = input([ ...
        'Select action: (1) load existing dataset + replot (overwrite plots), ', ...
        '(2) rerun + overwrite dataset and plots, ', ...
        '(3) just error out: '], 's');
    userChoice = str2double(strtrim(userInput));
end
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

function modelNames = local_default_model_names(nModels)
% LOCAL_DEFAULT_MODEL_NAMES Provide fallback names if metadata is missing.
baseNames = {'position dot1', 'position dot2', 'direction dot1', 'direction dot2', ...
    'predicted position dot1', 'predicted position dot2', ...
    'predicted direction dot1', 'predicted direction dot2'};
if nModels <= numel(baseNames)
    modelNames = baseNames(1:nModels);
    return;
end
modelNames = cell(1, nModels);
for iModel = 1:nModels
    modelNames{iModel} = sprintf('model %d', iModel);
end
end
