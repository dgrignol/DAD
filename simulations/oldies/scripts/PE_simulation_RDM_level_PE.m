% PE_simulation.m
%
% Purpose:
%   Run counterfactual and prediction-error dRSA simulations on dot-motion
%   inputs by combining
%   position-only dRSA (as in simulations/scripts/PIPELINE_simulation_posOnly.m)
%   with direction-based models (as in simulations/scripts/PIPELINE_simulation.m).
%   The script loads center-relative dot-path data (dot 1 + dot 2, visual
%   degrees) for nondeviant and deviant conditions, prepares concatenated
%   inputs, and executes dRSA for both position and direction data. For the
%   deviant condition, it optionally adds "predicted deviant" position models
%   built from MovDot_SubXX_predicted.mat (deviant-only paths without
%   deviation). When running the deviant condition, it also adds predicted
%   direction models derived from the predicted paths, yielding 8 models
%   total. For deviant trials, a dedicated Prediction Error (PE) branch is
%   computed at the RDM level:
%       RDM_PE = RDM_observed_deviant - RDM_predicted_smooth
%   and is plotted against three neural targets (PE, predicted, observed)
%   using 2x4 matrix grids (row 1: position data, row 2: direction data).
%
% Example usage (from simulations/ in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   PE_simulation;
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts');
%   PE_simulation;
%
% Configuration (edit at the start of the script):
%   - participantNumber selects MovDot_SubXX.mat and the derived simulation
%     inputs (MovDot_SubXX_nondeviant/deviant.mat).
%   - inputConditions selects which conditions to run (default: both).
%   - shuffleDot2Trials decouples dot2 trial order from dot1 (diagnostic).
%   - centerRelativeRectSize overrides Cfg.rectSize when recentering.
%   - dRSA matrix figures are always generated in both colorbar modes:
%     shared [0 1] row-wise colorbars and per-subplot colorbars.
%   - textScaleFactor scales figure text (axes/titles/labels) before export.
%   - dRSAtypeToRun chooses which dRSA variant to execute ('PCR' or 'corr').
%   - if a matching results .mat already exists, the script prints a
%     parameter summary and prompts whether to reuse/replot or fully rerun.
%   - resampleSubsamples toggles trial subsampling across iterations to
%     reduce memory use; resampleIterations, resampleSubsampleCount, and
%     resampleWithReplacement control the resampling behavior.
%   - suppressDispText silences function-level status text (0 = show, 1 = suppress).
%   - cutPostDev keeps only the post-deviant portion of each trial; deviantOnset
%     (fraction of trial duration, 0-1) sets where the cut begins.
%   - output filenames include subject/condition, dRSA type, and (only when
%     enabled) shuffle/resampling tags.
%   - output directories are organized by subject/condition/dRSAtype with
%     subfolders for results, matrices, and diagonals; deviant predicted
%     models are stored under matrices/predicted (never "wannabe"), and
%     PE matrices are stored under matrices/PE. When
%     cutPostDev is enabled, outputs are stored under subXX/postDev.
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
% Outputs (organized):
%   - simulations/output/subXX/paths/paths_subXX_allConditions_all_conditions.png
%   - simulations/output/subXX/<condition>/<dRSAtype>/results/
%       dRSA_<subject>_<condition>_<dRSAtype>_[<shuffle>][_<resample>]_results.mat
%   - simulations/output/subXX/<condition>/<dRSAtype>/diagonal/
%       dRSA_<subject>_<condition>_<dRSAtype>_[<shuffle>][_<resample>]_diagonal.png
%   - simulations/output/subXX/<condition>/<dRSAtype>/matrices/
%       nondeviant: commCbar/ and sepCbar/
%       deviant: base/commCbar, base/sepCbar,
%                predicted/commCbar, predicted/sepCbar,
%                PE/commCbar, PE/sepCbar
%   - simulations/output/subXX/postDev/... (same structure, when cutPostDev = true)
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
participantNumber = 86;
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
dRSAtypeToRun = 'corr'; % 'PCR' or 'corr'
suppressDispText = 1; % 1 to suppress function-level status text, 0 to show.
validConditions = {'nondeviant', 'deviant'};
if ~all(ismember(inputConditions, validConditions))
    error('inputConditions must be drawn from: %s', strjoin(validConditions, ', '));
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
    'Title', sprintf('Dot paths (sub%02d, all conditions - Tot trial n: %02d)', participantNumber, size(pathsByCondition.nondeviant.dot1,1)));

%% Run counterfactual dRSA per condition
% Data flow: condition list -> per-condition load -> position/direction dRSA outputs.
for iCond = 1:numel(inputConditions)
    inputCondition = inputConditions{iCond};
    usePredictedModels = strcmp(inputCondition, 'deviant');
    reuseExistingResults = false;
    loadedResults = struct();
    cutFrame = [];
    peModelNames = {'PE position dot1', 'PE position dot2', ...
        'PE direction dot1', 'PE direction dot2'};
    peTargetTags = {'neuralPE', 'neuralPredicted', 'neuralObserved'};
    peRowNamesByTarget = { ...
        {'Data: PE position dot1', 'Data: PE direction dot1'}, ...
        {'Data: predicted position dot1', 'Data: predicted direction dot1'}, ...
        {'Data: observed position dot1', 'Data: observed direction dot1'}};
    dRSA_PE_neuralPE_position = [];
    dRSA_PE_neuralPE_direction = [];
    dRSA_PE_neuralPredicted_position = [];
    dRSA_PE_neuralPredicted_direction = [];
    dRSA_PE_neuralObserved_position = [];
    dRSA_PE_neuralObserved_direction = [];
    paramsPE = struct();
    peAutocorrBorder = [];
    peSubsamples = [];
    peUsedPostDevOnly = false;
    peCutFrame = [];

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
    % Data flow: model-group flag -> concrete matrix folders for export.
    if usePredictedModels
        baseMatricesDir = fullfile(matricesOutputDir, 'base');
        predictedMatricesDir = fullfile(matricesOutputDir, 'predicted');
        peMatricesDir = fullfile(matricesOutputDir, 'PE');
        matricesBaseCommDir = fullfile(baseMatricesDir, 'commCbar');
        matricesBaseSepDir = fullfile(baseMatricesDir, 'sepCbar');
        matricesPredCommDir = fullfile(predictedMatricesDir, 'commCbar');
        matricesPredSepDir = fullfile(predictedMatricesDir, 'sepCbar');
        matricesPECommDir = fullfile(peMatricesDir, 'commCbar');
        matricesPESepDir = fullfile(peMatricesDir, 'sepCbar');
    else
        matricesBaseCommDir = fullfile(matricesOutputDir, 'commCbar');
        matricesBaseSepDir = fullfile(matricesOutputDir, 'sepCbar');
        matricesPredCommDir = '';
        matricesPredSepDir = '';
        matricesPECommDir = '';
        matricesPESepDir = '';
    end
    local_ensure_dir(resultsOutputDir);
    local_ensure_dir(diagonalOutputDir);
    local_ensure_dir(matricesBaseCommDir);
    local_ensure_dir(matricesBaseSepDir);
    if usePredictedModels
        local_ensure_dir(matricesPredCommDir);
        local_ensure_dir(matricesPredSepDir);
        local_ensure_dir(matricesPECommDir);
        local_ensure_dir(matricesPESepDir);
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
        userChoice = local_prompt_existing_results_action();
        if userChoice == 1
            loadedResults = load(resultsFileToCheck);
            requiredVars = {'dRSA_position', 'dRSA_direction', ...
                'dRSA_diagonal_position', 'dRSA_diagonal_direction', 'params'};
            if usePredictedModels
                requiredVars = [requiredVars, { ...
                    'dRSA_PE_neuralPE_position', 'dRSA_PE_neuralPE_direction', ...
                    'dRSA_PE_neuralPredicted_position', 'dRSA_PE_neuralPredicted_direction', ...
                    'dRSA_PE_neuralObserved_position', 'dRSA_PE_neuralObserved_direction'}];
            end
            missingVars = requiredVars(~isfield(loadedResults, requiredVars));
            if isempty(missingVars)
                reuseExistingResults = true;
            else
                warning(['Existing file is missing required variables (%s). ' ...
                    'Falling back to full re-run.'], strjoin(missingVars, ', '));
            end
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
        if usePredictedModels
            dRSA_PE_neuralPE_position = loadedResults.dRSA_PE_neuralPE_position;
            dRSA_PE_neuralPE_direction = loadedResults.dRSA_PE_neuralPE_direction;
            dRSA_PE_neuralPredicted_position = loadedResults.dRSA_PE_neuralPredicted_position;
            dRSA_PE_neuralPredicted_direction = loadedResults.dRSA_PE_neuralPredicted_direction;
            dRSA_PE_neuralObserved_position = loadedResults.dRSA_PE_neuralObserved_position;
            dRSA_PE_neuralObserved_direction = loadedResults.dRSA_PE_neuralObserved_direction;
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peModelNames')
                peModelNames = loadedResults.repro.peModelNames;
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'paramsPE')
                paramsPE = loadedResults.repro.paramsPE;
            else
                paramsPE = paramsCore;
                paramsPE.modelNames = peModelNames;
                paramsPE.modelToTest = 1:numel(peModelNames);
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peAutocorrBorder')
                peAutocorrBorder = loadedResults.repro.peAutocorrBorder;
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peUsedPostDevOnly')
                peUsedPostDevOnly = loadedResults.repro.peUsedPostDevOnly;
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peCutFrame')
                peCutFrame = loadedResults.repro.peCutFrame;
            end
            if isfield(loadedResults, 'repro') && isfield(loadedResults.repro, 'peSubsamples')
                peSubsamples = loadedResults.repro.peSubsamples;
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
    % Data flow: dot positions -> frame-to-frame displacement -> angle -> unit vector features.
    dot1Dx = diff(dot1Trials(:, 1, :), 1, 3); % frame-to-frame x displacement (trials x 1 x time-1)
    dot1Dy = diff(dot1Trials(:, 2, :), 1, 3); % frame-to-frame y displacement (trials x 1 x time-1)
    dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx); % pad with first displacement so length matches original time
    dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy); % pad with first displacement so length matches original time
    dot1Angle = atan2(dot1Dy, dot1Dx); % per-timepoint direction angle in radians
    dot1Direction = cat(2, cos(dot1Angle), sin(dot1Angle)); % convert angles to unit vectors [cos; sin]
    modelDirectionDot1 = dRSA_concatenate(dot1Direction, [], plotConcat, ...
        'suppressDispText', suppressDispText); % flatten trial/time into dRSA-ready direction model

    dot2Dx = diff(dot2Trials(:, 1, :), 1, 3); % frame-to-frame x displacement for dot 2 (trials x 1 x time-1)
    dot2Dy = diff(dot2Trials(:, 2, :), 1, 3); % frame-to-frame y displacement for dot 2 (trials x 1 x time-1)
    dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx); % pad with first displacement so length matches original time
    dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy); % pad with first displacement so length matches original time
    dot2Angle = atan2(dot2Dy, dot2Dx); % per-timepoint direction angle in radians
    dot2Direction = cat(2, cos(dot2Angle), sin(dot2Angle)); % convert angles to unit vectors [cos; sin]
    modelDirectionDot2 = dRSA_concatenate(dot2Direction, [], plotConcat, ...
        'suppressDispText', suppressDispText); % flatten trial/time into dRSA-ready direction model
    dataDirection = modelDirectionDot1; % use dot 1 direction as the observed data series
    if usePredictedModels
        predictedDot1Dx = diff(dot1TrialsPredicted(:, 1, :), 1, 3);
        predictedDot1Dy = diff(dot1TrialsPredicted(:, 2, :), 1, 3);
        predictedDot1Dx = cat(3, predictedDot1Dx(:, :, 1), predictedDot1Dx);
        predictedDot1Dy = cat(3, predictedDot1Dy(:, :, 1), predictedDot1Dy);
        predictedDot1Angle = atan2(predictedDot1Dy, predictedDot1Dx);
        predictedDot1Direction = cat(2, cos(predictedDot1Angle), sin(predictedDot1Angle));
        modelDirectionDot1Predicted = dRSA_concatenate(predictedDot1Direction, [], plotConcat, ...
            'suppressDispText', suppressDispText);

        predictedDot2Dx = diff(dot2TrialsPredicted(:, 1, :), 1, 3);
        predictedDot2Dy = diff(dot2TrialsPredicted(:, 2, :), 1, 3);
        predictedDot2Dx = cat(3, predictedDot2Dx(:, :, 1), predictedDot2Dx);
        predictedDot2Dy = cat(3, predictedDot2Dy(:, :, 1), predictedDot2Dy);
        predictedDot2Angle = atan2(predictedDot2Dy, predictedDot2Dx);
        predictedDot2Direction = cat(2, cos(predictedDot2Angle), sin(predictedDot2Angle));
        modelDirectionDot2Predicted = dRSA_concatenate(predictedDot2Direction, [], plotConcat, ...
            'suppressDispText', suppressDispText);
    end

    %% Build trial-locked subsamples
    % Data flow: trial length -> trigger mask -> subsamples tensor for dRSA.
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
    % In case we do the PCR, it is better to calculate the border outside,
    % because otherwise we would need to recalculate it for each iteration.
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

    % For the PCR
    params.dRSAtype = dRSAtypeToRun; % 'PCR' or 'corr'
    params.PCR.AdditionalPCA = 1;
    params.PCR.Method = 'ExplainedVar';
    params.PCR.Methodfactor = 0.1;

    % Regress out autocorrelation
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

    %% Run dedicated Prediction Error dRSA branch (deviant only)
    % Data flow: observed/predicted paths -> observed/predicted RDMs ->
    % PE RDMs (observed - predicted) -> dRSA with three neural targets.
    if usePredictedModels
        % To avoid the zero-valued pre-deviance segment in PE, this branch
        % uses post-deviance windows unless the whole script already does.
        dot1TrialsPEObserved = dot1Trials;
        dot2TrialsPEObserved = dot2Trials;
        dot1TrialsPEPredicted = dot1TrialsPredicted;
        dot2TrialsPEPredicted = dot2TrialsPredicted;
        if ~cutPostDev
            [dot1TrialsPEObserved, peCutFrame] = local_cut_postdeviant(dot1TrialsPEObserved, deviantOnset);
            [dot2TrialsPEObserved, ~] = local_cut_postdeviant(dot2TrialsPEObserved, deviantOnset);
            [dot1TrialsPEPredicted, ~] = local_cut_postdeviant(dot1TrialsPEPredicted, deviantOnset);
            [dot2TrialsPEPredicted, ~] = local_cut_postdeviant(dot2TrialsPEPredicted, deviantOnset);
            peUsedPostDevOnly = true;
        else
            peCutFrame = cutFrame;
            peUsedPostDevOnly = true;
        end

        % Prepare observed/predicted feature-time series for PE-model RDMs.
        peObservedPositionDot1 = dRSA_concatenate(dot1TrialsPEObserved, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        peObservedPositionDot2 = dRSA_concatenate(dot2TrialsPEObserved, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        pePredictedPositionDot1 = dRSA_concatenate(dot1TrialsPEPredicted, [], plotConcat, ...
            'suppressDispText', suppressDispText);
        pePredictedPositionDot2 = dRSA_concatenate(dot2TrialsPEPredicted, [], plotConcat, ...
            'suppressDispText', suppressDispText);

        peObservedDirectionDot1 = dRSA_concatenate( ...
            local_paths_to_direction_features(dot1TrialsPEObserved), [], plotConcat, ...
            'suppressDispText', suppressDispText);
        peObservedDirectionDot2 = dRSA_concatenate( ...
            local_paths_to_direction_features(dot2TrialsPEObserved), [], plotConcat, ...
            'suppressDispText', suppressDispText);
        pePredictedDirectionDot1 = dRSA_concatenate( ...
            local_paths_to_direction_features(dot1TrialsPEPredicted), [], plotConcat, ...
            'suppressDispText', suppressDispText);
        pePredictedDirectionDot2 = dRSA_concatenate( ...
            local_paths_to_direction_features(dot2TrialsPEPredicted), [], plotConcat, ...
            'suppressDispText', suppressDispText);

        % Build PE-specific trial-locked subsamples (post-deviance windows).
        peTrialLen = size(dot1TrialsPEObserved, 3);
        peTotalTime = size(peObservedPositionDot1, 2);
        peMaskSubsampling = true(1, peTotalTime);
        if mod(peTotalTime, peTrialLen) ~= 0
            error('PE total time (%d) is not a multiple of PE trial length (%d).', ...
                peTotalTime, peTrialLen);
        end
        peTrialStarts = 1:peTrialLen:peTotalTime;
        peMaskTrigger = false(1, peTotalTime);
        peMaskTrigger(peTrialStarts) = true;

        peOpt.PreTrigger = 0;
        peOpt.PostTrigger = peTrialLen - 1;
        peOpt.spacing = 0;
        peOpt.nSubSamples = numel(peTrialStarts);
        peOpt.nIter = 1;
        peOpt.checkRepetition = 0;
        peSubsamples = dRSA_triggered_subsampling(peMaskSubsampling, peMaskTrigger, peOpt);

        if resampleSubsamples
            peBaseSubsamples = peSubsamples(:, :, 1);
            peBaseCount = size(peBaseSubsamples, 1);
            if resampleSubsampleCount > peBaseCount
                error('PE resampleSubsampleCount (%d) exceeds available subsamples (%d).', ...
                    resampleSubsampleCount, peBaseCount);
            end
            if resampleWithReplacement
                peResamplePercent = 100 * (resampleSubsampleCount / peBaseCount);
                peSubsamples = dRSA_resample_subsamples(peBaseSubsamples, ...
                    resampleIterations, peResamplePercent);
            else
                peSubsamples = zeros(resampleSubsampleCount, ...
                    size(peBaseSubsamples, 2), resampleIterations);
                for iIter = 1:resampleIterations
                    peDrawIdx = randperm(peBaseCount, resampleSubsampleCount);
                    peSubsamples(:, :, iIter) = peBaseSubsamples(peDrawIdx, :);
                end
            end
        end

        % PE dRSA settings: 4 PE models (pos dot1/2, dir dot1/2).
        paramsPE = params;
        paramsPE.modelNames = peModelNames;
        paramsPE.modelToTest = 1:numel(peModelNames);
        paramsPE.modelDistMeasure = {'euclidean', 'euclidean', 'cosine', 'cosine'};
        paramsPE.neuralDistMeasure = paramsPE.modelDistMeasure{1};
        paramsPE.nIter = size(peSubsamples, 3);
        if ~isfield(paramsPE.PCR, 'Method')
            paramsPE.PCR.Method = 'ExplainedVar';
        end
        if ~isfield(paramsPE.PCR, 'Methodfactor')
            paramsPE.PCR.Methodfactor = 0.1;
        end
        peAvgHalfWindowSamples = floor((peTrialLen - 1) / 2);
        paramsPE.AverageTime = peAvgHalfWindowSamples / paramsPE.fs;
        paramsPE.modeltoRegressout = cell(1, numel(peModelNames));
        for iModel = 1:numel(peModelNames)
            paramsPE.modeltoRegressout{iModel} = setdiff(1:numel(peModelNames), iModel);
        end

        peModelRDMByIter = cell(paramsPE.nIter, 1);
        peNeuralPEPositionRDMByIter = cell(paramsPE.nIter, 1);
        peNeuralPEDirectionRDMByIter = cell(paramsPE.nIter, 1);
        peNeuralPredictedPositionRDMByIter = cell(paramsPE.nIter, 1);
        peNeuralPredictedDirectionRDMByIter = cell(paramsPE.nIter, 1);
        peNeuralObservedPositionRDMByIter = cell(paramsPE.nIter, 1);
        peNeuralObservedDirectionRDMByIter = cell(paramsPE.nIter, 1);

        for iIter = 1:paramsPE.nIter
            CurrPESubsamples = peSubsamples(:, :, iIter);

            observedPositionDot1RDM = dRSA_computeRDM( ...
                peObservedPositionDot1, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{1});
            observedPositionDot2RDM = dRSA_computeRDM( ...
                peObservedPositionDot2, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{2});
            observedDirectionDot1RDM = dRSA_computeRDM( ...
                peObservedDirectionDot1, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{3});
            observedDirectionDot2RDM = dRSA_computeRDM( ...
                peObservedDirectionDot2, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{4});

            predictedPositionDot1RDM = dRSA_computeRDM( ...
                pePredictedPositionDot1, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{1});
            predictedPositionDot2RDM = dRSA_computeRDM( ...
                pePredictedPositionDot2, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{2});
            predictedDirectionDot1RDM = dRSA_computeRDM( ...
                pePredictedDirectionDot1, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{3});
            predictedDirectionDot2RDM = dRSA_computeRDM( ...
                pePredictedDirectionDot2, paramsPE, CurrPESubsamples, paramsPE.modelDistMeasure{4});

            peModelRDMByIter{iIter} = { ...
                observedPositionDot1RDM - predictedPositionDot1RDM, ...
                observedPositionDot2RDM - predictedPositionDot2RDM, ...
                observedDirectionDot1RDM - predictedDirectionDot1RDM, ...
                observedDirectionDot2RDM - predictedDirectionDot2RDM};

            peNeuralPEPositionRDMByIter{iIter} = peModelRDMByIter{iIter}{1};
            peNeuralPEDirectionRDMByIter{iIter} = peModelRDMByIter{iIter}{3};
            peNeuralPredictedPositionRDMByIter{iIter} = predictedPositionDot1RDM;
            peNeuralPredictedDirectionRDMByIter{iIter} = predictedDirectionDot1RDM;
            peNeuralObservedPositionRDMByIter{iIter} = observedPositionDot1RDM;
            peNeuralObservedDirectionRDMByIter{iIter} = observedDirectionDot1RDM;
        end

        % For PCR, compute autocorrelation borders directly on PE model RDMs.
        if strcmp(paramsPE.dRSAtype, 'PCR') && paramsPE.PCR.RegressAutocor
            peAutocorrBorder = local_compute_autocorr_border_from_model_rdms( ...
                peModelRDMByIter, paramsPE, suppressDispText);
        elseif strcmp(paramsPE.dRSAtype, 'PCR')
            peAutocorrBorder = nan(numel(peModelNames), 1);
        else
            peAutocorrBorder = [];
        end

        [dRSA_PE_neuralPE_position, dRSA_PE_neuralPE_direction] = ...
            local_run_drsa_from_rdm_sets( ...
            peModelRDMByIter, peNeuralPEPositionRDMByIter, ...
            peNeuralPEDirectionRDMByIter, paramsPE, peAutocorrBorder);
        [dRSA_PE_neuralPredicted_position, dRSA_PE_neuralPredicted_direction] = ...
            local_run_drsa_from_rdm_sets( ...
            peModelRDMByIter, peNeuralPredictedPositionRDMByIter, ...
            peNeuralPredictedDirectionRDMByIter, paramsPE, peAutocorrBorder);
        [dRSA_PE_neuralObserved_position, dRSA_PE_neuralObserved_direction] = ...
            local_run_drsa_from_rdm_sets( ...
            peModelRDMByIter, peNeuralObservedPositionRDMByIter, ...
            peNeuralObservedDirectionRDMByIter, paramsPE, peAutocorrBorder);
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

    %% Plot dedicated PE matrices (deviant-only)
    % Data flow: PE model dRSA + neural-target choice -> separate 2x4 grids.
    figMatricesPECommon = gobjects(0, 1);
    figMatricesPESeparate = gobjects(0, 1);
    if usePredictedModels
        pePositionTargets = { ...
            dRSA_PE_neuralPE_position, ...
            dRSA_PE_neuralPredicted_position, ...
            dRSA_PE_neuralObserved_position};
        peDirectionTargets = { ...
            dRSA_PE_neuralPE_direction, ...
            dRSA_PE_neuralPredicted_direction, ...
            dRSA_PE_neuralObserved_direction};
        peTimeSeconds = (0:size(pePositionTargets{1}, 1)-1) / paramsPE.fs;
        peModelIdx = 1:numel(peModelNames);
        for iPETarget = 1:numel(peTargetTags)
            peDRAAll = {pePositionTargets{iPETarget}, peDirectionTargets{iPETarget}};
            figCommonTmp = local_plot_dRSA_matrices( ...
                peDRAAll, peModelNames, peRowNamesByTarget{iPETarget}, peTimeSeconds, ...
                {'base'}, peModelIdx, [], true, cbarLimits, 'on');
            figSeparateTmp = local_plot_dRSA_matrices( ...
                peDRAAll, peModelNames, peRowNamesByTarget{iPETarget}, peTimeSeconds, ...
                {'base'}, peModelIdx, [], false, cbarLimits, 'on');
            figMatricesPECommon(iPETarget, 1) = figCommonTmp(1);
            figMatricesPESeparate(iPETarget, 1) = figSeparateTmp(1);
            set(figMatricesPECommon(iPETarget), 'Name', ...
                sprintf('dRSA PE matrices (%s, common cbar)', peTargetTags{iPETarget}));
            set(figMatricesPESeparate(iPETarget), 'Name', ...
                sprintf('dRSA PE matrices (%s, separate cbar)', peTargetTags{iPETarget}));
        end
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
        repro.paramsPE = paramsPE;
        repro.peAutocorrBorder = peAutocorrBorder;
        repro.peUsedPostDevOnly = peUsedPostDevOnly;
        repro.peCutFrame = peCutFrame;
        repro.peSubsamples = peSubsamples;

        save(resultsFile, ...
            'dRSA_position', 'dRSA_direction', ...
            'dRSA_PE_neuralPE_position', 'dRSA_PE_neuralPE_direction', ...
            'dRSA_PE_neuralPredicted_position', 'dRSA_PE_neuralPredicted_direction', ...
            'dRSA_PE_neuralObserved_position', 'dRSA_PE_neuralObserved_direction', ...
            'dRSA_diagonal_position', 'dRSA_diagonal_direction', ...
            'params', 'paramsDiagonal', 'repro');
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
    for iFigPrep = 1:numel(figMatricesPECommon)
        local_scale_figure_text(figMatricesPECommon(iFigPrep), textScaleFactor);
        local_prepare_figure_for_export(figMatricesPECommon(iFigPrep));
    end
    for iFigPrep = 1:numel(figMatricesPESeparate)
        local_scale_figure_text(figMatricesPESeparate(iFigPrep), textScaleFactor);
        local_prepare_figure_for_export(figMatricesPESeparate(iFigPrep));
    end
    % Keep lag plots at their native size/text for readability consistency.
    if iCond == 1
        local_scale_figure_text(figAllPaths, textScaleFactor);
        local_prepare_figure_for_export(figAllPaths);
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
        for iPETarget = 1:numel(peTargetTags)
            local_safe_print(figMatricesPECommon(iPETarget), ...
                fullfile(matricesPECommDir, [resultsBase '_matrices_PE_' ...
                peTargetTags{iPETarget} '_commCbar.png']));
            local_safe_print(figMatricesPESeparate(iPETarget), ...
                fullfile(matricesPESepDir, [resultsBase '_matrices_PE_' ...
                peTargetTags{iPETarget} '_sepCbar.png']));
        end
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
function directionFeatures = local_paths_to_direction_features(paths)
% LOCAL_PATHS_TO_DIRECTION_FEATURES Convert x/y paths into unit-vector directions.
%
% Usage example:
%   directionFeatures = local_paths_to_direction_features(dot1Trials);
%
% Inputs:
%   paths : trials x 2 x time center-relative x/y trajectories.
%
% Output:
%   directionFeatures : trials x 2 x time unit vectors [cos(theta); sin(theta)].
%
% Data flow:
%   x/y path samples -> frame-wise displacement -> angle -> unit vectors.
if isempty(paths)
    directionFeatures = [];
    return;
end
dx = diff(paths(:, 1, :), 1, 3);
dy = diff(paths(:, 2, :), 1, 3);
dx = cat(3, dx(:, :, 1), dx); % keep original trial length by front padding.
dy = cat(3, dy(:, :, 1), dy);
ang = atan2(dy, dx);
directionFeatures = cat(2, cos(ang), sin(ang));
end

function [dRSAPosition, dRSADirection] = local_run_drsa_from_rdm_sets( ...
        modelRDMByIter, neuralPositionRDMByIter, neuralDirectionRDMByIter, ...
        params, autocorrBorder)
% LOCAL_RUN_DRSA_FROM_RDM_SETS Run dRSA from precomputed RDM collections.
%
% Usage example:
%   [dPos, dDir] = local_run_drsa_from_rdm_sets(modelRDMByIter, ...
%       neuralPosRDMByIter, neuralDirRDMByIter, paramsPE, peAutocorrBorder);
%
% Inputs:
%   modelRDMByIter           : nIterx1 cell, each entry is 1xnModels cell
%                              of model RDMs in features x time format.
%   neuralPositionRDMByIter  : nIterx1 cell of neural RDMs (features x time).
%   neuralDirectionRDMByIter : nIterx1 cell of neural RDMs (features x time).
%   params                   : dRSA parameter struct (modelToTest, dRSAtype, PCR...).
%   autocorrBorder           : autocorrelation borders for PCR (or []).
%
% Outputs:
%   dRSAPosition : time x time x nModels dRSA matrix (position neural target).
%   dRSADirection: time x time x nModels dRSA matrix (direction neural target).
%
% Data flow:
%   RDM sets per iteration -> corr/PCR dRSA per iteration -> iteration mean.
nIter = numel(modelRDMByIter);
if nIter == 0
    dRSAPosition = [];
    dRSADirection = [];
    return;
end

firstModelIdx = params.modelToTest(1);
nTime = size(modelRDMByIter{1}{firstModelIdx}, 2);
nModels = numel(params.modelToTest);
dRSAPositionIter = zeros(nIter, nTime, nTime, nModels);
dRSADirectionIter = zeros(nIter, nTime, nTime, nModels);

for iIter = 1:nIter
    currModelRDM = modelRDMByIter{iIter};
    currNeuralPosRDM = neuralPositionRDMByIter{iIter};
    currNeuralDirRDM = neuralDirectionRDMByIter{iIter};
    if strcmp(params.dRSAtype, 'corr')
        dRSAPositionIter(iIter, :, :, :) = dRSA_corr(currModelRDM, currNeuralPosRDM, params);
        dRSADirectionIter(iIter, :, :, :) = dRSA_corr(currModelRDM, currNeuralDirRDM, params);
    elseif strcmp(params.dRSAtype, 'PCR')
        dRSAPositionIter(iIter, :, :, :) = dRSA_PCR(currModelRDM, currNeuralPosRDM, ...
            autocorrBorder, params);
        dRSADirectionIter(iIter, :, :, :) = dRSA_PCR(currModelRDM, currNeuralDirRDM, ...
            autocorrBorder, params);
    else
        error('Unsupported dRSAtype for PE RDM branch: %s', params.dRSAtype);
    end
end

dRSAPosition = mean(dRSAPositionIter, 1);
dRSAPosition = reshape(dRSAPosition, size(dRSAPosition, 2), ...
    size(dRSAPosition, 3), size(dRSAPosition, 4));
dRSADirection = mean(dRSADirectionIter, 1);
dRSADirection = reshape(dRSADirection, size(dRSADirection, 2), ...
    size(dRSADirection, 3), size(dRSADirection, 4));
end

function border = local_compute_autocorr_border_from_model_rdms( ...
        modelRDMByIter, params, suppressDispText)
% LOCAL_COMPUTE_AUTOCORR_BORDER_FROM_MODEL_RDMS Estimate PCR borders from RDMs.
%
% Usage example:
%   border = local_compute_autocorr_border_from_model_rdms( ...
%       peModelRDMByIter, paramsPE, 1);
%
% Inputs:
%   modelRDMByIter  : nIterx1 cell, each entry is 1xnModels cell of
%                     model RDMs in features x time format.
%   params          : dRSA params with .Var, .AverageTime, .fs, .modelToTest.
%   suppressDispText: 0/1 print suppression flag.
%
% Output:
%   border          : nModels x 1 autocorrelation border in samples.
%
% Data flow:
%   per-model RDM autocorr (corr) -> diagonal average -> thresholded border.
nIter = numel(modelRDMByIter);
if nIter == 0
    border = [];
    return;
end

modelToTest = params.modelToTest;
nModels = numel(modelToTest);
nTime = size(modelRDMByIter{1}{modelToTest(1)}, 2);
autocorrIter = zeros(nIter, nTime, nTime, nModels);
for iIter = 1:nIter
    for iModelPos = 1:nModels
        iModel = modelToTest(iModelPos);
        currModelRDM = modelRDMByIter{iIter}{iModel};
        autocorrIter(iIter, :, :, iModelPos) = corr(currModelRDM, currModelRDM);
    end
end
autocorrMean = squeeze(mean(autocorrIter, 1));
if nModels == 1
    autocorrMean = reshape(autocorrMean, nTime, nTime, 1);
end

paramsForAverage = params;
paramsForAverage.modelToTest = 1:nModels;
averagedAutocorr = dRSA_average(autocorrMean, paramsForAverage);
lag0Idx = (size(averagedAutocorr, 2) + 1) / 2;
regValue = sqrt(params.Var);
border = NaN(nModels, 1);

for iModelPos = 1:nModels
    peak = max(averagedAutocorr(iModelPos, :));
    leftSide = squeeze(averagedAutocorr(iModelPos, 1:lag0Idx));
    rightSide = squeeze(averagedAutocorr(iModelPos, lag0Idx:end));
    leftBorder = find(leftSide < regValue * peak, 1, 'last');
    rightBorder = length(rightSide) - find(rightSide < regValue * peak, 1, 'first');
    if isempty(leftBorder) || isempty(rightBorder)
        border(iModelPos) = NaN;
        if ~suppressDispText
            warning('PE autocorr border not found for model index %d.', modelToTest(iModelPos));
        end
        continue;
    end
    border(iModelPos) = ceil(nanmean([leftBorder rightBorder]));
end
end

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

function userChoice = local_prompt_existing_results_action()
% LOCAL_PROMPT_EXISTING_RESULTS_ACTION Ask whether to reuse or rerun.
% Output: 1 = load existing and replot, 2 = rerun and overwrite all.
userChoice = 0;
while ~ismember(userChoice, [1, 2])
    userInput = input([ ...
        'Select action: (1) load existing dataset + replot (overwrite plots), ', ...
        '(2) rerun + overwrite dataset and plots: '], 's');
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
