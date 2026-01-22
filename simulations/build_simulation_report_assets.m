% build_simulation_report_assets.m
%
% Purpose:
%   Compute simulation results and build grouped figures for the report.
%   Saves all outputs in simulations/output so the Report Generator script
%   can assemble a PDF without recomputing dRSA.
%
% Example usage (from repo root in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   run('simulations/build_simulation_report_assets.m');
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/build_simulation_report_assets.m');
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat (raw dot-motion paths)
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat (generated)
%
% Outputs:
%   - simulations/output/subXX/report_assets_subXX.mat
%   - simulations/output/subXX/report_figures/*.png
%
% Key assumptions:
%   - Dot paths are trials × 2 × time in visual degrees.
%   - Direction models use unit vectors derived from frame-to-frame angles.
%   - dRSA functions in simulations/functions are available on the path.

%% Resolve paths and add dependencies
% Data flow: script location -> simulations dir -> repo root -> addpath.
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptPath = which('build_simulation_report_assets.m');
end
simDir = fileparts(scriptPath);
repoRoot = fileparts(simDir);
addpath(simDir);
addpath(fullfile(simDir, 'functions'));
addpath(fullfile(simDir, 'debug'));
addpath(repoRoot);

%% Report configuration
% Data flow: configuration -> input selection -> report assets.
% Options:
%   - useSavedResults: true to reuse cached dRSA results if available.
%   - saveResults: true to save computed results for future re-plotting.
%   - includePCR: true to compute and plot PCR-based dRSA outputs.
useSavedResults = false;
saveResults = true;
includePCR = true;
participantNumber = 98;
conditions = {'nondeviant', 'deviant'};
% Use a single condition for the 2x4 dRSA matrices to avoid a 2x8 layout.
drsaReportCondition = conditions{1};
simulationInputDir = fullfile(simDir, 'input');
subjectLabel = sprintf('sub%02d', participantNumber);
outputDir = fullfile(simDir, 'output', subjectLabel);
figureOutputDir = fullfile(outputDir, 'report_figures');
if ~exist(figureOutputDir, 'dir')
    mkdir(figureOutputDir);
end
assetsFile = fullfile(outputDir, sprintf('report_assets_%s.mat', subjectLabel));

% Figure typography settings (increase for report readability).
fontSizes = struct();
fontSizes.title = 36;
fontSizes.label = 32;
fontSizes.tick = 28;
fontSizes.legend = 28;
fontSizes.colorbar = 24;

params = struct();
params.nIter = 1;
params.fs = 120;
params.modelToTest = [1 2 3 4];
params.Var = 0.1;
params.modelDistMeasure = {'euclidean', 'euclidean', 'cosine', 'cosine'};
params.neuralDistMeasure = 'euclidean';
params.dRSAtype = 'corr';
params.modeltoRegressout = {[2 3 4] [1 3 4] [1 2 4] [1 2 3]};

paramsPCR = params;
paramsPCR.dRSAtype = 'PCR';

%% Suppress figure windows while saving report assets
% Data flow: default figure visibility -> hidden figure generation.
defaultFigureVisible = get(0, 'DefaultFigureVisible');
cleanupVisibility = onCleanup(@() set(0, 'DefaultFigureVisible', defaultFigureVisible));
set(0, 'DefaultFigureVisible', 'off');

%% Ensure condition-specific inputs exist
% Data flow: raw MovDot input -> condition-specific dot paths.
inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
if ~exist(simulationInputDir, 'dir')
    mkdir(simulationInputDir);
end
for iCond = 1:numel(conditions)
    conditionLabel = conditions{iCond};
    conditionFile = fullfile(simulationInputDir, ...
        sprintf('MovDot_Sub%02d_%s.mat', participantNumber, conditionLabel));
    if ~exist(conditionFile, 'file')
        build_movdot_simulation_inputs(inputFile, ...
            'OutputDir', simulationInputDir, ...
            'MoveDotScript', fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m'));
        break;
    end
end

%% Compute or load cached results
% Data flow: cached MAT file -> results struct or fresh computations.
if useSavedResults
    if ~exist(assetsFile, 'file')
        error(['Missing cached results file: %s. Set useSavedResults = false ', ...
            'to compute results the first time.'], assetsFile);
    end
    loaded = load(assetsFile, 'results', 'participantNumber', 'conditions', 'includePCR');
    results = loaded.results;
    % Enforce exact match when using cached results.
    if loaded.participantNumber ~= participantNumber || ~isequal(loaded.conditions, conditions)
        error(['Cached results were computed for different settings. ', ...
            'Set useSavedResults = false to recompute with current settings.']);
    end
    if ~isfield(loaded, 'includePCR') || loaded.includePCR ~= includePCR
        error(['Cached results were computed with a different includePCR setting. ', ...
            'Set useSavedResults = false to recompute with current settings.']);
    end
else
    results = compute_results(simulationInputDir, participantNumber, conditions, params, paramsPCR, includePCR);
end

if saveResults
    save(assetsFile, 'results', 'participantNumber', 'conditions', 'includePCR');
end

%% Build grouped figures and save to disk
% Data flow: results -> grouped figures -> PNGs in output directory.
fig = plot_grouped_paths(results, conditions, fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('dot_paths_grouped_%s.png', subjectLabel)));

rectSize = load_rect_size(inputFile);
fig = plot_condition_position_distribution(results, conditions, fontSizes, subjectLabel, rectSize);
save_figure(fig, fullfile(figureOutputDir, ...
    sprintf('dot_paths_position_distribution_%s.png', subjectLabel)));

fig = plot_grouped_distances(results, conditions, fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('distance_matrices_grouped_%s.png', subjectLabel)));

fig = plot_grouped_drsa_combined(results, drsaReportCondition, 'corr', fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('drsa_matrices_corr_%s.png', subjectLabel)));

if includePCR
    fig = plot_grouped_drsa_combined(results, drsaReportCondition, 'pcr', fontSizes, subjectLabel);
    save_figure(fig, fullfile(figureOutputDir, sprintf('drsa_matrices_pcr_%s.png', subjectLabel)));
end

fig = plot_all_lagged(results, conditions, fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('lagged_drsa_%s.png', subjectLabel)));

%% Save report assets metadata
% Data flow: parameters + paths -> MAT file for report generator.
save(assetsFile, 'results', 'participantNumber', 'conditions', 'params', ...
    'paramsPCR', 'fontSizes', 'figureOutputDir', 'includePCR');

%% Restore figure visibility (ensures plots are re-enabled after the script)
% Data flow: onCleanup handle -> restore DefaultFigureVisible.
clear cleanupVisibility

%% Local helper functions
% Data flow: encapsulate loading, dRSA computation, and plotting utilities.

function [dot1Paths, dot2Paths] = load_condition_paths(simInputDir, participantNumber, conditionLabel)
% load_condition_paths
%
% Purpose:
%   Load condition-specific dot path arrays for a participant.
%
% Inputs:
%   simInputDir       : directory containing MovDot_SubXX_condition.mat
%   participantNumber : numeric participant ID
%   conditionLabel    : 'nondeviant' or 'deviant'
%
% Outputs:
%   dot1Paths : trials × 2 × time array for dot 1
%   dot2Paths : trials × 2 × time array for dot 2
%
% Key assumptions:
%   - The MAT-file contains dot1GreenPaths and dot2YellowPaths variables.

    inputFile = fullfile(simInputDir, ...
        sprintf('MovDot_Sub%02d_%s.mat', participantNumber, lower(conditionLabel)));
    if ~exist(inputFile, 'file')
        error('Condition input file not found: %s', inputFile);
    end
    simData = load(inputFile);
    dot1Paths = simData.dot1GreenPaths;
    dot2Paths = simData.dot2YellowPaths;
end

function [dRSA, dRSA_diagonal, diagTimeVec, modelNames] = run_drsa(dot1Paths, dot2Paths, paramsIn)
% run_drsa
%
% Purpose:
%   Compute dRSA matrices and their averaged diagonals for position and
%   direction models using dot 1 position as data.
%
% Inputs:
%   dot1Paths : trials × 2 × time array for dot 1
%   dot2Paths : trials × 2 × time array for dot 2
%   paramsIn  : struct of dRSA parameters (model distances, fs, etc.)
%
% Outputs:
%   dRSA          : time × time × model dRSA matrices
%   dRSA_diagonal : model × time diagonal averages
%   diagTimeVec   : lag vector in samples for plotting
%   modelNames    : 1x4 cell array of model names
%
% Key assumptions:
%   - dot paths share the same number of trials and time samples.
%   - Trial length divides the concatenated time series evenly.

    % Prepare concatenated data and models.
    dataPosition = dRSA_concatenate(dot1Paths);

    % Position models (per dot).
    modelPositionDot1 = dataPosition;
    modelPositionDot2 = dRSA_concatenate(dot2Paths);

    % Direction models (per dot).
    % Data flow: positions -> displacement -> angle -> unit vector features.
    dot1Dx = diff(dot1Paths(:, 1, :), 1, 3);
    dot1Dy = diff(dot1Paths(:, 2, :), 1, 3);
    dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx);
    dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy);
    dot1Angle = atan2(dot1Dy, dot1Dx);
    dot1Direction = cat(2, cos(dot1Angle), sin(dot1Angle));
    modelDirectionDot1 = dRSA_concatenate(dot1Direction);

    dot2Dx = diff(dot2Paths(:, 1, :), 1, 3);
    dot2Dy = diff(dot2Paths(:, 2, :), 1, 3);
    dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx);
    dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy);
    dot2Angle = atan2(dot2Dy, dot2Dx);
    dot2Direction = cat(2, cos(dot2Angle), sin(dot2Angle));
    modelDirectionDot2 = dRSA_concatenate(dot2Direction);

    model = {modelPositionDot1, modelPositionDot2, modelDirectionDot1, modelDirectionDot2};
    modelNames = {'position dot1', 'position dot2', 'direction dot1', 'direction dot2'};
    Y = dataPosition;

    % Build trigger-aligned subsamples (trial-locked).
    trialLen = size(dot1Paths, 3);
    totalTime = size(dataPosition, 2);
    if mod(totalTime, trialLen) ~= 0
        error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
    end
    trialStarts = 1:trialLen:totalTime;
    maskSubsampling = true(1, totalTime);
    maskTrigger = false(1, totalTime);
    maskTrigger(trialStarts) = true;

    opt.PreTrigger = 0;
    opt.PostTrigger = trialLen - 1;
    opt.spacing = 0;
    opt.nSubSamples = numel(trialStarts);
    opt.nIter = 1;
    opt.checkRepetition = 0;

    subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);

    % Copy params and set average window based on trial length.
    params = paramsIn;
    avgHalfWindowSamples = floor((trialLen - 1) / 2);
    params.AverageTime = avgHalfWindowSamples / params.fs;

    % Compute autocorrelation borders and dRSA matrices.
    Autocorrborder = dRSA_border(model, subsamples, params);
    dRSA_Iter = [];
    for iIter = 1:params.nIter
        CurrSubsamples = subsamples(:, :, iIter);
        dRSAma = dRSA_coreFunction(Y, model, params, ...
            'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);
        dRSA_Iter(iIter, :, :, :) = dRSAma;
    end

    dRSA = mean(dRSA_Iter, 1);
    dRSA = reshape(dRSA, size(dRSA, 2), size(dRSA, 3), size(dRSA, 4));

    % Average across the diagonal for lagged plots.
    diagParams = params;
    diagParams.AverageTime = 2;
    diagParams.fs = params.fs;
    dRSA_diagonal = dRSA_average(dRSA, diagParams);
    diagTimeVec = (1:size(dRSA_diagonal, 2)) - ceil(size(dRSA_diagonal, 2) / 2);
end

function results = compute_results(simulationInputDir, participantNumber, conditions, params, paramsPCR, includePCR)
% compute_results
%
% Purpose:
%   Compute per-condition dRSA outputs for corr and optionally PCR, and cache results.
%
% Inputs:
%   simulationInputDir : directory with condition-specific dot inputs
%   participantNumber : numeric participant ID
%   conditions        : cell array of condition labels
%   params            : dRSA parameters for corr
%   paramsPCR         : dRSA parameters for PCR
%   includePCR        : true to compute PCR outputs
%
% Outputs:
%   results : struct containing dot paths and dRSA outputs per condition

    results = struct();
    for iCond = 1:numel(conditions)
        conditionLabel = conditions{iCond};
        [dot1Paths, dot2Paths] = load_condition_paths(simulationInputDir, ...
            participantNumber, conditionLabel);

        % Compute dRSA matrices (corr; position data from dot 1).
        [dRSA_position, dRSA_diagonal_position, diagTimeVec_position, modelNames] = ...
            run_drsa(dot1Paths, dot2Paths, params);

        % Compute dRSA matrices (corr; direction data from dot 1).
        [dRSA_direction, dRSA_diagonal_direction, diagTimeVec_direction] = ...
            run_drsa_direction(dot1Paths, dot2Paths, params);

        % Compute dRSA matrices (PCR; position and direction data from dot 1).
        dRSA_position_pcr = [];
        dRSA_diagonal_position_pcr = [];
        diagTimeVec_position_pcr = [];
        dRSA_direction_pcr = [];
        dRSA_diagonal_direction_pcr = [];
        diagTimeVec_direction_pcr = [];
        if includePCR
            [dRSA_position_pcr, dRSA_diagonal_position_pcr, diagTimeVec_position_pcr] = ...
                run_drsa(dot1Paths, dot2Paths, paramsPCR);
            [dRSA_direction_pcr, dRSA_diagonal_direction_pcr, diagTimeVec_direction_pcr] = ...
                run_drsa_direction(dot1Paths, dot2Paths, paramsPCR);
        end

        results.(conditionLabel).dot1Paths = dot1Paths;
        results.(conditionLabel).dot2Paths = dot2Paths;
        results.(conditionLabel).modelNames = modelNames;
        results.(conditionLabel).dRSA_position = dRSA_position;
        results.(conditionLabel).dRSA_direction = dRSA_direction;
        results.(conditionLabel).dRSA_diagonal_position = dRSA_diagonal_position;
        results.(conditionLabel).diagTimeVec_position = diagTimeVec_position;
        results.(conditionLabel).dRSA_diagonal_direction = dRSA_diagonal_direction;
        results.(conditionLabel).diagTimeVec_direction = diagTimeVec_direction;
        results.(conditionLabel).dRSA_position_pcr = dRSA_position_pcr;
        results.(conditionLabel).dRSA_direction_pcr = dRSA_direction_pcr;
        results.(conditionLabel).dRSA_diagonal_position_pcr = dRSA_diagonal_position_pcr;
        results.(conditionLabel).diagTimeVec_position_pcr = diagTimeVec_position_pcr;
        results.(conditionLabel).dRSA_diagonal_direction_pcr = dRSA_diagonal_direction_pcr;
        results.(conditionLabel).diagTimeVec_direction_pcr = diagTimeVec_direction_pcr;
    end
end

function figHandle = plot_condition_position_distribution(results, conditions, fontSizes, subjectLabel, rectSize)
% plot_condition_position_distribution
%
% Purpose:
%   Plot mean distance-from-center profiles for dot 1 and dot 2, averaged
%   across conditions, to mirror plot_position_distribution in the report.
%
% Inputs:
%   results    : struct with per-condition dot paths
%   conditions : cell array of condition labels
%   fontSizes  : struct of typography sizes
%   rectSize   : 1x2 [width height] in visual degrees (optional)
%
% Outputs:
%   figHandle : figure handle for the position distribution plot

    nCond = numel(conditions);
    meanSumDot1 = [];
    meanSumDot2 = [];
    for iCond = 1:nCond
        conditionLabel = conditions{iCond};
        dot1 = results.(conditionLabel).dot1Paths;
        dot2 = results.(conditionLabel).dot2Paths;
        meanDot1 = compute_mean_distance_from_center(dot1, rectSize, [0 0]);
        meanDot2 = compute_mean_distance_from_center(dot2, rectSize, [0 0]);
        if iCond == 1
            meanSumDot1 = meanDot1;
            meanSumDot2 = meanDot2;
        else
            meanSumDot1 = meanSumDot1 + meanDot1;
            meanSumDot2 = meanSumDot2 + meanDot2;
        end
    end
    meanDot1 = meanSumDot1 ./ nCond;
    meanDot2 = meanSumDot2 ./ nCond;

    figHandle = figure('Name', 'Mean distance from center', ...
        'NumberTitle', 'off', 'Visible', 'off');
    figHandle.Position = [100 100 1800 700];
    sgtitle(sprintf('Mean distance from center (mean across conditions, %s)', subjectLabel), ...
        'FontSize', fontSizes.title + 2);

    subplot(1, 2, 1);
    plot(meanDot1, 'LineWidth', 1.8);
    grid on;
    title('Dot 1 distance from center', 'FontSize', fontSizes.title);
    xlabel('Time (samples)', 'FontSize', fontSizes.label);
    ylabel('Mean distance from center (visual degrees)', 'FontSize', fontSizes.label);
    set(gca, 'FontSize', fontSizes.tick);

    subplot(1, 2, 2);
    plot(meanDot2, 'LineWidth', 1.8);
    grid on;
    title('Dot 2 distance from center', 'FontSize', fontSizes.title);
    xlabel('Time (samples)', 'FontSize', fontSizes.label);
    ylabel('Mean distance from center (visual degrees)', 'FontSize', fontSizes.label);
    set(gca, 'FontSize', fontSizes.tick);
end

function figHandle = plot_grouped_paths(results, conditions, fontSizes, subjectLabel)
% plot_grouped_paths
%
% Purpose:
%   Plot dot paths grouped by dot (rows) and condition (columns).
%
% Inputs:
%   results    : struct with per-condition dot paths
%   conditions : cell array of condition labels
%   fontSizes  : struct of typography sizes
%
% Outputs:
%   figHandle : figure handle for the grouped dot-path plot

    nCond = numel(conditions);
    axisLimits = compute_paths_limits(results, conditions);
    figHandle = figure('Name', 'Dot paths grouped', 'NumberTitle', 'off', 'Visible', 'off');
    figHandle.Position = [100 100 2200 900];
    sgtitle(sprintf('Dot paths (%s)', subjectLabel), 'FontSize', fontSizes.title + 2);
    lastAx = [];
    for iCond = 1:nCond
        conditionLabel = conditions{iCond};
        dot1 = results.(conditionLabel).dot1Paths;
        dot2 = results.(conditionLabel).dot2Paths;

        subplot(2, nCond, iCond);
        plot_paths_axes(gca, dot1, sprintf('Dot 1 paths (%s)', conditionLabel), ...
            fontSizes, axisLimits, false);

        subplot(2, nCond, nCond + iCond);
        plot_paths_axes(gca, dot2, sprintf('Dot 2 paths (%s)', conditionLabel), ...
            fontSizes, axisLimits, false);
        lastAx = gca;
    end
    if ~isempty(lastAx)
        add_shared_colorbar(lastAx, fontSizes);
    end
end

function figHandle = plot_grouped_distances(results, conditions, fontSizes, subjectLabel)
% plot_grouped_distances
%
% Purpose:
%   Plot position-time distance matrices grouped by dot and condition.
%
% Inputs:
%   results    : struct with per-condition dot paths
%   conditions : cell array of condition labels
%   fontSizes  : struct of typography sizes
%
% Outputs:
%   figHandle : figure handle for the grouped distance matrix plot

    nCond = numel(conditions);
    figHandle = figure('Name', 'Distance matrices grouped', 'NumberTitle', 'off', 'Visible', 'off');
    figHandle.Position = [100 100 2200 900];
    sgtitle(sprintf('Position-time distance matrices (%s)', subjectLabel), ...
        'FontSize', fontSizes.title + 2);
    for iCond = 1:nCond
        conditionLabel = conditions{iCond};
        dot1 = results.(conditionLabel).dot1Paths;
        dot2 = results.(conditionLabel).dot2Paths;

        subplot(2, nCond, iCond);
        plot_distance_axes(gca, dot1, sprintf('Dot 1 distance (%s)', conditionLabel), fontSizes);

        subplot(2, nCond, nCond + iCond);
        plot_distance_axes(gca, dot2, sprintf('Dot 2 distance (%s)', conditionLabel), fontSizes);
    end
end

function figHandle = plot_grouped_drsa_combined(results, conditionLabel, drsaType, fontSizes, subjectLabel)
% plot_grouped_drsa_combined
%
% Purpose:
%   Plot dRSA matrices in a single figure with two rows:
%   top row = dot 1 position data, bottom row = dot 1 direction data.
%
% Inputs:
%   results        : struct with per-condition dRSA outputs
%   conditionLabel : condition label to plot (e.g., 'nondeviant')
%   drsaType       : 'corr' or 'pcr' to select the dRSA variant
%   fontSizes      : struct of typography sizes
%
% Outputs:
%   figHandle : figure handle for the combined dRSA plot

    nModels = 4;
    figHandle = figure('Name', 'dRSA matrices combined', ...
        'NumberTitle', 'off', 'Visible', 'off');
    figHandle.Position = [100 100 2200 900];
    sgtitle(sprintf('dRSA matrices (%s, %s; %s)', ...
        upper(drsaType), subjectLabel, format_condition_label(conditionLabel)), ...
        'FontSize', fontSizes.title + 2);

    if ~isfield(results, conditionLabel)
        error('Condition label not found in results: %s', conditionLabel);
    end
    rowDataLabels = {'data: position dot1', 'data: direction dot1'};
    modelNames = results.(conditionLabel).modelNames;
    for iRow = 1:2
        if iRow == 1
            dRSA = select_drsa(results.(conditionLabel), drsaType, false);
        else
            dRSA = select_drsa(results.(conditionLabel), drsaType, true);
        end
        for iModel = 1:nModels
            subplot(2, nModels, (iRow - 1) * nModels + iModel);
            ax = gca;
            plot_drsa_axes(ax, dRSA, iModel, modelNames{iModel}, false, false, fontSizes);
            if iRow == 2
                xlabel(ax, sprintf('Time in %s (samples)', modelNames{iModel}), ...
                    'FontSize', fontSizes.label);
            end
            if iModel == 1
                ylabel(ax, sprintf('Time in %s (samples)', rowDataLabels{iRow}), ...
                    'FontSize', fontSizes.label);
            end
        end
    end
end

function figHandle = plot_all_lagged(results, conditions, fontSizes, subjectLabel)
% plot_all_lagged
%
% Purpose:
%   Plot all lagged dRSA curves in one figure with two panels:
%   - position data (dot 1)
%   - direction data (dot 1)
%
% Inputs:
%   results    : struct with per-condition dRSA diagonals
%   conditions : cell array of condition labels
%   fontSizes  : struct of typography sizes
%
% Outputs:
%   figHandle : figure handle for the combined lagged plots

    figHandle = figure('Name', 'Lagged dRSA combined', ...
        'NumberTitle', 'off', 'Visible', 'off');
    figHandle.Position = [100 100 2000 700];
    sgtitle(sprintf('Lagged dRSA (%s)', subjectLabel), 'FontSize', fontSizes.title + 2);

    % Panel 1: position data (dot 1).
    subplot(1, 2, 1);
    plot_lagged_within_combined(gca, results, conditions, false, fontSizes);

    % Panel 2: direction data (dot 1).
    subplot(1, 2, 2);
    plot_lagged_within_combined(gca, results, conditions, true, fontSizes);
end

function plot_paths_axes(ax, paths, titleText, fontSizes, axisLimits, showColorbar)
% plot_paths_axes
%
% Purpose:
%   Plot dot paths into the provided axes with time-graded colors.
%
% Inputs:
%   ax        : target axes handle
%   paths     : trials × 2 × time array of dot positions
%   titleText : subplot title string
%   fontSizes : struct of typography sizes
%   axisLimits : 1x4 [xmin xmax ymin ymax] for consistent axes
%   showColorbar: true to add a colorbar on this axes

    nTrials = size(paths, 1);
    nTime = size(paths, 3);
    xVals = reshape(paths(:, 1, :), [], 1);
    yVals = reshape(paths(:, 2, :), [], 1);
    timeIdx = repmat(1:nTime, nTrials, 1);
    timeIdx = timeIdx(:);
    cmap = parula(nTime);
    pointColors = cmap(timeIdx, :);

    scatter(ax, xVals, yVals, 6, pointColors, 'filled');
    axis(ax, 'equal');
    if ~isempty(axisLimits)
        axis(ax, axisLimits);
    end
    colormap(ax, cmap);
    caxis(ax, [1 nTime]);
    if showColorbar
        cb = colorbar(ax);
        cb.Label.String = 'Time (samples)';
        cb.FontSize = fontSizes.colorbar;
    end
    title(ax, titleText, 'FontSize', fontSizes.title);
    xlabel(ax, 'X (visual degrees)', 'FontSize', fontSizes.label);
    ylabel(ax, 'Y (visual degrees)', 'FontSize', fontSizes.label);
    set(ax, 'FontSize', fontSizes.tick);
end

function plot_distance_axes(ax, paths, titleText, fontSizes)
% plot_distance_axes
%
% Purpose:
%   Plot the time-by-time distance matrix for dot paths into provided axes.
%
% Inputs:
%   ax        : target axes handle
%   paths     : trials × 2 × time array of dot positions
%   titleText : subplot title string
%   fontSizes : struct of typography sizes

    distMean = compute_distance_matrix(paths);
    imagesc(ax, distMean);
    set(ax, 'YDir', 'normal');
    axis(ax, 'image');
    title(ax, titleText, 'FontSize', fontSizes.title);
    xlabel(ax, 'Time (samples)', 'FontSize', fontSizes.label);
    ylabel(ax, 'Time (samples)', 'FontSize', fontSizes.label);
    cb = colorbar(ax);
    cb.FontSize = fontSizes.colorbar;
    set(ax, 'FontSize', fontSizes.tick);
end

function plot_drsa_axes(ax, dRSA, modelIdx, titleText, showYLabel, showXLabel, fontSizes)
% plot_drsa_axes
%
% Purpose:
%   Plot a single dRSA matrix into the provided axes.
%
% Inputs:
%   ax         : target axes handle
%   dRSA       : time × time × model dRSA matrices
%   modelIdx   : model index to plot
%   titleText  : subplot title string
%   showYLabel : true to show y-axis label
%   showXLabel : true to show x-axis label
%   fontSizes  : struct of typography sizes

    nTime = size(dRSA, 1);
    imagesc(ax, dRSA(:, :, modelIdx));
    set(ax, 'YDir', 'normal');
    axis(ax, 'image');
    title(ax, titleText, 'FontSize', fontSizes.title);
    if showXLabel
        xlabel(ax, 'Time in model (samples)', 'FontSize', fontSizes.label);
    end
    if showYLabel
        ylabel(ax, 'Time in data (samples)', 'FontSize', fontSizes.label);
    end
    cb = colorbar(ax, 'eastoutside');
    cb.FontSize = fontSizes.colorbar;
    cb.LineWidth = 0.5;
    set(ax, 'FontSize', fontSizes.tick);
    hold(ax, 'on');
    plot(ax, 1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
    hold(ax, 'off');
end

function plot_lagged_within_combined(ax, results, conditions, useSwapped, fontSizes)
% plot_lagged_within_combined
%
% Purpose:
%   Plot within-condition lagged curves in a single axes with legends.
%
% Inputs:
%   ax         : target axes handle
%   results    : struct with per-condition dRSA diagonals
%   conditions : cell array of condition labels
%   useSwapped : true to use direction-data dRSA diagonals
%   fontSizes  : struct of typography sizes

    hold(ax, 'on');
    labels = {};
    for iCond = 1:numel(conditions)
        conditionLabel = conditions{iCond};
        if useSwapped
            diagVals = results.(conditionLabel).dRSA_diagonal_direction;
            tVec = results.(conditionLabel).diagTimeVec_direction;
            labelSuffix = 'direction data (dot 1)';
        else
            diagVals = results.(conditionLabel).dRSA_diagonal_position;
            tVec = results.(conditionLabel).diagTimeVec_position;
            labelSuffix = 'position data (dot 1)';
        end
        plot(ax, tVec, diagVals(1, :), 'LineWidth', 1.4);
        labels{end + 1} = sprintf('%s %s', format_condition_label(conditionLabel), 'position dot1');
        plot(ax, tVec, diagVals(2, :), 'LineWidth', 1.4);
        labels{end + 1} = sprintf('%s %s', format_condition_label(conditionLabel), 'position dot2');
        plot(ax, tVec, diagVals(3, :), 'LineWidth', 1.4);
        labels{end + 1} = sprintf('%s %s', format_condition_label(conditionLabel), 'direction dot1');
        plot(ax, tVec, diagVals(4, :), 'LineWidth', 1.4);
        labels{end + 1} = sprintf('%s %s', format_condition_label(conditionLabel), 'direction dot2');
    end
    hold(ax, 'off');
    xlabel(ax, 'Lag (samples)', 'FontSize', fontSizes.label);
    ylabel(ax, 'dRSA (average across diagonal)', 'FontSize', fontSizes.label);
    title(ax, sprintf('Lagged dRSA (%s)', labelSuffix), 'FontSize', fontSizes.title);
    legend(ax, labels, 'Location', 'best', 'FontSize', fontSizes.legend);
    set(ax, 'FontSize', fontSizes.tick);
    grid(ax, 'on');
end

function distMean = compute_distance_matrix(paths)
% compute_distance_matrix
%
% Purpose:
%   Compute mean time-by-time distance matrix across trials.
%
% Inputs:
%   paths : trials × 2 × time array of dot positions
%
% Outputs:
%   distMean : time × time matrix of mean Euclidean distances

    nTrials = size(paths, 1);
    nTime = size(paths, 3);
    distSum = zeros(nTime, nTime);
    for iTrial = 1:nTrials
        pos = squeeze(paths(iTrial, :, :))';
        distSum = distSum + pdist2(pos, pos);
    end
    distMean = distSum ./ nTrials;
end

function axisLimits = compute_paths_limits(results, conditions)
% compute_paths_limits
%
% Purpose:
%   Compute a shared axis range for dot-path plots across conditions/dots.
%
% Inputs:
%   results    : struct with per-condition dot paths
%   conditions : cell array of condition labels
%
% Outputs:
%   axisLimits : 1x4 [xmin xmax ymin ymax] with a small padding

    xAll = [];
    yAll = [];
    for iCond = 1:numel(conditions)
        conditionLabel = conditions{iCond};
        dot1 = results.(conditionLabel).dot1Paths;
        dot2 = results.(conditionLabel).dot2Paths;
        xAll = [xAll; dot1(:, 1, :); dot2(:, 1, :)];
        yAll = [yAll; dot1(:, 2, :); dot2(:, 2, :)];
    end
    xAll = xAll(:);
    yAll = yAll(:);
    xMin = min(xAll);
    xMax = max(xAll);
    yMin = min(yAll);
    yMax = max(yAll);
    xPad = 0.05 * (xMax - xMin);
    yPad = 0.05 * (yMax - yMin);
    axisLimits = [xMin - xPad, xMax + xPad, yMin - yPad, yMax + yPad];
end

function label = format_condition_label(conditionLabel)
% format_condition_label
%
% Purpose:
%   Convert condition labels to report-friendly names.
%
% Inputs:
%   conditionLabel : condition string used in filenames
%
% Outputs:
%   label : formatted label for plots

    if strcmpi(conditionLabel, 'nondeviant')
        label = 'non-deviant';
    elseif strcmpi(conditionLabel, 'deviant')
        label = 'deviant';
    else
        label = conditionLabel;
    end
end

function dRSA = select_drsa(conditionStruct, drsaType, useSwapped)
% select_drsa
%
% Purpose:
%   Select the correct dRSA matrix based on dRSA type and data source.
%
% Inputs:
%   conditionStruct : struct containing per-condition dRSA outputs
%   drsaType        : 'corr' or 'pcr'
%   useSwapped      : true to use direction data from dot 1
%
% Outputs:
%   dRSA : time × time × model dRSA matrices

    if strcmpi(drsaType, 'pcr')
        if useSwapped
            dRSA = conditionStruct.dRSA_direction_pcr;
        else
            dRSA = conditionStruct.dRSA_position_pcr;
        end
    else
        if useSwapped
            dRSA = conditionStruct.dRSA_direction;
        else
            dRSA = conditionStruct.dRSA_position;
        end
    end
end

function [dRSA, dRSA_diagonal, diagTimeVec] = run_drsa_direction(dot1Paths, dot2Paths, paramsIn)
% run_drsa_direction
%
% Purpose:
%   Compute dRSA matrices and their diagonals using dot 1 direction as data.
%
% Inputs:
%   dot1Paths : trials × 2 × time array for dot 1
%   dot2Paths : trials × 2 × time array for dot 2
%   paramsIn  : struct of dRSA parameters (model distances, fs, etc.)
%
% Outputs:
%   dRSA          : time × time × model dRSA matrices
%   dRSA_diagonal : model × time diagonal averages
%   diagTimeVec   : lag vector in samples for plotting
%
% Key assumptions:
%   - Direction is derived from frame-to-frame displacement.

    % Build direction data for dot 1 (unit vectors).
    dot1Dx = diff(dot1Paths(:, 1, :), 1, 3);
    dot1Dy = diff(dot1Paths(:, 2, :), 1, 3);
    dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx);
    dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy);
    dot1Angle = atan2(dot1Dy, dot1Dx);
    dot1Direction = cat(2, cos(dot1Angle), sin(dot1Angle));
    dataDirection = dRSA_concatenate(dot1Direction);

    % Recompute position and direction models to align with data.
    dataPosition = dRSA_concatenate(dot1Paths);
    modelPositionDot1 = dataPosition;
    modelPositionDot2 = dRSA_concatenate(dot2Paths);

    dot2Dx = diff(dot2Paths(:, 1, :), 1, 3);
    dot2Dy = diff(dot2Paths(:, 2, :), 1, 3);
    dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx);
    dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy);
    dot2Angle = atan2(dot2Dy, dot2Dx);
    dot2Direction = cat(2, cos(dot2Angle), sin(dot2Angle));
    modelDirectionDot2 = dRSA_concatenate(dot2Direction);

    modelDirectionDot1 = dRSA_concatenate(dot1Direction);
    model = {modelPositionDot1, modelPositionDot2, modelDirectionDot1, modelDirectionDot2};

    % Build trigger-aligned subsamples (trial-locked).
    trialLen = size(dot1Paths, 3);
    totalTime = size(dataDirection, 2);
    if mod(totalTime, trialLen) ~= 0
        error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
    end
    trialStarts = 1:trialLen:totalTime;
    maskSubsampling = true(1, totalTime);
    maskTrigger = false(1, totalTime);
    maskTrigger(trialStarts) = true;

    opt.PreTrigger = 0;
    opt.PostTrigger = trialLen - 1;
    opt.spacing = 0;
    opt.nSubSamples = numel(trialStarts);
    opt.nIter = 1;
    opt.checkRepetition = 0;

    subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);

    % Copy params and set average window based on trial length.
    params = paramsIn;
    avgHalfWindowSamples = floor((trialLen - 1) / 2);
    params.AverageTime = avgHalfWindowSamples / params.fs;

    % Compute autocorrelation borders and dRSA matrices (direction data).
    Autocorrborder = dRSA_border(model, subsamples, params);
    dRSA_Iter = [];
    for iIter = 1:params.nIter
        CurrSubsamples = subsamples(:, :, iIter);
        dRSAma = dRSA_coreFunction(dataDirection, model, params, ...
            'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);
        dRSA_Iter(iIter, :, :, :) = dRSAma;
    end

    dRSA = mean(dRSA_Iter, 1);
    dRSA = reshape(dRSA, size(dRSA, 2), size(dRSA, 3), size(dRSA, 4));

    % Average across the diagonal for lagged plots.
    diagParams = params;
    diagParams.AverageTime = 2;
    diagParams.fs = params.fs;
    dRSA_diagonal = dRSA_average(dRSA, diagParams);
    diagTimeVec = (1:size(dRSA_diagonal, 2)) - ceil(size(dRSA_diagonal, 2) / 2);
end

function meanDist = compute_mean_distance_from_center(paths, rectSize, rectOrigin)
% compute_mean_distance_from_center
%
% Purpose:
%   Compute mean distance from the stimulus center at each timepoint.
%
% Inputs:
%   paths      : trials × 2 × time array of dot positions
%   rectSize   : 1x2 [width height] in visual degrees (optional)
%   rectOrigin : 1x2 [x0 y0] origin of the rect (default: [0 0])
%
% Outputs:
%   meanDist : 1 × time vector of mean distances

    if nargin < 3
        rectOrigin = [0 0];
    end

    if isempty(rectSize)
        xVals = reshape(paths(:, 1, :), [], 1);
        yVals = reshape(paths(:, 2, :), [], 1);
        xVals = xVals(~isnan(xVals));
        yVals = yVals(~isnan(yVals));
        if isempty(xVals) || isempty(yVals)
            error('Cannot infer bounds: all positions are NaN.');
        end
        rectOrigin = [min(xVals) min(yVals)];
        rectSize = [max(xVals) - min(xVals), max(yVals) - min(yVals)];
    end
    center = rectOrigin + rectSize / 2;

    nTime = size(paths, 3);
    meanDist = zeros(1, nTime);
    for iTime = 1:nTime
        positions = squeeze(paths(:, :, iTime)); % trials x 2
        deltas = positions - center;
        distances = sqrt(sum(deltas.^2, 2));
        validMask = ~any(isnan(positions), 2);
        if any(validMask)
            meanDist(iTime) = mean(distances(validMask));
        else
            meanDist(iTime) = NaN;
        end
    end
end

function rectSize = load_rect_size(inputFile)
% load_rect_size
%
% Purpose:
%   Load rectSize from the experiment input file if present.
%
% Inputs:
%   inputFile : path to MovDot_SubXX.mat
%
% Outputs:
%   rectSize : 1x2 [width height] or [] if unavailable

    rectSize = [];
    if exist(inputFile, 'file')
        cfgData = load(inputFile, 'Cfg');
        if isfield(cfgData, 'Cfg') && isfield(cfgData.Cfg, 'rectSize')
            rectSize = cfgData.Cfg.rectSize;
        end
    end
end

function add_shared_colorbar(ax, fontSizes)
% add_shared_colorbar
%
% Purpose:
%   Add a single colorbar to a representative axes for dot-path plots.
%
% Inputs:
%   ax        : axes handle to attach the colorbar
%   fontSizes : struct of typography sizes

    cb = colorbar(ax);
    cb.Label.String = 'Time (samples)';
    cb.FontSize = fontSizes.colorbar;
end

function save_figure(figHandle, outFile)
% save_figure
%
% Purpose:
%   Save a figure to disk and close it to free resources.
%
% Inputs:
%   figHandle : handle to the MATLAB figure
%   outFile   : full path to the output PNG file

    % Use on-screen sizing for consistent export.
    set(figHandle, 'PaperPositionMode', 'auto');
    print(figHandle, outFile, '-dpng', '-r200');
    close(figHandle);
end
