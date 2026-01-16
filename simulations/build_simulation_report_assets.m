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
useSavedResults = true;
saveResults = true;
participantNumber = 98;
conditions = {'nondeviant', 'deviant'};
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
params.modelToTest = [1 2];
params.Var = 0.1;
params.modelDistMeasure = {'euclidean', 'euclidean'};
params.neuralDistMeasure = 'euclidean';
params.dRSAtype = 'corr';
params.modeltoRegressout = {2 1};

paramsPCR = params;
paramsPCR.dRSAtype = 'PCR';

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
    loaded = load(assetsFile, 'results', 'participantNumber', 'conditions');
    results = loaded.results;
    % Enforce exact match when using cached results.
    if loaded.participantNumber ~= participantNumber || ~isequal(loaded.conditions, conditions)
        error(['Cached results were computed for different settings. ', ...
            'Set useSavedResults = false to recompute with current settings.']);
    end
else
    results = compute_results(simulationInputDir, participantNumber, conditions, params, paramsPCR);
end

if saveResults
    save(assetsFile, 'results', 'participantNumber', 'conditions');
end

%% Build grouped figures and save to disk
% Data flow: results -> grouped figures -> PNGs in output directory.
fig = plot_grouped_paths(results, conditions, fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('dot_paths_grouped_%s.png', subjectLabel)));

fig = plot_grouped_distances(results, conditions, fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('distance_matrices_grouped_%s.png', subjectLabel)));

fig = plot_grouped_drsa_combined(results, conditions, 'corr', fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('drsa_matrices_corr_%s.png', subjectLabel)));

fig = plot_grouped_drsa_combined(results, conditions, 'pcr', fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('drsa_matrices_pcr_%s.png', subjectLabel)));

fig = plot_all_lagged(results, conditions, fontSizes, subjectLabel);
save_figure(fig, fullfile(figureOutputDir, sprintf('lagged_drsa_%s.png', subjectLabel)));

%% Save report assets metadata
% Data flow: parameters + paths -> MAT file for report generator.
save(assetsFile, 'results', 'participantNumber', 'conditions', 'params', ...
    'paramsPCR', 'fontSizes', 'figureOutputDir');

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

function [dRSA, dRSA_diagonal, diagTimeVec] = run_drsa(dot1Paths, dot2Paths, paramsIn)
% run_drsa
%
% Purpose:
%   Compute dRSA matrices and their averaged diagonals for dot 1 and dot 2.
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
%   - dot paths share the same number of trials and time samples.
%   - Trial length divides the concatenated time series evenly.

    % Prepare concatenated data and models.
    data = dRSA_concatenate(dot1Paths);
    model1 = data;
    model2 = dRSA_concatenate(dot2Paths);
    model = {model1, model2};
    Y = data;

    % Build trigger-aligned subsamples (trial-locked).
    trialLen = size(dot1Paths, 3);
    totalTime = size(data, 2);
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

function results = compute_results(simulationInputDir, participantNumber, conditions, params, paramsPCR)
% compute_results
%
% Purpose:
%   Compute per-condition dRSA outputs for corr and PCR, and cache results.
%
% Inputs:
%   simulationInputDir : directory with condition-specific dot inputs
%   participantNumber : numeric participant ID
%   conditions        : cell array of condition labels
%   params            : dRSA parameters for corr
%   paramsPCR         : dRSA parameters for PCR
%
% Outputs:
%   results : struct containing dot paths and dRSA outputs per condition

    results = struct();
    for iCond = 1:numel(conditions)
        conditionLabel = conditions{iCond};
        [dot1Paths, dot2Paths] = load_condition_paths(simulationInputDir, ...
            participantNumber, conditionLabel);

        % Compute dRSA matrices (corr; dot 1 as data).
        [dRSA, dRSA_diagonal, diagTimeVec] = run_drsa(dot1Paths, dot2Paths, params);

        % Compute dRSA matrices (corr; dot 2 as data).
        [dRSA_swapped, dRSA_diagonal_swapped, diagTimeVec_swapped] = ...
            run_drsa(dot2Paths, dot1Paths, params);

        % Compute dRSA matrices (PCR; dot 1 as data).
        [dRSA_pcr, dRSA_diagonal_pcr, diagTimeVec_pcr] = ...
            run_drsa(dot1Paths, dot2Paths, paramsPCR);

        % Compute dRSA matrices (PCR; dot 2 as data).
        [dRSA_swapped_pcr, dRSA_diagonal_swapped_pcr, diagTimeVec_swapped_pcr] = ...
            run_drsa(dot2Paths, dot1Paths, paramsPCR);

        results.(conditionLabel).dot1Paths = dot1Paths;
        results.(conditionLabel).dot2Paths = dot2Paths;
        results.(conditionLabel).dRSA = dRSA;
        results.(conditionLabel).dRSA_swapped = dRSA_swapped;
        results.(conditionLabel).dRSA_diagonal = dRSA_diagonal;
        results.(conditionLabel).diagTimeVec = diagTimeVec;
        results.(conditionLabel).dRSA_diagonal_swapped = dRSA_diagonal_swapped;
        results.(conditionLabel).diagTimeVec_swapped = diagTimeVec_swapped;
        results.(conditionLabel).dRSA_pcr = dRSA_pcr;
        results.(conditionLabel).dRSA_swapped_pcr = dRSA_swapped_pcr;
        results.(conditionLabel).dRSA_diagonal_pcr = dRSA_diagonal_pcr;
        results.(conditionLabel).diagTimeVec_pcr = diagTimeVec_pcr;
        results.(conditionLabel).dRSA_diagonal_swapped_pcr = dRSA_diagonal_swapped_pcr;
        results.(conditionLabel).diagTimeVec_swapped_pcr = diagTimeVec_swapped_pcr;
    end
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

function figHandle = plot_grouped_drsa_combined(results, conditions, drsaType, fontSizes, subjectLabel)
% plot_grouped_drsa_combined
%
% Purpose:
%   Plot dRSA matrices in a single figure with two rows:
%   top row = dot 1 as data, bottom row = dot 2 as data.
%
% Inputs:
%   results    : struct with per-condition dRSA outputs
%   conditions : cell array of condition labels
%   drsaType   : 'corr' or 'pcr' to select the dRSA variant
%   fontSizes  : struct of typography sizes
%
% Outputs:
%   figHandle : figure handle for the combined dRSA plot

    nCond = numel(conditions);
    nModels = 2;
    nCols = nCond * nModels;
    figHandle = figure('Name', 'dRSA matrices combined', ...
        'NumberTitle', 'off', 'Visible', 'off');
    figHandle.Position = [100 100 2400 900];
    sgtitle(sprintf('dRSA matrices (%s, %s)', upper(drsaType), subjectLabel), ...
        'FontSize', fontSizes.title + 2);

    for iRow = 1:2
        for iCond = 1:nCond
            conditionLabel = conditions{iCond};
            if iRow == 1
                dRSA = select_drsa(results.(conditionLabel), drsaType, false);
                rowLabel = 'dot 1 as data';
            else
                dRSA = select_drsa(results.(conditionLabel), drsaType, true);
                rowLabel = 'dot 2 as data';
            end
            for iModel = 1:nModels
                colIndex = (iCond - 1) * nModels + iModel;
                subplot(2, nCols, (iRow - 1) * nCols + colIndex);
                ax = gca;
                plot_drsa_axes(ax, dRSA, iModel, ...
                    sprintf('Model %d (%s)', iModel, conditionLabel), ...
                    iCond == 1, iRow == 2, fontSizes);
                if iCond == 1 && iModel == 1
                    ylabel(ax, sprintf('%s\nTime (samples)', rowLabel));
                end
            end
        end
    end
end

function figHandle = plot_all_lagged(results, conditions, fontSizes, subjectLabel)
% plot_all_lagged
%
% Purpose:
%   Plot all lagged dRSA curves in one figure with two panels:
%   - dot 1 as data
%   - dot 2 as data
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

    % Panel 1: dot 1 as data.
    subplot(1, 2, 1);
    plot_lagged_within_combined(gca, results, conditions, false, fontSizes);

    % Panel 2: dot 2 as data.
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
        xlabel(ax, 'Time (samples)', 'FontSize', fontSizes.label);
    end
    if showYLabel
        ylabel(ax, 'Time (samples)', 'FontSize', fontSizes.label);
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
%   useSwapped : true to use swapped dRSA diagonals
%   fontSizes  : struct of typography sizes

    hold(ax, 'on');
    labels = {};
    for iCond = 1:numel(conditions)
        conditionLabel = conditions{iCond};
        if useSwapped
            diagVals = results.(conditionLabel).dRSA_diagonal_swapped;
            tVec = results.(conditionLabel).diagTimeVec_swapped;
            labelSuffix = 'dot 2 as data';
        else
            diagVals = results.(conditionLabel).dRSA_diagonal;
            tVec = results.(conditionLabel).diagTimeVec;
            labelSuffix = 'dot 1 as data';
        end
        plot(ax, tVec, diagVals(1, :), 'LineWidth', 1.4);
        labels{end + 1} = sprintf('%s model 1', format_condition_label(conditionLabel));
        plot(ax, tVec, diagVals(2, :), 'LineWidth', 1.4);
        labels{end + 1} = sprintf('%s model 2', format_condition_label(conditionLabel));
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
%   useSwapped      : true to use dot 2 as data
%
% Outputs:
%   dRSA : time × time × model dRSA matrices

    if strcmpi(drsaType, 'pcr')
        if useSwapped
            dRSA = conditionStruct.dRSA_swapped_pcr;
        else
            dRSA = conditionStruct.dRSA_pcr;
        end
    else
        if useSwapped
            dRSA = conditionStruct.dRSA_swapped;
        else
            dRSA = conditionStruct.dRSA;
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
