% pipeline_recursive.m
%
% Purpose:
%   Run dRSA repeatedly on random trial subsets (size = subsetTrialCount)
%   to inspect variability across subsets. Each run uses params.nIter = 1
%   and re-computes inputs, models, triggers, and dRSA matrices for that
%   subset only. No outputs are saved to disk.
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
%   - nRecursiveRuns controls how many dRSA runs to execute (N times).
%   - subsetTrialCount sets the number of trials per run (unique per subset).
%   - shuffleDot2Trials decouples dot2 trial order from dot1 (diagnostic).
%   - rngSeed controls reproducibility of subset sampling.
%
% Inputs:
%   - simulations/input/MovDot_SubXX_{nondeviant|deviant}.mat with:
%       dot1GreenPaths (trials × 2 × time) and dot2YellowPaths (same shape)
%
% Outputs:
%   - None saved. Figures are displayed only for inspection.
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
participantNumber = 93;
inputCondition = 'nondeviant'; % 'nondeviant' or 'deviant' to match output files.
inputCondition = lower(inputCondition);
shuffleDot2Trials = false; % true to shuffle dot2 trials relative to dot1 before dRSA

% Recursive dRSA configuration (N runs with unique trials per subset).
nRecursiveRuns = 4; % N times to run dRSA on different trial subsets.
subsetTrialCount = 20; % number of trials per run (unique within each subset).

% Random seed for subset sampling (set numeric for reproducibility).
% Use 'shuffle' to re-seed from the current time.
rngSeed = 'shuffle';

validConditions = {'nondeviant', 'deviant'};
if ~ismember(inputCondition, validConditions)
    error('inputCondition must be one of: %s', strjoin(validConditions, ', '));
end
simulationInputDir = fullfile(simDir, 'input');

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

%% Preallocate a results struct for inspection in the workspace
% Data flow: per-run dRSA outputs -> results struct array (no file saving).
results = struct( ...
    'subsetIdx', [], ...
    'dRSA_position', [], ...
    'dRSA_direction', [], ...
    'modelNames', [], ...
    'paramsCore', [], ...
    'trialLen', []);
results = repmat(results, nRecursiveRuns, 1);

%% Run dRSA recursively across random trial subsets
% Data flow: subset indices -> per-subset data/models -> dRSA outputs + plots.
for runIdx = 1:nRecursiveRuns
    % Select a random subset of trials (unique within the subset).
    subsetIdx = randperm(nTrialsTotal, subsetTrialCount);
    dot1Trials = dot1TrialsAll(subsetIdx, :, :);
    dot2Trials = dot2TrialsAll(subsetIdx, :, :);

    % Plot subset paths for quick sanity checks (no file output).
    figPaths = figure('Name', sprintf('Run %d: dot paths', runIdx), 'NumberTitle', 'off');
    tiledlayout(figPaths, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    plot_paths_no_save(dot1Trials, 'Title', sprintf('Run %d: Dot 1 paths', runIdx));
    nexttile;
    plot_paths_no_save(dot2Trials, 'Title', sprintf('Run %d: Dot 2 paths', runIdx));

    %% Prepare data, models, and trial-locked subsamples
    % Data flow: trial arrays -> per-dot position + direction models -> concatenated time series.
    % Position models use raw x/y paths; direction models use per-timepoint angles.
    dataPosition = dRSA_concatenate(dot1Trials); % position data from dot 1

    % Position models (per dot).
    modelPositionDot1 = dataPosition;
    modelPositionDot2 = dRSA_concatenate(dot2Trials);

    % Direction models (per dot).
    % Data flow: dot positions -> frame-to-frame displacement -> angle -> unit vectors.
    dot1Dx = diff(dot1Trials(:, 1, :), 1, 3); % frame-to-frame x displacement (trials x 1 x time-1)
    dot1Dy = diff(dot1Trials(:, 2, :), 1, 3); % frame-to-frame y displacement (trials x 1 x time-1)
    dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx); % pad so length matches original time
    dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy); % pad so length matches original time
    dot1Angle = atan2(dot1Dy, dot1Dx); % per-timepoint direction angle in radians
    dot1Direction = cat(2, cos(dot1Angle), sin(dot1Angle)); % unit vectors [cos; sin]
    modelDirectionDot1 = dRSA_concatenate(dot1Direction);

    dot2Dx = diff(dot2Trials(:, 1, :), 1, 3); % frame-to-frame x displacement for dot 2
    dot2Dy = diff(dot2Trials(:, 2, :), 1, 3); % frame-to-frame y displacement for dot 2
    dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx); % pad so length matches original time
    dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy); % pad so length matches original time
    dot2Angle = atan2(dot2Dy, dot2Dx); % per-timepoint direction angle in radians
    dot2Direction = cat(2, cos(dot2Angle), sin(dot2Angle)); % unit vectors [cos; sin]
    modelDirectionDot2 = dRSA_concatenate(dot2Direction);

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
    subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);

    %% Configure models and core dRSA parameters
    % Data flow: models + subsamples -> dRSA params -> autocorrelation border.
    model = {modelPositionDot1, modelPositionDot2, modelDirectionDot1, modelDirectionDot2};
    modelNames = {'position dot1', 'position dot2', 'direction dot1', 'direction dot2'};

    params.nIter = 1; % force single-iteration dRSA per subset
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
        Autocorrborder = dRSA_border(model, subsamples, params);
    end


    %% Run dRSA for position data (single iteration)
    % Data flow: position data + models -> dRSA matrices (time x time x model).
    % Match neural distance to the position model (model 1) when data are positions.
    params.neuralDistMeasure = params.modelDistMeasure{1};
    dRSA_position = dRSA_coreFunction(dataPosition, model, params, ...
        'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', Autocorrborder);

    %% Run dRSA for direction data (single iteration)
    % Data flow: direction data + models -> dRSA matrices (time x time x model).
    % Match neural distance to the direction model (model 3) when data are directions.
    params.neuralDistMeasure = params.modelDistMeasure{3};
    dRSA_direction = dRSA_coreFunction(dataDirection, model, params, ...
        'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', Autocorrborder);

    %% Plot dRSA matrices for this subset
    % Data flow: dRSA per data type -> grid of heatmaps for inspection.
    nModels = size(dRSA_position, 3);
    nTime = size(dRSA_position, 1);
    dRSAAll = {dRSA_position, dRSA_direction};
    rowNames = {'Data: position dot1', 'Data: direction dot1'};

    figMatrices = figure('Name', sprintf('Run %d: dRSA matrices', runIdx), 'NumberTitle', 'off');
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
            xlabel(sprintf('Time in %s (samples)', modelNames{iModel}));
            ylabel(sprintf('Time in %s (samples)', rowNames{iRow}));
            hold on;
            plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
            hold off;
            if iModel == 1
                rowAxes(iRow) = gca;
            end
        end
    end
    for iRow = 1:numel(rowAxes)
        if isgraphics(rowAxes(iRow))
            text(rowAxes(iRow), -0.35, 0.5, rowNames{iRow}, ...
                'Units', 'normalized', 'Rotation', 90, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end

    %% Store run outputs in the workspace
    % Data flow: per-run outputs -> results struct (no file saving).
    results(runIdx).subsetIdx = subsetIdx;
    results(runIdx).dRSA_position = dRSA_position;
    results(runIdx).dRSA_direction = dRSA_direction;
    results(runIdx).modelNames = modelNames;
    results(runIdx).paramsCore = params;
    results(runIdx).trialLen = trialLen;

    % Ensure figures are rendered before moving to the next run.
    drawnow;
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
    parser.parse(paths, varargin{:});

    opts = parser.Results;

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
