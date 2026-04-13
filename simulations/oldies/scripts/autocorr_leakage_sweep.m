%% Position Autocorrelation Leakage Sweep (pathScale x curvatureScale)
% Script: autocorr_leakage_sweep.m
%
% Purpose:
%   Quantify how position-model autocorrelation leakage changes when path
%   spatial extent (pathScale) and turning-rate dynamics (curvatureScale)
%   are varied in a controlled, analysis-side sweep.
%
% Why this script exists:
%   - dRSA timing precision is limited when position autocorrelation stays
%     high at off-diagonal lags.
%   - Before changing stimulus generation code, this script provides a
%     reproducible parameter map showing which combinations increase or
%     decrease leakage width.
%
% Core data flow:
%   1) Load condition-specific dot paths for one subject.
%   2) Transform each trial with (pathScale, curvatureScale):
%      - pathScale rescales per-frame displacement magnitude,
%      - curvatureScale rescales frame-to-frame heading increments.
%   3) Recompute position autocorrelation via dRSA (position vs itself).
%   4) Extract leakage width (first positive lag where autocorr < threshold).
%   5) Save CSV + Markdown + MAT + heatmap for quick comparison.
%
% Example usage (from repo root in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   run('simulations/scripts/autocorr_leakage_sweep.m');
%
% Example usage (non-interactive terminal run):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch "run('simulations/scripts/autocorr_leakage_sweep.m')"
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat
%   - simulations/functions/* dRSA utilities
%
% Outputs:
%   - simulations/output/subXX/autocorr_sweep/position_autocorr_sweep_*.csv
%   - simulations/output/subXX/autocorr_sweep/position_autocorr_sweep_*.md
%   - simulations/output/subXX/autocorr_sweep/position_autocorr_sweep_*.mat
%   - simulations/output/subXX/autocorr_sweep/position_autocorr_sweep_*.png
%
% Key assumptions:
%   - Input paths are trials x 2 x time in visual degrees.
%   - Sweep transformations are analysis-side diagnostics and do not enforce
%     screen-bound feasibility after scaling.
%   - Autocorrelation is computed with corr-style dRSA (position model vs
%     itself), then summarized with lag-diagonal means.

clear all
close all
clc

%% Resolve paths and add dependencies
% Data flow: script location -> repo paths -> MATLAB path setup.
scriptDir = fileparts(mfilename('fullpath'));
simDir = fileparts(scriptDir);
repoRoot = fileparts(simDir);
addpath(scriptDir);
addpath(simDir);
addpath(fullfile(simDir, 'functions'));
addpath(fullfile(simDir, 'debug'));
addpath(repoRoot);

%% Sweep configuration
% Data flow: user-editable config -> file selection + metric extraction.
participantNumber = 98;
inputCondition = 'nondeviant'; % 'nondeviant' or 'deviant'
useCenterRelative = true; % true matches current simulation scripts.
plotSampleRateHz = 120; % keep aligned with simulation defaults.

% Leakage metric: first positive lag where autocorr drops below threshold.
autocorrDropThreshold = 0.20;

% Parameter grid to test.
pathScales = [0.50, 0.75, 1.00];
curvatureScales = [0.50, 0.75, 1.00, 1.25, 1.50];

% Console verbosity for helper functions.
suppressDispText = 1;

%% Build/load condition-specific simulation input
% Data flow: MovDot_SubXX -> condition file -> dot1 paths tensor.
inputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
if ~exist(inputFile, 'file')
    error('Input stimulus file not found: %s', inputFile);
end

simulationInputDir = fullfile(simDir, 'input');
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
build_movdot_simulation_inputs(inputFile, ...
    'OutputDir', simulationInputDir, ...
    'MoveDotScript', moveDotScript, ...
    'suppressDispText', suppressDispText);

conditionFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_%s.mat', participantNumber, lower(inputCondition)));
if ~exist(conditionFile, 'file')
    error('Condition input file not found: %s', conditionFile);
end
conditionData = load(conditionFile);

if useCenterRelative && isfield(conditionData, 'dot1GreenPathsCenterRelative')
    baseDot1Paths = conditionData.dot1GreenPathsCenterRelative;
    pathFrame = 'center-relative';
elseif isfield(conditionData, 'dot1GreenPaths')
    baseDot1Paths = conditionData.dot1GreenPaths;
    pathFrame = 'absolute-rect';
else
    error('dot1 path fields not found in %s', conditionFile);
end

if ~isnumeric(baseDot1Paths) || ndims(baseDot1Paths) ~= 3 || size(baseDot1Paths, 2) ~= 2
    error('dot1 paths must be trials x 2 x time. Found size: %s', mat2str(size(baseDot1Paths)));
end
baseDot1Paths = double(baseDot1Paths);

%% Prepare output folders and run metadata
% Data flow: participant + condition -> deterministic output paths.
subjectOutputDir = fullfile(simDir, 'output', sprintf('sub%02d', participantNumber));
analysisOutputDir = fullfile(subjectOutputDir, 'autocorr_sweep');
if ~exist(analysisOutputDir, 'dir')
    mkdir(analysisOutputDir);
end

timestampTag = datestr(now, 'yyyymmddTHHMMSS');
runTag = sprintf('sub%02d_%s_%s', participantNumber, lower(inputCondition), timestampTag);

%% Sweep loop over pathScale x curvatureScale
% Data flow: transformed paths -> autocorr matrix -> leakage metrics table.
nPath = numel(pathScales);
nCurv = numel(curvatureScales);
dropLagSecMatrix = nan(nPath, nCurv);
dropLagSampleMatrix = nan(nPath, nCurv);
meanAbsOffDiagMatrix = nan(nPath, nCurv);
lagCurves = cell(nPath, nCurv);
lagSamplesRef = [];

resultRows = repmat(struct( ...
    'participantNumber', participantNumber, ...
    'inputCondition', string(inputCondition), ...
    'pathCoordinateFrame', string(pathFrame), ...
    'pathScale', nan, ...
    'curvatureScale', nan, ...
    'dropThreshold', autocorrDropThreshold, ...
    'dropLagSamples', nan, ...
    'dropLagSeconds', nan, ...
    'meanAbsOffDiagCorr', nan), nPath * nCurv, 1);

rowIdx = 0;
for iPath = 1:nPath
    for iCurv = 1:nCurv
        rowIdx = rowIdx + 1;
        currPathScale = pathScales(iPath);
        currCurvScale = curvatureScales(iCurv);

        transformedPaths = local_transform_paths(baseDot1Paths, currPathScale, currCurvScale);
        [autocorrMat, lagSamples, lagCorr] = local_compute_position_autocorr( ...
            transformedPaths, plotSampleRateHz, suppressDispText);

        [dropLagSamples, dropLagSeconds, meanAbsOffDiagCorr] = ...
            local_extract_leakage_metrics(lagSamples, lagCorr, autocorrDropThreshold, plotSampleRateHz);

        dropLagSecMatrix(iPath, iCurv) = dropLagSeconds;
        dropLagSampleMatrix(iPath, iCurv) = dropLagSamples;
        meanAbsOffDiagMatrix(iPath, iCurv) = meanAbsOffDiagCorr;
        lagCurves{iPath, iCurv} = lagCorr;
        lagSamplesRef = lagSamples;

        resultRows(rowIdx).pathScale = currPathScale;
        resultRows(rowIdx).curvatureScale = currCurvScale;
        resultRows(rowIdx).dropLagSamples = dropLagSamples;
        resultRows(rowIdx).dropLagSeconds = dropLagSeconds;
        resultRows(rowIdx).meanAbsOffDiagCorr = meanAbsOffDiagCorr;

        fprintf('Sweep (%d/%d): pathScale=%.2f, curvatureScale=%.2f -> dropLag=%.4fs\n', ...
            rowIdx, nPath * nCurv, currPathScale, currCurvScale, dropLagSeconds);

        % Keep one representative matrix for sanity checks and report context.
        if iPath == 1 && iCurv == 1
            baselineAutocorrMatrix = autocorrMat; %#ok<NASGU>
        end
    end
end

resultsTable = struct2table(resultRows);
resultsTable = sortrows(resultsTable, {'dropLagSeconds', 'meanAbsOffDiagCorr', 'pathScale', 'curvatureScale'}, ...
    {'ascend', 'ascend', 'ascend', 'ascend'});

%% Write outputs (CSV, Markdown, MAT)
% Data flow: results table + matrices -> persistent artifacts for review.
csvFile = fullfile(analysisOutputDir, sprintf('position_autocorr_sweep_%s.csv', runTag));
writetable(resultsTable, csvFile);

mdFile = fullfile(analysisOutputDir, sprintf('position_autocorr_sweep_%s.md', runTag));
local_write_markdown_report(mdFile, resultsTable, ...
    participantNumber, inputCondition, pathFrame, ...
    pathScales, curvatureScales, autocorrDropThreshold);

matFile = fullfile(analysisOutputDir, sprintf('position_autocorr_sweep_%s.mat', runTag));
save(matFile, 'resultsTable', 'dropLagSecMatrix', 'dropLagSampleMatrix', ...
    'meanAbsOffDiagMatrix', 'pathScales', 'curvatureScales', ...
    'autocorrDropThreshold', 'lagCurves', 'lagSamplesRef', ...
    'participantNumber', 'inputCondition', 'pathFrame', 'plotSampleRateHz');

%% Plot and save leakage heatmap
% Data flow: drop-lag matrix -> heatmap figure for quick parameter scanning.
figHeat = figure('Name', 'Position autocorr leakage sweep', 'NumberTitle', 'off');
imagesc(curvatureScales, pathScales, dropLagSecMatrix);
set(gca, 'YDir', 'normal');
xlabel('Curvature scale');
ylabel('Path scale');
title(sprintf('Lag where autocorr < %.2f (s) | sub%02d %s', ...
    autocorrDropThreshold, participantNumber, lower(inputCondition)));
cb = colorbar;
cb.Label.String = 'Leakage width (s)';

for iPath = 1:nPath
    for iCurv = 1:nCurv
        cellVal = dropLagSecMatrix(iPath, iCurv);
        if isnan(cellVal)
            label = 'NaN';
        else
            label = sprintf('%.3f', cellVal);
        end
        text(curvatureScales(iCurv), pathScales(iPath), label, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9);
    end
end

pngFile = fullfile(analysisOutputDir, sprintf('position_autocorr_sweep_%s.png', runTag));
print(figHeat, pngFile, '-dpng', '-r300');

%% Print concise run summary
% Data flow: sorted table -> top rows for terminal-level sanity checks.
fprintf('\nTop 5 settings with smallest leakage width:\n');
disp(resultsTable(1:min(5, height(resultsTable)), ...
    {'pathScale', 'curvatureScale', 'dropLagSeconds', 'dropLagSamples', 'meanAbsOffDiagCorr'}));

fprintf('Saved:\n- %s\n- %s\n- %s\n- %s\n', csvFile, mdFile, matFile, pngFile);

%% Local helper functions
% Data flow: keep transformation, dRSA, metrics, and reporting logic explicit.

function transformedPaths = local_transform_paths(basePaths, pathScale, curvatureScale)
% local_transform_paths
%
% Purpose:
%   Apply analysis-side kinematic scaling to each trial:
%   - pathScale rescales per-step displacement magnitude,
%   - curvatureScale rescales heading increments (turn rate).
%
% Inputs:
%   basePaths       : trials x 2 x time (x/y positions).
%   pathScale       : scalar multiplier on per-frame displacement magnitude.
%   curvatureScale  : scalar multiplier on frame-to-frame heading deltas.
%
% Output:
%   transformedPaths: trials x 2 x time transformed trajectories.

if ~isscalar(pathScale) || pathScale <= 0 || ~isfinite(pathScale)
    error('pathScale must be a finite scalar > 0.');
end
if ~isscalar(curvatureScale) || curvatureScale <= 0 || ~isfinite(curvatureScale)
    error('curvatureScale must be a finite scalar > 0.');
end

[nTrials, ~, nTime] = size(basePaths);
transformedPaths = zeros(size(basePaths));

for iTrial = 1:nTrials
    trialXY = [squeeze(basePaths(iTrial, 1, :)), squeeze(basePaths(iTrial, 2, :))];
    transformedXY = zeros(nTime, 2);
    transformedXY(1, :) = trialXY(1, :);

    if nTime > 1
        stepVec = diff(trialXY, 1, 1);
        stepLen = sqrt(sum(stepVec .^ 2, 2));

        heading = zeros(nTime - 1, 1);
        for iStep = 1:(nTime - 1)
            if stepLen(iStep) <= eps
                if iStep == 1
                    heading(iStep) = 0;
                else
                    heading(iStep) = heading(iStep - 1);
                end
            else
                heading(iStep) = atan2(stepVec(iStep, 2), stepVec(iStep, 1));
            end
        end

        turnDelta = [0; local_wrap_to_pi(diff(heading))];
        scaledHeading = zeros(nTime - 1, 1);
        scaledHeading(1) = heading(1);
        for iStep = 2:(nTime - 1)
            % Scale turning increments while preserving temporal order.
            scaledHeading(iStep) = scaledHeading(iStep - 1) + curvatureScale * turnDelta(iStep);
        end

        scaledStepLen = pathScale * stepLen;
        for iStep = 1:(nTime - 1)
            if scaledStepLen(iStep) <= eps
                transformedXY(iStep + 1, :) = transformedXY(iStep, :);
            else
                stepDelta = scaledStepLen(iStep) * [cos(scaledHeading(iStep)), sin(scaledHeading(iStep))];
                transformedXY(iStep + 1, :) = transformedXY(iStep, :) + stepDelta;
            end
        end
    end

    transformedPaths(iTrial, 1, :) = transformedXY(:, 1);
    transformedPaths(iTrial, 2, :) = transformedXY(:, 2);
end
end

function wrapped = local_wrap_to_pi(angleValues)
% local_wrap_to_pi
%
% Purpose:
%   Wrap input angles (radians) to the interval [-pi, pi).
%
% Input:
%   angleValues : numeric vector/matrix of angles in radians.
%
% Output:
%   wrapped     : angles wrapped to [-pi, pi).
wrapped = mod(angleValues + pi, 2 * pi) - pi;
end

function [autocorrMat, lagSamples, lagCorr] = local_compute_position_autocorr(dotPaths, fs, suppressDispText)
% local_compute_position_autocorr
%
% Purpose:
%   Compute corr-style position autocorrelation matrix using the same dRSA
%   building blocks used in the simulation pipeline.
%
% Inputs:
%   dotPaths        : trials x 2 x time position tensor.
%   fs              : sampling rate (Hz).
%   suppressDispText: 0/1 passed into helper functions.
%
% Outputs:
%   autocorrMat     : time x time autocorrelation matrix.
%   lagSamples      : lag vector in samples for diagonal extraction.
%   lagCorr         : mean correlation per lag (diagonal average).

dataPosition = dRSA_concatenate(dotPaths, [], 0, 'suppressDispText', suppressDispText);
model = {dataPosition};

trialLen = size(dotPaths, 3);
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
subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt, ...
    'suppressDispText', suppressDispText);

params = struct();
params.nIter = 1;
params.fs = fs;
params.AverageTime = floor((trialLen - 1) / 2) / fs;
params.modelToTest = 1;
params.Var = 0.1;
% In dRSA autocorr mode (model vs itself), core expects a single distance
% label string, not a 1x1 cell array.
params.modelDistMeasure = 'euclidean';
params.neuralDistMeasure = 'euclidean';
params.dRSAtype = 'corr';

autocorrMat = dRSA_coreFunction(dataPosition, model, params, ...
    'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', [], ...
    'suppressDispText', suppressDispText);
autocorrMat = reshape(autocorrMat, size(autocorrMat, 1), size(autocorrMat, 2));

nTime = size(autocorrMat, 1);
lagSamples = -(nTime - 1):(nTime - 1);
lagCorr = nan(size(lagSamples));
for iLag = 1:numel(lagSamples)
    lagCorr(iLag) = mean(diag(autocorrMat, lagSamples(iLag)), 'omitnan');
end
end

function [dropLagSamples, dropLagSeconds, meanAbsOffDiagCorr] = local_extract_leakage_metrics( ...
        lagSamples, lagCorr, dropThreshold, fs)
% local_extract_leakage_metrics
%
% Purpose:
%   Convert a lag-correlation profile into compact leakage metrics.
%
% Inputs:
%   lagSamples    : lag vector in samples.
%   lagCorr       : mean correlation per lag.
%   dropThreshold : threshold used to define leakage width.
%   fs            : sampling rate (Hz) for seconds conversion.
%
% Outputs:
%   dropLagSamples    : first non-negative lag where corr < threshold.
%   dropLagSeconds    : same value in seconds.
%   meanAbsOffDiagCorr: mean absolute correlation for |lag|>0.

if ~isscalar(dropThreshold) || ~isfinite(dropThreshold)
    error('dropThreshold must be a finite scalar.');
end

offDiagMask = lagSamples ~= 0;
meanAbsOffDiagCorr = mean(abs(lagCorr(offDiagMask)), 'omitnan');

positiveMask = lagSamples >= 0;
positiveLags = lagSamples(positiveMask);
positiveCorr = lagCorr(positiveMask);
dropIdx = find(positiveCorr < dropThreshold, 1, 'first');

if isempty(dropIdx)
    dropLagSamples = NaN;
    dropLagSeconds = NaN;
else
    dropLagSamples = positiveLags(dropIdx);
    dropLagSeconds = dropLagSamples / fs;
end
end

function local_write_markdown_report(reportPath, resultsTable, ...
        participantNumber, inputCondition, pathFrame, ...
        pathScales, curvatureScales, dropThreshold)
% local_write_markdown_report
%
% Purpose:
%   Write a lightweight Markdown summary to accompany the CSV output.
%
% Inputs:
%   reportPath         : markdown output path.
%   resultsTable       : sorted sweep results table.
%   participantNumber  : subject identifier.
%   inputCondition     : condition label.
%   pathFrame          : coordinate frame label.
%   pathScales         : tested path-scale vector.
%   curvatureScales    : tested curvature-scale vector.
%   dropThreshold      : leakage threshold used for width extraction.

fid = fopen(reportPath, 'w');
if fid == -1
    error('Could not open markdown report for writing: %s', reportPath);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '# Position Autocorrelation Leakage Sweep\n\n');
fprintf(fid, '- Subject: `sub%02d`\n', participantNumber);
fprintf(fid, '- Condition: `%s`\n', lower(inputCondition));
fprintf(fid, '- Path frame: `%s`\n', pathFrame);
fprintf(fid, '- Drop threshold: `%.3f`\n', dropThreshold);
fprintf(fid, '- pathScales: `%s`\n', mat2str(pathScales));
fprintf(fid, '- curvatureScales: `%s`\n\n', mat2str(curvatureScales));

fprintf(fid, '## Best settings (lowest leakage width)\n\n');
nRows = min(10, height(resultsTable));
fprintf(fid, '| Rank | pathScale | curvatureScale | dropLagSeconds | dropLagSamples | meanAbsOffDiagCorr |\n');
fprintf(fid, '| --- | --- | --- | --- | --- | --- |\n');
for iRow = 1:nRows
    fprintf(fid, '| %d | %.2f | %.2f | %.6f | %.0f | %.6f |\n', ...
        iRow, ...
        resultsTable.pathScale(iRow), ...
        resultsTable.curvatureScale(iRow), ...
        resultsTable.dropLagSeconds(iRow), ...
        resultsTable.dropLagSamples(iRow), ...
        resultsTable.meanAbsOffDiagCorr(iRow));
end
fprintf(fid, '\n');
end
