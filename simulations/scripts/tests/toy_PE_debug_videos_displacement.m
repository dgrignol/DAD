% toy_PE_debug_videos_displacement.m
%
% Purpose:
%   Generate synchronized debug videos for post-deviance PE inspection in
%   displacement-aware stimuli (for example v16), where path visualization
%   recentering should preserve the observed jump at deviant onset.
%
% What this script produces:
%   - 3 path videos (observed, predicted, PE) with one fixed color per trial.
%   - 3 RDM videos (observed, predicted, PE) with dissimilarity colorbars.
%   - 1 combined panel video (2x3):
%       row 1 = paths [observed, predicted, PE]
%       row 2 = RDMs  [observed, predicted, PE]
%   - 1 combined panel video (3x3, rate view):
%       row 1 = paths [observed, predicted, PE]
%       row 2 = RDMs  [observed, predicted, PE]
%       row 3 = RDM change-rate curves [observed, predicted, PE]
%
% Why this script exists:
%   The static plots in toy_PE_debug.m are useful snapshots, but they do not
%   show how trial geometry and RDM structure evolve frame-by-frame after
%   deviance. This script makes that temporal evolution explicit.
%
% Data flow summary:
%   1) Load (or build) observed/predicted deviant center-relative paths.
%   2) Derive deviant onset frame and define a pre-deviant anchor frame.
%   3) Build full-trial PE paths as observed - predicted.
%   4) Optionally recenter each full trial at the anchor frame (default:
%      end of first half, i.e., cutFrame-1), then cut to post-deviance for
%      path visualization so displacement remains visible at sample 1.
%   5) Compute frame-wise trial RDMs (default: raw post-deviance paths,
%      independent of recentering used for path visualization).
%   6) Compute frame-wise RDM change rate (mean absolute off-diagonal delta).
%   7) Render and save 6 individual videos + 2 combined panel videos.
%
% Usage example (from this folder):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/tests');
%   toy_PE_debug_videos_displacement;
%
% Usage example (from anywhere):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/tests');
%   toy_PE_debug_videos_displacement;
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat
%   - experiment/input_files/MovDot_SubXX_predicted.mat
%   - simulations/input/MovDot_SubXX_deviant.mat
%   - simulations/input/MovDot_SubXX_predicted_deviant.mat
%
% Outputs:
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_paths_observed_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_paths_predicted_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_paths_pe_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_rdm_observed_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_rdm_predicted_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_rdm_pe_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_panel_*.mp4
%   - simulations/scripts/tests/derivatives/subXX/videos/toy_PE_debug_panel3x3_rate_*.mp4
%
% Assumptions:
%   - Center-relative path field dot1GreenPathsCenterRelative exists.
%   - Observed and predicted deviant paths are shape-aligned.
%   - Trial count is >= 2 for meaningful RDM construction with pdist.
%   - Default anchor uses a global deviant frame from deviantOnset; this is
%     exact when onset is fixed and approximate when onset variance exists.

clear all
close all

%% Resolve repository paths and add dependencies
% Data flow: script location -> repo root -> helper paths for simulation IO.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));
moveDotScript = fullfile(repoRoot, 'experiment', 'MoveDot1_experiment_vX.m');
addpath(fullfile(repoRoot, 'simulations'));
addpath(fullfile(repoRoot, 'simulations', 'debug'));
addpath(fullfile(repoRoot, 'simulations', 'functions'));

%% User configuration
% Data flow: config -> loading/cropping/recentering -> frame and video setup.
participantNumber = 75;
deviantOnset = 0.5;
enforcePostDevCut = true;
recenterPathsForDisplay = true; % true recenters path videos only.
recenterAnchorMode = 'preDeviantEnd'; % default keeps deviant displacement visible after cut.
suppressDispText = 1;

frameStride = 1; % use >1 to reduce output duration and render time.
plotSampleRateHz = 120;
videoFrameRate = 20;
videoQuality = 95;
figureVisibility = 'off'; % set to 'on' to preview while rendering.
showProgressPrint = true; % true prints stage/frame progress to console.
progressPrintEveryFrames = 20; % frame print cadence when showProgressPrint=true.

pathDistMeasure = 'euclidean';
rdmColormapName = 'parula';
useRecenteredPathsForRDM = false; % default false keeps RDM/rate on raw post-deviance paths.
rdmChangeMetric = 'meanAbsOffDiag'; % currently supported: meanAbsOffDiag

saveIndividualVideos = true;
savePanelVideo2x3 = true;
savePanelVideo3x3Rate = true;

testsDerivativesDir = fullfile(scriptDir, 'derivatives');
local_ensure_dir(testsDerivativesDir);
if isstring(rdmChangeMetric) && isscalar(rdmChangeMetric)
    rdmChangeMetric = char(rdmChangeMetric);
end
rdmChangeMetric = lower(strtrim(rdmChangeMetric));
if ~(isscalar(useRecenteredPathsForRDM) && ...
        (islogical(useRecenteredPathsForRDM) || isnumeric(useRecenteredPathsForRDM)))
    error('useRecenteredPathsForRDM must be a scalar logical/numeric value.');
end
if logical(useRecenteredPathsForRDM)
    rdmSourceMode = 'sameaspaths';
else
    rdmSourceMode = 'rawpostdev';
end
if ~strcmp(rdmChangeMetric, 'meanabsoffdiag')
    error('rdmChangeMetric currently supports only ''meanAbsOffDiag''.');
end
if ~(isscalar(showProgressPrint) && ...
        (islogical(showProgressPrint) || isnumeric(showProgressPrint)))
    error('showProgressPrint must be a scalar logical/numeric value.');
end
if ~isscalar(progressPrintEveryFrames) || progressPrintEveryFrames < 1 || ...
        floor(progressPrintEveryFrames) ~= progressPrintEveryFrames
    error('progressPrintEveryFrames must be a positive integer.');
end
if isstring(recenterAnchorMode) && isscalar(recenterAnchorMode)
    recenterAnchorMode = char(recenterAnchorMode);
end
recenterAnchorMode = lower(strtrim(recenterAnchorMode));
if ~any(strcmp(recenterAnchorMode, {'predeviantend', 'deviantsample'}))
    error('recenterAnchorMode must be ''preDeviantEnd'' or ''deviantSample''.');
end

if showProgressPrint
    fprintf('[toy_PE_debug_videos_displacement] Starting run for sub%02d\n', participantNumber);
end

%% Resolve input files and build missing simulation inputs
% Data flow: subject id -> expected input files -> optional simulation-input build.
simulationInputDir = fullfile(repoRoot, 'simulations', 'input');
observedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d.mat', participantNumber));
predictedInputFile = fullfile(repoRoot, 'experiment', 'input_files', ...
    sprintf('MovDot_Sub%02d_predicted.mat', participantNumber));
deviantSimulationFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_deviant.mat', participantNumber));
predictedSimulationFile = fullfile(simulationInputDir, ...
    sprintf('MovDot_Sub%02d_predicted_deviant.mat', participantNumber));

if ~isfile(observedInputFile)
    error('Observed input file not found: %s', observedInputFile);
end
if ~isfile(predictedInputFile)
    error('Predicted input file not found: %s', predictedInputFile);
end

if ~isfile(deviantSimulationFile)
    if showProgressPrint
        fprintf('[toy_PE_debug_videos_displacement] Building missing observed deviant simulation input...\n');
    end
    build_movdot_simulation_inputs(observedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'suppressDispText', suppressDispText);
end
if ~isfile(predictedSimulationFile)
    if showProgressPrint
        fprintf('[toy_PE_debug_videos_displacement] Building missing predicted deviant simulation input...\n');
    end
    build_movdot_simulation_inputs(predictedInputFile, ...
        'OutputDir', simulationInputDir, ...
        'MoveDotScript', moveDotScript, ...
        'AllowMissingConditions', true, ...
        'suppressDispText', suppressDispText);
end

%% Load observed/predicted deviant dot1 paths and validate alignment
% Data flow: .mat files -> center-relative dot1 tensors -> shape validation.
observedData = load(deviantSimulationFile);
predictedData = load(predictedSimulationFile);
requiredField = 'dot1GreenPathsCenterRelative';
if ~isfield(observedData, requiredField)
    error('Observed deviant file lacks %s: %s', requiredField, deviantSimulationFile);
end
if ~isfield(predictedData, requiredField)
    error('Predicted deviant file lacks %s: %s', requiredField, predictedSimulationFile);
end

observedDot1Full = observedData.dot1GreenPathsCenterRelative;
predictedDot1Full = predictedData.dot1GreenPathsCenterRelative;
if ~isequal(size(observedDot1Full), size(predictedDot1Full))
    error(['Observed/predicted dot1 shapes must match. ', ...
        'Observed=%s Predicted=%s'], ...
        mat2str(size(observedDot1Full)), mat2str(size(predictedDot1Full)));
end
if size(observedDot1Full, 1) < 2
    error('Need at least 2 trials for frame-wise RDM videos.');
end

%% Resolve deviant frame and build full-trial PE paths
% Data flow: full trial lengths + onset fraction -> cut frame -> full PE tensor.
fullTrialFrames = size(observedDot1Full, 3);
cutFrame = local_compute_cut_frame(fullTrialFrames, deviantOnset);
peDot1Full = observedDot1Full - predictedDot1Full;

%% Build path-visualization tensors with displacement-aware recentering
% Data flow: full paths -> optional anchor-frame recenter -> optional post-dev cut.
observedDot1PathFull = observedDot1Full;
predictedDot1PathFull = predictedDot1Full;
peDot1PathFull = peDot1Full;

anchorFrame = local_resolve_anchor_frame(fullTrialFrames, cutFrame, recenterAnchorMode);
if recenterPathsForDisplay
    if showProgressPrint
        fprintf(['[toy_PE_debug_videos_displacement] Recentering path-visualization tensors ', ...
            'at frame %d (%s).\n'], anchorFrame, recenterAnchorMode);
    end
    observedDot1PathFull = local_recenter_paths_at_frame(observedDot1PathFull, anchorFrame);
    predictedDot1PathFull = local_recenter_paths_at_frame(predictedDot1PathFull, anchorFrame);
    peDot1PathFull = local_recenter_paths_at_frame(peDot1PathFull, anchorFrame);
    local_assert_anchor_at_origin( ...
        {observedDot1PathFull, predictedDot1PathFull, peDot1PathFull}, ...
        anchorFrame, 1e-10);
end

if enforcePostDevCut
    observedDot1Path = local_cut_from_frame(observedDot1PathFull, cutFrame);
    predictedDot1Path = local_cut_from_frame(predictedDot1PathFull, cutFrame);
    peDot1Path = local_cut_from_frame(peDot1PathFull, cutFrame);
    observedDot1RawWindow = local_cut_from_frame(observedDot1Full, cutFrame);
    predictedDot1RawWindow = local_cut_from_frame(predictedDot1Full, cutFrame);
    peDot1RawWindow = local_cut_from_frame(peDot1Full, cutFrame);
else
    observedDot1Path = observedDot1PathFull;
    predictedDot1Path = predictedDot1PathFull;
    peDot1Path = peDot1PathFull;
    observedDot1RawWindow = observedDot1Full;
    predictedDot1RawWindow = predictedDot1Full;
    peDot1RawWindow = peDot1Full;
end

%% Build frame schedule and choose path source for RDM computation
% Data flow: post-deviance tensors + RDM source option -> path and RDM datasets.
nTime = size(observedDot1Path, 3);
frameIndices = 1:frameStride:nTime;
if frameIndices(end) ~= nTime
    frameIndices = [frameIndices, nTime];
end
timeSeconds = (frameIndices - 1) ./ plotSampleRateHz;

datasetsPaths = {observedDot1Path, predictedDot1Path, peDot1Path};
datasetLabels = {'Observed', 'Predicted', 'PE'};
pathAxisLimits = local_compute_global_axis_limits(datasetsPaths);
trialColors = local_trial_colors(size(observedDot1Path, 1));

if strcmp(rdmSourceMode, 'sameaspaths')
    % Optional branch: tie RDM/rate to the recentered path visualization tensors.
    datasetsRdmSource = datasetsPaths;
else
    % Default branch: keep RDM/rate based on original post-deviance tensors.
    datasetsRdmSource = {observedDot1RawWindow, predictedDot1RawWindow, peDot1RawWindow};
end

rdmObserved = local_compute_frame_rdms(datasetsRdmSource{1}, frameIndices, pathDistMeasure);
rdmPredicted = local_compute_frame_rdms(datasetsRdmSource{2}, frameIndices, pathDistMeasure);
rdmPE = local_compute_frame_rdms(datasetsRdmSource{3}, frameIndices, pathDistMeasure);
datasetsRdm = {rdmObserved, rdmPredicted, rdmPE};
rdmColorLimits = local_compute_rdm_limits(datasetsRdm);
rdmChangeRates = local_compute_rdm_change_rates(datasetsRdm, rdmChangeMetric);
rdmRateYLimits = local_compute_rate_axis_limits(rdmChangeRates);

%% Prepare output folder and filename tags
% Data flow: subject/crop config -> deterministic output filenames.
subjectLabel = sprintf('sub%02d', participantNumber);
videoOutputDir = fullfile(testsDerivativesDir, subjectLabel, 'videos');
local_ensure_dir(videoOutputDir);
if enforcePostDevCut
    runTag = sprintf('postDevFrom_%d', cutFrame);
else
    runTag = 'fullTrial';
end
if recenterPathsForDisplay
    centerTag = local_centering_tag_for_filename(recenterAnchorMode);
    centerDisplayTag = local_centering_tag_for_display(recenterAnchorMode, anchorFrame);
else
    centerTag = 'rawOrigin';
    centerDisplayTag = 'raw-origin paths';
end
fileTag = sprintf('%s_%s', runTag, centerTag);

%% Render six individual videos
% Data flow: datasets + frame-wise plotting functions -> individual MP4 files.
savedVideoPaths = {};
if saveIndividualVideos
    if showProgressPrint
        fprintf('[toy_PE_debug_videos_displacement] Rendering individual videos...\n');
    end
    pathVideoNames = {'paths_observed', 'paths_predicted', 'paths_pe'};
    rdmVideoNames = {'rdm_observed', 'rdm_predicted', 'rdm_pe'};

    for iData = 1:3
        pathOutBase = fullfile(videoOutputDir, ...
            sprintf('toy_PE_debug_%s_%s_%s', pathVideoNames{iData}, subjectLabel, fileTag));
        pathTitle = sprintf('%s paths (dot1)', datasetLabels{iData});
        pathVideoPath = local_write_paths_video( ...
            datasetsPaths{iData}, trialColors, frameIndices, timeSeconds, ...
            pathAxisLimits, pathOutBase, pathTitle, figureVisibility, ...
            videoFrameRate, videoQuality, showProgressPrint, ...
            progressPrintEveryFrames, sprintf('%s paths', lower(datasetLabels{iData})));
        savedVideoPaths{end + 1} = pathVideoPath; %#ok<AGROW>

        rdmOutBase = fullfile(videoOutputDir, ...
            sprintf('toy_PE_debug_%s_%s_%s', rdmVideoNames{iData}, subjectLabel, fileTag));
        rdmTitle = sprintf('%s RDM', datasetLabels{iData});
        rdmVideoPath = local_write_rdm_video( ...
            datasetsRdm{iData}, timeSeconds, rdmColorLimits, rdmOutBase, rdmTitle, ...
            figureVisibility, videoFrameRate, videoQuality, rdmColormapName, ...
            pathDistMeasure, showProgressPrint, progressPrintEveryFrames, ...
            sprintf('%s rdm', lower(datasetLabels{iData})));
        savedVideoPaths{end + 1} = rdmVideoPath; %#ok<AGROW>
    end
end

%% Render combined 2x3 panel video
% Data flow: all datasets -> synchronized panel frame updates -> one MP4 file.
if savePanelVideo2x3
    if showProgressPrint
        fprintf('[toy_PE_debug_videos_displacement] Rendering combined 2x3 panel video...\n');
    end
    panelOutBase = fullfile(videoOutputDir, ...
        sprintf('toy_PE_debug_panel_%s_%s', subjectLabel, fileTag));
    panelVideoPath = local_write_panel_video( ...
        datasetsPaths, datasetsRdm, datasetLabels, trialColors, frameIndices, ...
        timeSeconds, pathAxisLimits, rdmColorLimits, panelOutBase, ...
        figureVisibility, videoFrameRate, videoQuality, rdmColormapName, ...
        participantNumber, cutFrame, enforcePostDevCut, centerDisplayTag, ...
        rdmSourceMode, showProgressPrint, progressPrintEveryFrames, 'panel 2x3');
    savedVideoPaths{end + 1} = panelVideoPath; %#ok<AGROW>
end

%% Render combined 3x3 panel video with RDM rate row
% Data flow: path + RDM + RDM-rate series -> synchronized 3x3 panel video.
if savePanelVideo3x3Rate
    if showProgressPrint
        fprintf('[toy_PE_debug_videos_displacement] Rendering combined 3x3 rate panel video...\n');
    end
    panelRateOutBase = fullfile(videoOutputDir, ...
        sprintf('toy_PE_debug_panel3x3_rate_%s_%s', subjectLabel, fileTag));
    panelRateVideoPath = local_write_panel_video_with_rate( ...
        datasetsPaths, datasetsRdm, rdmChangeRates, rdmRateYLimits, datasetLabels, ...
        trialColors, frameIndices, timeSeconds, pathAxisLimits, rdmColorLimits, ...
        panelRateOutBase, figureVisibility, videoFrameRate, videoQuality, ...
        rdmColormapName, participantNumber, cutFrame, enforcePostDevCut, ...
        centerDisplayTag, rdmSourceMode, rdmChangeMetric, showProgressPrint, ...
        progressPrintEveryFrames, 'panel 3x3 rate');
    savedVideoPaths{end + 1} = panelRateVideoPath; %#ok<AGROW>
end

%% Report saved outputs
% Data flow: accumulated output path list -> compact terminal summary.
fprintf('\nSaved video outputs (%d):\n', numel(savedVideoPaths));
for iPath = 1:numel(savedVideoPaths)
    fprintf('  %s\n', savedVideoPaths{iPath});
end

%% Local helpers
function outPath = local_write_paths_video( ...
        paths, trialColors, frameIndices, timeSeconds, axisLimits, outBasePath, ...
        titleBase, figVisibility, frameRate, quality, showProgressPrint, ...
        progressPrintEveryFrames, progressLabel)
% LOCAL_WRITE_PATHS_VIDEO Render one paths-only video.
%
% Inputs:
%   paths        : trials x 2 x time path tensor.
%   trialColors  : nTrials x 3 RGB matrix.
%   frameIndices : vector of source frame indices to render.
%   timeSeconds  : vector of timestamps aligned to frameIndices.
%   axisLimits   : [xmin xmax ymin ymax] used for fixed limits.
%   outBasePath  : output path without extension.
%   titleBase    : panel title text.
%   figVisibility: 'on'/'off'.
%   frameRate    : video frame rate.
%   quality      : quality for MPEG-4 when supported.
%   showProgressPrint      : true to print progress updates.
%   progressPrintEveryFrames: print cadence in frames.
%   progressLabel          : short label for progress lines.
%
% Output:
%   outPath : written video path (extension may be mp4 or avi).

[writerObj, outPath] = local_open_video_writer(outBasePath, frameRate, quality);
fig = figure('Visible', figVisibility, 'Color', 'w', ...
    'Name', titleBase, 'NumberTitle', 'off', 'Position', [100, 100, 950, 760]);
ax = axes('Parent', fig);
nFrames = numel(frameIndices);
if showProgressPrint
    fprintf('[toy_PE_debug_videos_displacement] %s -> %s (%d frames)\n', ...
        progressLabel, outPath, nFrames);
end

for iFrame = 1:nFrames
    currFrame = frameIndices(iFrame);
    cla(ax);
    local_plot_paths_frame(ax, paths, currFrame, trialColors, axisLimits);
    title(ax, sprintf('%s | frame %d | t=%.3fs', titleBase, currFrame, timeSeconds(iFrame)));
    frameImage = getframe(fig);
    writeVideo(writerObj, frameImage);
    local_print_video_progress(progressLabel, iFrame, nFrames, ...
        progressPrintEveryFrames, showProgressPrint);
end

close(writerObj);
close(fig);
end

function outPath = local_write_rdm_video( ...
        rdms, timeSeconds, caxisLimits, outBasePath, titleBase, ...
        figVisibility, frameRate, quality, rdmColormapName, distMeasureLabel, ...
        showProgressPrint, progressPrintEveryFrames, progressLabel)
% LOCAL_WRITE_RDM_VIDEO Render one RDM-only video.
%
% Inputs:
%   rdms            : nTrials x nTrials x nFrames cube.
%   timeSeconds     : 1 x nFrames timestamps.
%   caxisLimits     : [min max] dissimilarity limits.
%   outBasePath     : output path without extension.
%   titleBase       : panel title text.
%   figVisibility   : 'on'/'off'.
%   frameRate       : video frame rate.
%   quality         : quality for MPEG-4 when supported.
%   rdmColormapName : MATLAB colormap name string.
%   distMeasureLabel: distance metric label used in colorbar title.
%   showProgressPrint      : true to print progress updates.
%   progressPrintEveryFrames: print cadence in frames.
%   progressLabel          : short label for progress lines.
%
% Output:
%   outPath : written video path (extension may be mp4 or avi).

[writerObj, outPath] = local_open_video_writer(outBasePath, frameRate, quality);
fig = figure('Visible', figVisibility, 'Color', 'w', ...
    'Name', titleBase, 'NumberTitle', 'off', 'Position', [120, 120, 900, 760]);
ax = axes('Parent', fig);
imagesc(ax, rdms(:, :, 1));
set(ax, 'YDir', 'normal');
axis(ax, 'image');
colormap(ax, rdmColormapName);
caxis(ax, caxisLimits);
cb = colorbar(ax);
cb.Label.String = sprintf('Dissimilarity (%s)', distMeasureLabel);
xlabel(ax, 'Trial index');
ylabel(ax, 'Trial index');
nFrames = size(rdms, 3);
if showProgressPrint
    fprintf('[toy_PE_debug_videos_displacement] %s -> %s (%d frames)\n', ...
        progressLabel, outPath, nFrames);
end

for iFrame = 1:nFrames
    set(get(ax, 'Children'), 'CData', rdms(:, :, iFrame));
    title(ax, sprintf('%s | t=%.3fs', titleBase, timeSeconds(iFrame)));
    frameImage = getframe(fig);
    writeVideo(writerObj, frameImage);
    local_print_video_progress(progressLabel, iFrame, nFrames, ...
        progressPrintEveryFrames, showProgressPrint);
end

close(writerObj);
close(fig);
end

function outPath = local_write_panel_video( ...
        datasetsPaths, datasetsRdm, datasetLabels, trialColors, frameIndices, ...
        timeSeconds, pathAxisLimits, rdmColorLimits, outBasePath, ...
        figVisibility, frameRate, quality, rdmColormapName, ...
        participantNumber, cutFrame, enforcePostDevCut, centerDisplayTag, ...
        rdmSourceMode, showProgressPrint, progressPrintEveryFrames, progressLabel)
% LOCAL_WRITE_PANEL_VIDEO Render synchronized 2x3 panel video.
%
% Panel layout:
%   Row 1: observed/predicted/PE paths
%   Row 2: observed/predicted/PE RDMs
%
% Inputs:
%   datasetsPaths      : 1x3 cell of path tensors.
%   datasetsRdm        : 1x3 cell of RDM cubes.
%   datasetLabels      : 1x3 labels.
%   trialColors        : nTrials x 3 RGB matrix.
%   frameIndices       : source frame indices.
%   timeSeconds        : timestamps aligned to frameIndices.
%   pathAxisLimits     : [xmin xmax ymin ymax].
%   rdmColorLimits     : [min max] for all RDM panels.
%   outBasePath        : output path without extension.
%   figVisibility      : 'on'/'off'.
%   frameRate/quality  : writer options.
%   rdmColormapName    : colormap name.
%   participantNumber  : subject id.
%   cutFrame           : post-dev frame index if applicable.
%   enforcePostDevCut  : boolean for title annotation.
%   centerDisplayTag   : text label for path recentering mode in title.
%   rdmSourceMode      : 'rawpostdev' or 'sameaspaths'.
%   showProgressPrint      : true to print progress updates.
%   progressPrintEveryFrames: print cadence in frames.
%   progressLabel          : short label for progress lines.
%
% Output:
%   outPath : written video path.

[writerObj, outPath] = local_open_video_writer(outBasePath, frameRate, quality);
fig = figure('Visible', figVisibility, 'Color', 'w', ...
    'Name', 'Toy PE debug panel video', 'NumberTitle', 'off', ...
    'Position', [60, 60, 1700, 980]);
tl = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nFrames = numel(frameIndices);
if showProgressPrint
    fprintf('[toy_PE_debug_videos_displacement] %s -> %s (%d frames)\n', ...
        progressLabel, outPath, nFrames);
end

pathAxes = gobjects(1, 3);
rdmAxes = gobjects(1, 3);
rdmImg = gobjects(1, 3);

for iCol = 1:3
    pathAxes(iCol) = nexttile(tl, iCol);
    rdmAxes(iCol) = nexttile(tl, 3 + iCol);

    % Initialize RDM panels once, then update CData each frame.
    rdmImg(iCol) = imagesc(rdmAxes(iCol), datasetsRdm{iCol}(:, :, 1));
    set(rdmAxes(iCol), 'YDir', 'normal');
    axis(rdmAxes(iCol), 'image');
    colormap(rdmAxes(iCol), rdmColormapName);
    caxis(rdmAxes(iCol), rdmColorLimits);
    colorbar(rdmAxes(iCol));
    xlabel(rdmAxes(iCol), 'Trial');
    ylabel(rdmAxes(iCol), 'Trial');
end

for iFrame = 1:nFrames
    currFrame = frameIndices(iFrame);
    for iCol = 1:3
        cla(pathAxes(iCol));
        local_plot_paths_frame(pathAxes(iCol), datasetsPaths{iCol}, currFrame, ...
            trialColors, pathAxisLimits);
        title(pathAxes(iCol), sprintf('%s paths', datasetLabels{iCol}));
        set(rdmImg(iCol), 'CData', datasetsRdm{iCol}(:, :, iFrame));
        title(rdmAxes(iCol), sprintf('%s RDM', datasetLabels{iCol}));
    end

    if enforcePostDevCut
        devTag = sprintf('postDev frame %d', cutFrame);
    else
        devTag = 'full trial';
    end
    title(tl, sprintf('Sub%02d | %s | %s | RDM source: %s | frame %d | t=%.3fs', ...
        participantNumber, devTag, centerDisplayTag, rdmSourceMode, currFrame, ...
        timeSeconds(iFrame)));

    frameImage = getframe(fig);
    writeVideo(writerObj, frameImage);
    local_print_video_progress(progressLabel, iFrame, nFrames, ...
        progressPrintEveryFrames, showProgressPrint);
end

close(writerObj);
close(fig);
end

function outPath = local_write_panel_video_with_rate( ...
        datasetsPaths, datasetsRdm, rdmChangeRates, rdmRateYLimits, datasetLabels, ...
        trialColors, frameIndices, timeSeconds, pathAxisLimits, rdmColorLimits, ...
        outBasePath, figVisibility, frameRate, quality, rdmColormapName, ...
        participantNumber, cutFrame, enforcePostDevCut, centerDisplayTag, ...
        rdmSourceMode, rdmChangeMetric, showProgressPrint, ...
        progressPrintEveryFrames, progressLabel)
% LOCAL_WRITE_PANEL_VIDEO_WITH_RATE Render synchronized 3x3 panel video.
%
% Panel layout:
%   Row 1: observed/predicted/PE paths
%   Row 2: observed/predicted/PE RDMs
%   Row 3: observed/predicted/PE RDM change-rate curves
%
% Inputs:
%   datasetsPaths      : 1x3 cell of path tensors.
%   datasetsRdm        : 1x3 cell of RDM cubes.
%   rdmChangeRates     : 1x3 cell of per-frame RDM change rates.
%   rdmRateYLimits     : [min max] limits shared by all rate panels.
%   datasetLabels      : 1x3 labels.
%   trialColors        : nTrials x 3 RGB matrix.
%   frameIndices       : source frame indices.
%   timeSeconds        : timestamps aligned to frameIndices.
%   pathAxisLimits     : [xmin xmax ymin ymax].
%   rdmColorLimits     : [min max] for all RDM panels.
%   outBasePath        : output path without extension.
%   figVisibility      : 'on'/'off'.
%   frameRate/quality  : writer options.
%   rdmColormapName    : colormap name.
%   participantNumber  : subject id.
%   cutFrame           : post-dev frame index if applicable.
%   enforcePostDevCut  : boolean for title annotation.
%   centerDisplayTag   : text label for path recentering mode in title.
%   rdmSourceMode      : 'rawpostdev' or 'sameaspaths'.
%   rdmChangeMetric    : change-rate metric label.
%   showProgressPrint      : true to print progress updates.
%   progressPrintEveryFrames: print cadence in frames.
%   progressLabel          : short label for progress lines.
%
% Output:
%   outPath : written video path.

[writerObj, outPath] = local_open_video_writer(outBasePath, frameRate, quality);
fig = figure('Visible', figVisibility, 'Color', 'w', ...
    'Name', 'Toy PE debug panel video (3x3 rate)', 'NumberTitle', 'off', ...
    'Position', [40, 40, 1780, 1280]);
tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nFrames = numel(frameIndices);
if showProgressPrint
    fprintf('[toy_PE_debug_videos_displacement] %s -> %s (%d frames)\n', ...
        progressLabel, outPath, nFrames);
end

pathAxes = gobjects(1, 3);
rdmAxes = gobjects(1, 3);
rateAxes = gobjects(1, 3);
rdmImg = gobjects(1, 3);
rateCursor = gobjects(1, 3);
rateMarker = gobjects(1, 3);

for iCol = 1:3
    pathAxes(iCol) = nexttile(tl, iCol);
    rdmAxes(iCol) = nexttile(tl, 3 + iCol);
    rateAxes(iCol) = nexttile(tl, 6 + iCol);

    % Initialize RDM panels once, then update CData every frame.
    rdmImg(iCol) = imagesc(rdmAxes(iCol), datasetsRdm{iCol}(:, :, 1));
    set(rdmAxes(iCol), 'YDir', 'normal');
    axis(rdmAxes(iCol), 'image');
    colormap(rdmAxes(iCol), rdmColormapName);
    caxis(rdmAxes(iCol), rdmColorLimits);
    colorbar(rdmAxes(iCol));
    xlabel(rdmAxes(iCol), 'Trial');
    ylabel(rdmAxes(iCol), 'Trial');
    title(rdmAxes(iCol), sprintf('%s RDM', datasetLabels{iCol}));

    % Initialize rate panel with full curve + moving cursor.
    plot(rateAxes(iCol), timeSeconds, rdmChangeRates{iCol}, ...
        '-', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.3);
    hold(rateAxes(iCol), 'on');
    rateCursor(iCol) = line(rateAxes(iCol), ...
        [timeSeconds(1), timeSeconds(1)], rdmRateYLimits, ...
        'Color', [0.85 0.10 0.10], 'LineWidth', 1.2, 'LineStyle', '--');
    initialRate = local_safe_rate_value(rdmChangeRates{iCol}, 1, rdmRateYLimits(1));
    rateMarker(iCol) = plot(rateAxes(iCol), timeSeconds(1), initialRate, ...
        'o', 'Color', [0.85 0.10 0.10], ...
        'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerSize', 4);
    hold(rateAxes(iCol), 'off');
    grid(rateAxes(iCol), 'on');
    box(rateAxes(iCol), 'on');
    xlabel(rateAxes(iCol), 'Time (s)');
    ylabel(rateAxes(iCol), 'RDM change rate');
    xlim(rateAxes(iCol), [timeSeconds(1), timeSeconds(end)]);
    ylim(rateAxes(iCol), rdmRateYLimits);
    title(rateAxes(iCol), sprintf('%s RDM change', datasetLabels{iCol}));
end

for iFrame = 1:nFrames
    currFrame = frameIndices(iFrame);
    currTime = timeSeconds(iFrame);

    for iCol = 1:3
        cla(pathAxes(iCol));
        local_plot_paths_frame(pathAxes(iCol), datasetsPaths{iCol}, currFrame, ...
            trialColors, pathAxisLimits);
        title(pathAxes(iCol), sprintf('%s paths', datasetLabels{iCol}));

        set(rdmImg(iCol), 'CData', datasetsRdm{iCol}(:, :, iFrame));
        set(rateCursor(iCol), 'XData', [currTime, currTime], 'YData', rdmRateYLimits);
        currRate = local_safe_rate_value(rdmChangeRates{iCol}, iFrame, rdmRateYLimits(1));
        set(rateMarker(iCol), 'XData', currTime, 'YData', currRate);
    end

    if enforcePostDevCut
        devTag = sprintf('postDev frame %d', cutFrame);
    else
        devTag = 'full trial';
    end
    title(tl, sprintf(['Sub%02d | %s | %s | RDM source: %s | metric: %s | ', ...
        'frame %d | t=%.3fs'], ...
        participantNumber, devTag, centerDisplayTag, rdmSourceMode, rdmChangeMetric, ...
        currFrame, currTime));

    frameImage = getframe(fig);
    writeVideo(writerObj, frameImage);
    local_print_video_progress(progressLabel, iFrame, nFrames, ...
        progressPrintEveryFrames, showProgressPrint);
end

close(writerObj);
close(fig);
end

function local_plot_paths_frame(ax, paths, frameIdx, trialColors, axisLimits)
% LOCAL_PLOT_PATHS_FRAME Plot cumulative paths up to one frame.
%
% Inputs:
%   ax         : target axes handle.
%   paths      : trials x 2 x time tensor.
%   frameIdx   : scalar frame index to draw up to.
%   trialColors: nTrials x 3 RGB matrix.
%   axisLimits : [xmin xmax ymin ymax].
%
% Data flow:
%   per-trial samples 1:frameIdx -> colored lines + current sample marker.

hold(ax, 'on');
nTrials = size(paths, 1);
for iTrial = 1:nTrials
    xVals = squeeze(paths(iTrial, 1, 1:frameIdx));
    yVals = squeeze(paths(iTrial, 2, 1:frameIdx));
    plot(ax, xVals, yVals, '-', 'Color', trialColors(iTrial, :), 'LineWidth', 1.0);
    plot(ax, xVals(end), yVals(end), 'o', ...
        'Color', trialColors(iTrial, :), 'MarkerFaceColor', trialColors(iTrial, :), ...
        'MarkerSize', 3);
end
xlabel(ax, 'X (visual degrees)');
ylabel(ax, 'Y (visual degrees)');
grid(ax, 'on');
box(ax, 'on');
axis(ax, 'equal');
xlim(ax, axisLimits(1:2));
ylim(ax, axisLimits(3:4));
hold(ax, 'off');
end

function rdms = local_compute_frame_rdms(paths, frameIndices, distMeasure)
% LOCAL_COMPUTE_FRAME_RDMS Compute trial RDM per selected frame.
%
% Inputs:
%   paths       : trials x 2 x time tensor.
%   frameIndices: frame indices to evaluate.
%   distMeasure : pdist metric string.
%
% Output:
%   rdms : nTrials x nTrials x nFrames dissimilarity cube.
%
% Data flow:
%   trial coordinates at each frame -> pdist -> squareform RDM.

nTrials = size(paths, 1);
nFrames = numel(frameIndices);
rdms = zeros(nTrials, nTrials, nFrames);
for iFrame = 1:nFrames
    currFrame = frameIndices(iFrame);
    frameSamples = squeeze(paths(:, :, currFrame));
    frameDistances = pdist(frameSamples, distMeasure);
    rdms(:, :, iFrame) = squareform(frameDistances);
end
end

function limits = local_compute_rdm_limits(datasetsRdm)
% LOCAL_COMPUTE_RDM_LIMITS Compute shared colorbar limits for RDM videos.
%
% Input:
%   datasetsRdm : cell array of RDM cubes.
%
% Output:
%   limits : [min max] dissimilarity range.
%
% Data flow:
%   all RDM values -> global min/max -> robust fallback if flat.

allVals = [];
for iData = 1:numel(datasetsRdm)
    allVals = [allVals; datasetsRdm{iData}(:)]; %#ok<AGROW>
end
minVal = min(allVals);
maxVal = max(allVals);
if ~(isfinite(minVal) && isfinite(maxVal))
    minVal = 0;
    maxVal = 1;
end
if abs(maxVal - minVal) < eps
    maxVal = minVal + 1;
end
limits = [minVal, maxVal];
end

function rdmChangeRates = local_compute_rdm_change_rates(datasetsRdm, metricName)
% LOCAL_COMPUTE_RDM_CHANGE_RATES Compute per-frame RDM change-rate vectors.
%
% Inputs:
%   datasetsRdm : cell array of RDM cubes (nTrials x nTrials x nFrames).
%   metricName  : rate metric label (currently 'meanabsoffdiag').
%
% Output:
%   rdmChangeRates : cell array, each entry is [1 x nFrames] rate vector.
%
% Data flow:
%   RDM cube -> frame-to-frame deltas -> scalar rate at each frame.

rdmChangeRates = cell(size(datasetsRdm));
for iData = 1:numel(datasetsRdm)
    currRdm = datasetsRdm{iData};
    nTrials = size(currRdm, 1);
    nFrames = size(currRdm, 3);
    offDiagMask = ~eye(nTrials);

    currRate = zeros(1, nFrames);
    for iFrame = 2:nFrames
        deltaRdm = currRdm(:, :, iFrame) - currRdm(:, :, iFrame - 1);
        if strcmp(metricName, 'meanabsoffdiag')
            currRate(iFrame) = mean(abs(deltaRdm(offDiagMask)));
        else
            error('Unsupported RDM change metric: %s', metricName);
        end
    end
    rdmChangeRates{iData} = currRate;
end
end

function limits = local_compute_rate_axis_limits(rateSeriesCell)
% LOCAL_COMPUTE_RATE_AXIS_LIMITS Compute shared y-limits for rate plots.
%
% Input:
%   rateSeriesCell : cell array of per-frame rate vectors.
%
% Output:
%   limits : [ymin ymax] limits shared across all rate panels.
%
% Data flow:
%   all rate values -> min/max -> padded axis limits.

allRates = [];
for iSeries = 1:numel(rateSeriesCell)
    allRates = [allRates, rateSeriesCell{iSeries}]; %#ok<AGROW>
end
allRates = allRates(isfinite(allRates));
if isempty(allRates)
    limits = [0, 1];
    return;
end

yMin = min(allRates);
yMax = max(allRates);
if yMax <= yMin
    yMax = yMin + 1;
end

if yMin >= 0
    yMinPlot = 0;
else
    yMinPlot = yMin;
end
yPad = 0.08 * (yMax - yMinPlot);
limits = [yMinPlot - yPad, yMax + yPad];
end

function valueOut = local_safe_rate_value(rateSeries, frameIdx, fallbackValue)
% LOCAL_SAFE_RATE_VALUE Return finite rate value with fallback guard.
%
% Inputs:
%   rateSeries    : per-frame rate vector.
%   frameIdx      : current frame index.
%   fallbackValue : value used if rate is NaN/Inf.
%
% Output:
%   valueOut : finite scalar value for marker plotting.
valueOut = rateSeries(frameIdx);
if ~isfinite(valueOut)
    valueOut = fallbackValue;
end
end

function local_print_video_progress( ...
        progressLabel, currFrame, totalFrames, printEveryFrames, showProgressPrint)
% LOCAL_PRINT_VIDEO_PROGRESS Emit periodic progress lines for long renders.
%
% Inputs:
%   progressLabel      : short label identifying the current video.
%   currFrame          : current frame index (1-based).
%   totalFrames        : total number of frames in the render loop.
%   printEveryFrames   : cadence in frames for progress updates.
%   showProgressPrint  : true to enable progress output.
%
% Data flow:
%   frame counters -> threshold check -> concise status line.
if ~showProgressPrint
    return;
end
if currFrame == 1 || currFrame == totalFrames || mod(currFrame, printEveryFrames) == 0
    fprintf('[toy_PE_debug_videos_displacement] %s progress: %d/%d (%.1f%%)\n', ...
        progressLabel, currFrame, totalFrames, 100 * currFrame / totalFrames);
end
end

function limits = local_compute_global_axis_limits(datasetsPaths)
% LOCAL_COMPUTE_GLOBAL_AXIS_LIMITS Compute fixed xy limits for path videos.
%
% Input:
%   datasetsPaths : cell array of paths tensors (trials x 2 x time).
%
% Output:
%   limits : [xmin xmax ymin ymax] padded and shared across datasets.
%
% Data flow:
%   all coordinates -> global min/max -> padding for stable framing.

xAll = [];
yAll = [];
for iData = 1:numel(datasetsPaths)
    curr = datasetsPaths{iData};
    xAll = [xAll; reshape(curr(:, 1, :), [], 1)]; %#ok<AGROW>
    yAll = [yAll; reshape(curr(:, 2, :), [], 1)]; %#ok<AGROW>
end
xMin = min(xAll);
xMax = max(xAll);
yMin = min(yAll);
yMax = max(yAll);
xSpan = max(xMax - xMin, eps);
ySpan = max(yMax - yMin, eps);
padX = 0.05 * xSpan;
padY = 0.05 * ySpan;
limits = [xMin - padX, xMax + padX, yMin - padY, yMax + padY];
end

function colors = local_trial_colors(nTrials)
% LOCAL_TRIAL_COLORS Build stable per-trial colors for path rendering.
%
% Input:
%   nTrials : number of trials.
%
% Output:
%   colors  : nTrials x 3 RGB matrix.
if nTrials <= 0
    colors = zeros(0, 3);
    return;
end
colors = hsv(nTrials);
end

function pathsCentered = local_recenter_paths_at_frame(paths, anchorFrame)
% LOCAL_RECENTER_PATHS_AT_FRAME Shift each trial so anchorFrame is [0, 0].
%
% Usage example:
%   pathsCentered = local_recenter_paths_at_frame(paths, 120);
%
% Inputs:
%   paths       : trials x 2 x time path tensor.
%   anchorFrame : 1-based frame index used as trial-wise origin.
%
% Output:
%   pathsCentered : recentered path tensor.
%
% Data flow:
%   per-trial anchor at anchorFrame -> subtract anchor from all frames.
if isempty(paths)
    pathsCentered = paths;
    return;
end
nTime = size(paths, 3);
if ~isscalar(anchorFrame) || floor(anchorFrame) ~= anchorFrame || ...
        anchorFrame < 1 || anchorFrame > nTime
    error('anchorFrame must be an integer in [1, %d].', nTime);
end
anchor = paths(:, :, anchorFrame);
pathsCentered = paths - repmat(anchor, 1, 1, nTime);
end

function local_assert_anchor_at_origin(datasets, anchorFrame, tolerance)
% LOCAL_ASSERT_ANCHOR_AT_ORIGIN Validate origin at anchorFrame after recentering.
%
% Inputs:
%   datasets    : cell array of recentered path tensors.
%   anchorFrame : frame that should be [0, 0] for every trial.
%   tolerance   : max absolute allowed residual.
for iData = 1:numel(datasets)
    originVals = datasets{iData}(:, :, anchorFrame);
    maxAbsOrigin = max(abs(originVals(:)));
    assert(maxAbsOrigin <= tolerance, ...
        ['Anchor recentering failed for dataset %d at frame %d ', ...
        '(max abs %.3g > %.3g).'], ...
        iData, anchorFrame, maxAbsOrigin, tolerance);
end
end

function [writerObj, outPath] = local_open_video_writer(outBasePath, frameRate, quality)
% LOCAL_OPEN_VIDEO_WRITER Open video writer with MPEG-4 fallback handling.
%
% Inputs:
%   outBasePath : output path without extension.
%   frameRate   : desired frames per second.
%   quality     : quality value for MPEG-4 when supported.
%
% Outputs:
%   writerObj : opened VideoWriter object.
%   outPath   : chosen output file path (.mp4 preferred, .avi fallback).
%
% Data flow:
%   preferred MP4 writer -> fallback AVI writer if unavailable.
[outDir, outName, ~] = fileparts(outBasePath);
local_ensure_dir(outDir);
mp4Path = fullfile(outDir, [outName '.mp4']);
aviPath = fullfile(outDir, [outName '.avi']);

try
    writerObj = VideoWriter(mp4Path, 'MPEG-4');
    outPath = mp4Path;
catch
    warning('MPEG-4 writer unavailable; falling back to Motion JPEG AVI for %s', outName);
    writerObj = VideoWriter(aviPath, 'Motion JPEG AVI');
    outPath = aviPath;
end

writerObj.FrameRate = frameRate;
if isprop(writerObj, 'Quality')
    writerObj.Quality = quality;
end
open(writerObj);
end

function cutFrame = local_compute_cut_frame(nTime, deviantOnset)
% LOCAL_COMPUTE_CUT_FRAME Convert deviantOnset fraction to a 1-based frame.
%
% Inputs:
%   nTime       : number of frames in full trajectory.
%   deviantOnset: onset fraction in [0, 1].
%
% Output:
%   cutFrame : onset frame index clamped to [1, nTime].
if ~isscalar(nTime) || floor(nTime) ~= nTime || nTime < 1
    error('nTime must be a positive integer.');
end
if ~isscalar(deviantOnset) || deviantOnset < 0 || deviantOnset > 1
    error('deviantOnset must be a scalar between 0 and 1.');
end
cutFrame = round(deviantOnset * nTime);
cutFrame = max(1, min(nTime, cutFrame));
end

function anchorFrame = local_resolve_anchor_frame(nTime, cutFrame, recenterAnchorMode)
% LOCAL_RESOLVE_ANCHOR_FRAME Select the recenter frame for path visualization.
%
% Inputs:
%   nTime             : full-trial length.
%   cutFrame          : post-deviance onset frame.
%   recenterAnchorMode: 'predeviantend' or 'deviantsample'.
%
% Output:
%   anchorFrame : frame used as [0, 0] after recentering.
%
% Notes:
%   - 'predeviantend' is displacement-aware: anchor at cutFrame-1 so the
%     first post-deviant sample preserves the observed jump magnitude.
if ~isscalar(cutFrame) || floor(cutFrame) ~= cutFrame || cutFrame < 1 || cutFrame > nTime
    error('cutFrame must be an integer in [1, %d].', nTime);
end
switch recenterAnchorMode
    case 'predeviantend'
        anchorFrame = max(cutFrame - 1, 1);
    case 'deviantsample'
        anchorFrame = cutFrame;
    otherwise
        error('Unsupported recenterAnchorMode: %s', recenterAnchorMode);
end
end

function pathsCut = local_cut_from_frame(paths, cutFrame)
% LOCAL_CUT_FROM_FRAME Slice tensor from cutFrame to end.
%
% Inputs:
%   paths    : trials x features x time tensor.
%   cutFrame : 1-based start frame for slice.
%
% Output:
%   pathsCut : trials x features x (time-cutFrame+1) tensor.
if isempty(paths)
    pathsCut = paths;
    return;
end
nTime = size(paths, 3);
if ~isscalar(cutFrame) || floor(cutFrame) ~= cutFrame || cutFrame < 1 || cutFrame > nTime
    error('cutFrame must be an integer in [1, %d].', nTime);
end
pathsCut = paths(:, :, cutFrame:end);
end

function centerTag = local_centering_tag_for_filename(recenterAnchorMode)
% LOCAL_CENTERING_TAG_FOR_FILENAME Build compact recenter tag for filenames.
switch recenterAnchorMode
    case 'predeviantend'
        centerTag = 'preDevAnchor';
    case 'deviantsample'
        centerTag = 'devSampleAnchor';
    otherwise
        error('Unsupported recenterAnchorMode: %s', recenterAnchorMode);
end
end

function centerDisplayTag = local_centering_tag_for_display(recenterAnchorMode, anchorFrame)
% LOCAL_CENTERING_TAG_FOR_DISPLAY Build verbose recenter tag for titles.
switch recenterAnchorMode
    case 'predeviantend'
        centerDisplayTag = sprintf('pre-dev anchor (frame %d)', anchorFrame);
    case 'deviantsample'
        centerDisplayTag = sprintf('deviant-sample anchor (frame %d)', anchorFrame);
    otherwise
        error('Unsupported recenterAnchorMode: %s', recenterAnchorMode);
end
end

function local_ensure_dir(dirPath)
% LOCAL_ENSURE_DIR Create directory if missing (including parent dirs).
%
% Input:
%   dirPath : target directory path.
if isempty(dirPath)
    return;
end
if ~exist(dirPath, 'dir')
    mkdir(dirPath);
end
end
