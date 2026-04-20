%% Plot curvature-limit paths with fixed directionChange = 0
% Script: plot_baseline_curvature_limit_paths.m
%
% Purpose:
%   Visualize how trajectories diverge after the deviance frame when using
%   the four limit values from the curvature windows reported in methods:
%     - -0.8 deg/frame
%     - -0.3755 deg/frame
%     -  0.3755 deg/frame
%     -  0.8 deg/frame
%
%   The visualization keeps:
%     - directionChange = 0 (no extra turn impulse),
%     - one shared initial direction,
%     - one shared pre-deviance trajectory segment drawn as dashed black.
%
%   To keep a common pre-deviance segment (as in the turn-limit figure), the
%   curvature limit values above are applied from deviance onward, while the
%   pre-deviance segment uses one fixed baselineCurvatureDeg reference value.
%
% Usage example (interactive):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper');
%   run('plot_baseline_curvature_limit_paths.m');
%
% Usage example (non-interactive, preferred):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/plot_baseline_curvature_limit_paths.m')"
%
% Usage example (override parameters before run):
%   baselineCurvatureDeg = 0.3755;
%   curvatureLimitDegPerFrame = [-0.8, -0.3755, 0.3755, 0.8];
%   figVisibility = 'off';
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/plot_baseline_curvature_limit_paths.m');
%
% Inputs (optional workspace overrides):
%   fps                      : sampling rate in Hz (default 120)
%   trialDurationSec         : trial duration in seconds (default 2.67)
%   dotSpeedDegPerSec        : dot speed in deg/s (default 3.73)
%   devianceFrame            : 1-based deviance frame index (default 130)
%   initialDirectionDeg      : initial heading in degrees (default 0)
%   baselineCurvatureDeg     : shared pre-deviance curvature in deg/frame
%                              (default 0.3755)
%   curvatureLimitDegPerFrame: curvature values applied from deviance onward
%                              (default [-0.8 -0.3755 0.3755 0.8])
%   directionChangeDeg       : extra turn impulse at deviance (default 0)
%   figVisibility            : 'on' or 'off' (default 'off')
%
% Outputs:
%   - Figure with one shared dashed-black pre-deviance segment and multiple
%     post-deviance curvature-limit continuations.
%   - PNG saved to:
%       paper/figures/baseline_curvature_limit_paths.png
%
% Key assumptions:
%   - Curvature is interpreted as signed heading increment in deg/frame.
%   - Curvature change is applied from onset step (devianceFrame - 1) onward.
%   - No arena-bound or fixation-collision correction is applied here because
%     the goal is geometric comparison of curvature limit effects.

%% Resolve defaults and validate caller overrides
% Data flow: optional workspace inputs -> validated numeric controls.
clearvars -except fps trialDurationSec dotSpeedDegPerSec devianceFrame initialDirectionDeg baselineCurvatureDeg curvatureLimitDegPerFrame directionChangeDeg figVisibility;
close all;

if ~exist('fps', 'var') || isempty(fps)
    fps = 120;
end
if ~exist('trialDurationSec', 'var') || isempty(trialDurationSec)
    trialDurationSec = 2.67;
end
if ~exist('dotSpeedDegPerSec', 'var') || isempty(dotSpeedDegPerSec)
    dotSpeedDegPerSec = 3.73;
end
if ~exist('devianceFrame', 'var') || isempty(devianceFrame)
    devianceFrame = 130;
end
if ~exist('initialDirectionDeg', 'var') || isempty(initialDirectionDeg)
    initialDirectionDeg = 0;
end
if ~exist('baselineCurvatureDeg', 'var') || isempty(baselineCurvatureDeg)
    baselineCurvatureDeg = 0.3755;
end
if ~exist('curvatureLimitDegPerFrame', 'var') || isempty(curvatureLimitDegPerFrame)
    curvatureLimitDegPerFrame = [-0.8, -0.3755, 0.3755, 0.8];
end
if ~exist('directionChangeDeg', 'var') || isempty(directionChangeDeg)
    directionChangeDeg = 0;
end
if ~exist('figVisibility', 'var') || isempty(figVisibility)
    figVisibility = 'off';
end

fps = double(fps);
trialDurationSec = double(trialDurationSec);
dotSpeedDegPerSec = double(dotSpeedDegPerSec);
devianceFrame = round(double(devianceFrame));
initialDirectionDeg = double(initialDirectionDeg);
baselineCurvatureDeg = double(baselineCurvatureDeg);
curvatureLimitDegPerFrame = double(curvatureLimitDegPerFrame(:)');
directionChangeDeg = double(directionChangeDeg);
figVisibility = lower(char(figVisibility));

if fps <= 0
    error('fps must be > 0.');
end
if trialDurationSec <= 0
    error('trialDurationSec must be > 0.');
end
if dotSpeedDegPerSec <= 0
    error('dotSpeedDegPerSec must be > 0.');
end
if ~ismember(figVisibility, {'on', 'off'})
    error('figVisibility must be ''on'' or ''off''.');
end
if numel(curvatureLimitDegPerFrame) < 1
    error('curvatureLimitDegPerFrame must contain at least one value.');
end

framesPerTrial = round(trialDurationSec * fps);
if framesPerTrial < 2
    error('framesPerTrial must be >= 2. Check trialDurationSec and fps.');
end
if devianceFrame < 2 || devianceFrame > framesPerTrial
    error('devianceFrame must be in [2, framesPerTrial].');
end

dotSpeedDegPerFrame = dotSpeedDegPerSec / fps;
nSteps = framesPerTrial - 1;
devianceStep = devianceFrame - 1;

%% Build shared nondeviant reference path
% Data flow: fixed initial direction + zero turn + fixed baseline curvature -> baseline XY.
turnNondeviantDeg = zeros(nSteps, 1);
curvatureNondeviantDeg = repmat(baselineCurvatureDeg, nSteps, 1);
directionNondeviantDeg = local_integrate_direction( ...
    initialDirectionDeg, turnNondeviantDeg, curvatureNondeviantDeg);
xyNondeviant = local_integrate_xy(directionNondeviantDeg, dotSpeedDegPerFrame);

%% Build deviant paths for curvature limit values
% Data flow: shared prefix controls + per-variant post-deviance curvature -> shifted deviant XY.
nVariants = numel(curvatureLimitDegPerFrame);
xyDeviants = zeros(framesPerTrial, 2, nVariants);

for iVariant = 1:nVariants
    postCurvatureDeg = curvatureLimitDegPerFrame(iVariant);

    turnDeviantDeg = zeros(nSteps, 1);
    turnDeviantDeg(devianceStep) = directionChangeDeg;

    curvatureDeviantDeg = curvatureNondeviantDeg;
    curvatureDeviantDeg(devianceStep:end) = postCurvatureDeg;

    directionDeviantRawDeg = local_integrate_direction( ...
        initialDirectionDeg, turnDeviantDeg, curvatureDeviantDeg);
    xyDeviantRaw = local_integrate_xy(directionDeviantRawDeg, dotSpeedDegPerFrame);

    % Match the generator-style shared-prefix rule, then translate suffix to
    % preserve continuity at the deviance frame.
    xyDeviant = xyDeviantRaw;
    xyDeviant(1:devianceFrame, :) = xyNondeviant(1:devianceFrame, :);
    suffixRaw = xyDeviantRaw(devianceFrame:end, :);
    suffixShift = xyDeviant(devianceFrame, :) - suffixRaw(1, :);
    xyDeviant(devianceFrame:end, :) = suffixRaw + suffixShift;

    xyDeviants(:, :, iVariant) = xyDeviant;
end

%% Render and save figure
% Data flow: baseline + deviant tensors -> dashed-prefix styled figure -> PNG file.
scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir, 'figures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
outputPng = fullfile(outputDir, 'baseline_curvature_limit_paths.png');

fig = figure('Name', 'Baseline curvature limit paths', ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Visible', figVisibility);
hold on;

preIdx = 1:devianceFrame;
postIdx = devianceFrame:framesPerTrial;

% Draw the shared pre-deviance path as requested (dashed black).
plot(xyNondeviant(preIdx, 1), xyNondeviant(preIdx, 2), ...
    'k--', 'LineWidth', 2.0, ...
    'DisplayName', 'shared pre-deviance path');

pathColors = [ ...
    0.00, 0.45, 0.74; ...
    0.85, 0.33, 0.10; ...
    0.47, 0.67, 0.19; ...
    0.49, 0.18, 0.56];

for iVariant = 1:nVariants
    colorIdx = mod(iVariant - 1, size(pathColors, 1)) + 1;
    curvatureLabel = sprintf('post-onset curvature = %+g deg/frame', curvatureLimitDegPerFrame(iVariant));
    plot(xyDeviants(postIdx, 1, iVariant), xyDeviants(postIdx, 2, iVariant), ...
        '-', 'Color', pathColors(colorIdx, :), 'LineWidth', 1.8, ...
        'DisplayName', curvatureLabel);
end

% Draw nondeviant continuation after variants so it remains visible even
% when a variant exactly overlaps it (e.g., +0.3755 deg/frame).
plot(xyNondeviant(postIdx, 1), xyNondeviant(postIdx, 2), ...
    'k-', 'LineWidth', 2.6, ...
    'DisplayName', 'nondeviant post-deviance (turn = 0)');

plot(xyNondeviant(devianceFrame, 1), xyNondeviant(devianceFrame, 2), ...
    'ko', 'MarkerFaceColor', [1 1 1], 'MarkerSize', 7, ...
    'DisplayName', sprintf('deviance frame %d', devianceFrame));

axis equal;
grid on;
box on;
xlabel('x (deg)');
ylabel('y (deg)');
title(sprintf([ ...
    'Curvature-window limit effects at deviance (directionChange = %.1f deg, pre curvature = %.4f deg/frame)'], ...
    directionChangeDeg, baselineCurvatureDeg));
legend('Location', 'bestoutside');

print(fig, outputPng, '-dpng', '-r300');
fprintf('Saved figure: %s\n', outputPng);

%% Local helpers
function directionsDeg = local_integrate_direction(initialDirectionDeg, turnDeg, curvatureDeg)
% LOCAL_INTEGRATE_DIRECTION Integrate per-step turn+curvature into frame directions.
nFrames = numel(turnDeg) + 1;
directionsDeg = zeros(nFrames, 1);
directionsDeg(1) = initialDirectionDeg;
if nFrames > 1
    directionsDeg(2:end) = initialDirectionDeg + cumsum(turnDeg + curvatureDeg, 1);
end
end

function xy = local_integrate_xy(directionDeg, dotSpeedDegPerFrame)
% LOCAL_INTEGRATE_XY Integrate frame directions into relative x/y trajectory.
nFrames = numel(directionDeg);
xy = zeros(nFrames, 2);
if nFrames > 1
    steps = [cosd(directionDeg(2:end)), sind(directionDeg(2:end))] * dotSpeedDegPerFrame;
    xy(2:end, :) = cumsum(steps, 1);
end
end
