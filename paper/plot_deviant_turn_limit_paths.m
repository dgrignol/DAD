%% Plot deviant-turn limit paths with fixed directionChange = 0
% Script: plot_deviant_turn_limit_paths.m
%
% Purpose:
%   Visualize how one shared baseline trajectory changes after the deviance
%   frame when the deviant signed-turn value is fixed to the four boundary
%   values used in the methods text:
%     - -81 deg
%     - -10 deg
%     -  10 deg
%     -  81 deg
%
%   The visualization keeps:
%     - directionChange = 0 (no extra nondeviant direction change),
%     - one shared initial direction,
%     - one shared baseline curvature,
%     - one shared pre-deviance path segment.
%
%   Post-deviance path construction mirrors the V28 generation logic:
%     1) integrate a raw deviant path with the chosen onset turn,
%     2) force identity with nondeviant trajectory through deviance frame,
%     3) translate the post-deviance suffix to keep positional continuity.
%
% Usage example (interactive):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper');
%   run('plot_deviant_turn_limit_paths.m');
%
% Usage example (non-interactive, preferred):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/plot_deviant_turn_limit_paths.m')"
%
% Usage example (override parameters before run):
%   baselineCurvatureDeg = -0.3755;
%   figVisibility = 'off';
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/plot_deviant_turn_limit_paths.m');
%
% Inputs (optional workspace overrides):
%   fps                 : sampling rate in Hz (default 120)
%   trialDurationSec    : trial duration in seconds (default 2.67)
%   dotSpeedDegPerSec   : dot speed in deg/s (default 3.73)
%   devianceFrame       : 1-based deviance frame index (default 130)
%   initialDirectionDeg : initial heading in degrees (default 0)
%   baselineCurvatureDeg: constant baseline curvature in deg/frame
%                         (default 0.3755, one lower bound in allowed window)
%   directionChangeDeg  : extra direction change added at deviance (default 0)
%   deviantTurnDegValues: deviant signed-turn values in deg (default [-81 -10 10 81])
%   figVisibility       : 'on' or 'off' (default 'off')
%
% Outputs:
%   - Figure with a shared dashed-black pre-deviance segment and post-deviance
%     segments for nondeviant + 4 deviant variants.
%   - PNG saved to:
%       paper/figures/deviant_turn_limit_paths.png
%
% Key assumptions:
%   - Curvature is interpreted as signed heading increment in deg/frame.
%   - Turn perturbation is applied at the onset step (devianceFrame - 1).
%   - No arena-bound or fixation-collision correction is applied here because
%     the goal is geometric comparison of turn limits on a shared path.

%% Resolve defaults and validate caller overrides
% Data flow: optional workspace inputs -> validated numeric controls.
clearvars -except fps trialDurationSec dotSpeedDegPerSec devianceFrame initialDirectionDeg baselineCurvatureDeg directionChangeDeg deviantTurnDegValues figVisibility;
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
if ~exist('directionChangeDeg', 'var') || isempty(directionChangeDeg)
    directionChangeDeg = 0;
end
if ~exist('deviantTurnDegValues', 'var') || isempty(deviantTurnDegValues)
    deviantTurnDegValues = [-81, -10, 10, 81];
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
directionChangeDeg = double(directionChangeDeg);
deviantTurnDegValues = double(deviantTurnDegValues(:)');
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
if numel(deviantTurnDegValues) < 1
    error('deviantTurnDegValues must contain at least one value.');
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
% Data flow: fixed initial direction + zero turn + constant curvature -> baseline XY.
turnNondeviantDeg = zeros(nSteps, 1);
curvatureNondeviantDeg = repmat(baselineCurvatureDeg, nSteps, 1);
directionNondeviantDeg = local_integrate_direction( ...
    initialDirectionDeg, turnNondeviantDeg, curvatureNondeviantDeg);
xyNondeviant = local_integrate_xy(directionNondeviantDeg, dotSpeedDegPerFrame);

%% Build deviant paths for the requested signed-turn values
% Data flow: one turn value at onset -> raw deviant XY -> forced shared prefix + translated suffix.
nVariants = numel(deviantTurnDegValues);
xyDeviants = zeros(framesPerTrial, 2, nVariants);

for iVariant = 1:nVariants
    onsetTurnDeg = directionChangeDeg + deviantTurnDegValues(iVariant);

    turnDeviantDeg = zeros(nSteps, 1);
    turnDeviantDeg(devianceStep) = onsetTurnDeg;
    curvatureDeviantDeg = curvatureNondeviantDeg;

    directionDeviantRawDeg = local_integrate_direction( ...
        initialDirectionDeg, turnDeviantDeg, curvatureDeviantDeg);
    xyDeviantRaw = local_integrate_xy(directionDeviantRawDeg, dotSpeedDegPerFrame);

    % Match the generator's explicit pre-deviance identity rule, then shift
    % the suffix so continuity is preserved at the deviance frame.
    xyDeviant = xyDeviantRaw;
    xyDeviant(1:devianceFrame, :) = xyNondeviant(1:devianceFrame, :);
    suffixRaw = xyDeviantRaw(devianceFrame:end, :);
    suffixShift = xyDeviant(devianceFrame, :) - suffixRaw(1, :);
    xyDeviant(devianceFrame:end, :) = suffixRaw + suffixShift;

    xyDeviants(:, :, iVariant) = xyDeviant;
end

%% Render and save figure
% Data flow: baseline + deviant path tensors -> styled overlay figure -> PNG file.
scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir, 'figures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
outputPng = fullfile(outputDir, 'deviant_turn_limit_paths.png');

fig = figure('Name', 'Deviant turn limit paths', ...
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

% Draw the nondeviant continuation after deviance as reference.
plot(xyNondeviant(postIdx, 1), xyNondeviant(postIdx, 2), ...
    'k-', 'LineWidth', 2.0, ...
    'DisplayName', 'nondeviant post-deviance');

pathColors = [ ...
    0.00, 0.45, 0.74; ...
    0.85, 0.33, 0.10; ...
    0.47, 0.67, 0.19; ...
    0.49, 0.18, 0.56];

for iVariant = 1:nVariants
    colorIdx = mod(iVariant - 1, size(pathColors, 1)) + 1;
    turnLabel = sprintf('deviant turn = %+g deg', deviantTurnDegValues(iVariant));
    % Plot only post-deviance segments because the prefix is shared and
    % already shown as dashed black.
    plot(xyDeviants(postIdx, 1, iVariant), xyDeviants(postIdx, 2, iVariant), ...
        '-', 'Color', pathColors(colorIdx, :), 'LineWidth', 1.8, ...
        'DisplayName', turnLabel);
end

plot(xyNondeviant(devianceFrame, 1), xyNondeviant(devianceFrame, 2), ...
    'ko', 'MarkerFaceColor', [1 1 1], 'MarkerSize', 7, ...
    'DisplayName', sprintf('deviance frame %d', devianceFrame));

axis equal;
grid on;
box on;
xlabel('x (deg)');
ylabel('y (deg)');
title(sprintf([ ...
    'Path divergence at deviance (directionChange = %.1f deg, baseline curvature = %.4f deg/frame)'], ...
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
