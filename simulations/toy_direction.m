% toy_direction.m
%
% Purpose:
%   Build a tiny, fully transparent direction model from toy dot paths and
%   visualize every intermediate step, then compute an RDM using cosine
%   similarity between direction vectors.
%
% Example usage (from simulations/ in MATLAB):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   toy_direction;
%
% Example usage (from anywhere in MATLAB):
%   addpath('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations');
%   toy_direction;
%
% Inputs:
%   - None (toy data are defined inside the script).
%
% Outputs:
%   - Figures showing: raw paths, dx/dy, direction angles, unit vectors,
%     trial-time concatenation order, cosine similarity, and RDM.
%
% Key assumptions:
%   - Trials are shaped as (trials x 2 x time), with x in index 1, y in index 2.
%   - Direction is inferred from frame-to-frame displacement.
%   - Cosine similarity is computed on unit vectors (cos(theta), sin(theta)).

clear variables
close all
clc

%% Define tiny toy paths
% Data flow: manual x/y sequences -> dot1Trials and dot2Trials arrays.
nTrials = 2;
nTime = 5;

% Dot 1 trials: simple right-then-up and up-then-left trajectories.
dot1Trial1 = [0 1 2 2 2; 0 0 0 1 2];   % right then up
dot1Trial2 = [0 0 0 -1 -2; 0 1 2 2 2]; % up then left
dot1Trials = zeros(nTrials, 2, nTime);
dot1Trials(1, :, :) = dot1Trial1;
dot1Trials(2, :, :) = dot1Trial2;

% Dot 2 trials: diagonals to show different directions.
dot2Trial1 = [0 1 2 3 4; 0 -1 -2 -3 -4]; % down-right
dot2Trial2 = [0 -1 -2 -3 -4; 0 -1 -2 -3 -4]; % down-left
dot2Trials = zeros(nTrials, 2, nTime);
dot2Trials(1, :, :) = dot2Trial1;
dot2Trials(2, :, :) = dot2Trial2;

%% Visualize raw paths with timepoint labels
% Data flow: dot paths -> per-trial path plots with labeled timepoints.
figure('Name', 'Toy Paths (Dot 1 and Dot 2)');
for iTrial = 1:nTrials
    subplot(2, nTrials, iTrial);
    plot_trial_path_with_labels(squeeze(dot1Trials(iTrial, :, :)), ...
        sprintf('Dot 1 trial %d', iTrial));
    subplot(2, nTrials, nTrials + iTrial);
    plot_trial_path_with_labels(squeeze(dot2Trials(iTrial, :, :)), ...
        sprintf('Dot 2 trial %d', iTrial));
end

%% Compute frame-to-frame displacement (dx, dy)
% Data flow: positions -> diff over time -> padded to match original length.
dot1Dx = diff(dot1Trials(:, 1, :), 1, 3);
dot1Dy = diff(dot1Trials(:, 2, :), 1, 3);
dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx); % pad with first displacement
dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy); % pad with first displacement

dot2Dx = diff(dot2Trials(:, 1, :), 1, 3);
dot2Dy = diff(dot2Trials(:, 2, :), 1, 3);
dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx); % pad with first displacement
dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy); % pad with first displacement

%% Visualize dx and dy over time
% Data flow: dx/dy tensors -> per-trial time series plots.
figure('Name', 'Displacement Over Time (dx and dy)');
plot_displacement_series(dot1Dx, dot1Dy, 'Dot 1', 1);
plot_displacement_series(dot2Dx, dot2Dy, 'Dot 2', 2);

%% Convert displacement into direction angles
% Data flow: dx/dy -> angle in radians -> angle in degrees for readability.
dot1AngleRad = atan2(dot1Dy, dot1Dx);
dot2AngleRad = atan2(dot2Dy, dot2Dx);
dot1AngleDeg = rad2deg(dot1AngleRad);
dot2AngleDeg = rad2deg(dot2AngleRad);

%% Visualize direction angles over time
% Data flow: angle matrices -> per-trial angle time series plots.
figure('Name', 'Direction Angles Over Time');
plot_angle_series(dot1AngleDeg, 'Dot 1', 1);
plot_angle_series(dot2AngleDeg, 'Dot 2', 2);

%% Convert angles to unit direction vectors
% Data flow: angles -> unit vectors [cos(theta), sin(theta)] per timepoint.
dot1Direction = cat(2, cos(dot1AngleRad), sin(dot1AngleRad));
dot2Direction = cat(2, cos(dot2AngleRad), sin(dot2AngleRad));

%% Visualize unit vectors overlaid on paths
% Data flow: positions + dx/dy -> quiver plot for each trial.
figure('Name', 'Direction Vectors Overlaid on Paths');
for iTrial = 1:nTrials
    subplot(2, nTrials, iTrial);
    plot_path_with_vectors(squeeze(dot1Trials(iTrial, :, :)), ...
        squeeze(dot1Dx(iTrial, :, :)), squeeze(dot1Dy(iTrial, :, :)), ...
        sprintf('Dot 1 trial %d', iTrial));
    subplot(2, nTrials, nTrials + iTrial);
    plot_path_with_vectors(squeeze(dot2Trials(iTrial, :, :)), ...
        squeeze(dot2Dx(iTrial, :, :)), squeeze(dot2Dy(iTrial, :, :)), ...
        sprintf('Dot 2 trial %d', iTrial));
end

%% Concatenate trials and time into a sample-by-feature matrix
% Data flow: (trial x 2 x time) -> (samples x 2) with samples ordered by trial then time.
dot1Samples = reshape(permute(dot1Direction, [1 3 2]), nTrials * nTime, 2);
dot2Samples = reshape(permute(dot2Direction, [1 3 2]), nTrials * nTime, 2);

% Build labels so each sample can be traced back to (trial, time).
sampleLabels = cell(nTrials * nTime, 1);
labelIdx = 1;
for iTrial = 1:nTrials
    for iTime = 1:nTime
        sampleLabels{labelIdx} = sprintf('tr%d_t%d', iTrial, iTime);
        labelIdx = labelIdx + 1;
    end
end

%% Visualize concatenation order
% Data flow: sample index -> trial x time grid for interpretability.
sampleIndexGrid = reshape(1:(nTrials * nTime), nTime, nTrials)'; % trials x time
figure('Name', 'Sample Index Order (Trials x Time)');
imagesc(sampleIndexGrid);
axis image
colormap(parula)
colorbar
title('Sample index order (rows = trials, cols = time)')
xlabel('Time index')
ylabel('Trial index')
set(gca, 'XTick', 1:nTime, 'YTick', 1:nTrials);

%% Compute cosine similarity and RDM for dot 1 directions
% Data flow: samples -> cosine similarity matrix -> dissimilarity (RDM).
[dot1CosSim, dot1Rdm] = cosine_similarity_rdm(dot1Samples);

figure('Name', 'Dot 1 Cosine Similarity');
imagesc(dot1CosSim);
axis image
colormap(parula)
colorbar
title('Dot 1 cosine similarity (samples x samples)')
set(gca, 'XTick', 1:(nTrials * nTime), 'YTick', 1:(nTrials * nTime));
set(gca, 'XTickLabel', sampleLabels, 'YTickLabel', sampleLabels, ...
    'TickLabelInterpreter', 'none');
xtickangle(45);

figure('Name', 'Dot 1 RDM (1 - Cosine Similarity)');
imagesc(dot1Rdm);
axis image
colormap(parula)
colorbar
title('Dot 1 RDM = 1 - cosine similarity')
set(gca, 'XTick', 1:(nTrials * nTime), 'YTick', 1:(nTrials * nTime));
set(gca, 'XTickLabel', sampleLabels, 'YTickLabel', sampleLabels, ...
    'TickLabelInterpreter', 'none');
xtickangle(45);

%% Compute cosine similarity and RDM for dot 2 directions
% Data flow: samples -> cosine similarity matrix -> dissimilarity (RDM).
[dot2CosSim, dot2Rdm] = cosine_similarity_rdm(dot2Samples);

figure('Name', 'Dot 2 Cosine Similarity');
imagesc(dot2CosSim);
axis image
colormap(parula)
colorbar
title('Dot 2 cosine similarity (samples x samples)')
set(gca, 'XTick', 1:(nTrials * nTime), 'YTick', 1:(nTrials * nTime));
set(gca, 'XTickLabel', sampleLabels, 'YTickLabel', sampleLabels, ...
    'TickLabelInterpreter', 'none');
xtickangle(45);

figure('Name', 'Dot 2 RDM (1 - Cosine Similarity)');
imagesc(dot2Rdm);
axis image
colormap(parula)
colorbar
title('Dot 2 RDM = 1 - cosine similarity')
set(gca, 'XTick', 1:(nTrials * nTime), 'YTick', 1:(nTrials * nTime));
set(gca, 'XTickLabel', sampleLabels, 'YTickLabel', sampleLabels, ...
    'TickLabelInterpreter', 'none');
xtickangle(45);

%% Helper functions (local to this script)
% Data flow: inputs -> visualization or similarity outputs for clarity.

function plot_trial_path_with_labels(trialXY, plotTitle)
% plot_trial_path_with_labels Plot a single trial path with labeled timepoints.
% Inputs:
%   trialXY: 2 x time matrix of [x; y] positions.
%   plotTitle: title string for the axes.
plot(trialXY(1, :), trialXY(2, :), '-o', 'LineWidth', 1.5);
axis equal
grid on
title(plotTitle)
xlabel('x')
ylabel('y')
for iTime = 1:size(trialXY, 2)
    text(trialXY(1, iTime), trialXY(2, iTime), ...
        sprintf('t%d', iTime), 'VerticalAlignment', 'bottom');
end
end

function plot_displacement_series(dx, dy, labelPrefix, rowIdx)
% plot_displacement_series Plot dx and dy time series for each trial.
% Inputs:
%   dx, dy: trials x 1 x time displacement matrices.
%   labelPrefix: string to identify the dot in titles.
%   rowIdx: row index for the 2x2 subplot grid (1 or 2).
nTrialsLocal = size(dx, 1);
nTimeLocal = size(dx, 3);
t = 1:nTimeLocal;

subplot(2, 2, (rowIdx - 1) * 2 + 1);
hold on
for iTrial = 1:nTrialsLocal
    plot(t, squeeze(dx(iTrial, 1, :)), '-o', 'LineWidth', 1.2);
end
hold off
grid on
title([labelPrefix ' dx over time'])
xlabel('time')
ylabel('dx')

subplot(2, 2, (rowIdx - 1) * 2 + 2);
hold on
for iTrial = 1:nTrialsLocal
    plot(t, squeeze(dy(iTrial, 1, :)), '-o', 'LineWidth', 1.2);
end
hold off
grid on
title([labelPrefix ' dy over time'])
xlabel('time')
ylabel('dy')
end

function plot_angle_series(angleDeg, labelPrefix, rowIdx)
% plot_angle_series Plot direction angle in degrees for each trial.
% Inputs:
%   angleDeg: trials x 1 x time matrix of angles in degrees.
%   labelPrefix: string to identify the dot in titles.
%   rowIdx: row index for the 2x1 subplot grid (1 or 2).
nTrialsLocal = size(angleDeg, 1);
nTimeLocal = size(angleDeg, 3);
t = 1:nTimeLocal;

subplot(2, 1, rowIdx);
hold on
for iTrial = 1:nTrialsLocal
    plot(t, squeeze(angleDeg(iTrial, 1, :)), '-o', 'LineWidth', 1.2);
end
hold off
grid on
title([labelPrefix ' angle (deg) over time'])
xlabel('time')
ylabel('angle (deg)')
end

function plot_path_with_vectors(trialXY, dx, dy, plotTitle)
% plot_path_with_vectors Overlay displacement vectors on a trial path.
% Inputs:
%   trialXY: 2 x time matrix of [x; y] positions.
%   dx, dy: 1 x time displacement vectors.
%   plotTitle: title string for the axes.
% Note:
%   Force dx/dy to row vectors so quiver sees consistent sizes.
dx = dx(:).';
dy = dy(:).';
plot(trialXY(1, :), trialXY(2, :), '-o', 'LineWidth', 1.5);
hold on
quiver(trialXY(1, :), trialXY(2, :), dx, dy, 0, 'LineWidth', 1.2);
hold off
axis equal
grid on
title(plotTitle)
xlabel('x')
ylabel('y')
end

function [cosSim, rdm] = cosine_similarity_rdm(samples)
% cosine_similarity_rdm Compute cosine similarity and RDM (1 - similarity).
% Inputs:
%   samples: nSamples x nFeatures matrix (rows are samples).
% Outputs:
%   cosSim: nSamples x nSamples cosine similarity matrix.
%   rdm: nSamples x nSamples dissimilarity matrix (1 - cosine similarity).
norms = sqrt(sum(samples.^2, 2));
norms = max(norms, eps); % avoid divide-by-zero if any sample is zero
cosSim = (samples * samples') ./ (norms * norms');
rdm = 1 - cosSim;
end
