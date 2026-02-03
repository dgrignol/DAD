% compute_avg_path_distance.m
%
% Purpose:
%   Compute the average Euclidean distance between two dot paths across
%   trials and timepoints, with NaN-safe summaries for missing samples.
%
% Example usage:
%   [meanDistance, meanByTrial, meanByTime] = ...
%       compute_avg_path_distance(dot1Paths, dot2Paths);
%
% Inputs:
%   dot1Paths: trials x 2 x time array of [x;y] positions (visual degrees).
%   dot2Paths: trials x 2 x time array of [x;y] positions (same shape).
%
% Outputs:
%   meanDistance: scalar average distance across all trials and timepoints.
%   meanDistanceByTrial: trials x 1 mean distance (averaged over timepoints).
%   meanDistanceByTime: 1 x time mean distance (averaged across trials).
%   distanceMatrix: trials x 1 x time Euclidean distances per sample.
%
% Assumptions:
%   - Input paths are aligned in trial/time and share the same shape.
%   - Missing samples are encoded as NaN and excluded from means.
function [meanDistance, meanDistanceByTrial, meanDistanceByTime, distanceMatrix] = ...
    compute_avg_path_distance(dot1Paths, dot2Paths)

%% Validate inputs and shapes
% Data flow: raw arrays -> size checks -> error for mismatched path geometry.
if ~isequal(size(dot1Paths), size(dot2Paths))
    error('Dot path arrays must have identical sizes.');
end
if size(dot1Paths, 2) ~= 2
    error('Dot path arrays must be trials x 2 x time.');
end

%% Compute per-sample Euclidean distances
% Data flow: x/y differences -> squared distance -> sqrt to recover Euclidean norm.
dx = dot1Paths(:, 1, :) - dot2Paths(:, 1, :);
dy = dot1Paths(:, 2, :) - dot2Paths(:, 2, :);
distanceMatrix = sqrt(dx.^2 + dy.^2);

%% Summarize distances across trials and time
% Data flow: distanceMatrix -> mean over time, mean over trials, global mean.
meanDistanceByTrial = squeeze(mean(distanceMatrix, 3, 'omitnan'));
meanDistanceByTime = squeeze(mean(distanceMatrix, 1, 'omitnan'));
meanDistance = mean(distanceMatrix, 'all', 'omitnan');
end
