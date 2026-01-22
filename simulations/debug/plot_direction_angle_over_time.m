function figHandle = plot_direction_angle_over_time(paths, varargin)
% plot_direction_angle_over_time
%
% Purpose:
%   Plot per-trial direction angle over time from dot position paths.
%   Direction is computed from frame-to-frame displacement; the first
%   timepoint uses the t1->t2 displacement to infer direction.
%
% Example usage:
%   % paths: trials x 2 x time (visual degrees)
%   fig = plot_direction_angle_over_time(dot1GreenPaths, ...
%       'DotLabel', 'dot1');
%
% Example usage (plot fewer trials):
%   fig = plot_direction_angle_over_time(dot2YellowPaths, ...
%       'DotLabel', 'dot2', 'NumTrialsToPlot', 3);
%
% Inputs:
%   paths : numeric array, trials x 2 x time (features are [x, y])
%
% Name/value options:
%   'DotLabel'        : label used in the title (default: 'dot')
%   'NumTrialsToPlot' : number of trials to show (default: 5)
%   'Unwrap'          : true to unwrap angles over time (default: true)
%
% Outputs:
%   figHandle : handle to the created figure
%
% Key assumptions:
%   - paths are in visual degrees, sampled uniformly in time.

    %% Input validation and options
    % Data flow: raw inputs -> validated options for plotting.
    parser = inputParser;
    parser.addRequired('paths', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addParameter('DotLabel', 'dot', @(x) ischar(x) || isstring(x));
    parser.addParameter('NumTrialsToPlot', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('Unwrap', true, @(x) islogical(x) && isscalar(x));
    parser.parse(paths, varargin{:});
    opts = parser.Results;

    %% Compute direction angles (radians)
    % Data flow: positions -> displacement -> angle per trial/time.
    dx = diff(paths(:, 1, :), 1, 3);
    dy = diff(paths(:, 2, :), 1, 3);
    dx = cat(3, dx(:, :, 1), dx);
    dy = cat(3, dy(:, :, 1), dy);
    angles = atan2(dy, dx);
    angles = squeeze(angles); % trials x time
    if opts.Unwrap
        angles = unwrap(angles, [], 2);
    end
    anglesDeg = rad2deg(angles);

    %% Plot mean and example trials
    % Data flow: per-trial angles -> mean + selected trials -> figure.
    nTrials = size(anglesDeg, 1);
    nPlot = min(opts.NumTrialsToPlot, nTrials);
    figHandle = figure('Name', 'Direction angle over time', 'NumberTitle', 'off');
    ax = axes(figHandle);
    hold(ax, 'on');
    plot(ax, mean(anglesDeg, 1, 'omitnan'), 'k-', 'LineWidth', 2);
    for iTrial = 1:nPlot
        plot(ax, anglesDeg(iTrial, :), 'LineWidth', 1);
    end
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'Time (samples)');
    ylabel(ax, 'Direction (degrees)');
    title(ax, sprintf('Direction angle over time (%s)', char(opts.DotLabel)));
    legend(ax, ['mean', arrayfun(@(x) sprintf('trial %d', x), ...
        1:nPlot, 'UniformOutput', false)], 'Location', 'best');
end
