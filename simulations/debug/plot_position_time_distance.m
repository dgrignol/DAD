function figHandle = plot_position_time_distance(paths, varargin)
% plot_position_time_distance
%
% Purpose:
%   Visualize a time-by-time distance matrix computed from dot positions.
%   Each cell (t, t') is the mean Euclidean distance between positions at
%   time t and time t' across trials. Central "squares" indicate time
%   windows where positions are more tightly clustered across trials.
%
% Example usage:
%   % paths: trials × 2 × time
%   fig = plot_position_time_distance(dot1GreenPaths, 'Title', 'Dot 1 distance');
%
% Example usage (custom colormap):
%   fig = plot_position_time_distance(dot1GreenPaths, 'Colormap', 'hot');
%
% Inputs:
%   paths : numeric array, trials × 2 × time (features are [x, y])
%
% Name/value options:
%   'Title'    : figure title string (default: '')
%   'Colormap' : colormap name or n×3 array (default: parula)
%   'AxisEqual': true/false to enforce equal axes (default: true)
%
% Outputs:
%   figHandle : handle to the created figure
%
% Key assumptions:
%   - paths is trials × 2 × time, where features are [x, y].
%   - Distances are computed in the same units as the input (visual degrees).

    %% Input validation and options
    parser = inputParser;
    parser.addRequired('paths', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addParameter('Title', '', @(x) ischar(x) || isstring(x));
    parser.addParameter('Colormap', parula, @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && size(x, 2) == 3));
    parser.addParameter('AxisEqual', true, @(x) islogical(x) && isscalar(x));
    parser.parse(paths, varargin{:});

    opts = parser.Results;

    %% Compute mean distance matrix across trials
    % Data flow: trial positions -> per-trial distance matrices -> mean.
    nTrials = size(paths, 1);
    nTime = size(paths, 3);
    distSum = zeros(nTime, nTime);

    for iTrial = 1:nTrials
        pos = squeeze(paths(iTrial, :, :))'; % time × 2
        distSum = distSum + pdist2(pos, pos);
    end

    distMean = distSum ./ nTrials;

    %% Plot distance matrix
    figHandle = figure('Name', 'Time-by-time distance', 'NumberTitle', 'off');
    ax = axes(figHandle);
    imagesc(ax, distMean);
    set(ax, 'YDir', 'normal');
    colormap(ax, opts.Colormap);
    colorbar(ax);
    title(ax, char(opts.Title));
    xlabel(ax, 'Time (samples)');
    ylabel(ax, 'Time (samples)');

    if opts.AxisEqual
        axis(ax, 'image');
    end

    %% Save plot to simulations/debug
    debugDir = fullfile(fileparts(mfilename('fullpath')), 'debug');
    if ~exist(debugDir, 'dir')
        mkdir(debugDir);
    end
    safeTitle = regexprep(char(opts.Title), '[^A-Za-z0-9_-]+', '_');
    if isempty(strtrim(safeTitle))
        safeTitle = 'plot_position_time_distance';
    end
    outFile = fullfile(debugDir, [safeTitle '.png']);
    print(figHandle, outFile, '-dpng', '-r300');
end
