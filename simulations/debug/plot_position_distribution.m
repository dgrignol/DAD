function figHandle = plot_position_distribution(paths, varargin)
% plot_position_distribution
%
% Purpose:
%   Plot the mean distance from the stimulus center at each timepoint for a
%   set of dot paths (trials x 2 x time). The center is defined by the
%   stimulus borders from stimuli_generation_v5.m (Cfg.rectSize) when
%   provided, or inferred from data bounds as a fallback.
%
% Example usage:
%   % dotPaths: trials x 2 x time (visual degrees)
%   fig = plot_position_distribution(dot1GreenPaths, ...
%       'Title', 'Dot 1 distance from center', ...
%       'RectSize', [10 10]);
%
% Example usage (with inferred bounds):
%   fig = plot_position_distribution(dot2YellowPaths, ...
%       'Title', 'Dot 2 distance from center');
%
% Inputs:
%   paths : numeric array, trials x 2 x time (features are [x, y])
%
% Name/value options:
%   'Title'      : figure title string (default: '')
%   'RectSize'   : 1x2 [width height] in visual degrees (default: [])
%   'RectOrigin' : 1x2 [x0 y0] origin of the rect (default: [0 0])
%   'LineColor'  : 1x3 RGB for the line (default: [0 0 0])
%   'LineWidth'  : line width (default: 1.5)
%
% Outputs:
%   figHandle : handle to the created figure
%
% Key assumptions:
%   - Coordinates are in visual degrees, with borders defined by
%     Config.dotRectSize saved as Cfg.rectSize.
%   - If RectSize is empty, bounds are inferred from data and may under-
%     estimate the true borders if paths do not reach edges.

    %% Input validation and options
    parser = inputParser;
    parser.addRequired('paths', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addParameter('Title', '', @(x) ischar(x) || isstring(x));
    parser.addParameter('RectSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    parser.addParameter('RectOrigin', [0 0], @(x) isnumeric(x) && numel(x) == 2);
    parser.addParameter('LineColor', [0 0 0], @(x) isnumeric(x) && numel(x) == 3);
    parser.addParameter('LineWidth', 1.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.parse(paths, varargin{:});

    opts = parser.Results;

    %% Resolve rect bounds and center
    % Data flow: rect size/origin (or inferred bounds) -> center coordinate.
    if isempty(opts.RectSize)
        % Fallback when Cfg.rectSize is unavailable: infer bounds from data.
        xVals = reshape(paths(:, 1, :), [], 1);
        yVals = reshape(paths(:, 2, :), [], 1);
        xVals = xVals(~isnan(xVals));
        yVals = yVals(~isnan(yVals));
        if isempty(xVals) || isempty(yVals)
            error('Cannot infer bounds: all positions are NaN.');
        end
        rectOrigin = [min(xVals) min(yVals)];
        rectSize = [max(xVals) - min(xVals), max(yVals) - min(yVals)];
    else
        rectSize = reshape(opts.RectSize, 1, 2);
        rectOrigin = reshape(opts.RectOrigin, 1, 2);
    end
    center = rectOrigin + rectSize / 2;

    %% Compute mean distance per timepoint
    % Data flow: trial positions -> distance to center -> mean across trials.
    nTime = size(paths, 3);
    meanDist = zeros(1, nTime);

    for iTime = 1:nTime
        positions = squeeze(paths(:, :, iTime)); % trials x 2
        deltas = positions - center;
        distances = sqrt(sum(deltas.^2, 2));
        validMask = ~any(isnan(positions), 2);
        if any(validMask)
            meanDist(iTime) = mean(distances(validMask));
        else
            meanDist(iTime) = NaN;
        end
    end

    %% Plot distance profile
    figHandle = figure('Name', 'Mean distance from center', 'NumberTitle', 'off');
    ax = axes(figHandle);
    plot(ax, meanDist, 'Color', opts.LineColor, 'LineWidth', opts.LineWidth);
    grid(ax, 'on');
    xlabel(ax, 'Time (samples)');
    ylabel(ax, 'Mean distance from center (visual degrees)');
    title(ax, char(opts.Title));

    %% Save plot to simulations/debug
    debugDir = fullfile(fileparts(mfilename('fullpath')), 'debug');
    if ~exist(debugDir, 'dir')
        mkdir(debugDir);
    end
    safeTitle = regexprep(char(opts.Title), '[^A-Za-z0-9_-]+', '_');
    if isempty(strtrim(safeTitle))
        safeTitle = 'plot_position_distribution';
    end
    outFile = fullfile(debugDir, [safeTitle '.png']);
    print(figHandle, outFile, '-dpng', '-r300');
end
