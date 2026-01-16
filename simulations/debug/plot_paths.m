function figHandle = plot_paths(paths, varargin)
% plot_paths
%
% Purpose:
%   Plot 2D dot paths (trials × features × time) using a time-graded color
%   ramp, so early samples appear cool and late samples appear warm.
%
% Example usage:
%   % dotPaths: trials × 2 × time (visual degrees)
%   fig = plot_paths(dot1GreenPaths, 'Title', 'Dot 1 paths');
%
% Example usage (with transparency):
%   fig = plot_paths(dot1GreenPaths, 'UseAlpha', true, 'AlphaValue', 0.05);
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
%
% Key assumptions:
%   - paths is trials × 2 × time, where features are [x, y].
%   - Time dimension is shared across trials.

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
    figHandle = figure('Name', 'Dot paths', 'NumberTitle', 'off');
    ax = axes(figHandle);
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

    %% Save plot to simulations/debug
    debugDir = fullfile(fileparts(mfilename('fullpath')), 'debug');
    if ~exist(debugDir, 'dir')
        mkdir(debugDir);
    end
    safeTitle = regexprep(char(opts.Title), '[^A-Za-z0-9_-]+', '_');
    if isempty(strtrim(safeTitle))
        safeTitle = 'plot_paths';
    end
    outFile = fullfile(debugDir, [safeTitle '.png']);
    print(figHandle, outFile, '-dpng', '-r300');
end
