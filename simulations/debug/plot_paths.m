function figHandle = plot_paths(paths, varargin)
% plot_paths
%
% Purpose:
%   Plot 2D dot paths (trials × features × time) using a time-graded color
%   ramp, so early samples appear cool and late samples appear warm. The
%   colorbar is expressed in seconds using the provided sampling rate.
%
% Example usage (basic):
%   % dotPaths: trials × 2 × time (visual degrees)
%   fig = plot_paths(dot1GreenPaths, 'Title', 'Dot 1 paths', 'SampleRateHz', 120);
%
% Example usage (with transparency):
%   fig = plot_paths(dot1GreenPaths, 'UseAlpha', true, 'AlphaValue', 0.05, ...
%       'SampleRateHz', 120);
%
% Example usage (subplot integration):
%   fig = figure;
%   ax = subplot(2, 2, 1);
%   plot_paths(dot1GreenPaths, 'ParentAxes', ax, 'SavePlot', false, ...
%       'Title', 'Dot 1', 'SampleRateHz', 120);
%
% Inputs:
%   paths : numeric array, trials × 2 × time (features are [x, y])
%
% Name/value options:
%   'Title'        : figure title string (default: '')
%   'Colormap'     : colormap name or n×3 array (default: parula)
%   'MarkerSize'   : scatter marker size (default: 8)
%   'UseAlpha'     : true/false to enable transparency (default: false)
%   'AlphaValue'   : alpha for marker face/edge when UseAlpha = true (default: 0.1)
%   'AxisEqual'    : true/false to enforce equal axes (default: true)
%   'AxisLimits'   : 1×4 [xmin xmax ymin ymax] (default: [])
%   'SampleRateHz' : sampling rate (Hz) for seconds colorbar (default: 1)
%   'ParentAxes'   : axes handle to plot into (default: [])
%   'SavePlot'     : true/false to save PNG into simulations/debug (default: true)
%
% Outputs:
%   figHandle : handle to the created figure (or parent figure if provided)
%
% Key assumptions:
%   - paths is trials × 2 × time, where features are [x, y].
%   - Time dimension is shared across trials.
%   - SampleRateHz reflects the sampling rate used to interpret seconds.

    %% Input validation and options
    parser = inputParser;
    parser.addRequired('paths', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addParameter('Title', '', @(x) ischar(x) || isstring(x));
    parser.addParameter('Colormap', parula, @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && size(x, 2) == 3));
    parser.addParameter('MarkerSize', 8, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('UseAlpha', false, @(x) islogical(x) && isscalar(x));
    parser.addParameter('AlphaValue', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    parser.addParameter('AxisEqual', true, @(x) islogical(x) && isscalar(x));
    parser.addParameter('AxisLimits', [-5 5 -5 5], @(x) isempty(x) || (isnumeric(x) && numel(x) == 4));
    parser.addParameter('SampleRateHz', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('ParentAxes', [], @(x) isempty(x) || isgraphics(x, 'axes'));
    parser.addParameter('SavePlot', true, @(x) islogical(x) && isscalar(x));
    parser.parse(paths, varargin{:});

    opts = parser.Results;

    %% Flatten trials into a single scatter set
    % Data flow: trials × 2 × time -> x/y vectors + time indices.
    nTrials = size(paths, 1);
    nTime = size(paths, 3);
    timeSeconds = (0:(nTime - 1)) / opts.SampleRateHz;

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
    % Data flow: x/y/time vectors -> scatter -> colorbar in seconds.
    if isempty(opts.ParentAxes)
        figHandle = figure('Name', 'Dot paths', 'NumberTitle', 'off');
        ax = axes(figHandle);
    else
        ax = opts.ParentAxes;
        figHandle = ancestor(ax, 'figure');
    end
    hold(ax, 'on');

    if opts.UseAlpha
        scatter(ax, xVals, yVals, opts.MarkerSize, pointColors, 'filled', ...
            'MarkerFaceAlpha', opts.AlphaValue, 'MarkerEdgeAlpha', opts.AlphaValue);
    else
        scatter(ax, xVals, yVals, opts.MarkerSize, pointColors, 'filled');
    end

    colormap(ax, cmap);
    caxis(ax, [timeSeconds(1) timeSeconds(end)]);
    cb = colorbar(ax);
    cb.Label.String = 'Time (s)';
    cb.Ticks = linspace(timeSeconds(1), timeSeconds(end), 5);
    cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);

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
    % Data flow: figure -> PNG on disk (optional).
    if opts.SavePlot
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
end
