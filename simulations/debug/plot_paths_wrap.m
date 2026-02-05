function figHandle = plot_paths_wrap(dot1Nondeviant, dot2Nondeviant, dot1Deviant, dot2Deviant, varargin)
% plot_paths_wrap
%
% Purpose:
%   Plot dot-motion paths for all conditions on a single figure with
%   four subplots (2x2). Columns correspond to dot1/dot2, and rows
%   correspond to nondeviant/deviant conditions. Each subplot uses a
%   time-graded color ramp with a seconds-based colorbar.
%
% Example usage (arrays already loaded):
%   fig = plot_paths_wrap(dot1Nondev, dot2Nondev, dot1Dev, dot2Dev, ...
%       'SampleCount', 20, 'SampleRateHz', 120, 'Title', 'All conditions');
%
% Example usage (with transparency):
%   fig = plot_paths_wrap(dot1Nondev, dot2Nondev, dot1Dev, dot2Dev, ...
%       'UseAlpha', true, 'AlphaValue', 0.05, 'SampleRateHz', 120);
%
% Inputs:
%   dot1Nondeviant : numeric array, trials × 2 × time (dot1, nondeviant)
%   dot2Nondeviant : numeric array, trials × 2 × time (dot2, nondeviant)
%   dot1Deviant    : numeric array, trials × 2 × time (dot1, deviant)
%   dot2Deviant    : numeric array, trials × 2 × time (dot2, deviant)
%
% Name/value options:
%   'Title'        : figure title string (default: 'Dot paths (all conditions)')
%   'Colormap'     : colormap name or n×3 array (default: parula)
%   'MarkerSize'   : scatter marker size (default: 8)
%   'UseAlpha'     : true/false to enable transparency (default: false)
%   'AlphaValue'   : alpha for marker face/edge when UseAlpha = true (default: 0.1)
%   'AxisEqual'    : true/false to enforce equal axes (default: true)
%   'AxisLimits'   : 1×4 [xmin xmax ymin ymax] (default: [])
%   'SampleRateHz' : sampling rate (Hz) for seconds colorbar (default: 1)
%   'SampleCount'  : number of paths to sample per condition (default: 20)
%
% Outputs:
%   figHandle : handle to the created figure
%
% Key assumptions:
%   - Inputs are trials × 2 × time in visual degrees, with shared time length.
%   - Each condition is sampled independently (with replacement).

    %% Input validation and options
    parser = inputParser;
    parser.addRequired('dot1Nondeviant', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addRequired('dot2Nondeviant', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addRequired('dot1Deviant', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addRequired('dot2Deviant', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addParameter('Title', 'Dot paths (all conditions)', @(x) ischar(x) || isstring(x));
    parser.addParameter('Colormap', parula, @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && size(x, 2) == 3));
    parser.addParameter('MarkerSize', 8, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('UseAlpha', false, @(x) islogical(x) && isscalar(x));
    parser.addParameter('AlphaValue', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    parser.addParameter('AxisEqual', true, @(x) islogical(x) && isscalar(x));
    parser.addParameter('AxisLimits', [-5 5 -5 5], @(x) isempty(x) || (isnumeric(x) && numel(x) == 4));
    parser.addParameter('SampleRateHz', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.addParameter('SampleCount', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.parse(dot1Nondeviant, dot2Nondeviant, dot1Deviant, dot2Deviant, varargin{:});

    opts = parser.Results;

    %% Sample paths independently per condition
    % Data flow: full condition arrays -> random trial indices -> sampled paths.
    dot1NondevSample = sample_paths(dot1Nondeviant, opts.SampleCount);
    dot2NondevSample = sample_paths(dot2Nondeviant, opts.SampleCount);
    dot1DevSample = sample_paths(dot1Deviant, opts.SampleCount);
    dot2DevSample = sample_paths(dot2Deviant, opts.SampleCount);

    %% Create the 2x2 subplot layout
    % Data flow: sampled paths -> plot_paths (axes) -> combined figure.
    figHandle = figure('Name', char(opts.Title), 'NumberTitle', 'off');
    sgtitle(figHandle, char(opts.Title));

    ax11 = subplot(2, 2, 1, 'Parent', figHandle);
    plot_paths(dot1NondevSample, 'ParentAxes', ax11, 'SavePlot', false, ...
        'Title', 'Nondeviant: dot1', 'Colormap', opts.Colormap, ...
        'MarkerSize', opts.MarkerSize, 'UseAlpha', opts.UseAlpha, ...
        'AlphaValue', opts.AlphaValue, 'AxisEqual', opts.AxisEqual, ...
        'AxisLimits', opts.AxisLimits, 'SampleRateHz', opts.SampleRateHz);

    ax12 = subplot(2, 2, 2, 'Parent', figHandle);
    plot_paths(dot2NondevSample, 'ParentAxes', ax12, 'SavePlot', false, ...
        'Title', 'Nondeviant: dot2', 'Colormap', opts.Colormap, ...
        'MarkerSize', opts.MarkerSize, 'UseAlpha', opts.UseAlpha, ...
        'AlphaValue', opts.AlphaValue, 'AxisEqual', opts.AxisEqual, ...
        'AxisLimits', opts.AxisLimits, 'SampleRateHz', opts.SampleRateHz);

    ax21 = subplot(2, 2, 3, 'Parent', figHandle);
    plot_paths(dot1DevSample, 'ParentAxes', ax21, 'SavePlot', false, ...
        'Title', 'Deviant: dot1', 'Colormap', opts.Colormap, ...
        'MarkerSize', opts.MarkerSize, 'UseAlpha', opts.UseAlpha, ...
        'AlphaValue', opts.AlphaValue, 'AxisEqual', opts.AxisEqual, ...
        'AxisLimits', opts.AxisLimits, 'SampleRateHz', opts.SampleRateHz);

    ax22 = subplot(2, 2, 4, 'Parent', figHandle);
    plot_paths(dot2DevSample, 'ParentAxes', ax22, 'SavePlot', false, ...
        'Title', 'Deviant: dot2', 'Colormap', opts.Colormap, ...
        'MarkerSize', opts.MarkerSize, 'UseAlpha', opts.UseAlpha, ...
        'AlphaValue', opts.AlphaValue, 'AxisEqual', opts.AxisEqual, ...
        'AxisLimits', opts.AxisLimits, 'SampleRateHz', opts.SampleRateHz);

    %% Save combined plot to simulations/debug
    % Data flow: combined figure -> PNG on disk.
    debugDir = fullfile(fileparts(mfilename('fullpath')), 'debug');
    if ~exist(debugDir, 'dir')
        mkdir(debugDir);
    end
    safeTitle = regexprep(char(opts.Title), '[^A-Za-z0-9_-]+', '_');
    if isempty(strtrim(safeTitle))
        safeTitle = 'plot_paths_wrap';
    end
    outFile = fullfile(debugDir, [safeTitle '.png']);
    print(figHandle, outFile, '-dpng', '-r300');
end

%% Local helper: sample trials with replacement
function sampledPaths = sample_paths(paths, sampleCount)
% sample_paths
%
% Purpose:
%   Draw a random subset of trials (with replacement) for plotting.
%
% Inputs:
%   paths       : trials × 2 × time array
%   sampleCount : number of trials to sample
%
% Outputs:
%   sampledPaths : sampleCount × 2 × time array

    nTrials = size(paths, 1);
    nPick = min(sampleCount, nTrials);
    drawIdx = randsample(nTrials, nPick, true);
    sampledPaths = paths(drawIdx, :, :);
end
