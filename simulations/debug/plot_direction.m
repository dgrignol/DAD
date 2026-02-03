function figHandle = plot_direction(dotDirection, varargin)
% plot_direction
%
% Purpose:
%   Plot a single trial from a dot-direction array (trials x features x time)
%   using two stacked subplots: cosine (feature 1) and sine (feature 2).
%   Intended as a quick diagnostic for direction timecourses.
%
% Example usage:
%   % dot2Direction: trials x 2 x time
%   fig = plot_direction(dot2Direction);
%
% Example usage (explicit trial selection):
%   fig = plot_direction(dot2Direction, 'TrialIndex', 3, 'Title', 'Dot 2 direction');
%
% Inputs:
%   dotDirection : numeric array, trials x 2 x time (features are [cos, sin])
%
% Name/value options:
%   'TrialIndex' : scalar trial index (default: [], uses random trial)
%   'Title'      : figure title string (default: 'Dot direction (random trial)')
%   'LineWidth'  : line width (default: 1.5)
%
% Outputs:
%   figHandle : handle to the created figure
%
% Key assumptions:
%   - Feature 1 is cosine and feature 2 is sine of the direction angle.
%   - Timepoints are sampled uniformly, plotted as 1..N samples.

    %% Validate inputs and parse options
    % Data flow: raw dotDirection -> validated shape -> plotting options.
    parser = inputParser;
    parser.addRequired('dotDirection', ...
        @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) >= 2);
    parser.addParameter('TrialIndex', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x, 1) == 0));
    parser.addParameter('Title', 'Dot direction (random trial)', ...
        @(x) ischar(x) || isstring(x));
    parser.addParameter('LineWidth', 1.5, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    parser.parse(dotDirection, varargin{:});
    opts = parser.Results;

    %% Select the trial to plot
    % Data flow: trial selection -> cosine/sine time series.
    trialCount = size(dotDirection, 1);
    if isempty(opts.TrialIndex)
        trialIdx = randi(trialCount);
    else
        trialIdx = opts.TrialIndex;
    end
    if trialIdx > trialCount
        error('TrialIndex (%d) exceeds trial count (%d).', trialIdx, trialCount);
    end

    nTime = size(dotDirection, 3);
    timeVec = 1:nTime;
    cosSeries = squeeze(dotDirection(trialIdx, 1, :));
    sinSeries = squeeze(dotDirection(trialIdx, 2, :));

    %% Plot cosine and sine features over time
    % Data flow: time vector + feature series -> stacked diagnostic plot.
    figHandle = figure('Name', char(opts.Title), 'NumberTitle', 'off');
    tickCandidates = [80 160 240 320];
    tickCandidates = tickCandidates(tickCandidates <= nTime);
    if isempty(tickCandidates)
        tickCandidates = unique(round(linspace(1, nTime, min(4, nTime))));
    end

    subplot(2, 1, 1);
    plot(timeVec, cosSeries, 'LineWidth', opts.LineWidth);
    xlim([1 nTime]);
    xticks(tickCandidates);
    xlabel('Time (samples)');
    ylabel('Cosine');
    title('Feature 1 (cos)');
    grid on;

    subplot(2, 1, 2);
    plot(timeVec, sinSeries, 'LineWidth', opts.LineWidth);
    xlim([1 nTime]);
    xticks(tickCandidates);
    xlabel('Time (samples)');
    ylabel('Sine');
    title('Feature 2 (sin)');
    grid on;
end
