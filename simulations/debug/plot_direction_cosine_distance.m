function figHandle = plot_direction_cosine_distance(dotDirection, varargin)
% plot_direction_cosine_distance
%
% Purpose:
%   Plot a time-by-time cosine distance matrix for a single trial from a
%   dot-direction array (trials x features x time). This mirrors the
%   cosine distance used by dRSA when params.modelDistMeasure = 'cosine'.
%
% Example usage:
%   % dot2Direction: trials x 2 x time
%   fig = plot_direction_cosine_distance(dot2Direction);
%
% Example usage (explicit trial selection):
%   fig = plot_direction_cosine_distance(dot2Direction, ...
%       'TrialIndex', 4, ...
%       'Title', 'Dot 2 cosine distance (trial 4)');
%
% Inputs:
%   dotDirection : numeric array, trials x 2 x time (features are [cos, sin])
%
% Name/value options:
%   'TrialIndex' : scalar trial index (default: [], uses random trial)
%   'Title'      : figure title string (default: 'Cosine distance (random trial)')
%
% Outputs:
%   figHandle : handle to the created figure
%
% Key assumptions:
%   - Feature 1 is cosine and feature 2 is sine of the direction angle.
%   - Cosine distance matches MATLAB pdist(..., 'cosine') = 1 - cosine similarity.

    %% Validate inputs and parse options
    % Data flow: raw dotDirection -> validated shape -> plotting options.
    parser = inputParser;
    parser.addRequired('dotDirection', ...
        @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) >= 2);
    parser.addParameter('TrialIndex', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x, 1) == 0));
    parser.addParameter('Title', 'Cosine distance (random trial)', ...
        @(x) ischar(x) || isstring(x));
    parser.parse(dotDirection, varargin{:});
    opts = parser.Results;

    %% Select the trial to plot
    % Data flow: trial selection -> time-by-feature samples.
    trialCount = size(dotDirection, 1);
    if isempty(opts.TrialIndex)
        trialIdx = randi(trialCount);
    else
        trialIdx = opts.TrialIndex;
    end
    if trialIdx > trialCount
        error('TrialIndex (%d) exceeds trial count (%d).', trialIdx, trialCount);
    end

    trialSeries = squeeze(dotDirection(trialIdx, 1:2, :)); % 2 x time
    samples = trialSeries'; % time x features for pdist

    %% Compute cosine distance matrix
    % Data flow: time samples -> pdist cosine -> squareform matrix.
    distVec = pdist(samples, 'cosine');
    distMat = squareform(distVec);

    %% Plot the cosine distance matrix
    % Data flow: distance matrix -> heatmap with time ticks.
    figHandle = figure('Name', char(opts.Title), 'NumberTitle', 'off');
    imagesc(distMat);
    axis square;
    axis xy;
    colorbar;
    title(char(opts.Title));
    xlabel('Time (samples)');
    ylabel('Time (samples)');

    nTime = size(distMat, 1);
    tickCandidates = [80 160 240 320];
    tickCandidates = tickCandidates(tickCandidates <= nTime);
    if isempty(tickCandidates)
        tickCandidates = unique(round(linspace(1, nTime, min(4, nTime))));
    end
    xticks(tickCandidates);
    yticks(tickCandidates);
end
