function [figHandle, diag] = plot_direction_dRSA_by_turn(dot1Trials, dot2Trials, params, varargin)
% plot_direction_dRSA_by_turn
%
% Purpose:
%   Split dot1 trials by angular-velocity sign (CCW vs CW) and recompute
%   the direction dRSA map for each group. This diagnostic tests whether
%   mixing opposite turn directions produces the periodic stripe pattern.
%
% Example usage:
%   % dot1Trials/dot2Trials: trials x 2 x time (visual degrees)
%   fig = plot_direction_dRSA_by_turn(dot1Trials, dot2Trials, params);
%
% Example usage (custom title and minimum trial count):
%   fig = plot_direction_dRSA_by_turn(dot1Trials, dot2Trials, params, ...
%       'TitlePrefix', 'Dot 1 direction dRSA', ...
%       'MinTrialsPerGroup', 4);
%
% Inputs:
%   dot1Trials : numeric array, trials x 2 x time (dot 1 positions)
%   dot2Trials : numeric array, trials x 2 x time (dot 2 positions)
%   params     : struct with dRSA parameters (uses neuralDistMeasure, dRSAtype)
%
% Name/value options:
%   'TitlePrefix'       : prefix for subplot titles (default: 'dRSA direction dot1')
%   'MinTrialsPerGroup' : minimum trials per group (default: 3)
%   'SlopeEps'          : threshold for near-zero angular velocity (default: 0)
%
% Outputs:
%   figHandle : handle to the created figure (empty if a group is too small)
%   diag      : struct with indices, slope estimates, and group counts
%
% Key assumptions:
%   - Direction angles are approximately linear in time per trial.
%   - Positive angular velocity corresponds to CCW in atan2 convention.
%   - dRSA maps are computed with cosine distance for direction models.

    %% Input validation and options
    % Data flow: raw inputs -> validated sizes -> diagnostic options.
    parser = inputParser;
    parser.addRequired('dot1Trials', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addRequired('dot2Trials', @(x) isnumeric(x) && ndims(x) == 3 && size(x, 2) == 2);
    parser.addRequired('params', @isstruct);
    parser.addParameter('TitlePrefix', 'dRSA direction dot1', @(x) ischar(x) || isstring(x));
    parser.addParameter('MinTrialsPerGroup', 3, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 2 && mod(x, 1) == 0);
    parser.addParameter('SlopeEps', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parser.parse(dot1Trials, dot2Trials, params, varargin{:});
    opts = parser.Results;

    if size(dot1Trials, 1) ~= size(dot2Trials, 1) || size(dot1Trials, 3) ~= size(dot2Trials, 3)
        error('dot1Trials and dot2Trials must have matching trial and time dimensions.');
    end
    if isfield(params, 'dRSAtype') && ~strcmp(params.dRSAtype, 'corr')
        error('plot_direction_dRSA_by_turn currently supports params.dRSAtype = ''corr'' only.');
    end

    %% Estimate per-trial angular velocity sign
    % Data flow: positions -> displacement -> unwrapped angles -> mean slope -> sign.
    dx = diff(dot1Trials(:, 1, :), 1, 3);
    dy = diff(dot1Trials(:, 2, :), 1, 3);
    dx = cat(3, dx(:, :, 1), dx);
    dy = cat(3, dy(:, :, 1), dy);
    angles = atan2(dy, dx); % trials x 1 x time
    angles = squeeze(angles); % trials x time
    angles = unwrap(angles, [], 2);
    meanStep = mean(diff(angles, 1, 2), 2, 'omitnan');
    meanStep(abs(meanStep) <= opts.SlopeEps) = 0;

    ccwMask = meanStep > 0;
    cwMask = meanStep < 0;
    zeroMask = meanStep == 0;

    diag = struct();
    diag.meanStep = meanStep;
    diag.ccwIdx = find(ccwMask);
    diag.cwIdx = find(cwMask);
    diag.zeroIdx = find(zeroMask);
    diag.counts = [numel(diag.ccwIdx), numel(diag.cwIdx)];

    if diag.counts(1) < opts.MinTrialsPerGroup || diag.counts(2) < opts.MinTrialsPerGroup
        warning('Not enough trials per group (CCW=%d, CW=%d).', diag.counts(1), diag.counts(2));
        figHandle = [];
        return;
    end

    %% Compute dRSA maps for CCW and CW trial groups
    % Data flow: trial subsets -> direction models -> subsamples -> dRSA maps.
    dRSA_ccw = compute_direction_dRSA(dot1Trials(ccwMask, :, :), dot2Trials(ccwMask, :, :), params);
    dRSA_cw = compute_direction_dRSA(dot1Trials(cwMask, :, :), dot2Trials(cwMask, :, :), params);

    %% Plot dRSA maps with shared color limits
    % Data flow: CCW/CW dRSA matrices -> side-by-side heatmaps.
    figHandle = figure('Name', 'dRSA direction dot1 by turn', 'NumberTitle', 'off');
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    cMin = min([min(dRSA_ccw(:)), min(dRSA_cw(:))]);
    cMax = max([max(dRSA_ccw(:)), max(dRSA_cw(:))]);
    if cMin == cMax
        cMin = cMin - 1;
        cMax = cMax + 1;
    end

    tickCandidates = [80 160 240 320];
    nTime = size(dRSA_ccw, 1);
    tickCandidates = tickCandidates(tickCandidates <= nTime);
    if isempty(tickCandidates)
        tickCandidates = unique(round(linspace(1, nTime, min(4, nTime))));
    end

    nexttile;
    imagesc(dRSA_ccw);
    set(gca, 'YDir', 'normal');
    axis image;
    caxis([cMin cMax]);
    colorbar;
    title(sprintf('%s CCW (n=%d)', char(opts.TitlePrefix), diag.counts(1)));
    xlabel('Time in direction dot1 (samples)');
    ylabel('Time in direction dot1 (samples)');
    xticks(tickCandidates);
    yticks(tickCandidates);
    hold on;
    plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
    hold off;

    nexttile;
    imagesc(dRSA_cw);
    set(gca, 'YDir', 'normal');
    axis image;
    caxis([cMin cMax]);
    colorbar;
    title(sprintf('%s CW (n=%d)', char(opts.TitlePrefix), diag.counts(2)));
    xlabel('Time in direction dot1 (samples)');
    ylabel('Time in direction dot1 (samples)');
    xticks(tickCandidates);
    yticks(tickCandidates);
    hold on;
    plot(1:nTime, 1:nTime, 'w-', 'LineWidth', 1);
    hold off;
end

function dRSAmat = compute_direction_dRSA(dot1Trials, dot2Trials, params)
% compute_direction_dRSA
% Helper that recomputes the direction dRSA map for a subset of trials.

    %% Build direction models from trial subsets
    % Data flow: positions -> displacement -> angles -> direction vectors.
    trialLen = size(dot1Trials, 3);

    dot1Dx = diff(dot1Trials(:, 1, :), 1, 3);
    dot1Dy = diff(dot1Trials(:, 2, :), 1, 3);
    dot1Dx = cat(3, dot1Dx(:, :, 1), dot1Dx);
    dot1Dy = cat(3, dot1Dy(:, :, 1), dot1Dy);
    dot1Angle = atan2(dot1Dy, dot1Dx);
    dot1Direction = cat(2, cos(dot1Angle), sin(dot1Angle));

    dot2Dx = diff(dot2Trials(:, 1, :), 1, 3);
    dot2Dy = diff(dot2Trials(:, 2, :), 1, 3);
    dot2Dx = cat(3, dot2Dx(:, :, 1), dot2Dx);
    dot2Dy = cat(3, dot2Dy(:, :, 1), dot2Dy);
    dot2Angle = atan2(dot2Dy, dot2Dx);
    dot2Direction = cat(2, cos(dot2Angle), sin(dot2Angle));

    %% Concatenate trials and prepare subsamples
    % Data flow: trial arrays -> concatenated series -> trial-locked subsamples.
    modelDirectionDot1 = dRSA_concatenate(dot1Direction, [], 0, 'suppressDispText', 1);
    modelDirectionDot2 = dRSA_concatenate(dot2Direction, [], 0, 'suppressDispText', 1);
    dataDirection = modelDirectionDot1;

    totalTime = size(dataDirection, 2);
    if mod(totalTime, trialLen) ~= 0
        error('Total time (%d) is not a multiple of trial length (%d).', totalTime, trialLen);
    end
    trialStarts = 1:trialLen:totalTime;
    maskSubsampling = true(1, totalTime);
    maskTrigger = false(1, totalTime);
    maskTrigger(trialStarts) = true;

    opt.PreTrigger = 0;
    opt.PostTrigger = trialLen - 1;
    opt.spacing = 0;
    opt.nSubSamples = numel(trialStarts);
    opt.nIter = 1;
    opt.checkRepetition = 0;
    subsamples = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt);

    %% Run dRSA for direction data with direction models
    % Data flow: direction data/models + subsamples -> dRSA time-by-time map.
    paramsGroup = params;
    paramsGroup.nIter = 1;
    paramsGroup.modelToTest = [1 2];
    paramsGroup.modelDistMeasure = {'cosine', 'cosine'};
    model = {modelDirectionDot1, modelDirectionDot2};

    CurrSubsamples = subsamples(:, :, 1);
    dRSAmatFull = dRSA_coreFunction(dataDirection, model, paramsGroup, ...
        'CurrSubsamples', CurrSubsamples, 'suppressDispText', 1);
    dRSAmat = dRSAmatFull(:, :, 1);
end
