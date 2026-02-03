function [SSIndices, SSIndicesPlot, usedIterations] = dRSA_triggered_subsampling(maskSubsampling, maskTrigger, opt, varargin)
% dRSA_triggered_subsampling
% Time-lock subsamples to specific trigger points (e.g. word onsets).
%
% INPUTS:
%   maskSubsampling : 1×N logical vector (true = available)
%   maskTrigger     : 1×N logical vector (true = trigger point)
%   opt.PreTrigger  : # timepoints before trigger to include
%   opt.PostTrigger : # timepoints after trigger to include
%   opt.nSubSamples : # subsamples per iteration
%   opt.nIter       : # iterations
%   opt.spacing     : optional minimal spacing between trigger windows (default = 0)
%   opt.checkRepetition : optional (default = true)
%
% OUTPUTS:
%   SSIndices      : [nSubSamples × SubSampleDur × nIter] matrix of indices
%   SSIndicesPlot  : [nIter × N] logical mask of selected points
%   usedIterations : number of unique valid iterations
%
% Name/value options:
%   'suppressDispText' : suppress console output (default = 0)
%
% -------------------------------------------------------------------------
% Based on: dRSA_random_subsampling
% -------------------------------------------------------------------------

%% --- Parameter setup ---
if nargin < 3
    error('You must provide maskSubsampling, maskTrigger, and opt.');
end
if ~islogical(maskSubsampling) || ~islogical(maskTrigger)
    error('Both maskSubsampling and maskTrigger must be logical.');
end
if ~isfield(opt, 'checkRepetition'), opt.checkRepetition = true; end
if ~isfield(opt, 'spacing'), opt.spacing = 0; end

% Data flow: name/value options -> suppressDispText -> status output gating.
parser = inputParser;
parser.addParameter('suppressDispText', 0, ...
    @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
parser.parse(varargin{:});
suppressDispText = logical(parser.Results.suppressDispText);

PreTrigger  = opt.PreTrigger;
PostTrigger = opt.PostTrigger;
nSubSamples = opt.nSubSamples;
nIter       = opt.nIter;
spacing     = opt.spacing;
checkRepetition = opt.checkRepetition;

SubSampleDur = PreTrigger + PostTrigger + 1;
N = numel(maskSubsampling);

rng('shuffle'); % randomize RNG seed

%% --- Step 1: Identify all trigger points
allTriggers = find(maskTrigger);
nTriggersTotal = numel(allTriggers);
if ~suppressDispText
    fprintf('\nFound %d total trigger points.\n', nTriggersTotal);
end

%% --- Step 2: Check which triggers are usable (fully within available data)
validTriggers = [];
for t = allTriggers(:)'
    startIdx = t - PreTrigger;
    endIdx   = t + PostTrigger;
    if startIdx < 1 || endIdx > N
        continue; % outside data range
    end
    if all(maskSubsampling(startIdx:endIdx))
        validTriggers(end+1) = t; %#ok<AGROW>
    end
end

nValid = numel(validTriggers);
if ~suppressDispText
    fprintf('%d trigger points are valid for subsampling (%.1f%% of total).\n', ...
        nValid, 100 * nValid / nTriggersTotal);
end

if nValid == 0
    error('No valid triggers found: check masks or Pre/PostTrigger parameters.');
end

%% --- Step 3: Initialize outputs
SSIndices = nan(nSubSamples, SubSampleDur, nIter);
SSIndicesPlot = false(nIter, N);
usedIterations = 0;

% Optional repetition tracker
if checkRepetition
    seen = containers.Map('KeyType','char','ValueType','logical');
    repeatedCount = 0;
    repeatedSubsequent = 0;
end

%% --- Step 4: Generate iterations
for iter = 1:nIter
    success = false;
    attempt = 0;
    maxAttempts = 5000;

    while ~success && attempt < maxAttempts
        attempt = attempt + 1;

        availableTriggers = validTriggers;
        selectedTriggers = zeros(1, nSubSamples);

        for s = 1:nSubSamples
            if isempty(availableTriggers)
                break;
            end
            tIdx = availableTriggers(randi(numel(availableTriggers)));
            selectedTriggers(s) = tIdx;

            % Remove triggers too close to this one (based on window+spacing)
            excludeRange = (tIdx - (PreTrigger + spacing)) : (tIdx + PostTrigger + spacing);
            excludeRange = excludeRange(excludeRange >= 1 & excludeRange <= N);
            availableTriggers = setdiff(availableTriggers, excludeRange);
        end

        % --- Check if we got enough subsamples
        if all(selectedTriggers)
            selectedTriggers = sort(selectedTriggers);

            % --- Repetition check
            if checkRepetition
                key = mat2str(selectedTriggers);
                if isKey(seen, key)
                    repeatedCount = repeatedCount + 1;
                    repeatedSubsequent = repeatedSubsequent + 1;
                    if repeatedSubsequent >= 100
                        warning('Iteration %d: too many repeated combinations, stopping early.', iter);
                        break;
                    end
                    continue;
                else
                    seen(key) = true;
                    repeatedSubsequent = 0;
                end
            end

            % --- Compute indices
            for s = 1:nSubSamples
                t = selectedTriggers(s);
                SSIndices(s, :, iter) = (t - PreTrigger):(t + PostTrigger);
            end

            % --- Visualization mask
            idxAll = SSIndices(:, :, iter);
            SSIndicesPlot(iter, idxAll(:)) = true;

            success = true;
            usedIterations = usedIterations + 1;
        end
    end

    if ~success
        warning('Iteration %d: could not generate valid set after %d attempts.', iter, maxAttempts);
        break;
    end
end

if ~suppressDispText
    fprintf('\nGenerated %d unique subsampling iterations (requested %d)\n', ...
        usedIterations, nIter);
end

%% --- Step 5: Summary
if ~suppressDispText
    fprintf('Subsample duration: %d points (from %d before to %d after trigger)\n', ...
        SubSampleDur, PreTrigger, PostTrigger);
    fprintf('Average trigger spacing constraint: %d points\n', spacing);
end

%% --- Step 6: Post-hoc repetition check (if checkRepetition = false)
if ~checkRepetition
    allStarts = nan(nSubSamples, nIter);
    for i = 1:nIter
        if any(isnan(SSIndices(1,:,i)))
            continue
        end
        allStarts(:, i) = squeeze(SSIndices(:,ceil(SubSampleDur/2),i)); % center index per subsample
    end
    [~, uniqueIdx] = unique(allStarts.', 'rows');
    nRepetitions = nIter - numel(uniqueIdx);
    if nRepetitions > 0
        warning('Detected %d repeated iterations (%.2f%%) after sampling without repetition check.', ...
            nRepetitions, 100 * nRepetitions / nIter);
    else
        if ~suppressDispText
            fprintf('No repeated iterations detected in post-hoc check.\n');
        end
    end
end

end
