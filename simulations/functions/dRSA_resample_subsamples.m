function subsamplesOut = dRSA_resample_subsamples(subsamplesIn, nIter, percentKeep)
% dRSA_resample_subsamples
%
% Purpose:
%   Resample subsamples with replacement across iterations to build a
%   N × A × Y tensor, where A is a percentage of the original Z subsamples.
%
% Example usage:
%   % subsamples: Z × Y (or Z × Y × 1) from dRSA_triggered_subsampling
%   subsamplesOut = dRSA_resample_subsamples(subsamples, 50, 50);
%
% Inputs:
%   subsamplesIn : Z × Y or Z × Y × 1 numeric array of subsample indices
%   nIter        : number of resampling iterations (N)
%   percentKeep  : percentage of Z to keep per iteration (B, 0-100]
%
% Output:
%   subsamplesOut : A × Y × N numeric array, where A = round(B/100 * Z)
%
% Key assumptions:
%   - subsamplesIn contains integer indices into the time axis.
%   - Resampling is done with replacement across the Z subsamples.
%   - If subsamplesIn is 3D with a singleton 3rd dimension, it is squeezed.

    %% Input validation and reshape
    if nargin < 3
        error('Provide subsamplesIn, nIter, and percentKeep.');
    end
    if ~isnumeric(subsamplesIn)
        error('subsamplesIn must be numeric.');
    end
    if ~isscalar(nIter) || nIter < 1 || nIter ~= round(nIter)
        error('nIter must be a positive integer.');
    end
    if ~isscalar(percentKeep) || percentKeep <= 0 || percentKeep > 100
        error('percentKeep must be in (0, 100].');
    end

    if ndims(subsamplesIn) == 3 && size(subsamplesIn, 3) == 1
        subsamplesIn = subsamplesIn(:, :, 1);
    end
    if ndims(subsamplesIn) ~= 2
        error('subsamplesIn must be Z × Y (or Z × Y × 1).');
    end

    %% Determine output size and resample
    [zCount, yCount] = size(subsamplesIn);
    aCount = max(1, round((percentKeep / 100) * zCount));
    subsamplesOut = zeros(aCount, yCount, nIter);

    % Resample Z rows with replacement per iteration.
    for iIter = 1:nIter
        drawIdx = randi(zCount, [aCount, 1]);
        subsamplesOut(:, :, iIter) = subsamplesIn(drawIdx, :);
    end
end
