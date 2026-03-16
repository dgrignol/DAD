function RDM = dRSA_computeRDM (data, params, CurrSubsamples, DistMeasure)

% function for calculating RDMS

%% INPUT
%data in format features x time
%DistMeasure: type of dist measure
% CurrSubsamples: a matrix of nSubsamples*subSampleDuration
%params


%% OUTPUT
% matrix in format events x time





%% calculate RDM
nTimePoints = size(CurrSubsamples,2);
nFeatures = size(CurrSubsamples,1) * (size(CurrSubsamples,1) - 1) / 2;
RDM = zeros(nFeatures, nTimePoints);

for iTime = 1:nTimePoints
    if  strcmp(DistMeasure, 'correlation')
        RDM(:,iTime) = dRSA_fastpdist (data(:, CurrSubsamples(:,iTime))', DistMeasure);
    else
        RDM(:,iTime) = pdist(data(:, CurrSubsamples(:,iTime))', DistMeasure);
    end
end

%% Optional debug plot: inspect one time-sample square RDM
% Data flow: vector-form RDM(:, t) -> square trial-by-trial matrix via squareform.
if isfield(params, 'debug') && isstruct(params.debug) && ...
        isfield(params.debug, 'plotRDM') && params.debug.plotRDM

    % Default target is the 100th sample unless overridden by params.debug.sampleIdx.
    if isfield(params.debug, 'sampleIdx')
        debugSampleIdx = params.debug.sampleIdx;
    else
        debugSampleIdx = 100;
    end

    if debugSampleIdx >= 1 && debugSampleIdx <= nTimePoints
        rdmVec = RDM(:, debugSampleIdx);      % condensed pdist vector
        rdmSq  = squareform(rdmVec);          % [nSubsamples x nSubsamples]

        figure('Name', sprintf('DEBUG RDM t=%d', debugSampleIdx), 'NumberTitle', 'off');
        imagesc(rdmSq);
        axis square;
        set(gca, 'YDir', 'normal');
        colorbar;
        title(sprintf('DEBUG RDM at sample %d (%s)', debugSampleIdx, DistMeasure));
        xlabel('Subsample / trial');
        ylabel('Subsample / trial');
    else
        warning('debug.sampleIdx=%d is out of range [1, %d].', debugSampleIdx, nTimePoints);
    end
    
    % Save debug figure only when requested
    if isfield(params.debug, 'saveSnapshoot') && params.debug.saveSnapshoot

        outDir = fullfile(pwd, 'RDMs_snapshots');
        if ~exist(outDir, 'dir'); mkdir(outDir); end

        ts = datestr(now, 'yyyymmdd_HHMMSS_FFF'); % unique timestamp
        baseName = sprintf('RDM_t%03d_%s_%s', debugSampleIdx, DistMeasure, ts);
        outFile = fullfile(outDir, [baseName '.png']);

        % Extra protection against same-name collisions
        k = 1;
        while exist(outFile, 'file')
            outFile = fullfile(outDir, sprintf('%s_%03d.png', baseName, k));
            k = k + 1;
        end

        exportgraphics(gcf, outFile, 'Resolution', 300);
    end

end





%% Normalization

if isfield(params, 'normalization') || strcmp(params.dRSAtype, 'PCR')

    %BEFORE DOING PCR:  all the RDMS should have been normalized here and centered
    if ~isfield(params, 'normalization') && strcmp(params.dRSAtype, 'PCR')
        params.normalization = 'Standardize';  %add a default version if we do PCR
    end

    if strcmp(params.normalization, 'Standardize')
        RDM = dRSA_standardizeRDM (RDM);
    elseif strcmp(params.normalization, 'Rescale')
        RDM = dRSA_rescaleRDM (RDM);
    else
        error('Unknown Normalization method');
    end


    %maybe again call here a new function for normalizing and centering
    %this funciton can be used in dRSA_PCR

end

end