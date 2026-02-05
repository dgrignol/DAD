function border  = dRSA_border (model, subsamples, params, varargin)

%% DRSA_BORDER Compute regression borders from model autocorrelation.
% Usage example:
%   border = dRSA_border(model, subsamples, params, ...
%       'suppressDispText', 1, 'plotAutocorr', 'on');
%
% Inputs:
%   model: 1*nModels cell array (also used for regress out models,
%          automatically regresses out other models).
%   subsamples: nSubsamples*subSampleDuration*nIter numeric array.
%   params:
%     .modelToTest = [1 4 5] or cell {1 4 5}
%     .AverageTime = 2; % in s
%     .Var % variance threshold factor
%     .AverageTime % used for calculating the border for that interval
%     .modelNames (optional) cell array of model labels for diagnostics
%
% Name/value options:
%   'suppressDispText' : suppress console output (default = 0)
%   'plotAutocorr'     : plot autocorrelation with borders and mRSA matrix
%                        subplots (default = 'off')
%
% Output:
%   border: In samples the distance from which we regress out our own model;
%           models we are not interested in are NaN.

%% Preparation
% Data flow: name/value options -> local flags -> console/plot gating.
parser = inputParser;
parser.addParameter('suppressDispText', 0, ...
    @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
parser.addParameter('plotAutocorr', 'off', ...
    @(x) (ischar(x) || isstring(x) || isscalar(x)));
parser.parse(varargin{:});
suppressDispText = logical(parser.Results.suppressDispText);
if ischar(parser.Results.plotAutocorr) || isstring(parser.Results.plotAutocorr)
    plotAutocorr = strcmpi(parser.Results.plotAutocorr, 'on');
else
    plotAutocorr = logical(parser.Results.plotAutocorr);
end

function [x_line, y_line] = local_diag_line(nLag, offset)
% LOCAL_DIAG_LINE Return line coords for a diagonal offset in a square matrix.
% offset > 0: above main diagonal (y = x + offset).
% offset < 0: below main diagonal (y = x + offset).
if offset >= 0
    x_line = 1:(nLag - offset);
    y_line = (1 + offset):nLag;
else
    x_line = (1 - offset):nLag;
    y_line = 1:(nLag + offset);
end
end
border = NaN(length(model), 1);
regvalue = sqrt(params.Var ); %apparently, explains 10% of variance ?
%We identify 

%if cell, transform to array
if iscell(params.modelToTest)
    params.modelToTest = [params.modelToTest{:}];  % flatten cell into numeric array
end
%% Autocorrelation
% Data flow: per model/iteration -> dRSA_coreFunction -> dRSA_Iter -> mRSA.
% First we calculate the autocorrelation for all the models.
params.modelDistMeasursaved = params.modelDistMeasure; % save it for later 


for iModel = params.modelToTest
    
    for iIter = 1:params.nIter
        CurrSubsamples = subsamples(:,:,iIter); % subsamples is nSubsamples*subSampleDuration*iterations
        params.dRSAtype = 'corr';
        
        params.modelDistMeasure = params.modelDistMeasursaved {iModel}; 
        dRSAma = dRSA_coreFunction(model{iModel}, {model{iModel}}, params, ...
            'CurrSubsamples', CurrSubsamples, 'suppressDispText', suppressDispText); % get autocorrelation
        dRSA_Iter(iIter,:,:,:) = dRSAma;
    end
    
    dRSA = mean(dRSA_Iter ,1); % 1 = Fmean across first dim
    dRSA = reshape(dRSA, size(dRSA,2), size(dRSA,3));  %remove singleton dimensions
    
    mRSA(:, :, iModel) = dRSA;  %we leave those we do not want empty!
    
end % of iModel

%% Average across time
% Data flow: mRSA -> dRSA_average -> Averaged_Autocorr (model x time).
Averaged_Autocorr = dRSA_average(mRSA, params); 

%% Identify borders
% Data flow: Averaged_Autocorr -> borders per model -> border output.

%find where Lag = 0
% Use the middle sample as lag-0 reference for border extraction/plotting.
lag0Idx = (size(Averaged_Autocorr, 2) + 1) / 2;

%% Optional shared color limits for autocorr plots
% Data flow: mRSA + modelToTest -> global min/max -> shared caxis for comparability.
sharedCaxis = [];
if plotAutocorr
    modelList = params.modelToTest(:)';
    sharedCaxis = [min(mRSA(:, :, modelList), [], 'all', 'omitnan'), ...
        max(mRSA(:, :, modelList), [], 'all', 'omitnan')];
end

for iModel = params.modelToTest
    peak = max(Averaged_Autocorr(iModel,:));  %the highest autocorrelation for this model

    %Define the left and right side of the peak
    LeftSide = squeeze(Averaged_Autocorr(iModel, 1:lag0Idx));
    RightSide = squeeze(Averaged_Autocorr(iModel, lag0Idx:end));
    
    %find the regression borders where less than xx% of variance has been explained
    LeftBorder = find((LeftSide) < regvalue*peak,1, 'last');  %flip to find the first from the left side
    RightBorder = length(RightSide) - find((RightSide) < regvalue*peak, 1, 'first');
    
    border(iModel, :) = ceil(nanmean([LeftBorder RightBorder ])); %take the average
    
     if isnan(border(iModel,:))
        if ~suppressDispText
            fprintf('ERROR: For Model "%d" the regression border could not be calculated! \n', iModel);
        end
    end
    if ~suppressDispText
        % Data flow: borders + lag span -> seconds + window metrics for logging.
        leftLagSamples = LeftBorder - lag0Idx;
        rightIdx = lag0Idx + (length(RightSide) - RightBorder) - 1;
        rightLagSamples = rightIdx - lag0Idx;
        leftLagSec = leftLagSamples / params.fs;
        rightLagSec = rightLagSamples / params.fs;
        windowSizeSamples = rightLagSamples - leftLagSamples;
        totalWindowSamples = size(Averaged_Autocorr, 2) - 1;
        windowSizeSec = windowSizeSamples / params.fs;
        percTotWind = 100 * (windowSizeSamples / totalWindowSamples);
        modelLabel = '';
        if isfield(params, 'modelNames') ...
                && numel(params.modelNames) >= iModel ...
                && ~isempty(params.modelNames{iModel})
            modelLabel = sprintf(' (%s)', params.modelNames{iModel});
        end
        fprintf('dRSA_border: model %d%s LeftLag=%.6g RightLag=%.6g (s), ', ...
            iModel, modelLabel, leftLagSec, rightLagSec);
        fprintf('WindowSize=%.6g (s), PercTotWind=%.2f%%, ', ...
            windowSizeSec, percTotWind);
        % peak is often ~1 in autocorr outputs, so yline may be constant across models.
        fprintf('corrThreshold=%.6g (autocorr threshold), peak=%.6g\n', regvalue * peak, peak);
    end

    % Optional visualization for diagnostics (compatible with older scripts).
    if plotAutocorr
        % Data flow: matrix indices -> seconds axis starting at 0 -> plotted diagnostics.
        nLag = size(mRSA, 1);
        lag_idx = (1:size(Averaged_Autocorr, 2)) - lag0Idx;
        lag_idx_sec = lag_idx / params.fs;
        figure;
        plotTitleSuffix = '';
        if isfield(params, 'modelNames') ...
                && numel(params.modelNames) >= iModel ...
                && ~isempty(params.modelNames{iModel})
            plotTitleSuffix = sprintf(' (%s)', params.modelNames{iModel});
        end
        sgtitle(sprintf('Diagnostics for model %d%s', iModel, plotTitleSuffix));
        % Data flow: mRSA(:,:,iModel) -> imagesc subplot for quick inspection.
        subplot(1, 2, 1);
        time_samples = 0:(nLag - 1);
        time_seconds = time_samples / params.fs;
        imagesc(time_seconds, time_seconds, mRSA(:, :, iModel));
        axis tight;
        colorbar;
        if ~isempty(sharedCaxis)
            caxis(sharedCaxis);
        end
        set(gca, 'YDir', 'normal');
        hold on;
        % Diagonal overlays: lag borders are diagonal offsets in the matrix.
        diag_offsets = [0];
        if ~isempty(LeftBorder)
            diag_offsets(end+1) = LeftBorder - lag0Idx;
        end
        if ~isempty(RightBorder)
            right_idx = lag0Idx + (length(RightSide) - RightBorder) - 1;
            diag_offsets(end+1) = right_idx - lag0Idx;
        end
        for iOff = 1:numel(diag_offsets)
            offset = diag_offsets(iOff);
            [x_line, y_line] = local_diag_line(nLag, offset);
            x_line = time_seconds(x_line);
            y_line = time_seconds(y_line);
            if offset == 0
                line(x_line, y_line, 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);
            elseif offset < 0
                line(x_line, y_line, 'LineStyle', '--', 'Color', [0.8 0.2 0.2]);
            else
                line(x_line, y_line, 'LineStyle', '--', 'Color', [0.2 0.2 0.8]);
            end
        end
        hold off;
        title('mRSA with borders');
        xlabel('Model time (s)');
        ylabel('Model time (s)');

        % Data flow: Averaged_Autocorr -> autocorr plot with borders.
        subplot(1, 2, 2);
        plot(lag_idx_sec, Averaged_Autocorr(iModel, :), 'k', 'LineWidth', 1.5);
        hold on;
        yline(regvalue * peak, '--', 'Color', [0.5 0.5 0.5]);
        xline(0, ':', 'Color', [0.2 0.2 0.2]);
        if ~isempty(LeftBorder)
            xline((LeftBorder - lag0Idx) / params.fs, '--', 'Color', [0.8 0.2 0.2]);
        end
        if ~isempty(RightBorder)
            right_idx = lag0Idx + (length(RightSide) - RightBorder) - 1;
            xline((right_idx - lag0Idx) / params.fs, '--', 'Color', [0.2 0.2 0.8]);
        end
        hold off;
        title('Autocorr with borders');
        xlabel('Lag (s, centered at Lag0)');
        ylabel('Autocorrelation');
    end
end






end
