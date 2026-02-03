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
[~, Lag0] = max(Averaged_Autocorr, [], 2);


for iModel = params.modelToTest
    peak = max(Averaged_Autocorr(iModel,:));  %the highest autocorrelation for this model

    %Define the left and right side of the peak
    LeftSide = squeeze(Averaged_Autocorr(iModel,1:Lag0(iModel)));
    RightSide = squeeze(Averaged_Autocorr(iModel,Lag0(iModel):end));
    
    %find the regression borders where less than xx% of variance has been explained
    LeftBorder = find((LeftSide) < regvalue*peak,1, 'last');  %flip to find the first from the left side
    RightBorder = length(RightSide) - find((RightSide) < regvalue*peak, 1, 'first');
    
    border(iModel, :) = ceil(nanmean([LeftBorder RightBorder ])); %take the average
    
     if isnan(border(iModel,:))
        if ~suppressDispText
            fprintf('ERROR: For Model "%d" the regression border could not be calculated! \n', iModel);
        end
    end

    % Optional visualization for diagnostics (compatible with older scripts).
    if plotAutocorr
        time_idx = (1:size(Averaged_Autocorr, 2)) - Lag0(iModel);
        figure;
        % Data flow: mRSA(:,:,iModel) -> imagesc subplot for quick inspection.
        subplot(1, 2, 1);
        imagesc(mRSA(:, :, iModel));
        axis tight;
        colorbar;
        set(gca, 'YDir', 'normal');
        hold on;
        % Diagonal overlays: lag borders are diagonal offsets in the matrix.
        nLag = size(mRSA, 1);
        diag_offsets = [0];
        if ~isempty(LeftBorder)
            diag_offsets(end+1) = LeftBorder - Lag0(iModel);
        end
        if ~isempty(RightBorder)
            right_idx = Lag0(iModel) + (length(RightSide) - RightBorder) - 1;
            diag_offsets(end+1) = right_idx - Lag0(iModel);
        end
        for iOff = 1:numel(diag_offsets)
            offset = diag_offsets(iOff);
            [x_line, y_line] = local_diag_line(nLag, offset);
            if offset == 0
                line(x_line, y_line, 'LineStyle', ':', 'Color', [0.2 0.2 0.2]);
            elseif offset < 0
                line(x_line, y_line, 'LineStyle', '--', 'Color', [0.8 0.2 0.2]);
            else
                line(x_line, y_line, 'LineStyle', '--', 'Color', [0.2 0.2 0.8]);
            end
        end
        hold off;
        title(sprintf('mRSA with borders (model %d)', iModel));
        xlabel('Lag (samples)');
        ylabel('Lag (samples)');

        % Data flow: Averaged_Autocorr -> autocorr plot with borders.
        subplot(1, 2, 2);
        plot(time_idx, Averaged_Autocorr(iModel, :), 'k', 'LineWidth', 1.5);
        hold on;
        yline(regvalue * peak, '--', 'Color', [0.5 0.5 0.5]);
        xline(0, ':', 'Color', [0.2 0.2 0.2]);
        if ~isempty(LeftBorder)
            xline(LeftBorder - Lag0(iModel), '--', 'Color', [0.8 0.2 0.2]);
        end
        if ~isempty(RightBorder)
            right_idx = Lag0(iModel) + (length(RightSide) - RightBorder) - 1;
            xline(right_idx - Lag0(iModel), '--', 'Color', [0.2 0.2 0.8]);
        end
        hold off;
        title(sprintf('Autocorr with borders (model %d)', iModel));
        xlabel('Lag (samples, centered at 0)');
        ylabel('Autocorrelation');
    end
end






end
