
% PIPELINE
% preparatory steps (done by individ user):
% - smoothing data etc
% - concatenate runs (e.g. for movie)
% - average across runs/repetitions
% - create mask (e.g. on TPs with too many NaN repetitions)

data = rand(5,20,10000); 
model = data(:,1:4,:);

%% pepare data and models
%[data mask] = dRSA_concatenate(data,mask) %see documentation with help dRSA_concatenate
[data mask] = dRSA_concatenate(data,mask); %see documentation with help dRSA_concatenate
[model] = dRSA_concatenate(model); %see documentation with help dRSA_concatenate


% optional module: create subsamples
  maskSubsampling = logical(mask.mask); % 1×N logical vector (true = unavailable, false = available) (ALLOW ALSO NUMERICAL)
  opt.SubSampleDur = 500; % : number of points per subsample (ALLOW ALSO SECONDS?)
  opt.spacing  = 10; %   : required spacing between subsamples (in points)
  opt.nSubSamples = 10; % : number of subsamples per iteration
  opt.nIter   = 50; %    : number of iterations (repetitions)
  opt.checkRepetition = 0; % (optional): true/false (default = true)

subsamples = dRSA_random_subsampling(maskSubsampling, opt);
% include illustration
% mask is maskTypes*time
% output: nSubsamples*subSampleDuration*iterations
% for later: option to provide predefined time points (e.g. for predictable vs. unpredictable time points)


%% Simulations



%% dRSA

Y = model.*(-1); % data;
model = {model};
model{2} = model{1};

% wrapper for many subsamples (most cases, except e.g. Ayman)


% In case we do the PCR, it is better to calculate the border outside, because otherwise we would need to recalculate it
% for each iteration
params.nIter = 10;  %how many Iterations?
params.AverageTime = 2; %in s
params.fs = 100; %framerate or how many samples fit into 1 second
params.modelToTest = [1 2];  %array of models to test
params.Var = 0.1; % how much variance? 
params.modelDistMeasure = {'euclidean', 'cosine'};

params.dRSAtype =  'PCR';% 'corr'; %

Autocorrborder = [];
if ~strcmp(params.dRSAtype, 'corr') % Autocorrelation not used with 'corr' type.
    Autocorrborder = dRSA_border(model, subsamples, params);
end

%For the PCR
params.modeltoRegressout = {2, 1}; % {[4 5] [1 3] [1 2]};


%the other stuff use default values. Can be changed, see documentation of PCR function

dRSA_Iter = [];

for iIter = 1:params.nIter
    CurrSubsamples = subsamples(:,:,iIter); % subsamples is nSubsamples*subSampleDuration*iterations
    dRSAma = dRSA_coreFunction(Y,model, params, ...
        'CurrSubsamples', CurrSubsamples, 'Autocorrborder', Autocorrborder);%
    % Y is features*time or a finished RDM of subsamples x subsamples x time
    % models: 1*nModels cell arrray (also used for regress out models, automatically regresses out other models)
    
    
    %In this funciton we have default values
    % we also have varagin: If we want, we can add the autocorrelation and already Subsamples at the end, but we dont need to
    
    
    dRSA_Iter(iIter,:,:,:) = dRSAma;
end  %of nIter

% average dRSAmats (do by summing current with previous summed)  
 % ----------------------------------QUESTION:  include this into averaging function?
dRSA = mean(dRSA_Iter ,1); % 1 = Fmean across first dim
dRSA = reshape(dRSA, size(dRSA,2), size(dRSA,3), size(dRSA,4));


% - module for averaging across Time
params.AverageTime = 2; %in s, how much should be left and right of the zerolag middle? 
params.fs = 100; %framerate or how many samples fit into 1 second
dRSA_diagonal = dRSA_average(dRSA, params);



%% stats and plots

% - module for lag plot
figure;
plot(dRSA_diagonal(1,:));