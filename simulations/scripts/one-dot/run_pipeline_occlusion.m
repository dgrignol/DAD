% run_pipeline_occlusion.m
%
% Purpose:
%   Batch runner for one-dot occlusion simulations using
%   PE_simulation_diff_1Dot_occlusion_v1.m across all supported PCR
%   regression strategies.
%
% Usage example (edit participants first, then run):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/one-dot');
%   run_pipeline_occlusion
%
% Usage example (non-interactive):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/simulations/scripts/one-dot'); run('run_pipeline_occlusion.m');"
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX.mat
%   - experiment/input_files/MovDot_SubXX_predicted.mat
%
% Outputs:
%   - simulations/output/SubXX_oneDot_occlusion/<condition>/PCR/...
%   for each configured participant and PCR strategy.
%
% Key assumptions:
%   - The occlusion script accepts:
%       inputConditions = {'occluded_nondeviant','occluded_deviant'}
%   - existingResultsAction = 2 is desired for batch reruns (overwrite).

%% Configure participants and strategies
% Data flow: user-edited lists -> nested loop -> one simulation run per pair.
if ~exist('participants', 'var') || isempty(participants)
    participants = 70;
end
if ~exist('strategies', 'var') || isempty(strategies)
    strategies = { ...
        'baseline_pcr_border', ...
        'ridge_full_autocorr', ...
        'ridge_tapered_autocorr', ...
        'ar1_prewhite_ridge'};
end

%% Resolve target simulation script
% Data flow: current file location -> sibling occlusion simulation script path.
thisDir = fileparts(mfilename('fullpath'));
scriptPath = fullfile(thisDir, 'PE_simulation_diff_1Dot_occlusion_v1.m');
if ~isfile(scriptPath)
    error('Target simulation script not found: %s', scriptPath);
end

%% Run batch
% Data flow: participants x strategies -> parameterized simulation invocations.
for iP = 1:numel(participants)
    for iS = 1:numel(strategies)
        run_one_dot_occlusion_with_explicit_params(participants(iP), strategies{iS}, scriptPath);
    end
end

%% Local worker
function run_one_dot_occlusion_with_explicit_params(pid, strategy, scriptPath)
% RUN_ONE_DOT_OCCLUSION_WITH_EXPLICIT_PARAMS Run one participant/strategy pair.
participantNumber = pid;
dRSAtypeToRun = 'PCR'; %#ok<NASGU>
inputConditions = {'occluded_nondeviant', 'occluded_deviant'}; %#ok<NASGU>
existingResultsAction = 2; %#ok<NASGU>
suppressDispText = 0; %#ok<NASGU>
plotSampleRateHz = 120; %#ok<NASGU>
cutPostDev = false; %#ok<NASGU>
cutPostDevPlot = false; %#ok<NASGU>

params = struct();
params.PCR = struct();
params.PCR.RegressStrategy = strategy;

if strcmp(strategy, 'ar1_prewhite_ridge')
    params.PCR.RidgeLambdaFactor = 0.04;
else
    params.PCR.RidgeLambdaFactor = 0.06;
end

params.PCR.TaperSigmaFactor = 0.33;
params.PCR.AutoPenaltyStrength = 2;
params.PCR.NormalizeAutoPenaltyPerRow = 1;
params.PCR.StandardizePredictors = 1;
params.PCR.AR1Clip = 0.99;
assert(isstruct(params) && isfield(params, 'PCR')); % explicit use for analyzer clarity.

fprintf('Running participant %02d with strategy %s...\n', participantNumber, strategy);
run(scriptPath);
end
