% 07_dRSA.m
%
% Purpose:
%   Run real-neural dRSA for one-dot trajectories using identity-averaged MEG
%   exports and subject-matched observed/predicted trajectory models.
%
% Data-flow summary:
%   1) Load neural tensors and channel labels for:
%      always_visible, occluded_nondeviant, occluded_deviant.
%   2) Check that channel/time dimensions and channel order are consistent
%      across conditions on the exported canonical channel axis.
%   3) Remove channels that are all-NaN over identities x time in any
%      condition, and print removed channel labels to console.
%   4) Load observed and predicted trajectory model files for the same subject.
%   5) Validate sequence order (must be identity 1..20), then build:
%      observed position/direction and predicted position/direction.
%   6) Build PE models in occluded_deviant as observed - predicted.
%   7) Run dRSA (corr or PCR) with neural euclidean and neural cosine rows.
%   8) Save matrix figures and MAT results for each condition.
%
% Usage examples:
%   % Example 1: run with defaults (sub-01, ses-01, corr, 120 Hz).
%   eval(fileread('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/analysis/07_dRSA.m'));
%
%   % Example 2: override selected settings before running.
%   subj_str = '01';
%   dRSAtypeToRun = 'corr'; % or 'PCR'
%   plotSampleRateHz = 120;
%   eval(fileread('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/analysis/07_dRSA.m'));
%
% Inputs:
%   Neural MAT files:
%   /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/analysis/
%   sub-XX/ses-YY/drsa/sub-XX_ses-YY_task-Attention_condition-<condition>_desc-identityAveraged_neural.mat
%
%   Observed/predicted model files:
%   /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/input_files/subXX/
%   MovDot_SubXX_V28_rescueTraject.mat
%   MovDot_SubXX_V28_rescueTraject_predicted.mat
%
% Required variables and dimensions:
%   - neural: [20 x nChannels x nTime] numeric
%   - channel_labels: [nChannels x 1] or [1 x nChannels], same order in all
%     conditions; may include bad/missing channels exported as all-NaN
%   - xySeqs / xySeqsPredicted entries: condition_label, sequence, xy
%   - sequence order must be exactly 1..20 per condition
%   - xy must have size [nTime x 2] to match neural time axis
%
% Outputs:
%   Workspace variables:
%     - drsa_input.<condition>.neural
%     - drsa_input.<condition>.channel_labels
%     - model_input.<condition>.observed_position
%     - model_input.<condition>.observed_direction
%     - model_input.occluded_deviant.predicted_position
%     - model_input.occluded_deviant.predicted_direction
%     - drsa_output.<condition>.dRSA_neural_correlation
%
%   Saved files:
%   /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/analysis/
%     sub-XX/ses-YY/drsa/matlab_real_neural/07_dRSA/<condition>/<dRSAtype>/
%       matrices/*.png
%       results/*.mat
%
% Assumptions:
%   - Identity codes in neural exports map directly to sequence IDs 1..20 in
%     the model files.
%   - Neural exports start at t=0 and use 120 Hz downsampled timing.
%   - This script is intentionally strict: no fallback model file search.
%
% Notes:
%   - PE models are defined as observed - predicted, matching the simulation
%     convention.
%   - The script is written as a compact, inspectable pipeline, with explicit
%     checks at each handoff between neural and model spaces.

%% Reset workspace state while preserving optional runtime overrides
clearvars -except subj_str dRSAtypeToRun plotSampleRateHz
dRSAtypeToRun = 'PCR';
%% Runtime defaults and canonical input locations
if ~exist('subj_str', 'var') || isempty(subj_str)
    subj_str = '01'; % default subject
end

ses_str = '01'; % DAD MEG currently uses ses-01 only
if ~exist('dRSAtypeToRun', 'var') || isempty(dRSAtypeToRun)
    dRSAtypeToRun = 'corr'; % valid: 'corr' or 'PCR'
end
if ~exist('plotSampleRateHz', 'var') || isempty(plotSampleRateHz)
    plotSampleRateHz = 120; % used for matrix axis labeling
end

script_dir = '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/analysis';
meg_root = fileparts(fileparts(script_dir)); % .../MEG
drsa_dir = fullfile(meg_root, 'derivatives', 'analysis', ...
    sprintf('sub-%s', subj_str), sprintf('ses-%s', ses_str), 'drsa');

conditions = {'always_visible', 'occluded_nondeviant', 'occluded_deviant'};
drsa_input = struct();

%% Load neural tensors from step-06 exports (one file per condition)
for i_cond = 1:numel(conditions)
    cond = conditions{i_cond};
    % File naming is strict to avoid accidental subject/version mismatch.
    mat_file = fullfile(drsa_dir, sprintf( ...
        'sub-%s_ses-%s_task-Attention_condition-%s_desc-identityAveraged_neural.mat', ...
        subj_str, ses_str, cond));
    if ~isfile(mat_file)
        error('Missing input file for condition "%s": %s', cond, mat_file);
    end

    S = load(mat_file, 'neural', 'channel_labels');
    if ~isfield(S, 'neural') || ~isfield(S, 'channel_labels')
        error('File %s must contain both variables: neural, channel_labels.', mat_file);
    end
    if ~isnumeric(S.neural) || ndims(S.neural) ~= 3 || size(S.neural, 1) ~= 20
        error('Invalid neural tensor in %s. Expected numeric [20 x nChannels x nTime].', mat_file);
    end

    drsa_input.(cond) = S;
end

%% Cross-condition consistency checks for neural channel/time axes
ref_cond = conditions{1};
ref_size = size(drsa_input.(ref_cond).neural);
ref_labels = drsa_input.(ref_cond).channel_labels(:);

for i_cond = 2:numel(conditions)
    cond = conditions{i_cond};
    curr_size = size(drsa_input.(cond).neural);
    if ~isequal(curr_size(2:3), ref_size(2:3))
        error(['Channel/time mismatch between %s and %s. ' ...
            '%s: [%d x %d], %s: [%d x %d].'], ...
            ref_cond, cond, ref_cond, ref_size(2), ref_size(3), cond, curr_size(2), curr_size(3));
    end

    curr_labels = drsa_input.(cond).channel_labels(:);
    if numel(curr_labels) ~= numel(ref_labels) || ~isequal(curr_labels, ref_labels)
        error('Channel label mismatch between %s and %s.', ref_cond, cond);
    end
end

%% Remove all-NaN channels from canonical export axis before dRSA
%
% Export step 06 now keeps a fixed canonical channel list and writes bad/missing
% channels as NaN. To preserve historical dRSA behavior, we remove channels that
% are fully NaN across identities and time, then run dRSA on remaining channels.
nan_mask_by_condition = false(numel(conditions), ref_size(2));
for i_cond = 1:numel(conditions)
    cond = conditions{i_cond};
    neural_cond = drsa_input.(cond).neural;
    nan_mask_by_condition(i_cond, :) = squeeze(all(all(isnan(neural_cond), 3), 1));
end

remove_nan_mask = any(nan_mask_by_condition, 1);
keep_mask = ~remove_nan_mask;
removed_channel_labels = cellstr(string(ref_labels(remove_nan_mask)));
kept_channel_labels = cellstr(string(ref_labels(keep_mask)));

fprintf(['Loaded canonical dRSA inputs for sub-%s ses-%s before NaN removal: ' ...
    '%d conditions, %d identities, %d channels, %d timepoints.\n'], ...
    subj_str, ses_str, numel(conditions), ref_size(1), ref_size(2), ref_size(3));
fprintf('All-NaN channels removed before dRSA: %d\n', numel(removed_channel_labels));
if ~isempty(removed_channel_labels)
    fprintf('Removed channel labels (all-NaN in >=1 condition):\n');
    for i_ch = 1:numel(removed_channel_labels)
        fprintf('  %s\n', removed_channel_labels{i_ch});
    end
else
    fprintf('Removed channel labels (all-NaN in >=1 condition): none\n');
end

if ~any(keep_mask)
    error(['All channels are all-NaN after export. Check preprocessing step 06 ' ...
        'and bad-channel marking before running dRSA.']);
end

for i_cond = 1:numel(conditions)
    cond = conditions{i_cond};
    drsa_input.(cond).neural = drsa_input.(cond).neural(:, keep_mask, :);
    drsa_input.(cond).channel_labels = kept_channel_labels(:);
end

ref_size = size(drsa_input.(ref_cond).neural);
fprintf('Channels retained for dRSA after NaN removal: %d\n', ref_size(2));

%% Load subject-specific one-dot model files (strict filenames, no fallback)
dad_root = fileparts(meg_root); % .../DAD
model_dir = fullfile(dad_root, 'experiment', 'input_files', sprintf('sub%s', subj_str));
if ~isfolder(model_dir)
    error('Missing subject model directory: %s', model_dir);
end

obs_file = fullfile(model_dir, sprintf('MovDot_Sub%s_V28_rescueTraject.mat', subj_str));
pred_file = fullfile(model_dir, sprintf('MovDot_Sub%s_V28_rescueTraject_predicted.mat', subj_str));
if ~isfile(obs_file)
    error('Missing observed model file: %s', obs_file);
end
if ~isfile(pred_file)
    error('Missing predicted model file: %s', pred_file);
end

%% Parse observed models and convert to identity-ordered position/direction
Obs = load(obs_file, 'xySeqs', 'Cfg');
if ~isfield(Obs, 'xySeqs') || isempty(Obs.xySeqs)
    error('Observed model file must contain non-empty xySeqs: %s', obs_file);
end
if ~isfield(Obs, 'Cfg') || ~isfield(Obs.Cfg, 'rectSize')
    error('Observed model file must contain Cfg.rectSize: %s', obs_file);
end
obs_rect_size = double(Obs.Cfg.rectSize(:)');
if numel(obs_rect_size) ~= 2
    error('Invalid observed Cfg.rectSize in %s. Expected [width height].', obs_file);
end
obs_center_shift = reshape(obs_rect_size ./ 2, [1, 2, 1]);

Pred = load(pred_file, 'xySeqsPredicted', 'Cfg');
if ~isfield(Pred, 'xySeqsPredicted') || isempty(Pred.xySeqsPredicted)
    error('Predicted model file must contain non-empty xySeqsPredicted: %s', pred_file);
end

identity_axis = 1:20;
n_time_neural = ref_size(3);
model_input = struct();

obs_trials = Obs.xySeqs(:);
if ~isfield(obs_trials, 'condition_label') || ~isfield(obs_trials, 'sequence') || ~isfield(obs_trials, 'xy')
    error('Observed xySeqs entries must contain condition_label, sequence, and xy.');
end
obs_labels = lower(string({obs_trials.condition_label}));

for i_cond = 1:numel(conditions)
    cond = conditions{i_cond};
    % Keep only trials matching the current condition label.
    cond_trials = obs_trials(obs_labels == cond);
    if isempty(cond_trials)
        error('No observed trials found for condition %s in %s.', cond, obs_file);
    end

    % Tiny guard: file order must already be identity 1..20.
    seq_values = round(double([cond_trials.sequence]));
    if ~isequal(seq_values(:)', identity_axis)
        error(['Observed sequence order is not 1..20 for %s. ' ...
            'Current order: %s'], cond, mat2str(seq_values));
    end

    % Fill [identity x feature(2) x time] observed position tensor.
    observed_position = zeros(numel(identity_axis), 2, n_time_neural);
    for i_id = 1:numel(identity_axis)
        xy = double(cond_trials(i_id).xy);
        if size(xy, 2) ~= 2
            error('Expected one-dot xy with 2 columns for condition %s sequence %d.', cond, i_id);
        end
        if size(xy, 1) ~= n_time_neural
            error(['Observed model timepoints mismatch for %s sequence %d. ' ...
                'Model=%d, neural=%d.'], cond, i_id, size(xy, 1), n_time_neural);
        end
        observed_position(i_id, 1, :) = reshape(xy(:, 1), 1, 1, []);
        observed_position(i_id, 2, :) = reshape(xy(:, 2), 1, 1, []);
    end

    % Re-center to stimulus-centered coordinates using half screen size.
    observed_position = bsxfun(@minus, observed_position, obs_center_shift);

    % Direction feature: angle from framewise displacement, encoded as cos/sin.
    dx = diff(observed_position(:, 1, :), 1, 3);
    dy = diff(observed_position(:, 2, :), 1, 3);
    % Prepend first derivative to preserve original time length.
    dx = cat(3, dx(:, :, 1), dx);
    dy = cat(3, dy(:, :, 1), dy);
    observed_direction = cat(2, cos(atan2(dy, dx)), sin(atan2(dy, dx)));

    model_input.(cond).observed_position = observed_position;
    model_input.(cond).observed_direction = observed_direction;
end

%% Parse predicted model for occluded_deviant and build position/direction
pred_trials = Pred.xySeqsPredicted(:);
if ~isfield(pred_trials, 'condition_label') || ~isfield(pred_trials, 'sequence') || ~isfield(pred_trials, 'xy')
    error('Predicted xySeqsPredicted entries must contain condition_label, sequence, and xy.');
end
pred_labels = lower(string({pred_trials.condition_label}));
pred_mask = pred_labels == "occluded_deviant_predicted" | pred_labels == "occluded_deviant";
pred_trials = pred_trials(pred_mask);
if isempty(pred_trials)
    error('No predicted occluded_deviant trials found in %s.', pred_file);
end

% Tiny guard: predicted identities must also be in 1..20 order.
pred_seq_values = round(double([pred_trials.sequence]));
if ~isequal(pred_seq_values(:)', identity_axis)
    error(['Predicted sequence order is not 1..20 for occluded_deviant. ' ...
        'Current order: %s'], mat2str(pred_seq_values));
end

% Fill [identity x feature(2) x time] predicted position tensor.
predicted_position = zeros(numel(identity_axis), 2, n_time_neural);
for i_id = 1:numel(identity_axis)
    xy = double(pred_trials(i_id).xy);
    if size(xy, 2) ~= 2
        error('Expected one-dot predicted xy with 2 columns for sequence %d.', i_id);
    end
    if size(xy, 1) ~= n_time_neural
        error(['Predicted model timepoints mismatch for occluded_deviant sequence %d. ' ...
            'Model=%d, neural=%d.'], i_id, size(xy, 1), n_time_neural);
    end
    predicted_position(i_id, 1, :) = reshape(xy(:, 1), 1, 1, []);
    predicted_position(i_id, 2, :) = reshape(xy(:, 2), 1, 1, []);
end

% Use predicted rectSize if provided; otherwise reuse observed rectSize.
if isfield(Pred, 'Cfg') && isfield(Pred.Cfg, 'rectSize')
    pred_rect_size = double(Pred.Cfg.rectSize(:)');
else
    pred_rect_size = obs_rect_size;
end
if numel(pred_rect_size) ~= 2
    error('Invalid predicted Cfg.rectSize in %s. Expected [width height].', pred_file);
end
pred_center_shift = reshape(pred_rect_size ./ 2, [1, 2, 1]);
predicted_position = bsxfun(@minus, predicted_position, pred_center_shift);

% Predicted direction is built with the same transform as observed.
dx = diff(predicted_position(:, 1, :), 1, 3);
dy = diff(predicted_position(:, 2, :), 1, 3);
dx = cat(3, dx(:, :, 1), dx);
dy = cat(3, dy(:, :, 1), dy);
predicted_direction = cat(2, cos(atan2(dy, dx)), sin(atan2(dy, dx)));

model_input.occluded_deviant.predicted_position = predicted_position;
model_input.occluded_deviant.predicted_direction = predicted_direction;

fprintf(['Loaded model trajectories for sub-%s from:\n  observed: %s\n  predicted: %s\n' ...
    'Aligned to %d identities and %d timepoints.\n'], ...
    subj_str, obs_file, pred_file, numel(identity_axis), n_time_neural);

%% Add and validate dRSA helper functions from simulations/functions
drsa_func_dir = fullfile(dad_root, 'simulations', 'functions');
if ~isfolder(drsa_func_dir)
    error('Missing dRSA helper directory: %s', drsa_func_dir);
end
addpath(drsa_func_dir);

required_funcs = {'dRSA_concatenate', 'dRSA_coreFunction', 'dRSA_triggered_subsampling'};
for i_req = 1:numel(required_funcs)
    if exist(required_funcs{i_req}, 'file') ~= 2
        error('Missing required function on path: %s', required_funcs{i_req});
    end
end

if strcmpi(dRSAtypeToRun, 'PCR')
    if exist('dRSA_border', 'file') ~= 2
        error('PCR mode requires dRSA_border.m on path.');
    end
else
    dRSAtypeToRun = 'corr';
end

%% Run dRSA condition-by-condition and save figures/results
output_root = fullfile(drsa_dir, 'matlab_real_neural', '07_dRSA');
drsa_output = struct();

for i_cond = 1:numel(conditions)
    cond = conditions{i_cond};
    neural_tensor = drsa_input.(cond).neural;

    observed_position = model_input.(cond).observed_position;
    observed_direction = model_input.(cond).observed_direction;

    model_cell = cell(1, 2);
    model_names = {'observed position dot1', 'observed direction dot1'};
    model_dist = {'euclidean', 'cosine'};
    model_cell{1} = dRSA_concatenate(observed_position, [], 0, 'suppressDispText', 1);
    model_cell{2} = dRSA_concatenate(observed_direction, [], 0, 'suppressDispText', 1);

    if strcmp(cond, 'occluded_deviant')
        % Only occluded_deviant has explicit predicted and PE model terms.
        predicted_position = model_input.occluded_deviant.predicted_position;
        predicted_direction = model_input.occluded_deviant.predicted_direction;
        % PE definition follows simulation script convention.
        pe_position = observed_position - predicted_position;
        pe_direction = observed_direction - predicted_direction;

        model_cell{3} = dRSA_concatenate(predicted_position, [], 0, 'suppressDispText', 1);
        model_cell{4} = dRSA_concatenate(predicted_direction, [], 0, 'suppressDispText', 1);
        model_cell{5} = dRSA_concatenate(pe_position, [], 0, 'suppressDispText', 1);
        model_cell{6} = dRSA_concatenate(pe_direction, [], 0, 'suppressDispText', 1);

        model_names = {'observed position dot1', 'observed direction dot1', ...
            'predicted position dot1', 'predicted direction dot1', ...
            'PE position dot1', 'PE direction dot1'};
        model_dist = {'euclidean', 'cosine', 'euclidean', 'cosine', 'euclidean', 'cosine'};
    end

    % Concatenate identity-by-time to the dRSA expected [feature x long-time].
    neural_concat = dRSA_concatenate(neural_tensor, [], 0, 'suppressDispText', 1);
    trial_len = size(neural_tensor, 3);
    total_time = size(neural_concat, 2);

    % Triggered subsampling: one contiguous window per identity epoch.
    trial_starts = 1:trial_len:total_time;
    mask_subsampling = true(1, total_time);
    mask_trigger = false(1, total_time);
    mask_trigger(trial_starts) = true;

    opt = struct();
    opt.PreTrigger = 0;
    opt.PostTrigger = trial_len - 1;
    opt.spacing = 0;
    opt.nSubSamples = numel(trial_starts);
    opt.nIter = 1;
    opt.checkRepetition = 0;
    subsamples = dRSA_triggered_subsampling(mask_subsampling, mask_trigger, opt, 'suppressDispText', 1);

    % Common dRSA parameters, including leave-other-models-out regressors.
    params = struct();
    params.dRSAtype = dRSAtypeToRun;
    params.modelNames = model_names;
    params.modelToTest = 1:numel(model_cell);
    params.modelDistMeasure = model_dist;
    params.neuralDistMeasure = 'euclidean';
    params.modeltoRegressout = cell(1, numel(model_cell));
    for i_model = 1:numel(model_cell)
        params.modeltoRegressout{i_model} = setdiff(1:numel(model_cell), i_model);
    end
    params.nIter = size(subsamples, 3);
    params.fs = plotSampleRateHz;
    params.AverageTime = floor((trial_len - 1) / 2) / params.fs;
    params.Var = 0.1;
    params.PCR = struct();
    params.PCR.RegressStrategy = 'baseline_pcr_border';
    params.PCR.RegressAutocor = 1;
    params.PCR.RessModel = 1;
    params.PCR.AdditionalPCA = 1;

    % PCR computes autocorrelation borders; corr leaves this empty.
    autocorr_border = [];
    if strcmpi(params.dRSAtype, 'PCR')
        autocorr_border = dRSA_border(model_cell, subsamples, params, ...
            'suppressDispText', 1, 'plotAutocorr', 'off');
    end

    % Neural RDM distance is fixed to correlation.
    params.neuralDistMeasure = 'correlation';
    dRSA_neural_correlation = dRSA_coreFunction(neural_concat, model_cell, params, ...
        'CurrSubsamples', subsamples(:, :, 1), ...
        'Autocorrborder', autocorr_border, ...
        'suppressDispText', 1);

    drsa_output.(cond).dRSA_neural_correlation = dRSA_neural_correlation;
    drsa_output.(cond).model_names = model_names;
    drsa_output.(cond).model_dist = model_dist;

    % Create condition-specific output folders.
    cond_root = fullfile(output_root, cond, lower(dRSAtypeToRun));
    matrices_dir = fullfile(cond_root, 'matrices');
    results_dir = fullfile(cond_root, 'results');
    if ~isfolder(matrices_dir), mkdir(matrices_dir); end
    if ~isfolder(results_dir), mkdir(results_dir); end

    time_seconds = (0:size(dRSA_neural_correlation, 1)-1) ./ plotSampleRateHz;

    % Plot a 1 x nModels grid with neural correlation dRSA matrices.
    fig = figure('Name', sprintf('dRSA matrices (%s, %s)', cond, dRSAtypeToRun), ...
        'NumberTitle', 'off', 'Visible', 'on');
    tiledlayout(1, numel(model_names), 'Padding', 'compact', 'TileSpacing', 'compact');
    for i_model = 1:numel(model_names)
        nexttile(i_model);
        imagesc(time_seconds, time_seconds, dRSA_neural_correlation(:, :, i_model));
        set(gca, 'YDir', 'normal');
        axis image;
        colorbar;
        title(model_names{i_model}, 'Interpreter', 'none');
        xlabel('neural time (s)');
        ylabel('model time (s)');
        hold on;
        plot(time_seconds, time_seconds, 'w-', 'LineWidth', 1);
        hold off;
    end

    % Save matrix figure and numeric result payload for reproducibility.
    fig_file = fullfile(matrices_dir, sprintf('dRSA_%s_%s_matrices.png', cond, lower(dRSAtypeToRun)));
    exportgraphics(fig, fig_file, 'Resolution', 300, 'BackgroundColor', 'white');

    results_file = fullfile(results_dir, sprintf('dRSA_%s_%s_results.mat', cond, lower(dRSAtypeToRun)));
    save(results_file, 'dRSA_neural_correlation', 'model_names', 'model_dist', 'params', ...
        'obs_file', 'pred_file', 'subj_str', 'ses_str');

    fprintf('Saved dRSA outputs for %s (%s):\n  %s\n  %s\n', ...
        cond, dRSAtypeToRun, fig_file, results_file);
end

fprintf('Done. All dRSA matrices were generated under:\n%s\n', output_root);
