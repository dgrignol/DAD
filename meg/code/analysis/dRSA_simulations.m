% dRSA_simulations.m
%
% Purpose:
%   Run simulation-style dRSA analyses for one participant using the DAD
%   one-dot trajectories and identity-averaged MEG exports, then write all
%   outputs under the participant DRSA derivatives "simulations" subfolder.
%
%   The script computes two analysis families:
%   1) Model-vs-model simulations per condition:
%      - position model with position model (autocorrelation-like matrix)
%      - position model with direction model
%   2) Neural-vs-neural simulations per condition:
%      - neural data with itself (autocorrelation-like matrix)
%      - run independently for each MEGIN sensor group (no cross-group pairs)
%        wholeBrain, occipital, parietal, frontal, temporal.
%
% Data-flow summary:
%   1) Load one-dot trajectory models (observed; optional predicted branch for
%      occluded_deviant only if requested).
%   2) For each condition, build position and direction model tensors aligned
%      to the neural trial/time axes.
%   3) Build trial-locked subsamples once per condition.
%   4) Run model-vs-model dRSA and save matrices/results.
%   5) Load neural data for the same condition, remove all-NaN channels,
%      optionally smooth, then run neural-vs-neural dRSA per sensor group.
%   6) Save matrices/results to:
%      <drsa_dir>/simulations/...
%
% Usage examples:
%   % Example 1: default run for sub-01/ses-01 and all three conditions.
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/analysis/dRSA_simulations.m');
%
%   % Example 2: run one condition only with PCR model simulations.
%   subj_str = '01';
%   ses_str = '01';
%   conditions_to_run = {'always_visible'};
%   dRSAtypeToRun = 'PCR';
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/analysis/dRSA_simulations.m');
%
%   % Example 3: sandbox test (redirect outputs away from real derivatives).
%   output_root_override = '/tmp/my_drsa_sim_sandbox';
%   run('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/analysis/dRSA_simulations.m');
%
% Inputs:
%   - Neural MAT files (identity-averaged):
%     sub-XX_ses-YY_task-Attention_condition-<condition>_desc-identityAveraged_neural.mat
%     with variables:
%       neural         [nIdentity x nChannels x nTime]
%       channel_labels [nChannels x 1] or [1 x nChannels]
%   - Trajectory model files:
%       MovDot_SubXX_V28_rescueTraject.mat
%       MovDot_SubXX_V28_rescueTraject_predicted.mat
%
% Outputs:
%   - PNG matrices and MAT result files in:
%     /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/derivatives/
%       analysis/sub-XX/ses-YY/drsa/simulations/
%         model_vs_model/<condition>/<dRSAtype>/{matrices,results}
%         neural_vs_neural/<condition>/corr/{matrices,results}
%
% Assumptions:
%   - Trajectory sequence labels map to neural identities (typically 1..20).
%   - Neural and model time axes match for each condition.
%   - All neural-vs-neural runs requested here are autocorrelation-only
%     (group with itself), not cross-group comparisons.

%% Reset workspace while preserving runtime overrides
clearvars -except subj_str ses_str conditions_to_run dRSAtypeToRun plotSampleRateHz ...
    pcrStrategy doSmoothNeural smoothWindowSamples model_path_source ...
    sensor_groups_to_run output_root_override

%% Runtime defaults and configuration
if ~exist('subj_str', 'var') || isempty(subj_str), subj_str = '01'; end
if ~exist('ses_str', 'var') || isempty(ses_str), ses_str = '01'; end
if ~exist('conditions_to_run', 'var') || isempty(conditions_to_run)
    conditions_to_run = {'always_visible', 'occluded_nondeviant', 'occluded_deviant'};
end
if ischar(conditions_to_run) || isstring(conditions_to_run)
    conditions_to_run = cellstr(conditions_to_run);
end
conditions_to_run = cellfun(@(c) lower(strtrim(char(c))), conditions_to_run, 'UniformOutput', false);
if ~exist('dRSAtypeToRun', 'var') || isempty(dRSAtypeToRun), dRSAtypeToRun = 'corr'; end
if ~exist('plotSampleRateHz', 'var') || isempty(plotSampleRateHz), plotSampleRateHz = 120; end
if ~exist('pcrStrategy', 'var') || isempty(pcrStrategy), pcrStrategy = 'baseline_pcr_border'; end
if ~exist('doSmoothNeural', 'var') || isempty(doSmoothNeural), doSmoothNeural = true; end
if ~exist('smoothWindowSamples', 'var') || isempty(smoothWindowSamples), smoothWindowSamples = 15; end
if ~exist('model_path_source', 'var') || isempty(model_path_source), model_path_source = 'observed'; end
if ~exist('sensor_groups_to_run', 'var') || isempty(sensor_groups_to_run)
    sensor_groups_to_run = {'wholeBrain', 'occipital', 'parietal', 'frontal', 'temporal'};
end
if ischar(sensor_groups_to_run) || isstring(sensor_groups_to_run)
    sensor_groups_to_run = cellstr(sensor_groups_to_run);
end

valid_conditions = {'always_visible', 'occluded_nondeviant', 'occluded_deviant'};
if ~all(ismember(conditions_to_run, valid_conditions))
    error('conditions_to_run must be drawn from: %s', strjoin(valid_conditions, ', '));
end
valid_drsa_types = {'corr', 'PCR'};
if ~any(strcmpi(dRSAtypeToRun, valid_drsa_types))
    error('dRSAtypeToRun must be one of: %s', strjoin(valid_drsa_types, ', '));
end
model_path_source = lower(strtrim(char(model_path_source)));
if ~any(strcmp(model_path_source, {'observed', 'predicted'}))
    error('model_path_source must be either "observed" or "predicted".');
end

%% Resolve paths and dependencies
% Data flow: root paths -> input files -> output folders -> helper functions.
dad_root = '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD';
meg_root = fullfile(dad_root, 'MEG');
drsa_dir = fullfile(meg_root, 'derivatives', 'analysis', ...
    sprintf('sub-%s', subj_str), sprintf('ses-%s', ses_str), 'drsa');
model_dir = fullfile(dad_root, 'experiment', 'input_files', sprintf('sub%s', subj_str));

if exist('output_root_override', 'var') && ~isempty(output_root_override)
    base_out_dir = char(output_root_override);
else
    base_out_dir = drsa_dir;
end

obs_file = fullfile(model_dir, sprintf('MovDot_Sub%s_V28_rescueTraject.mat', subj_str));
pred_file = fullfile(model_dir, sprintf('MovDot_Sub%s_V28_rescueTraject_predicted.mat', subj_str));
if ~isfile(obs_file)
    error('Observed model file not found: %s', obs_file);
end
if ~isfile(pred_file)
    error('Predicted model file not found: %s', pred_file);
end

addpath(fullfile(dad_root, 'simulations', 'functions'));

sensor_group_lists = local_get_megin_sensor_groups();
[sensor_group_keys, sensor_group_tags] = local_normalize_sensor_groups(sensor_groups_to_run);
for i_group = 1:numel(sensor_group_keys)
    if ~isfield(sensor_group_lists, sensor_group_keys{i_group})
        error('Unsupported sensor group "%s" after normalization.', sensor_group_keys{i_group});
    end
end

%% Load trajectory sources once
% Data flow: model source files -> reusable trial structs for each condition.
Obs = load(obs_file, 'xySeqs', 'Cfg');
if ~isfield(Obs, 'xySeqs') || isempty(Obs.xySeqs)
    error('Observed file %s must contain non-empty xySeqs.', obs_file);
end
if ~isfield(Obs, 'Cfg') || ~isfield(Obs.Cfg, 'rectSize')
    error('Observed file %s must contain Cfg.rectSize.', obs_file);
end
obs_center_shift = reshape(double(Obs.Cfg.rectSize(:) ./ 2), [1, 2, 1]);

Pred = load(pred_file, 'xySeqsPredicted', 'Cfg');
if ~isfield(Pred, 'xySeqsPredicted') || isempty(Pred.xySeqsPredicted)
    error('Predicted file %s must contain non-empty xySeqsPredicted.', pred_file);
end
if isfield(Pred, 'Cfg') && isfield(Pred.Cfg, 'rectSize')
    pred_center_shift = reshape(double(Pred.Cfg.rectSize(:) ./ 2), [1, 2, 1]);
else
    pred_center_shift = obs_center_shift;
end

fprintf('Running dRSA simulations for sub-%s ses-%s\n', subj_str, ses_str);
fprintf('Output root: %s\n', base_out_dir);

%% Condition loop: model simulations + neural autocorrelation per sensor group
% Data flow: condition -> model tensors + neural tensors -> dRSA outputs.
for i_cond = 1:numel(conditions_to_run)
    condition_str = conditions_to_run{i_cond};

    %% Load and validate neural input for current condition
    neural_file = fullfile(drsa_dir, sprintf( ...
        'sub-%s_ses-%s_task-Attention_condition-%s_desc-identityAveraged_neural.mat', ...
        subj_str, ses_str, condition_str));
    if ~isfile(neural_file)
        error('Missing neural input for condition %s: %s', condition_str, neural_file);
    end

    N = load(neural_file, 'neural', 'channel_labels');
    if ~isfield(N, 'neural') || ~isfield(N, 'channel_labels')
        error('Neural file %s must contain neural and channel_labels.', neural_file);
    end
    neural = N.neural;
    channel_labels = string(N.channel_labels(:));
    if ndims(neural) ~= 3
        error('neural must be 3D [identity x channels x time] in %s.', neural_file);
    end
    if numel(channel_labels) ~= size(neural, 2)
        error('channel_labels (%d) do not match neural channels (%d) in %s.', ...
            numel(channel_labels), size(neural, 2), neural_file);
    end

    % Remove bad channels exported as full-NaN before any group selection.
    nan_mask = squeeze(all(all(isnan(neural), 3), 1));
    neural = neural(:, ~nan_mask, :);
    channel_labels = channel_labels(~nan_mask);
    if size(neural, 2) == 0
        error('No channels remain after all-NaN removal for condition %s.', condition_str);
    end

    if doSmoothNeural
        neural = smoothdata(neural, 3, 'gaussian', smoothWindowSamples);
    end

    n_identity = size(neural, 1);
    n_time = size(neural, 3);

    %% Build condition-matched model position tensor
    % Data flow: selected trajectory source -> sorted condition trials -> position tensor.
    [model_trials, model_source_used, model_center_shift] = local_get_model_trials_for_condition( ...
        condition_str, model_path_source, Obs, Pred, obs_center_shift, pred_center_shift);
    [model_trials_sorted, sequence_values] = local_sort_and_validate_sequences(model_trials, n_identity);

    model_position = zeros(n_identity, 2, n_time);
    for i_id = 1:n_identity
        xy = double(model_trials_sorted(i_id).xy);
        if size(xy, 2) ~= 2
            error('Condition %s sequence %d is not one-dot (xy has %d columns).', ...
                condition_str, i_id, size(xy, 2));
        end
        if size(xy, 1) ~= n_time
            error('Model/neural time mismatch in %s sequence %d: model=%d neural=%d.', ...
                condition_str, i_id, size(xy, 1), n_time);
        end
        model_position(i_id, 1, :) = reshape(xy(:, 1), 1, 1, []);
        model_position(i_id, 2, :) = reshape(xy(:, 2), 1, 1, []);
    end
    model_position = bsxfun(@minus, model_position, model_center_shift);
    model_direction = local_positions_to_direction(model_position);

    %% Build concatenated model streams and trial-locked subsamples
    % Data flow: trial tensors -> feature x time streams + trigger subsamples.
    model_pos_concat = dRSA_concatenate(model_position, [], 0, 'suppressDispText', 1);
    model_dir_concat = dRSA_concatenate(model_direction, [], 0, 'suppressDispText', 1);
    [subsamples, trigger_opt] = local_build_trial_subsamples(n_time, size(model_pos_concat, 2));

    %% Model-vs-model simulations
    % Data flow: model streams -> dRSA matrices -> condition-specific outputs.
    model_out_root = fullfile(base_out_dir, 'simulations', 'model_vs_model', ...
        condition_str, lower(dRSAtypeToRun));
    model_matrices_dir = fullfile(model_out_root, 'matrices');
    model_results_dir = fullfile(model_out_root, 'results');
    if ~isfolder(model_matrices_dir), mkdir(model_matrices_dir); end
    if ~isfolder(model_results_dir), mkdir(model_results_dir); end

    model_params_self = local_build_drsa_params(lower(dRSAtypeToRun), plotSampleRateHz, n_time, pcrStrategy);
    model_params_self.modelNames = {'position model'};
    model_params_self.modelToTest = 1;
    model_params_self.modeltoRegressout = {[]};
    if strcmpi(model_params_self.dRSAtype, 'corr')
        % Char is required here so dRSA_coreFunction autocorr shortcut can call pdist directly.
        model_params_self.modelDistMeasure = 'euclidean';
    else
        model_params_self.modelDistMeasure = {'euclidean'};
    end
    model_params_self.neuralDistMeasure = 'euclidean';

    model_border_self = local_compute_autocorr_border(model_params_self, {model_pos_concat}, subsamples);
    dRSA_model_pos_pos = dRSA_coreFunction(model_pos_concat, {model_pos_concat}, model_params_self, ...
        'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', model_border_self, 'suppressDispText', 1);

    model_params_cross = local_build_drsa_params(lower(dRSAtypeToRun), plotSampleRateHz, n_time, pcrStrategy);
    model_params_cross.modelNames = {'direction model'};
    model_params_cross.modelToTest = 1;
    model_params_cross.modeltoRegressout = {[]};
    model_params_cross.modelDistMeasure = {'cosine'};
    model_params_cross.neuralDistMeasure = 'euclidean';

    model_border_cross = local_compute_autocorr_border(model_params_cross, {model_dir_concat}, subsamples);
    dRSA_model_pos_dir = dRSA_coreFunction(model_pos_concat, {model_dir_concat}, model_params_cross, ...
        'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', model_border_cross, 'suppressDispText', 1);

    model_base_tag = sprintf('dRSA_simulations_%s_%s', condition_str, lower(dRSAtypeToRun));
    model_fig_pos_pos = fullfile(model_matrices_dir, sprintf('%s_position_with_position.png', model_base_tag));
    model_fig_pos_dir = fullfile(model_matrices_dir, sprintf('%s_position_with_direction.png', model_base_tag));
    local_export_matrix_figure(dRSA_model_pos_pos(:, :, 1), plotSampleRateHz, ...
        sprintf('Model simulation: position with position (%s)', condition_str), model_fig_pos_pos);
    local_export_matrix_figure(dRSA_model_pos_dir(:, :, 1), plotSampleRateHz, ...
        sprintf('Model simulation: position with direction (%s)', condition_str), model_fig_pos_dir);

    model_res_file = fullfile(model_results_dir, sprintf('%s_model_vs_model.mat', model_base_tag));
    save(model_res_file, ...
        'dRSA_model_pos_pos', 'dRSA_model_pos_dir', ...
        'model_params_self', 'model_params_cross', ...
        'condition_str', 'subj_str', 'ses_str', ...
        'model_path_source', 'model_source_used', 'sequence_values', ...
        'obs_file', 'pred_file', 'neural_file', ...
        'model_fig_pos_pos', 'model_fig_pos_dir', ...
        'trigger_opt');

    %% Neural-vs-neural autocorrelation per sensor group
    % Data flow: cleaned neural tensor -> sensor-group slices -> per-group dRSA.
    neural_out_root = fullfile(base_out_dir, 'simulations', 'neural_vs_neural', condition_str, 'corr');
    neural_matrices_dir = fullfile(neural_out_root, 'matrices');
    neural_results_dir = fullfile(neural_out_root, 'results');
    if ~isfolder(neural_matrices_dir), mkdir(neural_matrices_dir); end
    if ~isfolder(neural_results_dir), mkdir(neural_results_dir); end

    for i_group = 1:numel(sensor_group_keys)
        sensor_group_key = sensor_group_keys{i_group};
        sensor_group_tag = sensor_group_tags{i_group};

        target_labels = string(sensor_group_lists.(sensor_group_key)(:));
        selected_mask = ismember(channel_labels, target_labels);
        selected_labels = channel_labels(selected_mask);
        if ~any(selected_mask)
            error('No channels selected for group %s in condition %s.', sensor_group_tag, condition_str);
        end

        neural_group = neural(:, selected_mask, :);
        neural_concat = dRSA_concatenate(neural_group, [], 0, 'suppressDispText', 1);

        neural_params = local_build_drsa_params('corr', plotSampleRateHz, n_time, pcrStrategy);
        neural_params.modelNames = {sprintf('neural %s', sensor_group_tag)};
        neural_params.modelToTest = 1;
        neural_params.modeltoRegressout = {[]};
        neural_params.modelDistMeasure = 'correlation';
        neural_params.neuralDistMeasure = 'correlation';

        dRSA_neural_self = dRSA_coreFunction(neural_concat, {neural_concat}, neural_params, ...
            'CurrSubsamples', subsamples(:, :, 1), 'Autocorrborder', [], 'suppressDispText', 1);

        neural_base_tag = sprintf('dRSA_simulations_%s_corr_neural_with_neural_%s', ...
            condition_str, sensor_group_tag);
        neural_fig_file = fullfile(neural_matrices_dir, sprintf('%s.png', neural_base_tag));
        local_export_matrix_figure(dRSA_neural_self(:, :, 1), plotSampleRateHz, ...
            sprintf('Neural autocorr: %s with %s (%s)', sensor_group_tag, sensor_group_tag, condition_str), ...
            neural_fig_file);

        neural_res_file = fullfile(neural_results_dir, sprintf('%s.mat', neural_base_tag));
        save(neural_res_file, ...
            'dRSA_neural_self', 'neural_params', ...
            'subj_str', 'ses_str', 'condition_str', ...
            'sensor_group_key', 'sensor_group_tag', ...
            'selected_labels', 'neural_file', 'neural_fig_file', 'trigger_opt');

        fprintf(['Saved neural autocorrelation (%s, %s):\n  %s\n  %s\n'], ...
            condition_str, sensor_group_tag, neural_fig_file, neural_res_file);
    end

    fprintf(['Completed condition %s using %s model paths. ' ...
        'Model results: %s\n'], condition_str, model_source_used, model_res_file);
end

fprintf('All requested dRSA simulations completed.\n');

%% Local helpers
function [trials, model_source_used, center_shift] = local_get_model_trials_for_condition( ...
    condition_str, model_path_source, Obs, Pred, obs_center_shift, pred_center_shift)
% Select observed/predicted trajectory trials for a condition.
% Predicted source is only valid for occluded_deviant in this dataset.
    if strcmp(model_path_source, 'predicted') && strcmp(condition_str, 'occluded_deviant')
        pred_trials = Pred.xySeqsPredicted(:);
        pred_labels = lower(string({pred_trials.condition_label}));
        mask = pred_labels == "occluded_deviant_predicted" | pred_labels == "occluded_deviant";
        trials = pred_trials(mask);
        model_source_used = 'predicted';
        center_shift = pred_center_shift;
        return;
    end

    obs_trials = Obs.xySeqs(:);
    obs_labels = lower(string({obs_trials.condition_label}));
    trials = obs_trials(obs_labels == condition_str);
    model_source_used = 'observed';
    center_shift = obs_center_shift;
end

function [sorted_trials, seq_values_sorted] = local_sort_and_validate_sequences(trials, n_identity)
% Enforce one trial per identity and sorted sequence ordering.
    if isempty(trials)
        error('No model trials found for requested condition.');
    end
    if ~isfield(trials, 'sequence') || ~isfield(trials, 'xy')
        error('Model trials must contain sequence and xy fields.');
    end

    seq_values = round(double([trials.sequence]));
    [seq_values_sorted, order] = sort(seq_values);
    expected = 1:n_identity;
    if numel(seq_values_sorted) ~= n_identity || ~isequal(seq_values_sorted(:)', expected)
        error('Model sequence values must be exactly 1:%d. Found: %s', ...
            n_identity, mat2str(seq_values_sorted));
    end
    sorted_trials = trials(order);
end

function direction_tensor = local_positions_to_direction(position_tensor)
% Convert position [identity x 2 x time] to direction [identity x 2 x time]
% using framewise displacement encoded as cos(theta), sin(theta).
    dx = diff(position_tensor(:, 1, :), 1, 3);
    dy = diff(position_tensor(:, 2, :), 1, 3);
    dx = cat(3, dx(:, :, 1), dx);
    dy = cat(3, dy(:, :, 1), dy);
    direction_tensor = cat(2, cos(atan2(dy, dx)), sin(atan2(dy, dx)));
end

function [subsamples, opt] = local_build_trial_subsamples(trial_len, total_time)
% Build triggered subsamples that align with trial starts.
    trial_starts = 1:trial_len:total_time;
    if trial_starts(end) + trial_len - 1 > total_time
        error('Trial concatenation mismatch: trial_len=%d total_time=%d.', trial_len, total_time);
    end
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
end

function params = local_build_drsa_params(drsa_type, fs_hz, trial_len, pcr_strategy)
% Build shared dRSA parameter struct with explicit defaults.
    params = struct();
    params.dRSAtype = drsa_type;
    params.nIter = 1;
    params.fs = fs_hz;
    params.AverageTime = floor((trial_len - 1) / 2) / params.fs;
    params.Var = 0.1;
    params.PCR = struct();
    params.PCR.RegressStrategy = pcr_strategy;
    params.PCR.RegressAutocor = 1;
    params.PCR.RessModel = 1;
    params.PCR.AdditionalPCA = 1;
end

function border = local_compute_autocorr_border(params, models_cell, subsamples)
% Compute PCR autocorr border only for baseline PCR strategy.
    border = [];
    if strcmpi(params.dRSAtype, 'PCR')
        if strcmp(params.PCR.RegressStrategy, 'baseline_pcr_border')
            border = dRSA_border(models_cell, subsamples, params, ...
                'suppressDispText', 1, 'plotAutocorr', 'on');
        else
            border = nan(numel(models_cell), 1);
        end
    end
end

function local_export_matrix_figure(matrix_2d, fs_hz, title_text, out_file)
% Export one dRSA matrix figure with shared axis/time formatting.
    t = (0:size(matrix_2d, 1)-1) ./ fs_hz;
    fig = figure('Visible', 'off', 'Color', 'w', 'Name', title_text, 'NumberTitle', 'off');
    imagesc(t, t, matrix_2d);
    set(gca, 'YDir', 'normal');
    axis image;
    colorbar;
    xlabel('neural time (s)');
    ylabel('model time (s)');
    title(title_text, 'Interpreter', 'none');
    hold on;
    plot(t, t, 'w-', 'LineWidth', 1);
    hold off;
    exportgraphics(fig, out_file, 'Resolution', 300, 'BackgroundColor', 'white');
    close(fig);
end

function [keys, tags] = local_normalize_sensor_groups(sensor_groups_requested)
% Normalize sensor-group names and keep stable lowercase tags for outputs.
    keys = cell(size(sensor_groups_requested));
    tags = cell(size(sensor_groups_requested));
    for i = 1:numel(sensor_groups_requested)
        raw = lower(strtrim(char(sensor_groups_requested{i})));
        if strcmp(raw, 'wholebrain') || strcmp(raw, 'whole_brain')
            key = 'wholebrain';
        elseif strcmp(raw, 'occipital')
            key = 'occipital';
        elseif strcmp(raw, 'parietal') || strcmp(raw, 'pariatel')
            key = 'parietal';
        elseif strcmp(raw, 'frontal')
            key = 'frontal';
        elseif strcmp(raw, 'temporal')
            key = 'temporal';
        else
            error('Unknown sensor group "%s".', sensor_groups_requested{i});
        end
        keys{i} = key;
        tags{i} = key;
    end
end

function sensor_groups = local_get_megin_sensor_groups()
% Hard-coded MEGIN/Neuromag Vectorview lobe selections.
% These labels match the naming used in exported channel_labels.
    sensor_groups = struct();
    sensor_groups.frontal = {'MEG0121', 'MEG0122', 'MEG0123', 'MEG0311', 'MEG0312', 'MEG0313', 'MEG0321', 'MEG0322', 'MEG0323', 'MEG0331', 'MEG0332', 'MEG0333', 'MEG0341', 'MEG0342', 'MEG0343', 'MEG0511', 'MEG0512', 'MEG0513', 'MEG0521', 'MEG0522', 'MEG0523', 'MEG0531', 'MEG0532', 'MEG0533', 'MEG0541', 'MEG0542', 'MEG0543', 'MEG0611', 'MEG0612', 'MEG0613', 'MEG0621', 'MEG0622', 'MEG0623', 'MEG0641', 'MEG0642', 'MEG0643', 'MEG0811', 'MEG0812', 'MEG0813', 'MEG0821', 'MEG0822', 'MEG0823', 'MEG0911', 'MEG0912', 'MEG0913', 'MEG0921', 'MEG0922', 'MEG0923', 'MEG0931', 'MEG0932', 'MEG0933', 'MEG0941', 'MEG0942', 'MEG0943', 'MEG1011', 'MEG1012', 'MEG1013', 'MEG1021', 'MEG1022', 'MEG1023', 'MEG1031', 'MEG1032', 'MEG1033', 'MEG1211', 'MEG1212', 'MEG1213', 'MEG1221', 'MEG1222', 'MEG1223', 'MEG1231', 'MEG1232', 'MEG1233', 'MEG1241', 'MEG1242', 'MEG1243', 'MEG1411', 'MEG1412', 'MEG1413'};
    sensor_groups.temporal = {'MEG0111', 'MEG0112', 'MEG0113', 'MEG0131', 'MEG0132', 'MEG0133', 'MEG0141', 'MEG0142', 'MEG0143', 'MEG0211', 'MEG0212', 'MEG0213', 'MEG0221', 'MEG0222', 'MEG0223', 'MEG0231', 'MEG0232', 'MEG0233', 'MEG0241', 'MEG0242', 'MEG0243', 'MEG1311', 'MEG1312', 'MEG1313', 'MEG1321', 'MEG1322', 'MEG1323', 'MEG1331', 'MEG1332', 'MEG1333', 'MEG1341', 'MEG1342', 'MEG1343', 'MEG1421', 'MEG1422', 'MEG1423', 'MEG1431', 'MEG1432', 'MEG1433', 'MEG1441', 'MEG1442', 'MEG1443', 'MEG1511', 'MEG1512', 'MEG1513', 'MEG1521', 'MEG1522', 'MEG1523', 'MEG1531', 'MEG1532', 'MEG1533', 'MEG1541', 'MEG1542', 'MEG1543', 'MEG1611', 'MEG1612', 'MEG1613', 'MEG1621', 'MEG1622', 'MEG1623', 'MEG2411', 'MEG2412', 'MEG2413', 'MEG2421', 'MEG2422', 'MEG2423', 'MEG2611', 'MEG2612', 'MEG2613', 'MEG2621', 'MEG2622', 'MEG2623', 'MEG2631', 'MEG2632', 'MEG2633', 'MEG2641', 'MEG2642', 'MEG2643'};
    sensor_groups.parietal = {'MEG0411', 'MEG0412', 'MEG0413', 'MEG0421', 'MEG0422', 'MEG0423', 'MEG0431', 'MEG0432', 'MEG0433', 'MEG0441', 'MEG0442', 'MEG0443', 'MEG0631', 'MEG0632', 'MEG0633', 'MEG0711', 'MEG0712', 'MEG0713', 'MEG0721', 'MEG0722', 'MEG0723', 'MEG0731', 'MEG0732', 'MEG0733', 'MEG0741', 'MEG0742', 'MEG0743', 'MEG1041', 'MEG1042', 'MEG1043', 'MEG1111', 'MEG1112', 'MEG1113', 'MEG1121', 'MEG1122', 'MEG1123', 'MEG1131', 'MEG1132', 'MEG1133', 'MEG1141', 'MEG1142', 'MEG1143', 'MEG1631', 'MEG1632', 'MEG1633', 'MEG1811', 'MEG1812', 'MEG1813', 'MEG1821', 'MEG1822', 'MEG1823', 'MEG1831', 'MEG1832', 'MEG1833', 'MEG1841', 'MEG1842', 'MEG1843', 'MEG2011', 'MEG2012', 'MEG2013', 'MEG2021', 'MEG2022', 'MEG2023', 'MEG2211', 'MEG2212', 'MEG2213', 'MEG2221', 'MEG2222', 'MEG2223', 'MEG2231', 'MEG2232', 'MEG2233', 'MEG2241', 'MEG2242', 'MEG2243', 'MEG2441', 'MEG2442', 'MEG2443'};
    sensor_groups.occipital = {'MEG1641', 'MEG1642', 'MEG1643', 'MEG1711', 'MEG1712', 'MEG1713', 'MEG1721', 'MEG1722', 'MEG1723', 'MEG1731', 'MEG1732', 'MEG1733', 'MEG1741', 'MEG1742', 'MEG1743', 'MEG1911', 'MEG1912', 'MEG1913', 'MEG1921', 'MEG1922', 'MEG1923', 'MEG1931', 'MEG1932', 'MEG1933', 'MEG1941', 'MEG1942', 'MEG1943', 'MEG2031', 'MEG2032', 'MEG2033', 'MEG2041', 'MEG2042', 'MEG2043', 'MEG2111', 'MEG2112', 'MEG2113', 'MEG2121', 'MEG2122', 'MEG2123', 'MEG2131', 'MEG2132', 'MEG2133', 'MEG2141', 'MEG2142', 'MEG2143', 'MEG2311', 'MEG2312', 'MEG2313', 'MEG2321', 'MEG2322', 'MEG2323', 'MEG2331', 'MEG2332', 'MEG2333', 'MEG2341', 'MEG2342', 'MEG2343', 'MEG2431', 'MEG2432', 'MEG2433', 'MEG2511', 'MEG2512', 'MEG2513', 'MEG2521', 'MEG2522', 'MEG2523', 'MEG2531', 'MEG2532', 'MEG2533', 'MEG2541', 'MEG2542', 'MEG2543'};

    % Whole brain is the union of all lobe lists.
    sensor_groups.wholebrain = unique([ ...
        sensor_groups.frontal, ...
        sensor_groups.temporal, ...
        sensor_groups.parietal, ...
        sensor_groups.occipital], 'stable');
end
