% A_dRSA_alt.m
%
% Purpose:
%   Run a single real-neural dRSA analysis (neural vs one-dot position model)
%   for one condition, with a tunable MEGIN sensor group selection.
%
% Data-flow summary:
%   1) Load identity-averaged neural tensor + canonical MEG channel labels.
%   2) Load observed or predicted one-dot trajectory model for the condition.
%   3) Remove all-NaN channels from the canonical export axis.
%   4) Apply a MEGIN sensor-group filter:
%      wholeBrain | occipital | parietal | frontal | temporal.
%   5) Build model position tensor, concatenate trial-wise inputs, and run dRSA.
%   6) Save matrix PNG + MAT result with a filename tag that includes the
%      selected sensor group.
%
% Usage examples:
%   % Example 1: default run (sub-01, always_visible, corr, wholeBrain).
%   eval(fileread('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/analysis/A_dRSA_alt.m'));
%
%   % Example 2: occipital-only run.
%   subj_str = '01';
%   condition_str = 'always_visible';
%   dRSAtypeToRun = 'corr'; % or 'PCR'
%   sensor_group = 'occipital'; % wholeBrain|occipital|parietal|frontal|temporal
%   eval(fileread('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/analysis/A_dRSA_alt.m'));
%
% Inputs:
%   - Neural MAT (identity-averaged): neural [20 x nChannels x nTime]
%   - Channel labels: channel_labels [nChannels x 1] or [1 x nChannels]
%   - Observed model MAT: xySeqs + Cfg.rectSize
%   - Predicted model MAT (used for occluded_deviant + predicted path_type)
%
% Outputs:
%   - dRSA matrix figure (PNG)
%   - dRSA result MAT
%   - Filenames include sensor tag, e.g.:
%     A_dRSA_alt_always_visible_corr_occipital.png
%
% Assumptions:
%   - Model trials are already ordered by identity 1..20.
%   - Neural and model time axes are already aligned.
%   - Export step 06 used canonical channel labels, with bad channels written
%     as all-NaN (removed before dRSA).
%   - MEGIN/Vectorview lobe groups follow the standard Neuromag selections.

clearvars -except subj_str ses_str condition_str dRSAtypeToRun plotSampleRateHz pcrStrategy doSmoothNeural smoothWindowSamples path_type sensor_group output_root_override
smoothWindowSamples = 20;
%% Runtime defaults
if ~exist('subj_str', 'var') || isempty(subj_str), subj_str = '01'; end
if ~exist('ses_str', 'var') || isempty(ses_str), ses_str = '01'; end
if ~exist('condition_str', 'var') || isempty(condition_str), condition_str = 'always_visible'; end
if ~exist('path_type', 'var') || isempty(path_type), path_type = 'observed'; end
if ~exist('dRSAtypeToRun', 'var') || isempty(dRSAtypeToRun), dRSAtypeToRun = 'corr'; end
if ~exist('plotSampleRateHz', 'var') || isempty(plotSampleRateHz), plotSampleRateHz = 120; end
if ~exist('pcrStrategy', 'var') || isempty(pcrStrategy), pcrStrategy = 'baseline_pcr_border'; end %'ridge_full_autocorr'; end
if ~exist('doSmoothNeural', 'var') || isempty(doSmoothNeural), doSmoothNeural = true; end
if ~exist('smoothWindowSamples', 'var') || isempty(smoothWindowSamples), smoothWindowSamples = 15; end
if ~exist('sensor_group', 'var') || isempty(sensor_group), sensor_group = 'wholeBrain'; end

% Canonicalize region selector while keeping a stable output tag string.
sensor_group_key = lower(strtrim(char(sensor_group)));
if strcmp(sensor_group_key, 'wholebrain') || strcmp(sensor_group_key, 'whole_brain')
    sensor_group_key = 'wholebrain';
elseif strcmp(sensor_group_key, 'frontal')
    sensor_group_key = 'frontal';
elseif strcmp(sensor_group_key, 'temporal')
    sensor_group_key = 'temporal';
elseif strcmp(sensor_group_key, 'parietal') || strcmp(sensor_group_key, 'pariatel')
    sensor_group_key = 'parietal';
elseif strcmp(sensor_group_key, 'occipital')
    sensor_group_key = 'occipital';
else
    error(['Invalid sensor_group: %s. Valid options are wholeBrain, occipital, ' ...
        'parietal, frontal, temporal.'], sensor_group);
end
sensor_group_tag = sensor_group_key;

%% Paths
dad_root = '/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD';
meg_root = fullfile(dad_root, 'MEG');
drsa_dir = fullfile(meg_root, 'derivatives', 'analysis', ...
    sprintf('sub-%s', subj_str), sprintf('ses-%s', ses_str), 'drsa');
model_dir = fullfile(dad_root, 'experiment', 'input_files', sprintf('sub%s', subj_str));

neural_file = fullfile(drsa_dir, sprintf( ...
    'sub-%s_ses-%s_task-Attention_condition-%s_desc-identityAveraged_neural.mat', ...
    subj_str, ses_str, condition_str));
obs_file = fullfile(model_dir, sprintf('MovDot_Sub%s_V28_rescueTraject.mat', subj_str));
pred_file = fullfile(model_dir, sprintf('MovDot_Sub%s_V28_rescueTraject_predicted.mat', subj_str));

%% Load neural
N = load(neural_file, 'neural', 'channel_labels');
neural = N.neural;
channel_labels = string(N.channel_labels(:));

%% if occluded deviant and var path_type is 'predicted', position model will be from predicted not observed
use_predicted_model = strcmpi(condition_str, 'occluded_deviant') && strcmpi(path_type, 'predicted');
if use_predicted_model
    M = load(pred_file, 'xySeqsPredicted', 'Cfg');
    trials = M.xySeqsPredicted( ...
        strcmpi({M.xySeqsPredicted.condition_label}, 'occluded_deviant_predicted') | ...
        strcmpi({M.xySeqsPredicted.condition_label}, 'occluded_deviant'));
    model_source = 'predicted';
else
    M = load(obs_file, 'xySeqs', 'Cfg');
    trials = M.xySeqs(strcmpi({M.xySeqs.condition_label}, condition_str));
    model_source = 'observed';
end

%% Remove all-NaN channels from canonical export axis
if numel(channel_labels) ~= size(neural, 2)
    error('channel_labels count (%d) does not match neural channels (%d).', ...
        numel(channel_labels), size(neural, 2));
end
nan_mask = squeeze(all(all(isnan(neural), 3), 1));
removed_labels = cellstr(channel_labels(nan_mask));
if any(nan_mask)
    fprintf('All-NaN channels removed before A_dRSA_alt: %d\n', numel(removed_labels));
    fprintf('Removed channel labels:\n');
    for i_ch = 1:numel(removed_labels)
        fprintf('  %s\n', removed_labels{i_ch});
    end
else
    fprintf('All-NaN channels removed before A_dRSA_alt: 0\n');
end
neural = neural(:, ~nan_mask, :);
channel_labels = channel_labels(~nan_mask);
if size(neural, 2) == 0
    error('No channels remain after removing all-NaN channels in %s.', neural_file);
end

%% Apply MEGIN sensor-group selection after NaN cleanup
%
% These channel lists are hard-coded from the standard Neuromag/MEGIN
% Vectorview lobe selections (same family as your attached layout).
sensor_groups = struct();
sensor_groups.wholebrain = cellstr(channel_labels);
sensor_groups.frontal = {'MEG0121', 'MEG0122', 'MEG0123', 'MEG0311', 'MEG0312', 'MEG0313', 'MEG0321', 'MEG0322', 'MEG0323', 'MEG0331', 'MEG0332', 'MEG0333', 'MEG0341', 'MEG0342', 'MEG0343', 'MEG0511', 'MEG0512', 'MEG0513', 'MEG0521', 'MEG0522', 'MEG0523', 'MEG0531', 'MEG0532', 'MEG0533', 'MEG0541', 'MEG0542', 'MEG0543', 'MEG0611', 'MEG0612', 'MEG0613', 'MEG0621', 'MEG0622', 'MEG0623', 'MEG0641', 'MEG0642', 'MEG0643', 'MEG0811', 'MEG0812', 'MEG0813', 'MEG0821', 'MEG0822', 'MEG0823', 'MEG0911', 'MEG0912', 'MEG0913', 'MEG0921', 'MEG0922', 'MEG0923', 'MEG0931', 'MEG0932', 'MEG0933', 'MEG0941', 'MEG0942', 'MEG0943', 'MEG1011', 'MEG1012', 'MEG1013', 'MEG1021', 'MEG1022', 'MEG1023', 'MEG1031', 'MEG1032', 'MEG1033', 'MEG1211', 'MEG1212', 'MEG1213', 'MEG1221', 'MEG1222', 'MEG1223', 'MEG1231', 'MEG1232', 'MEG1233', 'MEG1241', 'MEG1242', 'MEG1243', 'MEG1411', 'MEG1412', 'MEG1413'};
sensor_groups.temporal = {'MEG0111', 'MEG0112', 'MEG0113', 'MEG0131', 'MEG0132', 'MEG0133', 'MEG0141', 'MEG0142', 'MEG0143', 'MEG0211', 'MEG0212', 'MEG0213', 'MEG0221', 'MEG0222', 'MEG0223', 'MEG0231', 'MEG0232', 'MEG0233', 'MEG0241', 'MEG0242', 'MEG0243', 'MEG1311', 'MEG1312', 'MEG1313', 'MEG1321', 'MEG1322', 'MEG1323', 'MEG1331', 'MEG1332', 'MEG1333', 'MEG1341', 'MEG1342', 'MEG1343', 'MEG1421', 'MEG1422', 'MEG1423', 'MEG1431', 'MEG1432', 'MEG1433', 'MEG1441', 'MEG1442', 'MEG1443', 'MEG1511', 'MEG1512', 'MEG1513', 'MEG1521', 'MEG1522', 'MEG1523', 'MEG1531', 'MEG1532', 'MEG1533', 'MEG1541', 'MEG1542', 'MEG1543', 'MEG1611', 'MEG1612', 'MEG1613', 'MEG1621', 'MEG1622', 'MEG1623', 'MEG2411', 'MEG2412', 'MEG2413', 'MEG2421', 'MEG2422', 'MEG2423', 'MEG2611', 'MEG2612', 'MEG2613', 'MEG2621', 'MEG2622', 'MEG2623', 'MEG2631', 'MEG2632', 'MEG2633', 'MEG2641', 'MEG2642', 'MEG2643'};
sensor_groups.parietal = {'MEG0411', 'MEG0412', 'MEG0413', 'MEG0421', 'MEG0422', 'MEG0423', 'MEG0431', 'MEG0432', 'MEG0433', 'MEG0441', 'MEG0442', 'MEG0443', 'MEG0631', 'MEG0632', 'MEG0633', 'MEG0711', 'MEG0712', 'MEG0713', 'MEG0721', 'MEG0722', 'MEG0723', 'MEG0731', 'MEG0732', 'MEG0733', 'MEG0741', 'MEG0742', 'MEG0743', 'MEG1041', 'MEG1042', 'MEG1043', 'MEG1111', 'MEG1112', 'MEG1113', 'MEG1121', 'MEG1122', 'MEG1123', 'MEG1131', 'MEG1132', 'MEG1133', 'MEG1141', 'MEG1142', 'MEG1143', 'MEG1631', 'MEG1632', 'MEG1633', 'MEG1811', 'MEG1812', 'MEG1813', 'MEG1821', 'MEG1822', 'MEG1823', 'MEG1831', 'MEG1832', 'MEG1833', 'MEG1841', 'MEG1842', 'MEG1843', 'MEG2011', 'MEG2012', 'MEG2013', 'MEG2021', 'MEG2022', 'MEG2023', 'MEG2211', 'MEG2212', 'MEG2213', 'MEG2221', 'MEG2222', 'MEG2223', 'MEG2231', 'MEG2232', 'MEG2233', 'MEG2241', 'MEG2242', 'MEG2243', 'MEG2441', 'MEG2442', 'MEG2443'};
sensor_groups.occipital = {'MEG1641', 'MEG1642', 'MEG1643', 'MEG1711', 'MEG1712', 'MEG1713', 'MEG1721', 'MEG1722', 'MEG1723', 'MEG1731', 'MEG1732', 'MEG1733', 'MEG1741', 'MEG1742', 'MEG1743', 'MEG1911', 'MEG1912', 'MEG1913', 'MEG1921', 'MEG1922', 'MEG1923', 'MEG1931', 'MEG1932', 'MEG1933', 'MEG1941', 'MEG1942', 'MEG1943', 'MEG2031', 'MEG2032', 'MEG2033', 'MEG2041', 'MEG2042', 'MEG2043', 'MEG2111', 'MEG2112', 'MEG2113', 'MEG2121', 'MEG2122', 'MEG2123', 'MEG2131', 'MEG2132', 'MEG2133', 'MEG2141', 'MEG2142', 'MEG2143', 'MEG2311', 'MEG2312', 'MEG2313', 'MEG2321', 'MEG2322', 'MEG2323', 'MEG2331', 'MEG2332', 'MEG2333', 'MEG2341', 'MEG2342', 'MEG2343', 'MEG2431', 'MEG2432', 'MEG2433', 'MEG2511', 'MEG2512', 'MEG2513', 'MEG2521', 'MEG2522', 'MEG2523', 'MEG2531', 'MEG2532', 'MEG2533', 'MEG2541', 'MEG2542', 'MEG2543'};

target_labels = string(sensor_groups.(sensor_group_key)(:));
selected_mask = ismember(channel_labels, target_labels);
selected_channel_labels = channel_labels(selected_mask);
missing_group_labels = setdiff(target_labels, channel_labels, 'stable');

fprintf('Sensor group requested: %s\n', sensor_group_tag);
fprintf('Channels retained after group filter: %d\n', nnz(selected_mask));
if ~isempty(missing_group_labels)
    fprintf('Group labels missing after NaN-removal / export differences: %d\n', numel(missing_group_labels));
end

neural = neural(:, selected_mask, :);
channel_labels = selected_channel_labels;
if size(neural, 2) == 0
    error('No channels remain after applying sensor_group=%s.', sensor_group_tag);
end

%% Build position tensor [identity x 2 x time] from selected model source
n_id = size(neural, 1);
n_time = size(neural, 3);
model_center_shift = reshape(double(M.Cfg.rectSize(:) ./ 2), [1, 2, 1]);
model_position = zeros(n_id, 2, n_time);
for i_id = 1:n_id
    xy = double(trials(i_id).xy);
    model_position(i_id, 1, :) = reshape(xy(:, 1), 1, 1, []);
    model_position(i_id, 2, :) = reshape(xy(:, 2), 1, 1, []);
end
model_position = bsxfun(@minus, model_position, model_center_shift);

%% Add dRSA helpers
addpath(fullfile(dad_root, 'simulations', 'functions'));

%% let's try smoothing
if doSmoothNeural
    neural = smoothdata(neural, 3, 'gaussian', smoothWindowSamples);
end

%% Prepare concatenated inputs and triggered subsampling
neural_concat = dRSA_concatenate(neural, [], 0, 'suppressDispText', 1);
model_pos = dRSA_concatenate(model_position, [], 0, 'suppressDispText', 1);

trial_len = size(neural, 3);
total_time = size(neural_concat, 2);
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

%% dRSA parameters (single model: observed position)
params = struct();
params.dRSAtype = dRSAtypeToRun;
params.modelNames = {'position model'};
params.modelToTest = 1;
params.modelDistMeasure = {'euclidean'};
params.neuralDistMeasure = 'correlation';
params.modeltoRegressout = {[]};
params.nIter = 1;
params.fs = plotSampleRateHz;
params.AverageTime = floor((trial_len - 1) / 2) / params.fs;
params.Var = 0.1;
params.PCR = struct();
params.PCR.RegressStrategy = pcrStrategy;
params.PCR.RegressAutocor = 1;
params.PCR.RessModel = 1;
params.PCR.AdditionalPCA = 1;

autocorr_border = [];
if strcmpi(params.dRSAtype, 'PCR')
    autocorr_border = dRSA_border({model_pos}, subsamples, params, 'suppressDispText', 1, 'plotAutocorr', 'on');
end

%% Run one dRSA matrix
dRSA_neural_position = dRSA_coreFunction(neural_concat, {model_pos}, params, ...
    'CurrSubsamples', subsamples(:, :, 1), ...
    'Autocorrborder', autocorr_border, ...
    'suppressDispText', 1);

%% Save outputs
% Allow sandbox validation runs to redirect outputs without touching project
% derivatives. When unset, behavior is unchanged (writes under drsa_dir).
if exist('output_root_override', 'var') && ~isempty(output_root_override)
    base_out_dir = char(output_root_override);
else
    base_out_dir = drsa_dir;
end

if strcmpi(condition_str, 'occluded_deviant')
    source_tag = lower(model_source);
    out_root = fullfile(base_out_dir, 'matlab_real_neural', 'A_dRSA_alt', ...
        condition_str, source_tag, lower(dRSAtypeToRun));
    file_tag = sprintf('A_dRSA_alt_%s_%s_%s_%s', condition_str, source_tag, lower(dRSAtypeToRun), sensor_group_tag);
else
    out_root = fullfile(base_out_dir, 'matlab_real_neural', 'A_dRSA_alt', ...
        condition_str, lower(dRSAtypeToRun));
    file_tag = sprintf('A_dRSA_alt_%s_%s_%s', condition_str, lower(dRSAtypeToRun), sensor_group_tag);
end
matrices_dir = fullfile(out_root, 'matrices'); if ~isfolder(matrices_dir), mkdir(matrices_dir); end
results_dir = fullfile(out_root, 'results'); if ~isfolder(results_dir), mkdir(results_dir); end

t = (0:size(dRSA_neural_position, 1)-1) ./ plotSampleRateHz;
fig = figure('Name', sprintf('A_dRSA_alt (%s, %s)', condition_str, dRSAtypeToRun), 'NumberTitle', 'off', 'Visible', 'on');
imagesc(t, t, dRSA_neural_position(:, :, 1));
set(gca, 'YDir', 'normal');
axis image;
%caxis([0 1]);
colorbar;
xlabel(sprintf('neural time (s) \n RDM measure %s', params.neuralDistMeasure));
ylabel(sprintf('model time (s) \n RDM measure %s', params.modelDistMeasure{1}));
title(sprintf(' %s \n %s \n %s (%s) - Smoothing %d ', strrep(condition_str, "_", " "), [model_source ' ' params.modelNames{1}], dRSAtypeToRun, sensor_group_tag, smoothWindowSamples), 'Interpreter', 'none');
hold on;
plot(t, t, 'w-', 'LineWidth', 1);
if contains(lower(condition_str), 'occluded')
    y_occ_complete = (double(trials(1).occlusion_complete_frame) - 1) ./ plotSampleRateHz;
    y_reapp_start = (double(trials(1).reappearance_start_frame) - 1) ./ plotSampleRateHz;
    yline(y_occ_complete, '--', 'LineWidth', 2);
    yline(y_reapp_start, '--', 'LineWidth', 2);
end
hold off;

fig_file = fullfile(matrices_dir, sprintf('%s.png', file_tag));
exportgraphics(fig, fig_file, 'Resolution', 300, 'BackgroundColor', 'white');

res_file = fullfile(results_dir, sprintf('%s.mat', file_tag));
save(res_file, 'dRSA_neural_position', 'params', 'subj_str', 'ses_str', 'condition_str', ...
    'path_type', 'model_source', 'sensor_group_tag', 'channel_labels', ...
    'neural_file', 'obs_file', 'pred_file', 'fig_file');

fprintf('Saved:\n  %s\n  %s\n', fig_file, res_file);
