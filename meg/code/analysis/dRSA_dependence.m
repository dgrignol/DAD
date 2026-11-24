function dRSA_dependence(handoff_path, output_dir, conditions, sensor_groups, smooth_neural, smooth_window, model_sources)
% DRSA_DEPENDENCE Write descriptive model and MEG RDM correlation matrices.
%
% Usage:
%   addpath('meg/code/analysis');
%   dRSA_dependence('/tmp/dad/07_handoff.json', '/tmp/dad/staging', ...
%       {'occluded_deviant'}, {'occipital'}, true, 15);
%   In this integration, run_drsa_dependence_cohort.py supplies these
%   arguments from a validated Stage-06 handoff.
%
% Inputs: one subject's existing neural exports and trajectory MAT, an
% existing staging directory, condition/group lists and neural smoothing
% settings matching A_dRSA_alt, and optional model-source list. Output: one v7 MAT per condition with raw
% time-by-time Pearson correlations, indexed [row time, column time], and
% associated identity, time, channel and trajectory metadata.

    % Legacy callers keep the original observed/predicted matrix inventory.
    if nargin < 7, model_sources = {'observed', 'predicted'}; end
    handoff = drsa_stage_utils.load_handoff(handoff_path);
    drsa_stage_utils.add_helper_paths(handoff);
    identities = drsa_stage_utils.vectorize_numeric(handoff.identity_values);
    time_axis_sec = drsa_stage_utils.vectorize_numeric(handoff.time_axis_sec);
    n_time = numel(time_axis_sec);
    params = drsa_stage_utils.build_drsa_params('corr', handoff.sampling_rate_hz, n_time, '');
    [subsamples, ~] = drsa_stage_utils.build_trial_subsamples(n_time, numel(identities) * n_time);
    subsamples = subsamples(:, :, 1);
    observed = drsa_stage_utils.load_models(handoff);
    megin_groups = drsa_stage_utils.get_megin_sensor_groups();

    for i = 1:numel(conditions)
        condition = conditions{i};
        fprintf('RDM dependence | %s | %s\n', handoff.subject_bids_id, condition);

        % Build each model once with the exact A_dRSA_alt feature construction.
        % Cross-model panels use matched condensed identity-pair RDM entries;
        % MATLAB corr stores first-argument time on rows, second on columns.
        sources = {'observed'};
        if strcmp(condition, 'occluded_deviant'), sources = model_sources; end
        rdms = struct(); model_metadata = struct(); source_conditions = struct();
        matrices = struct(); shared_prefix_frames = [];
        for m = 1:numel(sources)
            source = sources{m};
            [tensor, ~, source_condition, feature] = drsa_build_alt_model( ...
                condition, source, observed, identities, n_time);
            rdms.(source) = tensor_rdm(tensor, params, subsamples, 'euclidean');
            matrices.([source '_self']) = corr(rdms.(source), rdms.(source));
            source_conditions.(source) = source_condition;
            model_metadata.(source) = struct('feature', feature, ...
                'feature_count', size(tensor, 2), 'rdm_distance', 'euclidean', ...
                'label', drsa_stage_utils.stimulus_model_label(condition, source, feature));
        end
        if strcmp(condition, 'occluded_deviant')
            matrices.observed_vs_predicted = corr(rdms.observed, rdms.predicted);
            [trials, ~, ~, ~] = drsa_stage_utils.model_trials_for_condition(condition, 'observed', observed);
            shared_prefix_frames = min([trials.deviance_frame]);
            for m = 1:numel(sources)
                source = sources{m};
                if ismember(source, {'observed', 'predicted'}), continue; end
                matrices.(['observed_vs_' source]) = corr(rdms.observed, rdms.(source));
                matrices.(['predicted_vs_' source]) = corr(rdms.predicted, rdms.(source));
            end
        end

        % Use the same channel exclusion, temporal smoothing and metric as alt.
        [neural, labels, ~] = drsa_stage_utils.load_neural_condition(handoff, condition);
        [neural, labels, removed_labels] = drsa_stage_utils.remove_all_nan_channels(neural, labels);
        if smooth_neural
            neural = smoothdata(neural, 3, 'gaussian', smooth_window);
        end
        selected_channel_labels = struct(); normalization_info = struct();
        for g = 1:numel(sensor_groups)
            [key, tag] = drsa_stage_utils.normalize_sensor_group(sensor_groups{g});
            selected = ismember(labels, string(megin_groups.(key)(:)));
            if ~any(selected)
                error('No channels remain for sensor group %s in %s.', tag, condition);
            end
            % Reuse the exact production transform after selecting this ROI.
            % Model RDMs above remain unchanged by sensor-noise normalization.
            [neural_group, retained_labels, normalization_info.(tag)] = ...
                drsa_stage_utils.apply_noise_normalization(neural(:, selected, :), labels(selected), handoff, 'alt', key);
            selected_channel_labels.(tag) = cellstr(retained_labels);
            neural_rdm = tensor_rdm(neural_group, params, subsamples, 'correlation');
            matrices.(['neural_' tag]) = corr(neural_rdm, neural_rdm);
        end

        save(fullfile(output_dir, [condition '.mat']), 'matrices', 'condition', ...
            'time_axis_sec', 'identities', 'source_conditions', 'shared_prefix_frames', ...
            'selected_channel_labels', 'removed_labels', 'model_metadata', 'normalization_info', '-v7');
    end
end

function rdm = tensor_rdm(tensor, params, subsamples, metric)
    % Flatten identity epochs just as Stage 07 does before its RDM calculation.
    concatenated = dRSA_concatenate(tensor, [], 0, 'suppressDispText', 1);
    rdm = dRSA_computeRDM(concatenated, params, subsamples, metric);
end
