function [model_tensor, model_source, model_source_condition, model_feature] = ...
        drsa_build_alt_model( ...
        condition_str, path_type, Obs, identity_values, n_time)
% DRSA_BUILD_ALT_MODEL Shared, numerically unchanged A_dRSA_alt model builder.
%
% Usage:
%   [tensor, source, condition, feature] = drsa_build_alt_model( ...
%       'occluded_deviant', 'prediction_error', Obs, 1:40, 320);
%   Called by A_dRSA_alt and dRSA_dependence after validated handoff loading.
%   Euclidean RDMs are computed by dRSA_computeRDM after concatenation.
%   No feature rescaling, masking, or additional smoothing is applied here.
%
% Inputs:
%   condition_str   : condition requested by Stage 07.
%   path_type       : observed, predicted, prediction_error, or pe_* source.
%   Obs             : observed trajectory MAT payload with xySeqs and Cfg.
%   identity_values : identity axis exported by Stage 06.
%   n_time          : neural/model sample count.
%
% Outputs:
%   model_tensor           : identity x feature x time dRSA model tensor.
%   model_source           : stable source string used in output folders.
%   model_source_condition : readable provenance for the source condition(s).
%   model_feature          : feature label passed to report label helpers.
%
% Data flow:
%   path_type -> observed/predicted trial selection -> tensor construction ->
%   optional PE-variant transform.

if local_is_pe_model_source(path_type)
    if ~strcmp(condition_str, 'occluded_deviant')
        error('%s path_type is valid only for occluded_deviant; got %s.', ...
            path_type, condition_str);
    end

    [observed_trials, ~, observed_center_shift, ~] = ...
        drsa_stage_utils.model_trials_for_condition(condition_str, 'observed', Obs);
    observed_sorted = drsa_stage_utils.sort_and_validate_sequences( ...
        observed_trials, identity_values);
    [observed_position, ~] = drsa_stage_utils.build_position_direction( ...
        observed_sorted, identity_values, n_time, observed_center_shift);

    [predicted_trials, ~, predicted_center_shift, predicted_source_condition] = ...
        drsa_stage_utils.model_trials_for_condition(condition_str, 'predicted', Obs);
    predicted_sorted = drsa_stage_utils.sort_and_validate_sequences( ...
        predicted_trials, identity_values);
    [predicted_position, ~] = drsa_stage_utils.build_position_direction( ...
        predicted_sorted, identity_values, n_time, predicted_center_shift);

    [model_tensor, model_feature] = local_build_pe_variant_tensor( ...
        path_type, observed_position, predicted_position, ...
        observed_sorted, predicted_sorted, Obs.Cfg, n_time);
    model_source = path_type;
    model_source_condition = sprintf('%s-minus-%s', ...
        condition_str, predicted_source_condition);
else
    [trials, model_source, center_shift, model_source_condition] = ...
        drsa_stage_utils.model_trials_for_condition(condition_str, path_type, Obs);
    sorted_trials = drsa_stage_utils.sort_and_validate_sequences(trials, identity_values);
    [model_tensor, ~] = drsa_stage_utils.build_position_direction( ...
        sorted_trials, identity_values, n_time, center_shift);
    model_feature = 'position';
end
end

function tf = local_is_pe_model_source(path_type)
% LOCAL_IS_PE_MODEL_SOURCE Return true for legacy and explicit PE sources.
path_type = lower(strtrim(char(path_type)));
tf = strcmp(path_type, 'prediction_error') || strncmp(path_type, 'pe_', 3);
end

function [model_tensor, model_feature] = local_build_pe_variant_tensor( ...
        path_type, observed_position, predicted_position, ...
        observed_trials, predicted_trials, Cfg, n_time)
% LOCAL_BUILD_PE_VARIANT_TENSOR Build the requested PE model tensor.
%
% The legacy prediction_error source remains the raw signed displacement. New
% variants use deterministic transforms so Stage 07 outputs are reproducible
% and do not depend on a Monte Carlo particle count.
residual = observed_position - predicted_position;

switch lower(strtrim(char(path_type)))
    case 'prediction_error'
        model_tensor = residual;
        model_feature = 'position';
    case 'pe_prior_unconstrained_50'
        model_tensor = local_pe_prior_unconstrained_50(residual, Cfg);
        model_feature = 'position';
    case 'pe_prior_constrained_50'
        model_tensor = local_pe_prior_constrained_50(residual, Cfg);
        model_feature = 'position';
    case 'pe_precision_weighted'
        model_tensor = local_pe_precision_weighted(residual, Cfg);
        model_feature = 'position';
    case 'pe_kalman_innovation'
        model_tensor = local_pe_kalman_innovation(observed_position, predicted_position);
        model_feature = 'position';
    case 'pe_latent_curvature'
        observed_curvature = local_curvature_tensor(observed_trials, n_time);
        predicted_curvature = local_curvature_tensor(predicted_trials, n_time);
        model_tensor = observed_curvature - predicted_curvature;
        model_feature = 'latent_curvature';
    case 'pe_belief_update'
        model_tensor = local_pe_belief_update(residual, Cfg);
        model_feature = 'belief_update';
    otherwise
        error('Unsupported PE path_type: %s', path_type);
end
end

function model_tensor = local_pe_prior_unconstrained_50(residual, Cfg)
% LOCAL_PE_PRIOR_UNCONSTRAINED_50 50% known path plus broad uniform branch.
%
% The known branch is a Gaussian around the expected nondeviant path. The
% unknown branch is uniform over the motion rectangle. The output keeps signed
% residual direction but gates it by the posterior responsibility of the known
% branch and appends a surprise channel.
sigma = local_spatial_sigma(Cfg);
rect_area = local_motion_rect_area(Cfg);
mag_sq = sum(residual .^ 2, 2);
known_density = exp(-0.5 .* mag_sq ./ (sigma .^ 2)) ./ (2 .* pi .* sigma .^ 2);
uniform_density = 1 ./ rect_area;
mixture_density = 0.5 .* known_density + 0.5 .* uniform_density;
responsibility = (0.5 .* known_density) ./ max(mixture_density, realmin);
weighted_residual = bsxfun(@times, residual, responsibility);
surprise = -log(max(mixture_density, realmin));
model_tensor = cat(2, weighted_residual, surprise);
end

function model_tensor = local_pe_prior_constrained_50(residual, Cfg)
% LOCAL_PE_PRIOR_CONSTRAINED_50 50% known path plus empirical deviant branch.
%
% The constrained branch is estimated from the same generated stimulus set:
% at each time point, the possible deviant displacement is approximated by the
% leave-one-out residual distribution across sequence identities. This is a
% deterministic proxy for knowing the experiment's turn/curvature constraints.
loo_mean = local_leave_one_out_mean(residual);
mixture_expected_residual = 0.5 .* loo_mean;
centered_residual = residual - mixture_expected_residual;
sigma = local_timewise_feature_scale(residual, Cfg);
standardized = centered_residual ./ sigma;
surprise = sqrt(sum(standardized .^ 2, 2));
model_tensor = cat(2, centered_residual, surprise);
end

function model_tensor = local_pe_precision_weighted(residual, Cfg)
% LOCAL_PE_PRECISION_WEIGHTED Scale signed PE by empirical predictive precision.
sigma = local_timewise_feature_scale(residual, Cfg);
model_tensor = residual ./ sigma;
end

function model_tensor = local_pe_kalman_innovation(observed_position, predicted_position)
% LOCAL_PE_KALMAN_INNOVATION Simple online innovation with fixed Kalman gain.
%
% The state follows the expected trajectory velocity and is corrected toward
% each observed sample. This produces a transient innovation that decays as the
% observer updates after reappearance.
gain = 0.35;
n_identity = size(observed_position, 1);
n_time = size(observed_position, 3);
model_tensor = zeros(size(observed_position));
for i_id = 1:n_identity
    state_xy = reshape(predicted_position(i_id, :, 1), 1, []);
    for i_time = 1:n_time
        observed_xy = reshape(observed_position(i_id, :, i_time), 1, []);
        innovation_xy = observed_xy - state_xy;
        model_tensor(i_id, :, i_time) = reshape(innovation_xy, 1, 2, 1);
        if i_time < n_time
            expected_step = reshape( ...
                predicted_position(i_id, :, i_time + 1) - predicted_position(i_id, :, i_time), ...
                1, []);
            state_xy = state_xy + expected_step + gain .* innovation_xy;
        end
    end
end
end

function model_tensor = local_pe_belief_update(residual, Cfg)
% LOCAL_PE_BELIEF_UPDATE Convert residual magnitude into latent-cause evidence.
%
% Channel 1 is an approximate posterior that the active path is deviant.
% Channel 2 is the positive frame-to-frame posterior update, which highlights
% when the residual first becomes informative.
sigma = local_spatial_sigma(Cfg);
mag = sqrt(sum(residual .^ 2, 2));
posterior_deviant = 1 - exp(-0.5 .* (mag ./ sigma) .^ 2);
if size(posterior_deviant, 3) < 2
    positive_update = posterior_deviant;
else
    delta = diff(posterior_deviant, 1, 3);
    positive_update = cat(3, posterior_deviant(:, :, 1), max(delta, 0));
end
model_tensor = cat(2, posterior_deviant, positive_update);
end

function loo_mean = local_leave_one_out_mean(tensor)
% LOCAL_LEAVE_ONE_OUT_MEAN Mean across identities excluding the current row.
n_identity = size(tensor, 1);
if n_identity < 2
    loo_mean = zeros(size(tensor));
    return;
end
total = sum(tensor, 1);
loo_mean = bsxfun(@rdivide, bsxfun(@minus, total, tensor), n_identity - 1);
end

function sigma = local_timewise_feature_scale(tensor, Cfg)
% LOCAL_TIMEWISE_FEATURE_SCALE Estimate per-time residual spread with a floor.
scale_floor = local_spatial_sigma(Cfg);
sigma = std(tensor, 0, 1);
sigma = max(sigma, scale_floor);
sigma = repmat(sigma, [size(tensor, 1), 1, 1]);
end

function sigma = local_spatial_sigma(Cfg)
% LOCAL_SPATIAL_SIGMA Resolve a conservative spatial uncertainty floor.
sigma = 0.25;
if isfield(Cfg, 'dot_w') && ~isempty(Cfg.dot_w)
    sigma = max(sigma, double(Cfg.dot_w) / 2);
elseif isfield(Cfg, 'dotWidth') && ~isempty(Cfg.dotWidth)
    sigma = max(sigma, double(Cfg.dotWidth) / 2);
end
if ~isfinite(sigma) || sigma <= 0
    sigma = 0.25;
end
end

function rect_area = local_motion_rect_area(Cfg)
% LOCAL_MOTION_RECT_AREA Return the visual-degree area used by the uniform branch.
rect_area = 100;
if isfield(Cfg, 'rectSize') && numel(Cfg.rectSize) >= 2
    rect_size = double(Cfg.rectSize(:)');
    rect_area = rect_size(1) .* rect_size(2);
end
if ~isfinite(rect_area) || rect_area <= 0
    rect_area = 100;
end
end

function curvature = local_curvature_tensor(sorted_trials, n_time)
% LOCAL_CURVATURE_TENSOR Load saved curvyness or derive signed turn rate.
n_identity = numel(sorted_trials);
curvature = zeros(n_identity, 1, n_time);
for i_id = 1:n_identity
    if isfield(sorted_trials(i_id), 'curvyness') && ...
            ~isempty(sorted_trials(i_id).curvyness)
        values = double(sorted_trials(i_id).curvyness(:));
    else
        values = local_derive_curvature_from_xy(sorted_trials(i_id).xy);
    end
    if numel(values) ~= n_time
        error('Curvature length mismatch for identity index %d: model=%d neural=%d.', ...
            i_id, numel(values), n_time);
    end
    curvature(i_id, 1, :) = reshape(values, 1, 1, []);
end
end

function curvature = local_derive_curvature_from_xy(xy)
% LOCAL_DERIVE_CURVATURE_FROM_XY Fallback signed turn-rate from path samples.
xy = double(xy);
if size(xy, 1) < 3
    curvature = zeros(size(xy, 1), 1);
    return;
end
dx = diff(xy(:, 1), 1, 1);
dy = diff(xy(:, 2), 1, 1);
angles_deg = rad2deg(unwrap(atan2(dy, dx)));
turn_rate = [0; diff(angles_deg)];
curvature = [turn_rate(1); turn_rate(:)];
curvature = curvature(1:size(xy, 1));
end
