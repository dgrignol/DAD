%% Stimuli Generation V23 (config-driven one-dot occlusion, no source subject files)
% Script: stimuli_generation_V23_NotWiredOccluder.m
% Author: Dami (V23 fixation-collision modes by Codex)
%
% Purpose:
%   Generate a complete one-dot occlusion dataset directly from config
%   parameters, without loading any source subject trajectory MAT files.
%
%   This script keeps the fixed-frame occlusion geometry behavior introduced
%   in V20 while replacing source-file transforms with de novo trajectory
%   synthesis (v17-style controls):
%     - one baseline nondeviant trajectory per sequence,
%     - one deviant trajectory per sequence (turn/curvature change),
%     - one predicted deviant baseline per sequence.
%
%   From those generated trajectories, the script exports:
%     1) always_visible
%     2) occluded_nondeviant
%     3) occluded_deviant
%   plus occluded_deviant_predicted in MovDot_SubXX_predicted.mat.
%
% Optional central-fixation collision handling:
%   - configured through Config_stimuli_generation_V23_NotWiredOccluder.fixationCollisionMode.
%   - modes:
%       * 'off'   : no special handling.
%       * 'retry' : reject colliding candidates and resample.
%       * 'move'  : minimally translate colliding paired trajectories out of
%                   the fixation exclusion zone when feasible, otherwise retry.
%   - exclusion radius is controlled by
%     Config_stimuli_generation_V23_NotWiredOccluder.fixationExclusionRadiusDeg.
%
% Fixed occlusion timeline:
%   - deviance frame fixed at 130 (default, configurable via Config_stimuli_generation_V23_NotWiredOccluder)
%   - first fully occluded frame fixed at 130
%   - fully occluded through frame 190 (inclusive)
%   - nominal reappearance search starts at frame 191
%     (actual reappearance metadata is computed geometrically)
%
% Occlusion geometry model (notwired_circle default):
%   - Start occluder:
%       * uses one fixed circle (center+radius) per trial
%       * occlusion_start_frame is inferred geometrically as the first frame
%         where this full circle is not fully visible anymore
%       * tangent at contact is perpendicular to path direction
%       * center side is chosen from local curvature sign
%       * radius is solved from exact frame constraints so full occlusion
%         lasts through fixedOcclusionEndFrame and reappearance starts after
%         the fixed interval (same timing semantics as v12/v22).
%   - Deviant reanchor occluder:
%       * same radius as start occluder
%       * at reappearance, center is recomputed from direction sampled at
%         fixedOcclusionEndFrame - 1.
%
% Trigger/event metadata compatibility:
%   The exported trial fields are aligned with:
%     - experiment/MoveDot1_experiment_occlusion_v1.m
%     - experiment/trigger_codes.md
%
% Usage example (interactive from experiment/):
%   addpath('lib');
%   stimuli_generation_V23_NotWiredOccluder;
%
% Usage example (set fixed-frame overrides before running):
%   addpath('lib');
%   fixedDevianceFrame = 130;
%   fixedOcclusionEndFrame = 190;
%   overwriteExisting = true;
%   stimuli_generation_V23_NotWiredOccluder;
%
% Usage example (enable fixation-cross avoidance in config and run):
%   addpath('lib');
%   % In Config_stimuli_generation_V23_NotWiredOccluder:
%   %   fixationCollisionMode = 'move'; % or 'retry' / 'off'
%   stimuli_generation_V23_NotWiredOccluder;
%
% Inputs:
%   - Config_stimuli_generation_V23_NotWiredOccluder class on MATLAB path.
%   - Interactive target subject ID (always prompted, no default target).
%
% Outputs:
%   - input_files/MovDot_SubXX_V23_NotWiredOccluder.mat
%       with xySeqs containing always_visible / occluded_nondeviant /
%       occluded_deviant and full occlusion metadata fields.
%   - input_files/MovDot_SubXX_V23_NotWiredOccluder_predicted.mat
%       with xySeqsPredicted containing occluded_deviant_predicted trials
%       and matching metadata.
%
% Key assumptions:
%   - One-dot only: xy is frames x 2 ([x y]) in visual degrees.
%   - Direction variance contains at least one nondeviant value (0) and one
%     deviant value (>0); the first deviant value is used for generation.
%   - Trial pairing is sequence-index based and one-to-one across conditions.
%   - Output directories are relative to experiment/ unless overridden by
%     Config_stimuli_generation_V23_NotWiredOccluder.

addpath('lib/');

%% Setup and runtime parameter defaults
% Data flow: caller overrides + Config_stimuli_generation_V23_NotWiredOccluder constants -> validated controls.
clearvars -except overwriteExisting fixedDevianceFrame fixedOcclusionEndFrame targetSubjectID;
clc;

if ~(exist('fixedDevianceFrame', 'var') && ~isempty(fixedDevianceFrame))
    fixedDevianceFrame = Config_stimuli_generation_V23_NotWiredOccluder.fixedDevianceFrame;
end
if ~(exist('fixedOcclusionEndFrame', 'var') && ~isempty(fixedOcclusionEndFrame))
    fixedOcclusionEndFrame = Config_stimuli_generation_V23_NotWiredOccluder.fixedOcclusionEndFrame;
end
if ~(exist('overwriteExisting', 'var') && ~isempty(overwriteExisting))
    overwriteExisting = true;
end

if ~isscalar(fixedDevianceFrame) || ~isscalar(fixedOcclusionEndFrame) || ...
        isnan(fixedDevianceFrame) || isnan(fixedOcclusionEndFrame)
    error('fixedDevianceFrame and fixedOcclusionEndFrame must be numeric scalars.');
end
fixedDevianceFrame = round(fixedDevianceFrame);
fixedOcclusionEndFrame = round(fixedOcclusionEndFrame);
overwriteExisting = logical(overwriteExisting);

%% Resolve target subject ID (workspace override or interactive prompt)
% Data flow: optional caller-provided subject ID -> fallback dialog -> validated ID.
if ~(exist('targetSubjectID', 'var') && ~isempty(targetSubjectID))
    prompt = {'Target subject ID (new config-driven occlusion dataset):'};
    dlgTitle = 'V23 target subject';
    dims = [1 70];
    definput = {''}; % no default target by design.
    userInput = inputdlg(prompt, dlgTitle, dims, definput);
    if isempty(userInput)
        disp('Canceled by user.');
        return;
    end
    targetSubjectID = str2double(userInput{1});
end

if isnan(targetSubjectID) || targetSubjectID < 0
    error('Target subject ID must be a numeric value >= 0.');
end
targetSubjectID = round(targetSubjectID);

%% Resolve output paths and overwrite guard
% Data flow: targetSubjectID + Config_stimuli_generation_V23_NotWiredOccluder directories -> output files.
inputDir = Config_stimuli_generation_V23_NotWiredOccluder.inputDirectory;
if ~exist(inputDir, 'dir')
    mkdir(inputDir);
end

targetObservedFile = fullfile(inputDir, sprintf( ...
    Config_stimuli_generation_V23_NotWiredOccluder.observedFilePattern, targetSubjectID));
targetPredictedFile = fullfile(inputDir, sprintf( ...
    Config_stimuli_generation_V23_NotWiredOccluder.predictedFilePattern, targetSubjectID));

if isfile(targetObservedFile) || isfile(targetPredictedFile)
    doOverwrite = logical(overwriteExisting);
    if ~doOverwrite
        choice = '';
        while ~ismember(lower(choice), {'o', 'c'})
            choice = input(sprintf(['Target file(s) already exist for Sub%02d. ' ...
                '[O]verwrite or [C]ancel? '], targetSubjectID), 's');
        end
        if strcmpi(choice, 'c')
            disp('Canceled by user.');
            return;
        end
    end
end

%% Build generation controls from Config_stimuli_generation_V23_NotWiredOccluder
% Data flow: Config_stimuli_generation_V23_NotWiredOccluder constants -> trajectory synthesis controls.
stimulusTypeConfig = Config_stimuli_generation_V23_NotWiredOccluder.likelihood;
framesPerTrial = round(Config_stimuli_generation_V23_NotWiredOccluder.trialDuration * Config_stimuli_generation_V23_NotWiredOccluder.frameFrequency);
nTrials = Config_stimuli_generation_V23_NotWiredOccluder.trialsPerCondition;
fps = double(Config_stimuli_generation_V23_NotWiredOccluder.frameFrequency);
alphaFadeInFrames = max(1, round(0.25 * fps));
dotRadiusDeg = double(Config_stimuli_generation_V23_NotWiredOccluder.dotWidth) / 2;

if nTrials < 1
    error('Config_stimuli_generation_V23_NotWiredOccluder.trialsPerCondition must be >= 1.');
end
if framesPerTrial < 2
    error('Computed framesPerTrial (trialDuration * frameFrequency) must be >= 2.');
end
if fixedDevianceFrame < 1 || fixedDevianceFrame > framesPerTrial
    error('fixedDevianceFrame (%d) must be within [1, %d].', fixedDevianceFrame, framesPerTrial);
end
if fixedOcclusionEndFrame < fixedDevianceFrame || fixedOcclusionEndFrame > framesPerTrial
    error(['fixedOcclusionEndFrame (%d) must be within [%d, %d] and must not ' ...
        'precede fixedDevianceFrame.'], fixedOcclusionEndFrame, fixedDevianceFrame, framesPerTrial);
end

directionVariance = double(stimulusTypeConfig.directionVariance(:)');
if ~any(directionVariance == 0)
    error('Config_stimuli_generation_V23_NotWiredOccluder.likelihood.directionVariance must include nondeviant value 0.');
end
if ~any(directionVariance > 0)
    error('Config_stimuli_generation_V23_NotWiredOccluder.likelihood.directionVariance must include at least one deviant value > 0.');
end
deviantConditionCode = directionVariance(find(directionVariance > 0, 1, 'first'));
if numel(directionVariance(directionVariance > 0)) > 1
    warning('Multiple deviant directionVariance values found; using first deviant value %g.', ...
        deviantConditionCode);
end

if numel(stimulusTypeConfig.pathDuration) > 1
    warning(['Config_stimuli_generation_V23_NotWiredOccluder.likelihood.pathDuration has multiple values. ' ...
        'V23 uses the first value for one-dot occlusion generation.']);
end
pathDurationSec = double(stimulusTypeConfig.pathDuration(1));

generationCfg = struct();
generationCfg.framesPerTrial = framesPerTrial;
generationCfg.dotSpeedDegPerFrame = double(Config_stimuli_generation_V23_NotWiredOccluder.dotSpeedDegPerFrame);
generationCfg.minPosDeg = [Config_stimuli_generation_V23_NotWiredOccluder.dotWidth / 2, Config_stimuli_generation_V23_NotWiredOccluder.dotWidth / 2];
generationCfg.maxPosDeg = double(Config_stimuli_generation_V23_NotWiredOccluder.dotRectSize) - generationCfg.minPosDeg;
generationCfg.initialCurvatureWindows = local_parse_interval_windows( ...
    Config_stimuli_generation_V23_NotWiredOccluder.initialCurvatureWindows, 'Config_stimuli_generation_V23_NotWiredOccluder.initialCurvatureWindows', false, -inf, inf);
generationCfg.deviantCurvatureWindows = local_parse_interval_windows( ...
    Config_stimuli_generation_V23_NotWiredOccluder.deviantCurvatureWindows, 'Config_stimuli_generation_V23_NotWiredOccluder.deviantCurvatureWindows', false, -inf, inf);
generationCfg.deviantTurnWindowsDeg = local_parse_interval_windows( ...
    stimulusTypeConfig.deviantSignedTurnWindows, ...
    'Config_stimuli_generation_V23_NotWiredOccluder.likelihood.deviantSignedTurnWindows', true, -180, 180);
generationCfg.flipCurvatureOnDeviant = logical(Config_stimuli_generation_V23_NotWiredOccluder.flipCurvatureOnDeviant);
generationCfg.randomizeCurvatureOnDeviant = logical(Config_stimuli_generation_V23_NotWiredOccluder.randomizeCurvatureOnDeviant);
generationCfg.pathDurationSec = pathDurationSec;
generationCfg.maxAttemptsPerTrial = 60000;
generationCfg.fixationCollisionMode = lower(strtrim(char( ...
    Config_stimuli_generation_V23_NotWiredOccluder.fixationCollisionMode)));
generationCfg.fixationExclusionRadiusDeg = double( ...
    Config_stimuli_generation_V23_NotWiredOccluder.fixationExclusionRadiusDeg);
generationCfg.fixationCenterDeg = double(Config_stimuli_generation_V23_NotWiredOccluder.dotRectSize(:)') / 2;
generationCfg.fixationMovePaddingDeg = double( ...
    Config_stimuli_generation_V23_NotWiredOccluder.fixationMovePaddingDeg);
generationCfg.fixationMoveDirectionSamples = double( ...
    Config_stimuli_generation_V23_NotWiredOccluder.fixationMoveDirectionSamples);
generationCfg.fixationMoveShiftSamples = double( ...
    Config_stimuli_generation_V23_NotWiredOccluder.fixationMoveShiftSamples);

% Validate fixation-collision controls once, then enforce them per trial.
if numel(generationCfg.fixationCenterDeg) ~= 2
    error('Config_stimuli_generation_V23_NotWiredOccluder.dotRectSize must be [width height].');
end
if ~ismember(generationCfg.fixationCollisionMode, {'off', 'retry', 'move'})
    error(['Config_stimuli_generation_V23_NotWiredOccluder.fixationCollisionMode must be one of: ' ...
        '''off'', ''retry'', ''move''.']);
end
if ~isfinite(generationCfg.fixationExclusionRadiusDeg) || generationCfg.fixationExclusionRadiusDeg < 0
    error('Config_stimuli_generation_V23_NotWiredOccluder.fixationExclusionRadiusDeg must be finite and >= 0.');
end
if ~isfinite(generationCfg.fixationMovePaddingDeg) || generationCfg.fixationMovePaddingDeg < 0
    error('Config_stimuli_generation_V23_NotWiredOccluder.fixationMovePaddingDeg must be finite and >= 0.');
end
if generationCfg.fixationMoveDirectionSamples < 4 || ...
        floor(generationCfg.fixationMoveDirectionSamples) ~= generationCfg.fixationMoveDirectionSamples
    error('Config_stimuli_generation_V23_NotWiredOccluder.fixationMoveDirectionSamples must be an integer >= 4.');
end
if generationCfg.fixationMoveShiftSamples < 10 || ...
        floor(generationCfg.fixationMoveShiftSamples) ~= generationCfg.fixationMoveShiftSamples
    error('Config_stimuli_generation_V23_NotWiredOccluder.fixationMoveShiftSamples must be an integer >= 10.');
end

% Match v17 geometric floor logic to avoid near-straight paths that cannot
% fit inside the movement rectangle under fixed speed/frame constraints.
maxTurnRadiusDeg = min(double(Config_stimuli_generation_V23_NotWiredOccluder.dotRectSize) - double(Config_stimuli_generation_V23_NotWiredOccluder.dotWidth)) / 2;
if maxTurnRadiusDeg <= 0
    error('Config_stimuli_generation_V23_NotWiredOccluder.dotRectSize minus dotWidth must be > 0.');
end
generationCfg.minimumAbsCurvatureDeg = rad2deg(generationCfg.dotSpeedDegPerFrame / maxTurnRadiusDeg);

%% Seed RNG from target subject for reproducibility
% Data flow: targetSubjectID -> RNG state -> deterministic trial generation.
rng(targetSubjectID);
rngState = rng;

%% Generate base one-dot trajectories (nondeviant, deviant, predicted baseline)
% Data flow: generationCfg + fixedDevianceFrame -> paired base trajectories.
baseNondeviantCell = cell(nTrials, 1);
baseDeviantCell = cell(nTrials, 1);
basePredictedCell = cell(nTrials, 1);

for iTrial = 1:nTrials
    [trialNondev, trialDev, trialPred] = local_generate_trial_triplet( ...
        iTrial, generationCfg, fixedDevianceFrame, deviantConditionCode);
    baseNondeviantCell{iTrial} = trialNondev;
    baseDeviantCell{iTrial} = trialDev;
    basePredictedCell{iTrial} = trialPred;
end

baseNondeviant = vertcat(baseNondeviantCell{:});
baseDeviant = vertcat(baseDeviantCell{:});
basePredicted = vertcat(basePredictedCell{:});

%% Build observed output trials for occlusion paradigm
% Data flow: paired base trials -> condition metadata -> flattened xySeqs.
observedCell = cell(3 * nTrials, 1);
splicedObservedDeviant = cell(nTrials, 1);

for iTrial = 1:nTrials
    trialNondev = baseNondeviant(iTrial);
    trialDevRaw = baseDeviant(iTrial);

    devFrame = local_infer_deviance_frame(framesPerTrial, fixedDevianceFrame);
    timing = local_build_occlusion_timing( ...
        framesPerTrial, devFrame, fixedOcclusionEndFrame, alphaFadeInFrames);
    trialDevSpliced = local_build_spliced_occluded_deviant(trialNondev, trialDevRaw, devFrame);

    notwiredNondev = local_build_notwired_occluder_metadata( ...
        trialNondev.xy, devFrame, timing, dotRadiusDeg, fixedOcclusionEndFrame, false);
    notwiredDev = local_build_notwired_occluder_metadata( ...
        trialDevSpliced.xy, devFrame, timing, dotRadiusDeg, fixedOcclusionEndFrame, true);

    [alphaVisible, geomVisible] = local_build_visible_profiles(framesPerTrial);
    timingNondev = timing;
    timingNondev.occlusion_start_frame = notwiredNondev.occlusion_start_frame;
    [alphaOccNondev, geomOccNondev] = local_build_alpha_profiles(framesPerTrial, timingNondev);
    timingDev = timing;
    timingDev.occlusion_start_frame = notwiredDev.occlusion_start_frame;
    [alphaOccDev, geomOccDev] = local_build_alpha_profiles(framesPerTrial, timingDev);

    rowVisible = iTrial;
    rowOccNondev = nTrials + iTrial;
    rowOccDev = 2 * nTrials + iTrial;

    observedCell{rowVisible} = local_attach_occlusion_fields( ...
        trialNondev, -1, 'always_visible', false, ...
        timing, notwiredNondev, alphaVisible, geomVisible);

    observedCell{rowOccNondev} = local_attach_occlusion_fields( ...
        trialNondev, 0, 'occluded_nondeviant', true, ...
        timingNondev, notwiredNondev, alphaOccNondev, geomOccNondev);

    observedCell{rowOccDev} = local_attach_occlusion_fields( ...
        trialDevSpliced, deviantConditionCode, 'occluded_deviant', true, ...
        timingDev, notwiredDev, alphaOccDev, geomOccDev);

    splicedObservedDeviant{iTrial} = observedCell{rowOccDev};
end
xySeqs = vertcat(observedCell{:});

%% Build predicted output trials (occluded deviant predicted branch)
% Data flow: predicted baseline + observed deviant prefix -> xySeqsPredicted.
predictedCell = cell(nTrials, 1);
for iTrial = 1:nTrials
    trialObservedDev = splicedObservedDeviant{iTrial};
    trialPredictedDev = basePredicted(iTrial);

    devFrame = local_infer_deviance_frame(framesPerTrial, fixedDevianceFrame);
    timing = local_build_occlusion_timing( ...
        framesPerTrial, devFrame, fixedOcclusionEndFrame, alphaFadeInFrames);
    notwiredPred = local_build_notwired_occluder_metadata( ...
        trialPredictedDev.xy, devFrame, timing, dotRadiusDeg, fixedOcclusionEndFrame, true);
    timingPred = timing;
    timingPred.occlusion_start_frame = notwiredPred.occlusion_start_frame;
    [alphaOcc, geomOcc] = local_build_alpha_profiles(framesPerTrial, timingPred);

    % Modeling assumption retained from v18-v20:
    % predicted and observed are forced identical until reappearance starts.
    adjustedPred = trialPredictedDev;
    rs = timing.reappearance_start_frame;
    if rs > 1
        adjustedPred.xy(1:(rs - 1), :) = trialObservedDev.xy(1:(rs - 1), :);
    end

    predictedCell{iTrial} = local_attach_occlusion_fields( ...
        adjustedPred, deviantConditionCode, 'occluded_deviant_predicted', true, ...
        timingPred, notwiredPred, alphaOcc, geomOcc);
end
xySeqsPredicted = vertcat(predictedCell{:});

%% Build Cfg and reproducibility metadata
% Data flow: Config_stimuli_generation_V23_NotWiredOccluder snapshot + runtime values -> Cfg + repro outputs.
Cfg = struct();
Cfg.dpf = Config_stimuli_generation_V23_NotWiredOccluder.dotSpeedDegPerFrame;
Cfg.fps = Config_stimuli_generation_V23_NotWiredOccluder.frameFrequency;
Cfg.Stimulitype = Config_stimuli_generation_V23_NotWiredOccluder.stimulusType;
Cfg.dot_w = Config_stimuli_generation_V23_NotWiredOccluder.dotWidth;
Cfg.rectSize = Config_stimuli_generation_V23_NotWiredOccluder.dotRectSize;
Cfg.DirChange = stimulusTypeConfig.directionChange;

configProps = properties('Config_stimuli_generation_V23_NotWiredOccluder');
configSnapshot = struct();
for iProp = 1:numel(configProps)
    propName = configProps{iProp};
    configSnapshot.(propName) = Config_stimuli_generation_V23_NotWiredOccluder.(propName);
end

repro = struct();
repro.script = struct( ...
    'name', mfilename, ...
    'version', 'V23_fixationCollisionModes', ...
    'parameters', struct( ...
        'fixedDevianceFrame', fixedDevianceFrame, ...
        'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
        'alphaFadeInFrames', alphaFadeInFrames, ...
        'deviantConditionCode', deviantConditionCode, ...
        'pathDurationSec', pathDurationSec, ...
        'minimumAbsCurvatureDeg', generationCfg.minimumAbsCurvatureDeg, ...
        'maxAttemptsPerTrial', generationCfg.maxAttemptsPerTrial, ...
        'fixationCollisionMode', generationCfg.fixationCollisionMode, ...
        'fixationExclusionRadiusDeg', generationCfg.fixationExclusionRadiusDeg, ...
        'fixationMovePaddingDeg', generationCfg.fixationMovePaddingDeg, ...
        'fixationMoveDirectionSamples', generationCfg.fixationMoveDirectionSamples, ...
        'fixationMoveShiftSamples', generationCfg.fixationMoveShiftSamples));
repro.inputs = struct('targetSubjectID', targetSubjectID);
repro.rng = rngState;
repro.config = configSnapshot;
repro.v23_notwired_occlusion = struct( ...
    'script', mfilename, ...
    'targetSubjectID', targetSubjectID, ...
    'fixedDevianceFrame', fixedDevianceFrame, ...
    'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
    'nominalReappearanceFrame', fixedOcclusionEndFrame + 1, ...
    'dotRadiusDeg', dotRadiusDeg, ...
    'deviantConstructionRule', 'config_generated_with_nondeviant_prefix_locked_at_deviance', ...
    'notWiredRule', ['radius_solved_from_exact_trajectory_constraints;' ...
        'start_tangent_perpendicular_to_contact_direction;' ...
        'deviant_reanchor_uses_fixedOcclusionEndFrame_minus_one_direction'], ...
    'defaultModel', 'notwired_circle', ...
    'optionalModel', 'alpha', ...
    'conditionCodes', struct('always_visible', -1, 'occluded_nondeviant', 0, ...
        'occluded_deviant', deviantConditionCode));

%% Save outputs
% Data flow: generated structs -> MAT files in Config_stimuli_generation_V23_NotWiredOccluder.inputDirectory.
save(targetObservedFile, 'xySeqs', 'Cfg', 'repro');
save(targetPredictedFile, 'xySeqsPredicted', 'Cfg');

fprintf('Saved observed occlusion dataset: %s\n', targetObservedFile);
fprintf('Saved predicted occlusion dataset: %s\n', targetPredictedFile);
fprintf('Trials per condition: %d (always_visible / occluded_nondeviant / occluded_deviant).\n', nTrials);
fprintf(['Fixed occlusion timing: full occlusion starts at deviance=%d, ' ...
    'last required full occlusion=%d, nominal reappearance search start=%d.\n'], ...
    fixedDevianceFrame, fixedOcclusionEndFrame, fixedOcclusionEndFrame + 1);

%% Local helper functions
function [trialNondev, trialDev, trialPred] = local_generate_trial_triplet( ...
        sequenceIndex, generationCfg, devFrame, deviantConditionCode)
% LOCAL_GENERATE_TRIAL_TRIPLET Generate paired one-dot paths for one sequence.
%
% Data flow:
%   sampled direction/curvature/turn controls -> relative paths ->
%   feasible placement -> absolute one-dot trial structs.
for attempt = 1:generationCfg.maxAttemptsPerTrial
    initialDirectionDeg = Utils.RandAngleDegree();
    baselineCurvatureDeg = local_sample_from_windows(generationCfg.initialCurvatureWindows, 1);

    % Clamp tiny curvature magnitudes to the geometric feasibility floor.
    if abs(baselineCurvatureDeg) < generationCfg.minimumAbsCurvatureDeg
        if baselineCurvatureDeg == 0
            baselineCurvatureDeg = generationCfg.minimumAbsCurvatureDeg * Utils.RandPosNeg();
        else
            baselineCurvatureDeg = sign(baselineCurvatureDeg) * generationCfg.minimumAbsCurvatureDeg;
        end
    end

    nSteps = generationCfg.framesPerTrial - 1;
    turnNondevDeg = zeros(nSteps, 1);
    turnDeviantDeg = zeros(nSteps, 1);
    if nSteps > 0 && devFrame > 1
        sampledTurnDeg = local_sample_deviant_turn( ...
            generationCfg.deviantTurnWindowsDeg, deviantConditionCode);
        turnDeviantDeg(devFrame - 1) = sampledTurnDeg;
    end

    curvatureNondevDeg = repmat(baselineCurvatureDeg, nSteps, 1);
    curvatureDeviantDeg = curvatureNondevDeg;
    if nSteps > 0 && devFrame > 1
        onsetIndex = devFrame - 1;
        if generationCfg.randomizeCurvatureOnDeviant
            postCurvatureDeg = local_sample_from_windows(generationCfg.deviantCurvatureWindows, 1);
            curvatureDeviantDeg(onsetIndex:end) = postCurvatureDeg;
        elseif generationCfg.flipCurvatureOnDeviant
            curvatureDeviantDeg(onsetIndex:end) = -baselineCurvatureDeg;
        end
    end

    directionNondevDeg = local_integrate_direction( ...
        initialDirectionDeg, turnNondevDeg, curvatureNondevDeg);
    directionDeviantDegRaw = local_integrate_direction( ...
        initialDirectionDeg, turnDeviantDeg, curvatureDeviantDeg);

    relativeNondevXY = local_integrate_xy(directionNondevDeg, generationCfg.dotSpeedDegPerFrame);
    relativeDeviantRawXY = local_integrate_xy(directionDeviantDegRaw, generationCfg.dotSpeedDegPerFrame);

    % Explicitly enforce nondeviant-prefix identity through devFrame.
    relativeDeviantXY = relativeDeviantRawXY;
    relativeDeviantXY(1:devFrame, :) = relativeNondevXY(1:devFrame, :);
    if devFrame < generationCfg.framesPerTrial
        suffix = relativeDeviantRawXY(devFrame:end, :);
        shiftVec = relativeDeviantXY(devFrame, :) - suffix(1, :);
        relativeDeviantXY(devFrame:end, :) = suffix + shiftVec;
    end

    % Place both trajectories together so every frame fits the movement bounds.
    [isFeasible, startOffset] = local_sample_feasible_start( ...
        relativeNondevXY, relativeDeviantXY, generationCfg.minPosDeg, generationCfg.maxPosDeg);
    if ~isFeasible
        continue;
    end

    nondevXY = relativeNondevXY + startOffset;
    deviantXY = relativeDeviantXY + startOffset;
    predictedXY = relativeNondevXY + startOffset; % no-deviant baseline for deviant branch.

    if ~local_paths_within_bounds(nondevXY, generationCfg.minPosDeg, generationCfg.maxPosDeg) || ...
            ~local_paths_within_bounds(deviantXY, generationCfg.minPosDeg, generationCfg.maxPosDeg)
        continue;
    end

    % Optional fixation-zone handling strategy:
    %   - off: keep candidate as-is
    %   - retry: reject colliding candidate
    %   - move: minimally translate candidate out of exclusion zone, else retry
    collisionNondev = ~local_path_avoids_fixation_zone( ...
        nondevXY, generationCfg.fixationCenterDeg, generationCfg.fixationExclusionRadiusDeg);
    collisionDeviant = ~local_path_avoids_fixation_zone( ...
        deviantXY, generationCfg.fixationCenterDeg, generationCfg.fixationExclusionRadiusDeg);
    if collisionNondev || collisionDeviant
        switch generationCfg.fixationCollisionMode
            case 'off'
                % Keep the colliding candidate unchanged.
            case 'retry'
                continue;
            case 'move'
                [nondevMoved, deviantMoved, predictedMoved, moveSucceeded] = ...
                    local_translate_paths_out_of_fixation_zone( ...
                    nondevXY, deviantXY, predictedXY, ...
                    generationCfg.fixationCenterDeg, generationCfg.fixationExclusionRadiusDeg, ...
                    generationCfg.minPosDeg, generationCfg.maxPosDeg, ...
                    generationCfg.fixationMovePaddingDeg, ...
                    generationCfg.fixationMoveDirectionSamples, ...
                    generationCfg.fixationMoveShiftSamples);
                if ~moveSucceeded
                    continue;
                end
                nondevXY = nondevMoved;
                deviantXY = deviantMoved;
                predictedXY = predictedMoved;
            otherwise
                error('Unexpected fixationCollisionMode: %s', generationCfg.fixationCollisionMode);
        end
    end

    pathAllNondev = zeros(generationCfg.framesPerTrial, 1);
    pathAllNondev(1) = 1;
    pathAllDeviant = pathAllNondev;
    pathAllDeviant(devFrame) = 1;

    curvynessNondev = [baselineCurvatureDeg; curvatureNondevDeg];
    curvynessDeviant = [baselineCurvatureDeg; curvatureDeviantDeg];

    trialNondev = local_build_base_trial( ...
        0, generationCfg.pathDurationSec, sequenceIndex, nondevXY, ...
        pathAllNondev, curvynessNondev, directionNondevDeg);
    trialDev = local_build_base_trial( ...
        deviantConditionCode, generationCfg.pathDurationSec, sequenceIndex, deviantXY, ...
        pathAllDeviant, curvynessDeviant, directionDeviantDegRaw);
    trialPred = local_build_base_trial( ...
        deviantConditionCode, generationCfg.pathDurationSec, sequenceIndex, predictedXY, ...
        pathAllNondev, curvynessNondev, directionNondevDeg);
    return;
end

error(['Failed to generate a feasible trial after %d attempts (sequence %d). ' ...
    'Check Config_stimuli_generation_V23_NotWiredOccluder geometry and curvature settings.'], ...
    generationCfg.maxAttemptsPerTrial, sequenceIndex);
end

function trial = local_build_base_trial(conditionCode, pathDurationSec, sequenceIndex, ...
        xy, pathAll, curvyness, angleDirection)
% LOCAL_BUILD_BASE_TRIAL Build base one-dot trial struct before occlusion fields.
trial = struct();
trial.condition = conditionCode;
trial.PredictionRange = pathDurationSec;
trial.sequence = sequenceIndex;
trial.xy = xy;
trial.pathAll = pathAll;
trial.curvyness = curvyness;
trial.AngleDirection = angleDirection;
end

function directionsDeg = local_integrate_direction(initialDirectionDeg, turnDeg, curvatureDeg)
% LOCAL_INTEGRATE_DIRECTION Integrate per-step turn+curvature into frame directions.
nFrames = numel(turnDeg) + 1;
directionsDeg = zeros(nFrames, 1);
directionsDeg(1) = initialDirectionDeg;
if nFrames > 1
    directionsDeg(2:end) = initialDirectionDeg + cumsum(turnDeg + curvatureDeg, 1);
end
end

function xy = local_integrate_xy(directionDeg, dotSpeedDegPerFrame)
% LOCAL_INTEGRATE_XY Integrate frame directions into relative x/y trajectory.
nFrames = numel(directionDeg);
xy = zeros(nFrames, 2);
if nFrames > 1
    steps = [cosd(directionDeg(2:end)), sind(directionDeg(2:end))] * dotSpeedDegPerFrame;
    xy(2:end, :) = cumsum(steps, 1);
end
end

function [isFeasible, startOffset] = local_sample_feasible_start(nondevXY, devXY, minPosDeg, maxPosDeg)
% LOCAL_SAMPLE_FEASIBLE_START Sample one start offset that keeps both paths in bounds.
allX = [nondevXY(:, 1); devXY(:, 1)];
allY = [nondevXY(:, 2); devXY(:, 2)];

xRange = [minPosDeg(1) - min(allX), maxPosDeg(1) - max(allX)];
yRange = [minPosDeg(2) - min(allY), maxPosDeg(2) - max(allY)];

if xRange(1) > xRange(2) || yRange(1) > yRange(2)
    isFeasible = false;
    startOffset = [NaN, NaN];
    return;
end

startOffset = [ ...
    xRange(1) + rand(1) * (xRange(2) - xRange(1)), ...
    yRange(1) + rand(1) * (yRange(2) - yRange(1))];
isFeasible = true;
end

function inBounds = local_paths_within_bounds(xy, minPosDeg, maxPosDeg)
% LOCAL_PATHS_WITHIN_BOUNDS True when every point is inside [minPosDeg, maxPosDeg].
inBounds = all(xy(:, 1) >= minPosDeg(1) & xy(:, 1) <= maxPosDeg(1) & ...
    xy(:, 2) >= minPosDeg(2) & xy(:, 2) <= maxPosDeg(2));
end

function avoids = local_path_avoids_fixation_zone(xy, fixationCenterDeg, exclusionRadiusDeg)
% LOCAL_PATH_AVOIDS_FIXATION_ZONE True when all frames stay outside exclusion radius.
if exclusionRadiusDeg <= 0
    avoids = true;
    return;
end
distToFixation = sqrt(sum((double(xy(:, 1:2)) - fixationCenterDeg) .^ 2, 2));
avoids = all(distToFixation >= exclusionRadiusDeg);
end

function [nondevMoved, deviantMoved, predictedMoved, moveSucceeded] = ...
        local_translate_paths_out_of_fixation_zone( ...
        nondevXY, deviantXY, predictedXY, fixationCenterDeg, exclusionRadiusDeg, ...
        minPosDeg, maxPosDeg, paddingDeg, directionSamples, shiftSamples)
% LOCAL_TRANSLATE_PATHS_OUT_OF_FIXATION_ZONE
% Try a minimal rigid translation that makes both paired paths leave
% the fixation exclusion zone while remaining in movement bounds.
%
% Data flow:
%   colliding absolute trajectories -> candidate translation search ->
%   translated trajectories or failure (caller retries new sample).
nondevMoved = nondevXY;
deviantMoved = deviantXY;
predictedMoved = predictedXY;
moveSucceeded = false;

if exclusionRadiusDeg <= 0
    moveSucceeded = true;
    return;
end

if local_path_avoids_fixation_zone(nondevXY, fixationCenterDeg, exclusionRadiusDeg) && ...
        local_path_avoids_fixation_zone(deviantXY, fixationCenterDeg, exclusionRadiusDeg)
    moveSucceeded = true;
    return;
end

allXY = [double(nondevXY(:, 1:2)); double(deviantXY(:, 1:2))];
targetRadius = exclusionRadiusDeg + max(0, paddingDeg);
dirs = local_build_translation_directions(allXY, fixationCenterDeg, directionSamples);

bestFound = false;
bestShiftNorm = inf;
bestShiftVec = [0, 0];

for iDir = 1:size(dirs, 1)
    dirVec = dirs(iDir, :);
    dirNorm = norm(dirVec);
    if dirNorm < 1e-12
        continue;
    end
    unitDir = dirVec / dirNorm;

    sMax = local_max_shift_in_direction(allXY, unitDir, minPosDeg, maxPosDeg);
    if ~isfinite(sMax) || sMax <= 0
        continue;
    end

    sVals = linspace(0, sMax, shiftSamples + 1);
    candidateS = NaN;
    for iS = 2:numel(sVals)
        s = sVals(iS);
        shifted = allXY + s * unitDir;
        d = sqrt(sum((shifted - fixationCenterDeg) .^ 2, 2));
        if all(d >= targetRadius)
            candidateS = s;
            break;
        end
    end

    if ~isnan(candidateS) && candidateS < bestShiftNorm
        bestFound = true;
        bestShiftNorm = candidateS;
        bestShiftVec = candidateS * unitDir;
    end
end

if ~bestFound
    return;
end

nondevMoved = nondevXY + bestShiftVec;
deviantMoved = deviantXY + bestShiftVec;
predictedMoved = predictedXY + bestShiftVec;

if ~local_paths_within_bounds(nondevMoved, minPosDeg, maxPosDeg) || ...
        ~local_paths_within_bounds(deviantMoved, minPosDeg, maxPosDeg) || ...
        ~local_paths_within_bounds(predictedMoved, minPosDeg, maxPosDeg)
    return;
end
if ~local_path_avoids_fixation_zone(nondevMoved, fixationCenterDeg, exclusionRadiusDeg) || ...
        ~local_path_avoids_fixation_zone(deviantMoved, fixationCenterDeg, exclusionRadiusDeg)
    return;
end

moveSucceeded = true;
end

function dirs = local_build_translation_directions(allXY, fixationCenterDeg, directionSamples)
% LOCAL_BUILD_TRANSLATION_DIRECTIONS Build deterministic search directions.
theta = linspace(0, 2 * pi, directionSamples + 1)';
theta(end) = [];
uniformDirs = [cos(theta), sin(theta)];

centroidDir = mean(allXY, 1) - fixationCenterDeg;
if norm(centroidDir) < 1e-12
    centroidDir = [1, 0];
end

distToCenter = sqrt(sum((allXY - fixationCenterDeg) .^ 2, 2));
[~, idxSorted] = sort(distToCenter, 'ascend');
nLocal = min(8, numel(idxSorted));
localDirs = allXY(idxSorted(1:nLocal), :) - fixationCenterDeg;
validLocal = vecnorm(localDirs, 2, 2) > 1e-12;
localDirs = localDirs(validLocal, :);

dirs = [centroidDir; localDirs; uniformDirs];
end

function sMax = local_max_shift_in_direction(allXY, unitDir, minPosDeg, maxPosDeg)
% LOCAL_MAX_SHIFT_IN_DIRECTION Maximum feasible translation along one direction.
allX = allXY(:, 1);
allY = allXY(:, 2);

if abs(unitDir(1)) < 1e-12
    sX = inf;
elseif unitDir(1) > 0
    sX = (maxPosDeg(1) - max(allX)) / unitDir(1);
else
    sX = (minPosDeg(1) - min(allX)) / unitDir(1);
end

if abs(unitDir(2)) < 1e-12
    sY = inf;
elseif unitDir(2) > 0
    sY = (maxPosDeg(2) - max(allY)) / unitDir(2);
else
    sY = (minPosDeg(2) - min(allY)) / unitDir(2);
end

sMax = min([sX, sY]);
if ~isfinite(sMax)
    sMax = 0;
elseif sMax < 0
    sMax = 0;
end
end

function value = local_sample_deviant_turn(turnWindowsDeg, fallbackMagnitudeDeg)
% LOCAL_SAMPLE_DEVIANT_TURN Sample one signed deviant turn in degrees.
%
% Falls back to +-fallbackMagnitudeDeg if no explicit windows are available.
if isempty(turnWindowsDeg)
    value = abs(fallbackMagnitudeDeg) * Utils.RandPosNeg();
    return;
end
value = local_sample_from_windows(turnWindowsDeg, 1);
end

function samples = local_sample_from_windows(windows, sampleCount)
% LOCAL_SAMPLE_FROM_WINDOWS Sample uniformly from a union of numeric intervals.
if isempty(windows)
    error('Cannot sample from empty interval windows.');
end
if sampleCount < 1 || floor(sampleCount) ~= sampleCount
    error('sampleCount must be a positive integer.');
end

windowLengths = windows(:, 2) - windows(:, 1);
if any(windowLengths <= 0)
    error('Each interval window must satisfy max > min.');
end

cumLengths = cumsum(windowLengths);
totalLength = cumLengths(end);
r = rand(sampleCount, 1) * totalLength;
samples = zeros(sampleCount, 1);
for i = 1:sampleCount
    intervalIndex = find(r(i) <= cumLengths, 1, 'first');
    intervalStart = windows(intervalIndex, 1);
    intervalLength = windowLengths(intervalIndex);
    offsetInside = r(i) - (cumLengths(intervalIndex) - intervalLength);
    samples(i) = intervalStart + offsetInside;
end
end

function windows = local_parse_interval_windows(rawWindows, fieldName, ...
        allowEmpty, minAllowed, maxAllowed)
% LOCAL_PARSE_INTERVAL_WINDOWS Validate [min max; ...] interval-window inputs.
if isempty(rawWindows)
    if allowEmpty
        windows = zeros(0, 2);
        return;
    end
    error('%s must not be empty.', fieldName);
end

if ~isnumeric(rawWindows) || size(rawWindows, 2) ~= 2
    error('%s must be an Nx2 numeric matrix of [min max] intervals.', fieldName);
end
if any(~isfinite(rawWindows), 'all')
    error('%s must contain only finite values.', fieldName);
end
if any(rawWindows(:, 2) <= rawWindows(:, 1))
    error('%s intervals must satisfy max > min for every row.', fieldName);
end
if any(rawWindows(:, 1) < minAllowed) || any(rawWindows(:, 2) > maxAllowed)
    error('%s intervals must be within [%g, %g].', fieldName, minAllowed, maxAllowed);
end

windows = double(rawWindows);
end

function devFrame = local_infer_deviance_frame(nFrames, fixedDevianceFrame)
% LOCAL_INFER_DEVIANCE_FRAME Clamp fixed deviance frame to trial length.
devFrame = max(1, min(nFrames, round(fixedDevianceFrame)));
end

function trialSpliced = local_build_spliced_occluded_deviant(trialNondev, trialDev, devFrame)
% LOCAL_BUILD_SPLICED_OCCLUDED_DEVIANT Lock pre-deviance frames to nondeviant branch.
%
% Rule (same intent as V20):
%   - Frames 1..deviance are copied from nondeviant.
%   - Frames deviance..end come from deviant branch, translated so the
%     splice is exactly position-continuous at deviance.
trialSpliced = trialDev;

xyNondev = double(trialNondev.xy);
xyDeviant = double(trialDev.xy);
if size(xyNondev, 1) ~= size(xyDeviant, 1) || size(xyNondev, 2) ~= 2 || size(xyDeviant, 2) ~= 2
    error('Splice builder expects one-dot trajectories with equal frame count.');
end

nFrames = size(xyDeviant, 1);
spliceFrame = max(1, min(nFrames, round(devFrame)));
xyOut = xyNondev;
xyOut(1:spliceFrame, :) = xyNondev(1:spliceFrame, :);

if spliceFrame < nFrames
    sourcePost = xyDeviant(spliceFrame:end, :);
    shiftVec = xyOut(spliceFrame, :) - sourcePost(1, :);
    xyOut(spliceFrame:end, :) = sourcePost + shiftVec;
end

trialSpliced.xy = xyOut;
if isfield(trialSpliced, 'pathAll')
    trialSpliced.pathAll = trialNondev.pathAll;
    trialSpliced.pathAll(spliceFrame) = 1;
end
end

function timing = local_build_occlusion_timing(nFrames, devFrame, occlusionEndFrameInclusive, fadeInFrames)
% LOCAL_BUILD_OCCLUSION_TIMING Build fixed-frame V23 occlusion event indices.
occlusionStartFrame = max(1, devFrame - 1);
occlusionFullEndFrame = max(devFrame, min(nFrames, round(occlusionEndFrameInclusive)));
reappearanceStartFrame = min(nFrames, occlusionFullEndFrame + 1);
reappearanceEndFrame = min(nFrames, reappearanceStartFrame + max(1, fadeInFrames) - 1);

timing = struct();
timing.deviance_frame = devFrame;
timing.occlusion_start_frame = occlusionStartFrame;
timing.occlusion_complete_frame = devFrame;
timing.occlusion_full_end_frame = occlusionFullEndFrame;
timing.occlusion_end_frame = reappearanceStartFrame;
timing.reappearance_start_frame = reappearanceStartFrame;
timing.reappearance_end_frame = reappearanceEndFrame;
timing.occlusion_end_complete_frame = reappearanceEndFrame;
end

function notwired = local_build_notwired_occluder_metadata( ...
        xy, devFrame, timing, dotRadiusDeg, fixedOcclusionEndFrame, reanchorEnabled)
% LOCAL_BUILD_NOTWIRED_OCCLUDER_METADATA Build geometry for the not-wired circle occluder.
%
% Data flow:
%   one-dot trajectory + fixed occlusion timing -> start occluder geometry
%   solved from exact frame constraints -> optional deviant reanchor geometry.
nFrames = size(xy, 1);
contactFrame = max(1, min(nFrames, round(devFrame)));
hideStartFrame = contactFrame;
hideEndFrame = max(hideStartFrame, min(nFrames, round(timing.occlusion_end_frame) - 1));
resumeFrame = min(nFrames, hideEndFrame + 1);

contactXY = double(xy(contactFrame, 1:2));
startDirection = local_unit_direction_at_frame(xy, contactFrame);
curvatureSign = local_curvature_sign_at_frame(xy, contactFrame);
sideSign = local_side_sign_from_curvature(curvatureSign);

% Deviant-only reanchor:
% keep radius fixed and rotate center placement using direction sampled at
% fixedOcclusionEndFrame - 1 (requested by the user).
reanchorDirectionFrame = max(1, min(nFrames, round(fixedOcclusionEndFrame + ...
    Config_stimuli_generation_V23_NotWiredOccluder.notWiredReanchorDirectionFrameOffset)));
reanchorDirection = local_unit_direction_at_frame(xy, reanchorDirectionFrame);
% Anchor reanchor contact at the last fully hidden frame so the first
% visible frame can start at resumeFrame.
reanchorContactFrame = hideEndFrame;
reanchorContactXY = double(xy(reanchorContactFrame, 1:2));

% Solve one fixed not-wired circle geometry per trial:
%   - fully covers the dot at deviance (frame 130)
%   - never fully covers any frame before deviance
%   - keeps first reappearance frame visible at frame 191
% while staying as close as possible to the speed/curvature target offset.
targetCenterOffset = local_estimate_notwired_center_offset( ...
    xy, contactFrame, startDirection, dotRadiusDeg);
[centerOffset, radiusSolved, resumeConstraintSatisfied] = local_solve_notwired_geometry( ...
    xy, contactFrame, startDirection, sideSign, dotRadiusDeg, targetCenterOffset, ...
    resumeFrame, logical(reanchorEnabled), reanchorDirection, reanchorContactFrame);

startCenterXY = contactXY + sideSign .* startDirection .* centerOffset;
reanchorCenterXY = reanchorContactXY + sideSign .* reanchorDirection .* centerOffset;
usedFallbackSide = false;
occlusionStartFrame = local_find_notwired_occlusion_start_frame( ...
    xy, startCenterXY, radiusSolved, dotRadiusDeg, contactFrame);

notwired = struct();
notwired.start_center_xy = startCenterXY;
notwired.radius = radiusSolved;
notwired.side_sign = sideSign;
notwired.curvature_sign = curvatureSign;
notwired.start_contact_frame = contactFrame;
notwired.start_direction_xy = startDirection;
notwired.hide_start_frame = hideStartFrame;
notwired.hide_end_frame = hideEndFrame;
notwired.resume_frame = resumeFrame;
notwired.reanchor_enabled = logical(reanchorEnabled);
notwired.reanchor_center_xy = reanchorCenterXY;
notwired.reanchor_direction_xy = reanchorDirection;
notwired.reanchor_direction_frame = reanchorDirectionFrame;
notwired.reanchor_contact_frame = reanchorContactFrame;
notwired.used_fallback_side = logical(usedFallbackSide);
notwired.resume_constraint_satisfied = logical(resumeConstraintSatisfied);

% Keep legacy twocircle compatibility fields for downstream readers.
notwired.center_xy = startCenterXY;
notwired.radius_pre = max(1e-6, dotRadiusDeg);
notwired.radius_post = radiusSolved;
notwired.pre_contact_frame = contactFrame;
notwired.post_contact_frame = contactFrame;
notwired.post_deactivate_frame = nFrames;
notwired.occlusion_start_frame = occlusionStartFrame;
notwired.occlusion_complete_frame = timing.occlusion_complete_frame;
notwired.occlusion_end_frame = timing.occlusion_end_frame;
notwired.occlusion_end_complete_frame = timing.occlusion_end_complete_frame;
end

function occlusionStartFrame = local_find_notwired_occlusion_start_frame( ...
        xy, centerXY, radiusDeg, dotRadiusDeg, contactFrame)
% LOCAL_FIND_NOTWIRED_OCCLUSION_START_FRAME Infer first frame with partial/full overlap.
nFrames = size(xy, 1);
contactFrame = max(1, min(nFrames, round(contactFrame)));
occlusionStartFrame = contactFrame;

for iFrame = 1:contactFrame
    d = norm(double(xy(iFrame, 1:2)) - centerXY);
    if radiusDeg >= (d - dotRadiusDeg)
        occlusionStartFrame = iFrame;
        return;
    end
end
end

function centerOffset = local_estimate_notwired_center_offset( ...
        xy, contactFrame, startDirection, dotRadiusDeg)
% LOCAL_ESTIMATE_NOTWIRED_CENTER_OFFSET Estimate center offset from local speed/curvature.
%
% Rationale:
%   Use local path speed and local turning angle to set the center distance
%   along the start direction. This keeps the occluder placement explicitly
%   dependent on local speed and curvature, while capping extreme offsets
%   from near-zero turn angles that would otherwise make the occluder cover
%   almost the whole early path.
nFrames = size(xy, 1);
contactFrame = max(1, min(nFrames, round(contactFrame)));

if contactFrame > 1
    prevXY = double(xy(contactFrame - 1, 1:2));
else
    prevXY = double(xy(contactFrame, 1:2)) - startDirection;
end
currXY = double(xy(contactFrame, 1:2));
localStep = max(1e-8, norm(currXY - prevXY));

if contactFrame < nFrames
    nextDirection = local_unit_direction_at_frame(xy, contactFrame + 1);
else
    nextDirection = startDirection;
end
cosAngle = max(-1, min(1, dot(startDirection, nextDirection)));
turnRad = acos(cosAngle);

if turnRad < 1e-4
    curvatureRadius = 5 * localStep;
else
    curvatureRadius = localStep / turnRad;
end

% Base geometric scale from dot size and local displacement.
baseOffset = max(2 * dotRadiusDeg, 2 * localStep);

% Curvature-derived term with an upper cap to prevent very large circles
% when local turn is close to zero.
curvatureCap = 4 * dotRadiusDeg + 8 * localStep;
curvatureOffset = min(curvatureRadius, curvatureCap);

centerOffset = max(baseOffset, curvatureOffset);
end

function [centerOffset, radiusSolved, resumeConstraintSatisfied] = local_solve_notwired_geometry( ...
        xy, contactFrame, startDirection, sideSign, dotRadiusDeg, targetOffset, ...
        resumeFrame, reanchorEnabled, reanchorDirection, reanchorContactFrame)
% LOCAL_SOLVE_NOTWIRED_GEOMETRY Solve fixed-circle geometry with frame constraints.
%
% Data flow:
%   tangent direction + curvature-sign side + trajectory samples ->
%   1D center-offset search -> fixed circle radius.
%
% Constraints enforced by this solver:
%   1) full occlusion at contact/deviance frame exactly (radius = offset + dot radius)
%   2) no full occlusion before contact frame
%   3) frame resumeFrame is visible (for deviant trials this uses reanchor center).
%
% Selection rule:
%   Choose the feasible offset nearest to the speed/curvature target,
%   preferring offsets not larger than the target when available.
strictEps = 1e-6;
nFrames = size(xy, 1);
contactFrame = max(1, min(nFrames, round(contactFrame)));
resumeFrame = max(1, min(nFrames, round(resumeFrame)));
reanchorContactFrame = max(1, min(nFrames, round(reanchorContactFrame)));
targetOffset = max(1e-4, double(targetOffset));

contactXY = double(xy(contactFrame, 1:2));
resumeXY = double(xy(resumeFrame, 1:2));
reanchorContactXY = double(xy(reanchorContactFrame, 1:2));

if contactFrame > 1
    prevXY = double(xy(contactFrame - 1, 1:2));
else
    prevXY = contactXY - startDirection;
end
localStep = max(1e-5, norm(contactXY - prevXY));

% Keep the lower bound near zero so the trajectory-constrained solver can
% find feasible circles even when the maximal admissible offset is small.
minOffset = 1e-6;
maxOffset = max([targetOffset * 2.5, 12 * dotRadiusDeg, 30 * localStep]);
if maxOffset <= minOffset
    maxOffset = minOffset + max(10 * localStep, dotRadiusDeg);
end

candidateOffsets = linspace(minOffset, maxOffset, 5000)';
nOffsets = numel(candidateOffsets);

startCentersX = contactXY(1) + sideSign .* startDirection(1) .* candidateOffsets;
startCentersY = contactXY(2) + sideSign .* startDirection(2) .* candidateOffsets;

if contactFrame > 1
    preX = double(xy(1:(contactFrame - 1), 1));
    preY = double(xy(1:(contactFrame - 1), 2));
    dPre = sqrt((preX - startCentersX').^2 + (preY - startCentersY').^2);
    marginPreNoFull = min(dPre - candidateOffsets', [], 1)';
    prePartialOverlap = any((dPre - candidateOffsets') <= (2 * dotRadiusDeg), 1)';
else
    marginPreNoFull = inf(nOffsets, 1);
    prePartialOverlap = false(nOffsets, 1);
end

dResumeStart = sqrt((resumeXY(1) - startCentersX).^2 + (resumeXY(2) - startCentersY).^2);
marginResumeStart = dResumeStart - candidateOffsets;

if reanchorEnabled
    reanchorCentersX = reanchorContactXY(1) + sideSign .* reanchorDirection(1) .* candidateOffsets;
    reanchorCentersY = reanchorContactXY(2) + sideSign .* reanchorDirection(2) .* candidateOffsets;
    dResumeReanchor = sqrt((resumeXY(1) - reanchorCentersX).^2 + (resumeXY(2) - reanchorCentersY).^2);
    marginResume = dResumeReanchor - candidateOffsets;
else
    marginResume = marginResumeStart;
end

feasibleMask = (marginPreNoFull > strictEps) & (marginResume > strictEps);
preferredMask = feasibleMask & prePartialOverlap;

if any(preferredMask)
    candidateMask = preferredMask;
elseif any(feasibleMask)
    candidateMask = feasibleMask;
else
    candidateMask = true(nOffsets, 1);
end

offsetSubset = candidateOffsets(candidateMask);
subsetLeTarget = offsetSubset <= targetOffset;
if any(subsetLeTarget)
    centerOffset = offsetSubset(find(subsetLeTarget, 1, 'last'));
else
    [~, minIdx] = min(abs(offsetSubset - targetOffset));
    centerOffset = offsetSubset(minIdx);
end

radiusSolved = centerOffset + dotRadiusDeg;

startCenterXY = contactXY + sideSign .* startDirection .* centerOffset;
resumeMargin = norm(resumeXY - startCenterXY) - centerOffset;
if reanchorEnabled
    reanchorCenterXY = reanchorContactXY + sideSign .* reanchorDirection .* centerOffset;
    resumeMargin = norm(resumeXY - reanchorCenterXY) - centerOffset;
end
resumeConstraintSatisfied = (resumeMargin > strictEps);
end

function unitDirection = local_unit_direction_at_frame(xy, frameIdx)
% LOCAL_UNIT_DIRECTION_AT_FRAME Estimate unit path direction at one frame.
nFrames = size(xy, 1);
frameIdx = max(1, min(nFrames, round(frameIdx)));

if frameIdx > 1
    directionVec = double(xy(frameIdx, 1:2)) - double(xy(frameIdx - 1, 1:2));
elseif frameIdx < nFrames
    directionVec = double(xy(frameIdx + 1, 1:2)) - double(xy(frameIdx, 1:2));
else
    directionVec = [1, 0];
end

if norm(directionVec) <= 1e-12
    if frameIdx < nFrames
        directionVec = double(xy(min(nFrames, frameIdx + 1), 1:2)) - double(xy(frameIdx, 1:2));
    end
    if norm(directionVec) <= 1e-12 && frameIdx > 1
        directionVec = double(xy(frameIdx, 1:2)) - double(xy(frameIdx - 1, 1:2));
    end
end

if norm(directionVec) <= 1e-12
    unitDirection = [1, 0];
else
    unitDirection = directionVec ./ norm(directionVec);
end
end

function curvatureSign = local_curvature_sign_at_frame(xy, frameIdx)
% LOCAL_CURVATURE_SIGN_AT_FRAME Estimate signed local curvature from 3 consecutive frames.
nFrames = size(xy, 1);
frameIdx = max(2, min(nFrames - 1, round(frameIdx)));

if nFrames < 3
    curvatureSign = 1;
    return;
end

prevXY = double(xy(frameIdx - 1, 1:2));
currXY = double(xy(frameIdx, 1:2));
nextXY = double(xy(frameIdx + 1, 1:2));
v1 = currXY - prevXY;
v2 = nextXY - currXY;
crossZ = v1(1) * v2(2) - v1(2) * v2(1);
if abs(crossZ) <= 1e-12
    curvatureSign = 0;
else
    curvatureSign = sign(crossZ);
end
end

function sideSign = local_side_sign_from_curvature(curvatureSign)
% LOCAL_SIDE_SIGN_FROM_CURVATURE Map signed curvature to center-side choice.
if curvatureSign > 0
    sideSign = 1;
elseif curvatureSign < 0
    sideSign = -1;
else
    sideSign = 1;
end
end

function [alphaProfile, geomProfile] = local_build_visible_profiles(nFrames)
% LOCAL_BUILD_VISIBLE_PROFILES Visibility profiles for always-visible condition.
alphaProfile = ones(nFrames, 1);
geomProfile = ones(nFrames, 1);
end

function [alphaProfile, geomProfile] = local_build_alpha_profiles(nFrames, timing)
% LOCAL_BUILD_ALPHA_PROFILES Smooth alpha fallback profile for occluded conditions.
alphaProfile = ones(nFrames, 1);

fadeOutFrames = timing.occlusion_start_frame:timing.deviance_frame;
if ~isempty(fadeOutFrames)
    n = numel(fadeOutFrames);
    phase = (0:(n - 1)) ./ max(n - 1, 1);
    alphaProfile(fadeOutFrames) = 0.5 * (1 + cos(pi * phase));
end

if timing.reappearance_start_frame > (timing.deviance_frame + 1)
    alphaProfile((timing.deviance_frame + 1):(timing.reappearance_start_frame - 1)) = 0;
end

fadeInFrames = timing.reappearance_start_frame:timing.reappearance_end_frame;
if ~isempty(fadeInFrames)
    n = numel(fadeInFrames);
    phase = (0:(n - 1)) ./ max(n - 1, 1);
    alphaProfile(fadeInFrames) = 0.5 * (1 - cos(pi * phase));
end

alphaProfile = max(0, min(1, alphaProfile));
geomProfile = alphaProfile;
end

function outTrial = local_attach_occlusion_fields(inTrial, conditionCode, conditionLabel, ...
        occlusionEnabled, timing, notwired, alphaProfile, geomProfile)
% LOCAL_ATTACH_OCCLUSION_FIELDS Attach condition and occlusion metadata to a trial.
outTrial = inTrial;
outTrial.condition = conditionCode;
outTrial.condition_label = conditionLabel;
outTrial.occlusion_enabled = logical(occlusionEnabled);
outTrial.occlusion_model_default = 'notwired_circle';
outTrial.occlusion_model_options = {'notwired_circle', 'alpha'};

outTrial.deviance_frame = timing.deviance_frame;
outTrial.occlusion_start_frame = notwired.occlusion_start_frame;
outTrial.occlusion_complete_frame = notwired.occlusion_complete_frame;
outTrial.occlusion_end_frame = notwired.occlusion_end_frame;
outTrial.reappearance_start_frame = notwired.occlusion_end_frame;
outTrial.reappearance_end_frame = notwired.occlusion_end_complete_frame;
outTrial.occlusion_end_complete_frame = notwired.occlusion_end_complete_frame;

outTrial.notwired_center_start_xy = notwired.start_center_xy;
outTrial.notwired_radius = notwired.radius;
outTrial.notwired_side_sign = notwired.side_sign;
outTrial.notwired_curvature_sign = notwired.curvature_sign;
outTrial.notwired_direction_start_xy = notwired.start_direction_xy;
outTrial.notwired_start_contact_frame = notwired.start_contact_frame;
outTrial.notwired_hide_start_frame = notwired.hide_start_frame;
outTrial.notwired_hide_end_frame = notwired.hide_end_frame;
outTrial.notwired_resume_frame = notwired.resume_frame;
outTrial.notwired_reanchor_enabled = notwired.reanchor_enabled;
outTrial.notwired_reanchor_center_xy = notwired.reanchor_center_xy;
outTrial.notwired_reanchor_direction_xy = notwired.reanchor_direction_xy;
outTrial.notwired_reanchor_direction_frame = notwired.reanchor_direction_frame;
outTrial.notwired_reanchor_contact_frame = notwired.reanchor_contact_frame;
outTrial.notwired_used_fallback_side = notwired.used_fallback_side;

% Keep legacy twocircle metadata for compatibility with older analysis code.
outTrial.twocircle_center_xy = notwired.center_xy;
outTrial.twocircle_radius_pre = notwired.radius_pre;
outTrial.twocircle_radius_post = notwired.radius_post;
outTrial.twocircle_pre_contact_frame = notwired.pre_contact_frame;
outTrial.twocircle_post_contact_frame = notwired.post_contact_frame;
outTrial.twocircle_post_deactivate_frame = notwired.post_deactivate_frame;

outTrial.visibility_alpha = alphaProfile(:);
outTrial.visibility_geom = geomProfile(:);
end
