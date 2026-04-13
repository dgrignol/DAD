%% Stimuli Generation V24 PathBandOccluder (config-driven one-dot occlusion, no source subject files)
% Script: stimuli_generation_V24_PathBandOccluder.m
% Author: Dami (V24 fixation-collision modes by Codex)
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
%   plus occluded_deviant_predicted in MovDot_SubXX_V24_PathBandOccluder_predicted.mat.
%
% Optional central-fixation collision handling:
%   - configured through Config_stimuli_generation_V24_PathBandOccluder.fixationCollisionMode.
%   - modes:
%       * 'off'   : no special handling.
%       * 'retry' : reject colliding candidates and resample.
%       * 'move'  : minimally translate colliding paired trajectories out of
%                   the fixation exclusion zone when feasible, otherwise retry.
%   - exclusion radius is controlled by
%     Config_stimuli_generation_V24_PathBandOccluder.fixationExclusionRadiusDeg.
%
% Fixed occlusion timeline:
%   - deviance frame fixed at 130 (default, configurable via Config_stimuli_generation_V24_PathBandOccluder)
%   - first fully occluded frame fixed at 130
%   - fully occluded through frame 190 (inclusive)
%   - nominal reappearance search starts at frame 191
%     (actual reappearance metadata is computed geometrically)
%
% Occlusion geometry model (pathband default):
%   - pre band:
%       * follows nondeviant path segment from frame 130 to frame 190
%       * active for frames < 130
%       * used to produce partial occlusion onset based on speed and width
%   - post band:
%       * activates at frame 130
%       * follows active branch segment from frame 130 to frame 190
%         (for deviant trials this switches to deviant branch at deviance)
%       * remains active to trial end so reappearance is geometric
%   - width:
%       * controlled by Config_stimuli_generation_V24_PathBandOccluder.pathBandWidthDeg
%       * constrained to be strictly larger than dot diameter
%   - terminal style:
%       * controlled by Config_stimuli_generation_V24_PathBandOccluder.pathBandTerminalStyle
%       * 'round' keeps rounded terminal caps; 'straight' uses wall-like
%         straight terminals with configurable start backshift
%
% Trigger/event metadata compatibility:
%   The exported trial fields are aligned with:
%     - experiment/MoveDot1_experiment_occlusion_v1.m
%     - experiment/trigger_codes.md
%
% Usage example (interactive from experiment/):
%   addpath('lib');
%   stimuli_generation_V24_PathBandOccluder;
%
% Usage example (set fixed-frame overrides before running):
%   addpath('lib');
%   fixedDevianceFrame = 130;
%   fixedOcclusionEndFrame = 190;
%   overwriteExisting = true;
%   stimuli_generation_V24_PathBandOccluder;
%
% Usage example (force round terminal metadata for comparison):
%   addpath('lib');
%   pathBandTerminalStyle = 'round';
%   pathBandStraightBackshiftDotRadiusScale = 0.0;
%   stimuli_generation_V24_PathBandOccluder;
%
% Usage example (enable fixation-cross avoidance in config and run):
%   addpath('lib');
%   % In Config_stimuli_generation_V24_PathBandOccluder:
%   %   fixationCollisionMode = 'move'; % or 'retry' / 'off'
%   stimuli_generation_V24_PathBandOccluder;
%
% Inputs:
%   - Config_stimuli_generation_V24_PathBandOccluder class on MATLAB path.
%   - targetSubjectID optional workspace override; otherwise prompted.
%
% Outputs:
%   - input_files/MovDot_SubXX_V24_PathBandOccluder.mat
%       with xySeqs containing always_visible / occluded_nondeviant /
%       occluded_deviant and full occlusion metadata fields.
%   - input_files/MovDot_SubXX_V24_PathBandOccluder_predicted.mat
%       with xySeqsPredicted containing occluded_deviant_predicted trials
%       and matching metadata.
%
% Key assumptions:
%   - One-dot only: xy is frames x 2 ([x y]) in visual degrees.
%   - Direction variance contains at least one nondeviant value (0) and one
%     deviant value (>0); the first deviant value is used for generation.
%   - Trial pairing is sequence-index based and one-to-one across conditions.
%   - Output directories are relative to experiment/ unless overridden by
%     Config_stimuli_generation_V24_PathBandOccluder.

addpath('lib/');

%% Setup and runtime parameter defaults
% Data flow: caller overrides + Config_stimuli_generation_V24_PathBandOccluder constants -> validated controls.
clearvars -except overwriteExisting fixedDevianceFrame fixedOcclusionEndFrame targetSubjectID pathBandTerminalStyle pathBandStraightBackshiftDotRadiusScale;
clc;

if ~(exist('fixedDevianceFrame', 'var') && ~isempty(fixedDevianceFrame))
    fixedDevianceFrame = Config_stimuli_generation_V24_PathBandOccluder.fixedDevianceFrame;
end
if ~(exist('fixedOcclusionEndFrame', 'var') && ~isempty(fixedOcclusionEndFrame))
    fixedOcclusionEndFrame = Config_stimuli_generation_V24_PathBandOccluder.fixedOcclusionEndFrame;
end
if ~(exist('overwriteExisting', 'var') && ~isempty(overwriteExisting))
    overwriteExisting = true;
end
if ~(exist('pathBandTerminalStyle', 'var') && ~isempty(pathBandTerminalStyle))
    pathBandTerminalStyle = Config_stimuli_generation_V24_PathBandOccluder.pathBandTerminalStyle;
end
if ~(exist('pathBandStraightBackshiftDotRadiusScale', 'var') && ...
        ~isempty(pathBandStraightBackshiftDotRadiusScale))
    pathBandStraightBackshiftDotRadiusScale = ...
        Config_stimuli_generation_V24_PathBandOccluder.pathBandStraightBackshiftDotRadiusScale;
end

if ~isscalar(fixedDevianceFrame) || ~isscalar(fixedOcclusionEndFrame) || ...
        isnan(fixedDevianceFrame) || isnan(fixedOcclusionEndFrame)
    error('fixedDevianceFrame and fixedOcclusionEndFrame must be numeric scalars.');
end
fixedDevianceFrame = round(fixedDevianceFrame);
fixedOcclusionEndFrame = round(fixedOcclusionEndFrame);
overwriteExisting = logical(overwriteExisting);
pathBandTerminalStyle = lower(strtrim(char(pathBandTerminalStyle)));
pathBandStraightBackshiftDotRadiusScale = double(pathBandStraightBackshiftDotRadiusScale);
if ~ismember(pathBandTerminalStyle, {'round', 'straight'})
    error('pathBandTerminalStyle must be ''round'' or ''straight''.');
end
if ~isfinite(pathBandStraightBackshiftDotRadiusScale) || ...
        pathBandStraightBackshiftDotRadiusScale < 0
    error('pathBandStraightBackshiftDotRadiusScale must be a non-negative finite scalar.');
end

%% Resolve target subject ID (workspace override or interactive prompt)
% Data flow: optional caller-provided subject ID -> fallback dialog -> validated ID.
if ~(exist('targetSubjectID', 'var') && ~isempty(targetSubjectID))
    prompt = {'Target subject ID (new config-driven occlusion dataset):'};
    dlgTitle = 'V24 target subject';
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
% Data flow: targetSubjectID + Config_stimuli_generation_V24_PathBandOccluder directories -> output files.
inputDir = Config_stimuli_generation_V24_PathBandOccluder.inputDirectory;
if ~exist(inputDir, 'dir')
    mkdir(inputDir);
end

targetObservedFile = fullfile(inputDir, sprintf( ...
    Config_stimuli_generation_V24_PathBandOccluder.observedFilePattern, targetSubjectID));
targetPredictedFile = fullfile(inputDir, sprintf( ...
    Config_stimuli_generation_V24_PathBandOccluder.predictedFilePattern, targetSubjectID));

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

%% Build generation controls from Config_stimuli_generation_V24_PathBandOccluder
% Data flow: Config_stimuli_generation_V24_PathBandOccluder constants -> trajectory synthesis controls.
stimulusTypeConfig = Config_stimuli_generation_V24_PathBandOccluder.likelihood;
framesPerTrial = round(Config_stimuli_generation_V24_PathBandOccluder.trialDuration * Config_stimuli_generation_V24_PathBandOccluder.frameFrequency);
nTrials = Config_stimuli_generation_V24_PathBandOccluder.trialsPerCondition;
fps = double(Config_stimuli_generation_V24_PathBandOccluder.frameFrequency);
alphaFadeInFrames = max(1, round(0.25 * fps));
dotRadiusDeg = double(Config_stimuli_generation_V24_PathBandOccluder.dotWidth) / 2;
pathBandWidthDeg = double(Config_stimuli_generation_V24_PathBandOccluder.pathBandWidthDeg);
pathBandStyle = struct( ...
    'terminalStyle', pathBandTerminalStyle, ...
    'straightBackshiftDotRadiusScale', pathBandStraightBackshiftDotRadiusScale);

if nTrials < 1
    error('Config_stimuli_generation_V24_PathBandOccluder.trialsPerCondition must be >= 1.');
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
if ~isfinite(pathBandWidthDeg) || pathBandWidthDeg <= (2 * dotRadiusDeg)
    error(['Config_stimuli_generation_V24_PathBandOccluder.pathBandWidthDeg must be finite and ' ...
        'strictly larger than the dot diameter (%g deg).'], 2 * dotRadiusDeg);
end

directionVariance = double(stimulusTypeConfig.directionVariance(:)');
if ~any(directionVariance == 0)
    error('Config_stimuli_generation_V24_PathBandOccluder.likelihood.directionVariance must include nondeviant value 0.');
end
if ~any(directionVariance > 0)
    error('Config_stimuli_generation_V24_PathBandOccluder.likelihood.directionVariance must include at least one deviant value > 0.');
end
deviantConditionCode = directionVariance(find(directionVariance > 0, 1, 'first'));
if numel(directionVariance(directionVariance > 0)) > 1
    warning('Multiple deviant directionVariance values found; using first deviant value %g.', ...
        deviantConditionCode);
end

if numel(stimulusTypeConfig.pathDuration) > 1
    warning(['Config_stimuli_generation_V24_PathBandOccluder.likelihood.pathDuration has multiple values. ' ...
        'V24 uses the first value for one-dot occlusion generation.']);
end
pathDurationSec = double(stimulusTypeConfig.pathDuration(1));

generationCfg = struct();
generationCfg.framesPerTrial = framesPerTrial;
generationCfg.dotSpeedDegPerFrame = double(Config_stimuli_generation_V24_PathBandOccluder.dotSpeedDegPerFrame);
generationCfg.minPosDeg = [Config_stimuli_generation_V24_PathBandOccluder.dotWidth / 2, Config_stimuli_generation_V24_PathBandOccluder.dotWidth / 2];
generationCfg.maxPosDeg = double(Config_stimuli_generation_V24_PathBandOccluder.dotRectSize) - generationCfg.minPosDeg;
generationCfg.initialCurvatureWindows = local_parse_interval_windows( ...
    Config_stimuli_generation_V24_PathBandOccluder.initialCurvatureWindows, 'Config_stimuli_generation_V24_PathBandOccluder.initialCurvatureWindows', false, -inf, inf);
generationCfg.deviantCurvatureWindows = local_parse_interval_windows( ...
    Config_stimuli_generation_V24_PathBandOccluder.deviantCurvatureWindows, 'Config_stimuli_generation_V24_PathBandOccluder.deviantCurvatureWindows', false, -inf, inf);
generationCfg.deviantTurnWindowsDeg = local_parse_interval_windows( ...
    stimulusTypeConfig.deviantSignedTurnWindows, ...
    'Config_stimuli_generation_V24_PathBandOccluder.likelihood.deviantSignedTurnWindows', true, -180, 180);
generationCfg.flipCurvatureOnDeviant = logical(Config_stimuli_generation_V24_PathBandOccluder.flipCurvatureOnDeviant);
generationCfg.randomizeCurvatureOnDeviant = logical(Config_stimuli_generation_V24_PathBandOccluder.randomizeCurvatureOnDeviant);
generationCfg.pathDurationSec = pathDurationSec;
generationCfg.maxAttemptsPerTrial = 60000;
generationCfg.fixationCollisionMode = lower(strtrim(char( ...
    Config_stimuli_generation_V24_PathBandOccluder.fixationCollisionMode)));
generationCfg.fixationExclusionRadiusDeg = double( ...
    Config_stimuli_generation_V24_PathBandOccluder.fixationExclusionRadiusDeg);
generationCfg.fixationCenterDeg = double(Config_stimuli_generation_V24_PathBandOccluder.dotRectSize(:)') / 2;
generationCfg.fixationMovePaddingDeg = double( ...
    Config_stimuli_generation_V24_PathBandOccluder.fixationMovePaddingDeg);
generationCfg.fixationMoveDirectionSamples = double( ...
    Config_stimuli_generation_V24_PathBandOccluder.fixationMoveDirectionSamples);
generationCfg.fixationMoveShiftSamples = double( ...
    Config_stimuli_generation_V24_PathBandOccluder.fixationMoveShiftSamples);

% Validate fixation-collision controls once, then enforce them per trial.
if numel(generationCfg.fixationCenterDeg) ~= 2
    error('Config_stimuli_generation_V24_PathBandOccluder.dotRectSize must be [width height].');
end
if ~ismember(generationCfg.fixationCollisionMode, {'off', 'retry', 'move'})
    error(['Config_stimuli_generation_V24_PathBandOccluder.fixationCollisionMode must be one of: ' ...
        '''off'', ''retry'', ''move''.']);
end
if ~isfinite(generationCfg.fixationExclusionRadiusDeg) || generationCfg.fixationExclusionRadiusDeg < 0
    error('Config_stimuli_generation_V24_PathBandOccluder.fixationExclusionRadiusDeg must be finite and >= 0.');
end
if ~isfinite(generationCfg.fixationMovePaddingDeg) || generationCfg.fixationMovePaddingDeg < 0
    error('Config_stimuli_generation_V24_PathBandOccluder.fixationMovePaddingDeg must be finite and >= 0.');
end
if generationCfg.fixationMoveDirectionSamples < 4 || ...
        floor(generationCfg.fixationMoveDirectionSamples) ~= generationCfg.fixationMoveDirectionSamples
    error('Config_stimuli_generation_V24_PathBandOccluder.fixationMoveDirectionSamples must be an integer >= 4.');
end
if generationCfg.fixationMoveShiftSamples < 10 || ...
        floor(generationCfg.fixationMoveShiftSamples) ~= generationCfg.fixationMoveShiftSamples
    error('Config_stimuli_generation_V24_PathBandOccluder.fixationMoveShiftSamples must be an integer >= 10.');
end

% Match v17 geometric floor logic to avoid near-straight paths that cannot
% fit inside the movement rectangle under fixed speed/frame constraints.
maxTurnRadiusDeg = min(double(Config_stimuli_generation_V24_PathBandOccluder.dotRectSize) - double(Config_stimuli_generation_V24_PathBandOccluder.dotWidth)) / 2;
if maxTurnRadiusDeg <= 0
    error('Config_stimuli_generation_V24_PathBandOccluder.dotRectSize minus dotWidth must be > 0.');
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

    pathbandNondev = local_build_pathband_metadata( ...
        trialNondev.xy, trialNondev.xy, devFrame, timing, dotRadiusDeg, pathBandWidthDeg, pathBandStyle);
    pathbandDev = local_build_pathband_metadata( ...
        trialNondev.xy, trialDevSpliced.xy, devFrame, timing, dotRadiusDeg, pathBandWidthDeg, pathBandStyle);

    [alphaVisible, geomVisible] = local_build_visible_profiles(framesPerTrial);
    [alphaOcc, geomOcc] = local_build_alpha_profiles(framesPerTrial, timing);

    rowVisible = iTrial;
    rowOccNondev = nTrials + iTrial;
    rowOccDev = 2 * nTrials + iTrial;

    observedCell{rowVisible} = local_attach_occlusion_fields( ...
        trialNondev, -1, 'always_visible', false, ...
        timing, pathbandNondev, alphaVisible, geomVisible);

    observedCell{rowOccNondev} = local_attach_occlusion_fields( ...
        trialNondev, 0, 'occluded_nondeviant', true, ...
        timing, pathbandNondev, alphaOcc, geomOcc);

    observedCell{rowOccDev} = local_attach_occlusion_fields( ...
        trialDevSpliced, deviantConditionCode, 'occluded_deviant', true, ...
        timing, pathbandDev, alphaOcc, geomOcc);

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
    [alphaOcc, geomOcc] = local_build_alpha_profiles(framesPerTrial, timing);

    % Modeling assumption retained from v18-v20:
    % predicted and observed are forced identical until reappearance starts.
    adjustedPred = trialPredictedDev;
    rs = timing.reappearance_start_frame;
    if rs > 1
        adjustedPred.xy(1:(rs - 1), :) = trialObservedDev.xy(1:(rs - 1), :);
    end
    pathbandPred = local_build_pathband_metadata( ...
        adjustedPred.xy, adjustedPred.xy, devFrame, timing, dotRadiusDeg, pathBandWidthDeg, pathBandStyle);

    predictedCell{iTrial} = local_attach_occlusion_fields( ...
        adjustedPred, deviantConditionCode, 'occluded_deviant_predicted', true, ...
        timing, pathbandPred, alphaOcc, geomOcc);
end
xySeqsPredicted = vertcat(predictedCell{:});

%% Build Cfg and reproducibility metadata
% Data flow: Config_stimuli_generation_V24_PathBandOccluder snapshot + runtime values -> Cfg + repro outputs.
Cfg = struct();
Cfg.dpf = Config_stimuli_generation_V24_PathBandOccluder.dotSpeedDegPerFrame;
Cfg.fps = Config_stimuli_generation_V24_PathBandOccluder.frameFrequency;
Cfg.Stimulitype = Config_stimuli_generation_V24_PathBandOccluder.stimulusType;
Cfg.dot_w = Config_stimuli_generation_V24_PathBandOccluder.dotWidth;
Cfg.rectSize = Config_stimuli_generation_V24_PathBandOccluder.dotRectSize;
Cfg.DirChange = stimulusTypeConfig.directionChange;

configProps = properties('Config_stimuli_generation_V24_PathBandOccluder');
configSnapshot = struct();
for iProp = 1:numel(configProps)
    propName = configProps{iProp};
    configSnapshot.(propName) = Config_stimuli_generation_V24_PathBandOccluder.(propName);
end

repro = struct();
repro.script = struct( ...
    'name', mfilename, ...
    'version', 'V24_pathBandOccluder', ...
    'parameters', struct( ...
        'fixedDevianceFrame', fixedDevianceFrame, ...
        'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
        'alphaFadeInFrames', alphaFadeInFrames, ...
        'pathBandWidthDeg', pathBandWidthDeg, ...
        'pathBandTerminalStyle', pathBandStyle.terminalStyle, ...
        'pathBandStraightBackshiftDotRadiusScale', pathBandStyle.straightBackshiftDotRadiusScale, ...
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
repro.v24_pathband_occlusion = struct( ...
    'script', mfilename, ...
    'targetSubjectID', targetSubjectID, ...
    'fixedDevianceFrame', fixedDevianceFrame, ...
    'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
    'nominalReappearanceFrame', fixedOcclusionEndFrame + 1, ...
    'dotRadiusDeg', dotRadiusDeg, ...
    'pathBandWidthDeg', pathBandWidthDeg, ...
    'pathBandTerminalStyle', pathBandStyle.terminalStyle, ...
    'pathBandStraightBackshiftDotRadiusScale', pathBandStyle.straightBackshiftDotRadiusScale, ...
    'deviantConstructionRule', 'config_generated_with_nondeviant_prefix_locked_at_deviance', ...
    'pathbandRule', ['pre_band_follows_nondeviant_trajectory_from_deviance_to_full_occlusion_end;' ...
        'post_band_follows_active_branch_trajectory_from_deviance_to_full_occlusion_end;' ...
        'terminal_style_controls_start_end_caps_and_start_backshift;' ...
        'band_width_strictly_larger_than_dot_diameter'], ...
    'defaultModel', 'pathband', ...
    'optionalModel', 'alpha', ...
    'conditionCodes', struct('always_visible', -1, 'occluded_nondeviant', 0, ...
        'occluded_deviant', deviantConditionCode));

%% Save outputs
% Data flow: generated structs -> MAT files in Config_stimuli_generation_V24_PathBandOccluder.inputDirectory.
save(targetObservedFile, 'xySeqs', 'Cfg', 'repro');
save(targetPredictedFile, 'xySeqsPredicted', 'Cfg');

fprintf('Saved observed occlusion dataset: %s\n', targetObservedFile);
fprintf('Saved predicted occlusion dataset: %s\n', targetPredictedFile);
fprintf('Trials per condition: %d (always_visible / occluded_nondeviant / occluded_deviant).\n', nTrials);
fprintf(['Fixed occlusion timing: deviance=%d, first full occlusion=%d, ' ...
    'last required full occlusion=%d, nominal reappearance search start=%d.\n'], ...
    fixedDevianceFrame, fixedDevianceFrame, fixedOcclusionEndFrame, fixedOcclusionEndFrame + 1);

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
    'Check Config_stimuli_generation_V24_PathBandOccluder geometry and curvature settings.'], ...
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
% LOCAL_BUILD_OCCLUSION_TIMING Build fixed-frame V24 occlusion event indices.
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

function pathband = local_build_pathband_metadata( ...
        xyPre, xyPost, devFrame, timing, dotRadiusDeg, pathBandWidthDeg, pathBandStyle)
% LOCAL_BUILD_PATHBAND_METADATA Compute path-band geometry and occlusion events.
%
% Data flow:
%   pre-branch trajectory + post-branch trajectory + band width ->
%   framewise overlap tests -> occlusion event frame metadata.
%
% Geometry rule:
%   - Pre-deviance band follows the nondeviant trajectory segment from
%     deviance to fixed full-occlusion end (used only for frames < deviance).
%   - Post-deviance band follows the active branch (nondeviant or deviant)
%     over the same segment frame window.
%   - Band width is strictly greater than dot diameter.
%
% Timing rule enforced:
%   - First fully occluded frame is fixed to deviance frame by disallowing
%     full invisibility for frames before deviance.
if size(xyPre, 2) ~= 2 || size(xyPost, 2) ~= 2 || size(xyPre, 1) ~= size(xyPost, 1)
    error('Path-band metadata expects paired one-dot trajectories with identical frame count.');
end

nFrames = size(xyPost, 1);
devFrame = max(1, min(nFrames, round(devFrame)));
holdEndFrame = max(devFrame, min(nFrames, timing.occlusion_full_end_frame));

prePolyline = local_build_pathband_polyline(xyPre, devFrame, holdEndFrame);
postPolyline = local_build_pathband_polyline(xyPost, devFrame, holdEndFrame);

% Resolve style options and apply geometry transforms (e.g., straight-style
% start backshift) before evaluating visibility events.
styleOpts = local_resolve_pathband_style_options(pathBandStyle, dotRadiusDeg);
prePolyline = local_apply_polyline_style(prePolyline, styleOpts);
postPolyline = local_apply_polyline_style(postPolyline, styleOpts);

halfWidthDeg = max(1e-6, pathBandWidthDeg / 2);
fullInvisibleDistanceDeg = max(0, halfWidthDeg - dotRadiusDeg);
fullVisibleDistanceDeg = halfWidthDeg + dotRadiusDeg;

occlusionStartFrame = NaN;
occlusionCompleteFrame = NaN;
for iFrame = 1:nFrames
    dotCenter = double(xyPost(iFrame, 1:2));
    if iFrame < devFrame
        activePolyline = prePolyline;
        allowFullInvisibility = false;
    else
        activePolyline = postPolyline;
        allowFullInvisibility = true;
    end
    [isInvisible, isFullyVisible] = local_pathband_visibility_state( ...
        dotCenter, activePolyline, fullInvisibleDistanceDeg, ...
        fullVisibleDistanceDeg, allowFullInvisibility, styleOpts);
    if isnan(occlusionStartFrame) && ~isFullyVisible
        occlusionStartFrame = iFrame;
    end
    if isInvisible
        occlusionCompleteFrame = iFrame;
        break;
    end
end
if isnan(occlusionStartFrame)
    occlusionStartFrame = max(1, devFrame - 1);
end
if isnan(occlusionCompleteFrame)
    occlusionCompleteFrame = devFrame;
end

% Required behavior: first complete occlusion is the deviance frame.
occlusionCompleteFrame = max(devFrame, occlusionCompleteFrame);

occlusionEndFrame = nFrames;
searchStart = min(nFrames, holdEndFrame + 1);
for iFrame = searchStart:nFrames
    dotCenter = double(xyPost(iFrame, 1:2));
    [isInvisible, ~] = local_pathband_visibility_state( ...
        dotCenter, postPolyline, fullInvisibleDistanceDeg, ...
        fullVisibleDistanceDeg, true, styleOpts);
    if ~isInvisible
        occlusionEndFrame = iFrame;
        break;
    end
end

occlusionEndCompleteFrame = timing.reappearance_end_frame;
for iFrame = max(occlusionEndFrame, searchStart):nFrames
    dotCenter = double(xyPost(iFrame, 1:2));
    [~, isFullyVisible] = local_pathband_visibility_state( ...
        dotCenter, postPolyline, fullInvisibleDistanceDeg, ...
        fullVisibleDistanceDeg, true, styleOpts);
    if isFullyVisible
        occlusionEndCompleteFrame = iFrame;
        break;
    end
end

pathband = struct();
pathband.pre_xy = prePolyline;
pathband.post_xy = postPolyline;
pathband.width_deg = pathBandWidthDeg;
pathband.half_width_deg = halfWidthDeg;
pathband.pre_anchor_frame = devFrame;
pathband.post_anchor_frame = devFrame;
pathband.post_deactivate_frame = nFrames;
pathband.occlusion_start_frame = occlusionStartFrame;
pathband.occlusion_complete_frame = occlusionCompleteFrame;
pathband.occlusion_end_frame = occlusionEndFrame;
pathband.occlusion_end_complete_frame = occlusionEndCompleteFrame;
pathband.terminal_style = styleOpts.terminalStyle;
pathband.straight_backshift_dot_radius_scale = styleOpts.straightBackshiftDotRadiusScale;
end

function polyline = local_build_pathband_polyline(xy, startFrame, endFrame)
% LOCAL_BUILD_PATHBAND_POLYLINE Build a finite polyline segment for band occlusion.
nFrames = size(xy, 1);
startFrame = max(1, min(nFrames, round(startFrame)));
endFrame = max(startFrame, min(nFrames, round(endFrame)));
polyline = double(xy(startFrame:endFrame, 1:2));

finiteMask = isfinite(polyline(:, 1)) & isfinite(polyline(:, 2));
polyline = polyline(finiteMask, :);
if isempty(polyline)
    polyline = double(xy(startFrame, 1:2));
end

if size(polyline, 1) > 1
    segmentLen = sqrt(sum(diff(polyline, 1, 1) .^ 2, 2));
    keepMask = [true; (segmentLen > 1e-10)];
    polyline = polyline(keepMask, :);
end

if size(polyline, 1) < 2
    if endFrame < nFrames
        fallbackPoint = double(xy(endFrame + 1, 1:2));
        if any(~isfinite(fallbackPoint))
            fallbackPoint = polyline(1, :) + [1e-6, 0];
        end
    else
        fallbackPoint = polyline(1, :) + [1e-6, 0];
    end
    polyline = [polyline; fallbackPoint];
end
end

function [isInvisible, isFullyVisible] = local_pathband_visibility_state( ...
        dotCenter, polyline, fullInvisibleDistanceDeg, fullVisibleDistanceDeg, ...
        allowFullInvisibility, styleOpts)
% LOCAL_PATHBAND_VISIBILITY_STATE Determine full-invisible/full-visible states.
distanceToBandCenter = local_distance_point_to_polyline_style(dotCenter, polyline, styleOpts);
isInvisible = logical(allowFullInvisibility) && ...
    (distanceToBandCenter <= (fullInvisibleDistanceDeg + 1e-12));
isFullyVisible = distanceToBandCenter >= (fullVisibleDistanceDeg - 1e-12);
end

function styleOpts = local_resolve_pathband_style_options(pathBandStyle, dotRadiusDeg)
% LOCAL_RESOLVE_PATHBAND_STYLE_OPTIONS Normalize style controls for geometry checks.
styleOpts = struct();
styleOpts.terminalStyle = 'round';
styleOpts.drawStartCap = true;
styleOpts.drawEndCap = true;
styleOpts.startBackshiftDeg = 0;
styleOpts.straightBackshiftDotRadiusScale = 0;

if nargin < 1 || isempty(pathBandStyle)
    return;
end

if isstruct(pathBandStyle) && isfield(pathBandStyle, 'terminalStyle')
    terminalStyle = lower(strtrim(char(pathBandStyle.terminalStyle)));
else
    terminalStyle = lower(strtrim(char(pathBandStyle)));
end
if ~ismember(terminalStyle, {'round', 'straight'})
    error('Unsupported pathBand terminal style: %s', terminalStyle);
end

styleOpts.terminalStyle = terminalStyle;
if strcmp(terminalStyle, 'straight')
    backshiftScale = 0.0;
    if isstruct(pathBandStyle) && isfield(pathBandStyle, 'straightBackshiftDotRadiusScale')
        backshiftScale = double(pathBandStyle.straightBackshiftDotRadiusScale);
    end
    if ~isfinite(backshiftScale) || backshiftScale < 0
        backshiftScale = 0;
    end
    styleOpts.drawStartCap = false;
    styleOpts.drawEndCap = false;
    styleOpts.straightBackshiftDotRadiusScale = backshiftScale;
    styleOpts.startBackshiftDeg = max(0, backshiftScale * double(dotRadiusDeg));
end
end

function polylineOut = local_apply_polyline_style(polylineIn, styleOpts)
% LOCAL_APPLY_POLYLINE_STYLE Apply style-specific transforms to polyline geometry.
polylineOut = double(polylineIn);
if isempty(polylineOut) || size(polylineOut, 1) < 2
    return;
end
polylineOut = local_apply_start_backshift_to_polyline(polylineOut, styleOpts.startBackshiftDeg);
end

function polylineOut = local_apply_start_backshift_to_polyline(polylineIn, backshiftDeg)
% LOCAL_APPLY_START_BACKSHIFT_TO_POLYLINE Extend first segment backward by degrees.
polylineOut = double(polylineIn);
if backshiftDeg <= 0 || isempty(polylineOut) || size(polylineOut, 1) < 2
    return;
end

p1 = polylineOut(1, :);
for iSeg = 1:(size(polylineOut, 1) - 1)
    p2 = polylineOut(iSeg + 1, :);
    if ~(all(isfinite(p1)) && all(isfinite(p2)))
        continue;
    end
    v = p2 - p1;
    vNorm = norm(v);
    if vNorm <= 1e-12
        continue;
    end
    u = v ./ vNorm;
    polylineOut(1, :) = p1 - backshiftDeg .* u;
    return;
end
end

function minDistance = local_distance_point_to_polyline_style(pointXY, polylineXY, styleOpts)
% LOCAL_DISTANCE_POINT_TO_POLYLINE_STYLE Distance to style-aware polyline support.
%
% Terminal behavior:
%   - start terminal uses rounded cap only if styleOpts.drawStartCap=true.
%   - end terminal uses rounded cap only if styleOpts.drawEndCap=true.
% Internal joins are always connected through neighboring segments.
if size(polylineXY, 1) == 1
    if logical(styleOpts.drawStartCap) || logical(styleOpts.drawEndCap)
        minDistance = norm(pointXY - polylineXY(1, :));
    else
        minDistance = inf;
    end
    return;
end

nSeg = size(polylineXY, 1) - 1;
minDistance = inf;
for iSeg = 1:nSeg
    segStart = polylineXY(iSeg, :);
    segEnd = polylineXY(iSeg + 1, :);
    segVec = segEnd - segStart;
    segLen2 = sum(segVec .^ 2);
    if segLen2 <= 1e-12
        continue;
    end

    ptVec = pointXY - segStart;
    t = sum(ptVec .* segVec) / segLen2;
    if t >= 0 && t <= 1
        closestPoint = segStart + t .* segVec;
        d = norm(pointXY - closestPoint);
        if d < minDistance
            minDistance = d;
        end
    elseif t < 0
        allowEndpoint = (iSeg > 1) || logical(styleOpts.drawStartCap);
        if allowEndpoint
            d = norm(pointXY - segStart);
            if d < minDistance
                minDistance = d;
            end
        end
    else
        allowEndpoint = (iSeg < nSeg) || logical(styleOpts.drawEndCap);
        if allowEndpoint
            d = norm(pointXY - segEnd);
            if d < minDistance
                minDistance = d;
            end
        end
    end
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
        occlusionEnabled, timing, pathband, alphaProfile, geomProfile)
% LOCAL_ATTACH_OCCLUSION_FIELDS Attach condition and occlusion metadata to a trial.
outTrial = inTrial;
outTrial.condition = conditionCode;
outTrial.condition_label = conditionLabel;
outTrial.occlusion_enabled = logical(occlusionEnabled);
outTrial.occlusion_model_default = 'pathband';
outTrial.occlusion_model_options = {'pathband', 'alpha'};

outTrial.deviance_frame = timing.deviance_frame;
outTrial.occlusion_start_frame = pathband.occlusion_start_frame;
outTrial.occlusion_complete_frame = pathband.occlusion_complete_frame;
outTrial.occlusion_end_frame = pathband.occlusion_end_frame;
outTrial.reappearance_start_frame = pathband.occlusion_end_frame;
outTrial.reappearance_end_frame = pathband.occlusion_end_complete_frame;
outTrial.occlusion_end_complete_frame = pathband.occlusion_end_complete_frame;

% Path-band occluder fields used by v14 runtime.
outTrial.pathband_pre_xy = pathband.pre_xy;
outTrial.pathband_post_xy = pathband.post_xy;
outTrial.pathband_width_deg = pathband.width_deg;
outTrial.pathband_half_width_deg = pathband.half_width_deg;
outTrial.pathband_pre_anchor_frame = pathband.pre_anchor_frame;
outTrial.pathband_post_anchor_frame = pathband.post_anchor_frame;
outTrial.pathband_post_deactivate_frame = pathband.post_deactivate_frame;
outTrial.pathband_terminal_style = pathband.terminal_style;
outTrial.pathband_straight_backshift_dot_radius_scale = ...
    pathband.straight_backshift_dot_radius_scale;

outTrial.visibility_alpha = alphaProfile(:);
outTrial.visibility_geom = geomProfile(:);
end
