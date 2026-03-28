%% Stimuli Generation V21 (config-driven one-dot occlusion, no source subject files)
% Script: stimuli_generation_V21.m
% Author: Dami (V21 sandbox prototype by Codex)
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
% Fixed occlusion timeline:
%   - deviance frame fixed at 130 (default, configurable via Config_stimuli_generation_V21)
%   - first fully occluded frame fixed at 130
%   - fully occluded through frame 190 (inclusive)
%   - nominal reappearance search starts at frame 191
%     (actual reappearance metadata is computed geometrically)
%
% Occlusion geometry model (twocircle default):
%   - pre occluder:
%       * center = dot position at frame 130
%       * radius = moving-dot radius
%       * active for frames < 130
%   - post occluder:
%       * activates at frame 130
%       * center anchored at frame-130 position
%       * radius chosen to fully occlude all centers in 130..190
%       * remains active to trial end so reappearance is geometric
%
% Trigger/event metadata compatibility:
%   The exported trial fields are aligned with:
%     - experiment/MoveDot1_experiment_occlusion_v1.m
%     - experiment/trigger_codes.md
%
% Usage example (interactive from experiment/):
%   addpath('lib');
%   stimuli_generation_V21;
%
% Usage example (set fixed-frame overrides before running):
%   addpath('lib');
%   fixedDevianceFrame = 130;
%   fixedOcclusionEndFrame = 190;
%   overwriteExisting = true;
%   stimuli_generation_V21;
%
% Inputs:
%   - Config_stimuli_generation_V21 class on MATLAB path.
%   - Interactive target subject ID (always prompted, no default target).
%
% Outputs:
%   - input_files/MovDot_SubXX.mat
%       with xySeqs containing always_visible / occluded_nondeviant /
%       occluded_deviant and full occlusion metadata fields.
%   - input_files/MovDot_SubXX_predicted.mat
%       with xySeqsPredicted containing occluded_deviant_predicted trials
%       and matching metadata.
%
% Key assumptions:
%   - One-dot only: xy is frames x 2 ([x y]) in visual degrees.
%   - Direction variance contains at least one nondeviant value (0) and one
%     deviant value (>0); the first deviant value is used for generation.
%   - Trial pairing is sequence-index based and one-to-one across conditions.
%   - Output directories are relative to experiment/ unless overridden by
%     Config_stimuli_generation_V21.

addpath('lib/');

%% Setup and runtime parameter defaults
% Data flow: caller overrides + Config_stimuli_generation_V21 constants -> validated controls.
clearvars -except overwriteExisting fixedDevianceFrame fixedOcclusionEndFrame;
clc;

if ~(exist('fixedDevianceFrame', 'var') && ~isempty(fixedDevianceFrame))
    fixedDevianceFrame = Config_stimuli_generation_V21.fixedDevianceFrame;
end
if ~(exist('fixedOcclusionEndFrame', 'var') && ~isempty(fixedOcclusionEndFrame))
    fixedOcclusionEndFrame = Config_stimuli_generation_V21.fixedOcclusionEndFrame;
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

%% Interactive target subject prompt (always requested)
% Data flow: dialog input -> targetSubjectID -> RNG seed + output filenames.
prompt = {'Target subject ID (new config-driven occlusion dataset):'};
dlgTitle = 'V21 target subject';
dims = [1 70];
definput = {''}; % no default target by design.
userInput = inputdlg(prompt, dlgTitle, dims, definput);
if isempty(userInput)
    disp('Canceled by user.');
    return;
end
targetSubjectID = str2double(userInput{1});

if isnan(targetSubjectID) || targetSubjectID < 0
    error('Target subject ID must be a numeric value >= 0.');
end
targetSubjectID = round(targetSubjectID);

%% Resolve output paths and overwrite guard
% Data flow: targetSubjectID + Config_stimuli_generation_V21 directories -> output files.
inputDir = Config_stimuli_generation_V21.inputDirectory;
if ~exist(inputDir, 'dir')
    mkdir(inputDir);
end

targetObservedFile = fullfile(inputDir, sprintf('MovDot_Sub%02d.mat', targetSubjectID));
targetPredictedFile = fullfile(inputDir, sprintf('MovDot_Sub%02d_predicted.mat', targetSubjectID));

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

%% Build generation controls from Config_stimuli_generation_V21
% Data flow: Config_stimuli_generation_V21 constants -> trajectory synthesis controls.
stimulusTypeConfig = Config_stimuli_generation_V21.likelihood;
framesPerTrial = round(Config_stimuli_generation_V21.trialDuration * Config_stimuli_generation_V21.frameFrequency);
nTrials = Config_stimuli_generation_V21.trialsPerCondition;
fps = double(Config_stimuli_generation_V21.frameFrequency);
alphaFadeInFrames = max(1, round(0.25 * fps));
dotRadiusDeg = double(Config_stimuli_generation_V21.dotWidth) / 2;

if nTrials < 1
    error('Config_stimuli_generation_V21.trialsPerCondition must be >= 1.');
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
    error('Config_stimuli_generation_V21.likelihood.directionVariance must include nondeviant value 0.');
end
if ~any(directionVariance > 0)
    error('Config_stimuli_generation_V21.likelihood.directionVariance must include at least one deviant value > 0.');
end
deviantConditionCode = directionVariance(find(directionVariance > 0, 1, 'first'));
if numel(directionVariance(directionVariance > 0)) > 1
    warning('Multiple deviant directionVariance values found; using first deviant value %g.', ...
        deviantConditionCode);
end

if numel(stimulusTypeConfig.pathDuration) > 1
    warning(['Config_stimuli_generation_V21.likelihood.pathDuration has multiple values. ' ...
        'V21 uses the first value for one-dot occlusion generation.']);
end
pathDurationSec = double(stimulusTypeConfig.pathDuration(1));

generationCfg = struct();
generationCfg.framesPerTrial = framesPerTrial;
generationCfg.dotSpeedDegPerFrame = double(Config_stimuli_generation_V21.dotSpeedDegPerFrame);
generationCfg.minPosDeg = [Config_stimuli_generation_V21.dotWidth / 2, Config_stimuli_generation_V21.dotWidth / 2];
generationCfg.maxPosDeg = double(Config_stimuli_generation_V21.dotRectSize) - generationCfg.minPosDeg;
generationCfg.initialCurvatureWindows = local_parse_interval_windows( ...
    Config_stimuli_generation_V21.initialCurvatureWindows, 'Config_stimuli_generation_V21.initialCurvatureWindows', false, -inf, inf);
generationCfg.deviantCurvatureWindows = local_parse_interval_windows( ...
    Config_stimuli_generation_V21.deviantCurvatureWindows, 'Config_stimuli_generation_V21.deviantCurvatureWindows', false, -inf, inf);
generationCfg.deviantTurnWindowsDeg = local_parse_interval_windows( ...
    stimulusTypeConfig.deviantSignedTurnWindows, ...
    'Config_stimuli_generation_V21.likelihood.deviantSignedTurnWindows', true, -180, 180);
generationCfg.flipCurvatureOnDeviant = logical(Config_stimuli_generation_V21.flipCurvatureOnDeviant);
generationCfg.randomizeCurvatureOnDeviant = logical(Config_stimuli_generation_V21.randomizeCurvatureOnDeviant);
generationCfg.pathDurationSec = pathDurationSec;
generationCfg.maxAttemptsPerTrial = 60000;

% Match v17 geometric floor logic to avoid near-straight paths that cannot
% fit inside the movement rectangle under fixed speed/frame constraints.
maxTurnRadiusDeg = min(double(Config_stimuli_generation_V21.dotRectSize) - double(Config_stimuli_generation_V21.dotWidth)) / 2;
if maxTurnRadiusDeg <= 0
    error('Config_stimuli_generation_V21.dotRectSize minus dotWidth must be > 0.');
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

    twocircleNondev = local_build_twocircle_metadata(trialNondev.xy, devFrame, timing, dotRadiusDeg);
    twocircleDev = local_build_twocircle_metadata(trialDevSpliced.xy, devFrame, timing, dotRadiusDeg);

    [alphaVisible, geomVisible] = local_build_visible_profiles(framesPerTrial);
    [alphaOcc, geomOcc] = local_build_alpha_profiles(framesPerTrial, timing);

    rowVisible = iTrial;
    rowOccNondev = nTrials + iTrial;
    rowOccDev = 2 * nTrials + iTrial;

    observedCell{rowVisible} = local_attach_occlusion_fields( ...
        trialNondev, -1, 'always_visible', false, ...
        timing, twocircleNondev, alphaVisible, geomVisible);

    observedCell{rowOccNondev} = local_attach_occlusion_fields( ...
        trialNondev, 0, 'occluded_nondeviant', true, ...
        timing, twocircleNondev, alphaOcc, geomOcc);

    observedCell{rowOccDev} = local_attach_occlusion_fields( ...
        trialDevSpliced, deviantConditionCode, 'occluded_deviant', true, ...
        timing, twocircleDev, alphaOcc, geomOcc);

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
    twocirclePred = local_build_twocircle_metadata(trialPredictedDev.xy, devFrame, timing, dotRadiusDeg);
    [alphaOcc, geomOcc] = local_build_alpha_profiles(framesPerTrial, timing);

    % Modeling assumption retained from v18-v20:
    % predicted and observed are forced identical until reappearance starts.
    adjustedPred = trialPredictedDev;
    rs = timing.reappearance_start_frame;
    if rs > 1
        adjustedPred.xy(1:(rs - 1), :) = trialObservedDev.xy(1:(rs - 1), :);
    end

    predictedCell{iTrial} = local_attach_occlusion_fields( ...
        adjustedPred, deviantConditionCode, 'occluded_deviant_predicted', true, ...
        timing, twocirclePred, alphaOcc, geomOcc);
end
xySeqsPredicted = vertcat(predictedCell{:});

%% Build Cfg and reproducibility metadata
% Data flow: Config_stimuli_generation_V21 snapshot + runtime values -> Cfg + repro outputs.
Cfg = struct();
Cfg.dpf = Config_stimuli_generation_V21.dotSpeedDegPerFrame;
Cfg.fps = Config_stimuli_generation_V21.frameFrequency;
Cfg.Stimulitype = Config_stimuli_generation_V21.stimulusType;
Cfg.dot_w = Config_stimuli_generation_V21.dotWidth;
Cfg.rectSize = Config_stimuli_generation_V21.dotRectSize;
Cfg.DirChange = stimulusTypeConfig.directionChange;

configProps = properties('Config_stimuli_generation_V21');
configSnapshot = struct();
for iProp = 1:numel(configProps)
    propName = configProps{iProp};
    configSnapshot.(propName) = Config_stimuli_generation_V21.(propName);
end

repro = struct();
repro.script = struct( ...
    'name', mfilename, ...
    'version', 'V21_configDrivenOcclusion', ...
    'parameters', struct( ...
        'fixedDevianceFrame', fixedDevianceFrame, ...
        'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
        'alphaFadeInFrames', alphaFadeInFrames, ...
        'deviantConditionCode', deviantConditionCode, ...
        'pathDurationSec', pathDurationSec, ...
        'minimumAbsCurvatureDeg', generationCfg.minimumAbsCurvatureDeg, ...
        'maxAttemptsPerTrial', generationCfg.maxAttemptsPerTrial));
repro.inputs = struct('targetSubjectID', targetSubjectID);
repro.rng = rngState;
repro.config = configSnapshot;
repro.v21_occlusion = struct( ...
    'script', mfilename, ...
    'targetSubjectID', targetSubjectID, ...
    'fixedDevianceFrame', fixedDevianceFrame, ...
    'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
    'nominalReappearanceFrame', fixedOcclusionEndFrame + 1, ...
    'dotRadiusDeg', dotRadiusDeg, ...
    'deviantConstructionRule', 'config_generated_with_nondeviant_prefix_locked_at_deviance', ...
    'twocircleRule', 'pre_radius_equals_dot_radius;post_radius_covers_frames_130_to_190;post_active_until_trial_end', ...
    'defaultModel', 'twocircle', ...
    'optionalModel', 'alpha', ...
    'conditionCodes', struct('always_visible', -1, 'occluded_nondeviant', 0, ...
        'occluded_deviant', deviantConditionCode));

%% Save outputs
% Data flow: generated structs -> MAT files in Config_stimuli_generation_V21.inputDirectory.
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
    'Check Config_stimuli_generation_V21 geometry and curvature settings.'], ...
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
% LOCAL_BUILD_OCCLUSION_TIMING Build fixed-frame V21 occlusion event indices.
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

function twocircle = local_build_twocircle_metadata(xy, devFrame, timing, dotRadiusDeg)
% LOCAL_BUILD_TWOCIRCLE_METADATA Compute twocircle geometry and event frames.
%
% V21 fixed-frame geometry:
%   - pre circle: radius == dot radius at frame-130 center.
%   - post circle: radius covers every center in frames 130..190.
nFrames = size(xy, 1);
center = double(xy(devFrame, 1:2));
preFrame = devFrame;
postFrame = devFrame;
holdEndFrame = max(devFrame, min(nFrames, timing.occlusion_full_end_frame));

rPre = max(1e-6, dotRadiusDeg);

windowXY = double(xy(devFrame:holdEndFrame, 1:2));
distWindow = sqrt(sum((windowXY - center) .^ 2, 2));
rPost = max(1e-6, max(distWindow) + dotRadiusDeg);

preActive = false(nFrames, 1);
preActive(1:max(1, devFrame - 1)) = true;

postActive = false(nFrames, 1);
postActive(devFrame:nFrames) = true;
postDeactivateFrame = nFrames;

occlusionStartFrame = NaN;
occlusionCompleteFrame = NaN;
for iFrame = 1:nFrames
    [isInvisible, isFullyVisible] = local_twocircle_visibility_state( ...
        double(xy(iFrame, 1:2)), center, rPre, rPost, dotRadiusDeg, ...
        preActive(iFrame), postActive(iFrame));
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

occlusionEndFrame = nFrames;
searchStart = min(nFrames, holdEndFrame + 1);
for iFrame = searchStart:nFrames
    [isInvisible, ~] = local_twocircle_visibility_state( ...
        double(xy(iFrame, 1:2)), center, rPre, rPost, dotRadiusDeg, ...
        preActive(iFrame), postActive(iFrame));
    if ~isInvisible
        occlusionEndFrame = iFrame;
        break;
    end
end

occlusionEndCompleteFrame = timing.reappearance_end_frame;
for iFrame = max(occlusionEndFrame, searchStart):nFrames
    [~, isFullyVisible] = local_twocircle_visibility_state( ...
        double(xy(iFrame, 1:2)), center, rPre, rPost, dotRadiusDeg, ...
        preActive(iFrame), postActive(iFrame));
    if isFullyVisible
        occlusionEndCompleteFrame = iFrame;
        break;
    end
end

twocircle = struct();
twocircle.center_xy = center;
twocircle.radius_pre = rPre;
twocircle.radius_post = rPost;
twocircle.pre_contact_frame = preFrame;
twocircle.post_contact_frame = postFrame;
twocircle.post_deactivate_frame = postDeactivateFrame;
twocircle.occlusion_start_frame = occlusionStartFrame;
twocircle.occlusion_complete_frame = occlusionCompleteFrame;
twocircle.occlusion_end_frame = occlusionEndFrame;
twocircle.occlusion_end_complete_frame = occlusionEndCompleteFrame;
end

function [isInvisible, isFullyVisible] = local_twocircle_visibility_state( ...
        dotCenter, circleCenter, rPre, rPost, dotRadiusDeg, preActive, postActive)
% LOCAL_TWOCIRCLE_VISIBILITY_STATE Determine fully invisible / fully visible states.
activeRadii = [];
if preActive
    activeRadii(end + 1) = rPre;
end
if postActive
    activeRadii(end + 1) = rPost;
end

if isempty(activeRadii)
    isInvisible = false;
    isFullyVisible = true;
    return;
end

distCenter = norm(dotCenter - circleCenter);
isInvisible = any(distCenter <= (activeRadii - dotRadiusDeg));
isFullyVisible = all(distCenter >= (activeRadii + dotRadiusDeg));
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
        occlusionEnabled, timing, twocircle, alphaProfile, geomProfile)
% LOCAL_ATTACH_OCCLUSION_FIELDS Attach condition and occlusion metadata to a trial.
outTrial = inTrial;
outTrial.condition = conditionCode;
outTrial.condition_label = conditionLabel;
outTrial.occlusion_enabled = logical(occlusionEnabled);
outTrial.occlusion_model_default = 'twocircle';
outTrial.occlusion_model_options = {'twocircle', 'alpha'};

outTrial.deviance_frame = timing.deviance_frame;
outTrial.occlusion_start_frame = twocircle.occlusion_start_frame;
outTrial.occlusion_complete_frame = twocircle.occlusion_complete_frame;
outTrial.occlusion_end_frame = twocircle.occlusion_end_frame;
outTrial.reappearance_start_frame = twocircle.occlusion_end_frame;
outTrial.reappearance_end_frame = twocircle.occlusion_end_complete_frame;
outTrial.occlusion_end_complete_frame = twocircle.occlusion_end_complete_frame;

outTrial.twocircle_center_xy = twocircle.center_xy;
outTrial.twocircle_radius_pre = twocircle.radius_pre;
outTrial.twocircle_radius_post = twocircle.radius_post;
outTrial.twocircle_pre_contact_frame = twocircle.pre_contact_frame;
outTrial.twocircle_post_contact_frame = twocircle.post_contact_frame;
outTrial.twocircle_post_deactivate_frame = twocircle.post_deactivate_frame;

outTrial.visibility_alpha = alphaProfile(:);
outTrial.visibility_geom = geomProfile(:);
end
