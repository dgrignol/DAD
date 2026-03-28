%% Stimuli Generation V20 (fixed-frame occlusion window with matched pre-deviance branch)
% Script: stimuli_generation_V20.m
% Author: Dami (V20)
%
% Purpose:
%   Build a three-condition one-dot occlusion dataset from an existing
%   one-dot source dataset while enforcing a fixed occlusion timeline:
%     - deviance frame fixed at 130
%     - first fully occluded frame fixed at 130
%     - fully occluded through frame 190 (inclusive)
%     - nominal reappearance search starts at frame 191
%       (actual reappearance frame is computed geometrically)
%
%   Conditions generated:
%     1) always_visible
%     2) occluded_nondeviant
%     3) occluded_deviant
%
% Occlusion model:
%   - default model remains "twocircle" (alpha fallback retained).
%   - pre-deviance occluder:
%       * center at xy(frame 130)
%       * radius equal to moving-dot radius
%       * active only for frames < 130
%   - post-deviance occluder:
%       * activates at frame 130
%       * center also anchored at xy(frame 130)
%       * radius computed to fully occlude every frame in 130..190
%       * remains active until trial end so reappearance is geometric
%         (gradual when trajectory exits the occluder boundary)
%
% Trigger-related metadata (saved per trial):
%   - occlusion_start_frame
%   - occlusion_complete_frame (first fully invisible frame)
%   - occlusion_end_frame      (first frame of reappearance)
%   - occlusion_end_complete_frame (first fully visible frame after reappearance)
%
% Usage example (interactive from experiment/):
%   addpath('lib');
%   stimuli_generation_V20;
%
% Default parameter values when variables are not pre-set:
%   sourceSubjectID = 73
%   targetSubjectID = prompted every run (no default)
%   fixedDevianceFrame = 130
%   fixedOcclusionEndFrame = 190
%   overwriteExisting = true
%
% Usage example (non-interactive from repo root):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    addpath('lib'); sourceSubjectID=73; targetSubjectID=99; ...
%    fixedDevianceFrame=130; fixedOcclusionEndFrame=190; ...
%    overwriteExisting=true; stimuli_generation_V20;"
%
% Inputs:
%   - Config class (for input/output directories)
%   - Source observed file:  input_files/MovDot_SubXX.mat
%   - Source predicted file: input_files/MovDot_SubXX_predicted.mat
%
% Outputs:
%   - Observed target file:  input_files/MovDot_SubYY.mat
%   - Predicted target file: input_files/MovDot_SubYY_predicted.mat
%   Both outputs include condition labels and occlusion metadata.
%
% Key assumptions:
%   - Source files are one-dot trajectories (xy columns [x y]).
%   - Source observed contains both nondeviant and deviant trials.
%   - Source predicted contains deviant predicted baseline trajectories.
%   - Trial pairing uses sequence sorting and one-to-one alignment.
%
% V20 deviant branch rule:
%   - occluded_deviant is explicitly built so its pre-deviance segment is
%     identical to the paired nondeviant trajectory.
%   - post-deviance comes from the deviant source branch, translated so the
%     splice is position-continuous at frame 130.

addpath('lib/');

%% Setup and user parameters
% Data flow: caller-provided vars (or defaults) -> build controls.
clearvars -except sourceSubjectID targetSubjectID overwriteExisting fixedDevianceFrame fixedOcclusionEndFrame;
clc;

% Default values when variables are not provided by the caller.
if ~(exist('sourceSubjectID', 'var') && ~isempty(sourceSubjectID))
    sourceSubjectID = 73;
end
if ~(exist('fixedDevianceFrame', 'var') && ~isempty(fixedDevianceFrame))
    fixedDevianceFrame = 130;
end
if ~(exist('fixedOcclusionEndFrame', 'var') && ~isempty(fixedOcclusionEndFrame))
    fixedOcclusionEndFrame = 190;
end
if ~(exist('overwriteExisting', 'var') && ~isempty(overwriteExisting))
    overwriteExisting = true;
end

% Safety gate: always request an explicit target subject ID.
prompt = {'Target subject ID (new occlusion dataset):'};
dlgTitle = 'V20 target subject';
dims = [1 70];
definput = {''};
userInput = inputdlg(prompt, dlgTitle, dims, definput);
if isempty(userInput)
    disp('Canceled by user.');
    return;
end
targetSubjectID = str2double(userInput{1});

if any(isnan([sourceSubjectID, targetSubjectID]))
    error('Source and target subject IDs must be numeric.');
end
if sourceSubjectID < 0 || targetSubjectID < 0
    error('Subject IDs must be >= 0.');
end
if sourceSubjectID == targetSubjectID
    warning('Source and target IDs are equal (%d). Target files may overwrite source IDs.', targetSubjectID);
end

sourceSubjectID = round(sourceSubjectID);
targetSubjectID = round(targetSubjectID);
if ~isscalar(fixedDevianceFrame) || ~isscalar(fixedOcclusionEndFrame) || ...
        isnan(fixedDevianceFrame) || isnan(fixedOcclusionEndFrame)
    error('fixedDevianceFrame and fixedOcclusionEndFrame must be numeric scalars.');
end
fixedDevianceFrame = round(fixedDevianceFrame);
fixedOcclusionEndFrame = round(fixedOcclusionEndFrame);
overwriteExisting = logical(overwriteExisting);

%% Resolve file paths and guard outputs
% Data flow: subject IDs + Config directories -> absolute source/target filenames.
inputDir = Config.inputDirectory;
if ~exist(inputDir, 'dir')
    mkdir(inputDir);
end

sourceObservedFile = fullfile(inputDir, sprintf('MovDot_Sub%02d.mat', sourceSubjectID));
sourcePredictedFile = fullfile(inputDir, sprintf('MovDot_Sub%02d_predicted.mat', sourceSubjectID));
targetObservedFile = fullfile(inputDir, sprintf('MovDot_Sub%02d.mat', targetSubjectID));
targetPredictedFile = fullfile(inputDir, sprintf('MovDot_Sub%02d_predicted.mat', targetSubjectID));

if ~isfile(sourceObservedFile)
    error('Source observed file not found: %s', sourceObservedFile);
end
if ~isfile(sourcePredictedFile)
    error('Source predicted file not found: %s', sourcePredictedFile);
end

if isfile(targetObservedFile) || isfile(targetPredictedFile)
    doOverwrite = false;
    if exist('overwriteExisting', 'var') && ~isempty(overwriteExisting)
        doOverwrite = logical(overwriteExisting);
    end
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

%% Load source datasets and validate one-dot format
% Data flow: source MATs -> flat trial arrays -> one-dot condition pools.
obsData = load(sourceObservedFile, 'xySeqs', 'Cfg', 'repro');
predData = load(sourcePredictedFile, 'xySeqsPredicted', 'Cfg');

if ~isfield(obsData, 'xySeqs') || ~isfield(obsData, 'Cfg')
    error('Observed source must contain xySeqs and Cfg: %s', sourceObservedFile);
end
if ~isfield(predData, 'xySeqsPredicted')
    error('Predicted source must contain xySeqsPredicted: %s', sourcePredictedFile);
end

allObserved = obsData.xySeqs(:);
allPredicted = predData.xySeqsPredicted(:);
if isempty(allObserved)
    error('Observed source contains no trials.');
end
if isempty(allPredicted)
    error('Predicted source contains no trials.');
end

obsValid = arrayfun(@(s) isnumeric(s.xy) && ~isempty(s.xy) && size(s.xy, 2) >= 2, allObserved);
predValid = arrayfun(@(s) isnumeric(s.xy) && ~isempty(s.xy) && size(s.xy, 2) >= 2, allPredicted);
allObserved = allObserved(obsValid);
allPredicted = allPredicted(predValid);

obsCols = unique(arrayfun(@(s) size(s.xy, 2), allObserved));
predCols = unique(arrayfun(@(s) size(s.xy, 2), allPredicted));
if any(obsCols ~= 2)
    error('Observed source must be one-dot xy=[x y]. Detected columns: %s', mat2str(obsCols));
end
if any(predCols ~= 2)
    error('Predicted source must be one-dot xy=[x y]. Detected columns: %s', mat2str(predCols));
end

observedNondeviant = allObserved([allObserved.condition] == 0);
observedDeviant = allObserved([allObserved.condition] > 0);
predictedDeviant = allPredicted([allPredicted.condition] > 0);

if isempty(observedNondeviant) || isempty(observedDeviant) || isempty(predictedDeviant)
    error('Failed to split source trials into nondeviant/deviant pools.');
end
if numel(observedNondeviant) ~= numel(observedDeviant)
    error('Observed source must have equal nondeviant and deviant trial counts.');
end
if numel(predictedDeviant) ~= numel(observedDeviant)
    error('Predicted deviant trial count (%d) must equal observed deviant count (%d).', ...
        numel(predictedDeviant), numel(observedDeviant));
end

% Deterministic alignment by sequence index when available.
if isfield(observedNondeviant, 'sequence')
    [~, idx] = sort([observedNondeviant.sequence]);
    observedNondeviant = observedNondeviant(idx);
end
if isfield(observedDeviant, 'sequence')
    [~, idx] = sort([observedDeviant.sequence]);
    observedDeviant = observedDeviant(idx);
end
if isfield(predictedDeviant, 'sequence')
    [~, idx] = sort([predictedDeviant.sequence]);
    predictedDeviant = predictedDeviant(idx);
end

nTrials = numel(observedNondeviant);
framesPerTrial = size(observedNondeviant(1).xy, 1);
fps = double(obsData.Cfg.fps);

if fixedDevianceFrame < 1 || fixedDevianceFrame > framesPerTrial
    error('fixedDevianceFrame (%d) must be within [1, %d].', fixedDevianceFrame, framesPerTrial);
end
if fixedOcclusionEndFrame < fixedDevianceFrame || fixedOcclusionEndFrame > framesPerTrial
    error(['fixedOcclusionEndFrame (%d) must be within [%d, %d] and must not ' ...
        'precede fixedDevianceFrame.'], fixedOcclusionEndFrame, fixedDevianceFrame, framesPerTrial);
end

% Alpha fallback uses a short fade-in window after reappearance.
alphaFadeInFrames = max(1, round(0.25 * fps));

if ~isfield(obsData.Cfg, 'dot_w')
    error('Source Cfg.dot_w is required to derive twocircle radii.');
end
dotRadiusDeg = double(obsData.Cfg.dot_w) / 2;

%% Build observed output trials for three conditions
% Data flow: paired source trials -> condition-specific metadata -> flattened xySeqs.
observedCell = cell(3 * nTrials, 1);
splicedObservedDeviant = cell(nTrials, 1);

for iTrial = 1:nTrials
    trialNondev = observedNondeviant(iTrial);
    trialDev = observedDeviant(iTrial);

    devFrame = local_infer_deviance_frame(framesPerTrial, fixedDevianceFrame);
    timing = local_build_occlusion_timing( ...
        framesPerTrial, devFrame, fixedOcclusionEndFrame, alphaFadeInFrames);
    trialDevSpliced = local_build_spliced_occluded_deviant(trialNondev, trialDev, devFrame);

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
        trialDevSpliced, 45, 'occluded_deviant', true, ...
        timing, twocircleDev, alphaOcc, geomOcc);
    splicedObservedDeviant{iTrial} = observedCell{rowOccDev};
end
xySeqs = vertcat(observedCell{:});

%% Build predicted output trials (occluded deviant only)
% Data flow: predicted source + spliced observed deviant prefix override -> xySeqsPredicted.
predictedCell = cell(nTrials, 1);
for iTrial = 1:nTrials
    trialObservedDev = splicedObservedDeviant{iTrial};
    trialPredictedDev = predictedDeviant(iTrial);

    devFrame = local_infer_deviance_frame(framesPerTrial, fixedDevianceFrame);
    timing = local_build_occlusion_timing( ...
        framesPerTrial, devFrame, fixedOcclusionEndFrame, alphaFadeInFrames);
    twocirclePred = local_build_twocircle_metadata(trialPredictedDev.xy, devFrame, timing, dotRadiusDeg);
    [alphaOcc, geomOcc] = local_build_alpha_profiles(framesPerTrial, timing);

    % Modeling assumption: predicted matches observed until reappearance starts.
    adjustedPred = trialPredictedDev;
    rs = timing.reappearance_start_frame;
    if rs > 1
        adjustedPred.xy(1:(rs - 1), :) = trialObservedDev.xy(1:(rs - 1), :);
    end

    predictedCell{iTrial} = local_attach_occlusion_fields( ...
        adjustedPred, 45, 'occluded_deviant_predicted', true, ...
        timing, twocirclePred, alphaOcc, geomOcc);
end
xySeqsPredicted = vertcat(predictedCell{:});

%% Save target datasets with reproducibility metadata
% Data flow: generated trials + source Cfg -> output files and repro snapshot.
Cfg = obsData.Cfg;

repro = struct();
if isfield(obsData, 'repro')
    repro = obsData.repro;
end
repro.v20_occlusion = struct( ...
    'script', mfilename, ...
    'sourceSubjectID', sourceSubjectID, ...
    'targetSubjectID', targetSubjectID, ...
    'sourceObservedFile', sourceObservedFile, ...
    'sourcePredictedFile', sourcePredictedFile, ...
    'fixedDevianceFrame', fixedDevianceFrame, ...
    'fixedOcclusionEndFrame', fixedOcclusionEndFrame, ...
    'nominalReappearanceFrame', fixedOcclusionEndFrame + 1, ...
    'alphaFadeInFrames', alphaFadeInFrames, ...
    'dotRadiusDeg', dotRadiusDeg, ...
    'deviantConstructionRule', 'nondeviant_prefix_plus_translated_deviant_suffix', ...
    'twocircleRule', 'pre_radius_equals_dot_radius;post_radius_covers_frames_130_to_190;post_active_until_trial_end', ...
    'defaultModel', 'twocircle', ...
    'optionalModel', 'alpha', ...
    'conditionCodes', struct('always_visible', -1, 'occluded_nondeviant', 0, 'occluded_deviant', 45));

save(targetObservedFile, 'xySeqs', 'Cfg', 'repro');
save(targetPredictedFile, 'xySeqsPredicted', 'Cfg');

fprintf('Saved observed occlusion dataset: %s\n', targetObservedFile);
fprintf('Saved predicted occlusion dataset: %s\n', targetPredictedFile);
fprintf('Trials per condition: %d (always_visible / occluded_nondeviant / occluded_deviant).\n', nTrials);
fprintf(['Fixed occlusion timing: deviance=%d, first full occlusion=%d, ' ...
    'last required full occlusion=%d. Reappearance is computed geometrically.\n'], ...
    fixedDevianceFrame, fixedDevianceFrame, fixedOcclusionEndFrame);

%% Local helper functions
function devFrame = local_infer_deviance_frame(nFrames, fixedDevianceFrame)
% LOCAL_INFER_DEVIANCE_FRAME Clamp fixed V20 deviance frame to trial length.
devFrame = max(1, min(nFrames, round(fixedDevianceFrame)));
end

function trialSpliced = local_build_spliced_occluded_deviant(trialNondev, trialDev, devFrame)
% LOCAL_BUILD_SPLICED_OCCLUDED_DEVIANT Build occluded-deviant with nondeviant prefix.
%
% Rule in V20:
%   - Frames 1..fixed_deviance_frame are copied from the nondeviant trial.
%   - The post-deviance branch is taken from the source deviant trajectory
%     starting at its own (source) deviance frame, then time-resampled to the
%     remaining V20 trial duration and translated so continuity holds exactly
%     at frame 130.
%   - This removes the second inherited deviance point and keeps one deviance
%     transition at the requested fixed frame.
trialSpliced = trialDev;

xyNondev = double(trialNondev.xy);
xyDeviant = double(trialDev.xy);
if size(xyNondev, 1) ~= size(xyDeviant, 1) || size(xyNondev, 2) ~= 2 || size(xyDeviant, 2) ~= 2
    error('Splice builder expects one-dot trajectories with equal frame count.');
end

nFrames = size(xyDeviant, 1);
spliceFrame = max(1, min(nFrames, round(devFrame)));
sourceDevFrame = local_infer_source_deviance_frame(trialDev, nFrames);
sourceDevFrame = max(1, min(nFrames, sourceDevFrame));

xyOut = xyNondev;
xyOut(1:spliceFrame, :) = xyNondev(1:spliceFrame, :);

if spliceFrame < nFrames
    % Resample the source post-deviance branch onto the V20 post-deviance
    % frame span to preserve total trial length while moving deviance onset.
    sourcePost = xyDeviant(sourceDevFrame:end, :);
    sourceCount = size(sourcePost, 1);
    targetCount = nFrames - spliceFrame + 1;

    if sourceCount < 2
        postResampled = repmat(sourcePost(1, :), targetCount, 1);
    else
        tSource = linspace(0, 1, sourceCount);
        tTarget = linspace(0, 1, targetCount);
        postResampled = zeros(targetCount, 2);
        postResampled(:, 1) = interp1(tSource, sourcePost(:, 1), tTarget, 'linear');
        postResampled(:, 2) = interp1(tSource, sourcePost(:, 2), tTarget, 'linear');
    end

    shiftVec = xyOut(spliceFrame, :) - postResampled(1, :);
    postAligned = postResampled + shiftVec;
    xyOut(spliceFrame:end, :) = postAligned;
end

trialSpliced.xy = xyOut;
end

function sourceDevFrame = local_infer_source_deviance_frame(trialStruct, nFrames)
% LOCAL_INFER_SOURCE_DEVIANCE_FRAME Infer source deviance onset from pathAll.
%
% Fallback rule:
%   If pathAll is missing/invalid, use midpoint.
sourceDevFrame = max(1, min(nFrames, round(0.5 * nFrames)));
if ~isfield(trialStruct, 'pathAll') || isempty(trialStruct.pathAll)
    return;
end

pathAll = trialStruct.pathAll(:);
if numel(pathAll) ~= nFrames
    return;
end

changeIdx = find(pathAll(2:end) > 0, 1, 'first');
if ~isempty(changeIdx)
    sourceDevFrame = changeIdx + 1;
end
end

function timing = local_build_occlusion_timing(nFrames, devFrame, occlusionEndFrameInclusive, fadeInFrames)
% LOCAL_BUILD_OCCLUSION_TIMING Build fixed-frame V20 occlusion event indices.
%
% Data flow:
%   fixed deviance + fixed occlusion end -> start/end/complete event frames.
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
% LOCAL_BUILD_TWOCIRCLE_METADATA Compute V20 twocircle geometry and event frames.
%
% Geometry in V20:
%   - pre circle: radius exactly equals dot radius, centered at frame 130 position.
%   - post circle: radius large enough to keep dot fully hidden from frame 130
%     through frame 190 inclusive, then remains active to trial end.
nFrames = size(xy, 1);
center = double(xy(devFrame, 1:2));
preFrame = devFrame;
postFrame = devFrame;
holdEndFrame = max(devFrame, min(nFrames, timing.occlusion_full_end_frame));

% Pre-occluder is explicitly dot-sized at deviance location.
rPre = max(1e-6, dotRadiusDeg);

% Post-occluder radius covers every dot center in the full-occlusion window.
windowXY = double(xy(devFrame:holdEndFrame, 1:2));
distWindow = sqrt(sum((windowXY - center) .^ 2, 2));
rPost = max(1e-6, max(distWindow) + dotRadiusDeg);

% Pre circle active strictly before deviance.
preActive = false(nFrames, 1);
preActive(1:max(1, devFrame - 1)) = true;

% Post circle stays active until trial end to support gradual reappearance.
postActive = false(nFrames, 1);
postActive(devFrame:nFrames) = true;
postDeactivateFrame = nFrames;

% Compute first partial occlusion and first full occlusion using geometric state.
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

% First reappearance frame (first frame that is not fully invisible after the
% enforced full-occlusion window) and first fully visible-again frame.
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
%
% Simplified analytic state test for a disk under active circle occluders.
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

% Invisible if dot disk is entirely inside any active circle occluder.
isInvisible = any(distCenter <= (activeRadii - dotRadiusDeg));

% Fully visible if dot disk is entirely outside all active circles.
isFullyVisible = all(distCenter >= (activeRadii + dotRadiusDeg));
end

function [alphaProfile, geomProfile] = local_build_visible_profiles(nFrames)
% LOCAL_BUILD_VISIBLE_PROFILES Visibility profiles for always-visible condition.
alphaProfile = ones(nFrames, 1);
geomProfile = ones(nFrames, 1);
end

function [alphaProfile, geomProfile] = local_build_alpha_profiles(nFrames, timing)
% LOCAL_BUILD_ALPHA_PROFILES Smooth alpha fallback profile for occluded conditions.
%
% Timeline:
%   - occlusion_start_frame -> deviance_frame : fade out 1 -> 0
%   - deviance_frame+1 -> reappearance_start_frame-1 : hold 0
%   - reappearance_start_frame -> reappearance_end_frame : fade in 0 -> 1

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
