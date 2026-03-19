%% Stimuli Generation v19 (occlusion paradigm with matched pre-deviance branch)
% Script: stimuli_generation_v19.m
% Author: Codex (v19)
%
% Purpose:
%   Build a new three-condition one-dot occlusion dataset from an existing
%   one-dot source dataset, while preserving trajectory identity across
%   conditions and adding explicit per-trial occlusion metadata.
%
%   This script is designed for the new occlusion paradigm where:
%     1) always_visible         : trajectory always visible
%     2) occluded_nondeviant    : same trajectory with invisible-object occlusion
%     3) occluded_deviant       : deviant trajectory with invisible-object occlusion
%
%   The default occlusion model is "twocircle" and "alpha" is retained as
%   an optional fallback model.
%
% Occlusion timing controls:
%   - leadSec  (X): occlusion starts at deviance - X seconds
%   - lagSec   (Y): reappearance starts at deviance + Y seconds
%   - defaults: X = 0.25, Y = 0.25
%
% Twocircle model metadata (saved per trial):
%   - two circles are centered on the deviance point
%   - pre circle is active only before deviance
%   - post circle is active from deviance onward
%   - pre/post radii are computed so disk-edge contact occurs at
%     deviance-X and deviance+Y respectively
%   - post circle includes a first-full-exit frame to prevent recurrent
%     disappearances if the path later re-enters the circle
%
% Trigger-related metadata (saved per trial):
%   - occlusion_start_frame
%   - occlusion_complete_frame (first fully invisible frame in twocircle model)
%   - occlusion_end_frame      (first frame of reappearance)
%   - occlusion_end_complete_frame (first fully visible frame after reappearance)
%
% Usage example (interactive from experiment/):
%   addpath('lib');
%   stimuli_generation_v19;
%
% Usage example (non-interactive from repo root):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); addpath('lib'); stimuli_generation_v19;"
%
% Inputs:
%   - Config class (for input/output directories)
%   - Source observed file:  input_files/MovDot_SubXX.mat
%   - Source predicted file: input_files/MovDot_SubXX_predicted.mat
%
% Outputs:
%   - Observed target file:
%       input_files/MovDot_SubYY.mat
%       with 3 condition labels and occlusion metadata
%   - Predicted target file:
%       input_files/MovDot_SubYY_predicted.mat
%       with occluded_deviant_predicted trials and matching metadata
%
% Key assumptions:
%   - Source files are one-dot trajectories (xy columns [x y]).
%   - Source observed contains both nondeviant and deviant trials.
%   - Source predicted contains deviant predicted baseline trajectories.
%   - Trial pairing uses sequence sorting and one-to-one alignment.
%
% v19 deviant branch rule:
%   - occluded_deviant is explicitly built so its pre-deviance segment is
%     identical to the paired nondeviant trajectory.
%   - Post-deviance comes from the deviant source branch, translated so the
%     splice is position-continuous at deviance.

addpath('lib/');

%% Setup and user parameters
% Data flow: user input -> source/target IDs + occlusion timing -> build controls.
clearvars -except sourceSubjectID targetSubjectID leadSec lagSec overwriteExisting;
clc;

if ~(exist('sourceSubjectID', 'var') && exist('targetSubjectID', 'var') && ...
        exist('leadSec', 'var') && exist('lagSec', 'var'))
    prompt = { ...
        'Source subject ID (existing one-dot dataset):', ...
        'Target subject ID (new occlusion dataset):', ...
        'Occlusion lead X in seconds (deviance - X):', ...
        'Occlusion lag Y in seconds (deviance + Y):'};
    dlgTitle = 'v19 occlusion dataset parameters';
    dims = [1 70];
    definput = {'73', '70', '0.25', '0.25'};
    userInput = inputdlg(prompt, dlgTitle, dims, definput);
    if isempty(userInput)
        disp('Canceled by user.');
        return;
    end

    sourceSubjectID = str2double(userInput{1});
    targetSubjectID = str2double(userInput{2});
    leadSec = str2double(userInput{3});
    lagSec = str2double(userInput{4});
end

if any(isnan([sourceSubjectID, targetSubjectID, leadSec, lagSec]))
    error('All input fields must be numeric.');
end
if sourceSubjectID < 0 || targetSubjectID < 0
    error('Subject IDs must be >= 0.');
end
if leadSec <= 0 || lagSec <= 0
    error('leadSec and lagSec must be > 0.');
end
if sourceSubjectID == targetSubjectID
    warning('Source and target IDs are equal (%d). Target files may overwrite source IDs.', targetSubjectID);
end

sourceSubjectID = round(sourceSubjectID);
targetSubjectID = round(targetSubjectID);

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
leadFrames = max(1, round(leadSec * fps));
lagFrames = max(1, round(lagSec * fps));
alphaFadeInFrames = leadFrames;

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

    devFrame = local_infer_deviance_frame(trialDev, framesPerTrial);
    timing = local_build_occlusion_timing(framesPerTrial, devFrame, leadFrames, lagFrames, alphaFadeInFrames);
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

    devFrame = local_infer_deviance_frame(trialObservedDev, framesPerTrial);
    timing = local_build_occlusion_timing(framesPerTrial, devFrame, leadFrames, lagFrames, alphaFadeInFrames);
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
repro.v19_occlusion = struct( ...
    'script', mfilename, ...
    'sourceSubjectID', sourceSubjectID, ...
    'targetSubjectID', targetSubjectID, ...
    'sourceObservedFile', sourceObservedFile, ...
    'sourcePredictedFile', sourcePredictedFile, ...
    'leadSec', leadSec, ...
    'lagSec', lagSec, ...
    'leadFrames', leadFrames, ...
    'lagFrames', lagFrames, ...
    'alphaFadeInFrames', alphaFadeInFrames, ...
    'dotRadiusDeg', dotRadiusDeg, ...
    'deviantConstructionRule', 'nondeviant_prefix_plus_translated_deviant_suffix', ...
    'defaultModel', 'twocircle', ...
    'optionalModel', 'alpha', ...
    'conditionCodes', struct('always_visible', -1, 'occluded_nondeviant', 0, 'occluded_deviant', 45));

save(targetObservedFile, 'xySeqs', 'Cfg', 'repro');
save(targetPredictedFile, 'xySeqsPredicted', 'Cfg');

fprintf('Saved observed occlusion dataset: %s\n', targetObservedFile);
fprintf('Saved predicted occlusion dataset: %s\n', targetPredictedFile);
fprintf('Trials per condition: %d (always_visible / occluded_nondeviant / occluded_deviant).\n', nTrials);
fprintf('Occlusion timing defaults used: deviance - %.3fs to deviance + %.3fs\n', leadSec, lagSec);

%% Local helper functions
function devFrame = local_infer_deviance_frame(trialStruct, nFrames)
% LOCAL_INFER_DEVIANCE_FRAME Infer deviance onset from pathAll markers.
%
% Rule:
%   Use first non-initial pathAll change if available; otherwise 50%% frame.
devFrame = max(1, min(nFrames, round(0.5 * nFrames)));
if ~isfield(trialStruct, 'pathAll') || isempty(trialStruct.pathAll)
    return;
end
pathAll = trialStruct.pathAll;
if size(pathAll, 1) ~= nFrames
    return;
end
changeIdx = find(pathAll(2:end, 1) > 0, 1, 'first');
if ~isempty(changeIdx)
    devFrame = changeIdx + 1;
end
end

function trialSpliced = local_build_spliced_occluded_deviant(trialNondev, trialDev, devFrame)
% LOCAL_BUILD_SPLICED_OCCLUDED_DEVIANT Build occluded-deviant with nondeviant prefix.
%
% Rule:
%   - Frames 1..deviance_frame are copied from the nondeviant trial.
%   - Frames after deviance come from the deviant branch, translated so the
%     splice is position-continuous at deviance.
trialSpliced = trialDev;

xyNondev = double(trialNondev.xy);
xyDeviant = double(trialDev.xy);
if size(xyNondev, 1) ~= size(xyDeviant, 1) || size(xyNondev, 2) ~= 2 || size(xyDeviant, 2) ~= 2
    error('Splice builder expects one-dot trajectories with equal frame count.');
end

nFrames = size(xyDeviant, 1);
spliceFrame = max(1, min(nFrames, round(devFrame)));

xyOut = xyDeviant;
xyOut(1:spliceFrame, :) = xyNondev(1:spliceFrame, :);

if spliceFrame < nFrames
    shiftVec = xyOut(spliceFrame, :) - xyDeviant(spliceFrame, :);
    xyOut((spliceFrame + 1):end, :) = xyDeviant((spliceFrame + 1):end, :) + shiftVec;
end

trialSpliced.xy = xyOut;
end

function timing = local_build_occlusion_timing(nFrames, devFrame, leadFrames, lagFrames, fadeInFrames)
% LOCAL_BUILD_OCCLUSION_TIMING Build frame indices for all occlusion events.
%
% Data flow:
%   deviance + user lead/lag -> start/end/complete event frames.
occlusionStartFrame = max(1, devFrame - leadFrames);
reappearanceStartFrame = min(nFrames, devFrame + lagFrames);
reappearanceEndFrame = min(nFrames, reappearanceStartFrame + fadeInFrames);

timing = struct();
timing.deviance_frame = devFrame;
timing.occlusion_start_frame = occlusionStartFrame;
timing.occlusion_complete_frame = devFrame; % default for alpha; twocircle may override.
timing.occlusion_end_frame = reappearanceStartFrame;
timing.reappearance_start_frame = reappearanceStartFrame;
timing.reappearance_end_frame = reappearanceEndFrame;
timing.occlusion_end_complete_frame = reappearanceEndFrame; % default for alpha; twocircle may override.
end

function twocircle = local_build_twocircle_metadata(xy, devFrame, timing, dotRadiusDeg)
% LOCAL_BUILD_TWOCIRCLE_METADATA Compute twocircle center/radii and event frames.
%
% Geometry:
%   pre contact:  d = R_pre + r_dot  at occlusion_start_frame
%   post contact: d = R_post - r_dot at reappearance_start_frame
nFrames = size(xy, 1);
center = double(xy(devFrame, 1:2));
preFrame = max(1, min(nFrames, timing.occlusion_start_frame));
postFrame = max(1, min(nFrames, timing.reappearance_start_frame));

dPre = norm(double(xy(preFrame, 1:2)) - center);
dPost = norm(double(xy(postFrame, 1:2)) - center);

rPre = max(1e-6, dPre - dotRadiusDeg);
rPost = max(1e-6, dPost + dotRadiusDeg);

% Pre circle active strictly before deviance.
preActive = false(nFrames, 1);
preActive(1:max(1, devFrame - 1)) = true;

% Post circle active from deviance onward until first full exit.
postActive = false(nFrames, 1);
postActive(devFrame:nFrames) = true;

% Find first full exit for post circle to prevent recurring re-occlusion.
postDeactivateFrame = nFrames;
for iFrame = devFrame:nFrames
    distCenter = norm(double(xy(iFrame, 1:2)) - center);
    if distCenter >= (rPost + dotRadiusDeg)
        postDeactivateFrame = iFrame;
        break;
    end
end
if postDeactivateFrame < nFrames
    postActive((postDeactivateFrame + 1):end) = false;
end

% Compute first fully invisible / first fully visible-again event frames.
% Fully invisible: entire dot inside active occluder circle.
occlusionCompleteFrame = NaN;
for iFrame = timing.occlusion_start_frame:nFrames
    [isInvisible, ~] = local_twocircle_visibility_state( ...
        double(xy(iFrame, 1:2)), center, rPre, rPost, dotRadiusDeg, ...
        preActive(iFrame), postActive(iFrame));
    if isInvisible
        occlusionCompleteFrame = iFrame;
        break;
    end
end
if isnan(occlusionCompleteFrame)
    occlusionCompleteFrame = timing.deviance_frame;
end

% First fully visible-again frame after reappearance starts.
occlusionEndCompleteFrame = timing.reappearance_end_frame;
for iFrame = timing.reappearance_start_frame:nFrames
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
twocircle.occlusion_complete_frame = occlusionCompleteFrame;
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
outTrial.occlusion_start_frame = timing.occlusion_start_frame;
outTrial.occlusion_complete_frame = twocircle.occlusion_complete_frame;
outTrial.occlusion_end_frame = timing.occlusion_end_frame;
outTrial.reappearance_start_frame = timing.reappearance_start_frame;
outTrial.reappearance_end_frame = timing.reappearance_end_frame;
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
