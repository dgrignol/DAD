% MoveDot1_experiment_occlusion_v18_rescueTraject.m
%
% Purpose:
%   Run the one-dot occlusion paradigm with three runs per block and
%   catch-trial logic planned upstream by:
%     CreateInputFiles_v21_rescueTraject.m
%
% Schedule model (from v21 TrialOrder/CatchPlan):
%   - Run 1: always_visible base trials + inserted type-1 catches.
%   - Run 2: mixed occluded_nondeviant + occluded_deviant base trials +
%            inserted type-2 catches.
%   - Run 3: mixed occluded_nondeviant + occluded_deviant base trials +
%            inserted type-2 catches.
%
% Timing-focused v18 rescueTraject additions:
%   - SkipSyncTests policy is no longer forced in MEG mode.
%   - Full flip diagnostics are stored per frame:
%       * VBLTimestamp
%       * StimulusOnsetTime
%       * FlipTimestamp
%       * Missed
%   - Online guardrails track missed flips and can warn or abort.
%   - Occluder draw geometry is precomputed per trial (including straight
%     terminal transforms and polygon assembly) to reduce per-frame load.
%   - EyeLink support is restored with optional fixation-break detection,
%     replay scheduling, and EyeLink trigger/event messages.
%   - EyeLink setup/calibration can be requested even when fake gaze replay
%     mode is enabled (trackeye=1, fakeEyeTracker=1, eyeTrackerDoCalibration=1).
%
% Trigger logic (MEG/DataPixx):
%   - Condition and occlusion triggers are preserved.
%   - Sequence-identity trigger fires once per trial at fixed frame offset.
%   - Catch question emits one mutually exclusive end-status trigger:
%       * 201 catch_question_correct
%       * 202 catch_question_incorrect
%       * 113 catch_question_timeout
%   - Eye-tracker-specific trigger codes are restored:
%       * 150 gaze break
%       * 151 replay start
%
% Usage example (interactive from experiment/):
%   MoveDot1_experiment_occlusion_v18_rescueTraject
%
% Usage example (non-interactive shell from repo root):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    iSub=66; iBlock=1; viewingDistanceMm=1000; ...
%    run('MoveDot1_experiment_occlusion_v18_rescueTraject.m');"
%
% Inputs:
%   - experiment/input_files/MovDot_SubXX_V28_rescueTraject.mat
%   - experiment/input_files/SubXX_TrialStruct_v21_rescueTraject.mat
%
% Workspace overrides (optional):
%   - iSub, iBlock, viewingDistanceMm
%   - dryRunValidateScheduleOnly, debug_runNumber, dryRunSimulatedBreakKeys
%   - question_text, debugTriggerMode, runColorCueEnabled
%   - practice_mode, practiceRunPercentageTrials
%   - ITIRangeSec (seconds, [min max], default [0.500 1.000] when unset)
%   - debugShowFrameInfoOverlay, debugShowFullPathOverlay
%   - pathBandEntranceStyle and straight-terminal tuning parameters
%   - debugSkipSyncTests (0/1; explicit development override for SkipSyncTests)
%   - disablePriorityBoost (0/1; practice safety guard for unsupported
%     PTB Priority/Mach MEX environments)
%   - flipMissedWarnThresholdCount / flipMissedAbortEnabled /
%     flipMissedAbortThresholdCount (runtime guardrail overrides)
%   - trackeye, ignoreEyeTracker, fakeEyeTracker, fakeGazeBreakRate
%   - enableFixationAbort, fixWindowDeg, fixBreakToleranceFrames
%   - fixationWarningDurationSec, replayLimit, allowInfiniteReplays
%   - eyeTrackerDoCalibration, eyeTrackerAllowBetweenBlockCalibration,
%     eyeTrackerSendTriggerMessages
%   - eyeTrackerImageTransferEnabled
%   - run1TransitionCountdownEnabled
%
% Outputs:
%   - experiment/output_files/MoveDot1_occlusion_v18_rescueTraject_SUBXX_BLOCKYY.mat
%   - experiment/output_files/debug_actual_triggers_occlusion_v18_rescueTraject_subXX_blockYY_runZZ.csv
%     (debug mode; columns: trigger, trial, frame, seconds, label, iti_duration_sec)
%
% Assumptions:
%   - Input trials include condition_label and pathband metadata fields.
%   - Input trajectories are one-dot [x y] in visual degrees.
%   - Generation fps in Cfg.fps matches runtime intent (nominally 120 Hz).
%   - EyeLink runtime is optional and enabled via config/workspace override.

%% Clear state and define runtime configuration
% Data flow: static conf values + user input -> runtime settings and file paths.
clearvars -except iSub iBlock viewingDistanceMm dryRunValidateScheduleOnly debug_runNumber dryRunSimulatedBreakKeys question_text debugTriggerMode runColorCueEnabled practice_mode practiceRunPercentageTrials ITIRangeSec debugShowFrameInfoOverlay debugShowFullPathOverlay pathBandEntranceStyle pathBandStraightBackshiftDotRadiusScale pathBandStraightEndForwardDotRadiusScale pathBandStraightPostStartWidthExtraDotRadiusScale pathBandStraightStartBackshiftPixelOffset debugSkipSyncTests disablePriorityBoost flipMissedWarnThresholdCount flipMissedAbortEnabled flipMissedAbortThresholdCount trackeye ignoreEyeTracker fakeEyeTracker fakeGazeBreakRate enableFixationAbort fixWindowDeg fixBreakToleranceFrames fixationWarningDurationSec replayLimit allowInfiniteReplays eyeTrackerDoCalibration eyeTrackerAllowBetweenBlockCalibration eyeTrackerSendTriggerMessages eyeTrackerImageTransferEnabled run1TransitionCountdownEnabled;
clc;
sca;
format shortG;
addpath('lib/');

if ~exist('dryRunValidateScheduleOnly', 'var') || isempty(dryRunValidateScheduleOnly)
    dryRunValidateScheduleOnly = false;
end
if ~exist('debug_runNumber', 'var')
    debug_runNumber = [];
end
if ~exist('debugTriggerMode', 'var') || isempty(debugTriggerMode)
    debugTriggerMode = 0;
end
if ~exist('debugShowFrameInfoOverlay', 'var') || isempty(debugShowFrameInfoOverlay)
    debugShowFrameInfoOverlay = false;
end
if ~exist('debugShowFullPathOverlay', 'var') || isempty(debugShowFullPathOverlay)
    debugShowFullPathOverlay = false;
end
if ~exist('debugSkipSyncTests', 'var')
    debugSkipSyncTests = [];
end
if ~exist('disablePriorityBoost', 'var') || isempty(disablePriorityBoost)
    disablePriorityBoost = true;
end
disablePriorityBoost = logical(disablePriorityBoost);
if ~exist('flipMissedWarnThresholdCount', 'var')
    flipMissedWarnThresholdCount = [];
end
if ~exist('flipMissedAbortEnabled', 'var')
    flipMissedAbortEnabled = [];
end
if ~exist('flipMissedAbortThresholdCount', 'var')
    flipMissedAbortThresholdCount = [];
end
if ~exist('trackeye', 'var')
    trackeye = [];
end
if ~exist('ignoreEyeTracker', 'var')
    ignoreEyeTracker = [];
end
if ~exist('fakeEyeTracker', 'var')
    fakeEyeTracker = [];
end
if ~exist('fakeGazeBreakRate', 'var')
    fakeGazeBreakRate = [];
end
if ~exist('enableFixationAbort', 'var')
    enableFixationAbort = [];
end
if ~exist('fixWindowDeg', 'var')
    fixWindowDeg = [];
end
if ~exist('fixBreakToleranceFrames', 'var')
    fixBreakToleranceFrames = [];
end
if ~exist('fixationWarningDurationSec', 'var')
    fixationWarningDurationSec = [];
end
if ~exist('replayLimit', 'var')
    replayLimit = [];
end
if ~exist('allowInfiniteReplays', 'var')
    allowInfiniteReplays = [];
end
if ~exist('eyeTrackerDoCalibration', 'var')
    eyeTrackerDoCalibration = [];
end
if ~exist('eyeTrackerAllowBetweenBlockCalibration', 'var')
    eyeTrackerAllowBetweenBlockCalibration = [];
end
if ~exist('eyeTrackerSendTriggerMessages', 'var')
    eyeTrackerSendTriggerMessages = [];
end
if ~exist('eyeTrackerImageTransferEnabled', 'var')
    eyeTrackerImageTransferEnabled = [];
end
if ~exist('run1TransitionCountdownEnabled', 'var')
    run1TransitionCountdownEnabled = [];
end
if ~exist('ITIRangeSec', 'var')
    ITIRangeSec = [];
end
debugShowFrameInfoOverlay = logical(debugShowFrameInfoOverlay);
debugShowFullPathOverlay = logical(debugShowFullPathOverlay);

scheduleCfg = Config_schedule_CreateInputV21_MoveDotV18_rescueTraject;
if scheduleCfg.runsPerBlock ~= 3
    error('Config_schedule_CreateInputV21_MoveDotV18_rescueTraject.runsPerBlock must be 3 for this runtime.');
end
if ~exist('runColorCueEnabled', 'var') || isempty(runColorCueEnabled)
    runColorCueEnabled = logical(scheduleCfg.runColorCueEnabled);
else
    runColorCueEnabled = logical(runColorCueEnabled);
end

Conf = struct();
Conf.MEG = 0;
Conf.debugTriggerMode = double(debugTriggerMode ~= 0);
Conf.debug = Conf.debugTriggerMode; % Backward-compatible alias used by legacy debug paths.
Conf.enableCatchTrials = 1;
Conf.occlusionMode = 'pathband'; % default: pathband; optional: alpha.
Conf.experimentName = 'MoveDot1_occlusion_v18_rescueTraject';
Conf.background = [0 0 0];
Conf.rectColor = [0 0 0];
Conf.rectBorderColor = [0 0 0];
Conf.dotColor = [0 255 0];
Conf.dotColorRun1Base = [0 255 0];
Conf.dotColorRuns23Base = [242 223 0];
Conf.runColorCueEnabled = logical(runColorCueEnabled);
Conf.fixColor = [180 180 180];
Conf.fixSizeDeg = 0.2;
Conf.showFixCrossBetweenTrials = 0; % kept for compatibility; not used in post-trial ITI.
% ITI policy (v18 jittered): one post-trial blank with fixation.
% Data flow: ITIRangeSec override/default -> precomputed per-run vector ->
% source-trial mapping (replay-safe) -> trial-specific WaitSecs duration.
Conf.ITIRangeSec = local_resolve_iti_range(ITIRangeSec, [0.500 1.000]);
Conf.postTrialBlank1Sec = Conf.ITIRangeSec(1); % fallback used only if trial ITI is invalid.
Conf.postTrialGridSec = 0.000;
Conf.postTrialBlank2Sec = 0.000;
% Legacy mask settings kept for easy restore:
% Conf.postTrialBlank1Sec = 0.150;
% Conf.postTrialGridSec = 0.150;
% Conf.postTrialBlank2Sec = 0.150;
% Conf.gridMaskBgColor = [45 45 45];
% Conf.gridMaskLineColor = [70 70 70];
% Conf.gridMaskSpacingPx = 24;
% Conf.gridMaskLineWidthPx = 1;
Conf.gridMaskBgColor = [45 45 45];
Conf.gridMaskLineColor = [70 70 70];
Conf.gridMaskSpacingPx = 24;
Conf.gridMaskLineWidthPx = 1;
Conf.sideMarginPx = 0;

% Eye-tracker controls (restored from legacy vX behavior, config-driven).
% Data flow:
%   schedule defaults -> workspace overrides -> derived runtime mode:
%   real EyeLink, fake EyeLink, or disabled.
Conf.trackeye = logical(local_schedule_get_or_default(scheduleCfg, 'eyeTrackerEnabled', false));
Conf.ignoreEyeTracker = logical(local_schedule_get_or_default(scheduleCfg, 'ignoreEyeTracker', false));
Conf.fakeEyeTracker = logical(local_schedule_get_or_default(scheduleCfg, 'fakeEyeTracker', false));
Conf.fakeGazeBreakRate = double(local_schedule_get_or_default(scheduleCfg, 'fakeGazeBreakRate', 0.05));
Conf.enableFixationAbort = logical(local_schedule_get_or_default(scheduleCfg, 'enableFixationAbort', true));
Conf.fixWindowDeg = double(local_schedule_get_or_default(scheduleCfg, 'fixWindowDeg', 7.0));
Conf.fixBreakToleranceFrames = round(double(local_schedule_get_or_default(scheduleCfg, 'fixBreakToleranceFrames', 5)));
Conf.fixationWarningDurationSec = double(local_schedule_get_or_default(scheduleCfg, 'fixationWarningDurationSec', 2.0));
Conf.replayLimit = double(local_schedule_get_or_default(scheduleCfg, 'replayLimit', 1));
Conf.allowInfiniteReplays = logical(local_schedule_get_or_default(scheduleCfg, 'allowInfiniteReplays', false));
Conf.eyeTrackerDoCalibration = logical(local_schedule_get_or_default(scheduleCfg, 'eyeTrackerDoCalibration', true));
Conf.eyeTrackerAllowBetweenBlockCalibration = logical(local_schedule_get_or_default( ...
    scheduleCfg, 'eyeTrackerAllowBetweenBlockCalibration', true));
Conf.eyeTrackerSendTriggerMessages = logical(local_schedule_get_or_default(scheduleCfg, 'eyeTrackerSendTriggerMessages', true));
Conf.eyeTrackerImageTransferEnabled = logical(local_schedule_get_or_default(scheduleCfg, 'eyeTrackerImageTransferEnabled', true));
Conf.eyeTrackerWaitSec = double(local_schedule_get_or_default(scheduleCfg, 'eyeTrackerWaitSec', 0.01));
Conf.eyeTrackerEdfPolicy = lower(strtrim(char(local_schedule_get_or_default( ...
    scheduleCfg, 'eyeTrackerEdfPolicy', 'per_block'))));
Conf.eyeTrackerReceiveFileAfterBlock = logical(local_schedule_get_or_default( ...
    scheduleCfg, 'eyeTrackerReceiveFileAfterBlock', true));
Conf.eyeTrackerFixBmpPath = char(local_schedule_get_or_default(scheduleCfg, 'eyeTrackerFixBmpPath', fullfile('.', 'fixation.bmp')));
if ~isempty(trackeye)
    Conf.trackeye = logical(trackeye);
end
if ~isempty(ignoreEyeTracker)
    Conf.ignoreEyeTracker = logical(ignoreEyeTracker);
end
if ~isempty(fakeEyeTracker)
    Conf.fakeEyeTracker = logical(fakeEyeTracker);
end
if ~isempty(fakeGazeBreakRate)
    Conf.fakeGazeBreakRate = double(fakeGazeBreakRate);
end
if ~isempty(enableFixationAbort)
    Conf.enableFixationAbort = logical(enableFixationAbort);
end
if ~isempty(fixWindowDeg)
    Conf.fixWindowDeg = double(fixWindowDeg);
end
if ~isempty(fixBreakToleranceFrames)
    Conf.fixBreakToleranceFrames = round(double(fixBreakToleranceFrames));
end
if ~isempty(fixationWarningDurationSec)
    Conf.fixationWarningDurationSec = double(fixationWarningDurationSec);
end
if ~isempty(replayLimit)
    Conf.replayLimit = double(replayLimit);
end
if ~isempty(allowInfiniteReplays)
    Conf.allowInfiniteReplays = logical(allowInfiniteReplays);
end
if ~isempty(eyeTrackerDoCalibration)
    Conf.eyeTrackerDoCalibration = logical(eyeTrackerDoCalibration);
end
if ~isempty(eyeTrackerAllowBetweenBlockCalibration)
    Conf.eyeTrackerAllowBetweenBlockCalibration = logical(eyeTrackerAllowBetweenBlockCalibration);
end
if ~isempty(eyeTrackerSendTriggerMessages)
    Conf.eyeTrackerSendTriggerMessages = logical(eyeTrackerSendTriggerMessages);
end
if ~isempty(eyeTrackerImageTransferEnabled)
    Conf.eyeTrackerImageTransferEnabled = logical(eyeTrackerImageTransferEnabled);
end
if Conf.ignoreEyeTracker
    Conf.fakeEyeTracker = false;
end
Conf.useFakeEyeTracker = Conf.fakeEyeTracker;
Conf.useEyelink = Conf.trackeye && ~Conf.ignoreEyeTracker && ~Conf.useFakeEyeTracker;
% Allow EyeLink setup/calibration even in fake-gaze mode when explicitly requested.
% Data flow:
%   runtime eye-tracker mode flags -> EyeLink session setup gate.
Conf.useEyelinkCalibrationWithFake = Conf.trackeye && ~Conf.ignoreEyeTracker && ...
    Conf.useFakeEyeTracker && Conf.eyeTrackerDoCalibration;
Conf.useEyelinkSession = Conf.useEyelink || Conf.useEyelinkCalibrationWithFake;
Conf.useGazeMonitoring = Conf.enableFixationAbort && (Conf.useEyelink || Conf.useFakeEyeTracker);
if ~isfinite(Conf.fakeGazeBreakRate) || Conf.fakeGazeBreakRate < 0 || Conf.fakeGazeBreakRate > 1
    error('fakeGazeBreakRate must be in [0, 1].');
end
if ~isfinite(Conf.fixWindowDeg) || Conf.fixWindowDeg <= 0
    error('fixWindowDeg must be finite and > 0.');
end
if ~isfinite(Conf.fixBreakToleranceFrames) || Conf.fixBreakToleranceFrames < 1
    error('fixBreakToleranceFrames must be >= 1.');
end
if ~isfinite(Conf.fixationWarningDurationSec) || Conf.fixationWarningDurationSec < 0
    error('fixationWarningDurationSec must be >= 0.');
end
if ~isfinite(Conf.replayLimit) || Conf.replayLimit < 0
    error('replayLimit must be >= 0.');
end
if ~ismember(Conf.eyeTrackerEdfPolicy, {'per_block'})
    error('Unsupported eyeTrackerEdfPolicy: %s', Conf.eyeTrackerEdfPolicy);
end
Conf.maxReplaysPerTrial = local_compute_max_replays_per_trial(Conf.replayLimit, Conf.allowInfiniteReplays);
if Conf.useEyelink
    fprintf('Eye tracker mode: REAL EyeLink enabled.\n');
elseif Conf.useFakeEyeTracker
    fprintf('Eye tracker mode: FAKE gaze samples enabled (rate=%.3f).\n', Conf.fakeGazeBreakRate);
    if Conf.useEyelinkCalibrationWithFake
        fprintf(['WARNING: fakeEyeTracker = 1, if this is a real run: stop, ' ...
            'disable fakeEyeTracker and re-run.\n']);
        fprintf('Eye tracker note: running EyeLink setup/calibration despite fake mode.\n');
    end
else
    fprintf('Eye tracker mode: disabled.\n');
end

% Debug visual overlays (disabled by default):
% - frame/timing text in the upper-left corner
% - full trajectory path line
Conf.debugShowFrameInfoOverlay = debugShowFrameInfoOverlay;
Conf.debugShowFullPathOverlay = debugShowFullPathOverlay;
Conf.debugOverlayTextColor = [240 240 240];
Conf.debugOverlayTextOffsetPx = [16 16];
Conf.debugPathLineColor = [95 95 95];
Conf.debugPathLineWidthPx = 2;

% Path-band terminal style controls (start and end edges).
% Data flow: schedule defaults -> optional workspace overrides -> validated
% runtime rendering behavior.
if ~exist('pathBandEntranceStyle', 'var') || isempty(pathBandEntranceStyle)
    pathBandEntranceStyle = scheduleCfg.pathBandEntranceStyle;
end
if ~exist('pathBandStraightBackshiftDotRadiusScale', 'var') || ...
        isempty(pathBandStraightBackshiftDotRadiusScale)
    pathBandStraightBackshiftDotRadiusScale = ...
        scheduleCfg.pathBandStraightBackshiftDotRadiusScale;
end
if ~exist('pathBandStraightEndForwardDotRadiusScale', 'var') || ...
        isempty(pathBandStraightEndForwardDotRadiusScale)
    if isfield(scheduleCfg, 'pathBandStraightEndForwardDotRadiusScale')
        pathBandStraightEndForwardDotRadiusScale = ...
            scheduleCfg.pathBandStraightEndForwardDotRadiusScale;
    else
        pathBandStraightEndForwardDotRadiusScale = 1.0;
    end
end
if ~exist('pathBandStraightPostStartWidthExtraDotRadiusScale', 'var') || ...
        isempty(pathBandStraightPostStartWidthExtraDotRadiusScale)
    if isfield(scheduleCfg, 'pathBandStraightPostStartWidthExtraDotRadiusScale')
        pathBandStraightPostStartWidthExtraDotRadiusScale = ...
            scheduleCfg.pathBandStraightPostStartWidthExtraDotRadiusScale;
    else
        pathBandStraightPostStartWidthExtraDotRadiusScale = 0.15;
    end
end
if ~exist('pathBandStraightStartBackshiftPixelOffset', 'var') || ...
        isempty(pathBandStraightStartBackshiftPixelOffset)
    if isfield(scheduleCfg, 'pathBandStraightStartBackshiftPixelOffset')
        pathBandStraightStartBackshiftPixelOffset = ...
            scheduleCfg.pathBandStraightStartBackshiftPixelOffset;
    else
        pathBandStraightStartBackshiftPixelOffset = -0.5;
    end
end
Conf.pathBandEntranceStyle = lower(strtrim(char(pathBandEntranceStyle)));
Conf.pathBandStraightBackshiftDotRadiusScale = double(pathBandStraightBackshiftDotRadiusScale);
Conf.pathBandStraightEndForwardDotRadiusScale = double(pathBandStraightEndForwardDotRadiusScale);
Conf.pathBandStraightPostStartWidthExtraDotRadiusScale = ...
    double(pathBandStraightPostStartWidthExtraDotRadiusScale);
Conf.pathBandStraightStartBackshiftPixelOffset = ...
    double(pathBandStraightStartBackshiftPixelOffset);
if ~ismember(Conf.pathBandEntranceStyle, {'round', 'straight'})
    error('pathBandEntranceStyle must be ''round'' or ''straight''.');
end
if ~isfinite(Conf.pathBandStraightBackshiftDotRadiusScale) || ...
        Conf.pathBandStraightBackshiftDotRadiusScale < 0
    error('pathBandStraightBackshiftDotRadiusScale must be a non-negative finite scalar.');
end
if ~isfinite(Conf.pathBandStraightEndForwardDotRadiusScale) || ...
        Conf.pathBandStraightEndForwardDotRadiusScale < 0
    error('pathBandStraightEndForwardDotRadiusScale must be a non-negative finite scalar.');
end
if ~isfinite(Conf.pathBandStraightPostStartWidthExtraDotRadiusScale) || ...
        Conf.pathBandStraightPostStartWidthExtraDotRadiusScale < 0
    error('pathBandStraightPostStartWidthExtraDotRadiusScale must be a non-negative finite scalar.');
end
if ~isfinite(Conf.pathBandStraightStartBackshiftPixelOffset)
    error('pathBandStraightStartBackshiftPixelOffset must be a finite scalar.');
end

if ~exist('question_text', 'var') || isempty(question_text)
    question_text = scheduleCfg.catchQuestionText;
end
Conf.catchPromptText = char(question_text);
Conf.catchResponseTimeoutSec = double(scheduleCfg.catchQuestionTimeoutSec);
Conf.catchResponseYesCode = double(scheduleCfg.catchResponseYesCode);
Conf.catchResponseNoCode = double(scheduleCfg.catchResponseNoCode);
Conf.catchType1DisappearRangeSec = double(scheduleCfg.catchType1DisappearRangeSec);
Conf.catchType1InvisibleDurationSec = double(scheduleCfg.catchType1InvisibleDurationSec);
Conf.catchType1ChangedPathProbability = double(scheduleCfg.catchType1ChangedPathProbability);
Conf.catchKeyYes = {'RightArrow', '8*', 'y', 'Y'};
Conf.catchKeyNo = {'LeftArrow', '1!', 'n', 'N'};
Conf.catchNoLabel = 'NO';
Conf.catchYesLabel = 'YES';
Conf.catchFeedbackDurationSec = 0.500;
Conf.catchFeedbackCorrectColor = [0 220 0];
Conf.catchFeedbackIncorrectColor = [220 0 0];
Conf.catchTimeoutText = 'Too slow';
Conf.catchQuestionTopYFrac = 0.12;
Conf.catchQuestionBottomYFrac = 0.88;
Conf.catchBottomLeftXFrac = 0.18;
Conf.catchBottomRightXFrac = 0.82;
Conf.megButtonMask = uint32(2^0 + 2^1 + 2^2 + 2^3);
Conf.megButtonYesValue = uint32(8); % blue
Conf.megButtonNoValue = uint32(1); % red
Conf.enableKeyboardDebugResponse = true;

% Start gate and run-transition messages (v18 rescueTraject).
% Data flow: schedule config constants -> runtime text, toggle, and timing controls.
Conf.startMessageEnabled = logical(scheduleCfg.startMessageEnabled);
Conf.startMessageText = char(scheduleCfg.startMessageText);
Conf.startMessageDurationSec = double(scheduleCfg.startMessageDurationSec);
Conf.startMessageAcceptKeyNames = {'1!', '8*', '1', '8', 'KP_1', 'KP_8', 'num_1', 'num_8'};
Conf.startMessageAcceptMegValues = uint32([1, 8]);
Conf.run1TransitionMessageEnabled = logical(scheduleCfg.run1TransitionMessageEnabled);
Conf.run1TransitionMessageText = char(scheduleCfg.run1TransitionMessageText);
Conf.run1TransitionMessageDurationSec = double(scheduleCfg.run1TransitionMessageDurationSec);
Conf.run1TransitionCountdownEnabled = logical(local_schedule_get_or_default( ...
    scheduleCfg, 'run1TransitionCountdownEnabled', true));
if ~isempty(run1TransitionCountdownEnabled)
    Conf.run1TransitionCountdownEnabled = logical(run1TransitionCountdownEnabled);
end
Conf.endOfBlockMessageEnabled = logical(scheduleCfg.endOfBlockMessageEnabled);
Conf.endOfBlockMessageTextTemplate = char(scheduleCfg.endOfBlockMessageTextTemplate);
Conf.endOfBlockMessageDurationSec = double(scheduleCfg.endOfBlockMessageDurationSec);
Conf.betweenBlockCalibrationChoiceEnabled = logical(local_schedule_get_or_default( ...
    scheduleCfg, 'betweenBlockCalibrationChoiceEnabled', true));
Conf.betweenBlockCalibrationChoiceText = char(local_schedule_get_or_default( ...
    scheduleCfg, 'betweenBlockCalibrationChoiceText', ...
    ['Calibration choice.\n\n' ...
     'Press Red / key 1 to continue without calibration.\n' ...
     'Press Blue / key 8 to run calibration now.\n\n' ...
     'Then the next block starts.']));
Conf.betweenBlockCalibrationChoiceDurationSec = double(local_schedule_get_or_default( ...
    scheduleCfg, 'betweenBlockCalibrationChoiceDurationSec', 0.0));
Conf.calibrationChoiceSkipMegValues = uint32(1);
Conf.calibrationChoiceRunMegValues = uint32(8);
Conf.calibrationChoiceSkipKeyNames = {'1!', '1', 'KP_1', 'num_1'};
Conf.calibrationChoiceRunKeyNames = {'8*', '8', 'KP_8', 'num_8'};
Conf.finalMessageEnabled = logical(scheduleCfg.finalMessageEnabled);
Conf.finalMessageText = char(scheduleCfg.finalMessageText);
Conf.finalMessageDurationSec = double(scheduleCfg.finalMessageDurationSec);
Conf.postMessageFixationEnabled = logical(scheduleCfg.postMessageFixationEnabled);
Conf.postMessageFixationDurationSec = double(scheduleCfg.postMessageFixationDurationSec);
Conf.messageResponseLockoutSec = double(scheduleCfg.messageResponseLockoutSec);
Conf.practiceRunPercentageTrials = double(scheduleCfg.practiceRunPercentageTrials);
if exist('practiceRunPercentageTrials', 'var') && ~isempty(practiceRunPercentageTrials)
    Conf.practiceRunPercentageTrials = double(practiceRunPercentageTrials);
end
if ~isfinite(Conf.practiceRunPercentageTrials) || Conf.practiceRunPercentageTrials <= 0 || ...
        Conf.practiceRunPercentageTrials > 100
    error('practiceRunPercentageTrials must be in (0, 100].');
end

% Timing diagnostics controls (v17):
% schedule defaults -> optional workspace overrides -> validated runtime.
Conf.skipSyncTestsWhenMEG = logical(scheduleCfg.skipSyncTestsWhenMEG);
Conf.skipSyncTestsInDebug = logical(scheduleCfg.skipSyncTestsInDebug);
Conf.flipMissedWarnThresholdCount = max(1, round(double(scheduleCfg.flipMissedWarnThresholdCount)));
Conf.flipMissedAbortEnabled = logical(scheduleCfg.flipMissedAbortEnabled);
Conf.flipMissedAbortThresholdCount = max(1, round(double(scheduleCfg.flipMissedAbortThresholdCount)));
Conf.flipMissedEpsilonSec = max(0, double(scheduleCfg.flipMissedEpsilonSec));

if ~isempty(flipMissedWarnThresholdCount)
    Conf.flipMissedWarnThresholdCount = max(1, round(double(flipMissedWarnThresholdCount)));
end
if ~isempty(flipMissedAbortEnabled)
    Conf.flipMissedAbortEnabled = logical(flipMissedAbortEnabled);
end
if ~isempty(flipMissedAbortThresholdCount)
    Conf.flipMissedAbortThresholdCount = max(1, round(double(flipMissedAbortThresholdCount)));
end

if isempty(debugSkipSyncTests)
    Conf.debugSkipSyncTestsOverride = [];
else
    Conf.debugSkipSyncTestsOverride = logical(debugSkipSyncTests);
end

% Trigger map for this occlusion paradigm.
Conf.trigger.conditionAlwaysVisible = 31;
Conf.trigger.conditionOccludedNondeviant = 41;
Conf.trigger.conditionOccludedDeviant = 51;
Conf.trigger.occlusionStart = 111;
Conf.trigger.occlusionComplete = 112;
Conf.trigger.occlusionEndStart = 114;
Conf.trigger.occlusionEndComplete = 115;
Conf.trigger.catch1Disappear = 116;
Conf.trigger.catch1Reappear = 117;
Conf.trigger.catchQuestionStart = 118;
Conf.trigger.catchQuestionEndCorrect = 201;
Conf.trigger.catchQuestionEndIncorrect = 202;
Conf.trigger.catchQuestionEndTimeout = 113;
Conf.trigger.gazeBreak = 150;
Conf.trigger.replayStart = 151;
Conf.trigger.escapeTermination = double(local_schedule_get_or_default( ...
    scheduleCfg, 'triggerEscTermination', 152));
Conf.sequenceIdentityTriggerFrame = 12;

if ~ismember(lower(Conf.occlusionMode), {'pathband', 'alpha'})
    error('Conf.occlusionMode must be ''pathband'' or ''alpha''.');
end

%% Collect user input and resolve files
% Data flow: dialog -> subject/block/display geometry -> input/output file names.
if ~(exist('iSub', 'var') && exist('iBlock', 'var') && exist('viewingDistanceMm', 'var') && exist('practice_mode', 'var'))
    prompt = {'Subject Number:', 'Block Number:', 'Viewing Distance (mm):', 'Practice mode (0/1):'};
    dlgtitle = 'Occlusion v18 rescueTraject (practice by run) input';
    dims = [1 50];
    definput = {'70', '1', '1000', '0'};
    userInput = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(userInput)
        disp('Canceled by user.');
        return;
    end

    iSub = str2double(userInput{1});
    iBlock = str2double(userInput{2});
    viewingDistanceMm = str2double(userInput{3});
    practice_mode = str2double(userInput{4});
end

if ~exist('practice_mode', 'var') || isempty(practice_mode)
    practice_mode = 0;
end

if any(isnan([iSub, iBlock, viewingDistanceMm, practice_mode])) || viewingDistanceMm <= 0
    error('Invalid subject/block/viewing-distance/practice input.');
end
iSub = round(iSub);
iBlock = round(iBlock);
practice_mode = round(practice_mode);
if ~ismember(practice_mode, [0, 1])
    error('Practice mode must be 0 or 1.');
end
Conf.practiceModeEnabled = logical(practice_mode);

if ~isempty(debug_runNumber)
    debug_runNumber = round(double(debug_runNumber));
    if ~isscalar(debug_runNumber) || ~ismember(debug_runNumber, 1:scheduleCfg.runsPerBlock)
        error('debug_runNumber must be empty or one of 1..%d.', scheduleCfg.runsPerBlock);
    end
end

rng(iSub * 100 + iBlock);

rootDir = '.';
inputDir = fullfile(rootDir, 'input_files');
outputDir = fullfile(rootDir, 'output_files');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

inputFile = fullfile(inputDir, sprintf(scheduleCfg.inputFilePattern, iSub));
trialStructFile = fullfile(inputDir, sprintf(scheduleCfg.outputFilePattern, iSub));
if ~isfile(inputFile)
    error('Input file not found: %s', inputFile);
end
if ~isfile(trialStructFile)
    error(['Required schedule file not found: %s\n' ...
        'Generate it with CreateInputFiles_v21_rescueTraject.m first.'], trialStructFile);
end

outputBase = sprintf('%s_SUB%02d_BLOCK%02d', Conf.experimentName, iSub, iBlock);
outputMat = fullfile(outputDir, [outputBase '.mat']);

if exist(outputMat, 'file') == 2
    error(['Output already exists; refusing to overwrite: %s\n' ...
        'if you intended to run from a specific block, set iBlock variable'], outputMat);
end

%% Load trials and validate metadata
% Data flow: input MAT -> one-dot trial structs -> condition grouping + metadata checks.
Dat = load(inputFile, 'xySeqs', 'Cfg');
if ~isfield(Dat, 'xySeqs') || ~isfield(Dat, 'Cfg')
    error('Input file must contain xySeqs and Cfg: %s', inputFile);
end

xyAll = Dat.xySeqs(:);
if isempty(xyAll)
    error('No trials found in xySeqs.');
end

requiredFields = { ...
    'xy', 'condition_label', 'deviance_frame', 'occlusion_enabled', ...
    'occlusion_start_frame', 'occlusion_complete_frame', ...
    'occlusion_end_frame', 'occlusion_end_complete_frame', ...
    'pathband_pre_xy', 'pathband_post_xy', 'pathband_width_deg'};

for iTrial = 1:numel(xyAll)
    for iField = 1:numel(requiredFields)
        if ~isfield(xyAll(iTrial), requiredFields{iField})
            error('Trial %d missing required field: %s', iTrial, requiredFields{iField});
        end
    end
    if size(xyAll(iTrial).xy, 2) ~= 2
        error('Trial %d xy must be one-dot [x y].', iTrial);
    end
end

framesPerTrial = size(xyAll(1).xy, 1);
for iTrial = 2:numel(xyAll)
    if size(xyAll(iTrial).xy, 1) ~= framesPerTrial
        error('All trials must have identical frame count in this runtime.');
    end
end

% Build per-source sequence trigger codes from input metadata.
% This is the value sent as identity trigger at frame Conf.sequenceIdentityTriggerFrame.
sequenceTriggerBySource = local_extract_source_sequence_trigger_codes(xyAll);
Conf.sequenceTriggerRange = 1:max(sequenceTriggerBySource);
local_validate_sequence_trigger_range(Conf.sequenceTriggerRange, Conf.trigger);
if ~isscalar(Conf.sequenceIdentityTriggerFrame) || Conf.sequenceIdentityTriggerFrame < 1 || ...
        Conf.sequenceIdentityTriggerFrame > framesPerTrial
    error('Conf.sequenceIdentityTriggerFrame must be within 1..%d frames.', framesPerTrial);
end

% Condition indexing by explicit labels.
labels = arrayfun(@(s) string(s.condition_label), xyAll, 'UniformOutput', true);
idxAlways = find(labels == "always_visible");
idxOccNondev = find(labels == "occluded_nondeviant");
idxOccDev = find(labels == "occluded_deviant");

if isempty(idxAlways) || isempty(idxOccNondev) || isempty(idxOccDev)
    error('Input must contain always_visible, occluded_nondeviant, occluded_deviant labels.');
end

nPerCond = min([numel(idxAlways), numel(idxOccNondev), numel(idxOccDev)]);
idxAlways = idxAlways(1:nPerCond);
idxOccNondev = idxOccNondev(1:nPerCond);
idxOccDev = idxOccDev(1:nPerCond);

% Required scheduling path:
%   Load TrialOrder/CatchPlan from v20 schedule file and resolve selected blocks/runs.
tsData = load(trialStructFile, 'TrialOrder', 'CatchPlan', 'Schedule');
if ~isfield(tsData, 'TrialOrder')
    error('TrialOrder missing in required schedule file: %s', trialStructFile);
end

if ndims(tsData.TrialOrder) ~= 3
    error('TrialOrder in %s must be 3D [block x trial x run] for v18 rescueTraject mode.', ...
        trialStructFile);
end

nBlocksAvailable = size(tsData.TrialOrder, 1);
runsPerBlockAvailable = size(tsData.TrialOrder, 3);
if iBlock < 1 || iBlock > nBlocksAvailable
    error('Block number %d is outside available range 1..%d.', iBlock, nBlocksAvailable);
end
if runsPerBlockAvailable ~= scheduleCfg.runsPerBlock
    error('TrialOrder run dimension (%d) does not match expected runsPerBlock (%d).', ...
        runsPerBlockAvailable, scheduleCfg.runsPerBlock);
end

if isfield(tsData, 'CatchPlan')
    catchPlan = local_normalize_catch_plan(tsData.CatchPlan, tsData.TrialOrder);
else
    warning('CatchPlan missing in %s; running with no catches.', trialStructFile);
    catchPlan = local_empty_catch_plan(tsData.TrialOrder);
end

if isfield(tsData, 'Schedule')
    if isfield(tsData.Schedule, 'numBlocks')
        scheduleBlocks = double(tsData.Schedule.numBlocks);
        if scheduleBlocks ~= nBlocksAvailable
            warning('Schedule.numBlocks (%d) does not match TrialOrder block size (%d). Using TrialOrder dimensions.', ...
                scheduleBlocks, nBlocksAvailable);
        end
    end
    if isfield(tsData.Schedule, 'runsPerBlock')
        scheduleRuns = double(tsData.Schedule.runsPerBlock);
        if scheduleRuns ~= runsPerBlockAvailable
            warning('Schedule.runsPerBlock (%d) does not match TrialOrder run size (%d). Using TrialOrder dimensions.', ...
                scheduleRuns, runsPerBlockAvailable);
        end
    end
end

if isempty(debug_runNumber)
    runsToExecute = 1:runsPerBlockAvailable;
else
    runsToExecute = debug_runNumber;
end

blocksToExecute = iBlock:nBlocksAvailable;

if isempty(debug_runNumber)
    fprintf('Using v21 TrialOrder from %s (start block %d, runs 1..%d).\n', ...
        trialStructFile, iBlock, runsPerBlockAvailable);
else
    fprintf('Using v21 TrialOrder from %s (start block %d, debug run %d only).\n', ...
        trialStructFile, iBlock, debug_runNumber);
end

if logical(Conf.practiceModeEnabled)
    fprintf(['Practice mode ENABLED: keeping %.2f%% of each run and forcing ' ...
        'selected trials to catches (run1->type1, runs2/3->type2).\n'], ...
        Conf.practiceRunPercentageTrials);
else
    fprintf('Practice mode disabled: running full scheduled trials.\n');
end

if ~exist('dryRunSimulatedBreakKeys', 'var') || isempty(dryRunSimulatedBreakKeys)
    dryRunSimulatedBreakKeys = {};
end

%% Dry-run validation path (no Psychtoolbox windows)
% Data flow: selected blocks/runs -> schedule checks + placeholder block saves.
if logical(dryRunValidateScheduleOnly)
    for iBlockIdx = 1:numel(blocksToExecute)
        currentBlock = blocksToExecute(iBlockIdx);
        outputBase = sprintf('%s_SUB%02d_BLOCK%02d', Conf.experimentName, iSub, currentBlock);
        outputMat = fullfile(outputDir, [outputBase '.mat']);
        if exist(outputMat, 'file') == 2
            error(['Output already exists; refusing to overwrite: %s\n' ...
                'if you intended to run from a specific block, set iBlock variable'], outputMat);
        end
        if Conf.useEyelink
            Conf = local_eyelink_open_block_file(Conf, iSub, currentBlock, outputDir);
        end

        runResults = repmat(struct( ...
            'subject_num', iSub, ...
            'block_num', currentBlock, ...
            'run_num', [], ...
            'trialOrder', [], ...
            'trialOrderBase', [], ...
            'trialIsReplay', [], ...
            'trialReplayInstance', [], ...
            'itiRangeSec', Conf.ITIRangeSec, ...
            'itiPrecomputeSeed', [], ...
            'itiPrecomputedSec', [], ...
            'itiUsedSec', [], ...
            'itiSourceTrialIndex', [], ...
            'output', [], ...
            'nTrialsPlanned', [], ...
            'nTrialsCompleted', [], ...
            'nCatchPlanned', [], ...
            'nCatchCompleted', [], ...
            'nReplaysQueued', 0, ...
            'aborted', false, ...
            'abort_reason', '', ...
            'debugTriggerLog', [], ...
            'debugTriggerTrial', [], ...
            'debugTriggerFrame', [], ...
            'debugTriggerSec', [], ...
            'debugTriggerLabel', {{}}, ...
            'debugTriggerFile', '', ...
            'eyeTrackerEdfFile', '', ...
            'eyeTrackerLocalEdfPath', ''), 1, numel(runsToExecute));

        for iRunIdx = 1:numel(runsToExecute)
            currentRun = runsToExecute(iRunIdx);
            trialOrderBase = local_resolve_block_run_order(tsData.TrialOrder, currentBlock, currentRun, numel(xyAll));
            if isempty(trialOrderBase)
                error('Invalid TrialOrder slice for block %d run %d in %s.', ...
                    currentBlock, currentRun, trialStructFile);
            end

            catchTypeBase = local_resolve_block_run_catch( ...
                catchPlan.typeCode, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), 0);
            catchExpectedBase = local_resolve_block_run_catch( ...
                catchPlan.expectedResponseCode, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), 0);
            catchBranchBase = local_resolve_block_run_catch( ...
                catchPlan.branchChangedPath, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);
            catchDisappearBase = local_resolve_block_run_catch( ...
                catchPlan.disappearFrame, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);
            catchReappearBase = local_resolve_block_run_catch( ...
                catchPlan.reappearFrame, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);
            catchAltSourceBase = local_resolve_block_run_catch( ...
                catchPlan.altSourceTrialIndex, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);

            [trialOrder, catchTypeVec] = local_apply_practice_mode_to_run( ...
                trialOrderBase, catchTypeBase, catchExpectedBase, catchBranchBase, ...
                catchDisappearBase, catchReappearBase, catchAltSourceBase, ...
                currentRun, xyAll, double(Dat.Cfg.fps), framesPerTrial, ...
                Conf.practiceModeEnabled, Conf.practiceRunPercentageTrials, Conf);

            runResults(iRunIdx).run_num = currentRun;
            runResults(iRunIdx).trialOrder = trialOrder;
            [itiPrecomputedSec, itiSeed] = local_precompute_run_iti_seconds( ...
                iSub, currentBlock, currentRun, numel(trialOrder), Conf.ITIRangeSec);
            itiSourceTrialIndex = 1:numel(trialOrder);
            itiUsedSec = itiPrecomputedSec;
            local_log_iti_precompute('Dry-run', iSub, currentBlock, currentRun, ...
                itiSeed, Conf.ITIRangeSec, itiPrecomputedSec);
            runResults(iRunIdx).itiRangeSec = Conf.ITIRangeSec;
            runResults(iRunIdx).itiPrecomputeSeed = itiSeed;
            runResults(iRunIdx).itiPrecomputedSec = itiPrecomputedSec;
            runResults(iRunIdx).itiUsedSec = itiUsedSec;
            runResults(iRunIdx).itiSourceTrialIndex = itiSourceTrialIndex;
            runResults(iRunIdx).output = struct([]);
            runResults(iRunIdx).nTrialsPlanned = numel(trialOrder);
            runResults(iRunIdx).nTrialsCompleted = numel(trialOrder);
            runResults(iRunIdx).nCatchPlanned = sum(catchTypeVec > 0);
            runResults(iRunIdx).nCatchCompleted = sum(catchTypeVec > 0);
            runDotColor = local_dot_color_for_run(Conf, iSub, currentRun);
            fprintf('Dry-run schedule validation mode: subject %02d block %d run %d -> %d trials.\n', ...
                iSub, currentBlock, currentRun, numel(trialOrder));
            fprintf('Dry-run run-color mapping: subject %02d run %d -> [%d %d %d] (enabled=%d)\n', ...
                iSub, currentRun, runDotColor(1), runDotColor(2), runDotColor(3), Conf.runColorCueEnabled);
        end

        save(outputMat, 'runResults', 'Conf', 'iSub', 'currentBlock', 'viewingDistanceMm', ...
            'inputFile', 'trialStructFile', 'runsToExecute', 'debug_runNumber', ...
            'dryRunValidateScheduleOnly', 'catchPlan');
        fprintf('Dry-run block output saved: %s\n', outputMat);

        if currentBlock < nBlocksAvailable
            action = local_resolve_dryrun_break_action(dryRunSimulatedBreakKeys, iBlockIdx);
            fprintf('Dry-run end-of-block action after block %d: %s\n', currentBlock, action);
            if strcmp(action, 'quit')
                fprintf('Dry-run session terminated by simulated q after block %d.\n', currentBlock);
                break;
            end
        else
            fprintf('Dry-run final message: You reached the end of the experiment! Thank you for participating.\n');
        end
    end
    return;
end

%% Screen setup and geometry transforms
% Data flow: display hardware + viewing distance + Cfg rect sizes -> pixel-space draw geometry.
try
    %% MEG/DataPixx setup (optional)
    % Data flow: Conf.MEG flag -> DataPixx initialization and trigger readiness.
    if Conf.MEG
        Datapixx('Open');
        Datapixx('SetVideoMode', 0);
        Datapixx('StopAllSchedules');
        Datapixx('EnableDinDebounce');
        Datapixx('SetDinLog');
        Datapixx('StopDinLog');
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
    end
    
    skipSyncTests = local_resolve_skip_sync_tests(Conf);
    Screen('Preference', 'SkipSyncTests', double(skipSyncTests));
    if skipSyncTests
        fprintf('PTB sync policy: SkipSyncTests=1 (debug/development mode).\n');
    else
        fprintf('PTB sync policy: SkipSyncTests=0 (timing-safe mode).\n');
    end
    screens = Screen('Screens');
    screenNumber = max(screens);
    [window, windowRect] = Screen('OpenWindow', screenNumber, Conf.background);
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    KbName('UnifyKeyNames');
    Conf.escapeKeyCode = KbName('ESCAPE');
    Conf.catchYesKeyCodes = local_keynames_to_codes(Conf.catchKeyYes);
    Conf.catchNoKeyCodes = local_keynames_to_codes(Conf.catchKeyNo);
    Conf.startMessageAcceptKeyCodes = local_keynames_to_codes(Conf.startMessageAcceptKeyNames);
    Conf.calibrationChoiceSkipKeyCodes = local_keynames_to_codes(Conf.calibrationChoiceSkipKeyNames);
    Conf.calibrationChoiceRunKeyCodes = local_keynames_to_codes(Conf.calibrationChoiceRunKeyNames);

    ifi = Screen('GetFlipInterval', window);
    actualFrameRate = 1 / ifi;
    Conf.refrate = Dat.Cfg.fps;
    if abs(actualFrameRate - Conf.refrate) > max(1, 0.05 * Conf.refrate)
        fprintf('WARNING: measured %.2f Hz differs from Cfg.fps %.2f Hz\n', actualFrameRate, Conf.refrate);
    end
    if actualFrameRate > 0
        playbackDurationScale = Conf.refrate / actualFrameRate;
        if abs(playbackDurationScale - 1) > 0.05
            fprintf(['NOTE: stimuli were generated for %.2f Hz but display is %.2f Hz.\n' ...
                '      Motion playback duration is scaled by ~%.2fx on this display.\n'], ...
                Conf.refrate, actualFrameRate, playbackDurationScale);
        end
    end

    [screenWidthPx, ~] = Screen('WindowSize', window);
    [monWidthMm, ~] = Screen('DisplaySize', screenNumber);
    ppd = pi * screenWidthPx / atan(monWidthMm / viewingDistanceMm / 2) / 360;
    Conf.fixWindowPx = Conf.fixWindowDeg * ppd;

    rectSizeDeg = double(Dat.Cfg.rectSize(:)');
    if numel(rectSizeDeg) ~= 2
        error('Cfg.rectSize must be [width height].');
    end

    rectSizePx = floor(rectSizeDeg * ppd);
    [xCenter, yCenter] = RectCenter(windowRect);
    rectCoords = [xCenter - rectSizePx(1)/2, yCenter - rectSizePx(2)/2, ...
        xCenter + rectSizePx(1)/2, yCenter + rectSizePx(2)/2];

    dotRadiusPx = max(2, floor((double(Dat.Cfg.dot_w) * ppd) / 2));
    dotSizePx = dotRadiusPx * 2;
    fixSizePx = max(4, floor(Conf.fixSizeDeg * ppd));

    % Initialize EyeLink session (optional) after the PTB window is ready.
    % Data flow:
    %   runtime Conf eye-tracker flags + PTB window geometry -> EyeLink
    %   session state used by trial recording and/or calibration plus cleanup.
    Conf.eyeTrackerInitialized = false;
    Conf.eyeTrackerEdfFile = '';
    Conf.eyeTrackerLocalEdfPath = '';
    if Conf.useEyelinkSession
        Conf = local_setup_eyelink_session(Conf, window, windowRect, screenNumber, iSub, iBlock, outputDir);
    end

    % Convert all trajectories to pixel coordinates once.
    stimuliPx = zeros(numel(xyAll), framesPerTrial, 2);
    for iTrial = 1:numel(xyAll)
        xy = double(xyAll(iTrial).xy);
        stimuliPx(iTrial, :, 1) = xy(:, 1) * ppd + rectCoords(1);
        stimuliPx(iTrial, :, 2) = xy(:, 2) * ppd + rectCoords(2);
    end


    %% Execute selected blocks
    % Data flow: selected blocks -> run loops -> per-block save + break/continue flow.
    HideCursor;
    % Priority policy:
    %   disablePriorityBoost=true is the practice-safe default. This avoids
    %   Mach time MEX crashes on unsupported/codesign-blocked PTB setups.
    if disablePriorityBoost
        fprintf('PTB priority policy: disabled (practice safety mode).\n');
        try
            Priority(0);
        catch
        end
    else
        try
            Priority(MaxPriority(window));
        catch MEPriority
            warning(['Priority boost failed; continuing at normal priority. ' ...
                'Error was: %s'], MEPriority.message);
            try
                Priority(0);
            catch
            end
        end
    end

    stopExperimentAfterBlock = false;
    stopReason = '';

    % Start-of-experiment gate: show welcome message and wait for 1/8.
    startAction = local_show_start_message_gate(window, Conf);
    if strcmp(startAction, 'abort')
        stopExperimentAfterBlock = true;
        stopReason = 'escape_on_start_message';
        local_emit_escape_termination_trigger(Conf, iBlock, 'escape_before_block');
    elseif strcmp(startAction, 'continue')
        abortedDuringFixation = local_show_post_message_fixation(window, Conf, xCenter, yCenter, fixSizePx);
        if abortedDuringFixation
            stopExperimentAfterBlock = true;
            stopReason = 'escape_during_post_message_fixation';
            local_emit_escape_termination_trigger(Conf, iBlock, 'escape_before_block');
        end
    end

    % Session-level catch accuracy accumulators (aggregated by run index).
    % Data flow:
    %   per-trial catch correctness -> per-run totals -> end-of-session summary.
    runCatchScored = zeros(1, runsPerBlockAvailable);
    runCatchAnswered = zeros(1, runsPerBlockAvailable);
    runCatchCorrect = zeros(1, runsPerBlockAvailable);

    for iBlockIdx = 1:numel(blocksToExecute)
        if stopExperimentAfterBlock
            break;
        end
        currentBlock = blocksToExecute(iBlockIdx);
        outputBase = sprintf('%s_SUB%02d_BLOCK%02d', Conf.experimentName, iSub, currentBlock);
        outputMat = fullfile(outputDir, [outputBase '.mat']);
        if exist(outputMat, 'file') == 2
            error(['Output already exists; refusing to overwrite: %s\n' ...
                'if you intended to run from a specific block, set iBlock variable'], outputMat);
        end

        runTrialOrders = cell(1, numel(runsToExecute));
        for iRunIdx = 1:numel(runsToExecute)
            currentRun = runsToExecute(iRunIdx);
            runTrialOrders{iRunIdx} = local_resolve_block_run_order( ...
                tsData.TrialOrder, currentBlock, currentRun, numel(xyAll));
            if isempty(runTrialOrders{iRunIdx})
                error('Invalid TrialOrder slice for block %d run %d in %s.', ...
                    currentBlock, currentRun, trialStructFile);
            end
        end

        runResults = repmat(struct( ...
            'subject_num', iSub, ...
            'block_num', currentBlock, ...
            'run_num', [], ...
            'trialOrder', [], ...
            'trialOrderBase', [], ...
            'trialIsReplay', [], ...
            'trialReplayInstance', [], ...
            'itiRangeSec', Conf.ITIRangeSec, ...
            'itiPrecomputeSeed', [], ...
            'itiPrecomputedSec', [], ...
            'itiUsedSec', [], ...
            'itiSourceTrialIndex', [], ...
            'output', [], ...
            'nTrialsPlanned', [], ...
            'nTrialsCompleted', [], ...
            'nCatchPlanned', [], ...
            'nCatchCompleted', [], ...
            'nReplaysQueued', 0, ...
            'aborted', false, ...
            'abort_reason', '', ...
            'debugTriggerLog', [], ...
            'debugTriggerTrial', [], ...
            'debugTriggerFrame', [], ...
            'debugTriggerSec', [], ...
            'debugTriggerLabel', {{}}, ...
            'debugTriggerFile', '', ...
            'eyeTrackerEdfFile', '', ...
            'eyeTrackerLocalEdfPath', ''), 1, numel(runsToExecute));

        % Precompute all run schedules and ITI vectors once at block start.
        % Data flow:
        %   TrialOrder/CatchPlan + practice policy -> run-specific base schedule
        %   and catch vectors -> deterministic ITI vector per run (before trial loop).
        runTrialOrdersPlanned = cell(1, numel(runsToExecute));
        runCatchTypeVec = cell(1, numel(runsToExecute));
        runCatchExpectedVec = cell(1, numel(runsToExecute));
        runCatchBranchVec = cell(1, numel(runsToExecute));
        runCatchDisappearVec = cell(1, numel(runsToExecute));
        runCatchReappearVec = cell(1, numel(runsToExecute));
        runCatchAltSourceVec = cell(1, numel(runsToExecute));
        runItiPrecomputedSec = cell(1, numel(runsToExecute));
        runItiPrecomputeSeed = nan(1, numel(runsToExecute));

        for iRunIdx = 1:numel(runsToExecute)
            currentRun = runsToExecute(iRunIdx);
            trialOrderBase = runTrialOrders{iRunIdx};
            catchTypeBase = local_resolve_block_run_catch( ...
                catchPlan.typeCode, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), 0);
            catchExpectedBase = local_resolve_block_run_catch( ...
                catchPlan.expectedResponseCode, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), 0);
            catchBranchBase = local_resolve_block_run_catch( ...
                catchPlan.branchChangedPath, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);
            catchDisappearBase = local_resolve_block_run_catch( ...
                catchPlan.disappearFrame, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);
            catchReappearBase = local_resolve_block_run_catch( ...
                catchPlan.reappearFrame, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);
            catchAltSourceBase = local_resolve_block_run_catch( ...
                catchPlan.altSourceTrialIndex, tsData.TrialOrder, currentBlock, currentRun, numel(xyAll), NaN);

            [trialOrder, catchTypeOut, catchExpectedOut, catchBranchOut, catchDisappearOut, catchReappearOut, catchAltSourceOut] = ...
                local_apply_practice_mode_to_run( ...
                trialOrderBase, catchTypeBase, catchExpectedBase, catchBranchBase, ...
                catchDisappearBase, catchReappearBase, catchAltSourceBase, ...
                currentRun, xyAll, double(Dat.Cfg.fps), framesPerTrial, ...
                Conf.practiceModeEnabled, Conf.practiceRunPercentageTrials, Conf);

            runTrialOrdersPlanned{iRunIdx} = trialOrder;
            runCatchTypeVec{iRunIdx} = catchTypeOut;
            runCatchExpectedVec{iRunIdx} = catchExpectedOut;
            runCatchBranchVec{iRunIdx} = catchBranchOut;
            runCatchDisappearVec{iRunIdx} = catchDisappearOut;
            runCatchReappearVec{iRunIdx} = catchReappearOut;
            runCatchAltSourceVec{iRunIdx} = catchAltSourceOut;
            [runItiPrecomputedSec{iRunIdx}, runItiPrecomputeSeed(iRunIdx)] = ...
                local_precompute_run_iti_seconds( ...
                iSub, currentBlock, currentRun, numel(trialOrder), Conf.ITIRangeSec);
            local_log_iti_precompute('Run-precompute', iSub, currentBlock, currentRun, ...
                runItiPrecomputeSeed(iRunIdx), Conf.ITIRangeSec, runItiPrecomputedSec{iRunIdx});
        end

        blockAborted = false;
        blockAbortReason = '';

        for iRunIdx = 1:numel(runsToExecute)
            % Show the transition message only when entering runs after run 1.
            if iRunIdx > 1
                previousRun = runsToExecute(iRunIdx - 1);
                if previousRun == 1 && logical(Conf.run1TransitionMessageEnabled)
                    abortedDuringTransition = local_show_timed_transition_message( ...
                        window, Conf, Conf.run1TransitionMessageText, Conf.run1TransitionMessageDurationSec);
                    if abortedDuringTransition
                        blockAborted = true;
                        blockAbortReason = 'escape_key';
                        break;
                    end

                    abortedDuringFixation = local_show_post_message_fixation(window, Conf, xCenter, yCenter, fixSizePx);
                    if abortedDuringFixation
                        blockAborted = true;
                        blockAbortReason = 'escape_key';
                        break;
                    end
                end
            end

            currentRun = runsToExecute(iRunIdx);
            trialOrderBase = runTrialOrders{iRunIdx};
            runDotColor = local_dot_color_for_run(Conf, iSub, currentRun);
            trialOrder = runTrialOrdersPlanned{iRunIdx};
            catchTypeVec = runCatchTypeVec{iRunIdx};
            catchExpectedVec = runCatchExpectedVec{iRunIdx};
            catchBranchVec = runCatchBranchVec{iRunIdx};
            catchDisappearVec = runCatchDisappearVec{iRunIdx};
            catchReappearVec = runCatchReappearVec{iRunIdx};
            catchAltSourceVec = runCatchAltSourceVec{iRunIdx};
            nTrials = numel(trialOrder);
            itiPrecomputedSec = runItiPrecomputedSec{iRunIdx};
            itiSeed = runItiPrecomputeSeed(iRunIdx);

            fprintf('Starting block %d run %d (%d trials).\n', currentBlock, currentRun, nTrials);

            % Run-local trigger traces.
            debugTriggerLog = [];
            debugTriggerTrial = [];
            debugTriggerFrame = [];
            debugTriggerSec = [];
            debugTriggerLabel = {};

            % Run-local trial schedule:
            % base schedule plus optional replay queue entries appended when
            % fixation breaks are detected.
            trialSchedule = trialOrder(:)';
            trialIsReplay = false(1, numel(trialSchedule));
            trialReplayInstance = zeros(1, numel(trialSchedule));
            itiSourceTrialIndexSchedule = 1:numel(trialSchedule);
            replayCountsBySource = zeros(1, numel(xyAll));
            replayStartTriggerSent = false;
            replayQueuedCount = 0;

            % Run-local output struct template.
            outputTemplate = struct( ...
                'subject_num', iSub, ...
                'block_num', currentBlock, ...
                'run_num', currentRun, ...
                'trial_num', [], ...
                'is_replay_trial', false, ...
                'replay_instance', 0, ...
                'iti_duration_sec', NaN, ...
                'iti_source_trial_index', NaN, ...
                'fixation_break', false, ...
                'fixation_break_frame', NaN, ...
                'fixation_break_time', NaN, ...
                'replay_queued', false, ...
                'source_index', [], ...
                'condition_label', '', ...
                'condition_code', [], ...
                'condition_trigger', [], ...
                'sequence_trigger_code', [], ...
                'sequence_trigger_frame', [], ...
                'dot_color_rgb', [], ...
                'occlusion_mode', '', ...
                'occlusion_enabled', [], ...
                'deviance_frame', [], ...
                'occlusion_start_frame', [], ...
                'occlusion_complete_frame', [], ...
                'occlusion_end_frame', [], ...
                'occlusion_end_complete_frame', [], ...
                'pathband_post_deactivate_frame', [], ...
                'catch_type_code', [], ...
                'catch_type_label', '', ...
                'catch_expected_response_code', [], ...
                'catch_branch_changed_path', [], ...
                'catch_disappear_frame', [], ...
                'catch_reappear_frame', [], ...
                'catch_alt_source_index', [], ...
                'catch_question_presented', false, ...
                'catch_question_start_time', [], ...
                'catch_question_end_time', [], ...
                'catch_response_code', [], ...
                'catch_response_label', '', ...
                'catch_response_rt_sec', [], ...
                'catch_response_correct', [], ...
                'catch_response_device', '', ...
                'catch_timed_out', false, ...
                'frame_flip_time', [], ...
                'frame_vbl_timestamp', [], ...
                'frame_stimulus_onset_time', [], ...
                'frame_flip_timestamp', [], ...
                'frame_missed', [], ...
                'frame_missed_count', 0, ...
                'frame_missed_frames', [], ...
                'frame_missed_values', [], ...
                'frame_duration_ms', [], ...
                'trial_start_time', [], ...
                'trial_end_time', []);
            output = repmat(outputTemplate, 1, numel(trialSchedule));

            abortRequested = false;
            abortReason = '';
            completedTrials = 0;
            completedCatches = 0;

            iTrial = 1;
            while iTrial <= numel(trialSchedule)
                if iTrial > numel(output)
                    output(iTrial) = outputTemplate;
                end
                srcIdx = trialSchedule(iTrial);
                isReplayTrial = trialIsReplay(iTrial);
                replayInstance = trialReplayInstance(iTrial);
                itiSourceTrialIndex = local_pick_from_vector( ...
                    itiSourceTrialIndexSchedule, iTrial, NaN);
                itiDurationSec = local_pick_from_vector( ...
                    itiPrecomputedSec, itiSourceTrialIndex, NaN);
                if ~isfinite(itiDurationSec) || itiDurationSec < 0
                    error(['Invalid ITI duration lookup for block %d run %d executed trial %d ' ...
                        '(source trial index %g).'], ...
                        currentBlock, currentRun, iTrial, itiSourceTrialIndex);
                end
                trial = xyAll(srcIdx);
                trialStimBase = squeeze(stimuliPx(srcIdx, :, :));

                condLabel = char(string(trial.condition_label));
                condTrigger = local_condition_trigger(condLabel, Conf.trigger);
                sequenceTriggerCode = sequenceTriggerBySource(srcIdx);
                sequenceTriggerFrame = round(double(Conf.sequenceIdentityTriggerFrame));

                if isReplayTrial
                    catchTypeCode = 0;
                    catchExpectedCode = 0;
                    catchBranchChanged = NaN;
                    catchDisappearFrame = NaN;
                    catchReappearFrame = NaN;
                    catchAltSourceIdx = NaN;
                else
                    baseTrialIdx = min(iTrial, numel(catchTypeVec));
                    catchTypeCode = local_pick_from_vector(catchTypeVec, baseTrialIdx, 0);
                    catchExpectedCode = local_pick_from_vector(catchExpectedVec, baseTrialIdx, 0);
                    catchBranchChanged = local_pick_from_vector(catchBranchVec, baseTrialIdx, NaN);
                    catchDisappearFrame = local_pick_from_vector(catchDisappearVec, baseTrialIdx, NaN);
                    catchReappearFrame = local_pick_from_vector(catchReappearVec, baseTrialIdx, NaN);
                    catchAltSourceIdx = local_pick_from_vector(catchAltSourceVec, baseTrialIdx, NaN);
                end

                isCatchType1 = (catchTypeCode == 1);
                isCatchType2 = (catchTypeCode == 2);
                if ~(isCatchType1 || isCatchType2)
                    catchTypeCode = 0;
                end

                % Read event frames from metadata and clamp to valid range.
                devFrame = local_clamp_frame(double(trial.deviance_frame), framesPerTrial);
                occStartFrame = local_clamp_frame(double(trial.occlusion_start_frame), framesPerTrial);
                occCompleteFrame = local_clamp_frame(double(trial.occlusion_complete_frame), framesPerTrial);
                occEndStartFrame = local_clamp_frame(double(trial.occlusion_end_frame), framesPerTrial);
                occEndCompleteFrame = local_clamp_frame(double(trial.occlusion_end_complete_frame), framesPerTrial);

                postDeactivateFrame = framesPerTrial;
                if isfield(trial, 'pathband_post_deactivate_frame') && ~isempty(trial.pathband_post_deactivate_frame)
                    postDeactivateFrame = local_clamp_frame(double(trial.pathband_post_deactivate_frame), framesPerTrial);
                end

                pathbandLineWidthPx = max(1, round(double(trial.pathband_width_deg) * ppd));
                pathbandPreCoords = local_polyline_to_px(double(trial.pathband_pre_xy), ppd, rectCoords);
                pathbandPostCoords = local_polyline_to_px(double(trial.pathband_post_xy), ppd, rectCoords);
                % Pre/post bands use separate draw options so straight-mode
                % entrance widening can be applied only where needed.
                pathbandPreDrawOpts = local_make_pathband_draw_options( ...
                    Conf, true, pathbandLineWidthPx, dotRadiusPx, false);
                pathbandPostDrawOpts = local_make_pathband_draw_options( ...
                    Conf, true, pathbandLineWidthPx, dotRadiusPx, true);
                pathOverlayDrawOpts = local_make_pathband_draw_options( ...
                    Conf, false, pathbandLineWidthPx, dotRadiusPx, false);
                % Build draw plans once per trial so style transforms and
                % straight-band polygons are not rebuilt every frame.
                pathbandPrePlan = local_precompute_polyline_draw_plan( ...
                    pathbandPreCoords, pathbandLineWidthPx, pathbandPreDrawOpts);
                pathbandPostPlan = local_precompute_polyline_draw_plan( ...
                    pathbandPostCoords, pathbandLineWidthPx, pathbandPostDrawOpts);

                trialStim = trialStimBase;
                dotVisibleMask = true(framesPerTrial, 1);
                if isCatchType1
                    altStim = trialStimBase;
                    if isfinite(catchAltSourceIdx) && catchAltSourceIdx >= 1 && catchAltSourceIdx <= numel(xyAll)
                        altStim = squeeze(stimuliPx(round(catchAltSourceIdx), :, :));
                    end
                    [trialStim, dotVisibleMask, catchDisappearFrame, catchReappearFrame] = ...
                        local_prepare_type1_catch_stimulus( ...
                        trialStimBase, altStim, catchBranchChanged, catchDisappearFrame, catchReappearFrame);
                end

                % Optional debug path overlay uses the same trajectory that
                % is fed to frame-by-frame drawing for this trial.
                pathOverlayCoords = [];
                pathOverlayPlan = struct([]);
                if logical(Conf.debugShowFullPathOverlay)
                    finitePathMask = isfinite(trialStim(:, 1)) & isfinite(trialStim(:, 2));
                    if nnz(finitePathMask) >= 2
                        pathOverlayCoords = local_polyline_px_to_drawlines_coords( ...
                            trialStim(finitePathMask, 1:2));
                        pathOverlayPlan = local_precompute_polyline_draw_plan( ...
                            pathOverlayCoords, Conf.debugPathLineWidthPx, pathOverlayDrawOpts);
                    end
                end

                % Store static trial metadata.
                output(iTrial).trial_num = iTrial;
                output(iTrial).is_replay_trial = logical(isReplayTrial);
                output(iTrial).replay_instance = replayInstance;
                output(iTrial).iti_duration_sec = itiDurationSec;
                output(iTrial).iti_source_trial_index = itiSourceTrialIndex;
                output(iTrial).fixation_break = false;
                output(iTrial).fixation_break_frame = NaN;
                output(iTrial).fixation_break_time = NaN;
                output(iTrial).replay_queued = false;
                output(iTrial).source_index = srcIdx;
                output(iTrial).condition_label = condLabel;
                output(iTrial).condition_code = trial.condition;
                output(iTrial).condition_trigger = condTrigger;
                output(iTrial).sequence_trigger_code = sequenceTriggerCode;
                output(iTrial).sequence_trigger_frame = sequenceTriggerFrame;
                output(iTrial).dot_color_rgb = runDotColor;
                output(iTrial).occlusion_mode = Conf.occlusionMode;
                output(iTrial).occlusion_enabled = logical(trial.occlusion_enabled);
                output(iTrial).deviance_frame = devFrame;
                output(iTrial).occlusion_start_frame = occStartFrame;
                output(iTrial).occlusion_complete_frame = occCompleteFrame;
                output(iTrial).occlusion_end_frame = occEndStartFrame;
                output(iTrial).occlusion_end_complete_frame = occEndCompleteFrame;
                output(iTrial).pathband_post_deactivate_frame = postDeactivateFrame;
                output(iTrial).catch_type_code = catchTypeCode;
                output(iTrial).catch_type_label = local_catch_type_label(catchTypeCode);
                output(iTrial).catch_expected_response_code = catchExpectedCode;
                output(iTrial).catch_branch_changed_path = catchBranchChanged;
                output(iTrial).catch_disappear_frame = catchDisappearFrame;
                output(iTrial).catch_reappear_frame = catchReappearFrame;
                output(iTrial).catch_alt_source_index = catchAltSourceIdx;

                % No pre-trial pause in v18 rescueTraject mode:
                % timing between trials is handled by an explicit post-trial interval.

                % Start EyeLink recording for this trial (real hardware mode).
                % Data flow:
                %   trial metadata + replay status -> EyeLink TRIALID/status
                %   messages and StartRecording.
                if Conf.useEyelink
                    local_eyelink_start_trial_recording(Conf, iSub, currentBlock, currentRun, ...
                        iTrial, srcIdx, isReplayTrial, replayInstance, screenNumber);
                end

                % Replay-start trigger (151): emitted once per run immediately
                % before the first replayed trial.
                if isReplayTrial && ~replayStartTriggerSent
                    [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                        local_emit_trigger(Conf, Conf.trigger.replayStart, iTrial, 0, ...
                        'replay_start', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                    replayStartTriggerSent = true;
                    WaitSecs(0.1);
                end

                % Trial start trigger (condition-specific).
                [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                    local_emit_trigger(Conf, condTrigger, iTrial, 1, ...
                    'trial_condition_onset', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);

                % Frame loop with fixed-timestep flipping.
                frameVblTimestamp = nan(framesPerTrial, 1);
                frameStimulusOnsetTime = nan(framesPerTrial, 1);
                frameFlipTimestamp = nan(framesPerTrial, 1);
                frameMissed = nan(framesPerTrial, 1);
                frameMissedFrames = zeros(framesPerTrial, 1);
                frameMissedValues = nan(framesPerTrial, 1);
                frameMissedCount = 0;
                gazeOutsideCounter = 0;
                gazeBreakTriggered = false;
                gazeBreakFrame = NaN;
                gazeBreakTime = NaN;
                [fakeBreakThisTrial, fakeBreakStartFrame, fakeBreakEndFrame] = ...
                    local_plan_fake_gaze_break(Conf, framesPerTrial);
                [vbl, trialStimulusOnsetTime, trialFlipTimestamp, trialMissed] = Screen('Flip', window);
                output(iTrial).trial_start_time = trialStimulusOnsetTime;
                trialFrameTimer = tic;

                for iFrame = 1:framesPerTrial
                    local_draw_arena(window, rectCoords, Conf);

                    xDot = trialStim(iFrame, 1);
                    yDot = trialStim(iFrame, 2);
                    dotVisibleThisFrame = dotVisibleMask(iFrame) && isfinite(xDot) && isfinite(yDot);

                    % Fixation-break monitoring (real EyeLink samples or fake schedule).
                    % If fixation is lost for enough consecutive frames, trigger 150,
                    % mark trial as broken, and queue one replay (configurable).
                    if Conf.useGazeMonitoring && ~gazeBreakTriggered
                        hasSample = false;
                        eyeWithinWindow = false;
                        if Conf.useFakeEyeTracker
                            hasSample = true;
                            if fakeBreakThisTrial && iFrame >= fakeBreakStartFrame && iFrame <= fakeBreakEndFrame
                                eyeWithinWindow = false;
                            else
                                eyeWithinWindow = true;
                            end
                        elseif Conf.useEyelink
                            [hasSample, eyeWithinWindow] = local_sample_eyelink_fixation(Conf, xCenter, yCenter);
                        end
                        if hasSample
                            if eyeWithinWindow
                                gazeOutsideCounter = 0;
                            else
                                gazeOutsideCounter = gazeOutsideCounter + 1;
                            end
                            if gazeOutsideCounter >= Conf.fixBreakToleranceFrames
                                gazeBreakTriggered = true;
                                gazeBreakFrame = iFrame;
                                gazeBreakTime = GetSecs;
                                output(iTrial).fixation_break = true;
                                output(iTrial).fixation_break_frame = gazeBreakFrame;
                                output(iTrial).fixation_break_time = gazeBreakTime - trialStimulusOnsetTime;
                                if Conf.useEyelink
                                    local_eyelink_safe_message(sprintf('GAZE_BREAK frame %d', gazeBreakFrame));
                                end
                                [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                    local_emit_trigger(Conf, Conf.trigger.gazeBreak, iTrial, iFrame, ...
                                    'gaze_break', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                                if replayCountsBySource(srcIdx) < Conf.maxReplaysPerTrial
                                    replayCountsBySource(srcIdx) = replayCountsBySource(srcIdx) + 1;
                                    trialSchedule(end + 1) = srcIdx; %#ok<AGROW>
                                    trialIsReplay(end + 1) = true; %#ok<AGROW>
                                    trialReplayInstance(end + 1) = replayCountsBySource(srcIdx); %#ok<AGROW>
                                    itiSourceTrialIndexSchedule(end + 1) = itiSourceTrialIndex; %#ok<AGROW>
                                    output(iTrial).replay_queued = true;
                                    replayQueuedCount = replayQueuedCount + 1;
                                end
                                break;
                            end
                        end
                    end

                    % Draw dot under selected occlusion mode.
                    if strcmpi(Conf.occlusionMode, 'alpha')
                        if dotVisibleThisFrame
                            alpha = 1.0;
                            if isfield(trial, 'visibility_alpha') && ~isempty(trial.visibility_alpha)
                                alpha = double(trial.visibility_alpha(iFrame));
                            end
                            alpha = max(0, min(1, alpha));
                            colorAlpha = local_blend_with_background(runDotColor, Conf.rectColor, alpha);
                            Screen('DrawDots', window, [xDot; yDot], dotSizePx, colorAlpha, [0 0], 1);
                        end
                    else
                        % Path-band mode: draw dot, then mask using the pre/post
                        % trajectory bands selected by frame index.
                        if dotVisibleThisFrame
                            Screen('DrawDots', window, [xDot; yDot], dotSizePx, runDotColor, [0 0], 1);
                        end

                        if logical(trial.occlusion_enabled) && ~isCatchType1
                            preActive = (iFrame < devFrame);
                            postActive = (iFrame >= devFrame) && (iFrame <= postDeactivateFrame);

                            % Pre-band is explicitly disabled at deviance transition.
                            if preActive
                                local_draw_polyline_plan(window, pathbandPrePlan, Conf.rectColor);
                            end
                            if postActive
                                local_draw_polyline_plan(window, pathbandPostPlan, Conf.rectColor);
                            end
                        end
                    end

                    % Optional persistent path overlay for debugging.
                    % Draw it above dot/occluder so the full trajectory
                    % remains visible even when the occluder overlaps it.
                    if logical(Conf.debugShowFullPathOverlay) && ~isempty(pathOverlayCoords)
                        local_draw_polyline_plan(window, pathOverlayPlan, Conf.debugPathLineColor);
                    end

                    % Draw fixation after dot/mask so it is never occluded.
                    local_draw_fixation(window, xCenter, yCenter, fixSizePx, Conf.fixColor);

                    % Optional frame/timing overlay in the upper-left corner.
                    if logical(Conf.debugShowFrameInfoOverlay)
                        elapsedSec = toc(trialFrameTimer);
                        overlayText = sprintf('trial %d/%d, frame %d/%d | t=%.3fs', ...
                            iTrial, numel(trialSchedule), iFrame, framesPerTrial, elapsedSec);
                        overlayX = windowRect(1) + round(Conf.debugOverlayTextOffsetPx(1));
                        overlayY = windowRect(2) + round(Conf.debugOverlayTextOffsetPx(2));
                        DrawFormattedText(window, overlayText, overlayX, overlayY, Conf.debugOverlayTextColor);
                    end

                    % Sequence-identity trigger: emitted once per trial at fixed frame offset.
                    if iFrame == sequenceTriggerFrame
                        [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                            local_emit_trigger(Conf, sequenceTriggerCode, iTrial, iFrame, ...
                            'trial_sequence_identity', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                    end

                    % Emit occlusion event triggers exactly at configured frames.
                    if isCatchType1
                        if isfinite(catchDisappearFrame) && iFrame == round(catchDisappearFrame)
                            [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                local_emit_trigger(Conf, Conf.trigger.catch1Disappear, iTrial, iFrame, ...
                                'catch1_disappear', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                        end
                        if isfinite(catchReappearFrame) && iFrame == round(catchReappearFrame)
                            [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                local_emit_trigger(Conf, Conf.trigger.catch1Reappear, iTrial, iFrame, ...
                                'catch1_reappear', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                        end
                    elseif logical(trial.occlusion_enabled)
                        if iFrame == occStartFrame
                            [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                local_emit_trigger(Conf, Conf.trigger.occlusionStart, iTrial, iFrame, ...
                                'occlusion_start', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                        end
                        if iFrame == occCompleteFrame
                            [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                local_emit_trigger(Conf, Conf.trigger.occlusionComplete, iTrial, iFrame, ...
                                'occlusion_complete', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                        end
                        if iFrame == occEndStartFrame
                            [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                local_emit_trigger(Conf, Conf.trigger.occlusionEndStart, iTrial, iFrame, ...
                                'occlusion_end_start', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                        end
                        if iFrame == occEndCompleteFrame
                            [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
                                local_emit_trigger(Conf, Conf.trigger.occlusionEndComplete, iTrial, iFrame, ...
                                'occlusion_end_complete', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                        end
                    end

                    if iFrame == 1
                        [flipVbl, flipStimulusOnset, flipTimestamp, flipMissed] = Screen('Flip', window);
                    else
                        requestedWhen = vbl + (1 - 0.5) * ifi;
                        [flipVbl, flipStimulusOnset, flipTimestamp, flipMissed] = ...
                            Screen('Flip', window, requestedWhen);
                    end
                    frameVblTimestamp(iFrame) = flipVbl;
                    frameStimulusOnsetTime(iFrame) = flipStimulusOnset;
                    frameFlipTimestamp(iFrame) = flipTimestamp;
                    frameMissed(iFrame) = flipMissed;
                    vbl = flipVbl;

                    if isfinite(flipMissed) && (flipMissed > Conf.flipMissedEpsilonSec)
                        frameMissedCount = frameMissedCount + 1;
                        frameMissedFrames(frameMissedCount) = iFrame;
                        frameMissedValues(frameMissedCount) = flipMissed;
                        if frameMissedCount == Conf.flipMissedWarnThresholdCount
                            warning(['Run %d trial %d reached missed-flip warning threshold (%d). ' ...
                                'Current frame=%d, missed=%.6f s.'], ...
                                currentRun, iTrial, Conf.flipMissedWarnThresholdCount, iFrame, flipMissed);
                        end
                        if Conf.flipMissedAbortEnabled && ...
                                frameMissedCount >= Conf.flipMissedAbortThresholdCount
                            abortRequested = true;
                            abortReason = 'flip_missed_abort_threshold';
                            fprintf(['Run %d aborted on trial %d: missed-flip count %d ' ...
                                'reached threshold %d.\n'], ...
                                currentRun, iTrial, frameMissedCount, Conf.flipMissedAbortThresholdCount);
                            break;
                        end
                    end

                    [keyDown, ~, keyCode] = KbCheck;
                    if keyDown && keyCode(Conf.escapeKeyCode)
                        abortRequested = true;
                        abortReason = 'escape_key';
                        KbReleaseWait;
                        break;
                    end
                end

                if Conf.useEyelink
                    local_eyelink_stop_trial_recording(Conf);
                end

                if gazeBreakTriggered
                    local_present_fixation_break_warning(window, Conf);
                end

                if abortRequested
                    if ~strcmp(abortReason, 'escape_key')
                        fprintf('Run %d interrupted (%s) during trial %d.\n', currentRun, abortReason, iTrial);
                    end
                    break;
                end

                output(iTrial).frame_flip_time = frameVblTimestamp;
                output(iTrial).frame_vbl_timestamp = frameVblTimestamp;
                output(iTrial).frame_stimulus_onset_time = frameStimulusOnsetTime;
                output(iTrial).frame_flip_timestamp = frameFlipTimestamp;
                output(iTrial).frame_missed = frameMissed;
                output(iTrial).frame_missed_count = frameMissedCount;
                output(iTrial).frame_missed_frames = frameMissedFrames(1:frameMissedCount);
                output(iTrial).frame_missed_values = frameMissedValues(1:frameMissedCount);
                validFlipMask = isfinite(frameVblTimestamp);
                if any(validFlipMask)
                    lastFlipIdx = find(validFlipMask, 1, 'last');
                    output(iTrial).trial_end_time = frameVblTimestamp(lastFlipIdx);
                else
                    output(iTrial).trial_end_time = NaN;
                end

                if framesPerTrial > 1
                    validFlip = frameVblTimestamp(isfinite(frameVblTimestamp));
                    if numel(validFlip) > 1
                        output(iTrial).frame_duration_ms = diff(validFlip) * 1000;
                    else
                        output(iTrial).frame_duration_ms = [];
                    end
                else
                    output(iTrial).frame_duration_ms = [];
                end

                if ~gazeBreakTriggered && Conf.enableCatchTrials && (catchTypeCode > 0)
                    [output(iTrial), debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel, abortFromQuestion] = ...
                        local_run_catch_question( ...
                        output(iTrial), window, Conf, iTrial, xCenter, yCenter, fixSizePx, ...
                        debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);
                    if abortFromQuestion
                        abortRequested = true;
                        abortReason = 'escape_key';
                        completedTrials = iTrial;
                        break;
                    end
                    completedCatches = completedCatches + 1;
                end

                completedTrials = iTrial;
                local_present_post_trial_interval(window, Conf, xCenter, yCenter, fixSizePx, itiDurationSec);
                iTrial = iTrial + 1;
            end

            output = output(1:completedTrials);
            executedTrialOrder = trialSchedule(1:completedTrials);
            itiSourceTrialIndexUsed = itiSourceTrialIndexSchedule(1:completedTrials);
            if isempty(output)
                itiUsedSec = [];
            else
                itiUsedSec = [output.iti_duration_sec];
            end

            runResults(iRunIdx).run_num = currentRun;
            runResults(iRunIdx).trialOrder = executedTrialOrder;
            runResults(iRunIdx).trialOrderBase = trialOrder;
            runResults(iRunIdx).trialIsReplay = trialIsReplay(1:completedTrials);
            runResults(iRunIdx).trialReplayInstance = trialReplayInstance(1:completedTrials);
            runResults(iRunIdx).itiRangeSec = Conf.ITIRangeSec;
            runResults(iRunIdx).itiPrecomputeSeed = itiSeed;
            runResults(iRunIdx).itiPrecomputedSec = itiPrecomputedSec;
            runResults(iRunIdx).itiUsedSec = itiUsedSec;
            runResults(iRunIdx).itiSourceTrialIndex = itiSourceTrialIndexUsed;
            runResults(iRunIdx).output = output;
            runResults(iRunIdx).nTrialsPlanned = nTrials;
            runResults(iRunIdx).nTrialsCompleted = completedTrials;
            runResults(iRunIdx).nCatchPlanned = sum(catchTypeVec > 0);
            runResults(iRunIdx).nCatchCompleted = completedCatches;
            runResults(iRunIdx).nReplaysQueued = replayQueuedCount;
            runResults(iRunIdx).aborted = abortRequested;
            runResults(iRunIdx).abort_reason = abortReason;
            runResults(iRunIdx).debugTriggerLog = debugTriggerLog;
            runResults(iRunIdx).debugTriggerTrial = debugTriggerTrial;
            runResults(iRunIdx).debugTriggerFrame = debugTriggerFrame;
            runResults(iRunIdx).debugTriggerSec = debugTriggerSec;
            runResults(iRunIdx).debugTriggerLabel = debugTriggerLabel;
            runResults(iRunIdx).eyeTrackerEdfFile = Conf.eyeTrackerEdfFile;
            runResults(iRunIdx).eyeTrackerLocalEdfPath = Conf.eyeTrackerLocalEdfPath;

            [scoredCount, answeredCount, correctCount] = local_count_run_catch_accuracy(output);
            runCatchScored(currentRun) = runCatchScored(currentRun) + scoredCount;
            runCatchAnswered(currentRun) = runCatchAnswered(currentRun) + answeredCount;
            runCatchCorrect(currentRun) = runCatchCorrect(currentRun) + correctCount;

            % Report measured frame pacing to help diagnose display-dependent lag.
            local_report_run_frame_timing(output, currentBlock, currentRun, actualFrameRate);

            if Conf.debug || ~isempty(debugTriggerLog)
                debugTriggerFile = fullfile(outputDir, sprintf( ...
                    'debug_actual_triggers_occlusion_v18_rescueTraject_sub%02d_block%02d_run%02d.csv', ...
                    iSub, currentBlock, currentRun));
                debugTriggerItiSec = nan(size(debugTriggerTrial(:)));
                if ~isempty(output)
                    itiByTrial = [output.iti_duration_sec];
                    validItiMask = isfinite(debugTriggerTrial(:)) & ...
                        debugTriggerTrial(:) >= 1 & debugTriggerTrial(:) <= numel(itiByTrial);
                    mappedTrialIdx = round(double(debugTriggerTrial(validItiMask)));
                    debugTriggerItiSec(validItiMask) = itiByTrial(mappedTrialIdx);
                end
                T = table(debugTriggerLog(:), debugTriggerTrial(:), debugTriggerFrame(:), ...
                    debugTriggerSec(:), string(debugTriggerLabel(:)), debugTriggerItiSec(:), ...
                    'VariableNames', {'trigger', 'trial', 'frame', 'seconds', 'label', 'iti_duration_sec'});
                writetable(T, debugTriggerFile);
                runResults(iRunIdx).debugTriggerFile = debugTriggerFile;
                fprintf('Debug trigger log saved: %s\n', debugTriggerFile);
            end

            if abortRequested
                blockAborted = true;
                blockAbortReason = abortReason;
                break;
            end
        end

        if Conf.useEyelink
            Conf = local_eyelink_finalize_block_file(Conf);
        end

        save(outputMat, 'runResults', 'Conf', 'iSub', 'currentBlock', 'viewingDistanceMm', ...
            'inputFile', 'trialStructFile', 'runsToExecute', 'debug_runNumber', 'catchPlan');
        fprintf('Saved block output: %s\n', outputMat);

        if blockAborted
            if strcmp(blockAbortReason, 'escape_key')
                local_emit_escape_termination_trigger(Conf, currentBlock, 'escape_during_block');
                fprintf('esc pressed during block %d\n', currentBlock);
            end
            stopExperimentAfterBlock = true;
            stopReason = blockAbortReason;
            break;
        end

        isLastBlock = (currentBlock == nBlocksAvailable);
        if ~isLastBlock
            remainingBlocks = nBlocksAvailable - currentBlock;
            breakAction = local_wait_end_of_block_action(window, Conf, currentBlock, remainingBlocks);
            if strcmp(breakAction, 'abort')
                stopExperimentAfterBlock = true;
                stopReason = 'escape_after_block_message';
                local_emit_escape_termination_trigger(Conf, currentBlock, 'escape_after_block');
                fprintf('esc pressed at the end of block %d\n', currentBlock);
                break;
            end

            if logical(Conf.betweenBlockCalibrationChoiceEnabled)
                calibrationAction = local_wait_between_block_calibration_action(window, Conf);
                if strcmp(calibrationAction, 'abort')
                    stopExperimentAfterBlock = true;
                    stopReason = 'escape_after_block_message';
                    local_emit_escape_termination_trigger(Conf, currentBlock, 'escape_after_block');
                    fprintf('esc pressed at the end of block %d\n', currentBlock);
                    break;
                elseif strcmp(calibrationAction, 'calibrate')
                    Conf = local_eyelink_run_calibration(Conf);
                end
            end

            if logical(Conf.endOfBlockMessageEnabled)
                abortedDuringFixation = local_show_post_message_fixation(window, Conf, xCenter, yCenter, fixSizePx);
                if abortedDuringFixation
                    stopExperimentAfterBlock = true;
                    stopReason = 'escape_during_post_message_fixation';
                    local_emit_escape_termination_trigger(Conf, currentBlock, 'escape_after_block');
                    fprintf('esc pressed at the end of block %d\n', currentBlock);
                    break;
                end
            end
        else
            finalAction = local_show_final_message(window, Conf);
            if strcmp(finalAction, 'abort')
                stopExperimentAfterBlock = true;
                stopReason = 'escape_on_final_message';
                local_emit_escape_termination_trigger(Conf, currentBlock, 'escape_termination');
                break;
            end

            if logical(Conf.finalMessageEnabled)
                abortedDuringFixation = local_show_post_message_fixation(window, Conf, xCenter, yCenter, fixSizePx);
                if abortedDuringFixation
                    stopExperimentAfterBlock = true;
                    stopReason = 'escape_during_post_message_fixation';
                    local_emit_escape_termination_trigger(Conf, currentBlock, 'escape_termination');
                end
                break;
            end
        end
    end

    if stopExperimentAfterBlock
        fprintf('Experiment ended early (%s).\n', stopReason);
    end

    local_print_run_accuracy_summary(runCatchCorrect, runCatchScored, runCatchAnswered);

    local_cleanup_runtime(Conf);
catch ME
    local_cleanup_runtime(Conf);
    rethrow(ME);
end

%% Local helper functions
function triggerValue = local_condition_trigger(conditionLabel, triggerMap)
% LOCAL_CONDITION_TRIGGER Map condition label to condition-specific onset trigger.
switch lower(strtrim(conditionLabel))
    case 'always_visible'
        triggerValue = triggerMap.conditionAlwaysVisible;
    case 'occluded_nondeviant'
        triggerValue = triggerMap.conditionOccludedNondeviant;
    case 'occluded_deviant'
        triggerValue = triggerMap.conditionOccludedDeviant;
    otherwise
        error('Unsupported condition label: %s', conditionLabel);
end
end

function skipSyncTests = local_resolve_skip_sync_tests(Conf)
% LOCAL_RESOLVE_SKIP_SYNC_TESTS Resolve PTB sync-test policy for this run.
%
% Policy:
%   1) explicit workspace override (debugSkipSyncTests) takes priority.
%   2) otherwise, MEG uses schedule default skipSyncTestsWhenMEG.
%   3) non-MEG can enable SkipSyncTests only in explicit debug modes.

if isfield(Conf, 'debugSkipSyncTestsOverride') && ~isempty(Conf.debugSkipSyncTestsOverride)
    skipSyncTests = logical(Conf.debugSkipSyncTestsOverride);
else
    if isfield(Conf, 'MEG') && logical(Conf.MEG)
        if isfield(Conf, 'skipSyncTestsWhenMEG')
            skipSyncTests = logical(Conf.skipSyncTestsWhenMEG);
        else
            skipSyncTests = false;
        end
    else
        debugLikeMode = logical(Conf.debugTriggerMode) || ...
            logical(Conf.debugShowFrameInfoOverlay) || logical(Conf.debugShowFullPathOverlay);
        skipSyncTests = logical(Conf.skipSyncTestsInDebug) && debugLikeMode;
    end
end

% Safety: do not allow accidental SkipSyncTests in MEG unless explicitly enabled.
if isfield(Conf, 'MEG') && logical(Conf.MEG) && ...
        isfield(Conf, 'skipSyncTestsWhenMEG') && ~logical(Conf.skipSyncTestsWhenMEG)
    skipSyncTests = false;
end
end

function value = local_schedule_get_or_default(scheduleCfg, fieldName, defaultValue)
% LOCAL_SCHEDULE_GET_OR_DEFAULT Read schedule config field with fallback.
if isprop(scheduleCfg, fieldName)
    value = scheduleCfg.(fieldName);
else
    value = defaultValue;
end
end

function maxReplays = local_compute_max_replays_per_trial(replayLimit, allowInfiniteReplays)
% LOCAL_COMPUTE_MAX_REPLAYS_PER_TRIAL Resolve replay cap used by gaze breaks.
if logical(allowInfiniteReplays)
    maxReplays = intmax('int32');
else
    replayLimit = max(0, round(double(replayLimit)));
    maxReplays = replayLimit;
end
end

function Conf = local_setup_eyelink_session(Conf, window, windowRect, screenNumber, iSub, iBlock, outputDir)
% LOCAL_SETUP_EYELINK_SESSION Initialize EyeLink setup/calibration session.
%
% Data flow:
%   runtime flags + PTB window geometry -> configured EyeLink session stored
%   in Conf fields for trial-level recording, optional calibration, and cleanup.
if ~isfield(Conf, 'useEyelinkSession') || ~Conf.useEyelinkSession
    return;
end
if exist('Eyelink', 'file') ~= 3 && exist('Eyelink', 'file') ~= 2
    error('EyeLink session requested, but Eyelink MEX is not available on MATLAB path.');
end
if EyelinkInit() ~= 1
    error('EyeLink initialization failed (EyelinkInit ~= 1).');
end
el = EyelinkInitDefaults(window);
el.foregroundcolour = 0;
el.backgroundcolour = double(Conf.background(:)');
Conf.eyeTrackerEyelinkDefaults = el;
Conf.eyeTrackerBlockFileOpen = false;
Conf.eyeTrackerEdfFile = '';
Conf.eyeTrackerLocalEdfPath = '';
Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
screenXpixels = windowRect(3) - windowRect(1);
screenYpixels = windowRect(4) - windowRect(2);
Eyelink('command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
Eyelink('command', 'calibration_type = HV9');
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');
Eyelink('command', 'saccade_velocity_threshold = 35');
Eyelink('command', 'saccade_acceleration_threshold = 9500');
if Eyelink('IsConnected') ~= 1
    error('EyeLink disconnected during setup.');
end
if Conf.eyeTrackerDoCalibration
    EyelinkDoTrackerSetup(Conf.eyeTrackerEyelinkDefaults);
end
if Conf.eyeTrackerImageTransferEnabled
    imgFile = Conf.eyeTrackerFixBmpPath;
    if exist(imgFile, 'file') == 2
        try
            imgInfo = imfinfo(imgFile);
            [screenW, screenH] = Screen('WindowSize', screenNumber);
            trackerX = round(screenW / 2 - imgInfo.Width / 2);
            trackerY = round(screenH / 2 - imgInfo.Height / 2);
            Eyelink('ImageTransfer', imgFile, 0, 0, imgInfo.Width, imgInfo.Height, trackerX, trackerY, 1);
        catch
            % Keep setup resilient when host image transfer is unsupported.
        end
    end
end
WaitSecs(max(0, double(Conf.eyeTrackerWaitSec)));
Conf.eyeTrackerInitialized = true;
end

function Conf = local_eyelink_open_block_file(Conf, iSub, currentBlock, outputDir)
% LOCAL_EYELINK_OPEN_BLOCK_FILE Open one EDF file for the current block.
if ~Conf.useEyelink || ~Conf.eyeTrackerInitialized || ...
        ~isfield(Conf, 'eyeTrackerBlockFileOpen') || ~Conf.eyeTrackerBlockFileOpen
    return;
end
if isfield(Conf, 'eyeTrackerBlockFileOpen') && Conf.eyeTrackerBlockFileOpen
    Conf = local_eyelink_finalize_block_file(Conf);
end

selectedBase = '';
selectedLocalPath = '';
for attemptIdx = 0:99
    candidateBase = local_make_eyelink_block_edf_basename(iSub, currentBlock, attemptIdx);
    candidateLocalPath = fullfile(outputDir, [candidateBase '.edf']);
    if exist(candidateLocalPath, 'file') ~= 2
        selectedBase = candidateBase;
        selectedLocalPath = candidateLocalPath;
        break;
    end
end
if isempty(selectedBase)
    error('Could not allocate unique EyeLink EDF basename for subject %d block %d.', iSub, currentBlock);
end

Conf.eyeTrackerEdfFile = [selectedBase '.edf'];
Conf.eyeTrackerLocalEdfPath = selectedLocalPath;
Eyelink('Openfile', Conf.eyeTrackerEdfFile);
Eyelink('command', sprintf(['add_file_preamble_text ' ...
    '''MoveDot occlusion v18 rescueTraject: subject %d ; block %d ; practice %d ; time %s'''], ...
    iSub, currentBlock, double(Conf.practiceModeEnabled), datestr(now, 'YYYYmmddHHMMSS')));
Conf.eyeTrackerBlockFileOpen = true;
end

function Conf = local_eyelink_finalize_block_file(Conf)
% LOCAL_EYELINK_FINALIZE_BLOCK_FILE Stop recording, receive EDF, and close file for this block.
if ~Conf.useEyelink || ~Conf.eyeTrackerInitialized || ...
        ~isfield(Conf, 'eyeTrackerBlockFileOpen') || ~Conf.eyeTrackerBlockFileOpen
    return;
end
if ~isfield(Conf, 'eyeTrackerBlockFileOpen') || ~Conf.eyeTrackerBlockFileOpen
    return;
end

try
    Eyelink('StopRecording');
catch
end
try
    Eyelink('Command', 'set_idle_mode');
catch
end
if isfield(Conf, 'eyeTrackerReceiveFileAfterBlock') && Conf.eyeTrackerReceiveFileAfterBlock && ...
        isfield(Conf, 'practiceModeEnabled') && ~Conf.practiceModeEnabled
    try
        if isfield(Conf, 'eyeTrackerEdfFile') && ~isempty(Conf.eyeTrackerEdfFile)
            fprintf('Receiving EyeLink EDF: %s\n', Conf.eyeTrackerEdfFile);
            status = Eyelink('ReceiveFile');
            if status > 0
                fprintf('EyeLink ReceiveFile status: %d\n', status);
            end
            if isfield(Conf, 'eyeTrackerLocalEdfPath') && ~isempty(Conf.eyeTrackerLocalEdfPath) && ...
                    exist(Conf.eyeTrackerEdfFile, 'file') == 2
                try
                    movefile(Conf.eyeTrackerEdfFile, Conf.eyeTrackerLocalEdfPath, 'f');
                catch
                    % Keep cleanup resilient if file move fails.
                end
            end
        end
    catch
        fprintf('Warning: failed to receive EyeLink EDF file.\n');
    end
end
try
    Eyelink('CloseFile');
catch
end
Conf.eyeTrackerBlockFileOpen = false;
end

function baseName = local_make_eyelink_block_edf_basename(iSub, iBlock, attemptIdx)
% LOCAL_MAKE_EYELINK_BLOCK_EDF_BASENAME Build <=8-char per-block EDF basename.
% Format: MBssbbaa, where ss=subject mod 100, bb=block mod 100, aa=attempt mod 100.
subjectPart = mod(round(double(iSub)), 100);
blockPart = mod(round(double(iBlock)), 100);
attemptPart = mod(round(double(attemptIdx)), 100);
baseName = sprintf('MB%02d%02d%02d', subjectPart, blockPart, attemptPart);
end

function local_eyelink_start_trial_recording(Conf, iSub, currentBlock, currentRun, ...
        iTrial, srcIdx, isReplayTrial, replayInstance, screenNumber)
% LOCAL_EYELINK_START_TRIAL_RECORDING Emit trial metadata and start recording.
if ~Conf.useEyelink || ~Conf.eyeTrackerInitialized
    return;
end
local_eyelink_safe_message(sprintf('TRIALID %d', iTrial));
try
    Eyelink('command', ['record_status_message "SUBJECT = %d ; RUN = %d ; ' ...
        'BLOCK = %d ; TRIAL = %d ; SOURCE = %d ; REPLAY = %d ; RINST = %d"'], ...
        iSub, currentRun, currentBlock, iTrial, srcIdx, double(isReplayTrial), replayInstance);
catch
    % Keep runtime robust if status message format is unsupported.
end
WaitSecs(max(0, double(Conf.eyeTrackerWaitSec)));
try
    Eyelink('Command', 'set_idle_mode');
catch
end
if Conf.eyeTrackerImageTransferEnabled
    imgFile = Conf.eyeTrackerFixBmpPath;
    if exist(imgFile, 'file') == 2
        try
            imgInfo = imfinfo(imgFile);
            [screenW, screenH] = Screen('WindowSize', screenNumber);
            trackerX = round(screenW / 2 - imgInfo.Width / 2);
            trackerY = round(screenH / 2 - imgInfo.Height / 2);
            Eyelink('ImageTransfer', imgFile, 0, 0, imgInfo.Width, imgInfo.Height, trackerX, trackerY, 1);
        catch
            % Non-fatal in runtime.
        end
    end
end
WaitSecs(max(0, double(Conf.eyeTrackerWaitSec)));
Eyelink('StartRecording', 1, 1, 1, 1);
WaitSecs(max(0, double(Conf.eyeTrackerWaitSec)));
if isReplayTrial
    local_eyelink_safe_message(sprintf('REPLAY_START source %d instance %d', srcIdx, replayInstance));
end
end

function local_eyelink_stop_trial_recording(Conf)
% LOCAL_EYELINK_STOP_TRIAL_RECORDING Stop EyeLink recording for one trial.
if ~Conf.useEyelink || ~Conf.eyeTrackerInitialized
    return;
end
try
    Eyelink('StopRecording');
catch
    % Keep stop best-effort.
end
end

function local_eyelink_safe_message(msg)
% LOCAL_EYELINK_SAFE_MESSAGE Best-effort EyeLink message emission.
try
    Eyelink('Message', char(msg));
catch
    % Keep messaging best-effort to avoid runtime aborts.
end
end

function [hasSample, eyeWithinWindow] = local_sample_eyelink_fixation(Conf, xCenter, yCenter)
% LOCAL_SAMPLE_EYELINK_FIXATION Read newest sample and test fixation window.
hasSample = false;
eyeWithinWindow = false;
try
    sample = Eyelink('NewestFloatSample');
catch
    sample = [];
end
if isempty(sample) || ~isstruct(sample)
    return;
end
gazeX = double(sample.gx);
gazeY = double(sample.gy);
pupilA = double(sample.pa);
validEyes = (~isnan(gazeX)) & (~isnan(gazeY)) & (pupilA > 0);
hasSample = true;
if any(validEyes)
    distFromFix = sqrt((gazeX - xCenter).^2 + (gazeY - yCenter).^2);
    eyeWithinWindow = any(distFromFix(validEyes) <= Conf.fixWindowPx);
end
end

function [fakeBreakThisTrial, fakeBreakStartFrame, fakeBreakEndFrame] = ...
        local_plan_fake_gaze_break(Conf, framesPerTrial)
% LOCAL_PLAN_FAKE_GAZE_BREAK Build one fake break schedule for debug mode.
fakeBreakThisTrial = false;
fakeBreakStartFrame = NaN;
fakeBreakEndFrame = NaN;
if ~Conf.useFakeEyeTracker || ~Conf.enableFixationAbort
    return;
end
if framesPerTrial < Conf.fixBreakToleranceFrames
    return;
end
fakeBreakThisTrial = rand < Conf.fakeGazeBreakRate;
if fakeBreakThisTrial
    lastStart = framesPerTrial - Conf.fixBreakToleranceFrames + 1;
    fakeBreakStartFrame = randi([1, lastStart]);
    fakeBreakEndFrame = fakeBreakStartFrame + Conf.fixBreakToleranceFrames - 1;
end
end

function local_present_fixation_break_warning(window, Conf)
% LOCAL_PRESENT_FIXATION_BREAK_WARNING Show warning screen after fixation break.
if Conf.fixationWarningDurationSec <= 0
    return;
end
Screen('FillRect', window, Conf.background);
DrawFormattedText(window, 'Please keep your fixation on the cross.', ...
    'center', 'center', Conf.fixColor);
Screen('Flip', window);
WaitSecs(max(0, double(Conf.fixationWarningDurationSec)));
end

function dotColor = local_dot_color_for_run(Conf, subjectNum, runNum)
% LOCAL_DOT_COLOR_FOR_RUN Resolve run-family color cue with odd/even counterbalancing.
%
% Rules when Conf.runColorCueEnabled is true:
%   odd subject:  run 1 -> dotColorRun1Base, runs 2/3 -> dotColorRuns23Base
%   even subject: run 1 -> dotColorRuns23Base, runs 2/3 -> dotColorRun1Base
%
% If disabled, use Conf.dotColor for all runs to preserve legacy behavior.
if ~isfield(Conf, 'runColorCueEnabled') || ~logical(Conf.runColorCueEnabled)
    dotColor = Conf.dotColor;
    return;
end

isOddSubject = mod(round(double(subjectNum)), 2) == 1;
if round(double(runNum)) == 1
    if isOddSubject
        dotColor = Conf.dotColorRun1Base;
    else
        dotColor = Conf.dotColorRuns23Base;
    end
else
    if isOddSubject
        dotColor = Conf.dotColorRuns23Base;
    else
        dotColor = Conf.dotColorRun1Base;
    end
end
end

function sequenceCodes = local_extract_source_sequence_trigger_codes(xyAll)
% LOCAL_EXTRACT_SOURCE_SEQUENCE_TRIGGER_CODES Read per-source sequence trigger values.
%
% Data flow:
%   xySeqs(source).sequence -> sequenceCodes(source) -> emitted trigger value.
sequenceCodes = nan(numel(xyAll), 1);
for iSrc = 1:numel(xyAll)
    if ~isfield(xyAll(iSrc), 'sequence') || isempty(xyAll(iSrc).sequence)
        error('Trial %d missing sequence field required for sequence-identity trigger.', iSrc);
    end
    seqValue = double(xyAll(iSrc).sequence);
    if ~isfinite(seqValue) || seqValue < 1 || mod(seqValue, 1) ~= 0
        error('Trial %d has invalid sequence value (%g); expected positive integer.', iSrc, seqValue);
    end
    sequenceCodes(iSrc) = seqValue;
end
end

function local_validate_sequence_trigger_range(sequenceRange, triggerMap)
% LOCAL_VALIDATE_SEQUENCE_TRIGGER_RANGE Ensure dynamic sequence range has no conflicts.
%
% Motivation:
%   Sequence-identity triggers occupy low values (typically 1..20). This
%   guard prevents overlap with fixed condition/event/catch trigger codes.
fixedValues = unique(cell2mat(struct2cell(triggerMap)));
fixedValues = fixedValues(isfinite(fixedValues));
overlap = intersect(unique(double(sequenceRange(:)')), unique(double(fixedValues(:)')));
if ~isempty(overlap)
    error(['Sequence trigger range %d..%d overlaps fixed trigger codes %s. ' ...
        'Lower trialsPerCondition or remap fixed trigger codes.'], ...
        min(sequenceRange), max(sequenceRange), mat2str(overlap));
end
end

function itiRangeSec = local_resolve_iti_range(itiRangeOverride, defaultRangeSec)
% LOCAL_RESOLVE_ITI_RANGE Validate runtime ITI jitter range in seconds.
%
% Data flow:
%   workspace override ITIRangeSec (or default) -> normalized [min max].
if nargin < 2 || isempty(defaultRangeSec)
    defaultRangeSec = [0.500 1.000];
end
if isempty(itiRangeOverride)
    itiRangeSec = double(defaultRangeSec(:)');
else
    itiRangeSec = double(itiRangeOverride(:)');
end
if numel(itiRangeSec) ~= 2 || any(~isfinite(itiRangeSec))
    error('ITIRangeSec must be a finite 1x2 numeric vector [min max] in seconds.');
end
if any(itiRangeSec < 0)
    error('ITIRangeSec values must be >= 0.');
end
if itiRangeSec(1) > itiRangeSec(2)
    error('ITIRangeSec must satisfy ITIRangeSec(1) <= ITIRangeSec(2).');
end
end

function [itiPrecomputedSec, itiSeed] = local_precompute_run_iti_seconds(iSub, blockNum, runNum, nBaseTrials, itiRangeSec)
% LOCAL_PRECOMPUTE_RUN_ITI_SECONDS Precompute replay-stable ITI values for one run.
%
% Data flow:
%   subject/block/run identifiers + range -> deterministic seed -> ITI vector.
nBaseTrials = round(double(nBaseTrials));
if ~isfinite(nBaseTrials) || nBaseTrials < 0
    error('nBaseTrials must be a finite integer >= 0.');
end
itiSeed = local_derive_iti_seed(iSub, blockNum, runNum);
if nBaseTrials == 0
    itiPrecomputedSec = zeros(1, 0);
    return;
end
itiMin = double(itiRangeSec(1));
itiMax = double(itiRangeSec(2));
stream = RandStream('mt19937ar', 'Seed', itiSeed);
itiPrecomputedSec = itiMin + (itiMax - itiMin) .* rand(stream, 1, nBaseTrials);
end

function itiSeed = local_derive_iti_seed(iSub, blockNum, runNum)
% LOCAL_DERIVE_ITI_SEED Derive deterministic ITI seed from participant identity.
%
% Assumption:
%   Participant ID is the primary seed source; block/run offsets keep seeds
%   deterministic per block/run for rescueTraject sessions.
subjectSeed = max(0, round(double(iSub)));
blockSeed = max(0, round(double(blockNum)));
runSeed = max(0, round(double(runNum)));
itiSeed = subjectSeed * 10000 + blockSeed * 100 + runSeed;
itiSeed = mod(itiSeed, 2^31 - 1);
if itiSeed <= 0
    itiSeed = max(1, subjectSeed + blockSeed + runSeed + 1);
end
end

function local_log_iti_precompute(prefixLabel, iSub, blockNum, runNum, itiSeed, itiRangeSec, itiPrecomputedSec)
% LOCAL_LOG_ITI_PRECOMPUTE Print ITI precompute diagnostics for dry-run/runtime logs.
prefixLabel = char(string(prefixLabel));
if isempty(itiPrecomputedSec)
    itiMin = NaN;
    itiMax = NaN;
else
    itiMin = min(itiPrecomputedSec);
    itiMax = max(itiPrecomputedSec);
end
fprintf(['%s ITI precompute: sub=%02d block=%d run=%d seed=%d ' ...
    'range=[%.3f %.3f] n=%d min=%.4f max=%.4f\n'], ...
    prefixLabel, iSub, blockNum, runNum, itiSeed, ...
    itiRangeSec(1), itiRangeSec(2), numel(itiPrecomputedSec), itiMin, itiMax);
fprintf('%s ITI vector (sec): %s\n', prefixLabel, mat2str(itiPrecomputedSec, 4));
end

function local_draw_arena(window, rectCoords, Conf)
% LOCAL_DRAW_ARENA Draw the arena border and interior rectangle.
marginRect = rectCoords + [-2 -2 2 2];
Screen('FillRect', window, Conf.rectBorderColor, marginRect);
Screen('FillRect', window, Conf.rectColor, rectCoords);
end

function local_draw_fixation(window, xCenter, yCenter, fixSizePx, fixColor)
% LOCAL_DRAW_FIXATION Draw a simple cross fixation marker.
lineWidth = max(1, floor(fixSizePx / 3));
xCenter = double(xCenter);
yCenter = double(yCenter);
fixSizePx = double(fixSizePx);
fixColor = double(fixColor(:)');
Screen('DrawLine', window, fixColor, xCenter - fixSizePx, yCenter, ...
    xCenter + fixSizePx, yCenter, lineWidth);
Screen('DrawLine', window, fixColor, xCenter, yCenter - fixSizePx, ...
    xCenter, yCenter + fixSizePx, lineWidth);
end

function local_present_post_trial_interval(window, Conf, xCenter, yCenter, fixSizePx, itiDurationSec)
% LOCAL_PRESENT_POST_TRIAL_INTERVAL Present one precomputed post-trial ITI.
%
% Sequence (current):
%   1) Blank (background only) with fixation for itiDurationSec.
%
% Data flow:
%   precomputed run-level ITI vector -> trial source-index mapping ->
%   trial-specific itiDurationSec -> immediate blank flip (dot removed) -> hold.
%
% Legacy sequence (disabled, kept for quick restore):
%   1) Blank (background only) for Conf.postTrialBlank1Sec
%   2) Low-luminance grid mask for Conf.postTrialGridSec
%   3) Blank (background only) for Conf.postTrialBlank2Sec
if ~(isfinite(itiDurationSec) && itiDurationSec >= 0)
    itiDurationSec = double(Conf.postTrialBlank1Sec);
end
local_draw_post_trial_blank(window, Conf, xCenter, yCenter, fixSizePx);
Screen('Flip', window);
WaitSecs(max(0, double(itiDurationSec)));
% local_draw_low_luminance_grid_mask(window, Conf, xCenter, yCenter, fixSizePx);
% Screen('Flip', window);
% WaitSecs(max(0, double(Conf.postTrialGridSec)));
% local_draw_post_trial_blank(window, Conf, xCenter, yCenter, fixSizePx);
% Screen('Flip', window);
% WaitSecs(max(0, double(Conf.postTrialBlank2Sec)));
end

function local_draw_post_trial_blank(window, Conf, xCenter, yCenter, fixSizePx)
% LOCAL_DRAW_POST_TRIAL_BLANK Draw an empty post-trial screen.
Screen('FillRect', window, Conf.background);
local_draw_fixation(window, xCenter, yCenter, fixSizePx, Conf.fixColor);
end

function local_draw_low_luminance_grid_mask(window, Conf, xCenter, yCenter, fixSizePx) %#ok<DEFNU>
% LOCAL_DRAW_LOW_LUMINANCE_GRID_MASK Draw a low-luminance fullscreen grid mask.
%
% Assumption:
%   The visual masking pattern is intentionally subtle to reduce abrupt
%   luminance transitions between trials during debugging and data collection.
%   This function is currently unused by active ITI flow (blank-only ITI).
Screen('FillRect', window, Conf.background);
windowRect = Screen('Rect', window);
Screen('FillRect', window, Conf.gridMaskBgColor, windowRect);

lineColor = uint8(round(double(Conf.gridMaskLineColor(:)')));
spacingPx = max(4, round(double(Conf.gridMaskSpacingPx)));
lineWidth = max(1, round(double(Conf.gridMaskLineWidthPx)));

xStart = floor(windowRect(1));
xEnd = ceil(windowRect(3));
yStart = floor(windowRect(2));
yEnd = ceil(windowRect(4));

for x = xStart:spacingPx:xEnd
    Screen('DrawLine', window, lineColor, x, yStart, x, yEnd, lineWidth);
end
for y = yStart:spacingPx:yEnd
    Screen('DrawLine', window, lineColor, xStart, y, xEnd, y, lineWidth);
end

% Keep fixation visible over the ITI mask.
local_draw_fixation(window, xCenter, yCenter, fixSizePx, Conf.fixColor);
end

function blendColor = local_blend_with_background(dotColor, bgColor, alpha)
% LOCAL_BLEND_WITH_BACKGROUND Blend dot color with background for alpha fallback.
blendColor = uint8(round(alpha .* double(dotColor) + (1 - alpha) .* double(bgColor)));
end

function frameIdx = local_clamp_frame(rawFrame, nFrames)
% LOCAL_CLAMP_FRAME Clamp frame indices to valid 1..nFrames range.
if isempty(rawFrame) || ~isfinite(rawFrame)
    frameIdx = 1;
else
    frameIdx = round(rawFrame);
    frameIdx = max(1, min(nFrames, frameIdx));
end
end

function trialOrder = local_resolve_block_run_order(rawTrialOrder, iBlock, iRun, nTrialsAvailable)
% LOCAL_RESOLVE_BLOCK_RUN_ORDER Return one block/run trial order from v21 TrialOrder.
%
% Inputs:
%   rawTrialOrder    : numeric TrialOrder array from MAT file.
%   iBlock           : requested block index (1-based).
%   iRun             : requested run index (1-based).
%   nTrialsAvailable : max valid trial index from xySeqs.
%
% Output:
%   trialOrder       : validated row vector of trial indices or [].
trialOrder = [];
if isempty(rawTrialOrder) || ~isnumeric(rawTrialOrder)
    return;
end

if ndims(rawTrialOrder) ~= 3
    return;
end

if iBlock < 1 || iBlock > size(rawTrialOrder, 1) || iRun < 1 || iRun > size(rawTrialOrder, 3)
    return;
end

candidate = squeeze(rawTrialOrder(iBlock, :, iRun));
candidate = candidate(:)';

if isempty(candidate)
    return;
end

candidate = round(double(candidate(:)'));
validMask = isfinite(candidate) & candidate >= 1 & candidate <= nTrialsAvailable;
candidate = candidate(validMask);
if isempty(candidate)
    return;
end

trialOrder = candidate;
end

function catchPlan = local_empty_catch_plan(rawTrialOrder)
% LOCAL_EMPTY_CATCH_PLAN Build an all-zero catch plan aligned to TrialOrder.
catchPlan = struct();
catchPlan.typeCode = zeros(size(rawTrialOrder));
catchPlan.expectedResponseCode = zeros(size(rawTrialOrder));
catchPlan.branchChangedPath = nan(size(rawTrialOrder));
catchPlan.disappearFrame = nan(size(rawTrialOrder));
catchPlan.reappearFrame = nan(size(rawTrialOrder));
catchPlan.altSourceTrialIndex = nan(size(rawTrialOrder));
end

function catchPlan = local_normalize_catch_plan(rawCatchPlan, rawTrialOrder)
% LOCAL_NORMALIZE_CATCH_PLAN Validate and normalize CatchPlan fields.
catchPlan = local_empty_catch_plan(rawTrialOrder);
if ~isstruct(rawCatchPlan)
    warning('CatchPlan is not a struct; using no-catch fallback.');
    return;
end

fieldNames = {'typeCode', 'expectedResponseCode', 'branchChangedPath', ...
    'disappearFrame', 'reappearFrame', 'altSourceTrialIndex'};
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    if isfield(rawCatchPlan, fieldName)
        candidate = rawCatchPlan.(fieldName);
        if isequal(size(candidate), size(rawTrialOrder))
            catchPlan.(fieldName) = double(candidate);
        else
            warning('CatchPlan.%s size mismatch; using fallback defaults.', fieldName);
        end
    end
end
end

function catchVec = local_resolve_block_run_catch(rawCatchPlane, rawTrialOrder, iBlock, iRun, nTrialsAvailable, defaultValue)
% LOCAL_RESOLVE_BLOCK_RUN_CATCH Resolve one catch metadata vector aligned to run order.
catchVec = [];
if isempty(rawCatchPlane) || isempty(rawTrialOrder)
    return;
end
if ndims(rawTrialOrder) ~= 3 || ndims(rawCatchPlane) ~= 3
    return;
end
if iBlock < 1 || iBlock > size(rawTrialOrder, 1) || iRun < 1 || iRun > size(rawTrialOrder, 3)
    return;
end

trialCandidate = squeeze(rawTrialOrder(iBlock, :, iRun));
catchCandidate = squeeze(rawCatchPlane(iBlock, :, iRun));
trialCandidate = double(trialCandidate(:)');
catchCandidate = double(catchCandidate(:)');

validMask = isfinite(trialCandidate) & trialCandidate >= 1 & trialCandidate <= nTrialsAvailable;
catchVec = catchCandidate(validMask);
if isempty(catchVec)
    return;
end
replaceMask = ~isfinite(catchVec);
if any(replaceMask)
    catchVec(replaceMask) = defaultValue;
end
end

function [trialOrderOut, catchTypeOut, catchExpectedOut, catchBranchOut, catchDisappearOut, catchReappearOut, catchAltSourceOut] = ...
        local_apply_practice_mode_to_run( ...
        trialOrderIn, catchTypeIn, catchExpectedIn, catchBranchIn, ...
        catchDisappearIn, catchReappearIn, catchAltSourceIn, ...
        currentRun, xyAll, inputFps, framesPerTrial, ...
        practiceModeEnabled, practiceRunPercentageTrials, Conf)
% LOCAL_APPLY_PRACTICE_MODE_TO_RUN Optionally trim run length and synthesize practice catches.
%
% Data flow:
%   scheduled trial/catch vectors + practice controls -> adjusted run vectors.
trialOrderOut = double(trialOrderIn(:)');
catchTypeOut = double(catchTypeIn(:)');
catchExpectedOut = double(catchExpectedIn(:)');
catchBranchOut = double(catchBranchIn(:)');
catchDisappearOut = double(catchDisappearIn(:)');
catchReappearOut = double(catchReappearIn(:)');
catchAltSourceOut = double(catchAltSourceIn(:)');

if isempty(trialOrderOut)
    return;
end

if ~logical(practiceModeEnabled)
    return;
end

nTrials = numel(trialOrderOut);
nKeep = round((double(practiceRunPercentageTrials) / 100) * nTrials);
nKeep = max(1, min(nTrials, nKeep));

if nKeep < nTrials
    selectedIdx = randperm(nTrials, nKeep);
    selectedIdx = sort(selectedIdx);
    trialOrderOut = trialOrderOut(selectedIdx);
end

if round(double(currentRun)) == 1
    % Practice rule for run 1 (always_visible): force type-1 catches.
    %
    % Critical constraint:
    %   For changed-path catches, align the hidden window so that the
    %   deviance frame occurs while the dot is invisible.
    catchTypeOut = ones(1, nKeep);
    catchBranchOut = zeros(1, nKeep);
    catchDisappearOut = nan(1, nKeep);
    catchReappearOut = nan(1, nKeep);
    catchAltSourceOut = nan(1, nKeep);
    catchExpectedOut = Conf.catchResponseNoCode * ones(1, nKeep);

    fps = max(1, round(double(inputFps)));
    disappearRangeFrames = max(1, round(double(Conf.catchType1DisappearRangeSec(:)') * fps));
    disappearMinFrame = max(1, min(disappearRangeFrames));
    disappearMaxFrame = min(max(1, framesPerTrial - 1), max(disappearRangeFrames));
    if disappearMaxFrame < disappearMinFrame
        disappearMinFrame = 1;
        disappearMaxFrame = max(1, framesPerTrial - 1);
    end
    invisibleFrames = max(1, round(double(Conf.catchType1InvisibleDurationSec) * fps));
    changeProb = max(0, min(1, double(Conf.catchType1ChangedPathProbability)));
    branchMask = rand(1, nKeep) < changeProb;

    for iKeep = 1:nKeep
        sourceIdx = trialOrderOut(iKeep);
        devianceTargetFrame = NaN;

        if branchMask(iKeep)
            altSourceIdx = local_pick_practice_alt_source(sourceIdx, trialOrderOut, xyAll);
            if isfinite(altSourceIdx) && altSourceIdx ~= sourceIdx
                catchBranchOut(iKeep) = 1;
                catchAltSourceOut(iKeep) = altSourceIdx;
                catchExpectedOut(iKeep) = Conf.catchResponseYesCode;
                devianceTargetFrame = local_pick_type1_practice_deviance_frame( ...
                    sourceIdx, altSourceIdx, xyAll, framesPerTrial);
            end
        end

        [disappearFrame, reappearFrame] = local_pick_type1_practice_hidden_window( ...
            disappearMinFrame, disappearMaxFrame, invisibleFrames, framesPerTrial, devianceTargetFrame);
        catchDisappearOut(iKeep) = disappearFrame;
        catchReappearOut(iKeep) = reappearFrame;
    end
else
    % Practice rule for runs 2/3 (occluded): force type-2 catches.
    catchTypeOut = 2 * ones(1, nKeep);
    catchBranchOut = nan(1, nKeep);
    catchDisappearOut = nan(1, nKeep);
    catchReappearOut = nan(1, nKeep);
    catchAltSourceOut = nan(1, nKeep);
    catchExpectedOut = zeros(1, nKeep);

    for iKeep = 1:nKeep
        catchExpectedOut(iKeep) = local_expected_code_for_type2_condition( ...
            trialOrderOut(iKeep), xyAll, Conf.catchResponseYesCode, Conf.catchResponseNoCode);
    end
end
end

function altSourceIdx = local_pick_practice_alt_source(sourceIdx, runTrialOrder, xyAll)
% LOCAL_PICK_PRACTICE_ALT_SOURCE Pick an alternate source index for run1 changed-path practice catches.
altSourceIdx = NaN;
if sourceIdx < 1 || sourceIdx > numel(xyAll)
    return;
end

candidates = [];
if isfield(xyAll(sourceIdx), 'sequence') && ~isempty(xyAll(sourceIdx).sequence)
    sourceSeq = double(xyAll(sourceIdx).sequence);
    for iTrial = 1:numel(xyAll)
        if iTrial == sourceIdx
            continue;
        end
        if ~isfield(xyAll(iTrial), 'sequence') || isempty(xyAll(iTrial).sequence)
            continue;
        end
        if double(xyAll(iTrial).sequence) ~= sourceSeq
            continue;
        end
        if ~isfield(xyAll(iTrial), 'condition_label')
            continue;
        end
        condLabel = strtrim(char(string(xyAll(iTrial).condition_label)));
        if strcmpi(condLabel, 'occluded_deviant')
            candidates(end + 1) = iTrial; %#ok<AGROW>
        end
    end
end

if isempty(candidates)
    runCandidates = runTrialOrder(runTrialOrder ~= sourceIdx);
    runCandidates = runCandidates(isfinite(runCandidates));
    candidates = double(runCandidates(:)');
end

if isempty(candidates)
    return;
end
altSourceIdx = candidates(randi(numel(candidates), 1, 1));
end

function devianceFrame = local_pick_type1_practice_deviance_frame(sourceIdx, altSourceIdx, xyAll, framesPerTrial)
% LOCAL_PICK_TYPE1_PRACTICE_DEVIANCE_FRAME Resolve deviance frame for changed-path run1 catches.
%
% Priority:
%   1) Alternate source trial deviance frame (preferred for changed path)
%   2) Original source trial deviance frame (fallback)
%
% Output:
%   NaN when no valid frame is available.
devianceFrame = local_read_trial_deviance_frame(altSourceIdx, xyAll, framesPerTrial);
if ~isfinite(devianceFrame)
    devianceFrame = local_read_trial_deviance_frame(sourceIdx, xyAll, framesPerTrial);
end
end

function devianceFrame = local_read_trial_deviance_frame(trialIdx, xyAll, framesPerTrial)
% LOCAL_READ_TRIAL_DEVIANCE_FRAME Read and clamp one trial deviance frame.
devianceFrame = NaN;
if ~isfinite(trialIdx)
    return;
end
trialIdx = round(double(trialIdx));
if trialIdx < 1 || trialIdx > numel(xyAll)
    return;
end
if ~isfield(xyAll(trialIdx), 'deviance_frame') || isempty(xyAll(trialIdx).deviance_frame)
    return;
end

rawFrame = double(xyAll(trialIdx).deviance_frame);
if ~isfinite(rawFrame)
    return;
end
devianceFrame = round(rawFrame);
devianceFrame = max(1, min(framesPerTrial, devianceFrame));
end

function [disappearFrame, reappearFrame] = local_pick_type1_practice_hidden_window( ...
        disappearMinFrame, disappearMaxFrame, invisibleFrames, framesPerTrial, devianceTargetFrame)
% LOCAL_PICK_TYPE1_PRACTICE_HIDDEN_WINDOW Pick hide/reappear frames for run1 type-1 practice catches.
%
% Rule for changed-path catches:
%   If a deviance target frame is provided, enforce:
%     disappearFrame <= devianceTargetFrame < reappearFrame
%   so the path change happens while the dot is invisible.
%
% Fallback:
%   If constraints cannot be met inside the configured disappearance range,
%   the function relaxes that range to keep the deviance hidden.
disappearFrame = randi([disappearMinFrame, disappearMaxFrame], 1, 1);

if isfinite(devianceTargetFrame)
    targetFrame = round(double(devianceTargetFrame));
    targetFrame = max(1, min(framesPerTrial - 1, targetFrame));

    inRangeLower = max(disappearMinFrame, targetFrame - invisibleFrames + 1);
    inRangeUpper = min(disappearMaxFrame, targetFrame);
    if inRangeLower <= inRangeUpper
        disappearFrame = randi([inRangeLower, inRangeUpper], 1, 1);
    else
        relaxedLower = max(1, targetFrame - invisibleFrames + 1);
        relaxedUpper = min(framesPerTrial - 1, targetFrame);
        if relaxedLower <= relaxedUpper
            disappearFrame = randi([relaxedLower, relaxedUpper], 1, 1);
        else
            disappearFrame = max(1, min(framesPerTrial - 1, targetFrame));
        end
    end
end

reappearFrame = min(framesPerTrial, disappearFrame + invisibleFrames);
if reappearFrame <= disappearFrame
    reappearFrame = min(framesPerTrial, disappearFrame + 1);
end

if isfinite(devianceTargetFrame)
    targetFrame = round(double(devianceTargetFrame));
    targetFrame = max(1, min(framesPerTrial - 1, targetFrame));
    if ~(disappearFrame <= targetFrame && targetFrame < reappearFrame)
        reappearFrame = min(framesPerTrial, max(reappearFrame, targetFrame + 1));
    end
end
end

function expectedCode = local_expected_code_for_type2_condition(sourceIdx, xyAll, yesCode, noCode)
% LOCAL_EXPECTED_CODE_FOR_TYPE2_CONDITION Resolve expected answer for type-2 practice catches.
expectedCode = noCode;
if sourceIdx < 1 || sourceIdx > numel(xyAll)
    return;
end
if ~isfield(xyAll(sourceIdx), 'condition_label')
    return;
end
condLabel = strtrim(char(string(xyAll(sourceIdx).condition_label)));
if strcmpi(condLabel, 'occluded_deviant')
    expectedCode = yesCode;
end
end

function [scoredCount, answeredCount, correctCount] = local_count_run_catch_accuracy(output)
% LOCAL_COUNT_RUN_CATCH_ACCURACY Count scored/answered/correct catches in one run output.
%
% Definition:
%   scoredCount   = catches with a defined expected response code (>0)
%   answeredCount = scored catches with a submitted response code (>0)
%   correctCount  = scored catches marked as correct
%
% Notes:
%   - Timeouts contribute to scoredCount but not answeredCount.
%   - Accuracy denominator uses scoredCount so timeouts reduce accuracy.
scoredCount = 0;
answeredCount = 0;
correctCount = 0;
if isempty(output)
    return;
end

expectedVec = nan(1, numel(output));
responseVec = nan(1, numel(output));
correctVec = nan(1, numel(output));
for iTrial = 1:numel(output)
    if isfield(output(iTrial), 'catch_expected_response_code')
        expectedVec(iTrial) = local_numeric_scalar_or_nan(output(iTrial).catch_expected_response_code);
    end
    if isfield(output(iTrial), 'catch_response_code')
        responseVec(iTrial) = local_numeric_scalar_or_nan(output(iTrial).catch_response_code);
    end
    if isfield(output(iTrial), 'catch_response_correct')
        correctVec(iTrial) = local_numeric_scalar_or_nan(output(iTrial).catch_response_correct);
    end
end

scoredMask = isfinite(expectedVec) & (expectedVec > 0);
scoredCount = sum(scoredMask);
if scoredCount == 0
    return;
end

answeredCount = sum(scoredMask & isfinite(responseVec) & responseVec > 0);
correctMask = scoredMask & isfinite(correctVec);
if any(correctMask)
    correctCount = sum(correctVec(correctMask) > 0.5);
end
end

function valueOut = local_numeric_scalar_or_nan(valueIn)
% LOCAL_NUMERIC_SCALAR_OR_NAN Read one numeric scalar from optional fields.
%
% Data flow:
%   struct field (possibly empty/vector) -> scalar numeric or NaN.
%   This prevents assignment-size errors when legacy/default fields are [].
valueOut = NaN;
if isempty(valueIn) || ~(isnumeric(valueIn) || islogical(valueIn))
    return;
end
valueIn = double(valueIn);
if isempty(valueIn)
    return;
end
valueOut = valueIn(1);
if ~isfinite(valueOut)
    valueOut = NaN;
end
end

function local_print_run_accuracy_summary(runCatchCorrect, runCatchScored, runCatchAnswered)
% LOCAL_PRINT_RUN_ACCURACY_SUMMARY Print catch accuracy aggregated by run index.
fprintf('\nCatch accuracy summary by run:\n');
for iRun = 1:numel(runCatchScored)
    scored = runCatchScored(iRun);
    correct = runCatchCorrect(iRun);
    answered = runCatchAnswered(iRun);
    timedOut = max(0, scored - answered);
    if scored > 0
        accuracyPct = 100 * (correct / scored);
        fprintf(['  Run %d: %d/%d correct (%.1f%%), answered %d, ' ...
            'timeouts %d.\n'], iRun, correct, scored, accuracyPct, answered, timedOut);
    else
        fprintf('  Run %d: no scored catch responses.\n', iRun);
    end
end
end

function value = local_pick_from_vector(vec, idx, defaultValue)
% LOCAL_PICK_FROM_VECTOR Safe scalar extraction with default fallback.
value = defaultValue;
if idx <= 0 || idx > numel(vec)
    return;
end
candidate = vec(idx);
if isfinite(candidate)
    value = candidate;
end
end

function label = local_catch_type_label(typeCode)
% LOCAL_CATCH_TYPE_LABEL Human-readable label for catch type code.
switch round(typeCode)
    case 1
        label = 'always_visible_disappear_reappear';
    case 2
        label = 'occlusion_question_only';
    otherwise
        label = 'none';
end
end

function [trialStim, visibleMask, disappearFrame, reappearFrame] = ...
        local_prepare_type1_catch_stimulus(baseStim, altStim, branchChangedRaw, disappearRaw, reappearRaw)
% LOCAL_PREPARE_TYPE1_CATCH_STIMULUS Build one type-1 catch stimulus trajectory.
nFrames = size(baseStim, 1);
trialStim = baseStim;
visibleMask = true(nFrames, 1);

disappearFrame = local_clamp_frame(disappearRaw, nFrames);
reappearFrame = local_clamp_frame(reappearRaw, nFrames);
if reappearFrame <= disappearFrame
    reappearFrame = min(nFrames, disappearFrame + 1);
end

hiddenStart = disappearFrame;
hiddenEnd = max(hiddenStart, reappearFrame - 1);
visibleMask(hiddenStart:hiddenEnd) = false;

branchChanged = isfinite(branchChangedRaw) && round(branchChangedRaw) == 1;
if branchChanged && isequal(size(altStim), size(baseStim))
    trialStim(reappearFrame:end, :) = altStim(reappearFrame:end, :);
end
end

function [trialOut, debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel, abortFromQuestion] = ...
        local_run_catch_question(trialOut, window, Conf, trialNum, xCenter, yCenter, fixSizePx, ...
        debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel)
% LOCAL_RUN_CATCH_QUESTION Show catch prompt and collect response (MEG + keyboard).
%
% Data flow:
%   prompt onset -> response polling loop -> one end-status trigger.
%   Exactly one of correct/incorrect/timeout is emitted per catch question.
abortFromQuestion = false;
trialOut.catch_question_presented = true;

questionStart = GetSecs;
trialOut.catch_question_start_time = questionStart;
trialOut.catch_response_code = 0;
trialOut.catch_response_label = 'none';
trialOut.catch_response_rt_sec = NaN;
trialOut.catch_response_device = '';
trialOut.catch_timed_out = false;

[debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
    local_emit_trigger(Conf, Conf.trigger.catchQuestionStart, trialNum, 0, ...
    'catch_question_start', debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);

local_start_meg_response_log(Conf);
cleanupObj = onCleanup(@() local_stop_meg_response_log(Conf));

deadline = questionStart + Conf.catchResponseTimeoutSec;
while GetSecs < deadline
    % Keep fixation continuously visible while showing top/bottom question text.
    local_draw_catch_question_screen(window, Conf, xCenter, yCenter, fixSizePx);
    Screen('Flip', window);

    [responseCode, responseDevice, responseTime, escapePressed] = local_poll_catch_response(Conf);
    if escapePressed
        abortFromQuestion = true;
        break;
    end
    if responseCode > 0
        trialOut.catch_response_code = responseCode;
        trialOut.catch_response_label = local_response_label(responseCode);
        trialOut.catch_response_rt_sec = responseTime - questionStart;
        trialOut.catch_response_device = responseDevice;
        break;
    end
end

if ~abortFromQuestion
    % Resolve timeout status first; this state is mutually exclusive with
    % answered catches and drives which single end-status trigger is sent.
    if trialOut.catch_response_code == 0
        trialOut.catch_timed_out = true;
    end

    expectedCode = trialOut.catch_expected_response_code;
    if isfinite(expectedCode) && expectedCode > 0 && trialOut.catch_response_code > 0
        trialOut.catch_response_correct = (round(expectedCode) == round(trialOut.catch_response_code));
    else
        trialOut.catch_response_correct = NaN;
    end

    if trialOut.catch_timed_out
        endStatusTrigger = Conf.trigger.catchQuestionEndTimeout;
        endStatusLabel = 'catch_question_timeout';
    elseif isfinite(trialOut.catch_response_correct) && logical(trialOut.catch_response_correct)
        endStatusTrigger = Conf.trigger.catchQuestionEndCorrect;
        endStatusLabel = 'catch_question_correct';
    else
        endStatusTrigger = Conf.trigger.catchQuestionEndIncorrect;
        endStatusLabel = 'catch_question_incorrect';
    end

    trialOut.catch_question_end_time = GetSecs;
    [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
        local_emit_trigger(Conf, endStatusTrigger, trialNum, 0, ...
        endStatusLabel, debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel);

    % Present explicit feedback immediately after the question phase.
    local_present_catch_feedback(window, Conf, xCenter, yCenter, fixSizePx, ...
        trialOut.catch_timed_out, trialOut.catch_response_correct);
end
end

function [responseCode, responseDevice, responseTime, escapePressed] = local_poll_catch_response(Conf)
% LOCAL_POLL_CATCH_RESPONSE Poll one response event from keyboard or MEG button box.
responseCode = 0;
responseDevice = '';
responseTime = NaN;
escapePressed = false;

if Conf.enableKeyboardDebugResponse
    [keyDown, secs, keyCode] = KbCheck;
    if keyDown
        if keyCode(Conf.escapeKeyCode)
            escapePressed = true;
            KbReleaseWait;
            return;
        end
        if any(keyCode(Conf.catchYesKeyCodes))
            responseCode = 1;
            responseDevice = 'keyboard';
            responseTime = secs;
            KbReleaseWait;
            return;
        end
        if any(keyCode(Conf.catchNoKeyCodes))
            responseCode = 2;
            responseDevice = 'keyboard';
            responseTime = secs;
            KbReleaseWait;
            return;
        end
    end
end

if Conf.MEG
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    if status.newLogFrames > 0 && mod(status.currentWriteFrame, 2)
        [data, ~] = Datapixx('ReadDinLog');
        if ~isempty(data)
            rawValue = bitand(uint32(data(end)), Conf.megButtonMask);
            if rawValue == Conf.megButtonYesValue
                responseCode = 1;
                responseDevice = 'meg_button_box';
                responseTime = GetSecs;
                return;
            elseif rawValue == Conf.megButtonNoValue
                responseCode = 2;
                responseDevice = 'meg_button_box';
                responseTime = GetSecs;
                return;
            end
        end
    end
end
end

function local_start_meg_response_log(Conf)
% LOCAL_START_MEG_RESPONSE_LOG Prepare DataPixx Din logger for response capture.
if ~Conf.MEG
    return;
end
Datapixx('EnableDinDebounce');
Datapixx('SetDinLog');
Datapixx('StartDinLog');
Datapixx('RegWrRd');
end

function local_stop_meg_response_log(Conf)
% LOCAL_STOP_MEG_RESPONSE_LOG Stop DataPixx Din logger after response capture.
if ~Conf.MEG
    return;
end
Datapixx('StopDinLog');
Datapixx('RegWrRd');
end

function keyCodes = local_keynames_to_codes(keyNames)
% LOCAL_KEYNAMES_TO_CODES Convert key-name list to numeric KbCheck codes.
keyCodes = [];
for iName = 1:numel(keyNames)
    try
        code = KbName(keyNames{iName});
        if isnumeric(code)
            keyCodes = [keyCodes, code(:)']; %#ok<AGROW>
        end
    catch
        % Ignore invalid key names so runtime remains robust.
    end
end
keyCodes = unique(keyCodes(keyCodes > 0));
end

function label = local_response_label(responseCode)
% LOCAL_RESPONSE_LABEL Human-readable response label.
switch round(responseCode)
    case 1
        label = 'yes';
    case 2
        label = 'no';
    otherwise
        label = 'none';
end
end

function xyPx = local_xy_to_px(xyDeg, ppd, rectCoords)
% LOCAL_XY_TO_PX Convert one [x y] point from visual degrees to screen pixels.
xyPx = [xyDeg(1) * ppd + rectCoords(1), xyDeg(2) * ppd + rectCoords(2)];
end

function coordsPx = local_polyline_to_px(xyDeg, ppd, rectCoords)
% LOCAL_POLYLINE_TO_PX Convert Nx2 degree polyline to DrawLines segment pairs.
%
% DrawLines expects point pairs [p1 p2 p3 p4 ...], where each adjacent pair
% is one segment. Therefore we expand a polyline of N points into 2*(N-1)
% coordinates: [p1 p2 p2 p3 ... p(N-1) pN].
coordsPx = [];
if isempty(xyDeg) || size(xyDeg, 2) ~= 2
    return;
end

xyDeg = double(xyDeg);
finiteMask = isfinite(xyDeg(:, 1)) & isfinite(xyDeg(:, 2));
xyDeg = xyDeg(finiteMask, :);
if size(xyDeg, 1) < 2
    return;
end

xPx = xyDeg(:, 1) .* ppd + rectCoords(1);
yPx = xyDeg(:, 2) .* ppd + rectCoords(2);
coordsPx = local_polyline_px_to_drawlines_coords([xPx, yPx]);
end

function coordsPx = local_polyline_px_to_drawlines_coords(xyPx)
% LOCAL_POLYLINE_PX_TO_DRAWLINES_COORDS Expand Nx2 pixel polyline to PTB pair format.
coordsPx = [];
if isempty(xyPx) || size(xyPx, 2) ~= 2
    return;
end

xyPx = double(xyPx);
finiteMask = isfinite(xyPx(:, 1)) & isfinite(xyPx(:, 2));
xyPx = xyPx(finiteMask, :);
if size(xyPx, 1) < 2
    return;
end

xStart = xyPx(1:(end - 1), 1)';
xEnd = xyPx(2:end, 1)';
yStart = xyPx(1:(end - 1), 2)';
yEnd = xyPx(2:end, 2)';

coordsPx = [reshape([xStart; xEnd], 1, []); ...
    reshape([yStart; yEnd], 1, [])];
end

function drawPlan = local_precompute_polyline_draw_plan(coordsPx, lineWidthPx, drawOpts)
% LOCAL_PRECOMPUTE_POLYLINE_DRAW_PLAN Prepare per-trial occluder draw geometry.
%
% Data flow:
%   raw paired coords + style opts -> one reusable plan consumed per frame.
%   For straight style, polygon vertices are built once to avoid rebuilding
%   tangent/normal geometry inside the frame loop.
drawPlan = struct( ...
    'valid', false, ...
    'useStraightBandPolygon', false, ...
    'preparedCoordsPx', [], ...
    'lineWidthPx', [], ...
    'preparedDrawOpts', struct(), ...
    'straightBandPolygonPx', [], ...
    'straightStartBoostPolygonPx', []);

if isempty(coordsPx) || size(coordsPx, 1) ~= 2
    return;
end
if nargin < 3 || isempty(drawOpts)
    drawOpts = struct();
end

lineWidthPx = double(lineWidthPx);
if ~isfinite(lineWidthPx) || lineWidthPx <= 0
    return;
end
lineWidthPx = max(1, lineWidthPx);

startBackshiftPx = 0;
endForwardShiftPx = 0;
startWidthExtraPx = 0;
drawStartCap = true;
drawEndCap = true;

if isfield(drawOpts, 'startBackshiftPx') && ~isempty(drawOpts.startBackshiftPx)
    startBackshiftPx = double(drawOpts.startBackshiftPx);
end
if isfield(drawOpts, 'endForwardShiftPx') && ~isempty(drawOpts.endForwardShiftPx)
    endForwardShiftPx = double(drawOpts.endForwardShiftPx);
end
if isfield(drawOpts, 'startWidthExtraPx') && ~isempty(drawOpts.startWidthExtraPx)
    startWidthExtraPx = double(drawOpts.startWidthExtraPx);
end
if isfield(drawOpts, 'drawStartCap') && ~isempty(drawOpts.drawStartCap)
    drawStartCap = logical(drawOpts.drawStartCap);
end
if isfield(drawOpts, 'drawEndCap') && ~isempty(drawOpts.drawEndCap)
    drawEndCap = logical(drawOpts.drawEndCap);
end

if ~isfinite(startBackshiftPx) || startBackshiftPx < 0
    startBackshiftPx = 0;
end
if ~isfinite(endForwardShiftPx) || endForwardShiftPx < 0
    endForwardShiftPx = 0;
end
if ~isfinite(startWidthExtraPx) || startWidthExtraPx < 0
    startWidthExtraPx = 0;
end

preparedCoordsPx = local_apply_polyline_start_backshift(double(coordsPx), startBackshiftPx);
preparedCoordsPx = local_apply_polyline_end_forward_shift(preparedCoordsPx, endForwardShiftPx);
if mod(size(preparedCoordsPx, 2), 2) ~= 0
    preparedCoordsPx = preparedCoordsPx(:, 1:(end - 1));
end
if size(preparedCoordsPx, 2) < 2
    return;
end

preparedDrawOpts = drawOpts;
preparedDrawOpts.startBackshiftPx = 0;
preparedDrawOpts.endForwardShiftPx = 0;

drawPlan.valid = true;
drawPlan.preparedCoordsPx = preparedCoordsPx;
drawPlan.lineWidthPx = lineWidthPx;
drawPlan.preparedDrawOpts = preparedDrawOpts;

if ~drawStartCap && ~drawEndCap
    [straightBandPolygonPx, straightBoostPolygonPx, ok] = ...
        local_build_straight_polygons_from_coords(preparedCoordsPx, lineWidthPx, startWidthExtraPx);
    if ok
        drawPlan.useStraightBandPolygon = true;
        drawPlan.straightBandPolygonPx = straightBandPolygonPx;
        drawPlan.straightStartBoostPolygonPx = straightBoostPolygonPx;
    end
end
end

function local_draw_polyline_plan(window, drawPlan, lineColor)
% LOCAL_DRAW_POLYLINE_PLAN Draw one precomputed polyline plan.
if nargin < 2 || isempty(drawPlan) || ~isstruct(drawPlan) || ...
        ~isfield(drawPlan, 'valid') || ~logical(drawPlan.valid)
    return;
end

if isfield(drawPlan, 'useStraightBandPolygon') && logical(drawPlan.useStraightBandPolygon)
    if isfield(drawPlan, 'straightBandPolygonPx') && size(drawPlan.straightBandPolygonPx, 1) >= 3
        Screen('FillPoly', window, lineColor, drawPlan.straightBandPolygonPx, 1);
    end
    if isfield(drawPlan, 'straightStartBoostPolygonPx') && size(drawPlan.straightStartBoostPolygonPx, 1) >= 3
        Screen('FillPoly', window, lineColor, drawPlan.straightStartBoostPolygonPx, 1);
    end
    return;
end

local_draw_polyline_segments(window, drawPlan.preparedCoordsPx, ...
    drawPlan.lineWidthPx, lineColor, drawPlan.preparedDrawOpts);
end

function [bandPolygonPx, boostPolygonPx, ok] = ...
        local_build_straight_polygons_from_coords(coordsPx, lineWidthPx, startWidthExtraPx)
% LOCAL_BUILD_STRAIGHT_POLYGONS_FROM_COORDS Build cached polygons for straight style.
ok = false;
bandPolygonPx = zeros(0, 2);
boostPolygonPx = zeros(0, 2);

if isempty(coordsPx) || size(coordsPx, 1) ~= 2
    return;
end
if mod(size(coordsPx, 2), 2) ~= 0
    coordsPx = coordsPx(:, 1:(end - 1));
end
if size(coordsPx, 2) < 2
    return;
end

polylinePx = local_drawlines_coords_to_polyline(coordsPx, size(coordsPx, 2));
if size(polylinePx, 1) < 2
    return;
end

halfWidthPx = max(1, double(lineWidthPx)) / 2;
bandPolygonPx = local_build_straight_band_polygon_vertices(polylinePx, halfWidthPx);
if isempty(bandPolygonPx)
    return;
end
boostPolygonPx = local_build_straight_start_width_boost_polygon( ...
    polylinePx, halfWidthPx, startWidthExtraPx);
ok = true;
end

function polygonPx = local_build_straight_band_polygon_vertices(polylinePx, halfWidthPx)
% LOCAL_BUILD_STRAIGHT_BAND_POLYGON_VERTICES Build filled-band polygon vertices.
polygonPx = zeros(0, 2);
if isempty(polylinePx) || size(polylinePx, 1) < 2 || ~isfinite(halfWidthPx) || halfWidthPx <= 0
    return;
end

nPts = size(polylinePx, 1);
tangent = zeros(nPts, 2);
for iPt = 1:nPts
    if iPt == 1
        v = polylinePx(2, :) - polylinePx(1, :);
    elseif iPt == nPts
        v = polylinePx(nPts, :) - polylinePx(nPts - 1, :);
    else
        vPrev = polylinePx(iPt, :) - polylinePx(iPt - 1, :);
        vNext = polylinePx(iPt + 1, :) - polylinePx(iPt, :);
        nPrev = norm(vPrev);
        nNext = norm(vNext);
        if nPrev <= 1e-9 && nNext <= 1e-9
            v = [0, 0];
        elseif nPrev <= 1e-9
            v = vNext;
        elseif nNext <= 1e-9
            v = vPrev;
        else
            uPrev = vPrev ./ nPrev;
            uNext = vNext ./ nNext;
            v = uPrev + uNext;
            if norm(v) <= 1e-9
                v = uNext;
            end
        end
    end

    vNorm = norm(v);
    if vNorm <= 1e-9
        if iPt > 1 && norm(tangent(iPt - 1, :)) > 0
            tangent(iPt, :) = tangent(iPt - 1, :);
        else
            tangent(iPt, :) = [1, 0];
        end
    else
        tangent(iPt, :) = v ./ vNorm;
    end
end

normal = [-tangent(:, 2), tangent(:, 1)];
leftEdge = polylinePx + halfWidthPx .* normal;
rightEdge = polylinePx - halfWidthPx .* normal;
polygonPx = [leftEdge; flipud(rightEdge)];

finiteMask = isfinite(polygonPx(:, 1)) & isfinite(polygonPx(:, 2));
polygonPx = polygonPx(finiteMask, :);
if size(polygonPx, 1) < 3
    polygonPx = zeros(0, 2);
    return;
end
polyArea = polyarea(polygonPx(:, 1), polygonPx(:, 2));
if ~isfinite(polyArea) || polyArea <= 1e-3
    polygonPx = zeros(0, 2);
end
end

function boostPoly = local_build_straight_start_width_boost_polygon( ...
        polylinePx, halfWidthPx, startWidthExtraPx)
% LOCAL_BUILD_STRAIGHT_START_WIDTH_BOOST_POLYGON Build optional start-wall boost.
boostPoly = zeros(0, 2);
if ~isfinite(startWidthExtraPx) || startWidthExtraPx <= 0
    return;
end
if isempty(polylinePx) || size(polylinePx, 1) < 2
    return;
end

p1 = polylinePx(1, :);
tangent = [];
for iPt = 2:size(polylinePx, 1)
    d = polylinePx(iPt, :) - p1;
    dNorm = norm(d);
    if dNorm > 1e-9
        tangent = d ./ dNorm;
        break;
    end
end
if isempty(tangent)
    return;
end

normal = [-tangent(2), tangent(1)];
depthPx = max(1, halfWidthPx);
wideHalf = halfWidthPx + startWidthExtraPx;
p2 = p1 + depthPx .* tangent;
boostPoly = [p1 + wideHalf .* normal; ...
    p1 - wideHalf .* normal; ...
    p2 - halfWidthPx .* normal; ...
    p2 + halfWidthPx .* normal];
end

function local_draw_polyline_segments(window, coordsPx, lineWidthPx, lineColor, drawOpts)
% LOCAL_DRAW_POLYLINE_SEGMENTS Draw paired 2xN coordinates as a thick filled band.
%
% Rationale:
%   On some GPUs/PTB backends, DrawLine/DrawLines reject large line widths
%   (e.g., >16 px) even when the occluder must be much wider to fully cover
%   the dot. To avoid hardware line-width caps, each segment is rendered as
%   a filled oriented quad, then round end-caps are added with FillOval.
%
% Terminal-style extension:
%   drawOpts can disable rounded start/end caps and backshift the first
%   segment start point to create straight "wall-like" terminal edges.
%   It can also extend the end terminal forward and add a local start-wall
%   width boost (perpendicular to motion) to seal entrance artifacts
%   without shifting wall timing along the path tangent.
if isempty(coordsPx) || size(coordsPx, 1) ~= 2
    return;
end

if nargin < 5 || isempty(drawOpts)
    drawOpts = struct();
end
drawStartCap = true;
drawEndCap = true;
drawInternalJoinCaps = true;
startBackshiftPx = 0;
endForwardShiftPx = 0;
startWidthExtraPx = 0;
if isfield(drawOpts, 'drawStartCap') && ~isempty(drawOpts.drawStartCap)
    drawStartCap = logical(drawOpts.drawStartCap);
end
if isfield(drawOpts, 'startBackshiftPx') && ~isempty(drawOpts.startBackshiftPx)
    startBackshiftPx = double(drawOpts.startBackshiftPx);
end
if isfield(drawOpts, 'endForwardShiftPx') && ~isempty(drawOpts.endForwardShiftPx)
    endForwardShiftPx = double(drawOpts.endForwardShiftPx);
end
if isfield(drawOpts, 'startWidthExtraPx') && ~isempty(drawOpts.startWidthExtraPx)
    startWidthExtraPx = double(drawOpts.startWidthExtraPx);
end
if isfield(drawOpts, 'drawEndCap') && ~isempty(drawOpts.drawEndCap)
    drawEndCap = logical(drawOpts.drawEndCap);
end
if isfield(drawOpts, 'drawInternalJoinCaps') && ~isempty(drawOpts.drawInternalJoinCaps)
    drawInternalJoinCaps = logical(drawOpts.drawInternalJoinCaps);
end
if ~isfinite(startBackshiftPx) || startBackshiftPx < 0
    startBackshiftPx = 0;
end
if ~isfinite(endForwardShiftPx) || endForwardShiftPx < 0
    endForwardShiftPx = 0;
end
if ~isfinite(startWidthExtraPx) || startWidthExtraPx < 0
    startWidthExtraPx = 0;
end
coordsPx = local_apply_polyline_start_backshift(coordsPx, startBackshiftPx);
coordsPx = local_apply_polyline_end_forward_shift(coordsPx, endForwardShiftPx);

nCols = size(coordsPx, 2);
if nCols < 2
    return;
end
if mod(nCols, 2) ~= 0
    % Drop trailing orphan point defensively.
    nCols = nCols - 1;
end
if nCols < 2
    return;
end

lineWidthPx = double(lineWidthPx);
if ~isfinite(lineWidthPx) || lineWidthPx <= 0
    return;
end
lineWidthPx = max(1, lineWidthPx);
halfWidthPx = lineWidthPx / 2;

lineColor = double(lineColor(:)');
if isempty(lineColor)
    lineColor = [255 255 255];
elseif numel(lineColor) == 1
    lineColor = repmat(lineColor, 1, 3);
elseif numel(lineColor) > 3
    lineColor = lineColor(1:3);
end

% Straight-terminal fast path: draw one continuous filled band polygon
% instead of many short quads. This avoids seam artifacts ("sunburst")
% near exits/entrances while preserving wall-like terminals.
if ~drawStartCap && ~drawEndCap
    polylinePx = local_drawlines_coords_to_polyline(coordsPx, nCols);
    if local_fill_straight_band_polygon( ...
            window, polylinePx, halfWidthPx, lineColor, startWidthExtraPx)
        return;
    end
end

% For continuous polylines in paired format [p1 p2 p2 p3 ...], draw the
% first vertex cap once (optional), then:
%   - internal join caps (optional, recommended for continuity),
%   - final endpoint cap (optional for straight terminal edge).
% This keeps joins closed while allowing straight wall-like terminals.
needStartCap = drawStartCap;
nSegments = floor(nCols / 2);
segEndCumLenPx = zeros(nSegments, 1);
totalPolylineLenPx = 0;
for iSegLen = 1:nSegments
    iColLen = 2 * iSegLen - 1;
    x1Len = coordsPx(1, iColLen);
    y1Len = coordsPx(2, iColLen);
    x2Len = coordsPx(1, iColLen + 1);
    y2Len = coordsPx(2, iColLen + 1);
    segLenPx = 0;
    if isfinite(x1Len) && isfinite(y1Len) && isfinite(x2Len) && isfinite(y2Len)
        segLenPx = hypot(x2Len - x1Len, y2Len - y1Len);
        if ~isfinite(segLenPx) || segLenPx <= eps
            segLenPx = 0;
        end
    end
    totalPolylineLenPx = totalPolylineLenPx + segLenPx;
    segEndCumLenPx(iSegLen) = totalPolylineLenPx;
end

% When terminal caps are disabled (straight style), keep a flat terminal
% zone of one half-width at each end by suppressing round join discs there.
terminalFlatZonePx = halfWidthPx;

for iSeg = 1:nSegments
    iCol = 2 * iSeg - 1;
    x1 = coordsPx(1, iCol);
    y1 = coordsPx(2, iCol);
    x2 = coordsPx(1, iCol + 1);
    y2 = coordsPx(2, iCol + 1);
    if ~(isfinite(x1) && isfinite(y1) && isfinite(x2) && isfinite(y2))
        continue;
    end

    dx = x2 - x1;
    dy = y2 - y1;
    segLen = hypot(dx, dy);
    if segLen <= eps
        if needStartCap
            local_fill_disc(window, x1, y1, halfWidthPx, lineColor);
            needStartCap = false;
        end
        if iSeg < nSegments
            if drawInternalJoinCaps
                drawJoinCap = local_should_draw_join_cap( ...
                    iSeg, drawStartCap, drawEndCap, terminalFlatZonePx, ...
                    segEndCumLenPx, totalPolylineLenPx);
                if drawJoinCap
                    local_fill_disc(window, x2, y2, halfWidthPx, lineColor);
                end
            end
        elseif drawEndCap
            local_fill_disc(window, x2, y2, halfWidthPx, lineColor);
        end
        continue;
    end

    % Build one convex quad expanded perpendicular to segment direction.
    nx = -dy / segLen;
    ny = dx / segLen;
    ox = nx * halfWidthPx;
    oy = ny * halfWidthPx;
    quad = [x1 + ox, y1 + oy; ...
        x1 - ox, y1 - oy; ...
        x2 - ox, y2 - oy; ...
        x2 + ox, y2 + oy];
    Screen('FillPoly', window, lineColor, quad, 1);

    if needStartCap
        local_fill_disc(window, x1, y1, halfWidthPx, lineColor);
        needStartCap = false;
    end
    if iSeg < nSegments
        if drawInternalJoinCaps
            drawJoinCap = local_should_draw_join_cap( ...
                iSeg, drawStartCap, drawEndCap, terminalFlatZonePx, ...
                segEndCumLenPx, totalPolylineLenPx);
            if drawJoinCap
                local_fill_disc(window, x2, y2, halfWidthPx, lineColor);
            end
        end
    elseif drawEndCap
        local_fill_disc(window, x2, y2, halfWidthPx, lineColor);
    end
end
end

function polylinePx = local_drawlines_coords_to_polyline(coordsPx, nCols)
% LOCAL_DRAWLINES_COORDS_TO_POLYLINE Convert DrawLines pair format to Nx2 polyline.
polylinePx = zeros(0, 2);
if isempty(coordsPx) || size(coordsPx, 1) ~= 2
    return;
end
if nargin < 2 || isempty(nCols)
    nCols = size(coordsPx, 2);
end
if nCols < 2
    return;
end
if mod(nCols, 2) ~= 0
    nCols = nCols - 1;
end
if nCols < 2
    return;
end

nSeg = nCols / 2;
polylinePx = zeros(nSeg + 1, 2);
polylinePx(1, :) = [coordsPx(1, 1), coordsPx(2, 1)];
for iSeg = 1:nSeg
    iEnd = 2 * iSeg;
    polylinePx(iSeg + 1, :) = [coordsPx(1, iEnd), coordsPx(2, iEnd)];
end

finiteMask = isfinite(polylinePx(:, 1)) & isfinite(polylinePx(:, 2));
polylinePx = polylinePx(finiteMask, :);
if size(polylinePx, 1) < 2
    polylinePx = zeros(0, 2);
    return;
end

% Remove duplicate consecutive points to keep tangent estimates stable.
keepMask = true(size(polylinePx, 1), 1);
for iPt = 2:size(polylinePx, 1)
    if norm(polylinePx(iPt, :) - polylinePx(iPt - 1, :)) <= 1e-6
        keepMask(iPt) = false;
    end
end
polylinePx = polylinePx(keepMask, :);
if size(polylinePx, 1) < 2
    polylinePx = zeros(0, 2);
end
end

function drawn = local_fill_straight_band_polygon( ...
        window, polylinePx, halfWidthPx, lineColor, startWidthExtraPx)
% LOCAL_FILL_STRAIGHT_BAND_POLYGON Draw continuous wall-terminated band.
%
% Geometry:
%   Build left/right offsets around the center polyline and fill the
%   closed polygon [left(1..N), right(N..1)].
%
% Optional start widening:
%   startWidthExtraPx enlarges the wall thickness only near the first
%   terminal edge (perpendicular to motion), without shifting the wall
%   position along the path tangent.
drawn = false;
if isempty(polylinePx) || size(polylinePx, 1) < 2 || ~isfinite(halfWidthPx) || halfWidthPx <= 0
    return;
end
if nargin < 5 || isempty(startWidthExtraPx)
    startWidthExtraPx = 0;
end

nPts = size(polylinePx, 1);
tangent = zeros(nPts, 2);
for iPt = 1:nPts
    if iPt == 1
        v = polylinePx(2, :) - polylinePx(1, :);
    elseif iPt == nPts
        v = polylinePx(nPts, :) - polylinePx(nPts - 1, :);
    else
        vPrev = polylinePx(iPt, :) - polylinePx(iPt - 1, :);
        vNext = polylinePx(iPt + 1, :) - polylinePx(iPt, :);
        nPrev = norm(vPrev);
        nNext = norm(vNext);
        if nPrev <= 1e-9 && nNext <= 1e-9
            v = [0, 0];
        elseif nPrev <= 1e-9
            v = vNext;
        elseif nNext <= 1e-9
            v = vPrev;
        else
            uPrev = vPrev ./ nPrev;
            uNext = vNext ./ nNext;
            v = uPrev + uNext;
            if norm(v) <= 1e-9
                v = uNext;
            end
        end
    end

    vNorm = norm(v);
    if vNorm <= 1e-9
        if iPt > 1 && norm(tangent(iPt - 1, :)) > 0
            tangent(iPt, :) = tangent(iPt - 1, :);
        else
            tangent(iPt, :) = [1, 0];
        end
    else
        tangent(iPt, :) = v ./ vNorm;
    end
end

normal = [-tangent(:, 2), tangent(:, 1)];
leftEdge = polylinePx + halfWidthPx .* normal;
rightEdge = polylinePx - halfWidthPx .* normal;

polygonPx = [leftEdge; flipud(rightEdge)];
finiteMask = isfinite(polygonPx(:, 1)) & isfinite(polygonPx(:, 2));
polygonPx = polygonPx(finiteMask, :);
if size(polygonPx, 1) < 3
    return;
end

polyArea = polyarea(polygonPx(:, 1), polygonPx(:, 2));
if ~isfinite(polyArea) || polyArea <= 1e-3
    return;
end

Screen('FillPoly', window, lineColor, polygonPx, 1);
local_fill_straight_start_width_boost( ...
    window, polylinePx, halfWidthPx, startWidthExtraPx, lineColor);
drawn = true;
end

function local_fill_straight_start_width_boost( ...
        window, polylinePx, halfWidthPx, startWidthExtraPx, lineColor)
% LOCAL_FILL_STRAIGHT_START_WIDTH_BOOST Widen only the first terminal wall.
%
% Data flow:
%   first polyline tangent + requested width boost -> short tapered strip.
%   The strip starts exactly at the existing terminal wall location, so the
%   path-wise collision timing is unchanged while perpendicular coverage
%   near entrance is increased.
if ~isfinite(startWidthExtraPx) || startWidthExtraPx <= 0
    return;
end
if isempty(polylinePx) || size(polylinePx, 1) < 2
    return;
end

p1 = polylinePx(1, :);
tangent = [];
for iPt = 2:size(polylinePx, 1)
    d = polylinePx(iPt, :) - p1;
    dNorm = norm(d);
    if dNorm > 1e-9
        tangent = d ./ dNorm;
        break;
    end
end
if isempty(tangent)
    return;
end

normal = [-tangent(2), tangent(1)];
depthPx = max(1, halfWidthPx);
wideHalf = halfWidthPx + startWidthExtraPx;
p2 = p1 + depthPx .* tangent;

boostPoly = [p1 + wideHalf .* normal; ...
    p1 - wideHalf .* normal; ...
    p2 - halfWidthPx .* normal; ...
    p2 + halfWidthPx .* normal];
Screen('FillPoly', window, lineColor, boostPoly, 1);
end

function drawJoinCap = local_should_draw_join_cap( ...
        iSeg, drawStartCap, drawEndCap, terminalFlatZonePx, ...
        segEndCumLenPx, totalPolylineLenPx)
% LOCAL_SHOULD_DRAW_JOIN_CAP Keep straight terminal zones free of round joins.
drawJoinCap = true;
if ~isfinite(iSeg) || iSeg < 1 || iSeg > numel(segEndCumLenPx)
    return;
end

joinDistFromStartPx = segEndCumLenPx(iSeg);
if ~isfinite(joinDistFromStartPx)
    joinDistFromStartPx = inf;
end
if ~isfinite(totalPolylineLenPx)
    totalPolylineLenPx = joinDistFromStartPx;
end
joinDistFromEndPx = max(0, totalPolylineLenPx - joinDistFromStartPx);

terminalEpsPx = 1e-6;
if ~drawStartCap && (joinDistFromStartPx <= (terminalFlatZonePx + terminalEpsPx))
    drawJoinCap = false;
    return;
end
if ~drawEndCap && (joinDistFromEndPx <= (terminalFlatZonePx + terminalEpsPx))
    drawJoinCap = false;
end
end

function drawOpts = local_make_pathband_draw_options( ...
        Conf, isOccluderBand, lineWidthPx, dotRadiusPx, isPostBand)
% LOCAL_MAKE_PATHBAND_DRAW_OPTIONS Resolve terminal-style draw options.
%
% Data flow:
%   validated runtime Conf + trial geometry metadata -> drawing flags.
%   Geometry shifts are generated in stimuli_generation_V28_rescueTraject
%   and serialized into pathband_*_xy. Runtime still applies small terminal
%   compensations in straight mode to preserve frame-locked occlusion
%   timing and avoid entrance leakage artifacts.
drawOpts = struct();
drawOpts.drawStartCap = true;
drawOpts.drawEndCap = true;
drawOpts.drawInternalJoinCaps = true;
drawOpts.startBackshiftPx = 0;
drawOpts.endForwardShiftPx = 0;
drawOpts.startWidthExtraPx = 0;

if nargin < 3 || isempty(lineWidthPx)
    lineWidthPx = 1;
end
if nargin < 4 || isempty(dotRadiusPx)
    dotRadiusPx = 0;
end
if nargin < 5 || isempty(isPostBand)
    isPostBand = false;
end

if ~logical(isOccluderBand)
    return;
end

style = 'round';
if isfield(Conf, 'pathBandEntranceStyle') && ~isempty(Conf.pathBandEntranceStyle)
    style = lower(strtrim(char(Conf.pathBandEntranceStyle)));
end
switch style
    case 'round'
        return;
    case 'straight'
        drawOpts.drawStartCap = false;
        drawOpts.drawEndCap = false;
        drawOpts.drawInternalJoinCaps = true;

        halfLineWidthPx = 0.5 * max(1, double(lineWidthPx));
        dotRadiusPx = max(0, double(dotRadiusPx));

        backshiftScale = 1;
        if isfield(Conf, 'pathBandStraightBackshiftDotRadiusScale') && ...
                ~isempty(Conf.pathBandStraightBackshiftDotRadiusScale)
            backshiftScale = double(Conf.pathBandStraightBackshiftDotRadiusScale);
        end
        if ~isfinite(backshiftScale) || backshiftScale < 0
            backshiftScale = 1;
        end
        drawOpts.startBackshiftPx = halfLineWidthPx + backshiftScale * dotRadiusPx;

        % Pixel-domain correction for display-dependent quantization:
        % values < 0 slightly reduce the straight start-wall backshift to
        % avoid one-frame-early full occlusion on some ppd/rounding combos.
        if isfield(Conf, 'pathBandStraightStartBackshiftPixelOffset') && ...
                ~isempty(Conf.pathBandStraightStartBackshiftPixelOffset)
            drawOpts.startBackshiftPx = drawOpts.startBackshiftPx + ...
                double(Conf.pathBandStraightStartBackshiftPixelOffset);
        end
        if ~isfinite(drawOpts.startBackshiftPx)
            drawOpts.startBackshiftPx = 0;
        end
        drawOpts.startBackshiftPx = max(0, drawOpts.startBackshiftPx);

        endForwardScale = 1;
        if isfield(Conf, 'pathBandStraightEndForwardDotRadiusScale') && ...
                ~isempty(Conf.pathBandStraightEndForwardDotRadiusScale)
            endForwardScale = double(Conf.pathBandStraightEndForwardDotRadiusScale);
        end
        if ~isfinite(endForwardScale) || endForwardScale < 0
            endForwardScale = 1;
        end
        drawOpts.endForwardShiftPx = endForwardScale * dotRadiusPx;

        % Increase only perpendicular coverage at post-band entrance
        % without shifting wall position along the motion tangent.
        if logical(isPostBand)
            widthExtraScale = 0;
            if isfield(Conf, 'pathBandStraightPostStartWidthExtraDotRadiusScale') && ...
                    ~isempty(Conf.pathBandStraightPostStartWidthExtraDotRadiusScale)
                widthExtraScale = double(Conf.pathBandStraightPostStartWidthExtraDotRadiusScale);
            end
            if ~isfinite(widthExtraScale) || widthExtraScale < 0
                widthExtraScale = 0;
            end
            drawOpts.startWidthExtraPx = widthExtraScale * dotRadiusPx;
        end
    otherwise
        error('Unsupported pathBandEntranceStyle: %s', style);
end
end

function coordsPx = local_apply_polyline_start_backshift(coordsPx, backshiftPx)
% LOCAL_APPLY_POLYLINE_START_BACKSHIFT Extend first segment backwards by pixels.
if backshiftPx <= 0 || isempty(coordsPx) || size(coordsPx, 1) ~= 2 || size(coordsPx, 2) < 2
    return;
end

x1 = coordsPx(1, 1);
y1 = coordsPx(2, 1);
nCols = size(coordsPx, 2);
for iCol = 2:2:nCols
    x2 = coordsPx(1, iCol);
    y2 = coordsPx(2, iCol);
    if ~(isfinite(x1) && isfinite(y1) && isfinite(x2) && isfinite(y2))
        continue;
    end

    dx = x2 - x1;
    dy = y2 - y1;
    segLen = hypot(dx, dy);
    if segLen <= eps
        continue;
    end

    ux = dx / segLen;
    uy = dy / segLen;
    coordsPx(1, 1) = x1 - ux * backshiftPx;
    coordsPx(2, 1) = y1 - uy * backshiftPx;
    return;
end
end

function coordsPx = local_apply_polyline_end_forward_shift(coordsPx, forwardShiftPx)
% LOCAL_APPLY_POLYLINE_END_FORWARD_SHIFT Extend last segment forward by pixels.
if forwardShiftPx <= 0 || isempty(coordsPx) || size(coordsPx, 1) ~= 2 || size(coordsPx, 2) < 2
    return;
end

nCols = size(coordsPx, 2);
x2 = coordsPx(1, nCols);
y2 = coordsPx(2, nCols);
for iColStart = (nCols - 1):-2:1
    x1 = coordsPx(1, iColStart);
    y1 = coordsPx(2, iColStart);
    if ~(isfinite(x1) && isfinite(y1) && isfinite(x2) && isfinite(y2))
        continue;
    end

    dx = x2 - x1;
    dy = y2 - y1;
    segLen = hypot(dx, dy);
    if segLen <= eps
        continue;
    end

    ux = dx / segLen;
    uy = dy / segLen;
    coordsPx(1, nCols) = x2 + ux * forwardShiftPx;
    coordsPx(2, nCols) = y2 + uy * forwardShiftPx;
    return;
end
end

function local_fill_disc(window, xCenter, yCenter, radiusPx, color)
% LOCAL_FILL_DISC Draw a circular cap used to join adjacent thick segments.
if ~isfinite(xCenter) || ~isfinite(yCenter) || ~isfinite(radiusPx) || radiusPx <= 0
    return;
end
capRect = [xCenter - radiusPx, yCenter - radiusPx, ...
    xCenter + radiusPx, yCenter + radiusPx];
Screen('FillOval', window, color, capRect);
end

function [debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel] = ...
        local_emit_trigger(Conf, triggerValue, trialNum, frameNum, label, ...
        debugTriggerLog, debugTriggerTrial, debugTriggerFrame, debugTriggerSec, debugTriggerLabel)
% LOCAL_EMIT_TRIGGER Emit one trigger (MEG optional) and append debug logs.
if Conf.MEG
    Datapixx('StopDoutSchedule');
    pulse = [1 0] .* triggerValue;
    Datapixx('WriteDoutBuffer', pulse);
    Datapixx('SetDoutSchedule', 1.0 / Conf.refrate, 1000, 2);
    Datapixx('StartDoutSchedule');
    Datapixx('RegWr');
end
if isfield(Conf, 'useEyelink') && Conf.useEyelink && ...
        isfield(Conf, 'eyeTrackerSendTriggerMessages') && Conf.eyeTrackerSendTriggerMessages
    local_eyelink_safe_message(sprintf('TRIGGER %d trial %d frame %d label %s', ...
        round(double(triggerValue)), round(double(trialNum)), round(double(frameNum)), char(label)));
end
triggerTimeSec = GetSecs;

debugTriggerLog(end + 1, 1) = triggerValue;
debugTriggerTrial(end + 1, 1) = trialNum;
debugTriggerFrame(end + 1, 1) = frameNum;
debugTriggerSec(end + 1, 1) = triggerTimeSec;
debugTriggerLabel{end + 1, 1} = label;

if local_is_trigger_debug_enabled(Conf)
    fprintf('[trigger t=%.6f] value=%d trial=%d frame=%d label=%s\n', ...
        triggerTimeSec, triggerValue, trialNum, frameNum, label);
end
end

function local_emit_escape_termination_trigger(Conf, blockNum, label)
% LOCAL_EMIT_ESCAPE_TERMINATION_TRIGGER Emit the dedicated ESC-termination trigger.
if ~isfield(Conf, 'trigger') || ~isfield(Conf.trigger, 'escapeTermination')
    return;
end
[~, ~, ~, ~, ~] = local_emit_trigger(Conf, Conf.trigger.escapeTermination, ...
    round(double(blockNum)), 0, char(label), [], [], [], [], {});
end

function local_draw_catch_question_screen(window, Conf, xCenter, yCenter, fixSizePx)
% LOCAL_DRAW_CATCH_QUESTION_SCREEN Draw the catch-question layout.
%
% Layout:
%   - fixation at center
%   - configurable question text near top
%   - NO (left) and YES (right) near bottom
Screen('FillRect', window, Conf.background);
local_draw_fixation(window, xCenter, yCenter, fixSizePx, Conf.fixColor);

windowRect = Screen('Rect', window);
screenW = windowRect(3) - windowRect(1);
screenH = windowRect(4) - windowRect(2);

questionY = windowRect(2) + round(double(Conf.catchQuestionTopYFrac) * screenH);
bottomY = windowRect(2) + round(double(Conf.catchQuestionBottomYFrac) * screenH);
noX = windowRect(1) + round(double(Conf.catchBottomLeftXFrac) * screenW);
yesX = windowRect(1) + round(double(Conf.catchBottomRightXFrac) * screenW);

DrawFormattedText(window, Conf.catchPromptText, 'center', questionY, Conf.fixColor);
DrawFormattedText(window, Conf.catchNoLabel, noX, bottomY, Conf.fixColor);
DrawFormattedText(window, Conf.catchYesLabel, yesX, bottomY, Conf.fixColor);
end

function local_present_catch_feedback(window, Conf, xCenter, yCenter, fixSizePx, timedOut, responseCorrect)
% LOCAL_PRESENT_CATCH_FEEDBACK Show post-question feedback for catch trials.
Screen('FillRect', window, Conf.background);

if timedOut
    % Timeout feedback is text-only (no fixation).
    DrawFormattedText(window, Conf.catchTimeoutText, 'center', 'center', Conf.fixColor);
else
    if isfinite(responseCorrect)
        if responseCorrect
            feedbackColor = Conf.catchFeedbackCorrectColor;
        else
            feedbackColor = Conf.catchFeedbackIncorrectColor;
        end
    else
        feedbackColor = Conf.fixColor;
    end
    local_draw_fixation(window, xCenter, yCenter, fixSizePx, feedbackColor);
end

Screen('Flip', window);
WaitSecs(max(0, double(Conf.catchFeedbackDurationSec)));
end

function enabled = local_is_trigger_debug_enabled(Conf)
% LOCAL_IS_TRIGGER_DEBUG_ENABLED Resolve trigger-debug print mode.
enabled = false;
if isfield(Conf, 'debugTriggerMode')
    enabled = (double(Conf.debugTriggerMode) ~= 0);
elseif isfield(Conf, 'debug')
    enabled = (double(Conf.debug) ~= 0);
end
end

function local_report_run_frame_timing(output, blockNum, runNum, actualFrameRate)
% LOCAL_REPORT_RUN_FRAME_TIMING Print measured frame pacing statistics.
%
% Purpose:
%   Provide a compact timing report for diagnosing apparent lag on displays
%   that differ from the stimulus-generation refresh rate.
durCells = {output.frame_duration_ms};
if isempty(durCells)
    return;
end
nonEmptyMask = ~cellfun(@isempty, durCells);
if ~any(nonEmptyMask)
    return;
end

allDurMs = vertcat(durCells{nonEmptyMask});
allDurMs = allDurMs(isfinite(allDurMs) & allDurMs > 0);
if isempty(allDurMs)
    return;
end

if actualFrameRate > 0
    nominalMs = 1000 / actualFrameRate;
else
    nominalMs = NaN;
end

medianMs = median(allDurMs);
maxMs = max(allDurMs);
jitterMs = std(allDurMs);
sortedDur = sort(allDurMs);
p95Idx = max(1, min(numel(sortedDur), ceil(0.95 * numel(sortedDur))));
p95Ms = sortedDur(p95Idx);

fprintf(['Timing summary block %d run %d: nominal %.2f ms, median %.2f ms, ' ...
    'p95 %.2f ms, max %.2f ms, std %.2f ms.\n'], ...
    blockNum, runNum, nominalMs, medianMs, p95Ms, maxMs, jitterMs);
end

function action = local_show_start_message_gate(window, Conf)
% LOCAL_SHOW_START_MESSAGE_GATE Show the configurable start screen and gate run start.
%
% Behavior:
%   - If disabled, continue immediately.
%   - If duration <= 0, wait indefinitely for keyboard/button-box 1/8.
%   - If duration > 0, continue on key press or timeout.
action = 'continue';
if ~isfield(Conf, 'startMessageEnabled') || ~logical(Conf.startMessageEnabled)
    return;
end

waitDurationSec = double(Conf.startMessageDurationSec);
if isfinite(waitDurationSec) && waitDurationSec > 0
    deadline = GetSecs + waitDurationSec;
else
    deadline = inf;
end
acceptAfter = GetSecs + max(0, double(Conf.messageResponseLockoutSec));

% Data flow: start-message state -> MEG log + keyboard polling -> gate action.
local_start_meg_response_log(Conf);
cleanupObj = onCleanup(@() local_stop_meg_response_log(Conf));

while true
    local_draw_text_message_screen(window, Conf, Conf.startMessageText);
    Screen('Flip', window);

    [accepted, escapePressed] = local_poll_start_message_response(Conf, acceptAfter);
    if escapePressed
        action = 'abort';
        return;
    end
    if accepted
        return;
    end
    if GetSecs >= deadline
        return;
    end
end
end

function [accepted, escapePressed] = local_poll_start_message_response(Conf, acceptAfter)
% LOCAL_POLL_START_MESSAGE_RESPONSE Poll one response for the start gate.
%
% Response mapping:
%   - keyboard 1/8
%   - MEG button-box value 1/8
accepted = false;
escapePressed = false;
canAccept = GetSecs >= acceptAfter;

if Conf.enableKeyboardDebugResponse
    [keyDown, ~, keyCode] = KbCheck;
    if keyDown
        if keyCode(Conf.escapeKeyCode)
            escapePressed = true;
            KbReleaseWait;
            return;
        end
        if canAccept && any(keyCode(Conf.startMessageAcceptKeyCodes))
            accepted = true;
            KbReleaseWait;
            return;
        end
    end
end

if Conf.MEG
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    if status.newLogFrames > 0 && mod(status.currentWriteFrame, 2)
        [data, ~] = Datapixx('ReadDinLog');
        if ~isempty(data)
            rawValue = bitand(uint32(data(end)), Conf.megButtonMask);
            if canAccept && any(rawValue == Conf.startMessageAcceptMegValues)
                accepted = true;
                return;
            end
        end
    end
end
end

function abortRequested = local_show_timed_transition_message(window, Conf, messageText, durationSec)
% LOCAL_SHOW_TIMED_TRANSITION_MESSAGE Show an informational message for a fixed duration.
%
% Behavior:
%   - No response is required.
%   - ESC can still abort the session for safety.
abortRequested = false;
showDurationSec = max(0, double(durationSec));
if showDurationSec == 0
    return;
end

deadline = GetSecs + showDurationSec;
while GetSecs < deadline
    messageToDraw = char(messageText);
    if isfield(Conf, 'run1TransitionCountdownEnabled') && logical(Conf.run1TransitionCountdownEnabled)
        remainingSec = max(0, ceil(deadline - GetSecs));
        messageToDraw = sprintf('%s\n\nSecond task starts in: %d s', messageToDraw, remainingSec);
    end
    local_draw_text_message_screen(window, Conf, messageToDraw);
    Screen('Flip', window);

    [keyDown, ~, keyCode] = KbCheck;
    if keyDown && keyCode(Conf.escapeKeyCode)
        abortRequested = true;
        KbReleaseWait;
        return;
    end
end
end

function abortRequested = local_show_post_message_fixation(window, Conf, xCenter, yCenter, fixSizePx)
% LOCAL_SHOW_POST_MESSAGE_FIXATION Show fixation-only screen after message screens.
abortRequested = false;
if ~isfield(Conf, 'postMessageFixationEnabled') || ~logical(Conf.postMessageFixationEnabled)
    return;
end

showDurationSec = max(0, double(Conf.postMessageFixationDurationSec));
if showDurationSec == 0
    return;
end

deadline = GetSecs + showDurationSec;
while GetSecs < deadline
    Screen('FillRect', window, Conf.background);
    local_draw_fixation(window, xCenter, yCenter, fixSizePx, Conf.fixColor);
    Screen('Flip', window);

    [keyDown, ~, keyCode] = KbCheck;
    if keyDown && keyCode(Conf.escapeKeyCode)
        abortRequested = true;
        KbReleaseWait;
        return;
    end
end
end

function local_draw_text_message_screen(window, Conf, messageText)
% LOCAL_DRAW_TEXT_MESSAGE_SCREEN Draw a centered full-screen instruction message.
Screen('FillRect', window, Conf.background);
DrawFormattedText(window, char(messageText), 'center', 'center', Conf.fixColor);
end

function action = local_wait_end_of_block_action(window, Conf, currentBlock, remainingBlocks)
% LOCAL_WAIT_END_OF_BLOCK_ACTION Show block-end message and wait for 1/8.
action = 'continue';
if ~isfield(Conf, 'endOfBlockMessageEnabled') || ~logical(Conf.endOfBlockMessageEnabled)
    return;
end

% Resolve block-end text from config format string with a safety fallback.
msgTemplate = char(Conf.endOfBlockMessageTextTemplate);
fallbackTemplate = ['Block %d ended.\n' ...
    'Remaining blocks: %d.\n\n' ...
    'Press Red or Blue to start the next block.'];
try
    if isempty(strfind(msgTemplate, '%')) %#ok<STREMP>
        % Template was likely pre-formatted incorrectly; restore full message.
        msg = sprintf(fallbackTemplate, currentBlock, remainingBlocks);
    else
        msg = sprintf(msgTemplate, currentBlock, remainingBlocks);
    end
catch
    % Keep runtime resilient to malformed custom templates.
    msg = sprintf(fallbackTemplate, currentBlock, remainingBlocks);
end

waitDurationSec = double(Conf.endOfBlockMessageDurationSec);
if isfinite(waitDurationSec) && waitDurationSec > 0
    deadline = GetSecs + waitDurationSec;
else
    deadline = inf;
end
acceptAfter = GetSecs + max(0, double(Conf.messageResponseLockoutSec));

local_start_meg_response_log(Conf);
cleanupObj = onCleanup(@() local_stop_meg_response_log(Conf));

while true
    local_draw_text_message_screen(window, Conf, msg);
    Screen('Flip', window);

    [accepted, escapePressed] = local_poll_start_message_response(Conf, acceptAfter);
    if escapePressed
        action = 'abort';
        return;
    end
    if accepted
        return;
    end
    if GetSecs >= deadline
        return;
    end
end
end

function action = local_wait_between_block_calibration_action(window, Conf)
% LOCAL_WAIT_BETWEEN_BLOCK_CALIBRATION_ACTION Ask whether to recalibrate before next block.
%
% Response mapping:
%   - skip calibration: red / key 1
%   - run calibration now: blue / key 8
%   - ESC: abort session
action = 'skip';
if ~isfield(Conf, 'betweenBlockCalibrationChoiceEnabled') || ...
        ~logical(Conf.betweenBlockCalibrationChoiceEnabled)
    return;
end
if ~isfield(Conf, 'eyeTrackerAllowBetweenBlockCalibration') || ...
        ~logical(Conf.eyeTrackerAllowBetweenBlockCalibration) || ...
        ~isfield(Conf, 'useEyelinkSession') || ~logical(Conf.useEyelinkSession)
    return;
end

waitDurationSec = double(Conf.betweenBlockCalibrationChoiceDurationSec);
if isfinite(waitDurationSec) && waitDurationSec > 0
    deadline = GetSecs + waitDurationSec;
else
    deadline = inf;
end
acceptAfter = GetSecs + max(0, double(Conf.messageResponseLockoutSec));

local_start_meg_response_log(Conf);
cleanupObj = onCleanup(@() local_stop_meg_response_log(Conf)); %#ok<NASGU>

while true
    local_draw_text_message_screen(window, Conf, Conf.betweenBlockCalibrationChoiceText);
    Screen('Flip', window);

    canAccept = GetSecs >= acceptAfter;
    if Conf.enableKeyboardDebugResponse
        [keyDown, ~, keyCode] = KbCheck;
        if keyDown
            if keyCode(Conf.escapeKeyCode)
                action = 'abort';
                KbReleaseWait;
                return;
            end
            if canAccept && any(keyCode(Conf.calibrationChoiceSkipKeyCodes))
                action = 'skip';
                KbReleaseWait;
                return;
            end
            if canAccept && any(keyCode(Conf.calibrationChoiceRunKeyCodes))
                action = 'calibrate';
                KbReleaseWait;
                return;
            end
        end
    end

    if Conf.MEG
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        if status.newLogFrames > 0 && mod(status.currentWriteFrame, 2)
            [data, ~] = Datapixx('ReadDinLog');
            if ~isempty(data)
                rawValue = bitand(uint32(data(end)), Conf.megButtonMask);
                if canAccept && any(rawValue == Conf.calibrationChoiceSkipMegValues)
                    action = 'skip';
                    return;
                end
                if canAccept && any(rawValue == Conf.calibrationChoiceRunMegValues)
                    action = 'calibrate';
                    return;
                end
            end
        end
    end

    if GetSecs >= deadline
        action = 'skip';
        return;
    end
end
end

function Conf = local_eyelink_run_calibration(Conf)
% LOCAL_EYELINK_RUN_CALIBRATION Run EyeLink calibration on participant request.
if ~isfield(Conf, 'useEyelinkSession') || ~Conf.useEyelinkSession || ...
        ~isfield(Conf, 'eyeTrackerInitialized') || ~Conf.eyeTrackerInitialized
    return;
end
if ~isfield(Conf, 'eyeTrackerEyelinkDefaults') || isempty(Conf.eyeTrackerEyelinkDefaults)
    return;
end
try
    EyelinkDoTrackerSetup(Conf.eyeTrackerEyelinkDefaults);
catch
    % Keep runtime robust if calibration is unavailable on current host.
end
WaitSecs(max(0, double(Conf.eyeTrackerWaitSec)));
end

function action = local_show_final_message(window, Conf)
% LOCAL_SHOW_FINAL_MESSAGE Show the final completion message.
action = 'continue';
if ~isfield(Conf, 'finalMessageEnabled') || ~logical(Conf.finalMessageEnabled)
    return;
end

msg = char(Conf.finalMessageText);
showDurationSec = max(0, double(Conf.finalMessageDurationSec));
if showDurationSec == 0
    local_draw_text_message_screen(window, Conf, msg);
    Screen('Flip', window);
    return;
end

deadline = GetSecs + showDurationSec;
while GetSecs < deadline
    local_draw_text_message_screen(window, Conf, msg);
    Screen('Flip', window);

    [keyDown, ~, keyCode] = KbCheck;
    if keyDown && keyCode(Conf.escapeKeyCode)
        action = 'abort';
        KbReleaseWait;
        return;
    end
end
end

function action = local_resolve_dryrun_break_action(simKeys, blockStepIdx)
% LOCAL_RESOLVE_DRYRUN_BREAK_ACTION Resolve simulated SPACE/Q action in dry-run tests.
action = 'continue';
if isempty(simKeys)
    return;
end

if ischar(simKeys) || isstring(simKeys)
    key = lower(char(simKeys));
elseif iscell(simKeys) && blockStepIdx <= numel(simKeys)
    key = lower(char(simKeys{blockStepIdx}));
else
    key = 'space';
end

if strcmp(key, 'q') || strcmp(key, 'quit') || strcmp(key, 'escape')
    action = 'quit';
end
end

function local_cleanup_runtime(Conf)
% LOCAL_CLEANUP_RUNTIME Best-effort cleanup for Psychtoolbox, EyeLink, DataPixx.
% EyeLink teardown runs for both:
%  - real EyeLink acquisition mode, and
%  - fake-gaze mode with explicit EyeLink calibration setup.
eyelinkSessionActive = false;
if isfield(Conf, 'useEyelinkSession')
    eyelinkSessionActive = logical(Conf.useEyelinkSession);
elseif isfield(Conf, 'useEyelink')
    eyelinkSessionActive = logical(Conf.useEyelink);
end
if eyelinkSessionActive
    if isfield(Conf, 'useEyelink') && Conf.useEyelink
        try
            Conf = local_eyelink_finalize_block_file(Conf);
        catch
        end
    end
    try
        Eyelink('StopRecording');
    catch
    end
    try
        Eyelink('Command', 'set_idle_mode');
    catch
    end
    try
        Eyelink('CloseFile');
    catch
    end
    try
        Eyelink('ShutDown');
    catch
    end
end
if isfield(Conf, 'MEG') && Conf.MEG
    try
        Datapixx('StopDinLog');
        Datapixx('SetDoutValues', 0);
        Datapixx('RegWrRd');
        Datapixx('Close');
    catch
        % Keep cleanup best-effort; do not raise nested cleanup errors.
    end
end

try
    Priority(0);
catch
end
try
    ShowCursor;
catch
end
try
    sca;
catch
end
end
