% Config_runtime_v17_blockResume.m
%
% Purpose:
%   Define workspace variables consumed by:
%     MoveDot1_experiment_occlusion_v17_blockResume.m
%   so the experiment can be started non-interactively with a reproducible
%   runtime preset (subject/block/mode/debug/eye-tracker/occluder settings).
%
% How it works:
%   - This file is a SCRIPT (not a class and not a function).
%   - Running it populates variables in the caller workspace.
%   - The experiment script keeps these variables via its
%     clearvars -except list and uses them as runtime overrides.
%
% Usage example (inside experiment/):
%   run('lib/Config_runtime_v17_blockResume.m');
%   run('MoveDot1_experiment_occlusion_v17_blockResume.m');
%
% Usage example (MATLAB batch from repo root):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    run('lib/Config_runtime_v17_blockResume.m'); ...
%    run('MoveDot1_experiment_occlusion_v17_blockResume.m');"
%
% Notes:
%   - Set optional fields to [] to use schedule defaults from:
%     Config_schedule_CreateInputV20_MoveDotV17_blockResume.
%   - This preset does NOT set Conf.MEG directly because the current
%     experiment script does not expose a workspace override for Conf.MEG.
%
% Inputs:
%   None (manual edits to values below).
%
% Outputs:
%   Variables in workspace:
%     iSub, iBlock, viewingDistanceMm, practice_mode, ...
%     ITIRangeSec, and all optional override variables supported by v17
%     block-resume.

%% Core session identity and run mode
% Data flow: these values define participant/session selection and suppress
% interactive dialogs inside the experiment script.
iSub = 50; % Accepted: finite integer subject ID. Default when unset in workspace: input dialog opens (dialog default value is 70). Example: 202.
iBlock = 1; % Accepted: finite integer block index in [1..numBlocks]. Default when unset: input dialog opens (dialog default value is 1).
viewingDistanceMm = 1000; % Accepted: finite scalar > 0 (mm). Default when unset: input dialog opens (dialog default value is 1000).
practice_mode = 0; % Accepted: {0,1}. Default when unset/empty: 0 (full mode). 1 enables practice-by-run trimming.

%% Execution mode controls
% Data flow: dry-run controls bypass PTB drawing and validate schedule logic.
dryRunValidateScheduleOnly = false; % Accepted: logical scalar. Default when unset: false. true -> no PTB window/triggers; schedule validation only.
debug_runNumber = []; % Accepted: [] or integer run index in [1..3]. Default when unset: []. [] runs all runs in each block.
dryRunSimulatedBreakKeys = {}; % Accepted: cell array of strings (dry-run only). Default when unset/empty: {}. Example: {'space','q'}.

%% Catch-question and debug display controls
% Data flow: these fields override default prompt/debug behavior in runtime.
question_text = 'Did the trajectory change while occluded?'; % Accepted: char/string scalar. Default when unset/empty: schedule value ('Has the dot changed its course?').
debugTriggerMode = 0; % Accepted: {0,1} (or logical). Default when unset: 0. 1 prints trigger emits in console.
runColorCueEnabled = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.runColorCueEnabled (current default: true).
practiceRunPercentageTrials = []; % Accepted: [] or finite scalar in (0,100]. Default when []/unset: scheduleCfg.practiceRunPercentageTrials (current default: 15).
debugShowFrameInfoOverlay = false; % Accepted: logical scalar. Default when unset: false.
debugShowFullPathOverlay = false; % Accepted: logical scalar. Default when unset: false.
run1TransitionCountdownEnabled = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.run1TransitionCountdownEnabled (current default: true).

%% Post-trial ITI jitter controls
% Data flow: ITIRangeSec override -> v17 runtime ITI precompute at run start
% -> trial-specific ITI hold (replay-safe via source-trial mapping).
ITIRangeSec = [0.500, 1.000]; % Accepted: [] or finite 1x2 [min max] sec with min<=max and both >=0. Default when unset/[] in runtime script: [0.500, 1.000].

%% Path-band occluder geometry controls
% Data flow: these values tune straight-terminal entrance/exit behavior.
% Keep [] to use defaults defined in the experiment script.
pathBandEntranceStyle = 'straight'; % Accepted: 'straight' or 'round'. Default when unset/empty: scheduleCfg.pathBandEntranceStyle (current default: 'straight').
pathBandStraightBackshiftDotRadiusScale = []; % Accepted: [] or finite scalar >= 0. Default when []/unset: scheduleCfg.pathBandStraightBackshiftDotRadiusScale (current default: 0.0).
pathBandStraightEndForwardDotRadiusScale = []; % Accepted: [] or finite scalar >= 0. Default when []/unset: scheduleCfg.pathBandStraightEndForwardDotRadiusScale (current default: 1.0).
pathBandStraightPostStartWidthExtraDotRadiusScale = []; % Accepted: [] or finite scalar >= 0. Default when []/unset: scheduleCfg.pathBandStraightPostStartWidthExtraDotRadiusScale (current default: 0.15).
pathBandStraightStartBackshiftPixelOffset = []; % Accepted: [] or finite scalar (pixels, signed). Default when []/unset: scheduleCfg.pathBandStraightStartBackshiftPixelOffset (current default: -0.5).

%% PTB sync policy override and flip-miss guardrails
% Data flow: optional timing guard overrides; [] keeps schedule defaults.
debugSkipSyncTests = []; % Accepted: [] or {0,1}/logical. Default when []/unset: policy mode (MEG->schedule skipSyncTestsWhenMEG=false; debug non-MEG may use skipSyncTestsInDebug=true).
flipMissedWarnThresholdCount = []; % Accepted: [] or integer >= 1. Default when []/unset: scheduleCfg.flipMissedWarnThresholdCount (current default: 1).
flipMissedAbortEnabled = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.flipMissedAbortEnabled (current default: false).
flipMissedAbortThresholdCount = []; % Accepted: [] or integer >= 1. Default when []/unset: scheduleCfg.flipMissedAbortThresholdCount (current default: 5).

%% Eye-tracker runtime controls
% Data flow: overrides for EyeLink/fake/disabled gaze monitoring behavior.
trackeye = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.eyeTrackerEnabled (current default: false).
ignoreEyeTracker = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.ignoreEyeTracker (current default: false).
fakeEyeTracker = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.fakeEyeTracker (current default: false).
fakeGazeBreakRate = []; % Accepted: [] or finite scalar in [0,1]. Default when []/unset: scheduleCfg.fakeGazeBreakRate (current default: 0.05).
enableFixationAbort = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.enableFixationAbort (current default: true).
fixWindowDeg = []; % Accepted: [] or finite scalar > 0 (deg). Default when []/unset: scheduleCfg.fixWindowDeg (current default: 7.0).
fixBreakToleranceFrames = []; % Accepted: [] or integer >= 1. Default when []/unset: scheduleCfg.fixBreakToleranceFrames (current default: 5).
fixationWarningDurationSec = []; % Accepted: [] or finite scalar >= 0 (sec). Default when []/unset: scheduleCfg.fixationWarningDurationSec (current default: 2.0).
replayLimit = []; % Accepted: [] or integer >= 0. Default when []/unset: scheduleCfg.replayLimit (current default: 1).
allowInfiniteReplays = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.allowInfiniteReplays (current default: false).
eyeTrackerDoCalibration = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.eyeTrackerDoCalibration (current default: true).
eyeTrackerAllowBetweenBlockCalibration = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.eyeTrackerAllowBetweenBlockCalibration (current default: true).
eyeTrackerSendTriggerMessages = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.eyeTrackerSendTriggerMessages (current default: true).
eyeTrackerImageTransferEnabled = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.eyeTrackerImageTransferEnabled (current default: true).
                                    % Behavior details: if unsupported or file missing, code fails softly (try/catch), so recording continues.
                                    % Practical use: set true when you want host-side visual context; set false to minimize host-command overhead.
                                    % Example: eyeTrackerImageTransferEnabled = false;  % lean runtime in MEG sessions.

%% Example presets (optional)
% Data flow: quick copy/paste profiles for common scenarios.
% Example A: hardware run with EyeLink and strict timing checks
% trackeye = true; ignoreEyeTracker = false; fakeEyeTracker = false;
% debugSkipSyncTests = 0; flipMissedAbortEnabled = true; flipMissedAbortThresholdCount = 3;
%
% Example B: development run on laptop without EyeLink
% trackeye = false; ignoreEyeTracker = true; fakeEyeTracker = false;
% debugShowFrameInfoOverlay = true; debugShowFullPathOverlay = true; debugSkipSyncTests = 1;

%% Preset summary print
% Data flow: lightweight audit line confirming active runtime preset.
itiRangeForPrint = ITIRangeSec;
if isempty(itiRangeForPrint)
    itiRangeForPrint = [0.500, 1.000];
end
fprintf(['Loaded runtime preset for MoveDot1 v17 eyeTrackerReplay blockResume: ' ...
    'SUB=%02d startBlock=%02d practice=%d dryRun=%d debugRun=%s ITI=[%.3f %.3f] entranceStyle=%s\n'], ...
    iSub, iBlock, practice_mode, double(dryRunValidateScheduleOnly), ...
    mat2str(debug_runNumber), itiRangeForPrint(1), itiRangeForPrint(2), char(pathBandEntranceStyle));
