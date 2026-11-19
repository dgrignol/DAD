% Config_runtime_v24_allTrialGuard.m
%
% Purpose:
%   Define workspace variables consumed by:
%     MoveDot1_experiment_occlusion_v24_allTrialGuard.m
%   so the experiment can be started non-interactively with a reproducible
%   runtime preset (subject/block/mode/debug/eye-tracker/occluder,
%   scheduled-onset, display-validation, and photodiode settings).
%
% How it works:
%   - This file is a SCRIPT (not a class and not a function).
%   - Running it populates variables in the caller workspace.
%   - The experiment script keeps these variables via its
%     clearvars -except list and uses them as runtime overrides.
%
% Usage example (inside experiment/):
%   run('lib/Config_runtime_v24_allTrialGuard.m');
%   run('MoveDot1_experiment_occlusion_v24_allTrialGuard.m');
%
% Usage example (MATLAB batch from repo root):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    run('lib/Config_runtime_v24_allTrialGuard.m'); ...
%    run('MoveDot1_experiment_occlusion_v24_allTrialGuard.m');"
%
% Notes:
%   - Set optional fields to [] to use schedule defaults from:
%     Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.
%   - MEG mode is controlled by workspace variable MEG (0/1).
%   - This preset measures and validates the active display mode; it never
%     configures the operating system or PROPixx refresh mode.
%
% Inputs:
%   None (manual edits to values below).
%
% Outputs:
%   Variables in workspace:
%     iSub, iBlock, viewingDistanceMm, practice_mode, MEG, ...
%     ITIRangeSec, photodiodePatchEnabled, photodiodePatchSizePx,
%     photodiodePatchPositionPx, photodiode colors/pulse widths, and all
%     scheduled-onset, all-trial guard, PTB sync, and warm-up controls used by v24.
%   - The invisible renderer preflight remains run-boundary-only, but the
%     3 ms lead check applies to every retained 71/81/91 onset. Up to six
%     low-patch rollover frames are allowed; exhaustion is saved and the trial
%     continues so this safeguard never terminates participant acquisition.

%% Core session identity and run mode
% Data flow: these values define participant/session selection and suppress
% interactive dialogs inside the experiment script.
iSub = 50; % Accepted: finite integer subject ID. Default when unset in workspace: input dialog opens (dialog default value is 70). Example: 202.
iBlock = 1; % Accepted: finite integer block index in [1..numBlocks]. Default when unset: input dialog opens (dialog default value is 1).
viewingDistanceMm = 1000; % Accepted: finite scalar > 0 (mm). Default when unset: input dialog opens (dialog default value is 1000).
practice_mode = 0; % Accepted: {0,1}. Default when unset/empty: 0 (full mode). 1 enables practice-by-run trimming.
MEG = true; % Accepted: [] or logical/{0,1}. Default when unset/[] in runtime script: 0 (disabled). true enables DataPixx/MEG trigger path.

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
% Data flow: ITIRangeSec override -> v23 runtime ITI precompute at run start
% -> trial-specific ITI hold (replay-safe via source-trial mapping).
ITIRangeSec = [1.000, 1.500]; % Accepted: [] or finite 1x2 [min max] sec with min<=max and both >=0. Default when unset/[] in runtime script: [0.500, 1.000].

%% ITI visual separation style
% Data flow: runtime ITI visual-mode overrides -> post-trial ITI renderer.
% Supported:
%   - 'screen_gray': full screen AND arena rectangle turn gray during ITI
%     only during the configured ITI visual-separation window.
%   - 'fixation_gray': only fixation turns gray during the configured
%     ITI visual-separation window (background stays black).
%   - 'fixation_short': fixation-only cue during ITI with start delay and
%     fixed cue duration (no tail guard and no full-screen color change).
%   - 'off': no ITI visual-separation cue.
itiVisualMode = 'fixation_short'; % Accepted: [] or {'screen_gray','fixation_gray','fixation_short','off'}. Default when unset/[] in runtime script: 'screen_gray'.
itiGreyColor = [90, 90, 90]; % Accepted: [] or finite RGB triplet in [0,255]. Default when unset/[]: [90 90 90]. Also used as fixation cue color in 'fixation_short'.
itiVisualSeparationDur = [0.300, 0.700]; % Accepted: [] or finite 1x2 [startSec endSec] with 0<=start<=end. Default when unset/[]: [0.300, 0.700].
itiFixationShortDelaySec = 0.100; % Accepted: [] or finite scalar >= 0 (sec). Used only when itiVisualMode='fixation_short'. Default when unset/[]: 0.100.
itiFixationShortDurationSec = 0.300; % Accepted: [] or finite scalar >= 0 (sec). Used only when itiVisualMode='fixation_short'. Default when unset/[]: 0.300.

%% Catch-question layout controls
% Data flow: catch-question layout fractions -> question/YES/NO draw positions.
catchQuestionTopYFrac = 0.42; % Accepted: [] or finite scalar in (0,1). Vertical fraction for question text.
catchQuestionBottomYFrac = 0.58; % Accepted: [] or finite scalar in (0,1). Vertical fraction for YES/NO labels.
catchBottomLeftXFrac = 0.38; % Accepted: [] or finite scalar in (0,1). Horizontal fraction for NO label.
catchBottomRightXFrac = 0.62; % Accepted: [] or finite scalar in (0,1). Horizontal fraction for YES label.

%% Mini blink-break controls
% Data flow: block trial counter -> optional timed blink-break message.
miniBlinkBreakEnabled = true; % Accepted: [] or logical. true enables automatic blink-break screens.
miniBlinkBreakEveryNTrials = 10; % Accepted: [] or integer >= 1. Candidate break interval in executed trials.
miniBlinkBreakDurationSec = 5.0; % Accepted: [] or finite scalar >= 0. Countdown duration in seconds.

%% Path-band occluder geometry controls
% Data flow: these values tune straight-terminal entrance/exit behavior.
% Keep [] to use defaults defined in the experiment script.
pathBandEntranceStyle = 'straight'; % Accepted: 'straight' or 'round'. Default when unset/empty: scheduleCfg.pathBandEntranceStyle (current default: 'straight').
pathBandStraightBackshiftDotRadiusScale = []; % Accepted: [] or finite scalar >= 0. Default when []/unset: scheduleCfg.pathBandStraightBackshiftDotRadiusScale (current default: 0.0).
pathBandStraightEndForwardDotRadiusScale = []; % Accepted: [] or finite scalar >= 0. Default when []/unset: scheduleCfg.pathBandStraightEndForwardDotRadiusScale (current default: 1.0).
pathBandStraightPostStartWidthExtraDotRadiusScale = []; % Accepted: [] or finite scalar >= 0. Default when []/unset: scheduleCfg.pathBandStraightPostStartWidthExtraDotRadiusScale (current default: 0.15).
pathBandStraightStartBackshiftPixelOffset = []; % Accepted: [] or finite scalar (pixels, signed). Default when []/unset: scheduleCfg.pathBandStraightStartBackshiftPixelOffset (current default: -0.5).

%% Photodiode patch and pulse controls
% Data flow: values below -> v23 runtime validation -> absolute patch
% rectangle after the PTB window opens -> final draw call before every flip.
% Position is an [x y] pixel offset from the window's top-left corner.
% Condition triggers 71/81/91 use the trial pulse width; visual event
% triggers 111/112/114/115/116/117/118 use the event pulse width.
photodiodePatchEnabled = []; % Accepted: [] or logical. Default when []/unset: follows MEG (enabled in MEG mode, disabled otherwise).
photodiodePatchSizePx = [50 50]; % Accepted: finite [width height] values >= 1 pixel.
photodiodePatchPositionPx = [0 0]; % Accepted: finite non-negative [x y] pixel offset from window top-left.
photodiodePatchHighColor = [255 255 255]; % Accepted: finite RGB triplet in [0,255].
photodiodePatchLowColor = [0 0 0]; % Accepted: finite RGB triplet in [0,255], different from high color.
photodiodeTrialPulseFrames = 4; % Accepted: integer >= 1; used for condition/trial-onset markers.
photodiodeEventPulseFrames = 2; % Accepted: integer >= 1; used for marked within-trial/question events.

%% Scheduled trial-onset and display-validation controls
% Data flow:
%   measured PTB IFI + existing randomized ITI -> continuous scheduled
%   fixation frames -> predetermined first-stimulus target frame.
%
% The warm-up runs once per MATLAB experiment launch, immediately before
% the first retained trial marker. With the defaults it uses ten cycles of
% four high-patch frames plus one explicit low frame. At 120 Hz this takes
% about 0.42 s, excluding the single non-retained anchor flip.
scheduledOnsetWarmupEnabled = []; % Accepted: [] or logical. [] follows MEG mode; true enables one session-level warm-up.
scheduledOnsetWarmupCycles = 10; % Accepted: integer >= 0. Ten cycles are approximately 0.42 s at 120 Hz.
scheduledOnsetWarmupTriggerCode = 153; % Accepted: integer in [1,255], collision-checked against fixed and actual sequence codes.
scheduledOnsetReanchorFrames = 10; % Accepted: integer >= 1. Silent baseline frames used only when no randomized ITI precedes an onset.
scheduledOnsetDeadlineFraction = 0.5; % Accepted: finite scalar in (0,1). Matches the validated continuous-loop test.
scheduledOnsetBoundaryPreflightEnabled = true; % Accepted: logical. Invisibly renders and erases frame 1 once before each run's first retained onset.
scheduledOnsetBoundaryMinArmLeadSec = 0.003; % Accepted: finite scalar > 0. Legacy name retained; checked before every retained 71/81/91 onset.
scheduledOnsetBoundaryMaxRetries = 6; % Accepted: integer >= 0. Legacy name retained; maximum low-patch rollover frames before recording exhaustion and continuing.
displayRefreshToleranceHz = 0.50; % Accepted: finite scalar >= 0. Allowed deviation from the input stimulus rate.
displayRefreshAbortEnabled = true; % Accepted: logical. true stops before acquisition when measured PTB/DataPixx refresh is outside tolerance.

%% PTB sync policy override and flip-miss guardrails
% Data flow: the lab-tested SyncTestSettings are applied before OpenWindow,
% saved in Conf, and restored during cleanup. SkipSyncTests remains disabled
% for acquisition; these thresholds keep PTB's synchronization checks active.
ptbSyncTestSettingsEnabled = true; % Lab default: apply the accepted PTB synchronization thresholds.
ptbSyncTestMaxStddevSec = 0.0005; % Lab-tested maximum accepted refresh-interval standard deviation (sec).
ptbSyncTestMinSamples = 50; % Lab-tested minimum valid timing samples.
ptbSyncTestMaxDeviationFraction = 0.1; % Lab-tested fractional deviation from nominal refresh interval.
ptbSyncTestMaxDurationSec = 5; % Lab-tested maximum duration of one sync-test run (sec).
debugSkipSyncTests = []; % Accepted: [] or {0,1}/logical. Default when []/unset: policy mode (MEG->schedule skipSyncTestsWhenMEG=false; debug non-MEG may use skipSyncTestsInDebug=true).
flipMissedWarnThresholdCount = []; % Accepted: [] or integer >= 1. Default when []/unset: scheduleCfg.flipMissedWarnThresholdCount (current default: 1).
flipMissedAbortEnabled = []; % Accepted: [] or logical. Default when []/unset: scheduleCfg.flipMissedAbortEnabled (current default: false).
flipMissedAbortThresholdCount = []; % Accepted: [] or integer >= 1. Default when []/unset: scheduleCfg.flipMissedAbortThresholdCount (current default: 5).

%% Eye-tracker runtime controls
% Data flow: overrides for EyeLink/fake/disabled gaze monitoring behavior.
trackeye = true; % Lab default: real EyeLink gaze stream enabled.
ignoreEyeTracker = false; % Lab default: do not bypass eye tracking.
fakeEyeTracker = false; % Lab default: fake gaze disabled in MEG sessions.
fakeGazeBreakRate = []; % Accepted: [] or finite scalar in [0,1]. Default when []/unset: scheduleCfg.fakeGazeBreakRate (current default: 0.05).
enableFixationAbort = true; % Lab default: enable fixation-break replay behavior.
fixWindowDeg = 3.0; % Lab default fixation threshold in degrees. Edit this value to change the threshold.
fixBreakToleranceFrames = []; % Accepted: [] or integer >= 1. Default when []/unset: scheduleCfg.fixBreakToleranceFrames (current default: 5).
fixationWarningDurationSec = []; % Accepted: [] or finite scalar >= 0 (sec). Default when []/unset: scheduleCfg.fixationWarningDurationSec (current default: 2.0).
replayLimit = 1; % Lab default: at most one replay per source trial.
allowInfiniteReplays = false; % Lab default: keep replay behavior bounded.
eyeTrackerDoCalibration = true; % Lab default: run calibration at session start.
eyeTrackerAllowBetweenBlockCalibration = true; % Lab default: allow recalibration between blocks.
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
fprintf(['Loaded runtime preset for MoveDot1 v24 all-trial guard: ' ...
    'SUB=%02d startBlock=%02d practice=%d dryRun=%d debugRun=%s ' ...
    'ITI=[%.3f %.3f] entranceStyle=%s diodeEnabled=%s size=%s position=%s ' ...
    'warmupCycles=%d warmupCode=%d guardLead=%.4fs retries=%d ' ...
    'ptbSyncSettings=%d refreshTolerance=%.3fHz\n'], ...
    iSub, iBlock, practice_mode, double(dryRunValidateScheduleOnly), ...
    mat2str(debug_runNumber), itiRangeForPrint(1), itiRangeForPrint(2), ...
    char(pathBandEntranceStyle), mat2str(photodiodePatchEnabled), ...
    mat2str(photodiodePatchSizePx), mat2str(photodiodePatchPositionPx), ...
    scheduledOnsetWarmupCycles, scheduledOnsetWarmupTriggerCode, ...
    scheduledOnsetBoundaryMinArmLeadSec, scheduledOnsetBoundaryMaxRetries, ...
    double(ptbSyncTestSettingsEnabled), displayRefreshToleranceHz);
