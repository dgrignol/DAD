% runExperiment.m
%
% Purpose:
%   One-command launcher for the current real experiment stack:
%     1) stimuli_generation_V30_rescueTraject.m
%     2) CreateInputFiles_v23_rescueTraject.m
%     3) MoveDot1_experiment_occlusion_v20_rescueTraject.m
%
%   The launcher is configured for real acquisition mode:
%     - MEG/DataPixx enabled
%     - real EyeLink enabled (no fake eye tracker)
%     - EyeLink calibration enabled
%
%   Running this script should be enough to prepare inputs and start the
%   experiment flow without manual per-step calls.
%
% Usage example (interactive from experiment/):
%   run('runExperiment.m');
%
% Usage example (batch from shell):
%   arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    run('runExperiment.m');"
%
% Optional workspace overrides (set before running this script):
%   - doGenerateStimuli (logical, default true)
%   - doGenerateSchedule (logical, default true)
%   - overwriteExistingInputs (logical, default false)
%   - showPathPlots (logical, default false)
%   - randomSeed (integer, default iSub*100 + 1)
%   - numBlocks (integer >= 1, default schedule config value)
%   - runtimeConfigScript (path, default 'lib/Config_runtime_v20_rescueTraject.m')
%
% Inputs:
%   - lib/Config_runtime_v20_rescueTraject.m (runtime preset)
%   - lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m
%   - lib/Config_stimuli_generation_V30_rescueTraject.m
%
% Outputs:
%   - input_files/MovDot_SubXX_V30_rescueTraject*.mat
%   - input_files/SubXX_TrialStruct_v23_rescueTraject.mat
%   - output_files/MoveDot1_occlusion_v20_rescueTraject_SUBXX_BLOCKYY.mat
%   - output_files/MBssbbaa.edf (one EDF per completed block)
%
% Key assumptions:
%   - Psychtoolbox, DataPixx, and EyeLink toolboxes are installed on the
%     acquisition machine.
%   - This script is launched from MATLAB with access to experiment/lib.
%   - Subject/block values are provided by the runtime config preset.

%% Launcher toggles and defaults
% Data flow: optional caller overrides -> deterministic launcher behavior.
if ~exist('doGenerateStimuli', 'var') || isempty(doGenerateStimuli)
    doGenerateStimuli = true;
end
if ~exist('doGenerateSchedule', 'var') || isempty(doGenerateSchedule)
    doGenerateSchedule = true;
end
if ~exist('overwriteExistingInputs', 'var') || isempty(overwriteExistingInputs)
    overwriteExistingInputs = false;
end
if ~exist('showPathPlots', 'var') || isempty(showPathPlots)
    showPathPlots = false;
end
if ~exist('randomSeed', 'var')
    randomSeed = [];
end
if ~exist('numBlocks', 'var')
    numBlocks = [];
end
if ~exist('runtimeConfigScript', 'var') || isempty(runtimeConfigScript)
    runtimeConfigScript = fullfile('lib', 'Config_runtime_v20_rescueTraject.m');
end

%% Resolve launcher location and load runtime preset
% Data flow: script location + runtime preset file -> session identity and
% runtime override variables in workspace.
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end
addpath('lib');

if exist(runtimeConfigScript, 'file') ~= 2
    error('Runtime config script not found: %s', runtimeConfigScript);
end
run(runtimeConfigScript);

%% Validate required runtime identity and derive schedule defaults
% Data flow: runtime preset + schedule defaults -> complete generation and
% runtime input set (subject, block, seed, block count).
if ~exist('iSub', 'var') || isempty(iSub) || ~isfinite(iSub)
    error('Runtime preset must define finite iSub.');
end
if ~exist('iBlock', 'var') || isempty(iBlock) || ~isfinite(iBlock)
    error('Runtime preset must define finite iBlock.');
end
if ~exist('viewingDistanceMm', 'var') || isempty(viewingDistanceMm) || ~isfinite(viewingDistanceMm)
    error('Runtime preset must define finite viewingDistanceMm.');
end

iSub = round(double(iSub));
iBlock = round(double(iBlock));
viewingDistanceMm = double(viewingDistanceMm);

scheduleCfg = Config_schedule_CreateInputV23_MoveDotV20_rescueTraject;
if isempty(randomSeed)
    randomSeed = iSub * 100 + 1;
end
if isempty(numBlocks)
    numBlocks = double(scheduleCfg.numBlocks);
end
randomSeed = round(double(randomSeed));
numBlocks = round(double(numBlocks));
if numBlocks < 1
    error('numBlocks must be >= 1.');
end

%% Force real acquisition profile (MEG + real EyeLink)
% Data flow: launcher policy -> workspace overrides consumed by runtime.
MEG = true;
trackeye = true;
ignoreEyeTracker = false;
fakeEyeTracker = false;
eyeTrackerDoCalibration = true;
eyeTrackerAllowBetweenBlockCalibration = true;
dryRunValidateScheduleOnly = false;
debugSkipSyncTests = 0;

fprintf(['Real acquisition profile: MEG=%d trackeye=%d ignoreEyeTracker=%d ' ...
    'fakeEyeTracker=%d calibration=%d\n'], ...
    double(MEG), double(trackeye), double(ignoreEyeTracker), ...
    double(fakeEyeTracker), double(eyeTrackerDoCalibration));

%% Snapshot launcher state across script-level clearvars calls
% Data flow: launcher workspace -> temp MAT -> reloaded before each step.
launcherCfg = struct();
launcherCfg.doGenerateStimuli = logical(doGenerateStimuli);
launcherCfg.doGenerateSchedule = logical(doGenerateSchedule);
launcherCfg.overwriteExistingInputs = logical(overwriteExistingInputs);
launcherCfg.showPathPlots = logical(showPathPlots);
launcherCfg.iSub = iSub;
launcherCfg.iBlock = iBlock;
launcherCfg.viewingDistanceMm = viewingDistanceMm;
launcherCfg.randomSeed = randomSeed;
launcherCfg.numBlocks = numBlocks;

% Runtime overrides passed to MoveDot1_experiment_occlusion_v20_rescueTraject.m
launcherCfg.MEG = logical(MEG);
launcherCfg.practice_mode = practice_mode;
launcherCfg.dryRunValidateScheduleOnly = logical(dryRunValidateScheduleOnly);
launcherCfg.debug_runNumber = debug_runNumber;
launcherCfg.dryRunSimulatedBreakKeys = dryRunSimulatedBreakKeys;
launcherCfg.question_text = question_text;
launcherCfg.debugTriggerMode = debugTriggerMode;
launcherCfg.runColorCueEnabled = runColorCueEnabled;
launcherCfg.practiceRunPercentageTrials = practiceRunPercentageTrials;
launcherCfg.ITIRangeSec = ITIRangeSec;
launcherCfg.itiVisualMode = itiVisualMode;
launcherCfg.itiGreyColor = itiGreyColor;
launcherCfg.itiVisualSeparationDur = itiVisualSeparationDur;
launcherCfg.catchQuestionTopYFrac = catchQuestionTopYFrac;
launcherCfg.catchQuestionBottomYFrac = catchQuestionBottomYFrac;
launcherCfg.catchBottomLeftXFrac = catchBottomLeftXFrac;
launcherCfg.catchBottomRightXFrac = catchBottomRightXFrac;
launcherCfg.miniBlinkBreakEnabled = miniBlinkBreakEnabled;
launcherCfg.miniBlinkBreakEveryNTrials = miniBlinkBreakEveryNTrials;
launcherCfg.miniBlinkBreakDurationSec = miniBlinkBreakDurationSec;
launcherCfg.pathBandEntranceStyle = pathBandEntranceStyle;
launcherCfg.pathBandStraightBackshiftDotRadiusScale = pathBandStraightBackshiftDotRadiusScale;
launcherCfg.pathBandStraightEndForwardDotRadiusScale = pathBandStraightEndForwardDotRadiusScale;
launcherCfg.pathBandStraightPostStartWidthExtraDotRadiusScale = pathBandStraightPostStartWidthExtraDotRadiusScale;
launcherCfg.pathBandStraightStartBackshiftPixelOffset = pathBandStraightStartBackshiftPixelOffset;
launcherCfg.debugSkipSyncTests = debugSkipSyncTests;
launcherCfg.flipMissedWarnThresholdCount = flipMissedWarnThresholdCount;
launcherCfg.flipMissedAbortEnabled = flipMissedAbortEnabled;
launcherCfg.flipMissedAbortThresholdCount = flipMissedAbortThresholdCount;
launcherCfg.trackeye = logical(trackeye);
launcherCfg.ignoreEyeTracker = logical(ignoreEyeTracker);
launcherCfg.fakeEyeTracker = logical(fakeEyeTracker);
launcherCfg.fakeGazeBreakRate = fakeGazeBreakRate;
launcherCfg.enableFixationAbort = enableFixationAbort;
launcherCfg.fixWindowDeg = fixWindowDeg;
launcherCfg.fixBreakToleranceFrames = fixBreakToleranceFrames;
launcherCfg.fixationWarningDurationSec = fixationWarningDurationSec;
launcherCfg.replayLimit = replayLimit;
launcherCfg.allowInfiniteReplays = allowInfiniteReplays;
launcherCfg.eyeTrackerDoCalibration = logical(eyeTrackerDoCalibration);
launcherCfg.eyeTrackerAllowBetweenBlockCalibration = logical(eyeTrackerAllowBetweenBlockCalibration);
launcherCfg.eyeTrackerSendTriggerMessages = eyeTrackerSendTriggerMessages;
launcherCfg.eyeTrackerImageTransferEnabled = eyeTrackerImageTransferEnabled;
launcherCfg.run1TransitionCountdownEnabled = run1TransitionCountdownEnabled;

launcherCfgPath = fullfile(tempdir, sprintf('runExperiment_cfg_pid%d.mat', feature('getpid')));
save(launcherCfgPath, 'launcherCfg');

%% Step 1: Generate V30 trajectories (optional)
% Data flow: subject id + generation flags -> MovDot_SubXX_V30_rescueTraject MAT files.
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;
if launcherCfg.doGenerateStimuli
    targetSubjectID = launcherCfg.iSub;
    overwriteExisting = launcherCfg.overwriteExistingInputs;
    showPathPlots = launcherCfg.showPathPlots;
    run('stimuli_generation_V30_rescueTraject.m');
end

%% Step 2: Generate v23 schedule/catch plan (optional)
% Data flow: subject + seed + block count -> SubXX_TrialStruct_v23_rescueTraject.mat.
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;
if launcherCfg.doGenerateSchedule
    iSub = launcherCfg.iSub; %#ok<NASGU>  % consumed by CreateInputFiles script workspace
    randomSeed = launcherCfg.randomSeed;
    numBlocks = launcherCfg.numBlocks;
    overwriteExisting = launcherCfg.overwriteExistingInputs;
    run('CreateInputFiles_v23_rescueTraject.m');
end

%% Step 3: Launch runtime experiment (v20)
% Data flow: runtime overrides + prepared input files -> experiment execution.
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;

iSub = launcherCfg.iSub;
iBlock = launcherCfg.iBlock;
viewingDistanceMm = launcherCfg.viewingDistanceMm;
MEG = launcherCfg.MEG;
practice_mode = launcherCfg.practice_mode;
dryRunValidateScheduleOnly = launcherCfg.dryRunValidateScheduleOnly;
debug_runNumber = launcherCfg.debug_runNumber;
dryRunSimulatedBreakKeys = launcherCfg.dryRunSimulatedBreakKeys;
question_text = launcherCfg.question_text;
debugTriggerMode = launcherCfg.debugTriggerMode;
runColorCueEnabled = launcherCfg.runColorCueEnabled;
practiceRunPercentageTrials = launcherCfg.practiceRunPercentageTrials;
ITIRangeSec = launcherCfg.ITIRangeSec;
itiVisualMode = launcherCfg.itiVisualMode;
itiGreyColor = launcherCfg.itiGreyColor;
itiVisualSeparationDur = launcherCfg.itiVisualSeparationDur;
catchQuestionTopYFrac = launcherCfg.catchQuestionTopYFrac;
catchQuestionBottomYFrac = launcherCfg.catchQuestionBottomYFrac;
catchBottomLeftXFrac = launcherCfg.catchBottomLeftXFrac;
catchBottomRightXFrac = launcherCfg.catchBottomRightXFrac;
miniBlinkBreakEnabled = launcherCfg.miniBlinkBreakEnabled;
miniBlinkBreakEveryNTrials = launcherCfg.miniBlinkBreakEveryNTrials;
miniBlinkBreakDurationSec = launcherCfg.miniBlinkBreakDurationSec;
pathBandEntranceStyle = launcherCfg.pathBandEntranceStyle;
pathBandStraightBackshiftDotRadiusScale = launcherCfg.pathBandStraightBackshiftDotRadiusScale;
pathBandStraightEndForwardDotRadiusScale = launcherCfg.pathBandStraightEndForwardDotRadiusScale;
pathBandStraightPostStartWidthExtraDotRadiusScale = launcherCfg.pathBandStraightPostStartWidthExtraDotRadiusScale;
pathBandStraightStartBackshiftPixelOffset = launcherCfg.pathBandStraightStartBackshiftPixelOffset;
debugSkipSyncTests = launcherCfg.debugSkipSyncTests;
flipMissedWarnThresholdCount = launcherCfg.flipMissedWarnThresholdCount;
flipMissedAbortEnabled = launcherCfg.flipMissedAbortEnabled;
flipMissedAbortThresholdCount = launcherCfg.flipMissedAbortThresholdCount;
trackeye = launcherCfg.trackeye;
ignoreEyeTracker = launcherCfg.ignoreEyeTracker;
fakeEyeTracker = launcherCfg.fakeEyeTracker;
fakeGazeBreakRate = launcherCfg.fakeGazeBreakRate;
enableFixationAbort = launcherCfg.enableFixationAbort;
fixWindowDeg = launcherCfg.fixWindowDeg;
fixBreakToleranceFrames = launcherCfg.fixBreakToleranceFrames;
fixationWarningDurationSec = launcherCfg.fixationWarningDurationSec;
replayLimit = launcherCfg.replayLimit;
allowInfiniteReplays = launcherCfg.allowInfiniteReplays;
eyeTrackerDoCalibration = launcherCfg.eyeTrackerDoCalibration;
eyeTrackerAllowBetweenBlockCalibration = launcherCfg.eyeTrackerAllowBetweenBlockCalibration;
eyeTrackerSendTriggerMessages = launcherCfg.eyeTrackerSendTriggerMessages;
eyeTrackerImageTransferEnabled = launcherCfg.eyeTrackerImageTransferEnabled;
run1TransitionCountdownEnabled = launcherCfg.run1TransitionCountdownEnabled;

run('MoveDot1_experiment_occlusion_v20_rescueTraject.m');

%% Cleanup temporary launcher snapshot
% Data flow: best-effort deletion of temp state file.
launcherCfgPath = fullfile(tempdir, sprintf('runExperiment_cfg_pid%d.mat', feature('getpid')));
if exist(launcherCfgPath, 'file') == 2
    delete(launcherCfgPath);
end
