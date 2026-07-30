% runTestExperiment.m
%
% Purpose:
%   End-to-end local test launcher for the current V30/v23/v21 experiment
%   stack using a laptop-safe runtime profile:
%     - no 11MEG
%     - fake eye-tracker mode (no EyeLink hardware required)
%     - participant 99
%
% What this script does:
%   1) Generates subject-specific trajectories (V30) for participant 99.
%   2) Builds trial schedule and catch plan (v23) for participant 99.
%   3) Launches the v21 photodiode runtime with fake eye-tracker enabled.
%
% Usage example (from experiment/):
%   run('runTestExperiment.m');
%11
% Usage example (absolute path, MATLAB batch):
%   arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    run('runTestExperiment.m');"
%
% Inputs:
%   None required. Edit the toggles below if you want to skip generation or
%   scheduling and only launch the runtime.
%
% Outputs:
%   - input_files/MovDot_Sub99_V30_rescueTraject*.mat
%   - input_files/Sub99_TrialStruct_v23_rescueTraject.mat
%   - output_files/MoveDot1_occlusion_v21_photodiode_SUB99_BLOCKXX.mat
%     (for each completed block)
%
% Key assumptions:
%   - Script is run from experiment/ (it self-corrects with scriptDir).
%   - Psychtoolbox is installed for runtime display.
%   - Fake eye-tracker mode is sufficient for logic/UI testing.

%% Test toggles and subject/session identity
% Data flow: local toggles -> generation/schedule/runtime script calls.
% NOTE:
%   Downstream scripts use clearvars because they are script-based, so
%   control variables from this launcher may be removed between steps.
%   Keep toggles as literal booleans in each step and re-assign per-step
%   inputs right before each run('...') call.
if ~exist('doGenerateStimuli', 'var')
    doGenerateStimuli = true;
end
if ~exist('doGenerateSchedule', 'var')
    doGenerateSchedule = true;
end

% Participant/session identity for test runs.
% Keep explicit literals below as well, because called scripts clear caller
% workspace variables.
testSub = 98;
testStartBlock = 1;
testViewingDistanceMm = 1000;

% Keep full mode by default so all run transitions and block-end messaging
% can be tested in a realistic flow.
practice_mode = 0;

%% Runtime behavior profile (safe local test)
% Data flow: workspace overrides -> v21 runtime Conf.* fields.
trackeye = true;
ignoreEyeTracker = false;
fakeEyeTracker = true;
eyeTrackerDoCalibration = false;

% Keep gaze-break logic available, but disable random fake breaks for a
% predictable test trajectory flow.
enableFixationAbort = true;
fakeGazeBreakRate = 0.00;

% New timing/visual settings under test.
ITIRangeSec = [1.000, 1.500];
itiVisualMode = 'fixation_short';
itiGreyColor = [70, 70, 70];
itiVisualSeparationDur = [0.300, 0.700];
miniBlinkBreakEnabled = true;
miniBlinkBreakEveryNTrials = 10;
miniBlinkBreakDurationSec = 5.0;

% Centered catch-question defaults under test.
catchQuestionTopYFrac = 0.42;
catchQuestionBottomYFrac = 0.58;
catchBottomLeftXFrac = 0.38;
catchBottomRightXFrac = 0.62;

% Optional runtime toggles.
run1TransitionCountdownEnabled = true;
debugTriggerMode = 0;
debugSkipSyncTests = 1; % Local test convenience on non-calibrated displays.
if ~exist('runtimeDryRunValidateScheduleOnly', 'var')
    runtimeDryRunValidateScheduleOnly = false; % Set true to validate schedule without opening PTB windows.
end

overwriteExisting = true;
if ~exist('deleteExistingRuntimeOutput', 'var')
    deleteExistingRuntimeOutput = true;
end

%% Snapshot launcher options across script-level clearvars calls
% Data flow: launcher workspace -> temporary MAT snapshot -> rehydrated
% options before each step. This guards against child scripts that execute
% clearvars -except ... and would otherwise remove control flags.
launcherCfg = struct();
launcherCfg.doGenerateStimuli = logical(doGenerateStimuli);
launcherCfg.doGenerateSchedule = logical(doGenerateSchedule);
launcherCfg.testSub = testSub;
launcherCfg.testStartBlock = testStartBlock;
launcherCfg.testViewingDistanceMm = testViewingDistanceMm;
launcherCfg.practice_mode = practice_mode;
launcherCfg.trackeye = trackeye;
launcherCfg.ignoreEyeTracker = ignoreEyeTracker;
launcherCfg.fakeEyeTracker = fakeEyeTracker;
launcherCfg.eyeTrackerDoCalibration = eyeTrackerDoCalibration;
launcherCfg.enableFixationAbort = enableFixationAbort;
launcherCfg.fakeGazeBreakRate = fakeGazeBreakRate;
launcherCfg.ITIRangeSec = ITIRangeSec;
launcherCfg.itiVisualMode = itiVisualMode;
launcherCfg.itiGreyColor = itiGreyColor;
launcherCfg.itiVisualSeparationDur = itiVisualSeparationDur;
launcherCfg.miniBlinkBreakEnabled = miniBlinkBreakEnabled;
launcherCfg.miniBlinkBreakEveryNTrials = miniBlinkBreakEveryNTrials;
launcherCfg.miniBlinkBreakDurationSec = miniBlinkBreakDurationSec;
launcherCfg.catchQuestionTopYFrac = catchQuestionTopYFrac;
launcherCfg.catchQuestionBottomYFrac = catchQuestionBottomYFrac;
launcherCfg.catchBottomLeftXFrac = catchBottomLeftXFrac;
launcherCfg.catchBottomRightXFrac = catchBottomRightXFrac;
launcherCfg.run1TransitionCountdownEnabled = run1TransitionCountdownEnabled;
launcherCfg.debugTriggerMode = debugTriggerMode;
launcherCfg.debugSkipSyncTests = debugSkipSyncTests;
launcherCfg.runtimeDryRunValidateScheduleOnly = logical(runtimeDryRunValidateScheduleOnly);
launcherCfg.deleteExistingRuntimeOutput = logical(deleteExistingRuntimeOutput);
launcherCfg.overwriteExisting = logical(overwriteExisting);
launcherCfg.randomSeed = 9901;
launcherCfg.numBlocks = 1;

launcherCfgPath = fullfile(tempdir, sprintf('runTestExperiment_cfg_pid%d.mat', feature('getpid')));
save(launcherCfgPath, 'launcherCfg');

%% Resolve script directory and enter experiment root
% Data flow: script location -> deterministic cwd for relative paths.
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end
addpath('lib');

%% Step 1: Generate V30 stimuli for participant 99 (optional)
% Data flow: targetSubjectID + overwriteExisting -> input_files/*.mat
launcherCfgPath = fullfile(tempdir, sprintf('runTestExperiment_cfg_pid%d.mat', feature('getpid')));
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;
if launcherCfg.doGenerateStimuli
    targetSubjectID = launcherCfg.testSub;
    showPathPlots = false;
    overwriteExisting = launcherCfg.overwriteExisting;
    run('stimuli_generation_V30_rescueTraject.m');
end

%% Step 2: Create v23 trial schedule for participant 99 (optional)
% Data flow: iSub/randomSeed/numBlocks -> TrialStruct + CatchPlan MAT.
launcherCfgPath = fullfile(tempdir, sprintf('runTestExperiment_cfg_pid%d.mat', feature('getpid')));
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;
if launcherCfg.doGenerateSchedule
    iSub = launcherCfg.testSub;
    randomSeed = launcherCfg.randomSeed;
    numBlocks = launcherCfg.numBlocks;
    overwriteExisting = launcherCfg.overwriteExisting;
    run('CreateInputFiles_v23_rescueTraject.m');
end

%% Step 3: Launch v21 photodiode runtime experiment
% Data flow: workspace overrides above -> runtime script.
launcherCfgPath = fullfile(tempdir, sprintf('runTestExperiment_cfg_pid%d.mat', feature('getpid')));
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;
iSub = launcherCfg.testSub;
iBlock = launcherCfg.testStartBlock;
viewingDistanceMm = launcherCfg.testViewingDistanceMm;
practice_mode = launcherCfg.practice_mode;
trackeye = launcherCfg.trackeye;
ignoreEyeTracker = launcherCfg.ignoreEyeTracker;
fakeEyeTracker = launcherCfg.fakeEyeTracker;
eyeTrackerDoCalibration = launcherCfg.eyeTrackerDoCalibration;
enableFixationAbort = launcherCfg.enableFixationAbort;
fakeGazeBreakRate = launcherCfg.fakeGazeBreakRate;
ITIRangeSec = launcherCfg.ITIRangeSec;
itiVisualMode = launcherCfg.itiVisualMode;
itiGreyColor = launcherCfg.itiGreyColor;
itiVisualSeparationDur = launcherCfg.itiVisualSeparationDur;
miniBlinkBreakEnabled = launcherCfg.miniBlinkBreakEnabled;
miniBlinkBreakEveryNTrials = launcherCfg.miniBlinkBreakEveryNTrials;
miniBlinkBreakDurationSec = launcherCfg.miniBlinkBreakDurationSec;
catchQuestionTopYFrac = launcherCfg.catchQuestionTopYFrac;
catchQuestionBottomYFrac = launcherCfg.catchQuestionBottomYFrac;
catchBottomLeftXFrac = launcherCfg.catchBottomLeftXFrac;
catchBottomRightXFrac = launcherCfg.catchBottomRightXFrac;
run1TransitionCountdownEnabled = launcherCfg.run1TransitionCountdownEnabled;
debugTriggerMode = launcherCfg.debugTriggerMode;
debugSkipSyncTests = launcherCfg.debugSkipSyncTests;
dryRunValidateScheduleOnly = logical(launcherCfg.runtimeDryRunValidateScheduleOnly);
overwriteExisting = launcherCfg.overwriteExisting;
deleteExistingRuntimeOutput = launcherCfg.deleteExistingRuntimeOutput;

if deleteExistingRuntimeOutput
    outputMat = fullfile('output_files', sprintf( ...
        'MoveDot1_occlusion_v21_photodiode_SUB%02d_BLOCK%02d.mat', ...
        iSub, iBlock));
    if exist(outputMat, 'file') == 2
        delete(outputMat);
    end
end
run('MoveDot1_experiment_occlusion_v21_photodiode.m');

% Best-effort cleanup of temporary launcher state.
launcherCfgPath = fullfile(tempdir, sprintf('runTestExperiment_cfg_pid%d.mat', feature('getpid')));
if exist(launcherCfgPath, 'file') == 2
    delete(launcherCfgPath);
end
