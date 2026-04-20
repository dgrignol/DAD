%% run_practice_occlusion_v1
% Script: run_practice_occlusion_v1.m
%
% Purpose:
%   Run an isolated practice session for the one-dot occlusion experiment
%   without MEG hardware, using one shared input dataset for all participants.
%
%   Practice structure:
%   - Run 1: catch-only practice trials derived from always_visible trials.
%   - Run 2: catch-only practice trials derived from occluded trials.
%   - Run 3 is intentionally skipped in practice.
%
%   This script keeps all inputs/outputs inside experiment/practice so it
%   does not alter real-run experiment files.
%
% What this script does:
%   1) Generates one shared practice input/schedule (if missing) using the
%      latest V28/v21 generation scripts in practice-local paths.
%   2) Asks participant number and whether to use EyeLink.
%   3) Runs practice run 1 and run 2 with independently tunable trial counts.
%   4) Shows a clear transition message between run 1 and run 2.
%   5) Lets the operator repeat practice with a key/button-style choice.
%   6) Saves one consolidated MAT file with catch accuracy and RT metrics.
%   7) Prints summary metrics by catch type for run 1 and run 2.
%
% Usage example (interactive from experiment/practice):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/practice');
%   run_practice_occlusion_v1;
%
% Usage example (non-interactive dry-run test, no PTB window):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/practice'); ...
%    practiceParticipantNumber=1; practiceUseEyeTracker=false; ...
%    practiceDryRunValidateScheduleOnly=true; practiceAutoRepeatCount=1; ...
%    practiceRun1TrialCount=4; practiceRun2TrialCount=4; ...
%    run_practice_occlusion_v1;"
%
% Inputs (workspace overrides, optional):
%   - practiceParticipantNumber        : positive integer participant ID.
%   - practiceUseEyeTracker            : logical, true enables EyeLink.
%   - practiceRun1TrialCount           : integer >=1, target practice trials for run 1.
%   - practiceRun2TrialCount           : integer >=1, target practice trials for run 2.
%   - practiceViewingDistanceMm        : positive scalar, default 1000.
%   - practiceDryRunValidateScheduleOnly: logical, true skips PTB rendering.
%   - practiceTargetRefreshHz          : numeric scalar or []. When [] the
%     launcher detects refresh from system info and derives target fps.
%   - practiceStrictSafeMode           : logical, default true. When true,
%     non-dry runs are blocked unless strict safety preflight passes.
%   - practiceSkipSyncTests            : logical, default true in practice.
%   - practiceDisablePriorityBoost     : logical, default true in practice.
%   - practiceAutoContinueBetweenRuns  : logical, default false; when true,
%     skip the run1->run2 wait-for-key gate (useful for automated testing).
%   - practiceAutoRepeatCount          : integer >=1 for non-interactive repeats.
%
% Outputs:
%   - practice/output_files/Practice_SubXX_YYYYmmdd_HHMMSS.mat
%     containing one struct `practiceResults` with run metrics and summaries.
%
% Assumptions:
%   - Parent scripts exist one level up:
%       MoveDot1_experiment_occlusion_v18_rescueTraject.m (practice-local copy)
%       stimuli_generation_V28_rescueTraject.m (practice-local copy)
%       CreateInputFiles_v21_rescueTraject.m (practice-local copy)
%   - Psychtoolbox/EyeLink availability follows your local setup.

function run_practice_occlusion_v1

%% Practice parameters and defaults
% Data flow: fixed practice constants + optional workspace overrides ->
% runtime controls for generation and execution.
practiceSharedInputSubjectId = local_get_base_or_default('practiceSharedInputSubjectId', 91);
practiceSharedRandomSeed = local_get_base_or_default('practiceSharedRandomSeed', 9101);
practiceRun1TrialCount = local_get_base_or_default('practiceRun1TrialCount', 8);
practiceRun2TrialCount = local_get_base_or_default('practiceRun2TrialCount', 8);
practiceViewingDistanceMm = local_get_base_or_default('practiceViewingDistanceMm', 1000);
practiceDryRunValidateScheduleOnly = local_get_base_or_default('practiceDryRunValidateScheduleOnly', false);
practiceTargetRefreshHz = local_get_base_or_default('practiceTargetRefreshHz', []);
practiceStrictSafeMode = local_get_base_or_default('practiceStrictSafeMode', true);
practiceSkipSyncTests = local_get_base_or_default('practiceSkipSyncTests', true);
practiceDisablePriorityBoost = local_get_base_or_default('practiceDisablePriorityBoost', true);
practiceAutoContinueBetweenRuns = local_get_base_or_default('practiceAutoContinueBetweenRuns', false);
practiceAutoRepeatCount = local_get_base_or_default('practiceAutoRepeatCount', []);
practiceParticipantNumber = local_get_base_or_default('practiceParticipantNumber', []);
practiceUseEyeTracker = local_get_base_or_default('practiceUseEyeTracker', []);

practiceRun1TrialCount = round(double(practiceRun1TrialCount));
practiceRun2TrialCount = round(double(practiceRun2TrialCount));
practiceViewingDistanceMm = double(practiceViewingDistanceMm);
practiceDryRunValidateScheduleOnly = logical(practiceDryRunValidateScheduleOnly);
if ~isempty(practiceTargetRefreshHz)
    practiceTargetRefreshHz = double(practiceTargetRefreshHz);
end
practiceStrictSafeMode = logical(practiceStrictSafeMode);
practiceSkipSyncTests = logical(practiceSkipSyncTests);
practiceDisablePriorityBoost = logical(practiceDisablePriorityBoost);
practiceAutoContinueBetweenRuns = logical(practiceAutoContinueBetweenRuns);

if practiceRun1TrialCount < 1 || practiceRun2TrialCount < 1
    error('practiceRun1TrialCount and practiceRun2TrialCount must be >= 1.');
end
if ~isfinite(practiceViewingDistanceMm) || practiceViewingDistanceMm <= 0
    error('practiceViewingDistanceMm must be a positive finite scalar.');
end
if ~isempty(practiceAutoRepeatCount)
    practiceAutoRepeatCount = round(double(practiceAutoRepeatCount));
    if ~isfinite(practiceAutoRepeatCount) || practiceAutoRepeatCount < 1
        error('practiceAutoRepeatCount must be empty or an integer >= 1.');
    end
end

%% Resolve script paths and prep environment
% Data flow: script location -> practice folder paths and MATLAB path.
thisScriptPath = mfilename('fullpath');
if isempty(thisScriptPath)
    error('Could not resolve script path for run_practice_occlusion_v1.m.');
end
practiceDir = fileparts(thisScriptPath);
libDir = fullfile(practiceDir, 'lib');
inputDir = fullfile(practiceDir, 'input_files');
outputDir = fullfile(practiceDir, 'output_files');
rawOutputDir = fullfile(outputDir, 'raw_runtime');

if ~exist(inputDir, 'dir')
    mkdir(inputDir);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(rawOutputDir, 'dir')
    mkdir(rawOutputDir);
end

addpath(libDir);

%% Collect participant-level runtime choices
% Data flow: operator input -> participant/session control variables.
if isempty(practiceParticipantNumber)
    practiceParticipantNumber = local_prompt_numeric( ...
        'Participant number for PRACTICE output file:', 1);
end
practiceParticipantNumber = round(double(practiceParticipantNumber));
if ~isfinite(practiceParticipantNumber) || practiceParticipantNumber < 1
    error('practiceParticipantNumber must be a positive integer.');
end

if isempty(practiceUseEyeTracker)
    practiceUseEyeTracker = local_prompt_yes_no( ...
        'Use EyeLink for practice (1 = yes, 0 = no)?', false);
end
practiceUseEyeTracker = logical(practiceUseEyeTracker);

%% Resolve display refresh and target generation fps
% Data flow: optional override + system display info -> target fps for
% input generation/caching.
[detectedRefreshHz, targetFrameFrequencyHz, refreshDetectionMethod] = ...
    local_resolve_practice_refresh_hz(practiceTargetRefreshHz);

%% Ensure shared practice inputs exist (generate if missing)
% Data flow: shared subject ID + detected target fps + local scripts ->
% fps-specific cached MAT inputs, then canonical runtime input files.
inputSelection = local_ensure_shared_practice_inputs( ...
    practiceDir, practiceSharedInputSubjectId, practiceSharedRandomSeed, targetFrameFrequencyHz);

%% Resolve run base lengths and map trial counts to practice percentages
% Data flow: TrialOrder sizes + requested counts -> exact percentage values
% consumed by v18 runtime practice-mode trial trimming.
trialStructPath = inputSelection.canonicalSchedulePath;
baseRun1Length = local_count_valid_run_trials(trialStructPath, 1, 1);
baseRun2Length = local_count_valid_run_trials(trialStructPath, 1, 2);

if practiceRun1TrialCount > baseRun1Length
    error('Requested run 1 practice trials (%d) exceed available run length (%d).', ...
        practiceRun1TrialCount, baseRun1Length);
end
if practiceRun2TrialCount > baseRun2Length
    error('Requested run 2 practice trials (%d) exceed available run length (%d).', ...
        practiceRun2TrialCount, baseRun2Length);
end

practiceRun1Pct = local_trials_to_percentage(practiceRun1TrialCount, baseRun1Length);
practiceRun2Pct = local_trials_to_percentage(practiceRun2TrialCount, baseRun2Length);

fprintf('\nPractice setup:\n');
fprintf('  Shared input subject: %d\n', practiceSharedInputSubjectId);
fprintf('  Participant number:   %d\n', practiceParticipantNumber);
fprintf('  Eye tracker:          %d\n', double(practiceUseEyeTracker));
fprintf('  Detected refresh:     %.2f Hz (%s)\n', detectedRefreshHz, refreshDetectionMethod);
fprintf('  Target generation fps:%d Hz\n', targetFrameFrequencyHz);
fprintf('  Active input file:    %s\n', inputSelection.canonicalObservedPath);
fprintf('  Run 1 target trials:  %d/%d (%.4f%%)\n', ...
    practiceRun1TrialCount, baseRun1Length, practiceRun1Pct);
fprintf('  Run 2 target trials:  %d/%d (%.4f%%)\n', ...
    practiceRun2TrialCount, baseRun2Length, practiceRun2Pct);
fprintf('  Dry-run mode:         %d\n\n', double(practiceDryRunValidateScheduleOnly));
fprintf('  Strict safe mode:     %d\n', double(practiceStrictSafeMode));
fprintf('  Skip sync tests:      %d\n\n', double(practiceSkipSyncTests));
fprintf('  Disable priority:     %d\n', double(practiceDisablePriorityBoost));
fprintf('  Auto-continue R1->R2: %d\n\n', double(practiceAutoContinueBetweenRuns));

% Section: strict preflight safety gate for non-dry runs.
% Data flow: launcher safety flags + local PTB environment checks -> allow
% or block full-screen runtime execution.
local_run_strict_safety_preflight( ...
    practiceStrictSafeMode, ...
    practiceDryRunValidateScheduleOnly, ...
    practiceSkipSyncTests, ...
    practiceDisablePriorityBoost);

%% Run practice attempts (run 1 then run 2) with optional repeat
% Data flow: participant controls + run-specific overrides -> per-attempt
% runtime files -> consolidated attempt summaries.
attemptSummaries = struct([]);
attemptIndex = 0;
continuePractice = true;

while continuePractice
    attemptIndex = attemptIndex + 1;
    fprintf('\n=== PRACTICE ATTEMPT %d ===\n', attemptIndex);

    run1Raw = local_execute_runtime_run( ...
        practiceDir, rawOutputDir, ...
        practiceSharedInputSubjectId, practiceViewingDistanceMm, ...
        1, practiceRun1Pct, practiceUseEyeTracker, ...
        true, practiceDryRunValidateScheduleOnly, practiceSkipSyncTests, ...
        practiceDisablePriorityBoost, attemptIndex);

    local_show_run1_to_run2_message( ...
        practiceDryRunValidateScheduleOnly, practiceAutoContinueBetweenRuns);

    run2Raw = local_execute_runtime_run( ...
        practiceDir, rawOutputDir, ...
        practiceSharedInputSubjectId, practiceViewingDistanceMm, ...
        2, practiceRun2Pct, practiceUseEyeTracker, ...
        false, practiceDryRunValidateScheduleOnly, practiceSkipSyncTests, ...
        practiceDisablePriorityBoost, attemptIndex);

    attemptSummary = local_build_attempt_summary( ...
        attemptIndex, run1Raw, run2Raw, practiceDryRunValidateScheduleOnly);
    attemptSummaries = local_append_attempt(attemptSummaries, attemptSummary);

    local_print_attempt_summary(attemptSummary);

    if ~isempty(practiceAutoRepeatCount)
        continuePractice = (attemptIndex < practiceAutoRepeatCount);
    else
        continuePractice = local_prompt_repeat_choice(practiceDryRunValidateScheduleOnly);
    end
end

%% Save one consolidated participant practice file
% Data flow: all attempt summaries + participant metadata -> single MAT file.
saveStamp = datestr(now, 'yyyymmdd_HHMMSS');
summaryFilePath = fullfile(outputDir, sprintf('Practice_Sub%02d_%s.mat', ...
    practiceParticipantNumber, saveStamp));

practiceResults = struct();
practiceResults.participantNumber = practiceParticipantNumber;
practiceResults.useEyeTracker = practiceUseEyeTracker;
practiceResults.sharedInputSubjectId = practiceSharedInputSubjectId;
practiceResults.sharedRandomSeed = practiceSharedRandomSeed;
practiceResults.practiceRun1TrialCount = practiceRun1TrialCount;
practiceResults.practiceRun2TrialCount = practiceRun2TrialCount;
practiceResults.practiceRun1Percentage = practiceRun1Pct;
practiceResults.practiceRun2Percentage = practiceRun2Pct;
practiceResults.viewingDistanceMm = practiceViewingDistanceMm;
practiceResults.detectedRefreshHz = detectedRefreshHz;
practiceResults.refreshDetectionMethod = refreshDetectionMethod;
practiceResults.targetFrameFrequencyHz = targetFrameFrequencyHz;
practiceResults.inputSelection = inputSelection;
practiceResults.dryRunValidateScheduleOnly = practiceDryRunValidateScheduleOnly;
practiceResults.strictSafeMode = practiceStrictSafeMode;
practiceResults.skipSyncTests = practiceSkipSyncTests;
practiceResults.disablePriorityBoost = practiceDisablePriorityBoost;
practiceResults.autoContinueBetweenRuns = practiceAutoContinueBetweenRuns;
practiceResults.attempts = attemptSummaries;
practiceResults.createdAtIso8601 = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

save(summaryFilePath, 'practiceResults', '-v7');
fprintf('\nSaved consolidated practice summary: %s\n', summaryFilePath);

end

%% Local helper functions
function value = local_prompt_numeric(promptText, defaultValue)
% LOCAL_PROMPT_NUMERIC Request one numeric value from the operator.
defaultStr = sprintf('%d', round(defaultValue));
answer = inputdlg({promptText}, 'Practice input', [1 70], {defaultStr});
if isempty(answer)
    error('Canceled by user.');
end
value = str2double(answer{1});
if ~isfinite(value)
    error('Invalid numeric input for: %s', promptText);
end
end

function flag = local_prompt_yes_no(promptText, defaultNo)
% LOCAL_PROMPT_YES_NO Request a yes/no logical choice with 1/0 entry.
if logical(defaultNo)
    defaultStr = '0';
else
    defaultStr = '1';
end
answer = inputdlg({promptText}, 'Practice input', [1 70], {defaultStr});
if isempty(answer)
    error('Canceled by user.');
end
value = str2double(answer{1});
if ~isfinite(value) || ~ismember(round(value), [0, 1])
    error('Invalid yes/no input. Use 1 for yes, 0 for no.');
end
flag = logical(round(value));
end

function inputSelection = local_ensure_shared_practice_inputs(practiceDir, sharedSubId, sharedSeed, targetFpsHz)
% LOCAL_ENSURE_SHARED_PRACTICE_INPUTS Select/generate fps-matched shared inputs.
%
% Data flow:
%   target fps -> fps-tagged cache files -> canonical runtime input files.
%   If cache is missing, generate canonical files at target fps and cache them.
addpath(fullfile(practiceDir, 'lib'));
scheduleCfg = Config_schedule_CreateInputV21_MoveDotV18_rescueTraject;

targetSubjectID = round(double(sharedSubId));
targetFpsHz = round(double(targetFpsHz));
fpsTag = sprintf('%03dHz', targetFpsHz);

canonicalObservedPath = fullfile(practiceDir, 'input_files', ...
    sprintf(scheduleCfg.inputFilePattern, targetSubjectID));
canonicalPredictedPath = fullfile(practiceDir, 'input_files', ...
    sprintf('MovDot_Sub%02d_V28_rescueTraject_predicted.mat', targetSubjectID));
canonicalSchedulePath = fullfile(practiceDir, 'input_files', ...
    sprintf(scheduleCfg.outputFilePattern, targetSubjectID));

cachedObservedPath = fullfile(practiceDir, 'input_files', ...
    sprintf('MovDot_Sub%02d_V28_rescueTraject_%s.mat', targetSubjectID, fpsTag));
cachedPredictedPath = fullfile(practiceDir, 'input_files', ...
    sprintf('MovDot_Sub%02d_V28_rescueTraject_predicted_%s.mat', targetSubjectID, fpsTag));
cachedSchedulePath = fullfile(practiceDir, 'input_files', ...
    sprintf('Sub%02d_TrialStruct_v21_rescueTraject_%s.mat', targetSubjectID, fpsTag));

haveCachedInput = isfile(cachedObservedPath) && ...
    isfile(cachedPredictedPath) && isfile(cachedSchedulePath) && ...
    local_input_file_matches_fps(cachedObservedPath, targetFpsHz);

if haveCachedInput
    % Use cached fps-specific files by syncing them into canonical names
    % consumed by runtime and schedule scripts.
    copyfile(cachedObservedPath, canonicalObservedPath, 'f');
    copyfile(cachedPredictedPath, canonicalPredictedPath, 'f');
    copyfile(cachedSchedulePath, canonicalSchedulePath, 'f');
    fprintf('Using cached practice inputs for %d Hz (%s).\n', targetFpsHz, fpsTag);
else
    fprintf(['Practice inputs for %d Hz not found (or stale). ' ...
        'Generating in %s ...\n'], targetFpsHz, fullfile(practiceDir, 'input_files'));

    practiceDirEsc = local_escape_single_quotes(practiceDir);
    stimScriptEsc = local_escape_single_quotes(fullfile(practiceDir, 'stimuli_generation_V28_rescueTraject.m'));

    % Section: generate shared V28 stimulus file with target fps override.
    % Data flow: target fps + shared subject -> canonical observed/predicted files.
    overwriteExisting = true;
    baseFpsHz = double(Config_stimuli_generation_V28_rescueTraject.frameFrequency);
    baseFixedDevFrame = double(Config_stimuli_generation_V28_rescueTraject.fixedDevianceFrame);
    baseFixedOccEndFrame = double(Config_stimuli_generation_V28_rescueTraject.fixedOcclusionEndFrame);
    scaledFramesPerTrial = round(double(Config_stimuli_generation_V28_rescueTraject.trialDuration) * targetFpsHz);
    scaledFixedDevFrame = max(1, round(baseFixedDevFrame * (targetFpsHz / baseFpsHz)));
    scaledFixedOccEndFrame = max(scaledFixedDevFrame, ...
        round(baseFixedOccEndFrame * (targetFpsHz / baseFpsHz)));
    scaledFixedOccEndFrame = min(scaledFramesPerTrial, scaledFixedOccEndFrame);
    cmdStim = sprintf([ ...
        'cd(''%s''); ' ...
        'targetSubjectID=%d; overwriteExisting=%d; ' ...
        'frameFrequencyOverride=%.10f; fixedDevianceFrame=%d; fixedOcclusionEndFrame=%d; ' ...
        'run(''%s'');'], ...
        practiceDirEsc, targetSubjectID, double(overwriteExisting), ...
        double(targetFpsHz), scaledFixedDevFrame, scaledFixedOccEndFrame, stimScriptEsc);
    evalin('base', cmdStim);

    if ~local_input_file_matches_fps(canonicalObservedPath, targetFpsHz)
        actualFps = local_read_input_cfg_fps(canonicalObservedPath);
        error('Generated canonical input fps mismatch. Expected %d Hz but found %.3f Hz in %s.', ...
            targetFpsHz, actualFps, canonicalObservedPath);
    end

    % Section: generate schedule file aligned to canonical input.
    % Data flow: canonical observed input + schedule script -> canonical TrialStruct.
    iSub = targetSubjectID;
    randomSeed = round(double(sharedSeed));
    numBlocks = 1;
    overwriteExisting = true;
    createScriptEsc = local_escape_single_quotes(fullfile(practiceDir, 'CreateInputFiles_v21_rescueTraject.m'));
    cmdCreate = sprintf([ ...
        'cd(''%s''); ' ...
        'iSub=%d; randomSeed=%d; numBlocks=%d; overwriteExisting=%d; ' ...
        'run(''%s'');'], ...
        practiceDirEsc, iSub, randomSeed, numBlocks, double(overwriteExisting), createScriptEsc);
    evalin('base', cmdCreate);

    if ~isfile(canonicalObservedPath) || ~isfile(canonicalPredictedPath) || ~isfile(canonicalSchedulePath)
        error('Failed to generate one or more canonical practice input files for %d Hz.', targetFpsHz);
    end

    % Cache generated canonical files under fps-tagged names.
    copyfile(canonicalObservedPath, cachedObservedPath, 'f');
    copyfile(canonicalPredictedPath, cachedPredictedPath, 'f');
    copyfile(canonicalSchedulePath, cachedSchedulePath, 'f');
    fprintf('Generated and cached practice inputs for %d Hz (%s).\n', targetFpsHz, fpsTag);
end

inputSelection = struct();
inputSelection.targetFpsHz = targetFpsHz;
inputSelection.fpsTag = fpsTag;
inputSelection.cachedObservedPath = cachedObservedPath;
inputSelection.cachedPredictedPath = cachedPredictedPath;
inputSelection.cachedSchedulePath = cachedSchedulePath;
inputSelection.canonicalObservedPath = canonicalObservedPath;
inputSelection.canonicalPredictedPath = canonicalPredictedPath;
inputSelection.canonicalSchedulePath = canonicalSchedulePath;
end

function [detectedRefreshHz, targetFpsHz, methodLabel] = local_resolve_practice_refresh_hz(overrideRefreshHz)
% LOCAL_RESOLVE_PRACTICE_REFRESH_HZ Resolve display refresh and target fps.
if ~isempty(overrideRefreshHz)
    detectedRefreshHz = double(overrideRefreshHz);
    methodLabel = 'override';
else
    [detectedRefreshHz, methodLabel] = local_detect_system_refresh_hz();
end

if ~isfinite(detectedRefreshHz) || detectedRefreshHz <= 0
    detectedRefreshHz = 60.0;
    methodLabel = 'fallback_default_60hz';
end

% Practice generation uses integer frame rates.
targetFpsHz = round(detectedRefreshHz);
targetFpsHz = max(30, min(240, targetFpsHz));
end

function [refreshHz, methodLabel] = local_detect_system_refresh_hz()
% LOCAL_DETECT_SYSTEM_REFRESH_HZ Detect nominal display refresh from system info.
% macOS: use system_profiler first (stable nominal refresh like 60/120 Hz).
if ismac
    [status, cmdOut] = system('system_profiler SPDisplaysDataType 2>/dev/null');
    if status == 0
        tokens = regexp(cmdOut, 'Refresh Rate:\s*([0-9]+(?:\.[0-9]+)?)\s*Hz', 'tokens', 'once');
        if ~isempty(tokens)
            refreshHz = str2double(tokens{1});
            if isfinite(refreshHz) && refreshHz > 0
                methodLabel = 'system_profiler';
                return;
            end
        end
    end
end

% Java fallback.
try
    ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
    devices = ge.getScreenDevices();
    if ~isempty(devices)
        mode = devices(1).getDisplayMode();
        rr = double(mode.getRefreshRate());
        if isfinite(rr) && rr > 1
            refreshHz = rr;
            methodLabel = 'java_display_mode';
            return;
        end
    end
catch
end

refreshHz = 60.0;
methodLabel = 'fallback_default_60hz';
end

function isMatch = local_input_file_matches_fps(inputPath, targetFpsHz)
% LOCAL_INPUT_FILE_MATCHES_FPS True when input MAT Cfg.fps matches target.
isMatch = false;
if ~isfile(inputPath)
    return;
end
actualFps = local_read_input_cfg_fps(inputPath);
if isfinite(actualFps) && abs(actualFps - double(targetFpsHz)) <= 0.5
    isMatch = true;
end
end

function fps = local_read_input_cfg_fps(inputPath)
% LOCAL_READ_INPUT_CFG_FPS Read Cfg.fps from one generated input file.
fps = NaN;
if ~isfile(inputPath)
    return;
end
data = load(inputPath, 'Cfg');
if ~isfield(data, 'Cfg') || ~isfield(data.Cfg, 'fps')
    return;
end
fps = double(data.Cfg.fps);
if isempty(fps) || ~isfinite(fps)
    fps = NaN;
else
    fps = fps(1);
end
end

function nTrials = local_count_valid_run_trials(trialStructPath, blockIndex, runIndex)
% LOCAL_COUNT_VALID_RUN_TRIALS Count valid source trials in one TrialOrder slice.
tsData = load(trialStructPath, 'TrialOrder');
if ~isfield(tsData, 'TrialOrder') || ndims(tsData.TrialOrder) ~= 3
    error('Invalid TrialOrder in %s', trialStructPath);
end
if blockIndex < 1 || blockIndex > size(tsData.TrialOrder, 1)
    error('Block index %d out of range.', blockIndex);
end
if runIndex < 1 || runIndex > size(tsData.TrialOrder, 3)
    error('Run index %d out of range.', runIndex);
end

runVec = squeeze(tsData.TrialOrder(blockIndex, :, runIndex));
runVec = double(runVec(:)');
validMask = isfinite(runVec) & runVec >= 1;
nTrials = sum(validMask);
if nTrials < 1
    error('Resolved empty run for block %d run %d in %s.', blockIndex, runIndex, trialStructPath);
end
end

function pct = local_trials_to_percentage(targetTrialCount, baseTrialCount)
% LOCAL_TRIALS_TO_PERCENTAGE Convert an absolute trial target to runtime percentage.
pct = 100 * (double(targetTrialCount) / double(baseTrialCount));
pct = max(0.001, min(100.0, pct));
end

function runRaw = local_execute_runtime_run( ...
        practiceDir, rawOutputDir, ...
        sharedSubId, viewingDistanceMm, runNumber, runPracticePct, ...
        useEyeTracker, doCalibration, dryRunValidate, skipSyncTests, ...
        disablePriorityBoost, attemptIndex)
% LOCAL_EXECUTE_RUNTIME_RUN Execute one runtime run via v18 script.
%
% Data flow:
%   run-specific overrides -> parent v18 runtime output MAT -> archived raw MAT.
runtimeOutputPath = fullfile(practiceDir, 'output_files', ...
    sprintf('MoveDot1_occlusion_v18_rescueTraject_SUB%02d_BLOCK%02d.mat', sharedSubId, 1));

% Move stale runtime output aside to keep this call deterministic.
if isfile(runtimeOutputPath)
    stalePath = fullfile(rawOutputDir, sprintf('stale_before_run%d_attempt%02d_%s.mat', ...
        runNumber, attemptIndex, datestr(now, 'yyyymmdd_HHMMSSFFF')));
    movefile(runtimeOutputPath, stalePath);
end

% Section: runtime overrides for a non-MEG, catch-only practice run.
% Data flow: wrapper controls -> base-workspace command -> v18 runtime.
practiceDirEsc = local_escape_single_quotes(practiceDir);
runtimeScriptEsc = local_escape_single_quotes(fullfile(practiceDir, 'MoveDot1_experiment_occlusion_v18_rescueTraject.m'));
cmdRuntime = sprintf([ ...
    'cd(''%s''); ' ...
    'iSub=%d; iBlock=1; viewingDistanceMm=%.10f; ' ...
    'practice_mode=1; practiceRunPercentageTrials=%.10f; debug_runNumber=%d; ' ...
    'dryRunValidateScheduleOnly=%d; dryRunSimulatedBreakKeys={}; ' ...
    'debugSkipSyncTests=%d; ' ...
    'disablePriorityBoost=%d; ' ...
    'trackeye=%d; ignoreEyeTracker=%d; fakeEyeTracker=0; ' ...
    'enableFixationAbort=%d; replayLimit=1; allowInfiniteReplays=0; ' ...
    'eyeTrackerDoCalibration=%d; eyeTrackerAllowBetweenBlockCalibration=0; ' ...
    'run1TransitionCountdownEnabled=0; ' ...
    'run(''%s'');'], ...
    practiceDirEsc, ...
    round(double(sharedSubId)), double(viewingDistanceMm), double(runPracticePct), round(double(runNumber)), ...
    double(logical(dryRunValidate)), ...
    double(logical(skipSyncTests)), ...
    double(logical(disablePriorityBoost)), ...
    double(logical(useEyeTracker)), double(~logical(useEyeTracker)), ...
    double(logical(useEyeTracker)), ...
    double(logical(doCalibration) && logical(useEyeTracker)), ...
    runtimeScriptEsc);
try
    evalin('base', cmdRuntime);
catch ME
    % Failsafe cleanup path to avoid leaving a full-screen PTB window open
    % when runtime throws before its own cleanup completes.
    local_force_ptb_reset();
    rethrow(ME);
end

if ~isfile(runtimeOutputPath)
    error('Runtime output not found after run %d: %s', runNumber, runtimeOutputPath);
end

rawFileName = sprintf('raw_attempt%02d_run%d_%s.mat', ...
    attemptIndex, runNumber, datestr(now, 'yyyymmdd_HHMMSSFFF'));
rawPath = fullfile(rawOutputDir, rawFileName);
movefile(runtimeOutputPath, rawPath);

loaded = load(rawPath);
runRaw = struct();
runRaw.rawPath = rawPath;
runRaw.loaded = loaded;
runRaw.runNumber = runNumber;
runRaw.practicePct = runPracticePct;
end

function local_show_run1_to_run2_message(dryRunValidate, autoContinue)
% LOCAL_SHOW_RUN1_TO_RUN2_MESSAGE Print explicit separation/instructions.
fprintf('\n--- END OF RUN 1 ---\n');
fprintf('Run 2 will now start.\n');
fprintf('Instruction reminder for run 2:\n');
fprintf('  Keep central fixation and answer the catch question quickly/accurately.\n');
fprintf('Press key/button 1 or 8 to continue to run 2.\n');
if ~logical(dryRunValidate) && ~logical(autoContinue)
    local_wait_for_1_or_8();
elseif ~logical(dryRunValidate) && logical(autoContinue)
    fprintf('Auto-continue enabled: proceeding to run 2 without key press.\n');
end
end

function attemptsOut = local_append_attempt(attemptsIn, attemptSummary)
% LOCAL_APPEND_ATTEMPT Append one attempt summary to struct array.
if isempty(attemptsIn)
    attemptsOut = attemptSummary;
else
    attemptsOut = [attemptsIn, attemptSummary]; %#ok<AGROW>
end
end

function attemptSummary = local_build_attempt_summary(attemptIndex, run1Raw, run2Raw, dryRunValidate)
% LOCAL_BUILD_ATTEMPT_SUMMARY Build run-level and catch-type summary metrics.
attemptSummary = struct();
attemptSummary.attemptIndex = attemptIndex;
attemptSummary.isDryRun = logical(dryRunValidate);
attemptSummary.run1 = local_extract_run_metrics(run1Raw);
attemptSummary.run2 = local_extract_run_metrics(run2Raw);
attemptSummary.createdAtIso8601 = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
end

function metrics = local_extract_run_metrics(runRaw)
% LOCAL_EXTRACT_RUN_METRICS Collect catch responses and summary stats from one run.
metrics = struct();
metrics.rawPath = runRaw.rawPath;
metrics.runNumber = runRaw.runNumber;
metrics.practicePercentage = runRaw.practicePct;
metrics.summaryByCatchType = struct([]);
metrics.catchTrials = struct([]);
metrics.nTrialsCompleted = NaN;

loaded = runRaw.loaded;
if ~isfield(loaded, 'runResults') || isempty(loaded.runResults)
    return;
end

rr = loaded.runResults(1);
if isfield(rr, 'nTrialsCompleted') && ~isempty(rr.nTrialsCompleted)
    metrics.nTrialsCompleted = double(rr.nTrialsCompleted);
end
if ~isfield(rr, 'output') || isempty(rr.output)
    return;
end

out = rr.output;
catchMask = false(1, numel(out));
for iTrial = 1:numel(out)
    catchType = local_field_numeric_or_nan(out(iTrial), 'catch_type_code');
    catchMask(iTrial) = isfinite(catchType) && catchType > 0;
end

catchIdx = find(catchMask);
if isempty(catchIdx)
    return;
end

catchRows = repmat(struct( ...
    'trialIndex', NaN, ...
    'catchTypeCode', NaN, ...
    'catchTypeLabel', '', ...
    'expectedCode', NaN, ...
    'responseCode', NaN, ...
    'responseCorrect', NaN, ...
    'responseRtSec', NaN, ...
    'timedOut', false), 1, numel(catchIdx));

for iRow = 1:numel(catchIdx)
    iTrial = catchIdx(iRow);
    catchTypeCode = local_field_numeric_or_nan(out(iTrial), 'catch_type_code');
    catchRows(iRow).trialIndex = iTrial;
    catchRows(iRow).catchTypeCode = catchTypeCode;
    catchRows(iRow).catchTypeLabel = local_catch_type_label(catchTypeCode);
    catchRows(iRow).expectedCode = local_field_numeric_or_nan(out(iTrial), 'catch_expected_response_code');
    catchRows(iRow).responseCode = local_field_numeric_or_nan(out(iTrial), 'catch_response_code');
    catchRows(iRow).responseCorrect = local_field_numeric_or_nan(out(iTrial), 'catch_response_correct');
    catchRows(iRow).responseRtSec = local_field_numeric_or_nan(out(iTrial), 'catch_response_rt_sec');
    catchRows(iRow).timedOut = local_field_logical(out(iTrial), 'catch_timed_out');
end

metrics.catchTrials = catchRows;
metrics.summaryByCatchType = local_summarize_by_catch_type(catchRows);
end

function summaryStruct = local_summarize_by_catch_type(catchRows)
% LOCAL_SUMMARIZE_BY_CATCH_TYPE Compute accuracy/RT summaries per catch type.
if isempty(catchRows)
    summaryStruct = struct([]);
    return;
end

typeCodes = [catchRows.catchTypeCode];
uniqueTypes = unique(typeCodes(isfinite(typeCodes) & typeCodes > 0));
summaryStruct = repmat(struct( ...
    'catchTypeCode', NaN, ...
    'catchTypeLabel', '', ...
    'nCatchTrials', 0, ...
    'nScored', 0, ...
    'nAnswered', 0, ...
    'nCorrect', 0, ...
    'nTimedOut', 0, ...
    'accuracyPct', NaN, ...
    'meanRtSec', NaN, ...
    'medianRtSec', NaN), 1, numel(uniqueTypes));

for iType = 1:numel(uniqueTypes)
    tc = uniqueTypes(iType);
    mask = (typeCodes == tc);
    rows = catchRows(mask);

    expected = [rows.expectedCode];
    response = [rows.responseCode];
    correct = [rows.responseCorrect];
    rts = [rows.responseRtSec];
    timeout = logical([rows.timedOut]);

    scoredMask = isfinite(expected) & expected > 0;
    answeredMask = scoredMask & isfinite(response) & response > 0;
    correctMask = scoredMask & isfinite(correct) & correct > 0.5;
    rtMask = answeredMask & isfinite(rts) & rts >= 0;

    nScored = sum(scoredMask);
    nAnswered = sum(answeredMask);
    nCorrect = sum(correctMask);

    summaryStruct(iType).catchTypeCode = tc;
    summaryStruct(iType).catchTypeLabel = local_catch_type_label(tc);
    summaryStruct(iType).nCatchTrials = numel(rows);
    summaryStruct(iType).nScored = nScored;
    summaryStruct(iType).nAnswered = nAnswered;
    summaryStruct(iType).nCorrect = nCorrect;
    summaryStruct(iType).nTimedOut = sum(timeout & scoredMask);

    if nScored > 0
        summaryStruct(iType).accuracyPct = 100 * (double(nCorrect) / double(nScored));
    end
    if any(rtMask)
        summaryStruct(iType).meanRtSec = mean(rts(rtMask));
        summaryStruct(iType).medianRtSec = median(rts(rtMask));
    end
end
end

function local_print_attempt_summary(attemptSummary)
% LOCAL_PRINT_ATTEMPT_SUMMARY Print compact run/catch-type summary lines.
fprintf('\nPractice attempt %d summary:\n', attemptSummary.attemptIndex);
local_print_run_summary_line(attemptSummary.run1, 'Run 1');
local_print_run_summary_line(attemptSummary.run2, 'Run 2');
end

function local_print_run_summary_line(runMetrics, runLabel)
% LOCAL_PRINT_RUN_SUMMARY_LINE Print one run summary by catch type.
if isempty(runMetrics.summaryByCatchType)
    fprintf('  %s: no catch summary available.\n', runLabel);
    return;
end
for iRow = 1:numel(runMetrics.summaryByCatchType)
    row = runMetrics.summaryByCatchType(iRow);
    fprintf(['  %s [%s]: %d/%d correct (%.1f%%), answered=%d, ' ...
        'timeouts=%d, meanRT=%.3fs, medianRT=%.3fs\n'], ...
        runLabel, row.catchTypeLabel, row.nCorrect, row.nScored, ...
        row.accuracyPct, row.nAnswered, row.nTimedOut, ...
        row.meanRtSec, row.medianRtSec);
end
end

function continuePractice = local_prompt_repeat_choice(dryRunValidate)
% LOCAL_PROMPT_REPEAT_CHOICE Ask whether to repeat another practice attempt.
if logical(dryRunValidate)
    continuePractice = false;
    return;
end

fprintf('\nRepeat practice? Press 8 to repeat, 1 to finish and save.\n');
choiceIsRepeat = local_wait_for_1_or_8();
continuePractice = logical(choiceIsRepeat);
end

function isKey8 = local_wait_for_1_or_8()
% LOCAL_WAIT_FOR_1_OR_8 Wait for keyboard/button-style key 1 or 8.
%
% Return:
%   - true  when key 8 is pressed
%   - false when key 1 is pressed
isKey8 = false;

% Fallback for environments where PTB keyboard polling is unavailable.
if exist('KbCheck', 'file') ~= 2 || exist('KbName', 'file') ~= 2
    value = input('Type 8 to continue/repeat, 1 to stop: ');
    if isfinite(value) && round(value) == 8
        isKey8 = true;
    else
        isKey8 = false;
    end
    return;
end

KbName('UnifyKeyNames');
keysFor1 = {'1!', '1', 'KP_1', 'num_1'};
keysFor8 = {'8*', '8', 'KP_8', 'num_8'};
keyCodes1 = local_keynames_to_codes(keysFor1);
keyCodes8 = local_keynames_to_codes(keysFor8);

% Some PTB builds/layouts do not expose keypad aliases (e.g., KP_1/KP_8).
% If mapping fails, fall back to a terminal prompt instead of crashing.
if isempty(keyCodes1) || isempty(keyCodes8)
    warning(['Could not resolve PTB keycodes for 1/8 on this system. ' ...
        'Falling back to terminal input for continue/repeat choice.']);
    value = input('Type 8 to continue/repeat, 1 to stop: ');
    if isfinite(value) && round(value) == 8
        isKey8 = true;
    else
        isKey8 = false;
    end
    return;
end

while true
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        if any(keyCode(keyCodes8))
            isKey8 = true;
            KbReleaseWait;
            return;
        end
        if any(keyCode(keyCodes1))
            isKey8 = false;
            KbReleaseWait;
            return;
        end
        KbReleaseWait;
    end
    WaitSecs(0.01);
end
end

function keyCodes = local_keynames_to_codes(keyNames)
% LOCAL_KEYNAMES_TO_CODES Convert key-name cellstr to valid key codes.
%
% Note:
%   KbName throws on unknown key labels in some PTB builds. We treat those
%   labels as unavailable and continue scanning alternatives.
keyCodes = [];
for iKey = 1:numel(keyNames)
    try
        code = KbName(keyNames{iKey});
    catch
        continue;
    end
    if isempty(code)
        continue;
    end
    if ischar(code) || isstring(code)
        continue;
    end
    if numel(code) > 1
        code = code(1);
    end
    if isnumeric(code) && isfinite(code)
        keyCodes(end + 1) = round(double(code)); %#ok<AGROW>
    end
end
keyCodes = unique(keyCodes);
end

function value = local_field_numeric_or_nan(s, fieldName)
% LOCAL_FIELD_NUMERIC_OR_NAN Read one numeric scalar field or NaN fallback.
value = NaN;
if ~isstruct(s) || ~isfield(s, fieldName)
    return;
end
candidate = s.(fieldName);
if isempty(candidate) || ~(isnumeric(candidate) || islogical(candidate))
    return;
end
candidate = double(candidate(:));
if isempty(candidate)
    return;
end
value = candidate(1);
if ~isfinite(value)
    value = NaN;
end
end

function value = local_field_logical(s, fieldName)
% LOCAL_FIELD_LOGICAL Read one logical-like scalar field with false fallback.
value = false;
if ~isstruct(s) || ~isfield(s, fieldName)
    return;
end
candidate = s.(fieldName);
if isempty(candidate)
    return;
end
if islogical(candidate)
    value = logical(candidate(1));
elseif isnumeric(candidate)
    value = isfinite(candidate(1)) && candidate(1) ~= 0;
end
end

function label = local_catch_type_label(typeCode)
% LOCAL_CATCH_TYPE_LABEL Human-readable label for catch type code.
switch round(typeCode)
    case 1
        label = 'type1_disappear_reappear';
    case 2
        label = 'type2_occlusion_question';
    otherwise
        label = 'unknown_or_none';
end
end

function escaped = local_escape_single_quotes(rawPath)
% LOCAL_ESCAPE_SINGLE_QUOTES Escape single quotes for evalin command strings.
escaped = strrep(char(rawPath), '''', '''''');
end

function local_run_strict_safety_preflight( ...
        strictSafeMode, dryRunValidate, skipSyncTests, disablePriorityBoost)
% LOCAL_RUN_STRICT_SAFETY_PREFLIGHT Block unsafe non-dry runs early.
%
% Safety policy:
%   - Applies only when strictSafeMode=true and dryRunValidate=false.
%   - On macOS, strict mode blocks non-dry PTB fullscreen runs entirely to
%     avoid potential black-screen lockups on unsupported PTB setups.
%   - Requires skipSyncTests=1 and disablePriorityBoost=1.
%   - Probes PTB priority MEX availability before opening any PTB window.
if ~logical(strictSafeMode) || logical(dryRunValidate)
    return;
end

if ismac
    error(sprintf(['Strict safe mode blocked launch before PTB window open.\n' ...
        'Non-dry fullscreen PTB runs are disabled in strict mode on macOS for safety.\n' ...
        'Use practiceDryRunValidateScheduleOnly=true for safe testing, or set\n' ...
        'practiceStrictSafeMode=false only if you accept risk.']));
end

if ~logical(skipSyncTests)
    error(['Strict safe mode blocked launch: practiceSkipSyncTests must be true ' ...
        'for non-dry runs on this setup.']);
end
if ~logical(disablePriorityBoost)
    error(['Strict safe mode blocked launch: practiceDisablePriorityBoost must be true ' ...
        'to avoid Priority/MEX crashes.']);
end

[probeOk, probeId, probeMessage] = local_probe_ptb_priority_mex();
if ~probeOk
    error(['Strict safe mode blocked launch before PTB window open.\n' ...
        'Priority/MEX preflight failed (%s): %s\n' ...
        'Use dry-run mode, or disable strict safe mode only if you accept risk.\n' ...
        'Suggested safe test override: practiceDryRunValidateScheduleOnly=true'], ...
        probeId, probeMessage);
end
end

function [ok, errId, errMessage] = local_probe_ptb_priority_mex()
% LOCAL_PROBE_PTB_PRIORITY_MEX Probe PTB Priority path without opening a window.
ok = true;
errId = '';
errMessage = '';

if exist('Priority', 'file') ~= 2
    ok = false;
    errId = 'Practice:PriorityFunctionMissing';
    errMessage = 'Psychtoolbox Priority function is not available on path.';
    return;
end

% Probe the same code path family used by Priority(MaxPriority(window)):
% a temporary non-zero priority request followed by reset.
try
    Priority(1);
    Priority(0);
catch ME
    ok = false;
    errId = char(ME.identifier);
    errMessage = char(ME.message);
    try
        Priority(0);
    catch
    end
end
end

function local_force_ptb_reset()
% LOCAL_FORCE_PTB_RESET Best-effort PTB cleanup on wrapper-level failures.
try
    evalin('base', 'try, Priority(0); catch, end;');
catch
end
try
    evalin('base', 'try, ShowCursor; catch, end;');
catch
end
try
    evalin('base', 'try, sca; catch, end;');
catch
end
end

function value = local_get_base_or_default(varName, defaultValue)
% LOCAL_GET_BASE_OR_DEFAULT Read optional override from base workspace.
value = defaultValue;
if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
    value = evalin('base', varName);
end
end
