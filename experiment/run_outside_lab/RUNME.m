%% RUNME.m
% Script: RUNME.m
%
% Purpose:
%   One-command launcher for sharing this experiment outside the lab.
%   The script runs the full pipeline in this folder copy only:
%     1) Generate one-dot trajectories (V30) at the detected display refresh.
%     2) Build schedule/catch-plan inputs (v23).
%     3) Launch runtime presentation (v20) in laptop-safe profile.
%
% Why t his launcher exists:
%   - Lab members can try a local test with only MATLAB + Psychtoolbox.
%   - No edits are required in the original lab scripts.
%   - Refresh-rate-dependent generation is handled automatically per
%   machine. Unless impossible to detect, then a user input is required.
%
% Usage example (interactive):
%   1) Open MATLAB.
%   2) cd into this folder (run_outside_lab).
%   3) run('RUNME.m')
%
% Usage example (non-interactive shell):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/absolute/path/to/run_outside_lab'); ...
%    viewingDistanceMm=600; practice_mode=1; ...
%    run('RUNME.m');"
%
% Optional workspace defaults (set before run):
%   - subjectStartID           : first auto-assignable subject id (default 50)
%   - iBlock                   : start block index (default 1)
%   - numBlocks                : number of blocks to schedule (default 1)
%   - viewingDistanceMm        : viewing distance in mm (default 600)
%   - practice_mode            : 0 full run, 1 shortened practice flow
%   - debugSkipSyncTests       : 1 enables PTB SkipSyncTests (default 1)
%   - generationFrameFrequencyHz : manual fps override for generation
%   - dryRunValidateScheduleOnly : true validates schedule without PTB display
%   - allowInteractivePrompts  : true allows dialogs when values are missing
%
% Subject numbering policy:
%   RUNME auto-selects the first subject id with no existing files,
%   scanning from subjectStartID upward (50, 51, 52, ... by default).
%
% Startup dialog policy:
%   RUNME always shows a setup popup before starting. Fields are pre-filled
%   from workspace values when present.
%
% Outputs:
%   - input_files/MovDot_SubXX_V30_rescueTraject*.mat
%   - input_files/SubXX_TrialStruct_v23_rescueTraject.mat
%   - output_files/MoveDot1_occlusion_v20_rescueTraject_SUBXX_BLOCKYY.mat
%
% Key assumptions:
%   - MATLAB + Psychtoolbox are available.
%   - This script is run from this copied folder (it self-corrects via scriptDir).
%   - The copied v20 runtime is used in non-MEG, non-EyeLink mode.

%% Resolve launcher location and ensure local copy takes precedence
% Data flow: launcher path -> deterministic working directory + local lib path.
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end
addpath(fullfile(scriptDir, 'lib'), '-begin');

if exist('input_files', 'dir') ~= 7
    mkdir('input_files');
end
if exist('output_files', 'dir') ~= 7
    mkdir('output_files');
end

% Reset any open PTB windows from previous interrupted runs.
sca;

%% Resolve run parameters (always-show setup dialog with workspace-prefill)
% Data flow: caller defaults + setup dialog -> validated launcher config.
if ~exist('allowInteractivePrompts', 'var') || isempty(allowInteractivePrompts)
    allowInteractivePrompts = true;
end
allowInteractivePrompts = logical(allowInteractivePrompts);
% Always show the startup setup popup, even if a previous run left
% allowInteractivePrompts=false in the base workspace.
allowInteractivePrompts = true;

if ~exist('iBlock', 'var') || isempty(iBlock)
    iBlock = 1;
end
if ~exist('numBlocks', 'var') || isempty(numBlocks)
    numBlocks = 1;
end
if ~exist('viewingDistanceMm', 'var') || isempty(viewingDistanceMm)
    viewingDistanceMm = NaN;
end
if ~exist('practice_mode', 'var') || isempty(practice_mode)
    practice_mode = 1;
end
if ~exist('debugSkipSyncTests', 'var') || isempty(debugSkipSyncTests)
    debugSkipSyncTests = 1;
end
if ~exist('dryRunValidateScheduleOnly', 'var') || isempty(dryRunValidateScheduleOnly)
    dryRunValidateScheduleOnly = false;
end
if ~exist('generationFrameFrequencyHz', 'var')
    generationFrameFrequencyHz = [];
end
if ~exist('subjectStartID', 'var') || isempty(subjectStartID)
    subjectStartID = 50;
end

defaultStartID = local_default_numeric_string(subjectStartID, 50);
defaultDistance = local_default_numeric_string(viewingDistanceMm, 600);
defaultPractice = local_default_numeric_string(practice_mode, 1);
defaultBlocks = local_default_numeric_string(numBlocks, 1);
defaultBlockStart = local_default_numeric_string(iBlock, 1);

prompt = {
    'Start subject ID (auto-select first free >= this value):', ...
    'Viewing Distance (mm):', ...
    'Practice mode (0 full / 1 practice):', ...
    'Number of blocks to schedule:', ...
    'Start block index:'};
dlgTitle = 'RUNME setup (outside-lab copy)';
dims = [1 80];
answer = inputdlg(prompt, dlgTitle, dims, ...
    {defaultStartID, defaultDistance, defaultPractice, defaultBlocks, defaultBlockStart});

if isempty(answer)
    fprintf('RUNME canceled by user.\n');
    return;
end

subjectStartID = str2double(answer{1});
viewingDistanceMm = str2double(answer{2});
practice_mode = str2double(answer{3});
numBlocks = str2double(answer{4});
iBlock = str2double(answer{5});

% Final validation and canonical numeric casting.
iBlock = round(double(iBlock));
numBlocks = round(double(numBlocks));
viewingDistanceMm = double(viewingDistanceMm);
practice_mode = round(double(practice_mode));
subjectStartID = round(double(subjectStartID));
debugSkipSyncTests = logical(round(double(debugSkipSyncTests)) ~= 0);
dryRunValidateScheduleOnly = logical(dryRunValidateScheduleOnly);

if ~isfinite(iBlock) || iBlock < 1
    error('iBlock must be a finite integer >= 1.');
end
if ~isfinite(numBlocks) || numBlocks < 1
    error('numBlocks must be a finite integer >= 1.');
end
if ~isfinite(viewingDistanceMm) || viewingDistanceMm <= 0
    error('viewingDistanceMm must be a finite scalar > 0.');
end
if ~ismember(practice_mode, [0, 1])
    error('practice_mode must be 0 or 1.');
end
if ~isfinite(subjectStartID) || subjectStartID < 0
    error('subjectStartID must be a finite integer >= 0.');
end

%% Auto-assign first free subject id and schedule seed
% Data flow: subjectStartID + on-disk artifacts -> first available subject id -> run seed.
%
% Policy:
%   Starting from subjectStartID (default 50), pick the first subject id
%   for which no run_outside_lab input/output artifacts exist.
iSub = local_pick_next_available_subject(subjectStartID);
randomSeed = round(double(iSub) * 100 + 1);
fprintf('RUNME: auto-selected subject %d (search start %d).\n', iSub, subjectStartID);

%% Resolve generation refresh rate with PTB auto-detection + fallback prompt
% Data flow: display query -> generation fps override -> stimuli generation timing.
resolvedFrameRateHz = NaN;
if ~isempty(generationFrameFrequencyHz)
    resolvedFrameRateHz = double(generationFrameFrequencyHz);
end

if ~(isfinite(resolvedFrameRateHz) && resolvedFrameRateHz > 0)
    resolvedFrameRateHz = local_detect_refresh_rate_hz(debugSkipSyncTests);
end

if ~(isfinite(resolvedFrameRateHz) && resolvedFrameRateHz > 0)
    if ~allowInteractivePrompts
        error(['Could not auto-detect display refresh rate and prompts are disabled. ' ...
            'Set generationFrameFrequencyHz manually.']);
    end
    resolvedFrameRateHz = local_prompt_refresh_rate_hz(60);
end

fprintf('RUNME: using generation refresh %.3f Hz\n', resolvedFrameRateHz);

%% Snapshot launcher state across script-level clearvars calls
% Data flow: launcher settings -> temporary MAT file -> reloaded per step.
launcherCfg = struct();
launcherCfg.iSub = iSub;
launcherCfg.iBlock = iBlock;
launcherCfg.numBlocks = numBlocks;
launcherCfg.viewingDistanceMm = viewingDistanceMm;
launcherCfg.practice_mode = practice_mode;
launcherCfg.randomSeed = randomSeed;
launcherCfg.debugSkipSyncTests = debugSkipSyncTests;
launcherCfg.dryRunValidateScheduleOnly = dryRunValidateScheduleOnly;
launcherCfg.allowInteractivePrompts = false; % force deterministic non-dialog child scripts.
launcherCfg.generationFrameFrequencyHz = resolvedFrameRateHz;

launcherCfgPath = fullfile(tempdir, sprintf('RUNME_outside_lab_cfg_pid%d.mat', feature('getpid')));
save(launcherCfgPath, 'launcherCfg');

%% Step 1: Generate stimuli at resolved refresh rate
% Data flow: launcherCfg + refresh rate -> input_files/MovDot_SubXX*.mat
% Ensure copied classdefs are reloaded from this folder path.
clear Config_stimuli_generation_V30_rescueTraject Config_schedule_CreateInputV23_MoveDotV20_rescueTraject;
launcherCfgPath = fullfile(tempdir, sprintf('RUNME_outside_lab_cfg_pid%d.mat', feature('getpid')));
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;

targetSubjectID = launcherCfg.iSub;
iSub = launcherCfg.iSub;
overwriteExisting = true;
allowInteractivePrompts = launcherCfg.allowInteractivePrompts;
generationFrameFrequencyHz = launcherCfg.generationFrameFrequencyHz;
autoScaleFixedFramesFromConfig = true;
showPathPlots = false;
run('stimuli_generation_V30_rescueTraject.m');

%% Step 2: Build schedule and catch plan
% Data flow: launcherCfg schedule parameters -> input_files/SubXX_TrialStruct*.mat
launcherCfgPath = fullfile(tempdir, sprintf('RUNME_outside_lab_cfg_pid%d.mat', feature('getpid')));
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;

iSub = launcherCfg.iSub;
randomSeed = launcherCfg.randomSeed;
numBlocks = launcherCfg.numBlocks;
overwriteExisting = true;
allowInteractivePrompts = launcherCfg.allowInteractivePrompts;
run('CreateInputFiles_v23_rescueTraject.m');

%% Step 3: Launch runtime in laptop-safe profile
% Data flow: launcherCfg runtime settings + generated files -> experiment output MAT.
launcherCfgPath = fullfile(tempdir, sprintf('RUNME_outside_lab_cfg_pid%d.mat', feature('getpid')));
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;

iSub = launcherCfg.iSub;
iBlock = launcherCfg.iBlock;
viewingDistanceMm = launcherCfg.viewingDistanceMm;
practice_mode = launcherCfg.practice_mode;
allowInteractivePrompts = launcherCfg.allowInteractivePrompts;

% Outside-lab defaults: no MEG and no EyeLink dependencies.
MEG = false;
trackeye = false;
ignoreEyeTracker = true;
fakeEyeTracker = false;
eyeTrackerDoCalibration = false;
eyeTrackerAllowBetweenBlockCalibration = false;
eyeTrackerSendTriggerMessages = false;
eyeTrackerImageTransferEnabled = false;

% Keep runtime robust on unknown consumer displays.
debugSkipSyncTests = launcherCfg.debugSkipSyncTests;
dryRunValidateScheduleOnly = launcherCfg.dryRunValidateScheduleOnly;

% Keep ITI and catch-screen geometry aligned with current test profile.
ITIRangeSec = [1.000, 1.500];
itiVisualMode = 'fixation_short';
itiGreyColor = [70, 70, 70];
itiVisualSeparationDur = [0.300, 0.700];
miniBlinkBreakEnabled = true;
miniBlinkBreakEveryNTrials = 10;
miniBlinkBreakDurationSec = 5.0;
catchQuestionTopYFrac = 0.42;
catchQuestionBottomYFrac = 0.58;
catchBottomLeftXFrac = 0.38;
catchBottomRightXFrac = 0.62;
run1TransitionCountdownEnabled = true;
debugTriggerMode = 0;
question_text = 'Has the dot changed its course?';

displayRunMsg = sprintf(['RUNME launch: Sub%02d Block%02d | viewingDistanceMm=%.1f | ' ...
    'practice=%d | seed=%d | fps=%.3f Hz\n'], ...
    iSub, iBlock, viewingDistanceMm, practice_mode, launcherCfg.randomSeed, launcherCfg.generationFrameFrequencyHz);
fprintf('%s', displayRunMsg);

run('MoveDot1_experiment_occlusion_v20_rescueTraject.m');

%% Cleanup
% Data flow: temporary launcher state file -> best-effort delete.
launcherCfgPath = fullfile(tempdir, sprintf('RUNME_outside_lab_cfg_pid%d.mat', feature('getpid')));
if exist(launcherCfgPath, 'file') == 2
    delete(launcherCfgPath);
end

%% Local helper functions
function valStr = local_default_numeric_string(rawVal, fallbackVal)
% LOCAL_DEFAULT_NUMERIC_STRING Build a display-safe default numeric string.
if nargin < 2
    fallbackVal = NaN;
end
if isempty(rawVal)
    rawVal = NaN;
end
rawVal = double(rawVal);
if isscalar(rawVal) && isfinite(rawVal)
    valStr = sprintf('%g', rawVal);
elseif isscalar(fallbackVal) && isfinite(double(fallbackVal))
    valStr = sprintf('%g', double(fallbackVal));
else
    valStr = '';
end
end

function resolvedHz = local_prompt_refresh_rate_hz(defaultHz)
% LOCAL_PROMPT_REFRESH_RATE_HZ Ask user for refresh-rate fallback.
if nargin < 1 || isempty(defaultHz)
    defaultHz = 60;
end
answer = inputdlg({'Enter display refresh rate (Hz):'}, ...
    'Refresh-rate fallback', [1 70], {sprintf('%.3f', double(defaultHz))});
if isempty(answer)
    error('Refresh-rate prompt canceled by user.');
end
resolvedHz = str2double(answer{1});
if ~isscalar(resolvedHz) || ~isfinite(resolvedHz) || resolvedHz <= 0
    error('Refresh rate must be a positive numeric value in Hz.');
end
end

function detectedHz = local_detect_refresh_rate_hz(skipSyncTests)
% LOCAL_DETECT_REFRESH_RATE_HZ Detect active display refresh via PTB.
%
% Detection strategy:
%   1) Read nominal display rate from Screen('NominalFrameRate').
%   2) Open a temporary PTB window and refine with GetFlipInterval.
%   3) Keep best finite estimate in [20, 500] Hz.
detectedHz = NaN;
window = [];
oldSkip = [];

try
    screens = Screen('Screens');
    if isempty(screens)
        return;
    end
    screenNumber = max(screens);

    nominalHz = double(Screen('NominalFrameRate', screenNumber));
    if isfinite(nominalHz) && nominalHz >= 20 && nominalHz <= 500
        detectedHz = nominalHz;
    end

    oldSkip = Screen('Preference', 'SkipSyncTests');
    Screen('Preference', 'SkipSyncTests', double(logical(skipSyncTests)));

    [window, ~] = Screen('OpenWindow', screenNumber, 0);
    ifi = double(Screen('GetFlipInterval', window));
    measuredHz = 1 / ifi;
    if isfinite(measuredHz) && measuredHz >= 20 && measuredHz <= 500
        detectedHz = measuredHz;
    end
catch detectErr
    fprintf('RUNME warning: refresh detection failed (%s).\n', detectErr.message);
end

if ~isempty(window)
    Screen('Close', window);
end
if ~isempty(oldSkip)
    Screen('Preference', 'SkipSyncTests', oldSkip);
end
sca;
end

function subjectID = local_pick_next_available_subject(startID)
% LOCAL_PICK_NEXT_AVAILABLE_SUBJECT Return first subject id with no local artifacts.
%
% Data flow:
%   start subject id -> per-id artifact existence checks -> first free id.
startID = round(double(startID));
if ~isscalar(startID) || ~isfinite(startID) || startID < 0
    error('startID must be a finite scalar >= 0.');
end

maxSearchID = startID + 100000;
for candidateID = startID:maxSearchID
    if ~local_subject_has_existing_artifacts(candidateID)
        subjectID = candidateID;
        return;
    end
end

error('Could not find an available subject id in range [%d, %d].', startID, maxSearchID);
end

function occupied = local_subject_has_existing_artifacts(iSub)
% LOCAL_SUBJECT_HAS_EXISTING_ARTIFACTS Check whether this subject id is already used.
%
% Data flow:
%   subject id -> expected input/output filenames -> occupied flag.
iSub = round(double(iSub));
occupied = false;

observedFile = fullfile('input_files', sprintf('MovDot_Sub%02d_V30_rescueTraject.mat', iSub));
trialStructFile = fullfile('input_files', sprintf('Sub%02d_TrialStruct_v23_rescueTraject.mat', iSub));

if exist(observedFile, 'file') == 2 || exist(trialStructFile, 'file') == 2
    occupied = true;
    return;
end

runtimePattern = fullfile('output_files', ...
    sprintf('MoveDot1_occlusion_v20_rescueTraject_SUB%02d_BLOCK*.mat', iSub));
runtimeMatches = dir(runtimePattern);
debugCsvPattern = fullfile('output_files', ...
    sprintf('debug_actual_triggers_occlusion_v20_rescueTraject_sub%02d_block*_run*.csv', iSub));
debugCsvMatches = dir(debugCsvPattern);
if ~isempty(runtimeMatches) || ~isempty(debugCsvMatches)
    occupied = true;
end
end
