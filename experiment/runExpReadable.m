% runExpReadable.m
%
% Purpose:
%   Minimal readable launcher that executes the full real experiment stack
%   in explicit order:
%     1) runtime config preset load
%     2) stimuli_generation_V30_rescueTraject.m
%     3) CreateInputFiles_v23_rescueTraject.m
%     4) MoveDot1_experiment_occlusion_v20_rescueTraject.m
%
% Why this launcher exists:
%   - Keep one place where subject/block and overwrite policy are edited.
%   - Keep subject identity consistent across generation + runtime.
%   - Avoid interactive popups during prepared runs.
%   - Survive script-level clearvars inside called scripts.
%
% Usage (interactive from experiment/):
%   run('runExpReadable.m')
%
% Usage (non-interactive shell, MATLAB R2020a on Apple Silicon):
%   arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    run('runExpReadable.m');"
%
% Inputs:
%   - lib/Config_runtime_v20_rescueTraject.m
%   - lib/Config_stimuli_generation_V30_rescueTraject.m
%   - lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m
%
% Outputs:
%   - input_files/MovDot_SubXX_V30_rescueTraject.mat
%   - input_files/SubXX_TrialStruct_v23_rescueTraject.mat
%   - output_files/MoveDot1_occlusion_v20_rescueTraject_SUBXX_BLOCKYY.mat
%
% Key assumptions:
%   - Called scripts use relative paths from experiment/.
%   - Called scripts may clear workspace variables (clearvars -except ...).

%% Resolve launcher location and pin current working directory
% Data flow: script location -> cwd for all relative input/output paths.
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end

%% Launcher targets and overwrite policy
% Data flow: launcher edits -> launcherCfg checkpoint -> restored workspace.
targetSub = 95;
targetBlock = 1;
randomSeed = targetSub;
overwriteStimuliInput = true;
overwriteScheduleInput = true;
fixation_window_allowed = 3; % in visual degree
allowInteractivePrompts = false;

launcherCfg = struct( ...
    'targetSub', targetSub, ...
    'targetBlock', targetBlock, ...
    'randomSeed', randomSeed, ...
    'overwriteStimuliInput', logical(overwriteStimuliInput), ...
    'overwriteScheduleInput', logical(overwriteScheduleInput), ...
    'fixationWindowDeg', double(fixation_window_allowed), ...
    'allowInteractivePrompts', logical(allowInteractivePrompts));

launcherCfgPath = fullfile(tempdir, ...
    sprintf('runExpReadable_cfg_pid%d.mat', feature('getpid')));
save(launcherCfgPath, 'launcherCfg');

%% Step 1: Stimuli generation (V30)
% Data flow: targetSub + overwrite flag -> MovDot_SubXX_V30_rescueTraject.mat
run('lib/Config_runtime_v20_rescueTraject.m');
addpath('lib');
cfgStim = Config_stimuli_generation_V30_rescueTraject;  % explicit dependency touchpoint
iSub = launcherCfg.targetSub;
targetSubjectID = iSub;  % consumed by called script workspace
overwriteExisting = launcherCfg.overwriteStimuliInput; %#ok<NASGU>  % consumed by called script workspace
allowInteractivePrompts = launcherCfg.allowInteractivePrompts; %#ok<NASGU>  % consumed by called script workspace
stimObservedFile = fullfile(cfgStim.inputDirectory, ...
    sprintf(cfgStim.observedFilePattern, launcherCfg.targetSub));
stimFilesExist = isfile(stimObservedFile);
if stimFilesExist && ~launcherCfg.overwriteStimuliInput
    fprintf(['Step 1 skipped: existing stimuli files found for Sub%02d and ' ...
        'overwriteStimuliInput=false.\n'], launcherCfg.targetSub);
else
    run('stimuli_generation_V30_rescueTraject.m');
end

% Stimulus generation script clears most workspace variables; restore launcher state.
if ~exist('launcherCfgPath', 'var') || isempty(launcherCfgPath)
    launcherCfgPath = fullfile(tempdir, ...
        sprintf('runExpReadable_cfg_pid%d.mat', feature('getpid')));
end
if exist(launcherCfgPath, 'file') ~= 2
    error('Launcher checkpoint file not found after stimuli generation: %s', launcherCfgPath);
end
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;

%% Step 2: Schedule/input creation (v23)
% Data flow: targetSub + randomSeed + overwrite flag -> SubXX_TrialStruct_v23*.mat
cfgInput = Config_schedule_CreateInputV23_MoveDotV20_rescueTraject;
iSub = launcherCfg.targetSub; %#ok<NASGU>  % consumed by called script workspace
randomSeed = launcherCfg.randomSeed;  % consumed by called script workspace
numBlocks = cfgInput.numBlocks;  % consumed by called script workspace
overwriteExisting = launcherCfg.overwriteScheduleInput;  % consumed by called script workspace
allowInteractivePrompts = launcherCfg.allowInteractivePrompts; %#ok<NASGU>  % consumed by called script workspace
scheduleFile = fullfile('.', 'input_files', ...
    sprintf(cfgInput.outputFilePattern, launcherCfg.targetSub));
if isfile(scheduleFile) && ~launcherCfg.overwriteScheduleInput
    fprintf(['Step 2 skipped: existing schedule file found for Sub%02d and ' ...
        'overwriteScheduleInput=false.\n'], launcherCfg.targetSub);
else
    run('CreateInputFiles_v23_rescueTraject.m');
end

% Input creation script also clears workspace; restore launcher state again.
if ~exist('launcherCfgPath', 'var') || isempty(launcherCfgPath)
    launcherCfgPath = fullfile(tempdir, ...
        sprintf('runExpReadable_cfg_pid%d.mat', feature('getpid')));
end
if exist(launcherCfgPath, 'file') ~= 2
    error('Launcher checkpoint file not found after input creation: %s', launcherCfgPath);
end
cfgData = load(launcherCfgPath, 'launcherCfg');
launcherCfg = cfgData.launcherCfg;

%% Step 2b: Verify generated files before runtime launch
% Data flow: targetSub + config naming patterns -> strict existence checks.
scheduleCfg = Config_schedule_CreateInputV23_MoveDotV20_rescueTraject;
expectedInputFile = fullfile('.', 'input_files', ...
    sprintf(scheduleCfg.inputFilePattern, launcherCfg.targetSub));
expectedTrialStructFile = fullfile('.', 'input_files', ...
    sprintf(scheduleCfg.outputFilePattern, launcherCfg.targetSub));
if exist(expectedInputFile, 'file') ~= 2
    error(['Expected generated input file not found after step 1: %s\n' ...
        'Current directory: %s'], expectedInputFile, pwd);
end
if exist(expectedTrialStructFile, 'file') ~= 2
    error(['Expected generated schedule file not found after step 2: %s\n' ...
        'Current directory: %s\n' ...
        'Check that CreateInputFiles_v23_rescueTraject.m ran for the same subject.'], ...
        expectedTrialStructFile, pwd);
end

%% Step 3: Runtime launch (MoveDot1 experiment)
% Data flow: targetSub/targetBlock + runtime overrides -> live block execution.
run('lib/Config_runtime_v20_rescueTraject.m');
iSub = launcherCfg.targetSub;
iBlock = launcherCfg.targetBlock;
fixWindowDeg = launcherCfg.fixationWindowDeg;  % consumed by called script workspace
allowInteractivePrompts = launcherCfg.allowInteractivePrompts;  % consumed by called script workspace
MEG = true;
trackeye = true;
ignoreEyeTracker = false;
fakeEyeTracker = false;
eyeTrackerDoCalibration = true;
eyeTrackerAllowBetweenBlockCalibration = true;
dryRunValidateScheduleOnly = false;
debugSkipSyncTests = 0;
run('MoveDot1_experiment_occlusion_v20_rescueTraject.m');

%% Cleanup
% Data flow: best-effort temp checkpoint cleanup at normal launcher exit.
if exist(launcherCfgPath, 'file') == 2
    delete(launcherCfgPath);
end
