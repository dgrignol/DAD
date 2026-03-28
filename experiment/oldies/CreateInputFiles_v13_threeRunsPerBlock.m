% CreateInputFiles_v13_threeRunsPerBlock.m
%
% Purpose:
%   Build an additive, block-aware TrialStruct artifact for the one-dot
%   occlusion paradigm, with fixed "3 runs per block" scheduling:
%     - Run 1: all always_visible trials (shuffled)
%     - Run 2: random half of occluded_nondeviant + random half of
%              occluded_deviant (shuffled together)
%     - Run 3: remaining half of occluded_nondeviant + remaining half of
%              occluded_deviant (shuffled together)
%
% Why this v13 file exists:
%   - Keep existing scheduling scripts untouched.
%   - Provide a versioned artifact and versioned logic for the requested
%     multi-block, three-run schedule.
%
% Usage example (interactive):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment');
%   CreateInputFiles_v13_threeRunsPerBlock
%
% Usage example (non-interactive):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    iSub=66; randomSeed=6601; numBlocks=2; run('CreateInputFiles_v13_threeRunsPerBlock.m');"
%
% Inputs:
%   - input_files/MovDot_SubXX.mat with:
%       * xySeqs
%       * Cfg
%   - Config_occlusion_schedule_v1 constants:
%       * numBlocks
%       * runsPerBlock (must be 3 for this script)
%       * run2FractionPerOccluded (must be 0.5 for this fixed design)
%
% Workspace overrides (optional):
%   - iSub
%   - randomSeed
%   - numBlocks
%   - overwriteExisting
%
% Outputs:
%   - input_files/SubXX_TrialStruct_v13_threeRunsPerBlock.mat with:
%       * TrialStruct  (nRuns x nTrialsFlattened)
%       * Catch        (no-catch defaults)
%       * TrialOrder   (numBlocks x runLength x nRuns)
%       * BlockOrder   (1 x numBlocks)
%       * Schedule     (explicit schedule metadata)
%
% Key assumptions:
%   - Input labels include exactly:
%       always_visible, occluded_nondeviant, occluded_deviant
%   - Occluded condition counts are equal and even.
%   - For equal run length:
%       count(always_visible) == count(occluded_nondeviant) ==
%       count(occluded_deviant)

%% Setup and input collection
% Data flow: dialog/workspace overrides + schedule config -> file names and controls.
clearvars -except iSub randomSeed numBlocks overwriteExisting;
clc;
addpath('lib/');

scheduleCfg = Config_occlusion_schedule_v1;
if scheduleCfg.runsPerBlock ~= 3
    error('Config_occlusion_schedule_v1.runsPerBlock must be 3 for this script.');
end
if abs(scheduleCfg.run2FractionPerOccluded - 0.5) > eps
    error('Config_occlusion_schedule_v1.run2FractionPerOccluded must be 0.5.');
end

if ~(exist('iSub', 'var') && exist('randomSeed', 'var'))
    prompt = {'Subject Number:', 'Random seed:'};
    dlgtitle = 'CreateInputFiles v13 (three runs per block)';
    dims = [1 60];
    definput = {'66', '6601'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(answer)
        disp('Canceled by user.');
        return;
    end
    iSub = str2double(answer{1});
    randomSeed = str2double(answer{2});
end

if ~exist('numBlocks', 'var') || isempty(numBlocks)
    numBlocks = double(scheduleCfg.numBlocks);
end

if any(isnan([iSub, randomSeed, numBlocks]))
    error('Subject, random seed, and numBlocks must be numeric.');
end

iSub = round(iSub);
randomSeed = round(randomSeed);
numBlocks = round(numBlocks);
if iSub < 0 || numBlocks < 1
    error('Subject must be >= 0 and numBlocks must be >= 1.');
end

inputDir = fullfile('.', 'input_files');
inputFile = fullfile(inputDir, sprintf(scheduleCfg.inputFilePattern, iSub));
outputFile = fullfile(inputDir, sprintf(scheduleCfg.outputFilePattern, iSub));

if ~isfile(inputFile)
    error('Input file not found: %s', inputFile);
end

if exist(outputFile, 'file') == 2
    doOverwrite = false;
    if exist('overwriteExisting', 'var') && ~isempty(overwriteExisting)
        doOverwrite = logical(overwriteExisting);
    end
    if ~doOverwrite
        error('Output already exists; refusing overwrite: %s', outputFile);
    end
end

%% Load source trajectories and validate required labels
% Data flow: MovDot_SubXX.mat -> xySeqs labels -> condition index pools.
Dat = load(inputFile, 'xySeqs', 'Cfg');
if ~isfield(Dat, 'xySeqs') || isempty(Dat.xySeqs)
    error('Input file must contain non-empty xySeqs: %s', inputFile);
end
if ~isfield(Dat, 'Cfg')
    error('Input file must contain Cfg: %s', inputFile);
end

xyAll = Dat.xySeqs(:);
if ~isfield(xyAll, 'xy')
    error('xySeqs entries must contain .xy');
end
if ~isfield(xyAll, 'condition_label')
    error('xySeqs entries must contain .condition_label');
end

labels = arrayfun(@(s) string(s.condition_label), xyAll, 'UniformOutput', true);
idxAlways = find(labels == "always_visible");
idxNondev = find(labels == "occluded_nondeviant");
idxDev = find(labels == "occluded_deviant");

if isempty(idxAlways) || isempty(idxNondev) || isempty(idxDev)
    error(['Expected condition labels always_visible, occluded_nondeviant, and ' ...
        'occluded_deviant in %s.'], inputFile);
end

nAlways = numel(idxAlways);
nNondev = numel(idxNondev);
nDev = numel(idxDev);

if nNondev ~= nDev
    error('occluded_nondeviant (%d) and occluded_deviant (%d) counts must match.', nNondev, nDev);
end

if mod(nNondev, 2) ~= 0
    error('Occluded condition count must be even. Found %d.', nNondev);
end

if logical(scheduleCfg.enforceEqualRunLengths) && (nAlways ~= nNondev)
    error('For equal run lengths, always_visible (%d) must match each occluded count (%d).', ...
        nAlways, nNondev);
end

nOccPerCondition = nNondev;
nOccHalf = nOccPerCondition / 2;
runLength = nAlways;
nRuns = scheduleCfg.runsPerBlock;

if runLength ~= 2 * nOccHalf
    error('Run length mismatch: run1 length=%d but run2/run3 length=%d. Use matched counts or disable enforceEqualRunLengths in config.', ...
        runLength, 2 * nOccHalf);
end

%% Build block-aware TrialOrder with fixed three-run design
% Data flow: condition pools + block seed offsets -> TrialOrder(block, trial, run).
TrialOrder = zeros(numBlocks, runLength, nRuns);
BlockOrder = 1:numBlocks;

for iBlock = 1:numBlocks
    % Seed each block deterministically to keep reproducibility while
    % preventing identical shuffles across blocks.
    rng(randomSeed + iBlock - 1);

    alwaysPool = idxAlways(randperm(numel(idxAlways)));
    nondevPool = idxNondev(randperm(numel(idxNondev)));
    devPool = idxDev(randperm(numel(idxDev)));

    run1Order = alwaysPool(:)';

    run2Pool = [nondevPool(1:nOccHalf), devPool(1:nOccHalf)];
    run3Pool = [nondevPool((nOccHalf + 1):end), devPool((nOccHalf + 1):end)];

    run2Order = run2Pool(randperm(numel(run2Pool)));
    run3Order = run3Pool(randperm(numel(run3Pool)));

    TrialOrder(iBlock, :, 1) = run1Order;
    TrialOrder(iBlock, :, 2) = run2Order;
    TrialOrder(iBlock, :, 3) = run3Order;
end

%% Build flattened TrialStruct for compatibility consumers
% Data flow: per-run concatenated block orders -> TrialStruct metadata rows.
nTrialsFlattened = numBlocks * runLength;
trialTemplate = struct( ...
    'Type', [], ...
    'Start', [], ...
    'Condition', [], ...
    'PathDuration', [], ...
    'VideoNumber', [], ...
    'Congruence', [], ...
    'IncongruentOffSet', [], ...
    'SimulatedPath', {{}}, ...
    'SimulatedNewAngle', {{}}, ...
    'End', [], ...
    'ConditionLabel', '', ...
    'SourceTrialIndex', [], ...
    'SequenceID', [], ...
    'OcclusionEnabled', false, ...
    'RunInBlock', [], ...
    'BlockIndex', []);

TrialStruct = repmat(trialTemplate, nRuns, nTrialsFlattened);

for iRun = 1:nRuns
    runOrderFlat = local_flatten_run_order(TrialOrder, iRun);
    for iSeq = 1:nTrialsFlattened
        sourceIdx = runOrderFlat(iSeq);
        sourceTrial = xyAll(sourceIdx);

        TrialStruct(iRun, iSeq).Type = [];
        TrialStruct(iRun, iSeq).Start = [];
        TrialStruct(iRun, iSeq).End = [];
        TrialStruct(iRun, iSeq).Congruence = [];
        TrialStruct(iRun, iSeq).IncongruentOffSet = [];
        TrialStruct(iRun, iSeq).SimulatedPath = cell(0, 1);
        TrialStruct(iRun, iSeq).SimulatedNewAngle = cell(0, 1);

        TrialStruct(iRun, iSeq).Condition = local_condition_code(sourceTrial);
        TrialStruct(iRun, iSeq).ConditionLabel = local_condition_label(sourceTrial);
        TrialStruct(iRun, iSeq).SourceTrialIndex = sourceIdx;
        TrialStruct(iRun, iSeq).SequenceID = local_sequence_id(sourceTrial, sourceIdx);
        TrialStruct(iRun, iSeq).OcclusionEnabled = local_occlusion_enabled(sourceTrial);
        TrialStruct(iRun, iSeq).RunInBlock = iRun;
        TrialStruct(iRun, iSeq).BlockIndex = ceil(iSeq / runLength);

        if isfield(sourceTrial, 'pathDuration') && ~isempty(sourceTrial.pathDuration)
            TrialStruct(iRun, iSeq).PathDuration = double(sourceTrial.pathDuration);
        else
            if isfield(Dat.Cfg, 'fps') && Dat.Cfg.fps > 0
                TrialStruct(iRun, iSeq).PathDuration = size(sourceTrial.xy, 1) / double(Dat.Cfg.fps);
            else
                TrialStruct(iRun, iSeq).PathDuration = NaN;
            end
        end

        if isfield(sourceTrial, 'sequence') && ~isempty(sourceTrial.sequence)
            TrialStruct(iRun, iSeq).VideoNumber = double(sourceTrial.sequence);
        else
            TrialStruct(iRun, iSeq).VideoNumber = double(sourceIdx);
        end
    end
end

%% Build no-catch Catch struct
% Data flow: compatibility defaults -> catch container with all zeros.
Catch = struct();
Catch.OcclAngleDisplacement = 0;
Catch.PlacingNewPos = 0;
Catch.MaxOverlap = 0;
Catch.nCatchPerSeq = 0;
Catch.CatchRatioFix = 0;
Catch.CatchStart = 0;
Catch.CatchDistance = 0;
Catch.BouncingEnd = 0;
Catch.FixDuration = 0;
Catch.FixResponseDuration = 0;
Catch.FixWaitingDuration = 0;
Catch.FixTotalDuration = 0;
Catch.OcclDuration = 0;
Catch.OcclVideoDuration = 0;
Catch.OcclResponseDuration = 0;
Catch.OcclWaitingDuration = 0;
Catch.OcclTotalDuration = 0;
Catch.CatchRatioOcc = 0;
Catch.OcclDeviance = 0;
Catch.OcclDistance = 0;
Catch.CatchDuration = [0 0];
if isfield(Dat.Cfg, 'Stimulitype')
    Catch.Stimulitype = Dat.Cfg.Stimulitype;
else
    Catch.Stimulitype = 0;
end
Catch.BouncingStart = 0;
Catch.nCatchPerRun = 0;
Catch.nCatchTotal = 0;

%% Build explicit schedule metadata
% Data flow: resolved counts + controls -> Schedule struct for runtime and audit.
Schedule = struct();
Schedule.version = 'v13_three_runs_per_block';
Schedule.numBlocks = double(numBlocks);
Schedule.runsPerBlock = double(nRuns);
Schedule.runLength = double(runLength);
Schedule.randomSeed = double(randomSeed);
Schedule.nAlwaysVisible = double(nAlways);
Schedule.nOccludedNondeviant = double(nNondev);
Schedule.nOccludedDeviant = double(nDev);
Schedule.run2PerOccludedCondition = double(nOccHalf);
Schedule.run3PerOccludedCondition = double(nOccHalf);
Schedule.design = struct( ...
    'run1', 'all always_visible shuffled', ...
    'run2', 'random half nondeviant + random half deviant shuffled', ...
    'run3', 'remaining half nondeviant + remaining half deviant shuffled');

%% Save output
% Data flow: generated compatibility structs + schedule metadata -> MAT artifact.
save(outputFile, 'TrialStruct', 'Catch', 'TrialOrder', 'BlockOrder', 'Schedule');

fprintf('Saved v13 three-runs-per-block trial struct: %s\n', outputFile);
fprintf(['Subject %02d | blocks=%d | runs/block=%d | run length=%d | ' ...
    'always=%d nondev=%d dev=%d\n'], ...
    iSub, numBlocks, nRuns, runLength, nAlways, nNondev, nDev);

%% Local helpers
function runOrderFlat = local_flatten_run_order(trialOrder, iRun)
% LOCAL_FLATTEN_RUN_ORDER Concatenate one run across all blocks.
runOrderFlat = [];
for iBlock = 1:size(trialOrder, 1)
    blockOrder = squeeze(trialOrder(iBlock, :, iRun));
    runOrderFlat = [runOrderFlat, blockOrder(:)']; %#ok<AGROW>
end
end

function condCode = local_condition_code(trial)
% LOCAL_CONDITION_CODE Map condition label to legacy numeric condition code.
if isfield(trial, 'condition') && ~isempty(trial.condition)
    condCode = double(trial.condition);
    return;
end
label = local_condition_label(trial);
switch label
    case 'always_visible'
        condCode = -1;
    case 'occluded_nondeviant'
        condCode = 0;
    case 'occluded_deviant'
        condCode = 45;
    otherwise
        condCode = NaN;
end
end

function label = local_condition_label(trial)
% LOCAL_CONDITION_LABEL Return lowercase condition label or empty string.
if isfield(trial, 'condition_label') && ~isempty(trial.condition_label)
    label = lower(char(trial.condition_label));
else
    label = '';
end
end

function seqId = local_sequence_id(trial, fallbackIdx)
% LOCAL_SEQUENCE_ID Read source sequence field if present, else fallback index.
if isfield(trial, 'sequence') && ~isempty(trial.sequence)
    seqId = double(trial.sequence);
else
    seqId = double(fallbackIdx);
end
end

function tf = local_occlusion_enabled(trial)
% LOCAL_OCCLUSION_ENABLED Read occlusion_enabled with false fallback.
if isfield(trial, 'occlusion_enabled') && ~isempty(trial.occlusion_enabled)
    tf = logical(trial.occlusion_enabled);
else
    tf = false;
end
end
