% CreateInputFiles_v14_threeRunsPerBlock_catch.m
%
% Purpose:
%   Build an additive, block-aware TrialStruct artifact for the one-dot
%   occlusion paradigm with fixed three-runs-per-block scheduling and
%   explicit catch-trial insertion metadata.
%
% Base schedule model (before catch insertion):
%   - Run 1: all always_visible trials (shuffled)
%   - Run 2: random half of occluded_nondeviant + random half of
%            occluded_deviant (shuffled together)
%   - Run 3: remaining half of occluded_nondeviant + remaining half of
%            occluded_deviant (shuffled together)
%
% Catch schedule model (v14 additions):
%   - Catch type 1 (always_visible disappear/reappear) is inserted in run 1.
%   - Catch type 2 (normal occlusion trial + end question) is inserted
%     across runs 2 and 3.
%   - Catch counts are percentage-based and rounded:
%       nCatchType1 = round(catchRateType1Run1 * nRun1BaseTrials)
%       nCatchType2 = round(catchRateType2Runs23 * (nRun2BaseTrials + nRun3BaseTrials))
%
% Why this v14 file exists:
%   - Keep existing v13 scheduling scripts untouched.
%   - Add catch planning upstream (input creation layer) so runtime can
%     consume deterministic per-slot catch metadata.
%
% Usage example (interactive):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment');
%   CreateInputFiles_v14_threeRunsPerBlock_catch
%
% Usage example (non-interactive):
%   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
%   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
%    iSub=66; randomSeed=6601; numBlocks=1; ...
%    run('CreateInputFiles_v14_threeRunsPerBlock_catch.m');"
%
% Inputs:
%   - input_files/MovDot_SubXX.mat with:
%       * xySeqs
%       * Cfg
%   - Config_schedule_CreateInputV14_MoveDotV8 constants:
%       * schedule controls (numBlocks, runsPerBlock, run2 split)
%       * catch controls (rates, disappear range, gap, prompt settings)
%
% Workspace overrides (optional):
%   - iSub
%   - randomSeed
%   - numBlocks
%   - overwriteExisting
%
% Outputs:
%   - input_files/SubXX_TrialStruct_v14_threeRunsPerBlock_catch.mat with:
%       * TrialStruct  (nRuns x nTrialsFlattened; padded slots allowed)
%       * Catch        (legacy compatibility struct + v14 catch counts)
%       * TrialOrder   (numBlocks x maxRunLength x nRuns; NaN padded)
%       * CatchPlan    (per-slot catch metadata aligned to TrialOrder dims)
%       * BlockOrder   (1 x numBlocks)
%       * Schedule     (explicit schedule + catch settings)
%
% Key assumptions:
%   - Input labels include exactly:
%       always_visible, occluded_nondeviant, occluded_deviant
%   - Base occluded condition counts are equal and even.
%   - Per-trial sequence IDs are available or can be approximated.

%% Setup and input collection
% Data flow: dialog/workspace overrides + config -> file names and controls.
clearvars -except iSub randomSeed numBlocks overwriteExisting;
clc;
addpath('lib/');

scheduleCfg = Config_schedule_CreateInputV14_MoveDotV8;
if scheduleCfg.runsPerBlock ~= 3
    error('Config_schedule_CreateInputV14_MoveDotV8.runsPerBlock must be 3 for this script.');
end
if abs(scheduleCfg.run2FractionPerOccluded - 0.5) > eps
    error('Config_schedule_CreateInputV14_MoveDotV8.run2FractionPerOccluded must be 0.5.');
end

if ~(exist('iSub', 'var') && exist('randomSeed', 'var'))
    prompt = {'Subject Number:', 'Random seed:'};
    dlgtitle = 'CreateInputFiles v14 (three runs per block + catches)';
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

if logical(scheduleCfg.enforceEqualBaseRunLengths) && (nAlways ~= nNondev)
    error(['For equal base run lengths, always_visible (%d) must match each ' ...
        'occluded count (%d).'], nAlways, nNondev);
end

framesPerTrial = size(xyAll(1).xy, 1);
fps = double(Dat.Cfg.fps);
if ~isfinite(fps) || fps <= 0
    error('Cfg.fps must be positive.');
end

nOccPerCondition = nNondev;
nOccHalf = nOccPerCondition / 2;
baseRunLength = nAlways;
nRuns = scheduleCfg.runsPerBlock;

if baseRunLength ~= 2 * nOccHalf
    error(['Base run length mismatch: run1=%d but run2/run3 base=%d. ' ...
        'Use matched condition counts for v14.'], baseRunLength, 2 * nOccHalf);
end

%% Build sequence->deviant index map for type-1 plausible-path catches
% Data flow: source trial metadata -> lookup table used by catch insertion.
seqToDevIdx = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:numel(idxDev)
    trialIdx = idxDev(k);
    seqId = local_sequence_id(xyAll(trialIdx), trialIdx);
    if ~isKey(seqToDevIdx, seqId)
        seqToDevIdx(seqId) = double(trialIdx);
    end
end

%% Build block-aware TrialOrder + CatchPlan (variable run lengths, NaN padded later)
% Data flow: condition pools + catch controls -> per-block run orders/catch slots.
BlockOrder = 1:numBlocks;

runOrders = cell(numBlocks, nRuns);
runCatchType = cell(numBlocks, nRuns);
runCatchExpected = cell(numBlocks, nRuns);
runCatchBranchChanged = cell(numBlocks, nRuns);
runCatchDisappearFrame = cell(numBlocks, nRuns);
runCatchReappearFrame = cell(numBlocks, nRuns);
runCatchAltSourceIdx = cell(numBlocks, nRuns);
runLengthsByBlockRun = zeros(numBlocks, nRuns);

catchCountType1Total = 0;
catchCountType2Total = 0;

for iBlock = 1:numBlocks
    % Seed each block deterministically to keep reproducibility while
    % preventing identical shuffles across blocks.
    rng(randomSeed + iBlock - 1);

    alwaysPool = idxAlways(randperm(numel(idxAlways)));
    nondevPool = idxNondev(randperm(numel(idxNondev)));
    devPool = idxDev(randperm(numel(idxDev)));

    baseRun1 = alwaysPool(:)';
    baseRun2Pool = [nondevPool(1:nOccHalf), devPool(1:nOccHalf)];
    baseRun3Pool = [nondevPool((nOccHalf + 1):end), devPool((nOccHalf + 1):end)];
    baseRun2 = baseRun2Pool(randperm(numel(baseRun2Pool)));
    baseRun3 = baseRun3Pool(randperm(numel(baseRun3Pool)));

    % Determine catch counts for this block.
    nCatchType1 = round(double(scheduleCfg.catchRateType1Run1) * numel(baseRun1));
    nCatchType2 = round(double(scheduleCfg.catchRateType2Runs23) * (numel(baseRun2) + numel(baseRun3)));
    nCatchType1 = max(0, nCatchType1);
    nCatchType2 = max(0, nCatchType2);

    catchSpecsType1 = local_make_type1_catch_specs( ...
        nCatchType1, baseRun1, xyAll, seqToDevIdx, fps, framesPerTrial, scheduleCfg);

    [run1Order, run1CatchMeta] = local_insert_catches_into_run(baseRun1, catchSpecsType1);

    % Type-2 catches are sampled from the combined base run2+run3 pools.
    catchSpecsType2Combined = local_make_type2_catch_specs(nCatchType2, baseRun2, baseRun3, xyAll, scheduleCfg);
    catchSpecsType2Run2 = catchSpecsType2Combined([catchSpecsType2Combined.targetRun] == 2);
    catchSpecsType2Run3 = catchSpecsType2Combined([catchSpecsType2Combined.targetRun] == 3);

    [run2Order, run2CatchMeta] = local_insert_catches_into_run(baseRun2, catchSpecsType2Run2);
    [run3Order, run3CatchMeta] = local_insert_catches_into_run(baseRun3, catchSpecsType2Run3);

    runOrders{iBlock, 1} = run1Order;
    runOrders{iBlock, 2} = run2Order;
    runOrders{iBlock, 3} = run3Order;

    runCatchType{iBlock, 1} = run1CatchMeta.typeCode;
    runCatchType{iBlock, 2} = run2CatchMeta.typeCode;
    runCatchType{iBlock, 3} = run3CatchMeta.typeCode;

    runCatchExpected{iBlock, 1} = run1CatchMeta.expectedResponse;
    runCatchExpected{iBlock, 2} = run2CatchMeta.expectedResponse;
    runCatchExpected{iBlock, 3} = run3CatchMeta.expectedResponse;

    runCatchBranchChanged{iBlock, 1} = run1CatchMeta.branchChanged;
    runCatchBranchChanged{iBlock, 2} = run2CatchMeta.branchChanged;
    runCatchBranchChanged{iBlock, 3} = run3CatchMeta.branchChanged;

    runCatchDisappearFrame{iBlock, 1} = run1CatchMeta.disappearFrame;
    runCatchDisappearFrame{iBlock, 2} = run2CatchMeta.disappearFrame;
    runCatchDisappearFrame{iBlock, 3} = run3CatchMeta.disappearFrame;

    runCatchReappearFrame{iBlock, 1} = run1CatchMeta.reappearFrame;
    runCatchReappearFrame{iBlock, 2} = run2CatchMeta.reappearFrame;
    runCatchReappearFrame{iBlock, 3} = run3CatchMeta.reappearFrame;

    runCatchAltSourceIdx{iBlock, 1} = run1CatchMeta.altSourceIdx;
    runCatchAltSourceIdx{iBlock, 2} = run2CatchMeta.altSourceIdx;
    runCatchAltSourceIdx{iBlock, 3} = run3CatchMeta.altSourceIdx;

    runLengthsByBlockRun(iBlock, 1) = numel(run1Order);
    runLengthsByBlockRun(iBlock, 2) = numel(run2Order);
    runLengthsByBlockRun(iBlock, 3) = numel(run3Order);

    catchCountType1Total = catchCountType1Total + sum(run1CatchMeta.typeCode == 1);
    catchCountType2Total = catchCountType2Total + sum(run2CatchMeta.typeCode == 2) + sum(run3CatchMeta.typeCode == 2);
end

maxRunLength = max(runLengthsByBlockRun(:));

%% Materialize padded TrialOrder and padded CatchPlan arrays
% Data flow: variable-length per-block/per-run vectors -> fixed-size MAT arrays.
TrialOrder = nan(numBlocks, maxRunLength, nRuns);
CatchPlan = struct();
CatchPlan.typeCode = zeros(numBlocks, maxRunLength, nRuns);
CatchPlan.expectedResponseCode = zeros(numBlocks, maxRunLength, nRuns);
CatchPlan.branchChangedPath = nan(numBlocks, maxRunLength, nRuns);
CatchPlan.disappearFrame = nan(numBlocks, maxRunLength, nRuns);
CatchPlan.reappearFrame = nan(numBlocks, maxRunLength, nRuns);
CatchPlan.altSourceTrialIndex = nan(numBlocks, maxRunLength, nRuns);

for iBlock = 1:numBlocks
    for iRun = 1:nRuns
        thisOrder = runOrders{iBlock, iRun};
        thisLen = numel(thisOrder);
        TrialOrder(iBlock, 1:thisLen, iRun) = thisOrder;

        CatchPlan.typeCode(iBlock, 1:thisLen, iRun) = runCatchType{iBlock, iRun};
        CatchPlan.expectedResponseCode(iBlock, 1:thisLen, iRun) = runCatchExpected{iBlock, iRun};
        CatchPlan.branchChangedPath(iBlock, 1:thisLen, iRun) = runCatchBranchChanged{iBlock, iRun};
        CatchPlan.disappearFrame(iBlock, 1:thisLen, iRun) = runCatchDisappearFrame{iBlock, iRun};
        CatchPlan.reappearFrame(iBlock, 1:thisLen, iRun) = runCatchReappearFrame{iBlock, iRun};
        CatchPlan.altSourceTrialIndex(iBlock, 1:thisLen, iRun) = runCatchAltSourceIdx{iBlock, iRun};
    end
end
CatchPlan.typeLabelMap = {'none', 'always_visible_disappear_reappear', 'occlusion_question_only'};
CatchPlan.responseCodeMap = struct('none', 0, 'yes', scheduleCfg.catchResponseYesCode, 'no', scheduleCfg.catchResponseNoCode);

%% Build flattened TrialStruct for compatibility consumers
% Data flow: padded TrialOrder + CatchPlan -> flattened rows grouped by run.
nTrialsFlattened = numBlocks * maxRunLength;
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
    'BlockIndex', [], ...
    'CatchTypeCode', [], ...
    'CatchExpectedResponse', [], ...
    'CatchBranchChangedPath', [], ...
    'CatchDisappearFrame', [], ...
    'CatchReappearFrame', [], ...
    'CatchAltSourceIndex', []);

TrialStruct = repmat(trialTemplate, nRuns, nTrialsFlattened);

for iRun = 1:nRuns
    runOrderFlat = local_flatten_run_order(TrialOrder, iRun);
    catchTypeFlat = local_flatten_run_order(CatchPlan.typeCode, iRun);
    catchExpectedFlat = local_flatten_run_order(CatchPlan.expectedResponseCode, iRun);
    catchBranchFlat = local_flatten_run_order(CatchPlan.branchChangedPath, iRun);
    catchDisappearFlat = local_flatten_run_order(CatchPlan.disappearFrame, iRun);
    catchReappearFlat = local_flatten_run_order(CatchPlan.reappearFrame, iRun);
    catchAltIdxFlat = local_flatten_run_order(CatchPlan.altSourceTrialIndex, iRun);

    for iSeq = 1:nTrialsFlattened
        sourceIdx = runOrderFlat(iSeq);
        blockIdx = ceil(iSeq / maxRunLength);

        TrialStruct(iRun, iSeq).Type = [];
        TrialStruct(iRun, iSeq).Start = [];
        TrialStruct(iRun, iSeq).End = [];
        TrialStruct(iRun, iSeq).Congruence = [];
        TrialStruct(iRun, iSeq).IncongruentOffSet = [];
        TrialStruct(iRun, iSeq).SimulatedPath = cell(0, 1);
        TrialStruct(iRun, iSeq).SimulatedNewAngle = cell(0, 1);
        TrialStruct(iRun, iSeq).RunInBlock = iRun;
        TrialStruct(iRun, iSeq).BlockIndex = blockIdx;

        TrialStruct(iRun, iSeq).CatchTypeCode = catchTypeFlat(iSeq);
        TrialStruct(iRun, iSeq).CatchExpectedResponse = catchExpectedFlat(iSeq);
        TrialStruct(iRun, iSeq).CatchBranchChangedPath = catchBranchFlat(iSeq);
        TrialStruct(iRun, iSeq).CatchDisappearFrame = catchDisappearFlat(iSeq);
        TrialStruct(iRun, iSeq).CatchReappearFrame = catchReappearFlat(iSeq);
        TrialStruct(iRun, iSeq).CatchAltSourceIndex = catchAltIdxFlat(iSeq);

        if ~isfinite(sourceIdx) || sourceIdx < 1 || sourceIdx > numel(xyAll)
            TrialStruct(iRun, iSeq).Condition = NaN;
            TrialStruct(iRun, iSeq).ConditionLabel = '';
            TrialStruct(iRun, iSeq).SourceTrialIndex = NaN;
            TrialStruct(iRun, iSeq).SequenceID = NaN;
            TrialStruct(iRun, iSeq).OcclusionEnabled = false;
            TrialStruct(iRun, iSeq).PathDuration = NaN;
            TrialStruct(iRun, iSeq).VideoNumber = NaN;
            continue;
        end

        sourceIdx = round(sourceIdx);
        sourceTrial = xyAll(sourceIdx);

        TrialStruct(iRun, iSeq).Condition = local_condition_code(sourceTrial);
        TrialStruct(iRun, iSeq).ConditionLabel = local_condition_label(sourceTrial);
        TrialStruct(iRun, iSeq).SourceTrialIndex = sourceIdx;
        TrialStruct(iRun, iSeq).SequenceID = local_sequence_id(sourceTrial, sourceIdx);
        TrialStruct(iRun, iSeq).OcclusionEnabled = local_occlusion_enabled(sourceTrial);

        if isfield(sourceTrial, 'pathDuration') && ~isempty(sourceTrial.pathDuration)
            TrialStruct(iRun, iSeq).PathDuration = double(sourceTrial.pathDuration);
        else
            TrialStruct(iRun, iSeq).PathDuration = size(sourceTrial.xy, 1) / fps;
        end

        if isfield(sourceTrial, 'sequence') && ~isempty(sourceTrial.sequence)
            TrialStruct(iRun, iSeq).VideoNumber = double(sourceTrial.sequence);
        else
            TrialStruct(iRun, iSeq).VideoNumber = double(sourceIdx);
        end
    end
end

%% Build compatibility Catch struct
% Data flow: legacy fields + v14 counters -> downstream-compatible container.
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

% Explicit v14 catch counters.
Catch.V14 = struct();
Catch.V14.nCatchType1Total = catchCountType1Total;
Catch.V14.nCatchType2Total = catchCountType2Total;
Catch.V14.nCatchTotal = catchCountType1Total + catchCountType2Total;
Catch.V14.catchRateType1Run1 = double(scheduleCfg.catchRateType1Run1);
Catch.V14.catchRateType2Runs23 = double(scheduleCfg.catchRateType2Runs23);

%% Build explicit schedule metadata
% Data flow: resolved counts + controls -> Schedule struct for runtime and audit.
Schedule = struct();
Schedule.version = 'v14_three_runs_per_block_with_catches';
Schedule.numBlocks = double(numBlocks);
Schedule.runsPerBlock = double(nRuns);
Schedule.baseRunLength = double(baseRunLength);
Schedule.maxRunLength = double(maxRunLength);
Schedule.runLengthsByBlockRun = double(runLengthsByBlockRun);
Schedule.randomSeed = double(randomSeed);
Schedule.framesPerTrial = double(framesPerTrial);
Schedule.fps = double(fps);
Schedule.nAlwaysVisible = double(nAlways);
Schedule.nOccludedNondeviant = double(nNondev);
Schedule.nOccludedDeviant = double(nDev);
Schedule.run2PerOccludedCondition = double(nOccHalf);
Schedule.run3PerOccludedCondition = double(nOccHalf);
Schedule.catchSettings = struct( ...
    'catchRateType1Run1', double(scheduleCfg.catchRateType1Run1), ...
    'catchRateType2Runs23', double(scheduleCfg.catchRateType2Runs23), ...
    'type1DisappearRangeSec', double(scheduleCfg.catchType1DisappearRangeSec), ...
    'type1InvisibleDurationSec', double(scheduleCfg.catchType1InvisibleDurationSec), ...
    'type1ChangedPathProbability', double(scheduleCfg.catchType1ChangedPathProbability), ...
    'questionText', char(scheduleCfg.catchQuestionText), ...
    'questionTimeoutSec', double(scheduleCfg.catchQuestionTimeoutSec), ...
    'responseYesCode', double(scheduleCfg.catchResponseYesCode), ...
    'responseNoCode', double(scheduleCfg.catchResponseNoCode));
Schedule.catchCounts = struct( ...
    'type1Total', double(catchCountType1Total), ...
    'type2Total', double(catchCountType2Total), ...
    'total', double(catchCountType1Total + catchCountType2Total));
Schedule.design = struct( ...
    'run1', 'always_visible + inserted type1 catches', ...
    'run2', 'occluded pools + subset of inserted type2 catches', ...
    'run3', 'occluded pools + subset of inserted type2 catches');

%% Save output
% Data flow: generated structs + metadata -> MAT artifact.
save(outputFile, 'TrialStruct', 'Catch', 'TrialOrder', 'CatchPlan', 'BlockOrder', 'Schedule');

fprintf('Saved v14 three-runs-per-block catch trial struct: %s\n', outputFile);
fprintf(['Subject %02d | blocks=%d | runs/block=%d | base run=%d | max run=%d | ' ...
    'type1 catches=%d | type2 catches=%d\n'], ...
    iSub, numBlocks, nRuns, baseRunLength, maxRunLength, catchCountType1Total, catchCountType2Total);

%% Local helpers
function runOrderFlat = local_flatten_run_order(trialOrder, iRun)
% LOCAL_FLATTEN_RUN_ORDER Concatenate one run across all blocks.
runOrderFlat = [];
for iBlock = 1:size(trialOrder, 1)
    blockOrder = squeeze(trialOrder(iBlock, :, iRun));
    runOrderFlat = [runOrderFlat, blockOrder(:)']; %#ok<AGROW>
end
end

function specs = local_make_type1_catch_specs(nCatch, baseRun1, xyAll, seqToDevIdx, fps, framesPerTrial, scheduleCfg)
% LOCAL_MAKE_TYPE1_CATCH_SPECS Build type-1 catch metadata for run 1.
specTemplate = struct( ...
    'typeCode', 1, ...
    'expectedResponse', scheduleCfg.catchResponseNoCode, ...
    'branchChanged', 0, ...
    'disappearFrame', NaN, ...
    'reappearFrame', NaN, ...
    'altSourceIdx', NaN, ...
    'targetRun', 1, ...
    'sourceIdx', NaN);
specs = repmat(specTemplate, 1, nCatch);
if nCatch == 0
    return;
end

sourceCandidates = baseRun1(:)';
if isempty(sourceCandidates)
    specs = specs([]); % no valid source slots
    return;
end

sampleWithReplacement = nCatch > numel(sourceCandidates);
if sampleWithReplacement
    pickIdx = randi(numel(sourceCandidates), [1, nCatch]);
else
    pickIdx = randperm(numel(sourceCandidates), nCatch);
end
chosenSources = sourceCandidates(pickIdx);

for iCatch = 1:nCatch
    srcIdx = chosenSources(iCatch);
    srcTrial = xyAll(srcIdx);
    seqId = local_sequence_id(srcTrial, srcIdx);

    branchChanged = rand < double(scheduleCfg.catchType1ChangedPathProbability);
    altIdx = srcIdx;
    if branchChanged && isKey(seqToDevIdx, seqId)
        altIdx = seqToDevIdx(seqId);
    else
        branchChanged = false;
    end

    devFrameAlt = local_trial_deviance_frame(xyAll(altIdx), framesPerTrial);
    disappearFrame = local_sample_disappear_frame( ...
        scheduleCfg.catchType1DisappearRangeSec, ...
        scheduleCfg.catchType1InvisibleDurationSec, ...
        fps, ...
        framesPerTrial, ...
        branchChanged, ...
        devFrameAlt);
    gapFrames = max(1, round(double(scheduleCfg.catchType1InvisibleDurationSec) * fps));
    reappearFrame = min(framesPerTrial, disappearFrame + gapFrames);

    specs(iCatch).typeCode = 1;
    specs(iCatch).expectedResponse = local_expected_response_for_type1(branchChanged, scheduleCfg);
    specs(iCatch).branchChanged = double(branchChanged);
    specs(iCatch).disappearFrame = double(disappearFrame);
    specs(iCatch).reappearFrame = double(reappearFrame);
    specs(iCatch).altSourceIdx = double(altIdx);
    specs(iCatch).targetRun = 1;
    specs(iCatch).sourceIdx = double(srcIdx);
end
end

function specs = local_make_type2_catch_specs(nCatch, baseRun2, baseRun3, xyAll, scheduleCfg)
% LOCAL_MAKE_TYPE2_CATCH_SPECS Build type-2 catch metadata across runs 2 and 3.
specTemplate = struct( ...
    'typeCode', 2, ...
    'expectedResponse', scheduleCfg.catchResponseNoCode, ...
    'branchChanged', NaN, ...
    'disappearFrame', NaN, ...
    'reappearFrame', NaN, ...
    'altSourceIdx', NaN, ...
    'targetRun', NaN, ...
    'sourceIdx', NaN);
specs = repmat(specTemplate, 1, nCatch);
if nCatch == 0
    return;
end

combined = struct('targetRun', {}, 'sourceIdx', {});
for k = 1:numel(baseRun2)
    combined(end + 1) = struct('targetRun', 2, 'sourceIdx', double(baseRun2(k))); %#ok<AGROW>
end
for k = 1:numel(baseRun3)
    combined(end + 1) = struct('targetRun', 3, 'sourceIdx', double(baseRun3(k))); %#ok<AGROW>
end

if isempty(combined)
    specs = specs([]);
    return;
end

nCatch = min(nCatch, numel(combined));
pick = randperm(numel(combined), nCatch);
chosen = combined(pick);

for iCatch = 1:nCatch
    srcIdx = chosen(iCatch).sourceIdx;
    srcTrial = xyAll(srcIdx);
    label = local_condition_label(srcTrial);

    expected = scheduleCfg.catchResponseNoCode;
    if strcmp(label, 'occluded_deviant')
        expected = scheduleCfg.catchResponseYesCode;
    end

    specs(iCatch).typeCode = 2;
    specs(iCatch).expectedResponse = double(expected);
    specs(iCatch).branchChanged = NaN;
    specs(iCatch).disappearFrame = NaN;
    specs(iCatch).reappearFrame = NaN;
    specs(iCatch).altSourceIdx = NaN;
    specs(iCatch).targetRun = double(chosen(iCatch).targetRun);
    specs(iCatch).sourceIdx = double(srcIdx);
end
end

function [augmentedOrder, meta] = local_insert_catches_into_run(baseOrder, catchSpecs)
% LOCAL_INSERT_CATCHES_INTO_RUN Insert catch trials into a base run order.
%
% Data flow:
%   base non-catch order + catch specs -> one augmented run with per-slot
%   metadata aligned to the augmented order.
nBase = numel(baseOrder);
nCatch = numel(catchSpecs);
nFinal = nBase + nCatch;

augmentedOrder = nan(1, nFinal);
meta = struct();
meta.typeCode = zeros(1, nFinal);
meta.expectedResponse = zeros(1, nFinal);
meta.branchChanged = nan(1, nFinal);
meta.disappearFrame = nan(1, nFinal);
meta.reappearFrame = nan(1, nFinal);
meta.altSourceIdx = nan(1, nFinal);

if nCatch > 0
    catchSlots = sort(randperm(nFinal, nCatch));
    catchSpecs = catchSpecs(randperm(nCatch));
else
    catchSlots = [];
end

iBase = 1;
iCatch = 1;
for iPos = 1:nFinal
    isCatchSlot = iCatch <= nCatch && iPos == catchSlots(iCatch);
    if isCatchSlot
        spec = catchSpecs(iCatch);
        augmentedOrder(iPos) = spec.sourceIdx;
        meta.typeCode(iPos) = spec.typeCode;
        meta.expectedResponse(iPos) = spec.expectedResponse;
        meta.branchChanged(iPos) = spec.branchChanged;
        meta.disappearFrame(iPos) = spec.disappearFrame;
        meta.reappearFrame(iPos) = spec.reappearFrame;
        meta.altSourceIdx(iPos) = spec.altSourceIdx;
        iCatch = iCatch + 1;
    else
        augmentedOrder(iPos) = baseOrder(iBase);
        iBase = iBase + 1;
    end
end
end

function expected = local_expected_response_for_type1(branchChanged, scheduleCfg)
% LOCAL_EXPECTED_RESPONSE_FOR_TYPE1 Map type-1 branch to yes/no code.
if branchChanged
    expected = double(scheduleCfg.catchResponseYesCode);
else
    expected = double(scheduleCfg.catchResponseNoCode);
end
end

function frameIdx = local_sample_disappear_frame(rangeSec, invisibleSec, fps, framesPerTrial, enforceHiddenDeviance, devFrame)
% LOCAL_SAMPLE_DISAPPEAR_FRAME Sample type-1 disappearance frame.
%
% For changed-path catches, try to enforce that deviance happens while the
% dot is hidden by sampling from the overlap between:
%   configured range and [devFrame-gap+1, devFrame].
rangeFrames = sort(round(double(rangeSec(:)') * fps));
if numel(rangeFrames) ~= 2
    rangeFrames = [1, max(1, framesPerTrial - 1)];
end
minFrame = max(2, rangeFrames(1));
maxFrame = min(framesPerTrial - 1, rangeFrames(2));
if maxFrame < minFrame
    maxFrame = minFrame;
end

candidateFrames = minFrame:maxFrame;
if enforceHiddenDeviance
    gapFrames = max(1, round(double(invisibleSec) * fps));
    devWindow = (devFrame - gapFrames + 1):devFrame;
    constrained = intersect(candidateFrames, devWindow);
    if ~isempty(constrained)
        candidateFrames = constrained;
    end
end

frameIdx = candidateFrames(randi(numel(candidateFrames)));
frameIdx = max(2, min(framesPerTrial - 1, frameIdx));
end

function devFrame = local_trial_deviance_frame(trial, framesPerTrial)
% LOCAL_TRIAL_DEVIANCE_FRAME Return clamped deviance frame for one trial.
if isfield(trial, 'deviance_frame') && ~isempty(trial.deviance_frame)
    devFrame = round(double(trial.deviance_frame));
else
    devFrame = round(framesPerTrial / 2);
end
devFrame = max(1, min(framesPerTrial, devFrame));
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
