%% Stimuli Generation (uniform starts, boundary-safe placement, dRSA-proxy trial gate)
% Script: stimuli_generation_v15_experimentalPathScale.m
% Author: Marisa (original), Ayman Hatoum (v5), updated by Dami (v6-v8), v9-v15 updates by Codex
%
% Purpose:
%   Generate dot-motion stimuli with uniform starting positions while keeping
%   the rectangular bounds strictly off-limits. This version removes
%   segmented sub-paths and center-biased start points, and instead:
%     - generates trajectories in relative coordinates first,
%     - solves feasible start-position ranges that keep paths in bounds, and
%     - avoids boundary-hit rejection bias from whole-trial resampling.
%   In practice this keeps the v13 path family while reducing hidden
%   selection pressure that can couple position and direction models in dRSA.
%
%   When a deviant occurs (likelihood stimulus), curvature can be modulated
%   from that frame onward for the deviant path (via Config) by either:
%     - flipping curvature sign (legacy behavior), or
%     - sampling a new curvature value in a configurable +/- range.
%
%   v14 keeps curvature smooth and predictable by design:
%     - one constant curvature value per dot at trial start,
%     - optional deviant-only curvature change at deviant onset,
%     - no additional within-trial curvature updates.
%
%   To reduce residual dRSA position-direction cross-correlation without
%   changing trajectory smoothness, v14 adds a dRSA-proxy-aware gate:
%     - generate a candidate trial using the same constant-curvature rules,
%     - compute a proxy of the target dRSA cross-matrix by correlating
%       position-RDM columns (euclidean) against direction-RDM columns
%       (cosine) over sampled frames in the accepted-trial bank plus
%       candidate trial,
%     - accept candidates only when the proxy score improves (or meets an
%       optional absolute cap).
%   This preserves path smoothness while shaping the trial ensemble directly
%   toward lower position-dot1 vs direction-dot1 dRSA coupling.
%
%   v15 keeps the same constant-curvature rules and deviant-point behavior as
%   v14, but adds a local `pathScale` multiplier on step size so path spatial
%   extent can be reduced without changing frames-per-trial or deviant timing.
%   By default, `pathScale = 0.5` produces half-length paths.
%
%   Optionally plots a random subset of paths per condition after generation.
%   This version also saves the no-deviant baseline paths for deviant
%   conditions to a separate analysis-only file.
%
%   Output matches the structure expected by:
%     - experiment/MoveDot1_experiment_vX.m
%     - simulations/build_movdot_simulation_inputs.m
%
% Example usage (from repo root in MATLAB):
%   addpath('experiment');
%   stimuli_generation_v15_experimentalPathScale;
%   % Follow the dialog prompts for viewing distance and subject ID.
%
% Example usage (custom working directory):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   addpath('experiment');
%   stimuli_generation_v15_experimentalPathScale;
%
% Inputs:
%   - Config (class/struct on MATLAB path) with screen, dot, and timing params.
%   - Utils helper class on MATLAB path.
%   - Dialog values for viewing distance and subject ID.
%
% Outputs:
%   - Saves xySeqs, Cfg, and repro to Config.inputDirectory using Config.stimuliFileName.
%   - repro captures Config, script settings, and RNG state needed to reproduce paths.
%   - Cfg.dpf stores the effective (scaled) per-frame displacement and
%     Cfg.pathScale stores the local spatial scaling factor.
%   - xySeqs(...).xy is frames x 4 with columns [x1 y1 x2 y2] in visual degrees.
%   - Saves xySeqsPredicted and Cfg to input_files/MovDot_SubXX_predicted.mat.
%     The "predicted" file stores no-deviant baseline paths for deviant conditions
%     only (analysis use; not presented to participants).
%   - Optional per-condition path plots after generation using time colors.
%
% Key assumptions:
%   - Coordinates are in visual degrees in the rectangle [0..rectSize].
%   - Deviant logic, curviness, and step size are preserved from v5.
%   - Dots must remain at least Config.minDistanceBetweenDots apart.
%   - Curvature modulation affects only the deviant path after deviant onset.
%   - If both deviant-curvature options are enabled, randomization takes
%     precedence over sign flipping.
%   - `pathScale` applies only to per-frame displacement, so it changes path
%     spatial extent while preserving frame count and deviant timing.
%   - Curvature stays constant within each trial except optional deviant-point
%     modulation in deviant conditions.
%   - No-deviant baselines reuse the same start points and pre-deviant
%     curvature, but remove deviant-only curvature modulation and angle changes.
%   - Boundary checks are satisfied by construction via feasible placement
%     (instead of rejecting boundary-violating starts).
%   - The dRSA-proxy gate rejects whole candidate trials; it does not edit
%     trajectories frame-by-frame.

addpath('lib/');

%% Setup environment
% Data flow: user dialog -> viewing distance/subject ID -> RNG seed.
clear all;
clc;
ShowCursor;
renderPreview = false; % set true to draw stimuli during generation
plotPathsAtEnd = true; % set true to plot path samples after generation
plotPathsPerCondition = 20; % number of paths to plot per condition (randomly sampled)
plotPathConditions = [0 45]; % directionVariance values to plot as conditions
pathScale = 0.5; % local path-length scale on dotSpeedDegPerFrame (0.5 -> half-length paths)
flipCurvatureOnDeviant = Config.flipCurvatureOnDeviant; % set true to flip curvature sign at deviant onset
randomizeCurvatureOnDeviant = Config.randomizeCurvatureOnDeviant; % set true to sample a new curvature at deviant onset
deviantCurvatureRange = Config.deviantCurvatureRange; % sampled post-deviant curvature range is [-range, +range]
enforceCurvatureFeasibilityFloor = true; % set true to avoid near-straight trajectories that cannot fit bounds
enforceMinDistanceOnNoDevBaseline = true; % set true to apply min-distance checks to the analysis-only no-deviant path too
maxAttemptsPerTrial = 60000; % hard stop to avoid infinite loops under incompatible parameter sets
enableDrsaProxyGate = true; % set true to enable dRSA-style position-vs-direction proxy gating
applyDrsaProxyGateToNondeviantOnly = true; % set true to apply gate only in nondeviant likelihood condition(s)
drsaProxyGateMinAcceptedTrials = 12; % gate starts only after this many accepted trials in current condition
drsaProxyGateFrameStride = 4; % evaluate proxy every N frames to control compute cost
drsaProxyGateExcludeMainDiagLags = 1; % ignore |lag|<=N sampled bins when averaging proxy matrix
drsaProxyGateMinScoreImprovement = 0.0005; % candidate must improve proxy score by at least this amount
drsaProxyGateMaxMeanAbsCorr = inf; % optional absolute cap (set finite to enforce hard maximum)
if deviantCurvatureRange < 0
    error('Config.deviantCurvatureRange must be >= 0.');
end
if ~isscalar(pathScale) || ~isfinite(pathScale) || pathScale <= 0
    error('pathScale must be a finite scalar > 0.');
end
if drsaProxyGateMinAcceptedTrials < 0 || ...
        floor(drsaProxyGateMinAcceptedTrials) ~= drsaProxyGateMinAcceptedTrials
    error('drsaProxyGateMinAcceptedTrials must be a non-negative integer.');
end
if drsaProxyGateFrameStride < 1 || floor(drsaProxyGateFrameStride) ~= drsaProxyGateFrameStride
    error('drsaProxyGateFrameStride must be a positive integer.');
end
if drsaProxyGateExcludeMainDiagLags < 0 || ...
        floor(drsaProxyGateExcludeMainDiagLags) ~= drsaProxyGateExcludeMainDiagLags
    error('drsaProxyGateExcludeMainDiagLags must be a non-negative integer.');
end
if drsaProxyGateMinScoreImprovement < 0
    error('drsaProxyGateMinScoreImprovement must be >= 0.');
end
if ~(isinf(drsaProxyGateMaxMeanAbsCorr) || drsaProxyGateMaxMeanAbsCorr >= 0)
    error('drsaProxyGateMaxMeanAbsCorr must be >= 0 or inf.');
end

userDoubles = Utils.GetUserDoubles(Config.dialogTitle, Config.dialogDimensions, ...
    Config.dialogPrompts, Config.dialogDefaults);
viewingDistance = userDoubles(1);
subjectID = userDoubles(2);
% Centralized scaled step used by all trajectory updates in v15.
scaledDotSpeedDegPerFrame = Config.dotSpeedDegPerFrame * pathScale;

% Seed RNG for reproducibility (per subject)
rng(subjectID);
% Snapshot RNG state after seeding to reproduce path generation later.
rngState = rng;

% Ensure output folders exist before saving
if ~exist(Config.inputDirectory, 'dir')
    mkdir(Config.inputDirectory);
end
if ~exist(Config.outputDirectory, 'dir')
    mkdir(Config.outputDirectory);
end

%% Output file guard
% Data flow: subject ID -> output files -> user choice -> generate/plot path.
outputFile = [Config.inputDirectory sprintf(Config.stimuliFileName, subjectID)];
predictedFile = [Config.inputDirectory sprintf('MovDot_Sub%02d_predicted.mat', subjectID)];
skipGeneration = false;
skipSave = false;
if exist(outputFile, 'file')
    userChoice = '';
    while ~ismember(lower(userChoice), {'o', 'p', 'c'})
        userChoice = input(sprintf('Output exists (%s). [O]verwrite, [P]lot existing, or [C]ancel? ', outputFile), 's');
    end
    if strcmpi(userChoice, 'c')
        disp('Canceled by user.');
        return;
    elseif strcmpi(userChoice, 'p')
        loadedData = load(outputFile, 'xySeqs', 'Cfg');
        xySeqs = loadedData.xySeqs;
        Cfg = loadedData.Cfg;
        framesPerTrial = size(xySeqs(1, 1, 1).xy, 1);
        plotPathsAtEnd = true;
        skipGeneration = true;
        skipSave = true;
    end
end

%% Screen and geometry setup
% Data flow: screen size + viewing distance -> pixels-per-degree -> dot rect.
if renderPreview
    Screen('Preference', 'SkipSyncTests', 1);
    PsychDebugWindowConfiguration([], .02);

    monitorIndex = max(Screen('Screens'));
    monitorPixelSize = Screen('Rect', monitorIndex);
    [monitorMilliWidth, monitorMilliHeight] = Screen('DisplaySize', monitorIndex);

    pixelsPerDegrees = monitorPixelSize(3) - monitorPixelSize(1);
    pixelsPerDegrees = pixelsPerDegrees / Utils.ComputeVisualAngleDegrees(viewingDistance, monitorMilliWidth);

    win = Screen('OpenWindow', monitorIndex, Config.screenBackgroundColor, Config.screenRect);

    screenCenter = floor((Config.screenRect(3:4) - Config.screenRect(1:2))/2);
    screenCenterDegrees = screenCenter / pixelsPerDegrees;

    dotRectDegrees = [screenCenterDegrees(1)-Config.dotRectSize(1)/2 screenCenterDegrees(2)-Config.dotRectSize(2)/2 ...
        screenCenterDegrees(1)+Config.dotRectSize(1)/2 screenCenterDegrees(2)+Config.dotRectSize(2)/2];

    dotSizePixels = Config.dotWidth * pixelsPerDegrees;
    [minSmooth, maxSmooth] = Screen('DrawDots', win);
    dotSizePixels = min(max(dotSizePixels, minSmooth), maxSmooth);
end

%% Stimulus-type configuration
% Data flow: Config.stimulusType -> per-condition parameters.
stimulusTypeConfig = struct();
if Config.stimulusType == Utils.likelihood
    stimulusTypeConfig = Config.likelihood;
elseif Config.stimulusType == Utils.path_duration
    stimulusTypeConfig = Config.path_duration;
elseif Config.stimulusType == Utils.path_duration_norm
    stimulusTypeConfig = Config.path_duration_norm;
end

%% Main generation loop
% Data flow: condition params -> trial loops -> relative paths -> feasible placement -> xySeqs outputs.
if ~skipGeneration
    framesPerTrial = round(Config.trialDuration * Config.frameFrequency);
    minPos = [Config.dotWidth/2, Config.dotWidth/2];
    maxPos = Config.dotRectSize - minPos;
    xySeqsPredicted = struct();
    generationAttemptStats = struct( ...
        'totalAttempts', 0, ...
        'acceptedTrials', 0, ...
        'rangeFailures', 0, ...
        'minDistanceFailures', 0, ...
        'drsaProxyFailures', 0, ...
        'maxAttemptsSingleTrial', 0);

    % Geometry-derived floor for constant-curvature trajectories.
    % Data flow: dot rectangle + speed -> min curvature magnitude that avoids straight-line overflow.
    maxTurnRadiusDeg = min(Config.dotRectSize - Config.dotWidth) / 2;
    if maxTurnRadiusDeg <= 0
        error('Dot rectangle minus dot width must be > 0 to generate paths.');
    end
    curvatureFeasibilityFloorDeg = rad2deg(scaledDotSpeedDegPerFrame / maxTurnRadiusDeg);
    effectiveCurvatureFloorDeg = min(curvatureFeasibilityFloorDeg, abs(Config.curvFactor));
    if enforceCurvatureFeasibilityFloor && abs(Config.curvFactor) < curvatureFeasibilityFloorDeg
        warning(['Config.curvFactor (%.4f) is below geometry feasibility floor (%.4f). ' ...
            'Using %.4f as effective floor to avoid impossible requests.'], ...
            abs(Config.curvFactor), curvatureFeasibilityFloorDeg, effectiveCurvatureFloorDeg);
    end

    for dirVarIndex = 1:length(stimulusTypeConfig.directionVariance)
        dirVar = stimulusTypeConfig.directionVariance(dirVarIndex);
        isDeviantCondition = (Config.stimulusType == Utils.likelihood && dirVar ~= 0);

        for pathDurIndex = 1:length(stimulusTypeConfig.pathDuration)
            pathDuration = stimulusTypeConfig.pathDuration(pathDurIndex);
            trialDeviances = linspace(dirVar, 360 - dirVar, Config.trialsPerCondition);
            isNondeviantCondition = (Config.stimulusType ~= Utils.likelihood || dirVar == 0);
            useDrsaProxyGate = enableDrsaProxyGate && ...
                (~applyDrsaProxyGateToNondeviantOnly || isNondeviantCondition);
            drsaProxyFrameIndices = 1:drsaProxyGateFrameStride:framesPerTrial;
            if drsaProxyFrameIndices(end) ~= framesPerTrial
                drsaProxyFrameIndices(end + 1) = framesPerTrial;
            end
            acceptedDot1Positions = zeros(Config.trialsPerCondition, framesPerTrial, 2);
            acceptedDot1Directions = zeros(Config.trialsPerCondition, framesPerTrial, 2);
            acceptedTrialCount = 0;

            totalSpread = zeros(Config.numXGrids, Config.numYGrids, 2);

            for trialPerCondIndex = 1:Config.trialsPerCondition
                % Constraint-driven resampling: regenerate when placement,
                % min-distance, or dRSA-proxy gate criteria fail.
                % For the dRSA proxy gate, only candidates that improve the
                % cross-model proxy are accepted once enough trials exist.
                trialValid = false;
                attemptsThisTrial = 0;
                currentDrsaProxyScore = NaN;
                if useDrsaProxyGate && acceptedTrialCount >= drsaProxyGateMinAcceptedTrials
                    currentDrsaProxyScore = local_compute_drsa_proxy_score( ...
                        acceptedDot1Positions(1:acceptedTrialCount, :, :), ...
                        acceptedDot1Directions(1:acceptedTrialCount, :, :), ...
                        drsaProxyFrameIndices, drsaProxyGateExcludeMainDiagLags);
                end
                while ~trialValid
                    attemptsThisTrial = attemptsThisTrial + 1;
                    if attemptsThisTrial > maxAttemptsPerTrial
                        error(['Max attempts reached (%d) in condition dirVar=%g, trial=%d. ' ...
                            'Check curvature/min-distance feasibility.'], ...
                            maxAttemptsPerTrial, dirVar, trialPerCondIndex);
                    end
                    generationAttemptStats.totalAttempts = generationAttemptStats.totalAttempts + 1;
                    trialValid = true;

                    % Per-trial parameters: initial directions and baseline curvature.
                    directionAngle = [Utils.RandAngleDegree(), Utils.RandAngleDegree()];
                    directionAngleNoDev = directionAngle;
                    curvynessFactor = Config.curvFactor * [ ...
                        Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand), ...
                        Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand)];

                    % Optional anti-bias safeguard: avoid near-zero curvature that cannot fit in-bounds.
                    if enforceCurvatureFeasibilityFloor
                        for dotIndex = 1:2
                            curvSign = sign(curvynessFactor(dotIndex));
                            if curvSign == 0
                                curvSign = Utils.RandPosNeg();
                            end
                            curvMag = abs(curvynessFactor(dotIndex));
                            if curvMag < effectiveCurvatureFloorDeg
                                if Config.isCurvFactorRand && abs(Config.curvFactor) > effectiveCurvatureFloorDeg
                                    curvMag = effectiveCurvatureFloorDeg + rand(1) * ...
                                        (abs(Config.curvFactor) - effectiveCurvatureFloorDeg);
                                else
                                    curvMag = effectiveCurvatureFloorDeg;
                                end
                                curvynessFactor(dotIndex) = curvSign * curvMag;
                            end
                        end
                    end

                    % Vectorized path generation: build per-frame angles from relative origin.
                    numSteps = framesPerTrial - 1;
                    devFrameOnset = [];
                    directionAngleChange = zeros(max(numSteps, 0), 1);
                    if numSteps > 0
                        if Config.stimulusType == Utils.likelihood
                            if dirVar ~= 0
                                devOnsetVariance = rand(1) * Utils.RandPosNeg() * ...
                                    Config.deviantOnsetVariance * Config.deviantOnset;
                                devFrameOnset = round((Config.deviantOnset + devOnsetVariance) * framesPerTrial);
                                devFrameOnset = min(max(devFrameOnset, 2), framesPerTrial);
                                directionAngleChange(devFrameOnset-1) = trialDeviances(trialPerCondIndex);
                            end
                        else
                            if Config.stimulusType == Utils.path_duration
                                directionAngleChange(1) = Utils.RandAngleDegree();
                            elseif Config.stimulusType == Utils.path_duration_norm
                                directionAngleChange(1) = normrnd(0, stimulusTypeConfig.directionChange);
                            end
                        end
                    end

                    directionAngleChange = repmat(directionAngleChange, 1, 2);
                    directionAngleChangeNoDev = directionAngleChange;
                    if Config.stimulusType == Utils.likelihood && dirVar ~= 0
                        directionAngleChangeNoDev(:) = 0; % baseline represents the no-deviant path
                    end

                    curvynessPerStep = repmat(curvynessFactor, max(numSteps, 0), 1);
                    curvynessPerStepNoDev = curvynessPerStep;
                    if ~isempty(devFrameOnset) && numSteps > 0
                        flipIndex = devFrameOnset - 1;
                        numFlipSteps = numSteps - flipIndex + 1;

                        % Deviant curvature mode precedence:
                        % randomizeCurvatureOnDeviant > flipCurvatureOnDeviant.
                        if randomizeCurvatureOnDeviant
                            deviantCurvynessFactor = (2 * rand(1, 2) - 1) * deviantCurvatureRange;
                            curvynessPerStep(flipIndex:end, :) = repmat(deviantCurvynessFactor, numFlipSteps, 1);
                        elseif flipCurvatureOnDeviant
                            curvynessPerStep(flipIndex:end, :) = repmat(-curvynessFactor, numFlipSteps, 1);
                        end
                    end

                    allPathsDirectionAngle = zeros(framesPerTrial, 2);
                    allPathsDirectionAngleNoDev = zeros(framesPerTrial, 2);
                    allPathsDirectionAngle(1, :) = directionAngle;
                    allPathsDirectionAngleNoDev(1, :) = directionAngleNoDev;
                    if numSteps > 0
                        allPathsDirectionAngle(2:end, :) = directionAngle + ...
                            cumsum(directionAngleChange + curvynessPerStep, 1);
                        allPathsDirectionAngleNoDev(2:end, :) = directionAngleNoDev + ...
                            cumsum(directionAngleChangeNoDev + curvynessPerStepNoDev, 1);
                    end

                    % Data flow: angles -> direction vectors -> relative (origin-centered) positions.
                    directionVectors = [cosd(allPathsDirectionAngle(:, 1)), sind(allPathsDirectionAngle(:, 1)), ...
                        cosd(allPathsDirectionAngle(:, 2)), sind(allPathsDirectionAngle(:, 2))];
                    dummyDirectionVectors = [cosd(allPathsDirectionAngleNoDev(:, 1)), sind(allPathsDirectionAngleNoDev(:, 1)), ...
                        cosd(allPathsDirectionAngleNoDev(:, 2)), sind(allPathsDirectionAngleNoDev(:, 2))];

                    relativePathsFrameDotXY = zeros(framesPerTrial, 4);
                    relativeDummyPathsFrameDotXY = zeros(framesPerTrial, 4);
                    if numSteps > 0
                        stepVectors = directionVectors(2:end, :) * scaledDotSpeedDegPerFrame;
                        dummyStepVectors = dummyDirectionVectors(2:end, :) * scaledDotSpeedDegPerFrame;
                        relativePathsFrameDotXY(2:end, :) = cumsum(stepVectors, 1);
                        relativeDummyPathsFrameDotXY(2:end, :) = cumsum(dummyStepVectors, 1);
                    end

                    % Boundary-safe start placement.
                    % Data flow: relative paths -> feasible start ranges -> absolute paths in-bounds.
                    dot1X = [relativePathsFrameDotXY(:, 1); relativeDummyPathsFrameDotXY(:, 1)];
                    dot1Y = [relativePathsFrameDotXY(:, 2); relativeDummyPathsFrameDotXY(:, 2)];
                    dot2X = [relativePathsFrameDotXY(:, 3); relativeDummyPathsFrameDotXY(:, 3)];
                    dot2Y = [relativePathsFrameDotXY(:, 4); relativeDummyPathsFrameDotXY(:, 4)];

                    dot1XRange = [minPos(1) - min(dot1X), maxPos(1) - max(dot1X)];
                    dot1YRange = [minPos(2) - min(dot1Y), maxPos(2) - max(dot1Y)];
                    dot2XRange = [minPos(1) - min(dot2X), maxPos(1) - max(dot2X)];
                    dot2YRange = [minPos(2) - min(dot2Y), maxPos(2) - max(dot2Y)];

                    if dot1XRange(1) > dot1XRange(2) || dot1YRange(1) > dot1YRange(2) || ...
                            dot2XRange(1) > dot2XRange(2) || dot2YRange(1) > dot2YRange(2)
                        generationAttemptStats.rangeFailures = generationAttemptStats.rangeFailures + 1;
                        trialValid = false;
                        continue;
                    end

                    startDot1 = [ ...
                        dot1XRange(1) + rand(1) * (dot1XRange(2) - dot1XRange(1)), ...
                        dot1YRange(1) + rand(1) * (dot1YRange(2) - dot1YRange(1))];
                    startDot2 = [ ...
                        dot2XRange(1) + rand(1) * (dot2XRange(2) - dot2XRange(1)), ...
                        dot2YRange(1) + rand(1) * (dot2YRange(2) - dot2YRange(1))];

                    allPathsFrameDotXY = relativePathsFrameDotXY + [startDot1 startDot2];
                    dummyPathsFrameDotXY = relativeDummyPathsFrameDotXY + [startDot1 startDot2];

                    % Min-distance constraints across full trajectories.
                    interDotDeltaObserved = allPathsFrameDotXY(:, 1:2) - allPathsFrameDotXY(:, 3:4);
                    interDotDistObserved = sqrt(sum(interDotDeltaObserved.^2, 2));
                    distTooCloseObserved = interDotDistObserved < Config.minDistanceBetweenDots;

                    distTooCloseNoDev = false(size(distTooCloseObserved));
                    if enforceMinDistanceOnNoDevBaseline
                        interDotDeltaNoDev = dummyPathsFrameDotXY(:, 1:2) - dummyPathsFrameDotXY(:, 3:4);
                        interDotDistNoDev = sqrt(sum(interDotDeltaNoDev.^2, 2));
                        distTooCloseNoDev = interDotDistNoDev < Config.minDistanceBetweenDots;
                    end
                    if any(distTooCloseObserved | distTooCloseNoDev)
                        generationAttemptStats.minDistanceFailures = generationAttemptStats.minDistanceFailures + 1;
                        trialValid = false;
                        continue;
                    end

                    % Optional dRSA-proxy-aware gate (dot1, nondeviant by default).
                    % Data flow: accepted trial bank + candidate trial -> dRSA-style
                    % position-vs-direction proxy score -> accept/reject.
                    if useDrsaProxyGate && acceptedTrialCount >= drsaProxyGateMinAcceptedTrials
                        candidateBankPositions = cat(1, ...
                            acceptedDot1Positions(1:acceptedTrialCount, :, :), ...
                            reshape(allPathsFrameDotXY(:, 1:2), [1, framesPerTrial, 2]));
                        candidateBankDirections = cat(1, ...
                            acceptedDot1Directions(1:acceptedTrialCount, :, :), ...
                            reshape(directionVectors(:, 1:2), [1, framesPerTrial, 2]));
                        candidateDrsaProxyScore = local_compute_drsa_proxy_score( ...
                            candidateBankPositions, candidateBankDirections, ...
                            drsaProxyFrameIndices, drsaProxyGateExcludeMainDiagLags);
                        if ~isfinite(candidateDrsaProxyScore)
                            candidateDrsaProxyScore = inf;
                        end

                        scoreImprovement = currentDrsaProxyScore - candidateDrsaProxyScore;
                        relativePass = ~isfinite(currentDrsaProxyScore) || ...
                            scoreImprovement >= drsaProxyGateMinScoreImprovement;
                        absolutePass = isinf(drsaProxyGateMaxMeanAbsCorr) || ...
                            candidateDrsaProxyScore <= drsaProxyGateMaxMeanAbsCorr;
                        if ~(relativePass && absolutePass)
                            generationAttemptStats.drsaProxyFailures = ...
                                generationAttemptStats.drsaProxyFailures + 1;
                            trialValid = false;
                            continue;
                        end
                    end

                    % Per-trial buffers from vectorized outputs.
                    allPathsCurvyness = zeros(framesPerTrial, 2);
                    allPathsCurvyness(1, :) = curvynessFactor;
                    if numSteps > 0
                        allPathsCurvyness(2:end, :) = curvynessPerStep;
                    end
                    allPathsStartingPoint = zeros(framesPerTrial, 2);
                    allPathsStartingPoint(1, :) = [1, 1];
                    if numSteps > 0
                        changeMask = any(directionAngleChange, 2);
                        allPathsStartingPoint(2:end, :) = repmat(changeMask, 1, 2);
                    end

                    % Analysis-only buffers for no-deviant baselines (deviant conditions only).
                    % Data flow: no-deviant angle changes -> per-frame curvyness/start flags.
                    if isDeviantCondition
                        allPathsCurvynessNoDev = zeros(framesPerTrial, 2);
                        allPathsCurvynessNoDev(1, :) = curvynessFactor;
                        if numSteps > 0
                            allPathsCurvynessNoDev(2:end, :) = curvynessPerStepNoDev;
                        end

                        allPathsStartingPointNoDev = zeros(framesPerTrial, 2);
                        allPathsStartingPointNoDev(1, :) = [1, 1];
                        if numSteps > 0
                            changeMaskNoDev = any(directionAngleChangeNoDev, 2);
                            allPathsStartingPointNoDev(2:end, :) = repmat(changeMaskNoDev, 1, 2);
                        end
                    end

                    % Grid occupancy bookkeeping and optional preview rendering.
                    gridCount = zeros(Config.numXGrids, Config.numYGrids, framesPerTrial, 2);
                    for frameIndex = 1:framesPerTrial
                        gridCount(:, :, frameIndex, 1) = Utils.GetFrameInGrid( ...
                            Config.numXGrids, Config.xGridSize, Config.numYGrids, Config.yGridSize, ...
                            allPathsFrameDotXY(frameIndex, 1:2));
                        gridCount(:, :, frameIndex, 2) = Utils.GetFrameInGrid( ...
                            Config.numXGrids, Config.xGridSize, Config.numYGrids, Config.yGridSize, ...
                            allPathsFrameDotXY(frameIndex, 3:4));

                        if renderPreview
                            % Preview drawing (optional, matches v5 behavior).
                            xyDrawDot1(1,1) = allPathsFrameDotXY(frameIndex, 1) * pixelsPerDegrees + dotRectDegrees(1) * pixelsPerDegrees;
                            xyDrawDot1(1,2) = allPathsFrameDotXY(frameIndex, 2) * pixelsPerDegrees + dotRectDegrees(2) * pixelsPerDegrees;
                            xyDrawDot2(1,1) = allPathsFrameDotXY(frameIndex, 3) * pixelsPerDegrees + dotRectDegrees(1) * pixelsPerDegrees;
                            xyDrawDot2(1,2) = allPathsFrameDotXY(frameIndex, 4) * pixelsPerDegrees + dotRectDegrees(2) * pixelsPerDegrees;

                            Screen('FillRect', win, Config.dotRectColor, dotRectDegrees * pixelsPerDegrees);
                            Screen('DrawDots', win, xyDrawDot1, dotSizePixels, Config.dot1Color*255, [0 0], 1);
                            Screen('DrawDots', win, xyDrawDot2, dotSizePixels, Config.dot2Color*255, [0 0], 1);

                            Screen('TextSize', win, 40);
                            Screen('TextFont', win, 'Arial');
                            DrawFormattedText(win, '+', 'center', 'center', Config.crossDefaultColor*255);
                            Screen('Flip', win, [], 1);
                        end
                    end
                end
                generationAttemptStats.acceptedTrials = generationAttemptStats.acceptedTrials + 1;
                generationAttemptStats.maxAttemptsSingleTrial = max( ...
                    generationAttemptStats.maxAttemptsSingleTrial, attemptsThisTrial);
                acceptedTrialCount = acceptedTrialCount + 1;
                acceptedDot1Positions(acceptedTrialCount, :, :) = allPathsFrameDotXY(:, 1:2);
                acceptedDot1Directions(acceptedTrialCount, :, :) = directionVectors(:, 1:2);

                % Store trial outputs (match v5 struct layout).
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).condition          = dirVar;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).PredictionRange    = pathDuration;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).sequence           = trialPerCondIndex;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).xy                 = allPathsFrameDotXY;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).pathAll            = allPathsStartingPoint;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).curvyness          = allPathsCurvyness;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).AngleDirection     = allPathsDirectionAngle;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).gridCounter        = gridCount;
                gridSum = squeeze(sum(gridCount, 3));
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).gridTot            = gridSum;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).tolleranceGrid     = Config.spreadoutTolerance;
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).borderCounter      = zeros(1, framesPerTrial);
                xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).borderTot          = 0;

                % Analysis-only storage: no-deviant baseline for deviant conditions.
                % Data flow: dummy path -> xySeqsPredicted output struct.
                if isDeviantCondition
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).condition       = dirVar;
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).PredictionRange = pathDuration;
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).sequence        = trialPerCondIndex;
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).xy              = dummyPathsFrameDotXY;
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).pathAll         = allPathsStartingPointNoDev;
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).curvyness       = allPathsCurvynessNoDev;
                    xySeqsPredicted(dirVarIndex, pathDurIndex, trialPerCondIndex).AngleDirection  = allPathsDirectionAngleNoDev;
                end

                totalSpread = totalSpread + gridSum;
                fprintf("End of trial %d within condition %d\n", trialPerCondIndex, dirVar);
            end

            xySeqs(dirVarIndex, pathDurIndex).totalSpread = totalSpread;
            fprintf("End of condition [%d] = %d\n", dirVarIndex, dirVar);
        end
    end

    if generationAttemptStats.acceptedTrials > 0
        generationAttemptStats.meanAttemptsPerAcceptedTrial = ...
            generationAttemptStats.totalAttempts / generationAttemptStats.acceptedTrials;
    else
        generationAttemptStats.meanAttemptsPerAcceptedTrial = NaN;
    end
end

%% Optional path plots
% Data flow: xySeqs -> sampled trials -> time-graded scatter in 2x2 grid.
if plotPathsAtEnd
    rng('shuffle');
    plotStream = RandStream.getGlobalStream;
    axisLimits = [0 Config.dotRectSize(1) 0 Config.dotRectSize(2)];
    for pathDurIndex = 1:length(stimulusTypeConfig.pathDuration)
        pathDuration = stimulusTypeConfig.pathDuration(pathDurIndex);
        figure('Color', [1 1 1]);

        for condIndex = 1:length(plotPathConditions)
            condValue = plotPathConditions(condIndex);
            dirVarIndex = find(stimulusTypeConfig.directionVariance == condValue, 1);
            for dotIndex = 1:2
                subplotIndex = (condIndex - 1) * 2 + dotIndex;
                ax = subplot(2, 2, subplotIndex);
                hold(ax, 'on');

                if isempty(dirVarIndex)
                    text(ax, 0.5, 0.5, 'Condition not found', 'HorizontalAlignment', 'center');
                    axis(ax, 'off');
                    hold(ax, 'off');
                    continue;
                end

                availablePaths = size(xySeqs, 3);
                numToPlot = min(plotPathsPerCondition, availablePaths);
                if numToPlot == 0
                    text(ax, 0.5, 0.5, 'No paths available', 'HorizontalAlignment', 'center');
                    axis(ax, 'off');
                    hold(ax, 'off');
                    continue;
                end
                sampleIdx = randperm(plotStream, availablePaths, numToPlot);

                % Data flow: selected trials -> [x y time] -> time-colored scatter.
                dotPaths = zeros(numToPlot, 2, framesPerTrial);
                for sampleIndex = 1:numToPlot
                    xy = xySeqs(dirVarIndex, pathDurIndex, sampleIdx(sampleIndex)).xy;
                    if dotIndex == 1
                        dotPaths(sampleIndex, 1, :) = xy(:, 1);
                        dotPaths(sampleIndex, 2, :) = xy(:, 2);
                    else
                        dotPaths(sampleIndex, 1, :) = xy(:, 3);
                        dotPaths(sampleIndex, 2, :) = xy(:, 4);
                    end
                end

                xVals = reshape(dotPaths(:, 1, :), [], 1);
                yVals = reshape(dotPaths(:, 2, :), [], 1);
                timeIdx = repmat(1:framesPerTrial, numToPlot, 1);
                timeIdx = timeIdx(:);
                cmap = parula(framesPerTrial);
                pointColors = cmap(timeIdx, :);

                scatter(ax, xVals, yVals, 8, pointColors, 'filled');
                colormap(ax, cmap);
                caxis(ax, [1 framesPerTrial]);
                cb = colorbar(ax);
                cb.Label.String = 'Time (frames)';

                axis(ax, 'equal');
                axis(ax, axisLimits);
                xlabel(ax, 'x (deg)');
                ylabel(ax, 'y (deg)');
                title(ax, sprintf('dirVar=%g, dot=%d', condValue, dotIndex));
                hold(ax, 'off');
            end
        end
        if exist('sgtitle', 'file') == 2
            sgtitle(sprintf('pathDuration=%g (n=%d)', pathDuration, min(plotPathsPerCondition, size(xySeqs, 3))));
        end
    end
end

%% Save outputs (xySeqs + Cfg, plus no-deviant baselines)
% Data flow: Config + script settings + RNG seed -> repro + Cfg -> .mat on disk.
if ~skipSave
    %% Reproducibility snapshot
    % Data flow: Config constants + script params + inputs -> repro struct.
    configProps = properties('Config');
    configSnapshot = struct();
    for propIndex = 1:numel(configProps)
        propName = configProps{propIndex};
        configSnapshot.(propName) = Config.(propName); % capture constant Config values
    end

    repro = struct();
    repro.script = struct( ...
        'name', mfilename, ...
        'version', 'v15_experimentalPathScale', ...
        'parameters', struct( ...
            'renderPreview', renderPreview, ...
            'plotPathsAtEnd', plotPathsAtEnd, ...
            'plotPathsPerCondition', plotPathsPerCondition, ...
            'plotPathConditions', plotPathConditions, ...
            'pathScale', pathScale, ...
            'scaledDotSpeedDegPerFrame', scaledDotSpeedDegPerFrame, ...
            'flipCurvatureOnDeviant', flipCurvatureOnDeviant, ...
            'randomizeCurvatureOnDeviant', randomizeCurvatureOnDeviant, ...
            'deviantCurvatureRange', deviantCurvatureRange, ...
            'enforceCurvatureFeasibilityFloor', enforceCurvatureFeasibilityFloor, ...
            'enforceMinDistanceOnNoDevBaseline', enforceMinDistanceOnNoDevBaseline, ...
            'enableDrsaProxyGate', enableDrsaProxyGate, ...
            'applyDrsaProxyGateToNondeviantOnly', applyDrsaProxyGateToNondeviantOnly, ...
            'drsaProxyGateMinAcceptedTrials', drsaProxyGateMinAcceptedTrials, ...
            'drsaProxyGateFrameStride', drsaProxyGateFrameStride, ...
            'drsaProxyGateExcludeMainDiagLags', drsaProxyGateExcludeMainDiagLags, ...
            'drsaProxyGateMinScoreImprovement', drsaProxyGateMinScoreImprovement, ...
            'drsaProxyGateMaxMeanAbsCorr', drsaProxyGateMaxMeanAbsCorr, ...
            'maxAttemptsPerTrial', maxAttemptsPerTrial));
    repro.inputs = struct('subjectID', subjectID, 'viewingDistance', viewingDistance);
    repro.rng = rngState;
    repro.config = configSnapshot;
    repro.stimulusTypeConfig = stimulusTypeConfig;
    savedDrsaProxyFrameIndices = 1:drsaProxyGateFrameStride:framesPerTrial;
    if savedDrsaProxyFrameIndices(end) ~= framesPerTrial
        savedDrsaProxyFrameIndices(end + 1) = framesPerTrial;
    end
    repro.derived = struct( ...
        'framesPerTrial', framesPerTrial, ...
        'minPos', minPos, ...
        'maxPos', maxPos, ...
        'pathScale', pathScale, ...
        'scaledDotSpeedDegPerFrame', scaledDotSpeedDegPerFrame, ...
        'drsaProxyFrameIndices', savedDrsaProxyFrameIndices, ...
        'curvatureFeasibilityFloorDeg', curvatureFeasibilityFloorDeg, ...
        'effectiveCurvatureFloorDeg', effectiveCurvatureFloorDeg, ...
        'generationAttemptStats', generationAttemptStats);

    Cfg = struct();
    Cfg.dpf = scaledDotSpeedDegPerFrame;
    Cfg.pathScale = pathScale;
    Cfg.fps = Config.frameFrequency;
    Cfg.Stimulitype = Config.stimulusType;
    Cfg.dot_w = Config.dotWidth;
    Cfg.rectSize = Config.dotRectSize;
    Cfg.DirChange = stimulusTypeConfig.directionChange;

    xySeqs = repmat(xySeqs, 1, Config.trialRepetetion, 1);
    save([Config.inputDirectory sprintf(Config.stimuliFileName, subjectID)], 'xySeqs', 'Cfg', 'repro');

    % Analysis-only output: no-deviant baselines for deviant conditions.
    xySeqsPredicted = repmat(xySeqsPredicted, 1, Config.trialRepetetion, 1);
    save(predictedFile, 'xySeqsPredicted', 'Cfg');
end

% Close all onscreens and offscreens
if renderPreview
    sca;
end

%% Local helpers (dRSA-proxy scoring)
function meanAbsCorrScore = local_compute_drsa_proxy_score( ...
        dot1Positions, dot1Directions, frameIndices, excludeDiagLagBins)
% LOCAL_COMPUTE_DRSA_PROXY_SCORE Compute a dRSA-style proxy score.
%
% Purpose:
%   Compute the acceptance-gate proxy used in v14/v15 by mirroring the
%   position-dot1 vs direction-dot1 corr dRSA branch at reduced time
%   resolution:
%     - position RDM columns use euclidean pair distances across trials,
%     - direction RDM columns use cosine pair distances across trials,
%     - score is mean(abs(corr(directionRDM, positionRDM))).
%
% Example usage:
%   score = local_compute_drsa_proxy_score( ...
%       dot1PosBank, dot1DirBank, 1:4:320, 1);
%
% Inputs:
%   - dot1Positions: [nTrials x nFrames x 2] dot1 observed positions.
%   - dot1Directions: [nTrials x nFrames x 2] dot1 direction unit vectors.
%   - frameIndices: sampled frames used to build proxy RDM columns.
%   - excludeDiagLagBins: ignore |lag| <= this many sampled bins in score.
%
% Output:
%   - meanAbsCorrScore: scalar proxy score (lower is better).
%
% Data flow:
%   sampled frames -> position/direction RDM columns -> cross-time columnwise
%   correlation matrix -> abs + lag mask average -> scalar proxy score.
    meanAbsCorrScore = 0;
    nTrials = size(dot1Positions, 1);
    if nTrials < 3 || isempty(frameIndices)
        return;
    end

    nFramesToScore = numel(frameIndices);
    trialPairs = nchoosek(1:nTrials, 2);
    pairI = trialPairs(:, 1);
    pairJ = trialPairs(:, 2);
    nPairs = size(trialPairs, 1);
    rdmPosition = zeros(nPairs, nFramesToScore);
    rdmDirection = zeros(nPairs, nFramesToScore);

    for frameListIndex = 1:nFramesToScore
        frameIndex = frameIndices(frameListIndex);
        posFrame = squeeze(dot1Positions(:, frameIndex, :));
        dirFrame = squeeze(dot1Directions(:, frameIndex, :));
        rdmPosition(:, frameListIndex) = local_pairwise_euclidean_from_pairs(posFrame, pairI, pairJ);
        rdmDirection(:, frameListIndex) = local_pairwise_cosine_from_pairs(dirFrame, pairI, pairJ);
    end

    corrMatrix = local_columnwise_correlation(rdmDirection, rdmPosition);
    if excludeDiagLagBins > 0
        [rowIdx, colIdx] = ndgrid(1:nFramesToScore, 1:nFramesToScore);
        validMask = abs(rowIdx - colIdx) > excludeDiagLagBins;
    else
        validMask = true(size(corrMatrix));
    end
    scoreValues = abs(corrMatrix(validMask));
    scoreValues = scoreValues(~isnan(scoreValues));
    if ~isempty(scoreValues)
        meanAbsCorrScore = mean(scoreValues);
    end
end

function distValues = local_pairwise_euclidean_from_pairs(pointsXY, pairI, pairJ)
% LOCAL_PAIRWISE_EUCLIDEAN_FROM_PAIRS Euclidean distances for selected trial pairs.
%
% Inputs:
%   - pointsXY: [nTrials x 2] coordinates for one frame.
%   - pairI/pairJ: pair index vectors in pdist condensed order.
%
% Output:
%   - distValues: [nPairs x 1] euclidean distances.
    deltas = pointsXY(pairI, :) - pointsXY(pairJ, :);
    distValues = sqrt(sum(deltas.^2, 2));
end

function distValues = local_pairwise_cosine_from_pairs(vectorsXY, pairI, pairJ)
% LOCAL_PAIRWISE_COSINE_FROM_PAIRS Cosine distances for selected trial pairs.
%
% Inputs:
%   - vectorsXY: [nTrials x 2] vectors for one frame.
%   - pairI/pairJ: pair index vectors in pdist condensed order.
%
% Output:
%   - distValues: [nPairs x 1] cosine distances (1 - cosine similarity).
    vecA = vectorsXY(pairI, :);
    vecB = vectorsXY(pairJ, :);
    dotAB = sum(vecA .* vecB, 2);
    normA = sqrt(sum(vecA.^2, 2));
    normB = sqrt(sum(vecB.^2, 2));
    denom = normA .* normB;
    cosSim = zeros(size(dotAB));
    validMask = denom > eps;
    cosSim(validMask) = dotAB(validMask) ./ denom(validMask);
    cosSim = max(min(cosSim, 1), -1);
    distValues = 1 - cosSim;
end

function corrMatrix = local_columnwise_correlation(matrixX, matrixY)
% LOCAL_COLUMNWISE_CORRELATION Compute corr(X, Y) without toolbox dependency.
%
% Inputs:
%   - matrixX: [nObs x nColsX]
%   - matrixY: [nObs x nColsY]
%
% Output:
%   - corrMatrix: [nColsX x nColsY] columnwise Pearson correlations.
    centeredX = bsxfun(@minus, matrixX, mean(matrixX, 1));
    centeredY = bsxfun(@minus, matrixY, mean(matrixY, 1));
    normX = sqrt(sum(centeredX.^2, 1));
    normY = sqrt(sum(centeredY.^2, 1));
    denom = normX' * normY;
    corrMatrix = (centeredX' * centeredY) ./ denom;
    corrMatrix(denom <= eps) = NaN;
end
