%% Stimuli Generation (deviant displacement + turn/curvature controls, boundary-safe, dRSA-proxy gate)
% Script: stimuli_generation_v16_Displacement.m
% Author: Marisa (original), Ayman Hatoum (v5), updated by Dami (v6-v8), v9-v16 updates by Dami
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
%   v16 keeps curvature smooth and predictable by design:
%     - one constant curvature value per dot at trial start,
%     - optional deviant-only curvature change at deviant onset,
%     - no additional within-trial curvature updates.
%   In v16, initial and post-deviant curvature can be constrained with
%   explicit windows from Config:
%     - Config.initialCurvatureWindows
%     - Config.deviantCurvatureWindows
%   Each is an [N x 2] interval list in deg/frame, so curvature sampling can
%   enforce signed thresholds directly (for example:
%   [-0.8, -0.3755] U [0.3755, 0.8]).
%   Deviant turn angles can be configured in two ways:
%     - legacy: generated from directionVariance as linspace(dirVar, 360-dirVar),
%     - explicit signed windows from Config.deviantSignedTurnWindows
%       (e.g., [-180 -45; 45 180]) for direct threshold-style control.
%
%   v16 adds an optional deviant displacement layer, applied on top of the
%   existing turn/curvature deviant logic:
%     - at deviant onset, the observed path can jump to a new post-onset
%       start point sampled in an annular-sector region (relative to local
%       pre-deviant heading),
%     - post-onset frames preserve the generated shape and are translated by
%       that sampled displacement vector.
%   This uses two primary controls:
%     - Config.likelihood.deviantDisplacementRadiusRangeDeg = [innerRadius, outerRadius]
%     - Config.likelihood.deviantDisplacementAngleWindowsDeg = [angleMin angleMax; ...]
%       in signed degrees (e.g., [-140 -40; 40 140]).
%   A displacement mode toggle can override/disable displacement geometry:
%     - Config.likelihood.deviantDisplacementMode = 'off' | 'constrained' | 'freeStart'
%       where:
%         'off' disables displacement,
%         'constrained' uses configured radius + angle windows,
%         'freeStart' ignores configured radius + angle windows and forces
%         effective controls to full-circle angles [-180, 0; 0, 180] and
%         a screen-scale radius [0, max(Config.dotRectSize)].
%   Turn and curvature deviant manipulations remain available and can be
%   combined with displacement or disabled independently through existing
%   Config settings.
%
%   To reduce residual dRSA position-direction cross-correlation without
%   changing trajectory smoothness, v16 keeps the v14 dRSA-proxy-aware gate:
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
%   stimuli_generation_v16_Displacement;
%   % Follow the dialog prompts for viewing distance and subject ID.
%
% Example usage (custom working directory):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   addpath('experiment');
%   stimuli_generation_v16_Displacement;
%
% Inputs:
%   - Config (class/struct on MATLAB path) with screen, dot, and timing params.
%   - Utils helper class on MATLAB path.
%   - Dialog values for viewing distance and subject ID.
%
% Outputs:
%   - Saves xySeqs, Cfg, and repro to Config.inputDirectory using Config.stimuliFileName.
%   - repro captures Config, script settings, and RNG state needed to reproduce paths.
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
%   - Curvature stays constant within each trial except optional deviant-point
%     modulation in deviant conditions.
%   - Deviant displacement, when enabled, shifts only observed-path frames
%     from deviant onset onward; no-deviant baselines stay unshifted.
%   - Displacement is configured via Config.likelihood and can be toggled
%     by deviantDisplacementMode ('off'|'constrained'|'freeStart').
%   - In 'constrained', displacement can also be disabled by setting
%     radius to [0, 0] or angle windows to [].
%   - In 'freeStart', configured displacement radius/angle windows are
%     ignored and replaced by full-circle controls.
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
flipCurvatureOnDeviant = Config.flipCurvatureOnDeviant; % set true to flip curvature sign at deviant onset
randomizeCurvatureOnDeviant = Config.randomizeCurvatureOnDeviant; % set true to sample a new curvature at deviant onset
deviantCurvatureRange = Config.deviantCurvatureRange; % sampled post-deviant curvature range is [-range, +range]
% Explicit curvature windows for v16.
initialCurvatureWindows = local_parse_interval_windows(Config.initialCurvatureWindows, ...
    'Config.initialCurvatureWindows', false, -inf, inf);
deviantCurvatureWindows = local_parse_interval_windows(Config.deviantCurvatureWindows, ...
    'Config.deviantCurvatureWindows', false, -inf, inf);
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

% Deviant displacement controls for v16.
% Data flow: Config.likelihood settings -> validated controls -> effective mode-specific controls.
% Minimum control surface:
%   0) displacement mode ('off'|'constrained'|'freeStart').
%   1) radius range of the annulus [inner, outer] in visual degrees.
%   2) signed angle windows (degrees) relative to pre-deviant heading.
deviantDisplacementModeRaw = local_get_struct_field_or_default( ...
    stimulusTypeConfig, 'deviantDisplacementMode', ...
    local_get_struct_field_or_default(Config.likelihood, ...
    'deviantDisplacementMode', 'constrained'));
deviantDisplacementMode = local_parse_deviant_displacement_mode( ...
    deviantDisplacementModeRaw, 'deviantDisplacementMode');
deviantDisplacementRadiusRangeRaw = local_get_struct_field_or_default( ...
    stimulusTypeConfig, 'deviantDisplacementRadiusRangeDeg', ...
    local_get_struct_field_or_default(Config.likelihood, ...
    'deviantDisplacementRadiusRangeDeg', [0, 0]));
deviantDisplacementAngleWindowsRaw = local_get_struct_field_or_default( ...
    stimulusTypeConfig, 'deviantDisplacementAngleWindowsDeg', ...
    local_get_struct_field_or_default(Config.likelihood, ...
    'deviantDisplacementAngleWindowsDeg', []));
deviantDisplacementRadiusRangeDegConfigured = local_parse_radius_range( ...
    deviantDisplacementRadiusRangeRaw, 'deviantDisplacementRadiusRangeDeg');
deviantDisplacementAngleWindowsDegConfigured = local_parse_interval_windows( ...
    deviantDisplacementAngleWindowsRaw, ...
    'deviantDisplacementAngleWindowsDeg', true, -180, 180);
deviantDisplacementRadiusRangeDeg = deviantDisplacementRadiusRangeDegConfigured;
deviantDisplacementAngleWindowsDeg = deviantDisplacementAngleWindowsDegConfigured;
if strcmp(deviantDisplacementMode, 'off')
    % Explicit mode override: disable displacement regardless of geometry.
    deviantDisplacementRadiusRangeDeg = [0, 0];
elseif strcmp(deviantDisplacementMode, 'freeStart')
    % Radius/angle-unconstrained mode with feasible board-scale radius.
    deviantDisplacementRadiusRangeDeg = [0, max(Config.dotRectSize)];
    deviantDisplacementAngleWindowsDeg = [-180, 0; 0, 180];
end
enableDeviantDisplacement = local_is_deviant_displacement_enabled( ...
    deviantDisplacementRadiusRangeDeg, deviantDisplacementAngleWindowsDeg);

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
    curvatureFeasibilityFloorDeg = rad2deg(Config.dotSpeedDegPerFrame / maxTurnRadiusDeg);
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
            % Build per-trial deviant turn schedule (degrees).
            % Data flow: Config likelihood settings -> turn vector -> deviant frame angle update.
            trialDeviances = local_build_likelihood_deviant_turns(stimulusTypeConfig, ...
                dirVar, Config.trialsPerCondition);
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
                    % Baseline curvature sampling for v16:
                    % Data flow: Config.initialCurvatureWindows -> per-dot constant curvature.
                    curvynessFactor = local_sample_from_windows(initialCurvatureWindows, 2);

                    % Optional anti-bias safeguard: avoid near-zero curvature that cannot fit in-bounds.
                    % Data flow: sampled baseline curvature -> floor check -> floor-safe curvature.
                    if enforceCurvatureFeasibilityFloor
                        for dotIndex = 1:2
                            curvSign = sign(curvynessFactor(dotIndex));
                            if curvSign == 0
                                curvSign = Utils.RandPosNeg();
                            end
                            curvMag = abs(curvynessFactor(dotIndex));
                            if curvMag < effectiveCurvatureFloorDeg
                                curvMag = effectiveCurvatureFloorDeg;
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
                            % Post-deviant curvature sampling for v16:
                            % Data flow: Config.deviantCurvatureWindows -> deviant-path curvature.
                            deviantCurvynessFactor = local_sample_from_windows(deviantCurvatureWindows, 2);
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
                        stepVectors = directionVectors(2:end, :) * Config.dotSpeedDegPerFrame;
                        dummyStepVectors = dummyDirectionVectors(2:end, :) * Config.dotSpeedDegPerFrame;
                        relativePathsFrameDotXY(2:end, :) = cumsum(stepVectors, 1);
                        relativeDummyPathsFrameDotXY(2:end, :) = cumsum(dummyStepVectors, 1);
                    end

                    % Optional deviant displacement (v16): translate observed-path
                    % frames from deviant onset onward by an annular-sector sample.
                    % Data flow: dev onset + pre-deviant heading + annular-sector
                    % params -> per-dot displacement -> shifted observed relative path.
                    displacementOnsetMask = false(max(numSteps, 0), 1);
                    if enableDeviantDisplacement && ~isempty(devFrameOnset) && numSteps > 0
                        postDeviantFrameCount = framesPerTrial - devFrameOnset + 1;
                        preDeviantHeadingDeg = allPathsDirectionAngle(devFrameOnset - 1, :);
                        displacementDot1 = local_sample_deviant_displacement( ...
                            deviantDisplacementRadiusRangeDeg, ...
                            deviantDisplacementAngleWindowsDeg, ...
                            preDeviantHeadingDeg(1));
                        displacementDot2 = local_sample_deviant_displacement( ...
                            deviantDisplacementRadiusRangeDeg, ...
                            deviantDisplacementAngleWindowsDeg, ...
                            preDeviantHeadingDeg(2));

                        relativePathsFrameDotXY(devFrameOnset:end, 1:2) = ...
                            relativePathsFrameDotXY(devFrameOnset:end, 1:2) + ...
                            repmat(displacementDot1, postDeviantFrameCount, 1);
                        relativePathsFrameDotXY(devFrameOnset:end, 3:4) = ...
                            relativePathsFrameDotXY(devFrameOnset:end, 3:4) + ...
                            repmat(displacementDot2, postDeviantFrameCount, 1);
                        displacementOnsetMask(devFrameOnset - 1) = true;
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
                        % Mark deviant onset when either an angle-change or a
                        % displacement jump is applied on the observed path.
                        changeMask = any(directionAngleChange, 2) | displacementOnsetMask;
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
        'version', 'v16_Displacement', ...
        'parameters', struct( ...
            'renderPreview', renderPreview, ...
            'plotPathsAtEnd', plotPathsAtEnd, ...
            'plotPathsPerCondition', plotPathsPerCondition, ...
            'plotPathConditions', plotPathConditions, ...
            'initialCurvatureWindows', initialCurvatureWindows, ...
            'deviantCurvatureWindows', deviantCurvatureWindows, ...
            'flipCurvatureOnDeviant', flipCurvatureOnDeviant, ...
            'randomizeCurvatureOnDeviant', randomizeCurvatureOnDeviant, ...
            'deviantCurvatureRange', deviantCurvatureRange, ...
            'deviantDisplacementMode', deviantDisplacementMode, ...
            'deviantDisplacementRadiusRangeDegConfigured', deviantDisplacementRadiusRangeDegConfigured, ...
            'deviantDisplacementAngleWindowsDegConfigured', deviantDisplacementAngleWindowsDegConfigured, ...
            'deviantDisplacementRadiusRangeDeg', deviantDisplacementRadiusRangeDeg, ...
            'deviantDisplacementAngleWindowsDeg', deviantDisplacementAngleWindowsDeg, ...
            'enableDeviantDisplacement', enableDeviantDisplacement, ...
            'deviantSignedTurnWindows', local_get_struct_field_or_default(stimulusTypeConfig, ...
                'deviantSignedTurnWindows', []), ...
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
        'drsaProxyFrameIndices', savedDrsaProxyFrameIndices, ...
        'curvatureFeasibilityFloorDeg', curvatureFeasibilityFloorDeg, ...
        'effectiveCurvatureFloorDeg', effectiveCurvatureFloorDeg, ...
        'generationAttemptStats', generationAttemptStats);

    Cfg = struct();
    Cfg.dpf = Config.dotSpeedDegPerFrame;
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

%% Local helpers (interval-window sampling)
function mode = local_parse_deviant_displacement_mode(rawMode, variableName)
% LOCAL_PARSE_DEVIANT_DISPLACEMENT_MODE Validate deviant displacement mode.
%
% Purpose:
%   Parse displacement mode while accepting mixed-case user strings and
%   canonicalizing to one of the supported values.
%
% Example usage:
%   mode = local_parse_deviant_displacement_mode('freeStart', ...
%       'deviantDisplacementMode');
%
% Inputs:
%   - rawMode: char/string scalar in {'off','constrained','freeStart'}.
%   - variableName: label used in error messages.
%
% Output:
%   - mode: canonical mode string ('off'|'constrained'|'freeStart').
%
% Assumptions:
%   - mode names are case-insensitive at input.
    if isstring(rawMode)
        if numel(rawMode) ~= 1
            error('%s must be a scalar string/char value.', variableName);
        end
        rawMode = char(rawMode);
    end
    if ~ischar(rawMode)
        error('%s must be a string/char value.', variableName);
    end

    normalizedMode = lower(strtrim(rawMode));
    if strcmp(normalizedMode, 'off')
        mode = 'off';
    elseif strcmp(normalizedMode, 'constrained')
        mode = 'constrained';
    elseif strcmp(normalizedMode, 'freestart')
        mode = 'freeStart';
    else
        error(['%s must be one of: ''off'', ''constrained'', ''freeStart''.'], ...
            variableName);
    end
end

function radiusRange = local_parse_radius_range(rawRange, variableName)
% LOCAL_PARSE_RADIUS_RANGE Validate [innerRadius, outerRadius] annulus bounds.
%
% Purpose:
%   Validate the deviant displacement annulus radius controls.
%
% Example usage:
%   r = local_parse_radius_range([0.45, 1.25], ...
%       'deviantDisplacementRadiusRangeDeg');
%
% Inputs:
%   - rawRange: numeric vector with exactly two values [inner, outer].
%   - variableName: label used in error messages.
%
% Output:
%   - radiusRange: [1 x 2] validated range in degrees.
%
% Assumptions:
%   - 0 <= inner <= outer.
    if ~isnumeric(rawRange) || numel(rawRange) ~= 2
        error('%s must be a numeric [inner, outer] vector.', variableName);
    end

    radiusRange = reshape(rawRange, 1, 2);
    if any(~isfinite(radiusRange))
        error('%s contains non-finite values.', variableName);
    end
    if any(radiusRange < 0)
        error('%s values must be >= 0.', variableName);
    end
    if radiusRange(1) > radiusRange(2)
        error('%s must satisfy inner <= outer.', variableName);
    end
end

function enabled = local_is_deviant_displacement_enabled(radiusRange, angleWindows)
% LOCAL_IS_DEVIANT_DISPLACEMENT_ENABLED Decide if deviant displacement is active.
%
% Data flow:
%   radius and angle controls -> activation boolean.
    enabled = ~isempty(angleWindows) && radiusRange(2) > 0;
end

function parsedWindows = local_parse_interval_windows(rawWindows, variableName, allowEmpty, minBound, maxBound)
% LOCAL_PARSE_INTERVAL_WINDOWS Validate and normalize numeric interval windows.
%
% Purpose:
%   Normalize user-configured windows to an [N x 2] numeric matrix and
%   validate interval semantics for curvature/angle sampling.
%
% Example usage:
%   w = local_parse_interval_windows(Config.initialCurvatureWindows, ...
%       'Config.initialCurvatureWindows', false, -inf, inf);
%
% Inputs:
%   - rawWindows: [N x 2] matrix or vector [min1 max1 min2 max2 ...].
%   - variableName: label used in error messages.
%   - allowEmpty: true allows [] (returns []).
%   - minBound/maxBound: optional global bounds for interval endpoints.
%
% Output:
%   - parsedWindows: sorted [N x 2] non-overlapping intervals.
    if isempty(rawWindows)
        if allowEmpty
            parsedWindows = [];
            return;
        end
        error('%s must not be empty.', variableName);
    end
    if ~isnumeric(rawWindows)
        error('%s must be numeric.', variableName);
    end

    if isvector(rawWindows)
        if mod(numel(rawWindows), 2) ~= 0
            error('%s vector form requires an even number of values.', variableName);
        end
        parsedWindows = reshape(rawWindows, 2, [])';
    elseif size(rawWindows, 2) == 2
        parsedWindows = rawWindows;
    else
        error('%s must be Nx2 or a vector of paired bounds.', variableName);
    end

    if any(~isfinite(parsedWindows(:)))
        error('%s contains non-finite values.', variableName);
    end
    if any(parsedWindows(:, 1) >= parsedWindows(:, 2))
        error('%s requires [min, max] with min < max for each row.', variableName);
    end
    if any(parsedWindows(:) < minBound | parsedWindows(:) > maxBound)
        error('%s bounds must stay within [%g, %g].', variableName, minBound, maxBound);
    end

    parsedWindows = sortrows(parsedWindows, 1);
    if size(parsedWindows, 1) > 1
        if any(parsedWindows(2:end, 1) < parsedWindows(1:end-1, 2))
            error('%s rows must not overlap.', variableName);
        end
    end
end

function sampledValues = local_sample_from_windows(intervalWindows, nValues)
% LOCAL_SAMPLE_FROM_WINDOWS Sample uniformly across a union of intervals.
%
% Purpose:
%   Draw nValues samples from the union of disjoint intervals, with sampling
%   probability proportional to interval width.
%
% Example usage:
%   kappa = local_sample_from_windows([-0.8 -0.3755; 0.3755 0.8], 2);
%
% Inputs:
%   - intervalWindows: [N x 2] sorted disjoint [min, max] intervals.
%   - nValues: number of values to draw.
%
% Output:
%   - sampledValues: [1 x nValues] sampled values.
    if nValues < 1 || floor(nValues) ~= nValues
        error('nValues must be a positive integer.');
    end
    if isempty(intervalWindows)
        error('intervalWindows must not be empty.');
    end

    intervalWidths = intervalWindows(:, 2) - intervalWindows(:, 1);
    totalWidth = sum(intervalWidths);
    if totalWidth <= 0
        error('intervalWindows must have positive total width.');
    end

    cumulativeWidths = [0; cumsum(intervalWidths)];
    sampledValues = zeros(1, nValues);
    for sampleIndex = 1:nValues
        draw = rand(1) * totalWidth;
        windowIndex = find(draw <= cumulativeWidths(2:end), 1, 'first');
        if isempty(windowIndex)
            windowIndex = size(intervalWindows, 1);
        end
        offset = draw - cumulativeWidths(windowIndex);
        sampledValues(sampleIndex) = intervalWindows(windowIndex, 1) + offset;
    end
end

function displacementXY = local_sample_deviant_displacement(radiusRangeDeg, angleWindowsDeg, referenceHeadingDeg)
% LOCAL_SAMPLE_DEVIANT_DISPLACEMENT Sample one displacement vector in an annular sector.
%
% Purpose:
%   Draw a post-deviant displacement vector using annular-sector controls.
%   Angle windows are signed and interpreted relative to reference heading.
%
% Example usage:
%   d = local_sample_deviant_displacement([0.45 1.25], [-140 -40; 40 140], 90);
%
% Inputs:
%   - radiusRangeDeg: [inner, outer] annulus radii in degrees.
%   - angleWindowsDeg: [N x 2] signed relative angle windows in degrees.
%   - referenceHeadingDeg: heading used as 0 deg for relative windows.
%
% Output:
%   - displacementXY: [1 x 2] translation vector [dx, dy] in visual degrees.
%
% Data flow:
%   radius-range + angle windows + local heading -> sampled polar offset -> XY vector.
    if isempty(angleWindowsDeg)
        error('angleWindowsDeg must not be empty when displacement sampling is enabled.');
    end

    innerRadius = radiusRangeDeg(1);
    outerRadius = radiusRangeDeg(2);
    if outerRadius == innerRadius
        sampledRadius = outerRadius;
    else
        % Uniform area density in the annulus (not uniform in radius).
        sampledRadius = sqrt(rand(1) * (outerRadius^2 - innerRadius^2) + innerRadius^2);
    end

    sampledRelativeAngleDeg = local_sample_from_windows(angleWindowsDeg, 1);
    absoluteAngleDeg = referenceHeadingDeg + sampledRelativeAngleDeg;
    displacementXY = sampledRadius * [cosd(absoluteAngleDeg), sind(absoluteAngleDeg)];
end

%% Local helpers (deviant turn scheduling)
function trialDeviances = local_build_likelihood_deviant_turns(stimulusTypeConfig, ...
        dirVar, nTrials)
% LOCAL_BUILD_LIKELIHOOD_DEVIANT_TURNS Build per-trial deviant turn angles.
%
% Purpose:
%   Build the deviant turn values (degrees) used at the deviant onset frame
%   for likelihood stimuli. This helper supports two behaviors:
%     1) Legacy behavior (default): linspace(dirVar, 360-dirVar, nTrials).
%     2) Explicit signed-turn windows from
%        stimulusTypeConfig.deviantSignedTurnWindows, e.g.:
%           [-180 -45; 45 180]
%        for threshold-style control in signed degrees.
%
% Example usage:
%   turns = local_build_likelihood_deviant_turns(stimulusTypeConfig, 45, 50);
%
% Inputs:
%   - stimulusTypeConfig: likelihood config struct.
%   - dirVar: current directionVariance condition value.
%   - nTrials: number of trials to generate in this condition.
%
% Output:
%   - trialDeviances: [1 x nTrials] turn angles in degrees.
%
% Data flow:
%   Config likelihood fields -> validated windows (optional) ->
%   deterministic coverage across allowed ranges -> per-trial turn vector.
    if nTrials < 1 || floor(nTrials) ~= nTrials
        error('nTrials must be a positive integer.');
    end

    trialDeviances = zeros(1, nTrials);

    % Nondeviant conditions never inject a turn at deviant onset.
    if dirVar == 0
        return;
    end

    signedTurnWindows = local_get_struct_field_or_default(stimulusTypeConfig, ...
        'deviantSignedTurnWindows', []);
    if isempty(signedTurnWindows)
        % Backward-compatible behavior used in earlier versions.
        trialDeviances = linspace(dirVar, 360 - dirVar, nTrials);
        return;
    end

    turnWindows = local_parse_signed_turn_windows(signedTurnWindows);
    windowWidths = turnWindows(:, 2) - turnWindows(:, 1);
    totalWidth = sum(windowWidths);
    cumulativeWidth = [0; cumsum(windowWidths)];

    % Deterministic, evenly spaced coverage over the union of allowed windows.
    % Trial order is still randomized later by TrialOrder in CreateInputFiles.
    samplePositions = ((1:nTrials)' - 0.5) / nTrials * totalWidth;
    for trialIndex = 1:nTrials
        windowIndex = find(samplePositions(trialIndex) <= cumulativeWidth(2:end), 1, 'first');
        offsetWithinWindow = samplePositions(trialIndex) - cumulativeWidth(windowIndex);
        trialDeviances(trialIndex) = turnWindows(windowIndex, 1) + offsetWithinWindow;
    end
end

function turnWindows = local_parse_signed_turn_windows(rawTurnWindows)
% LOCAL_PARSE_SIGNED_TURN_WINDOWS Validate and normalize signed turn windows.
%
% Purpose:
%   Normalize user-configured signed windows to an [N x 2] numeric matrix
%   and validate domain constraints for likelihood deviant turns.
%
% Inputs:
%   - rawTurnWindows: either [N x 2] matrix or 1D vector with pairwise
%     boundaries [min1, max1, min2, max2, ...].
%
% Output:
%   - turnWindows: sorted [N x 2] intervals in signed degrees.
%
% Assumptions:
%   - Signed angle domain is [-180, 180].
%   - Each row is [min, max] with min < max.
%   - Windows do not overlap.
    if ~isnumeric(rawTurnWindows) || isempty(rawTurnWindows)
        error(['deviantSignedTurnWindows must be numeric and non-empty when provided. ' ...
            'Use [] to disable explicit signed windows.']);
    end

    if isvector(rawTurnWindows)
        if mod(numel(rawTurnWindows), 2) ~= 0
            error('Vector deviantSignedTurnWindows must have an even number of elements.');
        end
        turnWindows = reshape(rawTurnWindows, 2, [])';
    elseif size(rawTurnWindows, 2) == 2
        turnWindows = rawTurnWindows;
    else
        error('deviantSignedTurnWindows must be Nx2 or a vector of paired bounds.');
    end

    if any(~isfinite(turnWindows(:)))
        error('deviantSignedTurnWindows contains non-finite values.');
    end
    if any(turnWindows(:, 1) >= turnWindows(:, 2))
        error('Each deviantSignedTurnWindows row must satisfy min < max.');
    end
    if any(turnWindows(:) < -180 | turnWindows(:) > 180)
        error('deviantSignedTurnWindows bounds must stay within [-180, 180].');
    end

    turnWindows = sortrows(turnWindows, 1);
    if size(turnWindows, 1) > 1
        if any(turnWindows(2:end, 1) < turnWindows(1:end-1, 2))
            error('deviantSignedTurnWindows rows must not overlap.');
        end
    end
end

function value = local_get_struct_field_or_default(structValue, fieldName, defaultValue)
% LOCAL_GET_STRUCT_FIELD_OR_DEFAULT Read an optional struct field safely.
    value = defaultValue;
    if isstruct(structValue) && isfield(structValue, fieldName)
        value = structValue.(fieldName);
    end
end

%% Local helpers (dRSA-proxy scoring)
function meanAbsCorrScore = local_compute_drsa_proxy_score( ...
        dot1Positions, dot1Directions, frameIndices, excludeDiagLagBins)
% LOCAL_COMPUTE_DRSA_PROXY_SCORE Compute a dRSA-style proxy score.
%
% Purpose:
%   Compute the acceptance-gate proxy used in v14 by mirroring the
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
