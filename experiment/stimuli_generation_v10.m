%% Stimuli Generation (uniform starts, boundary rejection, curvature flip on deviant)
% Script: stimuli_generation_v10.m
% Author: Marisa (original), Ayman Hatoum (v5), updated by Dami (v6-v8), v9-v10 updates by Codex
%
% Purpose:
%   Generate dot-motion stimuli with uniform starting positions while keeping
%   the rectangular bounds strictly off-limits. This version removes
%   segmented sub-paths and center-biased start points, and instead:
%     - draws initial positions uniformly within the allowed rectangle,
%     - advances positions with constant step size, and
%     - rejects any path that would cross the boundary.
%   When a deviant occurs (likelihood stimulus), the curvature sign can
%   flip from that frame onward for the deviant path (configurable below).
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
%   stimuli_generation_v10;
%   % Follow the dialog prompts for viewing distance and subject ID.
%
% Example usage (custom working directory):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   addpath('experiment');
%   stimuli_generation_v10;
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
%   - Saves xySeqsWannabeDev and Cfg to input_files/MovDot_SubXX_wannabeDev.mat.
%     The "wannabe" file stores no-deviant baseline paths for deviant conditions
%     only (analysis use; not presented to participants).
%   - Optional per-condition path plots after generation using time colors.
%
% Key assumptions:
%   - Coordinates are in visual degrees in the rectangle [0..rectSize].
%   - Deviant logic, curviness, and step size are preserved from v5.
%   - Dots must remain at least Config.minDistanceBetweenDots apart.
%   - If enabled, curvature flips only for the deviant path (baseline keeps
%     constant curvature for boundary validation).
%   - No-deviant baselines reuse the same start points and curvature, but
%     remove deviant angle changes.
%   - Rejection sampling introduces boundary-related selection effects.

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
flipCurvatureOnDeviant = false; % set true to flip curvature sign at deviant onset

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
wannabeFile = [Config.inputDirectory sprintf('MovDot_Sub%02d_wannabeDev.mat', subjectID)];
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
% Data flow: condition params -> trial loops -> xySeqs + xySeqsWannabeDev outputs.
if ~skipGeneration
    framesPerTrial = round(Config.trialDuration * Config.frameFrequency);
    minPos = [Config.dotWidth/2, Config.dotWidth/2];
    maxPos = Config.dotRectSize - minPos;
    xySeqsWannabeDev = struct();

    for dirVarIndex = 1:length(stimulusTypeConfig.directionVariance)
        dirVar = stimulusTypeConfig.directionVariance(dirVarIndex);
        isDeviantCondition = (Config.stimulusType == Utils.likelihood && dirVar ~= 0);

        for pathDurIndex = 1:length(stimulusTypeConfig.pathDuration)
            pathDuration = stimulusTypeConfig.pathDuration(pathDurIndex);
            trialDeviances = linspace(dirVar, 360 - dirVar, Config.trialsPerCondition);

            totalSpread = zeros(Config.numXGrids, Config.numYGrids, 2);

            for trialPerCondIndex = 1:Config.trialsPerCondition
                % Rejection sampling: keep regenerating until a full trial stays in-bounds.
                trialValid = false;
                while ~trialValid
                    trialValid = true;

                % Per-trial parameters: directions, curvyness, and start points.
                directionAngle = [Utils.RandAngleDegree(), Utils.RandAngleDegree()];
                directionAngleNoDev = directionAngle;
                curvynessFactor = Config.curvFactor * [ ...
                    Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand), ...
                    Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand)];

                % Uniform starting positions inside bounds.
                % Data flow: random start -> distance check -> accept/reject trial.
                startDot1 = rand(1, 2) .* (maxPos - minPos) + minPos;
                startDot2 = rand(1, 2) .* (maxPos - minPos) + minPos;
                frameDotXY = [startDot1 startDot2];
                if norm(startDot1 - startDot2) < Config.minDistanceBetweenDots
                    trialValid = false;
                    continue;
                end

                % Vectorized path generation: build per-frame angles and positions.
                % Data flow: deviant onset -> optional curvature flip -> cumulative angles.
                numSteps = framesPerTrial - 1;
                devFrameOnset = [];
                directionAngleChange = zeros(max(numSteps, 0), 1);
                if numSteps > 0
                    if Config.stimulusType == Utils.likelihood
                        if dirVar ~= 0
                            devOnsetVariance = rand(1) * Utils.RandPosNeg() * Config.deviantOnsetVariance * Config.deviantOnset;
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

                % Data flow: per-step angle changes + curvature -> cumulative angles.
                directionAngleChange = repmat(directionAngleChange, 1, 2);
                directionAngleChangeNoDev = directionAngleChange;
                if Config.stimulusType == Utils.likelihood && dirVar ~= 0
                    directionAngleChangeNoDev(:) = 0; % baseline represents the no-deviant path
                end

                curvynessPerStep = repmat(curvynessFactor, max(numSteps, 0), 1);
                curvynessPerStepNoDev = curvynessPerStep;
                if flipCurvatureOnDeviant && ~isempty(devFrameOnset) && numSteps > 0
                    flipIndex = devFrameOnset - 1;
                    numFlipSteps = numSteps - flipIndex + 1;
                    curvynessPerStep(flipIndex:end, :) = repmat(-curvynessFactor, numFlipSteps, 1);
                end

                allPathsDirectionAngle = zeros(framesPerTrial, 2);
                allPathsDirectionAngleNoDev = zeros(framesPerTrial, 2);
                allPathsDirectionAngle(1, :) = directionAngle;
                allPathsDirectionAngleNoDev(1, :) = directionAngleNoDev;
                if numSteps > 0
                    allPathsDirectionAngle(2:end, :) = directionAngle + cumsum(directionAngleChange + curvynessPerStep, 1);
                    allPathsDirectionAngleNoDev(2:end, :) = directionAngleNoDev + cumsum(directionAngleChangeNoDev + curvynessPerStepNoDev, 1);
                end

                % Data flow: angles -> direction vectors -> cumulative positions.
                directionVectors = [cosd(allPathsDirectionAngle(:, 1)), sind(allPathsDirectionAngle(:, 1)), ...
                    cosd(allPathsDirectionAngle(:, 2)), sind(allPathsDirectionAngle(:, 2))];
                dummyDirectionVectors = [cosd(allPathsDirectionAngleNoDev(:, 1)), sind(allPathsDirectionAngleNoDev(:, 1)), ...
                    cosd(allPathsDirectionAngleNoDev(:, 2)), sind(allPathsDirectionAngleNoDev(:, 2))];

                allPathsFrameDotXY = zeros(framesPerTrial, 4);
                dummyPathsFrameDotXY = zeros(framesPerTrial, 4);
                allPathsFrameDotXY(1, :) = frameDotXY;
                dummyPathsFrameDotXY(1, :) = frameDotXY;
                if numSteps > 0
                    stepVectors = directionVectors(2:end, :) * Config.dotSpeedDegPerFrame;
                    dummyStepVectors = dummyDirectionVectors(2:end, :) * Config.dotSpeedDegPerFrame;
                    allPathsFrameDotXY(2:end, :) = frameDotXY + cumsum(stepVectors, 1);
                    dummyPathsFrameDotXY(2:end, :) = frameDotXY + cumsum(dummyStepVectors, 1);
                end

                % Constraint checks: distance and boundary tests across frames.
                if numSteps > 0
                    trialViolation = zeros(numSteps, 1);
                    interDotDelta = allPathsFrameDotXY(2:end, 1:2) - allPathsFrameDotXY(2:end, 3:4);
                    interDotDist = sqrt(sum(interDotDelta.^2, 2));
                    distTooClose = interDotDist < Config.minDistanceBetweenDots;
                    leftHit = any(allPathsFrameDotXY(2:end, [1 3]) < minPos(1), 2) ...
                        | any(dummyPathsFrameDotXY(2:end, [1 3]) < minPos(1), 2);
                    topHit = any(allPathsFrameDotXY(2:end, [2 4]) < minPos(2), 2) ...
                        | any(dummyPathsFrameDotXY(2:end, [2 4]) < minPos(2), 2);
                    rightHit = any(allPathsFrameDotXY(2:end, [1 3]) > maxPos(1), 2) ...
                        | any(dummyPathsFrameDotXY(2:end, [1 3]) > maxPos(1), 2);
                    bottomHit = any(allPathsFrameDotXY(2:end, [2 4]) > maxPos(2), 2) ...
                        | any(dummyPathsFrameDotXY(2:end, [2 4]) > maxPos(2), 2);

                    trialViolation(distTooClose) = 1;
                    trialViolation(trialViolation == 0 & leftHit) = 2;
                    trialViolation(trialViolation == 0 & topHit) = 3;
                    trialViolation(trialViolation == 0 & rightHit) = 4;
                    trialViolation(trialViolation == 0 & bottomHit) = 5;

                    firstViolation = find(trialViolation > 0, 1, "first");
                    if ~isempty(firstViolation)
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
                % Data flow: dummy path -> xySeqsWannabeDev output struct.
                if isDeviantCondition
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).condition       = dirVar;
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).PredictionRange = pathDuration;
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).sequence        = trialPerCondIndex;
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).xy              = dummyPathsFrameDotXY;
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).pathAll         = allPathsStartingPointNoDev;
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).curvyness       = allPathsCurvynessNoDev;
                    xySeqsWannabeDev(dirVarIndex, pathDurIndex, trialPerCondIndex).AngleDirection  = allPathsDirectionAngleNoDev;
                end

                totalSpread = totalSpread + gridSum;
                fprintf("End of trial %d within condition %d\n", trialPerCondIndex, dirVar);
            end

            xySeqs(dirVarIndex, pathDurIndex).totalSpread = totalSpread;
            fprintf("End of condition [%d] = %d\n", dirVarIndex, dirVar);
        end
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
        'version', 'v10', ...
        'parameters', struct( ...
            'renderPreview', renderPreview, ...
            'plotPathsAtEnd', plotPathsAtEnd, ...
            'plotPathsPerCondition', plotPathsPerCondition, ...
            'plotPathConditions', plotPathConditions, ...
            'flipCurvatureOnDeviant', flipCurvatureOnDeviant));
    repro.inputs = struct('subjectID', subjectID, 'viewingDistance', viewingDistance);
    repro.rng = rngState;
    repro.config = configSnapshot;
    repro.stimulusTypeConfig = stimulusTypeConfig;
    repro.derived = struct('framesPerTrial', framesPerTrial, 'minPos', minPos, 'maxPos', maxPos);

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
    xySeqsWannabeDev = repmat(xySeqsWannabeDev, 1, Config.trialRepetetion, 1);
    save(wannabeFile, 'xySeqsWannabeDev', 'Cfg');
end

% Close all onscreens and offscreens
if renderPreview
    sca;
end
