%% Stimuli Generation (uniform starts, boundary rejection)
% Script: stimuli_generation_v6.m
% Author: Marisa (original), Ayman Hatoum (v5), updated by Dami (v6)
%
% Purpose:
%   Generate dot-motion stimuli with uniform starting positions while keeping
%   the rectangular bounds strictly off-limits. This version removes
%   segmented sub-paths and center-biased start points, and instead:
%     - draws initial positions uniformly within the allowed rectangle,
%     - advances positions with constant step size, and
%     - rejects any path that would cross the boundary.
%
%   Output matches the structure expected by:
%     - experiment/MoveDot1_experiment_vX.m
%     - simulations/build_movdot_simulation_inputs.m
%
% Example usage (from repo root in MATLAB):
%   addpath('experiment');
%   stimuli_generation_v6;
%   % Follow the dialog prompts for viewing distance and subject ID.
%
% Example usage (custom working directory):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   addpath('experiment');
%   stimuli_generation_v6;
%
% Inputs:
%   - Config (class/struct on MATLAB path) with screen, dot, and timing params.
%   - Utils helper class on MATLAB path.
%   - Dialog values for viewing distance and subject ID.
%
% Outputs:
%   - Saves xySeqs and Cfg to Config.inputDirectory using Config.stimuliFileName.
%   - xySeqs(...).xy is frames x 4 with columns [x1 y1 x2 y2] in visual degrees.
%
% Key assumptions:
%   - Coordinates are in visual degrees in the rectangle [0..rectSize].
%   - Deviant logic, curviness, and step size are preserved from v5.
%   - Rejection sampling introduces boundary-related selection effects.

addpath('lib/');

%% Setup environment
% Data flow: user dialog -> viewing distance/subject ID -> RNG seed.
clear all;
clc;
ShowCursor;
renderPreview = false; % set true to draw stimuli during generation

userDoubles = Utils.GetUserDoubles(Config.dialogTitle, Config.dialogDimensions, ...
    Config.dialogPrompts, Config.dialogDefaults);
viewingDistance = userDoubles(1);
subjectID = userDoubles(2);

% Seed RNG for reproducibility (per subject)
rng(subjectID);

% Ensure output folders exist before saving
if ~exist(Config.inputDirectory, 'dir')
    mkdir(Config.inputDirectory);
end
if ~exist(Config.outputDirectory, 'dir')
    mkdir(Config.outputDirectory);
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
% Data flow: condition params -> trial loops -> xySeqs output struct.
framesPerTrial = round(Config.trialDuration * Config.frameFrequency);
minPos = [Config.dotWidth/2, Config.dotWidth/2];
maxPos = Config.dotRectSize - minPos;

for dirVarIndex = 1:length(stimulusTypeConfig.directionVariance)
    dirVar = stimulusTypeConfig.directionVariance(dirVarIndex);

    for pathDurIndex = 1:length(stimulusTypeConfig.pathDuration)
        pathDuration = stimulusTypeConfig.pathDuration(pathDurIndex);
        trialDeviances = linspace(dirVar, 360 - dirVar, Config.trialsPerCondition);

        totalSpread = zeros(Config.numXGrids, Config.numYGrids, 2);

        for trialPerCondIndex = 1:Config.trialsPerCondition
            % Rejection sampling: keep regenerating until a full trial stays in-bounds.
            trialValid = false;
            while ~trialValid
                trialValid = true;

                % Per-trial buffers
                allPathsFrameDotXY = zeros(framesPerTrial, 4);
                allPathsStartingPoint = zeros(framesPerTrial, 2);
                allPathsDirectionAngle = zeros(framesPerTrial, 2);
                allPathsCurvyness = zeros(framesPerTrial, 2);
                gridCount = zeros(Config.numXGrids, Config.numYGrids, framesPerTrial, 2);

                % Initial directions (one per dot)
                directionAngle = [Utils.RandAngleDegree(), Utils.RandAngleDegree()];
                directionAngleNoDev = directionAngle;

                % Per-trial curvyness (constant over the trial)
                curvynessFactor = Config.curvFactor * [ ...
                    Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand), ...
                    Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand)];

                % Uniform starting positions inside bounds
                startDot1 = rand(1, 2) .* (maxPos - minPos) + minPos;
                startDot2 = rand(1, 2) .* (maxPos - minPos) + minPos;
                frameDotXY = [startDot1 startDot2];
                dummyFrameDotXY = frameDotXY; % baseline path for deviant logic

                if Config.stimulusType == Utils.likelihood
                    devOnsetVariance = rand(1) * Utils.RandPosNeg() * Config.deviantOnsetVariance * Config.deviantOnset;
                    devFrameOnset = round((Config.deviantOnset + devOnsetVariance) * framesPerTrial);
                    devFrameOnset = min(max(devFrameOnset, 2), framesPerTrial);
                end

                for frameIndex = 1:framesPerTrial
                    if frameIndex == 1
                        allPathsFrameDotXY(frameIndex, :) = frameDotXY;
                        allPathsStartingPoint(frameIndex, :) = [1, 1];
                        allPathsDirectionAngle(frameIndex, :) = directionAngle;
                        allPathsCurvyness(frameIndex, :) = curvynessFactor;
                    else
                        % Compute per-frame direction changes based on stimulus type.
                        if Config.stimulusType == Utils.likelihood
                            if frameIndex == devFrameOnset && dirVar ~= 0
                                directionAngleChange = trialDeviances(trialPerCondIndex);
                            else
                                directionAngleChange = 0;
                            end
                        else
                            if frameIndex == 2
                                if Config.stimulusType == Utils.path_duration
                                    directionAngleChange = Utils.RandAngleDegree();
                                elseif Config.stimulusType == Utils.path_duration_norm
                                    directionAngleChange = normrnd(0, stimulusTypeConfig.directionChange);
                                end
                            else
                                directionAngleChange = 0;
                            end
                        end

                        % Update real vs baseline angles, then compute candidate positions.
                        directionAngleChangeNoDev = directionAngleChange;
                        if Config.stimulusType == Utils.likelihood && dirVar ~= 0
                            directionAngleChangeNoDev = 0; % baseline represents the no-deviant path
                        end
                        directionAngle = directionAngle + curvynessFactor + directionAngleChange;
                        directionAngleNoDev = directionAngleNoDev + curvynessFactor + directionAngleChangeNoDev;

                        directionVector = Utils.GetDirectionVector(directionAngle, length(directionAngle));
                        nextDotXY = allPathsFrameDotXY(frameIndex-1, :) + directionVector * Config.dotSpeedDegPerFrame;

                        dummyDirectionVector = Utils.GetDirectionVector(directionAngleNoDev, length(directionAngleNoDev));
                        dummyDotXY = dummyFrameDotXY + dummyDirectionVector * Config.dotSpeedDegPerFrame;

                        % Reject the trial if any dot crosses boundaries.
                        if any(nextDotXY([1 3]) < minPos(1)) || any(nextDotXY([1 3]) > maxPos(1)) || ...
                                any(nextDotXY([2 4]) < minPos(2)) || any(nextDotXY([2 4]) > maxPos(2)) || ...
                                any(dummyDotXY([1 3]) < minPos(1)) || any(dummyDotXY([1 3]) > maxPos(1)) || ...
                                any(dummyDotXY([2 4]) < minPos(2)) || any(dummyDotXY([2 4]) > maxPos(2))
                            trialValid = false;
                            break;
                        end

                        allPathsFrameDotXY(frameIndex, :) = nextDotXY;
                        dummyFrameDotXY = dummyDotXY;
                        allPathsDirectionAngle(frameIndex, :) = directionAngle;
                        allPathsCurvyness(frameIndex, :) = curvynessFactor;

                        if any(directionAngleChange)
                            allPathsStartingPoint(frameIndex, :) = [1, 1];
                        else
                            allPathsStartingPoint(frameIndex, :) = [0, 0];
                        end
                    end

                % Grid occupancy bookkeeping (no longer used for acceptance).
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

            totalSpread = totalSpread + gridSum;
            fprintf("End of trial %d within condition %d\n", trialPerCondIndex, dirVar);
        end

        xySeqs(dirVarIndex, pathDurIndex).totalSpread = totalSpread;
        fprintf("End of condition [%d] = %d\n", dirVarIndex, dirVar);
    end
end

%% Save outputs (xySeqs + Cfg)
% Data flow: Config -> Cfg struct -> .mat on disk.
Cfg = struct();
Cfg.dpf = Config.dotSpeedDegPerFrame;
Cfg.fps = Config.frameFrequency;
Cfg.Stimulitype = Config.stimulusType;
Cfg.dot_w = Config.dotWidth;
Cfg.rectSize = Config.dotRectSize;
Cfg.DirChange = stimulusTypeConfig.directionChange;

xySeqs = repmat(xySeqs, 1, Config.trialRepetetion, 1);
save([Config.inputDirectory sprintf(Config.stimuliFileName, subjectID)], 'xySeqs', 'Cfg');

% Close all onscreens and offscreens
if renderPreview
    sca;
end
