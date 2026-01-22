%% Stimuli Generation (uniform starts, boundary-safe placement)
% Script: stimuli_generation_v7.m
% Author: Marisa (original), Ayman Hatoum (v5), updated by Dami (v6, v7)
%
% Purpose:
%   Generate dot-motion stimuli with uniform starting positions while keeping
%   the rectangular bounds strictly off-limits. This version removes
%   segmented sub-paths and center-biased start points, and instead:
%     - draws initial positions uniformly within the feasible rectangle,
%     - advances positions with constant step size, and
%     - places paths so they fit the boundary without per-frame rejection.
%
%   v7 update:
%     - curvature (curvynessFactor) is randomized independently for each dot
%       on every trial, while staying constant within the trial,
%     - curvature magnitude is bounded by geometric limits derived from the
%       trial duration and the dot rectangle,
%     - start positions are chosen after path generation so the full path
%       fits within bounds (faster than rejection sampling),
%     - plots the generated dot paths with time-graded colors at the end.
%
%   Output matches the structure expected by:
%     - experiment/MoveDot1_experiment_vX.m
%     - simulations/build_movdot_simulation_inputs.m
%
% Example usage (from repo root in MATLAB):
%   addpath('experiment');
%   stimuli_generation_v7;
%   % Follow the dialog prompts for viewing distance and subject ID.
%
% Example usage (custom working directory):
%   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD');
%   addpath('experiment');
%   stimuli_generation_v7;
%
% Inputs:
%   - Config (class/struct on MATLAB path) with screen, dot, and timing params.
%   - Utils helper class on MATLAB path.
%   - Dialog values for viewing distance and subject ID.
%
% Outputs:
%   - Saves xySeqs and Cfg to Config.inputDirectory using Config.stimuliFileName.
%   - xySeqs(...).xy is frames x 4 with columns [x1 y1 x2 y2] in visual degrees.
%   - Saves path plots to simulations/debug as PNGs.
%
% Key assumptions:
%   - Coordinates are in visual degrees in the rectangle [0..rectSize].
%   - Deviant logic and step size are preserved from v6.
%   - Curvyness uses randomized valence and bounded magnitude per dot per trial.
%   - Start positions are solved from the generated relative paths to avoid
%     per-frame boundary rejection.
%   - Path placement removes boundary-related rejection bias.

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
% Data flow: condition params -> trial loops -> bounded curvature + path fit -> xySeqs output struct.
framesPerTrial = round(Config.trialDuration * Config.frameFrequency);
minPos = [Config.dotWidth/2, Config.dotWidth/2];
maxPos = Config.dotRectSize - minPos;
% Curvature bounds: min fits largest circle in rect, max yields full circle per trial.
dotRectInner = Config.dotRectSize - Config.dotWidth;
maxRadius = min(dotRectInner) / 2;
curvMinDeg = rad2deg(Config.dotSpeedDegPerFrame / maxRadius);
curvMaxDeg = 360 / framesPerTrial;
if curvMinDeg > curvMaxDeg
    warning('Curvature bounds invalid (curvMin > curvMax). Using curvMax for both.');
    curvMinDeg = curvMaxDeg;
end

for dirVarIndex = 1:length(stimulusTypeConfig.directionVariance)
    dirVar = stimulusTypeConfig.directionVariance(dirVarIndex);

    for pathDurIndex = 1:length(stimulusTypeConfig.pathDuration)
        pathDuration = stimulusTypeConfig.pathDuration(pathDurIndex);
        trialDeviances = linspace(dirVar, 360 - dirVar, Config.trialsPerCondition);

        totalSpread = zeros(Config.numXGrids, Config.numYGrids, 2);

        for trialPerCondIndex = 1:Config.trialsPerCondition
            % Fast generation: build paths first, then place them inside bounds.
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

                % Per-trial per-dot curvyness (constant within trial, bounded per dot)
                % Data flow: random sign + bounded magnitude -> degrees per frame.
                curvynessFactor = [ ...
                    Utils.RandPosNeg() * (curvMinDeg + rand(1) * (curvMaxDeg - curvMinDeg)), ...
                    Utils.RandPosNeg() * (curvMinDeg + rand(1) * (curvMaxDeg - curvMinDeg))];

                if Config.stimulusType == Utils.likelihood
                    devOnsetVariance = rand(1) * Utils.RandPosNeg() * Config.deviantOnsetVariance * Config.deviantOnset;
                    devFrameOnset = round((Config.deviantOnset + devOnsetVariance) * framesPerTrial);
                    devFrameOnset = min(max(devFrameOnset, 2), framesPerTrial);
                end

                % Build relative paths from origin (dot start at [0,0]).
                relPaths = zeros(framesPerTrial, 4);
                relPathsNoDev = zeros(framesPerTrial, 4);

                for frameIndex = 1:framesPerTrial
                    if frameIndex == 1
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

                        % Update real vs baseline angles, then compute relative positions.
                        directionAngleChangeNoDev = directionAngleChange;
                        if Config.stimulusType == Utils.likelihood && dirVar ~= 0
                            directionAngleChangeNoDev = 0; % baseline represents the no-deviant path
                        end
                        directionAngle = directionAngle + curvynessFactor + directionAngleChange;
                        directionAngleNoDev = directionAngleNoDev + curvynessFactor + directionAngleChangeNoDev;

                        directionVector = Utils.GetDirectionVector(directionAngle, length(directionAngle));
                        dummyDirectionVector = Utils.GetDirectionVector(directionAngleNoDev, length(directionAngleNoDev));
                        relPaths(frameIndex, :) = relPaths(frameIndex-1, :) + directionVector * Config.dotSpeedDegPerFrame;
                        relPathsNoDev(frameIndex, :) = relPathsNoDev(frameIndex-1, :) + dummyDirectionVector * Config.dotSpeedDegPerFrame;
                        allPathsDirectionAngle(frameIndex, :) = directionAngle;
                        allPathsCurvyness(frameIndex, :) = curvynessFactor;

                        if any(directionAngleChange)
                            allPathsStartingPoint(frameIndex, :) = [1, 1];
                        else
                            allPathsStartingPoint(frameIndex, :) = [0, 0];
                        end
                    end
                end

                % Choose start positions so both real and baseline paths fit inside bounds.
                dot1X = [relPaths(:, 1); relPathsNoDev(:, 1)];
                dot1Y = [relPaths(:, 2); relPathsNoDev(:, 2)];
                dot2X = [relPaths(:, 3); relPathsNoDev(:, 3)];
                dot2Y = [relPaths(:, 4); relPathsNoDev(:, 4)];

                dot1XMin = min(dot1X);
                dot1XMax = max(dot1X);
                dot1YMin = min(dot1Y);
                dot1YMax = max(dot1Y);
                dot2XMin = min(dot2X);
                dot2XMax = max(dot2X);
                dot2YMin = min(dot2Y);
                dot2YMax = max(dot2Y);

                dot1XRange = [minPos(1) - dot1XMin, maxPos(1) - dot1XMax];
                dot1YRange = [minPos(2) - dot1YMin, maxPos(2) - dot1YMax];
                dot2XRange = [minPos(1) - dot2XMin, maxPos(1) - dot2XMax];
                dot2YRange = [minPos(2) - dot2YMin, maxPos(2) - dot2YMax];

                if dot1XRange(1) > dot1XRange(2) || dot1YRange(1) > dot1YRange(2) || ...
                        dot2XRange(1) > dot2XRange(2) || dot2YRange(1) > dot2YRange(2)
                    trialValid = false;
                    continue;
                end

                startDot1 = [ ...
                    dot1XRange(1) + rand(1) * (dot1XRange(2) - dot1XRange(1)), ...
                    dot1YRange(1) + rand(1) * (dot1YRange(2) - dot1YRange(1))];
                startDot2 = [ ...
                    dot2XRange(1) + rand(1) * (dot2XRange(2) - dot2XRange(1)), ...
                    dot2YRange(1) + rand(1) * (dot2YRange(2) - dot2YRange(1))];

                frameDotXY = relPaths + [startDot1 startDot2];
                dummyFrameDotXY = relPathsNoDev + [startDot1 startDot2];

                % Grid occupancy bookkeeping (no longer used for acceptance).
                for frameIndex = 1:framesPerTrial
                    allPathsFrameDotXY(frameIndex, :) = frameDotXY(frameIndex, :);
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

%% Plot generated paths (time-graded scatter, separated by condition and dot)
% Data flow: xySeqs (dirVar x pathDur x trials) -> dot paths arrays -> plots.
[nDirVar, nPathDur, ~] = size(xySeqs);
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
debugDir = fullfile(repoRoot, 'simulations', 'debug');
if ~exist(debugDir, 'dir')
    mkdir(debugDir);
end

for dirVarIndex = 1:nDirVar
    for pathDurIndex = 1:nPathDur
        seqs = squeeze(xySeqs(dirVarIndex, pathDurIndex, :));
        hasXY = arrayfun(@(seq) isfield(seq, 'xy') && ~isempty(seq.xy), seqs);
        seqs = seqs(hasXY);

        if isempty(seqs)
            warning('No xy sequences for condition dirVar=%d, pathDurIndex=%d.', dirVarIndex, pathDurIndex);
            continue;
        end

        nTime = size(seqs(1).xy, 1);
        valid = arrayfun(@(seq) size(seq.xy, 1) == nTime, seqs);
        seqs = seqs(valid);

        if isempty(seqs)
            warning('No consistent-length xy sequences for condition dirVar=%d, pathDurIndex=%d.', ...
                dirVarIndex, pathDurIndex);
            continue;
        end

        nSeqs = numel(seqs);
        dot1Paths = zeros(nSeqs, 2, nTime);
        dot2Paths = zeros(nSeqs, 2, nTime);

        for seqIndex = 1:nSeqs
            seqXY = seqs(seqIndex).xy;
            dot1Paths(seqIndex, 1, :) = seqXY(:, 1);
            dot1Paths(seqIndex, 2, :) = seqXY(:, 2);
            dot2Paths(seqIndex, 1, :) = seqXY(:, 3);
            dot2Paths(seqIndex, 2, :) = seqXY(:, 4);
        end

        % Shared plotting parameters (match plot_paths.m style).
        cmap = parula(nTime);
        timeIdx = repmat(1:nTime, nSeqs, 1);
        timeIdx = timeIdx(:);

        conditionLabel = sprintf('dirVar=%g pathDur=%g', ...
            seqs(1).condition, seqs(1).PredictionRange);
        safeLabel = regexprep(conditionLabel, '[^A-Za-z0-9_-]+', '_');

        % Dot 1 plot
        xVals = reshape(dot1Paths(:, 1, :), [], 1);
        yVals = reshape(dot1Paths(:, 2, :), [], 1);
        pointColors = cmap(timeIdx, :);

        figDot1 = figure('Name', ['Dot 1 paths - ' conditionLabel], 'NumberTitle', 'off');
        axDot1 = axes(figDot1);
        hold(axDot1, 'on');
        scatter(axDot1, xVals, yVals, 8, pointColors, 'filled');
        colormap(axDot1, cmap);
        caxis(axDot1, [1 nTime]);
        cbDot1 = colorbar(axDot1);
        cbDot1.Label.String = 'Time (samples)';
        xlabel(axDot1, 'X (visual degrees)');
        ylabel(axDot1, 'Y (visual degrees)');
        title(axDot1, ['Dot 1 paths - ' conditionLabel]);
        axis(axDot1, 'equal');
        hold(axDot1, 'off');

        % Dot 2 plot
        xVals = reshape(dot2Paths(:, 1, :), [], 1);
        yVals = reshape(dot2Paths(:, 2, :), [], 1);
        pointColors = cmap(timeIdx, :);

        figDot2 = figure('Name', ['Dot 2 paths - ' conditionLabel], 'NumberTitle', 'off');
        axDot2 = axes(figDot2);
        hold(axDot2, 'on');
        scatter(axDot2, xVals, yVals, 8, pointColors, 'filled');
        colormap(axDot2, cmap);
        caxis(axDot2, [1 nTime]);
        cbDot2 = colorbar(axDot2);
        cbDot2.Label.String = 'Time (samples)';
        xlabel(axDot2, 'X (visual degrees)');
        ylabel(axDot2, 'Y (visual degrees)');
        title(axDot2, ['Dot 2 paths - ' conditionLabel]);
        axis(axDot2, 'equal');
        hold(axDot2, 'off');

        % Save figures alongside simulations/debug outputs.
        print(figDot1, fullfile(debugDir, ...
            sprintf('stimuli_generation_v7_dot1_%s.png', safeLabel)), '-dpng', '-r300');
        print(figDot2, fullfile(debugDir, ...
            sprintf('stimuli_generation_v7_dot2_%s.png', safeLabel)), '-dpng', '-r300');
    end
end
