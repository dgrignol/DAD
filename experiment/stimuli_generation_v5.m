%% Stimuli Generation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Ayman Hatoum
%

addpath('lib/');

%setup environment
clear all;
clc;
ShowCursor;
Screen('Preference', 'SkipSyncTests', 1);

%shows transparent if needed
PsychDebugWindowConfiguration([], .02);

%block the command prompt inout while the app is running
%ListenChar(2);

userDoubles = Utils.GetUserDoubles(Config.dialogTitle, Config.dialogDimensions, ...
    Config.dialogPrompts, Config.dialogDefaults);

viewingDistance = userDoubles(1);
subjectID = userDoubles(2);

%get monitor pixels and millis
monitorIndex = max(Screen('Screens'));
monitorPixelSize = Screen('Rect', monitorIndex);
[monitorMilliWidth, monitorMilliHeight] = Screen('DisplaySize', monitorIndex);

%divide the pixels width by the visual angle degrees
pixelsPerDegrees = monitorPixelSize(3) - monitorPixelSize(1); %monitor pixels width
pixelsPerDegrees = pixelsPerDegrees / Utils.ComputeVisualAngleDegrees(viewingDistance, monitorMilliWidth);

%their code computes millis per degree but not used

%setup win
win = Screen('OpenWindow', monitorIndex, Config.screenBackgroundColor, Config.screenRect);

screenCenter = floor((Config.screenRect(3:4) - Config.screenRect(1:2))/2); %in pixels: take deltaX and deltaY divide by two
screenCenterDegrees = screenCenter / pixelsPerDegrees;

%dot rect coordinates in degrees
dotRectDegrees = [screenCenterDegrees(1)-Config.dotRectSize(1)/2 screenCenterDegrees(2)-Config.dotRectSize(2)/2 ...
    screenCenterDegrees(1)+Config.dotRectSize(1)/2 screenCenterDegrees(2)+Config.dotRectSize(2)/2];

dotRectCenterDegrees = floor((dotRectDegrees(3:4) - dotRectDegrees(1:2))/2);

%randomly assign +/- to the curvyness factor, for both dots
curvynessFactor = [Config.curvFactor*Utils.RandPosNeg(), Config.curvFactor*Utils.RandPosNeg()];

%compute size of dot in pixels and the best fit without aliasing
dotSizePixels = Config.dotWidth*pixelsPerDegrees;
[minSmooth, maxSmooth] = Screen('DrawDots', win);
dotSizePixels = min(max(dotSizePixels, minSmooth), maxSmooth);

%set the stimulus type config struct
stimulusTypeConfig = struct();
if Config.stimulusType == Utils.likelihood
    stimulusTypeConfig = Config.likelihood;
elseif Config.stimulusType == Utils.path_duration
    stimulusTypeConfig = Config.path_duration;
elseif Config.stimulusType == Utils.path_duration_norm
    stimulusTypeConfig = Config.path_duration_norm;
end

%recreating their loops

%new comment: seems to not make sense that it accounts for jitter since
%already that is counted when you compute the number of frames
extraNumPaths = 0; %this one they add it for counting for the jitter, to be checked and then moved to the Config

%dynamic variables they defined - to be checked if can be improved
spreadCount(length(stimulusTypeConfig.directionVariance), ...
    length(stimulusTypeConfig.pathDuration), ...
    Config.trialsPerCondition).gridCount(:,:,:) = zeros(Config.numXGrids, Config.numYGrids, 2);

for dirVarIndex = 1:length(stimulusTypeConfig.directionVariance)    %so this one define the 1st factor

    dirVar = stimulusTypeConfig.directionVariance(dirVarIndex);

    for pathDurIndex = 1:length(stimulusTypeConfig.pathDuration)    %this one define a second factor but just for the design
        
        pathDuration = stimulusTypeConfig.pathDuration(pathDurIndex);
        
        %get the possible amount of paths within the trial
        numPaths = (Config.trialDuration/pathDuration) + extraNumPaths; 

        trialDeviances = linspace(dirVar, 360 - dirVar, Config.trialsPerCondition);

        outOfLimitSpread = 1;
        while outOfLimitSpread

        for trialPerCondIndex = 1:Config.trialsPerCondition

            % repetition = 0; %to check what its serves, in thier code rep
            % repetition = repetition + 1;

            frameCount = 0; %frame counter

            %all trials per condition paramters of all generated paths -->
            %these could be added to a struct
            allPathsFrameDotXY = []; %its length would be all the xy's off all paths within all trials within a condition
            allPathsStartingPoint = []; %should have the same length of above, and 1 if new path 0 otherwise, ie continuing path --> might be thought of to save indeces instead of a logic array
            allPathsCurvyness = []; %has the factor of curvynes for each xy point
            allPathsDirectionAngle = []; %should have the direction angle of each point xy

            gridCount = [];

            %get a random start angle of direction for both dots
            directionAngle = [Utils.RandAngleDegree(), Utils.RandAngleDegree()];

            pathIndex = 1;
            % for pathIndex = 1:numPaths
            while pathIndex <= numPaths
                
                %compute the curvyness factor for each path, check the
                %function definition for details
                curvynessFactor = curvynessFactor*Utils.ComputeCurvyness(Config.isCurvValenceRand, Config.isCurvFactorRand);

                %get the coordinates of the start point of the path
                if ~Config.isConnectedPaths || pathIndex == 1
                    %if paths are not connected OR if it is the first path
                    % generate a start point
                    %around the center randomly by adding a random offset
                    %to the center coordinates which in term is focused by
                    %a factor around the center
                    frameDotXY = dotRectCenterDegrees + Utils.ComputeOffsetCoordinates(Config.focusAroundCenterFactor, Config.dotRectSize);
                    %add another starting point for the second dot
                    frameDotXY = [frameDotXY dotRectCenterDegrees + Utils.ComputeOffsetCoordinates(Config.focusAroundCenterFactor, Config.dotRectSize)];
                else
                    %if paths should be connected, take the previous
                    %coordinates (ie the last frame)
                    frameDotXY = allPathsFrameDotXY(frameCount, :);
                end

                %both dots should have same path duration so no need to
                %compute twice

                %compute the path variance which they call jitter
                pathDurationVar = rand(1)*Utils.RandPosNeg()*Config.pathDurationVariance*pathDuration;
    
                numFramesPerPath = round((pathDuration + pathDurationVar)*Config.frameFrequency);
               
                if Config.stimulusType == Utils.likelihood  %to be change to reflect the stimulus later on
                    devOnsetVariance = rand(1)*Utils.RandPosNeg()*Config.deviantOnsetVariance*Config.deviantOnset;
                    devFrameOnset = round((Config.deviantOnset + devOnsetVariance)*numFramesPerPath);
                end

                %so here they start the frame loop but above already they
                %got the first point which corresponds to the first frame
                %!!!!! --> can be improved

                outOfBoxPath = 0; %flag to check if the path is to be repeated
                dotsTooClose = 0; %flag to check if the path is to be repeated

                for frameIndex = 1:numFramesPerPath

                    frameCount = frameCount + 1; %why need a frame count apart from index??? --> seems that can be easily replaced with frame index apart from the part before the loop, Updated: since they use more than one path, and the frameIndex is per path, so framCount have all paths

                    gridCount(:, :, frameCount, 1) = zeros(Config.numXGrids, Config.numYGrids);
                    gridCount(:, :, frameCount, 2) = zeros(Config.numXGrids, Config.numYGrids);

                    if frameIndex == 1

                        %check if the dots are too close, and repeat path
                        %if so
                        if norm(frameDotXY(1:2) - frameDotXY(3:4)) < Config.minDistanceBetweenDots
                            disp("dots too close from their starting point, repeating path...");
                            dotsTooClose = 1;
                            break;
                        end

                        allPathsFrameDotXY(frameCount, :) = frameDotXY; %first dot
                        allPathsStartingPoint(frameCount, :) = [1, 1]; %when it starts
                        allPathsDirectionAngle(frameCount, :) = directionAngle;
                        allPathsCurvyness(frameCount, :) = curvynessFactor;

                    else
                                                
                        %get info for rest of the frames
                        %compute the direction change with respect to the
                        %stimulus type
                        if Config.stimulusType == Utils.likelihood

                            %for this type the dir angle change happens on
                            %each frame and gets picked from that normal
                            %distro define below -- this was before for
                            %thiers

                            %Now it is the deviant onset
                            if frameIndex == devFrameOnset && dirVar ~= 0
                                directionAngleChange = trialDeviances(trialPerCondIndex);%[normrnd(0, dirVar), normrnd(0, dirVar)];     
                                fprintf("directionAngleChange: [%d %d]\n", round(directionAngleChange), round(directionAngleChange));
                            else
                                directionAngleChange = 0;
                            end
                        else
                            %for the two other types, only dir change is
                            %aplied on the second frame and then remains
                            %fixed, so zero dir change
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

                        directionAngle = directionAngle + curvynessFactor + directionAngleChange; %in their code it is written that this is not necessary here

                        directionVector = Utils.GetDirectionVector(directionAngle, length(directionAngle));
                        nextDotXY = allPathsFrameDotXY(frameCount-1,:) + directionVector*Config.dotSpeedDegPerFrame;
                        
                        %check if any of the dots hits the boundaries, and
                        %flag for repeating path if so
                        if nextDotXY(1) < round(Config.dotWidth/2, 6) || nextDotXY(3) < round(Config.dotWidth/2, 6)
                            disp("left side hit, repeating path...");
                            outOfBoxPath = 1;
                            break;
                        elseif nextDotXY(2) < round(Config.dotWidth/2, 6) || nextDotXY(4) < round(Config.dotWidth/2, 6)
                            disp("top side hit, repeating path...");
                            outOfBoxPath = 1;
                            break;
                        elseif nextDotXY(1) > Config.dotRectSize(1) - round(Config.dotWidth/2, 6) ...
                                || nextDotXY(3) > Config.dotRectSize(1) - round(Config.dotWidth/2, 6)
                            disp("right side hit, repeating path...");
                            outOfBoxPath = 1;
                            break;
                        elseif nextDotXY(2) > Config.dotRectSize(2) - round(Config.dotWidth/2, 6) ...
                                || nextDotXY(4) > Config.dotRectSize(2) - round(Config.dotWidth/2, 6)
                            disp("bottom side hit, repeating path...");
                            outOfBoxPath = 1;
                            break;
                        end

                        %check if the dots are too close, and repeat path
                        %if so
                        if norm(nextDotXY(1:2) - nextDotXY(3:4)) < Config.minDistanceBetweenDots
                            disp("dots too close, repeating path...");
                            dotsTooClose = 1;
                            break;
                        end

                        allPathsFrameDotXY(frameCount,:) = nextDotXY;
                        allPathsDirectionAngle(frameCount, :) = directionAngle;
                        allPathsCurvyness(frameCount, :) = curvynessFactor;

                        %this one allows to differentiate when the deviant
                        %onset occurs, in theory will also be ones in case
                        %of stimulusType != 0 and frameIndex == 2 since
                        %there is an angle change there
                        if any(directionAngleChange)
                            allPathsStartingPoint(frameCount, :) = [1, 1];
                        else
                            allPathsStartingPoint(frameCount, :) = [0, 0];
                        end
                        

                    end %end else frameIndex == 1

                    %check on which grid does the current frame dot fall,
                    %for both dots
                    gridCount(:, :, frameCount, 1) = Utils.GetFrameInGrid( ...
                        Config.numXGrids, Config.xGridSize, Config.numYGrids, Config.yGridSize, ...
                        allPathsFrameDotXY(frameCount, 1:2));
                    gridCount(:, :, frameCount, 2) = Utils.GetFrameInGrid( ...
                        Config.numXGrids, Config.xGridSize, Config.numYGrids, Config.yGridSize, ...
                        allPathsFrameDotXY(frameCount, 3:4));  

                    %% Transform for each screen visual angle into pixel --> also their code to be reviewed
                    xyDrawDot1(1,1) = allPathsFrameDotXY(frameCount,1)*pixelsPerDegrees + dotRectDegrees(1)*pixelsPerDegrees;
                    xyDrawDot1(1,2) = allPathsFrameDotXY(frameCount,2)*pixelsPerDegrees + dotRectDegrees(2)*pixelsPerDegrees;

                    xyDrawDot2(1,1) = allPathsFrameDotXY(frameCount,3)*pixelsPerDegrees + dotRectDegrees(1)*pixelsPerDegrees;
                    xyDrawDot2(1,2) = allPathsFrameDotXY(frameCount,4)*pixelsPerDegrees + dotRectDegrees(2)*pixelsPerDegrees;

                    Screen('FillRect', win, Config.dotRectColor, dotRectDegrees*pixelsPerDegrees);
                   
                    Screen('DrawDots', win, xyDrawDot1, dotSizePixels, ...
                        Config.dot1Color*255, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    Screen('DrawDots', win, xyDrawDot2, dotSizePixels, ...
                        Config.dot2Color*255, [0 0], 1); % the last two variables are the 'center' and 'dot form'

                    % fix cross
                    Screen('TextSize', win, 40);
                    Screen('TextFont', win, 'Arial');
                    DrawFormattedText(win, '+', 'center', 'center', Config.crossDefaultColor*255);

                    vbl = Screen('Flip', win, [], 1);

                end %end loop of numFramesPerPath

                if outOfBoxPath || dotsTooClose
                                        
                    % pathIndex = pathIndex - 1;
                    frameCount = frameCount - frameIndex;
                    %clean what was saved
                    allPathsFrameDotXY((frameCount+1):end, :) = [];
                    allPathsStartingPoint((frameCount+1):end, :) = [];
                    allPathsDirectionAngle((frameCount+1):end, :) = [];
                    allPathsCurvyness((frameCount+1):end, :) = [];
                else
                    pathIndex = pathIndex + 1;
                end

            end %end loop for numPaths, the while loop
            
            spreadCount(dirVarIndex, pathDurIndex, trialPerCondIndex).gridCount(:,:,:) = sum(gridCount, 3);
            
            %mantaining the same nomenclature of marisas' since this serves
            %as an input to the other script
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).condition          = dirVar;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).PredictionRange    = pathDuration;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).sequence           = trialPerCondIndex;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).xy                 = allPathsFrameDotXY;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).pathAll            = allPathsStartingPoint;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).curvyness          = allPathsCurvyness;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).AngleDirection     = allPathsDirectionAngle;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).gridCounter        = gridCount;
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).gridTot            = spreadCount(dirVarIndex, pathDurIndex, trialPerCondIndex).gridCount(:,:,:);
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).tolleranceGrid     = Config.spreadoutTolerance;
  
            %added for the sake if integration with marisas code
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).borderCounter = zeros(1, length(allPathsFrameDotXY));
            xySeqs(dirVarIndex, pathDurIndex, trialPerCondIndex).borderTot     = 0;

            fprintf("End of trial %d within condition %d (loop trialsPerCondition)\n", trialPerCondIndex, dirVar);
        end %end loop of trialsPerCondition
        
        % disp("End of path duration " + pathDuration);

        totalSpread = zeros(Config.numXGrids, Config.numYGrids);
        for spreadIndex = 1:Config.trialsPerCondition
            totalSpread = totalSpread + spreadCount(dirVarIndex, pathDurIndex, spreadIndex).gridCount;
        end

        if (all(totalSpread <= Config.spreadoutToleranceInterval(2)*Config.trialsPerCondition, "all")) ...
            && (all(totalSpread >= Config.spreadoutToleranceInterval(1)*Config.trialsPerCondition, "all"))
                disp("Spreadout is within the tolerance interval");
                outOfLimitSpread = 0;
        else  
            disp("Spreadout is not within the tolerance interval!!.. repeating");
        end

        end %end while loop of outOfLimitSpread
      
    end %end loop of pathDuration

    xySeqs(dirVarIndex, pathDurIndex).totalSpread = totalSpread;
    fprintf("End of condition [%d] = %d (loop directionVariance)\n", dirVarIndex, dirVar);
end %end loop of directionVariance

%for integration the Cfg struct is needed in their scripts to be removed
%later on and use the Config class instead
Cfg = struct();
Cfg.dpf = Config.dotSpeedDegPerFrame;
Cfg.fps = Config.frameFrequency;
Cfg.Stimulitype = Config.stimulusType;
Cfg.dot_w = Config.dotWidth;
Cfg.rectSize = Config.dotRectSize;
Cfg.DirChange = stimulusTypeConfig.directionChange; %this one only needed in the case of stimulusType!=0

%replicate the trial the amount of repetetions
xySeqs = repmat(xySeqs, 1, Config.trialRepetetion, 1);

save([Config.inputDirectory sprintf(Config.stimuliFileName, subjectID)], 'xySeqs', 'Cfg');
%close all onscreens and offscreens
sca;

%unblock the command prompt inout
%ListenChar(0);

