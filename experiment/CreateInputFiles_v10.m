
%Creates Catch Trial File that is used in other Matlab Script.
% Here we specify the Parameters of the Catch Trials

clear all;
clc;

addpath('lib/'); %AH

Catch.OcclAngleDisplacement = 60; %in degrees

Practice = input('Is this a practice ? [1/0]');
%Practice = 1;

%if Practice == false
    iSub = input('Number of Test Subject: ');
%else
 %   iSub = 2;
%end
Filename = sprintf('MovDot_Sub%02d', iSub); % 'xyRange_v7 - 42P'; %Name of XY Position File

% 0 = random angle; 1 = 90° angle, 2 = radius based on distance between
% disappearing and reappearing position with a ° angle
%Catch.PlacingNewPos = input('Occlusion Task: Should the deviated position start at a fixed position (1) or not (0)?');
Catch.PlacingNewPos = 2;

%% Loading trials
%savedir = 'C:\Users\annik\Desktop\uni\Rovereto\Corsi\AAA_Tirocinio\code\MovImage\InputFiles';
%savedir = 'InputFiles/'; % '\\cimec-storage5\morwur\Projects\MARISA\InputFiles';
savedir = fullfile('.', 'input_files'); % avoid creating ".\input_files" on Unix
if ~exist(savedir, 'dir')
    mkdir(savedir);
end

Dat = load (fullfile(savedir, Filename));

%% Settings
nRuns = 10;
numBlocks = 2; %AH: number of blocks, one per dot

%AH: change to see if it passes
Catch.MaxOverlap = 0.6; %if not a practice run, how much are the catch trials allowed to overlap?

%AH: change the catch per seq
Catch.nCatchPerSeq = 0.125;%1.25; %*16 needs to be a whole number
maxCatchTrials = ceil(Catch.nCatchPerSeq);

%30% are Fixation Cross trials

%AH: putting fixation catch trials to zero
% Catch.CatchRatioFix = 0.3; %how many percentage of Catch Trials are fixation trials?
Catch.CatchRatioFix = 0;

Catch.CatchStart = 0.2; %in s, minimal time after which catch trial is allowed to start

Catch.CatchDistance = 4; %Abstand zwischen Starts von Catch Trials
Catch.BouncingEnd = 0.5; %Duration after Bouncing where Trial is allowed to start

%Fixation Catch Trial
Catch.FixDuration = 0.3; %ins s
Catch.FixResponseDuration = 1; %in s
Catch.FixWaitingDuration = 0.5; %in s
Catch.FixTotalDuration = Catch.FixDuration + Catch.FixResponseDuration+ Catch.FixWaitingDuration;

%AH: change the durations

%Occlusion Catch Trial
Catch.OcclDuration = 0.2; %duration of disappearing
Catch.OcclVideoDuration = 0.3;%0.7; % video afterwards
Catch.OcclResponseDuration = 0.5;%1; %how long to respond
Catch.OcclWaitingDuration = 0.2;%0.5; %how long until trial ends
Catch.OcclTotalDuration = Catch.OcclDuration + Catch.OcclVideoDuration + Catch.OcclResponseDuration + Catch.OcclWaitingDuration;
Catch.CatchRatioOcc = 1-Catch.CatchRatioFix;
Catch.OcclDeviance  = 0.5; %Temporal Deviance in s in Time task
Catch.OcclDistance = round((Dat.Cfg.dpf*Catch.OcclDuration*Dat.Cfg.fps*1));  %Spacial Distance in Visual Angle in Space Task %

Catch.CatchDuration = [Catch.FixTotalDuration, Catch.OcclTotalDuration]; %1 Fixation, 2 Occlusion


%$$$$
%Dat.Cfg.dot_speed = 50;
Catch.Stimulitype = Dat.Cfg.Stimulitype;


if Catch.Stimulitype == 0 %Likelihood
    Catch.BouncingStart = Catch.OcclVideoDuration + Catch.OcclDuration;
else
    Catch.Stimulitype == 1 %PathChange
    Catch.BouncingStart = Catch.OcclVideoDuration + Catch.OcclDuration;
end

%For practice purposes:
if Practice == true
    Catch.MaxOverlap = 1;
    %AH: per seq in practice can be 0.7
    Catch.nCatchPerSeq = 0.7;%2.75; %or 3.25 or 2.5 or 2.75 or 4.0
    maxCatchTrials = ceil(Catch.nCatchPerSeq);
end

% tito = [Dat.xySeqs Dat.xySeqs]; %change XYAll and that's it, also  Dat.xySeqs for getting nSeq correct or actually change xySeq in origin sice used in the experiment one as well

% how many catch trials
nSeq = numel(Dat.xySeqs);
nFrame = size(Dat.xySeqs(1,1).xy, 1); %AH: this one is fixed for them while for me no, so moved to loop

Catch.nCatchPerRun = round(nSeq*Catch.nCatchPerSeq);
Catch.nCatchTotal = Catch.nCatchPerRun* nRuns;

if Practice == true
    %AH: same as below changes 
    % Numberof3 =  (Catch.nCatchPerRun - nSeq)/ (maxCatchTrials - 1);
    % Numberof1 = nSeq - Numberof3;
    % CatchVector = [repmat(1, 1, Numberof1) repmat(maxCatchTrials, 1, Numberof3)];

    Numberof0 = abs(round(nSeq*Catch.nCatchPerSeq) - nSeq);
    Numberof1 = nSeq - Numberof0;
    CatchVector = [repmat(1, 1, Numberof1) repmat(0, 1, Numberof0)];
    
else
    Numberof2 = round(nSeq*Catch.nCatchPerSeq) - nSeq;
    
    %AH: for them numberof2 is which seq should have two catch trials, for
    %me it should be less than 1 the number of catch per sequences (ie
    %trials) so some trial will have zero catches
    Numberof0 = abs(round(nSeq*Catch.nCatchPerSeq) - nSeq);

    %this was theirs
    % Numberof1 = nSeq - Numberof2;
    % CatchVector = [repmat(1, 1, Numberof1) repmat(2, 1, Numberof2)];

    Numberof1 = nSeq - Numberof0;
    CatchVector = [repmat(1, 1, Numberof1) repmat(0, 1, Numberof0)];
end


%Create the Pool of Catch Trials
CatchTrialPool = [repmat(1, 1, round(Catch.CatchRatioFix*nSeq*Catch.nCatchPerSeq)) repmat(2,  1, round(Catch.CatchRatioOcc*nSeq*Catch.nCatchPerSeq))  ]; %1 = Fixation; 2 = Occlusion
CatchOcclusionCongruence = [repmat(0, 1, floor(Catch.CatchRatioFix*nSeq*Catch.nCatchPerSeq)), repmat([1,1],  1, ceil(Catch.CatchRatioOcc*nSeq*Catch.nCatchPerSeq*0.5 ))]; %0 = Fixation, 1= Congruent, 2 = Incongruent

DrawStuff = struct('Condition', [], ...
                    'OldPosition',[],...
                    'CorrectPosition',[],...
                    'NewPosition', []);

%Shuffle the Order of CatchTrials
for iRun = 1:nRuns
    ShuffledCatchMatrix = Shuffle(CatchVector);
    shuffler = Shuffle(1:(Catch.nCatchPerSeq*nSeq)); %shuffel so that congruence is shuffled same as Catch Trials
    ii=0;
    for iSeq=1:nSeq
        for i = 1:ShuffledCatchMatrix(iSeq)
            ii=ii+1;
            ShuffledCatchTrialPool{iRun, iSeq}(i)= CatchTrialPool(shuffler(ii)); %shuffel which kind of catch trial we use
            ShuffledCatchOcclusionCongruence{iRun, iSeq}(i) = CatchOcclusionCongruence(shuffler(ii)); %Also shuffle if occlusion trials are congruent or incongruent
        end
    end
end

%AH: so above has an error that dynamically allocates
%ShuffledCatchTrialPool, so if nothing in the last trial, the cell will not
%reach 80 as dimension and below throws error

%For the Occlusion Trial with a time deviation:
UncongruentStart = cell(nSeq, 1);      %for Congruence, half of the trials have a time lag of 1; the other half of -1
for iSeq = 1:nSeq
    Uncong = (length(find([ShuffledCatchOcclusionCongruence{:,iSeq}] == 2)));
    UncongruentStart{iSeq} = Shuffle([repmat(1, 1, floor(Uncong/2)) repmat(-1, 1, ceil(Uncong/2))]);
end

clear Uncong;
%allSubjectsallRuns = [];


%% Start of Creation

XYAll = Dat.xySeqs(:);

nRestart = 0;

TrialStruct = struct('Type', cell(nRuns, nSeq), ...
    'Start', cell(nRuns, nSeq), ...
    'Condition', cell(nRuns, nSeq), ...
    'PathDuration', cell(nRuns, nSeq), ...
    'VideoNumber', cell(nRuns, nSeq), ...
    'Congruence', cell(nRuns, nSeq), ...
    'IncongruentOffSet', cell(nRuns, nSeq), ...
    'SimulatedPath', cell(nRuns, nSeq), ...
    'SimulatedNewAngle', cell(nRuns, nSeq), ...
    'End', cell(nRuns, nSeq));
TrialStruct(iRun, iSeq).SimulatedPath = cell(1);
TrialStruct(iRun, iSeq).SimulatedNewAngle = cell(1);


%Shuffle the Order of Trials
for iRun = 1:nRuns

    %AH: add another trial order for each block
    for blockIndex = 1:numBlocks
        TrialOrder(blockIndex, :, iRun) = Shuffle(1:nSeq);
    end

    %AH: shuffle order of blocks as well
    BlockOrder(iRun,:) = Shuffle(1:numBlocks);
end

for iSeq = 1:nSeq
    constraints = 0; %while loop until all conditions are met
    CurrentCon = XYAll(iSeq).condition;
    CurrentPathDuration = XYAll(iSeq).PredictionRange;
    i = 0;
    s = 0;

    nFrame = size(XYAll(iSeq).xy, 1); %AH


    while constraints == 0
        i = i+1;
        allRuns = zeros(nRuns, nFrame);
        constraints = 1;
        c = 1;

        for iRun = 1:nRuns
            allIdx = [];
            NoStartBouncing = [];

            if nRestart == 30 %stop if too long
                disp('No solution found, trying again');
                iRun = 1;
                i = 1;
                nRestart = 1;
                Catch.MaxOverlap = Catch.MaxOverlap + 0.1;

                if Catch.MaxOverlap >= 0.5
                    disp('STOP');
                    break
                end

            elseif i == 4000 %maybe problems with placing the catch trials for sequence from the start of all 10 Runs --> restart again
                iRun = 1;
                i = 1;
                disp('restart')
                nRestart = nRestart +1;
                disp(nRestart);
                constraints  = 0;
            end


            nCatch = length(ShuffledCatchTrialPool{iRun, iSeq});  %how many catch trials per sequence?

            %AH: there is a deadlock down in the while var==0, so to avoid
            %it if no catch in this seq continue
            if nCatch == 0
                continue;
            end

            %where does CatchTrial not allowed to start because of
            %bouncing?
            start_bouncing = find(XYAll(iSeq).borderCounter ~= 0) - round(Dat.Cfg.fps*Catch.BouncingStart); %where does bouncing start?
            end_bouncing = find(XYAll(iSeq).borderCounter ~= 0) + round(Dat.Cfg.fps*Catch.BouncingEnd) ; %where does it end


            if Catch.Stimulitype == 0  %likelihood
                % indices = find(XYAll(iSeq).pathAll)'; 
                
                %AH: adding (:,1), for taking one dot for now should be
                %added where necessary
                
                %AH: so the indeces for now, all start from the begining of
                %path and end at its end since only one path unlike theirs

                %AH: updating the indices according to the change in
                %stimuli to add the dev onset as new path
                indices = find(XYAll(iSeq).pathAll(:,1))';
                indices = [indices, length(XYAll(iSeq).pathAll)];
                start_index = indices(1:end-1) + Dat.Cfg.fps*(Catch.CatchStart); %based on this, first possible start position
                end_index = indices(2:end) - Dat.Cfg.fps*Catch.OcclDuration -1; %last possible start position is where the Occlusion can still can occur within the path
                for iPath = 1:numel(start_index)
                    allIdx = [allIdx, start_index(iPath):end_index(iPath)];
                end

            else %for Path Duration
                start_index = Dat.Cfg.fps*(max(Catch.CatchStart, Catch.OcclDeviance))+1; %take the maximum amount
                end_index = nFrame - Dat.Cfg.fps*Catch.CatchDuration(2);
                allIdx = start_index:end_index;
            end


            for iBounce = 1:length(start_bouncing)
                NoStartBouncing = [NoStartBouncing, start_bouncing(iBounce):end_bouncing(iBounce)];
            end

            allIdx(find(allIdx >= (nFrame - Dat.Cfg.fps*Catch.CatchDuration(2)))) = []; %remove all startframes that are too close to end
            %AH: remove the response and waiting time
            % in the other code there is a frame index that gets counted
            % and does allow to go within all states
            % allIdx(find(allIdx >= (nFrame - Dat.Cfg.fps*(Catch.OcclDuration + Catch.OcclVideoDuration) ))) = [];

            allIdx(ismember(allIdx, unique(NoStartBouncing))) = [];
            counter = 0;

            startpos = [];
            var = 0;
            while var == 0
                for iCatch = 1:nCatch
                    catchindexx = ceil(ceil((iCatch-1)*((length(allIdx)+nCatch)/nCatch)) + floor((length(allIdx)+nCatch)/nCatch)*rand(1,1));
                    counter = counter +1;
                    if length(startpos) == nCatch
                        startpos = sort(startpos);
                        var = 1;
                    elseif catchindexx < length(allIdx)
                        minDistSatisfied = all(abs(startpos - allIdx(catchindexx)) >= Catch.CatchDistance * Dat.Cfg.fps);
                        if nCatch > 1 && minDistSatisfied == true
                            startpos = [startpos, allIdx(catchindexx)];
                        elseif nCatch == 1
                            startpos = [allIdx(catchindexx)];
                            var = 1;
                        end
                    else
                        var = 0;
                        if counter > 500000 %maybe problems with placing the catch trials for sequence from the start of all 10 Runs --> restart again
                            counter = 0;
                            disp('restart placing the CatchTrials')
                            startpos = [];
                            var = 0;
                        end
                    end %end of placing catchtrials
                end %end of for loop
            end %end of while loop


            % startpos = sort(randsample(allIdx, nCatch));


            %if only one catch trial, it can be placed anywhere, even at the end, %because there is stil enough time at the end.
            if nCatch > 1 && min(diff(startpos)) < Catch.CatchDistance*Dat.Cfg.fps %minimum distance between Catch Trials
                if mod(i, 200) == 0
                    disp('More Distance needed')
                end
                constraints = 0;
            end


            for iCatch = 1:nCatch
                StartCatch = startpos(iCatch) ; %where does it start?
                %and where does it end?
                %CatchEnd = StartCatch +  round(Catch.CatchDuration(ShuffledCatchTrialPool{iRun,iSeq}(iCatch))*Dat.Cfg.fps);
                CatchEnd =  StartCatch + round(Catch.CatchDuration(ShuffledCatchTrialPool{iRun,iSeq}(iCatch))*Dat.Cfg.fps); %where does Catch Trial end


                % Determine the type of catch trial (fixation or occlusion)
                %1 = Fixation; 2 = Occlusion
                % and the congruency (videostart or not)
                TrialStruct(iRun, iSeq).Type(iCatch) = ShuffledCatchTrialPool{iRun,iSeq}(iCatch);
                TrialStruct(iRun, iSeq).Congruence(iCatch) = ShuffledCatchOcclusionCongruence{iRun, iSeq}(iCatch) ;
                TrialStruct(iRun, iSeq).Condition(iCatch) = CurrentCon;
                TrialStruct(iRun, iSeq).PathDuration(iCatch) = CurrentPathDuration;
                TrialStruct(iRun, iSeq).VideoNumber(iCatch) = iSeq;


                if TrialStruct(iRun, iSeq).Congruence(iCatch) == 2
                    TrialStruct(iRun, iSeq).IncongruentOffSet(iCatch) = UncongruentStart{iSeq}(c); %-1 or 1
                    c = c+1;

                    %% Simulating Trial for Space Task and Angle Task
                    StartPosition = XYAll(iSeq).xy(StartCatch + Dat.Cfg.fps*Catch.OcclDuration,:); %where does Dot reappear

                    if Catch.PlacingNewPos == 1 %if we want a fixed Position
                        %Needs to be fixed first
                        % DevianceDirections = [deg2rad(XYAll(iSeq).AngleDirection((StartCatch+Catch.OcclVideoDuration * Dat.Cfg.fps)) + 90); deg2rad(XYAll(iSeq).AngleDirection((StartCatch+Catch.OcclVideoDuration * Dat.Cfg.fps))- 90)];
                        % NewPositions = StartPosition + [cos(DevianceDirections), sin(DevianceDirections)] .* Catch.OcclDistance;

                    elseif Catch.PlacingNewPos == 0 %if we want a random possible start
                        DevianceDirections = rand(8000, 1) * 2 * pi;  %Create 8000 new positions
                        NewPositions = StartPosition + [cos(DevianceDirections), sin(DevianceDirections)] .* Catch.OcclDistance;

                    elseif Catch.PlacingNewPos == 2 
                        pair = [XYAll(iSeq).xy(StartCatch,:); StartPosition ];

                        % DistanceOccl = pdist(pair, "euclidean");
                        % %AH:first was like this

                        %AH: add (:,1:2) to take just the first dot pair,
                        %and add (:,3:4) also the second dot in other
                        %dimension, or actually no need since they have
                        %same speed so one euclidean will do for both
                        % DistanceOccl = [pdist(pair(:,1:2), "euclidean"), pdist(pair(:,3:4), "euclidean")];
                        DistanceOccl = pdist(pair(:,1:2), "euclidean");

                        distVector = StartPosition - XYAll(iSeq).xy(StartCatch,:);% vector of direcion from point where it disappear to where it should reappear 
                        %distVersor =  distVector/norm(distVector);% the versor from disappearing to reappering point (versor -> norm = 1); in rad
                        
                        %AH:first was like this
                        % Direction = atan2(distVector(2), distVector(1)); % atan2(y, x) returns the angle in radians

                        Direction = [atan2(distVector(2), distVector(1)), atan2(distVector(4), distVector(3))];

                       DevianceDirections1 = [Direction + deg2rad(Catch.OcclAngleDisplacement)]; %in rad
                       DevianceDirections2 = [Direction - deg2rad(Catch.OcclAngleDisplacement)]; %in rad
                        
                       %AH: also here the cos sin should be managed
                       %accordingly, was like this before
                       % NewPosition1 = XYAll(iSeq).xy(StartCatch,:) + [cos(DevianceDirections1), sin(DevianceDirections1)] .* DistanceOccl;
                       % NewPosition2 = XYAll(iSeq).xy(StartCatch,:) + [cos(DevianceDirections2), sin(DevianceDirections2)] .* DistanceOccl;
                       
                       NewPosition1 = XYAll(iSeq).xy(StartCatch,:) + [cos(DevianceDirections1(1)), sin(DevianceDirections1(1)), cos(DevianceDirections1(2)), sin(DevianceDirections1(2))] .* DistanceOccl;
                       NewPosition2 = XYAll(iSeq).xy(StartCatch,:) + [cos(DevianceDirections2(1)), sin(DevianceDirections2(1)), cos(DevianceDirections2(2)), sin(DevianceDirections2(2))] .* DistanceOccl;
                       NewPositions = [NewPosition1; NewPosition2 ];

                    end

                    %AH: also here chane indeces for checking boundaries
                    %with the box for both points (:,1:2:3) instead of
                    %(:,1) and (:,2:2:4) instead of (:,2)
                    ValidIndices = NewPositions(:,1:2:3) >=  round(Dat.Cfg.dot_w/2,6) &  ...
                        NewPositions(:,2:2:4) >= round(Dat.Cfg.dot_w/2,6) & ...
                        NewPositions(:,1:2:3) <= Dat.Cfg.rectSize(1)-round(Dat.Cfg.dot_w/2,6) & ...
                        NewPositions(:,2:2:4) <= Dat.Cfg.rectSize(2)-round(Dat.Cfg.dot_w/2,6) ;% if dot is touching rect border

                    jj = 0;

                    if max(ValidIndices) ~= 1
                        if mod(i, 2000) == 0
                            disp('New Position needed')
                        end

                        if Catch.PlacingNewPos == 1
                            if mod(i, 2000) == 0
                                disp('Start Position would be outside of Boundaries')
                            end
                            constraints = 0;
                        elseif Catch.PlacingNewPos == 0
                            constraints = 0; %if we want a random direction on a circle
                        end
                        break;


                    else
                        %AH: also here, i have to find the valid position
                        %which not necessarly are on the same row for both
                        % ValidPositionIndex = Shuffle(find(ValidIndices == 1));

                        ValidPositionIndexDot1 = Shuffle(find(ValidIndices(:,1)));
                        ValidPositionIndexDot2 = Shuffle(find(ValidIndices(:,2)));
                        NewPosition = [NewPositions(ValidPositionIndexDot1(1), 1:2), ...
                            NewPositions(ValidPositionIndexDot2(1), 3:4)];

                        xyReplace = [];
                        xySimulated = [];
                        jj = 0;

                    end


                    % OcclVideoDuration was before, should be
                    % OcclDuration, since the strating point after
                    % occlusion i suppose
                    for iMov = (StartCatch+Catch.OcclDuration * Dat.Cfg.fps):(StartCatch+(Catch.OcclDuration + Catch.OcclVideoDuration)* Dat.Cfg.fps)
                        currentPath = XYAll(iSeq).pathAll(iMov,:); %AH: add (iMov,:) instead of (iMov)
                        jj = jj+1;

                        if Catch.Stimulitype == 0
                            NewDir = 0;
                            mu = 0;%normrnd(0, CurrentCon); %AH: commented this since we need the catch to mantain direction
                        else
                            if currentPath == 1
                                if  Catch.Stimulitype == 1
                                    NewDir = rand(1,1)*360;
                                elseif Catch.Stimulitype == 2
                                    NewDir = normrnd(0, Dat.Cfg.DirChange);
                                end

                            elseif currentPath == 0
                                NewDir = 0;
                            end
                            mu = 0;

                        end

                        % OcclVideoDuration was before, should be
                        % OcclDuration, since the strating point after
                        % occlusion i suppose
                        if iMov == StartCatch+Catch.OcclDuration * Dat.Cfg.fps
                            xyReplace(jj,:) = NewPosition; %for the spatial task
                            AngleDirection = XYAll(iSeq).AngleDirection(iMov,:); %AH: add (iMov,:) instead of (iMov)

                            xySimulated(jj,:) = StartPosition; %for the angle task
                            NewAngle = XYAll(iSeq).AngleDirection(iMov,:) + rand(1,1)*360; %AH: add (iMov,:) instead of (iMov)

                        else
                            %AH: add (iMov,:) instead of (iMov)
                            AngleDirection = AngleDirection +  XYAll(iSeq).curvyness(iMov,:) + NewDir + mu; %Calculate new simulated Position of Dot
                            NewAngle = NewAngle +  XYAll(iSeq).curvyness(iMov,:) + NewDir + mu;


                            % if xyReplace(jj-1,1) < (Dat.Cfg.ppd*Dat.Cfg.sizeDotInPixel)/2 || xyReplace(jj-1,2) < (Dat.Cfg.ppd*Dat.Cfg.sizeDotInPixel)/2 || xyReplace(jj-1,1) > Dat.Cfg.rectSize(1)-(Dat.Cfg.ppd*Dat.Cfg.sizeDotInPixel)/2 || xyReplace(jj-1,2) > Dat.Cfg.rectSize(2)-(Dat.Cfg.ppd*Dat.Cfg.sizeDotInPixel)/2 % if dot is touching rect border
                            %     disp (' pleessszzz');
                            % end

                            % [xyReplace(jj,:), borderCounter(jj), AngleDirection] = createNewPoint(xyReplace(jj-1,:), Dat.Cfg.dpf, AngleDirection,  Dat.Cfg.rectSize ,  Dat.Cfg);
                            % [xySimulated(jj,:), borderCounterSimulated(jj), NewAngle] = createNewPoint(xySimulated(jj-1,:), Dat.Cfg.dpf, NewAngle,  Dat.Cfg.rectSize,  Dat.Cfg);

                                                    
                            directionVector = Utils.GetDirectionVector(AngleDirection, length(AngleDirection));
                            xyReplace(jj,:) = xyReplace(jj-1,:) + directionVector*Dat.Cfg.dpf;
                            directionVector = Utils.GetDirectionVector(NewAngle, length(NewAngle));
                            xySimulated(jj,:) = xySimulated(jj-1,:) + directionVector*Dat.Cfg.dpf;
                            
                            NextPositions = [xyReplace(jj,:); xySimulated(jj,:)];
                            ValidNextPositions = NextPositions(:,1:2:3) >=  round(Dat.Cfg.dot_w/2,6) &  ...
                                NextPositions(:,2:2:4) >= round(Dat.Cfg.dot_w/2,6) & ...
                                NextPositions(:,1:2:3) <= Dat.Cfg.rectSize(1)-round(Dat.Cfg.dot_w/2,6) & ...
                                NextPositions(:,2:2:4) <= Dat.Cfg.rectSize(2)-round(Dat.Cfg.dot_w/2,6) ;% if next dots are touching rect border
                                                
                            if ~all(ValidNextPositions, 'all')
                                disp("the NextPositions hit the boundaries!");
                                constraints = 0;
                                break;
                            end
                            
                        end
                    end
                    TrialStruct(iRun,iSeq).SimulatedPath{iCatch} = xyReplace; %for Task with Spatial Deviation
                    TrialStruct(iRun, iSeq).SimulatedNewAngle{iCatch} = xySimulated; %For Task with a new Angle


                else


                    xyReplace = [];
                    jj = 0;

                    yTweak = 0.04;
                    for iMov = (StartCatch):(StartCatch+(Catch.OcclDuration)* Dat.Cfg.fps)
                        jj = jj+1;

                        yTweak = -yTweak;

                        originalCoords = XYAll(iSeq).xy(iMov, :);
                        xyReplace(jj,:) = originalCoords + [0 yTweak 0 yTweak];


                        ValidNextPositions = xyReplace(jj,1:2:3) >=  round(Dat.Cfg.dot_w/2,6) &  ...
                                xyReplace(jj,2:2:4) >= round(Dat.Cfg.dot_w/2,6) & ...
                                xyReplace(jj,1:2:3) <= Dat.Cfg.rectSize(1)-round(Dat.Cfg.dot_w/2,6) & ...
                                xyReplace(jj,2:2:4) <= Dat.Cfg.rectSize(2)-round(Dat.Cfg.dot_w/2,6) ;% if next dots are touching rect border
                                                
                        if ~all(ValidNextPositions, 'all')
                            disp("the NextPositions hit the boundaries!");
                            constraints = 0;
                            break;
                        end
                    end

                    TrialStruct(iRun, iSeq).IncongruentOffSet(iCatch) = 0;
                    TrialStruct(iRun, iSeq).SimulatedPath{iCatch} = xyReplace;
                    TrialStruct(iRun, iSeq).SimulatedNewAngle{iCatch} = 0;
                end

                %Record the Start and End of Trials
                TrialStruct(iRun, iSeq).Start(iCatch) = StartCatch;
                TrialStruct(iRun, iSeq).End(iCatch) = CatchEnd;
                allRuns(iRun,StartCatch:CatchEnd) = 1;
            end
        end
        % allSubjectsallRuns = [allSubjectsallRuns, allRuns] ;

        %end

        %AH: was before greter or equal

        if max(mean(allRuns))> Catch.MaxOverlap %0.2 not possible
            constraints=0;
        end

        if mod(i,2000) == 0
            %clc;
            disp('iSeq:')
            disp(iSeq);
            disp('iRun:')
            disp(iRun);
            disp(i);
            % imagesc(allRuns);
        end
    end
end



filename1 = sprintf('Sub%02d_TrialStruct',iSub);
%filename3 = 'Practice_TrialStruct';
filename3 = convertStringsToChars(append('Practice_', Filename));


if Practice == true
    Practice_Struct = fullfile(savedir, [filename3 '.mat']);
    
        save(Practice_Struct, 'TrialStruct', 'Catch', 'TrialOrder', 'BlockOrder');
        disp(['Variables saved to file: ', Practice_Struct]);
    
     
else
    fullPath_struct = fullfile(savedir, [filename1 '.mat']);
    
    if exist(fullPath_struct, 'file') == 2
        % File already exists, do not overwrite
        disp('File already exists. Not overwriting.');
    else
        % File does not exist, proceed with saving
        save(fullPath_struct, 'TrialStruct', 'Catch', 'TrialOrder', 'BlockOrder');
        disp(['Variables saved to file: ', fullPath_struct]);
    end
    
end



%clear all

