%% Clear everything

clear all;
clc;
sca;
format shortG

%% Values for MEG, EyeTracker, PracticeTrials, Debugging etc.
Conf.MEG =0;
Conf.trackeye = 0;
Conf.practice = 1;
Conf.debug = 0;
Conf.trialAmountPractice = 5;%16;
Conf.practiceDemoIndex = 10;

Conf.experiment_name = 'MoveDot1'; %% Experiment Data
if Conf.practice
    Conf.experiment_name = 'practice';
end

desiredFrameRate = 60; % independent of screen resolution
Conf.Task  = 1;


if Conf.MEG
    Conf.trackeye = true;
end
if Conf.trackeye
    Screen('Preference',  'SkipSyncTests', 0);
end


%% Input Data

% Create a dialog window for user input
prompt = {'Subject Number:', 'Real Block Number:', 'Viewing Distance (mm):'};
dlgtitle = 'User Input';
dims = [1 50]; % dimensions of the input boxes

% Set default values for subject number and practice run
definput = {'20', '1', '1000'}; % default values

user_input = inputdlg(prompt, dlgtitle, dims, definput);

% Convert the entered values
iSub = str2double(user_input{1});
iRun = str2double(user_input{2});
viewing_distance_mm = str2double(user_input{3});

% Check if the entered values are valid
if isnan(iSub) || isnan(Conf.practice) || isnan(Conf.Task) || isnan(iRun) || isnan(Conf.debug) || isnan(viewing_distance_mm) || viewing_distance_mm <= 0
    error('Invalid input. Please enter valid values.');
end

NameFile = sprintf('MovDot_Sub%02d',iSub); %'MovDot_Pilot_Sub1';

% Assign the entered values to the configuration variables
Conf.display.dist = viewing_distance_mm;


%% Load the stimulus structure

rootdir = '.';

savedir = [rootdir '/output_files'];
Inputdir = [rootdir '/input_files'];

Dat = load(fullfile(Inputdir, NameFile));


%% Load Trial Infos, including catch trials

%For CatchTrials:
if Conf.practice == true
    load (fullfile (Inputdir,['Practice_' NameFile] ));
    
else %For normal experiment:
    load (sprintf ('%s/Sub%02d_TrialStruct.mat', Inputdir,iSub));
end


%% Colors & Lines : For Videos

%AH: this one they use it to print Missed and also lines in practice??
Conf.colors = [27 158 119;217 95 2;117 112 179;231,41,138;102,166,30];%RGB values of colors that are colorblind and grayscale friendly

%AH: all of these should come from config

Conf.background  = 30; %Color of Screen
Conf.RectColor = [0 0 0]; % Color of Rectangle in which Dot moves;

Conf.ColorOval = [180, 180, 180]; %Color of the Oval

Conf.CrossColor = [10 10 10]; %Color of Fixation Cross
Conf.CrossColorChange = [50 50 250]; %[190 230 230]; %Color when Fixation Cross changes Color

%white and black are later

%% Experiment Settings

%AH: here also from config

Conf.refrate = Dat.Cfg.fps; %set to 60 because Marisa's Laptop can't handle more, also annika laptop
Conf.rectSize_visualangle   = Dat.Cfg.rectSize;

Conf.display.rect = []; %full screen if we are not debugging
Conf.debugRect = [1 1 800 800] ; %size of debug window

%AH: check this one if zero what happens
Conf.resetTime = 0;0.2; % how much time do we want to go back after a catch trial (in sec)

%% Behavioral Part

KbName('UnifyKeyNames');
Key.spaceKey = KbName('space');
Key.escapeKey = KbName('ESCAPE');
Key.LeftKey = KbName('LeftArrow');
Key.RightKey = KbName('RightArrow');
Key.EnterKey = KbName('Return');
Key.V = KbName('v');
Key.C = KbName('c');

RestrictKeysForKbCheck([Key.escapeKey Key.spaceKey Key.LeftKey Key.RightKey Key.EnterKey Key.V Key.C]);

%% General

% PSYCHTOOLBOX
%PsychDefaultSetup(2);
% if Conf.MEG
%     win = Screen('OpenWindow',1);
% end
Screen('Preference', 'SkipSyncTests', 1);%do screen sync + size tests

% For real experiment hide cursor
if Conf.debug == false
    %     HideCursor;
    % PsychDebugWindowConfiguration
elseif Conf.debug == true
    Conf.display.rect = Conf.debugRect;
    PsychDebugWindowConfiguration([],.8); %for debugging show transparent
end

if Conf.MEG==1
    screenNumber = 1; %in the lab, we use the left screen for the stimulus computer
else
    screens = Screen('Screens');
    screenNumber = max(screens); %for the laptop
end

Conf.screenSize = Screen('Rect', screenNumber); %shows the Pixel Number

[Conf.mon_width, Conf.mon_height] = Screen('DisplaySize', screenNumber); %in mm, show me the width and height of screen, for visual angle

Conf.ppd = pi * (Conf.screenSize(3)-Conf.screenSize(1)) / atan(Conf.mon_width/Conf.display.dist/2) / 360; % pixels per degree %% transforming pixelsize into visual angles
Conf.mpd =  2*Conf.display.dist*tan(deg2rad(1)/2); % millimeters per degree angle


%% Transforming from Visual Angle into Pixel

Conf.rectSize = floor(Conf.rectSize_visualangle * Conf.ppd);
Conf.sizeDotInPixel = floor(Dat.Cfg.dot_w*Conf.ppd);

Conf.FixSize = floor(0.2*Conf.ppd); %Pixels of Fixation Cross %before: (0.24697

%Information about presented Text
Conf.Textsize = round(Conf.FixSize*1.5);
Conf.Text = 'Arial';

%AH settings
Conf.warningQstMarkColor = Conf.colors(2,:);%[255 0 0];
Conf.normalQstMarkColor = 100;%[255 0 0];
Conf.showQstMarkOnOcclusion = 0;
Conf.qstMarkSizeFactor = 1.5;

Conf.jitterCatches = 1;

Conf.DotColor1 = [0 255 0]; %Color of Dot1
Conf.DotColor2 = [255 0 0]; %Color of Dot2

Conf.showFixCrossBetweenTrials = 1;

Conf.blockConditions = [1, 2]; %1: GREEN, 2: RED
Conf.block1Msg = ['New Block \n\n Pay attention to the GREEN dot \n\n'];
Conf.block2Msg = ['New Block \n\n Pay attention to the RED dot \n\n'];
Conf.blockMsgColor = 100;

%% Fixation Cross / Bullseye Fixation Cross

Conf.xCoords = [-Conf.FixSize Conf.FixSize 0 0];
Conf.yCoords = [0 0 -Conf.FixSize Conf.FixSize];
Conf.allCoords = [Conf.xCoords; Conf.yCoords]; %Coordinates of Fixation Cross
Conf.lineWidthPix = Conf.FixSize/3; %lineWidth of Fixation Cross

Conf.RadiusOut = Conf.FixSize; % the bigger Oval around Fixation Cross
Conf.RadiusIn = Conf.lineWidthPix /1.5; %the small oval within Fixation Cross


%% Open the experiment windowcenter
[window, windowRect] = Screen('OpenWindow', screenNumber, Conf.background, Conf.display.rect);

% Get screen info
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[xCenter, yCenter] = RectCenter(windowRect);% get screen center

Conf.black = BlackIndex(screenNumber); %black color
Conf.white = WhiteIndex(screenNumber); %white color

Conf.rectCoords  = [xCenter-Conf.rectSize(1)/2 yCenter-Conf.rectSize(2)/2 xCenter+Conf.rectSize(1)/2 yCenter+Conf.rectSize(2)/2];
%AH: settings
Conf.marginRect = Conf.rectCoords + [-10 -10 10 10];

if ~Conf.MEG
    [minsmooth,maxsmooth] = Screen('DrawDots', window,[0 0]);
    Conf.sizeDotInPixel = min(max(Conf.sizeDotInPixel, minsmooth), maxsmooth);
    clearvars minsmooth maxsmooth %to clean up workspace
end

[screenXpixels, screenYpixels] = Screen('WindowSize', window);%get screen size Pixels --> for EyeTracker


%Starting Experiment
DrawFormattedText(window,...
    'Loading videos...',...
    'center', yCenter-Conf.rectSize(2)*(3/8), Conf.black);
Screen('Flip', window);


%%
% Set PTB to top priority so no other running processes on this PC interfer
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);


%% Load the information of Stimuli and XY Position of Dot

nSeq = numel(Dat.xySeqs);
nFrame = size(Dat.xySeqs(1,1).xy, 1);   %AH: this one is not fixed for me
Conf.TotalnFrames =  nFrame;
xyAll = Dat.xySeqs (:);


if Dat.Cfg.Stimulitype == 0
    nCond = length(unique([xyAll.condition]));
    
else
    nCond = length(unique([xyAll.PredictionRange]));
end



%AH: comment this one and leave dynamic for now since nFrame is not fixed
% stimuli = zeros(nSeq, nFrame, 2); % Assuming XY positions are 2D
condMatrix = zeros(4, nSeq);

for iSeq = 1:nSeq
    %condMatrix(:, iSeq) = [xyAll(iSeq).condition; xyAll(iSeq).sequence; iSeq];
    condMatrix(:, iSeq) = [xyAll(iSeq).condition; xyAll(iSeq).sequence; iSeq; xyAll(iSeq).PredictionRange];
    
    nFrame = size(xyAll(iSeq).xy, 1);
    for iframe = 1:nFrame
        %stimuli (iSeq, iframe, :)  = xyAll(iSeq).xy (iframe, :);
        stimuli (iSeq, iframe, 1)  = xyAll(iSeq).xy (iframe, 1)*Conf.ppd  +Conf.rectCoords(1);  %From Visual Angle into Pixels
        stimuli (iSeq, iframe, 2)  = xyAll(iSeq).xy (iframe, 2)*Conf.ppd   +Conf.rectCoords(2); %Also, adjust the Position of the Dot according to the rectanlge it's in
        %AH: second dot 
        stimuli (iSeq, iframe, 3)  = xyAll(iSeq).xy (iframe, 3)*Conf.ppd  +Conf.rectCoords(1);  %From Visual Angle into Pixels
        stimuli (iSeq, iframe, 4)  = xyAll(iSeq).xy (iframe, 4)*Conf.ppd   +Conf.rectCoords(2); %Also, adjust the Position of the Dot according to the rectanlge it's in
    end
end


%AH: MOVED down

% Randomise the order of Videos presented

%AH: so the trial order is generated else where and is used as the random
%indeces of the order, a run is made of 16 trials in her case

%AH: our case, we have 20 trials per cond, her case was 4 trial per 4 conds, so
%we have 40, she has 16

% trialOrder = TrialOrder(iRun,:);        %AH: is loaded from --> load (sprintf ('%s/Sub%02d_TrialStruct.mat', Inputdir,iSub));
% condMatrixShuffled = condMatrix(:, trialOrder); %shuffel the conditions of matrix
% clear condMatrix

%AH: this one something related to datapix
TriggerValues = [1:16, 17, 21:36, 41:56, 60, 99, 101, 113, 201, 202];


%% Set up Eye Link

if Conf.trackeye
    % INITIALIZE EYELINK
    % It is better not to send too many Eyelink commands to the eye-tracker in a row. For this reason, between them, we wait for a short time, here defined.
    elk.wait = 0.01;
    
    % This code initializes the connection with the eyelink: if something fails, it exit program with error
    if EyelinkInit()~= 1
        error('Eyelink disconnected !!!');
    end;
    
    % We need to provide Eyelink with details about the graphics environment and perform some initializations. The initialization information is returned in a
    % structure that also contains useful defaults and control codes (e.g. tracker state bit and Eyelink key values). The structure, moreover, act asn an handle
    % for subsequent commands, like "windowHandle" for Psychtoolbox.
    elk.el = EyelinkInitDefaults(window);
    
    % Here we create the name for the eyelink datafile. Data gathered from the eye tracker are saved on the eye-tracking PC in a file. Data from all users are
    % saved in the same folder and the folder is routinely cleaned up without any advice. So, be sure to copy your data after the experiment and choose an
    % unique name for the datafile (containing date/time, subject number etc...). It has to be less than 8 characters long.
    
    if iSub < 10
        Eyefilename = ['mbMD0' num2str(iSub) 'b' num2str(iRun) '.edf'];
    else
        Eyefilename = ['mbMD' num2str(iSub) 'b' num2str(iRun) '.edf'];
    end
    
    %Eyefilename = Eyefilename.edf;
    elk.edfFile = sprintf(Eyefilename);		% Create file name
    Eyelink('Openfile', elk.edfFile);									% Open the file to the eye-tracker
    
    % Writing a short preamble to the file helps if the name became not that informative ;-)
    Eyelink('command', sprintf('add_file_preamble_text ''Marisa Birk Project MovImage: subject %d ; block %d ; practice %d ; time %s''', iSub, iRun, Conf.practice, datestr(now, 'YYYYmmddhhMM')));
    
    % Setting the eye-tracker so as to record GAZE of  LEFT and RIGHT eyes, together with pupil AREA
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    
    % Setting the proper recording resolution, proper calibration type, as well as the data file content
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
    
    % Setting the proper calibration type. Usually we use 9 points calibration. For a long range mount also 13 points (HV13) is a good (longer) calibration.
    Eyelink('command', 'calibration_type = HV9');
    
    % Setting the proper data file content
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');
    
    % Setting link data (used for gaze cursor, optional)
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');
    
    % Saccade detection thresholds (optional)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    
    % Now make sure that we are still connected to the Eyelink ... otherwise throw error
    if Eyelink('IsConnected')~=1
        error('Eyelink disconnected !!!');
    end
    
    % EYELINK CALIBRATION
    % This code allow the EyeLink software to take control of your psychtoolbox screen. This means that at this point you will see participant eyes as recorded
    % by the Eye-tracker camera on the MEG whiteboard, a condition essential for setting up the camera. After setting up the camera you can perform calibration
    % and validation at this step.
    
    % Some calibration parameters
    elk.el.foregroundcolour = 0;
    elk.el.backgroundcolour = [Conf.background Conf.background Conf.background];
    
    
    % Give eye-tracker control of the screen for camera setup and calibration, until you exit back to psychtoolbox by pressing ESC
    EyelinkDoTrackerSetup(elk.el);
end

if Conf.MEG
    % INITIALIZE DATAPIXX
    Datapixx('Open');					% Open DataPixx
    
    Datapixx('SetVideoMode', 0);		% This set video mode to normal passthrought, no stereo mode. C24, Straight passthrough from DVI 8-bit RGB to VGA RGB.
    % In this configuration luminance is linear with RGB (see our wiki).
    
    Datapixx('StopAllSchedules');		% Stop all schedules (audio waveforms, triggers etc...)
    
    Datapixx('SetDoutValues', 0);		% Set digital output to zero, as required to prepare for triggering
    
    Datapixx('EnableDinDebounce');		% Enable response debouncing. This is required to prune out spurious button presses after a real response
    
    Datapixx('SetDinLog');				% Clear digital input logger, i.e: clear old responses in the register
    Datapixx('StopDinLog');				% Stop running response logger
    
    Datapixx('RegWrRd');				% So far, no real changes occurred on the physical blue box devices. This command synchronize local and remote registers
    % in a read/write mode and immediately. Only now, the blue box status is as determined by the above initializations.
    
    responseButtonsMask = 2^0 + 2^1 + 2^2 + 2^3;	% Values of response buttons are stored
    % %in a cumbersome binary way. This is a binary mask useful to
    % transform them in decimal human-readable values.
    % %In particular, red = 1, blue = 8, geen = X and yellow =
    % X. It works. Just believe it. I do, I am a true believer. Neo is the one.
    
end




%% Practice Trial

if Conf.practice == true    
    
    nSeq = Conf.trialAmountPractice;   %reduce the number of trials to practice number

    if Conf.MEG
        Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
        
        % Reset and fire up the response logger
        Datapixx('SetDinLog');
        Datapixx('StartDinLog');
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');
        
        while status.newLogFrames == 0  %1
            Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
            Screen('TextSize', window, Conf.Textsize);
            Screen('TextFont',window, Conf.Text);
            DrawFormattedText(window,...
                'Thank you for participating in this experiment',...
                'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
            DrawFormattedText(window,'- Press space to start the instructions -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
            
            Screen('Flip', window);
            
            Datapixx('RegWrRd');                        % Commit changes to/from DP
            status = Datapixx('GetDinStatus');			% Get response logger status
        end
        
        WaitSecs(.5);
        Datapixx('StopDinLog');
        
        Datapixx('EnableDinDebounce');
        Datapixx('SetDinLog');
        Datapixx('StartDinLog');
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        
        while status.newLogFrames == 0 %2
            
            Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
            Screen('TextSize', window, Conf.Textsize);
            Screen('TextFont',window, Conf.Text);
            DrawFormattedText(window,...
                'Each trial starts with a fixation cross \n\n Whenever you see a fixation cross, please fixate it',...
                'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
            DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut], Conf.RadiusOut*2 )
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
            Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
            
            Screen('Flip', window);
            
            Datapixx('RegWrRd');                        % Commit changes to/from DP
            status = Datapixx('GetDinStatus');			% Get response logger status
            
        end
        
        WaitSecs(.5);
        Datapixx('StopDinLog');
        
        Datapixx('EnableDinDebounce');
        Datapixx('SetDinLog');
        Datapixx('StartDinLog');
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        
        while status.newLogFrames == 0 %3
            
            Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
            Screen('TextSize', window, Conf.Textsize);
            Screen('TextFont',window, Conf.Text);
            DrawFormattedText(window,...
                ['Next you will see videos of two moving dots of distinct colors \n\n At the beginning of each block \n\n' ...
                ' You will be instructed to pay attention to one of the moving dots \n\n' ...
                'Pay attention to the instructed dot \n\n But keep fixating the cross while doing so'],...
                'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
            DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
            % AH: stimuli(1,10,1:2) instead of stimuli(1,10,:)
            Screen('DrawDots', window, squeeze(stimuli(2,Conf.practiceDemoIndex,1:2))', Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
            Screen('DrawDots', window, squeeze(stimuli(2,Conf.practiceDemoIndex,3:4))', Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1); % the last two variables are the 'center' and 'dot form'
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
            Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
            
            Screen('Flip', window);
            
            Datapixx('RegWrRd');                        % Commit changes to/from DP
            status = Datapixx('GetDinStatus');			% Get response logger status
            
        end
        
        %AH: if fixation trials exist, show their instructions
        if Catch.CatchRatioFix > 0

        WaitSecs(.5);
        Datapixx('StopDinLog');
        
        Datapixx('EnableDinDebounce');
        Datapixx('SetDinLog');
        Datapixx('StartDinLog');
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        
        while status.newLogFrames == 0 %4
            
            Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
            Screen('TextSize', window, Conf.Textsize);
            Screen('TextFont',window, Conf.Text);
            DrawFormattedText(window,...
                'Sometimes the fixation cross \n\n shortly changes its color \n\n You have to indicate the color change \n\n by pressing the left arrow key',...
                'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
            DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
            % AH: stimuli(1,10,1:2) instead of stimuli(1,10,:)
            Screen('DrawDots', window, squeeze(stimuli(2,Conf.practiceDemoIndex,1:2))', Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
            Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColorChange, [xCenter yCenter], 2);
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
            
            Screen('Flip', window);
            
            Datapixx('RegWrRd');                        % Commit changes to/from DP
            status = Datapixx('GetDinStatus');			% Get response logger status
            
        end

        end
        
        WaitSecs(.5);
        Datapixx('StopDinLog');
        
        Datapixx('EnableDinDebounce');
        Datapixx('SetDinLog');
        Datapixx('StartDinLog');
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        
        while status.newLogFrames == 0 %5
            
            Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
            Screen('TextSize', window, Conf.Textsize);
            Screen('TextFont',window, Conf.Text);
            DrawFormattedText(window,...
                'Sometimes the dots will \n\n move shortly in a tweaky manner. \n\n \n\n  You have to indicate \n\n  whenever the instructed dot \n\n shows such a tweak. \n\n \n\n By pressing the left arrow key',...
                'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
            DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
            
            Screen('Flip', window);
            
            Datapixx('RegWrRd');                        % Commit changes to/from DP
            status = Datapixx('GetDinStatus');			% Get response logger status
            
        end
        
        WaitSecs(.5);
        Datapixx('StopDinLog');
        
        Datapixx('EnableDinDebounce');
        Datapixx('SetDinLog');
        Datapixx('StartDinLog');
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Please remember to relax and fixate the cross',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- When ready, \n\n inform the experimenter -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('Flip', window);
        KbStrokeWait;
        
        
    else %if not within MEG
        %1
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Thank you for participating in this experiment',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- Press space to start the instructions -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('Flip', window);
        KbStrokeWait;
        
        %2
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Each trial starts with a fixation cross \n\n Whenever you see a fixation cross, please fixate it',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
        Screen('Flip', window);
        KbStrokeWait;
        
        %3
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            ['Next you will see videos of two moving dots of distinct colors \n\n At the beginning of each block \n\n ' ...
            'You will be instructed to pay attention to one of the moving dots \n\n' ...
            'Pay attention to the instructed dot \n\n But keep fixating the cross while doing so'],...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        % AH: stimuli(1,10,1:2) instead of stimuli(1,10,:)
        Screen('DrawDots', window, squeeze(stimuli(2,Conf.practiceDemoIndex,1:2))', Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
        Screen('DrawDots', window, squeeze(stimuli(2,Conf.practiceDemoIndex,3:4))', Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1); % the last two variables are the 'center' and 'dot form'
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
        Screen('Flip', window);
        KbStrokeWait;
        
        %AH: if fixation trials exist, show their instructions
        if Catch.CatchRatioFix > 0

        %4
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Sometimes the fixation cross \n\n shortly changes its color \n\n You have to indicate the color change \n\n by pressing the left arrow key',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        % AH: stimuli(1,10,1:2) instead of stimuli(1,10,:)
        Screen('DrawDots', window, squeeze(stimuli(2,Conf.practiceDemoIndex,1:2))', Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
        
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColorChange, [xCenter yCenter], 2);
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
        Screen('Flip', window);
        KbStrokeWait;
        
        end

        %5
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Sometimes the dots will \n\n move shortly in a tweaky manner. \n\n \n\n  You have to indicate \n\n  whenever the instructed dot \n\n shows such a tweak. \n\n \n\n By pressing the left arrow key',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- Press space to continue -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('Flip', window);
        KbStrokeWait;
        
        % 6
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Please remember to relax and fixate the cross',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- When ready, \n\n press space to start with some practice trials -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('Flip', window);
        KbStrokeWait;
    end
    
else %no Instructions needed if not a practice trial
    %Screen('DrawTexture', window, emptyTex);
    
    if Conf.MEG
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'I will now measure \n\n the position of your head \n\n Please do not move your head \n\n anymore from now on \n\n Please remember \n\n to relax and fixate the cross',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- When ready, inform the experimenter -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('Flip', window);
        KbStrokeWait;
        
    else
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        Screen('TextSize', window, Conf.Textsize);
        Screen('TextFont',window, Conf.Text);
        DrawFormattedText(window,...
            'Please remember \n\n to relax and fixate the cross',...
            'center', yCenter-Conf.rectSize(2)*(3/8), Conf.white);
        DrawFormattedText(window,'- When ready, \n\n press space to start with the next block -','center', yCenter+Conf.rectSize(2)*(3/8), Conf.white);
        Screen('Flip', window);
        KbStrokeWait;
    end
    
end

%% Output File

%Nr = condMatrixShuffled(3,:);

%Subject Nr
% Number of Block / of Rrun
%Nr of Trial of Test subject
% Number of Video per Condition (1-4)
%Condition (Predictability)
% PathDuration (Predictability)
%Video_nr: 1-16 video
%numcatch = number of catch trials
%startcatch: framenumber in which catchtrial starts
%endcatch: first frame after catch trial ends
%catch_type: 1 = fixation; 2 = occlusion
%congruence: 0 = fixation, 1 = congruent, 2 = incongruent
%Response: 1 = left arrow key; 2= right arrow key / blue button vs red
%button
%Correct_response: 1 = correct; 0 = false
%RT: Reaction Time
%testTiming: ~ 30s long
%testITI: how long is ITI

clear output;
output = struct(...
    'subject_num', num2cell(iSub * ones(1, nSeq)), ...
    'run_num', num2cell(iRun * ones(1, nSeq)), ...    %AH: was block_num
    'block_cond', [], ...
    'trial_num', [], ...
    'sequence', [], ...    %AH: 1:nSeq since nSeq in practice is changes and first was defined before it
    'condition', [], ...
    'PathDuration', [], ...
    'predictability',[], ...
    'video_nr', [], ...
    'numcatch', [], ...
    'startcatch', [], ...
    'endcatch', [], ...
    'catch_type', [], ...
    'congruence', [], ...
    'lagtime', [], ...
    'response', [], ...
    'correct_response', [], ...
    'RT', [], ...
    'testTiming', [], ...
    'testITI', [], ...
    'FrameDuration', [], ...
    'SequenceStartTime', [], ...
    'SequenceEndTime', [], ...
    'CatchStartTime', [], ...
    'CatchEndTime', [], ...
     'XYPositionPerFrame', [] );

%% Starting Experiment

DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(3)], 'center', 'center', Conf.black);
Screen('Flip', window);
WaitSecs(1);
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(2)], 'center', 'center', Conf.black);
Screen('Flip', window);
WaitSecs(1);
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(1)], 'center', 'center', Conf.black);
Screen('Flip', window);
WaitSecs(1);
%
%

numBlocks = size(BlockOrder, 2);
output = repmat(output, numBlocks, 1); %AH: second block output(2,:)

blockConditionOrder = BlockOrder(iRun, :);

for blockIndex = 1:numBlocks

% For Feedback at end
correctCountFixAll = 0;
correctCountOccAll = 0;
correctCountFixCond = zeros(1,nCond);
correctCountOccCond = zeros(1,nCond);
CompletedFixVideosPerCond = zeros(1,nCond);  %if stopped inbetween
CompletedOccVideosPerCond = zeros(1,nCond);  %if stopped inbetween

blockCondition = blockConditionOrder(blockIndex);

trialOrder = TrialOrder(blockIndex, :, iRun);
condMatrixShuffled = condMatrix(:, trialOrder); %shuffel the conditions of matrix


%% Create the ITI for between Trials, also the trigger  values between 0 and 255 for MEG and eyetracker

%AH: so here they creat ITI for the 16 trial 1.8 secs plus rand, and then
%assign it according to the trial order

ITI = 0.5 + .04 * rand(1, length(condMatrixShuffled));%Random ITI  %AH: length(condMatrixShuffled) instead of nSeq --> make ITI have all possible confMatrix indeces even if in practice, since in practive nSeq is made smaller above
% Video Number = 1:16
% 17 --> end of video
% +20 = Fixation Cross Catch Trial Starts
% + 40 = Occlusion Catch Trial Starts
% 60 = End of Catch Trial

% uniqueCond = unique(condMatrixShuffled(1, :)); %Get the condition values
% mapping = containers.Map(uniqueCond, 1:length(uniqueCond)); % Dictionary for condition --> 1:4
% mappedValues = cellfun(@(x) mapping(x), num2cell(condMatrixShuffled(1, :))); % Have 1:4 instead of condition values

%AH moved within loop
% for iTrial = 1:nSeq
%     iSeq = condMatrixShuffled(3,iTrial);
%     output(blockIndex, iTrial).ITI = ITI(iSeq);     %AH: reverse order of indexing
% end



%AH: new block
if blockCondition == 1
    Screen('FillRect', window, Conf.DotColor1, Conf.marginRect);
    Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
    DrawFormattedText(window, Conf.block1Msg, 'center', 'center', Conf.DotColor1);
    marginColor = Conf.DotColor1;
else
    Screen('FillRect', window, Conf.DotColor2, Conf.marginRect);
    Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
    DrawFormattedText(window, Conf.block2Msg, 'center', 'center', Conf.DotColor2);
    marginColor = Conf.DotColor2;
end 

Screen('Flip', window);
WaitSecs(5);

 StartTimeX = GetSecs;

for iTrial = 1:nSeq
    %AH: first was the else clause by itself above this
    Screen('FillRect', window, marginColor, Conf.marginRect);
    Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
    if Conf.showFixCrossBetweenTrials == 1
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
    else       
        DrawFormattedText(window, 'New Trial','center', 'center', 100);
    end
    
    ITIstart = Screen('Flip', window);
    iSeq = condMatrixShuffled(3,iTrial); %iSeq is the random video number
    
    %Fill output File
    output(blockIndex, iTrial).numcatch  = length(TrialStruct(iRun, iSeq).Start);
    output(blockIndex, iTrial).catch_type = [TrialStruct(iRun, iSeq).Type];
    output(blockIndex, iTrial).startcatch = [TrialStruct(iRun, iSeq).Start];
    output(blockIndex, iTrial).endcatch = [TrialStruct(iRun, iSeq).End];
    output(blockIndex, iTrial).congruence = [TrialStruct(iRun, iSeq).Congruence];
    output(blockIndex, iTrial).trial_num = iTrial;
    output(blockIndex, iTrial).lagtime = [TrialStruct(iRun, iSeq).IncongruentOffSet];
    %AH: mod(iSeq - 1, nCond) + 1 instead of mod(iSeq - 1, 4) + 1, also
    %down
    output(blockIndex, iTrial).predictability = mod(iSeq - 1, nCond) + 1; %the condition of the video
    
    %AH:
    output(blockIndex, iTrial).block_cond = blockCondition;
    output(blockIndex, iTrial).ITI = ITI(iSeq);
    output(blockIndex, iTrial).sequence = condMatrixShuffled(2, iSeq);
    output(blockIndex, iTrial).condition = condMatrixShuffled(1, iSeq);
    output(blockIndex, iTrial).PathDuration = condMatrixShuffled(4, iSeq);
    output(blockIndex, iTrial).video_nr = condMatrixShuffled(3, iSeq);
    


    currentStim = squeeze(stimuli(iSeq,:,:));  %the video
    %AH: remove the zero elements
    currentStim = currentStim(all(currentStim~=0, 2),:);
    currentCond = mod(iSeq - 1, nCond) + 1;
    
    if Conf.trackeye == true
        % EYELINK RECORDING
        
        Eyelink('Message', 'TRIALID %d', iTrial);
        Eyelink('command', 'record_status_message "SUBJECT = %d ; RUN = %d ; BLOCK = %d ; TRIAL = %d"', iSub, iRun, blockIndex, iTrial);
        
        % As specified before, it is better not to send to many Eyelink commands to the eye-tracker in a row.  So, after the command, we wait the pre-set time
        WaitSecs(elk.wait);
        
        % Here we start recording eyelink data (left/right gaze and pupil size), preceded by a short pause
        Eyelink('Command', 'set_idle_mode');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %transfer image to host
        [width, height] = Screen('WindowSize', screenNumber);
        imgfile= 'fixation.bmp';
        transferimginfo=imfinfo(imgfile);
        
        %fprintf('img file name is %s\n',transferimginfo.Filename);
        
        % image file should be 24bit or 32bit bitmap
        % parameters of ImageTransfer:
        % imagePath, xPosition, yPosition, width, height, trackerXPosition, trackerYPosition, xferoptions
        transferStatus =  Eyelink('ImageTransfer',transferimginfo.Filename,0,0,transferimginfo.Width,transferimginfo.Height,width/2-transferimginfo.Width/2 ,height/2-transferimginfo.Height/2,1);
        if transferStatus ~= 0
            fprintf('*****Image transfer Failed*****-------\n');
        end
        
        % WaitSecs(0.1); %AH:not needed?
        
        
        %%
        WaitSecs(elk.wait);
        Eyelink('StartRecording', 1, 1, 1, 1);
        %%
        WaitSecs(elk.wait);
        
    end
    
    
    disco = 0; %for photo diode
    respTime = 0; %better to initialize now
    responseStart = 0;
    XYPosition = [];
    
    
    
    % Continue with fixation cross in between trials (ITI)

    %AH: first was the else clause by itself above this
    Screen('FillRect', window, marginColor, Conf.marginRect);
    Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
    if Conf.showFixCrossBetweenTrials == 1
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
    else       
        DrawFormattedText(window, 'New Trial','center', 'center', 100);
    end

    
    ITImiddle = Screen('Flip', window);
    WaitSecs(output(blockIndex, iTrial).ITI - (ITImiddle - ITIstart));
    % Note that use  of WaitSecs is not very accurate,
    % but for the ITI it does not matter that much. It's randomly jittered anyway.
    % What is crucial is the timing of sequence onset and each frame in the sequen ce
    
    frameOnsets = (0:(length(currentStim) + 900))/desiredFrameRate; %extra frames for the catch Trial
    
    catchframe = 0; %Counter for within a Catch Trial
    catchcount = 1; %Counter which Catch Trial we are in
    iframe = 1;
    
    framecounter = 0; %for the photodiode
    framelog = []; %For Catch Trials, have we already been inside?
    
    
    testITI2 = GetSecs;
    output(blockIndex, iTrial).testITI = testITI2-ITIstart;
     output(blockIndex, iTrial).SequenceStartTime = testITI2- StartTimeX;

    
    %for iframe = 1: length(currentStim)
    %time = GetSecs
    while iframe <= length(currentStim)
        framecounter = framecounter + 1;
        
        response = -1;
        if Conf.MEG
            Datapixx('RegWrRd');						% Commit changes to/from DP
            status = Datapixx('GetDinStatus');			% Get response logger status
            
            if status.newLogFrames > 0	&& mod(status.currentWriteFrame,2)	% We've got new data in response buffer !!! --> = 1=
                [data, time] = Datapixx('ReadDinLog');	% Read data in
                
                response = bitand(data(end), responseButtonsMask);
                if response == 1
                    %AH: make any response of buttun as 1, the congruent
                    response = 1;%2; %incongruent = red button
                elseif response == 8
                    response = 1; %congruent = blue button
                else
                    disp(response);
                end
                %responses(iframe) = response;
                
                
                if response == 1 | response == 2
                    % Eyelink triggering
                    Eyelink('Message', sprintf('Response %d', response));
                    
                    % Send RT/response trigger (response value + 200)
                    Datapixx('EnableDinDebounce');
                    Datapixx('StopDoutSchedule');
                    triggerPulse = [1 0] .* (response+200);
                    Datapixx('WriteDoutBuffer', triggerPulse);
                    Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                    Datapixx('StartDoutSchedule');
                    Datapixx('RegWr');
                    % disp(sprintf('button press, trigger %d',triggerPulse(1)));
                    if ~ ismember(triggerPulse(1),TriggerValues)
                        disp(sprintf('1XXXXXXXXX %d, vpix output: %d, mask: %d',triggerPulse(1),data(end),responseButtonsMask));
                        status;
                    end
                else
                    disp(response);
                end
            end
        end
        
        %Background
        Screen('FillRect', window, marginColor, Conf.marginRect);
        Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
        
        if Conf.MEG
            %photo diode
            disco = 1 - disco; %during duration of video let square for photodiode flicker between black and white
            %discocolor = disco*255;
            Screen('FillRect', window, disco*255, [0 0 50 50]);
        end
        
        
        if any(iframe >= TrialStruct(iRun,iSeq).Start & iframe < TrialStruct(iRun,iSeq).End) && ~ ismember(iframe, framelog) %Catch Trial
            catchframe = catchframe+1;
            
            if catchframe == 1
                s = 0; %For Reacting --> it only records the first time the button was pressed
                v = 0; %In case the trial was missed --> to record the "missed" response, but only during first time
                sp = 0; %For the spatial Deviation Task
                
                CatchStart = GetSecs;
                output(blockIndex, iTrial).CatchStartTime = [output(blockIndex, iTrial).CatchStartTime, CatchStart - StartTimeX];
                
                if Conf.MEG
                    if TrialStruct(iRun,iSeq).Type(catchcount) == 1 %fixation
                        triggerPulse = [1 0] .* 100;  % AH: trigger-100
                    elseif TrialStruct(iRun,iSeq).Type(catchcount) == 2 % occlusion
                        triggerPulse = [1 0] .* 100;  % AH: trigger-100
                    end
                    
                    %Datapixx triggering
                    Datapixx('StopDoutSchedule');
                    Datapixx('WriteDoutBuffer', triggerPulse);
                    Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                    Datapixx('StartDoutSchedule');
                    %disp(sprintf('catch trial start, trigger %d',triggerPulse(1)));
                    if ~ ismember(triggerPulse(1),TriggerValues)
                        disp(sprintf('2XXXXXXXXX %d',triggerPulse(1)))
                    end
                    % Eyelink triggering
                    Eyelink('Message', sprintf('Video onset %d', triggerPulse(1)));
                    
                    Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                    
                    % Reset and fire up the response logger
                    Datapixx('SetDinLog');
                    Datapixx('StartDinLog');
                end
            end
            
            
            if TrialStruct(iRun,iSeq).Type(catchcount) == 1 %fixation
                if catchframe < Catch.FixDuration*Conf.refrate                 %Fixation Cross Color Change
                    
                    Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    % AH: draw second dot
                    Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                    XYPosition(framecounter, :) = currentStim(iframe,:);
                    
                    
                    Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                    Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColorChange, [xCenter yCenter], 2);
                    Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    
                    
                    if Conf.MEG && response>0 && s == 0 %For MEG, the response is recorded at the beginning of the loop
                        s = 1;
                        
                        output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, response];
                        
                        respTime = time(end);	%
                        
                        %in output
                        output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                        
                        % Eyelink triggering
                        Eyelink('Message', sprintf('Response %d', response));
                        
                        % Send RT/response trigger (response value + 128)
                        %Datapixx('EnableDinDebounce');
                        %Datapixx('StopDoutSchedule');
                        
                        %triggerPulse = [1 0] .* 101; %correct
                        %Datapixx('WriteDoutBuffer', triggerPulse);
                        
                        %% ???
                        %Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                        %Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                        
                        %Datapixx('StartDoutSchedule');
                        %Datapixx('RegWr');
                        %disp(sprintf('fixation response, trigger (correctness) %d, RT: %0.1f',triggerPulse(1),respTime - responseStart));
                        
                    end
                    
                    if Conf.MEG == false
                        
                        [keyIsDown,respTime,keyCode] = KbCheck; %not  outside loop
                        
                        if keyIsDown == true && keyCode(Key.spaceKey) == 0 &&  s == 0  %only first time when key is pressed, enter here %if pressed, record RT
                            
                            output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                            s = 1;
                            
                            if keyCode(Key.LeftKey) == true
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 1];
                                response = 1;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                                
                            elseif keyCode(Key.RightKey) == true
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 2];
                                response = 2;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                            end
                        end
                    end
                    
                    
                    %Response Period of Fixation Cross
                elseif (Catch.FixDuration*Conf.refrate) <= catchframe && catchframe < ((Catch.FixDuration + Catch.FixResponseDuration)*Conf.refrate)
                    Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                    Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                    Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    
                    Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    % AH: draw second dot
                    Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                    XYPosition(framecounter, :) = currentStim(iframe,:);
                    
                    
                    
                    if Conf.MEG && response>0 && s == 0
                        s = 1;
                        
                        output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, response];
                        respTime = time(end);	%
                        
                        %in output
                        
                        output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                        
                        % Eyelink triggering
                        Eyelink('Message', sprintf('Response %d', response));
                        
                        % Send RT/response trigger (response value + 128)
                        %Datapixx('EnableDinDebounce');
                        %Datapixx('StopDoutSchedule');
                        
                        %triggerPulse = [1 0] .* 101; %correct
                        %Datapixx('WriteDoutBuffer', triggerPulse);
                        
                        %% ??
                        %Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                        % Datapixx('StartDoutSchedule');
                        % Datapixx('RegWr');
                        %disp(sprintf('fixation response, trigger (correctness) %d, RT: %0.1f',triggerPulse(1), respTime - responseStart));
                        
                    end
                    
                    if Conf.MEG == false
                        
                        [keyIsDown,respTime,keyCode] = KbCheck;
                        
                        if keyIsDown == true && keyCode(Key.spaceKey) == 0 &&  s == 0  %only first time when key is pressed, enter here %if pressed, record RT
                            
                            output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                            s = 1;
                            
                            if keyCode(Key.LeftKey ) == true
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 1];
                                response = 1;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                            elseif keyCode(Key.RightKey) == true
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 2];
                                %return
                                response = 2;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                            end
                        end
                    end
                    
                    % Waiting Period of Fixation Cross
                elseif ((Catch.FixDuration + Catch.FixResponseDuration)*Conf.refrate) <= catchframe && catchframe < ((Catch.FixDuration + Catch.FixResponseDuration + Catch.FixWaitingDuration)*Conf.refrate)
                    %Fixation Cross
                    % Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                    % Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                    % Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    %Screen('DrawDots', window, currentStim(iframe,:), Conf.sizeDotInPixel, Conf.DotColor, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    
                    XYPosition(framecounter, :) = NaN;
                    
                    
                    if s == 0 %if no response has been given
                        % Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        % Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.colors(2,:), [xCenter yCenter], 2);
                        % Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                        %
                        % DrawFormattedText(window,  'Missed',  'center', yCenter-Conf.rectSize(2)*(1/10), Conf.white);
                        
                        % AH: since when size is changed within the qst
                        % mark thing it remains and we want missed to be
                        % as text size
                        Screen('TextSize', window, Conf.Textsize);
                        Screen('TextFont',window, Conf.Text);

                        DrawFormattedText(window,  'Missed',  'center', 'center', Conf.colors(2,:));
                        
                        if v == 0 %the first frame of missed response, send a trigger
                            output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 0];
                            output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT, NaN];
                            v = 1;
                            
                            if Conf.MEG
                                % Send RT/response trigger (response value + 128)
                                Datapixx('EnableDinDebounce');
                                Datapixx('StopDoutSchedule');
                                triggerPulse = [1 0] .* 113;
                                Datapixx('WriteDoutBuffer', triggerPulse);
                                Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                                Datapixx('StartDoutSchedule');
                                %disp(sprintf('missed, trigger %d',triggerPulse(1)));
                                if ~ ismember(triggerPulse(1),TriggerValues)
                                    disp(sprintf('3XXXXXXXXX %d',triggerPulse(1)))
                                end
                                Datapixx('RegWr');
                            end
                        end
                    else
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    end
                    
                elseif  catchframe == ((Catch.FixDuration + Catch.FixResponseDuration  + Catch.FixWaitingDuration)*Conf.refrate) %last Catch Trial frame
                    
                    XYPosition(framecounter, :) = NaN;
                    
                    if s == 0
                        % AH: since when size is changed within the qst
                        % mark thing it remains and we want missed to be
                        % as text size
                        Screen('TextSize', window, Conf.Textsize);
                        Screen('TextFont',window, Conf.Text);

                        DrawFormattedText(window,  'Missed',  'center', 'center', Conf.colors(2,:));
                    else
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    end
                    
                    
                    
                    framelog = [framelog, TrialStruct(iRun,iSeq).Start(catchcount):(TrialStruct(iRun,iSeq).End(catchcount) - 1)];
                    iframe = TrialStruct(iRun,iSeq).Start(catchcount)-round(Conf.resetTime*Conf.refrate);
                    
                    catchcount = catchcount+1; %Next Catch Trial
                    catchframe = 0; %Reset the Catchframe Counter
                    
                    CatchEnd = GetSecs;
                    output(blockIndex, iTrial).CatchEndTime = [output(blockIndex, iTrial).CatchEndTime, CatchEnd - StartTimeX];
             
                    
                    if Conf.MEG
                        triggerPulse = [1 0] .* 101;  % AH: trigger-101
                        
                        %Datapixx triggering
                        Datapixx('StopDoutSchedule');
                        Datapixx('WriteDoutBuffer', triggerPulse);
                        Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                        Datapixx('StartDoutSchedule');
                        %disp(sprintf('catch trial end, trigger %d',triggerPulse(1)));
                        if ~ ismember(triggerPulse(1),TriggerValues)
                            disp(sprintf('4XXXXXXXXX %d',triggerPulse(1)))
                        end
                        % Eyelink triggering
                        Eyelink('Message', sprintf('Video onset %d', triggerPulse(1)));
                        
                        Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                        
                        % Reset and fire up the response logger
                        Datapixx('SetDinLog');
                        Datapixx('StartDinLog');
                    end
                end
                
                
            elseif TrialStruct(iRun,iSeq).Type(catchcount) == 2 %oclusion
                
                if catchframe < Conf.refrate*Catch.OcclDuration %Dot disappears

                    Screen('FillRect', window, marginColor, Conf.marginRect);
                    Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
                    
                    %Draw dot in red during occlusion:
                    %Screen('DrawDots', window, currentStim(iframe,:), Conf.sizeDotInPixel, [200, 0, 0 ], [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    
                    if Conf.showQstMarkOnOcclusion == 0
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    else
                        Screen('TextSize', window, Conf.Textsize*Conf.qstMarkSizeFactor);
                        Screen('TextFont',window, Conf.Text);
                        DrawFormattedText(window,  '?',  'center', 'center', Conf.normalQstMarkColor);
                    end

                    sp = sp+1; %Counter for the simulated Dotpaths with different position
                    
                    spacial_NewPos(1) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 1) *Conf.ppd  + Conf.rectCoords(1);
                    spacial_NewPos(2) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 2) *Conf.ppd + Conf.rectCoords(2);
                    %AH: second point
                    spacial_NewPos(3) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 3) *Conf.ppd  + Conf.rectCoords(1);
                    spacial_NewPos(4) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 4) *Conf.ppd + Conf.rectCoords(2);

                    if blockCondition == 1
                        Screen('DrawDots', window,  spacial_NewPos(1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                        
                        % AH: draw dot 2 normally not jitery
                        Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                            
                    % AH: draw second dot
                    elseif blockCondition == 2
                        Screen('DrawDots', window, spacial_NewPos(3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);

                        % AH: draw dot 1 normally not jitery
                        Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1);
                    end

                    XYPosition(framecounter, :) = spacial_NewPos;

                    % XYPosition(framecounter, :) = NaN;

                    %AH: here is when the dot first disappears, so we
                    %prompt a question mark
                    
                    
                    
                    %Video
                elseif  Conf.refrate*Catch.OcclDuration <= catchframe && catchframe < (Conf.refrate*(Catch.OcclDuration + Catch.OcclVideoDuration ))
                    if catchframe == Conf.refrate*Catch.OcclDuration
                        %responseStart = GetSecs;
                        FrozenFrame = iframe;
                    end
                    
                    if Conf.showQstMarkOnOcclusion == 0
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )          
                    else
                        Screen('TextSize', window, Conf.Textsize*Conf.qstMarkSizeFactor);
                        Screen('TextFont',window, Conf.Text);
                        DrawFormattedText(window,  '?',  'center', 'center', Conf.normalQstMarkColor);
                    end

                    %AH: here is when the dot comes back after occlusion, should we
                    %prompt a question mark also??

                    
                    if TrialStruct(iRun, iSeq).Congruence(catchcount) == 1 %Congruent
                        if Conf.Task == 2
                            Screen('DrawDots', window, currentStim(FrozenFrame,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                            % AH: draw second dot
                            Screen('DrawDots', window, currentStim(FrozenFrame,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                        else %if dot position is not frozen
                            Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                            % AH: draw second dot
                            Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                        end
                        
                        XYPosition(framecounter, :) = currentStim(iframe,:);
                        
                        
                    elseif TrialStruct(iRun, iSeq).Congruence(catchcount) == 2 %Incongruent
                        %For Debugging: Draw where dot reappers
                        % Screen('DrawDots', window, currentStim(FrozenFrame,:), Conf.sizeDotInPixel, [0, 0, 200], [0 0], 1); % the last two variables are the 'center' and 'dot form'
                        
                        if Conf.Task == 1 %Spatial Deviation
                            sp = sp+1; %Counter for the simulated Dotpaths with different position
                            
                            spacial_NewPos(1) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 1) *Conf.ppd  + Conf.rectCoords(1);
                            spacial_NewPos(2) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 2) *Conf.ppd + Conf.rectCoords(2);
                            %AH: second point
                            spacial_NewPos(3) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 3) *Conf.ppd  + Conf.rectCoords(1);
                            spacial_NewPos(4) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(sp, 4) *Conf.ppd + Conf.rectCoords(2);

                            Screen('DrawDots', window,  spacial_NewPos(1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                            % AH: draw second dot
                            Screen('DrawDots', window, spacial_NewPos(3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                            XYPosition(framecounter, :) = spacial_NewPos;
                            
                            
                        elseif Conf.Task == 0 %Temporal deviation
                            randomstart = TrialStruct(iRun, iSeq).IncongruentOffSet(catchcount);
                            random_offset = (Catch.OcclDeviance*Conf.refrate)*randomstart;
                            Screen('DrawDots', window, currentStim(iframe+random_offset,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                            % AH: draw second dot
                            Screen('DrawDots', window, currentStim(iframe+random_offset,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                            XYPosition(framecounter, :) = currentStim(iframe+random_offset,:);
                            
                            
                        elseif Conf.Task == 2 %Task where Dot is frozen
                            spacial_NewPos(1) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(1, 1) *Conf.ppd + Conf.rectCoords(1);
                            spacial_NewPos(2) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(1, 2) *Conf.ppd + Conf.rectCoords(2);
                            %AH: second point
                            spacial_NewPos(3) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(1, 3) *Conf.ppd + Conf.rectCoords(1);
                            spacial_NewPos(4) = TrialStruct(iRun, iSeq).SimulatedPath{catchcount}(1, 4) *Conf.ppd + Conf.rectCoords(2);

                            Screen('DrawDots', window,  spacial_NewPos(1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                            % AH: draw second dot
                            Screen('DrawDots', window, spacial_NewPos(3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                            XYPosition(framecounter, :) = currentStim(iframe,:);
                            
                            
                        elseif Conf.Task == 3 %Task with new angle
                            sp = sp+1;
                            spacial_NewPos(1) = TrialStruct(iRun, iSeq).SimulatedNewAngle{catchcount}(sp, 1)*Conf.ppd + Conf.rectCoords(1);
                            spacial_NewPos(2) = TrialStruct(iRun, iSeq).SimulatedNewAngle{catchcount}(sp, 2)*Conf.ppd + Conf.rectCoords(2);
                            %AH: second point
                            spacial_NewPos(3) = TrialStruct(iRun, iSeq).SimulatedNewAngle{catchcount}(sp, 3)*Conf.ppd + Conf.rectCoords(1);
                            spacial_NewPos(4) = TrialStruct(iRun, iSeq).SimulatedNewAngle{catchcount}(sp, 4)*Conf.ppd + Conf.rectCoords(2);

                            Screen('DrawDots', window,  spacial_NewPos(1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                            % AH: draw second dot
                            Screen('DrawDots', window, spacial_NewPos(3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                            XYPosition(framecounter, :) = spacial_NewPos;
                            
                        end
                    end
                    
                    if Conf.MEG && response>0 && s == 0
                        s = 1;
                        
                        output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, response];
                        
                        respTime = time(end);	%
                        
                        %in output
                        
                        output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                        
                        %                         if output(blockIndex, iTrial).response(catchcount) == output(blockIndex, iTrial).congruence(catchcount)
                        %                             triggerPulse = [1 0] .* 101; % correct response
                        %                         else
                        %                             triggerPulse = [1 0] .* 99; %false response
                        %                         end
                        
                        %% ??
                        % Send RT/response trigger (response value + 128)
                        %Datapixx('EnableDinDebounce');
                        %Datapixx('StopDoutSchedule');
                        
                        %Datapixx('WriteDoutBuffer', triggerPulse);
                        
                        % Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                        %Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                        
                        %Datapixx('StartDoutSchedule');
                        %Datapixx('RegWr');
                        %disp(sprintf('Occlusion, trigger (correctness) %d, RT: %0.1f',triggerPulse(1), respTime - responseStart));
                        
                    end
                    
                    if Conf.MEG == false
                        
                        [keyIsDown,respTime,keyCode] = KbCheck;
                        
                        if keyIsDown == true && keyCode(Key.spaceKey) == 0 &&  s == 0 %only first time when key is pressed, enter here %if pressed, record RT
                            
                            output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                            s = 1;
                            
                            if keyCode(Key.LeftKey )
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 1];
                                response = 1;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                            elseif keyCode(Key.RightKey)
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 2];
                                response = 2;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                                %return
                            end
                        end
                    end
                    
                    %Response
                elseif (Conf.refrate*(Catch.OcclDuration + Catch.OcclVideoDuration )) <= catchframe && catchframe < (Conf.refrate*(Catch.OcclDuration + Catch.OcclVideoDuration + Catch.OcclResponseDuration ))
                    
                    if Conf.showQstMarkOnOcclusion == 0 %AH: this is marisas case
                        % DrawFormattedText(window,  '?',  'center', 'center', Conf.normalQstMarkColor);
                    else
                        Screen('TextSize', window, Conf.Textsize*Conf.qstMarkSizeFactor);
                        Screen('TextFont',window, Conf.Text);
                        DrawFormattedText(window,  '?',  'center', 'center', Conf.warningQstMarkColor);
                    end

                    if Conf.jitterCatches
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                    
                        Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                        % AH: draw second dot
                        Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                        XYPosition(framecounter, :) = currentStim(iframe,:);
                    else
                    
                        XYPosition(framecounter, :) = NaN;

                    end
                    
                    if Conf.MEG && response>0 && s == 0
                        s = 1;
                        
                        output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, response];
                        respTime = time(end);	%
                        
                        %in output
                        
                        output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                        
                        %                         if output(blockIndex, iTrial).response(catchcount) == output(blockIndex, iTrial).congruence(catchcount)
                        %                             triggerPulse = [1 0] .* 101; % correct response
                        %                         else
                        %                             triggerPulse = [1 0] .* 99; %false response
                        %                         end
                        
                        
                        % Send RT/response trigger (response value + 128)
                        %                         Datapixx('EnableDinDebounce');
                        %                         Datapixx('StopDoutSchedule');
                        %
                        %                         Datapixx('WriteDoutBuffer', triggerPulse);
                        %
                        %                         Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                        %                         %Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                        %
                        %                         Datapixx('StartDoutSchedule');
                        %                         Datapixx('RegWr');
                        %                         disp(sprintf('Occlusion, trigger (correctness) %d, RT: %0.1f',triggerPulse(1), respTime - responseStart));
                        
                    end
                    
                    if Conf.MEG == false
                        
                         [keyIsDown,respTime,keyCode] = KbCheck;
                        
                        if keyIsDown == true && keyCode(Key.spaceKey) == 0 &&  s == 0 %only first time when key is pressed, enter here %if pressed, record RT
                            
                            output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT,  respTime - responseStart];
                            s = 1;
                            
                            if keyCode(Key.LeftKey) == true
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 1];
                                response = 1;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                            elseif keyCode(Key.RightKey) == true
                                output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 2];
                                response = 2;
                                if Conf.trackeye
                                    Eyelink('Message', sprintf('Response %d', response));
                                end
                                %return
                            else
                                s = 0;
                            end
                        end
                    end
                    
                    
                    %Waiting
                elseif (Conf.refrate*(Catch.OcclDuration + Catch.OcclVideoDuration + Catch.OcclResponseDuration )) <= catchframe && catchframe < (Conf.refrate*(Catch.OcclDuration + Catch.OcclVideoDuration + Catch.OcclResponseDuration + Catch.OcclWaitingDuration ))
                    %DrawFormattedText(window,  'Waiting',  'center', yCenter/3, Conf.white);
                    
                    % XYPosition(framecounter, :) = NaN;
                    
                    
                    if s == 1 %if button was pressed
                        
                        if Conf.practice == true  %Feedback given in practice trials
                            if output(blockIndex, iTrial).response(catchcount) == output(blockIndex, iTrial).congruence(catchcount)
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                                Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix,  Conf.colors(1,:), [xCenter yCenter], 2);
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                                
                            else
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                                Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix,  Conf.colors(2,:), [xCenter yCenter], 2);
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                                
                            end
                            
                        else %if not practice trial, show no Feedback
                            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                            Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                        end
                        
                        
                    elseif s == 0
                        
                        % % AH: since when size is changed within the qst
                        % % mark thing it remains and we want missed to be
                        % % as text size
                        % Screen('TextSize', window, Conf.Textsize);
                        % Screen('TextFont',window, Conf.Text);
                        % 
                        % DrawFormattedText(window,  'Missed',  'center', 'center', Conf.colors(2,:));
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )

                        if v == 0
                            output(blockIndex, iTrial).response = [output(blockIndex, iTrial).response, 0];
                            output(blockIndex, iTrial).RT = [output(blockIndex, iTrial).RT, NaN];
                            
                            v = 1;
                            
                            if Conf.MEG
                                % % Send RT/response trigger (response value + 128)
                                % Datapixx('EnableDinDebounce');
                                % Datapixx('StopDoutSchedule');
                                % triggerPulse = [1 0] .* 113;
                                % Datapixx('WriteDoutBuffer', triggerPulse);
                                % Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                                % %Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                                % 
                                % Datapixx('StartDoutSchedule');
                                % %disp(sprintf('missed, trigger %d',triggerPulse(1)));
                                % if ~ ismember(triggerPulse(1),TriggerValues)
                                %     disp(sprintf('5XXXXXXXXX %d',triggerPulse(1)))
                                % end
                                % Datapixx('RegWr');
                                
                            end
                        end
                    end

                    Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    % AH: draw second dot
                    Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                    XYPosition(framecounter, :) = currentStim(iframe,:);
                    
                elseif catchframe ==  Conf.refrate*(Catch.OcclDuration + Catch.OcclVideoDuration + Catch.OcclResponseDuration + Catch.OcclWaitingDuration ) %last Catch Trial frame
                    
                    % XYPosition(framecounter, :) = NaN;
                    
                    if s == 1 %if button was pressed
                        if Conf.practice == true  %Feedback given in practice trials
                            if output(blockIndex, iTrial).response(catchcount) == output(blockIndex, iTrial).congruence(catchcount)
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                                Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix,  Conf.colors(1,:), [xCenter yCenter], 2);
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                                
                            else
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                                Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix,  Conf.colors(2,:), [xCenter yCenter], 2);
                                Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                            end
                            
                        else %if not practice trial, show no Feedback
                            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                            Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
                        end
                        
                        
                    elseif s == 0
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
                        Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
                        Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )

                        % DrawFormattedText(window,  'Missed',  'center', yCenter-Conf.rectSize(2)*(1/10), Conf.white);
                        
                        % % AH: since when size is changed within the qst
                        % % mark thing it remains and we want missed to be
                        % % as text size
                        % Screen('TextSize', window, Conf.Textsize);
                        % Screen('TextFont',window, Conf.Text);
                        % 
                        % DrawFormattedText(window,  'Missed',  'center', 'center', Conf.colors(2,:));
                        
                    end

                    Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
                    % AH: draw second dot
                    Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
                    XYPosition(framecounter, :) = currentStim(iframe,:);

                    
                    framelog = [framelog, TrialStruct(iRun,iSeq).Start(catchcount):(TrialStruct(iRun,iSeq).End(catchcount)- 1) ];
                    
                    %AH: no change in iFrame needed
                    % iframe = TrialStruct(iRun,iSeq).Start(catchcount) - round(Conf.resetTime*Conf.refrate);
                    
                    
                    catchcount = catchcount+1;
                    catchframe = 0;
                    
                    CatchEnd = GetSecs;
                    output(blockIndex, iTrial).CatchEndTime = [output(blockIndex, iTrial).CatchEndTime, CatchEnd - StartTimeX];
                    
                    if Conf.MEG
                        triggerPulse = [1 0] .* 101;  % AH: trigger-101
                        %Datapixx triggering
                        Datapixx('StopDoutSchedule');
                        Datapixx('WriteDoutBuffer', triggerPulse);
                        Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                        Datapixx('StartDoutSchedule');
                        %disp(sprintf('last catch trial frame, trigger %d',triggerPulse(1)));
                        if ~ ismember(triggerPulse(1),TriggerValues)
                            disp(sprintf('6XXXXXXXXX %d',triggerPulse(1)))
                        end
                        % Eyelink triggering
                        Eyelink('Message', sprintf('Video onset %d', triggerPulse(1)));
                        
                        Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                        
                        % Reset and fire up the response logger
                        Datapixx('SetDinLog');
                        Datapixx('StartDinLog');
                    end
                end
            end
        else %If not a catch trial
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusOut, yCenter - Conf.RadiusOut, xCenter+Conf.RadiusOut, yCenter+ Conf.RadiusOut]) %Conf.RadiusOut*2
            Screen('DrawLines', window, Conf.allCoords, Conf.lineWidthPix, Conf.CrossColor, [xCenter yCenter], 2);
            Screen('FillOval', window, Conf.ColorOval, [xCenter - Conf.RadiusIn, yCenter - Conf.RadiusIn, xCenter+Conf.RadiusIn, yCenter+ Conf.RadiusIn], Conf.RadiusIn*2 )
            
            Screen('DrawDots', window, currentStim(iframe,1:2), Conf.sizeDotInPixel, Conf.DotColor1, [0 0], 1); % the last two variables are the 'center' and 'dot form'
            % AH: draw second dot
            Screen('DrawDots', window, currentStim(iframe,3:4), Conf.sizeDotInPixel, Conf.DotColor2, [0 0], 1);
            XYPosition(framecounter, :) = currentStim(iframe,:);
            
            
            if Conf.MEG %Triggering
                if iframe == 1 %First Frame

                    % should be a value between 1 and 80: 1-40 = attend
                    % red, 41-80: attend green
                    stimulus =  condMatrixShuffled(2, iTrial);
                    if condMatrixShuffled(1, iTrial)==40
                        stimulus = stimulus+20;
                    end
                    if blockCondition==2
                        stimulus = stimulus+40;
                    end
                    
                    if output(blockIndex, iTrial).numcatch         
                        stimulus = 102; % AH: trigger-102 catch trial
                    end
          
                    triggerPulse = [1 0] .* stimulus; % AH: trigger-stimulus

                    %Datapixx triggering
                    Datapixx('StopDoutSchedule');
                    Datapixx('WriteDoutBuffer', triggerPulse);
                    Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx) %$$ change trigger to video_nr
                    Datapixx('StartDoutSchedule');
                    
                    % Eyelink triggering
                    Eyelink('Message', sprintf('Video onset %d', triggerPulse(1)));
                    
                    % White square for photodiode
                    %Screen('FillRect', window, [1 1 1], [0 0 30 30]);
                    
                    Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                    
                    % Reset and fire up the response logger
                    Datapixx('SetDinLog');
                    Datapixx('StartDinLog');
                    %disp(sprintf('first frame, trigger %d',triggerPulse(1)));
                    if ~ ismember(triggerPulse(1),TriggerValues)
                        disp(sprintf('7XXXXXXXXX %d',triggerPulse(1)))
                    end
                elseif iframe == length(currentStim) %last frame of video
                    if output(blockIndex, iTrial).numcatch         
                        lastFrameTrigger = 103; % AH: trigger-103 catch trial
                    else
                        lastFrameTrigger = 81; %AH: trigger-81 as last trigger
                    end

                    triggerPulse = [1 0] .* lastFrameTrigger; %AH: trigger-lastFrameTrigger as last trigger
                    %Datapixx triggering
                    Datapixx('StopDoutSchedule');
                    Datapixx('WriteDoutBuffer', triggerPulse);
                    Datapixx('SetDoutSchedule', 1.0/Conf.refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                    Datapixx('StartDoutSchedule');
                    
                    % Eyelink triggering
                    Eyelink('Message', sprintf('Video onset %d', triggerPulse(1)));
                    
                    Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                    
                    % Reset and fire up the response logger
                    Datapixx('SetDinLog');
                    Datapixx('StartDinLog');
                    %disp(sprintf('last frame, trigger %d',triggerPulse(1)));
                    if ~ ismember(triggerPulse(1),TriggerValues)
                        disp(sprintf('8XXXXXXXXX %d',triggerPulse(1)))
                    end
                end
            end
        end
        
        output(blockIndex, iTrial).XYPositionPerFrame(framecounter, :) = XYPosition(framecounter, :);
        
        
        timePassed = GetSecs - testITI2;
        while timePassed <= frameOnsets(framecounter)
            timePassed = GetSecs - testITI2;
            
             if Conf.MEG
                Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                %Datapixx('SetDinLog');
                %Datapixx('StartDinLog');
                
                Datapixx('RegWrRd');						% Commit changes to/from DP
                status = Datapixx('GetDinStatus');			% Get response logger status
                
            else
                [keyIsDown,respTime,keyCode] = KbCheck; % keep collecting responses
            end
            
        end
        
        %disp(sprintf('trial: %d; frame: %d; current: %0.5f s; desired onset: %0.5f s',iTrial, iframe, timePassed, frameOnsets(iframe)));
        
                
        %AH: is you have Nan in correct_response meand it is a catch and a
        %no response so missed, if no cath it is empty
        if any(output(blockIndex, iTrial).response == 0) && iframe == size(currentStim,1)

            % AH: since when size is changed within the qst
            % mark thing it remains and we want missed to be
            % as text size
            Screen('FillRect', window, marginColor, Conf.marginRect);
            Screen('FillRect', window, Conf.RectColor, Conf.rectCoords);
            Screen('TextSize', window, Conf.Textsize);
            Screen('TextFont',window, Conf.Text);

            DrawFormattedText(window,  'Missed',  'center', 'center', Conf.colors(2,:));
                    
            %Now send all instructions to the Datapixx box, as close as possible to the actual screen flip in PTB
        
            if Conf.MEG
                Datapixx('RegWrVideoSync');
            end
            
            vbl = Screen('Flip', window);
            flipTime(framecounter) = vbl;

            WaitSecs(0.3);
 
        else
                    
            %Now send all instructions to the Datapixx box, as close as possible to the actual screen flip in PTB
        
            if Conf.MEG
                Datapixx('RegWrVideoSync');
            end
            
            vbl = Screen('Flip', window);
            flipTime(framecounter) = vbl;
        end


        
        if framecounter > 1
            output(blockIndex, iTrial).FrameDuration = [output(blockIndex, iTrial).FrameDuration,  diff([flipTime(framecounter - 1) flipTime(framecounter)])*1000]; %add the frame time in ms
        else
            output(blockIndex, iTrial).FrameDuration = [output(blockIndex, iTrial).FrameDuration, diff([testITI2 flipTime(framecounter)])*1000]; %add the frame time in ms
        end
        
        
        if catchframe == Conf.refrate*Catch.OcclDuration && TrialStruct(iRun,iSeq).Type(catchcount) == 2
            if Conf.MEG
                Datapixx('SetMarker');
                Datapixx('RegWrRd');
                responseStart = Datapixx('GetMarker');
            else
                responseStart = GetSecs;
            end
        elseif  catchframe == 1 && TrialStruct(iRun,iSeq).Type(catchcount) == 1
            if Conf.MEG
                Datapixx('SetMarker');
                Datapixx('RegWrRd');
                responseStart = Datapixx('GetMarker');
            else
                responseStart = GetSecs;
            end
        end
        
        [keyIsDown,respTime,keyCode] = KbCheck; %to stop experiment
        if keyCode(Key.escapeKey)
            sca;
            disp('The escape key has been pressed');
            return
        end
        
        %% Test Timing (Ask Ingmar why)
        %temporary timing test
        if iframe == 1
            test1 = GetSecs;
        end
        if iframe == size(currentStim,1) % ceil(30*Conf.refrate) %30s long)
            test2 = GetSecs;
            output(blockIndex, iTrial).testTiming = test2-test1;
        end
        
        iframe = iframe + 1; %next frame
    end  % End of Frames
    
    
    % Catch Trials: Correct Response
    %Congruence: 0 = Fixation,  1= congruent, 2 = incongruent
    %Response: 0 = missed, 1 = leftkey; 2 = rightkey
    %Maybe put it outside of for loop in separate for loop?
    
    if Conf.trackeye
        % Stop eyelink
        Eyelink('StopRecording');
    end
    
    for i = 1:length(output(blockIndex, iTrial).response)
        response = output(blockIndex, iTrial).response(i);
        congruence = output(blockIndex, iTrial).congruence(i);
        catchtype = output(blockIndex, iTrial).catch_type(i);
        
        if response == 1
            if congruence == 1
                output(blockIndex, iTrial).correct_response(i) = 1;
                correctCountOccAll =  correctCountOccAll + 1;
                correctCountOccCond(currentCond) = correctCountOccCond(currentCond) + 1;
                
            elseif congruence == 0
                output(blockIndex, iTrial).correct_response(i) = 1;
                correctCountFixAll = correctCountFixAll + 1;
                correctCountFixCond(currentCond) = correctCountFixCond(currentCond) + 1;
                
            else
                output(blockIndex, iTrial).correct_response(i) = 0;
            end
            
        elseif response == 2
            if congruence == 2
                output(blockIndex, iTrial).correct_response(i) = 1;
                correctCountOccAll =  correctCountOccAll + 1;
                correctCountOccCond(currentCond) = correctCountOccCond(currentCond) + 1;
                
            elseif congruence == 0
                output(blockIndex, iTrial).correct_response(i) = 1;
                correctCountFixAll = correctCountFixAll + 1;
                correctCountFixCond(currentCond) = correctCountFixCond(currentCond) + 1;
                
            else
                output(blockIndex, iTrial).correct_response(i) = 0;
            end
            
        elseif response == 0
             output(blockIndex, iTrial).correct_response(i) = NaN;
            
        else
            output(blockIndex, iTrial).correct_response(i) = NaN;
            
        end

        if catchtype == 1
            CompletedFixVideosPerCond(currentCond) = CompletedFixVideosPerCond(currentCond) + 1;
        elseif catchtype == 2
            CompletedOccVideosPerCond(currentCond) =  CompletedOccVideosPerCond(currentCond) + 1;
        end
    end
    
EndTime = GetSecs;
output(blockIndex, iTrial).SequenceEndTime = EndTime- StartTimeX;
    
end %end of iTrial


% FEEDBACK AFTER EACH BLOCK
FixPercentCorrectAll = (correctCountFixAll / sum([output(blockIndex,:).catch_type] == 1)) * 100;
OccPercentCorrectAll = (correctCountOccAll / sum([output(blockIndex,:).catch_type] == 2)) * 100;


disp(sprintf('Fixation Cross Correct Total: %0.3f',FixPercentCorrectAll));
disp(sprintf('Occlusion Correct Total: %0.3f',OccPercentCorrectAll));

fprintf('\n');
fprintf('\n');

disp('Fixation Cross per Condition:')
for i = 1:nCond
    CorrectFix = (correctCountFixCond(i) / CompletedFixVideosPerCond(i)) * 100;
    if isnan(CorrectFix)
        CorrectFix = 0;
    end
    fprintf('Condition %d Correct: %0.3f\n', i, CorrectFix);
end

fprintf('\n');
fprintf('\n');
disp('Occlusion per Condition:')
for i = 1:nCond
    CorrectFix = (correctCountOccCond(i) / CompletedOccVideosPerCond(i)) * 100;
    if isnan(CorrectFix)
        CorrectFix = 0;
    end
    fprintf('Condition %d Correct: %0.3f\n', i, CorrectFix);
end

end %block loop

% DrawFormattedText(window,['You had ' num2str(FixPercentCorrect) '% correct in the Fixation Task \n\n and '...
%     num2str(OccPercentCorrect) '% correct in the Occlusion Task'],'center', yCenter/3, Conf.black);
DrawFormattedText(window,'Please tell the experimenter that you finished this block \n\n and take a small break','center', 'center', Conf.black);
Screen('Flip', window);
% sound(sin(2 * pi * 500 * linspace(0, 0.5, 0.5 * 22050)), 22050); %create a beeping sound or half a second length

KbStrokeWait;



if Conf.trackeye
    if Conf.practice == 0
        try
            fprintf('Receiving data file ''%s''\n', elk.edfFile );
            status = Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2 == exist(elk.edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', elk.edfFile, savedir );
            end
        catch
            fprintf('Problem receiving data file ''%s''\n', elk.edfFile );
        end
    end
    
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    Eyelink('ShutDown');
end

if Conf.MEG
    % Close DataPixx
    Datapixx('StopAllSchedules');		% Stop all schedules
    Datapixx('SetDoutValues', 0);       % Reset triggers to zero
    Datapixx('StopDinLog');             % Stop response buttons recording
    Datapixx('Close');                  % Close DataPixx
end



% Psychtoolbox
Priority(0);
sca



%% New Output of CatchTrials

CatchOutput = struct(...
    'Subject_num', [], ...
    'iRun', [], ...
    'Condition', [], ...
    'PathDuration', [], ...
    'Start', [], ...
    'End', [], ...
    'Type', [], ...
    'Lagtime', [],...
    'Congruence', [], ...
    'Correctness', [], ...
    'Reaction_Time', [], ...
    'VideoNr', [], ...
    'Response', []);


counter = 1;

for iBlock = 1:numBlocks    
for iSeq = 1:nSeq
    for iCatch = 1:output(iBlock, iSeq).numcatch
        CatchOutput(counter).Subject_num = output(iBlock, iSeq).subject_num;
        CatchOutput(counter).iRun = output(iBlock, iSeq).run_num;%AH: was block_num
        CatchOutput(counter).block_cond = output(iBlock, iSeq).block_cond;
        CatchOutput(counter).Condition = output(iBlock, iSeq).condition;
        CatchOutput(counter).PathDuration = output(iBlock, iSeq).PathDuration;
        CatchOutput(counter).Task = Conf.Task;
        CatchOutput(counter).Start =  output(iBlock, iSeq).startcatch(iCatch);
        CatchOutput(counter).End = output(iBlock, iSeq).endcatch(iCatch);
        CatchOutput(counter).Type = output(iBlock, iSeq).catch_type(iCatch);
        CatchOutput(counter).Lagtime = output(iBlock, iSeq).lagtime(iCatch);
        CatchOutput(counter).Congruence = output(iBlock, iSeq).congruence(iCatch);
        CatchOutput(counter).Correctness = output(iBlock, iSeq).correct_response(iCatch);
        CatchOutput(counter).Reaction_Time = output(iBlock, iSeq).RT(iCatch);
        CatchOutput(counter).VideoNr = output(iBlock, iSeq).video_nr;
        CatchOutput(counter).Response = output(iBlock, iSeq).response;
        counter = counter +1; %For Counting
    end
end
end



%% Analysis


%t = date('now', 'Format', 'yyyy-MM-dd''_THHmmss');
fileName = sprintf('%s_SUB%02d_RUN%02d',Conf.experiment_name,iSub,iRun); % sprintf('%s_OutputFiles.mat', strrep(datestr(t), ':', '_'));
filePath = fullfile(savedir, fileName);

if exist(filePath, 'file') == 2
    % File already exists, do not overwrite
    disp('File already exists. Not overwriting.');
else
    % File does not exist, proceed with saving
    save(filePath, 'CatchOutput', 'Catch', 'output', 'Conf');
    disp(['Variables saved to file: ', filePath]);
end
