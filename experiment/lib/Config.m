classdef Config
    properties(Constant)
    
        %general
        trialsPerCondition = 50;%20 %10000;  %nSeq in their code
        trialDuration = 2.67;       %in seconds
        trialRepetetion = 1;
        %general - pc settings
        frameFrequency = 120;   %frames per senconds
        %in their code they define the inter frame interval, which is not used
        framesPerTrial = round(Config.trialDuration*Config.frameFrequency);

        stimulusType = Utils.likelihood;
        % stimulusType = Utils.path_duration;

        %curvyness
        curvFactor = 0.8;%0.4%0.9 %1.05 %new parameter for 120Hz. Old parameter 2.1;            %added to the direction angle in case of constant #HERE 0.45
        isCurvValenceRand = true;   %if set the curvyness valence (ie +/- will be random from one path to the other), else will alternate
        isCurvFactorRand = true;    %if set the curvFactor will be multiplied by a rand before added

        %dot settings
        dotRectSize = [10, 10]; %[14, 14]                 %rect size where the dot can move in visual degrees
        dotWidth = 0.51442;%Config.dotRectSize(1)/28;    %in degrees
        dotSpeedDegPerSec = 3.73;                                              %visual angle degrees per second
        dotSpeedDegPerFrame = Config.dotSpeedDegPerSec / Config.frameFrequency; %visual angle degrees per frame
        %their code computes the dot speed pixels per frame but not used
        dotRectColor = [10 10 10];
        pathDurationVariance = 0;
        isConnectedPaths = false;               %the end of previous path is used for the new one if flag is set

        minDistanceBetweenDots = Config.dotWidth*2; %in degrees, how much should be apart from each other the dots
        %dot settings - deviant
        deviantOnset = 0.5; %midway the path duration
        deviantOnsetVariance = 0;
        flipCurvatureOnDeviant = true; %set true to flip curvature sign at deviant onset
        randomizeCurvatureOnDeviant = true; %set true to sample a new post-deviant curvature value (deviant path only) #HERE flase
        deviantCurvatureRange = 0.35; %post-deviant curvature is sampled uniformly in [-range, +range]
        % Explicit curvature windows (deg/frame) for stimuli_generation_v15.m.
        % Format: one row per allowed interval [minCurvature, maxCurvature].
        % These defaults enforce |curvature| >= 0.3755 and <= 0.8.
        initialCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        deviantCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        % Optional explicit deviant-turn windows (signed degrees) for
        % likelihood stimulus generation in stimuli_generation_v14.m.
        % Format: one row per allowed interval [minSignedTurn, maxSignedTurn].
        % Example below enforces |turn| >= 45 deg and allows both CW/CCW:
        %   [-180, -45] U [45, 180]
        % Set to [] to keep legacy generation from directionVariance.
        deviantSignedTurnWindows = [-80, -20; 20, 80];
        % Optional deviant displacement region for stimuli_generation_v16_Displacement.m.
        % Radius is [innerRadius, outerRadius] in visual degrees.
        % Angles are signed degrees relative to pre-deviant heading.
        % Displacement mode selects how radius/angle controls are interpreted:
        %   - 'off': disable displacement.
        %   - 'constrained': use radius + angle windows as configured.
        %   - 'freeStart': ignore configured radius and angle windows by
        %     forcing full-circle angles [-180, 0; 0, 180] and a
        %     screen-scale radius [0, max(dotRectSize)].
        % Set radius to [0, 0] or angle windows to [] to disable displacement
        % when using 'constrained'.
        deviantDisplacementMode = 'freeStart';
        deviantDisplacementRadiusRangeDeg = [0.5, 1.8];
        deviantDisplacementAngleWindowsDeg = [-140, -40; 40, 140];
        % dRSA proxy-gate controls for stimuli_generation_v14+/v17 style
        % acceptance filtering.
        % - enableDrsaProxyGate toggles gate use in all conditions.
        % - applyDrsaProxyGateToNondeviantOnly limits the gate to
        %   nondeviant conditions when enabled.
        enableDrsaProxyGate = true;
        applyDrsaProxyGateToNondeviantOnly = false;

        numXGrids = 2;              %num of grids in the rect of the dot
        numYGrids = 2;
        spreadoutTolerance = 0.5;   %the accepted spreadout unequality
        
        xGridSize = Config.dotRectSize(1)/Config.numXGrids; %single X grid size of the dotRect in degrees
        yGridSize = Config.dotRectSize(2)/Config.numYGrids; %single Y grid size of the dotRect in degrees
        
        %get the perfect spread first and then add and subtract the margin
        %tolerance to create the acceptance interval
        perfectSpreadout = Config.framesPerTrial/(Config.numXGrids*Config.numYGrids);
        spreadoutToleranceInterval = [Config.perfectSpreadout - Config.perfectSpreadout*Config.spreadoutTolerance ...
            Config.perfectSpreadout + Config.perfectSpreadout*Config.spreadoutTolerance];

        %this one might be changed or ommited, they use it to avoid extrem start positions
        %Update: they use it keep the start point of the path more around
        %the center
        focusAroundCenterFactor = 0.5;%0.2; %was 0.3 before

        %dot settings - likelihood
        %was before 'pathDuration', 1.0, ... 'directionVariance', [0, 10, 20, 30]
        likelihood = struct('pathDuration', 2.6667, ...    %in seconds
            'directionVariance', [0, 45], ...   %represents the predictability
            'directionChange', 30, ... %not used in this type
            'deviantSignedTurnWindows', Config.deviantSignedTurnWindows, ... %explicit signed deviant turns (optional)
            'deviantDisplacementMode', Config.deviantDisplacementMode, ... %off|constrained|freeStart
            'deviantDisplacementRadiusRangeDeg', Config.deviantDisplacementRadiusRangeDeg, ... %annulus bounds (optional)
            'deviantDisplacementAngleWindowsDeg', Config.deviantDisplacementAngleWindowsDeg); %signed angle windows (optional)

        %dot settings - path_duration
        % was before [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.6, 2, 2.5, 3.2]
        path_duration = struct('pathDuration', [2.6667], ...    %in seconds
            'directionVariance', 0, ...   %represents the predictability
            'directionChange', 30); %not used in this type

        %dot settings - path_duration_norm
        path_duration_norm = struct('pathDuration', [0.1333, 0.6667, 1.2000, 1.6000], ...    %in seconds
            'directionVariance', 0, ...   %represents the predictability
            'directionChange', 100); 

        %dialogs
        dialogTitle = 'Generation parameters';
        dialogPrompts = {'Enter the viewing distance (in millimeters):', ...
            'Enter the subject ID:'};
        dialogDimensions = [1 100];
        dialogDefaults = {'1000', '1'};

        %screen settings
        screenRect = [800 50 2300 1550]; %was before [0 0 800 800]
        screenBackgroundColor = [128 128 128];

        %behavioral settings
        crossDefaultColor = [0 1 0];%[0 0 1]
        crossChangeColor = [0.5 0 1];
        crossTargetColor = [0 1 0];

        dot1Color = [1 0 0];        %red
        dot2Color = [0.9 0.9 0.9];  %whitish

        %save/load settings
        inputDirectory = 'input_files/';
        outputDirectory = 'output_files/'; 
        stimuliFileName = 'MovDot_Sub%02d.mat';
    end
end
