classdef Config_stimuli_generation_V21
    properties(Constant)

        % CONFIG_STIMULI_GENERATION_V21
        %
        % Purpose:
        %   Centralize one-dot trajectory generation controls used by:
        %     - experiment/stimuli_generation_V21.m
        %
        % Usage example (interactive):
        %   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment');
        %   addpath('lib/');
        %   stimuli_generation_V21;
        %
        % Usage example (non-interactive):
        %   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
        %   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
        %    addpath('lib/'); run('stimuli_generation_V21.m');"
        %
        % Key assumptions:
        %   - One-dot generation only.
        %   - Fixed-frame occlusion defaults match V21 expectations.
        %   - Output location points to experiment/input_files by default.

        % General one-dot generation controls.
        trialsPerCondition = 20;
        trialDuration = 2.67; % seconds

        frameFrequency = 120; % Hz

        % Keep one-dot generation in likelihood mode to match the current
        % occlusion paradigm assumptions.
        stimulusType = Utils.likelihood;

        % Dot geometry and kinematics (visual degrees).
        dotRectSize = [10, 10];
        dotWidth = 0.51442;
        dotSpeedDegPerSec = 3.73;
        dotSpeedDegPerFrame = Config_stimuli_generation_V21.dotSpeedDegPerSec / Config_stimuli_generation_V21.frameFrequency;

        % Deviance and curvature controls used by the config-driven
        % one-dot trajectory synthesis in stimuli_generation_V21.m.
        flipCurvatureOnDeviant = false;
        randomizeCurvatureOnDeviant = true;
        %deviantCurvatureRange = 0.35;
        initialCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        deviantCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        deviantSignedTurnWindows = [-81, 1; 1, 81];

        % Fixed-frame occlusion timing controls for the V21 paradigm.
        fixedDevianceFrame = 130;
        fixedOcclusionEndFrame = 190;

        % Single stimulus family kept for one-dot occlusion generation.
        likelihood = struct( ...
            'pathDuration', Config_stimuli_generation_V21.trialDuration, ...
            'directionVariance', [0, 45], ...
            'directionChange', 30, ...
            'deviantSignedTurnWindows', Config_stimuli_generation_V21.deviantSignedTurnWindows);

        % Save/load settings.
        inputDirectory = 'input_files/';
    end
end
