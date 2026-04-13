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

        % Number of unique sequence IDs generated per condition.
        % Example: 20 -> sequence labels 1..20 in each condition family.
        trialsPerCondition = 20;
        % Trial duration in seconds used to derive frame count.
        % Example: 2.67 s at 120 Hz -> 320 frames/trial.
        trialDuration = 2.67; % seconds

        % Display/generation frame rate assumption.
        % Example: keep 120 Hz for MEG presentation timing alignment.
        frameFrequency = 120; % Hz

        % Stimulus family selector (kept in likelihood mode for occlusion).
        % Example: Utils.likelihood enables [0,45] condition coding.
        stimulusType = Utils.likelihood;

        % Motion arena size in visual degrees [width height].
        % Example: [10, 10] means a 10x10 deg motion window.
        dotRectSize = [10, 10];
        % Dot diameter in visual degrees.
        % Example: 0.51442 deg (used for bounds and visibility geometry).
        dotWidth = 0.51442;
        % Dot speed in visual degrees per second.
        % Example: 3.73 deg/s at 120 Hz -> ~0.0311 deg/frame.
        dotSpeedDegPerSec = 3.73;
        % Derived framewise displacement (deg/frame).
        % Example: consumed directly by path integrator each frame.
        dotSpeedDegPerFrame = Config_stimuli_generation_V21.dotSpeedDegPerSec / Config_stimuli_generation_V21.frameFrequency;

        % Deviance and curvature controls used by the config-driven
        % one-dot trajectory synthesis in stimuli_generation_V21.m.
        % If true, invert curvature sign after deviance onset.
        % Example: false keeps sign unchanged unless randomization is enabled.
        flipCurvatureOnDeviant = false;
        % If true, resample post-onset curvature from deviant windows.
        % Example: true introduces post-onset curvature variability.
        randomizeCurvatureOnDeviant = true;
        %deviantCurvatureRange = 0.35;
        % Allowed baseline curvature intervals (deg/frame-equivalent units).
        % Example: excludes near-zero curvature to avoid nearly straight paths.
        initialCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        % Allowed post-deviance curvature intervals when randomization is on.
        % Example: match baseline windows to keep magnitude regimes balanced.
        deviantCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        % Allowed signed deviant-turn windows in degrees at onset.
        % Example: [-81,1;1,81] excludes exactly-zero signed turn.
        deviantSignedTurnWindows = [-81, 1; 1, 81];

        % Fixed-frame occlusion timing controls for the V21 paradigm.
        % Fixed frame index where deviance and full occlusion begin.
        % Example: 130 at 120 Hz is ~1.08 s after trial start.
        fixedDevianceFrame = 130;
        % Last frame guaranteed to be fully occluded (inclusive).
        % Example: 190 keeps full occlusion through ~1.58 s.
        fixedOcclusionEndFrame = 190;

        % Single stimulus family kept for one-dot occlusion generation.
        % Likelihood-mode condition definition consumed by generator.
        % Example: directionVariance [0,45] -> nondeviant + deviant condition.
        likelihood = struct( ...
            'pathDuration', Config_stimuli_generation_V21.trialDuration, ...
            'directionVariance', [0, 45], ...
            'directionChange', 30, ...
            'deviantSignedTurnWindows', Config_stimuli_generation_V21.deviantSignedTurnWindows);

        % Output folder (relative to experiment/ unless absolute path is used).
        % Example: 'input_files/' writes MovDot_SubXX*.mat there.
        inputDirectory = 'input_files/';
    end
end
