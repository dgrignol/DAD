classdef Config_stimuli_generation_V26_eyeTrackerReplay
    properties(Constant)

        % CONFIG_STIMULI_GENERATION_V26_PATHBANDOCCLUDER
        %
        % Purpose:
        %   Centralize one-dot trajectory generation controls used by:
        %     - experiment/stimuli_generation_V26_eyeTrackerReplay.m
        %
        % Usage example (interactive):
        %   cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment');
        %   addpath('lib/');
        %   stimuli_generation_V26_eyeTrackerReplay;
        %
        % Usage example (non-interactive):
        %   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
        %   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
        %    addpath('lib/'); run('stimuli_generation_V26_eyeTrackerReplay.m');"
        %
        % Key assumptions:
        %   - One-dot generation only.
        %   - Fixed-frame occlusion defaults keep first full occlusion at
        %     fixedDevianceFrame.
        %   - Path-band width is strictly larger than the dot diameter so
        %     the dot is fully covered at the deviance frame.

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
        dotSpeedDegPerFrame = Config_stimuli_generation_V26_eyeTrackerReplay.dotSpeedDegPerSec / ...
            Config_stimuli_generation_V26_eyeTrackerReplay.frameFrequency;

        % Path-band occluder thickness controls.
        % Width multiplier applied to dot diameter.
        % Example: 1.10 -> band width is 10% larger than dot diameter.
        pathBandWidthMultiplier = 1.10;
        % Absolute safety margin (deg) to guarantee width > dot diameter.
        % Example: 0.005 deg keeps strict inequality even if multiplier=1.
        pathBandMinMarginDeg = 0.005;
        % Final occluder width used by path-band geometry (deg).
        % Example: max(dotWidth*multiplier, dotWidth+margin).
        pathBandWidthDeg = max( ...
            Config_stimuli_generation_V26_eyeTrackerReplay.dotWidth * ...
            Config_stimuli_generation_V26_eyeTrackerReplay.pathBandWidthMultiplier, ...
            Config_stimuli_generation_V26_eyeTrackerReplay.dotWidth + ...
            Config_stimuli_generation_V26_eyeTrackerReplay.pathBandMinMarginDeg);

        % Path-band terminal style used by geometry-based visibility metadata.
        % Supported values:
        %   'round'    -> rounded start/end caps at band terminals.
        %   'straight' -> straight wall-like terminal edges at both ends.
        % Example: 'straight' aligns metadata with runtime straight terminals.
        pathBandTerminalStyle = 'straight';
        % Backshift amount (in dot-radius units) for straight terminal style.
        %
        % Timing note:
        %   Keep this at 0.0 by default to preserve first complete occlusion
        %   at fixedDevianceFrame. Positive values can shift complete
        %   occlusion earlier.
        %
        % Example: 0.0 -> no backward shift (default).
        pathBandStraightBackshiftDotRadiusScale = 0.0;

        % Optional central-fixation collision handling for generated paths.
        % Modes:
        %   - 'off'   : no fixation-zone handling.
        %   - 'retry' : reject colliding candidates and resample.
        %   - 'move'  : try to translate colliding candidates minimally out
        %               of the exclusion zone before falling back to retry.
        % Collision handling mode around central fixation zone.
        % Example: 'move' (default here), alternatives: 'retry', 'off'.
        fixationCollisionMode = 'move';
        % Approximate fixation cross full span in visual degrees.
        % Example: if runtime has Conf.fixSizeDeg = 0.2, set this to ~0.4.
        fixationCrossSizeDegApprox = 0.4;
        % Extra conservative clearance around fixation, beyond
        % cross half-span + dot radius, to prevent near-miss trajectories.
        % Example: 0.10 keeps the dot visibly farther from the cross center.
        fixationSafetyMarginDeg = 0.10;
        % Exclusion radius for the dot center around fixation.
        % Example: cross half-span + dot radius + safety margin.
        fixationExclusionRadiusDeg = ...
            (Config_stimuli_generation_V26_eyeTrackerReplay.fixationCrossSizeDegApprox / 2) + ...
            (Config_stimuli_generation_V26_eyeTrackerReplay.dotWidth / 2) + ...
            Config_stimuli_generation_V26_eyeTrackerReplay.fixationSafetyMarginDeg;
        % Extra safety padding added to the move-mode target radius.
        % Example: 0.01 deg ensures a small clearance margin.
        fixationMovePaddingDeg = 0.01;
        % Number of direction candidates tested in move mode.
        % Example: 32 gives dense angular search with low runtime cost.
        fixationMoveDirectionSamples = 32;
        % Number of shift magnitudes sampled per candidate direction.
        % Example: 320 supports fine-grained minimal-shift search.
        fixationMoveShiftSamples = 320;

        % Deviance and curvature controls used by config-driven synthesis.
        % If true, invert curvature sign after deviance onset.
        % Example: false keeps sign unchanged unless randomization is enabled.
        flipCurvatureOnDeviant = false;
        % If true, resample post-onset curvature from deviant windows.
        % Example: true introduces post-onset curvature variability.
        randomizeCurvatureOnDeviant = true;
        % Allowed baseline curvature intervals (deg/frame-equivalent units).
        % Example: excludes near-zero curvature to avoid nearly straight paths.
        initialCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        % Allowed post-deviance curvature intervals when randomization is on.
        % Example: match baseline windows to keep magnitude regimes balanced.
        deviantCurvatureWindows = [-0.8, -0.3755; 0.3755, 0.8];
        % Allowed signed deviant-turn windows in degrees at onset.
        % Example: [-81,-10;10,81] excludes exactly-zero signed turn.
        deviantSignedTurnWindows = [-81, -10; 10, 81];

        % Fixed-frame occlusion timing controls for the V26 paradigm.
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
            'pathDuration', Config_stimuli_generation_V26_eyeTrackerReplay.trialDuration, ...
            'directionVariance', [0, 45], ...
            'directionChange', 30, ...
            'deviantSignedTurnWindows', ...
            Config_stimuli_generation_V26_eyeTrackerReplay.deviantSignedTurnWindows);

        % Output folder (relative to experiment/ unless absolute path is used).
        % Example: 'input_files/' writes MovDot_SubXX*.mat there.
        inputDirectory = 'input_files/';

        % Versioned output filenames for the path-band occluder branch.
        % Example: subject 66 -> MovDot_Sub66_V26_eyeTrackerReplay.mat
        observedFilePattern = 'MovDot_Sub%02d_V26_eyeTrackerReplay.mat';
        % Example: subject 66 -> MovDot_Sub66_V26_eyeTrackerReplay_predicted.mat
        predictedFilePattern = 'MovDot_Sub%02d_V26_eyeTrackerReplay_predicted.mat';
    end
end
