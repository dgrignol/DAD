classdef Config_schedule_CreateInputV21_MoveDotV18_rescueTraject
    properties(Constant)

        % CONFIG_SCHEDULE_CREATEINPUTV21_MOVEDOTV18_RESCUETRAJECT
        %
        % Purpose:
        %   Centralize block/run schedule controls and catch-trial controls
        %   for:
        %     - experiment/CreateInputFiles_v21_rescueTraject.m
        %     - experiment/MoveDot1_experiment_occlusion_v18_rescueTraject.m
        %
        % Usage example (interactive):
        %   In MATLAB, open this file and set:
        %       Config_schedule_CreateInputV21_MoveDotV18_rescueTraject.numBlocks
        %   then run:
        %       CreateInputFiles_v21_rescueTraject
        %
        % Usage example (non-interactive):
        %   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
        %   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
        %    iSub=66; randomSeed=6601; ...
        %    run('CreateInputFiles_v21_rescueTraject.m');"
        %
        % Scheduling model:
        %   - Exactly 3 runs per block.
        %   - Base (non-catch) schedule:
        %       Run 1 -> always_visible only (shuffled)
        %       Run 2 -> first random half of each occluded condition
        %       Run 3 -> remaining random half of each occluded condition
        %   - Catch trial insertion:
        %       Type 1 catches are inserted into run 1 only.
        %       Type 2 catches are inserted across runs 2 and 3.
        %
        % Key assumptions:
        %   - Input datasets contain condition labels:
        %       always_visible, occluded_nondeviant, occluded_deviant
        %   - Base occluded condition counts are equal and even.
        %   - Catch rates are fractions in [0, 1].

        % Number of blocks generated in TrialOrder.
        % Example: 2 -> block 1 and block 2 each with runs 1..3.
        numBlocks = 1;

        % Runs per block (kept fixed for this exp version family).
        % Example: 3 -> run1(always_visible), run2/3(occluded mix).
        runsPerBlock = 3; % DO NOT CHANGE THIS

        % Fraction of each occluded condition sent to run 2.
        % Example: 0.5 splits 20 trials as 10 in run 2, 10 in run 3.
        run2FractionPerOccluded = 0.5;

        % Enforce exact equal source condition counts before catch insertion.
        % For the base design this implies:
        %   count(always_visible) == count(occluded_nondeviant) ==
        %   count(occluded_deviant), with occluded counts even.
        % Strict base-balance check before catch insertion.
        % Example: true errors if condition counts are mismatched.
        enforceEqualBaseRunLengths = true;

        % Catch rates per catch type.
        % Example with 20 base trials per condition and 0.10:
        %   run1 type-1 catches = round(20 * 0.10) = 2
        %   runs2+3 type-2 catches total = round((20+20) * 0.10) = 4
        % Catch fraction for run-1 (always_visible) base trials.
        % Example: 0.10 with 20 base trials -> 2 type-1 catches.
        catchRateType1Run1 = 0.00;
        % Catch fraction across pooled runs 2+3 occluded base trials.
        % Example: 0.10 with 40 pooled base trials -> 4 type-2 catches total.
        catchRateType2Runs23 = 0.00;

        % Catch type 1 timing controls.
        % Allowed disappearance onset range for type-1 catches.
        % Example: [0.30,1.00] seconds from trial start.
        catchType1DisappearRangeSec = [0.30, 1.00];
        % Invisible duration before reappearance in type-1 catches.
        % Example: 0.50 -> dot hidden for 500 ms.
        catchType1InvisibleDurationSec = 0.50;
        % Probability that a type-1 catch reappears on changed path.
        % Example: 0.50 -> half changed, half unchanged on average.
        catchType1ChangedPathProbability = 0.50;

        % Response prompt controls (shared by both catch types).
        % Prompt shown for catch decisions.
        % Example: customize to local language or shorter wording.
        catchQuestionText = 'Has the dot changed its course?';
        % Response window for catch question.
        % Example: 4.0 -> timeout after 4 seconds.
        catchQuestionTimeoutSec = 4.0;
        % Logical code used by scheduler for YES correctness checks.
        % Example: 1 means "changed".
        catchResponseYesCode = 1;
        % Logical code used by scheduler for NO correctness checks.
        % Example: 2 means "not changed".
        catchResponseNoCode = 2;

        % Practice mode: percentage of scheduled trials to keep per run
        % when practice mode is enabled from the runtime input dialog.
        % Example:
        %   100 -> full run length (all selected trials become catch trials)
        %    15 -> keep 15%% of each run (all selected trials are catch).
        practiceRunPercentageTrials = 15;

        % Start-of-experiment message shown before the first run starts.
        % Set duration <= 0 to require an explicit keyboard/button-box
        % response (1 or 8) with no timeout.
        % Enable pre-run message before run 1.
        % Example: true shows welcome/instruction screen.
        startMessageEnabled = false;
        % Text shown in the pre-run message.
        % Example: includes button-box instructions.
        startMessageText = sprintf(['Welcome!\n\n' ...
            'Press Red/Blue on the button box\n' ...
            'or key 1/8 on keyboard\n' ...
            'to start the experiment.']);
        % Message duration (seconds). <=0 means wait for response.
        % Example: 0.0 -> require explicit response (1 or 8).
        startMessageDurationSec = 0.0;

        % Mid-experiment transition message shown after run 1 and before
        % the following run (typically run 2).
        % Enable transition message between run 1 and run 2.
        % Example: true adds a short rest prompt after run 1.
        run1TransitionMessageEnabled = false;
        % Text shown in the run1->run2 transition screen.
        % Example: remind participant to rest before continuing.
        run1TransitionMessageText = sprintf(['End of first task.\n\n' ...
            'You can rest your eyes.\n' ...
            'Second task will start soon.']);
        % Transition message duration in seconds.
        % Example: 5.0 -> auto-continue after 5 s.
        run1TransitionMessageDurationSec = 0.0;
        % If true, append a visible seconds countdown to the run1->run2 message.
        % Example: true shows "Second task starts in: 5 s ... 1 s".
        run1TransitionCountdownEnabled = false;

        % End-of-block message (non-final blocks).
        % Runtime fills the two %d placeholders with:
        %   1) ended block number
        %   2) remaining block count
        % Enable message shown between non-final blocks.
        % Example: true shows "block ended / remaining blocks" prompt.
        endOfBlockMessageEnabled = false;
        % End-of-block template with placeholders:
        %   first %d = completed block index, second %d = remaining blocks.
        % IMPORTANT:
        %   Keep this as a raw format template (do not wrap in sprintf here).
        %   Runtime injects block values with sprintf at display time.
        endOfBlockMessageTextTemplate = ['Block %d ended.\n' ...
            'Remaining blocks: %d.\n\n' ...
            'Press Red or Blue to start the next block.'];
        % End-of-block message duration. <=0 waits for response.
        % Example: 0.0 -> participant must press button to continue.
        endOfBlockMessageDurationSec = 0.0;
        % If true, show a second screen after the end-of-block continue gate:
        % skip calibration or run calibration now, then continue.
        betweenBlockCalibrationChoiceEnabled = false;
        % Between-block calibration choice text.
        % Response mapping:
        %   - NO CALIBRATION: Red / key 1
        %   - CALIBRATE NOW: Blue / key 8
        betweenBlockCalibrationChoiceText = sprintf(['Calibration choice.\n\n' ...
            'Press Red / key 1 to continue without calibration.\n' ...
            'Press Blue / key 8 to run calibration now.\n\n' ...
            'Then the next block starts.']);
        % Optional timeout for calibration-choice screen (seconds).
        % Use <=0 to wait until a valid key/button is received.
        betweenBlockCalibrationChoiceDurationSec = 0.0;

        % End-of-experiment message shown after the last block.
        % Enable final message after last block.
        % Example: true shows closing text before quit.
        finalMessageEnabled = false;
        % Closing text at experiment end.
        % Example: replace with lab-standard thank-you wording.
        finalMessageText = 'you completed the experiment bye thanks';
        % Final message duration in seconds.
        % Example: 5.0 -> auto-close after 5 s.
        finalMessageDurationSec = 0.0;

        % Shared controls for message transitions.
        % Add fixation period immediately after message screens.
        % Example: true enforces short recentering interval.
        postMessageFixationEnabled = false;
        % Duration of post-message fixation interval.
        % Example: 3.0 seconds.
        postMessageFixationDurationSec = 0.0;

        % For message screens that require button confirmation (except
        % catch trials), wait this long before accepting inputs.
        % Ignore button responses during initial lockout period.
        % Example: 1.0 prevents accidental carry-over presses.
        messageResponseLockoutSec = 0.0;

        % Run-color cue switch for experiment v18.
        % If true, run 1 and runs 2/3 use different dot colors with
        % odd/even subject counterbalancing in runtime.
        % Enable run-family color cue metadata used by runtime.
        % Example: true -> run 1 color differs from runs 2/3.
        runColorCueEnabled = true;

        % Path-band terminal rendering style for runtime occluder drawing.
        % Supported values:
        %   'round'    -> rounded start/end caps at both occluder terminals.
        %   'straight' -> wall-like straight start/end at both terminals.
        % Example: 'straight' for wall-like terminal edges at deviance and reappearance.
        pathBandEntranceStyle = 'straight';

        % Backshift amount (in dot-radius units) applied when
        % pathBandEntranceStyle='straight'.
        %
        % Timing note:
        %   The default is 0.0 so complete occlusion stays locked to the
        %   deviance frame (fixedDevianceFrame). Positive values move the
        %   straight terminal backward along the path and can make complete
        %   occlusion happen earlier.
        %
        % Example:
        %   0.0 -> no backshift (default, timing-safe).
        %   0.5 -> half-radius backward shift.
        pathBandStraightBackshiftDotRadiusScale = 0.0;

        % Forward extension amount (in dot-radius units) applied at
        % straight occluder exits to keep the last full-occlusion frame
        % aligned with fixedOcclusionEndFrame.
        %
        % Example:
        %   1.0 -> extend the wall by one dot radius (default).
        pathBandStraightEndForwardDotRadiusScale = 1.0;

        % Perpendicular width boost (in dot-radius units) applied only at
        % the post-occluder entrance wall. This increases local vertical
        % coverage without shifting entrance timing along path tangent.
        %
        % Example:
        %   0.15 -> small entrance thickening (default).
        pathBandStraightPostStartWidthExtraDotRadiusScale = 0.15;

        % Pixel correction added to straight start backshift after
        % dot-radius scaling. This compensates display-dependent
        % quantization so first full occlusion stays aligned with
        % fixedDevianceFrame.
        %
        % Example:
        %   -0.5 -> reduce start backshift by half pixel (default).
        pathBandStraightStartBackshiftPixelOffset = -0.5;

        % Psychtoolbox sync-test policy for runtime presentation.
        %
        % Rationale:
        %   - MEG sessions should not force SkipSyncTests.
        %   - Development/debug sessions can opt in to SkipSyncTests to run
        %     on non-calibrated displays.
        %
        % Example:
        %   skipSyncTestsWhenMEG = false keeps MEG strict by default.
        skipSyncTestsWhenMEG = false;
        % When true, runtime enables SkipSyncTests automatically in
        % non-MEG sessions only when a debug visual/trigger mode is active.
        % Example: true keeps development convenient while preserving MEG strictness.
        skipSyncTestsInDebug = true;

        % Flip-miss diagnostics guardrails (runtime frame pacing).
        % Count threshold that triggers a warning message during a trial.
        % Example: 1 warns on first missed flip.
        flipMissedWarnThresholdCount = 1;
        % If true, abort trial loop when missed-flip count reaches threshold.
        % Example: false keeps warn-only default behavior.
        flipMissedAbortEnabled = false;
        % Abort threshold used only when flipMissedAbortEnabled is true.
        % Example: 5 aborts after five missed flips in one trial.
        flipMissedAbortThresholdCount = 5;
        % Positive epsilon for classifying missed flips.
        % Example: 1e-6 ignores tiny numerical noise around zero.
        flipMissedEpsilonSec = 1e-6;

        % Eye-tracker controls for MoveDot1_experiment_occlusion_v18_rescueTraject.
        % Master switch for eye-tracker support.
        % Example: true enables EyeLink setup and trial recording.
        eyeTrackerEnabled = false;
        % If true, disables EyeLink I/O and fixation-break logic even when
        % eyeTrackerEnabled is true.
        ignoreEyeTracker = false;
        % If true, simulate gaze samples to exercise fixation-break and
        % replay logic without EyeLink hardware.
        fakeEyeTracker = false;
        % Per-trial probability of a simulated fixation break in fake mode.
        % Example: 0.05 means 5%% of trials schedule a fake break.
        fakeGazeBreakRate = 0.05;
        % Enable fixation-break detection and replay scheduling (trigger 150/151).
        enableFixationAbort = true;
        % Fixation acceptance radius in visual degrees.
        % Example: 7 deg from fixation center.
        fixWindowDeg = 7.0;
        % Consecutive outside-window frames required to flag a fixation break.
        fixBreakToleranceFrames = 5;
        % Warning duration after a fixation break, in seconds.
        fixationWarningDurationSec = 2.0;
        % Replay cap per source trial after fixation break.
        % Example: 1 means at most one replay of that source trial.
        replayLimit = 1;
        % If true, ignore replayLimit and allow unlimited replay queueing.
        allowInfiniteReplays = false;
        % If true, run EyeLink tracker setup/calibration before the session.
        eyeTrackerDoCalibration = true;
        % If true, mirror experiment trigger events as EyeLink messages.
        eyeTrackerSendTriggerMessages = true;
        % If true, transfer fixation bitmap to EyeLink host display.
        eyeTrackerImageTransferEnabled = true;
        % Small pause between EyeLink setup/recording commands.
        eyeTrackerWaitSec = 0.01;
        % EyeLink EDF file policy.
        % Allowed values:
        %   - 'per_block': one EDF per completed block (recommended).
        % The runtime opens and receives one EDF for each block before showing
        % the block-end rest/continue message.
        eyeTrackerEdfPolicy = 'per_block';
        % If true, attempt EDF retrieval immediately after each block.
        eyeTrackerReceiveFileAfterBlock = true;
        % If true and eye tracking is enabled, allow between-block calibration.
        eyeTrackerAllowBetweenBlockCalibration = false;
        % Path to fixation bitmap used for EyeLink host image transfer.
        eyeTrackerFixBmpPath = './fixation.bmp';
        % Dedicated trigger code for termination by ESC.
        % Keep this outside sequence range and outside existing fixed triggers.
        triggerEscTermination = 152;

        % Input/output naming for additive v21 rescueTraject schedule artifacts.
        % Input stimulus filename pattern.
        % Example: subject 66 -> MovDot_Sub66_V28_rescueTraject.mat
        inputFilePattern = 'MovDot_Sub%02d_V28_rescueTraject.mat';
        % Output schedule filename pattern.
        % Example: subject 66 -> Sub66_TrialStruct_v21_rescueTraject.mat
        outputFilePattern = 'Sub%02d_TrialStruct_v21_rescueTraject.mat';
    end
end
