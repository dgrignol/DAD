classdef Config_schedule_CreateInputV14_MoveDotV8
    properties(Constant)

        % CONFIG_SCHEDULE_CREATEINPUTV14_MOVEDOTV8
        %
        % Purpose:
        %   Centralize block/run schedule controls and catch-trial controls
        %   for:
        %     - experiment/CreateInputFiles_v14_threeRunsPerBlock_catch.m
        %     - experiment/MoveDot1_experiment_occlusion_v8_runColorCueMessages.m
        %     - experiment/MoveDot1_experiment_occlusion_v5_sequenceTriggers.m
        %
        % Usage example (interactive):
        %   In MATLAB, open this file and set:
        %       Config_schedule_CreateInputV14_MoveDotV8.numBlocks
        %   then run:
        %       CreateInputFiles_v14_threeRunsPerBlock_catch
        %
        % Usage example (non-interactive):
        %   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
        %   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
        %    iSub=66; randomSeed=6601; run('CreateInputFiles_v14_threeRunsPerBlock_catch.m');"
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
        numBlocks = 2;

        % Runs per block (kept fixed for this family).
        % Example: 3 -> run1(always_visible), run2/3(occluded mix).
        runsPerBlock = 3;

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
        catchRateType1Run1 = 0.10;
        % Catch fraction across pooled runs 2+3 occluded base trials.
        % Example: 0.10 with 40 pooled base trials -> 4 type-2 catches total.
        catchRateType2Runs23 = 0.10;

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

        % Start-of-experiment message shown before the first run starts.
        % Set duration <= 0 to require an explicit keyboard/button-box
        % response (1 or 8) with no timeout.
        % Enable pre-run message before run 1.
        % Example: true shows welcome/instruction screen.
        startMessageEnabled = true;
        % Text shown in the pre-run message.
        % Example: includes button-box instructions.
        startMessageText = sprintf(['Welcome!\n\n' ...
            'Press Red or Blue on the button box\n' ...
            'to start the experiment.']);
        % Message duration (seconds). <=0 means wait for response.
        % Example: 0.0 -> require explicit response (1 or 8).
        startMessageDurationSec = 0.0;

        % Mid-experiment transition message shown after run 1 and before
        % the following run (typically run 2).
        % Enable transition message between run 1 and run 2.
        % Example: true adds a short rest prompt after run 1.
        run1TransitionMessageEnabled = true;
        % Text shown in the run1->run2 transition screen.
        % Example: remind participant to rest before continuing.
        run1TransitionMessageText = sprintf(['End of first task.\n\n' ...
            'You ca rest your eyes.\n' ...
            'Second task will start soon.']);
        % Transition message duration in seconds.
        % Example: 5.0 -> auto-continue after 5 s.
        run1TransitionMessageDurationSec = 5.0;

        % End-of-block message (non-final blocks).
        % Runtime fills the two %d placeholders with:
        %   1) ended block number
        %   2) remaining block count
        % Enable message shown between non-final blocks.
        % Example: true shows "block ended / remaining blocks" prompt.
        endOfBlockMessageEnabled = true;
        % End-of-block template with placeholders:
        %   first %d = completed block index, second %d = remaining blocks.
        endOfBlockMessageTextTemplate = sprintf(['Block %d ended.\n' ...
            'Remaining blocks: %d.\n\n' ...
            'Press Red or Blue to start the next block.']);
        % End-of-block message duration. <=0 waits for response.
        % Example: 0.0 -> participant must press button to continue.
        endOfBlockMessageDurationSec = 0.0;

        % End-of-experiment message shown after the last block.
        % Enable final message after last block.
        % Example: true shows closing text before quit.
        finalMessageEnabled = true;
        % Closing text at experiment end.
        % Example: replace with lab-standard thank-you wording.
        finalMessageText = 'you completed the experiment bye thanks';
        % Final message duration in seconds.
        % Example: 5.0 -> auto-close after 5 s.
        finalMessageDurationSec = 5.0;

        % Shared controls for message transitions.
        % Add fixation period immediately after message screens.
        % Example: true enforces short recentering interval.
        postMessageFixationEnabled = true;
        % Duration of post-message fixation interval.
        % Example: 3.0 seconds.
        postMessageFixationDurationSec = 3.0;

        % For message screens that require button confirmation (except
        % catch trials), wait this long before accepting inputs.
        % Ignore button responses during initial lockout period.
        % Example: 1.0 prevents accidental carry-over presses.
        messageResponseLockoutSec = 1.0;

        % Run-color cue switch for experiment v8.
        % If true, run 1 and runs 2/3 use different dot colors with
        % odd/even subject counterbalancing in runtime.
        % Enable run-family color cue metadata used by runtime.
        % Example: true -> run 1 color differs from runs 2/3.
        runColorCueEnabled = true;

        % Input/output naming for additive v14 schedule artifacts.
        % Input stimulus filename pattern.
        % Example: subject 66 -> MovDot_Sub66.mat
        inputFilePattern = 'MovDot_Sub%02d.mat';
        % Output schedule filename pattern.
        % Example: subject 66 -> Sub66_TrialStruct_v14_threeRunsPerBlock_catch.mat
        outputFilePattern = 'Sub%02d_TrialStruct_v14_threeRunsPerBlock_catch.mat';
    end
end
