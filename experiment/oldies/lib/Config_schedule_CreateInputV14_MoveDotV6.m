classdef Config_schedule_CreateInputV14_MoveDotV6
    properties(Constant)

        % CONFIG_SCHEDULE_CREATEINPUTV14_MOVEDOTV6
        %
        % Purpose:
        %   Centralize block/run schedule controls and catch-trial controls
        %   for:
        %     - experiment/CreateInputFiles_v14_threeRunsPerBlock_catch.m
        %     - experiment/MoveDot1_experiment_occlusion_v6_runColorCue.m
        %     - experiment/MoveDot1_experiment_occlusion_v5_sequenceTriggers.m
        %
        % Usage example (interactive):
        %   In MATLAB, open this file and set:
        %       Config_schedule_CreateInputV14_MoveDotV6.numBlocks
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
        numBlocks = 1;

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

        % Run-color cue switch for experiment v6.
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
