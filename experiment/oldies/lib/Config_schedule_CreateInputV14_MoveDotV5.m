classdef Config_schedule_CreateInputV14_MoveDotV5
    properties(Constant)

        % CONFIG_SCHEDULE_CREATEINPUTV14_MOVEDOTV5
        %
        % Purpose:
        %   Centralize block/run schedule controls and catch-trial controls
        %   for:
        %     - experiment/CreateInputFiles_v14_threeRunsPerBlock_catch.m
        %     - experiment/MoveDot1_experiment_occlusion_v5_sequenceTriggers.m
        %     - experiment/MoveDot1_experiment_occlusion_v4_catchTrials.m
        %
        % Usage example (interactive):
        %   In MATLAB, open this file and set:
        %       Config_schedule_CreateInputV14_MoveDotV5.numBlocks
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

        % Number of blocks to schedule.
        numBlocks = 1;

        % Fixed number of runs per block for this schedule family.
        runsPerBlock = 3;

        % Fraction of each occluded condition assigned to run 2.
        % Run 3 receives the remaining fraction.
        run2FractionPerOccluded = 0.5;

        % Enforce exact equal source condition counts before catch insertion.
        % For the base design this implies:
        %   count(always_visible) == count(occluded_nondeviant) ==
        %   count(occluded_deviant), with occluded counts even.
        enforceEqualBaseRunLengths = true;

        % Catch rates per catch type.
        % Example with 20 base trials per condition and 0.10:
        %   run1 type-1 catches = round(20 * 0.10) = 2
        %   runs2+3 type-2 catches total = round((20+20) * 0.10) = 4
        catchRateType1Run1 = 0.10;
        catchRateType2Runs23 = 0.10;

        % Catch type 1 timing controls.
        catchType1DisappearRangeSec = [0.30, 1.00];
        catchType1InvisibleDurationSec = 0.50;
        catchType1ChangedPathProbability = 0.50;

        % Response prompt controls (shared by both catch types).
        catchQuestionText = 'Has the dot changed its course?';
        catchQuestionTimeoutSec = 4.0;
        catchResponseYesCode = 1;
        catchResponseNoCode = 2;

        % Input/output naming for additive v14 schedule artifacts.
        inputFilePattern = 'MovDot_Sub%02d.mat';
        outputFilePattern = 'Sub%02d_TrialStruct_v14_threeRunsPerBlock_catch.mat';
    end
end
