classdef Config_occlusion_schedule_v1
    properties(Constant)

        % CONFIG_OCCLUSION_SCHEDULE_V1
        %
        % Purpose:
        %   Centralize schedule controls for the additive v13 trial-order
        %   builder and v2 occlusion runtime. This keeps existing files
        %   untouched while exposing one place where block/run scheduling is set.
        %
        % Usage example (interactive):
        %   In MATLAB, open this file and set:
        %       Config_occlusion_schedule_v1.numBlocks
        %   then run:
        %       CreateInputFiles_v13_threeRunsPerBlock
        %
        % Usage example (non-interactive):
        %   /Applications/MATLAB_R2020a.app/bin/matlab -batch ...
        %   "cd('/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment'); ...
        %    iSub=66; randomSeed=6601; run('CreateInputFiles_v13_threeRunsPerBlock.m');"
        %
        % Scheduling model:
        %   - Exactly 3 runs per block.
        %   - Run 1: all always_visible trials (shuffled).
        %   - Run 2: random half of occluded_nondeviant + random half of
        %            occluded_deviant, then shuffled together.
        %   - Run 3: remaining half from the same two occluded conditions,
        %            then shuffled together.
        %
        % Key assumptions:
        %   - Input datasets contain condition labels:
        %       always_visible, occluded_nondeviant, occluded_deviant
        %   - Occluded condition counts are equal and even.

        % Number of blocks to schedule.
        numBlocks = 1;

        % Fixed number of runs per block for this schedule family.
        runsPerBlock = 3;

        % Fraction of each occluded condition assigned to run 2.
        % Run 3 receives the remaining fraction.
        run2FractionPerOccluded = 0.5;

        % Enforce exact equal run lengths across all runs.
        % For the current design this implies:
        %   count(always_visible) == count(occluded_nondeviant) ==
        %   count(occluded_deviant), with occluded counts even.
        enforceEqualRunLengths = true;

        % Input/output naming for additive v13 schedule artifacts.
        inputFilePattern = 'MovDot_Sub%02d.mat';
        outputFilePattern = 'Sub%02d_TrialStruct_v13_threeRunsPerBlock.mat';
    end
end
