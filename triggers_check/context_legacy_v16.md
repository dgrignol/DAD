# Legacy trigger-check context: v16

> Historical reference for the v16 trigger investigation. The experiment,
> configuration, and data paths below describe that investigation's scope.
> Use the [current preprocessing roadmap](../meg/docs/COHORT_PREPROCESSING_REVIEW_GUIDE.md)
> for the cohort workflow and resolve each participant's actual runtime from
> their saved sources and profile.

This historical investigation contains MEG data from an experiment, together
with its related input files and experiment outputs.

## Historical experiment references

The experiment references used for this investigation were:

- `stimuli_generation_V26_PathBandOccluder.m`
- `CreateInputFiles_v19_threeRunsPerBlock_catch_PathBandOccluder.m`
- `MoveDot1_experiment_occlusion_v16_PathBandOccluder.m`
- `trigger_codes_occlusion_v7.md`

These files are located in:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/oldies`

## Historical configuration references

The configuration references used for this investigation were:

- `Config_stimuli_generation_V26_PathBandOccluder.m`
- `Config_schedule_CreateInputV19_MoveDotV16_PathBandOccluder.m`

These files are located in:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/oldies/lib`

## Historical investigation goal

The goal was to check trigger correspondence across three sources:

- Expected triggers: derived from input files and experiment code.
- Reported triggers: derived from `/experiment/output_files/debug_actual_triggers_*.csv`.
- Observed triggers: written on the trigger channel in the MEG FIF files at `/data/*.fif`.

This file preserves historical context only.
