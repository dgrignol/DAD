# Project Context

This project contains MEG data from an experiment, together with the related input files and experiment outputs.

## Canonical Experiment Files
The experiment files to use as reference are:

- `stimuli_generation_V26_PathBandOccluder.m`
- `CreateInputFiles_v19_threeRunsPerBlock_catch_PathBandOccluder.m`
- `MoveDot1_experiment_occlusion_v16_PathBandOccluder.m`
- `trigger_codes_occlusion_v7.md`

These files are located in:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/oldies`

## Canonical Config Files
The config files to use as reference are:

- `Config_stimuli_generation_V26_PathBandOccluder.m`
- `Config_schedule_CreateInputV19_MoveDotV16_PathBandOccluder.m`

These files are located in:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/oldies/lib`

## Project Goal
The goal of this project is to check trigger correspondence across three sources:

- Expected triggers: derived from input files and experiment code.
- Reported triggers: derived from `/experiment/output_files/debug_actual_triggers_*.csv`.
- Observed triggers: written on the trigger channel in the MEG FIF files at `/data/*.fif`.

For now, this file defines the context only.
