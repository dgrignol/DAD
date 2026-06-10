# Multisubject MEG Refactor Prompts

Use the chats in order. Paste the Common Context followed by the prompt for
the selected stage.

## Common Context

```text
Work in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD and follow AGENTS.md.
Do not revert unrelated working-tree changes.

Raw data root:
  /Volumes/MORWUR/Projects/DAMIANO/Attention/data

The complete raw-data tree is read-only. Never create, modify, rename, or
delete anything there. Every output must pass the shared raw-tree guard in:
  meg/code/preprocessing_basic/pipeline_config.py

Use this derivative root unless a configuration override is supplied:
  /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives
Former derivative root, retained only for migration context:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives

Infrastructure already implemented:
- pipeline_config.py is the single importable configuration API. It merges
  YAML, performs stage-specific validation, resolves derivative templates
  relative to project.derivatives_root, and protects the raw tree.
- dataset_inventory.py is reusable optional preflight/QC. Its JSON/Markdown
  reports are not scientific inputs consumed by later stages.
- run_pipeline.py is the orchestration entrypoint. It writes one resolved
  config snapshot, runs preflight, logs subprocess execution, and currently
  registers only trigger_extraction.
- 00_config.py and 01_inventory_and_load.py were removed. Do not recreate them
  or reintroduce numbered configuration/inventory stages.
- 02_trigger_extraction_custom.py imports pipeline_config.py directly.

Verified subject profiles:
- sub01_150426: stimuli V28, schedule v21, runtime v18; conditions 31/41/51;
  identities 1:20; 11 completed behavioral blocks; EOG001/EOG002 present;
  saved eye-tracker data lost, although online gaze/replay logic worked.
- sub02_130526 and sub03_200526: V29/v22/v19; conditions 31/41/51;
  identities 1:20; 10 completed blocks; no EOG channels; saved EDF eye tracking.
- sub04_220526 through sub07_260527: V30/v23/v20; conditions 71/81/91;
  identities 1:40; 5 completed blocks; no EOG channels; saved EDF eye tracking.
  Sub04 block02/run03 terminated after 18 of 47 trials.

Use trials actually presented and saved, including replay and ESC information.
Do not use generated/prepared trial or block counts as expected data counts.

MEG is organized by acquisition run, while behavioral output and EDF files are
organized by experiment block. Never zip those lists positionally.
- Sub01 has seven FIFs, irregular run numbering, and some FIFs contain multiple
  behavioral blocks.
- Subjects 2-7 have five task FIFs.

Always load subject-specific trajectory files from each external input_files
folder. Never hardcode V28 or assume exactly 20 identities.

ITI changed: sub01 0.5-1.0 s; sub02-07 1.0-1.5 s, with visual-content
differences. Trigger 150 records online gaze breaks independently of whether an
EDF file was saved.

Eye/ocular data:
- sub01: EOG exists in FIF; saved eye-tracker data is unavailable.
- sub02-07: no EOG channels in FIF; EDF files exist per behavioral block.

Current baseline facts:
- 04_epoching.py uses baseline=None and linear detrending.
- sanity_checks/01_dot_side_decoding.py applies -0.4 to -0.2 s.
- XX_ERF_manual.py defaults to -0.1 to 0 s.

Implement the requested stage, not only an audit. Add tests and extensive
script documentation as required by AGENTS.md. Keep the existing infrastructure
separation: configuration, discovery/preflight, orchestration, then numbered
scientific transformations.
```

## Chat 1: External Dataset Profiles And Run Manifest

```text
Extend the existing pipeline_config.py, dataset_inventory.py, run_pipeline.py,
and YAML configuration for the external multisubject dataset. Do not redesign
or replace these modules.

Make subject selection accept 01-07 and safely discover exactly one external
folder matching sub<NN>_*. Fail on zero or ambiguous matches. Discover:
- task FIFs under meg/sub<NN>_run*.fif;
- input_files trajectory and TrialStruct files;
- output_files behavioral MAT and debug trigger CSV files;
- eyetracker EDF files;
- optional empty-room recordings separately from task runs.

Create subject profiles for all seven participants containing folder/version
information, condition triggers, identity range, ocular modalities, completed
behavioral blocks, and known exceptions. Store facts once and derive stage
configuration from them.

Create a durable source manifest that preserves acquisition-run filenames and
source paths without claiming that MEG run equals behavioral block. Include a
place for later trigger-derived block/run associations. Do not make the
optional inventory report a required downstream input; later stages should be
able to reproduce discovery from config, while saved manifests provide
provenance and lineage.

Change project.data_root to the external read-only root and keep derivatives
under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives. Add tests for all
subject-folder matches, sub01 irregular run names, sub07 empty-room exclusion,
sub07 date mismatch, ambiguous folders, missing inputs, and raw-tree output
rejection.

Update MULTISUBJECT_MIGRATION_AUDIT.md and PIPELINE_SPEC_basic.md. Run
compileall, pytest, config validation, and run_pipeline.py --dry-run for
sub01-sub07. Do not run expensive preprocessing.
```

## Chat 2: Trigger Extraction Across Experiment Versions

```text
Read the migration audit and use the profiles/manifest foundation from Chat 1.
Refactor 02_trigger_extraction_custom.py for all seven participants.

Condition codes, identity range, runtime version, expected CSV pattern, and
required trigger policy must come from the subject profile:
- v18/v19: 31/41/51 and identities 1:20;
- v20: 71/81/91 and identities 1:40.

Audit how debug_actual_triggers CSV files relate to behavioral blocks, replays,
and MEG acquisition runs. Do not concatenate every CSV blindly. Observed STI101
events are primary evidence; reported CSVs are alignment/QC evidence.

Preserve source FIF, acquisition run, behavioral block when determinable,
replay status, ESC events, and raw sample frame in machine-readable outputs.
Use actual saved/presented trials rather than prepared schedule totals.

Add synthetic parser tests for both trigger generations, identity 40, replay
151, gaze break 150, ESC 152, transition pulses, multi-block FIF files, and a
partial final block. Register/update the trigger stage in run_pipeline.py
without making inventory JSON a required input.

Run targeted read-only validation on representative v18, v19, and v20 FIFs and
debug CSVs. Update the migration audit with exact verified counts and any
unresolved alignment cases.
```

## Chat 3: Maxwell Filtering, Artefacts And ICA

```text
Read the migration audit and refactor 03_A_MaxFilt.py,
03_B_ArtefactAnnotation.py, 03_C_ICA.py, and 03_C_ICA_auto.py into
configuration-driven noninteractive CLI stages using pipeline_config.py and the
shared source-lineage conventions.

Remove hardcoded sub01 paths, block-based filename substitutions,
unconditional overwrite behavior, and interactive debug code from normal
execution. Preserve the source acquisition-run identity through every
derivative filename and metadata record.

Evaluate the current policy of deriving bad channels from the first run and
applying them to all runs. Make the selected policy explicit/configurable and
report per-run differences. Keep Maxwell calibration/crosstalk resources
configurable and outside subject raw-data assumptions.

Implement ocular handling by available modality:
- sub01 has EOG but no saved eye-tracker file;
- sub02-07 have no EOG and must not call EOG-only functions;
- EDF files do not map one-to-one to MEG runs;
- trigger-150 trial rejection is separate from physiological artefact removal.

Make manual ICA exclusions subject-configured and noninteractive in batch mode.
EOG diagnostics must skip cleanly when channels are absent. MEGNet dependency
or inference failures must have explicit QC status and must not silently switch
methods.

Register these stages in run_pipeline.py only after their CLIs are stable.
Add lightweight header/synthetic tests and targeted smoke tests; avoid full
expensive processing unless required. Update documentation and the audit.
```

## Chat 4: Epoching, Export And Sanity Checks

```text
Refactor 04_epoching.py, 06_export_neural_for_drsa.py, XX_ERF_manual.py, and
sanity_checks/01_dot_side_decoding.py using pipeline_config.py, subject
profiles, and source lineage.

Replace hardcoded 31/41/51, V28, identities 1:20, sub01 paths, and
cleaned-to-raw filename replacement. Support dynamic identity counts of 20 or
40 and discover the exact subject-specific trajectory model.

Audit trial exclusion logic for trigger-150 gaze breaks, replays, ESC
termination, catch triggers and both catch types, run boundaries, multi-block
FIF files, and partial runs. State which exclusions change neural data and
which are report-only.

Evaluate ITI changes against every pre-stimulus window and baseline operation.
Because epoching currently uses baseline=None, do not claim it performs
baseline correction. Decide and document whether the -1 s epoch start and
downstream -0.4:-0.2 or -0.1:0 baselines are valid for sub01's 0.5 s minimum
ITI. Do not silently change the scientific baseline policy.

Export identity count, model path/version, event map, trial/exclusion counts,
and provenance into MAT/JSON metadata. Add tests for 20- and 40-identity
exports, missing identities, and lineage matching. Register stable stages in
run_pipeline.py and update the audit.
```

## Chat 5: MATLAB dRSA And Final Orchestration

```text
Refactor A_dRSA.m, A_dRSA_alt.m, and dRSA_simulations.m for external
subject-specific model files and dynamic identity counts. Invoke MATLAB only
through /Applications/MATLAB_R2020a.app/bin/matlab.

Remove V28 and exactly-20 assumptions. Sort and validate model trials by
sequence rather than MAT-file order. Preserve strict checks that neural/model
identity sets, time axes, conditions, and channel axes match.

Treat 20 versus 40 identities as a scientific comparability issue. Within-
subject dRSA may use each subject's complete identity set. Do not silently
truncate sub04-07, duplicate identities, or concatenate identity axes across
subjects. If matched RDM sampling is needed, implement it only as an explicit
sensitivity-analysis option with documented subsampling and fixed seeds.

Complete run_pipeline.py so it can select subjects 01-07 and stages, and
supports dry-run, resume, overwrite, continue-on-error, and ICA variant
selection. Invoke each stage in a separate process, retain commands/status/
durations, and keep MATLAB analysis optional. The inventory remains automatic
preflight/QC, not a scientific dependency.

Run dry-run validation for all subjects, compileall, pytest, and MATLAB
syntax/smoke checks that avoid full recomputation. Reconcile
MULTISUBJECT_MIGRATION_AUDIT.md and PIPELINE_SPEC_basic.md with the final
implemented pipeline.
```
