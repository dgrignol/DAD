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
  default, site, subject, and optional tuning YAML; performs stage-specific
  validation; resolves derivative templates relative to
  project.derivatives_root; and protects the raw tree.
- dataset_inventory.py is reusable optional preflight/QC. Its JSON/Markdown
  reports are not scientific inputs consumed by later stages.
- run_pipeline.py is the orchestration entrypoint. It writes one resolved
  config snapshot and durable source manifest, runs preflight, executes each
  stage in a separate process, and records commands/status/durations.
- Stable registered stages are trigger_extraction, maxwell_filter,
  artifact_annotation, and ica_manual. The runner already supports stage
  selection, dry-run, review/apply passes, resume, overwrite,
  continue-on-error, interactive choices, and site/tuning configuration.
- 03_C_ICA_auto.py is deprecated experimental code. It is intentionally absent
  from run_pipeline.py and production requirements. Do not restore automatic
  ICA selection or use it as a fallback.
- 00_config.py and 01_inventory_and_load.py were removed. Do not recreate them
  or reintroduce numbered configuration/inventory stages.
- 02_trigger_extraction_custom.py imports pipeline_config.py directly.
- The canonical trigger products are 02_trigger_events.csv/.json,
  02_trigger_alignment.csv, 02_trigger_source_associations.json, and the
  trigger-derived associations retained in source_manifest.json.
- 03_A, 03_B, and manual 03_C are configuration-driven noninteractive CLIs.
  They preserve the complete source acquisition-run stem in derivative names
  and record derivative lineage in source_manifest.json.
- Stage 03 has explicit review and apply passes. Review creates QC evidence and
  a reusable prefit ICA model. Apply requires reviewed
  03_tuning_choices.yaml; ica.manual_exclude=null means unreviewed, while []
  means reviewed with no exclusions.
- Cluster execution uses meg/params/pipeline_cluster.yaml and meg/cluster.
  MNE 1.10.2 is required to expose the embedded Polhemus head-shape points.
  Noninteractive ICA diagnostics must use the matplotlib MNE browser backend,
  not Qt/X11.

Verified subject profiles:
- sub01_150426: stimuli V28, schedule v21, runtime v18; conditions 31/41/51;
  identities 1:20; 11 completed behavioral blocks; EOG001/EOG002 present;
  saved eye-tracker data lost, although online gaze/replay logic worked.
- sub02_130526 and sub03_200526: V29/v22/v19; conditions 31/41/51;
  identities 1:20; 10 completed blocks; no EOG channels; saved EDF eye tracking.
- sub04_220526, sub05_260526, sub06_260526, and sub07_270526:
  V30/v23/v20; conditions 71/81/91; identities 1:40; 5 completed blocks;
  no EOG channels; saved EDF eye tracking. Sub04 block02/run03 terminated after
  18 of 47 trials. The raw sub07 date 270526 is authoritative.

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
- sub02-07 also contain MISC001-MISC006 in every checked task FIF. Read-only
  signal checks show that MISC001-MISC003 carry the three analog signals sent
  from the left eye tracker, while MISC004-MISC006 are recorded but effectively
  unwired/near-flat because the second-eye cable was unavailable.
- Preserve these MISC channels and their acquisition timing. Do not treat them
  as EOG or discard them before a future eye-tracker synchronization analysis.
  The complete binocular data remain in the EDF files; the left-eye MISC traces
  can later be used to estimate EDF-to-MEG timing offset and drift.
- EDF-to-MEG synchronization is deferred to a separate future task. It should
  verify the exact MISC-to-EDF signal mapping, polarity and scale, then use
  robust windowed alignment (for example cross-correlation plus drift checks)
  rather than assuming file starts or behavioral blocks are already aligned.

Stage-03 scientific and lineage facts:
- Maxwell bad channels are detected per acquisition run. The selected policy
  is configurable as reference_run, union, or per_run; union is the default
  and per-run differences are reported. Calibration and cross-talk files are
  configured project resources.
- Maxwell filtering uses every run's static dev_head_t and maps all runs to the
  actual transform medoid: sub01 run02, sub02 run05, sub03 run03, sub04 run04,
  sub05 run02, sub06 run03, and sub07 run03.
- Maxwell reports include bad-channel detection PSDs and matched before/after
  Maxwell PSDs for every acquisition run.
- 03_B defaults to FLUX muscle threshold 10. Sub01 EOG is diagnostics-only;
  subjects 02-07 skip EOG-only functions. Trigger 150 is not handled by 03_B.
- Manual ICA fits MEG-only data at 200 Hz, 1-40 Hz, FastICA, 30 components,
  random state 97. EOG scores are diagnostic only.
- Rebuilt annotated runs share the common Maxwell destination. ICA retains
  on_mismatch="warn" only for defensive compatibility with historical
  derivatives; new production derivatives should have matching dev_head_t.
- Downstream scientific stages must consume lineage-matched manual-ICA apply
  outputs, <exact-source-stem>_ica.fif. Do not substitute review-only SSS or
  annotation files and do not reconstruct filenames by replacing block text.
- Before starting a downstream chat, verify each selected subject has a
  completed 03_C apply report, explicit ICA exclusions (possibly []), cleaned
  FIFs for every acquisition run, and completed source-manifest lineage. If
  Stage 03 is still at review_required or a cluster recovery is incomplete,
  stop and report that prerequisite instead of bypassing it.

Current baseline facts:
- The current legacy 04_epoching.py hardcodes tmin=-1.0, while
  pipeline_default.yaml currently specifies tmin=-0.5. This discrepancy must
  be resolved explicitly during Chat 4.
- 04_epoching.py uses baseline=None and linear detrending.
- sanity_checks/01_dot_side_decoding.py applies -0.4 to -0.2 s. (Should be fine in all ITI cases: both for sub01 0.5-1.0 s and sub02-07 1.0-1.5 s).
- XX_ERF_manual.py defaults to -0.1 to 0 s.

Implement the requested stage. Add tests and extensive script documentation as
required by AGENTS.md. Keep the code as readable and simple as possible.
Preserve the existing infrastructure separation: configuration,
discovery/preflight, orchestration, then numbered scientific transformations.
Do not refactor or rerun completed upstream stages unless a demonstrated
lineage or compatibility defect requires it.
```

## DONE -> Chat 1: External Dataset Profiles And Run Manifest

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

## DONE -> Chat 2: Trigger Extraction Across Experiment Versions

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

## DONE (implementation; participant review/apply must satisfy the gate above) -> Chat 3: Maxwell Filtering, Artefacts And ICA

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

## DONE (implementation; participant epoch/export execution remains separate) -> Chat 4: Epoching, Export And Sanity Checks

```text
First verify the Stage-03 completion gate from the Common Context for every
selected subject. The scientific input is the lineage-matched manual-ICA apply
output <exact-source-stem>_ica.fif. Do not consume review-only SSS/annotation
files, do not use the deprecated automatic ICA branch, and do not rerun Stage
03 unless a concrete upstream defect is found.

Refactor 04_epoching.py and 06_export_neural_for_drsa.py into documented,
configuration-driven noninteractive CLIs using pipeline_config.py, subject
profiles, the canonical trigger products, and source_manifest.json lineage.
Refactor XX_ERF_manual.py and sanity_checks/01_dot_side_decoding.py to use the
same resolved paths/event definitions; they may remain explicitly optional
inspection/sanity workflows and must not become required scientific inputs.

Replace hardcoded 31/41/51, V28, identities 1:20, sub01 paths, and
cleaned-to-raw filename replacement. Remove the manual/auto ICA variant
switch: production input is manual ICA only. Support dynamic identity counts
of 20 or 40 and discover the exact subject-specific observed trajectory file
from each external input_files folder. Use its sequence-matched nondeviant
condition as the expected continuation.

Preserve acquisition-run identity until events and metadata have been mapped
into one concatenated epoch sample frame. Match cleaned FIFs to raw acquisition
runs through source lineage, never sorted position or filename substitution.
Carry source FIF, acquisition-run label, raw sample frame, behavioral
block/run/trial when aligned, replay status, gaze-break/ESC evidence, identity,
condition, and exclusion reasons into epoch metadata and machine-readable
reports.

Audit trial exclusion logic for trigger-150 gaze breaks, replays, ESC
termination, catch triggers and both catch types, run boundaries, multi-block
FIF files, transition-resolved trigger rows, and partial runs. Use only trials
actually presented and saved. Define trial intervals within their acquisition
run so the last trial cannot extend into the next run. State which exclusions
remove epochs, which annotations are handled by reject_by_annotation, and
which checks are report-only. Never conflate trigger-150 trial rejection with
03_B/03_C physiological cleaning.

Evaluate ITI changes against every pre-stimulus window and baseline operation.
Reconcile the legacy tmin=-1.0 with the configured tmin=-0.5; do not preserve
either value accidentally. Because epoching uses baseline=None, do not claim
it performs baseline correction. The known downstream -0.4:-0.2 and -0.1:0
windows fit inside a -0.5 s epoch and are considered compatible with the
verified ITIs, but document the reasoning and test boundary behavior. Do not
silently change the scientific baseline or detrending policy.

Export identity count, model path/version, event map, trial/exclusion counts,
source lineage, ordered channel labels, time axis, sampling rate, bad/missing
channel representation, and exact inclusion policy into MAT/JSON metadata.
Do not synthesize missing identities or silently truncate 40 to 20. If an
identity has no surviving epoch in a required condition, fail explicitly or
use a separately configured and documented missing-data policy.

Keep the canonical MEG channel axis stable across subjects and conditions.
Preserve upstream bad-channel metadata; if bad/missing channels are exported
as NaN, record that policy explicitly and test it. MISC/EOG/STI channels must
not enter the dRSA neural tensor, but their omission from the export must not
modify the source FIFs.

Add tests for 20- and 40-identity subjects, both condition-code generations,
sub01 irregular runs, sub04 partial/ESC termination, replay and trigger-150
handling, catch types, run-boundary isolation, missing identities, channel-axis
consistency, raw-to-cleaned lineage matching, and guarded output paths. Add
targeted synthetic/small smoke tests rather than full participant epoching.

Register epoching and export in run_pipeline.py only after their CLIs and
dependency contracts are stable. Preserve all existing runner behavior and
keep ERF/decoding optional unless there is a clear reason to expose them as
explicit non-default stages. Update MULTISUBJECT_MIGRATION_AUDIT.md and
PIPELINE_SPEC_basic.md.
```

## DONE -> Chat 5: MATLAB dRSA And Final Orchestration

```text
First verify that Chat 4 produced lineage-valid epochs and neural exports for
the selected subjects, with explicit identity count, condition map, model
paths, channel labels, and time-axis metadata. Do not reconstruct missing
metadata inside MATLAB.

Refactor meg/code/analysis/A_dRSA.m, A_dRSA_alt.m, and dRSA_simulations.m for
the external subject-specific model files and dynamic identity counts. Invoke
MATLAB only through /Applications/MATLAB_R2020a.app/bin/matlab, using
`arch -x86_64` if MATLAB R2020a requires Intel startup on Apple Silicon.

Remove hardcoded project-root case differences, V28 filenames, and exactly-20
tensor checks. Resolve the observed model path from exported subject metadata
or another single validated handoff generated by Python. Build the expected
occluded-deviant model from the matching `occluded_nondeviant` condition in
that same MAT. Sort and validate model trials by sequence rather than MAT-file order.
Preserve strict checks that neural/model identity sets, time axes, conditions,
channel axes, channel-label order, and subject/version provenance match.

Treat 20 versus 40 identities as a scientific comparability issue. Within-
subject dRSA may use each subject's complete identity set. Do not silently
truncate sub04-07, duplicate identities, or concatenate identity axes across
subjects. If matched RDM sampling is needed, implement it only as an explicit
sensitivity-analysis option with documented subsampling and fixed seeds.

Extend the existing runner rather than redesigning it. It already supports
stable stage selection, dry-run, resume, overwrite, continue-on-error,
review/apply controls, separate subprocesses, and command/status/duration
logging. Register the stable Chat-4 epoch/export stages and an optional MATLAB
analysis handoff only after their contracts are tested. Do not add ICA variant
selection: production uses manual ICA and automatic ICA remains deprecated.

If multisubject execution is added, keep one resolved subject/config/manifest
per subprocess and per run log. A thin subjects 01-07 wrapper is preferable to
mixing subjects inside one scientific stage process. Preserve automatic
inventory preflight as QC, never as a downstream scientific input.

Keep real-neural dRSA, alternate sensor-group analyses, and simulations
optional and independently selectable. Write outputs under the configured
derivative root with subject/session/method/model-source/sensor-group metadata.
Do not overwrite existing scientific results without an explicit flag.

Run dry-run validation for all subjects, compileall, pytest, and MATLAB
syntax/smoke checks in temporary output roots that avoid full recomputation.
Test at least one 20-identity and one 40-identity fixture, shuffled MAT trial
order, missing/duplicate identities, neural/model time mismatch, channel-order
mismatch, and deterministic matched-sampling behavior if implemented.
Reconcile MULTISUBJECT_MIGRATION_AUDIT.md and PIPELINE_SPEC_basic.md with the
final implemented pipeline.
```

## Chat 6: Cluster Dot-Side Decoding Sanity Checks

```text
First verify that Chat 4 produced lineage-valid epochs for the selected
subjects and that Chat 5 did not change the Stage-04 epoching contract. This
is an optional sanity-check workflow, not a production scientific dependency.
Do not modify raw data, do not modify Stage-04 epochs, do not modify dRSA
exports, and do not register this as a required run_pipeline.py stage unless
there is a clear infrastructure reason.

Port meg/code/sanity_checks/01_dot_side_decoding.py to the DIPSEN cluster for
one-subject-per-array-task execution across subjects 01-07. Reuse the existing
cluster conventions in meg/cluster: pipeline_cluster.yaml for site paths,
headless matplotlib/MNE settings, micromamba environment activation, guarded
raw-tree output checks, cluster_logs, and explicit per-subject validation.
Create a thin SLURM wrapper, for example meg/cluster/submit_dot_side_decoding.sh,
that runs the sanity check for one subject per task and exposes condition,
decoder window, trial averaging, random seed, resume/rebuild behavior, and
resource settings through documented environment variables. Keep the wrapper
independent from Stage 04/06/07 submission; it should consume completed epochs,
not rerun epoching.

Update deployment so the minimal cluster payload includes
meg/code/sanity_checks/, the new dot-side decoding wrapper, and any helper
files required by that script. Keep deployment outside the raw data tree and do
not add rsync deletion behavior. Update meg/cluster/README.md with deployment,
submission, monitoring, retry, rebuild, and output-inspection instructions for
the dot-side sanity check.

Change the decoding algorithm from independent sample-by-sample decoding to
sliding-window decoding. Expose a tuneable sample-count width parameter, for
example --decode-window-samples, default 10 samples. Interpret this value in
the resampled decoding data grid after --decode-sfreq has been applied. For
each output time point, train/test the classifier on features from the centered
sample window around that time point, flattening channel x window samples into
one feature vector per trial or averaged pseudo-trial. At epoch boundaries, use
a documented policy such as clipping the window to available samples or
emitting only centers with complete windows; record the exact policy, window
sample counts, effective duration in milliseconds, center time, start/end
sample indices, and start/end times in the scores CSV, summary JSON, and HTML
report. Keep the default classifier behavior otherwise comparable unless a
change is justified and tested.

Add a tuneable same-identity trial averaging parameter before classification,
for example --trials-per-average, default 1. Averaging must happen only within
the selected condition and within the same trial_sequence_identity; never mix
different identities, conditions, subjects, or sessions. Use a deterministic
seeded shuffle within each identity, then form non-overlapping groups of
trials_per_average epochs and average the MEG data across each complete group.
By default, drop incomplete remainder groups so every pseudo-trial has the same
number of source trials and therefore a comparable noise reduction. Record all
dropped remainder counts by identity, the source epoch indices for every
pseudo-trial, and the random seed in machine-readable outputs. If you expose an
optional keep-remainder mode, it must be explicit, documented as lower-SNR for
the final partial group, represented in metadata, and covered by tests.

Ensure labels remain model-derived from the exact subject-specific observed
trajectory MAT file and the configured profile identity range. After trial
averaging, each pseudo-trial inherits the identity of its source group and uses
the same dynamic dot-side label time series for that identity. Validate that
the model identity set matches the subject profile, including 20-identity
subjects 01-03 and 40-identity subjects 04-07; do not truncate, duplicate, or
synthesize identities. If trial averaging removes all usable pseudo-trials for
an identity or produces too few samples for stratified decoding at any time
window, fail with a clear message or skip only that window with explicit NaNs
and a documented reason.

Preserve the existing optional ERF/report behavior, including the -0.4 to
-0.2 s ERF baseline on copied inspection epochs only. Do not claim production
epoching performs baseline correction. Keep MISC/EOG/STI channels out of the
decoder feature matrix, preserve source channel metadata, and continue to
write the three-panel PNG, scores CSV, summary JSON, and MNE HTML report under
the configured derivatives sanity-check tree.

Add focused tests for the new decoding helpers and cluster plumbing: 10-sample
default windows, even and odd window widths, edge-window policy, same-identity
trial averaging, deterministic seed behavior, dropped remainder metadata,
20- and 40-identity model label validation, too-few-class/window handling,
new wrapper bash syntax, and deployment payload inclusion. Use synthetic small
arrays/epochs rather than full participant decoding for tests. Update
MULTISUBJECT_MIGRATION_AUDIT.md and PIPELINE_SPEC_basic.md with the optional
cluster sanity-check contract and the actual verification results.

Run lightweight verification before participant jobs: python -m compileall -q
meg/code/preprocessing_basic meg/code/sanity_checks meg/tests; targeted pytest
for sanity-check helpers and cluster wrappers; bash -n for all cluster
wrappers; pipeline_config.py --validate-only with pipeline_cluster.yaml for
subjects 01-07; and at least one local or dry-run invocation that validates the
dot-side sanity-check inputs for a 20-identity and a 40-identity subject
without running all participants. Document the exact sbatch command for running
subjects 01-07 once the wrapper is ready.
```
