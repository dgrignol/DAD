# Multisubject Migration Audit

## Scope

This audit records the implemented configuration, discovery, provenance, and
orchestration foundation for the external seven-subject dataset.

- Read-only raw root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/data`
- Default derivative root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives`
- Former derivative root, retained only as migration context:
  `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives`
- Supported subject selectors: `01` through `07`

The raw tree was inspected read-only. No file or directory was created,
modified, renamed, or deleted there.

## Implemented Foundation

### Configuration

`pipeline_config.py` remains the single importable configuration API.

1. Default and optional subject YAML are recursively merged.
2. `subject.id` is normalized and restricted to `01`-`07`.
3. One profile is selected from `subject_profiles`.
4. Stage-facing trigger, epoch condition, identity, EOG, and EDF settings are
   derived from that profile.
5. Exactly one external `subNN_*` directory is required.
6. Derivative templates resolve below `project.derivatives_root`.
7. Every infrastructure output passes the shared raw-tree guard.

Subject YAML files are selectors only. Verified facts are stored once in
`pipeline_default.yaml`.

### Discovery

`dataset_inventory.py` discovers these categories independently:

- task MEG acquisition FIFs,
- date-matched empty-room FIFs,
- presented trajectory MAT,
- TrialStruct MAT,
- saved behavioral block MAT files,
- debug trigger CSV files,
- EDF files.

Task FIFs use `meg/subNN_run*.fif`. Empty-room files are searched separately at
the raw root using the date suffix from the uniquely discovered subject folder.

The inventory is optional QC. Later stages can reproduce discovery directly
from configuration.

### Durable Source Manifest

`run_pipeline.py` always writes `source_manifest.json` beside the resolved
configuration and run log.

The manifest preserves each source filename, absolute source path,
raw-root-relative path, and file size. MEG acquisition-run labels are recorded
from MEG filenames; behavioral block/run labels are recorded only from their
own files. The manifest contains an empty
`trigger_derived_associations.meg_acquisition_run_to_behavioral_segments` list
for later evidence-based mapping.

No MEG/behavior/EDF lists are zipped positionally.

## Verified Profiles

| Subject | Folder | Versions | Conditions | Identities | Saved blocks | Ocular data | Known exception |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 01 | `sub01_150426` | V28 / v21 / v18 | 31/41/51 | 1-20 | 1-11 | EOG001/EOG002; no saved EDF | Seven irregular acquisition FIFs; several contain multiple behavioral blocks |
| 02 | `sub02_130526` | V29 / v22 / v19 | 31/41/51 | 1-20 | 1-10 | EDF; MISC001-MISC006; no EOG | None |
| 03 | `sub03_200526` | V29 / v22 / v19 | 31/41/51 | 1-20 | 1-10 | EDF; MISC001-MISC006; no EOG | None |
| 04 | `sub04_220526` | V30 / v23 / v20 | 71/81/91 | 1-40 | 1-5 | EDF; MISC001-MISC006; no EOG | block02/run03 stopped after 18 of 47 trials |
| 05 | `sub05_260526` | V30 / v23 / v20 | 71/81/91 | 1-40 | 1-5 | EDF; MISC001-MISC006; no EOG | None |
| 06 | `sub06_260526` | V30 / v23 / v20 | 71/81/91 | 1-40 | 1-5 | EDF; MISC001-MISC006; no EOG | None |
| 07 | `sub07_270526` | V30 / v23 / v20 | 71/81/91 | 1-40 | 1-5 | EDF; MISC001-MISC006; no EOG | Prior migration label was `sub07_260527`; raw folder/FIF date `270526` is authoritative |

For subjects 02-07, MISC001-MISC003 are retained as left-eye analog timing
signals. MISC004-MISC006 are also retained even though currently near-flat.
EDF-to-MEG synchronization remains deferred.

## Read-Only Source Counts

| Subject | Task FIF | Presented trajectory | TrialStruct | Behavioral MAT | Trigger CSV | EDF | Date-matched empty room |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 | 7 | 1 | 1 | 11 | 33 | 0 | 0 |
| 02 | 5 | 1 | 1 | 10 | 30 | 10 | 0 |
| 03 | 5 | 1 | 1 | 10 | 30 | 10 | 1 |
| 04 | 5 | 1 | 1 | 5 | 15 | 5 | 0 |
| 05 | 5 | 1 | 1 | 5 | 15 | 5 | 1 |
| 06 | 5 | 1 | 1 | 5 | 15 | 5 | 1 |
| 07 | 5 | 1 | 1 | 5 | 15 | 5 | 1 |

These counts describe files actually present. They are not expected trial
counts. Scientific trial selection must use saved behavioral/debug data and
observed triggers, including replay, gaze-break code 150, and ESC code 152.

## Trigger Extraction Refactor

`02_trigger_extraction_custom.py` is now a documented numbered entry point over
the reusable `trigger_extraction.py` implementation.

1. Runtime version, exact debug CSV pattern, condition codes, identity range,
   required codes, optional codes, and 150/151/152 meanings come from the
   selected subject profile.
2. FIF files are ordered and recorded as MEG acquisition runs. Debug CSV files
   are parsed as unique behavioral block/run segments. Those lists are never
   zipped positionally.
3. Every non-zero STI101 window is reconstructed in its source FIF sample
   frame. Reported CSV rows guide transition resolution and alignment but never
   create a scientific event.
4. Segment alignment attaches behavioral block/run/trial, replay state, and
   reported provenance only where sequence evidence supports it. Only exact
   `MATCH` or transition-resolved `MATCH_WARNING` rows establish a durable
   acquisition-run-to-behavioral-segment association.
5. Partial and ESC-terminated segments use only rows actually saved in their
   CSV. No generated schedule or prepared block total is used.

Machine-readable outputs are:

- `02_trigger_events.csv` and `02_trigger_events.json`,
- `02_trigger_alignment.csv`,
- `02_transition_warnings.csv`,
- `02_trigger_source_associations.json`,
- the trigger stage JSON/Markdown report,
- a trigger-derived update to `source_manifest.json` when the runner supplies
  that manifest.

The legacy deep report bundle is also generated by default under
`02_triggers_extraction_custum_report_deeper`. It preserves the historical
Markdown, alignment CSV, anomaly CSV, anomaly plots, and transition-pulse
plots, but renders them from the canonical segment-aware alignment rather than
running a second comparator with legacy block-encoded FIF assumptions.

The event output preserves source FIF, acquisition run, raw sample frame,
behavioral block/run/trial when aligned, replay status, gaze-break code 150,
and ESC code 152.

### Verified v18 Evidence

The retained read-only sub01 extraction artifacts generated on April 24, 2026
provide the available full-dataset v18 reference:

- 7 acquisition FIFs,
- 33 CSV segments = 11 behavioral blocks x 3 within-block runs,
- 3,626 reported rows,
- 3,637 reconstructed raw trigger windows,
- 3,624 resolved observed events,
- 13 transition pairs, all resolved from reported trigger evidence,
- alignment: 3,611 `MATCH`, 13 `MATCH_WARNING`, 2 `REPORTED_ONLY`,
  0 `MISMATCH`, and 0 `OBSERVED_ONLY`.

The two reported-only rows are the first condition/identity pair in
block01/run01: codes 31 and 9. The first observed event corresponds to the next
presented trial. This remains an explicit alignment exception rather than a
synthesized MEG event.

Exact observed v18 code counts:

```text
1:36, 2:41, 3:34, 4:37, 5:42, 6:38, 7:40, 8:34, 9:37, 10:35,
11:35, 12:36, 13:38, 14:39, 15:35, 16:37, 17:35, 18:36, 19:35,
20:36, 31:242, 41:253, 51:244, 111:488, 112:488, 114:487,
115:487, 116:22, 117:22, 118:65, 150:14, 151:11, 201:60, 202:5
```

Codes 113 and 152 were absent in the sub01 reference, so they are optional
profile codes rather than universal required-code failures.

### Verified v19 Evidence

The complete sub02 v19 dataset was read from the mounted external volume on
June 10, 2026, with outputs redirected to `/private/tmp`:

- 5 acquisition FIFs,
- 30 CSV segments = 10 behavioral blocks x 3 within-block runs,
- 3,637 reported rows,
- 3,647 reconstructed raw trigger windows,
- 3,637 resolved observed events,
- 10 transition pairs, all resolved from reported trigger evidence,
- alignment: 3,627 `MATCH`, 10 `MATCH_WARNING`, and no mismatch,
  reported-only, or observed-only rows,
- 49 gaze-break events, 22 replay-start events, and no ESC event.

Exact observed v19 code counts:

```text
1:39, 2:39, 3:36, 4:38, 5:35, 6:34, 7:41, 8:37, 9:35, 10:35,
11:41, 12:35, 13:37, 14:39, 15:36, 16:34, 17:32, 18:37, 19:35,
20:35, 31:244, 41:249, 51:245, 111:468, 112:468, 114:468,
115:468, 116:30, 117:30, 118:83, 150:49, 151:22, 201:74, 202:9
```

### Verified v20 Evidence

The complete sub04 v20 dataset was read from the mounted external volume on
June 10, 2026, with outputs redirected to `/private/tmp`:

- 5 acquisition FIFs,
- 15 CSV segments = 5 behavioral blocks x 3 within-block runs,
- 3,305 reported rows,
- 3,319 reconstructed raw trigger windows,
- 3,306 resolved observed events,
- 13 transition pairs, all resolved from reported trigger evidence,
- alignment: 3,292 `MATCH`, 13 `MATCH_WARNING`, 1 `OBSERVED_ONLY`,
  and no mismatch or reported-only rows,
- 26 gaze-break events, 5 replay-start events, and 1 ESC event.

Exact observed v20 code counts:

```text
1:15, 2:16, 3:17, 4:14, 5:16, 6:19, 7:21, 8:18, 9:16, 10:16,
11:16, 12:18, 13:15, 14:17, 15:17, 16:18, 17:17, 18:16, 19:17,
20:14, 21:14, 22:17, 23:18, 24:19, 25:20, 26:17, 27:16, 28:17,
29:16, 30:17, 31:18, 32:16, 33:19, 34:18, 35:16, 36:17, 37:17,
38:17, 39:14, 40:15, 71:238, 81:219, 91:222, 111:426, 112:426,
114:426, 115:426, 116:30, 117:30, 118:80, 150:26, 151:5,
152:1, 201:70, 202:10
```

The sole unresolved CSV alignment case is the known partial sub04
block02/run03:

- its CSV contains 85 actual saved rows and ends at trial label 19 after
  trials 1-18 were completed,
- all 85 saved rows match acquisition run 2,
- STI101 then records ESC 152 at raw sample frame 1,058,774,
- no matching 152 row was saved in the CSV,
- the event remains `OBSERVED_ONLY` but is assigned behavioral block 2/run 3
  from the preceding exact match in the same FIF,
- the segment termination status is `observed_esc_after_saved_rows`.

This is retained as an explicit acquisition-vs-runtime evidence difference;
the ESC event is not dropped or synthesized.

Synthetic tests additionally cover symmetric transition orders, identity 40,
replay 151, gaze break 150, CSV-saved ESC 152, observed-only ESC 152,
multi-block FIFs, partial segments, and exact-match-only source associations.

### All-Subject Production Run

The registered trigger stage was run through `run_pipeline.py` for all seven
subjects on June 10, 2026. Outputs were written under the configured external
derivative root. Each run completed with `inventory_status=ok`, return code 0,
a resolved trigger-derived source-manifest association, and `aligned` status
for every acquisition FIF.

| Subject | Runtime | FIF | CSV segments | Reported rows | Observed events | MATCH | MATCH_WARNING | REPORTED_ONLY | OBSERVED_ONLY | Deep anomaly plots | Deep transition plots |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 | v18 | 7 | 33 | 3,626 | 3,624 | 3,611 | 13 | 2 | 0 | 2 | 13 |
| 02 | v19 | 5 | 30 | 3,637 | 3,637 | 3,627 | 10 | 0 | 0 | 0 | 10 |
| 03 | v19 | 5 | 30 | 3,515 | 3,515 | 3,496 | 19 | 0 | 0 | 0 | 19 |
| 04 | v20 | 5 | 15 | 3,305 | 3,306 | 3,292 | 13 | 0 | 1 | 1 | 13 |
| 05 | v20 | 5 | 15 | 3,514 | 3,514 | 3,497 | 17 | 0 | 0 | 0 | 17 |
| 06 | v20 | 5 | 15 | 3,540 | 3,540 | 3,526 | 14 | 0 | 0 | 0 | 14 |
| 07 | v20 | 5 | 15 | 3,584 | 3,584 | 3,562 | 22 | 0 | 0 | 0 | 22 |

All 108 transition-warning plots and all 3 unmatched-row anomaly plots were
generated without plotting errors. The legacy deep Markdown for each subject
is located at:

```text
/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-XX/ses-01/02_triggers_extraction_custum_report_deeper/subXX_all_blocks_reported_vs_observed.md
```

The unmatched alignment cases are limited to the two known sub01
`REPORTED_ONLY` rows at the start of block01/run01 and the known sub04
`OBSERVED_ONLY` ESC 152 after the saved rows of block02/run03. Every
`MATCH_WARNING` row is a resolved short STI101 transition backed by the saved
reported sequence; none remains an unresolved transition.

## Safety and Failure Behavior

- Zero `subNN_*` folder matches fail.
- Multiple `subNN_*` folder matches fail as ambiguous.
- Missing task FIF, presented trajectory, TrialStruct, behavioral MAT, or
  trigger CSV inputs fail preflight.
- Missing EDF fails when the selected profile expects saved EDF data.
- Missing EDF for sub01 is valid and records the eye-model skip reason.
- Any output inside the configured raw root is rejected before writing.
- Empty-room recordings cannot enter the task-run list.

## Verification

The migration foundation is covered by tests for:

- all seven subject-folder matches,
- sub01 irregular run filenames,
- sub07 empty-room exclusion,
- sub07 date-label discrepancy,
- zero and ambiguous folders,
- missing required inputs,
- raw-tree output rejection,
- profile-derived condition and identity settings,
- manifest separation and pending trigger-derived associations.

Verification completed on June 10, 2026:

- `python -m compileall -q meg/code/preprocessing_basic meg/tests`: passed.
- `python -m pytest meg/tests -q`: `38 passed`.
- `pipeline_config.py --validate-only`: passed for sub01-sub07.
- `run_pipeline.py --dry-run`: completed with `inventory=ok` for sub01-sub07.
- `git diff --check`: passed.

## Stage 07 MATLAB dRSA Refactor

The MATLAB dRSA layer is now separate from Stage 06. Stage 06 writes only the
identity-averaged neural tensor contract; Stage 07 consumes those exports and
creates MATLAB dRSA results, PNG figures, and figure-bearing reports.

### Handoff and orchestration

`07_run_matlab_drsa.py` is the optional Python orchestration CLI for MATLAB
dRSA. It validates completed `export_drsa` lineage, reads every Stage-06
condition sidecar, discovers the subject-specific observed trajectory MAT
through the resolved profile, records the nondeviant expected-source
condition, then writes:

- `07_matlab_drsa_handoff.json`
- `07_matlab_drsa_report.json`
- `07_matlab_drsa_report.html`

under `derivatives/reports/sub-XX/ses-01/`.

The handoff records the subject/session, identity axis, conditions, time axis,
sampling rate, canonical channel labels, Stage-06 MAT/JSON paths, model paths,
and derivative output directories. MATLAB scripts consume that handoff instead
of reconstructing paths or assuming V28/20 identities.

`run_pipeline.py` registers `matlab_drsa` as an optional downstream stage. It
inherits resume/overwrite/dry-run behavior and remains independent from
epoching and neural export. `meg/cluster/submit_drsa_analysis.sh` submits one
subject per SLURM array task. `meg/cluster/deploy_review_to_cluster.sh` now
deploys `meg/code/analysis/`, `07_run_matlab_drsa.py`, the Stage-07 wrapper,
and `simulations/functions/`.

### MATLAB scripts

The refactored MATLAB layer contains:

- `drsa_stage_utils.m`: shared handoff parsing, model loading, sequence
  sorting, dynamic identity validation, sensor-group selection, NaN-channel
  handling, and headless PNG export.
- `A_dRSA.m`: main real-neural dRSA for all exported conditions. It preserves
  20- and 40-identity subjects and writes condition PNG/MAT outputs under
  `analysis/sub-XX/ses-01/drsa/matlab_real_neural/07_dRSA/`.
- `A_dRSA_alt.m`: optional single-model/sensor-group real-neural dRSA.
- `dRSA_simulations.m`: optional model-vs-model and neural-vs-neural
  simulation-style outputs.

All figures are exported headlessly with `Visible='off'` and `exportgraphics`.
The subject HTML report embeds all discovered Stage-07 PNGs and links the
machine-readable manifest.

### Stage-7 verification

Final Stage-07 verification on June 16, 2026:

- `python -m py_compile` for the new Stage-07 CLI and touched runner/config
  modules: passed.
- `bash -n` for deployment and Stage-07 cluster wrappers: passed.
- MATLAB R2020a `checkcode` through `arch -x86_64` for all four analysis files:
  zero messages.
- `python -m pytest meg/tests -q`: `94 passed`.
- `07_run_matlab_drsa.py --validate-only`: passed for subjects 01-07 against
  the completed Stage-06 derivatives.
- Representative full trigger extraction passed for v18 sub01, v19 sub02,
  and v20 sub04 with temporary derivatives before the production run.
- Production trigger extraction passed through `run_pipeline.py` for
  sub01-sub07 under
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives`.

The dry-runs wrote resolved configuration, source manifest, and run log files
only under the configured derivative root. Expensive preprocessing was not
run during those dry-runs; the later production execution ran the registered
trigger stage only.

## Stage 03 Configuration Refactor

On June 11, 2026, 03_A, 03_B, and the production manual 03_C branch were
refactored into import-safe configuration-driven CLIs.

### Shared behavior

- Optional tuning YAML is merged after the subject selector.
- Every derivative preserves the complete raw source stem and acquisition-run
  label.
- FIF headers, JSON reports, and `source_manifest.json` record source lineage.
- Manifest refreshes preserve trigger-derived associations and derivative
  lineage.
- Existing scientific derivatives require explicit overwrite permission.
- Interrupted production 03 stages support guarded resume; existing outputs
  are skipped only after stage-specific metadata validation.
- Interactive choices are written to guarded `03_tuning_choices.yaml` for a
  later noninteractive apply pass.
- Machine paths can be supplied through a distinct site-config layer. The
  merge order is default, site, subject, then tuning.
- `meg/cluster/submit_preprocessing_review.sh` runs one noninteractive review
  job per subject as a SLURM array and never invokes automatic ICA.
- `render_stage03_mne_reports.py` converts completed stage JSON and PNG
  artifacts into portable, self-contained MNE reports without reopening FIF
  files or rerunning scientific processing. When available,
  `03_tuning_choices.yaml` is embedded so report evidence and accepted
  decisions remain visible together.

### 03_A Maxwell filtering

- Bad channels are detected independently for every acquisition run.
- The default `union` policy applies every channel detected in any run to all
  runs. `reference_run`, `per_run`, and global/per-run manual overrides remain
  available.
- Every run uses its static `dev_head_t`. A common actual-run destination is
  selected as the transform medoid over the subject's digitization cloud:
  sub01 run02, sub02 run05, sub03 run03, sub04 run04, sub05 run02, sub06
  run03, and sub07 run03.
- Reports include matched raw-versus-Maxwell PSD plots for every acquisition
  run in addition to bad-channel detection PSDs.
- Cross-talk and calibration files are explicit project resources rather than
  subject raw-data assumptions.
- Resume also validates the destination run and destination-transform
  fingerprint; native-position historical derivatives cannot be reused.

### 03_B physiological artifacts

- The noninteractive default retains the FLUX muscle threshold of 10.
- Interactive mode replots each run's muscle z-score and threshold until the
  operator accepts a value.
- Sub01 EOG defaults to diagnostics-only.
- Subjects 02-07 skip EOG-only functions because no EOG channels exist.
- MISC001-MISC006 are preserved and never treated as EOG.
- Trigger 150 remains outside physiological artifact annotation.
- Resume validates lineage, threshold, channels, bad-channel metadata, and
  FIF/CSV annotation-count consistency.

### 03_C manual ICA

- Production fitting retains FastICA, 30 components, seed 97, 200 Hz, and
  1-40 Hz fit preprocessing.
- Prefit cache reuse requires a matching mount-independent input/parameter
  fingerprint; legacy absolute-path metadata is normalized by filename.
- Review reports include source traces, topographies, component properties,
  and EOG diagnostics where available.
- Apply mode requires an explicit exclusion list; null means unreviewed and an
  empty list means reviewed with no rejection.
- The hardcoded ICA000 blink workaround was removed.
- Apply resume validates the post-rejection model, exclusions, lineage,
  channels, bads, samples, sampling frequency, and preserved MISC channels.

### Automatic ICA status

`03_C_ICA_auto.py` is retained only as deprecated experimental research code.
It exits before loading MEGNet dependencies during normal execution, is absent
from `run_pipeline.py`, and is not used for participant processing. Optional
dependencies are isolated in `requirements-ica-auto-experimental.txt`.

### Verification additions

Tests cover acquisition-run naming, bad-channel policy selection, tuning-file
merges, raw-tree rejection, modality-aware EOG handling, MISC/bad metadata
preservation, manual exclusion parsing, automatic-ICA deprecation, runner
command construction, and read-only headers for all 37 task FIFs.

## Stage 04 Epoching And dRSA Export Refactor

On June 15, 2026, the Stage-03 completion gate was rechecked for subjects
01-07. Every subject had:

- a completed `03_C_ICA_apply_report.json`;
- explicit manual exclusions;
- `ica_manual` apply/completed source-manifest lineage;
- one `<exact-source-stem>_ica.fif` for every acquisition run.

The manifest paths were produced on the cluster under `/mnt/storage`, while
the same files are locally mounted under `/Volumes/MORWUR`. Downstream
resolution therefore validates the recorded basename/source stem and resolves
the local file from the configured derivative root rather than trusting or
rewriting an absolute mount prefix.

### Shared downstream contracts

`downstream_stage_utils.py` now provides:

- the manual-ICA apply completion gate;
- raw-to-cleaned matching by exact source stem;
- exact observed trajectory discovery and sequence-matched expected-path selection;
- canonical trigger CSV validation;
- per-run trial intervals and identity matching;
- explicit gaze, catch, ESC, replay, transition, and exclusion metadata;
- missing-identity policy validation.

### Epoching

`04_epoching.py` is a configuration-driven noninteractive CLI. It uses only
manual-ICA apply outputs and the canonical trigger event table. Runs are
epoched independently and concatenated after metadata mapping.

The legacy `tmin=-1.0` discrepancy was resolved in favor of the reviewed
configured value `-0.5 s`. Production epoching retains `tmax=2.667 s`,
`baseline=None`, linear detrending, amplitude rejection, and
`reject_by_annotation`.

Default trial removals are trigger 150, both catch families, ESC termination,
and unresolved transition evidence. Replay trials remain included and
labelled. The final trial interval is bounded by its acquisition-run end.

Outputs include epochs FIF, surviving metadata CSV, complete trial-decision
CSV, JSON/HTML reports, and epoching source-manifest lineage.

### dRSA export

`06_export_neural_for_drsa.py` is a configuration-driven noninteractive CLI.
It supports 20- and 40-identity profiles and writes one MAT tensor plus JSON
sidecar per condition.

The channel axis is the fixed repository-owned MEG list. Upstream bad channels
and missing channels are represented as NaN and documented. MISC, EOG, STI,
and other non-MEG channels are omitted without modifying source data. Missing
required identity cells fail by default; the optional `nan` policy must be
configured explicitly.

MAT/JSON metadata includes model paths/version, event definitions, identity
and trial counts, inclusion policy, source lineage, channel labels, time axis,
sampling rate, and bad/missing channel representation.

### Optional workflows and runner

`XX_ERF_manual.py` and `sanity_checks/01_dot_side_decoding.py` now use resolved
epoch paths, subject event definitions, profile identity ranges, and exact
subject-specific trajectory models. They remain optional and are not
registered as required stages.

`run_pipeline.py` registers stable `epoching` and `export_drsa` stages while
preserving preprocessing review/apply, dry-run, resume, overwrite, and
continue-on-error behavior.

### Dot-side decoding cluster sanity check

`sanity_checks/01_dot_side_decoding.py` now validates completed Stage-04
epoching lineage before consuming epochs, remains outside `run_pipeline.py`,
and writes only optional sanity-check derivatives. The decoder uses
model-derived dynamic dot-side labels from the exact subject-specific
trajectory MAT file and validates the configured profile identity range, so
20-identity subjects 01-03 and 40-identity subjects 04-07 are handled without
truncation, duplication, or synthesized identities.

The optional workflow now has four explicit analysis identities:
`always_visible_observed`, `occluded_nondeviant_observed`,
`occluded_deviant_observed`, and `occluded_deviant_predicted`. Predicted
deviant decoding keeps `occluded_deviant` neural epochs while using the
matching `occluded_nondeviant` trajectory as the expected dot-side label path.
Condition, model source, model-label condition, and analysis key are recorded,
and output paths include both condition and source to prevent collisions.

Decoding changed from independent sample-by-sample features to centered
sliding windows on the resampled decoding grid. The default window is
10 samples after `--decode-sfreq`; odd widths are symmetric, even widths place
the extra sample after the center, and epoch-edge windows are clipped to
available samples. The scores CSV and summary JSON record requested/effective
window sample counts, duration in milliseconds, center time, start/end sample
indices, start/end times, edge clipping, class counts, and explicit skip
reasons.

`--trials-per-average` adds same-identity pseudo-trial averaging before
classification. Averaging never mixes identities, conditions, subjects, or
sessions. Numeric values use complete non-overlapping groups after seeded
within-identity shuffling; incomplete remainder groups are dropped and recorded
by identity, with source epoch indices stored for every pseudo-trial. The
literal value `all_trials_per_identity` averages every retained source trial
for an identity into one identity-level pseudo-trial. An identity with no
complete group is dropped and explicitly recorded without aborting usable
identities; the run fails only if every identity is underfilled. Subject and
group reports include training/testing diagnostics for classifier type, CV
folds, averaging mode, pseudo-trial counts, class balance, skip reasons, and
whether decoding is fully or partially well-defined. At group level, time
points with any subject-level NaN are excluded from the statistic and left as
plot gaps so cluster tests do not bridge skipped windows.

`meg/cluster/submit_dot_side_decoding.sh` runs one subject per SLURM array task
and executes all four standard analyses sequentially by default. It reuses the
cluster conventions for `pipeline_cluster.yaml`, micromamba activation,
headless matplotlib/MNE settings, `cluster_logs`, validation, and resume/rebuild
behavior. Deployment includes `meg/code/sanity_checks/` and the dot-side
wrapper without adding rsync deletion behavior.

`sanity_checks/02_dot_side_group_report.py` and
`meg/cluster/submit_dot_side_group_report.sh` build the optional group report
after subject-level sanity-check outputs exist. The report validates matched
time grids and decoding settings across subjects, uses balanced accuracy as
the default metric, tests accuracy minus 0.5 with a one-sample sign-flip
cluster permutation over the time axis, and annotates the group plot with the
cluster-forming p threshold and equivalent t cutoff. Display bands default to
95% confidence intervals; SEM is available explicitly with `--interval sem`.
The combined HTML report emits all four group results and cluster tables first,
then subjects 01-07. Occluded participant and group decoding plots mark the
configured first complete occlusion and first reappearance times. Each subject
section also begins with source-trial counts by identity for all three neural
conditions before its four decoding figures; duplicate observed/predicted
deviant summaries are checked for count consistency rather than plotted twice.

### Chat-4 verification

Focused tests cover both condition-code generations, 20/40 identities, sub01
irregular run lineage, sub04 partial/ESC termination, trigger 150, replay,
both catch types, resolved and unresolved transitions, acquisition-run
boundaries, missing identities, fixed channel order, bad-channel NaNs,
non-MEG omission, raw-tree guards, and manifest updates.

A small synthetic cleaned FIF completed epoching and dRSA MAT export. No full
participant epoching or export was run as part of this refactor.

Final verification on June 15, 2026:

- `python -m compileall -q meg/code/preprocessing_basic meg/code/sanity_checks meg/tests`: passed.
- `python -m pytest meg/tests -q`: `85 passed`.
- `pipeline_config.py --validate-only`: passed for subjects 01-07.
- `04_epoching.py --validate-only`: passed the real Stage-03 lineage/cleaned-file gate for subjects 01-07.
- `run_pipeline.py --stages epoching export_drsa --dry-run`: completed with `inventory_status=ok` for subjects 01-07.
- Real observed trajectory label construction produced 20 identities for sub01 V28 and 40 identities for sub04 V30.
- `git diff --check`: passed.

### Chat-6 verification

Final verification on June 16, 2026:

- `python -m compileall -q meg/code/preprocessing_basic meg/code/sanity_checks meg/tests`: passed.
- `python -m pytest meg/tests/test_dot_side_decoding.py meg/tests/test_dot_side_group_report.py meg/tests/test_cluster_wrappers.py -q`: `15 passed`.
- `bash -n meg/cluster/*.sh`: passed for all cluster wrappers.
- Dot-side `--validate-only --print-summary`: passed for sub01
  `always_visible` as the 20-identity profile and sub04 `always_visible` as
  the 40-identity profile.
- Local `pipeline_config.py --validate-only` with default paths passed for
  subjects 01-07.
- Cluster-site `pipeline_config.py --site-config
  meg/params/pipeline_cluster.yaml --validate-only` could not complete on the
  laptop because `/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data`
  is not mounted locally. Run this validation on DIPSEN before
  participant array submission.
- A full local sub05 `always_visible` run completed with 201 source trials,
  201 decoder trials, 320 decoded windows, 10-sample sliding windows, no
  skipped windows, peak unbalanced accuracy 0.7861 at 0.1917 s, and peak
  balanced accuracy 0.7000 at 0.1750 s. Outputs were written under
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/sub-05/ses-01/dot_side_decoding/always_visible/`.
- Local subject-level dot-side outputs were generated or resumed for subjects
  01-07. The local group report completed for `always_visible` with balanced
  accuracy, 95% CI, cluster-forming `p=0.05`, exact sign flips (`n=128`), and
  significant above-chance clusters at 0.000-0.483 s (p=0.0078),
  0.533-0.617 s (p=0.0313), 0.875-0.992 s (p=0.0313), and
  1.217-1.283 s (p=0.0313).

Four-analysis extension verification on June 22, 2026:

- Focused dot-side and wrapper tests: `30 passed`.
- Complete MEG test suite: `140 passed` with six existing synthetic-FIF naming
  warnings.
- Compileall, dot-side wrapper `bash -n`, and targeted `git diff --check`:
  passed.
- Read-only real-input validation passed all four analyses for sub01 and sub04,
  covering 20- and 40-identity profiles and predicted deviant labels.
- Synthetic combined-report generation verified four group analyses before all
  participant-specific sections.
- Existing participant decoding derivatives were not regenerated during this
  implementation pass.
