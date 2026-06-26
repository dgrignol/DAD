# MEG Preprocessing Basic Pipeline Specification

## Scope

This document describes the implemented external multisubject infrastructure
and production stages through dRSA neural export.

- Code root:
  `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic`
- Read-only raw root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/data`
- Default derivative root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives`
- Subjects: `sub-01` through `sub-07`
- Session: `ses-01`

The raw tree is never an output location. Every infrastructure and scientific
stage output passes `pipeline_config.assert_output_path_outside_raw()` before
directory creation or writing.

## Implemented Stages

| Order | Script | Production status | Role |
| --- | --- | --- | --- |
| Infra | `pipeline_config.py` | Active | Merge default, subject, and optional tuning YAML; validate stage settings; resolve paths; guard outputs. |
| Infra | `dataset_inventory.py` | Active | Optional read-only preflight and durable source manifest. |
| Runner | `run_pipeline.py` | Active through dRSA export | Write provenance, run preflight, and execute stable CLIs in separate processes. |
| 02 | `02_trigger_extraction_custom.py` | Active | Reconstruct STI101 events and evidence-based behavioral-segment associations. |
| 03_A | `03_A_MaxFilt.py` | Active | Per-run bad-channel detection, explicit policy selection, Maxwell filtering, and QC. |
| 03_B | `03_B_ArtefactAnnotation.py` | Active | Muscle annotation and modality-aware EOG diagnostics/annotation. |
| 03_C | `03_C_ICA.py` | Active | Production manual ICA review and apply workflow. |
| 03_C auto | `03_C_ICA_auto.py` | Deprecated experimental | Retained for future research only; blocked by default and absent from the runner. |
| 04 | `04_epoching.py` | Active | Lineage-matched, run-isolated trial selection and concatenated condition epoching. |
| 06 | `06_export_neural_for_drsa.py` | Active | Identity-averaged MAT/JSON export on a fixed MEG channel axis. |
| Optional | `XX_ERF_manual.py` | Optional inspection | Interactive ERF inspection using resolved epochs and event definitions. |
| Optional | `sanity_checks/01_dot_side_decoding.py` | Optional sanity check | Configuration-driven sliding-window dot-side decoding for observed and expected/predicted trajectory sources. |
| Optional | `sanity_checks/03_dot_side_time_generalization.py` | Optional sanity check | Train-time x test-time dot-side decoding matrix for observed and expected/predicted trajectory sources. |

## Configuration Layers

Configuration is resolved in this order:

1. `meg/params/pipeline_default.yaml`
2. optional site paths such as `meg/params/pipeline_cluster.yaml`
3. subject selector such as `meg/params/sub-04.yaml`
4. optional accepted tuning YAML supplied with `--tuning-config`

Interactive stages update a guarded YAML file, defaulting to:

```text
<derivatives reports dir>/03_tuning_choices.yaml
```

That file contains only the selected subject and accepted decisions, for
example:

```yaml
subject:
  id: "04"
maxwell_filter:
  bad_channel_policy: union
  reference_acquisition_run: 1
  manual_bad_channels:
    add: []
    remove: []
    by_run: {}
  destination:
    policy: central_run
    acquisition_run: null
artifact_annotation:
  muscle_threshold_by_run:
    run01: 10.0
    run02: 8.0
ica:
  manual_exclude: [0, 5]
```

`ica.manual_exclude: null` means unreviewed. `manual_exclude: []` means the ICA
was reviewed and no components should be removed.

## Review And Apply

Production stages support two passes.

### Cluster site override

Machine-specific paths are a separate configuration layer:

```bash
python meg/code/preprocessing_basic/run_pipeline.py \
  --config meg/params/pipeline_default.yaml \
  --site-config meg/params/pipeline_cluster.yaml \
  --subject-config meg/params/sub-02.yaml \
  --stages maxwell_filter artifact_annotation ica_manual \
  --processing-pass review --resume
```

The merge order is default -> site -> subject -> tuning. The cluster override
contains paths only; scientific settings remain in `pipeline_default.yaml`.
The documented SGE array workflow is in `meg/cluster/README.md`.

### Review pass

```bash
python meg/code/preprocessing_basic/run_pipeline.py \
  --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-01.yaml \
  --stages maxwell_filter artifact_annotation ica_manual \
  --processing-pass review
```

- 03_A and 03_B create review derivatives using configured or FLUX defaults.
- 03_C fits or reuses the ICA model and writes component diagnostics.
- A missing ICA exclusion list produces `review_required`, not an implicit
  cleaning decision.

Add `--interactive` to prompt for each bad-channel policy, each acquisition
run's muscle threshold, and the ICA exclusion list. Plots are saved before
every prompt. `--no-display` retains terminal prompts but suppresses GUI
windows, allowing an operator to inspect the saved files separately.

### Apply pass

```bash
python meg/code/preprocessing_basic/run_pipeline.py \
  --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-01.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/03_tuning_choices.yaml \
  --stages maxwell_filter artifact_annotation ica_manual \
  --processing-pass apply --overwrite
```

The explicit `--overwrite` is required to replace first-pass scientific
derivatives. Pass-specific reports remain separate.

### Portable report rendering

The production stage CLIs retain lightweight HTML generation so report writing
does not become a dependency of scientific processing. When review processing
runs on the cluster, those HTML files can contain cluster-absolute PNG links
that do not resolve from the local mount.

`render_stage03_mne_reports.py` rebuilds the presentation layer only. It reads
the completed pass-specific JSON and PNG files, embeds the figures in
self-contained `mne.Report` HTML, and never opens or modifies FIF derivatives:

```bash
python meg/code/preprocessing_basic/render_stage03_mne_reports.py \
  --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
  --subjects 01 02 03 04 05 06 07 \
  --processing-pass review
```

Outputs use the suffix `_mne_report.html`, leaving the original stage HTML
unchanged. When a subject's `03_tuning_choices.yaml` exists, every rendered
report includes its complete contents plus the section relevant to that stage.
This keeps report evidence and the decisions intended for the apply pass
visible together.

### Interrupted runs

SMB or other network-volume interruptions can leave only some acquisition runs
completed. All production 03 stages support resume without replacing completed
derivatives:

```bash
python meg/code/preprocessing_basic/run_pipeline.py \
  --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-01.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/03_tuning_choices.yaml \
  --stages maxwell_filter artifact_annotation ica_manual \
  --processing-pass review --resume
```

Maxwell resume recomputes bad-channel detection evidence, then reads every
existing SSS FIF header. It skips a file only when the stored source filename,
acquisition-run label, stage identity, and applied bad-channel list all match
the current resolved decision.

Artifact-annotation resume validates source/stage lineage, the accepted muscle
threshold, channel identity, bad-channel metadata, and FIF/CSV annotation-count
consistency. ICA apply resume validates the post-rejection model, exclusion
list, source/stage lineage, channel identity, bad channels, sample count,
sampling frequency, and MISC channels. `--resume` and `--overwrite` are
mutually exclusive.

## Acquisition-Run Lineage

MEG FIFs are acquisition runs. Behavioral MAT/CSV and EDF files are organized
by experiment block and are never zipped positionally with MEG files.

Every 03-stage derivative preserves the complete source stem:

```text
sub01_run02_block02_03_150426.fif
sub01_run02_block02_03_150426_raw_sss.fif
sub01_run02_block02_03_150426_ann.fif
sub01_run02_block02_03_150426_ica.fif
```

Each output also records the source filename and acquisition-run label in:

- its FIF `info["description"]`,
- the pass-specific JSON report,
- `source_manifest.json` under `derivative_lineage`.

Refreshing the source manifest preserves existing trigger-derived associations
and derivative-lineage sections.

## 03_A Maxwell Filtering

Inputs are discovered with the configured task FIF glob. Facility calibration
and cross-talk files are explicit project resources:

```text
meg/data/fwmegmaxfilter/ct_sparse_triux2.fif
meg/data/fwmegmaxfilter/sss_cal_3182_250410.dat
```

Bad channels are detected independently for every acquisition run using
`find_bad_channels_maxwell`. The report always includes:

- raw-header bads,
- automatically noisy and flat channels,
- missing and extra channels versus the reference run,
- final applied bad channels.

Configurable policies:

- `union`: apply the union of every run's detections to all runs. This is the
  production default.
- `reference_run`: apply one reference run's list to all runs.
- `per_run`: apply each acquisition run's own detections.

Global and per-run manual additions/removals are supported. Maxwell filtering
reattaches the applied list to every output header. Guarded resume validates
that header before treating an existing derivative as completed.

Every raw run does contain a static `dev_head_t`. The default
`destination.policy: central_run` chooses an actual acquisition run as the
common Maxwell destination. Selection is the medoid minimizing summed pairwise
RMS displacement of the subject's measured digitization points under each
run's transform. This incorporates translation and rotation without creating
an artificial averaged transform. `acquisition_run` remains `null` because the
destination is always selected from the measured transform centrality.

The detected central runs are:

| Subject | Common destination |
| --- | --- |
| sub-01 | run02 |
| sub-02 | run05 |
| sub-03 | run03 |
| sub-04 | run04 |
| sub-05 | run02 |
| sub-06 | run03 |
| sub-07 | run03 |

Each run report contains both its bad-channel detection PSD and a matched
raw-versus-Maxwell PSD comparison for magnetometers and gradiometers.

## 03_B Physiological Artifact Annotation

Default muscle parameters retain the local FLUX workflow:

- threshold: `10`
- channel type: magnetometer
- filter: `110-140 Hz`
- minimum good segment: `0.2 s`

Interactive mode displays/saves the z-score trace with the current threshold
and repeats until the operator accepts a value. Thresholds may differ by
acquisition run. The accepted value is also stored in the output FIF
description so resume can reject a derivative created with a different
threshold.

EOG behavior:

- `auto` resolves to `diagnostics_only` when configured EOG channels exist.
- Sub01 therefore receives EOG diagnostics but no default EOG-derived blink
  rejection.
- Subjects 02-07 contain no EOG and skip EOG-only functions cleanly.
- An explicit `annotate` mode remains available if scientifically approved.

EDF files are not consumed by this stage. MISC001-MISC006 are retained with
their acquisition timing and are never relabeled as EOG. Trigger 150 is not a
physiological artifact annotation; its trial-level handling remains downstream.

## 03_C Manual ICA

The production ICA fit copy uses:

- MEG channels only,
- resampling to `200 Hz`,
- band-pass `1-40 Hz`,
- FastICA,
- 30 components,
- random state `97`.

New annotated acquisition runs inherit the common Maxwell destination, so
their `dev_head_t` values should agree. Concatenation retains
`on_mismatch="warn"` as defensive compatibility for historical derivatives,
but a rebuilt production run should not rely on it. The prefit model is
reusable only when a fingerprint of stable input filenames, sizes,
modification times, and fit parameters matches. Absolute mount paths are
excluded so shared derivatives remain compatible between local
`/Volumes/MORWUR` and cluster `/mnt/storage` access.

Review outputs include:

- ICA source traces,
- component topography pages,
- per-component property plots,
- EOG score rankings when EOG channels exist.

EOG scores are diagnostic and never auto-applied. Subjects without EOG skip
the diagnostic. The former hardcoded ICA000 blink annotation is removed.

The apply pass requires an explicit exclusion list, saves a post-rejection
model, applies it to every original annotated run, and preserves bad-channel
and MISC metadata. Interrupted apply passes can use guarded resume without
recomputing already validated cleaned runs.

## 04 Epoching

Epoching first enforces the Stage-03 completion gate:

- `source_manifest.json` must contain `ica_manual` with `processing_pass=apply`
  and `status=completed`;
- ICA exclusions must be explicit, where `[]` is valid;
- every configured raw acquisition run must have exactly one completed lineage
  record and one local `<exact-source-stem>_ica.fif`;
- the manifest basename must match the raw source stem exactly. Cluster mount
  prefixes are treated as provenance and relocated through the configured
  derivative root.

Each cleaned run is read and epoched independently. Trial intervals are
`[condition onset, next condition onset)` and the final interval ends at that
run's sample boundary, so evidence cannot leak into the next acquisition run.
Resolved transition-warning rows are accepted; unresolved transition evidence
is excluded by default.

Default trial-level removals are:

- trigger 150 gaze-break evidence;
- catch type 1 evidence: 116/117;
- catch type 2/question evidence: 113/118/201/202;
- ESC termination evidence: 152;
- unresolved transition evidence.

Replay trials are included because they were actually presented and saved.
Replay status remains in epoch metadata. Physiological `BAD_*` annotations and
amplitude thresholds are handled separately by MNE epoch rejection. Trigger
150 is never reinterpreted as a physiological annotation.

The production epoch settings are:

- `tmin=-0.5 s`, intentionally replacing the legacy `-1.0 s`;
- `tmax=2.667 s`;
- `baseline=None`, so epoching performs no baseline correction;
- `detrend=1`;
- `reject_by_annotation=true`.

Surviving per-run epochs are concatenated only after source/run metadata has
been attached. Metadata includes source FIF, acquisition run, raw and
concatenated sample frames, aligned behavioral block/run/trial, replay state,
transition evidence, identity, condition, trigger evidence, exclusion reasons,
and annotation policy. The complete trial-decision CSV also retains trials
removed before or during MNE epoch creation.

## 06 dRSA Neural Export

The export stage consumes only completed step-04 epochs and epoching lineage.
Identity values come from the subject profile: 1-20 for subjects 01-03 and
1-40 for subjects 04-07. One exact observed trajectory MAT is discovered from
each external subject's `input_files` folder and recorded with the profile
stimulus version. The expected occluded-deviant model uses the matching
`occluded_nondeviant` path from that same MAT, keyed by sequence identity.

The output tensor contract is:

```text
identity x canonical MEG channel x time
```

The canonical channel order comes from the repository-owned
`drsa_canonical_meg_channels.txt` and is never initialized or changed during
export. Channels marked bad upstream or absent from an epochs file remain on
the axis as NaN. MISC, EOG, STI, and other non-MEG channels are omitted without
modifying the source epochs or FIFs.

The default missing-identity policy is `error`: a required identity with no
surviving epoch in a condition stops export. A separately configured `nan`
policy is available and leaves the identity row explicitly missing. No 40-ID
profile is truncated to 20 and no identity is synthesized silently.

Each condition MAT and JSON sidecar records identity values and counts,
trajectory paths/version, event map, source lineage, inclusion policy, ordered
channel labels, time axis, sampling rate, bad/missing channel labels and NaN
policy. Aggregate JSON/HTML reports and source-manifest lineage are also
written.

## 07 MATLAB dRSA

Stage 07 is optional and independently rerunnable. It consumes completed
Stage-06 MAT/JSON exports and the completed `export_drsa` source-manifest
lineage. It does not modify epochs, neural exports, or raw data.

`07_run_matlab_drsa.py` validates the Stage-06 products and writes a
subject-level handoff:

```text
derivatives/reports/sub-XX/ses-01/07_matlab_drsa_handoff.json
```

That handoff is the only path/provenance contract consumed by MATLAB. It
contains Stage-06 neural paths, sidecar paths, trajectory model paths, identity
values, time axis, sampling rate, condition definitions, channel labels, and
output directories. MATLAB scripts must not hardcode V28, exactly 20
identities, local-only roots, or model filenames.

Registered optional analysis families are:

- `real`: `A_dRSA.m`, the main real-neural dRSA. Outputs are written under
  `analysis/sub-XX/ses-01/drsa/matlab_real_neural/07_dRSA/`.
- `alt`: `A_dRSA_alt.m`, single-model/sensor-group exploratory dRSA. For
  `occluded_deviant`, `--model-sources` can request `observed`, `predicted`,
  and `prediction_error`; the prediction-error model is actual deviant
  position minus expected/predicted nondeviant position.
- `simulations`: `dRSA_simulations.m`, model-vs-model and neural-vs-neural
  simulation-style checks.

Each family writes PNG matrix figures and MAT result files under the subject's
configured `drsa` derivative directory. The Python stage scans those outputs
and writes:

```text
derivatives/reports/sub-XX/ses-01/07_matlab_drsa_report.json
derivatives/reports/sub-XX/ses-01/07_matlab_drsa_report.html
```

The HTML report embeds discovered PNG figures. Missing or failed MATLAB work
is represented in the JSON/HTML status rather than being silently treated as
completed.

Add `--smoothing-plots` (or `DAD_DRSA_SMOOTHING_PLOTS=1`) to write
display-only `_smoothed.png` variants next to the raw matrix PNGs. The
underlying MAT result is still the unsmoothed dRSA matrix; smoothing is for
figure readability and uses `plot_smoothing_window` samples.

After all subject Stage-07 outputs exist, a one-time group report can be built
from the numeric MAT results:

```bash
python meg/code/preprocessing_basic/07_run_matlab_drsa.py \
  --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-01.yaml \
  --group-average-only \
  --group-subjects 01 02 03 04 05 06 07 \
  --smoothing-plots
```

Group outputs are written under:

```text
derivatives/analysis/group/ses-01/drsa/group_average/
derivatives/reports/group/ses-01/07_matlab_drsa_group_report.json
derivatives/reports/group/ses-01/07_matlab_drsa_group_report.html
derivatives/reports/group/ses-01/07_matlab_drsa_group_report_smoothed.json
derivatives/reports/group/ses-01/07_matlab_drsa_group_report_smoothed.html
```

Raw group figures average unsmoothed subject matrices. When smoothing is
enabled, the group report also averages the per-subject smoothed matrices into
separate `_smoothed.png` figures. The raw and smoothed reports are separate:
each report starts with group-average figures and then includes
participant-by-participant dRSA figures organized by subject and condition.

Cluster execution uses:

```bash
qsub -v DAD_MATLAB_BIN=/absolute/path/to/MATLAB/R2020a/bin/matlab \
  -t 1-7 meg/cluster/submit_drsa_analysis.sh
```

`DAD_MATLAB_BIN` must be an absolute executable path. Local MATLAB R2020a on
Apple Silicon may require `arch -x86_64`, which the Stage-07 CLI exposes as
`--use-arch-x86-64`.

## Deprecated Automatic ICA

`03_C_ICA_auto.py` is retained for future evaluation but:

- exits with a deprecation error by default,
- is absent from `SUPPORTED_STAGES`,
- is not installed through production requirements,
- is not run on participants,
- never acts as a fallback when manual ICA fails.

Its optional research dependencies are listed separately in
`meg/requirements-ica-auto-experimental.txt`.

## Ocular And Trial Boundaries

- Sub01: EOG001/EOG002 in FIF; saved EDF unavailable.
- Subjects 02-07: no EOG; saved EDF exists per behavioral block.
- Subjects 02-07: MISC001-MISC003 carry left-eye analog signals;
  MISC004-MISC006 remain preserved despite being near-flat.
- EDF-to-MEG synchronization is deferred and must not assume matching file
  starts or positional block/run correspondence.
- Trigger 150 records online gaze breaks independently of EDF availability.
- Trigger-150 trial rejection is separate from 03_B/03_C physiological
  cleaning.

## Baseline Facts

- `04_epoching.py`: `tmin=-0.5 s`, `baseline=None`, linear detrending.
- `sanity_checks/01_dot_side_decoding.py`: `-0.4` to `-0.2 s`.
- `XX_ERF_manual.py`: default `-0.1` to `0 s`.

Both optional baseline windows fit inside the epoch and are compatible with
the verified 0.5-1.0 s sub01 ITI and 1.0-1.5 s subject 02-07 ITIs. Baseline
correction occurs only in those optional copied inspection objects, not in
production epoching.

## Optional Dot-Side Decoding Sanity Check

`sanity_checks/01_dot_side_decoding.py` is an optional diagnostic workflow and
is not registered as a required `run_pipeline.py` stage. It consumes completed
Stage-04 epochs plus `source_manifest.json`, validates that the epoching
lineage is completed, validates the configured condition map and
`baseline=None` epoching contract, and then reads the exact subject-specific
trajectory MAT file through the resolved configuration.

Four standard analyses are supported and namespaced independently:

- `always_visible_observed`;
- `occluded_nondeviant_observed`;
- `occluded_deviant_observed`;
- `occluded_deviant_predicted`.

Predicted deviant decoding uses `occluded_deviant` neural epochs but derives
the dynamic dot-side labels from the matching `occluded_nondeviant` trajectory
in the same MAT file, following the established expected-path convention used
by dRSA. Condition, model source, model-label condition, and analysis key are
recorded in every subject summary. Outputs are separated under
`<condition>/<model-source>` so observed and predicted deviant results cannot
overwrite each other.

The decoder uses only MEG gradiometer features selected from posterior sensors;
MISC, EOG, STI, and other non-MEG channels are not feature inputs. The ERF
inspection branch uses copied MEG epochs with a local `-0.4` to `-0.2 s`
baseline and does not alter production epochs.

Decoding uses centered sliding windows on the resampled decoding grid. The
default is `--decode-window-samples 10` after `--decode-sfreq 120.0`; odd
widths are symmetric, even widths place the extra sample after the center, and
epoch-edge windows are clipped to available samples. The scores CSV records
center time, start/end sample indices, start/end times, effective sample count,
effective duration in milliseconds, edge-clipping state, class counts, skip
reason, and unbalanced/balanced accuracies.

For occluded conditions, participant decoding panels and group accuracy plots
mark the first completely occluded frame and first reappearance frame. The
boundaries are resolved from `matlab_drsa.occlusion_timing` using the shared
1-based-frame conversion `(frame - 1) / frame_frequency_hz`; the default
frames 130 and 191 therefore appear at 1.075 s and 1.583 s.

`--trials-per-average` optionally averages source trials before
classification, but only within the selected condition and the same
`trial_sequence_identity`. Numeric values use fixed-count pseudo-trials: source
trials are shuffled deterministically with `--random-seed`, complete
non-overlapping groups are averaged, and incomplete remainder groups are
dropped. If one identity has fewer trials than the requested group size, that
identity is fully dropped and recorded while other identities continue. The
literal value `all_trials_per_identity` instead averages every retained source
trial for an identity into one identity-level pseudo-trial. The run fails only
when no complete pseudo-trials remain. The summary JSON records averaging mode,
dropped identities/counts, source epoch indices for every retained pseudo-trial,
and training/testing diagnostics for the classifier.

Labels are always model-derived from the exact subject trajectory file and
the configured profile identity range. The validation path supports both
20-identity subjects 01-03 and 40-identity subjects 04-07 and fails on model
identity mismatches rather than truncating or synthesizing identities.

Cluster execution is provided by:

```bash
qsub -t 1-7 meg/cluster/submit_dot_side_decoding.sh
```

By default the wrapper runs all four standard analyses sequentially inside
each participant task. It exposes analysis selection, decode time grid,
decode sampling rate, window width, trial averaging, seed, CV folds,
validate-only, resume, rebuild, and output-directory overrides through
documented environment variables. It uses `meg/params/pipeline_cluster.yaml`,
micromamba activation, headless matplotlib/MNE settings, `cluster_logs`, and
the shared raw-tree output guard.

`sanity_checks/02_dot_side_group_report.py` builds one optional combined report
after subject-level outputs exist. For each analysis it reads every subject's
summary JSON and scores CSV, validates that condition, model source,
model-label condition, time axis, decode window, sampling rate, fold count,
seed, and trial-averaging settings match, and tests the selected metric against
chance over time.
When subject-level decoding marks a time point as NaN because class labels are
not well-defined for the configured stratified folds, group statistics exclude
that time point unless every included subject has a finite score. Excluded time
points remain gaps in the group plot, so temporal clusters cannot bridge
skipped windows.

The default group metric is
`accuracy_balanced_undersampled_sets`, preserving the subject-level
undersampling strategy as the primary sanity-check readout. The default
statistical test is a one-sample sign-flip cluster permutation on
`accuracy - 0.5`, using a directional greater-than-chance alternative,
cluster-forming `p=0.05`, and cluster-level alpha `0.05`. With subjects 01-07
the sign-flip distribution is exact by default. The group figure annotates the
cluster-forming p threshold and equivalent t cutoff, marks significant
clusters, and prints cluster p-values on the plot.

The default group display band is a 95% confidence interval around the mean.
SEM is available explicitly with `--interval sem`; the report records which
interval type was used. The HTML report contains all four group figures and
cluster tables first. Only after the complete group section does it include
subjects 01-07. Each subject section begins with a three-condition source-trial
count figure by trajectory identity, followed by all requested three-panel
analyses. Each group and participant result includes a training/testing
diagnostic table covering the classifier, stratified CV folds, averaging mode,
pseudo-trial/identity counts, class-balance summaries, skipped-window reasons,
and an assessment of whether the decoded time points are fully or only
partially well-defined.
Counts are taken before optional pseudo-trial averaging; observed and predicted
deviant summaries must agree because they consume the same neural epochs.
The group wrapper is:

```bash
qsub meg/cluster/submit_dot_side_group_report.sh
```

The 1D group report can optionally smooth each subject's finite accuracy
segments before group statistics and plotting with `--smooth-sigma-samples`.
The smoothing is Gaussian, is applied before the sign-flip cluster test, and is
never allowed to cross NaN skipped-window gaps. At the default 120 Hz decoding
grid, `--smooth-sigma-samples 1.0` corresponds to one decoded sample of
Gaussian sigma, about 8.3 ms, with an approximate full-width at half maximum of
2.355 samples or 19.6 ms. The cluster wrapper exposes the same setting through
`DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES`.

`sanity_checks/03_dot_side_time_generalization.py` is the matrix companion to
the line-decoding sanity check. It preserves the same four analysis keys,
trajectory-label rules, `--decode-window-samples` sliding-window definition,
and same-identity trial averaging options. For each train-time center it builds
stratified folds from the train-time labels once, fits the regularized LDA
inside each fold, and reuses those folds/models across every test-time center.
The primary output matrix stores `accuracy_balanced_undersampled_sets` with
rows as train time and columns as test time.

Train rows are skipped when train-time labels are single-class or the minority
class cannot support the configured `--n-splits`. Test cells are NaN when the
held-out test labels cannot be balanced in every reused fold. Because
time-generalization requires train and test feature vectors to have the same
length, cells crossing different edge-clipped window widths are also explicit
NaNs with a `feature_dimension_mismatch` reason. The NPZ bundle is the
reload-friendly matrix format; the companion long CSV provides one row per
train/test cell for auditing.

Participant outputs are written under:

```text
derivatives/sanity_checks/sub-XX/ses-01/dot_side_time_generalization/<condition>/<model-source>/
```

Each participant bundle contains a summary JSON, matrix NPZ, long matrix CSV,
heatmap PNG, and MNE HTML report. Occluded-condition heatmaps mark first
complete occlusion and first reappearance on both train-time and test-time
axes.

`sanity_checks/04_dot_side_time_generalization_group_report.py` builds the
group companion. It validates that included subjects share matrix shape, time
axis, decoding settings, random seed, and trial-averaging settings. Group
statistics use a conservative all-subject finite-cell mask, so clusters cannot
bridge NaN cells from skipped participant train/test combinations. The
inferential test is a one-sample sign-flip cluster permutation over the 2D
train-time x test-time lattice using 4-neighbor adjacency; the default
cluster-forming threshold is `p=0.05`, cluster alpha is `0.05`, and subjects
01-07 use exact sign flips by default. Group plots show the valid-cell group
mean matrix with significant 2D clusters shaded and outlined in red.
Detailed 2D cluster tables are omitted from HTML reports by default; add
`--include-cluster-table` to the group-report script, or set
`DAD_DOT_SIDE_TG_INCLUDE_CLUSTER_TABLE=1` in the cluster wrapper, when the
per-cluster train/test ranges, masses, and p-values should be shown inline.
Optional `--smooth-sigma-samples` applies isotropic Gaussian smoothing to each
subject matrix before group statistics and plotting. Smoothing is constrained
to finite connected train/test components, so skipped cells remain NaN and 2D
clusters cannot bridge those gaps. The time-generalization group wrapper
exposes `DAD_DOT_SIDE_TG_SMOOTH_SIGMA_SAMPLES`, falling back to the shared
`DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES` when unset.

Cluster wrappers are:

```bash
qsub -t 1-7 meg/cluster/submit_dot_side_time_generalization.sh
qsub meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

## Verification Contract

Required lightweight checks before participant processing:

- compile all preprocessing scripts and tests,
- run all pytest tests,
- validate configuration for subjects 01-07,
- run pipeline dry-run for subjects 01-07,
- inspect all raw FIF headers read-only,
- run synthetic MISC/bad-channel/ICA preservation smoke tests,
- run representative short real-data Maxwell/artifact/ICA smoke tests.
- run synthetic run-isolated epoching and dRSA export smoke tests,
- validate Stage-07 handoff generation for subjects 01-07 after Stage 06,
- run MATLAB `checkcode` on Stage-07 analysis scripts before cluster
  submission,
- run dot-side line and time-generalization sanity-check helper and wrapper tests,
- run `bash -n` for all cluster wrappers,
- validate the dot-side sanity-check handoff for at least one 20-identity and
  one 40-identity subject before participant array jobs.

Full stage processing is run only after those checks pass.

Chat-6 dot-side verification on June 16, 2026:

- `python -m compileall -q meg/code/preprocessing_basic meg/code/sanity_checks meg/tests`: passed.
- `python -m pytest meg/tests/test_dot_side_decoding.py meg/tests/test_dot_side_group_report.py meg/tests/test_cluster_wrappers.py -q`: `15 passed`.
- `bash -n meg/cluster/*.sh`: passed for all cluster wrappers.
- `01_dot_side_decoding.py --validate-only --print-summary`: passed for sub01
  as a 20-identity subject and sub04 as a 40-identity subject.
- Local `pipeline_config.py --validate-only` with default paths passed for
  subjects 01-07.
- `pipeline_config.py --site-config meg/params/pipeline_cluster.yaml
  --validate-only` cannot complete on the laptop because the cluster raw root
  `/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data` is not mounted
  locally; this validation must be run on the CIMeC cluster before array
  submission.
- A local full subject run completed for sub05 `always_visible` with
  `--decode-window-samples 10`, `--trials-per-average 1`, and seed `42`,
  writing the three-panel PNG, scores CSV, summary JSON, and MNE report under
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/sub-05/ses-01/dot_side_decoding/always_visible/`.
- Local subject-level dot-side outputs were generated or resumed for subjects
  01-07. The local group report completed for `always_visible` with balanced
  accuracy, 95% CI, cluster-forming `p=0.05`, exact sign flips (`n=128`),
  and four significant above-chance clusters. Outputs were written under
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/dot_side_decoding/always_visible/`
  and
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_decoding/always_visible/`.

Four-analysis decoding extension verification on June 22, 2026:

- `python -m compileall -q meg/code/preprocessing_basic meg/code/sanity_checks meg/tests`: passed.
- `python -m pytest meg/tests/test_dot_side_decoding.py meg/tests/test_dot_side_group_report.py meg/tests/test_cluster_wrappers.py -q`: `30 passed`.
- `python -m pytest meg/tests -q`: `140 passed`, with six pre-existing MNE
  filename-convention warnings from synthetic FIF fixtures.
- `bash -n` passed for both dot-side cluster wrappers, and targeted
  `git diff --check` passed.
- Read-only `01_dot_side_decoding.py --validate-only` passed all four standard
  analyses for sub01 (20 identities) and sub04 (40 identities), including
  deviant neural epochs with predicted/nondeviant model labels.
- A synthetic two-subject run built all four group figures, combined JSON, and
  embedded HTML and verified that the complete group-analysis section precedes
  every participant-specific section.
- No participant decoding outputs or existing derivatives were regenerated as
  part of this implementation verification.
