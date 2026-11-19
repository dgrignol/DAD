# Archived DAD MEG Pipeline SGE Workflow

> **Historical reference only:** this file preserves the former SGE `qsub`
> workflow. The active wrappers now use SLURM; do not run the commands in
> this archive. See [`README.md`](README.md) for current DIPSEN instructions.

## Purpose

This file documents deployment of DAD MEG preprocessing, downstream export,
optional MATLAB dRSA, and optional sanity checks. The job workflow below is
the historical CIMeC Sun Grid Engine setup, with one subject per array task.

The automatic MEGNet branch is deprecated and is never submitted.

## Files

- `submit_preprocessing_stage1.sh`: cohort source/trigger review and the separate P24 acquisition audit.
- `submit_preprocessing_review.sh`: Stage-03 review array job.
- `submit_preprocessing_apply.sh`: Stage-03 apply array job.
- `submit_epoch_export.sh`: Stage-04 epoching then Stage-06 export array job.
- `submit_drsa_analysis.sh`: optional Stage-07 MATLAB dRSA array job.
- `submit_drsa_group_average.sh`: optional Stage-07 numeric group-average report job.
- `submit_drsa_line_summaries.sh`: optional Stage-08 dRSA line-summary report job.
- `submit_drsa_stats.sh`: optional Stage-10 dRSA matrix and line-summary stats job.
- `submit_dot_side_decoding.sh`: optional dot-side sanity-check array job.
- `submit_dot_side_group_report.sh`: optional dot-side group report job.
- `submit_dot_side_time_generalization.sh`: optional dot-side matrix decoding array job.
- `submit_dot_side_time_generalization_group_report.sh`: optional dot-side matrix group report job.
- `../params/pipeline_cluster.yaml`: cluster-only path overrides.
- `../params/pipeline_default.yaml`: shared scientific defaults.
- `../params/sub-XX.yaml`: subject selectors.

## Current DIPSEN connection

The tested login command and the retained project directory are:

```bash
ssh damiano.grignolio@hpc.dipsen.unitn.it
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
```

Run `meg/cluster/deploy_review_to_cluster.sh` from the local DAD directory to
copy the review payload to this host and directory. Compute jobs require the
SLURM workflow described by DIPSEN; the DAD job wrappers are not ported yet.

## One-Time Deployment

The repository currently lives on the laptop, while the cluster project files
are under `/mnt/storage/tier2/morwur`. Only a small review-pass payload is
required:

- production Python files from `meg/code/preprocessing_basic` through export;
- optional dRSA post-processing and sanity-check files;
- `meg/params/pipeline_default.yaml`, `pipeline_cluster.yaml`, and all 39
  subject selectors;
- `meg/cluster`;
- the Maxwell cross-talk FIF and calibration DAT files.

Raw data and derivatives already live in their separate shared trees. Tests,
Git metadata, project documentation, local derivatives, and deprecated
automatic ICA dependencies are not required.

From the local DAD directory, use the deployment helper:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
meg/cluster/deploy_review_to_cluster.sh
```

The helper preserves the `DAD/meg/...` layout, does not use deletion flags,
and refuses deployment into the raw `data` tree.

This destination is outside the read-only raw-data tree. Do not deploy code
inside `Attention/data`.

## Legacy SGE job workflow: sources and triggers

The commands and job-log behavior in this section are for the old SGE cluster.
They are not valid on DIPSEN, where jobs must use SLURM after the DAD wrappers
have been ported.

After deployment, submit participants 08-39 from the cluster head node. The
existing 01-07 derivatives are left alone. Grid Engine writes job logs under
`cluster_logs`, and the Python pipeline writes reports under
`Attention/derivatives/reports/sub-XX/ses-01`.

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
qsub -tc 4 -t 8-39 meg/cluster/submit_preprocessing_stage1.sh
```

Each task writes `dataset_inventory.json/.md`, `source_manifest.json`, and
the Stage-02 reports under `derivatives/reports/sub-XX/ses-01`. The canonical
`02_trigger_events.csv` and `02_trigger_alignment.csv` live under
`derivatives/preprocessing/sub-XX/ses-01/events`.
The inventory includes the experimenter's session/run notes as **unverified
review clues** and reports the acquisition high-pass, low-pass, and channel
list from each FIF header. The trigger table keeps digital STI101 samples;
it does not yet shift epochs to photodiode time. The legacy deep trigger plot
bundle is deferred at this gate to keep the cohort job focused.

P24 additionally runs `p24_acquisition_adapter.py`. It measures the 50%
MISC008 light-rise lag only in run 5, reports a run-5 lag distribution, and
labels runs 1-4 as estimates based on that distribution. Its
`p24_visual_timing_proposal.csv` and `p24_acquisition_audit.json` are review
artifacts, not accepted corrections. The irreversible 0.1 Hz online
high-pass in runs 1-4 is listed in the audit. Review P24's photodiode and
filter choices before any Stage-03/04 application.

At this gate, compare inventory warnings, event/code counts, alignment
anomalies, acquisition-run/behavioral-block associations, and session notes
against the recordings. Compare P24's accepted run-5 lag distribution and
flagged events independently before deciding whether and how to use the
earlier-run estimates. Once those decisions are recorded, prepare the next
stage for bad channels, muscle annotations, and ICA review.

## Legacy SGE Python environment notes

The following shell and environment notes were written for the old SGE
cluster. For current DIPSEN login, use the new host and project path below;
do not SSH directly to a compute node. The DAD submission wrappers still need
SLURM support before cluster jobs can be run.

The micromamba setup and package-pinning commands that follow are inherited
from the old cluster and have not yet been verified on DIPSEN.

```bash
ssh -Y damiano.grignolio@hpc.dipsen.unitn.it
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
```

On DIPSEN, request compute resources through SLURM according to the official
cluster documentation; the old `compute-0-3.local` SSH command is no longer
appropriate.

Initialize micromamba in the new SSH shell. Environments below
`$HOME/micromamba/envs` are in the shared home directory, so an environment
created from the head node is available on compute nodes:

```bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$($HOME/bin/micromamba shell hook -s bash)"
```

For a new environment, install the pinned conda-forge binaries:

```bash
micromamba create -y -f meg/cluster/environment.yml
micromamba activate dad-meg
```

If `dad-meg` already exists, repair or update it in place:

```bash
micromamba install -y -n dad-meg --override-channels -c conda-forge \
  "python=3.11" \
  "pyyaml=6.0" \
  "numpy=1.26" \
  "scipy=1.13" \
  "mne=1.10.2" \
  "matplotlib=3.8" \
  "pandas=2.2" \
  "scikit-learn=1.5"
micromamba activate dad-meg
```

Do not rerun an unrestricted `pip install -r requirements.txt`. On this
cluster, pip may select a NumPy source distribution requiring a newer compiler
than GCC 4.8.5. `environment.yml` uses compiled conda-forge packages. The
requirements file is retained only as a pinned binary-wheel fallback.

Verify the environment and cluster paths while still on the compute node:

```bash
python -c "import mne, numpy, scipy; print(mne.__version__, numpy.__version__, scipy.__version__)"

python meg/code/preprocessing_basic/pipeline_config.py \
  --config meg/params/pipeline_default.yaml \
  --site-config meg/params/pipeline_cluster.yaml \
  --subject-config meg/params/sub-02.yaml \
  --require-section maxwell_filter \
  --require-section artifact_annotation \
  --require-section ica \
  --validate-only --print-summary

exit
```

MNE 1.10.2 is required for this dataset. MNE 1.8 and 1.9 expose only the
cardinal and HPI points from these FIF headers, causing `origin: auto` to fail
with `No digitization points found`. MNE 1.10.2 correctly exposes the embedded
Polhemus extra points used for the Maxwell sphere fit.

If an older environment was already created, repair it on a compute node:

```bash
micromamba install -y -n dad-meg --override-channels -c conda-forge "mne=1.10.2"
micromamba activate dad-meg
python -c "import mne; assert mne.__version__ == '1.10.2'; print(mne.__version__)"
```

The summary must show:

```text
Raw subject directory: /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data/sub02_130526
Derivative root: /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives
```

## Submit The Review Pass

On the cluster head node:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
qsub -t 2-7 meg/cluster/submit_preprocessing_review.sh
```

Use `2-7` to run the subjects not already reviewed locally. To include or
resume sub01 too:

```bash
qsub -t 1-7 meg/cluster/submit_preprocessing_review.sh
```

The script normally uses guarded `--resume`. If a subject already has
`03_tuning_choices.yaml`, those accepted choices are reused. Otherwise the
first pass uses:

- bad-channel policy: `union`;
- common Maxwell destination: central actual acquisition run;
- muscle threshold: `10`;
- ICA exclusions: `null`, meaning review still required.

After a scientific change that invalidates existing derivatives, explicitly
rebuild all selected stages:

```bash
qsub -v DAD_REBUILD=1 -t 1-7 \
  meg/cluster/submit_preprocessing_review.sh
```

`DAD_REBUILD=1` uses `--overwrite`; it is intentionally different from
recovery/resume. Existing subject tuning YAML is still applied, so review each
file before rebuilding. In particular, an explicit
`maxwell_filter.bad_channel_policy` in tuning overrides the new union default.

To recover after a later-stage infrastructure failure without repeating
Maxwell bad-channel detection, select a colon-separated stage subset:

```bash
qsub -v DAD_STAGES=artifact_annotation:ica_manual -t 1 \
  meg/cluster/submit_preprocessing_review.sh

qsub -v DAD_STAGES=ica_manual -t 2-7 \
  meg/cluster/submit_preprocessing_review.sh
```

Only `maxwell_filter`, `artifact_annotation`, and `ica_manual` are accepted.
Dependencies must already exist when an earlier stage is omitted.

## Monitor

```bash
qstat
qstat -j JOB_ID
tail -f cluster_logs/DAD_meg_review.JOB_ID.TASK_ID.log
```

Cancel a job only when required:

```bash
qdel JOB_ID
```

Each successful log ends with the subject report and figure directories.

## Inspect The First-Pass Outputs

Reports:

```text
/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/reports/sub-XX/ses-01/
```

Decision plots:

```text
/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/preprocessing/sub-XX/ses-01/figures/
```

Inspect, in order:

1. `03_A_MaxFilt_review_report.html`, every detection PSD, and every
   `*_before_after_psd.png`. Confirm the selected central destination, union
   bad-channel list, and raw-versus-Maxwell attenuation.
2. `03_B_ArtefactAnnotation_review_report.html` and every
   `*_muscle_zscore.png`. Decide whether threshold `10` cleanly separates
   bursts from background for each acquisition run.
3. `03_C_ICA_review_report.html`, source traces, topography pages, and all
   component property plots. EOG ranking is diagnostic only. Choose exclusions
   from time course, topography, spectrum, and segment consistency.

For subjects 02-07, EOG-only diagnostics must be reported as skipped. Their
MISC001-MISC006 channels must remain present. Trigger 150 is not handled by
these physiological cleaning stages.

### Render portable MNE reports locally

The stage CLIs write lightweight HTML that links its figures. Reports produced
on the cluster contain cluster-absolute image paths, so those links do not
resolve through the local `/Volumes/MORWUR` mount. The completed JSON metadata
and PNG figures can be converted locally without opening FIF files or rerunning
any scientific processing:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD

python meg/code/preprocessing_basic/render_stage03_mne_reports.py \
  --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
  --subjects 01 02 03 04 05 06 07 \
  --processing-pass review
```

This creates self-contained reports beside the original HTML:

```text
03_A_MaxFilt_review_mne_report.html
03_B_ArtefactAnnotation_review_mne_report.html
03_C_ICA_review_mne_report.html
```

Every PNG is embedded in the HTML, so the reports remain usable across mount
layouts. If `03_tuning_choices.yaml` is present for a subject, each report also
contains the complete tuning file and highlights the choices relevant to that
stage. Re-run with `--overwrite` only when intentionally replacing an earlier
rendered MNE report.

## Record Decisions

Create or edit:

```text
.../derivatives/reports/sub-XX/ses-01/03_tuning_choices.yaml
```

Example:

```yaml
subject:
  id: "02"

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
    run03: 10.0
    run04: 10.0
    run05: 9.0

ica:
  manual_exclude: [0, 5]
```

Use `manual_exclude: []` when the ICA was reviewed and no component should be
removed. Do not leave it `null` for the apply pass.

## Run The Apply Pass

Launch apply only after every submitted subject has a reviewed tuning file:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
qsub -t 1-7 meg/cluster/submit_preprocessing_apply.sh
```

The apply wrapper requires `03_tuning_choices.yaml`, sets
`--processing-pass apply --overwrite`, rebuilds Maxwell and annotation
derivatives from the accepted choices, and creates final `*_ica.fif` files.
Missing tuning or `ica.manual_exclude: null` fails explicitly.

The apply wrapper resolves its delegated review wrapper from
`DAD_PROJECT_ROOT` (defaulting to the documented cluster project path). This
is required because Grid Engine executes a temporary spooled copy of the
submitted apply script, so the runtime script directory is not the repository.

After completion, recreate portable local MNE reports for the apply pass:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD

python meg/code/preprocessing_basic/render_stage03_mne_reports.py \
  --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
  --subjects 01 02 03 04 05 06 07 \
  --processing-pass apply --overwrite
```

## Run Stage 04 Then Stage 06

Resync the cluster after changing epoching/export/MATLAB dRSA/dRSA
post-processing/sanity-check code. The
deployment helper includes `04_epoching.py`, `06_export_neural_for_drsa.py`,
`07_run_matlab_drsa.py`, the dRSA post-processing scripts, `meg/code/analysis`,
`meg/code/sanity_checks`, `simulations/functions`, `downstream_stage_utils.py`,
the fixed channel list, current configuration, and the downstream submission
wrappers:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
meg/cluster/deploy_review_to_cluster.sh
```

Then submit one task per subject from the cluster head node:

```bash
ssh -Y damiano.grignolio@hpc.dipsen.unitn.it
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
qsub -t 1-7 meg/cluster/submit_epoch_export.sh
```

Each task validates the completed manual-ICA apply gate, runs `epoching`, and
only after success runs `export_drsa`. Default mode uses `--resume`, which is
also appropriate for retrying a Stage-06 failure after Stage 04 completed:

```bash
qsub -t 4 meg/cluster/submit_epoch_export.sh
```

To intentionally replace existing Stage-04 and Stage-06 derivatives without
rerunning Stage 03:

```bash
qsub -v DAD_REBUILD=1 -t 1-7 \
  meg/cluster/submit_epoch_export.sh
```

Monitor jobs and logs:

```bash
qstat
tail -f cluster_logs/DAD_epoch_export.JOB_ID.TASK_ID.log
```

Successful outputs are written under:

```text
.../derivatives/analysis/sub-XX/ses-01/meg/
.../derivatives/analysis/sub-XX/ses-01/drsa/
.../derivatives/reports/sub-XX/ses-01/
```

## Run Stage 07 MATLAB dRSA

Stage 07 consumes the completed Stage-06 `drsa` MAT/JSON exports and writes
MATLAB dRSA PNG/MAT outputs plus JSON/HTML reports. The wrapper first honors
`DAD_MATLAB_BIN`, then checks the validated CIMeC installation at
`/state/partition1/MATLAB/R2020a/bin/matlab`, common install paths, `matlab` on
`PATH`, and environment modules. Submit the array with:

```bash
ssh -Y damiano.grignolio@hpc.dipsen.unitn.it
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
qsub -t 1-7 meg/cluster/submit_drsa_analysis.sh
```

The job keeps the JVM enabled because MATLAB R2020a disables `figure` under
`-nojvm`. It requests 24 GB and launches MATLAB with headless software OpenGL;
both `figure` and `exportgraphics` were validated on a CIMeC compute node. If
auto-discovery fails, resubmit explicitly with:

```bash
qsub -v DAD_MATLAB_BIN=/state/partition1/MATLAB/R2020a/bin/matlab \
  -t 1-7 meg/cluster/submit_drsa_analysis.sh
```

Default mode runs the configured real-neural dRSA analysis with `--resume`.
To intentionally replace existing Stage-07 outputs:

```bash
qsub -v DAD_REBUILD=1 \
  -t 1-7 meg/cluster/submit_drsa_analysis.sh
```

Outputs are written under:

```text
.../derivatives/analysis/sub-XX/ses-01/drsa/matlab_real_neural/
.../derivatives/reports/sub-XX/ses-01/07_matlab_drsa_report.html
.../derivatives/reports/sub-XX/ses-01/07_matlab_drsa_report.json
```

To request smoothed display figures in the subject jobs:

```bash
qsub -v DAD_DRSA_SMOOTHING_PLOTS=1 \
  -t 1-7 meg/cluster/submit_drsa_analysis.sh
```

To run the optional sensor-group alt dRSA for the occluded-deviant
prediction-error position model:

```bash
qsub -v DAD_DRSA_ANALYSES=alt,DAD_DRSA_CONDITIONS=occluded_deviant,DAD_DRSA_SENSOR_GROUPS=wholeBrain,DAD_DRSA_MODEL_SOURCES=prediction_error \
  -t 1-7 meg/cluster/submit_drsa_analysis.sh
```

Build the numeric group-average report once after all subject jobs finish:

```bash
qsub -hold_jid <stage07_job_id> meg/cluster/submit_drsa_group_average.sh
```

The group-average wrapper always runs `--group-average-only` and
`--print-summary`. It defaults to subjects `01..07`, the standard cluster
configuration files, and smoothed group figures. Override the useful knobs with
environment variables when needed:

```bash
export DAD_DRSA_GROUP_CONFIG=/path/to/pipeline_default.yaml
export DAD_DRSA_GROUP_SITE_CONFIG=/path/to/pipeline_cluster.yaml
export DAD_DRSA_GROUP_SUBJECT_CONFIG=/path/to/sub-01.yaml
export DAD_DRSA_GROUP_SUBJECTS=01,02,03,04,05,06,07
export DAD_DRSA_GROUP_SMOOTHING_PLOTS=1
export DAD_DRSA_GROUP_PLOT_SMOOTHING_WINDOW=9

qsub -v DAD_DRSA_GROUP_CONFIG,DAD_DRSA_GROUP_SITE_CONFIG,DAD_DRSA_GROUP_SUBJECT_CONFIG,DAD_DRSA_GROUP_SUBJECTS,DAD_DRSA_GROUP_SMOOTHING_PLOTS,DAD_DRSA_GROUP_PLOT_SMOOTHING_WINDOW \
  -hold_jid <stage07_job_id> meg/cluster/submit_drsa_group_average.sh
```

The group report is written to:

```text
.../derivatives/analysis/group/ses-01/drsa/group_average/
.../derivatives/reports/group/ses-01/07_matlab_drsa_group_report.html
.../derivatives/reports/group/ses-01/07_matlab_drsa_group_report.json
.../derivatives/reports/group/ses-01/07_matlab_drsa_group_report_smoothed.html
.../derivatives/reports/group/ses-01/07_matlab_drsa_group_report_smoothed.json
```

The raw and smoothed HTML reports both start with group-average figures and
then collect participant figures under subject sections with condition
subsections.

## Run Optional dRSA Line, Surface, And Stats Reports

These are optional post-processing jobs. They consume completed Stage-07 MAT
files and do not rerun MATLAB, dRSA export, epoching, or preprocessing.

Build Stage-08 group line summaries:

```bash
qsub -hold_jid <stage07_job_id> meg/cluster/submit_drsa_line_summaries.sh
```

The Stage-09 sliding lag-surface script is report-oriented and can be run
directly on the cluster after deployment:

```bash
python meg/code/preprocessing_basic/09_plot_drsa_sliding_lag_surfaces.py \
  --config meg/params/pipeline_default.yaml \
  --site-config meg/params/pipeline_cluster.yaml \
  --subjects 01 02 03 04 05 06 07 \
  --group-only --resume --print-summary
```

Build Stage-10 dRSA stats as a single group job:

```bash
qsub meg/cluster/submit_drsa_stats.sh
```

The Stage-10 wrapper defaults to subjects `01..07`, both matrix and line-summary
targets, position plus prediction-error model selectors, all available
conditions and sensor groups, analysis families `A_dRSA,A_dRSA_alt`, the actual
observed-position paired contrast `occluded_deviant - occluded_nondeviant`,
significant-only cluster tables in the HTML report, `1024` sign-flip
permutations, cluster-forming `p=0.05`, cluster p significance `0.05`, and a
two-sided test. Simulation families remain opt-in. Use rebuild mode only when
intentionally replacing existing stats outputs:

```bash
qsub -v DAD_REBUILD=1 meg/cluster/submit_drsa_stats.sh
```

Restrict the stats job when inspecting a focused family or sensor group:

```bash
qsub -v DAD_DRSA_STATS_TARGETS=line-summaries,DAD_DRSA_STATS_ANALYSIS_FAMILIES=A_dRSA_alt,DAD_DRSA_STATS_SENSOR_GROUPS=occipital \
  meg/cluster/submit_drsa_stats.sh
```

Opt in to simulation stats only when that model-only control is specifically
needed:

```bash
qsub -v DAD_DRSA_STATS_ANALYSIS_FAMILIES=A_dRSA,A_dRSA_alt,simulations_model_vs_model \
  meg/cluster/submit_drsa_stats.sh
```

Disable the default deviant-minus-nondeviant paired sensitivity analysis:

```bash
qsub -v DAD_DRSA_STATS_PAIRED_CONTRASTS=none \
  meg/cluster/submit_drsa_stats.sh
```

Show every suprathreshold cluster in the HTML tables instead of only
significant clusters:

```bash
qsub -v DAD_DRSA_STATS_REPORT_CLUSTER_TABLES=all \
  meg/cluster/submit_drsa_stats.sh
```

Outputs are written below:

```text
.../derivatives/analysis/group/ses-01/drsa/stats/
.../derivatives/reports/group/ses-01/10_drsa_stats_report.html
.../derivatives/reports/group/ses-01/10_drsa_stats_report_embedded.html
```

When the MORWUR volume is mounted locally, inspect the report at:

```text
/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/10_drsa_stats_report.html
```

## Run Optional Dot-Side Decoding Sanity Check

This is an optional sanity-check workflow. It consumes completed Stage-04
epochs and the existing `source_manifest.json`; it does not rerun epoching,
dRSA export, MATLAB dRSA, or any preprocessing stage.

First validate a representative 20-identity and 40-identity subject without
writing outputs:

```bash
ssh -Y damiano.grignolio@hpc.dipsen.unitn.it
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
qsub -v DAD_VALIDATE_ONLY=1 -t 1,4 \
  meg/cluster/submit_dot_side_decoding.sh
```

Submit all subjects with the four standard analyses, 120 Hz decoding grid,
10-sample clipped sliding windows, no trial averaging beyond single trials,
5-fold stratified CV, and seed `42`:

1. `always_visible_observed`
2. `occluded_nondeviant_observed`
3. `occluded_deviant_observed`
4. `occluded_deviant_predicted`

```bash
qsub -t 1-7 meg/cluster/submit_dot_side_decoding.sh
```

Run only selected analyses by exporting the comma-separated specification
before `qsub`. Exporting avoids Grid Engine interpreting the comma as a
separator between environment-variable assignments:

```bash
export DAD_DOT_SIDE_ANALYSES="occluded_deviant_observed,occluded_deviant_predicted"
qsub -v DAD_DOT_SIDE_ANALYSES,DAD_DECODE_WINDOW_SAMPLES=15,DAD_TRIALS_PER_AVERAGE=2,DAD_RANDOM_SEED=97 \
  -t 1-7 meg/cluster/submit_dot_side_decoding.sh
```

Use one identity-level pseudo-trial per retained identity by setting the same
variable to the literal value `all_trials_per_identity`:

```bash
qsub -v DAD_REBUILD=1,DAD_TRIALS_PER_AVERAGE=all_trials_per_identity \
  -t 1-7 meg/cluster/submit_dot_side_decoding.sh
```

Retry a failed participant with guarded resume:

```bash
qsub -t 5 meg/cluster/submit_dot_side_decoding.sh
```

Intentionally replace existing sanity-check outputs:

```bash
qsub -v DAD_REBUILD=1 -t 1-7 \
  meg/cluster/submit_dot_side_decoding.sh
```

Resource settings can be overridden with standard `qsub` options when wider
windows or trial averaging increase runtime:

```bash
qsub -pe smp 8 -l h_vmem=32G -t 1-7 \
  meg/cluster/submit_dot_side_decoding.sh
```

Monitor jobs and logs:

```bash
qstat
tail -f cluster_logs/DAD_dot_side.JOB_ID.TASK_ID.log
```

Outputs are written under:

```text
.../derivatives/sanity_checks/sub-XX/ses-01/dot_side_decoding/<condition>/<model-source>/
```

Inspect, in order:

1. `<condition>_<model-source>_dot_side_decoding_summary.json`: Stage-04 handoff validation,
   identity range, pseudo-trial counts, dropped remainders, window policy,
   training/testing diagnostics, skipped-window reasons, and peak/mean accuracies.
2. `<condition>_<model-source>_dot_side_decoding_scores.csv`: one row per decoded center
   time, including start/end sample indices, start/end times, effective
   window width, class counts, skip reason, and accuracy values.
3. `<condition>_<model-source>_dot_side_decoding_3panel.png`: ERFs, ipsi-minus-contra, and
   decoding curves. Occluded-condition decoding panels mark the first complete
   occlusion and first reappearance frames from shared configuration.
4. `<condition>_<model-source>_dot_side_decoding_report.html`: embedded summary,
   result figure, training/testing assessment, and MNE layout figures.

### Build the group report

After all requested subject-level jobs finish, build the combined report as a
single non-array job. To rebuild subjects and automatically wait before the
group aggregation, submit both jobs with a Grid Engine dependency:

```bash
decode_job="$(qsub -terse -v DAD_REBUILD=1 -t 1-7 meg/cluster/submit_dot_side_decoding.sh)"
decode_job_id="${decode_job%%.*}"
qsub -hold_jid "$decode_job_id" -v DAD_REBUILD=1 \
  meg/cluster/submit_dot_side_group_report.sh
```

The combined report runs the same group analysis separately for all four
condition/model-source pairs. It uses the balanced accuracy column
`accuracy_balanced_undersampled_sets`, tests accuracy minus chance (`0.5`) with
a one-sample sign-flip cluster permutation over the time axis, uses
cluster-forming `p=0.05`, cluster alpha `0.05`, and a directional
greater-than-chance alternative. The group figure annotates the
cluster-forming p threshold and the equivalent t cutoff. Display bands are
95% confidence intervals by default.

Use SEM display bands instead:

```bash
qsub -v DAD_DOT_SIDE_INTERVAL=sem \
  meg/cluster/submit_dot_side_group_report.sh
```

Change the cluster-forming p threshold or run a two-sided test:

```bash
qsub -v DAD_CLUSTER_FORMING_P=0.01,DAD_CLUSTER_TAIL=two-sided \
  meg/cluster/submit_dot_side_group_report.sh
```

Smooth each subject's accuracy series before group statistics and plotting.
The default `0.0` disables smoothing; `1.0` applies a light Gaussian smoothing
kernel with sigma equal to one decoded sample:

```bash
qsub -v DAD_REBUILD=1,DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES=1.0 \
  meg/cluster/submit_dot_side_group_report.sh
```

Validate group inputs before writing outputs:

```bash
qsub -v DAD_VALIDATE_ONLY=1 \
  meg/cluster/submit_dot_side_group_report.sh
```

Intentionally replace an existing group report:

```bash
qsub -v DAD_REBUILD=1 \
  meg/cluster/submit_dot_side_group_report.sh
```

Group figures and the combined report are written under:

```text
.../derivatives/sanity_checks/group/ses-01/dot_side_decoding/<condition>/<model-source>/
.../derivatives/reports/group/ses-01/dot_side_decoding/combined/
```

The combined embedded report is:

```text
.../dot_side_decoding/combined/dot_side_all_analyses_group_report_embedded.html
```

It contains:

1. all four group mean analyses, individual traces, intervals, chance lines,
   cluster-permutation tables, and configured occlusion/reappearance markers
   on the three occluded analyses;
2. only after all group analyses, subjects 01-07 with each subject's four
   three-panel figures plus a three-condition source-trial count figure by
   trajectory identity before trial averaging.

### Run Optional Dot-Side Time-Generalization Check

The time-generalization companion trains the dot-side classifier at each
sliding-window center and tests it at every other center, producing a
train-time x test-time matrix for the same four analysis keys:

1. `always_visible_observed`
2. `occluded_nondeviant_observed`
3. `occluded_deviant_observed`
4. `occluded_deviant_predicted`

Validate representative inputs without writing outputs:

```bash
qsub -v DAD_VALIDATE_ONLY=1 -t 1,4 \
  meg/cluster/submit_dot_side_time_generalization.sh
```

Run all subjects with the default 120 Hz grid, 10-sample windows, 5-fold
stratified CV, and seed `42`:

```bash
qsub -t 1-7 meg/cluster/submit_dot_side_time_generalization.sh
```

Use identity-level pseudo-trials or selected analyses:

```bash
export DAD_DOT_SIDE_TG_ANALYSES="occluded_deviant_observed,occluded_deviant_predicted"
qsub -v DAD_DOT_SIDE_TG_ANALYSES,DAD_TRIALS_PER_AVERAGE=all_trials_per_identity \
  -t 1-7 meg/cluster/submit_dot_side_time_generalization.sh
```

The matrix decoder is heavier than the 1D line decoder because each valid
train window is evaluated against every valid test window. The default wrapper
requests 4 slots and 32G. If larger windows or single-trial decoding run
slowly, request more memory explicitly:

```bash
qsub -pe smp 8 -l h_vmem=48G -t 1-7 \
  meg/cluster/submit_dot_side_time_generalization.sh
```

Participant outputs are written under:

```text
.../derivatives/sanity_checks/sub-XX/ses-01/dot_side_time_generalization/<condition>/<model-source>/
```

Inspect, in order:

1. `<condition>_<model-source>_dot_side_time_generalization_summary.json`
2. `<condition>_<model-source>_dot_side_time_generalization_matrix.npz`
3. `<condition>_<model-source>_dot_side_time_generalization_matrix.csv`
4. `<condition>_<model-source>_dot_side_time_generalization_heatmap.png`
5. `<condition>_<model-source>_dot_side_time_generalization_report.html`

Build the group time-generalization report after the participant array
finishes. With a dependency:

```bash
tg_job="$(qsub -terse -v DAD_REBUILD=1 -t 1-7 meg/cluster/submit_dot_side_time_generalization.sh)"
tg_job_id="${tg_job%%.*}"
qsub -hold_jid "$tg_job_id" -v DAD_REBUILD=1 \
  meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

The group report uses `accuracy_balanced_undersampled_sets` by default and
tests `accuracy - 0.5` with a 2D sign-flip cluster permutation over the
train-time x test-time lattice. A group cell is included only when every
included subject has a finite value, so clusters cannot bridge NaN cells from
skipped participant train/test windows. Significant 2D clusters are shaded and
outlined in red on the group mean matrix.

The detailed 2D cluster table is omitted from the HTML report by default so
the group report stays compact. Regenerate with the explicit opt-in variable
when the per-cluster train/test ranges, masses, and p-values should be shown
inline:

```bash
qsub -v DAD_REBUILD=1,DAD_DOT_SIDE_TG_INCLUDE_CLUSTER_TABLE=1 \
  meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

Change the cluster-forming threshold or run a two-sided test:

```bash
qsub -v DAD_CLUSTER_FORMING_P=0.01,DAD_CLUSTER_TAIL=two-sided \
  meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

Smooth each subject's matrix before group statistics and plotting. The
time-generalization wrapper accepts `DAD_DOT_SIDE_TG_SMOOTH_SIGMA_SAMPLES` and
falls back to the shared `DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES` if the
time-generalization-specific variable is unset:

```bash
qsub -v DAD_REBUILD=1,DAD_DOT_SIDE_TG_SMOOTH_SIGMA_SAMPLES=1.0 \
  meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

Group outputs are written under:

```text
.../derivatives/sanity_checks/group/ses-01/dot_side_time_generalization/<condition>/<model-source>/
.../derivatives/reports/group/ses-01/dot_side_time_generalization/combined/
```

The combined embedded report is:

```text
.../dot_side_time_generalization/combined/dot_side_time_generalization_all_analyses_group_report_embedded.html
```

### Local dot-side commands

Run all subject-level dot-side outputs locally from the laptop:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
analyses=(
  "always_visible_observed"
  "occluded_nondeviant_observed"
  "occluded_deviant_observed"
  "occluded_deviant_predicted"
)
for sub in 01 02 03 04 05 06 07; do
  for analysis in "${analyses[@]}"; do
    condition="${analysis%_*}"
    model_source="${analysis##*_}"
    python meg/code/sanity_checks/01_dot_side_decoding.py \
      --config meg/params/pipeline_default.yaml \
      --subject-config "meg/params/sub-${sub}.yaml" \
      --subject "$sub" \
      --condition "$condition" \
      --model-source "$model_source" \
      --decode-window-samples 10 \
      --trials-per-average 1 \
      --random-seed 42 \
      --overwrite
  done
done
```

Build the local group report:

```bash
python meg/code/sanity_checks/02_dot_side_group_report.py \
  --config meg/params/pipeline_default.yaml \
  --analyses all \
  --group-subjects 01 02 03 04 05 06 07 \
  --metric accuracy_balanced_undersampled_sets \
  --interval ci95 \
  --cluster-forming-p 0.05 \
  --cluster-alpha 0.05 \
  --tail greater \
  --overwrite \
  --print-summary
```

Use `--interval sem` to draw SEM bands instead of the default 95% CI.

Open the portable combined report from the laptop after generation:

```bash
open /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/dot_side_decoding/combined/dot_side_all_analyses_group_report_embedded.html
```
