# Poster analysis: settings, readiness and cluster commands

Updated 2026-11-17 after user-approved consolidation of P21/P23/P24 into
`noica-20261117-v2`. All 30 participants now use the same run folder and
ordinary participant selectors. No analysis job was submitted by Codex.

## Current cohort and inputs

P04-P39 excluding P10, P15, P29, P30, P34 and P35: **30 participants**.
P01-P03 remain historical pilots with different experiment versions.

The completed recovery epochs and exports for P21, P23 and P24 were moved into
`processing_runs/noica-20261117-v2/analysis/sub-XX` and matching `reports/sub-XX`
folders. Large FIF, CSV, MAT and sidecar files were not rewritten or duplicated.
Their saved epoch/export scientific contracts remain unchanged. Report paths,
the storage label and ordinary participant selectors were updated and validated.
The local and shared cluster selectors already match; no deployment is required
for this consolidation. Historical recovery selectors remain available.

The participant settings used to create these recovered epochs remain explicit:

| Participant | Executed recovery setting now used by the ordinary selector |
| --- | --- |
| P21 | Two specifically approved blink-annotation end trims before epoching. |
| P23 | Amplitude screen disabled; annotation and trial-policy rejection remain enabled. |
| P24 | EDF-derived six-channel eye reconstruction for all five acquisitions, with consistent channel ordering. |
| P38 | Original reviewed exact-source filtering exception, loaded from upstream tuning. |

Earlier incomplete main-run outputs and original recovery reports/selectors are
archived at:
`/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/previous_runs/cohort-consolidation-20261117-v1/`.
Its `consolidation.json` records all moves and before/after validation. The
migration procedure is in `meg/sandbox/consolidate_recovery_20261117/`.

## Existing job outcomes and rerun

| Job | Outcome before consolidation |
| --- | --- |
| 11869 | 27 participant dRSA jobs completed; P21/P23/P24 failed because the old main folders lacked their recovered inputs. |
| 11870 | Group average: stuck at `DependencyNeverSatisfied`. |
| 11871 | Group statistics: stuck at `DependencyNeverSatisfied`. |
| 11872 | ERF: failed because the original P24 epochs were absent. |
| 11873 | Peak dRSA: stuck at `DependencyNeverSatisfied`. |

The 27 completed participants have 108 verified finite 320 x 320 dRSA matrices
and 324 matching output bindings. Their raw actual/predicted matrices match
exactly through model-time frame 130. The rerun uses `DAD_REBUILD=0`: completed
participants are reused and P21/P23/P24 are computed from the recovered inputs.
The full 30-participant array provides a new common dependency for group jobs.
ERF starts independently because all required epochs are now in the main run.

Run the cancellation command below once to remove the three obsolete waiting
jobs. Then copy the entire submission block from `cd` to `--print-summary'`.
The ERF baseline is the user's final choice **-0.5 to 0 seconds**.

## Recommended initial settings

| Item | Setting | Alternatives |
| --- | --- | --- |
| MATLAB analysis | `alt` = `A_dRSA_alt.m` | `real` = `A_dRSA.m`, for the broader multi-model position/direction analysis |
| dRSA measure | `corr`: Pearson correlation of neural and position RDMs | `PCR` is a different analysis and is not the poster replication |
| Sensors | Occipital; current alt code uses good magnetometers and gradiometers together | Other lobes, whole brain, or a separately implemented sensor-type comparison |
| Models | Observed in all three conditions; predicted in occluded-deviant only | Add `prediction_error` or registered `pe_*` sources later |
| Neural sampling/smoothing | Existing 120 Hz exports; Gaussian neural smoothing over 15 samples, as in the pilot | Changing this requires separate outputs and changes the analysis |
| Matrix display | Raw plus Gaussian display smoothing, window 9 samples | Raw MAT matrices remain available |
| Peak curve | Raw saved matrices; 25 ms brain windows, 25 ms steps, lag search +/-1 s; minimum 3 valid samples per lag; participant peaks first, then group mean +/- SEM | Also generate 200 ms windows as a separately named sensitivity figure |
| ERF | Occipital magnetometers, three source conditions; RMS/GFP overlay; baseline -0.5 to 0 s; analysis -0.5 to 2.667 s, 200 Hz | Gradiometers or other sensor groups |
| Statistics | 10,000 permutations, seed 42; two-sided, forming p=.05, cluster alpha=.05 | Omit peak cluster bars for a descriptive plot |

The peak-window history is inconsistent: the supplied poster caption states
25 ms, but the archived SVG with permutation bars records 200 ms. The archive
also contains descriptive 25/50/100/200 ms variants. The commands below expose
the window explicitly rather than claiming an exact reconstruction of every
manually assembled poster element. Peak cluster statistics use Gaussian sigma
1 in units of the 25 ms output grid; plotted means are not smoothed by that
option. Each of the six pairwise contrasts is corrected over time separately,
not jointly over all six contrasts. Stage-10 matrix tests are separate from
these peak-curve comparisons.

These analyses consume the existing epochs/exports. ICA stays off and epoch
timing, annotations, epoch baseline/detrending and rejection settings are
preserved. The user-selected ERF baseline is now -0.5 to 0 s, replacing the
poster's -0.3 to -0.1 s interval. It is applied during ERF analysis and fits
inside the saved epochs; no epoching/export rerun is needed for this change.

## Actual versus predicted equality

Read-only checks on all 30 actual trajectory files found exactly zero
coordinate difference and zero Euclidean position-RDM difference through
frame 130. The first different frame is 131 for every participant.
At 120 Hz, the exported frame times are `(frame-1)/120`, so frame 130 is
1.075 s and frame 131 is 1.08333 s. Reappearance frame 191 is 1.58333 s.
The poster uses rounded frame/120 descriptions; keep the established export
time axis and handoff timing.

For the same occluded-deviant neural data, raw correlation dRSA must therefore
match in the first 130 **model-time rows**, across all neural times. Archived
pilot P01, P04 and P07 matrices satisfy that exactly (maximum difference 0).
Their first 130 neural-time columns differ, which is allowed because those
columns include post-deviance model-time rows.

The poster brain-time peak searches +/-1 s of lag. Thus an early brain time
can match a later model time where actual and predicted paths already differ.
The peak curves need not coincide before occlusion. To test the shared-prefix
invariant, compare raw matched model-time rows, not maxima across all lags.
Display smoothing can also blend rows across the divergence boundary. No
change to the dRSA computation was justified by these equality checks.

## Deployment status

The analysis launchers were already deployed for jobs 11869/11872. The three
updated canonical participant selectors have also been written and checked in
the shared cluster project. The submission block below uses those selectors;
it unsets an old recovery-selector override if still present in the terminal.
No new Python/MATLAB analysis code is required for this consolidation.
After the move, all 30 ordinary-selector dRSA input contracts and all 30 ERF
headers passed. The 27 existing dRSA resume contracts also match the deployed
cluster settings; their 108 matrices and 324 output bindings remain valid.

```bash
scancel 11870 11871 11873
```

## Full-cohort submission after consolidation

Run this block on the cluster login node, in Bash. It submits participant dRSA
as a parallel array. Group averaging, matrix statistics and peak plotting wait
for the entire array to succeed. ERF aggregation is independent once all epochs
exist. The settings deliberately override the wrappers' broader defaults.
Every job explicitly requests 4 CPUs per task. MATLAB/numerical-library thread
limits follow that allocation; the peak job sets matching limits explicitly.
Memory requests retain the existing wrapper defaults (24G participant dRSA,
12G group averaging, 16G statistics/ERF) and 8G for peak plotting. Group stages
join participant results; extra CPUs do not make every aggregation step parallel.

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
unset DAD_SUBJECT_CONFIG_DIR
export DAD_PROJECT_ROOT="$PWD"
export DAD_DERIVATIVES_ROOT=/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export DAD_EPOCH_PREPARATION=inline_no_ica
export DAD_PROCESSING_RUN=noica-20261117-v2
export DAD_REBUILD=0 DAD_ALLOW_MISSING=0 DAD_VALIDATE_ONLY=0
export DAD_POSTER_SUBJECTS="04 05 06 07 08 09 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 31 32 33 36 37 38 39"
export DAD_DRSA_GROUP_SUBJECTS="${DAD_POSTER_SUBJECTS// /,}"
export DAD_DOT_ERF_GROUP_SUBJECTS="$DAD_DRSA_GROUP_SUBJECTS"
export DAD_DRSA_STATS_SUBJECTS="$DAD_DRSA_GROUP_SUBJECTS"
export DAD_DRSA_ANALYSES=alt
export DAD_DRSA_CONDITIONS=always_visible,occluded_nondeviant,occluded_deviant
export DAD_DRSA_SENSOR_GROUPS=occipital
export DAD_DRSA_MODEL_SOURCES=observed,predicted
export DAD_DRSA_SMOOTHING_PLOTS=1
export DAD_DRSA_PLOT_SMOOTHING_WINDOW=9
export DAD_DRSA_GROUP_SMOOTHING_PLOTS=1
export DAD_DRSA_GROUP_PLOT_SMOOTHING_WINDOW=9

drsa_job=$(sbatch --parsable --cpus-per-task=4 --array=4-9,11-14,16-28,31-33,36-39 \
    meg/cluster/submit_drsa_analysis.sh)
sbatch --cpus-per-task=4 --dependency=afterok:"$drsa_job" meg/cluster/submit_drsa_group_average.sh

export DAD_DRSA_STATS_TARGETS=matrices
export DAD_DRSA_STATS_ANALYSIS_FAMILIES=A_dRSA_alt
export DAD_DRSA_STATS_SENSOR_GROUPS=occipital
export DAD_DRSA_STATS_MODEL_SOURCES=observed,predicted
export DAD_DRSA_STATS_N_PERMUTATIONS=10000
export DAD_DRSA_STATS_SEED=42
export DAD_DRSA_STATS_CLUSTER_ALPHA=0.05
export DAD_DRSA_STATS_SIGNIFICANCE_ALPHA=0.05
export DAD_DRSA_STATS_TAIL=two-sided
export DAD_DRSA_STATS_INCLUDE_SMOOTHED=0
export DAD_DRSA_STATS_PAIRED_CONTRASTS=none
sbatch --cpus-per-task=4 --dependency=afterok:"$drsa_job" meg/cluster/submit_drsa_stats.sh

export DAD_DOT_ERF_SENSOR_GROUPS=occipital
export DAD_DOT_ERF_SENSOR_TYPE=mag
export DAD_DOT_ERF_BASELINE=-0.500,-0.000
export DAD_DOT_ERF_TMIN=-0.5 DAD_DOT_ERF_TMAX=2.667
export DAD_DOT_ERF_ANALYSIS_SFREQ=200
export DAD_DOT_ERF_GROUP_OVERLAY_BAND=sem
export DAD_N_PERMUTATIONS=10000 DAD_RANDOM_SEED=42
export DAD_CLUSTER_FORMING_P=0.05 DAD_CLUSTER_ALPHA=0.05
export DAD_CLUSTER_TAIL=two-sided
sbatch --cpus-per-task=4 --array=1 meg/cluster/submit_dot_erf_analysis.sh

export DAD_POSTER_PEAK_WINDOW=0.025
sbatch --dependency=afterok:"$drsa_job" --job-name=DAD_drsa_peak \
    --cpus-per-task=4 --mem=8G --time=02:00:00 --export=ALL \
    --output=cluster_logs/%x.%j.log --error=cluster_logs/%x.%j.err \
    --wrap='export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4; exec "$HOME/bin/micromamba" run -n dad-meg python "$DAD_DERIVATIVES_ROOT/pilot/poster_figures/export_drsa_poster_figures.py" --project-root "$DAD_PROJECT_ROOT" --derivatives-root "$DAD_DERIVATIVES_ROOT" --site-config "$DAD_PROJECT_ROOT/meg/params/pipeline_cluster.yaml" --subjects $DAD_POSTER_SUBJECTS --figures brain-time-peak --sensor-group occipital --brain-window-sec "$DAD_POSTER_PEAK_WINDOW" --brain-step-sec 0.025 --lag-window-sec 1 --min-valid-samples 3 --band sem --cluster-permutation --cluster-n-permutations 10000 --random-seed 42 --output-dir "$DAD_DERIVATIVES_ROOT/processing_runs/$DAD_PROCESSING_RUN/poster_figures/cohort30/brain_time_peak" --formats png svg pdf --print-summary'
```

To reproduce the archived statistical SVG's window instead, set
`DAD_POSTER_PEAK_WINDOW=0.200` before submitting the peak job. Its output name
contains the window, so both versions can coexist. Omit `--cluster-permutation`
for a descriptive plot. Do not change the processing-run label to rename a
plot: that would select different input epochs/exports.

The array and all three group subject lists cover the same 30 participants,
including P21, P23, P24 and P38. No reduced-cohort submission is used here.

## Output locations on the Mac

All current outputs are under:

`/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/processing_runs/noica-20261117-v2/`

| Output | Path below that root |
| --- | --- |
| Subject neural MAT exports | `analysis/sub-XX/ses-01/drsa/*desc-identityAveraged_neural.mat` |
| Subject dRSA matrices and results | `analysis/sub-XX/ses-01/drsa/matlab_real_neural/A_dRSA_alt/` |
| Four group position dRSA matrices | `analysis/group/ses-01/drsa/group_average/matlab_real_neural/A_dRSA_alt/` |
| Group dRSA report | `reports/group/ses-01/07_matlab_drsa_group_report_smoothed_embedded.html` |
| Matrix statistical figures/masks | `analysis/group/ses-01/drsa/stats/` |
| Matrix statistics report | `reports/group/ses-01/10_drsa_stats_report_embedded.html` |
| ERF condition overlay PNG | `sanity_checks/group/ses-01/dot_erf/dot_erf_mag_occipital_condition_overlay_gfp_rms.png` |
| ERF arrays/contrasts | `sanity_checks/group/ses-01/dot_erf/` |
| ERF report | `reports/group/ses-01/dot_erf/dot_erf_mag_occipital_group_report_embedded.html` |
| Brain-time peak PNG/SVG/PDF | `poster_figures/cohort30/brain_time_peak/` |

The four matrix branches are `always_visible/corr`,
`occluded_nondeviant/corr`, `occluded_deviant/observed/corr`, and
`occluded_deviant/predicted/corr`; look inside each `matrices` directory.
The peak script is the original archived pilot exporter, explicitly redirected
to current inputs and a new output directory. Stage 09 itself makes
model-time-window summaries; it does not need to be run first for this plot.

On the cluster replace `/Volumes/MORWUR` with `/mnt/storage/tier2/morwur`.
Poster layout/colors/captions were edited manually; these commands reproduce
the analysis families and generate new figures, not the assembled poster PDF.
