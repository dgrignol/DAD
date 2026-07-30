# 30 Jul 2026 18:02

## Add v21 photodiode experiment runtime

### `experiment/CreateInputFiles_v23_rescueTraject.m`
- Preserve `launcherCfgPath` across the schedule script's workspace cleanup so `runExperiment.m` can reload its temporary launcher snapshot before starting the runtime.
- Update the active message-flow dependency documentation from the retired v20 runtime to the v21 photodiode runtime.

### `experiment/MoveDot1_experiment_occlusion_v20_rescueTraject.m` -> `experiment/MoveDot1_experiment_occlusion_v21_photodiode.m`
- Rename the active runtime and its output/debug identities from v20 rescueTraject to v21 photodiode while retaining the established schedule, EyeLink, replay, ITI, and display behavior.
- Move condition-trigger emission for codes `71`, `81`, and `91` to the fully drawn first-frame presentation boundary immediately before the photodiode patch draw and display flip.
- Add validated, configurable patch enablement, size, top-left-relative position, high/low colors, four-frame trial pulses, and two-frame marked-event pulses for trigger codes `111`, `112`, `114`, `115`, `116`, `117`, and `118`.
- Route experiment-controlled screen flips through a photodiode-aware wrapper, keep the patch low outside requested pulses, prevent overlapping requests from collapsing the required low transition, and force a final low frame during cleanup.
- Extend runtime configuration persistence, diagnostics, output naming, and catch-question presentation so trigger logging and photodiode frames remain aligned.

### `experiment/lib/Config_runtime_v20_rescueTraject.m` -> `experiment/lib/Config_runtime_v21_photodiode.m`
- Rename the development runtime preset for v21 and expose the patch geometry, colors, enablement policy, and trial/event pulse widths as documented workspace overrides.
- Keep the non-lab EyeLink and fixation controls delegated to schedule defaults while printing the active photodiode settings in the preset audit line.

### `experiment/lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m`
- Update active runtime and eye-tracker comments to identify the v21 photodiode consumer while preserving the existing v23 schedule configuration and compatibility-oriented filename.

### `experiment/runExpReadable.m`
- Switch both runtime-preset loads and the final direct launch from v20 to v21.
- Update documented inputs and output filenames so the readable launcher matches the active photodiode stack.

### `experiment/runExperiment.m`
- Default to the v21 config and runtime, validate an overridable `runtimeExperimentScript`, and retain its path in the temporary launcher snapshot.
- Carry all photodiode patch and pulse overrides across generator/schedule scripts before launching the selected runtime.
- Update launcher documentation and output naming for the V30/v23/v21 stack.

### `experiment/runTestExperiment.m`
- Point the laptop-safe end-to-end launcher and optional output cleanup at the v21 runtime/output identity.
- Update the documented test stack while retaining fake-eye-tracker and non-MEG defaults.

### `experiment/run_MEG_lab/CreateInputFiles_v23_rescueTraject.m`
- Mirror the active schedule dependency documentation and preserve `launcherCfgPath` through workspace cleanup so the standalone lab launcher can continue after schedule generation.

### `experiment/run_MEG_lab/MoveDot1_experiment_occlusion_v20_rescueTraject.m` -> `experiment/run_MEG_lab/MoveDot1_experiment_occlusion_v21_photodiode.m`
- Mirror the v21 photodiode runtime rename and timing implementation in the standalone MEG-lab deployment copy.
- Keep the deployed runtime byte-identical to the development runtime so condition/event pulses, low-state cleanup, and trigger timing follow the same code path in both environments.

### `experiment/run_MEG_lab/lib/Config_runtime_v20_rescueTraject.m` -> `experiment/run_MEG_lab/lib/Config_runtime_v21_photodiode.m`
- Mirror the v21 patch and pulse controls in the lab preset.
- Retain lab-specific real-EyeLink defaults, the `3.0 deg` fixation window, bounded replay policy, and calibration behavior.

### `experiment/run_MEG_lab/lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m`
- Update active runtime and eye-tracker comments to the v21 photodiode name while retaining the mirrored v23 schedule defaults.

### `experiment/run_MEG_lab/runExpReadable.m`
- Switch the standalone lab launcher to the v21 runtime/config/output identity.
- Preserve the lab-specific readable-launch defaults, including the documented `3.0 deg` fixation preset.

### `experiment/run_MEG_lab/runExperiment.m`
- Mirror the configurable v21 runtime selection and photodiode override handoff used by the development launcher.
- Preserve the launcher snapshot across generation and schedule steps before starting acquisition.

### `experiment/run_MEG_lab/runTestExperiment.m`
- Mirror the v21 test-launch and output-cleanup paths in the standalone lab copy.

### `experiment/run_MEG_lab/stimuli_generation_V30_rescueTraject.m`
- Preserve `launcherCfgPath` during generator workspace cleanup so the lab launcher can reload its configuration for schedule generation and runtime launch.
- Update the active runtime compatibility comment to v21.

### `experiment/stimuli_generation_V30_rescueTraject.m`
- Preserve `launcherCfgPath` during generator workspace cleanup for the multi-step launcher.
- Update exported trigger/event metadata documentation to reference the active v21 runtime.

Suggested commit message:
`Add v21 photodiode experiment runtime`

# 29 Jun 2026 12:09

## Add dot ERF sanity-check analysis

### `meg/code/sanity_checks/05_dot_erf_analysis.py`
- Add a WORK_IN_PROGRESS dot-condition ERF sanity-check workflow that consumes completed Stage-04 epoch FIF files and source manifests without modifying decoding or dRSA outputs.
- Support read-only validation and full report generation for requested group subjects, one MEG sensor type per run, and one or more Stage-07 MEGIN sensor groups filtered to common good channels.
- Apply baseline correction inside the script, crop/resample the analysis window, compute per-subject condition ERFs, and construct `occluded_combined` as an equal-weight average of nondeviant and deviant ERFs rather than a pooled trial average.
- Resolve guarded subject and group output paths for ERF arrays, trial-count diagnostics, imbalance diagnostics, cluster-stat arrays/JSON, linked HTML reports, embedded HTML reports, and diagnostic figures.
- Add paired ERF contrasts for `occluded_deviant - occluded_nondeviant` and `occluded_combined - always_visible`.
- Run MNE spatio-temporal one-sample cluster tests over subject-level channel x time contrasts, plus a signed channel-average 1-D temporal cluster test for each contrast.
- Add optional trial-count equalized sensitivity analyses when within-subject source-count imbalance exceeds the configured threshold.
- Generate group-first figures for selected sensor layout, condition butterfly/topomap plots, source-condition RMS/GFP overlays, trial counts, contrast butterflies, cluster heatmaps, and channel-average cluster highlights.
- Write report-ready metadata that records units, scaling, baseline, occlusion timing, sensor-selection policy, cluster parameters, topomap latencies, imbalance policy, skipped inputs, and the WIP interpretation warning.

Suggested commit message:
`Add dot ERF sanity-check analysis`

# 29 Jun 2026 12:08

## Add Stage-10 dRSA group statistics workflow

### `meg/code/preprocessing_basic/10_drsa_stats.py`
- Add an optional Stage-10 group statistics CLI that reads existing Stage-07 dRSA result MAT files without rerunning MATLAB, epoching, export, or preprocessing.
- Assemble matched subject stacks for 2-D matrix tests and Stage-08-compatible 1-D lag-summary tests, preserving the convention `lag = neural time - model time`.
- Add raw and explicitly labeled `smoothed_display` input variants, with default model selection focused on position-like and prediction-error matrices from `A_dRSA` and `A_dRSA_alt`.
- Implement one-sample sign-flip cluster permutation tests with an automatic MNE backend for fully finite arrays and a NumPy/SciPy fallback for NaN-containing support.
- Add the default paired actual-position contrast `occluded_deviant - occluded_nondeviant`, first differencing matrices within subject before matrix or lag-summary aggregation.
- Write guarded statistical outputs below the group derivative tree, including NPZ arrays, per-test JSON, cluster tables, summary JSON/CSV, linked HTML, and embedded HTML reports.
- Include report organization by input variant, condition/contrast, target, sensor group, model source, and analysis family, with significant-only cluster tables by default.

### `meg/tests/test_drsa_stats.py`
- Add synthetic unit coverage for positive 2-D and 1-D cluster detection with the NumPy backend.
- Verify paired matrix preparation drops subjects missing one side of the deviant/nondeviant contrast.
- Verify paired line summaries are computed from the within-subject difference matrix before lag averaging, preserving asymmetric NaN support.
- Cover MNE cluster-mask normalization, automatic fallback to NumPy for zero-variance features, CLI defaults, significant-only report tables, and safe handling of NaNs.
- Add an end-to-end CLI test on a fake derivative tree that writes stats outputs, linked report HTML, and embedded report images without MATLAB.

### `meg/cluster/submit_drsa_stats.sh`
- Add a one-task SGE wrapper for the Stage-10 dRSA stats report using existing Stage-07 result MAT files.
- Expose cluster controls for subjects, conditions, model selectors, model sources, sensor groups, analysis families, matrix/line targets, line intervals, lag window, lag mode, smoothed-display inclusion, permutation count, cluster thresholds, tail, backend, paired contrasts, and cluster-table display.
- Activate the shared micromamba environment, configure headless plotting/cache paths, validate required deployment files, print core dependency availability, and run the Stage-10 Python script in resume or rebuild mode.

### `meg/cluster/deploy_review_to_cluster.sh`
- Extend the deployment payload to include the Stage-09 lag-surface script, the new Stage-10 stats script, the Stage-10 wrapper, related docs, and dRSA post-processing tests.
- Update the script header to describe optional dRSA post-processing as part of the cluster payload.

### `meg/cluster/README.md`
- Register `submit_drsa_stats.sh` in the available cluster wrapper list.
- Clarify that deployment now includes optional dRSA post-processing code in addition to sanity-check code.
- Add a runbook section for optional Stage-08 line summaries, Stage-09 lag-surface reports, and Stage-10 matrix/line-summary statistics.
- Document Stage-10 defaults, focused qsub overrides, simulation-family opt-in, paired-contrast disabling, cluster-table verbosity, and expected report/output locations.

Suggested commit message:
`Add Stage-10 dRSA group statistics workflow`

# 26 Jun 2026 19:09

## Add SAMBA poster results reasoning note

### `paper/samba_poster_results_reasoning.md`
- Add a decision document for the SAMBA poster that argues for a cautious framing around dynamic updating of moving-object position representations after occlusion.
- Summarize the current pilot evidence from occipital dRSA lag summaries, dot-side decoding, and time-generalization sanity checks, emphasizing actual-trajectory information after reappearance over a strong prediction-error claim.
- Record the local methods, pipeline, report, result-presentation, and literature sources used to ground the poster reasoning.
- Define poster-ready framing elements: title directions, intro logic, research question, working hypotheses, result-panel priorities, discussion language, and a claim hierarchy.
- Identify current limitations, including small `N=7`, pilot-stack heterogeneity, display-only smoothing, descriptive dRSA interpretation, and the need for subject-level contrasts.
- Add ready-to-run follow-up prompts for direct dRSA contrasts, group ERF reappearance contrasts, direction/curvature dRSA, actual-vs-predicted dot-side decoding contrasts, prediction-error model comparison, and robustness/confound checks.
- End with a poster-safe bottom-line claim that dynamic position coding appears to align rapidly with the actually reappearing trajectory while direct prediction-error evidence remains to be quantified.

Suggested commit message:
`Add SAMBA poster results reasoning note`

# 26 Jun 2026 19:08

## Add dRSA sliding lag surface summaries

### `meg/code/preprocessing_basic/09_plot_drsa_sliding_lag_surfaces.py`
- Add an optional Stage-09 post-processing script that reads existing Stage-07 dRSA result MAT files without rerunning MATLAB or recomputing dRSA.
- Compute lag-by-model-time surfaces with the Stage-08 convention `lag = neural time - model time`, using model-anchored windows that restrict model-time rows while allowing neural-time samples anywhere inside matrix bounds.
- Target the initial occipital `A_dRSA_alt` position summaries for actual visible, actual nondeviant, actual deviant, expected/predicted deviant, and prediction-error deviant models.
- Write static heatmaps, standalone/in-report Plotly 3-D waterfall HTML, and two peak-timecourse comparison figures with shared report value limits.
- Support subject and group scopes, explicit subject lists, lag/window/step controls, minimum finite sample counts, shared value-limit modes, peak-latency masking, `--resume`, `--overwrite`, and `--print-summary`.
- Guard output paths against raw-tree writes and generate linked plus embedded HTML reports with occlusion markers, skipped-input diagnostics, and interpretation-risk notes.

### `meg/docs/DRSA_SLIDING_LAG_SURFACES.md`
- Add an interpretation guide for the exploratory Stage-09 lag surfaces, including the matrix convention, five initial target summaries, heatmap/waterfall outputs, and peak-timecourse extraction.
- Document descriptive-use caveats around overlapping windows, temporal autocorrelation, edge support, prediction-error undefined support, model-anchored sampling, and peak-lag instability.

### `meg/tests/test_drsa_sliding_lag_surfaces.py`
- Add unit coverage that verifies model-anchored windows average diagonal lag samples while permitting neural columns outside the local model-time window.
- Test peak magnitude/lag extraction and relative peak-latency threshold calculation on synthetic surfaces with missing rows.
- Verify target selection returns exactly the five intended occipital position summaries.
- Add an end-to-end CLI test with synthetic Stage-07-like MAT files that builds a group report, five heatmaps, five Plotly waterfalls, two peak-timecourse PNGs, and embedded HTML images.

Suggested commit message:
`Add dRSA sliding lag surface summaries`

# 26 Jun 2026 19:07

## Add cluster report controls and dRSA group-average wrapper

### `meg/cluster/README.md`
- Register `submit_drsa_group_average.sh` in the cluster script list as the optional Stage-07 numeric group-average report job.
- Replace the manual group-average Python command with a `qsub -hold_jid` workflow for running the report after Stage-07 subject jobs finish.
- Document the group-average wrapper's config, site-config, representative subject-config, subject-list, smoothing, and smoothing-window environment overrides.
- Add the explicit opt-in command for including detailed 2D dot-side time-generalization cluster tables in the HTML report.

### `meg/cluster/deploy_review_to_cluster.sh`
- Include `meg/cluster/submit_drsa_group_average.sh` in the minimal cluster deployment payload so group-average aggregation is available after upload.

### `meg/cluster/submit_dot_side_time_generalization_group_report.sh`
- Add `DAD_DOT_SIDE_TG_INCLUDE_CLUSTER_TABLE` as a validated `0`/`1` wrapper control.
- Forward `--include-cluster-table` to the Python group-report script only when the opt-in variable is enabled.
- Update submission examples, environment-variable documentation, and runtime logging to show whether detailed 2D cluster tables are included.

### `meg/cluster/submit_drsa_group_average.sh`
- Add a one-task SGE wrapper that builds the optional Stage-07 numeric group-average dRSA report from existing subject-level MATLAB result MAT files.
- Activate the shared micromamba environment, configure headless plotting/cache paths, validate required configs and Python helpers, and print MAT-reader versions before aggregation.
- Run `07_run_matlab_drsa.py --group-average-only --group-subjects ... --print-summary`, with documented controls for subject selection, smoothed group figures, and plot smoothing window.

### `meg/code/sanity_checks/04_dot_side_time_generalization_group_report.py`
- Add `--include-cluster-table` to keep the detailed per-cluster HTML table opt-in while preserving cluster overlays and JSON metadata by default.
- Record `include_cluster_table` in validate-only and completed machine reports.
- Include a summary line that states whether the cluster table was included or omitted in the generated HTML.

### `meg/docs/PIPELINE_SPEC_basic.md`
- Document that detailed 2D cluster tables are omitted from HTML reports by default and can be restored with the script flag or cluster wrapper environment variable.

### `meg/tests/test_cluster_wrappers.py`
- Add the new group-average wrapper to shell syntax coverage and deployment-payload checks.
- Verify the dot-side time-generalization group wrapper exposes the cluster-table opt-in control.
- Add coverage that the Stage-07 group-average wrapper exposes config, subject, smoothing, `--group-average-only`, `--group-subjects`, and `--print-summary` controls.

### `meg/tests/test_dot_side_time_generalization.py`
- Update the group-report HTML test to assert detailed 2D cluster tables are omitted by default.
- Add an explicit `--include-cluster-table` run that verifies the machine-report flag and embedded HTML table heading are present when requested.

Suggested commit message:
`Add cluster report controls and dRSA group-average wrapper`

# 26 Jun 2026 19:06

## Preserve undefined dRSA smoothing support

### `meg/code/analysis/drsa_stage_utils.m`
- Keep a finite-support mask before MATLAB display smoothing and restore original non-finite cells to `NaN` afterward.
- Document that undefined cells define analysis support, so smoothing must not make impossible model/neural regions visible in exported matrices.

### `meg/code/preprocessing_basic/07_run_matlab_drsa.py`
- Replace the local Python group-display smoothing implementation with the shared `drsa_result_utils.smooth_matrix_for_plot` helper.
- Keep the Stage-07 wrapper function so group-average plotting continues to use the same call site while inheriting the shared NaN-support policy.

### `meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py`
- Clarify in the script header that `--include-smoothed` is display-only and cannot create values outside a matrix's originally valid support.

### `meg/code/preprocessing_basic/drsa_result_utils.py`
- Extend the shared smoothing helper so cells that were non-finite before smoothing remain `NaN` after normalized Gaussian smoothing.
- Update the module and function docs to state that prediction-error masks and other undefined supports are preserved for display summaries and group averages.

### `meg/tests/test_drsa_line_summaries.py`
- Add a regression test for prediction-error matrices whose pre-occlusion model-time rows are undefined.
- Verify that smoothing preserves the empty pre-occlusion support and that a model-anchored lag curve over that interval remains all `NaN`.

### `meg/tests/test_matlab_drsa_stage.py`
- Add a Stage-07 regression test that exercises the Python group-display smoother on a matrix with undefined top rows.
- Assert that smoothed group-display matrices keep the original undefined rows non-finite while retaining finite values in valid rows.

Suggested commit message:
`Preserve undefined dRSA support during smoothing`

# 24 Jun 2026 20:48

## Add group smoothing for dot-side cluster reports

### `meg/code/sanity_checks/02_dot_side_group_report.py`
- Add `--smooth-sigma-samples` to optionally smooth each subject's finite 1D accuracy segments before group statistics and plotting.
- Preserve skipped-window NaNs as hard boundaries by smoothing only within contiguous finite subject segments.
- Record smoothing settings, raw subject values, and tested/plotted subject values in the group report payload and diagnostics.

### `meg/code/sanity_checks/04_dot_side_time_generalization_group_report.py`
- Add `--smooth-sigma-samples` to optionally smooth each subject's train-time x test-time matrix before 2D cluster permutation and group plotting.
- Smooth only within finite 4-neighbor connected components using normalized Gaussian convolution so skipped cells cannot be bridged by smoothing.
- Record smoothing settings, raw matrices, and tested/plotted matrices in the group report payload and diagnostics.

### `meg/cluster/submit_dot_side_group_report.sh`
- Expose 1D group smoothing through `DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES` and pass it to the Python group-report script.

### `meg/cluster/submit_dot_side_time_generalization_group_report.sh`
- Expose matrix smoothing through `DAD_DOT_SIDE_TG_SMOOTH_SIGMA_SAMPLES`, falling back to `DAD_DOT_SIDE_SMOOTH_SIGMA_SAMPLES` when the time-generalization-specific variable is unset.

### `meg/cluster/README.md`
- Document qsub commands for regenerating the 1D and 2D group reports with smoothed subject values.

### `meg/docs/PIPELINE_SPEC_basic.md`
- Document the smoothing option, NaN-preserving policy, cluster-stat timing, and the default-grid interpretation of sigma in samples and milliseconds.

### `meg/tests/test_dot_side_group_report.py`
- Add a synthetic test that verifies 1D smoothing changes finite values without leaking across NaN gaps.

### `meg/tests/test_dot_side_time_generalization.py`
- Add a synthetic test that verifies 2D smoothing preserves NaN barriers and does not spread signal across disconnected matrix components.

### `meg/tests/test_cluster_wrappers.py`
- Extend wrapper tests to cover the new smoothing environment variables and `--smooth-sigma-samples` forwarding.

Suggested commit message:
`Add group smoothing for dot-side reports`

# 24 Jun 2026 19:02

## Add dot-side time-generalization decoding

### `meg/code/sanity_checks/03_dot_side_time_generalization.py`
- Add the participant-level train-time x test-time dot-side decoding workflow as a separate sanity-check script rather than expanding the existing 1D line decoder.
- Reuse the established Stage-04 validation, model-condition resolution, dynamic dot-side labels, posterior gradiometer selection, same-identity trial averaging, and regularized LDA behavior from the existing dot-side decoder.
- Build one stratified fold assignment per train time, reuse it across all test times, and store explicit NaNs for unsupported train rows, unbalanceable test cells, and train/test feature-dimension mismatches from edge-clipped windows.
- Export summary JSON, compressed NPZ matrices, long-form matrix CSV, participant heatmap PNG, and MNE HTML report with occlusion/reappearance markers on both axes.

### `meg/code/sanity_checks/04_dot_side_time_generalization_group_report.py`
- Add the group companion report for completed participant time-generalization outputs.
- Load and validate subject matrix bundles, enforce shared time axes and decoding settings, and compute group-valid cells only where every included subject has a finite value.
- Implement a 2D one-sample sign-flip cluster permutation over the train-time x test-time lattice with 4-neighbor adjacency, exact sign flips for small groups, configurable cluster thresholds, and NaN gaps that prevent cluster bridging.
- Write group mean matrix NPZ files plus group-first linked and embedded HTML reports, with significant 2D clusters shaded and outlined on the group heatmaps before participant-specific sections.

### `meg/cluster/submit_dot_side_time_generalization.sh`
- Add an SGE array wrapper that runs one participant per task and sequentially evaluates the four standard dot-side analysis keys.
- Expose decode grid, window width, same-identity trial averaging, random seed, fold count, validate-only, resume, rebuild, selected-analysis, and optional output-root controls through environment variables.

### `meg/cluster/submit_dot_side_time_generalization_group_report.sh`
- Add a one-task SGE aggregation wrapper for the time-generalization group report.
- Expose subject selection, metric selection, missing-subject handling, rebuild/resume behavior, and 2D cluster-permutation settings through documented environment variables.

### `meg/cluster/deploy_review_to_cluster.sh`
- Include the new time-generalization subject and group wrappers in the minimal cluster deployment payload.

### `meg/cluster/README.md`
- Document the cluster runbook for validating, submitting, monitoring, and aggregating dot-side time-generalization jobs.
- Record the output directory layout, matrix artifact names, default resource expectations, dependency-submission example, and 2D cluster-permutation interpretation.

### `meg/docs/PIPELINE_SPEC_basic.md`
- Register the time-generalization script as an optional sanity-check workflow.
- Document the matrix data flow, NPZ/CSV output contract, train/test NaN policies, edge-window feature policy, group all-subject finite-cell mask, 2D cluster statistics, and qsub wrappers.

### `meg/tests/test_dot_side_time_generalization.py`
- Add synthetic helper tests for train-row skipping, test-cell skipping, feature-dimension mismatch handling, 2D cluster detection, NaN gaps that split clusters, and group-first HTML generation.

### `meg/tests/test_cluster_wrappers.py`
- Extend wrapper and deployment tests so the new time-generalization qsub scripts parse, expose the expected controls, and are included in the cluster payload.

Suggested commit message:
`Add dot-side time-generalization decoding`

# 24 Jun 2026 18:33

## Add pilot abstract

### `paper/abstract.md`
- Add a concise pilot abstract describing the moving-object occlusion MEG study, dynamic RSA approach, observed/expected/prediction-error model comparison, and preliminary seven-participant findings.
- Summarize the current interpretation that position representations are weakened or transformed during occlusion and that post-reappearance activity follows observed sensory evidence more strongly than a separable prediction-error code.

Suggested commit message:
`Add pilot abstract`

# 24 Jun 2026 18:32

## Document the multisubject MEG pipeline state

### `meg/docs/MULTISUBJECT_MIGRATION_AUDIT.md`
- Add the full migration audit for the external seven-subject MEG dataset.
- Record verified subject profiles, read-only source counts, trigger extraction evidence, Stage-03 preprocessing policy, epoch/export contracts, MATLAB dRSA orchestration, cluster execution, dot-side sanity checks, and remaining caveats.

### `meg/docs/MULTISUBJECT_REFACTOR_PROMPTS.md`
- Update the implementation prompts to reflect completed chats, current pipeline contracts, manual-ICA production status, downstream lineage gates, cluster execution, dot-side sanity workflows, and Stage-08 line-summary work.
- Clarify future-work boundaries so optional sanity checks and reports do not become required scientific dependencies.

### `meg/docs/PIPELINE_SPEC_basic.md`
- Reconcile the pipeline specification with the new configuration/profile model, source manifest, Stage-03 review/apply flow, epoching/export gates, MATLAB handoff, cluster wrappers, and optional sanity workflows.
- Update data-flow, output, validation, and execution details for the current multisubject pipeline.

### `meg/docs/TOOD_AND_DISCUSS.md`
- Add follow-up discussion notes about report interpretability and position decoding.
- Update the line-noise discussion with MNE spectrum-fit filtering as an option.

### `meg/docs/TRIGGERS_EVENTS.md`
- Update the trigger/event contract for profile-specific condition codes, canonical event fields, transition-warning preservation, deep report outputs, observed-only context inference, and downstream epoching mappings.

### `meg/docs/deviations/DEVIATION_03B_SUB01_SKIP_EOG_BLINK.md`
- Replace the older sub01 blink-skip description with the current diagnostics-only EOG policy and explicit manual ICA component-review requirement.

### `meg/docs/deviations/DEVIATION_03C_ICA_DEV_HEAD_T.md`
- Mark the old native-position and ICA000 blink-workaround decisions as superseded.
- Document the current central Maxwell destination policy and explicit manual ICA exclusion model.

### `meg/docs/meg_preprocessing_pipeline.mmd`
- Update the Mermaid source to show the resolved configuration flow, source manifest, Stage-03 review/apply decisions, manual-ICA lineage gate, epoching, dRSA export, optional ERF/decoding checks, and deprecated auto-ICA branch.

Suggested commit message:
`Document multisubject MEG pipeline state`

# 24 Jun 2026 18:31

## Add CIMeC cluster orchestration

### `CIMeC_cluster_docs.md`
- Add the user-facing CIMeC cluster reference covering VPN/SSH access, compute-node usage, qsub basics, environment activation, and MATLAB submission examples.

### `meg/cluster/README.md`
- Add the DAD-specific CIMeC cluster runbook for deployment, micromamba environment creation, preprocessing review/apply, epoch/export, MATLAB dRSA, line summaries, and monitoring.

### `meg/cluster/deploy_review_to_cluster.sh`
- Add a guarded rsync deployment helper that copies the minimal DAD MEG processing payload to the cluster destination outside the raw data tree.
- Include preprocessing, downstream, MATLAB analysis, sanity-check, parameter, cluster, calibration, and simulation helper files.

### `meg/cluster/environment.yml`
- Add the conda-forge `dad-meg` environment definition with pinned Python, MNE, numerical, plotting, and machine-learning packages for cluster execution.

### `meg/cluster/requirements.txt`
- Add a binary-wheel-only pip fallback for the same headless runtime dependencies.

### `meg/cluster/submit_dot_side_decoding.sh`
- Add an SGE array wrapper for optional per-subject dot-side decoding sanity checks with documented decoding, averaging, seed, validation, and rebuild controls.

### `meg/cluster/submit_dot_side_group_report.sh`
- Add a one-task SGE wrapper for combined group dot-side reports, cluster-permutation settings, missing-subject handling, and validation-only runs.

### `meg/cluster/submit_drsa_analysis.sh`
- Add an SGE array wrapper for optional Stage-07 MATLAB dRSA with MATLAB binary discovery, model-source controls, analysis selection, and memory notes.

### `meg/cluster/submit_drsa_line_summaries.sh`
- Add a one-task SGE wrapper for Stage-08 dRSA line-summary reports with subject, condition, sensor-group, smoothing, shading, and interval-lag controls.

### `meg/cluster/submit_epoch_export.sh`
- Add an SGE array wrapper for Stage-04 epoching followed by Stage-06 neural export, with resume/rebuild behavior and source-manifest prerequisites.

### `meg/cluster/submit_preprocessing_apply.sh`
- Add an apply-pass wrapper that delegates to the shared review wrapper with explicit apply and rebuild settings.

### `meg/cluster/submit_preprocessing_review.sh`
- Add the main SGE array wrapper for Stage-03 review, including selected-stage recovery, resume/rebuild behavior, tuning reuse, and headless runtime setup.

### `meg/code/preprocessing_basic/run_pipeline.py`
- Register stable preprocessing, epoching, export, and optional MATLAB dRSA stages.
- Forward site/tuning configs, source manifests, review/apply controls, overwrite/resume flags, and deep-trigger-report options to stage subprocesses.
- Refresh source manifests while preserving trigger-derived associations and derivative lineage.

### `meg/params/pipeline_cluster.yaml`
- Add CIMeC path overrides for project, raw data, derivatives, docs, params, Maxwell resources, and MATLAB binary selection without changing scientific defaults.

### `meg/tests/test_cluster_wrappers.py`
- Cover wrapper syntax, deployment payload contents, dot-side controls, group-report controls, Stage-08 controls, MATLAB model-source controls, and CIMeC MATLAB resource discovery.

### `meg/tests/test_run_pipeline.py`
- Cover source-manifest forwarding, site-config forwarding, deep-report opt-out, review controls, downstream noninteractive controls, resume forwarding, and manifest preservation.

Suggested commit message:
`Add CIMeC cluster orchestration`

# 24 Jun 2026 18:30

## Add dot-side decoding sanity reports

### `meg/code/sanity_checks/01_dot_side_decoding.py`
- Refactor the dot-side sanity check to resolve configuration, subject paths, model MATs, and raw-tree output guards through the shared MEG pipeline infrastructure.
- Add observed/predicted model-source selection, with predicted occluded-deviant labels derived from the matched `occluded_nondeviant` trajectory in the same observed MAT.
- Add centered sliding-window decoding, deterministic same-identity trial averaging, all-trials-per-identity averaging, occlusion timing markers, richer summary JSON, and validate-only/resume/overwrite controls.

### `meg/code/sanity_checks/02_dot_side_group_report.py`
- Add a combined group report builder for completed per-subject dot-side outputs.
- Implement standard analysis registry, group figures, source-trial count plots, linked and embedded HTML reports, and one-sample sign-flip cluster permutation testing.
- Support missing-subject handling, CI/SEM display bands, configurable metrics, and occlusion timing markers.

### `meg/tests/test_dot_side_decoding.py`
- Cover sliding-window edge clipping, odd/even window centering, same-identity averaging, all-trials-per-identity mode, 20/40-identity model validation, predicted deviant model selection, output paths, skipped windows, and diagnostics.

### `meg/tests/test_dot_side_group_report.py`
- Cover cluster permutation behavior, CI/SEM intervals, NaN gaps, standard analysis keys, cross-analysis setting validation, report ordering, missing-subject handling, time-axis mismatch rejection, combined report generation, and occlusion markers.

Suggested commit message:
`Add dot-side decoding sanity reports`

# 24 Jun 2026 18:29

## Add MATLAB dRSA handoff and line summaries

### `meg/code/analysis/A_dRSA.m`
- Refactor real-neural dRSA to consume the Python-generated Stage-07 handoff instead of hardcoded subject paths and V28/predicted MAT filenames.
- Validate shared identity, time, and channel axes from Stage-06 sidecars; remove all-NaN channels before analysis.
- Build observed, expected/predicted, and prediction-error models from the shared observed trajectory MAT.

### `meg/code/analysis/A_dRSA_alt.m`
- Refactor alternate sensor-group dRSA to use the same handoff, model-selection, dynamic identity, and output conventions as the main real-neural script.
- Keep sensor-group analyses selectable without reconstructing paths inside MATLAB.

### `meg/code/analysis/dRSA_simulations.m`
- Refactor model-vs-model and neural-vs-neural simulations to use handoff-driven paths, model sources, sensor groups, and smoothed display outputs.
- Preserve expected/predicted occluded-deviant handling through the matched nondeviant trajectory condition.

### `meg/code/analysis/drsa_stage_utils.m`
- Add shared MATLAB utilities for loading handoffs, validating Stage-06 neural sidecars, selecting trajectory trials, enforcing expected-path correspondence, building position/direction tensors, constructing dRSA parameters, exporting figures, smoothing matrices, and normalizing sensor groups.

### `meg/code/preprocessing_basic/07_run_matlab_drsa.py`
- Add the optional Stage-07 Python orchestration layer that validates Stage-06 exports, writes a MATLAB handoff JSON, launches selected MATLAB analysis families, records reports, supports resume/overwrite, and updates MATLAB dRSA lineage.
- Include local/cluster MATLAB binary validation, optional `arch -x86_64` startup, model-source selection, smoothed plot options, and group-average report generation.

### `meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py`
- Add Stage-08 summary plotting from existing Stage-07 result MAT files without rerunning MATLAB.
- Compute interval lag curves, group averages, optional display smoothing, and optional shaded uncertainty bands with linked and embedded HTML reports.

### `meg/code/preprocessing_basic/drsa_result_utils.py`
- Add shared Python helpers for reading MATLAB v7/v7.3 dRSA results, normalizing result metadata, deriving occlusion timing, smoothing display matrices, and embedding report images.

### `meg/tests/test_drsa_line_summaries.py`
- Cover lag-curve conventions, occlusion-frame conversion, interval support modes, shaded bands, NaN handling, result-label inference, CLI execution on synthetic derivatives, and embedded report generation.

### `meg/tests/test_matlab_drsa_stage.py`
- Cover Stage-07 handoff validation, expected output path coverage, MATLAB batch command construction, observed/predicted model-source selection, and group-average report behavior.

Suggested commit message:
`Add MATLAB dRSA handoff and line summaries`

# 24 Jun 2026 18:28

## Add lineage-gated epoching and dRSA export

### `meg/code/preprocessing_basic/04_epoching.py`
- Refactor epoching into a noninteractive CLI that requires completed manual-ICA apply lineage from `source_manifest.json`.
- Build run-isolated trial decisions from canonical step-02 events, preserving source FIF, acquisition run, raw sample frame, behavioral context, replay state, identity, condition, and exclusion evidence.
- Write concatenated epochs, surviving metadata, full trial decisions, reports, and epoching lineage while documenting the -0.5 to 2.667 s window, no-baseline policy, detrending, and annotation rejection behavior.

### `meg/code/preprocessing_basic/06_export_neural_for_drsa.py`
- Refactor dRSA neural export into a configuration-driven CLI that consumes completed Stage-04 epochs and lineage.
- Export one identity-averaged MAT tensor and JSON sidecar per condition on a fixed canonical MEG channel axis, with bad or absent channels represented as NaN.
- Validate identity coverage, condition mappings, trajectory-model discovery, time axis, epoch report status, and source-manifest lineage before writing outputs.

### `meg/code/preprocessing_basic/XX_ERF_manual.py`
- Update the optional ERF inspection workflow to resolve the configured subject, epochs path, and condition definitions through `pipeline_config.py`.
- Keep the script interactive and report-free while documenting the inspection-only -0.1 to 0 s baseline.

### `meg/code/preprocessing_basic/downstream_stage_utils.py`
- Add shared downstream contracts for manual-ICA apply gating, trajectory-model discovery, canonical trigger row loading, run-grouped event rows, trial-decision construction, exclusion policies, and identity coverage.
- Enforce source-stem lineage matching instead of sorted-position or filename-substitution assumptions.

### `meg/code/preprocessing_basic/drsa_canonical_meg_channels.txt`
- Clarify that the file is a repository-owned fixed Neuromag channel axis and is not auto-generated by export.

### `meg/tests/test_downstream_stage_utils.py`
- Cover observed-only trajectory discovery, v18/v20 trigger generations, gaze/catch/ESC exclusions, identity-code collisions, replay handling, run-boundary isolation, partial ESC trials, identity coverage, and lineage matching.

### `meg/tests/test_epoch_export_smoke.py`
- Add a synthetic smoke test that exercises epoching and dRSA export contracts without full participant preprocessing.

Suggested commit message:
`Add lineage-gated epoching and dRSA export`

# 24 Jun 2026 18:27

## Add lineage-aware stage 03 preprocessing

### `meg/code/preprocessing_basic/03_A_MaxFilt.py`
- Refactor Maxwell filtering into a configuration-driven CLI with review/apply passes, explicit overwrite/resume behavior, and source-manifest lineage updates.
- Discover raw acquisition runs from the resolved profile and preserve exact source stems in `_raw_sss.fif` outputs.
- Add MNE digitization compatibility checks, per-run bad-channel detection, configurable bad-channel policies, central-run destination selection, before/after QC, and guarded report outputs.

### `meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py`
- Refactor artifact annotation into a configuration-driven CLI that consumes lineage-matched SSS outputs and writes `_ann.fif`/`_ann.csv` derivatives.
- Preserve MISC channels, separate EOG diagnostics from annotation policy, and default missing-EOG subjects to disabled EOG-only processing.
- Add review/apply reporting, threshold tuning support, resume validation, and source-manifest lineage updates.

### `meg/code/preprocessing_basic/03_C_ICA.py`
- Refactor manual ICA into the production stage-03 ICA branch with explicit component exclusions, compatible prefit cache reuse, review/apply passes, and guarded cleaned FIF outputs.
- Fit MEG-only data with configured resampling/filtering while keeping EOG scores diagnostic and preserving original non-MEG channels in cleaned run outputs.
- Validate resume compatibility through lineage, exclusions, channel identity, bad channels, sample counts, and MISC-channel preservation.

### `meg/code/preprocessing_basic/03_C_ICA_auto.py`
- Mark the automatic MEGNet branch as deprecated experimental code.
- Add a production gate so normal execution exits before optional MEGNet dependencies are imported unless a developer-only escape hatch is used.

### `meg/code/preprocessing_basic/preprocessing_stage_utils.py`
- Add shared acquisition-run discovery, run-label parsing, derivative-path construction, guarded output helpers, JSON/text/HTML writers, tuning-file merging, and source-manifest stage updates.
- Centralize lineage records so 03_A, 03_B, and 03_C use the same source/stage provenance contract.

### `meg/code/preprocessing_basic/render_stage03_mne_reports.py`
- Add a report-only utility that renders self-contained MNE HTML reports from existing stage-03 JSON metadata, PNG figures, and optional tuning YAML.
- Validate report inputs and enforce raw-tree output guards without opening raw FIF files or recomputing scientific derivatives.

### `meg/requirements-ica-auto-experimental.txt`
- Move MEGNet/ONNX automatic-ICA dependencies into an explicit optional file for the deprecated research branch.

### `meg/requirements.txt`
- Add runtime dependencies needed by the production refactor and remove automatic-ICA packages from the production environment.

### `meg/tests/test_preprocessing_stage_policies.py`
- Cover Maxwell bad-channel policies, digitization errors, central destination selection, resume validation, EOG mode resolution, manual ICA parsing, mount-independent ICA cache fingerprints, and auto-ICA deprecation.

### `meg/tests/test_preprocessing_stage_smoke.py`
- Add synthetic smoke tests for artifact annotation, manual ICA apply behavior, headless ICA diagnostics, and Maxwell before/after QC using existing FIF outputs.

### `meg/tests/test_preprocessing_stage_utils.py`
- Cover run-label parsing, exact source-stem derivative names, tuning-file merges, raw-tree rejection, NumPy JSON normalization, and manifest-stage preservation.

### `meg/tests/test_real_preprocessing_headers.py`
- Add mounted-dataset header checks for expected ocular modalities across real subject task FIFs, skipped automatically when the raw dataset is unavailable.

### `meg/tests/test_render_stage03_mne_reports.py`
- Cover MNE report rendering from synthetic stage JSON/PNG inputs and raw-tree output refusal.

Suggested commit message:
`Add lineage-aware stage 03 preprocessing`

# 24 Jun 2026 18:26

## Add multisubject config and trigger extraction foundation

### `meg/code/preprocessing_basic/02_trigger_extraction_custom.py`
- Convert the numbered trigger stage into a documented wrapper around reusable trigger extraction logic.
- Keep observed STI101 pulses as the scientific event source and route reported CSV rows through the new alignment/QC implementation.

### `meg/code/preprocessing_basic/dataset_inventory.py`
- Expand read-only discovery from a subject-01 preflight into explicit multisubject source discovery for task FIFs, empty-room recordings, trajectory MATs, TrialStruct MATs, behavioral MATs, debug CSVs, and EDF files.
- Add durable source-manifest construction that preserves source filenames and keeps MEG acquisition runs, behavioral blocks, trigger CSVs, and EDF files separate until evidence-based association.
- Update inventory reports and validation to reflect the external raw tree and profile-specific EDF expectations.

### `meg/code/preprocessing_basic/deep_trigger_report.py`
- Add a legacy-compatible deep trigger report renderer driven by canonical step-02 events, alignments, transition warnings, and segment summaries.
- Preserve familiar alignment/anomaly CSV and plot outputs while avoiding a second independent FIF/CSV alignment pass.

### `meg/code/preprocessing_basic/pipeline_config.py`
- Add supported subject normalization for participants 01-07 and a verified subject-profile selection layer.
- Merge default, site, subject, and tuning YAMLs in a defined order, with selected-profile facts deriving trigger ranges, condition maps, EOG settings, and discovery globs.
- Resolve the unique external `subNN_*` raw folder, derivative paths, subject input/output folders, and raw-root guards from the selected profile.

### `meg/code/preprocessing_basic/trigger_extraction.py`
- Add the reusable implementation for non-zero STI101 window reconstruction, transition-pulse resolution, saved debug CSV parsing, segment alignment, replay/gaze/ESC metadata, source associations, and canonical output writing.
- Preserve raw sample frames and source FIF identity for downstream run-isolated epoching.
- Treat reported rows as alignment evidence only; missing reported rows do not create MEG events.

### `meg/params/pipeline_default.yaml`
- Replace subject-01-specific defaults with shared pipeline settings plus verified subject profiles for participants 01-07.
- Store per-subject versions, condition codes, identity ranges, trigger policies, ocular modalities, completed blocks, and known exceptions once in `subject_profiles`.
- Add derived stage defaults for trigger parsing, artifact handling, ICA, epoching, dRSA export, MATLAB dRSA, and eye-model skip policy.

### `meg/params/sub-01.yaml`
- Reduce the subject file to a selector for participant 01, with scientific facts centralized in `pipeline_default.yaml`.

### `meg/params/sub-02.yaml`
- Add a selector file for participant 02.

### `meg/params/sub-03.yaml`
- Add a selector file for participant 03.

### `meg/params/sub-04.yaml`
- Add a selector file for participant 04.

### `meg/params/sub-05.yaml`
- Add a selector file for participant 05.

### `meg/params/sub-06.yaml`
- Add a selector file for participant 06.

### `meg/params/sub-07.yaml`
- Add a selector file for participant 07.

### `meg/tests/test_deep_trigger_report.py`
- Cover canonical-to-legacy deep report generation and raw-tree output rejection.

### `meg/tests/test_multisubject_discovery.py`
- Cover exact subject-folder discovery, profile-derived stage configuration, sub01 irregular run preservation, sub07 date handling, empty-room separation, inventory validation, and source-manifest modality separation.

### `meg/tests/test_no_edf_behavior.py`
- Update the no-EDF test to the external subject-folder layout and profile-driven EDF expectation.

### `meg/tests/test_pipeline_config.py`
- Extend configuration tests for profile completeness, site/tuning merge precedence, selector normalization, tuning overrides, derivative rooting, and raw-tree guards.

### `meg/tests/test_trigger_parser.py`
- Add coverage for v18/v20 CSV parsing, v20 identity-code collisions, replay/gaze/ESC metadata, multi-block acquisition FIF alignment, exact-match-only source associations, and observed-only ESC context inheritance.

Suggested commit message:
`Add multisubject config and trigger foundation`

# 24 Jun 2026 18:25

## Use nondeviant paths as expected continuation

### `experiment/lib/Config_stimuli_generation_V30_rescueTraject.m`
- Remove the separate `predictedFilePattern` contract from the main V30 stimulus configuration.
- Keep the observed MAT as the single canonical generated stimulus file for downstream experiment and analysis code.

### `experiment/practice/lib/Config_stimuli_generation_V28_rescueTraject.m`
- Remove the practice V28 predicted-file pattern so practice generation follows the same single-file expected-continuation contract.

### `experiment/practice/run_practice_occlusion_v1.m`
- Stop checking, copying, caching, and returning FPS-specific predicted MAT files.
- Treat the canonical observed MAT and schedule MAT as the complete practice input set.

### `experiment/practice/stimuli_generation_V28_rescueTraject.m`
- Remove the predicted baseline branch, `xySeqsPredicted` construction, and separate predicted MAT save.
- Document that the sequence-matched `occluded_nondeviant` trial is the expected continuation for occluded-deviant analyses.
- Simplify trajectory generation and fixation-zone translation helpers from triplets to nondeviant/deviant pairs.

### `experiment/runExpReadable.m`
- Update launcher documentation and skip logic so V30 generation is considered complete when the observed MAT exists.

### `experiment/run_MEG_lab/lib/Config_stimuli_generation_V30_rescueTraject.m`
- Remove the lab-copy predicted-file pattern to keep deployment parity with the main V30 generator.

### `experiment/run_MEG_lab/runExpReadable.m`
- Update MEG-lab launch checks to require only the observed V30 MAT before skipping stimulus generation.

### `experiment/run_MEG_lab/stimuli_generation_V30_rescueTraject.m`
- Remove predicted trajectory generation and saving from the MEG-lab V30 copy.
- Use paired nondeviant/deviant base trajectories and the observed condition set as the downstream model source.

### `experiment/run_outside_lab/RUNME.m`
- Update unused-subject detection so outside-lab subject slots are occupied by observed stimulus or TrialStruct files, not by the removed predicted MAT.

### `experiment/run_outside_lab/lib/Config_stimuli_generation_V29_rescueTraject.m`
- Remove the outside-lab V29 predicted-file pattern.

### `experiment/run_outside_lab/lib/Config_stimuli_generation_V30_rescueTraject.m`
- Remove the outside-lab V30 predicted-file pattern.

### `experiment/run_outside_lab/stimuli_generation_V29_rescueTraject.m`
- Remove V29 predicted branch generation and save logic.
- Retain one observed MAT containing the condition-matched nondeviant continuation.

### `experiment/run_outside_lab/stimuli_generation_V30_rescueTraject.m`
- Remove V30 predicted branch generation and save logic from the outside-lab copy.
- Keep the condition set in one MAT as the complete model source for downstream expected-versus-observed comparisons.

### `experiment/stimuli_generation_V30_rescueTraject.m`
- Remove the main V30 predicted MAT output and predicted trajectory helper flow.
- Preserve sequence-matched shared paths in `always_visible` and `occluded_nondeviant`, with occluded-deviant sharing the prefix through the deviance frame.

### `experiment/stimuli_generation_versions.md`
- Document the single-file expected-continuation contract for V28, V29, and V30 rescueTraject outputs.
- State that standalone expected-path MAT files are intentionally omitted.

### `simulations/functions/build_movdot_simulation_inputs.m`
- Drop `xySeqsPredicted` autodetection and predicted-specific condition handling.
- Require `xySeqs` as the model source so simulation inputs follow the new observed-MAT contract.

### `simulations/scripts/one-dot/PE_simulation_diff_1Dot_occlusion_v2.m`
- Resolve only one compatible observed stimulus MAT.
- Load sequence-matched `occluded_nondeviant` trials as the expected/predicted continuation and validate sequence identity alignment.
- Record expected-source metadata in reproducibility output instead of predicted-file provenance.

### `simulations/scripts/one-dot/plot_occlusion_paths_1Dot_v2.m`
- Update resolver logic to select compatible MAT files by `xySeqs` presence rather than `_predicted` naming.

Suggested commit message:
`Use nondeviant paths as expected continuation`

# 24 Jun 2026 18:24

## Ignore PowerPoint artifacts

### `.gitignore`
- Add `*.pptx` to the generated/binary artifact ignore block so PowerPoint decks do not appear as normal untracked commit candidates.
- This intentionally keeps `results_pilot.pptx` out of the working-tree commit plan unless it is force-added later as a deliberate binary artifact.

Suggested commit message:
`Ignore PowerPoint artifacts`

# 10 Jun 2026 15:04

## Update current methods and preserve pilot protocols

### `experiment/timing.md` (deleted)
- Remove the obsolete static timing assessment for the V14/V24 path-band stack.
- Its completed hardening items, including strict MEG sync policy, full flip diagnostics, missed-frame guardrails, and precomputed occluder geometry, are implemented and documented in the current V20 runtime and methods.

### `paper/methods.md`
- Update the primary methods to the V30 generation, v23 scheduling, V20 runtime, and v10 trigger map used for subjects 04-07.
- Add participant-specific presented-stack, completed-trial, replay, catch, trigger-range, ITI, fixation-window, color-cue, and eye-recording differences for subjects 01-07.
- Describe 40-path generation, zero-inclusive deviant-turn sampling with effective-deviance guards, 15% catch insertion, centered prompts, 1000-1500 ms ITIs, fixation cues, blink breaks, message cushions, and current timing estimates.
- Clarify that reported counts come from saved/presented data rather than prepared schedules and that subject 01 retained online replay behavior despite loss of offline eye-tracker files.

### `paper/methods_pilot_sub01_02_03/methods_subjects_01_02_03.md`
- Preserve the earlier V29/v22/v19 protocol as a dedicated reference for analyses of subjects 01-03.
- Document the 20-path schedule, legacy trigger codes, catch counts, visual timing, message flow, and run-color behavior separately from the current V30/v23/v20 methods.

Suggested commit message:
`Update methods and archive pilot protocols`

# 10 Jun 2026 15:03

## Add the MEG preprocessing, trigger-QA, and dRSA pipeline

### `meg/code/analysis/A_dRSA.m`
- Add real-neural dRSA for observed/predicted one-dot models, prediction-error models, channel consistency checks, NaN-channel removal, correlation/PCR execution, and MAT/figure output.

### `meg/code/analysis/A_dRSA_alt.m`
- Add a configurable single-condition dRSA variant with observed/predicted model selection and whole-brain or lobe-specific MEGIN sensor groups.

### `meg/code/analysis/dRSA_simulations.m`
- Add model-vs-model and neural-vs-neural simulation analyses across conditions and sensor groups with redirectable output roots.

### `meg/code/preprocessing_basic/02_trigger_extraction_custom.py`
- Add config-driven FIF trigger reconstruction, expected-trigger-first transition resolution, machine-readable event outputs, QA summaries, and integrated deep reported-vs-observed confrontation reporting.
- Preserve both global concatenated timing and source-FIF timing fields for downstream epoching and auditability.

### `meg/code/preprocessing_basic/03_A_MaxFilt.py`
- Add first-look QC, Maxwell bad-channel detection, calibration/cross-talk checks, SSS/Maxwell filtering, PSD figures, HTML reporting, and bad-channel propagation.

### `meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py`
- Add blink/muscle annotation, subject-01 EOG exception handling, annotation CSV/FIF output, QC figures, reports, and bad-channel preservation.

### `meg/code/preprocessing_basic/03_C_ICA.py`
- Add the manual ICA branch with FastICA fitting, subject-specific component exclusion, blink-window handling, cleaned FIF/model output, figures, and HTML reporting.

### `meg/code/preprocessing_basic/03_C_ICA_auto.py`
- Add an automatic MEGNet ICA-label branch with probability-threshold exclusions for eye and heart components, saved models/labels, cleaned FIF output, and reports.

### `meg/code/preprocessing_basic/04_epoching.py`
- Add trigger-table-to-annotation reconstruction, cross-run concatenation, condition epoching, sequence-identity metadata assignment, integrity validation, saved epochs, and report figures.

### `meg/code/preprocessing_basic/06_export_neural_for_drsa.py`
- Add condition/identity averaging, canonical-channel alignment, NaN placeholders for bad or absent channels, 120 Hz export, MAT output, and dRSA-preparation QC reporting.

### `meg/code/preprocessing_basic/XX_ERF_manual.py`
- Add an interactive IPython ERF inspection workflow for condition averages, trial images, topo traces, and topomaps without writing derivatives.

### `meg/code/preprocessing_basic/dataset_inventory.py`
- Add read-only dataset discovery and FIF-header preflight with optional JSON/Markdown inventories and explicit no-EDF skip reporting.

### `meg/code/preprocessing_basic/drsa_canonical_meg_channels.txt`
- Add the stable ordered MEG channel axis used to keep dRSA tensor dimensions consistent across subjects and preprocessing outcomes.

### `meg/code/preprocessing_basic/pipeline_config.py`
- Add recursive YAML loading/merging, required-section validation, path resolution, provenance snapshots, and guards that reject outputs inside the raw-data tree.

### `meg/code/preprocessing_basic/run_pipeline.py`
- Add the orchestration entry point for configuration resolution, automatic preflight, run logging, dry runs, and the currently registered trigger-extraction stage.

### `meg/code/sanity_checks/01_dot_side_decoding.py`
- Add time-resolved left/right dot-side decoding with balanced and unbalanced classifiers, contra/ipsi posterior ERFs, summary data, and MNE report output.

### `meg/docs/MULTISUBJECT_REFACTOR_PROMPTS.md`
- Record staged implementation prompts and shared constraints for extending the currently subject-01-oriented pipeline to heterogeneous subjects 01-07.

### `meg/docs/PIPELINE_SPEC_basic.md`
- Document implemented infrastructure and scientific stages, FLUX anchors, data flow, outputs, and the editable Mermaid pipeline diagram.

### `meg/docs/RUN_LOG.md`
- Record configuration/inventory refactoring, derivative-root migration, test baselines, and prior trigger-extraction execution evidence.

### `meg/docs/TOOD_AND_DISCUSS.md`
- Capture unresolved preprocessing, cleaning, head-shift, decoding, ERF, and eye-tracking discussion items.

### `meg/docs/TRIGGERS_EVENTS.md`
- Define the trigger/event contract, expected fixed and dynamic codes, transition-pulse policy, event fields, QA requirements, and downstream epoching interface.

### `meg/docs/deviations/DEVIATION_03B_SUB01_SKIP_EOG_BLINK.md`
- Document the subject-01 deviation that skips unreliable EOG blink annotation in step 03_B and defers blink handling to ICA.

### `meg/docs/deviations/DEVIATION_03C_ICA_DEV_HEAD_T.md`
- Document the ICA fitting decision around device-to-head transforms, cross-run application, risks, alternatives, and validation evidence.

### `meg/docs/meg_preprocessing_pipeline.mmd`
- Add the editable Mermaid source for the preprocessing, trigger, ICA, epoching, and dRSA export workflow.

### `meg/docs/render_meg_preprocessing_pipeline.sh`
- Add a location-independent Mermaid CLI renderer for ignored SVG and PDF pipeline exports.

### `meg/params/pipeline_default.yaml`
- Add global project/layout paths, FIF/trigger/eye discovery rules, trigger-parser settings, required-code coverage, stage defaults, and runtime policies.

### `meg/params/sub-01.yaml`
- Add subject-01 path, naming, discovery, and overwrite/save-intermediate overrides.

### `meg/requirements-dev.txt`
- Add pytest and coverage dependencies on top of the runtime requirements.

### `meg/requirements.txt`
- Define the YAML, numerical, MNE, plotting, Qt, pandas, scikit-learn, ICA-label, and ONNX runtime dependency ranges.

### `meg/tests/test_no_edf_behavior.py`
- Test that missing eye files produce an explicit, valid eye-model skip rather than a pipeline failure when configured to skip.

### `meg/tests/test_pipeline_config.py`
- Test derivative-template rooting, raw-tree write protection, and configuration path behavior.

### `meg/tests/test_trigger_parser.py`
- Test symmetric expected-trigger-first collapse, width fallback, and unresolved deterministic fallback for contiguous trigger transitions.

### `triggers_check/compare_reported_vs_observed_triggers_subject.py`
- Extend FIF discovery to `_meg.fif`, legacy `.fif`, and gzipped FIF names.
- Resolve transition pairs symmetrically against the expected trigger stream, then use width and deterministic fallbacks while reporting decision source, confidence, order, expected code, and unresolved state.
- Add configurable expected-stream lookahead and expose the richer transition metadata in plots and Markdown reports.

### Validation
- `python3 -m pytest meg/tests -q`: 6 passed.
- `python3 -m compileall -q meg/code triggers_check/compare_reported_vs_observed_triggers_subject.py`: passed.
- `bash -n meg/docs/render_meg_preprocessing_pipeline.sh`: passed.
- MATLAB `checkcode` parsed all three dRSA scripts; two had zero messages and `dRSA_simulations.m` had one advisory message.

Suggested commit message:
`Add MEG preprocessing and dRSA pipeline`

# 10 Jun 2026 15:01

## Add a standalone laptop-safe experiment bundle

### `experiment/run_outside_lab/CreateInputFiles_v22_rescueTraject.m`
- Include the compatibility v22 scheduler for V29/v19 datasets and previously generated outside-lab subjects.

### `experiment/run_outside_lab/CreateInputFiles_v23_rescueTraject.m`
- Include the current v23 three-run schedule and 15% catch-plan builder for V30/v20 runs.

### `experiment/run_outside_lab/MoveDot1_experiment_occlusion_v19_rescueTraject.m`
- Preserve the prior v19 laptop runtime for compatibility with V29/v22 inputs and historical outside-lab sessions.

### `experiment/run_outside_lab/MoveDot1_experiment_occlusion_v20_rescueTraject.m`
- Add the current laptop-adapted V20 runtime with keyboard arrow mappings, centered catch labels, non-MEG operation, and no dependency on EyeLink hardware.

### `experiment/run_outside_lab/RUNME.m`
- Add the primary one-command laptop launcher with startup dialog, automatic unused subject selection, display refresh-rate detection, local directory creation, and non-MEG runtime defaults.
- Scale V30 timing metadata to the detected display refresh rate and keep generation, scheduling, and runtime self-contained in the copied folder.

### `experiment/run_outside_lab/fixation.bmp`
- Include the confirmed fixation bitmap used by the copied runtime.
- Binary warning: this required bitmap cannot be line-reviewed. The user explicitly confirmed the asset is intentional; to unstage it before commit: `git reset HEAD -- experiment/run_outside_lab/fixation.bmp`.

### `experiment/run_outside_lab/lib/Config_runtime_v19_rescueTraject.m`
- Include the compatibility v19 runtime preset for prior outside-lab datasets.

### `experiment/run_outside_lab/lib/Config_runtime_v20_rescueTraject.m`
- Include the current V20 preset surface for ITI cues, catch layout, blink breaks, path-band geometry, and optional fake gaze.

### `experiment/run_outside_lab/lib/Config_schedule_CreateInputV22_MoveDotV19_rescueTraject.m`
- Include schedule/catch/message constants and V29/v22/v19 filename contracts.

### `experiment/run_outside_lab/lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m`
- Adapt v23 message gates and calibration-choice prompts to LEFT/RIGHT keyboard controls while retaining the current catch and block design.

### `experiment/run_outside_lab/lib/Config_stimuli_generation_V29_rescueTraject.m`
- Include V29 generation parameters for compatibility with existing outside-lab subject files.

### `experiment/run_outside_lab/lib/Config_stimuli_generation_V30_rescueTraject.m`
- Include the current 40-path V30 generation parameters.

### `experiment/run_outside_lab/lib/Utils.m`
- Package the local utility class required by both generation lineages.

### `experiment/run_outside_lab/lib/createNewPoint.m`
- Package the local boundary and bounce helper required by trajectory synthesis.

### `experiment/run_outside_lab/lib/plot_paths.m`
- Include trajectory preview rendering without requiring parent-repository helpers.

### `experiment/run_outside_lab/stimuli_generation_V29_rescueTraject.m`
- Preserve V29 config-driven generation for compatibility runs.

### `experiment/run_outside_lab/stimuli_generation_V30_rescueTraject.m`
- Extend V30 generation with optional target-display refresh-rate scaling for frame count, per-frame speed, and fixed event indices.
- Resolve plotting helpers locally when the full repository simulation tree is unavailable.

Suggested commit message:
`Add standalone laptop experiment bundle`

# 10 Jun 2026 15:00

## Add a standalone V30/v23/v20 MEG-lab experiment bundle

### `experiment/run_MEG_lab/CreateInputFiles_v23_rescueTraject.m`
- Package the canonical v23 schedule/catch-plan builder inside the self-contained lab directory.

### `experiment/run_MEG_lab/MEG_lab_procedure.md`
- Include the acquisition-room setup, participant preparation, EyeLink, head-position, response mapping, shutdown, and troubleshooting procedure with lab-copy-specific wording.

### `experiment/run_MEG_lab/MoveDot1_experiment_occlusion_v20_rescueTraject.m`
- Package the complete V20 runtime with MEG/DataPixx triggers, EyeLink replay handling, timing diagnostics, ITI cues, blink breaks, catch scoring, and block-level saves.

### `experiment/run_MEG_lab/README_RUN_MEG_LAB.md`
- Document interactive and Intel-architecture MATLAB launch commands, active V30/v23/v20 defaults, block-resume editing, and skip-if-existing behavior.

### `experiment/run_MEG_lab/fixation.bmp`
- Include the confirmed runtime fixation bitmap required for EyeLink image transfer in the standalone copy.
- Binary warning: this required bitmap cannot be line-reviewed. The user explicitly confirmed the asset is intentional; to unstage it before commit: `git reset HEAD -- experiment/run_MEG_lab/fixation.bmp`.

### `experiment/run_MEG_lab/lib/Config_runtime_v20_rescueTraject.m`
- Specialize the runtime preset for real lab acquisition: EyeLink enabled, fake gaze disabled, 3-degree fixation threshold, bounded replay, and calibration enabled.

### `experiment/run_MEG_lab/lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m`
- Include the matching v23 schedule, catch, message, path-band, timing, and EyeLink policy constants.

### `experiment/run_MEG_lab/lib/Config_stimuli_generation_V30_rescueTraject.m`
- Include the matching 40-path V30 stimulus-generation constants and fixed occlusion geometry.

### `experiment/run_MEG_lab/lib/Utils.m`
- Package the legacy motion-generation utility class required by the local configs and generator.

### `experiment/run_MEG_lab/lib/createNewPoint.m`
- Package the local boundary-intersection and bounce helper used by trajectory generation.

### `experiment/run_MEG_lab/lib/plot_paths.m`
- Include the reusable trajectory visualization helper so path previews do not depend on files outside the lab bundle.

### `experiment/run_MEG_lab/participant_instruction_script.md`
- Include the participant-facing fixation, catch-response, blink-break, and comfort instructions alongside the executable bundle.

### `experiment/run_MEG_lab/runExpReadable.m`
- Add the primary lab launcher with subject/block controls, non-overwriting generation defaults, and a 3-degree runtime fixation default from the lab config.

### `experiment/run_MEG_lab/runExperiment.m`
- Include the full configurable acquisition launcher as an alternative entry point.

### `experiment/run_MEG_lab/runTestExperiment.m`
- Include the fake-eye-tracker test launcher for pre-session logic and display checks.

### `experiment/run_MEG_lab/stimuli_generation_V30_rescueTraject.m`
- Package config-driven observed/predicted trajectory generation and condition path plotting.

### `experiment/run_MEG_lab/trigger_codes_occlusion_v10_rescueTraject.md`
- Include the V20 condition, sequence, occlusion, catch, gaze/replay, and ESC trigger reference.

Suggested commit message:
`Add standalone MEG lab experiment bundle`

# 10 Jun 2026 14:59

## Refresh experiment procedures, participant guidance, and work tracking

### `experiment/MEG_lab_procedure.md`
- Update acquisition, preparation, consent, camera, projector, EyeLink, head-position, and participant-response instructions for the current lab workflow.
- Clarify the red/blue response mapping, fixation and blink-break guidance, optional EOG/ECG checks, and alternative frontal-sensor checks.
- Move trigger/response recovery steps into a troubleshooting section and mark unresolved lab-specific values with explicit `^^` update placeholders.

### `experiment/TODO.md`
- Change the active experiment inventory to the V30 generation, v23 scheduling, and v20 runtime stack.
- Record completed work for 40 unique paths, 15% catch rates, centered questions, fixation cues, blink breaks, block accuracy, EyeLink calibration/EDF handling, replay policy, and save/abort behavior.
- Keep current participant-session follow-ups visible and correct the replay description to reference v23-generated `TrialOrder`.

### `experiment/participant_instruction_script.md`
- Add a concise experimenter-facing participant script covering fixation, catch questions, response mapping, blink timing, block structure, head stability, and discomfort reporting.

Suggested commit message:
`Refresh experiment procedures and participant guidance`

# 10 Jun 2026 14:57

## Upgrade the canonical experiment stack to V30/v23/v20

### `experiment/CreateInputFiles_v21_rescueTraject.m` -> `experiment/CreateInputFiles_v23_rescueTraject.m`
- Advance the schedule builder to the V30/v23/v20 file contract and preserve deterministic three-run, block-aware scheduling.
- Generate percentage-based type-1 and type-2 catch plans, padded `TrialOrder`/`CatchPlan` arrays, compatibility `TrialStruct` data, and explicit schedule metadata.
- Add non-interactive subject/seed aliases, overwrite protection, validation of balanced condition pools, and detailed input/output documentation.

### `experiment/MoveDot1_experiment_occlusion_v18_rescueTraject.m` -> `experiment/MoveDot1_experiment_occlusion_v20_rescueTraject.m`
- Upgrade the runtime to consume V30 trajectories and v23 schedule/catch-plan artifacts.
- Add configurable jittered ITIs, fixation-only separation cues, centered catch prompts, timed blink breaks, message-screen cushions, and block-level catch-accuracy reporting.
- Restore configurable EyeLink acquisition, per-block EDF receipt, fixation-break replay handling, calibration debounce, and runtime `i`/`o` fixation-abort controls.
- Add strict MEG sync-test policy, full Psychtoolbox flip diagnostics, missed-flip guardrails, and precomputed straight-terminal path-band geometry.
- Update trigger behavior to the v10 map: condition onsets `71/81/91`, mutually exclusive catch results `201/202/113`, gaze/replay `150/151`, and ESC termination `152`.

### `experiment/fixation.bmp`
- Replace the EyeLink/runtime fixation bitmap with the intended updated 1024x768 acquisition asset shared by the canonical and standalone experiment bundles.
- Binary warning: bitmap changes are normally reviewed carefully because they cannot be line-diffed. The user explicitly confirmed this asset change is intentional, so it remains staged. To reverse that decision before commit: `git reset HEAD -- experiment/fixation.bmp`.
- No `.gitignore` rule is appropriate because this bitmap is a required runtime input rather than generated output.

### `experiment/lib/Config_runtime_v18_rescueTraject.m` -> `experiment/lib/Config_runtime_v20_rescueTraject.m`
- Expand the runtime preset for V20 with documented MEG, EyeLink, ITI, catch-layout, blink-break, path-band, sync-test, and flip-guard controls.
- Default the operational preset to 1000-1500 ms ITIs with a 100 ms delayed, 300 ms fixation-only cue while allowing schedule defaults through empty overrides.

### `experiment/lib/Config_schedule_CreateInputV21_MoveDotV18_rescueTraject.m` -> `experiment/lib/Config_schedule_CreateInputV23_MoveDotV20_rescueTraject.m`
- Update file naming and dependencies to the V30/v23/v20 stack.
- Set 40 base trials per condition, 15% catch rates, three-run block structure, message/calibration flow, straight path-band terminal controls, timing guardrails, and per-block EyeLink EDF policy.

### `experiment/lib/Config_stimuli_generation_V28_rescueTraject.m` -> `experiment/lib/Config_stimuli_generation_V30_rescueTraject.m`
- Raise generation to 40 unique sequence identities per condition and centralize V30 motion, fixation exclusion, curvature, deviance, and path-band geometry parameters.
- Keep fixed full-occlusion timing anchored at frame 130 through frame 190 under the 120 Hz MEG profile.

### `experiment/lib/plot_paths.m`
- Add a reusable path visualization helper with caller-provided axes, shared limits, time-based coloring, configurable save behavior, and input validation.

### `experiment/runExpReadable.m`
- Add a readable, checkpoint-backed launcher that keeps subject, block, overwrite, and fixation settings stable across child scripts that call `clearvars`.
- Run generation, schedule creation, input verification, and real MEG/EyeLink runtime in explicit order with skip-if-existing behavior.

### `experiment/runExperiment.m`
- Replace the legacy single-runtime call with a documented end-to-end acquisition launcher for V30 generation, v23 schedule creation, and v20 execution.
- Preserve runtime state through a temporary MAT checkpoint and force the intended real DataPixx/MEG and EyeLink acquisition profile.

### `experiment/runTestExperiment.m`
- Add a laptop-safe end-to-end launcher using fake gaze, configurable dry-run validation, local sync-test bypass, one-block scheduling, and cleanup of test outputs/state.
- Exercise the current ITI cue, blink-break, catch-layout, and message flow without requiring EyeLink hardware.

### `experiment/stimuli_generation_V28_rescueTraject.m` -> `experiment/stimuli_generation_V30_rescueTraject.m`
- Replace source-subject transforms with config-driven generation of paired nondeviant, deviant, and predicted trajectories.
- Enforce shared pre-deviance paths, effective-deviance checks, fixation-zone collision handling, fixed-frame path-band occlusion metadata, and deterministic subject seeding.
- Save observed/predicted MAT artifacts plus an optional three-condition path visualization generated through `plot_paths`.

### `experiment/trigger_codes_occlusion_v9_rescueTraject.md` -> `experiment/trigger_codes_occlusion_v10_rescueTraject.md`
- Document the V20 trigger map, 40-value sequence range, mutually exclusive catch-question outcomes, per-block EDF policy, and standardized ESC logging.

Suggested commit message:
`Upgrade experiment stack to V30 v23 and v20`

# 10 Jun 2026 14:56

## Ignore generated MEG, experiment, and Python artifacts

### `.gitignore`
- Ignore Python bytecode and cache directories repository-wide through `**/__pycache__/`, `*.py[cod]`, and `.pytest_cache/`.
- Ignore generated input/output trees for both standalone experiment bundles so subject schedules, recordings, trigger CSVs, and preview images cannot be staged accidentally.
- Exclude the MEG raw-data, derivatives, and temporary-work trees, including nested repositories and large FIF intermediates under `meg/tmp/`.
- Exclude generated PDF and SVG exports of the MEG pipeline diagram while retaining its editable Mermaid source and render script.

### `__pycache__/inspect_fif_report.cpython-313.pyc` (deleted)
- Remove the tracked CPython bytecode artifact; future regenerated copies are covered by the repository-wide cache rules.

### `__pycache__/trigger_events.cpython-313.pyc` (deleted)
- Remove the tracked trigger-tool bytecode artifact and prevent Python interpreter output from re-entering version control.

Suggested commit message:
`Ignore generated MEG and Python artifacts`

# 20 Apr 2026 11:13

## Stage rescueTraject migration, practice runner, trigger-check expansion, and artifact hygiene

### `.gitignore`
- Add default-ignore rules for generated/binary artifacts: `*.mat`, `*.png`, `*.csv`, `*.html`, `*.fif`.
- Add folder-level ignore rules for practice generated outputs: `experiment/practice/input_files/` and `experiment/practice/output_files/`.
- Add ignore rules for `triggers_check/oldTriggerChecks_asBakcup/*` while explicitly allowing the archived source script `triggers_check/oldTriggerChecks_asBakcup/compare_reported_vs_observed_triggers.py`.
- Keep ignore strategy aligned with the hygiene request so `paper/figures/*.png` stays out of the index unless explicitly forced.

### `experiment/CreateInputFiles_v20_threeRunsPerBlock_catch_blockResume.m` -> `experiment/CreateInputFiles_v21_rescueTraject.m`
- Rename the scheduler/input-builder entry point to the rescueTraject lineage and expand script-header documentation with purpose, data-flow assumptions, and MATLAB usage examples.
- Keep deterministic three-runs-per-block partitioning while documenting catch insertion policy and v21 output structures (`TrialStruct`, `CatchPlan`, `Schedule`, `TrialOrder`).

### `experiment/MEG_lab_procedure.md`
- Add an operator note about the end-of-block head-position popup and acceptable displacement range guidance (`<40/60 mm`).

### `experiment/MoveDot1_experiment_occlusion_v17_blockResume.m` -> `experiment/MoveDot1_experiment_occlusion_v18_rescueTraject.m`
- Rename the runtime experiment script to the v18 rescueTraject variant and update top-level documentation for restored EyeLink flow, replay triggers, flip diagnostics, and runtime override surface.
- Preserve condition/catch trigger behavior while documenting timing guardrails and calibration/recalibration controls.

### `experiment/TODO.md`
- Rewrite the TODO header to reflect the current active rescueTraject script/config set.
- Add open tasks for block-end accuracy printouts, EyeLink recalibration escape flow, and participant feedback about catch-trial density.
- Mark practice launcher work as done and reformat completed/open sections for the current workflow snapshot.

### `experiment/decreasing_speed.md` (deleted)
- Remove the legacy speed-scaling note file that referenced older config/script names no longer used in the active rescueTraject flow.

### `experiment/lib/Config_runtime_v17_blockResume.m` -> `experiment/lib/Config_runtime_v18_rescueTraject.m`
- Rename runtime preset script to match v18 rescueTraject and refresh configuration documentation blocks, defaults, and accepted override semantics.
- Keep preset-driven launch flow (config script then runtime script) explicit for batch and interactive operation.

### `experiment/lib/Config_schedule_CreateInputV20_MoveDotV17_blockResume.m` -> `experiment/lib/Config_schedule_CreateInputV21_MoveDotV18_rescueTraject.m`
- Rename schedule config class to v21/v18 rescueTraject naming and keep schedule/catch/message/timing policy centralized in one constant config surface.
- Preserve fixed three-run assumptions and catch parametrization while aligning usage docs to the new script family.

### `experiment/lib/Config_stimuli_generation_V28_rescueTraject.m`
- Add the dedicated V28 rescueTraject generation config class in active `experiment/lib` with fixed-frame occlusion controls, turn-window limits, and rescue-related parameters.
- Document one-dot assumptions, path-band geometry controls, and fixation-collision handling modes used by the new generator.

### `experiment/practice/CreateInputFiles_v21_rescueTraject.m`
- Add a practice-local copy of the v21 input-builder to keep practice generation isolated from main experiment directories.
- Keep same catch/schedule model but scoped to `experiment/practice` paths.

### `experiment/practice/MoveDot1_experiment_occlusion_v18_rescueTraject.m`
- Add a practice-local runtime copy for non-MEG practice execution with the same rescueTraject runtime behavior surface.
- Maintain trigger/debug/eye-tracker compatibility while constraining file IO to the practice tree.

### `experiment/practice/fixation.bmp`
- Add a fixed bitmap asset used by the practice-local EyeLink/visual setup path (`eyeTrackerFixBmpPath`).
- Keep it versioned as a required runtime asset for the isolated practice launcher.

### `experiment/practice/lib/Config_schedule_CreateInputV21_MoveDotV18_rescueTraject.m`
- Add practice-local schedule config with run/catch/message defaults tuned for practice execution.
- Keep class-level constants self-contained so practice does not depend on parent config files.

### `experiment/lib/Config_stimuli_generation_V27_blockResume.m` -> `experiment/practice/lib/Config_stimuli_generation_V28_rescueTraject.m`
- Move the prior active config lineage into the practice-local tree as part of the practice package.
- Update naming/content toward V28 rescueTraject controls inside practice scope.

### `experiment/practice/lib/Utils.m`
- Add practice-local `Utils` support class to decouple practice scripts from parent-path dependency resolution.

### `experiment/practice/run_practice_occlusion_v1.m`
- Add a dedicated top-level practice orchestrator with extensive script docs, safety preflight, repeat-loop controls, run1/run2 trial-count selection, and consolidated participant metrics output.
- Implement practice input caching/generation and runtime execution flow fully within `experiment/practice`.

### `experiment/stimuli_generation_V27_blockResume.m` -> `experiment/practice/stimuli_generation_V28_rescueTraject.m`
- Move the previous active generator lineage into the practice-local package and align naming/docs with V28 rescueTraject.
- Keep generation behavior available for practice input regeneration without mutating the active experiment tree.

### `experiment/stimuli_generation_V28_rescueTraject.m`
- Add the active config-driven V28 rescueTraject generator in `experiment/` with explicit turn-magnitude rescue fallback (`turnRescueScaleGrid`) before candidate rejection.
- Preserve fixed-frame occlusion timing, shared pre-deviance branch behavior, and metadata exports needed by runtime/analysis.

### `experiment/stimuli_generation_versions.md`
- Add a new V28_rescueTraject section documenting added turn-magnitude rescue behavior and per-trial rescue metadata fields.
- Update the “which version to use” guidance to include V28_rescueTraject for better CW/CCW balance under shared placement constraints.

### `experiment/trigger_codes_occlusion_v8_blockResume.md` -> `experiment/trigger_codes_occlusion_v9_rescueTraject.md`
- Rename trigger-code documentation file to v9 rescueTraject and keep trigger map aligned to the renamed runtime lineage.

### `paper/methods.md`
- Update methods text for the new deviant turn bounds (`[-81,-15] U [15,81]`) and add a curvature-unit clarification (`deg/frame`).
- Note schedule flexibility wording for total block count.

### `paper/plot_baseline_curvature_limit_paths.m`
- Add a standalone plotting script that visualizes post-deviance trajectory outcomes for curvature-window boundary values while enforcing a shared pre-deviance segment.
- Include documented workspace overrides, validation, and deterministic PNG export path.

### `paper/plot_deviant_turn_limit_paths.m`
- Add a complementary plotting script for deviant-turn boundary values with shared-prefix enforcement and translated suffix continuity.
- Keep plotting/export behavior parallel to the baseline-curvature plot helper for reproducible methods figures.

### `reports/sub04/sub04_report.html` (deleted)
- Remove tracked generated report artifact from `reports/sub04` to reduce repository noise from export outputs.

### `reports/sub04/sub04_run01_report.html` (deleted)
- Remove tracked run-level generated report artifact from `reports/sub04`.

### `triggers_check/compare_reported_vs_observed_triggers.md`
- Add usage and flag reference documentation for the single-run confrontation script, including explicit-path mode and auto-resolution behavior.

### `triggers_check/compare_reported_vs_observed_triggers.py`
- Generalize script language from V17-only to version-agnostic derivative handling.
- Add wildcard-based reported-CSV auto-resolution with optional `--variant` filtering and deterministic candidate ranking by parsed version and mtime.
- Refresh CLI help text and examples to rescueTraject naming and update transition warning plot-window default documentation.

### `triggers_check/compare_reported_vs_observed_triggers_subject.md`
- Add subject-level workflow documentation covering concatenated multi-run confrontation, memory model, output bundle, and scroller usage.

### `triggers_check/compare_reported_vs_observed_triggers_subject.py`
- Add a subject-level confrontation driver that discovers all subject CSV runs, maps block/run coverage to FIFs, concatenates reported/observed timelines, and emits one global report bundle.
- Keep memory footprint bounded by reading one FIF trigger channel at a time and storing compact window metadata.

### `oldTriggerChecks_asBakcup/compare_reported_vs_observed_triggers.py` -> `triggers_check/oldTriggerChecks_asBakcup/compare_reported_vs_observed_triggers.py`
- Move archived trigger-check script into the `triggers_check/oldTriggerChecks_asBakcup` location.
- Preserve file contents unchanged (`R100`) to keep archival provenance while consolidating trigger-check assets under one tree.

Suggested commit message:
`Migrate rescueTraject experiment flow and add practice tooling`

# 13 Apr 2026 12:27

## Stop tracking archived `oldies` source trees while keeping files on disk

### `.gitignore`
- Add a repository-wide archive rule `**/oldies/` so any `oldies` directory is treated as archival and excluded from future tracking by default.
- Keep existing artifact-specific ignore entries in place; this new rule broadens scope to archived source trees (not just generated files inside them).

### `experiment/oldies/**` and `simulations/oldies/**` (index-only removals)
- Untrack all currently tracked files under `experiment/oldies` and `simulations/oldies` via `git rm -r --cached`.
- This change removes archived sources from git history going forward while intentionally preserving the files on local disk for offline reference.
- Verified representative files still exist locally after untracking:
  `experiment/oldies/CreateInputFiles_v10.m`, `simulations/oldies/scripts/PE_simulation_diff.m`.

# 13 Apr 2026 12:24

## Ignore archived generated artifacts in experiment and trigger-check trees

### `.gitignore`
- Add folder-level ignore rules for archived/generated artifact trees so large legacy `.mat/.csv/.png` outputs and trigger-check report products are kept out of version control.
- Newly ignored paths:
  `experiment/oldies/input_files/`, `experiment/oldies/output_files/`,
  `oldTriggerChecks_asBakcup/experiment/`, `oldTriggerChecks_asBakcup/reports/`,
  `triggers_check/experiment/`, `triggers_check/reports/`, `triggers_check/__pycache__/`.
- Keep the ignore strategy directory-scoped (instead of broad extension globs) to avoid masking source files elsewhere in the repository.

# 13 Apr 2026 12:17

## Add one-dot occlusion v2 simulation scripts

### `simulations/scripts/one-dot/PE_simulation_diff_1Dot_occlusion_v2.m`
- Add a new one-dot occlusion PE simulation entry point (v2) with extensive script-level and section-level documentation, explicit usage examples, and clarified data-flow assumptions.
- Keep the canonical PE definition (`observed - predicted`) explicit while restricting the analysis surface to one-dot occlusion conditions and adding stricter input/version resolution behavior.
- Extend configuration controls for existing-results policy and PCR regression strategy handling so non-interactive runs can choose deterministic overwrite/reuse behavior.

### `simulations/scripts/one-dot/plot_occlusion_paths_1Dot_v2.m`
- Add a new minimal one-dot path-inspection helper (v2) that resolves observed input files by version priority and samples aligned trial identities across `always_visible`, `occluded_nondeviant`, and `occluded_deviant`.
- Render shared-axis three-panel path comparisons and save reproducible PNG outputs under the subject-specific one-dot occlusion output tree.

# 13 Apr 2026 12:16

## Archive simulations scripts tree into `simulations/oldies`

### `simulations/scripts/*.m` -> `simulations/oldies/scripts/*.m` (root scripts)
- Move the root simulation pipeline/debug scripts to `simulations/oldies/scripts` as pure path relocations (`R100`), preserving file content and behavior unchanged.
- Moved files:
  `PE_simulation_RDM_level_PE.m`, `PE_simulation_diff.m`, `PIPELINE_simulation.m`, `PIPELINE_simulation_posOnly.m`, `autocorr_leakage_sweep.m`, `barebone_pipeline.m`, `counterfactual_simulation.m`, `pipeline_recursive.m`, `pipeline_recursive_backup.m`, `toy_PE_diff.m`, `toy_direction.m`.

### `simulations/scripts/one-dot/*.m` -> `simulations/oldies/scripts/one-dot/*.m`
- Archive legacy one-dot simulation entry points and pipeline launchers under `simulations/oldies/scripts/one-dot` without modifying implementation details.
- Moved files:
  `PE_simulation_diff_1Dot.m`, `PE_simulation_diff_1Dot_occlusion_v1.m`, `plot_occlusion_paths_1Dot.m`, `run_pipeline.m`, `run_pipeline_occlusion.m`.

### `simulations/scripts/step_by_step/*.m` -> `simulations/oldies/scripts/step_by_step/*.m`
- Move step-by-step teaching/debug pipeline scripts to the simulations archive tree to keep the active `simulations/scripts` path focused on current workflows.
- Moved files:
  `pipeline_base.m`, `pipeline_base_nondeviant_dot1_minimal.m`.

### `simulations/scripts/tests/*.m` -> `simulations/oldies/scripts/tests/*.m`
- Move toy/debug test scripts to `simulations/oldies/scripts/tests` as archival snapshots.
- Moved files:
  `toy_PE_debug.m`, `toy_PE_debug_videos.m`, `toy_PE_debug_videos_displacement.m`.

# 13 Apr 2026 12:13

## Archive legacy experiment lineage under oldies and clean tracked trigger/output artifacts

### `experiment/oldies/CreateInputFiles_v14_threeRunsPerBlock_catch.m`, `experiment/oldies/CreateInputFiles_v15_threeRunsPerBlock_catch.m`, `experiment/oldies/CreateInputFiles_v16_threeRunsPerBlock_catch_NotWiredOccluder.m`, `experiment/oldies/CreateInputFiles_v17_threeRunsPerBlock_catch_PathBandOccluder.m`, `experiment/oldies/CreateInputFiles_v18_threeRunsPerBlock_catch_PathBandOccluder.m`, `experiment/oldies/CreateInputFiles_v19_threeRunsPerBlock_catch_eyeTrackerReplay.m`
- Add archived CreateInputFiles script snapshots (v14-v19) into `experiment/oldies` so prior schedule and catch-planning iterations remain recoverable without staying in the active root.

### `experiment/oldies/MoveDot1_experiment_occlusion_v8_runColorCueMessages.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v9_addedPracticeMode.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v10_practiceByRun.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v11_practiceByRun.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v12_practiceByRun.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v13_NotWiredOccluder.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v14_PathBandOccluder.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v15_PathBandOccluder.m`, `experiment/oldies/MoveDot1_experiment_occlusion_v16_eyeTrackerReplay.m`, `experiment/oldies/MoveDot1_experiment_vX copy.m`
- Add archived runtime script snapshots for the v8-v16 occlusion family plus the generic `MoveDot1_experiment_vX` copy, consolidating historical runtime variants under `experiment/oldies`.

### `experiment/oldies/lib/Config_runtime_v16_eyeTrackerReplay.m`, `experiment/oldies/lib/Config_schedule_CreateInputV14_MoveDotV6.m`, `experiment/oldies/lib/Config_schedule_CreateInputV14_MoveDotV7.m`, `experiment/oldies/lib/Config_schedule_CreateInputV14_MoveDotV8.m`, `experiment/oldies/lib/Config_schedule_CreateInputV14_MoveDotV9.m`, `experiment/oldies/lib/Config_schedule_CreateInputV15_MoveDotV9.m`, `experiment/oldies/lib/Config_schedule_CreateInputV15_MoveDotV10.m`, `experiment/oldies/lib/Config_schedule_CreateInputV15_MoveDotV11.m`, `experiment/oldies/lib/Config_schedule_CreateInputV15_MoveDotV12.m`, `experiment/oldies/lib/Config_schedule_CreateInputV16_MoveDotV13_NotWiredOccluder.m`, `experiment/oldies/lib/Config_schedule_CreateInputV17_MoveDotV14_PathBandOccluder.m`, `experiment/oldies/lib/Config_schedule_CreateInputV18_MoveDotV15_PathBandOccluder.m`, `experiment/oldies/lib/Config_schedule_CreateInputV19_MoveDotV16_eyeTrackerReplay.m`
- Add archived schedule/runtime config lineage files to keep historical parameter surfaces and run-structure assumptions versioned alongside their archived scripts.

### `experiment/oldies/lib/Config_stimuli_generation_V21.m`, `experiment/oldies/lib/Config_stimuli_generation_V23_NotWiredOccluder.m`, `experiment/oldies/lib/Config_stimuli_generation_V24_PathBandOccluder.m`, `experiment/oldies/lib/Config_stimuli_generation_V25_PathBandOccluder.m`, `experiment/oldies/lib/Config_stimuli_generation_V26_eyeTrackerReplay.m`, `experiment/oldies/input_files/Config_stimuli_generation_V22.m`
- Add archived generation config snapshots across V21-V26 plus the V22 input-file config snapshot for reproducible regeneration/debugging of old datasets.

### `experiment/oldies/stimuli_generation_V22.m`, `experiment/oldies/stimuli_generation_V23_NotWiredOccluder.m`, `experiment/oldies/stimuli_generation_V24_PathBandOccluder.m`, `experiment/oldies/stimuli_generation_V25_PathBandOccluder.m`, `experiment/oldies/stimuli_generation_V26_eyeTrackerReplay.m`
- Add archived generator script snapshots so prior trajectory-generation branches remain available in the `oldies` tree after active-path consolidation.

### `experiment/oldies/trigger_codes_occlusion_v6.md`, `experiment/oldies/trigger_codes_occlusion_v7_eyeTrackerReplay.md`
- Add archived trigger-code maps for older runtime versions to preserve historical trigger semantics in the same archive tree as the corresponding scripts/configs.

### `derivatives/triggers/sub04/*.csv`, `derivatives/triggers/sub98/*.csv`, `derivatives/triggers/sub99/*.csv` (deleted)
- Remove tracked trigger-check CSV artifacts for archived subjects/runs from `derivatives/triggers`, keeping trigger-data cleanup separate from active script sources.
- CSV paths are intentionally not enumerated one-by-one here to keep notes concise.

### `experiment/output_files/practice_SUB99_RUN01.mat` (deleted binary)
- Remove a tracked output artifact from `experiment/output_files` as part of output-tree cleanup.
- Warning (binary extension present in staged set): if this `.mat` deletion was unintended, unstage it with:
  `git reset HEAD -- experiment/output_files/practice_SUB99_RUN01.mat`
- `.gitignore` guidance: no additional rule needed because `experiment/output_files/` is already ignored.

# 13 Apr 2026 12:12

## Add trigger confrontation tooling for reported-vs-observed checks

### `triggers_check/compare_reported_vs_observed_triggers.py`
- Add the V17-focused trigger confrontation CLI that aligns reported MATLAB debug CSV events against observed FIF stim-channel windows using semi-global sequence alignment.
- Include anomaly extraction/reporting, transition-pulse collapse warnings, markdown summary generation, per-anomaly plots, and optional interactive scroller support for rapid trigger-channel inspection.
- Keep script-level documentation explicit for expected inputs, output artifacts, and default analysis windows so the workflow can run non-interactively from subject/block/run selectors.

### `oldTriggerChecks_asBakcup/compare_reported_vs_observed_triggers.py`
- Add an archived trigger-confrontation script snapshot preserving the earlier tool behavior for historical comparison and backward-reference debugging.

### `triggers_check/context.md`
- Add a local context manifest describing canonical legacy experiment/config references and the three-way trigger comparison goal (expected, reported, observed) used by the trigger-check workspace.

### `triggers_check/AGENTS.md`
- Add workspace-local coding and MATLAB invocation instructions aligned with repository AGENTS policies for documentation/commenting and absolute-path MATLAB execution.

# 13 Apr 2026 12:12

## Update experiment documentation and add simulation/methods notes

### `.gitignore`
- Add a root `*.gdoc` ignore rule so Google Docs link-export sidecar files are not accidentally tracked in future documentation updates.

### `experiment/MEG_lab_procedure.md`
- Add an explicit pre-session check item for confirming Polhemus power state in the "Before subject" preparation checklist.
- Keep checklist formatting aligned with the existing operations flow used during acquisition setup.

### `experiment/TODO.md`
- Replace the short legacy TODO header with a dated, structured list that separates active priorities from completed items and older deferred tasks.
- Document recent decisions and implementation outcomes for EyeLink calibration flow, replay behavior, abort messaging, save policy, and display color/mask adjustments.

### `experiment/timing.md`
- Add a dedicated timing-assessment note for the path-band runtime with MEG-oriented interpretation, current risks, and a prioritized hardening checklist.
- Record scope limits (static inspection vs live capture), trigger/flip alignment caveats, and acceptance criteria for production timing validation.

### `paper/methods.md`
- Add a methods draft covering V27 block-resume stimulus generation, occlusion geometry/timing anchors, schedule and catch-trial policies, message flow, and trigger-code mapping.
- Include explicit parameter values and assumptions used by the current runtime/generator stack so manuscript text stays synchronized with implementation.

### `simulations/report/PCR_full_ridge_timeline_explainer_Sub52_V27.md`
- Add a detailed explainer of the one-dot occlusion PCR `ridge_full_autocorr` strategy, including design-matrix construction, penalty behavior, and timeline-window interpretation for standard vs PE branches.
- Document metadata-backed frame/timeline values from the referenced Sub52 run so result interpretation remains reproducible.

# 13 Apr 2026 12:11

## Add block-resume v17/v20/v27 experiment pipeline and retire legacy active files

### `AGENTS.md`
- Remove the extra trigger-commentary alignment rule so repository-level coding guidance stays focused on documentation/commenting and local MATLAB invocation policy.

### `experiment/CreateInputFiles_v14_threeRunsPerBlock_catch.m` -> `experiment/CreateInputFiles_v20_threeRunsPerBlock_catch_blockResume.m`
- Promote the input-builder script to the v20 block-resume variant and align script naming with the new runtime/generation family.
- Keep deterministic three-runs-per-block scheduling while integrating explicit catch-planning metadata (including run-scoped catch type placement and expected response fields) into the generated TrialStruct artifacts.
- Expand script-level and section-level documentation so non-interactive/batch usage and workspace override behavior are explicit.

### `experiment/MoveDot1_experiment_occlusion_v17_blockResume.m`
- Add the new active v17 runtime entry point for the one-dot occlusion paradigm with block-resume support, catch handling, and schedule compatibility with the v20 TrialStruct/CatchPlan format.
- Restore and document eye-tracker related controls (including optional calibration gates, fixation-break detection, and replay scheduling) while keeping debug-mode and hardware-mode behavior configurable via workspace overrides.
- Extend timing instrumentation and runtime guardrails (flip diagnostics, missed-frame counters, and optional abort thresholds) to support MEG-focused timing validation workflows.
- Keep trigger emissions aligned to the new v8 trigger map, including explicit catch-question status closure behavior and replay/gaze-break markers.

### `experiment/MoveDot1_experiment_occlusion_v8_runColorCueMessages.m` (deleted)
- Remove the previous active v8 runtime script from the active experiment root now that v17 block-resume is the canonical runtime entry point.

### `experiment/lib/Config_runtime_v17_blockResume.m`
- Add a runtime preset script that exposes subject/session controls and all major v17 overrides (ITI jitter, debug overlays, eye-tracker modes, calibration choices, and path-band geometry tuning) in one reproducible launcher surface.
- Document accepted value ranges and default fallback behavior for each override so operators can run deterministic presets in interactive and batch sessions.

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV6.m` (deleted)
- Remove the old V14/V6 schedule config from active lib after migration to the new v20/v17 block-resume schedule class.

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV7.m` (deleted)
- Remove the old V14/V7 schedule config from active lib to avoid parallel active schedule definitions after v20 adoption.

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV8.m` (deleted)
- Remove the old V14/V8 schedule config from active lib and consolidate schedule ownership in the v20/v17 block-resume class.

### `experiment/lib/Config_schedule_CreateInputV20_MoveDotV17_blockResume.m`
- Add the canonical schedule/config class for the new runtime/generator chain, including fixed 3-run partitioning, catch-rate controls, message flow defaults, eye-tracker policy defaults, and timing guard parameters.
- Keep run-family color-cue and transition messaging controls co-located with catch/planning controls so runtime and input creation consume one shared schedule contract.

### `experiment/lib/Config_stimuli_generation_V21.m` (deleted)
- Remove the previous V21 generation config from active lib and retire the superseded config lineage from the active path.

### `experiment/lib/Config_stimuli_generation_V22.m` -> `experiment/lib/Config_stimuli_generation_V27_blockResume.m`
- Promote the generation config to the V27 block-resume variant and keep one-dot synthesis controls, fixed-frame occlusion anchors, and path-band geometry policy in a dedicated class.
- Preserve fixation-collision handling controls and deviance-curvature windows while aligning output naming and assumptions with the new runtime/input-builder family.

### `experiment/stimuli_generation_V22.m` -> `experiment/stimuli_generation_V27_blockResume.m`
- Promote the generator script to V27 block-resume naming and update script-level documentation to reflect the de-novo config-driven synthesis flow and exported metadata contracts.
- Keep fixed-frame occlusion timing and path-band metadata generation explicit, including terminal-style controls and compatibility with the v17 runtime.

### `experiment/trigger_codes_occlusion_v6.md` (deleted)
- Remove the outdated v6 trigger reference from the active root in favor of the new block-resume trigger map.

### `experiment/trigger_codes_occlusion_v8_blockResume.md`
- Add the block-resume trigger-code table documenting condition onsets, sequence identity range, occlusion event markers, catch question lifecycle, gaze/replay events, and dedicated ESC termination code.
- Clarify runtime semantics for mutually exclusive catch-question close triggers and preserve collision-free trigger-space expectations with dynamic sequence IDs.

# 28 Mar 2026 22:17

## Add V22 generator with fixation-collision controls

### `experiment/stimuli_generation_V22.m`
- Add a new active V22 config-driven one-dot occlusion generator that extends the V21 generation flow with explicit fixation-zone collision handling modes (`off`, `retry`, `move`).
- Keep fixed-frame occlusion timing and output metadata contracts aligned with the existing occlusion runtime/simulation pipeline while introducing shape-preserving translation fallback logic for fixation collisions.

### `experiment/lib/Config_stimuli_generation_V22.m`
- Add the dedicated V22 config class with fixation-collision parameters (`fixationCollisionMode`, exclusion radius, padding, direction/shift sampling controls) plus core one-dot motion and occlusion defaults.

### `experiment/lib/Config_stimuli_generation_V21.m`
- Expand parameter-level comments/examples to clarify generation assumptions and units for V21 controls without changing the underlying default numeric behavior.

### `experiment/stimuli_generation_versions.md`
- Extend the version guide with a V22 section documenting fixation-collision behavior, retained V21 features, and MATLAB usage guidance for selecting collision mode.

# 28 Mar 2026 22:17

## Introduce v8 run-color runtime and schedule config variants

### `experiment/MoveDot1_experiment_occlusion_v8_runColorCueMessages.m`
- Add a new active v8 occlusion runtime variant that keeps the run-color-cue flow and catch-trial architecture while advancing the runtime naming/version lineage beyond the archived v6/v7 scripts.
- Keep trigger and schedule assumptions consistent with the existing one-dot occlusion family (condition onsets, sequence identity pulse, occlusion event triggers, and catch-response triggers).

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV7.m`
- Add a V7 schedule/config snapshot for the three-runs-per-block catch pipeline so v7 runtime assumptions remain explicit and reproducible.
- Set defaults/documentation fields for block count, catch rates, question timing, and run-color cue behavior in one place.

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV8.m`
- Add the active v8 schedule/config class used by the new runtime, preserving the same schedule/catch contract while versioning config ownership to v8.

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV6.m`
- Expand section-level comments and parameter-level examples to clarify schedule/catch data flow and runtime assumptions without changing the configured numeric behavior.

# 28 Mar 2026 22:16

## Archive v6/V21 runtime lineage under oldies

### `experiment/MoveDot1_experiment_occlusion_v6_runColorCue.m` -> `experiment/oldies/MoveDot1_experiment_occlusion_v6_runColorCue.m`
- Move the v6 run-color-cue runtime out of the active experiment root into `oldies` as a pure archive relocation, preserving the script content for reproducibility/back-reference.

### `experiment/oldies/MoveDot1_experiment_occlusion_v7_runColorCueMessages.m`
- Add the v7 run-color-cue-with-messages runtime variant into the archive tree so intermediate runtime iterations remain recoverable without crowding the active experiment folder.

### `experiment/stimuli_generation_V21.m` -> `experiment/oldies/stimuli_generation_V21.m`
- Move the V21 config-driven one-dot generator from active root to `oldies`, keeping the full generator implementation available while preparing the active path for V22.

# 28 Mar 2026 15:23

## Experiment folder archival and active v6/V21 layout

### `experiment/CreateInputFiles_v14_threeRunsPerBlock_catch.m`
- Add the v14 catch-aware TrialStruct builder at the active experiment root, now wired to `Config_schedule_CreateInputV14_MoveDotV6` for schedule, catch-rate, and question timing controls.
- Preserve the three-runs-per-block data flow (always-visible in run 1; mixed occlusion conditions in runs 2/3) while adding deterministic `CatchPlan` metadata aligned to `TrialOrder` slots.

### `experiment/MoveDot1_experiment_occlusion_v6_runColorCue.m`
- Add the v6 runtime as the active occlusion experiment entry point, including block-aware scheduling, catch-trial question flow, sequence-identity trigger support, and run-family color cue counterbalancing.
- Keep trigger semantics aligned with the occlusion event map (`occlusion_start`, `occlusion_complete`, `occlusion_end_start`, `occlusion_end_complete`) and catch-response events used by DataPixx/MEG workflows.

### `experiment/lib/Config_schedule_CreateInputV14_MoveDotV6.m`
- Add the active v6 schedule config class for run partitioning, catch-rate controls, and run-color cue toggles used by both input-file creation and runtime execution.

### `experiment/lib/Config_stimuli_generation_V21.m`
- Add a dedicated V21 generator config class that isolates one-dot occlusion synthesis parameters, fixed-frame occlusion timing defaults, and output naming controls from the legacy multi-version `Config` class.

### `experiment/runExperiment.m`
- Add a small convenience launcher that pre-sets subject/block/debug run variables and calls the selected occlusion runtime from an absolute local path.

### `experiment/stimuli_generation_V21.m`
- Add the active V21 config-driven generator that synthesizes one-dot trajectories without source subject MAT transforms while preserving fixed-frame occlusion geometry and metadata compatibility expectations.
- Export observed condition sets plus predicted occluded-deviant outputs in the current one-dot format used by downstream simulation tooling.

### `experiment/trigger_codes_occlusion_v6.md`
- Add the v6 trigger reference table, including condition onsets, sequence-identity range, occlusion event pulses, catch-question lifecycle, and response/timeout events.

### `experiment/stimuli_generation_versions.md`
- Extend the version guide with V20/V21 capability notes, fixed-frame occlusion timing rules, and usage guidance for choosing transformed-source versus config-driven generation paths.

## Legacy experiment archive moves (`experiment` -> `experiment/oldies`)

### Legacy CreateInput script moves
- Move `experiment/CreateInputFiles_v10.m` to `experiment/oldies/CreateInputFiles_v10.m` without behavior changes.
- Move `experiment/CreateInputFiles_v11.m` to `experiment/oldies/CreateInputFiles_v11.m` without behavior changes.
- Move `experiment/CreateInputFiles_v12_temporary.m` to `experiment/oldies/CreateInputFiles_v12_temporary.m`, keeping compatibility-helper logic while updating internal script-name references.
- Add `experiment/oldies/CreateInputFiles_v13_threeRunsPerBlock.m` to preserve the prior v13 schedule builder in the archive area.

### Legacy MoveDot runtime moves
- Move `experiment/MoveDot1_experiment_vX.m` to `experiment/oldies/MoveDot1_experiment_vX.m` as a pure archive relocation.
- Move `experiment/MoveDot1_experiment_vX_occlusion_v1.m` to `experiment/oldies/MoveDot1_experiment_occlusion_v1.m` and align script-name references in comments.
- Add archived copies of intermediate occlusion runtime variants:
  `experiment/oldies/MoveDot1_experiment_occlusion_v2_threeRunsPerBlock.m`,
  `experiment/oldies/MoveDot1_experiment_occlusion_v3_blocksBreak.m`,
  `experiment/oldies/MoveDot1_experiment_occlusion_v4_catchTrials.m`,
  `experiment/oldies/MoveDot1_experiment_occlusion_v5_sequenceTriggers.m`.

### Legacy config class moves and archival
- Move `experiment/lib/Config.m` to `experiment/oldies/lib/Config.m`, keeping legacy generator constants available under the archive tree.
- Add `experiment/oldies/lib/Config_occlusion_schedule_v1.m` and `experiment/oldies/lib/Config_schedule_CreateInputV14_MoveDotV5.m` to retain earlier scheduling config versions alongside archived runtimes.

### Legacy stimulus generator moves
- Move legacy generators from the active root to `experiment/oldies/`:
  `stimuli_generation_v05.m`, `stimuli_generation_v06.m`, `stimuli_generation_v07.m`, `stimuli_generation_v08.m`, `stimuli_generation_v09.m`, `stimuli_generation_v10.m`, `stimuli_generation_v12.m`, `stimuli_generation_v13.m`, `stimuli_generation_v14.m`, `stimuli_generation_v15.m`, `stimuli_generation_v15_experimentalPathScale.m`, `stimuli_generation_v16_Displacement.m`, `stimuli_generation_v17.m`, `stimuli_generation_v18.m`, and `stimuli_generation_v19.m`.
- Add `experiment/oldies/stimuli_generation_V20.m` as the archived fixed-frame V20 generator.
- Normalize authorship comments in moved v09-v19 scripts to the maintained naming convention.

### Legacy trigger doc moves
- Move `experiment/trigger_codes.md` to `experiment/oldies/trigger_codes.md` with updated runtime-name references.
- Add archived trigger maps `experiment/oldies/trigger_codes_occlusion_v4.md` and `experiment/oldies/trigger_codes_occlusion_v5.md` next to their archived runtime variants.

### `experiment/triggers_report_sub03_copy_withNotes.pdf` (deleted)
- Remove the tracked trigger-report PDF from the active `experiment/` root as part of the folder cleanup; an unstaged copy remains in `experiment/oldies/` for local archival.

## One-dot occlusion simulation/report refresh

### `simulations/scripts/one-dot/PE_simulation_diff_1Dot_occlusion_v1.m`
- Add matrix-guide overlays keyed to occlusion event timing (full disappearance and first reappearance) and propagate resolved event metadata into reproducibility outputs.
- Keep guide-time conversion explicit for both standard and PE-cut analysis windows so plotted markers stay aligned after frame trimming.

### `simulations/scripts/one-dot/run_pipeline_occlusion.m`
- Extend the batch runner to branch cleanly across `dRSAtypeToRun` values (`PCR` vs `corr`), normalize strategy inputs, and avoid redundant per-strategy reruns for the `corr` branch.
- Update participant default and runtime logging so batch invocations print both dRSA branch and strategy context per run.

### `simulations/report/PE_simulation_diff_1Dot.md`
- Append a generated Subject 67 section (`<!-- SUBJECT:67:START/END -->`) with matrix embeds and path checks across observed/predicted/PE comparisons and corr/PCR strategy variants.
- Preserve report structure and section conventions so future scripted subject refreshes can replace bounded subject blocks in place.

# 19 Mar 2026 15:52

## Occlusion PE report refresh

### `simulations/report/PE_simulation_diff_1Dot.md`
- Append a generated Subject 70 section (`<!-- SUBJECT:70:START/END -->`) for the one-dot occlusion outputs under `Sub70_oneDot_occlusion`.
- Add matrix embeds and path-check lines for `occluded_deviant` across Observed/Predicted/PE comparisons, including correlation and PCR strategy variants (legacy, ARC1, Ridge Full, Ridge Tapered) with both common and separate colorbar layouts.
- Add matrix embeds and path-check lines for `occluded_nondeviant` observed-only comparisons across the same corr/PCR strategy set and colorbar layout variants.
- Preserve the existing report structure so section-level updates remain script-replaceable for future subject refreshes.

# 19 Mar 2026 15:52

## One-dot occlusion simulation pipeline

### `simulations/scripts/one-dot/PE_simulation_diff_1Dot_occlusion_v1.m`
- Add a dedicated one-dot occlusion PE simulation entry point that runs only `occluded_nondeviant` and `occluded_deviant`, with strict one-dot input validation and extensive runtime/parameter documentation.
- Keep canonical PE construction explicit (`observed - predicted`) and run deviant-only PE analyses against three neural targets (`neuralPE`, `neuralPredicted`, `neuralObserved`).
- Add post-reappearance enforcement for occluded-deviant PE streams, deriving cut points from per-trial metadata with controlled fallback behavior when metadata is missing.
- Support both `corr` and `PCR` dRSA modes and expose the full PCR regression-strategy surface (`baseline_pcr_border`, `ridge_full_autocorr`, `ridge_tapered_autocorr`, `ar1_prewhite_ridge`) with reproducibility-oriented defaults.
- Organize outputs under `simulations/output/SubXX_oneDot_occlusion/...`, including strategy-tagged filenames to avoid collisions with legacy PCR outputs.

### `simulations/scripts/one-dot/plot_occlusion_paths_1Dot.m`
- Add a lightweight plotting helper to compare shared trial identities across `always_visible`, `occluded_nondeviant`, and `occluded_deviant` with shared axis limits.
- Load one-dot center-relative paths from `MovDot_SubXX.mat`, align samples by `sequence` when available, and fall back to index-based sampling when no full sequence intersection exists.
- Save reproducible path-comparison figures into the occlusion output tree for quick geometry sanity checks.

### `simulations/scripts/one-dot/run_pipeline_occlusion.m`
- Add a batch runner that loops participants and PCR strategies, then calls `PE_simulation_diff_1Dot_occlusion_v1.m` with explicit per-run parameter initialization.
- Set default batch behavior to rerun/overwrite existing outputs (`existingResultsAction = 2`) for non-interactive pipeline use.

### `simulations/scripts/one-dot/run_pipeline.m`
- Extend the one-dot strategy list to include `baseline_pcr_border` alongside ridge and AR1 strategies so baseline PCR outputs are regenerated in the same batch workflow.

# 19 Mar 2026 15:51

## Occlusion TrialStruct compatibility helper

### `experiment/CreateInputFiles_v12_temporary.m`
- Add a temporary converter that builds `SubXX_TrialStruct.mat` from occlusion-era `MovDot_SubXX.mat` inputs without reintroducing catch-trial logic.
- Preserve occlusion condition traceability by copying condition label, source trial index, sequence identity, and occlusion-enabled metadata into each `TrialStruct` row.
- Build balanced run orders across `always_visible`, `occluded_nondeviant`, and `occluded_deviant` pools by truncating to the minimum shared count per condition.
- Emit a no-catch `Catch` struct and legacy-compatible `TrialOrder`/`BlockOrder` containers so older workflows that require `SubXX_TrialStruct.mat` can run unchanged.
- Add numeric input validation and overwrite guards (including optional non-interactive `overwriteExisting`) to keep batch runs deterministic.

# 19 Mar 2026 15:51

## One-dot occlusion experiment core

### `experiment/stimuli_generation_v18.m`
- Add a v18 occlusion dataset builder that transforms an existing one-dot source into three aligned condition sets: `always_visible`, `occluded_nondeviant`, and `occluded_deviant`.
- Add configurable deviance-centered timing controls (`deviance-X`, `deviance+Y`) and emit per-trial twocircle geometry plus alpha fallback visibility profiles.
- Export trigger-aligned occlusion event frames (`start`, `complete`, `end_start`, `end_complete`) and construct predicted occluded-deviant trials that match observed trajectories until reappearance before diverging.
- Keep strict one-dot validation and deterministic trial pairing so generated observed/predicted outputs remain aligned for downstream runtime and simulations.

### `experiment/stimuli_generation_v19.m`
- Add a v19 generator variant that enforces a matched pre-deviance branch for occluded deviant trials and splices in the deviant suffix with translation for position-continuous continuity at deviance.
- Preserve the v18 output contract (three observed conditions, predicted branch, twocircle/alpha metadata, and event-frame fields) while changing only the deviant-branch construction rule.
- Add non-interactive overwrite control (`overwriteExisting`) so batch regeneration can run without prompt loops when target files already exist.

### `experiment/MoveDot1_experiment_occlusion_v1.m`
- Add a dedicated one-dot occlusion runtime that consumes `MovDot_SubXX.mat` metadata from v18/v19, validates required occlusion fields, and enforces one-dot trajectory shape.
- Implement condition-specific trial-onset triggers (`31`, `41`, `51`) and occlusion event triggers (`111`, `112`, `114`, `115`) with per-frame logging for debug/replay workflows.
- Support two rendering paths (`twocircle` default, `alpha` fallback), balancing trial scheduling via `SubXX_TrialStruct.mat` when available and an internal balanced fallback otherwise.
- Save run outputs and optional debug trigger CSVs under `experiment/output_files` with explicit no-overwrite guards.

### `experiment/trigger_codes.md`
- Add the canonical trigger-code table for `MoveDot1_experiment_occlusion_v1.m`, including condition onset pulses and the four occlusion event pulses.
- Clarify that event frames come from generator metadata and that condition and event trigger emissions are independent during a trial.

### `experiment/stimuli_generation_versions.md`
- Extend the version guide with v18/v19 selection guidance and capability summaries, including occlusion timing controls, twocircle metadata, event-frame exports, and v19 matched pre-deviance behavior.
- Document how v18/v19 differ from earlier versions so users can choose between de-novo generation and occlusion-paradigm transformation workflows.

# 17 Mar 2026 15:48

## Report tracking and ignore scope

### `.gitignore`
- Narrow the ignore scope from `simulations/report/` to `simulations/report/old/` so report deliverables under `simulations/report/` can be tracked while the archived `old` subtree stays ignored.
- Add a root-level `/$MDN` rule to prevent committing the accidental markdown-dump artifact created in the repository root.
- Add `simulations/report/**/*.pdf` to keep generated report exports out of version control while retaining markdown report sources.
- Keep the existing generated-output excludes (`simulations/output/`, debug artifacts, and archived PDFs under `simulations/report_old/**/*.pdf`) in place.

### `simulations/report/PE_simulation_diff_1Dot.md`
- Add a generated one-dot PE matrix report for Sub71, Sub72, and Sub73 with a top index and anchor links for deviant, nondeviant, and post-deviant sections.
- Enumerate matrix embeds across Observed/Predicted/PE comparisons and corr/PCR variants (legacy, ARC1, Ridge Full, Ridge Tapered), including per-image path-check strings that mirror the expected output tree.
- Provide a text-first report source that can be diffed and regenerated independently from the PDF export.

### `simulations/report/PE_simulation_diff_1Dot_backup.md`
- Add a backup markdown copy of the same generated report content as `PE_simulation_diff_1Dot.md` to preserve a second text snapshot in the report folder.
- Keep the backup in markdown form so it remains inspectable and diffable in the same way as the primary report source.

## Legacy script removal

### `simulations/report_old/generate_PIPELINE_simulation_report.m` (deleted)
- Remove the tracked MATLAB Report Generator script from the legacy `simulations/report_old` path.
- Retire this archived tracked copy now that report tracking is being redirected through the active `simulations/report` area and its scoped ignore rules.

# 16 Mar 2026 15:44

## Report generator archive

### `simulations/report/generate_PIPELINE_simulation_report.m` -> `simulations/report_old/generate_PIPELINE_simulation_report.m`
- Move the MATLAB Report Generator entry point out of the active `simulations/report` tree into `simulations/report_old` without changing its implementation.
- Preserve the script source in the archive location so older report-generation workflows remain recoverable while signaling that the file is no longer part of the active reporting path.

# 16 Mar 2026 15:44

## Diagnostic and tutorial scripts

### `simulations/scripts/autocorr_leakage_sweep.m`
- Add a dedicated analysis-side sweep that rescales path extent and curvature dynamics, recomputes position autocorrelation, and summarizes how leakage width changes across the tested parameter grid.
- Save CSV, Markdown, MAT, and heatmap outputs so stimulus-shape hypotheses can be compared quantitatively before changing generation code.

### `simulations/scripts/barebone_pipeline.m`
- Add a compact real-stimulus dRSA demo that runs both `corr` and `PCR` on dot-position streams using the original concatenate/subsample/dRSA structure.
- Keep the script focused on minimal end-to-end data flow so the core dRSA API can be inspected without the complexity of the full PE pipelines.

### `simulations/scripts/step_by_step/pipeline_base.m`
- Add a stripped-down step-by-step baseline script for inspecting the dRSA building blocks in isolation.
- Keep it separate from the real-data minimal examples so the introductory walkthrough remains easy to modify during debugging.

### `simulations/scripts/step_by_step/pipeline_base_nondeviant_dot1_minimal.m`
- Add a documented minimal real-data example that extracts nondeviant dot1 trials, builds position and direction streams, and runs exactly two small dRSA tests for direct inspection.
- Support both one-dot and two-dot input formats so the walkthrough remains usable across the newer stimulus-generation variants.

### `simulations/scripts/tests/toy_PE_debug.m`
- Add a compact PE debugging script for inspecting observed, predicted, and PE streams around the deviant segment without running the full production pipeline.
- Keep the script focused on geometric and representational sanity checks that help diagnose post-deviant behavior before larger batch runs.

### `simulations/scripts/tests/toy_PE_debug_videos.m`
- Add a video-oriented PE debug script that synchronizes path evolution and matrix-level diagnostics across post-deviant time for richer visual inspection.
- Separate the video export workflow from the still-image debug script so heavier derivative generation remains optional.

### `simulations/scripts/tests/toy_PE_debug_videos_displacement.m`
- Add a displacement-aware variant of the PE video debugger so observed-path jump manipulations can be inspected with the same synchronized visual tooling.
- Keep the displacement-specific logic isolated from the generic PE video debugger to avoid conflating the two test cases.

# 16 Mar 2026 15:43

## Presentation source

### `control-center/mock_presentation/README.md`
- Add a source-level README that explains the reproducible slide-deck workflow, the expected figure outputs, and the acceptance checks for rebuilding the presentation assets and PowerPoint export.
- Document the relationship between the markdown deck, figure-generator script, and generated outputs so the presentation can be refreshed from source rather than edited only in binary form.

### `control-center/mock_presentation/build_presentation.sh`
- Add a one-command build script that regenerates figure assets for a selected subject and renders the markdown deck to a PowerPoint file with `pandoc`.
- Resolve paths relative to the script directory so the build can be launched from any working directory without hand-editing output paths.

### `control-center/mock_presentation/generate_expected_figures.py`
- Add the figure-generation script that builds both real trajectory-derived assets and synthetic expected-pattern illustrations for the MEG meeting deck.
- Keep the figure generation separate from the markdown slides so presentation assets can be refreshed independently when the source data or expected-pattern sketches change.

### `control-center/mock_presentation/meg_meeting_slides.md`
- Add the markdown source for the meeting presentation, covering literature framing, hypotheses, design, trigger mapping, planned dRSA analysis, and the expected result figures.
- Reference the local `figures/` assets directly so the deck content remains editable in text form while the binary presentation output stays reproducible.

# 16 Mar 2026 15:43

## Planning documents

### `control-center/experiment_blueprint.md`
- Expand the blueprint from a short work guide into a broader design document that separates non-deviant and deviant hypotheses, adds concrete interpretation examples, and records multiple open design choices that still need to be locked before data collection.
- Rename the primary manipulation framing from "predictability" to explicit non-deviant versus deviant trial structure, and connect that framing to the current Config defaults and analysis constraints.
- Flesh out the "minimal items to lock" section with more concrete hypothesis/contrast definitions, including matched post-onset deviant comparisons and optional context/attention extensions.

### `control-center/objectives.md`
- Add a standalone project checklist that separates experiment definition, presentation preparation, implementation, piloting, and full data collection milestones.
- Include a long-form "Deep Research" prompt that formalizes the current theory-audit request, expected terminology, and the boundary conditions for rewriting the hypothesis set without redesigning the paradigm.

# 16 Mar 2026 15:42

## Stimulus generation controls and version expansion

### `experiment/lib/Config.m`
- Update the current generation defaults toward a smaller local batch size and a higher baseline curvature, while also switching the deviant-curvature path to random post-onset curvature sampling.
- Add explicit interval-based controls for baseline curvature, post-deviant curvature, signed deviant turns, and deviant displacement geometry so newer generators can constrain path families directly from Config instead of hard-coding ranges in each script.
- Extend the `likelihood` config struct with the new signed-turn and displacement fields so existing stimulus-type selection code can access the added controls through the same configuration object.

### `experiment/stimuli_generation_v14.m`
- Extend v14 to support explicit signed deviant-turn windows for likelihood-mode generation while preserving the existing proxy-gated, boundary-safe generation structure.
- Keep the saved reproducibility metadata aligned with the new turn-window controls so generated inputs can be audited against the configured deviant schedule.

### `experiment/stimuli_generation_v15.m`
- Add a new boundary-safe two-dot generator that introduces explicit initial/deviant curvature windows on top of the v14 signed-turn and dRSA-proxy-gate logic.
- Preserve the predicted no-deviant export path and the existing output contract so downstream experiment and simulation code can consume v15 outputs without loader changes.

### `experiment/stimuli_generation_v15_experimentalPathScale.m`
- Add an experimental v15 variant that shortens path geometry through a local path-scaling manipulation while keeping the same broad generation pipeline and saved outputs.
- Keep the variant separate from the main v15 script so experimental geometry changes do not silently alter the canonical generator.

### `experiment/stimuli_generation_v16_Displacement.m`
- Add a v16 generator that layers deviant-onset displacement control onto the existing turn and curvature manipulations, with configurable annular-sector sampling and mode-based overrides.
- Preserve the dRSA-proxy gate and no-deviant predicted export so displacement can be tested without breaking the analysis-oriented output structure.

### `experiment/stimuli_generation_v17.m`
- Add a one-dot generator that carries forward the v16 deviant controls and proxy gate while collapsing the path representation to a single observed dot.
- Document and export the one-dot output shape explicitly so downstream one-dot analyses can rely on a stable `frames x 2` coordinate contract.

### `experiment/stimuli_generation_versions.md`
- Expand the input-file audit table with the newer curvature, turn, displacement, gate, and onset fields so subject-level generation settings can be compared from the saved inputs alone.
- Add capability notes for v15, v15 experimental path scaling, v16 displacement, and v17 one-dot generation to explain when each generator should be used and what each adds relative to the base behavior.

# 16 Mar 2026 15:41

## dRSA PCR backends and one-dot workflow

### `simulations/functions/dRSA_PCR.m`
- Replace the legacy PCR-only implementation with a strategy-driven autocorrelation-removal backend that supports the existing PCR+border path plus `ridge_full_autocorr`, `ridge_tapered_autocorr`, and `ar1_prewhite_ridge`.
- Add helper blocks for explicit autocorrelation-regressor construction, lag-weighted ridge penalties, pooled AR(1) estimation/prewhitening, and first-predictor extraction so the target dRSA regressor is estimated consistently across strategies.
- Expand the function header and section comments to document the strategy options, parameter semantics, and the row-wise regression flow used when converting model/neural RDM streams into dRSA matrices.

### `simulations/functions/dRSA_coreFunction.m`
- Extend the PCR parameter contract with strategy-selection and regularization controls (`RegressStrategy`, `RidgeLambdaFactor`, `TaperSigmaFactor`, `StandardizePredictors`, `AR1Clip`, `AutoPenaltyStrength`, `NormalizeAutoPenaltyPerRow`).
- Validate the new strategy names and guard the weighted-penalty settings so downstream PCR calls fail early on unsupported configurations instead of silently reverting to legacy behavior.

### `simulations/functions/dRSA_computeRDM.m`
- Add an optional debug path that reconstructs and plots one square RDM sample from the condensed vector form, with configurable sample index and optional PNG export.
- Keep the debug branch gated behind `params.debug` so standard pipeline runs preserve the original computation path unless debugging is explicitly enabled.

### `simulations/scripts/PE_simulation_diff.m`
- Broaden the script from corr-only usage to a `dRSAtypeToRun` switch that accepts both `corr` and `PCR`, while keeping the shared simulation and export path intact.
- Add top-level debug parameters for the new `dRSA_computeRDM` snapshot path and keep the PCR parameter struct populated so the updated core functions can be exercised from the legacy PE workflow.
- Update the script commentary to reflect the new branch selection semantics and keep the existing-results rerun path compatible with the expanded dRSA modes.

### `simulations/scripts/one-dot/PE_simulation_diff_1Dot.m`
- Add a dedicated one-dot PE simulation pipeline with full script-level documentation, explicit usage examples, and section comments that trace data flow from input loading through dRSA/PE computation to organized output export.
- Support both `corr` and `PCR` execution, including the new PCR regression strategies, deviant-only analysis controls, and separate plotting control via `cutPostDevPlot` so analysis cropping and visualization cropping can diverge cleanly.
- Persist the expanded reproducibility metadata and strategy tags into output naming so one-dot comparison runs can coexist without overwriting each other.

### `simulations/scripts/one-dot/run_pipeline.m`
- Add a documented helper that batch-runs the one-dot PCR workflow across participants and regression strategies using explicit parameter overrides.
- Keep the helper in function form so the script-level `clearvars -except` in the main pipeline does not wipe the outer batch-loop variables during repeated runs.

## Documentation

### `simulations/documentation/dRSA_autocorr_regression_strategies.md`
- Add a focused write-up of the available PCR autocorrelation-removal strategies, their expected artifact behavior, and the rationale for choosing ridge-based implementations over lasso/elastic-net variants in this repo.
- Tie the conceptual explanation back to the concrete code paths so the strategy behavior can be audited against the current implementation.

### `simulations/documentation/dRSA_type_PCR_explained.md`
- Add a beginner-oriented, code-referenced walkthrough of `dRSAtype='PCR'` versus `dRSAtype='corr'`, including purpose, data flow, and the main parameter switches that affect behavior.
- Include follow-the-code guidance so readers can trace the high-level conceptual differences down to the specific simulation and function entry points.

# 16 Mar 2026 15:41

## Repo hygiene

### `.gitignore`
- Add ignore rules for generated mock-presentation outputs so source slide materials can be tracked without pulling `.pptx` exports, figure renders, or temporary Office lock files into commits.
- Ignore dRSA test snapshot and derivative export folders under `simulations/scripts/tests` to keep local debugging runs from polluting the working tree with generated PNG and MP4 artifacts.
- Ignore archived report PDFs under `simulations/report_old` and local shell history files so only source documents and report-generation scripts remain visible to git.

# 25 Feb 2026 05:21

## Agent workflow documentation

### `AGENTS.md`
- Add a local-tooling section that documents the absolute MATLAB R2020a binary path on this machine and makes absolute-path invocation the canonical command style.
- Provide non-interactive MATLAB execution examples for standard runs and Apple Silicon fallback (`arch -x86_64`) so scripted checks use a deterministic launch path.
- Keep the new guidance adjacent to existing documentation/commenting requirements so agent runs follow repo-specific constraints in one place.

# 25 Feb 2026 01:09

## Control-center docs and repo hygiene

### `.gitignore`
- Add `control-center/presentation/` to keep generated slides, figure artifacts, and temporary office-lock files out of version control by default.
- Preserve the existing ignore rules for other experiment/simulation generated outputs.

### `control-center/experiment_blueprint.md`
- Add a comprehensive experiment work guide that links scientific hypotheses to concrete implementation files, trigger references, stimulus constraints, and analysis commitments.
- Document decision points, quality gates, and canonical source-of-truth paths so future design changes can be traced directly to code/config files.

### `control-center/experiment_description.md`
- Add an early-stage brainstorming companion note that captures theoretical motivation and open design questions while pointing readers to the blueprint as the maintained canonical document.

# 25 Feb 2026 01:09

## Hypothesis testing script defaults

### `simulations/scripts/PE_simulation_diff.m`
- Change the local default `participantNumber` from 87 to 80 so quick runs target the Sub80 dataset without manual edits.
- Keep all other participant/configuration defaults unchanged.

# 25 Feb 2026 01:08

## Simulation plotting semantics

### `simulations/scripts/PE_simulation_RDM_level_PE.m`
- Correct dRSA matrix axis labels so x now reports neural/data time and y reports model time, matching the `corr(modelRDM, neuralRDM)` convention used by the computation.
- Add an inline clarification comment at the plotting point to keep future matrix-layout edits aligned with the intended axis semantics.

### `simulations/scripts/PE_simulation_diff.m`
- Apply the same axis-label correction to both the main dRSA matrix plots and the PE-target matrix plots so figure annotations match the underlying matrix dimensions.
- Add matching inline comments documenting the x/y convention in both plotting blocks to reduce interpretation mistakes when comparing figures across scripts.

### `simulations/scripts/PIPELINE_simulation.m`
- Swap the dRSA x/y label assignment to align axis text with neural-time on x and model-time on y.
- Add the same matrix-orientation comment used in the other scripts for consistent plotting semantics.

### `simulations/scripts/PIPELINE_simulation_posOnly.m`
- Correct the sample-based dRSA axis labels so x refers to neural/data time samples and y refers to model-time samples.
- Add an inline orientation note to keep the pos-only pipeline aligned with the full pipeline's dRSA interpretation.

### `simulations/scripts/counterfactual_simulation.m`
- Update dRSA matrix labels to the neural-time-x/model-time-y convention and document the convention inline.

### `simulations/scripts/pipeline_recursive.m`
- Correct dRSA axis labels in both per-condition and averaged matrix plotting blocks to prevent swapped-time interpretation in recursive outputs.
- Add inline comments in both blocks describing the `corr(modelRDM, neuralRDM)` axis mapping.

### `simulations/scripts/pipeline_recursive_backup.m`
- Apply the same dRSA axis-label correction and inline orientation comment in the backup recursive plotting path.

# 25 Feb 2026 01:08

## Experiment configuration and stimuli generation

### `experiment/lib/Config.m`
- Reduce `trialsPerCondition` from 1000 to 100 to shorten local generation/runtime loops while preserving the existing parameter structure.

### `experiment/stimuli_generation_v13.m`
- Add a new v13 generator that keeps the v12 output contract (`xySeqs`, `xySeqsPredicted`, `Cfg`, `repro`) while replacing whole-trial boundary rejection with relative-path generation plus feasible start-position placement.
- Add geometry-aware curvature-floor protection, optional no-deviant-baseline min-distance enforcement, and hard attempt limits so incompatible parameter sets fail explicitly instead of looping indefinitely.
- Expand script-level documentation and section comments with usage examples, data-flow notes, assumptions, and explicit output semantics aligned with experiment/runtime and simulation input builders.

### `experiment/stimuli_generation_v14.m`
- Add a new v14 generator that builds on v13 and introduces an optional dRSA-proxy-aware trial gate to reduce position-vs-direction cross-model coupling by accepting only candidates that improve a configurable proxy score.
- Keep within-trial curvature constant (except optional deviant-point modulation), expose proxy-gate tuning parameters, and persist the gate settings/failure counts into reproducibility metadata.
- Preserve the same saved artifacts (`xySeqs`, predicted deviant baselines, `Cfg`, `repro`) so downstream experiment and simulation workflows continue to load the expected structure.

## Stimuli generation documentation

### `experiment/stimuli_generation_versions.md`
- Extend the input-file parameter table with `randomizeCurvatureOnDeviant` and `deviantCurvatureRange`, and add additional subject rows (including Sub80-Sub86) so stored inputs can be audited against deviant-curvature settings.
- Add a dedicated table parameter guide section documenting path-level effects for each tracked parameter, with concrete examples to align documentation with generation behavior.
- Document new capability guidance for selecting v13 (boundary-safe placement and feasibility safeguards) vs v14 (dRSA-proxy-gated selection for reduced residual position-direction coupling).

# 24 Feb 2026 10:59

## Simulation folder reorganization

### `simulations/compute_avg_path_distance.m` -> `simulations/functions/compute_avg_path_distance.m`
- Move the path-distance helper into `simulations/functions` to align with the function/script split.
- Keep function signature and implementation unchanged so existing callers preserve identical outputs and assumptions.

### `simulations/PIPELINE_simulation.m` -> `simulations/scripts/PIPELINE_simulation.m`
- Move the full pipeline script into `simulations/scripts` and update top-level usage examples to the new location.
- Patch dependency/data path resolution to `scriptDir -> simDir -> repoRoot`, then use `simDir` for `input`, `output`, `functions`, and `debug` paths so runtime behavior remains stable after relocation.
- Keep `compute_avg_path_distance(...)` usage unchanged while resolving it via `addpath(fullfile(simDir, 'functions'))`.

### `simulations/PIPELINE_simulation_posOnly.m` -> `simulations/scripts/PIPELINE_simulation_posOnly.m`
- Move the position-only pipeline script into `simulations/scripts` and update header usage examples accordingly.
- Apply the same `scriptDir -> simDir -> repoRoot` path bootstrap and route `input/output/functions/debug` lookups through `simDir` to preserve previous behavior from the new folder.

## Cross-script reference cleanup

### `simulations/scripts/counterfactual_simulation.m`
- Update documentation references to the new pipeline paths under `simulations/scripts`.
- Update addpath example to use the scripts folder location.

### `simulations/scripts/PE_simulation_RDM_level_PE.m`
- Update documentation references to `simulations/scripts/PIPELINE_simulation*.m`.
- Update addpath example to `simulations/scripts`.

### `simulations/scripts/toy_direction.m`
- Update addpath example to `simulations/scripts` to match the new script location.

# 20 Feb 2026 11:53

## Simulation folder reorganization

### `simulations/build_movdot_simulation_inputs.m` -> `simulations/functions/build_movdot_simulation_inputs.m`
- Move the simulation-input builder into `simulations/functions` and keep it callable by name from all simulation scripts.
- Update internal path bootstrap (`scriptDir` -> `simDir` -> `repoRoot`) so default output still targets `simulations/input` and repository-relative path resolution remains unchanged after the move.
- Refresh header usage examples so they point to the new function location.

### `simulations/build_simulation_report_assets.m` -> `simulations/functions/build_simulation_report_assets.m`
- Move report-asset generation logic into `simulations/functions` as requested.
- Adjust startup path wiring to add the moved script directory plus `simulations` and `simulations/debug`, preserving access to dRSA/debug helpers and stable repo-root resolution.
- Update usage examples in the header to the new path.

### `simulations/generate_PIPELINE_simulation_report.m` -> `simulations/report/generate_PIPELINE_simulation_report.m`
- Move report generation into `simulations/report` and keep output under `simulations/report/subXX`.
- Fix moved-script path resolution by setting `simDir` to the parent `simulations` folder and re-adding `simulations/functions` so dependencies remain discoverable.
- Update help text/error guidance to reference `simulations/functions/build_simulation_report_assets.m`.

### `simulations/hypothesis_testing/PE_simulation_RDM_level_PE.m` -> `simulations/scripts/PE_simulation_RDM_level_PE.m`
- Relocate the RDM-level PE hypothesis-testing script under the new `simulations/scripts` folder with no behavior changes.

### `simulations/hypothesis_testing/PE_simulation_diff.m` -> `simulations/scripts/PE_simulation_diff.m`
- Relocate the signal-level PE simulation script into `simulations/scripts`.
- Update usage examples to replace `simulations/hypothesis_testing` with `simulations/scripts`.

### `simulations/hypothesis_testing/counterfactual_simulation.m` -> `simulations/scripts/counterfactual_simulation.m`
- Relocate the counterfactual simulation script into `simulations/scripts` with logic unchanged.

### `simulations/pipeline_recursive.m` -> `simulations/scripts/pipeline_recursive.m`
- Move the recursive pipeline script into `simulations/scripts`.
- Patch moved-script dependency wiring to use `simDir` for `input`, `output`, `functions`, and `debug` folders so runs still resolve the same resources after relocation.

### `simulations/scripts/pipeline_recursive_backup.m`
- Add the backup recursive pipeline script under `simulations/scripts` and align its path bootstrap with the moved-folder structure (`scriptDir` -> `simDir` -> `repoRoot`).

### `simulations/hypothesis_testing/toy_PE_diff.m` -> `simulations/scripts/toy_PE_diff.m`
- Relocate the toy PE visualization script into `simulations/scripts`.
- Update usage examples to the new folder while keeping execution logic unchanged.

### `simulations/toy_direction.m` -> `simulations/scripts/toy_direction.m`
- Relocate the toy direction modeling script into `simulations/scripts` without changing computation behavior.

# 20 Feb 2026 11:20

## Hypothesis testing scripts

### `simulations/hypothesis_testing/PE_simulation_diff.m`
- Add a full signal-level PE pipeline where PE is computed as observed minus predicted streams before dRSA (`PE_signal = observed - predicted`), keeping standard counterfactual outputs in the same run.
- Add deviant-only PE model families (position/direction, dot1/dot2) evaluated against three neural targets (`neuralPE`, `neuralPredicted`, `neuralObserved`) with dedicated matrix exports.
- Enforce post-deviance handling for PE analyses even when global `cutPostDev` is disabled, and track PE-specific cut diagnostics in reproducibility metadata.
- Add non-interactive existing-results behavior control (`existingResultsAction`) for batch rerun/reuse/error workflows.

### `simulations/hypothesis_testing/PE_simulation_RDM_level_PE.m`
- Add an RDM-level PE analysis variant that complements the signal-level workflow and expands hypothesis-testing coverage with a separate, documented script entry point.
- Preserve the organized output conventions and reproducibility-oriented script structure used across the updated simulation stack.

### `simulations/hypothesis_testing/toy_PE_diff.m`
- Add a compact visualization utility that overlays observed, predicted, centered-PE, origin-anchored PE, and PE vector fields for a small deviant trial subset.
- Include optional post-deviance cropping and deterministic export under `simulations/output/subXX/paths` to support quick geometric sanity checks of PE behavior.

# 20 Feb 2026 11:20

## Hypothesis testing scripts

### `simulations/hypothesis_testing/counterfactual_simulation.m`
- Add optional post-deviance truncation (`cutPostDev`, `deviantOnset`) and apply it consistently to condition-level path previews, observed trial streams, and deviant predicted streams.
- Move exports to structured output folders with explicit matrix subgrouping: nondeviant (`commCbar`, `sepCbar`) and deviant (`base` vs `predicted`, each with `commCbar`/`sepCbar`).
- Add backward-compatible reuse checks that can load legacy flat results files while writing new exports to the organized directory layout.
- Persist new repro metadata (`cutPostDev`, `deviantOnset`, `cutFrame`) and include these fields in existing-results parameter matching diagnostics.
- Add local safety helpers for directory creation and guarded figure export to avoid crashes when optional figure handles are invalid.

# 20 Feb 2026 11:20

## Simulation pipeline scripts

### `simulations/PIPELINE_simulation.m`
- Reorganize outputs into a nested subject/condition/dRSA folder tree (`results`, `matrices`, `diagonal`) and store all-condition path overviews under `subXX/paths`.
- Normalize parsed condition labels by stripping `run_default` tokens so output filenames/folders remain stable across generated input naming variants.
- Add a local directory-creation helper and route every save/print call through the new structured output paths.

### `simulations/PIPELINE_simulation_posOnly.m`
- Re-enable the previously disabled save block and gate it behind `saveOutputs` so runs can be compute-only or export-enabled without editing code blocks.
- Apply the same structured output layout and condition-label cleanup used by the full pipeline, including consistent optional shuffle/resample tags.
- Persist reproducibility metadata to the organized results path while keeping all-condition path exports in subject-level `paths` folders.

## Simulation input building and diagnostics

### `simulations/build_movdot_simulation_inputs.m`
- Cast dot-path arrays and center shifts to `double` before center-relative recentering to avoid mixed-class arithmetic failures in `bsxfun`.
- Document the recentering type-conversion assumption directly in the function header assumptions list.

# 20 Feb 2026 11:19

## Experiment configuration and stimuli generation

### `experiment/lib/Config.m`
- Enable `isCurvFactorRand` by default so per-trial curvature magnitude can vary without script-level edits.
- Add `randomizeCurvatureOnDeviant` and `deviantCurvatureRange` to support configurable post-deviant curvature resampling in the canonical config path.

### `experiment/stimuli_generation_v12.m`
- Add a new v12 generator that keeps the v10 data contract (`xySeqs`, `xySeqsPredicted`, `Cfg`, `repro`) while introducing a random post-deviant curvature mode.
- Implement explicit deviant-curvature precedence (`randomizeCurvatureOnDeviant` over `flipCurvatureOnDeviant`) and keep no-deviant predicted baselines free of deviant-only curvature/direction changes.
- Expand script-level documentation and section comments with usage examples, data-flow notes, assumptions, and output semantics aligned to the experiment/simulation pipeline.

## Stimuli generation documentation

### `experiment/stimuli_generation_versions.md`
- Extend the input-file parameter table with `isCurvFactorRand` so stored stimuli can be compared by both curvature magnitude and curvature-randomization behavior.
- Add the v12 capability profile and update the parameter-effects section to document random post-deviant curvature controls (`randomizeCurvatureOnDeviant`, `deviantCurvatureRange`).

# 05 Feb 2026 16:06

## Experiment configuration and stimuli generation

### `experiment/lib/Config.m`
- Raise `trialsPerCondition` to 1000 and update the curvature defaults (`curvFactor`, `isCurvValenceRand`) for the current pilot assumptions.
- Change the dialog default viewing distance to 1000 mm to match the new setup.

### `experiment/stimuli_generation_v10.m`
- Rename the analysis-only baseline output from “wannabe” to “predicted,” including variable names and the saved `MovDot_SubXX_predicted.mat` file.
- Update inline documentation to match the predicted baseline naming and data flow.

## Stimuli generation documentation

### `experiment/stimuli_generation_versions.md`
- Add a `dotWidth` column to the input-file parameter table and reorder the rows while inserting Sub87 metadata.
- Update the capability notes to reference the predicted baseline export filename and struct (`xySeqsPredicted`).

## Simulation input building and diagnostics

### `simulations/build_movdot_simulation_inputs.m`
- Switch auto-detection from `xySeqsWannabeDev` to `xySeqsPredicted`, including console messaging and error text.
- Treat predicted deviant inputs as the special-case single-condition sources (condition 45) and keep the empty-trial filtering notes consistent.

### `simulations/debug/plot_paths.m`
- Add seconds-aware colorbar support via `SampleRateHz`, plus optional `ParentAxes`/`SavePlot` controls for subplot integration.
- Default the axis limits to `[-5 5 -5 5]`, update examples, and emit formatted second ticks on the colorbar.

### `simulations/debug/plot_paths_wrap.m`
- Add a helper to plot all conditions in a 2x2 grid with time-graded colorbars, configurable sampling, and a single combined export.
- Use `plot_paths` in parent axes mode to keep subplot styling and saving consistent.

## Simulation pipeline scripts

### `simulations/PIPELINE_simulation.m`
- Load center-relative paths for both conditions, plot all four dot/condition paths via `plot_paths_wrap`, and route the overview plot to a single `_all_conditions` output.
- Silence concatenation display output, switch to seconds-based axes/ticks for dRSA matrices and diagonals, and mark the midpoint with `M`.
- Compact output filename tags to only include shuffle/resample when enabled and set the default `dRSAtype` to `PCR`.
- Pass `modelNames` into dRSA border diagnostics and standardize the sampling rate reference via `plotSampleRateHz`.

### `simulations/PIPELINE_simulation_posOnly.m`
- Move to center-relative paths with optional `centerRelativeRectSize`, predicted-devian model inputs, and a shared all-conditions path plot.
- Replace “wannabe” nomenclature with “predicted” throughout the deviant-only model workflow and metadata labels.
- Comment out the output save block while preserving the new compact tag logic and predicted-model repro fields for future re-enable.

## dRSA helper functions

### `simulations/functions/dRSA_border.m`
- Anchor lag-0 to the midpoint sample, compute borders and diagnostics in seconds, and include model labels in console summaries.
- Add shared color limits for autocorr heatmaps and update plots to use seconds-based axes and titles.

## Hypothesis testing scripts

### `simulations/hypothesis_testing/counterfactual_simulation.m`
- Add a counterfactual dRSA pipeline that combines position and direction models (including predicted variants) with reuse prompts for existing results.
- Generate both common and per-subplot colorbar matrix figures, scale text for export, and save consolidated path plots for all conditions.
- Include local helpers for plot layout, figure scaling, and parameter match reporting to keep runs reproducible.

# 03 Feb 2026 11:43

## Experiment configuration and stimuli generation

### `experiment/lib/Config.m`
- Adjust the default `curvFactor` to 0.5 while retaining the prior value as an inline comment.
- Add `flipCurvatureOnDeviant` under deviant settings to centralize the curvature flip toggle.

### `experiment/stimuli_generation_v10.m`
- Read `flipCurvatureOnDeviant` from `Config` instead of a script-local literal.
- Update the header note to point to the Config-driven toggle.

## Stimuli generation documentation

### `experiment/stimuli_generation_versions.md`
- Add an input file parameter summary table with a source priority of `inferred_repro` → `repro` → `Cfg`.
- Populate summary rows for Sub04/05 (missing fields) and Sub88–98 with key parameters including curvature and direction variance.

## Simulation pipeline scripts

### `simulations/PIPELINE_simulation.m`
- Switch the default `participantNumber` to 88 for local runs.

# 03 Feb 2026 10:45

## Simulation pipeline scripts

### `simulations/PIPELINE_simulation.m` (renamed from `simulations/messy_PIPELINE_simulation_experimental.m`)
- Rename the experimental pipeline script to the cleaned `PIPELINE_` name without changing its contents.

### `simulations/PIPELINE_simulation_posOnly.m` (renamed from `simulations/messy_PIPELINE_simulation_posOnly.m`)
- Rename the position-only pipeline script to the cleaned `PIPELINE_` name without changing its contents.

### `simulations/messy_PIPELINE_simulation.m`
- Remove the legacy pipeline script now that the cleaned pipeline naming is in place.

# 03 Feb 2026 10:28

## Repository hygiene

### `.gitignore`
- Ignore editor backup `*~` files to keep temporary artifacts out of version control.

## Simulation input building and experimental pipeline

### `simulations/build_movdot_simulation_inputs.m`
- Remove the `RecenterToFixation` toggle and always emit center-relative paths alongside absolute paths.
- Require a valid rect size for center-relative outputs, using `RectSize` overrides or the stimulus `Cfg.rectSize`.

### `simulations/messy_PIPELINE_simulation_experimental.m`
- Always load center-relative paths and error with guidance if they are missing from inputs.
- Sample a fixed number of paths for overview plots to keep figure output compact.

# 29 Jan 2026 16:13

## Experiment configuration and stimuli generation

### `experiment/lib/Config.m`
- Reduce `trialsPerCondition` to 20 for faster runs while preserving the original value as an inline comment.
- Adjust default curvature settings (`curvFactor` and `isCurvValenceRand`) to match the latest pilot assumptions.

### `experiment/stimuli_generation_v05.m` (renamed from `experiment/stimuli_generation_v5.m`)
- Rename to the zero-padded version label while keeping the v05 generation logic intact.
- Add RNG seeding + state snapshot notes and expand the header documentation to spell out outputs and reproducibility.

### `experiment/stimuli_generation_v06.m` (renamed from `experiment/stimuli_generation_v6.m`)
- Rename to the zero-padded version label while preserving the uniform-start, boundary-rejection logic.
- Expand the header block with reproducibility notes and clarify the output struct fields.

### `experiment/stimuli_generation_v07.m` (renamed from `experiment/stimuli_generation_v7.m`)
- Rename to the zero-padded version label and keep the fit-to-bounds path placement workflow.
- Document the bounded curvature logic and reproducibility metadata in the script header.

### `experiment/stimuli_generation_v08.m` (renamed from `experiment/stimuli_generation_v8.m`)
- Rename to the zero-padded version label and retain the vectorized, min-distance-checked generation path.
- Update header documentation to capture RNG state and output details.

### `experiment/stimuli_generation_v09.m` (renamed from `experiment/stimuli_generation_v9.m`)
- Rename to the zero-padded version label and keep the curvature-flip-on-deviant option.
- Document reproducibility metadata and deviant curvature handling in the script header.

### `experiment/stimuli_generation_v10.m`
- Add a new stimuli generator that preserves v09 behavior while exporting no-deviant baseline paths as `MovDot_SubXX_wannabeDev.mat`.
- Capture RNG state and script settings in the saved repro metadata for analysis reproducibility.

### `experiment/stimuli_generation_versions.md`
- Add a version guide describing which `stimuli_generation_vXX.m` script to use for each path capability.
- Summarize how each version changes bounds handling, curvature logic, and deviant behavior.

### `experiment/CreateInputFiles_v11.m`
- Add a new catch-trial input generator with usage docs, jitter scaling controls, and boundary guards for occlusion paths.
- Centralize catch trial timing, jitter scaling modes, and boundary checks in a script-level configuration block.

## Simulation input building and utilities

### `simulations/build_movdot_simulation_inputs.m`
- Add support for `xySeqsWannabeDev` with auto-detection or `XySeqsField` overrides, plus missing-condition tolerance.
- Introduce `suppressDispText` and empty-trial filtering to keep batch runs quiet and robust.
- Add optional center-relative outputs (`RecenterToFixation`, `RectSize`) and persist the recentering metadata.

### `simulations/compute_avg_path_distance.m`
- Add a helper to compute per-sample dot distances plus mean summaries across trials/time with NaN handling.

## Simulation diagnostics (debug helpers)

### `simulations/debug/plot_direction.m`
- Add a quick diagnostic plot for cosine/sine direction timecourses from a random or specified trial.

### `simulations/debug/plot_direction_cosine_distance.m`
- Add a cosine-distance heatmap diagnostic for direction time series, matching dRSA cosine metrics.

### `simulations/debug/plot_direction_dRSA_by_turn.m`
- Add a diagnostic that splits trials by CW/CCW turn sign and recomputes direction dRSA maps per group.

## dRSA helper functions

### `simulations/functions/dRSA_border.m`
- Add `suppressDispText` and `plotAutocorr` options plus a local diagonal helper for plotting regression borders.
- Gate console output and add optional autocorrelation visualizations to inspect border placement.

### `simulations/functions/dRSA_concatenate.m`
- Add `suppressDispText` parsing and route console output through a quiet flag for batch usage.
- Allow `plotConcat` to be passed as a name/value fallback and propagate quiet mode to reshape warnings.

### `simulations/functions/dRSA_coreFunction.m`
- Add `suppressDispText` handling and use the configured `params.modelDistMeasure` for autocorrelation runs.

### `simulations/functions/dRSA_triggered_subsampling.m`
- Add `suppressDispText` to silence trigger/subsample status prints during batch runs.

## Simulation pipeline scripts

### `simulations/messy_PIPELINE_simulation.m`
- Add quiet-mode wiring, compute dot-to-dot distance summaries, and save path overview figures alongside dRSA outputs.
- Add direction diagnostics (timecourse, cosine distance, angle plots) and align neural distance metrics with model type.
- Skip autocorr border computation for corr-only runs, gating it behind `dRSAtype` and quiet options.

### `simulations/messy_PIPELINE_simulation_experimental.m`
- Add an experimental pipeline that can switch between absolute and center-relative paths with recentering support.
- Integrate dot distance summaries, direction diagnostics, and the same dRSA workflow as the main pipeline.

### `simulations/messy_PIPELINE_simulation_posOnly.m`
- Add a position-only pipeline that supports deviant and wannabe-deviant inputs with aligned metadata handling.
- Preserve the trial-locked subsampling workflow and debug plotting for position-only analyses.

### `simulations/pipeline_recursive.m`
- Add top-k subset tracking with duplicate filtering, checkpointing, and scoring by mean/max absolute correlation.
- Support resuming scans, plotting multiple best subsets, and averaging dRSA matrices across top runs.
- Append a Markdown report section with plots and ranked subset metadata for reproducibility.

# 22 Jan 2026 16:36

## Repository hygiene

### `.gitignore`
- Ignore the local `simulations/pipeline_recursive_backup.m` copy so the backup stays out of version control.

# 22 Jan 2026 16:35

## Simulation pipeline scripts

### `simulations/pipeline_recursive.m`
- Switch the recursive scan to score random subsets by max dot1-vs-dot2 position correlation, then reuse the lowest-correlation subset for a single full dRSA run.
- Store per-iteration subset indices and max-correlation scores in `scanResults`, then capture final outputs in `bestResults` with the selected subset metadata.
- Limit path and dRSA matrix plotting to the selected subset to keep the scan loop non-graphical and faster for batch usage.
- Update script-level documentation to reflect the new subset-selection workflow and output variables.

## dRSA functions

### `simulations/functions/dRSA_concatenate.m`
- Add a `plotConcat` input (default on) to suppress diagnostic figures during batch runs without changing concatenation behavior.
- Validate the new plotting flag and document the new call signatures, including empty-mask usage when only toggling plotting.
- Gate the plotting block behind `plotConcat` and keep the existing figure layout intact when enabled.

# 22 Jan 2026 15:36

## Simulation pipeline scripts

### `simulations/pipeline_recursive.m`
- Add a recursive dRSA pipeline that rebuilds simulation inputs and samples random trial subsets to inspect run-to-run variability without writing outputs.
- Configure per-run plotting for dot paths and dRSA matrices, with position and direction models run separately and labeled by model/time axes.
- Include trial-locked trigger subsampling over full trials, with guardrails for subset size and trial-length alignment, and store per-run results in a workspace struct.
- Embed a local plotting helper that mirrors the debug visualization style while keeping figures in-memory only.

# 16 Jan 2026 22:34

## Simulation pipeline scripts

### `simulations/build_movdot_simulation_inputs.m`
- Resolve stimulus/output paths relative to the repo, load MovDot_SubXX stimuli, and split trials into deviant vs. non-deviant inputs.
- Save dot1/dot2 path arrays with frame-count validation plus metadata (units, dot colors, condition definition) in the output .mat files.

### `simulations/build_simulation_report_assets.m`
- Configure participant inputs, ensure condition-specific dot paths exist, and load cached results or compute dRSA matrices per condition.
- Generate grouped figures and persist report assets/parameters in `simulations/output/subXX` for the report generator.

### `simulations/generate_PIPELINE_simulation_report.m`
- Load precomputed report assets and assemble a PDF with parameter summaries and figure sections via MATLAB Report Generator.
- Validate required assets and write the final report into `simulations/report/subXX`.

### `simulations/messy_PIPELINE_simulation.m`
- Build condition-specific inputs, load dot paths, and run trial-locked triggered subsampling for dRSA computations.
- Plot dot-path and distance diagnostics, then save dRSA matrices/diagonal summaries and PNG outputs.

## Simulation debug plotting

### `simulations/debug/plot_paths.m`
- Render time-graded scatter plots for dot trajectories with optional alpha/axis controls.
- Auto-sanitize the title into a filename and save PNGs into `simulations/debug`.

### `simulations/debug/plot_position_time_distance.m`
- Compute time-by-time mean distance matrices across trials and visualize them as heatmaps.
- Save diagnostics to `simulations/debug` with sanitized titles for repeatability.

## dRSA functions

### `simulations/functions/dRSA_PCR.m`
- Implement PCR-based dRSA with optional autocorrelation and model-regression steps per timepoint.
- Support extra PCA preprocessing and multiple component-selection modes for large models.

### `simulations/functions/dRSA_average.m`
- Average dRSA diagonals across a symmetric time window defined by `AverageTime` and `fs`, handling out-of-range samples as NaNs.

### `simulations/functions/dRSA_border.m`
- Compute autocorrelation regression borders by running dRSA autocorrelations, averaging diagonals, and thresholding by variance.

### `simulations/functions/dRSA_computeRDM.m`
- Build time-resolved RDMs via pdist (or correlation fast-path), with optional normalization for PCR workflows.

### `simulations/functions/dRSA_concatenate.m`
- Concatenate trial/segment data into a features×time matrix from arrays, cell matrices, or .mat paths.
- Generate and append a "Concatenation mask" that flags boundaries, with diagnostic plotting.

### `simulations/functions/dRSA_coreFunction.m`
- Centralize dRSA execution by filling defaults, validating subsamples, and computing neural/model RDMs as needed.
- Dispatch to correlation or PCR engines, including an autocorrelation-only path when requested.

### `simulations/functions/dRSA_corr.m`
- Compute correlation-based dRSA matrices for each requested model against the neural RDM.

### `simulations/functions/dRSA_fastpdist.m`
- Provide a vectorized correlation-distance implementation for fast RDM construction.

### `simulations/functions/dRSA_random_subsampling.m`
- Generate random non-overlapping subsamples from an availability mask with spacing and repetition checks.
- Return indices plus a visualization mask and emit warnings when repetition rates are high.

### `simulations/functions/dRSA_resample_subsamples.m`
- Resample existing subsample sets with replacement to build iteration tensors at a configurable keep percentage.

### `simulations/functions/dRSA_rescaleRDM.m`
- Rescale RDM values to the [0, 1] interval and mean-center each column.

### `simulations/functions/dRSA_standardizeRDM.m`
- Mean-center and z-score RDM values to normalize variance for downstream analysis.

### `simulations/functions/dRSA_subsampling_diagnostics.m`
- Report mask statistics, iteration overlap (Jaccard), coverage, and start-point entropy with summary plots.

### `simulations/functions/dRSA_trialwise_subsampling.m`
- Sample subsamples trial-by-trial, limiting to one window per trial with repetition control and visualization masks.

### `simulations/functions/dRSA_triggered_subsampling.m`
- Time-lock subsamples to trigger points with pre/post windows, spacing constraints, and repetition tracking.

# 16 Jan 2026 22:25

## Experiment documentation

### `experiment/MoveDot1_experiment_vX.m`
- Indent the top banner line and adjust spacing in the output path comment; no runtime behavior changes.

## Trigger analysis tooling

### `inspect_fif_report.py`
- Update replay expectations so each replay-start (151) opens a replay window until the next block's first expected trigger, and tag gaze breaks (150) plus the missing 81 as GazeBreak in the match column.
- Build replay windows for every 151 using expected block order and an end-of-recording fallback so matching pauses per replay block.
- Handle gaze breaks, replay starts, and response triggers (201) explicitly during alignment so they stay visible without advancing expected matching.
- Add a gaze-break count to replay summaries and add HTML match styling for RESPONSE and GazeBreak rows.

## Reports

### `reports/sub98/sub98_report.html`
- Regenerate the report with updated trigger matching so GazeBreak/REPLAY/RESPONSE labels and replay summary counts are included.
- Refresh MNE section IDs and the footer timestamp from the new output.

# 14 Jan 2026 12:50

## Replay trigger semantics and reporting

### `experiment/MoveDot1_experiment_vX.m`
- Update trigger documentation so replay-start (151) is described as firing once per block before the first replay trial, keeping the 100 ms pause note.
- Move the replay-start guard into the per-block loop so each block can emit one 151 before its first replayed trial.
- Comment out the MEG photodiode flicker block in the video playback section.

### `experiment/trigger_codes.md`
- Clarify that trigger 151 is emitted once per block, immediately before the first replayed trial in that block, with the 100 ms pause preserved.

## Trigger analysis tooling

### `expected_triggers.py`
- Fix run-count detection by using the first axis of `TrialStruct` so 2D run-by-trial arrays no longer inflate the run count.

### `inspect_fif_report.py`
- Add replay-rule documentation and a usage line in the script header, aligning expectations with trigger-code definitions.
- Introduce replay validation checks (150-to-81 violations, 151 ordering/counts) and surface them in a new HTML summary panel.
- Mark triggers between the first 151 and the first expected block-2 trigger as REPLAY so expected matching pauses during replay segments.
- Add a block column to the HTML trigger table and color-code REPLAY matches for clarity.
- Centralize missing-row detection and reuse it for event extraction and trial numbering.

## Reports

### `reports/sub98/sub98_report.html` (new)
- New trigger report for subject 98 with replay-aware matching and summary panels.

# 13 Jan 2026 11:57

## Eye-tracker testing modes and gaze-break handling

### `experiment/MoveDot1_experiment_vX.m`
- Extend script-level documentation with gaze-break (150) details, usage steps,
  and testing mode behavior for ignore vs fake tracking, aligned to
  `experiment/trigger_codes.md`.
- Add `Conf.IgnoreEyeTracker`, `Conf.FakeEyeTracker`, and
  `Conf.fakeGazeBreakRate`, then derive `useEyelink`, `useFakeEyeTracker`, and
  `useGazeMonitoring` so real vs simulated gaze is resolved in one place.
- Gate EyeLink setup, recording, and messaging on `useEyelink` to avoid
  hardware I/O during testing while preserving MEG-trigger flow.
- Schedule per-trial fake gaze breaks with a random onset that respects
  `Conf.fixBreakToleranceFrames`, generate out-of-window samples during the
  break, and feed them into the existing fixation-break logic/output.
- Preserve inline trigger-code annotations for gaze-break (150) and replay
  start (151) so the testing path stays consistent with the experiment
  definitions.

# 13 Jan 2026 11:25

## Experiment runtime trigger behavior

### `experiment/MoveDot1_experiment_vX.m`
- Add script-level documentation with usage, inputs/outputs, and trigger logic assumptions aligned to `experiment/trigger_codes.md`.
- Introduce a run-level guard so the replay-start trigger (151) is emitted once before the first replay trial, then pause 100 ms to keep separation from subsequent trial-start pulses.

## Trigger reference updates

### `experiment/trigger_codes.md`
- Clarify that trigger 151 fires once before the first replayed trial begins and includes a 100 ms pause after the pulse.

# 9 Jan 2026 16:23

## Documentation and workflow guidance

### `AGENTS.md`
- Add a requirement to include extensive script-level documentation with usage
  examples, and tighten wording around section-level comments describing data flow.

### `trigger_pipeline.md`
- Update expected/actual trigger output paths to use `derivatives/triggers/subXX`
  and clarify that the actual-trigger CSV is now unified per subject.
- Add multi-run usage examples for expected-trigger generation and report
  creation, including default plotting behavior.
- Document the detailed matching procedure used by `inspect_fif_report.py`,
  covering overlap resolution, expected alignment, missing rows, and extra match
  columns.

## Trigger analysis tooling

### `expected_triggers.py`
- Expand the module docstring with inputs/outputs, run selection behavior, and
  usage examples aligned to `experiment/MoveDot1_experiment_vX.m` and
  `experiment/trigger_codes.md`.
- Make `--run` optional so omission generates all runs found in `TrialStruct`,
  while re-seeding per run to keep outputs identical to standalone invocations.
- Add `resolve_runs` and update CSV output/logging to include run IDs explicitly.

### `inspect_fif_report.py`
- Add multi-run processing when `--subject` is set and `--fif` is default:
  discover all `subXX_runYY.fif` files, require expected block CSVs per run, and
  generate a single report with per-run sections.
- Replace per-run/per-block CSV outputs with one unified detected-trigger CSV
  containing `trial`, `trial_exp`, `run`, and `block` columns populated from
  in-memory alignment results.
- Introduce helpers for run-specific path resolution, CSV parsing, and block
  assignment; external match CSVs now filter by run and skip missing rows to
  avoid false matches.
- Add trial counters based on trial-start triggers (1-80/102) for both detected
  and expected sequences, and surface them in the HTML tables.
- Update default report naming for multi-run output and restrict `--browse-raw`
  to the first run with a clear log message.

### `plot_trigger_chan.py`
- Expand script documentation with inputs/outputs and usage examples, and switch
  the default missing-trigger source to the unified CSV when `--subject` is set.
- Allow `--plot-missing` to accept CSV or report HTML, with subject/run-aware
  path resolution and robust parsing based on table headers.
- Default to absolute time origin when missing markers are plotted and skip
  overlays cleanly when sources are absent.

## Reports

### `reports/sub04/sub04_report.html` (new)
- New multi-run MNE report for subject 04 with per-run raw, sensor layout, and
  trigger-channel sections plus trigger tables that include trial/trial_exp
  columns.
- Footer indicates creation on 2026-01-09 11:33:03 with MNE-Python 1.10.1.

## Python bytecode artifacts

### `__pycache__/inspect_fif_report.cpython-313.pyc`
- Updated Python 3.13 bytecode cache generated by the trigger-report script; no
  textual diff available.

# 8 Jan 2025 15:45

## Changes made to the experiment code and documentation

### `experiment/MoveDot1_experiment_vX.m`
- Delay the `desiredFrameRate` assignment until stimulus FPS is loaded, then set
  `Conf.refrate` from `Dat.Cfg.fps` and mirror it into `desiredFrameRate`.
- Add a safety guard (for non-practice runs) that checks for existing output
  files (`*_SUBxx_RUNyy` and `*_firstBlock`, both with and without `.mat`) and
  aborts the run with a detailed error if any exist, preventing accidental
  overwrites.
- Precompute catch-phase durations in whole frames (fixation and occlusion
  phases) to avoid repeated floating-point comparisons and to keep frame
  boundaries consistent.
- Replace all catch-phase comparisons with the precomputed frame counts,
  including the "last frame" checks and occlusion boundary logic.
- Add a refresh-rate sanity check immediately after opening the Psychtoolbox
  window: compute the measured Hz from `GetFlipInterval` and warn if it differs
  from the expected rate by more than 1 Hz or 5%.
- Update the EyeLink file preamble text to the new project label
  ("Damiano Attention Dot Project").

### `experiment/lib/Config.m`
- Increase `frameFrequency` from 60 to 120 to reflect the higher display
  refresh rate expected by the experiment configuration.

### `experiment/trigger_codes.md`
- Rewrite the trial-onset mapping to make the block/condition logic explicit:
  blockCondition 1 uses 1-20 (condition 0) or 21-40 (condition 45), and
  blockCondition 2 uses 41-60 or 61-80.
- Clarify catch start/end semantics as "glitch" boundaries rather than generic
  start/end wording.
- Expand response trigger notes: 201 is expected; 202 is flagged as unexpected
  and is not in `TriggerValues`.
- Add notes about how `condition` maps to `directionVariance` and how
  `blockCondition` is sourced from `BlockOrder`.

### `experiment/TODO.md`
- Add a new "SAVING ISSUE" checklist with concrete failure modes and mitigation
  ideas (backup save, saving order, catch output indexing risks, teardown
  ordering, and partial save on ESC).
- Add a task to update scripts for the new data folder layout (notably
  `inspect_fif_report.py` for `data/subXX/MEG/subXX_runYY.fif`).
- Add tasks to build a trigger-fix script for subject 99 simulation issues and
  to improve trigger checking for replays (using 150/151 and summary stats).
- Mark two items in the "FROM DAVIDE" section as completed.

## Trigger analysis tooling

### `trigger_events.py` (new)
- Introduce `find_events_with_overlaps`, a helper that wraps
  `mne.find_events(..., output="step")` and adds non-zero-to-non-zero steps so
  overlapping triggers are not lost.
- Implement step merging, mask handling, and duplicate filtering to keep the
  output compatible with MNE event arrays while preserving overlaps.

### `inspect_fif_report.py`
- Switch event detection from `mne.find_events` to
  `find_events_with_overlaps`, enabling detection of overlapping trigger
  transitions.
- Update default FIF selection: when `--subject` is provided and `--fif` is the
  default, the script now targets
  `data/subXX/MEG/subXX_runYY.fif` (with `--run`).
- Expand expected-trigger loading to auto-discover block files via a glob
  pattern, track per-block ranges, and maintain block-order concatenation.
- Add overlap resolution logic that uses the expected sequence to decide which
  overlapping trigger to close or start, and advance the expected index
  accordingly.
- Add automatic CSV output routing: when using the default output path and
  `--subject` is set, write per-block outputs if expected block files exist;
  otherwise write a single per-run output.
- Add an HTML "Missing triggers" summary that counts and reports missing
  expected triggers (e.g., "101 -> X missing out of Y") beneath the table.

### `plot_trigger_chan.py` (new)
- Provide a CLI tool for interactive, scrollable visualization of the trigger
  channel with a slider and mouse-wheel scrolling.
- Support detached plotting (`--no-block`) by spawning a child Python process.
- Add `--plot-missing` to read a report HTML, infer missing-trigger midpoints,
  and overlay them in red on the plot.
- Allow selecting the time origin (`raw` vs `absolute`) and handle clipping of
  missing points to the data range.

## Reports

### `reports/dRSA_Att_test00_report.html` (new)
- New MNE report for the `dRSA_Att_test00` data with a "Triggers" section,
  including an embedded trigger-channel plot and a trigger table.
- Includes a "Check passed: all triggers allowed" status.
- Footer indicates creation on 2026-01-07 16:06:05 with MNE-Python 1.11.0.

### `reports/sub04_run01_report.html` (new)
- New MNE report for subject 04, run 01 with trigger-channel visualization and
  a trigger table that includes a `match` column against expected triggers.
- The table shows multiple 201 entries with match=0 and several 101 entries
  marked "MISSING".
- Includes a missing-trigger summary line (101 missing out of 20) and a
  "Check passed: all triggers allowed" status.
- Footer indicates creation on 2026-01-08 10:55:22 with MNE-Python 1.10.1.

### `reports/sub99_report.html`
- Update the report title to reference `sub99_run01.fif` instead of
  `raw_sub99.fif`.
- Trigger table now includes two match columns: `match` and `match debug1`,
  consistent with comparing against multiple actual-trigger sources.
- Footer timestamp updated to 2026-01-07 16:06:40 (MNE-Python 1.11.0).

## Python bytecode artifacts

### `__pycache__/inspect_fif_report.cpython-313.pyc`
### `__pycache__/trigger_events.cpython-313.pyc`
- Binary `.pyc` files generated by Python 3.13 when running the scripts above.
  These are compiled caches with no textual diff.
# 9 Jan 2026 18:28

## Experiment planning updates

### `experiment/TODO.md`
- Replace the initial "Questions" header with the active checklist and add a
  task to verify non-deviant paths stay within constraints.
- Move the practice-response timing, drift/saccade handling, and false-positive
  feedback items into DONE so the active list reflects remaining work.

## Stimulus generation logic

### `experiment/stimuli_generation_v5.m`
- Add script-level documentation covering purpose, usage, inputs/outputs, and
  assumptions, including the baseline no-deviant path behavior.
- Track a baseline no-deviant direction/position alongside the real path and
  reset it when starting a new path segment.
- Update per-frame angle and position updates to keep baseline angles free of
  deviants for likelihood trials, then compute baseline positions in parallel.
- Expand boundary checks to fail paths if either real or baseline trajectories
  exit bounds, and advance the baseline position after successful frames.

## Commit notes

### `commits_notes.md`
- Add a new entry documenting the staged experiment checklist and
  stimuli-generation changes.

# 12 Jan 2026 19:34

## Catch-trial generation

### `experiment/CreateInputFiles_v10.m`
- Add script-level documentation with usage, inputs/outputs, and catch-jitter
  assumptions for the input/output file workflow.
- Precompute the dot-rect center and max distance to support distance-based
  jitter scaling with an `eps` guard against divide-by-zero.
- Introduce configurable jitter parameters (`Catch.JitterBase`,
  `Catch.JitterSlope`, `Catch.JitterScaling`) to control occlusion catch
  amplitude and scaling mode (`off`, `linear`, `log`, `exp`, with `on` as
  linear).
- Replace the fixed alternating y-offset with per-dot, perpendicular jitter
  driven by `AngleDirection + 90`, alternating sign each frame and scaling by
  distance from the rect center before boundary validation.

## Experiment runtime defaults

### `experiment/MoveDot1_experiment_vX.m`
- Flip the default run-mode flags to a practice configuration: disable MEG,
  eye-tracking, and debug, while enabling practice mode.

### `experiment/MoveDot1_experiment_vX.m~` (new)
- Add a full snapshot of the experiment script (editor backup) mirroring the
  current defaults and logic.

## Stimulus configuration and generation

### `experiment/lib/Config.m`
- Retune curvature settings for 120 Hz by halving `curvFactor` to 1.05 and
  randomizing the curviness valence per path (`isCurvValenceRand = true`).

### `experiment/stimuli_generation_v5.m`
- Expand assumptions to note per-dot curviness valence randomization when
  enabled in `Config`.
- Update per-path curvyness updates to compute signed factors independently for
  each dot, preserving per-dot randomness across paths.

# 13 Jan 2026 10:33

## Fixation monitoring

### `experiment/MoveDot1_experiment_vX.m`
- Reduce the fixation window radius from 10 to 7 degrees by updating
  `Conf.fixWindowDeg`, tightening the eye-tracking break threshold.

## Repository cleanup

### `experiment/MoveDot1_experiment_vX.m~` (deleted)
- Remove the editor backup copy of the experiment script to avoid tracking a
  redundant snapshot.

## Commit notes

### `commits_notes.md`
- Add a new entry covering the fixation window update and backup-file removal.
# 13 Jan 2026 19:45

## MEG video sync gating

### `experiment/MoveDot1_experiment_vX.m`
- Add a per-frame `needsVideoSync` flag that is set only when a trigger pulse is queued, tying `RegWrVideoSync` to actual trigger activity.
- Gate all `Datapixx('RegWrVideoSync')` calls on `Conf.MEG && needsVideoSync` to avoid extra sync writes on frames without scheduled triggers.
- Clarify the trigger scheduling data flow with inline comments so MEG sync behavior stays aligned with trigger emission.
# 21 Jul 2026 11:17

## Expand dRSA PE models and group reporting

### `meg/cluster/submit_drsa_analysis.sh`
- Extend the documented Stage-07 model-source override to include explicit `pe_*` prediction-error variants alongside observed, predicted, and legacy prediction-error sources.
- Clarify that matching occluded-nondeviant trajectories provide both the expected continuation and the baseline used to construct prediction-error models.

### `meg/cluster/submit_drsa_line_summaries.sh`
- Change the Stage-08 interval lag-curve default from no shading to across-subject SEM shading while preserving explicit controls for unshaded and within-subject diagonal-SD displays.
- Update wrapper examples and environment-variable documentation so cluster runs state the new poster-readable group-band default.

### `meg/code/analysis/A_dRSA_alt.m`
- Register six explicit prediction-error variants: 50/50 unconstrained-prior, 50/50 constrained-prior, precision-weighted, Kalman-innovation, latent-curvature, and belief-update models.
- Refactor alternate-model construction into local helpers that preserve the legacy observed, predicted, and raw signed prediction-error behavior while assigning each new variant a distinct source folder, feature label, provenance string, and MAT metadata field.
- Build deterministic PE tensors from actual-minus-expected trajectories using mixture responsibility, leave-one-out constrained priors, timewise precision scaling, online innovation updates, curvature differences, or latent-cause evidence as appropriate.
- Add guarded spatial-scale, motion-area, leave-one-out, curvature, and feature-scale helpers so malformed or sparse configuration values retain conservative finite fallbacks.

### `meg/code/analysis/drsa_stage_utils.m`
- Add human-readable stimulus-model labels for every registered PE source, including feature-specific labels for latent curvature and belief updates.
- Keep legacy prediction-error and neural/model label behavior unchanged for existing outputs.

### `meg/code/preprocessing_basic/07_run_matlab_drsa.py`
- Register the explicit PE model-source names in the Stage-07 CLI and configuration flow, including the existing `pe` alias for legacy `prediction_error`.
- Restrict all PE models to occluded-deviant `A_dRSA_alt` batches and exclude them from simulation batches, preventing unsupported or misleading duplicated outputs.
- Update examples, validation messages, batch-empty errors, and expected-output planning to describe the expanded model-source contract.

### `meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py`
- Make across-subject SEM the default interval lag-curve shading mode and retain `none` and within-subject diagonal-SD as explicit alternatives.
- Include explicit `pe_*` records in report plot selection so new prediction-error models are not filtered from group summaries.
- Clarify that raw Stage-07 matrices may already contain MATLAB neural smoothing and that Python `smoothed_display` is an additional display-only transform.

### `meg/code/preprocessing_basic/09_plot_drsa_sliding_lag_surfaces.py`
- Expand the default Stage-09 target registry from the five core position summaries to include all six explicit PE variants.
- Add the new variants to the occluded-deviant peak-timecourse comparison while preserving stable target ordering and shared report scaling.
- Generalize documentation and selection descriptions so reports adapt to whichever registered PE outputs are present.

### `meg/code/preprocessing_basic/10_drsa_stats.py`
- Document explicit PE variants as supported model and model-source selectors for Stage-10 statistical reports.
- Clarify in report metadata that raw inputs are Stage-07 MAT matrices and may already reflect configured MATLAB neural smoothing.

### `meg/code/preprocessing_basic/drsa_result_utils.py`
- Detect arbitrary `pe_*` model sources directly from result-directory components before applying legacy observed/predicted inference.
- Map each registered PE source to a stable report label so Stage-08, Stage-09, and Stage-10 consumers display self-describing model names.

### `meg/code/preprocessing_basic/pipeline_config.py`
- Accept all registered explicit PE sources in `matlab_drsa.model_source` and `matlab_drsa.model_sources` validation.
- Expand validation errors to distinguish the legacy prediction-error source from registered `pe_*` variants.

### `meg/docs/DRSA_SLIDING_LAG_SURFACES.md`
- Document the six additional occluded-deviant PE targets and their inclusion in the model-comparison peak timecourse.
- Generalize interpretation guidance from a fixed five-target report to a variable target set with shared color scaling.

### `meg/docs/PIPELINE_SPEC_basic.md`
- Extend the Stage-07 alternate-analysis specification with explicit `pe_*` sources, separate output folders, and preserved legacy actual-minus-expected prediction-error semantics.

### `meg/tests/test_cluster_wrappers.py`
- Verify the Stage-08 cluster wrapper defaults to across-subject SEM shading.
- Verify the Stage-07 wrapper exposes an explicit constrained-prior PE source in its documented model-source controls.

### `meg/tests/test_drsa_line_summaries.py`
- Add coverage for the new SEM shading default and corresponding output/report layout.
- Verify `pe_*` result paths produce stable source metadata and self-describing stimulus-model labels.
- Update the end-to-end fake-derivative report expectations for the new default shading directory and summary text.

### `meg/tests/test_drsa_sliding_lag_surfaces.py`
- Expand target-selection fixtures and assertions to cover every explicit PE variant.
- Generate synthetic Stage-07 results for the enlarged target registry and derive expected image/Plotly counts from the configured target and comparison collections.

### `meg/tests/test_matlab_drsa_stage.py`
- Verify expected-output planning includes explicit PE folders only for occluded-deviant alternate analyses.
- Verify MATLAB batch construction excludes PE sources from simulations and nondeviant conditions while passing the selected source into `A_dRSA_alt.m`.

### `meg/tests/test_pipeline_config.py`
- Add a complete MATLAB dRSA configuration fixture demonstrating that legacy and explicit PE model-source selections pass shared configuration validation.

Suggested commit message:
`Expand dRSA PE models and group reporting`
# 21 Jul 2026 11:18

## Add sensor-group dot-side decoding

### `meg/cluster/submit_dot_side_decoding.sh`
- Expand the SGE array from one task per participant to one task per participant x requested decoder sensor group, with the default seven-subject/four-region configuration mapped across 28 tasks.
- Add validated `DAD_DOT_SIDE_SENSOR_GROUPS` and `DAD_DECODER_SENSOR_TYPE` controls, including compatibility aliases for whole-brain and the historical parietal misspelling.
- Pass the selected Stage-07 sensor group and MEG channel type into each participant analysis, nest custom and default outputs by group, and print group-specific completion paths.
- Change the default same-identity pseudo-trial averaging request from one trial to two trials while preserving explicit overrides and the all-trials mode.

### `meg/cluster/submit_dot_side_group_report.sh`
- Add a validated sensor-group list and run the group-report command once per requested Stage-07 region instead of building a single unscoped aggregate.
- Pass the selected group into the Python report, retain validate/rebuild/resume and missing-subject controls per group, and report the new group-scoped output layout.
- Change the default group uncertainty band from 95% CI to SEM for poster-readable summaries while preserving the `DAD_DOT_SIDE_INTERVAL` override.

### `meg/code/sanity_checks/01_dot_side_decoding.py`
- Add `--sensor-group` and `--decoder-sensor-type` controls while preserving the historical posterior-gradiometer decoder as the backward-compatible default.
- Parse the canonical MEGIN/Neuromag lobe channel definitions from the shared Stage-07 MATLAB utility, normalize group aliases, filter requested groups to available magnetometer or gradiometer channels, and preserve MNE channel order.
- Record the complete decoder-selection policy, available/requested/selected channel counts, dropped channels, group name, sensor type, and selected channel list in validation and completed summary metadata.
- Scope non-posterior output directories by sensor group and update figure/report titles so results cannot be mistaken for the historical posterior feature set.
- Validate unsupported group/type combinations and require at least two selected channels before decoding.

### `meg/code/sanity_checks/02_dot_side_group_report.py`
- Add a one-group-per-run CLI contract with normalized group names, group-scoped participant input paths, group figures, combined JSON, and linked/embedded HTML reports while preserving historical posterior paths.
- Validate saved participant summaries against the requested sensor group and enforce consistent decoder group, channel type, and channel count across analyses before aggregation.
- Carry sensor-group provenance through validate-only results, completed report metadata, resume checks, analysis records, figure titles, and report summaries.
- Use analysis-specific colors for group means, uncertainty bands, and significant-cluster markers, and omit individual-subject traces from the group accuracy panel for a cleaner aggregate display.

### `meg/tests/test_dot_side_decoding.py`
- Add a synthetic MNE channel-selection test confirming Stage-07 group membership is filtered to the requested sensor type while retaining channel order and selection metadata.

### `meg/tests/test_dot_side_group_report.py`
- Extend fake participant bundles with sensor-group, decoder-type, and decoder-channel-count metadata.
- Verify Stage-07 participant and group report paths are nested and named by sensor group without changing historical posterior paths.
- Verify aggregation rejects a participant summary whose recorded group does not match the requested report group.

Suggested commit message:
`Add sensor-group dot-side decoding`
# 21 Jul 2026 11:19

## Improve dot ERF cluster execution

### `meg/cluster/submit_dot_erf_analysis.sh`
- Add a fully documented SGE array wrapper that runs the group dot-condition ERF workflow once per requested Stage-07 sensor group, with default tasks for occipital, parietal, frontal, and temporal magnetometers.
- Expose subject selection, sensor type/groups, baseline, time window, analysis sampling rate, topomap times, group uncertainty band, cluster thresholds, permutation count, random seed, cluster-table verbosity, local epoch cache, missing-subject handling, validation, rebuild, and resume controls through environment variables.
- Accept colon-separated multi-value qsub overrides as well as comma-separated shell defaults so SGE environment assignment parsing does not split subject, sensor-group, baseline, or topomap lists incorrectly.
- Validate task-to-group mapping and scalar choices before activating the shared micromamba environment, configure isolated headless plotting/cache directories, verify required deployment files, and print the resolved command and scientific settings before execution.
- Forward one selected group to `05_dot_erf_analysis.py` per array task and report the canonical group ERF and HTML output roots on completion.

### `meg/code/sanity_checks/05_dot_erf_analysis.py`
- Change the default ERF baseline to -0.300..-0.100 seconds and the direct CLI analysis start to -1.0 seconds; the cluster wrapper retains -0.5 seconds for compatibility with existing Stage-04 epoch derivatives.
- Add `--local-epoch-cache-dir` for manually copied same-named Stage-04 FIF files when large reads from SMB or other network volumes are unreliable.
- Keep source-manifest validation and recorded input lineage anchored to the canonical derivative while routing only MNE reads through the optional local cache path.
- Validate that the cache directory and per-subject cached FIF files exist, propagate the chosen MNE read path through scanned subject headers, and record it separately in subject summaries.
- Preload each subject's epoch FIF once before condition slicing to avoid repeated lazy seeks through externally backed files and the associated low-level read failures.

Suggested commit message:
`Improve dot ERF cluster execution`
# 21 Jul 2026 11:20

## Ignore generated paper and poster derivatives

### `.gitignore`
- Ignore the complete `paper/poster/` artifact tree, including generated SVG/PDF exports and the PNG/PPTX files already covered by global binary-artifact rules.
- Ignore top-level `paper/dot_erf_*.svg` exports so alternate ERF render locations do not reappear as untracked files.
- Keep authored Markdown and MATLAB files elsewhere under `paper/` eligible for version control.

Suggested commit message:
`Ignore generated paper and poster derivatives`
