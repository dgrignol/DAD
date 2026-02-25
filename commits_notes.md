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
