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
