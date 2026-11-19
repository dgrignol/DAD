# DAD analysis todo

Editable task list for the moving-dot MEG analyses. Each prompt is intended for a new chat in this project, so it only repeats task-specific context.

**Priority:** P0 = immediate; P1 = next; P2 = later.  
**Difficulty:** 1/5 = straightforward; 5/5 = substantial analysis or implementation.  
**Status:** change `TODO` to `IN PROGRESS` or `DONE` as work proceeds.

## 1. [TODO] P0 — Model autocorrelation and cross-correlation · Difficulty 4/5

**Suggested chat setup:** GPT-6.2 Astra · High reasoning

**Certain context**
- The current dRSA plots are produced by `meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py` and `09_plot_drsa_sliding_lag_surfaces.py`.
- `meg/docs/POSTER_ANALYSIS_RUNBOOK_20261117.md` records that actual and predicted trajectories match through frame 130 and first differ at frame 131.

**Inference to verify**
- The useful diagnostic is likely lagged correlation between time-resolved model RDMs, using matched pairwise entries. Confirm this from the model construction before choosing the estimator.

**Prompt for a new chat**
> In this DAD project, compute temporal autocorrelation for each model and lagged cross-correlation between model pairs, then add aligned diagnostic panels to the existing dRSA plots. First inspect how the models/RDMs are constructed and confirm the correlation estimator and shared support. Treat the frame-130 shared prefix and frame-131 divergence as an explicit check. Keep the current dRSA results intact.

## 2. [TODO] P1 — Recheck contrasts after the dependence diagnostics · Difficulty 3/5

**Suggested chat setup:** GPT-6.2 Astra · High reasoning

**Certain context**
- `meg/code/preprocessing_basic/10_drsa_stats.py` supports matrix and line-summary statistics and paired contrasts.
- The runbook contains a submission block with paired contrasts set to `none`; inspect the actual job settings and saved outputs before assuming which contrasts were run.

**Inference to verify**
- Model dependence may change how confidently a contrast can be interpreted as unique model information. Reassess this only after seeing the autocorrelation/cross-correlation results.

**Prompt for a new chat**
> In this DAD project, after the model dependence diagnostics are available, audit and re-evaluate the existing dRSA contrast results. Check the actual job configuration and outputs first, then relate the contrasts to the measured model dependence. Recompute only if needed, preserving the current results and clearly distinguishing shared from unique model information.

## 3. [TODO] P1 — Reduce statistical job runtime · Difficulty 3/5

**Suggested chat setup:** GPT-6.2 Sol · High reasoning

**Certain context**
- In the “Optimise code for faster cluster jobs” chat, Job 11946 was reported as requesting 4 CPUs and 16 GB, with about 394 MB peak memory and roughly 1.15 CPU cores used on average.
- Relevant files: `meg/code/preprocessing_basic/10_drsa_stats.py` and `meg/cluster/submit_drsa_stats.sh`.

**Inference to verify**
- The main gain is likely bounded parallelism across independent tests rather than more memory. Confirm the current bottleneck and resource behavior by profiling before changing the implementation.

**Prompt for a new chat**
> In this DAD project, profile and improve the runtime of the Stage-10 dRSA statistics job. Use the Job 11946 profile from “Optimise code for faster cluster jobs” as context, verify the current MNE and NumPy paths, and benchmark safe parallelization within the allocated CPUs. Keep each group-level test intact, control nested threads, and compare outputs as well as runtime.

## 4. [TODO] P2 — Improve the ERF analysis and discussion · Difficulty 4/5

**Suggested chat setup:** GPT-6.2 Astra · High reasoning

**Certain context**
- The “ERF Univariate MEG Sensor Analysis” chat discussed reappearance-locked deviant-minus-nondeviant ERFs, sensors contralateral to the actual reappearance side, physical stimulus differences, and the danger of assuming a larger signed deviant response.
- Current analysis: `meg/code/sanity_checks/05_dot_erf_analysis.py`; settings and outputs are summarized in `meg/docs/POSTER_ANALYSIS_RUNBOOK_20261117.md`.

**Inference to verify**
- The prior discussion’s proposed posterior ROI and latency windows should be checked against the current sensor types, trial structure, and literature before changing the analysis.

**Prompt for a new chat**
> In this DAD project, review and improve the deviant-versus-nondeviant ERF analysis using the “ERF Univariate MEG Sensor Analysis” discussion and relevant literature. Inspect the current script and outputs, verify sensor selection and latency choices, and account for actual reappearance side and physical stimulus differences. Separate confirmatory choices from exploratory analyses; do not assume the deviant ERF must be larger.

## 5. [TODO] P2 — Test regression controls for model dependence · Difficulty 5/5

**Suggested chat setup:** GPT-6.2 Astra · xhigh reasoning

**Certain context**
- The current poster dRSA analysis uses `corr`; the runbook describes PCR as a different analysis.
- Existing methods are documented in `simulations/documentation/dRSA_autocorr_regression_strategies.md` and `simulations/documentation/dRSA_type_PCR_explained.md`.

**Inference to verify**
- Regression may be useful as a complementary estimate of unique model contribution, but it should not replace the current correlation analysis without validation. Check which documented methods can be applied to the current MEG pipeline.

**Prompt for a new chat**
> In this DAD project, evaluate methods for controlling model autocorrelation and overlap with other models, using the existing PCR/ridge documentation as a starting point. First establish which methods apply to the current MEG dRSA pipeline and what each would estimate. Compare suitable approaches with validation, keep the current `corr` analysis as a reference, and save any new analysis separately.

## Inbox — new ideas to organize

Add new items here as plain bullets. Ask for review to assign priority, difficulty, suggested model/reasoning, and a project-aware prompt.

- [ ] New idea:
