# SAMBA Poster Reasoning From Pilot MEG Results

Date: 2026-06-25

Purpose: reason through the current pilot results and decide what is strongest
to show on the SAMBA poster. This is not a finished poster text. It is a
decision document: what the results seem to say, what can be claimed safely,
which figures are worth using, and which follow-up analyses would most improve
the poster.

## Technical Summary

The strongest current story is not "we see a deviant-specific prediction-error
representation." The stronger and safer story is that the pipeline can recover
time-resolved position information for moving dots, and that after occlusion the
MEG signal is more consistent with the actually observed trajectory than with a
sustained hidden non-deviant prediction.

In the post-occlusion interval, the occipital dRSA lag summary shows strong
position-model matches for both occluded non-deviant trials and occluded deviant
trials using the actual deviant path model. The deviant-predicted/non-deviant
continuation model is visibly weaker after reappearance. The dot-side decoding
sanity check points in the same direction: post-reappearance dot-side decoding
is significant for occluded non-deviant and occluded deviant actual labels, but
not for the deviant-predicted labels.

This supports a poster framing around dynamic updating of object-position
representations during occlusion and reappearance. The predictive-coding angle
should be phrased cautiously: a predictive-coding account would motivate looking
for a prediction-error response at deviant reappearance, but the current dRSA
results alone mostly show the representational content after updating. A
group-level ERF/evoked contrast around reappearance is therefore one of the
highest-value missing analyses.

## Sources Inspected

Local experiment and methods context:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/methods.md`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/abstract.md`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/control-center/experiment_description.md`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/control-center/experiment_blueprint.md`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/docs/PIPELINE_SPEC_basic.md`

Result sources:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/results_pilot.pptx`
- `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/07_matlab_drsa_group_report_smoothed.json`
- `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/dot_side_decoding/combined/dot_side_all_analyses_group_report.json`
- `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/dot_side_time_generalization/combined/dot_side_time_generalization_all_analyses_group_report.json`
- Analysis code in `/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/analysis`
- Sanity-check code in `/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/sanity_checks`

Literature anchors:

- de Vries and Wurm, 2023, dynamic RSA for lagged and predictive
  representations of naturalistic dynamic input:
  https://www.nature.com/articles/s41467-023-39355-y
- Hogendoorn and Burkitt, 2018, time-resolved EEG evidence for anticipatory
  object-position coding during predictable motion:
  https://psychologicalsciences.unimelb.edu.au/__data/assets/pdf_file/0004/3142147/2018-Hogendoorn-and-Burkitt.pdf
- Teichmann et al., 2022, MEG decoding of object representations during dynamic
  occlusion:
  https://www.researchgate.net/publication/360212225_The_nature_of_neural_object_representations_during_dynamic_occlusion
- Walsh et al., 2020, review of neurophysiological evidence for predictive
  processing:
  https://pubmed.ncbi.nlm.nih.gov/32147856/
- Jiang, Summerfield, and Egner, 2016, prediction error across object features:
  https://scholars.duke.edu/individual/pub1157610

## Experiment Context That Matters For The Poster

Use the current `paper/methods.md` as the main methods source. The
control-center documents are useful for the theoretical motivation, but they do
not perfectly describe the actual pilot.

Actual pilot design:

- N = 7 MEG participants.
- One moving dot, not the older two-dot attention design.
- Conditions:
  - `always_visible`: dot remains visible.
  - `occluded_nondeviant`: dot is temporarily occluded and continues on the
    expected path.
  - `occluded_deviant`: dot is temporarily occluded and changes direction/path
    while hidden.
- Trial duration: 2.67 s, 120 Hz, 320 frames.
- Deviance and first full occlusion: frame 130, about 1.083 s.
- Full occlusion through frame 190.
- Reappearance starts around frame 191, about 1.592 s.
- Deviant and non-deviant trajectories share the same initial path through the
  deviance frame, then diverge in deviant trials.
- Subjects were not perfectly homogeneous in stimulus/runtime stack:
  sub01-sub03 used earlier pilot stacks and 20 sequence identities; sub04-sub07
  used the later rescueTraject stack and 40 sequence identities. This is a
  poster limitation and a reason to include subject-level or leave-one-out
  robustness if possible.

Analysis context:

- Stage 06 exports identity-averaged neural tensors for dRSA.
- Stage 07 MATLAB dRSA compares time-varying neural RDMs with time-varying
  stimulus model RDMs.
- `A_dRSA_alt` is the most useful current poster result family because it gives
  single-model sensor-group matrices, including occipital position dRSA.
- The displayed group dRSA matrices are "average of smoothed" display figures;
  smoothing is for visualization, not a replacement for the underlying MAT
  results.
- Dot-side decoding and time-generalization are sanity checks. They verify that
  left/right position information is present in the MEG signal, but they do not
  replace the dRSA analysis.

## Poster Framing

### Recommended Title Direction

Best current title style:

> Tracking how the brain updates moving-object position representations during
> occlusion

More predictive-coding-forward but still cautious:

> Predictions, sensory evidence, and updating during occluded motion

I would avoid a title that claims "prediction error coding" unless the ERF or PE
model analysis is completed and supports it.

### Intro Logic

1. Dynamic vision requires the brain to continuously estimate where moving
   objects are and where they will be next.
2. Occlusion separates visible sensory input from predicted continuation: while
   the object is hidden, only an internal estimate can track its expected path.
3. A hidden deviant path tests what happens when the prediction is violated at
   reappearance.
4. Dynamic RSA lets us ask whether MEG activity is more similar to the actual
   trajectory, the expected trajectory, or a prediction-error model across time.

### Main Research Question

> During temporary occlusion of a moving object, does MEG activity reflect the
> expected continuation, the actually observed trajectory after reappearance, or
> a prediction-error response to mismatch?

### Working Hypotheses To Put On The Poster

Prediction/trajectory-maintenance hypothesis:

- During occlusion, position representations should align with the expected
  continuation of the trajectory.
- In deviant trials, if the predicted trajectory is maintained through
  reappearance, MEG should initially match the expected non-deviant continuation.

Sensory-updating hypothesis:

- After reappearance, MEG should rapidly align with the actually observed path.
- Deviant and non-deviant actual-path representations may look similar if both
  are dominated by recovered visible position information.

Predictive-coding/prediction-error hypothesis:

- Deviant reappearance should evoke a mismatch response because the sensory
  input violates the expected continuation.
- That response may appear more clearly in ERFs or response amplitude/timing than
  in dRSA content, because dRSA asks which trajectory information is represented,
  not necessarily whether the response is stronger.

## Main Results

### Result 1: Position Information Is Recoverable In The MEG Signal

The dot-side decoding sanity check shows robust left/right position information.
This is important because it validates the basic signal before interpreting
dRSA.

Group dot-side decoding, balanced undersampled metric:

- Always visible observed labels: maximum group mean accuracy 0.811 at 0.158 s;
  significant cluster 0.050-0.658 s, p = 0.0078.
- Occluded non-deviant observed labels: maximum group mean accuracy 0.840 at
  1.750 s; significant clusters 0.000-0.633 s, p = 0.0078;
  1.008-1.308 s, p = 0.0234; and 1.608-2.150 s, p = 0.0078.
- Occluded deviant observed labels: maximum group mean accuracy 0.821 at
  1.758 s; significant clusters 0.000-0.733 s, p = 0.0078 and
  1.700-2.333 s, p = 0.0078.
- Occluded deviant predicted labels: maximum group mean accuracy 0.820 at
  0.183 s; significant early cluster 0.000-0.733 s, p = 0.0078, but no
  significant post-reappearance cluster.

Interpretation:

- The decoding checks support the same broad conclusion as the dRSA: after
  reappearance, the signal carries actual-position information in both
  occluded conditions.
- The absence of a significant post-reappearance cluster for the deviant
  predicted labels weakens a simple sustained-predicted-path story.

Poster use:

- Use this as a small validation/supporting panel, not the central result.
- If space is tight, include only the occluded non-deviant and occluded deviant
  observed/predicted curves in a compact panel or report the result in text.

Figure paths:

![Always visible dot-side decoding](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_decoding/always_visible/observed/always_visible_observed_dot_side_group_accuracy.png)

![Occluded non-deviant dot-side decoding](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_decoding/occluded_nondeviant/observed/occluded_nondeviant_observed_dot_side_group_accuracy.png)

![Occluded deviant actual dot-side decoding](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_decoding/occluded_deviant/observed/occluded_deviant_observed_dot_side_group_accuracy.png)

![Occluded deviant predicted dot-side decoding](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_decoding/occluded_deviant/predicted/occluded_deviant_predicted_dot_side_group_accuracy.png)

### Result 2: dRSA Shows Strong Actual-Trajectory Position Structure After Reappearance

The most poster-relevant dRSA figures are the occipital `A_dRSA_alt` position
matrices.

Use these as the main result panel:

![Always visible occipital position dRSA](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/group_average/matlab_real_neural/A_dRSA_alt/always_visible/corr/matrices/A_dRSA_alt_always_visible_corr_occipital_dRSA_neural_position_groupAverage_smoothed.png)

![Occluded non-deviant occipital position dRSA](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/group_average/matlab_real_neural/A_dRSA_alt/occluded_nondeviant/corr/matrices/A_dRSA_alt_occluded_nondeviant_corr_occipital_dRSA_neural_position_groupAverage_smoothed.png)

![Occluded deviant actual-path occipital position dRSA](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/group_average/matlab_real_neural/A_dRSA_alt/occluded_deviant/observed/corr/matrices/A_dRSA_alt_occluded_deviant_observed_corr_occipital_dRSA_neural_position_groupAverage_smoothed.png)

![Occluded deviant predicted-path occipital position dRSA](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/group_average/matlab_real_neural/A_dRSA_alt/occluded_deviant/predicted/corr/matrices/A_dRSA_alt_occluded_deviant_predicted_corr_occipital_dRSA_neural_position_groupAverage_smoothed.png)

Interpretation:

- The always-visible condition shows that continuous visible motion produces
  structured position-model matching.
- In both occluded non-deviant and occluded deviant actual-path models, the
  post-reappearance period has a strong vertical band around neural time
  ~1.7-1.9 s.
- The deviant predicted-path matrix does not show a comparably strong
  post-reappearance structure.
- This looks more like rapid representational updating to the visible/actual
  path than like sustained representation of the expected non-deviant path after
  the deviant is revealed.

Critical caveat:

- The matrices are descriptive group averages. The poster can show them as the
  main visual pattern, but the direct statistical contrast between actual
  deviant, actual non-deviant, and predicted deviant models should be treated as
  a high-priority next analysis.

### Result 3: Lag Curves Make The Post-Occlusion Pattern Easier To Explain

For the poster, the lag-curve summary may be more readable than the full dRSA
matrices because it condenses neural time minus model time.

Use the post-occlusion lag curve as the central result:

![Post-occlusion occipital dRSA lag curve](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/line_summaries/interval_lag_curves/group_smoothed-display_a-drsa-alt_occipital_corr_stimulus-model-comparison_drsa-neural-position_matrix_post-occlusion_lag_curve.png)

Supporting lag curves:

![Pre-occlusion occipital dRSA lag curve](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/line_summaries/interval_lag_curves/group_smoothed-display_a-drsa-alt_occipital_corr_stimulus-model-comparison_drsa-neural-position_matrix_pre-occlusion_lag_curve.png)

![During-occlusion occipital dRSA lag curve](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/line_summaries/interval_lag_curves/group_smoothed-display_a-drsa-alt_occipital_corr_stimulus-model-comparison_drsa-neural-position_matrix_during-occlusion_lag_curve.png)

Interpretation:

- Pre-occlusion: all occluded conditions share the same visible trajectory until
  deviance/occlusion, so differences here should not be over-interpreted.
- During occlusion: the curves overlap substantially, consistent with a
  weaker or less separable hidden-position representation.
- Post-occlusion: actual non-deviant and actual deviant paths are both strong,
  while the deviant predicted-path model is clearly lower.

Best poster wording:

> After occlusion, occipital dRSA was strongest for the actually visible path in
> both non-deviant and deviant trials. The expected non-deviant continuation in
> deviant trials showed weaker post-reappearance structure, suggesting rapid
> updating toward sensory evidence rather than sustained representation of the
> violated prediction.

Avoid saying:

> Deviant trials have a stronger prediction-error representation.

The current evidence does not establish that yet.

### Result 4: Time-Generalization Supports Stable Position Information, But Is Secondary

The dot-side time-generalization results show broad above-chance train-test
generalization in the early visible period and useful structure in the occluded
conditions.

Group figures:

![Always visible dot-side time-generalization](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_time_generalization/always_visible/observed/always_visible_observed_dot_side_time_generalization_group_matrix.png)

![Occluded non-deviant dot-side time-generalization](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_time_generalization/occluded_nondeviant/observed/occluded_nondeviant_observed_dot_side_time_generalization_group_matrix.png)

![Occluded deviant actual dot-side time-generalization](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_time_generalization/occluded_deviant/observed/occluded_deviant_observed_dot_side_time_generalization_group_matrix.png)

![Occluded deviant predicted dot-side time-generalization](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/sanity_checks/group/ses-01/dot_side_time_generalization/occluded_deviant/predicted/occluded_deviant_predicted_dot_side_time_generalization_group_matrix.png)

Poster use:

- Keep time-generalization off the main poster unless there is extra space.
- It is useful as a robustness/sanity-check result, but it does not answer the
  prediction-vs-updating question as directly as the dRSA model comparison.

## What To Put On The Poster

### Intro / Background Panel

Suggested text:

> Predictive accounts of perception propose that the brain continuously compares
> incoming sensory input with expected future states. Moving-object occlusion is
> a strong test case: while an object is hidden, its trajectory can only be
> inferred from recent motion; when it reappears, the sensory input can confirm
> or violate that prediction. We used MEG and dynamic RSA to ask whether neural
> activity tracks the actual trajectory, the expected continuation, or their
> mismatch.

### Research Question Panel

Suggested text:

> When a moving object disappears and reappears, does MEG activity reflect
> the predicted path, the actual path, or the prediction error caused by a hidden
> trajectory change?

Show a simple schematic with:

- visible dot path,
- occluder/path-band interval,
- non-deviant continuation,
- deviant hidden turn and reappearance,
- analysis windows: pre-occlusion, during occlusion, post-occlusion.

Use or adapt existing schematic/path figures if helpful:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/figures/baseline_curvature_limit_paths.png`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/figures/deviant_turn_limit_paths.png`

Those figures are more stimulus-generation illustrations than final task
schematics, so a cleaner custom schematic would probably be better for the
poster.

### Methods Panel

Keep it compact:

- N = 7 pilot MEG participants.
- One dot moved along smooth curved paths for 2.67 s.
- Conditions: always visible, occluded non-deviant, occluded deviant.
- Deviance and complete occlusion started at ~1.083 s; reappearance started at
  ~1.592 s.
- dRSA compared MEG RDMs with trajectory models at all neural-time by
  model-time combinations.
- Primary analysis shown: occipital sensor-group dRSA with position models.
- Sanity checks: dot-side decoding and time-generalization.

### Results Panel

Recommended layout:

1. Small validation result:
   - dot-side decoding confirms recoverable position information.
   - One compact line plot or text summary.

2. Main dRSA matrix result:
   - occluded non-deviant actual path,
   - occluded deviant actual path,
   - occluded deviant predicted path.
   - Always-visible can be a reference panel if space allows.

3. Main summary result:
   - post-occlusion lag curve comparing actual non-deviant, actual deviant, and
     deviant predicted-path models.

4. Optional:
   - during-occlusion lag curve if you want to emphasize the weak/ambiguous
     hidden-position representation.

Suggested results headline:

> After reappearance, MEG position structure follows the actual path more than
> the violated expected continuation.

More cautious version:

> Pilot dRSA suggests rapid post-occlusion updating toward sensory evidence.

### Discussion Panel

Suggested text:

> The pilot data show robust position-related MEG structure for moving dots and
> a clear recovery of actual-path information after occlusion. In deviant trials,
> the expected non-deviant continuation does not show a strong post-reappearance
> dRSA pattern, suggesting that the representation is rapidly updated once
> sensory evidence becomes available. This is consistent with dynamic updating
> accounts and with prior work showing that position can be represented during
> occlusion only to a limited extent.

Predictive-coding wording:

> A predictive-coding interpretation remains plausible but requires care. A
> deviant reappearance should produce a prediction-error response, but that
> response may be expressed as ERF amplitude/timing rather than as a stable
> trajectory-content pattern in dRSA. The next critical test is therefore a
> direct deviant-vs-nondeviant ERF contrast around reappearance and a
> subject-level contrast of actual-path, predicted-path, and prediction-error
> dRSA models.

Limitations to state:

- Pilot N = 7.
- Participants used slightly different stimulus/runtime stacks.
- Current poster figures are mostly descriptive group averages.
- Direct actual-vs-predicted and deviant-vs-nondeviant statistics should be
  added before making strong inferential claims.
- ERF/prediction-error evidence has not yet been aggregated at group level.

## Critical Evaluation Against The Literature

Dynamic RSA literature supports the method, but also warns against over-reading
latency. de Vries and Wurm define dRSA as a way to compare neural and model
representations across neural-time and model-time. Negative lags are especially
strong evidence for prediction, but positive lags do not prove purely bottom-up
processing because top-down expectations can sharpen or speed up post-stimulus
representations. For this poster, that means the dRSA lag curves are good for
asking "which trajectory model fits the MEG pattern?" but less definitive for
claiming a specific computational mechanism.

The moving-object prediction literature makes a sustained predictive-position
account plausible. Hogendoorn and Burkitt found anticipatory position coding
when motion was predictable. If the DAD task generated a strong hidden-path
prediction, one might expect the deviant trials to retain some match to the
expected non-deviant continuation during occlusion or immediately after
reappearance. The current post-reappearance predicted-path dRSA and predicted
dot-side decoding do not show that strongly.

The dynamic-occlusion literature makes the current result less surprising.
Teichmann et al. report that object position information can be represented
during occlusion to a limited extent, while richer identity features are not
maintained. The DAD pilot fits a similar pattern: position information is
recoverable and returns after reappearance, but hidden-maintained trajectory
content is not clearly dominant.

Predictive processing literature motivates a PE analysis but does not let us
call the current dRSA result a PE result. In a predictive-coding account,
unexpected sensory input should produce prediction error. However, prediction
error can appear as response magnitude, timing, or distributed model mismatch,
not necessarily as stronger decoding of a PE trajectory model. Therefore, the
poster should say "motivates testing prediction error" rather than "shows
prediction error" unless the ERF/PE analyses are completed.

## Recommended Poster Claim Hierarchy

Strongest claim:

- The pipeline recovers dynamic position information from MEG during moving-dot
  perception.

Strong but still descriptive claim:

- After occlusion, both non-deviant and deviant trials show strong dRSA matches
  to the actual post-reappearance path.

Important negative/contrastive claim:

- In deviant trials, the expected non-deviant continuation is weaker after
  reappearance than the actual deviant path.

Cautious interpretation:

- The pilot pattern favors rapid sensory updating over a sustained explicit
  representation of the violated prediction.

Do not claim yet:

- A distinct neural prediction-error representation.
- A statistically established deviant-vs-nondeviant dRSA difference.
- A definitive predictive-coding result.

## High-Value Next Analyses And Ready Prompts

These are ordered by expected value for the poster.

### 1. Direct Subject-Level dRSA Contrasts For The Main Poster Claim

Why it is high value:

- It would turn the main visual result into a statistical statement.
- It directly tests actual deviant vs actual non-deviant vs predicted
  non-deviant continuation.
- It would give a clean poster-ready summary figure with subject-level points.

Prompt for a new chat:

```text
We are in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD. Follow AGENTS.md. Create and run a subject-level group contrast analysis for the pilot dRSA results on the MORWUR volume.

Context:
- Analysis code is in /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/analysis and /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic.
- Current group figures are in /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa.
- The current poster reasoning document is /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper/samba_poster_results_reasoning.md.
- Use the Stage-07 A_dRSA_alt occipital position MAT outputs for subjects sub-01 through sub-07:
  always_visible observed, occluded_nondeviant observed, occluded_deviant observed, and occluded_deviant predicted.

Task:
1. Inspect the existing MAT structures and group report JSON to identify the correct variables, time axes, model-time axes, and occlusion timing.
2. Write a documented Python or MATLAB script under meg/code/analysis or meg/code/sanity_checks, following existing project style, that extracts subject-level dRSA values from defined intervals:
   - pre-occlusion,
   - during occlusion,
   - post-occlusion,
   and from lag bands centered on zero/predictive/lagged ranges.
3. Run paired subject-level contrasts:
   - occluded_deviant actual vs occluded_deviant predicted in the post-occlusion interval,
   - occluded_deviant actual vs occluded_nondeviant actual in the post-occlusion interval,
   - occluded_nondeviant actual vs occluded_deviant predicted in the post-occlusion interval.
4. Use a small-N appropriate nonparametric sign-flip/permutation or bootstrap approach. Report exact p-values where feasible.
5. Create poster-ready figures:
   - subject-level paired dots with group mean/CI,
   - post-occlusion lag curves with significance/contrast annotation,
   - optional compact matrix difference maps.
6. Save outputs under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/poster_contrasts/ without overwriting existing results unless explicitly requested.
7. Update or create a concise Markdown report in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper summarizing methods, outputs, and interpretation.

Deliver the script path, output figure paths, and a short interpretation suitable for updating the SAMBA poster.
```

### 2. Group ERF Contrast Around Reappearance

Why it is high value:

- This directly addresses the predictive-coding question you raised.
- Predictive coding predicts that unexpected reappearance should generate a
  stronger mismatch/prediction-error response, even if dRSA content rapidly
  updates to the actual path.
- If the ERF differs but dRSA does not, that is an interesting poster story:
  response magnitude signals surprise, representational geometry tracks updated
  position.

Prompt for a new chat:

```text
We are in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD. Follow AGENTS.md. Create and run a group ERF analysis for the pilot MEG occlusion experiment.

Goal:
Test whether occluded_deviant and occluded_nondeviant trials differ in ERF amplitude/topography around dot reappearance, and compare this with the dRSA result that actual-path representations look similar after reappearance.

Context:
- Actual pilot methods are in paper/methods.md.
- Pipeline details are in meg/docs/PIPELINE_SPEC_basic.md.
- Preprocessed epochs and reports are under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives.
- Existing sanity-check code 01_dot_side_decoding.py already builds ERF/lateralization/topography diagnostics; reuse its helper logic where possible instead of rewriting everything.

Task:
1. Find the final cleaned/epoched data used by the dot-side decoding and dRSA exports for sub-01 through sub-07.
2. Write a documented Python script, consistent with the existing sanity_checks style, that loads the relevant epochs and computes evoked responses for:
   - occluded_nondeviant,
   - occluded_deviant,
   time-locked to trial onset but analyzed around reappearance (~1.592 s) and/or transformed to reappearance-relative time.
3. Use sensible baseline choices and explicitly document them. Include at least one baseline that does not contaminate the reappearance response.
4. Compute group-level deviant minus nondeviant ERFs for occipital sensors and whole-head MEG, with:
   - sensor-time cluster permutation if feasible,
   - simpler ROI time-course statistics if cluster testing is too heavy.
5. Create poster-ready figures:
   - occipital ERF time courses with reappearance marked,
   - deviant-minus-nondeviant difference wave,
   - topography/topomap at the strongest post-reappearance window.
6. Save outputs under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/erf_reappearance_contrast/ and /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/erf_reappearance_contrast/.
7. Create a Markdown report in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper explaining whether the ERF supports a prediction-error/mismatch interpretation.

Deliver the script path, output figures, statistics, and a short poster-ready interpretation.
```

### 3. Direction/Curvature dRSA For The Deviant Manipulation

Why it is high value:

- The deviant manipulation is a direction/path change, so a position-only model
  may miss part of the theoretically relevant effect.
- Direction dRSA may better capture whether the brain represents the new
  deviant trajectory parameter or the expected continuation.

Prompt for a new chat:

```text
We are in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD. Follow AGENTS.md. Extend or reuse the dRSA analysis to test direction/curvature-related representations in the pilot MEG data.

Goal:
Determine whether the hidden deviant turn is better captured by direction/curvature models than by position-only dRSA, especially after reappearance.

Context:
- MATLAB analysis scripts are in /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/analysis.
- A_dRSA.m already includes position and direction models in the main real-neural dRSA.
- A_dRSA_alt.m currently gives useful single-model sensor-group outputs for position; inspect whether it can be extended to direction or whether main A_dRSA outputs are sufficient.
- Group reports exist under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01 and /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa.

Task:
1. Inspect existing direction-model outputs for subjects sub-01 through sub-07.
2. If occipital sensor-group direction outputs do not already exist, modify or add a documented script to compute A_dRSA_alt-style direction dRSA for:
   - occluded_nondeviant observed,
   - occluded_deviant observed,
   - occluded_deviant predicted,
   - occluded_deviant prediction_error if supported.
3. Run the analysis on the available subjects.
4. Build group average matrices and lag-curve summaries analogous to the current position results.
5. Quantify post-reappearance actual deviant direction vs predicted direction and actual deviant vs nondeviant direction contrasts.
6. Save outputs under the existing dRSA derivative tree using a clear direction-model folder/name, without overwriting existing position results.
7. Create a short Markdown report in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper with figure paths and poster interpretation.

Deliver script changes, output paths, and whether direction/curvature adds a stronger deviant-specific result than position.
```

### 4. Direct Actual-Vs-Predicted Dot-Side Decoding Contrast

Why it is high value:

- Current decoding reports test each curve against chance separately.
- The poster claim would be stronger if there is a direct subject-level
  contrast showing actual deviant labels outperform predicted labels after
  reappearance.

Prompt for a new chat:

```text
We are in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD. Follow AGENTS.md. Create and run a direct group contrast for dot-side decoding actual vs predicted labels.

Context:
- Existing scripts are in /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/code/sanity_checks:
  01_dot_side_decoding.py, 02_dot_side_group_report.py, 03_dot_side_time_generalization.py, and 04_dot_side_time_generalization_group_report.py.
- Existing group report JSON:
  /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/dot_side_decoding/combined/dot_side_all_analyses_group_report.json
- Current separate results suggest post-reappearance decoding for occluded_deviant observed labels but not for occluded_deviant predicted labels.

Task:
1. Inspect per-subject dot-side decoding score CSVs referenced by the group JSON.
2. Write a documented group contrast script that loads per-subject accuracy time courses and computes paired contrasts:
   - occluded_deviant_observed minus occluded_deviant_predicted,
   - occluded_deviant_observed minus occluded_nondeviant_observed,
   - occluded_nondeviant_observed minus occluded_deviant_predicted.
3. Focus on the post-reappearance interval, but also report the full time-course cluster test.
4. Use small-N appropriate sign-flip cluster permutation testing and clearly mark occlusion/reappearance times.
5. Create poster-ready contrast line plots and a compact summary table.
6. Save outputs under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/dot_side_decoding/contrasts/ and matching analysis outputs if needed.
7. Create a Markdown report in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper with results and interpretation.

Deliver the script path, output figure paths, statistics, and whether the decoding contrast supports the dRSA interpretation.
```

### 5. Prediction-Error Model Comparison

Why it is high value:

- It tests the strongest predictive-coding interpretation directly.
- Current prediction-error dRSA figures exist, but the poster should not claim PE
  without a quantified comparison to actual and predicted models.

Prompt for a new chat:

```text
We are in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD. Follow AGENTS.md. Quantify prediction-error model evidence in the pilot dRSA results.

Goal:
Test whether prediction-error position and/or direction models explain MEG dRSA structure in occluded_deviant trials beyond actual-path and predicted-path models.

Context:
- A_dRSA.m includes prediction-error position/direction models for occluded_deviant.
- A_dRSA_alt.m supports path_type = prediction_error for position sensor-group checks.
- Group dRSA outputs are under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa.
- Current poster reasoning cautions that PE is not yet established.

Task:
1. Inventory existing subject-level PE dRSA outputs for wholebrain and occipital sensors.
2. Extract subject-level post-reappearance PE dRSA values and lag curves.
3. Compare PE model strength against:
   - actual deviant position,
   - expected/predicted position,
   - actual deviant direction if available,
   - expected/predicted direction if available.
4. Use paired permutation/sign-flip statistics and effect sizes.
5. Make poster-ready figures:
   - PE vs actual vs predicted subject-level paired plots,
   - PE lag curve,
   - optional PE matrix.
6. If PE evidence is weak or diffuse, state that explicitly and explain whether this supports rapid updating rather than separable PE coding.
7. Save outputs under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/analysis/group/ses-01/drsa/prediction_error_model_comparison/.
8. Create a Markdown report in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper.

Deliver script path, output paths, statistics, and a poster-safe interpretation of PE evidence.
```

### 6. Robustness And Confound Checks

Why it is high value:

- N = 7 is small and the pilot stacks changed across subjects.
- A leave-one-subject-out and early-vs-late-stack check would protect the poster
  from a result driven by one subject or by the 20-vs-40 identity change.
- Eye/gaze checks would help address whether post-reappearance position
  decoding is contaminated by systematic eye position.

Prompt for a new chat:

```text
We are in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD. Follow AGENTS.md. Create robustness and confound checks for the pilot poster results.

Goal:
Test whether the main post-occlusion dRSA/decoding interpretation is robust across subjects, pilot stacks, and available ocular/gaze information.

Context:
- Pilot methods and subject stack differences are documented in paper/methods.md.
- Main dRSA and decoding group reports are under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01.
- sub01-sub03 used earlier stacks and fewer sequence identities; sub04-sub07 used the later rescueTraject stack.
- sub01 eye-tracker data were lost but EOG was available; sub02-sub07 use saved eye-tracker data.

Task:
1. For the main dRSA post-occlusion metric and dot-side decoding post-reappearance metric, compute:
   - leave-one-subject-out group curves/metrics,
   - sub01-sub03 vs sub04-sub07 descriptive comparisons,
   - subject-level effect-size plots.
2. Inspect available eye-tracking/EOG outputs and determine whether a simple ocular confound check can be performed without rebuilding preprocessing.
3. If feasible, compute condition-wise gaze/ocular summaries around reappearance and test whether deviant/nondeviant differ.
4. Create figures that show whether the main result is robust and whether ocular differences plausibly explain it.
5. Save outputs under /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/ses-01/poster_robustness/.
6. Create a Markdown report in /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/paper with a clear recommendation: safe for poster, needs caveat, or should not be claimed.

Deliver script path, output paths, robustness results, and a short caveat paragraph for the poster.
```

## Bottom Line For The Current Poster

If the poster had to be finalized now, I would center it on this message:

> Pilot MEG+dRSA results show dynamic position coding for moving objects and
> suggest that, after occlusion, neural representations rapidly align with the
> actually reappearing trajectory rather than the violated expected
> continuation. This provides a promising route for testing prediction and
> prediction-error accounts of dynamic vision, but direct PE/ERF and
> subject-level contrast analyses are still needed.

That is interesting, coherent with the existing figures, and not over-claimed.
