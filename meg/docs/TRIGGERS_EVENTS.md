# TRIGGERS_EVENTS.md

## Purpose
This document defines the canonical trigger and event handling contract for the
MEG preprocessing pipeline in this workspace.

It is the reference for:

- trigger extraction behavior,
- transition-trigger handling,
- event table fields,
- quality checks and required-code coverage,
- outputs consumed by later preprocessing/epoching steps.

## Scope
Applies to:

- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py`
- `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/04_epoching.py`
- step companion reports and subject summary reports

## Primary Source References
Behavior is grounded in:

- `experiment/MoveDot1_experiment_occlusion_v18_rescueTraject.m`
- `experiment/trigger_codes_occlusion_v9_rescueTraject.md`
- the selected multisubject profile and its exact runtime CSV pattern,
- retained read-only reported-vs-observed validation artifacts.

## Trigger Code Contract
The parser must support and report the profile-selected condition generation:

- v18/v19: `31, 41, 51` and identities `1..20`,
- v20: `71, 81, 91` and identities `1..40`,
- `111, 112, 114, 115`
- `116, 117, 118`
- `201, 202, 113`
- `150, 151, 152`

### Canonical Meaning Table
| Code | Label | Meaning |
| --- | --- | --- |
| 31 | `trial_condition_onset_always_visible` | Condition onset for always-visible trial. |
| 41 | `trial_condition_onset_occluded_nondeviant` | Condition onset for occluded nondeviant trial. |
| 51 | `trial_condition_onset_occluded_deviant` | Condition onset for occluded deviant trial. |
| 71 | `trial_condition_onset_always_visible` | v20 condition onset for always-visible trial. |
| 81 | `trial_condition_onset_occluded_nondeviant` | v20 condition onset for occluded nondeviant trial. |
| 91 | `trial_condition_onset_occluded_deviant` | v20 condition onset for occluded deviant trial. |
| 1..N | `trial_sequence_identity` | Per-trial sequence ID sent at fixed frame offset. |
| 111 | `occlusion_start` | First frame where occlusion starts. |
| 112 | `occlusion_complete` | First frame where dot is fully occluded. |
| 114 | `occlusion_end_start` | First frame where reappearance starts. |
| 115 | `occlusion_end_complete` | First frame where reappearance is complete. |
| 116 | `catch1_disappear` | Type-1 catch disappearance event. |
| 117 | `catch1_reappear` | Type-1 catch reappearance event. |
| 118 | `catch_question_start` | Catch question start event. |
| 201 | `catch_question_correct` | Catch question ended with correct response. |
| 202 | `catch_question_incorrect` | Catch question ended with incorrect response. |
| 113 | `catch_question_timeout` | Catch question timed out. |
| 150 | `gaze_break` | Fixation/gaze break event. |
| 151 | `replay_start` | Replay start event (after gaze break). |
| 152 | `esc_termination` | ESC-based termination event. |

## Transition Handling Contract
Observed trigger streams can include short digital-line transition pulses that
must be handled explicitly.

### Transition Pulse Pattern
Known hardware pattern family:

- `0 -> X -> Y -> 0`

where either order can occur:

- `0 -> A -> B -> 0`
- `0 -> B -> A -> 0`

and one of the two non-zero codes is the intended event while the other is a
line-update artifact.

### Resolution Rule (Default)
The parser must resolve contiguous transition pairs by keeping the intended
code and dropping the non-intended code.

Primary rule (expected-trigger-first):

1. if expected trigger information is available at that position, keep the code
   that matches expected and drop the other code,
2. apply this symmetrically for both orders (`A->B` and `B->A`),
3. do not rely on pulse width when expected trigger identity resolves the pair.

Fallback rule (only when expected trigger is missing or ambiguous):

1. if exactly one code has width `<= transition_max_samples` (default `1`),
   drop the short code and keep the other code,
2. if both (or neither) satisfy the short-width criterion, apply deterministic
   fallback (keep first, drop second) and mark unresolved warning.

### Reporting Rule
Transition-resolution actions are **not silently discarded**. They must be
reported with at least:

- dropped code,
- kept code,
- expected code (if available),
- decision source (`expected_trigger`, `width_fallback`, `deterministic_fallback`),
- transition order pattern (`X_then_Y`),
- dropped width in samples,
- dropped start/end sample/time,
- source FIF path,
- unresolved flag and confidence score.

### Alignment Status Rule
When reported-vs-observed alignment is performed, transition-resolved cases are
represented as warning matches (`MATCH_WARNING`) rather than hard mismatches
when the kept code matches the expected/reported code, for both orders:

- `0 -> reported -> unexpected -> 0`
- `0 -> unexpected -> reported -> 0`

The built-in segment alignment uses the same reported-trigger-first symmetric
rule, with width and deterministic fallbacks only when reported identity is
unavailable.

## Window Preservation Contract
The parser does not apply a generic time-based debounce. Every reconstructed
non-zero STI101 window is retained except the explicitly documented contiguous
non-zero transition pair, whose kept/dropped evidence is written to the
transition warning table. Raw-window and resolved-event counts remain in the
stage report.

## Required Outputs
Step `02_trigger_extraction_custom.py` must write both machine-readable and
human-readable outputs.

### Machine-Readable Output
A canonical event table (e.g., CSV/JSON/Parquet) with at least:

- `event_index`
- `trigger_code`
- `event_label`
- `sample_start`
- `sample_end`
- `time_start_sec`
- `time_end_sec`
- `width_samples`
- `source_fif`
- `source_fif_path`
- `acquisition_run`
- `raw_sample_frame`
- `behavioral_block`, `behavioral_run`, `behavioral_trial`
- `behavioral_context_source`
- `reported_source_csv`, `reported_row_in_csv`
- `reported_alignment_status`, `reported_alignment_confidence`
- `replay_status`, `is_replay_trial`
- `is_gaze_break`, `is_escape_termination`
- `is_transition_warning` (boolean)
- `transition_partner_code` (nullable)
- `expected_trigger_code` (nullable)
- `transition_decision_source` (nullable string)
- `transition_resolution_confidence` (nullable float)
- `unresolved_transition` (boolean)

### Human-Readable Output
A concise summary with at least:

- total detected events,
- count per trigger code,
- required-code presence/absence table,
- transition warnings summary,
- any mismatch or unresolved conditions.

The stage also generates the legacy deep report bundle by default under:

`02_triggers_extraction_custum_report_deeper/`

The bundle retains the historical subject filenames:

- `subXX_all_blocks_reported_vs_observed.md`,
- `subXX_all_blocks_reported_vs_observed.csv`,
- `subXX_all_blocks_reported_vs_observed_anomalies.csv`,
- `subXX_all_blocks_reported_vs_observed_plots/`,
- `subXX_all_blocks_reported_vs_observed_plots/transition_pulses/`.

Unlike the retired standalone comparator, the renderer consumes the canonical
step-02 alignment and source associations. It does not infer behavioral blocks
from acquisition FIF filenames or perform a second independent alignment.

### Observed-Only Context Rule
An observed STI101 event without a saved CSV row remains `OBSERVED_ONLY`.
Behavioral block/run may be inferred only from exact `MATCH` or
`MATCH_WARNING` neighbors in the same source FIF:

1. matching preceding and following segment labels,
2. an ESC 152 following an exact segment match,
3. an unambiguous leading or trailing event within one FIF.

The inference source is stored in `behavioral_context_source`. It does not
turn the event into a reported match or increase exact-match association
counts.

## Required Validation Checks
At minimum, trigger step validation must report:

1. `stim_channel_present` in all analyzed FIF files,
2. profile-required subject-code coverage (present/missing),
3. dynamic sequence code range sanity (`1..N`),
4. counts of transition-resolution warnings (including collapsed pulses),
5. event table row count before/after post-processing,
6. reported alignment status counts and segment-level source associations.

## Downstream Epoching Link
Epoching condition onsets must come from the selected profile. The current
generation mappings are 31/41/51 for v18/v19 and 71/81/91 for v20.

Any change to these mapping assumptions must be updated here and in:

- params files,
- step companion reports,
- subject summary report entries.
