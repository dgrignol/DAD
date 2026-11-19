# MEG Preprocessing Basic Pipeline Specification

Last reviewed: **2026-11-12**.

## Current continuous-in-memory epoch route (2026-11-12)

The user approved accepted 4A/4B annotated inputs -> continuous spectrum-fit
50/100 Hz (10 s) in memory -> unchanged Stage 5 epochs/export, omitting ICA now.
Future selected exact-model ICA applies before that same continuous line removal.
No permanent continuous output is required. `DAD_EPOCH_PREPARATION=inline_no_ica`
(or `inline_ica`) requires a new `DAD_PROCESSING_RUN` label; epochs, exports,
ERF and MATLAB outputs/reports share that isolated namespace. Derived manifests
retain canonical ancestry and actual preparation; input/model/timing/operator
and artifact bindings guard resume. The legacy saved-ICA route remains available.
Accepted Stage 2/3 eligibility/timing, native sampling/auxiliary/annotations and
all established epoch settings remain required. See the executable
[epoching handoff](INLINE_EPOCHING_HANDOFF_20261112.md); no cohort epochs or
application QC are claimed by implementation/testing.

## Parallel cohort execution contract

User-directed default from 2026-11-03: schedule all eligible independent
participants together using disjoint processes or participant array tasks.
Bound active workers by allocated CPUs/memory and concrete shared resources;
limit nested BLAS/OpenMP threads. Keep prerequisite barriers and shared final
publication ordered. No serial cohort loop is the default when participant work
can be isolated. Expose worker controls and per-participant/run progress.

`participant_parallel.py` provides spawned workers, allocation checks, stable
result ordering and one numerical-library thread per child. The ocular importer
joins whole-table read-only preflight and fresh-input phases before parallel
participant saves. Annotation counting and participant HTML rendering join
workers before their single cohort publisher. Current standalone triage/report
wrappers request 36 CPUs/64 GB; --jobs or DAD_PARTICIPANT_WORKERS may reduce
concurrency within the actual allocation. Current scientific checks, source/model
binding, manual precedence and output immutability are unchanged. Existing job
10590 is not restarted/resized by this execution-policy change.

## Scope

This is the technical reference for the external multisubject infrastructure,
preprocessing, and downstream analysis interfaces. It describes available code
and required output contracts. Scientific acceptance and case-specific rules
belong to the
[shared decisions record](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md);
commands, output locations, and review procedure belong to the
[cohort review guide](COHORT_PREPROCESSING_REVIEW_GUIDE.md).

- Code root:
  `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic`
- Read-only raw root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/data`
- Default derivative root:
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives`
- Subjects registered for preflight and trigger extraction: `sub-01` through `sub-39`
- Current analysis cohort: P04-P14 and P16-P39 (35). P01-P03 remain historical
  pilots; P15 is excluded for its participant-specific remaining coverage,
  with all data/provenance retained. No general numeric exclusion rule follows.
- Session: `ses-01`

**Current cohort status:** Stages 1 (sources/digital triggers) and 2 (trial
accounting under the recorded policy) were executed and accepted on 2026-10-19.
The immutable Stage 3 evidence begins with `s3-v001`; the coordinated
`s3-v002` is now published for all 36 participants with the separate revised
eligibility application, P24 donor filter and Task-2 onset method. Original
Stage 2 and `s3-v001` evidence remain unchanged. Stage 3 scientific acceptance
is recorded separately in the 36 hash-bound review records; epoching still waits
for accepted cleaning lineage and production prerequisites. See the [application/API contract](STAGE3_REVIEW_APPLICATION.md) and
[Stage 3 contract](STAGE3_VISUAL_TIMING.md) for exact keys, hashes and handoff.
`p24_acquisition_adapter.py`
remains a deferred prototype outside Stage 1; its old proposals are not accepted
participant corrections.

Stage numbers describe cohort workflow and review gates. Numeric script prefixes
(`02`, `03_A`, `04`, etc.) are legacy compatibility identifiers, not workflow
stage numbers. Stage 1 uses script `02`; Stage 2 prepares trials; Stage 3 is
visual timing; Stage 4 cleaning uses `03_A`-`03_C`; Stage 5 epoching uses `04`.
The runner's `--stages` values are descriptive selectors such as `maxwell_filter`
and `epoching`. Stages 2 and 3 use standalone interfaces rather than runner
selectors. Available code does not establish execution or scientific acceptance.

Keep existing filenames, numbered report paths, `03_tuning_choices.yaml`, receipt
keys such as `stage4_review.03_A_MaxFilt`, FIF provenance tokens, and manifest/hash
bindings unchanged. Historical records and existing messages may still use
"Stage-03" for cleaning or "Stage-04" for epoching; these refer to legacy script
IDs. The workflow names and crosswalk below govern current instructions.

The raw tree is never an output location. Every infrastructure and scientific
stage output passes `pipeline_config.assert_output_path_outside_raw()` before
directory creation or writing.

## Stage interfaces and implementation status

| Workflow stage | Runner selector / standalone interface | Existing entrypoint(s) | Implementation status | Scientific acceptance for the current cohort |
| --- | --- | --- | --- | --- |
| 1. Sources, acquisition headers and digital triggers | `trigger_extraction` plus runner preflight/inventory | `dataset_inventory.py`, `02_trigger_extraction_custom.py` (implementation: `trigger_extraction.py`) | Implemented; executed | Accepted 2026-10-19. Deep trigger reports remain enabled. |
| 2. Trial preparation | Standalone CLI and Stage 2 wrappers | `trial_preparation.py` | Implemented; executed | Accounting under the recorded policy accepted 2026-10-19; original ledger and hashes preserved. |
| 3. Cohort visual timing | Standalone CLI and Stage 3 wrappers/application workflow | `visual_timing.py`, `visual_timing_application.py`, `visual_timing_task2.py` | Implemented; `s3-v002` published | Exact timing/eligibility scope accepted 2026-10-21 for 36 participants through separate hash-bound reviews. |
| 4A. Channels, SQUID review, Maxwell and filter settings | `maxwell_filter` plus SQUID preparation/review helper | `03_A_MaxFilt.py`, `review_squid_jumps.py` | Baseline executed; SQUID extension implemented | Baseline accepted 2026-10-23 for all 185 runs. SQUID implementation approved 2026-10-28; candidate findings/new exclusions remain unaccepted. Final temporal filtering occurs in legacy 4C apply or the approved inline epoch-preparation route. |
| 4B. Artifact intervals | `artifact_annotation` plus source-specific ocular review helpers | `03_B_ArtefactAnnotation.py` | Implemented; current masks operative | Current-mask processing policy accepted 2026-11-09 for 35 participants/180 acquisitions, preserving 18 genuine manual and 162 automatic runs. Detailed interval reinspection is deferred; final ICA/cleaned acceptance remains separate. |
| 4C. ICA and cleaning completion | `ica_manual` review/apply | `03_C_ICA.py` | 35 final models and paired continuous/task review books prepared; P19 numerical repair published | Final component choices/application QC/cleaned-output acceptance pending; model-bound exclusions and accepted 4A/4B remain required. |
| 5. Epoching and preprocessing handoff | `epoching` | `04_epoching.py` | Implemented; saved-ICA and inline continuous preparation; cohort execution pending | Approved no-ICA inline route requires accepted 4A/4B, Stage 2/3 and preservation/Stage 5 QC. ICA-enabled routes retain exact-model Stage 4C guards. |

Infrastructure, deferred code and optional downstream interfaces have no additional
cohort stage number. Analysis-specific exports are outside the preprocessing
endpoint. See the sections below for technical contracts and the shared decisions
record for accepted constraints.

| Interface | Implementation status | Role / scope |
| --- | --- | --- |
| `pipeline_config.py` | Implemented | Merge default, subject, site and optional tuning YAML; validate settings; resolve paths; guard outputs. |
| `run_pipeline.py` | Implemented | Provenance, preflight and subprocess dispatch through export; optional MATLAB dRSA registered. Stages 2 and 3 remain standalone. |
| `p24_acquisition_adapter.py` | Deferred prototype | Outside Stage 1 and the accepted Stage 3 method; old detector defaults/proposals are not accepted corrections. |
| `03_C_ICA_auto.py` | Deprecated experimental | Blocked by default and absent from production dispatch. |
| `06_export_neural_for_drsa.py` (`export_drsa`) | Implemented downstream export | Identity-averaged MAT/JSON export on a fixed MEG channel axis. Separate from Stage 5 acceptance. |
| `07_run_matlab_drsa.py` (`matlab_drsa`) | Optional downstream analysis | Validate completed neural exports and prepare/run the MATLAB dRSA handoff. |
| `XX_ERF_manual.py` | Optional inspection | Interactive ERF inspection using resolved epochs and event definitions. |
| `sanity_checks/01_dot_side_decoding.py` | Optional sanity check | Sliding-window dot-side decoding for observed and expected/predicted trajectory sources. |
| `sanity_checks/03_dot_side_time_generalization.py` | Optional sanity check | Train-time x test-time dot-side decoding for observed and expected/predicted trajectory sources. |

## Stage 2 trial preparation requirements

Prepare a trial-decision ledger from the accepted Stage 1 evidence before
epoch creation. It must preserve source FIFs, original sample frames,
acquisition/block/trial identity, trigger evidence, eligibility, and exclusion
reasons. Carry forward the accepted source-specific exclusions and incomplete
recording boundaries without modifying Stage 1 event/alignment tables.

The exact P30 source/sample keys and exclusion reason are maintained in
[the accepted decision](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md#p30-unreported-incomplete-sequences).
Both excluded candidates must remain visible in the ledger and produce no
epochs. P37's unavailable tail and incomplete final trial, P29's separate partial
recording and block identities, and other partial-session boundaries must also
survive trial preparation.

Implementation is available in `trial_preparation.py`, with dedicated
`submit_preprocessing_stage2.sh` and `submit_preprocessing_stage2_summary.sh`
SLURM wrappers. Execution and scientific acceptance remain separate gates.
See the cohort guide for exact commands and output paths.

Stage 2 acceptance covers ledger accounting, count/identity coverage, selection
effects under the recorded policy, and application of accepted Stage 1 rules.
Those decisions remain accepted unless changed inputs or contradictory evidence
require review. Historical output reports retain pending-review text and their
verified hashes; the decisions record holds scientific acceptance. Any new
accounting inconsistency still requires investigation even with zero
`review_required` rows. Documentation-only confirmations do not require
regenerating those accepted bundles.

Each recorded condition onset owns a stable participant/source/sample key.
Intervals stop at acquisition, behavioral-run, replay-start and warm-up
boundaries. Resolved warning onsets remain candidates; unresolved onsets also
remain accounted for rather than disappearing. Exact accepted exclusions from
subject YAMLs are matched before identity validation. An exclusion with a
missing source/sample or mismatched condition code stops preparation. P37's
accepted final-onset key additionally requires canonical unavailable-tail
evidence in the same source-constrained behavioral/replay context.
P04/P29/P38 terminal-onset keys also require their exact reviewed ESC marker
inside the candidate interval. Those interrupted-trial exclusions therefore
survive changes to the generic ESC policy. P29's unpaired pre-warm-up B6 ESC
has no condition-onset key and creates no trial. There are six accepted
source-specific candidate exclusions across P04/P29/P30/P37/P38.

The trial builder is shared with step 04. Stage 2 uses its non-strict reporting
mode to retain ambiguous identities as review-required, never eligible. Ledger
rows retain event/sample/alignment references, accepted reasons, policy reasons,
and original digital timing. Missing behavioral tail contexts have a separate
table with no invented sample. Acquisition headers are read without loading
signals and provide actual `first_samp` and exclusive end coordinates. Coverage
of the configured future digital epoch window is checked independently of
visual timing and does not establish full stimulus presentation by itself.

Dedicated Stage 2 output bundles refuse existing files by default; explicit
overwrite affects only the enumerated new Stage 2 paths. The completion JSON is
written last and records configuration, decisions-record reference/hash,
implementation, input and output fingerprints, counts and pending acceptance.
The scientific fingerprint excludes orchestration logs/timestamps and added
cleaning lineage while retaining trigger tables, source identities/sizes,
associations, and acquisition headers. Cohort aggregation requires all 36
matching completed bundles, current policies, and intact output/input hashes.

`epoching.require_stage2_ledger: true` is the production default. Step 04 must
load the completed ledger with matching policy/evidence and intact outputs,
reject any remaining review-required candidates, and preserve all exclusions.
Later cleaning must preserve the raw recording sample coordinates/sfreq; step
04 checks these against the ledger. The final epoch ledger updates inclusion
and MNE drop outcomes without modifying the original Stage 2 accounting.
Later visual-timing integration requires an explicit reviewed handoff retaining
the original digital source/sample keys.

## Stage 3 visual timing: requirements and evidence

### Implementation boundary and future handoff

The [roadmap](COHORT_PREPROCESSING_REVIEW_GUIDE.md#next-stage-prepare-cohort-visual-timing)
defines the accepted Stage 3 scope and the remaining Stage 4/5 workflow. The [decided constraints](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md#visual-timing-constraints-and-directions)
govern missing/unusable diode cases, including P11's own later-block donors and
P24's run-5-only donor. They do not select numerical estimates or QC parameters.

Current `04_epoching.py` reads the verified Stage 2 ledger and preserves its
original digital `raw_sample_frame` keys. For P04-P39 production,
`epoching.require_stage3_timing: true` requires an accepted hash-bound timing
review record (`--timing-review` or `epoching.timing_review_path`). The reviewed
overlay supplies a separate `epoch_event_sample_raw` for MNE events. Pilot and
unreviewed proposals cannot feed epoching. `run_pipeline.py` has no registered
cohort timing stage; dedicated standalone wrappers preserve prior evidence.
The standalone P24 adapter writes prototype proposals and must not be treated
as an integrated cohort correction or run during a read-only preflight.
Its preliminary active proposal/audit files were intentionally removed in the
approved Stage 1 cleanup; archived copies are historical. Their absence is not
a missing accepted Stage 1 output.

The candidate [Stage 3 interface](STAGE3_VISUAL_TIMING.md) retains participant,
exact source FIF and original absolute raw onset sample as the stable trial key.
Keep Stage 1 events/alignment, Stage 2 bundles and their hashes intact. It writes separate
versioned outputs with original/proposed samples and times, actual sampling
rate and `first_samp`, signed lag and units, measured/estimated/fallback status,
donor provenance, QC and uncertainty. For a measured edge,
`lag_samples = diode_edge_sample_raw - digital_sample_raw`; proposed timing adds
the accepted signed lag to the original digital anchor. Candidate estimate conversion
uses target-native sampling and nearest-integer rounding with ties away from
zero, retaining residuals; this rule remains subject to Stage 3 review.
Never fabricate an observed pulse or
restore an excluded/unrecorded trial. Recheck affected epoch-window coverage
through the reviewed handoff; do not rewrite the accepted accounting snapshot.

### Published eligibility and exclusion application — `s3-v002`

The immutable `s3-v001` bundles and accepted Stage 2 trial ledger remain the
original evidence snapshot. The revised policy and review decisions are applied
in the separate hash-bound application and published with Task-2 estimates in
`s3-v002`. The supported consumer validates application hashes and mandatory
scientific review scopes. Do not edit Stage 2 eligibility, rewrite an `s3-v001`
bundle, or treat this proposal as formal Stage 3 acceptance.

The feasibility audit joined the saved behavioral `runResults` MAT outputs to
Stage 2 and Stage 3 using exact block, behavioral run, trial, replay, condition,
source FIF and onset sample. The saved outputs preserve `catch_type_code` with
trial identity, while Stage 3 `contexts.json` and the Stage 2 ledger do not store
that classification. Across P04-P39, the MATs contain 2,191 type-2 outputs for
motion conditions 81/91. Five are P37 block 1/run 3 trials 4, 8, 13, 16 and 21
in the logged unavailable tail of `sub37_run01_051026.fif`; they have no recorded
condition onset and no Stage 2 or Stage 3 trial key. The other 2,186 have exact
Stage 2 and Stage 3 keys. Of these, 2,081 were excluded only by the old
`catch_type_2` rule and 105 retain additional independent exclusion reasons.
Their existing maps have proposed windows inside the recordings and no stored
trial-end overlap; 1,819 onset measurements are marked measured, 343 have no
recorded diode event and 24 are marked confirmed unusable. The published application yields 2,079 newly eligible catches and six newly
excluded original trials: 23,963 policy-eligible trials total (21,890 original
+ 2,079 catches - 6 exclusions). Of the newly eligible catches, 1,718 onset
assignments are measured and 361 estimated; all proposed windows are in-bounds
with no trial-end overlap. These remain proposal counts, not scientific
acceptance.

The logged policy makes type-2 catches in occlusion conditions eligible in
principle for motion analysis using the ordinary condition onset, codes 81/91,
subject to all other applicable exclusions and review. Type-1
disappearance/reappearance catches remain excluded. Derive the class only from
saved `catch_type_code`, condition and exact trial identity; Stage 2 marker
fields and code 118 alone are not classification. Keep code 118 question timing
separate from motion onset: retain and review it as lower-priority behavioral
timing evidence, never substitute or pool it with motion-onset lags. Its
behavioral endpoint and clock alignment are still undefined.

| Participant | Exact recorded trial | Review effect | Original Stage 2 status and timing/donor handling |
| --- | --- | --- | --- |
| P24 | Run 5 / block 5 / behavioral run 2 / trial 51, replay=true; code 81 at 582257 in `sub24_run05_100926.fif` | Exclude the whole trial from participant analysis | Exclude all five measured trial events from P24 donor inputs: 81@582257, 111@583362, 112@583496, 114@584004, 115@584138. Recompute their code-specific estimates; expected donor counts are 49 to 48 for code 81 and 90 to 89 separately for codes 111/112/114/115. Preserve all five measured rows. |
| P15 | Run 1 / block 1 / behavioral run 3 / trial 57, replay=true; code 81 at 1074390 in `sub15_run01_030926.fif` | Exclude the whole trial from participant analysis | All timings are measured; no donor estimate applies. |
| P28 | Run 5 / block 5 / behavioral run 2 / trial 13; code 91 at 352984 in `sub28_run05_150926.fif` | Exclude the whole trial from participant analysis | All timings are measured; no donor estimate applies. |
| P29 | Run 1 / block 1 / behavioral run 3 / trial 10, non-replay; code 81 at 1500707 in `sub29_run01_150926.fif` | Exclude the whole trial for its confirmed sequence irregularity, independent of catch type | It was excluded by the original catch rule; it stays out after type-2 inclusion. All timings remain measured evidence; no donor estimate applies. |
| P29 | Same source/run/block; trial 12, non-replay; code 81 at 1519036 | Exclude the whole trial for its confirmed sequence irregularity | Originally eligible. All timings remain measured evidence; no donor estimate applies. |
| P29 | Same source/run/block; trial 13, non-replay; code 81 at 1523905 | Exclude the whole trial for its confirmed sequence irregularity | Originally eligible. All timings remain measured evidence; no donor estimate applies. |
| P33 | Run 4 / block 4 / behavioral run 3 / trial 39, non-replay; code 81 at 811230 in `sub33_run04_220926.fif` | Exclude the whole trial for its confirmed sequence irregularity | Originally eligible. All timings remain measured evidence; no donor estimate applies. |
| P38 | Run 2 / block 2 / behavioral run 2 / trial 24, non-replay; code 81 at 622800 in `sub38_run02_051026.fif` | Exclude the whole trial for its confirmed sequence irregularity, independent of catch type | It was excluded by the original catch rule; it stays out after type-2 inclusion. All timings remain measured evidence; no donor estimate applies. |

P29's decision also retains the 37 currently eligible always-visible trials in
`sub29_run02_150926_halfBlock.fif`, each with its own measured code-71 onset
(18-20 ms; median 19 ms). Its six type-1 catches and interrupted final trial
44 (`71@269558`, then `152@270938`) remain excluded. Preserve this source and
block identity separately from the other acquisitions, including block 6; the
decision adds no averaging or weighting rule. P29's trials 10/12/13 remain
sequence-excluded. Under the old Stage 2 policy, trials 12/13 would reduce its
636 eligible condition onsets to 634; trial 10 was already catch-excluded, so
this historical count does not change for it. Neither 634 nor any old summary is
a final denominator under the revised policy.

The user accepted a future positive estimated onset correction for P04-P07 and
P17: use the median, not the mean, of usable genuinely measured condition-onset
lags from other participants and their donor blocks. Preliminary v24/v25
medians are 19 ms separately for codes 71/81/91. The original `s3-v001`
digital-fallback rows remain unchanged; apply the accepted method only in a
future supported version after freezing and verifying donor keys, exclusions,
weighting, dispersion, native-sample conversion, uncertainty, coverage and
provenance. P04-P07's v20 trigger-before-draw/unscheduled-flip flow differs
from the donor runtimes; transfer bias and participant transfer remain unknown.
P17 used v25 but participant/display transfer remains unknown and its
irreversible 0.1 Hz high-pass limitation remains. The user requested a later
v20 optical bench measurement for P04-P07 that preserves trigger/draw/flip
order and records instrumentation effects. It has not been performed and is
not a retrospective calibration. Neither estimate is a participant-specific
measurement or evidence that +19 ms is closer than zero.

P38's short first recording is handled as an ordinary block with its own source
identity. In `sub38_run01_051026.fif`, 16 of 20 code-71 onset rows are eligible;
all 16 are measured and their proposed windows are inside the recording. The
other three have gaze-break exclusions and one has gaze-break, ESC and reviewed
interruption exclusions. Keep the smaller eligible count visible; apply no
special timing method or weighting. Retain this source separately from block 6.
P38's remaining measured onset assignments are accepted event by event. Its
13 ms code-71 edge is specific to acquisition 5/block 5/behavioral run 1/trial
2 (`71@68884`, diode edge at 68897); retain each other accepted onset at its
own measured edge. Do not use 13 ms as a group correction. P29 likewise retains
measured onset assignments event by event. These case decisions are included in the formal Stage 3 cohort acceptance and
hash-bound review records. P33's remaining 601 original-eligible measured onset
assignments and all 2,079 newly eligible catch assignments were accepted on
2026-10-21; no general one-refresh cutoff was accepted.

### Accepted versioned application and handoff

The accepted `s3-v002` overlay preserves all Stage 2 keys and original
`included_by_trial_policy` values and hashes, as well as all `s3-v001`
measurements, proposals, maps and unavailable-tail evidence. It separately
records `final_analysis_eligible` based on saved MAT `catch_type_code`, exact
trial identity, the restricted type-2 condition-81/91 policy and eight exact
sequence decisions. The cohort summarizer validates 36 matching final bundles.

The final total is 23,963 trials (21,890 original + 2,079 newly eligible
catches - six newly excluded original trials). Five P37 MAT-only tail outputs
remain unrecorded. New catch onsets are 1,718 measured and 361 estimated; all
proposed windows are in bounds with no trial-end overlap. The donor table has
22,120 unique QC-usable measured events and ordinary per-event medians of 19 ms
for codes 71/81/91 (n=7,341/7,422/7,357). P24's five reviewed events are absent
from donor inputs; affected counts are 48 for code 81 and 89 for codes
111/112/114/115, and medians remain 19 ms. Codes 116/117 have six donors each
and no estimate. Cross-participant/runtime and optical patch-to-dot transfer
error remain unknown.

Thirty-six accepted timing JSONs in `s3-v002/reviews/` bind exact manifest and
application hashes, accepted assignment bases and mandatory scientific scopes.
`04_epoching.py` validates these records through `load_reviewed_timing()`,
automatically discovering the adjacent bound `application.json` unless an
explicit `--eligibility-application` path is passed. Epoching still waits for
accepted cleaning lineage. Step 06 and downstream products were not run; Stage
4 cleaning is independent of the accepted trial-timing decision.

### Runtime and evidence distinctions

- **Runtime facts:** P04-P07 used v20. Its condition trigger precedes first-frame
  preparation/drawing and an unscheduled first flip. V21 moves the trigger after
  drawing but still uses an unscheduled first flip. V22 schedules onset, with
  boundary issues investigated in v23/v24; v24 adds the all-trial guard, and v25
  adds dot-offset blanking. Resolve the actual saved runtime per recording:
  P08's accepted block-1 v24 / later v25 distinction must survive this audit.
- **Historical measured evidence:** the timing-test reports below describe
  particular recordings, refresh settings and event classes. They motivate QC;
  their measured lags are not automatically participant calibrations. Trial,
  mid-trial, question and warm-up pulses must not be pooled by assumption.
- **Historical estimate:** Appendix B of the original audit gives approximately
  14.3 +/- 2.7 ms for v20 condition onsets using an assumed extra pre-flip cost.
  This is a counterfactual, not a measured or accepted P04-P07 calibration.
- **Proposals:** the old 50% onset definition, sensitivity thresholds, pooling,
  outlier limits and fallback procedures remain Stage 3 review items. Diode-patch
  light onset is not automatically the moving dot's exact luminance onset.
- **Confirmed versus reported:** use the decisions record for confirmed
  acquisition facts and decided donor directions. P09 run-2 detachment remains
  a session-note lead; P11 blocks 1-2 are now confirmed unusable. Inspect donor
  signal quality without reopening the decided P11 donor restriction.

### Investigation and code catalog

| Reference | Technical use and scope |
| --- | --- |
| [Original audit and Appendix B](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/photodiode_timing_report.md), [detector](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/analyze_photodiode_timing.py), [metadata](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/analysis_metadata.json) | Threshold sensitivity, raw waveform QC, event-class differences, clock interpretation and the explicitly unaccepted v20 counterfactual. Follow recording provenance. |
| [v22 audit](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/v22_scheduled_onset_20260731/photodiode_timing_report.md), [v23 boundary audit](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/260820_sub97/analysis/v23_boundary_guard_timing_report.md), [v24 retest](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/v24_all_trial_guard_resume_retest_20260820/v24_resume_retest_timing_report.md) | Scheduled-onset, run-start, guard/rollover and stop/resume evidence; accompanying analysis metadata identify the test inputs. |
| [Diode-only experiment timeline](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/experiment_timing_photodiode_only/photodiode_experiment_timing_report.md), [analysis](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/tmp/photodiode_timing_audit/analyze_photodiode_experiment_timeline.py) | Distinguish physical within-trial intervals from event-locking shifts. |
| [Legacy v20 source](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/run_outside_lab/MoveDot1_experiment_occlusion_v20_rescueTraject.m), [v21 source](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/MoveDot1_experiment_occlusion_v21_photodiode.m), [v24 source](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/MoveDot1_experiment_occlusion_v24_allTrialGuard.m), [v25 source](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/MoveDot1_experiment_occlusion_v25_dotOffsetBlank.m), [trigger dictionary](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/experiment/trigger_codes_occlusion_v10_rescueTraject.md) | Compare trigger/flip placement and pulse families; validate reference copies against the participant's saved version/profile rather than assuming today's launcher. |
| [Preprocessing code](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic), [deferred P24 adapter](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/p24_acquisition_adapter.py) | Canonical extraction in `02_trigger_extraction_custom.py`; ledger/validator in `trial_preparation.py`; shared trial logic in `downstream_stage_utils.py`; digital-only consumer in `04_epoching.py`. |
| [Default profiles](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml), [cluster site paths](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_cluster.yaml), [participant selectors](/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params) | Resolve exact source/version maps, accepted exceptions and local/cluster roots. |

### Source evidence and representative products

- [Session Log](</Volumes/MORWUR/Projects/DAMIANO/Attention/data/Notes_from_session/Demographics & Sessions log - Session Log.csv>)
  and [Runs Details](</Volumes/MORWUR/Projects/DAMIANO/Attention/data/Notes_from_session/Demographics & Sessions log - Runs Details.csv>):
  contextual experimenter notes. [P24 wrong settings](/Volumes/MORWUR/Projects/DAMIANO/Attention/data/Notes_from_session/wrong_settings_P24_first4blocks.jpeg)
  and [correct settings](/Volumes/MORWUR/Projects/DAMIANO/Attention/data/Notes_from_session/right_settings.jpeg)
  document its acquisition change; verify against saved headers.
- [P04 inventory](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/dataset_inventory.json),
  [P11 inventory](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-11/ses-01/dataset_inventory.json)
  and [P24 inventory](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-24/ses-01/dataset_inventory.json):
  representative channel/filter/runtime context; use equivalent inventories and
  source manifests for every recording. Raw inputs remain under
  `/Volumes/MORWUR/Projects/DAMIANO/Attention/data` and are read-only.
- [P24 canonical events](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/preprocessing/sub-24/ses-01/events/02_trigger_events.csv),
  [P11 ledger](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/preprocessing/sub-11/ses-01/trial_preparation/stage2_trial_decisions.csv),
  [P24 ledger](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/preprocessing/sub-24/ses-01/trial_preparation/stage2_trial_decisions.csv)
  and [cohort summary](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/group/trial_preparation/stage2_cohort_summary.md):
  accepted anchor/accounting examples. These are input references, not Stage 3
  correction outputs. Detailed case rules remain in the decisions record.

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

## Preliminary Annotation-Retention Triage

Approved 2026-10-29 as a prioritization workflow before complete-source interval
review throughout the cohort. It is not an additional numbered scientific gate
and grants no Stage 4 acceptance or production epoching permission.

| Interface | Inputs and behavior | New outputs |
| --- | --- | --- |
| `review_ocular_components.py` | Accepted 4A/current muscle evidence and compatible cached frontal ICA. Walk exact acquisitions, show maps and native sources, save blink/eye selections independently. Zero-based numbers/names, NONE, DEFER and QUIT are explicit. | Atomic model/source-bound choice JSONs and per-participant CSV under `annotation_triage/component_selections/`. These are not production tuning or manual-review receipts. |
| `build_ocular_component_reports.py` | Batch-render existing cached models/sources into ONE offline participant HTML, with internal run selection and navigable traces. Preserve/list saved manual/component work. | `annotation_triage/component_reports/LABEL/sub-XX.html`, bound manifest, note template and shared index. Internal hashed run assets support generation/resume only. |
| `import_ocular_component_choices.py` | Validate a per-run human decision CSV against exact report/model/source identities. Default read-only; --apply persists only compatible new choices after whole-table validation. | Existing component-selection JSON/CSV format; identical choices retained, conflicts rejected, no annotations/acceptance created. |
| `annotation_triage.py` | Preserve validated manual reviews; otherwise use existing detectors/current settings plus automatic muscle/source BADs. Verify original Stage 2 and accepted Stage 3 application/timing. | Immutable hashed masks under `annotation_triage/proposals/`; new report snapshots under `annotation_triage/reports/`. |
| `render_annotation_triage_report.py` | Read saved triage tables only and reconcile every plotted cell with trial outcomes. Embed PNG figures in cohort/participant HTML. Optional baseline comparisons require identical trial identities, eligibility and native epoch coordinates. | New immutable `annotation_triage/visual_reports/LABEL/`: self-contained index/participant HTML, identity bars, condition/run/coverage/reason figures, derived summaries and SHA-256 provenance. |
| `manual_eye_annotations.py --triage-proposals MASK.json` | Validate the exact counted mask/current inputs. Open those intervals without detector reseeding, EDF candidates or SQUID proposals. Refuse any existing manual work; reopen that work normally instead. | Only explicit whole-source SAVE writes ordinary manual annotations/receipts, including triage provenance, automatic-label deletion vetoes and an explicit EDF opt-out. |

The visual renderer does not reopen MEG/ICA or run detectors. After annotation
edits, regenerate the count snapshot first; rendering old tables cannot update
the scientific counts. `submit_annotation_retention_report.sh` supports explicit
full-cohort recount plus render without repeating the completed choice import,
or render-only mode chained after a successful count job. The user submits jobs.
Each participant book has three identity 1-40 barplots on a shared count scale,
run retention, nonadditive rejection causes, exclusive outcomes and exact counts.
Missing inputs remain unranked, true zero survivors remain zero, and disabled
detectors/NONE remain explicit incomplete ocular scoring coverage. No exclusion
threshold, new artifact detector, annotation change or scientific gate is added.

Preserve P04/P05 all-run reviews and P06 runs 1-3. P06 runs 4-5 are included in
automatic processing when unreviewed. Protection is based on current files and
receipts, not a hard-coded list. Orphan/stale manual work is an input error,
never permission to replace it. Reviewed muscle sets replace automatic muscle,
including zero; legacy ocular-only reviews retain existing automatic muscle.
EDF candidates and SQUID proposals are neither loaded nor promoted by triage.
Existing source BADs remain. No SSS, final ICA, temporal policy or sample grid
changes. Compatible detector masks are reusable only with the same source,
annotations, selection/model/cache and detector settings; content hashes detect
corruption. A manual SAVE during detection prevents publishing stale auto work.

Annotation overlap matches MNE's strict overlap on the rounded native epoch
slice, with the inclusive last epoch sample and nonzero first_samp. Exact legacy
BAD_filter_edge remains informational. Keep all overlapping BAD categories and
count each rejected trial once; separate reason counts can overlap. No amplitude
rejection, baseline correction or detrending runs in this metadata-only step.
Every expected condition/identity cell is emitted, including zeros. Zero surviving
epochs is a valid result. Missing/stale inputs are explicitly unranked, while
missing/deferred choices, NONE and unscored source stretches remain visible as
incomplete ocular coverage rather than artifact-free periods.

Reports include full trial metadata, condition counts/percentages, run/condition
counts, all 120 condition/identity cells, overlapping reasons, errors, review
priority and a static count heatmap. Rank by minimum condition count; missing
identity cells break ties. Report snapshots cannot overwrite earlier evidence.
Users can rerun selected participants or the cohort after targeted manual SAVE;
unchanged detector masks are reused. Counts are conditional on the current mask
and do not establish clean signals or a participant-replacement threshold.

The EDF opt-out is source-specific and requires an actual full manual SAVE from
the triage workflow. Later manual reopening retains it; 4B skips that source's
EDF preparation/reconfirmation while retaining all ordinary annotation, lineage
and independent acceptance checks. Historical manual receipts are unchanged.
Production 4A/4B/4C gates, accepted timing and Stage 5 remain as specified below.

**Offline component-report addition, approved 2026-10-30:** one self-contained
participant HTML embeds all exported run data and mag/grad maps, without CDN,
server or Python runtime. Run selector, pan, time slider/jump, zoom and display
gain operate only on embedded arrays. Each default 200-Hz display bin retains
the largest absolute original native sample and its exact within-bin offset.
Any nonfinite acquired sample makes its bin explicitly unscored; final-bin
padding is never treated as missing acquisition. Values use per-component
maximum-absolute native normalization and compact float32/uint16 encoding.
These copies are diagnostic only; native cached traces and detector inputs stay
unchanged. No new fitting, classifier or artifact annotation runs during export.

Versioned participant manifests bind HTML content, viewer/template identity,
source/native evidence, component count and exact cache/model hashes. New
component selections alone do not invalidate a read-only export. Existing
books remain immutable; changed source/viewer/settings require a new label.
Report generation reuses compatible internal assets and refreshes the shared
index under a short lock, including array jobs. Input errors remain errors, not
empty successful trace reports.

The index refreshes itself every 30 seconds; participant books do not refresh.
Index publication prefers atomic replacement. For cross-client shared-storage
PermissionError on replacing an existing file, only the locked mutable
index.html/decisions_template.csv products may be updated in place with flush/
fsync/truncation. Their ownership/ACLs are preserved. This fallback is explicitly
forbidden for immutable books, manifests and assets. A reader may briefly see
an updating index; participant data stay independently immutable. --index-only
rebuilds navigation from completed small manifests without reading source FIFs,
rerendering books, fitting, annotation changes or scientific acceptance.

Human decisions remain per acquisition, never participant-wide component indices.
The importer rejects duplicates, partial choices, unquoted excess CSV columns,
out-of-range indices, wrong source/model/report and existing manual work.
NONE is explicit per group; DEFER applies to both; unfilled rows remain unreviewed.
Complete preflight precedes all writes, with another input check before commit
and ordinary atomic per-source conflict checks. Identical existing decisions are
preserved byte-for-byte. Site-config text may differ across cluster/macOS mounts:
native evidence must match and current effective source/mask/fit settings are
independently validated against the exact model. No stale-model waiver occurs.

## Review And Apply

Stage 4 has three scientific review points, in order: 4A channels/SSS/filtering,
4B muscle/manual ocular intervals, and 4C final ICA/components/applied cleaning.
The user runs cohort jobs. User and Codex inspect independently, reconcile and
record acceptance before downstream work. Execution success is not acceptance.

Keep one subject tuning file at `reports/sub-XX/ses-01/03_tuning_choices.yaml`.
`review_cleaning.py` records explicit reconciled decisions there, binding the
current report evidence, configuration and (for ICA) model contents. Existing
large FIFs use portable filename/size/mtime fingerprints; small calibration,
annotation and model files use SHA-256. Downstream stages reject missing or
changed acceptance. Formal scientific decisions belong in shared storage.

The cluster review wrapper defaults to **maxwell_filter only**. Set
`DAD_STAGES=artifact_annotation` for 4B or `DAD_STAGES=ica_manual` for 4C.
Do not submit the historical combined review/apply chain. Applying ICA selects
only `ica_manual`: it reuses accepted SSS/annotation data and the exact model,
and cannot silently refit. See the cluster runbook for runnable commands.

A provisional 4B pass generates muscle scores and edge-exposure previews even
when eyes are unreviewed. Review each exact source locally with
`manual_eye_annotations.py`, then rebuild 4B explicitly with --overwrite.
A reviewed empty ocular file is valid; a missing receipt/file is unreviewed.
A reviewed empty ICA exclusion list is valid; null/missing is unreviewed.

Review-only processing creates scientific intermediate derivatives and reports.
`--resume` validates their source/settings/content lineage; `--overwrite` is
explicit rebuilding. An upstream change invalidates affected downstream work.
Portable `render_stage03_mne_reports.py` can re-render JSON/PNG evidence without
changing FIFs. Each 4A/4B/4C script calls it automatically for its primary
`_report.html`: figures are embedded, with explanations of the processing,
inspection criteria and limits. Its `--primary-report --overwrite` options
repair older HTML directly from saved evidence without rerunning preprocessing;
the default report-only output remains `_mne_report.html`.
For 4A, `maxwell_qc.py` refreshes an initial individual-channel raw PSD with
sensor-location insets, noisy/flat channel-score heatmaps and one reproducible
random 20-second native raw MEG trace window per exact source.
Red `[D]` marks this run's detected channels; amber `[A]` marks applied-only
channels from harmonization/manual choices. The final matched raw/SSS PSD shows
individual retained channels, excluding applied bads from both sides with
matching y scales per sensor type. PSDs use the same first up-to-120-second
native samples, 1-60 Hz and Welch settings; no diagnostic filter is applied.
The duplicate detection-average PSD is omitted when the new plots are present.
Saved scores and original 4A JSON are reused without changing classifications.
Plot metadata records exact window samples and binds images to the stage report;
changed reports/lists require a plot refresh, not automatic scientific acceptance.
The cluster wrapper's `DAD_MAXWELL_QC_ONLY=1` path invokes this utility directly
without running preprocessing or modifying accepted decisions/data. This plotting
code lives in `meg/code/preprocessing_basic/maxwell_qc.py`; processing orchestration
remains in `03_A_MaxFilt.py`, and HTML composition in `render_stage03_mne_reports.py`.
Stage 4's runner writes `stage4_resolved_config.json` and
`stage4_pipeline_run.json`, preserving Stage 1 inventory/execution snapshots and
the manifest's complete accepted source projection; only derivative lineage is
appended. Do not use ordinary runner dry-run as a read-only preflight: it writes
its execution/configuration snapshots, even though it does not run stages.

## Acquisition-Run Lineage

MEG FIFs are acquisition runs. Behavioral MAT/CSV and EDF files are organized
by experiment block and are never zipped positionally with MEG files.

Every step-03 derivative preserves the complete source stem. Historical pilot
example:

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

<a id="03_a-maxwell-filtering"></a>

## Stage 4A: Maxwell Filtering

Entry point: `03_A_MaxFilt.py`; runner selector: `maxwell_filter`.

Facility inputs now reside at:

```text
/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/meg_calibration_files/fwmegmaxfilter/ct_sparse_triux2.fif
/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/meg_calibration_files/fwmegmaxfilter/sss_cal_3182_250410.dat
```

The workspace deployment copy uses the same `meg/meg_calibration_files` layout;
cluster paths use `/mnt/storage/tier2/morwur`. Relocation preserved contents.
Pilot `meg/data/sub01`, the ZIP and raw `Attention/data` are separate resources.

Detect bad channels independently per acquisition, preserving header bads and
noisy/flat evidence. Save per-window score arrays and per-run raw/SSS PSDs.
Review flagged channels and suitable unflagged controls across the run. The
starting union policy harmonizes the participant's bad-channel basis; explicit
reference/per-run policies and global/per-run add/remove overrides remain.
Manual overrides can defeat harmonization, so final ICA requires agreeing
channel/bad-channel bases rather than silently accepting concatenation warnings.

### SQUID screening and review within Stage 4A

Implementation approved 2026-10-28. Workflow Stage 4A remains one scientific
review gate; `03_A` is a legacy script identifier. Order: channel detection and
acceptance -> SQUID candidates/review -> SSS. The original accepted cohort
participant unions are the screening exclusions; restored good channels remain
eligible. New channel decisions default to union; an explicit run-only exception
must still satisfy the final-ICA channel-basis contract.

`squid_jumps.py` implements short median smoothing (default 9 ms), signed first
differences and absolute robust channel scores. Detector v2 floors the smoothed
derivative's 1.4826 MAD with the raw derivative's 1.4826 MAD: median filtering
can produce a tiny positive scale, not just an exact-zero plateau. This prevents
ordinary high-frequency/quantized noise receiving inflated scores. A deterministic bounded first pass estimates
scale; the second scores every native difference, with chunk halos and no
differencing across EDGE/acquisition gaps. Default score threshold 10, merge gap
0.1 s and +/-0.5 s review padding are adjustable screening parameters, not SQUID
diagnoses or measured recovery durations. Zero-scale/nonfinite channels are
reported as unscored. Baseline changes provide context without requiring
a permanent step. Magnetometers/gradiometers retain independent channel scales.
Within-channel groups and cross-channel peak groups are anchored and bounded
by the merge gap, preventing transitive chaining across an entire recording.
The algorithm identifier is part of proposal identity; v1 pilot outputs are
preserved and are not reused by the corrected detector.

The top-level `squid_jumps` YAML section separates detector/UI settings from
operative Maxwell settings. Preparation uses accepted reports/saved channel
evidence and writes only native candidate JSON/CSV under the preprocessing
`squid_jumps/` directory. Content-addressed files are immutable; channel-basis
snapshots retain configuration/source fingerprints and the accepted evidence.
The mutable channel-basis pointer can support review iteration while rebuilt
SSS awaits acceptance; it grants no downstream scientific acceptance. Fresh
channels need explicit interactive acceptance before SQUID review.

`review_squid_jumps.py --mode review --source EXACT_FILENAME` opens native,
unfiltered traces. Display-only window mean centering and sensor-type scales
computed from that displayed window prevent remote drift/channel offsets from
flattening active traces; within-window steps and source samples remain unchanged.
New review sessions show only the selected candidate's flagged sensors by default;
ALL displays all screened sensors, and VIEW/NEXT/PREV restores event focus. The
title and terminal identify exact flagged channels. Overlapping proposal shading
is shown once, preserving original candidate IDs/evidence and the existing
retained-interval union on SAVE. VIEW retains edits; SAVE records reviewer,
observations, dismissed candidates and retained intervals, including explicit
zero. DISCARD changes nothing. Reopening compatible work never reseeds deleted
candidates. `--reconcile` explicitly reopens changed-context decisions as
unreviewed native-coordinate starting points; new candidates are listed, not
seeded over prior manual decisions. The review supports dismissal, editable
`BAD_squid_jumps`, and channel exclusion through `raw.info['bads']`. BAD_chan is
a channel action, not a time annotation. Channel promotion resolves unchanged
seeded windows only when all contributing sensors were excluded; manually
adjusted/added windows remain deliberate decisions.

`--mode inspect --source EXACT_FILENAME --candidate ID` is a separate read-only
viewer that can run beside an older editable browser with unsaved changes. It
reads exactly one intact current-algorithm/source proposal and bounded native
samples. It never calls preparation, channel-basis persistence or review saving.
Each flagged sensor has independent physical y scales for six-second context,
a 300 ms raw/median close-up and median-filtered first differences. `--context`
adds available colocated sensors clearly labeled as controls, including accepted
bad controls when present. Event-level scores retain their cached definition;
short-window plots do not invent channel-specific scores or classify cause.
DETAIL opens the same plots within new review sessions; `--candidate ID` also
selects the starting candidate for a new editable review. All source samples,
candidate evidence, existing manual annotations and pending review edits remain
unchanged. A muscle-like proposal can be dismissed from SQUID review without
changing Stage 4B muscle decisions. Neither display nor detector score establishes
a SQUID diagnosis or automatic channel exclusion.

All review JSONs are immutable. Only the exact source's `squid_jumps.reviews`
entry in the existing `03_tuning_choices.yaml` is merged using the locked,
latest-map atomic writer. No SQUID operation writes `manual_eye/`, ocular/muscle
receipts or ICA/timing products. Preserve completed P04 all-run and P05 run-1
manual BAD reviews, plus ongoing P05 run-2 work (user status 2026-10-28).

Application validates every exact-source review before SSS. Intervals use
half-open absolute native samples, including first_samp. Only nonempty operative
interval/channel actions enter SSS configuration/signatures; empty reviews and
detector-only tuning do not change SSS identity. Compatible accepted outputs
and reports remain untouched, with a separate SQUID review summary. Changed
actions require explicit affected rebuild/reconciliation: the final-temporal-only
legacy compatibility exception never accepts changed SQUID treatment. Preserve
prior 4A reports/scores/HTML and changed SSS under
`reports/.../history/squid_jumps/<prior-report-sha256>/`; a changed large SSS is
moved before replacement. Resume validates intervals, native/auxiliary samples,
source/settings/resources and final bads; unchanged runs retain their SSS stamps.

SSS explicitly skips `('edge','bad_acq_skip','BAD_squid_jumps')`. MNE 1.10's
plain-SSS buffer check fails when a skipped window splits a full recording;
`maxwell_filter_with_jumps` therefore filters each good stretch independently
with st_duration=None and the common destination/calibration. It restores those
MEG samples onto the original full grid, retaining original annotations and
auxiliary channels. Skipped MEG samples remain acquired/untransformed and BAD;
they must not be interpreted as usable SSS data. No samples are deleted.

Stage 4A owns this category. Stage 4B preserves it from the accepted SSS source
once; copies in old editable eye snapshots are ignored as context, preventing
duplication or replacement. GUI edits to SQUID context in 4B are not saved;
revise treatment through 4A. Severe-period masks, supplementary/final ICA and
final spectrum-fit operations exclude it, and Stage 5 rejects overlapping
epochs while retaining original Stage 2/3 eligibility/timing evidence. Changed
SSS invalidates affected 4B bindings/caches and participant ICA; preserve prior
manual intervals as candidates and explicitly reconfirm them, rather than
deleting or automatically reapproving reviews. Unchanged SSS retains reviews.

No additional final high-pass/low-pass, destination, tSSS, acquisition-limit,
ICA-setting, timing/eligibility or epoch-window change is introduced. Preparation,
tests and implementation approval do not submit jobs or establish acceptance.

References: [FieldTrip detector](https://github.com/fieldtrip/fieldtrip/blob/master/ft_artifact_jump.m),
[FieldTrip artifact tutorial](https://www.fieldtriptoolbox.org/tutorial/preproc/automatic_artifact_rejection/),
[MNE SSS/skip semantics](https://mne.tools/1.10/generated/mne.preprocessing.maxwell_filter.html).

**Retain DAD central_run:** the common destination is the actual acquisition
minimizing summed pairwise RMS displacement of one measured digitization cloud
under the inverse static transforms. It includes translation and rotation.
Reports retain the selected source, all pairwise distances/summed scores and
transform fingerprint. Missing/invalid transforms fail; numerical ties within
1e-6 mm use the lowest acquisition number and are reported.

The previously proposed head-origin alternative was probably a first-pass
oversimplification of DAD's existing central_run approach. It is a different
algorithm and is not used. The selected method is DAD's implementation choice,
not a criterion attributed to FLUX. See the implementation and [MNE destination
API](https://mne.tools/1.10/generated/mne.preprocessing.maxwell_filter.html).

**SSS only:** `maxwell_filter.filter.st_duration_sec: null`; non-null is rejected.
The dataset has per-run HPI, no continuous cHPI. FLUX section 2.4.2 recommends
SSS unless reasons such as excessive movement justify tSSS, discussing continuous
HPI for optimal tSSS operation and its parameter choices. This is a dataset
policy, not a claim that software cannot execute tSSS without cHPI.
[FLUX paper](https://doi.org/10.1016/j.neuroimage.2022.119047).
Time-weighted average head position remains a future lab option requiring a
suitable continuous trajectory; it is not selected here. Duration-weighting a
few static transforms does not establish that missing trajectory.

**Final temporal decision updated 2026-10-28:** no additional high-pass or
low-pass on final continuous MEG. Retain
`notch_filter(freqs=[50,100], method='spectrum_fit', filter_length='10s')` after
ICA. The 10-second window remains subject to evidence review. MISC/STI remain
untouched. P17 and P24 runs 1-4 retain their irreversible 0.1 Hz acquisition
high-pass. The 1-40 Hz filter belongs only to ICA fitting copies.

The intended corrected motion-onset epochs span -0.5 to +2.667 s. The user
expects activity at 0.01 Hz and below to contribute limited within-epoch
variation and prefers retaining low-frequency content over the former long FIR
support. Short epochs do not remove drift; offsets/slopes and trial-order effects
still need QC. This introduces no Stage 5 baseline or detrending decision.

Keep SSS/annotated data broadband for muscle detection and independent ICA fit
preparation. The former 0.01 Hz FIR had 330001 taps at 1000 Hz (165 s
half-support). Its saved exposure reports remain historical information and
never reject otherwise eligible epochs. A new high-pass margin is unnecessary
when no high-pass is applied. Keep incomplete-epoch and reviewed BAD artifact
rules; review spectrum-fit boundaries/short stretches separately. Do not bridge
BAD periods in either line removal or ICA-fit preparation.

**Runtime implemented 2026-10-28:** `highpass_hz: null` omits the additional
high-pass independently of spectrum-fit. Previews have no FIR design/margin,
and `review_cleaning.py --stage 4B` takes no margin argument for this policy.
Spectrum-fit remains a separate operation on MEG only, with BAD/EDGE boundaries
splitting processing; boundary and short-stretch QC remains required.

The SSS configuration projection excludes only named final temporal settings;
all Maxwell/resource/non-temporal settings remain bound. Existing 4A receipts
are validated read-only against the exact known historical temporal projection,
after receipt/report, raw-source and output-stamp checks. Original receipts,
reasons, reports, channel unions and destinations remain unchanged. No 4A rerun
or new acceptance is required. New 4B configuration/evidence binds final policy;
its report records both historical accepted and effective SSS configuration
hashes. On temporal-policy refresh, superseded JSON/embedded HTML previews are
preserved under `history/final_filter/<report-sha256>/` beside the 4B report.
ICA fit identity remains independent; final application additionally binds the
independent-spectrum-fit operation contract, rejecting obsolete output reuse.

**P38 clarification approved 2026-10-28:** spectrum-fit also runs for
`sub38_run01_051026.fif` (137 s). The older full-filter exception was motivated
by the long 0.01 Hz FIR. Preserve its original reason in the existing tuning
`preprocessing.final_filter_skip_sources` map; the historical key now protects
only a future additional high-pass. The exact-source preview retains a warning
to revisit the short recording before such a change. It never bypasses current
line removal, SSS or the separately filtered ICA-fit copies. Scientific 4B/4C
and applied-output acceptance remain pending; see shared PREPROCESSING_DECISIONS.md.

<a id="03_b-physiological-artifact-annotation"></a>

## Current-mask policy adoption (2026-11-07)

The working cohort is P04-P14/P16-P39 (35); P15 remains readable for recovery
and immutable 36-participant provenance, but is excluded from current analysis
and Stage 4/5 wrapper defaults. No general minimum-count threshold is introduced.

`stage4b_adoption.py` supplies an explicit-policy alternative to exhaustive
manual SAVE for the current pass. Its immutable whole-cohort plan validates exact
retention masks and current source/muscle/manual/model/timing inputs; independent
participant processes rebuild standard 4B FIF/CSV/report/lineage, archiving prior
working products. A joined output audit checks native preservation, actual saved
annotation/CSV equality and every trial outcome before authorized policy-specific
processing receipts are recorded. `artifact_annotation.mask_adoption` and per-run
`annotation_adoption` are distinct from unchanged `manual_eye_reviews`.
Saved replacements, reviewed-empty muscle and deletion provenance always win.
Changed inputs fail closed; no EDF/SQUID proposal generation or SSS rerun occurs.
`require_review` revalidates the frozen adoption bindings at downstream entry.
See [commands and evidence contract](STAGE4_CURRENT_MASK_HANDOFF.md).

This explicitly supersedes the exhaustive per-acquisition manual-review gate in
the historical procedures below for this pass. It does not accept final ICA
components, applied cleaning or clean continuous data. Final fitting/settings and
all native/Stage 2/3 contracts remain unchanged. A current-model offline component
book supplements ordinary 4C reports; incomplete properties block component
acceptance. P10 native application QC includes long matched windows and slow
spectra; its retention decision stays open. Stage 5 waits for all final gates.

## Stage 4B: Physiological Artifact Annotation

Entry point: `03_B_ArtefactAnnotation.py`; runner selector: `artifact_annotation`.

Retain FLUX muscle starting settings: mag, 110-140 Hz, threshold 10, minimum good
gap 0.2 s. Inspect full score traces, flagged intervals, nearby unflagged periods
and ordinary controls before choosing per-run thresholds. Score NPZ and PNG files,
combined annotated FIF and CSV remain exact-source products.

Manual review uses available/usable MISC eye analog channels with frontal MEG.
MISC008 is photodiode, not eye data. Preserve original channel names/types; do
not assume all six eye channels work or relabel them EOG. P17 and P24 runs 1-4
have no MISC. The local browser requires a source-specific usability/meaning note,
usable-channel declaration, and explicit SAVE after the entire source is reviewed.
It saves per-source manual-eye -annot.fif files with coordinate/input/hash receipts
in the existing tuning YAML. Missing and reviewed-empty are distinct. Existing
annotations are preserved, including nonzero first_samp and missing meas_date.
Manual BAD_sensor_fault/BAD_movement/BAD_other labels are retained as well; automatic
BAD_muscle remains owned by Stage 4B. EOG helpers are legacy and never executed here.

**EDF-assisted guidance, approved 2026-10-23:** `edf_blink_candidates.py` preserves
original-clock blink endpoints, TRIALID occurrences and recording START/STOP bounds
through the pinned eyelinkio 0.3.0 EDFAPI adapter. It compares accepted canonical
MEG onset/identity/frame evidence within behavioral blocks; filenames alone never
establish correspondence. Exact trial sequences or unique neighbouring contexts
permit candidate placement. Duplicate/missing/mismatched markers, unresolved
transitions, >20 ms local marker-interval disagreement and unknown recording
bounds are flagged, with no automatic repair or ordinal shift. The tolerance is
a guidance QC parameter, not a scientific lag estimate or an epoch criterion.
Blink endpoints map from each trial's own recorded trigger in native MEG samples;
no global offset, fixed ~70 ms correction or Stage 3 motion-onset shift is used.

The helper saves candidate -annot.fif, original-event CSV and mapping JSON in an
exact-source `edf_blink_candidates/` folder beside existing preprocessing outputs.
`EDF_blink_candidate` and `EDF_blink_long_candidate` never mask ICA-fit samples or
reject epochs. >0.5 s is a review flag, not an automatic exclusion; overlapping
binocular bouts merge only within a trial, retaining original records in the CSV.
Raw EDF inputs, canonical events, source associations, settings and SSS coordinates
bind the sidecars and 4B evidence. Guarded resume reuses compatible muscle outputs.

The local raw browser displays candidates with original MISC/frontal channels and
temporary normalized EDF gaze/pupil MISC channels (arbitrary units). Nearest measured samples are
placed on the native grid separately inside valid recording/trial bounds; NaNs,
coverage=0 and informational EDF_no_guidance spans expose gaps. No measured zero
signal is invented, and missing eye values stay NaN even within recorded coverage.
These display channels never enter saved MEG FIFs, ICA or analyses. Optional
--no-edf-traces omits samples. The operator edits/deletes/moves/resizes candidates
and adds missing intervals. Only complete-source SAVE explicitly promotes retained
candidates (including long ones) to BAD_blink and writes the manual-eye receipt.
An unchanged proposal fingerprint prevents deleted candidates from reappearing;
new/changed guidance requires reconfirmation before 4B acceptance. A reviewed zero
is distinct from unavailable EDF guidance or missing review. No cohort processing
or scientific acceptance is implied by implementing this addition.

The 2026-10-26 approved supplementary workflow prepares source-specific caches
on the cluster with `blink_annotation.py` (existing wrapper mode
`DAD_BLINK_AID_ONLY=1`). It reads accepted SSS and current provisional muscle
evidence, selecting ALL good MNE Left-frontal/Right-frontal mag/grad sensors.
Existing isolated 1-40 Hz/200 Hz fit copies feed an ICA capped at ten components
and numerical rank; final 4C n_components=None is unchanged. No fresh manual
severe-period pass is required. Known muscle/severe BADs are excluded, while
ocular examples remain for this supplementary fit only. Unmarked artifacts may
still spoil separation; convergence/fit coverage and manual component inspection
are required. No clear component is not evidence of no ocular artifacts.

The model, native-grid source FIF and metadata live in `blink_aid/`, separate
from final ICA. Source/mask/channels/fit settings bind cache reuse; model hashes,
source stamps and native coordinates are checked. NaN gaps expose excluded,
too-short and filter-edge stretches without generating epoch exclusions. The
per-participant preparation report records failures without accepting cleaning.
Only affected caches need rebuilding; detector settings do not invalidate ICA.

Local `manual_eye_annotations.py --preliminary-ica` LOADS these caches; it never
fits locally. Maps and source traces precede two component/name prompts, blink
first and eye-movement second (`0` = `ICA000`, the first component). Each accepts
comma-separated numbers/names or explicit NONE; invalid/blank answers retry.
The same component may serve both groups. NONE skips only that group's detector,
never establishes artifact absence, and does not disable manual annotation.
Blink sources provide editable auto proposals: per-stretch
detection-only 10 Hz low-pass of the cached 1-40 Hz signals, actual filter guards,
median/MAD absolute scores, SciPy peak detection, merged overlapping intervals.
Starting parameters: five robust standard deviations, 0.25 s peak spacing and
approximate +/-0.25 s windows. Broad peaks (>0.5 s half-prominence width) are
flagged, not discarded. These are review settings, not validated cohort choices
or exact eyelid-duration estimates. Flat/missing/unscored/no-peak results remain
unreviewed. Original MISC identities and EDF gap/coverage rules are unchanged.

Eye proposals use an independent duration method (median-duration-v1) rather
than the blink peak detector. Optional artifact_annotation.eye_movement_auto
settings start with median_window_sec=0.1, threshold_std=4, min_duration_sec=0.1,
merge_gap_sec=0.1 and padding_sec=0.05. Round smoothing to an odd centered native
sample window; trim its half-window support at every finite-stretch boundary.
Per-component median/MAD absolute scores identify continuous excursions meeting
minimum duration before short-gap merging and padding clipped to the same
scored stretch. Union overlapping eye-component proposals without bridging NaN
gaps. Record effective smoothing support, native half-open sample endpoints
including first_samp, component diagnostics and unscored spans. Samples, date,
sampling and recording boundaries are unchanged. Ordinary samples must dominate
robust scaling; flat/unscorable cases require manual review. The cached 1-40 Hz
filter attenuates sustained position plateaus, so this aid cannot recover long
eye position changes or measure complete movement durations. These are starting
review settings, not cohort-validated thresholds or an automatic classifier.

Existing all-frontal-v1 caches and the blink detector/settings/fingerprint inputs
are unchanged. Missing eye selection in an older receipt means NONE in memory,
including P04; no load-time normalization/migration writes occur or imply review.
Keep independent blink_auto and eye_movement_auto provenance, settings, scores,
fingerprints, provisional_intervals_raw and removed_spans_raw. Automatic ocular
categories may overlap; one group never suppresses the other's category.

The combined raw view includes all available MISC, EDF, frontal MEG, selected
normalized AID traces, separate BLINK_score/EYE_MOVEMENT_score channels (one means
that group's threshold; capped at two for display) and annotations. Shared
components are displayed once. Source/AID display normalization uses maximum absolute amplitude
so rare large peaks remain visible; cached/detection amplitudes are unchanged. A score
overview shows each threshold and blink broad-peak flags. `BAD_blink_auto` and
`BAD_eye_movement_auto` live only on this disposable review object until SAVE.
BACK reopens both selections without refitting, retaining edits and unchanged
groups; changed nonempty groups ask separately before regeneration. NONE removes
only that group's untouched session proposals. If regeneration is declined,
selected_components stays distinct from the components/fingerprint actually
used for existing scores/proposals. REDETECT asks BLINK or EYE, then adjusts only
that group's threshold. Only its untouched session proposals are replaced;
manual edits/deletions and previously saved work
are retained. Remembered proposal veto spans never become artifact/epoch rules.
Identical reviewed proposals are not reseeded on reopening. MISC usefulness and
observations are requested after inspection, immediately before persistence.
SAVE records new selections only then; manual-only saves retain prior detector
provenance/vetoes, while missing reviews remain distinct from reviewed-empty.
Older P04 work needs no migration or blanket refit for this addition. Existing
source/muscle/severe-mask changes still invalidate only affected ICA caches.

Default local timing guidance uses the accepted Stage 3 `s3-v002` review (or
explicit --timing-review/configured epoching.timing_review_path), loaded through
the existing read-only Stage 2/3 validators. Accepted corrected motion-onset
samples become event lines; green EPOCH_preview annotations span -0.5 to +2.667 s.
Green guides show accepted eligible timed trials, before later artifact/amplitude rejection.
Recording coordinates must match exactly. Incomplete guides are clipped only
for display and counted; usual incomplete-epoch rules remain. Press a and toggle
EPOCH_preview via Qt Select Visible or Matplotlib visibility checkboxes. Guide
edits are discarded at browser close and their annotations are excluded from
SAVE; neither guides nor existing filter-edge exposure create rejection rules.
--no-timing-marks disables both guides and event lines, without loading timing.
--browser-backend matplotlib provides thin default traces (0.5 points); Qt stays
the default. A scoped color palette supports green on MNE 1.10/1.11 without new
dependencies or private browser mutation. No blink-aid cache rebuild is needed.
Corrected onset width defaults to 2 independently of traces; --onset-linewidth 3
can make it stronger. A small compatibility helper styles existing event artists
and preserves their width across scrolling on both backends; timing is unchanged.
Markers are now fully opaque black by default, above annotation shading. Use
--onset-color blue or --onset-color "#0057B8" for another opaque named/RGB hex
color; labels follow that color. Opacity is fixed at one, without a transparency
option. Both backends retain the style on scrolling; no cluster rerun is needed.
Native code-150 gaze_break markers are orange, above all plot layers and share
--onset-linewidth (default 2). Exact canonical/ledger evidence assigns the owning
trial, including replay identity; ambiguous/missing matches are reported, and
accepted-eligibility contradictions stop for reconciliation. Already excluded
gaze-break trials get separate pale-orange EPOCH_invalid_preview intended windows,
using their hash-bound stored Stage 3 correction, never an invented shift. Missing
usable timing is reported. They do not receive green guides. Both preview types
toggle independently and are removed on close/SAVE; --no-timing-marks hides all
these guides/markers. No BAD interval or new scientific exclusion is introduced.
Every review opening pre-registers BAD_blink, BAD_eye_movement and BAD_segment in
the native annotation selector. They are UI descriptions only, without dummy
intervals; existing labels and marked intervals remain intact. Unused labels do
not exist in the saved annotation FIF and cannot affect fitting/epoch rejection.
They are offered again on reopening, including BACK/REDETECT. Manual BAD_blink
stays distinct from BAD_blink_auto; BAD_segment denotes unclear/non-prototypical
artifacts. The existing explicit-SAVE and reviewed-zero policies remain in force.

SAVE explicitly reviews all retained intervals and stores model/component,
detector settings/peaks and edits/proposal veto provenance in the existing manual
receipt. Retained auto labels then follow usual BAD rejection rules; no batch
output consumes them before reviewed SAVE and a 4B rebuild. DISCARDALL leaves
existing files/receipts unchanged; blanks/typos retry. Explicit
`--overwrite --recover-annotations` still recovers a failed receipt write as
unreviewed work without reseeding EDF candidates. NumPy metadata are normalized
before YAML serialization. Independent review/reconciliation precedes 4B
acceptance; no final ICA/removal or Stage 5 work is performed by this aid.

The local view also performs muscle QC. SAVE stores every BAD-prefixed category,
including BAD_muscle, in the existing manual file. A source/settings-bound
muscle_review receipt distinguishes the reviewed full muscle replacement from
legacy ocular-only reviews, including reviewed-empty muscle sets. Reopening and
4B rebuilding replace the muscle category with this snapshot; they never union
deleted automatic detections back in. Legacy reviews retain automatic muscle
annotations until explicitly reviewed under the new workflow. Other source
annotations and all native/auxiliary preservation checks remain protected.
Automatic scores/counts remain diagnostics and do not use reviewed muscle labels
as their own normalization mask. Recheck rebuilt diagnostics before acceptance.
Changed muscle detection settings invalidate replacement reuse; changed effective
fit masks invalidate only the affected preliminary ICA caches. No raw/source
FIF, trial eligibility, filter-edge rejection policy or final ICA decision changes.

EDF is optional supplementary visual context only. Its recording is discontinuous;
no reliable global continuous EDF-to-MEG synchronization is established. Trigger
matches can support exact trial-segment context. The apparent recording start
about 1.17 s before code 71 is provisional and the endpoint is unknown. Display
unrecorded intervals as gaps, retain segment identities and expose uncertain
alignment; never display concatenated trials as a continuous measured signal.
No EDF synchronization study or new EDF conversion dependency is included.

<a id="03_c-manual-ica"></a>

## Stage 4C: Manual ICA and Cleaning Completion

Entry point: `03_C_ICA.py`; runner selector: `ica_manual`.

Only accepted 4B annotations can feed final fitting. Routine ocular intervals
remain excluded, along with **every BAD-prefixed category**. Filtering uses
`skip_by_annotation=('edge','bad')`, not MNE's narrower default. Prepare each good
stretch separately: native-rate MEG 1-40 Hz FIR; remove full half-support at each
edge (normally 1.65 s); polyphase resample that stretch to 200 Hz; trim a further
0.1 s at each resampling edge. Require at least 1 s of retained data. These are
fit-only losses, never edits to trial eligibility. Concatenate only afterward,
with transform/channel mismatch errors. Source absolute sample ranges and actual
fit sample counts/durations are saved and reconciled with ICA n_samples_.

Use MNE `ICA(method='fastica', n_components=None, random_state=97, max_iter='auto')`.
None uses MNE's near-total-variance selection with numerical/rank safeguards; it
is not a claim of one component per sensor. There is currently no clear evidence
favoring a smaller fixed decomposition for DAD. Log n_components_, PCA variance,
iterations, sample count and input summaries. No group example component indices
are transferred. See [MNE ICA](https://mne.tools/1.10/generated/mne.preprocessing.ICA.html).

P19's numerical reconciliation is documented separately in
[P19_ICA_NUMERICAL_RECONCILIATION.md](P19_ICA_NUMERICAL_RECONCILIATION.md).
After an unsuccessful higher-cap retry, source-bound diagnosis confirmed near
two-step recurrence. A deterministic restart of the saved seed-97 solution
converged in 114 restart-phase iterations, retaining exact PCA/whitening, all
69 directions/601,558 samples, original contrast/tolerance/cap and configuration.
This explicit initializer is optimizer history, not a fresh random draw.
Original weights/evidence are archived; all new ICs need human review. The
restart metadata appears in both reports and incompatible caches refuse silent
default refitting. Require successful publication and the joined current
readiness snapshot; no component/application/cleaned acceptance is inferred.

Review all component maps/properties/source traces with native MISC/frontal context.
The [final component-review guide](FINAL_ICA_COMPONENT_REVIEW_GUIDE.md) defines
the participant-level procedure. `03_C_ICA_review_components.html` is the human
review book with per-plot model-pinned commands; `03_C_ICA_review_report.html`
records processing/fit evidence without duplicating the figures. Raw-input
properties are arbitrary 2 s fit segments, never stimulus ERFs. The report-only
`review_final_ica.py` helper restores numeric axes and full concatenated display
coverage without refitting; legacy first-stretch panels remain labelled until
regenerated. Its separate TASK diagnostic uses accepted Stage 2/3 motion-onset
samples and the configured epoch window but creates no production epochs.
Interface v3 displays both native continuous examples and accepted-motion-onset
trial images/condition means by default for every IC. A single chunked native
projection per acquisition supplies all components; no per-IC source rereads,
refitting, native temporal filtering or exported epochs are introduced. The
native manifest records trial identities/drop reasons and exact BAD windows.
Missing continuous/task panels leave v3 diagnostic coverage incomplete. The
report-only wrapper's native mode reuses stamp-validated full-fit properties.
The accepted list is bound to exact model SHA-256; explicit [] is reviewed empty.
Changing annotations, input data or fit preparation requires refitting/review.
Apply never silently refits. Compatible cache reuse is portable across mount
paths and validates model contents as well as input/fit fingerprints.

Apply to each original-resolution annotated run separately. Save matched native
ICA-only before/after plots before final temporal processing, then spectrum-fit
line processing with every BAD/EDGE interval separating operations and no final
high-pass (2026-10-28 decision and implementation above), including spectrum-fit
for P38 run 1. Preserve its warning about a future high-pass. Historical `filter_edge_exposure`
intervals remain informational. No filtering operation deletes samples.
Preserve sampling rate, n_times, first_samp, measurement origin, channels/order,
bad metadata, original annotations, exact STI and MISC within FIF storage precision.
Preservation checks run on new and resumed outputs. Post-rejection models and
cleaned files must match the prefit model, exclusions, inputs and final settings.

Final ocular annotations are not removed after ICA. Under later reject_by_annotation,
BAD_blink/BAD_eye_movement/BAD_muscle and other separately reviewed BAD artifacts
reject overlapping epochs. Filter-edge exposure does not: even an exact legacy
`BAD_filter_edge` label is renamed in memory before Stage 5 creates epochs.
Incomplete epochs retain their usual boundary rejection. Trigger-150 handling
remains trial accounting, not annotation.
Stage 4 completion requires independently reviewed applied outputs and a shared
acceptance entry; successful jobs/implementation tests do not establish completion.

<a id="04-epoching"></a>

## Stage 5: Epoching

Entry point: `04_epoching.py`; runner selector: `epoching`.

Epoching first enforces the Stage 4 cleaning completion gate:

- `source_manifest.json` must contain `ica_manual` with `processing_pass=apply`
  and `status=completed`;
- ICA exclusions must be explicit, where `[]` is valid;
- every configured raw acquisition run must have exactly one completed lineage
  record and one local `<exact-source-stem>_ica.fif`;
- the manifest basename must match the raw source stem exactly. Cluster mount
  prefixes are treated as provenance and relocated through the configured
  derivative root.

Each cleaned run is read and epoched independently. Production selection comes
from the verified Stage 2 ledger. Trial intervals run from condition onset to
the next onset or an earlier behavioral-run/replay/warm-up boundary; the final
interval ends no later than the recording boundary. Evidence cannot leak into
the next acquisition or behavioral run.
Resolved transition-warning rows are accepted; unresolved transition evidence
is excluded by default.

Default trial-level removals are:

- trigger 150 gaze-break evidence;
- catch type 1 evidence: 116/117;
- catch type 2/question evidence: 113/118/201/202;
- ESC termination evidence: 152;
- unresolved transition evidence.

The trial-decision ledger must honor the accepted exclusions and recording
boundaries in the [Stage 2 requirements](#stage-2-trial-preparation-requirements),
retaining excluded candidates and their reasons before creating epochs.

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

<a id="06-drsa-neural-export"></a>

## Downstream: dRSA Neural Export

Entry point: `06_export_neural_for_drsa.py`; runner selector: `export_drsa`.

This downstream export consumes only completed Stage 5 epochs and epoching lineage.
Identity values come from the subject profile. The pilot profiles use 1-20
for subjects 01-03 and 1-40 for subjects 04-07; resolve each current-cohort
profile independently. One exact observed trajectory MAT is discovered from
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

<a id="07-matlab-drsa"></a>

## Optional Downstream: MATLAB dRSA

Script `07` is optional and independently rerunnable. It consumes completed
script `06` MAT/JSON exports and the completed `export_drsa` source-manifest
lineage. It does not modify epochs, neural exports, or raw data.

`07_run_matlab_drsa.py` validates the script `06` products and writes a
subject-level handoff:

```text
derivatives/reports/sub-XX/ses-01/07_matlab_drsa_handoff.json
```

That handoff is the only path/provenance contract consumed by MATLAB. It
contains script `06` neural paths, sidecar paths, trajectory model paths, identity
values, time axis, sampling rate, condition definitions, channel labels, and
output directories. MATLAB scripts must not hardcode V28, exactly 20
identities, local-only roots, or model filenames.

Registered optional analysis families are:

- `real`: `A_dRSA.m`, the main real-neural dRSA. Outputs are written under
  `analysis/sub-XX/ses-01/drsa/matlab_real_neural/07_dRSA/`.
- `alt`: `A_dRSA_alt.m`, single-model/sensor-group exploratory dRSA. For
  `occluded_deviant`, `--model-sources` can request `observed`, `predicted`,
  `prediction_error`, and explicit `pe_*` variants. The legacy
  `prediction_error` model is actual deviant position minus
  expected/predicted nondeviant position. The `pe_*` variants keep that
  source unchanged and write separate output folders.
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

After the selected subjects' step-07 outputs exist, a group report can be built
from the numeric MAT results. The commands in this section retain the historical
P01-P07 pilot selection; choose the analysis cohort explicitly for later runs.

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

Historical pilot cluster example:

```bash
DAD_MATLAB_BIN=/absolute/path/to/matlab sbatch --array=1-7%4 \
  meg/cluster/submit_drsa_analysis.sh
```

`DAD_MATLAB_BIN` must be an absolute executable path. Local MATLAB R2020a on
Apple Silicon may require `arch -x86_64`, which the script `07` CLI exposes as
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

Pilot channel coverage: sub-01 has EOG001/EOG002 in FIF and no saved EDF;
subjects 02-07 have no EOG and have saved EDF per behavioral block. For those
subjects, MISC001-MISC003 carry left-eye analog signals and MISC004-MISC006
are retained despite being near-flat. Use the current cohort's inventories
for its channel and EDF coverage.

- EDF-to-MEG synchronization is deferred and must not assume matching file
  starts or positional block/run correspondence.
- Trigger 150 records online gaze breaks independently of EDF availability.
- Trigger-150 trial rejection is separate from Stage 4B/4C physiological
  cleaning.

## Baseline Facts

- `04_epoching.py`: `tmin=-0.5 s`, `baseline=None`, linear detrending.
- `sanity_checks/01_dot_side_decoding.py`: `-0.4` to `-0.2 s`.
- `XX_ERF_manual.py`: default `-0.1` to `0 s`.

Both optional baseline windows fit inside the epoch. Their ITI compatibility
was verified in the pilot (0.5-1.0 s for sub-01; 1.0-1.5 s for subjects 02-07).
Assess the current cohort against its saved timing. Optional baseline
correction operates on copied inspection objects; production epoching uses
`baseline=None`.

## Optional Dot-Side Decoding Sanity Check

`sanity_checks/01_dot_side_decoding.py` is an optional diagnostic workflow and
is not registered as a required `run_pipeline.py` stage. It consumes completed
Stage 5 epochs plus `source_manifest.json`, validates that the epoching
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

Labels are model-derived from the exact subject trajectory file and configured
profile identity range. The validation path was checked with the pilot's
20-identity and 40-identity profiles. It fails on model identity mismatches
rather than truncating or synthesizing identities.

Historical pilot cluster example:

```bash
sbatch --array=1-7%4 meg/cluster/submit_dot_side_decoding.sh
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
cluster-forming `p=0.05`, and cluster-level alpha `0.05`. The seven-participant
pilot used exact sign flips; resolve the permutation settings for the selected
cohort. The group figure annotates the
cluster-forming p threshold and equivalent t cutoff, marks significant
clusters, and prints cluster p-values on the plot.

The default group display band is a 95% confidence interval around the mean.
SEM is available explicitly with `--interval sem`; the report records which
interval type was used. The HTML report contains all four group figures and
cluster tables first, followed by participant sections. Each participant
section begins with a three-condition source-trial
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
sbatch meg/cluster/submit_dot_side_group_report.sh
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
cluster-forming threshold is `p=0.05` and cluster alpha is `0.05`. The
seven-participant pilot used exact sign flips. Group plots show the valid-cell
group mean matrix with significant 2D clusters shaded and outlined in red.
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

Historical pilot dependency example:

```bash
tg_job_id="$(sbatch --parsable --array=1-7%4 \
  meg/cluster/submit_dot_side_time_generalization.sh | cut -d';' -f1)"
sbatch --dependency="afterok:${tg_job_id}" \
  meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

## Verification Contract

### Stage 4A SQUID implementation, 2026-10-28

252 targeted synthetic checks passed on local MNE 1.10.1, covering the new
detector/review path, actual SSS on synthetic coils with skipped artifacts,
Stage 4A prepare/review/apply/resume/archival, unchanged accepted SSS/report
preservation, native annotation/epoch timing, channel scope, review-zero/deletion
persistence and explicit reconciliation. Existing Stage 4B manual-save/recovery,
muscle/ICA, Maxwell QC, report, configuration, tuning-lock and wrapper checks
also passed. A stale plotting smoke test was updated to use the already-current
maxwell_qc.save_psds interface. Wrapper construction used harmless runtime stubs;
no scheduler submission occurred. Shell syntax checks passed.

No cohort recording, P04/P05 manual annotation/tuning file, open P05 run-2
review, accepted timing bundle, SSS output or model was changed by validation.
Scientific candidate review/performance, deployment and any affected-output
reconciliation/acceptance remain user-run work. Stage 3 documents require no
method/schema/hash change: this addition is later artifact handling only.

### Final temporal policy implementation, 2026-10-28

237 targeted Stage 4/configuration/annotation-save/EDF/report/wrapper checks
passed on both MNE 1.11.0 and 1.10.1 using temporary synthetic inputs only.
All 68 preprocessing/test Python sources compiled in memory, and the review,
apply and deployment wrappers passed shell syntax checks. Coverage includes independent
50/100 Hz spectrum-fit with no high-pass or margin, P38's retained future-HP
warning, native/acquisition-bandwidth preservation, BAD/EDGE isolation, archived
preview bytes, margin-free receipt creation, stale downstream evidence/output
rejection and unchanged-model reuse after final-only changes. Existing save-lock,
atomic replacement, recovery and no-reseed checks also pass. Synthetic naming,
missing sensor-layout and deliberately non-physiological ICA warnings remain.

The updated prerequisite checker validated all 36 existing 4A acceptances
read-only against current layered settings/resource contents and saved source/
SSS statistics; no cohort signal samples were loaded. Original 4A records and
Stage 1-3 evidence, P04/P05 annotation files/receipts and Stage 5 policy were
preserved. Cluster command construction was exercised with harmless runtime
stubs; no scheduler job, cohort processing, epoching or analysis was run.
4B/4C/applied-output scientific acceptance remains pending. Refresh/review
instructions are in the cluster runbook.

Use these checks when implementing or preparing the relevant pipeline step.
Keep completed cohort data review in the decisions record; rerun affected
checks when code, configuration, or inputs change.

- Compile all preprocessing scripts/tests and run the full pytest suite.
- Validate configurations and dry runs for the participants selected for that
  step, using the paths available on the execution machine. The cluster raw
  root is not mounted on the laptop.
- Check wrapper syntax with `bash -n` before submission.
- For cleaning changes, check MISC/bad-channel/ICA preservation with synthetic
  data and representative short real-data examples; inspect relevant headers.
- For epoching/export changes, check run isolation, trial exclusions, identity
  coverage, fixed channel order, and source lineage.
- For MATLAB dRSA changes, validate completed step-06 handoffs and run MATLAB
  `checkcode` before submission.
- For decoding changes, test line/time-generalization helpers and wrappers,
  including both supported identity counts and invalid model handoffs.

Complete the relevant implementation checks before processing that step;
scientific acceptance remains a separate review gate.

### Stage 4 implementation verification, 2026-10-22

- Full local suite: **378 passed** (MNE 1.10.1 / Python 3.13.3). Warnings include
  established FIF naming conventions, one deliberately non-physiological
  synthetic ICA convergence warning, and an existing plotting deprecation.
- All 36 cohort configurations validate; all preprocessing/test sources compile;
  all cluster shell wrappers pass syntax checks. New CLI argument handling and
  staged acceptance/application/resume were exercised with temporary FIFs.
- Regression checks cover all BAD categories, missing measurement dates,
  nonzero/large first_samp, reviewed-empty eyes/components, annotation endpoint
  serialization, auxiliary preservation, stale configuration/evidence/model
  rejection, retained medoid ties/invalid transforms, and isolated resampling.
- Changing huge artifacts inside excluded intervals leaves actual ICA fit
  samples identical. A historical synthetic example exercises the former 0.01 Hz
  FIR and spectrum-fit 50/100 Hz removal with explicit edge annotations.
- Bounded real-data checks used only 20-40 s windows of P04 run 1, P17 run 1,
  P24 run 5 and P38 run 1. Protected fit copies retained source coordinates.
  One P04 window exercised the actual SSS solver and FIF round-trip, preserving
  STI/MISC and native coordinates. Masks used in these checks are synthetic
  validation inputs, not cohort artifact decisions.
- Native context and before/after figures were generated; report rendering
  covers current manual ocular evidence and reuses exact-model review component
  figures in apply reports. The preliminary-source projection was tested with
  interactive windows replaced; hands-on Qt browser usability still needs the
  operator's first local review.
- Calibration contents match their pre-relocation hashes; pilot data/ZIP remain.
  All 36 accepted Stage 3 manifest/application review hash bindings still match.
  No cohort cleaning output, epoch, export, analysis or job was produced. Existing
  Stage 1-3 derivatives remain untouched. Cluster runtime execution and scientific
  cleaning acceptance remain future user-run review work.

### Historical pilot verification

The June 2026 checks covered P01-P07 and established the pilot interfaces.
Their participant outputs are now in the pilot archive described in the
[review guide](COHORT_PREPROCESSING_REVIEW_GUIDE.md#archives).

| Date | Recorded verification |
| --- | --- |
| 2026-06-16 | Compilation, cluster-wrapper syntax, and 15 targeted dot-side tests passed. Local configuration validation covered P01-P07; decoding handoffs covered 20- and 40-identity profiles. Subject/group pilot outputs were generated. |
| 2026-06-22 | Compilation, wrapper syntax, 30 targeted tests, and the 140-test suite passed (six existing synthetic-FIF filename warnings). All four decoding handoffs passed for P01/P04; a synthetic group report verified group-first ordering. Existing participant outputs were not regenerated. |

These dated checks describe pilot development. Current cohort acceptance is
recorded in the shared decisions document.
