# DAD cohort preprocessing: roadmap and review gates

For participant-level final component decisions use the [final ICA review guide](FINAL_ICA_COMPONENT_REVIEW_GUIDE.md). The component book is the human interface; the standard HTML records processing evidence.

2026-11-10 update: array 11197 passed the joined 35-participant/180-acquisition
report audit: 2,548 full-fit property panels plus 5,096 default native continuous
AND task-locked panels, with unchanged models and matched sidecars/output stamps.
The saved current-mask policy receipt already accepts Stage 4B for final fitting;
final component and cleaned-output acceptance remain pending. See the final ICA
guide for commands, interpretation and the separately authorized P19 numerical
reconciliation. Report completion is not component acceptance.

2026-11-11: P19's 10,000-step retry still failed; source-bound diagnosis
confirmed an optimizer near two-cycle. A documented deterministic restart of
the original seed-97 solution converged in 114 restart-phase iterations, keeping
all 69 directions/samples and the original PCA, contrast and tolerance. See
[P19 numerical reconciliation](P19_ICA_NUMERICAL_RECONCILIATION.md). Require its
completed publication record and the joined final-ica-readiness-20261111-v1
snapshot before the component handoff. All new ICs require human review; no old
indices or scientific acceptance transfer. Local unvalidated checkpoints now
preserve computation across share interruptions without relaxing validation.

Last reviewed: **2026-11-11**. Working cohort: **P04-P14 and P16-P39 (35)**;
P01-P03 and P15 are excluded from current analyses. P15 data/provenance are retained.
Raw `Attention/data/` remains read-only.

**Endpoint:** reviewed cleaned continuous FIFs and analysis-ready epochs, with
metadata, a final trial/rejection ledger, coverage/QC reports, and source lineage.
Analysis-specific exports and scientific analyses follow separately.

## Records and current status

- [PREPROCESSING_DECISIONS.md](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md):
  authoritative accepted choices, case rules, rationale, and pending decisions.
- [PIPELINE_SPEC_basic.md](PIPELINE_SPEC_basic.md): technical interfaces,
  parameters, outputs, and validation requirements.
- [Cluster runbook](../cluster/README.md): deployment, commands, and rerun behavior.

Stage numbers below refer to workflow/review gates. Script prefixes (`02`,
`03_A`-`03_C`, `04`, `06`) are legacy compatibility identifiers; retain existing
filenames, report paths and review keys. Use the specification's
[stage-to-script crosswalk](PIPELINE_SPEC_basic.md#stage-interfaces-and-implementation-status)
for dispatch and acceptance status. The naming clarification does not reopen
accepted decisions or require regenerating evidence.

Edit the decisions record in shared storage; no competing workspace copy is
maintained. Edit this guide and the specification in the workspace, then
synchronize their shared-storage copies. Deployment preserves the decisions record.

Execution, implementation validation, and scientific acceptance are separate.
The current decisions record explicitly accepts **Stages 1 and 2 on 2026-10-19**;
this was verified against [Reconcile stage 2 outputs](codex://threads/01a154b1-f77e-7532-8f7e-b0f95089ec70).
Stage 2's identity accounting is resolved. Historical reports retain their
run-time pending-review text and hashes; acceptance is recorded separately.
The immutable Stage 3 evidence begins with `s3-v001`; the coordinated
`s3-v002` was accepted on 2026-10-21 for all 36 participants, with separate
hash-bound review records for timing assignments, revised eligibility and
mandatory scopes. See the [application/API contract](STAGE3_REVIEW_APPLICATION.md)
before Stage 5. Cleaning lineage and ordinary epoching prerequisites remain.

## Current annotation-policy handoff (2026-11-07 decision)

Proceed with current annotations and defer further detailed interval inspection
until later findings justify it. Preserve every saved manual review dynamically,
including edits/deletions and reviewed-empty muscle sets. Other runs use the exact
current automatic ocular/muscle masks and existing source BADs. P15 is excluded
for insufficient remaining coverage (76/75/79 survivors and ten empty identity
cells); this participant-specific decision establishes no general threshold.
P10 stays in processing with its approximately 0.33 Hz concern open for native
before/after final ICA QC.

The [current-mask handoff](STAGE4_CURRENT_MASK_HANDOFF.md) is the active procedure:
whole-cohort validation -> separate, auditable Stage 4B adoption -> joined
reconciliation and explicit-policy processing acceptance -> final ICA fitting.
Automatic adoption never fabricates a manual SAVE or marks an untouched source
manually reviewed. This policy supersedes older requirements below for every
run to undergo exhaustive manual interval review or fresh EDF/SQUID review before
this pass. Source/model checks, accepted SSS, native coordinates and Stage 2/3
application evidence remain protected. All historical snapshots remain intact.

Final component decisions, application QC and cleaned-output acceptance remain
explicit later steps. No epoching follows merely from fit/job completion.
The new final-model component book supports efficient per-component navigation;
preliminary frontal ICA indices never transfer to final exclusions.

## Participant execution policy

From 2026-11-03, run all eligible independent participants in parallel by
default, bounded by allocated CPU/memory and concrete shared-resource limits.
Use disjoint processes or scheduler array tasks; preserve true prerequisite
barriers and a single final cohort publisher. Whole-table component-choice
validation must finish before any bulk import save. Independent annotation
counts and participant HTML books use process workers; shared summaries follow
completion. Worker limits, allocations and run progress must be explicit.
Existing jobs retain their original allocation/code; never restart or replace
manual work merely to change concurrency. See the cluster runbook's
[parallel execution policy](../cluster/README.md#parallel-participant-execution).

## Stage 1: sources, acquisition headers, and digital triggers

**Executed and accepted.** Purpose/prerequisites: reconcile raw sources, session
context, acquisition settings, and recorded STI101 against behavioral evidence.
Outputs: inventories, source maps, canonical event/alignment tables, and deep QC.
Review closed source/trigger discrepancies and recording boundaries; acceptance
and case rules remain in the [Stage 1 outcome](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md#stage-1-gate-outcome).

## Stage 2: trial preparation

**Executed and accepted for accounting under the recorded policy.** Prerequisite:
accepted Stage 1 evidence. Outputs: source/sample-keyed candidate ledgers,
recording/unavailable-tail coverage, eligibility and condition/identity/replay
counts. Review reconciled counts, selection effects, and application of accepted
cases; see the [Stage 2 outcome](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md#stage-2-trial-accounting-gate-outcome).
This acceptance does not settle visual timing, cleaning, final participant/run
inclusion, or broader epoch-rejection policy.

## Remaining roadmap

Stage 4 is one cleaning stage with three evidence-dependent review points.
Implementation or a completed job does not establish current-cohort scientific
acceptance.

| Stage | Purpose and prerequisites | Main outputs | New human review and acceptance criteria |
| --- | --- | --- | --- |
| **3. Cohort visual timing — accepted 2026-10-21** | Apply accepted `s3-v002` assignments and final eligibility to later epoching. | Hash-bound timing/application review records; original Stage 2 and `s3-v001` preserved. | Recheck lineage and accepted hashes at Stage 5; review only if inputs or contradictory evidence change. |
| **4A. Channels, SQUID review, Maxwell and filtering** | Prepare comparable continuous signals from accepted recordings/headers; independent of timing. Screen SQUID candidates after accepted channels and before SSS. | Per-run channel scores/lists, SQUID proposals/decisions, destination/filter settings, Maxwell FIFs, before/after QC. | Review flagged/unflagged channels and SQUID candidates; choose BAD_squid_jumps windows or channel exclusions, including reviewed zero. Accept resulting SSS application and existing acquisition constraints. |
| **4B. Artifact intervals** | Identify muscle and other unusable intervals on accepted 4A outputs. | Score traces, annotations, annotated FIFs, available ocular diagnostics. | Check flagged and nearby unflagged intervals; decide thresholds/annotation policy. Accept supported per-run annotations with preserved channel metadata. |
| **4C. ICA and cleaning completion** | Separate/reject physiological components using accepted 4A/4B inputs. | Models, component maps/properties/traces, explicit exclusions, cleaned FIFs and QC. | Review components and residual artifacts; a reviewed empty list is valid. Accept justified exclusions and before/after QC, preserving channels and recording coordinates. |
| **5. Epoching and preprocessing handoff** | Combine accepted trial accounting, revised eligibility, accepted timing and cleaning; finalize inclusion/rejection rules. | Epochs, metadata, final trial/drop ledger, condition/identity/channel coverage and lineage. | Apply the type-2 policy and eight exact whole-trial exclusions through a supported application without rewriting Stage 2. Reconcile coverage/balance and downstream exports. Accept applied timing/rules, traceable losses, justified channel handling and no unexplained exclusions. This closes preprocessing. |

## Dependencies and review boundaries

- **Keep trial preparation before timing:** it establishes recorded versus
  unavailable candidates independently of visual-lag assumptions. Timing adds a
  mapping while preserving original source/sample keys.
- **Assess diode signals on raw recordings before transformations.** Timing
  and continuous cleaning are independent branches; both must be accepted before
  final epochs. The working order is timing, then cleaning for manageable review.
- **Keep 4A/4B/4C reviews distinct:** upstream choices change downstream evidence
  and fits. Keep coding, tests, fitting/refitting, application and routine
  validation inside their parent stage. Applying an accepted ICA list does not
  create another approval gate by itself.
- **Preserve original trial accounting:** the separate application retains the
  Stage 2 ledger, eligibility values and hashes, then record revised final
  analysis eligibility separately. Reconcile every exact saved catch class,
  independent exclusion, unavailable tail, timing assignment and coverage row.
- **Keep review decisions exact:** the eight confirmed whole-trial exclusions
  are P24, P15, P28, P29 trials 10/12/13, P33 trial 39 and P38 trial 24.
  Only P24 removes measured events from donor inputs. Preserve every measured
  row, and apply no donor changes for the other seven decisions.
- **Keep P29 recordings separate:** retain 37 currently eligible always-visible
  trials from its partial half-block using their own measured code-71 onsets;
  the six type-1 catches and interrupted final trial remain excluded. Preserve
  the other source/block identities, including block 6. No averaging or
  weighting rule follows from this decision.
- **Retain P38's short block:** 16 of 20 code-71 onset rows from
  `sub38_run01_051026.fif` are eligible and measured, with their windows inside
  the recording. Keep the four existing gaze/termination/interruption
  exclusions and source identity; process these 16 like other eligible events.
  The user accepts P38's 13 ms onset only as its own measured edge; do not apply
  it to other events or infer a weighting rule from the block's smaller count.
- **Keep prompt timing distinct:** code 118 is a later catch-question prompt.
  Retain it for separate lower-priority behavioral review; do not substitute it
  for the motion onset or pool its lag distribution with motion events. The
  behavioral endpoint/clock is undefined.
- **Stage 5 prerequisite:** Stage 3 acceptance is recorded with exact timing
  and application hashes/scopes. Keep epoching blocked until accepted Stage 4
  cleaning lineage and the remaining production gates are ready.
- **Remaining choices:** per-run cleaning decisions, final filter-edge exposure,
  final epoch settings, rejection and participant/run inclusion. No EDF synchronization
  study is required for Stage 4. Session notes are context, not automatic exclusions.

## Accepted Stage 3 handoff and remaining roadmap

The `s3-v002` motion-onset and revised eligibility scope was accepted on
2026-10-21. The cohort has 23,963 policy-eligible trials (21,890 original +
2,079 newly eligible catches - six newly excluded original trials), after all
eight exact sequence decisions. Five P37 saved MAT outputs in the unavailable
tail remain unrecorded. The 2,079 new catch onsets are in-bounds without
trial-end overlap: 1,718 measured and 361 estimated. Donor medians for codes
71/81/91 are 19 ms using an ordinary unweighted median over 22,120 unique
QC-usable measured events (n=7,341/7,422/7,357). P24's five donor events were
removed; affected counts are 48 for code 81 and 89 for 111/112/114/115, while
medians remain 19 ms. Codes 116/117 have six donors each and no estimate.
P04-P07/P17 use estimated +19 ms onsets with transfer error unknown. P33's
601 remaining original-eligible measured onsets were accepted. Code 118 prompt
timing remains deferred behavioral work with undefined endpoint/clock.

Thirty-six accepted review JSONs under `s3-v002/reviews/` bind exact participant
manifest hashes, application hashes, accepted assignment bases and required
review scopes. `04_epoching.py` automatically discovers the bound
`application.json` beside each bundle. `--eligibility-application` may supply
an explicit path. Production may use the existing per-subject
`epoching.timing_review_path` or `--timing-review` after accepted cleaning
lineage is ready. Stage 4 cleaning remains independent; no cleaning, epoching,
step-06 export or analysis was run for this handoff.

## Preliminary annotation retention and review prioritization

This describes the earlier triage preparation. The 2026-11-07 current-mask
policy above supersedes its exhaustive manual-SAVE prerequisite for this pass.

On 2026-10-29 the user adopted a preliminary annotation-retention pass before
completing detailed interval review throughout the cohort. This supplies the
review order; final Stage 4 acceptance and production epoching retain their
existing gates. Component choices and automatic masks are explicitly provisional.

Preserve the completed P04/P05 reviews and P06 runs 1-3. P06 runs 4-5 are ordinary
unreviewed inputs: the user stopped their manual review for now. Every execution
checks current exact-source receipts/files, so later manual saves also take
precedence. Never replace manual deletions or reviewed-empty muscle sets with
automatic detections. Legacy ocular-only receipts keep existing automatic muscle.

1. Prefer the offline participant books approved 2026-10-30:
   `build_ocular_component_reports.py` makes ONE self-contained HTML per
   participant with internal run navigation, maps and scrollable/zoomable traces.
   Generate in advance locally or with the report-only cluster wrapper. Existing
   manual work/component choices are preserved and listed. Record choices and
   notes per RUN, send them in chat or use a copy of the decision template.
   `import_ocular_component_choices.py` validates the complete table against
   exact reports and current models before --apply saves new choices. No
   missing field is inferred as NONE; conflicts never overwrite prior work.
   The previous local component-only viewer, `review_ocular_components.py`,
   remains available. Inspect
   cached maps/native sources and save independent zero-based blink/eye choices,
   explicit NONE, or DEFER. Valid completed choices resume automatically; manual
   work is skipped. Missing/stale caches need targeted preparation, not local fits.
2. Run `annotation_triage.py`. It combines preserved manual decisions or current
   automatic ocular/muscle/source BAD intervals with verified Stage 2 and accepted
   `s3-v002` timing/eligibility. EDF/SQUID proposals contribute nothing. Existing
   source BAD annotations remain operative; filter-edge exposure remains informational.
3. Inspect the report under `derivatives/annotation_triage/reports/`. The separate
   `render_annotation_triage_report.py` renders these saved tables as self-contained
   embedded-figure HTML in `annotation_triage/visual_reports/LABEL/`: a cohort
   overview and one participant book, including three identity barplots, run
   retention, overlapping causes, exact counts and ocular coverage. Input-error
   participants are unranked. Optional before/after comparisons require unchanged
   trial accounting and epoch coordinates. Prioritize
   the lowest surviving condition count, checking identity gaps, retention rates,
   incomplete windows and ocular scoring coverage. Input errors remain unranked.
   Missing choices, NONE and NaN/short stretches never mean artifact absence.
4. Open the exact counted mask with `manual_eye_annotations.py --triage-proposals`.
   Perform ordinary whole-source inspection/editing/SAVE. Existing manual work
   is reopened normally instead. Recompute selected participants or regenerate
   the cohort snapshot; compatible automatic detector masks are reused. Rendering
   an old count snapshot alone does not incorporate newer annotations. The user
   can run `submit_annotation_retention_report.sh` with `DAD_TRIAGE_RECOUNT=1`
   and a new label to recompute the cohort and render HTML in one job, without
   repeating the completed component-choice import. Preserve all earlier snapshots.

Counts precede final ICA, line treatment and amplitude rejection, so they are
conditional on annotation quality and cannot alone establish participant exclusion
or replacement. Retain consistent artifact criteria and later signal-quality QC.
All triage writes have a dedicated namespace; production annotations, tuning
receipts, accepted SSS and Stage 2/3 evidence remain inputs. Complete manual SAVE
from a triage mask records an explicit EDF opt-out, also retained on later
reopening and 4B rebuilding. It does not accept 4B or final ICA automatically.
See the [runbook](../cluster/README.md#preliminary-annotation-retention-triage)
for local commands, outputs and targeted review usage.

The offline books use compact peak-preserving display bins with exact native
peak offsets, per-component maximum-absolute normalization and explicit gaps.
They change no scientific source sample or detector input. Viewer controls run
entirely from embedded data, so navigating to another run reads no FIF. An
immutable per-participant manifest binds the exact report, model and source;
the importer revalidates those inputs, including cross-mount cluster/local
usage. These are component-choice reports, not full manual artifact review,
ICA component-removal decisions or cohort scientific acceptance.

## Next stage: Stage 4 cleaning review

Implementation/settings approved 2026-10-22; cohort cleaning is not yet accepted.
The bounded preflight found 185 acquisitions, 42.75 h, native 1000 Hz. P25/P29/
P35/P37/P38 have six acquisitions; partial source boundaries remain unchanged.
All current static transforms are available. P17 and P24 runs 1-4 have no MISC
and retain irreversible 0.1 Hz acquisition high-pass; all other recordings have
six eye MISC channels, with later recordings also carrying MISC008 photodiode.
Channel presence does not establish usable analog eye information.

The selected settings are DAD central_run (actual-run digitization-cloud medoid),
SSS only, no additional final high-pass (user decision 2026-10-28), spectrum-fit
50/100 Hz removal, no final 40 Hz low-pass, and final ICA n_components=None.
The no-high-pass runtime/lineage update is implemented and synthetically tested.
Spectrum-fit also runs for P38 run 1 under the user's 2026-10-28 clarification.
The annotation aid does not
change final fitting: routine ocular intervals and all other BAD categories
remain excluded. Preliminary frontal ICA is supplementary ocular annotation
context only; EOG is legacy. EDF is discontinuous optional trial-segment context,
not reliably globally aligned, and needs explicit unrecorded gaps. No EDF study.

The approved 2026-10-23 addition generates trial-mapped EDF blink candidates in
separate sidecars during a guarded 4B resume, reusing unchanged muscle outputs.
Candidates are not BAD and exclude nothing. The local browser loads editable
candidates plus normalized, temporary EDF gaze/pupil traces with NaN gaps and
explicit coverage, alongside MISC/frontal MEG. Complete-source SAVE explicitly
promotes retained candidates to BAD_blink; long events (>0.5 s), ambiguous trials
and uncovered periods need inspection. Changed guidance requires reconfirmation;
deleted candidates are not reintroduced when unchanged work is reopened. No
global offset, fixed lag correction, continuous signal reconstruction or new ICA
analysis is performed. See the cluster runbook for dependency/refresh/review commands.

The 2026-10-26 approved revision caches the preliminary aid on the cluster, using
ALL good frontal mag/grad sensors and at most ten components. Existing known
muscle/severe BADs protect the fit; ocular examples remain available for this
supplementary fit only, with no fresh manual severe-artifact pass. Locally load
maps/source traces, choose blink components followed by eye-movement components
(`0` = `ICA000`, first), and inspect both categories beside EDF/MISC/frontal
signals. Each prompt accepts comma-separated numbers/names or explicit NONE;
blank/invalid answers retry. NONE disables only its detector, not manual review
or the possibility of artifacts. A component may be selected for both groups.
Blink proposals remain BAD_blink_auto. Independent duration-based eye proposals
use BAD_eye_movement_auto with their own scores/settings/provenance and vetoes.
Five robust standard deviations, 0.25 s peak spacing and +/-0.25 s windows are
review starting points. Broad/flat/unscored cases and unflagged periods still
need inspection; an unclear component is not evidence of no blinks. No final ICA
settings/removal policy changes. Source-specific caches check lineage before reuse.
The eye aid uses 100 ms median smoothing, threshold 4, 100 ms minimum excursions/
merge gap and 50 ms padding within native finite stretches. It does not reuse
blink peak windows. Its cached 1-40 Hz filtering limits long position plateaus;
unscored gaps and unflagged periods still need inspection. See the specification
for assumptions/coordinates and the runbook for parameters/review commands.
Compatible older blink-only caches are reused; absent eye selection means NONE
in memory without migrating receipts or establishing whole-source review.

The local browser now defaults to accepted Stage 3 corrected onset lines plus
green EPOCH_preview annotations (-0.5 to +2.667 s). Their records are checked
through the existing Stage 2/3 validators; only accepted timed trials receive
guides, before later artifact/amplitude rejection. Toggle the label with a ->
Select Visible in Qt (visibility checkboxes in Matplotlib). Preview edits are
discarded on close and never saved as artifacts or used to reject epochs.
--no-timing-marks disables all timing guides; --browser-backend matplotlib gives
thin default traces. Incomplete windows are counted and clipped only for display;
their usual Stage 5 rule and separately reviewed artifacts remain. Filter-edge
exposure stays informational. These display changes require no cluster rerun.

Code 150 gaze_break is an opaque orange marker at its original recorded sample,
above all plot layers and sharing --onset-linewidth. Exact ledger evidence links
its owning acquisition/trial/replay. The already excluded trial gets a separate
pale-orange EPOCH_invalid_preview using stored hash-bound Stage 3 timing for the
intended -0.5/+2.667 s window, without a green guide. Matching/timing gaps are
reported; no new correction or exclusion is inferred. Both preview types toggle
independently, are never saved as artifacts, and disappear with --no-timing-marks.

This browser also reviews muscle labels: SAVE retains edited BAD_muscle alongside
all other BAD categories. A new muscle-review field binds the complete replacement
set to its source/settings, including reviewed zero. Reopening and rebuilding 4B
must honor deletions/resizing/additions without restoring automatic intervals.
Legacy ocular-only receipts retain automatic muscle labels and do not imply zero.
Review rebuilt automatic score diagnostics and applied annotation counts; changed
detector settings require reconciliation, and changed muscle fit masks require
only affected preliminary-ICA cache rebuilds. Other preservation checks remain.

BACK reopens both selections on the same model; unchanged groups retain results
and changed nonempty groups ask separately before regeneration. Declining keeps
the prior scores/provenance distinct from a new pending selection. REDETECT asks
BLINK or EYE and adjusts only that group's threshold. Manual edits/deletions and previously saved work
are preserved. MISC usefulness/observations are asked AFTER viewing. SAVE explicitly
reviews retained auto intervals and only then saves BAD labels and model/detector
provenance; DISCARDALL leaves previous files/receipts unchanged. Blank/typo actions
retry. A failed receipt write remains recoverable with --overwrite
--recover-annotations, without reseeding deleted EDF proposals. See the runbook
for preparation-only pilot/cohort commands; independent review and 4B acceptance
remain necessary before any final 4C fitting.

The group reference compared was commit
e584e5b49fa2da4e8492319baa4198910db8c039. Borrow its simple processing sequence,
retaining DAD layout, accepted decisions and preservation checks. The previously
proposed head-origin alternative was probably a first-pass oversimplification
of DAD central_run and is not selected. TWA head position remains a future lab
option requiring continuous position data unavailable here.

Stage 4A was accepted on 2026-10-23 for P04-P39 (36 participants, 185 runs).
All participant-union channel sets and central_run destinations were confirmed;
36 existing tuning files hold report/configuration-bound acceptance receipts.
The authoritative shared decision record contains per-participant distances and
the historical final-filter edge policy/P38 warning and its later clarification.
4B/4C remain unaccepted; no 4A refresh/reacceptance is required for final policy.

**SQUID extension approved 2026-10-28:** this belongs inside workflow Stage 4A
(legacy script `03_A_MaxFilt.py`), after accepted channel decisions and before
SSS. The accepted participant unions are skipped during screening; no automatic
channel re-review is required. Proposals are non-BAD until explicit review;
confirmed intervals use exactly `BAD_squid_jumps`, while channel exclusions use
bad-channel metadata and retain union as the default. A run-only exception is
explicit and must preserve the final-ICA channel-basis contract. Preparation and
local review do not change existing SSS or 4B files and do not create a new
numbered scientific gate. The no-rerun exception above concerns final temporal
policy only: changed operative SQUID treatment requires affected 4A rebuild/QC
and renewed acceptance before dependent work. Reviewed zero/unchanged SSS keeps
the accepted reports and manual review bindings.

Preserve completed manual BAD reviews for **all P04 runs and P05 run 1**.
P05 run 2 was open for manual review when implementation was approved; do not
interrupt, reload or replace that process or its annotation/tuning files.
Immutable SQUID sidecars and locked exact-source tuning merges keep these
independent reviews separate. After a real SSS change, retain saved 4B intervals
as candidates for source-specific revalidation, including manual deletions and
reviewed-empty results; do not reset/reseed or automatically reapprove them.
Stage 4B displays Stage 4A SQUID context, preserving it once from the SSS source;
changes to that category return through Stage 4A. All BAD masks and later epoch
rejection honor it. User-submitted processing/review steps are in the
[runbook](../cluster/README.md#stage-4a-squid-candidates-before-sss).

SQUID review usability now includes a focused flagged-sensor view with explicit
candidate/channel labels and local physical scales. Overlapping preview shading
is displayed once. DETAIL provides per-sensor context/edge/difference evidence;
the independent `--mode inspect --candidate ID` viewer can run beside an already
open editable session without discarding unsaved work. It writes no decisions
or preparation artifacts. Keep muscle-like flags separate from confirmed SQUID
decisions; dismissal here does not alter reviewed Stage 4B BAD intervals.

1. For the accepted cohort, reuse the 4A channel/SSS evidence to prepare and
   review SQUID candidates. Finish those decisions before scheduling the pending
   4B refresh. Candidate/display/detector changes and reviewed-zero/unchanged SSS
   need no SSS rerun or repeat manual BAD review. Actual interval/channel changes
   require archived, affected 4A rebuild/QC/acceptance and explicit revalidation
   of saved 4B work only where its SSS changed. A union channel addition can
   affect every run. Never automatically rebind or reseed existing manual work.
   See the [coordinated refresh](../cluster/README.md#coordinate-the-pending-4b-refresh-with-squid-review).
2. Once current 4A is settled, user runs provisional 4B muscle evidence and the optional cached frontal blink
   aid. Locally choose components, inspect/edit proposals and review every source's eyes.
   Rebuild 4B with reviewed eyes/thresholds, inspect flagged/unflagged intervals,
   exact coordinates and final-filter edge exposure, then reconcile acceptance.
3. User runs 4C final fitting only after accepted 4B. Inspect components and
   native MISC/frontal context, record model-bound exclusions (reviewed empty
   is valid), then user applies that exact model to the native runs. Inspect
   ICA-only before/after and final filtering/preservation before acceptance.

For 4A channel review, the initial per-channel raw PSD includes sensor-location
insets and highlights bad-channel curves. Noisy/flat score heatmaps identify this run's detected
channels in red `[D]` and applied-only harmonization/manual choices in amber
`[A]`. Inspect these beside the reproducible random 20-second raw MEG trace
window and the final per-channel raw/SSS PSD, which excludes applied bads from
both sides and uses matching scales/windows. Threshold crossings are evidence,
not new channel decisions; missing scores are not zeros. The random window
can miss intermittent faults, so use the full-run heatmaps to guide further
inspection. `maxwell_qc.py` / `DAD_MAXWELL_QC_ONLY=1` refresh these displays from
saved evidence without rerunning detection/SSS or editing scientific outputs.

The user superseded the 0.01 Hz final FIR on 2026-10-28: omit the additional
high-pass for the roughly 3 s corrected-onset analysis. Slow activity is expected
to vary little within an epoch, but short epochs do not establish negligible
drift; later QC must inspect offsets/slopes and trial-order effects without
implicitly selecting baseline correction or detrending. The old 165 s exposure
reports remain historical information and never reject otherwise eligible
epochs. No high-pass margin is needed; review the
remaining 10 s spectrum-fit operation at BAD boundaries/short stretches.
P38 run 1 (`sub38_run01_051026.fif`, 137 s) receives spectrum-fit too, as clarified
by the user on 2026-10-28. Its original exact-source tuning reason is preserved
as a warning/guard for any future additional high-pass. Incomplete epochs and
reviewed BAD artifacts retain their usual rules. Final-policy previews and 4B
evidence must be refreshed before new acceptance; omit the obsolete margin
argument. Compatible annotation FIF/CSV files can be reused with guarded resume;
changed manual masks need an explicit affected 4B rebuild after review saves
finish. Superseded exposure reports are archived automatically. P04/P05 saved
manual FIFs/receipts and the locked atomic save/recovery behavior are preserved.
Fit preparation
is independent: 1-40 Hz on separate good stretches, edge guards, then separate
200 Hz resampling, with source/sample maps and actual fit count reconciliation.

See the [specification](PIPELINE_SPEC_basic.md#review-and-apply) for exact input,
output, masking and cache contracts, and the [runbook](../cluster/README.md#stage-4-cleaning-separate-4a-4b-and-4c-reviews)
for commands. Stage 4 runner snapshots preserve accepted Stage 1 evidence;
source-manifest changes are additive cleaning lineage only.

Any SSS/source/bad-channel/destination change invalidates affected annotations,
participant ICA and cleaned outputs. Annotation/mask/fit-setting changes require
participant refit/component review and new cleaned outputs. Exclusion-only
changes reuse a compatible prefit model. Final-filter-only changes regenerate
final outputs/QC without refitting unchanged fit inputs. Report-only changes
regenerate presentation only. Keep old manual intervals as candidates for
revalidation after changed SSS rather than deleting reviewed work.

Completion requires reconciled acceptance of applied continuous recordings,
exact cleaning lineage, documented limitations and a Stage 5 handoff. Preserve
Stage 2 and both Stage 3 versions/applications/reviews, source maps and boundaries.
No trial eligibility edits or timing shifts affect continuous FIF coordinates.
Later epochs use -0.5 to +2.667 s around corrected motion onset. Whether incorrect
catch answers should exclude trials remains a separate Stage 5 decision.

## Review practice and optional downstream work

Each handoff states purpose, prerequisites, command, outputs, inspection criteria,
and the new decision. We inspect independently, reconcile findings, and record
acceptance in the decisions file. Later stages verify accepted rules without
renewed approval unless inputs or contradictory evidence change. Revise this
roadmap when findings change dependencies; detail only the next stage.

Step `06` dRSA export is an analysis-specific data handoff. MATLAB dRSA, ERF
analyses, decoding/time-generalization and group/statistical reports are optional
downstream analyses with their own prerequisites. Diagnostic plots needed for
preprocessing QC remain within their parent stage.

## Archives

Historical pilot products and the first Stage 1 attempt remain under
`derivatives/pilot/` and `derivatives/previous_runs/stage1_2026-10-15/`, with archive
manifests. They provide context; historical cleaning choices are not transferred
automatically to the current cohort.
