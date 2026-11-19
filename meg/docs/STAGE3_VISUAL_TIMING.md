# Stage 3: cohort visual timing implementation and review

The cohort Stage 3 review was accepted on 2026-10-21 for the exact `s3-v002`
proposals and separately hash-bound eligibility applications. The original
`s3-v001` bundles and accepted Stage 2 ledgers remain immutable. Thirty-six
accepted review JSONs in `s3-v002/reviews/` bind each participant manifest hash,
application hash, assignment bases and mandatory scopes. Scientific acceptance
is recorded separately from the proposal manifests. The shared
`PREPROCESSING_DECISIONS.md` is the formal acceptance authority.

The accepted scope covers motion-onset timing (codes 71/81/91), the restricted
type-2 catch eligibility rule for conditions 81/91, and eight exact
participant-analysis exclusions. It does not define prompt timing for code 118,
provide P24 estimates for codes 116/117, or resolve transfer error from measured
donors to unmeasured participants. Stage 4 cleaning remains independent;
epoching and downstream exports were not run as part of the review.

## Scope and preservation

`visual_timing.py` is standalone; it never calls `run_pipeline.py`, the P24
prototype, cleaning, ICA or epoch creation. Raw `Attention/data` is read-only.
Stage 1 events/alignment/source maps and accepted Stage 2 bundles remain intact.
The completed Stage 2 loader verifies scientific evidence, policy and output
hashes before timing. Exact accepted exclusions remain excluded, and unavailable
behavioral tails acquire no samples or candidates. Stable identity is participant
+ exact source FIF + original absolute raw onset sample, including `first_samp`.

The original Stage 3 producer and `s3-v001` maps preserve the Stage 2
eligibility snapshot. The separate `visual_timing_application.py` creates a
hash-bound overlay without rewriting either source. Task 2 adds the reviewed
onset method and publishes one coherent `s3-v002` for all 36 participants.
`04_epoching.py` consumes the reviewed bundle plus its bound application only
after loading the accepted timing JSON and mandatory scientific scopes. Do not
edit Stage 2 or `s3-v001`, or treat proposal manifests as acceptance records.

Validate current inputs and headers without loading signals or writing outputs:

```bash
python meg/code/preprocessing_basic/visual_timing.py \
  --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-11.yaml \
  --version s3-input-validation --validate-only
```

Saved MAT Conf/runResults are read as contextual evidence. Exact block/run/trial/
replay keys attach guard, boundary, missed-flip and dot-offset diagnostics;
condition-code disagreement stops processing. Behavioral/PTB times never replace
recorded observations. Configured stimuli V30, schedule v23 and runtime versions
are separate. P04-P07 used v20; P08 block 1 used v24, later blocks v25; P09-P39
used v25. A runtime/profile mismatch stops processing. No current launcher is
used to infer a participant's historical runtime.

## Candidate detector and scientific limitations

`meg/params/stage3_visual_timing.yaml` is separate from the hashed Stage 2 policy.
Every value remains a proposal. The implementation reads 250-ms raw windows
(-75 to +175 ms) of STI101/MISC008 at native rate. Ordinary baseline is the
median at -50:-5 ms; the candidate search is -20:+80 ms. The bilateral search
exposes unexpectedly early/late responses rather than forcing positive lags.

Recording-level polarity comes from ordinary condition-onset controls in those
same windows. Ambiguous polarity abstains; it is not chosen independently for
each event. A sustained provisional response identifies the first rise; its
first 8 ms define response amplitude. First fractional crossings at 10/25/50/
75/90% remain on the native integer grid. The 50% crossing is the candidate
primary edge. Negative lags are retained and may reveal baseline contamination.

Candidate flags are amplitude <5 mV, baseline SNR <8, 10-90% span >3 ms,
ordinary baseline variation >20% of amplitude, pulse-envelope deviation >4 ms,
multiple rises, missing fall, truncated windows or nonfinite data. Flags retain
candidate samples but do not enter automatic measurement assignments/donors.
These thresholds are not accepted participant QC criteria. Rounded levels,
dynamic range, amplitude, noise and waveform plots support signal-usability
review; channel presence or patch enablement alone never establishes usability.
Clipping/saturation and marginal traces require visual review; no ADC-rail
calibration is inferred from a participant waveform.

Warm-up 153 uses a post-anchor low-trough baseline and a subsequent rise because
the previous plateau can remain high at the digital edge. It is diagnostic only.
The bounded pilot exposed the same contamination at the first retained onset
after warm-up, including false selection of the preceding warm-up transition.
For condition onsets only, accepted canonical history must identify a preceding
153 within 150 ms before the explicit `post_warmup_low_trough` method is used.
That method uses the 10th percentile at 0:30 ms and a subsequent rising edge
searched at 0:80 ms, excluding the previous physical pulse. It is recorded
per event and remains a candidate definition for independent review. Ordinary
events retain the bilateral search; arbitrary bad baselines are not repaired.
Expected pulse envelopes use the saved refresh/pulse settings, not burst counts.
Condition 71/81/91, mid-trial 111/112/114/115/116/117, question 118 and warm-up
are kept distinct. No global lag is applied. Identity, response, gaze, replay
and termination codes do not receive fabricated diode observations.

Recording/code summaries expose distributions, tails and linear drift as a
diagnostic, never an automatic correction. Saved boundary/rollover/restart fields
allow stratified inspection. Genuine timing tails should stay measured; detector
failures are different. Substantial drift, discontinuities, polarity changes or
new failures require a revised version/segment method after review, not silent
pooling or exclusion. The initial implementation does not automatically fit
change points or extrapolate drift into missing recordings.

`physical_intervals.csv` compares measured same-context diode endpoints with
saved motion-frame differences; no missing endpoints are reconstructed.
Question UI frame indices are not treated as motion frames. No dedicated
optical dot-offset endpoint is inferred from v25's saved PTB offset diagnostics.

Patch onset is not the moving dot's exact luminance onset: spatial scanout,
display response and patch-to-dot transfer remain limitations. Historical v21
and v22/v23/v24 test results motivate QC but are not participant calibrations.
Appendix B's v20 counterfactual is never a default correction or donor.

## Restricted estimates and fallbacks

- **P04-P07:** The accepted `s3-v002` proposal assigns estimated +19 ms
  motion-onset shifts derived from the ordinary median of unique QC-usable
  measured events for codes 71/81/91. The assignments are labelled estimated,
  never measured. These participants used v20 and have no diode measurements.
  Cross-runtime transfer from v24/v25 to v20, participant transfer and
  patch-to-dot error remain unknown; 19 ms is not evidence that this estimate
  is closer than zero for an individual onset.
- **P17:** The accepted `s3-v002` proposal assigns estimated +19 ms motion
  onsets from the same measured-event median. All five P17 contexts use v25,
  matching later donor runtime, but participant/display transfer remains
  unknown. Preserve the irreversible 0.1 Hz online high-pass and lack of P17
  photodiode measurements. The estimate is not a measured or trial-specific
  physical lag.
- **P11 blocks 1-2:** Confirmed unusable; only its own measured later blocks
  starting at block 3 may donate. Early patch position [193,0] differs from
  donor blocks 3-5 [245,208]. Spatial-transfer uncertainty stays unknown.
- **P24 acquisitions 1-4:** Only measured acquisition-5 events may donate.
  Preserve irreversible 0.1 Hz online high-pass/no-diode limitations. Run-5
  measurements cannot validate transfer across that settings change.
- **P09 acquisition 2:** Suspected detachment is a lead, not a confirmed
  unusable recording. Assess individual/temporal signal coverage and neighboring
  controls.
- Other unusable/flagged edges remain unresolved. No cohort mean replaces
  event-specific measured assignments.

The present `s3-v001` P04-P07/P17 rows retain their original 0-shift digital
fallback as historical evidence. The accepted method direction requires a new
supported estimator version, recalculated estimates/maps/reports and coverage
review before timing acceptance or epoching. For P04-P07, a separate requested
future bench measurement should reproduce the historical v20 trigger/draw/flip
order with added optical marker instrumentation and document any instrumentation
or hardware effects. That bench work is not performed and cannot retroactively
calibrate historical participants.

Estimates are condition/code-specific. Candidate minimum is 10 usable donors
per code; this is an explicit proposal, not the old adapter's 20/80% rule.
P11's location is the median of donor-block medians; P24 has one donor block.
The event-weighted median is also reported so weighting effects are visible.
Donor eligibility ordinarily uses waveform quality, independently of
participant-analysis selection. Accepted Stage 2 exclusions remain excluded
regardless of their diagnostic signal quality. Any additional donor-input
exclusion must be an explicit, separately recorded decision: the planned P24
application removes all five events from its named whole trial, while P15/P28
participant-analysis exclusions do not remove measured donor evidence.

Deterministic hierarchical bootstrap resamples donor blocks, behavioral-run
clusters and events (2,000 resamples, seed 20261020). Reports distinguish the
location CI, observed donor 5-95% spread, SD, held-out donor-block errors and
unknown context transfer. With one donor block, validation is within acquisition
only, with behavioral-run holdout errors reported. Sparse pilot donors cannot
establish full calibration. A narrow CI of the
median is not a target-trial prediction interval and does not remove jitter.
Condition pooling, alternative weighting, outliers and transfer assumptions
remain independent-review decisions; changing them requires a new version.

## Versioned outputs

Root: `Attention/derivatives/visual_timing/VERSION/sub-XX/ses-01/`.
Version identifiers cannot traverse directories. A bundle directory is reserved
exclusively; there is no overwrite. `--resume` requires matching code, method,
inputs and intact completed outputs. Incomplete bundles remain inspectable and
cannot masquerade as completed; choose a new version for a fresh attempt.

| File | Contract |
| --- | --- |
| `manifest.json` | Written last; schema, mode, version, code/method/input and all output hashes; scientific status is recorded separately in hash-bound review JSONs. |
| `contexts.json` | Saved runtime, patch/refresh settings and exact trial diagnostic joins. |
| `event_measurements.csv` | Recorded marked events, original anchor, candidate edge, signed lag, thresholds, waveform QC and measurement status. |
| `event_timing_proposals.csv` | Separate measured/estimated/fallback/unresolved assignments; warm-up diagnostic only. |
| `trial_timing_map.csv` | Every original Stage 2 key and exclusion, proposed sample, applied shift, assignment basis/status, uncertainty and window coverage. |
| `estimates.json` | Exact donor event keys, estimator, bootstrap and transfer limitations. |
| `recording_summary.csv` | Per-source/code usability coverage, acquisition facts, lag distributions and drift diagnostics. |
| `physical_intervals.csv` | Physical same-trial intervals and available frame-based comparisons. |
| `delay_violin_with_dots.png`, `lag_by_trigger.png` | Lag consistency and magnitude by trigger factor and exact recording/runtime; measured dots, median/percentiles, mean +/- sample SD, counts and flagged candidates. |
| `lag_over_time.png`, `signal_usability.png`, `threshold_sensitivity.png`, `measurement_coverage.png`, `context_comparisons.png` | Recording-separated stability, signal quality, threshold sensitivity, measured/failed counts and saved-context comparisons. |
| `physical_intervals.png`, `epoch_coverage.png`, donor figures where applicable | Physical sequence sanity, original/proposed epoch coverage, restricted donor locations/spread/confidence and held-out errors. |
| `physical_intervals_violin.png`, `interval_residual_vs_duration.png`, `interval_residual_vs_time.png`, `interval_residual_by_pair.png` | Additional full-range interval distributions/points, duration/time diagnostics and categorical endpoint-pair distributions; the original interval/bar figure is retained. |
| `interval_details.csv`, `interval_statistics.csv`, `interval_transition_statistics.csv`, `interval_tail_pairs.csv` | Verified endpoint joins, source/pair means and sample SD, distinct-trial counts and an exact-key inspection inventory; no new exclusions or assignments. |
| `waveforms_recording_*.png`, `diagnostics.png`, `lag_sensitivity.png` | Native waveform gallery; the last two names retain compatibility as aliases of the first gallery page and categorical lag bars. |
| `waveform_windows.npz`, `diagnostic_index.csv` | Reinspectable selected native STI101/diode windows; all participant pilot reads are cached. |
| `report.html`, `report.md` | Primary self-contained offline review report and Markdown companion; every figure has axes/units, derivation, observed findings, decision guidance and an example. HTML embeds figures and exact evidence downloads. |
| `lag_statistics.csv`, `correction_summary.csv`, `signal_statistics.csv`, `context_statistics.csv`, `donor_summary.csv` | Inspectable counts, means, sample SD, medians/percentiles and actual proposed shifts; missing dispersion is explicit. |
| `coverage_summary.csv`, `coverage_review.csv`, `event_qc_review.csv`, `tail_events.csv`, `timing_sanity_checks.json`, `report_summary.json` | Review tables, complete QC/tail lists, assignment/boundary checks and machine-readable report explanations/provenance. |
| `code_visual_timing.py.txt`, `code_visual_timing_utils.py.txt`, `code_visual_timing_report.py.txt`, `method.yaml` | Frozen implementation and method alongside their fingerprints. |

For each measurement, `lag_samples = diode_edge_sample_raw - digital_sample_raw`.
Assignments add the signed shift. Recording-relative seconds subtract `first_samp`
before dividing by sfreq. Donor estimates are converted to target-native samples,
rounding to nearest integer with ties away from zero; residuals are retained.
Digital fallback keeps physical lag blank rather than claiming measured zero.

Coverage uses inclusive epoch endpoints and the recording's exclusive end, with
the currently hashed -0.5:+2.667 s window. Supplied detector/donor uncertainty is
checked separately; it does not bound unknown transfer error. Boundary overlap
and coverage changes are findings for Stage 5, not edits to Stage 2 eligibility.

### Participant reports and report-only rebuilding

The reporting upgrade was approved on 2026-10-21. Open a participant's
`report.html` to review lag consistency across trigger factors and recordings,
the quantity and basis of actual proposed corrections, QC, donors/uncertainty
and timing-dependent coverage. The report is a single portable HTML file with
embedded figures and downloadable tables; no web server or network is needed.
`report.md` and PNG/CSV companions remain convenient alternatives.

All trigger codes are factors. Bar height is measured mean lag; whiskers are
sample SD (n-1 denominator), with medians and n also shown. Violin/strip plots
retain individual measured and flagged edges. Recording/runtime groups are
kept separate on comparable scales. Sparse/constant groups have no invented
density and n=1 has no SD. Measured trials use their own edges; these descriptive
centers do not silently become corrections. Donor location CIs, donor spread,
held-out errors and unknown transfer uncertainty are explained separately.
Counts/coverage bars do not receive meaningless SD whiskers.

The interval additions were approved on 2026-10-21 with the original
`physical_intervals.png` bar figure retained. An additional
`physical_intervals_violin.png` shows the same recording groups/order, every
measured pair, a median and mean +/- sample SD. Sparse/constant groups retain
points without a fabricated density. Three further figures show residual
(observed minus expected, ms) against expected duration, native recording time
and categorical endpoint trigger pair. All finite tails remain on full-range
scales; a central zoom and a new pair-coverage panel were not selected.

These views check sequence consistency, not absolute visual latency. Observed
intervals use the two saved native diode-edge samples. Expected intervals use
the saved motion-frame difference and nominal `Conf.refrate`; measured refresh
and behavioral/PTB timestamps do not replace that clock. Report generation
verifies exact endpoint identities, measured status, same source/trial context
and native-sample/frame arithmetic. Inconsistent evidence fails explicitly.
Several pairs can share one trial/onset: pair SD is descriptive spread, not
independent-trial uncertainty. Trigger definitions, actual worked examples,
per-source/pair statistics and exact largest-residual keys accompany the fully
explained figures. The largest-residual inventory is not a rejection rule.

The lag-distribution x limits include group means +/- two sample SD and every
finite measured/flagged lag, with a small native-grid-aware margin. Zero is not
forced into a positive-only view and no finite tail is clipped. Recordings with
no finite lag appear in a compact notice and remain in the QC/coverage tables.
Each coloured dot is one measured marked event; an open diamond is the median
within that exact recording/runtime/trigger group (middle sorted value, or mean
of the middle two for even n). Category-axis jitter only separates overlapping
dots; it carries no timing/order meaning. HTML figures have a CSS-only full-
resolution view with horizontal scrolling for readable labels in narrow panels.

Every figure includes what its axes/marks measure, source fields/calculation,
the actual participant findings, how they inform a decision and a clearly
labelled example. Native waveform selection retains all pilot windows; full
proposals cache up to 16 representative windows per recording, prioritizing
event-family controls, post-warm-up context, QC cases and tails. The gallery
index makes this selection explicit; it is not a random sample or proof of
whole-recording usability. Complete measurement/QC tables remain available.

New timing jobs create the richer report automatically. Existing immutable
bundles can be presented without remeasuring signals or editing timing maps:

```bash
python meg/code/preprocessing_basic/visual_timing_report.py \
  --bundle /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/visual_timing/s3-pilot-20261020-v003/sub-11/ses-01 \
  --output /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/visual_timing/NEW_REPORT_VERSION/sub-11/ses-01 \
  --raw-root /Volumes/MORWUR/Projects/DAMIANO/Attention/data
```

The report-only command verifies all source output hashes before/after rendering,
records the exact source manifest hash, copies scientific table bytes intact and
uses only retained native caches. It refuses an existing output directory or
outputs outside `derivatives/reports/visual_timing`. It never opens a raw FIF or
calls a detector/estimator. A presentation report cannot substitute for an
accepted timing-map bundle. Both pilot and full proposal reports retain the
independent-review gate and participant-specific decided donor restrictions.

The reporting upgrade passed the full 301-test local suite on 2026-10-21.
The [four cached pilot reports](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/visual_timing/s3-review-20261021-v004/index.html)
are available for P08/P09/P11/P24. Their
[validation record](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/visual_timing/s3-review-20261021-v004/report_validation.md)
documents exact output/source/embedded-table hash checks and all 48 unchanged
native cache windows, with zero new raw-signal reads. These reports retain their
pilot limitations and cannot establish full-cohort timing or donor calibration.

The approved interval additions passed the full 317-test local suite on
2026-10-21. The [updated full-cohort participant reports](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/visual_timing/s3-review-20261021-v005/index.html)
present P04-P39 from the existing `s3-v001` bundles. They retain the bar figure
and add the four selected views with complete explanations and exact interval
tables. Report-only publication preserves scientific-table bytes and source
hashes; it does not remeasure signals, alter timing assignments or approve
scientific choices. The [validation record](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/visual_timing/s3-review-20261021-v005/report_validation.md)
records cohort coverage, output/source checks and the report-only boundary.

## Bounded validation and cluster execution

The fixed participant pilot selects before reading signals: P08 acquisitions
1/2 (6+6); P09 acquisitions 1/2/3 (3+6+3); P11 acquisitions 1-5 (1+1+4+3+3);
P24 acquisition 5 (12). It covers recording positions, condition codes and
available diagnostic families without exceeding 12 windows per participant or
48 participant windows total. P04/P17 and partial-boundary checks are header/
table-only. Historical test windows are a separate small validation set.

```bash
python meg/code/preprocessing_basic/visual_timing.py \
  --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-09.yaml --version s3-pilot-v001 --pilot
```

Pilot bundles are explicitly incomplete and can never feed epoching or a full
cohort summary. Inspect traces and independent edge calculations before using
the candidate method for a full cohort proposal. See the cluster runbook for
user-submitted validation/pilot/cohort/summary commands. Resource requests are
starting bounds; inspect actual pilot elapsed time and scheduler accounting.

Candidate refinements can reuse an intact native pilot cache without reopening
participant signal samples. Keep the source version and create a new version:

```bash
python meg/code/preprocessing_basic/visual_timing.py \
  --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-09.yaml \
  --version s3-pilot-20261020-v002 --pilot \
  --pilot-cache-from s3-pilot-20261020-v001
```

Cache reuse verifies all source-bundle hashes, accepted inputs and native window
coordinates, records the source manifest hash and cannot expand the signal
scope. A matching resume must retain the same cache provenance argument.
The manifest records software versions, actual signal reads/native windows
assessed, elapsed time and peak process memory. Cohort aggregation requires
matching implementation and software versions across full proposal bundles.

### Implementation validation completed on 2026-10-20

The [validation and handoff record](/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/visual_timing/s3-validation-20261020-v002/validation/implementation_validation.md)
links the immutable evidence bundles and preservation check. The full local
test suite passed 285 tests (8 synthetic-FIF filename warnings); all 50 Python
sources compiled in memory and all 18 cluster wrappers passed shell syntax.
All 36 participants passed accepted-input/header verification: 185 recordings,
26,467 original candidates, 21,890 eligible and 4,577 excluded. Stage 1 scientific
inputs and accepted Stage 2 report/output hashes remained unchanged.

The original bounded reads were 48 participant and 14 historical 250-ms windows.
All refinements reused those native caches. Final participant pilot
`s3-pilot-20261020-v003` has 46 usable measured edges and the two confirmed
unusable P11 early-block signals. The independent maximum-slope comparison
agrees within 1 ms for all 46 usable edges, with no digital-anchor/code mismatch.
All 14 available historical comparisons reproduce archived lag measurements
exactly; missing v21/old 100-Hz sources were not remeasured.

The six sampled P09 acquisition-2 transitions were measurable despite attenuation
(about 0.0128-0.0183 V; 18-20 ms candidate lags). This does not settle whole-run
usability or the detachment lead. The bounded P11/P24 donors were insufficient
per condition, so no numerical early-block estimate was produced. Genuine
early/low-tail edges remain visible, including the context-qualified post-warmup
method refinement. These checks validate implementation on a bounded sample;
they do not accept the detector, QC, fallbacks, transfer assumptions or cohort
timing. No cohort jobs, cleaning or epochs were run.

## Reviewed epoching handoff

`04_epoching.py` accepts `--timing-review /absolute/path/review.json` or
`epoching.timing_review_path` in a subject/tuning configuration. P04-P39
production defaults to `require_stage3_timing: true`. This additional setting is
outside the Stage 2 policy projection, so it does not invalidate accepted bundles.
The runner has no Stage 3 registration; timing uses its dedicated wrappers.

After independent reconciliation, a human-authorized review record must identify:

```json
{
  "schema_version": "1.0",
  "status": "accepted",
  "reviewer": "named independent reviewer",
  "reviewed_utc": "actual review timestamp",
  "decision_reference": "PREPROCESSING_DECISIONS.md#actual-accepted-entry",
  "bundle_manifest": "relative/or/absolute/path/to/manifest.json",
  "bundle_manifest_sha256": "actual accepted manifest hash",
  "accepted_assignment_bases": ["measured", "estimated", "digital_fallback"]
}
```

This is a schema example, not an acceptance record. Only include assignment
bases actually reviewed for that participant. No measurement job creates it.
The separate application record is bound to the timing bundle. The review JSON
must also accept its exact application hash and all mandatory review scopes;
`application.json` beside the bundle is discovered automatically unless an
explicit `--eligibility-application` path is supplied. The loader verifies
bundle/output hashes, Stage 2 report identity, complete keys,
coordinates/conditions, unchanged original eligibility plus the separately
accepted final application, shift arithmetic, exact usable
measured edges, decided donor restrictions, event order and collisions. Original
`raw_sample_frame` remains untouched; `epoch_event_sample_raw` is the MNE event.
Epoch metadata/report retain timing review and manifest hashes. Epoch resume
cannot skip an absent/different reviewed timing gate. Cleaning must still pass
its separate lineage/coordinate checks before any epochs can be created.

Changing the epoch window remains an explicit future contract change: neither
the Stage 2 hash policy nor the reviewed timing window is silently relaxed.

<a id="planned-eligibility-and-exclusion-application-not-implemented"></a>

## Eligibility and exclusion application: accepted cohort handoff

The accepted Stage 2 ledger and immutable `s3-v001` bundles remain unchanged.
The parent-controlled review log records both a restricted catch-policy
revision and exact sequence-quality decisions. The separate application,
P24-only donor recomputation and hash-bound consumer are implemented; Task 2
added the cross-participant measured-median onset method and published one
coherent `s3-v002` for all 36 participants. See [the application/API
contract](STAGE3_REVIEW_APPLICATION.md). Original Stage 2 eligibility,
`s3-v001`, and every measured event row remain intact. Final overlays preserve
original flags and add `final_analysis_eligible`. The reconciled policy total is
23,963: 21,890 originally eligible + 2,079 newly eligible type-2 catches - six
newly excluded original trials, after all eight exact sequence decisions. The user accepted the motion-onset and revised eligibility scopes on 2026-10-21.
Thirty-six separate review records bind exact participant manifest/application
hashes and mandatory scopes; participant bundle contents remain immutable.

### Feasibility audit from preserved evidence

The saved behavioral `runResults` MAT outputs preserve `catch_type_code` and
condition/block/run/trial/replay identity. The Stage 3 contexts omit the catch
type; Stage 2 stores event-code evidence, not the saved class. A read-only join
across P04-P39 found 2,191 saved type-2 outputs for motion conditions 81/91.
Five are MAT-only outputs for P37 block 1/run 3 trials 4, 8, 13, 16 and 21, in
the logged unavailable tail of `sub37_run01_051026.fif`; they have no recorded
onset, Stage 2 key or timing map. The other 2,186 join exactly to Stage 2 trials
and `s3-v001` maps. Of those, 2,081 were excluded only by the original
`catch_type_2` rule; 105 retain additional independent exclusion reasons. All
2,186 map rows have a proposed window inside the recording and no trial-end
overlap. Their timing assignments include 1,819 measured onsets, 343 with no
recorded diode event and 24 marked confirmed unusable. These counts establish a
viable evidence path, not revised eligibility, timing acceptance, or the final
cohort denominator.

### Published `s3-v002` reconciliation

The final application retains 2,186 recorded type-2 keys: 2,081 lose only the
old `catch_type_2` reason, while 105 keep independent exclusions. Of the 2,081,
2,079 become newly eligible and two remain excluded by exact sequence review.
Four more sequence decisions exclude previously eligible trials, for six newly
excluded original trials overall. Five P37 MAT-only tail outputs remain
unrecorded and receive no timing assignment or epoch. All 2,079 newly eligible
catch onsets have proposed samples within the recording and no trial-end
overlap: 1,718 measured and 361 estimated; none is digitally substituted.
Proposal-wide checks found no out-of-bounds windows, trial-end overlaps, or
measured-edge mismatches.

The donor table contains 22,120 unique QC-usable measured events, separate
from participant-analysis eligibility. Codes 71/81/91 have medians of 19 ms
and donor counts 7,341 / 7,422 / 7,357. The method is the ordinary median over
unique QC-usable measured cohort events, once per event, without block weighting
or bootstrap. P04-P07/P17 use estimated +19 ms for eligible 71/81/91 onsets
(counts 624/661/660/660/661); transfer error from historical v20 remains
unknown. P24's five reviewed events are absent from the pooled donor table.
Its affected donor counts changed 49→48 for code 81 and 90→89 for each of
111/112/114/115; all five medians remain 19 ms. Codes 116/117 retain six
donors each and no numerical estimate.

P33's sequence exclusion leaves 601 remaining original-eligible onset rows,
all measured. The user accepted those 601 assignments and the 2,079 newly
eligible catch assignments on 2026-10-21; the review records identify the
exact scopes and hashes.

### Revised catch policy and question timing

Type-2 catches in occlusion conditions are eligible in principle for participant
motion analysis using the ordinary condition onset (81/91), subject to all other
exclusions and review. Type-1 disappearance/reappearance catches remain
excluded. Determine the class from saved `catch_type_code`, condition and exact
trial identity. Do not infer class from Stage 2 catch-marker evidence or code
118 alone.

Code 118 marks the later question prompt. Keep its measured, estimated or
unresolved timing as separate, lower-priority evidence for possible behavioral
analysis. Do not use it as motion onset or pool it with motion-onset
photodiode lags. The behavioral endpoint and clock alignment have not been
defined, and no prompt estimate may be invented.

### Exact whole-trial decisions

The following confirmed user decisions apply to participant analysis. Preserve
all recorded measurements. Only P24 changes donor membership.

| Participant | Exact key and source | Participant-analysis decision | Donor input |
| --- | --- | --- | --- |
| P24 | Run 5 / block 5 / behavioral run 2 / trial 51 (replay); code 81 at 582257, `sub24_run05_100926.fif` | Exclude the whole trial | Exclude all five measured events: 81@582257, 111@583362, 112@583496, 114@584004, 115@584138. Recompute estimates for codes 81/111/112/114/115; expected counts 49 to 48 for 81 and 90 to 89 for each remaining code. |
| P15 | Run 1 / block 1 / behavioral run 3 / trial 57, replay=true; code 81 at 1074390, `sub15_run01_030926.fif` | Exclude the whole trial | None; timings are measured. |
| P28 | Run 5 / block 5 / behavioral run 2 / trial 13; code 91 at 352984, `sub28_run05_150926.fif` | Exclude the whole trial | None; timings are measured. |
| P29 | Run 1 / block 1 / behavioral run 3 / trial 10, non-replay; code 81 at 1500707, `sub29_run01_150926.fif` | Exclude for confirmed sequence irregularity, independently of catch type | None; timings remain measured evidence. It was already excluded under the old catch policy. |
| P29 | Same run/block/source; trial 12, non-replay; code 81 at 1519036 | Exclude for confirmed sequence irregularity | None; timings remain measured evidence. |
| P29 | Same run/block/source; trial 13, non-replay; code 81 at 1523905 | Exclude for confirmed sequence irregularity | None; timings remain measured evidence. |
| P33 | Run 4 / block 4 / behavioral run 3 / trial 39, non-replay; code 81 at 811230, `sub33_run04_220926.fif` | Exclude for confirmed sequence irregularity | None; timings remain measured evidence. |
| P38 | Run 2 / block 2 / behavioral run 2 / trial 24, non-replay; code 81 at 622800, `sub38_run02_051026.fif` | Exclude for confirmed sequence irregularity, independently of catch type | None; timings remain measured evidence. It was already excluded under the old catch policy. |

P29's trial 10 onset-to-111 interval is 1209 ms optically versus 941.67 ms
expected; saved display timing is 1208.64 ms. Trial 12's 112-to-114 interval is
776 ms versus 508.33 ms expected, with 775.16 ms in saved display timing. Trial
13's onset-to-111 interval is 1134 ms versus 941.67 ms expected, with 1133.73 ms
in saved display timing. P33's excess is about 9-9.67 ms; the saved frame
80-to-81 update is 17.104 ms before occlusion. P38 trial 24's excess is
26-26.67 ms, matching saved display updates of about 21.59 and 20.05 ms into
frames 86/87. After inspecting both native traces, the user judged these
trial-specific excesses too large and reconfirmed the P33 trial-39 and P38
trial-24 exclusions. Causes and perceptual effects remain unknown. No blanket
one-refresh cutoff is established. Keep local diode onset lags separate from
these within-trial sequence checks.

P29 also retains all 37 currently eligible always-visible trials in
`sub29_run02_150926_halfBlock.fif`, each using its own measured code-71 onset
(18-20 ms; median 19 ms). Keep six type-1 catches and the interrupted final
trial 44 (`71@269558`, `152@270938`) excluded. Preserve its source and block
identity separately from the other acquisitions, including block 6; this
retention accepts no averaging or weighting rule. P29's other retained onset
assignments remain event-wise measured, without median replacement.

P38's short first recording is handled as an ordinary block with its own source
identity. In `sub38_run01_051026.fif`, 16 of 20 code-71 onset rows are eligible;
all 16 are measured and their proposed windows are inside the recording. The
other three have gaze-break exclusions and one has gaze-break, ESC and reviewed
interruption exclusions. Keep the smaller eligible count visible; apply no
special timing method or weighting. Retain this source separately from block 6.
P38's remaining measured onset assignments are accepted event by event. Its
13 ms code-71 edge is specific to acquisition 5/block 5/behavioral run 1/trial
2 (`71@68884`, diode sample 68897); keep the other accepted events at their own
measured lags. No 13 ms group shift is accepted. Native traces for P33 trial
39 and P38 trial 24 were presented, and the user reconfirmed both whole-trial
exclusions after inspection. P33's remaining 601 original-eligible measured onset assignments were accepted
with the cohort on 2026-10-21. The case decisions are included in the formal
cohort acceptance and hash-bound review records.

### Accepted cohort method, scope and remaining limitations

The accepted `s3-v002` donor table contains 22,120 unique QC-usable measured
events. The estimator is the ordinary median over unique events, without block
weighting or bootstrap. Codes 71/81/91 have n=7,341/7,422/7,357, each with a
19 ms median. P04-P07/P17 use estimated +19 ms assignments for their eligible
motion onsets (624/661/660/660/661 rows). P11's early blocks use its own usable
later measured blocks; P24's unmeasured runs use only its run-5 donors.
Transfer error remains unknown, particularly across P04-P07's v20 runtime.
Donor spread does not quantify target transfer error. The requested v20 optical
bench measurement is still future work and cannot calibrate historical
participants.

The accepted application yields 23,963 eligible trials (21,890 original +
2,079 newly eligible type-2 catches - six newly excluded original trials). Of
2,186 recorded type-2 keys, 105 retain independent exclusions; five additional
P37 MAT outputs in the unavailable tail have no recorded onset or assignment.
New catch assignments are 1,718 measured and 361 estimated. All proposed
windows are recording-contained without trial-end overlap. P33 retains 601
original-eligible measured onsets after its exact trial-39 exclusion. The user
accepted the new catch and P33 scopes on 2026-10-21.

P24's five events remain measured evidence but are removed from pooled donor
inputs. Its affected donor counts are 48 for code 81 and 89 for each of
111/112/114/115; all affected medians remain 19 ms. Codes 116/117 still have
six donors each and no numerical estimate. P29 retains 37 eligible measured
onsets in its partial half-block. P38 retains 16 eligible measured onsets in
its short first recording; its accepted 13 ms onset remains event-specific.
Code 118 remains deferred, lower-priority behavioral timing; its endpoint and
clock are undefined and it is never substituted for motion onset.

The user accepted the exact `s3-v002` timing/application scope on 2026-10-21.
Separate review records in `s3-v002/reviews/` bind all 36 participant manifests,
application hashes, eligible assignment bases and mandatory scopes. Original
Stage 2 and `s3-v001` evidence remain unchanged. The real `load_reviewed_timing`
consumer passed P04-P10 and targeted P11/P17/P24/P29/P33/P38 checks; the prior
one-pass integrity reconciliation covered all 36 bundles. No cleaning,
epoching, step-06 export or analysis was run. Stage 4 cleaning remains
independent; Stage 5 may use the accepted records only after accepted cleaning
lineage and ordinary production prerequisites are ready.

## Completed review and supported Stage 5 gate

The independent review reconciled original keys and unavailable boundaries,
reviewed raw trace usability and sequence tails, checked P09 signal usability,
P11/P24 donor restrictions, donor provenance and transfer limits, and verified
proposed coverage. Accepted records bind exact bundle/application hashes and
scope. The Stage 5 consumer still rechecks Stage 2 lineage, keys, arithmetic,
measured edges, application hashes, assignment bases and accepted scopes at
runtime. Cleaning lineage and epoch configuration remain separate prerequisites.
