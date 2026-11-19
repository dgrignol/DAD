# Stage 3 review application: implementation and accepted s3-v002 handoff

The Task-1 application and Task-2 onset method are implemented and published
as one coherent `s3-v002` for all 36 participants. Scientific acceptance was
recorded separately on 2026-10-21 in the formal decisions file and 36
hash-bound review records. The implementation flow below documents how the
version was produced; do not rerun it to recreate the accepted application.

## Run and output boundary

```bash
PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python3 -B \
  meg/code/preprocessing_basic/visual_timing_application.py \
  --config meg/params/pipeline_default.yaml --source-version s3-v001 \
  --target-version s3-v002 --workers 4 --validate-only

PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python3 -B \
  meg/code/preprocessing_basic/visual_timing_application.py \
  --config meg/params/pipeline_default.yaml --source-version s3-v001 \
  --target-version s3-v002 --workers 4
```

The default new area is
`derivatives/visual_timing/s3-v002-task1-staging`. Existing destinations are
refused. The CLI reads saved MAT metadata and existing measured-event tables;
it does not read new waveform samples, run the detector, refresh a manifest,
submit jobs, or run cleaning, epoching, exports or analysis. P24's full report
is rebuilt from existing caches and filtered estimates; the other Task-1
reports provide application counts/evidence for the final Task-2 report pass.

## Exact policy and independent donor membership

`visual_timing_application.policy_definition()` records the one supported
versioned policy and all eight exact sequence decisions from the working log:
P24, P15, P28, P29 trials 10/12/13, P33 and P38. It does not introduce a generic
sequence-error cutoff. Every original trial key/source/sample remains present.

The saved MAT `catch_type_code=2` qualifies only conditions 81/91, joined by
block, behavioral run, trial, replay and condition, with saved debug CSV identity
checked across Windows/POSIX path spellings. Only `catch_type_2` is waived;
gaze, interruption, type-1, unavailable-tail and other independent reasons
remain. Sequence exclusions override catch inclusion. A question marker alone
does not establish catch class. Empty planned-run MAT placeholders do not
create trials. Accepted interrupted onsets absent from completed MAT output
retain their original disposition with explicit class-unavailable evidence.

Only P24 changes donors: remove **all five** measured keys from
`sub24_run05_100926.fif`: `81@582257`, `111@583362`, `112@583496`,
`114@584004`, `115@584138`. Their measurement rows remain unchanged. Recompute
dependent P24 estimates/maps/event proposals. Expected counts are 49 to 48 for
81, and 90 to 89 for each of 111/112/114/115. P24 116/117 retain six donors each,
an insufficient-donor status and no numerical estimate. P11's own blocks 3+
restriction and estimates are unchanged. **Never filter any donor pool by final
analysis eligibility**; the other seven sequence decisions are analysis only.

P29's short source keeps its 37 eligible always-visible onsets, six type-1
exclusions and exact interrupted final onset. P38's short source keeps 16 of
20 onsets and all four independent exclusions. Source/block identities,
including block 6, remain separate. Five P37 MAT-only type-2 outputs retain no
recorded onset or trial key.

## Immutable fields and staged products

| Field/product | Meaning |
|---|---|
| `included_by_trial_policy`, `exclusion_reasons` | Original Stage 2 policy values, preserved in timing maps and later epoch metadata. |
| `final_analysis_eligible`, `analysis_exclusion_reasons` | Computational selection after the application; scientific acceptance is bound separately in each participant review JSON. |
| `eligibility_overlay.csv` | One row per original key, exact identity, saved class, original/final eligibility, waived/retained reasons and sequence decision. |
| `saved_trial_metadata.csv` | Frozen saved class/context inventory, with MAT filename/hash and no fabricated signal times. |
| `application.json` | Policy, source/Stage-2/MAT/log hashes, overlay hashes, independent donor exclusions, scientific review requirements and before/after summary. |
| `review_notes_snapshot.md` | Frozen input snapshot; later authorized edits to live `s3-v001/review_notes.md` do not invalidate this application. |
| `new_catch_timing_qc.csv` | Newly included onsets with measured/estimated status, flags, edge span and coverage; acceptance is recorded separately. |
| `analysis_counts.csv` | Original/final counts and coverage by source, condition, identity, replay, block and joint condition/identity/source factors. |
| `unrecorded_saved_type2.csv` | Accepted unavailable-tail saved outputs, explicitly without an onset/key. |
| `donor_exclusions.csv`, `p24_estimate_changes.csv` | Exact donor removals and P24 numerical/count sensitivity; not an implicit analysis filter. |
| `task1_handoff.json`, `participant_counts.csv` | Whole-cohort checkpoint, per-participant manifest hashes, counts, API and final-version reservation. |

The source working log is input only. Its frozen content/hash is retained;
the live log and formal decisions record are not edited by Task 1. Existing
Stage 2 verification remains intact. Later formal acceptance updates are not
treated as edits to accepted Stage 2 evidence.

## Completed Task 2 handoff

```python
from visual_timing_application import load_task1_stage
manifest, contents = load_task1_stage(staging_root, "04")
```

The loader requires the complete sealed 36-participant checkpoint and verifies
the participant manifest and all staged output hashes. `contents` maps filenames
to immutable input bytes. Do not edit staging in place. Keep `application.json`,
its hash, metadata/overlay, source snapshots and original policy columns intact.
Read timing maps as CSV strings and preserve original coordinate spellings;
write newly computed shifts/proposed samples as integers on native coordinates.

Task 2 completed these steps before publication:

1. Apply only its supported onset method for P04-P07/P17. Read genuinely measured
   donor events from the original source catalog; use a separate donor artifact,
   not appended or relabelled `event_measurements.csv` rows. Exclude only the five
   P24 donor keys, independently of analysis eligibility. Freeze exact member
   keys, source hashes, event weighting, code-specific medians and uncertainty.
2. Extend the estimated-assignment validator explicitly for that method and its
   five targets. The final loader validates the accepted method while retaining
   event-specific labels and unknown transfer error.
3. Rebuild maps/proposals and all final participant reports/QC/coverage using
   `final_analysis_eligible`; retain original eligibility/reasons separately.
   Keep codes 118 and other families separate, and no P24 116/117 numerical value.
4. Publish all 36 to the fresh final `s3-v002` using `visual_timing.publish`.
   Final mode is `cohort_proposal`, `task1_staging=False`. Preserve
   `review_application_sha256`, original Stage-2/source hashes, and extend the
   coherent input/method/implementation contract for Task 2. The staging source
   contract is a parent reference, not a reason to retain a stale final contract.

## Supported consumer and acceptance gate

`04_epoching.py` uses `analysis_eligible()` for candidate events, MNE outcome
accounting and final count/reason reports. Original Stage 2 flags and reasons
remain in metadata. It accepts optional `--eligibility-application` naming the
same `application.json` already bound to the reviewed bundle; it cannot attach
an arbitrary override or bypass review. Resume rechecks the complete gate.

In addition to the existing accepted timing record fields, a revised bundle
requires **human-authorized** `accepted_eligibility_application_sha256` matching
its application and `accepted_review_scopes` containing every requirement in
`application.json`. Requirements include final eligibility, newly included
catch onsets where applicable, and P33's remaining original onset scope.
Staging/pilot modes, changed original coordinates/eligibility, missing keys,
changed MAT class evidence, modified overlays, excluded donor members,
unaccepted bases/scopes, invalid shifts, collisions and order changes fail.

No actual participant epoching or downstream analysis has been run for this
implementation. Scientific reconciliation accepted the new catches, P33, the
Task-2 donor/method results and reviewed coverage on 2026-10-21. Transfer error
remains unknown. Stage 4 cleaning and downstream condition/identity balance
checks remain separate prerequisites.

### Task-2 release hardening completed

`load_applied_eligibility()` derives mandatory scopes from the recomputed
overlay and participant rather than trusting only declared
`scientific_review_requirements`: `final_analysis_eligibility` always,
`newly_included_catch_onsets` when catches are added, and
`p33_remaining_original_onsets` for P33. The focused consumer suite passed 48
tests before the accepted release. The supported loader also passed the
read-only accepted-record checks documented in
`s3-v002/independent_reconciliation_summary.json`.
