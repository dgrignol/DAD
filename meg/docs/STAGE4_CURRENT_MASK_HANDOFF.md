# Current-mask adoption and final ICA handoff

For participant-level final component decisions use the [final ICA review guide](FINAL_ICA_COMPONENT_REVIEW_GUIDE.md). The component book is the human interface; the standard HTML records processing evidence.

Updated 2026-11-11; implements Damiano's explicit 2026-11-07 decision.
The authoritative scientific record remains
[PREPROCESSING_DECISIONS.md](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md).

## Current scope and evidence

The working cohort is **P04-P14 and P16-P39: 35 participants, 180 acquisitions**.
P15 is excluded from current analyses for its insufficient remaining coverage:
76 always-visible, 75 occluded-nondeviant, 79 occluded-deviant annotation survivors,
with ten empty condition-by-identity cells. This is a participant-specific
decision, not a general numerical exclusion rule. Preserve its data, all earlier
36-participant snapshots and all Stage 1-3 evidence for possible recovery.

The input is the immutable `annotation_triage/reports/retention-20261107-v1`
snapshot (36 original participants, 185 acquisitions, zero input errors) and its
exact per-run mask links. Its embedded-figure companion remains
`annotation_triage/visual_reports/retention-20261107-v1/index.html`.
Current saved manual work is discovered dynamically: at preparation, P04/P05 all
runs, P06 runs 1-3 and P26 all five runs contain 18 genuine saved reviews.
The 35-participant pass therefore has 18 preserved-manual and 162 automatic runs.

Further exhaustive manual interval inspection is deferred until later findings
justify it. Adopt current automatic ocular/muscle masks and existing source BADs
for processing; never represent this as manual SAVE. Manual files and exact-source
receipts retain precedence, including edited/deleted ocular proposals, complete
muscle replacements, reviewed-empty muscle and deletion veto provenance. Legacy
ocular-only reviews retain their existing automatic muscle category. A saved file
with a missing/stale receipt fails validation; it is never replaced automatically.

Read-only preparation completed on 2026-11-09: all 35 participants/180 runs
validated, with 18 preserved manual reviews and 162 automatic masks. There are
19,991 conditional annotation survivors in the working cohort before FIF storage;
every participant's timing/application and survivor counts match the Nov 7
snapshot. No manual save superseded that snapshot. The immutable plan is
`stage4b_adoption/current-masks-20261109-v1/plan.json`, with content identity
`2074b2b80f24bcce90fd8cccf4dee843011cc59fa32e786c94dfaacbd40cb231`.
The adjacent `preparation_verification.json` records input-only status. Production
adoption and joined policy acceptance subsequently completed. The saved
`stage4b_adoption/current-masks-20261109-v1/policy_acceptance.json`, generated
2026-11-09 at 11:46 UTC, records all 35 participants as
`accepted_for_final_ica_fit` under this exact plan. It explicitly leaves final
ICA decisions and cleaned-output acceptance pending. Its operative total is
19,988 annotation survivors. The three differences from the pre-storage total
are recorded `BAD_blink_auto` endpoint-overlap changes in P36 runs 2, 3 and 4;
the receipt identifies each exact trial, with no new padding/threshold rule.
The 35 final prefit models
now exist; use the final component-review guide for their current diagnostics
and the documented [P19 numerical reconciliation](P19_ICA_NUMERICAL_RECONCILIATION.md). Require the completed publication and joined current readiness snapshot before review. The commands below document
the completed adoption path, not a request to repeat it.

## Implementation and audit boundary

`stage4b_adoption.py plan` validates the entire working cohort in independent
processes before publishing one immutable portable plan. It rechecks accepted 4A,
current SSS/native coordinates, muscle basis, manual files/receipts, automatic
component selections and cached source/model identities, detector settings and
accepted Stage 2/3 eligibility/timing/application evidence. It reruns no detector.
Newer valid manual saves supersede snapshot masks during planning and are recorded;
changes after planning require a new plan. Missing/deferred selections are errors.
Unscored detector periods remain recorded limitations, not evidence of no artifacts.

`apply` runs independently per participant. It renames original annotation FIF/CSV
outputs into `preprocessing/sub-XX/ses-01/history/mask_adoption/<plan-sha256>/`,
preserving their modification times and avoiding a second large copy. Prior report
JSON/HTML and tuning are copied there. The source SSS, manual files, old triage
snapshots and calibration/model resources remain unchanged. Retry uses the same
archived basis and accepts only outputs bearing the exact planned mask token.
The stage retains existing muscle score PNG/NPZ diagnostics without recomputing
or interpreting them as the operative manual muscle set.

Standard Stage 4B outputs are rebuilt as `_ann.fif`, `_ann.csv`,
`03_B_ArtefactAnnotation_review_report.json/html` and source-manifest lineage.
The tuning YAML gains a separate `artifact_annotation.mask_adoption` pointer.
Per-run `manual_eye_review` remains genuine reviewed or explicitly unreviewed;
`annotation_adoption` separately records method, policy and mask SHA-256.
Original BAD/non-BAD boundaries and source metadata remain protected, except the
explicitly editable muscle category under a valid reviewed replacement.
Native sample coordinates, channels/bads, STI101 and MISC are checked after saving.
Per-trial eligibility, incomplete windows, all overlapping BAD causes and survivors
are checked against an independent public-MNE annotation-only FIF round-trip; equal
totals alone are insufficient. Endpoint storage precision can change strict
overlap exactly at an epoch boundary. Original proposal masks remain intact;
`serialization_outcome_changes` records every such before/operative trial delta.
No new interval padding or numerical exclusion criterion is introduced. Refreshed
operative condition, identity and run-condition retention CSVs accompany the
joined cohort receipt, which is published last.

`reconcile` reopens all participant outputs, checks CSV/FIF and frozen source,
model, manual and Stage 2/3 bindings, then joins every result. Only after this
whole-cohort validation barrier can `--accept-policy` record the user's already
explicit **annotation-policy acceptance for processing** in `stage4_review`.
This is distinct from exhaustive manual review and clean-data acceptance. The
single immutable `policy_acceptance.json` is published after participant receipts.
Downstream consumers recheck the adoption bindings and require the joined
cohort publication through `require_review`.
Successful array tasks alone cannot satisfy the gate.

## Exact commands

No scheduler jobs are submitted by preparation. Run deployment from the local DAD
workspace after reviewing this handoff:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
meg/cluster/deploy_review_to_cluster.sh
```

The shared derivatives plan is `stage4b_adoption/current-masks-20261109-v1/plan.json`.
Its validation status is recorded in the adjacent README. If preparation did not
complete, or any input/manual save changed, use a **new label** and submit ONE
whole-cohort planning job first (the example below uses a new cluster label):

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
export DAD_ADOPTION_LABEL=current-masks-20261109-cluster-v1
DAD_ADOPTION_MODE=plan sbatch --cpus-per-task=35 --mem=64G \
  meg/cluster/submit_stage4b_adoption.sh
```

After a validated plan exists, use its actual label. These commands submit the
whole working cohort, then one dependent joined reconciliation, then final fitting
only after policy processing acceptance. No further interval-review confirmation
is required. The shell variables are ordinary job identifiers; inspect job logs
if any dependency remains unsatisfied.

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
export DAD_ADOPTION_LABEL=current-masks-20261109-v1
adopt_job=$(DAD_ADOPTION_MODE=apply sbatch --parsable --array=4-14,16-39 \
  meg/cluster/submit_stage4b_adoption.sh)
adopt_job=${adopt_job%%;*}
audit_job=$(DAD_ADOPTION_MODE=reconcile DAD_ADOPTION_ACCEPT=1 \
  sbatch --parsable --dependency="afterok:${adopt_job}" \
  --cpus-per-task=35 --mem=64G meg/cluster/submit_stage4b_adoption.sh)
audit_job=${audit_job%%;*}
DAD_STAGES=ica_manual DAD_PROCESSING_PASS=review DAD_REBUILD=0 \
  sbatch --dependency="afterok:${audit_job}" --array=4-14,16-39 \
  meg/cluster/submit_preprocessing_review.sh
```

Each adoption task requests 1 CPU/8 GB; final fitting defaults to 4 CPUs/24 GB per
participant. The scheduler bounds concurrent tasks by available allocations.
Use `--cpus-per-task`, `--mem` and, if required by actual available resources,
an array concurrency suffix such as `4-14,16-39%8`. All 35 remain submitted together.
Cohort plan/reconciliation exposes `DAD_PARTICIPANT_WORKERS` (default allocated
CPUs); child numerical-library threads are one. Genuine shared publication and
the whole-cohort validation barrier are ordered. Run progress is printed to
`cluster_logs/DAD_stage4b_adopt.*` and `DAD_meg_review.*`.

Reusing an existing plan label for planning or reconciliation fails intentionally;
apply resumes only its own compatible outputs. A failed input audit requires
reconciling the actual changed source/review and a new plan, not disabling checks.
Do not rerun ordinary automatic 4B detection over an adopted mask plan.

## Final fit outputs and next scientific decision

Accepted 4A SSS is reused unchanged. This pass generates no EDF guidance or SQUID
proposals. Final continuous processing remains spectrum-fit 50/100 Hz removal
(10 s window), including P38 run 1, with **no additional high-pass or low-pass**.
The original irreversible acquisition high-pass constraints remain.

Final ICA prepares each admissible native stretch independently at 1-40 Hz,
trims full FIR support, resamples separate fit copies to 200 Hz and concatenates
afterward. Every operative BAD category, including routine ocular intervals, is
excluded. FastICA uses `n_components=None`, seed 97 and `max_iter='auto'`.
No preliminary per-run annotation-aid index becomes a final-model exclusion.

Per participant, inspect:

- `preprocessing/sub-XX/ses-01/ica_models/sub-XX_ses-01_ica_prefit.fif` and metadata:
  input fingerprint, actual component count, PCA variance, fit samples/maps,
  iterations and whether the iteration limit was reached.
- `reports/sub-XX/ses-01/03_C_ICA_review_components.html`: one self-contained
  current-model book with component selection and keyboard navigation, maps,
  properties, source traces and per-run native frontal/MISC context.
- `03_C_ICA_review_report.json/html`: standard stage evidence and native fit
  interval/sample accounting. Incomplete current-model properties cannot be
  accepted using stale figures from an older decomposition.

The **next review decision** is the exclusion list for each exact final model,
including an explicitly reviewed empty list where justified. Reconcile component
evidence and convergence before recording it with `review_cleaning.py --stage 4C
--exclude ...` (no example indices are decisions). Do not apply an unreviewed model.

After explicit final component decisions, the existing apply wrapper operates on
the original native runs. It produces `_ica.fif`, a post-rejection model, matched
ICA-only comparisons, final line-removal PSDs and native preservation evidence.
P10 stays in the processing cohort; its approximately 0.33 Hz concern remains
open. Its application QC additionally compares up to three matched native 60 s
windows and 0.1-2 Hz spectra per acquisition, before line removal and without
high-pass/detrending. Assess residual contamination and task-signal preservation;
fit-only 1-40 Hz diagnostics cannot settle that question.

Application QC, P10's later retention decision and **cleaned-output scientific
acceptance remain explicit later work**. No epoching command is part of this
handoff. Stage 5 remains blocked until accepted final component decisions,
successful reviewed application/preservation evidence, cleaned-output acceptance
and the already accepted Stage 2/3 eligibility/timing prerequisites are reconciled.
