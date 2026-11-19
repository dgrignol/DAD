# Continuous preparation during epoching - 2026-11-12

**2026-11-17 update:** production job 11699 completed P37 and stopped elsewhere.
The approved short-fragment guard is now implemented. Use the current
[retry command](SHORT_FRAGMENT_RETRY_20261117.md), with a new run label and P10
pending its separate all-rejected-recording diagnostic. The commands and local
preflight status below describe the earlier v1 handoff.

The approved current route is accepted native SSS/annotated recordings ->
continuous spectrum-fit line removal in memory -> epochs -> neural export.
**ICA is disabled for this run.** It does not change saved component choices.
No new continuous `_ica.fif` is written. Future runs can enable the exact saved
ICA model/choices before the same line removal, with a different output label.
The user submits every cluster job.

## Next executable step

From the cluster DAD root, with the prepared runtime and existing `dad-meg`
environment (MNE 1.10.2), export this pair in the shell used for all submissions:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
export DAD_EPOCH_PREPARATION=inline_no_ica
export DAD_PROCESSING_RUN=noica-20261112-v1
DAD_VALIDATE_ONLY=1 sbatch --array=4-14,16-39 --cpus-per-task=1 --mem=24G \
  --time=1-00:00:00 meg/cluster/submit_epoch_export.sh
```

This first array validates all 35 participants independently, including exact
accepted 4A/4B sources, saved manual annotations/deletions, Stage 2 ledger and
accepted `s3-v002` Stage 3 onset/eligibility bindings. It reads headers and small
evidence files; it produces no scientific derivatives. Check that every array
task succeeds. Missing ICA receipts do not block this explicitly no-ICA route.

After that array succeeds, create epochs and exports:

```bash
DAD_VALIDATE_ONLY=0 sbatch --array=4-14,16-39 --cpus-per-task=1 --mem=24G \
  --time=1-00:00:00 meg/cluster/submit_epoch_export.sh
```

Both arrays include the full eligible cohort without a participant throttle;
SLURM schedules them within available resources. One participant per process,
one numerical thread at this allocation. Acquisitions remain ordered within a
participant to bound memory. Export starts only after that participant's guarded
epoch stage succeeds. P15 is excluded; P32 remains included. The 35 participants
have 180 acquisitions. No SSS, annotation detection, timing estimation, ICA fit,
model save or continuous-output stage is rerun.

## Settings and space

- Continuous: native sampling, channels, transforms, acquisition boundaries and
  annotations preserved; all STI101 and MISC samples checked. Spectrum-fit at
  50/100 Hz, 10 s window, independently within good stretches bounded by BAD/EDGE
  annotations. P38 run 1 receives this operation too. No additional high-pass,
  final low-pass or resampling. Existing acquisition high-pass limits remain.
- Epochs: accepted corrected-onset anchors; -0.5 to +2.667 s, no baseline,
  linear detrend (`1`), BAD rejection, grad `5e-10` and mag `4e-12` thresholds.
  Accepted Stage 2/3 eligibility and trial identities stay fixed. Final-filter
  edge exposure remains informational under the existing exact-label rule.
- Export: existing 0/+2.667 s window at 120 Hz, identities 1-40 within each
  condition; fixed canonical MEG order, bad/missing channels as NaN. No extra
  trial selection. Missing required identity coverage remains a visible failure.
- MORWUR currently has about 233.7 GiB free. Avoiding continuous copies saves
  approximately 182.5 GiB. Allow up to approximately 76.8 GiB for native epochs,
  plus exports, analysis and small QC. Retain earlier inputs and outputs.

The continuous operation occurs before cutting epochs, so each short trial is
not filtered independently. In-memory filtering reproduces the former continuous
operation; omitting an intermediate float32 FIF save can change only numerical
rounding relative to that saved-file reference.

## Output and restart contract

All downstream outputs for this selection are under:

```text
derivatives/processing_runs/noica-20261112-v1/
  analysis/sub-XX/ses-01/meg/       epochs, metadata, all trial decisions
  analysis/sub-XX/ses-01/drsa/      exports and subject MATLAB outputs
  reports/sub-XX/ses-01/           derived manifest, epoch/export/MATLAB reports
  reports/sub-XX/ses-01/epoch_input_qc/  per-acquisition line spectra; P10 slow QC
  sanity_checks/.../dot_erf/       subject and group ERF outputs
  reports/group/.../dot_erf/       ERF group review
```

Canonical upstream manifests, source recordings, models, masks, tuning and timing
are read-only. A derived manifest records canonical ancestry and actual epoch
preparation; it never claims a completed saved-ICA application stage. Contracts
bind native split-file stamps/headers/annotations, settings, timing, trigger CSV,
continuous/epoch implementation and actual output FIF parts/small sidecars.
Default resume validates the completed contract and outputs before reuse.

A label identifies one immutable epoch variant. Changing ICA mode, model,
exclusions, source/masks, timing or epoch settings requires a NEW label.
`DAD_REBUILD=1` does not overwrite labelled epoch variants. After an incomplete
labelled write, preserve that evidence and use a new label. Exports and MATLAB
also check their input/settings/output bindings; ERF resume checks cohort epoch
contracts and analysis settings. Keep both environment variables set identically
when running later stages. The legacy `saved_ica` route remains the default when
no explicit inline mode/run selection is provided.

## Join and inspect before analyses

After all epoch/export tasks finish, run the cohort audit (separate participant
processes, one numerical thread each; joined summary only after all finish):

```bash
DAD_QC_PHASE=epochs DAD_QC_TASK_SIGNAL=0 DAD_QC_LABEL=noica-epoch-audit-20261112-v1 \
  sbatch --cpus-per-task=8 --mem=24G --time=1-00:00:00 \
  meg/cluster/submit_ica_application_qc.sh
```

The existing audit wrapper's name is historical: `phase=epochs` checks the chosen
run's epochs/exports and does not apply ICA. It reconciles exact final trial keys,
source/onset/condition/identity fields, every extra MNE drop and per-cell counts
against the accepted annotation-survivor footprint. It reads actual export MATs
and checks channel order, identity/time axes, condition counts and NaN masks.
Audit outputs use a new immutable review label under canonical
`reports/ica_review_updates/`; scientific acceptance is not inferred.

Inspect the epoch report's line spectra and boundary/short-stretch summaries.
P10 has additional native long-window traces and sub-Hz PSDs; retention remains
open. With ICA disabled, remaining ocular/cardiac activity and extra amplitude
rejections are expected QC questions, not reasons to silently change thresholds,
annotations, identities or participant membership.

After the cohort audit succeeds, prepare the existing ERF and real dRSA analyses:

```bash
export DAD_DOT_ERF_GROUP_SUBJECTS=04,05,06,07,08,09,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39
DAD_N_PERMUTATIONS=10000 DAD_RANDOM_SEED=42 DAD_ALLOW_MISSING=0 \
  sbatch --array=1-4 --cpus-per-task=4 --mem=16G --time=2-00:00:00 \
  meg/cluster/submit_dot_erf_analysis.sh
DAD_DRSA_ANALYSES=real DAD_DRSA_SENSOR_GROUPS=wholeBrain DAD_DRSA_MODEL_SOURCES=observed \
  sbatch --array=4-14,16-39 --cpus-per-task=4 --mem=24G --time=1-00:00:00 \
  meg/cluster/submit_drsa_analysis.sh
```

ERF retains existing mag sensors, occipital/parietal/frontal/temporal groups,
-0.300/-0.100 s baseline, -0.5/+2.667 s analysis window and 200 Hz analysis rate.
Only the user-approved permutation choice changes: 10,000 with seed 42 (rather
than exhaustive sign enumeration). Existing ERF interpretation cautions remain.
The MATLAB stage retains its existing analysis/baseline/model settings; MATLAB
threads are bounded by allocated CPUs. Group dRSA publication follows completion
and joining of the selected subject runs. Keep the same run environment, then:

```bash
DAD_DRSA_GROUP_SUBJECTS="$DAD_DOT_ERF_GROUP_SUBJECTS" \
  sbatch --cpus-per-task=2 --mem=12G --time=06:00:00 \
  meg/cluster/submit_drsa_group_average.sh
```

Do not mix processing-run directories; verify all 35 subjects are represented.

## Later ICA-enabled run

Before using upstream ICA review/recording tools, unset `DAD_EPOCH_PREPARATION`
and `DAD_PROCESSING_RUN` so those tools address canonical ICA reports.
After returning to the ICA stage review, finish the exact-model guarded Stage 4C
recording for any still unrecorded participants. See the separate
[ICA review status](ICA_APPLICATION_HANDOFF_20261112.md). Do not reapprove the
already authorized lists; genuinely new choices require their own authorization.
Then use the same epoch/export workflow with:

```bash
export DAD_EPOCH_PREPARATION=inline_ica
export DAD_PROCESSING_RUN=ica-reviewed-NEW-LABEL
DAD_VALIDATE_ONLY=1 sbatch --array=4-14,16-39 --cpus-per-task=1 --mem=24G \
  meg/cluster/submit_epoch_export.sh
# After all validation tasks succeed:
DAD_VALIDATE_ONLY=0 sbatch --array=4-14,16-39 --cpus-per-task=1 --mem=24G \
  meg/cluster/submit_epoch_export.sh
```

This requires the exact accepted model-bound receipt, current model/metadata,
source/cache fingerprint and explicit valid exclusions; no refit fallback exists.
It applies selected ICA spatial rejection on the native continuous recording,
then the same boundary-respecting line removal, then epochs. ICA-only before/after
figures are added to the scoped QC. No permanent continuous ICA copy is needed.
Both no-ICA and ICA runs remain available for comparison.

## Implementation checks

The model-compatible local environment passed 210 targeted checks, including
real synthetic spectrum-fit and saved-model application, operation order,
BAD/EDGE isolation, native coordinates/auxiliary preservation, trial identities,
float32-FIF reference parity, export MAT counts/axes and changed-contract guards.
Eight participants (P04-P11) passed read-only native-input/Stage 2/3 validation on
the Mac shared mount. The local full-cohort attempt was then stopped because
SMB header reads were very slow; **27 participants remain unvalidated by that
attempt**. The cluster preflight above is the full-cohort execution gate.
Runtime deployment/previous bytes and test logs are preserved in
`reports/ica_review_updates/inline-epoch-preparation-20261112-v1/`; partial local
preflight evidence is under `inline-epoch-input-preflight-20261112-v1/`.
No scientific cohort derivatives were produced and no jobs were submitted.
