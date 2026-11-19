# Epoch/export retry after the short-fragment fix - 2026-11-17

The approved fix leaves intervals of eight native samples or fewer untouched by
line-noise removal. Other intervals receive the same continuous 50/100 Hz,
10-second spectrum-fit operation. Empty intervals are also bypassed. Samples,
annotations, timing, epoch selection and ICA choices are preserved. ICA remains
off for the current run; future selected ICA still runs before line removal.
No new continuous FIF files are saved.

The complete check covered 35 participants and 180 acquisitions: 233 fragments,
1-8 samples each, totalling 0.894 seconds. None intersects any of the 19,988 full
epoch windows surviving the accepted annotations. The check used the cluster's
MNE 1.10.2, NumPy 1.26.4 and SciPy 1.13.1, and verified the exact current native
annotation hashes, sample coordinates and accepted Stage 2/3 epoch footprint.
All fragment decisions are automatic; individual fragment inspection is not
required. Each acquisition's epoch report records the skipped sample ranges.

## Submit the retry

Runtime files are copied to the shared cluster project; another deployment is
not needed. The user submits the jobs. From the cluster DAD root:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
export DAD_EPOCH_PREPARATION=inline_no_ica
export DAD_PROCESSING_RUN=noica-20261117-v2
DAD_VALIDATE_ONLY=0 DAD_REBUILD=0 sbatch --array=4-9,11-14,16-39 \
  --cpus-per-task=1 --mem=24G --time=1-00:00:00 \
  meg/cluster/submit_epoch_export.sh
```

This submits all 34 ready participants in parallel with disjoint outputs and one
numerical thread per task. The existing wrapper validates accepted inputs and
timing before processing, then exports only after epoching succeeds. P37 reruns
in this version so the new batch has one consistent processing namespace; its
completed v1 epochs and three MAT files remain preserved. The earlier failed
attempts and their QC also remain in `noica-20261112-v1`.

**P10 remains in the 35-participant cohort but is pending a separate rejection
diagnostic.** Its first acquisition lost every candidate epoch under the existing
rejection settings. The short-fragment guard does not change that result, alter
thresholds or exclude the participant. P15 remains excluded by the prior cohort
decision. Group ERF/dRSA awaits resolution of the pending P10 coverage.

## Expected output and next check

Successful outputs are under:

```text
derivatives/processing_runs/noica-20261117-v2/
  analysis/sub-XX/ses-01/meg/   native epochs, metadata, complete trial decisions
  analysis/sub-XX/ses-01/drsa/  three identity-averaged neural MAT/JSON exports
  reports/sub-XX/ses-01/       epoch/export reports, derived manifest and small QC
```

This batch should yield 34 participant epoch sets and 102 neural MAT files.
Existing epoch, baseline, detrend, rejection, export, identity and channel-axis
settings are unchanged. Successful filtering does not establish final amplitude
rejection or complete identity coverage; the existing automated checks remain.
After completion, check all participant reports, trial counts and exported MATs.
Use the [main handoff](INLINE_EPOCHING_HANDOFF_20261112.md) for subsequent QC and
ERF/dRSA steps with this same run label after P10 is resolved.

## Evidence and implementation

The fix is confined to `epoch_input_preparation.py`. It uses MNE's own sample
mask and the same public array notch operation used by `Raw.notch_filter`.
Ordinary recordings without unsupported fragments still call the original
helper. The model-bound `continuous_cleaning.py` has not been modified.

Audit results, source backups, tests and deployment checks are preserved under
`derivatives/reports/ica_review_updates/short-fragment-fix-20261117-v1/`.
The archived sandbox tests passed 86 cases; release tests exercise the actual
production helper, both ICA modes, boundary handling, epoch identity and MAT
export. No cohort jobs are submitted by this implementation.
