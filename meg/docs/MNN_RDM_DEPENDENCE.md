# Expanded RDM dependence in the MNN namespace

The Stage-07/08/09/10 jobs 13012--13015 do not produce the separate expanded
RDM dependence report. This addition submits its own dependent workflow without
changing those running jobs or their frozen code. It does not yet replace the
Stage-07 HTML report or integrate diagnostics into the default main launcher.
That presentation/orchestration decision is reserved for a separate review.

## Inputs and unchanged model estimand

Source epochs/exports remain `noica-20261117-v2`. Brain-model matrices,
completed execution reports and the exact saved MNN transforms come from
`noica-20261117-v2-drsa-mnn-v1`. The accepted cohort is exactly the same 30
participants recorded in `SENSOR_NOISE_MNN.md`. Five sensor groups and all nine
model sources are selected in their valid conditions.

The reviewed expanded renderer, validators and model builder were adapted from
`dgrignol/dad-cluster-analysis` revision
`a6cdb3cdf70e54b4169c7f38f659d8cf50d94e43`. That checkout and the ERF worktree
are read-only during this integration. The extracted model builder is compared
with the MNN branch's existing A_dRSA_alt builder, rather than assuming that
matching function names establish equivalent scientific definitions.

Model-only Euclidean RDM self/cross correlations do not change with sensor MNN.
MEG correlation RDM self-correlations must use exactly the same fixed transform
as the brain-model matrices displayed alongside them. The diagnostic MATLAB
code applies the producer's saved ROI transform after its unchanged 15-sample
smoothing/ROI selection; it never refits normalization. It saves normalization
metadata with every neural self panel. Python binds the transform bytes and
contract, rejects raw/MNN mixtures, validates all requested brain matrices,
and records the producer analysis contract.

Undefined model Pearson correlations remain NaN, including diagonals. The
130-frame observed/predicted shared prefix is retained as a validation check.
Model definitions, signed lags, averaging and SEM conventions are unchanged.

## Report layout and resources

The expanded report retains 15 overviews and 35 focused PE pages (50 PNG/NPZ
pairs), with +/-2 s lags, 200-ms sliding row windows, 25-ms steps and at least
three finite samples per cell. Summary statistics are computed per participant
before equal-weight finite-participant reduction; SEM requires at least two.
The existing color palette and layout are reused. The report identifies MNN
explicitly on its neural panels and in its provenance.

The standalone cohort runner supports `--neural-normalization mnn` (default)
and explicit `none`. MNN requires compatible production inputs and a saved
transform. No original input manifest, Stage-07 result or previous diagnostic
report is modified. All outputs use a fresh child of the MNN result namespace.

Example on DIPSEN, from a separate frozen diagnostic code snapshot:

```sh
python meg/cluster/submit_mnn_rdm_dependence.py \
  --run-dir /home/USER/runs/dad-mnn-dependence-20261121 \
  --production-run-dir /home/USER/runs/dad-mnn-20261121 \
  --brain-results-root /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/processing_runs/noica-20261117-v2-drsa-mnn-v1 \
  --diagnostics-root /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/processing_runs/noica-20261117-v2-drsa-mnn-v1/rdm_dependence_expanded_mnn_20261121 \
  --dependency 13012
```

Dependencies are the complete MNN producer array -> four-worker cohort preflight
-> all 30 independent participant tasks -> four-worker group join/render.
Participant tasks request 2 CPUs and 24 GB and bound MATLAB/library threads to
the allocation. Group work requests 4 CPUs and 16 GB; numerical workers use one
thread each. Shared publication is serial after all inputs and summary stacks
join. Controls for workers, task CPUs, concurrency, memory and time are exposed.
Private per-job caches avoid interfering with other MATLAB/plotting jobs.
The final page is `group/lag-2s/ses-01/index.html` under the diagnostic root.

## Validation

Synthetic regressions cover lag orientation, prefix equality, finite support,
SEM, immutable artifacts, missing participants, PE metadata, MNN/raw mismatch,
and producer/preflight/array/group submission order. A real P04 R2020a check
compares all nine model-builder outputs exactly with the existing branch's
model code. Model self-correlations agree with the original reference within
2.9e-15; MNN occipital neural self-correlations agree with the approved sandbox
within 3.8e-8 (the public MNE baseline-read path). All five ROI axes and MNN
metadata pass validation. These are implementation checks, not significance
or mechanistic claims from descriptive correlation panels.

## Submission-limit fallback used on 21 November

The scheduler accepted preflight 13080, then rejected the participant array
because the account had reached its 200 submitted-task limit. That pending
preflight had done no work. Only this own pending preflight is replaced; no
existing producer, decoding, ERF, summary or statistics job is altered.

`run_mnn_rdm_allocation.py` runs the complete diagnostic workflow within one
allocation, dependent on the entire producer array 13012. Request 60 CPUs,
192 GB and four hours: validate with four workers, compute all 30 participants
with 30 independent processes and two threads each, then join and render with
four group workers. The measured approximately 4-GB participant footprint
budgets about 120 GB before headroom; 192 GB is an aggregate allocation rather
than 30 separate MATLAB startup limits. Each pool worker has private MATLAB
and plotting caches. All eligible participants enter the pool together.
Shared publication remains ordered after all participant processes complete.

The control directory preserves the failed submission record and the final
single-job ID, resource settings, code revision and phase progress. This is a
scheduler/resource change only; inputs, matrices, normalization, summaries and
output paths are identical to the array workflow.
