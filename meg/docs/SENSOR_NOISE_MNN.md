# Sensor-space multivariate noise normalization

New Stage-07 sensor analyses default to `--neural-normalization mnn`.
`--neural-normalization none` reproduces the previous metric explicitly.
Stages 08--10 accept the same selector. `DAD_DRSA_NORMALIZATION` supplies it to
cluster jobs. Existing result files remain readable as historical results;
recomputing an old MATLAB handoff requires generating an explicit new handoff.
No raw processing, Maxwell filtering, epoching or Stage-06 export is changed.

## Data path and estimator

Stage 07 validates the accepted input contracts, reads retained Stage-04
prestimulus epochs through public MNE, and retains channels finite across all
three Stage-06 conditions. The original tensors remain in T and T/m.
One fixed transform per participant, ROI and temporal-preprocessing profile
is fitted to pooled background activity across conditions and acquisition runs.
It is reused for every identity, neural time and stimulus model. Covariance is
estimated within the selected ROI; whole-head whitening is not followed by
selecting a sensor subset.

The `sensor_noise_mnn_v1` estimator resamples negative-time epochs to 120 Hz
with FFT edge padding, matches the target neural Gaussian smoothing (15 samples
for A_dRSA_alt; no such smoothing for A_dRSA), and uses -0.4 through -0.2 s
inclusive. Smoothing support excludes stimulus onset. Channel means are removed
for covariance estimation; no new baseline subtraction is applied to task data.
All retained epochs contribute, including unequal identity repetition counts.

Write D for channel noise SDs, z = D^-1 x and R for their noise correlation
matrix. Retain eigenmodes above 1e-6 of R's largest eigenvalue. Shrink retained
eigenvalues toward their mean; select alpha from
`0, 1e-6, 1e-5, 1e-4, .001, .01, .05, .1, .2, .5, 1`
using acquisition-run-held-out Gaussian noise likelihood. Apply the symmetric
retained inverse square root of R to z. Then run the unchanged 1-Pearson neural
RDM and existing dRSA code. Whitened Pearson has a fixed spatial-centering
convention; it is not a rotation-invariant Mahalanobis angle. Crossnobis answers
a different question and is not introduced here.

`noise_normalization/noise_transforms.mat` stores labels, SDs, covariance,
eigenvalues, rank, alpha and symmetric transforms. Its JSON sidecar records
source file bindings, code hashes, physical units, preprocessing and software,
and run-held-out likelihood/calibration diagnostics. Holdouts select alpha;
these diagnostics are not unbiased nested-cross-validation performance scores.
Rank is an effective numerical rule, not a declaration of physical SSS rank.
Prestimulus activity models ongoing background, not purely instrumental noise.
Unequal trial counts and angular-distance measurement bias remain limitations.

## Isolation and resuming

Source epochs, exports, sidecars and their manifest stay read-only. Default
labelled outputs are sibling processing-run directories named
`SOURCE-drsa-mnn-v1` or `SOURCE-drsa-none-v1`, outside legacy recursive searches.
`DAD_DRSA_ANALYSIS_LABEL` selects a fresh simple label when settings change.
Stage-07 handoffs, reports and derived manifests live in the variant directory.
Resume verifies export evidence, normalization artifacts, settings and helper
hashes. Changed inputs/settings require a fresh namespace. Stage 08--10 readers
reject normalization mixtures before participant grouping. Simulations retain
their own model-only estimand.

Sensor groups are unchanged, now stored once in `megin_sensor_groups.json`
for Python fitting and MATLAB selection. The real-neural whole-head A_dRSA
profile uses its original complete retained MEG axis; alt uses its original
predefined groups. Stage-09 reports use separate paths for additional sensor
groups so a full expansion cannot replace the occipital report.

## Approved production expansion

`meg/cluster/run_drsa_mnn_cohort.py` prepares selectors from saved accepted epoch
settings, performs whole-cohort read-only preflight in parallel, and submits
participant and reporting arrays followed by dependent group/statistics jobs.

```sh
python meg/cluster/run_drsa_mnn_cohort.py prepare --run-dir /absolute/run
python meg/cluster/run_drsa_mnn_cohort.py preflight --run-dir /absolute/run --workers 4
python meg/cluster/run_drsa_mnn_cohort.py submit --run-dir /absolute/run --max-concurrent 30 --cpus 4 --memory 24G --hours 24
```

The user-approved cohort is exactly:
`04 05 06 07 08 09 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 31 32 33 36 37 38 39`.
There are 30 participant report folders and one group folder. No additional
export-ready participant is silently added. All five sensor groups and nine
model sources from v2-drsa-expansions are selected explicitly, yielding 55
condition/model/ROI matrices per participant. Observed models cover all three
conditions; predicted and seven PE sources apply to occluded-deviant. The
launcher verifies the whole cohort and MNN provenance before group publication.
Numerical-library threads are capped at one; MATLAB is capped at the allocated
CPU count. Expose concurrency, CPU, memory and duration controls at submission.

## Validation

Synthetic checks cover unit conversion, sensor permutation, supported rank,
whitening and output isolation. Existing Stage-07/08/09/10 tests preserve legacy
calculations by explicitly selecting none. Fresh R2020a entrypoint calculations
on P04 reproduce the saved unnormalized matrix to 8.9e-16 and the approved MNN
sandbox to approximately 1.2e-12. The full cluster path through public MNE
baseline reads agrees with the sandbox dRSA matrix within 3.3e-8. Full production cohort preflight must pass
before submission. These checks validate implementation, not new scientific
acceptance of other conditions or evidence for prediction from lag peaks.

See [the integration validation record](MNN_INTEGRATION_VALIDATION_20261124.md)
for the coordinated repository update and exact MATLAB comparisons.
