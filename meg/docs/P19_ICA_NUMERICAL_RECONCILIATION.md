# P19 final ICA numerical reconciliation

Date: 2026-11-11. This records numerical evidence and optimizer provenance,
not human component decisions or cleaned-data acceptance. Damiano requested
autonomous investigation and repair of the nonconvergence warning. No
participant exclusion, mask edit, rank reduction or tolerance relaxation is used.

## Evidence and cause

The original final model's SHA-256 is
`8785158b3aa3578af6994d20a356c8e426b5ebf9f73fb92852bdb0485ca806c0`.
It fitted 601,558 samples and retained 69 PCA directions selected by
`n_components=None`. The original cluster stderr explicitly contains sklearn's
ConvergenceWarning at the 1,000-iteration cap. A matched-stack 10,000-iteration
attempt still warned after 4,108.7 seconds. Its weights were not saved before a
network interruption; that attempt is negative convergence evidence, not an
available replacement model.

Source-bound diagnosis used MNE 1.10.2, NumPy 1.26.4, SciPy 1.13.1 and sklearn
1.5.2. This matches the logged/pinned numerical packages, not the Linux platform
or BLAS implementation. Full source/model/accepted-policy validation passed
before and after diagnosis. Whitening covariance error was 1.90e-14 and
unmixing orthogonality error 1.80e-13. The retained-versus-discarded PCA boundary
was approximately 0.04394 versus 2.76e-8 in explained variance.

One unchanged parallel/logcosh update had a stopping residual of 0.40767;
two updates returned within 6.65e-6 of the original directions, versus the
original tolerance of 1e-4. Permutation matching did not remove the one-step
discrepancy. This confirms a near two-cycle, not simple IC label exchange.
The update is approximately a reflection in the old IC coordinates. A
non-Gaussian direction (excess kurtosis about 1.98) is distributed across several
old ICs. The earlier suggestion that weak Gaussian directions alone explain
the failure was not established. Within-good-piece spectra did not identify
a coherent narrowband pair. These observations identify an optimization trap;
they do not identify its physiological source or select artifact exclusions.

## Repair and exact interpretation

Compute the ordinary FastICA proposal from the saved whitened unmixing matrix
W and the original whitened fit data. Express it in the old IC coordinates,
symmetrize that curvature matrix and obtain its orthonormal eigenvectors V.
Initialize the restart with `V.T @ W`; run **standard sklearn FastICA** with
parallel updates, logcosh, the same 1e-4 tolerance and 1,000-step phase cap.
The observer used for progress/checkpoints returns sklearn's update unchanged;
a real-optimizer test checks exact equality to the unobserved implementation.

This is a deterministic restart of the original seed-97 solution. It passes
`random_state=97` with an explicit initializer; sklearn draws no fresh random
matrix in that case. Do not describe it as a fresh random seed-97 MNE fit.
The original MNE PCA/whitener/channel basis is retained exactly and actual
sklearn output is represented in that basis with MNE's usual mixing inverse.
The original phase used 1,000 iterations; the **restart phase used 114**.
The 10,000-step negative attempt is a separate attempt, not this trajectory.

The saved candidate's SHA-256 is
`e53af044aecfa61c06895a276591284f9f6602c921d491c2114ca97ce90b58d0`.
It generated no convergence warning. A further unchanged update of the saved
weights had residual 8.20941e-5, below 1e-4. Rebuilding protected fit preparation
matched every interval/filter/sample field; only the cluster/Mac input mount
prefix is normalized for comparison. Every original source binding remains
checked. All 69 directions and all 601,558 samples remain, with identical PCA,
means and whitening. Method, contrast, tolerance, configured cap, seed ancestry,
BAD masks and 1-40 Hz / 200 Hz preparation are unchanged.

A native 10 s no-exclusion reconstruction had MEG relative error 4.23e-15;
STI101/MISC samples were exactly equal (including NaNs). This is a numerical
matrix/application sanity check, not full applied-cleaning QC or task-signal
preservation acceptance. No cleaned continuous file or production epoch is
created by the numerical reconciliation.

Display QA of current IC000 illustrates a later review question: its fit-copy
standard deviations are approximately 0.036, 0.038, 0.037, 1.598 and 1.569 AU
across runs 1-5. There are no exactly zero fitted activations. Its pale early
image rows therefore reflect much smaller activation, not missing trials or
blanked data; the shared-scale native continuous panels also show the run
difference. Inspect this nonstationarity during component review. It identifies
neither a physiological cause nor a component/participant exclusion by itself.

## Auditable artifacts and publication gate

Shared candidate namespace:
`derivatives/ica_convergence/sub-19/curvature-restart-20261111-v1/`.
It contains the original and candidate FIFs, full original fingerprint,
initializer/final matrix files, diagnosis/probe results, script snapshots and
SHA-256 bindings, exact serialized-weight checks and a comparison of both
decompositions. Matches are descriptive; no old IC index is transferred.
The best one-to-one absolute activation correlations have median 0.719 and
minimum 0.266; 45 of 69 are below 0.9. Thus this is a materially different
decomposition, not merely relabeling/sign changes or a repair of seven indices.
Those descriptive numbers are not acceptance thresholds. Review all current ICs.

`result.json` records completed numerical and before/after source validation.
`publication.json` exists only after transactional installation, complete
regeneration of all component evidence with refitting forbidden, and repeated
ordinary validation. Its `published` and `diagnostic_trials_unchanged` fields
must both be true before treating this as the operative review model.

Publication archives the old model/metadata/report/book/figures/tuning/manifest
under `reports/sub-19/ses-01/history/ica_review/p19-curvature-restart-20261111-v1/`.
It leaves tuning unchanged and records explicit `optimizer_restart` history in
the new metadata and both review reports. Incompatible restart caches refuse
silent default refitting. Only P19's model changes; the other 34 hashes remain
bound to the immutable array-11197 audit.

The first publication attempt stopped at the generated-output stable-positive
timestamp gate on SMB, before rendering. It restored the original model/evidence.
The retry explicitly finalizes timestamps on the newly written model/metadata,
checks exact model bytes and still requires stable positive stamps. It reuses
the complete same-candidate archive only after inventory/stamp checks (plus
stage-specific manifest checks that preserve unrelated lineage/notes). Source
timestamps and all source/model validations remain untouched. Rollback and
changed-archive rejection have real filesystem tests. Attempt logs are retained.

The joined readiness snapshot is
`reports/ica_review_updates/final-ica-readiness-20261111-v1/summary.json`.
Require its successful existence before using the cohort handoff. Its worker
join checks all current model hashes, all property/continuous/task sidecars and
counts, bound evidence stamps and exact operative retained-trial metadata.
It does not reread all raw headers or large PNG/HTML bytes. The ordinary P19
publication validation does reread source headers; other participants retain
their earlier validated unchanged model/source evidence.

## Next human decision

Review all 69 **current-model** P19 ICs using the
[final ICA guide](FINAL_ICA_COMPONENT_REVIEW_GUIDE.md), alongside the other
working participants. Numerical convergence does not determine whether an IC
is neural, artifact or mixed. Explicit model-bound component decisions,
application QC, P10's native slow-contamination assessment, cleaned-output
acceptance and epoching prerequisites remain later gates.

Primary implementation references:
[sklearn 1.5.2 FastICA update/stopping code](https://github.com/scikit-learn/scikit-learn/blob/1.5.2/sklearn/decomposition/_fastica.py)
and [FastICA initialization/API](https://scikit-learn.org/1.5/modules/generated/sklearn.decomposition.FastICA.html).
The reflection/curvature diagnosis and initializer are project-specific numerical
observations and a documented deterministic restart, not a new physiological
classification method or a claim that this is a documented MNE automatic repair.
