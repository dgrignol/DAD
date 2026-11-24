# MNN integration validation - 24 November 2026

## Accepted behavior and provenance

Sensor dRSA defaults to one fixed participant/ROI/profile MNN transform followed
by the existing Pearson RDM. Explicit `--neural-normalization none` or
`DAD_DRSA_NORMALIZATION=none` selects the legacy metric. Source FIFs, SSS,
epoching and Stage-06 exports are unchanged. Derived analyses use separate
namespaces; altered code/settings require a fresh namespace, not a rewrite of
an accepted run. MNN is fitted from pooled prestimulus activity, independent of
condition-specific task amplitudes. Its limitations remain those documented in
[SENSOR_NOISE_MNN.md](SENSOR_NOISE_MNN.md).

The source DAD branch contains the reviewed feature commits `51dda15`,
`9d621dd`, `d90c6fc`, `fc4aabf` and `5e96286`, based on `82aca23`.
The cluster integration starts from dad-cluster-analysis main `a6cdb3c` and
ports their delta while retaining the existing parallel Stage-10 implementation,
shared PE model builder, diagnostic CLI and group-only renderer. No unrelated
repository history is merged between the two repositories.

The standard deployment payload now includes the MNN and diagnostic Python
helpers, launchers, report tests and group-report policy. The cluster expansion
wrappers propagate source run, output label and normalization independently;
readers receive the same output label as the producer. Historical diagnostics
select none explicitly. These are preparations for future deployments; existing
frozen numerical snapshots are not changed.

## Integration checks

- DAD: 104 relevant Python checks passed (66 dRSA/diagnostic tests and 38 cluster
  wrapper tests).
- dad-cluster-analysis: 99 relevant Python checks passed, including serial vs
  parallel statistics, expansion dispatch, explicit none/MNN selection,
  normalization provenance, immutable diagnostic publication and group-only
  linked/embedded reports. A further 20 affected checks passed after the final
  diagnostic-root and CLI-example reconciliation; these overlap the 99.
- Changed Python files parse; changed shell files pass `bash -n`; code whitespace
  checks pass. Cluster tests retain its existing parallel-equivalence coverage.
- Fresh local MATLAB R2020a entry points, using archived read-only P04 input
  copies, reproduce both saved always-visible occipital dRSA matrices (MNN and
  none) with maximum absolute difference **0** in both repositories.
- The observed model self matrix and all five MNN neural self matrices match
  their saved diagnostic references with maximum absolute difference **0**.
  Every one of the nine occluded-deviant model builders matches the reference
  tensor, source condition and feature metadata exactly.

All newly generated MATLAB outputs and detailed numerical evidence are outside
the project and its production derivatives. This is a representative integration
check, not a new full-cohort calculation or an additional scientific inference.
It reuses the validated noise artifact; it does not repeat the noise fit.

## Completed production results retained

The existing `noica-20261117-v2-drsa-mnn-v1` run remains the accepted output
namespace. Jobs 13012-13015 and 13101 completed successfully. The prior full-run
audit verified exactly the approved 30 participants, 1,650 brain-model MATs,
group summaries/statistics and the 50 PNG/NPZ RDM-dependence pairs. The Stage-07
group linked and embedded reports contain 55 group figures per variant;
participant figures remain in participant reports. Production arrays, PNGs,
JSON, source manifests and frozen runtime hashes are retained unchanged.

## Scope for subsequent work

RDM dependence remains a separate diagnostic workflow. Replacing Stage-07
reports or adding these reports to the default main launcher requires the
separate workflow-coverage review; it is not part of this integration. A merge
records reviewed code for future use; it does not deploy, submit jobs, rerun
statistics or retroactively change the code provenance of completed results.
