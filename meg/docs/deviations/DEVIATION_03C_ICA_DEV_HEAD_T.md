# Deviation Log: 03_C ICA `dev_head_t` mismatch across runs

Date: 2026-04-23

This note documents a FLUX deviation according to the template in
`/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/docs/PIPELINE_SPEC.md`.

- Step:
  `03_C_ICA.py` (FLUX ICA sequence, concatenation before ICA fit).

- FLUX baseline behavior:
  In FLUX `MNEPython/ICA.ipynb`, run-01 and run-02 are concatenated directly via
  `mne.io.concatenate_raws([raw_resmpl_all, raw_resmpl])` after resampling and
  1-40 Hz filtering, with no explicit `on_mismatch` handling and no explicit
  head-position realignment step in the ICA notebook.

- Local change:
  In local script
  `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA.py`,
  concatenation is set to:
  `mne.io.concatenate_raws([...], on_mismatch="warn")`.

- Why deviation is needed for this dataset/task:
  Local subject-01 run files differ in `info['dev_head_t']` (head position
  between runs), which causes MNE to raise by default during concatenation.
  Without this change, ICA fitting on concatenated runs crashes before fitting.
  Estimated mismatch in this session (relative to block01):
  - translation: ~2.1-3.1 mm
  - rotation: ~2.5-4.9 deg

- Expected impact on preprocessing/analysis:
  ICA can be fit on all concatenated runs (FLUX-like workflow preserved).
  However, because runs are not explicitly transformed to a common destination,
  cross-run sensor-to-head alignment is imperfect. This can slightly affect the
  decomposition basis compared with a fully realigned workflow.

- Risk or validation concern:
  Potentially suboptimal ICA decomposition due to mixed head positions across
  runs. Must verify component interpretability (ocular/cardiac components) via
  source traces, topographies, and before/after artefact attenuation plots.

- Validation evidence path (report/log/figure):
  - Script handling:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA.py`
    (concatenation with `on_mismatch="warn"`).
  - Runtime evidence:
    previous traceback from `03_C_ICA.py` before this fix (MNE
    `ValueError: raws[1].info['dev_head_t'] differs`).
  - Quantified mismatch evidence:
    subject-01 `dev_head_t` comparison against block01 yielded
    translation ~2.1-3.1 mm and rotation ~2.5-4.9 deg across runs.
  - Post-fix QC outputs (after rerun):
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/03_C_ICA_report.html`
    and ICA figures under
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/figures/`.

## Possible Solutions Considered

1. Keep current approach (`on_mismatch="warn"`), concatenate, fit one ICA.
- Pros: minimal change, closest to FLUX logic, fast to proceed.
- Cons: does not explicitly solve cross-run head-position differences.

2. Re-run Maxwell filtering with a common head-position destination for all runs,
   then concatenate and fit ICA.
- Pros: best geometric consistency across runs; recommended long-term fix.
- Cons: requires rerunning SSS/maxfilter outputs and downstream steps.

3. Fit ICA separately per run (no concatenation across mismatched `dev_head_t`).
- Pros: avoids head-position mismatch in fitting.
- Cons: deviates more from FLUX two-run concatenated ICA pattern and increases
  handling complexity (multiple ICA models).

## Superseded Decision

Date superseded: 2026-06-12

The temporary native-position workflow is no longer the production decision.
Production 03_A uses each run's static `dev_head_t` as its source transform and
maps all runs to the most central actual run. Centrality is the transform
medoid measured by RMS displacement over the subject's digitization point
cloud. The destinations are sub01 run02, sub02 run05, sub03 run03, sub04 run04,
sub05 run02, sub06 run03, and sub07 run03.

`on_mismatch="warn"` remains in 03_C only as defensive compatibility for
historical derivatives. Rebuilt stage-03 derivatives should share one
destination transform and must not depend on that warning path.

---

## Superseded 03_C Deviation: sub-01 blink annotation source (ICA000)

Date: 2026-04-23

- Step:
  `03_C_ICA.py` (ICA fitting, inspection, and application).

- FLUX baseline behavior:
  Blink annotations are expected to come from the pre-ICA artefact-annotation
  step (EOG-based labels in annotated raw files).

- Historical local change:
  The legacy implementation generated blink annotations from ICA component 0.

- Why deviation is needed for this dataset/task:
  Subject 01 EOG-based blink labels were not sufficiently aligned with observed
  blink-like transients in ICA source traces; subject-specific ICA000 labeling
  is used as an operational workaround.

- Current production status:
  This workaround was removed in the configuration-driven refactor. Component
  order is arbitrary, so production ICA never assumes ICA000 is ocular.
  Sub01 EOG is diagnostics-only and ocular cleaning uses explicit manual
  component exclusions.

- Risk or validation concern:
  The legacy approach could label non-ocular transients. Removing it avoids
  index-based false assumptions but requires careful manual review.

- Validation evidence path (report/log/figure):
  - Script logic:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA.py`
  - Current review reports and component diagnostics under:
    `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/`
