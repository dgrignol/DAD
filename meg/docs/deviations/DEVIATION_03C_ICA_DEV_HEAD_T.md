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

## Current Decision

For now, keep the current local solution:

- `on_mismatch="warn"` in `03_C_ICA.py` concatenation.

This is accepted as an explicit temporary deviation to preserve pipeline
continuity and FLUX-style structure, with later upgrade path to a common
Maxwell destination if needed.

---

## Additional 03_C Deviation: sub-01 blink annotation source (ICA000)

Date: 2026-04-23

- Step:
  `03_C_ICA.py` (ICA fitting, inspection, and application).

- FLUX baseline behavior:
  Blink annotations are expected to come from the pre-ICA artefact-annotation
  step (EOG-based labels in annotated raw files).

- Local change:
  For subject `sub-01` only, blink annotations are generated in `03_C` from
  ICA component 0 (`ICA000`) and injected into:
  1) the concatenated source-inspection raw used by `ica.plot_sources`, and
  2) each saved `_ica.fif` output file.

- Why deviation is needed for this dataset/task:
  Subject 01 EOG-based blink labels were not sufficiently aligned with observed
  blink-like transients in ICA source traces; subject-specific ICA000 labeling
  is used as an operational workaround.

- Expected impact on preprocessing/analysis:
  Sub-01 blink labels in `_ica.fif` outputs now reflect ICA000 event detection
  instead of EOG001-based annotations from `03_B`.

- Risk or validation concern:
  ICA000 can include non-ocular transients; ICA000-based blink labels may include
  false positives relative to pure blink physiology.

- Validation evidence path (report/log/figure):
  - Script logic:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA.py`
  - Report log entries under:
    "Blink annotation source" and "ICA-applied outputs" in
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/03_C_ICA_report.html`
