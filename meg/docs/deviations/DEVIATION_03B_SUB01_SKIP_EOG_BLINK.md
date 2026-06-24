# Deviation Log: 03_B skip EOG blink annotation for sub-01

Date: 2026-04-23

This note documents a FLUX deviation according to the template in
`/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/docs/PIPELINE_SPEC.md`.

- Step:
  `03_B_ArtefactAnnotation.py` (blink + muscle annotation before ICA).

- FLUX baseline behavior:
  Blink annotations are created from EOG channel events (in this local script,
  `find_eog_events(..., ch_name="EOG001")`) and added to the `_ann.fif` files.

- Local change:
  For subject `sub-01`, EOG is retained as diagnostics-only by default.
  `find_eog_events` counts and EOG trace plots are reported, but no EOG-derived
  `BAD_blink` annotations are added unless an explicit future configuration
  changes `artifact_annotation.eog_mode` to `annotate`.

- Why deviation is needed for this dataset/task:
  In subject 01, EOG-based blink labels were found to be unreliable for the
  observed ocular-like artefacts in ICA component traces (several clear blink-like
  transients were not tagged by EOG001-based blink annotations).

- Expected impact on preprocessing/analysis:
  `03_B` contributes existing and muscle annotations for sub01 while preserving
  EOG diagnostic evidence. Blink-like ICA components are handled through the
  explicit manual ICA exclusion list rather than adding detector-derived blink
  windows.

- Risk or validation concern:
  No separate blink annotations are available for epoch-level
  `reject_by_annotation`. Ocular removal depends on the manually reviewed ICA
  decomposition. EOG plots and component properties must therefore be reviewed
  together.

- Validation evidence path (report/log/figure):
  - Script logic:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py`
  - Runtime output:
    pass-specific EOG status and plots in
    `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/03_B_ArtefactAnnotation_review_report.html`

## Possible Solutions Considered

1. Keep EOG001 blink annotation for all subjects.
- Pros: closest FLUX parity.
- Cons: misses blink-like artefacts observed in sub-01 ICA traces.

2. Use a combined EOG001/EOG002 policy in 03_B for sub-01.
- Pros: still EOG-based.
- Cons: requires additional detector tuning and validation work.

3. Keep EOG diagnostics, skip default EOG annotations, and use manually
   reviewed ICA exclusions.
- Pros: preserves evidence without trusting a weak detector or a hardcoded
  component index.
- Cons: requires explicit component review.

## Current Decision

Current production decision:

- `03_B` uses diagnostics-only EOG handling for sub01.
- `03_C` does not create ICA-derived blink annotations.
- Ocular component rejection is explicit in `ica.manual_exclude`.
