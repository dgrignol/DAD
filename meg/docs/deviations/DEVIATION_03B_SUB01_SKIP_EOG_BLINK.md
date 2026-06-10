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
  For subject `sub-01` only, EOG-based blink annotation is intentionally skipped
  in `03_B_ArtefactAnnotation.py`. The blink-only QC figure is also skipped and
  replaced with an explicit warning in the report/log.

- Why deviation is needed for this dataset/task:
  In subject 01, EOG-based blink labels were found to be unreliable for the
  observed ocular-like artefacts in ICA component traces (several clear blink-like
  transients were not tagged by EOG001-based blink annotations).

- Expected impact on preprocessing/analysis:
  `03_B` now contributes only event + muscle annotations for sub-01. Blink labels
  for sub-01 are deferred to step 03_C, where they are derived from ICA component
  0 and reattached to outputs.

- Risk or validation concern:
  ICA000-based blink detection can capture non-blink transients if component 0 is
  not purely ocular. This is accepted as a temporary subject-specific workaround.

- Validation evidence path (report/log/figure):
  - Script logic:
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py`
  - Runtime output:
    warning text under section "Blink-only QC figure" in
    `/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/03_B_ArtefactAnnotation_report.html`

## Possible Solutions Considered

1. Keep EOG001 blink annotation for all subjects.
- Pros: closest FLUX parity.
- Cons: misses blink-like artefacts observed in sub-01 ICA traces.

2. Use a combined EOG001/EOG002 policy in 03_B for sub-01.
- Pros: still EOG-based.
- Cons: requires additional detector tuning and validation work.

3. Skip 03_B blink annotation for sub-01 and defer to ICA in 03_C.
- Pros: minimal immediate change that reflects observed sub-01 data behavior.
- Cons: introduces cross-step dependency and a subject-specific branch.

## Current Decision

For now, keep this subject-specific deviation:

- `03_B` skips EOG blink annotation for sub-01.
- Blink annotations are added in `03_C` from ICA component 0 for sub-01.
