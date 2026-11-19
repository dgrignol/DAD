# P24 EDF auxiliary reconstruction - 2026-11-17

The user requested EDF-derived eye channels in **all five acquisitions**, zeros
where EDF samples are unavailable, and omission of MISC008 from the new run-5
epoch input. Original acquisition/SSS/annotated FIFs, EDFs, accepted annotations,
visual timing and `noica-20261117-v2` outputs remain unchanged. No jobs have been
submitted, no deployment has been performed, and no new epochs have been made.

## Prepared inputs

The immutable auxiliary manifest is at:

`derivatives/auxiliary_reconstruction/p24-edf-20261117-v1/manifest.json`

It binds five new compressed arrays, individual validity masks, the original
native coordinates/annotations/file stamps, the canonical trigger CSV, accepted
trial-local EDF mapping files and reconstruction implementation. The six fields
are **left pupil/X/Y, right pupil/X/Y** in MISC001-006. Pupil remains in EDF area
units; gaze remains in pixels. Run 5 is replaced too, so units are consistent.

Only uniquely matched recorded trial segments receive measured nearest EDF
samples. Invalid samples, missing eyes, recording-off gaps, unmatched trials and
out-of-coverage periods are zero. Separate validity masks distinguish filled
zeros from measured zero coordinates. No interpolation bridges a gap; no clock
drift, extra lag or physiological fit is applied to the output arrays. These
channels do not introduce new BAD annotations or rejection decisions.

The opt-in runtime adapter supplies all six fields in every acquisition, omits
MISC008 and orders every run to **325 channels, including all 306 MEG channels**.
It verifies retained native samples, triggers, annotations and acquisition
settings. Runs 1-4 retain their acquired 0.1 Hz high-pass; run 5 retains 0 Hz.
The original run-5 photodiode stays available for the accepted timing evidence.

Use the new selector `meg/params/sub-24-edf-aux.yaml`, with the original reviewed
tuning and Stage-3 timing. Its new label is `noica-p24-edf-20261117-v1`.
An exported `DAD_PROCESSING_RUN` overrides the selector, so select that same new
label explicitly or unset the variable when using the selector. The ordinary
array wrapper still chooses `sub-24.yaml`; it does not automatically opt in.

## Validation and comparison report

All five native sources, accepted inputs/timing and the reconstruction passed
the actual Stage-04 **validation-only** route. Evidence is saved locally at:

`meg/sandbox/p24_edf_auxiliary_review_20261117/preflight.log`

The comparison report is:

`meg/sandbox/p24_edf_auxiliary_review_20261117/index.html`

It reports every measured correlation, including disagreement. On the stable
covered core, original analog versus EDF correlations are:

| Field | Left | Right |
| --- | ---: | ---: |
| Pupil | 0.9476 | 0.9999 |
| Gaze X | 0.8006 | 0.9924 |
| Gaze Y | 0.6886 | 0.9927 |

The stable diagnostic core removes existing BAD intervals, 100 ms around EDF
recording edges, and EDF blink intervals padded by 100 ms. About 355 seconds
per eye remain. These are report-only masks, not modifications of annotations
or reconstructed samples. The report separately includes all valid coverage,
which has weaker correlations. A labelled -2 ms analog-delay sensitivity is
diagnostic only and does not alter the reconstruction or MEG timing.

Near-perfect agreement is supported for the **right eye's stable covered
periods**, not for all six channels or all coverage. The source of the remaining
left-eye disagreement has not been established. No scientific acceptance of
new eye-based analyses follows from this auxiliary reconstruction.

## Subsequent execution

Deployment must include `p24_auxiliary_reconstruction.py`, the updated
`epoch_input_preparation.py` and `pipeline_config.py`, and the new selector.
The deployment manifest includes these files. The new arrays already reside in
the shared derivatives tree. Existing completed participants can retain their
verified old epochs; a later cohort handoff must explicitly bind P24's separate
recovery run rather than silently changing the cohort's global input label.

P24 still needs user-submitted epoching/export in the new label after the
amplitude-rejection policy is settled. Its failed v2 calculation saved no
epochs, so continuous line removal and epoch construction must run again.
SSS, annotations and accepted timing can be reused. The auxiliary substitution
does not change the current amplitude thresholds or predict final identity
coverage after epoching.
