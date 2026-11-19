# ICA review status - 2026-11-12

The component choices are **provisional**. Damiano has not manually inspected
the proposals and needs to return to inspect them. Keep this reminder in ICA
review documents and receipts; downstream outputs need no repeated labels.
Component authorization, application QC and final scientific acceptance remain
separate. Deferred ICs remain retained, with their questions open.

The original frozen authorization covers 34 participants and 29 zero-based
exclusions. Separately recovered P32 has explicit user authorization for
**[5, 34]**, retaining its other 69 ICs, including 35 deferred components.
Together these cover 35 participants, 180 acquisitions and 31 exclusions.
P15 remains excluded from processing, with its data/provenance retained.

Current recording state:

- P32's guarded Stage 4C recording succeeded, with independent authorization,
  source recovery, complete proposal and before/after snapshots under
  `derivatives/reports/ica_review_updates/p32-source-recovery-20261112-v1/`
  and `p32-ica-recording-20261112-v1/`.
- The original 34-participant recorder was stopped before its whole-table write
  barrier. Its last preserved state was 24/34 validated and **0/34 recorded**.
  The stopped worker processes have been terminated; do not resume that stale
  process tree or infer completion from partial files. Evidence remains under
  `provisional-ica-recording-20261112-v1/` and `...-v2/`.
- The frozen original authorization/choices remain under
  `provisional-ica-authorization-20261112-v1/`. P32's separate authorization
  supplements the original null entry; it is not an inferred empty list.
- Future ICA-enabled production requires finishing current exact-model validation
  for the entire applicable table before any bulk recording writes, through
  `record_provisional_ica_choices.py` and the guarded Stage 4C workflow. Preserve
  existing receipts/conflicts and snapshots. No need to ask the user to authorize
  the already accepted lists again. New choices remain separate decisions.

**The next step has changed: do not produce continuous ICA derivatives now.**
The user approved proceeding without component rejection and preparing line-noise
removal in memory before epoch cutting. Future selected ICA rejection will be
inserted immediately before that same line removal. See the executable
[continuous-preparation/epoching handoff](INLINE_EPOCHING_HANDOFF_20261112.md).
That route preserves the saved models and this review work, while avoiding about
182.5 GiB of permanent continuous copies.

The previous full continuous-application handoff was archived before replacement
under `derivatives/reports/ica_review_updates/inline-epoch-preparation-20261112-v1/source_before/`.
No cohort ICA application, production epochs or cluster job was executed by Codex
as part of this implementation. Actual new-run QC and scientific review follow
user-submitted processing. P10's approximately 0.33 Hz residual and participant
retention remain open, particularly until a later ICA-enabled comparison.
