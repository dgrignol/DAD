# Final ICA: participant-by-participant component review

Updated 2026-11-10. This is the human review guide for the 35-participant working
cohort P04-P14/P16-P39. P15 remains excluded by the participant-specific coverage
decision; preserve its data. P10 remains in processing with its slow-contamination
concern open. Authoritative choices remain in
[PREPROCESSING_DECISIONS.md](/Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/meg/docs/PREPROCESSING_DECISIONS.md).

## Which report to open

Open `reports/sub-XX/ses-01/03_C_ICA_review_components.html` for human review.
It contains the guide, component navigation, maps/properties, source/context
plots and a model-pinned command under every figure. Local commands use the Mac
workspace; expandable cluster commands use the cluster site configuration.
Interface v3 shows BOTH native continuous examples and task-locked trial images
with condition means by default, under every final IC. They are visible figures,
not just commands to generate optional plots. Each run is projected once for all
ICs. The shared native manifest records retained trial identities, native sample
bounds, BAD categories, drop reasons and the accepted timing gate.

`03_C_ICA_review_report.html` is the processing/evidence report: model identity,
fit duration/sample maps, convergence, configuration and tuning. Interface v2
links to the component book instead of duplicating its figures. These are two
views of ONE model, not competing decompositions. Earlier versions are archived
under `reports/sub-XX/ses-01/history/ica_review/LABEL/` before replacement.

## Precise references and how to use each view

The [MNE plot_properties API](https://mne.tools/1.10/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA.plot_properties)
defines the five panels, segment/epoch rejection, PSD standard deviation and
95% parametric confidence band. The [MNE ICA tutorial](https://mne.tools/1.10/auto_tutorials/preprocessing/40_artifact_correction_ica.html#selecting-ica-components-manually)
demonstrates source traces, maps, properties and physiological-reference views
together; it also shows Raw and event-averaged sources separately.
[ICLabel: Telling Components Apart](https://labeling.ucsd.edu/tutorial/labels),
particularly Brain/Eye/Muscle/Heart and their example captions, explains why
temporal morphology and event-related images complement spatial/spectral clues.
It is an EEG tutorial: do not import its scalp/dipole residual-variance thresholds
or classifier labels into this MEG/SSS analysis.

For DAD, use that evidence as follows; these are review instructions, not
automatic classification rules:

| View | Specific question | Interpretation limit |
| --- | --- | --- |
| Mixing map (mag/grad) | Is the spatial pattern compatible with the source hypothesis and the native frontal/MISC evidence? | A dipolar-looking map alone does not identify a neural source. |
| Fit PSD (1-40 Hz) | Is there a reproducible oscillatory peak or a broad spectral shape worth investigating? | This fit-only PSD cannot assess P10's 0.33 Hz issue or muscle/line activity above 40 Hz. |
| Arbitrary 2 s image and variance | Are fluctuations persistent, intermittent, or dominated by a few segments? Investigate variance extremes in source traces. | Segment mean is not an ERF; high variance is not an automatic exclusion, and variance-drop labels are not trial retention. |
| Native continuous traces | Do repeated blinks, steps, bursts, periodic waveforms or drift coincide with labelled BADs/reference traces? Check ordinary controls too. | Displayed 20 s examples are deterministic, not a representative sample. Scroll all acquisitions when needed. |
| Task-locked trial image and mean | Does activation recur at accepted motion-onset latencies, with jitter or condition dependence? Inspect individual trials as well as the mean. | Both artifacts and neural signals can be task-related; absent means do not rule out neural activity. |

Continuous examples include BAD periods, which fit data exclude. Their window
means are subtracted only for display, without variance normalization or native
filtering; each IC uses a common y-scale across its run panels. Task plots use
the configured BAD rejection, baseline and detrend and remain pre-cleaning
diagnostics, not final Stage 5 epochs. Each source interval retains native time;
task zero is accepted motion onset, never the arbitrary segment start.

The joined `native-task-retention-audit-20261110-v1` snapshot checks all 35
participants: exactly 19,988 retained trial identities, all condition/run/identity
cells, the accepted timing gates and native -0.5/+2.667 s bounds agree with the
operative Stage 4B footprint, including the three recorded P36 serialization
deltas. These are pre-cleaning IC diagnostics; the audit creates no final epochs
and establishes no neural/artifact classification or scientific acceptance.
The reproducible metadata audit command is:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 meg/code/preprocessing_basic/audit_ica_task_diagnostics.py \
  --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
  --label NEW-UNUSED-LABEL --jobs 2
```

This supplements the report workflow's source/header validation; it does not
replace that gate. Reuse the existing snapshot unless inputs/reports change.

## P19 numerical reconciliation

The original fit hit its 1,000-iteration cap with an explicit ConvergenceWarning.
A matched-stack 10,000-iteration attempt still failed. Source-bound diagnosis
confirmed near two-step recurrence: one update changed directions by 0.40767,
whereas two returned within 6.65e-6; the original tolerance is 1e-4. Whitening
and the PCA boundary were healthy. This points to an optimizer cycle; it does
not identify a physiological artifact or justify excluding any IC.

The documented deterministic restart of the saved seed-97 solution converged
in 114 restart-phase iterations, with a further-update residual of 8.20941e-5.
It keeps all 69 directions, all 601,558 fitted samples, exact original PCA and
whitening, parallel/logcosh FastICA and the original tolerance/cap. It uses an
explicit initialization derived from update-curvature eigenvectors, not a fresh
random draw. The configuration fingerprint remains unchanged; exact model
hashes and explicit optimizer_restart history record the different trajectory.
All 69 NEW ICs require human review; no old indices/choices transfer.

See [P19 numerical reconciliation](P19_ICA_NUMERICAL_RECONCILIATION.md) for exact
provenance, limitations, artifacts and archive. Require the candidate directory's
publication.json to report published=true and diagnostic_trials_unchanged=true.
The joined cohort readiness snapshot is
reports/ica_review_updates/final-ica-readiness-20261111-v1/summary.json.
These gates establish numerical/software readiness only; final component and
cleaned-output acceptance remain explicit later decisions.

The share disconnected before the negative 10,000-step attempt saved weights.
Its completed warning remains evidence; no nonexistent model is reused.
New attempts preserve actual weights locally as explicitly unvalidated
checkpoints before network revalidation. MNE 1.10.2 resets its loaded max_iter
attribute to a default; convergence accounting reads the persisted fit_params
cap instead. That I/O bookkeeping issue did not cause the original warning.

## What the property panel actually represents

The original pipeline called `ICA.plot_properties(fit_raw, picks=component)`.
For Raw input, MNE creates arbitrary nonoverlapping **2 s segments**, not task
epochs. With 200 Hz fit copies, their sample times are 0 through 1.995 s. Rows
are segments; the black curve is their mean, and the grey region is a 95%
parametric confidence interval. The spectrum's grey region is one standard
deviation across segments. This mean is not a stimulus ERF; absence of a mean
response here says nothing about task locking. IC activation has arbitrary
units and polarity.

The old MNE title, "Segment image and ERP/ERF", was misleading for Raw input.
MNE also clears the image-axis tick locator; because the mean axis shares it,
both can lose numeric x ticks. The plotting helper now restores numeric ticks
and labels the panel "2 s fit segments (not stimulus-locked)". The variance
panel no longer implies that the old "Dropped segments: 0%" measured trial or
annotation retention; MNE's percentage omitted annotation/boundary drops.

An additional **MNE 1.10.x display bug** was reproduced: ICA.get_sources() on
concatenated Raw data can retain only the first piece's metadata duration while
its data array contains all source samples. P08's model fitted **595,205 samples**
at 200 Hz, but its first retained piece has **10,231 samples (51.155 s)**. That
produced the screenshot's **25** complete 2 s segments. Fitting and model weights
were unaffected. Full source-browser transformation and native single-run
context are distinct code paths.

The corrected display uses one public RawArray containing all fit samples and
the original join annotations before plotting. The protected fit preparation,
model files, numerical settings, native data and exclusions are unchanged.
Segments crossing BAD/EDGE joins are omitted from property displays. No refit
is needed for this correction.
Read-only verification on the saved P08 model reconstructed all 595,205 fit
samples and plotted 1,366 retained segments, with numeric 0-1.995 s ticks.

Reports not yet regenerated explicitly say `legacy_first_fit_stretch` and
provide both an exact first-stretch reproduction command and a corrected
full-fit command. `full_fit` identifies regenerated panels. The original
1-40 Hz/200 Hz fit preparation is rebuilt only as a disposable display copy.
The runtime used for display and the fitted model's runtime are recorded.
Some earlier reports (including P04) have no saved property panels despite
requesting them. Their books explicitly remain incomplete. Run the report-only
regeneration below before final participant decisions; it both fills missing
panels and corrects first-stretch displays. Processing HTML can document this
incompleteness without granting component acceptance.

Primary references: [MNE 1.10 ICA API](https://mne.tools/1.10/generated/mne.preprocessing.ICA.html),
[MNE 1.10.2 property-plot implementation](https://github.com/mne-tools/mne-python/blob/v1.10.2/mne/viz/ica.py),
[source projection implementation](https://github.com/mne-tools/mne-python/blob/v1.10.2/mne/preprocessing/ica.py).

## Four different time coordinates

| View | Zero / time axis | Can it establish task locking? |
| --- | --- | --- |
| Property image and segment mean | Start of each arbitrary 2 s fit segment; 0-1.995 s | No |
| Fit-source browser | Compressed concatenated fit time; excluded gaps removed | No |
| Native source/context | Acquisition seconds, absolute native sample divided by native sampling rate | Shows original timing, but is not averaged around task events |
| TASK diagnostic | Accepted Stage 3 corrected motion-onset event, with configured -0.5 to +2.667 s | Shows task-locked IC activation; task locking alone does not classify neural versus artifact |

The TASK command joins the unchanged accepted Stage 2 ledger to the hash-bound
Stage 3 timing/application. Its event is `epoch_event_sample_raw` for condition
codes 71/81/91, not an uncorrected digital trial-start pulse, warmup or catch
question. Estimated motion-onset assignments remain estimates with their
documented uncertainty. It uses the configured native epoch window, baseline
(currently none) and detrend (currently linear, `1`), together with BAD rejection.
Original recording boundaries and incomplete windows retain MNE's rules.

These TASK plots use **uncleaned native IC projections**, before final ICA
removal and line processing. Sensor-amplitude rejection cannot be transferred
to arbitrary IC units; it is not applied. They are transient diagnostics, with
trial/source identities and drop reasons, not saved Stage 5 epochs or final
ERFs. `--include-bad` is a clearly labelled diagnostic comparison only; it changes
no eligibility, annotation or production rule. Routine ocular intervals are
otherwise rejected here, so lack of an ocular response is not proof of absence.

## What the native context figure shows

The figure projects the final model onto one native, unfiltered **5 s** window
without removing components. Its panels contain all IC projections, selected
good frontal magnetometers and original MISC eye channels (excluding MISC008,
the photodiode). Empty/flat MISC is not magnified into evidence of usability.
New frontal trace panels subtract each channel's five-second window mean for
display in pT, so large native DC offsets do not hide its variation. This is
recorded in the figure sidecar; it changes neither projection input nor native
data and introduces no temporal filter. Older panels show the original DC levels.
Their commands retain those DC levels; new centered panels explicitly carry
`--center-frontal`, so commands reproduce the pictured display operation.

The default choice is deterministic and component-independent: select the first
`BAD_blink*` or `BAD_eye_movement*` annotation, start 2 s before it, and clamp the
window to the recording. Without an ocular annotation, start at
`min(20 s, recording midpoint)`, then clamp. This is neither the strongest
artifact nor a representative statistical sample. Every overlapping BAD category,
its original onset/duration and exact native window bounds appear below the
plot. Newly rendered plots also label the shaded categories in their legend.
`--category BAD_muscle` selects the first example of that exact category.

Each IC heatmap row is centered and scaled by its own standard deviation within
that window. Colours show temporal shape, not comparable amplitude or dominance
across components/windows. Co-occurrence with a BAD interval is a useful lead,
not a BAD-enrichment statistic or an automatic reason for removal. Use full-run
source traces, several events, ordinary controls and independent spatial/spectral
evidence to investigate it.

## Review procedure

1. Verify participant, exact prefit SHA-256, current source/annotation bindings,
   component count, fitted versus mapped samples and convergence. The current
   P08 model has 75 ICs and 73 iterations. P19's reconciled model has 69 ICs
   and 114 restart-phase iterations; verify its explicit restart provenance and
   successful publication/readiness evidence above. Its original capped model
   and warning remain archived. Numerical convergence never selects exclusions.
2. Inspect every zero-based final-model IC. Start with maps and spectra, then
   temporal morphology and repeated native examples. In the interactive property
   window, `T` cycles sensor types and `L` toggles log scales. A spatially plausible
   map, large variance or high-frequency power alone is insufficient to decide.
3. For artifact candidates, compare source traces with frontal/MISC context,
   labelled ocular/muscle periods and ordinary controls across acquisitions.
   Consider blink/eye, muscle, cardiac, sensor and mixed-source explanations.
   Do not turn absent/unusable physiological references into negative evidence.
4. Inspect task-locked diagnostics for uncertain or apparently neural/mixed ICs.
   A task-locked response supports investigation, but eye/muscle activity can also
   be task-related. A component may mix signals. These plots alone select no
   exclusion and do not establish preservation of task signals after removal.
5. Maintain a table with IC index, keep/exclude/defer, proposed source, confidence,
   evidence (plot/run/native time) and unresolved questions. Inspect all ICs;
   defer uncertain cases. No preliminary frontal-ICA index transfers here.
6. Discuss a proposed participant-level list. Only after explicit approval record
   model-bound exclusions with `review_cleaning.py --stage 4C --exclude ...`.
   `--exclude` alone is an explicitly reviewed empty list. No plot command writes
   a decision. Browser selections exist only in memory and are printed as unsaved
   marks when the source browser closes.
7. Final application, matched ICA-only before/after QC, final line-removal and
   preservation checks follow. P10's approximately 0.33 Hz concern requires native
   before/after long windows; 1-40 Hz fit plots cannot settle it. Cleaned-output
   acceptance and Stage 5 remain later explicit steps.

## Independent commands and report refresh

Every image's command pins its exact model SHA-256. The following shorter examples
resolve and validate the current P08 model; commands printed in HTML also carry
`--expected-model SHA256`. Run them in the existing MNE environment:

```bash
python3 /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/review_final_ica.py --subject 08 --plot sources --source sub08_run01_210826.fif
python3 /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/review_final_ica.py --subject 08 --plot raw --source sub08_run01_210826.fif
python3 /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/review_final_ica.py --subject 08 --plot properties --component 0 --scope full
python3 /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/review_final_ica.py --subject 08 --plot task --component 0
python3 /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/review_final_ica.py --subject 08 --plot continuous --component 0
python3 /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/review_final_ica.py --subject 08 --plot context --source sub08_run01_210826.fif --category BAD_muscle --center-frontal
```

Source browsers scroll in time and across ICs. Add `--component N` for one IC.
The separate raw browser provides frontal/MISC context. Independent plots can
use `--save /tmp/name.png`; they cannot save production models, cleaned data or
epochs. Reconstructing full fit views reads native data and can take time.
Choose a new diagnostic filename: an existing PNG or JSON sidecar is preserved.

Array 11096 was verified on 2026-11-10: the joined audit checks all 35 models,
180 runs and 2,548 full-fit property figures, with unchanged model hashes and
matched output evidence. Those version-2 books did not yet include default
task-locked/continuous per-IC figures. Array **11197** subsequently completed the
paired views for all 35 participants/180 acquisitions. Its joined audit verifies
2,548 full-fit property panels and 5,096 native continuous/task panels, with
unchanged model hashes, matching sidecars and bound output stamps. The audit
records that large HTML hashes were not reread; their bound size/mtime was
checked. This is software-output readiness, not scientific component acceptance.
See `reports/ica_review_updates/array-11197-audit-20261110-v1/summary.json`.

The commands below document the completed submission. Do not repeat its immutable
label; choose a NEW label only if another report refresh is required. Native-only
mode reuses validated full-fit properties (the user submits any future jobs):

```bash
/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/cluster/deploy_review_to_cluster.sh
```

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
DAD_ICA_REPORT_LABEL=continuous-task-20261110-v1 \
  sbatch --array=4-14,16-39 meg/cluster/submit_final_ica_reports.sh
```

Default resources are 4 CPUs/24 GB per participant; scheduler resources bound
concurrency. Override `--cpus-per-task`, `--mem`, or an array cap when required.
Default `DAD_ICA_REPORT_MODE=native` uses `--native-only`, validates every cached
fit picture's stamp and full-sample sidecar, and adds both native plot types.
Use `DAD_ICA_REPORT_MODE=full` only when all fit properties also need rebuilding.
The wrapper runs no fitting/application/epoching command and writes no decisions.
Join the array before cohort review. Metadata-only HTML updates also support
parallel `review_final_ica.py --refresh --subjects ... --label NEW --jobs N`.
Earlier reports are archived. Existing final-model decisions cause refusal to
replace their bound evidence; use independent read-only plots for those cases.
