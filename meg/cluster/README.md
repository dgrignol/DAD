# DAD MEG Pipeline on DIPSEN HPC

This runbook covers DAD's SLURM submission wrappers on the DIPSEN cluster. The
wrappers read raw MEG data from
`/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data` and write
derivatives under `/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives`.
The raw tree is read-only and must never be used as an output location.

The current analysis cohort is **P04-P14 and P16-P39 (35 participants)**.
P15 is excluded by the explicit 2026-11-07 coverage decision, with data/provenance
retained; P01-P03 are excluded because of experiment-version differences. See the
[cohort review guide](../docs/COHORT_PREPROCESSING_REVIEW_GUIDE.md) for each
stage's purpose, expected derivatives, inspection checklist, and review gate.

Stage numbers refer to workflow/review gates. Numeric script prefixes and report
names remain legacy compatibility identifiers: Stage 1 uses script `02`, Stage 4
uses `03_A`-`03_C`, and Stage 5 uses `04`; Stage 3 is visual timing. See the
[stage-to-script crosswalk](../docs/PIPELINE_SPEC_basic.md#stage-interfaces-and-implementation-status)
for selectors and acceptance status. Keep existing commands, filenames, tuning
keys and lineage identifiers unchanged.

## Current epoch/export retry (2026-11-17)

The [short-fragment retry](../docs/SHORT_FRAGMENT_RETRY_20261117.md) fixes the
millisecond-interval line-removal failure and provides the current cluster
command. It leaves <=8-sample islands unfiltered, with zero overlap with the
accepted annotation-surviving epoch footprint. Use the new processing-run label;
P10 remains pending its separate rejection diagnostic. P37's v1 outputs remain.

## Current executable route (2026-11-12)

The user approved proceeding without ICA rejection: accepted annotated continuous
recordings -> 50/100 Hz spectrum-fit in memory -> epochs/export. Future reviewed
ICA applies immediately before that same line removal. No permanent continuous
ICA copies are required. Use a distinct processing-run label and the
[epoching handoff](../docs/INLINE_EPOCHING_HANDOFF_20261112.md) for the full-cohort
read-only preflight, production array, QC and ERF/dRSA commands. The user submits
all jobs. The older continuous-ICA application commands below are not the current
next step. Clear the two downstream run environment variables before upstream
review/recording tools.

## Earlier Stage 4 handoff: current masks, then final ICA

Use [STAGE4_CURRENT_MASK_HANDOFF.md](../docs/STAGE4_CURRENT_MASK_HANDOFF.md)
for exact deployment, Stage 4B adoption/reconciliation and dependent final-fit
commands for `--array=4-14,16-39`. The explicit 2026-11-07 policy defers exhaustive
manual interval review. `submit_stage4b_adoption.sh` adopts validated current masks
with manual precedence and separate automatic-adoption provenance, then joins
all output audits before processing acceptance. No manual SAVE is manufactured.
Accepted SSS is reused; this pass adds no EDF/SQUID proposals. This active handoff
supersedes historical all-run manual-review and SQUID prerequisites below.

The final fit still uses 1-40 Hz, separate 200 Hz fit copies, FastICA,
`n_components=None`, seed 97 and all operative BAD categories excluded. Final
continuous processing has spectrum-fit 50/100 Hz removal, including P38 run 1,
and no additional high-pass/low-pass. Inspect the final-model component navigator;
component decisions, application QC, P10 slow-contamination review and cleaned
output acceptance remain later explicit work. No Stage 5 job is authorized here.
Historical cohort commands elsewhere in this runbook describe earlier work;
current Stage 4/5 defaults exclude P15 and reject an explicit P15 task.

## Parallel participant execution

**User instruction, 2026-11-03:** run all independent participants in parallel
whenever the operation permits it. Submit the whole eligible cohort together;
do not silently serialize participant processing. Use disjoint participant
processes or scheduler array tasks. Shared cohort publication and true prerequisite
barriers run only after the participant tasks finish. Manual work, scientific
validation, review gates and immutable output protection retain their scope.

Participant preprocessing arrays now default to `--array=4-39`, without the
historical `%4` concurrency throttle. The scheduler still enforces CPU, memory,
account and available-node limits. Explicit throttles remain available when
required by a concrete resource limit; record the reason instead of silently
using a serial cohort loop. Dependencies within a participant, interactive human
review, licensed-resource constraints and shared aggregate writes may require
ordering or a lower worker limit.

The choice importer, annotation counter and embedded HTML renderer now use
separate spawned participant processes. The triage and retention-report wrappers
request **36 CPUs / 64 GB** by default, allowing all 36 participants to run
concurrently once scheduled. Each child uses one numerical-library thread;
cohort plots/manifests have one final writer. The importer completes whole-table
validation AND a fresh-input barrier before starting any participant save worker.
Runs within one participant remain ordered to protect its shared selection CSV.
Use `--jobs N` for direct Python calls, or `DAD_PARTICIPANT_WORKERS=N` in wrappers.
Workers cannot exceed the allocated CPU count. For a deliberately smaller
allocation, for example:

```bash
DAD_PARTICIPANT_WORKERS=8 DAD_TRIAGE_LABEL=NEW_SNAPSHOT DAD_TRIAGE_RECOUNT=1 \
sbatch --cpus-per-task=8 --mem=16G meg/cluster/submit_annotation_retention_report.sh
```

Default worker counts follow the ACTUAL allocation, including older queued jobs.
Changing code/default requests does not resize or restart an existing job. Job
10590 keeps its original sequential process and allocation; completed masks and
manual work remain intact. New jobs print worker counts and per-run progress.
Parallelism removes the serial cohort bottleneck; shared storage and scheduler
availability still determine achieved throughput.

## Access and deployment

The tested SSH login is:

```bash
ssh damiano.grignolio@hpc.dipsen.unitn.it
```

The remote DAD project path remains
`/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD`. From the local
project directory, deploy the review payload with:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
meg/cluster/deploy_review_to_cluster.sh
```

The helper defaults to the tested host and remote path, preserves the
`DAD/meg/...` layout, and refuses deployment into the raw `Attention/data`
tree. Raw data and derivative data stay in their existing shared locations.

DIPSEN uses SLURM. Connect only to the login host; request compute resources
with `sbatch` (or `srun`/`salloc` for short debugging). Direct SSH to compute
nodes is disabled. The default CPU partition is `all-compute-cpu`; the DAD jobs
do not need the GPU partition. SLURM's normal queue has a maximum wall time of
three days. See the [DIPSEN HPC guide](https://hpc.dipsen.unitn.it/) for account
and cluster details.

## Stage 2: prepare trial accounting

Stage 1 is accepted. Stage 2 prepares trial decisions and recording coverage;
it does not create epochs or start cleaning. Deploy the updated payload with
the command above, then submit from the remote DAD root:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
stage2_job=$(sbatch --parsable --array=4-39 meg/cluster/submit_preprocessing_stage2.sh)
echo "Stage 2 participant job: $stage2_job"
sbatch --dependency="afterok:${stage2_job}" meg/cluster/submit_preprocessing_stage2_summary.sh
```

The dependent summary requires all 36 participant jobs to succeed. Inspect
`reports/group/trial_preparation/stage2_cohort_summary.md` under the derivative
root first, then the participant ledgers/reports. The
[Stage 2 handoff](../docs/COHORT_PREPROCESSING_REVIEW_GUIDE.md#stage-2-trial-preparation)
documents every output and the independent-review checklist.

Existing Stage 2 outputs are protected by default. Explicit
`DAD_STAGE2_OVERWRITE=1` replaces only selected Stage 2 bundles; there is no
implicit resume. Use an explicit reduced array for any deliberate regeneration.
`DAD_STAGE2_VALIDATE_ONLY=1` performs subject accounting without derivative
writes. These wrappers call `trial_preparation.py` directly, preserving Stage 1
tables, inventory, source manifest, and execution provenance. They never submit
another job internally.

## Stage 3: accepted visual timing

Stages 1 and 2 are scientifically accepted. Stage 3 `s3-v002` was accepted on
2026-10-21 for all 36 participants with revised eligibility, P24 donor filtering
and the measured-median onset method. Thirty-six separate accepted review
records under `/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/visual_timing/s3-v002/reviews/`
bind exact participant manifest hashes, application hashes, eligible assignment
bases and mandatory review scopes. The accepted Stage 2 ledger and original
`s3-v001` bundles remain unchanged. Do not rerun Stage 3 measurement or the
metadata application for this accepted handoff.

The cohort has 23,963 policy-eligible trials (21,890 original + 2,079 newly
eligible type-2 catches - six newly excluded original trials), after all eight
exact sequence decisions. Five P37 MAT outputs in the unavailable recording
tail remain unrecorded. The 2,079 new catch onsets are in-bounds without
trial-end overlap: 1,718 measured and 361 estimated. The ordinary unweighted
medians are 19 ms for codes 71/81/91 from 22,120 unique QC-usable measured
events. P24's five donor events are excluded; counts decrease 49→48 for code 81
and 90→89 for codes 111/112/114/115, while affected medians remain 19 ms.
Codes 116/117 have six donors each and no estimate. P33's 601 remaining
original-eligible measured onsets were accepted; no general sequence cutoff was
adopted. P29 retains 37 eligible measured events in its partial half-block;
P38 retains 16 eligible measured events in its short first block, with the 13 ms
edge event-specific. Code 118 prompt timing remains deferred behavioral work.

See the [Stage 3 contract](../docs/STAGE3_VISUAL_TIMING.md) and [formal
acceptance](../docs/PREPROCESSING_DECISIONS.md#stage-3-gate-outcome) for exact
keys, donor methods and limitations. Stage 4 cleaning is independent. Stage 5
still requires accepted cleaning lineage and ordinary production prerequisites.

## Prepare the shared environment

The wrappers expect `micromamba` at `$HOME/bin/micromamba` and an environment
named `dad-meg`, unless `MICROMAMBA_BIN` or `DAD_MAMBA_ENV` is set before
submission. If the environment is not present on DIPSEN, create it from a
short interactive compute allocation rather than running package installs on
the login node:

```bash
srun --qos=interactive --partition=all-compute-cpu --nodes=1 --ntasks=1 \
  --cpus-per-task=2 --mem=4G --time=02:00:00 --pty bash -l
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$("$HOME/bin/micromamba" shell hook -s bash)"
micromamba create -y -f meg/cluster/environment.yml
micromamba activate dad-meg
```

Exit the interactive shell after setup. The wrappers activate the environment
on their allocated compute nodes. The environment path is in shared home, so
it is visible to compute jobs. If your micromamba executable is elsewhere,
submit with an explicit exported path, for example:

```bash
MICROMAMBA_BIN=/absolute/path/to/micromamba sbatch --array=4-39 \
  meg/cluster/submit_preprocessing_stage1.sh
```

## SLURM job behavior

Run submissions from the remote DAD project root. Create `cluster_logs` before
the first `sbatch`: SLURM opens its output file before the wrapper starts.
Array tasks are throttled to at most four running tasks by default. The
wrappers request CPU, memory, and wall time directly in `#SBATCH` lines;
command-line options can override them. The wall times are conservative
starting estimates; if a task reaches its limit, increase it with `--time` up
to DIPSEN's three-day normal-QOS maximum. The wrappers use the default CPU
partition.

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs

# Check queued/running jobs and inspect a submitted job.
squeue --me
scontrol show job <job_id>

# Cancel a job only if needed.
scancel <job_id>
```

Each task writes stdout to a `.log` file and stderr to a matching `.err` file.
Array names use
`cluster_logs/<job-name>.<array-job-id>_<task-id>.<log-or-err>`; single-job
names use `cluster_logs/<job-name>.<job-id>.<log-or-err>`. SLURM expands these
names before the wrapper runs.

Environment settings can be passed by exporting them or by placing them on
the `sbatch` command. For example:

```bash
DAD_REBUILD=1 sbatch --array=4-39 \
  meg/cluster/submit_preprocessing_review.sh
```

The wrappers use `--export=ALL`, so shell variables are inherited by the job.
Quote values containing spaces; comma-separated values such as analysis lists
are ordinary shell environment values under SLURM.

## Stage 1: source, trigger, and acquisition review

Pilot derivatives are preserved under `derivatives/pilot`. The first P08-P39
Stage-1 attempt is preserved separately under
`derivatives/previous_runs/stage1_2026-10-15`. Both archives have verified move
manifests. The fresh array covers P04-P39.

This gate inventories sources and session notes, extracts STI101 triggers,
and writes both canonical alignment outputs and the full deep report with
anomaly/transition plots. It does not apply photodiode timing corrections or
run Maxwell filtering, muscle review, or ICA. P24's adapter is no longer called:
visual-timing measurements and decisions are handled in the separate cohort
Stage 3 gate, accepted on 2026-10-21.

The P04-P39 digital Stage-1 outputs are already complete. Deploying this scope
change and removing only the two preliminary active P24 timing files requires
**no job resubmission**. The existing digital outputs remain valid.

For a future fresh cohort execution, the command remains:

```bash
sbatch --array=4-39 meg/cluster/submit_preprocessing_stage1.sh
```

If P24's digital outputs later need regeneration, P24 alone is sufficient:

```bash
sbatch --array=24 meg/cluster/submit_preprocessing_stage1.sh
```

That optional rerun rewrites P24's inventory/provenance, canonical trigger
files/reports, execution log, and deep plots. Trigger extraction is not
protected by `--resume`, and it writes existing outputs without requiring
`--overwrite`. Do not run it when other existing derivatives must be preserved.

Review each participant's inventory, source manifest, event counts, alignment
reports, and experimenter session notes together with the recordings. Treat
the notes as review clues, not as ground truth. Preserve P17's 0.1 Hz
acquisition high-pass and lack of photodiode, P24's early acquisition
differences, and P38's six blocks with a partial first block as context from
the shared `meg/docs/PREPROCESSING_DECISIONS.md` record. Retain the guide's
Stage-1 acceptance checklist; defer P24 calibration and all timing fallback
decisions to the cohort-wide visual-timing gate.

P37's reviewed acquisition/block map prevents its truncated first recording
from shifting later acquisitions' alignment. Its missing logged suffix remains
`REPORTED_ONLY`, with `recording_tail_unavailable` context in the canonical
and deep CSVs. Check the final incomplete trial and the separate earlier
missing gaze-break marker before accepting the usable portion.

## Preliminary annotation-retention triage

User-approved 2026-10-29: obtain a preliminary retention ranking before detailed
manual annotation review throughout the cohort. Preserve P04/P05 and P06 runs
1-3; P06 runs 4-5 remain unreviewed and are included in the automatic workflow.
Newly saved manual work is detected dynamically and takes precedence. No EDF or
SQUID candidates are used here, and this workflow rebuilds no accepted SSS.

### Recommended: one offline report per participant

Approved 2026-10-30 to remove the interactive loading delay. Build maps and
navigable traces in advance, then open ONE participant HTML in your normal web
browser. Use its run selector/previous-next buttons to move between acquisitions.
Viewing needs no Python process, FIF read, network asset, server or new package.
The older interactive viewer remains available below.

The real P14 pilot is under:

```text
/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/annotation_triage/component_reports/ocular-reports-20261030-v1/index.html
```

Reports list existing manual work and saved choices without changing them.
Only runs still pending at export receive new trace/map data. New component
choices saved during generation do not invalidate a read-only trace export.
Scientific source/model/mask changes still require reconciliation. No automatic
classification or scientific acceptance occurs during export.

Generate further participant books locally (background batch, not interactive):

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/build_ocular_component_reports.py \
  --subjects {15..39} --label ocular-reports-20261030-v1 --resume
```

The preferred large batch runs on cluster compute nodes after deploying the
updated payload. From the cluster DAD root, the USER submits:

```bash
DAD_COMPONENT_REPORT_LABEL=ocular-reports-20261030-v1 \
  sbatch --array=15-39 meg/cluster/submit_ocular_component_reports.sh
```

Each task makes one participant book and updates the shared index under a lock.
The index page reloads every 30 seconds; participant books never auto-reload.
Publication normally uses atomic replacement. If shared-storage permissions
deny replacing an existing Mac/cluster-created index, a locked in-place fallback
updates only index.html and the generated blank template, preserving their ACLs.
Immutable participant reports/manifests/assets never use this fallback.
No nested jobs are submitted. A common label is mandatory. Existing complete
books are reused only with compatible inputs/viewer/settings; use a new label
when those change. Internal per-run assets support interrupted generation but
are not separate reports for the user to open. The local P14 pilot need not be
regenerated on the cluster; start the remaining array at 15.

If reports exist but the navigation index is stale, repair navigation locally
without rerunning report generation:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/build_ocular_component_reports.py --subjects 14 --label ocular-reports-20261030-v1 --index-only
```

The 2026-10-30 batch 10384 hit a shared-storage permission error replacing the
Mac-created index after participant books had already been published. Navigation
was repaired locally from the completed manifests; report/FIF regeneration and
job resubmission were unnecessary. Already-running older scripts keep their
loaded publication code; the updated publisher applies to newly started runs.

Every book has magnetometer/gradiometer maps, zero-based component labels, a
time slider, horizontal drag, time jump, zoom and gain controls. Arrow keys
scroll time when no form field is focused. Gray bins contain unscored native
samples; they are never joined across a gap. The compact display retains the
largest absolute native sample per bin with its exact time offset, at a default
200 bins/s. Original native sources are unchanged and remain detector inputs.
Display reduction/normalization is diagnostic only, not scientific resampling.

Keep decisions PER RUN. Copy the report's note row, or fill a COPY of its
decision template and send it in chat:

```text
P15 R02 | blink: [zero-based indices or NONE] | eye: [indices or NONE] | notes: ...
```

Multiple indices and a component selected for both groups are allowed. NONE
disables one detector, never establishes no artifacts. DEFER defers a run.
Reports do not save notes or scientific decisions. The generated templates are
blank references; keep your completed notes separately so a batch index refresh
cannot overwrite them. CSV imports require explicit choices in both groups;
blank rows stay unfinished. Quote comma-separated indices in CSV.

After converting the human notes to the template CSV, validate the ENTIRE table:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/import_ocular_component_choices.py \
  --reports-dir /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/annotation_triage/component_reports/ocular-reports-20261030-v1 \
  --decisions /absolute/path/to/completed_notes.csv
```

Only after successful validation, repeat with --apply. The importer checks
the exact report/manifest, current native evidence/model, fitted component
bounds, both groups, duplicates and existing manual/choice precedence before
writing anything. It permits cluster-to-macOS path relocation while revalidating
the current effective fit/source/mask settings. Identical existing choices are
kept verbatim; conflicting choices require explicit reconciliation, not a bulk
overwrite. Per-source atomic saves retain the existing component-choice format.
Then run annotation_triage.py for the selected participants to create provisional
annotations and retention counts, as described below. EDF/SQUID proposals and
production 4A/4B/4C/epoching gates retain their existing scope.

For the completed 2026-11-03 notes, a converted 131-row CSV is saved separately
as `annotation_triage/review_notes/ocular_component_choices_completed_20261103.csv`.
The original text notes remain unchanged, including the P35 run-6 uncertainty.
Run the following ONE compute job (not an array) from the cluster DAD root:

```bash
mkdir -p cluster_logs
DAD_COMPONENT_REPORT_LABEL=ocular-reports-20261030-v1 \
DAD_TRIAGE_DECISIONS=/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/annotation_triage/review_notes/ocular_component_choices_completed_20261103.csv \
DAD_TRIAGE_LABEL=retention-20261103-v1 \
sbatch meg/cluster/submit_annotation_triage.sh
```

This validates the WHOLE decision table before importing choices, then generates
provisional masks and a single P04-P39 retention snapshot. Earlier 41 component
selections and 13 complete manual reviews also contribute; current manual work
retains precedence. A failed import stops the job before counting. Logs use
`cluster_logs/DAD_annotation_triage.JOBID.log` and `.err`. The resulting
`annotation_triage/reports/retention-20261103-v1/report.md` and
`review_priority.csv` identify the first participants for targeted review.
These are annotation-only window counts before final ICA/amplitude rejection.
Use a new retention label if an earlier attempt already published that folder.
Only the user submits the job. The wrapper performs no production acceptance,
ICA fitting/removal or epoch creation.

### Embedded HTML retention reports

The original counter writes identity/condition/run CSV tables, a text summary
and one cohort condition heatmap. For the richer visual report requested on
2026-11-03, `render_annotation_triage_report.py` reads those saved tables only.
It reconciles counts with the trial ledger and creates self-contained HTML with
every figure embedded, one cohort overview and one book per participant:

- Three identity barplots per participant, identities 1-40 on x and trial counts
  on y. Blue survivors plus pale-red loss equal the eligible count; all three
  conditions share a y scale. Exact counts, true zero cells and denominators are
  retained.
- Cohort surviving-count and retention-rate heatmaps in review-priority order;
  all 120 condition/identity cells, run retention, overlapping rejection causes,
  exclusive trial-accounting tables and ocular scoring coverage.
- Optional previous-snapshot points on identity bars and a cohort change
  heatmap/CSV. Comparison requires unchanged source/trial identities, accepted
  eligibility and native epoch coordinates. Input-error participants stay
  unranked, with diagnostic error pages rather than zero-survivor plots.

The HTML files need no external image files, JavaScript, server or network.
PNG exports, condition summaries and SHA-256 provenance are also saved. Original
snapshots and annotations remain unchanged. To prepare the initial books after
job 10590 succeeds, the user can submit this ONE dependent job now:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
DAD_TRIAGE_LABEL=retention-20261103-v1 \
sbatch --dependency=afterok:10590 meg/cluster/submit_annotation_retention_report.sh
```

Open the mounted output
`/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/annotation_triage/visual_reports/retention-20261103-v1/index.html`.
The job log is `cluster_logs/DAD_retention_report.JOBID.log`.

After saving targeted manual annotation reviews, count the CURRENT full cohort
again with a NEW snapshot label, then render and compare automatically. This
does not repeat the completed notes import:

```bash
DAD_TRIAGE_RECOUNT=1 \
DAD_TRIAGE_LABEL=retention-20261103-v2 \
DAD_TRIAGE_COMPARE_TO=retention-20261103-v1 \
sbatch meg/cluster/submit_annotation_retention_report.sh
```

Re-rendering an old snapshot alone cannot reflect newer annotation edits. To
change only the presentation, keep `DAD_TRIAGE_RECOUNT=0` (default) and select
a new `DAD_TRIAGE_VISUAL_LABEL`; no MEG data or detectors are read. Local plotting
after a count job finishes is also possible:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/render_annotation_triage_report.py \
  --report-dir /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/annotation_triage/reports/retention-20261103-v1 \
  --label retention-20261103-v1
```

Every count/visual snapshot is immutable. Choose a new corresponding label if
it already exists. Counts remain provisional and precede final ICA/amplitude
rejection; no participant-exclusion threshold is introduced.

### Existing interactive component viewer

Run the bulk viewer locally from the workspace DAD root. Its default cohort is
P04-P39; existing valid manual reviews and compatible component choices are
skipped. The first command is strictly read-only and opens no GUI:

```bash
cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/review_ocular_components.py --check-only
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/review_ocular_components.py
```

For each recording inspect maps and scrollable native cached sources, close the
source browser, enter blink and eye-movement components separately, then SAVE.
`0`/`ICA000` is the first component; comma-separated indices/names and shared
components are allowed. NONE disables just that detector; DEFER revisits the
recording later. QUIT retains earlier saves. Optional observations are recorded.
Choices are saved immediately in model/source-bound JSONs and a per-subject
`component_selections.csv` under `derivatives/annotation_triage/component_selections/`.
No whole-source interval review or final ICA exclusion is implied.

Use `--subjects 06` or `--source sub06_run04_260526.fif` to narrow the queue.
`--revisit` explicitly replaces prior choices, including stale choices after a
new cache; completed manual work still cannot be replaced. Cache problems print
targeted cluster preparation advice and require reconciliation. Preparation
continues to use the existing `DAD_BLINK_AID_ONLY=1` mode; no fitting occurs locally.

After choices, generate a new preliminary report locally:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/annotation_triage.py
```

This reads accepted Stage 2/3 eligibility/onsets and current annotations without
creating epochs. Existing source BADs and reviewed manual replacements are
preserved; otherwise existing detectors propose blink/eye intervals alongside
automatic muscle. Settings stay provisional. Missing selections and unscored
periods are explicit, not evidence of clean data. Outputs include:

- `report.md` / `report.json`, `review_priority.csv`, `condition_retention.png`;
- `trial_retention.csv` with source/trial keys, all BAD overlaps and exact mask paths;
- `condition_retention.csv`, `run_condition_retention.csv`, `identity_retention.csv`;
- `overlapping_rejection_reasons.csv` and `errors.csv` (reason counts can overlap).

All reports live in new folders under `derivatives/annotation_triage/reports/`.
An optional `--label NAME` selects a new immutable folder. Rank by the lowest
surviving condition count, then identity gaps; compare retention percentages and
ocular coverage too. Participants with input errors are unranked; their partial
tables are diagnostic. Zero survivors is a valid result. These counts precede
final ICA/line treatment/amplitude rejection and do not decide participant replacement.

Each trial points to its content-hashed provisional JSON. For targeted review,
use its exact path and source (example placeholders must be replaced):

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python -B meg/code/preprocessing_basic/manual_eye_annotations.py \
  --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-06.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-06/ses-01/03_tuning_choices.yaml \
  --source sub06_run04_260526.fif --basis misc_and_frontal \
  --triage-proposals /absolute/path/from/trial_retention.csv.json
```

Inspect the whole source, including unflagged/unscored periods; edit and SAVE
normally. Imported masks never overwrite existing manual work: reopen that work
with its ordinary command and `--overwrite`. A triage-based complete manual SAVE
records source-specific EDF opt-out and detector/deletion provenance. Reopening
retains it, and later 4B rebuilds omit that source's EDF guidance. Earlier reviews
and independent 4B/4C acceptance remain intact. Do not combine --triage-proposals
with --preliminary-ica or --recover-annotations.

Rerun `annotation_triage.py --subjects 06` for an affected participant, or rerun
the default cohort command for an updated comparative ranking. Compatible
detector masks are reused after current input validation. Older snapshots stay
unchanged. Deployment includes both new tools; a cluster report invocation uses
`python .../annotation_triage.py --site-config meg/params/pipeline_cluster.yaml`
inside an ordinary compute allocation. No new scheduler wrapper or job is
submitted by these tools. Production epoching still requires the accepted
cleaning lineage described below.

## Stage 4 cleaning: separate 4A, 4B and 4C reviews

Stages 1-3 are accepted. The user submits these jobs; user and Codex inspect
independently, reconcile and record scientific decisions. The approved 2026-10-22
implementation/settings do not accept current-cohort cleaning outputs.

Calibration files now use `meg/meg_calibration_files/fwmegmaxfilter` in the
workspace deployment copy and shared/cluster DAD tree. Contents are unchanged.
Deploy the updated payload with the existing deployment helper before running.
It preserves the shared decisions record and never deploys raw/derivative data.

### 4A: per-run channels and SSS

From the remote DAD root, the review wrapper now defaults to 4A only:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
mkdir -p cluster_logs
sbatch --array=4-39 meg/cluster/submit_preprocessing_review.sh
```

Inspect `reports/sub-XX/ses-01/03_A_MaxFilt_review_report.json/html`, per-run bad
channel score NPZ files and matching raw/SSS PSD PNGs. Review flagged/unflagged
examples across each acquisition. Retain DAD central_run and SSS only. Record
any explicit manual additions/removals in 03_tuning_choices.yaml; regenerate 4A
with `DAD_REBUILD=1` for changed choices, then reconcile the new evidence.

The primary HTML report embeds its figures and includes processing explanations,
figure limitations and review criteria. It can be copied/opened independently of
the cluster figure paths. To repair reports from jobs run before this reporting
fix, render the existing JSON and PNGs only; do not resubmit preprocessing:

```bash
python meg/code/preprocessing_basic/render_stage03_mne_reports.py \
  --derivatives-root /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives \
  --raw-root /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/data \
  --subjects {04..39} --stages maxwell_filter --processing-pass review \
  --primary-report --overwrite
```

This command changes HTML only. It opens no FIFs and preserves JSON evidence,
channel decisions, review receipts, annotations and models. On macOS, use the
corresponding `/Volumes/MORWUR/Projects/DAMIANO/Attention/` roots. The same renderer
is called automatically by future 4A/4B/4C processing; omit `--primary-report` to
write a separate `_mne_report.html` instead.

To replace the old detection-average PSD with per-run noisy/flat score heatmaps
and highlighted raw traces, the user submits a plot-only refresh:

```bash
DAD_MAXWELL_QC_ONLY=1 sbatch --array=4-39 meg/cluster/submit_preprocessing_review.sh
```

This reuses saved scores and decisions. Each recording has an initial raw PSD
with one curve per sensor and sensor-location insets, then heatmaps and one
reproducible random 20-second raw trace window. The last PSD figure compares
individual retained channels before/after SSS, excluding applied bads on both
sides with matching y scales. PSDs read the same first up-to-120-second raw/SSS
windows (1-60 Hz, identical Welch settings). It writes PNGs,
`03_A_MaxFilt_review_plot_metadata.json` and HTML only.
Detected channels in this exact run are red `[D]`; applied-only channels from
harmonization/manual choices are amber `[A]`. These labels do not change bad
channel decisions. Heatmaps distinguish missing scores and threshold crossings;
a random raw window can miss an intermittent fault. The exact native sample
interval is recorded, and all MEG channels share physical scales by sensor type.
Raw overview PSDs highlight detected bads with named red dashed curves and
applied-only bads in amber dash-dot curves. Plotting does not filter, resample
or modify data. Ordinary 4A processing now
produces these plots automatically. Optional `DAD_QC_WINDOW_SEC` and `DAD_QC_SEED`
change display sampling only; defaults are 20 and 97.

Local equivalent, useful for one participant's layout review:

```bash
python meg/code/preprocessing_basic/maxwell_qc.py \
  --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
  --raw-root /Volumes/MORWUR/Projects/DAMIANO/Attention/data \
  --subjects 04 --window-sec 20 --seed 97 --overwrite
```

After reconciliation, record each participant's exact 4A acceptance locally.
Example for P04 (the reason must describe the actual reconciled finding):

```bash
python meg/code/preprocessing_basic/review_cleaning.py \
  --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/03_tuning_choices.yaml \
  --stage 4A --reviewer 'User and Codex' --reason 'Reconciled channels and per-run SSS evidence'
```

Do not run acceptance commands before inspection. Formal acceptance belongs in
the shared PREPROCESSING_DECISIONS.md. Configuration/report/source changes make
a receipt stale. Stage 4 execution keeps Stage 1 manifests/inventories/logs intact
apart from additive derivative lineage; Stage 4 execution snapshots have their
own stage4_ prefix.

### Stage 4A SQUID candidates before SSS

Implementation approved 2026-10-28. This is part of workflow 4A, using the legacy
`03_A_MaxFilt.py` ID; Stage 3 remains visual timing. Existing accepted automatic
participant-union channels are skipped. Installing/deploying the extension does
not rewrite SSS, manual BAD annotations, tuning receipts or accepted timing.
Preserve P04 all-run and P05 run-1 manual reviews, and finish the currently open
P05 run-2 review normally. No restart/reload is required for that open process.
The commands below are user-submitted procedures, not jobs run during implementation.

After deploying the updated payload, prepare candidates only for selected
participants. For example, from the cluster DAD root:

```bash
DAD_SQUID_PREPARE_ONLY=1 sbatch --array=4%1 meg/cluster/submit_preprocessing_review.sh
```

`DAD_SOURCE=EXACT_SOURCE.fif` optionally limits preparation to one acquisition.
This mode reuses accepted channel lists/scores, reads native MEG and writes
immutable JSON/CSV under `preprocessing/sub-XX/ses-01/squid_jumps/`. It runs no
channel detector, SSS, muscle processing, ICA or annotation/tuning save. Existing
proposals are reused after identity/content verification; detector changes create
new proposals rather than overwriting old decisions. The alternative local
preparation command is the review command below with `--mode prepare` and no
`--source`, or `03_A_MaxFilt.py ... --squid-prepare-only`.

The initial P04 run-1 pilot (job 10078) exposed inflated scores and transitive
merging into a recording-long v1 proposal. Do not SAVE that proposal as BAD.
Detector v2 floors the smoothed derivative MAD with the raw derivative MAD and
bounds channel/cross-channel peak groups. Its new algorithm identity creates
fresh candidates while retaining the v1 JSON/CSV. Deploy v2 before preparing
further sources; installing it changes no reviewed intervals or accepted SSS.

Review one exact source locally, without opening or replacing its eye-review
file. Example for P04; substitute the actual source:

```bash
python meg/code/preprocessing_basic/review_squid_jumps.py \
  --mode review --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-04.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/03_tuning_choices.yaml \
  --source sub04_run01_220526.fif
```

Inspect the listed candidates and nearby unflagged context. New sessions show
the selected candidate's FLAGGED sensors only, with the candidate ID/time in
the title and exact names in the terminal. ALL shows all screened sensors;
VIEW/NEXT/PREV returns to focused review. Native traces use display-only window
mean centering and physical mag/grad scales computed from that displayed window,
so remote drift or channel offsets do not flatten active traces. No filtering or source-sample change occurs;
within-window steps remain visible. Delete/resize `possible_squid_jumps`, or add the exact label
`BAD_squid_jumps`. Mark a channel by clicking its name; this is bad-channel
metadata, not a time annotation. The default channel scope is participant union;
choose run explicitly only with a justified exception and compatible final ICA
bases. Extend windows to actual recovery; a lasting shift may warrant channel
exclusion or a longer interval. A +/-0.5 s proposal is not an established
artifact duration. Inspect listed unscored channels separately.

Overlapping proposal shading is merged for display, so one muscle burst does
not obscure the traces with stacked fills or require many identical deletions.
The original candidate IDs/JSON/CSV and detector sensitivity remain unchanged.
SAVE already merges retained overlapping windows; display consolidation does
not itself confirm or reject any candidate.

VIEW (optionally followed by a candidate number) reopens current edits without
reseeding proposals. SAVE asks for scope, reviewer, observations and explicit
confirmation, then saves immutable review JSON and atomically merges only this
source's `squid_jumps.reviews` entry into the existing tuning file. Reviewed zero
is explicit. DISCARD preserves all previous saved work. Other sources and active
eye saves retain the locked latest-map merge behavior. After changed detector
settings/channel evidence, `--reconcile` loads old decisions as unreviewed native
starting intervals, lists new candidates without overwriting manual edits or
deletions, and requires renewed inspection/SAVE. Original decisions remain intact.

For precise read-only inspection, new review sessions accept `DETAIL 20` after
closing the editable browser. Each flagged sensor gets its own labeled physical
y axes: six seconds of native context, a 300 ms edge close-up with the detector's
median trace, and median-filtered first differences. These short windows never
re-estimate a channel score or infer an artifact cause. An absent lasting step
does not establish good data, and a detector score does not prove a SQUID jump.

An already open review runs its original code and may contain unsaved edits.
Leave it open; do not restart/discard it merely to get the new display. Open an
independent inspector in a SECOND local terminal, for example for candidate 20
(672.999 s, MEG0141 and MEG0241 in P04 run 1):

```bash
MPLBACKEND=QtAgg /opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python \
  meg/code/preprocessing_basic/review_squid_jumps.py \
  --mode inspect --config meg/params/pipeline_default.yaml \
  --subject-config meg/params/sub-04.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/03_tuning_choices.yaml \
  --source sub04_run01_220526.fif --candidate 20
```

`--context` also shows available colocated sensors explicitly marked as context,
not flagged by this event. Inspection reads one existing current proposal and
bounded raw samples, creates no candidate/channel-basis/annotation/tuning file,
and never saves review decisions. Close the inspector and return to the original
review to make edits. `--candidate 20` can also start a NEW editable review at
that candidate. No detector rerun, cluster submission, SSS change or migration
of saved P04/P05 BAD work is needed for these display improvements. Muscle-like
proposals may be dismissed from SQUID review; existing Stage 4B muscle decisions
remain independent and preserved.

Review every configured source before invoking 4A application. For a new
participant without accepted channels, run `03_A_MaxFilt.py ... --interactive`
first: accept the channel policy; it prepares SQUID proposals and stops before
SSS while reviews are missing. This deliberate review-required result is not
an instruction to disable the gate. `squid_jumps.enabled: true` is the new 4A
default; old configurations lacking the section retain their previous behavior.

After review, the user may run 4A with guarded resume to retain compatible SSS:

```bash
sbatch --array=4%1 meg/cluster/submit_preprocessing_review.sh
```

If operative windows/channels changed, resume fails safely. Rebuild only the
explicit affected participant list with `DAD_REBUILD=1`. Changed SSS files and
old 4A report/score/HTML evidence are preserved under the report directory's
`history/squid_jumps/<prior-report-sha256>/` before replacement; unchanged runs
retain their SSS stamps. Missing/stale SQUID reviews stop before SSS. Compatible
accepted SSS and reports remain unchanged for reviewed-zero/unchanged actions;
a separate SQUID summary records those reviews. Prepare/review alone changes
none of these SSS/4B files.

Confirmed windows are skipped by SSS. The MNE-1.10-compatible implementation runs
SSS separately on good stretches and restores them to the original native grid.
Skipped samples remain acquired/untransformed and BAD; recording length,
first_samp and auxiliary channels remain preserved. Reconcile new before/after
evidence and record affected 4A acceptance with the existing review command.
The compatibility exception for final temporal policy never accepts changed
SQUID treatment.

Actual SSS changes invalidate affected 4B bindings/caches and participant ICA.
Preserve previous eye/muscle BAD intervals as candidates for explicit revalidation
against current accepted SSS; do not delete, reset/reseed or auto-approve them.
Unchanged SSS requires no repeat manual review. Stage 4B preserves the current
4A-owned SQUID category once; browser edits to its context are discarded rather
than written into eye snapshots. Revise SQUID treatment through 4A. All severe
BAD masks and Stage 5 epoch rejection honor `BAD_squid_jumps`. This adds no
timing/eligibility rewrite, final high-pass, ICA setting or epoch-window change.

#### Coordinate the pending 4B refresh with SQUID review

The earlier no-high-pass handoff proposed refreshing P04's saved manual work
directly in 4B. With the approved SQUID extension, finish/reconcile the selected
participant's SQUID reviews first, then choose the affected refresh below. This
avoids rebuilding 4B against SSS that may immediately change. Candidate-only
preparation and browser/inspector improvements do not require any SSS/4B rebuild
and never establish a confirmed artifact; muscle-like proposals can be dismissed.

- **Reviewed zero or unchanged operative SSS:** preserve accepted 4A reports,
  source stamps and existing manual reviews. Refresh only the affected 4B
  evidence; use guarded resume when FIF/CSV signatures still match, or the
  explicit 4B rebuild when saved manual masks/thresholds have changed.
- **Confirmed windows or new channel exclusions that change SSS:** preserve
  historical 4A evidence/SSS and saved manual files/receipts, apply/QC/reaccept
  only genuinely changed 4A treatment, then explicitly revalidate the affected
  saved 4B intervals before the 4B refresh. Native-coordinate intervals, edits,
  deletions and reviewed-empty choices remain the starting work, not a reason
  to redo everything. Do not automatically rebind old receipts to changed SSS.
  A participant-union channel addition can change every run; a source-local
  interval does not itself require repeating unchanged sources' manual reviews.

Changing SSS also invalidates affected supplementary caches and the participant's
final ICA input fingerprint. Rebuild/review those models only when required;
candidate/UI-only changes do not invalidate models. Keep 4C fitting and final
application deferred until the respective current 4A/4B/component gates pass.
The final temporal policy still uses no high-pass/low-pass, no FIR margin, and
10 s spectrum-fit at 50/100 Hz on every source including P38 run 1. All BAD
categories, including BAD_squid_jumps, protect fitting and final line removal.

Finish open ocular/muscle reviews and SAVE prompts normally. Do not close/reload
them to get SQUID display changes or schedule a rebuild of their inputs while
they remain open. P04 runs 1-2's older ocular-only receipts preserve completed
ocular work; absence of muscle-review fields is not reviewed-zero muscle.
Reconcile that remaining coverage separately before 4B acceptance.

### 4B: provisional muscle evidence, local ocular review, then annotations

For the current pass use the [policy adoption handoff](../docs/STAGE4_CURRENT_MASK_HANDOFF.md).
The historical manual-entry procedure below remains available for later targeted
review; exhaustive new inspection and fresh EDF/SQUID work are deferred.

Stage 4A closed on 2026-10-23: 36 report/configuration-bound tuning receipts
cover P04-P39 and all 185 runs. No 4A rebuild was needed for that acceptance.
Continue existing 4B review work; separately apply the new SQUID extension above
before changed SSS is used, with affected revalidation as documented.
Neither 4B nor 4C has been accepted. Shared tuning retains P38 run 1's historical
high-pass warning; current spectrum-fit runs on every recording. Calibration
resources and accepted earlier stages are intact.

After accepted 4A, generate muscle evidence:

```bash
DAD_STAGES=artifact_annotation sbatch --array=4-39 meg/cluster/submit_preprocessing_review.sh
```

Unreviewed eyes produce provisional review_required reports. Muscle starting
settings are threshold 10, mag, 110-140 Hz and 0.2 s minimum good gap. Inspect
flagged and neighbouring unflagged intervals before accepting thresholds.

#### EDF-assisted guidance (approved 2026-10-23)

The existing exploration is reused for trial-local placement, without repeating
its ICA, peak/lag scans or synchronization analyses. `artifact_annotation.edf_blink_candidates`
is enabled in the existing YAML. Before the first updated job, add its one pinned
reader dependency to the cluster environment (inside your usual allocation):

```bash
micromamba run -n dad-meg python -m pip install --only-binary=:all: eyelinkio==0.3.0
micromamba run -n dad-meg python -c 'from eyelinkio.edf import _edf2py; print("EDFAPI reader ready")'
```

The updated environment/requirements files also include this reader. It bundles
EDFAPI; its small version-specific adapter preserves original timestamps and
recording START/STOP boundaries that the public reader compresses. No global
offset, drift fit, Stage 3 shift, or fixed physiological correction is applied.
The 20 ms within-trial marker tolerance is a mismatch guard, never a correction
or epoch rule. Trial identity/onset/frame sequences must match; unique neighbouring
contexts can recover clear matches after missing trials. Ambiguous/mismatched
trials, unknown recording bounds and out-of-recording events are flagged, not
silently assigned to another trial. Filename block numbers only narrow input
files; accepted source associations and trial markers establish correspondence.

Use the ordinary 4B command above with guarded resume: it reuses unchanged muscle
FIF/CSV/scores and adds exact-source sidecars under `preprocessing/sub-XX/ses-01/edf_blink_candidates/`:

- `*-annot.fif`: `EDF_blink_candidate` / `EDF_blink_long_candidate`, neither BAD.
- `*_blinks.csv`: original eye events/endpoints, mapped samples and flags.
- `*_mapping.json`: input fingerprints, matched/unmatched trials and recorded coverage.

Overlapping binocular bouts merge only within the same trial. Durations over
0.5 s are retained with a long-event review flag; they are not automatically
rejected or converted to artifacts. Missing/unmatched EDF coverage is not evidence
of no artifacts. Candidate labels affect neither ICA fitting nor epoch rejection.
New/changed guidance requires complete-source ocular review before 4B acceptance.

Locally inspect each exact acquisition, using usable MISC and frontal MEG.
Example (repeat with every source, including partial/additional acquisitions):

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/dRSA/bin/python meg/code/preprocessing_basic/manual_eye_annotations.py \
  --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/03_tuning_choices.yaml \
  --source sub04_run01_220526.fif --basis edf_misc_and_frontal
```

The local dRSA environment already used by the exploration has the reader and Qt
browser. Candidates and available EDF vertical-gaze/pupil traces load automatically.
EDF traces use the same trial anchors on the native MEG grid, with nearest measured
samples only: no interpolation across recording/missing-sample gaps. NaNs remain
visible gaps; `EDF_coverage` (0/1) and `EDF_no_guidance` identify unavailable periods.
Coverage=1 means an original sample time is represented, not that every eye value
is valid. Each `*_norm` trace is median-centred and scaled by its 95th-percentile
absolute deviation for display; Y means vertical gaze and P means pupil, with
L/R naming the eye. Unavailable eye traces are omitted. Browser units for these
extra EDF MISC channels are arbitrary units (AU), not microvolts; original MEG/MISC
channels keep their identities/units. All extra channels exist only in the disposable browser. Use `--no-edf-traces`
to show candidate labels alone. This comparison supports annotation, not a global
alignment claim. Source/proposal changes stop stale reuse; the reports expose limits.
An explicit recorded display gain keeps normalized traces and the coverage step
visible even when long gaps would make automatic scaling zero. It changes only
synthetic EDF channels, never original MISC samples. Compare waveform/timing
context rather than interpreting these display amplitudes as calibrated analog units.

Delete, move, resize or create intervals in the raw browser. SAVE explicitly
confirms the entire source and all retained candidates, including long ones,
have been reviewed, then converts retained candidates to `BAD_blink`. Delete
uncertain candidates first; add reviewed BAD intervals where EDF provides no
guidance. After closing the browser, type exactly `SAVE` or `DISCARDALL`;
blank answers and typos such as `SAFE` repeat the prompt. `DISCARDALL` discards
only this session's edits, leaving any previously saved file/receipt unchanged.
Deleting all candidates and explicitly saving is a valid reviewed-empty outcome.
Reopening unchanged reviewed guidance
does not resurrect deleted candidates. Original proposals remain separate from
the edited manual-eye file and its fingerprint-bound receipt. After manual SAVE,
rebuild 4B below to apply the newly reviewed BAD intervals; resume must reject
changed manual decisions. 4A never needs rebuilding for this addition.

#### Cached frontal ICA and automatic blink proposals (approved 2026-10-26)

Deploy the updated payload first. From the cluster DAD directory, pilot one
exact recording without rerunning SSS, muscle detection or final ICA:

```bash
DAD_BLINK_AID_ONLY=1 DAD_SOURCE=sub04_run01_220526.fif sbatch --array=4 meg/cluster/submit_preprocessing_review.sh
```

After inspecting the pilot, prepare P04-P39; compatible caches resume:

```bash
DAD_BLINK_AID_ONLY=1 sbatch --array=4-39 meg/cluster/submit_preprocessing_review.sh
```

The standalone `blink_annotation.py` reads accepted 4A and current provisional
4B muscle evidence. It uses ALL good MNE Left-frontal/Right-frontal mag/grad
channels, up to ten components (limited by numerical rank), and existing
1-40 Hz/200 Hz fit-copy settings. It excludes known muscle/severe BAD intervals;
routine ocular examples remain available for this auxiliary fit only. No new
manual severe-artifact pass is required. Severe unmarked artifacts can still
spoil the fit; convergence/rank/coverage are recorded, and no usable component
is a valid reason to continue manually or prepare a targeted refit. Ten PCA
dimensions do not guarantee a clean blink component. Final 4C remains separate
with n_components=None and every reviewed BAD category excluded.

Each exact source gets `blink_aid/*-ica.fif`, `*_sources_raw.fif` and
`*_blink_aid.json`. Cached source traces preserve native sampling, first_samp,
sample count and date, with NaN gaps for excluded/short stretches and FIR guards.
`reports/sub-XX/ses-01/03_B_blink_aid_preparation.json` reports per-run preparation
or failures, never acceptance. A model/source/mask/settings mismatch stops stale
reuse. Source or fit-mask/settings changes require only affected cache rebuilds;
component or detector-threshold changes require local re-detection only. Use
`DAD_BLINK_AID_ONLY=1 DAD_REBUILD=1 DAD_SOURCE=<exact filename> sbatch --array=<participant> ...`
for a targeted rebuild. Regenerating a cache never edits reviewed manual work.

Locally append `--preliminary-ica` to the exact-source manual command above.
It loads the cache, opens component maps and native source traces, and prompts
first for blink components, then eye-movement components. `0` or `ICA000` selects
the first, `1` or `ICA001` the second; comma-separated numbers/names may be mixed.
Type exactly `NONE` for either unclear group; blank/invalid answers retry. The
same component can be selected for both groups. NONE skips only that detector,
never establishes absence of artifacts, and leaves manual review available.
It never fits locally.
Sources use arbitrary display units; temporary source/AID traces are normalized
by their maximum absolute amplitude so rare large blinks are not hidden by
browser clipping. Detection uses unmodified cached source amplitudes.

Detection low-passes the already 1-40 Hz sources at 10 Hz, separately across
finite stretches, trims its actual additional FIR support, and uses absolute
median/MAD scores with SciPy find_peaks. Starting settings in the existing
`artifact_annotation.blink_auto` block are five robust standard deviations,
0.25 s minimum peak spacing and +/-0.25 s candidate windows. These are review
starting points, not cohort acceptance or measured eyelid durations. Overlapping
proposals merge; broad peaks (>0.5 s half-prominence width) are flagged on the
score overview and printed for inspection/resizing, never silently discarded.
Flat/no-component/no-peak results remain unreviewed, not proof of no blinks.

The combined browser shows all original MISC (including MISC008 photodiode
context when present), EDF gaps/coverage, frontal MEG, selected AID traces and
editable `BAD_blink_auto` and `BAD_eye_movement_auto` proposals. `BLINK_score=1`
and `EYE_MOVEMENT_score=1` are their independent detector thresholds
(the browser score is capped at 2 solely for display);
the separate overview displays scores and the threshold line. Inspect the whole
recording, including unflagged and unscored periods. Candidate BAD labels exist
only on the disposable view until explicit SAVE; preparation produces no BAD
blink file usable by the batch pipeline. EDF remains discontinuous supplementary
context with no new synchronization or EOG workflow.

Corrected motion-onset lines and green `EPOCH_preview` annotations are shown
by default (approved 2026-10-26). Each preview spans the accepted Stage 3 window,
-0.5 to +2.667 seconds around its corrected onset. The existing read-only Stage
2/3 validators check ledger, review, timing bundle and eligibility application
hashes before display. The default review is
`<derivatives_root>/visual_timing/s3-v002/reviews/sub-XX_review.json`;
`--timing-review /path/to/sub-XX_review.json` or `epoching.timing_review_path`
can select an explicit accepted record. Verification may take time on shared
storage, especially when a participant uses the cohort donor evidence.

Green guides require an accepted epoch onset. Gaze-break trials additionally
have the separate invalid guides described below; no correction is invented for
unassigned trials. The terminal reports counts, including
incomplete windows, whose guides are clipped to this exact recording for display.
These are prospective windows before artifact/amplitude rejection, not a promise
that an epoch will survive Stage 5. Original STI101, native samples, first_samp,
recording boundaries and trial eligibility are unchanged.

Press `a` to open annotation controls. In Qt, choose **Select Visible** and toggle
`EPOCH_preview` or `EPOCH_invalid_preview`; in Matplotlib use its annotation-label
visibility checkboxes. Keep both labels unchanged. Guides can be resized/deleted in the temporary view,
but their edits are discarded when that browser closes; reopening restores the
accepted coordinates. They never enter the saved artifact file or any rejection
policy. `--no-timing-marks` disables both preview types and all timing-event lines
and skips their timing-evidence load. Existing filter-edge exposure stays
informational, and usual incomplete-epoch/reviewed-artifact rules still apply.

For thinner traces, append the already supported `--browser-backend matplotlib`
to the local annotation command. Its default trace width is 0.5 points on the
installed MNE versions. Qt remains the default. No custom trace-width mechanism
or dependency upgrade is needed. Green guides use the public annotation-color
option on newer MNE, or a scoped Matplotlib color palette on MNE 1.10/1.11; other
BAD annotations retain red. Changing these local display options does not
invalidate the cluster blink-aid caches or require a cluster rerun.

Corrected onset lines use width 2 by default, separately from the thin signal
traces. Append `--onset-linewidth 3` (or another finite positive value) to the
local command for stronger markers. Units are Matplotlib points or Qt logical
pixels. Installed MNE fixes its default event width, so a small display helper
styles the existing event artists and preserves the width when scrolling. It
does not change timing, traces, annotations, saved artifacts or cluster caches.

Onset lines are fully opaque black by default (approved 2026-10-27). Qt draws
them above annotation fills, so green epoch/red artifact shading cannot tint
them; signal traces keep their existing style. Use `--onset-color blue` or
`--onset-color "#0057B8"` to choose a named/RGB hex color. Transparent colors
are rejected; there is no opacity option. The event labels use the same color.
Color, width and opacity persist across scrolling in both backends. For example:

```bash
# Append to the usual exact-source local review command:
--onset-linewidth 3 --onset-color "#0057B8"
```

Restart the local helper to load the updated display defaults; no cluster rerun
or annotation reprocessing is needed.

Code 150 (`gaze_break`) is now displayed as an opaque orange vertical line above
traces, annotations and the other plot layers (approved 2026-10-28). It uses the
same `--onset-linewidth` value/default 2; there is no separate width control.
Its sample is the original native STI101 sample, without a visual-lag shift.
The canonical trigger table and exact ledger evidence tuple identify its owning
trial within the same acquisition/replay attempt; ambiguous/unmatched markers
are reported rather than assigned to a nearby trial. Contradictions with accepted
gaze-break eligibility stop the display for reconciliation.

Each matched already-excluded gaze-break trial gets a pale-orange
`EPOCH_invalid_preview` around its existing stored corrected onset, with the same
-0.5/+2.667 s intended window. It does not also get a green `EPOCH_preview`.
The display reads the excluded trial's hash-bound Stage 3 entry, validates its
native coordinates and shift arithmetic, and reports missing usable timing.
These are counterfactual intended windows, not usable epochs. Other excluded
trials keep their existing display behavior. The ordinary onset-color option
does not recolor orange gaze-break markers. Both preview types are independently
toggleable and discarded on browser close/SAVE. No BAD interval is created from
code 150, and no accepted eligibility/timing, artifact decisions or cluster
cache changes. Restart the local helper; no cluster rerun is required.

Every browser opening now offers `BAD_blink`, `BAD_eye_movement` and `BAD_segment`
without needing Add Description (approved 2026-10-27). Use `BAD_blink` for manual
blinks, keeping them distinct from `BAD_blink_auto`; use `BAD_segment` for unclear
or non-prototypical artifacts. Existing labels/intervals are retained and choices
are not duplicated. Qt starts on the manual `BAD_blink` choice; Matplotlib keeps
its usual sorted selector. BACK/REDETECT browser reopenings offer the same labels.

These are UI choices, with no placeholder/zero-duration annotations. If you never
use `BAD_segment`, that description is absent from the saved FIF, rather than
stored with an empty interval list. Reopening offers it again automatically.
SAVE writes only actual reviewed BAD intervals; unused descriptions change no
fit/epoch mask or lineage. Restart the helper to load the added choices; no
cluster rerun or annotation-file migration is needed.

After closing the browser, `BACK` reopens both selections on the same model.
Unchanged groups keep their results; changed nonempty groups ask separately
whether to regenerate. NONE disables only its group and removes only its
untouched session proposals. Declining regeneration keeps the previous scores/
proposal provenance separate from the newly displayed selection. `REDETECT`
asks `BLINK` or `EYE`, then adjusts only that group's robust threshold.
Regeneration replaces only that group's untouched proposals made in this session;
manual edits, deleted intervals and previous saved annotations remain protected.
Deleted/moved proposal spans are remembered as proposal vetoes, not artifact or
epoch exclusions. Reopening unchanged reviewed proposals does not reseed them.
MISC usefulness/meanings and observations are requested AFTER inspection when
SAVE is selected. SAVE marks retained auto intervals as reviewed artifacts,
saves their model/detector provenance in the existing receipt, and gives them
the usual final-ICA/epoch rejection behavior. DISCARDALL leaves saved files and
receipts unchanged; blank/typo actions retry. Use --overwrite to reopen saved
work. Full-source review and the independent 4B reconciliation gate still apply.

The same browser is also muscle QC (approved 2026-10-26). You can move, resize,
delete or add `BAD_muscle` intervals as well as ocular/other BAD labels. SAVE now
retains every BAD category in the existing exact-source manual annotation file.
When current provisional 4B evidence is available, its settings/source lineage
is checked and SAVE records the complete reviewed muscle set as a replacement,
including an explicit zero when all muscle intervals were deleted. Reopening
restores that set without resurrecting the automatic intervals. Keep the exact
`BAD_muscle` description for muscle segments; other BAD labels remain distinct.

The next explicit 4B rebuild applies this reviewed muscle set instead of adding
automatic detections back to it. Scores and automatic detection counts remain
diagnostic evidence; reviewed muscle intervals are not fed back into their own
score normalization. The report identifies `muscle_review_mode` and the applied
count. Revised ocular/severe masks can still change diagnostic scores, so inspect
the rebuilt scores/annotations/filter preview before the independent 4B gate.
Native/auxiliary samples and all other source-annotation preservation checks
remain in place. No Stage 5 decision or filter-edge rejection rule is introduced.

Older ocular-only receipts do not imply reviewed-zero muscle: automatic muscle
labels remain active until a new complete-source SAVE records muscle QC. Muscle
edits discarded by an older script were never stored and need to be redone.
Changing muscle detector settings/thresholds after muscle QC requires explicit
reconciliation of the affected source's muscle review before reuse; stale
replacements fail validation. Changed effective muscle masks also invalidate
the affected preliminary ICA cache. Rebuild that cache with the targeted
`DAD_BLINK_AID_ONLY=1 DAD_REBUILD=1 DAD_SOURCE=<exact filename>` command before
reusing it; unchanged masks keep compatible caches. DISCARDALL still changes
no saved files. Restart the local script to load this new SAVE behavior.

#### Separate eye-movement duration aid (approved 2026-10-28)

The eye detector reuses the existing `blink_aid/` sources without changing fit
identity, blink settings or the historical blink fingerprint formula. Its
`artifact_annotation.eye_movement_auto` settings are independent: 100 ms centered
median smoothing, absolute median/MAD threshold 4, minimum continuous excursion
100 ms, maximum merge gap 100 ms and 50 ms padding. These are provisional review
settings, not validated cohort thresholds. The median preserves step-like
changes and suppresses brief spikes; contiguous excursions retain duration
instead of becoming fixed blink windows. Ordinary activity should dominate the
finite samples used for robust normalization. The already cached 1-40 Hz filter
attenuates long gaze plateaus, so this cannot recover sustained eye position or
measure complete movement durations. Inspect unflagged periods manually.

Process each selected component and finite stretch separately, round the median
window to an odd native sample count, and trim half-window support at each
boundary. Apply minimum duration before gap merging, then pad/clip within that
same scored stretch. Merge overlapping proposals across selected eye components;
never bridge NaN gaps or recording boundaries. Keep source sampling/first_samp/
measurement date unchanged; record half-open absolute native sample endpoints,
effective smoothing support, diagnostics and unscored spans in the eye receipt.
No additional resampling, global alignment, classifier or dependency is needed.

The groups have distinct scores, settings, fingerprints, provisional intervals
and deletion/veto histories. Their automatic categories may overlap, including
when using the same component. Existing saved/edited intervals are protected;
only untouched session proposals for the chosen group are disposable. Manual
BAD_blink, BAD_eye_movement and BAD_segment remain available without dummy
intervals, and complete/reviewed-empty muscle replacements still survive SAVE
and 4B rebuilding. EDF and green/orange timing guides retain their existing
display-only contracts and onset styling controls.

Older blink-only receipts need no migration: missing eye selection means NONE
in memory, not corruption or a completed review. Loading/DISCARDALL writes
nothing. Only explicit complete-source SAVE stores new eye selection/detection
provenance under `eye_movement_auto`, separately from `blink_auto`. A selection
declined at BACK is recorded separately from the detector actually run; compact
earlier detector provenance remains when replacing its settings/selection.
Manual-only SAVE also retains previous detector provenance and deletion vetoes.
Reviewed-empty sets remain distinct from missing reviews.

This addition requires no cohort refit, SSS/muscle rerun or cache migration.
Current source/severe/muscle mask changes can independently require the existing
targeted cache rebuild; never bypass stale-cache checks to reopen supplementary
ICA. Finish any currently open older review normally, then restart the helper
to get the two prompts. Reviewed saves still require the ordinary affected 4B
rebuild and independent reconciliation/acceptance before any final 4C work.

References: [MNE ICA](https://mne.tools/stable/generated/mne.preprocessing.ICA.html),
[MNE blink detection](https://mne.tools/stable/generated/mne.preprocessing.find_eog_events.html),
[SciPy peak detection](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html),
[BLINKER](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2017.00012/full).
BLINKER supplies methodological context; its EEG validation is not a DAD MEG
performance claim. No classifier, dependency, epoch-exclusion margin or final
component-removal decision is introduced by this aid.

If an earlier SAVE wrote the annotation file but failed to write its tuning
receipt, or overlapping old review processes overwrote its receipt, append
`--overwrite --recover-annotations` to the exact-source command, optionally with
`--preliminary-ica`. This explicitly recovers the existing annotations as an
unreviewed starting point without reseeding EDF or ICA blink candidates. Saved
muscle edits replace provisional muscle intervals during recovery when the
saved file contains BAD_muscle; an orphan file without muscle labels cannot
establish reviewed zero. Current
EDF traces/provenance still load; recheck the full source/current guidance before
SAVE. Recovery itself is not acceptance. NumPy metadata now serialize through
the existing report conversion before tuning YAML is written. Since 2026-10-28,
SAVE reloads the latest tuning map and merges only the current source's receipt
under a short lock, then atomically replaces the YAML. Different source reviews
may overlap; avoid editing the same recording in two processes. After SAVE,
finish the MISC/observations prompts and wait for "Saved reviewed ocular
annotations (...)". SAVE alone precedes file writing. Already open processes
still run the older code: save/close them before starting any fixed-version
session. The small hidden .03_tuning_choices.yaml.lock is a write lock, not a
review record.

Update muscle thresholds in the existing tuning file as needed. Rebuild 4B only:

```bash
DAD_STAGES=artifact_annotation DAD_REBUILD=1 sbatch --array=4-39 meg/cluster/submit_preprocessing_review.sh
```

Inspect merged FIF/CSV annotations, full muscle scores, eye review receipts and
`final_filter_preview` in each run report. On 2026-10-28 the user superseded the
0.01 Hz final high-pass: retain low-frequency content for the roughly 3 s
onset-locked analysis, with residual drift QC still required. Keep spectrum-fit
50/100 Hz removal with its 10 s window, no final low-pass and no resampling.
P38 run 1 (`sub38_run01_051026.fif`, 137 s) receives spectrum-fit too: the user
clarified on 2026-10-28 that the short-run concern applied to the long high-pass.
Its original tuning exception remains a warning/guard for a future additional
high-pass and never bypasses current line removal. The fit copies remain separate.

Runtime/lineage migration is implemented and synthetically verified. Null
high-pass, independent spectrum-fit and margin-free 4B receipt handling are
supported. Omit --filter-edge-margin-sec for this policy. The 36 accepted 4A
receipts/reports remain unchanged: read-only legacy compatibility validates
their original configuration and source/output evidence, while new downstream
reports record the separate effective SSS identity. Do not rerun/reaccept 4A
for a final-temporal-only change; changed SQUID treatment follows the separate
[coordinated refresh](#coordinate-the-pending-4b-refresh-with-squid-review) above.

After deploying the updated workspace payload and resolving the SQUID branch
above, refresh only selected 4B evidence
with the following commands from the cluster DAD directory. The wrapper defaults
to guarded --resume, reusing compatible annotation FIF/CSV data and refreshing
previews/evidence. Before replacement, old temporal-policy JSON/embedded HTML
reports are archived under the participant report's `history/final_filter/`.
Use an explicit affected list; `4` below is an example after P04's manual saves
and SQUID decisions are reconciled and its SSS prerequisites are current, not a
cohort submission or an instruction to interrupt active review:

```bash
cd /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/DAD
DAD_AFFECTED_SUBJECTS=4
DAD_STAGES=artifact_annotation DAD_REBUILD=0 sbatch --array="${DAD_AFFECTED_SUBJECTS}%1" meg/cluster/submit_preprocessing_review.sh
```

If manual intervals/thresholds changed since the annotated FIF/CSV was built,
resume deliberately fails its existing signature check. After completing SAVE
and its prompts, use the explicit affected annotation rebuild instead:

```bash
DAD_STAGES=artifact_annotation DAD_REBUILD=1 sbatch --array="${DAD_AFFECTED_SUBJECTS}%1" meg/cluster/submit_preprocessing_review.sh
```

Neither command changes manual annotation FIFs or their exact-source receipts.
Keep the locked latest-map merge/atomic replacement and recovery workflow above.
Review every source and spectrum-fit boundary/short-stretch evidence, reconcile
4B, then record its receipt without a margin. No new receipt is recorded by this
implementation. The former 165 s exposure remains historical information.

`filter_edge_exposure` annotations never reject otherwise eligible epochs.
Stage 5 also normalizes the exact legacy `BAD_filter_edge` label in memory;
incomplete epochs and separately reviewed BAD artifacts retain their rules.
Spectrum-fit still needs boundary/short-stretch review and must not bridge BAD
periods. Formal 4B acceptance waits for independent inspection/reconciliation.

### 4C: protected final fit, model-specific decisions and application

For existing fitted models use [FINAL_ICA_COMPONENT_REVIEW_GUIDE.md](../docs/FINAL_ICA_COMPONENT_REVIEW_GUIDE.md). Array 11197 passed the joined 35-participant/180-acquisition audit: 2,548 full-fit property panels plus 5,096 default continuous/task panels, with unchanged models and matched sidecars/output stamps. `submit_final_ica_reports.sh` defaults to native-only mode, reusing validated fit properties. Native data are projected once per acquisition for all components. Set `DAD_ICA_REPORT_MODE=full` to also rebuild fit properties. No fitting, application or decisions occur; earlier HTML/changed figures are archived. P19's authorized numerical convergence reconciliation is separate from report completion.

P19 numerical repair is documented in
[P19_ICA_NUMERICAL_RECONCILIATION.md](../docs/P19_ICA_NUMERICAL_RECONCILIATION.md):
a near two-cycle, unsuccessful higher-cap retry, then a deterministic restart
of the saved seed-97 solution with 114 restart-phase iterations. It changes
the model weights, not fit configuration/masks/rank, and transfers no indices.
Require completed publication and the current joined readiness snapshot before
component review. Incompatible restart caches refuse silent default refitting.
No repeat cluster report array is needed for an already published local repair
on the shared derivatives. Deploy current code before future cluster operations.

The final joined audit (after all participant writes finish) is reproducible
with a NEW label; it does not fit/apply/accept data:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 meg/code/preprocessing_basic/audit_final_ica_readiness.py \
  --derivatives-root /mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives \
  --label NEW-UNUSED-LABEL --jobs 4
```

It joins four bounded participant processes before publication, verifying
current model hashes, all panel counts/sidecars and bound file stamps, plus
exact operative retained-trial identities and native timing bounds. It does
not reread all raw headers or large PNG/HTML hashes and grants no acceptance.


After accepted 4B only:

```bash
DAD_STAGES=ica_manual sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_review.sh
```

Fit copies omit every BAD category, including routine ocular intervals. Each good
stretch receives native 1-40 Hz filtering, full filter-edge guards and separate
200 Hz resampling before concatenation. The production decomposition uses MNE
method=fastica, n_components=None and random_state=97. Inspect actual dimension,
convergence, fit sample/duration accounting, maps, properties and source traces
with MISC/frontal context. No EOG scoring runs.

After independent component reconciliation, record exact-model exclusions:

```bash
python meg/code/preprocessing_basic/review_cleaning.py \
  --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
  --tuning-config /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/03_tuning_choices.yaml \
  --stage 4C --exclude --reviewer 'User and Codex' --reason 'Reconciled component evidence; no components selected'
```

Here --exclude with no numbers is an explicit reviewed empty list. Supply actual
reviewed indices when needed; no example component indices are cohort decisions.

Apply only the reviewed model, reusing accepted upstream derivatives:

```bash
sbatch --array=4-14,16-39 meg/cluster/submit_preprocessing_apply.sh
```

The apply wrapper selects ica_manual only. It cannot refit on changed inputs.
Inspect ICA-only before/after PNGs and final temporal-filter/edge/preservation
records, reconcile findings, and record applied cleaning acceptance separately.
The final outputs retain native rate/count/first_samp, source boundaries, channels,
STI101/MISC and annotations. The selected final processing is spectrum-fit
50/100 Hz removal only, with no additional high-pass or low-pass, including
P38 run 1. P17/P24 acquisition constraints and the future-high-pass warning remain.

For a final-filter-only refresh after the required 4B/4C reviews, use the affected
list with the existing apply wrapper:

```bash
sbatch --array="${DAD_AFFECTED_SUBJECTS}%1" meg/cluster/submit_preprocessing_apply.sh
```

It selects only ica_manual apply and reuses the exact compatible prefit model.
Final settings/operation-contract changes invalidate obsolete cleaned-output
reuse. If annotation masks changed, first refit/review only the affected participant
models through ordinary 4C; do not reuse changed fit inputs. Command construction
is tested with harmless runtime stubs and temporary synthetic CLI data, with no
scheduler submission or cohort processing during implementation validation.

For interruption recovery, use the existing runner with the same site/subject/
tuning configuration and `--stages ica_manual --processing-pass apply --resume`.
Changed upstream inputs require explicit affected-stage rebuild/refit/review;
changed exclusions alone reuse a compatible prefit model and rebuild final outputs.
Report-only rendering does not require scientific recomputation. Stage 5 remains
blocked until applied Stage 4 outputs and lineage have been scientifically accepted.

## Stage 5 epoching and subsequent script-06 export

`submit_epoch_export.sh` runs Stage 5 through `04_epoching.py`, then the downstream
`06_export_neural_for_drsa.py` export. The export is outside the preprocessing
endpoint and has its own prerequisites; the combined wrapper does not establish
scientific acceptance of either output.

Submit only after Stage 4 cleaning apply outputs are complete and their lineage is
accepted. The wrapper accepts P04-P39. Set `epoching.timing_review_path` to each
participant's accepted Stage 3 review record in the `epoching` block of its
subject selector before using this wrapper. For subject `sub-XX.yaml`, set
`timing_review_path` to
`/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives/visual_timing/s3-v002/reviews/sub-XX_review.json`.
For `s3-v002`, the review record accepts the exact adjacent `application.json`
hash and all required review scopes. The consumer discovers that application
automatically; `--eligibility-application` is an optional explicit path.
The cohort review guide's timing,
missing-data, and channel-consistency checks must be completed first.

**Current prerequisite:** Stage 3 acceptance is complete. The consumer verifies
the accepted timing JSON and the exact application hash/scopes; the bound
`application.json` is discovered beside the bundle unless the optional
`--eligibility-application` path is supplied. Stage 4 cleaning apply outputs and
lineage must be accepted before submission. The command below runs production;
it does not create acceptance or exclusions.

```bash
sbatch --array=4-39 meg/cluster/submit_epoch_export.sh
```

To resume selected supported participants, use an array selector such as
`--array=4`.
The default is guarded resume. Use `DAD_REBUILD=1` only when intentionally
replacing existing epoching/export outputs:

```bash
DAD_REBUILD=1 sbatch --array=24 meg/cluster/submit_epoch_export.sh
```

## Optional downstream and sanity-check jobs

These wrappers operate on their documented prerequisite derivatives. Check
each wrapper's header for all runtime controls and expected outputs.

| Script | Purpose |
| --- | --- |
| `submit_drsa_analysis.sh` | MATLAB dRSA per participant (script `07`) |
| `submit_drsa_group_average.sh` | Numeric dRSA group report |
| `submit_drsa_line_summaries.sh` | dRSA line-summary report |
| `submit_drsa_stats.sh` | Group dRSA matrix and line-summary statistics |
| `submit_dot_side_decoding.sh` | Per-subject dot-side time decoding |
| `submit_dot_side_group_report.sh` | Dot-side group report |
| `submit_dot_side_time_generalization.sh` | Per-subject train-time × test-time decoding |
| `submit_dot_side_time_generalization_group_report.sh` | Time-generalization group report |
| `submit_dot_erf_analysis.sh` | Group dot-condition ERF analysis |

Example submissions:

```bash
sbatch --array=4-39%4 meg/cluster/submit_drsa_analysis.sh
sbatch meg/cluster/submit_drsa_group_average.sh
sbatch meg/cluster/submit_drsa_line_summaries.sh
sbatch meg/cluster/submit_drsa_stats.sh
sbatch --array=1-144%4 meg/cluster/submit_dot_side_decoding.sh
sbatch meg/cluster/submit_dot_side_group_report.sh
sbatch --array=4-39%4 meg/cluster/submit_dot_side_time_generalization.sh
sbatch meg/cluster/submit_dot_side_time_generalization_group_report.sh
sbatch --array=1-4%4 meg/cluster/submit_dot_erf_analysis.sh
```

For dependent group reports, submit the report only after its subject job
succeeds. For example:

```bash
array_job_id="$(sbatch --parsable --array=4-39%4 \
  meg/cluster/submit_dot_side_time_generalization.sh | cut -d';' -f1)"
sbatch --dependency="afterok:${array_job_id}" \
  meg/cluster/submit_dot_side_time_generalization_group_report.sh
```

The MATLAB wrapper prefers the DIPSEN `matlab/R2020a` module for compatibility
and falls back to `matlab/R2025b` or another available MATLAB module. Set
`DAD_MATLAB_BIN=/absolute/path/to/matlab` before `sbatch` if you need a specific
executable.

## Former SGE guide

The previous runbook is preserved in [`README_SGE_legacy.md`](README_SGE_legacy.md)
for historical reference. Its `qsub` commands do not apply to these SLURM
wrappers.
