#!/usr/bin/env python3
"""
04_epoching.py

Purpose
Adapt the FLUX MNE-Python `Conditions.ipynb` step to the local DAD MEG
dataset with minimal changes, using robust trigger information from
`02_trigger_extraction_custom.py` when direct event reconstruction from raw
stim channels is unreliable.

Current implementation status
- This file is being implemented cell-by-cell to mirror FLUX structure.
- Current implementation includes:
  1) FLUX-style event preview from first cleaned run,
  2) all-trigger annotation injection in-memory from step-02 CSV,
  3) annotation round-trip validation via `events_from_annotations`,
  4) concatenation across runs, and
  5) condition-specific epoching for 31/41/51.
  6) propagation of upstream bad-channel metadata into the saved epochs header.

Usage example
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/preprocessing_basic/04_epoching.py

Planned inputs
- ICA-cleaned FIF files from step 03_C:
  .../derivatives/preprocessing/sub-01/ses-01/sss/*_ica.fif
  or
  .../derivatives/preprocessing/sub-01/ses-01/sss/*_ica_auto.fif
- Trigger events table from step 02:
  .../derivatives/preprocessing/sub-01/ses-01/events/02_trigger_events.csv

Planned outputs
- Epoched FIF:
  .../derivatives/analysis/sub-01/ses-01/meg/sub-01_ses-01_task-Attention_run-all_epo.fif
  (preserves upstream bad-channel metadata in `epochs.info['bads']`)
- HTML report:
  .../derivatives/reports/sub-01/ses-01/04_epoching_report.html
"""

# %%
# -------------------------------------------------------------------------
# Imports (FLUX preparation cell equivalent)
# -------------------------------------------------------------------------
import html
import os
import os.path as op
import csv
from bisect import bisect_left
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd

# %%
# -------------------------------------------------------------------------
# Subject/session, condition codes, and epoch-window settings
# -------------------------------------------------------------------------
subj_str = "01"
ses_str = "01"

# Choose cleaned input branch from step 03_C:
# - "manual" -> *_ica.fif
# - "auto"   -> *_ica_auto.fif
ica_variant = "manual"

# Condition mapping requested for this dataset.
event_id = {
    "always_visible": 31,
    "occluded_nondeviant": 41,
    "occluded_deviant": 51,
}

# Trial-level exclusion rule: reject condition trials containing gaze-break
# trigger(s) between trial onset and the next condition-onset.
gaze_break_code = 150

# Catch-trial exclusion for always-visible condition (31):
# reject only condition-31 trials containing trigger 116 in
# [trial_onset, next_condition_onset) (last trial: [trial_onset, run_end]).
catch_trial_reject_enabled = True
catch_code = 116
catch_applies_to_code = 31

# Epoch time window.
tmin = -1.0
tmax = 2.667

# I/O performance setting:
# keep run loading lazy during per-run event bookkeeping to reduce startup time
# and avoid reading full files before epoch extraction.
preload_runs_for_event_bookkeeping = False

# %%
# -------------------------------------------------------------------------
# Paths and output layout (FLUX set-paths cell adapted to local structure)
# -------------------------------------------------------------------------
project_root = "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG"

preproc_dir = op.join(
    project_root,
    "derivatives",
    "preprocessing",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
)
sss_dir = op.join(preproc_dir, "sss")
events_csv = op.join(preproc_dir, "events", "02_trigger_events.csv")

analysis_dir = op.join(
    project_root,
    "derivatives",
    "analysis",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
    "meg",
)
figure_output_dir = op.join(preproc_dir, "figures")
out_report = op.join(
    project_root,
    "derivatives",
    "reports",
    f"sub-{subj_str}",
    f"ses-{ses_str}",
    "04_epoching_report.html",
)
out_epochs = op.join(
    analysis_dir,
    f"sub-{subj_str}_ses-{ses_str}_task-Attention_run-all_epo.fif",
)
out_epoch_metadata_csv = out_epochs.replace("_epo.fif", "_epoch_metadata.csv")

os.makedirs(analysis_dir, exist_ok=True)
os.makedirs(figure_output_dir, exist_ok=True)
os.makedirs(op.dirname(out_report), exist_ok=True)

# Input cleaned-file suffix depends on selected ICA branch.
input_suffix = "_ica_auto.fif" if ica_variant == "auto" else "_ica.fif"
ica_files = sorted(
    op.join(sss_dir, fname)
    for fname in os.listdir(sss_dir)
    if fname.endswith(input_suffix)
)
if not ica_files:
    raise FileNotFoundError(
        f"No cleaned ICA files found in {sss_dir} with suffix {input_suffix}. "
        "Run 03_C_ICA.py or 03_C_ICA_auto.py first."
    )
if not op.exists(events_csv):
    raise FileNotFoundError(
        f"Required trigger events file not found: {events_csv}. "
        "Run 02_trigger_extraction_custom.py first."
    )

print(f"Selected ICA variant: {ica_variant}")
print(f"Input cleaned files ({len(ica_files)}):")
for pth in ica_files:
    print(f"  - {pth}")
print(f"Input trigger events CSV: {events_csv}")
print(f"Output epochs file: {out_epochs}")
print(f"Output report file: {out_report}")

# %%
# -------------------------------------------------------------------------
# Shared report/figure helpers (same style as other 03_* scripts)
# -------------------------------------------------------------------------
report = mne.Report(title=f"sub-{subj_str} | 04 Epoching")


def savefig_logged(fig, filename, **kwargs):
    """
    Save a matplotlib figure into the shared figure output folder and log path.
    """
    out_path = op.join(figure_output_dir, filename)
    fig.savefig(out_path, **kwargs)
    print(f"Saved figure to: {out_path}")
    return out_path


def log_to_console_and_report(text, *, title, section, tags=("log", "qc")):
    """
    Print text to console and store the same text in a report text section.
    """
    print(text)
    report.add_html(
        f"<pre>{html.escape(text)}</pre>",
        title=title,
        section=section,
        tags=tags,
    )


log_to_console_and_report(
    "FLUX Conditions adaptation: setup complete.\n"
    f"ICA branch={ica_variant}\n"
    f"event_id={event_id}\n"
    f"epoch_window=[{tmin}, {tmax}] sec\n"
    f"gaze_break_code={gaze_break_code}\n"
    f"catch_trial_reject_enabled={catch_trial_reject_enabled}\n"
    f"catch_code={catch_code}\n"
    f"catch_applies_to_code={catch_applies_to_code}\n"
    f"preload_runs_for_event_bookkeeping={preload_runs_for_event_bookkeeping}\n"
    f"n_input_files={len(ica_files)}\n"
    f"events_csv={events_csv}\n"
    f"out_epochs={out_epochs}",
    title="Setup summary",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# FLUX "Reading the events" adaptation: parse step-02 CSV with strict contract
# -------------------------------------------------------------------------
def cleaned_fname_to_source_meg_fname(cleaned_path):
    """
    Map cleaned run basename to source raw basename used by step-02 events CSV.

    Examples
    - SUB01_block01_ica_auto.fif -> SUB01_block01_meg.fif
    - SUB01_block01_ica.fif      -> SUB01_block01_meg.fif
    """
    base = op.basename(cleaned_path)
    if base.endswith("_ica_auto.fif"):
        return base.replace("_ica_auto.fif", "_meg.fif")
    if base.endswith("_ica.fif"):
        return base.replace("_ica.fif", "_meg.fif")
    raise ValueError(f"Unsupported cleaned filename pattern: {base}")


def csv_bool(value):
    """
    Convert CSV bool-like text from step-02 to Python bool.
    """
    return str(value).strip().lower() == "true"


required_event_csv_columns = {
    "trigger_code",
    "event_label",
    "sample_start_raw",
    "sample_end_raw",
    "time_start_raw_sec",
    "time_end_raw_sec",
    "source_fif",
    "is_transition_warning",
    "unresolved_transition",
}

event_rows_all = []
with open(events_csv, "r", encoding="utf-8") as fobj:
    reader = csv.DictReader(fobj)
    missing_cols = required_event_csv_columns.difference(set(reader.fieldnames or []))
    if missing_cols:
        raise RuntimeError(
            "02_trigger_events.csv missing required columns for epoching contract: "
            f"{sorted(missing_cols)}. Re-run 02_trigger_extraction_custom.py."
        )
    for row in reader:
        event_rows_all.append(row)

# Keep resolved events only (all trigger codes). Resolved transition-warning
# rows are retained because step 02 has already recovered the expected trigger
# and discarded the unexpected transition partner.
event_rows_resolved = [
    row
    for row in event_rows_all
    if not csv_bool(row["unresolved_transition"])
]
if not event_rows_resolved:
    raise RuntimeError("No resolved events found in 02_trigger_events.csv.")

rows_resolved_by_source = defaultdict(list)
for row in event_rows_resolved:
    rows_resolved_by_source[row["source_fif"]].append(row)
for source_fif_name in rows_resolved_by_source:
    rows_resolved_by_source[source_fif_name].sort(key=lambda row: int(row["sample_start_raw"]))

resolved_code_counts = Counter(int(row["trigger_code"]) for row in event_rows_resolved)
n_transition_warning_resolved = sum(csv_bool(row["is_transition_warning"]) for row in event_rows_resolved)
log_to_console_and_report(
    "Loaded step-02 event table with raw-frame fields.\n"
    f"n_rows_total={len(event_rows_all)}\n"
    f"n_rows_resolved={len(event_rows_resolved)}\n"
    f"n_transition_warning_resolved={n_transition_warning_resolved}\n"
    f"n_source_files={len(rows_resolved_by_source)}\n"
    f"resolved_code_counts={dict(sorted(resolved_code_counts.items()))}",
    title="Event CSV contract check",
    section="Events extraction",
)


def inject_annotations_and_validate(raw_obj, rows_source, source_fif_name):
    """
    Inject all resolved trigger annotations from CSV and validate round-trip.

    Validation target
    Compare `(sample_start_raw, event_label)` pairs from CSV against
    `events_from_annotations(..., event_id='auto')` output on the same raw.
    """
    sfreq = float(raw_obj.info["sfreq"])
    sample_min = int(raw_obj.first_samp)
    sample_max = int(raw_obj.first_samp + raw_obj.n_times - 1)

    onset = []
    duration = []
    description = []
    csv_pairs = []
    out_of_range_rows = 0
    annotation_orig_time = raw_obj.annotations.orig_time

    for row in rows_source:
        sample_start_raw = int(row["sample_start_raw"])
        label = str(row["event_label"])
        if sample_start_raw < sample_min or sample_start_raw > sample_max:
            out_of_range_rows += 1
            continue
        onset_rel = (sample_start_raw - raw_obj.first_samp) / sfreq
        # Requested epoching frame is run-relative (sample - first_samp). When
        # existing annotations use an absolute `orig_time`, shift by first_time
        # so the new annotations live in the same absolute frame.
        onset_store = onset_rel + raw_obj.first_time if annotation_orig_time is not None else onset_rel
        onset.append(onset_store)
        duration.append(0.0)
        description.append(label)
        csv_pairs.append((sample_start_raw, label))

    if out_of_range_rows:
        raise RuntimeError(
            f"Found {out_of_range_rows} out-of-range annotation rows for {source_fif_name}. "
            "Raw-frame sample contract is violated."
        )

    annotations = mne.Annotations(
        onset=onset,
        duration=duration,
        description=description,
        orig_time=annotation_orig_time,
    )
    raw_obj.set_annotations(raw_obj.annotations + annotations)

    events_auto, event_id_auto = mne.events_from_annotations(raw_obj, event_id="auto")
    id_to_label = {value: key for key, value in event_id_auto.items()}
    csv_labels = set(description)
    auto_pairs = [
        (int(sample), id_to_label[int(code)])
        for sample, _prev, code in events_auto
        if id_to_label[int(code)] in csv_labels
    ]

    csv_counter = Counter(csv_pairs)
    auto_counter = Counter(auto_pairs)
    if csv_counter != auto_counter:
        missing = list((csv_counter - auto_counter).items())[:5]
        extra = list((auto_counter - csv_counter).items())[:5]
        raise RuntimeError(
            "events_from_annotations round-trip mismatch for "
            f"{source_fif_name}. Missing={missing}; Extra={extra}"
        )

    return {
        "n_csv_pairs": len(csv_pairs),
        "n_auto_pairs": len(auto_pairs),
        "n_unique_labels": len(csv_labels),
    }


# FLUX equivalent loads run-01 for illustration; keep same intent locally.
raw_preview = mne.io.read_raw_fif(ica_files[0], preload=False, verbose=True)
source_fif_preview = cleaned_fname_to_source_meg_fname(ica_files[0])
rows_preview_all = rows_resolved_by_source.get(source_fif_preview, [])
if not rows_preview_all:
    raise RuntimeError(
        "No resolved events found for first cleaned run source file: "
        f"{source_fif_preview}"
    )

# Inject all trigger annotations for the preview file and validate against CSV.
preview_validation = inject_annotations_and_validate(
    raw_preview,
    rows_preview_all,
    source_fif_preview,
)

# Build condition events for epoching preview from raw-frame samples.
target_codes = set(event_id.values())
event_rows_preview = [
    row for row in rows_preview_all if int(row["trigger_code"]) in target_codes
]
preview_events = np.array(
    [
        [int(row["sample_start_raw"]), 0, int(row["trigger_code"])]
        for row in event_rows_preview
    ],
    dtype=int,
)
if preview_events.size == 0:
    raise RuntimeError(
        "No preview condition events found for first cleaned run. "
        f"source_fif={source_fif_preview}, event_codes={sorted(target_codes)}"
    )
preview_events = preview_events[np.argsort(preview_events[:, 0])]

preview_counts = {
    name: int(np.sum(preview_events[:, 2] == code))
    for name, code in event_id.items()
}
log_to_console_and_report(
    "Loaded first cleaned run for FLUX-style event illustration.\n"
    f"raw_preview={ica_files[0]}\n"
    f"source_fif_preview={source_fif_preview}\n"
    f"n_preview_events={len(preview_events)}\n"
    f"preview_validation={preview_validation}\n"
    f"event_counts={preview_counts}\n"
    f"first_preview_event={preview_events[0].tolist()}",
    title="Preview run and events",
    section="Events extraction",
)

# %%
# -------------------------------------------------------------------------
# FLUX plotting equivalent 1: preview raw segment from first cleaned run
# -------------------------------------------------------------------------
# FLUX notebook uses `raw.plot(start=50)` for visual trigger/context checks.
sti_ch_names = [ch for ch in raw_preview.ch_names if ch.upper().startswith("STI")]
if not sti_ch_names:
    raise RuntimeError(
        "No STI channels found in preview run. Cannot generate STI preview plot."
    )
sti_ch_names = sti_ch_names[:6]

fig_raw_preview = raw_preview.plot(start=50, picks=sti_ch_names, show=True)
report.add_figure(
    fig=fig_raw_preview,
    title=f"Raw preview ({op.basename(ica_files[0])})",
    caption=(
        "FLUX-equivalent first-look plot from the first cleaned run "
        f"(first STI channels: {', '.join(sti_ch_names)})."
    ),
    section="Events extraction",
    tags=("raw", "qc"),
)
savefig_logged(
    fig_raw_preview,
    "04_epoching_preview_raw.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig_raw_preview)

# %%
# -------------------------------------------------------------------------
# FLUX plotting equivalent 2: first N trigger events stem plot
# -------------------------------------------------------------------------
# Keep FLUX convention of plotting the first ~80 event rows as sample/code stem.
first_n_events = 80
n_plot = min(first_n_events, len(preview_events))

fig_events_stem, ax = plt.subplots(figsize=(10, 4))
markerline, stemlines, baseline = ax.stem(
    preview_events[:n_plot, 0],
    preview_events[:n_plot, 2],
)
plt.setp(markerline, markersize=4)
plt.setp(stemlines, linewidth=1.0)
plt.setp(baseline, linewidth=0.8)
ax.set_xlabel("samples")
ax.set_ylabel("Trigger value")
ax.set_title(f"First {n_plot} preview events (codes 31/41/51)")

report.add_figure(
    fig=fig_events_stem,
    title="Event stem plot (first N events)",
    caption=(
        "FLUX-equivalent event-array visualization. "
        "X-axis is sample index, Y-axis is trigger code."
    ),
    section="Events extraction",
    tags=("events", "qc"),
)
savefig_logged(
    fig_events_stem,
    "04_epoching_preview_events_stem.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig_events_stem)

log_to_console_and_report(
    "Generated FLUX-equivalent preview plots.\n"
    f"raw_plot_start_sec=50\n"
    f"stem_plot_n_events={n_plot}",
    title="Preview plots",
    section="Events extraction",
)

# %%
# -------------------------------------------------------------------------
# FLUX conditions block adaptation: per-run event arrays + concatenation
# -------------------------------------------------------------------------
def identity_row_for_trial(row, idx, rows_run_cond, rows_run_all, sample_max, source_fif_name):
    """
    Find the one trial-sequence identity trigger belonging to one condition row.

    The interval boundary uses all resolved condition rows, including resolved
    transition warnings, so a recovered expected condition trigger owns only the
    identity inside its own trial interval.
    """
    start_sample = int(row["sample_start_raw"])
    end_sample_exclusive = (
        int(rows_run_cond[idx + 1]["sample_start_raw"])
        if idx + 1 < len(rows_run_cond)
        else sample_max + 1
    )
    identity_rows = [
        candidate
        for candidate in rows_run_all
        if start_sample < int(candidate["sample_start_raw"]) < end_sample_exclusive
        and str(candidate["event_label"]) == "trial_sequence_identity"
        and 1 <= int(candidate["trigger_code"]) <= 20
        and not csv_bool(candidate["is_transition_warning"])
    ]
    if len(identity_rows) != 1:
        raise RuntimeError(
            "Ambiguous trial_sequence_identity mapping for "
            f"{source_fif_name}: condition_sample={start_sample}, "
            f"condition_code={row['trigger_code']}, "
            f"interval_end={end_sample_exclusive}, "
            f"identity_candidates={[candidate['sample_start_raw'] for candidate in identity_rows]}"
        )
    return identity_rows[0]


def build_epoch_metadata_row(row, identity_row, source_fif_name):
    """
    Build one metadata record aligned with one condition event row.
    """
    return {
        "trial_sequence_identity": int(identity_row["trigger_code"]),
        "condition_name": str(row["event_label"]).replace("trial_condition_onset_", ""),
        "condition_code": int(row["trigger_code"]),
        "source_fif": source_fif_name,
        "sample_start_raw": int(row["sample_start_raw"]),
        "sample_end_raw": int(row["sample_end_raw"]),
        "time_start_raw_sec": float(row["time_start_raw_sec"]),
        "time_end_raw_sec": float(row["time_end_raw_sec"]),
        "is_transition_warning": csv_bool(row["is_transition_warning"]),
        "transition_partner_code": str(row.get("transition_partner_code", "")),
        "expected_trigger_code": str(row.get("expected_trigger_code", "")),
        "transition_decision_source": str(row.get("transition_decision_source", "")),
        "transition_resolution_confidence": str(row.get("transition_resolution_confidence", "")),
        "identity_sample_start_raw": int(identity_row["sample_start_raw"]),
    }


raw_list = []
events_list = []
epoch_metadata_rows = []
per_run_summary = []
global_counts_before_gaze = Counter()
global_counts_after_gaze = Counter()
global_counts_rejected_gaze = Counter()
global_trials_before_gaze = 0
global_trials_after_gaze = 0
global_counts_before_catch = Counter()
global_counts_after_catch = Counter()
global_counts_rejected_catch = Counter()
global_counts_would_reject_catch = Counter()
global_trials_before_catch = 0
global_trials_after_catch = 0

for clean_path in ica_files:
    raw_run = mne.io.read_raw_fif(
        clean_path,
        preload=preload_runs_for_event_bookkeeping,
        verbose=True,
    )
    source_fif_run = cleaned_fname_to_source_meg_fname(clean_path)
    rows_run_all = rows_resolved_by_source.get(source_fif_run, [])
    if not rows_run_all:
        raise RuntimeError(f"No resolved events found for source file {source_fif_run}")

    run_validation = inject_annotations_and_validate(raw_run, rows_run_all, source_fif_run)

    rows_run_cond = [row for row in rows_run_all if int(row["trigger_code"]) in target_codes]
    if not rows_run_cond:
        raise RuntimeError(
            f"No condition events (31/41/51) found for {source_fif_run}. "
            "Cannot continue epoching."
        )

    sample_min = int(raw_run.first_samp)
    sample_max = int(raw_run.first_samp + raw_run.n_times - 1)
    cond_samples = np.array([int(row["sample_start_raw"]) for row in rows_run_cond], dtype=int)
    if cond_samples.min() < sample_min or cond_samples.max() > sample_max:
        raise RuntimeError(
            f"Event sample out of range in {source_fif_run}: "
            f"valid=[{sample_min}, {sample_max}], "
            f"observed=[{cond_samples.min()}, {cond_samples.max()}]"
        )

    # Trial-level gaze-break exclusion:
    # reject trial i if a gaze-break occurs in [trial_i_onset, trial_{i+1}_onset).
    # For the last trial in a run, use [trial_last_onset, run_end].
    gaze_break_samples = sorted(
        int(row["sample_start_raw"])
        for row in rows_run_all
        if int(row["trigger_code"]) == gaze_break_code
    )
    run_events_after_gaze_rows = []
    run_events_rows = []
    run_metadata_rows = []
    rejected_by_gaze_rows = []
    rejected_by_catch_rows = []
    catch_would_reject_rows = []
    catch_samples = sorted(
        int(row["sample_start_raw"])
        for row in rows_run_all
        if int(row["trigger_code"]) == catch_code
    )
    for idx, row in enumerate(rows_run_cond):
        start_sample = int(row["sample_start_raw"])
        trial_code = int(row["trigger_code"])
        end_sample_exclusive = (
            int(rows_run_cond[idx + 1]["sample_start_raw"])
            if idx + 1 < len(rows_run_cond)
            else sample_max + 1
        )
        gb_idx = bisect_left(gaze_break_samples, start_sample)
        has_gaze_break = (
            gb_idx < len(gaze_break_samples)
            and gaze_break_samples[gb_idx] < end_sample_exclusive
        )
        if has_gaze_break:
            rejected_by_gaze_rows.append(row)
            continue

        run_events_after_gaze_rows.append(row)

        has_catch = False
        if trial_code == catch_applies_to_code:
            catch_idx = bisect_left(catch_samples, start_sample)
            has_catch = (
                catch_idx < len(catch_samples)
                and catch_samples[catch_idx] < end_sample_exclusive
            )
        if has_catch:
            catch_would_reject_rows.append(row)
            if catch_trial_reject_enabled:
                rejected_by_catch_rows.append(row)
                continue

        identity_row = identity_row_for_trial(
            row,
            idx,
            rows_run_cond,
            rows_run_all,
            sample_max,
            source_fif_run,
        )
        run_events_rows.append(row)
        run_metadata_rows.append(build_epoch_metadata_row(row, identity_row, source_fif_run))

    if not run_events_rows:
        raise RuntimeError(
            "All condition events were rejected by gaze/catch rules in "
            f"{source_fif_run}. "
            "Cannot continue epoching."
        )

    run_events = np.array(
        [
            [int(row["sample_start_raw"]), 0, int(row["trigger_code"])]
            for row in run_events_rows
        ],
        dtype=int,
    )
    run_events = run_events[np.argsort(run_events[:, 0])]

    raw_list.append(raw_run)
    events_list.append(run_events)
    epoch_metadata_rows.extend(run_metadata_rows)

    cond_counts_before = Counter(int(row["trigger_code"]) for row in rows_run_cond)
    cond_counts_after_gaze = Counter(int(row["trigger_code"]) for row in run_events_after_gaze_rows)
    cond_counts_rejected_gaze = Counter(int(row["trigger_code"]) for row in rejected_by_gaze_rows)
    cond_counts_after = Counter(int(row["trigger_code"]) for row in run_events_rows)
    cond_counts_rejected_catch = Counter(int(row["trigger_code"]) for row in rejected_by_catch_rows)
    cond_counts_would_reject_catch = Counter(int(row["trigger_code"]) for row in catch_would_reject_rows)
    global_counts_before_gaze.update(cond_counts_before)
    global_counts_after_gaze.update(cond_counts_after_gaze)
    global_counts_rejected_gaze.update(cond_counts_rejected_gaze)
    global_trials_before_gaze += len(rows_run_cond)
    global_trials_after_gaze += len(run_events_after_gaze_rows)
    global_counts_before_catch.update(cond_counts_after_gaze)
    global_counts_after_catch.update(cond_counts_after)
    global_counts_rejected_catch.update(cond_counts_rejected_catch)
    global_counts_would_reject_catch.update(cond_counts_would_reject_catch)
    global_trials_before_catch += len(run_events_after_gaze_rows)
    global_trials_after_catch += len(run_events_rows)

    per_run_summary.append(
        {
            "clean_file": op.basename(clean_path),
            "source_file": source_fif_run,
            "n_all_resolved_events": len(rows_run_all),
            "n_condition_events_before_gaze": len(rows_run_cond),
            "n_condition_events_after_gaze": len(run_events_after_gaze_rows),
            "n_condition_events_rejected_gaze": len(rejected_by_gaze_rows),
            "n_condition_events_before_catch": len(run_events_after_gaze_rows),
            "n_condition_events_after_catch": len(run_events_rows),
            "n_condition_events_rejected_catch": len(rejected_by_catch_rows),
            "n_condition_events_would_reject_catch": len(catch_would_reject_rows),
            "n_gaze_break_events": len(gaze_break_samples),
            "n_catch_events": len(catch_samples),
            "condition_counts_before_gaze": dict(cond_counts_before),
            "condition_counts_after_gaze": dict(cond_counts_after_gaze),
            "condition_counts_rejected_gaze": dict(cond_counts_rejected_gaze),
            "condition_counts_before_catch": dict(cond_counts_after_gaze),
            "condition_counts_after_catch": dict(cond_counts_after),
            "condition_counts_rejected_catch": dict(cond_counts_rejected_catch),
            "condition_counts_would_reject_catch": dict(cond_counts_would_reject_catch),
            "annotation_validation": run_validation,
        }
    )

log_to_console_and_report(
    "Per-run event extraction summary:\n"
    + "\n".join(
        f"  - {item['clean_file']} ({item['source_file']}): "
        f"resolved={item['n_all_resolved_events']}, "
        f"condition_before_gaze={item['n_condition_events_before_gaze']}, "
        f"condition_after_gaze={item['n_condition_events_after_gaze']}, "
        f"rejected_gaze={item['n_condition_events_rejected_gaze']}, "
        f"n_gaze_break_events={item['n_gaze_break_events']}, "
        f"counts_before_gaze={item['condition_counts_before_gaze']}, "
        f"counts_after_gaze={item['condition_counts_after_gaze']}, "
        f"counts_rejected_gaze={item['condition_counts_rejected_gaze']}, "
        f"condition_before_catch={item['n_condition_events_before_catch']}, "
        f"condition_after_catch={item['n_condition_events_after_catch']}, "
        f"rejected_catch={item['n_condition_events_rejected_catch']}, "
        f"would_reject_catch={item['n_condition_events_would_reject_catch']}, "
        f"n_catch_events={item['n_catch_events']}, "
        f"counts_before_catch={item['condition_counts_before_catch']}, "
        f"counts_after_catch={item['condition_counts_after_catch']}, "
        f"counts_rejected_catch={item['condition_counts_rejected_catch']}, "
        f"counts_would_reject_catch={item['condition_counts_would_reject_catch']}, "
        f"annot_check={item['annotation_validation']}"
        for item in per_run_summary
    ),
    title="Per-run event arrays",
    section="Events extraction",
)

log_to_console_and_report(
    "Gaze-break trial filter summary.\n"
    f"gaze_break_code={gaze_break_code}\n"
    f"interval_rule=[trial_onset, next_trial_onset)\n"
    f"n_trials_before_gaze={global_trials_before_gaze}\n"
    f"n_trials_after_gaze={global_trials_after_gaze}\n"
    f"n_trials_rejected_gaze={global_trials_before_gaze - global_trials_after_gaze}\n"
    f"counts_before_gaze={dict(sorted(global_counts_before_gaze.items()))}\n"
    f"counts_after_gaze={dict(sorted(global_counts_after_gaze.items()))}\n"
    f"counts_rejected_gaze={dict(sorted(global_counts_rejected_gaze.items()))}",
    title="Gaze-break filter",
    section="Events extraction",
)

log_to_console_and_report(
    "Catch trial filter summary.\n"
    f"catch_trial_reject_enabled={catch_trial_reject_enabled}\n"
    f"catch_code={catch_code}\n"
    f"catch_applies_to_code={catch_applies_to_code}\n"
    f"interval_rule=[trial_onset, next_trial_onset)\n"
    f"n_trials_before_catch={global_trials_before_catch}\n"
    f"n_trials_would_reject_catch={sum(global_counts_would_reject_catch.values())}\n"
    f"n_trials_rejected_catch={sum(global_counts_rejected_catch.values())}\n"
    f"n_trials_after_catch={global_trials_after_catch}\n"
    f"counts_before_catch={dict(sorted(global_counts_before_catch.items()))}\n"
    f"counts_would_reject_catch={dict(sorted(global_counts_would_reject_catch.items()))}\n"
    f"counts_rejected_catch={dict(sorted(global_counts_rejected_catch.items()))}\n"
    f"counts_after_catch={dict(sorted(global_counts_after_catch.items()))}\n"
    f"filter_mode={'ON (trials removed)' if catch_trial_reject_enabled else 'OFF (report-only, no trial removed)'}",
    title="Catch filter",
    section="Events extraction",
)

raw_all, events = mne.concatenate_raws(raw_list, events_list=events_list, on_mismatch="warn")
del raw_list
# Keep one explicit copy of propagated bad channels so epochs saving does not
# depend on implicit MNE internals.
raw_all_bads = list(dict.fromkeys(str(ch) for ch in raw_all.info.get("bads", [])))

if len(epoch_metadata_rows) != len(events):
    raise RuntimeError(
        "Epoch metadata row count does not match concatenated event count: "
        f"metadata={len(epoch_metadata_rows)}, events={len(events)}"
    )
event_codes_after_concat = [int(code) for _sample, _prev, code in events]
metadata_codes = [int(row["condition_code"]) for row in epoch_metadata_rows]
if event_codes_after_concat != metadata_codes:
    raise RuntimeError("Epoch metadata condition-code order does not match events order.")

log_to_console_and_report(
    "Concatenated runs for epoching.\n"
    f"n_channels={raw_all.info['nchan']}\n"
    f"n_times={raw_all.n_times}\n"
    f"sfreq={raw_all.info['sfreq']}\n"
    f"n_events={len(events)}\n"
    f"propagated_bads={raw_all_bads}",
    title="Concatenation summary",
    section="Events extraction",
)

# %%
# -------------------------------------------------------------------------
# FLUX rejection and epoching block adapted to requested event_ids/window
# -------------------------------------------------------------------------
# Optional blink-window-center trial rejection (ON/OFF).
# Rule:
# - Compute blink center sample from annotation window center.
# - Reject trial if center is in [trial_onset, next_trial_onset).
# - For last trial in a run, use [trial_onset, run_end].
blink_trial_reject_enabled = False
blink_annotation_keyword = "blink"

events_for_epoching = events.copy()
event_sort_order = np.argsort(events_for_epoching[:, 0])
events_for_epoching = events_for_epoching[event_sort_order]
metadata_for_epoching = [epoch_metadata_rows[idx] for idx in event_sort_order]

# Identify run boundaries in concatenated raw (inserted by MNE as boundary annotations).
boundary_samples = sorted(
    {
        int(raw_all.time_as_index(onset, use_rounding=True)[0])
        for onset, desc in zip(raw_all.annotations.onset, raw_all.annotations.description)
        if "boundary" in str(desc).lower()
    }
)
raw_end_exclusive = int(raw_all.first_samp + raw_all.n_times)

# Extract blink-window centers from annotations (center = onset + duration/2).
blink_center_samples = sorted(
    int(raw_all.time_as_index(onset + duration / 2.0, use_rounding=True)[0])
    for onset, duration, desc in zip(
        raw_all.annotations.onset,
        raw_all.annotations.duration,
        raw_all.annotations.description,
    )
    if blink_annotation_keyword.lower() in str(desc).lower()
)

blink_keep_mask = np.ones(len(events_for_epoching), dtype=bool)
blink_rejected_rows = []
for idx_trial in range(len(events_for_epoching)):
    start_sample = int(events_for_epoching[idx_trial, 0])
    code = int(events_for_epoching[idx_trial, 2])
    next_trial_start = (
        int(events_for_epoching[idx_trial + 1, 0])
        if idx_trial + 1 < len(events_for_epoching)
        else raw_end_exclusive
    )

    # Limit trial window to current run end, so the last trial in each run
    # follows [trial_onset, run_end] instead of crossing into next run.
    boundary_idx = bisect_left(boundary_samples, start_sample + 1)
    run_end_exclusive = (
        boundary_samples[boundary_idx]
        if boundary_idx < len(boundary_samples)
        else raw_end_exclusive
    )
    end_exclusive = min(next_trial_start, run_end_exclusive)

    center_idx = bisect_left(blink_center_samples, start_sample)
    has_blink_center = (
        center_idx < len(blink_center_samples)
        and blink_center_samples[center_idx] < end_exclusive
    )
    if has_blink_center:
        blink_keep_mask[idx_trial] = False
        blink_rejected_rows.append((start_sample, code))

blink_counts_before = Counter(int(code) for code in events_for_epoching[:, 2])
blink_counts_rejected = Counter(code for _sample, code in blink_rejected_rows)
if blink_trial_reject_enabled:
    events_for_epoching = events_for_epoching[blink_keep_mask]
    metadata_for_epoching = [
        row for row, keep in zip(metadata_for_epoching, blink_keep_mask) if keep
    ]
    if len(events_for_epoching) == 0:
        raise RuntimeError(
            "Blink-center trial rejection removed all condition events. "
            "Disable blink_trial_reject_enabled or inspect blink annotations."
        )
blink_counts_after = Counter(int(code) for code in events_for_epoching[:, 2])

log_to_console_and_report(
    "Blink-center trial rejection summary.\n"
    f"blink_trial_reject_enabled={blink_trial_reject_enabled}\n"
    f"blink_annotation_keyword={blink_annotation_keyword}\n"
    "interval_rule=[trial_onset, next_trial_onset), "
    "last trial in run [trial_onset, run_end]\n"
    f"n_boundary_markers={len(boundary_samples)}\n"
    f"n_blink_windows_detected={len(blink_center_samples)}\n"
    f"n_trials_before_blink_filter={int(len(events))}\n"
    f"n_trials_rejected_by_blink_center={int(np.sum(~blink_keep_mask))}\n"
    f"n_trials_after_blink_filter={int(len(events_for_epoching))}\n"
    f"counts_before_blink_filter={dict(sorted(blink_counts_before.items()))}\n"
    f"counts_rejected_by_blink_center={dict(sorted(blink_counts_rejected.items()))}\n"
    f"counts_after_blink_filter={dict(sorted(blink_counts_after.items()))}\n"
    f"filter_mode={'ON (trials removed)' if blink_trial_reject_enabled else 'OFF (report-only, no trial removed)'}",
    title="Blink-center filter",
    section="Epoching",
)

reject = dict(
    grad=5000e-13,  # T/m
    mag=4e-12,      # T
)

epochs = mne.Epochs(
    raw_all,
    events_for_epoching,
    event_id=event_id,
    tmin=tmin,
    tmax=tmax,
    baseline=None,
    proj=False,
    picks="all",
    detrend=1,
    reject=reject,
    reject_by_annotation=True,
    metadata=pd.DataFrame(metadata_for_epoching),
    preload=True,
    verbose=True,
)
# Preserve upstream bad channels explicitly in epochs metadata/header.
epochs.info["bads"] = list(raw_all_bads)

epoch_counts = {cond_name: len(epochs[cond_name]) for cond_name in event_id}
log_to_console_and_report(
    "Epoching completed.\n"
    f"tmin={tmin}, tmax={tmax}\n"
    f"event_id={event_id}\n"
    f"n_epochs_total={len(epochs)}\n"
    f"n_epochs_by_condition={epoch_counts}\n"
    f"metadata_columns={list(epochs.metadata.columns) if epochs.metadata is not None else []}\n"
    f"n_marked_bad_channels={len(epochs.info.get('bads', []))}\n"
    f"marked_bad_channels={epochs.info.get('bads', [])}",
    title="Epoching summary",
    section="Epoching",
)

counts_before_selected = {
    code: int(global_counts_before_gaze.get(code, 0))
    for code in sorted(event_id.values())
}
counts_after_gaze_selected = {
    code: int(global_counts_after_gaze.get(code, 0))
    for code in sorted(event_id.values())
}
counts_after_catch_selected = {
    code: int(global_counts_after_catch.get(code, 0))
    for code in sorted(event_id.values())
}
final_counts_by_code = {
    code: int(len(epochs[cond_name]))
    for cond_name, code in sorted(event_id.items(), key=lambda item: item[1])
}
log_to_console_and_report(
    "Validation run completed successfully:\n\n"
    f"Trials before gaze filter: {global_trials_before_gaze}\n"
    f"Trials after gaze filter: {global_trials_after_gaze}\n"
    f"Rejected by gaze rule: {global_trials_before_gaze - global_trials_after_gaze}\n"
    f"Counts before gaze filter: {counts_before_selected}\n"
    f"Counts after gaze filter: {counts_after_gaze_selected}\n"
    f"Catch filter enabled: {catch_trial_reject_enabled}\n"
    f"Trials before catch filter: {global_trials_before_catch}\n"
    f"Trials would reject by catch: {sum(global_counts_would_reject_catch.values())}\n"
    f"Trials rejected by catch: {sum(global_counts_rejected_catch.values())}\n"
    f"Trials after catch filter: {global_trials_after_catch}\n"
    f"Counts after catch filter: {counts_after_catch_selected}\n"
    "Final epochs after amplitude rejection: "
    f"{len(epochs)} ({', '.join(f'{code}:{final_counts_by_code[code]}' for code in sorted(final_counts_by_code))})",
    title="Validation summary",
    section="Epoching",
)

# %%
# -------------------------------------------------------------------------
# FLUX "Plotting the trials" adaptation
# -------------------------------------------------------------------------
# FLUX cell equivalents:
# - epochs['cue_left'].plot(n_epochs=10, picks=['grad'])
# - epochs['cue_left'].plot(n_epochs=1, picks=['stim'])
# - epochs['cue_left'].filter(0.0,30).crop(-1,2.667).plot_image(picks=['MEG2343'], ...)
plot_trials_condition = "always_visible"
plot_trials_n_epochs_grad = 10
plot_trials_n_epochs_stim = 1
plot_trials_image_channel = "MEG2343"

if plot_trials_condition not in epochs.event_id:
    raise RuntimeError(
        f"Plot condition '{plot_trials_condition}' not in epochs.event_id={epochs.event_id}."
    )

# 1) Trial browser plot for gradiometers (non-interactive, saved + report).
fig_trials_grad = epochs[plot_trials_condition].plot(
    n_epochs=plot_trials_n_epochs_grad,
    picks=["grad"],
    show=False,
)
# Interactive debug version (uncomment when needed):
# mne.viz.set_browser_backend("qt")
# epochs[plot_trials_condition].plot(
#     n_epochs=plot_trials_n_epochs_grad,
#     picks=["grad"],
#     show=True,
#     block=True,
# )
report.add_figure(
    fig=fig_trials_grad,
    title=f"Trials ({plot_trials_condition}) - gradiometers",
    caption=(
        "FLUX-equivalent trial browser for gradiometers. "
        f"Condition={plot_trials_condition}, n_epochs={plot_trials_n_epochs_grad}."
    ),
    section="Epoching",
    tags=("epochs", "qc"),
)
savefig_logged(fig_trials_grad, "04_epoching_trials_grad.png", dpi=150, bbox_inches="tight")
plt.close(fig_trials_grad)

# 2) Trial browser plot for stimulus channels (non-interactive, saved + report).
fig_trials_stim = epochs[plot_trials_condition].plot(
    n_epochs=plot_trials_n_epochs_stim,
    picks=["stim"],
    show=False,
)
# Interactive debug version (uncomment when needed):
# mne.viz.set_browser_backend("qt")
# epochs[plot_trials_condition].plot(
#     n_epochs=plot_trials_n_epochs_stim,
#     picks=["stim"],
#     show=True,
#     block=True,
# )
report.add_figure(
    fig=fig_trials_stim,
    title=f"Trials ({plot_trials_condition}) - stimulus channels",
    caption=(
        "FLUX-equivalent trial browser for trigger channels. "
        f"Condition={plot_trials_condition}, n_epochs={plot_trials_n_epochs_stim}."
    ),
    section="Epoching",
    tags=("epochs", "qc"),
)
savefig_logged(fig_trials_stim, "04_epoching_trials_stim.png", dpi=150, bbox_inches="tight")
plt.close(fig_trials_stim)

# 3) Trial image plot on representative channel (non-interactive, saved + report).
if plot_trials_image_channel not in epochs.ch_names:
    grad_pick = mne.pick_types(epochs.info, meg="grad")
    if len(grad_pick) == 0:
        raise RuntimeError(
            f"Channel '{plot_trials_image_channel}' not found and no grad channel available."
        )
    plot_trials_image_channel = epochs.ch_names[grad_pick[0]]

epochs_img = epochs[plot_trials_condition].copy().filter(0.0, 30).crop(-1, 2.667)
fig_trials_image = epochs_img.plot_image(
    picks=[plot_trials_image_channel],
    vmin=-150,
    vmax=150,
    show=False,
)
if isinstance(fig_trials_image, list):
    fig_trials_image = fig_trials_image[0]
# Interactive debug version (uncomment when needed):
# epochs[plot_trials_condition].copy().filter(0.0, 30).crop(-1, 2.667).plot_image(
#     picks=[plot_trials_image_channel],
#     vmin=-150,
#     vmax=150,
#     show=True,
# )
report.add_figure(
    fig=fig_trials_image,
    title=f"Trial image ({plot_trials_condition}) - {plot_trials_image_channel}",
    caption=(
        "FLUX-equivalent image plot (low-pass 30 Hz, crop -0.1..0.4 s). "
        f"Condition={plot_trials_condition}, channel={plot_trials_image_channel}."
    ),
    section="Epoching",
    tags=("epochs", "qc"),
)
savefig_logged(fig_trials_image, "04_epoching_trials_image.png", dpi=150, bbox_inches="tight")
plt.close(fig_trials_image)

log_to_console_and_report(
    "Generated FLUX-equivalent 'Plotting the trials' outputs.\n"
    f"plot_trials_condition={plot_trials_condition}\n"
    f"plot_trials_n_epochs_grad={plot_trials_n_epochs_grad}\n"
    f"plot_trials_n_epochs_stim={plot_trials_n_epochs_stim}\n"
    f"plot_trials_image_channel={plot_trials_image_channel}\n"
    "saved_figures=['04_epoching_trials_grad.png', "
    "'04_epoching_trials_stim.png', '04_epoching_trials_image.png']",
    title="Trial plots",
    section="Epoching",
)

fig_drop_log = epochs.plot_drop_log(show=False)
report.add_figure(
    fig=fig_drop_log,
    title="Epoch drop log",
    section="Epoching",
    tags=("epochs", "qc"),
)
savefig_logged(fig_drop_log, "04_epoching_drop_log.png", dpi=150, bbox_inches="tight")
plt.close(fig_drop_log)

# %%
# -------------------------------------------------------------------------
# Save outputs: epochs FIF + HTML report
# -------------------------------------------------------------------------
epochs.save(out_epochs, overwrite=True)
if epochs.metadata is None or "trial_sequence_identity" not in epochs.metadata:
    raise RuntimeError("Epoch metadata was not saved with trial_sequence_identity.")
epochs.metadata.to_csv(out_epoch_metadata_csv, index=False)
log_to_console_and_report(
    "Saved epoch output.\n"
    f"epochs_fif={out_epochs}\n"
    f"epoch_metadata_csv={out_epoch_metadata_csv}",
    title="Epoch output",
    section="Outputs",
)

report.save(out_report, overwrite=True, open_browser=False)
print(f"Saved HTML report to: {out_report}")
