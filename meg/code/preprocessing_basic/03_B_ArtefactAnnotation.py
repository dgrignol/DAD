#!/usr/bin/env python3
"""
03_B_ArtefactAnnotation.py

Purpose
Adapt the FLUX MNE-Python ArtefactAnnotation workflow to the local DAD MEG
dataset with minimal changes. This script:
1) loads Maxwell-filtered FIF files for sub-01/ses-01,
2) detects eye-blink artefacts from EOG001 in a consistent sample frame
   (except subject 01: blink annotation is intentionally skipped here),
3) detects muscle artefacts from magnetometers (110-140 Hz z-score),
4) combines blink/muscle annotations with existing annotations,
5) saves annotated FIF + CSV files while preserving input bad-channel metadata,
6) saves static QC figures, and
7) writes an HTML MNE report with figures and key print logs.

Usage example
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py

Inputs
- Maxwell-filtered FIF files:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/sss/*_raw_sss.fif

Outputs
- Annotated FIF files next to SSS files:
  .../sss/*_ann.fif
  (`raw.info['bads']` preserved from corresponding `*_raw_sss.fif` input)
- Annotation CSV files:
  .../sss/*_ann.csv
- QC figures:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/figures/
- HTML report:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/03_B_ArtefactAnnotation_report.html
"""

# %%
# -------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------
import html
import os
import os.path as op

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import annotate_muscle_zscore

# %%
# -------------------------------------------------------------------------
# Subject/session selection and output layout
# -------------------------------------------------------------------------
subj_str = "01"
ses_str = "01"
# Subject-specific exception: for sub-01 we do not trust EOG blink labels.
# Blink annotations are intentionally deferred to ICA component-0 logic in 03_C.
skip_eog_blink_for_subject = subj_str == "01"

sss_input_dir = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"preprocessing/sub-{subj_str}/ses-{ses_str}/sss"
)
figure_output_dir = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"preprocessing/sub-{subj_str}/ses-{ses_str}/figures"
)
out_report = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"reports/sub-{subj_str}/ses-{ses_str}/03_B_ArtefactAnnotation_report.html"
)

os.makedirs(figure_output_dir, exist_ok=True)
os.makedirs(op.dirname(out_report), exist_ok=True)


# -------------------------------------------------------------------------
# Shared helpers for figure saving and sectioned print logging to report
# -------------------------------------------------------------------------
def savefig_logged(fig, filename, **kwargs):
    """
    Save a matplotlib figure into the shared figure folder and print the path.
    """
    out_path = op.join(figure_output_dir, filename)
    fig.savefig(out_path, **kwargs)
    print(f"Saved figure to: {out_path}")


report = mne.Report(title=f"sub-{subj_str} | 03_B ArtefactAnnotation") # initialize html report with title


def log_to_console_and_report(text, *, title, section, tags=("log", "qc")):
    """
    Print text to console and store the same text in an MNE report section.
    """
    print(text)
    report.add_html(
        f"<pre>{html.escape(text)}</pre>",
        title=title,
        section=section,
        tags=tags,
    )


def build_blink_annotations(raw, *, ch_name="EOG001"):
    """
    Detect blink events on one raw object and return aligned annotations.

    Data-flow note
    `find_eog_events` returns event sample indices in the same sample frame as
    `raw.first_samp`. To keep blink annotations aligned after `raw.set_annotations`
    (with `orig_time=raw.info['meas_date']`), we keep this absolute-sample
    frame and do *not* subtract `raw._first_time`.

    Parameters
    ----------
    raw : mne.io.BaseRaw
        Raw object containing EOG channels.
    ch_name : str
        EOG channel used for blink detection.

    Returns
    -------
    annotations_blink : mne.Annotations
        Blink annotations centered around each detected event.
    n_blinks : int
        Number of detected blink events.
    n_out_of_range : int
        Number of blink onset samples outside [first_samp, last_samp].
    """
    eog_events = mne.preprocessing.find_eog_events(raw, ch_name=ch_name)
    n_blinks = len(eog_events)
    onset = eog_events[:, 0] / raw.info["sfreq"] - 0.25
    duration = np.repeat(0.5, n_blinks)
    description = ["blink"] * n_blinks
    annotations_blink = mne.Annotations(
        onset,
        duration,
        description,
        raw.info["meas_date"],
    )
    # Validate in sample space (absolute frame), not in raw.times (relative frame).
    onset_samples = np.round(onset * raw.info["sfreq"]).astype(int)
    sample_min = raw.first_samp
    sample_max = raw.first_samp + raw.n_times - 1
    n_out_of_range = int(np.sum((onset_samples < sample_min) | (onset_samples > sample_max)))
    return annotations_blink, n_blinks, n_out_of_range


# %%
# -------------------------------------------------------------------------
# FLUX "Set the path and load the FIF-files" section adapted to local SSS files
# -------------------------------------------------------------------------
sss_files = sorted(
    op.join(sss_input_dir, fname)
    for fname in os.listdir(sss_input_dir)
    if fname.endswith("_raw_sss.fif")
)
if not sss_files:
    raise FileNotFoundError(
        f"No SSS FIF files found under: {sss_input_dir}. Run 03_A_MaxFilt.py first."
    )

raw1_path = sss_files[0]
ann1_fif = raw1_path.replace("_raw_sss.fif", "_ann.fif")
ann1_csv = ann1_fif.replace(".fif", ".csv")

log_to_console_and_report(
    "Input SSS files discovered:\n"
    + "\n".join(f"  - {p}" for p in sss_files)
    + f"\n\nFirst-file outputs:\n  - {ann1_fif}\n  - {ann1_csv}",
    title="Input/output overview",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# FLUX run-01 equivalent: load first SSS file with preload=True
# -------------------------------------------------------------------------
raw1 = mne.io.read_raw_fif(raw1_path, preload=True, verbose=True)
raw1_input_bads = list(dict.fromkeys(str(ch) for ch in raw1.info.get("bads", [])))
log_to_console_and_report(
    f"Loaded first file:\n  - {raw1_path}\nRaw summary: {raw1}",
    title="Loaded first file",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# FLUX "Identifying eye blinks": detect blinks from EOG001
# -------------------------------------------------------------------------
if skip_eog_blink_for_subject:
    annotations_blink = mne.Annotations([], [], [], raw1.info["meas_date"])
    n_blinks = 0
    n_out_of_range = 0
else:
    annotations_blink, n_blinks, n_out_of_range = build_blink_annotations(raw1, ch_name="EOG001")

log_to_console_and_report(
    (
        "WARNING: blink annotation skipped in 03_B for sub-01 because "
        "EOG-based blink labels are unreliable in this subject. "
        "Blink annotation is deferred to ICA component-0 logic in 03_C.\n"
        if skip_eog_blink_for_subject
        else ""
    )
    + f"Detected blink events (block01 equivalent): {n_blinks}\n"
    + f"Blink onset out-of-range count: {n_out_of_range}",
    title="Blink detection",
    section="Block01 annotation",
)

# %%
# -------------------------------------------------------------------------
# Blink-only QC: frontal MEG + vEOG scrollable view with highlighted blinks
# -------------------------------------------------------------------------
if skip_eog_blink_for_subject:
    log_to_console_and_report(
        "WARNING: blink-only QC figure is intentionally omitted for sub-01 in 03_B. "
        "Reason: EOG-based blink annotation is disabled for this subject and replaced "
        "by ICA component-0 blink annotation in 03_C.",
        title="Blink-only QC figure",
        section="Block01 annotation",
    )
else:
    # Build a blink-only view right after blink detection for targeted QC.
    raw_blink_qc = raw1.copy()
    raw_blink_qc.set_annotations(annotations_blink)

    # Use user-specified frontal channels (FLUX adaptation override for this dataset).
    frontal_meg_requested = [
    #    "MEG0922",
        "MEG0923",
    #    "MEG0921",
    #    "MEG1212",
        "MEG1213",
    #    "MEG1211",
    #    "MEG0512",
        "MEG0513",
    #    "MEG0511",
    #    "MEG0312",
    #    "MEG0311",
        "MEG0313",
    ]
    frontal_meg = [ch for ch in frontal_meg_requested if ch in raw_blink_qc.ch_names]
    missing_frontal_meg = [ch for ch in frontal_meg_requested if ch not in raw_blink_qc.ch_names]
    if missing_frontal_meg:
        print(f"WARNING: requested frontal channels not found: {missing_frontal_meg}")
    if not frontal_meg:
        frontal_meg = [ch for ch in raw_blink_qc.ch_names if ch.startswith("MEG")][:12]
    picks_blink_qc = frontal_meg + (["EOG001"] if "EOG001" in raw_blink_qc.ch_names else [])

    start_blink_qc = max(0.0, float(annotations_blink.onset[0]) - 1.0) if n_blinks else 0.0
    fig_blink_qc = raw_blink_qc.plot(
        picks=picks_blink_qc,
        n_channels=len(picks_blink_qc),
        duration=12.0,
        start=start_blink_qc,
        proj=False,
        scalings=dict(eog=500e-6),
        show=True,
    )

    # Recolor blink annotation spans to yellow with 0.4 alpha for visual salience.
    for key in list(fig_blink_qc.mne.annotation_segment_colors):
        if str(key) == "blink":
            fig_blink_qc.mne.annotation_segment_colors[key] = (255/255, 184/255, 28/255, 0.4)
    fig_blink_qc._redraw(annotations=True)

    report.add_figure(
        fig=fig_blink_qc,
        title="Blink-only frontal + vEOG view (block01)",
        caption=f"Blinks={n_blinks}; start={start_blink_qc:.2f}s; picks={picks_blink_qc}",
        section="Block01 annotation",
        tags=("blink", "eog", "qc"),
    )
    savefig_logged(
        fig_blink_qc,
        "03_B_ArtefactAnnotation_block01_blink_frontal_eog.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig_blink_qc)
    log_to_console_and_report(
        "Blink-only QC figure generated with frontal MEG + EOG001.\n"
        f"n_blinks={n_blinks}\n"
        f"start={start_blink_qc:.2f}s\n"
        f"picks={picks_blink_qc}",
        title="Blink-only QC figure",
        section="Block01 annotation",
    )

# %%
# -------------------------------------------------------------------------
# FLUX "Finding muscle artefacts": annotate high-frequency muscle bursts
# -------------------------------------------------------------------------
threshold_muscle = 10
annotations_muscle, scores_muscle = annotate_muscle_zscore(
    raw1,
    ch_type="mag",
    threshold=threshold_muscle,
    min_length_good=0.2,
    filter_freq=[110, 140],
)
log_to_console_and_report(
    "Muscle annotation parameters:\n"
    f"  - threshold_muscle={threshold_muscle}\n"
    "  - ch_type=mag\n"
    "  - min_length_good=0.2\n"
    "  - filter_freq=[110, 140]\n"
    f"Detected muscle segments: {len(annotations_muscle)}",
    title="Muscle detection",
    section="Block01 annotation",
)

# %%
# -------------------------------------------------------------------------
# FLUX z-score threshold check figure
# -------------------------------------------------------------------------
fig_muscle_scores, ax = plt.subplots(figsize=(10, 4))
ax.plot(raw1.times, scores_muscle)
ax.axhline(y=threshold_muscle, color="r")
ax.set(
    xlabel="time, (s)",
    ylabel="zscore",
    title=f"Muscle activity (threshold = {threshold_muscle})",
)
report.add_figure(
    fig=fig_muscle_scores,
    title="Muscle z-score over time (block01)",
    section="Block01 annotation",
    tags=("muscle", "qc"),
)
savefig_logged(
    fig_muscle_scores,
    "03_B_ArtefactAnnotation_block01_muscle_zscore.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig_muscle_scores)

# %%
# -------------------------------------------------------------------------
# FLUX "Include annotations in dataset and inspect": combine annotations
# -------------------------------------------------------------------------
annotations_event = raw1.annotations
raw1.set_annotations(annotations_event + annotations_blink + annotations_muscle)
log_to_console_and_report(
    f"Combined annotations in first file: {len(raw1.annotations)} total segments.",
    title="Combined annotations",
    section="Block01 annotation",
)

# %%
# -------------------------------------------------------------------------
# FLUX inspection plot: raw browser view with annotations
# -------------------------------------------------------------------------
# Keep FLUX-like browser defaults here; avoid expensive global auto-scaling pass.
fig_annotated_raw = raw1.plot(start=50, duration=10, proj=False, show=True)
report.add_figure(
    fig=fig_annotated_raw,
    title="Annotated raw first-look (block01)",
    section="Block01 annotation",
    tags=("raw", "qc"),
)
savefig_logged(
    fig_annotated_raw,
    "03_B_ArtefactAnnotation_block01_annotated_raw.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig_annotated_raw)

# %%
# -------------------------------------------------------------------------
# FLUX ocular-focused inspection plot: EOG-only channels
# -------------------------------------------------------------------------
# Keep FLUX channel-type assignment step for parity; these channels are usually
# already typed as EOG in this dataset, so this is idempotent.
raw1.set_channel_types({"EOG001": "eog"})
raw1.set_channel_types({"EOG002": "eog"})

eog_picks = mne.pick_types(raw1.info, meg=False, eog=True)
scl = dict(eog=500e-6)
fig_eog = raw1.plot(order=eog_picks, scalings=scl, start=50, duration=10, show=True)
report.add_figure(
    fig=fig_eog,
    title="EOG channels with blink annotations (block01)",
    section="Block01 annotation",
    tags=("eog", "qc"),
)
savefig_logged(
    fig_eog,
    "03_B_ArtefactAnnotation_block01_eog_raw.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig_eog)

# %%
# -------------------------------------------------------------------------
# FLUX save step for first file
# -------------------------------------------------------------------------
# Keep upstream bad-channel metadata intact in the annotated output file.
raw1.info["bads"] = list(raw1_input_bads)
raw1.save(ann1_fif, overwrite=True)
raw1.annotations.save(ann1_csv, overwrite=True)
log_to_console_and_report(
    "Saved first annotated files:\n"
    f"  - {ann1_fif}\n"
    f"  - {ann1_csv}\n"
    f"  - propagated_bads={raw1_input_bads}",
    title="Saved first outputs",
    section="Save outputs",
)

# %%
# -------------------------------------------------------------------------
# FLUX run-02 block generalized to all remaining SSS files for subject 01
# -------------------------------------------------------------------------
per_file_summary = []
for raw_path in sss_files[1:]:
    print(f"Annotating additional file: {raw_path}")
    raw_run = mne.io.read_raw_fif(raw_path, preload=True, verbose=True)
    raw_run_input_bads = list(dict.fromkeys(str(ch) for ch in raw_run.info.get("bads", [])))

    # Blinks: same subject-specific policy as run-01.
    if skip_eog_blink_for_subject:
        annotations_blink = mne.Annotations([], [], [], raw_run.info["meas_date"])
        n_blinks = 0
        n_out_of_range = 0
    else:
        annotations_blink, n_blinks, n_out_of_range = build_blink_annotations(
            raw_run,
            ch_name="EOG001",
        )

    # Muscle (same FLUX parameters).
    annotations_muscle, scores_muscle_run = annotate_muscle_zscore(
        raw_run,
        ch_type="mag",
        threshold=threshold_muscle,
        min_length_good=0.2,
        filter_freq=[110, 140],
    )

    # Per-file muscle z-score QC figure (saved and embedded in report).
    run_label = op.basename(raw_path).replace("_raw_sss.fif", "")
    fig_muscle_scores_run, ax = plt.subplots(figsize=(10, 4))
    ax.plot(raw_run.times, scores_muscle_run)
    ax.axhline(y=threshold_muscle, color="r")
    ax.set(
        xlabel="time, (s)",
        ylabel="zscore",
        title=f"Muscle activity ({run_label}, threshold = {threshold_muscle})",
    )
    report.add_figure(
        fig=fig_muscle_scores_run,
        title=f"Muscle z-score over time ({run_label})",
        section="Additional files annotation",
        tags=("muscle", "qc"),
    )
    savefig_logged(
        fig_muscle_scores_run,
        f"03_B_ArtefactAnnotation_{run_label}_muscle_zscore.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig_muscle_scores_run)

    # Combine and save.
    raw_run.set_annotations(raw_run.annotations + annotations_blink + annotations_muscle)
    # Preserve upstream bad-channel metadata exactly as received.
    raw_run.info["bads"] = list(raw_run_input_bads)

    out_ann_fif = raw_path.replace("_raw_sss.fif", "_ann.fif")
    out_ann_csv = out_ann_fif.replace(".fif", ".csv")
    raw_run.save(out_ann_fif, overwrite=True)
    raw_run.annotations.save(out_ann_csv, overwrite=True)

    per_file_summary.append(
        {
            "file": op.basename(raw_path),
            "blinks": n_blinks,
            "blink_out_of_range": n_out_of_range,
            "muscle_segments": len(annotations_muscle),
            "total_annotations": len(raw_run.annotations),
            "propagated_bads": list(raw_run_input_bads),
            "ann_fif": out_ann_fif,
            "ann_csv": out_ann_csv,
        }
    )

if not per_file_summary:
    log_to_console_and_report(
        "No additional SSS files found beyond the first file.",
        title="Additional files",
        section="Save outputs",
    )
else:
    summary_lines = ["Additional-file annotation summary:"]
    for item in per_file_summary:
        summary_lines.append(
            "  "
            f"{item['file']}: blinks={item['blinks']}, "
            f"blink_out_of_range={item['blink_out_of_range']}, "
            f"muscle_segments={item['muscle_segments']}, "
            f"total_annotations={item['total_annotations']}, "
            f"propagated_bads={item['propagated_bads']}"
        )
        summary_lines.append(f"    saved_fif={item['ann_fif']}")
        summary_lines.append(f"    saved_csv={item['ann_csv']}")
    log_to_console_and_report(
        "\n".join(summary_lines),
        title="Additional-file summary",
        section="Save outputs",
    )

# %%
# -------------------------------------------------------------------------
# Final report save
# -------------------------------------------------------------------------
report.save(out_report, overwrite=True, open_browser=False)
print(f"Saved HTML report to: {out_report}")
