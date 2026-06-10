#!/usr/bin/env python3
"""
03_C_ICA_auto.py

Purpose
Provide an automatic ICA-cleaning alternative to `03_C_ICA.py` using MEGNet
(`mne-icalabel`) as the primary component-labeling method for MEG.

This script follows the local `preprocessing_basic` flow and keeps the key data
flow of the FLUX ICA stage, while replacing manual component rejection with a
transparent rule-based automatic policy.

High-level flow
1) load annotation-ready FIF files (`*_ann.fif`) from step 03_B,
2) build an ICA-fit copy (MEG-only, 1-100 Hz band-pass, 50/60 Hz notch,
   250 Hz resampling),
3) fit ICA with Infomax (or load cached pre-fit ICA if available),
4) label components with MEGNet (`mne_icalabel`),
5) auto-exclude only targeted artefact labels above threshold,
6) apply ICA to the original annotated raws (not the fit copy),
7) preserve input bad-channel metadata in each cleaned output header, and
8) save cleaned outputs, ICA caches, labeling CSV, figures, and HTML report.

Usage example
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA_auto.py

Inputs
- Annotated FIF files from step 03_B:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/sss/*_ann.fif

Outputs
- ICA-cleaned FIF files:
  .../sss/*_ica_auto.fif
  (`raw.info['bads']` preserved from corresponding `*_ann.fif` input)
- Saved ICA model files:
  .../ica_models/sub-01_ses-01_ica_auto_prefit.fif
  .../ica_models/sub-01_ses-01_ica_auto_postreject.fif
- Component-label summary CSV:
  .../ica_models/sub-01_ses-01_ica_auto_megnet_labels.csv
- Static figures:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/figures/03_C_ICA_auto_*.png
- HTML report:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/03_C_ICA_auto_report.html

Automatic exclusion policy
- Target labels: `eye blink`, `eye movement`, `heart beat`.
- A component is auto-excluded only if:
  predicted_label in target_labels AND predicted_label_probability >= threshold.
- `brain/other` is never auto-excluded by this script.

Decision on prior subject-01 blink workaround
- The manual subject-specific hard-coded ICA-index workaround used in
  `03_C_ICA.py` is intentionally NOT reused here.
- In this auto script, blink/eye rejection is driven only by MEGNet labels and
  the probability threshold; no hard-coded component index is applied.
- For subject 01 only, blink windows are additionally annotated from the
  MEGNet-identified `eye blink` component (if available).
"""

# %%
# -------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------
import csv
import html
import os
import os.path as op
from collections import Counter

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import ICA, read_ica

# Fail early with explicit install guidance if MEGNet runtime deps are missing.
try:
    import onnxruntime as _ort  # noqa: F401
except Exception as exc:
    raise RuntimeError(
        "Missing dependency: onnxruntime is required for MEGNet inference. "
        "Install with: pip install onnxruntime"
    ) from exc

try:
    from mne_icalabel import label_components
    from mne_icalabel.megnet import megnet_label_components
except Exception as exc:
    raise RuntimeError(
        "Missing dependency: mne-icalabel is required for MEGNet labeling. "
        "Install with: pip install mne-icalabel"
    ) from exc

# %%
# -------------------------------------------------------------------------
# Subject/session selection and key editable constants
# -------------------------------------------------------------------------
subj_str = "01"
ses_str = "01"

# ICA-fit preparation settings aligned with MEGNet expectations.
FIT_L_FREQ = 1.0
FIT_H_FREQ = 100.0
FIT_RESAMPLE_SFREQ = 250.0
FIT_DEFAULT_NOTCH_FREQ = 50.0
FIT_VALID_NOTCH_FREQS = (50.0, 60.0)

# ICA fitting settings (MEGNet guidance prefers infomax).
ICA_METHOD = "infomax"
ICA_RANDOM_STATE = 97
ICA_N_COMPONENTS = 20
ICA_MAX_ITER = "auto"

# Cache control for expensive ICA fitting.
REUSE_PREFIT_ICA = True

# Automatic exclusion policy settings.
AUTO_EXCLUDE_LABELS = ("eye blink", "eye movement", "heart beat")
AUTO_EXCLUDE_PROBA_THRESHOLD = 0.80
# Subject-specific behavior: add blink annotations from MEGNet `eye blink` IC.
USE_MEGNET_BLINK_ANNOTATION_FOR_SUBJECT = subj_str == "01"

# Expected MEGNet class order from mne_icalabel.megnet docs.
MEGNET_CLASS_ORDER = ("brain/other", "eye movement", "heart beat", "eye blink")

# %%
# -------------------------------------------------------------------------
# Paths and output targets
# -------------------------------------------------------------------------
sss_dir = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"preprocessing/sub-{subj_str}/ses-{ses_str}/sss"
)
figure_output_dir = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"preprocessing/sub-{subj_str}/ses-{ses_str}/figures"
)
ica_model_output_dir = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"preprocessing/sub-{subj_str}/ses-{ses_str}/ica_models"
)
out_report = (
    "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/"
    f"reports/sub-{subj_str}/ses-{ses_str}/03_C_ICA_auto_report.html"
)
ica_prefit_file = op.join(
    ica_model_output_dir,
    f"sub-{subj_str}_ses-{ses_str}_ica_auto_prefit.fif",
)
ica_postreject_file = op.join(
    ica_model_output_dir,
    f"sub-{subj_str}_ses-{ses_str}_ica_auto_postreject.fif",
)
labels_summary_csv = op.join(
    ica_model_output_dir,
    f"sub-{subj_str}_ses-{ses_str}_ica_auto_megnet_labels.csv",
)

os.makedirs(figure_output_dir, exist_ok=True)
os.makedirs(ica_model_output_dir, exist_ok=True)
os.makedirs(op.dirname(out_report), exist_ok=True)


# -------------------------------------------------------------------------
# Shared helper functions for figure/report/CSV logging
# -------------------------------------------------------------------------
def savefig_logged(fig, filename, **kwargs):
    """
    Save a matplotlib figure into the shared figure output folder and log path.
    """
    out_path = op.join(figure_output_dir, filename)
    fig.savefig(out_path, **kwargs)
    print(f"Saved figure to: {out_path}")
    return out_path


report = mne.Report(title=f"sub-{subj_str} | 03_C ICA auto (MEGNet)")


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


def infer_notch_frequency(raw, *, fallback=FIT_DEFAULT_NOTCH_FREQ):
    """
    Infer notch frequency from raw.info['line_freq'], fallback to a default.
    """
    line_freq = raw.info.get("line_freq")
    if line_freq in FIT_VALID_NOTCH_FREQS:
        return float(line_freq), "raw.info['line_freq']"
    return float(fallback), "fallback_default"


def prob_col_name(class_label):
    """
    Convert class label to stable CSV column name.
    """
    return "prob_" + class_label.replace("/", "_").replace(" ", "_")


def write_label_summary_csv(csv_path, rows):
    """
    Write per-component MEGNet labeling summary to CSV.
    """
    if not rows:
        raise RuntimeError("No MEGNet label rows available to write CSV summary.")
    fieldnames = list(rows[0].keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved label summary CSV to: {csv_path}")


def build_component_blink_annotations(ica_model, raw_fit_like, component_idx):
    """
    Build blink annotations by running blink-event detection on one ICA source.

    Parameters
    ----------
    ica_model : mne.preprocessing.ICA
        Fitted ICA model.
    raw_fit_like : mne.io.BaseRaw
        Raw instance preprocessed like the ICA fit copy (MEG, 1-100, notch,
        resampled to 250 Hz).
    component_idx : int
        ICA component index used as proxy blink channel.

    Returns
    -------
    annotations_blink : mne.Annotations
        Blink annotations centered around detected events.
    n_blinks : int
        Number of detected blink-like events.
    """
    sources = ica_model.get_sources(raw_fit_like)
    component_name = sources.ch_names[component_idx]
    blink_events = mne.preprocessing.find_eog_events(
        sources,
        ch_name=component_name,
        verbose="ERROR",
    )
    n_blinks = len(blink_events)
    onset = blink_events[:, 0] / raw_fit_like.info["sfreq"] - 0.25
    duration = np.repeat(0.5, n_blinks)
    description = ["blink"] * n_blinks
    annotations_blink = mne.Annotations(
        onset,
        duration,
        description,
        raw_fit_like.info["meas_date"],
    )
    return annotations_blink, n_blinks


# %%
# -------------------------------------------------------------------------
# Input discovery: load annotation FIFs and derive auto-ICA output filenames
# -------------------------------------------------------------------------
ann_files = sorted(
    op.join(sss_dir, fname)
    for fname in os.listdir(sss_dir)
    if fname.endswith("_ann.fif")
)
if not ann_files:
    raise FileNotFoundError(
        f"No input *_ann.fif files found under: {sss_dir}. Run 03_B first."
    )

ica_auto_files = [ann_path.replace("_ann.fif", "_ica_auto.fif") for ann_path in ann_files]

log_to_console_and_report(
    "Discovered inputs:\n"
    + "\n".join(f"  - {p}" for p in ann_files)
    + "\n\nPlanned cleaned outputs:\n"
    + "\n".join(f"  - {p}" for p in ica_auto_files)
    + "\n\nICA model outputs:\n"
    + f"  - prefit: {ica_prefit_file}\n"
    + f"  - postreject: {ica_postreject_file}\n"
    + f"  - label summary: {labels_summary_csv}\n"
    + "\nPolicy notes:\n"
    + "  - Manual subject-specific blink workaround is NOT used here.\n"
    + f"  - Auto-exclude labels: {AUTO_EXCLUDE_LABELS}\n"
    + f"  - Probability threshold: {AUTO_EXCLUDE_PROBA_THRESHOLD:.2f}",
    title="Input/output overview",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# ICA-fit copy preparation: MEG-only, 1-100 Hz, notch 50/60, resample 250 Hz
# -------------------------------------------------------------------------
raw_fit_all = None
fit_input_summaries = []
selected_notch_freq = None
selected_notch_source = None

for ann_path in ann_files:
    raw = mne.io.read_raw_fif(ann_path, preload=True, verbose=True)
    if selected_notch_freq is None:
        selected_notch_freq, selected_notch_source = infer_notch_frequency(raw)

    # Keep the fit copy strictly MEG-only and MEGNet-compatible preprocessing.
    raw_fit = raw.copy().pick("meg")
    raw_fit.filter(FIT_L_FREQ, FIT_H_FREQ, verbose=True)
    raw_fit.notch_filter([selected_notch_freq], verbose=True)
    raw_fit.resample(FIT_RESAMPLE_SFREQ)

    if raw_fit_all is None:
        raw_fit_all = mne.io.concatenate_raws([raw_fit])
    else:
        # Runs may differ in dev_head_t; continue with explicit warning handling.
        raw_fit_all = mne.io.concatenate_raws(
            [raw_fit_all, raw_fit],
            on_mismatch="warn",
        )

    fit_input_summaries.append(
        f"{op.basename(ann_path)}: n_channels={raw_fit.info['nchan']}, "
        f"n_times={raw_fit.n_times}, sfreq={raw_fit.info['sfreq']}"
    )

log_to_console_and_report(
    "ICA fit input summary:\n"
    + "\n".join(f"  - {line}" for line in fit_input_summaries)
    + "\n\nFit-copy preprocessing settings:\n"
    + f"  - filter: {FIT_L_FREQ:.1f}-{FIT_H_FREQ:.1f} Hz\n"
    + f"  - notch: {selected_notch_freq:.1f} Hz ({selected_notch_source})\n"
    + f"  - resample: {FIT_RESAMPLE_SFREQ:.1f} Hz\n"
    + f"  - concatenated n_channels={raw_fit_all.info['nchan']}, "
    f"n_times={raw_fit_all.n_times}, sfreq={raw_fit_all.info['sfreq']}",
    title="ICA fit input summary",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# ICA fit/caching: load cached prefit model or fit new Infomax decomposition
# -------------------------------------------------------------------------
if REUSE_PREFIT_ICA and op.exists(ica_prefit_file):
    ica_cached = read_ica(ica_prefit_file, verbose=True)
    if int(ica_cached.n_components_) == int(ICA_N_COMPONENTS):
        ica = ica_cached
        log_to_console_and_report(
            "Loaded cached pre-fit ICA model.\n"
            f"  - {ica_prefit_file}\n"
            f"method={ica.method}\n"
            f"n_components_={ica.n_components_}",
            title="ICA fit summary",
            section="ICA decomposition",
        )
    else:
        log_to_console_and_report(
            "Cached pre-fit ICA ignored because component count does not match "
            "current configuration.\n"
            f"  - cache_file={ica_prefit_file}\n"
            f"  - cache_n_components={ica_cached.n_components_}\n"
            f"  - requested_n_components={ICA_N_COMPONENTS}\n"
            "Proceeding with a fresh fit.",
            title="ICA fit summary",
            section="ICA decomposition",
        )
        REUSE_PREFIT_ICA = False

if not (REUSE_PREFIT_ICA and op.exists(ica_prefit_file) and "ica" in locals()):
    try:
        ica = ICA(
            method=ICA_METHOD,
            random_state=ICA_RANDOM_STATE,
            n_components=ICA_N_COMPONENTS,
            max_iter=ICA_MAX_ITER,
            verbose=True,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Unable to initialize ICA with method='{ICA_METHOD}'. "
            "Check local MNE version and ICA settings."
        ) from exc

    ica.fit(raw_fit_all, verbose=True)
    ica.save(ica_prefit_file, overwrite=True)
    log_to_console_and_report(
        "Fitted ICA model and saved pre-fit cache.\n"
        f"method={ica.method}\n"
        f"n_components_={ica.n_components_}\n"
        f"random_state={ICA_RANDOM_STATE}\n"
        f"saved_prefit={ica_prefit_file}",
        title="ICA fit summary",
        section="ICA decomposition",
    )

# %%
# -------------------------------------------------------------------------
# ICA visual diagnostics: source traces and component maps (if feasible)
# -------------------------------------------------------------------------
try:
    # Use show=False to avoid blocking GUI main loop in unattended script runs.
    fig_sources = ica.plot_sources(raw_fit_all, title="ICA source traces (fit copy)", show=False)
    report.add_figure(
        fig=fig_sources,
        title="ICA source traces (fit copy)",
        section="ICA decomposition",
        tags=("ica", "qc"),
    )
    savefig_logged(fig_sources, "03_C_ICA_auto_sources.png", dpi=150, bbox_inches="tight")
    plt.close(fig_sources)
except Exception as exc:
    log_to_console_and_report(
        f"WARNING: could not render ICA source traces figure: {exc}",
        title="ICA source figure warning",
        section="ICA decomposition",
    )

try:
    # Use show=False to keep plotting non-blocking in script mode.
    fig_components = ica.plot_components(show=False)
    if not isinstance(fig_components, list):
        fig_components = [fig_components]
    for idx, fig_comp in enumerate(fig_components, start=1):
        report.add_figure(
            fig=fig_comp,
            title=f"ICA component topographies ({idx}/{len(fig_components)})",
            section="ICA decomposition",
            tags=("ica", "qc"),
        )
        savefig_logged(
            fig_comp,
            f"03_C_ICA_auto_components_{idx:02d}.png",
            dpi=150,
            bbox_inches="tight",
        )
        plt.close(fig_comp)
except Exception as exc:
    log_to_console_and_report(
        f"WARNING: could not render ICA component topography figures: {exc}",
        title="ICA component figure warning",
        section="ICA decomposition",
    )

# %%
# -------------------------------------------------------------------------
# MEGNet labeling: stable label API + full per-class probabilities
# -------------------------------------------------------------------------
label_result = label_components(raw_fit_all, ica, method="megnet")
labels = list(label_result["labels"])
predicted_label_proba = np.asarray(label_result["y_pred_proba"], dtype=float)

# `label_components` provides labels + selected class probability; we also keep
# full class probabilities for transparent per-class reporting in CSV.
proba_matrix = np.asarray(megnet_label_components(raw_fit_all, ica), dtype=float)
if proba_matrix.ndim != 2 or proba_matrix.shape[0] != len(labels):
    raise RuntimeError(
        "Unexpected MEGNet probability matrix shape. "
        f"Received shape={proba_matrix.shape}, n_labels={len(labels)}."
    )

if proba_matrix.shape[1] == len(MEGNET_CLASS_ORDER):
    class_order = list(MEGNET_CLASS_ORDER)
else:
    # Defensive fallback if class dimensionality changes in future releases.
    class_order = [f"class_{idx}" for idx in range(proba_matrix.shape[1])]

rows = []
auto_exclude = []
for comp_idx, label in enumerate(labels):
    pred_prob = float(predicted_label_proba[comp_idx])
    is_target_label = label in AUTO_EXCLUDE_LABELS
    meets_threshold = pred_prob >= AUTO_EXCLUDE_PROBA_THRESHOLD
    auto_flag = bool(is_target_label and meets_threshold)
    if auto_flag:
        auto_exclude.append(comp_idx)

    if auto_flag:
        exclusion_rule = (
            "label_in_target_set_and_predicted_label_probability_at_or_above_threshold"
        )
    elif not is_target_label:
        exclusion_rule = "label_not_in_target_set"
    else:
        exclusion_rule = "predicted_label_probability_below_threshold"

    row = {
        "component": comp_idx,
        "predicted_label": label,
        "predicted_label_probability": pred_prob,
        "auto_excluded": int(auto_flag),
        "exclusion_rule": exclusion_rule,
        "threshold_used": AUTO_EXCLUDE_PROBA_THRESHOLD,
    }
    for col_idx, class_name in enumerate(class_order):
        row[prob_col_name(class_name)] = float(proba_matrix[comp_idx, col_idx])
    rows.append(row)

write_label_summary_csv(labels_summary_csv, rows)

label_counts = Counter(labels)
ica.exclude = sorted(auto_exclude)

# Subject-specific blink-window annotation from MEGNet-labeled blink component.
blink_component_idx = None
if USE_MEGNET_BLINK_ANNOTATION_FOR_SUBJECT:
    blink_candidates = [idx for idx, label in enumerate(labels) if label == "eye blink"]
    if blink_candidates:
        blink_component_idx = int(
            max(blink_candidates, key=lambda idx: float(predicted_label_proba[idx]))
        )
        try:
            annotations_blink_ref, n_blinks_ref = build_component_blink_annotations(
                ica,
                raw_fit_all,
                blink_component_idx,
            )
            log_to_console_and_report(
                "Subject-specific blink annotation enabled (sub-01).\n"
                f"MEGNet blink component selected: ICA{blink_component_idx:03d}\n"
                f"Reference blink windows detected on fit-copy data: {n_blinks_ref}",
                title="Blink annotation from MEGNet",
                section="ICA labeling",
            )
            del annotations_blink_ref
        except Exception as exc:
            log_to_console_and_report(
                "WARNING: could not build reference blink annotations from "
                f"ICA{blink_component_idx:03d}: {exc}\n"
                "Blink-window injection into outputs will be skipped.",
                title="Blink annotation from MEGNet",
                section="ICA labeling",
            )
            blink_component_idx = None
    else:
        log_to_console_and_report(
            "Subject-specific blink annotation enabled, but no component was "
            "labeled as `eye blink` by MEGNet. Blink-window injection skipped.",
            title="Blink annotation from MEGNet",
            section="ICA labeling",
        )

ica.save(ica_postreject_file, overwrite=True)

log_to_console_and_report(
    "MEGNet labeling summary:\n"
    + "\n".join(f"  - {label}: {count}" for label, count in sorted(label_counts.items()))
    + "\n\nAuto-exclusion policy:\n"
    + f"  - target_labels={AUTO_EXCLUDE_LABELS}\n"
    + f"  - probability_threshold={AUTO_EXCLUDE_PROBA_THRESHOLD:.2f}\n"
    + f"  - excluded_components={ica.exclude}\n"
    + f"  - blink_component_for_sub01={blink_component_idx}\n"
    + "\nSaved outputs:\n"
    + f"  - labels_csv={labels_summary_csv}\n"
    + f"  - postreject_ica={ica_postreject_file}",
    title="MEGNet labeling summary",
    section="ICA labeling",
)

# %%
# -------------------------------------------------------------------------
# Apply ICA to original annotated raws and save cleaned *_ica_auto.fif outputs
# -------------------------------------------------------------------------
processed_pairs = []
for ann_path, out_path in zip(ann_files, ica_auto_files, strict=True):
    raw_ica = mne.io.read_raw_fif(ann_path, preload=True, verbose=True)
    input_bads = list(dict.fromkeys(str(ch) for ch in raw_ica.info.get("bads", [])))
    # For sub-01 only, add blink windows from MEGNet blink component before save.
    if blink_component_idx is not None:
        raw_fit_like_run = raw_ica.copy().pick("meg")
        raw_fit_like_run.filter(FIT_L_FREQ, FIT_H_FREQ, verbose=False)
        raw_fit_like_run.notch_filter([selected_notch_freq], verbose=False)
        raw_fit_like_run.resample(FIT_RESAMPLE_SFREQ)
        annotations_blink_run, n_blinks_run = build_component_blink_annotations(
            ica,
            raw_fit_like_run,
            blink_component_idx,
        )
        raw_ica.set_annotations(raw_ica.annotations + annotations_blink_run)
        print(
            f"Added {n_blinks_run} blink annotations from ICA{blink_component_idx:03d} "
            f"to {op.basename(ann_path)}"
        )
    ica.apply(raw_ica)
    # Preserve upstream bad-channel metadata after ICA projection.
    raw_ica.info["bads"] = list(input_bads)
    raw_ica.save(out_path, overwrite=True)
    processed_pairs.append((ann_path, out_path, input_bads))

log_to_console_and_report(
    "Applied ICA to original annotated raws:\n"
    + "\n".join(
        f"  - {src} -> {dst} (propagated_bads={bads})"
        for src, dst, bads in processed_pairs
    ),
    title="ICA-applied outputs",
    section="Apply ICA",
)

# %%
# -------------------------------------------------------------------------
# Before/after comparison figure pair (FLUX-style channel-level inspection)
# -------------------------------------------------------------------------
compare_chs = ["MEG0311", "MEG0121", "MEG1211", "MEG1411"]
raw_before = mne.io.read_raw_fif(ann_files[0], preload=True, verbose=False)
raw_after = mne.io.read_raw_fif(ica_auto_files[0], preload=True, verbose=False)

available_compare_chs = [ch for ch in compare_chs if ch in raw_before.ch_names]
if not available_compare_chs:
    available_compare_chs = [ch for ch in raw_before.ch_names if ch.startswith("MEG")][:4]
chan_idxs = [raw_before.ch_names.index(ch) for ch in available_compare_chs]

fig_before = raw_before.plot(order=chan_idxs, duration=5, show=False)
report.add_figure(
    fig=fig_before,
    title=f"Before ICA auto ({op.basename(ann_files[0])})",
    caption=f"Channels: {available_compare_chs}",
    section="Before/after ICA",
    tags=("raw", "qc"),
)
savefig_logged(fig_before, "03_C_ICA_auto_before_projection.png", dpi=150, bbox_inches="tight")
plt.close(fig_before)

fig_after = raw_after.plot(order=chan_idxs, duration=5, show=False)
report.add_figure(
    fig=fig_after,
    title=f"After ICA auto ({op.basename(ica_auto_files[0])})",
    caption=f"Channels: {available_compare_chs}",
    section="Before/after ICA",
    tags=("raw", "qc"),
)
savefig_logged(fig_after, "03_C_ICA_auto_after_projection.png", dpi=150, bbox_inches="tight")
plt.close(fig_after)

log_to_console_and_report(
    "Before/after comparison summary:\n"
    + f"  - before={ann_files[0]}\n"
    + f"  - after={ica_auto_files[0]}\n"
    + f"  - channels={available_compare_chs}",
    title="Before/after comparison",
    section="Before/after ICA",
)

# %%
# -------------------------------------------------------------------------
# Final report save and output path summary
# -------------------------------------------------------------------------
output_lines = [
    "Final output paths:",
    f"  - report={out_report}",
    f"  - prefit_ica={ica_prefit_file}",
    f"  - postreject_ica={ica_postreject_file}",
    f"  - labels_csv={labels_summary_csv}",
]
output_lines.extend(f"  - cleaned_fif={path}" for path in ica_auto_files)
log_to_console_and_report(
    "\n".join(output_lines),
    title="Output paths",
    section="Outputs",
)

report.save(out_report, overwrite=True, open_browser=False)
print(f"Saved HTML report to: {out_report}")
