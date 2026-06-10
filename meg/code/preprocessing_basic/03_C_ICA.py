#!/usr/bin/env python3
"""
03_C_ICA.py

Purpose
Adapt the FLUX MNE-Python ICA workflow to the local DAD MEG dataset with
minimal changes. This script follows the FLUX ICA sequence:
1) load annotation-ready FIF files (`*_ann.fif`),
2) resample MEG data to 200 Hz and band-pass filter 1-40 Hz for ICA fitting,
3) fit FastICA with 30 components (random_state=97),
4) inspect ICA sources/components and score component-EOG coupling,
5) exclude artefact components using a manual subject-specific list (or prompt
   interactively when no subject list is configured),
6) apply ICA to each annotated file,
7) preserve input bad-channel metadata in each cleaned output header, and
8) save ICA-cleaned outputs (`*_ica.fif`), figures, and HTML report.

Usage example
python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA.py

Inputs
- Annotated FIF files from step 03_B:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/sss/*_ann.fif

Outputs
- ICA-cleaned FIF files:
  .../sss/*_ica.fif
  (`raw.info['bads']` preserved from corresponding `*_ann.fif` input)
- Saved ICA model files:
  .../ica_models/sub-01_ses-01_ica_prefit.fif
  .../ica_models/sub-01_ses-01_ica_postreject.fif
- Static figures:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/figures/
- HTML report:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/reports/sub-01/ses-01/03_C_ICA_report.html

Key FLUX assumptions kept
- ICA fit is performed on MEG-only data after 200 Hz resampling and 1-40 Hz
  filtering to improve decomposition efficiency.
- ICA projection is applied to the original annotated raw files (not the
  filtered/resampled fit copy).

Local adaptation for debugging/QA
- The script now computes per-file/per-EOG ICA correlation scores and reports
  ranked ocular-candidate components. This makes it explicit that ICA index
  ordering is arbitrary and should not be interpreted as "component 1 = blinks"
  without channel-based scoring.
- Subject-specific exception for sub-01: blink annotations are generated from
  ICA component 0 (ICA000) instead of EOG-based blink labels from step 03_B.

Local adaptation note
- !!!IMPORTANT!!! runs can have different dev_head_t (different head
        # positions between blocks). Keep FLUX-style concatenated ICA fitting,
        # but allow this mismatch explicitly so fitting can proceed.
- The script caches ICA decomposition just after fitting (pre-rejection) and
  saves a second ICA file after assigning excluded components (post-rejection)
  so later reruns can skip expensive ICA fitting.
- If no manual exclusion list is configured for the current subject, the script
  saves an interim report and exits, so components can be chosen and the script
  rerun with `reuse_prefit_ica = True`.
"""

# %%
# -------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------
import html
import os
import os.path as op
from collections import defaultdict

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import ICA, read_ica

# %%
# -------------------------------------------------------------------------
# Subject/session selection and shared output paths
# -------------------------------------------------------------------------
subj_str = "01"
ses_str = "01"

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
    f"reports/sub-{subj_str}/ses-{ses_str}/03_C_ICA_report.html"
)
ica_prefit_file = op.join(
    ica_model_output_dir,
    f"sub-{subj_str}_ses-{ses_str}_ica_prefit.fif",
)
ica_postreject_file = op.join(
    ica_model_output_dir,
    f"sub-{subj_str}_ses-{ses_str}_ica_postreject.fif",
)
# If a pre-fit ICA cache exists, load it to skip expensive re-fitting.
reuse_prefit_ica = False
# Manual exclusion lists by subject.
manual_exclude_by_subject = {"01": [0, 5]}  # e.g. {"01": [1, 5]} or {"01": []} to explicitly set no exclusions; if subject not listed, script prompts for manual selection.
manual_exclude_default = []
# EOG diagnostic threshold (used only for ranking/visual diagnostics).
eog_guided_score_threshold = 0.6
eog_guided_max_components = 2
# Subject-specific exception: use ICA000-based blink annotation for sub-01.
use_ica0_blink_for_subject = subj_str == "01"

os.makedirs(figure_output_dir, exist_ok=True)
os.makedirs(ica_model_output_dir, exist_ok=True)
os.makedirs(op.dirname(out_report), exist_ok=True)


# -------------------------------------------------------------------------
# Shared helper functions for consistent figure/report logging
# -------------------------------------------------------------------------
def savefig_logged(fig, filename, **kwargs):
    """
    Save a matplotlib figure into the shared figure output folder and log path.
    """
    out_path = op.join(figure_output_dir, filename)
    fig.savefig(out_path, **kwargs)
    print(f"Saved figure to: {out_path}")


report = mne.Report(title=f"sub-{subj_str} | 03_C ICA")


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


def summarize_eog_component_scores(ica_model, ann_paths, *, top_n=5):
    """
    Compute component-vs-EOG score summaries on each annotated run.

    Parameters
    ----------
    ica_model : mne.preprocessing.ICA
        Fitted ICA instance to score.
    ann_paths : list[str]
        Paths to `*_ann.fif` files containing MEG + EOG channels.
    top_n : int
        Number of top absolute-score components reported per EOG channel.

    Returns
    -------
    summary_text : str
        Human-readable multiline summary for logs/report.
    ranked_components : list[tuple[int, float, int]]
        Component ranking as `(component_index, mean_abs_score, n_scores)`.
    auto_candidates : list[int]
        Components whose mean absolute score exceeds configured threshold.
    """
    # Collect scores per component across files/channels for robust ranking.
    component_scores = defaultdict(list)
    summary_lines = ["EOG diagnostics per file/channel (absolute-score ranking):"]
    for ann_path in ann_paths:
        raw_diag = mne.io.read_raw_fif(ann_path, preload=True, verbose=False)
        eog_names = [
            ch_name
            for ch_name, ch_type in zip(raw_diag.ch_names, raw_diag.get_channel_types(), strict=True)
            if ch_type == "eog"
        ]
        if not eog_names:
            summary_lines.append(f"  - {op.basename(ann_path)}: no EOG channels found.")
            continue

        summary_lines.append(f"  - {op.basename(ann_path)}")
        for eog_name in eog_names:
            inds, scores = ica_model.find_bads_eog(raw_diag, ch_name=eog_name, verbose=False)
            top_idx = np.argsort(np.abs(scores))[::-1][:top_n]
            top_tokens = []
            for comp_idx in top_idx:
                abs_score = float(np.abs(scores[comp_idx]))
                component_scores[int(comp_idx)].append(abs_score)
                top_tokens.append(f"ICA{int(comp_idx):03d}:{scores[comp_idx]:+.3f}")
            summary_lines.append(
                f"    * {eog_name} suggested={list(map(int, inds))[:top_n]} "
                f"top={', '.join(top_tokens)}"
            )

    ranked_components = sorted(
        (
            (comp_idx, float(np.mean(abs_vals)), len(abs_vals))
            for comp_idx, abs_vals in component_scores.items()
        ),
        key=lambda row: row[1],
        reverse=True,
    )

    ranked_tokens = [
        f"ICA{comp_idx:03d}: mean|score|={mean_abs:.3f} (n={n_vals})"
        for comp_idx, mean_abs, n_vals in ranked_components[:10]
    ]
    summary_lines.append("Aggregate EOG score ranking (top 10):")
    if ranked_tokens:
        summary_lines.extend(f"    * {token}" for token in ranked_tokens)
    else:
        summary_lines.append("    * No valid EOG scores available.")

    auto_candidates = [
        comp_idx
        for comp_idx, mean_abs, _n_vals in ranked_components
        if mean_abs >= eog_guided_score_threshold
    ][:eog_guided_max_components]
    summary_lines.append(
        "Auto-candidate components above threshold "
        f"{eog_guided_score_threshold:.2f}: {auto_candidates}"
    )
    return "\n".join(summary_lines), ranked_components, auto_candidates


def build_ica0_blink_annotations(ica_model, raw_for_ica):
    """
    Detect blink-like events from ICA component 0 on ICA-fit-like data.

    Parameters
    ----------
    ica_model : mne.preprocessing.ICA
        Fitted ICA model.
    raw_for_ica : mne.io.BaseRaw
        Raw object already prepared like ICA fit input (MEG, 200 Hz, 1-40 Hz).

    Returns
    -------
    annotations_blink : mne.Annotations
        Blink annotations with description "blink".
    n_blinks : int
        Number of detected blink-like events.
    """
    ica_sources = ica_model.get_sources(raw_for_ica)
    component_name = ica_sources.ch_names[0]
    ica_events = mne.preprocessing.find_eog_events(
        ica_sources,
        ch_name=component_name,
        verbose="ERROR",
    )
    n_blinks = len(ica_events)
    onset = ica_events[:, 0] / raw_for_ica.info["sfreq"] - 0.25
    duration = np.repeat(0.5, n_blinks)
    description = ["blink"] * n_blinks
    annotations_blink = mne.Annotations(
        onset,
        duration,
        description,
        raw_for_ica.info["meas_date"],
    )
    return annotations_blink, n_blinks


# %%
# -------------------------------------------------------------------------
# FLUX "Set local paths": discover input annotation files and derive output ica files
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

ica_files = [ann_path.replace("_ann.fif", "_ica.fif") for ann_path in ann_files]

log_to_console_and_report(
    "Discovered ICA inputs:\n"
    + "\n".join(f"  - {p}" for p in ann_files)
    + "\n\nPlanned ICA outputs:\n"
    + "\n".join(f"  - {p}" for p in ica_files)
    + "\n\nICA model cache files:\n"
    + f"  - prefit: {ica_prefit_file}\n"
    + f"  - postreject: {ica_postreject_file}",
    title="Input/output overview",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# FLUX resampling/filtering block: prepare concatenated MEG data for ICA fitting
# -------------------------------------------------------------------------
raw_resmpl_all = None
fit_input_summaries = []

for ann_path in ann_files:
    raw = mne.io.read_raw_fif(ann_path, preload=True, verbose=True)
    raw_resmpl = raw.copy().pick("meg")
    raw_resmpl.resample(200)  # FLUX: down-sample to 200 Hz before ICA.
    raw_resmpl.filter(1, 40)  # FLUX: 1-40 Hz band-pass for ICA decomposition.

    if raw_resmpl_all is None:
        raw_resmpl_all = mne.io.concatenate_raws([raw_resmpl])
    else:
        # Local adaptation: runs can have different dev_head_t (different head
        # positions between blocks). Keep FLUX-style concatenated ICA fitting,
        # but allow this mismatch explicitly so fitting can proceed.
        raw_resmpl_all = mne.io.concatenate_raws(
            [raw_resmpl_all, raw_resmpl],
            on_mismatch="warn",
        )

    fit_input_summaries.append(
        f"{op.basename(ann_path)}: n_times={raw_resmpl.n_times}, sfreq={raw_resmpl.info['sfreq']}"
    )

del raw_resmpl
log_to_console_and_report(
    "ICA fit input summary:\n"
    + "\n".join(f"  - {line}" for line in fit_input_summaries)
    + f"\nConcatenated MEG data: n_channels={raw_resmpl_all.info['nchan']}, "
    f"n_times={raw_resmpl_all.n_times}, sfreq={raw_resmpl_all.info['sfreq']}",
    title="ICA fit input",
    section="Preparation",
)

# %%
# -------------------------------------------------------------------------
# FLUX ICA fitting block: fastICA on concatenated preprocessed MEG data
# -------------------------------------------------------------------------
if reuse_prefit_ica and op.exists(ica_prefit_file):
    ica = read_ica(ica_prefit_file, verbose=True)
    log_to_console_and_report(
        "Loaded cached pre-fit ICA model.\n"
        f"File: {ica_prefit_file}\n"
        f"method={ica.method}\n"
        f"n_components_={ica.n_components_}",
        title="ICA fit",
        section="ICA decomposition",
    )
else:
    ica = ICA(
        method="fastica",
        random_state=97,
        n_components=30,
        verbose=True,
    )
    ica.fit(raw_resmpl_all, verbose=True)
    # Save decomposition before exclusions so future reruns can skip fitting.
    ica.save(ica_prefit_file, overwrite=True)
    log_to_console_and_report(
        f"ICA fitted successfully.\nmethod={ica.method}\n"
        f"n_components_={ica.n_components_}\nrandom_state=97\n"
        f"Saved pre-fit ICA model: {ica_prefit_file}",
        title="ICA fit",
        section="ICA decomposition",
    )

# %%
# -------------------------------------------------------------------------
# Subject-specific blink-annotation block: inject ICA000 blink labels for sub-01
# -------------------------------------------------------------------------
if use_ica0_blink_for_subject:
    ica0_blink_annotations_all, n_ica0_blinks_all = build_ica0_blink_annotations(ica, raw_resmpl_all)
    raw_resmpl_all.set_annotations(raw_resmpl_all.annotations + ica0_blink_annotations_all)
    log_to_console_and_report(
        "Subject-specific blink policy active (sub-01): using ICA000-based blink "
        "annotation in 03_C source inspection/apply outputs.\n"
        f"Detected ICA000 blink-like events on concatenated fit input: {n_ica0_blinks_all}",
        title="Blink annotation source",
        section="ICA decomposition",
    )
else:
    log_to_console_and_report(
        "Blink policy: using existing blink annotations from 03_B (EOG-based).",
        title="Blink annotation source",
        section="ICA decomposition",
    )

# %%
# -------------------------------------------------------------------------
# FLUX component-inspection block: source traces and topographies
# -------------------------------------------------------------------------
fig_sources = ica.plot_sources(raw_resmpl_all, title="ICA", show=True)
report.add_figure(
    fig=fig_sources,
    title="ICA source traces",
    caption=(
        "Component order is arbitrary. Interpret ocular/cardiac components "
        "using channel-based scores (see EOG-guided diagnostics), not index."
    ),
    section="ICA decomposition",
    tags=("ica", "qc"),
)
savefig_logged(fig_sources, "03_C_ICA_sources.png", dpi=150, bbox_inches="tight")
plt.close(fig_sources)

fig_components = ica.plot_components(show=True)
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
        f"03_C_ICA_components_{idx:02d}.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig_comp)

# %%
# -------------------------------------------------------------------------
# EOG-guided diagnostics block: score components against EOG channels
# -------------------------------------------------------------------------
eog_summary_text, eog_ranked_components, eog_auto_candidates = summarize_eog_component_scores(
    ica,
    ann_files,
)
log_to_console_and_report(
    eog_summary_text,
    title="EOG-guided component diagnostics",
    section="ICA decomposition",
)

# %%
# -------------------------------------------------------------------------
# Exclusion block: manual subject-specific list
# -------------------------------------------------------------------------
has_subject_exclude = subj_str in manual_exclude_by_subject
manual_exclude_list = list(manual_exclude_by_subject.get(subj_str, manual_exclude_default))
exclude_list = list(manual_exclude_list)

if not has_subject_exclude:
    report.save(out_report, overwrite=True, open_browser=False)
    raise SystemExit(
        "No manual ICA components are configured for rejection.\n"
        f"Choose components to reject based on report:\n  {out_report}\n"
        "Then set manual_exclude_by_subject for this subject, set "
        "`reuse_prefit_ica = True`, and re-run 03_C_ICA.py."
    )

valid_exclude = [idx for idx in exclude_list if idx < ica.n_components_]
dropped_exclude = [idx for idx in exclude_list if idx >= ica.n_components_]
merged_exclude = valid_exclude
exclude_strategy = "Manual exclusion list only."

ica.exclude = merged_exclude
ica.save(ica_postreject_file, overwrite=True)

log_to_console_and_report(
    f"Exclude strategy: {exclude_strategy}\n"
    f"Configured/manual exclude list: {exclude_list}\n"
    f"ICA exclude set (valid indices): {merged_exclude}\n"
    f"EOG diagnostic candidates (not auto-applied): {eog_auto_candidates}\n"
    f"Dropped indices not available in this fit: {dropped_exclude}\n"
    f"Saved post-rejection ICA model: {ica_postreject_file}",
    title="ICA exclude list",
    section="ICA decomposition",
)

# %%
# -------------------------------------------------------------------------
# FLUX projection/apply block: apply ICA to each annotated file and save
# -------------------------------------------------------------------------
processed_pairs = []
for ann_path, ica_path in zip(ann_files, ica_files):
    raw_ica = mne.io.read_raw_fif(ann_path, preload=True, verbose=True)
    input_bads = list(dict.fromkeys(str(ch) for ch in raw_ica.info.get("bads", [])))

    # Sub-01 exception: replace missing/unreliable EOG blink labels with ICA000 labels.
    ica0_blink_count = 0
    if use_ica0_blink_for_subject:
        raw_ica_for_blink = raw_ica.copy().pick("meg")
        raw_ica_for_blink.resample(200)
        raw_ica_for_blink.filter(1, 40)
        ica0_blink_annotations, ica0_blink_count = build_ica0_blink_annotations(
            ica,
            raw_ica_for_blink,
        )
        raw_ica.set_annotations(raw_ica.annotations + ica0_blink_annotations)

    ica.apply(raw_ica)
    # Preserve upstream bad-channel metadata after ICA projection.
    raw_ica.info["bads"] = list(input_bads)
    raw_ica.save(ica_path, overwrite=True)
    processed_pairs.append((ann_path, ica_path, ica0_blink_count, input_bads))

log_to_console_and_report(
    "Applied ICA and saved files:\n"
    + "\n".join(
        f"  - {src} -> {dst}"
        + (f" (ICA000 blinks added: {n_blinks})" if use_ica0_blink_for_subject else "")
        + f" (propagated_bads={bads})"
        for src, dst, n_blinks, bads in processed_pairs
    ),
    title="ICA-applied outputs",
    section="Apply ICA",
)

# %%
# -------------------------------------------------------------------------
# FLUX before/after comparison block: visualize artefact attenuation
# -------------------------------------------------------------------------
# FLUX reference channel set for frontal inspection (MEG mags).
compare_chs = [];#["MEG0311", "MEG0121", "MEG1211", "MEG1411"]
compare_ann = ann_files[0]
compare_ica = ica_files[0]

raw_before = mne.io.read_raw_fif(compare_ann, preload=True, verbose=False)
raw_after = mne.io.read_raw_fif(compare_ica, preload=True, verbose=False)

available_compare_chs = [ch for ch in compare_chs if ch in raw_before.ch_names]
missing_compare_chs = [ch for ch in compare_chs if ch not in raw_before.ch_names]
if not available_compare_chs:
    available_compare_chs = [ch for ch in raw_before.ch_names if ch.startswith("MEG")][:4]
    missing_compare_chs = compare_chs
chan_idxs = [raw_before.ch_names.index(ch) for ch in available_compare_chs]

if missing_compare_chs:
    print(f"WARNING: missing FLUX comparison channels: {missing_compare_chs}")

fig_before = raw_before.plot(order=chan_idxs, duration=5, show=True)
report.add_figure(
    fig=fig_before,
    title=f"Before ICA ({op.basename(compare_ann)})",
    caption=f"Channels: {available_compare_chs}",
    section="Before/after ICA",
    tags=("raw", "qc"),
)
savefig_logged(fig_before, "03_C_ICA_before_projection.png", dpi=150, bbox_inches="tight")
plt.close(fig_before)

fig_after = raw_after.plot(order=chan_idxs, duration=5, show=True)
report.add_figure(
    fig=fig_after,
    title=f"After ICA ({op.basename(compare_ica)})",
    caption=f"Channels: {available_compare_chs}",
    section="Before/after ICA",
    tags=("raw", "qc"),
)
savefig_logged(fig_after, "03_C_ICA_after_projection.png", dpi=150, bbox_inches="tight")
plt.close(fig_after)

log_to_console_and_report(
    f"Before/after ICA comparison used file:\n  - before: {compare_ann}\n  - after: {compare_ica}\n"
    f"Channels: {available_compare_chs}",
    title="Before/after comparison",
    section="Before/after ICA",
)



# %%
# -------------------------------------------------------------------------
# Optional: visualize EOG scores with a bar plot for top-ranked components
# -------------------------------------------------------------------------
if eog_ranked_components:
    top_comps, top_scores, _n_vals = zip(*eog_ranked_components[:10])
    fig_eog_scores, ax = plt.subplots(figsize=(8, 4))
    ax.bar(
        [f"ICA{comp_idx:03d}" for comp_idx in top_comps],
        top_scores,
        color="skyblue",
        edgecolor="black",
    )
    ax.axhline(eog_guided_score_threshold, color="red", linestyle="--", label="Auto-candidate threshold")
    ax.set_ylabel("Mean absolute EOG score")
    ax.set_title("Top-ranked ICA components by EOG score")
    ax.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()
    report.add_figure(
        fig=fig_eog_scores,
        title="EOG score ranking of ICA components",
        caption=(
            "Bar plot of mean absolute EOG scores for top-ranked components. "
            f"Red dashed line indicates diagnostic threshold ({eog_guided_score_threshold:.2f})."
        ),
        section="EOG diagnostics",
        tags=("ica", "qc", "eog"),
    )
    savefig_logged(fig_eog_scores, "03_C_ICA_eog_scores.png", dpi=150, bbox_inches="tight")
    plt.close(fig_eog_scores)
else:
    log_to_console_and_report(
        "No valid EOG scores available to plot.",
        title="EOG score ranking",
        section="EOG diagnostics",
    )


"""
# %%
# -------------------------------------------------------------------------
# Optional: Scrollable plot for seeing why we can't rely on EOG for annotation (hopefully only for subj 01)
# -------------------------------------------------------------------------

import mne
from mne.preprocessing import read_ica

# Paths
ica_path = "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/ica_models/sub-01_ses-01_ica_prefit.fif"
raw_path = "/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/preprocessing/sub-01/ses-01/sss/SUB01_block01_ann.fif"

# Window to inspect
start_sec = 97.0
duration_sec = 5.0

# Load data
raw = mne.io.read_raw_fif(raw_path, preload=True)
ica = read_ica(ica_path)

# Recreate ICA-fit preprocessing on MEG
raw_fit = raw.copy().pick("meg")
raw_fit.resample(200)
raw_fit.filter(1.0, 40.0)

# ICA sources
src = ica.get_sources(raw_fit)

# Keep first and fifth components (ICA000 and ICA004)
src_2 = src.copy().pick([src.ch_names[0], src.ch_names[4]])

# EOG channels
eog_chs = [ch for ch in ("EOG001", "EOG002") if ch in raw.ch_names]
if not eog_chs:
    raise RuntimeError("No EOG channels found in raw file.")

eog = raw.copy().pick(eog_chs)
eog.resample(200)
eog.filter(1.0, 40.0, picks="all")

# Combine ICA components + EOG in one scrollable browser
combo = src_2.copy().add_channels([eog], force_update_info=True)

combo.plot(
    start=start_sec,
    duration=duration_sec,
    n_channels=len(combo.ch_names),
    scalings="auto",
    title="ICA000 + ICA004 + EOG001/EOG002",
    block=True,
)
"""


# %%
# -------------------------------------------------------------------------
# Final report save
# -------------------------------------------------------------------------
report.save(out_report, overwrite=True, open_browser=False)
print(f"Saved HTML report to: {out_report}")

# %%
