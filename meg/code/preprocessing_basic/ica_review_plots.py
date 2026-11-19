#!/usr/bin/env python3
"""Shared display helpers for final ICA review; never fit or persist ICA decisions.

Usage inside 03_C_ICA.py or review_final_ica.py
-------------------------------------------
    display = continuous_display_copy(fit_raw)
    figures = model.plot_properties(display, picks=[0], show=False)
    describe_properties(figures[0], display.info['sfreq'])
    path, selection = native_context(model, native, record, figure_dir)

MNE 1.10.x get_sources() collapses concatenated Raw bookkeeping in an order
that can report only the first piece's duration although its source data array
contains all samples. Rebuild a SINGLE public RawArray for DISPLAY, retaining
all data, coordinates, annotations and joins. This changes neither fit inputs,
model weights, protected preparation code, native derivatives nor exclusions.

Properties use arbitrary two-second segments, NOT stimulus epochs. The helper
restores shared-axis ticks and replaces the misleading ERP/ERF title. Native
context projects the final model onto one original unfiltered five-second
window; annotations, selection policy and display normalization are explicit.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import mne
import numpy as np


# -----------------------------------------------------------------------------
# Full-length display objects and honest segment/property labels
# -----------------------------------------------------------------------------

def attach_annotations(raw, annotations, *, globalize=False):
    """Attach acquisition-time annotations once, optionally projecting sensor BADs globally.

    IC activation channels have no one-to-one native sensor identity. Diagnostic
    IC epochs therefore use the workflow's whole-interval BAD mask, retaining
    labels/times but omitting sensor-specific channel names on the disposable copy.
    Native objects and their original annotation channel assignments are untouched.
    """
    annotations = annotations.copy()
    if globalize:
        annotations = mne.Annotations(annotations.onset, annotations.duration, annotations.description,
            orig_time=annotations.orig_time, extras=annotations.extras)
    if annotations.orig_time is None:
        annotations.onset -= raw.first_time
    raw.set_annotations(annotations, emit_warning=False)
    return raw


def continuous_display_copy(raw):
    """Return one RawArray with every input sample and the exact original joins."""
    result = mne.io.RawArray(raw.get_data(), raw.info.copy(), first_samp=raw.first_samp,
                             verbose="ERROR")
    attach_annotations(result, raw.annotations)
    if result.n_times != raw.n_times or result.first_samp != raw.first_samp:
        raise RuntimeError("ICA display copy changed sample coordinates.")
    return result


def describe_properties(figure, sfreq):
    """Label non-task segments and return actual plotted sample/segment counts."""
    axes = {ax.get_label(): ax for ax in figure.axes}
    image, average = axes.get("image"), axes.get("erp")
    if image is None or average is None:
        raise RuntimeError("Unexpected MNE property layout; cannot verify the time axis.")
    stop = float(average.get_xlim()[1])
    ticks = [value for value in (0., .5, 1., 1.5) if value < stop] + [stop]
    # Setting an empty locator on the image also empties its SHARED mean axis.
    # Restore the common locator, then hide only the upper tick labels.
    average.set_xticks(ticks, [f"{value:g}" for value in ticks])
    average.tick_params(axis="x", labelbottom=True)
    image.tick_params(axis="x", labelbottom=False)
    image.set_title("2 s fit segments (not stimulus-locked)", fontsize=10)
    average.set_xlabel("Time within arbitrary segment (s)")
    average.set_ylabel("IC (AU)")
    variance = axes.get("variance")
    if variance is not None:
        variance.set_title("Variance of retained display segments", fontsize=10)
    segments = int(image.images[0].get_array().shape[0]) if image.images else None
    return {"kind": "fixed_length_fit_segments", "stimulus_locked": False,
        "segment_duration_sec": 2., "last_sample_sec": stop,
        "sfreq_hz": float(sfreq), "segments_plotted": segments,
        "mean_band": "95_percent_parametric_confidence_interval",
        "spectrum_band": "one_standard_deviation_across_segments"}


# -----------------------------------------------------------------------------
# Native examples: deterministic selection, labelled BADs and exact source time
# -----------------------------------------------------------------------------

def context_selection(raw, record, category=None):
    """Choose the first ocular interval, or a requested BAD label, without ranking ICs.

    Default starts two seconds before the first BAD_blink/BAD_eye_movement;
    absent ocular labels, use min(20 s, recording midpoint). Clamp the five-second
    window to the recording. Return metadata and an unpreloaded cropped Raw.
    """
    bad = [a for a in raw.annotations if a["description"].lower().startswith("bad")]
    selected = [a for a in bad if (a["description"].casefold() == category.casefold()
        if category else a["description"].lower().startswith(("bad_blink", "bad_eye_movement")))]
    if category and not selected:
        raise ValueError(f"No {category} annotation in this recording.")
    anchor = selected[0] if selected else None
    start = anchor["onset"] - raw.first_time - 2. if anchor else min(20., raw.times[-1] / 2)
    start = max(0., min(start, max(0., raw.times[-1] - 5.)))
    preview = raw.copy().crop(start, min(start + 5., raw.times[-1]))
    left = preview.first_samp / raw.info["sfreq"]
    right = (preview.first_samp + preview.n_times) / raw.info["sfreq"]
    spans = [{"description": str(a["description"]), "onset_acquisition_sec": float(a["onset"]),
        "duration_sec": float(a["duration"])} for a in bad
        if a["onset"] < right and a["onset"] + a["duration"] >= left]
    metadata = {"source_filename": record.source_filename,
        "selection_policy": "first_requested_BAD" if category else "first_ocular_else_fixed_fallback",
        "anchor_category": str(anchor["description"]) if anchor else None,
        "anchor_onset_acquisition_sec": float(anchor["onset"]) if anchor else None,
        "window_start_sample_raw": int(preview.first_samp),
        "window_stop_sample_raw_exclusive": int(preview.first_samp + preview.n_times),
        "window_start_recording_sec": float(preview.first_time - raw.first_time),
        "window_start_acquisition_sec": float(left), "window_stop_acquisition_sec": float(right),
        "bad_annotations": spans,
        "component_display": "mean_centered_and_unit_SD_independently_in_this_window",
        "signal_preparation": "native_unfiltered_annotation_input; no ICA removal",
        "interpretation": "qualitative temporal context, not BAD enrichment or component ranking"}
    return metadata, preview


def native_context(model, raw, record, figure_dir, *, category=None, output=None, show=False, center_frontal=True):
    """Plot native ICs, display-centered frontal mags and original MISC for one window.

    Subtract each frontal trace's window mean for display only: native SSS DC
    offsets can otherwise hide its temporal variation. No temporal filtering,
    model input change or native sample modification is performed.
    """
    metadata, preview = context_selection(raw, record, category)
    preview.load_data()
    sources = model.get_sources(continuous_display_copy(preview))
    if sources.n_times != preview.n_times:
        raise RuntimeError("Native ICA source projection was truncated.")
    data = sources.get_data()
    scale = np.std(data, axis=1, keepdims=True)
    scaled = (data - np.mean(data, axis=1, keepdims=True)) / np.maximum(scale, np.finfo(float).tiny)
    times = (np.arange(preview.n_times) + preview.first_samp) / preview.info["sfreq"]
    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True,
        gridspec_kw={"height_ratios": [3, 1, 1]})
    axes[0].imshow(scaled, origin="lower", aspect="auto", interpolation="nearest",
        cmap="RdBu_r", vmin=-4, vmax=4, extent=[times[0], times[-1], -.5, model.n_components_ - .5])
    axes[0].set_ylabel("Final IC index; within-window display z-score")
    axes[0].set_yticks(np.arange(0, model.n_components_, max(1, model.n_components_ // 15)))
    frontal = [ch for ch, kind in zip(preview.ch_names, preview.get_channel_types())
        if kind == "mag" and ch not in preview.info["bads"] and
        ch.startswith(("MEG011", "MEG012", "MEG031", "MEG041"))][:4]
    misc = [ch for ch, kind in zip(preview.ch_names, preview.get_channel_types()) if kind == "misc" and ch != "MISC008"]
    for ch in frontal:
        trace = preview.get_data(picks=[ch])[0]
        axes[1].plot(times, (trace - np.mean(trace) if center_frontal else trace) * 1e12, label=ch)
    for ch in misc:
        axes[2].plot(times, preview.get_data(picks=[ch])[0], label=ch)
    metadata["frontal_magnetometer_display"] = ("per_channel_window_mean_subtracted; pT; no temporal filtering"
        if center_frontal else "original_acquisition_levels_with_DC; pT; no temporal filtering")
    axes[1].set_ylabel("Window-centered frontal mag (pT)" if center_frontal else "Good frontal mag (pT)")
    axes[2].set_ylabel("MISC acquisition values")
    for ax, names in ((axes[1], frontal), (axes[2], misc)):
        if names:
            ax.legend(fontsize=7, loc="upper right")
        else:
            ax.text(.01, .5, "No channels available for this context", transform=ax.transAxes)
    labels = sorted({a["description"] for a in metadata["bad_annotations"]})
    colors = {label: plt.get_cmap("tab10")(i % 10) for i, label in enumerate(labels)}
    for annotation in metadata["bad_annotations"]:
        for ax in axes:
            ax.axvspan(annotation["onset_acquisition_sec"], annotation["onset_acquisition_sec"] + annotation["duration_sec"],
                       color=colors[annotation["description"]], alpha=.16)
    if labels:
        axes[0].legend(handles=[Patch(facecolor=colors[label], alpha=.35, label=label) for label in labels],
                       loc="upper right", fontsize=7)
    axes[-1].set_xlabel("Acquisition seconds (absolute native sample / sfreq)")
    anchor = metadata["anchor_category"] or "fixed fallback; no ocular label"
    fig.suptitle(record.source_filename + " | first example: " + anchor +
        "\nFinal-model native context; z-scores do not rank component strength", fontsize=11)
    fig.tight_layout()
    path = Path(output) if output else Path(figure_dir) / f"03_C_ICA_review_{record.source_stem}_native_context.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    if show:
        plt.show(block=True)
    plt.close(fig)
    return path, metadata
