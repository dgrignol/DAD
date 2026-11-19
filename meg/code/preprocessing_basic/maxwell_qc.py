#!/usr/bin/env python3
"""Plot per-run Maxwell channel evidence without repeating preprocessing.

Usage
-----
Refresh P04's existing 4A plots and embedded report on the local mount:

    python meg/code/preprocessing_basic/maxwell_qc.py \
      --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
      --raw-root /Volumes/MORWUR/Projects/DAMIANO/Attention/data \
      --subjects 04 --window-sec 20 --seed 97 --overwrite

Use --subjects {04..39} for the cohort. On the cluster, substitute its mount
roots. This command submits no jobs and never runs detection, SSS or filtering.

Inputs / outputs
----------------
Existing 4A JSON + per-run score NPZ + exact raw/SSS sources -> raw per-channel
PSD with sensor locations -> noisy/flat heatmaps -> native 20-second MEG traces
-> matched per-channel raw/SSS PSDs excluding applied bads -> plot metadata and
a rebuilt self-contained primary HTML report. Original
4A JSON, decisions, scores, source manifests and FIFs are never overwritten.
03_A_MaxFilt.py also calls these plotting functions for future processing runs.

Scientific / display assumptions
-------------------------------
Labels come from the recorded per-run decision; threshold crossings never
silently reclassify a sensor. Red [D] means detected in this run (including
header bads); amber [A] means applied only through harmonization/manual choices.
NaN scores are unscored, not zero. Noisy scores exceed their own channel/window
limits; flat scores fall below theirs. Columns represent saved evaluated windows,
including any gaps, rather than an invented continuous time axis.

Trace windows are selected uniformly on the native sample grid with a stable
seed derived from the source filename and --seed. All MEG sensors are shown,
separating mag/grad and retaining physical units. Only a constant median per
trace is removed for display; there is no filtering, resampling or per-channel
amplitude normalization. A random window may miss an intermittent fault; use
the full-recording heatmaps to choose additional inspections when needed.

PSDs use the same first up-to-120 seconds, Welch windows of up to 2000 samples
with 50% overlap, and 1-60 Hz. Each curve is a separate sensor, colored by its
location using MNE's public Spectrum.plot API. First-plot bads are also named
and overlaid red/amber with dashed lines. The final raw/SSS panels exclude the
same applied bad-channel list and share y limits within each sensor type. PSDs
are power per Hz in dB, not amplitudes. These views cannot resolve the 0.01 Hz
transition or inspect 100 Hz line noise. No diagnostic filter is applied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import mne
import numpy as np

from continuous_cleaning import file_hash, file_stamp
from pipeline_config import normalize_subject_id
from preprocessing_stage_utils import guarded_path, write_json


DETECTED_COLOR = "#b2182b"
APPLIED_COLOR = "#b77900"


# -----------------------------------------------------------------------------
# Score arrays and deterministic native-sample window selection
# -----------------------------------------------------------------------------

def load_scores(path):
    """Read numeric/string NPZ arrays without allowing pickled Python objects."""
    with np.load(path, allow_pickle=False) as saved:
        return {key: saved[key] for key in saved.files}


def score_values(scores, channel_type, criterion):
    """Return channel names, scores, limits and exact per-cell exceedance mask.

    MNE stores per-channel limits as (channels, 1); broadcasting also supports
    per-window limits. Missing/header-excluded values remain masked in plots.
    """
    names = np.asarray(scores["ch_names"]).astype(str)
    subset = np.asarray(scores["ch_types"]) == channel_type
    values = np.asarray(scores[f"scores_{criterion}"], dtype=float)[subset]
    limits = np.broadcast_to(np.asarray(scores[f"limits_{criterion}"], dtype=float),
                             np.asarray(scores[f"scores_{criterion}"]).shape)[subset]
    valid = np.isfinite(values) & np.isfinite(limits)
    crossed = valid & ((values > limits) if criterion == "noisy" else (values < limits))
    return names[subset], values, limits, crossed


def choose_window(raw, source_filename, duration_sec=20., seed=97):
    """Choose a repeatable, uniformly random half-open native sample interval."""
    if not np.isfinite(duration_sec) or duration_sec <= 0:
        raise ValueError("Trace duration must be positive and finite.")
    count = min(raw.n_times, max(1, int(round(duration_sec * raw.info["sfreq"]))))
    stable_seed = int.from_bytes(hashlib.sha256(f"{seed}:{source_filename}".encode()).digest()[:8], "big")
    start = int(np.random.default_rng(stable_seed).integers(raw.n_times - count + 1))
    return {"start_index": start, "stop_index": start + count,
            "start_sec": start / raw.info["sfreq"], "stop_sec": (start + count) / raw.info["sfreq"],
            "first_samp": int(raw.first_samp), "sfreq_hz": float(raw.info["sfreq"]),
            "native_start_sample": int(raw.first_samp + start),
            "native_stop_sample_exclusive": int(raw.first_samp + start + count), "seed": int(seed)}


# -----------------------------------------------------------------------------
# Per-channel PSDs: raw overview first, matched retained-channel SSS view last
# -----------------------------------------------------------------------------

def save_psds(*, raw, sss, source_filename, detected_bads, applied_bads,
              prefix, raw_root, show=False):
    """Save two PSD figures using MNE's individual-channel/sensor-inset layout.

    Inputs are native raw/SSS readers for the same acquisition. Only the first
    up-to-120 seconds are read. One raw spectrum supplies both the all-channel
    overview and the retained-channel comparison, avoiding duplicate estimation.
    Returns figure metadata and the exact sampling/channel choices for review.
    Neither input's samples, bad list or annotations are modified.
    """
    for label, source in (("raw", raw), ("SSS", sss)):
        if source.n_times < 4:
            raise ValueError(f"{label} recording is too short for a PSD.")
    if (raw.n_times, raw.first_samp, raw.info["sfreq"], raw.ch_names) != (
            sss.n_times, sss.first_samp, sss.info["sfreq"], sss.ch_names):
        raise ValueError("Raw/SSS PSD inputs have different native coordinates/channels.")
    stop = min(raw.n_times, int(round(120 * raw.info["sfreq"])))
    n_fft = min(2000, stop)
    settings = dict(method="welch", fmin=1., fmax=min(60., raw.info["sfreq"] / 2),
                    tmin=0., tmax=(stop - 1) / raw.info["sfreq"],
                    n_fft=n_fft, n_overlap=n_fft // 2, exclude=(), verbose="ERROR")
    all_names = [raw.ch_names[i] for i in mne.pick_types(raw.info, meg=True, exclude=[])]
    good_names = [name for name in all_names if name not in applied_bads]
    if not good_names:
        raise ValueError("No retained MEG channels remain for the raw/SSS comparison.")
    raw_spectrum = raw.compute_psd(picks=all_names, **settings)
    sss_spectrum = sss.compute_psd(picks=good_names, **settings)
    # Display copies use spatial colors for every location. Recorded decisions
    # are highlighted explicitly below; no Raw.info bad-channel list is edited.
    raw_spectrum.info["bads"] = []
    sss_spectrum.info["bads"] = []
    channel_types = [kind for kind in ("mag", "grad") if kind in raw.get_channel_types()]
    overview, raw_axes = plt.subplots(len(channel_types), 1, figsize=(15, 5 * len(channel_types)), squeeze=False)
    comparison, pair_axes = plt.subplots(len(channel_types), 2,
        figsize=(16, 5 * len(channel_types)), squeeze=False)
    for row, kind in enumerate(channel_types):
        names = [raw.ch_names[i] for i in mne.pick_types(raw.info, meg=kind, exclude=[])]
        good = [name for name in names if name in good_names]
        axis = raw_axes[row, 0]
        raw_spectrum.plot(picks=names, average=False, spatial_colors=True,
                          exclude=(), alpha=.4, axes=axis, show=False)
        # MNE handles spatial colors and the sensor inset. Overlay explicitly
        # selected channel PSDs in physical power units to make the bads clear.
        values = raw_spectrum.get_data(picks=names)
        scale = 1e15 if kind == "mag" else 1e13  # fT or fT/cm before squaring.
        for name, value in zip(names, values):
            if name in set(detected_bads) | set(applied_bads):
                detected = name in detected_bads
                power_db = 10 * np.log10(np.maximum(value * scale ** 2, np.finfo(float).tiny))
                axis.plot(raw_spectrum.freqs, power_db,
                    color=DETECTED_COLOR if detected else APPLIED_COLOR,
                    linestyle="--" if detected else "-.", linewidth=1.3,
                    label=name + (" [D]" if detected else " [A]"), zorder=4)
        if axis.get_legend_handles_labels()[0]:
            axis.legend(loc="upper left", fontsize=7, ncol=3)
        axis.set_title(f"Raw {kind}: all {len(names)} individual channels")
        axis.set_xlabel("Frequency (Hz)")
        for column, (spectrum, title) in enumerate(((raw_spectrum, "Before SSS"), (sss_spectrum, "After SSS"))):
            target = pair_axes[row, column]
            if good:
                spectrum.plot(picks=good, average=False, spatial_colors=True,
                              exclude=(), alpha=.5, axes=target, show=False)
            else:
                target.text(.5, .5, "No retained channels of this type", transform=target.transAxes, ha="center")
            target.set_title(f"{title}: {len(good)} individual {kind} channels")
            target.set_xlabel("Frequency (Hz)")
        # Comparable axes prevent autoscaling from disguising the size of a
        # spectral change. Sensor-location colors are stable across both sides.
        limits = [axis.get_ylim() for axis in pair_axes[row]]
        for axis in pair_axes[row]:
            axis.set_ylim(min(limit[0] for limit in limits), max(limit[1] for limit in limits))
    duration = stop / raw.info["sfreq"]
    overview.suptitle(f"{source_filename} | raw PSD, first {duration:.3f} s\n"
        "Each curve is a sensor; inset shows location colors. Dashed red [D]: detected here; dash-dot amber [A]: applied only.", fontsize=11)
    comparison.suptitle(f"{source_filename} | matched raw/SSS PSD, first {duration:.3f} s\n"
        f"Same retained channels and Welch settings; applied bads excluded on both sides: {', '.join(applied_bads) or '(none)'}. "
        "Matched scales within mag/grad.", fontsize=10, wrap=True)
    plots = []
    for figure, suffix, title, section in (
        (overview, "raw_channels_psd", "Raw PSD: every channel and sensor locations", "Raw PSD and sensor locations"),
        (comparison, "before_after_psd", "Retained-channel PSD: before and after SSS", "Before and after Maxwell filtering"),
    ):
        # MNE's location insets are positioned relative to their parent axes;
        # explicit margins avoid tight_layout warnings and clipped axis labels.
        figure.subplots_adjust(left=.07, right=.98, bottom=.07, top=.88,
                               hspace=.35, wspace=.22)
        path = guarded_path(f"{prefix}_{suffix}.png", raw_root)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=110)
        if show:
            plt.show(block=True)
        plt.close(figure)
        plots.append({"filename": path.name, "title": title, "section": section})
    return plots, {"start_index": 0, "stop_index": stop, "duration_sec": duration,
        "native_start_sample": int(raw.first_samp), "native_stop_sample_exclusive": int(raw.first_samp + stop),
        "fmin_hz": settings["fmin"], "fmax_hz": settings["fmax"], "n_fft": n_fft,
        "n_overlap": n_fft // 2, "raw_channels": all_names, "comparison_channels": good_names}


# -----------------------------------------------------------------------------
# Heatmaps: all scores and threshold crossings, with explicit selected rows
# -----------------------------------------------------------------------------

def save_heatmaps(*, scores, source_filename, detected_bads, applied_bads,
                  prefix, raw_root, show=False):
    """Save MNE-style two-panel noisy/flat heatmaps for each available MEG type.

    The right panel masks non-crossings using each cell's actual limit, rather
    than merely clipping colors at a global threshold. Red/amber row outlines
    and labels identify recorded channel decisions on both panels. Sparse
    background labels keep dense arrays readable; every selected row is named.
    """
    plots = []
    detected = set(detected_bads)
    applied_only = set(applied_bads or []) - detected
    bins = np.asarray(scores["bins"])
    if bins.ndim != 2 or bins.shape[1] != 2 or not len(bins):
        raise ValueError("Maxwell score bins must contain start/stop pairs.")
    for channel_type in ("mag", "grad"):
        for criterion in ("noisy", "flat"):
            names, values, limits, crossed = score_values(scores, channel_type, criterion)
            if not len(names):
                continue
            if values.shape[1] != len(bins):
                raise ValueError("Score columns do not match the evaluated windows.")
            # Flat scores are amplitudes; convert both scores and limits into
            # physical display units. Noisy scores remain dimensionless.
            multiplier = (1e15 if channel_type == "mag" else 1e13) if criterion == "flat" else 1.
            values, limits = values * multiplier, limits * multiplier
            unit = ("fT" if channel_type == "mag" else "fT/cm") if criterion == "flat" else "z score"
            figure, axes = plt.subplots(1, 2, figsize=(15, max(7, len(names) * .100)),
                                        layout="constrained", sharey=True)
            selected = {i for i, name in enumerate(names) if name in detected | applied_only}
            # Omit nearby background labels so selected names cannot collide
            # with them. Every selected row still has its own label/outline.
            background_ticks = {i for i in range(0, len(names), max(1, len(names) // 45))
                                if all(abs(i - chosen) > 2 for chosen in selected)}
            ticks = sorted(background_ticks | selected)
            for axis, thresholded in zip(axes, (False, True)):
                display = np.ma.masked_invalid(values)
                if thresholded:
                    display = np.ma.masked_where(~crossed, display)
                cmap = plt.get_cmap("Reds" if criterion == "noisy" else "Reds_r").copy()
                cmap.set_bad("#e4e4e4" if not thresholded else "white")
                kwargs = {}
                if thresholded:
                    finite_limits = limits[np.isfinite(limits)]
                    if finite_limits.size:
                        kwargs = {"vmin": float(np.min(finite_limits))} if criterion == "noisy" else {"vmin": 0., "vmax": float(np.max(finite_limits))}
                if display.count() == 0:
                    kwargs.update(vmin=0., vmax=1.)
                image = axis.imshow(display, aspect="auto", interpolation="nearest", cmap=cmap, **kwargs)
                figure.colorbar(image, ax=axis, label=unit, shrink=.75)
                axis.set_title("All scores" if not thresholded else ("Scores > own limit" if criterion == "noisy" else "Scores < own limit"))
                x_ticks = np.unique(np.linspace(0, len(bins) - 1, min(12, len(bins))).astype(int))
                axis.set_xticks(x_ticks, [f"{bins[i, 0]:.1f}-{bins[i, 1]:.1f}" for i in x_ticks], rotation=55, ha="right", fontsize=8)
                axis.set_xlabel("Evaluated window (seconds from recording start)")
                axis.set_yticks(ticks, [names[i] + (" [D]" if names[i] in detected else " [A]" if names[i] in applied_only else "") for i in ticks], fontsize=7)
                axis.tick_params(axis="y", labelleft=True)
                for tick, i in zip(axis.get_yticklabels(), ticks):
                    if i in selected:
                        color = DETECTED_COLOR if names[i] in detected else APPLIED_COLOR
                        tick.set_color(color); tick.set_fontweight("bold")
                        axis.add_patch(Rectangle((-.5, i - .5), len(bins), 1, fill=False, edgecolor=color, linewidth=.8))
                if thresholded and not crossed.any():
                    axis.text(.5, .5, "No finite scores crossed their limits", transform=axis.transAxes,
                              ha="center", bbox={"facecolor": "white", "alpha": .9})
            axes[0].set_ylabel("MEG channel (original channel order)")
            figure.suptitle(f"{source_filename} | {criterion} {channel_type}\n"
                "Red [D]: detected in this run; amber [A]: applied only by policy/manual choice.\n"
                "Blank right cells did not cross; gray left cells are unscored. Crossings are not new decisions.", fontsize=11)
            path = guarded_path(f"{prefix}_{criterion}_{channel_type}_heatmap.png", raw_root)
            path.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(path, dpi=110)
            if show:
                plt.show(block=True)
            plt.close(figure)
            plots.append({"filename": path.name, "title": f"{criterion.capitalize()} {channel_type}: scores and limit crossings", "section": "Bad-channel score heatmaps"})
    return plots


# -----------------------------------------------------------------------------
# Raw traces: a bounded reproducible random window, never a processed derivative
# -----------------------------------------------------------------------------

def save_traces(*, raw, source_filename, detected_bads, applied_bads, prefix,
                raw_root, duration_sec=20., seed=97, show=False):
    """Overlay all raw MEG channels and identify each detected channel in a legend.

    Red-family colors distinguish detected traces individually. Amber traces
    were applied but not detected here; other channels are gray context. All
    traces within a sensor type share a physical y-axis scale, without clipping.
    Data are read only for the requested window and are not written back.
    """
    window = choose_window(raw, source_filename, duration_sec, seed)
    figure, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=True, layout="constrained")
    detected = set(detected_bads)
    applied_only = set(applied_bads or []) - detected
    for axis, channel_type in zip(axes, ("mag", "grad")):
        picks = mne.pick_types(raw.info, meg=channel_type, eeg=False, exclude=[])
        if not len(picks):
            axis.text(.5, .5, f"No {channel_type} channels", transform=axis.transAxes, ha="center")
            continue
        data = raw.get_data(picks=picks, start=window["start_index"], stop=window["stop_index"])
        data = (data - np.median(data, axis=1, keepdims=True)) * (1e15 if channel_type == "mag" else 1e13)
        times = (window["start_index"] + np.arange(data.shape[1])) / raw.info["sfreq"]
        names = [raw.ch_names[i] for i in picks]
        for i, name in enumerate(names):
            if name not in detected:
                axis.plot(times, data[i], color=APPLIED_COLOR if name in applied_only else "0.55",
                          alpha=.65 if name in applied_only else .20, linewidth=.55)
        selected = [(i, name) for i, name in enumerate(names) if name in detected]
        colors = plt.get_cmap("Reds")(np.linspace(.55, .95, max(1, len(selected))))
        for (i, name), color in zip(selected, colors):
            retained = " (not applied)" if applied_bads is not None and name not in applied_bads else ""
            axis.plot(times, data[i], color=color, linewidth=.85, label=name + retained, zorder=3)
        handles = [Line2D([], [], color="0.55", label="Other raw MEG channels")]
        if applied_only.intersection(names):
            handles.append(Line2D([], [], color=APPLIED_COLOR, label="Applied only: " + ", ".join(sorted(applied_only.intersection(names)))))
        lines, labels = axis.get_legend_handles_labels()
        axis.legend(handles=handles + lines, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
        axis.set_ylabel("Magnetometers (fT)" if channel_type == "mag" else "Gradiometers (fT/cm)")
        axis.set_title(f"{channel_type}: {len(selected)} detected in this run, {len(picks)} channels shown")
        axis.grid(alpha=.2)
    axes[-1].set_xlabel("Seconds from recording start (native samples)")
    figure.suptitle(f"{source_filename} | random raw window {window['start_sec']:.3f}-{window['stop_sec']:.3f} s\n"
        f"Native samples [{window['native_start_sample']}, {window['native_stop_sample_exclusive']}) | seed {seed}\n"
        "Warm-colored named traces: detected here. Constant medians removed for display only; no filtering or amplitude normalization.", fontsize=11)
    path = guarded_path(f"{prefix}_bad_channel_traces.png", raw_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=110)
    if show:
        plt.show(block=True)
    plt.close(figure)
    return {"filename": path.name, "title": f"Raw channels: random {window['stop_sec'] - window['start_sec']:.1f} s window", "section": "Bad-channel raw traces"}, window


def build_run_qc(*, raw, source_path, sss_path, scores, scores_path, detected_bads, applied_bads,
                 prefix, raw_root, duration_sec=20., seed=97, show=False):
    """Generate one run's images and portable metadata without editing its data."""
    source_path = Path(source_path)
    absent = (set(detected_bads) | set(applied_bads)) - set(raw.ch_names)
    if absent:
        raise ValueError(f"QC decision channels absent from source: {sorted(absent)}")
    if not set(np.asarray(scores["ch_names"]).astype(str)).issubset(raw.ch_names):
        raise ValueError("Scored channels do not match this recording.")
    sss = mne.io.read_raw_fif(sss_path, preload=False, verbose="ERROR")
    psd_plots, psd_settings = save_psds(raw=raw, sss=sss, source_filename=source_path.name,
        detected_bads=detected_bads, applied_bads=applied_bads, prefix=prefix, raw_root=raw_root, show=show)
    plots = save_heatmaps(scores=scores, source_filename=source_path.name,
        detected_bads=detected_bads, applied_bads=applied_bads, prefix=prefix, raw_root=raw_root, show=show)
    trace, window = save_traces(raw=raw, source_filename=source_path.name,
        detected_bads=detected_bads, applied_bads=applied_bads, prefix=prefix, raw_root=raw_root,
        duration_sec=duration_sec, seed=seed, show=show)
    return {"source_filename": source_path.name, "source_stamp": file_stamp(source_path),
            "score_sha256": file_hash(scores_path), "detected_bad_channels": list(detected_bads),
            "applied_bad_channels": list(applied_bads), "window": window, "psd": psd_settings,
            "plots": psd_plots[:1] + plots + [trace] + psd_plots[1:]}


# -----------------------------------------------------------------------------
# Refresh saved evidence: verify exact sources, write only plots/plot metadata/HTML
# -----------------------------------------------------------------------------

def refresh_subject(*, derivatives_root, raw_root, subject_id, processing_pass="review",
                    duration_sec=20., seed=97, overwrite=False):
    """Refresh an existing participant report; fail on changed raw source stamps."""
    subject = "sub-" + normalize_subject_id(subject_id)
    report_dir = Path(derivatives_root) / "reports" / subject / "ses-01"
    figure_dir = Path(derivatives_root) / "preprocessing" / subject / "ses-01/figures"
    stem = f"03_A_MaxFilt_{processing_pass}"
    report_path = report_dir / f"{stem}_report.json"
    payload = json.loads(report_path.read_text())
    if payload["subject_bids_id"] != subject or payload["processing_pass"] != processing_pass:
        raise ValueError("Report subject/pass does not match the requested QC refresh.")
    metadata_path = guarded_path(report_dir / f"{stem}_plot_metadata.json", raw_root)
    if metadata_path.exists() and not overwrite:
        raise FileExistsError(f"Plot metadata exists; pass --overwrite: {metadata_path}")
    sources = {row["stamp"]["name"]: row for row in payload["source_files"]}
    outputs = {row["stamp"]["name"]: row for row in payload["evidence_files"]}
    runs = []
    for run in payload["runs"]:
        source = sources[run["source_filename"]]
        path = Path(raw_root) / source["relative_path"]
        if file_stamp(path) != source["stamp"]:
            raise RuntimeError(f"Raw source changed since 4A: {path}")
        # Resolve the SSS acquisition through portable saved lineage, never an
        # absolute path from another machine's /mnt or /Volumes mount point.
        output = outputs[Path(run["output_path"]).name]
        sss_path = Path(derivatives_root) / output["relative_path"]
        if file_stamp(sss_path) != output["stamp"]:
            raise RuntimeError(f"SSS source changed since 4A: {sss_path}")
        prefix = figure_dir / f"{stem}_{path.stem}"
        scores_path = report_dir / f"{stem}_{path.stem}_bad_channel_scores.npz"
        raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
        runs.append(build_run_qc(raw=raw, source_path=path, sss_path=sss_path, scores=load_scores(scores_path),
            scores_path=scores_path, detected_bads=run["detected_bad_channels"],
            applied_bads=run["applied_bad_channels"], prefix=prefix, raw_root=raw_root,
            duration_sec=duration_sec, seed=seed))
        print(f"Plotted {subject}: {path.name}", flush=True)
    write_json(metadata_path, {"stage_report_sha256": file_hash(report_path), "runs": runs}, raw_root=raw_root)
    from render_stage03_mne_reports import render_stage_report
    return render_stage_report(derivatives_root=Path(derivatives_root), raw_root=Path(raw_root),
        subject_id=subject_id, processing_pass=processing_pass, stage_name="maxwell_filter",
        overwrite=True, primary_report=True)


def main(argv=None):
    """Refresh requested participants from existing scores and bounded raw windows."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--derivatives-root", required=True, type=Path)
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--subjects", nargs="+", required=True)
    parser.add_argument("--processing-pass", choices=("review", "apply"), default="review")
    parser.add_argument("--window-sec", type=float, default=20.)
    parser.add_argument("--seed", type=int, default=97)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    for subject in args.subjects:
        print(refresh_subject(derivatives_root=args.derivatives_root, raw_root=args.raw_root,
            subject_id=subject, processing_pass=args.processing_pass, duration_sec=args.window_sec,
            seed=args.seed, overwrite=args.overwrite), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
