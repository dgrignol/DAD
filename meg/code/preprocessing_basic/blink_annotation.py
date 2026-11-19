#!/usr/bin/env python3
"""Prepare frontal ICA and independent blink/eye-movement aids for Stage 4B.

Usage
-----
    # Headless preparation, normally through the existing cluster wrapper:
    python meg/code/preprocessing_basic/blink_annotation.py \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml --resume
    # Add --source sub04_run01_220526.fif to prepare just one exact acquisition.
    # Local review then uses manual_eye_annotations.py --preliminary-ica.
    # Inspect maps/sources, choose blink components, then eye components.
    # Each selection accepts 0 / ICA000 (first component), or explicit NONE.

Data flow: accepted SSS + current muscle/known severe annotations -> disposable
1-40 Hz/200 Hz fit copies -> ICA on ALL good left/right frontal mag/grad sensors
(at most ten components) -> model, native-grid sources and a small JSON sidecar.
The native source cache has NaN gaps for excluded/too-short stretches and filter
guards. No original MEG/MISC channel, recording coordinate or annotation is edited.
No final ICA, component removal, epoching or scientific acceptance occurs here.

Local detection uses selected components only, median/MAD normalization, a
detection-only low-pass at 10 Hz (sources are already 1-40 Hz), either polarity,
and SciPy peak finding. Approximate +/-0.25 s intervals are editable proposals,
not measured eyelid durations. BAD_blink_auto is used ONLY in the disposable
review object until explicit complete-source SAVE. Known muscle/severe BAD
intervals remain excluded; routine ocular examples are retained for this aid
only. Final Stage 4C still excludes every reviewed BAD category.

Independent eye proposals use detect_eye_movements on the SAME cached sources:
100 ms centered median smoothing, absolute median/MAD score >=4 for >=100 ms,
100 ms maximum within-stretch merge gap and 50 ms padding. These optional
artifact_annotation.eye_movement_auto parameters are review starting points.
Round the smoothing window to an odd native sample count, discard its half-
window support at every finite-stretch boundary, and never bridge NaN gaps.
Returned intervals index native samples; metadata uses absolute half-open
sample endpoints including first_samp. No resampling/time shift is introduced.
The 1 Hz cache high-pass attenuates long gaze plateaus, so sustained position
and full event duration cannot be inferred. Ordinary activity should occupy
most scored samples for robust scaling. BAD_eye_movement_auto, its settings,
scores and veto spans stay separate from blink proposals, even for the same
selected component. NONE disables only its group and never proves absence.
Both automatic BAD labels remain provisional until explicit whole-source SAVE.
Existing blink settings, digest inputs and all-frontal-v1 cache identities stay
unchanged; missing legacy eye selection requires no file/receipt migration.

References: MNE ICA and find_eog_events documentation; scipy.signal.find_peaks
and scipy.stats.median_abs_deviation; BLINKER (doi:10.3389/fnins.2017.00012).
The defaults below are DAD review starting points, not validated cohort choices.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mne
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import median_filter
from scipy.stats import median_abs_deviation

from continuous_cleaning import (annotation_rows, coordinates, digest, file_hash,
    file_stamp, filter_design, good_stretches, muscle_parameters, prepare_fit_copy,
    require_review, reviewed_muscle_override)
from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import (derivative_path_for_run, discover_raw_run_records,
    prepare_output_path, stage_output_dirs, utc_now, write_json)


AUTO_LABEL = "BAD_blink_auto"
EYE_AUTO_LABEL = "BAD_eye_movement_auto"
EYE_METHOD = "median-duration-v1"
DEFAULTS = {"n_components": 10, "threshold_std": 5., "min_peak_distance_sec": .25,
            "tmin_sec": -.25, "tmax_sec": .25, "long_peak_width_sec": .5}
EYE_DEFAULTS = {"median_window_sec": .1, "threshold_std": 4.,
                "min_duration_sec": .1, "merge_gap_sec": .1, "padding_sec": .05}


# -----------------------------------------------------------------------------
# Exact-source input, frontal sensor selection and compatible cache identity
# -----------------------------------------------------------------------------

def blink_settings(cfg):
    """Read the existing artifact section; never change the final ICA dimension."""
    return {**DEFAULTS, **cfg.get("artifact_annotation", {}).get("blink_auto", {})}


def eye_movement_settings(cfg):
    """Read optional eye aid settings without changing blink/ICA cache identity."""
    return {**EYE_DEFAULTS, **cfg.get("artifact_annotation", {}).get("eye_movement_auto", {})}


def frontal_channel_groups(raw):
    """Match both shipped MNE Vectorview frontal regions, including mag and grad."""
    return {side: [ch.replace(" ", "") for ch in mne.read_vectorview_selection(side)
                   if ch.replace(" ", "") in raw.ch_names]
            for side in ("Left-frontal", "Right-frontal")}


def aid_mask(raw):
    """Copy annotations, keeping ocular examples for the supplementary fit only."""
    result = raw.copy()
    remove = [i for i, label in enumerate(result.annotations.description)
              if label.lower().startswith(("bad_blink", "bad_eye_movement"))
              or label.lower() == "bad_filter_edge"]
    result.annotations.delete(remove)
    return result


def load_aid_input(cfg, record):
    """Load SSS headers plus existing, source-bound muscle/severe annotations.

    This does not require completed ocular review. A previous reviewed manual
    file contributes known severe intervals, but its ocular labels are removed
    by aid_mask. Missing or stale provisional muscle evidence requires a 4B
    refresh rather than silently fitting without the expected protection.
    A reviewed complete muscle set replaces provisional muscle labels, including
    reviewed zero; changed effective masks invalidate the usual cache identity.
    """
    dirs = stage_output_dirs(cfg)
    source = derivative_path_for_run(record, dirs["sss"], "_raw_sss.fif")
    raw = mne.io.read_raw_fif(source, preload=False, verbose="ERROR")
    report = json.loads((dirs["reports"] / "03_B_ArtefactAnnotation_review_report.json").read_text())
    rows = [r for r in report["runs"] if r["source_filename"] == record.source_filename]
    if len(rows) != 1 or rows[0].get("input_stamp") != file_stamp(source):
        raise RuntimeError("Missing/stale source-specific 4B muscle evidence; refresh that evidence first.")
    settings = cfg["artifact_annotation"]
    by_run = settings.get("muscle_threshold_by_run", {})
    run = record.acquisition_run_label
    threshold = next((float(by_run[key]) for key in (run, str(run), f"{run:02d}", f"run{run:02d}")
                      if key in by_run), float(settings["muscle_threshold"]))
    if rows[0]["muscle_threshold"] != threshold:
        raise RuntimeError("Muscle threshold changed; refresh the affected 4B evidence before aid preparation.")
    annotated = derivative_path_for_run(record, dirs["sss"], "_ann.fif")
    evidence = next((r for r in report.get("evidence_files", [])
                     if Path(r["relative_path"]).name == annotated.name), None)
    if not evidence or file_stamp(annotated) != evidence["stamp"]:
        raise RuntimeError("Provisional annotation evidence changed; refresh its report before preparation.")
    provisional = mne.io.read_raw_fif(annotated, preload=False, verbose="ERROR")
    used_review = rows[0].get("manual_eye_review", {})
    used_review = used_review if used_review.get("status") == "reviewed" else None
    # This is the existing 03_B artifact_signature contract. New blink settings
    # do not change muscle evidence, but changed muscle filters/settings do.
    signature = digest({"input": file_stamp(source), **muscle_parameters(settings, threshold),
        "eye_review": used_review})
    if f"stage4_input_sha256={signature}" not in str(provisional.info.get("description")):
        raise RuntimeError("Muscle processing settings changed; refresh affected 4B evidence before aid preparation.")
    ann = provisional.annotations
    if ann.orig_time != raw.annotations.orig_time:
        raise RuntimeError("Provisional annotations have a different acquisition origin.")
    # Append only BAD context; original source annotations already remain on raw.
    seen = {(float(a["onset"]), float(a["duration"]), a["description"]) for a in raw.annotations}
    manual = Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "manual_eye" / f"{record.source_stem}-annot.fif"
    receipt = settings.get("manual_eye_reviews", {}).get(record.source_filename, {})
    additions = [ann]
    if receipt.get("status") == "reviewed":
        if not manual.is_file() or receipt.get("sss_input") != file_stamp(source) or receipt.get("annotation_sha256") != file_hash(manual):
            raise RuntimeError("Previous manual review changed; reconcile it before reusing severe-period context.")
        reviewed = mne.read_annotations(manual)
        if reviewed_muscle_override(receipt, settings, threshold):
            # Replacement must precede deduplication; otherwise deleted muscle
            # periods would remain in the old automatic fit mask.
            raw.annotations.delete([i for i, label in enumerate(raw.annotations.description) if label.lower() == "bad_muscle"])
            additions = [ann[[i for i, label in enumerate(ann.description) if label.lower() != "bad_muscle"]]]
            seen = {(float(a["onset"]), float(a["duration"]), a["description"]) for a in raw.annotations}
        additions.append(reviewed)
    for annotations in additions:
        if annotations.orig_time != raw.annotations.orig_time:
            raise RuntimeError("Known annotations have a different acquisition origin.")
        for a in annotations:
            key = (float(a["onset"]), float(a["duration"]), a["description"])
            if a["description"].lower().startswith("bad") and key not in seen:
                raw.annotations.append(*key)
                seen.add(key)
    return aid_mask(raw), source


def cache_paths(cfg, record):
    """Keep this auxiliary model separate from the final Stage 4C models."""
    folder = Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "blink_aid"
    return {"model": folder / f"{record.source_stem}-ica.fif",
            "sources": folder / f"{record.source_stem}_sources_raw.fif",
            "metadata": folder / f"{record.source_stem}_blink_aid.json"}


def cache_identity(cfg, raw, source):
    """Bind fit inputs; detector thresholds deliberately do not invalidate ICA."""
    frontal = {ch for group in frontal_channel_groups(raw).values() for ch in group}
    kinds = dict(zip(raw.ch_names, raw.get_channel_types()))
    channels = [ch for ch in raw.ch_names if ch in frontal and kinds[ch] in ("mag", "grad")
                and ch not in raw.info["bads"]]
    if len(channels) < 2:
        raise RuntimeError("Too few good frontal MEG channels for the annotation aid.")
    fit = {key: value for key, value in cfg["ica"].items() if key.startswith("fit_")}
    return {"preparation": "all-frontal-v1", "sss_input": file_stamp(source),
            "coordinates": coordinates(raw), "mask": annotation_rows(aid_mask(raw)),
            "channels": channels, "fit_settings": fit,
            "n_components_max": int(blink_settings(cfg)["n_components"]),
            "method": "fastica", "random_state": 97, "max_iter": "auto"}


def read_cache(cfg, record, identity):
    """Validate the saved model, source grid and current input/mask before reuse."""
    paths = cache_paths(cfg, record)
    if not all(p.is_file() for p in paths.values()):
        raise RuntimeError("Blink-aid cache is missing/incomplete; prepare it on the cluster first.")
    metadata = json.loads(paths["metadata"].read_text())
    if metadata["fingerprint"] != digest(identity):
        raise RuntimeError("Blink-aid cache is stale: source, sensors, fit mask or fit settings changed; prepare only the affected cache again.")
    if metadata["model_sha256"] != file_hash(paths["model"]) or metadata["sources_stamp"] != file_stamp(paths["sources"]):
        raise RuntimeError("Blink-aid cache files changed after preparation.")
    sources = mne.io.read_raw_fif(paths["sources"], preload=False, verbose="ERROR")
    expected = identity["coordinates"]
    if any(coordinates(sources)[key] != expected[key] for key in ("sfreq", "n_times", "first_samp", "meas_date")):
        raise RuntimeError("Cached ICA sources no longer match the native recording grid.")
    model = mne.preprocessing.read_ica(paths["model"], verbose="ERROR")
    if model.ch_names != identity["channels"] or sources.ch_names != [f"ICA{i:03}" for i in range(model.n_components_)]:
        raise RuntimeError("Cached ICA channel/component basis changed.")
    return {"model": model, "native_sources": sources, "metadata": metadata}


def load_cached_aid(cfg, record):
    """Read compatible preparation outputs locally, without fitting/filtering MEG."""
    raw, source = load_aid_input(cfg, record)
    return read_cache(cfg, record, cache_identity(cfg, raw, source))


# -----------------------------------------------------------------------------
# Headless preparation; fit copies and native-grid projections stay separate
# -----------------------------------------------------------------------------

def prepare_aid(cfg, record, *, overwrite=False):
    """Prepare one cache, or reuse a fully compatible one; write metadata last."""
    raw, source_path = load_aid_input(cfg, record)
    identity = cache_identity(cfg, raw, source_path)
    paths = cache_paths(cfg, record)
    if any(p.exists() for p in paths.values()) and not overwrite:
        return {**read_cache(cfg, record, identity)["metadata"], "execution_status": "reused"}
    raw_root = cfg["resolved_paths"]["data_root"]
    for path in paths.values():
        prepare_output_path(path, raw_root=raw_root, overwrite=overwrite)
    frontal = raw.copy().pick(identity["channels"])
    fit, summary = prepare_fit_copy(frontal, cfg["ica"])
    rank = int(sum(mne.compute_rank(fit, rank=None, verbose="ERROR").values()))
    n_components = min(identity["n_components_max"], rank)
    if n_components < 2:
        raise RuntimeError("Frontal fit rank is below two; use manual ocular review for this run.")
    model = mne.preprocessing.ICA(n_components=n_components, method="fastica", random_state=97, max_iter="auto")
    model.fit(fit, reject_by_annotation=True)
    # Project separate native-rate stretches, trimming their full FIR support.
    # NaN gaps are unscored context, never new BAD/epoch-exclusion annotations.
    guard = int(np.ceil(summary["filter"]["half_support_sec"] * raw.info["sfreq"]))
    data = np.full((n_components, raw.n_times), np.nan)
    intervals = []
    for start, stop in good_stretches(frontal):
        if stop - start <= 2 * guard:
            continue
        segment = frontal.copy().crop(start / raw.info["sfreq"], (stop - 1) / raw.info["sfreq"]).load_data()
        segment.filter(cfg["ica"]["fit_l_freq_hz"], cfg["ica"]["fit_h_freq_hz"],
                       skip_by_annotation=("edge", "bad"), verbose="ERROR")
        projected = model.get_sources(segment).get_data()
        data[:, start + guard:stop - guard] = projected[:, guard:segment.n_times - guard]
        intervals.append({"start_raw": int(raw.first_samp + start + guard),
                          "stop_raw_exclusive": int(raw.first_samp + stop - guard)})
    info = mne.create_info([f"ICA{i:03}" for i in range(n_components)], raw.info["sfreq"], "misc")
    sources = mne.io.RawArray(data, info, first_samp=raw.first_samp, verbose="ERROR")
    sources.set_meas_date(raw.info["meas_date"])
    model.save(paths["model"], overwrite=overwrite)
    sources.save(paths["sources"], fmt="single", overwrite=overwrite, verbose="ERROR")
    metadata = {"purpose": "supplementary_annotation_only", "status": "prepared_unreviewed",
        "source_filename": record.source_filename, "generated_utc": utc_now(),
        "fingerprint": digest(identity), "identity": identity, "mne_version": mne.__version__,
        "model_sha256": file_hash(paths["model"]), "sources_stamp": file_stamp(paths["sources"]),
        "n_components_fitted": int(model.n_components_), "estimated_rank": rank,
        "converged": bool(model.n_iter_ < 1000), "fit_summary": summary,
        "projection_intervals": intervals, "unscored_samples": int(np.isnan(data[0]).sum()),
        "execution_status": "prepared"}
    write_json(paths["metadata"], metadata, raw_root=raw_root, overwrite=overwrite)
    return metadata


# -----------------------------------------------------------------------------
# Local peak proposals, with independent stretch filtering and exact coordinates
# -----------------------------------------------------------------------------

def finite_stretches(valid):
    """Return contiguous half-open spans; never concatenate gaps for detection."""
    cuts = np.flatnonzero(np.diff(np.r_[False, valid, False])).tolist()
    return list(zip(cuts[::2], cuts[1::2]))


def detect_blinks(sources, components, settings):
    """Return native-grid candidate intervals, robust scores and review metadata.

    Polarity does not matter. Each selected source is centered/scaled using
    finite samples after the additional low-pass guard. Peak width is a review
    flag only; no long event is silently discarded. Empty/flat sources produce
    an explicit diagnostic, not evidence that the recording contains no blinks.
    """
    sfreq = float(sources.info["sfreq"])
    scores = np.full(sources.n_times, np.nan)
    if not components:
        return [], scores, {"status": "no_component_selected", "candidate_count": 0, "settings": settings}
    design = filter_design(sfreq, None, 10.)
    guard = int(np.ceil(design["half_support_sec"] * sfreq))
    peaks_found, diagnostics = [], []
    for component in components:
        signal = sources.get_data(picks=[f"ICA{component:03}"])[0]
        filtered = np.full_like(signal, np.nan)
        for start, stop in finite_stretches(np.isfinite(signal)):
            if stop - start <= 2 * guard:
                continue
            piece = mne.filter.filter_data(signal[start:stop], sfreq, None, 10., verbose="ERROR")
            filtered[start + guard:stop - guard] = piece[guard:len(piece) - guard]
        valid = np.isfinite(filtered)
        scale = float(median_abs_deviation(filtered[valid], scale="normal")) if valid.any() else 0.
        if not np.isfinite(scale) or scale <= np.finfo(float).eps:
            diagnostics.append({"component": int(component), "status": "flat_or_unscorable"})
            continue
        score = np.abs((filtered - np.median(filtered[valid])) / scale)
        scores = np.fmax(scores, score)
        for start, stop in finite_stretches(np.isfinite(score)):
            peaks, properties = find_peaks(score[start:stop], height=settings["threshold_std"],
                distance=max(1, int(round(settings["min_peak_distance_sec"] * sfreq))), width=(None, None))
            for i, peak in enumerate(peaks + start):
                left = max(start, int(peak + round(settings["tmin_sec"] * sfreq)))
                right = min(stop, int(peak + round(settings["tmax_sec"] * sfreq)) + 1)
                peaks_found.append({"start": left, "stop": right, "peak_raw": int(peak + sources.first_samp),
                    "component": int(component), "score": float(score[peak]),
                    "long_event": bool(properties["widths"][i] / sfreq > settings["long_peak_width_sec"])})
        diagnostics.append({"component": int(component), "robust_std": scale, "status": "scored"})
    # Union duplicates only when their intervals overlap; the original peak
    # rows stay in the receipt for reviewing long events and multiple sources.
    intervals = []
    for row in sorted(peaks_found, key=lambda r: r["start"]):
        if intervals and row["start"] < intervals[-1][1]:
            intervals[-1][1] = max(intervals[-1][1], row["stop"])
        else:
            intervals.append([row["start"], row["stop"]])
    return intervals, scores, {"status": "proposals_for_review" if peaks_found else "no_peaks_review_still_required",
        "candidate_count": len(intervals), "settings": dict(settings), "components": list(components),
        "peaks": peaks_found, "diagnostics": diagnostics, "detection_filter": design,
        "unscored_samples": int(np.isnan(scores).sum())}


# -----------------------------------------------------------------------------
# Local duration proposals for step-like or sustained component excursions
# -----------------------------------------------------------------------------

def detect_eye_movements(sources, components, settings):
    """Return half-open native-index intervals, scores and eye-only provenance.

    A centered odd-window median preserves steps while suppressing short spikes.
    Center/scale each component over its finite, fully supported smoothed samples
    using median/MAD. Require a continuous absolute-score excursion of at least
    min_duration_sec BEFORE merging short gaps and padding within that same
    scored stretch. Merge overlapping component proposals, never NaN gaps.

    ICA sources already have a 1-40 Hz bandpass: long position plateaus/drift may
    be absent, and threshold durations are not measured eye-movement durations.
    Robust normalization assumes ordinary activity occupies most scored samples.
    Category follows the operator's component choice; mixed components can also
    propose blinks. Every candidate requires visual review.
    Flat/unscorable components and NONE provide no evidence of artifact absence.
    Source samples/annotations are never changed or resampled. Metadata endpoints
    include first_samp; returned intervals index the unchanged native source grid.
    """
    sfreq = float(sources.info["sfreq"])
    window = max(1, int(round(settings["median_window_sec"] * sfreq)))
    window += int(window % 2 == 0)
    guard = window // 2
    minimum = max(1, int(np.ceil(settings["min_duration_sec"] * sfreq)))
    gap = int(np.floor(settings["merge_gap_sec"] * sfreq))
    padding = int(round(settings["padding_sec"] * sfreq))
    scores = np.full(sources.n_times, np.nan)
    proposals, events, diagnostics = [], [], []
    for component in components:
        signal = sources.get_data(picks=[f"ICA{component:03}"])[0]
        smoothed = np.full_like(signal, np.nan)
        for start, stop in finite_stretches(np.isfinite(signal)):
            if stop - start <= 2 * guard:
                continue
            piece = median_filter(signal[start:stop], size=window, mode="nearest")
            # Discard every sample whose centered window would need padding.
            smoothed[start + guard:stop - guard] = piece[guard:len(piece) - guard] if guard else piece
        valid = np.isfinite(smoothed)
        scale = float(median_abs_deviation(smoothed[valid], scale="normal")) if valid.any() else 0.
        if not np.isfinite(scale) or scale <= np.finfo(float).eps:
            diagnostics.append({"component": int(component), "status": "flat_or_unscorable"})
            continue
        score = np.abs((smoothed - np.median(smoothed[valid])) / scale)
        scores = np.fmax(scores, score)
        for start, stop in finite_stretches(np.isfinite(score)):
            bouts = []
            for left, right in finite_stretches(score[start:stop] >= settings["threshold_std"]):
                left, right = start + left, start + right
                if right - left < minimum:
                    continue
                if bouts and left - bouts[-1][1] <= gap:
                    bouts[-1][1] = right
                else:
                    bouts.append([left, right])
            for left, right in bouts:
                padded = [max(start, left - padding), min(stop, right + padding)]
                proposals.append(padded)
                events.append({"component": int(component),
                    "start_raw": int(left + sources.first_samp),
                    "stop_raw_exclusive": int(right + sources.first_samp),
                    "proposal_start_raw": int(padded[0] + sources.first_samp),
                    "proposal_stop_raw_exclusive": int(padded[1] + sources.first_samp),
                    "max_score": float(np.max(score[left:right]))})
        diagnostics.append({"component": int(component), "robust_std": scale, "status": "scored"})
    intervals = []
    for start, stop in sorted(proposals):
        if intervals and start < intervals[-1][1]:
            intervals[-1][1] = max(intervals[-1][1], stop)
        else:
            intervals.append([start, stop])
    return intervals, scores, {"method": EYE_METHOD,
        "status": ("no_component_selected" if not components else
                   "proposals_for_review" if intervals else "no_excursions_review_still_required"),
        "candidate_count": len(intervals), "settings": dict(settings), "components": list(components),
        "events": events, "diagnostics": diagnostics, "median_window_samples": window,
        "smoothing_guard_samples": guard, "sfreq": sfreq, "first_samp": int(sources.first_samp),
        "unscored_samples": int(np.isnan(scores).sum()),
        "unscored_intervals_raw": [[int(a + sources.first_samp), int(b + sources.first_samp)]
                                   for a, b in finite_stretches(np.isnan(scores))]}


def annotation_key(raw, annotation):
    """Identify an edited interval by native sample endpoints and description."""
    start = int(round(annotation["onset"] * raw.info["sfreq"]))
    stop = start + int(round(annotation["duration"] * raw.info["sfreq"]))
    return (start, stop, str(annotation["description"]))


def remember_deletions(display, state, *, label=AUTO_LABEL):
    """Remember only this label's deleted/moved spans, without editing annotations."""
    present = {annotation_key(display, a) for a in display.annotations}
    removed = state.setdefault("removed_spans_raw", [])
    for key in state.get("previous_auto", []):
        span = list(key[:2])
        if key not in present and span not in removed:
            removed.append(span)
    state["previous_auto"] = [annotation_key(display, a) for a in display.annotations if a["description"] == label]


def replace_proposals(display, intervals, state, *, label=AUTO_LABEL):
    """Replace untouched session proposals, retaining manual work and deletions.

    Previously saved annotations are never considered disposable session
    proposals. Edited generated intervals no longer match their original key,
    so they survive regeneration. Remembered deletions suppress overlapping new
    proposals; the operator may still draw an interval back manually.
    ``label`` defaults to the historical blink label. The other automatic eye
    category may overlap; manual/severe BAD overlap protection remains in force.
    """
    remember_deletions(display, state, label=label)
    generated = {tuple(key) for key in state.get("generated", [])}
    indices = [i for i, a in enumerate(display.annotations) if annotation_key(display, a) in generated]
    display.annotations.delete(indices)
    # Keep the original blink overlap policy for manual/severe BAD intervals.
    # The other automatic ocular group never suppresses this group's category.
    other = EYE_AUTO_LABEL if label == AUTO_LABEL else AUTO_LABEL
    existing = [(key[0], key[1]) for a in display.annotations
                if a["description"].lower().startswith("bad") and a["description"] != other
                for key in [annotation_key(display, a)]]
    generated = []
    for start, stop in intervals:
        absolute = [int(start + display.first_samp), int(stop + display.first_samp)]
        if any(absolute[0] < b and absolute[1] > a for a, b in existing + state["removed_spans_raw"]):
            continue
        display.annotations.append(absolute[0] / display.info["sfreq"],
                                   (stop - start) / display.info["sfreq"], label)
        generated.append((*absolute, label))
    state["generated"] = generated
    state["previous_auto"] = [annotation_key(display, a) for a in display.annotations if a["description"] == label]


# -----------------------------------------------------------------------------
# Cluster entry point: prepare only, with per-run failure reporting and no GUI
# -----------------------------------------------------------------------------

def main(argv=None):
    """Prepare selected participant acquisitions; never write review receipts."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for flag in ("config", "subject-config", "tuning-config"):
        parser.add_argument("--" + flag, required=True)
    parser.add_argument("--site-config")
    parser.add_argument("--source", help="One exact filename; otherwise prepare all configured acquisitions.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true")
    mode.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    cfg, *_ = load_and_resolve_config(args.config, args.subject_config, args.tuning_config,
        site_config_path=args.site_config, required_sections=("artifact_annotation", "ica"))
    dirs = stage_output_dirs(cfg)
    require_review(cfg, "03_A_MaxFilt", dirs["reports"])
    records = discover_raw_run_records(cfg)
    if args.source:
        records = [r for r in records if r.source_filename == args.source]
        if len(records) != 1:
            raise ValueError("--source must match one exact configured acquisition filename.")
    results = []
    for record in records:
        try:
            result = prepare_aid(cfg, record, overwrite=args.overwrite)
            results.append({"source_filename": record.source_filename, "status": result["execution_status"],
                "fingerprint": result["fingerprint"], "n_components": result["n_components_fitted"],
                "converged": result["converged"]})
        except Exception as exc:
            # Keep processing independent runs, while ensuring the job reports
            # failure rather than misrepresenting an unavailable aid as no blinks.
            results.append({"source_filename": record.source_filename, "status": "failed", "error": str(exc)})
        print(results[-1], flush=True)
    failed = any(row["status"] == "failed" for row in results)
    report = {"stage": "4B_supplementary_blink_aid", "scientific_acceptance": False,
              "generated_utc": utc_now(), "runs": results, "status": "failed" if failed else "prepared_unreviewed"}
    write_json(dirs["reports"] / "03_B_blink_aid_preparation.json", report,
               raw_root=cfg["resolved_paths"]["data_root"], overwrite=True)
    return 2 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
