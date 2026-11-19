"""Stage 4B EyeLink guidance: trial-local candidates and disposable EDF traces.

Usage
-----
The 4B runner calls prepare_candidates once per exact MEG source. Re-run with:
    DAD_STAGES=artifact_annotation sbatch --array=4-39%4 meg/cluster/submit_preprocessing_review.sh
Then open manual_eye_annotations.py with the usual exact --source and tuning.
EDF traces are shown by default; --no-edf-traces displays candidates alone.

Inputs are the accepted source manifest/canonical trigger CSV, original EDFs,
and accepted SSS coordinates. Output sidecars contain non-BAD annotations, a
blink CSV and trial/mapping provenance JSON. No MEG signals are changed here.
Only explicit manual SAVE promotes retained candidates to BAD_blink elsewhere.

The exploration in tmp_test_align_EDF_FIF supplies the trial-local approach,
not a clock calibration or ICA model. EDF recording-off gaps must survive:
eyelinkio.read_edf compresses its time axis, so this small pinned EDFAPI adapter
reads original timestamps and recording START/STOP records instead. Samples
are read only for a local display, never for batch candidate generation.
Mapping assumes equal elapsed seconds within each short matched trial; it is
guidance, not proven continuous synchronization. No global/physiological offset,
drift estimate, Stage 3 timing correction or waveform alignment is applied.
"""

from __future__ import annotations

import csv
import ctypes
import json
import re
from collections import Counter
from pathlib import Path

import mne
import numpy as np

from continuous_cleaning import coordinates, digest, file_hash, file_stamp
from preprocessing_stage_utils import guarded_path, utc_now

ONSET_CODES = {71, 81, 91}
CANDIDATE_LABELS = ("EDF_blink_candidate", "EDF_blink_long_candidate")
TRIGGER = re.compile(r"^TRIGGER\s+(\d+)\s+trial\s+(\d+)\s+frame\s+(\d+)\s+label\s+(.*)$")
TRIAL = re.compile(r"^TRIALID\s+(\d+)\s*$")


# -----------------------------------------------------------------------------
# Original-clock EDF decoding; keep this version-specific adapter in one place
# -----------------------------------------------------------------------------

def reader_api():
    """Validate the pinned reader/native runtime without opening a recording."""
    from importlib.metadata import version, PackageNotFoundError
    try:
        if version("eyelinkio") != "0.3.0":
            raise RuntimeError("Stage 4B EDF guidance requires eyelinkio==0.3.0.")
        from eyelinkio.edf import _edf2py as api, _defines as defs
    except (PackageNotFoundError, ImportError, OSError, AssertionError) as exc:
        raise RuntimeError("EDF reader unavailable: install eyelinkio==0.3.0 and check its EDFAPI import before 4B.") from exc
    return api, defs


def read_edf_records(path, *, samples=False):
    """Read messages, blink endpoints and recording blocks in original seconds.

    eyelinkio 0.3.0 bundles the native EDFAPI library on Linux/macOS. Its small
    internal ctypes adapter is pinned and isolated here because the public
    reader discards recording gaps. Events-only opens do not load gaze samples.
    With samples=True return [seconds, gaze_y_L, pupil_L, gaze_y_R, pupil_R].
    Missing fields/eyes and zero (invalid) pupil sizes become NaN, not zeros.
    """
    api, defs = reader_api()
    error = ctypes.c_int(0)
    handle = api.edf_open_file(str(Path(path).resolve()).encode(), 2, 1, int(samples), ctypes.byref(error))
    if not handle or error.value:
        if handle:
            api.edf_close_file(handle)
        raise RuntimeError(f"Cannot open EDF {Path(path).name}: EDFAPI error {error.value}")
    result = {"messages": [], "blinks": [], "blocks": [], "samples": []}
    active = None
    try:
        while True:
            kind = api.edf_get_next_data(handle)
            if kind == defs.event_constants["NO_PENDING_ITEMS"]:
                break
            if kind not in (24, 4, 30, 200):
                continue
            data = api.edf_get_float_data(handle).contents
            if kind == 24:  # MESSAGEEVENT: copy bytes before the next SDK call.
                event = data.fe
                message = ctypes.string_at(ctypes.byref(event.message[0]), event.message.contents.len + 2)[2:]
                result["messages"].append({"time": event.sttime / 1000.,
                                           "text": message.rstrip(b"\0").decode("utf8", errors="replace").strip()})
            elif kind == 4:  # ENDBLINK includes the original start and end.
                event = data.fe
                result["blinks"].append({"start": event.sttime / 1000., "stop": event.entime / 1000.,
                                          "eye": int(event.eye)})
            elif kind == 30:  # RECORDING_INFO: preserve real recording-off gaps.
                rec = data.rec
                if rec.state and active is None:
                    active = {"start": rec.time / 1000., "stop": None, "sfreq": float(rec.sample_rate)}
                    result["blocks"].append(active)
                elif not rec.state and active is not None:
                    active["stop"] = rec.time / 1000.
                    active = None
            elif samples:
                sample = data.fs
                # Half-millisecond timestamps are an EDF format flag, not a lag correction.
                row = [(sample.time + (0.5 if sample.flags & defs.SAMPLE_ADD_OFFSET else 0.)) / 1000.]
                for eye, eye_flag in [(0, defs.SAMPLE_LEFT), (1, defs.SAMPLE_RIGHT)]:
                    gaze = float(sample.gy[eye]) if sample.flags & eye_flag and sample.flags & defs.SAMPLE_GAZEXY else np.nan
                    pupil = float(sample.pa[eye]) if sample.flags & eye_flag and sample.flags & defs.SAMPLE_PUPILSIZE else np.nan
                    row.extend([gaze if -32768 < gaze < 1e8 else np.nan,
                                pupil if 0 < pupil < 1e8 else np.nan])
                result["samples"].append(row)
    finally:
        api.edf_close_file(handle)
    if samples:
        result["samples"] = np.asarray(result["samples"], dtype=float).reshape(-1, 5)
    return result


# -----------------------------------------------------------------------------
# Trial correspondence: accepted MEG associations, local markers, no clock fit
# -----------------------------------------------------------------------------

def edf_trials(records):
    """Split messages by TRIALID occurrence; repeated numeric trial IDs are valid."""
    trials = []
    for message in records["messages"]:
        trial = TRIAL.match(message["text"])
        marker = TRIGGER.match(message["text"])
        if trial:
            trials.append({"ordinal": len(trials), "trial_id": int(trial[1]),
                           "start": message["time"], "markers": []})
        elif marker and trials:
            trials[-1]["markers"].append({"time": message["time"], "code": int(marker[1]),
                                          "trial_id": int(marker[2]), "frame": int(marker[3]), "label": marker[4]})
    for index, trial in enumerate(trials):
        trial["stop"] = trials[index + 1]["start"] if index + 1 < len(trials) else float("inf")
        onsets = [m for m in trial["markers"] if m["code"] in ONSET_CODES and m["label"] == "trial_condition_onset"]
        trial["anchor"] = onsets[0] if len(onsets) == 1 else None
    return trials


def meg_trials(event_rows):
    """Split canonical observed triggers without editing accepted associations."""
    trials = []
    for row in event_rows:
        marker = {"code": int(row["trigger_code"]), "frame": int(row.get("behavioral_frame") or -1),
                  "sample": int(row["raw_sample_frame"]), "label": row.get("reported_label", ""),
                  "unresolved": str(row.get("unresolved_transition", "False")).lower() == "true"}
        if marker["code"] in ONSET_CODES:
            trials.append({"anchor": marker, "source": row["source_fif"],
                           "block": int(row.get("behavioral_block") or -1),
                           "trial_id": int(row.get("behavioral_trial") or -1), "markers": []})
        if trials and row["source_fif"] == trials[-1]["source"]:
            trials[-1]["markers"].append(marker)
    return trials


def trial_signature(trial, identity_range):
    """Use onset/identity/frame evidence, not filename order, to align trials."""
    anchor = trial.get("anchor")
    if anchor is None:
        return None
    markers = [m for m in trial["markers"] if m["code"] in ONSET_CODES
               or identity_range[0] <= m["code"] <= identity_range[1]]
    if len(markers) < 2:
        return None
    return (trial["trial_id"], tuple((m["code"], m["frame"]) for m in markers))


def match_trials(eyes, megs, identity_range, sfreq, tolerance):
    """Match exact sequences, or uniquely anchored neighbouring trial contexts.

    After a mismatch, only a unique three-trial context in both streams is used.
    Missing/repeated trials cannot shift subsequent ordinal pairings. Full EDF
    marker keys must be an unambiguous ordered subset of the matched MEG trial;
    extra MEG events remain reported. Interval disagreement with the local onset
    flags a trial; it is never compensated by an offset or clock-rate estimate.
    """
    left = [trial_signature(t, identity_range) for t in eyes]
    right = [trial_signature(t, identity_range) for t in megs]
    if left == right and all(s is not None for s in left):
        proposed = dict(enumerate(range(len(eyes))))
    else:
        def contexts(sequence):
            return [(sequence[i - 1] if i else "START", value,
                     sequence[i + 1] if i + 1 < len(sequence) else "END")
                    for i, value in enumerate(sequence)]
        a, b = contexts(left), contexts(right)
        ca, cb = Counter(a), Counter(b)
        proposed = {i: b.index(key) for i, key in enumerate(a)
                    if left[i] is not None and ca[key] == cb[key] == 1}
        # Crossed unique contexts are contradictory evidence, not permission to reorder.
        proposed = {i: j for i, j in proposed.items()
                    if all((k - i) * (v - j) > 0 for k, v in proposed.items() if k != i)}
    result = []
    for i, eye in enumerate(eyes):
        row = {"edf_trial_ordinal": i, "trial_id": eye["trial_id"], "status": "ambiguous_or_unmatched"}
        if i in proposed:
            meg = megs[proposed[i]]
            keys = [(m["code"], m["frame"]) for m in meg["markers"]]
            positions = []; good = eye["anchor"] is not None
            for marker in eye["markers"]:
                candidates = [j for j, key in enumerate(keys) if key == (marker["code"], marker["frame"])]
                if len(candidates) != 1 or marker["trial_id"] != eye["trial_id"]:
                    good = False; break
                positions.append(candidates[0])
                delta = ((meg["markers"][candidates[0]]["sample"] - meg["anchor"]["sample"]) / sfreq
                         - (marker["time"] - eye["anchor"]["time"]))
                if abs(delta) > tolerance or meg["markers"][candidates[0]]["unresolved"]:
                    good = False
            good = good and bool(positions) and positions == sorted(set(positions))
            row.update(status="matched" if good else "mismatched_markers", source=meg["source"],
                       edf_anchor_sec=eye["anchor"]["time"], meg_anchor_sample=meg["anchor"]["sample"],
                       extra_meg_marker_count=len(keys) - len(positions))
        result.append(row)
    for index, meg in enumerate(megs):
        if index not in proposed.values():
            result.append({"edf_trial_ordinal": None, "trial_id": meg["trial_id"],
                           "source": meg["source"], "meg_anchor_sample": meg["anchor"]["sample"],
                           "status": "meg_trial_without_unique_edf_match"})
    return result


def map_blinks(records, eyes, matches, source, sfreq, first_samp, n_times, long_sec):
    """Map both endpoints within one matched recorded trial; never bridge gaps."""
    rows = []; segments = []
    for trial, match in zip(eyes, matches):
        if match.get("source") != source or match["status"] != "matched":
            continue
        # START/STOP must be known, contain the anchor, and stay in this trial.
        blocks = [b for b in records["blocks"] if b["stop"] is not None
                  and trial["start"] <= b["start"] <= match["edf_anchor_sec"] < b["stop"] <= trial["stop"]]
        if len(blocks) != 1:
            match["status"] = "unknown_or_ambiguous_recording_bounds"
            continue
        block = blocks[0]
        segments.append({**match, **block})
    for number, blink in enumerate(records["blinks"]):
        trial_indices = [i for i, trial in enumerate(eyes) if trial["start"] <= blink["start"] < trial["stop"]]
        row = {"event_id": number, "eye": blink["eye"], "edf_start_sec": blink["start"],
               "edf_stop_sec": blink["stop"], "duration_sec": blink["stop"] - blink["start"],
               "long_event": blink["stop"] - blink["start"] > long_sec, "status": "outside_matched_recording"}
        if blink["stop"] <= blink["start"]:
            row["status"] = "invalid_edf_interval"
            rows.append(row)
            continue
        if len(trial_indices) == 1:
            index = trial_indices[0]; match = matches[index]
            row.update(edf_trial_ordinal=index, trial_id=eyes[index]["trial_id"])
            if match.get("source") not in (None, source):
                continue
            segments_here = [s for s in segments if s["edf_trial_ordinal"] == index
                             and s["start"] <= blink["start"] < blink["stop"] <= s["stop"]]
            if len(segments_here) == 1:
                seg = segments_here[0]
                a, b = [int(seg["meg_anchor_sample"] + np.rint((t - seg["edf_anchor_sec"]) * sfreq))
                        for t in [blink["start"], blink["stop"]]]
                if first_samp <= a < b <= first_samp + n_times:
                    row.update(status="candidate", sample_start_raw=a, sample_stop_raw=b)
                else:
                    row["status"] = "outside_meg_recording"
            elif match["status"] != "matched":
                row["status"] = match["status"]
        rows.append(row)
    return rows, segments


# -----------------------------------------------------------------------------
# Focused sidecars: independent from muscle FIFs and human annotation files
# -----------------------------------------------------------------------------

def candidate_paths(cfg, record):
    """Use the existing preprocessing tree and complete raw source stem."""
    folder = Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "edf_blink_candidates"
    return {"annotations": folder / f"{record.source_stem}-annot.fif",
            "csv": folder / f"{record.source_stem}_blinks.csv", "json": folder / f"{record.source_stem}_mapping.json"}


def load_context(cfg):
    """Read saved Stage 1 evidence once; absence is explicit, never reviewed zero."""
    root = Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"])
    events_path = root / "events/02_trigger_events.csv"
    manifest_path = Path(cfg["resolved_paths"]["derivatives_reports_dir"]) / "source_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    events = []
    if events_path.exists():
        with events_path.open(newline="") as handle:
            events = list(csv.DictReader(handle))
    return {"manifest": manifest, "events": events, "events_path": events_path,
            "events_sha256": file_hash(events_path) if events_path.exists() else None}


def candidate_input(cfg, record, raw, sss_path, context):
    """Fingerprint original EDFs, canonical triggers and exact SSS coordinates."""
    manifest = context["manifest"]
    associations = manifest.get("trigger_derived_associations", {}).get("meg_acquisition_run_to_behavioral_segments", [])
    blocks = {int(s["behavioral_block"]) for a in associations if a["source_fif"] == record.source_filename
              for s in a["behavioral_segments"]}
    files = []
    for row in manifest.get("sources", {}).get("eyetracker_edf_files", []):
        # Filename narrows possible block files; the trial/marker checks establish correspondence.
        name = row["source_filename"]
        if any(name.lower() == f'MB{int(cfg["subject"]["id"]):02d}{block:02d}00.edf'.lower() for block in blocks):
            path = Path(cfg["resolved_paths"]["data_root"]) / row["raw_relative_path"]
            files.append({"relative_path": row["raw_relative_path"], "stamp": file_stamp(path)})
    signature = {"sss_input": file_stamp(sss_path), "coordinates": coordinates(raw),
                 "edf_inputs": files, "events_sha256": context["events_sha256"],
                 "association_sha256": digest(associations), "settings": cfg["artifact_annotation"]["edf_blink_candidates"],
                 "implementation_sha256": file_hash(__file__)}
    return signature, blocks


def prepare_candidates(cfg, record, raw, sss_path, context, cache):
    """Generate/cache guidance only; unchanged muscle outputs need no rebuild."""
    settings = cfg["artifact_annotation"].get("edf_blink_candidates", {})
    if not settings.get("enabled", False):
        return {"status": "disabled"}
    signature, blocks = candidate_input(cfg, record, raw, sss_path, context)
    fingerprint = digest(signature); paths = candidate_paths(cfg, record)
    if paths["json"].exists():
        saved = json.loads(paths["json"].read_text())
        if saved.get("fingerprint") == fingerprint and all(paths[k].exists() and file_hash(paths[k]) == saved[k + "_sha256"]
                                                        for k in ("annotations", "csv")):
            return {**saved, "paths": {k: str(p) for k, p in paths.items()}}
    rows = []; trials = []; segments = []
    for source in signature["edf_inputs"]:
        path = Path(cfg["resolved_paths"]["data_root"]) / source["relative_path"]
        if sum(Path(s["relative_path"]).stem[4:6] == path.stem[4:6] for s in signature["edf_inputs"]) != 1:
            trials.append({"edf_relative_path": source["relative_path"], "status": "ambiguous_edf_block_files"})
            continue
        if str(path) not in cache:
            cache[str(path)] = read_edf_records(path)
        records = cache[str(path)]; eyes = edf_trials(records)
        block = int(path.stem[4:6])
        events = [r for r in context["events"] if int(r.get("behavioral_block") or -1) == block]
        megs = meg_trials(events)
        matches = match_trials(eyes, megs, cfg["subject_profile"]["identity_range"], raw.info["sfreq"], settings["marker_tolerance_sec"])
        blink_rows, mapped_segments = map_blinks(records, eyes, matches, record.source_filename,
            raw.info["sfreq"], raw.first_samp, raw.n_times, settings["long_event_sec"])
        for collection in [blink_rows, matches, mapped_segments]:
            for row in collection:
                row["edf_relative_path"] = source["relative_path"]
        rows.extend(blink_rows); trials.extend(matches); segments.extend(mapped_segments)
    # Merge overlapping binocular bouts within the same recorded trial only;
    # keep every original eye event and endpoint in the CSV for inspection.
    bouts = []
    for row in sorted((r for r in rows if r["status"] == "candidate"), key=lambda r: r["sample_start_raw"]):
        key = (row["edf_relative_path"], row["edf_trial_ordinal"])
        if bouts and bouts[-1]["key"] == key and row["sample_start_raw"] <= bouts[-1]["stop"]:
            bouts[-1]["stop"] = max(bouts[-1]["stop"], row["sample_stop_raw"])
            bouts[-1]["long"] |= row["long_event"]
        else:
            bouts.append({"key": key, "start": row["sample_start_raw"], "stop": row["sample_stop_raw"], "long": row["long_event"]})
    sfreq = raw.info["sfreq"]
    ann = mne.Annotations([b["start"] / sfreq for b in bouts], [(b["stop"] - b["start"]) / sfreq for b in bouts],
                          [CANDIDATE_LABELS[int(b["long"] or (b["stop"] - b["start"]) / sfreq > settings["long_event_sec"])] for b in bouts],
                          orig_time=raw.annotations.orig_time)
    for path in paths.values():
        guarded_path(path, cfg["resolved_paths"]["data_root"]).parent.mkdir(parents=True, exist_ok=True)
    ann.save(paths["annotations"], overwrite=True)
    columns = ["edf_relative_path", "event_id", "eye", "edf_trial_ordinal", "trial_id", "edf_start_sec", "edf_stop_sec", "duration_sec", "long_event", "status", "sample_start_raw", "sample_stop_raw"]
    with paths["csv"].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns); writer.writeheader(); writer.writerows(rows)
    result = {**signature, "fingerprint": fingerprint, "generated_utc": utc_now(),
              "status": "candidates_for_manual_review" if segments else "unavailable_or_unmatched",
              "candidate_count": len(ann), "long_candidate_count": sum(ann.description == CANDIDATE_LABELS[1]),
              "matched_trial_count": len(segments),
              "flagged_trial_count": sum(t["status"] != "matched" and t.get("source") in (None, record.source_filename) for t in trials),
              "trials": trials, "segments": segments,
              "annotations_sha256": file_hash(paths["annotations"]), "csv_sha256": file_hash(paths["csv"])}
    result["long_candidate_count"] = int(result["long_candidate_count"])
    paths["json"].write_text(json.dumps(result, indent=2))
    return {**result, "paths": {k: str(p) for k, p in paths.items()}}


def load_candidates(cfg, record, raw, sss_path):
    """Require a current 4B sidecar before displaying/promoting any candidates."""
    path = candidate_paths(cfg, record)["json"]
    if not path.exists():
        raise RuntimeError("EDF candidate evidence missing; rerun Stage 4B before ocular review.")
    saved = json.loads(path.read_text()); context = load_context(cfg)
    signature, _ = candidate_input(cfg, record, raw, sss_path, context)
    paths = candidate_paths(cfg, record)
    if saved["fingerprint"] != digest(signature) or any(file_hash(paths[k]) != saved[k + "_sha256"] for k in ("annotations", "csv")):
        raise RuntimeError("EDF candidate evidence changed; rerun Stage 4B before ocular review.")
    return saved, mne.read_annotations(paths["annotations"])


def promote_candidates_after_review(raw, *, confirmed):
    """Promote edited, retained candidates only after explicit human confirmation."""
    if confirmed is not True:
        raise ValueError("Explicit complete-source review is required before promoting EDF candidates.")
    mapping = {label: "BAD_blink" for label in CANDIDATE_LABELS if label in raw.annotations.description}
    if mapping:
        raw.annotations.rename(mapping)


# -----------------------------------------------------------------------------
# Browser-only samples: nearest measured sample inside each valid trial/block
# -----------------------------------------------------------------------------

def display_samples(raw, proposal, records_by_file):
    """Return normalized pupil/gaze plus coverage on the native MEG grid.

    NaNs break unrecorded/missing intervals. Nearest-neighbour placement is only
    for display and never crosses a recording boundary or fills sample gaps.
    Robust normalization is per trace, not a conversion of native EDF units.
    """
    traces = np.full((4, raw.n_times), np.nan); coverage = np.zeros(raw.n_times)
    for segment in proposal["segments"]:
        records = records_by_file[segment["edf_relative_path"]]
        samples = records["samples"]
        samples = samples[(samples[:, 0] >= segment["start"]) & (samples[:, 0] < segment["stop"])]
        if not len(samples):
            continue
        a, b = [int(segment["meg_anchor_sample"] + np.rint((t - segment["edf_anchor_sec"]) * raw.info["sfreq"]) - raw.first_samp)
                for t in (segment["start"], segment["stop"])]
        indices = np.arange(max(0, a), min(raw.n_times, b))
        times = segment["edf_anchor_sec"] + (indices + raw.first_samp - segment["meg_anchor_sample"]) / raw.info["sfreq"]
        right = np.clip(np.searchsorted(samples[:, 0], times), 0, len(samples) - 1)
        left = np.maximum(0, right - 1)
        nearest = np.where(abs(samples[left, 0] - times) < abs(samples[right, 0] - times), left, right)
        valid = abs(samples[nearest, 0] - times) <= 0.51 / segment["sfreq"]
        traces[:, indices[valid]] = samples[nearest[valid], 1:].T
        coverage[indices[valid]] = 1.
    scales = []
    for trace in traces:
        valid = np.isfinite(trace)
        center = float(np.median(trace[valid])) if valid.any() else 0.
        scale = float(np.percentile(abs(trace[valid] - center), 95)) if valid.any() else 1.
        scale = scale if scale > 0 else 1.
        trace[valid] = (trace[valid] - center) / scale
        scales.append({"median": center, "p95_absolute_deviation": scale})
    return np.vstack([traces, coverage]), scales


def add_edf_display(display, cfg, proposal):
    """Add temporary EDF MISC channels in arbitrary units; originals stay intact."""
    records = {s["edf_relative_path"]: None for s in proposal["segments"]}
    for relative in records:
        records[relative] = read_edf_records(Path(cfg["resolved_paths"]["data_root"]) / relative, samples=True)
    data, scales = display_samples(display, proposal, records)
    # Long zero/NaN gaps can make MNE's automatic MISC IQR equal zero. Use one
    # explicit display gain and enlarge only the synthetic traces to match it.
    # Original analog MISC values are never normalized or changed.
    display.load_data()  # Load this disposable channel subset once, not once per trace.
    gain = 1.
    original_misc = [ch for ch, kind in zip(display.ch_names, display.get_channel_types()) if kind == "misc"]
    if original_misc:
        values = display.get_data(picks=original_misc)[:, ::max(1, display.n_times // 10000)]
        for trace in values:
            finite = trace[np.isfinite(trace)]
            if len(finite):
                gain = max(gain, float(np.percentile(abs(finite - np.median(finite)), 95)))
    # Short labels fit the raw browser: Y=vertical gaze, P=pupil, L/R=eye.
    names = ["EDF_Y_L_norm", "EDF_P_L_norm", "EDF_Y_R_norm", "EDF_P_R_norm", "EDF_coverage"]
    available = [i for i in range(4) if np.isfinite(data[i]).any()] + [4]
    shown_names = [names[i] for i in available]
    # MISC supplies arbitrary browser units. BIO would incorrectly label these
    # dimensionless display traces as microvolts in MNE's public Raw browser.
    info = mne.create_info(shown_names, display.info["sfreq"], ["misc"] * len(shown_names))
    extra = mne.io.RawArray(data[available] * gain, info, first_samp=display.first_samp, verbose="ERROR")
    extra.set_meas_date(display.info["meas_date"])
    display.load_data().add_channels([extra], force_update_info=True)
    # These labels are also informational. Coverage=0 means no EDF guidance,
    # not absence of artifacts. NaN signal gaps remain visible in the browser.
    missing = data[-1] == 0
    cuts = np.diff(np.r_[False, missing, False].astype(int))
    starts, stops = np.flatnonzero(cuts == 1), np.flatnonzero(cuts == -1)
    display.annotations.append(starts / display.info["sfreq"] + display.first_time,
                               (stops - starts) / display.info["sfreq"], ["EDF_no_guidance"] * len(starts))
    return {"purpose": "discontinuous_visual_guidance_only", "channels": shown_names,
            "normalization": dict(zip(names[:4], scales)), "display_gain": gain}
