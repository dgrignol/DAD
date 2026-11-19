#!/usr/bin/env python3
"""Stage 4A SQUID candidate detection, review storage and SSS application.

Usage
-----
Preparation from already accepted channel decisions (no SSS or eye-file writes)::

    python meg/code/preprocessing_basic/review_squid_jumps.py --mode prepare \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml

The companion CLI reviews one exact acquisition using MNE's native browser.
03_A_MaxFilt.py calls these helpers between accepted channels and SSS. Inputs
are native raw FIFs, the final accepted bad-channel basis and detector settings.
Outputs are content-addressed candidate JSON/CSV and immutable reviewed JSON,
referenced by exact-source entries in the existing tuning YAML. This module
never opens or writes manual_eye annotation files, models or timing bundles.

Data flow: per-channel median-filtered first differences -> robust channel
scores -> merged time candidates -> explicit manual decisions -> BAD_squid_jumps
and optional channel exclusions. Scores are screening aids, not hardware
diagnoses. A review window is editable and is not a measured recovery duration.
Only operative exclusions enter SSS identities; empty reviews change no samples.
All intervals use half-open absolute native samples, including first_samp.
"""

from __future__ import annotations

import copy
import csv
import json
import math
import os
from pathlib import Path
import tempfile

import mne
import numpy as np
from scipy.ndimage import median_filter

from continuous_cleaning import coordinates, digest, file_hash, file_stamp, require_review, stage_configuration
from preprocessing_stage_utils import guarded_path, update_tuning_file, utc_now

LABEL = "BAD_squid_jumps"
CANDIDATE_LABEL = "possible_squid_jumps"
ALGORITHM = "median_absdiff_raw_mad_floor_peak_groups_v2"
DEFAULTS = {"enabled": False, "median_window_sec": 0.009, "threshold": 10.,
            "merge_gap_sec": 0.1, "review_padding_sec": 0.5, "block_sec": 30.}


# -----------------------------------------------------------------------------
# Scientific settings and operative decisions (independent of review UI state)
# -----------------------------------------------------------------------------

def parameters(cfg):
    """Validate configurable detector values; missing settings retain old CLIs."""
    result = {key: cfg.get("squid_jumps", {}).get(key, value) for key, value in DEFAULTS.items()}
    if not isinstance(result["enabled"], bool):
        raise ValueError("squid_jumps.enabled must be boolean.")
    for key in DEFAULTS.keys() - {"enabled"}:
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"squid_jumps.{key} must be finite numeric.")
        if value < 0 or (key in ("threshold", "block_sec") and value == 0):
            raise ValueError(f"Invalid squid_jumps.{key}.")
    return result


def application(cfg):
    """Return only nonempty saved actions, excluding thresholds/GUI/provenance.

    Receipt contents are validated against immutable JSON before application.
    Keeping empty reviews out of this projection preserves exact legacy SSS
    configuration/signatures and existing downstream manual review bindings.
    """
    result = {}
    for source, receipt in sorted(cfg.get("squid_jumps", {}).get("reviews", {}).items()):
        action = receipt.get("application", {})
        if action.get("intervals") or action.get("bad_channels"):
            result[source] = {"intervals": action.get("intervals", []),
                              "bad_channels": action.get("bad_channels", []),
                              "channel_scope": action.get("channel_scope", "union")}
    return result


def run_action(cfg, source):
    """Return this source's intervals and source/global channel additions."""
    actions = application(cfg)
    channels = set()
    for name, action in actions.items():
        if name == source or action["channel_scope"] == "union":
            channels.update(action["bad_channels"])
    return {"intervals": actions.get(source, {}).get("intervals", []),
            "bad_channels": sorted(channels)}


def attach_intervals(raw, intervals):
    """Attach global BAD windows directly in acquisition coordinates once.

    Native samples are half-open. No sample is removed, zeroed or interpolated.
    orig_time=None still uses raw.first_time for directly appended annotations.
    """
    existing = [(round(a["onset"] * raw.info["sfreq"]),
                 round((a["onset"] + a["duration"]) * raw.info["sfreq"]))
                for a in raw.annotations if a["description"] == LABEL]
    for start, stop in intervals:
        if not raw.first_samp <= start < stop <= raw.first_samp + raw.n_times:
            raise ValueError("SQUID interval is outside the exact native recording.")
        if (start, stop) not in existing:
            raw.annotations.append(start / raw.info["sfreq"], (stop-start) / raw.info["sfreq"], LABEL)
            existing.append((start, stop))


def without_owned_context(annotations):
    """Remove copied SQUID context from a Stage 4B editable snapshot.

    The authoritative intervals always come from the current accepted SSS
    source. Older eye snapshots can carry stale copies after revalidation;
    these copies must neither duplicate nor replace Stage 4A decisions.
    """
    return annotations[[i for i, label in enumerate(annotations.description)
                        if str(label).lower() != LABEL.lower()]].copy()


def maxwell_filter_with_jumps(raw, **kwargs):
    """SSS good stretches independently, retaining skipped native samples.

    MNE 1.10 derives a full-recording buffer for st_duration=None, then rejects
    that buffer when skip_by_annotation shortens the longest good stretch.
    SSS-only has no temporal fit, so filtering isolated stretches implements
    the documented skip semantics without inventing tSSS or cropping outputs.
    Metadata comes from the first filtered stretch; original annotations and
    every auxiliary sample remain on the original full native grid. Skipped
    MEG samples remain acquired, untransformed and BAD, never usable as SSS data.
    """
    from mne.annotations import _annotations_starts_stops
    if kwargs.get("st_duration") is not None:
        raise ValueError("DAD SQUID skipping supports SSS only (st_duration=None).")
    skips = ("edge", "bad_acq_skip", LABEL)
    starts, stops = _annotations_starts_stops(raw, skips, invert=True)
    bounds = [(int(a), int(b)) for a, b in zip(starts, stops) if b > a]
    if bounds == [(0, raw.n_times)]:
        return mne.preprocessing.maxwell_filter(raw, **{**kwargs, "skip_by_annotation": skips})
    if not bounds:
        raise ValueError("No good samples remain for SSS after SQUID/acquisition exclusions.")
    result = raw.copy().load_data()
    picks = mne.pick_types(raw.info, meg=True, ref_meg=True, exclude=[])
    template = None
    for start, stop in bounds:
        segment = raw.copy().crop(start/raw.info["sfreq"], (stop-1)/raw.info["sfreq"])
        # Boundaries have already been respected. Disposable segment labels
        # can be clipped to zero length by crop; none is an output annotation.
        segment.set_annotations(None)
        filtered = mne.preprocessing.maxwell_filter(segment, **{**kwargs, "skip_by_annotation": skips})
        if filtered.ch_names != raw.ch_names:
            raise RuntimeError("SSS changed the channel basis unexpectedly.")
        result._data[picks, start:stop] = filtered.get_data(picks=picks)
        if template is None:
            template = filtered.info.copy()
    result.info = template
    return result


def root(cfg):
    """Guard a separate derivative area; manual_eye files are never outputs."""
    return guarded_path(Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "squid_jumps",
                        cfg["resolved_paths"]["data_root"])


def channel_basis_configuration(cfg):
    """Bind detector exclusions to channel/SSS settings before new SQUID actions."""
    baseline = copy.deepcopy(cfg)
    baseline.pop("squid_jumps", None)
    return stage_configuration(baseline, "03_A_MaxFilt")


def accepted_report(cfg):
    """Validate accepted current or historical 4A evidence, read-only.

    New SQUID actions may intentionally differ from the accepted baseline. The
    fallback validates that original baseline without accepting the new actions;
    downstream require_review never uses this fallback for changed SSS inputs.
    """
    report_dir = Path(cfg["resolved_paths"]["derivatives_reports_dir"])
    try:
        require_review(cfg, "03_A_MaxFilt", report_dir)
    except RuntimeError:
        baseline = copy.deepcopy(cfg)
        baseline.pop("squid_jumps", None)
        require_review(baseline, "03_A_MaxFilt", report_dir)
    return json.loads((report_dir / "03_A_MaxFilt_review_report.json").read_text())


def read_channel_basis(cfg, records):
    """Validate saved channel-only evidence while a changed SSS awaits review.

    This is not an SSS acceptance bypass: only detector channel exclusions may
    use it. Downstream processing still checks the actual 4A report and receipt.
    """
    basis = json.loads((root(cfg) / "channel_basis.json").read_text())
    if basis.get("content_sha256") != digest({k:v for k,v in basis.items() if k != "content_sha256"}):
        raise RuntimeError("SQUID channel evidence changed/corrupt.")
    if basis["configuration"] != channel_basis_configuration(cfg):
        raise RuntimeError("SQUID channel configuration changed; reconcile accepted channels first.")
    if basis["sources"] != {r.source_filename: file_stamp(r.source_path) for r in records}:
        raise RuntimeError("SQUID channel source evidence changed.")
    if not basis.get("channel_reviewed"):
        raise RuntimeError("Accept fresh channels with 03_A_MaxFilt.py --interactive before SQUID review.")
    return basis


def remember_channel_basis(cfg, records, channels, *, reviewed, report=None):
    """Snapshot channel-only evidence without rewriting any acceptance record."""
    payload = {"configuration": channel_basis_configuration(cfg),
               "sources": {r.source_filename: file_stamp(r.source_path) for r in records},
               "channels": channels, "channel_reviewed": bool(reviewed)}
    if report is not None:
        payload["accepted_channel_report"] = report
    try:
        existing = read_channel_basis(cfg, records)
    except (FileNotFoundError, RuntimeError):
        existing = None
    if existing and existing["channels"] == channels and (report is None or existing.get("accepted_channel_report")):
        return existing  # Do not erase previous interactive channel acceptance.
    payload["content_sha256"] = digest(payload)
    directory = root(cfg); directory.mkdir(parents=True, exist_ok=True)
    immutable = directory / f"channel-basis-{payload['content_sha256']}.json"
    encoded = json.dumps(payload, indent=2, allow_nan=False)
    if not immutable.exists():
        with immutable.open("x") as handle:
            handle.write(encoded)
    # Atomically update only this SQUID channel pointer; immutable history and
    # the shared tuning/eye files are completely independent of this write.
    with tempfile.NamedTemporaryFile(mode="w", dir=directory, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(encoded)
    try:
        os.replace(temporary, directory / "channel_basis.json")
    finally:
        temporary.unlink(missing_ok=True)
    return payload


# -----------------------------------------------------------------------------
# Bounded two-pass detector; halos preserve candidates at chunk boundaries
# -----------------------------------------------------------------------------

def _segments(raw, channels, settings):
    """Yield core differences with filter halos, never bridging acquisition gaps."""
    from mne.annotations import _annotations_starts_stops
    sfreq = float(raw.info["sfreq"])
    width = max(1, int(round(settings["median_window_sec"] * sfreq)))
    width += 1 - width % 2
    halo = width // 2
    block = max(width + 1, int(round(settings["block_sec"] * sfreq)))
    starts, stops = _annotations_starts_stops(raw, ("edge", "bad_acq_skip"), invert=True)
    for seg_start, seg_stop in zip(starts, stops):
        for start in range(int(seg_start), int(seg_stop)-1, block):
            stop = min(start + block, int(seg_stop)-1)
            read_start, read_stop = max(int(seg_start), start-halo), min(int(seg_stop), stop+halo+1)
            # Read each on-disk FIF block once for all retained sensors; separate
            # channel reads would repeatedly decode the same interleaved data.
            data = raw.get_data(picks=channels, start=read_start, stop=read_stop)
            smooth = median_filter(data, size=(1, width), mode="nearest")
            offset = start-read_start
            yield start, np.diff(smooth, axis=1)[:,offset:offset+stop-start], np.diff(data, axis=1)[:,offset:offset+stop-start]


def detect(raw, excluded, settings):
    """Return channel-preserving event candidates and unscored-channel reasons.

    The first pass estimates each channel's derivative distribution using a
    deterministic bounded subsample. The second scores every native difference.
    Raw-derivative MAD floors the smoothed-derivative scale: median filtering
    can leave a small positive MAD on quantized/high-frequency recordings,
    which must not inflate ordinary fluctuations into artifacts. Channel and
    cross-channel groups are bounded by the merge gap rather than chained
    indefinitely. Flat/nonfinite channels are reported as unscored.
    """
    settings = parameters({"squid_jumps": settings})
    picks = [raw.ch_names[p] for p in mne.pick_types(raw.info, meg=True, exclude=[])
             if raw.ch_names[p] not in excluded]
    sfreq = float(raw.info["sfreq"])
    hits, unscored = [], []
    stride = max(1, int(np.ceil(raw.n_times / 200000)))
    sampled, unsmoothed = [], []
    if picks:
        for start, diff, raw_diff in _segments(raw, picks, settings):
            # Align sampling to the native recording grid. Copy only retained
            # statistics so small strided views cannot retain full read blocks.
            offset = (-start) % stride
            sampled.append(diff[:,offset::stride].copy())
            unsmoothed.append(raw_diff[:,offset::stride].copy())
    values_by_channel = np.concatenate(sampled, axis=1) if sampled else np.empty((len(picks),0))
    raw_by_channel = np.concatenate(unsmoothed, axis=1) if unsmoothed else np.empty((len(picks),0))
    del sampled, unsmoothed
    centers, scales = {}, {}
    for row, channel in enumerate(picks):
        values = values_by_channel[row]
        values = values[np.isfinite(values)]
        if not values.size:
            unscored.append({"channel": channel, "reason": "no_finite_differences"})
            continue
        center = float(np.median(values))
        scale = 1.4826 * float(np.median(np.abs(values-center)))
        fallback = raw_by_channel[row]
        fallback = fallback[np.isfinite(fallback)]
        raw_scale = 1.4826 * float(np.median(np.abs(fallback-np.median(fallback)))) if fallback.size else 0.
        # A merely nonzero filtered MAD is insufficient: the P04 pilot had
        # filtered/raw scale ratios near 0.06, giving widespread false flags.
        scale = max(scale, raw_scale)
        if not np.isfinite(scale) or scale <= 0:
            unscored.append({"channel": channel, "reason": "zero_robust_scale"})
            continue
        centers[channel], scales[channel] = center, scale
    del values_by_channel, raw_by_channel
    # Collapse within-channel crossings across block boundaries, retaining the
    # strongest edge and support before the final cross-channel merge.
    currents = {ch: None for ch in centers}
    gap = max(1, int(round(settings["merge_gap_sec"] * sfreq)))
    for start, differences, _ in (_segments(raw, picks, settings) if centers else []):
        for row, channel in enumerate(picks):
            if channel not in centers:
                continue
            current = currents[channel]
            scores = np.abs(differences[row]-centers[channel]) / scales[channel]
            indices = np.flatnonzero(np.isfinite(scores) & (scores >= settings["threshold"]))
            for index in indices:
                sample = int(raw.first_samp + start + index + 1)
                score = float(scores[index])
                if current is None or sample-current["first_sample"] > gap:
                    if current is not None:
                        hits.append(current)
                    current = {"channel": channel, "first_sample": sample, "last_sample": sample,
                               "peak_sample": sample, "score": score}
                else:
                    current["last_sample"] = sample
                    if score > current["score"]:
                        current.update(peak_sample=sample, score=score)
            currents[channel] = current
    for current in currents.values():
        if current is not None:
            hits.append(current)
    merged = []
    gap = int(round(settings["merge_gap_sec"] * sfreq))
    pad = int(round(settings["review_padding_sec"] * sfreq))
    # Group around an anchored peak, not the advancing end of the previous
    # group's support. Transitive cross-sensor chaining can otherwise merge
    # unrelated fluctuations across an entire recording.
    anchor = None
    for hit in sorted(hits, key=lambda row: row["peak_sample"]):
        if anchor is None or hit["peak_sample"]-anchor > gap:
            anchor = hit["peak_sample"]
            merged.append({"first_sample": hit["first_sample"], "last_sample": hit["last_sample"],
                           "peak_sample": hit["peak_sample"], "score": hit["score"], "channels": [hit["channel"]]})
        else:
            event = merged[-1]
            event["first_sample"] = min(event["first_sample"], hit["first_sample"])
            event["last_sample"] = max(event["last_sample"], hit["last_sample"])
            event["channels"] = sorted(set(event["channels"] + [hit["channel"]]))
            if hit["score"] > event["score"]:
                event.update(peak_sample=hit["peak_sample"], score=hit["score"])
    for i, event in enumerate(merged):
        event.update(id=i+1, start_sample=int(max(raw.first_samp, event["first_sample"]-pad)),
                     stop_sample=int(min(raw.first_samp+raw.n_times, event["last_sample"]+pad+1)),
                     time_sec=(event["peak_sample"]-raw.first_samp)/sfreq)
        # Baseline change is context only, never a gate that loses decaying jumps.
        peak = event["peak_sample"]-raw.first_samp
        radius = max(1, round(.05*sfreq))
        before = raw.get_data(picks=event["channels"], start=max(0, peak-radius), stop=peak)
        after = raw.get_data(picks=event["channels"], start=peak, stop=min(raw.n_times, peak+radius))
        changes = [float(np.median(after[j])-np.median(before[j])) for j in range(len(event["channels"]))]
        event["baseline_change"] = {ch: value if math.isfinite(value) else None
                                    for ch, value in zip(event["channels"], changes)}
    return {"screened_channels": picks, "excluded_channels": sorted(set(excluded)),
            "unscored_channels": unscored, "candidates": merged}


# -----------------------------------------------------------------------------
# Immutable proposals and explicit review receipts; no manual-eye file writes
# -----------------------------------------------------------------------------

def _identity(cfg, record, raw, excluded):
    """Bind one proposal to native source, accepted channels and detector settings."""
    return {"algorithm": ALGORITHM, "source": file_stamp(record.source_path),
            "coordinates": coordinates(raw), "excluded_channels": sorted(excluded),
            "channel_basis_configuration": channel_basis_configuration(cfg),
            "parameters": {k: v for k, v in parameters(cfg).items() if k not in ("enabled", "block_sec")}}


def proposal_path(cfg, record, identity):
    """Resolve an immutable, exact-source candidate filename outside raw data."""
    return root(cfg) / f"{record.source_stem}-{digest(identity)}-candidates.json"


def prepare(cfg, record, excluded):
    """Reuse or generate one immutable proposal; accepted science stays read-only."""
    raw = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
    identity = _identity(cfg, record, raw, excluded)
    path = proposal_path(cfg, record, identity)
    if path.is_file():
        payload = json.loads(path.read_text())
        if payload.get("identity") != identity or payload.get("content_sha256") != digest({k:v for k,v in payload.items() if k != "content_sha256"}):
            raise RuntimeError("SQUID proposal changed/corrupt; preserve it and prepare a fresh validated proposal.")
        return path, payload
    payload = {"identity": identity, "fingerprint": digest(identity), "source_filename": record.source_filename,
               "generated_utc": utc_now(), **detect(raw, excluded, parameters(cfg))}
    payload["content_sha256"] = digest(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation never replaces previous proposals or manual decisions.
    with path.open("x") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
    with path.with_suffix(".csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "time_sec", "peak_sample", "score", "channels", "start_sample", "stop_sample"])
        writer.writeheader()
        for row in payload["candidates"]:
            writer.writerow({key: row[key] for key in writer.fieldnames})
    return path, payload


def load_review(cfg, record, excluded):
    """Require matching source/settings/proposals and intact reviewed actions."""
    receipt = cfg.get("squid_jumps", {}).get("reviews", {}).get(record.source_filename)
    if not receipt or receipt.get("status") != "reviewed":
        raise RuntimeError(f"SQUID review required for {record.source_filename}; prepare/review before SSS.")
    raw = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
    identity = _identity(cfg, record, raw, excluded)
    if receipt.get("fingerprint") != digest(identity):
        raise RuntimeError(f"SQUID source/settings/channel basis changed: {record.source_filename}; reconcile review.")
    path = root(cfg) / receipt["review_file"]
    if Path(receipt["review_file"]).name != receipt["review_file"] or file_hash(path) != receipt["review_sha256"]:
        raise RuntimeError("SQUID review file changed or is not an exact-source local artifact.")
    saved = json.loads(path.read_text())
    if saved.get("application") != receipt.get("application") or saved.get("fingerprint") != receipt["fingerprint"]:
        raise RuntimeError("SQUID saved review and tuning receipt disagree.")
    candidate_file = proposal_path(cfg, record, identity)
    if file_hash(candidate_file) != saved.get("proposal_sha256"):
        raise RuntimeError("Reviewed SQUID proposals changed.")
    validate_actions(raw, saved["application"], set(excluded))
    return saved


def validate_actions(raw, action, excluded):
    """Validate exact native bounds and MEG additions without implicit scope."""
    if action.get("channel_scope") not in ("union", "run"):
        raise ValueError("SQUID channel scope must be union or explicitly run.")
    meg = {raw.ch_names[p] for p in mne.pick_types(raw.info, meg=True, exclude=[])}
    if not set(action.get("bad_channels", [])).issubset(meg-excluded):
        raise ValueError("New bad channels must be screened MEG channels.")
    for pair in action.get("intervals", []):
        if len(pair) != 2 or any(type(sample) is not int for sample in pair):
            raise ValueError("SQUID intervals require integer half-open native sample pairs.")
    attach_intervals(raw.copy(), action.get("intervals", []))


def save_review(cfg, record, excluded, proposal, proposal_file, action, *, reviewer, reason, dismissed):
    """Write immutable decisions and merge only this source under the YAML lock.

    The existing update_tuning_file reloads the latest map atomically, preserving
    overlapping P04/P05 eye-review saves. No saved BAD annotation is ever touched.
    """
    raw = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
    if _identity(cfg, record, raw, excluded) != proposal["identity"]:
        raise RuntimeError("Raw source or screening basis changed during review.")
    validate_actions(raw, action, set(excluded))
    if not reviewer.strip() or not reason.strip():
        raise ValueError("Explicit reviewer and review observations are required.")
    saved = {"status": "reviewed", "fingerprint": proposal["fingerprint"],
             "source_filename": record.source_filename, "application": action,
             "proposal_sha256": file_hash(proposal_file), "reviewed_utc": utc_now(),
             "reviewer": reviewer, "reason": reason, "dismissed_candidate_ids": sorted(dismissed),
             "screened_channels": proposal["screened_channels"], "excluded_channels": sorted(excluded),
             "unscored_channels": proposal["unscored_channels"]}
    path = root(cfg) / f"{record.source_stem}-{digest(saved)}-review.json"
    if not path.exists():
        with path.open("x") as handle:
            json.dump(saved, handle, indent=2, allow_nan=False)
    receipt = {"status": "reviewed", "fingerprint": saved["fingerprint"], "review_file": path.name,
               "review_sha256": file_hash(path), "application": action}
    update_tuning_file(path=cfg["squid_jumps"]["tuning_path"], raw_root=cfg["resolved_paths"]["data_root"],
        subject_id=cfg["subject"]["id"], section="squid_jumps", mapping_key="reviews",
        values={record.source_filename: receipt})
    return saved
