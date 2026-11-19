#!/usr/bin/env python3
"""Prepare and apply P24's explicit EDF auxiliary-channel reconstruction.

Usage::

    PYTHONPATH=/path/to/eyelinkio/site python p24_auxiliary_reconstruction.py \
      --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
      --evidence-log /Volumes/MORWUR/Projects/DAMIANO/Attention/DAD/cluster_logs/DAD_epoch_export.11793_24.log \
      --output /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/auxiliary_reconstruction/p24-edf-20261117-v1 \
      --workers 4

The builder reads original-clock EDF samples, accepted trial-local mapping
sidecars and the failed run's frozen native-input contract. It creates six
EDF-native-value arrays for all five acquisitions, a validity mask and a manifest.
MISC order is left pupil/x/y, right pupil/x/y, verified against run-5 recordings.
Coordinates are pixels and pupil sizes are EDF area values; these are NOT
claimed to be reconstructed analog voltages. Run 5's acquired MISC001-006 are
replaced only in the new in-memory variant, making units consistent across runs.
The original acquired run-5 FIF remains immutable and available for comparison.

Only uniquely matched, recording-on segments receive nearest measured EDF
samples (maximum half-sample error). Invalid eyes, missing samples, recording
gaps, unmatched trials and times outside coverage receive explicit zeros, with
the separate validity mask distinguishing zero fill from a measured zero.
There is no interpolation across gaps, clock-rate fit, physiological alignment,
new BAD annotation or MEG timing change. This is an auxiliary reconstruction,
not new evidence for changing accepted annotations or visual-onset timing.

The runtime adapter is opt-in via epoching.auxiliary_reconstruction_manifest,
a derivative-root-relative path. It is restricted to labelled inline_no_ica
P24 runs, validates every bound input before processing, supplies MISC001-006
from EDF in every run and omits MISC008 in run 5 in memory. It checks every retained sample,
annotation and acquisition setting, then reorders all five runs to 325 common
channels. Source FIFs, EDFs, masks, models and old processing runs are immutable.
No epoching, export, deployment or scheduler submission is performed by this
CLI. Independent acquisitions use separate processes and disjoint outputs;
publish the shared manifest only after all workers finish successfully.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import ctypes
import json
import multiprocessing
from pathlib import Path

import mne
import numpy as np

from continuous_cleaning import annotation_rows, coordinates, digest, file_hash, file_stamp

CHANNELS = [f"MISC00{i}" for i in range(1, 7)]
FIELDS = ["pupil_left", "gaze_x_left", "gaze_y_left", "pupil_right", "gaze_x_right", "gaze_y_right"]
UNITS = ["EDF pupil area", "pixels", "pixels", "EDF pupil area", "pixels", "pixels"]
SCHEMA = "p24-edf-auxiliary-v1"


# -----------------------------------------------------------------------------
# Decode original timestamps and map only reviewed trial/recording segments.
# -----------------------------------------------------------------------------

def read_six_fields(path):
    """Read six valid-eye EDF fields without compressing recording-off gaps."""
    from edf_blink_candidates import reader_api
    api, defs = reader_api()
    error = ctypes.c_int()
    handle = api.edf_open_file(str(Path(path).resolve()).encode(), 2, 0, 1, ctypes.byref(error))
    if not handle or error.value:
        if handle: api.edf_close_file(handle)
        raise RuntimeError(f"Cannot read original EDF: {path}")
    rows = []
    try:
        while True:
            kind = api.edf_get_next_data(handle)
            if kind == defs.event_constants["NO_PENDING_ITEMS"]: break
            if kind != 200: continue
            sample = api.edf_get_float_data(handle).contents.fs
            row = [(sample.time + (.5 if sample.flags & defs.SAMPLE_ADD_OFFSET else 0.)) / 1000.]
            for eye, flag in [(0, defs.SAMPLE_LEFT), (1, defs.SAMPLE_RIGHT)]:
                pupil = float(sample.pa[eye]) if sample.flags & flag and sample.flags & defs.SAMPLE_PUPILSIZE else np.nan
                valid = np.isfinite(pupil) and 0 < pupil < 1e8
                row.append(pupil if valid else np.nan)
                for field in (sample.gx, sample.gy):
                    value = float(field[eye]) if valid and sample.flags & defs.SAMPLE_GAZEXY else np.nan
                    row.append(value if -32768 < value < 1e8 else np.nan)
            rows.append(row)
    finally:
        api.edf_close_file(handle)
    data = np.asarray(rows, dtype=float).reshape(-1, 7)
    if not len(data) or np.any(np.diff(data[:, 0]) <= 0):
        raise RuntimeError("EDF sample timestamps must be present and strictly increasing")
    return data


def map_samples(samples, segments, *, first_samp, n_times, sfreq):
    """Place measured samples on native MEG coordinates; preserve every gap."""
    data = np.zeros((6, n_times), dtype=np.float32)
    valid = np.zeros((6, n_times), dtype=bool)
    assigned = np.zeros(n_times, dtype=bool)
    for segment in segments:
        if segment["status"] != "matched" or segment["stop"] is None:
            raise RuntimeError("Unreviewed or open EDF segment cannot be reconstructed")
        bounds = [int(segment["meg_anchor_sample"] + np.rint((t-segment["edf_anchor_sec"])*sfreq)) - first_samp
                  for t in (segment["start"], segment["stop"])]
        indices = np.arange(max(0, bounds[0]), min(n_times, bounds[1]))
        times = segment["edf_anchor_sec"] + (indices+first_samp-segment["meg_anchor_sample"])/sfreq
        right = np.clip(np.searchsorted(samples[:, 0], times), 0, len(samples)-1)
        left = np.maximum(0, right-1)
        nearest = np.where(abs(samples[left, 0]-times) < abs(samples[right, 0]-times), left, right)
        measured = (abs(samples[nearest, 0]-times) <= .51/segment["sfreq"])
        measured &= (samples[nearest, 0] >= segment["start"]) & (samples[nearest, 0] < segment["stop"])
        indices, nearest = indices[measured], nearest[measured]
        if assigned[indices].any(): raise RuntimeError("Mapped EDF recording segments overlap")
        values = samples[nearest, 1:].T
        good = np.isfinite(values)
        data[:, indices] = np.where(good, values, 0.)
        valid[:, indices] = good
        assigned[indices] = True
    return data, valid


def _build_one(task):
    """Write one disjoint native-grid auxiliary array and return its binding."""
    native, mapping, edf_path, output = task
    print(f"P24 {native['source_filename']}: reconstructing EDF auxiliary samples", flush=True)
    samples = read_six_fields(edf_path)
    coord = native["coordinates"]
    data, valid = map_samples(samples, mapping["segments"], first_samp=coord["first_samp"],
        n_times=coord["n_times"], sfreq=coord["sfreq"])
    target = Path(output) / (Path(native["source_filename"]).stem + "_edf_aux.npz")
    with target.open("xb") as handle:
        np.savez_compressed(handle, data=data, valid=valid)
    coverage = [float(v.mean()) for v in valid]
    print(f"P24 {native['source_filename']}: saved; valid fractions {coverage}", flush=True)
    return {"source_filename": native["source_filename"], "array": target.name,
        "array_sha256": file_hash(target), "valid_fraction_by_channel": dict(zip(CHANNELS, coverage)),
        "zero_fill_samples_by_channel": dict(zip(CHANNELS, [int((~v).sum()) for v in valid])),
        "matched_segments": len(mapping["segments"]), "fields": FIELDS, "units": UNITS}


# -----------------------------------------------------------------------------
# Bind the explicit opt-in plan, without changing existing run contracts.
# -----------------------------------------------------------------------------

def plan_path(cfg):
    """Resolve only an explicit derivative-root-relative reconstruction plan."""
    value = cfg["epoching"].get("auxiliary_reconstruction_manifest")
    if value is None: return None
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise RuntimeError("Auxiliary reconstruction manifest must be relative to derivatives_root")
    return Path(cfg["resolved_paths"]["derivatives_root"]) / relative


def validate_plan(cfg, inputs):
    """Validate all plan/array/source bindings before any participant processing."""
    path = plan_path(cfg)
    if path is None: return None
    if (str(cfg["subject"]["id"]).zfill(2) != "24" or cfg["epoching"].get("preparation_mode") != "inline_no_ica"
            or not cfg["epoching"].get("processing_run")):
        raise RuntimeError("EDF auxiliary reconstruction requires an explicitly labelled P24 no-ICA variant")
    plan = json.loads(path.read_text())
    if plan.get("schema") != SCHEMA or plan.get("native_inputs") != inputs:
        raise RuntimeError("P24 reconstruction source coordinates, annotations or file bindings differ")
    if plan["channels"] != CHANNELS or plan["removed_channels"] != ["MISC008"]:
        raise RuntimeError("Unexpected P24 auxiliary channel policy")
    if plan.get("implementation_sha256") != file_hash(__file__):
        raise RuntimeError("P24 auxiliary reconstruction implementation changed")
    for row in plan["reconstructed_runs"]:
        name = Path(row["array"])
        if name.name != str(name) or file_hash(path.parent/name) != row["array_sha256"]:
            raise RuntimeError("P24 auxiliary reconstruction array changed")
    return {"manifest_sha256": file_hash(path), "definition": plan}


def apply_plan(raw, cfg, source_filename, evidence):
    """Apply the authorized auxiliary change in memory; verify retained samples."""
    path = plan_path(cfg)
    if path is None: return None
    if file_hash(path) != evidence["manifest_sha256"]:
        raise RuntimeError("Auxiliary plan changed after participant validation")
    plan = evidence["definition"]
    before_annotations = annotation_rows(raw)
    before_settings = {k: raw.info[k] for k in ("sfreq", "highpass", "lowpass", "bads", "meas_date")}
    retained = [n for n in raw.ch_names if n not in CHANNELS + ["MISC008"]]
    # Digest all retained samples, not only MEG, to verify that replacement of
    # the explicitly authorized eye channels never changes a trigger or sensor.
    import hashlib
    def retained_hash():
        h = hashlib.sha256()
        for a in range(0, raw.n_times, 10000):
            h.update(np.ascontiguousarray(raw.get_data(picks=retained, start=a, stop=a+10000)).tobytes())
        return h.hexdigest()
    before_hash = retained_hash()
    row = next((r for r in plan["reconstructed_runs"] if r["source_filename"] == source_filename), None)
    omitted = [n for n in raw.ch_names if n in CHANNELS + ["MISC008"]]
    if row is not None:
        if any(n in raw.ch_names for n in CHANNELS) and source_filename != plan["reference_acquired_auxiliary_source"]:
            raise RuntimeError("Unexpected acquired eye channels outside P24 run 5")
        target = path.parent/row["array"]
        if file_hash(target) != row["array_sha256"]: raise RuntimeError("Auxiliary array changed")
        with np.load(target, allow_pickle=False) as payload:
            data, valid = payload["data"], payload["valid"]
        if data.shape != (6, raw.n_times) or valid.shape != data.shape or not np.isfinite(data).all() or np.any(data[~valid] != 0):
            raise RuntimeError("Invalid auxiliary values, shape or zero-fill mask")
        if omitted: raw.drop_channels(omitted)
        auxiliary = mne.io.RawArray(data.astype(float), mne.create_info(CHANNELS, raw.info["sfreq"], "misc"),
            first_samp=raw.first_samp, verbose="ERROR")
        auxiliary.set_meas_date(raw.info["meas_date"])
        raw.add_channels([auxiliary], force_update_info=True)
    else:
        raise RuntimeError("P24 source is absent from the auxiliary policy")
    if "MISC008" in raw.ch_names: raw.drop_channels(["MISC008"])
    raw.reorder_channels(plan["final_channel_order"])
    if retained_hash() != before_hash or annotation_rows(raw) != before_annotations:
        raise RuntimeError("A retained sample or annotation changed during auxiliary reconstruction")
    if any(raw.info[k] != v for k,v in before_settings.items()):
        raise RuntimeError("Acquisition settings changed during auxiliary reconstruction")
    return {"status": "applied_in_memory", "retained_samples_verified": True,
        "coordinates_after": coordinates(raw), "EDF_supplied_channels": CHANNELS,
        "replaced_or_omitted_acquired_channels": omitted,
        "channel_units": dict(zip(CHANNELS, UNITS)),
        "zero_fill_samples_by_channel": row["zero_fill_samples_by_channel"],
        "manifest_sha256": evidence["manifest_sha256"]}


# -----------------------------------------------------------------------------
# Validate the entire input table before launching independent builders.
# -----------------------------------------------------------------------------

def main():
    """Create immutable auxiliary arrays; never invoke epoching or a scheduler."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derivatives-root", type=Path, required=True)
    parser.add_argument("--evidence-log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1: parser.error("--workers must be positive")
    decoder = json.JSONDecoder(); text = args.evidence_log.read_text(); reports = []
    for offset, char in enumerate(text):
        if char != "{": continue
        try: obj, _ = decoder.raw_decode(text[offset:])
        except ValueError: continue
        if isinstance(obj, dict) and obj.get("stage") == "epoching" and obj.get("subject_bids_id") == "sub-24" and "preparation_contract" in obj: reports.append(obj)
    if len(reports) != 1: raise RuntimeError("Require one exact P24 epoch input contract in the evidence log")
    native = reports[0]["preparation_contract"]["definition"]["preparation"]["inputs"]
    if len(native) != 5: raise RuntimeError("Require the five accepted P24 acquisitions")
    maps = []; tasks = []
    events = args.derivatives_root/"preprocessing/sub-24/ses-01/events/02_trigger_events.csv"
    for item in native:
        source = item["source_filename"]
        folder = args.derivatives_root/"preprocessing/sub-24/ses-01"
        for part in item["parts"]:
            if file_stamp(folder/"sss"/part["name"]) != part: raise RuntimeError("Native FIF changed since epoch preflight")
        mapping_path = folder/"edf_blink_candidates"/(Path(source).stem+"_mapping.json")
        mapping = json.loads(mapping_path.read_text())
        if mapping["coordinates"] != item["coordinates"] or mapping["events_sha256"] != file_hash(events):
            raise RuntimeError("EDF mapping coordinates or canonical events changed")
        if len(mapping["edf_inputs"]) != 1: raise RuntimeError("Require one uniquely associated EDF per acquisition")
        edf = mapping["edf_inputs"][0]
        edf_path = args.derivatives_root.parent/"data"/edf["relative_path"]
        if file_stamp(edf_path) != edf["stamp"]: raise RuntimeError("Original EDF changed")
        maps.append({"source_filename": source, "mapping_sha256": file_hash(mapping_path), "edf": edf,
            "unmatched_trial_count": sum(t["status"] != "matched" for t in mapping["trials"])})
        tasks.append((item, mapping, str(edf_path), str(args.output)))
    final_channels = [n for n in native[-1]["coordinates"]["channels"] if n != "MISC008"]
    if len(final_channels) != 325 or any(set(x["coordinates"]["channels"]) | set(CHANNELS) != set(final_channels) for x in native[:-1]):
        raise RuntimeError("Unexpected native P24 channel inventories")
    if args.output.resolve().is_relative_to((args.derivatives_root.parent/"data").resolve()):
        raise RuntimeError("Refusing to write inside raw data")
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        for future in as_completed([pool.submit(_build_one, t) for t in tasks]): rows.append(future.result())
    rows.sort(key=lambda r:r["source_filename"])
    manifest = {"schema": SCHEMA, "subject": "sub-24", "units_policy": "EDF native pupil-area and pixel units in all five runs; original acquired voltages remain in source FIFs",
        "implementation_sha256": file_hash(__file__), "native_inputs": native,
        "channels": CHANNELS, "removed_channels": ["MISC008"], "final_channel_order": final_channels,
        "reference_acquired_auxiliary_source": native[-1]["source_filename"], "mapping_inputs": maps,
        "canonical_events_sha256": file_hash(events), "reconstructed_runs": rows,
        "missing_policy": "zero with separate validity mask; no interpolation across unrecorded/unmatched intervals"}
    with (args.output/"manifest.json").open("x") as handle: json.dump(manifest, handle, indent=2)
    print(f"P24 auxiliary inputs prepared: {args.output/'manifest.json'}", flush=True)


if __name__ == "__main__": main()
