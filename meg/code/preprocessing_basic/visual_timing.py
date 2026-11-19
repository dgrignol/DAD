#!/usr/bin/env python3
"""Measure/propose Stage 3 cohort visual timing without changing Stage 1/2.

Usage from the DAD root (local paths are the pipeline defaults)::

    python meg/code/preprocessing_basic/visual_timing.py \
        --config meg/params/pipeline_default.yaml \
        --subject-config meg/params/sub-11.yaml --version s3-pilot-v001 --pilot

    python meg/code/preprocessing_basic/visual_timing.py \
        --config meg/params/pipeline_default.yaml \
        --site-config meg/params/pipeline_cluster.yaml \
        --subject-config meg/params/sub-11.yaml --version s3-v001

``--validate-only`` checks accepted inputs and FIF headers without loading any
signals or creating outputs. ``--pilot`` reads exactly the selected short
windows (at most 12 per supported participant P08/P09/P11/P24). Full cohort
jobs are user-submitted. ``--summarize-cohort`` requires all 36 matching full
proposal bundles; it never accepts their science. ``--resume`` verifies a
completed, matching immutable bundle. There is deliberately no overwrite flag.

Data flow: verified Stage 2 + canonical events + raw headers/saved Conf ->
short native STI101/diode windows -> edge/QC measurements -> restricted P11/P24
estimates and explicit digital fallbacks -> immutable versioned proposal.
Raw data, events/alignment, accepted ledgers and source manifests are inputs.
All output paths are guarded and live below derivatives/visual_timing/VERSION.
Each participant bundle includes an explained offline ``report.html`` with
categorical trigger-lag distributions, mean +/- SD summaries, donor evidence,
QC/sanity checks and embedded exact tables. Existing bundles can be presented
without signal rereads using ``visual_timing_report.py`` in a new report folder.
No cleaning, epochs, source-manifest refresh or numerical scientific acceptance
occurs here. A separate reviewed record is necessary for later epoching.
"""

from __future__ import annotations

import argparse
import io
import json
import re
import platform
import resource
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import mne
import numpy as np
import yaml
from scipy.io import loadmat

from downstream_stage_utils import csv_bool, load_trigger_event_rows
from pipeline_config import ANALYSIS_SUBJECT_IDS, load_and_resolve_config
from preprocessing_stage_utils import discover_raw_run_records, guarded_path, utc_now
from trial_preparation import (csv_text, digest, evidence_paths, file_hash,
                               load_prepared_decisions, read_csv, stage2_paths, validate_stage1)
from visual_timing_utils import (CONDITION_CODES, DEFAULT_METHOD, SCHEMA_VERSION,
                                 VISUAL_CODES, donor_allowed, estimate_lag, event_family,
                                 infer_polarity, measure_edge, round_half_away,
                                 verify_bundle, window_coverage)
from visual_timing_report import build_report_artifacts


DEFAULT_METHOD_PATH = Path(__file__).resolve().parents[2] / "params/stage3_visual_timing.yaml"
CONTEXT_FIELDS = ("experimentName", "refrate", "measuredIfiSec", "measuredRefreshHz",
                  "photodiodePatchEnabled", "photodiodePatchPositionPx", "photodiodePatchSizePx",
                  "photodiodeTrialPulseFrames", "photodiodeEventPulseFrames",
                  "photodiodeMarkedTriggerCodes", "scheduledOnsetGuardScope",
                  "scheduledOnsetGuardContinueAfterExhaustion", "ptbAppliedSkipSyncTests")
TRIAL_FIELDS = ("condition_trigger", "sequence_trigger_code", "trial_completed",
                "trial_onset_schedule_type", "trial_onset_boundary_trial",
                "trial_onset_guard_retry_count", "trial_onset_guard_exhausted",
                "trial_onset_guard_continued_after_exhaustion", "trial_onset_guard_final_lead_sec",
                "trial_onset_scheduling_lead_sec", "trial_onset_missed", "frame_missed_count",
                "dot_offset_presented", "dot_offset_missed", "dot_endpoint_hold_sec")


# -----------------------------------------------------------------------------
# Input fingerprints, exact source/Conf resolution and read-only validation
# -----------------------------------------------------------------------------


def load_method(path: str | Path) -> dict:
    """Require the documented detector contract; all values remain proposals."""
    method = yaml.safe_load(Path(path).read_text())
    if not isinstance(method, dict) or set(method) != set(DEFAULT_METHOD):
        raise ValueError("Stage 3 method config has missing/unknown fields.")
    if (method["threshold_fractions"] != [.1, .25, .5, .75, .9]
            or method["primary_fraction"] != .5 or method["approval"] != "proposal_only"
            or method["rounding"] != DEFAULT_METHOD["rounding"]
            or method["estimator"] != DEFAULT_METHOD["estimator"]):
        raise ValueError("Unsupported Stage 3 method/approval contract.")
    for key in ("first_response_peak_ms", "sustained_ms", "minimum_amplitude_v", "minimum_snr",
                "maximum_sensitivity_ms", "pulse_tolerance_ms", "minimum_donor_events_per_code",
                "bootstrap_resamples", "post_warmup_max_gap_ms"):
        if not isinstance(method[key], (int, float)) or not np.isfinite(method[key]) or method[key] <= 0:
            raise ValueError(f"Invalid Stage 3 method value: {key}")
    for key in ("window_ms", "baseline_ms", "search_ms", "post_warmup_trough_ms"):
        pair = method[key]
        if not isinstance(pair, list) or len(pair) != 2 or not np.isfinite(pair).all() or pair[0] >= pair[1]:
            raise ValueError(f"Invalid detector interval: {key}")
    for key in ("minimum_donor_events_per_code", "bootstrap_resamples", "random_seed"):
        if type(method[key]) is not int or method[key] < 0:
            raise ValueError(f"Integer method setting required: {key}")
    if method["minimum_donor_events_per_code"] < 2:
        raise ValueError("At least two donor events are needed for variability diagnostics.")
    if method["window_ms"][1] - method["window_ms"][0] > 250:
        raise ValueError("Stage 3 short-window method is bounded to 250 ms.")
    if not (method["window_ms"][0] <= method["baseline_ms"][0] < method["baseline_ms"][1] <= 0
            and method["window_ms"][0] <= method["search_ms"][0] <= 0
            and 80 <= method["search_ms"][1] < method["window_ms"][1]):
        raise ValueError("Detector baseline/search must fit the native waveform window.")
    if not 0 <= method["post_warmup_trough_ms"][0] < method["post_warmup_trough_ms"][1] <= 30:
        raise ValueError("Post-warm-up low trough must fit the documented 0:30-ms context window.")
    return method


def version_root(cfg, version):
    """Resolve a safe, non-traversing immutable version directory."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", version) or ".." in version:
        raise ValueError("Version must be a simple non-traversing identifier.")
    root = Path(cfg["resolved_paths"]["derivatives_root"]) / "visual_timing" / version
    guarded_path(root, cfg["resolved_paths"]["data_root"])
    return root


def authority_path(cfg):
    """Read the single shared decisions record, never create a workspace copy."""
    p = Path(cfg["project"]["docs_root"]) / "PREPROCESSING_DECISIONS.md"
    if not p.is_file():
        p = Path(cfg["resolved_paths"]["data_root"]).parent / "DAD/meg/docs/PREPROCESSING_DECISIONS.md"
    text = p.read_text()
    if "Accepted on 2026-10-19" not in text or "Stage 2 trial-accounting gate outcome" not in text:
        raise RuntimeError("Authoritative accepted Stage 1/2 record unavailable.")
    return p


def read_contexts(cfg):
    """Resolve saved runtime/patch settings by exact saved behavioral block.

    Saved MAT files are small read-only evidence. Their hashes and Conf values
    are recorded; filenames and experimentName must agree with each profile's
    runtime_by_block, not today's launcher or a schedule/stimulus version.
    """
    contexts, hashes = {}, {}
    folder = Path(cfg["resolved_paths"]["subject_data_dir"]) / "output_files"
    subject = cfg["subject"]["id"]
    expected = cfg["subject_profile"]
    for block in expected["completed_behavioral_blocks"]:
        runtime = expected.get("runtime_by_block", {}).get(str(block),
                  expected.get("runtime_by_block", {}).get(block, expected["versions"]["runtime"]))
        matches = sorted(folder.glob(f"MoveDot1_occlusion_*_SUB{subject}_BLOCK{block:02d}.mat"))
        if len(matches) != 1 or not (f"_{runtime}_" in matches[0].name or
                                    runtime == "v20" and "_v20_rescueTraject_" in matches[0].name):
            raise RuntimeError(f"Ambiguous/wrong saved runtime for {subject} block {block}: {matches}")
        p = matches[0]
        saved = loadmat(p, variable_names=["Conf", "runResults"], squeeze_me=True, struct_as_record=False)
        c = saved.get("Conf")
        if c is None or runtime not in str(getattr(c, "experimentName", "")):
            raise RuntimeError(f"Saved Conf runtime mismatch: {p}")
        values = {}
        for key in CONTEXT_FIELDS:
            value = getattr(c, key, None)
            values[key] = value.tolist() if isinstance(value, np.ndarray) else value.item() if isinstance(value, np.generic) else value
        contexts[str(block)] = dict(values, runtime=runtime, source_mat=p.name,
                                    stimuli=expected["versions"]["stimuli"], schedule=expected["versions"]["schedule"])
        # Join saved context by exact block/run/trial/replay identity, never by
        # array position. PTB/guard fields remain explanatory metadata only.
        trials = {}
        for run in np.atleast_1d(saved.get("runResults", [])):
            for trial in np.atleast_1d(getattr(run, "output", [])):
                if not hasattr(trial, "trial_num"):
                    continue
                replay = int(bool(getattr(trial, "is_replay_trial", False)))
                trial_key = f"{int(run.run_num)}:{int(trial.trial_num)}:{replay}"
                if trial_key in trials:
                    raise RuntimeError(f"Ambiguous saved trial diagnostics: {p.name}/{trial_key}")
                diagnostics = {}
                for key in TRIAL_FIELDS:
                    v = getattr(trial, key, None)
                    if isinstance(v, np.ndarray):
                        v = v.item() if v.size == 1 else None
                    if isinstance(v, np.generic):
                        v = v.item()
                    if isinstance(v, float) and not np.isfinite(v):
                        v = None
                    diagnostics[key] = v
                trials[trial_key] = diagnostics
        contexts[str(block)]["trial_diagnostics"] = trials
        hashes[p.name] = file_hash(p)
    return contexts, hashes


def validate_inputs(cfg):
    """Verify accepted tables, policies, mandatory exclusions and actual headers.

    No signals or derivatives are read/created by the header portion. The
    original Stage 2 report remains authoritative for the immutable accounting
    snapshot; later decisions-document edits do not invalidate its old hash.
    """
    if cfg["subject"]["id"] not in ANALYSIS_SUBJECT_IDS:
        raise ValueError("Stage 3 cohort is P04-P39.")
    authority = authority_path(cfg)
    inputs = validate_stage1(cfg, evidence_paths(cfg))
    ledger = load_prepared_decisions(cfg)
    keys = [r["trial_key"] for r in ledger]
    if len(keys) != len(set(keys)):
        raise RuntimeError("Duplicate original Stage 2 trial keys.")
    for accepted in cfg["epoching"].get("accepted_trial_exclusions", []):
        rows = [r for r in ledger if r["source_fif"] == accepted["source_fif"]
                and r["raw_sample_frame"] == accepted["raw_sample_frame"]]
        if (len(rows) != 1 or rows[0]["included_by_trial_policy"] or
                rows[0]["condition_code"] != accepted["condition_code"] or
                rows[0]["accepted_exclusion_reason"] != accepted["reason"]):
            raise RuntimeError("Accepted source/sample exclusion changed.")
    events = load_trigger_event_rows(evidence_paths(cfg)["events"])
    sources = discover_raw_run_records(cfg)
    expected_sources = {r["source_fif"] for r in ledger}
    if expected_sources != {s.source_filename for s in sources} or expected_sources != {e["source_fif"] for e in events}:
        raise RuntimeError("Stage 3 source set differs from accepted trial/event coverage.")
    coverage = {r["source_fif"]: r for r in read_csv(stage2_paths(cfg)["coverage"])}
    headers = []
    for source in sources:
        saved = coverage[source.source_filename]
        with mne.io.read_raw_fif(source.source_path, preload=False, verbose="ERROR") as raw:
            if (raw.first_samp != int(saved["first_sample_raw"]) or
                    raw.first_samp + raw.n_times != int(saved["end_sample_raw_exclusive"]) or
                    raw.info["sfreq"] != float(saved["sfreq_hz"]) or
                    Path(source.source_path).stat().st_size != int(saved["source_size_bytes"])):
                raise RuntimeError(f"Raw acquisition coordinates/size changed: {source.source_filename}")
            headers.append(dict(saved, diode_channel_present="MISC008" in raw.ch_names))
    contexts, mat_hashes = read_contexts(cfg)
    for event in events:
        block = event.get("behavioral_block")
        if block and str(block) not in contexts:
            raise RuntimeError("An event has no resolved saved-runtime context.")
    signature = {"stage2_report_sha256": file_hash(stage2_paths(cfg)["report"]),
                 "stage2_output_sha256": json.loads(stage2_paths(cfg)["report"].read_text())["output_sha256"],
                 "saved_mat_sha256": mat_hashes, "decisions_record_sha256": file_hash(authority)}
    return ledger, events, sources, headers, contexts, signature


# -----------------------------------------------------------------------------
# Bounded pilot selection and native-channel measurement on original anchors
# -----------------------------------------------------------------------------


def event_key(subject, event):
    """Keep event identity distinct from its owning condition-trial identity."""
    return f"sub-{subject}:{event['source_fif']}:{event['raw_sample_frame']}:{event['trigger_code']}"


def select_pilot(events, subject):
    """Select <=12 fixed windows per participant before looking at diode data.

    P09 emphasizes acquisition 2 with controls. P11 includes confirmed unusable
    early blocks plus all three donor blocks. P24 reads acquisition 5 only.
    Other pilot participants are forbidden, keeping four tasks <=48 windows.
    """
    allocation = {"08": {1: 6, 2: 6}, "09": {1: 3, 2: 6, 3: 3},
                  "11": {1: 1, 2: 1, 3: 4, 4: 3, 5: 3}, "24": {5: 12}}
    if subject not in allocation:
        raise ValueError("Signal pilot is bounded to P08/P09/P11/P24; use --validate-only for other cases.")
    selected = []
    for run, count in allocation[subject].items():
        rows = sorted([e for e in events if int(e["acquisition_run"]) == run
                       and int(e["trigger_code"]) in VISUAL_CODES], key=lambda r: int(r["raw_sample_frame"]))
        chosen = []
        conditions = [e for e in rows if int(e["trigger_code"]) in CONDITION_CODES]
        # Three positions per source, then code coverage and diagnostic families.
        candidates = ([conditions[i] for i in sorted(set([0, len(conditions)//2, len(conditions)-1]))]
                      if conditions else [])
        for code in (71, 81, 91, 153, 111, 118, 112, 114, 115, 116, 117):
            group = [e for e in rows if int(e["trigger_code"]) == code]
            if group:
                candidates.append(group[0])
        candidates += [rows[i] for i in np.linspace(0, max(0, len(rows)-1), count, dtype=int)] if rows else []
        for e in candidates + rows:
            if e not in chosen:
                chosen.append(e)
            if len(chosen) == count:
                break
        selected.extend(chosen)
    return selected


def confirmed_unusable(subject, event):
    """Confirmed facts constrain assignments, irrespective of channel presence."""
    return subject == "11" and int(event.get("behavioral_block") or 0) in (1, 2)


def measure_events(cfg, events, sources, contexts, method, *, pilot, cached_windows=None):
    """Read only two short native channel windows per selected observed event.

    Polarity controls come from these same selected windows, so pilot inference
    performs no hidden extra reads. No production edge can leave its source
    recording. Confirmed-unusable P11 early windows are diagnostics only.
    """
    subject = cfg["subject"]["id"]
    selected = select_pilot(events, subject) if pilot else [e for e in events if int(e["trigger_code"]) in VISUAL_CODES]
    by_source = defaultdict(list)
    for e in selected:
        by_source[e["source_fif"]].append(e)
    measurements, waveforms, reads = [], {}, 0
    history = {}
    for source in sources:
        previous = None
        for event in sorted([e for e in events if e["source_fif"] == source.source_filename
                             and int(e["trigger_code"]) in VISUAL_CODES], key=lambda e:int(e["raw_sample_frame"])):
            history[event_key(subject,event)] = previous
            previous = event
    for source in sources:
        rows = by_source[source.source_filename]
        if not rows:
            continue
        with mne.io.read_raw_fif(source.source_path, preload=False, verbose="ERROR") as raw:
            sfreq = float(raw.info["sfreq"])
            before, after = [round_half_away(ms * sfreq / 1000) for ms in method["window_ms"]]
            windows = {}
            for e in rows:
                key = event_key(subject, e)
                center = int(e["raw_sample_frame"]) - int(raw.first_samp)
                if "MISC008" not in raw.ch_names or confirmed_unusable(subject, e) and not pilot:
                    continue
                start, stop = max(0, center + before), min(raw.n_times, center + after)
                if cached_windows is None:
                    values = raw.get_data(picks=["STI101", "MISC008"], start=start, stop=stop)
                    reads += 1
                else:
                    if key not in cached_windows:
                        raise RuntimeError(f"Pilot cache is missing a selected native window: {key}")
                    values, cached_anchor, cached_sfreq = cached_windows[key]
                    if cached_anchor != center-start or cached_sfreq != sfreq or values.shape != (2,stop-start):
                        raise RuntimeError("Cached native window coordinates differ from accepted source.")
                windows[key] = (values, center - start)
            controls = [(windows[event_key(subject, e)][0][1], windows[event_key(subject, e)][1])
                        for e in rows if int(e["trigger_code"]) in CONDITION_CODES
                        and event_key(subject, e) in windows and not confirmed_unusable(subject, e)]
            polarity, polarity_evidence = infer_polarity(controls, sfreq, method)
            for e in rows:
                key, code = event_key(subject, e), int(e["trigger_code"])
                context = contexts.get(str(e.get("behavioral_block", "")), {})
                frames = int(context.get("photodiodeTrialPulseFrames") or 4) if code in CONDITION_CODES | {153} else int(context.get("photodiodeEventPulseFrames") or 2)
                refresh = float(context.get("measuredRefreshHz") or context.get("refrate") or 120)
                row = {"event_key": key, "subject_bids_id": f"sub-{subject}",
                       "source_fif": e["source_fif"], "event_index": e["event_index"],
                       "digital_sample_raw": int(e["raw_sample_frame"]), "trigger_code": code,
                       "acquisition_run": int(e["acquisition_run"]), "sfreq_hz": sfreq,
                       "first_samp": int(raw.first_samp), "family": event_family(code),
                       "behavioral_block": e.get("behavioral_block", ""),
                       "behavioral_run": e.get("behavioral_run", ""),
                       "behavioral_trial": e.get("behavioral_trial", ""),
                       "behavioral_frame": e.get("behavioral_frame", ""),
                       "reported_source_csv": Path(e.get("reported_source_csv", "")).name,
                       "is_replay_trial": e.get("is_replay_trial", ""),
                       "runtime": context.get("runtime", "unassigned_control"),
                       "patch_position_px": json.dumps(context.get("photodiodePatchPositionPx")),
                       "measurement_status": "not_recorded", "qc_reasons": "diode_channel_absent",
                       "polarity_evidence": json.dumps(polarity_evidence, sort_keys=True),
                       "signed_lag_definition": "diode_edge_minus_original_STI101"}
                trial_id = f"{e.get('behavioral_run','')}:{e.get('behavioral_trial','')}:{int(csv_bool(e.get('is_replay_trial',False)))}"
                diagnostics = context.get("trial_diagnostics", {}).get(trial_id, {})
                row.update({f"saved_{k}": v for k, v in diagnostics.items()})
                row["saved_diagnostics_join"] = "exact_block_run_trial_replay" if diagnostics else "not_available"
                if code in CONDITION_CODES and diagnostics.get("condition_trigger") not in (None, code):
                    raise RuntimeError(f"Saved condition and accepted event disagree: {key}")
                previous = history.get(key)
                gap = ((row["digital_sample_raw"]-int(previous["raw_sample_frame"]))*1000/sfreq) if previous else None
                post_warmup = bool(code in CONDITION_CODES and previous and int(previous["trigger_code"])==153
                                   and 0 < gap <= method["post_warmup_max_gap_ms"])
                row.update(preceding_visual_code=int(previous["trigger_code"]) if previous else None,
                           preceding_visual_gap_ms=gap, post_warmup_context=post_warmup)
                if key in windows:
                    values, anchor = windows[key]
                    waveforms[key] = (values, anchor, sfreq)
                    if confirmed_unusable(subject, e):
                        row.update(measurement_status="confirmed_unusable", qc_reasons="P11_blocks_1_2_confirmed_unusable")
                    elif polarity is None:
                        row.update(measurement_status="unusable", qc_reasons="ambiguous_recording_polarity")
                    else:
                        row.update(measure_edge(values[1], sfreq, anchor, polarity, frames, refresh, method,
                                                warmup=code == 153, post_warmup=post_warmup))
                        if row["edge_offset_samples"] is not None:
                            row["diode_edge_sample_raw"] = row["digital_sample_raw"] + row["edge_offset_samples"]
                    row["anchor_recorded_STI101_value"] = int(round(values[0, anchor]))
                elif confirmed_unusable(subject, e):
                    row.update(measurement_status="confirmed_unusable", qc_reasons="P11_blocks_1_2_confirmed_unusable")
                measurements.append(row)
    if pilot and reads > 12:
        raise RuntimeError("Pilot exceeded the fixed 12-window participant bound.")
    return measurements, waveforms, reads


def event_proposals(subject, measurements, estimates, method):
    """Add separate event assignments; warm-up remains a diagnostic control.

    P11/P24 non-condition estimates follow the same own-participant/code-only
    donor restrictions. No extra signal reads occur. Missing/flagged edges in
    other participants stay unresolved rather than using a global correction.
    """
    proposals = []
    for event in measurements:
        row = dict(event, assignment_basis="unresolved", proposed_event_sample_raw=None,
                   applied_shift_samples=None, lag_estimate_ms=None, scientific_approval="pending")
        if event["trigger_code"] == 153:
            row["assignment_basis"] = "diagnostic_only"
        elif event["measurement_status"] == "measured":
            row.update(assignment_basis="measured", proposed_event_sample_raw=event["diode_edge_sample_raw"],
                       applied_shift_samples=event["edge_offset_samples"], lag_estimate_ms=event["lag_ms"])
        elif subject in {"04", "05", "06", "07", "17"}:
            row.update(assignment_basis="digital_fallback", proposed_event_sample_raw=event["digital_sample_raw"],
                       applied_shift_samples=0, uncertainty_kind="physical_lag_unknown")
        else:
            target = dict(event, condition_code=event["trigger_code"])
            donors = [d for d in measurements if donor_allowed(subject, target, d)]
            if donors:
                eid = f"P{subject}_code{event['trigger_code']}"
                if eid not in estimates:
                    estimates[eid] = estimate_lag(donors, method)
                estimate = estimates[eid]
                row["donor_estimate_id"] = eid
                if estimate["status"] == "estimated":
                    shift = round_half_away(estimate["lag_estimate_ms"] * event["sfreq_hz"] / 1000)
                    row.update(assignment_basis="estimated", proposed_event_sample_raw=event["digital_sample_raw"]+shift,
                               applied_shift_samples=shift, lag_estimate_ms=estimate["lag_estimate_ms"],
                               uncertainty_kind="donor_location_ci_and_observed_percentiles;transfer_unknown")
        proposals.append(row)
    return proposals


# -----------------------------------------------------------------------------
# Restricted timing proposals, native coverage and recording-level diagnostics
# -----------------------------------------------------------------------------


def build_map(cfg, ledger, events, measurements, method):
    """Retain all trial keys; never resurrect an excluded/unrecorded candidate."""
    subject = cfg["subject"]["id"]
    measured = {r["event_key"]: r for r in measurements}
    maps, estimates = [], {}
    for original in ledger:
        row = {k: original[k] for k in ("trial_key", "subject_bids_id", "source_fif", "raw_sample_frame",
               "condition_code", "acquisition_run", "behavioral_block", "behavioral_run", "sfreq_hz",
               "recording_first_sample_raw", "recording_end_sample_raw_exclusive", "included_by_trial_policy",
               "exclusion_reasons", "accepted_exclusion_reason", "trial_interval_end_raw_exclusive")}
        row.update(online_highpass_hz=original.get("online_highpass_hz"),
                   online_lowpass_hz=original.get("online_lowpass_hz"))
        row.update(assignment_status="review_required", assignment_basis="unresolved",
                   proposed_sample_raw=None, applied_shift_samples=None, lag_estimate_ms=None,
                   uncertainty_kind="unknown", donor_estimate_id="", scientific_approval="pending")
        key = f"sub-{subject}:{original['source_fif']}:{original['raw_sample_frame']}:{original['condition_code']}"
        edge = measured.get(key)
        uncertainty = None
        if edge and edge["measurement_status"] == "measured":
            shift = edge["edge_offset_samples"]
            lag = edge["lag_ms"]
            row.update(assignment_basis="measured", assignment_status="proposed", applied_shift_samples=shift,
                       proposed_sample_raw=int(original["raw_sample_frame"]) + shift, lag_estimate_ms=lag,
                       uncertainty_kind=edge["uncertainty_kind"], measurement_event_key=key)
            lo = edge["edge_10_offset_samples"] * 1000 / float(original["sfreq_hz"])
            hi = edge["edge_90_offset_samples"] * 1000 / float(original["sfreq_hz"])
            uncertainty = (lo - 1000 / float(original["sfreq_hz"]), hi)
        elif subject in {"04", "05", "06", "07", "17"}:
            row.update(assignment_basis="digital_fallback", assignment_status="proposed", applied_shift_samples=0,
                       proposed_sample_raw=int(original["raw_sample_frame"]),
                       uncertainty_kind="physical_lag_unknown;zero_applied_shift_is_not_zero_latency")
        elif subject == "11" and int(original["behavioral_block"]) in (1, 2) or subject == "24" and int(original["acquisition_run"]) < 5:
            eid = f"P{subject}_code{original['condition_code']}"
            if eid not in estimates:
                donors = [d for d in measurements if donor_allowed(subject, original, d)]
                estimates[eid] = estimate_lag(donors, method)
                estimates[eid]["target_context_limitation"] = (
                    "P11 patch position changes from [193,0] to [245,208];spatial_transfer_unknown"
                    if subject == "11" else "P24 earlier acquisition settings differ;run5_only;cross_run_transfer_unknown")
            estimate = estimates[eid]
            row["donor_estimate_id"] = eid
            if estimate["status"] == "estimated":
                lag = estimate["lag_estimate_ms"]
                samples = lag * float(original["sfreq_hz"]) / 1000
                shift = round_half_away(samples)
                row.update(assignment_basis="estimated", assignment_status="proposed", applied_shift_samples=shift,
                           proposed_sample_raw=int(original["raw_sample_frame"]) + shift, lag_estimate_ms=lag,
                           rounding_residual_samples=shift-samples,
                           uncertainty_kind="donor_location_ci_and_observed_percentiles;transfer_unknown")
                uncertainty = tuple(estimate["donor_observed_p05_p95_ms"])
        if row["proposed_sample_raw"] is not None:
            row.update(window_coverage(original, row["proposed_sample_raw"], cfg["epoching"]["tmin_sec"],
                                       cfg["epoching"]["tmax_sec"], uncertainty))
            row["proposed_recording_time_sec"] = (row["proposed_sample_raw"] - int(original["recording_first_sample_raw"])) / float(original["sfreq_hz"])
        row["original_recording_time_sec"] = (int(original["raw_sample_frame"]) - int(original["recording_first_sample_raw"])) / float(original["sfreq_hz"])
        maps.append(row)
    return maps, estimates


def recording_summaries(headers, measurements):
    """Summarize each code separately, retaining tails and diagnostic trends."""
    summaries = []
    for header in headers:
        rows = [r for r in measurements if r["source_fif"] == header["source_fif"]]
        for code in sorted({r["trigger_code"] for r in rows}) or [None]:
            group = [r for r in rows if r["trigger_code"] == code]
            usable = [r for r in group if r["measurement_status"] == "measured"]
            values = np.asarray([r["lag_ms"] for r in usable], float)
            x = np.asarray([(r["digital_sample_raw"] - r["first_samp"]) / r["sfreq_hz"] for r in usable])
            row = {"source_fif": header["source_fif"], "acquisition_run": header["acquisition_run"],
                   "diode_channel_present": header["diode_channel_present"], "trigger_code": code,
                   "n_assessed": len(group), "n_measured": len(usable),
                   "measurement_status_counts": json.dumps(dict(Counter(r["measurement_status"] for r in group))),
                   "usability": "not_recorded" if not header["diode_channel_present"] else
                       "review_required" if not usable else "measured_subset_requires_review",
                   "n_flagged": sum(r["measurement_status"] == "flagged" for r in group)}
            row.update({k: header[k] for k in ("first_sample_raw", "end_sample_raw_exclusive",
                       "sfreq_hz", "online_highpass_hz", "online_lowpass_hz")})
            if len(values):
                row.update(mean_lag_ms=float(np.mean(values)), median_lag_ms=float(np.median(values)),
                           sd_lag_ms=float(np.std(values, ddof=1)) if len(values)>1 else None,
                           min_lag_ms=float(min(values)), max_lag_ms=float(max(values)),
                           p05_lag_ms=float(np.quantile(values,.05)), p95_lag_ms=float(np.quantile(values,.95)),
                           scaled_mad_ms=float(1.4826*np.median(np.abs(values-np.median(values)))))
                if len(values) >= 3 and np.ptp(x):
                    row["linear_drift_ms_per_hour_diagnostic_only"] = float(np.polyfit(x-x.mean(), values, 1)[0]*3600)
            amplitudes=[float(r["amplitude_v"]) for r in group if r.get("amplitude_v") is not None]
            if amplitudes:
                row.update(min_amplitude_v=min(amplitudes),median_amplitude_v=float(np.median(amplitudes)),
                           max_amplitude_v=max(amplitudes))
            summaries.append(row)
    return summaries


def physical_intervals(measurements, contexts):
    """Compare within-trial diode-only intervals with saved frame differences.

    Both endpoints must be measured in the same exact source/context. Expected
    frame counts are context, never replacements for missing physical markers.
    There is no inferred optical dot-offset endpoint.
    """
    groups = defaultdict(list)
    for row in measurements:
        if row["measurement_status"] == "measured" and row["trigger_code"] != 153:
            key = tuple(str(row.get(k, "")) for k in ("source_fif", "behavioral_block", "behavioral_run",
                  "behavioral_trial", "is_replay_trial", "reported_source_csv"))
            groups[key].append(row)
    output = []
    for group in groups.values():
        starts = [r for r in group if r["trigger_code"] in CONDITION_CODES]
        pairs = [(starts[0], r) for r in group if r["trigger_code"] not in CONDITION_CODES] if len(starts) == 1 else []
        for a, b in ((111, 112), (112, 114), (114, 115), (116, 117)):
            left, right = [r for r in group if r["trigger_code"] == a], [r for r in group if r["trigger_code"] == b]
            if len(left) == len(right) == 1:
                pairs.append((left[0], right[0]))
        for left, right in pairs:
            if not left["behavioral_frame"] or not right["behavioral_frame"]:
                continue
            refresh = float(contexts[str(left["behavioral_block"])]["refrate"])
            observed = (right["diode_edge_sample_raw"]-left["diode_edge_sample_raw"])*1000/left["sfreq_hz"]
            expected = (int(right["behavioral_frame"])-int(left["behavioral_frame"]))*1000/refresh
            # Question frames belong to a separate UI loop, not the motion
            # frame axis; retain measured spacing without a false expectation.
            comparable = right["trigger_code"] != 118
            output.append({"left_event_key":left["event_key"],"right_event_key":right["event_key"],
                           "observed_diode_interval_ms":observed,
                           "expected_frame_interval_ms":expected if comparable else None,
                           "interval_error_ms":observed-expected if comparable else None,
                           "expectation_basis":"saved_motion_frames" if comparable else "question_UI_frame_axis_not_comparable"})
    return output


# -----------------------------------------------------------------------------
# Guarded immutable publication, plots and cohort aggregation
# -----------------------------------------------------------------------------


def publish(folder, contents, report, raw_root):
    """Reserve a new directory exclusively, then write completion metadata last."""
    guarded_path(folder, raw_root)
    folder = Path(folder)
    folder.parent.mkdir(parents=True, exist_ok=True)
    folder.mkdir(exist_ok=False)
    report["output_sha256"] = {}
    for name, content in contents.items():
        if Path(name).name != name:
            raise ValueError("Bundle output must be a direct child file.")
        path = folder / name
        guarded_path(path, raw_root)
        with path.open("xb") as f:
            f.write(content if isinstance(content, bytes) else content.encode())
        report["output_sha256"][name] = file_hash(path)
    with (folder / "manifest.json").open("x") as f:
        json.dump(report, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def diagnostic_images(measurements, waveforms):
    """Show actual native traces plus lag/sensitivity trends and visible failures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [r for r in measurements if r["event_key"] in waveforms]
    # Pilot displays every read. Full proposals show flags, tails and controls.
    order = sorted(rows, key=lambda r: (r["measurement_status"] == "measured", r.get("lag_ms") or -999))
    selected = order[:6]
    for i in np.linspace(0, max(0,len(rows)-1), 6, dtype=int) if rows else []:
        if rows[i] not in selected:
            selected.append(rows[i])
    if len(rows) <= 12:
        # Cache every bounded pilot read for independent/offline inspection;
        # repeated detector checks need not read participant FIFs again.
        selected = rows
    nrows=max(1,int(np.ceil(len(selected)/3)))
    fig, axes = plt.subplots(nrows, 3, figsize=(15, 3*nrows), squeeze=False)
    arrays, index = {}, []
    for ax, row in zip(axes.flat, selected):
        values, anchor, sf = waveforms[row["event_key"]]
        t = (np.arange(values.shape[1])-anchor)*1000/sf
        ax.plot(t, values[1]*1000, lw=1, label="raw diode mV")
        ax.axvline(0, color="k", ls="--", lw=.8)
        if row.get("lag_ms") is not None:
            ax.axvline(row["lag_ms"], color="red", ls=":")
        ax.set_title(f"{row['source_fif']}\n{row['trigger_code']} {row['measurement_status']} {row.get('qc_reasons','')}", fontsize=7)
        ax.set_xlabel("ms relative to original STI101")
        arrays[f"window_{len(index):03d}"] = np.vstack([t,values])
        index.append({"array": f"window_{len(index):03d}", "event_key": row["event_key"]})
    for ax in list(axes.flat)[len(selected):]:
        ax.axis("off")
    fig.suptitle("Stage 3 proposal QC: original raw diode, digital anchor, candidate 50% edge")
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf,format="png",dpi=130); plt.close(fig)
    out = io.BytesIO(); np.savez_compressed(out,**arrays)
    fig, axes = plt.subplots(2,2,figsize=(12,8))
    valid=[r for r in measurements if r.get("lag_ms") is not None]
    codes = sorted({r["trigger_code"] for r in valid})
    for position, code in enumerate(codes):
        group=[r for r in valid if r["trigger_code"]==code]
        axes[0,0].scatter([(r["digital_sample_raw"]-r["first_samp"])/r["sfreq_hz"] for r in group],
                          [r["lag_ms"] for r in group],s=8,label=str(code))
        # Codes are factors: equal categorical spacing, never numeric distance.
        usable = [r["lag_ms"] for r in group if r["measurement_status"] == "measured"]
        if usable:
            mean = float(np.mean(usable))
            axes[0,1].bar(position, mean, width=.65)
            if len(usable)>1:
                axes[0,1].errorbar(position,mean,yerr=float(np.std(usable,ddof=1)),fmt="none",ecolor="#222222",capsize=3)
            axes[0,1].annotate(f"n={len(usable)}",(position,mean),xytext=(0,6),textcoords="offset points",ha="center",fontsize=7)
    axes[0,0].set(xlabel="time within source recording (s)",ylabel="signed lag (ms)")
    if valid: axes[0,0].legend(fontsize=6,ncol=4)
    axes[0,1].set(xticks=range(len(codes)),xticklabels=[str(c) for c in codes],
                  xlabel="Trigger code (categorical)",ylabel="Measured lag mean +/- sample SD (ms)")
    axes[1,0].scatter([r["lag_ms"] for r in valid],[r["threshold_span_ms"] for r in valid],s=8)
    axes[1,0].set(xlabel="lag (ms)",ylabel="10-90% threshold span (ms)")
    counts=Counter(r["measurement_status"] for r in measurements)
    axes[1,1].bar(list(counts),list(counts.values()));axes[1,1].tick_params(axis="x",labelrotation=25)
    axes[1,1].set_ylabel("assessed events (including failures)")
    fig.tight_layout();summary=io.BytesIO();fig.savefig(summary,format="png",dpi=130);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for source in sorted({r["source_fif"] for r in measurements}):
        group=[r for r in measurements if r["source_fif"]==source and r.get("amplitude_v") is not None]
        x=[(r["digital_sample_raw"]-r["first_samp"])/r["sfreq_hz"] for r in group]
        axes[0].scatter(x,[r["amplitude_v"]*1000 for r in group],s=10,label=source)
        axes[1].scatter(x,[r.get("baseline_noise_v",0)*1000 for r in group],s=10)
        axes[2].scatter(x,[r.get("threshold_span_ms",np.nan) for r in group],s=10)
    axes[0].set_ylabel("response amplitude (mV)")
    axes[1].set_ylabel("baseline noise (mV)")
    axes[2].set(xlabel="time within each source recording (s); sources overlaid",ylabel="10-90% span (ms)")
    if measurements:axes[0].legend(fontsize=7)
    fig.suptitle("Signal usability over recording time: inspect attenuation and dropout")
    fig.tight_layout();signal=io.BytesIO();fig.savefig(signal,format="png",dpi=130);plt.close(fig)
    return {"diagnostics.png":buf.getvalue(),"lag_sensitivity.png":summary.getvalue(),
            "signal_usability.png":signal.getvalue(),
            "waveform_windows.npz":out.getvalue(),"diagnostic_index.csv":csv_text(index)}


def run_subject(cfg, version, method, *, pilot=False, validate_only=False, resume=False, pilot_cache_from=None):
    """Produce one proposal or bounded pilot; accepted evidence remains read-only."""
    start=time.monotonic()
    folder=version_root(cfg,version)/cfg["layout"]["subject_bids_id"]/cfg["layout"]["session_bids_id"]
    ledger,events,sources,headers,contexts,signature=validate_inputs(cfg)
    implementation={p.name:file_hash(p) for p in [Path(__file__),Path(__file__).with_name("visual_timing_utils.py"),
                                                Path(__file__).with_name("visual_timing_report.py")]}
    contract=dict(signature,method=method,implementation_sha256=implementation)
    cached_windows, cache_provenance = None, None
    if pilot_cache_from:
        if not pilot or pilot_cache_from == version:
            raise ValueError("--pilot-cache-from requires pilot mode and a different source version.")
        previous_folder=version_root(cfg,pilot_cache_from)/cfg["layout"]["subject_bids_id"]/cfg["layout"]["session_bids_id"]
        previous=verify_bundle(previous_folder)
        if previous["mode"]!="pilot" or previous["method"]["window_ms"]!=method["window_ms"] or any(previous[k]!=v for k,v in signature.items()):
            raise RuntimeError("Pilot cache has different accepted inputs or native waveform window.")
        cached_windows={}
        with np.load(previous_folder/"waveform_windows.npz",allow_pickle=False) as cache:
            for entry in read_csv(previous_folder/"diagnostic_index.csv"):
                values=cache[entry["array"]]
                anchor=int(np.flatnonzero(np.isclose(values[0],0))[0])
                sfreq=float(1000/(values[0,1]-values[0,0]))
                cached_windows[entry["event_key"]]=(values[1:].copy(),anchor,sfreq)
        cache_provenance={"version":pilot_cache_from,"manifest_sha256":file_hash(previous_folder/"manifest.json")}
        contract["pilot_cache"]=cache_provenance
    fingerprint=digest(contract)
    if folder.exists():
        previous=verify_bundle(folder)
        if not resume or previous["input_contract_sha256"]!=fingerprint or previous["mode"]!=("pilot" if pilot else "cohort_proposal"):
            raise RuntimeError("Existing immutable Stage 3 bundle; choose a new version or matching --resume.")
        return previous
    report={"schema_version":SCHEMA_VERSION,"version":version,"subject_bids_id":cfg["layout"]["subject_bids_id"],
            "status":"validated" if validate_only else "completed","mode":"pilot" if pilot else "cohort_proposal",
            "scientific_acceptance":"pending_independent_review","generated_utc":utc_now(),
            "input_contract_sha256":fingerprint,"input_contract":contract,**signature,
            "method":method,"implementation_sha256":implementation,
            "epoch_window_sec":[cfg["epoching"]["tmin_sec"],cfg["epoching"]["tmax_sec"]],
            "n_recordings":len(headers),"n_original_trial_keys":len(ledger)}
    import scipy
    report["software_versions"]={"python":platform.python_version(),"mne":mne.__version__,
                                 "numpy":np.__version__,"scipy":scipy.__version__,"pyyaml":yaml.__version__}
    report["pilot_cache"]=cache_provenance
    if validate_only:
        return report
    measurements,waveforms,reads=measure_events(cfg,events,sources,contexts,method,pilot=pilot,cached_windows=cached_windows)
    maps,estimates=build_map(cfg,ledger,events,measurements,method)
    proposals=event_proposals(cfg["subject"]["id"],measurements,estimates,method)
    # Detect a concurrent upstream change before publishing any new proposal.
    load_prepared_decisions(cfg)
    if file_hash(stage2_paths(cfg)["report"])!=signature["stage2_report_sha256"]:
        raise RuntimeError("Stage 2 report changed during timing measurement.")
    contexts_now,hashes_now=read_contexts(cfg)
    if hashes_now!=signature["saved_mat_sha256"] or contexts_now!=contexts:
        raise RuntimeError("Saved runtime context changed during timing measurement.")
    report.update(n_signal_windows_read=reads,n_assessed_visual_events=len(measurements),
                  measurement_status_counts=dict(Counter(r["measurement_status"] for r in measurements)),
                  assignment_counts=dict(Counter(r["assignment_basis"] for r in maps)),
                  eligible_assignment_counts=dict(Counter(r["assignment_basis"] for r in maps if r["included_by_trial_policy"])),
                  elapsed_seconds=time.monotonic()-start)
    report["n_native_windows_assessed"]=len(waveforms)
    report["peak_process_rss_mib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(2**20 if sys.platform=="darwin" else 1024)
    recording_rows = recording_summaries(headers,measurements)
    intervals = physical_intervals(measurements,contexts)
    # Reporting consumes observations/proposals; it cannot change the detector,
    # donor estimates, stable trial keys, eligibility or assigned sample shifts.
    artifacts = build_report_artifacts(measurements,maps,estimates,recording_rows,intervals,
        contexts,report,waveforms,source_identity={"bundle_manifest":"manifest.json",
            "raw_signal_windows_read":reads,"native_windows_assessed":len(waveforms)},proposals=proposals)
    contents={"event_measurements.csv":csv_text(measurements,("event_key","measurement_status")),
              "event_timing_proposals.csv":csv_text(proposals,("event_key","assignment_basis")),
              "trial_timing_map.csv":csv_text(maps),"recording_summary.csv":csv_text(recording_rows),
              "contexts.json":json.dumps(contexts,indent=2,sort_keys=True),
              "estimates.json":json.dumps(estimates,indent=2,sort_keys=True,allow_nan=False),
              "physical_intervals.csv":csv_text(intervals,
                  ("left_event_key","right_event_key","observed_diode_interval_ms","expected_frame_interval_ms","interval_error_ms")),
              }
    # Authoritative scientific CSV bytes win over any presentation companions.
    contents={**artifacts,**contents}
    # Freeze the actual scientific source alongside its hashes. This remains
    # reproducible even when the workspace has uncommitted development changes.
    for name,expected in implementation.items():
        source=Path(__file__).with_name(name)
        if file_hash(source)!=expected:
            raise RuntimeError("Timing implementation changed while the proposal was running.")
        contents[f"code_{name}.txt"]=source.read_bytes()
    contents["method.yaml"]=yaml.safe_dump(method,sort_keys=True)
    report["elapsed_seconds"]=time.monotonic()-start
    report["peak_process_rss_mib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(2**20 if sys.platform=="darwin" else 1024)
    publish(folder,contents,report,cfg["resolved_paths"]["data_root"])
    return report


def summarize_cohort(cfg,version,method,*,validate_only=False):
    """Require every participant's intact full proposal and matching method/code.

    Reload current Stage 2 evidence for each participant and require matching
    original trial counts. No signals are read and no scientific acceptance is
    inferred from a successful summary job.
    """
    root=version_root(cfg,version)
    reports=[]
    for subject in ANALYSIS_SUBJECT_IDS:
        current,_,_,_=load_and_resolve_config(cfg["_config_sources"]["default"],
            str(Path(cfg["project"]["params_root"])/f"sub-{subject}.yaml"),
            site_config_path=cfg["_config_sources"].get("site"),required_sections=("epoching",))
        ledger=load_prepared_decisions(current)
        folder=root/f"sub-{subject}"/current["layout"]["session_bids_id"]
        report=verify_bundle(folder)
        if (report["mode"]!="cohort_proposal" or report["method"]!=method or report["version"]!=version
            or report["stage2_report_sha256"]!=file_hash(stage2_paths(current)["report"])
            or report["n_original_trial_keys"]!=len(ledger)
            or reports and (report["implementation_sha256"]!=reports[0]["implementation_sha256"]
                           or report.get("software_versions")!=reports[0].get("software_versions"))):
            raise RuntimeError(f"Incomplete/stale/mismatched cohort proposal: {subject}")
        reports.append(report)
    result={"schema_version":SCHEMA_VERSION,"status":"completed","version":version,
            "scientific_acceptance":"pending_independent_review","n_participants":len(reports),
            "n_original_trial_keys":sum(r["n_original_trial_keys"] for r in reports),
            "participants":reports,"generated_utc":utc_now()}
    if not validate_only:
        lines=["# Stage 3 cohort visual timing","","P04-P39; scientific acceptance pending independent review.","",
               "Open each participant's self-contained HTML report for explained figures, actual proposed shifts, QC/donors and the independent-review checklist.","",
               "| Participant report | Recordings | Original keys | Eligible assignment counts |","| --- | ---: | ---: | --- |"]
        lines.extend(f"| [{r['subject_bids_id']}](../{r['subject_bids_id']}/ses-01/report.html) | {r['n_recordings']} | {r['n_original_trial_keys']} | {json.dumps(r['eligible_assignment_counts'],sort_keys=True)} |" for r in reports)
        publish(root/"group",{"cohort_report.md":"\n".join(lines)+"\n"},result,cfg["resolved_paths"]["data_root"])
    return result


def main():
    """Resolve configuration without invoking the runner or creating inventories."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",required=True)
    parser.add_argument("--site-config")
    parser.add_argument("--subject-config")
    parser.add_argument("--method-config",default=str(DEFAULT_METHOD_PATH))
    parser.add_argument("--version",required=True)
    parser.add_argument("--pilot",action="store_true")
    parser.add_argument("--pilot-cache-from",help="Reassess an intact bounded pilot cache; no participant signal rereads.")
    parser.add_argument("--validate-only",action="store_true")
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--summarize-cohort",action="store_true")
    args=parser.parse_args()
    if args.summarize_cohort and (args.pilot or args.resume or args.pilot_cache_from):
        parser.error("Cohort summary does not support pilot/resume.")
    if not args.summarize_cohort and not args.subject_config:
        parser.error("Participant operation requires --subject-config.")
    cfg,warnings,_,_=load_and_resolve_config(args.config,args.subject_config,
        site_config_path=args.site_config,required_sections=("epoching",))
    if warnings:
        raise RuntimeError(f"Configuration warnings require investigation: {warnings}")
    method=load_method(args.method_config)
    result=(summarize_cohort(cfg,args.version,method,validate_only=args.validate_only) if args.summarize_cohort else
            run_subject(cfg,args.version,method,pilot=args.pilot,validate_only=args.validate_only,resume=args.resume,
                        pilot_cache_from=args.pilot_cache_from))
    print(json.dumps({k:result[k] for k in ("status","version","scientific_acceptance","n_original_trial_keys")},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
