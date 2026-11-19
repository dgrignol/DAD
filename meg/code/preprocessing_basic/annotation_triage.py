#!/usr/bin/env python3
"""Estimate annotation-only epoch retention without accepting cleaning outputs.

Usage (from the DAD root)
-------------------------
    python meg/code/preprocessing_basic/annotation_triage.py --subjects 06 07 --jobs 2
    python meg/code/preprocessing_basic/annotation_triage.py --subjects 04 05 06 \
      --label first-retention
    # Cluster: add --site-config meg/params/pipeline_cluster.yaml.
    # First collect choices with review_ocular_components.py. Missing choices
    # remain explicit incomplete coverage; they never mean no ocular artifacts.

Execution: independent participants run in spawned processes, defaulting to the
allocated SLURM CPU count (or local CPUs). --jobs can reduce concurrency. Each
worker writes only its own masks and logs run-level progress. Cohort tables and
the summary plot are published by the parent after all workers finish.

Inputs: accepted 4A SSS, report-bound provisional muscle evidence, exact-source
manual reviews, compatible frontal ICA caches/selections, and accepted Stage 2/3
eligibility/timing. Reviewed manual intervals take precedence over proposals;
reviewed muscle replacements include empty sets. EDF and SQUID candidates are
never read or promoted. Existing source BAD intervals remain operative.

Data flow: validated source annotations + optional ocular detectors -> immutable
provisional JSON masks -> corrected epoch-window overlap accounting -> trial,
condition, identity and cohort reports. Only new files under
derivatives/annotation_triage are written. No signal processing outside the
existing cached-source detectors, epoch creation, amplitude rejection, ICA
component removal, acceptance receipt or production derivative is performed.
The masks may seed explicit targeted manual review; component selection alone
is not whole-source interval review. Counts are conditional on the current
annotations, not a claim of clean data. Short/NaN/disabled detection coverage is
reported separately. Epoch endpoints and annotation overlap follow MNE's native
sample grid, including nonzero first_samp and inclusive final epoch samples.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
from collections import Counter
from pathlib import Path
from uuid import uuid4

import mne
import numpy as np

from blink_annotation import (AUTO_LABEL, EYE_AUTO_LABEL, EYE_METHOD, blink_settings,
    detect_blinks, detect_eye_movements, eye_movement_settings, finite_stretches,
    load_cached_aid)
from continuous_cleaning import (annotation_rows, coordinates, digest, file_hash,
    file_stamp, require_review, reviewed_muscle_override)
from pipeline_config import load_and_resolve_config
from participant_parallel import add_parallel_argument, map_participants, participant_workers
from preprocessing_stage_utils import (derivative_path_for_run, discover_raw_run_records,
    guarded_path, stage_output_dirs, utc_now, write_json, write_text)
from trial_preparation import load_prepared_decisions
from visual_timing_utils import analysis_eligible, load_reviewed_timing

SCHEMA = "annotation-triage-v1"
PURPOSE = "preliminary_annotation_retention_only"


# -----------------------------------------------------------------------------
# Shared discovery and protected exact-source inputs; never run the pipeline
# -----------------------------------------------------------------------------

def subject_config(args, subject):
    """Resolve current subject tuning read-only, including saves made since startup."""
    selector = Path(args.subject_config_dir) / f"sub-{int(subject):02}.yaml"
    cfg, *_ = load_and_resolve_config(args.config, selector,
        site_config_path=args.site_config, required_sections=("ica", "artifact_annotation", "epoching"))
    tuning = Path(cfg["resolved_paths"].get("upstream_derivatives_reports_dir",
        cfg["resolved_paths"]["derivatives_reports_dir"])) / "03_tuning_choices.yaml"
    if not tuning.is_file():
        raise RuntimeError(f"Missing accepted participant tuning: {tuning}")
    cfg, *_ = load_and_resolve_config(args.config, selector, tuning,
        site_config_path=args.site_config, required_sections=("ica", "artifact_annotation", "epoching"))
    return cfg


def triage_root(cfg):
    """Use a dedicated namespace; no caller chooses a production output path."""
    base = Path(cfg["resolved_paths"]["derivatives_root"]).resolve()
    root = guarded_path(base / "annotation_triage", cfg["resolved_paths"]["data_root"])
    if root != base / "annotation_triage":
        raise RuntimeError("The triage namespace must not redirect through a symlink.")
    return root


def triage_path(cfg, path):
    """Reject descendant symlinks escaping the dedicated output namespace."""
    resolved = guarded_path(path, cfg["resolved_paths"]["data_root"])
    if not resolved.is_relative_to(triage_root(cfg)):
        raise RuntimeError("Triage output must remain inside derivatives/annotation_triage.")
    return resolved


def selection_path(cfg, record):
    """One resumable component-choice sidecar, separate from scientific tuning."""
    return triage_path(cfg, triage_root(cfg) / "component_selections" / cfg["layout"]["subject_bids_id"] / f"{record.source_stem}.json")


def has_manual_work(cfg, record):
    """Also protect orphan/stale saves: an invalid review is never replaced by auto."""
    from manual_eye_annotations import manual_path
    reviews = cfg["artifact_annotation"].get("manual_eye_reviews", {})
    return record.source_filename in reviews or manual_path(cfg, record).exists()


def load_context(cfg, record):
    """Compose current annotations in memory with existing validation/precedence.

    Return an unpreloaded source header, its path and a portable mask binding.
    Manual files/receipts must agree, even when the older 4B report predates SAVE.
    Legacy ocular-only reviews retain automatic muscle; explicit reviewed muscle
    sets replace it. Existing SSS/source BADs and acquisition origins are retained.
    """
    from manual_eye_annotations import load_manual, load_muscle_context
    dirs = stage_output_dirs(cfg)
    source = derivative_path_for_run(record, dirs["sss"], "_raw_sss.fif")
    raw = mne.io.read_raw_fif(source, preload=False, verbose="ERROR")
    muscle, basis = load_muscle_context(cfg, record, raw, source, dirs)
    if muscle is None:
        raise RuntimeError("Missing source-bound muscle evidence; prepare the affected 4B evidence first.")
    raw.annotations.delete([i for i, a in enumerate(raw.annotations) if a["description"].lower() == "bad_muscle"])
    raw.annotations.append(muscle.onset, muscle.duration, muscle.description)
    receipt = None
    if has_manual_work(cfg, record):
        manual, receipt = load_manual(cfg, record, raw, source)
        if reviewed_muscle_override(receipt, cfg["artifact_annotation"], basis["threshold"]):
            raw.annotations.delete([i for i, a in enumerate(raw.annotations) if a["description"].lower() == "bad_muscle"])
        seen = {(a["onset"], a["duration"], a["description"]) for a in raw.annotations}
        for a in manual:
            key = (a["onset"], a["duration"], a["description"])
            if key not in seen:
                raw.annotations.append(*key)
                seen.add(key)
    binding = {"subject_bids_id": cfg["layout"]["subject_bids_id"],
        "source_filename": record.source_filename, "acquisition_run": record.acquisition_run_label,
        "sss_input": file_stamp(source), "coordinates": coordinates(raw),
        "annotation_sha256": digest(annotation_rows(raw)), "muscle_basis": basis,
        "manual_receipt_sha256": digest(receipt) if receipt is not None else None}
    return raw, source, binding


def component_binding(record, cache):
    """Bind human indices to the exact source, cached model and source projection."""
    meta = cache["metadata"]
    return {"source_filename": record.source_filename, "acquisition_run": record.acquisition_run_label,
            "cache_fingerprint": meta["fingerprint"], "model_sha256": meta["model_sha256"],
            "sources_stamp": meta["sources_stamp"]}


def validate_selection(selection, record, cache):
    """Reject stale/wrong-source/out-of-range choices, including bool-as-int values."""
    if (selection.get("schema") != SCHEMA or selection.get("purpose") != "ocular_component_selection_only"
            or selection.get("status") != "selected" or selection.get("binding") != component_binding(record, cache)):
        raise RuntimeError("Component selection is stale or belongs to another source/model; revisit this recording.")
    for key in ("blink_components", "eye_movement_components"):
        values = selection.get(key)
        if (not isinstance(values, list) or any(type(i) is not int or not 0 <= i < cache["model"].n_components_ for i in values)
                or values != sorted(set(values))):
            raise ValueError(f"Invalid {key}; use displayed zero-based components or explicit NONE.")
    return selection


# -----------------------------------------------------------------------------
# Provisional detectors and immutable masks, independent of manual SAVE receipts
# -----------------------------------------------------------------------------

def mask_rows(raw):
    """Serialize rejecting intervals in acquisition seconds; filter edges are informational."""
    return [{"onset": float(a["onset"]), "duration": float(a["duration"]), "description": str(a["description"])}
            for a in raw.annotations if a["description"].lower().startswith("bad")
            and a["description"].lower() != "bad_filter_edge"]


def cached_mask(cfg, record, input_sha256):
    """Reuse an intact detector result only after all current inputs were validated."""
    folder = triage_root(cfg) / "proposals" / cfg["layout"]["subject_bids_id"] / record.source_stem
    for path in sorted(folder.glob("*.json")):
        saved = json.loads(path.read_text())
        if saved.get("input_sha256") == input_sha256:
            if digest({k: v for k, v in saved.items() if k != "payload_sha256"}) != saved.get("payload_sha256"):
                raise RuntimeError(f"Cached triage mask content changed: {path}")
            return saved
    return None


def provisional_mask(cfg, record):
    """Preserve manual work, or run existing ocular detectors on selected cached sources.

    Missing/deferred selections retain known BADs but carry incomplete coverage.
    Selected NONE disables that detector, never establishes artifact absence.
    No selection is loaded at all for an already manually reviewed recording.
    """
    raw, _, binding = load_context(cfg, record)
    result = {"schema": SCHEMA, "purpose": PURPOSE, "input_binding": binding,
              "method": "preserved_manual" if binding["manual_receipt_sha256"] else "provisional_automatic",
              "detectors": {}, "selection_sha256": None,
              "settings": {"blink": blink_settings(cfg), "eye_movement": eye_movement_settings(cfg)}}
    if result["method"] == "preserved_manual":
        result["coverage_status"] = "manual_review_preserved"
    else:
        path = selection_path(cfg, record)
        selection = json.loads(path.read_text()) if path.exists() else None
        result["coverage_status"] = "selection_missing"
        if selection is not None:
            if selection.get("status") == "deferred":
                if (selection.get("schema") != SCHEMA or selection.get("purpose") != "ocular_component_selection_only"
                        or selection.get("binding", {}).get("source_filename") != record.source_filename):
                    raise RuntimeError("Deferred selection belongs to another source or schema.")
                result["coverage_status"] = "selection_deferred"
            else:
                cache = load_cached_aid(cfg, record)
                validate_selection(selection, record, cache)
                result["selection_sha256"] = digest(selection)
                result["component_binding"] = component_binding(record, cache)
                result["components"] = {key: selection[key] for key in ("blink_components", "eye_movement_components")}
                result["coverage_status"] = "selected_components_only"
                result["input_sha256"] = digest(result)
                previous = cached_mask(cfg, record, result["input_sha256"])
                if previous is not None:
                    return raw, previous
                sources = cache["native_sources"]
                for key, label, detector in (("blink", AUTO_LABEL, detect_blinks),
                                             ("eye_movement", EYE_AUTO_LABEL, detect_eye_movements)):
                    components = selection[f"{key}_components"]
                    settings = result["settings"][key]
                    intervals, scores, metadata = detector(sources, components, settings)
                    metadata["unscored_intervals_raw"] = [[int(a + raw.first_samp), int(b + raw.first_samp)]
                        for a, b in finite_stretches(~np.isfinite(scores))]
                    metadata["provisional_intervals_raw"] = [[int(a + raw.first_samp), int(b + raw.first_samp)] for a, b in intervals]
                    identity = {"model": cache["metadata"]["model_sha256"], "components": components, "settings": settings}
                    if key == "eye_movement":
                        identity["method"] = EYE_METHOD
                    metadata.update(fingerprint=digest(identity), model_sha256=cache["metadata"]["model_sha256"],
                                    selected_components=components, removed_spans_raw=[])
                    result["detectors"][key] = metadata
                    for a, b in intervals:
                        raw.annotations.append((a + raw.first_samp) / raw.info["sfreq"], (b - a) / raw.info["sfreq"], label)
    result["annotations"] = mask_rows(raw)
    result["payload_sha256"] = digest(result)
    return raw, result


def save_mask(cfg, record, payload):
    """Write a content-addressed mask once, detecting corrupt existing sidecars."""
    path = triage_path(cfg, triage_root(cfg) / "proposals" / cfg["layout"]["subject_bids_id"] / record.source_stem / f'{payload["payload_sha256"]}.json')
    if path.exists():
        if json.loads(path.read_text()) != payload:
            raise RuntimeError(f"Existing immutable triage mask differs: {path}")
    else:
        write_json(path, payload, raw_root=cfg["resolved_paths"]["data_root"], overwrite=False)
    return path


def load_review_mask(cfg, record, path):
    """Verify and reconstruct an exact provisional mask for targeted manual review.

    No detector reruns or candidate reseeding occurs. A later manual SAVE always
    wins: refuse to load a provisional set over any existing manual work.
    """
    if has_manual_work(cfg, record):
        raise RuntimeError("Manual work already exists for this source; reopen that review instead of triage proposals.")
    payload = json.loads(Path(path).read_text())
    expected = payload.get("payload_sha256")
    if digest({k: v for k, v in payload.items() if k != "payload_sha256"}) != expected:
        raise RuntimeError("Triage proposal content hash mismatch.")
    raw, _, binding = load_context(cfg, record)
    if (payload.get("schema") != SCHEMA or payload.get("purpose") != PURPOSE
            or payload.get("input_binding") != binding or payload.get("method") != "provisional_automatic"):
        raise RuntimeError("Triage source/muscle inputs changed; regenerate the affected retention report.")
    settings = {"blink": blink_settings(cfg), "eye_movement": eye_movement_settings(cfg)}
    if payload.get("settings") != settings:
        raise RuntimeError("Triage detector settings changed; regenerate the affected report.")
    selected = selection_path(cfg, record)
    if payload.get("selection_sha256") is not None:
        selection = json.loads(selected.read_text())
        cache = load_cached_aid(cfg, record)
        validate_selection(selection, record, cache)
        if digest(selection) != payload["selection_sha256"]:
            raise RuntimeError("Component selections changed after triage; regenerate the affected report.")
    elif selected.exists() and json.loads(selected.read_text()).get("status") == "selected":
        raise RuntimeError("New component choices exist; regenerate the affected retention report.")
    return payload


def seed_review_mask(raw, payload):
    """Replace BAD context only in a disposable browser with the verified exact mask."""
    raw.annotations.delete([i for i, a in enumerate(raw.annotations)
                            if a["description"].lower().startswith("bad") and a["description"].lower() != "bad_filter_edge"])
    for a in payload["annotations"]:
        raw.annotations.append(a["onset"], a["duration"], a["description"])


# -----------------------------------------------------------------------------
# Metadata-only retention accounting; all overlapping causes remain inspectable
# -----------------------------------------------------------------------------

def count_trials(rows, raw, payload, epoch):
    """Evaluate eligibility, complete native windows, BAD overlaps and scoring gaps.

    Match MNE's strict overlap against the half-open sample slice, not a
    duration-percentage approximation. Even zero-duration BAD annotations inside
    the slice have MNE's rejection semantics. All BAD causes are retained, rather
    than only MNE's first matching label. No data samples are loaded.
    """
    sfreq, first, end = float(raw.info["sfreq"]), int(raw.first_samp), int(raw.first_samp + raw.n_times)
    output = []
    for original in rows:
        if (float(original["sfreq_hz"]) != sfreq or int(original["recording_first_sample_raw"]) != first
                or int(original["recording_end_sample_raw_exclusive"]) != end):
            raise RuntimeError("Accepted trial coordinates differ from the exact SSS recording.")
        row = dict(original)
        eligible = analysis_eligible(row)
        row.update(analysis_eligible=eligible, annotation_survives=False, incomplete_window=False,
                   bad_categories="", unscored_groups="", mask_method=payload["method"],
                   ocular_coverage_status=payload["coverage_status"])
        if eligible:
            onset = int(row["epoch_event_sample_raw"])
            left = onset + int(round(float(epoch["tmin_sec"]) * sfreq))
            stop = onset + int(round(float(epoch["tmax_sec"]) * sfreq)) + 1
            row.update(epoch_start_raw=left, epoch_stop_raw_exclusive=stop,
                       incomplete_window=left < first or stop > end)
            labels = sorted({a["description"] for a in payload["annotations"]
                if a["onset"] - raw.first_time < (stop - first) / sfreq
                and a["onset"] - raw.first_time + a["duration"] > (left - first) / sfreq})
            row["bad_categories"] = ";".join(labels)
            row["annotation_survives"] = not row["incomplete_window"] and not labels
            gaps = []
            if payload["method"] != "preserved_manual":
                for key in ("blink", "eye_movement"):
                    detector = payload["detectors"].get(key)
                    if not detector or any(a < stop and b > left for a, b in detector["unscored_intervals_raw"]):
                        gaps.append(key)
            row["unscored_groups"] = ";".join(gaps)
        row["triage_outcome"] = ("excluded_by_accepted_trial_policy" if not eligible else
            "incomplete_window" if row["incomplete_window"] else
            "overlaps_BAD" if row["bad_categories"] else "annotation_survivor")
        output.append(row)
    return output


def summarize(rows, cfg, sources):
    """Return condition/run/identity cells, including every expected zero-count cell."""
    conditions, identities, per_run = [], [], []
    low, high = cfg["subject_profile"]["identity_range"]
    for name, code in cfg["epoching"]["event_id_map"].items():
        subset = [r for r in rows if int(r["condition_code"]) == int(code)]
        eligible = sum(r["analysis_eligible"] for r in subset)
        survivors = sum(r["annotation_survives"] for r in subset)
        conditions.append({"condition": name, "condition_code": code, "eligible": eligible,
            "annotation_survivors": survivors, "retention_percent": 100 * survivors / eligible if eligible else None,
            "incomplete_windows": sum(r["analysis_eligible"] and r["incomplete_window"] for r in subset),
            "windows_overlapping_BAD": sum(r["analysis_eligible"] and bool(r["bad_categories"]) for r in subset),
            "survivors_with_unscored_ocular": sum(r["annotation_survives"] and bool(r["unscored_groups"]) for r in subset)})
        for identity in range(int(low), int(high) + 1):
            cell = [r for r in subset if r.get("trial_sequence_identity") == identity]
            identities.append({"condition": name, "identity": identity,
                "eligible": sum(r["analysis_eligible"] for r in cell),
                "annotation_survivors": sum(r["annotation_survives"] for r in cell)})
        for source in sources:
            cell = [r for r in subset if r["source_fif"] == source]
            cell_eligible = sum(r["analysis_eligible"] for r in cell)
            cell_survivors = sum(r["annotation_survives"] for r in cell)
            per_run.append({"source_fif": source, "condition": name,
                "eligible": cell_eligible, "annotation_survivors": cell_survivors,
                "retention_percent": 100 * cell_survivors / cell_eligible if cell_eligible else None})
    return conditions, identities, per_run


def csv_content(rows):
    """Export heterogeneous ledgers deterministically, retaining accepted metadata."""
    buffer = io.StringIO()
    fields = sorted({key for row in rows for key in row})
    if fields:
        writer = csv.DictWriter(buffer, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return buffer.getvalue()


def common_arguments(parser):
    """Share discovery flags between the bulk GUI and noninteractive report tool."""
    parser.add_argument("--config", default="meg/params/pipeline_default.yaml")
    parser.add_argument("--site-config")
    parser.add_argument("--subject-config-dir", default="meg/params")
    parser.add_argument("--subjects", nargs="+", type=int, choices=range(4, 40), default=list(range(4, 40)))


def calculate_participant(task):
    """Validate/count one participant; only its source-specific masks are written."""
    args, subject, expected_root = task
    all_trials, condition_rows, identity_rows, run_rows = [], [], [], []
    subject_id = f"sub-{subject:02}"
    print(f"{subject_id}: validating accepted inputs and calculating annotation retention", flush=True)
    subject_errors, trials, masks = [], [], []
    try:
        cfg = subject_config(args, subject)
        if not cfg["epoching"]["reject_by_annotation"]:
            raise RuntimeError("Annotation triage requires the accepted reject_by_annotation=true policy.")
        if triage_root(cfg) != expected_root:
            raise RuntimeError("All selected subjects must use the same derivatives root.")
        require_review(cfg, "03_A_MaxFilt", stage_output_dirs(cfg)["reports"])
        review = (cfg["epoching"].get("timing_review_path") or Path(cfg["resolved_paths"]["derivatives_root"]) /
                  "visual_timing/s3-v002/reviews" / f"{subject_id}_review.json")
        prepared, gate = load_reviewed_timing(cfg, load_prepared_decisions(cfg), review)
        records = discover_raw_run_records(cfg)
        if {r.source_filename for r in records} != {r["source_fif"] for r in prepared}:
            raise RuntimeError("Configured recordings and accepted ledger source coverage differ.")
        for record in records:
            try:
                cfg = subject_config(args, subject)
                raw, mask = provisional_mask(cfg, record)
                selected = [row for row in prepared if row["source_fif"] == record.source_filename]
                counted = count_trials(selected, raw, mask, cfg["epoching"])
                # A SAVE during detection must not publish stale automatic work.
                latest = subject_config(args, subject)
                _, _, latest_binding = load_context(latest, record)
                if latest_binding != mask["input_binding"]:
                    raise RuntimeError("Manual/source inputs changed during counting; rerun this subject.")
                path = save_mask(cfg, record, mask)
                for row in counted:
                    row["subject_bids_id"] = subject_id
                    row["triage_proposals"] = str(path)
                trials.extend(counted)
                masks.append({"source_fif": record.source_filename, "path": str(path),
                              "method": mask["method"], "coverage_status": mask["coverage_status"]})
                print(f"{subject_id} {record.source_filename}: mask and counts ready", flush=True)
            except (OSError, ValueError, RuntimeError, KeyError) as error:
                subject_errors.append({"subject_bids_id": subject_id, "source_fif": record.source_filename, "error": str(error)})
                print(f"{subject_id} {record.source_filename}: INPUT ERROR {error}", flush=True)
        conditions, identities, runs = summarize(trials, cfg, [r.source_filename for r in records])
        for table, destination in ((conditions, condition_rows), (identities, identity_rows), (runs, run_rows)):
            destination.extend({"subject_bids_id": subject_id,
                "input_status": "incomplete_inputs" if subject_errors else "counted", **row} for row in table)
        summary = {"subject_bids_id": subject_id, "status": "incomplete_inputs" if subject_errors else "counted",
            "minimum_condition_survivors": min((r["annotation_survivors"] for r in conditions), default=0),
            "total_annotation_survivors": sum(r["annotation_survives"] for r in trials),
            "missing_identity_cells": sum(r["annotation_survivors"] == 0 for r in identities),
            "runs_missing_selections": sum(r["coverage_status"] in ("selection_missing", "selection_deferred") for r in masks),
            "survivors_with_unscored_ocular": sum(r["annotation_survives"] and bool(r["unscored_groups"]) for r in trials),
            "masks": masks, "timing_gate": gate}
        all_trials.extend(trials)
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        subject_errors.append({"subject_bids_id": subject_id, "source_fif": "", "error": str(error)})
        summary = {"subject_bids_id": subject_id, "status": "incomplete_inputs"}
    print(f"{subject_id}: finished {len(masks)} run(s), {len(subject_errors)} input error(s)", flush=True)
    return summary, all_trials, condition_rows, identity_rows, run_rows, subject_errors


def run_report(args):
    """Build a new report snapshot; errors remain visible and are never ranked as zeros."""
    first_cfg = subject_config(args, args.subjects[0])
    label = args.label or "retention-" + utc_now().replace(":", "").replace("-", "") + "-" + uuid4().hex[:6]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label) or label in (".", ".."):
        raise ValueError("--label must be a simple folder name.")
    output = triage_path(first_cfg, triage_root(first_cfg) / "reports" / label)
    if output.exists():
        raise FileExistsError(f"Report snapshots are immutable; choose a new --label: {output}")
    raw_root = first_cfg["resolved_paths"]["data_root"]
    subjects, all_trials, condition_rows, identity_rows, run_rows, errors = [], [], [], [], [], []
    selected_subjects = sorted(set(args.subjects))
    workers = participant_workers(getattr(args, "jobs", None), len(selected_subjects))
    print(f"Counting {len(selected_subjects)} participants with {workers} participant worker(s)", flush=True)
    tasks = [(args, subject, triage_root(first_cfg)) for subject in selected_subjects]
    results = map_participants(calculate_participant, tasks, workers)
    # Shared tables are published once, after all independent participants finish.
    for summary, trials, conditions, identities, runs, subject_errors in results:
        subjects.append(summary)
        all_trials.extend(trials)
        condition_rows.extend(conditions)
        identity_rows.extend(identities)
        run_rows.extend(runs)
        errors.extend(subject_errors)
    ranked = sorted([r for r in subjects if r["status"] == "counted"],
                    key=lambda r: (r["minimum_condition_survivors"], -r["missing_identity_cells"], r["subject_bids_id"]))
    for rank, row in enumerate(ranked, 1):
        row["review_priority"] = rank
    payload = {"schema": SCHEMA, "purpose": PURPOSE, "generated_utc": utc_now(),
        "execution": {"participant_workers": workers},
        "status": "incomplete_inputs" if errors else "completed", "subjects": subjects,
        "errors": errors, "interpretation": "Annotation-only conditional counts; no amplitude rejection or claim of clean data."}
    write_json(output / "report.json", payload, raw_root=raw_root, overwrite=False)
    reasons = Counter()
    for row in all_trials:
        if row["analysis_eligible"]:
            for reason in filter(None, row["bad_categories"].split(";")):
                reasons[row["subject_bids_id"], row["source_fif"], row["condition_code"], reason] += 1
            if row["incomplete_window"]:
                reasons[row["subject_bids_id"], row["source_fif"], row["condition_code"], "incomplete_window"] += 1
    reason_rows = [{"subject_bids_id": subject, "source_fif": source, "condition_code": condition,
        "reason": reason, "affected_eligible_windows": count}
        for (subject, source, condition, reason), count in sorted(reasons.items())]
    for name, rows in (("trial_retention", all_trials), ("condition_retention", condition_rows),
                       ("identity_retention", identity_rows), ("run_condition_retention", run_rows),
                       ("overlapping_rejection_reasons", reason_rows), ("errors", errors)):
        write_text(output / f"{name}.csv", csv_content(rows), raw_root=raw_root, overwrite=False)
    export = [{k: v for k, v in row.items() if k not in ("masks", "timing_gate")} for row in ranked + [r for r in subjects if r["status"] != "counted"]]
    write_text(output / "review_priority.csv", csv_content(export), raw_root=raw_root, overwrite=False)
    lines = ["# Preliminary annotation retention", "", payload["interpretation"],
        "", "Rank by lowest surviving condition count; missing identity cells break ties.",
        "Missing selections and unscored periods indicate incomplete ocular screening.",
        "Reason counts overlap and must not be added to estimate unique rejected trials.",
        "Participants with input errors are unranked; their partial tables are diagnostic only.", "",
        "| Priority | Participant | Weakest condition | Missing identity cells | Missing selections | Survivors with unscored eyes |",
        "| --- | --- | ---: | ---: | ---: | ---: |"]
    for row in ranked:
        lines.append(f'| {row["review_priority"]} | {row["subject_bids_id"]} | {row["minimum_condition_survivors"]} | {row["missing_identity_cells"]} | {row["runs_missing_selections"]} | {row["survivors_with_unscored_ocular"]} |')
    lines.extend(["", "Each trial CSV row links its exact triage_proposals JSON. Pass that path to",
                  "manual_eye_annotations.py --triage-proposals PATH for targeted review.",
                  "Existing manual work must instead be reopened normally with --overwrite.", ""])
    for error in errors:
        lines.append(f'- INPUT ERROR {error["subject_bids_id"]} {error["source_fif"]}: {error["error"]}')
    write_text(output / "report.md", "\n".join(lines), raw_root=raw_root, overwrite=False)
    make_plot(output, condition_rows, ranked)
    print(f"Saved preliminary report: {output / 'report.md'}", flush=True)
    return 1 if errors else 0


def make_plot(output, conditions, ranked):
    """Render a scientific static count heatmap with participant/condition labels."""
    if not ranked:
        return
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.ticker import MaxNLocator
    names = sorted({r["condition"] for r in conditions})
    values = np.array([[next(r["annotation_survivors"] for r in conditions
        if r["subject_bids_id"] == subject["subject_bids_id"] and r["condition"] == name)
        for name in names] for subject in ranked])
    fig = Figure(figsize=(8, max(3, .3 * len(ranked) + 2)))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    artist = ax.imshow(values, aspect="auto", cmap="viridis", vmin=0, vmax=max(1, int(values.max())))
    ax.set(xticks=range(len(names)), xticklabels=names, yticks=range(len(ranked)),
           yticklabels=[r["subject_bids_id"] for r in ranked], title="Preliminary annotation survivors (review priority order)")
    for i, row in enumerate(values):
        for j, value in enumerate(row):
            ax.text(j, i, str(value), ha="center", va="center",
                    color="white" if float(artist.norm(value)) < .55 else "black", fontsize=8)
    colorbar = fig.colorbar(artist, ax=ax, label="Annotation survivors")
    colorbar.locator = MaxNLocator(nbins=5, integer=True)
    colorbar.update_ticks()
    fig.tight_layout()
    fig.savefig(output / "condition_retention.png", dpi=160)


def main(argv=None):
    """Generate one new snapshot; repeat selected subjects after targeted manual SAVE."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common_arguments(parser)
    add_parallel_argument(parser)
    parser.add_argument("--label", help="New immutable report-folder name; default is a unique timestamp.")
    return run_report(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
