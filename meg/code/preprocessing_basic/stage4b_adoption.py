#!/usr/bin/env python3
"""Adopt current annotation masks under Damiano's explicit 2026-11-07 policy.

Usage from the DAD root (cluster: add --site-config meg/params/pipeline_cluster.yaml)
--------------------------------------------------------------------------------
    python meg/code/preprocessing_basic/stage4b_adoption.py plan \
      --snapshot retention-20261107-v1 --label current-masks-20261109-v1 --jobs 8
    python meg/code/preprocessing_basic/stage4b_adoption.py apply \
      --label current-masks-20261109-v1 --subject 04
    python meg/code/preprocessing_basic/stage4b_adoption.py reconcile \
      --label current-masks-20261109-v1 --accept-policy --jobs 8

Inputs: exact immutable retention masks, accepted 4A, Stage 2/3 accounting and
timing/application evidence, current manual receipts/files, and validated cached
ocular models/selections. Plan validates the ENTIRE working cohort before any
production write. Later manual saves supersede the snapshot during planning;
after planning, changed inputs fail closed and require a new plan. No detector,
EDF/SQUID preparation or source/model validation bypass is introduced.

Data flow: parallel read-only validation/counts -> immutable portable plan ->
separate participant adoption processes -> joined independent reconciliation ->
optional policy-specific 4B processing acceptance -> one cohort publication.
Apply writes standard Stage 4B FIF/CSV/report/manifest outputs, preserving old
outputs by rename, all manual files/receipts and all non-BAD source context.
Source BAD categories remain, except an explicitly reviewed muscle replacement.
Reviewed-empty muscle and ocular edits/deletions retain their saved provenance.
Automatic adoption is NEVER a manual SAVE or a claim of manually reviewed data.

All numerical threads are bounded by the caller; cohort workers use one thread
each. --jobs exposes concurrency; cluster array tasks isolate participant writes.
Fit/model decisions, application QC, cleaned-output acceptance and epoching remain
later explicit steps. P15 data/provenance are retained but never adopted here.
"""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import fcntl
import importlib.util
import json
from pathlib import Path
import re
import shutil
import tempfile

import mne
import numpy as np

import annotation_triage as triage
from blink_annotation import blink_settings, eye_movement_settings, cache_paths, load_cached_aid
from continuous_cleaning import (annotation_rows, bind_report, check_preservation,
    coordinates, digest, file_hash, file_stamp, final_filter_configuration,
    final_filter_summary, require_review, reviewed_muscle_override, stage_configuration)
from manual_eye_annotations import load_manual, manual_path
from participant_parallel import add_parallel_argument, map_participants, participant_workers
from preprocessing_stage_utils import (derivative_path_for_run, discover_raw_run_records,
    guarded_path, stage_output_dirs, update_source_manifest_stage, update_tuning_file,
    utc_now, write_json, write_text)
from render_stage03_mne_reports import render_stage_report
from trial_preparation import load_prepared_decisions
from visual_timing_utils import load_reviewed_timing

POLICY = "current_annotations_defer_detailed_review_2026-11-07"
STAGE = "03_B_ArtefactAnnotation"
COHORT = tuple(i for i in range(4, 40) if i != 15)
TRIAL_FIELDS = ("trial_key", "source_fif", "analysis_eligible", "annotation_survives",
    "incomplete_window", "bad_categories", "unscored_groups")


# -----------------------------------------------------------------------------
# Portable immutable identities and current scientific prerequisites
# -----------------------------------------------------------------------------

def root(cfg):
    """Resolve the dedicated adoption namespace outside source data."""
    return guarded_path(Path(cfg["resolved_paths"]["derivatives_root"]) / "stage4b_adoption",
                        cfg["resolved_paths"]["data_root"])


def relative(cfg, path):
    """Store shared resources relative to derivatives, independently of mounts."""
    return str(Path(path).resolve().relative_to(Path(cfg["resolved_paths"]["derivatives_root"]).resolve()))


def resource(cfg, path):
    """Refuse escaping/symlinked resource paths in an adoption plan."""
    base = Path(cfg["resolved_paths"]["derivatives_root"]).resolve()
    result = (base / path).resolve()
    if not result.is_relative_to(base):
        raise ValueError("Adoption resources must remain inside derivatives.")
    return result


def portable(value):
    """Remove display-only absolute paths from already hash-bound timing records."""
    if isinstance(value, dict):
        return {k: portable(v) for k, v in value.items() if not k.endswith("_path")}
    if isinstance(value, list):
        return [portable(v) for v in value]
    return value


def base_configuration(cfg):
    """Bind all annotation/final settings without the newly added adoption pointer."""
    plain = deepcopy(cfg)
    plain["artifact_annotation"].pop("mask_adoption", None)
    return stage_configuration(plain, STAGE, include_edge=False)


def timing(cfg):
    """Revalidate unchanged accepted trial eligibility and Stage 3 application."""
    review = (cfg["epoching"].get("timing_review_path") or
        Path(cfg["resolved_paths"]["derivatives_root"]) / "visual_timing/s3-v002/reviews" /
        f'{cfg["layout"]["subject_bids_id"]}_review.json')
    rows, gate = load_reviewed_timing(cfg, load_prepared_decisions(cfg), review)
    return rows, {"gate": portable(gate), "prepared_rows_sha256": digest(portable(rows))}


def protected(cfg, path, *, small=True, basis=False):
    """Freeze small content hashes and large FIF stat identities explicitly."""
    return {"relative_path": relative(cfg, path), "stamp": file_stamp(path),
            "sha256": file_hash(path) if small else None, "basis": basis}


def check_payload(mask):
    """Verify a content-addressed proposal before interpreting its scientific rows."""
    if (mask.get("schema") != triage.SCHEMA or mask.get("purpose") != triage.PURPOSE or
        digest({k: v for k, v in mask.items() if k != "payload_sha256"}) != mask.get("payload_sha256")):
        raise RuntimeError("Invalid or changed immutable triage mask.")
    for row in mask["annotations"]:
        if (not row["description"].lower().startswith("bad") or
            row["description"].lower() == "bad_filter_edge" or
            not np.isfinite(row["onset"]) or not np.isfinite(row["duration"]) or row["duration"] < 0):
            raise RuntimeError("Invalid BAD interval in the adoption plan.")


def trial_footprint(rows):
    """Compare every trial outcome/cause, rather than only equal cohort totals."""
    return [{k: row[k] for k in TRIAL_FIELDS} for row in rows]


def validate_participant(task):
    """Read/validate one whole participant; no production or manual files are written."""
    args, subject, snapshot_row = task
    cfg = triage.subject_config(args, subject)
    dirs = stage_output_dirs(cfg)
    require_review(cfg, "03_A_MaxFilt", dirs["reports"])
    prepared, gate = timing(cfg)
    if portable(snapshot_row["timing_gate"]) != gate["gate"]:
        raise RuntimeError("Accepted Stage 3 timing/application differs from the retention snapshot.")
    records = discover_raw_run_records(cfg)
    if {r.source_filename for r in records} != {r["source_fif"] for r in prepared}:
        raise RuntimeError("Accepted ledger and acquisition coverage differ.")
    if not cfg["epoching"]["reject_by_annotation"]:
        raise RuntimeError("Adoption requires reject_by_annotation=true.")
    snapshots = {row["source_fif"]: row for row in snapshot_row["masks"]}
    if set(snapshots) != {r.source_filename for r in records}:
        raise RuntimeError("Snapshot and current acquisition coverage differ.")
    result = {"subject": subject, "base_configuration_sha256": base_configuration(cfg),
        "timing": gate, "runs": [], "protected_files": [], "expected_trials": []}
    for record in records:
        print(f'{cfg["layout"]["subject_bids_id"]} {record.source_filename}: validating mask/manual/model inputs', flush=True)
        old = snapshots[record.source_filename]
        suffix = old["path"].split("/annotation_triage/", 1)[1]
        proposal = resource(cfg, "annotation_triage/" + suffix)
        mask = json.loads(proposal.read_text())
        check_payload(mask)
        raw, source, binding = triage.load_context(cfg, record)
        manual_receipt = None
        changed_manual = False
        if triage.has_manual_work(cfg, record):
            _, manual_receipt = load_manual(cfg, record, raw, source)
            # This branch runs no detector and discovers reviews saved since the
            # retention snapshot. Full saved replacement/deletion semantics win.
            raw, current = triage.provisional_mask(cfg, record)
            changed_manual = current["payload_sha256"] != mask["payload_sha256"]
            mask = current
            result["protected_files"].append(protected(cfg, manual_path(cfg, record)))
        else:
            if (mask["method"] != "provisional_automatic" or mask["input_binding"] != binding or
                mask["coverage_status"] != "selected_components_only" or
                mask["settings"] != {"blink": blink_settings(cfg), "eye_movement": eye_movement_settings(cfg)}):
                raise RuntimeError(f"Stale/incomplete mask for {record.source_filename}; create a new retention snapshot.")
            cache = load_cached_aid(cfg, record)
            selection = json.loads(triage.selection_path(cfg, record).read_text())
            triage.validate_selection(selection, record, cache)
            if (mask["selection_sha256"] != digest(selection) or
                mask["component_binding"] != triage.component_binding(record, cache)):
                raise RuntimeError("Mask component/source/model identity changed.")
            for key, path in cache_paths(cfg, record).items():
                result["protected_files"].append(protected(cfg, path, small=key != "sources"))
            result["protected_files"].append(protected(cfg, triage.selection_path(cfg, record)))
        result["protected_files"].extend([protected(cfg, proposal), protected(cfg, source, small=False)])
        basis = derivative_path_for_run(record, dirs["sss"], "_ann.fif")
        result["protected_files"].append(protected(cfg, basis, small=False, basis=True))
        selected = [row for row in prepared if row["source_fif"] == record.source_filename]
        counted = triage.count_trials(selected, raw, mask, cfg["epoching"])
        result["expected_trials"].extend(trial_footprint(counted))
        result["runs"].append({"source_filename": record.source_filename,
            "mask": mask, "manual_eye_review": manual_receipt or {"status": "unreviewed"},
            "manual_superseded_snapshot": changed_manual, "snapshot_mask": relative(cfg, proposal)})
    # A SAVE during validation must not freeze stale automatic or manual inputs.
    latest = triage.subject_config(args, subject)
    if base_configuration(latest) != result["base_configuration_sha256"]:
        raise RuntimeError("Annotation configuration/manual receipts changed while planning.")
    return result


def read_plan(cfg, label):
    """Require a complete immutable plan, for precisely the authorized cohort."""
    path = root(cfg) / label / "plan.json"
    plan = json.loads(path.read_text())
    if (plan.get("policy") != POLICY or plan.get("cohort") != list(COHORT) or
        digest({k: v for k, v in plan.items() if k != "plan_sha256"}) != plan.get("plan_sha256")):
        raise RuntimeError("Invalid/changed adoption plan or working cohort.")
    return path, plan


def history_path(cfg, plan, path):
    """Preserve earlier working outputs with native stat identities, without copying GBs."""
    return (Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "history/mask_adoption" /
            plan["plan_sha256"] / Path(path).name)


def validate_frozen(cfg, plan, participant):
    """Check current inputs even after original annotated inputs have been archived."""
    if base_configuration(cfg) != participant["base_configuration_sha256"]:
        raise RuntimeError("Current settings/manual receipts differ from the frozen plan; replan.")
    require_review(cfg, "03_A_MaxFilt", stage_output_dirs(cfg)["reports"])
    prepared, gate = timing(cfg)
    if gate != participant["timing"]:
        raise RuntimeError("Accepted Stage 2/3 evidence changed after planning.")
    if file_hash(resource(cfg, plan["snapshot"]["relative_path"])) != plan["snapshot"]["sha256"]:
        raise RuntimeError("Retention snapshot changed after planning.")
    for item in participant["protected_files"]:
        path = resource(cfg, item["relative_path"])
        archived = history_path(cfg, plan, path)
        if item["basis"] and archived.exists():
            path = archived
        if file_stamp(path) != item["stamp"] or item["sha256"] and file_hash(path) != item["sha256"]:
            raise RuntimeError(f"Adoption input changed: {path}")
    for row in participant["runs"]:
        check_payload(row["mask"])
    return prepared


def adoption_pointer(cfg, path):
    """Keep policy provenance separate from every manual-review receipt."""
    return {"policy": POLICY, "plan_relative_path": relative(cfg, path), "plan_file_sha256": file_hash(path)}


def require_operative_adoption(cfg, payload, *, require_publication=False):
    """Downstream guard: acceptance remains bound to current source/manual/model evidence."""
    pointer = cfg["artifact_annotation"].get("mask_adoption")
    if not pointer:
        if payload.get("annotation_review_policy") == POLICY:
            raise RuntimeError("Adopted Stage 4B report is missing its policy binding.")
        return
    path = resource(cfg, pointer["plan_relative_path"])
    if pointer["policy"] != POLICY or file_hash(path) != pointer["plan_file_sha256"]:
        raise RuntimeError("Stage 4B adoption plan changed.")
    _, plan = read_plan(cfg, path.parent.name)
    participant = next(row for row in plan["participants"] if row["subject"] == int(cfg["subject"]["id"]))
    validate_frozen(cfg, plan, participant)
    if payload.get("adoption_plan_sha256") != plan["plan_sha256"] or payload.get("annotation_review_policy") != POLICY:
        raise RuntimeError("Stage 4B report belongs to a different annotation policy/plan.")
    if require_publication:
        publication = path.parent / "policy_acceptance.json"
        if not publication.is_file():
            raise RuntimeError("Whole-cohort policy acceptance publication is missing; final ICA must wait.")
        joined = json.loads(publication.read_text())
        rows = joined.get("participants", [])
        if (joined.get("status") != "accepted_for_final_ica_fit" or joined.get("policy") != POLICY or
            joined.get("adoption_plan_sha256") != plan["plan_sha256"] or
            sorted(row["subject"] for row in rows) != list(COHORT) or
            not any(row["subject"] == participant["subject"] and row["evidence_sha256"] == payload["evidence_sha256"] for row in rows)):
            raise RuntimeError("Whole-cohort policy publication differs from the current Stage 4B evidence.")


# -----------------------------------------------------------------------------
# Standard Stage 4B outputs: no detector rerun and no fabricated manual SAVE
# -----------------------------------------------------------------------------

def compose_annotations(raw, mask):
    """Replace only BAD rows with the validated complete operative mask.

    MNE attached onsets already include first_time. Appending after deletion,
    rather than set_annotations(), avoids applying that offset twice when the
    recording has no measurement date. Non-BAD boundaries/context retain extras.
    """
    remaining = list(mask["annotations"])
    removed = []
    for index, annotation in enumerate(raw.annotations):
        if not annotation["description"].lower().startswith("bad") or annotation["description"].lower() == "bad_filter_edge":
            continue
        row = {k: annotation[k] for k in ("onset", "duration", "description")}
        if row in remaining:
            remaining.remove(row)  # Keep source channel/extras metadata in place.
        else:
            removed.append(index)
    raw.annotations.delete(removed)
    for row in remaining:
        raw.annotations.append(row["onset"], row["duration"], row["description"])
    return raw


def stage4b_module():
    """Reuse existing artifact signatures/lineage interfaces without running main."""
    spec = importlib.util.spec_from_file_location("dad_stage4b", Path(__file__).with_name("03_B_ArtefactAnnotation.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def serialized_mask(raw, mask):
    """Predict the public MNE FIF annotation representation without signal IO.

    FIF stores endpoints as float32, and Raw loading crops/rebases annotations.
    A native-grid endpoint exactly on an epoch boundary can consequently change
    strict overlap after saving. Round-trip the small annotation object through
    MNE itself, then reproduce Raw loading with the public setter. This predicts
    the operative representation independently of the new production FIF; it
    never rounds manual decisions to a new grid or introduces an artifact margin.
    """
    rows = mask["annotations"]
    annotations = mne.Annotations([row["onset"] for row in rows],
        [row["duration"] for row in rows], [row["description"] for row in rows],
        orig_time=raw.annotations.orig_time)
    with tempfile.TemporaryDirectory(prefix="dad-mask-serialization-") as folder:
        path = Path(folder) / "mask-annot.fif"
        annotations.save(path, verbose="ERROR")
        stored = mne.read_annotations(path)
    if stored.orig_time is None:
        stored.onset -= raw.first_time
    preview = raw.copy().set_annotations(stored, emit_warning=False)
    return {**mask, "annotations": triage.mask_rows(preview)}


def verify_saved(cfg, record, entry, path, prepared):
    """Reconcile saved masks against an independently predicted FIF round-trip."""
    saved = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
    if coordinates(saved) != entry["mask"]["input_binding"]["coordinates"]:
        raise RuntimeError("Adoption changed native recording coordinates.")
    actual = {**entry["mask"], "annotations": triage.mask_rows(saved)}
    predicted = serialized_mask(saved, entry["mask"])
    canonical = lambda rows: sorted(rows, key=lambda a: (a["onset"], a["duration"], a["description"]))
    # Equal-time labels may retain a different stable order when original
    # source context is preserved. Compare the complete interval multiset;
    # the actual saved ordering remains bound separately by annotation_sha256.
    if canonical(actual["annotations"]) != canonical(predicted["annotations"]):
        raise RuntimeError("Saved operative BAD rows differ from the exact planned FIF serialization.")
    selected = [row for row in prepared if row["source_fif"] == record.source_filename]
    counted = triage.count_trials(selected, saved, actual, cfg["epoching"])
    expected = triage.count_trials(selected, saved, predicted, cfg["epoching"])
    original = trial_footprint(triage.count_trials(selected, saved, entry["mask"], cfg["epoching"]))
    footprint = trial_footprint(counted)
    if footprint != trial_footprint(expected):
        raise RuntimeError("Saved trial outcomes differ from independently predicted operative annotations.")
    changes = [{"trial_key": old["trial_key"], "before_serialization": old, "operative": new}
               for old, new in zip(original, footprint) if old != new]
    return saved, footprint, original, changes


def apply_participant(args):
    """Produce one participant's disjoint Stage 4B derivatives; never accept them."""
    cfg = triage.subject_config(args, args.subject)
    plan_path, plan = read_plan(cfg, args.label)
    participant = next(row for row in plan["participants"] if row["subject"] == args.subject)
    lock_path = root(cfg) / args.label / f"sub-{args.subject:02}.lock"
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        prepared = validate_frozen(cfg, plan, participant)
        dirs = stage_output_dirs(cfg)
        reports = dirs["reports"]
        module = stage4b_module()
        records = discover_raw_run_records(cfg)
        entries = {row["source_filename"]: row for row in participant["runs"]}
        if set(entries) != {r.source_filename for r in records}:
            raise RuntimeError("Acquisition coverage changed since planning.")
        # Archive all former products before any replacement. Rename preserves
        # basis stamps and uses no second GB-scale copy. Retry uses these originals.
        for record in records:
            for suffix in ("_ann.fif", "_ann.csv"):
                path = derivative_path_for_run(record, dirs["sss"], suffix)
                destination = history_path(cfg, plan, path)
                destination.parent.mkdir(parents=True, exist_ok=True)
                if path.exists() and not destination.exists():
                    path.rename(destination)
        for name in (f"{STAGE}_review_report.json", f"{STAGE}_review_report.html", "03_tuning_choices.yaml"):
            path = reports / name
            destination = history_path(cfg, plan, path)
            if path.exists() and not destination.exists():
                shutil.copy2(path, destination)
        run_reports, lineage, outcomes, original_outcomes = [], [], [], []
        for record in records:
            entry = entries[record.source_filename]
            mask = entry["mask"]
            source = derivative_path_for_run(record, dirs["sss"], "_raw_sss.fif")
            output = derivative_path_for_run(record, dirs["sss"], "_ann.fif")
            csv_path = derivative_path_for_run(record, dirs["sss"], "_ann.csv")
            threshold = float(mask["input_binding"]["muscle_basis"]["threshold"])
            receipt = entry["manual_eye_review"]
            reviewed = receipt if receipt.get("status") == "reviewed" else None
            signature = module.artifact_signature(source, threshold, cfg["artifact_annotation"], reviewed)
            # Keep the established muscle signature so ordinary later review can
            # read its context. A second, independent token binds mask adoption.
            token = f'annotation_adoption_sha256={mask["payload_sha256"]}'
            print(f'sub-{args.subject:02} {record.source_filename}: adopting {mask["method"]}', flush=True)
            raw = mne.io.read_raw_fif(source, preload=False, verbose="ERROR")
            comparison = raw.copy()
            if reviewed_muscle_override(reviewed, cfg["artifact_annotation"], threshold):
                comparison.annotations.delete([i for i, a in enumerate(comparison.annotations) if a["description"].lower() == "bad_muscle"])
            if output.exists():
                header = mne.io.read_raw_fif(output, preload=False, verbose="ERROR")
                if token not in str(header.info.get("description")):
                    raise RuntimeError("Existing output is not from this plan; refusing to replace it.")
            else:
                compose_annotations(raw, mask)
                raw.info["description"] = (f"DAD source={record.source_filename}; acquisition_run={record.acquisition_run_label}; "
                    f"stage={STAGE}; muscle_threshold={threshold:g}; stage4_input_sha256={signature}; {token}")
                check_preservation(comparison, raw, annotations=False, samples=False)
                raw.save(output, overwrite=False, verbose="ERROR")
            saved, counted, original, changes = verify_saved(cfg, record, entry, output, prepared)
            preservation = check_preservation(comparison, saved, annotations=False)
            saved.annotations.save(csv_path, overwrite=True)
            outcomes.extend(counted)
            original_outcomes.extend(original)
            report = {"source_filename": record.source_filename, "acquisition_run_label": record.acquisition_run_label,
                "execution_status": "completed", "input_path": str(source), "input_stamp": file_stamp(source),
                "output_fif": str(output), "output_csv": str(csv_path), "muscle_threshold": threshold,
                "manual_eye_review": receipt, "annotation_adoption": {"policy": POLICY,
                    "method": mask["method"], "mask_sha256": mask["payload_sha256"],
                    "coverage_status": mask["coverage_status"], "manual_superseded_snapshot": entry["manual_superseded_snapshot"]},
                "edf_blink_candidates": {"status": "disabled", "review_policy": POLICY},
                "muscle_annotation_count": int(sum(a.lower() == "bad_muscle" for a in saved.annotations.description)),
                "muscle_review_mode": "reviewed_replacement" if reviewed_muscle_override(reviewed, cfg["artifact_annotation"], threshold) else "preserved_current_automatic",
                "eog_status": "legacy_not_in_stage4",
                "misc_channels_preserved": [ch for ch, kind in zip(saved.ch_names, saved.get_channel_types()) if kind == "misc"],
                "annotation_sha256": digest(annotation_rows(saved)), "stage4_input_sha256": signature,
                "serialization_outcome_changes": changes,
                "serialization_contract": "public_mne_annotation_fif_roundtrip_then_raw_setter",
                "preservation": preservation, "total_annotation_count": len(saved.annotations),
                "final_filter_preview": final_filter_summary(saved, cfg["preprocessing"], source_filename=record.source_filename)}
            run_reports.append(report)
            lineage.append(module.build_annotation_lineage_record(record=record, input_path=source,
                           run_report=report, accepted_threshold=threshold))
        if original_outcomes != participant["expected_trials"]:
            raise RuntimeError("Original mask/trial accounting differs from the validated plan.")
        latest = triage.subject_config(args, args.subject)
        validate_frozen(latest, plan, participant)
        pointer = adoption_pointer(cfg, plan_path)
        cfg["artifact_annotation"]["mask_adoption"] = pointer
        sss = require_review(cfg, "03_A_MaxFilt", reports)
        payload = {"stage": "artifact_annotation", "generated_utc": utc_now(), "processing_pass": "review",
            "subject_bids_id": cfg["layout"]["subject_bids_id"], "status": "completed",
            "annotation_review_policy": POLICY, "adoption_plan_sha256": plan["plan_sha256"],
            "mne_version": mne.__version__,
            "scientific_acceptance": "policy_processing_acceptance_pending_cohort_reconciliation",
            "configuration_sha256": stage_configuration(cfg, STAGE, include_edge=False),
            "accepted_sss_lineage": {
                "accepted_configuration_sha256": sss["configuration_sha256"],
                "evidence_sha256": sss["evidence_sha256"],
                "effective_sss_configuration_sha256": stage_configuration(cfg, "03_A_MaxFilt"),
                "configuration_contract": "sss_independent_of_final_temporal_policy_2026-10-28"},
            "runs": run_reports}
        bind_report(payload, [plan_path] + [Path(row[key]) for row in run_reports for key in ("output_fif", "output_csv")],
            {"policy": POLICY, "plan": plan["plan_sha256"], "runs": run_reports}, reports,
            cfg["resolved_paths"]["derivatives_root"])
        write_json(reports / f"{STAGE}_review_report.json", payload, raw_root=cfg["resolved_paths"]["data_root"])
        update_tuning_file(path=reports / "03_tuning_choices.yaml", raw_root=cfg["resolved_paths"]["data_root"],
            subject_id=cfg["subject"]["id"], section="artifact_annotation", values={"mask_adoption": pointer})
        update_source_manifest_stage(manifest_path=reports / "source_manifest.json",
            raw_root=cfg["resolved_paths"]["data_root"], stage_name="artifact_annotation", processing_pass="review",
            records=lineage, status="completed", qc={"annotation_review_policy": POLICY, "plan_sha256": plan["plan_sha256"]})
        render_stage_report(derivatives_root=Path(cfg["resolved_paths"]["derivatives_root"]),
            raw_root=cfg["resolved_paths"]["data_root"], subject_id=cfg["layout"]["subject_bids_id"],
            processing_pass="review", stage_name="artifact_annotation", overwrite=True, primary_report=True)
        print(f'sub-{args.subject:02}: Stage 4B outputs ready; policy acceptance awaits cohort reconciliation', flush=True)


# -----------------------------------------------------------------------------
# Independent output audit, cohort join, then policy acceptance and publication
# -----------------------------------------------------------------------------

def reconcile_participant(task):
    """Reopen saved FIF/CSV and exact input evidence before any cohort acceptance."""
    args, subject = task
    cfg = triage.subject_config(args, subject)
    _, plan = read_plan(cfg, args.label)
    participant = next(row for row in plan["participants"] if row["subject"] == subject)
    prepared = validate_frozen(cfg, plan, participant)
    reports = stage_output_dirs(cfg)["reports"]
    payload = json.loads((reports / f"{STAGE}_review_report.json").read_text())
    require_operative_adoption(cfg, payload)
    # Use the real prerequisite checker with an in-memory candidate receipt;
    # this validates every saved report-bound file without writing acceptance.
    trial = deepcopy(cfg)
    trial.setdefault("stage4_review", {})[STAGE] = {"status": "accepted",
        "configuration_sha256": stage_configuration(cfg, STAGE), "evidence_sha256": payload["evidence_sha256"]}
    require_review(trial, STAGE, reports)
    entries = {row["source_filename"]: row for row in participant["runs"]}
    rows = {row["source_filename"]: row for row in payload["runs"]}
    outcomes, original_outcomes, changes = [], [], []
    module = stage4b_module()
    for record in discover_raw_run_records(cfg):
        print(f"sub-{subject:02} {record.source_filename}: reconciling operative FIF/CSV and trial accounting", flush=True)
        row = rows[record.source_filename]
        output = derivative_path_for_run(record, stage_output_dirs(cfg)["sss"], "_ann.fif")
        saved, counted, original, delta = verify_saved(cfg, record, entries[record.source_filename], output, prepared)
        if digest(annotation_rows(saved)) != row["annotation_sha256"] or not row["preservation"].get("auxiliary_samples_verified"):
            raise RuntimeError("Missing or changed native preservation evidence.")
        source = derivative_path_for_run(record, stage_output_dirs(cfg)["sss"], "_raw_sss.fif")
        module._validate_resumable_annotation(record=record, input_path=source, output_fif=output,
            output_csv=output.with_suffix(".csv"), threshold=row["muscle_threshold"],
            artifact_cfg=cfg["artifact_annotation"], raw_root=cfg["resolved_paths"]["data_root"],
            expected_signature=row["stage4_input_sha256"], reviewed_muscle=row["muscle_review_mode"] == "reviewed_replacement")
        outcomes.extend(counted)
        original_outcomes.extend(original)
        changes.extend(delta)
        if delta != row.get("serialization_outcome_changes", []):
            raise RuntimeError("Recorded FIF serialization changes differ from independent reconciliation.")
    if original_outcomes != participant["expected_trials"]:
        raise RuntimeError("Original annotation/trial accounting differs from the adopted plan.")
    lookup = {row["trial_key"]: row for row in prepared}
    full_rows = [{**lookup[row["trial_key"]], **row} for row in outcomes]
    conditions, identities, run_counts = triage.summarize(full_rows, cfg, list(rows))
    return {"subject": subject, "evidence_sha256": payload["evidence_sha256"],
        "configuration_sha256": stage_configuration(cfg, STAGE),
        "runs": len(rows), "annotation_survivors": sum(r["annotation_survives"] for r in outcomes),
        "serialization_outcome_changes": changes, "condition_retention": conditions,
        "identity_retention": identities, "run_condition_retention": run_counts}


def accept_participant(task):
    """Record explicitly authorized policy processing acceptance, never manual review."""
    args, result = task
    cfg = triage.subject_config(args, result["subject"])
    # Recheck after the cohort validation barrier and immediately before receipt.
    payload = json.loads((stage_output_dirs(cfg)["reports"] / f"{STAGE}_review_report.json").read_text())
    require_operative_adoption(cfg, payload)
    if (payload["evidence_sha256"] != result["evidence_sha256"] or
        stage_configuration(cfg, STAGE) != result["configuration_sha256"]):
        raise RuntimeError("Evidence changed after the cohort audit.")
    trial = deepcopy(cfg)
    trial.setdefault("stage4_review", {})[STAGE] = {"status": "accepted",
        "configuration_sha256": result["configuration_sha256"], "evidence_sha256": result["evidence_sha256"]}
    require_review(trial, STAGE, stage_output_dirs(cfg)["reports"])
    receipt = {"status": "accepted", "acceptance_kind": "explicit_annotation_policy_for_processing",
        "policy": POLICY, "reviewer": "Damiano (explicit policy decision 2026-11-07)",
        "reason": "Adopt current validated masks; defer further detailed interval inspection. Final ICA and cleaned-output QC remain open.",
        "reviewed_utc": utc_now(), "configuration_sha256": result["configuration_sha256"],
        "evidence_sha256": result["evidence_sha256"]}
    update_tuning_file(path=stage_output_dirs(cfg)["reports"] / "03_tuning_choices.yaml",
        raw_root=cfg["resolved_paths"]["data_root"], subject_id=cfg["subject"]["id"],
        section="stage4_review", values={STAGE: receipt})
    return result


def main(argv=None):
    """Keep parallel participant phases and shared immutable publication ordered."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("plan", "apply", "reconcile"))
    parser.add_argument("--config", default="meg/params/pipeline_default.yaml")
    parser.add_argument("--site-config")
    parser.add_argument("--subject-config-dir", default="meg/params")
    parser.add_argument("--snapshot", default="retention-20261107-v1")
    parser.add_argument("--label", required=True)
    parser.add_argument("--subject", type=int, choices=COHORT)
    parser.add_argument("--accept-policy", action="store_true")
    add_parallel_argument(parser)
    args = parser.parse_args(argv)
    for label in (args.label, args.snapshot):
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", label) or label in (".", ".."):
            parser.error("Labels must be simple folder names.")
    cfg = triage.subject_config(args, COHORT[0])
    target = root(cfg) / args.label
    if args.mode == "apply":
        if args.subject is None:
            parser.error("apply requires one --subject; submit the full working-cohort array.")
        apply_participant(args)
        return 0
    workers = participant_workers(args.jobs, len(COHORT))
    print(f'{args.mode}: {len(COHORT)} participants, {workers} worker(s); P15 excluded, data retained', flush=True)
    if args.mode == "plan":
        if target.exists():
            raise FileExistsError("Immutable plan already exists; choose a new label.")
        snapshot_path = resource(cfg, f"annotation_triage/reports/{args.snapshot}/report.json")
        snapshot = json.loads(snapshot_path.read_text())
        if snapshot.get("status") != "completed" or snapshot.get("errors") or len(snapshot["subjects"]) != 36:
            raise RuntimeError("Expected the complete error-free original-cohort snapshot.")
        subjects = {row["subject_bids_id"]: row for row in snapshot["subjects"]}
        results = map_participants(validate_participant,
            [(args, subject, subjects[f"sub-{subject:02}"]) for subject in COHORT], workers)
        # No shared publication (or production write) before the whole cohort joins.
        payload = {"schema": "stage4b-current-mask-adoption-v1", "policy": POLICY,
            "generated_utc": utc_now(), "cohort": list(COHORT), "excluded_participants": {"15": "participant_specific_insufficient_trial_coverage"},
            "snapshot": protected(cfg, snapshot_path), "participants": results,
            "execution": {"participant_workers": workers}}
        payload["plan_sha256"] = digest(payload)
        write_json(target / "plan.json", payload, raw_root=cfg["resolved_paths"]["data_root"], overwrite=False)
        methods = Counter(row["mask"]["method"] for participant in results for row in participant["runs"])
        write_text(target / "README.md", f"# Stage 4B policy adoption plan\n\nPolicy: {POLICY}\n\n"
            f"35 participants; {sum(methods.values())} runs. Methods: {dict(methods)}.\n\n"
            "Plan validation is complete. Production adoption and cohort reconciliation remain pending.\n"
            "Automatic adoption is distinct from manual review. No final ICA decisions or clean-data acceptance.\n",
            raw_root=cfg["resolved_paths"]["data_root"], overwrite=False)
        print(f"Validated immutable adoption plan: {target / 'plan.json'}", flush=True)
    else:
        publication = target / ("policy_acceptance.json" if args.accept_policy else "reconciliation.json")
        if publication.exists():
            raise FileExistsError(f"Immutable reconciliation already exists: {publication}")
        results = map_participants(reconcile_participant, [(args, subject) for subject in COHORT], workers)
        if args.accept_policy:
            results = map_participants(accept_participant, [(args, result) for result in results], workers)
        _, plan = read_plan(cfg, args.label)
        # Publish refreshed operative retention tables only after all participant
        # audits (and requested receipts) join. The cohort receipt is LAST, so
        # interrupted table publication cannot appear ready to downstream fits.
        for field in ("condition_retention", "identity_retention", "run_condition_retention"):
            table = [{"subject_bids_id": f'sub-{result["subject"]:02}', **row}
                     for result in results for row in result[field]]
            path = target / f'{publication.stem}_{field}.csv'
            content = triage.csv_content(table)
            if path.exists():
                if path.read_bytes() != content.encode("utf-8"):
                    raise RuntimeError(f"An immutable operative table differs: {path}")
            else:
                write_text(path, content, raw_root=cfg["resolved_paths"]["data_root"], overwrite=False)
        write_json(publication, {"policy": POLICY, "adoption_plan_sha256": plan["plan_sha256"], "generated_utc": utc_now(),
            "status": "accepted_for_final_ica_fit" if args.accept_policy else "reconciled_not_accepted",
            "participants": results, "final_ica_decisions": "pending", "cleaned_output_acceptance": "pending"},
            raw_root=cfg["resolved_paths"]["data_root"], overwrite=False)
        print(f"Joined cohort reconciliation: {publication}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
