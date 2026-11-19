#!/usr/bin/env python3
"""Prepare the Stage 2 cohort trial accounting without creating epochs.

Purpose and usage
-----------------
Read accepted Stage 1 tables, manifests, deep adjudications, session clues, and
raw FIF *headers*. Retain every recorded condition candidate in a decision
ledger, including excluded and ambiguous candidates. Behavioral-only entries
in a reviewed unavailable tail are a separate coverage table, never events.

From the deployed DAD root on a compute node::

    python meg/code/preprocessing_basic/trial_preparation.py \
      --config meg/params/pipeline_default.yaml \
      --site-config meg/params/pipeline_cluster.yaml \
      --subject-config meg/params/sub-30.yaml

Use ``--validate-only`` to build/check the accounting entirely in memory.
After all P04-P39 bundles exist, run the same command with
``--summarize-cohort`` (no subject selector is needed). ``--overwrite`` permits
replacement of this stage's explicitly named outputs only; it is never implicit.

Inputs, outputs, and assumptions
-------------------------------
Stage 1 evidence -> per-recording boundaries and condition intervals -> shared
trial selection -> exact accepted exclusions -> coverage and eligibility QC.
Outputs are new ``trial_preparation/`` bundles under each subject's preprocessing
and reports directories; the cohort report is under reports/group. The source
manifest, Stage 1 tables, raw files, and authoritative decisions are read-only.

Eligibility is provisional under the recorded configuration. Physiological
artifact rejection, participant inclusion, and visual-timing shifts have not
been evaluated. Session notes are contextual evidence, not automatic exclusions.
Original absolute sample coordinates include ``raw.first_samp``; sample counts
alone are not recording-end coordinates. Stage 04 consumes the same ledger and
checks that its scientific policy, Stage 1 inputs, and output hashes still agree.

Report checklists request acceptance of the new trial accounting and verification
that accepted Stage 1 decisions were applied. They reserve renewed scientific
review for changed inputs or new contradictory evidence. Previously published
bundles retain their recorded report text and hashes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping

import mne

from downstream_stage_utils import build_trial_decisions, csv_bool, load_trigger_event_rows
from pipeline_config import ANALYSIS_SUBJECT_IDS, load_and_resolve_config
from preprocessing_stage_utils import discover_raw_run_records, guarded_path, utc_now


SCHEMA_VERSION = "1.0"
TAIL_REASON = "recording_tail_unavailable"
RESOLVED_REPORTED_ONLY = {
    "RESOLVED_RECORDING_TAIL_UNAVAILABLE", "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK",
}


# -----------------------------------------------------------------------------
# Canonical paths, portable fingerprints, and guarded output transactions
# -----------------------------------------------------------------------------


def stage2_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Resolve a new dedicated bundle without touching Stage 1 paths."""
    paths = cfg["resolved_paths"]
    pre = Path(paths["derivatives_preprocessing_dir"]) / "trial_preparation"
    reports = Path(paths["derivatives_reports_dir"]) / "trial_preparation"
    return {"ledger": pre / "stage2_trial_decisions.csv",
            "coverage": pre / "stage2_recording_coverage.csv",
            "unavailable": pre / "stage2_unavailable_behavioral_trials.csv",
            "counts": reports / "stage2_trial_counts.csv",
            "identities": reports / "stage2_identity_counts.csv",
            "markdown": reports / "stage2_trial_preparation.md",
            "report": reports / "stage2_trial_preparation.json"}


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a small evidence table, retaining blank and original values."""
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def digest(value: Any) -> str:
    """Hash a JSON-compatible scientific contract independently of mount paths."""
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def file_hash(path: Path) -> str:
    """Hash small evidence/output files; never hash raw MEG signal files."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def policy_contract(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """Select the settings that determine trial accounting, excluding later QC."""
    epoch = cfg["epoching"]
    return {"subject": cfg["layout"]["subject_bids_id"],
            "session": cfg["layout"]["session_bids_id"],
            "conditions": epoch["event_id_map"],
            "identities": cfg["subject_profile"]["identity_range"],
            "trial_exclusion": epoch["trial_exclusion"],
            "accepted_exclusions": epoch.get("accepted_trial_exclusions", []),
            "epoch_window_sec": [epoch["tmin_sec"], epoch["tmax_sec"]]}


def evidence_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Resolve accepted tables locally rather than following historical mounts."""
    paths = cfg["resolved_paths"]
    reports = Path(paths["derivatives_reports_dir"])
    events = Path(paths["derivatives_preprocessing_dir"]) / "events"
    subject = cfg["subject"]["id"]
    deep = reports / "02_triggers_extraction_custum_report_deeper"
    return {"manifest": reports / "source_manifest.json",
            "inventory": reports / "dataset_inventory.json",
            "pipeline_run": reports / "pipeline_run.json",
            "events": events / "02_trigger_events.csv",
            "alignment": events / "02_trigger_alignment.csv",
            "warnings": events / "02_transition_warnings.csv",
            "associations": events / "02_trigger_source_associations.json",
            "deep": deep / f"sub{subject}_all_blocks_reported_vs_observed.csv"}


def evidence_fingerprints(paths: Mapping[str, Path]) -> dict[str, str]:
    """Allow later manifest lineage additions while detecting changed evidence.

    The manifest's complete Stage 1 projection is immutable. Later cleaning
    adds only derivative_lineage and must not invalidate this trial accounting.
    """
    result = {key: file_hash(path) for key, path in paths.items() if key != "manifest"}
    manifest = json.loads(paths["manifest"].read_text())
    projection = {key: value for key, value in manifest.items() if key != "derivative_lineage"}
    result["manifest_stage1_projection"] = digest(projection)
    return result


def scientific_fingerprints(paths: Mapping[str, Path]) -> dict[str, str]:
    """Ignore orchestration timestamps/logs that later stages legitimately update.

    Exact trigger tables stay immutable. Source identities, sizes, associations,
    and acquisition headers must also agree, even if a later runner refreshes
    the manifest's generated time, config paths, or derivative lineage.
    """
    result = {key: file_hash(paths[key]) for key in
              ("events", "alignment", "warnings", "associations", "deep")}
    manifest = json.loads(paths["manifest"].read_text())
    result["manifest_sources"] = digest({
        key: [{k: v for k, v in row.items() if k != "source_path"} for row in rows]
        for key, rows in manifest["sources"].items()})
    result["manifest_associations"] = digest([
        {key: value for key, value in row.items() if key != "source_fif_path"}
        for row in manifest["trigger_derived_associations"]["meg_acquisition_run_to_behavioral_segments"]])
    inventory = json.loads(paths["inventory"].read_text())
    result["acquisition_headers"] = digest([
        {key: value for key, value in row.items() if key != "path"}
        for row in inventory["files"]["meg_fif"]])
    return result


def csv_text(rows: list[dict[str, Any]], empty_fields: tuple[str, ...] = ()) -> str:
    """Serialize a review table with headers even when it has no rows."""
    fields = list(dict.fromkeys(key for row in rows for key in row)) or list(empty_fields)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def write_bundle(paths: Mapping[str, Path], texts: dict[str, str], report: dict[str, Any],
                 *, raw_root: str, overwrite: bool) -> None:
    """Precheck every destination, then publish the completion report last.

    Refuse any existing output by default, including a partially written bundle.
    Exclusive creation also prevents two normal jobs from silently overwriting
    each other. Explicit overwrite is confined to these enumerated Stage 2 files.
    """
    for path in paths.values():
        guarded_path(path, raw_root)
        if path.is_symlink() or path.parent.is_symlink():
            raise ValueError(f"Refusing a symlink in a Stage 2 output bundle: {path}")
        if path.exists() and not overwrite:
            raise FileExistsError(f"Stage 2 output already exists; no files changed: {path}")
    report["output_sha256"] = {
        key: hashlib.sha256(value.encode()).hexdigest() for key, value in texts.items()
    }
    texts["report"] = json.dumps(report, indent=2, sort_keys=True) + "\n"
    for key, value in texts.items():
        path = paths[key]
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w" if overwrite else "x", encoding="utf-8", newline="") as handle:
            handle.write(value)


# -----------------------------------------------------------------------------
# Accepted evidence validation and per-trial accounting
# -----------------------------------------------------------------------------


def validate_stage1(cfg: Mapping[str, Any], paths: Mapping[str, Path]) -> dict[str, Any]:
    """Require complete accepted inputs and no newly unresolved Stage 1 cases."""
    inputs = {key: json.loads(paths[key].read_text())
              for key in ("manifest", "inventory", "pipeline_run")}
    subject = cfg["layout"]["subject_bids_id"]
    if (inputs["inventory"].get("status") != "ok"
        or inputs["pipeline_run"].get("status") != "completed"
        or inputs["manifest"].get("warnings")
        or inputs["inventory"].get("validation", {}).get("errors")
        or inputs["inventory"].get("validation", {}).get("warnings")):
        raise RuntimeError("Stage 1 completion/validation no longer matches the accepted gate.")
    for name in ("inventory", "manifest"):
        if inputs[name].get("subject_bids_id") != subject:
            raise RuntimeError(f"Subject mismatch in Stage 1 {name}.")
        for key in ("condition_triggers", "identity_range"):
            if inputs[name]["subject_profile"][key] != cfg["subject_profile"][key]:
                raise RuntimeError(f"Stage 1 profile mismatch: {key}")
    if any(csv_bool(row["unresolved_transition"]) for row in read_csv(paths["warnings"])):
        raise RuntimeError("New unresolved transition evidence requires Stage 1 review.")
    alignment = read_csv(paths["alignment"])
    deep = {int(row["align_row"]): row for row in read_csv(paths["deep"])}
    for row in alignment:
        match = deep.get(int(row["alignment_index"]))
        if not match or match["status"] != row["status"]:
            raise RuntimeError("Canonical alignment and accepted deep report disagree.")
        if row["status"] == "REPORTED_ONLY" and match["stage1_adjudication"] not in RESOLVED_REPORTED_ONLY:
            raise RuntimeError("A reported-only case requires fresh Stage 1 review.")
    inputs["alignment"] = alignment
    inputs["deep"] = deep
    return inputs


def context_key(row: Mapping[str, Any]) -> tuple[str, ...]:
    """Match a behavioral trial including source CSV and replay identity."""
    return tuple(str(row.get(key, "")) for key in
                 ("behavioral_block", "behavioral_run", "behavioral_trial")) + (
        Path(str(row.get("reported_source_csv", ""))).name,
        str(row.get("reported_replay_status", row.get("replay_status", ""))),
    )


def summarize_counts(ledger: list[dict[str, Any]], cfg: Mapping[str, Any]) -> tuple[list, list]:
    """Count each eligibility state and retain zero eligible identity cells."""
    grouped: dict[tuple, Counter] = defaultdict(Counter)
    for row in ledger:
        key = tuple(row.get(name, "") for name in
                    ("source_fif", "acquisition_run", "behavioral_block", "behavioral_run", "replay_status", "condition_name"))
        grouped[key][row["eligibility_status"]] += 1
    counts = []
    subject = cfg["layout"]["subject_bids_id"]
    for key, values in sorted(grouped.items(), key=lambda item: str(item[0])):
        counts.append(dict(subject_bids_id=subject,
                           **dict(zip(("source_fif", "acquisition_run", "behavioral_block", "behavioral_run", "replay_status", "condition_name"), key)),
                           recorded_candidates=sum(values.values()),
                           eligible=values["eligible_under_current_policy"],
                           excluded=values["excluded"], review_required=values["review_required"]))
    low, high = cfg["subject_profile"]["identity_range"]
    identities = []
    for condition in cfg["epoching"]["event_id_map"]:
        for identity in range(low, high + 1):
            rows = [r for r in ledger if r["condition_name"] == condition
                    and r["trial_sequence_identity"] == identity]
            identities.append({"subject_bids_id": subject, "condition_name": condition,
                               "sequence_identity": identity, "recorded_candidates": len(rows),
                               "eligible": sum(r["included_by_trial_policy"] for r in rows),
                               "excluded": sum(r["eligibility_status"] == "excluded" for r in rows),
                               "review_required": sum(r["eligibility_status"] == "review_required" for r in rows)})
    return counts, identities


def prepare_subject(cfg: Mapping[str, Any], *, overwrite: bool = False,
                    validate_only: bool = False) -> dict[str, Any]:
    """Build the complete subject ledger; optionally publish only new outputs.

    No signal data are loaded. Header sizes and sample coordinates are checked
    against accepted inventory/manifest records. Unexpected evidence is never
    converted to an automatic scientific acceptance.
    """
    if cfg["subject"]["id"] not in ANALYSIS_SUBJECT_IDS:
        raise ValueError("Stage 2 cohort is P04-P39; pilot participants are excluded.")
    paths, outputs = evidence_paths(cfg), stage2_paths(cfg)
    fingerprints = evidence_fingerprints(paths)
    inputs = validate_stage1(cfg, paths)
    authority = Path(cfg["project"]["docs_root"]) / "PREPROCESSING_DECISIONS.md"
    if not authority.is_file():
        # The local checkout deliberately has no competing decisions copy.
        authority = Path(cfg["resolved_paths"]["data_root"]).parent / "DAD/meg/docs/PREPROCESSING_DECISIONS.md"
    authority_text = authority.read_text(encoding="utf-8")
    if "Accepted on 2026-10-19 for Stage 2 trial preparation" not in authority_text:
        raise RuntimeError("The authoritative Stage 1 acceptance record was not found.")
    events = load_trigger_event_rows(paths["events"])
    # The two reviewed adjacent-pulse exceptions remain two observed events.
    # Do not let a later canonical-table change silently undo their gaze breaks.
    for override in cfg.get("trigger_parser", {}).get("reviewed_keep_both_transitions", []):
        pair = []
        for prefix in ("first", "second"):
            found = [row for row in events if row["source_fif"] == override["source_fif"]
                     and int(row["sample_start_raw"]) == override[f"{prefix}_sample_raw"]
                     and int(row["trigger_code"]) == override[f"{prefix}_code"]]
            if (len(found) != 1 or found[0]["reported_alignment_status"] != "MATCH_WARNING"
                or found[0]["transition_decision_source"] != "reviewed_keep_both"):
                raise RuntimeError("A reviewed keep-both transition is missing or changed.")
            pair.append(found[0])
        if context_key(pair[0]) != context_key(pair[1]):
            raise RuntimeError("Reviewed keep-both events no longer share one trial context.")
    records = discover_raw_run_records(cfg)
    source_names = {record.source_filename for record in records}
    manifest_records = inputs["manifest"]["sources"]["task_meg_acquisition_runs"]
    if source_names != {record["source_filename"] for record in manifest_records}:
        raise RuntimeError("Raw source set differs from the accepted manifest.")
    if {row["source_fif"] for row in events} != source_names:
        raise RuntimeError("Canonical event sources differ from the accepted recording set.")
    associations = inputs["manifest"]["trigger_derived_associations"]["meg_acquisition_run_to_behavioral_segments"]
    if ({item["source_fif"] for item in associations} != source_names
        or any(item["association_status"] != "aligned" for item in associations)):
        raise RuntimeError("Accepted event-aligned source associations are incomplete.")
    accepted = cfg["epoching"].get("accepted_trial_exclusions", [])
    if any(item["source_fif"] not in source_names for item in accepted):
        raise RuntimeError("An accepted exclusion source is missing; refusing silent omission.")
    headers = {row["name"]: row for row in inputs["inventory"]["files"]["meg_fif"]}
    expected_sizes = {row["source_filename"]: row["size_bytes"] for row in manifest_records}
    alignment_by_event = {str(row["observed_event_index"]): row for row in inputs["alignment"]
                          if row.get("observed_event_index")}
    tails: dict[tuple, list] = defaultdict(list)
    reported_only: dict[tuple, list] = defaultdict(list)
    for row in inputs["alignment"]:
        if row["status"] == "REPORTED_ONLY":
            reported_only[context_key(row)].append(row)
        if row.get("reported_only_reason") == TAIL_REASON:
            if (row["status"] != "REPORTED_ONLY" or row.get("raw_sample_frame")
                or int(row["constrained_acquisition_run"]) not in
                    cfg["subject_profile"].get("truncated_acquisition_runs", [])):
                raise RuntimeError("Unavailable-tail evidence lacks the reviewed recording constraint.")
            tails[context_key(row)].append(row)

    ledger, coverage = [], []
    low, high = cfg["subject_profile"]["identity_range"]
    for record in records:
        source = record.source_filename
        run_events = [row for row in events if row["source_fif"] == source]
        # first_samp is essential: P37's cutoff is not n_samples in absolute
        # coordinates. Header-only access leaves all MEG samples on disk.
        raw = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
        try:
            first, end = int(raw.first_samp), int(raw.first_samp + raw.n_times)
            sfreq = float(raw.info["sfreq"])
            online_highpass, online_lowpass = float(raw.info["highpass"]), float(raw.info["lowpass"])
            n_samples = int(raw.n_times)
        finally:
            raw.close()
        if (n_samples != headers[source]["n_samples"] or sfreq != headers[source]["sfreq_hz"]
            or Path(record.source_path).stat().st_size != expected_sizes[source]
            or online_highpass != headers[source]["online_highpass_hz"]
            or online_lowpass != headers[source]["online_lowpass_hz"]):
            raise RuntimeError(f"Recording header/size changed since Stage 1: {source}")
        if any(not first <= int(row["sample_start_raw"]) < int(row["sample_end_raw"]) <= end
               for row in run_events):
            raise RuntimeError(f"Canonical events lie outside recording coordinates: {source}")
        decisions = build_trial_decisions(
            run_events, source_fif=source, run_end_sample_exclusive=end,
            event_id_map=cfg["epoching"]["event_id_map"], identity_range=range(low, high + 1),
            exclusion_policy=cfg["epoching"]["trial_exclusion"],
            accepted_exclusions=accepted, strict=False,
        )
        owned_events = set()
        for row in decisions:
            row["subject_bids_id"] = cfg["layout"]["subject_bids_id"]
            row["session_bids_id"] = cfg["layout"]["session_bids_id"]
            row["trial_key"] = f"{row['subject_bids_id']}:{source}:{row['raw_sample_frame']}"
            row["recording_first_sample_raw"] = first
            row["recording_end_sample_raw_exclusive"] = end
            row["sfreq_hz"] = sfreq
            row["online_highpass_hz"] = online_highpass
            row["online_lowpass_hz"] = online_lowpass
            row["timing_basis"] = "original_digital_STI101_unshifted"
            row["visual_timing_review"] = "pending_cohort_gate"
            indices = row["evidence_event_indices"].split(",")
            owned_events.update(indices)
            row["evidence_alignment_rows"] = ",".join(
                str(alignment_by_event[index]["alignment_index"])
                for index in indices if index in alignment_by_event)
            on_alignment = alignment_by_event.get(row["condition_event_index"], {})
            log_only_rows = reported_only.get(context_key(on_alignment), [])
            row["reported_only_alignment_rows"] = ",".join(r["alignment_index"] for r in log_only_rows)
            adjudications = [inputs["deep"][int(r["alignment_index"])] for r in log_only_rows]
            row["stage1_adjudications"] = ",".join(a["stage1_adjudication"] for a in adjudications)
            row["adjudication_partner_alignment_rows"] = ",".join(a["adjudication_partner_align_row"] for a in adjudications)
            missing = tails.get(context_key(on_alignment), [])
            row["unavailable_tail_alignment_rows"] = ",".join(r["alignment_index"] for r in missing)
            row["recording_tail_intersects_trial"] = bool(missing)
            matched_rule = next((item for item in accepted if item["source_fif"] == source
                                 and item["raw_sample_frame"] == row["raw_sample_frame"]), None)
            if matched_rule and matched_rule.get("require_recording_tail_evidence") and not missing:
                raise RuntimeError("Reviewed incomplete-trial exclusion lost its unavailable-tail evidence.")
            reasons = row["exclusion_reasons"].split(",") if row["exclusion_reasons"] else []
            reviews = row["review_reasons"].split(",") if row["review_reasons"] else []
            if (row["reported_alignment_status"] not in {"MATCH", "MATCH_WARNING"}
                and not row["accepted_exclusion_reason"]):
                reviews.append("unmatched_condition_onset")
            if any(a["stage1_adjudication"] == "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK" for a in adjudications):
                if "gaze_break" not in reasons:
                    reviews.append("gaze_adjudication_requires_excluded_trial")
            if missing and not row["accepted_exclusion_reason"]:
                reviews.append("unreviewed_recording_tail_trial")
            # This is digital-time coverage of the configured future epoch
            # window, not proof of full stimulus presentation or final timing.
            epoch_first = row["raw_sample_frame"] + round(cfg["epoching"]["tmin_sec"] * sfreq)
            epoch_last = row["raw_sample_frame"] + round(cfg["epoching"]["tmax_sec"] * sfreq)
            row["digital_epoch_window_inside_recording"] = first <= epoch_first and epoch_last < end
            if not row["digital_epoch_window_inside_recording"]:
                reasons.append("epoch_window_outside_recording")
            row["exclusion_reasons"] = ",".join(reasons)
            row["review_reasons"] = ",".join(reviews)
            row["included_by_trial_policy"] = not reasons and not reviews
            row["eligibility_status"] = ("review_required" if reviews else
                                         "excluded" if reasons else "eligible_under_current_policy")
            row["final_epoch_inclusion"] = "not_evaluated"
        orphan = [row for row in run_events if row["event_index"] not in owned_events]
        coverage.append({"subject_bids_id": cfg["layout"]["subject_bids_id"], "source_fif": source,
                         "acquisition_run": record.acquisition_run_label, "first_sample_raw": first,
                         "end_sample_raw_exclusive": end, "n_samples": n_samples, "sfreq_hz": sfreq,
                         "online_highpass_hz": online_highpass, "online_lowpass_hz": online_lowpass,
                         "source_size_bytes": expected_sizes[source], "recorded_events": len(run_events),
                         "recorded_candidates": len(decisions), "unowned_event_indices": ",".join(r["event_index"] for r in orphan),
                         "unowned_event_codes": ",".join(r["trigger_code"] for r in orphan),
                         "unowned_samples_raw": ",".join(r["sample_start_raw"] for r in orphan),
                         "behavioral_segments": json.dumps(next(a["behavioral_segments"] for a in associations if a["source_fif"] == source)),
                         "reviewed_truncated_recording": record.acquisition_run_label in cfg["subject_profile"].get("truncated_acquisition_runs", [])})
        ledger.extend(decisions)

    # Keep log-only trial contexts outside the scientific candidate ledger.
    unavailable = []
    for key, rows in sorted(tails.items()):
        run = int(rows[0]["constrained_acquisition_run"])
        source = next(record.source_filename for record in records if record.acquisition_run_label == run)
        recorded = [r for r in ledger if r["source_fif"] == source and
                    (str(r["behavioral_block"]), str(r["behavioral_run"]), str(r["behavioral_trial"])) == key[:3]
                    and Path(r["reported_source_csv"]).name == key[3] and r["replay_status"] == key[4]]
        unavailable.append({"subject_bids_id": cfg["layout"]["subject_bids_id"], "constrained_source_fif": source,
                            "acquisition_run": run, "behavioral_block": key[0], "behavioral_run": key[1],
                            "behavioral_trial": key[2], "reported_source_csv": key[3], "replay_status": key[4],
                            "raw_condition_onset_sample": "", "recorded_candidate_present": bool(recorded),
                            "alignment_rows": ",".join(r["alignment_index"] for r in rows),
                            "unavailable_logged_events": len(rows), "disposition": "stage1_recording_tail_unavailable"})
    if len({row["trial_key"] for row in ledger}) != len(ledger):
        raise RuntimeError("Trial keys are not unique.")
    if sum(bool(row["accepted_exclusion_reason"]) for row in ledger) != len(accepted):
        raise RuntimeError("Not every accepted source/sample exclusion matched exactly once.")
    # Detect changes during the run before publishing a report based on mixed inputs.
    if evidence_fingerprints(paths) != fingerprints:
        raise RuntimeError("Stage 1 evidence changed during Stage 2 preparation.")
    counts, identities = summarize_counts(ledger, cfg)
    status_counts = dict(Counter(row["eligibility_status"] for row in ledger))
    reason_counts = dict(Counter(reason for row in ledger for reason in row["exclusion_reasons"].split(",") if reason))
    report = {"schema_version": SCHEMA_VERSION, "status": "completed", "scientific_acceptance": "pending_independent_review",
              "generated_utc": utc_now(), "subject_bids_id": cfg["layout"]["subject_bids_id"],
              "session_bids_id": cfg["layout"]["session_bids_id"], "policy_contract": policy_contract(cfg),
              "stage1_input_sha256": fingerprints, "decisions_record": str(authority),
              "scientific_input_sha256": scientific_fingerprints(paths),
              "decisions_record_sha256": file_hash(authority), "implementation_sha256": file_hash(Path(__file__)),
              "shared_trial_builder_sha256": file_hash(Path(__file__).with_name("downstream_stage_utils.py")),
              "config_sources": cfg.get("_config_sources", {}), "n_recordings": len(records),
              "n_recorded_candidates": len(ledger), "eligibility_counts": status_counts,
              "exclusion_reason_counts": reason_counts, "accepted_exclusions_matched": len(accepted),
              "n_unavailable_logged_events": sum(len(rows) for rows in tails.values()),
              "n_unavailable_trial_contexts": len(unavailable),
              "zero_eligible_condition_identity_cells": [row for row in identities if not row["eligible"]],
              "session_review_clues": inputs["inventory"].get("session_review_clues", {}),
              "known_acquisition_exceptions": cfg["subject_profile"].get("known_exceptions", []),
              "recording_coverage": coverage, "outputs": {key: str(path) for key, path in outputs.items()}}
    markdown = render_subject_report(report, counts)
    if not validate_only:
        write_bundle(outputs, {"ledger": csv_text(ledger), "coverage": csv_text(coverage),
                               "unavailable": csv_text(unavailable, ("subject_bids_id", "constrained_source_fif", "alignment_rows", "disposition")),
                               "counts": csv_text(counts), "identities": csv_text(identities),
                               "markdown": markdown}, report,
                     raw_root=cfg["resolved_paths"]["data_root"], overwrite=overwrite)
    return report


# -----------------------------------------------------------------------------
# Human reports, cohort aggregation, and the future epoching handoff
# -----------------------------------------------------------------------------


def render_subject_report(report: Mapping[str, Any], counts: list) -> str:
    """Render saved accounting and new review tasks from report/count mappings.

    Return Markdown without reading evidence or changing decisions. Accepted
    Stage 1 dispositions are checked for correct ledger application; only new
    accounting issues or contradictory evidence require renewed review.
    """
    lines = [f"# Stage 2 trial preparation: {report['subject_bids_id']}", "",
             "Scientific acceptance: pending independent review. Original digital timing; no epochs created.", "",
             f"Recorded candidates: {report['n_recorded_candidates']}; recordings: {report['n_recordings']}.",
             f"Eligibility: `{json.dumps(report['eligibility_counts'], sort_keys=True)}`.",
             f"Reasons (may overlap): `{json.dumps(report['exclusion_reason_counts'], sort_keys=True)}`.",
             f"Accepted source/sample exclusions matched: {report['accepted_exclusions_matched']}.",
             f"Unavailable logged events: {report['n_unavailable_logged_events']} in {report['n_unavailable_trial_contexts']} behavioral contexts.", "",
             "## Counts by recording, block, and condition", "",
             "| FIF | Block/run | Replay | Condition | Recorded | Eligible | Excluded | Review |",
             "| --- | --- | --- | --- | ---: | ---: | ---: | ---: |"]
    for row in counts:
        lines.append(f"| {row['source_fif']} | {row['behavioral_block']}/{row['behavioral_run']} | {row['replay_status']} | {row['condition_name']} | {row['recorded_candidates']} | {row['eligible']} | {row['excluded']} | {row['review_required']} |")
    lines.extend(["", "## Decision before the next command", "",
                  "- Accept the new trial accounting: reconcile ledger rows and count tables by source, condition, identity, replay status, and eligibility.",
                  "- Review selection effects under the recorded policy, including missing identity coverage and future epoch-window exclusions; reason counts may overlap.",
                  "- Confirm accepted Stage 1 rules were applied using automatic checks and targeted ledger spot-checks. Reopen a decision only for changed inputs or new contradictory evidence.",
                  "- Investigate new review_required rows and any independently found count, identity, lineage, or policy-application inconsistency, including when no rows are flagged.",
                  "- Eligible means candidate under the recorded policy; final participant/epoch inclusion, visual timing, and physiological rejection remain pending.",
                  "- Session notes below are contextual clues, not automatic exclusions.", "",
                  "## Session and acquisition context", "", "```json",
                  json.dumps({"known_exceptions": report["known_acquisition_exceptions"],
                              "session_review_clues": report["session_review_clues"]}, indent=2), "```", ""])
    return "\n".join(lines)


def load_prepared_decisions(cfg: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Load a completed, unchanged Stage 2 ledger for production epoching.

    Reject stale policy/evidence, edited or partial output bundles, and any new
    review-required rows. Independent scientific acceptance remains a human gate.
    """
    paths = stage2_paths(cfg)
    report = json.loads(paths["report"].read_text())
    if report.get("schema_version") != SCHEMA_VERSION or report.get("status") != "completed":
        raise RuntimeError("Stage 2 bundle is incomplete or uses an unsupported schema.")
    if report["policy_contract"] != policy_contract(cfg):
        raise RuntimeError("Trial policy changed; regenerate and review Stage 2 explicitly.")
    if report["scientific_input_sha256"] != scientific_fingerprints(evidence_paths(cfg)):
        raise RuntimeError("Stage 1 evidence changed after Stage 2; review is required.")
    if set(report["output_sha256"]) != set(paths).difference({"report"}):
        raise RuntimeError("Stage 2 completion record lacks the required outputs.")
    for key, expected in report["output_sha256"].items():
        if file_hash(paths[key]) != expected:
            raise RuntimeError(f"Stage 2 output changed or is incomplete: {key}")
    rows = read_csv(paths["ledger"])
    if any(row["eligibility_status"] == "review_required" for row in rows):
        raise RuntimeError("Trial ledger still contains review-required candidates.")
    # CSV values need their original scalar types for MNE metadata and for
    # _mark_mne_outcomes, which indexes rows by absolute raw sample frame.
    integer_fields = {"acquisition_run", "raw_sample_frame", "sample_start_raw", "trial_interval_end_raw_exclusive",
                      "condition_code", "trial_sequence_identity", "identity_sample_start_raw",
                      "behavioral_block", "behavioral_run", "behavioral_trial", "unresolved_transition_count"}
    boolean_fields = {"included_by_trial_policy", "is_replay_trial", "is_transition_warning",
                      "gaze_break_evidence", "catch_type_1_evidence", "catch_type_2_evidence", "esc_termination_evidence"}
    for row in rows:
        for key in integer_fields:
            row[key] = int(row[key]) if row[key] else None
        for key in boolean_fields:
            row[key] = csv_bool(row[key])
    if len(rows) != report["n_recorded_candidates"]:
        raise RuntimeError("Stage 2 candidate count disagrees with the completion report.")
    return rows


def summarize_cohort(cfg: Mapping[str, Any], *, overwrite: bool = False,
                     validate_only: bool = False) -> dict[str, Any]:
    """Validate all 36 bundles and render the cohort accounting review summary.

    Read bundles under the configured derivative root and return the combined
    report. Write dedicated summary outputs unless validate_only is true;
    overwrite must be explicit. The checklist reviews new accounting results
    while preserving accepted Stage 1 dispositions.
    """
    root = Path(cfg["resolved_paths"]["derivatives_root"])
    reports, counts, identities = [], [], []
    for subject in ANALYSIS_SUBJECT_IDS:
        sub_cfg = dict(cfg, resolved_paths=dict(cfg["resolved_paths"]))
        # Only paths are needed here. No other participant's scientific policy
        # is inferred from the selector used to locate the derivative root.
        sub_cfg["resolved_paths"].update(
            derivatives_preprocessing_dir=str(root / "preprocessing" / f"sub-{subject}" / "ses-01"),
            derivatives_reports_dir=str(root / "reports" / f"sub-{subject}" / "ses-01"))
        paths = stage2_paths(sub_cfg)
        report = json.loads(paths["report"].read_text())
        if (report.get("status") != "completed" or report.get("schema_version") != SCHEMA_VERSION
            or report["subject_bids_id"] != f"sub-{subject}" or report["session_bids_id"] != "ses-01"):
            raise RuntimeError(f"Incomplete/mismatched Stage 2 subject bundle: {subject}")
        if set(report["output_sha256"]) != set(paths).difference({"report"}):
            raise RuntimeError(f"Incomplete Stage 2 completion record: {subject}")
        sources = cfg["_config_sources"]
        subject_cfg, warnings, _, _ = load_and_resolve_config(
            sources["default"], Path(cfg["project"]["params_root"]) / f"sub-{subject}.yaml",
            site_config_path=sources.get("site"), required_sections=("epoching",))
        if warnings or report["policy_contract"] != policy_contract(subject_cfg):
            raise RuntimeError(f"Current trial policy differs from subject {subject} bundle.")
        for key, expected in report["output_sha256"].items():
            if file_hash(paths[key]) != expected:
                raise RuntimeError(f"Incomplete/changed Stage 2 output: {subject}/{key}")
        if report["scientific_input_sha256"] != scientific_fingerprints(evidence_paths(dict(
            sub_cfg, subject={"id": subject}))):
            raise RuntimeError(f"Stage 1 inputs changed after subject {subject} preparation.")
        reports.append(report)
        counts.extend(read_csv(paths["counts"]))
        identities.extend(read_csv(paths["identities"]))
    totals = Counter()
    for report in reports:
        totals.update(report["eligibility_counts"])
    report = {"schema_version": SCHEMA_VERSION, "status": "completed", "scientific_acceptance": "pending_independent_review",
              "generated_utc": utc_now(), "n_subjects": len(reports),
              "eligibility_counts": dict(totals), "participants": reports}
    lines = ["# Stage 2 cohort trial preparation", "", "P04-P39; independent scientific review pending.", "",
             f"Totals: `{json.dumps(dict(totals), sort_keys=True)}`.", "",
             "| Participant | Recordings | Candidates | Eligible | Excluded | Review | Accepted exclusions | Unavailable logged events |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for item in reports:
        values = item["eligibility_counts"]
        lines.append(f"| {item['subject_bids_id']} | {item['n_recordings']} | {item['n_recorded_candidates']} | {values.get('eligible_under_current_policy', 0)} | {values.get('excluded', 0)} | {values.get('review_required', 0)} | {item['accepted_exclusions_matched']} | {item['n_unavailable_logged_events']} |")
    lines.extend(["", "## Decision before the next command", "",
                  "Accept the new trial accounting: reconcile ledger/count tables, condition/identity/replay coverage, selection reasons, and future epoch-window effects. Confirm accepted Stage 1 rules were applied through automatic checks and targeted ledger spot-checks; reopen a decision only for changed inputs or new contradictory evidence. Resolve new accounting inconsistencies even when no review-required rows are flagged. Timing, cleaning, and final inclusion decisions remain pending.", ""])
    out = root / "reports" / "group" / "trial_preparation"
    paths = {"report": out / "stage2_cohort_summary.json", "markdown": out / "stage2_cohort_summary.md",
             "counts": out / "stage2_trial_counts.csv", "identities": out / "stage2_identity_counts.csv"}
    if not validate_only:
        write_bundle(paths, {"markdown": "\n".join(lines), "counts": csv_text(counts),
                             "identities": csv_text(identities)}, report,
                     raw_root=cfg["resolved_paths"]["data_root"], overwrite=overwrite)
    return report


def main() -> int:
    """Resolve configuration and execute one bounded Stage 2 operation."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config")
    parser.add_argument("--subject-config")
    parser.add_argument("--summarize-cohort", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    cfg, warnings, _, _ = load_and_resolve_config(args.config, args.subject_config,
        site_config_path=args.site_config, required_sections=("epoching",))
    if warnings:
        raise RuntimeError(f"Configuration requires review: {warnings}")
    function = summarize_cohort if args.summarize_cohort else prepare_subject
    report = function(cfg, overwrite=args.overwrite, validate_only=args.validate_only)
    print(json.dumps({key: report[key] for key in ("status", "eligibility_counts", "scientific_acceptance")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
