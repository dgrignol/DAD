#!/usr/bin/env python3
"""Validate a frozen cohort authorization, then record guarded provisional ICA choices.

Usage (activate the model's existing MNE environment first)::

    python record_provisional_ica_choices.py --authorization /path/authorization.json \
        --output /path/NEW-UNUSED-DIRECTORY --jobs 4
    # Add --record to write decisions AFTER the whole eligible table validates.
    # --validate-subject 32 only retries P32 source validation; never writes choices.
    # Cluster validation additionally accepts --site-config meg/params/pipeline_cluster.yaml.

Inputs: adjacent authorization.json/choices.json, frozen complete proposal tables,
current exact models/reports/native annotated headers, operative 4A/4B receipts
and Stage 2/3 diagnostic-retention evidence. Null exclusions remain unresolved;
explicit empty lists are valid only for a completed authorized table. Model
indices never transfer. No classification or user visual review is inferred.

Outputs: a new immutable validation/recording directory with participant progress,
before/after tuning snapshots, ordinary Stage 4C command logs and joined summaries.
All eligible participants are submitted together to bounded separate processes
with one numerical thread per worker. No tuning writes begin until all validation
workers join successfully. Recording uses review_cleaning.py's ordinary guards,
then adds structured provisional provenance. Existing decisions are preserved and
conflicts stop the operation. The script never applies ICA, submits jobs, edits
annotations, creates epochs or grants cleaned-output/scientific acceptance.

Validation or recording failures leave explicit partial evidence and snapshots;
there is no automatic rollback over potentially concurrent manual work. Use a
new output directory for every attempt. Inspect a partial recording before reuse.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import shutil
import traceback

import mne

from audit_ica_task_diagnostics import audit_participant
from continuous_cleaning import file_hash, require_review, stage_configuration
from participant_parallel import map_participants, participant_workers
from pipeline_config import WORKING_SUBJECT_IDS, read_yaml_mapping
from preprocessing_stage_utils import stage_output_dirs, update_tuning_file, utc_now
from review_final_ica import validate_review
import review_cleaning

PROJECT = Path(__file__).resolve().parents[3]
STAGE = "03_C_ICA"
KIND = "user_provisional_component_processing_authorization"


# -----------------------------------------------------------------------------
# Frozen whole-table checks precede source workers and every production write.
# -----------------------------------------------------------------------------

def require(condition, message):
    """Raise a persistent guard (also when Python assertions are disabled)."""
    if not condition:
        raise RuntimeError(message)


def load_authorization(path):
    """Check exact cohort, null/empty semantics and adjacent table agreement."""
    auth = json.loads(path.read_text())
    choices = json.loads(path.with_name("choices.json").read_text())
    rows = auth["participants"]
    require(choices["authorization_file"] == path.name, "Choices reference a different authorization.")
    require(choices["participants"] == rows, "Frozen choices and authorization disagree.")
    subjects = [row["subject"] for row in rows]
    require(subjects == list(map(int, WORKING_SUBJECT_IDS)) == auth["working_cohort"], "Working cohort changed or duplicated.")
    require(auth["status"] == "provisional_component_processing_authorization", "Wrong authorization kind.")
    require(auth["user_visual_review_completed"] is False and auth["final_scientific_acceptance"] is False,
            "This importer requires honest provisional provenance.")
    blocked = [row["subject"] for row in rows if row["exclusions"] is None]
    eligible = [row for row in rows if row["exclusions"] is not None]
    require(blocked == auth["blocked_participants"], "Null choices do not match pending participants.")
    require(len(eligible) == auth["authorized_participant_count"], "Authorized participant count changed.")
    require(sum(len(row["exclusions"]) for row in eligible) == auth["authorized_exclusion_count"],
            "Authorized exclusion count changed.")
    for row in eligible:
        values = row["exclusions"]
        require(values == sorted(set(values)) and all(type(i) is int and 0 <= i < row["component_count"] for i in values),
                f"P{row['subject']:02}: invalid or repeated component index.")
        require(not set(values).intersection(row["deferred_components_retained"]), "A deferred component would be removed.")
    return auth, eligible


def validate_proposal(row, proposal):
    """Require complete zero-based coverage and exact exclude/defer agreement."""
    require(int(str(proposal["subject"]).lstrip("Pp")) == row["subject"] and
            proposal.get("model_sha256", proposal.get("current_model_sha256")) == row["model_sha256"],
            "Proposal belongs to another participant/model.")
    # The frozen batches predate one common schema. Normalize only their explicit
    # field aliases; unknown/missing fields still fail, never imply keep/empty.
    tables = [proposal[key] for key in ("rows", "evidence", "proposals", "components")
              if isinstance(proposal.get(key), list)]
    require(len(tables) == 1, "Expected exactly one component proposal table.")
    items = [{"component": item.get("component", item.get("component_zero_based")),
              "proposed": item.get("proposed", item.get("proposal", item.get("proposed_action")))}
             for item in tables[0]]
    require(sorted(item["component"] for item in items) == list(range(row["component_count"])),
            "Proposal does not cover every exact-model component once.")
    require(all(item["proposed"] in ("keep", "exclude", "defer") for item in items), "Unknown component choice.")
    require(sorted(item["component"] for item in items if item["proposed"] == "exclude") == row["exclusions"],
            "Proposal exclusions differ from frozen authorization.")
    require(sorted(item["component"] for item in items if item["proposed"] == "defer") == row["deferred_components_retained"],
            "Deferred components differ from frozen authorization.")
    if "proposed_exclusions" in proposal:
        require(proposal["proposed_exclusions"] == row["exclusions"], "Proposal summary differs from its rows.")
    if "all_components_reviewed" in proposal:
        require(proposal["all_components_reviewed"] == row["component_count"], "Proposal is incomplete.")


# -----------------------------------------------------------------------------
# Independent read-only participant checks include native headers and ordinary
# prerequisites. A caught failure is evidence, never an accepted empty choice.
# -----------------------------------------------------------------------------

def validate_participant(task):
    """Return a small auditable result and preserve failures without aborting peers."""
    args, row, readiness, output = task
    subject = row["subject"]
    try:
        cfg, report, records, inputs, model, metadata = validate_review(args, subject)
        root = Path(cfg["resolved_paths"]["derivatives_root"])
        dirs = stage_output_dirs(cfg)
        report_path = dirs["reports"] / "03_C_ICA_review_report.json"
        require(file_hash(report_path) == row["current_report_sha256"], "Frozen report identity changed.")
        require(report["evidence_sha256"] == row["current_report_evidence_sha256"], "Frozen report binding changed.")
        require(report["prefit_model_sha256"] == row["model_sha256"] == readiness["model_sha256"], "Frozen model changed.")
        require(report["n_components_fitted"] == row["component_count"] == readiness["components"], "Component count changed.")
        fit = report["fit_execution"]
        require(fit["n_samples_fitted"] == fit["mapped_fit_samples"] and not fit["iteration_limit_reached"],
                "Fit samples/convergence are not ready.")
        require(report["component_review_coverage"]["complete"], "Component evidence is incomplete.")
        trial = copy.deepcopy(cfg)
        trial.setdefault("stage4_review", {})[STAGE] = {"status": "accepted", "evidence_sha256": report["evidence_sha256"],
            "configuration_sha256": stage_configuration(cfg, STAGE)}
        require_review(trial, STAGE, dirs["reports"])
        tuning = dirs["reports"] / "03_tuning_choices.yaml"
        saved = read_yaml_mapping(tuning)
        receipt = saved.get("stage4_review", {}).get(STAGE)
        existing = saved.get("ica", {})
        if receipt or existing.get("manual_exclude") is not None or existing.get("reviewed_model_sha256"):
            require(receipt and receipt.get("acceptance_kind") == KIND and
                    existing.get("manual_exclude") == row["exclusions"] and
                    existing.get("reviewed_model_sha256") == row["model_sha256"] and
                    receipt.get("authorization_sha256") == args.authorization_sha256,
                    "Existing component decision conflicts or has different provenance; preserved.")
            require_review(cfg, STAGE, dirs["reports"])
        proposal_path = root / row["proposal_evidence"]
        require(file_hash(proposal_path) == row["proposal_evidence_sha256"], "Frozen proposal changed.")
        if row["exclusions"] is not None:
            validate_proposal(row, json.loads(proposal_path.read_text()))
        plan_path = root / "stage4b_adoption/current-masks-20261109-v1/plan.json"
        plan = json.loads(plan_path.read_text())
        policy = json.loads(plan_path.with_name("policy_acceptance.json").read_text())
        require(plan["plan_sha256"] == policy["adoption_plan_sha256"], "Adoption plan/receipt changed.")
        retained = audit_participant((root, next(x for x in plan["participants"] if x["subject"] == subject),
            next(x for x in policy["participants"] if x["subject"] == subject)))
        require(retained["retained_trial_keys_sha256"] == readiness["retained_trial_keys_sha256"], "Trial identities changed.")
        # A read-only import must not weaken the production runtime/cache guard.
        require(mne.__version__ == report["fingerprint"]["fit"]["mne_version"], "Use the saved model's MNE version.")
        if subject == 19:
            record = report["convergence_reconciliation"]
            publication = json.loads((root / "ica_convergence/sub-19" / Path(record["result"]).parent.name / "publication.json").read_text())
            require(publication["published"] and publication["diagnostic_trials_unchanged"] and
                    publication["new_model_sha256"] == row["model_sha256"], "P19 restart publication is not ready.")
        outputs = list(dirs["sss"].glob("*_ica*.fif")) + list(dirs["reports"].glob("03_C_ICA_apply*"))
        outputs += list(dirs["figures"].glob("03_C_ICA_apply*")) + list(dirs["ica_models"].glob("*postreject*"))
        result = {"subject": subject, "ok": True, "model_sha256": row["model_sha256"],
            "exclusions": row["exclusions"], "deferred_components_retained": row["deferred_components_retained"],
            "source_headers_validated": len(inputs), "retention": retained,
            "tuning": str(tuning), "tuning_sha256_before": file_hash(tuning),
            "already_recorded": bool(receipt), "prior_application_outputs": [str(x) for x in outputs],
            "epoching": cfg["epoching"], "preprocessing": cfg["preprocessing"]}
        print(f"P{subject:02}: VALIDATED; {len(inputs)} native headers; exclude={row['exclusions']}", flush=True)
    except Exception as error:
        result = {"subject": subject, "ok": False, "error": str(error), "traceback": traceback.format_exc()}
        print(f"P{subject:02}: VALIDATION FAILED: {error}", flush=True)
    (output / f"P{subject:02}_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


# -----------------------------------------------------------------------------
# A joined validation barrier and complete before snapshots precede recording.
# Each participant uses the ordinary Stage 4C writer; peers have disjoint YAMLs.
# -----------------------------------------------------------------------------

def record_participant(task):
    """Run existing guarded Stage 4C recording and attach truthful machine provenance."""
    args, row, checked, output = task
    subject = row["subject"]
    tuning = Path(checked["tuning"])
    result = {"subject": subject, "ok": False}
    try:
        require(file_hash(tuning) == checked["tuning_sha256_before"], "Tuning changed after cohort validation; preserved.")
        if not checked["already_recorded"]:
            reason = (f"User provisional acceptance on 2026-11-12 WITHOUT completed user visual review. "
                f"Frozen authorization SHA256={args.authorization_sha256}; exact model={row['model_sha256']}; "
                f"retain deferred ICs {row['deferred_components_retained']}. For exploratory processing only; "
                "cleaned-output acceptance and final scientific acceptance remain pending.")
            command = ["--config", str(args.config), "--subject-config", str(args.subject_config_dir / f"sub-{subject:02}.yaml"),
                "--tuning-config", str(tuning), "--stage", "4C", "--reviewer", "Damiano (provisional authorization; no user visual review)",
                "--reason", reason, "--exclude", *map(str, row["exclusions"])]
            if args.site_config:
                command += ["--site-config", str(args.site_config)]
            (output / f"P{subject:02}_record_command.json").write_text(json.dumps(command, indent=2) + "\n")
            review_cleaning.main(command)
            receipt = read_yaml_mapping(tuning)["stage4_review"][STAGE]
            receipt.update(acceptance_kind=KIND, authorization_sha256=args.authorization_sha256,
                authorization_relative_path=args.authorization_relative_path,
                choices_sha256=args.choices_sha256, proposal_evidence_sha256=row["proposal_evidence_sha256"],
                user_visual_review_completed=False, final_scientific_acceptance=False,
                cleaned_output_acceptance=False, deferred_components_retained=row["deferred_components_retained"])
            cfg = __import__("annotation_triage").subject_config(args, subject)
            update_tuning_file(path=tuning, raw_root=cfg["resolved_paths"]["data_root"], subject_id=subject,
                section="stage4_review", values={STAGE: receipt})
        cfg = __import__("annotation_triage").subject_config(args, subject)
        require_review(cfg, STAGE, stage_output_dirs(cfg)["reports"])
        after = read_yaml_mapping(tuning)
        before = read_yaml_mapping(output / "before" / f"sub-{subject:02}.yaml")
        stripped = copy.deepcopy(after)
        # Ignore only the two intended decision fields. A concurrent/unexpected
        # change to any fit setting must remain visible to this preservation QC.
        for key in ("manual_exclude", "reviewed_model_sha256"):
            if key in before.get("ica", {}):
                stripped["ica"][key] = before["ica"][key]
            else:
                stripped["ica"].pop(key, None)
        if not stripped.get("ica") and "ica" not in before:
            stripped.pop("ica", None)
        stripped["stage4_review"].pop(STAGE, None)
        expected = copy.deepcopy(before)
        expected["stage4_review"].pop(STAGE, None)
        require(stripped == expected, "Unrelated tuning content changed during recording.")
        require(after["ica"]["manual_exclude"] == row["exclusions"] and
                after["ica"]["reviewed_model_sha256"] == row["model_sha256"], "Recorded exact-model choices disagree.")
        shutil.copy2(tuning, output / "after" / f"sub-{subject:02}.yaml")
        result.update(ok=True, tuning_sha256_after=file_hash(tuning), receipt=after["stage4_review"][STAGE])
        print(f"P{subject:02}: RECORDED provisional exact-model exclusions {row['exclusions']}", flush=True)
    except Exception as error:
        result.update(error=str(error), traceback=traceback.format_exc())
        print(f"P{subject:02}: RECORDING FAILED: {error}; inspect before snapshot", flush=True)
    (output / f"P{subject:02}_recording.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    """Validate frozen tables, join source workers, snapshot, record, then publish."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=PROJECT / "meg/params/pipeline_default.yaml")
    parser.add_argument("--subject-config-dir", type=Path, default=PROJECT / "meg/params")
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--validate-subject", type=int, help="Read-only recovery of one pending participant, never recording.")
    args = parser.parse_args()
    auth, eligible = load_authorization(args.authorization)
    args.authorization_sha256 = file_hash(args.authorization)
    args.choices_sha256 = file_hash(args.authorization.with_name("choices.json"))
    # Resolve the same shared derivative root through either mount convention.
    cfg = __import__("annotation_triage").subject_config(args, eligible[0]["subject"])
    root = Path(cfg["resolved_paths"]["derivatives_root"])
    args.authorization_relative_path = str(args.authorization.resolve().relative_to(root.resolve()))
    ready_path = root / "reports/ica_review_updates/final-ica-readiness-20261111-v1/summary.json"
    require(file_hash(ready_path) == auth["readiness_sha256"], "Frozen readiness changed.")
    ready = {row["subject"]: row for row in json.loads(ready_path.read_text())["participants"]}
    if args.validate_subject is not None:
        require(not args.record, "Single-subject recovery cannot write decisions.")
        eligible = [row for row in auth["participants"] if row["subject"] == args.validate_subject]
        require(len(eligible) == 1, "Recovery participant is outside the working cohort.")
    output = args.output.resolve()
    require(output.is_relative_to(root.resolve() / "reports/ica_review_updates"), "Use the derivative ICA evidence namespace.")
    output.mkdir(parents=False, exist_ok=False)
    workers = participant_workers(args.jobs, len(eligible))
    results = map_participants(validate_participant, [(args, row, ready[row["subject"]], output) for row in eligible], workers)
    summary = {"generated_utc": utc_now(), "authorization_sha256": args.authorization_sha256,
        "choices_sha256": args.choices_sha256, "resources": {"participant_processes": workers, "numerical_threads_per_worker": 1},
        "participants": results, "all_validated": all(row["ok"] for row in results),
        "user_visual_review_completed": False, "cleaned_output_acceptance": False, "final_scientific_acceptance": False,
        "ica_applied": False, "cluster_jobs_submitted": False}
    (output / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    require(summary["all_validated"], "Whole-table validation failed; NO bulk decision writes started.")
    if args.record:
        (output / "before").mkdir()
        (output / "after").mkdir()
        for row in results:
            path = Path(row["tuning"])
            require(file_hash(path) == row["tuning_sha256_before"], "Tuning changed before snapshots; no writes started.")
            shutil.copy2(path, output / "before" / f"sub-{row['subject']:02}.yaml")
        summary["recordings"] = map_participants(record_participant,
            [(args, row, checked, output) for row, checked in zip(eligible, results)], workers)
        summary["all_recorded"] = all(row["ok"] for row in summary["recordings"])
        (output / "recording.json").write_text(json.dumps(summary, indent=2) + "\n")
        require(summary["all_recorded"], "Partial recording; preserve snapshots and inspect participant results.")
    print("JOINED", output, flush=True)


if __name__ == "__main__":
    main()
