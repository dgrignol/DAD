#!/usr/bin/env python3
"""Validate and import per-RUN choices made from offline participant reports.

Usage from the workspace
------------------------
    python meg/code/preprocessing_basic/import_ocular_component_choices.py \
      --reports-dir /derivatives/annotation_triage/component_reports/LABEL \
      --decisions /path/to/completed_notes.csv
    # Default validates the COMPLETE table without writing choices.
    # Repeat with --apply to persist validated choices in the existing format.
    # Then run annotation_triage.py for those participants to generate masks/counts.

Independent participants use spawned processes (--jobs; default allocated CPUs).
Whole-table validation and recheck barriers precede every save worker. Runs stay
ordered within each participant; outputs and locks are disjoint between subjects.

CSV columns: participant, run, source_filename, blink, eye_movement, notes.
Participant accepts P14/14/sub-14; run accepts 2/R02/run02. The exact source column
may be blank when participant/run uniquely identifies the manifest acquisition.
Selections accept zero-based indices or ICA names separated by commas, or NONE.
Use DEFER for both groups to defer the run. Fully blank choices are unfinished;
one missing group is an error. Notes never imply a selection or reviewed zero.

Data flow: human CSV + exact report manifest/HTML -> whole-table syntax/identity
validation -> current source/cache validation -> explicit atomic component-choice
saves. Report/model bindings, current inputs, fitted ranges and manual precedence
are checked before ANY new choice is written. No existing manual work is changed.
Identical existing choices are kept verbatim; conflicting existing choices fail.
Use the ordinary explicit review workflow to reconcile a correction instead of
blindly replacing earlier work. Concurrent writes are checked again on SAVE.
This script creates no annotations or scientific acceptance. Only the existing
annotation_triage component-selection namespace is written with --apply.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

from annotation_triage import (SCHEMA, common_arguments, component_binding,
    has_manual_work, selection_path, subject_config, validate_selection)
from blink_annotation import load_cached_aid
from build_ocular_component_reports import (input_fingerprints, load_report_manifest)
from continuous_cleaning import require_review
from preprocessing_stage_utils import discover_raw_run_records, stage_output_dirs, utc_now
from review_ocular_components import parse_components, save_selection
from participant_parallel import add_parallel_argument, map_participants, participant_workers


# -----------------------------------------------------------------------------
# Explicit note syntax; never infer NONE from blanks, prose or another run
# -----------------------------------------------------------------------------

def numeric_label(value, prefix):
    """Normalize participant/run labels without guessing condition or source identity."""
    pattern = r"(?:P|SUB-?)?(\d+)" if prefix == "participant" else r"(?:R|RUN)?(\d+)"
    match = re.fullmatch(pattern, str(value).strip(), re.IGNORECASE)
    if not match:
        raise ValueError(f"Invalid {prefix}: {value!r}")
    return int(match.group(1))


def parse_table(path):
    """Parse all rows and report missing/duplicate/partial decisions before I/O-heavy checks."""
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"participant", "run", "blink", "eye_movement", "notes"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CSV requires participant, run, blink, eye_movement and notes columns.")
        rows, seen = [], set()
        for line, original in enumerate(reader, 2):
            if original.get(None):
                raise ValueError(f"Line {line}: too many CSV fields; quote comma-separated selections.")
            row = {key: (value or "").strip() for key, value in original.items() if key is not None}
            if not any(row.values()):
                continue
            subject = numeric_label(row["participant"], "participant")
            run = numeric_label(row["run"], "run")
            if not 4 <= subject <= 39 or run <= 0:
                raise ValueError(f"Line {line}: participant/run outside the current cohort.")
            if (subject, run) in seen:
                raise ValueError(f"Line {line}: duplicate participant/run.")
            seen.add((subject, run))
            if not row["blink"] and not row["eye_movement"]:
                continue
            if not row["blink"] or not row["eye_movement"]:
                raise ValueError(f"Line {line}: both groups need explicit indices, NONE or DEFER.")
            if (row["blink"].upper() == "DEFER") != (row["eye_movement"].upper() == "DEFER"):
                raise ValueError(f"Line {line}: use DEFER in both groups to defer the run.")
            # Check general syntax now; actual fitted bounds are checked against the model later.
            for key in ("blink", "eye_movement"):
                choice = parse_components(row[key], 1000000)
                if choice == "QUIT":
                    raise ValueError(f"Line {line}: QUIT is not a saved decision.")
            rows.append({**row, "subject": subject, "acquisition_run": run, "line": line})
    return rows


def equivalent(previous, proposed):
    """Permit idempotent reruns while retaining original save times/provenance verbatim."""
    return all(previous.get(key) == proposed.get(key) for key in
               ("status", "binding", "blink_components", "eye_movement_components", "observations"))


# -----------------------------------------------------------------------------
# Whole-table preflight; exact report/current-native validation precedes writes
# -----------------------------------------------------------------------------

def prepare_import_subject(task):
    """Validate one participant read-only; no save occurs inside this phase."""
    args, rows = task
    planned, configs, manifests, records_by_subject = [], {}, {}, {}
    for row in rows:
        subject = row["subject"]
        if subject not in configs:
            cfg = subject_config(args, subject)
            require_review(cfg, "03_A_MaxFilt", stage_output_dirs(cfg)["reports"])
            manifest_path = Path(args.reports_dir) / f"sub-{subject:02}_manifest.json"
            manifests[subject] = load_report_manifest(manifest_path)
            configs[subject] = cfg
            records_by_subject[subject] = discover_raw_run_records(cfg)
        cfg, manifest = configs[subject], manifests[subject]
        if manifest["subject_bids_id"] != cfg["layout"]["subject_bids_id"]:
            raise RuntimeError("Wrong participant in report manifest.")
        matches = [r for r in records_by_subject[subject] if r.acquisition_run_label == row["acquisition_run"]]
        if len(matches) != 1:
            raise ValueError(f"Line {row['line']}: no unique configured acquisition.")
        record = matches[0]
        if row.get("source_filename") and row["source_filename"] != record.source_filename:
            raise ValueError(f"Line {row['line']}: exact source filename disagrees.")
        report_runs = [r for r in manifest["runs"] if r["source_filename"] == record.source_filename]
        if len(report_runs) != 1 or report_runs[0]["status"] != "ready":
            raise RuntimeError("This run was not exported as pending component-review data.")
        report_run = report_runs[0]
        if has_manual_work(cfg, record):
            raise RuntimeError(f"Manual work takes precedence: {record.source_filename}")
        path = selection_path(cfg, record)
        previous = json.loads(path.read_text()) if path.exists() else None
        # The source-specific choice file is an expected output of a prior
        # import. Other inputs must still match the exact reviewed report.
        current = input_fingerprints(cfg, record)
        expected = report_run["inputs"]
        # Reports may be built on the cluster and imported on macOS. Textual
        # site/selector paths are not model identity. Native evidence must match;
        # load_cached_aid below independently validates CURRENT effective fit,
        # source, mask and sensor settings against the exact exported model.
        portable = lambda values: {k: v for k, v in values.items()
                                   if k != "choices" and not k.startswith("config_")}
        if portable(current) != portable(expected):
            raise RuntimeError(f"Inputs changed since report generation: {record.source_filename}")
        print(f"Validating notes against current cached model: {record.source_filename}", flush=True)
        cache = load_cached_aid(cfg, record)
        binding = component_binding(record, cache)
        if binding != report_run["binding"] or cache["model"].n_components_ != report_run["n_components"]:
            raise RuntimeError(f"Report model/source binding changed: {record.source_filename}")
        blink = parse_components(row["blink"], cache["model"].n_components_)
        eye = parse_components(row["eye_movement"], cache["model"].n_components_)
        payload = {"schema": SCHEMA, "purpose": "ocular_component_selection_only",
            "status": "deferred" if blink == "DEFER" else "selected",
            "binding": binding, "observations": row["notes"], "reviewer": args.reviewer,
            "saved_utc": utc_now(), "report_review": {"report_id": manifest["report_id"],
                "manifest_sha256": manifest["manifest_sha256"], "report_sha256": manifest["report_sha256"]}}
        if payload["status"] == "selected":
            payload.update(blink_components=blink, eye_movement_components=eye)
            validate_selection(payload, record, cache)
        if previous is not None and not equivalent(previous, payload):
            raise RuntimeError(f"Existing choices conflict; reconcile explicitly: {record.source_filename}")
        if input_fingerprints(cfg, record) != current:
            raise RuntimeError("Inputs changed during note validation; retry the batch.")
        planned.append({"subject": subject, "line": row["line"], "cfg": cfg, "record": record, "payload": payload,
                        "previous": previous, "input_snapshot": current})
    return planned


def prepare_import(args, rows):
    """Join parallel participant preflights before permitting ANY bulk save."""
    groups = {subject: [r for r in rows if r["subject"] == subject]
              for subject in sorted({r["subject"] for r in rows})}
    workers = participant_workers(getattr(args, "jobs", None), len(groups))
    print(f"Validating {len(groups)} participants with {workers} participant worker(s)", flush=True)
    batches = map_participants(prepare_import_subject, [(args, group) for group in groups.values()], workers)
    return sorted([item for batch in batches for item in batch], key=lambda item: item["line"])


def recheck_import_subject(task):
    """Read-only validation barrier after the potentially long model preflight."""
    args, items = task
    for item in items:
        latest = subject_config(args, item["subject"])
        if has_manual_work(latest, item["record"]) or input_fingerprints(latest, item["record"]) != item["input_snapshot"]:
            raise RuntimeError("Inputs changed during batch validation; no new choices were saved.")


def apply_import_subject(task):
    """Save one disjoint participant after the WHOLE table passes both barriers."""
    args, items = task
    for item in items:
        if item["previous"] is not None:
            continue
        latest = subject_config(args, item["subject"])
        if has_manual_work(latest, item["record"]) or input_fingerprints(latest, item["record"]) != item["input_snapshot"]:
            raise RuntimeError("Inputs changed immediately before SAVE; remaining rows were not written.")
        # Runs of one participant stay ordered because their convenience CSV is
        # shared. Other participants have disjoint JSON/CSV outputs and locks.
        path = save_selection(item["cfg"], item["record"], item["payload"], previous=None)
        print("Saved:", path, flush=True)


def run_import(args):
    """Join both whole-table validation phases, then save independent participants."""
    rows = parse_table(args.decisions)
    planned = prepare_import(args, rows)
    print(f"Validated {len(planned)} explicit run decision(s); {sum(p['previous'] is not None for p in planned)} already saved.", flush=True)
    if not args.apply:
        print("Validation only. Repeat with --apply to save these component choices.")
        return 0
    groups = [[item for item in planned if item["subject"] == subject]
              for subject in sorted({item["subject"] for item in planned})]
    workers = participant_workers(getattr(args, "jobs", None), len(groups))
    tasks = [(args, group) for group in groups]
    # No save worker is started until EVERY participant passes the fresh-input
    # barrier. Each per-source save also checks for a later concurrent edit.
    map_participants(recheck_import_subject, tasks, workers)
    map_participants(apply_import_subject, tasks, workers)
    print("Component choices saved. Run annotation_triage.py for these subjects to generate provisional masks and retention counts.")
    return 0


def main(argv=None):
    """Default to read-only validation; --apply is the explicit persistence action."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common_arguments(parser)
    add_parallel_argument(parser)
    parser.add_argument("--reports-dir", required=True)
    parser.add_argument("--decisions", required=True)
    parser.add_argument("--reviewer", default="User")
    parser.add_argument("--apply", action="store_true")
    return run_import(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
