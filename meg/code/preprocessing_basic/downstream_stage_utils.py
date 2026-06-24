#!/usr/bin/env python3
"""
downstream_stage_utils.py

Purpose
Provide shared, importable contracts for DAD MEG epoching, dRSA export, and
optional sanity workflows. The module validates the completed manual-ICA
lineage, discovers the subject-specific trajectory model, and converts canonical
trigger rows into run-isolated trial decisions with explicit exclusion
evidence.

Usage examples
1) Resolve the exact manual-ICA apply outputs recorded in a source manifest:

   from downstream_stage_utils import resolve_manual_ica_apply_runs
   runs, gate = resolve_manual_ica_apply_runs(cfg, source_manifest_path)

2) Build trial decisions from the canonical trigger CSV:

   rows = load_trigger_event_rows(events_csv)
   decisions = build_trial_decisions(
       rows,
       source_fif="sub04_run03_220526.fif",
       run_end_sample_exclusive=500000,
       event_id_map=cfg["epoching"]["event_id_map"],
       identity_range=cfg["subject_profile"]["identity_range"],
       exclusion_policy=cfg["epoching"]["trial_exclusion"],
   )

Data flow
resolved configuration + source_manifest.json -> exact raw acquisition records
-> manual-ICA apply records -> local lineage-matched cleaned FIFs.

canonical trigger rows -> per-run condition intervals -> trial evidence and
profile-specific numeric identity -> inclusion/exclusion decision + epoch
metadata.

Safety and scientific assumptions
- The source manifest must report a completed ``ica_manual`` apply pass with an
  explicit exclusion list and one completed record per raw acquisition run.
- Cluster paths stored in the manifest are provenance. Local files are resolved
  from the configured derivative root and exact source stem, then checked
  against the recorded output basename.
- Trial intervals never cross acquisition-run boundaries.
- Identity detection uses the configured numeric profile range rather than the
  stored event label. This is required for v20 identity 31, whose numeric value
  overlaps the legacy v18 always-visible condition code.
- A trial that terminates through an enabled gaze-break, catch, ESC, replay, or
  unresolved-transition exclusion may legitimately have no identity trigger.
  Included trials still require exactly one resolved identity.
- Resolved transition-warning rows are retained; unresolved transition rows are
  reportable evidence and can exclude their owning trial by policy.
- Trigger-150 rejection is a trial-level behavioral rule and is independent of
  physiological annotations handled by MNE ``reject_by_annotation``.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from preprocessing_stage_utils import (
    AcquisitionRunRecord,
    derivative_path_for_run,
    discover_raw_run_records,
    stage_output_dirs,
)


REQUIRED_TRIGGER_COLUMNS = {
    "trigger_code",
    "event_label",
    "sample_start_raw",
    "sample_end_raw",
    "raw_sample_frame",
    "source_fif",
    "acquisition_run",
    "is_transition_warning",
    "unresolved_transition",
    "behavioral_block",
    "behavioral_run",
    "behavioral_trial",
    "behavioral_context_source",
    "replay_status",
    "is_replay_trial",
    "is_gaze_break",
    "is_escape_termination",
}


@dataclass(frozen=True)
class ManualIcaApplyRun:
    """Bind one raw acquisition run to its lineage-matched cleaned FIF."""

    source: AcquisitionRunRecord
    cleaned_path: str
    manifest_output_path: str
    manifest_status: str


def csv_bool(value: Any) -> bool:
    """Parse bool-like CSV/JSON values without treating arbitrary text as true."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def optional_int(value: Any) -> int | None:
    """Return an integer for populated values and ``None`` for empty fields."""
    text = str(value).strip()
    if not text:
        return None
    return int(float(text))


def _require_one_file(root: Path, pattern: str, label: str) -> Path:
    """Discover exactly one file below ``root`` and fail on missing/ambiguity."""
    matches = sorted(path.resolve() for path in root.glob(pattern) if path.is_file())
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one {label} matching {pattern!r} under {root}; "
            f"found {len(matches)}: {[str(path) for path in matches]}"
        )
    return matches[0]


def discover_trajectory_models(cfg: Mapping[str, Any]) -> dict[str, str]:
    """
    Discover the exact subject-specific observed trajectory model.

    The observed MAT contains the shared path identities for all conditions.
    The expected continuation used for occluded-deviant dRSA is therefore the
    matching ``occluded_nondeviant`` trial in this same file; no separate
    predicted trajectory file is required.

    Returns
    -------
    dict
        Absolute ``observed`` MAT path, the expected-model source condition,
        and the documented stimulus version from the selected subject profile.
    """
    input_dir = Path(cfg["resolved_paths"]["subject_input_files_dir"]).resolve()
    discovery = cfg["data_discovery"]
    observed_pattern = Path(str(discovery["trajectory_mat_glob"])).name
    return {
        "observed": str(
            _require_one_file(input_dir, observed_pattern, "observed trajectory model")
        ),
        "expected_source_condition": "occluded_nondeviant",
        "stimuli_version": str(
            cfg["subject_profile"].get("versions", {}).get("stimuli", "")
        ),
    }


def load_trigger_event_rows(path: str | Path) -> list[dict[str, str]]:
    """Read and validate the canonical step-02 trigger event table."""
    events_path = Path(path).expanduser().resolve()
    if not events_path.is_file():
        raise FileNotFoundError(f"Canonical trigger events CSV not found: {events_path}")
    with events_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_TRIGGER_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                "Canonical trigger events CSV is missing required columns: "
                f"{sorted(missing)}"
            )
        return [dict(row) for row in reader]


def rows_by_source_fif(
    rows: Iterable[Mapping[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    """Group trigger rows by source FIF and order each run by raw sample."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        source_fif = str(row["source_fif"])
        grouped.setdefault(source_fif, []).append(dict(row))
    for source_fif in grouped:
        grouped[source_fif].sort(key=lambda item: int(item["sample_start_raw"]))
    return grouped


def resolve_manual_ica_apply_runs(
    cfg: Mapping[str, Any],
    source_manifest_path: str | Path,
) -> tuple[list[ManualIcaApplyRun], dict[str, Any]]:
    """
    Enforce the Stage-03 completion gate and resolve cleaned run paths.

    The local derivative path is derived from the raw source record and the
    configured derivative root. The manifest's recorded basename must match
    exactly, preventing sorted-position matching and filename replacement.
    """
    manifest_path = Path(source_manifest_path).expanduser().resolve()
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Source manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    stage = (
        manifest.get("derivative_lineage", {})
        .get("stages", {})
        .get("ica_manual")
    )
    if not isinstance(stage, Mapping):
        raise RuntimeError("Source manifest has no `ica_manual` derivative lineage.")
    if stage.get("processing_pass") != "apply" or stage.get("status") != "completed":
        raise RuntimeError(
            "Stage-03 completion gate failed: `ica_manual` must have "
            "processing_pass='apply' and status='completed'."
        )
    manual_exclude = stage.get("qc", {}).get("manual_exclude")
    if not isinstance(manual_exclude, list) or any(
        not isinstance(value, int) for value in manual_exclude
    ):
        raise RuntimeError(
            "Stage-03 completion gate failed: ICA exclusions must be explicit "
            "(an empty list is valid)."
        )

    raw_records = discover_raw_run_records(cfg)
    manifest_records = stage.get("records")
    if not isinstance(manifest_records, list):
        raise RuntimeError("ICA lineage records must be a list.")
    by_stem: dict[str, Mapping[str, Any]] = {}
    for record in manifest_records:
        if not isinstance(record, Mapping):
            raise RuntimeError("ICA lineage contains a non-mapping record.")
        source_stem = str(record.get("source_stem", ""))
        if not source_stem or source_stem in by_stem:
            raise RuntimeError(
                "ICA lineage source stems must be populated and unique: "
                f"{source_stem!r}"
            )
        by_stem[source_stem] = record

    sss_dir = stage_output_dirs(cfg)["sss"]
    resolved_runs: list[ManualIcaApplyRun] = []
    for source in raw_records:
        record = by_stem.get(source.source_stem)
        if record is None:
            raise RuntimeError(
                "Stage-03 completion gate failed: no ICA apply record for "
                f"{source.source_filename}."
            )
        output_path_value = record.get("output_path")
        if record.get("status") != "completed" or not output_path_value:
            raise RuntimeError(
                "Stage-03 completion gate failed: incomplete ICA record for "
                f"{source.source_filename}: status={record.get('status')!r}, "
                f"output_path={output_path_value!r}"
            )
        manifest_output = Path(str(output_path_value))
        expected_local = derivative_path_for_run(source, sss_dir, "_ica.fif")
        if manifest_output.name != expected_local.name:
            raise RuntimeError(
                "ICA lineage output basename does not preserve the exact source "
                f"stem: recorded={manifest_output.name}, expected={expected_local.name}"
            )
        if not expected_local.is_file():
            raise FileNotFoundError(
                "Lineage-matched manual-ICA output is missing at the configured "
                f"derivative root: {expected_local}"
            )
        resolved_runs.append(
            ManualIcaApplyRun(
                source=source,
                cleaned_path=str(expected_local),
                manifest_output_path=str(output_path_value),
                manifest_status=str(record["status"]),
            )
        )

    extra_stems = sorted(set(by_stem).difference(item.source_stem for item in raw_records))
    if extra_stems:
        raise RuntimeError(
            "ICA lineage contains records that do not match configured raw runs: "
            f"{extra_stems}"
        )
    gate = {
        "status": "passed",
        "processing_pass": "apply",
        "manual_exclude": list(manual_exclude),
        "n_raw_runs": len(raw_records),
        "n_cleaned_runs": len(resolved_runs),
        "source_manifest": str(manifest_path),
        "runs": [
            {
                **asdict(item.source),
                "cleaned_path": item.cleaned_path,
                "manifest_output_path": item.manifest_output_path,
            }
            for item in resolved_runs
        ],
    }
    return resolved_runs, gate


def _trial_value(row: Mapping[str, Any], key: str) -> Any:
    """Normalize selected trigger metadata fields for CSV/JSON/MNE metadata."""
    if key in {"behavioral_block", "behavioral_run", "behavioral_trial"}:
        return optional_int(row.get(key, ""))
    if key in {
        "is_transition_warning",
        "is_replay_trial",
        "is_gaze_break",
        "is_escape_termination",
    }:
        return csv_bool(row.get(key, False))
    return str(row.get(key, ""))


def build_trial_decisions(
    rows: Iterable[Mapping[str, Any]],
    *,
    source_fif: str,
    run_end_sample_exclusive: int,
    event_id_map: Mapping[str, int],
    identity_range: Iterable[int],
    exclusion_policy: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """
    Convert one run's canonical events into explicit epoch trial decisions.

    Each condition onset owns ``[onset, next condition onset)`` within the same
    acquisition run; the last interval ends at the run boundary supplied by
    the caller. Included trials require exactly one resolved identity trigger.
    An explicitly excluded partial trial may have no identity because gaze or
    ESC termination can occur before the identity frame is presented.
    """
    run_rows = sorted(
        (dict(row) for row in rows if str(row["source_fif"]) == source_fif),
        key=lambda item: int(item["sample_start_raw"]),
    )
    if not run_rows:
        raise RuntimeError(f"No canonical trigger rows found for {source_fif}.")

    condition_by_code = {int(code): name for name, code in event_id_map.items()}
    identity_values = list(identity_range)
    if not identity_values:
        raise ValueError("Identity range must contain at least one value.")
    identity_min = min(identity_values)
    identity_max = max(identity_values)

    # Unresolved condition rows cannot define trustworthy epochs. Resolved
    # transition-warning condition rows remain valid scientific events.
    condition_rows = [
        row
        for row in run_rows
        if int(row["trigger_code"]) in condition_by_code
        and not csv_bool(row["unresolved_transition"])
    ]
    if not condition_rows:
        raise RuntimeError(f"No resolved condition onsets found for {source_fif}.")

    configured_codes = {
        "gaze_break": set(int(value) for value in exclusion_policy["gaze_break_codes"]),
        "catch_type_1": set(
            int(value) for value in exclusion_policy["catch_type_1_codes"]
        ),
        "catch_type_2": set(
            int(value) for value in exclusion_policy["catch_type_2_codes"]
        ),
        "esc_termination": set(
            int(value) for value in exclusion_policy["esc_termination_codes"]
        ),
    }
    decisions: list[dict[str, Any]] = []
    for index, condition_row in enumerate(condition_rows):
        start = int(condition_row["sample_start_raw"])
        end = (
            int(condition_rows[index + 1]["sample_start_raw"])
            if index + 1 < len(condition_rows)
            else int(run_end_sample_exclusive)
        )
        if end <= start:
            raise RuntimeError(
                f"Invalid trial interval in {source_fif}: start={start}, end={end}"
            )
        interval_rows = [
            row
            for row in run_rows
            if start <= int(row["sample_start_raw"]) < end
        ]
        # Use the profile's numeric identity axis, not event_label. In v20 the
        # valid identity 31 overlaps the legacy v18 condition code and older
        # canonical CSVs can therefore carry the wrong descriptive label.
        identity_rows = [
            row
            for row in interval_rows
            if identity_min <= int(row["trigger_code"]) <= identity_max
            and not csv_bool(row["unresolved_transition"])
        ]
        evidence_codes = [int(row["trigger_code"]) for row in interval_rows]
        unresolved_count = sum(
            csv_bool(row["unresolved_transition"]) for row in interval_rows
        )
        replay = csv_bool(condition_row.get("is_replay_trial", False)) or (
            str(condition_row.get("replay_status", "")).lower() == "replay"
        )

        evidence = {
            name: sorted(code_set.intersection(evidence_codes))
            for name, code_set in configured_codes.items()
        }
        exclusion_reasons: list[str] = []
        for name in ("gaze_break", "catch_type_1", "catch_type_2", "esc_termination"):
            if evidence[name] and bool(exclusion_policy[f"exclude_{name}"]):
                exclusion_reasons.append(name)
        if replay and bool(exclusion_policy["exclude_replay_trials"]):
            exclusion_reasons.append("replay")
        if unresolved_count and bool(
            exclusion_policy["exclude_unresolved_transitions"]
        ):
            exclusion_reasons.append("unresolved_transition")

        # Multiple identities are ambiguous even when another trigger excludes
        # the trial. A missing identity is acceptable only for an already
        # excluded partial trial; candidates sent to MNE always have one.
        if len(identity_rows) > 1 or (
            not identity_rows and not exclusion_reasons
        ):
            raise RuntimeError(
                "Expected one resolved trial identity inside the run-isolated "
                f"interval for {source_fif} at sample {start}; found "
                f"{[(row['trigger_code'], row['sample_start_raw']) for row in identity_rows]}"
            )
        identity_row = identity_rows[0] if identity_rows else None
        decision = {
            "source_fif": source_fif,
            "source_fif_path": str(condition_row.get("source_fif_path", "")),
            "acquisition_run": int(condition_row["acquisition_run"]),
            "raw_sample_frame": start,
            "sample_start_raw": start,
            "trial_interval_end_raw_exclusive": end,
            "condition_name": condition_by_code[int(condition_row["trigger_code"])],
            "condition_code": int(condition_row["trigger_code"]),
            "trial_sequence_identity": (
                int(identity_row["trigger_code"]) if identity_row else None
            ),
            "identity_sample_start_raw": (
                int(identity_row["sample_start_raw"]) if identity_row else None
            ),
            "identity_mapping_status": (
                "resolved_numeric_profile_range"
                if identity_row
                else "missing_excluded_before_identity"
            ),
            "behavioral_block": _trial_value(condition_row, "behavioral_block"),
            "behavioral_run": _trial_value(condition_row, "behavioral_run"),
            "behavioral_trial": _trial_value(condition_row, "behavioral_trial"),
            "behavioral_context_source": _trial_value(
                condition_row,
                "behavioral_context_source",
            ),
            "replay_status": _trial_value(condition_row, "replay_status"),
            "is_replay_trial": replay,
            "is_transition_warning": _trial_value(
                condition_row,
                "is_transition_warning",
            ),
            "transition_partner_code": str(
                condition_row.get("transition_partner_code", "")
            ),
            "expected_trigger_code": str(
                condition_row.get("expected_trigger_code", "")
            ),
            "transition_decision_source": str(
                condition_row.get("transition_decision_source", "")
            ),
            "transition_resolution_confidence": str(
                condition_row.get("transition_resolution_confidence", "")
            ),
            "gaze_break_evidence": bool(evidence["gaze_break"]),
            "catch_type_1_evidence": bool(evidence["catch_type_1"]),
            "catch_type_2_evidence": bool(evidence["catch_type_2"]),
            "esc_termination_evidence": bool(evidence["esc_termination"]),
            "unresolved_transition_count": int(unresolved_count),
            "evidence_trigger_codes": ",".join(str(code) for code in evidence_codes),
            "exclusion_reasons": ",".join(exclusion_reasons),
            "included_by_trial_policy": not exclusion_reasons,
            "annotation_rejection_policy": "handled_by_mne_reject_by_annotation",
        }
        decisions.append(decision)
    return decisions


def identity_coverage(
    observed: Iterable[int],
    required: Iterable[int],
    *,
    missing_policy: str,
) -> tuple[dict[int, int], list[int]]:
    """
    Count condition epochs by identity and enforce the missing-data policy.

    ``missing_policy='error'`` raises when any required identity is absent;
    ``'nan'`` returns the missing list so callers can retain an explicit NaN
    tensor row. Observed identities outside the required axis always fail.
    """
    observed_values = [int(value) for value in observed]
    required_values = [int(value) for value in required]
    unexpected = sorted(set(observed_values).difference(required_values))
    if unexpected:
        raise RuntimeError(
            f"Observed identities are outside the configured axis: {unexpected}"
        )
    counts = {
        identity: observed_values.count(identity) for identity in required_values
    }
    missing = [identity for identity, count in counts.items() if count == 0]
    if missing_policy not in {"error", "nan"}:
        raise ValueError(
            f"Unsupported missing identity policy: {missing_policy!r}"
        )
    if missing and missing_policy == "error":
        raise RuntimeError(
            "Required identities have no surviving epochs: "
            f"{missing}"
        )
    return counts, missing


__all__ = [
    "ManualIcaApplyRun",
    "build_trial_decisions",
    "csv_bool",
    "discover_trajectory_models",
    "identity_coverage",
    "load_trigger_event_rows",
    "optional_int",
    "resolve_manual_ica_apply_runs",
    "rows_by_source_fif",
]
