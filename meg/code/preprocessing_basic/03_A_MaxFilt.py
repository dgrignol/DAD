#!/usr/bin/env python3
"""
03_A_MaxFilt.py

Purpose
Run configuration-driven bad-channel review and Maxwell filtering for every
task MEG acquisition run belonging to one DAD participant. Every run is mapped
from its static ``dev_head_t`` to one central actual run's head position.

Usage examples
1) Noninteractive review with FLUX-aligned defaults:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_A_MaxFilt.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --processing-pass review --print-summary

2) Interactive policy review. QC plots are displayed when a GUI backend is
   available; use ``--no-display`` to inspect the saved PNGs before answering
   the terminal prompts:

   python .../03_A_MaxFilt.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --processing-pass review \
     --interactive --decisions-out .../03_tuning_choices.yaml

3) Apply accepted tuning and intentionally replace review derivatives:

   python .../03_A_MaxFilt.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --tuning-config .../03_tuning_choices.yaml \
     --processing-pass apply --overwrite

4) Resume after an interrupted network-volume write. Existing outputs are
   read and validated against source lineage and the selected bad-channel list
   before they are skipped:

   python .../03_A_MaxFilt.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --tuning-config .../03_tuning_choices.yaml \
     --processing-pass review --resume

Inputs
- Layered configuration; ``--site-config`` overrides machine paths without
  changing scientific settings.
- Raw task FIFs discovered from ``data_discovery.meg_fif_glob`` below the
  unique external ``subNN_*`` folder.
- Facility cross-talk and calibration resources configured in
  ``maxwell_filter``. They are not inferred from a subject raw-data folder.

Outputs
- ``<exact-source-stem>_raw_sss.fif`` for every acquisition run.
- Pass-specific JSON and self-contained MNE HTML QC reports, with embedded
  PNG figures, explanations and scientific review criteria. Per-run noisy/flat
  score heatmaps replace the old detection PSD; raw traces highlight detected
  channels in a reproducible random 20-second window. Use maxwell_qc.py to
  refresh these displays from existing scores without rerunning detection/SSS.
- Optional accepted tuning YAML for reuse by ``run_pipeline.py``.
- Derivative lineage in the shared source manifest when supplied.

Data flow
raw acquisition runs -> independent per-run Maxwell bad-channel detection ->
central static destination selection -> configured/interactive bad-channel
policy -> per-run Maxwell filtering -> before/after QC ->
lineage-preserving SSS FIFs and reports.

SQUID extension (workflow Stage 4A, legacy script ID 03_A)
Use --squid-prepare-only to prepare candidates from accepted channel evidence
without changing SSS, reports, tuning or existing 4B manual reviews. Review with
review_squid_jumps.py --mode review --source EXACT_SOURCE before the next SSS
pass. Confirmed BAD_squid_jumps windows are explicitly skipped by SSS; channel
actions default to participant union. Missing reviews stop enabled new 4A
processing before SSS. Existing accepted 4A/4B work is never migrated merely
by installing this extension. --overwrite archives changed SSS/evidence first;
compatible SSS files are reused, retaining their manual-eye review bindings.

Scientific assumptions
- Final temporal settings belong to 4B/application, independently of SSS.
  Existing 2026-10-23 report/receipt identities are validated read-only through
  the explicit legacy projection; this policy update requires no 4A rerun.
- Acquisition-run identity comes only from each FIF's ``runNN`` token.
- The default ``union`` policy applies every channel detected in any run to all
  runs. ``reference_run`` and ``per_run`` remain explicit alternatives.
- The common destination is an actual acquisition run, selected as the medoid
  of pairwise RMS displacement over the subject's digitization point cloud.
- Existing raw-header bad channels and configured manual additions/removals are
  retained in the decision record.

Safety
The complete raw tree is read-only. Every output path passes the shared
``assert_output_path_outside_raw`` guard before directories or files are
created. Existing outputs are rejected unless overwrite is explicitly enabled.
Resume mode never trusts an existing FIF solely because its filename matches:
the stored source/stage description and applied bad-channel list must agree
with the current decision.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Mapping

import mne
import numpy as np
from mne.io.constants import FIFF

from continuous_cleaning import (stage_configuration, bind_report, check_preservation, digest, file_hash, file_stamp, filter_design)

from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import (
    AcquisitionRunRecord,
    derivative_path_for_run,
    discover_raw_run_records,
    ensure_output_directory,
    guarded_path,
    lineage_record,
    prepare_output_path,
    stage_output_dirs,
    update_source_manifest_stage,
    update_tuning_file,
    utc_now,
    write_json,
)

from render_stage03_mne_reports import render_stage_report
from maxwell_qc import build_run_qc, load_scores, save_heatmaps
from squid_jumps import (LABEL as SQUID_LABEL, accepted_report, application as squid_application,
    attach_intervals, channel_basis_configuration, load_review as load_squid_review,
    parameters as squid_parameters, prepare as prepare_squid, root as squid_root,
    run_action as squid_run_action, maxwell_filter_with_jumps, remember_channel_basis,
    read_channel_basis)


STAGE_NAME = "maxwell_filter"
OUTPUT_SUFFIX = "_raw_sss.fif"


# -----------------------------------------------------------------------------
# Digitization compatibility checks
# -----------------------------------------------------------------------------


def digitization_summary(info: mne.Info) -> dict[str, int]:
    """
    Count digitization points relevant to automatic Maxwell-origin fitting.

    DAD task FIFs contain Polhemus head-shape points as ``extra`` points.
    MNE 1.8 and 1.9 omit those points while reading these files, even though
    cardinal and HPI points remain visible. Reporting each kind makes that
    reader-version failure distinguishable from genuinely absent digitization.
    """
    counts = {
        "cardinal": 0,
        "hpi": 0,
        "eeg": 0,
        "extra": 0,
        "total": 0,
    }
    kind_to_name = {
        int(FIFF.FIFFV_POINT_CARDINAL): "cardinal",
        int(FIFF.FIFFV_POINT_HPI): "hpi",
        int(FIFF.FIFFV_POINT_EEG): "eeg",
        int(FIFF.FIFFV_POINT_EXTRA): "extra",
    }
    for point in info.get("dig") or []:
        counts["total"] += 1
        name = kind_to_name.get(int(point["kind"]))
        if name is not None:
            counts[name] += 1
    return counts


def validate_auto_origin_digitization(
    info: mne.Info,
    *,
    origin: str | list[float] | tuple[float, float, float],
    source_filename: str,
) -> dict[str, int]:
    """
    Validate that MNE exposed enough head-shape points for ``origin="auto"``.

    A numerical origin does not require a sphere fit and therefore bypasses
    this validation. No numerical fallback is selected automatically: changing
    the Maxwell origin is a scientific decision, not an error-recovery policy.
    """
    summary = digitization_summary(info)
    if origin == "auto" and summary["extra"] + summary["eeg"] < 4:
        raise RuntimeError(
            f"{source_filename}: origin='auto' requires at least four extra/EEG "
            "digitization points, but this MNE reader exposed "
            f"{summary}. DAD task FIFs require MNE >= 1.10.0 to read their "
            f"embedded Polhemus head-shape points; found MNE {mne.__version__}. "
            "Upgrade the environment rather than silently changing the origin."
        )
    return summary


# -----------------------------------------------------------------------------
# Common static destination selection
# -----------------------------------------------------------------------------


def _centrality_head_points(info: mne.Info) -> np.ndarray:
    """
    Return actual digitization points in head coordinates for transform scoring.

    One common point cloud is transformed through every run's static
    ``dev_head_t``. This evaluates translation and rotation on the participant's
    measured head geometry rather than comparing transform coefficients.
    """
    included_kinds = {
        int(FIFF.FIFFV_POINT_CARDINAL),
        int(FIFF.FIFFV_POINT_HPI),
        int(FIFF.FIFFV_POINT_EEG),
        int(FIFF.FIFFV_POINT_EXTRA),
    }
    points = [
        np.asarray(point["r"], dtype=float)
        for point in info.get("dig") or []
        if int(point["coord_frame"]) == int(FIFF.FIFFV_COORD_HEAD)
        and int(point["kind"]) in included_kinds
    ]
    if len(points) < 4:
        raise RuntimeError(
            "Central-run destination selection requires at least four "
            f"digitization points in head coordinates; found {len(points)}."
        )
    return np.vstack(points)


def select_central_destination_run(
    dev_head_t_by_run: Mapping[int, Any],
    *,
    head_points: np.ndarray,
) -> tuple[int, dict[str, Any]]:
    """
    Select the actual run minimizing summed pairwise head-point displacement.

    ``dev_head_t`` maps device coordinates to head coordinates. Its inverse
    maps one shared measured head point cloud into each run's device frame.
    Pairwise RMS displacement therefore captures both translation and rotation.
    """
    if not dev_head_t_by_run:
        raise ValueError("No static dev_head_t transforms were supplied.")
    points = np.asarray(head_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError(
            "Head points must have shape (n_points, 3); "
            f"received {points.shape}."
        )

    device_points: dict[int, np.ndarray] = {}
    for run, transform in dev_head_t_by_run.items():
        if transform is None:
            raise RuntimeError(f"run{int(run):02d} has no static dev_head_t.")
        matrix = np.asarray(transform["trans"], dtype=float)
        if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
            raise RuntimeError(f"run{int(run):02d} has an invalid static transform.")
        if (int(transform.get("from", 1)) != 1 or int(transform.get("to", 4)) != 4
                or not np.allclose(matrix[3], [0, 0, 0, 1])
                or not np.allclose(matrix[:3, :3].T @ matrix[:3, :3], np.eye(3), atol=1e-5)
                or not np.isclose(np.linalg.det(matrix[:3, :3]), 1., atol=1e-5)):
            raise RuntimeError(f"run{int(run):02d} has invalid device-to-head frame/rotation.")
        head_to_device = np.linalg.inv(matrix)
        homogeneous = np.column_stack([points, np.ones(len(points))])
        device_points[int(run)] = (head_to_device @ homogeneous.T).T[:, :3]

    pairwise_mm: dict[str, dict[str, float]] = {}
    score_mm: dict[int, float] = {}
    for run_a in sorted(device_points):
        distances = {}
        for run_b in sorted(device_points):
            delta = device_points[run_a] - device_points[run_b]
            rms_m = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
            distances[f"run{run_b:02d}"] = rms_m * 1000.0
        pairwise_mm[f"run{run_a:02d}"] = distances
        score_mm[run_a] = float(sum(distances.values()))

    minimum = min(score_mm.values())
    tied = [run for run, score in score_mm.items() if abs(score - minimum) <= 1e-6]
    central_run = min(tied)
    return central_run, {
        "metric": "summed_pairwise_rms_digitization_displacement_mm",
        "digitization_point_count": int(len(points)),
        "tie_candidates": sorted(tied), "tie_tolerance_mm": 1e-6,
        "rationale": "Measured actual-run medoid includes translation and rotation of head geometry; no cHPI trajectory is available.",
        "references": ["PIPELINE_SPEC_basic.md#03_a-maxwell-filtering", "https://mne.tools/1.10/generated/mne.preprocessing.maxwell_filter.html"],
        "score_by_run_mm": {
            f"run{run:02d}": score_mm[run] for run in sorted(score_mm)
        },
        "pairwise_rms_mm": pairwise_mm,
    }


def resolve_destination(
    *,
    destination_cfg: Mapping[str, Any],
    dev_head_t_by_run: Mapping[int, Any],
    head_points: np.ndarray,
) -> tuple[int, Any, dict[str, Any]]:
    """Resolve the central actual-run destination from static transforms."""
    central_run, centrality = select_central_destination_run(
        dev_head_t_by_run,
        head_points=head_points,
    )
    policy = str(destination_cfg["policy"])
    configured_run = destination_cfg["acquisition_run"]
    if policy != "central_run" or configured_run is not None:
        raise ValueError(
            "Maxwell destination must use policy='central_run' with "
            "acquisition_run=null."
        )
    return central_run, dev_head_t_by_run[central_run], {
        "policy": policy,
        "configured_acquisition_run": configured_run,
        "selected_acquisition_run": central_run,
        "central_acquisition_run": central_run,
        **centrality,
    }


def destination_transform_fingerprint(transform: Any) -> str:
    """Return a stable fingerprint for resume and lineage validation."""
    matrix = np.asarray(transform["trans"], dtype=float)
    canonical = json.dumps(
        np.round(matrix, decimals=12).tolist(),
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(canonical).hexdigest()


# -----------------------------------------------------------------------------
# Pure bad-channel policy helpers
# -----------------------------------------------------------------------------


def _deduplicate(values: list[str]) -> list[str]:
    """Return channel names once while preserving first occurrence order."""
    return list(dict.fromkeys(str(value) for value in values))


def _run_override(
    manual_cfg: Mapping[str, Any],
    acquisition_run: int,
) -> tuple[list[str], list[str]]:
    """
    Return configured per-run manual bad-channel additions and removals.

    Run keys may be written as ``1``, ``"1"``, ``"01"``, or ``"run01"``.
    """
    by_run = manual_cfg.get("by_run", {})
    if not isinstance(by_run, Mapping):
        return [], []
    candidate_keys = (
        acquisition_run,
        str(acquisition_run),
        f"{acquisition_run:02d}",
        f"run{acquisition_run:02d}",
    )
    value = next((by_run[key] for key in candidate_keys if key in by_run), None)
    if value is None:
        return [], []
    if isinstance(value, list):
        return _deduplicate(value), []
    if not isinstance(value, Mapping):
        raise ValueError(
            "Per-run manual bad-channel settings must be a list or mapping; "
            f"run={acquisition_run}, value={value!r}"
        )
    return (
        _deduplicate(list(value.get("add", []))),
        _deduplicate(list(value.get("remove", []))),
    )


def select_bad_channels_by_policy(
    detected_by_run: Mapping[int, list[str]],
    *,
    policy: str,
    reference_acquisition_run: int | None,
    manual_cfg: Mapping[str, Any],
) -> tuple[dict[int, list[str]], int]:
    """
    Select applied bad channels for each run under one explicit policy.

    Parameters
    ----------
    detected_by_run
        Per-run automatic detections, including pre-existing raw-header bads.
    policy
        ``reference_run``, ``union``, or ``per_run``.
    reference_acquisition_run
        Reference run for ``reference_run``. ``None`` selects the smallest
        acquisition-run label.
    manual_cfg
        Mapping with global ``add``/``remove`` lists and optional ``by_run``.

    Returns
    -------
    selected_by_run, resolved_reference_run
        Applied channel lists and the resolved reference-run label.
    """
    if not detected_by_run:
        raise ValueError("No per-run bad-channel detections were supplied.")
    run_labels = sorted(int(run) for run in detected_by_run)
    resolved_reference = (
        min(run_labels)
        if reference_acquisition_run is None
        else int(reference_acquisition_run)
    )
    if resolved_reference not in detected_by_run:
        raise ValueError(
            "Configured reference acquisition run is not present: "
            f"{resolved_reference}; available={run_labels}"
        )

    if policy == "reference_run":
        base_by_run = {
            run: list(detected_by_run[resolved_reference]) for run in run_labels
        }
    elif policy == "union":
        union = sorted(
            {
                channel
                for channels in detected_by_run.values()
                for channel in channels
            }
        )
        base_by_run = {run: list(union) for run in run_labels}
    elif policy == "per_run":
        base_by_run = {
            run: list(detected_by_run[run]) for run in run_labels
        }
    else:
        raise ValueError(f"Unsupported bad-channel policy: {policy!r}")

    global_add = _deduplicate(list(manual_cfg.get("add", [])))
    global_remove = set(_deduplicate(list(manual_cfg.get("remove", []))))
    selected_by_run: dict[int, list[str]] = {}
    for run in run_labels:
        run_add, run_remove = _run_override(manual_cfg, run)
        removed = global_remove | set(run_remove)
        selected_by_run[run] = [
            channel
            for channel in _deduplicate(
                [*base_by_run[run], *global_add, *run_add]
            )
            if channel not in removed
        ]
    return selected_by_run, resolved_reference


# -----------------------------------------------------------------------------
# QC plotting and interactive review
# -----------------------------------------------------------------------------


def _parse_channel_prompt(text: str) -> list[str]:
    """Parse comma/space separated channel names from an interactive prompt."""
    return _deduplicate(
        [token.strip() for token in text.replace(",", " ").split() if token.strip()]
    )


def interactive_bad_channel_review(
    *,
    detected_by_run: Mapping[int, list[str]],
    default_policy: str,
    default_reference_run: int | None,
    default_manual_cfg: Mapping[str, Any],
) -> tuple[str, int | None, dict[str, Any]]:
    """
    Prompt until the operator accepts a policy and manual channel adjustments.

    Saved QC figures are generated before this function is called. The prompt
    prints every per-run list so ``--no-display`` remains fully usable.
    """
    print("\nPer-run automatic bad-channel detections:")
    for run, channels in sorted(detected_by_run.items()):
        print(f"  run{run:02d}: {channels}")

    policy = default_policy
    reference_run = default_reference_run
    manual_cfg = {
        "add": list(default_manual_cfg.get("add", [])),
        "remove": list(default_manual_cfg.get("remove", [])),
        "by_run": dict(default_manual_cfg.get("by_run", {})),
    }
    while True:
        entered_policy = input(
            "Bad-channel policy "
            f"[reference_run/union/per_run] (default {policy}): "
        ).strip()
        if entered_policy:
            if entered_policy not in ("reference_run", "union", "per_run"):
                print("Invalid policy; choose reference_run, union, or per_run.")
                continue
            policy = entered_policy
        if policy == "reference_run":
            default_ref = (
                min(detected_by_run)
                if reference_run is None
                else int(reference_run)
            )
            entered_ref = input(
                f"Reference acquisition run (default {default_ref}): "
            ).strip()
            reference_run = int(entered_ref) if entered_ref else default_ref

        entered_add = input(
            "Global manual bad-channel additions "
            f"(default {manual_cfg['add']}, blank keeps): "
        ).strip()
        if entered_add:
            manual_cfg["add"] = _parse_channel_prompt(entered_add)
        entered_remove = input(
            "Channels to remove from the selected bad list "
            f"(default {manual_cfg['remove']}, blank keeps): "
        ).strip()
        if entered_remove:
            manual_cfg["remove"] = _parse_channel_prompt(entered_remove)

        selected, resolved_reference = select_bad_channels_by_policy(
            detected_by_run,
            policy=policy,
            reference_acquisition_run=reference_run,
            manual_cfg=manual_cfg,
        )
        print("\nProposed applied bad channels:")
        for run, channels in sorted(selected.items()):
            print(f"  run{run:02d}: {channels}")
        answer = input("Accept this bad-channel decision? [Y/n]: ").strip().lower()
        if answer in ("", "y", "yes"):
            return policy, resolved_reference, manual_cfg


# -----------------------------------------------------------------------------
# Stage execution
# -----------------------------------------------------------------------------


def _detect_bad_channels(
    record: AcquisitionRunRecord,
    *,
    cfg: Mapping[str, Any],
    calibration_file: Path,
    crosstalk_file: Path,
) -> tuple[list[str], dict[str, Any], mne.io.BaseRaw]:
    """Detect noisy/flat channels for one run and retain header provenance."""
    raw = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
    input_bads = _deduplicate(list(raw.info.get("bads", [])))
    detection_raw = raw.copy()
    detection_raw.info["bads"] = []
    find_cfg = cfg["maxwell_filter"]["find_bad_channels"]
    dig_summary = validate_auto_origin_digitization(
        detection_raw.info,
        origin=find_cfg["origin"],
        source_filename=record.source_filename,
    )
    noisy, flat, scores = mne.preprocessing.find_bad_channels_maxwell(
        detection_raw,
        limit=float(find_cfg["limit"]),
        duration=float(find_cfg["duration_sec"]),
        min_count=int(find_cfg["min_count"]),
        origin=find_cfg["origin"],
        h_freq=float(find_cfg["h_freq_hz"]),
        calibration=str(calibration_file),
        cross_talk=str(crosstalk_file),
        return_scores=True,
        verbose="WARNING",
    )
    detected = _deduplicate([*input_bads, *noisy, *flat])
    return (
        detected,
        {
            "input_header_bads": input_bads,
            "auto_noisy_channels": list(noisy),
            "auto_flat_channels": list(flat),
            "detected_bad_channels": detected,
            "digitization_points": dig_summary,
            "scores": scores,
        },
        raw,
    )


def maxwell_signature(record, cfg, calibration_file, crosstalk_file, destination_fingerprint, applied_bads):
    """Bind SSS to source stat, resource contents and all operative decisions."""
    settings = dict(cfg["maxwell_filter"])
    settings.pop("calibration_file", None)
    settings.pop("crosstalk_file", None)
    payload = {"source": file_stamp(record.source_path), "settings": settings,
                   "calibration": file_hash(calibration_file), "crosstalk": file_hash(crosstalk_file),
                   "destination": destination_fingerprint, "bads": applied_bads}
    action = squid_run_action(cfg, record.source_filename)
    if action["intervals"] or action["bad_channels"]:
        payload["squid_jumps"] = action
    return digest(payload)


def _maxwell_filter_one_run(
    *,
    record: AcquisitionRunRecord,
    applied_bads: list[str],
    output_path: Path,
    cfg: Mapping[str, Any],
    calibration_file: Path,
    crosstalk_file: Path,
    destination_run: int,
    destination_transform: Any,
    destination_fingerprint: str,
    overwrite: bool,
) -> None:
    """Apply Maxwell filtering to one full run at the common destination."""
    raw = mne.io.read_raw_fif(
        record.source_path,
        preload=True,
        verbose="WARNING",
    )
    missing = sorted(set(applied_bads) - set(raw.ch_names))
    if missing:
        raise ValueError(
            f"Configured bad channels are absent from {record.source_filename}: "
            f"{missing}"
        )
    raw.info["bads"] = list(applied_bads)
    attach_intervals(raw, squid_run_action(cfg, record.source_filename)["intervals"])
    raw.fix_mag_coil_types()
    filter_cfg = cfg["maxwell_filter"]["filter"]
    validate_auto_origin_digitization(
        raw.info,
        origin=filter_cfg["origin"],
        source_filename=record.source_filename,
    )
    filtered = maxwell_filter_with_jumps(
        raw,
        origin=filter_cfg["origin"],
        calibration=str(calibration_file),
        cross_talk=str(crosstalk_file),
        st_duration=filter_cfg.get("st_duration_sec"),
        st_correlation=float(filter_cfg["st_correlation"]),
        destination=destination_transform,
        skip_by_annotation=("edge", "bad_acq_skip", SQUID_LABEL),
        verbose="WARNING",
    )
    # Maxwell filtering may rewrite Info; reattach the explicit applied policy.
    filtered.info["bads"] = list(applied_bads)
    filtered.info["description"] = (
        f"DAD source={record.source_filename}; "
        f"acquisition_run={record.acquisition_run_label}; "
        f"stage=03_A_MaxFilt; destination_run={destination_run}; "
        f"destination_transform_sha256={destination_fingerprint}"
    )
    signature = maxwell_signature(record, cfg, calibration_file, crosstalk_file,
                                  destination_fingerprint, applied_bads)
    filtered.info["description"] += f"; stage4_input_sha256={signature}"
    check_preservation(raw, filtered)
    filtered.save(output_path, overwrite=overwrite, verbose="WARNING")


def _validate_resumable_output(
    *,
    path: Path,
    record: AcquisitionRunRecord,
    applied_bads: list[str],
    destination_run: int,
    destination_transform: Any,
    destination_fingerprint: str,
    raw_root: str | Path,
    expected_signature: str | None = None,
    squid_intervals: list | None = None,
) -> dict[str, Any]:
    """
    Validate an existing Maxwell derivative before resume skips computation.

    A filename match is insufficient because a previous run may have used a
    different bad-channel policy. The FIF header must identify the same raw
    acquisition run and carry exactly the currently selected bad channels.
    """
    resolved = guarded_path(path, raw_root)
    try:
        existing = mne.io.read_raw_fif(
            resolved,
            preload=False,
            verbose="ERROR",
        )
    except Exception as exc:  # pylint: disable=broad-except
        raise RuntimeError(
            f"Cannot resume from unreadable Maxwell derivative {resolved}: {exc}"
        ) from exc
    if expected_signature is not None:
        if f"stage4_input_sha256={expected_signature}" not in str(existing.info.get("description")):
            raise RuntimeError("Cannot resume SSS: upstream source/settings/resources changed.")
        source = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
        source.info["bads"] = list(applied_bads)
        if squid_intervals:
            # FIF stores annotation endpoints at float32 precision. Compare new
            # SQUID endpoints on the native sample grid, while every original
            # non-SQUID annotation still passes the existing exact contract.
            expected = source.copy()
            attach_intervals(expected, squid_intervals)
            def squid_samples(raw):
                return sorted((round(a["onset"]*raw.info["sfreq"]),
                    round((a["onset"]+a["duration"])*raw.info["sfreq"]))
                    for a in raw.annotations if a["description"] == SQUID_LABEL)
            if squid_samples(expected) != squid_samples(existing):
                raise RuntimeError("Cannot resume SSS: reviewed SQUID intervals changed native samples.")
            compared = existing.copy()
            for raw in (source, compared):
                raw.annotations.delete(np.flatnonzero(raw.annotations.description == SQUID_LABEL))
            check_preservation(source, compared)
        else:
            check_preservation(source, existing)
    description = str(existing.info.get("description") or "")
    required_tokens = (
        f"source={record.source_filename}",
        f"acquisition_run={record.acquisition_run_label}",
        "stage=03_A_MaxFilt",
        f"destination_run={destination_run}",
        f"destination_transform_sha256={destination_fingerprint}",
    )
    missing_tokens = [token for token in required_tokens if token not in description]
    if missing_tokens:
        raise RuntimeError(
            "Cannot resume from Maxwell derivative with mismatched lineage: "
            f"{resolved}; missing description tokens={missing_tokens}"
        )
    existing_bads = list(existing.info.get("bads", []))
    if existing_bads != list(applied_bads):
        raise RuntimeError(
            "Cannot resume from Maxwell derivative created with different bad "
            f"channels: {resolved}; existing={existing_bads}, "
            f"selected={applied_bads}"
        )
    existing_transform = existing.info.get("dev_head_t")
    if existing_transform is None or not np.allclose(
        np.asarray(existing_transform["trans"], dtype=float),
        np.asarray(destination_transform["trans"], dtype=float),
        rtol=0.0,
        atol=1e-9,
    ):
        raise RuntimeError(
            "Cannot resume from Maxwell derivative with a different common "
            f"destination transform: {resolved}"
        )
    if existing.n_times <= 0:
        raise RuntimeError(
            f"Cannot resume from empty Maxwell derivative: {resolved}"
        )
    return {
        "n_times": int(existing.n_times),
        "sfreq_hz": float(existing.info["sfreq"]),
        "description": description,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the noninteractive-by-default stage CLI."""
    parser = argparse.ArgumentParser(
        description="Run DAD per-acquisition-run Maxwell filtering."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--source-manifest", default=None)
    parser.add_argument(
        "--processing-pass",
        choices=("review", "apply"),
        default="review",
    )
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument(
        "--no-display",
        action="store_true",
        help="Save decision plots but do not open GUI windows.",
    )
    parser.add_argument(
        "--decisions-out",
        default=None,
        help="Accepted interactive choices YAML; defaults in the report directory.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--squid-prepare-only", action="store_true",
        help="Prepare native SQUID candidates from accepted channels; never run SSS or change 4B files.")
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Validate and skip existing lineage-matched derivatives; fail on "
            "any mismatch. Mutually exclusive with --overwrite."
        ),
    )
    parser.add_argument("--print-summary", action="store_true")
    return parser


def archive_squid_rebuild(cfg, *, output_path=None):
    """Preserve old 4A evidence and changed SSS before explicit replacement.

    Small report/score artifacts are copied once under their content hash.
    A changed large SSS file is moved within the derivative filesystem. Manual
    annotations/tuning are never outputs here and retain their original bytes.
    """
    report_dir = Path(cfg["resolved_paths"]["derivatives_reports_dir"])
    report = report_dir / "03_A_MaxFilt_review_report.json"
    if not report.exists():
        return
    archive = ensure_output_directory(report_dir / "history" / "squid_jumps" / file_hash(report),
                                      cfg["resolved_paths"]["data_root"])
    for path in report_dir.glob("03_A_MaxFilt_*"):
        target = guarded_path(archive / path.name, cfg["resolved_paths"]["data_root"])
        if path.is_file() and not target.exists():
            shutil.copy2(path, target)
    if output_path is not None and Path(output_path).exists():
        target = guarded_path(archive / Path(output_path).name, cfg["resolved_paths"]["data_root"])
        if target.exists():
            raise RuntimeError(f"Archived SSS already exists; reconcile before replacement: {target}")
        Path(output_path).replace(target)


def main(argv: list[str] | None = None) -> int:
    """Run bad-channel review, Maxwell filtering, QC, and lineage recording."""
    args = build_arg_parser().parse_args(argv)
    try:
        cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            args.tuning_config,
            site_config_path=args.site_config,
            required_sections=("maxwell_filter",),
        )
        raw_root = cfg["resolved_paths"]["data_root"]
        squid_cfg = squid_parameters(cfg)
        if args.squid_prepare_only:
            # This branch deliberately precedes every scientific-output write.
            from review_squid_jumps import channel_lists
            records = discover_raw_run_records(cfg)
            channels = channel_lists(cfg, records)
            for record in records:
                path, proposal = prepare_squid(cfg, record, channels[record.source_filename])
                print(f"{record.source_filename}: {len(proposal['candidates'])} SQUID candidates: {path}")
            return 0
        overwrite = bool(args.overwrite or cfg["runtime"]["overwrite_existing"])
        if args.resume and overwrite:
            raise ValueError("--resume and --overwrite are mutually exclusive.")
        maxwell_cfg = cfg["maxwell_filter"]
        if maxwell_cfg["filter"].get("st_duration_sec") is not None:
            raise ValueError("DAD Stage 4 uses SSS only: st_duration_sec must be null.")
        if not maxwell_cfg["enabled"]:
            raise RuntimeError("Maxwell filtering is disabled in configuration.")

        calibration_file = Path(
            maxwell_cfg["calibration_file"]
        ).expanduser().resolve()
        crosstalk_file = Path(
            maxwell_cfg["crosstalk_file"]
        ).expanduser().resolve()
        for resource in (calibration_file, crosstalk_file):
            if not resource.is_file():
                raise FileNotFoundError(
                    f"Required Maxwell resource does not exist: {resource}"
                )

        output_dirs = stage_output_dirs(cfg)
        sss_dir = ensure_output_directory(output_dirs["sss"], raw_root)
        figure_dir = ensure_output_directory(output_dirs["figures"], raw_root)
        report_dir = ensure_output_directory(output_dirs["reports"], raw_root)
        subject_id = cfg["subject"]["id"]
        decisions_out = (
            Path(args.decisions_out).expanduser().resolve()
            if args.decisions_out
            else report_dir / "03_tuning_choices.yaml"
        )
        records = discover_raw_run_records(cfg)

        # Reuse already accepted channel detection and scores. Do not rerun or
        # overwrite accepted evidence just to investigate intermittent jumps.
        accepted = None
        accepted_current = False
        if (squid_cfg["enabled"] or squid_application(cfg)) and cfg.get("stage4_review", {}).get("03_A_MaxFilt", {}).get("status") == "accepted":
            try:
                accepted = accepted_report(cfg)
                accepted_current = True
            except RuntimeError:
                # A rebuilt SSS can await reconciliation while its original
                # accepted channel basis is still usable for SQUID iteration.
                accepted = read_channel_basis(cfg, records).get("accepted_channel_report")
                if accepted is None:
                    raise RuntimeError("Accepted channel evidence unavailable; reconcile Stage 4A first.")
        accepted_rows = {r["source_filename"]: r for r in accepted["runs"]} if accepted else {}

        detected_by_run: dict[int, list[str]] = {}
        detection_details: dict[int, dict[str, Any]] = {}
        raw_headers_by_run: dict[int, mne.io.BaseRaw] = {}
        dev_head_t_by_run: dict[int, Any] = {}
        for record in records:
            if accepted:
                row = accepted_rows[record.source_filename]
                raw_header = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
                detected = list(row["detected_bad_channels"])
                details = {key: row[key] for key in ("input_header_bads", "auto_noisy_channels",
                    "auto_flat_channels", "detected_bad_channels", "digitization_points", "bad_channel_scores_path")}
                # Score paths are mount-portable by exact basename.
                details["bad_channel_scores_path"] = str(report_dir / Path(details["bad_channel_scores_path"]).name)
            else:
                detected, details, raw_header = _detect_bad_channels(
                    record, cfg=cfg, calibration_file=calibration_file, crosstalk_file=crosstalk_file)
            detected_by_run[record.acquisition_run_label] = detected
            raw_headers_by_run[record.acquisition_run_label] = raw_header
            dev_head_t_by_run[record.acquisition_run_label] = (
                raw_header.info.get("dev_head_t")
            )
            if accepted:
                scores_path = Path(details["bad_channel_scores_path"])
            else:
                scores_path = prepare_output_path(report_dir / f"03_A_MaxFilt_{args.processing_pass}_{record.source_stem}_bad_channel_scores.npz",
                                                  raw_root=raw_root, overwrite=True)
                np.savez_compressed(scores_path, **details.pop("scores"))
            details["bad_channel_scores_path"] = str(scores_path)
            detection_details[record.acquisition_run_label] = details
            # Interactive choices can use the recorded detector scores before
            # harmonization. Final report plots are rebuilt with applied labels.
            if args.interactive and not accepted:
                save_heatmaps(scores=load_scores(scores_path), source_filename=record.source_filename,
                    detected_bads=detected, applied_bads=None,
                    prefix=figure_dir / f"03_A_MaxFilt_{args.processing_pass}_{record.source_stem}",
                    raw_root=raw_root, show=not args.no_display)

        first_run = min(raw_headers_by_run)
        destination_run, destination_transform, destination_details = (
            resolve_destination(
                destination_cfg=maxwell_cfg["destination"],
                dev_head_t_by_run=dev_head_t_by_run,
                head_points=_centrality_head_points(
                    raw_headers_by_run[first_run].info
                ),
            )
        )
        destination_fingerprint = destination_transform_fingerprint(
            destination_transform
        )
        destination_details["transform"] = np.asarray(
            destination_transform["trans"],
            dtype=float,
        ).tolist()
        destination_details["transform_sha256"] = destination_fingerprint

        policy = str(maxwell_cfg["bad_channel_policy"])
        reference_run = maxwell_cfg["reference_acquisition_run"]
        manual_cfg = dict(maxwell_cfg["manual_bad_channels"])
        if args.interactive and not accepted:
            policy, reference_run, manual_cfg = interactive_bad_channel_review(
                detected_by_run=detected_by_run,
                default_policy=policy,
                default_reference_run=reference_run,
                default_manual_cfg=manual_cfg,
            )
            update_tuning_file(
                path=decisions_out,
                raw_root=raw_root,
                subject_id=subject_id,
                section="maxwell_filter",
                values={
                    "bad_channel_policy": policy,
                    "reference_acquisition_run": reference_run,
                    "manual_bad_channels": manual_cfg,
                    "destination": dict(maxwell_cfg["destination"]),
                },
            )

        cfg["maxwell_filter"].update(bad_channel_policy=policy, reference_acquisition_run=reference_run,
                                     manual_bad_channels=manual_cfg)
        selected_by_run, resolved_reference = select_bad_channels_by_policy(
            detected_by_run,
            policy=policy,
            reference_acquisition_run=reference_run,
            manual_cfg=manual_cfg,
        )

        # All reviews are validated before any SSS output can be replaced.
        # Proposals and fresh channel evidence live outside accepted reports.
        squid_reviews, screening_by_run = {}, {r: list(chs) for r, chs in selected_by_run.items()}
        if squid_cfg["enabled"] or squid_application(cfg):
            unknown = set(cfg.get("squid_jumps", {}).get("reviews", {})) - {r.source_filename for r in records}
            if unknown:
                raise RuntimeError(f"SQUID reviews contain unconfigured sources: {sorted(unknown)}")
            if accepted:
                screening_by_run = {r.acquisition_run_label: list(accepted_rows[r.source_filename]
                    .get("squid_jumps", {}).get("screening_bad_channels",
                    accepted_rows[r.source_filename]["applied_bad_channels"])) for r in records}
            remember_channel_basis(cfg, records,
                {r.source_filename: screening_by_run[r.acquisition_run_label] for r in records},
                reviewed=bool(accepted or args.interactive), report=accepted)
            pending = []
            for record in records:
                excluded = screening_by_run[record.acquisition_run_label]
                prepare_squid(cfg, record, excluded)
                try:
                    squid_reviews[record.source_filename] = load_squid_review(cfg, record, excluded)
                except (RuntimeError, FileNotFoundError) as exc:
                    pending.append(str(exc))
            if pending:
                raise RuntimeError("SQUID review required before SSS. Use review_squid_jumps.py --mode review "
                                   "with the existing tuning file. " + " | ".join(pending))
            for record in records:
                run = record.acquisition_run_label
                additions = squid_run_action(cfg, record.source_filename)["bad_channels"]
                selected_by_run[run] = _deduplicate([*screening_by_run[run], *additions])
            # Reviewed zero (or unchanged operative decisions) must not replace
            # accepted reports, SSS stamps or existing manual-eye bindings.
            compatible = bool(accepted_current)
            for record in (records if compatible else []):
                path = derivative_path_for_run(record, sss_dir, OUTPUT_SUFFIX)
                header = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
                expected = maxwell_signature(record, cfg, calibration_file, crosstalk_file,
                    destination_fingerprint, selected_by_run[record.acquisition_run_label])
                compatible = compatible and f"stage4_input_sha256={expected}" in str(header.info.get("description"))
            if compatible:
                summary = {"status": "reviewed_unchanged_sss", "reviews": cfg["squid_jumps"].get("reviews", {}),
                           "application": squid_application(cfg)}
                path = squid_root(cfg) / f"{digest(summary)}-review-summary.json"
                if not path.exists():
                    write_json(path, summary, raw_root=raw_root)
                print("SQUID review complete; operative SSS unchanged. Accepted 4A and all 4B manual reviews retained.")
                return 0
            if accepted:
                archive_squid_rebuild(cfg)

        lineage_records: list[dict[str, Any]] = []
        run_reports: list[dict[str, Any]] = []
        plot_runs: list[dict[str, Any]] = []
        reference_set = set(detected_by_run[resolved_reference])
        for record in records:
            run = record.acquisition_run_label
            requested_output = derivative_path_for_run(
                record,
                sss_dir,
                OUTPUT_SUFFIX,
            )
            execution_status = "completed"
            resume_validation = None
            expected = maxwell_signature(record, cfg, calibration_file, crosstalk_file,
                                         destination_fingerprint, selected_by_run[run])
            unchanged = False
            if squid_reviews and requested_output.exists():
                header = mne.io.read_raw_fif(requested_output, preload=False, verbose="ERROR")
                unchanged = f"stage4_input_sha256={expected}" in str(header.info.get("description"))
            if (args.resume or unchanged) and requested_output.exists():
                resume_validation = _validate_resumable_output(
                    path=requested_output,
                    record=record,
                    applied_bads=selected_by_run[run],
                    expected_signature=maxwell_signature(record, cfg, calibration_file, crosstalk_file,
                                                         destination_fingerprint, selected_by_run[run]),
                    squid_intervals=squid_run_action(cfg, record.source_filename)["intervals"],
                    destination_run=destination_run,
                    destination_transform=destination_transform,
                    destination_fingerprint=destination_fingerprint,
                    raw_root=raw_root,
                )
                output_path = requested_output.resolve()
                execution_status = "resumed_existing"
            else:
                if squid_reviews and requested_output.exists() and overwrite:
                    archive_squid_rebuild(cfg, output_path=requested_output)
                output_path = prepare_output_path(
                    requested_output,
                    raw_root=raw_root,
                    overwrite=overwrite,
                )
                _maxwell_filter_one_run(
                    record=record,
                    applied_bads=selected_by_run[run],
                    output_path=output_path,
                    cfg=cfg,
                    calibration_file=calibration_file,
                    crosstalk_file=crosstalk_file,
                    destination_run=destination_run,
                    destination_transform=destination_transform,
                    destination_fingerprint=destination_fingerprint,
                    overwrite=overwrite,
                )
            comparison_path = figure_dir / (
                f"03_A_MaxFilt_{args.processing_pass}_{record.source_stem}_before_after_psd.png"
            )
            # Scores cover the whole run; raw trace inspection reads only one
            # reproducible 20-second window, including the detected channels.
            scores_path = detection_details[run]["bad_channel_scores_path"]
            plot_runs.append(build_run_qc(raw=raw_headers_by_run[run], source_path=record.source_path, sss_path=output_path,
                scores=load_scores(scores_path), scores_path=scores_path,
                detected_bads=detected_by_run[run], applied_bads=selected_by_run[run],
                prefix=figure_dir / f"03_A_MaxFilt_{args.processing_pass}_{record.source_stem}",
                raw_root=raw_root, show=bool(args.interactive and not args.no_display)))
            detected_set = set(detected_by_run[run])
            run_report = {
                "source_filename": record.source_filename,
                "acquisition_run_label": run,
                "execution_status": execution_status,
                "resume_validation": resume_validation,
                **detection_details[run],
                "missing_vs_reference": sorted(reference_set - detected_set),
                "extra_vs_reference": sorted(detected_set - reference_set),
                "applied_bad_channels": selected_by_run[run],
                "squid_jumps": {"screening_bad_channels": screening_by_run[run],
                    "review": squid_reviews.get(record.source_filename),
                    "application": squid_run_action(cfg, record.source_filename)},
                "destination_acquisition_run": destination_run,
                "destination_transform_sha256": destination_fingerprint,
                "before_after_qc_path": str(comparison_path),
                "output_path": str(output_path),
            }
            run_reports.append(run_report)
            lineage_records.append(
                lineage_record(
                    source=record,
                    input_path=record.source_path,
                    output_path=output_path,
                    status=execution_status,
                    metadata={
                        "bad_channel_policy": policy,
                        "reference_acquisition_run": resolved_reference,
                        "destination_acquisition_run": destination_run,
                        "destination_transform_sha256": destination_fingerprint,
                        "detected_bad_channels": detected_by_run[run],
                        "applied_bad_channels": selected_by_run[run],
                        "execution_status": execution_status,
                    },
                )
            )

        report_payload = {
            "stage": STAGE_NAME,
            "generated_utc": utc_now(),
            "processing_pass": args.processing_pass,
            "interactive": bool(args.interactive),
            "resume": bool(args.resume),
            "subject_bids_id": cfg["layout"]["subject_bids_id"],
            "warnings": warnings,
            "status": "completed",
            "scientific_acceptance": "separate_explicit_stage4_review",
            "mne_version": mne.__version__,
            "calibration_file": str(calibration_file),
            "crosstalk_file": str(crosstalk_file),
            "bad_channel_policy": policy,
            "reference_acquisition_run": resolved_reference,
            "manual_bad_channels": manual_cfg,
            "destination": destination_details,
            "runs": run_reports,
        }
        highpass = cfg["preprocessing"].get("highpass_hz")
        report_payload["final_filter_design"] = filter_design(1000., highpass) if highpass is not None else None
        report_payload["source_files"] = [{"relative_path": r.raw_relative_path, "stamp": file_stamp(r.source_path)} for r in records]
        report_payload["configuration_sha256"] = stage_configuration(cfg, "03_A_MaxFilt", include_edge=False)
        bind_report(report_payload, [r["output_path"] for r in run_reports],
                    {"signatures": [maxwell_signature(r, cfg, calibration_file, crosstalk_file, destination_fingerprint, selected_by_run[r.acquisition_run_label]) for r in records],
                     "sss_configuration_sha256": stage_configuration(cfg, "03_A_MaxFilt")}, report_dir, cfg["resolved_paths"]["derivatives_root"])
        report_json = write_json(
            report_dir
            / f"03_A_MaxFilt_{args.processing_pass}_report.json",
            report_payload,
            raw_root=raw_root,
        )
        # Plot-only metadata binds these displays to the exact saved decisions;
        # it does not introduce a second scientific configuration/acceptance.
        write_json(report_dir / f"03_A_MaxFilt_{args.processing_pass}_plot_metadata.json",
                   {"stage_report_sha256": file_hash(report_json), "runs": plot_runs}, raw_root=raw_root)
        # Render saved JSON/PNGs into the primary self-contained report. This
        # presentation step opens no FIFs and makes no scientific decisions.
        report_html = render_stage_report(
            derivatives_root=Path(cfg["resolved_paths"]["derivatives_root"]),
            raw_root=raw_root,
            subject_id=cfg["layout"]["subject_bids_id"],
            processing_pass=args.processing_pass,
            stage_name=STAGE_NAME,
            overwrite=True,
            primary_report=True,
        )
        manifest_status = update_source_manifest_stage(
            manifest_path=args.source_manifest,
            raw_root=raw_root,
            stage_name=STAGE_NAME,
            processing_pass=args.processing_pass,
            records=lineage_records,
            status="completed",
            qc={
                "report_json": str(report_json),
                "report_html": str(report_html),
                "bad_channel_policy": policy,
                "destination_acquisition_run": destination_run,
                "destination_transform_sha256": destination_fingerprint,
            },
        )
        if args.print_summary:
            print(f"Maxwell filtering: completed ({len(records)} runs)")
            print(f"Bad-channel policy: {policy}")
            print(f"Reference acquisition run: {resolved_reference}")
            print(f"Destination acquisition run: {destination_run}")
            print(f"JSON report: {report_json}")
            print(f"HTML report: {report_html}")
            print(f"Source manifest: {manifest_status}")
            if args.interactive:
                print(f"Accepted tuning: {decisions_out}")
        return 0
    except Exception as exc:  # pylint: disable=broad-except
        print(f"03_A MaxFilt error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
