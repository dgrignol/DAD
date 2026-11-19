#!/usr/bin/env python3
"""
06_export_neural_for_drsa.py

Processing runs
Export the same DAD_EPOCH_PREPARATION and DAD_PROCESSING_RUN values used for
Stage 04. Analysis/reports use the matching processing_runs/LABEL directory.
Labelled resume verifies input/settings contracts and artifact bindings.

Purpose
Export condition-specific DAD MEG tensors for dRSA from the configured step-04
epochs. Epochs are averaged by subject-specific trajectory identity on a fixed
repository-owned MEG channel axis, with bad or absent channels represented as
NaN and documented in MAT and JSON metadata.

Usage examples
1) Export subject 01:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/06_export_neural_for_drsa.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --source-manifest /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-01/ses-01/source_manifest.json

2) Validate the 40-identity and channel-axis contracts without loading data:

   python .../06_export_neural_for_drsa.py --config .../pipeline_default.yaml \
     --subject-config .../sub-04.yaml --source-manifest .../source_manifest.json \
     --validate-only --print-summary

Inputs
- Step-04 epochs FIF with complete epoch metadata.
- Step-04 machine-readable report.
- Fixed canonical MEG channel-label file.
- Exact observed trajectory MAT discovered from the configured external
  subject ``input_files`` directory. The matching ``occluded_nondeviant``
  condition supplies the expected continuation for occluded-deviant dRSA.
- Source manifest with completed epoching lineage.

Outputs
- One MAT neural tensor and JSON sidecar per condition.
- Aggregate JSON and HTML export reports.
- An ``export_drsa`` derivative-lineage section in the source manifest.

Tensor contract
- Shape: ``identity x canonical_meg_channel x time``.
- Identity axis is 1..20 or 1..40 according to the subject profile.
- No identities are synthesized or truncated. The default missing-data policy
  is an explicit error if any required condition/identity cell has no epoch.
- The canonical channel list is never initialized or modified by this stage.
- MISC, EOG, STI, and all other non-MEG channels are excluded from the tensor
  without modifying the source epochs or FIF files.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any, Mapping

import mne
import numpy as np
from scipy.io import savemat

from downstream_stage_utils import discover_trajectory_models, identity_coverage
from pipeline_config import load_and_resolve_config
from continuous_cleaning import file_hash, file_stamp
from epoch_input_preparation import (export_contract, processing_manifest,
    verify_epoch_artifacts, verify_export_artifacts)
from preprocessing_stage_utils import (
    ensure_output_directory,
    guarded_path,
    prepare_output_path,
    update_source_manifest_stage,
    utc_now,
    write_json,
    write_simple_html_report,
)


# -----------------------------------------------------------------------------
# CLI and path helpers
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the noninteractive dRSA export parser."""
    parser = argparse.ArgumentParser(
        description="Export identity-averaged DAD MEG tensors for dRSA."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-summary", action="store_true")
    return parser


def export_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return canonical epoch input, dRSA directory, and report paths."""
    resolved = cfg["resolved_paths"]
    subject = cfg["layout"]["subject_bids_id"]
    session = cfg["layout"]["session_bids_id"]
    analysis = Path(resolved["derivatives_analysis_dir"])
    reports = Path(resolved["derivatives_reports_dir"])
    project_root = Path(resolved["project_root"])
    epoch_stem = f"{subject}_{session}_task-Attention_run-all"
    canonical_value = Path(str(cfg["export"]["canonical_meg_channels_file"]))
    canonical_path = (
        canonical_value.resolve()
        if canonical_value.is_absolute()
        else (project_root / canonical_value).resolve()
    )
    return {
        "epochs_fif": analysis / "meg" / f"{epoch_stem}_epo.fif",
        "epoch_report_json": reports / "04_epoching_report.json",
        "output_dir": analysis / "drsa",
        "report_json": reports / "06_export_neural_for_drsa_report.json",
        "report_html": reports / "06_export_neural_for_drsa_report.html",
        "canonical_channels": canonical_path,
    }


def load_canonical_meg_channels(path: str | Path) -> list[str]:
    """Load a non-empty, duplicate-free fixed channel axis."""
    canonical_path = Path(path).expanduser().resolve()
    if not canonical_path.is_file():
        raise FileNotFoundError(
            "Configured canonical MEG channel file does not exist: "
            f"{canonical_path}"
        )
    channels = [
        line.strip()
        for line in canonical_path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    duplicates = sorted(
        name for name, count in Counter(channels).items() if count > 1
    )
    if not channels or duplicates:
        raise RuntimeError(
            "Canonical MEG channel axis must be non-empty and unique: "
            f"n_channels={len(channels)}, duplicates={duplicates}"
        )
    if any(not name.startswith("MEG") for name in channels):
        raise RuntimeError(
            "Canonical dRSA channel axis may contain MEG channels only."
        )
    return channels


def _mat_string_column(values: list[str]) -> np.ndarray:
    """Create a MATLAB-friendly string column with stable ordering."""
    return np.asarray(values, dtype=object).reshape(-1, 1)


def _condition_output_paths(
    cfg: Mapping[str, Any],
    output_dir: Path,
    condition: str,
) -> tuple[Path, Path]:
    """Build one condition's MAT and JSON sidecar paths."""
    subject = cfg["layout"]["subject_bids_id"]
    session = cfg["layout"]["session_bids_id"]
    stem = (
        f"{subject}_{session}_task-Attention_condition-{condition}_"
        "desc-identityAveraged_neural"
    )
    return output_dir / f"{stem}.mat", output_dir / f"{stem}.json"


# -----------------------------------------------------------------------------
# Export workflow
# -----------------------------------------------------------------------------


def run_export(
    cfg: Mapping[str, Any],
    *,
    source_manifest: str | Path,
    overwrite: bool,
    validate_only: bool = False,
    print_summary: bool = False,
) -> dict[str, Any]:
    """
    Validate the epoch contract and write one identity tensor per condition.

    Returns a complete aggregate report suitable for source-manifest QC
    references and downstream reproducibility checks.
    """
    raw_root = cfg["resolved_paths"]["data_root"]
    export_cfg = cfg["export"]
    paths = export_paths(cfg)
    canonical_channels = load_canonical_meg_channels(paths["canonical_channels"])
    trajectory_models = discover_trajectory_models(cfg)
    identity_min, identity_max = cfg["subject_profile"]["identity_range"]
    identities = list(range(int(identity_min), int(identity_max) + 1))
    event_id = {
        str(name): int(code)
        for name, code in cfg["epoching"]["event_id_map"].items()
    }

    if not paths["epochs_fif"].is_file():
        raise FileNotFoundError(f"Step-04 epochs FIF not found: {paths['epochs_fif']}")
    if not paths["epoch_report_json"].is_file():
        raise FileNotFoundError(
            f"Step-04 machine report not found: {paths['epoch_report_json']}"
        )
    epoch_report = verify_epoch_artifacts(cfg)
    if epoch_report.get("status") != "completed":
        raise RuntimeError("Step-04 report does not have status='completed'.")

    manifest_path = processing_manifest(cfg, source_manifest)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    epoch_lineage = (
        manifest.get("derivative_lineage", {})
        .get("stages", {})
        .get("epoching")
    )
    if not isinstance(epoch_lineage, Mapping) or epoch_lineage.get("status") != "completed":
        raise RuntimeError(
            "Source manifest must contain completed epoching lineage before export."
        )

    for key in ("output_dir", "report_json", "report_html"):
        guarded_path(paths[key], raw_root)
    condition_paths = {
        condition: _condition_output_paths(cfg, paths["output_dir"], condition)
        for condition in event_id
    }
    for mat_path, json_path in condition_paths.values():
        guarded_path(mat_path, raw_root)
        guarded_path(json_path, raw_root)

    setup = {
        "stage": "export_drsa",
        "epoch_contract_sha256": epoch_report.get("preparation_contract", {}).get("sha256"),
        "export_contract_sha256": export_contract(cfg, epoch_report, canonical_channels, trajectory_models),
        "status": "validated" if validate_only else "running",
        "generated_utc": utc_now(),
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "identity_values": identities,
        "identity_count": len(identities),
        "event_id_map": event_id,
        "trajectory_models": trajectory_models,
        "epoch_input": str(paths["epochs_fif"]),
        "epoch_report": str(paths["epoch_report_json"]),
        "source_manifest": str(manifest_path),
        "canonical_meg_channels_file": str(paths["canonical_channels"]),
        "n_canonical_meg_channels": len(canonical_channels),
        "export_window_sec": [
            float(export_cfg["tmin_sec"]),
            float(export_cfg["tmax_sec"]),
        ],
        "sampling_rate_hz": float(export_cfg["sampling_rate_hz"]),
        "missing_identity_policy": str(export_cfg["missing_identity_policy"]),
        "bad_or_missing_channel_policy": str(
            export_cfg["bad_or_missing_channel_policy"]
        ),
        "inclusion_policy": {
            "epoch_stage_trial_exclusion": epoch_report["trial_exclusion_policy"],
            "epoch_stage_annotation_rejection": epoch_report["annotation_policy"],
            "export": (
                "All surviving epochs are averaged within condition and "
                "trial_sequence_identity. No additional trial rejection occurs."
            ),
        },
        "condition_outputs": {
            condition: {"mat": str(values[0]), "json": str(values[1])}
            for condition, values in condition_paths.items()
        },
    }
    if print_summary:
        print(json.dumps(setup, indent=2, sort_keys=True))
    if validate_only:
        return setup

    if cfg["epoching"].get("processing_run") and (paths["report_json"].exists()
            or any(path.exists() for pair in condition_paths.values() for path in pair)):
        raise RuntimeError("Labelled exports already exist; resume matching outputs or use a new run label.")

    ensure_output_directory(paths["output_dir"], raw_root)
    ensure_output_directory(paths["report_json"].parent, raw_root)
    prepare_output_path(
        paths["report_json"],
        raw_root=raw_root,
        overwrite=overwrite,
    )
    prepare_output_path(
        paths["report_html"],
        raw_root=raw_root,
        overwrite=overwrite,
    )
    for mat_path, json_path in condition_paths.values():
        prepare_output_path(mat_path, raw_root=raw_root, overwrite=overwrite)
        prepare_output_path(json_path, raw_root=raw_root, overwrite=overwrite)

    epochs = mne.read_epochs(paths["epochs_fif"], preload=False, verbose=False)
    if epochs.metadata is None:
        raise RuntimeError("Epochs file has no metadata.")
    required_metadata = {
        "trial_sequence_identity",
        "condition_name",
        "source_fif",
        "acquisition_run",
        "raw_sample_frame",
        "exclusion_reasons",
        "is_replay_trial",
    }
    missing_metadata = required_metadata.difference(epochs.metadata.columns)
    if missing_metadata:
        raise RuntimeError(
            "Epoch metadata is missing export contract columns: "
            f"{sorted(missing_metadata)}"
        )
    if {
        name: int(code) for name, code in epochs.event_id.items()
    } != event_id:
        raise RuntimeError(
            f"Epoch event map mismatch: expected={event_id}, observed={epochs.event_id}"
        )
    if float(export_cfg["tmin_sec"]) < float(epochs.tmin) or float(
        export_cfg["tmax_sec"]
    ) > float(epochs.tmax):
        raise RuntimeError(
            "Configured export window lies outside epoch boundaries: "
            f"export={[export_cfg['tmin_sec'], export_cfg['tmax_sec']]}, "
            f"epochs={[epochs.tmin, epochs.tmax]}"
        )

    meg_picks = mne.pick_types(
        epochs.info,
        meg=True,
        eeg=False,
        eog=False,
        stim=False,
        misc=False,
        exclude=[],
    )
    available_meg = [epochs.ch_names[index] for index in meg_picks]
    unexpected_meg = sorted(set(available_meg).difference(canonical_channels))
    if unexpected_meg:
        raise RuntimeError(
            "Epochs contain MEG channels absent from the fixed canonical axis: "
            f"{unexpected_meg}"
        )
    present_canonical = [
        channel for channel in canonical_channels if channel in available_meg
    ]
    missing_canonical = [
        channel for channel in canonical_channels if channel not in available_meg
    ]
    marked_bad = [
        channel
        for channel in canonical_channels
        if channel in epochs.info.get("bads", [])
    ]
    signal_channels = [
        channel for channel in present_canonical if channel not in marked_bad
    ]
    if not signal_channels:
        raise RuntimeError("No usable canonical MEG channels remain for export.")

    condition_reports: list[dict[str, Any]] = []
    for condition, condition_code in event_id.items():
        condition_epochs = epochs[condition].copy()
        identities_by_epoch = condition_epochs.metadata[
            "trial_sequence_identity"
        ].to_numpy(dtype=int)
        try:
            counts, missing_identities = identity_coverage(
                identities_by_epoch,
                identities,
                missing_policy=str(export_cfg["missing_identity_policy"]),
            )
        except RuntimeError as exc:
            raise RuntimeError(f"Condition {condition}: {exc}") from exc

        condition_epochs.load_data()
        condition_epochs.pick(present_canonical)
        condition_epochs.reorder_channels(present_canonical)
        condition_epochs.crop(
            tmin=float(export_cfg["tmin_sec"]),
            tmax=float(export_cfg["tmax_sec"]),
        )
        condition_epochs.resample(
            float(export_cfg["sampling_rate_hz"]),
            npad="auto",
        )
        data = condition_epochs.get_data(copy=False)
        tensor = np.full(
            (
                len(identities),
                len(canonical_channels),
                data.shape[2],
            ),
            np.nan,
            dtype=np.float64,
        )
        present_index = {
            name: index for index, name in enumerate(condition_epochs.ch_names)
        }
        canonical_index = {
            name: index for index, name in enumerate(canonical_channels)
        }
        for identity_axis, identity in enumerate(identities):
            epoch_indices = np.flatnonzero(identities_by_epoch == identity)
            if len(epoch_indices) == 0:
                continue
            identity_mean = data[epoch_indices].mean(axis=0)
            for channel in signal_channels:
                tensor[
                    identity_axis,
                    canonical_index[channel],
                    :,
                ] = identity_mean[present_index[channel], :]

        mat_path, sidecar_path = condition_paths[condition]
        condition_report = {
            "condition_name": condition,
            "condition_code": condition_code,
            "identity_values": identities,
            "counts_by_identity": counts,
            "missing_identities": missing_identities,
            "tensor_shape": [int(value) for value in tensor.shape],
            "ordered_channel_labels": canonical_channels,
            "time_axis_sec": [float(value) for value in condition_epochs.times],
            "sampling_rate_hz": float(condition_epochs.info["sfreq"]),
            "marked_bad_channels": marked_bad,
            "missing_canonical_channels": missing_canonical,
            "nan_channel_labels": sorted(set(marked_bad + missing_canonical)),
            "bad_or_missing_channel_policy": export_cfg[
                "bad_or_missing_channel_policy"
            ],
            "excluded_non_meg_channel_types": ["misc", "eog", "stim", "ecg", "eeg"],
            "trajectory_models": trajectory_models,
            "event_id_map": event_id,
            "source_lineage": {
                "epochs_fif": str(paths["epochs_fif"]),
                "epoch_report": str(paths["epoch_report_json"]),
                "source_manifest": str(manifest_path),
            },
            "inclusion_policy": setup["inclusion_policy"],
            "mat_path": str(mat_path),
            "json_path": str(sidecar_path),
        }
        savemat(
            mat_path,
            {
                "neural": tensor,
                "channel_labels": _mat_string_column(canonical_channels),
                "time_axis_sec": np.asarray(condition_epochs.times, dtype=float),
                "sampling_rate_hz": np.asarray(
                    [[condition_epochs.info["sfreq"]]],
                    dtype=float,
                ),
                "identity_values": np.asarray(identities, dtype=int).reshape(-1, 1),
                "condition_name": np.asarray([[condition]], dtype=object),
                "condition_code": np.asarray([[condition_code]], dtype=int),
                "bad_channel_labels": _mat_string_column(marked_bad),
                "missing_channel_labels": _mat_string_column(missing_canonical),
                "bad_or_missing_channel_policy": np.asarray(
                    [[export_cfg["bad_or_missing_channel_policy"]]],
                    dtype=object,
                ),
                "observed_model_path": np.asarray(
                    [[trajectory_models["observed"]]],
                    dtype=object,
                ),
                "expected_model_source_condition": np.asarray(
                    [[trajectory_models["expected_source_condition"]]],
                    dtype=object,
                ),
            },
            do_compression=True,
        )
        write_json(
            sidecar_path,
            condition_report,
            raw_root=raw_root,
            overwrite=True,
        )
        condition_reports.append(condition_report)

    report = {
        **setup,
        "status": "completed",
        "finished_utc": utc_now(),
        "artifact_bindings": [],
        "channel_axis": {
            "ordered_channel_labels": canonical_channels,
            "n_channels": len(canonical_channels),
            "present_channels": present_canonical,
            "missing_channels": missing_canonical,
            "marked_bad_channels": marked_bad,
            "signal_channels": signal_channels,
            "non_meg_channels_omitted_without_source_modification": True,
        },
        "conditions": condition_reports,
    }
    # Bind all condition files after shared storage exposes stable timestamps.
    from review_final_ica import settle_review_outputs
    output_files = [path for pair in condition_paths.values() for path in pair]
    settle_review_outputs(output_files)
    report["artifact_bindings"] = [{"stamp": file_stamp(path),
        "sha256": file_hash(path) if path.suffix == ".json" else None} for path in output_files]
    write_json(paths["report_json"], report, raw_root=raw_root, overwrite=True)
    write_simple_html_report(
        paths["report_html"],
        title=f"{cfg['layout']['subject_bids_id']} | 06 dRSA Export",
        summary_lines=[
            "Stage status: completed",
            f"Identity count: {len(identities)}",
            f"Canonical MEG channels: {len(canonical_channels)}",
            f"Marked bad channels exported as NaN: {marked_bad}",
            f"Missing channels exported as NaN: {missing_canonical}",
            f"Conditions: {list(event_id)}",
            f"Machine report: {paths['report_json']}",
        ],
        image_paths=[],
        raw_root=raw_root,
    )
    manifest_status = update_source_manifest_stage(
        manifest_path=manifest_path,
        raw_root=raw_root,
        stage_name="export_drsa",
        processing_pass="apply",
        records=[
            {
                "input_path": str(paths["epochs_fif"]),
                "output_path": item["mat_path"],
                "status": "completed",
                "condition_name": item["condition_name"],
                "condition_code": item["condition_code"],
                "tensor_shape": item["tensor_shape"],
            }
            for item in condition_reports
        ],
        status="completed",
        qc={
            "report_json": str(paths["report_json"]),
            "report_html": str(paths["report_html"]),
            "identity_count": len(identities),
            "canonical_meg_channels_file": str(paths["canonical_channels"]),
            "bad_or_missing_channel_policy": export_cfg[
                "bad_or_missing_channel_policy"
            ],
        },
    )
    report["source_manifest_update"] = manifest_status
    write_json(paths["report_json"], report, raw_root=raw_root, overwrite=True)
    if print_summary:
        print(
            json.dumps(
                {
                    "identity_count": len(identities),
                    "conditions": [
                        {
                            "name": item["condition_name"],
                            "tensor_shape": item["tensor_shape"],
                            "missing_identities": item["missing_identities"],
                        }
                        for item in condition_reports
                    ],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return report


def main() -> int:
    """Resolve configuration and execute the dRSA export stage."""
    args = build_arg_parser().parse_args()
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        args.subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("epoching", "export"),
    )
    if warnings and args.print_summary:
        print(json.dumps({"configuration_warnings": warnings}, indent=2))
    paths = export_paths(cfg)
    if args.resume and paths["report_json"].is_file():
        if cfg["epoching"].get("processing_run"):
            current = run_export(cfg, source_manifest=args.source_manifest, overwrite=False, validate_only=True)
            previous = verify_export_artifacts(cfg)
            if previous["export_contract_sha256"] != current["export_contract_sha256"]:
                raise RuntimeError("Export inputs/settings changed; use a new processing-run label.")
        print(f"Resumed existing dRSA export report: {paths['report_json']}")
        return 0
    run_export(
        cfg,
        source_manifest=args.source_manifest,
        overwrite=bool(args.overwrite or cfg["runtime"]["overwrite_existing"]),
        validate_only=args.validate_only,
        print_summary=args.print_summary,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
