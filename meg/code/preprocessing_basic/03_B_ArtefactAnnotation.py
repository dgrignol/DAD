#!/usr/bin/env python3
"""
03_B_ArtefactAnnotation.py

Purpose
Detect and annotate physiological MEG artifacts for each lineage-matched
Maxwell-filtered acquisition run. Muscle detection follows the FLUX workflow.
EOG processing is selected from the subject's actual FIF modality and never
assumes that EDF files correspond one-to-one with MEG runs.

Usage examples
1) Noninteractive review using the FLUX muscle threshold of 10:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --processing-pass review --print-summary

2) Interactive threshold review:

   python .../03_B_ArtefactAnnotation.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --processing-pass review \
     --interactive --decisions-out .../03_tuning_choices.yaml

   For every acquisition run, the muscle z-score trace and threshold line are
   saved and optionally displayed. The operator may change the threshold,
   redraw the plot, and repeat until accepting it.

3) Apply reviewed thresholds:

   python .../03_B_ArtefactAnnotation.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --tuning-config .../03_tuning_choices.yaml \
     --processing-pass apply --overwrite

4) Resume an interrupted pass without replacing validated completed runs:

   python .../03_B_ArtefactAnnotation.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --tuning-config .../03_tuning_choices.yaml \
     --processing-pass review --resume

Inputs
- Layered configuration; ``--site-config`` overrides machine paths without
  changing scientific settings.
- ``<exact-source-stem>_raw_sss.fif`` from 03_A for each configured raw run.

Outputs
- ``<exact-source-stem>_ann.fif`` and ``<exact-source-stem>_ann.csv``.
- Pass-specific muscle/EOG QC PNGs and JSON/HTML reports.
- Optional accepted tuning YAML and source-manifest derivative lineage.

Data flow
raw source lineage -> expected SSS derivative -> muscle score/annotations ->
optional EOG diagnostics or annotations -> original annotations + physiological
annotations -> lineage-preserving annotated FIF and CSV.

Ocular policy
- ``auto`` becomes ``diagnostics_only`` when configured EOG channels are
  present, and ``disabled`` otherwise.
- Sub01 therefore receives EOG diagnostics but no default EOG-derived rejection
  annotations because its EOG quality is known to be limited.
- Subjects 02-07 have no EOG channels and never call EOG-only MNE functions.
- MISC001-MISC006 are preserved unchanged as acquisition-timed eye-tracker
  analog signals. They are not treated as EOG.
- Trigger 150 trial rejection is a later event/epoching decision and is not
  performed by this physiological-artifact stage.

Safety
All outputs pass the shared raw-tree guard. Existing derivatives require
explicit overwrite permission. Resume validates source/stage lineage, the
accepted muscle threshold, channel identity, bad channels, and annotation CSV
consistency before an existing run is skipped.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import annotate_muscle_zscore

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
    write_simple_html_report,
)


STAGE_NAME = "artifact_annotation"
INPUT_SUFFIX = "_raw_sss.fif"
OUTPUT_FIF_SUFFIX = "_ann.fif"
OUTPUT_CSV_SUFFIX = "_ann.csv"


# -----------------------------------------------------------------------------
# Configuration and modality helpers
# -----------------------------------------------------------------------------


def threshold_for_run(
    artifact_cfg: Mapping[str, Any],
    acquisition_run: int,
) -> float:
    """Resolve a run-specific muscle threshold with a global fallback."""
    by_run = artifact_cfg.get("muscle_threshold_by_run", {})
    candidate_keys = (
        acquisition_run,
        str(acquisition_run),
        f"{acquisition_run:02d}",
        f"run{acquisition_run:02d}",
    )
    if isinstance(by_run, Mapping):
        for key in candidate_keys:
            if key in by_run:
                return float(by_run[key])
    return float(artifact_cfg["muscle_threshold"])


def resolve_eog_mode(
    configured_mode: str,
    configured_channels: list[str],
    raw: mne.io.BaseRaw,
) -> tuple[str, list[str], list[str]]:
    """
    Resolve EOG behavior from configuration and channels actually present.

    Returns
    -------
    resolved_mode, available_channels, missing_channels
        ``auto`` resolves to ``diagnostics_only`` only when at least one
        configured EOG channel is present; otherwise it resolves to disabled.
    """
    available = [
        channel for channel in configured_channels if channel in raw.ch_names
    ]
    missing = [
        channel for channel in configured_channels if channel not in raw.ch_names
    ]
    if configured_mode == "auto":
        return (
            "diagnostics_only" if available else "disabled",
            available,
            missing,
        )
    if configured_mode in ("diagnostics_only", "annotate") and not available:
        return "disabled", available, missing
    return configured_mode, available, missing


def build_blink_annotations(
    raw: mne.io.BaseRaw,
    *,
    channel: str,
    tmin_sec: float,
    tmax_sec: float,
) -> tuple[mne.Annotations, np.ndarray]:
    """
    Detect EOG events and return annotations in the raw absolute sample frame.

    ``find_eog_events`` returns samples including ``raw.first_samp``. With
    ``orig_time=meas_date``, annotation onsets therefore retain that absolute
    frame rather than subtracting the run's first sample.
    """
    events = mne.preprocessing.find_eog_events(
        raw,
        ch_name=channel,
        verbose="ERROR",
    )
    onsets = events[:, 0] / raw.info["sfreq"] + float(tmin_sec)
    duration = float(tmax_sec) - float(tmin_sec)
    annotations = mne.Annotations(
        onset=onsets,
        duration=np.repeat(duration, len(events)),
        description=["BAD_blink"] * len(events),
        orig_time=raw.info["meas_date"],
    )
    return annotations, events


# -----------------------------------------------------------------------------
# Plotting and interactive threshold selection
# -----------------------------------------------------------------------------


def save_muscle_score_plot(
    *,
    record: AcquisitionRunRecord,
    raw: mne.io.BaseRaw,
    scores: np.ndarray,
    threshold: float,
    output_path: Path,
    show: bool,
) -> Path:
    """Save/display the FLUX-style muscle z-score decision plot."""
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(raw.times, scores, linewidth=0.6)
    ax.axhline(y=threshold, color="r")
    ax.set(
        xlabel="time (s)",
        ylabel="z-score",
        title=(
            f"{record.source_filename} | Muscle activity "
            f"(threshold = {threshold:g})"
        ),
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show(block=True)
    plt.close(fig)
    return output_path


def interactive_threshold_review(
    *,
    record: AcquisitionRunRecord,
    raw: mne.io.BaseRaw,
    scores: np.ndarray,
    initial_threshold: float,
    figure_path: Path,
    show: bool,
) -> float:
    """
    Replot and prompt until one run's muscle threshold is accepted.

    The z-score trace does not depend on the threshold, so only the decision
    line needs to be redrawn while the operator iterates.
    """
    threshold = float(initial_threshold)
    while True:
        save_muscle_score_plot(
            record=record,
            raw=raw,
            scores=scores,
            threshold=threshold,
            output_path=figure_path,
            show=show,
        )
        print(
            f"Muscle threshold plot for {record.source_filename}: {figure_path}"
        )
        answer = input(
            f"Accept muscle threshold {threshold:g}? "
            "[Y/n or enter a new numeric threshold]: "
        ).strip()
        if answer.lower() in ("", "y", "yes"):
            return threshold
        if answer.lower() in ("n", "no"):
            answer = input("New muscle threshold: ").strip()
        try:
            threshold = float(answer)
        except ValueError:
            print(f"Invalid numeric threshold: {answer!r}")


def save_eog_diagnostic_plot(
    *,
    record: AcquisitionRunRecord,
    raw: mne.io.BaseRaw,
    eog_channels: list[str],
    events: np.ndarray,
    output_path: Path,
    show: bool,
) -> Path:
    """Save/display a two-minute EOG trace with detected event markers."""
    duration = min(120.0, max(0.0, float(raw.times[-1])))
    preview = raw.copy().pick(eog_channels).crop(0.0, duration).load_data()
    data = preview.get_data()
    fig, axes = plt.subplots(
        len(eog_channels),
        1,
        figsize=(12, max(3, 2.5 * len(eog_channels))),
        sharex=True,
        squeeze=False,
    )
    event_times = (events[:, 0] - raw.first_samp) / raw.info["sfreq"]
    for index, channel in enumerate(eog_channels):
        axis = axes[index, 0]
        axis.plot(preview.times, data[index] * 1e6, linewidth=0.7)
        for event_time in event_times:
            if 0.0 <= event_time <= duration:
                axis.axvline(event_time, color="tab:red", alpha=0.25)
        axis.set_ylabel(f"{channel} (uV)")
        axis.grid(alpha=0.2)
    axes[-1, 0].set_xlabel("time (s)")
    fig.suptitle(
        f"{record.source_filename} | EOG diagnostics "
        f"({len(events)} detected events)"
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show(block=True)
    plt.close(fig)
    return output_path


# -----------------------------------------------------------------------------
# Stage execution
# -----------------------------------------------------------------------------


def _annotate_one_run(
    *,
    record: AcquisitionRunRecord,
    input_path: Path,
    output_fif: Path,
    output_csv: Path,
    artifact_cfg: Mapping[str, Any],
    initial_threshold: float,
    interactive: bool,
    show: bool,
    muscle_figure_path: Path,
    eog_figure_path: Path,
    overwrite: bool,
) -> tuple[dict[str, Any], float, list[Path]]:
    """Annotate one run and return its QC summary and accepted threshold."""
    raw = mne.io.read_raw_fif(
        input_path,
        preload=True,
        verbose="WARNING",
    )
    input_bads = list(raw.info.get("bads", []))
    original_annotations = raw.annotations.copy()
    misc_before = [
        channel
        for channel, kind in zip(
            raw.ch_names,
            raw.get_channel_types(),
            strict=True,
        )
        if kind == "misc"
    ]

    muscle_kwargs = {
        "ch_type": artifact_cfg["muscle_channel_type"],
        "min_length_good": float(
            artifact_cfg["muscle_min_length_good_sec"]
        ),
        "filter_freq": tuple(
            float(value) for value in artifact_cfg["muscle_filter_freq_hz"]
        ),
    }
    muscle_annotations, muscle_scores = annotate_muscle_zscore(
        raw,
        threshold=float(initial_threshold),
        verbose="WARNING",
        **muscle_kwargs,
    )
    accepted_threshold = float(initial_threshold)
    if interactive:
        accepted_threshold = interactive_threshold_review(
            record=record,
            raw=raw,
            scores=muscle_scores,
            initial_threshold=initial_threshold,
            figure_path=muscle_figure_path,
            show=show,
        )
        if accepted_threshold != float(initial_threshold):
            muscle_annotations, muscle_scores = annotate_muscle_zscore(
                raw,
                threshold=accepted_threshold,
                verbose="WARNING",
                **muscle_kwargs,
            )
    save_muscle_score_plot(
        record=record,
        raw=raw,
        scores=muscle_scores,
        threshold=accepted_threshold,
        output_path=muscle_figure_path,
        show=False,
    )

    configured_eog = list(artifact_cfg["eog_channel_names"])
    eog_mode, available_eog, missing_eog = resolve_eog_mode(
        str(artifact_cfg["eog_mode"]),
        configured_eog,
        raw,
    )
    blink_annotations = mne.Annotations(
        [],
        [],
        [],
        orig_time=raw.info["meas_date"],
    )
    blink_events = np.empty((0, 3), dtype=int)
    eog_status = "skipped_no_eog_channels"
    qc_images = [muscle_figure_path]
    if eog_mode in ("diagnostics_only", "annotate"):
        requested_blink_channel = artifact_cfg.get("blink_channel")
        blink_channel = (
            requested_blink_channel
            if requested_blink_channel in available_eog
            else available_eog[0]
        )
        blink_annotations, blink_events = build_blink_annotations(
            raw,
            channel=blink_channel,
            tmin_sec=float(artifact_cfg["blink_tmin_sec"]),
            tmax_sec=float(artifact_cfg["blink_tmax_sec"]),
        )
        save_eog_diagnostic_plot(
            record=record,
            raw=raw,
            eog_channels=available_eog,
            events=blink_events,
            output_path=eog_figure_path,
            show=show,
        )
        qc_images.append(eog_figure_path)
        eog_status = (
            "annotated" if eog_mode == "annotate" else "diagnostics_only"
        )

    combined = original_annotations + muscle_annotations
    if eog_mode == "annotate":
        combined = combined + blink_annotations
    raw.set_annotations(combined)
    raw.info["bads"] = list(input_bads)
    raw.info["description"] = (
        f"DAD source={record.source_filename}; "
        f"acquisition_run={record.acquisition_run_label}; "
        "stage=03_B_ArtefactAnnotation; "
        f"muscle_threshold={accepted_threshold:g}"
    )
    raw.save(output_fif, overwrite=overwrite, verbose="WARNING")
    raw.annotations.save(output_csv, overwrite=overwrite)

    misc_after = [
        channel
        for channel, kind in zip(
            raw.ch_names,
            raw.get_channel_types(),
            strict=True,
        )
        if kind == "misc"
    ]
    if misc_after != misc_before:
        raise RuntimeError(
            "MISC channel identity changed during artifact annotation: "
            f"before={misc_before}, after={misc_after}"
        )
    return (
        {
            "source_filename": record.source_filename,
            "acquisition_run_label": record.acquisition_run_label,
            "execution_status": "completed",
            "input_path": str(input_path),
            "output_fif": str(output_fif),
            "output_csv": str(output_csv),
            "input_bad_channels": input_bads,
            "original_annotation_count": len(original_annotations),
            "muscle_threshold": accepted_threshold,
            "muscle_annotation_count": len(muscle_annotations),
            "eog_mode_configured": artifact_cfg["eog_mode"],
            "eog_mode_resolved": eog_mode,
            "eog_status": eog_status,
            "available_eog_channels": available_eog,
            "missing_configured_eog_channels": missing_eog,
            "detected_eog_event_count": len(blink_events),
            "blink_annotations_applied": eog_mode == "annotate",
            "misc_channels_preserved": misc_after,
            "trigger_150_handling": "not_part_of_this_stage",
            "total_annotation_count": len(raw.annotations),
        },
        accepted_threshold,
        qc_images,
    )


def _validate_resumable_annotation(
    *,
    record: AcquisitionRunRecord,
    input_path: Path,
    output_fif: Path,
    output_csv: Path,
    threshold: float,
    artifact_cfg: Mapping[str, Any],
    raw_root: str | Path,
) -> tuple[dict[str, Any], list[Path]]:
    """
    Validate and summarize one existing annotated derivative for resume.

    Both FIF and CSV must describe the same annotation set. Input/output
    channel identity and bad-channel metadata are compared so preservation of
    MISC timing channels is checked rather than assumed.
    """
    resolved_fif = guarded_path(output_fif, raw_root)
    resolved_csv = guarded_path(output_csv, raw_root)
    if not resolved_csv.is_file():
        raise RuntimeError(
            "Cannot resume artifact annotation because its annotation CSV is "
            f"missing: {resolved_csv}"
        )
    try:
        source = mne.io.read_raw_fif(
            input_path,
            preload=False,
            verbose="ERROR",
        )
        existing = mne.io.read_raw_fif(
            resolved_fif,
            preload=False,
            verbose="ERROR",
        )
        with resolved_csv.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as csv_handle:
            csv_annotation_count = sum(1 for _row in csv.DictReader(csv_handle))
    except Exception as exc:  # pylint: disable=broad-except
        raise RuntimeError(
            f"Cannot resume from unreadable artifact derivative {resolved_fif}: "
            f"{exc}"
        ) from exc

    description = str(existing.info.get("description") or "")
    required_tokens = (
        f"source={record.source_filename}",
        f"acquisition_run={record.acquisition_run_label}",
        "stage=03_B_ArtefactAnnotation",
        f"muscle_threshold={threshold:g}",
    )
    missing_tokens = [token for token in required_tokens if token not in description]
    if missing_tokens:
        raise RuntimeError(
            "Cannot resume from artifact derivative with mismatched lineage or "
            f"threshold: {resolved_fif}; missing={missing_tokens}"
        )
    if existing.ch_names != source.ch_names:
        raise RuntimeError(
            "Cannot resume because artifact annotation changed channel identity: "
            f"{resolved_fif}"
        )
    if list(existing.info.get("bads", [])) != list(source.info.get("bads", [])):
        raise RuntimeError(
            "Cannot resume because artifact annotation changed bad-channel "
            f"metadata: {resolved_fif}"
        )
    if csv_annotation_count != len(existing.annotations):
        raise RuntimeError(
            "Cannot resume because FIF and CSV annotation counts differ: "
            f"{resolved_fif}={len(existing.annotations)}, "
            f"{resolved_csv}={csv_annotation_count}"
        )

    channel_types = existing.get_channel_types()
    misc_channels = [
        channel
        for channel, kind in zip(existing.ch_names, channel_types, strict=True)
        if kind == "misc"
    ]
    configured_eog = list(artifact_cfg["eog_channel_names"])
    eog_mode, available_eog, missing_eog = resolve_eog_mode(
        str(artifact_cfg["eog_mode"]),
        configured_eog,
        existing,
    )
    descriptions = [str(value) for value in existing.annotations.description]
    muscle_count = sum(value == "BAD_muscle" for value in descriptions)
    eog_annotation_count = sum(
        value.lower().startswith("bad_blink") for value in descriptions
    )
    return (
        {
            "source_filename": record.source_filename,
            "acquisition_run_label": record.acquisition_run_label,
            "execution_status": "resumed_existing",
            "input_path": str(input_path),
            "output_fif": str(resolved_fif),
            "output_csv": str(resolved_csv),
            "input_bad_channels": list(source.info.get("bads", [])),
            "original_annotation_count": None,
            "muscle_threshold": float(threshold),
            "muscle_annotation_count": muscle_count,
            "eog_mode_configured": artifact_cfg["eog_mode"],
            "eog_mode_resolved": eog_mode,
            "eog_status": (
                "resumed_annotated"
                if eog_mode == "annotate"
                else (
                    "resumed_diagnostics_only"
                    if eog_mode == "diagnostics_only"
                    else "skipped_no_eog_channels"
                )
            ),
            "available_eog_channels": available_eog,
            "missing_configured_eog_channels": missing_eog,
            "detected_eog_event_count": (
                eog_annotation_count if eog_mode == "annotate" else None
            ),
            "blink_annotations_applied": eog_mode == "annotate",
            "misc_channels_preserved": misc_channels,
            "trigger_150_handling": "not_part_of_this_stage",
            "total_annotation_count": len(existing.annotations),
        },
        [],
    )


# -----------------------------------------------------------------------------
# Lineage construction
# -----------------------------------------------------------------------------


def build_annotation_lineage_record(
    *,
    record: AcquisitionRunRecord,
    input_path: Path,
    run_report: Mapping[str, Any],
    accepted_threshold: float,
) -> dict[str, Any]:
    """
    Build consistent lineage for either a completed or resumed annotation run.

    Output paths come from ``run_report`` because both execution branches
    produce that common contract. This avoids referring to local output
    variables that exist only when a derivative is created during this call.
    """
    return lineage_record(
        source=record,
        input_path=input_path,
        output_path=run_report["output_fif"],
        status=run_report["execution_status"],
        metadata={
            "annotation_csv": str(run_report["output_csv"]),
            "muscle_threshold": accepted_threshold,
            "muscle_annotation_count": run_report["muscle_annotation_count"],
            "eog_status": run_report["eog_status"],
            "misc_channels_preserved": run_report["misc_channels_preserved"],
            "execution_status": run_report["execution_status"],
        },
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the artifact-annotation CLI."""
    parser = argparse.ArgumentParser(
        description="Annotate DAD MEG muscle and optional EOG artifacts."
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
    parser.add_argument("--no-display", action="store_true")
    parser.add_argument("--decisions-out", default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Validate and skip existing lineage-matched FIF/CSV pairs; fail on "
            "any mismatch. Mutually exclusive with --overwrite."
        ),
    )
    parser.add_argument("--print-summary", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run modality-aware artifact annotation for all acquisition runs."""
    args = build_arg_parser().parse_args(argv)
    try:
        cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            args.tuning_config,
            site_config_path=args.site_config,
            required_sections=("artifact_annotation",),
        )
        artifact_cfg = cfg["artifact_annotation"]
        if not artifact_cfg["enabled"]:
            raise RuntimeError("Artifact annotation is disabled in configuration.")
        raw_root = cfg["resolved_paths"]["data_root"]
        overwrite = bool(args.overwrite or cfg["runtime"]["overwrite_existing"])
        if args.resume and overwrite:
            raise ValueError("--resume and --overwrite are mutually exclusive.")
        output_dirs = stage_output_dirs(cfg)
        sss_dir = ensure_output_directory(output_dirs["sss"], raw_root)
        figure_dir = ensure_output_directory(output_dirs["figures"], raw_root)
        report_dir = ensure_output_directory(output_dirs["reports"], raw_root)
        decisions_out = (
            Path(args.decisions_out).expanduser().resolve()
            if args.decisions_out
            else report_dir / "03_tuning_choices.yaml"
        )

        raw_records = discover_raw_run_records(cfg)
        run_reports: list[dict[str, Any]] = []
        lineage_records: list[dict[str, Any]] = []
        qc_images: list[Path] = []
        accepted_by_run: dict[str, float] = {}
        for record in raw_records:
            input_path = derivative_path_for_run(record, sss_dir, INPUT_SUFFIX)
            if not input_path.is_file():
                raise FileNotFoundError(
                    "Expected lineage-matched Maxwell derivative is missing: "
                    f"{input_path}. Run 03_A_MaxFilt.py first."
                )
            requested_fif = derivative_path_for_run(
                record,
                sss_dir,
                OUTPUT_FIF_SUFFIX,
            )
            requested_csv = derivative_path_for_run(
                record,
                sss_dir,
                OUTPUT_CSV_SUFFIX,
            )
            muscle_figure = prepare_output_path(
                figure_dir
                / (
                    f"03_B_ArtefactAnnotation_{args.processing_pass}_"
                    f"{record.source_stem}_muscle_zscore.png"
                ),
                raw_root=raw_root,
                overwrite=True,
            )
            eog_figure = prepare_output_path(
                figure_dir
                / (
                    f"03_B_ArtefactAnnotation_{args.processing_pass}_"
                    f"{record.source_stem}_eog_diagnostics.png"
                ),
                raw_root=raw_root,
                overwrite=True,
            )
            initial_threshold = threshold_for_run(
                artifact_cfg,
                record.acquisition_run_label,
            )
            if args.resume and (requested_fif.exists() or requested_csv.exists()):
                if not requested_fif.is_file():
                    raise RuntimeError(
                        "Cannot resume artifact annotation because CSV exists "
                        f"without its FIF: {requested_csv}"
                    )
                run_report, run_images = _validate_resumable_annotation(
                    record=record,
                    input_path=input_path,
                    output_fif=requested_fif,
                    output_csv=requested_csv,
                    threshold=initial_threshold,
                    artifact_cfg=artifact_cfg,
                    raw_root=raw_root,
                )
                accepted_threshold = initial_threshold
                for candidate in (muscle_figure, eog_figure):
                    if candidate.is_file():
                        run_images.append(candidate)
            else:
                output_fif = prepare_output_path(
                    requested_fif,
                    raw_root=raw_root,
                    overwrite=overwrite,
                )
                output_csv = prepare_output_path(
                    requested_csv,
                    raw_root=raw_root,
                    overwrite=overwrite,
                )
                run_report, accepted_threshold, run_images = _annotate_one_run(
                    record=record,
                    input_path=input_path,
                    output_fif=output_fif,
                    output_csv=output_csv,
                    artifact_cfg=artifact_cfg,
                    initial_threshold=initial_threshold,
                    interactive=bool(args.interactive),
                    show=bool(args.interactive and not args.no_display),
                    muscle_figure_path=muscle_figure,
                    eog_figure_path=eog_figure,
                    overwrite=overwrite,
                )
            run_reports.append(run_report)
            qc_images.extend(run_images)
            accepted_by_run[f"run{record.acquisition_run_label:02d}"] = (
                accepted_threshold
            )
            lineage_records.append(
                build_annotation_lineage_record(
                    record=record,
                    input_path=input_path,
                    run_report=run_report,
                    accepted_threshold=accepted_threshold,
                )
            )

        if args.interactive:
            update_tuning_file(
                path=decisions_out,
                raw_root=raw_root,
                subject_id=cfg["subject"]["id"],
                section="artifact_annotation",
                values={"muscle_threshold_by_run": accepted_by_run},
            )

        report_payload = {
            "stage": STAGE_NAME,
            "generated_utc": utc_now(),
            "processing_pass": args.processing_pass,
            "interactive": bool(args.interactive),
            "resume": bool(args.resume),
            "subject_bids_id": cfg["layout"]["subject_bids_id"],
            "status": "completed",
            "warnings": warnings,
            "default_muscle_threshold": artifact_cfg["muscle_threshold"],
            "accepted_muscle_threshold_by_run": accepted_by_run,
            "eog_policy": (
                "Use FIF EOG only when present; EDF is not aligned or consumed."
            ),
            "trigger_150_policy": (
                "No trigger-150 rejection occurs in physiological annotation."
            ),
            "runs": run_reports,
        }
        report_json = write_json(
            report_dir
            / f"03_B_ArtefactAnnotation_{args.processing_pass}_report.json",
            report_payload,
            raw_root=raw_root,
        )
        summary_lines = [
            "status: completed",
            f"processing_pass: {args.processing_pass}",
            f"subject: {cfg['layout']['subject_bids_id']}",
            f"thresholds: {accepted_by_run}",
            "trigger 150: not handled by this stage",
            "",
            "Per-run results:",
        ]
        summary_lines.extend(
            (
                f"run{item['acquisition_run_label']:02d} "
                f"{item['source_filename']}: "
                f"execution={item['execution_status']}; "
                f"threshold={item['muscle_threshold']}; "
                f"muscle_segments={item['muscle_annotation_count']}; "
                f"eog_status={item['eog_status']}; "
                f"eog_events={item['detected_eog_event_count']}; "
                f"MISC={item['misc_channels_preserved']}"
            )
            for item in run_reports
        )
        report_html = write_simple_html_report(
            report_dir
            / f"03_B_ArtefactAnnotation_{args.processing_pass}_report.html",
            title=(
                f"{cfg['layout']['subject_bids_id']} | 03_B Artifacts "
                f"({args.processing_pass})"
            ),
            summary_lines=summary_lines,
            image_paths=qc_images,
            raw_root=raw_root,
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
                "accepted_muscle_threshold_by_run": accepted_by_run,
            },
        )
        if args.print_summary:
            print(f"Artifact annotation: completed ({len(raw_records)} runs)")
            print(f"Thresholds: {accepted_by_run}")
            print(f"JSON report: {report_json}")
            print(f"HTML report: {report_html}")
            print(f"Source manifest: {manifest_status}")
            if args.interactive:
                print(f"Accepted tuning: {decisions_out}")
        return 0
    except Exception as exc:  # pylint: disable=broad-except
        print(f"03_B ArtefactAnnotation error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
