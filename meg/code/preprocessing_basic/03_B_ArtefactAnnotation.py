#!/usr/bin/env python3
"""Stage 4B: per-acquisition muscle and manual ocular annotation.

Usage: python meg/code/preprocessing_basic/03_B_ArtefactAnnotation.py \
  --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
  --tuning-config /path/to/03_tuning_choices.yaml --processing-pass review

Requires accepted 4A evidence. A first pass may create provisional muscle scores
and annotation FIF/CSV files before manual eyes are reviewed. Use the local
manual_eye_annotations.py browser for each exact source, then explicitly rebuild
4B with --overwrite. Missing review differs from reviewed zero. Final ICA is
blocked until 4B acceptance. Threshold 10, mag, 110-140 Hz and minimum good gap
0.2 s are review starting points, not automatically accepted cohort settings.
Local SAVE can also review BAD_muscle: the saved complete muscle set replaces
automatic labels on rebuild, including reviewed zero. Automatic scores remain
diagnostic evidence; reviewed muscle masks are not fed back into their own
normalization. Legacy eye-only files retain automatic muscle detection. Changed
muscle settings require explicit reconciliation of existing muscle reviews.

Data flow: exact SSS source -> muscle scores plus original/manual BAD intervals
-> native-coordinate annotated FIF/CSV and final-filter edge-exposure preview.
Stage 4A-owned BAD_squid_jumps is preserved from the accepted SSS source once;
copied SQUID context in older manual snapshots never replaces or duplicates it.
With unchanged annotations, --resume refreshes the final-policy preview while
reusing FIF/CSV data. Superseded previews are archived as small JSON/HTML records
before report replacement. No high-pass margin is needed for spectrum-fit only.
EDF guidance is a separate non-BAD candidate annotation/CSV/mapping sidecar.
Guarded resume adds it without recomputing compatible muscle outputs. Original
recording gaps and mismatched trials remain explicit; only manual review/save
can promote edited candidates to artifact labels. Source-specific candidate
fingerprints must be reconciled before 4B acceptance. Complete manual reviews
started from --triage-proposals explicitly opt out of EDF guidance per source;
those validated receipts skip EDF preparation/reconfirmation, while all ordinary
review gates remain. Historical receipts retain their prior behavior.
No clock/lag study runs.
MISC/STI samples, native sampling, first_samp and acquisition boundaries remain.
EOG helpers below are legacy imports only; no Stage 4 EOG detection is executed.
No EDF synchronization, timing shifts, trial eligibility edits or epochs occur.
Reports embed PNGs and explain their interpretation and review limits.
Reports/tuning/lineage remain inside the existing derivative layout. Raw inputs
are read-only; --resume validates exact evidence and --overwrite is explicit.
For the explicitly authorized 2026-11-07 current-mask pass, use
stage4b_adoption.py instead of this detector/manual-only entry point. It adopts
validated exact masks, preserves saved manual work and records separate policy
provenance. Existing manual-eye receipts are never manufactured or relabelled.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import annotate_muscle_zscore

from continuous_cleaning import (stage_configuration, annotation_rows, bind_report, check_preservation, digest,
    file_hash, file_stamp, final_filter_configuration, final_filter_summary, muscle_parameters, require_review, reviewed_muscle_override)
from manual_eye_annotations import load_manual
from edf_blink_candidates import load_context, prepare_candidates, reader_api

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


STAGE_NAME = "artifact_annotation"
INPUT_SUFFIX = "_raw_sss.fif"
OUTPUT_FIF_SUFFIX = "_ann.fif"
OUTPUT_CSV_SUFFIX = "_ann.csv"


def preserve_filter_preview_history(payload, report_path, raw_root):
    """Archive superseded temporal-policy reports before refreshing 4B evidence.

    Only small report bytes are copied; annotations, scores and receipts are
    untouched. The content-addressed directory is immutable and old history
    references survive subsequent annotation-only rebuilds. The original FIR
    exposure remains historical information, never an epoch rejection rule.
    """
    report_path = guarded_path(report_path, raw_root)
    if not report_path.is_file():
        return
    previous = json.loads(report_path.read_text())
    history = list(previous.get("historical_final_filter_reports", []))
    def policies(report):
        return [final_filter_configuration(row.get("final_filter_preview", {}).get("settings", {}))
                for row in report.get("runs", [])]
    if policies(previous) != policies(payload):
        fingerprint = file_hash(report_path)
        directory = guarded_path(report_path.parent / "history/final_filter" / fingerprint, raw_root)
        directory.mkdir(parents=True, exist_ok=True)
        entry = {"json_sha256": fingerprint, "status": "historical_information_only"}
        for source in (report_path, report_path.with_suffix(".html")):
            if not source.is_file():
                continue
            target = guarded_path(directory / source.name, raw_root)
            data = source.read_bytes()
            try:
                with target.open("xb") as handle:
                    handle.write(data)
            except FileExistsError:
                if target.read_bytes() != data:
                    raise RuntimeError(f"Historical filter report changed: {target}")
            entry[source.suffix[1:]] = str(target.relative_to(report_path.parent))
        if entry not in history:
            history.append(entry)
    if history:
        payload["historical_final_filter_reports"] = history


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
    manual_annotations=None,
    manual_receipt=None,
    final_settings=None,
) -> tuple[dict[str, Any], float, list[Path]]:
    """Apply reviewed artifact edits, using automatic muscle labels as proposals.

    A missing manual review is allowed only to generate provisional 4B evidence;
    the CLI gate prevents accepting it or fitting final ICA. Existing acquisition
    annotations are preserved. Direct helper calls support synthetic smoke tests.
    A reviewed muscle replacement supersedes the whole BAD_muscle category,
    including source/automatic labels; other acquisition annotations remain.
    Automatic scores are still regenerated for inspection, independently of the
    final manually reviewed muscle set. A missing marker means legacy behavior.
    """
    muscle_override = reviewed_muscle_override(manual_receipt, artifact_cfg, initial_threshold)
    if muscle_override and manual_annotations is None:
        raise RuntimeError("Reviewed muscle replacement is missing its annotation file.")
    if muscle_override and interactive:
        raise RuntimeError("Change/reconcile the muscle threshold before applying a reviewed muscle replacement.")
    raw = mne.io.read_raw_fif(input_path, preload=True, verbose="WARNING")
    # Stage 4A SQUID decisions enter through the accepted SSS source only.
    # Editable eye/muscle snapshots cannot replace or duplicate these windows.
    if manual_annotations is not None:
        from squid_jumps import without_owned_context
        manual_annotations = without_owned_context(manual_annotations)
    source = raw.copy()
    original_annotations = raw.annotations.copy()
    if manual_annotations is not None:
        if manual_annotations.orig_time != raw.annotations.orig_time:
            raise ValueError("Manual annotations have a different acquisition origin.")
        # Exclude already reviewed severe/ocular BAD periods from muscle-score
        # normalization too; provisional scores must be rechecked after review.
        # Do not mask out the reviewed muscle intervals while measuring muscle
        # scores: that would make the diagnostic depend on its own final labels.
        score_context = manual_annotations[[i for i, label in enumerate(manual_annotations.description)
                                            if label.lower() != "bad_muscle"]]
        raw.annotations.append(score_context.onset, score_context.duration, score_context.description)
    muscle_kwargs = {"ch_type": artifact_cfg["muscle_channel_type"],
                     "min_length_good": float(artifact_cfg["muscle_min_length_good_sec"]),
                     "filter_freq": tuple(artifact_cfg["muscle_filter_freq_hz"])}
    muscle_annotations, muscle_scores = annotate_muscle_zscore(raw, threshold=initial_threshold,
                                                              verbose="WARNING", **muscle_kwargs)
    accepted_threshold = float(initial_threshold)
    if interactive:
        accepted_threshold = interactive_threshold_review(record=record, raw=raw, scores=muscle_scores,
            initial_threshold=initial_threshold, figure_path=muscle_figure_path, show=show)
        if accepted_threshold != float(initial_threshold):
            muscle_annotations, muscle_scores = annotate_muscle_zscore(raw, threshold=accepted_threshold,
                                                                       verbose="WARNING", **muscle_kwargs)
    save_muscle_score_plot(record=record, raw=raw, scores=muscle_scores, threshold=accepted_threshold,
                           output_path=muscle_figure_path, show=False)
    score_path = muscle_figure_path.with_suffix(".npz")
    np.savez_compressed(score_path, scores=muscle_scores, sfreq=raw.info["sfreq"], first_samp=raw.first_samp)
    # Muscle helper onsets are relative to raw start when meas_date is absent;
    # already attached acquisition/manual onsets include first_time. Append in
    # one common acquisition frame instead of adding first_time twice.
    muscle_onsets = muscle_annotations.onset + (raw.first_time if raw.info["meas_date"] is None else 0.)
    raw.set_annotations(None)
    preserved = (original_annotations[[i for i, label in enumerate(original_annotations.description)
                                      if label.lower() != "bad_muscle"]] if muscle_override else original_annotations)
    if muscle_override:
        # The comparison copy omits only the explicitly editable category.
        # Preserve every other source annotation and every native/auxiliary
        # sample check; the SSS source file itself is never edited.
        source.annotations.delete([i for i, label in enumerate(source.annotations.description) if label.lower() == "bad_muscle"])
    raw.annotations.append(preserved.onset, preserved.duration, preserved.description,
                           ch_names=preserved.ch_names, extras=preserved.extras)
    if not muscle_override:
        raw.annotations.append(muscle_onsets, muscle_annotations.duration, muscle_annotations.description)
    if manual_annotations is not None:
        if manual_annotations.orig_time != raw.annotations.orig_time:
            raise ValueError("Manual eye annotations have a different acquisition origin.")
        raw.annotations.append(manual_annotations.onset, manual_annotations.duration, manual_annotations.description)
    signature = artifact_signature(input_path, accepted_threshold, artifact_cfg, manual_receipt)
    raw.info["description"] = (f"DAD source={record.source_filename}; acquisition_run={record.acquisition_run_label}; "
                              f"stage=03_B_ArtefactAnnotation; muscle_threshold={accepted_threshold:g}; "
                              f"stage4_input_sha256={signature}")
    check_preservation(source, raw, annotations=False)
    raw.save(output_fif, overwrite=overwrite, verbose="WARNING")
    # FIF stores annotation endpoints as float32. Serialize the actual saved
    # times into CSV so resume compares one representation even at large
    # nonzero first_samp, rather than mistaking rounding for changed intervals.
    saved = mne.io.read_raw_fif(output_fif, preload=False, verbose="ERROR")
    check_preservation(source, saved, annotations=False)
    saved.annotations.save(output_csv, overwrite=overwrite)
    misc = [ch for ch, kind in zip(raw.ch_names, raw.get_channel_types()) if kind == "misc"]
    report = {"source_filename": record.source_filename, "acquisition_run_label": record.acquisition_run_label,
              "execution_status": "completed", "input_path": str(input_path), "input_stamp": file_stamp(input_path), "output_fif": str(output_fif),
              "output_csv": str(output_csv), "input_bad_channels": raw.info["bads"],
              "original_annotation_count": len(original_annotations), "muscle_threshold": accepted_threshold,
              "muscle_annotation_count": int(sum(label.lower() == "bad_muscle" for label in raw.annotations.description)),
              "automatic_muscle_annotation_count": len(muscle_annotations),
              "muscle_review_mode": "reviewed_replacement" if muscle_override else "automatic_with_manual_additions",
              "muscle_score_path": str(score_path),
              "manual_eye_review": manual_receipt or {"status": "unreviewed"},
              "eog_status": "legacy_not_in_stage4", "eog_mode_configured": "disabled",
              "eog_mode_resolved": "disabled", "available_eog_channels": [], "missing_configured_eog_channels": [],
              "detected_eog_event_count": 0, "blink_annotations_applied": manual_annotations is not None,
              "misc_channels_preserved": misc, "trigger_150_handling": "not_part_of_this_stage",
              "total_annotation_count": len(raw.annotations), "annotation_sha256": digest(annotation_rows(saved)),
              "stage4_input_sha256": signature}
    if final_settings is not None:
        report["final_filter_preview"] = final_filter_summary(raw, final_settings, source_filename=record.source_filename)
    return report, accepted_threshold, [muscle_figure_path]


def artifact_signature(input_path, threshold, settings, receipt):
    """Bind annotations to SSS content identity, all muscle settings and eye review."""
    return digest({"input": file_stamp(input_path), **muscle_parameters(settings, threshold),
                   "eye_review": receipt})


def _validate_resumable_annotation(
    *,
    record: AcquisitionRunRecord,
    input_path: Path,
    output_fif: Path,
    output_csv: Path,
    threshold: float,
    artifact_cfg: Mapping[str, Any],
    raw_root: str | Path,
    expected_signature=None,
    reviewed_muscle=False,
) -> tuple[dict[str, Any], list[Path]]:
    """
    Validate and summarize one existing annotated derivative for resume.

    Both FIF and CSV must describe the same annotation set. Input/output
    channel identity and bad-channel metadata are compared so preservation of
    MISC timing channels is checked rather than assumed.
    A caller-validated muscle replacement exempts only BAD_muscle from original
    annotation retention; every other preservation and signature check remains.
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

    if reviewed_muscle:
        source.annotations.delete([i for i, label in enumerate(source.annotations.description) if label.lower() == "bad_muscle"])
    check_preservation(source, existing, annotations=False)
    if expected_signature and f"stage4_input_sha256={expected_signature}" not in str(existing.info.get("description")):
        raise RuntimeError("Cannot resume artifact annotation: upstream input/settings/eye review changed.")
    # MNE cannot read an empty annotation CSV, and nonempty CSV readers may
    # rebase orig_time to the first row. Compare the actual serialized table.
    import io
    expected_csv = existing.annotations.to_data_frame().to_csv(index=False)
    expected_rows = list(csv.DictReader(io.StringIO(expected_csv)))
    with resolved_csv.open(newline="", encoding="utf-8") as handle:
        actual_rows = list(csv.DictReader(handle))
    if len(actual_rows) != len(expected_rows):
        raise RuntimeError("FIF and CSV annotation contents differ.")
    for expected, actual in zip(expected_rows, actual_rows):
        if (expected["onset"] != actual["onset"] or expected["description"] != actual["description"]
                or not np.isclose(float(expected["duration"]), float(actual["duration"]), atol=1e-6, rtol=0)):
            raise RuntimeError("FIF and CSV annotation contents differ.")
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
    eog_mode, available_eog, missing_eog = "disabled", [], []
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
            "input_stamp": file_stamp(input_path),
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
        description="Annotate DAD MEG muscle and reviewed manual ocular intervals."
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
        sss_receipt = require_review(cfg, "03_A_MaxFilt", cfg["resolved_paths"]["derivatives_reports_dir"])
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
        # Candidate sidecars do not change the muscle/SSS FIF. A guarded resume
        # can reuse muscle scores while adding fresh EDF guidance and provenance.
        # A completed targeted triage review explicitly omits EDF guidance.
        # Preserve legacy reviews unchanged; only those opt-out receipts avoid
        # EDF preparation and changed-guidance reconfirmation for their source.
        reviews = artifact_cfg.get("manual_eye_reviews", {})
        edf_needed = any(reviews.get(r.source_filename, {}).get("edf_guidance_policy") != "not_used_annotation_triage"
                         for r in raw_records)
        edf_context = load_context(cfg) if edf_needed and artifact_cfg.get("edf_blink_candidates", {}).get("enabled", False) else None
        if edf_context and edf_context["manifest"].get("sources", {}).get("eyetracker_edf_files"):
            reader_api()  # Fail on missing EDFAPI before any expensive FIF validation.
        edf_cache = {}
        edf_evidence = []
        edf_sources = {}
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
            manual_annotations, manual_receipt = None, None
            try:
                header = mne.io.read_raw_fif(input_path, preload=False, verbose="ERROR")
                manual_annotations, manual_receipt = load_manual(cfg, record, header, input_path)
            except RuntimeError as exc:
                # Only genuinely missing reviews can produce provisional reports.
                if args.processing_pass == "apply" or "Ocular review missing" not in str(exc):
                    raise
            initial_threshold = threshold_for_run(
                artifact_cfg,
                record.acquisition_run_label,
            )
            muscle_override = reviewed_muscle_override(manual_receipt, artifact_cfg, initial_threshold)
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
                    expected_signature=artifact_signature(input_path, initial_threshold, artifact_cfg, manual_receipt),
                    reviewed_muscle=muscle_override,
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
                    overwrite=overwrite, manual_annotations=manual_annotations,
                    manual_receipt=manual_receipt, final_settings=cfg["preprocessing"],
                )
            if "final_filter_preview" not in run_report:
                completed_header = mne.io.read_raw_fif(requested_fif, preload=False, verbose="ERROR")
                run_report["final_filter_preview"] = final_filter_summary(completed_header, cfg["preprocessing"], source_filename=record.source_filename)
                run_report["manual_eye_review"] = manual_receipt or {"status": "unreviewed"}
            run_report["muscle_review_mode"] = "reviewed_replacement" if muscle_override else "automatic_with_manual_additions"
            guidance = ({"status": "disabled", "review_policy": "not_used_annotation_triage"}
                        if manual_receipt and manual_receipt.get("edf_guidance_policy") == "not_used_annotation_triage"
                        else prepare_candidates(cfg, record, header, input_path, edf_context, edf_cache))
            run_report["edf_blink_candidates"] = {
                k: guidance[k] for k in ("status", "fingerprint", "candidate_count", "long_candidate_count",
                                        "matched_trial_count", "flagged_trial_count", "paths", "review_policy") if k in guidance}
            # Even a valid previous ocular review needs reconfirmation when EDF
            # guidance is newly introduced or changed. Unavailable guidance is
            # not reviewed zero; the operator must still inspect the full run.
            if guidance["status"] != "disabled":
                run_report["edf_blink_candidates"]["review_status"] = (
                    "reviewed" if manual_receipt and manual_receipt.get("edf_candidates_fingerprint") == guidance["fingerprint"] else "review_required")
                edf_evidence.extend(guidance["paths"].values())
                for source in guidance["edf_inputs"]:
                    edf_sources[source["relative_path"]] = source
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
            artifact_cfg["muscle_threshold_by_run"] = accepted_by_run
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
            "proposed_muscle_threshold_by_run": accepted_by_run,
            "scientific_acceptance": "separate_explicit_stage4_review",
            "eog_policy": (
                "EOG is legacy and excluded from Stage 4. MISC/frontal MEG support ocular review."
            ),
            "trigger_150_policy": (
                "No trigger-150 rejection occurs in physiological annotation."
            ),
            "runs": run_reports,
        }
        report_payload["status"] = "completed" if all(
            r["manual_eye_review"].get("status") == "reviewed"
            and r["edf_blink_candidates"].get("review_status", "reviewed") == "reviewed" for r in run_reports) else "review_required"
        report_payload["configuration_sha256"] = stage_configuration(cfg, "03_B_ArtefactAnnotation", include_edge=False)
        # Retain the original acceptance identity and make the new SSS-only
        # projection explicit in downstream evidence, without migrating YAML.
        report_payload["accepted_sss_lineage"] = {
            "accepted_configuration_sha256": sss_receipt["configuration_sha256"],
            "evidence_sha256": sss_receipt["evidence_sha256"],
            "effective_sss_configuration_sha256": stage_configuration(cfg, "03_A_MaxFilt"),
            "configuration_contract": "sss_independent_of_final_temporal_policy_2026-10-28"}
        if edf_context is not None:
            report_payload["source_files"] = list(edf_sources.values())
            if "trigger_derived_associations" in edf_context["manifest"]:
                report_payload["source_association_sha256"] = digest(edf_context["manifest"][
                    "trigger_derived_associations"].get("meg_acquisition_run_to_behavioral_segments", []))
            if edf_context["events_path"].exists():
                edf_evidence.append(edf_context["events_path"])
        bind_report(report_payload, [r["output_fif"] for r in run_reports] + [r["output_csv"] for r in run_reports] + edf_evidence,
                    {"thresholds": accepted_by_run, "eye_reviews": [r["manual_eye_review"] for r in run_reports],
                     "edf_guidance": [r["edf_blink_candidates"] for r in run_reports],
                     "final_filter": final_filter_configuration(cfg["preprocessing"], include_edge=False),
                     "accepted_sss_lineage": report_payload["accepted_sss_lineage"]},
                    report_dir, cfg["resolved_paths"]["derivatives_root"])
        preserve_filter_preview_history(report_payload,
            report_dir / f"03_B_ArtefactAnnotation_{args.processing_pass}_report.json", raw_root)
        report_json = write_json(
            report_dir
            / f"03_B_ArtefactAnnotation_{args.processing_pass}_report.json",
            report_payload,
            raw_root=raw_root,
        )
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
            status=report_payload["status"],
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
