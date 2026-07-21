#!/usr/bin/env python3
"""
01_dot_side_decoding.py

Purpose
Run a MEG sanity check that decodes dot horizontal side (left vs right half of
screen) with centered sliding windows for one neural condition and trajectory
model source, using either the historical posterior-gradiometer feature set or
a named Stage-07 sensor group, and generate a 3-panel figure:
1) baseline-corrected contra/ipsi ERFs from posterior ROIs (auto-selected or
   explicitly overridden by user-provided magnetometer channel names),
2) ipsi-minus-contra waveform,
3) time-resolved sliding-window decoding accuracy curves for unbalanced and
   class-balanced (undersampled) decoding. Occluded-condition decoding panels
   mark the first completely occluded frame and first reappearance frame from
   the shared pipeline timing configuration.

This script is intended for the DAD MEG sanity-check workflow and writes all
outputs into a dedicated sanity-check derivatives folder.

Usage examples
1) Default run for sub-01, ses-01, always_visible:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml

2) Run another subject/session:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --subject 02 --session 01 --condition always_visible

3) Validate a completed Stage-04 handoff without decoding or writing outputs:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-05.yaml \
     --subject 05 --validate-only --print-summary

4) Run 10-sample sliding-window decoding with same-identity averaging:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --subject 05 --condition always_visible \
     --decode-window-samples 10 --trials-per-average 2 --random-seed 97 \
     --overwrite

5) Average all retained repetitions of each identity before decoding:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --subject 05 --condition always_visible \
     --decode-window-samples 10 --trials-per-average all_trials_per_identity \
     --overwrite

6) Decode occluded-deviant epochs against their expected/predicted path:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --subject 05 --condition occluded_deviant --model-source predicted \
     --decode-window-samples 10 --trials-per-average 2 --overwrite

7) Decode a Stage-07 lobe sensor group for poster/group summaries:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --subject 05 --condition occluded_deviant --model-source observed \
     --sensor-group occipital --decoder-sensor-type grad \
     --decode-window-samples 10 --trials-per-average 2 --overwrite

8) Explicit input/output overrides:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --epochs-fif /abs/path/to/sub-01_ses-01_task-Attention_run-all_epo.fif \
     --model-mat /abs/path/to/MovDot_SubNN_VXX_rescueTraject.mat \
     --out-dir /abs/path/to/output/folder

9) Force specific ERF ROIs:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/01_dot_side_decoding.py \
     --left-roi-mag-channels MEG1731 \
     --right-roi-mag-channels MEG2241

Inputs
- Epoched FIF containing the requested condition and metadata column
  `trial_sequence_identity`.
- Model MAT with `xySeqs` entries containing `condition_label`, `sequence`,
  and `xy`, plus `Cfg.rectSize` for center-threshold definition.

Outputs
- Figure PNG: three-panel sanity-check plot.
- Scores CSV: per-timepoint unbalanced and balanced decoding accuracies.
- Summary JSON: settings, channel sets, sliding-window metadata, trial
  averaging source-index metadata, balancing strategy, and metrics.
- MNE HTML report: summary text + embedded three-panel figure + whole-sensor
  ERF layout figures for magnetometers and RMS-merged planar gradiometers.

Key assumptions
- `trial_sequence_identity` spans the configured 1..20 or 1..40 profile range
  and maps to the exact subject-specific observed trajectory model.
- Predicted decoding is valid only for `occluded_deviant`; it uses deviant
  neural epochs and the matching `occluded_nondeviant` trajectory from the
  same observed MAT file as the expected continuation.
- Dot side label is `left` when model x < `Cfg.rectSize(1)/2`, else `right`.
- `--sensor-group posterior` preserves the historical automatic posterior
  gradiometer selection. Stage-07 groups use the shared MATLAB lobe channel
  definitions and then filter to `--decoder-sensor-type`.
- Sliding windows are interpreted after `--decode-sfreq` resampling. Edge
  windows are clipped and every effective window is recorded.
- Trial averaging is deterministic, same-identity only, and drops incomplete
  remainder groups in fixed-count mode so pseudo-trials have comparable SNR.
  The literal value `all_trials_per_identity` instead creates one evoked
  pseudo-trial per retained identity.
- Balanced decoding uses random undersampling (majority -> minority count)
  inside each CV train/test fold; unbalanced decoding uses original fold data.
"""

from __future__ import annotations

import argparse
from collections import Counter
import html
import json
from pathlib import Path
import re
import sys
from typing import Any, Mapping

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
import scipy.io as sio
from mne.io.constants import FIFF
from mne.transforms import apply_trans


PREPROCESSING_DIR = Path(__file__).resolve().parents[1] / "preprocessing_basic"
sys.path.insert(0, str(PREPROCESSING_DIR))

from downstream_stage_utils import discover_trajectory_models  # noqa: E402
from drsa_result_utils import occlusion_timing_from_config  # noqa: E402
from pipeline_config import (  # noqa: E402
    assert_output_path_outside_raw,
    load_and_resolve_config,
)


MODEL_SOURCES = ("observed", "predicted")
ALL_TRIALS_PER_IDENTITY = "all_trials_per_identity"
FIXED_COUNT_AVERAGING_MODE = "fixed_count"
ALL_TRIALS_AVERAGING_MODE = "all_trials_per_identity"
DECODER_SENSOR_TYPES = ("grad", "mag")
STAGE07_SENSOR_GROUPS = ("wholebrain", "occipital", "parietal", "frontal", "temporal")
DECODING_SENSOR_GROUPS = ("posterior", *STAGE07_SENSOR_GROUPS)
MEGIN_SENSOR_GROUP_RE = re.compile(r"sensor_groups\.(\w+)\s*=\s*\{([^}]*)\};")
MEG_CHANNEL_RE = re.compile(r"'(MEG\d{4})'")
DECODER_MODEL_NAME = "Regularized linear discriminant analysis"
DECODER_MODEL_DETAILS = (
    "Train-fold z-scored binary LDA with a pooled covariance estimate and "
    "ridge regularization reg=1e-2. The Woodbury equivalent is used when the "
    "windowed feature dimension exceeds the available training samples."
)


# -----------------------------------------------------------------------------
# CLI and path resolution
# -----------------------------------------------------------------------------

def parse_trials_per_average(value: str | int) -> int | str:
    """
    Parse the trial-averaging request from CLI or tests.

    Accepted inputs are positive integers for fixed-count same-identity
    pseudo-trial groups, or the literal ``all_trials_per_identity`` for one
    evoked pseudo-trial per retained trajectory identity.
    """
    text = str(value).strip()
    if text == ALL_TRIALS_PER_IDENTITY:
        return ALL_TRIALS_PER_IDENTITY
    try:
        parsed = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--trials-per-average must be a positive integer or "
            f"{ALL_TRIALS_PER_IDENTITY!r}."
        ) from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("--trials-per-average must be at least 1.")
    return parsed


def normalize_trial_averaging_request(value: str | int) -> tuple[str, int | None]:
    """
    Return ``(mode, fixed_count)`` for a parsed trial-averaging request.

    ``fixed_count`` is an integer only for fixed-count mode; all-trials mode
    records ``None`` because the per-identity source count determines group
    size separately for every identity.
    """
    parsed = parse_trials_per_average(value)
    if parsed == ALL_TRIALS_PER_IDENTITY:
        return ALL_TRIALS_AVERAGING_MODE, None
    return FIXED_COUNT_AVERAGING_MODE, int(parsed)


def format_trial_averaging_request(value: str | int | None) -> str:
    """Create a readable representation of a trial-averaging request."""
    if value is None:
        return ALL_TRIALS_PER_IDENTITY
    if value == ALL_TRIALS_PER_IDENTITY:
        return ALL_TRIALS_PER_IDENTITY
    return str(int(value))


def build_arg_parser() -> argparse.ArgumentParser:
    """Create CLI parser for dot-side decoding sanity checks."""
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Dot-side decoding sanity check")
    parser.add_argument(
        "--config",
        default=str(project_root / "params" / "pipeline_default.yaml"),
        help="Shared pipeline default YAML.",
    )
    parser.add_argument("--site-config", default=None)
    parser.add_argument(
        "--subject-config",
        default=None,
        help="Subject selector YAML; defaults from --subject.",
    )
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--subject", type=str, default="01", help="Subject id without 'sub-' prefix.")
    parser.add_argument("--session", type=str, default="01", help="Session id without 'ses-' prefix.")
    parser.add_argument(
        "--condition",
        type=str,
        default="always_visible",
        help="Condition name in epochs.event_id and model condition_label.",
    )
    parser.add_argument(
        "--model-source",
        choices=MODEL_SOURCES,
        default="observed",
        help=(
            "Trajectory used to derive dynamic dot-side labels. Predicted is "
            "valid only for occluded_deviant and uses the matching "
            "occluded_nondeviant trajectory from the same MAT file."
        ),
    )
    parser.add_argument(
        "--sensor-group",
        default="posterior",
        help=(
            "Decoder sensor group. 'posterior' preserves the historical "
            "auto-selected posterior gradiometer decoder; Stage-07 groups "
            "include wholebrain, occipital, parietal, frontal, and temporal."
        ),
    )
    parser.add_argument(
        "--decoder-sensor-type",
        choices=DECODER_SENSOR_TYPES,
        default="grad",
        help="MEG sensor type to use as decoder features for Stage-07 groups.",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="MEG project root (default inferred from script location).",
    )
    parser.add_argument(
        "--dad-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="DAD root containing experiment input files.",
    )
    parser.add_argument(
        "--epochs-fif",
        type=Path,
        default=None,
        help="Optional explicit epochs FIF path. If omitted, derived from subject/session.",
    )
    parser.add_argument(
        "--model-mat",
        type=Path,
        default=None,
        help=(
            "Optional explicit trajectory MAT path. The subject-specific "
            "observed MAT contains both observed and expected trajectories."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Optional explicit output directory. If omitted, uses derivatives/sanity_checks.",
    )
    parser.add_argument(
        "--decode-sfreq",
        type=float,
        default=120.0,
        help="Sampling frequency used for decoding after resampling.",
    )
    parser.add_argument(
        "--decode-tmin",
        type=float,
        default=0.0,
        help="First decoded center time in seconds after condition onset.",
    )
    parser.add_argument(
        "--decode-tmax",
        type=float,
        default=2.667,
        help="Last decoded center time in seconds after condition onset.",
    )
    parser.add_argument(
        "--decode-window-samples",
        type=int,
        default=10,
        help=(
            "Centered sliding-window width in samples after --decode-sfreq "
            "resampling. Edge windows are clipped to available samples."
        ),
    )
    parser.add_argument(
        "--trials-per-average",
        type=parse_trials_per_average,
        default=1,
        help=(
            "Number of same-identity source trials to average into each "
            "pseudo-trial before classification. Incomplete remainder groups "
            "are dropped. Use the literal all_trials_per_identity to average "
            "all retained repetitions of each identity into one pseudo-trial."
        ),
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
        help="Seed for same-identity trial shuffling and decoding splits.",
    )
    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Stratified CV folds for time-resolved decoding.",
    )
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=None,
        help=(
            "Optional source_manifest.json used to validate that Stage 04 "
            "epoching completed for the consumed epochs file."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse existing sanity-check outputs when all expected files exist.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Intentionally replace existing sanity-check outputs.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs, lineage, and output contract without decoding.",
    )
    parser.add_argument(
        "--print-summary",
        action="store_true",
        help="Print resolved inputs/outputs during validation.",
    )
    parser.add_argument(
        "--left-roi-mag-channels",
        type=str,
        default=None,
        help=(
            "Comma-separated left ROI magnetometer channels for ERFs "
            "(example: MEG1731 or MEG1731,MEG1711). "
            "If omitted, channels are auto-selected."
        ),
    )
    parser.add_argument(
        "--right-roi-mag-channels",
        type=str,
        default=None,
        help=(
            "Comma-separated right ROI magnetometer channels for ERFs "
            "(example: MEG2241 or MEG2241,MEG2231). "
            "If omitted, channels are auto-selected."
        ),
    )
    return parser


def resolve_input_output_paths(
    args: argparse.Namespace,
) -> tuple[Path, Path, Path, Path, dict]:
    """
    Resolve and validate epochs/model inputs and output directory.

    Data flow
    - Build default paths from subject/session identifiers.
    - Respect explicit CLI path overrides when provided.
    - Validate required input files exist before any analysis runs.
    """
    project_root = Path(args.project_root).expanduser().resolve()
    subject_config = (
        Path(args.subject_config).expanduser().resolve()
        if args.subject_config
        else project_root / "params" / f"sub-{int(args.subject):02d}.yaml"
    )
    cfg, _warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("epoching",),
    )
    subject = str(cfg["subject"]["id"])
    session = str(cfg["layout"]["session_bids_id"]).removeprefix("ses-")
    resolved_analysis = Path(cfg["resolved_paths"]["derivatives_analysis_dir"])
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    sensor_group = normalize_sensor_group(args.sensor_group)

    default_epochs = (
        resolved_analysis
        / "meg"
        / f"sub-{subject}_ses-{session}_task-Attention_run-all_epo.fif"
    )
    default_model = Path(discover_trajectory_models(cfg)["observed"])
    default_base = (
        derivatives_root
        / "sanity_checks"
        / f"sub-{subject}"
        / f"ses-{session}"
        / "dot_side_decoding"
    )
    if sensor_group == "posterior":
        default_out = default_base / str(args.condition) / str(args.model_source)
    else:
        default_out = default_base / sensor_group / str(args.condition) / str(args.model_source)
    default_source_manifest = (
        Path(cfg["resolved_paths"]["derivatives_reports_dir"])
        / "source_manifest.json"
    )

    epochs_fif = Path(args.epochs_fif).expanduser().resolve() if args.epochs_fif else default_epochs
    model_mat = Path(args.model_mat).expanduser().resolve() if args.model_mat else default_model
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else default_out
    source_manifest = (
        Path(args.source_manifest).expanduser().resolve()
        if args.source_manifest
        else default_source_manifest
    )

    if not epochs_fif.exists():
        raise FileNotFoundError(f"Epochs FIF not found: {epochs_fif}")
    if not model_mat.exists():
        raise FileNotFoundError(f"Model MAT not found: {model_mat}")
    if not source_manifest.exists():
        raise FileNotFoundError(f"Source manifest not found: {source_manifest}")

    assert_output_path_outside_raw(out_dir, cfg["resolved_paths"]["data_root"])
    return epochs_fif, model_mat, out_dir, source_manifest, cfg


def analysis_key(condition: str, model_source: str) -> str:
    """Return the stable condition/model-source identifier used in outputs."""
    return f"{condition}_{model_source}"


def normalize_sensor_group(value: str) -> str:
    """Normalize decoder sensor-group names used on the CLI and in paths."""
    text = str(value).strip().lower().replace("_", "")
    if text in {"wholebrain", "whole"}:
        return "wholebrain"
    if text == "pariatel":
        return "parietal"
    if text in DECODING_SENSOR_GROUPS:
        return text
    raise ValueError(
        f"Unknown --sensor-group {value!r}; expected one of "
        f"{', '.join(DECODING_SENSOR_GROUPS)}."
    )


def stage07_sensor_group_definition_path() -> Path:
    """Return the MATLAB utility containing the canonical Stage-07 groups."""
    return Path(__file__).resolve().parents[1] / "analysis" / "drsa_stage_utils.m"


def load_stage07_sensor_groups(path: Path | None = None) -> dict[str, list[str]]:
    """
    Load Stage-07 MEGIN/Neuromag sensor groups from the MATLAB utility.

    The decoder uses the same lobe definitions as dRSA/ERF exports. Each group
    contains mixed magnetometer and planar gradiometer names; feature selection
    later filters that list to ``--decoder-sensor-type``.
    """
    definition_path = path or stage07_sensor_group_definition_path()
    text = definition_path.read_text(encoding="utf-8")
    groups: dict[str, list[str]] = {}
    for raw_group, body in MEGIN_SENSOR_GROUP_RE.findall(text):
        group = normalize_sensor_group(raw_group)
        if group == "posterior":
            continue
        channels = MEG_CHANNEL_RE.findall(body)
        if channels:
            groups[group] = channels
    for group in ("frontal", "temporal", "parietal", "occipital"):
        if group not in groups:
            raise RuntimeError(
                f"Stage-07 sensor group {group!r} was not found in {definition_path}."
            )
    wholebrain: list[str] = []
    for group in ("frontal", "temporal", "parietal", "occipital"):
        for channel in groups[group]:
            if channel not in wholebrain:
                wholebrain.append(channel)
    groups["wholebrain"] = wholebrain
    return groups


def resolve_model_condition(condition: str, model_source: str) -> str:
    """
    Resolve which trajectory condition supplies labels for one neural condition.

    Observed analyses use the matching condition. Predicted decoding is
    intentionally restricted to occluded-deviant neural epochs and uses the
    corresponding occluded-nondeviant path as the expected continuation.
    """
    if model_source not in MODEL_SOURCES:
        raise ValueError(
            f"Unsupported model source {model_source!r}; expected one of {MODEL_SOURCES}."
        )
    if model_source == "predicted":
        if condition != "occluded_deviant":
            raise ValueError(
                "Predicted model decoding is valid only for condition "
                "'occluded_deviant'."
            )
        return "occluded_nondeviant"
    return condition


def sanity_output_paths(
    out_dir: Path,
    condition: str,
    model_source: str,
) -> dict[str, Path]:
    """Return all expected dot-side sanity-check outputs."""
    stem = analysis_key(condition, model_source)
    return {
        "figure": out_dir / f"{stem}_dot_side_decoding_3panel.png",
        "scores_csv": out_dir / f"{stem}_dot_side_decoding_scores.csv",
        "summary_json": out_dir / f"{stem}_dot_side_decoding_summary.json",
        "report_html": out_dir / f"{stem}_dot_side_decoding_report.html",
    }


def validate_cli_args(args: argparse.Namespace) -> None:
    """Fail early for invalid command-line combinations and numeric settings."""
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if args.decode_sfreq <= 0:
        raise SystemExit("--decode-sfreq must be positive.")
    if args.decode_tmax <= args.decode_tmin:
        raise SystemExit("--decode-tmax must be greater than --decode-tmin.")
    if args.decode_window_samples < 1:
        raise SystemExit("--decode-window-samples must be at least 1.")
    try:
        normalize_trial_averaging_request(args.trials_per_average)
    except argparse.ArgumentTypeError as exc:
        raise SystemExit(str(exc)) from exc
    if args.n_splits < 2:
        raise SystemExit("--n-splits must be at least 2.")
    try:
        sensor_group = normalize_sensor_group(args.sensor_group)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if sensor_group == "posterior" and args.decoder_sensor_type != "grad":
        raise SystemExit("The historical posterior decoder supports only --decoder-sensor-type grad.")
    try:
        resolve_model_condition(args.condition, args.model_source)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc


def validate_stage04_epoching_contract(
    *,
    source_manifest: Path,
    epochs_fif: Path,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    """
    Validate that consumed epochs are recorded as completed Stage-04 outputs.

    The manifest stores cluster-absolute paths in some deployments, while the
    same derivatives can be mounted locally through ``/Volumes/MORWUR``. For
    that reason the check requires status and basename agreement, and records
    the manifest path for auditability instead of requiring identical prefixes.
    """
    manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    stages = (
        manifest.get("derivative_lineage", {})
        .get("stages", {})
    )
    epoching_stage = stages.get("epoching")
    if not epoching_stage:
        raise RuntimeError(
            f"Source manifest lacks derivative_lineage.stages.epoching: {source_manifest}"
        )
    if epoching_stage.get("status") != "completed":
        raise RuntimeError(
            "Stage-04 epoching is not completed in source_manifest: "
            f"{epoching_stage.get('status')!r}"
        )

    qc = epoching_stage.get("qc", {})
    manifest_epochs = Path(str(qc.get("epochs_fif", "")))
    if manifest_epochs.name != epochs_fif.name:
        raise RuntimeError(
            "Stage-04 manifest epochs basename does not match consumed epochs: "
            f"manifest={manifest_epochs.name!r}, consumed={epochs_fif.name!r}"
        )
    configured_event_map = dict(cfg["epoching"]["event_id_map"])
    manifest_event_map = dict(qc.get("event_id_map", {}))
    if manifest_event_map != configured_event_map:
        raise RuntimeError(
            "Stage-04 manifest event_id_map differs from resolved config: "
            f"manifest={manifest_event_map}, config={configured_event_map}"
        )
    if qc.get("baseline_sec", None) is not None:
        raise RuntimeError(
            "Stage-04 epoching manifest reports baseline correction, but the "
            "downstream contract requires baseline=None."
        )

    return {
        "source_manifest": str(source_manifest),
        "stage04_status": str(epoching_stage.get("status")),
        "manifest_epochs_fif": str(manifest_epochs),
        "consumed_epochs_fif": str(epochs_fif),
        "event_id_map": manifest_event_map,
        "baseline_sec": qc.get("baseline_sec", None),
        "n_final_epochs": qc.get("n_final_epochs"),
    }


def validate_decoding_input_contract(
    *,
    epochs: mne.Epochs,
    model_mat: Path,
    condition: str,
    model_source: str,
    cfg: dict[str, Any],
    n_times_target: int,
) -> dict[str, Any]:
    """
    Validate condition, metadata, and model identity coverage for decoding.

    This can run on non-preloaded epochs during ``--validate-only`` because it
    needs metadata and model labels, not full MEG sample arrays.
    """
    model_condition = resolve_model_condition(condition, model_source)
    configured_conditions = list(cfg["epoching"]["event_id_map"])
    if condition not in configured_conditions:
        raise RuntimeError(
            f"Condition {condition!r} not configured for this subject: "
            f"{configured_conditions}"
        )
    if condition not in epochs.event_id:
        raise RuntimeError(f"Condition {condition!r} not in epochs.event_id={epochs.event_id}.")

    epochs_condition = epochs[condition]
    if epochs_condition.metadata is None or "trial_sequence_identity" not in epochs_condition.metadata:
        raise RuntimeError("Epoch metadata missing required column 'trial_sequence_identity'.")

    identity_values = list(
        range(
            int(cfg["subject_profile"]["identity_range"][0]),
            int(cfg["subject_profile"]["identity_range"][1]) + 1,
        )
    )
    trial_identities = epochs_condition.metadata["trial_sequence_identity"].to_numpy(dtype=int)
    labels, _x_by_identity, center_x = build_dot_side_labels(
        model_mat=model_mat,
        condition=model_condition,
        trial_identities=trial_identities,
        identity_values=identity_values,
        n_times_target=n_times_target,
    )
    return {
        "condition": condition,
        "model_source": model_source,
        "model_condition": model_condition,
        "analysis_key": analysis_key(condition, model_source),
        "n_condition_epochs": int(len(epochs_condition)),
        "identity_range": [int(identity_values[0]), int(identity_values[-1])],
        "n_profile_identities": int(len(identity_values)),
        "observed_trial_identities": [
            int(value) for value in sorted(np.unique(trial_identities).tolist())
        ],
        "n_label_times_validated": int(labels.shape[1]),
        "center_x": float(center_x),
    }


# -----------------------------------------------------------------------------
# Geometry and channel selection
# -----------------------------------------------------------------------------

def _meg_head_coords(info: mne.Info, *, meg_kind: str) -> list[dict[str, float | int | str]]:
    """
    Return per-channel head-coordinate metadata for MEG channels.

    Data flow
    - Read MEG sensor positions from channel `loc` vectors.
    - Convert device coordinates to head coordinates with `dev_head_t`.
    - Keep only selected MEG kind (`mag` or `grad`).
    """
    if meg_kind not in {"mag", "grad"}:
        raise ValueError(f"Unsupported meg_kind={meg_kind!r}; expected 'mag' or 'grad'.")

    picks = mne.pick_types(info, meg=meg_kind, eeg=False, eog=False, stim=False, exclude=[])
    pick_set = set(int(i) for i in picks)
    dev_head_t = info.get("dev_head_t")

    rows: list[dict[str, float | int | str]] = []
    for idx, ch in enumerate(info["chs"]):
        if idx not in pick_set:
            continue
        xyz = np.array(ch["loc"][:3], dtype=float)
        if ch["coord_frame"] == FIFF.FIFFV_COORD_DEVICE and dev_head_t is not None:
            xyz = apply_trans(dev_head_t, xyz)
        rows.append(
            {
                "index": idx,
                "name": ch["ch_name"],
                "x": float(xyz[0]),
                "y": float(xyz[1]),
                "z": float(xyz[2]),
            }
        )
    return rows


def select_posterior_decoder_channels(info: mne.Info) -> list[str]:
    """
    Select posterior gradiometers used as decoder features.

    Rule
    - Keep grad channels with head-y below the 45th percentile to focus on
      occipito-parietal sensors while keeping enough channels for stability.
    """
    grad_rows = _meg_head_coords(info, meg_kind="grad")
    y_vals = np.array([row["y"] for row in grad_rows], dtype=float)
    y_threshold = float(np.quantile(y_vals, 0.45))
    selected = [row["name"] for row in grad_rows if row["y"] <= y_threshold]
    if len(selected) < 10:
        raise RuntimeError("Too few decoder channels selected; inspect sensor geometry.")
    return selected


def select_stage07_decoder_channels(
    info: mne.Info,
    *,
    sensor_group: str,
    decoder_sensor_type: str,
) -> tuple[list[str], dict[str, Any]]:
    """
    Select decoder features from a Stage-07 lobe sensor group.

    Data flow
    - Load the canonical Stage-07 mixed MEG channel list for the requested
      group.
    - Keep channels of ``decoder_sensor_type`` present in the current epochs.
    - Preserve MNE channel order so the feature matrix is stable across runs.
    """
    group = normalize_sensor_group(sensor_group)
    if group == "posterior":
        channels = select_posterior_decoder_channels(info)
        return channels, {
            "sensor_group": "posterior",
            "decoder_sensor_type": "grad",
            "selection_policy": (
                "Historical auto-selected posterior gradiometers: grad channels "
                "with head-y below the 45th percentile."
            ),
            "n_selected_channels": int(len(channels)),
            "selected_channels": channels,
        }

    if decoder_sensor_type not in DECODER_SENSOR_TYPES:
        raise ValueError(
            f"Unsupported decoder sensor type {decoder_sensor_type!r}; "
            f"expected one of {DECODER_SENSOR_TYPES}."
        )
    stage07_groups = load_stage07_sensor_groups()
    requested = stage07_groups[group]
    requested_set = set(requested)
    picks = mne.pick_types(
        info,
        meg=decoder_sensor_type,
        eeg=False,
        eog=False,
        stim=False,
        misc=False,
        exclude=[],
    )
    sensor_type_channels = [info["ch_names"][int(idx)] for idx in picks]
    sensor_type_set = set(sensor_type_channels)
    selected = [channel for channel in sensor_type_channels if channel in requested_set]
    dropped_not_sensor_type = [channel for channel in requested if channel not in sensor_type_set]
    if len(selected) < 2:
        raise RuntimeError(
            f"Too few {decoder_sensor_type} channels selected for Stage-07 "
            f"sensor group {group!r}: {selected}."
        )
    return selected, {
        "sensor_group": group,
        "decoder_sensor_type": decoder_sensor_type,
        "sensor_group_definition_file": str(stage07_sensor_group_definition_path()),
        "selection_policy": (
            "Stage-07 MEGIN/Neuromag lobe list filtered to the requested "
            "MNE sensor type and available channels in this epochs file."
        ),
        "n_stage07_group_channels": int(len(requested)),
        "stage07_group_channels": requested,
        "n_available_channels_for_sensor_type": int(len(sensor_type_channels)),
        "n_selected_channels": int(len(selected)),
        "selected_channels": selected,
        "n_dropped_stage07_channels_not_this_sensor_type_or_layout": int(
            len(dropped_not_sensor_type)
        ),
        "dropped_stage07_channels_not_this_sensor_type_or_layout": dropped_not_sensor_type,
    }


def select_left_right_mag_rois(info: mne.Info) -> tuple[list[str], list[str]]:
    """
    Select left/right posterior 3-channel magnetometer ROIs for ERFs.

    Rule
    - Keep only posterior magnetometers (lowest 20% head-y).
    - Split by hemisphere using sign of head-x.
    - Pick the 3 most posterior channels per side.
    """
    mag_rows = _meg_head_coords(info, meg_kind="mag")
    y_vals = np.array([row["y"] for row in mag_rows], dtype=float)
    y_threshold = float(np.quantile(y_vals, 0.20))

    left_candidates = [row for row in mag_rows if row["y"] <= y_threshold and row["x"] < 0]
    right_candidates = [row for row in mag_rows if row["y"] <= y_threshold and row["x"] > 0]
    left_sorted = sorted(left_candidates, key=lambda row: row["y"])
    right_sorted = sorted(right_candidates, key=lambda row: row["y"])

    left_roi = [row["name"] for row in left_sorted[:3]]
    right_roi = [row["name"] for row in right_sorted[:3]]

    if len(left_roi) < 3 or len(right_roi) < 3:
        raise RuntimeError("Could not identify 3-channel left/right posterior mag ROIs.")
    return left_roi, right_roi


def parse_mag_roi_override(raw_value: str | None) -> list[str] | None:
    """
    Parse comma-separated ROI override text into channel names.

    Data flow
    - Split on commas.
    - Strip surrounding whitespace from each token.
    - Drop empty entries and preserve user ordering.
    """
    if raw_value is None:
        return None
    channels = [token.strip() for token in str(raw_value).split(",")]
    channels = [ch for ch in channels if ch]
    return channels or None


def validate_mag_roi_override(info: mne.Info, channels: list[str], side_label: str) -> list[str]:
    """
    Validate a user-supplied ERF ROI channel list.

    Inputs
    - info: epochs info structure that defines available channel names/types.
    - channels: requested ROI channels in desired averaging order.
    - side_label: text label used in error messages ("left" or "right").

    Output
    - The same channel list if all names exist and are magnetometers.
    """
    if not channels:
        raise ValueError(f"{side_label} ROI override is empty after parsing.")

    existing_names = set(info["ch_names"])
    missing = [ch for ch in channels if ch not in existing_names]
    if missing:
        raise ValueError(
            f"{side_label} ROI contains missing channels: {missing}. "
            "Use channel names present in the loaded epochs file."
        )

    mag_picks = mne.pick_types(info, meg="mag", eeg=False, eog=False, stim=False, exclude=[])
    mag_names = {info["ch_names"][idx] for idx in mag_picks}
    non_mag = [ch for ch in channels if ch not in mag_names]
    if non_mag:
        raise ValueError(
            f"{side_label} ROI contains non-magnetometer channels: {non_mag}. "
            "ERF ROI overrides must use MEG magnetometer channels."
        )

    return channels


def resolve_erf_mag_rois(
    info: mne.Info,
    *,
    left_override: list[str] | None,
    right_override: list[str] | None,
) -> tuple[list[str], list[str]]:
    """
    Resolve ERF ROI channels from either auto-selection or explicit overrides.

    Rule
    - If no overrides are provided, use posterior 3-channel auto-selection per
      hemisphere.
    - If one or both overrides are provided, validate each supplied list and
      use it as-is for that hemisphere.
    """
    left_auto, right_auto = select_left_right_mag_rois(info)
    left_roi = left_auto if left_override is None else validate_mag_roi_override(info, left_override, "left")
    right_roi = right_auto if right_override is None else validate_mag_roi_override(info, right_override, "right")
    return left_roi, right_roi


# -----------------------------------------------------------------------------
# Label construction from model trajectories
# -----------------------------------------------------------------------------

def build_dot_side_labels(
    *,
    model_mat: Path,
    condition: str,
    trial_identities: np.ndarray,
    identity_values: list[int],
    n_times_target: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Build per-trial per-time binary labels from model x trajectories.

    Outputs
    - labels: shape (n_trials, n_times), values 0=right, 1=left.
    - x_by_identity: shape (n_identities, n_times), ordered by identity_values.
    - center_x: scalar x threshold defined as `Cfg.rectSize(1)/2`.
    """
    model = sio.loadmat(model_mat, simplify_cells=True)
    xyseqs = model["xySeqs"]
    cfg = model["Cfg"]

    trials = [row for row in xyseqs if str(row.get("condition_label", "")).lower() == condition.lower()]
    by_identity: dict[int, np.ndarray] = {}
    for row in trials:
        identity = int(row["sequence"])
        if identity in by_identity:
            raise RuntimeError(
                f"Duplicate identity {identity} for condition {condition!r} in model MAT."
            )
        by_identity[identity] = np.asarray(row["xy"], dtype=float)

    expected = set(int(value) for value in identity_values)
    found = set(by_identity.keys())
    if found != expected:
        raise RuntimeError(f"Model identity set mismatch: found={sorted(found)}, expected={sorted(expected)}")

    x_by_identity_full = np.stack(
        [by_identity[idx][:, 0] for idx in identity_values],
        axis=0,
    )
    n_times = min(int(n_times_target), int(x_by_identity_full.shape[1]))
    x_by_identity = x_by_identity_full[:, :n_times]

    rect_size = np.asarray(cfg["rectSize"], dtype=float).reshape(-1)
    if rect_size.size == 0:
        raise RuntimeError("Cfg.rectSize missing or empty in model MAT.")
    center_x = float(rect_size[0] / 2.0)

    labels_by_identity = (x_by_identity < center_x).astype(int)
    identity_to_axis = {
        identity: axis for axis, identity in enumerate(identity_values)
    }
    unexpected = sorted(set(int(value) for value in trial_identities).difference(expected))
    if unexpected:
        raise RuntimeError(
            "trial_sequence_identity outside configured identity values: "
            f"{unexpected}"
        )

    labels = labels_by_identity[
        np.asarray([identity_to_axis[int(value)] for value in trial_identities]),
        :,
    ]
    return labels, x_by_identity, center_x


# -----------------------------------------------------------------------------
# Decoding helpers
# -----------------------------------------------------------------------------

def build_centered_window_metadata(
    *,
    times: np.ndarray,
    window_samples: int,
    sfreq: float,
) -> list[dict[str, Any]]:
    """
    Build clipped centered sample-window metadata for every decoded time.

    Policy
    - The requested width is interpreted on the resampled decoding grid.
    - Odd widths are symmetric around the center sample.
    - Even widths include the extra sample after the center sample.
    - At epoch boundaries, windows are clipped to available samples rather than
      dropping centers. The exact effective width is recorded per row.
    """
    if window_samples < 1:
        raise ValueError("window_samples must be at least 1.")
    if times.ndim != 1 or times.size == 0:
        raise ValueError("times must be a non-empty 1D array.")

    n_times = int(times.size)
    samples_before = int((window_samples - 1) // 2)
    samples_after = int(window_samples - 1 - samples_before)
    rows: list[dict[str, Any]] = []
    for center_idx in range(n_times):
        start_idx = max(0, center_idx - samples_before)
        end_idx = min(n_times, center_idx + samples_after + 1)
        n_window = int(end_idx - start_idx)
        rows.append(
            {
                "center_sample": int(center_idx),
                "center_time_sec": float(times[center_idx]),
                "window_start_sample_inclusive": int(start_idx),
                "window_end_sample_exclusive": int(end_idx),
                "window_start_time_sec": float(times[start_idx]),
                "window_end_time_sec": float(times[end_idx - 1]),
                "requested_window_samples": int(window_samples),
                "effective_window_samples": n_window,
                "effective_window_ms": float(1000.0 * n_window / sfreq),
                "edge_policy": "clip_to_available_samples",
                "is_edge_clipped": bool(n_window != window_samples),
                "samples_before_center_requested": samples_before,
                "samples_after_center_requested": samples_after,
                "skip_reason": "",
            }
        )
    return rows


def average_trials_within_identity(
    data: np.ndarray,
    labels: np.ndarray,
    *,
    identities: np.ndarray,
    trials_per_average: int | str,
    random_seed: int,
    source_epoch_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """
    Average same-identity source-trial groups before classification.

    Data flow
    - Split condition epochs by ``trial_sequence_identity``.
    - Fixed-count mode deterministically shuffles source trials inside each
      identity, forms non-overlapping complete groups of
      ``trials_per_average``, and drops incomplete remainders.
    - ``all_trials_per_identity`` mode forms one group containing all retained
      source trials for each identity.
    - Average MEG data within each group; inherit the identity-specific dynamic
      dot-side label time series from the source trials.

    Returns
    -------
    averaged_data, averaged_labels, averaged_identities, metadata
        Pseudo-trial arrays plus machine-readable source-index and drop
        metadata. In fixed-count mode, remainder groups are dropped to keep
        pseudo-trial SNR comparable. An identity with fewer source trials than
        the requested group size is recorded as fully dropped; other identities
        continue. In all-trials mode, every identity with at least one retained
        source trial contributes exactly one pseudo-trial.
    """
    n_trials, _n_channels, n_times = data.shape
    if labels.shape != (n_trials, n_times):
        raise ValueError(f"labels shape mismatch: got {labels.shape}, expected {(n_trials, n_times)}")
    identities = np.asarray(identities, dtype=int)
    if identities.shape != (n_trials,):
        raise ValueError(f"identities shape mismatch: got {identities.shape}, expected {(n_trials,)}")
    try:
        averaging_mode, fixed_group_size = normalize_trial_averaging_request(
            trials_per_average
        )
    except argparse.ArgumentTypeError as exc:
        raise ValueError(str(exc)) from exc
    if source_epoch_indices is None:
        source_epoch_indices = np.arange(n_trials, dtype=int)
    source_epoch_indices = np.asarray(source_epoch_indices, dtype=int)
    if source_epoch_indices.shape != (n_trials,):
        raise ValueError(
            "source_epoch_indices shape mismatch: "
            f"got {source_epoch_indices.shape}, expected {(n_trials,)}"
        )

    rng = np.random.default_rng(random_seed)
    pseudo_data: list[np.ndarray] = []
    pseudo_labels: list[np.ndarray] = []
    pseudo_identities: list[int] = []
    pseudo_trials: list[dict[str, Any]] = []
    dropped_by_identity: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    pseudo_counts: dict[str, int] = {}
    dropped_identities: list[int] = []

    for identity in sorted(int(value) for value in np.unique(identities)):
        identity_rows = np.flatnonzero(identities == identity)
        source_counts[str(identity)] = int(identity_rows.size)
        if averaging_mode == ALL_TRIALS_AVERAGING_MODE:
            grouped_rows = [np.array(identity_rows, dtype=int)]
            dropped_by_identity[str(identity)] = 0
            pseudo_counts[str(identity)] = 1
        else:
            assert fixed_group_size is not None
            shuffled_rows = np.array(identity_rows, dtype=int)
            rng.shuffle(shuffled_rows)

            n_complete = int(shuffled_rows.size // fixed_group_size)
            remainder = int(shuffled_rows.size % fixed_group_size)
            dropped_by_identity[str(identity)] = remainder
            pseudo_counts[str(identity)] = n_complete
            if n_complete == 0:
                dropped_identities.append(identity)
                continue
            grouped_rows = [
                shuffled_rows[
                    group_idx * fixed_group_size : (group_idx + 1) * fixed_group_size
                ]
                for group_idx in range(n_complete)
            ]

        for group_idx, group_rows in enumerate(grouped_rows):
            group_labels = labels[group_rows]
            if not np.all(group_labels == group_labels[0]):
                raise RuntimeError(
                    f"Identity {identity} has inconsistent model-derived labels "
                    "across source trials."
                )
            pseudo_data.append(data[group_rows].mean(axis=0))
            pseudo_labels.append(group_labels[0].copy())
            pseudo_identities.append(identity)
            pseudo_trials.append(
                {
                    "pseudo_trial_index": int(len(pseudo_trials)),
                    "trial_sequence_identity": int(identity),
                    "identity_group_index": int(group_idx),
                    "source_row_indices_zero_based": [
                        int(value) for value in group_rows.tolist()
                    ],
                    "source_epoch_indices_zero_based": [
                        int(value) for value in source_epoch_indices[group_rows].tolist()
                    ],
                    "source_epoch_count": int(group_rows.size),
                }
            )

    if not pseudo_data:
        raise RuntimeError(
            "No complete same-identity pseudo-trial groups remain for "
            f"--trials-per-average={format_trial_averaging_request(trials_per_average)}."
        )

    metadata = {
        "mode": averaging_mode,
        "requested_value": format_trial_averaging_request(trials_per_average),
        "trials_per_average": fixed_group_size,
        "random_seed": int(random_seed),
        "drop_incomplete_remainders": averaging_mode == FIXED_COUNT_AVERAGING_MODE,
        "grouping_policy": (
            "One pseudo-trial per identity, averaging all retained source "
            "trials with that trial_sequence_identity."
            if averaging_mode == ALL_TRIALS_AVERAGING_MODE
            else (
                "Non-overlapping shuffled same-identity groups of "
                f"{fixed_group_size} source trials; incomplete remainders are dropped."
            )
        ),
        "n_source_trials": int(n_trials),
        "n_pseudo_trials": int(len(pseudo_trials)),
        "n_source_trials_by_identity": source_counts,
        "n_pseudo_trials_by_identity": pseudo_counts,
        "dropped_remainders_by_identity": dropped_by_identity,
        "n_dropped_remainder_trials_total": int(sum(dropped_by_identity.values())),
        "dropped_identities_no_complete_group": dropped_identities,
        "n_dropped_identities_no_complete_group": int(len(dropped_identities)),
        "retained_identity_values": sorted(set(pseudo_identities)),
        "pseudo_trials": pseudo_trials,
    }
    return (
        np.stack(pseudo_data, axis=0),
        np.stack(pseudo_labels, axis=0),
        np.asarray(pseudo_identities, dtype=int),
        metadata,
    )


def decode_sliding_windows(
    data: np.ndarray,
    labels: np.ndarray,
    *,
    times: np.ndarray,
    n_splits: int,
    window_samples: int,
    random_seed: int,
    sfreq: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """
    Decode dynamic dot-side labels with clipped centered sliding windows.

    Outputs
    - acc_unbalanced: mean accuracy using original fold class proportions.
    - acc_balanced: mean accuracy after random undersampling of majority class
      in each train/test fold subset.
    - window_rows: per-center metadata including effective window bounds,
      class counts, and explicit NaN skip reasons.
    """
    n_trials, _n_channels, n_times = data.shape
    if labels.shape != (n_trials, n_times):
        raise ValueError(f"labels shape mismatch: got {labels.shape}, expected {(n_trials, n_times)}")
    if times.shape != (n_times,):
        raise ValueError(f"times shape mismatch: got {times.shape}, expected {(n_times,)}")

    window_rows = build_centered_window_metadata(
        times=times,
        window_samples=window_samples,
        sfreq=sfreq,
    )
    acc_unbalanced = np.full(n_times, np.nan, dtype=float)
    acc_balanced = np.full(n_times, np.nan, dtype=float)

    for row in window_rows:
        t_idx = int(row["center_sample"])
        start_idx = int(row["window_start_sample_inclusive"])
        end_idx = int(row["window_end_sample_exclusive"])
        y_t = labels[:, t_idx]
        class_counts = {str(int(k)): int(v) for k, v in Counter(y_t.tolist()).items()}
        row["class_counts"] = class_counts
        if np.unique(y_t).size < 2:
            row["skip_reason"] = "single_class_labels"
            continue
        min_class_count = min(class_counts.values())
        if min_class_count < n_splits:
            row["skip_reason"] = (
                f"too_few_class_samples_for_{n_splits}_fold_stratified_cv"
            )
            continue

        x_t = data[:, :, start_idx:end_idx].reshape(n_trials, -1)
        split_indices = stratified_kfold_indices(
            y_t,
            n_splits=n_splits,
            random_state=random_seed,
        )
        fold_acc_unbalanced: list[float] = []
        fold_acc_balanced: list[float] = []

        for fold_idx, (train_idx, test_idx) in enumerate(split_indices):
            x_train = x_t[train_idx]
            y_train = y_t[train_idx]
            x_test = x_t[test_idx]
            y_test = y_t[test_idx]

            # Strategy 1: original class proportions.
            y_pred_unbalanced = lda_predict_regularized(x_train, y_train, x_test)
            fold_acc_unbalanced.append(float(accuracy(y_test, y_pred_unbalanced)))

            # Strategy 2: random undersampling to parity in train and test.
            rng = np.random.default_rng(
                (random_seed * 1000003) + 100000 + (t_idx * 997) + fold_idx
            )
            train_bal_idx = undersample_majority_to_minority_indices(y_train, rng=rng)
            test_bal_idx = undersample_majority_to_minority_indices(y_test, rng=rng)
            if train_bal_idx.size == 0 or test_bal_idx.size == 0:
                continue
            y_pred_balanced = lda_predict_regularized(
                x_train[train_bal_idx],
                y_train[train_bal_idx],
                x_test[test_bal_idx],
            )
            fold_acc_balanced.append(float(accuracy(y_test[test_bal_idx], y_pred_balanced)))

        if fold_acc_unbalanced:
            acc_unbalanced[t_idx] = float(np.mean(fold_acc_unbalanced))
        if fold_acc_balanced:
            acc_balanced[t_idx] = float(np.mean(fold_acc_balanced))
        if not fold_acc_unbalanced or not fold_acc_balanced:
            row["skip_reason"] = "no_valid_cv_folds_after_balancing"

    return acc_unbalanced, acc_balanced, window_rows


def summarize_training_testing_diagnostics(
    *,
    window_rows: list[dict[str, Any]],
    averaging_metadata: Mapping[str, Any],
    n_decoder_trials: int,
    n_decoder_identities: int,
    n_splits: int,
    window_samples: int,
    sfreq: float,
    accuracy_unbalanced: np.ndarray,
    accuracy_balanced: np.ndarray,
) -> dict[str, Any]:
    """
    Summarize the classification design and per-window CV viability.

    The returned dictionary is stored in the summary JSON and rendered in both
    participant and group reports. It is intentionally derived from the same
    ``window_rows`` consumed by the scores CSV so report text and machine
    diagnostics cannot drift apart.
    """
    n_windows = int(len(window_rows))
    skip_reasons = Counter(row.get("skip_reason", "") or "decoded" for row in window_rows)
    minority_counts: list[int] = []
    decoded_minority_counts: list[int] = []
    abs_imbalances: list[int] = []
    ratios: list[float] = []

    # Each row stores the binary class counts at that center sample. A missing
    # class has count zero, which is what determines whether stratified folds
    # are well-defined.
    for row in window_rows:
        counts_raw = row.get("class_counts", {})
        counts = {int(key): int(value) for key, value in dict(counts_raw).items()}
        count_0 = counts.get(0, 0)
        count_1 = counts.get(1, 0)
        minority = min(count_0, count_1)
        majority = max(count_0, count_1)
        minority_counts.append(int(minority))
        abs_imbalances.append(int(majority - minority))
        if minority > 0:
            ratios.append(float(majority / minority))
        if not row.get("skip_reason"):
            decoded_minority_counts.append(int(minority))

    n_decoded = int(skip_reasons.get("decoded", 0))
    valid_fraction = (float(n_decoded) / float(n_windows)) if n_windows else 0.0
    min_minority_all = min(minority_counts) if minority_counts else None
    min_minority_decoded = (
        min(decoded_minority_counts) if decoded_minority_counts else None
    )
    max_imbalance = max(abs_imbalances) if abs_imbalances else None
    max_ratio = max(ratios) if ratios else None

    if n_windows == 0 or n_decoded == 0:
        assessment_status = "not_well_defined"
        assessment_text = (
            "No time windows produced valid stratified binary decoding. Inspect "
            "skip reasons before interpreting accuracy values."
        )
    elif n_decoded == n_windows:
        assessment_status = "well_defined"
        assessment_text = (
            "Every decoded center had both classes represented with enough "
            "samples for the configured stratified folds."
        )
    else:
        assessment_status = "partially_defined"
        assessment_text = (
            "Some decoded centers lacked two classes or enough minority-class "
            "samples for the configured stratified folds; interpret continuous "
            "accuracy traces with those NaN windows in mind."
        )

    return {
        "classifier": {
            "name": DECODER_MODEL_NAME,
            "details": DECODER_MODEL_DETAILS,
            "classes": {"0": "right_or_center", "1": "left"},
        },
        "cross_validation": {
            "scheme": "stratified_k_fold",
            "n_splits": int(n_splits),
            "split_unit": "same-identity pseudo-trials",
            "skip_policy": (
                "A center sample is skipped when labels are single-class or "
                "the minority class has fewer samples than n_splits."
            ),
        },
        "feature_window": {
            "window_samples": int(window_samples),
            "nominal_window_ms": float(window_samples) / float(sfreq) * 1000.0,
            "edge_windows_are_clipped": True,
        },
        "evaluation": {
            "unbalanced_metric": "accuracy_unbalanced_train",
            "balanced_metric": "accuracy_balanced_undersampled_sets",
            "balancing_strategy": (
                "Random undersampling of majority class to minority count "
                "inside each CV train and test fold; no oversampling."
            ),
            "n_finite_unbalanced_windows": int(np.isfinite(accuracy_unbalanced).sum()),
            "n_finite_balanced_windows": int(np.isfinite(accuracy_balanced).sum()),
        },
        "sample_structure": {
            "averaging_mode": str(averaging_metadata.get("mode", FIXED_COUNT_AVERAGING_MODE)),
            "requested_value": str(
                averaging_metadata.get(
                    "requested_value",
                    format_trial_averaging_request(
                        averaging_metadata.get("trials_per_average", 1)
                    ),
                )
            ),
            "n_source_trials": int(averaging_metadata.get("n_source_trials", 0)),
            "n_decoder_trials": int(n_decoder_trials),
            "n_decoder_identities": int(n_decoder_identities),
            "n_dropped_identities_no_complete_group": int(
                averaging_metadata.get("n_dropped_identities_no_complete_group", 0)
            ),
            "n_dropped_remainder_trials_total": int(
                averaging_metadata.get("n_dropped_remainder_trials_total", 0)
            ),
        },
        "class_balance": {
            "n_windows": n_windows,
            "n_decoded_windows": n_decoded,
            "n_skipped_windows": int(n_windows - n_decoded),
            "decoded_window_fraction": valid_fraction,
            "skip_reasons": {str(key): int(value) for key, value in skip_reasons.items()},
            "min_minority_class_count_all_windows": (
                int(min_minority_all) if min_minority_all is not None else None
            ),
            "min_minority_class_count_decoded_windows": (
                int(min_minority_decoded) if min_minority_decoded is not None else None
            ),
            "max_abs_class_count_difference": (
                int(max_imbalance) if max_imbalance is not None else None
            ),
            "max_majority_to_minority_ratio": (
                float(max_ratio) if max_ratio is not None else None
            ),
        },
        "assessment": {
            "status": assessment_status,
            "text": assessment_text,
        },
    }


def _format_optional_number(value: Any, *, digits: int = 3) -> str:
    """Format optional numeric values for report tables."""
    if value is None:
        return "n/a"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def training_testing_diagnostics_html(summary: Mapping[str, Any]) -> str:
    """
    Render classifier, CV, class-balance, and assessment details as HTML.

    The report is intentionally concise but contains enough information to
    judge whether the decoded time points were well-defined under the selected
    averaging mode and number of stratified folds.
    """
    diagnostics = summary.get("training_testing_diagnostics")
    if not isinstance(diagnostics, Mapping):
        return (
            "<p>Training/testing diagnostics are unavailable for this summary. "
            "Rebuild the decoding output with the current script to populate "
            "this section.</p>"
        )

    classifier = diagnostics.get("classifier", {})
    cv = diagnostics.get("cross_validation", {})
    feature_window = diagnostics.get("feature_window", {})
    evaluation = diagnostics.get("evaluation", {})
    sample_structure = diagnostics.get("sample_structure", {})
    class_balance = diagnostics.get("class_balance", {})
    assessment = diagnostics.get("assessment", {})
    skip_reasons = class_balance.get("skip_reasons", {})
    skip_text = ", ".join(
        f"{html.escape(str(key))}: {int(value)}"
        for key, value in sorted(dict(skip_reasons).items())
    ) or "none"
    decoded_fraction = class_balance.get("decoded_window_fraction")
    decoded_percent = (
        "n/a"
        if decoded_fraction is None
        else f"{float(decoded_fraction) * 100.0:.1f}%"
    )

    rows = [
        ("Classifier", classifier.get("name", "n/a")),
        ("Classifier details", classifier.get("details", "n/a")),
        ("CV scheme", cv.get("scheme", "n/a")),
        ("CV folds", cv.get("n_splits", "n/a")),
        ("Split unit", cv.get("split_unit", "n/a")),
        ("Averaging mode", sample_structure.get("averaging_mode", "n/a")),
        ("Averaging request", sample_structure.get("requested_value", "n/a")),
        ("Source trials", sample_structure.get("n_source_trials", "n/a")),
        ("Decoder pseudo-trials", sample_structure.get("n_decoder_trials", "n/a")),
        ("Decoder identities", sample_structure.get("n_decoder_identities", "n/a")),
        ("Window samples", feature_window.get("window_samples", "n/a")),
        (
            "Nominal window ms",
            _format_optional_number(feature_window.get("nominal_window_ms"), digits=2),
        ),
        ("Unbalanced metric", evaluation.get("unbalanced_metric", "n/a")),
        ("Balanced metric", evaluation.get("balanced_metric", "n/a")),
        ("Balancing strategy", evaluation.get("balancing_strategy", "n/a")),
        (
            "Decoded windows",
            f"{class_balance.get('n_decoded_windows', 'n/a')} / "
            f"{class_balance.get('n_windows', 'n/a')} ({decoded_percent})",
        ),
        ("Skip reasons", skip_text),
        (
            "Minimum minority class count",
            class_balance.get("min_minority_class_count_all_windows", "n/a"),
        ),
        (
            "Minimum decoded-window minority count",
            class_balance.get("min_minority_class_count_decoded_windows", "n/a"),
        ),
        (
            "Worst majority/minority ratio",
            _format_optional_number(
                class_balance.get("max_majority_to_minority_ratio"),
                digits=2,
            ),
        ),
        ("Assessment", assessment.get("status", "n/a")),
    ]
    table_rows = "\n".join(
        "<tr><th>"
        + html.escape(str(label))
        + "</th><td>"
        + html.escape(str(value))
        + "</td></tr>"
        for label, value in rows
    )
    assessment_text = html.escape(str(assessment.get("text", "")))
    return (
        "<p>"
        + assessment_text
        + "</p><table><tbody>"
        + table_rows
        + "</tbody></table>"
    )


def stratified_kfold_indices(
    y: np.ndarray,
    *,
    n_splits: int,
    random_state: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    Build stratified K-fold train/test indices with NumPy only.

    Data flow
    - Shuffle each class index list.
    - Split class indices into `n_splits` chunks.
    - Combine class chunks per fold to preserve class ratios in test sets.
    """
    y = np.asarray(y, dtype=int)
    unique = np.unique(y)
    if unique.size != 2:
        raise ValueError("Stratified split expects binary labels.")

    rng = np.random.default_rng(random_state)
    per_class_folds: dict[int, list[np.ndarray]] = {}
    for cls in unique:
        cls_idx = np.flatnonzero(y == cls)
        rng.shuffle(cls_idx)
        per_class_folds[int(cls)] = list(np.array_split(cls_idx, n_splits))

    all_idx = np.arange(y.size, dtype=int)
    splits: list[tuple[np.ndarray, np.ndarray]] = []
    for fold in range(n_splits):
        test_idx = np.concatenate([per_class_folds[int(cls)][fold] for cls in unique])
        test_idx = np.sort(test_idx)
        train_mask = np.ones(y.size, dtype=bool)
        train_mask[test_idx] = False
        train_idx = all_idx[train_mask]
        splits.append((train_idx, test_idx))
    return splits


def accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute standard classification accuracy."""
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    if y_true.size == 0:
        return np.nan
    return float(np.mean(y_true == y_pred))


def undersample_majority_to_minority_indices(
    y: np.ndarray,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Balance binary labels by random undersampling.

    Strategy details
    - Find minority class count.
    - Randomly sample that many rows from both classes without replacement.
    - Shuffle selected indices for downstream model fitting.
    """
    y = np.asarray(y, dtype=int)
    classes = np.unique(y)
    if classes.size != 2:
        return np.array([], dtype=int)
    idx_a = np.flatnonzero(y == classes[0])
    idx_b = np.flatnonzero(y == classes[1])
    if idx_a.size == 0 or idx_b.size == 0:
        return np.array([], dtype=int)
    n_keep = min(idx_a.size, idx_b.size)
    sel_a = rng.choice(idx_a, size=n_keep, replace=False)
    sel_b = rng.choice(idx_b, size=n_keep, replace=False)
    balanced_idx = np.concatenate([sel_a, sel_b])
    rng.shuffle(balanced_idx)
    return balanced_idx


def lda_predict_regularized(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    *,
    reg: float = 1e-2,
) -> np.ndarray:
    """
    Predict binary labels with regularized LDA.

    Non-obvious implementation detail
    - Features are z-scored with train-only statistics.
    - Pooled covariance gets ridge regularization before inversion.
    """
    y_train = np.asarray(y_train, dtype=int)
    x_train = np.asarray(x_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)

    mu = x_train.mean(axis=0, keepdims=True)
    sigma = x_train.std(axis=0, keepdims=True)
    sigma[sigma == 0] = 1.0
    x_tr = (x_train - mu) / sigma
    x_te = (x_test - mu) / sigma

    x0 = x_tr[y_train == 0]
    x1 = x_tr[y_train == 1]
    if x0.size == 0 or x1.size == 0:
        raise ValueError("Regularized LDA requires both classes in the training set.")
    m0 = x0.mean(axis=0)
    m1 = x1.mean(axis=0)

    centered = np.vstack([x0 - m0, x1 - m1])
    denom = max(x0.shape[0] + x1.shape[0] - 2, 1)
    scaled_centered = centered / np.sqrt(float(denom))
    mean_delta = m1 - m0

    if mean_delta.size <= scaled_centered.shape[0]:
        pooled = scaled_centered.T @ scaled_centered
        pooled = pooled + reg * np.eye(pooled.shape[0], dtype=float)
        w = np.linalg.solve(pooled, mean_delta)
    else:
        # Woodbury form avoids repeatedly solving huge feature-space systems
        # after channel-by-window flattening.
        gram = scaled_centered @ scaled_centered.T
        middle = np.eye(gram.shape[0], dtype=float) + (gram / reg)
        rhs = (scaled_centered @ mean_delta) / reg
        correction = scaled_centered.T @ np.linalg.solve(middle, rhs)
        w = (mean_delta / reg) - (correction / reg)

    b = -0.5 * float((m1 + m0) @ w)
    score = x_te @ w + b
    return (score > 0.0).astype(int)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def build_three_panel_plot(
    *,
    times: np.ndarray,
    contra_erf_ft: np.ndarray,
    ipsi_erf_ft: np.ndarray,
    ipsi_minus_contra_ft: np.ndarray,
    accuracy_unbalanced: np.ndarray,
    accuracy_balanced: np.ndarray,
    left_roi: list[str],
    right_roi: list[str],
    n_decoder_channels: int,
    decoder_sensor_label: str,
    erf_baseline_window: tuple[float, float],
    condition_label: str,
    model_source: str,
    n_splits: int,
    occlusion_timing: dict[str, Any] | None = None,
) -> plt.Figure:
    """
    Create the participant ERF/lateralization/decoding figure.

    Occluded conditions receive labeled decoding-panel boundaries at the first
    completely occluded stimulus frame and the first reappearance frame. The
    timing mapping is expected to use seconds relative to condition onset.
    """
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
    ax_top, ax_mid, ax_bottom = axes

    # Top panel: baseline-corrected contra/ipsi ERFs.
    roi_count_label = f"L={len(left_roi)} mag, R={len(right_roi)} mag"
    ax_top.plot(times, contra_erf_ft, color="#1f77b4", lw=2, label=f"Contra posterior ROI ({roi_count_label})")
    ax_top.plot(times, ipsi_erf_ft, color="#d62728", lw=2, label=f"Ipsi posterior ROI ({roi_count_label})")
    ax_top.axvline(0.0, color="k", lw=1, ls="--", alpha=0.7)
    ax_top.set_ylabel("Amplitude (fT)")
    ax_top.set_title(
        f"{condition_label} ERFs (contra vs ipsi, posterior ROI means, "
        f"baseline {erf_baseline_window[0]:.1f}..{erf_baseline_window[1]:.1f} s)"
    )
    ax_top.legend(loc="upper right", frameon=True)

    # Middle panel: signed lateralization summary.
    ax_mid.plot(times, ipsi_minus_contra_ft, color="#9467bd", lw=2, label="Ipsi - Contra")
    ax_mid.axhline(0.0, color="k", lw=1, ls=":", alpha=0.8)
    ax_mid.axvline(0.0, color="k", lw=1, ls="--", alpha=0.7)
    ax_mid.set_ylabel("Amplitude (fT)")
    ax_mid.set_title("Ipsi - Contra")
    ax_mid.legend(loc="upper right", frameon=True)

    # Bottom panel: decoding accuracy comparison.
    ax_bottom.plot(
        times,
        accuracy_unbalanced,
        color="#2ca02c",
        lw=2,
        label=f"Accuracy (unbalanced train, {n_splits}-fold CV)",
    )
    ax_bottom.plot(
        times,
        accuracy_balanced,
        color="#ff7f0e",
        lw=2,
        label=f"Accuracy (balanced via random undersampling, {n_splits}-fold CV)",
    )
    ax_bottom.axhline(0.5, color="k", lw=1, ls=":", alpha=0.8, label="Chance = 0.5")
    ax_bottom.axvline(0.0, color="k", lw=1, ls="--", alpha=0.7)
    add_occlusion_boundary_lines(
        ax_bottom,
        condition=condition_label,
        occlusion_timing=occlusion_timing,
    )
    ax_bottom.set_xlabel("Time from condition onset (s)")
    ax_bottom.set_ylabel("Accuracy")
    ax_bottom.set_ylim(0.35, 1.0)
    ax_bottom.set_title(
        f"Dot side decoding from {decoder_sensor_label} "
        f"(n={n_decoder_channels}; model={model_source})"
    )
    ax_bottom.legend(loc="upper right", frameon=True)

    info_text = f"Left ROI: {', '.join(left_roi)}\nRight ROI: {', '.join(right_roi)}"
    ax_top.text(
        0.01,
        0.01,
        info_text,
        transform=ax_top.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"},
    )

    return fig


def add_occlusion_boundary_lines(
    ax: plt.Axes,
    *,
    condition: str,
    occlusion_timing: dict[str, Any] | None,
) -> None:
    """
    Add configured occlusion boundaries to one decoding axis when applicable.

    Always-visible plots remain unchanged. For an occluded condition, absence
    of timing metadata is an error because silently omitting the scientific
    event markers would make participant and group figures inconsistent.
    """
    if "occluded" not in str(condition).lower():
        return
    if occlusion_timing is None:
        raise ValueError(f"Occlusion timing is required for condition {condition!r}.")

    boundaries = (
        (
            "first_complete_occlusion_sec",
            "First complete occlusion",
            "#7f7f7f",
            "-.",
        ),
        ("reappearance_sec", "First reappearance", "#8c564b", "--"),
    )
    for key, label, color, linestyle in boundaries:
        if key not in occlusion_timing:
            raise ValueError(f"Occlusion timing lacks required field {key!r}.")
        ax.axvline(
            float(occlusion_timing[key]),
            color=color,
            lw=1.4,
            ls=linestyle,
            alpha=0.9,
            label=label,
        )


def build_sensor_layout_topo_figures(
    *,
    epochs_erf: mne.Epochs,
    labels: np.ndarray,
    condition_label: str,
) -> tuple[plt.Figure, plt.Figure]:
    """
    Build whole-sensor ipsi-minus-contra ERF layout plots.

    Data flow
    - Compute per-sensor ipsi-minus-contra time courses from baseline-corrected
      epochs using dynamic dot-side labels at each sample.
    - Build one evoked object for magnetometers and one for planar
      gradiometers.
    - Plot magnetometers with `plot_topo`.
    - Plot planar gradiometers with `merge_grads=True`, which combines each
      orthogonal grad pair using RMS magnitude (as in FLUX-style usage).

    Outputs
    - fig_mag_topo: ipsi-minus-contra ERFs over the full magnetometer layout.
    - fig_grad_topo_rms: ipsi-minus-contra ERFs over planar grad layout with
      RMS pair merge.
    """
    def _ipsi_minus_contra_evoked_for_kind(meg_kind: str) -> mne.EvokedArray:
        """
        Compute ipsi-minus-contra evoked traces for one MEG sensor kind.

        Non-obvious logic
        - "Ipsi" and "contra" are defined per channel hemisphere and per sample:
          left-hemisphere channels are ipsi when the dot is on the left, while
          right-hemisphere channels are ipsi when the dot is on the right.
        """
        rows = _meg_head_coords(epochs_erf.info, meg_kind=meg_kind)
        if not rows:
            raise RuntimeError(f"No {meg_kind} channels found for ipsi-contra topography.")

        channel_names = [str(row["name"]) for row in rows]
        x_coords = np.array([float(row["x"]) for row in rows], dtype=float)

        # Drop potential midline channels to avoid ambiguous ipsi/contra mapping.
        keep_mask = ~np.isclose(x_coords, 0.0, atol=1e-9)
        if not np.any(keep_mask):
            raise RuntimeError(f"No lateralized {meg_kind} channels remain after excluding x==0.")
        channel_names = [name for name, keep in zip(channel_names, keep_mask, strict=True) if keep]
        x_coords = x_coords[keep_mask]

        data = epochs_erf.get_data(picks=channel_names)  # (n_trials, n_channels, n_times)
        if labels.shape != (data.shape[0], data.shape[2]):
            raise RuntimeError(
                "Label matrix shape does not match ERF data shape for ipsi-contra topography."
            )
        dot_left_mask = labels.astype(bool)[:, np.newaxis, :]  # (n_trials, 1, n_times)
        channel_is_left = (x_coords < 0.0)[np.newaxis, :, np.newaxis]  # (1, n_channels, 1)

        # Per channel/time sample, choose ipsi and contra trial subsets.
        ipsi_mask = np.where(channel_is_left, dot_left_mask, ~dot_left_mask)
        contra_mask = ~ipsi_mask

        ipsi_count = ipsi_mask.sum(axis=0, dtype=float)
        contra_count = contra_mask.sum(axis=0, dtype=float)
        if np.any(ipsi_count == 0.0) or np.any(contra_count == 0.0):
            raise RuntimeError("Encountered empty ipsi/contra bin while building topography evokeds.")

        ipsi_sum = np.where(ipsi_mask, data, 0.0).sum(axis=0)
        contra_sum = np.where(contra_mask, data, 0.0).sum(axis=0)
        ipsi_mean = ipsi_sum / ipsi_count
        contra_mean = contra_sum / contra_count
        ipsi_minus_contra = ipsi_mean - contra_mean

        picks = mne.pick_channels(epochs_erf.ch_names, include=channel_names, ordered=True)
        picked_info = mne.pick_info(epochs_erf.info, picks)
        return mne.EvokedArray(
            ipsi_minus_contra,
            picked_info,
            tmin=float(epochs_erf.times[0]),
            nave=int(len(epochs_erf)),
            comment=f"ipsi_minus_contra_{meg_kind}",
        )

    evoked_mag_diff = _ipsi_minus_contra_evoked_for_kind("mag")
    evoked_grad_diff = _ipsi_minus_contra_evoked_for_kind("grad")

    # Magnetometer layout ERFs across all available mag channels.
    fig_mag_topo = evoked_mag_diff.plot_topo(
        title=f"{condition_label} | Magnetometers (ipsi - contra)",
        show=False,
    )

    # RMS-combined planar gradiometer layout (orthogonal pair merge).
    fig_grad_topo_rms = evoked_grad_diff.plot_topo(
        title=f"{condition_label} | Planar gradiometers (ipsi - contra, RMS pair merge)",
        merge_grads=True,
        show=False,
    )

    return fig_mag_topo, fig_grad_topo_rms


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

def main() -> None:
    """Run the end-to-end sanity-check decoding workflow."""
    args = build_arg_parser().parse_args()
    validate_cli_args(args)
    sensor_group = normalize_sensor_group(args.sensor_group)
    model_condition = resolve_model_condition(args.condition, args.model_source)
    current_analysis_key = analysis_key(args.condition, args.model_source)
    epochs_fif, model_mat, out_dir, source_manifest, cfg = resolve_input_output_paths(args)
    output_paths = sanity_output_paths(out_dir, args.condition, args.model_source)
    left_roi_override = parse_mag_roi_override(args.left_roi_mag_channels)
    right_roi_override = parse_mag_roi_override(args.right_roi_mag_channels)

    stage04_validation = validate_stage04_epoching_contract(
        source_manifest=source_manifest,
        epochs_fif=epochs_fif,
        cfg=cfg,
    )

    if args.resume and all(path.is_file() for path in output_paths.values()):
        print("Sanity-check decoding outputs already exist; resume mode reuses them.")
        for label, path in output_paths.items():
            print(f"{label}: {path}")
        return

    if not args.validate_only and not args.resume and not args.overwrite:
        existing_outputs = [str(path) for path in output_paths.values() if path.exists()]
        if existing_outputs:
            raise RuntimeError(
                "Refusing to overwrite existing sanity-check outputs without "
                f"--overwrite or --resume: {existing_outputs}"
            )

    # -------------------------------------------------------------------------
    # Load epochs and prepare separate decode/ERF branches.
    # -------------------------------------------------------------------------
    if args.validate_only:
        epochs_header = mne.read_epochs(epochs_fif, preload=False, verbose=False)
        validation = validate_decoding_input_contract(
            epochs=epochs_header,
            model_mat=model_mat,
            condition=args.condition,
            model_source=args.model_source,
            cfg=cfg,
            n_times_target=int(round((args.decode_tmax - args.decode_tmin) * args.decode_sfreq)) + 1,
        )
        decoder_channels, decoder_selection = select_stage07_decoder_channels(
            epochs_header[args.condition].info,
            sensor_group=sensor_group,
            decoder_sensor_type=str(args.decoder_sensor_type),
        )
        if args.print_summary:
            print(
                json.dumps(
                    {
                        "stage04_validation": stage04_validation,
                        "decoding_input_validation": validation,
                        "sensor_group": sensor_group,
                        "decoder_sensor_type": str(args.decoder_sensor_type),
                        "decoder_sensor_selection": decoder_selection,
                        "out_dir": str(out_dir),
                        "outputs": {key: str(value) for key, value in output_paths.items()},
                        "analysis_key": current_analysis_key,
                        "model_source": args.model_source,
                        "model_condition": model_condition,
                        "n_decoder_channels": int(len(decoder_channels)),
                        "decode_window_samples": int(args.decode_window_samples),
                        "trials_per_average": format_trial_averaging_request(
                            args.trials_per_average
                        ),
                        "random_seed": int(args.random_seed),
                    },
                    indent=2,
                )
            )
        else:
            print("Dot-side decoding input validation passed.")
        return

    epochs_all = mne.read_epochs(epochs_fif, preload=False, verbose=False)
    validate_decoding_input_contract(
        epochs=epochs_all,
        model_mat=model_mat,
        condition=args.condition,
        model_source=args.model_source,
        cfg=cfg,
        n_times_target=int(round((args.decode_tmax - args.decode_tmin) * args.decode_sfreq)) + 1,
    )

    epochs_condition = epochs_all[args.condition].copy().load_data()
    decoder_channels, decoder_selection = select_stage07_decoder_channels(
        epochs_condition.info,
        sensor_group=sensor_group,
        decoder_sensor_type=str(args.decoder_sensor_type),
    )
    left_roi, right_roi = resolve_erf_mag_rois(
        epochs_condition.info,
        left_override=left_roi_override,
        right_override=right_roi_override,
    )

    erf_baseline_window = (-0.4, -0.2)
    if (
        erf_baseline_window[0] < float(epochs_condition.tmin)
        or erf_baseline_window[1] > float(epochs_condition.tmax)
    ):
        raise RuntimeError(
            f"ERF baseline {erf_baseline_window} lies outside epoch range "
            f"[{epochs_condition.tmin}, {epochs_condition.tmax}]."
        )
    epochs_condition = epochs_condition.crop(
        tmin=min(erf_baseline_window[0], args.decode_tmin),
        tmax=args.decode_tmax,
    )

    epochs_decode = (
        epochs_condition.copy()
        .crop(tmin=args.decode_tmin, tmax=args.decode_tmax)
        .pick(decoder_channels)
        .resample(args.decode_sfreq)
    )
    erf_channel_picks = mne.pick_types(
        epochs_condition.info,
        meg=True,
        eeg=False,
        eog=False,
        stim=False,
        misc=False,
        exclude=[],
    )
    erf_channels = [epochs_condition.ch_names[int(idx)] for idx in erf_channel_picks]
    epochs_erf = (
        epochs_condition.copy()
        .pick(erf_channels)
        .resample(args.decode_sfreq)
        .apply_baseline(erf_baseline_window)
    )
    epochs_erf = epochs_erf.crop(
        tmin=args.decode_tmin,
        tmax=min(args.decode_tmax, float(epochs_erf.tmax)),
    )

    trial_identities = epochs_decode.metadata["trial_sequence_identity"].to_numpy(dtype=int)

    # -------------------------------------------------------------------------
    # Build model-derived dot-side labels and align time axes.
    # -------------------------------------------------------------------------
    labels, _x_by_identity, center_x = build_dot_side_labels(
        model_mat=model_mat,
        condition=model_condition,
        trial_identities=trial_identities,
        identity_values=list(
            range(
                int(cfg["subject_profile"]["identity_range"][0]),
                int(cfg["subject_profile"]["identity_range"][1]) + 1,
            )
        ),
        n_times_target=len(epochs_decode.times),
    )
    n_times_aligned = min(labels.shape[1], len(epochs_decode.times), len(epochs_erf.times))
    if n_times_aligned < 1:
        raise RuntimeError("No overlapping decoding/ERF/model time samples remain.")
    if n_times_aligned < len(epochs_decode.times):
        epochs_decode = epochs_decode.copy().crop(
            tmin=float(epochs_decode.times[0]),
            tmax=float(epochs_decode.times[n_times_aligned - 1]),
        )
    if n_times_aligned < len(epochs_erf.times):
        epochs_erf = epochs_erf.copy().crop(
            tmin=float(epochs_erf.times[0]),
            tmax=float(epochs_erf.times[n_times_aligned - 1]),
        )
    n_times_aligned = min(labels.shape[1], len(epochs_decode.times), len(epochs_erf.times))
    labels = labels[:, :n_times_aligned]
    times = epochs_decode.times[:n_times_aligned].copy()

    # -------------------------------------------------------------------------
    # Build contra/ipsi waveforms from original source trials.
    # -------------------------------------------------------------------------

    source_epoch_indices = np.asarray(epochs_decode.selection, dtype=int)
    decode_data_source = epochs_decode.get_data()[:, :, :n_times_aligned]
    left_data = epochs_erf.get_data(picks=left_roi)[:, :, :n_times_aligned]
    right_data = epochs_erf.get_data(picks=right_roi)[:, :, :n_times_aligned]

    # Dynamic contra/ipsi definition by dot side at each sample:
    # label=1 (dot left): contra=right ROI, ipsi=left ROI
    # label=0 (dot right): contra=left ROI, ipsi=right ROI
    left_roi_mean = left_data.mean(axis=1)
    right_roi_mean = right_data.mean(axis=1)
    dot_left_mask = labels.astype(bool)
    contra_samples = np.where(dot_left_mask, right_roi_mean, left_roi_mean)
    ipsi_samples = np.where(dot_left_mask, left_roi_mean, right_roi_mean)
    contra_erf_ft = contra_samples.mean(axis=0) * 1e15
    ipsi_erf_ft = ipsi_samples.mean(axis=0) * 1e15
    ipsi_minus_contra_ft = ipsi_erf_ft - contra_erf_ft

    # -------------------------------------------------------------------------
    # Build optional same-identity pseudo-trials and decode sliding windows.
    # -------------------------------------------------------------------------
    decode_data, decode_labels, decode_identities, averaging_metadata = average_trials_within_identity(
        decode_data_source,
        labels,
        identities=trial_identities,
        trials_per_average=args.trials_per_average,
        random_seed=args.random_seed,
        source_epoch_indices=source_epoch_indices,
    )
    accuracy_unbalanced, accuracy_balanced, window_rows = decode_sliding_windows(
        decode_data,
        decode_labels,
        times=times,
        n_splits=args.n_splits,
        window_samples=args.decode_window_samples,
        random_seed=args.random_seed,
        sfreq=float(epochs_decode.info["sfreq"]),
    )

    valid_unbalanced = np.isfinite(accuracy_unbalanced)
    valid_balanced = np.isfinite(accuracy_balanced)
    if not np.any(valid_unbalanced) or not np.any(valid_balanced):
        raise RuntimeError("No valid decoding accuracies were produced.")

    peak_idx_unbalanced = int(np.nanargmax(accuracy_unbalanced))
    peak_idx_balanced = int(np.nanargmax(accuracy_balanced))
    peak_time_unbalanced = float(times[peak_idx_unbalanced])
    peak_time_balanced = float(times[peak_idx_balanced])
    peak_accuracy_unbalanced = float(accuracy_unbalanced[peak_idx_unbalanced])
    peak_accuracy_balanced = float(accuracy_balanced[peak_idx_balanced])
    mean_accuracy_unbalanced = float(np.nanmean(accuracy_unbalanced))
    mean_accuracy_balanced = float(np.nanmean(accuracy_balanced))
    occlusion_timing = (
        occlusion_timing_from_config(cfg)
        if "occluded" in args.condition.lower()
        else None
    )

    # -------------------------------------------------------------------------
    # Export figure and tabular metrics.
    # -------------------------------------------------------------------------
    fig = build_three_panel_plot(
        times=times,
        contra_erf_ft=contra_erf_ft,
        ipsi_erf_ft=ipsi_erf_ft,
        ipsi_minus_contra_ft=ipsi_minus_contra_ft,
        accuracy_unbalanced=accuracy_unbalanced,
        accuracy_balanced=accuracy_balanced,
        left_roi=left_roi,
        right_roi=right_roi,
        n_decoder_channels=len(decoder_channels),
        decoder_sensor_label=(
            "posterior gradiometers"
            if sensor_group == "posterior"
            else f"Stage-07 {sensor_group} {args.decoder_sensor_type} channels"
        ),
        erf_baseline_window=erf_baseline_window,
        condition_label=args.condition,
        model_source=args.model_source,
        n_splits=args.n_splits,
        occlusion_timing=occlusion_timing,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_path = output_paths["figure"]
    fig.savefig(fig_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    # Build full-layout ERF topographies for report diagnostics.
    fig_mag_topo, fig_grad_topo_rms = build_sensor_layout_topo_figures(
        epochs_erf=epochs_erf,
        labels=labels,
        condition_label=args.condition,
    )

    scores_csv = output_paths["scores_csv"]
    scores_df = pd.DataFrame(window_rows)
    scores_df.insert(0, "time_sec", times)
    scores_df["accuracy_unbalanced_train"] = accuracy_unbalanced
    scores_df["accuracy_balanced_undersampled_sets"] = accuracy_balanced
    scores_df["left_fraction"] = decode_labels.mean(axis=0)
    scores_df.to_csv(scores_csv, index=False)

    skip_reasons = {
        str(reason): int(count)
        for reason, count in Counter(
            row.get("skip_reason", "") or "decoded" for row in window_rows
        ).items()
    }
    training_testing_diagnostics = summarize_training_testing_diagnostics(
        window_rows=window_rows,
        averaging_metadata=averaging_metadata,
        n_decoder_trials=int(decode_data.shape[0]),
        n_decoder_identities=int(np.unique(decode_identities).size),
        n_splits=int(args.n_splits),
        window_samples=int(args.decode_window_samples),
        sfreq=float(epochs_decode.info["sfreq"]),
        accuracy_unbalanced=accuracy_unbalanced,
        accuracy_balanced=accuracy_balanced,
    )

    summary = {
        "condition": args.condition,
        "model_source": args.model_source,
        "model_condition": model_condition,
        "analysis_key": current_analysis_key,
        "sensor_group": sensor_group,
        "decoder_sensor_type": str(args.decoder_sensor_type),
        "decoder_sensor_selection": decoder_selection,
        "subject": str(cfg["subject"]["id"]),
        "session": str(cfg["layout"]["session_bids_id"]).removeprefix("ses-"),
        "stage04_validation": stage04_validation,
        "event_id_map": dict(cfg["epoching"]["event_id_map"]),
        "identity_range": list(cfg["subject_profile"]["identity_range"]),
        "n_trials": int(decode_data.shape[0]),
        "n_source_trials": int(len(epochs_decode)),
        "n_decoder_trials": int(decode_data.shape[0]),
        "n_decoder_identities": int(np.unique(decode_identities).size),
        "n_times": int(len(times)),
        "sfreq_hz": float(epochs_decode.info["sfreq"]),
        "decode_tmin_sec": float(args.decode_tmin),
        "decode_tmax_sec": float(args.decode_tmax),
        "decode_window_samples": int(args.decode_window_samples),
        "decode_window_policy": (
            "Centered sliding windows on the resampled decoding grid; odd "
            "widths are symmetric, even widths place the extra sample after "
            "the center, and epoch-edge windows are clipped to available "
            "samples."
        ),
        "n_splits": int(args.n_splits),
        "random_seed": int(args.random_seed),
        "n_decoder_channels": int(len(decoder_channels)),
        "decoder_channels": decoder_channels,
        "left_roi_mag_channels": left_roi,
        "right_roi_mag_channels": right_roi,
        "erf_baseline_window_sec": [float(erf_baseline_window[0]), float(erf_baseline_window[1])],
        "occlusion_timing": occlusion_timing,
        "trial_averaging": averaging_metadata,
        "training_testing_diagnostics": training_testing_diagnostics,
        "window_skip_reasons": skip_reasons,
        "n_skipped_windows": int(sum(1 for row in window_rows if row.get("skip_reason"))),
        "n_decoded_windows": int(sum(1 for row in window_rows if not row.get("skip_reason"))),
        "balancing_strategy": (
            "Random undersampling of majority class to minority count inside each "
            "CV train and test fold. No oversampling is used."
        ),
        "center_x": center_x,
        "peak_accuracy_unbalanced_train": peak_accuracy_unbalanced,
        "peak_time_sec_unbalanced_train": peak_time_unbalanced,
        "mean_accuracy_unbalanced_train": mean_accuracy_unbalanced,
        "peak_accuracy_balanced_train": peak_accuracy_balanced,
        "peak_time_sec_balanced_train": peak_time_balanced,
        "mean_accuracy_balanced_train": mean_accuracy_balanced,
        "figure_path": str(fig_path),
        "scores_csv": str(scores_csv),
        "epochs_fif": str(epochs_fif),
        "model_mat": str(model_mat),
    }

    summary_json = output_paths["summary_json"]
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # -------------------------------------------------------------------------
    # Export MNE report with metrics and figure.
    # -------------------------------------------------------------------------
    report = mne.Report(
        title=(
            f"Dot-side decoding sanity check | sub-{args.subject} ses-{args.session} "
            f"| {args.condition} | {args.model_source} | {sensor_group}"
        )
    )

    report.add_html(
        "<pre>"
        + html.escape(
            "Dot-side decoding summary\n"
            f"condition={summary['condition']}\n"
            f"model_source={summary['model_source']}\n"
            f"model_condition={summary['model_condition']}\n"
            f"analysis_key={summary['analysis_key']}\n"
            f"sensor_group={summary['sensor_group']}\n"
            f"decoder_sensor_type={summary['decoder_sensor_type']}\n"
            f"subject=sub-{summary['subject']}\n"
            f"session=ses-{summary['session']}\n"
            f"n_trials={summary['n_trials']}\n"
            f"n_source_trials={summary['n_source_trials']}\n"
            f"trial_averaging_mode={summary['trial_averaging']['mode']}\n"
            f"trials_per_average={summary['trial_averaging']['requested_value']}\n"
            f"n_times={summary['n_times']} @ {summary['sfreq_hz']} Hz\n"
            f"decode_window_samples={summary['decode_window_samples']}\n"
            f"n_decoded_windows={summary['n_decoded_windows']}\n"
            f"n_skipped_windows={summary['n_skipped_windows']}\n"
            f"n_decoder_channels={summary['n_decoder_channels']}\n"
            f"left_roi={summary['left_roi_mag_channels']}\n"
            f"right_roi={summary['right_roi_mag_channels']}\n"
            f"erf_baseline_window_sec={summary['erf_baseline_window_sec']}\n"
            f"balancing_strategy={summary['balancing_strategy']}\n"
            f"window_skip_reasons={summary['window_skip_reasons']}\n"
            f"peak_accuracy_unbalanced_train={summary['peak_accuracy_unbalanced_train']:.4f} "
            f"at t={summary['peak_time_sec_unbalanced_train']:.3f} s\n"
            f"mean_accuracy_unbalanced_train={summary['mean_accuracy_unbalanced_train']:.4f}\n"
            f"peak_accuracy_balanced_train={summary['peak_accuracy_balanced_train']:.4f} "
            f"at t={summary['peak_time_sec_balanced_train']:.3f} s\n"
            f"mean_accuracy_balanced_train={summary['mean_accuracy_balanced_train']:.4f}\n"
            f"scores_csv={summary['scores_csv']}\n"
            f"figure_path={summary['figure_path']}"
        )
        + "</pre>",
        title="Decoding summary",
        section="sanity_check",
        tags=("summary", "decoding", "sanity_check"),
    )

    fig_for_report = plt.figure(figsize=(12, 8))
    img = plt.imread(fig_path)
    plt.imshow(img)
    plt.axis("off")
    report.add_figure(
        fig_for_report,
        title="ERF + lateralization + decoding accuracy",
        section="sanity_check",
        tags=("figure", "decoding", "sanity_check"),
    )
    plt.close(fig_for_report)
    report.add_html(
        training_testing_diagnostics_html(summary),
        title="Training/testing diagnostics",
        section="sanity_check",
        tags=("summary", "decoding", "training_testing", "sanity_check"),
    )
    report.add_figure(
        fig_mag_topo,
        title="ERF Topography | Magnetometers (ipsi - contra)",
        section="sanity_check",
        tags=("figure", "erf", "magnetometers", "sanity_check"),
    )
    plt.close(fig_mag_topo)
    report.add_figure(
        fig_grad_topo_rms,
        title="ERF Topography | Planar Gradiometers (ipsi - contra, RMS pair merge)",
        section="sanity_check",
        tags=("figure", "erf", "gradiometers", "rms", "sanity_check"),
    )
    plt.close(fig_grad_topo_rms)

    report_path = output_paths["report_html"]
    report.save(report_path, overwrite=True, open_browser=False)

    print("Sanity-check decoding complete.")
    print(f"Figure: {fig_path}")
    print(f"Scores CSV: {scores_csv}")
    print(f"Summary JSON: {summary_json}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
