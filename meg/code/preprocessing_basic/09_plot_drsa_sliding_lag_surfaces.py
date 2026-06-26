#!/usr/bin/env python3
"""
09_plot_drsa_sliding_lag_surfaces.py

Purpose
Generate exploratory sliding-window lag summaries from existing Stage-07 dRSA
result MAT files. This stage does not run MATLAB and does not recompute dRSA.
It reads saved dRSA matrices, restricts model-time rows to short overlapping
windows, and computes lag curves inside each window using the same
model-anchored lag convention as Stage 08:

    lag = neural time - model time
    sample = matrix[model_time_index, model_time_index + lag]

The first implementation is intentionally narrow and report-focused. It targets
the group or subject occipital A_dRSA_alt position summaries that are most
useful for inspecting how latency-like dRSA structure changes over model time:
actual visible, actual nondeviant, actual deviant, expected/predicted deviant,
and prediction-error deviant.

Usage examples
1) Build the default group report for subjects 01-07 from existing occipital
   Stage-07 A_dRSA_alt outputs:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/09_plot_drsa_sliding_lag_surfaces.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 \
     --group-only --overwrite --print-summary

2) Use a wider model-time window and coarser step while keeping the default
   +/-1 s lag search:

   python .../09_plot_drsa_sliding_lag_surfaces.py --config .../pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 --group-only \
     --model-window-sec 0.300 --model-step-sec 0.050 --overwrite

3) Restrict the lag range and require at least five finite diagonal samples
   before plotting one window/lag cell:

   python .../09_plot_drsa_sliding_lag_surfaces.py --config .../pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 --group-only \
     --lag-window-sec 0.5 --min-valid-samples 5 --overwrite

4) Build a subject-level report for one participant:

   python .../09_plot_drsa_sliding_lag_surfaces.py --config .../pipeline_default.yaml \
     --subject-config .../sub-04.yaml --subjects-only --overwrite --print-summary

Inputs
- Existing Stage-07 result MAT files below
  derivatives/analysis/sub-XX/ses-01/drsa/.
- pipeline_config.py-resolved derivative paths and matlab_drsa.occlusion_timing.

Outputs
- Subject figures and interactive HTML:
  derivatives/analysis/sub-XX/ses-01/drsa/sliding_lag_surfaces/
  derivatives/reports/sub-XX/ses-01/09_drsa_sliding_lag_surfaces_report.html
  derivatives/reports/sub-XX/ses-01/09_drsa_sliding_lag_surfaces_report_embedded.html
- Group figures and interactive HTML:
  derivatives/analysis/group/ses-01/drsa/sliding_lag_surfaces/
  derivatives/reports/group/ses-01/09_drsa_sliding_lag_surfaces_report.html
  derivatives/reports/group/ses-01/09_drsa_sliding_lag_surfaces_report_embedded.html

Safety and assumptions
- Raw data are never modified. Every output path passes the shared raw-tree
  guard before writing.
- Stage-07 matrices are interpreted exactly as Stage 07 plots them: x is neural
  time, y is model time, and arrays are indexed
  matrix[model_time_index, neural_time_index].
- Sliding windows are model-anchored: each window restricts only model-time
  rows. Neural-time columns may fall anywhere inside the full matrix when the
  requested lag points outside the window but remains inside matrix bounds.
- Heatmaps and 3-D waterfalls share one value range across all generated target
  summaries in the same report so the five plots are visually comparable.
- The interactive 3-D waterfalls are written as lightweight Plotly JavaScript
  HTML. They do not require the Python plotly package, but they load Plotly JS
  from a CDN unless a local/alternative URL is supplied.
- These plots are descriptive. Overlapping windows, temporal autocorrelation,
  edge support, unequal valid sample counts, and peak-lag instability must be
  considered before treating apparent ridges as statistical evidence.
"""

from __future__ import annotations

import argparse
import html
import json
import textwrap
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from drsa_result_utils import (
    CONDITIONS,
    analysis_family_from_path,
    condition_from_path,
    drsa_type_from_path,
    html_id,
    image_href,
    iter_matrix_slices,
    model_source_for_matrix_slice,
    model_time_axis_label,
    model_source_from_path,
    normalize_sensor_group,
    normalize_subject_label,
    occlusion_timing_from_config,
    read_drsa_result_file,
    sampling_step_seconds,
    sanitize_label,
    sensor_group_from_path,
    stimulus_model_label,
    time_axis_for_matrix,
    utc_now,
    write_embedded_html_copy,
)
from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import guarded_path


DEFAULT_SAMPLING_RATE_HZ = 120.0
DEFAULT_PLOTLY_JS_URL = "https://cdn.plot.ly/plotly-2.35.2.min.js"
PLOT_TYPES = (
    "lag_by_model_time_heatmap",
    "plotly_3d_waterfall",
    "sliding_window_peak_timecourse",
)
VALUE_LIMIT_MODES = ("auto-symmetric", "auto", "manual")
PEAK_LATENCY_MASK_MODES = ("relative", "absolute", "none")


@dataclass(frozen=True)
class TargetSummarySpec:
    """
    One condition/model summary that should become one heatmap and one waterfall.

    Data flow:
    Stage-07 matrix metadata -> target selector -> one sliding lag surface.
    """

    slug: str
    condition: str
    stimulus_model_label: str
    title: str


TARGET_SUMMARIES = (
    TargetSummarySpec(
        slug="actual-visible-position",
        condition="always_visible",
        stimulus_model_label="actual visible position",
        title="Actual visible position",
    ),
    TargetSummarySpec(
        slug="actual-nondeviant-position",
        condition="occluded_nondeviant",
        stimulus_model_label="actual nondeviant position",
        title="Actual nondeviant position",
    ),
    TargetSummarySpec(
        slug="actual-deviant-position",
        condition="occluded_deviant",
        stimulus_model_label="actual deviant position",
        title="Actual deviant position",
    ),
    TargetSummarySpec(
        slug="expected-predicted-position",
        condition="occluded_deviant",
        stimulus_model_label="expected/predicted position",
        title="Expected/predicted deviant position",
    ),
    TargetSummarySpec(
        slug="prediction-error-position",
        condition="occluded_deviant",
        stimulus_model_label="prediction error position",
        title="Prediction-error deviant position",
    ),
)

PEAK_TIMECOURSE_COMPARISONS = (
    {
        "slug": "actual-position-comparison",
        "title": "Actual-position peak timecourse",
        "target_slugs": (
            "actual-visible-position",
            "actual-nondeviant-position",
            "actual-deviant-position",
        ),
    },
    {
        "slug": "occluded-deviant-model-comparison",
        "title": "Occluded-deviant model peak timecourse",
        "target_slugs": (
            "actual-deviant-position",
            "expected-predicted-position",
            "prediction-error-position",
        ),
    },
)


@dataclass(frozen=True)
class MatrixRecord:
    """
    One 2-D dRSA matrix slice plus analysis metadata.

    Data flow:
    Stage-07 result MAT variable -> 2-D matrix slice -> sliding lag surface.
    """

    owner_scope: str
    owner_label: str
    condition: str
    analysis_family: str
    sensor_group: str
    drsa_type: str
    model_source: str
    stimulus_model_label: str
    model_time_axis_label: str
    variable: str
    slice_index: int
    slice_label: str
    matrix: np.ndarray
    time_axis_sec: np.ndarray
    relative_result: str
    source_paths: tuple[str, ...]
    subjects: tuple[str, ...] = field(default_factory=tuple)
    n_subjects: int = 1
    subject_records: tuple["MatrixRecord", ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class SlidingLagSurface:
    """
    Sliding-window lag summary for one matrix record.

    ``values`` and ``valid_counts`` are indexed as
    [model_window_index, lag_index]. The model-time coordinate is selected by
    CLI option and can be the window center, start, or end.
    """

    record: MatrixRecord
    target: TargetSummarySpec
    lag_samples: np.ndarray
    lag_sec: np.ndarray
    window_starts_sec: np.ndarray
    window_ends_sec: np.ndarray
    window_coords_sec: np.ndarray
    values: np.ndarray
    valid_counts: np.ndarray
    model_window_samples: int
    model_step_samples: int


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the Stage-09 sliding lag-surface command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate sliding-window lag-by-model-time dRSA heatmaps and "
            "interactive Plotly 3-D waterfalls from existing Stage-07 MAT files."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--subjects", nargs="+", default=None)
    parser.add_argument("--group-only", action="store_true")
    parser.add_argument("--subjects-only", action="store_true")
    parser.add_argument(
        "--lag-window-sec",
        type=float,
        default=1.0,
        help="Maximum absolute lag, in seconds, for lag searches.",
    )
    parser.add_argument(
        "--model-window-sec",
        type=float,
        default=0.200,
        help="Sliding model-time window size in seconds.",
    )
    parser.add_argument(
        "--model-step-sec",
        type=float,
        default=0.025,
        help="Step between consecutive model-time windows in seconds.",
    )
    parser.add_argument(
        "--window-coordinate",
        choices=("center", "start", "end"),
        default="center",
        help="Model-time coordinate used on the heatmap/waterfall depth axis.",
    )
    parser.add_argument(
        "--min-valid-samples",
        type=int,
        default=3,
        help=(
            "Minimum finite diagonal samples required for each "
            "window-lag cell. Cells below this count are plotted as missing."
        ),
    )
    parser.add_argument(
        "--sensor-group",
        default="occipital",
        help="Sensor group to plot. The first implementation defaults to occipital.",
    )
    parser.add_argument(
        "--analysis-family",
        default="A_dRSA_alt",
        help="Stage-07 analysis family to plot. Defaults to A_dRSA_alt.",
    )
    parser.add_argument(
        "--drsa-type",
        default="corr",
        help="dRSA metric folder to plot. Defaults to corr.",
    )
    parser.add_argument(
        "--variable",
        default="dRSA_neural_position",
        help="MATLAB matrix variable to plot. Defaults to dRSA_neural_position.",
    )
    parser.add_argument(
        "--value-limit-mode",
        choices=VALUE_LIMIT_MODES,
        default="auto-symmetric",
        help=(
            "Shared value-limit mode for all heatmaps and waterfall z axes: "
            "auto-symmetric, auto, or manual."
        ),
    )
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument(
        "--peak-latency-mask-mode",
        choices=PEAK_LATENCY_MASK_MODES,
        default="relative",
        help=(
            "How to mask peak-lag/latency timecourses when peak magnitude is "
            "weak: relative, absolute, or none."
        ),
    )
    parser.add_argument(
        "--peak-latency-min-magnitude-frac",
        type=float,
        default=0.20,
        help=(
            "Relative peak-latency mask threshold as a fraction of the maximum "
            "finite peak magnitude across report surfaces."
        ),
    )
    parser.add_argument(
        "--peak-latency-min-magnitude",
        type=float,
        default=None,
        help=(
            "Absolute peak-latency mask threshold. Required when "
            "--peak-latency-mask-mode absolute is used."
        ),
    )
    parser.add_argument(
        "--plotly-js-url",
        default=DEFAULT_PLOTLY_JS_URL,
        help=(
            "Plotly JavaScript URL for interactive waterfalls. This can be a "
            "CDN URL or a report-relative local JS path."
        ),
    )
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


def _analysis_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return subject and group Stage-09 input/output paths from config."""
    resolved = cfg["resolved_paths"]
    analysis_dir = Path(resolved["derivatives_analysis_dir"])
    reports_dir = Path(resolved["derivatives_reports_dir"])
    session = cfg["layout"]["session_bids_id"]
    group_analysis_dir = analysis_dir.parent.parent / "group" / session
    group_reports_dir = reports_dir.parent.parent / "group" / session
    return {
        "analysis_dir": analysis_dir,
        "reports_dir": reports_dir,
        "drsa_dir": analysis_dir / "drsa",
        "subject_surface_dir": analysis_dir / "drsa" / "sliding_lag_surfaces",
        "subject_report_html": reports_dir / "09_drsa_sliding_lag_surfaces_report.html",
        "group_analysis_dir": group_analysis_dir,
        "group_reports_dir": group_reports_dir,
        "group_drsa_dir": group_analysis_dir / "drsa",
        "group_surface_dir": group_analysis_dir / "drsa" / "sliding_lag_surfaces",
        "group_report_html": group_reports_dir / "09_drsa_sliding_lag_surfaces_report.html",
    }


def _load_cfg(
    args: argparse.Namespace,
    *,
    subject_config: str | Path | None = None,
) -> Mapping[str, Any]:
    """Load one resolved pipeline configuration for Stage-09 use."""
    cfg, warnings_list, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("export", "matlab_drsa"),
    )
    if warnings_list and args.print_summary:
        print(json.dumps({"configuration_warnings": warnings_list}, indent=2))
    return cfg


def _subject_config_for_id(args: argparse.Namespace, subject_label: str) -> Path | None:
    """Resolve the per-subject config file for a requested subject."""
    if args.subject_config and not args.subjects:
        return Path(args.subject_config)
    if args.subject_config and args.subjects and len(args.subjects) == 1:
        return Path(args.subject_config)
    subject_id = subject_label.replace("sub-", "")
    candidate = Path(args.config).expanduser().resolve().parent / f"sub-{subject_id}.yaml"
    return candidate if candidate.is_file() else None


def _subjects_from_args_or_config(
    args: argparse.Namespace,
    cfg: Mapping[str, Any],
) -> list[str]:
    """Resolve requested subject labels for subject or group operations."""
    if args.subjects:
        return [normalize_subject_label(value) for value in args.subjects]
    configured = cfg.get("matlab_drsa", {}).get("group_average_subjects")
    if args.group_only and configured:
        return [normalize_subject_label(value) for value in configured]
    return [str(cfg["layout"]["subject_bids_id"])]


def _sampling_rate_for_plots(cfg: Mapping[str, Any]) -> float:
    """Resolve the fallback sampling rate used when a MAT file lacks time_axis_sec."""
    return float(
        cfg.get("export", {}).get(
            "sampling_rate_hz",
            cfg.get("matlab_drsa", {}).get(
                "plot_sample_rate_hz",
                DEFAULT_SAMPLING_RATE_HZ,
            ),
        )
    )


def _should_skip_result_path(path: Path) -> bool:
    """Filter out non-result MAT files and generated plotting output folders."""
    text = str(path)
    return (
        "desc-identityAveraged_neural" in path.name
        or "/line_summaries/" in text
        or "line_summaries" in path.parts
        or "/sliding_lag_surfaces/" in text
        or "sliding_lag_surfaces" in path.parts
    )


def _read_subject_records(
    cfg: Mapping[str, Any],
    *,
    owner_label: str,
    sensor_group: str,
    analysis_family: str,
    drsa_type: str,
    variable: str,
) -> tuple[list[MatrixRecord], list[dict[str, str]]]:
    """
    Discover and read Stage-07 target matrices for one subject.

    Data flow:
    subject drsa folder -> result MAT files -> selected matrix variable ->
    MatrixRecord objects for target summary selection.
    """
    paths = _analysis_paths(cfg)
    sampling_rate = _sampling_rate_for_plots(cfg)
    drsa_dir = paths["drsa_dir"]
    normalized_sensor = normalize_sensor_group(sensor_group)
    records: list[MatrixRecord] = []
    skipped: list[dict[str, str]] = []
    if not drsa_dir.is_dir():
        skipped.append({"owner": owner_label, "reason": f"missing dRSA folder: {drsa_dir}"})
        return records, skipped

    for mat_path in sorted(drsa_dir.rglob("*.mat")):
        if _should_skip_result_path(mat_path):
            continue
        relative = mat_path.relative_to(drsa_dir)
        path_analysis_family = analysis_family_from_path(relative)
        path_sensor_group = sensor_group_from_path(relative)
        path_drsa_type = drsa_type_from_path(relative)
        if path_analysis_family != analysis_family:
            continue
        if path_sensor_group != normalized_sensor:
            continue
        if path_drsa_type != drsa_type:
            continue

        condition = condition_from_path(relative)
        path_model_source = model_source_from_path(relative)
        try:
            matrices, time_axis, result_metadata = read_drsa_result_file(mat_path)
        except RuntimeError as exc:
            skipped.append(
                {
                    "owner": owner_label,
                    "path": str(mat_path),
                    "reason": str(exc),
                }
            )
            continue
        if variable not in matrices:
            continue

        matrix_array = matrices[variable]
        model_names = (
            result_metadata.get("model_names")
            if variable == "dRSA_neural_correlation"
            else None
        )
        for slice_index, slice_label, matrix_2d in iter_matrix_slices(
            matrix_array,
            slice_labels=model_names,
        ):
            effective_model_source = model_source_for_matrix_slice(
                condition,
                path_model_source,
                variable,
                slice_label=slice_label,
                analysis_family=path_analysis_family,
            )
            model_label = stimulus_model_label(
                condition,
                path_model_source,
                variable,
                slice_label=slice_label,
                analysis_family=path_analysis_family,
            )
            axis_label = model_time_axis_label(
                condition,
                path_model_source,
                variable,
                slice_label=slice_label,
                analysis_family=path_analysis_family,
            )
            axis = time_axis_for_matrix(
                matrix_2d,
                time_axis=time_axis,
                sampling_rate_hz=sampling_rate,
            )
            records.append(
                MatrixRecord(
                    owner_scope="subject",
                    owner_label=owner_label,
                    condition=condition,
                    analysis_family=path_analysis_family,
                    sensor_group=path_sensor_group,
                    drsa_type=path_drsa_type,
                    model_source=effective_model_source,
                    stimulus_model_label=model_label,
                    model_time_axis_label=axis_label,
                    variable=variable,
                    slice_index=int(slice_index),
                    slice_label=slice_label,
                    matrix=np.asarray(matrix_2d, dtype=float),
                    time_axis_sec=np.asarray(axis, dtype=float),
                    relative_result=str(relative),
                    source_paths=(str(mat_path),),
                    subjects=(owner_label,),
                    n_subjects=1,
                )
            )
    return records, skipped


def _group_records(
    subject_records: Sequence[MatrixRecord],
) -> tuple[list[MatrixRecord], list[dict[str, str]]]:
    """
    Average matching matrix records across subjects.

    Records match when they share the same Stage-07 relative result path,
    variable, and slice. Shape mismatches are skipped with report metadata.
    """
    grouped: dict[tuple[str, str, int], list[MatrixRecord]] = {}
    for record in subject_records:
        key = (record.relative_result, record.variable, record.slice_index)
        grouped.setdefault(key, []).append(record)

    group_records: list[MatrixRecord] = []
    skipped: list[dict[str, str]] = []
    for key, items in sorted(grouped.items()):
        shapes = {tuple(item.matrix.shape) for item in items}
        if len(shapes) != 1:
            skipped.append(
                {
                    "relative_result": key[0],
                    "variable": key[1],
                    "reason": f"shape mismatch across subjects: {sorted(shapes)}",
                }
            )
            continue
        stack = np.stack([item.matrix for item in items], axis=0)
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            averaged = np.nanmean(stack, axis=0)
        first = items[0]
        subjects = tuple(item.owner_label for item in items)
        group_records.append(
            MatrixRecord(
                owner_scope="group",
                owner_label="group",
                condition=first.condition,
                analysis_family=first.analysis_family,
                sensor_group=first.sensor_group,
                drsa_type=first.drsa_type,
                model_source=first.model_source,
                stimulus_model_label=first.stimulus_model_label,
                model_time_axis_label=first.model_time_axis_label,
                variable=first.variable,
                slice_index=first.slice_index,
                slice_label=first.slice_label,
                matrix=averaged,
                time_axis_sec=first.time_axis_sec,
                relative_result=first.relative_result,
                source_paths=tuple(path for item in items for path in item.source_paths),
                subjects=subjects,
                n_subjects=len(items),
                subject_records=tuple(items),
            )
        )
    return group_records, skipped


def _target_spec_for_record(record: MatrixRecord) -> TargetSummarySpec | None:
    """Return the target summary spec matched by a record, if any."""
    for spec in TARGET_SUMMARIES:
        if (
            record.condition == spec.condition
            and record.stimulus_model_label == spec.stimulus_model_label
        ):
            return spec
    return None


def _target_records(
    records: Sequence[MatrixRecord],
) -> tuple[list[tuple[TargetSummarySpec, MatrixRecord]], list[dict[str, str]]]:
    """Select the five requested target records in stable report order."""
    by_slug: dict[str, list[MatrixRecord]] = {}
    for record in records:
        spec = _target_spec_for_record(record)
        if spec is not None:
            by_slug.setdefault(spec.slug, []).append(record)

    selected: list[tuple[TargetSummarySpec, MatrixRecord]] = []
    skipped: list[dict[str, str]] = []
    for spec in TARGET_SUMMARIES:
        matches = by_slug.get(spec.slug, [])
        if not matches:
            skipped.append(
                {
                    "target": spec.slug,
                    "condition": spec.condition,
                    "model": spec.stimulus_model_label,
                    "reason": "target matrix not found",
                }
            )
            continue
        matches = sorted(
            matches,
            key=lambda item: (
                item.relative_result,
                item.variable,
                item.slice_index,
                item.slice_label,
            ),
        )
        if len(matches) > 1:
            skipped.append(
                {
                    "target": spec.slug,
                    "reason": (
                        f"multiple matches found; using {matches[0].relative_result} "
                        f"slice {matches[0].slice_label}"
                    ),
                }
            )
        selected.append((spec, matches[0]))
    return selected, skipped


def _lag_samples(time_axis_sec: np.ndarray, lag_window_sec: float) -> np.ndarray:
    """Build integer lag offsets for the requested lag window."""
    dt = sampling_step_seconds(time_axis_sec, DEFAULT_SAMPLING_RATE_HZ)
    max_lag = int(round(float(lag_window_sec) / dt))
    return np.arange(-max_lag, max_lag + 1, dtype=int)


def _sample_count_from_seconds(
    time_axis_sec: np.ndarray,
    seconds: float,
    *,
    minimum: int = 1,
) -> int:
    """Convert a duration in seconds to a positive sample count."""
    dt = sampling_step_seconds(time_axis_sec, DEFAULT_SAMPLING_RATE_HZ)
    return max(int(minimum), int(round(float(seconds) / dt)))


def _nanmean_with_min_count(values: np.ndarray, min_valid_samples: int) -> tuple[float, int]:
    """Return a finite-sample mean and count, preserving NaN below threshold."""
    finite = np.isfinite(values)
    count = int(np.sum(finite))
    if count < int(min_valid_samples):
        return float("nan"), count
    return float(np.nanmean(values[finite])), count


def compute_sliding_lag_surface(
    matrix_2d: Any,
    time_axis_sec: Any,
    *,
    lag_window_sec: float,
    model_window_sec: float,
    model_step_sec: float,
    window_coordinate: str = "center",
    min_valid_samples: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    """
    Compute a model-anchored sliding-window lag surface.

    Returns
    -------
    lag_samples, lag_sec, window_coords_sec, values, valid_counts,
    model_window_samples, model_step_samples

    Data flow:
    matrix[model, neural] -> overlapping model-time row windows -> lag-indexed
    diagonal samples -> mean dRSA value and finite sample count for each
    window/lag cell.
    """
    if window_coordinate not in {"center", "start", "end"}:
        raise ValueError("window_coordinate must be center, start, or end.")
    if lag_window_sec <= 0:
        raise ValueError("lag_window_sec must be > 0.")
    if model_window_sec <= 0:
        raise ValueError("model_window_sec must be > 0.")
    if model_step_sec <= 0:
        raise ValueError("model_step_sec must be > 0.")
    if min_valid_samples < 1:
        raise ValueError("min_valid_samples must be >= 1.")

    matrix = np.asarray(matrix_2d, dtype=float)
    axis = np.asarray(time_axis_sec, dtype=float).squeeze()
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2-D matrix, found shape {matrix.shape}.")
    if axis.ndim != 1 or axis.size != matrix.shape[0]:
        raise ValueError(
            "time_axis_sec must be a 1-D array matching the model-time dimension."
        )

    n_model, n_neural = matrix.shape
    lag_samples = _lag_samples(axis, lag_window_sec)
    dt = sampling_step_seconds(axis, DEFAULT_SAMPLING_RATE_HZ)
    lag_sec = lag_samples.astype(float) * dt
    model_window_samples = _sample_count_from_seconds(axis, model_window_sec)
    model_step_samples = _sample_count_from_seconds(axis, model_step_sec)
    if model_window_samples > n_model:
        raise ValueError(
            f"model window has {model_window_samples} samples, but matrix has "
            f"only {n_model} model-time rows."
        )

    starts = np.arange(
        0,
        n_model - model_window_samples + 1,
        model_step_samples,
        dtype=int,
    )
    ends = starts + model_window_samples
    window_starts_sec = axis[starts]
    window_ends_sec = axis[ends - 1]
    if window_coordinate == "start":
        window_coords_sec = window_starts_sec
    elif window_coordinate == "end":
        window_coords_sec = window_ends_sec
    else:
        window_coords_sec = (window_starts_sec + window_ends_sec) / 2.0

    values = np.full((starts.size, lag_samples.size), np.nan, dtype=float)
    valid_counts = np.zeros((starts.size, lag_samples.size), dtype=int)
    for window_idx, (start, end) in enumerate(zip(starts, ends)):
        model_indices = np.arange(start, end, dtype=int)
        for lag_idx, lag in enumerate(lag_samples):
            neural_indices = model_indices + int(lag)
            in_bounds = (neural_indices >= 0) & (neural_indices < n_neural)
            if not np.any(in_bounds):
                continue
            samples = matrix[model_indices[in_bounds], neural_indices[in_bounds]]
            value, count = _nanmean_with_min_count(samples, min_valid_samples)
            values[window_idx, lag_idx] = value
            valid_counts[window_idx, lag_idx] = count

    return (
        lag_samples,
        lag_sec,
        window_coords_sec,
        values,
        valid_counts,
        model_window_samples,
        model_step_samples,
    )


def _surface_for_record(
    record: MatrixRecord,
    target: TargetSummarySpec,
    *,
    args: argparse.Namespace,
) -> SlidingLagSurface:
    """Compute the sliding lag surface and attach record/target metadata."""
    (
        lag_samples,
        lag_sec,
        window_coords_sec,
        values,
        valid_counts,
        model_window_samples,
        model_step_samples,
    ) = compute_sliding_lag_surface(
        record.matrix,
        record.time_axis_sec,
        lag_window_sec=args.lag_window_sec,
        model_window_sec=args.model_window_sec,
        model_step_sec=args.model_step_sec,
        window_coordinate=args.window_coordinate,
        min_valid_samples=args.min_valid_samples,
    )
    starts = np.arange(
        0,
        record.matrix.shape[0] - model_window_samples + 1,
        model_step_samples,
        dtype=int,
    )
    ends = starts + model_window_samples
    return SlidingLagSurface(
        record=record,
        target=target,
        lag_samples=lag_samples,
        lag_sec=lag_sec,
        window_starts_sec=record.time_axis_sec[starts],
        window_ends_sec=record.time_axis_sec[ends - 1],
        window_coords_sec=window_coords_sec,
        values=values,
        valid_counts=valid_counts,
        model_window_samples=model_window_samples,
        model_step_samples=model_step_samples,
    )


def _float_slug(value: float) -> str:
    """Return a compact filename token for a floating-point second value."""
    text = f"{float(value):.3f}".rstrip("0").rstrip(".")
    if text == "-0":
        text = "0"
    return text.replace("-", "neg").replace(".", "p")


def _surface_slug(surface: SlidingLagSurface, args: argparse.Namespace) -> str:
    """Build a stable output slug for one surface and parameter set."""
    record = surface.record
    parts = [
        record.owner_label,
        "raw",
        record.analysis_family,
        record.sensor_group,
        record.drsa_type,
        surface.target.slug,
        record.variable,
        f"win-{_float_slug(args.model_window_sec)}s",
        f"step-{_float_slug(args.model_step_sec)}s",
        f"lag-{_float_slug(args.lag_window_sec)}s",
        f"minvalid-{int(args.min_valid_samples)}",
    ]
    return "_".join(sanitize_label(part) for part in parts if part)


def _prepare_output_path(path: Path, *, raw_root: Path, args: argparse.Namespace) -> Path | None:
    """
    Guard one output path and apply --resume/--overwrite behavior.

    Returns ``None`` when an existing output should be reused in resume mode.
    """
    guarded = guarded_path(path, raw_root)
    if guarded.exists():
        if args.resume:
            return None
        if not args.overwrite:
            raise FileExistsError(
                f"Output already exists: {guarded}. Use --overwrite or --resume."
            )
    guarded.parent.mkdir(parents=True, exist_ok=True)
    return guarded


def _value_limits(
    surfaces: Sequence[SlidingLagSurface],
    *,
    mode: str,
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float]:
    """Resolve shared value limits for all generated plots in one report."""
    if mode == "manual":
        if vmin is None or vmax is None:
            raise ValueError("--value-limit-mode manual requires --vmin and --vmax.")
        if not vmin < vmax:
            raise ValueError("--vmin must be smaller than --vmax.")
        return float(vmin), float(vmax)

    finite_arrays = [
        surface.values[np.isfinite(surface.values)]
        for surface in surfaces
        if np.any(np.isfinite(surface.values))
    ]
    if not finite_arrays:
        return -1.0, 1.0
    values = np.concatenate(finite_arrays)
    data_min = float(np.nanmin(values))
    data_max = float(np.nanmax(values))
    if mode == "auto-symmetric":
        limit = max(abs(data_min), abs(data_max))
        if limit == 0.0:
            limit = 1.0
        return -limit, limit
    if data_min == data_max:
        padding = max(abs(data_min) * 0.05, 1e-6)
        return data_min - padding, data_max + padding
    return data_min, data_max


def extract_peak_timecourse(surface: SlidingLagSurface) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract peak magnitude and peak lag for each model-time window.

    Data flow:
    lag-by-model-time surface -> per-window argmax across lag -> peak magnitude
    and its lag. Empty/all-NaN windows remain NaN in both returned arrays.
    """
    peak_magnitudes = np.full(surface.values.shape[0], np.nan, dtype=float)
    peak_lags = np.full(surface.values.shape[0], np.nan, dtype=float)
    for row_idx, row in enumerate(surface.values):
        if np.any(np.isfinite(row)):
            peak_idx = int(np.nanargmax(row))
            peak_magnitudes[row_idx] = float(row[peak_idx])
            peak_lags[row_idx] = float(surface.lag_sec[peak_idx])
    return peak_magnitudes, peak_lags


def _peak_latency_threshold(
    surfaces: Sequence[SlidingLagSurface],
    *,
    args: argparse.Namespace,
) -> float | None:
    """Resolve the magnitude threshold used to mask peak-lag timecourses."""
    if args.peak_latency_mask_mode == "none":
        return None
    if args.peak_latency_mask_mode == "absolute":
        if args.peak_latency_min_magnitude is None:
            raise ValueError(
                "--peak-latency-mask-mode absolute requires "
                "--peak-latency-min-magnitude."
            )
        return float(args.peak_latency_min_magnitude)

    finite_peaks: list[np.ndarray] = []
    for surface in surfaces:
        peak_magnitudes, _peak_lags = extract_peak_timecourse(surface)
        finite = peak_magnitudes[np.isfinite(peak_magnitudes)]
        if finite.size:
            finite_peaks.append(finite)
    if not finite_peaks:
        return float("inf")
    max_peak = float(np.nanmax(np.concatenate(finite_peaks)))
    if max_peak <= 0.0:
        return float("inf")
    return max_peak * float(args.peak_latency_min_magnitude_frac)


def _compact_plot_title(parts: Iterable[Any], *, width: int = 90) -> str:
    """Return a wrapped Matplotlib title from plot-identifying fields."""
    text = " | ".join(str(part).strip() for part in parts if str(part).strip())
    if not text:
        return ""
    return "\n".join(
        textwrap.wrap(
            text,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )


def _plot_heatmap(
    surface: SlidingLagSurface,
    *,
    output_dir: Path,
    raw_root: Path,
    occlusion_timing: Mapping[str, Any],
    value_limits: tuple[float, float],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Write one static lag-by-model-time heatmap PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    heatmap_dir = output_dir / "heatmaps"
    out_path = heatmap_dir / f"{_surface_slug(surface, args)}_heatmap.png"
    write_target = _prepare_output_path(out_path, raw_root=raw_root, args=args)

    if write_target is not None:
        cmap = plt.get_cmap("coolwarm").copy()
        cmap.set_bad("#eeeeee")
        masked_values = np.ma.masked_invalid(surface.values)
        fig, ax = plt.subplots(figsize=(9.5, 6.2), constrained_layout=True)
        image = ax.imshow(
            masked_values,
            origin="lower",
            aspect="auto",
            extent=[
                float(surface.lag_sec[0]),
                float(surface.lag_sec[-1]),
                float(surface.window_coords_sec[0]),
                float(surface.window_coords_sec[-1]),
            ],
            vmin=value_limits[0],
            vmax=value_limits[1],
            cmap=cmap,
            interpolation="nearest",
        )
        ax.axvline(0.0, color="0.15", linestyle="--", linewidth=1.0)
        first = float(occlusion_timing["first_complete_occlusion_sec"])
        reappearance = float(occlusion_timing["reappearance_sec"])
        for boundary, label in (
            (first, "complete occlusion"),
            (reappearance, "reappearance"),
        ):
            if surface.window_coords_sec[0] <= boundary <= surface.window_coords_sec[-1]:
                ax.axhline(boundary, color="0.2", linestyle=":", linewidth=1.0)
                ax.text(
                    surface.lag_sec[0],
                    boundary,
                    f" {label}",
                    va="bottom",
                    ha="left",
                    fontsize=7,
                    color="0.2",
                )

        _peak_magnitudes, peak_lags = extract_peak_timecourse(surface)
        finite_peaks = np.isfinite(peak_lags)
        if np.any(finite_peaks):
            ax.scatter(
                peak_lags[finite_peaks],
                surface.window_coords_sec[finite_peaks],
                color="black",
                s=8,
                linewidths=0.0,
                alpha=0.65,
                label="descriptive peak lag points",
            )
            ax.legend(loc="upper right", fontsize=8, frameon=True)

        ax.set_xlabel("Lag (neural time - model time) [s]")
        ax.set_ylabel(f"Model-time window {args.window_coordinate} [s]")
        ax.set_title(
            _compact_plot_title(
                [
                    surface.record.owner_label,
                    surface.target.title,
                    surface.record.sensor_group,
                    surface.record.analysis_family,
                    surface.record.variable,
                ]
            ),
            fontsize=10,
        )
        cbar = fig.colorbar(image, ax=ax)
        cbar.set_label("Mean dRSA value")
        fig.savefig(write_target, dpi=180, facecolor="white")
        plt.close(fig)

    finite_counts = surface.valid_counts[surface.valid_counts > 0]
    return {
        "plot_type": "lag_by_model_time_heatmap",
        "figure": str(out_path),
        "owner": surface.record.owner_label,
        "target": surface.target.slug,
        "condition": surface.record.condition,
        "model_source": surface.record.model_source,
        "stimulus_model_label": surface.record.stimulus_model_label,
        "analysis_family": surface.record.analysis_family,
        "sensor_group": surface.record.sensor_group,
        "variable": surface.record.variable,
        "slice_label": surface.record.slice_label,
        "caption": (
            f"{surface.record.owner_label} | {surface.target.title} | "
            f"{surface.record.analysis_family} | {surface.record.sensor_group} | "
            f"{surface.record.variable} | window={args.model_window_sec:g}s, "
            f"step={args.model_step_sec:g}s, lag=+/-{args.lag_window_sec:g}s | "
            f"minimum finite samples per cell={args.min_valid_samples} | "
            f"valid-count range="
            f"{int(np.min(finite_counts)) if finite_counts.size else 0}-"
            f"{int(np.max(finite_counts)) if finite_counts.size else 0}"
        ),
    }


def _finite_or_none(value: Any) -> float | None:
    """Convert one numeric value to JSON-safe float/None."""
    numeric = float(value)
    if not np.isfinite(numeric):
        return None
    return numeric


def _plotly_data_for_surface(surface: SlidingLagSurface) -> list[dict[str, Any]]:
    """Build Plotly scatter3d traces for one sliding lag surface."""
    traces: list[dict[str, Any]] = []
    coord_min = float(np.nanmin(surface.window_coords_sec))
    coord_max = float(np.nanmax(surface.window_coords_sec))
    for row_idx, coord in enumerate(surface.window_coords_sec):
        row = surface.values[row_idx]
        if not np.any(np.isfinite(row)):
            continue
        coord_value = float(coord)
        hover_text = []
        for lag, value, count in zip(
            surface.lag_sec,
            row,
            surface.valid_counts[row_idx],
        ):
            approximate_neural = coord_value + float(lag)
            value_text = "missing" if not np.isfinite(value) else f"{float(value):.5g}"
            hover_text.append(
                "<br>".join(
                    [
                        html.escape(surface.target.title),
                        f"model window center/start/end: {coord_value:.4f} s",
                        (
                            "window span: "
                            f"{surface.window_starts_sec[row_idx]:.4f}-"
                            f"{surface.window_ends_sec[row_idx]:.4f} s"
                        ),
                        f"lag: {float(lag):.4f} s",
                        f"approx neural time: {approximate_neural:.4f} s",
                        f"dRSA value: {value_text}",
                        f"finite samples: {int(count)}",
                    ]
                )
            )
        traces.append(
            {
                "type": "scatter3d",
                "mode": "lines",
                "name": f"{coord_value:.3f} s",
                "showlegend": False,
                "x": [_finite_or_none(value) for value in surface.lag_sec],
                "y": [coord_value for _ in surface.lag_sec],
                "z": [_finite_or_none(value) for value in row],
                "hoverinfo": "text",
                "hovertext": hover_text,
                "line": {
                    "width": 3,
                    "color": [coord_value for _ in surface.lag_sec],
                    "colorscale": "Viridis",
                    "cmin": coord_min,
                    "cmax": coord_max,
                },
            }
        )
    if not traces:
        traces.append(
            {
                "type": "scatter3d",
                "mode": "markers",
                "name": "no finite windows",
                "x": [0.0],
                "y": [float(surface.window_coords_sec[0])],
                "z": [0.0],
                "marker": {"size": 1, "opacity": 0.0},
                "hoverinfo": "text",
                "hovertext": ["No finite sliding-window lag cells passed the mask."],
                "showlegend": False,
            }
        )
    return traces


def _plotly_layout_for_surface(
    surface: SlidingLagSurface,
    *,
    value_limits: tuple[float, float],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Build a Plotly layout shared by report-embedded and standalone plots."""
    title = (
        f"{surface.record.owner_label} | {surface.target.title} | "
        f"{surface.record.sensor_group} | {surface.record.analysis_family}"
    )
    return {
        "title": {"text": title, "font": {"size": 14}},
        "margin": {"l": 0, "r": 0, "b": 0, "t": 42},
        "scene": {
            "xaxis": {"title": "Lag (neural - model) [s]"},
            "yaxis": {"title": f"Model-time window {args.window_coordinate} [s]"},
            "zaxis": {
                "title": "Mean dRSA value",
                "range": [float(value_limits[0]), float(value_limits[1])],
            },
            "camera": {"eye": {"x": 1.7, "y": -1.8, "z": 1.05}},
        },
        "height": 640,
    }


def _plotly_div_html(entry: Mapping[str, Any]) -> str:
    """Return an in-report Plotly div and script block for one waterfall."""
    div_id = str(entry["div_id"])
    data_json = json.dumps(entry["plotly_data"], separators=(",", ":"))
    layout_json = json.dumps(entry["plotly_layout"], separators=(",", ":"))
    return "\n".join(
        [
            f"<div id=\"{html.escape(div_id)}\" class=\"plotly-chart\"></div>",
            "<script>",
            (
                "Plotly.newPlot("
                f"{json.dumps(div_id)}, {data_json}, {layout_json}, "
                "{responsive: true, displaylogo: false}"
                ");"
            ),
            "</script>",
        ]
    )


def _standalone_plotly_html(
    *,
    title: str,
    plotly_js_url: str,
    div_html: str,
) -> str:
    """Build a standalone HTML document for one interactive waterfall."""
    return "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            f"<script src=\"{html.escape(plotly_js_url)}\"></script>",
            "<style>",
            "body{font-family:system-ui,sans-serif;margin:1rem;}",
            ".plotly-chart{width:100%;height:680px;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            div_html,
            "</body></html>",
        ]
    )


def _write_waterfall(
    surface: SlidingLagSurface,
    *,
    output_dir: Path,
    raw_root: Path,
    value_limits: tuple[float, float],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Write one standalone Plotly 3-D waterfall HTML file and return report data."""
    waterfall_dir = output_dir / "plotly_waterfalls"
    out_path = waterfall_dir / f"{_surface_slug(surface, args)}_waterfall.html"
    write_target = _prepare_output_path(out_path, raw_root=raw_root, args=args)
    div_id = f"waterfall-{sanitize_label(_surface_slug(surface, args))}"
    plotly_data = _plotly_data_for_surface(surface)
    plotly_layout = _plotly_layout_for_surface(
        surface,
        value_limits=value_limits,
        args=args,
    )
    div_entry = {
        "div_id": div_id,
        "plotly_data": plotly_data,
        "plotly_layout": plotly_layout,
    }
    if write_target is not None:
        title = (
            f"{surface.record.owner_label} | {surface.target.title} | "
            "3-D sliding lag waterfall"
        )
        write_target.write_text(
            _standalone_plotly_html(
                title=title,
                plotly_js_url=args.plotly_js_url,
                div_html=_plotly_div_html(div_entry),
            ),
            encoding="utf-8",
        )

    return {
        "plot_type": "plotly_3d_waterfall",
        "html": str(out_path),
        "owner": surface.record.owner_label,
        "target": surface.target.slug,
        "condition": surface.record.condition,
        "model_source": surface.record.model_source,
        "stimulus_model_label": surface.record.stimulus_model_label,
        "analysis_family": surface.record.analysis_family,
        "sensor_group": surface.record.sensor_group,
        "variable": surface.record.variable,
        "slice_label": surface.record.slice_label,
        "caption": (
            f"{surface.record.owner_label} | {surface.target.title} | "
            f"{surface.record.analysis_family} | {surface.record.sensor_group} | "
            f"{surface.record.variable} | rotatable 3-D lines; x=lag, "
            f"depth=model-time window {args.window_coordinate}, z=mean dRSA."
        ),
        "div_id": div_id,
        "plotly_data": plotly_data,
        "plotly_layout": plotly_layout,
    }


def _comparison_slug(
    comparison: Mapping[str, Any],
    *,
    owner_label: str,
    args: argparse.Namespace,
) -> str:
    """Build a stable output slug for a peak-timecourse comparison figure."""
    parts = [
        owner_label,
        str(comparison["slug"]),
        f"win-{_float_slug(args.model_window_sec)}s",
        f"step-{_float_slug(args.model_step_sec)}s",
        f"lag-{_float_slug(args.lag_window_sec)}s",
        f"minvalid-{int(args.min_valid_samples)}",
    ]
    return "_".join(sanitize_label(part) for part in parts if part)


def _line_label_for_surface(surface: SlidingLagSurface) -> str:
    """Return a concise legend label for comparison timecourses."""
    return f"{surface.record.condition} | {surface.record.stimulus_model_label}"


def _add_occlusion_markers(
    ax: Any,
    *,
    occlusion_timing: Mapping[str, Any],
) -> None:
    """Add vertical occlusion-boundary markers to a timecourse axis."""
    first = float(occlusion_timing["first_complete_occlusion_sec"])
    reappearance = float(occlusion_timing["reappearance_sec"])
    for boundary, label in (
        (first, "complete occlusion"),
        (reappearance, "reappearance"),
    ):
        ax.axvline(boundary, color="0.35", linestyle=":", linewidth=1.0)
        ax.text(
            boundary,
            0.98,
            label,
            transform=ax.get_xaxis_transform(),
            rotation=90,
            va="top",
            ha="right",
            fontsize=7,
            color="0.3",
        )


def _plot_peak_timecourse_comparison(
    surfaces: Sequence[SlidingLagSurface],
    *,
    comparison: Mapping[str, Any],
    output_dir: Path,
    raw_root: Path,
    occlusion_timing: Mapping[str, Any],
    peak_latency_threshold: float | None,
    args: argparse.Namespace,
) -> dict[str, Any] | None:
    """
    Write one two-panel peak magnitude/latency comparison figure.

    Data flow:
    selected lag surfaces -> per-window peak magnitude and peak lag -> overlay
    comparison across conditions/models. Peak lag is masked when the peak
    magnitude is below the configured threshold.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    wanted = tuple(str(value) for value in comparison["target_slugs"])
    surface_by_slug = {surface.target.slug: surface for surface in surfaces}
    selected = [surface_by_slug[slug] for slug in wanted if slug in surface_by_slug]
    if not selected:
        return None

    peak_dir = output_dir / "peak_timecourses"
    out_path = peak_dir / (
        f"{_comparison_slug(comparison, owner_label=selected[0].record.owner_label, args=args)}"
        "_peak_timecourse.png"
    )
    write_target = _prepare_output_path(out_path, raw_root=raw_root, args=args)

    if write_target is not None:
        fig, (mag_ax, lag_ax) = plt.subplots(
            nrows=2,
            sharex=True,
            figsize=(10.5, 7.2),
            constrained_layout=True,
        )
        for surface in selected:
            peak_magnitudes, peak_lags = extract_peak_timecourse(surface)
            if peak_latency_threshold is None:
                masked_peak_lags = peak_lags
            else:
                masked_peak_lags = peak_lags.copy()
                below_threshold = (
                    ~np.isfinite(peak_magnitudes)
                    | (peak_magnitudes < float(peak_latency_threshold))
                )
                masked_peak_lags[below_threshold] = np.nan
            label = _line_label_for_surface(surface)
            mag_ax.plot(
                surface.window_coords_sec,
                peak_magnitudes,
                linewidth=1.6,
                label=label,
            )
            lag_ax.plot(
                surface.window_coords_sec,
                masked_peak_lags,
                linewidth=1.6,
                label=label,
            )

        for axis in (mag_ax, lag_ax):
            _add_occlusion_markers(axis, occlusion_timing=occlusion_timing)
            axis.grid(True, color="0.9", linewidth=0.8)
            axis.axhline(0.0, color="0.55", linestyle="--", linewidth=0.8)

        if peak_latency_threshold is not None:
            mag_ax.axhline(
                float(peak_latency_threshold),
                color="0.2",
                linestyle="--",
                linewidth=0.9,
                label="latency mask threshold",
            )
        mag_ax.set_ylabel("Peak dRSA magnitude")
        lag_ax.set_ylabel("Peak lag [s]")
        lag_ax.set_xlabel(f"Model-time window {args.window_coordinate} [s]")
        mag_ax.set_title(
            _compact_plot_title(
                [
                    selected[0].record.owner_label,
                    comparison["title"],
                    selected[0].record.sensor_group,
                    selected[0].record.analysis_family,
                    selected[0].record.variable,
                ]
            ),
            fontsize=10,
        )
        mag_ax.legend(loc="best", fontsize=8, frameon=True)
        fig.savefig(write_target, dpi=180, facecolor="white")
        plt.close(fig)

    return {
        "plot_type": "sliding_window_peak_timecourse",
        "figure": str(out_path),
        "owner": selected[0].record.owner_label,
        "comparison": str(comparison["slug"]),
        "targets": [surface.target.slug for surface in selected],
        "analysis_family": selected[0].record.analysis_family,
        "sensor_group": selected[0].record.sensor_group,
        "variable": selected[0].record.variable,
        "caption": (
            f"{selected[0].record.owner_label} | {comparison['title']} | "
            f"{selected[0].record.analysis_family} | {selected[0].record.sensor_group} | "
            f"{selected[0].record.variable} | top panel=max dRSA over lag per "
            f"model-time window; bottom panel=lag at that max. "
            f"Latency mask mode={args.peak_latency_mask_mode}; "
            f"threshold="
            f"{'none' if peak_latency_threshold is None else f'{peak_latency_threshold:.6g}'}."
        ),
    }


def _risk_notes() -> list[str]:
    """Return concise report-facing interpretation risks."""
    return [
        (
            "Overlapping windows share model-time rows, so neighboring rows in "
            "the heatmap/waterfall are not independent observations."
        ),
        (
            "Temporal autocorrelation in trajectories and neural patterns can "
            "produce smooth ridges even without a sharply timed neural delay."
        ),
        (
            "Large positive or negative lags have fewer available diagonal "
            "samples near matrix edges; cells below the configured finite-sample "
            "threshold are hidden."
        ),
        (
            "Model-anchored support restricts only model-time rows. Neural-time "
            "columns may fall before or after the displayed model window."
        ),
        (
            "Peak-lag traces are descriptive argmax overlays and should not be "
            "interpreted as inferential timing estimates without a later "
            "subject-level statistical procedure."
        ),
    ]


def _write_report(
    report_path: Path,
    *,
    title: str,
    summary_lines: Iterable[str],
    plot_entries: Sequence[Mapping[str, Any]],
    skipped_inputs: Sequence[Mapping[str, str]],
    raw_root: Path,
    plotly_js_url: str,
) -> Path:
    """
    Write the Stage-09 HTML report and embedded sibling copy.

    Heatmaps and peak timecourses are embedded as images. Waterfalls are
    embedded directly as Plotly divs and also linked to their standalone HTML
    files.
    """
    report = guarded_path(report_path, raw_root)
    report.parent.mkdir(parents=True, exist_ok=True)

    toc = "\n".join(
        f"<li><a href=\"#{html_id(plot_type)}\">{html.escape(plot_type)}</a></li>"
        for plot_type in PLOT_TYPES
    )
    blocks: list[str] = []
    for plot_type in PLOT_TYPES:
        type_entries = [item for item in plot_entries if item.get("plot_type") == plot_type]
        figures: list[str] = []
        for item in type_entries:
            caption = html.escape(
                " | ".join([str(item.get("caption", "")), f"plot_type={plot_type}"])
            )
            if plot_type in {
                "lag_by_model_time_heatmap",
                "sliding_window_peak_timecourse",
            }:
                href = image_href(report, str(item["figure"]))
                figures.append(
                    "<figure>"
                    f"<img src=\"{html.escape(href)}\"/>"
                    f"<figcaption>{caption}</figcaption>"
                    "</figure>"
                )
            else:
                linked = image_href(report, str(item["html"]))
                figures.append(
                    "<figure>"
                    + _plotly_div_html(item)
                    + (
                        f"<figcaption>{caption} | "
                        f"<a href=\"{html.escape(linked)}\">standalone interactive HTML</a>"
                        "</figcaption>"
                    )
                    + "</figure>"
                )
        if not figures:
            figures.append("<p>No figures for this plot type.</p>")
        blocks.append(
            f"<section id=\"{html_id(plot_type)}\"><h2>{html.escape(plot_type)}</h2>"
            + "\n".join(figures)
            + "</section>"
        )

    skipped_block = ""
    if skipped_inputs:
        rows = "\n".join(
            "<li>"
            + html.escape(
                " | ".join(f"{key}: {value}" for key, value in sorted(item.items()))
            )
            + "</li>"
            for item in skipped_inputs
        )
        skipped_block = f"<details><summary>Skipped inputs</summary><ul>{rows}</ul></details>"

    risk_rows = "\n".join(f"<li>{html.escape(note)}</li>" for note in _risk_notes())
    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            f"<script src=\"{html.escape(plotly_js_url)}\"></script>",
            "<style>",
            "body{font-family:system-ui,sans-serif;max-width:1680px;margin:2rem auto;padding:0 1rem;}",
            "nav,pre{background:#f5f5f5;padding:1rem;border-radius:6px;}",
            "figure{margin:1.5rem 0;padding-bottom:1rem;border-bottom:1px solid #ddd;}",
            "img{max-width:100%;height:auto;border:1px solid #ddd;}",
            ".plotly-chart{width:100%;height:680px;border:1px solid #ddd;}",
            "figcaption{font-size:.9rem;color:#333;margin-top:.4rem;}",
            "h2{margin-top:3rem;border-top:3px solid #222;padding-top:1rem;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<pre>{html.escape(chr(10).join(summary_lines))}</pre>",
            "<nav><strong>Plot types</strong><ul>",
            toc,
            "</ul></nav>",
            "<section><h2>Interpretation risks</h2><ul>",
            risk_rows,
            "</ul></section>",
            skipped_block,
            "\n".join(blocks),
            "</body></html>",
        ]
    )
    report.write_text(body, encoding="utf-8")
    write_embedded_html_copy(report, raw_root=raw_root)
    return report


def _plot_records(
    records: Sequence[MatrixRecord],
    *,
    output_dir: Path,
    report_html: Path,
    title: str,
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
    skipped_inputs: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    """Generate all Stage-09 plots and the HTML report for one owner scope."""
    raw_root = Path(cfg["resolved_paths"]["data_root"])
    output_dir = guarded_path(output_dir, raw_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    occlusion_timing = occlusion_timing_from_config(cfg)

    selected, target_skipped = _target_records(records)
    skipped: list[Mapping[str, str]] = [*skipped_inputs, *target_skipped]
    surfaces = [
        _surface_for_record(record, spec, args=args)
        for spec, record in selected
    ]
    value_limits = _value_limits(
        surfaces,
        mode=args.value_limit_mode,
        vmin=args.vmin,
        vmax=args.vmax,
    )
    peak_latency_threshold = _peak_latency_threshold(surfaces, args=args)

    plot_entries: list[dict[str, Any]] = []
    for surface in surfaces:
        plot_entries.append(
            _plot_heatmap(
                surface,
                output_dir=output_dir,
                raw_root=raw_root,
                occlusion_timing=occlusion_timing,
                value_limits=value_limits,
                args=args,
            )
        )
    for surface in surfaces:
        plot_entries.append(
            _write_waterfall(
                surface,
                output_dir=output_dir,
                raw_root=raw_root,
                value_limits=value_limits,
                args=args,
            )
        )
    for comparison in PEAK_TIMECOURSE_COMPARISONS:
        entry = _plot_peak_timecourse_comparison(
            surfaces,
            comparison=comparison,
            output_dir=output_dir,
            raw_root=raw_root,
            occlusion_timing=occlusion_timing,
            peak_latency_threshold=peak_latency_threshold,
            args=args,
        )
        if entry is not None:
            plot_entries.append(entry)

    first_sec = occlusion_timing["first_complete_occlusion_sec"]
    reapp_sec = occlusion_timing["reappearance_sec"]
    summary_lines = [
        "Stage status: completed",
        f"Generated UTC: {utc_now()}",
        "Lag convention: lag = neural time - model time",
        "Matrix interpretation: x = neural time, y = model time, indexed matrix[model_time, neural_time]",
        "Sliding-window support: model-anchored; windows restrict model-time rows only",
        (
            "Occlusion boundaries: "
            f"first_complete_occlusion frame {occlusion_timing['first_complete_occlusion_frame']} "
            f"= {first_sec:.6f} s; "
            f"reappearance frame {occlusion_timing['reappearance_frame']} "
            f"= {reapp_sec:.6f} s"
        ),
        f"Sensor group: {normalize_sensor_group(args.sensor_group)}",
        f"Analysis family: {args.analysis_family}",
        f"dRSA type: {args.drsa_type}",
        f"Variable: {args.variable}",
        f"Lag window: +/-{args.lag_window_sec:g} s",
        f"Model-time window: {args.model_window_sec:g} s",
        f"Model-time step: {args.model_step_sec:g} s",
        f"Window coordinate: {args.window_coordinate}",
        f"Minimum finite samples per cell: {args.min_valid_samples}",
        f"Shared value limits: {value_limits[0]:.6g} to {value_limits[1]:.6g}",
        f"Value-limit mode: {args.value_limit_mode}",
        f"Peak-latency mask mode: {args.peak_latency_mask_mode}",
        (
            "Peak-latency mask threshold: "
            + (
                "none"
                if peak_latency_threshold is None
                else f"{peak_latency_threshold:.6g}"
            )
        ),
        f"Target summaries found: {len(selected)} of {len(TARGET_SUMMARIES)}",
        f"Matrix records scanned: {len(records)}",
        f"Figures/interactive plots: {len(plot_entries)}",
        f"Skipped inputs or target notes: {len(skipped)}",
        f"Output directory: {output_dir}",
    ]

    _write_report(
        report_html,
        title=title,
        summary_lines=summary_lines,
        plot_entries=plot_entries,
        skipped_inputs=skipped,
        raw_root=raw_root,
        plotly_js_url=args.plotly_js_url,
    )
    return {
        "status": "completed",
        "record_count": len(records),
        "target_count": len(selected),
        "figure_count": len(plot_entries),
        "value_limits": value_limits,
        "output_dir": str(output_dir),
        "report_html": str(report_html),
        "report_embedded_html": str(report_html.with_name(f"{report_html.stem}_embedded.html")),
        "plot_entries": plot_entries,
        "skipped_inputs": list(skipped),
    }


def run_subject_scope(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Build Stage-09 outputs for one resolved subject configuration."""
    owner_label = str(cfg["layout"]["subject_bids_id"])
    paths = _analysis_paths(cfg)
    records, skipped = _read_subject_records(
        cfg,
        owner_label=owner_label,
        sensor_group=args.sensor_group,
        analysis_family=args.analysis_family,
        drsa_type=args.drsa_type,
        variable=args.variable,
    )
    if not records:
        raise RuntimeError(f"No matching Stage-07 result matrices were found for {owner_label}.")
    return _plot_records(
        records,
        output_dir=paths["subject_surface_dir"],
        report_html=paths["subject_report_html"],
        title=f"{owner_label} | 09 dRSA sliding lag surfaces",
        cfg=cfg,
        args=args,
        skipped_inputs=skipped,
    )


def run_group_scope(
    base_cfg: Mapping[str, Any],
    args: argparse.Namespace,
    subjects: Sequence[str],
) -> dict[str, Any]:
    """Build group-average Stage-09 outputs from existing subject MAT files."""
    subject_records: list[MatrixRecord] = []
    skipped: list[dict[str, str]] = []
    for subject in subjects:
        subject_cfg = _load_cfg(
            args,
            subject_config=_subject_config_for_id(args, subject),
        )
        records, subject_skipped = _read_subject_records(
            subject_cfg,
            owner_label=subject,
            sensor_group=args.sensor_group,
            analysis_family=args.analysis_family,
            drsa_type=args.drsa_type,
            variable=args.variable,
        )
        subject_records.extend(records)
        skipped.extend(subject_skipped)

    records, group_skipped = _group_records(subject_records)
    skipped.extend(group_skipped)
    if not records:
        raise RuntimeError(
            "No group-average sliding lag records were produced. "
            f"Skipped inputs: {skipped[:5]}"
        )
    paths = _analysis_paths(base_cfg)
    return _plot_records(
        records,
        output_dir=paths["group_surface_dir"],
        report_html=paths["group_report_html"],
        title=f"Group | {base_cfg['layout']['session_bids_id']} | 09 dRSA sliding lag surfaces",
        cfg=base_cfg,
        args=args,
        skipped_inputs=skipped,
    )


def main() -> int:
    """Resolve configuration, generate requested Stage-09 summaries, and exit."""
    args = build_arg_parser().parse_args()
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if args.group_only and args.subjects_only:
        raise SystemExit("--group-only and --subjects-only are mutually exclusive.")
    if args.lag_window_sec <= 0:
        raise SystemExit("--lag-window-sec must be > 0.")
    if args.model_window_sec <= 0:
        raise SystemExit("--model-window-sec must be > 0.")
    if args.model_step_sec <= 0:
        raise SystemExit("--model-step-sec must be > 0.")
    if args.min_valid_samples < 1:
        raise SystemExit("--min-valid-samples must be >= 1.")
    if args.value_limit_mode == "manual" and (args.vmin is None or args.vmax is None):
        raise SystemExit("--value-limit-mode manual requires --vmin and --vmax.")
    if args.value_limit_mode != "manual" and (args.vmin is not None or args.vmax is not None):
        raise SystemExit("--vmin/--vmax require --value-limit-mode manual.")
    if args.peak_latency_min_magnitude_frac < 0:
        raise SystemExit("--peak-latency-min-magnitude-frac must be >= 0.")
    if (
        args.peak_latency_mask_mode == "absolute"
        and args.peak_latency_min_magnitude is None
    ):
        raise SystemExit(
            "--peak-latency-mask-mode absolute requires "
            "--peak-latency-min-magnitude."
        )

    base_subject_config = args.subject_config
    if args.subjects:
        base_subject_config = _subject_config_for_id(args, normalize_subject_label(args.subjects[0]))
    base_cfg = _load_cfg(args, subject_config=base_subject_config)
    subjects = _subjects_from_args_or_config(args, base_cfg)

    reports: dict[str, Any] = {
        "stage": "drsa_sliding_lag_surfaces",
        "generated_utc": utc_now(),
        "subjects": subjects,
        "sensor_group": normalize_sensor_group(args.sensor_group),
        "analysis_family": args.analysis_family,
        "drsa_type": args.drsa_type,
        "variable": args.variable,
        "reports": {},
    }
    try:
        if not args.group_only:
            subject_reports = {}
            for subject in subjects if args.subjects else [base_cfg["layout"]["subject_bids_id"]]:
                subject_cfg = _load_cfg(
                    args,
                    subject_config=_subject_config_for_id(args, subject),
                )
                subject_reports[str(subject)] = run_subject_scope(subject_cfg, args)
            reports["reports"]["subjects"] = subject_reports

        if not args.subjects_only:
            reports["reports"]["group"] = run_group_scope(base_cfg, args, subjects)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"dRSA sliding lag-surface error: {exc}")
        return 2

    if args.print_summary:
        printable = {
            "stage": reports["stage"],
            "subjects": reports["subjects"],
            "sensor_group": reports["sensor_group"],
            "analysis_family": reports["analysis_family"],
            "drsa_type": reports["drsa_type"],
            "variable": reports["variable"],
            "reports": {
                scope: (
                    {
                        key: {
                            "record_count": value["record_count"],
                            "target_count": value["target_count"],
                            "figure_count": value["figure_count"],
                            "value_limits": value["value_limits"],
                            "output_dir": value["output_dir"],
                            "report_html": value["report_html"],
                            "report_embedded_html": value["report_embedded_html"],
                        }
                        for key, value in report.items()
                    }
                    if scope == "subjects"
                    else {
                        "record_count": report["record_count"],
                        "target_count": report["target_count"],
                        "figure_count": report["figure_count"],
                        "value_limits": report["value_limits"],
                        "output_dir": report["output_dir"],
                        "report_html": report["report_html"],
                        "report_embedded_html": report["report_embedded_html"],
                    }
                )
                for scope, report in reports["reports"].items()
            },
        }
        print(json.dumps(printable, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
