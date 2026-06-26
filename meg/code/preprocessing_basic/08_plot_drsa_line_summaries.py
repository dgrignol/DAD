#!/usr/bin/env python3
"""
08_plot_drsa_line_summaries.py

Purpose
Generate optional exploratory line-plot summaries from existing Stage-07 dRSA
result MAT files. This stage does not run MATLAB and does not recompute dRSA.
It consumes the saved Stage-07 matrices and writes lightweight interval lag
summaries under the configured derivatives tree:

Interval lag curves are computed by diagonal averaging inside occlusion
intervals: pre-occlusion, during occlusion, and post-occlusion. Every figure
carries an explicit stimulus/comparison model label, such as "actual deviant
position" or "expected/predicted position".

Usage examples
1) Build group-average line summaries from existing Stage-07 result MAT files:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 \
     --group-only --overwrite --print-summary

2) Build only smoothed-display group summaries for the configured subjects:

   python .../08_plot_drsa_line_summaries.py --config .../pipeline_default.yaml \
     --group-only --include-smoothed --overwrite --print-summary

3) Build subject-level summaries for one participant without touching MATLAB:

   python .../08_plot_drsa_line_summaries.py --config .../pipeline_default.yaml \
     --subject-config .../sub-04.yaml --subjects-only --resume --print-summary

4) Restrict to occipital sensor-group Stage-07 outputs and a +/-0.5 s lag range:

   python .../08_plot_drsa_line_summaries.py --config .../pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 --group-only \
     --sensor-groups occipital --lag-window-sec 0.5 --overwrite

5) Add descriptive within-subject diagonal-SD shading to group interval lag
   curves:

   python .../08_plot_drsa_line_summaries.py --config .../pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 --group-only \
     --shading within-subject-diagonal-sd --overwrite --print-summary

6) Use strict same-interval square support for interval lag curves, so both the
   model-time row and neural-time column must fall inside the named interval:

   python .../08_plot_drsa_line_summaries.py --config .../pipeline_default.yaml \
     --subjects 01 02 03 04 05 06 07 --group-only \
     --interval-lag-mode same-interval-square --overwrite --print-summary

Inputs
- Existing Stage-07 result MAT files below
  derivatives/analysis/sub-XX/ses-01/drsa/.
- pipeline_config.py-resolved derivative paths and matlab_drsa.occlusion_timing.

Outputs
- Subject figures:
  derivatives/analysis/sub-XX/ses-01/drsa/line_summaries/
  derivatives/reports/sub-XX/ses-01/08_drsa_line_summaries_report.html
  derivatives/reports/sub-XX/ses-01/08_drsa_line_summaries_report_embedded.html
- Group figures:
  derivatives/analysis/group/ses-01/drsa/line_summaries/
  derivatives/reports/group/ses-01/08_drsa_line_summaries_report.html
  derivatives/reports/group/ses-01/08_drsa_line_summaries_report_embedded.html

Safety and assumptions
- Raw data are never modified. Every output path passes the shared raw-tree
  guard before writing.
- Stage-07 matrix PNGs stay where they are. This script writes only to the new
  line_summaries folder and Stage-08 report filenames.
- The lag convention is preserved as neural time minus model time. Stage-07
  plots x = neural time and y = model time, so this script indexes each matrix
  as matrix[model_time_index, neural_time_index].
- ``--include-smoothed`` smooths matrices in Python only for display summaries;
  it does not require or create smoothed MAT files, and originally undefined
  cells remain NaN so display smoothing cannot create values outside a model's
  valid analysis support.
- ``--shading`` is optional and applies only to interval lag-curve plots. The
  default ``none`` preserves the original unshaded behavior and filenames.
- ``--interval-lag-mode`` controls the support for interval lag curves. The
  default ``model-anchored`` preserves the previous behavior: the interval mask
  restricts model-time rows, while neural-time columns may fall anywhere inside
  the full matrix. ``same-interval-square`` additionally requires the neural-time
  column to fall inside the same interval.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
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
    html_escape_join,
    html_id,
    image_href,
    iter_matrix_slices,
    model_source_for_matrix_slice,
    model_time_axis_label,
    model_source_from_path,
    normalize_sensor_group,
    normalize_subject_label,
    occlusion_intervals,
    occlusion_timing_from_config,
    read_drsa_result_file,
    sampling_step_seconds,
    sanitize_label,
    sensor_group_from_path,
    smooth_matrix_for_plot,
    stimulus_model_label,
    time_axis_for_matrix,
    utc_now,
    write_embedded_html_copy,
)
from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import guarded_path


DEFAULT_SAMPLING_RATE_HZ = 120.0
PLOT_TYPES = ("interval_lag_curve",)
SHADING_MODES = (
    "none",
    "within-subject-diagonal-sd",
    "across-subject-sem",
)
INTERVAL_LAG_MODES = (
    "model-anchored",
    "same-interval-square",
)


@dataclass(frozen=True)
class MatrixRecord:
    """
    One 2-D dRSA matrix slice plus its analysis metadata.

    Data flow:
    Stage-07 result MAT variable -> 2-D matrix slice -> lag summary plots.
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
    variant: str
    matrix: np.ndarray
    time_axis_sec: np.ndarray
    relative_result: str
    source_paths: tuple[str, ...]
    subjects: tuple[str, ...] = field(default_factory=tuple)
    n_subjects: int = 1
    subject_records: tuple["MatrixRecord", ...] = field(default_factory=tuple)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the Stage-08 command-line interface."""
    parser = argparse.ArgumentParser(
        description="Generate optional line summaries from existing dRSA MAT files."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--subjects", nargs="+", default=None)
    parser.add_argument("--group-only", action="store_true")
    parser.add_argument("--subjects-only", action="store_true")
    parser.add_argument(
        "--conditions",
        nargs="+",
        choices=CONDITIONS,
        default=None,
    )
    parser.add_argument(
        "--sensor-groups",
        nargs="+",
        default=None,
        help=(
            "Optional sensor-group filter. Defaults to all existing Stage-07 "
            "result MAT files."
        ),
    )
    parser.add_argument(
        "--lag-window-sec",
        type=float,
        default=1.0,
        help="Maximum absolute lag, in seconds, for lag searches.",
    )
    parser.add_argument(
        "--include-smoothed",
        action="store_true",
        help="Also compute display-only smoothed matrix summaries in Python.",
    )
    parser.add_argument(
        "--shading",
        choices=SHADING_MODES,
        default="none",
        help=(
            "Optional interval lag-curve shaded band: none, descriptive "
            "within-subject diagonal SD, or across-subject SEM."
        ),
    )
    parser.add_argument(
        "--interval-lag-mode",
        choices=INTERVAL_LAG_MODES,
        default="model-anchored",
        help=(
            "Interval lag-curve support: model-anchored restricts model-time "
            "rows only; same-interval-square also requires neural-time columns "
            "to be inside the same interval."
        ),
    )
    parser.add_argument(
        "--plot-smoothing-window",
        type=int,
        default=None,
        help="Gaussian smoothing window, in samples, for --include-smoothed.",
    )
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


def _analysis_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return subject and group Stage-08 output paths from resolved config."""
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
        "subject_line_dir": analysis_dir / "drsa" / "line_summaries",
        "subject_report_html": reports_dir / "08_drsa_line_summaries_report.html",
        "group_analysis_dir": group_analysis_dir,
        "group_reports_dir": group_reports_dir,
        "group_drsa_dir": group_analysis_dir / "drsa",
        "group_line_dir": group_analysis_dir / "drsa" / "line_summaries",
        "group_report_html": group_reports_dir / "08_drsa_line_summaries_report.html",
    }


def _load_cfg(
    args: argparse.Namespace,
    *,
    subject_config: str | Path | None = None,
) -> Mapping[str, Any]:
    """Load one resolved pipeline configuration for Stage-08 use."""
    cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("export", "matlab_drsa"),
    )
    if warnings and args.print_summary:
        print(json.dumps({"configuration_warnings": warnings}, indent=2))
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
    """Resolve the sample rate used when a MAT file lacks time_axis_sec."""
    return float(
        cfg.get("export", {}).get(
            "sampling_rate_hz",
            cfg.get("matlab_drsa", {}).get(
                "plot_sample_rate_hz",
                DEFAULT_SAMPLING_RATE_HZ,
            ),
        )
    )


def _selected_conditions(args: argparse.Namespace) -> set[str] | None:
    """Return the optional condition filter as a set."""
    return set(args.conditions) if args.conditions else None


def _selected_sensor_groups(args: argparse.Namespace) -> set[str] | None:
    """Return the optional sensor-group filter as normalized tags."""
    if not args.sensor_groups:
        return None
    return {normalize_sensor_group(value) for value in args.sensor_groups}


def _should_skip_result_path(path: Path) -> bool:
    """Filter out non-result MAT files and Stage-08 output folders."""
    name = path.name
    text = str(path)
    return (
        "desc-identityAveraged_neural" in name
        or "/line_summaries/" in text
        or "line_summaries" in path.parts
    )


def _read_subject_records(
    cfg: Mapping[str, Any],
    *,
    owner_label: str,
    include_smoothed: bool,
    smoothing_window_samples: int,
    conditions: set[str] | None,
    sensor_groups: set[str] | None,
) -> tuple[list[MatrixRecord], list[dict[str, str]]]:
    """
    Discover and read Stage-07 result matrices for one subject.

    Data flow:
    subject drsa folder -> result MAT files -> matrix variables -> per-slice
    MatrixRecord objects. Display-only smoothed variants are produced in memory
    when requested.
    """
    paths = _analysis_paths(cfg)
    sampling_rate = _sampling_rate_for_plots(cfg)
    drsa_dir = paths["drsa_dir"]
    records: list[MatrixRecord] = []
    skipped: list[dict[str, str]] = []
    if not drsa_dir.is_dir():
        skipped.append({"owner": owner_label, "reason": f"missing dRSA folder: {drsa_dir}"})
        return records, skipped

    for mat_path in sorted(drsa_dir.rglob("*.mat")):
        if _should_skip_result_path(mat_path):
            continue
        relative = mat_path.relative_to(drsa_dir)
        condition = condition_from_path(relative)
        sensor_group = sensor_group_from_path(relative)
        analysis_family = analysis_family_from_path(relative)
        path_model_source = model_source_from_path(relative)
        if conditions is not None and condition not in conditions:
            continue
        if sensor_groups is not None and sensor_group not in sensor_groups:
            continue
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
        if not matrices:
            skipped.append(
                {
                    "owner": owner_label,
                    "path": str(mat_path),
                    "reason": "no supported dRSA matrix variables found",
                }
            )
            continue

        for variable, matrix_array in matrices.items():
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
                    analysis_family=analysis_family,
                )
                model_label = stimulus_model_label(
                    condition,
                    path_model_source,
                    variable,
                    slice_label=slice_label,
                    analysis_family=analysis_family,
                )
                axis_label = model_time_axis_label(
                    condition,
                    path_model_source,
                    variable,
                    slice_label=slice_label,
                    analysis_family=analysis_family,
                )
                axis = time_axis_for_matrix(
                    matrix_2d,
                    time_axis=time_axis,
                    sampling_rate_hz=sampling_rate,
                )
                common = {
                    "owner_scope": "subject",
                    "owner_label": owner_label,
                    "condition": condition,
                    "analysis_family": analysis_family,
                    "sensor_group": sensor_group,
                    "drsa_type": drsa_type_from_path(relative),
                    "model_source": effective_model_source,
                    "stimulus_model_label": model_label,
                    "model_time_axis_label": axis_label,
                    "variable": variable,
                    "slice_index": int(slice_index),
                    "slice_label": slice_label,
                    "time_axis_sec": np.asarray(axis, dtype=float),
                    "relative_result": str(relative),
                    "source_paths": (str(mat_path),),
                    "subjects": (owner_label,),
                    "n_subjects": 1,
                }
                records.append(
                    MatrixRecord(
                        **common,
                        variant="raw",
                        matrix=np.asarray(matrix_2d, dtype=float),
                    )
                )
                if include_smoothed:
                    records.append(
                        MatrixRecord(
                            **common,
                            variant="smoothed_display",
                            matrix=smooth_matrix_for_plot(
                                matrix_2d,
                                smoothing_window_samples,
                            ),
                        )
                    )
    return records, skipped


def _group_records(
    subject_records: Sequence[MatrixRecord],
) -> tuple[list[MatrixRecord], list[dict[str, str]]]:
    """
    Average matching matrix records across subjects.

    Records match when they share the same Stage-07 relative result path,
    variable, slice, and raw/smoothed-display variant. Shape mismatches are
    skipped with explicit report metadata.
    """
    grouped: dict[tuple[str, str, int, str], list[MatrixRecord]] = {}
    for record in subject_records:
        key = (
            record.relative_result,
            record.variable,
            record.slice_index,
            record.variant,
        )
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
                variant=first.variant,
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


def _lag_samples(time_axis_sec: np.ndarray, lag_window_sec: float) -> np.ndarray:
    """Build integer lag offsets for the requested lag window."""
    dt = sampling_step_seconds(time_axis_sec, DEFAULT_SAMPLING_RATE_HZ)
    max_lag = int(round(float(lag_window_sec) / dt))
    return np.arange(-max_lag, max_lag + 1, dtype=int)


def _nanmean_or_nan(values: Sequence[float]) -> float:
    """Return nanmean while preserving NaN for empty/all-NaN samples."""
    if not values:
        return float("nan")
    array = np.asarray(values, dtype=float)
    finite = np.isfinite(array)
    if not np.any(finite):
        return float("nan")
    return float(np.nanmean(array))


def _nanstd_or_nan(values: Sequence[float]) -> float:
    """Return nanstd while preserving NaN for empty/all-NaN samples."""
    if not values:
        return float("nan")
    array = np.asarray(values, dtype=float)
    finite = np.isfinite(array)
    if not np.any(finite):
        return float("nan")
    return float(np.nanstd(array))


def _nanmean_axis0(values: Any) -> np.ndarray:
    """Return a NaN-safe column mean for stacked lag curves or bands."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"Expected a 2-D array, found shape {array.shape}.")
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(array, axis=0)


def compute_across_subject_sem(subject_curves: Any) -> np.ndarray:
    """
    Compute NaN-safe SEM across subject-level lag curves at each lag.

    The SEM uses the sample standard deviation (ddof=1) divided by the square
    root of the finite subject count at each lag. Lags with fewer than two
    finite subject values return NaN because between-subject variability cannot
    be estimated there.
    """
    curves = np.asarray(subject_curves, dtype=float)
    if curves.ndim != 2:
        raise ValueError(f"Expected a 2-D subject x lag array, found {curves.shape}.")
    finite_counts = np.sum(np.isfinite(curves), axis=0)
    sem = np.full(curves.shape[1], np.nan, dtype=float)
    estimable = finite_counts > 1
    if np.any(estimable):
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            std = np.nanstd(curves[:, estimable], axis=0, ddof=1)
        sem[estimable] = std / np.sqrt(finite_counts[estimable])
    return sem


def compute_interval_lag_curve_with_diagonal_sd(
    matrix_2d: Any,
    lag_samples: Sequence[int],
    model_time_mask: Sequence[bool],
    *,
    interval_lag_mode: str = "model-anchored",
    neural_time_mask: Sequence[bool] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute interval lag means and within-diagonal SD values.

    Lag convention
    ``lag = neural time - model time``. Stage 07 plots x = neural time and
    y = model time, so a sample at lag ``L`` is read as
    ``matrix[model_idx, model_idx + L]``.

    Interval modes
    ``model-anchored`` preserves the original Stage-08 behavior by applying the
    interval mask only to model-time rows. ``same-interval-square`` also applies
    the interval mask to neural-time columns, preventing the named interval from
    averaging diagonals that leave the interval on the neural-time axis.

    Returns
    -------
    lag_curve
        NaN-safe mean of the diagonal samples contributing to each lag.
    diagonal_sd
        NaN-safe standard deviation across the same contributing diagonal
        samples for each lag. This is descriptive within-matrix variability,
        not a confidence interval.
    """
    if interval_lag_mode not in INTERVAL_LAG_MODES:
        raise ValueError(
            f"Unknown interval_lag_mode {interval_lag_mode!r}; "
            f"expected one of {INTERVAL_LAG_MODES}."
        )
    matrix = np.asarray(matrix_2d, dtype=float)
    mask = np.asarray(model_time_mask, dtype=bool)
    n_model, n_neural = matrix.shape
    model_indices = np.flatnonzero(mask[:n_model])
    neural_mask = None
    if interval_lag_mode == "same-interval-square":
        neural_mask = (
            np.asarray(neural_time_mask, dtype=bool)
            if neural_time_mask is not None
            else mask
        )
    lag_values: list[float] = []
    lag_sds: list[float] = []
    for lag in lag_samples:
        samples: list[float] = []
        for model_idx in model_indices:
            neural_idx = int(model_idx + lag)
            if 0 <= neural_idx < n_neural:
                if neural_mask is not None:
                    if neural_idx >= len(neural_mask) or not bool(neural_mask[neural_idx]):
                        continue
                samples.append(float(matrix[model_idx, neural_idx]))
        lag_values.append(_nanmean_or_nan(samples))
        lag_sds.append(_nanstd_or_nan(samples))
    return np.asarray(lag_values, dtype=float), np.asarray(lag_sds, dtype=float)


def compute_interval_lag_curve(
    matrix_2d: Any,
    lag_samples: Sequence[int],
    model_time_mask: Sequence[bool],
    *,
    interval_lag_mode: str = "model-anchored",
    neural_time_mask: Sequence[bool] | None = None,
) -> np.ndarray:
    """
    Compute diagonal-average lag values inside one interval support.

    This compatibility wrapper preserves the original Stage-08 public helper
    while the shaded plot path also consumes the matching diagonal SD values.
    The default ``model-anchored`` mode applies the interval only to model-time
    rows; ``same-interval-square`` also requires neural-time columns to be inside
    the same interval.
    """
    curve, _diagonal_sd = compute_interval_lag_curve_with_diagonal_sd(
        matrix_2d,
        lag_samples,
        model_time_mask,
        interval_lag_mode=interval_lag_mode,
        neural_time_mask=neural_time_mask,
    )
    return curve


def _output_slug_for_record(record: MatrixRecord) -> str:
    """Build a compact output slug from one matrix record."""
    parts = [
        record.owner_label,
        record.variant,
        record.analysis_family,
        record.condition,
        record.sensor_group,
        record.drsa_type,
        record.model_source,
        record.stimulus_model_label,
        record.variable,
        record.slice_label,
    ]
    return "_".join(sanitize_label(part) for part in parts if part)


def _overlay_key(record: MatrixRecord) -> tuple[str, ...]:
    """Group records that should share interval lag overlay plots."""
    return (
        record.owner_scope,
        record.owner_label,
        record.variant,
        record.analysis_family,
        record.sensor_group,
        record.drsa_type,
        record.variable,
        record.slice_label,
    )


def _overlay_slug(records: Sequence[MatrixRecord], interval_name: str) -> str:
    """Build an output slug for one interval overlay plot."""
    first = records[0]
    labels = {record.stimulus_model_label for record in records}
    model_label_slug = (
        first.stimulus_model_label if len(labels) == 1 else "stimulus-model-comparison"
    )
    parts = [
        first.owner_label,
        first.variant,
        first.analysis_family,
        first.sensor_group,
        first.drsa_type,
        model_label_slug,
        first.variable,
        first.slice_label,
        interval_name,
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


def _shading_interpretation(shading_mode: str) -> str:
    """Return a report-ready interpretation of one interval shading mode."""
    if shading_mode == "within-subject-diagonal-sd":
        return (
            "within-subject-diagonal-sd: shaded band is +/- the mean "
            "subject-level SD across diagonal samples contributing to each lag; "
            "it is descriptive within-subject diagonal variability, not a "
            "statistical confidence interval."
        )
    if shading_mode == "across-subject-sem":
        return (
            "across-subject-sem: group shaded band is +/- the SEM across "
            "subject-level interval lag curves at each lag."
        )
    return "none: no shaded band; interval lag curves use the original unshaded behavior."


def _shading_legend_label(shading_mode: str) -> str:
    """Return a compact legend label for the shaded band."""
    if shading_mode == "within-subject-diagonal-sd":
        return "descriptive within-subject diagonal SD"
    if shading_mode == "across-subject-sem":
        return "across-subject SEM"
    return "no shading"


def _compact_plot_title(parts: Iterable[Any], *, width: int = 86) -> str:
    """
    Return a wrapped, compact Matplotlib title from plot-identifying fields.

    Data flow:
    metadata fields that identify the visual comparison -> short title inside
    the PNG. Complete condition/model/shading details remain in HTML captions.
    """
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


def _legend_label_for_record(record: MatrixRecord) -> str:
    """Return the curve legend label without duplicating band metadata."""
    return f"{record.condition} | {record.stimulus_model_label}"


def _wrapped_legend_label(label: str, *, width: int = 34) -> str:
    """Wrap long legend labels so the external legend panel stays readable."""
    return "\n".join(
        textwrap.wrap(
            str(label),
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    ) or str(label)


def _interval_lag_mode_interpretation(interval_lag_mode: str) -> str:
    """Return a report-ready description of interval lag support."""
    if interval_lag_mode == "same-interval-square":
        return (
            "same-interval-square: both model-time rows and neural-time columns "
            "must fall inside the named interval; lags outside the interval "
            "support are NaN."
        )
    return (
        "model-anchored: the named interval restricts model-time rows only; "
        "neural-time columns may fall anywhere inside the full matrix."
    )


def _interval_lag_summary_for_record(
    record: MatrixRecord,
    *,
    lag_samples: Sequence[int],
    interval_mask: Sequence[bool],
    interval_name: str,
    occlusion_timing: Mapping[str, Any],
    shading_mode: str,
    interval_lag_mode: str,
) -> tuple[np.ndarray, np.ndarray | None, str | None]:
    """
    Compute the line and optional shaded band for one interval lag curve.

    Data flow:
    subject matrices -> subject lag curves/diagonal SDs -> optional group mean
    line and group band. In ``none`` mode the group path intentionally preserves
    the original behavior by computing the curve from the averaged group matrix.
    """
    if shading_mode == "none":
        curve = compute_interval_lag_curve(
            record.matrix,
            lag_samples,
            interval_mask,
            interval_lag_mode=interval_lag_mode,
        )
        return curve, None, None

    if record.owner_scope == "group":
        if not record.subject_records:
            note = (
                f"{shading_mode} requested for {record.owner_label}, but "
                "subject-level matrices were unavailable; plotted the available "
                "group matrix without a shaded band."
            )
            warnings.warn(note, RuntimeWarning)
            curve = compute_interval_lag_curve(
                record.matrix,
                lag_samples,
                interval_mask,
                interval_lag_mode=interval_lag_mode,
            )
            return curve, None, note

        subject_curves: list[np.ndarray] = []
        subject_diagonal_sds: list[np.ndarray] = []
        for subject_record in record.subject_records:
            subject_mask = occlusion_intervals(
                subject_record.time_axis_sec,
                occlusion_timing,
            )[interval_name]
            curve, diagonal_sd = compute_interval_lag_curve_with_diagonal_sd(
                subject_record.matrix,
                lag_samples,
                subject_mask,
                interval_lag_mode=interval_lag_mode,
            )
            subject_curves.append(curve)
            subject_diagonal_sds.append(diagonal_sd)

        curve_stack = np.vstack(subject_curves)
        line_curve = _nanmean_axis0(curve_stack)
        if shading_mode == "within-subject-diagonal-sd":
            return line_curve, _nanmean_axis0(np.vstack(subject_diagonal_sds)), None
        return line_curve, compute_across_subject_sem(curve_stack), None

    curve, diagonal_sd = compute_interval_lag_curve_with_diagonal_sd(
        record.matrix,
        lag_samples,
        interval_mask,
        interval_lag_mode=interval_lag_mode,
    )
    if shading_mode == "within-subject-diagonal-sd":
        return curve, diagonal_sd, None

    note = (
        f"{shading_mode} requested for subject-level plot {record.owner_label}; "
        "across-subject SEM requires a group plot, so no shaded band was drawn."
    )
    warnings.warn(note, RuntimeWarning)
    return curve, None, note


def _plot_interval_lag_overlay(
    records: Sequence[MatrixRecord],
    *,
    interval_name: str,
    interval_masks: Sequence[np.ndarray],
    occlusion_timing: Mapping[str, Any],
    lag_window_sec: float,
    shading_mode: str,
    interval_lag_mode: str,
    output_dir: Path,
    raw_root: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Write one interval lag-curve overlay figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    first = records[0]
    interval_dir = output_dir / "interval_lag_curves"
    if interval_lag_mode != "model-anchored":
        interval_dir = interval_dir / f"mode-{sanitize_label(interval_lag_mode)}"
    if shading_mode != "none":
        interval_dir = interval_dir / f"shading-{sanitize_label(shading_mode)}"
    out_path = (
        interval_dir / f"{_overlay_slug(records, interval_name)}_lag_curve.png"
    )
    write_target = _prepare_output_path(out_path, raw_root=raw_root, args=args)
    conditions = [record.condition for record in records]
    plot_notes: list[str] = []

    if write_target is not None:
        fig, (ax, legend_ax) = plt.subplots(
            ncols=2,
            figsize=(10.5, 4.8),
            constrained_layout=True,
            gridspec_kw={"width_ratios": [1.0, 0.38]},
        )
        legend_ax.axis("off")
        for record, interval_mask in zip(records, interval_masks):
            line_label = _legend_label_for_record(record)
            lags = _lag_samples(record.time_axis_sec, lag_window_sec)
            lag_sec = lags * sampling_step_seconds(record.time_axis_sec, DEFAULT_SAMPLING_RATE_HZ)
            curve, band, note = _interval_lag_summary_for_record(
                record,
                lag_samples=lags,
                interval_mask=interval_mask,
                interval_name=interval_name,
                occlusion_timing=occlusion_timing,
                shading_mode=shading_mode,
                interval_lag_mode=interval_lag_mode,
            )
            if note:
                plot_notes.append(note)
            (line,) = ax.plot(lag_sec, curve, label=line_label)
            if band is not None:
                finite_band = np.isfinite(curve) & np.isfinite(band)
                if np.any(finite_band):
                    ax.fill_between(
                        lag_sec,
                        curve - band,
                        curve + band,
                        where=finite_band,
                        color=line.get_color(),
                        alpha=0.18,
                        linewidth=0.0,
                        label="_nolegend_",
                    )
            if np.any(np.isfinite(curve)):
                best_idx = int(np.nanargmax(curve))
                best_lag = float(lag_sec[best_idx])
                best_value = float(curve[best_idx])
                ax.scatter(
                    [best_lag],
                    [best_value],
                    color=line.get_color(),
                    s=22,
                    zorder=3,
                )
                ax.annotate(
                    (
                        f"{_wrapped_legend_label(line_label)}\n"
                        f"{best_lag:.3f}s, {best_value:.3g}"
                    ),
                    xy=(best_lag, best_value),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                    color=line.get_color(),
                )
        ax.axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
        ax.axhline(0.0, color="0.45", linestyle=":", linewidth=0.9)
        ax.set_xlim(-float(lag_window_sec), float(lag_window_sec))
        ax.set_xlabel("Lag (neural time - model time) [s]")
        model_labels = sorted({record.stimulus_model_label for record in records})
        y_label = (
            f"Mean dRSA value - {model_labels[0]}"
            if len(model_labels) == 1
            else "Mean dRSA value - stimulus model comparison"
        )
        ax.set_ylabel(y_label)
        ax.set_title(
            _compact_plot_title(
                [
                    first.owner_label,
                    first.variant,
                    interval_name,
                    first.analysis_family,
                    first.sensor_group,
                    f"{first.variable} {first.slice_label}",
                ]
            ),
            fontsize=10,
        )
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            legend_ax.legend(
                handles,
                [_wrapped_legend_label(label) for label in labels],
                loc="center left",
                fontsize=8,
                title="Curves",
                title_fontsize=8,
                frameon=True,
                borderaxespad=0.0,
                handlelength=2.4,
            )
        fig.savefig(write_target, dpi=180, facecolor="white")
        plt.close(fig)

    return {
        "plot_type": "interval_lag_curve",
        "figure": str(out_path),
        "owner": first.owner_label,
        "conditions": conditions,
        "model_sources": sorted({record.model_source for record in records}),
        "stimulus_model_labels": sorted(
            {record.stimulus_model_label for record in records}
        ),
        "analysis_family": first.analysis_family,
        "sensor_group": first.sensor_group,
        "variable": first.variable,
        "slice_label": first.slice_label,
        "variant": first.variant,
        "interval": interval_name,
        "condition_count": len(set(conditions)),
        "interval_lag_mode": interval_lag_mode,
        "caption": (
            f"{first.owner_label} | {first.variant} | {interval_name} | "
            f"{first.analysis_family} | {first.sensor_group} | {first.variable} "
            f"{first.slice_label} | conditions/models: "
            + "; ".join(
                f"{record.condition} = {record.stimulus_model_label}"
                for record in records
            )
            + f" | interval mode: {_interval_lag_mode_interpretation(interval_lag_mode)}"
            + f" | shading: {_shading_interpretation(shading_mode)}"
            + (f" | notes: {'; '.join(plot_notes)}" if plot_notes else "")
        ),
        "shading_mode": shading_mode,
        "shading_interpretation": _shading_interpretation(shading_mode),
        "interval_lag_mode_interpretation": _interval_lag_mode_interpretation(
            interval_lag_mode
        ),
        "notes": plot_notes,
    }


def _include_plot_entry_in_report(entry: Mapping[str, Any]) -> bool:
    """
    Decide whether one generated figure should be listed in the HTML report.

    Data flow:
    all generated plot metadata -> report-facing subset. Real-neural interval
    lag reports are kept focused on plots that overlay more than one condition.
    Prediction-error models are also retained even though they are inherently
    single-condition occluded-deviant summaries. Simulations remain unchanged
    because their comparison structure does not always map onto the three
    real-neural conditions.
    """
    if entry.get("plot_type") != "interval_lag_curve":
        return True
    analysis_family = str(entry.get("analysis_family", ""))
    if analysis_family.startswith("simulations"):
        return True
    model_sources = {
        str(source)
        for source in entry.get("model_sources", ())
        if str(source).strip()
    }
    if "prediction_error" in model_sources:
        return True
    conditions = {
        str(condition)
        for condition in entry.get("conditions", ())
        if str(condition).strip()
    }
    return len(conditions) > 1


def _write_report(
    report_path: Path,
    *,
    title: str,
    summary_lines: Iterable[str],
    plot_entries: Sequence[Mapping[str, Any]],
    skipped_inputs: Sequence[Mapping[str, str]],
    raw_root: Path,
) -> Path:
    """
    Write the static Stage-08 report and embedded sibling copy.

    The report groups figures by plot type and exposes metadata filters in the
    captions so users can navigate by subject/group, condition, analysis family,
    sensor group, variable, and raw/smoothed-display variant.
    """
    report = guarded_path(report_path, raw_root)
    report.parent.mkdir(parents=True, exist_ok=True)
    generated_entries = [dict(item) for item in plot_entries]
    entries = [item for item in generated_entries if _include_plot_entry_in_report(item)]
    omitted_entries = len(generated_entries) - len(entries)

    toc = "\n".join(
        f"<li><a href=\"#{html_id(plot_type)}\">{html.escape(plot_type)}</a></li>"
        for plot_type in PLOT_TYPES
    )
    blocks: list[str] = []
    for plot_type in PLOT_TYPES:
        type_entries = [item for item in entries if item.get("plot_type") == plot_type]
        figures: list[str] = []
        for item in type_entries:
            href = image_href(report, item["figure"])
            caption_bits = [
                item.get("caption", item["figure"]),
                f"plot_type={plot_type}",
            ]
            caption = " | ".join(str(bit) for bit in caption_bits)
            figures.append(
                "<figure>"
                f"<img src=\"{html.escape(href)}\"/>"
                f"<figcaption>{caption}</figcaption>"
                "</figure>"
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

    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            "<style>",
            "body{font-family:system-ui,sans-serif;max-width:1600px;margin:2rem auto;padding:0 1rem;}",
            "nav,pre{background:#f5f5f5;padding:1rem;border-radius:6px;}",
            "figure{margin:1.5rem 0;padding-bottom:1rem;border-bottom:1px solid #ddd;}",
            "img{max-width:100%;height:auto;border:1px solid #ddd;}",
            "figcaption{font-size:.9rem;color:#333;margin-top:.4rem;}",
            "h2{margin-top:3rem;border-top:3px solid #222;padding-top:1rem;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<pre>{html.escape(chr(10).join(summary_lines))}</pre>",
            "<nav><strong>Plot types</strong><ul>",
            toc,
            "</ul></nav>",
            (
                "<p><strong>Report filter:</strong> non-simulation interval lag "
                "curves with only one condition are omitted unless they are "
                "prediction-error summaries; simulation interval lag curves are "
                "unchanged. "
                f"Omitted generated figures: {omitted_entries}.</p>"
            ),
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
    """Generate all Stage-08 plots and the HTML report for one owner scope."""
    raw_root = Path(cfg["resolved_paths"]["data_root"])
    output_dir = guarded_path(output_dir, raw_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    occlusion_timing = occlusion_timing_from_config(cfg)

    plot_entries: list[dict[str, Any]] = []
    overlay_groups: dict[tuple[str, ...], list[MatrixRecord]] = {}
    for record in records:
        overlay_groups.setdefault(_overlay_key(record), []).append(record)

    # Family 1: interval lag curves. One figure per interval and shared context,
    # overlaying conditions when the Stage-07 output set supports that grouping.
    for group in overlay_groups.values():
        group = sorted(
            group,
            key=lambda item: (
                CONDITIONS.index(item.condition) if item.condition in CONDITIONS else 99,
                item.stimulus_model_label,
            ),
        )
        for interval_name in ("pre_occlusion", "during_occlusion", "post_occlusion"):
            masks = [
                occlusion_intervals(record.time_axis_sec, occlusion_timing)[interval_name]
                for record in group
            ]
            plot_entries.append(
                _plot_interval_lag_overlay(
                    group,
                    interval_name=interval_name,
                    interval_masks=masks,
                    occlusion_timing=occlusion_timing,
                    lag_window_sec=args.lag_window_sec,
                    shading_mode=args.shading,
                    interval_lag_mode=args.interval_lag_mode,
                    output_dir=output_dir,
                    raw_root=raw_root,
                    args=args,
                )
            )

    first_sec = occlusion_timing["first_complete_occlusion_sec"]
    reapp_sec = occlusion_timing["reappearance_sec"]
    summary_lines = [
        "Stage status: completed",
        f"Generated UTC: {utc_now()}",
        "Lag convention: lag = neural time - model time",
        "Matrix interpretation: x = neural time, y = model time, indexed matrix[model_time, neural_time]",
        (
            "Occlusion boundaries: "
            f"first_complete_occlusion frame {occlusion_timing['first_complete_occlusion_frame']} "
            f"= {first_sec:.6f} s; "
            f"reappearance frame {occlusion_timing['reappearance_frame']} "
            f"= {reapp_sec:.6f} s"
        ),
        (
            "Intervals on model time: "
            f"pre_occlusion [0, {first_sec:.6f}); "
            f"during_occlusion [{first_sec:.6f}, {reapp_sec:.6f}); "
            f"post_occlusion [{reapp_sec:.6f}, end]"
        ),
        f"Lag window: +/-{args.lag_window_sec:g} s",
        f"Interval lag-curve mode: {args.interval_lag_mode}",
        (
            "Interval lag-curve mode interpretation: "
            f"{_interval_lag_mode_interpretation(args.interval_lag_mode)}"
        ),
        f"Interval lag-curve shading mode: {args.shading}",
        f"Interval lag-curve shading interpretation: {_shading_interpretation(args.shading)}",
        (
            "Report figure selection: non-simulation interval lag curves include "
            "multi-condition overlays only; simulation interval lag curves remain "
            "unchanged."
        ),
        f"Matrix records: {len(records)}",
        f"Figures: {len(plot_entries)}",
        f"Skipped inputs: {len(skipped_inputs)}",
        f"Output directory: {output_dir}",
    ]
    variants = sorted({record.variant for record in records})
    summary_lines.append(
        "Variants: "
        + ", ".join(variants)
        + " (raw = unsmoothed Stage-07 matrices; smoothed_display = "
        "display-only Python smoothing when requested)"
    )
    model_labels = sorted({record.stimulus_model_label for record in records})
    summary_lines.append(f"Stimulus/comparison model labels: {', '.join(model_labels)}")

    _write_report(
        report_html,
        title=title,
        summary_lines=summary_lines,
        plot_entries=plot_entries,
        skipped_inputs=skipped_inputs,
        raw_root=raw_root,
    )
    return {
        "status": "completed",
        "record_count": len(records),
        "figure_count": len(plot_entries),
        "output_dir": str(output_dir),
        "report_html": str(report_html),
        "report_embedded_html": str(report_html.with_name(f"{report_html.stem}_embedded.html")),
        "plot_entries": plot_entries,
        "skipped_inputs": list(skipped_inputs),
    }


def run_subject_scope(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Build Stage-08 outputs for one resolved subject configuration."""
    owner_label = str(cfg["layout"]["subject_bids_id"])
    paths = _analysis_paths(cfg)
    smoothing_window = int(
        args.plot_smoothing_window
        or cfg.get("matlab_drsa", {}).get("plot_smoothing_window", 9)
    )
    records, skipped = _read_subject_records(
        cfg,
        owner_label=owner_label,
        include_smoothed=args.include_smoothed,
        smoothing_window_samples=smoothing_window,
        conditions=_selected_conditions(args),
        sensor_groups=_selected_sensor_groups(args),
    )
    if not records:
        raise RuntimeError(f"No Stage-07 result matrices were found for {owner_label}.")
    return _plot_records(
        records,
        output_dir=paths["subject_line_dir"],
        report_html=paths["subject_report_html"],
        title=f"{owner_label} | 08 dRSA line summaries",
        cfg=cfg,
        args=args,
        skipped_inputs=skipped,
    )


def run_group_scope(
    base_cfg: Mapping[str, Any],
    args: argparse.Namespace,
    subjects: Sequence[str],
) -> dict[str, Any]:
    """Build group-average Stage-08 outputs from existing subject MAT files."""
    smoothing_window = int(
        args.plot_smoothing_window
        or base_cfg.get("matlab_drsa", {}).get("plot_smoothing_window", 9)
    )
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
            include_smoothed=args.include_smoothed,
            smoothing_window_samples=smoothing_window,
            conditions=_selected_conditions(args),
            sensor_groups=_selected_sensor_groups(args),
        )
        subject_records.extend(records)
        skipped.extend(subject_skipped)

    records, group_skipped = _group_records(subject_records)
    skipped.extend(group_skipped)
    if not records:
        raise RuntimeError(
            "No group-average line-summary records were produced. "
            f"Skipped inputs: {skipped[:5]}"
        )
    paths = _analysis_paths(base_cfg)
    return _plot_records(
        records,
        output_dir=paths["group_line_dir"],
        report_html=paths["group_report_html"],
        title=f"Group | {base_cfg['layout']['session_bids_id']} | 08 dRSA line summaries",
        cfg=base_cfg,
        args=args,
        skipped_inputs=skipped,
    )


def main() -> int:
    """Resolve configuration, generate requested Stage-08 summaries, and exit."""
    args = build_arg_parser().parse_args()
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if args.group_only and args.subjects_only:
        raise SystemExit("--group-only and --subjects-only are mutually exclusive.")
    if args.lag_window_sec <= 0:
        raise SystemExit("--lag-window-sec must be > 0.")

    base_subject_config = args.subject_config
    if args.subjects:
        base_subject_config = _subject_config_for_id(args, normalize_subject_label(args.subjects[0]))
    base_cfg = _load_cfg(args, subject_config=base_subject_config)
    subjects = _subjects_from_args_or_config(args, base_cfg)

    reports: dict[str, Any] = {
        "stage": "drsa_line_summaries",
        "generated_utc": utc_now(),
        "subjects": subjects,
        "include_smoothed": bool(args.include_smoothed),
        "shading": args.shading,
        "interval_lag_mode": args.interval_lag_mode,
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
        print(f"dRSA line-summary error: {exc}")
        return 2

    if args.print_summary:
        printable = {
            "stage": reports["stage"],
            "subjects": reports["subjects"],
            "include_smoothed": reports["include_smoothed"],
            "shading": reports["shading"],
            "interval_lag_mode": reports["interval_lag_mode"],
            "reports": {
                scope: (
                    {
                        key: {
                            "record_count": value["record_count"],
                            "figure_count": value["figure_count"],
                            "output_dir": value["output_dir"],
                            "report_html": value["report_html"],
                            "report_embedded_html": value["report_embedded_html"],
                        }
                        for key, value in report.items()
                    }
                    if scope == "subjects"
                    else {
                        "record_count": report["record_count"],
                        "figure_count": report["figure_count"],
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
