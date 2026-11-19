#!/usr/bin/env python3
"""
10_drsa_stats.py

Purpose
Generate optional group-level statistical tests from existing Stage-07 dRSA
result MAT files and Stage-08-compatible line summaries. This script does not
run MATLAB, does not recompute dRSA, and does not modify raw data. It reads
saved subject-level dRSA matrices, assembles matched subject stacks, runs
cluster-permutation tests, and writes machine-readable statistics plus a static
HTML report:

- 2-D matrix tests use each dRSA matrix as matrix[model_time, neural_time].
  Stage-07 plots x = neural time and y = model time, so report figures keep the
  same interpretation.
- 1-D line-summary tests compute interval lag curves with the Stage-08 helper.
  Paired line-summary contrasts first build each subject's
  deviant-minus-nondeviant matrix, then aggregate that difference matrix into a
  lag curve so NaN support and interval masking match the 2-D paired contrast.
  The lag convention is:

      lag = neural time - model time
      sample = matrix[model_time_index, model_time_index + lag]

Usage examples
1) Run the default group stats for subjects 04-39 from the configured
   derivative tree:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/10_drsa_stats.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subjects {04..39} \
     --overwrite --print-summary

2) Run only occipital A_dRSA_alt line-summary stats with a +/-0.5 s lag range:

   python .../10_drsa_stats.py --config .../pipeline_default.yaml \
     --subjects {04..39} \
     --stats-targets line-summaries --analysis-families A_dRSA_alt \
     --sensor-groups occipital --line-lag-window-sec 0.5 --overwrite

3) Include display-smoothed matrices as explicitly labeled additional inputs:

   python .../10_drsa_stats.py --config .../pipeline_default.yaml \
     --subjects {04..39} \
     --include-smoothed --plot-smoothing-window 9 --overwrite

4) Run on the mounted MORWUR derivative tree with a site override:

   python .../10_drsa_stats.py --config .../pipeline_default.yaml \
     --site-config .../pipeline_cluster.yaml --subjects {04..39} \
     --resume --print-summary

Inputs
- Existing Stage-07 result MAT files below
  derivatives/analysis/sub-XX/ses-01/drsa/.
- pipeline_config.py-resolved data and derivative roots.
- Stage-08 lag-summary conventions from 08_plot_drsa_line_summaries.py.

Outputs
- Statistical arrays, figures, JSON, and CSV tables:
  derivatives/analysis/group/ses-01/drsa/stats/
- Linked and embedded reports:
  derivatives/reports/group/ses-01/10_drsa_stats_report.html
  derivatives/reports/group/ses-01/10_drsa_stats_report_embedded.html

Safety and assumptions
- Every output path is passed through the shared raw-tree guard.
- Subject-level Stage-07 MAT files are the only scientific inputs. PNG files are
  never scraped.
- ``--include-smoothed`` creates display-smoothed inputs in memory only and
  labels those tests as ``smoothed_display``. Smoothed-display tests are not
  independent of raw tests.
- The default model filter keeps position-like matrices plus prediction-error
  and explicit PE-variant matrices from A_dRSA and A_dRSA_alt. Simulations are
  excluded unless requested with ``--analysis-families``. Use ``--models``
  and/or ``--model-sources`` to override model selection.
- Cluster tests use one-sample sign-flip permutations. The default paired
  contrast tests actual observed position for occluded deviant minus occluded
  nondeviant. Paired contrasts are first differenced within subject, then
  tested against zero. Subjects lacking either side of a paired contrast are
  dropped and reported.
- MNE cluster functions are preferred when available and the data have complete
  finite support. A documented NumPy/SciPy sign-flip implementation is used for
  NaN-containing stacks or when MNE is unavailable.
Current default cohort is P04-P14/P16-P39 (2026-11-07 P15 coverage decision);
explicit historical inputs remain readable for provenance.
"""

from __future__ import annotations

import argparse
import csv
import html
import importlib.util
import itertools
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
    html_id,
    image_href,
    iter_matrix_slices,
    model_feature_label,
    model_source_for_matrix_slice,
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
from pipeline_config import WORKING_SUBJECT_IDS, load_and_resolve_config
from preprocessing_stage_utils import guarded_path


DEFAULT_SAMPLING_RATE_HZ = 120.0
DEFAULT_SUBJECTS = tuple(f"sub-{subject_id}" for subject_id in WORKING_SUBJECT_IDS)
STATS_TARGETS = ("matrices", "line-summaries")
TAIL_CHOICES = ("two-sided", "positive", "negative")
DEFAULT_LINE_INTERVALS = ("pre_occlusion", "during_occlusion", "post_occlusion")
VARIANT_RAW = "raw"
VARIANT_SMOOTHED = "smoothed_display"
DEFAULT_MODEL_SELECTORS = ("position", "prediction_error")
DEFAULT_ANALYSIS_FAMILIES = ("A_dRSA", "A_dRSA_alt")
PAIRED_CONTRASTS = ("none", "actual_deviant_position_vs_nondeviant")
REPORT_CLUSTER_TABLE_MODES = ("significant", "all", "none")


@dataclass(frozen=True)
class MatrixRecord:
    """
    One subject-level 2-D dRSA matrix slice plus selection metadata.

    Data flow:
    Stage-07 result MAT variable -> matrix slice -> group test stack.
    """

    subject: str
    condition: str
    analysis_family: str
    sensor_group: str
    drsa_type: str
    model_source: str
    model_feature: str
    stimulus_model_label: str
    variable: str
    slice_index: int
    slice_label: str
    variant: str
    matrix: np.ndarray
    time_axis_sec: np.ndarray
    relative_result: str
    source_path: str


@dataclass(frozen=True)
class ClusterInfo:
    """One cluster from a cluster-permutation result."""

    cluster_id: int
    mask: np.ndarray
    p_value: float
    statistic: float
    size: int
    sign: str


@dataclass(frozen=True)
class ClusterTestResult:
    """
    Complete result for one one-sample sign-flip cluster test.

    ``observations`` is indexed as subject x feature(s). A paired contrast is
    represented by within-subject differences before this object is created.
    """

    statistic_map: np.ndarray
    clusters: tuple[ClusterInfo, ...]
    threshold: float
    cluster_alpha: float
    significance_alpha: float
    tail: str
    n_permutations: int
    backend: str
    finite_observation_counts: np.ndarray
    null_max_statistics: np.ndarray


@dataclass(frozen=True)
class PreparedTest:
    """One matrix or line-summary statistical test ready for output/reporting."""

    slug: str
    stats_target: str
    test_type: str
    variant: str
    analysis_family: str
    sensor_group: str
    condition: str
    contrast: str
    model_source: str
    model_feature: str
    stimulus_model_label: str
    variable: str
    slice_label: str
    interval: str | None
    subjects: tuple[str, ...]
    dropped_subjects: tuple[str, ...]
    observations: np.ndarray
    mean_map: np.ndarray
    time_axis_sec: np.ndarray | None = None
    lag_sec: np.ndarray | None = None
    source_paths: tuple[str, ...] = field(default_factory=tuple)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the Stage-10 dRSA stats command-line interface."""
    stage08 = _load_stage08_module()
    parser = argparse.ArgumentParser(
        description="Run optional group-level cluster stats on existing dRSA outputs."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--subjects", nargs="+", default=None)
    parser.add_argument(
        "--conditions",
        nargs="+",
        choices=CONDITIONS,
        default=None,
        help="Condition filter. Defaults to all available conditions.",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help=(
            "Optional model/feature selector, e.g. position direction "
            "prediction_error pe_prior_constrained_50. Defaults to position "
            "plus prediction_error."
        ),
    )
    parser.add_argument(
        "--model-sources",
        nargs="+",
        default=None,
        help=(
            "Optional model-source filter, e.g. observed predicted "
            "prediction_error pe_prior_constrained_50."
        ),
    )
    parser.add_argument(
        "--sensor-groups",
        nargs="+",
        default=None,
        help="Optional sensor-group filter. Defaults to all available Stage-07 outputs.",
    )
    parser.add_argument(
        "--analysis-families",
        nargs="+",
        default=list(DEFAULT_ANALYSIS_FAMILIES),
        help="Optional analysis-family filter, e.g. A_dRSA A_dRSA_alt simulations_model_vs_model.",
    )
    parser.add_argument(
        "--stats-targets",
        nargs="+",
        choices=STATS_TARGETS,
        default=list(STATS_TARGETS),
        help="Statistical targets to generate.",
    )
    parser.add_argument(
        "--line-intervals",
        nargs="+",
        choices=DEFAULT_LINE_INTERVALS,
        default=list(DEFAULT_LINE_INTERVALS),
        help="Stage-08 interval summaries to test.",
    )
    parser.add_argument(
        "--line-lag-window-sec",
        type=float,
        default=1.0,
        help="Maximum absolute lag, in seconds, for line-summary tests.",
    )
    parser.add_argument(
        "--line-interval-lag-mode",
        choices=stage08.INTERVAL_LAG_MODES,
        default="model-anchored",
        help=(
            "Interval lag support. The default preserves Stage-08 behavior: "
            "model-time rows are restricted; neural-time columns may fall "
            "anywhere inside the full matrix."
        ),
    )
    parser.add_argument(
        "--include-smoothed",
        action="store_true",
        help="Also test display-smoothed matrices, clearly labeled as smoothed_display.",
    )
    parser.add_argument(
        "--paired-contrasts",
        nargs="+",
        choices=PAIRED_CONTRASTS,
        default=["actual_deviant_position_vs_nondeviant"],
        help=(
            "Optional paired contrasts. Defaults to the actual observed "
            "position contrast occluded_deviant - occluded_nondeviant. "
            "Use none to run only one-sample tests."
        ),
    )
    parser.add_argument(
        "--report-cluster-tables",
        choices=REPORT_CLUSTER_TABLE_MODES,
        default="significant",
        help=(
            "HTML report cluster table display. Machine-readable CSV/JSON "
            "always retain all clusters."
        ),
    )
    parser.add_argument(
        "--plot-smoothing-window",
        type=int,
        default=None,
        help="Gaussian smoothing window, in samples, for --include-smoothed.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=1024,
        help="Number of sign-flip permutations; exact enumeration is used when smaller.",
    )
    parser.add_argument(
        "--cluster-alpha",
        type=float,
        default=0.05,
        help="Cluster-forming p-value used to derive the t threshold.",
    )
    parser.add_argument(
        "--significance-alpha",
        type=float,
        default=0.05,
        help="Cluster p-value threshold used for report overlays.",
    )
    parser.add_argument("--tail", choices=TAIL_CHOICES, default="two-sided")
    parser.add_argument("--seed", type=int, default=97)
    parser.add_argument(
        "--stats-backend",
        choices=("auto", "mne", "numpy"),
        default="auto",
        help="Cluster-test backend. auto prefers MNE when finite support allows it.",
    )
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


_STAGE08_MODULE: Any | None = None


def _load_stage08_module() -> Any:
    """
    Load the numbered Stage-08 module so line stats use the same lag helper.

    Python cannot import ``08_plot_drsa_line_summaries`` with a normal import
    statement because the filename starts with a digit. The dynamic import keeps
    Stage-10 behavior tied to the tested Stage-08 implementation without moving
    shared code during this focused stats change.
    """
    global _STAGE08_MODULE  # pylint: disable=global-statement
    if _STAGE08_MODULE is not None:
        return _STAGE08_MODULE
    module_path = Path(__file__).resolve().with_name("08_plot_drsa_line_summaries.py")
    spec = importlib.util.spec_from_file_location("stage08_drsa_line_summaries", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load Stage-08 module from {module_path}.")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(spec.name, module)
    spec.loader.exec_module(module)
    _STAGE08_MODULE = module
    return module


def _analysis_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return subject input and group Stage-10 output paths from config."""
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
        "group_analysis_dir": group_analysis_dir,
        "group_reports_dir": group_reports_dir,
        "group_stats_dir": group_analysis_dir / "drsa" / "stats",
        "group_report_html": group_reports_dir / "10_drsa_stats_report.html",
    }


def _load_cfg(
    args: argparse.Namespace,
    *,
    subject_config: str | Path | None = None,
) -> Mapping[str, Any]:
    """Load one resolved pipeline configuration for Stage-10 use."""
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


def _subjects_from_args_or_config(args: argparse.Namespace, cfg: Mapping[str, Any]) -> list[str]:
    """Resolve group subjects; defaults to the configured group or 04..39."""
    if args.subjects:
        return [normalize_subject_label(value) for value in args.subjects]
    configured = cfg.get("matlab_drsa", {}).get("group_average_subjects")
    if configured:
        return [normalize_subject_label(value) for value in configured]
    return list(DEFAULT_SUBJECTS)


def _sampling_rate_for_plots(cfg: Mapping[str, Any]) -> float:
    """Resolve the fallback sample rate used when a MAT lacks time_axis_sec."""
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
    """Filter out non-result MAT files and generated post-processing folders."""
    text = str(path)
    return (
        "desc-identityAveraged_neural" in path.name
        or "/line_summaries/" in text
        or "line_summaries" in path.parts
        or "/sliding_lag_surfaces/" in text
        or "sliding_lag_surfaces" in path.parts
        or "/stats/" in text
        or "stats" in path.parts
        or "/group_average/" in text
        or "group_average" in path.parts
    )


def _selected_conditions(args: argparse.Namespace) -> set[str] | None:
    """Return the optional condition filter."""
    return set(args.conditions) if args.conditions else None


def _selected_sensor_groups(args: argparse.Namespace) -> set[str] | None:
    """Return the optional normalized sensor-group filter."""
    if not args.sensor_groups:
        return None
    return {normalize_sensor_group(value) for value in args.sensor_groups}


def _selected_analysis_families(args: argparse.Namespace) -> set[str] | None:
    """Return the optional analysis-family filter."""
    return set(args.analysis_families) if args.analysis_families else None


def _paired_deviant_nondeviant_enabled(args: argparse.Namespace) -> bool:
    """Return whether the actual-position deviant-minus-nondeviant contrast is enabled."""
    contrasts = set(args.paired_contrasts or ())
    return "actual_deviant_position_vs_nondeviant" in contrasts


def _record_matches_actual_position_contrast(record: MatrixRecord) -> bool:
    """Return whether one record belongs in the actual-position paired contrast."""
    return (
        record.condition in {"occluded_deviant", "occluded_nondeviant"}
        and record.model_source == "observed"
        and record.model_feature == "position"
    )


def _selector_text(values: Iterable[Any]) -> set[str]:
    """Normalize model selector text for permissive CLI matching."""
    return {str(value).strip().lower().replace("_", " ") for value in values}


def _record_matches_models(record: MatrixRecord, args: argparse.Namespace) -> bool:
    """
    Apply model/feature filters to one matrix record.

    The default follows the requested initial scope: position-like matrices and
    prediction-error matrices. Explicit ``--models`` selectors match model
    feature, model source, variable, slice label, or stimulus model label.
    """
    if args.model_sources:
        allowed_sources = _selector_text(args.model_sources)
        if record.model_source.lower().replace("_", " ") not in allowed_sources:
            return False

    selectors = _selector_text(args.models or DEFAULT_MODEL_SELECTORS)
    candidates = {
        record.model_feature.lower(),
        record.model_source.lower().replace("_", " "),
        record.model_source.lower(),
        record.variable.lower(),
        str(record.slice_label).lower(),
        record.stimulus_model_label.lower(),
    }
    candidates.update(record.stimulus_model_label.lower().replace("/", " ").split())
    candidates.update(record.stimulus_model_label.lower().replace("-", " ").split())
    if selectors.intersection(candidates):
        return True
    return any(selector in record.stimulus_model_label.lower() for selector in selectors)


def _read_subject_records(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
    *,
    subject: str,
) -> tuple[list[MatrixRecord], list[dict[str, str]]]:
    """
    Discover and read selected Stage-07 result matrices for one subject.

    Data flow:
    subject dRSA folder -> result MAT files -> matrix variables/slices -> raw
    and optional smoothed-display MatrixRecord objects.
    """
    paths = _analysis_paths(cfg)
    drsa_dir = paths["drsa_dir"]
    sampling_rate = _sampling_rate_for_plots(cfg)
    smoothing_window = int(
        args.plot_smoothing_window
        or cfg.get("matlab_drsa", {}).get("plot_smoothing_window", 9)
    )
    conditions = _selected_conditions(args)
    sensor_groups = _selected_sensor_groups(args)
    analysis_families = _selected_analysis_families(args)
    records: list[MatrixRecord] = []
    skipped: list[dict[str, str]] = []

    if not drsa_dir.is_dir():
        skipped.append({"subject": subject, "reason": f"missing dRSA folder: {drsa_dir}"})
        return records, skipped

    for mat_path in sorted(drsa_dir.rglob("*.mat")):
        if _should_skip_result_path(mat_path):
            continue
        relative = mat_path.relative_to(drsa_dir)
        condition = condition_from_path(relative)
        sensor_group = sensor_group_from_path(relative)
        analysis_family = analysis_family_from_path(relative)
        drsa_type = drsa_type_from_path(relative)
        path_model_source = model_source_from_path(relative)
        if conditions is not None and condition not in conditions:
            continue
        if sensor_groups is not None and sensor_group not in sensor_groups:
            continue
        if analysis_families is not None and analysis_family not in analysis_families:
            continue
        try:
            matrices, time_axis, metadata = read_drsa_result_file(mat_path)
        except RuntimeError as exc:
            skipped.append({"subject": subject, "path": str(mat_path), "reason": str(exc)})
            continue

        for variable, matrix_array in matrices.items():
            model_names = (
                metadata.get("model_names")
                if variable == "dRSA_neural_correlation"
                else None
            )
            for slice_index, slice_label, matrix_2d in iter_matrix_slices(
                matrix_array,
                slice_labels=model_names,
            ):
                model_source = model_source_for_matrix_slice(
                    condition,
                    path_model_source,
                    variable,
                    slice_label=slice_label,
                    analysis_family=analysis_family,
                )
                feature = model_feature_label(
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
                axis = time_axis_for_matrix(
                    matrix_2d,
                    time_axis=time_axis,
                    sampling_rate_hz=sampling_rate,
                )
                common = {
                    "subject": subject,
                    "condition": condition,
                    "analysis_family": analysis_family,
                    "sensor_group": sensor_group,
                    "drsa_type": drsa_type,
                    "model_source": model_source,
                    "model_feature": feature,
                    "stimulus_model_label": model_label,
                    "variable": variable,
                    "slice_index": int(slice_index),
                    "slice_label": str(slice_label),
                    "time_axis_sec": np.asarray(axis, dtype=float),
                    "relative_result": str(relative),
                    "source_path": str(mat_path),
                }
                raw_record = MatrixRecord(
                    **common,
                    variant=VARIANT_RAW,
                    matrix=np.asarray(matrix_2d, dtype=float),
                )
                if _record_matches_models(raw_record, args):
                    records.append(raw_record)
                if args.include_smoothed:
                    smooth_record = MatrixRecord(
                        **common,
                        variant=VARIANT_SMOOTHED,
                        matrix=smooth_matrix_for_plot(matrix_2d, smoothing_window),
                    )
                    if _record_matches_models(smooth_record, args):
                        records.append(smooth_record)

    return records, skipped


def _matrix_group_key(record: MatrixRecord) -> tuple[str, ...]:
    """Return the one-sample grouping key for matching matrix records."""
    return (
        record.variant,
        record.analysis_family,
        record.sensor_group,
        record.drsa_type,
        record.condition,
        record.model_source,
        record.model_feature,
        record.variable,
        record.slice_label,
    )


def _paired_group_key(record: MatrixRecord) -> tuple[str, ...]:
    """Return the paired-contrast grouping key, intentionally omitting condition."""
    return (
        record.variant,
        record.analysis_family,
        record.sensor_group,
        record.drsa_type,
        record.model_source,
        record.model_feature,
        record.variable,
    )


def _lag_samples(time_axis_sec: np.ndarray, lag_window_sec: float) -> np.ndarray:
    """Build integer lag offsets for the requested line-summary lag window."""
    dt = sampling_step_seconds(time_axis_sec, DEFAULT_SAMPLING_RATE_HZ)
    max_lag = int(round(float(lag_window_sec) / dt))
    return np.arange(-max_lag, max_lag + 1, dtype=int)


def _line_curve_for_matrix(
    matrix_2d: np.ndarray,
    time_axis_sec: np.ndarray,
    *,
    interval: str,
    occlusion_timing: Mapping[str, Any],
    lag_window_sec: float,
    interval_lag_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute one Stage-08-compatible interval lag curve for a matrix.

    Data flow:
    2-D matrix -> interval model-time mask -> integer lag samples -> Stage-08
    diagonal averaging helper. Paired line summaries pass an already
    differenced subject matrix here so the NaN support is established before
    interval/lag averaging.
    """
    stage08 = _load_stage08_module()
    axis = np.asarray(time_axis_sec, dtype=float)
    interval_mask = occlusion_intervals(axis, occlusion_timing)[interval]
    lags = _lag_samples(axis, lag_window_sec)
    lag_sec = lags * sampling_step_seconds(axis, DEFAULT_SAMPLING_RATE_HZ)
    curve = stage08.compute_interval_lag_curve(
        matrix_2d,
        lags,
        interval_mask,
        interval_lag_mode=interval_lag_mode,
    )
    return lag_sec.astype(float), np.asarray(curve, dtype=float)


def _line_curve_for_record(
    record: MatrixRecord,
    *,
    interval: str,
    occlusion_timing: Mapping[str, Any],
    lag_window_sec: float,
    interval_lag_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute one Stage-08-compatible interval lag curve for a matrix record."""
    return _line_curve_for_matrix(
        record.matrix,
        record.time_axis_sec,
        interval=interval,
        occlusion_timing=occlusion_timing,
        lag_window_sec=lag_window_sec,
        interval_lag_mode=interval_lag_mode,
    )


def _slug_for_parts(parts: Iterable[Any]) -> str:
    """Build a stable filename slug from metadata parts."""
    return "_".join(sanitize_label(part) for part in parts if str(part).strip())


def _records_by_subject(records: Sequence[MatrixRecord]) -> dict[str, MatrixRecord]:
    """Keep one deterministic record per subject for a matched test group."""
    by_subject: dict[str, MatrixRecord] = {}
    for record in sorted(records, key=lambda item: (item.subject, item.relative_result)):
        by_subject.setdefault(record.subject, record)
    return by_subject


def _prepare_matrix_one_sample_tests(
    records: Sequence[MatrixRecord],
    subjects: Sequence[str],
    *,
    skipped: list[dict[str, str]],
) -> list[PreparedTest]:
    """Assemble one-sample 2-D matrix tests from subject-level records."""
    grouped: dict[tuple[str, ...], list[MatrixRecord]] = {}
    for record in records:
        grouped.setdefault(_matrix_group_key(record), []).append(record)

    tests: list[PreparedTest] = []
    for key, items in sorted(grouped.items()):
        by_subject = _records_by_subject(items)
        matched_subjects = tuple(subject for subject in subjects if subject in by_subject)
        if len(matched_subjects) < 2:
            skipped.append({"target": "matrices", "key": repr(key), "reason": "fewer than two matched subjects"})
            continue
        shapes = {tuple(by_subject[subject].matrix.shape) for subject in matched_subjects}
        if len(shapes) != 1:
            skipped.append({"target": "matrices", "key": repr(key), "reason": f"shape mismatch: {sorted(shapes)}"})
            continue
        axes = [by_subject[subject].time_axis_sec for subject in matched_subjects]
        if any(axis.shape != axes[0].shape or not np.allclose(axis, axes[0]) for axis in axes[1:]):
            skipped.append({"target": "matrices", "key": repr(key), "reason": "time-axis mismatch"})
            continue

        stack = np.stack([by_subject[subject].matrix for subject in matched_subjects], axis=0)
        first = by_subject[matched_subjects[0]]
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_map = np.nanmean(stack, axis=0)
        tests.append(
            PreparedTest(
                slug=_slug_for_parts(("matrices", "onesample", *key)),
                stats_target="matrices",
                test_type="one_sample",
                variant=first.variant,
                analysis_family=first.analysis_family,
                sensor_group=first.sensor_group,
                condition=first.condition,
                contrast=f"{first.condition} vs 0",
                model_source=first.model_source,
                model_feature=first.model_feature,
                stimulus_model_label=first.stimulus_model_label,
                variable=first.variable,
                slice_label=first.slice_label,
                interval=None,
                subjects=matched_subjects,
                dropped_subjects=tuple(subject for subject in subjects if subject not in by_subject),
                observations=stack,
                mean_map=mean_map,
                time_axis_sec=first.time_axis_sec,
                source_paths=tuple(by_subject[subject].source_path for subject in matched_subjects),
            )
        )
    return tests


def _prepare_matrix_paired_tests(
    records: Sequence[MatrixRecord],
    subjects: Sequence[str],
    *,
    skipped: list[dict[str, str]],
) -> list[PreparedTest]:
    """Assemble actual-position paired deviant-minus-nondeviant 2-D matrix tests."""
    by_key_condition: dict[tuple[str, ...], dict[str, list[MatrixRecord]]] = {}
    for record in records:
        if not _record_matches_actual_position_contrast(record):
            continue
        by_key_condition.setdefault(_paired_group_key(record), {}).setdefault(
            record.condition,
            [],
        ).append(record)

    tests: list[PreparedTest] = []
    for key, condition_items in sorted(by_key_condition.items()):
        deviant = _records_by_subject(condition_items.get("occluded_deviant", []))
        nondeviant = _records_by_subject(condition_items.get("occluded_nondeviant", []))
        matched_subjects = tuple(
            subject for subject in subjects if subject in deviant and subject in nondeviant
        )
        dropped = tuple(subject for subject in subjects if subject not in matched_subjects)
        if len(matched_subjects) < 2:
            skipped.append({"target": "matrices", "key": repr(key), "reason": "paired contrast has fewer than two matched subjects"})
            continue
        shapes = {
            tuple(deviant[subject].matrix.shape)
            for subject in matched_subjects
        } | {
            tuple(nondeviant[subject].matrix.shape)
            for subject in matched_subjects
        }
        if len(shapes) != 1:
            skipped.append({"target": "matrices", "key": repr(key), "reason": f"paired shape mismatch: {sorted(shapes)}"})
            continue
        first = deviant[matched_subjects[0]]
        if any(
            not np.allclose(deviant[subject].time_axis_sec, first.time_axis_sec)
            or not np.allclose(nondeviant[subject].time_axis_sec, first.time_axis_sec)
            for subject in matched_subjects
        ):
            skipped.append({"target": "matrices", "key": repr(key), "reason": "paired time-axis mismatch"})
            continue

        stack = np.stack(
            [
                deviant[subject].matrix - nondeviant[subject].matrix
                for subject in matched_subjects
            ],
            axis=0,
        )
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_map = np.nanmean(stack, axis=0)
        tests.append(
            PreparedTest(
                slug=_slug_for_parts(("matrices", "actual-position-dev-minus-nondev", *key)),
                stats_target="matrices",
                test_type="paired",
                variant=first.variant,
                analysis_family=first.analysis_family,
                sensor_group=first.sensor_group,
                condition="occluded_deviant_minus_occluded_nondeviant",
                contrast="actual deviant position - actual nondeviant position",
                model_source=first.model_source,
                model_feature=first.model_feature,
                stimulus_model_label=first.stimulus_model_label,
                variable=first.variable,
                slice_label=first.slice_label,
                interval=None,
                subjects=matched_subjects,
                dropped_subjects=dropped,
                observations=stack,
                mean_map=mean_map,
                time_axis_sec=first.time_axis_sec,
                source_paths=tuple(
                    path
                    for subject in matched_subjects
                    for path in (deviant[subject].source_path, nondeviant[subject].source_path)
                ),
            )
        )
    return tests


def _prepare_line_one_sample_tests(
    records: Sequence[MatrixRecord],
    subjects: Sequence[str],
    args: argparse.Namespace,
    cfg: Mapping[str, Any],
    *,
    skipped: list[dict[str, str]],
) -> list[PreparedTest]:
    """Assemble one-sample 1-D interval lag-curve tests."""
    grouped: dict[tuple[str, ...], list[MatrixRecord]] = {}
    for record in records:
        grouped.setdefault(_matrix_group_key(record), []).append(record)
    timing = occlusion_timing_from_config(cfg)

    tests: list[PreparedTest] = []
    for key, items in sorted(grouped.items()):
        by_subject = _records_by_subject(items)
        matched_subjects = tuple(subject for subject in subjects if subject in by_subject)
        if len(matched_subjects) < 2:
            skipped.append({"target": "line-summaries", "key": repr(key), "reason": "fewer than two matched subjects"})
            continue
        for interval in args.line_intervals:
            curves: list[np.ndarray] = []
            lag_axis: np.ndarray | None = None
            axis_mismatch = False
            for subject in matched_subjects:
                lag_sec, curve = _line_curve_for_record(
                    by_subject[subject],
                    interval=interval,
                    occlusion_timing=timing,
                    lag_window_sec=args.line_lag_window_sec,
                    interval_lag_mode=args.line_interval_lag_mode,
                )
                if lag_axis is None:
                    lag_axis = lag_sec
                elif lag_sec.shape != lag_axis.shape or not np.allclose(lag_sec, lag_axis):
                    axis_mismatch = True
                    break
                curves.append(curve)
            if axis_mismatch or lag_axis is None:
                skipped.append({"target": "line-summaries", "key": repr(key), "interval": interval, "reason": "lag-axis mismatch"})
                continue
            stack = np.vstack(curves)
            first = by_subject[matched_subjects[0]]
            with np.errstate(invalid="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                mean_curve = np.nanmean(stack, axis=0)
            tests.append(
                PreparedTest(
                    slug=_slug_for_parts(("lines", interval, "onesample", *key)),
                    stats_target="line-summaries",
                    test_type="one_sample",
                    variant=first.variant,
                    analysis_family=first.analysis_family,
                    sensor_group=first.sensor_group,
                    condition=first.condition,
                    contrast=f"{first.condition} vs 0",
                    model_source=first.model_source,
                    model_feature=first.model_feature,
                    stimulus_model_label=first.stimulus_model_label,
                    variable=first.variable,
                    slice_label=first.slice_label,
                    interval=interval,
                    subjects=matched_subjects,
                    dropped_subjects=tuple(subject for subject in subjects if subject not in by_subject),
                    observations=stack,
                    mean_map=mean_curve,
                    lag_sec=lag_axis,
                    source_paths=tuple(by_subject[subject].source_path for subject in matched_subjects),
                )
            )
    return tests


def _prepare_line_paired_tests(
    records: Sequence[MatrixRecord],
    subjects: Sequence[str],
    args: argparse.Namespace,
    cfg: Mapping[str, Any],
    *,
    skipped: list[dict[str, str]],
) -> list[PreparedTest]:
    """
    Assemble actual-position paired deviant-minus-nondeviant 1-D lag-curve tests.

    Each subject is differenced at the matrix level before Stage-08 lag
    aggregation. This preserves the paired finite-sample support for NaN-heavy
    matrices and keeps line-summary contrasts consistent with the 2-D paired
    matrix tests.
    """
    # First group the two condition-specific records by all non-condition
    # metadata so each paired contrast compares like with like.
    by_key_condition: dict[tuple[str, ...], dict[str, list[MatrixRecord]]] = {}
    for record in records:
        if not _record_matches_actual_position_contrast(record):
            continue
        by_key_condition.setdefault(_paired_group_key(record), {}).setdefault(
            record.condition,
            [],
        ).append(record)
    timing = occlusion_timing_from_config(cfg)

    tests: list[PreparedTest] = []
    for key, condition_items in sorted(by_key_condition.items()):
        deviant = _records_by_subject(condition_items.get("occluded_deviant", []))
        nondeviant = _records_by_subject(condition_items.get("occluded_nondeviant", []))
        matched_subjects = tuple(
            subject for subject in subjects if subject in deviant and subject in nondeviant
        )
        dropped = tuple(subject for subject in subjects if subject not in matched_subjects)
        if len(matched_subjects) < 2:
            skipped.append({"target": "line-summaries", "key": repr(key), "reason": "paired contrast has fewer than two matched subjects"})
            continue
        for interval in args.line_intervals:
            curves: list[np.ndarray] = []
            lag_axis: np.ndarray | None = None
            axis_mismatch = False
            shape_mismatch = False
            for subject in matched_subjects:
                dev_record = deviant[subject]
                non_record = nondeviant[subject]
                dev_matrix = np.asarray(dev_record.matrix, dtype=float)
                non_matrix = np.asarray(non_record.matrix, dtype=float)
                if dev_matrix.shape != non_matrix.shape:
                    shape_mismatch = True
                    break
                if (
                    dev_record.time_axis_sec.shape != non_record.time_axis_sec.shape
                    or not np.allclose(dev_record.time_axis_sec, non_record.time_axis_sec)
                ):
                    axis_mismatch = True
                    break

                # Matrix-level differencing is intentionally before lag
                # aggregation; nanmean(curve_a) - nanmean(curve_b) is not
                # equivalent when condition matrices have different NaN support.
                lag_sec, curve = _line_curve_for_matrix(
                    dev_matrix - non_matrix,
                    dev_record.time_axis_sec,
                    interval=interval,
                    occlusion_timing=timing,
                    lag_window_sec=args.line_lag_window_sec,
                    interval_lag_mode=args.line_interval_lag_mode,
                )
                if lag_axis is None:
                    lag_axis = lag_sec
                elif lag_sec.shape != lag_axis.shape or not np.allclose(lag_sec, lag_axis):
                    axis_mismatch = True
                    break
                curves.append(curve)
            if shape_mismatch:
                skipped.append({"target": "line-summaries", "key": repr(key), "interval": interval, "reason": "paired shape mismatch"})
                continue
            if axis_mismatch or lag_axis is None:
                skipped.append({"target": "line-summaries", "key": repr(key), "interval": interval, "reason": "paired lag-axis mismatch"})
                continue
            stack = np.vstack(curves)
            first = deviant[matched_subjects[0]]
            with np.errstate(invalid="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                mean_curve = np.nanmean(stack, axis=0)
            tests.append(
                PreparedTest(
                    slug=_slug_for_parts(("lines", interval, "actual-position-dev-minus-nondev", *key)),
                    stats_target="line-summaries",
                    test_type="paired",
                    variant=first.variant,
                    analysis_family=first.analysis_family,
                    sensor_group=first.sensor_group,
                    condition="occluded_deviant_minus_occluded_nondeviant",
                    contrast="actual deviant position - actual nondeviant position",
                    model_source=first.model_source,
                    model_feature=first.model_feature,
                    stimulus_model_label=first.stimulus_model_label,
                    variable=first.variable,
                    slice_label=first.slice_label,
                    interval=interval,
                    subjects=matched_subjects,
                    dropped_subjects=dropped,
                    observations=stack,
                    mean_map=mean_curve,
                    lag_sec=lag_axis,
                    source_paths=tuple(
                        path
                        for subject in matched_subjects
                        for path in (deviant[subject].source_path, nondeviant[subject].source_path)
                    ),
                )
            )
    return tests


def _one_sample_t_map(observations: Any) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute a NaN-safe one-sample t map over the first axis.

    Features with fewer than two finite observations return t=0 so they cannot
    form clusters. Zero-variance nonzero effects are assigned a large finite
    signed statistic instead of infinity, which keeps plotting and JSON output
    stable while preserving suprathreshold direction.
    """
    data = np.asarray(observations, dtype=float)
    finite = np.isfinite(data)
    counts = np.sum(finite, axis=0)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        means = np.nanmean(data, axis=0)
        stds = np.nanstd(data, axis=0, ddof=1)
    sems = stds / np.sqrt(np.maximum(counts, 1))
    t_map = np.zeros(data.shape[1:], dtype=float)
    valid = (counts > 1) & np.isfinite(sems) & (sems > 0)
    t_map[valid] = means[valid] / sems[valid]
    zero_var = (counts > 1) & (~valid) & np.isfinite(means) & (means != 0)
    t_map[zero_var] = np.sign(means[zero_var]) * 1.0e12
    return t_map, counts


def _cluster_forming_threshold(
    *,
    n_observations: int,
    cluster_alpha: float,
    tail: str,
) -> float:
    """Convert a cluster-forming p-value to an absolute t threshold."""
    try:
        from scipy import stats
    except ImportError as exc:  # pragma: no cover - scipy is pinned for DAD.
        raise RuntimeError("Stage-10 dRSA stats require scipy for t thresholds.") from exc

    if not 0 < cluster_alpha < 1:
        raise ValueError("cluster_alpha must be between 0 and 1.")
    df = int(n_observations) - 1
    if df < 1:
        raise ValueError("At least two observations are required for a t threshold.")
    quantile = 1.0 - (cluster_alpha / 2.0 if tail == "two-sided" else cluster_alpha)
    return float(stats.t.ppf(quantile, df))


def _cluster_statistic(t_map: np.ndarray, mask: np.ndarray, sign: str, tail: str) -> float:
    """Return a cluster-mass statistic for one suprathreshold cluster."""
    values = np.asarray(t_map, dtype=float)[mask]
    if tail == "two-sided":
        return float(np.sum(np.abs(values)))
    if sign == "negative":
        return float(np.sum(-values))
    return float(np.sum(values))


def _find_clusters(t_map: np.ndarray, threshold: float, tail: str) -> list[tuple[np.ndarray, str, float]]:
    """Find connected suprathreshold clusters in a 1-D or 2-D t map."""
    try:
        from scipy import ndimage
    except ImportError as exc:  # pragma: no cover - scipy is pinned for DAD.
        raise RuntimeError("Stage-10 dRSA stats require scipy.ndimage for clustering.") from exc

    array = np.asarray(t_map, dtype=float)
    structure = ndimage.generate_binary_structure(array.ndim, 1)
    signed_masks: list[tuple[np.ndarray, str]] = []
    if tail in {"two-sided", "positive"}:
        signed_masks.append((array >= threshold, "positive"))
    if tail in {"two-sided", "negative"}:
        signed_masks.append((array <= -threshold, "negative"))

    clusters: list[tuple[np.ndarray, str, float]] = []
    for supra_mask, sign in signed_masks:
        labels, n_labels = ndimage.label(supra_mask, structure=structure)
        for label in range(1, n_labels + 1):
            mask = labels == label
            if not np.any(mask):
                continue
            clusters.append((mask, sign, _cluster_statistic(array, mask, sign, tail)))
    return clusters


def _sign_flip_matrix(
    rng: np.random.Generator,
    *,
    n_observations: int,
    n_permutations: int,
) -> np.ndarray:
    """Generate exact or random sign-flip matrices for permutation testing."""
    exact_total = 2 ** int(n_observations)
    if n_permutations >= exact_total and n_observations <= 16:
        flips = np.asarray(
            list(itertools.product((-1.0, 1.0), repeat=int(n_observations))),
            dtype=float,
        )
        return flips
    return rng.choice(
        np.asarray([-1.0, 1.0], dtype=float),
        size=(int(n_permutations), int(n_observations)),
        replace=True,
    )


def _cluster_permutation_numpy(
    observations: Any,
    *,
    n_permutations: int,
    cluster_alpha: float,
    significance_alpha: float,
    tail: str,
    seed: int,
) -> ClusterTestResult:
    """Run a NaN-safe NumPy/SciPy one-sample sign-flip cluster test."""
    data = np.asarray(observations, dtype=float)
    if data.ndim not in {2, 3}:
        raise ValueError(f"Expected observations x 1-D/2-D features, found {data.shape}.")
    n_observations = data.shape[0]
    threshold = _cluster_forming_threshold(
        n_observations=n_observations,
        cluster_alpha=cluster_alpha,
        tail=tail,
    )
    observed_t, counts = _one_sample_t_map(data)
    observed_clusters = _find_clusters(observed_t, threshold, tail)

    rng = np.random.default_rng(seed)
    flips = _sign_flip_matrix(
        rng,
        n_observations=n_observations,
        n_permutations=n_permutations,
    )
    null_max = np.zeros(flips.shape[0], dtype=float)
    reshape = (flips.shape[0], n_observations) + (1,) * (data.ndim - 1)
    # Looping keeps memory bounded for 2-D matrices while preserving exact
    # sign-flip behavior and NaN support.
    for perm_index, signs in enumerate(flips.reshape(reshape)):
        perm_t, _counts = _one_sample_t_map(data * signs)
        perm_clusters = _find_clusters(perm_t, threshold, tail)
        null_max[perm_index] = max((cluster[2] for cluster in perm_clusters), default=0.0)

    clusters: list[ClusterInfo] = []
    for cluster_id, (mask, sign, statistic) in enumerate(observed_clusters, start=1):
        p_value = (1.0 + float(np.sum(null_max >= statistic))) / (1.0 + float(null_max.size))
        clusters.append(
            ClusterInfo(
                cluster_id=cluster_id,
                mask=mask,
                p_value=float(p_value),
                statistic=float(statistic),
                size=int(np.sum(mask)),
                sign=sign,
            )
        )
    clusters.sort(key=lambda item: (item.p_value, -item.statistic, item.cluster_id))
    clusters = [
        ClusterInfo(index + 1, item.mask, item.p_value, item.statistic, item.size, item.sign)
        for index, item in enumerate(clusters)
    ]
    return ClusterTestResult(
        statistic_map=observed_t,
        clusters=tuple(clusters),
        threshold=threshold,
        cluster_alpha=cluster_alpha,
        significance_alpha=significance_alpha,
        tail=tail,
        n_permutations=int(flips.shape[0]),
        backend="numpy",
        finite_observation_counts=counts,
        null_max_statistics=null_max,
    )


def _cluster_permutation_mne(
    observations: Any,
    *,
    n_permutations: int,
    cluster_alpha: float,
    significance_alpha: float,
    tail: str,
    seed: int,
) -> ClusterTestResult:
    """Run MNE's one-sample cluster permutation when finite support permits."""
    try:
        from mne.stats import permutation_cluster_1samp_test
    except ImportError as exc:
        raise RuntimeError("MNE is not available.") from exc

    data = np.asarray(observations, dtype=float)
    if not np.isfinite(data).all():
        raise RuntimeError("MNE backend requires complete finite support.")
    threshold = _cluster_forming_threshold(
        n_observations=data.shape[0],
        cluster_alpha=cluster_alpha,
        tail=tail,
    )
    tail_int = {"positive": 1, "negative": -1, "two-sided": 0}[tail]
    t_obs, clusters, cluster_p_values, null = permutation_cluster_1samp_test(
        data,
        n_permutations=int(n_permutations),
        threshold=threshold,
        tail=tail_int,
        seed=int(seed),
        out_type="mask",
        n_jobs=1,
        verbose=False,
    )
    _t_map, counts = _one_sample_t_map(data)
    cluster_infos: list[ClusterInfo] = []
    for cluster_id, (mask, p_value) in enumerate(zip(clusters, cluster_p_values), start=1):
        mask_array = _mne_cluster_to_mask(mask, np.asarray(t_obs).shape)
        sign = "negative" if float(np.nanmean(t_obs[mask_array])) < 0 else "positive"
        stat = _cluster_statistic(np.asarray(t_obs), mask_array, sign, tail)
        cluster_infos.append(
            ClusterInfo(
                cluster_id=cluster_id,
                mask=mask_array,
                p_value=float(p_value),
                statistic=stat,
                size=int(np.sum(mask_array)),
                sign=sign,
            )
        )
    cluster_infos.sort(key=lambda item: (item.p_value, -item.statistic, item.cluster_id))
    cluster_infos = [
        ClusterInfo(index + 1, item.mask, item.p_value, item.statistic, item.size, item.sign)
        for index, item in enumerate(cluster_infos)
    ]
    return ClusterTestResult(
        statistic_map=np.asarray(t_obs, dtype=float),
        clusters=tuple(cluster_infos),
        threshold=threshold,
        cluster_alpha=cluster_alpha,
        significance_alpha=significance_alpha,
        tail=tail,
        n_permutations=int(n_permutations),
        backend="mne",
        finite_observation_counts=counts,
        null_max_statistics=np.asarray(null, dtype=float),
    )


def _mne_cluster_to_mask(cluster: Any, feature_shape: Sequence[int]) -> np.ndarray:
    """
    Normalize one MNE cluster object to a boolean feature mask.

    MNE can return boolean arrays for 2-D feature grids but tuple-of-slice or
    tuple-of-index clusters for 1-D features, even when ``out_type='mask'``.
    Stage-10 report code expects a mask matching the statistic map exactly.
    """
    shape = tuple(int(value) for value in feature_shape)
    if isinstance(cluster, tuple):
        mask = np.zeros(shape, dtype=bool)
        mask[cluster] = True
        return mask

    array = np.asarray(cluster)
    if array.dtype == bool:
        squeezed = np.squeeze(array)
        if squeezed.shape == shape:
            return squeezed.astype(bool, copy=False)
        if squeezed.size == int(np.prod(shape)):
            return squeezed.reshape(shape).astype(bool, copy=False)
    raise RuntimeError(
        f"Unsupported MNE cluster mask shape {array.shape}; expected {shape}."
    )


def _mne_auto_backend_is_safe(data: np.ndarray) -> bool:
    """
    Return whether the auto backend should use MNE for this observation stack.

    MNE's default one-sample t statistic emits infinities/NaNs when any feature
    has zero variance across subjects. The local NumPy backend handles those
    features deterministically, so ``auto`` uses MNE only for fully finite,
    non-constant feature grids.
    """
    if not np.isfinite(data).all():
        return False
    if data.shape[0] < 2:
        return False
    variances = np.var(data, axis=0, ddof=1)
    return bool(np.isfinite(variances).all() and np.all(variances > 0))


def cluster_permutation_1samp(
    observations: Any,
    *,
    n_permutations: int,
    cluster_alpha: float,
    significance_alpha: float,
    tail: str,
    seed: int,
    backend: str = "auto",
) -> ClusterTestResult:
    """
    Run a one-sample cluster-permutation test over 1-D or 2-D features.

    The MNE backend is preferred in ``auto`` mode when all observations are
    finite. NaN-containing arrays use the local NumPy/SciPy implementation
    because MNE's cluster functions do not accept partial finite support.
    """
    data = np.asarray(observations, dtype=float)
    if data.shape[0] < 2:
        raise ValueError("At least two observations are required.")
    if tail not in TAIL_CHOICES:
        raise ValueError(f"tail must be one of {TAIL_CHOICES}.")
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1.")

    use_mne = backend == "mne" or (backend == "auto" and _mne_auto_backend_is_safe(data))
    if use_mne:
        try:
            return _cluster_permutation_mne(
                data,
                n_permutations=n_permutations,
                cluster_alpha=cluster_alpha,
                significance_alpha=significance_alpha,
                tail=tail,
                seed=seed,
            )
        except RuntimeError:
            if backend == "mne":
                raise
    return _cluster_permutation_numpy(
        data,
        n_permutations=n_permutations,
        cluster_alpha=cluster_alpha,
        significance_alpha=significance_alpha,
        tail=tail,
        seed=seed,
    )


def _significant_mask(result: ClusterTestResult) -> np.ndarray:
    """Return the union of clusters meeting the report significance alpha."""
    mask = np.zeros_like(result.statistic_map, dtype=bool)
    for cluster in result.clusters:
        if cluster.p_value <= result.significance_alpha:
            mask |= cluster.mask
    return mask


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


def _json_safe(value: Any) -> Any:
    """Convert NumPy scalars/arrays into JSON-serializable Python values."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _write_json(path: Path, payload: Mapping[str, Any], *, raw_root: Path, args: argparse.Namespace) -> Path:
    """Write JSON unless resume mode reuses an existing file."""
    target = _prepare_output_path(path, raw_root=raw_root, args=args)
    if target is not None:
        target.write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True), encoding="utf-8")
    return path


def _cluster_extent_fields(
    cluster: ClusterInfo,
    *,
    test: PreparedTest,
) -> dict[str, Any]:
    """Describe cluster index and time/lag extents for CSV and HTML tables."""
    indices = np.argwhere(cluster.mask)
    fields: dict[str, Any] = {
        "cluster_id": cluster.cluster_id,
        "p_value": cluster.p_value,
        "cluster_statistic": cluster.statistic,
        "cluster_size": cluster.size,
        "cluster_sign": cluster.sign,
    }
    if indices.size == 0:
        return fields
    if test.stats_target == "matrices":
        rows = indices[:, 0]
        cols = indices[:, 1]
        axis = np.asarray(test.time_axis_sec, dtype=float)
        fields.update(
            {
                "model_index_min": int(rows.min()),
                "model_index_max": int(rows.max()),
                "neural_index_min": int(cols.min()),
                "neural_index_max": int(cols.max()),
                "model_time_min_sec": float(axis[rows.min()]),
                "model_time_max_sec": float(axis[rows.max()]),
                "neural_time_min_sec": float(axis[cols.min()]),
                "neural_time_max_sec": float(axis[cols.max()]),
            }
        )
    else:
        lag_indices = indices[:, 0]
        lag_axis = np.asarray(test.lag_sec, dtype=float)
        fields.update(
            {
                "lag_index_min": int(lag_indices.min()),
                "lag_index_max": int(lag_indices.max()),
                "lag_min_sec": float(lag_axis[lag_indices.min()]),
                "lag_max_sec": float(lag_axis[lag_indices.max()]),
            }
        )
    return fields


def _cluster_rows(test: PreparedTest, result: ClusterTestResult) -> list[dict[str, Any]]:
    """Build CSV/HTML cluster rows for one completed test."""
    rows = []
    for cluster in result.clusters:
        row = {
            "stats_target": test.stats_target,
            "test_type": test.test_type,
            "variant": test.variant,
            "analysis_family": test.analysis_family,
            "sensor_group": test.sensor_group,
            "condition": test.condition,
            "contrast": test.contrast,
            "model_source": test.model_source,
            "model_feature": test.model_feature,
            "stimulus_model_label": test.stimulus_model_label,
            "variable": test.variable,
            "slice_label": test.slice_label,
            "interval": test.interval or "",
            "n_subjects": len(test.subjects),
            "subjects": " ".join(test.subjects),
            "dropped_subjects": " ".join(test.dropped_subjects),
            "backend": result.backend,
            "tail": result.tail,
            "cluster_forming_alpha": result.cluster_alpha,
            "cluster_forming_t_threshold": result.threshold,
            "cluster_p_significance_alpha": result.significance_alpha,
            "n_permutations": result.n_permutations,
        }
        row.update(_cluster_extent_fields(cluster, test=test))
        rows.append(row)
    if not rows:
        rows.append(
            {
                "stats_target": test.stats_target,
                "test_type": test.test_type,
                "variant": test.variant,
                "analysis_family": test.analysis_family,
                "sensor_group": test.sensor_group,
                "condition": test.condition,
                "contrast": test.contrast,
                "model_source": test.model_source,
                "model_feature": test.model_feature,
                "stimulus_model_label": test.stimulus_model_label,
                "variable": test.variable,
                "slice_label": test.slice_label,
                "interval": test.interval or "",
                "n_subjects": len(test.subjects),
                "subjects": " ".join(test.subjects),
                "dropped_subjects": " ".join(test.dropped_subjects),
                "backend": result.backend,
                "tail": result.tail,
                "cluster_forming_alpha": result.cluster_alpha,
                "cluster_forming_t_threshold": result.threshold,
                "cluster_p_significance_alpha": result.significance_alpha,
                "n_permutations": result.n_permutations,
                "cluster_id": "",
                "p_value": "",
                "cluster_statistic": "",
                "cluster_size": 0,
                "cluster_sign": "",
            }
        )
    return rows


def _write_cluster_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    raw_root: Path,
    args: argparse.Namespace,
) -> Path:
    """Write one cluster table CSV unless resume mode reuses it."""
    target = _prepare_output_path(path, raw_root=raw_root, args=args)
    if target is None:
        return path
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with target.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def _save_npz(
    path: Path,
    *,
    test: PreparedTest,
    result: ClusterTestResult,
    raw_root: Path,
    args: argparse.Namespace,
) -> Path:
    """Save machine-readable arrays for one statistical test."""
    target = _prepare_output_path(path, raw_root=raw_root, args=args)
    if target is None:
        return path
    cluster_masks = (
        np.stack([cluster.mask for cluster in result.clusters], axis=0)
        if result.clusters
        else np.zeros((0,) + result.statistic_map.shape, dtype=bool)
    )
    payload = {
        "observations": test.observations,
        "mean_map": test.mean_map,
        "statistic_map": result.statistic_map,
        "finite_observation_counts": result.finite_observation_counts,
        "cluster_masks": cluster_masks,
        "cluster_p_values": np.asarray([cluster.p_value for cluster in result.clusters], dtype=float),
        "cluster_statistics": np.asarray([cluster.statistic for cluster in result.clusters], dtype=float),
        "null_max_statistics": result.null_max_statistics,
    }
    if test.time_axis_sec is not None:
        payload["time_axis_sec"] = test.time_axis_sec
    if test.lag_sec is not None:
        payload["lag_sec"] = test.lag_sec
    np.savez_compressed(target, **payload)
    return path


def _plot_matrix_result(
    path: Path,
    *,
    test: PreparedTest,
    result: ClusterTestResult,
    raw_root: Path,
    args: argparse.Namespace,
) -> Path:
    """Write a 2-D mean/stat/mask summary figure for one matrix test."""
    target = _prepare_output_path(path, raw_root=raw_root, args=args)
    if target is None:
        return path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    axis = np.asarray(test.time_axis_sec, dtype=float)
    extent = [float(axis[0]), float(axis[-1]), float(axis[0]), float(axis[-1])]
    sig_mask = _significant_mask(result)
    fig, axes = plt.subplots(ncols=3, figsize=(13.0, 4.2), constrained_layout=True)
    images = [
        axes[0].imshow(test.mean_map, origin="lower", aspect="auto", extent=extent),
        axes[1].imshow(result.statistic_map, origin="lower", aspect="auto", extent=extent, cmap="coolwarm"),
        axes[2].imshow(test.mean_map, origin="lower", aspect="auto", extent=extent),
    ]
    if np.any(sig_mask):
        axes[2].contour(
            axis,
            axis,
            sig_mask.astype(float),
            levels=[0.5],
            colors="black",
            linewidths=1.0,
        )
    titles = ["Mean input", "t statistic", f"p <= {result.significance_alpha:g} clusters"]
    for ax, title, image in zip(axes, titles, images):
        ax.set_title(title)
        ax.set_xlabel("Neural time (s)")
        ax.set_ylabel("Model time (s)")
        fig.colorbar(image, ax=ax, shrink=0.78)
    fig.suptitle(
        _compact_title(
            [
                test.variant,
                test.analysis_family,
                test.sensor_group,
                test.contrast,
                test.stimulus_model_label,
                test.variable,
            ]
        ),
        fontsize=10,
    )
    fig.savefig(target, dpi=180, facecolor="white")
    plt.close(fig)
    return path


def _plot_line_result(
    path: Path,
    *,
    test: PreparedTest,
    result: ClusterTestResult,
    raw_root: Path,
    args: argparse.Namespace,
) -> Path:
    """Write a 1-D lag-curve test figure with significant cluster bands."""
    target = _prepare_output_path(path, raw_root=raw_root, args=args)
    if target is None:
        return path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lag_sec = np.asarray(test.lag_sec, dtype=float)
    observations = np.asarray(test.observations, dtype=float)
    finite_counts = np.sum(np.isfinite(observations), axis=0)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sem = np.nanstd(observations, axis=0, ddof=1) / np.sqrt(np.maximum(finite_counts, 1))
    mean_curve = np.asarray(test.mean_map, dtype=float)
    sig_mask = _significant_mask(result)

    fig, (ax, stat_ax) = plt.subplots(
        nrows=2,
        figsize=(8.5, 5.5),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )
    ax.plot(lag_sec, mean_curve, color="tab:blue", linewidth=1.8, label="group mean")
    finite_band = np.isfinite(mean_curve) & np.isfinite(sem)
    if np.any(finite_band):
        ax.fill_between(
            lag_sec,
            mean_curve - sem,
            mean_curve + sem,
            where=finite_band,
            color="tab:blue",
            alpha=0.18,
            linewidth=0,
            label="SEM",
        )
    ax.axhline(0.0, color="0.4", linestyle=":", linewidth=0.9)
    ax.axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_ylabel("Mean dRSA value")
    ax.legend(loc="best", fontsize=8)
    ax.set_title(
        _compact_title(
            [
                test.variant,
                test.interval,
                test.analysis_family,
                test.sensor_group,
                test.contrast,
                test.stimulus_model_label,
            ]
        ),
        fontsize=10,
    )

    stat_ax.plot(lag_sec, result.statistic_map, color="tab:red", linewidth=1.4)
    stat_ax.axhline(result.threshold, color="0.35", linestyle=":", linewidth=0.9)
    stat_ax.axhline(-result.threshold, color="0.35", linestyle=":", linewidth=0.9)
    if np.any(sig_mask):
        ymin, ymax = stat_ax.get_ylim()
        for start, stop in _contiguous_true_runs(sig_mask):
            stat_ax.axvspan(
                lag_sec[start],
                lag_sec[stop],
                color="gold",
                alpha=0.28,
                linewidth=0,
            )
        stat_ax.set_ylim(ymin, ymax)
    stat_ax.set_xlabel("Lag (neural time - model time) [s]")
    stat_ax.set_ylabel("t statistic")
    fig.savefig(target, dpi=180, facecolor="white")
    plt.close(fig)
    return path


def _contiguous_true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return inclusive index spans for contiguous True values in a 1-D mask."""
    values = np.asarray(mask, dtype=bool).ravel()
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, value in enumerate(values):
        if value and start is None:
            start = index
        elif not value and start is not None:
            runs.append((start, index - 1))
            start = None
    if start is not None:
        runs.append((start, len(values) - 1))
    return runs


def _compact_title(parts: Iterable[Any], *, width: int = 95) -> str:
    """Return a wrapped, readable plot title from metadata parts."""
    text = " | ".join(str(part) for part in parts if str(part).strip())
    return "\n".join(textwrap.wrap(text, width=width)) if len(text) > width else text


def _test_metadata(test: PreparedTest, result: ClusterTestResult) -> dict[str, Any]:
    """Build JSON metadata for one completed statistical test."""
    return {
        "slug": test.slug,
        "stats_target": test.stats_target,
        "test_type": test.test_type,
        "variant": test.variant,
        "analysis_family": test.analysis_family,
        "sensor_group": test.sensor_group,
        "condition": test.condition,
        "contrast": test.contrast,
        "model_source": test.model_source,
        "model_feature": test.model_feature,
        "stimulus_model_label": test.stimulus_model_label,
        "variable": test.variable,
        "slice_label": test.slice_label,
        "interval": test.interval,
        "subjects": list(test.subjects),
        "dropped_subjects": list(test.dropped_subjects),
        "n_subjects": len(test.subjects),
        "source_paths": list(test.source_paths),
        "observation_shape": list(test.observations.shape),
        "backend": result.backend,
        "tail": result.tail,
        "cluster_forming_alpha": result.cluster_alpha,
        "cluster_forming_t_threshold": result.threshold,
        "cluster_p_significance_alpha": result.significance_alpha,
        "n_permutations": result.n_permutations,
        "cluster_count": len(result.clusters),
        "significant_cluster_count": sum(
            1 for cluster in result.clusters if cluster.p_value <= result.significance_alpha
        ),
        "finite_feature_count": int(np.sum(np.any(np.isfinite(test.observations), axis=0))),
    }


def _run_and_write_test(
    test: PreparedTest,
    *,
    output_dir: Path,
    raw_root: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Run one prepared test and write arrays, figure, JSON, and CSV outputs."""
    result = cluster_permutation_1samp(
        test.observations,
        n_permutations=args.n_permutations,
        cluster_alpha=args.cluster_alpha,
        significance_alpha=args.significance_alpha,
        tail=args.tail,
        seed=args.seed,
        backend=args.stats_backend,
    )
    target_dir = output_dir / test.stats_target / test.test_type
    figure_path = target_dir / "figures" / f"{test.slug}.png"
    npz_path = target_dir / "arrays" / f"{test.slug}.npz"
    json_path = target_dir / "json" / f"{test.slug}.json"
    csv_path = target_dir / "tables" / f"{test.slug}_clusters.csv"

    if test.stats_target == "matrices":
        _plot_matrix_result(figure_path, test=test, result=result, raw_root=raw_root, args=args)
    else:
        _plot_line_result(figure_path, test=test, result=result, raw_root=raw_root, args=args)
    _save_npz(npz_path, test=test, result=result, raw_root=raw_root, args=args)
    rows = _cluster_rows(test, result)
    _write_cluster_csv(csv_path, rows, raw_root=raw_root, args=args)
    _write_json(
        json_path,
        {
            "metadata": _test_metadata(test, result),
            "clusters": rows,
        },
        raw_root=raw_root,
        args=args,
    )

    return {
        "metadata": _test_metadata(test, result),
        "figure": str(figure_path),
        "npz": str(npz_path),
        "json": str(json_path),
        "cluster_csv": str(csv_path),
        "clusters": rows,
    }


def _write_summary_csv(
    path: Path,
    entries: Sequence[Mapping[str, Any]],
    *,
    raw_root: Path,
    args: argparse.Namespace,
) -> Path:
    """Write the combined machine-readable cluster summary CSV."""
    rows = [row for entry in entries for row in entry.get("clusters", [])]
    return _write_cluster_csv(path, rows, raw_root=raw_root, args=args)


def _write_report(
    report_path: Path,
    *,
    title: str,
    summary_lines: Sequence[str],
    entries: Sequence[Mapping[str, Any]],
    skipped_inputs: Sequence[Mapping[str, str]],
    raw_root: Path,
    report_cluster_tables: str,
) -> Path:
    """
    Write a static HTML report and embedded sibling copy.

    The report is organized by variant, condition, sensor group, and stats
    target, followed by the same hierarchy for contrasts. The nested index links
    to every rendered level and every individual test figure.
    """
    report = guarded_path(report_path, raw_root)
    report.parent.mkdir(parents=True, exist_ok=True)
    nav_items, sections = _organized_report_sections(
        report,
        entries=entries,
        report_cluster_tables=report_cluster_tables,
    )
    toc = _nav_tree_html(nav_items)

    skipped_block = ""
    if skipped_inputs:
        rows = "\n".join(
            "<li>"
            + html.escape(" | ".join(f"{key}: {value}" for key, value in sorted(item.items())))
            + "</li>"
            for item in skipped_inputs
        )
        skipped_block = f"<details><summary>Skipped inputs and tests</summary><ul>{rows}</ul></details>"

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
            "table{border-collapse:collapse;font-size:.82rem;margin:.7rem 0;}",
            "td,th{border:1px solid #ddd;padding:.25rem .4rem;text-align:left;}",
            "h2{margin-top:3rem;border-top:3px solid #222;padding-top:1rem;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<pre>{html.escape(chr(10).join(summary_lines))}</pre>",
            "<nav><strong>Index</strong>",
            toc,
            "</nav>",
            skipped_block,
            "\n".join(sections),
            "</body></html>",
        ]
    )
    report.write_text(body, encoding="utf-8")
    write_embedded_html_copy(report, raw_root=raw_root)
    return report


def _organized_report_sections(
    report: Path,
    *,
    entries: Sequence[Mapping[str, Any]],
    report_cluster_tables: str,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Build nested navigation and HTML sections for the Stage-10 report."""
    nav_items: list[dict[str, Any]] = []
    sections: list[str] = []
    variants = sorted(
        {str(entry["metadata"]["variant"]) for entry in entries},
        key=_variant_sort_key,
    )
    for variant in variants:
        variant_entries = [
            entry for entry in entries if str(entry["metadata"]["variant"]) == variant
        ]
        variant_anchor = html_id(f"variant-{variant}")
        variant_children: list[dict[str, Any]] = []
        variant_blocks: list[str] = []

        for condition in CONDITIONS:
            condition_entries = [
                entry
                for entry in variant_entries
                if entry["metadata"]["test_type"] == "one_sample"
                and entry["metadata"]["condition"] == condition
            ]
            if not condition_entries:
                continue
            label = _condition_label(condition)
            anchor = html_id(f"{variant}-condition-{condition}")
            children, html_block = _report_group_block(
                report,
                condition_entries,
                anchor=anchor,
                heading_level=3,
                title=label,
                report_cluster_tables=report_cluster_tables,
            )
            variant_children.append({"label": label, "anchor": anchor, "children": children})
            variant_blocks.append(html_block)

        contrast_entries = [
            entry
            for entry in variant_entries
            if entry["metadata"]["test_type"] != "one_sample"
            or entry["metadata"]["condition"] not in CONDITIONS
        ]
        if contrast_entries:
            anchor = html_id(f"{variant}-contrasts")
            children, html_block = _report_group_block(
                report,
                contrast_entries,
                anchor=anchor,
                heading_level=3,
                title="Contrasts",
                report_cluster_tables=report_cluster_tables,
            )
            variant_children.append({"label": "Contrasts", "anchor": anchor, "children": children})
            variant_blocks.append(html_block)

        nav_items.append(
            {
                "label": _variant_label(variant),
                "anchor": variant_anchor,
                "children": variant_children,
            }
        )
        sections.append(
            f"<section id=\"{variant_anchor}\"><h2>{html.escape(_variant_label(variant))}</h2>"
            + "\n".join(variant_blocks or ["<p>No tests for this variant.</p>"])
            + "</section>"
        )
    return nav_items, sections


def _report_group_block(
    report: Path,
    entries: Sequence[Mapping[str, Any]],
    *,
    anchor: str,
    heading_level: int,
    title: str,
    report_cluster_tables: str,
) -> tuple[list[dict[str, Any]], str]:
    """Render one condition or contrast group by sensor group and stats target."""
    nav_children: list[dict[str, Any]] = []
    blocks: list[str] = []
    sensor_groups = sorted(
        {str(entry["metadata"]["sensor_group"]) for entry in entries},
        key=_sensor_sort_key,
    )
    for sensor_group in sensor_groups:
        sensor_entries = [
            entry for entry in entries if str(entry["metadata"]["sensor_group"]) == sensor_group
        ]
        sensor_anchor = html_id(f"{anchor}-sensor-{sensor_group}")
        target_children: list[dict[str, Any]] = []
        target_blocks: list[str] = []
        for target in STATS_TARGETS:
            target_entries = [
                entry
                for entry in sensor_entries
                if str(entry["metadata"]["stats_target"]) == target
            ]
            if not target_entries:
                continue
            target_anchor = html_id(f"{sensor_anchor}-target-{target}")
            figure_children: list[dict[str, Any]] = []
            figures: list[str] = []
            for entry in sorted(target_entries, key=_entry_sort_key):
                figure_anchor = html_id(f"test-{entry['metadata']['slug']}")
                figure_children.append(
                    {
                        "label": _entry_nav_label(entry["metadata"]),
                        "anchor": figure_anchor,
                        "children": [],
                    }
                )
                figures.append(
                    _entry_figure_html(
                        report,
                        entry,
                        anchor=figure_anchor,
                        report_cluster_tables=report_cluster_tables,
                    )
                )
            target_children.append(
                {
                    "label": _target_label(target),
                    "anchor": target_anchor,
                    "children": figure_children,
                }
            )
            target_blocks.append(
                f"<section id=\"{target_anchor}\"><h{heading_level + 2}>{html.escape(_target_label(target))}</h{heading_level + 2}>"
                + "\n".join(figures)
                + "</section>"
            )
        nav_children.append(
            {
                "label": sensor_group,
                "anchor": sensor_anchor,
                "children": target_children,
            }
        )
        blocks.append(
            f"<section id=\"{sensor_anchor}\"><h{heading_level + 1}>Sensor group: {html.escape(sensor_group)}</h{heading_level + 1}>"
            + "\n".join(target_blocks or ["<p>No tests for this sensor group.</p>"])
            + "</section>"
        )
    html_block = (
        f"<section id=\"{anchor}\"><h{heading_level}>{html.escape(title)}</h{heading_level}>"
        + "\n".join(blocks or ["<p>No tests.</p>"])
        + "</section>"
    )
    return nav_children, html_block


def _entry_figure_html(
    report: Path,
    entry: Mapping[str, Any],
    *,
    anchor: str,
    report_cluster_tables: str,
) -> str:
    """Render one test figure and its significant-cluster table."""
    metadata = entry["metadata"]
    href = image_href(report, entry["figure"])
    significant = int(metadata["significant_cluster_count"])
    caption = (
        f"target={metadata['stats_target']} | test={metadata['test_type']} | "
        f"variant={metadata['variant']} | family={metadata['analysis_family']} | "
        f"sensor={metadata['sensor_group']} | condition/contrast={metadata['contrast']} | "
        f"model={metadata['stimulus_model_label']} | source={metadata['model_source']} | "
        f"feature={metadata['model_feature']} | variable={metadata['variable']} | "
        f"slice={metadata['slice_label']} | interval={metadata.get('interval') or 'n/a'} | "
        f"n={metadata['n_subjects']} | backend={metadata['backend']} | "
        f"tail={metadata['tail']} | permutations={metadata['n_permutations']} | "
        f"cluster-forming p={metadata['cluster_forming_alpha']} "
        f"(t={metadata['cluster_forming_t_threshold']:.3g}) | "
        f"significant clusters={significant}"
    )
    table_html = _cluster_table_html(
        entry.get("clusters", []),
        mode=report_cluster_tables,
        significance_alpha=float(metadata["cluster_p_significance_alpha"]),
    )
    return (
        f"<figure id=\"{anchor}\">"
        f"<img src=\"{html.escape(href)}\"/>"
        f"<figcaption>{html.escape(caption)}</figcaption>"
        f"{table_html}"
        "</figure>"
    )


def _nav_tree_html(items: Sequence[Mapping[str, Any]]) -> str:
    """Render a nested clickable index for report sections and figures."""
    if not items:
        return "<p>No report sections.</p>"
    rows = []
    for item in items:
        children = _nav_tree_html(item.get("children", [])) if item.get("children") else ""
        rows.append(
            "<li>"
            f"<a href=\"#{html.escape(str(item['anchor']))}\">{html.escape(str(item['label']))}</a>"
            f"{children}"
            "</li>"
        )
    return "<ul>" + "\n".join(rows) + "</ul>"


def _variant_sort_key(variant: str) -> tuple[int, str]:
    """Sort raw first, smoothed-display second, and any future variants last."""
    order = {VARIANT_RAW: 0, VARIANT_SMOOTHED: 1}
    return (order.get(variant, 99), variant)


def _variant_label(variant: str) -> str:
    """Return a reader-facing variant label."""
    if variant == VARIANT_RAW:
        return "Raw inputs"
    if variant == VARIANT_SMOOTHED:
        return "Smoothed-display inputs"
    return str(variant)


def _condition_label(condition: str) -> str:
    """Return a reader-facing condition heading."""
    labels = {
        "always_visible": "Always visible",
        "occluded_nondeviant": "Occluded nondeviant",
        "occluded_deviant": "Occluded deviant",
    }
    return labels.get(condition, condition)


def _sensor_sort_key(sensor_group: str) -> tuple[int, str]:
    """Sort common sensor groups in a stable anatomical/report order."""
    order = {
        "wholebrain": 0,
        "occipital": 1,
        "parietal": 2,
        "frontal": 3,
        "temporal": 4,
        "model": 90,
    }
    return (order.get(sensor_group, 50), sensor_group)


def _target_label(target: str) -> str:
    """Return a reader-facing stats target label."""
    return "Line summaries" if target == "line-summaries" else "Matrices"


def _entry_sort_key(entry: Mapping[str, Any]) -> tuple[Any, ...]:
    """Sort tests within one target subsection."""
    metadata = entry["metadata"]
    interval_order = {value: index for index, value in enumerate(DEFAULT_LINE_INTERVALS)}
    return (
        str(metadata.get("analysis_family", "")),
        str(metadata.get("model_source", "")),
        str(metadata.get("model_feature", "")),
        str(metadata.get("stimulus_model_label", "")),
        str(metadata.get("variable", "")),
        str(metadata.get("slice_label", "")),
        interval_order.get(str(metadata.get("interval", "")), 99),
        str(metadata.get("interval", "")),
    )


def _entry_nav_label(metadata: Mapping[str, Any]) -> str:
    """Return a concise label for one test figure in the clickable index."""
    parts = [
        metadata.get("analysis_family"),
        metadata.get("stimulus_model_label"),
        metadata.get("model_source"),
        metadata.get("interval") or metadata.get("stats_target"),
    ]
    return " | ".join(str(part) for part in parts if str(part).strip())


def _cluster_table_html(
    rows: Sequence[Mapping[str, Any]],
    *,
    mode: str,
    significance_alpha: float,
) -> str:
    """Render a compact per-figure cluster table for the selected HTML mode."""
    if mode == "none":
        return "<p>Cluster tables are hidden in this report; see the CSV/JSON outputs.</p>"
    cluster_rows = [row for row in rows if row.get("cluster_id")]
    if not cluster_rows:
        return "<p>No suprathreshold clusters.</p>"
    if mode == "significant":
        cluster_rows = [
            row
            for row in cluster_rows
            if row.get("p_value") not in {"", None}
            and float(row.get("p_value")) <= float(significance_alpha)
        ]
        if not cluster_rows:
            return (
                "<p>No significant clusters at the configured cluster p-value "
                f"threshold ({significance_alpha:g}).</p>"
            )
    visible_keys = [
        "cluster_id",
        "p_value",
        "cluster_size",
        "cluster_sign",
        "model_time_min_sec",
        "model_time_max_sec",
        "neural_time_min_sec",
        "neural_time_max_sec",
        "lag_min_sec",
        "lag_max_sec",
    ]
    header = "".join(f"<th>{html.escape(key)}</th>" for key in visible_keys)
    body_rows = []
    for row in cluster_rows:
        cells = []
        for key in visible_keys:
            value = row.get(key, "")
            if isinstance(value, float):
                value = f"{value:.4g}"
            cells.append(f"<td>{html.escape(str(value))}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return "<table><thead><tr>" + header + "</tr></thead><tbody>" + "".join(body_rows) + "</tbody></table>"


def run_group_scope(
    base_cfg: Mapping[str, Any],
    args: argparse.Namespace,
    subjects: Sequence[str],
) -> dict[str, Any]:
    """Read subject dRSA outputs, run selected stats, and write group reports."""
    raw_root = Path(base_cfg["resolved_paths"]["data_root"])
    paths = _analysis_paths(base_cfg)
    stats_dir = guarded_path(paths["group_stats_dir"], raw_root)
    stats_dir.mkdir(parents=True, exist_ok=True)

    all_records: list[MatrixRecord] = []
    skipped: list[dict[str, str]] = []
    for subject in subjects:
        subject_cfg = _load_cfg(args, subject_config=_subject_config_for_id(args, subject))
        records, subject_skipped = _read_subject_records(subject_cfg, args, subject=subject)
        all_records.extend(records)
        skipped.extend(subject_skipped)

    if not all_records:
        raise RuntimeError("No selected Stage-07 dRSA matrices were found.")

    prepared: list[PreparedTest] = []
    if "matrices" in args.stats_targets:
        prepared.extend(_prepare_matrix_one_sample_tests(all_records, subjects, skipped=skipped))
        if _paired_deviant_nondeviant_enabled(args):
            prepared.extend(_prepare_matrix_paired_tests(all_records, subjects, skipped=skipped))
    if "line-summaries" in args.stats_targets:
        prepared.extend(_prepare_line_one_sample_tests(all_records, subjects, args, base_cfg, skipped=skipped))
        if _paired_deviant_nondeviant_enabled(args):
            prepared.extend(_prepare_line_paired_tests(all_records, subjects, args, base_cfg, skipped=skipped))
    if not prepared:
        raise RuntimeError(f"No statistical tests were prepared. Skipped: {skipped[:5]}")

    entries = [
        _run_and_write_test(test, output_dir=stats_dir, raw_root=raw_root, args=args)
        for test in prepared
    ]
    summary_csv = _write_summary_csv(
        stats_dir / "10_drsa_stats_cluster_summary.csv",
        entries,
        raw_root=raw_root,
        args=args,
    )
    summary_json = _write_json(
        stats_dir / "10_drsa_stats_summary.json",
        {
            "stage": "drsa_stats",
            "generated_utc": utc_now(),
            "subjects": list(subjects),
            "record_count": len(all_records),
            "test_count": len(prepared),
            "analysis_families": list(args.analysis_families or []),
            "paired_contrasts": list(args.paired_contrasts or []),
            "report_cluster_tables": args.report_cluster_tables,
            "entries": [entry["metadata"] for entry in entries],
            "skipped_inputs": skipped,
            "summary_csv": str(summary_csv),
        },
        raw_root=raw_root,
        args=args,
    )

    first_sec = occlusion_timing_from_config(base_cfg)["first_complete_occlusion_sec"]
    reapp_sec = occlusion_timing_from_config(base_cfg)["reappearance_sec"]
    variants = sorted({record.variant for record in all_records})
    summary_lines = [
        "Stage status: completed",
        f"Generated UTC: {utc_now()}",
        f"Subjects: {' '.join(subjects)}",
        f"Input matrix records: {len(all_records)}",
        f"Statistical tests: {len(entries)}",
        f"Stats targets: {', '.join(args.stats_targets)}",
        f"Analysis families: {', '.join(args.analysis_families or [])}",
        f"Paired contrasts: {', '.join(args.paired_contrasts or [])}",
        f"HTML cluster tables: {args.report_cluster_tables}",
        f"Variants: {', '.join(variants)} (raw = Stage-07 MAT matrices, "
        "possibly MATLAB-neural-smoothed; smoothed_display only when requested)",
        "Matrix interpretation: x = neural time, y = model time, indexed matrix[model_time, neural_time]",
        "Lag convention: lag = neural time - model time",
        f"Line interval lag mode: {args.line_interval_lag_mode}",
        (
            "Model-anchored interval behavior: interval masks restrict model-time "
            "rows; neural-time columns may fall anywhere inside the full matrix."
        ),
        (
            "Occlusion intervals on model time: "
            f"pre_occlusion [0, {first_sec:.6f}); "
            f"during_occlusion [{first_sec:.6f}, {reapp_sec:.6f}); "
            f"post_occlusion [{reapp_sec:.6f}, end]"
        ),
        f"Line lag window: +/-{args.line_lag_window_sec:g} s",
        f"Cluster-forming alpha: {args.cluster_alpha:g}",
        f"Cluster p-value significance alpha: {args.significance_alpha:g}",
        f"Tail: {args.tail}",
        f"Requested permutations: {args.n_permutations}",
        f"Seed: {args.seed}",
        f"Stats backend: {args.stats_backend}",
        f"Machine summary JSON: {summary_json}",
        f"Machine summary CSV: {summary_csv}",
        f"Output directory: {stats_dir}",
        f"Skipped inputs/tests: {len(skipped)}",
    ]
    report_html = _write_report(
        paths["group_report_html"],
        title=f"Group | {base_cfg['layout']['session_bids_id']} | 10 dRSA stats",
        summary_lines=summary_lines,
        entries=entries,
        skipped_inputs=skipped,
        raw_root=raw_root,
        report_cluster_tables=args.report_cluster_tables,
    )
    return {
        "status": "completed",
        "record_count": len(all_records),
        "test_count": len(entries),
        "output_dir": str(stats_dir),
        "summary_json": str(summary_json),
        "summary_csv": str(summary_csv),
        "report_html": str(report_html),
        "report_embedded_html": str(report_html.with_name(f"{report_html.stem}_embedded.html")),
        "skipped_inputs": skipped,
    }


def main() -> int:
    """Resolve configuration, run requested Stage-10 stats, and exit."""
    args = build_arg_parser().parse_args()
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if args.line_lag_window_sec <= 0:
        raise SystemExit("--line-lag-window-sec must be > 0.")
    if args.n_permutations < 1:
        raise SystemExit("--n-permutations must be >= 1.")
    if not 0 < args.cluster_alpha < 1:
        raise SystemExit("--cluster-alpha must be between 0 and 1.")
    if not 0 < args.significance_alpha < 1:
        raise SystemExit("--significance-alpha must be between 0 and 1.")
    if "none" in set(args.paired_contrasts or []) and len(set(args.paired_contrasts or [])) > 1:
        raise SystemExit("--paired-contrasts none cannot be combined with other contrasts.")

    base_subject_config = args.subject_config
    if args.subjects:
        base_subject_config = _subject_config_for_id(args, normalize_subject_label(args.subjects[0]))
    base_cfg = _load_cfg(args, subject_config=base_subject_config)
    subjects = _subjects_from_args_or_config(args, base_cfg)

    try:
        report = run_group_scope(base_cfg, args, subjects)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"dRSA stats error: {exc}")
        return 2

    if args.print_summary:
        print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
