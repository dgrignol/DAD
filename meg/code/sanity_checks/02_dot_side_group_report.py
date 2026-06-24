#!/usr/bin/env python3
"""
02_dot_side_group_report.py

Purpose
Build a combined group-level report for the optional dot-side decoding sanity
check. The script consumes completed per-subject outputs from
``01_dot_side_decoding.py``, tests each requested neural-condition/model-source
analysis against chance with a 1D sign-flip cluster permutation test, writes
group figures and machine-readable JSON, and builds an HTML report containing
all group analyses before any participant-specific figures. Occluded-condition
accuracy plots mark the configured complete-occlusion and reappearance frames,
and each participant section starts with source-trial counts by identity for
all three neural conditions.

Usage examples
1) Build the default group report locally after subjects 01-07 have completed:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/02_dot_side_group_report.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --analyses all --group-subjects 01 02 03 04 05 06 07 \
     --overwrite

2) Validate inputs without writing outputs:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/02_dot_side_group_report.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --analyses all --group-subjects 01 02 03 04 05 06 07 \
     --validate-only --print-summary

3) Use SEM instead of the default 95% confidence interval display band:

   python .../02_dot_side_group_report.py --analyses all \
     --group-subjects 01 02 03 04 05 06 07 --interval sem --overwrite

4) Smooth each participant's time series before group statistics:

   python .../02_dot_side_group_report.py --analyses all \
     --group-subjects 01 02 03 04 05 06 07 \
     --smooth-sigma-samples 1.0 --overwrite

5) Run against cluster paths from the deployed DAD directory:

   python meg/code/sanity_checks/02_dot_side_group_report.py \
     --config meg/params/pipeline_default.yaml \
     --site-config meg/params/pipeline_cluster.yaml \
     --analyses all --group-subjects 01 02 03 04 05 06 07 \
     --overwrite

Inputs
- Per-subject summary JSON from ``01_dot_side_decoding.py``.
- Per-subject scores CSV referenced by each summary JSON.
- Per-subject three-panel PNG referenced by each summary JSON.

Outputs
- One group PNG per analysis under
  ``derivatives/sanity_checks/group/ses-01/dot_side_decoding/<condition>/<model-source>/``.
- Combined JSON and linked/embedded HTML reports under
  ``derivatives/reports/group/ses-01/dot_side_decoding/combined/``.
- One participant trial-count PNG per subject under the combined report's
  ``participant_trial_counts`` directory.

Statistical policy
- The default metric is ``accuracy_balanced_undersampled_sets`` to keep the
  current undersampling strategy as the group-level target.
- The default test is one-sample sign-flip cluster permutation on
  ``accuracy - 0.5`` with a directional greater-than-chance alternative.
- The default cluster-forming threshold is p=0.05; the threshold p and
  equivalent t cutoff are annotated directly on the group plot.
- The default display band is a 95% confidence interval. SEM is available with
  ``--interval sem``.
- Optional ``--smooth-sigma-samples`` applies Gaussian smoothing separately
  within each subject's finite contiguous time segments before group statistics
  and plotting. NaN windows remain NaN so clusters cannot bridge skipped
  windows.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
from itertools import product
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage, stats


PREPROCESSING_DIR = Path(__file__).resolve().parents[1] / "preprocessing_basic"
sys.path.insert(0, str(PREPROCESSING_DIR))

from pipeline_config import load_and_resolve_config  # noqa: E402
from drsa_result_utils import occlusion_timing_from_config  # noqa: E402
from preprocessing_stage_utils import guarded_path, utc_now, write_json  # noqa: E402


DEFAULT_SUBJECTS = ["01", "02", "03", "04", "05", "06", "07"]
DEFAULT_METRIC = "accuracy_balanced_undersampled_sets"
SUPPORTED_METRICS = (
    "accuracy_balanced_undersampled_sets",
    "accuracy_unbalanced_train",
)
STANDARD_ANALYSES = (
    {
        "key": "always_visible_observed",
        "condition": "always_visible",
        "model_source": "observed",
        "model_condition": "always_visible",
        "label": "Always visible | observed model",
    },
    {
        "key": "occluded_nondeviant_observed",
        "condition": "occluded_nondeviant",
        "model_source": "observed",
        "model_condition": "occluded_nondeviant",
        "label": "Occluded nondeviant | observed model",
    },
    {
        "key": "occluded_deviant_observed",
        "condition": "occluded_deviant",
        "model_source": "observed",
        "model_condition": "occluded_deviant",
        "label": "Occluded deviant | observed model",
    },
    {
        "key": "occluded_deviant_predicted",
        "condition": "occluded_deviant",
        "model_source": "predicted",
        "model_condition": "occluded_nondeviant",
        "label": "Occluded deviant | predicted model",
    },
)
ANALYSIS_BY_KEY = {str(item["key"]): item for item in STANDARD_ANALYSES}


# -----------------------------------------------------------------------------
# CLI and path resolution
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the group dot-side decoding report command-line parser."""
    meg_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Build group-level dot-side decoding sanity-check report."
    )
    parser.add_argument(
        "--config",
        default=str(meg_root / "params" / "pipeline_default.yaml"),
        help="Shared pipeline default YAML.",
    )
    parser.add_argument("--site-config", default=None)
    parser.add_argument(
        "--subject-config",
        default=None,
        help=(
            "Subject selector YAML used only to resolve shared paths. Defaults "
            "to sub-01."
        ),
    )
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument(
        "--analyses",
        nargs="+",
        default=None,
        help=(
            "Analysis keys to include, or 'all' for the four standard "
            "analyses. Defaults to all when neither --analyses nor "
            "--condition is supplied."
        ),
    )
    parser.add_argument(
        "--condition",
        default=None,
        help="Backward-compatible single-analysis neural condition selector.",
    )
    parser.add_argument(
        "--model-source",
        choices=("observed", "predicted"),
        default="observed",
        help="Model source used with the backward-compatible --condition mode.",
    )
    parser.add_argument(
        "--group-subjects",
        nargs="+",
        default=DEFAULT_SUBJECTS,
        help="Subjects to include, e.g. 01 02 sub-03.",
    )
    parser.add_argument(
        "--metric",
        choices=SUPPORTED_METRICS,
        default=DEFAULT_METRIC,
        help="Scores CSV column to test at group level.",
    )
    parser.add_argument(
        "--smooth-sigma-samples",
        type=float,
        default=0.0,
        help=(
            "Optional Gaussian smoothing sigma in decoded time samples, "
            "applied per subject before group statistics. 0 disables smoothing."
        ),
    )
    parser.add_argument(
        "--chance",
        type=float,
        default=0.5,
        help="Chance level to subtract before the one-sample test.",
    )
    parser.add_argument(
        "--cluster-forming-p",
        type=float,
        default=0.05,
        help="Pointwise p threshold used to form temporal clusters.",
    )
    parser.add_argument(
        "--cluster-alpha",
        type=float,
        default=0.05,
        help="Cluster-level p threshold used for significant cluster shading.",
    )
    parser.add_argument(
        "--tail",
        choices=("greater", "two-sided"),
        default="greater",
        help="Alternative hypothesis for the one-sample cluster test.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=0,
        help=(
            "Random sign-flip permutations. 0 uses exact enumeration when "
            "possible and otherwise 5000 random permutations."
        ),
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
        help="Seed for randomized sign-flip permutations.",
    )
    parser.add_argument(
        "--interval",
        choices=("ci95", "sem"),
        default="ci95",
        help="Display band around the group mean. Default is 95% CI.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Build from available subjects and list missing subjects instead of failing.",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


def normalize_subject_label(value: str) -> str:
    """Normalize subject input such as ``5`` or ``sub-05`` to ``sub-05``."""
    text = str(value).strip()
    if text.lower().startswith("sub-"):
        text = text[4:]
    return f"sub-{int(text):02d}"


def resolve_analysis_specs(args: argparse.Namespace) -> list[dict[str, str]]:
    """Resolve requested analysis keys into validated ordered specifications."""
    if args.analyses and args.condition:
        raise SystemExit("--analyses and --condition are mutually exclusive.")
    if args.analyses:
        requested = [str(value) for value in args.analyses]
    elif args.condition:
        requested = [f"{args.condition}_{args.model_source}"]
    else:
        requested = ["all"]

    if "all" in requested:
        if len(requested) != 1:
            raise SystemExit("The analysis selector 'all' cannot be combined with other keys.")
        return [dict(item) for item in STANDARD_ANALYSES]

    unknown = [key for key in requested if key not in ANALYSIS_BY_KEY]
    if unknown:
        raise SystemExit(
            f"Unsupported analysis keys: {unknown}. Expected one or more of "
            f"{list(ANALYSIS_BY_KEY)} or 'all'."
        )
    if len(set(requested)) != len(requested):
        raise SystemExit("Duplicate analysis keys are not allowed.")
    return [dict(ANALYSIS_BY_KEY[key]) for key in requested]


def resolve_config(args: argparse.Namespace) -> tuple[dict[str, Any], list[str]]:
    """Load configuration using sub-01 only as a shared path resolver."""
    meg_root = Path(__file__).resolve().parents[2]
    subject_config = (
        Path(args.subject_config).expanduser().resolve()
        if args.subject_config
        else meg_root / "params" / "sub-01.yaml"
    )
    cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("epoching",),
    )
    return cfg, warnings


def combined_output_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return output paths for the combined group report."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    report_dir = (
        derivatives_root
        / "reports"
        / "group"
        / session
        / "dot_side_decoding"
        / "combined"
    )
    stem = "dot_side_all_analyses_group"
    return {
        "report_dir": report_dir,
        "report_json": report_dir / f"{stem}_report.json",
        "report_html": report_dir / f"{stem}_report.html",
        "report_embedded_html": report_dir / f"{stem}_report_embedded.html",
    }


def group_figure_path(cfg: Mapping[str, Any], analysis: Mapping[str, str]) -> Path:
    """Return the group-figure path for one condition/model-source analysis."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    return (
        derivatives_root
        / "sanity_checks"
        / "group"
        / session
        / "dot_side_decoding"
        / analysis["condition"]
        / analysis["model_source"]
        / f"{analysis['key']}_dot_side_group_accuracy.png"
    )


def participant_trial_count_figure_path(
    cfg: Mapping[str, Any],
    subject: str,
) -> Path:
    """Return the combined report's trial-count figure path for one subject."""
    return (
        combined_output_paths(cfg)["report_dir"]
        / "participant_trial_counts"
        / f"{subject}_trial_counts_by_identity.png"
    )


def subject_output_paths(
    *,
    derivatives_root: Path,
    subject: str,
    session: str,
    condition: str,
    model_source: str = "observed",
) -> dict[str, Path]:
    """Return canonical subject-level dot-side output paths."""
    out_dir = (
        derivatives_root
        / "sanity_checks"
        / subject
        / session
        / "dot_side_decoding"
        / condition
        / model_source
    )
    stem = f"{condition}_{model_source}"
    return {
        "out_dir": out_dir,
        "summary_json": out_dir / f"{stem}_dot_side_decoding_summary.json",
        "scores_csv": out_dir / f"{stem}_dot_side_decoding_scores.csv",
        "figure": out_dir / f"{stem}_dot_side_decoding_3panel.png",
    }


def validate_cli_args(args: argparse.Namespace) -> None:
    """Validate basic CLI settings before reading any inputs."""
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if not 0 < args.cluster_forming_p < 1:
        raise SystemExit("--cluster-forming-p must be between 0 and 1.")
    if not 0 < args.cluster_alpha < 1:
        raise SystemExit("--cluster-alpha must be between 0 and 1.")
    if args.n_permutations < 0:
        raise SystemExit("--n-permutations must be non-negative.")
    if args.smooth_sigma_samples < 0:
        raise SystemExit("--smooth-sigma-samples must be non-negative.")


# -----------------------------------------------------------------------------
# Input loading and validation
# -----------------------------------------------------------------------------


def load_subject_score(
    *,
    derivatives_root: Path,
    subject: str,
    session: str,
    condition: str,
    model_source: str = "observed",
    metric: str,
) -> tuple[dict[str, Any] | None, dict[str, str] | None]:
    """
    Load one subject's dot-side summary and score series.

    Returns
    -------
    tuple
        ``(record, None)`` when loading succeeds, otherwise
        ``(None, skipped_input)`` describing the missing or invalid input.
    """
    paths = subject_output_paths(
        derivatives_root=derivatives_root,
        subject=subject,
        session=session,
        condition=condition,
        model_source=model_source,
    )
    summary_path = paths["summary_json"]
    if not summary_path.is_file():
        return None, {
            "subject": subject,
            "reason": f"missing summary JSON: {summary_path}",
        }
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("condition") != condition:
            raise RuntimeError(
                f"summary condition mismatch: {summary.get('condition')!r} != {condition!r}"
            )
        if summary.get("model_source") != model_source:
            raise RuntimeError(
                "summary model_source mismatch: "
                f"{summary.get('model_source')!r} != {model_source!r}"
            )
        scores_csv = Path(str(summary.get("scores_csv", paths["scores_csv"]))).expanduser().resolve()
        figure_path = Path(str(summary.get("figure_path", paths["figure"]))).expanduser().resolve()
        if not scores_csv.is_file():
            raise FileNotFoundError(f"missing scores CSV: {scores_csv}")
        if not figure_path.is_file():
            raise FileNotFoundError(f"missing participant figure: {figure_path}")
        scores = pd.read_csv(scores_csv)
        for column in ("time_sec", metric):
            if column not in scores.columns:
                raise RuntimeError(f"scores CSV lacks required column {column!r}")
        series = scores[metric].to_numpy(dtype=float)
        times = scores["time_sec"].to_numpy(dtype=float)
        if series.size == 0:
            raise RuntimeError("score series is empty")
        if not np.all(np.isfinite(times)):
            raise RuntimeError("time_sec contains non-finite values")
        record = {
            "subject": subject,
            "summary_path": str(summary_path),
            "scores_csv": str(scores_csv),
            "figure": str(figure_path),
            "summary": summary,
            "times": times,
            "values": series,
        }
        return record, None
    except Exception as exc:  # noqa: BLE001 - report the concrete skipped input.
        return None, {
            "subject": subject,
            "path": str(summary_path),
            "reason": str(exc),
        }


def trial_averaging_settings_from_summary(summary: Mapping[str, Any]) -> dict[str, Any]:
    """
    Extract the trial-averaging setting used for compatibility checks.

    Older summaries did not contain an explicit mode, so they are interpreted
    as fixed-count summaries with their recorded ``trials_per_average`` value.
    """
    trial_averaging = summary.get("trial_averaging", {})
    return {
        "mode": trial_averaging.get("mode", "fixed_count"),
        "trials_per_average": trial_averaging.get("trials_per_average"),
        "requested_value": trial_averaging.get(
            "requested_value",
            str(trial_averaging.get("trials_per_average")),
        ),
    }


def load_group_inputs(
    *,
    cfg: Mapping[str, Any],
    subjects: Iterable[str],
    condition: str,
    model_source: str = "observed",
    metric: str,
    allow_missing: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """
    Load and validate all requested subject-level dot-side outputs.

    Data flow
    - Derivative root + subject list -> canonical summary paths.
    - Summary JSON -> scores CSV and participant figure path.
    - Scores CSV -> time axis and selected metric series.
    - Cross-subject validation ensures all available subjects share one time
      axis and compatible decoding settings.
    """
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    records: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    for subject in subjects:
        record, skip = load_subject_score(
            derivatives_root=derivatives_root,
            subject=subject,
            session=session,
            condition=condition,
            model_source=model_source,
            metric=metric,
        )
        if record is not None:
            records.append(record)
        if skip is not None:
            skipped.append(skip)

    if skipped and not allow_missing:
        examples = "; ".join(item["reason"] for item in skipped[:3])
        raise RuntimeError(
            "Missing or invalid subject dot-side outputs. "
            f"Use --allow-missing to build a partial report. Examples: {examples}"
        )
    if len(records) < 2:
        raise RuntimeError(
            "At least two valid subject outputs are required for group statistics; "
            f"found {len(records)}."
        )

    reference = records[0]
    reference_times = reference["times"]
    reference_summary = reference["summary"]
    setting_keys = (
        "condition",
        "model_source",
        "model_condition",
        "analysis_key",
        "n_times",
        "sfreq_hz",
        "decode_window_samples",
        "decode_tmin_sec",
        "decode_tmax_sec",
        "n_splits",
        "random_seed",
    )
    for record in records[1:]:
        if record["times"].shape != reference_times.shape or not np.allclose(record["times"], reference_times):
            raise RuntimeError(
                "Subject time axis mismatch between "
                f"{reference['subject']} and {record['subject']}."
            )
        for key in setting_keys:
            if record["summary"].get(key) != reference_summary.get(key):
                raise RuntimeError(
                    f"Subject setting mismatch for {key}: "
                    f"{reference['subject']}={reference_summary.get(key)!r}, "
                    f"{record['subject']}={record['summary'].get(key)!r}"
                )
        ref_tpa = trial_averaging_settings_from_summary(reference_summary)
        cur_tpa = trial_averaging_settings_from_summary(record["summary"])
        if ref_tpa != cur_tpa:
            raise RuntimeError(
                "Subject trial-averaging setting mismatch: "
                f"{reference['subject']}={ref_tpa!r}, {record['subject']}={cur_tpa!r}"
            )

    return records, skipped


def validate_cross_analysis_settings(
    analysis_results: Iterable[Mapping[str, Any]],
) -> None:
    """Require one decoding configuration across all analyses in a report."""
    results = list(analysis_results)
    if len(results) < 2:
        return

    scalar_keys = (
        "sfreq_hz",
        "decode_window_samples",
        "decode_tmin_sec",
        "decode_tmax_sec",
        "n_splits",
        "random_seed",
    )

    def settings(result: Mapping[str, Any]) -> dict[str, Any]:
        records = list(result["participant_records"])
        if not records:
            raise RuntimeError(f"Analysis {result['key']} contains no participant records.")
        summary = records[0]["summary"]
        values = {key: summary.get(key) for key in scalar_keys}
        values["trial_averaging"] = trial_averaging_settings_from_summary(summary)
        return values

    reference = results[0]
    reference_settings = settings(reference)
    for result in results[1:]:
        current_settings = settings(result)
        if current_settings != reference_settings:
            raise RuntimeError(
                "Cross-analysis decoding setting mismatch between "
                f"{reference['key']}={reference_settings} and "
                f"{result['key']}={current_settings}."
            )


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------


def one_sample_t(values: np.ndarray) -> np.ndarray:
    """Compute one-sample t values along subjects for a subject x time array."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("values must be subject x time.")
    mean = np.nanmean(values, axis=0)
    sd = np.nanstd(values, axis=0, ddof=1)
    n = np.sum(np.isfinite(values), axis=0)
    sem = np.divide(sd, np.sqrt(n), out=np.full_like(sd, np.nan), where=n > 0)
    t_values = np.divide(mean, sem, out=np.zeros_like(mean), where=sem > 0)
    t_values[(sem == 0) & (mean > 0)] = np.inf
    t_values[(sem == 0) & (mean < 0)] = -np.inf
    return t_values


def cluster_threshold(
    *,
    n_subjects: int,
    cluster_forming_p: float,
    tail: str,
) -> float:
    """Return the t-value threshold for cluster formation."""
    df = n_subjects - 1
    if df < 1:
        raise ValueError("At least two subjects are required.")
    if tail == "greater":
        return float(stats.t.ppf(1.0 - cluster_forming_p, df=df))
    if tail == "two-sided":
        return float(stats.t.ppf(1.0 - (cluster_forming_p / 2.0), df=df))
    raise ValueError(f"Unsupported tail: {tail}")


def contiguous_true_regions(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return inclusive start/end index pairs for contiguous True regions."""
    mask = np.asarray(mask, dtype=bool)
    regions: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(mask):
        if value and start is None:
            start = idx
        elif not value and start is not None:
            regions.append((start, idx - 1))
            start = None
    if start is not None:
        regions.append((start, mask.size - 1))
    return regions


def cluster_masses(t_values: np.ndarray, threshold: float, tail: str) -> list[dict[str, Any]]:
    """Find temporal clusters and compute cluster masses from t values."""
    if tail == "greater":
        mask = t_values > threshold
        mass_source = t_values
    elif tail == "two-sided":
        mask = np.abs(t_values) > threshold
        mass_source = np.abs(t_values)
    else:
        raise ValueError(f"Unsupported tail: {tail}")
    clusters: list[dict[str, Any]] = []
    for start_idx, end_idx in contiguous_true_regions(mask):
        mass = float(np.nansum(mass_source[start_idx : end_idx + 1]))
        clusters.append(
            {
                "start_index": int(start_idx),
                "end_index": int(end_idx),
                "mass": mass,
                "n_timepoints": int(end_idx - start_idx + 1),
            }
        )
    return clusters


def finite_group_time_mask(values: np.ndarray) -> np.ndarray:
    """
    Return time points with finite scores for every included subject.

    Group sign-flip tests use one fixed subject set. When all-trials
    per-identity decoding skips a subject/time window, that time point is
    excluded for the group statistic rather than changing N across time.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("values must be subject x time.")
    return np.all(np.isfinite(values), axis=0)


def smooth_subject_series_1d(series: np.ndarray, *, sigma_samples: float) -> np.ndarray:
    """
    Smooth one subject's time series without crossing NaN gaps.

    Data flow
    - Identify finite contiguous regions in the subject series.
    - Apply Gaussian smoothing separately inside each region.
    - Restore the original NaN mask so skipped windows remain hard boundaries
      for group statistics and temporal clusters.
    """
    series = np.asarray(series, dtype=float)
    smoothed = np.full(series.shape, np.nan, dtype=float)
    if sigma_samples <= 0:
        return np.array(series, dtype=float, copy=True)
    for start_idx, end_idx in contiguous_true_regions(np.isfinite(series)):
        segment = series[start_idx : end_idx + 1]
        smoothed[start_idx : end_idx + 1] = ndimage.gaussian_filter1d(
            segment,
            sigma=float(sigma_samples),
            mode="nearest",
        )
    return smoothed


def smooth_group_values_1d(values: np.ndarray, *, sigma_samples: float) -> np.ndarray:
    """Apply NaN-preserving 1D Gaussian smoothing to every subject row."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("values must be subject x time.")
    if sigma_samples <= 0:
        return np.array(values, dtype=float, copy=True)
    return np.stack(
        [
            smooth_subject_series_1d(row, sigma_samples=float(sigma_samples))
            for row in values
        ],
        axis=0,
    )


def _sign_matrix(
    *,
    n_subjects: int,
    n_permutations: int,
    random_seed: int,
) -> tuple[np.ndarray, str]:
    """Return exact or randomized sign flips for one-sample permutation."""
    if n_permutations == 0 and n_subjects <= 15:
        signs = np.array(list(product([-1.0, 1.0], repeat=n_subjects)), dtype=float)
        return signs, "exact"
    rng = np.random.default_rng(random_seed)
    count = int(n_permutations or 5000)
    signs = rng.choice([-1.0, 1.0], size=(count, n_subjects))
    signs[0, :] = 1.0
    return signs.astype(float), "randomized"


def cluster_permutation_1d(
    values: np.ndarray,
    *,
    chance: float,
    cluster_forming_p: float,
    cluster_alpha: float,
    tail: str,
    n_permutations: int,
    random_seed: int,
) -> dict[str, Any]:
    """
    Test mean accuracy against chance using 1D sign-flip cluster permutation.

    ``values`` is subject x time accuracy. The tested quantity is
    ``values - chance``. Cluster p-values are computed against the permutation
    distribution of the maximum cluster mass. Time points with at least one
    non-finite subject value are excluded from cluster formation; invalid
    points remain NaN in the returned t series so clusters cannot bridge them.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("values must be subject x time.")
    if values.shape[0] < 2:
        raise ValueError("At least two subjects are required.")
    valid_time_mask = finite_group_time_mask(values)
    if not np.any(valid_time_mask):
        raise ValueError(
            "Group metric contains no time points with finite values for "
            "every included subject."
        )

    diffs_valid = values[:, valid_time_mask] - chance
    n_subjects = int(values.shape[0])
    threshold = cluster_threshold(
        n_subjects=n_subjects,
        cluster_forming_p=cluster_forming_p,
        tail=tail,
    )
    observed_t = np.full(values.shape[1], np.nan, dtype=float)
    observed_t[valid_time_mask] = one_sample_t(diffs_valid)
    observed_clusters = cluster_masses(observed_t, threshold, tail)

    signs, mode = _sign_matrix(
        n_subjects=n_subjects,
        n_permutations=n_permutations,
        random_seed=random_seed,
    )
    null_max_masses = np.zeros(signs.shape[0], dtype=float)
    for perm_idx, sign_row in enumerate(signs):
        permuted = diffs_valid * sign_row[:, np.newaxis]
        perm_t = np.full(values.shape[1], np.nan, dtype=float)
        perm_t[valid_time_mask] = one_sample_t(permuted)
        perm_clusters = cluster_masses(perm_t, threshold, tail)
        null_max_masses[perm_idx] = max(
            [cluster["mass"] for cluster in perm_clusters],
            default=0.0,
        )

    for cluster in observed_clusters:
        p_value = float(np.mean(null_max_masses >= float(cluster["mass"])))
        cluster["p_value"] = p_value
        cluster["significant"] = bool(p_value <= cluster_alpha)

    return {
        "chance": float(chance),
        "tail": tail,
        "cluster_forming_p": float(cluster_forming_p),
        "cluster_alpha": float(cluster_alpha),
        "t_threshold": float(threshold),
        "permutation_mode": mode,
        "n_permutations": int(signs.shape[0]),
        "observed_t": observed_t,
        "valid_time_mask": valid_time_mask,
        "n_valid_timepoints": int(valid_time_mask.sum()),
        "n_invalid_timepoints": int((~valid_time_mask).sum()),
        "clusters": observed_clusters,
        "null_max_cluster_masses": null_max_masses,
    }


def group_interval(values: np.ndarray, *, interval: str, ci_level: float = 0.95) -> dict[str, np.ndarray | str | float]:
    """Return mean and display interval for a subject x time score matrix."""
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    n_per_time = finite.sum(axis=0).astype(float)
    summed = np.nansum(values, axis=0)
    mean = np.divide(
        summed,
        n_per_time,
        out=np.full(values.shape[1], np.nan, dtype=float),
        where=n_per_time > 0,
    )
    centered = np.where(finite, values - mean[np.newaxis, :], 0.0)
    ss = np.sum(centered * centered, axis=0)
    variance = np.divide(
        ss,
        n_per_time - 1.0,
        out=np.full(values.shape[1], np.nan, dtype=float),
        where=n_per_time > 1,
    )
    sem = np.divide(
        np.sqrt(variance),
        np.sqrt(n_per_time),
        out=np.full(values.shape[1], np.nan, dtype=float),
        where=n_per_time > 1,
    )
    if interval == "sem":
        low = mean - sem
        high = mean + sem
        label = "mean +/- SEM"
    elif interval == "ci95":
        tcrit = stats.t.ppf(
            1.0 - ((1.0 - ci_level) / 2.0),
            df=np.maximum(n_per_time - 1.0, 1.0),
        )
        tcrit[n_per_time <= 1] = np.nan
        low = mean - (tcrit * sem)
        high = mean + (tcrit * sem)
        label = f"{int(ci_level * 100)}% CI"
    else:
        raise ValueError(f"Unsupported interval: {interval}")
    return {
        "mean": mean,
        "low": low,
        "high": high,
        "sem": sem,
        "n_per_time": n_per_time,
        "label": label,
        "ci_level": ci_level if interval == "ci95" else np.nan,
    }


# -----------------------------------------------------------------------------
# Plotting and report writing
# -----------------------------------------------------------------------------


def _time_for_cluster(times: np.ndarray, index: int) -> float:
    """Return a valid time value for a cluster index."""
    index = int(np.clip(index, 0, len(times) - 1))
    return float(times[index])


def add_occlusion_boundary_lines(
    ax: plt.Axes,
    *,
    condition: str,
    occlusion_timing: Mapping[str, Any],
) -> None:
    """Mark complete occlusion and reappearance on occluded accuracy plots."""
    if "occluded" not in str(condition).lower():
        return
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


def subject_trial_counts_by_condition(
    *,
    subject: str,
    analysis_results: Iterable[Mapping[str, Any]],
) -> tuple[list[int], dict[str, list[int]]]:
    """
    Collect one subject's source-trial counts for each neural condition.

    Counts come from each decoder summary's pre-averaging identity inventory.
    Missing identities are represented as zero. The predicted deviant result
    is not plotted twice, but its counts must match the observed deviant result
    because both analyses consume the same neural epochs.
    """
    identity_values: list[int] | None = None
    counts_by_condition: dict[str, list[int]] = {}

    for result in analysis_results:
        record = next(
            (
                item
                for item in result["participant_records"]
                if item["subject"] == subject
            ),
            None,
        )
        if record is None:
            continue
        summary = record["summary"]
        identity_range = summary.get("identity_range")
        if not isinstance(identity_range, list) or len(identity_range) != 2:
            raise RuntimeError(
                f"{subject} {result['key']} summary lacks a valid identity_range."
            )
        current_identities = list(
            range(int(identity_range[0]), int(identity_range[1]) + 1)
        )
        if identity_values is None:
            identity_values = current_identities
        elif current_identities != identity_values:
            raise RuntimeError(
                f"{subject} analyses disagree on the configured identity range."
            )

        source_counts = summary.get("trial_averaging", {}).get(
            "n_source_trials_by_identity"
        )
        if not isinstance(source_counts, Mapping):
            raise RuntimeError(
                f"{subject} {result['key']} summary lacks source trial counts by identity."
            )
        condition = str(result["condition"])
        current_counts = [
            int(source_counts.get(str(identity), source_counts.get(identity, 0)))
            for identity in current_identities
        ]
        if condition in counts_by_condition:
            if current_counts != counts_by_condition[condition]:
                raise RuntimeError(
                    f"{subject} duplicate analyses disagree on {condition} trial counts."
                )
        else:
            counts_by_condition[condition] = current_counts

    if identity_values is None:
        raise RuntimeError(f"No participant records are available for {subject}.")
    expected_conditions: list[str] = []
    for result in analysis_results:
        condition = str(result["condition"])
        if condition in counts_by_condition and condition not in expected_conditions:
            expected_conditions.append(condition)
    return identity_values, {
        condition: counts_by_condition[condition] for condition in expected_conditions
    }


def plot_subject_trial_counts(
    *,
    subject: str,
    identity_values: list[int],
    counts_by_condition: Mapping[str, list[int]],
    out_file: Path,
    raw_root: Path,
) -> Path:
    """Write one three-panel source-trial count figure for a participant."""
    out_path = guarded_path(out_file, raw_root)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "always_visible": "#4c78a8",
        "occluded_nondeviant": "#59a14f",
        "occluded_deviant": "#e15759",
    }
    max_count = max(
        max((int(value) for value in counts), default=0)
        for counts in counts_by_condition.values()
    )
    fig, axes = plt.subplots(
        len(counts_by_condition),
        1,
        figsize=(14, 8),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes_array = np.atleast_1d(axes)
    for ax, (condition, counts) in zip(
        axes_array,
        counts_by_condition.items(),
        strict=True,
    ):
        ax.bar(
            identity_values,
            counts,
            width=0.82,
            color=colors.get(condition, "#767676"),
            edgecolor="white",
            linewidth=0.35,
        )
        ax.set_title(f"{condition} | total source trials = {sum(counts)}")
        ax.set_ylabel("Trials")
        ax.set_ylim(0, max_count + 1)
        ax.grid(axis="y", color="0.88", lw=0.7)
    axes_array[-1].set_xlabel("Trajectory identity")
    axes_array[-1].set_xticks(identity_values)
    axes_array[-1].tick_params(axis="x", labelrotation=90, labelsize=8)
    fig.suptitle(f"{subject} | Source trial counts per identity before averaging")
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_group_accuracy(
    *,
    times: np.ndarray,
    values: np.ndarray,
    subjects: list[str],
    stats_result: Mapping[str, Any],
    interval_result: Mapping[str, Any],
    metric: str,
    condition: str,
    model_source: str,
    occlusion_timing: Mapping[str, Any],
    out_file: Path,
    raw_root: Path,
) -> Path:
    """Write the group accuracy figure with cluster-threshold annotations."""
    out_path = guarded_path(out_file, raw_root)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mean = np.asarray(interval_result["mean"], dtype=float)
    low = np.asarray(interval_result["low"], dtype=float)
    high = np.asarray(interval_result["high"], dtype=float)

    fig, ax = plt.subplots(figsize=(13, 6), constrained_layout=True)
    for _subject, series in zip(subjects, values, strict=True):
        ax.plot(times, series, color="0.65", lw=0.8, alpha=0.65, label="_nolegend_")

    ax.fill_between(
        times,
        low,
        high,
        color="#4c78a8",
        alpha=0.22,
        label=str(interval_result["label"]),
    )
    ax.plot(times, mean, color="#1f77b4", lw=2.2, label=f"Group mean ({metric})")
    chance = float(stats_result["chance"])
    ax.axhline(chance, color="k", lw=1.0, ls=":", label=f"Chance = {chance:g}")
    ax.axvline(0.0, color="k", lw=1.0, ls="--", alpha=0.7)
    add_occlusion_boundary_lines(
        ax,
        condition=condition,
        occlusion_timing=occlusion_timing,
    )

    y_min, y_max = ax.get_ylim()
    y_span = y_max - y_min
    if y_span <= 0:
        y_span = 1.0
    cluster_y = y_max - (0.08 * y_span)
    for cluster in stats_result["clusters"]:
        if not cluster["significant"]:
            continue
        start_time = _time_for_cluster(times, cluster["start_index"])
        end_time = _time_for_cluster(times, cluster["end_index"])
        color = "#d62728"
        alpha = 0.18
        ax.axvspan(start_time, end_time, color=color, alpha=alpha, lw=0)
        ax.plot([start_time, end_time], [cluster_y, cluster_y], color=color, lw=4)
        ax.text(
            (start_time + end_time) / 2.0,
            cluster_y + (0.02 * y_span),
            f"cluster p={cluster['p_value']:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
            color=color,
        )

    threshold_text = (
        f"Cluster-forming threshold: p={stats_result['cluster_forming_p']:.3g}, "
        f"t>{stats_result['t_threshold']:.3f}; "
        f"cluster alpha={stats_result['cluster_alpha']:.3g}; "
        f"{stats_result['permutation_mode']} sign flips n={stats_result['n_permutations']}"
    )
    ax.text(
        0.01,
        0.98,
        threshold_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "0.75"},
    )

    ax.set_title(
        f"Group dot-side decoding | {condition} | {model_source} model | "
        f"n={len(subjects)}"
    )
    ax.set_xlabel("Time from condition onset (s)")
    ax.set_ylabel("Accuracy")
    ax.set_ylim(min(0.35, float(np.nanmin(low)) - 0.03), max(1.0, float(np.nanmax(high)) + 0.03))
    ax.legend(loc="lower right", frameon=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _image_src(report_path: Path, image_path: str | Path, *, embed: bool) -> str:
    """Return a relative or data-URI image source for HTML reports."""
    resolved = Path(image_path).expanduser().resolve()
    if embed:
        mime = "image/png" if resolved.suffix.lower() == ".png" else "application/octet-stream"
        payload = base64.b64encode(resolved.read_bytes()).decode("ascii")
        return f"data:{mime};base64,{payload}"
    try:
        return str(resolved.relative_to(report_path.parent))
    except ValueError:
        return str(resolved)


def _figure_html(
    report_path: Path,
    *,
    image_path: str | Path,
    caption: str,
    embed: bool,
) -> str:
    """Build one HTML figure block."""
    src = _image_src(report_path, image_path, embed=embed)
    return (
        "<figure>"
        f"<img src=\"{html.escape(src)}\"/>"
        f"<figcaption>{html.escape(caption)}</figcaption>"
        "</figure>"
    )


def _cluster_table_html(clusters: Iterable[Mapping[str, Any]]) -> str:
    """Build the temporal-cluster table shared by every group analysis."""
    rows = "\n".join(
        "<tr>"
        f"<td>{cluster['start_time_sec']:.3f}</td>"
        f"<td>{cluster['end_time_sec']:.3f}</td>"
        f"<td>{cluster['n_timepoints']}</td>"
        f"<td>{cluster['mass']:.3f}</td>"
        f"<td>{cluster['p_value']:.4f}</td>"
        f"<td>{'yes' if cluster['significant'] else 'no'}</td>"
        "</tr>"
        for cluster in clusters
    ) or "<tr><td colspan=\"6\">No clusters formed.</td></tr>"
    return (
        "<table><thead><tr><th>start s</th><th>end s</th><th>n timepoints</th>"
        "<th>mass</th><th>p</th><th>significant</th></tr></thead><tbody>"
        f"{rows}</tbody></table>"
    )


def _skipped_inputs_html(skipped_inputs: Iterable[Mapping[str, str]]) -> str:
    """Build a compact disclosure for missing or invalid analysis inputs."""
    items = list(skipped_inputs)
    if not items:
        return ""
    rows = "\n".join(
        "<li>"
        + html.escape(
            " | ".join(f"{key}: {value}" for key, value in sorted(item.items()))
        )
        + "</li>"
        for item in items
    )
    return f"<details><summary>Skipped or missing inputs</summary><ul>{rows}</ul></details>"


def _format_optional_number(value: Any, *, digits: int = 3) -> str:
    """Format optional numeric values for compact HTML tables."""
    if value is None:
        return "n/a"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def _diagnostics_from_summary(summary: Mapping[str, Any]) -> Mapping[str, Any] | None:
    """Return the nested training/testing diagnostics when present."""
    diagnostics = summary.get("training_testing_diagnostics")
    if isinstance(diagnostics, Mapping):
        return diagnostics
    return None


def _subject_diagnostics_html(summary: Mapping[str, Any]) -> str:
    """
    Render one participant's training/testing diagnostics after its result.

    Missing diagnostics are tolerated so older subject outputs can still be
    listed, but the report asks for a rebuild before interpreting them.
    """
    diagnostics = _diagnostics_from_summary(summary)
    if diagnostics is None:
        return (
            "<h4>Training/testing diagnostics</h4>"
            "<p>Diagnostics are unavailable for this participant output. "
            "Rebuild subject decoding with the current script before using "
            "this report to assess classifier validity.</p>"
        )

    classifier = diagnostics.get("classifier", {})
    cv = diagnostics.get("cross_validation", {})
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
    return (
        "<h4>Training/testing diagnostics</h4><p>"
        + html.escape(str(assessment.get("text", "")))
        + "</p><table><tbody>"
        + table_rows
        + "</tbody></table>"
    )


def _analysis_diagnostics_html(result: Mapping[str, Any]) -> str:
    """
    Summarize training/testing diagnostics across participants for one analysis.

    The table intentionally keeps one row per participant so class-balance
    problems caused by identity loss or all-trials averaging remain visible.
    """
    records = list(result["participant_records"])
    if not records:
        return "<p>No valid participant records were available.</p>"

    first_diagnostics = _diagnostics_from_summary(records[0]["summary"])
    if first_diagnostics is None:
        return (
            "<p>Training/testing diagnostics are unavailable for these outputs. "
            "Rebuild subject decoding with the current script.</p>"
        )

    classifier = first_diagnostics.get("classifier", {})
    cv = first_diagnostics.get("cross_validation", {})
    evaluation = first_diagnostics.get("evaluation", {})
    smoothing = result.get("smoothing", {})
    sigma_samples = float(smoothing.get("sigma_samples", 0.0) or 0.0)
    common_rows = [
        ("Classifier", classifier.get("name", "n/a")),
        ("Classifier details", classifier.get("details", "n/a")),
        ("CV scheme", cv.get("scheme", "n/a")),
        ("CV folds", cv.get("n_splits", "n/a")),
        ("Split unit", cv.get("split_unit", "n/a")),
        ("Balanced metric", evaluation.get("balanced_metric", "n/a")),
        ("Balancing strategy", evaluation.get("balancing_strategy", "n/a")),
        (
            "Smoothing",
            (
                "none"
                if sigma_samples <= 0
                else f"Gaussian sigma={sigma_samples:g} samples before group statistics"
            ),
        ),
        (
            "Group-valid time points",
            f"{result.get('n_group_valid_timepoints', 'n/a')} valid, "
            f"{result.get('n_group_invalid_timepoints', 'n/a')} skipped",
        ),
    ]
    common_table = "\n".join(
        "<tr><th>"
        + html.escape(str(label))
        + "</th><td>"
        + html.escape(str(value))
        + "</td></tr>"
        for label, value in common_rows
    )

    participant_rows: list[str] = []
    for record in records:
        diagnostics = _diagnostics_from_summary(record["summary"])
        if diagnostics is None:
            participant_rows.append(
                f"<tr><td>{html.escape(record['subject'])}</td>"
                "<td colspan=\"8\">diagnostics unavailable</td></tr>"
            )
            continue
        sample_structure = diagnostics.get("sample_structure", {})
        class_balance = diagnostics.get("class_balance", {})
        assessment = diagnostics.get("assessment", {})
        decoded_fraction = class_balance.get("decoded_window_fraction")
        decoded_percent = (
            "n/a"
            if decoded_fraction is None
            else f"{float(decoded_fraction) * 100.0:.1f}%"
        )
        participant_rows.append(
            "<tr>"
            f"<td>{html.escape(record['subject'])}</td>"
            f"<td>{html.escape(str(sample_structure.get('averaging_mode', 'n/a')))}</td>"
            f"<td>{html.escape(str(sample_structure.get('requested_value', 'n/a')))}</td>"
            f"<td>{html.escape(str(sample_structure.get('n_decoder_trials', 'n/a')))}</td>"
            f"<td>{html.escape(str(sample_structure.get('n_decoder_identities', 'n/a')))}</td>"
            f"<td>{html.escape(str(class_balance.get('n_decoded_windows', 'n/a')))} / "
            f"{html.escape(str(class_balance.get('n_windows', 'n/a')))} "
            f"({html.escape(decoded_percent)})</td>"
            f"<td>{html.escape(str(class_balance.get('min_minority_class_count_all_windows', 'n/a')))}</td>"
            f"<td>{html.escape(_format_optional_number(class_balance.get('max_majority_to_minority_ratio'), digits=2))}</td>"
            f"<td>{html.escape(str(assessment.get('status', 'n/a')))}</td>"
            "</tr>"
        )

    return (
        "<table><tbody>"
        + common_table
        + "</tbody></table>"
        "<table><thead><tr><th>subject</th><th>averaging mode</th>"
        "<th>request</th><th>pseudo-trials</th><th>identities</th>"
        "<th>decoded windows</th><th>min minority</th><th>worst ratio</th>"
        "<th>assessment</th></tr></thead><tbody>"
        + "\n".join(participant_rows)
        + "</tbody></table>"
    )


def write_combined_html_report(
    path: Path,
    *,
    title: str,
    summary_lines: Iterable[str],
    subjects: list[str],
    analysis_results: list[Mapping[str, Any]],
    participant_trial_count_figures: Mapping[str, str],
    raw_root: Path,
    embed_images: bool,
) -> Path:
    """
    Write one report with every group result before participant-specific data.

    Data flow
    - Analysis results -> group figures, cluster tables, and missing-input notes.
    - The same results -> per-subject lookup across all requested analyses and
      a source-trial count figure covering every neural condition.
    - Group section is emitted completely before the participant section.
    """
    report_path = guarded_path(path, raw_root)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    analysis_links = "\n".join(
        f"<li><a href=\"#group-{html.escape(result['key'])}\">"
        f"{html.escape(result['label'])}</a></li>"
        for result in analysis_results
    )
    subject_links = "\n".join(
        f"<li><a href=\"#{html.escape(subject)}\">{html.escape(subject)}</a></li>"
        for subject in subjects
    )

    # Emit all group-level analyses first, including their inferential tables.
    group_sections: list[str] = []
    for result in analysis_results:
        group_sections.append(
            "".join(
                [
                    f"<article id=\"group-{html.escape(result['key'])}\">",
                    f"<h3>{html.escape(result['label'])}</h3>",
                    _figure_html(
                        report_path,
                        image_path=result["group_figure"],
                        caption=(
                            "Group mean accuracy with individual traces, display interval, "
                            "and cluster-permutation results."
                        ),
                        embed=embed_images,
                    ),
                    "<h4>Temporal clusters</h4>",
                    _cluster_table_html(result["clusters"]),
                    "<h4>Training/testing diagnostics</h4>",
                    _analysis_diagnostics_html(result),
                    _skipped_inputs_html(result["skipped_inputs"]),
                    "</article>",
                ]
            )
        )

    # Reindex each analysis by subject, then place all analyses under a subject.
    records_by_analysis = {
        result["key"]: {
            record["subject"]: record for record in result["participant_records"]
        }
        for result in analysis_results
    }
    participant_sections: list[str] = []
    for subject in subjects:
        analysis_blocks: list[str] = []
        for result in analysis_results:
            record = records_by_analysis[result["key"]].get(subject)
            if record is None:
                content = "<p>No valid participant output was available.</p>"
            else:
                summary = record["summary"]
                content = _figure_html(
                    report_path,
                    image_path=record["figure"],
                    caption=(
                        f"{subject} | {result['label']} | "
                        f"peak balanced={summary.get('peak_accuracy_balanced_train')} "
                        f"at {summary.get('peak_time_sec_balanced_train')} s"
                    ),
                    embed=embed_images,
                )
                content += _subject_diagnostics_html(summary)
            analysis_blocks.append(
                f"<article><h3>{html.escape(result['label'])}</h3>{content}</article>"
            )
        trial_count_path = participant_trial_count_figures.get(subject)
        trial_count_content = (
            _figure_html(
                report_path,
                image_path=trial_count_path,
                caption=(
                    f"{subject} | Source trial counts per trajectory identity for "
                    "always-visible, occluded nondeviant, and occluded deviant epochs."
                ),
                embed=embed_images,
            )
            if trial_count_path is not None
            else "<p>No participant trial-count figure was available.</p>"
        )
        participant_sections.append(
            f"<section id=\"{html.escape(subject)}\"><h2>{html.escape(subject)}</h2>"
            f"<article><h3>Trial counts per identity</h3>{trial_count_content}</article>"
            f"{''.join(analysis_blocks)}</section>"
        )

    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            "<style>",
            "body{font-family:system-ui,sans-serif;max-width:1280px;margin:2rem auto;padding:0 1rem;}",
            "nav,pre{background:#f5f5f5;padding:1rem;border-radius:6px;}",
            "nav{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:1rem;}",
            "figure{margin:1.5rem 0;padding-bottom:1rem;border-bottom:1px solid #ddd;}",
            "img{max-width:100%;height:auto;border:1px solid #ddd;}",
            "figcaption{font-size:.9rem;color:#333;margin-top:.4rem;}",
            "table{border-collapse:collapse;margin:1rem 0;width:100%;}",
            "td,th{border:1px solid #ddd;padding:.4rem;text-align:left;}",
            "h2{margin-top:3rem;border-top:3px solid #222;padding-top:1rem;}",
            "h3{margin-top:2rem;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<pre>{html.escape(chr(10).join(summary_lines))}</pre>",
            "<nav><div><strong>Group analyses</strong><ul>",
            analysis_links,
            "</ul></div><div><strong>Subjects</strong><ul>",
            subject_links,
            "</ul></div></nav>",
            "<section id=\"group-analyses\"><h2>Group analyses</h2>",
            "\n".join(group_sections),
            "</section>",
            "<section id=\"participants\"><h2>Participant-specific analyses</h2>",
            "\n".join(participant_sections),
            "</section>",
            "</body></html>",
        ]
    )
    report_path.write_text(body, encoding="utf-8")
    return report_path


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def build_group_report(cfg: Mapping[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    """Load all requested analyses and write one group-first combined report."""
    raw_root = Path(cfg["resolved_paths"]["data_root"])
    subjects = [normalize_subject_label(value) for value in args.group_subjects]
    analysis_specs = resolve_analysis_specs(args)
    output_paths = combined_output_paths(cfg)
    participant_count_paths = {
        subject: participant_trial_count_figure_path(cfg, subject)
        for subject in subjects
    }
    expected_outputs = [
        output_paths["report_json"],
        output_paths["report_html"],
        output_paths["report_embedded_html"],
        *[group_figure_path(cfg, spec) for spec in analysis_specs],
    ]

    if args.resume and all(path.is_file() for path in expected_outputs):
        existing_report = json.loads(
            output_paths["report_json"].read_text(encoding="utf-8")
        )
        requested_keys = [spec["key"] for spec in analysis_specs]
        if existing_report.get("analysis_keys") != requested_keys:
            raise RuntimeError(
                "Existing combined report analysis keys do not match the "
                f"requested selection: existing={existing_report.get('analysis_keys')}, "
                f"requested={requested_keys}. Use --overwrite to rebuild."
            )
        return {
            "stage": "dot_side_group_report",
            "status": "completed_resume",
            "analysis_keys": requested_keys,
            "report_files": {
                key: str(value)
                for key, value in output_paths.items()
                if key != "report_dir"
            },
        }
    if not args.validate_only and not args.overwrite:
        existing = [str(path) for path in expected_outputs if path.exists()]
        if existing:
            raise RuntimeError(
                "Refusing to overwrite existing combined report outputs without "
                f"--overwrite or --resume: {existing}"
            )

    # -------------------------------------------------------------------------
    # Compute each analysis independently with the same group-statistic policy.
    # -------------------------------------------------------------------------
    analysis_results: list[dict[str, Any]] = []
    occlusion_timing = occlusion_timing_from_config(cfg)
    for spec in analysis_specs:
        records, skipped = load_group_inputs(
            cfg=cfg,
            subjects=subjects,
            condition=spec["condition"],
            model_source=spec["model_source"],
            metric=args.metric,
            allow_missing=args.allow_missing,
        )
        times = records[0]["times"]
        raw_values = np.stack([record["values"] for record in records], axis=0)
        values = smooth_group_values_1d(
            raw_values,
            sigma_samples=float(args.smooth_sigma_samples),
        )
        valid_subjects = [record["subject"] for record in records]

        stats_result = cluster_permutation_1d(
            values,
            chance=args.chance,
            cluster_forming_p=args.cluster_forming_p,
            cluster_alpha=args.cluster_alpha,
            tail=args.tail,
            n_permutations=args.n_permutations,
            random_seed=args.random_seed,
        )
        valid_time_mask = np.asarray(stats_result["valid_time_mask"], dtype=bool)
        values_for_group_display = np.array(values, dtype=float, copy=True)
        values_for_group_display[:, ~valid_time_mask] = np.nan
        interval_result = group_interval(values_for_group_display, interval=args.interval)

        clusters_for_report: list[dict[str, Any]] = []
        for cluster in stats_result["clusters"]:
            item = dict(cluster)
            item["start_index"] = int(cluster["start_index"])
            item["end_index"] = int(cluster["end_index"])
            item["start_time_sec"] = float(times[item["start_index"]])
            item["end_time_sec"] = float(times[item["end_index"]])
            clusters_for_report.append(item)

        figure_path = group_figure_path(cfg, spec)
        if not args.validate_only:
            figure_path = plot_group_accuracy(
                times=times,
                values=values_for_group_display,
                subjects=valid_subjects,
                stats_result=stats_result,
                interval_result=interval_result,
                metric=args.metric,
                condition=spec["condition"],
                model_source=spec["model_source"],
                occlusion_timing=occlusion_timing,
                out_file=figure_path,
                raw_root=raw_root,
            )

        analysis_results.append(
            {
                **spec,
                "subjects_requested": subjects,
                "subjects_included": valid_subjects,
                "skipped_inputs": skipped,
                "times_sec": times,
                "group_mean": interval_result["mean"],
                "group_interval_low": interval_result["low"],
                "group_interval_high": interval_result["high"],
                "group_valid_time_mask": valid_time_mask,
                "n_group_valid_timepoints": int(stats_result["n_valid_timepoints"]),
                "n_group_invalid_timepoints": int(stats_result["n_invalid_timepoints"]),
                "smoothing": {
                    "sigma_samples": float(args.smooth_sigma_samples),
                    "applied_to": args.metric,
                    "applied_before": "group_statistics_and_plotting",
                    "nan_policy": "smooth finite contiguous subject segments only",
                },
                "subject_values_raw": {
                    record["subject"]: record["values"] for record in records
                },
                "subject_values": {
                    record["subject"]: values[idx]
                    for idx, record in enumerate(records)
                },
                "clusters": clusters_for_report,
                "cluster_forming_t_threshold": float(stats_result["t_threshold"]),
                "permutation_mode": stats_result["permutation_mode"],
                "n_permutations": int(stats_result["n_permutations"]),
                "display_interval_label": str(interval_result["label"]),
                "group_figure": str(figure_path),
                "participant_figures": [
                    {
                        "subject": record["subject"],
                        "figure": record["figure"],
                        "scores_csv": record["scores_csv"],
                        "summary_json": record["summary_path"],
                    }
                    for record in records
                ],
                "participant_records": records,
            }
        )

    validate_cross_analysis_settings(analysis_results)

    if args.validate_only:
        return {
            "stage": "dot_side_group_report",
            "status": "validated",
            "analysis_keys": [result["key"] for result in analysis_results],
            "metric": args.metric,
            "analyses": {
                result["key"]: {
                    "subjects_requested": result["subjects_requested"],
                    "subjects_included": result["subjects_included"],
                    "skipped_inputs": result["skipped_inputs"],
                    "n_times": int(np.asarray(result["times_sec"]).size),
                    "cluster_count": len(result["clusters"]),
                    "significant_cluster_count": int(
                        sum(1 for cluster in result["clusters"] if cluster["significant"])
                    ),
                }
                for result in analysis_results
            },
        }

    # Build one participant QC figure from pre-averaging source counts. This is
    # intentionally done after all analyses load so duplicate deviant results
    # can be checked for count consistency before plotting.
    participant_trial_count_figures: dict[str, str] = {}
    for subject in subjects:
        if not any(
            any(
                record["subject"] == subject
                for record in result["participant_records"]
            )
            for result in analysis_results
        ):
            continue
        identity_values, counts_by_condition = subject_trial_counts_by_condition(
            subject=subject,
            analysis_results=analysis_results,
        )
        count_path = plot_subject_trial_counts(
            subject=subject,
            identity_values=identity_values,
            counts_by_condition=counts_by_condition,
            out_file=participant_count_paths[subject],
            raw_root=raw_root,
        )
        participant_trial_count_figures[subject] = str(count_path)

    # -------------------------------------------------------------------------
    # Serialize one machine report and two HTML variants in group-first order.
    # -------------------------------------------------------------------------
    summary_lines = [
        "Stage status: completed",
        f"Analyses: {', '.join(result['key'] for result in analysis_results)}",
        f"Metric: {args.metric}",
        f"Subjects requested: {', '.join(subjects)}",
        f"Chance: {args.chance:g}",
        f"Tail: {args.tail}",
        f"Cluster-forming p: {args.cluster_forming_p:g}",
        f"Cluster alpha: {args.cluster_alpha:g}",
        f"Smoothing sigma samples: {args.smooth_sigma_samples:g}",
        f"Display interval: {args.interval}",
        f"Machine report: {output_paths['report_json']}",
    ]
    report = {
        "stage": "dot_side_group_report",
        "status": "completed",
        "generated_utc": utc_now(),
        "analysis_keys": [result["key"] for result in analysis_results],
        "metric": args.metric,
        "session_bids_id": str(cfg["layout"]["session_bids_id"]),
        "subjects_requested": subjects,
        "chance": float(args.chance),
        "tail": args.tail,
        "cluster_forming_p": float(args.cluster_forming_p),
        "cluster_alpha": float(args.cluster_alpha),
        "interval": str(args.interval),
        "smooth_sigma_samples": float(args.smooth_sigma_samples),
        "occlusion_timing": occlusion_timing,
        "participant_trial_count_figures": participant_trial_count_figures,
        "analyses": {
            result["key"]: {
                key: value
                for key, value in result.items()
                if key != "participant_records"
            }
            for result in analysis_results
        },
        "report_files": {
            "json": str(output_paths["report_json"]),
            "html": str(output_paths["report_html"]),
            "embedded_html": str(output_paths["report_embedded_html"]),
        },
    }
    write_json(output_paths["report_json"], report, raw_root=raw_root, overwrite=True)
    title = f"Group | {cfg['layout']['session_bids_id']} | Dot-side decoding | all analyses"
    write_combined_html_report(
        output_paths["report_html"],
        title=title,
        summary_lines=summary_lines,
        subjects=subjects,
        analysis_results=analysis_results,
        participant_trial_count_figures=participant_trial_count_figures,
        raw_root=raw_root,
        embed_images=False,
    )
    write_combined_html_report(
        output_paths["report_embedded_html"],
        title=title,
        summary_lines=summary_lines,
        subjects=subjects,
        analysis_results=analysis_results,
        participant_trial_count_figures=participant_trial_count_figures,
        raw_root=raw_root,
        embed_images=True,
    )
    return report


def main() -> int:
    """Resolve configuration and build or validate the group report."""
    args = build_arg_parser().parse_args()
    validate_cli_args(args)
    cfg, warnings = resolve_config(args)
    if warnings:
        for warning in warnings:
            print(f"WARNING: {warning}", file=sys.stderr)
    try:
        report = build_group_report(cfg, args)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if args.print_summary:
        printable = {key: value for key, value in report.items() if key != "analyses"}
        if "analyses" in report:
            printable["analyses"] = {
                key: {
                    "subjects_included": value.get("subjects_included", []),
                    "skipped_inputs": value.get("skipped_inputs", []),
                    "n_times": len(value.get("times_sec", [])),
                    "significant_cluster_count": sum(
                        1
                        for cluster in value.get("clusters", [])
                        if cluster.get("significant")
                    ),
                    "group_figure": value.get("group_figure"),
                }
                for key, value in report["analyses"].items()
            }
        print(json.dumps(printable, indent=2))
    return 0 if report.get("status") in {"validated", "completed", "completed_resume"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
