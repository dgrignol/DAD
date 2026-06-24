#!/usr/bin/env python3
"""
04_dot_side_time_generalization_group_report.py

Purpose
Build a group-level report for dot-side time-generalization sanity checks.
The script consumes participant outputs from
``03_dot_side_time_generalization.py``, computes per-analysis group mean
train-time x test-time matrices, tests accuracy against chance with a 2D
sign-flip cluster permutation, and writes a group-first HTML report followed
by participant-specific matrix sections.

Usage examples
1) Build the default report for all four standard analyses:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/04_dot_side_time_generalization_group_report.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --analyses all --group-subjects 01 02 03 04 05 06 07 --overwrite

2) Validate inputs without writing outputs:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/04_dot_side_time_generalization_group_report.py \
     --analyses all --group-subjects 01 02 03 04 05 06 07 \
     --validate-only --print-summary

3) Run against cluster paths from the deployed DAD directory:

   python meg/code/sanity_checks/04_dot_side_time_generalization_group_report.py \
     --config meg/params/pipeline_default.yaml \
     --site-config meg/params/pipeline_cluster.yaml \
     --analyses all --group-subjects 01 02 03 04 05 06 07 --overwrite

Inputs
- Per-subject summary JSON from ``03_dot_side_time_generalization.py``.
- Per-subject NPZ matrix bundle referenced by each summary JSON.
- Per-subject heatmap PNG referenced by each summary JSON.

Outputs
- One group heatmap PNG and one group NPZ per analysis under
  ``derivatives/sanity_checks/group/ses-01/dot_side_time_generalization``.
- Combined JSON and linked/embedded HTML reports under
  ``derivatives/reports/group/ses-01/dot_side_time_generalization/combined``.

Statistical policy
- The default metric is ``accuracy_balanced_undersampled_sets``.
- The tested quantity is ``accuracy - 0.5``.
- A cell is included in group statistics only when every included subject has
  a finite value. Invalid cells remain NaN, so 2D clusters cannot bridge
  train/test gaps caused by skipped participant cells.
- Clusters are connected components on the 2D train/test lattice with
  4-neighbor adjacency. The default test is directional greater-than-chance,
  cluster-forming p=0.05, cluster alpha=0.05, and exact sign flips for
  subjects 01-07.
"""

from __future__ import annotations

import argparse
import base64
from collections import deque
from itertools import product
import html
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
import numpy as np
from scipy import stats


def _load_sibling_module(module_name: str, path: Path) -> Any:
    """Import a sibling script whose filename is not a valid Python module id."""
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {module_name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


LINE_GROUP = _load_sibling_module(
    "dot_side_group_report_helpers",
    Path(__file__).with_name("02_dot_side_group_report.py"),
)
TIME_GEN = _load_sibling_module(
    "dot_side_time_generalization_helpers",
    Path(__file__).with_name("03_dot_side_time_generalization.py"),
)


DEFAULT_METRIC = TIME_GEN.BALANCED_METRIC
SUPPORTED_METRICS = (TIME_GEN.BALANCED_METRIC, TIME_GEN.UNBALANCED_METRIC)
STANDARD_ANALYSES = LINE_GROUP.STANDARD_ANALYSES
ANALYSIS_BY_KEY = LINE_GROUP.ANALYSIS_BY_KEY


# -----------------------------------------------------------------------------
# CLI and path resolution
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the group time-generalization report command-line parser."""
    meg_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Build group-level dot-side time-generalization report."
    )
    parser.add_argument(
        "--config",
        default=str(meg_root / "params" / "pipeline_default.yaml"),
        help="Shared pipeline default YAML.",
    )
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument(
        "--analyses",
        nargs="+",
        default=None,
        help="Analysis keys to include, or 'all' for the four standard analyses.",
    )
    parser.add_argument("--condition", default=None)
    parser.add_argument(
        "--model-source",
        choices=("observed", "predicted"),
        default="observed",
    )
    parser.add_argument(
        "--group-subjects",
        nargs="+",
        default=LINE_GROUP.DEFAULT_SUBJECTS,
        help="Subjects to include, e.g. 01 02 sub-03.",
    )
    parser.add_argument("--metric", choices=SUPPORTED_METRICS, default=DEFAULT_METRIC)
    parser.add_argument("--chance", type=float, default=0.5)
    parser.add_argument("--cluster-forming-p", type=float, default=0.05)
    parser.add_argument("--cluster-alpha", type=float, default=0.05)
    parser.add_argument("--tail", choices=("greater", "two-sided"), default="greater")
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=0,
        help="0 uses exact sign flips when feasible and otherwise 5000 random flips.",
    )
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


def validate_cli_args(args: argparse.Namespace) -> None:
    """Validate report options before reading participant outputs."""
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    if not 0 < args.cluster_forming_p < 1:
        raise SystemExit("--cluster-forming-p must be between 0 and 1.")
    if not 0 < args.cluster_alpha < 1:
        raise SystemExit("--cluster-alpha must be between 0 and 1.")
    if args.n_permutations < 0:
        raise SystemExit("--n-permutations must be non-negative.")


def combined_output_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return output paths for the combined group report."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    report_dir = (
        derivatives_root
        / "reports"
        / "group"
        / session
        / "dot_side_time_generalization"
        / "combined"
    )
    stem = "dot_side_time_generalization_all_analyses_group"
    return {
        "report_dir": report_dir,
        "report_json": report_dir / f"{stem}_report.json",
        "report_html": report_dir / f"{stem}_report.html",
        "report_embedded_html": report_dir / f"{stem}_report_embedded.html",
    }


def group_output_paths(
    cfg: Mapping[str, Any],
    analysis: Mapping[str, str],
) -> dict[str, Path]:
    """Return group-level matrix output paths for one analysis."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    out_dir = (
        derivatives_root
        / "sanity_checks"
        / "group"
        / session
        / "dot_side_time_generalization"
        / analysis["condition"]
        / analysis["model_source"]
    )
    stem = str(analysis["key"])
    return {
        "out_dir": out_dir,
        "figure": out_dir / f"{stem}_dot_side_time_generalization_group_matrix.png",
        "matrix_npz": out_dir / f"{stem}_dot_side_time_generalization_group_matrix.npz",
    }


def subject_output_paths(
    *,
    derivatives_root: Path,
    subject: str,
    session: str,
    condition: str,
    model_source: str,
) -> dict[str, Path]:
    """Return canonical participant time-generalization output paths."""
    out_dir = (
        derivatives_root
        / "sanity_checks"
        / subject
        / session
        / "dot_side_time_generalization"
        / condition
        / model_source
    )
    stem = f"{condition}_{model_source}"
    return {
        "out_dir": out_dir,
        "summary_json": out_dir / f"{stem}_dot_side_time_generalization_summary.json",
        "matrix_npz": out_dir / f"{stem}_dot_side_time_generalization_matrix.npz",
        "figure": out_dir / f"{stem}_dot_side_time_generalization_heatmap.png",
    }


# -----------------------------------------------------------------------------
# Input loading
# -----------------------------------------------------------------------------


def load_subject_matrix(
    *,
    derivatives_root: Path,
    subject: str,
    session: str,
    condition: str,
    model_source: str,
    metric: str,
) -> tuple[dict[str, Any] | None, dict[str, str] | None]:
    """
    Load one participant's time-generalization summary and matrix.

    Returns ``(record, None)`` when loading succeeds, otherwise
    ``(None, skipped_input)`` with a concrete reason for the report.
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
        matrix_npz = Path(str(summary.get("matrix_npz", paths["matrix_npz"]))).expanduser().resolve()
        figure_path = Path(str(summary.get("figure_path", paths["figure"]))).expanduser().resolve()
        if not matrix_npz.is_file():
            raise FileNotFoundError(f"missing matrix NPZ: {matrix_npz}")
        if not figure_path.is_file():
            raise FileNotFoundError(f"missing participant heatmap: {figure_path}")

        with np.load(matrix_npz) as npz:
            if metric not in npz:
                raise RuntimeError(f"matrix NPZ lacks metric {metric!r}")
            matrix = np.asarray(npz[metric], dtype=float)
            times = np.asarray(npz["times_sec"], dtype=float)
            cell_skip_reasons = np.asarray(npz["cell_skip_reasons"])

        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise RuntimeError(f"matrix must be square train x test; got {matrix.shape}")
        if times.ndim != 1 or times.size != matrix.shape[0]:
            raise RuntimeError("times_sec shape does not match matrix axes")
        if cell_skip_reasons.shape != matrix.shape:
            raise RuntimeError("cell_skip_reasons shape does not match matrix")
        if not np.all(np.isfinite(times)):
            raise RuntimeError("times_sec contains non-finite values")

        return {
            "subject": subject,
            "summary_path": str(summary_path),
            "matrix_npz": str(matrix_npz),
            "figure": str(figure_path),
            "summary": summary,
            "times": times,
            "matrix": matrix,
            "cell_skip_reasons": cell_skip_reasons,
        }, None
    except Exception as exc:  # noqa: BLE001 - report the concrete input error.
        return None, {
            "subject": subject,
            "path": str(summary_path),
            "reason": str(exc),
        }


def load_group_inputs(
    *,
    cfg: Mapping[str, Any],
    subjects: Iterable[str],
    condition: str,
    model_source: str,
    metric: str,
    allow_missing: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Load and validate participant matrices for one analysis."""
    derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
    session = str(cfg["layout"]["session_bids_id"])
    records: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    for subject in subjects:
        record, skip = load_subject_matrix(
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
            "Missing or invalid subject time-generalization outputs. "
            f"Use --allow-missing to build a partial report. Examples: {examples}"
        )
    if len(records) < 2:
        raise RuntimeError(
            "At least two valid subject outputs are required for group statistics; "
            f"found {len(records)}."
        )

    reference = records[0]
    reference_times = reference["times"]
    reference_shape = reference["matrix"].shape
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
    ref_tpa = LINE_GROUP.trial_averaging_settings_from_summary(reference_summary)

    for record in records[1:]:
        if record["matrix"].shape != reference_shape:
            raise RuntimeError(
                "Subject matrix shape mismatch between "
                f"{reference['subject']} and {record['subject']}."
            )
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
        cur_tpa = LINE_GROUP.trial_averaging_settings_from_summary(record["summary"])
        if cur_tpa != ref_tpa:
            raise RuntimeError(
                "Subject trial-averaging setting mismatch: "
                f"{reference['subject']}={ref_tpa!r}, {record['subject']}={cur_tpa!r}"
            )

    return records, skipped


def validate_cross_analysis_settings(
    analysis_results: Iterable[Mapping[str, Any]],
) -> None:
    """Require one decoding configuration across all analyses in a report."""
    LINE_GROUP.validate_cross_analysis_settings(analysis_results)


# -----------------------------------------------------------------------------
# 2D cluster permutation
# -----------------------------------------------------------------------------


def one_sample_t_flat(values: np.ndarray) -> np.ndarray:
    """Compute one-sample t values for a subjects x cells matrix."""
    values = np.asarray(values, dtype=float)
    mean = values.mean(axis=0)
    sd = values.std(axis=0, ddof=1)
    n = values.shape[0]
    sem = sd / np.sqrt(float(n))
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
    """Return the t-value threshold used to form 2D clusters."""
    df = n_subjects - 1
    if df < 1:
        raise ValueError("At least two subjects are required.")
    if tail == "greater":
        return float(stats.t.ppf(1.0 - cluster_forming_p, df=df))
    if tail == "two-sided":
        return float(stats.t.ppf(1.0 - (cluster_forming_p / 2.0), df=df))
    raise ValueError(f"Unsupported tail: {tail}")


def connected_components_2d(mask: np.ndarray) -> list[np.ndarray]:
    """Return 4-neighbor connected components as ``n_cells x 2`` coordinates."""
    mask = np.asarray(mask, dtype=bool)
    visited = np.zeros(mask.shape, dtype=bool)
    components: list[np.ndarray] = []
    n_rows, n_cols = mask.shape
    for row in range(n_rows):
        for col in range(n_cols):
            if visited[row, col] or not mask[row, col]:
                continue
            queue: deque[tuple[int, int]] = deque([(row, col)])
            visited[row, col] = True
            cells: list[tuple[int, int]] = []
            while queue:
                cur_row, cur_col = queue.popleft()
                cells.append((cur_row, cur_col))
                for nbr_row, nbr_col in (
                    (cur_row - 1, cur_col),
                    (cur_row + 1, cur_col),
                    (cur_row, cur_col - 1),
                    (cur_row, cur_col + 1),
                ):
                    if not (0 <= nbr_row < n_rows and 0 <= nbr_col < n_cols):
                        continue
                    if visited[nbr_row, nbr_col] or not mask[nbr_row, nbr_col]:
                        continue
                    visited[nbr_row, nbr_col] = True
                    queue.append((nbr_row, nbr_col))
            components.append(np.asarray(cells, dtype=int))
    return components


def cluster_masses_2d(
    t_values: np.ndarray,
    threshold: float,
    tail: str,
) -> list[dict[str, Any]]:
    """Find 2D clusters and compute mass from t values."""
    if tail == "greater":
        mask = t_values > threshold
        mass_source = t_values
    elif tail == "two-sided":
        mask = np.abs(t_values) > threshold
        mass_source = np.abs(t_values)
    else:
        raise ValueError(f"Unsupported tail: {tail}")
    mask = np.asarray(mask, dtype=bool)
    # Infinite t values can occur when every subject has the same positive
    # effect. The 1D dot-side cluster code treats those as valid supra-threshold
    # points, so only NaN gaps are removed from 2D cluster formation.
    mask[np.isnan(t_values)] = False

    clusters: list[dict[str, Any]] = []
    for component in connected_components_2d(mask):
        rows = component[:, 0]
        cols = component[:, 1]
        mass = float(np.nansum(mass_source[rows, cols]))
        clusters.append(
            {
                "train_start_index": int(rows.min()),
                "train_end_index": int(rows.max()),
                "test_start_index": int(cols.min()),
                "test_end_index": int(cols.max()),
                "n_cells": int(component.shape[0]),
                "mass": mass,
                "cell_indices": [(int(r), int(c)) for r, c in component.tolist()],
            }
        )
    return clusters


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


def cluster_permutation_2d(
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
    Test mean matrix against chance using 2D sign-flip cluster permutation.

    ``values`` is subjects x train_time x test_time. Only cells that are
    finite for every included subject enter the statistic.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError("values must be subject x train_time x test_time.")
    if values.shape[0] < 2:
        raise ValueError("At least two subjects are required.")

    valid_cell_mask = np.all(np.isfinite(values), axis=0)
    if not np.any(valid_cell_mask):
        raise ValueError(
            "Group metric contains no train/test cells with finite values for "
            "every included subject."
        )

    diffs_valid = (values[:, valid_cell_mask] - chance).astype(float)
    n_subjects = int(values.shape[0])
    threshold = cluster_threshold(
        n_subjects=n_subjects,
        cluster_forming_p=cluster_forming_p,
        tail=tail,
    )
    observed_t = np.full(values.shape[1:], np.nan, dtype=float)
    observed_t[valid_cell_mask] = one_sample_t_flat(diffs_valid)
    observed_clusters = cluster_masses_2d(observed_t, threshold, tail)

    signs, mode = _sign_matrix(
        n_subjects=n_subjects,
        n_permutations=n_permutations,
        random_seed=random_seed,
    )
    null_max_masses = np.zeros(signs.shape[0], dtype=float)
    for perm_idx, sign_row in enumerate(signs):
        permuted = diffs_valid * sign_row[:, np.newaxis]
        perm_t = np.full(values.shape[1:], np.nan, dtype=float)
        perm_t[valid_cell_mask] = one_sample_t_flat(permuted)
        perm_clusters = cluster_masses_2d(perm_t, threshold, tail)
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
        "valid_cell_mask": valid_cell_mask,
        "n_valid_cells": int(valid_cell_mask.sum()),
        "n_invalid_cells": int((~valid_cell_mask).sum()),
        "clusters": observed_clusters,
        "null_max_cluster_masses": null_max_masses,
    }


def group_matrix_summary(
    values: np.ndarray,
    *,
    valid_cell_mask: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return group mean/SEM matrices while preserving invalid-cell gaps."""
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values, axis=0)
    sem = np.nanstd(values, axis=0, ddof=1) / np.sqrt(float(values.shape[0]))
    mean = np.where(valid_cell_mask, mean, np.nan)
    sem = np.where(valid_cell_mask, sem, np.nan)
    return {"mean": mean, "sem": sem}


# -----------------------------------------------------------------------------
# Plotting and report writing
# -----------------------------------------------------------------------------


def significant_cluster_mask(
    clusters: Iterable[Mapping[str, Any]],
    shape: tuple[int, int],
) -> np.ndarray:
    """Return a boolean matrix covering significant 2D cluster cells."""
    mask = np.zeros(shape, dtype=bool)
    for cluster in clusters:
        if not cluster.get("significant"):
            continue
        for row, col in cluster.get("cell_indices", []):
            mask[int(row), int(col)] = True
    return mask


def plot_group_matrix(
    *,
    times: np.ndarray,
    group_mean: np.ndarray,
    stats_result: Mapping[str, Any],
    subjects: list[str],
    metric: str,
    condition: str,
    model_source: str,
    occlusion_timing: Mapping[str, Any],
    out_file: Path,
    raw_root: Path,
) -> Path:
    """Write the group mean matrix with 2D cluster overlays."""
    out_path = LINE_GROUP.guarded_path(out_file, raw_root)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    finite_values = group_mean[np.isfinite(group_mean)]
    if finite_values.size == 0:
        raise RuntimeError("Cannot plot a group matrix with no finite cells.")

    vmin = min(0.35, float(np.nanmin(finite_values)))
    vmax = max(0.75, float(np.nanmax(finite_values)))
    if vmin >= 0.5:
        vmin = 0.35
    if vmax <= 0.5:
        vmax = 0.75
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#d9d9d9")
    norm = TwoSlopeNorm(vmin=vmin, vcenter=float(stats_result["chance"]), vmax=vmax)

    fig, ax = plt.subplots(figsize=(10, 8.5), constrained_layout=True)
    im = ax.imshow(
        np.ma.masked_invalid(group_mean),
        origin="lower",
        interpolation="nearest",
        aspect="auto",
        extent=[float(times[0]), float(times[-1]), float(times[0]), float(times[-1])],
        cmap=cmap,
        norm=norm,
    )

    sig_mask = significant_cluster_mask(stats_result["clusters"], group_mean.shape)
    if np.any(sig_mask):
        overlay = np.ma.masked_where(~sig_mask, sig_mask.astype(float))
        ax.imshow(
            overlay,
            origin="lower",
            interpolation="nearest",
            aspect="auto",
            extent=[float(times[0]), float(times[-1]), float(times[0]), float(times[-1])],
            cmap=ListedColormap(["#d62728"]),
            alpha=0.18,
        )
        x_grid, y_grid = np.meshgrid(times, times)
        ax.contour(
            x_grid,
            y_grid,
            sig_mask.astype(float),
            levels=[0.5],
            colors="#d62728",
            linewidths=1.2,
        )
        for cluster in [c for c in stats_result["clusters"] if c.get("significant")][:8]:
            cells = np.asarray(cluster["cell_indices"], dtype=int)
            train_time = float(np.mean(times[cells[:, 0]]))
            test_time = float(np.mean(times[cells[:, 1]]))
            ax.text(
                test_time,
                train_time,
                f"p={cluster['p_value']:.3f}",
                color="#d62728",
                ha="center",
                va="center",
                fontsize=8,
                bbox={
                    "boxstyle": "round",
                    "facecolor": "white",
                    "alpha": 0.8,
                    "edgecolor": "#d62728",
                },
            )

    ax.plot([times[0], times[-1]], [times[0], times[-1]], color="white", lw=1.2, alpha=0.9)
    ax.axvline(0.0, color="k", lw=1.0, ls="--", alpha=0.7)
    ax.axhline(0.0, color="k", lw=1.0, ls="--", alpha=0.7)
    TIME_GEN.add_occlusion_boundary_lines_2d(
        ax,
        condition=condition,
        occlusion_timing=occlusion_timing,
    )
    threshold_text = (
        f"2D cluster threshold: p={stats_result['cluster_forming_p']:.3g}, "
        f"t>{stats_result['t_threshold']:.3f}; "
        f"alpha={stats_result['cluster_alpha']:.3g}; "
        f"{stats_result['permutation_mode']} sign flips n={stats_result['n_permutations']}; "
        f"valid cells={stats_result['n_valid_cells']}"
    )
    ax.text(
        0.01,
        0.99,
        threshold_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "0.75"},
    )
    ax.set_xlabel("Test time from condition onset (s)")
    ax.set_ylabel("Train time from condition onset (s)")
    ax.set_title(
        f"Group dot-side time generalization | {condition} | {model_source} model | "
        f"n={len(subjects)}"
    )
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="lower right", frameon=True)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"Group mean {metric}")
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def write_group_npz(
    path: Path,
    *,
    times: np.ndarray,
    subject_values: np.ndarray,
    group_mean: np.ndarray,
    group_sem: np.ndarray,
    stats_result: Mapping[str, Any],
    metric: str,
) -> Path:
    """Write the group matrix arrays and statistic grids."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        times_sec=np.asarray(times, dtype=float),
        subject_values=np.asarray(subject_values, dtype=float),
        group_mean=np.asarray(group_mean, dtype=float),
        group_sem=np.asarray(group_sem, dtype=float),
        observed_t=np.asarray(stats_result["observed_t"], dtype=float),
        valid_cell_mask=np.asarray(stats_result["valid_cell_mask"], dtype=bool),
        metric=np.asarray(metric),
    )
    return path


def _image_src(report_path: Path, image_path: str | Path, *, embed: bool) -> str:
    """Return a relative or data-URI image source for HTML reports."""
    resolved = Path(image_path).expanduser().resolve()
    if embed:
        payload = base64.b64encode(resolved.read_bytes()).decode("ascii")
        return f"data:image/png;base64,{payload}"
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


def _cluster_for_report(
    cluster: Mapping[str, Any],
    times: np.ndarray,
) -> dict[str, Any]:
    """Drop bulky cell lists while keeping interpretable cluster bounds."""
    return {
        "train_start_index": int(cluster["train_start_index"]),
        "train_end_index": int(cluster["train_end_index"]),
        "test_start_index": int(cluster["test_start_index"]),
        "test_end_index": int(cluster["test_end_index"]),
        "train_start_time_sec": float(times[int(cluster["train_start_index"])]),
        "train_end_time_sec": float(times[int(cluster["train_end_index"])]),
        "test_start_time_sec": float(times[int(cluster["test_start_index"])]),
        "test_end_time_sec": float(times[int(cluster["test_end_index"])]),
        "n_cells": int(cluster["n_cells"]),
        "mass": float(cluster["mass"]),
        "p_value": float(cluster["p_value"]),
        "significant": bool(cluster["significant"]),
    }


def _cluster_table_html(clusters: Iterable[Mapping[str, Any]]) -> str:
    """Build a 2D cluster table for one group analysis."""
    rows = "\n".join(
        "<tr>"
        f"<td>{cluster['train_start_time_sec']:.3f}..{cluster['train_end_time_sec']:.3f}</td>"
        f"<td>{cluster['test_start_time_sec']:.3f}..{cluster['test_end_time_sec']:.3f}</td>"
        f"<td>{cluster['n_cells']}</td>"
        f"<td>{cluster['mass']:.3f}</td>"
        f"<td>{cluster['p_value']:.4f}</td>"
        f"<td>{'yes' if cluster['significant'] else 'no'}</td>"
        "</tr>"
        for cluster in clusters
    ) or "<tr><td colspan=\"6\">No clusters formed.</td></tr>"
    return (
        "<table><thead><tr><th>train s</th><th>test s</th><th>n cells</th>"
        "<th>mass</th><th>p</th><th>significant</th></tr></thead><tbody>"
        f"{rows}</tbody></table>"
    )


def _skipped_inputs_html(skipped_inputs: Iterable[Mapping[str, str]]) -> str:
    """Build a compact disclosure for missing or invalid inputs."""
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


def _subject_diagnostics_html(summary: Mapping[str, Any]) -> str:
    """Render participant matrix diagnostics for the combined report."""
    return (
        "<h4>Training/testing diagnostics</h4>"
        + TIME_GEN.time_generalization_diagnostics_html(summary)
    )


def _analysis_diagnostics_html(result: Mapping[str, Any]) -> str:
    """Summarize group cell validity and participant diagnostics."""
    records = list(result["participant_records"])
    group_rows = [
        ("Subjects included", ", ".join(result["subjects_included"])),
        ("Metric", result["metric"]),
        (
            "Group-valid cells",
            f"{result['n_group_valid_cells']} valid, {result['n_group_invalid_cells']} invalid",
        ),
        ("Cluster-forming t", f"{result['cluster_forming_t_threshold']:.3f}"),
        ("Permutation mode", result["permutation_mode"]),
        ("Permutations", result["n_permutations"]),
    ]
    group_table = "\n".join(
        "<tr><th>"
        + html.escape(str(label))
        + "</th><td>"
        + html.escape(str(value))
        + "</td></tr>"
        for label, value in group_rows
    )

    participant_rows: list[str] = []
    for record in records:
        diagnostics = record["summary"].get("time_generalization_diagnostics", {})
        sample = diagnostics.get("sample_structure", {})
        cell_validity = diagnostics.get("cell_validity", {})
        assessment = diagnostics.get("assessment", {})
        frac = cell_validity.get("valid_balanced_cell_fraction")
        pct = "n/a" if frac is None else f"{float(frac) * 100.0:.1f}%"
        participant_rows.append(
            "<tr>"
            f"<td>{html.escape(record['subject'])}</td>"
            f"<td>{html.escape(str(sample.get('averaging_mode', 'n/a')))}</td>"
            f"<td>{html.escape(str(sample.get('requested_value', 'n/a')))}</td>"
            f"<td>{html.escape(str(sample.get('n_decoder_trials', 'n/a')))}</td>"
            f"<td>{html.escape(str(cell_validity.get('n_valid_balanced_cells', 'n/a')))} / "
            f"{html.escape(str(cell_validity.get('n_cells', 'n/a')))} ({html.escape(pct)})</td>"
            f"<td>{html.escape(str(assessment.get('status', 'n/a')))}</td>"
            "</tr>"
        )
    return (
        "<table><tbody>"
        + group_table
        + "</tbody></table>"
        "<table><thead><tr><th>subject</th><th>averaging mode</th>"
        "<th>request</th><th>pseudo-trials</th><th>valid cells</th>"
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
    raw_root: Path,
    embed_images: bool,
) -> Path:
    """Write one report with all group matrices before participant sections."""
    report_path = LINE_GROUP.guarded_path(path, raw_root)
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
                            "Group mean train-time x test-time matrix. Red "
                            "overlays mark significant 2D cluster-permutation "
                            "regions."
                        ),
                        embed=embed_images,
                    ),
                    "<h4>2D clusters</h4>",
                    _cluster_table_html(result["clusters"]),
                    "<h4>Diagnostics</h4>",
                    _analysis_diagnostics_html(result),
                    _skipped_inputs_html(result["skipped_inputs"]),
                    "</article>",
                ]
            )
        )

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
                        f"peak balanced={summary.get('peak_accuracy_balanced_train_test')} "
                        f"at train {summary.get('peak_train_time_sec_balanced')} s, "
                        f"test {summary.get('peak_test_time_sec_balanced')} s"
                    ),
                    embed=embed_images,
                )
                content += _subject_diagnostics_html(summary)
            analysis_blocks.append(
                f"<article><h3>{html.escape(result['label'])}</h3>{content}</article>"
            )
        participant_sections.append(
            f"<section id=\"{html.escape(subject)}\"><h2>{html.escape(subject)}</h2>"
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
            "<section id=\"participants\"><h2>Participant-specific matrices</h2>",
            "\n".join(participant_sections),
            "</section>",
            "</body></html>",
        ]
    )
    report_path.write_text(body, encoding="utf-8")
    return report_path


def _jsonify(value: Any) -> Any:
    """Convert NumPy-heavy report objects to JSON-serializable structures."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _jsonify(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(item) for item in value]
    return value


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def build_group_report(cfg: Mapping[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    """Load analyses, run 2D cluster permutation, and write the group report."""
    raw_root = Path(cfg["resolved_paths"]["data_root"])
    subjects = [LINE_GROUP.normalize_subject_label(value) for value in args.group_subjects]
    analysis_specs = LINE_GROUP.resolve_analysis_specs(args)
    output_paths = combined_output_paths(cfg)
    expected_outputs = [
        output_paths["report_json"],
        output_paths["report_html"],
        output_paths["report_embedded_html"],
        *[group_output_paths(cfg, spec)["figure"] for spec in analysis_specs],
        *[group_output_paths(cfg, spec)["matrix_npz"] for spec in analysis_specs],
    ]

    if args.resume and all(path.is_file() for path in expected_outputs):
        existing_report = json.loads(output_paths["report_json"].read_text(encoding="utf-8"))
        requested_keys = [spec["key"] for spec in analysis_specs]
        if existing_report.get("analysis_keys") != requested_keys:
            raise RuntimeError(
                "Existing time-generalization report analysis keys do not "
                f"match requested selection: existing={existing_report.get('analysis_keys')}, "
                f"requested={requested_keys}. Use --overwrite to rebuild."
            )
        return {
            "stage": "dot_side_time_generalization_group_report",
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
                "Refusing to overwrite existing time-generalization group "
                f"outputs without --overwrite or --resume: {existing}"
            )

    analysis_results: list[dict[str, Any]] = []
    occlusion_timing = LINE_GROUP.occlusion_timing_from_config(cfg)
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
        values = np.stack([record["matrix"] for record in records], axis=0)
        valid_subjects = [record["subject"] for record in records]
        stats_result = cluster_permutation_2d(
            values,
            chance=args.chance,
            cluster_forming_p=args.cluster_forming_p,
            cluster_alpha=args.cluster_alpha,
            tail=args.tail,
            n_permutations=args.n_permutations,
            random_seed=args.random_seed,
        )
        summary = group_matrix_summary(
            values,
            valid_cell_mask=np.asarray(stats_result["valid_cell_mask"], dtype=bool),
        )
        clusters_for_report = [
            _cluster_for_report(cluster, times) for cluster in stats_result["clusters"]
        ]
        paths = group_output_paths(cfg, spec)
        group_npz = paths["matrix_npz"]
        group_figure = paths["figure"]
        if not args.validate_only:
            group_npz = write_group_npz(
                group_npz,
                times=times,
                subject_values=values,
                group_mean=summary["mean"],
                group_sem=summary["sem"],
                stats_result=stats_result,
                metric=args.metric,
            )
            group_figure = plot_group_matrix(
                times=times,
                group_mean=summary["mean"],
                stats_result=stats_result,
                subjects=valid_subjects,
                metric=args.metric,
                condition=spec["condition"],
                model_source=spec["model_source"],
                occlusion_timing=occlusion_timing,
                out_file=group_figure,
                raw_root=raw_root,
            )

        analysis_results.append(
            {
                **spec,
                "metric": args.metric,
                "subjects_requested": subjects,
                "subjects_included": valid_subjects,
                "skipped_inputs": skipped,
                "times_sec": times,
                "group_mean": summary["mean"],
                "group_sem": summary["sem"],
                "group_valid_cell_mask": stats_result["valid_cell_mask"],
                "n_group_valid_cells": int(stats_result["n_valid_cells"]),
                "n_group_invalid_cells": int(stats_result["n_invalid_cells"]),
                "clusters": clusters_for_report,
                "clusters_with_cells": stats_result["clusters"],
                "cluster_forming_t_threshold": float(stats_result["t_threshold"]),
                "permutation_mode": str(stats_result["permutation_mode"]),
                "n_permutations": int(stats_result["n_permutations"]),
                "group_figure": str(group_figure),
                "group_matrix_npz": str(group_npz),
                "participant_figures": [
                    {
                        "subject": record["subject"],
                        "figure": record["figure"],
                        "matrix_npz": record["matrix_npz"],
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
            "stage": "dot_side_time_generalization_group_report",
            "status": "validated",
            "analysis_keys": [result["key"] for result in analysis_results],
            "metric": args.metric,
            "analyses": {
                result["key"]: {
                    "subjects_requested": result["subjects_requested"],
                    "subjects_included": result["subjects_included"],
                    "skipped_inputs": result["skipped_inputs"],
                    "n_times": int(np.asarray(result["times_sec"]).size),
                    "n_group_valid_cells": result["n_group_valid_cells"],
                    "n_group_invalid_cells": result["n_group_invalid_cells"],
                    "cluster_count": len(result["clusters"]),
                    "significant_cluster_count": int(
                        sum(1 for cluster in result["clusters"] if cluster["significant"])
                    ),
                }
                for result in analysis_results
            },
        }

    summary_lines = [
        "Stage status: completed",
        f"Analyses: {', '.join(result['key'] for result in analysis_results)}",
        f"Metric: {args.metric}",
        f"Subjects requested: {', '.join(subjects)}",
        f"Chance: {args.chance:g}",
        f"Tail: {args.tail}",
        f"2D cluster-forming p: {args.cluster_forming_p:g}",
        f"Cluster alpha: {args.cluster_alpha:g}",
        f"Machine report: {output_paths['report_json']}",
    ]
    report = {
        "stage": "dot_side_time_generalization_group_report",
        "status": "completed",
        "generated_utc": LINE_GROUP.utc_now(),
        "analysis_keys": [result["key"] for result in analysis_results],
        "metric": args.metric,
        "session_bids_id": str(cfg["layout"]["session_bids_id"]),
        "subjects_requested": subjects,
        "chance": float(args.chance),
        "tail": args.tail,
        "cluster_forming_p": float(args.cluster_forming_p),
        "cluster_alpha": float(args.cluster_alpha),
        "occlusion_timing": occlusion_timing,
        "analyses": {
            result["key"]: {
                key: value
                for key, value in result.items()
                if key not in {"participant_records", "clusters_with_cells"}
            }
            for result in analysis_results
        },
        "report_files": {
            "json": str(output_paths["report_json"]),
            "html": str(output_paths["report_html"]),
            "embedded_html": str(output_paths["report_embedded_html"]),
        },
    }
    output_paths["report_json"].parent.mkdir(parents=True, exist_ok=True)
    output_paths["report_json"].write_text(
        json.dumps(_jsonify(report), indent=2),
        encoding="utf-8",
    )
    title = (
        f"Group | {cfg['layout']['session_bids_id']} | "
        "Dot-side time generalization | all analyses"
    )
    write_combined_html_report(
        output_paths["report_html"],
        title=title,
        summary_lines=summary_lines,
        subjects=subjects,
        analysis_results=analysis_results,
        raw_root=raw_root,
        embed_images=False,
    )
    write_combined_html_report(
        output_paths["report_embedded_html"],
        title=title,
        summary_lines=summary_lines,
        subjects=subjects,
        analysis_results=analysis_results,
        raw_root=raw_root,
        embed_images=True,
    )
    return report


def main() -> int:
    """Resolve configuration and build or validate the group report."""
    args = build_arg_parser().parse_args()
    validate_cli_args(args)
    cfg, warnings = LINE_GROUP.resolve_config(args)
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
                    "n_group_valid_cells": value.get("n_group_valid_cells"),
                    "significant_cluster_count": sum(
                        1
                        for cluster in value.get("clusters", [])
                        if cluster.get("significant")
                    ),
                    "group_figure": value.get("group_figure"),
                }
                for key, value in report["analyses"].items()
            }
        print(json.dumps(_jsonify(printable), indent=2))
    return 0 if report.get("status") in {"validated", "completed", "completed_resume"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
