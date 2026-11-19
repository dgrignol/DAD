#!/usr/bin/env python3
"""
03_dot_side_time_generalization.py

Purpose
Run a MEG dot-side time-generalization sanity check for one participant,
neural condition, and trajectory model source. The workflow trains the same
regularized LDA classifier used by ``01_dot_side_decoding.py`` at each
centered sliding-window train time and tests that classifier at every other
test time, producing a train-time x test-time accuracy matrix.

Usage examples
1) Default run for sub-04, ses-01, always-visible observed labels:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/03_dot_side_time_generalization.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-04.yaml \
     --subject 01 --condition always_visible --model-source observed --overwrite

2) Run occluded-deviant epochs against their expected/predicted path:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/03_dot_side_time_generalization.py \
     --subject 05 --condition occluded_deviant --model-source predicted \
     --decode-window-samples 10 --trials-per-average all_trials_per_identity \
     --overwrite

3) Validate inputs and output paths without decoding:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/sanity_checks/03_dot_side_time_generalization.py \
     --subject 01 --condition occluded_nondeviant --validate-only --print-summary

Inputs
- Completed Stage-04 epochs FIF plus ``source_manifest.json``.
- Subject-specific observed trajectory MAT file containing the requested
  observed or expected/predicted trajectory condition.
- The same configuration files consumed by ``01_dot_side_decoding.py``.

Outputs
- NPZ matrix bundle with times, accuracy matrices, valid masks, fold counts,
  and skip-reason grids.
- Long-form matrix CSV with one row per train/test cell.
- PNG heatmap, with occlusion/reappearance markers on both train and test axes
  for occluded conditions.
- Summary JSON with settings, trial-averaging metadata, class-balance
  diagnostics, valid/invalid cell counts, and peak/mean accuracies.
- MNE HTML report with the heatmap and diagnostics.

Key assumptions
- Condition/model-source semantics, dynamic dot-side labels, and same-identity
  trial averaging are inherited from ``01_dot_side_decoding.py``.
- Cross-validation splits are defined once per train time from the train-time
  labels and then reused across every test time for that train row.
- Balanced decoding trains on a deterministic undersampled subset of each
  train fold and evaluates only when every held-out fold can be balanced at
  the test time. Cells that cannot support balanced evaluation are explicit
  NaNs and carry a skip reason.
- Time-generalization requires matching feature dimensions between train and
  test windows. Edge-clipped windows therefore only generalize to test windows
  with the same effective sample count; other cells are NaN with a
  ``feature_dimension_mismatch`` reason.
"""

from __future__ import annotations

import argparse
from collections import Counter
import html
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any, Mapping

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import mne
import numpy as np
import pandas as pd


def _load_sibling_module(module_name: str, path: Path) -> Any:
    """Import a sibling script whose filename is not a valid Python module id."""
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {module_name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


DOT_SIDE = _load_sibling_module(
    "dot_side_decoding_helpers",
    Path(__file__).with_name("01_dot_side_decoding.py"),
)


BALANCED_METRIC = "accuracy_balanced_undersampled_sets"
UNBALANCED_METRIC = "accuracy_unbalanced_train"
MATRIX_STAGE = "dot_side_time_generalization"


# -----------------------------------------------------------------------------
# CLI and path resolution
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the participant time-generalization command-line parser."""
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Dot-side train-time x test-time decoding sanity check."
    )
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
    parser.add_argument("--subject", type=str, default="04")
    parser.add_argument("--session", type=str, default="01")
    parser.add_argument("--condition", type=str, default="always_visible")
    parser.add_argument(
        "--model-source",
        choices=DOT_SIDE.MODEL_SOURCES,
        default="observed",
        help=(
            "Trajectory source for dynamic labels. Predicted is valid only for "
            "occluded_deviant and uses occluded_nondeviant trajectories."
        ),
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
    parser.add_argument("--epochs-fif", type=Path, default=None)
    parser.add_argument("--model-mat", type=Path, default=None)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Optional explicit output directory. If omitted, uses "
            "derivatives/sanity_checks/sub-XX/ses-YY/dot_side_time_generalization."
        ),
    )
    parser.add_argument("--decode-sfreq", type=float, default=120.0)
    parser.add_argument("--decode-tmin", type=float, default=0.0)
    parser.add_argument("--decode-tmax", type=float, default=2.667)
    parser.add_argument(
        "--decode-window-samples",
        type=int,
        default=10,
        help="Centered sliding-window width in samples after resampling.",
    )
    parser.add_argument(
        "--trials-per-average",
        type=DOT_SIDE.parse_trials_per_average,
        default=1,
        help=(
            "Positive integer same-identity averaging group size, or the "
            "literal all_trials_per_identity."
        ),
    )
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=None,
        help="Optional source_manifest.json used to validate Stage-04 lineage.",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-summary", action="store_true")
    return parser


def validate_cli_args(args: argparse.Namespace) -> None:
    """Fail early for invalid command-line combinations and numeric settings."""
    DOT_SIDE.validate_cli_args(args)


def time_generalization_output_paths(
    out_dir: Path,
    condition: str,
    model_source: str,
) -> dict[str, Path]:
    """Return all expected participant time-generalization outputs."""
    stem = DOT_SIDE.analysis_key(condition, model_source)
    return {
        "summary_json": out_dir / f"{stem}_dot_side_time_generalization_summary.json",
        "matrix_npz": out_dir / f"{stem}_dot_side_time_generalization_matrix.npz",
        "matrix_csv": out_dir / f"{stem}_dot_side_time_generalization_matrix.csv",
        "figure": out_dir / f"{stem}_dot_side_time_generalization_heatmap.png",
        "report_html": out_dir / f"{stem}_dot_side_time_generalization_report.html",
    }


def resolve_input_output_paths(
    args: argparse.Namespace,
) -> tuple[Path, Path, Path, Path, dict[str, Any]]:
    """
    Resolve epochs/model inputs plus the time-generalization output directory.

    The existing dot-side path resolver is reused for input validation and
    configuration loading. Its default output directory is then replaced with
    this workflow's separate ``dot_side_time_generalization`` namespace.
    """
    epochs_fif, model_mat, _line_out_dir, source_manifest, cfg = (
        DOT_SIDE.resolve_input_output_paths(args)
    )
    if args.out_dir is None:
        derivatives_root = Path(cfg["resolved_paths"]["derivatives_root"])
        subject = str(cfg["subject"]["id"])
        session = str(cfg["layout"]["session_bids_id"])
        out_dir = (
            derivatives_root
            / "sanity_checks"
            / f"sub-{subject}"
            / session
            / "dot_side_time_generalization"
            / str(args.condition)
            / str(args.model_source)
        )
    else:
        out_dir = Path(args.out_dir).expanduser().resolve()

    DOT_SIDE.assert_output_path_outside_raw(out_dir, cfg["resolved_paths"]["data_root"])
    return epochs_fif, model_mat, out_dir, source_manifest, cfg


# -----------------------------------------------------------------------------
# LDA helpers
# -----------------------------------------------------------------------------


def fit_regularized_lda(
    x_train: np.ndarray,
    y_train: np.ndarray,
    *,
    reg: float = 1e-2,
) -> dict[str, np.ndarray | float]:
    """
    Fit the same train-fold z-scored regularized binary LDA used in script 01.

    A lightweight fit/predict split is used here so one train-time model can be
    reused across all test times in a fold without refitting from test-time
    information.
    """
    y_train = np.asarray(y_train, dtype=int)
    x_train = np.asarray(x_train, dtype=float)
    x0 = x_train[y_train == 0]
    x1 = x_train[y_train == 1]
    if x0.size == 0 or x1.size == 0:
        raise ValueError("Regularized LDA requires both classes in training data.")

    mu = x_train.mean(axis=0, keepdims=True)
    sigma = x_train.std(axis=0, keepdims=True)
    sigma[sigma == 0.0] = 1.0
    x_tr = (x_train - mu) / sigma

    x0 = x_tr[y_train == 0]
    x1 = x_tr[y_train == 1]
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
        # Match the Woodbury form from script 01 for high-dimensional windows.
        gram = scaled_centered @ scaled_centered.T
        middle = np.eye(gram.shape[0], dtype=float) + (gram / reg)
        rhs = (scaled_centered @ mean_delta) / reg
        correction = scaled_centered.T @ np.linalg.solve(middle, rhs)
        w = (mean_delta / reg) - (correction / reg)

    b = -0.5 * float((m1 + m0) @ w)
    return {"mu": mu, "sigma": sigma, "w": w, "b": b}


def predict_regularized_lda(
    model: Mapping[str, np.ndarray | float],
    x_test: np.ndarray,
) -> np.ndarray:
    """Predict binary labels from a fitted regularized LDA model."""
    x_test = np.asarray(x_test, dtype=float)
    mu = np.asarray(model["mu"], dtype=float)
    sigma = np.asarray(model["sigma"], dtype=float)
    w = np.asarray(model["w"], dtype=float)
    b = float(model["b"])
    x_te = (x_test - mu) / sigma
    score = x_te @ w + b
    return (score > 0.0).astype(int)


# -----------------------------------------------------------------------------
# Time-generalization decoding
# -----------------------------------------------------------------------------


def _class_counts(y: np.ndarray) -> dict[str, int]:
    """Return stable JSON-ready binary class counts."""
    counts = Counter(np.asarray(y, dtype=int).tolist())
    return {str(int(key)): int(value) for key, value in sorted(counts.items())}


def _feature_matrices(
    data: np.ndarray,
    window_rows: list[dict[str, Any]],
) -> list[np.ndarray]:
    """Flatten each centered channel x sample window into trial feature rows."""
    n_trials = int(data.shape[0])
    matrices: list[np.ndarray] = []
    for row in window_rows:
        start_idx = int(row["window_start_sample_inclusive"])
        end_idx = int(row["window_end_sample_exclusive"])
        matrices.append(data[:, :, start_idx:end_idx].reshape(n_trials, -1))
    return matrices


def _train_balance_seed(random_seed: int, train_idx: int, fold_idx: int) -> int:
    """Return a deterministic seed for train-fold undersampling."""
    return int((random_seed * 1000003) + 300000 + (train_idx * 1009) + fold_idx)


def _test_balance_seed(
    random_seed: int,
    train_idx: int,
    test_idx: int,
    fold_idx: int,
) -> int:
    """Return a deterministic seed for held-out test-fold undersampling."""
    return int(
        (random_seed * 1000003)
        + 500000
        + (train_idx * 1009)
        + (test_idx * 997)
        + fold_idx
    )


def decode_time_generalization(
    data: np.ndarray,
    labels: np.ndarray,
    *,
    times: np.ndarray,
    n_splits: int,
    window_samples: int,
    random_seed: int,
    sfreq: float,
) -> dict[str, Any]:
    """
    Decode train-time x test-time dynamic dot-side labels.

    Data flow
    - Build the same centered window metadata used by the 1D decoder.
    - For each train time, validate train labels and build one stratified
      K-fold partition from those labels.
    - Fit one unbalanced and one balanced LDA model per train fold.
    - Reuse those train-fold models across every test time, marking cells NaN
      when test features have incompatible dimensionality or held-out test
      labels cannot be balanced in every fold.
    """
    n_trials, _n_channels, n_times = data.shape
    if labels.shape != (n_trials, n_times):
        raise ValueError(f"labels shape mismatch: got {labels.shape}, expected {(n_trials, n_times)}")
    if times.shape != (n_times,):
        raise ValueError(f"times shape mismatch: got {times.shape}, expected {(n_times,)}")

    window_rows = DOT_SIDE.build_centered_window_metadata(
        times=times,
        window_samples=window_samples,
        sfreq=sfreq,
    )
    features = _feature_matrices(data, window_rows)
    feature_dims = np.asarray([matrix.shape[1] for matrix in features], dtype=int)

    accuracy_unbalanced = np.full((n_times, n_times), np.nan, dtype=float)
    accuracy_balanced = np.full((n_times, n_times), np.nan, dtype=float)
    n_valid_unbalanced_folds = np.zeros((n_times, n_times), dtype=int)
    n_valid_balanced_folds = np.zeros((n_times, n_times), dtype=int)
    cell_skip_reasons = np.full((n_times, n_times), "", dtype=object)
    train_rows: list[dict[str, Any]] = []
    test_class_counts = [_class_counts(labels[:, idx]) for idx in range(n_times)]

    for train_idx, base_row in enumerate(window_rows):
        row = dict(base_row)
        row["train_time_sec"] = float(times[train_idx])
        row["train_feature_dim"] = int(feature_dims[train_idx])
        y_train_time = labels[:, train_idx]
        train_counts = _class_counts(y_train_time)
        row["train_class_counts"] = train_counts
        row["train_skip_reason"] = ""

        if np.unique(y_train_time).size < 2:
            row["train_skip_reason"] = "single_class_train_labels"
        else:
            min_class_count = min(train_counts.values())
            if min_class_count < n_splits:
                row["train_skip_reason"] = (
                    f"too_few_train_class_samples_for_{n_splits}_fold_stratified_cv"
                )

        if row["train_skip_reason"]:
            cell_skip_reasons[train_idx, :] = row["train_skip_reason"]
            row["n_valid_balanced_test_cells"] = 0
            row["n_invalid_balanced_test_cells"] = int(n_times)
            row["n_feature_dimension_matched_test_cells"] = 0
            train_rows.append(row)
            continue

        split_indices = DOT_SIDE.stratified_kfold_indices(
            y_train_time,
            n_splits=n_splits,
            random_state=random_seed,
        )
        unbalanced_models: list[dict[str, np.ndarray | float]] = []
        balanced_models: list[dict[str, np.ndarray | float]] = []
        usable_splits: list[tuple[np.ndarray, np.ndarray]] = []
        x_train_time = features[train_idx]

        for fold_idx, (train_fold_idx, test_fold_idx) in enumerate(split_indices):
            x_fold = x_train_time[train_fold_idx]
            y_fold = y_train_time[train_fold_idx]
            unbalanced_models.append(fit_regularized_lda(x_fold, y_fold))

            rng = np.random.default_rng(
                _train_balance_seed(random_seed, train_idx, fold_idx)
            )
            train_bal_idx = DOT_SIDE.undersample_majority_to_minority_indices(
                y_fold,
                rng=rng,
            )
            if train_bal_idx.size == 0:
                row["train_skip_reason"] = "no_valid_train_folds_after_balancing"
                break
            balanced_models.append(
                fit_regularized_lda(x_fold[train_bal_idx], y_fold[train_bal_idx])
            )
            usable_splits.append((train_fold_idx, test_fold_idx))

        if row["train_skip_reason"]:
            cell_skip_reasons[train_idx, :] = row["train_skip_reason"]
            row["n_valid_balanced_test_cells"] = 0
            row["n_invalid_balanced_test_cells"] = int(n_times)
            row["n_feature_dimension_matched_test_cells"] = 0
            train_rows.append(row)
            continue

        matched_dim_count = 0
        for test_idx in range(n_times):
            if feature_dims[test_idx] != feature_dims[train_idx]:
                cell_skip_reasons[train_idx, test_idx] = "feature_dimension_mismatch"
                continue
            matched_dim_count += 1

            x_test_time = features[test_idx]
            y_test_time = labels[:, test_idx]
            fold_acc_unbalanced: list[float] = []
            fold_acc_balanced: list[float] = []

            for fold_idx, (_train_fold_idx, test_fold_idx) in enumerate(usable_splits):
                y_test_fold = y_test_time[test_fold_idx]
                x_test_fold = x_test_time[test_fold_idx]

                y_pred_unbalanced = predict_regularized_lda(
                    unbalanced_models[fold_idx],
                    x_test_fold,
                )
                fold_acc_unbalanced.append(
                    float(DOT_SIDE.accuracy(y_test_fold, y_pred_unbalanced))
                )

                rng = np.random.default_rng(
                    _test_balance_seed(random_seed, train_idx, test_idx, fold_idx)
                )
                test_bal_idx = DOT_SIDE.undersample_majority_to_minority_indices(
                    y_test_fold,
                    rng=rng,
                )
                if test_bal_idx.size == 0:
                    continue
                y_pred_balanced = predict_regularized_lda(
                    balanced_models[fold_idx],
                    x_test_fold[test_bal_idx],
                )
                fold_acc_balanced.append(
                    float(DOT_SIDE.accuracy(y_test_fold[test_bal_idx], y_pred_balanced))
                )

            n_valid_unbalanced_folds[train_idx, test_idx] = len(fold_acc_unbalanced)
            n_valid_balanced_folds[train_idx, test_idx] = len(fold_acc_balanced)
            if len(fold_acc_unbalanced) == n_splits:
                accuracy_unbalanced[train_idx, test_idx] = float(
                    np.mean(fold_acc_unbalanced)
                )
            if len(fold_acc_balanced) == n_splits:
                accuracy_balanced[train_idx, test_idx] = float(
                    np.mean(fold_acc_balanced)
                )
            elif len(fold_acc_balanced) == 0:
                cell_skip_reasons[train_idx, test_idx] = (
                    "balanced_test_not_possible_any_fold"
                )
            else:
                cell_skip_reasons[train_idx, test_idx] = (
                    "balanced_test_not_possible_all_folds"
                )

        row["n_feature_dimension_matched_test_cells"] = int(matched_dim_count)
        row["n_valid_balanced_test_cells"] = int(
            np.isfinite(accuracy_balanced[train_idx]).sum()
        )
        row["n_invalid_balanced_test_cells"] = int(
            n_times - row["n_valid_balanced_test_cells"]
        )
        train_rows.append(row)

    return {
        "times": times,
        "window_rows": window_rows,
        "train_rows": train_rows,
        "test_class_counts": test_class_counts,
        "feature_dims": feature_dims,
        "accuracy_unbalanced": accuracy_unbalanced,
        "accuracy_balanced": accuracy_balanced,
        "n_valid_unbalanced_folds": n_valid_unbalanced_folds,
        "n_valid_balanced_folds": n_valid_balanced_folds,
        "cell_skip_reasons": cell_skip_reasons,
    }


def summarize_time_generalization_diagnostics(
    *,
    decode_result: Mapping[str, Any],
    averaging_metadata: Mapping[str, Any],
    n_decoder_trials: int,
    n_decoder_identities: int,
    n_splits: int,
    window_samples: int,
    sfreq: float,
) -> dict[str, Any]:
    """Summarize train/test validity, class balance, and NaN cell policy."""
    train_rows = list(decode_result["train_rows"])
    accuracy_balanced = np.asarray(decode_result["accuracy_balanced"], dtype=float)
    cell_skip_reasons = np.asarray(decode_result["cell_skip_reasons"], dtype=object)
    train_skip_reasons = Counter(
        row.get("train_skip_reason", "") or "decoded" for row in train_rows
    )
    cell_reason_counter = Counter(
        str(reason) if str(reason) else "decoded"
        for reason in cell_skip_reasons.reshape(-1)
    )

    train_minority_counts: list[int] = []
    for row in train_rows:
        counts = {int(k): int(v) for k, v in dict(row.get("train_class_counts", {})).items()}
        if counts:
            train_minority_counts.append(min(counts.get(0, 0), counts.get(1, 0)))

    n_train_times = int(len(train_rows))
    n_cells = int(accuracy_balanced.size)
    n_valid_cells = int(np.isfinite(accuracy_balanced).sum())
    valid_fraction = float(n_valid_cells / n_cells) if n_cells else 0.0
    n_valid_train_times = int(train_skip_reasons.get("decoded", 0))

    if n_cells == 0 or n_valid_cells == 0:
        assessment_status = "not_well_defined"
        assessment_text = (
            "No train/test cells produced finite balanced decoding accuracy. "
            "Inspect train-time and cell-level skip reasons before interpreting "
            "the matrix."
        )
    elif n_valid_cells == n_cells:
        assessment_status = "well_defined"
        assessment_text = (
            "Every train/test cell was valid under the configured stratified "
            "folds, feature-window policy, and balanced evaluation rule."
        )
    else:
        assessment_status = "partially_defined"
        assessment_text = (
            "Some train/test cells were NaN because train labels, window feature "
            "dimensions, or held-out test labels could not support the configured "
            "time-generalization evaluation."
        )

    return {
        "classifier": {
            "name": DOT_SIDE.DECODER_MODEL_NAME,
            "details": DOT_SIDE.DECODER_MODEL_DETAILS,
            "classes": {"0": "right_or_center", "1": "left"},
        },
        "cross_validation": {
            "scheme": "train_time_stratified_k_fold_reused_across_test_times",
            "n_splits": int(n_splits),
            "split_unit": "same-identity pseudo-trials",
            "fold_reuse_policy": (
                "For each train-time row, stratified folds are generated once "
                "from train-time labels and reused for every test-time column."
            ),
        },
        "feature_window": {
            "window_samples": int(window_samples),
            "nominal_window_ms": float(window_samples) / float(sfreq) * 1000.0,
            "edge_windows_are_clipped": True,
            "cross_time_feature_policy": (
                "Train/test cells require equal flattened feature dimensions; "
                "mismatched edge-clipped window widths are stored as NaN."
            ),
        },
        "evaluation": {
            "primary_metric": BALANCED_METRIC,
            "secondary_metric": UNBALANCED_METRIC,
            "balanced_cell_policy": (
                "A balanced train/test cell is finite only when every reused "
                "held-out fold contains both test-time classes after random "
                "undersampling."
            ),
            "n_finite_balanced_cells": n_valid_cells,
            "n_finite_unbalanced_cells": int(
                np.isfinite(np.asarray(decode_result["accuracy_unbalanced"])).sum()
            ),
        },
        "sample_structure": {
            "averaging_mode": str(averaging_metadata.get("mode", "fixed_count")),
            "requested_value": str(averaging_metadata.get("requested_value", "1")),
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
            "n_train_times": n_train_times,
            "n_valid_train_times": n_valid_train_times,
            "n_skipped_train_times": int(n_train_times - n_valid_train_times),
            "train_skip_reasons": {
                str(key): int(value) for key, value in sorted(train_skip_reasons.items())
            },
            "min_train_time_minority_class_count": (
                int(min(train_minority_counts)) if train_minority_counts else None
            ),
        },
        "cell_validity": {
            "n_cells": n_cells,
            "n_valid_balanced_cells": n_valid_cells,
            "n_invalid_balanced_cells": int(n_cells - n_valid_cells),
            "valid_balanced_cell_fraction": valid_fraction,
            "cell_skip_reasons": {
                str(key): int(value) for key, value in sorted(cell_reason_counter.items())
            },
        },
        "assessment": {
            "status": assessment_status,
            "text": assessment_text,
        },
    }


# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------


def add_occlusion_boundary_lines_2d(
    ax: plt.Axes,
    *,
    condition: str,
    occlusion_timing: Mapping[str, Any] | None,
) -> None:
    """Add vertical test-time and horizontal train-time event markers."""
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
        value = float(occlusion_timing[key])
        ax.axvline(value, color=color, lw=1.2, ls=linestyle, alpha=0.9, label=label)
        ax.axhline(value, color=color, lw=1.2, ls=linestyle, alpha=0.9, label="_nolegend_")


def plot_time_generalization_heatmap(
    *,
    times: np.ndarray,
    matrix: np.ndarray,
    condition: str,
    model_source: str,
    n_splits: int,
    n_decoder_channels: int,
    occlusion_timing: Mapping[str, Any] | None,
    out_file: Path,
) -> Path:
    """Write the participant train-time x test-time heatmap."""
    finite_values = matrix[np.isfinite(matrix)]
    if finite_values.size == 0:
        raise RuntimeError("Cannot plot a time-generalization matrix with no finite cells.")
    vmin = min(0.35, float(np.nanmin(finite_values)))
    vmax = max(0.75, float(np.nanmax(finite_values)))
    if vmin >= 0.5:
        vmin = 0.35
    if vmax <= 0.5:
        vmax = 0.75

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#d9d9d9")
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.5, vmax=vmax)

    fig, ax = plt.subplots(figsize=(9.5, 8), constrained_layout=True)
    im = ax.imshow(
        np.ma.masked_invalid(matrix),
        origin="lower",
        interpolation="nearest",
        aspect="auto",
        extent=[float(times[0]), float(times[-1]), float(times[0]), float(times[-1])],
        cmap=cmap,
        norm=norm,
    )
    ax.plot([times[0], times[-1]], [times[0], times[-1]], color="white", lw=1.2, alpha=0.9)
    ax.axvline(0.0, color="k", lw=1.0, ls="--", alpha=0.7)
    ax.axhline(0.0, color="k", lw=1.0, ls="--", alpha=0.7)
    add_occlusion_boundary_lines_2d(
        ax,
        condition=condition,
        occlusion_timing=occlusion_timing,
    )
    ax.set_xlabel("Test time from condition onset (s)")
    ax.set_ylabel("Train time from condition onset (s)")
    ax.set_title(
        "Dot-side time generalization | "
        f"{condition} | {model_source} model | {n_splits}-fold CV | "
        f"n={n_decoder_channels} gradiometers"
    )
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper right", frameon=True)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(BALANCED_METRIC)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_file


def write_matrix_csv(
    path: Path,
    *,
    times: np.ndarray,
    decode_result: Mapping[str, Any],
) -> Path:
    """Write one long-form CSV row per train/test cell."""
    window_rows = list(decode_result["window_rows"])
    train_rows = list(decode_result["train_rows"])
    test_class_counts = list(decode_result["test_class_counts"])
    acc_bal = np.asarray(decode_result["accuracy_balanced"], dtype=float)
    acc_un = np.asarray(decode_result["accuracy_unbalanced"], dtype=float)
    n_bal_folds = np.asarray(decode_result["n_valid_balanced_folds"], dtype=int)
    n_un_folds = np.asarray(decode_result["n_valid_unbalanced_folds"], dtype=int)
    reasons = np.asarray(decode_result["cell_skip_reasons"], dtype=object)
    feature_dims = np.asarray(decode_result["feature_dims"], dtype=int)

    records: list[dict[str, Any]] = []
    for train_idx, train_row in enumerate(train_rows):
        train_window = window_rows[train_idx]
        for test_idx, test_window in enumerate(window_rows):
            records.append(
                {
                    "train_index": int(train_idx),
                    "test_index": int(test_idx),
                    "train_time_sec": float(times[train_idx]),
                    "test_time_sec": float(times[test_idx]),
                    "train_window_start_time_sec": float(train_window["window_start_time_sec"]),
                    "train_window_end_time_sec": float(train_window["window_end_time_sec"]),
                    "test_window_start_time_sec": float(test_window["window_start_time_sec"]),
                    "test_window_end_time_sec": float(test_window["window_end_time_sec"]),
                    "train_effective_window_samples": int(train_window["effective_window_samples"]),
                    "test_effective_window_samples": int(test_window["effective_window_samples"]),
                    "train_feature_dim": int(feature_dims[train_idx]),
                    "test_feature_dim": int(feature_dims[test_idx]),
                    "train_class_counts": json.dumps(train_row.get("train_class_counts", {}), sort_keys=True),
                    "test_class_counts": json.dumps(test_class_counts[test_idx], sort_keys=True),
                    "train_skip_reason": str(train_row.get("train_skip_reason", "")),
                    "cell_skip_reason": str(reasons[train_idx, test_idx]),
                    "n_valid_unbalanced_folds": int(n_un_folds[train_idx, test_idx]),
                    "n_valid_balanced_folds": int(n_bal_folds[train_idx, test_idx]),
                    UNBALANCED_METRIC: acc_un[train_idx, test_idx],
                    BALANCED_METRIC: acc_bal[train_idx, test_idx],
                }
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame.from_records(records).to_csv(path, index=False)
    return path


def write_matrix_npz(
    path: Path,
    *,
    times: np.ndarray,
    decode_result: Mapping[str, Any],
) -> Path:
    """Write compressed NumPy arrays for reliable matrix reloads."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        times_sec=np.asarray(times, dtype=float),
        train_times_sec=np.asarray(times, dtype=float),
        test_times_sec=np.asarray(times, dtype=float),
        accuracy_balanced_undersampled_sets=np.asarray(
            decode_result["accuracy_balanced"],
            dtype=float,
        ),
        accuracy_unbalanced_train=np.asarray(
            decode_result["accuracy_unbalanced"],
            dtype=float,
        ),
        n_valid_balanced_folds=np.asarray(
            decode_result["n_valid_balanced_folds"],
            dtype=int,
        ),
        n_valid_unbalanced_folds=np.asarray(
            decode_result["n_valid_unbalanced_folds"],
            dtype=int,
        ),
        cell_skip_reasons=np.asarray(
            decode_result["cell_skip_reasons"],
            dtype="U96",
        ),
        train_skip_reasons=np.asarray(
            [
                str(row.get("train_skip_reason", ""))
                for row in decode_result["train_rows"]
            ],
            dtype="U96",
        ),
        effective_window_samples=np.asarray(
            [row["effective_window_samples"] for row in decode_result["window_rows"]],
            dtype=int,
        ),
        feature_dims=np.asarray(decode_result["feature_dims"], dtype=int),
    )
    return path


def _format_optional_number(value: Any, *, digits: int = 3) -> str:
    """Format optional numbers for compact report tables."""
    if value is None:
        return "n/a"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def time_generalization_diagnostics_html(summary: Mapping[str, Any]) -> str:
    """Render participant time-generalization diagnostics as HTML."""
    diagnostics = summary.get("time_generalization_diagnostics", {})
    classifier = diagnostics.get("classifier", {})
    cv = diagnostics.get("cross_validation", {})
    feature_window = diagnostics.get("feature_window", {})
    evaluation = diagnostics.get("evaluation", {})
    sample = diagnostics.get("sample_structure", {})
    class_balance = diagnostics.get("class_balance", {})
    cell_validity = diagnostics.get("cell_validity", {})
    assessment = diagnostics.get("assessment", {})
    cell_reasons = cell_validity.get("cell_skip_reasons", {})
    reason_text = ", ".join(
        f"{html.escape(str(key))}: {int(value)}"
        for key, value in sorted(dict(cell_reasons).items())
    ) or "none"
    valid_fraction = cell_validity.get("valid_balanced_cell_fraction")
    valid_percent = (
        "n/a" if valid_fraction is None else f"{float(valid_fraction) * 100.0:.1f}%"
    )
    rows = [
        ("Classifier", classifier.get("name", "n/a")),
        ("CV scheme", cv.get("scheme", "n/a")),
        ("CV folds", cv.get("n_splits", "n/a")),
        ("Fold reuse", cv.get("fold_reuse_policy", "n/a")),
        ("Averaging mode", sample.get("averaging_mode", "n/a")),
        ("Averaging request", sample.get("requested_value", "n/a")),
        ("Decoder pseudo-trials", sample.get("n_decoder_trials", "n/a")),
        ("Decoder identities", sample.get("n_decoder_identities", "n/a")),
        ("Window samples", feature_window.get("window_samples", "n/a")),
        (
            "Nominal window ms",
            _format_optional_number(feature_window.get("nominal_window_ms"), digits=2),
        ),
        ("Feature policy", feature_window.get("cross_time_feature_policy", "n/a")),
        ("Primary metric", evaluation.get("primary_metric", "n/a")),
        ("Balanced cell policy", evaluation.get("balanced_cell_policy", "n/a")),
        (
            "Valid train times",
            f"{class_balance.get('n_valid_train_times', 'n/a')} / "
            f"{class_balance.get('n_train_times', 'n/a')}",
        ),
        (
            "Valid matrix cells",
            f"{cell_validity.get('n_valid_balanced_cells', 'n/a')} / "
            f"{cell_validity.get('n_cells', 'n/a')} ({valid_percent})",
        ),
        ("Cell skip reasons", reason_text),
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
        "<p>"
        + html.escape(str(assessment.get("text", "")))
        + "</p><table><tbody>"
        + table_rows
        + "</tbody></table>"
    )


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def main() -> None:
    """Run the end-to-end participant time-generalization workflow."""
    args = build_arg_parser().parse_args()
    validate_cli_args(args)
    model_condition = DOT_SIDE.resolve_model_condition(args.condition, args.model_source)
    current_analysis_key = DOT_SIDE.analysis_key(args.condition, args.model_source)
    epochs_fif, model_mat, out_dir, source_manifest, cfg = resolve_input_output_paths(args)
    output_paths = time_generalization_output_paths(
        out_dir,
        args.condition,
        args.model_source,
    )

    stage04_validation = DOT_SIDE.validate_stage04_epoching_contract(
        source_manifest=source_manifest,
        epochs_fif=epochs_fif,
        cfg=cfg,
    )

    if args.resume and all(path.is_file() for path in output_paths.values()):
        print("Time-generalization outputs already exist; resume mode reuses them.")
        for label, path in output_paths.items():
            print(f"{label}: {path}")
        return

    if not args.validate_only and not args.resume and not args.overwrite:
        existing_outputs = [str(path) for path in output_paths.values() if path.exists()]
        if existing_outputs:
            raise RuntimeError(
                "Refusing to overwrite existing time-generalization outputs "
                f"without --overwrite or --resume: {existing_outputs}"
            )

    if args.validate_only:
        epochs_header = mne.read_epochs(epochs_fif, preload=False, verbose=False)
        validation = DOT_SIDE.validate_decoding_input_contract(
            epochs=epochs_header,
            model_mat=model_mat,
            condition=args.condition,
            model_source=args.model_source,
            cfg=cfg,
            n_times_target=int(round((args.decode_tmax - args.decode_tmin) * args.decode_sfreq)) + 1,
        )
        if args.print_summary:
            print(
                json.dumps(
                    {
                        "stage04_validation": stage04_validation,
                        "decoding_input_validation": validation,
                        "out_dir": str(out_dir),
                        "outputs": {
                            key: str(value) for key, value in output_paths.items()
                        },
                        "analysis_key": current_analysis_key,
                        "model_source": args.model_source,
                        "model_condition": model_condition,
                        "decode_window_samples": int(args.decode_window_samples),
                        "trials_per_average": DOT_SIDE.format_trial_averaging_request(
                            args.trials_per_average
                        ),
                        "random_seed": int(args.random_seed),
                    },
                    indent=2,
                )
            )
        else:
            print("Dot-side time-generalization input validation passed.")
        return

    # -------------------------------------------------------------------------
    # Load and align participant data with dynamic model-derived labels.
    # -------------------------------------------------------------------------
    epochs_all = mne.read_epochs(epochs_fif, preload=False, verbose=False)
    DOT_SIDE.validate_decoding_input_contract(
        epochs=epochs_all,
        model_mat=model_mat,
        condition=args.condition,
        model_source=args.model_source,
        cfg=cfg,
        n_times_target=int(round((args.decode_tmax - args.decode_tmin) * args.decode_sfreq)) + 1,
    )

    epochs_condition = epochs_all[args.condition].copy().load_data()
    decoder_channels = DOT_SIDE.select_posterior_decoder_channels(epochs_condition.info)
    epochs_decode = (
        epochs_condition.copy()
        .crop(tmin=args.decode_tmin, tmax=args.decode_tmax)
        .pick(decoder_channels)
        .resample(args.decode_sfreq)
    )
    trial_identities = epochs_decode.metadata["trial_sequence_identity"].to_numpy(dtype=int)
    identity_values = list(
        range(
            int(cfg["subject_profile"]["identity_range"][0]),
            int(cfg["subject_profile"]["identity_range"][1]) + 1,
        )
    )
    labels, _x_by_identity, center_x = DOT_SIDE.build_dot_side_labels(
        model_mat=model_mat,
        condition=model_condition,
        trial_identities=trial_identities,
        identity_values=identity_values,
        n_times_target=len(epochs_decode.times),
    )
    n_times_aligned = min(labels.shape[1], len(epochs_decode.times))
    if n_times_aligned < 1:
        raise RuntimeError("No overlapping decoding/model time samples remain.")
    if n_times_aligned < len(epochs_decode.times):
        epochs_decode = epochs_decode.copy().crop(
            tmin=float(epochs_decode.times[0]),
            tmax=float(epochs_decode.times[n_times_aligned - 1]),
        )
    labels = labels[:, :n_times_aligned]
    times = epochs_decode.times[:n_times_aligned].copy()
    decode_data_source = epochs_decode.get_data()[:, :, :n_times_aligned]
    source_epoch_indices = np.asarray(epochs_decode.selection, dtype=int)

    # -------------------------------------------------------------------------
    # Preserve same-identity averaging, then train/test every time pair.
    # -------------------------------------------------------------------------
    decode_data, decode_labels, decode_identities, averaging_metadata = (
        DOT_SIDE.average_trials_within_identity(
            decode_data_source,
            labels,
            identities=trial_identities,
            trials_per_average=args.trials_per_average,
            random_seed=args.random_seed,
            source_epoch_indices=source_epoch_indices,
        )
    )
    decode_result = decode_time_generalization(
        decode_data,
        decode_labels,
        times=times,
        n_splits=args.n_splits,
        window_samples=args.decode_window_samples,
        random_seed=args.random_seed,
        sfreq=float(epochs_decode.info["sfreq"]),
    )
    accuracy_balanced = np.asarray(decode_result["accuracy_balanced"], dtype=float)
    accuracy_unbalanced = np.asarray(decode_result["accuracy_unbalanced"], dtype=float)
    if not np.any(np.isfinite(accuracy_balanced)):
        raise RuntimeError("No finite balanced time-generalization cells were produced.")

    peak_flat_idx = int(np.nanargmax(accuracy_balanced))
    peak_train_idx, peak_test_idx = np.unravel_index(peak_flat_idx, accuracy_balanced.shape)
    occlusion_timing = (
        DOT_SIDE.occlusion_timing_from_config(cfg)
        if "occluded" in args.condition.lower()
        else None
    )
    diagnostics = summarize_time_generalization_diagnostics(
        decode_result=decode_result,
        averaging_metadata=averaging_metadata,
        n_decoder_trials=int(decode_data.shape[0]),
        n_decoder_identities=int(np.unique(decode_identities).size),
        n_splits=int(args.n_splits),
        window_samples=int(args.decode_window_samples),
        sfreq=float(epochs_decode.info["sfreq"]),
    )

    # -------------------------------------------------------------------------
    # Export matrix bundles, heatmap, and machine-readable summary.
    # -------------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_npz = write_matrix_npz(
        output_paths["matrix_npz"],
        times=times,
        decode_result=decode_result,
    )
    matrix_csv = write_matrix_csv(
        output_paths["matrix_csv"],
        times=times,
        decode_result=decode_result,
    )
    fig_path = plot_time_generalization_heatmap(
        times=times,
        matrix=accuracy_balanced,
        condition=args.condition,
        model_source=args.model_source,
        n_splits=int(args.n_splits),
        n_decoder_channels=len(decoder_channels),
        occlusion_timing=occlusion_timing,
        out_file=output_paths["figure"],
    )

    summary = {
        "stage": MATRIX_STAGE,
        "condition": args.condition,
        "model_source": args.model_source,
        "model_condition": model_condition,
        "analysis_key": current_analysis_key,
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
            "Centered sliding windows on the resampled decoding grid; epoch "
            "edges are clipped and train/test cells with mismatched flattened "
            "feature dimensions are stored as NaN."
        ),
        "n_splits": int(args.n_splits),
        "random_seed": int(args.random_seed),
        "n_decoder_channels": int(len(decoder_channels)),
        "decoder_channels": decoder_channels,
        "occlusion_timing": occlusion_timing,
        "trial_averaging": averaging_metadata,
        "time_generalization_diagnostics": diagnostics,
        "training_testing_diagnostics": diagnostics,
        "center_x": float(center_x),
        "primary_metric": BALANCED_METRIC,
        "secondary_metric": UNBALANCED_METRIC,
        "peak_accuracy_balanced_train_test": float(accuracy_balanced[peak_train_idx, peak_test_idx]),
        "peak_train_time_sec_balanced": float(times[peak_train_idx]),
        "peak_test_time_sec_balanced": float(times[peak_test_idx]),
        "mean_accuracy_balanced_train_test": float(np.nanmean(accuracy_balanced)),
        "mean_accuracy_unbalanced_train_test": float(np.nanmean(accuracy_unbalanced)),
        "matrix_npz": str(matrix_npz),
        "matrix_csv": str(matrix_csv),
        "figure_path": str(fig_path),
        "report_html": str(output_paths["report_html"]),
        "epochs_fif": str(epochs_fif),
        "model_mat": str(model_mat),
    }
    summary_json = output_paths["summary_json"]
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # -------------------------------------------------------------------------
    # Export MNE HTML report with the heatmap and diagnostics.
    # -------------------------------------------------------------------------
    report = mne.Report(
        title=(
            f"Dot-side time generalization | sub-{args.subject} ses-{args.session} "
            f"| {args.condition} | {args.model_source}"
        )
    )
    report.add_html(
        "<pre>"
        + html.escape(
            "Dot-side time-generalization summary\n"
            f"condition={summary['condition']}\n"
            f"model_source={summary['model_source']}\n"
            f"model_condition={summary['model_condition']}\n"
            f"analysis_key={summary['analysis_key']}\n"
            f"subject=sub-{summary['subject']}\n"
            f"session=ses-{summary['session']}\n"
            f"n_source_trials={summary['n_source_trials']}\n"
            f"n_decoder_trials={summary['n_decoder_trials']}\n"
            f"trial_averaging_mode={summary['trial_averaging']['mode']}\n"
            f"trials_per_average={summary['trial_averaging']['requested_value']}\n"
            f"n_times={summary['n_times']} @ {summary['sfreq_hz']} Hz\n"
            f"decode_window_samples={summary['decode_window_samples']}\n"
            f"n_decoder_channels={summary['n_decoder_channels']}\n"
            f"peak_balanced_accuracy={summary['peak_accuracy_balanced_train_test']:.4f} "
            f"at train={summary['peak_train_time_sec_balanced']:.3f} s, "
            f"test={summary['peak_test_time_sec_balanced']:.3f} s\n"
            f"mean_balanced_accuracy={summary['mean_accuracy_balanced_train_test']:.4f}\n"
            f"matrix_npz={summary['matrix_npz']}\n"
            f"matrix_csv={summary['matrix_csv']}\n"
            f"figure_path={summary['figure_path']}"
        )
        + "</pre>",
        title="Time-generalization summary",
        section="sanity_check",
        tags=("summary", "time_generalization", "sanity_check"),
    )
    fig_for_report = plt.figure(figsize=(10, 8))
    img = plt.imread(fig_path)
    plt.imshow(img)
    plt.axis("off")
    report.add_figure(
        fig_for_report,
        title="Train-time x test-time balanced accuracy",
        section="sanity_check",
        tags=("figure", "time_generalization", "sanity_check"),
    )
    plt.close(fig_for_report)
    report.add_html(
        time_generalization_diagnostics_html(summary),
        title="Training/testing diagnostics",
        section="sanity_check",
        tags=("summary", "time_generalization", "training_testing", "sanity_check"),
    )
    report.save(output_paths["report_html"], overwrite=True, open_browser=False)

    print("Dot-side time-generalization complete.")
    print(f"Summary JSON: {summary_json}")
    print(f"Matrix NPZ: {matrix_npz}")
    print(f"Matrix CSV: {matrix_csv}")
    print(f"Figure: {fig_path}")
    print(f"Report: {output_paths['report_html']}")


if __name__ == "__main__":
    main()
