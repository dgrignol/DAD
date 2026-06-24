"""
test_dot_side_time_generalization.py

Purpose
Validate the dot-side time-generalization helper contracts with synthetic
arrays and temporary subject outputs, without loading participant MEG data.

Usage
Run only these tests:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_dot_side_time_generalization.py -q
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SUBJECT_SCRIPT = (
    REPO_ROOT / "meg" / "code" / "sanity_checks" / "03_dot_side_time_generalization.py"
)
GROUP_SCRIPT = (
    REPO_ROOT
    / "meg"
    / "code"
    / "sanity_checks"
    / "04_dot_side_time_generalization_group_report.py"
)


def _load_script(name: str, path: Path):
    """Load a script module from a numeric filename."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


time_gen = _load_script("dot_side_time_generalization", SUBJECT_SCRIPT)
tg_group = _load_script("dot_side_time_generalization_group", GROUP_SCRIPT)


def test_time_generalization_marks_feature_dimension_mismatch():
    """Edge-clipped train/test windows with different widths must be NaN cells."""
    rng = np.random.default_rng(1)
    data = rng.normal(size=(8, 2, 5))
    labels = np.tile(np.array([[0], [0], [0], [0], [1], [1], [1], [1]]), (1, 5))
    result = time_gen.decode_time_generalization(
        data,
        labels,
        times=np.arange(5, dtype=float) / 10.0,
        n_splits=2,
        window_samples=3,
        random_seed=11,
        sfreq=10.0,
    )

    matrix = result["accuracy_balanced"]
    reasons = result["cell_skip_reasons"]
    assert np.isfinite(matrix[2, 2])
    assert np.isnan(matrix[2, 0])
    assert reasons[2, 0] == "feature_dimension_mismatch"
    assert result["n_valid_balanced_folds"][2, 2] == 2


def test_time_generalization_skips_invalid_train_rows_and_test_columns():
    """Train-time and test-time balance failures should be explicit NaNs."""
    rng = np.random.default_rng(2)
    data = rng.normal(size=(8, 2, 5))
    labels = np.tile(np.array([[0], [0], [0], [0], [1], [1], [1], [1]]), (1, 5))
    labels[:, 1] = 0
    labels[0, 1] = 1
    labels[:, 3] = 0

    result = time_gen.decode_time_generalization(
        data,
        labels,
        times=np.arange(5, dtype=float) / 10.0,
        n_splits=2,
        window_samples=1,
        random_seed=11,
        sfreq=10.0,
    )

    reasons = result["cell_skip_reasons"]
    assert set(reasons[1]) == {"too_few_train_class_samples_for_2_fold_stratified_cv"}
    assert np.all(np.isnan(result["accuracy_balanced"][1]))
    assert reasons[0, 3] == "balanced_test_not_possible_any_fold"
    assert np.isnan(result["accuracy_balanced"][0, 3])


def test_2d_cluster_permutation_recovers_positive_block():
    """A contiguous above-chance train/test block should form a 2D cluster."""
    values = np.full((7, 6, 6), 0.5, dtype=float)
    values[:, 2:4, 3:5] += np.array([0.12, 0.14, 0.13, 0.16, 0.15, 0.11, 0.17])[
        :, None, None
    ]
    result = tg_group.cluster_permutation_2d(
        values,
        chance=0.5,
        cluster_forming_p=0.05,
        cluster_alpha=0.05,
        tail="greater",
        n_permutations=0,
        random_seed=1,
    )

    assert result["permutation_mode"] == "exact"
    significant = [cluster for cluster in result["clusters"] if cluster["significant"]]
    assert significant
    cluster = significant[0]
    assert cluster["train_start_index"] == 2
    assert cluster["train_end_index"] == 3
    assert cluster["test_start_index"] == 3
    assert cluster["test_end_index"] == 4


def test_2d_cluster_permutation_nan_gap_prevents_bridge():
    """Invalid group cells should split clusters instead of bridging them."""
    values = np.full((4, 5, 7), 0.5, dtype=float)
    values[:, 2, 1:3] = 0.66
    values[:, 2, 4:6] = 0.67
    values[:, 2, 3] = 0.68
    values[0, 2, 3] = np.nan

    result = tg_group.cluster_permutation_2d(
        values,
        chance=0.5,
        cluster_forming_p=0.05,
        cluster_alpha=0.05,
        tail="greater",
        n_permutations=0,
        random_seed=1,
    )

    assert result["n_invalid_cells"] == 1
    assert np.isnan(result["observed_t"][2, 3])
    bounds = [
        (
            cluster["train_start_index"],
            cluster["train_end_index"],
            cluster["test_start_index"],
            cluster["test_end_index"],
        )
        for cluster in result["clusters"]
    ]
    assert (2, 2, 1, 2) in bounds
    assert (2, 2, 4, 5) in bounds


def _synthetic_tg_diagnostics(n_cells: int) -> dict[str, object]:
    """Return compact diagnostics for synthetic subject summaries."""
    return {
        "classifier": {
            "name": "Regularized linear discriminant analysis",
            "details": "Train-fold z-scored binary LDA with ridge regularization.",
        },
        "cross_validation": {
            "scheme": "train_time_stratified_k_fold_reused_across_test_times",
            "n_splits": 2,
            "split_unit": "same-identity pseudo-trials",
        },
        "feature_window": {
            "window_samples": 1,
            "nominal_window_ms": 8.33,
            "cross_time_feature_policy": "matching feature dimensions",
        },
        "evaluation": {
            "primary_metric": "accuracy_balanced_undersampled_sets",
            "balanced_cell_policy": "all folds must be balanceable",
        },
        "sample_structure": {
            "averaging_mode": "fixed_count",
            "requested_value": "1",
            "n_source_trials": 8,
            "n_decoder_trials": 8,
            "n_decoder_identities": 4,
        },
        "class_balance": {
            "n_train_times": 4,
            "n_valid_train_times": 4,
            "n_skipped_train_times": 0,
        },
        "cell_validity": {
            "n_cells": n_cells,
            "n_valid_balanced_cells": n_cells,
            "valid_balanced_cell_fraction": 1.0,
            "cell_skip_reasons": {"decoded": n_cells},
        },
        "assessment": {
            "status": "well_defined",
            "text": "Synthetic matrix is fully valid.",
        },
    }


def _write_subject_output(
    root: Path,
    subject: str,
    times: np.ndarray,
    matrix: np.ndarray,
) -> None:
    """Write a minimal time-generalization output bundle for loader tests."""
    out_dir = (
        root
        / "sanity_checks"
        / subject
        / "ses-01"
        / "dot_side_time_generalization"
        / "always_visible"
        / "observed"
    )
    out_dir.mkdir(parents=True)
    stem = "always_visible_observed"
    matrix_npz = out_dir / f"{stem}_dot_side_time_generalization_matrix.npz"
    figure = out_dir / f"{stem}_dot_side_time_generalization_heatmap.png"
    summary_json = out_dir / f"{stem}_dot_side_time_generalization_summary.json"
    np.savez_compressed(
        matrix_npz,
        times_sec=times,
        accuracy_balanced_undersampled_sets=matrix,
        accuracy_unbalanced_train=matrix + 0.01,
        cell_skip_reasons=np.full(matrix.shape, "", dtype="U8"),
    )
    figure.write_bytes(b"\x89PNG\r\n\x1a\n")
    summary_json.write_text(
        json.dumps(
            {
                "stage": "dot_side_time_generalization",
                "condition": "always_visible",
                "model_source": "observed",
                "model_condition": "always_visible",
                "analysis_key": stem,
                "subject": subject.removeprefix("sub-"),
                "session": "01",
                "n_times": int(times.size),
                "sfreq_hz": 120.0,
                "decode_window_samples": 1,
                "decode_tmin_sec": float(times[0]),
                "decode_tmax_sec": float(times[-1]),
                "n_splits": 2,
                "random_seed": 42,
                "trial_averaging": {
                    "mode": "fixed_count",
                    "requested_value": "1",
                    "trials_per_average": 1,
                },
                "time_generalization_diagnostics": _synthetic_tg_diagnostics(matrix.size),
                "matrix_npz": str(matrix_npz),
                "figure_path": str(figure),
                "peak_accuracy_balanced_train_test": float(np.nanmax(matrix)),
                "peak_train_time_sec_balanced": float(times[0]),
                "peak_test_time_sec_balanced": float(times[0]),
            }
        ),
        encoding="utf-8",
    )


def test_group_report_builds_matrix_report_and_group_first_html(tmp_path):
    """Synthetic participant matrices should produce a group-first report."""
    times = np.arange(4, dtype=float) / 120.0
    for subject_idx, subject in enumerate(("sub-01", "sub-02")):
        matrix = np.full((times.size, times.size), 0.56 + subject_idx * 0.02)
        matrix[1:3, 1:3] += 0.04
        _write_subject_output(tmp_path, subject, times, matrix)

    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    cfg = {
        "resolved_paths": {
            "derivatives_root": str(tmp_path),
            "data_root": str(raw_root),
        },
        "layout": {"session_bids_id": "ses-01"},
    }
    args = tg_group.build_arg_parser().parse_args(
        [
            "--analyses",
            "always_visible_observed",
            "--group-subjects",
            "01",
            "02",
            "--n-permutations",
            "8",
            "--overwrite",
        ]
    )
    report = tg_group.build_group_report(cfg, args)

    assert report["status"] == "completed"
    analysis = report["analyses"]["always_visible_observed"]
    assert Path(analysis["group_figure"]).is_file()
    assert Path(analysis["group_matrix_npz"]).is_file()
    embedded = Path(report["report_files"]["embedded_html"])
    html_text = embedded.read_text(encoding="utf-8")
    assert html_text.index('id="group-analyses"') < html_text.index(
        'id="participants"'
    )
    assert "2D clusters" in html_text
    assert "Participant-specific matrices" in html_text
