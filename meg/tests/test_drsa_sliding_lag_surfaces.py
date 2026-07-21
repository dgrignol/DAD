"""
test_drsa_sliding_lag_surfaces.py

Purpose
Validate the optional Stage-09 dRSA sliding lag-surface calculations and CLI
without running MATLAB. The tests use synthetic matrices with known lag support
and a tiny fake derivative tree with scipy-written MAT files.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_drsa_sliding_lag_surfaces.py -q
"""

from __future__ import annotations

import importlib.util
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat


REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = REPO_ROOT / "meg" / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))


def _load_stage09_module():
    """Load the numbered Stage-09 module for direct function tests."""
    module_path = CODE_DIR / "09_plot_drsa_sliding_lag_surfaces.py"
    spec = importlib.util.spec_from_file_location("stage09_drsa_sliding_lag", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STAGE09 = _load_stage09_module()


def _make_record(
    *,
    condition: str,
    label: str,
    model_source: str = "observed",
):
    """Create a minimal MatrixRecord for direct target-selection tests."""
    matrix = np.zeros((5, 5), dtype=float)
    axis = np.arange(5, dtype=float) / 10.0
    return STAGE09.MatrixRecord(
        owner_scope="subject",
        owner_label="sub-01",
        condition=condition,
        analysis_family="A_dRSA_alt",
        sensor_group="occipital",
        drsa_type="corr",
        model_source=model_source,
        stimulus_model_label=label,
        model_time_axis_label=f"model time (s) - {label}",
        variable="dRSA_neural_position",
        slice_index=0,
        slice_label="matrix",
        matrix=matrix,
        time_axis_sec=axis,
        relative_result=f"fake/{condition}/{model_source}.mat",
        source_paths=("/tmp/fake.mat",),
        subjects=("sub-01",),
    )


def _make_surface(values):
    """Create a minimal SlidingLagSurface around an explicit value array."""
    values = np.asarray(values, dtype=float)
    lags = np.array([-1, 0, 1], dtype=int)
    lag_sec = lags.astype(float) / 10.0
    coords = np.arange(values.shape[0], dtype=float) / 10.0
    record = _make_record(condition="always_visible", label="actual visible position")
    return STAGE09.SlidingLagSurface(
        record=record,
        target=STAGE09.TARGET_SUMMARIES[0],
        lag_samples=lags,
        lag_sec=lag_sec,
        window_starts_sec=coords,
        window_ends_sec=coords + 0.1,
        window_coords_sec=coords,
        values=values,
        valid_counts=np.where(np.isfinite(values), 1, 0),
        model_window_samples=2,
        model_step_samples=1,
    )


def test_sliding_lag_surface_uses_model_anchored_neural_support():
    """A window restricts model rows, while neural columns may leave the window."""
    matrix = np.full((5, 5), np.nan)
    matrix[0, 0] = 1.0
    matrix[1, 1] = 3.0
    matrix[0, 4] = 8.0
    axis = np.arange(5, dtype=float) / 10.0

    (
        lag_samples,
        _lag_sec,
        window_coords,
        values,
        valid_counts,
        window_samples,
        step_samples,
    ) = STAGE09.compute_sliding_lag_surface(
        matrix,
        axis,
        lag_window_sec=0.4,
        model_window_sec=0.2,
        model_step_sec=0.1,
        min_valid_samples=1,
    )

    assert window_samples == 2
    assert step_samples == 1
    np.testing.assert_allclose(window_coords[0], 0.05)
    zero_lag_index = int(np.where(lag_samples == 0)[0][0])
    plus_four_index = int(np.where(lag_samples == 4)[0][0])
    assert values[0, zero_lag_index] == 2.0
    assert valid_counts[0, zero_lag_index] == 2
    assert values[0, plus_four_index] == 8.0
    assert valid_counts[0, plus_four_index] == 1


def test_peak_timecourse_extracts_argmax_lag_and_relative_threshold():
    """Peak timecourses summarize each window's strongest averaged lag value."""
    surface = _make_surface(
        [
            [1.0, 4.0, 2.0],
            [np.nan, np.nan, np.nan],
            [2.0, 3.0, 5.0],
        ]
    )

    peak_magnitudes, peak_lags = STAGE09.extract_peak_timecourse(surface)

    np.testing.assert_allclose(peak_magnitudes[[0, 2]], [4.0, 5.0])
    assert np.isnan(peak_magnitudes[1])
    np.testing.assert_allclose(peak_lags[[0, 2]], [0.0, 0.1])
    assert np.isnan(peak_lags[1])

    args = argparse.Namespace(
        peak_latency_mask_mode="relative",
        peak_latency_min_magnitude_frac=0.2,
        peak_latency_min_magnitude=None,
    )
    assert STAGE09._peak_latency_threshold([surface], args=args) == 1.0


def test_target_selection_returns_the_requested_summaries():
    """The Stage-09 scope includes core position targets and explicit PE variants."""
    records = [
        _make_record(condition="always_visible", label="actual visible position"),
        _make_record(condition="occluded_nondeviant", label="actual nondeviant position"),
        _make_record(condition="occluded_deviant", label="actual deviant position"),
        _make_record(
            condition="occluded_deviant",
            label="expected/predicted position",
            model_source="predicted",
        ),
        _make_record(
            condition="occluded_deviant",
            label="prediction error position",
            model_source="prediction_error",
        ),
        _make_record(
            condition="occluded_deviant",
            label="50/50 unconstrained-prior PE position",
            model_source="pe_prior_unconstrained_50",
        ),
        _make_record(
            condition="occluded_deviant",
            label="50/50 constrained-prior PE position",
            model_source="pe_prior_constrained_50",
        ),
        _make_record(
            condition="occluded_deviant",
            label="precision-weighted PE position",
            model_source="pe_precision_weighted",
        ),
        _make_record(
            condition="occluded_deviant",
            label="Kalman innovation PE position",
            model_source="pe_kalman_innovation",
        ),
        _make_record(
            condition="occluded_deviant",
            label="latent curvature PE",
            model_source="pe_latent_curvature",
        ),
        _make_record(
            condition="occluded_deviant",
            label="belief-update PE",
            model_source="pe_belief_update",
        ),
        _make_record(condition="always_visible", label="actual visible direction"),
    ]

    selected, skipped = STAGE09._target_records(records)

    assert [spec.slug for spec, _record in selected] == [
        "actual-visible-position",
        "actual-nondeviant-position",
        "actual-deviant-position",
        "expected-predicted-position",
        "prediction-error-position",
        "pe-prior-unconstrained-50-position",
        "pe-prior-constrained-50-position",
        "pe-precision-weighted-position",
        "pe-kalman-innovation-position",
        "pe-latent-curvature",
        "pe-belief-update",
    ]
    assert not [item for item in skipped if item.get("reason") == "target matrix not found"]


def _write_fake_site_config(tmp_path: Path) -> Path:
    """Write a site override that points the pipeline at a fake derivative tree."""
    data_root = tmp_path / "data"
    derivatives_root = tmp_path / "derivatives"
    for folder in (
        data_root / "sub01_150426",
        data_root / "sub02_130526",
    ):
        (folder / "input_files").mkdir(parents=True, exist_ok=True)
    site_config = tmp_path / "site.yaml"
    site_config.write_text(
        "\n".join(
            [
                "project:",
                f"  root: \"{REPO_ROOT / 'meg'}\"",
                f"  data_root: \"{data_root}\"",
                f"  derivatives_root: \"{derivatives_root}\"",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return site_config


def _write_fake_result_mat(
    tmp_path: Path,
    subject: str,
    *,
    condition: str,
    model_source: str | None,
    value_offset: float,
) -> None:
    """Write one tiny Stage-07-like occipital A_dRSA_alt MAT file."""
    parts = [
        tmp_path,
        "derivatives",
        "analysis",
        subject,
        "ses-01",
        "drsa",
        "matlab_real_neural",
        "A_dRSA_alt",
        condition,
    ]
    if model_source is not None:
        parts.append(model_source)
    result_path = (
        Path(*parts)
        / "corr"
        / "results"
        / f"A_dRSA_alt_{condition}_{model_source + '_' if model_source else ''}corr_occipital.mat"
    )
    result_path.parent.mkdir(parents=True, exist_ok=True)
    matrix = np.zeros((8, 8, 1), dtype=float)
    for model_idx in range(8):
        for neural_idx in range(8):
            matrix[model_idx, neural_idx, 0] = value_offset + model_idx + neural_idx / 10.0
    savemat(
        result_path,
        {
            "dRSA_neural_position": matrix,
            "time_axis_sec": np.arange(8, dtype=float) / 10.0,
        },
    )


def _write_fake_target_set(tmp_path: Path, subject: str, value_offset: float) -> None:
    """Write all Stage-09 target matrices for one fake subject."""
    targets = [
        ("always_visible", None),
        ("occluded_nondeviant", None),
        ("occluded_deviant", "observed"),
        ("occluded_deviant", "predicted"),
        ("occluded_deviant", "prediction_error"),
        ("occluded_deviant", "pe_prior_unconstrained_50"),
        ("occluded_deviant", "pe_prior_constrained_50"),
        ("occluded_deviant", "pe_precision_weighted"),
        ("occluded_deviant", "pe_kalman_innovation"),
        ("occluded_deviant", "pe_latent_curvature"),
        ("occluded_deviant", "pe_belief_update"),
    ]
    for index, (condition, model_source) in enumerate(targets):
        _write_fake_result_mat(
            tmp_path,
            subject,
            condition=condition,
            model_source=model_source,
            value_offset=value_offset + float(index),
        )


def test_cli_builds_group_report_with_all_three_plot_types(tmp_path):
    """The CLI builds heatmaps, waterfalls, peak timecourses, and a report."""
    site_config = _write_fake_site_config(tmp_path)
    _write_fake_target_set(tmp_path, "sub-01", 10.0)
    _write_fake_target_set(tmp_path, "sub-02", 20.0)

    completed = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "09_plot_drsa_sliding_lag_surfaces.py"),
            "--config",
            str(REPO_ROOT / "meg" / "params" / "pipeline_default.yaml"),
            "--site-config",
            str(site_config),
            "--subjects",
            "01",
            "02",
            "--group-only",
            "--lag-window-sec",
            "0.2",
            "--model-window-sec",
            "0.2",
            "--model-step-sec",
            "0.1",
            "--min-valid-samples",
            "1",
            "--overwrite",
            "--print-summary",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    report = (
        tmp_path
        / "derivatives"
        / "reports"
        / "group"
        / "ses-01"
        / "09_drsa_sliding_lag_surfaces_report.html"
    )
    embedded = report.with_name("09_drsa_sliding_lag_surfaces_report_embedded.html")
    assert report.is_file()
    assert embedded.is_file()
    report_text = report.read_text(encoding="utf-8")
    embedded_text = embedded.read_text(encoding="utf-8")
    target_count = len(STAGE09.TARGET_SUMMARIES)
    comparison_count = len(STAGE09.PEAK_TIMECOURSE_COMPARISONS)
    assert f"Target summaries found: {target_count} of {target_count}" in report_text
    assert "Lag convention: lag = neural time - model time" in report_text
    assert "Peak-latency mask mode: relative" in report_text
    assert "sliding_window_peak_timecourse" in report_text
    assert "Overlapping windows share model-time rows" in report_text
    assert report_text.count("Plotly.newPlot") == target_count
    assert embedded_text.count("data:image/png;base64,") == target_count + comparison_count

    surface_dir = (
        tmp_path
        / "derivatives"
        / "analysis"
        / "group"
        / "ses-01"
        / "drsa"
        / "sliding_lag_surfaces"
    )
    assert len(sorted((surface_dir / "heatmaps").glob("*.png"))) == target_count
    assert len(sorted((surface_dir / "plotly_waterfalls").glob("*.html"))) == target_count
    assert len(sorted((surface_dir / "peak_timecourses").glob("*.png"))) == comparison_count
