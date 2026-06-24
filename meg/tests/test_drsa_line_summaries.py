"""
test_drsa_line_summaries.py

Purpose
Validate the optional Stage-08 dRSA line-summary calculations and CLI without
running MATLAB. The tests use synthetic matrices where the expected lag
behavior is known exactly, plus a tiny fake derivative tree with scipy-written
MAT files.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_drsa_line_summaries.py -q
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat


REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = REPO_ROOT / "meg" / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))

from drsa_result_utils import (  # noqa: E402
    frame_to_seconds,
    iter_matrix_slices,
    model_feature_label,
    model_source_for_matrix_slice,
    occlusion_timing_from_config,
    read_drsa_result_file,
)


def _load_stage08_module():
    """Load the numbered Stage-08 module for direct function tests."""
    module_path = CODE_DIR / "08_plot_drsa_line_summaries.py"
    spec = importlib.util.spec_from_file_location("stage08_drsa_line_summaries", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STAGE08 = _load_stage08_module()


def _make_record(
    matrix,
    *,
    owner_scope: str = "subject",
    owner_label: str = "sub-01",
    subject_records=(),
):
    """Create a minimal MatrixRecord for direct Stage-08 helper tests."""
    axis = np.arange(np.asarray(matrix).shape[0], dtype=float) / 10.0
    return STAGE08.MatrixRecord(
        owner_scope=owner_scope,
        owner_label=owner_label,
        condition="always_visible",
        analysis_family="A_dRSA_alt",
        sensor_group="occipital",
        drsa_type="corr",
        model_source="observed",
        stimulus_model_label="actual visible position",
        model_time_axis_label="model time (s) - actual visible position",
        variable="dRSA_neural_position",
        slice_index=0,
        slice_label="matrix",
        variant="raw",
        matrix=np.asarray(matrix, dtype=float),
        time_axis_sec=axis,
        relative_result="fake/results.mat",
        source_paths=("/tmp/fake/results.mat",),
        subjects=(owner_label,),
        n_subjects=1 if owner_scope == "subject" else len(subject_records),
        subject_records=tuple(subject_records),
    )


def _old_compute_lag_correlation(old_orientation_matrix, lag_window):
    """
    Reproduce the old helper's lag convention.

    The old matrix orientation was neural x model and used matrix.T before
    diagonal averaging. Stage 08 consumes current Stage-07 plotted matrices as
    model x neural, which is the transpose of this synthetic old matrix.
    """
    lags = np.arange(-lag_window, lag_window + 1)
    values = []
    for lag in lags:
        diag_vals = np.diagonal(old_orientation_matrix.T, offset=lag)
        values.append(np.nanmean(diag_vals))
    return lags, np.asarray(values)


def test_full_matrix_lag_curve_matches_old_compute_lag_correlation_convention():
    """Stage-08 current-orientation curves preserve the old lag definition."""
    old_orientation = np.arange(25, dtype=float).reshape(5, 5)
    current_stage07_orientation = old_orientation.T
    lags, expected = _old_compute_lag_correlation(old_orientation, 2)

    actual = STAGE08.compute_interval_lag_curve(
        current_stage07_orientation,
        lags,
        np.ones(5, dtype=bool),
    )

    np.testing.assert_allclose(actual, expected)


def test_interval_restricted_lag_curve_uses_only_model_time_rows_inside_interval():
    """Interval masks restrict diagonal samples by model-time row."""
    matrix = np.full((5, 5), np.nan)
    matrix[0, 0] = 100.0
    matrix[1, 1] = 1.0
    matrix[2, 2] = 3.0
    matrix[3, 3] = 200.0
    interval_mask = np.array([False, True, True, False, False])

    curve = STAGE08.compute_interval_lag_curve(matrix, [0], interval_mask)

    assert curve.tolist() == [2.0]


def test_same_interval_square_lag_curve_requires_neural_time_inside_interval():
    """Strict-square interval mode rejects diagonal samples outside neural interval."""
    matrix = np.full((5, 5), np.nan)
    matrix[1, 2] = 5.0
    matrix[2, 3] = 99.0
    matrix[1, 3] = 10.0
    matrix[2, 4] = 20.0
    interval_mask = np.array([False, True, True, False, False])

    model_anchored = STAGE08.compute_interval_lag_curve(
        matrix,
        [1, 2],
        interval_mask,
    )
    same_square = STAGE08.compute_interval_lag_curve(
        matrix,
        [1, 2],
        interval_mask,
        interval_lag_mode="same-interval-square",
    )

    np.testing.assert_allclose(model_anchored, [52.0, 15.0])
    assert same_square[0] == 5.0
    assert np.isnan(same_square[1])


def test_interval_lag_mode_default_is_model_anchored():
    """The CLI default preserves the previous model-anchored interval behavior."""
    args = STAGE08.build_arg_parser().parse_args(["--config", "fake.yaml"])

    assert args.interval_lag_mode == "model-anchored"
    assert not hasattr(args, "model_time_window_samples")


def test_model_time_window_cli_option_is_removed():
    """The removed peak-lag model-time window option is no longer accepted."""
    with pytest.raises(SystemExit):
        STAGE08.build_arg_parser().parse_args(
            ["--config", "fake.yaml", "--model-time-window-samples", "1"]
        )


def test_none_shading_preserves_existing_interval_curve_behavior():
    """Default none mode returns the existing curve and no band."""
    matrix = np.full((5, 5), np.nan)
    matrix[1, 1] = 1.0
    matrix[2, 2] = 3.0
    interval_mask = np.array([False, True, True, False, False])
    record = _make_record(matrix)

    expected = STAGE08.compute_interval_lag_curve(matrix, [0], interval_mask)
    curve, band, note = STAGE08._interval_lag_summary_for_record(
        record,
        lag_samples=[0],
        interval_mask=interval_mask,
        interval_name="pre_occlusion",
        occlusion_timing={"first_complete_occlusion_sec": 0.3, "reappearance_sec": 0.4},
        shading_mode="none",
        interval_lag_mode="model-anchored",
    )

    np.testing.assert_allclose(curve, expected)
    assert band is None
    assert note is None


def test_within_subject_diagonal_sd_matches_known_diagonal_samples():
    """Within-subject diagonal SD uses the same samples as the lag mean."""
    matrix = np.full((5, 5), np.nan)
    matrix[1, 1] = 1.0
    matrix[2, 2] = 3.0
    interval_mask = np.array([False, True, True, False, False])

    curve, diagonal_sd = STAGE08.compute_interval_lag_curve_with_diagonal_sd(
        matrix,
        [0],
        interval_mask,
    )

    np.testing.assert_allclose(curve, [2.0])
    np.testing.assert_allclose(diagonal_sd, [1.0])


def test_group_within_subject_shading_averages_subject_diagonal_sds():
    """Group within-subject shading averages subject-level diagonal SD bands."""
    first_matrix = np.full((5, 5), np.nan)
    first_matrix[1, 1] = 1.0
    first_matrix[2, 2] = 3.0
    second_matrix = np.full((5, 5), np.nan)
    second_matrix[1, 1] = 10.0
    second_matrix[2, 2] = 14.0
    first = _make_record(first_matrix, owner_label="sub-01")
    second = _make_record(second_matrix, owner_label="sub-02")
    group = _make_record(
        np.zeros((5, 5), dtype=float),
        owner_scope="group",
        owner_label="group",
        subject_records=(first, second),
    )
    timing = {"first_complete_occlusion_sec": 0.3, "reappearance_sec": 0.4}
    interval_mask = np.array([True, True, True, False, False])

    curve, band, note = STAGE08._interval_lag_summary_for_record(
        group,
        lag_samples=[0],
        interval_mask=interval_mask,
        interval_name="pre_occlusion",
        occlusion_timing=timing,
        shading_mode="within-subject-diagonal-sd",
        interval_lag_mode="model-anchored",
    )

    np.testing.assert_allclose(curve, [7.0])
    np.testing.assert_allclose(band, [1.5])
    assert note is None


def test_across_subject_sem_matches_known_subject_curves():
    """Across-subject SEM is computed over subject-level lag curves."""
    curves = np.array(
        [
            [1.0, 2.0],
            [3.0, 6.0],
            [5.0, 10.0],
        ]
    )

    sem = STAGE08.compute_across_subject_sem(curves)

    np.testing.assert_allclose(sem, [2.0 / np.sqrt(3.0), 4.0 / np.sqrt(3.0)])


def test_nan_values_are_handled_safely():
    """NaNs are ignored when possible and preserved when no sample is valid."""
    matrix = np.full((4, 4), np.nan)
    matrix[1, 2] = 5.0
    interval_mask = np.array([False, True, False, False])

    curve = STAGE08.compute_interval_lag_curve(matrix, [0, 1, 2], interval_mask)

    assert np.isnan(curve[0])
    assert curve[1] == 5.0
    assert np.isnan(curve[2])

    curve, diagonal_sd = STAGE08.compute_interval_lag_curve_with_diagonal_sd(
        matrix,
        [1, 2],
        interval_mask,
    )
    assert curve[0] == 5.0
    assert diagonal_sd[0] == 0.0
    assert np.isnan(curve[1])
    assert np.isnan(diagonal_sd[1])

    sem = STAGE08.compute_across_subject_sem(
        np.array(
            [
                [1.0, np.nan],
                [3.0, np.nan],
                [np.nan, np.nan],
            ]
        )
    )
    np.testing.assert_allclose(sem[0], 1.0)
    assert np.isnan(sem[1])


def test_occlusion_frame_to_second_conversion_is_one_based():
    """Occlusion frame conversion follows (frame - 1) / Hz."""
    cfg = {
        "matlab_drsa": {
            "occlusion_timing": {
                "frame_frequency_hz": 120.0,
                "first_complete_occlusion_frame": 130,
                "reappearance_frame": 191,
            }
        }
    }
    timing = occlusion_timing_from_config(cfg)

    assert frame_to_seconds(130, 120) == 129 / 120
    assert timing["first_complete_occlusion_sec"] == 129 / 120
    assert timing["reappearance_sec"] == 190 / 120


def test_stage07_model_names_replace_model_index_slice_labels(tmp_path):
    """A_dRSA model slices use MATLAB model_names instead of opaque modelXX labels."""
    result_path = tmp_path / "dRSA_occluded_deviant_corr_results.mat"
    matrix = np.zeros((4, 4, 6), dtype=float)
    model_names = np.array(
        [
            "actual deviant position",
            "actual deviant direction",
            "expected/predicted position",
            "expected/predicted direction",
            "prediction error position",
            "prediction error direction",
        ],
        dtype=object,
    )
    savemat(
        result_path,
        {
            "dRSA_neural_correlation": matrix,
            "time_axis_sec": np.arange(4, dtype=float) / 10.0,
            "model_names": model_names,
        },
    )

    matrices, _time_axis, metadata = read_drsa_result_file(result_path)
    slices = list(
        iter_matrix_slices(
            matrices["dRSA_neural_correlation"],
            slice_labels=metadata["model_names"],
        )
    )

    assert slices[3][1] == "expected/predicted direction"
    assert (
        model_source_for_matrix_slice(
            "occluded_deviant",
            "observed",
            "dRSA_neural_correlation",
            slice_label=slices[3][1],
            analysis_family="A_dRSA",
        )
        == "predicted"
    )
    assert (
        model_feature_label(
            "dRSA_neural_correlation",
            slice_label=slices[3][1],
            analysis_family="A_dRSA",
        )
        == "direction"
    )


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
    value_offset: float,
    *,
    condition: str = "always_visible",
) -> None:
    """Write one tiny condition-specific Stage-07-like MAT file for CLI tests."""
    result_path = (
        tmp_path
        / "derivatives"
        / "analysis"
        / subject
        / "ses-01"
        / "drsa"
        / "matlab_real_neural"
        / "A_dRSA_alt"
        / condition
        / "corr"
        / "results"
        / f"A_dRSA_alt_{condition}_corr_occipital.mat"
    )
    result_path.parent.mkdir(parents=True, exist_ok=True)
    matrix = np.zeros((6, 6, 1), dtype=float)
    for model_idx in range(5):
        matrix[model_idx, model_idx + 1, 0] = value_offset + model_idx
    savemat(
        result_path,
        {
            "dRSA_neural_position": matrix,
            "time_axis_sec": np.arange(6, dtype=float) / 10.0,
        },
    )


def _write_fake_condition_set(tmp_path: Path, subject: str, value_offset: float) -> None:
    """Write matched fake Stage-07 MAT files for all real-neural conditions."""
    for offset, condition in enumerate(
        ("always_visible", "occluded_nondeviant", "occluded_deviant")
    ):
        _write_fake_result_mat(
            tmp_path,
            subject,
            value_offset + float(offset),
            condition=condition,
        )


def _write_dummy_png(path: Path) -> None:
    """Write tiny image-like bytes for report embedding tests."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\n")


def test_report_omits_real_neural_single_condition_interval_plots(tmp_path):
    """Reports keep real-neural overlays focused while leaving simulations unchanged."""
    report = tmp_path / "reports" / "line_report.html"
    multi_condition_png = tmp_path / "figures" / "real_multi.png"
    real_single_png = tmp_path / "figures" / "real_single.png"
    simulation_single_png = tmp_path / "figures" / "simulation_single.png"
    for image_path in (multi_condition_png, real_single_png, simulation_single_png):
        _write_dummy_png(image_path)

    STAGE08._write_report(
        report,
        title="Synthetic report",
        summary_lines=["summary"],
        plot_entries=[
            {
                "plot_type": "interval_lag_curve",
                "figure": str(multi_condition_png),
                "analysis_family": "A_dRSA_alt",
                "conditions": [
                    "always_visible",
                    "occluded_nondeviant",
                    "occluded_deviant",
                ],
                "caption": "real multi-condition overlay",
            },
            {
                "plot_type": "interval_lag_curve",
                "figure": str(real_single_png),
                "analysis_family": "A_dRSA",
                "conditions": ["occluded_deviant"],
                "caption": "real single-condition plot",
            },
            {
                "plot_type": "interval_lag_curve",
                "figure": str(simulation_single_png),
                "analysis_family": "simulations_model_vs_model",
                "conditions": ["occluded_deviant"],
                "caption": "simulation single-condition plot",
            },
        ],
        skipped_inputs=[],
        raw_root=tmp_path / "raw",
    )

    report_text = report.read_text(encoding="utf-8")
    embedded_text = report.with_name("line_report_embedded.html").read_text(
        encoding="utf-8"
    )
    assert "real multi-condition overlay" in report_text
    assert "simulation single-condition plot" in report_text
    assert "real single-condition plot" not in report_text
    assert "Omitted generated figures: 1" in report_text
    assert embedded_text.count("data:image/png;base64,") == 2


def test_cli_runs_on_tiny_fake_derivative_tree_and_embeds_report_images(tmp_path):
    """The CLI can build group-only reports from existing MAT files without MATLAB."""
    site_config = _write_fake_site_config(tmp_path)
    _write_fake_condition_set(tmp_path, "sub-01", 10.0)
    _write_fake_condition_set(tmp_path, "sub-02", 20.0)

    completed = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "08_plot_drsa_line_summaries.py"),
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
            "--sensor-groups",
            "occipital",
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
        / "08_drsa_line_summaries_report.html"
    )
    embedded = report.with_name("08_drsa_line_summaries_report_embedded.html")
    assert report.is_file()
    assert embedded.is_file()
    assert "Lag convention: lag = neural time - model time" in report.read_text(
        encoding="utf-8"
    )
    assert "Interval lag-curve mode: model-anchored" in report.read_text(
        encoding="utf-8"
    )
    assert "actual visible position" in report.read_text(encoding="utf-8")
    assert "peak_lag_over_model_time" not in report.read_text(encoding="utf-8")
    assert "Peak-lag" not in report.read_text(encoding="utf-8")
    assert "data:image/png;base64," in embedded.read_text(encoding="utf-8")
    line_dir = (
        tmp_path
        / "derivatives"
        / "analysis"
        / "group"
        / "ses-01"
        / "drsa"
        / "line_summaries"
    )
    assert line_dir.is_dir()
    assert sorted((line_dir / "interval_lag_curves").glob("*.png"))
    assert not (line_dir / "peak_lag_over_model_time").exists()
    assert not sorted(line_dir.rglob("*.csv"))


def test_cli_with_within_subject_shading_reports_selected_mode(tmp_path):
    """The CLI can build shaded interval plots from a tiny fake derivative tree."""
    site_config = _write_fake_site_config(tmp_path)
    _write_fake_condition_set(tmp_path, "sub-01", 10.0)
    _write_fake_condition_set(tmp_path, "sub-02", 20.0)

    completed = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "08_plot_drsa_line_summaries.py"),
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
            "--sensor-groups",
            "occipital",
            "--shading",
            "within-subject-diagonal-sd",
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
        / "08_drsa_line_summaries_report.html"
    )
    report_text = report.read_text(encoding="utf-8")
    assert "Interval lag-curve shading mode: within-subject-diagonal-sd" in report_text
    assert "descriptive within-subject diagonal variability" in report_text
    shaded_dir = (
        tmp_path
        / "derivatives"
        / "analysis"
        / "group"
        / "ses-01"
        / "drsa"
        / "line_summaries"
        / "interval_lag_curves"
        / "shading-within-subject-diagonal-sd"
    )
    assert sorted(shaded_dir.glob("*.png"))
