"""
test_drsa_stats.py

Purpose
Validate the optional Stage-10 dRSA statistics calculations and CLI without
running MATLAB. The tests use synthetic subject x feature arrays and a tiny
fake derivative tree with scipy-written Stage-07-like MAT files.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_drsa_stats.py -q
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat


REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = REPO_ROOT / "meg" / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))


def _load_stage10_module():
    """Load the numbered Stage-10 module for direct function tests."""
    module_path = CODE_DIR / "10_drsa_stats.py"
    spec = importlib.util.spec_from_file_location("stage10_drsa_stats", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STAGE10 = _load_stage10_module()


def _known_cluster_observations(shape=(8, 6, 6)):
    """Create a subject stack with a positive 2-D cluster and weak background."""
    rng = np.random.default_rng(42)
    data = rng.normal(0.0, 0.03, size=shape)
    for subject_index in range(shape[0]):
        data[subject_index, 2:4, 3:5] += 0.9 + subject_index * 0.03
    return data


def test_2d_one_sample_cluster_detects_known_positive_cluster():
    """The NumPy backend detects a known positive 2-D cluster."""
    data = _known_cluster_observations()

    result = STAGE10.cluster_permutation_1samp(
        data,
        n_permutations=256,
        cluster_alpha=0.05,
        significance_alpha=0.05,
        tail="positive",
        seed=1,
        backend="numpy",
    )

    known = np.zeros((6, 6), dtype=bool)
    known[2:4, 3:5] = True
    assert result.clusters
    assert any(
        cluster.p_value <= 0.05 and np.any(cluster.mask & known)
        for cluster in result.clusters
    )


def test_2d_paired_contrast_detects_known_difference_and_drops_unmatched():
    """Paired matrix preparation drops subjects lacking either contrast side."""
    subjects = tuple(f"sub-{index:02d}" for index in range(1, 8))
    records = []
    axis = np.arange(6, dtype=float) / 10.0
    for subject_index, subject in enumerate(subjects[:6]):
        nondeviant = np.random.default_rng(subject_index).normal(0.0, 0.02, size=(6, 6))
        deviant = nondeviant.copy()
        deviant[2:4, 3:5] += 1.0 + subject_index * 0.02
        for condition, matrix in (
            ("occluded_nondeviant", nondeviant),
            ("occluded_deviant", deviant),
        ):
            records.append(
                STAGE10.MatrixRecord(
                    subject=subject,
                    condition=condition,
                    analysis_family="A_dRSA_alt",
                    sensor_group="occipital",
                    drsa_type="corr",
                    model_source="observed",
                    model_feature="position",
                    stimulus_model_label="actual position",
                    variable="dRSA_neural_position",
                    slice_index=0,
                    slice_label="matrix",
                    variant="raw",
                    matrix=matrix,
                    time_axis_sec=axis,
                    relative_result=f"{condition}/fake.mat",
                    source_path=f"/tmp/{subject}/{condition}.mat",
                )
            )
    records.append(
        STAGE10.MatrixRecord(
            subject="sub-07",
            condition="occluded_deviant",
            analysis_family="A_dRSA_alt",
            sensor_group="occipital",
            drsa_type="corr",
            model_source="observed",
            model_feature="position",
            stimulus_model_label="actual deviant position",
            variable="dRSA_neural_position",
            slice_index=0,
            slice_label="matrix",
            variant="raw",
            matrix=np.ones((6, 6), dtype=float),
            time_axis_sec=axis,
            relative_result="occluded_deviant/fake.mat",
            source_path="/tmp/sub-07/occluded_deviant.mat",
        )
    )

    skipped = []
    prepared = STAGE10._prepare_matrix_paired_tests(  # pylint: disable=protected-access
        records,
        subjects,
        skipped=skipped,
    )

    assert len(prepared) == 1
    assert prepared[0].subjects == subjects[:6]
    assert "sub-07" in prepared[0].dropped_subjects
    result = STAGE10.cluster_permutation_1samp(
        prepared[0].observations,
        n_permutations=128,
        cluster_alpha=0.05,
        significance_alpha=0.05,
        tail="positive",
        seed=2,
        backend="numpy",
    )
    assert any(cluster.p_value <= 0.05 for cluster in result.clusters)


def test_1d_paired_line_summary_uses_difference_matrix_before_lag_averaging():
    """Paired line prep computes lag curves from subject difference matrices."""
    subjects = ("sub-01", "sub-02")
    axis = np.arange(4, dtype=float) / 10.0
    records = []
    expected_zero_lag = []
    for subject, dev_diag, non_diag in (
        (
            "sub-01",
            (10.0, np.nan, np.nan, np.nan),
            (1.0, 100.0, np.nan, np.nan),
        ),
        (
            "sub-02",
            (20.0, np.nan, np.nan, np.nan),
            (2.0, 200.0, np.nan, np.nan),
        ),
    ):
        deviant = np.full((4, 4), np.nan, dtype=float)
        nondeviant = np.full((4, 4), np.nan, dtype=float)
        for index, (dev_value, non_value) in enumerate(zip(dev_diag, non_diag)):
            deviant[index, index] = dev_value
            nondeviant[index, index] = non_value
        expected_zero_lag.append(deviant[0, 0] - nondeviant[0, 0])

        # The asymmetric NaN support on the diagonal makes curve-minus-curve
        # differ from computing the curve on the paired difference matrix.
        for condition, matrix in (
            ("occluded_nondeviant", nondeviant),
            ("occluded_deviant", deviant),
        ):
            records.append(
                STAGE10.MatrixRecord(
                    subject=subject,
                    condition=condition,
                    analysis_family="A_dRSA_alt",
                    sensor_group="occipital",
                    drsa_type="corr",
                    model_source="observed",
                    model_feature="position",
                    stimulus_model_label="actual position",
                    variable="dRSA_neural_position",
                    slice_index=0,
                    slice_label="matrix",
                    variant="raw",
                    matrix=matrix,
                    time_axis_sec=axis,
                    relative_result=f"{condition}/fake.mat",
                    source_path=f"/tmp/{subject}/{condition}.mat",
                )
            )
    args = STAGE10.build_arg_parser().parse_args(
        [
            "--config",
            "fake.yaml",
            "--line-intervals",
            "pre_occlusion",
            "--line-lag-window-sec",
            "0.1",
        ]
    )

    skipped = []
    prepared = STAGE10._prepare_line_paired_tests(  # pylint: disable=protected-access
        records,
        subjects,
        args,
        {},
        skipped=skipped,
    )

    assert not skipped
    assert len(prepared) == 1
    zero_lag_index = int(np.flatnonzero(np.isclose(prepared[0].lag_sec, 0.0))[0])
    np.testing.assert_allclose(
        prepared[0].observations[:, zero_lag_index],
        np.asarray(expected_zero_lag, dtype=float),
    )
    np.testing.assert_allclose(
        prepared[0].mean_map[zero_lag_index],
        np.mean(expected_zero_lag),
    )


def test_1d_lag_curve_cluster_detects_known_lag_cluster():
    """The 1-D cluster test detects a known contiguous lag effect."""
    rng = np.random.default_rng(7)
    data = rng.normal(0.0, 0.03, size=(8, 11))
    data[:, 5:8] += np.linspace(0.8, 1.0, 8)[:, None]

    result = STAGE10.cluster_permutation_1samp(
        data,
        n_permutations=256,
        cluster_alpha=0.05,
        significance_alpha=0.05,
        tail="positive",
        seed=3,
        backend="numpy",
    )

    known = np.zeros(11, dtype=bool)
    known[5:8] = True
    assert any(
        cluster.p_value <= 0.05 and np.any(cluster.mask & known)
        for cluster in result.clusters
    )


def test_mne_1d_cluster_masks_are_normalized():
    """MNE 1-D slice clusters are converted to masks matching t_obs."""
    rng = np.random.default_rng(11)
    data = rng.normal(0.0, 0.05, size=(7, 31))
    data[:, 12:16] += 1.0

    result = STAGE10.cluster_permutation_1samp(
        data,
        n_permutations=32,
        cluster_alpha=0.2,
        significance_alpha=0.05,
        tail="positive",
        seed=12,
        backend="mne",
    )

    assert result.backend == "mne"
    assert result.clusters
    assert all(cluster.mask.shape == result.statistic_map.shape for cluster in result.clusters)


def test_auto_backend_uses_numpy_for_zero_variance_features():
    """Auto avoids MNE when constant features would create infinite t values."""
    data = np.ones((6, 5), dtype=float)

    result = STAGE10.cluster_permutation_1samp(
        data,
        n_permutations=32,
        cluster_alpha=0.2,
        significance_alpha=0.05,
        tail="positive",
        seed=13,
        backend="auto",
    )

    assert result.backend == "numpy"
    assert np.isfinite(result.statistic_map).all()


def test_cli_defaults_exclude_simulations_and_enable_actual_position_contrast():
    """Stage-10 defaults keep simulations out and enable the actual-position contrast."""
    args = STAGE10.build_arg_parser().parse_args(["--config", "fake.yaml"])

    assert args.analysis_families == ["A_dRSA", "A_dRSA_alt"]
    assert args.paired_contrasts == ["actual_deviant_position_vs_nondeviant"]
    assert args.report_cluster_tables == "significant"


def test_report_cluster_table_defaults_to_significant_rows_only():
    """HTML report tables show only significant clusters by default."""
    rows = [
        {
            "cluster_id": 1,
            "p_value": 0.01,
            "cluster_size": 3,
            "cluster_sign": "positive",
            "lag_min_sec": -0.1,
            "lag_max_sec": 0.1,
        },
        {
            "cluster_id": 2,
            "p_value": 0.2,
            "cluster_size": 4,
            "cluster_sign": "positive",
            "lag_min_sec": 0.2,
            "lag_max_sec": 0.3,
        },
    ]

    significant = STAGE10._cluster_table_html(  # pylint: disable=protected-access
        rows,
        mode="significant",
        significance_alpha=0.05,
    )
    all_rows = STAGE10._cluster_table_html(  # pylint: disable=protected-access
        rows,
        mode="all",
        significance_alpha=0.05,
    )

    assert "<td>1</td>" in significant
    assert "<td>2</td>" not in significant
    assert "<td>1</td>" in all_rows
    assert "<td>2</td>" in all_rows


def test_nan_values_are_handled_safely():
    """NaN-containing arrays use finite support and keep invalid features out."""
    data = _known_cluster_observations()
    data[:2, 0, 0] = np.nan
    data[:, 5, 5] = np.nan

    result = STAGE10.cluster_permutation_1samp(
        data,
        n_permutations=64,
        cluster_alpha=0.05,
        significance_alpha=0.05,
        tail="positive",
        seed=4,
        backend="numpy",
    )

    assert np.isfinite(result.statistic_map).all()
    assert result.finite_observation_counts[0, 0] == data.shape[0] - 2
    assert result.finite_observation_counts[5, 5] == 0


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
    value_offset: float,
) -> None:
    """Write one tiny Stage-07-like occipital A_dRSA_alt MAT file."""
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


def _write_fake_condition_set(tmp_path: Path, subject: str, value_offset: float) -> None:
    """Write matched fake Stage-07 MAT files for all real-neural conditions."""
    for offset, condition in enumerate(
        ("always_visible", "occluded_nondeviant", "occluded_deviant")
    ):
        _write_fake_result_mat(
            tmp_path,
            subject,
            condition=condition,
            value_offset=value_offset + float(offset),
        )


def test_cli_runs_on_tiny_fake_derivative_tree_and_embeds_report(tmp_path):
    """The CLI writes stats outputs and embedded HTML without MATLAB."""
    site_config = _write_fake_site_config(tmp_path)
    _write_fake_condition_set(tmp_path, "sub-01", 10.0)
    _write_fake_condition_set(tmp_path, "sub-02", 20.0)

    completed = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "10_drsa_stats.py"),
            "--neural-normalization", "none",
            "--config",
            str(REPO_ROOT / "meg" / "params" / "pipeline_default.yaml"),
            "--site-config",
            str(site_config),
            "--subjects",
            "01",
            "02",
            "--sensor-groups",
            "occipital",
            "--analysis-families",
            "A_dRSA_alt",
            "--stats-targets",
            "matrices",
            "line-summaries",
            "--line-intervals",
            "pre_occlusion",
            "--line-lag-window-sec",
            "0.2",
            "--include-smoothed",
            "--n-permutations",
            "16",
            "--cluster-alpha",
            "0.2",
            "--stats-backend",
            "numpy",
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
        / "10_drsa_stats_report.html"
    )
    embedded = report.with_name("10_drsa_stats_report_embedded.html")
    stats_dir = (
        tmp_path
        / "derivatives"
        / "analysis"
        / "group"
        / "ses-01"
        / "drsa"
        / "stats"
    )
    assert report.is_file()
    assert embedded.is_file()
    assert (stats_dir / "10_drsa_stats_summary.json").is_file()
    assert (stats_dir / "10_drsa_stats_cluster_summary.csv").is_file()
    assert sorted(stats_dir.rglob("*.npz"))
    report_text = report.read_text(encoding="utf-8")
    embedded_text = embedded.read_text(encoding="utf-8")
    assert "Lag convention: lag = neural time - model time" in report_text
    assert "Model-anchored interval behavior" in report_text
    assert "target=matrices" in report_text
    assert "target=line-summaries" in report_text
    assert "Paired contrasts: actual_deviant_position_vs_nondeviant" in report_text
    assert "<a href=\"#variant-raw\">Raw inputs</a>" in report_text
    assert "<a href=\"#variant-smoothed-display\">Smoothed-display inputs</a>" in report_text
    assert "<a href=\"#raw-condition-always-visible\">Always visible</a>" in report_text
    assert "<a href=\"#raw-condition-occluded-nondeviant\">Occluded nondeviant</a>" in report_text
    assert "<a href=\"#raw-condition-occluded-deviant\">Occluded deviant</a>" in report_text
    assert "<a href=\"#raw-contrasts\">Contrasts</a>" in report_text
    assert "actual deviant position - actual nondeviant position" in report_text
    assert "Sensor group: occipital" in report_text
    assert "data:image/png;base64," in embedded_text
