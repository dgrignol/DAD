"""
test_dot_side_group_report.py

Purpose
Validate the group-level dot-side decoding report helpers with synthetic
arrays and temporary subject outputs, without running participant decoding.

Usage
Run only these tests:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_dot_side_group_report.py -q
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "meg" / "code" / "sanity_checks" / "02_dot_side_group_report.py"

SPEC = importlib.util.spec_from_file_location("dot_side_group_report", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
group_report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(group_report)


def _synthetic_training_testing_diagnostics(
    *, n_decoder_trials: int = 8, n_decoder_identities: int = 4
) -> dict[str, object]:
    """Return compact decoder diagnostics for synthetic report fixtures."""
    return {
        "classifier": {
            "name": "Regularized linear discriminant analysis",
            "details": "Train-fold z-scored binary LDA with ridge regularization.",
        },
        "cross_validation": {
            "scheme": "stratified_k_fold",
            "n_splits": 5,
            "split_unit": "same-identity pseudo-trials",
        },
        "evaluation": {
            "unbalanced_metric": "accuracy_unbalanced_train",
            "balanced_metric": "accuracy_balanced_undersampled_sets",
            "balancing_strategy": "Random undersampling inside each CV fold.",
        },
        "sample_structure": {
            "averaging_mode": "fixed_count",
            "requested_value": "1",
            "n_source_trials": 17,
            "n_decoder_trials": n_decoder_trials,
            "n_decoder_identities": n_decoder_identities,
        },
        "class_balance": {
            "n_windows": 6,
            "n_decoded_windows": 6,
            "n_skipped_windows": 0,
            "decoded_window_fraction": 1.0,
            "skip_reasons": {"decoded": 6},
            "min_minority_class_count_all_windows": 5,
            "max_majority_to_minority_ratio": 1.2,
        },
        "assessment": {
            "status": "well_defined",
            "text": "Every decoded center had enough samples for stratified CV.",
        },
    }


def test_cluster_permutation_recovers_known_positive_cluster():
    """A strong contiguous above-chance interval should form a significant cluster."""
    values = np.full((7, 12), 0.5, dtype=float)
    values[:, 4:8] += np.array([0.12, 0.14, 0.13, 0.16, 0.15, 0.11, 0.17])[:, None]
    result = group_report.cluster_permutation_1d(
        values,
        chance=0.5,
        cluster_forming_p=0.05,
        cluster_alpha=0.05,
        tail="greater",
        n_permutations=0,
        random_seed=1,
    )

    assert result["permutation_mode"] == "exact"
    assert result["n_permutations"] == 128
    significant = [cluster for cluster in result["clusters"] if cluster["significant"]]
    assert significant
    assert significant[0]["start_index"] == 4
    assert significant[0]["end_index"] == 7
    assert significant[0]["p_value"] <= 0.05


def test_ci95_is_default_interval_and_sem_is_alternative():
    """The group display band supports 95% CI and SEM modes."""
    values = np.array(
        [
            [0.55, 0.60],
            [0.57, 0.62],
            [0.59, 0.64],
            [0.61, 0.66],
        ],
        dtype=float,
    )
    ci = group_report.group_interval(values, interval="ci95")
    sem = group_report.group_interval(values, interval="sem")

    assert ci["label"] == "95% CI"
    assert sem["label"] == "mean +/- SEM"
    np.testing.assert_allclose(ci["mean"], sem["mean"])
    assert np.all((ci["high"] - ci["mean"]) > (sem["high"] - sem["mean"]))


def test_cluster_permutation_tolerates_nan_timepoint_gaps():
    """Skipped subject windows must not abort or bridge temporal clusters."""
    values = np.full((4, 8), 0.5, dtype=float)
    values[:, 1:3] = 0.65
    values[:, 4:6] = 0.66
    values[:, 3] = 0.67
    values[0, 3] = np.nan

    result = group_report.cluster_permutation_1d(
        values,
        chance=0.5,
        cluster_forming_p=0.05,
        cluster_alpha=0.05,
        tail="greater",
        n_permutations=0,
        random_seed=1,
    )

    assert result["n_valid_timepoints"] == 7
    assert result["n_invalid_timepoints"] == 1
    assert np.isnan(result["observed_t"][3])
    assert [(cluster["start_index"], cluster["end_index"]) for cluster in result["clusters"]] == [
        (1, 2),
        (4, 5),
    ]


def test_1d_smoothing_preserves_nan_gaps():
    """Optional group smoothing must not leak signal across skipped windows."""
    values = np.array([[0.5, 0.5, 0.8, np.nan, 0.5, 0.5, 0.5]], dtype=float)
    smoothed = group_report.smooth_group_values_1d(values, sigma_samples=1.0)

    assert np.isnan(smoothed[0, 3])
    assert smoothed[0, 1] > values[0, 1]
    assert smoothed[0, 4] == pytest.approx(values[0, 4])
    np.testing.assert_allclose(
        group_report.smooth_group_values_1d(values, sigma_samples=0.0),
        values,
        equal_nan=True,
    )


def _write_subject_output(
    root: Path,
    subject: str,
    times: np.ndarray,
    values: np.ndarray,
    *,
    condition: str = "always_visible",
    model_source: str = "observed",
) -> None:
    """Write a minimal subject dot-side output bundle for input-loader tests."""
    out_dir = (
        root
        / "sanity_checks"
        / subject
        / "ses-01"
        / "dot_side_decoding"
        / condition
        / model_source
    )
    out_dir.mkdir(parents=True)
    stem = f"{condition}_{model_source}"
    scores_csv = out_dir / f"{stem}_dot_side_decoding_scores.csv"
    figure = out_dir / f"{stem}_dot_side_decoding_3panel.png"
    summary_json = out_dir / f"{stem}_dot_side_decoding_summary.json"
    pd.DataFrame(
        {
            "time_sec": times,
            "accuracy_balanced_undersampled_sets": values,
            "accuracy_unbalanced_train": values + 0.01,
        }
    ).to_csv(scores_csv, index=False)
    figure.write_bytes(b"\x89PNG\r\n\x1a\n")
    summary_json.write_text(
        json.dumps(
            {
                "condition": condition,
                "model_source": model_source,
                "model_condition": (
                    "occluded_nondeviant"
                    if model_source == "predicted"
                    else condition
                ),
                "analysis_key": stem,
                "subject": subject.removeprefix("sub-"),
                "session": "01",
                "identity_range": [1, 4],
                "n_times": int(times.size),
                "sfreq_hz": 120.0,
                "decode_window_samples": 10,
                "decode_tmin_sec": float(times[0]),
                "decode_tmax_sec": float(times[-1]),
                "n_splits": 5,
                "random_seed": 42,
                "trial_averaging": {
                    "mode": "fixed_count",
                    "requested_value": "1",
                    "trials_per_average": 1,
                    "n_source_trials_by_identity": {
                        "1": 5,
                        "2": 4,
                        "3": 5,
                        "4": 3,
                    },
                },
                "training_testing_diagnostics": _synthetic_training_testing_diagnostics(),
                "scores_csv": str(scores_csv),
                "figure_path": str(figure),
                "peak_accuracy_balanced_train": float(np.max(values)),
                "peak_time_sec_balanced_train": float(times[int(np.argmax(values))]),
            }
        ),
        encoding="utf-8",
    )


def test_standard_analysis_registry_contains_requested_four_analyses():
    """The combined report registry must contain the four scientific analyses."""
    args = group_report.build_arg_parser().parse_args(["--analyses", "all"])
    specs = group_report.resolve_analysis_specs(args)
    assert [spec["key"] for spec in specs] == [
        "always_visible_observed",
        "occluded_nondeviant_observed",
        "occluded_deviant_observed",
        "occluded_deviant_predicted",
    ]
    predicted = specs[-1]
    assert predicted["condition"] == "occluded_deviant"
    assert predicted["model_condition"] == "occluded_nondeviant"


def test_combined_report_rejects_cross_analysis_setting_mismatch():
    """A combined report must not mix incompatible decoding configurations."""
    base_summary = {
        "sfreq_hz": 120.0,
        "decode_window_samples": 10,
        "decode_tmin_sec": 0.0,
        "decode_tmax_sec": 2.667,
        "n_splits": 5,
        "random_seed": 42,
        "trial_averaging": {"trials_per_average": 1},
    }
    results = [
        {
            "key": "always_visible_observed",
            "participant_records": [{"summary": dict(base_summary)}],
        },
        {
            "key": "occluded_nondeviant_observed",
            "participant_records": [
                {
                    "summary": {
                        **base_summary,
                        "decode_window_samples": 12,
                    }
                }
            ],
        },
    ]
    with pytest.raises(RuntimeError, match="Cross-analysis decoding setting mismatch"):
        group_report.validate_cross_analysis_settings(results)


def test_combined_html_places_all_group_analyses_before_subjects(tmp_path):
    """No participant-specific section may appear before all group analyses."""
    figure = tmp_path / "figure.png"
    figure.write_bytes(b"\x89PNG\r\n\x1a\n")
    analysis_results = []
    for spec in group_report.STANDARD_ANALYSES:
        analysis_results.append(
            {
                **spec,
                "group_figure": str(figure),
                "clusters": [],
                "skipped_inputs": [],
                "participant_records": [
                    {
                        "subject": "sub-01",
                        "figure": str(figure),
                        "summary": {
                            "peak_accuracy_balanced_train": 0.6,
                            "peak_time_sec_balanced_train": 0.2,
                            "training_testing_diagnostics": (
                                _synthetic_training_testing_diagnostics()
                            ),
                        },
                    }
                ],
            }
        )

    report_path = group_report.write_combined_html_report(
        tmp_path / "combined.html",
        title="Combined test",
        summary_lines=["test"],
        subjects=["sub-01"],
        analysis_results=analysis_results,
        participant_trial_count_figures={"sub-01": str(figure)},
        raw_root=tmp_path / "raw",
        embed_images=False,
    )
    report_text = report_path.read_text(encoding="utf-8")
    assert report_text.index('id="group-analyses"') < report_text.index(
        'id="participants"'
    )
    assert "Trial counts per identity" in report_text
    assert "Training/testing diagnostics" in report_text
    assert "Regularized linear discriminant analysis" in report_text
    for spec in group_report.STANDARD_ANALYSES:
        assert report_text.index(f'id="group-{spec["key"]}"') < report_text.index(
            'id="participants"'
        )


def test_group_input_loader_reports_missing_subjects(tmp_path):
    """Missing requested subjects fail unless explicitly allowed."""
    times = np.arange(5, dtype=float) / 120.0
    _write_subject_output(tmp_path, "sub-01", times, np.full(5, 0.55))
    cfg = {
        "resolved_paths": {"derivatives_root": str(tmp_path)},
        "layout": {"session_bids_id": "ses-01"},
    }

    with pytest.raises(RuntimeError, match="Missing or invalid subject"):
        group_report.load_group_inputs(
            cfg=cfg,
            subjects=["sub-01", "sub-02"],
            condition="always_visible",
            metric="accuracy_balanced_undersampled_sets",
            allow_missing=False,
        )
    with pytest.raises(RuntimeError, match="At least two valid"):
        group_report.load_group_inputs(
            cfg=cfg,
            subjects=["sub-01", "sub-02"],
            condition="always_visible",
            metric="accuracy_balanced_undersampled_sets",
            allow_missing=True,
        )


def test_group_input_loader_rejects_time_axis_mismatch(tmp_path):
    """Available subjects must share the same decoded time grid."""
    _write_subject_output(tmp_path, "sub-01", np.arange(5, dtype=float), np.full(5, 0.55))
    _write_subject_output(tmp_path, "sub-02", np.arange(6, dtype=float), np.full(6, 0.56))
    cfg = {
        "resolved_paths": {"derivatives_root": str(tmp_path)},
        "layout": {"session_bids_id": "ses-01"},
    }

    with pytest.raises(RuntimeError, match="time axis mismatch"):
        group_report.load_group_inputs(
            cfg=cfg,
            subjects=["sub-01", "sub-02"],
            condition="always_visible",
            metric="accuracy_balanced_undersampled_sets",
            allow_missing=False,
        )


def test_combined_group_report_builds_all_analyses_and_embedded_html(tmp_path):
    """Synthetic subject bundles should produce one complete combined report."""
    times = np.arange(6, dtype=float) / 120.0
    for analysis_idx, spec in enumerate(group_report.STANDARD_ANALYSES):
        for subject_idx, subject in enumerate(("sub-01", "sub-02")):
            values = (
                np.linspace(0.53, 0.61, times.size)
                + (analysis_idx * 0.005)
                + (subject_idx * 0.01)
            )
            _write_subject_output(
                tmp_path,
                subject,
                times,
                values,
                condition=spec["condition"],
                model_source=spec["model_source"],
            )

    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    cfg = {
        "resolved_paths": {
            "derivatives_root": str(tmp_path),
            "data_root": str(raw_root),
        },
        "layout": {"session_bids_id": "ses-01"},
    }
    args = group_report.build_arg_parser().parse_args(
        [
            "--analyses",
            "all",
            "--group-subjects",
            "01",
            "02",
            "--overwrite",
        ]
    )
    report = group_report.build_group_report(cfg, args)

    assert report["status"] == "completed"
    assert report["analysis_keys"] == [
        spec["key"] for spec in group_report.STANDARD_ANALYSES
    ]
    embedded = Path(report["report_files"]["embedded_html"])
    assert embedded.is_file()
    html_text = embedded.read_text(encoding="utf-8")
    assert html_text.index('id="group-analyses"') < html_text.index(
        'id="participants"'
    )
    for analysis in report["analyses"].values():
        assert Path(analysis["group_figure"]).is_file()
    count_figures = report["participant_trial_count_figures"]
    assert set(count_figures) == {"sub-01", "sub-02"}
    assert all(Path(path).is_file() for path in count_figures.values())
    assert "Trial counts per identity" in html_text
    assert "Training/testing diagnostics" in html_text
    assert "well_defined" in html_text


def test_occluded_group_plot_markers_use_configured_seconds():
    """Group accuracy axes mark both configured occlusion boundaries."""
    fig, ax = group_report.plt.subplots()
    timing = {
        "first_complete_occlusion_sec": 129 / 120,
        "reappearance_sec": 190 / 120,
    }
    group_report.add_occlusion_boundary_lines(
        ax,
        condition="occluded_nondeviant",
        occlusion_timing=timing,
    )
    lines_by_label = {line.get_label(): line for line in ax.lines}
    assert lines_by_label["First complete occlusion"].get_xdata()[0] == pytest.approx(
        129 / 120
    )
    assert lines_by_label["First reappearance"].get_xdata()[0] == pytest.approx(
        190 / 120
    )
    group_report.plt.close(fig)
