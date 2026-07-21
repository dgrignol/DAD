"""
test_dot_side_decoding.py

Purpose
Validate helper contracts for the optional dot-side decoding sanity check
without loading participant MEG data or running the full report workflow.

Usage
Run only these tests:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_dot_side_decoding.py -q
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "meg" / "code" / "sanity_checks" / "01_dot_side_decoding.py"

SPEC = importlib.util.spec_from_file_location("dot_side_decoding", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
dot_side = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dot_side)


def test_default_ten_sample_windows_clip_edges():
    """The default 10-sample window is interpreted on the resampled grid."""
    times = np.arange(12, dtype=float) / 120.0
    rows = dot_side.build_centered_window_metadata(
        times=times,
        window_samples=10,
        sfreq=120.0,
    )

    assert rows[0]["window_start_sample_inclusive"] == 0
    assert rows[0]["window_end_sample_exclusive"] == 6
    assert rows[0]["effective_window_samples"] == 6
    assert rows[0]["is_edge_clipped"] is True
    assert rows[5]["window_start_sample_inclusive"] == 1
    assert rows[5]["window_end_sample_exclusive"] == 11
    assert rows[5]["effective_window_samples"] == 10
    assert rows[5]["effective_window_ms"] == pytest.approx(10 / 120 * 1000)


def test_occluded_participant_plot_markers_use_configured_seconds():
    """Participant decoding axes mark complete occlusion and reappearance."""
    fig, ax = dot_side.plt.subplots()
    dot_side.add_occlusion_boundary_lines(
        ax,
        condition="occluded_deviant",
        occlusion_timing={
            "first_complete_occlusion_sec": 129 / 120,
            "reappearance_sec": 190 / 120,
        },
    )
    lines_by_label = {line.get_label(): line for line in ax.lines}
    assert lines_by_label["First complete occlusion"].get_xdata()[0] == pytest.approx(
        129 / 120
    )
    assert lines_by_label["First reappearance"].get_xdata()[0] == pytest.approx(
        190 / 120
    )
    dot_side.plt.close(fig)


def test_odd_and_even_window_centering_policy():
    """Odd windows are symmetric; even windows place the extra sample after center."""
    times = np.arange(9, dtype=float)
    odd = dot_side.build_centered_window_metadata(
        times=times,
        window_samples=5,
        sfreq=1.0,
    )
    even = dot_side.build_centered_window_metadata(
        times=times,
        window_samples=4,
        sfreq=1.0,
    )

    assert odd[4]["window_start_sample_inclusive"] == 2
    assert odd[4]["window_end_sample_exclusive"] == 7
    assert odd[4]["samples_before_center_requested"] == 2
    assert odd[4]["samples_after_center_requested"] == 2
    assert even[4]["window_start_sample_inclusive"] == 3
    assert even[4]["window_end_sample_exclusive"] == 7
    assert even[4]["samples_before_center_requested"] == 1
    assert even[4]["samples_after_center_requested"] == 2


def test_same_identity_trial_averaging_drops_remainders_and_is_seeded():
    """Pseudo-trials average only complete shuffled groups within identity."""
    data = np.arange(7 * 2 * 3, dtype=float).reshape(7, 2, 3)
    labels = np.vstack(
        [
            np.zeros((3, 3), dtype=int),
            np.ones((4, 3), dtype=int),
        ]
    )
    identities = np.array([1, 1, 1, 2, 2, 2, 2])

    avg_a = dot_side.average_trials_within_identity(
        data,
        labels,
        identities=identities,
        trials_per_average=2,
        random_seed=97,
    )
    avg_b = dot_side.average_trials_within_identity(
        data,
        labels,
        identities=identities,
        trials_per_average=2,
        random_seed=97,
    )
    averaged_data, averaged_labels, averaged_identities, metadata = avg_a

    assert averaged_data.shape == (3, 2, 3)
    assert averaged_labels.shape == (3, 3)
    assert averaged_identities.tolist() == [1, 2, 2]
    assert metadata["dropped_remainders_by_identity"] == {"1": 1, "2": 0}
    assert metadata["n_dropped_remainder_trials_total"] == 1
    assert all(
        item["source_epoch_count"] == 2 for item in metadata["pseudo_trials"]
    )
    assert avg_b[3]["pseudo_trials"] == metadata["pseudo_trials"]
    np.testing.assert_allclose(avg_b[0], averaged_data)


def test_all_trials_per_identity_averaging_collapses_each_identity_once():
    """The all-trials mode creates one evoked pseudo-trial per identity."""
    data = np.arange(6 * 2 * 3, dtype=float).reshape(6, 2, 3)
    labels = np.vstack(
        [
            np.zeros((2, 3), dtype=int),
            np.ones((4, 3), dtype=int),
        ]
    )
    identities = np.array([1, 1, 2, 2, 2, 2])

    averaged_data, averaged_labels, averaged_identities, metadata = (
        dot_side.average_trials_within_identity(
            data,
            labels,
            identities=identities,
            trials_per_average="all_trials_per_identity",
            random_seed=42,
            source_epoch_indices=np.arange(10, 16),
        )
    )

    assert averaged_data.shape == (2, 2, 3)
    np.testing.assert_allclose(averaged_data[0], data[:2].mean(axis=0))
    np.testing.assert_allclose(averaged_data[1], data[2:].mean(axis=0))
    assert averaged_labels.tolist() == [[0, 0, 0], [1, 1, 1]]
    assert averaged_identities.tolist() == [1, 2]
    assert metadata["mode"] == "all_trials_per_identity"
    assert metadata["requested_value"] == "all_trials_per_identity"
    assert metadata["trials_per_average"] is None
    assert metadata["n_pseudo_trials_by_identity"] == {"1": 1, "2": 1}
    assert metadata["n_dropped_remainder_trials_total"] == 0
    assert metadata["pseudo_trials"][0]["source_epoch_indices_zero_based"] == [10, 11]
    assert metadata["pseudo_trials"][1]["source_epoch_count"] == 4


def test_trials_per_average_parser_accepts_integer_or_all_mode():
    """The CLI parser accepts fixed counts and the all-trials literal."""
    parser = dot_side.build_arg_parser()
    fixed = parser.parse_args(["--trials-per-average", "3"])
    all_mode = parser.parse_args(
        ["--trials-per-average", "all_trials_per_identity"]
    )

    assert fixed.trials_per_average == 3
    assert all_mode.trials_per_average == "all_trials_per_identity"
    with pytest.raises(SystemExit):
        parser.parse_args(["--trials-per-average", "0"])


def test_stage07_sensor_group_selection_filters_requested_sensor_type():
    """Stage-07 groups keep only requested channel types in MNE channel order."""
    info = dot_side.mne.create_info(
        ["MEG1641", "MEG1642", "MEG1643", "MEG0111"],
        sfreq=100.0,
        ch_types=["mag", "grad", "grad", "mag"],
    )

    selected, metadata = dot_side.select_stage07_decoder_channels(
        info,
        sensor_group="occipital",
        decoder_sensor_type="grad",
    )

    assert selected == ["MEG1642", "MEG1643"]
    assert metadata["sensor_group"] == "occipital"
    assert metadata["decoder_sensor_type"] == "grad"
    assert metadata["n_selected_channels"] == 2


def test_trial_averaging_records_identity_with_no_complete_group():
    """An underfilled identity is dropped without aborting other identities."""
    data = np.arange(5 * 2 * 3, dtype=float).reshape(5, 2, 3)
    labels = np.vstack(
        [
            np.zeros((1, 3), dtype=int),
            np.ones((4, 3), dtype=int),
        ]
    )
    identities = np.array([40, 1, 1, 1, 1])

    averaged_data, _labels, averaged_identities, metadata = (
        dot_side.average_trials_within_identity(
            data,
            labels,
            identities=identities,
            trials_per_average=2,
            random_seed=42,
        )
    )

    assert averaged_data.shape[0] == 2
    assert averaged_identities.tolist() == [1, 1]
    assert metadata["dropped_remainders_by_identity"]["40"] == 1
    assert metadata["dropped_identities_no_complete_group"] == [40]
    assert metadata["retained_identity_values"] == [1]


def test_model_identity_validation_supports_20_and_40_profiles(monkeypatch):
    """Model-derived labels must match the configured identity profile exactly."""
    def fake_loadmat_20(_path, simplify_cells):
        assert simplify_cells is True
        return {
            "Cfg": {"rectSize": [10.0, 8.0]},
            "xySeqs": [
                {
                    "condition_label": "always_visible",
                    "sequence": identity,
                    "xy": np.column_stack(
                        [
                            np.full(6, identity % 10, dtype=float),
                            np.zeros(6, dtype=float),
                        ]
                    ),
                }
                for identity in range(1, 21)
            ],
        }

    def fake_loadmat_40(_path, simplify_cells):
        base = fake_loadmat_20(_path, simplify_cells)
        base["xySeqs"].extend(
            [
                {
                    "condition_label": "always_visible",
                    "sequence": identity,
                    "xy": np.column_stack(
                        [
                            np.full(6, identity % 10, dtype=float),
                            np.zeros(6, dtype=float),
                        ]
                    ),
                }
                for identity in range(21, 41)
            ]
        )
        return base

    monkeypatch.setattr(dot_side.sio, "loadmat", fake_loadmat_20)
    labels_20, x_20, center_20 = dot_side.build_dot_side_labels(
        model_mat=Path("unused.mat"),
        condition="always_visible",
        trial_identities=np.array([1, 20]),
        identity_values=list(range(1, 21)),
        n_times_target=4,
    )
    assert labels_20.shape == (2, 4)
    assert x_20.shape == (20, 4)
    assert center_20 == 5.0
    with pytest.raises(RuntimeError, match="Model identity set mismatch"):
        dot_side.build_dot_side_labels(
            model_mat=Path("unused.mat"),
            condition="always_visible",
            trial_identities=np.array([1, 20]),
            identity_values=list(range(1, 41)),
            n_times_target=4,
        )

    monkeypatch.setattr(dot_side.sio, "loadmat", fake_loadmat_40)
    labels_40, x_40, _center_40 = dot_side.build_dot_side_labels(
        model_mat=Path("unused.mat"),
        condition="always_visible",
        trial_identities=np.array([1, 40]),
        identity_values=list(range(1, 41)),
        n_times_target=4,
    )
    assert labels_40.shape == (2, 4)
    assert x_40.shape == (40, 4)


def test_predicted_deviant_uses_nondeviant_model_condition(monkeypatch):
    """Predicted deviant labels come from nondeviant paths, not deviant paths."""
    assert (
        dot_side.resolve_model_condition("occluded_deviant", "predicted")
        == "occluded_nondeviant"
    )
    assert (
        dot_side.resolve_model_condition("occluded_deviant", "observed")
        == "occluded_deviant"
    )
    with pytest.raises(ValueError, match="valid only"):
        dot_side.resolve_model_condition("always_visible", "predicted")

    def fake_loadmat(_path, simplify_cells):
        assert simplify_cells is True
        return {
            "Cfg": {"rectSize": [10.0, 8.0]},
            "xySeqs": [
                {
                    "condition_label": condition,
                    "sequence": identity,
                    "xy": np.column_stack(
                        [
                            np.full(4, x_value, dtype=float),
                            np.zeros(4, dtype=float),
                        ]
                    ),
                }
                for condition, x_value in (
                    ("occluded_nondeviant", 2.0),
                    ("occluded_deviant", 8.0),
                )
                for identity in (1, 2)
            ],
        }

    monkeypatch.setattr(dot_side.sio, "loadmat", fake_loadmat)
    predicted_labels, _x, _center = dot_side.build_dot_side_labels(
        model_mat=Path("unused.mat"),
        condition=dot_side.resolve_model_condition(
            "occluded_deviant", "predicted"
        ),
        trial_identities=np.array([1, 2]),
        identity_values=[1, 2],
        n_times_target=4,
    )
    observed_labels, _x, _center = dot_side.build_dot_side_labels(
        model_mat=Path("unused.mat"),
        condition=dot_side.resolve_model_condition(
            "occluded_deviant", "observed"
        ),
        trial_identities=np.array([1, 2]),
        identity_values=[1, 2],
        n_times_target=4,
    )
    assert np.all(predicted_labels == 1)
    assert np.all(observed_labels == 0)


def test_output_paths_include_condition_and_model_source(tmp_path):
    """Observed and predicted deviant bundles must never share filenames."""
    observed = dot_side.sanity_output_paths(
        tmp_path / "observed", "occluded_deviant", "observed"
    )
    predicted = dot_side.sanity_output_paths(
        tmp_path / "predicted", "occluded_deviant", "predicted"
    )
    assert "occluded_deviant_observed" in observed["summary_json"].name
    assert "occluded_deviant_predicted" in predicted["summary_json"].name
    assert observed["summary_json"] != predicted["summary_json"]


def test_too_few_class_samples_are_explicit_nan_windows():
    """Windows with too few samples for stratified CV are skipped explicitly."""
    rng = np.random.default_rng(1)
    data = rng.normal(size=(4, 3, 6))
    labels = np.tile(np.array([[0], [0], [1], [1]]), (1, 6))
    accuracy_unbalanced, accuracy_balanced, rows = dot_side.decode_sliding_windows(
        data,
        labels,
        times=np.arange(6, dtype=float) / 10.0,
        n_splits=3,
        window_samples=3,
        random_seed=11,
        sfreq=10.0,
    )

    assert np.all(np.isnan(accuracy_unbalanced))
    assert np.all(np.isnan(accuracy_balanced))
    assert {
        row["skip_reason"] for row in rows
    } == {"too_few_class_samples_for_3_fold_stratified_cv"}


def test_training_testing_diagnostics_assesses_skipped_windows():
    """Diagnostics report class balance and whether CV was fully defined."""
    rng = np.random.default_rng(1)
    data = rng.normal(size=(8, 3, 6))
    labels = np.tile(np.array([[0], [0], [1], [1], [1], [1], [1], [1]]), (1, 6))
    acc_unbalanced, acc_balanced, rows = dot_side.decode_sliding_windows(
        data,
        labels,
        times=np.arange(6, dtype=float) / 10.0,
        n_splits=3,
        window_samples=3,
        random_seed=11,
        sfreq=10.0,
    )
    diagnostics = dot_side.summarize_training_testing_diagnostics(
        window_rows=rows,
        averaging_metadata={
            "mode": "all_trials_per_identity",
            "requested_value": "all_trials_per_identity",
            "n_source_trials": 40,
        },
        n_decoder_trials=8,
        n_decoder_identities=8,
        n_splits=3,
        window_samples=3,
        sfreq=10.0,
        accuracy_unbalanced=acc_unbalanced,
        accuracy_balanced=acc_balanced,
    )

    assert diagnostics["assessment"]["status"] == "not_well_defined"
    assert diagnostics["class_balance"]["n_skipped_windows"] == 6
    assert diagnostics["class_balance"]["min_minority_class_count_all_windows"] == 2
    assert diagnostics["sample_structure"]["averaging_mode"] == "all_trials_per_identity"
    html = dot_side.training_testing_diagnostics_html(
        {"training_testing_diagnostics": diagnostics}
    )
    assert "Regularized linear discriminant analysis" in html
    assert "too_few_class_samples_for_3_fold_stratified_cv" in html
