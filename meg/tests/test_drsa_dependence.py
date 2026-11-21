"""Regressions for the descriptive RDM correlation contract, lag and SEM.

Run with python -m pytest meg/tests/test_drsa_dependence.py from the repo root.
All synthetic result files are written under pytest's temporary directory.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

CODE = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE))
from drsa_dependence import DEVIANT, diagnostics_root, read_bundle, read_condition
from continuous_cleaning import file_hash
from plot_drsa_dependence import mean_and_support, summarize


def test_deviant_prefix_and_lag_orientation(tmp_path):
    rng = np.random.default_rng(19)
    axis = np.arange(60) / 120
    observed = rng.standard_normal((30, axis.size))
    predicted = observed.copy()
    predicted[:, 15:] = observed[:, 12:-3]
    observed_auto, predicted_auto = np.corrcoef(observed.T), np.corrcoef(predicted.T)
    cross = np.corrcoef(observed.T, predicted.T)[:axis.size, axis.size:]
    matrices = {
        "observed_self": observed_auto,
        "predicted_self": predicted_auto,
        "observed_vs_predicted": cross,
        "neural_occipital": observed_auto,
    }
    identities = list(range(1, 9))
    data = dict(
        matrices=matrices, condition=DEVIANT, time_axis_sec=axis, identities=identities,
        shared_prefix_frames=15,
    )
    path = tmp_path / "condition.mat"
    savemat(path, data)
    handoff = {"time_axis_sec": axis, "identity_values": identities}
    validated = read_condition(path, DEVIANT, ["occipital"], handoff)
    np.testing.assert_allclose(validated["matrices"]["observed_vs_predicted"], cross)

    late = axis >= axis[16]
    summaries = summarize([cross, cross], axis, {"late": late})
    values = summaries["curves"][0]
    assert summaries["lag_sec"][np.nanargmax(values)] == pytest.approx(3 / 120)
    assert summaries["curves_n"][0, np.nanargmax(values)] == 2

    # Transposing a signed comparison should flip the measured lag. A changed
    # prefix must be caught before a participant result is published.
    reversed_summary = summarize([cross.T], axis, {"late": late})
    assert reversed_summary["lag_sec"][np.nanargmax(reversed_summary["curves"][0])] == pytest.approx(-3 / 120)
    matrices["observed_vs_predicted"][0, 0] = 0
    savemat(path, data)
    with pytest.raises(ValueError, match="Shared prefix"):
        read_condition(path, DEVIANT, ["occipital"], handoff)


def test_isolated_output_and_artifact_verification(tmp_path):
    raw = tmp_path / "raw"
    derivative = tmp_path / "derivatives"
    paths = {
        "analysis_dir": derivative / "analysis" / "sub-04" / "ses-01",
        "reports_dir": derivative / "reports" / "sub-04" / "ses-01",
    }
    cfg = {"resolved_paths": {"data_root": str(raw)}}
    with pytest.raises(ValueError, match="outside the ordinary"):
        diagnostics_root(derivative / "analysis" / "group", cfg, paths)
    directory = diagnostics_root(derivative / "diagnostics", cfg, paths)
    directory.mkdir(parents=True)
    handoff = directory / "handoff.json"
    handoff.write_text("{}")
    result = directory / "always_visible.mat"
    result.write_bytes(b"test")
    report = {
        "status": "completed", "conditions": ["always_visible"], "contract_sha256": "settings",
        "sha256": {path.name: file_hash(path) for path in (handoff, result)},
    }
    (directory / "report.json").write_text(json.dumps(report))
    assert read_bundle(directory, contract="settings") == report
    with pytest.raises(ValueError, match="Different diagnostic inputs"):
        read_bundle(directory, contract="other")
    result.write_bytes(b"changed")
    with pytest.raises(ValueError, match="artifact changed"):
        read_bundle(directory, contract="settings")


def test_group_support_uses_finite_participants_only():
    mean, n = mean_and_support([[0, np.nan, np.nan], [1, 0.6, np.nan]])
    np.testing.assert_allclose(mean, [0.5, 0.6, np.nan])
    np.testing.assert_array_equal(n, [2, 1, 0])


def test_lag_window_and_sem_follow_finite_subject_support():
    axis = np.arange(70) / 25
    first = np.full((axis.size, axis.size), 0.2)
    second = np.full_like(first, 0.6)
    np.fill_diagonal(second[:, 1:], np.nan)
    for matrix in (first, second):
        np.fill_diagonal(matrix[:, 2:], np.nan)

    summaries = summarize([first, second], axis, {"full": axis >= 0}, lag_window_sec=2.0)
    np.testing.assert_allclose(summaries["lag_sec"][[0, -1]], [-2.0, 2.0])
    np.testing.assert_allclose(summaries["surface_lag_sec"], summaries["lag_sec"])
    zero = summaries["lag_sec"].size // 2
    bins = slice(zero, zero + 3)
    np.testing.assert_array_equal(summaries["curves_n"][0, bins], [2, 1, 0])
    np.testing.assert_allclose(summaries["curves"][0, bins], [0.4, 0.2, np.nan])
    np.testing.assert_allclose(summaries["curves_sem"][0, bins], [0.2, np.nan, np.nan])
