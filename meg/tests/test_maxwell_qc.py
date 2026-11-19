"""Test channel evidence plots, using only small synthetic data.

Usage: MPLBACKEND=Agg python -B -m pytest -p no:cacheprovider meg/tests/test_maxwell_qc.py -q

Check per-channel thresholds, missing scores, native sample coordinates,
reproducible random selection, bounded reads, highlights and stale-report guards.
No cohort data, SSS solver or scientific decisions are modified.
"""

import copy
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import maxwell_qc as qc
import render_stage03_mne_reports as renderer


def scores_example():
    """Use different limits to catch accidental global-threshold comparisons."""
    return {"ch_names": np.array(["MEG0111", "MEG0121", "MEG0131"]),
            "ch_types": np.array(["mag", "mag", "mag"]), "bins": np.array([[0., 5.], [8., 13.]]),
            "scores_noisy": np.array([[8., np.nan], [8., 10.], [np.nan, np.nan]]),
            "limits_noisy": np.array([[7.], [9.], [np.nan]]),
            "scores_flat": np.array([[.5e-15, 2e-15], [.5e-15, 2e-15], [np.nan, np.nan]]),
            "limits_flat": np.array([[1e-15], [.1e-15], [np.nan]])}


def test_each_cell_uses_its_own_limit_and_nan_is_unscored():
    """Transient crossings and header-excluded NaNs must not change labels."""
    scores = scores_example()
    _, values, limits, crossed = qc.score_values(scores, "mag", "noisy")
    np.testing.assert_array_equal(crossed, [[True, False], [False, True], [False, False]])
    assert np.isnan(values[2]).all() and np.isnan(limits[2]).all()
    _, _, _, crossed = qc.score_values(scores, "mag", "flat")
    np.testing.assert_array_equal(crossed, [[True, False], [False, False], [False, False]])


def test_random_window_is_repeatable_native_and_short_run_safe():
    """Time is recording-relative; native sample coordinates include first_samp."""
    raw = mne.io.RawArray(np.zeros((1, 6000)), mne.create_info(["MEG0111"], 100., ["mag"]), first_samp=1397000, verbose="ERROR")
    first = qc.choose_window(raw, "source.fif", 20., 97)
    assert first == qc.choose_window(raw, "source.fif", 20., 97)
    assert first != qc.choose_window(raw, "other_source.fif", 20., 97)
    assert first["stop_index"] - first["start_index"] == 2000
    assert first["native_start_sample"] == 1397000 + first["start_index"]
    assert first["stop_sec"] - first["start_sec"] == pytest.approx(20.)
    short = raw.copy().crop(tmax=2.)
    window = qc.choose_window(short, "source.fif", 20.)
    assert window["start_index"] == 0 and window["stop_index"] == short.n_times
    with pytest.raises(ValueError):
        qc.choose_window(raw, "source.fif", -1)


def test_heatmap_highlights_detected_and_applied_only_without_mutating_scores(tmp_path, monkeypatch):
    """Selected NaN rows remain visible and detected/applied sets stay separate."""
    captured = []
    monkeypatch.setattr(plt.Figure, "savefig", lambda figure, *a, **kw: captured.append(figure))
    scores = scores_example(); before = copy.deepcopy(scores)
    plots = qc.save_heatmaps(scores=scores, source_filename="source.fif",
        detected_bads=["MEG0131"], applied_bads=["MEG0121", "MEG0131"],
        prefix=tmp_path / "plots/source", raw_root=tmp_path / "raw")
    assert len(plots) == 2  # noisy/flat magnetometers; absent grads are skipped.
    for figure in captured:
        for axis in figure.axes[:2]:
            labels = {label.get_text(): label for label in axis.get_yticklabels()}
            assert labels["MEG0131 [D]"].get_color() == qc.DETECTED_COLOR
            assert labels["MEG0121 [A]"].get_color() == qc.APPLIED_COLOR
            assert len(axis.patches) == 2
    for key in scores:
        np.testing.assert_array_equal(scores[key], before[key])


def test_traces_read_only_twenty_seconds_and_keep_physical_amplitude(tmp_path, monkeypatch):
    """All sensors remain visible, including header bads, without normalization."""
    data = np.arange(2 * 3000, dtype=float).reshape(2, -1) * 1e-15
    raw = mne.io.RawArray(data, mne.create_info(["MEG0111", "MEG0121"], 100., ["mag", "mag"]), first_samp=4000, verbose="ERROR")
    raw.info["bads"] = ["MEG0111"]
    reads = []; captured = []; original = raw.get_data
    def bounded_read(*args, **kwargs):
        reads.append((kwargs["start"], kwargs["stop"]))
        return original(*args, **kwargs)
    monkeypatch.setattr(raw, "get_data", bounded_read)
    monkeypatch.setattr(plt.Figure, "savefig", lambda figure, *a, **kw: captured.append(figure))
    _, window = qc.save_traces(raw=raw, source_filename="source.fif", detected_bads=["MEG0111"],
        applied_bads=["MEG0111"], prefix=tmp_path / "plots/source", raw_root=tmp_path / "raw")
    assert reads == [(window["start_index"], window["stop_index"])]
    assert reads[0][1] - reads[0][0] == 2000
    axis = captured[0].axes[0]
    assert len(axis.lines) == 2
    flagged = next(line for line in axis.lines if line.get_label() == "MEG0111")
    segment = data[0, window["start_index"]:window["stop_index"]]
    np.testing.assert_allclose(flagged.get_ydata(), (segment - np.median(segment)) * 1e15)
    np.testing.assert_array_equal(raw._data, data)
    assert raw.info["bads"] == ["MEG0111"]


def test_renderer_replaces_detection_average_and_rejects_stale_highlights(tmp_path):
    """Presence of new metadata selects channel plots, never the old PSD."""
    for name in ["new_heatmap.png", "new_traces.png", "03_A_MaxFilt_review_source_before_after_psd.png"]:
        (tmp_path / name).write_bytes(b"present")
    run = {"source_filename": "source.fif", "acquisition_run_label": 1,
           "detected_bad_channels": ["MEG0111"], "applied_bad_channels": ["MEG0111", "MEG0121"]}
    metadata = {"runs": [{**run, "window": {"start_sec": 1., "stop_sec": 21.},
        "plots": [{"filename": name, "title": name, "section": "Channel evidence"}
                  for name in ["new_heatmap.png", "new_traces.png"]]}]}
    _, figures = renderer._maxwell_content({"runs": [run]}, tmp_path, "review", metadata)
    assert len(figures) == 3 and not any("detection_psd" in str(path) for path, _, _ in figures)
    metadata["runs"][0]["detected_bad_channels"] = []
    with pytest.raises(ValueError, match="Stale QC"):
        renderer._maxwell_content({"runs": [run]}, tmp_path, "review", metadata)


def test_plot_only_refresh_preserves_scientific_evidence_and_checks_source(tmp_path):
    """Run the real plot/report refresh and refuse changed acquisition inputs."""
    import json
    import os
    raw_root = tmp_path / "raw"
    source = raw_root / "sub04/meg/source_raw.fif"
    source.parent.mkdir(parents=True)
    scores = scores_example()
    raw = mne.io.RawArray(np.random.default_rng(9).normal(size=(3, 3000)) * 1e-13,
        mne.create_info(list(scores["ch_names"]), 100., ["mag"] * 3), first_samp=4000, verbose="ERROR")
    raw.save(source, overwrite=True, verbose="ERROR")
    derivatives = tmp_path / "derivatives"
    reports = derivatives / "reports/sub-04/ses-01"; reports.mkdir(parents=True)
    figures = derivatives / "preprocessing/sub-04/ses-01/figures"; figures.mkdir(parents=True)
    sss_path = derivatives / "preprocessing/sub-04/ses-01/sss/source_raw_sss.fif"
    sss_path.parent.mkdir()
    raw.copy().save(sss_path, overwrite=True, verbose="ERROR")
    stage = reports / "03_A_MaxFilt_review_report.json"
    stage.write_text(json.dumps({"subject_bids_id": "sub-04", "processing_pass": "review",
        "status": "completed", "scientific_acceptance": "separate_explicit_stage4_review",
        "source_files": [{"relative_path": str(source.relative_to(raw_root)), "stamp": qc.file_stamp(source)}],
        "evidence_files": [{"relative_path": str(sss_path.relative_to(derivatives)), "stamp": qc.file_stamp(sss_path)}],
        "runs": [{"source_filename": source.name, "acquisition_run_label": 1,
                  "output_path": str(sss_path),
                  "detected_bad_channels": ["MEG0131"], "applied_bad_channels": ["MEG0121", "MEG0131"]}]}))
    score_path = reports / "03_A_MaxFilt_review_source_raw_bad_channel_scores.npz"
    np.savez_compressed(score_path, **scores)
    stage_before = stage.read_bytes();score_before = score_path.read_bytes();raw_stamp = qc.file_stamp(source)
    output = qc.refresh_subject(derivatives_root=derivatives, raw_root=raw_root, subject_id="04")
    html = output.read_text()
    assert "data:image/png;base64," in html
    assert "Bad-channel raw traces" in html and "Bad-channel score heatmaps" in html
    assert "legacy detection-average PSD" not in html
    assert stage.read_bytes() == stage_before and score_path.read_bytes() == score_before
    assert qc.file_stamp(source) == raw_stamp
    assert html.index('alt="Run 01: Raw PSD:') < html.index('alt="Run 01: Noisy mag:')
    assert html.index('alt="Run 01: Retained-channel PSD:') > html.index('alt="Run 01: Raw channels:')
    metadata = reports / "03_A_MaxFilt_review_plot_metadata.json"
    before = metadata.read_bytes()
    # A changed SSS file also makes the plotted evidence stale, without writing
    # data, decisions or new metadata as a side effect of that validation error.
    sss_stat = sss_path.stat()
    os.utime(sss_path, ns=(sss_stat.st_atime_ns, sss_stat.st_mtime_ns + 1_000_000))
    with pytest.raises(RuntimeError, match="SSS source changed"):
        qc.refresh_subject(derivatives_root=derivatives, raw_root=raw_root, subject_id="04", overwrite=True)
    os.utime(sss_path, ns=(sss_stat.st_atime_ns, sss_stat.st_mtime_ns))
    stat = source.stat();os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
    with pytest.raises(RuntimeError, match="Raw source changed"):
        qc.refresh_subject(derivatives_root=derivatives, raw_root=raw_root, subject_id="04", overwrite=True)
    assert metadata.read_bytes() == before


def test_per_channel_psds_match_windows_exclude_bads_and_show_expected_gain(tmp_path, monkeypatch):
    """The final panels retain separate good sensors and comparable dB scales."""
    times = np.arange(30000) / 200.
    signals = np.array([np.sin(2*np.pi*8*times), 3*np.sin(2*np.pi*12*times),
                        50*np.sin(2*np.pi*20*times)]) * 1e-13
    # A broadband floor keeps this gain check away from numerical-zero bins
    # outside a mathematically pure sinusoid's support.
    signals += np.random.default_rng(8).normal(scale=1e-15, size=signals.shape)
    raw = mne.io.RawArray(signals, mne.create_info(
        ["MEG0111", "MEG0121", "MEG0131"], 200., ["mag"]*3), first_samp=4000, verbose="ERROR")
    raw.info["bads"] = ["MEG0131"]
    sss = mne.io.RawArray(signals * .1, raw.info.copy(), first_samp=4000, verbose="ERROR")
    figures = []
    monkeypatch.setattr(plt.Figure, "savefig", lambda figure, *a, **kw: figures.append(figure))
    before = raw.get_data().copy();sss_before = sss.get_data().copy()
    plots, settings = qc.save_psds(raw=raw, sss=sss, source_filename="source.fif",
        detected_bads=["MEG0131"], applied_bads=["MEG0121", "MEG0131"],
        prefix=tmp_path / "plots/source", raw_root=tmp_path / "raw")
    assert settings["duration_sec"] == 120.
    assert settings["stop_index"] == 24000
    assert settings["comparison_channels"] == ["MEG0111"]
    assert settings["raw_channels"] == raw.ch_names
    assert len(plots) == 2
    overview = figures[0].axes[0]
    assert overview.get_xlabel() == "Frequency (Hz)"
    assert {"MEG0121 [A]", "MEG0131 [D]"}.issubset(overview.get_legend_handles_labels()[1])
    left, right = figures[1].axes[:2]
    assert left.get_ylim() == right.get_ylim()
    assert left.get_xlabel() == right.get_xlabel() == "Frequency (Hz)"
    # Each retained sensor is a distinct curve; constant .1 amplitude gain
    # produces a -20 dB power difference on the same native Welch windows.
    left_curves = [line for line in left.lines if len(line.get_xdata()) > 2]
    right_curves = [line for line in right.lines if len(line.get_xdata()) > 2]
    assert len(left_curves) == len(right_curves) == 1
    np.testing.assert_allclose(right_curves[0].get_ydata() - left_curves[0].get_ydata(), -20., atol=.01)
    np.testing.assert_array_equal(raw.get_data(), before)
    np.testing.assert_array_equal(sss.get_data(), sss_before)
    assert raw.info["bads"] == sss.info["bads"] == ["MEG0131"]
