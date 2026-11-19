"""Focused EDF guidance tests, entirely on temporary/synthetic evidence.

Usage: python -B -m pytest -p no:cacheprovider meg/tests/test_edf_blink_candidates.py -q

These tests cover trial gaps, ambiguity, long/binocular events, native annotation
coordinates, sidecar reuse, explicit promotion and the disposable browser. No
cohort recording is read, processed, or written. SDK decoding is also checked
separately against eyelinkio's bundled generic fixtures, never DAD recordings.
"""

import copy
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import edf_blink_candidates as edf
from continuous_cleaning import good_stretches


def raw_example(date=False):
    """Make a native stream with nonzero first_samp and an unchanged MISC trace."""
    raw = mne.io.RawArray(np.zeros((2, 1200)), mne.create_info(["MEG0111", "MISC001"], 100., ["mag", "misc"]),
                          first_samp=1000, verbose="ERROR")
    if date:
        raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    return raw


def records_example():
    """Two measured EDF trials separated by a real, large recording-off gap."""
    messages = []
    for number, start in [(1, 100.), (2, 200.)]:
        messages.extend([{"time": start, "text": f"TRIALID {number}"},
                         {"time": start + .2, "text": f"TRIGGER 71 trial {number} frame 1 label trial_condition_onset"},
                         {"time": start + .4, "text": f"TRIGGER {number} trial {number} frame 2 label identity"}])
    return {"messages": messages, "blocks": [{"start": 100.1, "stop": 101.1, "sfreq": 100.},
                                               {"start": 200.1, "stop": 201.1, "sfreq": 100.}],
            "blinks": [{"start": 100.3, "stop": 100.5, "eye": 0},
                       {"start": 100.31, "stop": 100.51, "eye": 1},
                       {"start": 200.3, "stop": 201.0, "eye": 1},
                       {"start": 101.0, "stop": 200.2, "eye": 0}]}


def event_rows():
    """Accepted canonical trigger identities; elapsed EDF gaps are irrelevant."""
    rows = []
    for trial, sample in [(1, 1200), (2, 1800)]:
        for code, frame, delta in [(71, 1, 0), (trial, 2, 20)]:
            rows.append({"source_fif": "sub17_run01_test.fif", "trigger_code": code,
                         "raw_sample_frame": sample + delta, "behavioral_frame": frame,
                         "behavioral_trial": trial, "behavioral_block": 1,
                         "reported_label": "trial_condition_onset" if code == 71 else "identity"})
    return rows


def test_trial_local_mapping_never_bridges_gap_or_applies_a_fixed_offset():
    records = records_example(); eyes = edf.edf_trials(records)
    matches = edf.match_trials(eyes, edf.meg_trials(event_rows()), [1, 60], 100., .02)
    assert [r["status"] for r in matches] == ["matched", "matched"]
    rows, segments = edf.map_blinks(records, eyes, matches, "sub17_run01_test.fif", 100., 1000, 1200, .5)
    assert rows[0]["sample_start_raw"] == 1210
    assert rows[2]["sample_start_raw"] == 1810  # Independent trial anchor, not 11210.
    assert rows[2]["long_event"] and rows[2]["status"] == "candidate"
    assert rows[3]["status"] == "outside_matched_recording"
    assert len(segments) == 2


def test_mismatched_duplicate_markers_and_unknown_stops_are_flagged():
    records = records_example(); eyes = edf.edf_trials(records); megs = edf.meg_trials(event_rows())
    megs[0]["markers"].append(copy.deepcopy(megs[0]["markers"][1]))
    matches = edf.match_trials(eyes, megs, [1, 60], 100., .02)
    assert matches[0]["status"] in ("mismatched_markers", "ambiguous_or_unmatched")
    matches = edf.match_trials(eyes, edf.meg_trials(event_rows()), [1, 60], 100., .02)
    records["blocks"][1]["stop"] = None
    rows, segments = edf.map_blinks(records, eyes, matches, "sub17_run01_test.fif", 100., 1000, 1200, .5)
    assert len(segments) == 1 and segments[0]["edf_trial_ordinal"] == 0
    assert rows[2]["status"] != "candidate"
    assert matches[1]["status"] == "unknown_or_ambiguous_recording_bounds"
    # Timing disagreement is flagged, never repaired by fitting a lag.
    megs = edf.meg_trials(event_rows()); megs[0]["markers"][1]["sample"] += 10
    assert edf.match_trials(eyes, megs, [1, 60], 100., .02)[0]["status"] == "mismatched_markers"


def test_missing_trial_does_not_shift_later_unique_matches():
    eyes = []; megs = []
    for i in range(1, 9):
        markers = [{"code": 71, "frame": 1, "time": i * 10., "trial_id": i},
                   {"code": i, "frame": 2, "time": i * 10. + .2, "trial_id": i}]
        eyes.append({"trial_id": i, "markers": markers, "anchor": markers[0]})
        m = [{**r, "sample": int(r["time"] * 100), "unresolved": False} for r in markers]
        megs.append({"trial_id": i, "markers": m, "anchor": m[0], "source": "source.fif"})
    del megs[2]
    matches = edf.match_trials(eyes, megs, [1, 60], 100., .02)
    assert matches[2]["status"] == "ambiguous_or_unmatched"
    assert matches[6]["meg_anchor_sample"] == 7000
    repeated = [copy.deepcopy(eyes[0]) for _ in range(5)]
    assert not any(r["status"] == "matched" for r in edf.match_trials(repeated, megs, [1, 60], 100., .02))


@pytest.mark.parametrize("date", [False, True])
def test_candidates_are_non_excluding_and_promotion_requires_explicit_review(date):
    raw = raw_example(date)
    raw.set_annotations(mne.Annotations([2., 8.], [.2, .7], list(edf.CANDIDATE_LABELS)))
    assert good_stretches(raw) == [(0, raw.n_times)]
    events = np.array([[1200, 0, 1], [1800, 0, 1]])
    epochs = mne.Epochs(raw, events, {"trial": 1}, tmin=-.1, tmax=.3, baseline=None,
                       preload=True, reject_by_annotation=True, verbose="ERROR")
    assert len(epochs) == 2
    with pytest.raises(ValueError, match="Explicit"):
        edf.promote_candidates_after_review(raw, confirmed=False)
    # Simulate the browser: remove one candidate, move/resize the retained one.
    raw.annotations.delete(1)
    raw.annotations.onset[0] += .1; raw.annotations.duration[0] = .3
    edf.promote_candidates_after_review(raw, confirmed=True)
    assert raw.annotations.description.tolist() == ["BAD_blink"]
    assert raw.annotations.onset[0] == pytest.approx(raw.first_time + 2.1)
    epochs = mne.Epochs(raw, events, {"trial": 1}, tmin=-.1, tmax=.3, baseline=None,
                       preload=True, reject_by_annotation=True, verbose="ERROR")
    assert epochs.selection.tolist() == [1]


@pytest.mark.parametrize("date", [False, True])
def test_browser_samples_keep_nan_gaps_and_original_channels(date):
    raw = raw_example(date); before = raw.get_data().copy()
    records = records_example(); eyes = edf.edf_trials(records)
    matches = edf.match_trials(eyes, edf.meg_trials(event_rows()), [1, 60], 100., .02)
    _, segments = edf.map_blinks(records, eyes, matches, "sub17_run01_test.fif", 100., 1000, 1200, .5)
    for seg in segments:
        seg["edf_relative_path"] = "eye.edf"
    times = np.r_[np.arange(100.1, 101.1, .01), np.arange(200.1, 201.1, .01)]
    samples = np.column_stack([times, times, times * 2, times * 3, times * 4])
    # A missing measurement inside a trial must remain a gap as well.
    samples[20, 1:] = np.nan
    proposal = {"segments": segments}
    data, _ = edf.display_samples(raw, proposal, {"eye.edf": {"samples": samples}})
    assert np.isnan(data[:4, 400:700]).all()
    assert not data[-1, 400:700].any()
    assert np.isnan(data[:4, 210]).all()
    assert np.isfinite(data[:4, 200]).all() and data[-1, 200] == 1.
    np.testing.assert_array_equal(raw.get_data(), before)


@pytest.mark.parametrize("date", [False, True])
def test_sidecars_reuse_exact_evidence_and_preserve_native_onsets(tmp_path, monkeypatch, date):
    from preprocessing_stage_utils import AcquisitionRunRecord
    raw = raw_example(date); root = tmp_path / "raw"; root.mkdir()
    eye = root / "MB170100.edf"; eye.write_bytes(b"synthetic EDF input")
    sss = tmp_path / "sss.fif"; sss.write_bytes(b"synthetic SSS identity")
    derivative = tmp_path / "derivatives/preprocessing/sub-17/ses-01"
    cfg = {"subject": {"id": "17"}, "subject_profile": {"identity_range": [1, 60]},
           "resolved_paths": {"data_root": str(root), "derivatives_preprocessing_dir": str(derivative)},
           "artifact_annotation": {"edf_blink_candidates": {"enabled": True, "long_event_sec": .5, "marker_tolerance_sec": .02}}}
    record = AcquisitionRunRecord("raw", "sub17_run01_test.fif", "sub17_run01_test", 1, "raw")
    context = {"events": event_rows(), "events_sha256": "test", "manifest": {
        "sources": {"eyetracker_edf_files": [{"source_filename": eye.name, "raw_relative_path": eye.name}]},
        "trigger_derived_associations": {"meg_acquisition_run_to_behavioral_segments": [{
            "source_fif": record.source_filename, "behavioral_segments": [{"behavioral_block": 1}]}]}}}
    calls = []
    monkeypatch.setattr(edf, "read_edf_records", lambda path: calls.append(path) or records_example())
    first = edf.prepare_candidates(cfg, record, raw, sss, context, {})
    assert first["candidate_count"] == 2  # Overlapping binocular records merge.
    assert first["long_candidate_count"] == 1
    ann = mne.read_annotations(first["paths"]["annotations"])
    np.testing.assert_allclose(ann.onset, [12.1, 18.1], atol=1e-6)
    assert ann.orig_time == raw.annotations.orig_time
    assert not any(str(d).lower().startswith("bad") for d in ann.description)
    second = edf.prepare_candidates(cfg, record, raw, sss, context, {})
    assert first["fingerprint"] == second["fingerprint"] and len(calls) == 1
    eye.write_bytes(b"changed EDF identity")
    third = edf.prepare_candidates(cfg, record, raw, sss, context, {})
    assert third["fingerprint"] != first["fingerprint"] and len(calls) == 2


def test_matplotlib_browser_accepts_nan_guidance(tmp_path, monkeypatch):
    """Smoke-test the public MNE plot API actually used during local review."""
    raw = raw_example(); times = np.arange(100.1, 101.1, .01)
    sample_rows = np.column_stack([times, times, times, times, times])
    proposal = {"segments": [{"edf_relative_path": "eye.edf", "start": 100.1, "stop": 101.1,
        "edf_anchor_sec": 100.2, "meg_anchor_sample": 1200, "sfreq": 100.}]}
    monkeypatch.setattr(edf, "read_edf_records", lambda *a, **k: {"samples": sample_rows})
    before = raw.ch_names.copy()
    metadata = edf.add_edf_display(raw, {"resolved_paths": {"data_root": str(tmp_path)}}, proposal)
    assert raw.ch_names[:2] == before and raw.get_channel_types()[:2] == ["mag", "misc"]
    assert "EDF_no_guidance" in raw.annotations.description
    from manual_eye_annotations import plot_review_browser
    original_plot = raw.plot
    monkeypatch.setattr(raw, "plot", lambda **kwargs: original_plot(**{**kwargs, "show": False, "block": False}))
    with mne.viz.use_browser_backend("matplotlib"):
        figure = plot_review_browser(raw, ["MISC001"], metadata)
    assert all(name.startswith("EDF_") for name in raw.ch_names[:5])
    assert raw.ch_names[5] == "MISC001"
    figure.savefig(tmp_path / "edf_browser.png")
    import matplotlib.pyplot as plt
    plt.close(figure)
