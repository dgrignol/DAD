"""Synthetic checks for Stage 4B cached ICA and automatic blink proposals.

Usage: python -B -m pytest -p no:cacheprovider meg/tests/test_blink_annotation.py -q

Only temporary synthetic recordings are read or written. Tests cover all-frontal
mag/grad selection, the ten-component cap, native acquisition coordinates,
BAD-isolated preparation, cache lineage, peak polarity, duplicate intervals,
unscored gaps, and preservation of manual edits/deletions on re-detection.
These checks validate implementation, not detector accuracy on DAD recordings.
"""

from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import blink_annotation as blink
from continuous_cleaning import annotation_rows, coordinates, file_hash, file_stamp
from preprocessing_stage_utils import AcquisitionRunRecord


# -----------------------------------------------------------------------------
# Deterministic recordings with both frontal sensor types and auxiliary channels
# -----------------------------------------------------------------------------

def preparation_example(tmp_path, monkeypatch, *, date=True):
    """Make a source with eleven good frontal channels and a known bad interval."""
    names = [f"MEG{station}{axis}" for station in ("012", "031", "081", "091") for axis in ("1", "2", "3")]
    kinds = ["mag" if name.endswith("1") else "grad" for name in names]
    rng = np.random.default_rng(17)
    latent = rng.laplace(size=(12, 8000))
    data = rng.normal(size=(12, 12)) @ latent
    data *= np.array([1e-12 if kind == "mag" else 1e-10 for kind in kinds])[:, None]
    raw = mne.io.RawArray(np.vstack([data, np.zeros((2, 8000))]),
        mne.create_info(names + ["MISC001", "STI101"], 200., kinds + ["misc", "stim"]),
        first_samp=1397000, verbose="ERROR")
    if date:
        raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    raw.info["bads"] = ["MEG0122"]
    raw.set_annotations(mne.Annotations([10., 20.], [1., .3], ["BAD_sensor_fault", "BAD_blink"]))
    source = tmp_path / "native_sss.fif"; source.write_bytes(b"source_stat_fixture")
    cfg = {"ica": {"n_components": None, "fit_l_freq_hz": 1., "fit_h_freq_hz": 40.,
        "fit_resample_sfreq_hz": 100., "fit_resample_edge_sec": .1, "fit_min_good_sec": 1.},
        "artifact_annotation": {"blink_auto": dict(blink.DEFAULTS)},
        "resolved_paths": {"data_root": str(tmp_path / "raw"),
            "derivatives_preprocessing_dir": str(tmp_path / "preprocessing"),
            "derivatives_reports_dir": str(tmp_path / "reports")}}
    record = AcquisitionRunRecord(str(source), "sub99_run02_partial_test.fif", "sub99_run02_partial_test", 2, "source")
    monkeypatch.setattr(blink, "load_aid_input", lambda *a: (blink.aid_mask(raw), source))
    return cfg, raw, source, record


@pytest.mark.parametrize("date", [False, True])
def test_cache_roundtrip_uses_all_good_frontal_types_and_caps_only_the_aid(tmp_path, monkeypatch, date):
    cfg, raw, source, record = preparation_example(tmp_path, monkeypatch, date=date)
    original = coordinates(raw), annotation_rows(raw), raw.get_data().copy()
    result = blink.prepare_aid(cfg, record)
    assert result["n_components_fitted"] == 10
    assert len(result["identity"]["channels"]) == 11
    assert "MEG0122" not in result["identity"]["channels"]
    assert {kind for name, kind in zip(raw.ch_names, raw.get_channel_types()) if name in result["identity"]["channels"]} == {"mag", "grad"}
    assert "MISC001" not in result["identity"]["channels"]
    assert cfg["ica"]["n_components"] is None
    state = blink.load_cached_aid(cfg, record)
    sources = state["native_sources"]
    for key in ("first_samp", "n_times", "sfreq", "meas_date"):
        assert coordinates(sources)[key] == coordinates(raw)[key]
    values = sources.get_data()
    assert np.isnan(values[:, :330]).all()  # Full 1 Hz filter half-support at native 200 Hz.
    assert np.isnan(values[:, 2000:2200]).all()  # Known severe interval is never projected.
    assert np.isfinite(values[:, 4000:4060]).all()  # Ocular examples remain available for this aid.
    assert original[0] == coordinates(raw) and original[1] == annotation_rows(raw)
    np.testing.assert_array_equal(original[2], raw.get_data())
    assert blink.prepare_aid(cfg, record)["execution_status"] == "reused"
    cfg["artifact_annotation"]["blink_auto"]["threshold_std"] = 7.
    assert blink.load_cached_aid(cfg, record)["metadata"]["fingerprint"] == result["fingerprint"]
    # New detector tuning and reviewed ocular labels do not alter the legacy
    # fit identity or require preparation of an otherwise compatible cache.
    cfg["artifact_annotation"]["eye_movement_auto"] = {"threshold_std": 3.}
    raw.annotations.append(raw.first_time + 25., .5, blink.EYE_AUTO_LABEL)
    assert blink.load_cached_aid(cfg, record)["metadata"]["fingerprint"] == result["fingerprint"]
    cfg["ica"]["fit_h_freq_hz"] = 35.
    with pytest.raises(RuntimeError, match="stale"):
        blink.load_cached_aid(cfg, record)


def test_cache_rejects_changed_source_mask_or_cache_file(tmp_path, monkeypatch):
    cfg, raw, source, record = preparation_example(tmp_path, monkeypatch)
    blink.prepare_aid(cfg, record)
    raw.annotations.append(raw.first_time + 30., 1., "BAD_movement")
    with pytest.raises(RuntimeError, match="stale"):
        blink.load_cached_aid(cfg, record)
    raw.annotations.delete(len(raw.annotations) - 1)
    path = blink.cache_paths(cfg, record)["sources"]
    with path.open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(RuntimeError, match="files changed"):
        blink.load_cached_aid(cfg, record)


def test_cache_isolation_ignores_huge_samples_inside_a_known_bad_period(tmp_path, monkeypatch):
    cfg, raw, source, record = preparation_example(tmp_path, monkeypatch)
    blink.prepare_aid(cfg, record)
    first = blink.load_cached_aid(cfg, record)["native_sources"].get_data()
    raw._data[:12, 2000:2200] += 1e6
    blink.prepare_aid(cfg, record, overwrite=True)
    second = blink.load_cached_aid(cfg, record)["native_sources"].get_data()
    np.testing.assert_allclose(first, second, equal_nan=True, rtol=0, atol=0)


# -----------------------------------------------------------------------------
# Detection is scale/polarity invariant, gap-aware and explicitly provisional
# -----------------------------------------------------------------------------

def source_example():
    sfreq = 200.
    times = np.arange(8000) / sfreq
    signal = np.random.default_rng(5).normal(scale=.05, size=8000)
    for peak in (10., 25.):
        signal += 5. * np.exp(-.5 * ((times - peak) / .07) ** 2)
    data = np.vstack([signal, -2. * signal])
    data[:, 3000:3400] = np.nan
    return mne.io.RawArray(data, mne.create_info(["ICA000", "ICA001"], sfreq, "misc"),
                           first_samp=1397000, verbose="ERROR")


def test_detection_merges_duplicates_and_preserves_native_peak_samples():
    raw = source_example()
    intervals, scores, result = blink.detect_blinks(raw, [0, 1], dict(blink.DEFAULTS))
    assert len(intervals) == 2 and len(result["peaks"]) == 4
    assert result["status"] == "proposals_for_review"
    assert all(any(abs(p["peak_raw"] - raw.first_samp - t * 200) <= 2 for p in result["peaks"]) for t in (10, 25))
    assert np.isnan(scores[3000:3400]).all()
    assert np.isnan(scores[:100]).all()
    assert all(not (a < 3400 and b > 3000) for a, b in intervals)
    positive = blink.detect_blinks(raw, [0], dict(blink.DEFAULTS))[0]
    negative = blink.detect_blinks(raw, [1], dict(blink.DEFAULTS))[0]
    assert positive == negative == intervals
    assert raw.ch_names == ["ICA000", "ICA001"]


def test_empty_or_flat_detection_is_not_reviewed_zero():
    raw = source_example()
    assert blink.detect_blinks(raw, [], dict(blink.DEFAULTS))[2]["status"] == "no_component_selected"
    raw._data[:] = 0.
    intervals, scores, result = blink.detect_blinks(raw, [0], dict(blink.DEFAULTS))
    assert intervals == [] and result["diagnostics"][0]["status"] == "flat_or_unscorable"
    assert "review" in result["status"] and "reviewed" not in result


def test_unusually_broad_peaks_are_flagged_instead_of_discarded():
    raw = source_example()
    times = raw.times
    raw._data[0] += 15 * np.exp(-.5 * ((times - 33.) / .5) ** 2)
    _, _, result = blink.detect_blinks(raw, [0], dict(blink.DEFAULTS))
    assert any(p["long_event"] and abs(p["peak_raw"] - raw.first_samp - 6600) < 20 for p in result["peaks"])


# -----------------------------------------------------------------------------
# Regeneration changes proposals, never the operator's annotation decisions
# -----------------------------------------------------------------------------

def test_regeneration_preserves_deleted_resized_and_preexisting_annotations():
    raw = source_example()
    raw.set_annotations(mne.Annotations([30.], [.2], ["BAD_eye_movement"]))
    state = {}
    blink.replace_proposals(raw, [[1000, 1100], [2000, 2100], [4000, 4100]], state)
    # Delete the first, resize the second; leave the third untouched.
    raw.annotations.delete(0)
    raw.annotations.duration[0] = .7
    blink.replace_proposals(raw, [[1000, 1100], [2000, 2100], [4500, 4600]], state)
    keys = [blink.annotation_key(raw, a) for a in raw.annotations]
    assert not any(a == raw.first_samp + 1000 for a, b, label in keys)
    assert (raw.first_samp + 2000, raw.first_samp + 2140, blink.AUTO_LABEL) in keys
    assert (raw.first_samp + 4500, raw.first_samp + 4600, blink.AUTO_LABEL) in keys
    assert not any(a == raw.first_samp + 4000 for a, b, label in keys)
    assert "BAD_eye_movement" in raw.annotations.description
    assert len(state["removed_spans_raw"]) == 2


def test_reviewed_auto_intervals_are_not_removed_by_regeneration():
    raw = source_example()
    raw.annotations.append(raw.first_time + 5., .5, blink.AUTO_LABEL)
    state = {"previous_auto": [blink.annotation_key(raw, raw.annotations[0])]}
    blink.replace_proposals(raw, [[3000, 3100]], state)
    assert len(raw.annotations) == 2
    assert raw.annotations.onset[0] == raw.first_time + 5.


# -----------------------------------------------------------------------------
# Eye duration proposals: steps, native coordinates, gaps and independent edits
# -----------------------------------------------------------------------------

def eye_source_example():
    """Make two scaled/opposite sources with broad steps and a short spike."""
    rng = np.random.default_rng(71)
    signal = rng.normal(scale=.05, size=6000)
    signal[1200:1340] += 4.
    signal[4000:4240] -= 5.
    signal[600:602] += 20.  # Too brief for the duration aid.
    signal[2000:2400] = np.nan
    return mne.io.RawArray(np.vstack([signal, -2. * signal]),
        mne.create_info(["ICA000", "ICA001"], 200., "misc"), first_samp=1397000, verbose="ERROR")


def test_eye_duration_detection_handles_steps_polarity_and_native_gaps():
    raw = eye_source_example()
    before = raw.get_data().copy()
    settings = dict(blink.EYE_DEFAULTS)
    intervals, scores, result = blink.detect_eye_movements(raw, [0, 1], settings)
    assert len(intervals) == 2 and len(result["events"]) == 4
    assert intervals == blink.detect_eye_movements(raw, [0], settings)[0]
    assert intervals == blink.detect_eye_movements(raw, [1], settings)[0]
    assert all(not (a < 2400 and b > 2000) for a, b in intervals)
    assert np.isnan(scores[1990:2410]).all() and np.isnan(scores[:10]).all()
    assert all(b - a > 100 for a, b in intervals)  # Broad excursions retain duration.
    assert all(row["start_raw"] >= raw.first_samp for row in result["events"])
    assert result["median_window_samples"] == 21 and result["smoothing_guard_samples"] == 10
    assert [raw.first_samp + 1990, raw.first_samp + 2410] in result["unscored_intervals_raw"]
    assert result["method"] == blink.EYE_METHOD and result["components"] == [0, 1]
    np.testing.assert_array_equal(raw.get_data(), before)


def test_eye_duration_merges_only_supported_bouts_and_clips_at_gaps():
    raw = eye_source_example()
    raw._data[0, 1500:1540] += 5.
    # An 80 ms gap survives the 100 ms median's majority rule, but the
    # separate 100 ms merge setting should still join these supported bouts.
    raw._data[0, 1556:1596] += 5.
    settings = {**blink.EYE_DEFAULTS, "padding_sec": 0.}
    separated = blink.detect_eye_movements(raw, [0], {**settings, "merge_gap_sec": 0.})[0]
    merged = blink.detect_eye_movements(raw, [0], settings)[0]
    assert len(separated) == len(merged) + 1
    raw._data[0, 1960:2000] += 5.
    raw._data[0, 2400:2450] += 5.
    intervals, _, _ = blink.detect_eye_movements(raw, [0], {**settings, "merge_gap_sec": 3., "padding_sec": 1.})
    assert all(not (a < 2410 and b > 1990) for a, b in intervals)


@pytest.mark.parametrize("values", [0., np.nan])
def test_eye_empty_flat_and_unscorable_cases_never_claim_artifact_absence(values):
    raw = eye_source_example()
    assert blink.detect_eye_movements(raw, [], blink.EYE_DEFAULTS)[2]["status"] == "no_component_selected"
    raw._data[:] = values
    intervals, scores, result = blink.detect_eye_movements(raw, [1], blink.EYE_DEFAULTS)
    assert intervals == [] and np.isnan(scores).all()
    assert result["diagnostics"] == [{"component": 1, "status": "flat_or_unscorable"}]
    assert result["status"] == "no_excursions_review_still_required" and "reviewed" not in result


def test_automatic_categories_overlap_but_have_independent_vetoes_and_edits():
    raw = eye_source_example()
    blink_state, eye_state = {}, {}
    blink.replace_proposals(raw, [[1000, 1100]], blink_state)
    blink.replace_proposals(raw, [[1000, 1100]], eye_state, label=blink.EYE_AUTO_LABEL)
    assert set(raw.annotations.description) == {blink.AUTO_LABEL, blink.EYE_AUTO_LABEL}
    raw.annotations.delete([i for i, label in enumerate(raw.annotations.description) if label == blink.AUTO_LABEL])
    blink.replace_proposals(raw, [[1000, 1100]], blink_state)
    assert list(raw.annotations.description) == [blink.EYE_AUTO_LABEL]
    assert blink_state["removed_spans_raw"] and not eye_state["removed_spans_raw"]
    raw.annotations.duration[0] = .7
    blink.replace_proposals(raw, [[1000, 1100], [1400, 1500]], eye_state, label=blink.EYE_AUTO_LABEL)
    assert len(raw.annotations) == 2 and raw.annotations.duration[0] == .7
    # Disabling one group removes only untouched session proposals, not edits.
    blink.replace_proposals(raw, [], eye_state, label=blink.EYE_AUTO_LABEL)
    assert len(raw.annotations) == 1 and raw.annotations.duration[0] == .7


def test_batch_reports_failures_without_calling_a_final_stage(tmp_path, monkeypatch):
    cfg, raw, source, record = preparation_example(tmp_path, monkeypatch)
    other = AcquisitionRunRecord(str(source), "sub99_run03_test.fif", "sub99_run03_test", 3, "source")
    monkeypatch.setattr(blink, "load_and_resolve_config", lambda *a, **k: (cfg,))
    monkeypatch.setattr(blink, "require_review", lambda *a, **k: None)
    monkeypatch.setattr(blink, "discover_raw_run_records", lambda cfg: [record, other])
    seen = []

    def prepare(cfg, selected, **kwargs):
        seen.append(selected.source_filename)
        if selected is record:
            raise RuntimeError("Unusable fixture")
        return {"execution_status": "prepared", "fingerprint": "test", "n_components_fitted": 10, "converged": True}

    monkeypatch.setattr(blink, "prepare_aid", prepare)
    assert blink.main(["--config", "x", "--subject-config", "x", "--tuning-config", "x", "--resume"]) == 2
    report = json.loads((Path(cfg["resolved_paths"]["derivatives_reports_dir"]) / "03_B_blink_aid_preparation.json").read_text())
    assert len(seen) == 2 and report["status"] == "failed"
    assert report["scientific_acceptance"] is False
    assert [r["status"] for r in report["runs"]] == ["failed", "prepared"]


@pytest.mark.parametrize("manual_muscle", [None, "edited", "empty", "unchanged"])
def test_existing_muscle_context_is_source_bound_and_all_settings_are_checked(tmp_path, manual_muscle):
    """Exercise the real header-only input loader without any cohort recording."""
    from continuous_cleaning import digest
    preprocessing = tmp_path / "preprocessing"
    sss = preprocessing / "sss"; sss.mkdir(parents=True)
    reports = tmp_path / "reports"; reports.mkdir()
    cfg = {"artifact_annotation": {"muscle_threshold": 10., "muscle_threshold_by_run": {},
        "muscle_filter_freq_hz": [110., 140.], "muscle_channel_type": "mag", "muscle_min_length_good_sec": .2},
        "resolved_paths": {"derivatives_preprocessing_dir": str(preprocessing), "derivatives_reports_dir": str(reports)}}
    record = AcquisitionRunRecord("source", "sub99_run01_test.fif", "sub99_run01_test", 1, "source")
    raw = mne.io.RawArray(np.zeros((2, 1000)), mne.create_info(["MEG0121", "MEG0811"], 200., "mag"),
                           first_samp=2000, verbose="ERROR")
    source = sss / "sub99_run01_test_raw_sss.fif"; raw.save(source, overwrite=True, verbose="ERROR")
    signature = digest({"input": file_stamp(source), "threshold": 10.,
        "muscle_settings": {k: v for k, v in cfg["artifact_annotation"].items() if k.startswith("muscle_") and k != "muscle_threshold_by_run"},
        "eye_review": None})
    annotated = raw.copy()
    annotated.info["description"] = f"stage4_input_sha256={signature}"
    annotated.set_annotations(mne.Annotations([1., 3.], [.1, .1], ["BAD_muscle", "BAD_blink"]))
    output = sss / "sub99_run01_test_ann.fif"; annotated.save(output, overwrite=True, verbose="ERROR")
    report = {"runs": [{"source_filename": record.source_filename, "input_stamp": file_stamp(source), "muscle_threshold": 10.}],
        "evidence_files": [{"relative_path": str(output.relative_to(tmp_path)), "stamp": file_stamp(output)}]}
    (reports / "03_B_ArtefactAnnotation_review_report.json").write_text(json.dumps(report))
    loaded, _ = blink.load_aid_input(cfg, record)
    assert loaded.annotations.description.tolist() == ["BAD_muscle"]
    assert loaded.annotations.onset.tolist() == [raw.first_time + 1.]
    if manual_muscle is not None:
        from continuous_cleaning import muscle_parameters
        before = blink.annotation_rows(loaded)
        manual = mne.Annotations([], [], [], orig_time=raw.annotations.orig_time)
        if manual_muscle == "edited":
            manual.append(raw.first_time + 2., .4, "BAD_muscle")
        elif manual_muscle == "unchanged":
            current = mne.io.read_raw_fif(output, preload=False, verbose="ERROR").annotations
            manual = current[[i for i, label in enumerate(current.description) if label == "BAD_muscle"]]
        folder = preprocessing / "manual_eye"; folder.mkdir()
        path = folder / f"{record.source_stem}-annot.fif"; manual.save(path)
        cfg["artifact_annotation"]["manual_eye_reviews"] = {record.source_filename: {
            "status": "reviewed", "sss_input": file_stamp(source), "annotation_sha256": file_hash(path),
            "muscle_review": {"status": "reviewed", "settings_sha256": digest(muscle_parameters(cfg["artifact_annotation"], 10.))}}}
        changed, _ = blink.load_aid_input(cfg, record)
        if manual_muscle == "unchanged":
            assert blink.annotation_rows(changed) == before  # Merely saving must keep compatible aids.
        elif manual_muscle == "empty":
            assert blink.annotation_rows(changed) != before
            assert not len(changed.annotations)
        else:
            assert blink.annotation_rows(changed) != before
            np.testing.assert_allclose(changed.annotations.onset, [raw.first_time + 2.], atol=2e-6)
            np.testing.assert_allclose(changed.annotations.duration, [.4], atol=2e-6)
    cfg["artifact_annotation"]["muscle_filter_freq_hz"] = [100., 130.]
    with pytest.raises(RuntimeError, match="settings changed"):
        blink.load_aid_input(cfg, record)
