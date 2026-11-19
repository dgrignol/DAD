"""Synthetic regression checks for bulk ocular selection and annotation retention.

Usage
-----
    python -B -m pytest -p no:cacheprovider meg/tests/test_annotation_triage.py -q

Temporary native-grid recordings exercise real MNE annotation rejection parity,
manual precedence (including reviewed-empty muscle), cache-bound choice reuse,
resumable explicit selection, exact-mask targeted review and complete zero-count
reports. GUI interactions are mocked; static report rendering is real. No DAD
source, derivative, human review receipt or scheduler is used by these tests.
"""

from datetime import datetime, timezone
import argparse
from copy import deepcopy
import csv
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import annotation_triage as triage
import review_ocular_components as bulk
import manual_eye_annotations as manual
from continuous_cleaning import digest, muscle_parameters
from preprocessing_stage_utils import AcquisitionRunRecord


# -----------------------------------------------------------------------------
# Native coordinate fixtures and MNE parity, including empty survivors
# -----------------------------------------------------------------------------

def example(tmp_path, *, date=True):
    """Build a nonzero-origin recording, small selector and exact-source cache."""
    raw = mne.io.RawArray(np.zeros((2, 1000)), mne.create_info(["ICA000", "ICA001"], 100., "misc"),
                          first_samp=1397000, verbose="ERROR")
    if date:
        raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    cfg = {"layout": {"subject_bids_id": "sub-06"}, "subject": {"id": "06"},
        "resolved_paths": {"derivatives_root": str(tmp_path / "derivatives"), "data_root": str(tmp_path / "raw"),
            "derivatives_preprocessing_dir": str(tmp_path / "derivatives/preprocessing/sub-06/ses-01"),
            "derivatives_reports_dir": str(tmp_path / "derivatives/reports/sub-06/ses-01")},
        "artifact_annotation": {"muscle_threshold": 10., "muscle_filter_freq_hz": [110., 140.]},
        "subject_profile": {"identity_range": [1, 40]},
        "epoching": {"tmin_sec": -.1, "tmax_sec": .2, "reject_by_annotation": True,
                     "event_id_map": {"visible": 71, "occluded": 81, "deviant": 91}}}
    record = AcquisitionRunRecord(str(tmp_path / "raw/test.fif"), "sub06_run04_test.fif", "sub06_run04_test", 4, "test.fif")
    cache = {"model": SimpleNamespace(n_components_=2), "native_sources": raw,
        "metadata": {"fingerprint": "cache-identity", "model_sha256": "model-identity", "sources_stamp": {"size": 1}}}
    return cfg, raw, record, cache


def rows_for(raw, offsets):
    """Preserve digital source keys while supplying accepted corrected onsets."""
    return [{"trial_key": f"trial-{i}", "source_fif": "sub06_run04_test.fif", "sfreq_hz": raw.info["sfreq"],
        "recording_first_sample_raw": raw.first_samp,
        "recording_end_sample_raw_exclusive": raw.first_samp + raw.n_times,
        "condition_code": [71, 81, 91][i % 3], "trial_sequence_identity": i + 1,
        "included_by_trial_policy": False, "final_analysis_eligible": True,
        "epoch_event_sample_raw": raw.first_samp + offset} for i, offset in enumerate(offsets)]


def mask(raw, *, manual_review=True):
    """Represent current mask rows without modifying the source object."""
    return {"schema": triage.SCHEMA, "purpose": triage.PURPOSE,
        "method": "preserved_manual" if manual_review else "provisional_automatic",
        "coverage_status": "manual_review_preserved" if manual_review else "selection_missing",
        "annotations": triage.mask_rows(raw), "detectors": {}, "input_binding": {"source": "unchanged"}}


@pytest.mark.parametrize("date", [False, True])
def test_native_annotation_outcomes_match_real_mne_including_boundaries_and_zero_spans(tmp_path, date):
    cfg, raw, _, _ = example(tmp_path, date=date)
    raw.set_annotations(mne.Annotations([1., 1.1, 5., 0.], [.3, .1, 0., 9.],
        ["BAD_muscle", "BAD_blink_auto", "bad_zero_duration", "BAD_filter_edge"]))
    rows = rows_for(raw, [0, 100, 130, 150, 500, 600, 999])
    counted = triage.count_trials(rows, raw, mask(raw), cfg["epoching"])
    comparison = raw.copy()
    comparison.annotations.delete([i for i, a in enumerate(comparison.annotations) if a["description"] == "BAD_filter_edge"])
    events = np.array([[r["epoch_event_sample_raw"], 0, r["condition_code"]] for r in rows])
    epochs = mne.Epochs(comparison, events, tmin=-.1, tmax=.2, baseline=None, proj=False,
                        reject_by_annotation=True, preload=True, verbose="ERROR")
    assert [i for i, r in enumerate(counted) if r["annotation_survives"]] == list(epochs.selection)
    assert counted[0]["incomplete_window"] and counted[-1]["incomplete_window"]
    assert counted[1]["bad_categories"] == "BAD_blink_auto;BAD_muscle"
    assert all("filter_edge" not in r["bad_categories"] for r in counted)
    np.testing.assert_array_equal(raw.get_data(), np.zeros((2, 1000)))


def test_accepted_overlay_and_unscored_none_are_not_artifact_free(tmp_path):
    cfg, raw, _, _ = example(tmp_path)
    rows = rows_for(raw, [500, 600])
    rows[1]["included_by_trial_policy"] = True
    rows[1]["final_analysis_eligible"] = False
    counted = triage.count_trials(rows, raw, mask(raw, manual_review=False), cfg["epoching"])
    assert counted[0]["annotation_survives"]
    assert counted[0]["unscored_groups"] == "blink;eye_movement"
    assert counted[1]["triage_outcome"] == "excluded_by_accepted_trial_policy"
    wrong = deepcopy(rows); wrong[0]["recording_first_sample_raw"] += 1
    with pytest.raises(RuntimeError, match="coordinates"):
        triage.count_trials(wrong, raw, mask(raw), cfg["epoching"])


# -----------------------------------------------------------------------------
# Protected manual precedence, model bindings and exact-mask reopening
# -----------------------------------------------------------------------------

def test_reviewed_empty_muscle_replaces_auto_while_legacy_review_keeps_it(tmp_path, monkeypatch):
    cfg, raw, record, _ = example(tmp_path)
    source = tmp_path / "source_sss.fif"; source.write_bytes(b"source")
    monkeypatch.setattr(triage, "derivative_path_for_run", lambda *a: source)
    monkeypatch.setattr(mne.io, "read_raw_fif", lambda *a, **k: raw.copy())
    muscle = mne.Annotations([raw.first_time + 2.], [.5], ["BAD_muscle"], orig_time=raw.annotations.orig_time)
    # Match the source's acquisition annotation origin even for absent meas_date.
    monkeypatch.setattr(manual, "load_muscle_context", lambda *a: (muscle, {"threshold": 10.}))
    monkeypatch.setattr(triage, "has_manual_work", lambda *a: True)
    reviewed = mne.Annotations([raw.first_time + 4.], [.1], ["BAD_blink"], orig_time=raw.annotations.orig_time)
    receipt = {"status": "reviewed", "muscle_review": {"status": "reviewed",
        "settings_sha256": digest(muscle_parameters(cfg["artifact_annotation"], 10.)), "annotation_count": 0}}
    monkeypatch.setattr(manual, "load_manual", lambda *a: (reviewed, receipt))
    result, _, _ = triage.load_context(cfg, record)
    assert list(result.annotations.description) == ["BAD_blink"]
    receipt.pop("muscle_review")
    result, _, _ = triage.load_context(cfg, record)
    assert set(result.annotations.description) == {"BAD_muscle", "BAD_blink"}
    assert len(raw.annotations) == 0


def test_manual_masks_ignore_all_bulk_choices_and_do_not_load_cache(tmp_path, monkeypatch):
    cfg, raw, record, _ = example(tmp_path)
    raw.annotations.append(raw.first_time + 3, .5, "BAD_blink")
    monkeypatch.setattr(triage, "load_context", lambda *a: (raw.copy(), None, {"manual_receipt_sha256": "reviewed"}))
    monkeypatch.setattr(triage, "load_cached_aid", lambda *a: pytest.fail("Manual work must not load automatic cache"))
    _, result = triage.provisional_mask(cfg, record)
    assert result["method"] == "preserved_manual"
    assert result["detectors"] == {} and result["selection_sha256"] is None


def selected(record, cache):
    """Construct an explicit human selection; [] means deliberate NONE."""
    return {"schema": triage.SCHEMA, "purpose": "ocular_component_selection_only", "status": "selected",
            "binding": triage.component_binding(record, cache), "blink_components": [0], "eye_movement_components": []}


def test_selection_staleness_and_atomic_conflict_preserve_previous_choices(tmp_path):
    cfg, _, record, cache = example(tmp_path)
    choice = selected(record, cache)
    path = bulk.save_selection(cfg, record, choice)
    assert triage.validate_selection(choice, record, cache) == choice
    changed = deepcopy(cache); changed["metadata"]["model_sha256"] = "new-model"
    with pytest.raises(RuntimeError, match="stale"):
        triage.validate_selection(choice, record, changed)
    with pytest.raises(RuntimeError, match="another session"):
        bulk.save_selection(cfg, record, {**choice, "blink_components": [1]})
    assert json.loads(path.read_text()) == choice
    assert (path.parent / "component_selections.csv").is_file()
    bad = deepcopy(choice); bad["blink_components"] = [True]
    with pytest.raises(ValueError, match="Invalid"):
        triage.validate_selection(bad, record, cache)


def test_provisional_detection_none_gaps_and_hash_bound_targeted_reopening(tmp_path, monkeypatch):
    cfg, raw, record, cache = example(tmp_path)
    binding = {"manual_receipt_sha256": None, "source_filename": record.source_filename}
    monkeypatch.setattr(triage, "load_context", lambda *a: (raw.copy(), None, binding))
    monkeypatch.setattr(triage, "load_cached_aid", lambda *a: cache)
    monkeypatch.setattr(triage, "has_manual_work", lambda *a: False)
    choice = selected(record, cache)
    bulk.save_selection(cfg, record, choice)
    scores = np.ones(raw.n_times); scores[:50] = np.nan
    monkeypatch.setattr(triage, "detect_blinks", lambda *a: ([[100, 130]], scores,
        {"settings": triage.blink_settings(cfg), "status": "proposals_for_review"}))
    _, proposal = triage.provisional_mask(cfg, record)
    assert proposal["annotations"][0]["onset"] == (raw.first_samp + 100) / 100.
    assert proposal["detectors"]["blink"]["unscored_intervals_raw"] == [[raw.first_samp, raw.first_samp + 50]]
    assert proposal["detectors"]["eye_movement"]["unscored_intervals_raw"] == [[raw.first_samp, raw.first_samp + raw.n_times]]
    path = triage.save_mask(cfg, record, proposal)
    assert triage.load_review_mask(cfg, record, path) == proposal
    monkeypatch.setattr(triage, "detect_blinks", lambda *a: pytest.fail("Compatible detector results should be reused"))
    assert triage.provisional_mask(cfg, record)[1] == proposal
    seeded = raw.copy(); triage.seed_review_mask(seeded, proposal)
    assert triage.mask_rows(seeded) == proposal["annotations"]
    assert len(raw.annotations) == 0
    cfg["artifact_annotation"]["blink_auto"] = {"threshold_std": 7}
    with pytest.raises(RuntimeError, match="settings changed"):
        triage.load_review_mask(cfg, record, path)
    cfg["artifact_annotation"].pop("blink_auto")
    monkeypatch.setattr(triage, "has_manual_work", lambda *a: True)
    with pytest.raises(RuntimeError, match="Manual work already exists"):
        triage.load_review_mask(cfg, record, path)


# -----------------------------------------------------------------------------
# Resumable UI and end-to-end report snapshots without production writes
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("text,expected", [("0,ICA001,0", [0, 1]), ("NONE", []), ("defer", "DEFER"), ("QUIT", "QUIT")])
def test_component_input(text, expected):
    assert bulk.parse_components(text, 2) == expected


@pytest.mark.parametrize("text", ["", "-1", "2", "NONE,0", "ICAxyz"])
def test_component_input_rejects_invalid_answers(text):
    with pytest.raises(ValueError):
        bulk.parse_components(text, 2)


def bulk_fixture(tmp_path, monkeypatch):
    cfg, raw, record, cache = example(tmp_path)
    monkeypatch.setattr(bulk, "subject_config", lambda *a: cfg)
    monkeypatch.setattr(bulk, "require_review", lambda *a: None)
    monkeypatch.setattr(bulk, "discover_raw_run_records", lambda *a: [record])
    monkeypatch.setattr(bulk, "load_cached_aid", lambda *a: cache)
    monkeypatch.setattr(bulk, "has_manual_work", lambda *a: False)
    args = argparse.Namespace(subjects=[6], source=None, revisit=False, check_only=False,
                              reviewer="User", browser_backend="matplotlib")
    return cfg, raw, record, cache, args


def test_bulk_check_only_is_read_only_and_save_resume_skips_gui(tmp_path, monkeypatch):
    cfg, _, record, _, args = bulk_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(bulk, "show_components", lambda *a: None)
    args.check_only = True
    assert bulk.run_review(args) == 0
    assert not Path(cfg["resolved_paths"]["derivatives_root"]).exists()
    args.check_only = False
    answers = iter(["0", "NONE", "SAVE", "test note"])
    monkeypatch.setattr("builtins.input", lambda *a: next(answers))
    assert bulk.run_review(args) == 0
    saved = json.loads(triage.selection_path(cfg, record).read_text())
    assert saved["blink_components"] == [0] and saved["eye_movement_components"] == []
    monkeypatch.setattr(bulk, "show_components", lambda *a: pytest.fail("Resume must skip saved choices"))
    assert bulk.run_review(args) == 0
    monkeypatch.setattr(bulk, "has_manual_work", lambda *a: True)
    monkeypatch.setattr(bulk, "load_context", lambda *a: None)
    args.revisit = True
    assert bulk.run_review(args) == 0  # Manual still wins over --revisit.


def test_report_zero_survivors_is_valid_and_all_120_identity_cells_are_present(tmp_path, monkeypatch):
    cfg, raw, record, _ = example(tmp_path)
    raw.annotations.append(raw.first_time, 10., "BAD_muscle")
    proposal = mask(raw); proposal["payload_sha256"] = digest(proposal)
    monkeypatch.setattr(triage, "subject_config", lambda *a: cfg)
    monkeypatch.setattr(triage, "require_review", lambda *a: None)
    monkeypatch.setattr(triage, "discover_raw_run_records", lambda *a: [record])
    monkeypatch.setattr(triage, "load_prepared_decisions", lambda *a: rows_for(raw, [300, 400, 500]))
    monkeypatch.setattr(triage, "load_reviewed_timing", lambda c, r, p: (r, {"version": "accepted-test"}))
    monkeypatch.setattr(triage, "provisional_mask", lambda *a: (raw, proposal))
    monkeypatch.setattr(triage, "load_context", lambda *a: (raw, None, proposal["input_binding"]))
    args = argparse.Namespace(subjects=[6], label="zero-survivors")
    assert triage.run_report(args) == 0
    output = triage.triage_root(cfg) / "reports/zero-survivors"
    report = json.loads((output / "report.json").read_text())
    assert report["subjects"][0]["minimum_condition_survivors"] == 0
    assert report["subjects"][0]["missing_identity_cells"] == 120
    with (output / "identity_retention.csv").open() as handle:
        assert len(list(csv.DictReader(handle))) == 120
    assert (output / "condition_retention.png").stat().st_size > 1000
    with pytest.raises(FileExistsError):
        triage.run_report(args)
    assert not Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]).exists()


def test_report_input_failure_is_unranked_instead_of_zero_count(tmp_path, monkeypatch):
    cfg, _, _, _ = example(tmp_path)
    monkeypatch.setattr(triage, "subject_config", lambda *a: cfg)
    monkeypatch.setattr(triage, "require_review", lambda *a: (_ for _ in ()).throw(RuntimeError("stale source")))
    assert triage.run_report(argparse.Namespace(subjects=[6], label="error-case")) == 1
    report = json.loads((triage.triage_root(cfg) / "reports/error-case/report.json").read_text())
    assert report["subjects"][0] == {"subject_bids_id": "sub-06", "status": "incomplete_inputs"}
    assert report["errors"][0]["error"] == "stale source"
