"""Stage 4 scientific-boundary regression tests, using temporary synthetic data.

Usage: python -B -m pytest -p no:cacheprovider meg/tests/test_stage4_continuous.py -q

Large artifacts are deliberately varied inside every excluded category to test
that filtering/resampling cannot leak them into the actual ICA fitting samples.
Tests also exercise native acquisition coordinates, reviewed empty ocular files,
final-filter edge handling and changed-evidence rejection. No cohort derivative
is read or written by this module.
"""

import copy
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
from continuous_cleaning import (annotation_rows, bind_report, check_preservation,
    coordinates, digest, file_hash, file_stamp, filter_continuous, final_filter_summary,
    good_stretches, make_filter_edges_informational, prepare_fit_copy, require_review, stage_configuration)
from manual_eye_annotations import eye_annotations, load_manual, manual_path
from preprocessing_stage_utils import AcquisitionRunRecord


FIT = {"fit_l_freq_hz": 1., "fit_h_freq_hz": 40., "fit_resample_sfreq_hz": 200.,
       "fit_resample_edge_sec": .1, "fit_min_good_sec": 1.}


def raw_example(*, date=True, duration=60.):
    """Make MEG plus auxiliary data with a nonzero acquisition origin."""
    sfreq = 400.
    rng = np.random.default_rng(7)
    data = rng.normal(size=(5, int(duration * sfreq)))
    data[:3] *= 1e-12
    data[4] = np.arange(data.shape[1]) % 200
    info = mne.create_info(["MEG0111", "MEG0121", "MEG0311", "MISC001", "STI101"], sfreq,
                           ["mag", "mag", "mag", "misc", "stim"])
    raw = mne.io.RawArray(data, info, first_samp=4000, verbose="ERROR")
    if date:
        raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    return raw


@pytest.mark.parametrize("date", [True, False])
@pytest.mark.parametrize("historical_margin", [None, 165.])
@pytest.mark.parametrize("category", ["BAD_blink", "BAD_squid_jumps"])
def test_no_highpass_still_removes_lines_and_preserves_native_recording(date, historical_margin, category, monkeypatch):
    """No FIR call/margin/exposure; real spectrum-fit cannot bridge BAD/EDGE.

    Very large excluded artifacts are varied while every retained sample must
    stay identical. This exercises the actual MNE annotation-aware operation,
    including short good stretches and nonzero acquisition coordinates.
    """
    raw = raw_example(date=date)
    raw.set_annotations(mne.Annotations([20., 30., 30.25], [2., 0., .5],
                                       [category, "EDGE boundary", "BAD_muscle"]))
    with raw.info._unlock():
        raw.info["highpass"] = .1  # Irreversible P17/P24-style acquisition limit.
    times = np.arange(raw.n_times) / raw.info["sfreq"]
    for frequency in (50., 100.):
        raw._data[0] += 2e-11 * np.sin(2 * np.pi * frequency * times)
    before = raw.copy()
    damaged = raw.copy()
    damaged._data[:3, 8000:8800] += 1e9
    damaged._data[:3, 12100:12300] -= 1e12
    settings = {"highpass_hz": None, "filter_edge_margin_sec": historical_margin,
                "notch_freqs_hz": [50., 100.], "spectrum_fit_window": "10s"}
    def unexpected_filter(*args, **kwargs):
        raise AssertionError("A null high-pass must never call Raw.filter")
    monkeypatch.setattr(type(raw), "filter", unexpected_filter)
    summary = filter_continuous(raw, settings)
    filter_continuous(damaged, settings)
    for start, stop in good_stretches(before):
        np.testing.assert_array_equal(raw.get_data(start=start, stop=stop),
                                      damaged.get_data(start=start, stop=stop))
    check_preservation(before, raw)
    np.testing.assert_array_equal(raw.get_data()[:3, 8000:8800], before.get_data()[:3, 8000:8800])
    assert raw.info["highpass"] == .1
    assert summary["filter"] is None and summary["candidate_margin_sec"] is None
    assert summary["reviewed_margin_sec"] is None and summary["edge_exposed_good_seconds"] == 0
    assert summary["spectrum_fit"]["window"] == "10s"
    assert min(summary["good_stretch_durations_sec"]) == .25
    assert summary["epoch_exclusions_from_filter_edges"] == 0
    for frequency in (50., 100.):
        wave = np.exp(-2j * np.pi * frequency * times[800:7200])
        original = abs(np.mean(before.get_data(picks=[0])[0, 800:7200] * wave))
        cleaned = abs(np.mean(raw.get_data(picks=[0])[0, 800:7200] * wave))
        assert cleaned < .1 * original


def legacy_sss_fixture(tmp_path):
    """Create old-format accepted metadata with independent legacy hashing."""
    import yaml
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] / "params/pipeline_default.yaml").read_text())
    cfg["preprocessing"]["highpass_hz"] = .01
    resource = tmp_path / "resource"; resource.write_bytes(b"calibration")
    for key in ("calibration_file", "crosstalk_file"):
        cfg["maxwell_filter"][key] = str(resource)
    old_maxwell = dict(cfg["maxwell_filter"])
    for key in ("calibration_file", "crosstalk_file"):
        old_maxwell[key] = file_hash(resource)
    old_filters = dict(cfg["preprocessing"])
    old_filters.pop("filter_edge_margin_sec", None)
    old_filters.pop("final_filter_skip_sources", None)
    legacy = digest({"maxwell": old_maxwell, "continuous_filter": old_filters})
    reports = tmp_path / "reports/sub-04/ses-01"; reports.mkdir(parents=True)
    raw_root = tmp_path / "raw"; raw_root.mkdir()
    source = raw_root / "source.fif"; source.write_bytes(b"raw fingerprint only")
    output = tmp_path / "sss.fif"; output.write_bytes(b"SSS fingerprint only")
    payload = {"configuration_sha256": legacy,
               "source_files": [{"relative_path": source.name, "stamp": file_stamp(source)}]}
    bind_report(payload, [output], {"old_policy": .01}, reports, tmp_path)
    report = reports / "03_A_MaxFilt_review_report.json"; report.write_text(json.dumps(payload))
    receipt = {"status": "accepted", "configuration_sha256": legacy,
               "evidence_sha256": payload["evidence_sha256"], "reason": "Historical acceptance; retain verbatim"}
    cfg["stage4_review"] = {"03_A_MaxFilt": receipt}
    cfg["resolved_paths"] = {"data_root": str(raw_root)}
    cfg["preprocessing"]["highpass_hz"] = None
    return cfg, reports, report, source, output, resource


def test_accepted_legacy_sss_survives_final_filter_change_without_rewriting_history(tmp_path):
    cfg, reports, report, *_ = legacy_sss_fixture(tmp_path)
    original_report = report.read_bytes()
    original_receipt = copy.deepcopy(cfg["stage4_review"])
    original_a = stage_configuration(cfg, "03_A_MaxFilt")
    original_b = stage_configuration(cfg, "03_B_ArtefactAnnotation")
    original_c = stage_configuration(cfg, "03_C_ICA")
    require_review(cfg, "03_A_MaxFilt", reports)
    cfg["preprocessing"].update(notch_freqs_hz=[50.], spectrum_fit_window="2s")
    require_review(cfg, "03_A_MaxFilt", reports)
    assert stage_configuration(cfg, "03_A_MaxFilt") == original_a
    assert stage_configuration(cfg, "03_B_ArtefactAnnotation") != original_b
    assert stage_configuration(cfg, "03_C_ICA") == original_c
    assert report.read_bytes() == original_report
    assert cfg["stage4_review"] == original_receipt


@pytest.mark.parametrize("action", ["reviewed_zero", "interval", "union_channel", "run_channel"])
def test_final_policy_legacy_compatibility_does_not_accept_changed_squid_treatment(tmp_path, action):
    """Refresh ordering depends on operative SSS changes, not candidate/UI state.

    Empty SQUID review plus a final-only change preserves old 4A/4B input
    bindings. Any new interval or channel exclusion must fail that exception
    before downstream work; acceptance history stays byte-identical.
    """
    cfg, reports, report, *_ = legacy_sss_fixture(tmp_path)
    original_report = report.read_bytes()
    original_acceptance = copy.deepcopy(cfg["stage4_review"])
    applied = {"intervals": [], "bad_channels": [], "channel_scope": "union"}
    cfg["squid_jumps"] = {"enabled": True, "threshold": 8., "reviews": {
        "source.fif": {"status": "reviewed", "application": applied}}}
    cfg["preprocessing"]["highpass_hz"] = None
    if action == "interval":
        applied["intervals"] = [[14000, 14500]]
    elif action in ("union_channel", "run_channel"):
        applied["bad_channels"] = ["MEG0111"]
        applied["channel_scope"] = "union" if action == "union_channel" else "run"
    if action == "reviewed_zero":
        require_review(cfg, "03_A_MaxFilt", reports)
    else:
        with pytest.raises(RuntimeError, match="configuration changed"):
            require_review(cfg, "03_A_MaxFilt", reports)
    assert report.read_bytes() == original_report
    assert cfg["stage4_review"] == original_acceptance


@pytest.mark.parametrize("change", ["calibration", "destination", "bads", "sss_filter", "other_preprocessing",
                                          "raw", "output", "report_binding", "evidence_binding"])
def test_legacy_sss_migration_rejects_unrelated_changes(tmp_path, change):
    cfg, reports, report, source, output, resource = legacy_sss_fixture(tmp_path)
    if change == "calibration": resource.write_bytes(b"changed")
    elif change == "destination": cfg["maxwell_filter"]["destination"]["policy"] = "none"
    elif change == "bads": cfg["maxwell_filter"]["manual_bad_channels"]["add"] = ["MEG0111"]
    elif change == "sss_filter": cfg["maxwell_filter"]["filter"]["origin"] = [0., 0., .04]
    elif change == "other_preprocessing": cfg["preprocessing"]["preserve_misc_channels"] = False
    elif change == "raw": source.write_bytes(b"changed")
    elif change == "output": output.write_bytes(b"changed")
    elif change == "report_binding":
        payload = json.loads(report.read_text()); payload["configuration_sha256"] = "changed"
        report.write_text(json.dumps(payload))
    else: cfg["stage4_review"]["03_A_MaxFilt"]["evidence_sha256"] = "changed"
    with pytest.raises(RuntimeError):
        require_review(cfg, "03_A_MaxFilt", reports)


def test_filter_preview_refresh_archives_original_report_bytes(tmp_path):
    """Report refresh preserves old FIR exposure and carries history forward."""
    import importlib.util
    path = Path(__file__).resolve().parents[1] / "code/preprocessing_basic/03_B_ArtefactAnnotation.py"
    spec = importlib.util.spec_from_file_location("filter_history_test", path)
    stage = importlib.util.module_from_spec(spec); spec.loader.exec_module(stage)
    report = tmp_path / "review.json"
    old = {"runs": [{"final_filter_preview": {"settings": {"highpass_hz": .01}, "candidate_margin_sec": 165.}}]}
    original = json.dumps(old).encode(); report.write_bytes(original)
    report.with_suffix(".html").write_text("Historical embedded exposure report")
    current = {"runs": [{"final_filter_preview": {"settings": {"highpass_hz": None}}}]}
    stage.preserve_filter_preview_history(current, report, tmp_path / "raw")
    entry = current["historical_final_filter_reports"][0]
    assert (tmp_path / entry["json"]).read_bytes() == original
    assert (tmp_path / entry["html"]).read_text() == "Historical embedded exposure report"
    report.write_text(json.dumps(current))
    again = copy.deepcopy(current)
    stage.preserve_filter_preview_history(again, report, tmp_path / "raw")
    assert again == current


@pytest.mark.parametrize("date", [True, False])
@pytest.mark.parametrize("category", ["BAD_blink", "BAD_eye_movement", "BAD_eye_movement_auto", "BAD_muscle", "BAD_sensor_fault", "bad_jump", "BAD_squid_jumps"])
def test_bad_artifacts_cannot_enter_fit_data(date, category):
    raw = raw_example(date=date)
    raw.set_annotations(mne.Annotations([15., 35.], [2., 3.], [category, category]))
    damaged = raw.copy()
    damaged._data[:3, 6000:6800] += 1e9
    damaged._data[:3, 14000:15200] -= 1e12
    before = coordinates(raw)
    first, summary = prepare_fit_copy(raw, FIT)
    second, _ = prepare_fit_copy(damaged, FIT)
    np.testing.assert_array_equal(first.get_data(), second.get_data())
    assert first.n_times == sum(r["fit_samples"] for r in summary["intervals"])
    assert summary["filter"]["half_support_sec"] > 0
    assert coordinates(raw) == before
    assert all(b <= 6000 or a >= 6800 for a, b in good_stretches(raw))


def test_zero_duration_edge_splits_filtering():
    raw = raw_example()
    raw.set_annotations(mne.Annotations([20], [0], ["EDGE boundary"]))
    assert good_stretches(raw) == [(0, 8000), (8000, 24000)]


@pytest.mark.parametrize("date", [True, False])
@pytest.mark.parametrize("empty", [True, False])
def test_manual_eye_roundtrip_and_reviewed_zero(tmp_path, date, empty):
    raw = raw_example(date=date)
    raw.set_annotations(mne.Annotations([1, 2], [0.3, 0.5], ["BAD_muscle", "BAD_blink"]))
    ann = eye_annotations(raw)
    if empty:
        ann = ann[[]]
    source = tmp_path / "source_raw_sss.fif"
    raw.save(source, overwrite=True, verbose="ERROR")
    record = AcquisitionRunRecord("raw", "sub04_run01_date.fif", "sub04_run01_date", 1, "raw")
    cfg = {"resolved_paths": {"derivatives_preprocessing_dir": str(tmp_path)}, "artifact_annotation": {}}
    path = manual_path(cfg, record)
    path.parent.mkdir(parents=True)
    ann.save(path)
    receipt = {"status": "reviewed", "coordinates": coordinates(raw), "sss_input": file_stamp(source),
               "annotation_sha256": file_hash(path), "annotation_count": len(ann)}
    cfg["artifact_annotation"]["manual_eye_reviews"] = {record.source_filename: receipt}
    loaded, _ = load_manual(cfg, record, raw, source)
    assert len(loaded) == (0 if empty else 2)
    np.testing.assert_allclose(loaded.onset, ann.onset, atol=1e-6)
    if not empty:
        assert loaded.orig_time == raw.annotations.orig_time
    cfg["artifact_annotation"]["manual_eye_reviews"] = {}
    with pytest.raises(RuntimeError, match="review missing"):
        load_manual(cfg, record, raw, source)


def test_final_filter_requires_reviewed_exposure_not_epoch_loss():
    raw = raw_example()
    settings = {"highpass_hz": .01, "filter_edge_margin_sec": None}
    report = final_filter_summary(raw, settings)
    assert report["filter"]["half_support_sec"] == 165.
    assert report["seconds_outside_edge_exposure"] == 0
    assert report["epoch_exclusions_from_filter_edges"] == 0
    assert report["status"] == "filter_design_review_required"
    with pytest.raises(RuntimeError, match="Review"):
        filter_continuous(raw, settings)


@pytest.mark.parametrize("date", [True, False])
def test_final_filter_keeps_coordinates_and_auxiliary_samples(date):
    # Keep a bounded optional-FIR regression; production YAML omits this filter.
    raw = raw_example(date=date)
    raw.set_annotations(mne.Annotations([20], [1], ["BAD_sensor_fault"]))
    before = raw.copy()
    filter_continuous(raw, {"highpass_hz": 1., "filter_edge_margin_sec": 2., "spectrum_fit_window": "2s"})
    check_preservation(before, raw, annotations=False)
    assert "filter_edge_exposure" in raw.annotations.description
    assert "BAD_filter_edge" not in raw.annotations.description
    assert all(row in annotation_rows(raw)["rows"] for row in annotation_rows(before)["rows"])


@pytest.mark.parametrize("edge_name", ["filter_edge_exposure", "BAD_filter_edge", "bad_filter_edge"])
def test_edge_exposed_epochs_survive_but_artifacts_and_incomplete_epochs_do_not(edge_name):
    """Exercise MNE's actual rejection behavior with nonzero first_samp."""
    raw = raw_example()
    raw.set_annotations(mne.Annotations([10., 20.], [3., 3.], [edge_name, "BAD_sensor_fault"]))
    make_filter_edges_informational(raw)
    events = np.array([[raw.first_samp + int(t * raw.info["sfreq"]), 0, 1]
                       for t in [11., 21., 59.]])
    epochs = mne.Epochs(raw, events, event_id={"trial": 1}, tmin=-.5, tmax=2.667,
                        baseline=None, preload=True, reject_by_annotation=True, verbose="ERROR")
    assert epochs.selection.tolist() == [0]
    assert "BAD_sensor_fault" in epochs.drop_log[1]
    assert "TOO_SHORT" in epochs.drop_log[2]
    assert "BAD_sensor_fault" in raw.annotations.description


@pytest.mark.parametrize("highpass", [None, .01])
def test_p38_historical_skip_protects_only_a_future_highpass(highpass, monkeypatch):
    """P38 receives real spectrum-fit even with its untouched legacy tuning map."""
    raw = raw_example(duration=20.)
    times = np.arange(raw.n_times) / raw.info["sfreq"]
    raw._data[0] += 2e-11 * np.sin(2 * np.pi * 50 * times)
    before = raw.copy()
    settings = {"highpass_hz": highpass, "final_filter_skip_sources": {
        "sub38_run01_051026.fif": "User decision: do not apply final filtering to this short run."}}
    def unexpected_highpass(*args, **kwargs):
        raise AssertionError("The historical short-run warning must protect a future high-pass")
    monkeypatch.setattr(type(raw), "filter", unexpected_highpass)
    result = filter_continuous(raw, settings, source_filename="sub38_run01_051026.fif")
    check_preservation(before, raw)
    assert np.std(raw.get_data(picks=[0])) < .2 * np.std(before.get_data(picks=[0]))
    assert result["filter"] is None and result["source_highpass_warning"]
    assert result["spectrum_fit"]["status"] == "selected_for_every_source"
    other = final_filter_summary(raw, settings, source_filename="sub38_run02_051026.fif")
    assert (other["filter"] is None) == (highpass is None)
    with pytest.raises(ValueError, match="exact source_filename"):
        filter_continuous(raw, settings)


def test_final_filter_exception_binds_annotation_review_without_invalidating_sss(tmp_path):
    """A final-only decision changes 4B review identity, not completed SSS."""
    resource = tmp_path / "calibration"; resource.write_text("test resource")
    cfg = {"maxwell_filter": {"calibration_file": str(resource), "crosstalk_file": str(resource)},
           "artifact_annotation": {"muscle_threshold": 10},
           "preprocessing": {"highpass_hz": .01, "filter_edge_margin_sec": None}}
    before_a = stage_configuration(cfg, "03_A_MaxFilt")
    before_b = stage_configuration(cfg, "03_B_ArtefactAnnotation", include_edge=False)
    cfg["preprocessing"]["final_filter_skip_sources"] = {"sub38_run01_051026.fif": "User decision"}
    assert stage_configuration(cfg, "03_A_MaxFilt") == before_a
    assert stage_configuration(cfg, "03_B_ArtefactAnnotation", include_edge=False) != before_b


def test_acceptance_rejects_changed_evidence_and_configuration(tmp_path):
    reports = tmp_path / "reports/sub-04/ses-01"
    reports.mkdir(parents=True)
    data = tmp_path / "input.fif"
    data.write_bytes(b"source")
    cfg = {"ica": {"method": "fastica"}, "resolved_paths": {}}
    payload = {}
    bind_report(payload, [data], {"method": "fastica"}, reports, tmp_path)
    (reports / "03_C_ICA_review_report.json").write_text(json.dumps(payload))
    cfg["stage4_review"] = {"03_C_ICA": {"status": "accepted", "evidence_sha256": payload["evidence_sha256"],
                                           "configuration_sha256": stage_configuration(cfg, "03_C_ICA")}}
    require_review(cfg, "03_C_ICA", reports)
    cfg["ica"]["method"] = "infomax"
    with pytest.raises(RuntimeError, match="configuration changed"):
        require_review(cfg, "03_C_ICA", reports)
    cfg["ica"]["method"] = "fastica"
    data.write_bytes(b"changed")
    with pytest.raises(RuntimeError, match="evidence changed"):
        require_review(cfg, "03_C_ICA", reports)


def test_stage4_clis_separate_review_acceptance_and_application(tmp_path, monkeypatch):
    """Exercise actual stage/review CLIs with temporary FIFs and a mocked SSS solver.

    The expensive physical SSS solver is covered by a bounded real-data smoke
    check separately. This test verifies actual gate wiring, reviewed empty eyes,
    None ICA dimension selection, native application and stale-review rejection.
    The production null-high-pass policy exercises real margin-free 4B receipt
    creation and final spectrum-fit application on temporary synthetic FIFs.
    """
    import importlib.util
    import yaml
    import matplotlib.pyplot as plt
    from preprocessing_stage_utils import update_tuning_file
    import review_cleaning
    code = Path(__file__).resolve().parents[1] / "code/preprocessing_basic"
    modules = []
    for name in ["03_A_MaxFilt", "03_B_ArtefactAnnotation", "03_C_ICA"]:
        spec = importlib.util.spec_from_file_location(name, code / (name + ".py"))
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); modules.append(module)
    maxwell, artifacts, ica = modules
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "params/pipeline_default.yaml").read_text())
    raw_dir = tmp_path / "raw/sub04_220526/meg"; raw_dir.mkdir(parents=True)
    derivatives = tmp_path / "derivatives"
    cfg["project"].update(data_root=str(tmp_path / "raw"), derivatives_root=str(derivatives))
    cfg["preprocessing"].update(highpass_hz=None, spectrum_fit_window="2s")
    # This fixture tests the pre-existing muscle/ICA gates; EDF integration has
    # separate event/recording/sample fixtures and never reads cohort sources.
    cfg["artifact_annotation"]["edf_blink_candidates"]["enabled"] = False
    # Existing muscle/ICA compatibility is independent of the new SQUID gate;
    # its enabled prepare/review/apply path has dedicated native fixtures.
    cfg["squid_jumps"]["enabled"] = False
    for field in ["calibration_file", "crosstalk_file"]:
        resource = tmp_path / field; resource.write_text("synthetic resource")
        cfg["maxwell_filter"][field] = str(resource)
    config = tmp_path / "config.yaml"; config.write_text(yaml.safe_dump(cfg))
    source = raw_dir / "sub04_run01_220526.fif"
    raw = raw_example(duration=30.)
    raw.info["dev_head_t"] = mne.transforms.Transform("meg", "head", np.eye(4))
    raw.save(source, overwrite=True, verbose="ERROR")
    monkeypatch.setattr(maxwell, "validate_auto_origin_digitization", lambda *a, **k: {})
    monkeypatch.setattr(maxwell, "_centrality_head_points", lambda info: np.vstack([np.eye(3), np.zeros(3)]))
    def synthetic_scores(raw, **kwargs):
        """Return the full saved-score structure used by channel diagnostics."""
        picks = mne.pick_types(raw.info, meg=True, exclude=[])
        names = np.array([raw.ch_names[i] for i in picks])
        kinds = np.array([raw.get_channel_types()[i] for i in picks])
        return [], [], {"bins": np.array([[0., 5.]]), "ch_names": names, "ch_types": kinds,
            "scores_noisy": np.zeros((len(picks), 1)), "limits_noisy": np.full((len(picks), 1), 7.),
            "scores_flat": np.full((len(picks), 1), 1e-12), "limits_flat": np.full((len(picks), 1), 1e-15)}
    monkeypatch.setattr(mne.preprocessing, "find_bad_channels_maxwell", synthetic_scores)
    monkeypatch.setattr(mne.preprocessing, "maxwell_filter", lambda raw, **kwargs: raw.copy())
    common = ["--config", str(config), "--subject-config", str(root / "params/sub-04.yaml")]
    assert maxwell.main(common + ["--processing-pass", "review"]) == 0
    report_dir = derivatives / "reports/sub-04/ses-01"
    tuning = report_dir / "03_tuning_choices.yaml"
    accepted = common + ["--tuning-config", str(tuning), "--reviewer", "synthetic test", "--reason", "test-only evidence"]
    assert artifacts.main(common + ["--processing-pass", "review"]) == 2  # A still unreviewed.
    assert review_cleaning.main(accepted + ["--stage", "4A"]) == 0
    run = common + ["--tuning-config", str(tuning)]
    assert artifacts.main(run + ["--processing-pass", "review"]) == 0
    # Adding EDF sidecars must reuse the existing muscle FIF under --resume.
    # No EDF input is available in this synthetic fixture: report unavailable
    # guidance explicitly, rather than inventing reviewed-zero blinks.
    enabled = copy.deepcopy(cfg)
    enabled["artifact_annotation"]["edf_blink_candidates"]["enabled"] = True
    config.write_text(yaml.safe_dump(enabled))
    def no_muscle_recompute(*args, **kwargs):
        raise AssertionError("Candidate-only refresh must not recompute muscle scores")
    with monkeypatch.context() as scoped:
        scoped.setattr(artifacts, "annotate_muscle_zscore", no_muscle_recompute)
        assert artifacts.main(run + ["--processing-pass", "review", "--resume"]) == 0
    refreshed = json.loads((report_dir / "03_B_ArtefactAnnotation_review_report.json").read_text())
    assert refreshed["runs"][0]["execution_status"] == "resumed_existing"
    assert refreshed["runs"][0]["edf_blink_candidates"]["status"] == "unavailable_or_unmatched"
    assert refreshed["runs"][0]["edf_blink_candidates"]["review_status"] == "review_required"
    with pytest.raises(SystemExit):
        review_cleaning.main(accepted + ["--stage", "4B", "--filter-edge-margin-sec", "2"])
    config.write_text(yaml.safe_dump(cfg))
    # Save an explicit reviewed-empty ocular artifact set, without opening a GUI.
    from pipeline_config import load_and_resolve_config
    resolved, *_ = load_and_resolve_config(config, root / "params/sub-04.yaml", tuning)
    record = AcquisitionRunRecord(str(source), source.name, source.stem, 1, "sub04_220526/meg/" + source.name)
    sss = derivatives / "preprocessing/sub-04/ses-01/sss" / (source.stem + "_raw_sss.fif")
    header = mne.io.read_raw_fif(sss, preload=False, verbose="ERROR")
    ann_path = manual_path(resolved, record); ann_path.parent.mkdir(parents=True)
    mne.Annotations([], [], [], orig_time=header.annotations.orig_time).save(ann_path)
    receipt = {"status": "reviewed", "coordinates": coordinates(header), "sss_input": file_stamp(sss),
               "annotation_sha256": file_hash(ann_path), "annotation_count": 0, "basis": "synthetic"}
    update_tuning_file(path=tuning, raw_root=tmp_path / "raw", subject_id="04", section="artifact_annotation",
                       values={"manual_eye_reviews": {source.name: receipt}})
    assert artifacts.main(run + ["--processing-pass", "review", "--overwrite"]) == 0
    with pytest.raises(SystemExit):
        review_cleaning.main(accepted + ["--stage", "4B", "--filter-edge-margin-sec", "165"])
    assert review_cleaning.main(accepted + ["--stage", "4B"]) == 0
    assert yaml.safe_load(tuning.read_text()).get("preprocessing", {}).get("filter_edge_margin_sec") is None
    assert ica.main(run + ["--processing-pass", "apply"]) == 2  # No reviewed model yet.
    def synthetic_component_figures(**kwargs):
        """Supply tiny component PNGs; this gate test does not test ICA plots."""
        figure_dir = kwargs["figure_dir"]
        names = ["03_C_ICA_review_sources.png", "03_C_ICA_review_topographies_01.png"]
        names += [f"03_C_ICA_review_ICA{i:03d}_properties_01.png"
                  for i in range(int(kwargs["model"].n_components_))]
        paths = []
        for name in names:
            figure, axis = plt.subplots(figsize=(1, 1))
            axis.plot([0, 1], [0, 1])
            path = figure_dir / name
            figure.savefig(path, dpi=30)
            plt.close(figure)
            paths.append(path)
        return paths
    monkeypatch.setattr(ica, "generate_component_diagnostics", synthetic_component_figures)
    def synthetic_native_figures(cfg, model, records, inputs, figure_dir, model_sha):
        """Keep this model/application gate test independent of Stage 2/3 fixtures.

        Native extraction, exact accepted anchors and both real plot types are
        exercised in test_final_ica_review; here only diagnostic coverage matters.
        """
        from shutil import copy2
        paths = []
        for i in range(model.n_components_):
            for kind in ('continuous', 'task'):
                path = figure_dir / f'03_C_ICA_review_ICA{i:03}_{kind}_01.png'
                copy2(figure_dir / f'03_C_ICA_review_ICA{i:03}_properties_01.png', path)
                paths.append(path)
        manifest = {'components_complete': list(range(model.n_components_))}
        (figure_dir / '03_C_ICA_review_native_diagnostics.json').write_text(json.dumps(manifest))
        return paths, manifest
    import final_ica_native_diagnostics
    monkeypatch.setattr(final_ica_native_diagnostics, 'generate_native_diagnostics', synthetic_native_figures)
    assert ica.main(run + ["--processing-pass", "review"]) == 0
    assert review_cleaning.main(accepted + ["--stage", "4C", "--exclude"]) == 0
    assert ica.main(run + ["--processing-pass", "apply", "--overwrite"]) == 0
    assert ica.main(run + ["--processing-pass", "apply", "--resume"]) == 0
    final = sss.with_name(source.stem + "_ica.fif")
    check_preservation(header, mne.io.read_raw_fif(final, preload=False, verbose="ERROR"), annotations=False)
    # A final-only change requires refreshed 4B evidence/output application,
    # while the accepted SSS and the compatible participant model stay intact.
    sss_report = (report_dir / "03_A_MaxFilt_review_report.json").read_bytes()
    model_path = derivatives / "preprocessing/sub-04/ses-01/ica_models/sub-04_ses-01_ica_prefit.fif"
    model_bytes = model_path.read_bytes()
    annotation_stamp = file_stamp(sss.with_name(source.stem + "_ann.fif"))
    cfg["preprocessing"]["spectrum_fit_window"] = "4s"
    config.write_text(yaml.safe_dump(cfg))
    assert ica.main(run + ["--processing-pass", "apply", "--resume"]) == 2
    assert artifacts.main(run + ["--processing-pass", "review", "--resume"]) == 0
    assert review_cleaning.main(accepted + ["--stage", "4B"]) == 0
    assert ica.main(run + ["--processing-pass", "apply", "--resume"]) == 2  # Stale final application.
    def no_refit(*args, **kwargs):
        raise AssertionError("Final-only policy changes must reuse the compatible model")
    with monkeypatch.context() as scoped:
        scoped.setattr(ica.ICA, "fit", no_refit)
        assert ica.main(run + ["--processing-pass", "apply", "--overwrite"]) == 0
    assert model_path.read_bytes() == model_bytes
    assert file_stamp(sss.with_name(source.stem + "_ann.fif")) == annotation_stamp
    assert (report_dir / "03_A_MaxFilt_review_report.json").read_bytes() == sss_report
    changed = yaml.safe_load(tuning.read_text()); changed["ica"]["fit_l_freq_hz"] = 2.
    tuning.write_text(yaml.safe_dump(changed))
    assert ica.main(run + ["--processing-pass", "apply", "--resume"]) == 2


def test_runner_preserves_accepted_stage1_snapshots(tmp_path, monkeypatch):
    """Stage 4 adds its own execution snapshots and never refreshes Stage 1."""
    import yaml
    from types import SimpleNamespace
    import run_pipeline
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "params/pipeline_default.yaml").read_text())
    (tmp_path / "raw/sub04_220526").mkdir(parents=True)
    cfg["project"].update(data_root=str(tmp_path / "raw"), derivatives_root=str(tmp_path / "derivatives"))
    config = tmp_path / "config.yaml"; config.write_text(yaml.safe_dump(cfg))
    reports = tmp_path / "derivatives/reports/sub-04/ses-01"; reports.mkdir(parents=True)
    snapshots = [reports / name for name in ["source_manifest.json", "dataset_inventory.json", "pipeline_run.json", "resolved_config.json"]]
    for path in snapshots:
        path.write_text('{"accepted_stage1": "preserve exact bytes"}')
    before = [p.read_bytes() for p in snapshots]
    monkeypatch.setattr(sys, "argv", ["runner", "--config", str(config), "--subject-config", str(root / "params/sub-04.yaml"), "--stages", "maxwell_filter", "--dry-run"])
    assert run_pipeline.main() == 0
    assert [p.read_bytes() for p in snapshots] == before
    assert (reports / "stage4_pipeline_run.json").is_file()
    assert (reports / "stage4_resolved_config.json").is_file()


def test_preliminary_sources_stay_on_the_native_display_grid(monkeypatch):
    """Display/reselect real cached sources without a local fit or MEG mutation."""
    from manual_eye_annotations import add_preliminary_sources
    # Both real Vectorview frontal regions are required for the bilateral aid.
    original = raw_example(duration=30.)
    names = ["MEG0121", "MEG0311", "MEG0321", "MEG0811", "MEG0911", "MEG0921"]
    data = np.random.default_rng(21).normal(size=(6, original.n_times)) * 1e-12
    raw = mne.io.RawArray(data, mne.create_info(names, 400., "mag"), first_samp=4000, verbose="ERROR")
    display = raw.copy()
    fit, summary = prepare_fit_copy(raw, FIT)
    model = mne.preprocessing.ICA(n_components=4, random_state=97, max_iter="auto")
    model.fit(fit)
    sources = model.get_sources(raw.copy().filter(1., 40., verbose="ERROR"))
    monkeypatch.setattr(mne.preprocessing.ICA, "plot_components", lambda *a, **k: None)
    monkeypatch.setattr(mne.io.BaseRaw, "plot", lambda *a, **k: None)
    monkeypatch.setattr(mne.preprocessing.ICA, "fit", lambda *a, **k: pytest.fail("Local review must not fit ICA"))
    answers = iter(["ICA000", "NONE", "1", "NONE"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    state = {"model": model, "native_sources": sources, "metadata": {
        "identity": {"channels": names}, "fingerprint": "synthetic_cache", "model_sha256": "synthetic_model",
        "converged": True, "fit_summary": summary}}
    receipt = add_preliminary_sources(display, {"ica": FIT}, state=state)
    assert receipt["purpose"] == "supplementary_annotation_only"
    assert "AID000" in display.ch_names and "AID000" not in raw.ch_names
    model = state["model"]
    display.annotations.append(display.first_time + 2., .5, "BAD_blink")
    before = annotation_rows(display)
    # Revisit the same fitted decomposition; preserve edits and replace the aid.
    receipt = add_preliminary_sources(display, {"ica": FIT}, state=state)
    assert state["model"] is model
    assert "AID001" in display.ch_names and "AID000" not in display.ch_names
    assert receipt["selected_components"] == [1]
    assert annotation_rows(display) == before
    expected = sources.get_data(picks="ICA001")
    np.testing.assert_allclose(display.get_data(picks="AID001"), expected / np.max(abs(expected)))
    assert display.n_times == raw.n_times and display.first_samp == raw.first_samp
    assert display.info["sfreq"] == raw.info["sfreq"]
    np.testing.assert_array_equal(display.get_data(picks=raw.ch_names), raw.get_data())


def test_annotation_csv_matches_saved_fif_at_large_acquisition_times(tmp_path):
    """FIF float32 annotation endpoints must not cause false resume failures."""
    import importlib.util
    path = Path(__file__).resolve().parents[1] / "code/preprocessing_basic/03_B_ArtefactAnnotation.py"
    spec = importlib.util.spec_from_file_location("annotation_precision_test", path)
    stage = importlib.util.module_from_spec(spec); spec.loader.exec_module(stage)
    example = raw_example(date=False, duration=8.)
    raw = mne.io.RawArray(example.get_data(), example.info, first_samp=1397000, verbose="ERROR")
    raw.set_annotations(mne.Annotations([.75], [.1], ["acquisition_note"]))
    source = tmp_path / "source_raw_sss.fif"; raw.save(source, overwrite=True, verbose="ERROR")
    ocular = mne.Annotations([raw.first_time + 2.12345], [.234567], ["BAD_blink"])
    record = AcquisitionRunRecord(str(source), "source.fif", "source", 1, "source.fif")
    settings = {"muscle_channel_type": "mag", "muscle_min_length_good_sec": .2,
                "muscle_filter_freq_hz": [110., 140.], "eog_mode": "disabled", "eog_channel_names": []}
    output = tmp_path / "source_ann.fif"; csv_path = tmp_path / "source_ann.csv"
    stage._annotate_one_run(record=record, input_path=source, output_fif=output, output_csv=csv_path,
        artifact_cfg=settings, initial_threshold=10., interactive=False, show=False,
        muscle_figure_path=tmp_path / "score.png", eog_figure_path=tmp_path / "legacy.png",
        overwrite=True, manual_annotations=ocular)
    stage._validate_resumable_annotation(record=record, input_path=source, output_fif=output,
        output_csv=csv_path, threshold=10., artifact_cfg=settings, raw_root=tmp_path / "raw")


@pytest.mark.parametrize("date", [True, False])
@pytest.mark.parametrize("reviewed", ["edited", "empty", "legacy"])
def test_muscle_qc_replaces_automatic_intervals_including_reviewed_zero(tmp_path, monkeypatch, date, reviewed):
    """Exercise real FIF/CSV rebuilding; proposed labels must not resurrect."""
    import importlib.util
    from continuous_cleaning import digest, muscle_parameters
    path = Path(__file__).resolve().parents[1] / "code/preprocessing_basic/03_B_ArtefactAnnotation.py"
    spec = importlib.util.spec_from_file_location("muscle_qc_test", path)
    stage = importlib.util.module_from_spec(spec); spec.loader.exec_module(stage)
    raw = raw_example(date=date, duration=8.)
    raw.set_annotations(mne.Annotations([.2, .4], [.1, .1], ["acquisition_note", "BAD_muscle"]))
    source = tmp_path / "source_raw_sss.fif"; raw.save(source, overwrite=True, verbose="ERROR")
    record = AcquisitionRunRecord(str(source), "source.fif", "source", 1, "source.fif")
    settings = {"muscle_channel_type": "mag", "muscle_min_length_good_sec": .2,
                "muscle_filter_freq_hz": [110., 140.], "eog_mode": "disabled", "eog_channel_names": []}
    manual = mne.Annotations([raw.first_time + 5.], [.1], ["BAD_eye_movement"], orig_time=raw.annotations.orig_time)
    manual.append(raw.first_time + 5.5, .2, "BAD_eye_movement_auto")
    if reviewed == "edited":
        manual.append([raw.first_time + 1.5, raw.first_time + 6.], [1.2, .3], ["BAD_muscle"] * 2)
    receipt = {"status": "reviewed"}
    if reviewed != "legacy":
        receipt["muscle_review"] = {"status": "reviewed", "settings_sha256": digest(muscle_parameters(settings, 10.))}

    def proposals(data, **unused):
        # Scores must not mask themselves with newly reviewed muscle periods.
        assert sum(label == "BAD_muscle" for label in data.annotations.description) == 1
        offset = data.first_time if data.info["meas_date"] is not None else 0.
        return mne.Annotations([offset + 1., offset + 3.], [.1, .2], ["BAD_muscle"] * 2,
                               orig_time=data.annotations.orig_time), np.zeros(data.n_times)

    monkeypatch.setattr(stage, "annotate_muscle_zscore", proposals)
    output = tmp_path / "source_ann.fif"; csv_path = tmp_path / "source_ann.csv"
    report, _, _ = stage._annotate_one_run(record=record, input_path=source, output_fif=output, output_csv=csv_path,
        artifact_cfg=settings, initial_threshold=10., interactive=False, show=False,
        muscle_figure_path=tmp_path / "score.png", eog_figure_path=tmp_path / "legacy.png",
        overwrite=True, manual_annotations=manual, manual_receipt=receipt)
    result = mne.io.read_raw_fif(output, preload=False, verbose="ERROR")
    muscle = result.annotations[[i for i, label in enumerate(result.annotations.description) if label == "BAD_muscle"]]
    expected = {"edited": [1.5, 6.], "empty": [], "legacy": [.4, 1., 3.]}[reviewed]
    np.testing.assert_allclose(muscle.onset - raw.first_time, expected, atol=2e-6)
    assert {"acquisition_note", "BAD_eye_movement", "BAD_eye_movement_auto"}.issubset(result.annotations.description)
    assert report["automatic_muscle_annotation_count"] == 2
    assert report["muscle_annotation_count"] == len(expected)
    comparison = raw.copy()
    if reviewed != "legacy":
        comparison.annotations.delete([i for i, label in enumerate(comparison.annotations.description) if label == "BAD_muscle"])
    check_preservation(comparison, result, annotations=False)
    stage._validate_resumable_annotation(record=record, input_path=source, output_fif=output,
        output_csv=csv_path, threshold=10., artifact_cfg=settings, raw_root=tmp_path / "raw",
        expected_signature=stage.artifact_signature(source, 10., settings, receipt), reviewed_muscle=reviewed != "legacy")
    if reviewed != "legacy":
        settings["muscle_min_length_good_sec"] = .3
        with pytest.raises(RuntimeError, match="settings changed"):
            stage._annotate_one_run(record=record, input_path=source, output_fif=output, output_csv=csv_path,
                artifact_cfg=settings, initial_threshold=10., interactive=False, show=False,
                muscle_figure_path=tmp_path / "score.png", eog_figure_path=tmp_path / "legacy.png",
                overwrite=True, manual_annotations=manual, manual_receipt=receipt)


def test_historical_point01_filter_removes_lines_and_preserves_native_stream():
    """Keep historical FIR exposure testable without selecting it for production."""
    raw = raw_example(duration=700.)
    times = np.arange(raw.n_times) / raw.info["sfreq"]
    raw._data[0] += 2e-11 * (np.sin(2 * np.pi * 50 * times) + np.sin(2 * np.pi * 100 * times))
    raw.set_annotations(mne.Annotations([350.], [1.], ["BAD_sensor_fault"]))
    raw._data[:3, 140000:140400] += 1e3  # Must never be filtered into good neighbours.
    before = raw.copy()
    summary = filter_continuous(raw, {"highpass_hz": .01, "filter_edge_margin_sec": 165.,
                                     "notch_freqs_hz": [50., 100.], "notch_method": "spectrum_fit"})
    check_preservation(before, raw, annotations=False)
    assert summary["seconds_outside_edge_exposure"] > 0
    a, b = 166 * 400, 184 * 400
    for frequency in (50, 100):
        wave = np.exp(-2j * np.pi * frequency * times[a:b])
        original = abs(np.mean(before.get_data(picks=[0], start=a, stop=b)[0] * wave))
        cleaned = abs(np.mean(raw.get_data(picks=[0], start=a, stop=b)[0] * wave))
        assert cleaned < .1 * original
