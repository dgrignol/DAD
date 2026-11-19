"""Scientific-boundary tests for Stage 4A SQUID screening and review.

Usage: MPLBACKEND=Agg python -B -m pytest -p no:cacheprovider meg/tests/test_squid_jumps.py -q

All raw FIFs, review receipts, annotation sentinels and calibration placeholders
are temporary synthetic fixtures. No cohort data, active P05 browser, scheduler
or saved P04/P05 annotation is accessed. Tests cover detection rather than a
particular count of false positives, native time, immutable decisions, locked
sibling-save preservation, operative identities and real MNE SSS skipping.
"""

import copy
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import sys

import mne
import numpy as np
import pytest
import yaml

CODE = Path(__file__).resolve().parents[1] / "code/preprocessing_basic"
sys.path.insert(0, str(CODE))
import squid_jumps as squid
from continuous_cleaning import annotation_rows, coordinates, file_hash, good_stretches, stage_configuration
from manual_eye_annotations import eye_annotations
from preprocessing_stage_utils import AcquisitionRunRecord, update_tuning_file


def example(tmp_path, *, date=False):
    """Make source/auxiliary data and minimal portable science configuration."""
    rng = np.random.default_rng(42)
    data = rng.normal(scale=1e-13, size=(4, 10000))
    data[1] *= 100  # A separate physical sensor scale must not set its threshold.
    data[2] = np.arange(10000) % 200
    info = mne.create_info(["MEG0111", "MEG0122", "STI101", "MISC001"], 1000., ["mag", "grad", "stim", "misc"])
    raw = mne.io.RawArray(data, info, first_samp=12000, verbose="ERROR")
    if date:
        raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    directory = tmp_path / "raw"; directory.mkdir()
    source = directory / "sub04_run01_synthetic.fif"
    raw.save(source, overwrite=True, verbose="ERROR")
    resource = tmp_path / "calibration"; resource.write_bytes(b"synthetic")
    cfg = {"subject": {"id": "04"}, "maxwell_filter": {"calibration_file": str(resource), "crosstalk_file": str(resource)},
           "preprocessing": {}, "squid_jumps": {"enabled": True, "threshold": 10., "block_sec": .5,
           "tuning_path": str(tmp_path / "03_tuning_choices.yaml")},
           "resolved_paths": {"data_root": str(directory), "derivatives_preprocessing_dir": str(tmp_path / "preprocessing")}}
    record = AcquisitionRunRecord(str(source), source.name, source.stem, 1, source.name)
    return raw, cfg, record


@pytest.mark.parametrize("block", [.5, .997, 30.])
def test_persistent_and_decaying_steps_across_chunks_keep_single_sensor(tmp_path, block):
    raw, cfg, _ = example(tmp_path)
    raw._data[0, 1000:] += 1e-10
    raw._data[0, 4500:] += 1e-10*np.exp(-np.arange(5500)/100.)
    raw._data[1, 6000:] += 1e-6  # Accepted bad channel must never contribute.
    result = squid.detect(raw, ["MEG0122"], {**squid.parameters(cfg), "block_sec": block})
    peaks = [e["peak_sample"]-raw.first_samp for e in result["candidates"]]
    assert any(abs(p-1000) <= 4 for p in peaks)
    assert any(abs(p-4500) <= 4 for p in peaks)
    assert result["screened_channels"] == ["MEG0111"]
    assert all(e["channels"] == ["MEG0111"] for e in result["candidates"])


def test_true_acquisition_boundary_is_not_a_jump_and_flat_is_unscored(tmp_path):
    raw, cfg, _ = example(tmp_path)
    raw._data[0, 3000:] += 1e-9
    raw._data[1] = 0
    raw.set_annotations(mne.Annotations([3.], [0.], ["EDGE boundary"]))
    result = squid.detect(raw, [], squid.parameters(cfg))
    assert not any(abs(e["peak_sample"]-raw.first_samp-3000) < 10 for e in result["candidates"])
    assert result["unscored_channels"] == [{"channel": "MEG0122", "reason": "zero_robust_scale"}]


def test_small_positive_median_filter_scale_does_not_flag_quantized_noise(tmp_path):
    """Native high-frequency/quantized noise needs a raw-MAD scale floor.

    The first real pilot had positive filtered MADs around 6% of raw MAD;
    treating only exact-zero MAD as special inflated ordinary fluctuations.
    Retain true isolated steps without manufacturing broad noisy candidates.
    """
    raw, cfg, _ = example(tmp_path)
    rng = np.random.default_rng(83)
    t = np.arange(raw.n_times)/raw.info["sfreq"]
    noise = np.round((np.sin(2*np.pi*240*t) + .1*rng.normal(size=raw.n_times))/.01)*.01*1e-11
    raw._data[0] = noise
    result = squid.detect(raw, ["MEG0122"], squid.parameters(cfg))
    assert result["candidates"] == []
    raw._data[0,3000:] += 5e-10
    raw._data[0,7000:] -= 5e-10
    result = squid.detect(raw, ["MEG0122"], squid.parameters(cfg))
    assert any(abs(e["peak_sample"]-raw.first_samp-3000) <= 5 for e in result["candidates"])
    assert any(abs(e["peak_sample"]-raw.first_samp-7000) <= 5 for e in result["candidates"])


def test_dense_cross_channel_edges_never_chain_into_a_recording_long_event(tmp_path):
    """Staggered sensor events remain local even if their gaps overlap."""
    raw, cfg, _ = example(tmp_path)
    # Alternate edges every 80 ms between two sensors. Each channel's next
    # edge is outside 100 ms; naive cross-channel transitive merging chains
    # all of them into one long event even with correct channel thresholds.
    for i, sample in enumerate(range(1000,9000,80)):
        raw._data[i % 2,sample:] += 1e-8
    result = squid.detect(raw, [], squid.parameters(cfg))
    assert len(result["candidates"]) > 10
    assert max(e["last_sample"]-e["first_sample"] for e in result["candidates"]) <= 300
    assert all(e["stop_sample"]-e["start_sample"] <= 1301 for e in result["candidates"])


@pytest.mark.parametrize("date", [False, True])
def test_native_annotation_roundtrip_masks_epochs_without_changing_samples(tmp_path, date):
    raw, _, _ = example(tmp_path, date=date)
    original = raw.get_data().copy()
    squid.attach_intervals(raw, [[14000, 14500]])
    squid.attach_intervals(raw, [[14000, 14500]])
    assert list(raw.annotations.description) == [squid.LABEL]
    assert good_stretches(raw) == [(0, 2000), (2500, 10000)]
    target = tmp_path / "output_raw.fif"; raw.save(target, verbose="ERROR")
    saved = mne.io.read_raw_fif(target, preload=True, verbose="ERROR")
    assert coordinates(raw) == coordinates(saved)
    assert saved.annotations.onset[0] == pytest.approx(14.)
    events = np.array([[14050, 0, 1], [16000, 0, 1]])
    epochs = mne.Epochs(saved, events, tmin=0., tmax=.1, baseline=None, preload=True, verbose="ERROR")
    assert len(epochs) == 1 and squid.LABEL in epochs.drop_log[0]
    np.testing.assert_array_equal(raw.get_data(), original)


def test_immutable_review_and_sibling_manual_save_survive_stale_gui_config(tmp_path):
    raw, cfg, record = example(tmp_path)
    path, proposal = squid.prepare(cfg, record, [])
    # Model a P04/P05 eye save completed after the SQUID GUI loaded its cfg.
    sentinel = tmp_path / "manual_eye" / "p05run02-annot.fif"
    sentinel.parent.mkdir(); sentinel.write_bytes(b"active eye review: never touch")
    update_tuning_file(path=cfg["squid_jumps"]["tuning_path"], raw_root=cfg["resolved_paths"]["data_root"], subject_id="04",
        section="artifact_annotation", mapping_key="manual_eye_reviews", values={"p04_all": {"status": "reviewed"}, "p05run01": {"status": "reviewed"}})
    before = sentinel.read_bytes()
    action = {"intervals": [], "bad_channels": [], "channel_scope": "union"}
    saved = squid.save_review(cfg, record, [], proposal, path, action, reviewer="human", reason="Reviewed zero; flat controls inspected", dismissed=[])
    latest = yaml.safe_load(Path(cfg["squid_jumps"]["tuning_path"]).read_text())
    assert latest["artifact_annotation"]["manual_eye_reviews"] == {"p04_all": {"status": "reviewed"}, "p05run01": {"status": "reviewed"}}
    assert sentinel.read_bytes() == before
    cfg["squid_jumps"]["reviews"] = latest["squid_jumps"]["reviews"]
    assert squid.load_review(cfg, record, []) == saved
    old = root_files = sorted(squid.root(cfg).glob("*-review.json"))
    original = {p.name: p.read_bytes() for p in root_files}
    squid.save_review(cfg, record, [], proposal, path, action, reviewer="human", reason="Retained zero", dismissed=[])
    assert all((squid.root(cfg)/name).read_bytes() == data for name, data in original.items())
    assert len(list(squid.root(cfg).glob("*-review.json"))) == len(old)+1
    # A modified immutable artifact cannot be accepted by merely editing YAML.
    review_path = squid.root(cfg) / cfg["squid_jumps"]["reviews"][record.source_filename]["review_file"]
    review_path.write_text("{}")
    with pytest.raises(RuntimeError, match="changed"):
        squid.load_review(cfg, record, [])


def test_empty_reviews_and_detector_tuning_leave_sss_identity_unchanged(tmp_path):
    _, cfg, record = example(tmp_path)
    original = stage_configuration(cfg, "03_A_MaxFilt")
    cfg["squid_jumps"]["threshold"] = 8.
    cfg["squid_jumps"]["reviews"] = {record.source_filename: {"status": "reviewed", "application":
        {"intervals": [], "bad_channels": [], "channel_scope": "union"}}}
    assert stage_configuration(cfg, "03_A_MaxFilt") == original
    cfg["squid_jumps"]["reviews"][record.source_filename]["application"]["intervals"] = [[14000, 14500]]
    assert stage_configuration(cfg, "03_A_MaxFilt") != original
    assert squid.run_action(cfg, "another_run.fif") == {"intervals": [], "bad_channels": []}
    cfg["squid_jumps"]["reviews"][record.source_filename]["application"].update(bad_channels=["MEG0111"])
    assert squid.run_action(cfg, "another_run.fif")["bad_channels"] == ["MEG0111"]
    cfg["squid_jumps"]["reviews"][record.source_filename]["application"]["channel_scope"] = "run"
    assert squid.run_action(cfg, "another_run.fif")["bad_channels"] == []


def test_stage4b_snapshot_cannot_duplicate_or_replace_squid_owner(tmp_path):
    raw, _, _ = example(tmp_path)
    squid.attach_intervals(raw, [[14000, 14500]])
    raw.annotations.append(16., .2, "BAD_blink")
    raw.annotations.append(17., .3, "BAD_muscle")
    original = annotation_rows(raw)
    editable = eye_annotations(raw)
    assert set(editable.description) == {"BAD_blink", "BAD_muscle"}
    assert annotation_rows(raw) == original
    assert list(squid.without_owned_context(raw.annotations).description) == list(editable.description)


def test_real_mne_sss_never_mixes_skipped_artifact_into_good_samples():
    """Use real SSS on synthetic radial coils, not a mocked filter call."""
    nchan, samples = 40, 1000
    rng = np.random.default_rng(8)
    info = mne.create_info([f"MEG{i:03d}1" for i in range(nchan)], 200., "mag")
    # Arbitrary radial coils cover a sphere; low expansion orders keep this
    # synthetic geometry well determined without claiming a MEGIN calibration.
    for ch in info["chs"]:
        z = rng.normal(size=3); z /= np.linalg.norm(z)
        x = np.cross(z, [0., 0., 1.]); x /= np.linalg.norm(x)
        y = np.cross(z, x)
        # Vary radius: an exactly spherical radial array cannot distinguish
        # internal/external fields of the same angular order.
        ch["loc"] = np.r_[(.10+.04*rng.random())*z, x, y, z]
    raw = mne.io.RawArray(rng.normal(scale=1e-13, size=(nchan, samples)), info, first_samp=2000, verbose="ERROR")
    squid.attach_intervals(raw, [[2400, 2500]])
    damaged = raw.copy(); damaged._data[0, 400:500] += 1e-5
    kwargs = dict(origin=(0., 0., 0.), int_order=2, ext_order=1, regularize=None,
                  coord_frame="meg", skip_by_annotation=("edge", "bad_acq_skip", squid.LABEL), verbose="ERROR")
    first = squid.maxwell_filter_with_jumps(raw, **kwargs)
    second = squid.maxwell_filter_with_jumps(damaged, **kwargs)
    for start, stop in good_stretches(raw):
        np.testing.assert_array_equal(first.get_data(start=start, stop=stop), second.get_data(start=start, stop=stop))
    np.testing.assert_array_equal(second.get_data()[:,400:500], damaged.get_data()[:,400:500])
    assert coordinates(raw) == coordinates(first)


def test_browser_view_retains_deletion_and_saves_reviewed_zero(tmp_path, monkeypatch):
    """A deleted proposal must not be reseeded on VIEW or later reopening."""
    import review_squid_jumps as ui
    raw, cfg, record = example(tmp_path)
    raw._data[0, 5000:] += 1e-9
    raw.save(record.source_path, overwrite=True, verbose="ERROR")
    path, proposal = squid.prepare(cfg, record, ["MEG0122"])
    assert proposal["candidates"]
    views = []
    def browser(display, **kwargs):
        assert kwargs["remove_dc"] is True and all(value > 0 for value in kwargs["scalings"].values())
        views.append(sum(display.annotations.description == squid.CANDIDATE_LABEL))
        display.annotations.delete(np.flatnonzero(display.annotations.description == squid.CANDIDATE_LABEL))
    monkeypatch.setattr(mne.io.BaseRaw, "plot", browser)
    answers = iter(["VIEW", "SAVE", "operator", "Dismissed after raw inspection", "yes"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    saved = ui.review(cfg, record, ["MEG0122"], path, proposal)
    assert views[0] > 0 and views[1] == 0
    assert saved["application"]["intervals"] == []
    assert saved["dismissed_candidate_ids"]
    cfg["squid_jumps"]["reviews"] = yaml.safe_load(Path(cfg["squid_jumps"]["tuning_path"]).read_text())["squid_jumps"]["reviews"]
    answers = iter(["DISCARD"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    assert ui.review(cfg, record, ["MEG0122"], path, proposal) is None
    assert views[-1] == 0
    cfg["squid_jumps"]["threshold"] = 8.
    new_path, new_proposal = squid.prepare(cfg, record, ["MEG0122"])
    with pytest.raises(RuntimeError, match="reconcile"):
        ui.review(cfg, record, ["MEG0122"], new_path, new_proposal)
    answers = iter(["SAVE", "operator", "Rechecked new sensitivity; retain zero", "yes"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    reconciled = ui.review(cfg, record, ["MEG0122"], new_path, new_proposal, reconcile=True)
    assert views[-1] == 0 and reconciled["application"]["intervals"] == []
    assert reconciled["fingerprint"] != saved["fingerprint"]


def test_candidate_sensor_comes_first_without_losing_unflagged_controls(tmp_path):
    """Events on a later sensor cannot hide behind MNE's first channel page."""
    from review_squid_jumps import candidate_order
    raw, _, _ = example(tmp_path)
    proposal = {"screened_channels": ["MEG0111", "MEG0122"],
                "candidates": [{"channels": ["MEG0122"]}, {"channels": ["MEG0111"]}]}
    assert candidate_order(raw, proposal, 0) == [1, 0]
    assert candidate_order(raw, proposal, 1) == [0, 1]
    proposal["candidates"] = []
    assert candidate_order(raw, proposal, 0) == [0, 1]


def test_local_window_scale_ignores_remote_drift_and_acquired_dc(tmp_path):
    """Large offsets elsewhere cannot flatten the current candidate display."""
    from review_squid_jumps import window_scalings
    raw, _, _ = example(tmp_path)
    names = ["MEG0111", "MEG0122"]
    expected = window_scalings(raw,names,4.,2.)
    raw._data[0,:1000] += 1e-6
    raw._data[1,:1000] -= 1e-4
    raw._data[0] += 1e-9
    raw._data[1] -= 1e-7
    actual = window_scalings(raw,names,4.,2.)
    np.testing.assert_allclose(list(actual.values()),list(expected.values()),rtol=1e-8,atol=0)
    assert all(value > 0 for value in actual.values())


def test_display_proposal_union_preserves_original_ids_and_separate_periods():
    """One burst needs one editable shade, rather than stacked opaque spans."""
    from review_squid_jumps import proposal_windows
    proposal={"candidates":[{"id":1,"start_sample":100,"stop_sample":200},
                            {"id":2,"start_sample":150,"stop_sample":250},
                            {"id":3,"start_sample":400,"stop_sample":450}]}
    original=copy.deepcopy(proposal)
    assert proposal_windows(proposal)==[[100,250],[400,450]]
    assert proposal == original


def test_inspector_axes_label_flagged_sensors_and_leave_data_decisions_untouched(tmp_path):
    """Per-sensor physical close-ups cannot edit source samples or annotations."""
    from review_squid_jumps import inspection_figures
    import matplotlib.pyplot as plt
    raw, cfg, record = example(tmp_path)
    raw._data[0,5000:] += 1e-10
    raw.save(record.source_path,overwrite=True,verbose="ERROR")
    _, proposal = squid.prepare(cfg,record,[])
    original = raw.get_data().copy()
    original_annotations = annotation_rows(raw)
    fig, metrics = inspection_figures(raw,proposal,0)
    assert len(fig) == 1 and len(fig[0].axes) == 3*len(metrics)
    assert all(ch in fig[0].axes[3*i].get_title() for i,ch in enumerate(row["channel"] for row in metrics))
    assert "6 s context" in fig[0].axes[0].get_title()
    assert "300 ms" in fig[0].axes[1].get_title()
    assert "native sample" in fig[0].axes[2].get_ylabel()
    np.testing.assert_array_equal(raw.get_data(),original)
    assert annotation_rows(raw) == original_annotations
    for page in fig:
        plt.close(page)


def test_inspect_cli_reuses_candidates_without_any_preparation_or_review_write(tmp_path,monkeypatch):
    """A second inspector can coexist with a live unsaved editable review."""
    import review_squid_jumps as ui
    raw,cfg,record = example(tmp_path)
    raw._data[0,5000:] += 1e-10
    raw.save(record.source_path,overwrite=True,verbose="ERROR")
    _,proposal = squid.prepare(cfg,record,[])
    sentinel = tmp_path / "active-review.json"; sentinel.write_text("unsaved session is outside inspector")
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    monkeypatch.setattr(ui,"load_and_resolve_config",lambda *a,**k: (cfg,))
    monkeypatch.setattr(ui,"discover_raw_run_records",lambda _: [record])
    def cannot_write(*a,**k):
        raise AssertionError("Inspection must never prepare candidates, channel evidence or review receipts")
    monkeypatch.setattr(ui,"channel_lists",cannot_write)
    monkeypatch.setattr(ui,"prepare",cannot_write)
    monkeypatch.setattr(ui,"save_review",cannot_write)
    seen=[]
    monkeypatch.setattr(ui,"inspect_candidate",lambda raw,proposal,index,**kwargs: seen.append((index,kwargs)))
    assert ui.main(["--mode","inspect","--config","unused","--subject-config","unused",
        "--tuning-config",str(tmp_path/"03_tuning_choices.yaml"),"--source",record.source_filename,"--candidate","1"]) == 0
    assert seen == [(0,{"context":False})]
    assert before == {p:p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}


def test_enabled_stage4a_stops_before_sss_then_applies_reviewed_windows(tmp_path, monkeypatch):
    """Run the real 4A CLI path with synthetic I/O and a bounded SSS stand-in."""
    from pipeline_config import load_and_resolve_config
    from preprocessing_stage_utils import discover_raw_run_records
    spec = importlib.util.spec_from_file_location("squid_stage4a_cli", CODE / "03_A_MaxFilt.py")
    stage = importlib.util.module_from_spec(spec); spec.loader.exec_module(stage)
    cfg = yaml.safe_load((CODE.parents[1] / "params/pipeline_default.yaml").read_text())
    raw_dir = tmp_path / "raw/sub04_synthetic/meg"; raw_dir.mkdir(parents=True)
    cfg["project"].update(data_root=str(tmp_path / "raw"), derivatives_root=str(tmp_path / "derivatives"))
    for key in ("calibration_file", "crosstalk_file"):
        resource = tmp_path / key; resource.write_text("synthetic")
        cfg["maxwell_filter"][key] = str(resource)
    source = raw_dir / "sub04_run01_220526.fif"
    rng = np.random.default_rng(9)
    raw = mne.io.RawArray(rng.normal(scale=1e-13, size=(3,10000)),
        mne.create_info(["MEG0111", "MEG0121", "STI101"], 1000., ["mag", "mag", "stim"]), first_samp=1000, verbose="ERROR")
    # MNE 1.11 no longer supplies a synthetic device/head identity by default;
    # this fixture explicitly represents the valid static transforms in DAD.
    raw.info["dev_head_t"] = mne.transforms.Transform("meg", "head", np.eye(4))
    raw._data[0,5000:] += 1e-9
    raw.save(source, verbose="ERROR")
    config = tmp_path / "config.yaml"; config.write_text(yaml.safe_dump(cfg))
    subject = CODE.parents[1] / "params/sub-04.yaml"
    common = ["--config", str(config), "--subject-config", str(subject)]
    monkeypatch.setattr(stage, "validate_auto_origin_digitization", lambda *a, **k: {})
    monkeypatch.setattr(stage, "_centrality_head_points", lambda _: np.vstack([np.eye(3), np.zeros(3)]))
    def detect_channels(record, **kwargs):
        return ["MEG0121"], {"input_header_bads": [], "auto_noisy_channels": ["MEG0121"],
            "auto_flat_channels": [], "detected_bad_channels": ["MEG0121"], "digitization_points": {}, "scores": {}}, raw.copy()
    monkeypatch.setattr(stage, "_detect_bad_channels", detect_channels)
    monkeypatch.setattr(stage, "save_heatmaps", lambda **kwargs: None)
    filtered_bounds = []
    def solver(segment, **kwargs):
        filtered_bounds.append((segment.first_samp, segment.first_samp+segment.n_times))
        assert squid.LABEL in kwargs["skip_by_annotation"]
        result = segment.copy(); result._data[:2] *= .5
        return result
    monkeypatch.setattr(mne.preprocessing, "maxwell_filter", solver)
    answers = iter(["", "", "", "yes"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    assert stage.main(common + ["--interactive", "--no-display"]) == 2
    assert filtered_bounds == []  # Missing human decisions must precede SSS.
    resolved, *_ = load_and_resolve_config(config, subject, required_sections=("maxwell_filter",))
    tuning = Path(resolved["resolved_paths"]["derivatives_reports_dir"]) / "03_tuning_choices.yaml"
    resolved, *_ = load_and_resolve_config(config, subject, tuning, required_sections=("maxwell_filter",))
    resolved["squid_jumps"]["tuning_path"] = str(tuning)
    record = discover_raw_run_records(resolved)[0]
    path, proposal = squid.prepare(resolved, record, ["MEG0121"])
    squid.save_review(resolved, record, ["MEG0121"], proposal, path,
        {"intervals": [[5500,6500]], "bad_channels": [], "channel_scope": "union"},
        reviewer="human", reason="Synthetic jump window reviewed", dismissed=[])
    monkeypatch.setattr(stage, "build_run_qc", lambda **kwargs: {"source_filename": source.name, "plots": []})
    monkeypatch.setattr(stage, "render_stage_report", lambda **kwargs: tmp_path / "report.html")
    assert stage.main(common + ["--tuning-config", str(tuning)]) == 0
    assert filtered_bounds == [(1000,5500), (6500,11000)]
    output = Path(resolved["resolved_paths"]["derivatives_preprocessing_dir"]) / "sss" / f"{source.stem}_raw_sss.fif"
    saved = mne.io.read_raw_fif(output, preload=True, verbose="ERROR")
    assert list(saved.annotations.description) == [squid.LABEL]
    assert saved.info["bads"] == ["MEG0121"]
    np.testing.assert_allclose(saved.get_data()[:2,4500:5500], raw.get_data()[:2,4500:5500], rtol=1e-7)
    assert not (Path(resolved["resolved_paths"]["derivatives_preprocessing_dir"]) / "manual_eye").exists()
    # Synthetic scientific acceptance permits a no-op review/restart without
    # rewriting the accepted report, SSS stamp, or already reviewed 4B work.
    resolved, *_ = load_and_resolve_config(config, subject, tuning, required_sections=("maxwell_filter",))
    reports = Path(resolved["resolved_paths"]["derivatives_reports_dir"])
    report_path = reports / "03_A_MaxFilt_review_report.json"
    report = json.loads(report_path.read_text())
    update_tuning_file(path=tuning, raw_root=resolved["resolved_paths"]["data_root"], subject_id="04",
        section="stage4_review", values={"03_A_MaxFilt": {"status": "accepted",
            "configuration_sha256": report["configuration_sha256"], "evidence_sha256": report["evidence_sha256"]}})
    before_report, before_output = report_path.read_bytes(), squid.file_stamp(output)
    manual = Path(resolved["resolved_paths"]["derivatives_preprocessing_dir"]) / "manual_eye" / "already-reviewed-annot.fif"
    manual.parent.mkdir(); manual.write_bytes(b"P04/P05 saved manual BAD decisions")
    assert stage.main(common + ["--tuning-config", str(tuning), "--resume"]) == 0
    assert report_path.read_bytes() == before_report and squid.file_stamp(output) == before_output
    assert filtered_bounds == [(1000,5500), (6500,11000)]
    # A changed interval requires explicit rebuild; the old accepted evidence
    # and manually reviewed bytes survive while acceptance becomes stale.
    resolved["squid_jumps"]["tuning_path"] = str(tuning)
    squid.save_review(resolved, record, ["MEG0121"], proposal, path,
        {"intervals": [[5401,6602]], "bad_channels": [], "channel_scope": "union"},
        reviewer="human", reason="Adjusted to actual recovery", dismissed=[])
    assert stage.main(common + ["--tuning-config", str(tuning), "--resume"]) == 2
    assert squid.file_stamp(output) == before_output and report_path.read_bytes() == before_report
    assert stage.main(common + ["--tuning-config", str(tuning), "--overwrite"]) == 0
    assert manual.read_bytes() == b"P04/P05 saved manual BAD decisions"
    archived = list((reports / "history/squid_jumps").glob(f"*/{output.name}"))
    assert len(archived) == 1 and squid.file_stamp(archived[0]) == before_output
    # Resume compares newly saved, fractional-time endpoints on native samples.
    rebuilt_report = json.loads(report_path.read_text())
    entry = rebuilt_report["runs"][0]
    stage._validate_resumable_output(path=output, record=record, applied_bads=["MEG0121"],
        destination_run=1, destination_transform=raw.info["dev_head_t"],
        destination_fingerprint=rebuilt_report["destination"]["transform_sha256"],
        raw_root=resolved["resolved_paths"]["data_root"], squid_intervals=[[5401,6602]],
        expected_signature=mne.io.read_raw_fif(output, verbose="ERROR").info["description"].split("stage4_input_sha256=")[1])
