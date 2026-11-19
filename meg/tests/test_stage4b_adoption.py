"""Test policy adoption on temporary native FIFs without cohort data or jobs.

Usage
-----
    python -B -m pytest -p no:cacheprovider meg/tests/test_stage4b_adoption.py -q

Exercise real Stage 4B save/resume/reconciliation and prerequisite acceptance,
including nonzero first_samp, missing measurement date, source metadata, manual
reviewed-empty muscle, immutable history and deliberate stale-input failures.
Only accepted timing is stubbed; MNE FIF/CSV/preservation and policy gates run.
"""

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import os
import shlex
import subprocess
from types import SimpleNamespace

import mne
import numpy as np
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import stage4b_adoption as adoption
import annotation_triage as triage
from continuous_cleaning import (annotation_rows, coordinates, digest, file_hash,
    file_stamp, muscle_parameters, require_review, stage_configuration)
from manual_eye_annotations import manual_path
from preprocessing_stage_utils import AcquisitionRunRecord, stage_output_dirs


# -----------------------------------------------------------------------------
# Synthetic sources, exact policy plan, and real artifact-file evidence
# -----------------------------------------------------------------------------

def fixture(tmp_path, monkeypatch, *, manual_empty, date):
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] / "params/pipeline_default.yaml").read_text())
    base = tmp_path / "derivatives"
    cfg.update(subject={"id": "04"}, layout={"subject_bids_id": "sub-04", "session_bids_id": "ses-01"})
    cfg["subject_profile"] = cfg["subject_profiles"]["04"]
    cfg["resolved_paths"] = {"derivatives_root": str(base), "data_root": str(tmp_path / "raw"),
        "derivatives_preprocessing_dir": str(base / "preprocessing/sub-04/ses-01"),
        "derivatives_reports_dir": str(base / "reports/sub-04/ses-01")}
    for key in ("calibration_file", "crosstalk_file"):
        path = tmp_path / key
        path.write_text("test calibration")
        cfg["maxwell_filter"][key] = str(path)
    cfg["preprocessing"].update(highpass_hz=None, filter_edge_margin_sec=None)
    cfg["squid_jumps"]["enabled"] = False
    cfg["epoching"].update(tmin_sec=-.1, tmax_sec=.2, reject_by_annotation=True)
    dirs = stage_output_dirs(cfg)
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(97)
    data = rng.normal(size=(4, 6000)) * np.array([1e-12, 1e-12, 1., 1.])[:, None]
    data[2] = np.arange(6000) % 255
    raw = mne.io.RawArray(data, mne.create_info(["MEG0111", "MEG0121", "STI101", "MISC001"],
        500., ["mag", "mag", "stim", "misc"]), first_samp=1397000, verbose="ERROR")
    if date:
        raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    raw.set_annotations(mne.Annotations([1., 2., 0.], [.25, .5, 0.],
        ["BAD_source", "BAD_muscle", "EDGE boundary"], ch_names=[[], [], []],
        extras=[{"source": "immutable"}, {}, {}]))
    record = AcquisitionRunRecord(str(tmp_path / "raw/sub04_run01_test.fif"),
        "sub04_run01_test.fif", "sub04_run01_test", 1, "sub04_run01_test.fif")
    source = dirs["sss"] / (record.source_stem + "_raw_sss.fif")
    raw.save(source, verbose="ERROR")
    header = mne.io.read_raw_fif(source, verbose="ERROR")
    basis = dirs["sss"] / (record.source_stem + "_ann.fif")
    header.save(basis, verbose="ERROR")
    oldstamp = file_stamp(basis)
    header.annotations.save(basis.with_suffix(".csv"))
    tuning = dirs["reports"] / "03_tuning_choices.yaml"
    a_report = {"configuration_sha256": stage_configuration(cfg, "03_A_MaxFilt"),
        "evidence_sha256": "test-accepted-sss", "evidence_files": []}
    (dirs["reports"] / "03_A_MaxFilt_review_report.json").write_text(json.dumps(a_report))
    cfg["stage4_review"] = {"03_A_MaxFilt": {"status": "accepted", **{k: a_report[k] for k in ("configuration_sha256", "evidence_sha256")}}}
    receipt = None
    maskraw = header.copy()
    if manual_empty:
        maskraw.annotations.delete([i for i, a in enumerate(maskraw.annotations) if a["description"] == "BAD_muscle"])
        maskraw.annotations.append(maskraw.first_time + 4., .2, "BAD_blink")
        path = manual_path(cfg, record)
        path.parent.mkdir(parents=True)
        ann = maskraw.annotations[[i for i, a in enumerate(maskraw.annotations) if a["description"].startswith("BAD")]]
        ann.save(path)
        receipt = {"status": "reviewed", "coordinates": coordinates(header), "sss_input": file_stamp(source),
            "annotation_sha256": file_hash(path), "annotation_count": len(ann),
            "muscle_review": {"status": "reviewed", "annotation_count": 0,
                "settings_sha256": digest(muscle_parameters(cfg["artifact_annotation"], 10.))},
            "blink_auto": {"removed_spans_raw": [[1397500, 1397600]]}}
        cfg["artifact_annotation"]["manual_eye_reviews"] = {record.source_filename: receipt}
    else:
        maskraw.annotations.append(maskraw.first_time + 4., .2, "BAD_blink_auto")
    mask = {"schema": triage.SCHEMA, "purpose": triage.PURPOSE,
        "method": "preserved_manual" if manual_empty else "provisional_automatic",
        "coverage_status": "manual_review_preserved" if manual_empty else "selected_components_only",
        "input_binding": {"coordinates": coordinates(header), "muscle_basis": {"threshold": 10.}},
        "detectors": {}, "annotations": triage.mask_rows(maskraw)}
    mask["payload_sha256"] = digest(mask)
    rows = [{"trial_key": f"test-{i}", "source_fif": record.source_filename, "sfreq_hz": 500.,
        "recording_first_sample_raw": header.first_samp, "recording_end_sample_raw_exclusive": header.first_samp + header.n_times,
        "condition_code": 71, "trial_sequence_identity": i + 1, "included_by_trial_policy": True,
        "final_analysis_eligible": True, "epoch_event_sample_raw": header.first_samp + i * 500}
        for i in range(1, 8)]
    gate = {"gate": {"application_sha256": "test-accepted-timing"}, "prepared_rows_sha256": digest(rows)}
    snapshot = base / "annotation_triage/reports/test/report.json"
    snapshot.parent.mkdir(parents=True)
    snapshot.write_text('{"status":"completed"}')
    protected = [adoption.protected(cfg, source, small=False), adoption.protected(cfg, basis, small=False, basis=True)]
    if manual_empty:
        protected.append(adoption.protected(cfg, manual_path(cfg, record)))
    participant = {"subject": 4, "base_configuration_sha256": adoption.base_configuration(cfg), "timing": gate,
        "protected_files": protected, "runs": [{"source_filename": record.source_filename,
            "mask": mask, "manual_eye_review": receipt or {"status": "unreviewed"}, "manual_superseded_snapshot": False}],
        "expected_trials": adoption.trial_footprint(triage.count_trials(rows, header, mask, cfg["epoching"]))}
    plan = {"policy": adoption.POLICY, "cohort": list(adoption.COHORT), "participants": [participant],
        "snapshot": adoption.protected(cfg, snapshot)}
    plan["plan_sha256"] = digest(plan)
    target = adoption.root(cfg) / "test"
    target.mkdir(parents=True)
    (target / "plan.json").write_text(json.dumps(plan, indent=2))
    tuning.write_text(yaml.safe_dump({"stage4_review": cfg["stage4_review"], "artifact_annotation": {
        "manual_eye_reviews": cfg["artifact_annotation"].get("manual_eye_reviews", {})}}))

    def current(*args):
        result = deepcopy(cfg)
        saved = yaml.safe_load(tuning.read_text())
        result["artifact_annotation"].update(saved.get("artifact_annotation", {}))
        result["stage4_review"].update(saved.get("stage4_review", {}))
        return result

    monkeypatch.setattr(triage, "subject_config", current)
    monkeypatch.setattr(adoption, "timing", lambda cfg: (deepcopy(rows), deepcopy(gate)))
    monkeypatch.setattr(adoption, "discover_raw_run_records", lambda cfg: [record])
    # Presentation is tested separately; this test isolates real scientific
    # output saving and guards without needing a synthetic muscle-score PNG.
    monkeypatch.setattr(adoption, "render_stage_report", lambda **kwargs: None)
    return argparse.Namespace(subject=4, label="test"), cfg, current, plan, participant, record, oldstamp, tuning


@pytest.mark.parametrize("manual_empty", [False, True])
@pytest.mark.parametrize("date", [False, True])
def test_real_adoption_resume_reconcile_and_processing_acceptance(tmp_path, monkeypatch, manual_empty, date):
    args, cfg, current, plan, participant, record, oldstamp, tuning = fixture(tmp_path, monkeypatch, manual_empty=manual_empty, date=date)
    original_manual = deepcopy(cfg["artifact_annotation"].get("manual_eye_reviews", {}))
    adoption.apply_participant(args)
    basis = stage_output_dirs(cfg)["sss"] / (record.source_stem + "_ann.fif")
    archived = adoption.history_path(cfg, plan, basis)
    assert file_stamp(archived) == oldstamp
    assert basis.exists()
    saved = mne.io.read_raw_fif(basis, verbose="ERROR")
    assert any(a.get("extras", {}).get("source") == "immutable" for a in saved.annotations)
    assert ("BAD_muscle" in saved.annotations.description) is not manual_empty
    assert current()["artifact_annotation"].get("manual_eye_reviews", {}) == original_manual
    # Processing alone never writes an acceptance receipt.
    with pytest.raises(RuntimeError, match="unreviewed or changed"):
        require_review(current(), adoption.STAGE, stage_output_dirs(cfg)["reports"])
    adoption.apply_participant(args)  # Exact-plan interruption retry preserves history.
    result = adoption.reconcile_participant((args, 4))
    adoption.accept_participant((args, result))
    with pytest.raises(RuntimeError, match="Whole-cohort"):
        require_review(current(), adoption.STAGE, stage_output_dirs(cfg)["reports"])
    # Simulate the parent's joined publication using an explicit synthetic roster.
    path = adoption.root(cfg) / "test/plan.json"
    joined = [result if subject == 4 else {"subject": subject, "evidence_sha256": "synthetic-other-participant"}
              for subject in adoption.COHORT]
    (path.parent / "policy_acceptance.json").write_text(json.dumps({"policy": adoption.POLICY,
        "adoption_plan_sha256": plan["plan_sha256"], "status": "accepted_for_final_ica_fit", "participants": joined}))
    receipt = require_review(current(), adoption.STAGE, stage_output_dirs(cfg)["reports"])
    assert receipt["acceptance_kind"] == "explicit_annotation_policy_for_processing"
    report = json.loads((stage_output_dirs(cfg)["reports"] / f"{adoption.STAGE}_review_report.json").read_text())
    assert report["runs"][0]["manual_eye_review"]["status"] == ("reviewed" if manual_empty else "unreviewed")
    # A new SAVE/source/model change cannot inherit the old policy acceptance.
    source = stage_output_dirs(cfg)["sss"] / (record.source_stem + "_raw_sss.fif")
    source.touch()
    with pytest.raises(RuntimeError, match="changed"):
        require_review(current(), adoption.STAGE, stage_output_dirs(cfg)["reports"])


@pytest.mark.parametrize("changed", ["manual_file", "manual_receipt", "plan", "basis"])
def test_changed_inputs_stop_before_production_replacement(tmp_path, monkeypatch, changed):
    args, cfg, current, plan, participant, record, oldstamp, tuning = fixture(tmp_path, monkeypatch, manual_empty=True, date=False)
    if changed == "manual_file":
        manual_path(cfg, record).write_bytes(b"changed saved review")
    elif changed == "manual_receipt":
        data = yaml.safe_load(tuning.read_text())
        data["artifact_annotation"]["manual_eye_reviews"][record.source_filename]["annotation_count"] += 1
        tuning.write_text(yaml.safe_dump(data))
    elif changed == "plan":
        path = adoption.root(cfg) / "test/plan.json"
        data = json.loads(path.read_text())
        data["participants"][0]["runs"][0]["mask"]["annotations"] = []
        path.write_text(json.dumps(data))
    else:
        (stage_output_dirs(cfg)["sss"] / (record.source_stem + "_ann.fif")).touch()
    with pytest.raises(RuntimeError):
        adoption.apply_participant(args)
    basis = stage_output_dirs(cfg)["sss"] / (record.source_stem + "_ann.fif")
    assert not adoption.history_path(cfg, plan, basis).exists()


def test_missing_property_diagnostics_cannot_reuse_old_pngs(tmp_path):
    from final_ica_review_book import write_component_book
    from matplotlib import pyplot as plt
    figure, ax = plt.subplots()
    ax.plot([0, 1], [1, 0])
    path = tmp_path / "03_C_ICA_review_ICA000_properties_01.png"
    figure.savefig(path)
    plt.close(figure)
    payload = {"subject_bids_id": "sub-04", "n_components_fitted": 2,
        "prefit_model_sha256": "exact-final-model", "component_properties_requested": True}
    coverage = write_component_book(payload, [path], tmp_path / "book.html", tmp_path / "raw")
    assert not coverage["complete"] and coverage["missing_properties"] == [1]
    assert "exact-final-model" in (tmp_path / "book.html").read_text()


def test_plan_discovers_new_manual_save_without_reseeding_automatic_intervals(tmp_path, monkeypatch):
    args, cfg, current, plan, participant, record, oldstamp, tuning = fixture(tmp_path, monkeypatch, manual_empty=True, date=False)
    entry = participant["runs"][0]
    source = stage_output_dirs(cfg)["sss"] / (record.source_stem + "_raw_sss.fif")
    raw = adoption.compose_annotations(mne.io.read_raw_fif(source, verbose="ERROR"), entry["mask"])
    binding = {**entry["mask"]["input_binding"], "manual_receipt_sha256": digest(entry["manual_eye_review"])}
    monkeypatch.setattr(triage, "load_context", lambda *args: (raw.copy(), source, binding))
    monkeypatch.setattr(adoption, "load_cached_aid", lambda *args: pytest.fail("Manual review must not load preliminary selections/models"))
    proposal = Path(cfg["resolved_paths"]["derivatives_root"]) / "annotation_triage/proposals/sub-04/old.json"
    proposal.parent.mkdir(parents=True)
    proposal.write_text(json.dumps(entry["mask"]))
    result = adoption.validate_participant((args, 4, {"timing_gate": participant["timing"]["gate"], "masks": [{"source_fif": record.source_filename,
        "path": "/mnt/shared/annotation_triage/proposals/sub-04/old.json"}]}))
    assert result["runs"][0]["mask"]["method"] == "preserved_manual"
    assert result["runs"][0]["manual_superseded_snapshot"]
    assert not any(row["description"] == "BAD_muscle" for row in result["runs"][0]["mask"]["annotations"])
    assert result["runs"][0]["manual_eye_review"]["blink_auto"]["removed_spans_raw"] == [[1397500, 1397600]]


@pytest.mark.parametrize("mode,task,accept,success", [
    ("apply", "4", "1", True), ("apply", "15", "1", False),
    ("apply", None, "1", False), ("plan", None, "1", True),
    ("reconcile", None, "1", True), ("reconcile", None, "0", True),
    ("reconcile", "4", "1", False)])
def test_cluster_wrapper_dispatch_and_cohort_phase_barriers(tmp_path, mode, task, accept, success):
    project = tmp_path / "project"
    project.mkdir()
    binary = tmp_path / "bin"
    binary.mkdir()
    mamba = binary / "micromamba"
    mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n")
    mamba.chmod(0o755)
    runtime = binary / "python"
    runtime.write_text('#!/bin/bash\nprintf "%s\\n" "$*" > "$DAD_TEST_OUTPUT"\n')
    runtime.chmod(0o755)
    calls = tmp_path / "calls.txt"
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
        MICROMAMBA_BIN=str(mamba), DAD_PROJECT_ROOT=str(project), DAD_ADOPTION_LABEL="test",
        DAD_ADOPTION_MODE=mode, DAD_ADOPTION_ACCEPT=accept, SLURM_CPUS_PER_TASK="8",
        DAD_TEST_OUTPUT=str(calls), TMPDIR=str(tmp_path / "tmp"))
    env.pop("SLURM_ARRAY_TASK_ID", None)
    if task is not None:
        env["SLURM_ARRAY_TASK_ID"] = task
    wrapper = Path(__file__).resolve().parents[1] / "cluster/submit_stage4b_adoption.sh"
    result = subprocess.run(["bash", str(wrapper)], env=env, capture_output=True, text=True)
    assert (result.returncode == 0) is success, result.stdout + result.stderr
    if not success:
        assert not calls.exists()
        return
    tokens = shlex.split(calls.read_text())
    assert tokens[2] == mode
    assert ("--subject" in tokens) is (mode == "apply")
    assert ("--accept-policy" in tokens) is (mode == "reconcile" and accept == "1")
    if mode != "apply":
        assert tokens[tokens.index("--jobs") + 1] == "8"


def test_p10_slow_qc_keeps_native_signals_and_three_long_windows(tmp_path):
    import importlib.util
    spec = importlib.util.spec_from_file_location("p10_native_qc", Path(adoption.__file__).with_name("03_C_ICA.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sfreq = 200.
    times = np.arange(18000) / sfreq
    values = np.vstack([np.sin(2 * np.pi * times / 3) * 1e-12, np.cos(2 * np.pi * times / 3) * 1e-12])
    before = mne.io.RawArray(values, mne.create_info(["MEG0111", "MEG0121"], sfreq, "mag"), first_samp=10000, verbose="ERROR")
    after = before.copy()
    after._data *= .2
    record = SimpleNamespace(source_filename="sub10_run01_test.fif", source_stem="sub10_run01_test")
    evidence = module.save_p10_slow_qc(before, after, record, tmp_path)
    np.testing.assert_array_equal(before.get_data(), values)
    assert len(evidence["windows"]) == 3 and evidence["decision"] == "open"
    assert all(row["stop_sample_raw_exclusive"] - row["start_sample_raw"] == 12000 for row in evidence["windows"])
    assert Path(evidence["figure"]).is_file()


@pytest.mark.parametrize("date", [False, True])
def test_fif_boundary_rounding_is_predicted_and_reported_without_a_new_margin(tmp_path, monkeypatch, date):
    args, cfg, current, plan, participant, record, oldstamp, tuning = fixture(tmp_path, monkeypatch, manual_empty=False, date=date)
    entry = participant["runs"][0]
    # This interval starts exactly on trial 6's exclusive stop before FIF float32
    # serialization. Its shifted storage endpoint must remain visible in QC.
    first = entry["mask"]["input_binding"]["coordinates"]["first_samp"] / 500.
    entry["mask"]["annotations"].append({"onset": first + 6.202, "duration": .457, "description": "BAD_boundary"})
    entry["mask"]["annotations"].sort(key=lambda row: (row["onset"], row["duration"], row["description"]))
    entry["mask"]["payload_sha256"] = digest({k: v for k, v in entry["mask"].items() if k != "payload_sha256"})
    prepared, _ = adoption.timing(cfg)
    source = stage_output_dirs(cfg)["sss"] / (record.source_stem + "_raw_sss.fif")
    raw = mne.io.read_raw_fif(source, verbose="ERROR")
    participant["expected_trials"] = adoption.trial_footprint(triage.count_trials(prepared, raw, entry["mask"], cfg["epoching"]))
    plan["plan_sha256"] = digest({k: v for k, v in plan.items() if k != "plan_sha256"})
    path = adoption.root(cfg) / "test/plan.json"
    path.write_text(json.dumps(plan))
    adoption.apply_participant(args)
    result = adoption.reconcile_participant((args, 4))
    assert result["serialization_outcome_changes"]
    assert entry["mask"]["annotations"][-1]["duration"] == .457


@pytest.mark.parametrize("failure", [None, "audit", "table"])
def test_cohort_gate_is_published_only_after_all_audits_receipts_and_tables(tmp_path, monkeypatch, failure):
    fixture(tmp_path, monkeypatch, manual_empty=False, date=False)
    phases, writes = [], []
    results = [{"subject": subject, "evidence_sha256": "synthetic-audit", "configuration_sha256": "synthetic-settings",
        "condition_retention": [], "identity_retention": [], "run_condition_retention": []}
        for subject in adoption.COHORT]

    def phase(worker, tasks, workers):
        assert len(tasks) == len(adoption.COHORT)
        if worker is adoption.reconcile_participant:
            phases.append("all-participant-audits")
            if failure == "audit":
                raise RuntimeError("synthetic failed participant audit")
        else:
            assert phases == ["all-participant-audits"]
            phases.append("all-participant-receipts")
        return results

    def table(path, content, **kwargs):
        assert phases == ["all-participant-audits", "all-participant-receipts"]
        if failure == "table":
            raise RuntimeError("synthetic interrupted cohort table")
        writes.append("table")

    monkeypatch.setattr(adoption, "map_participants", phase)
    monkeypatch.setattr(adoption, "write_text", table)
    monkeypatch.setattr(adoption, "write_json", lambda *args, **kwargs: writes.append("cohort-gate"))
    command = ["reconcile", "--label", "test", "--accept-policy", "--jobs", "2"]
    if failure:
        with pytest.raises(RuntimeError, match="synthetic"):
            adoption.main(command)
        assert "cohort-gate" not in writes
        if failure == "audit":
            assert phases == ["all-participant-audits"]
    else:
        assert adoption.main(command) == 0
        assert writes == ["table", "table", "table", "cohort-gate"]


def test_equal_time_source_labels_keep_their_order_without_changing_mask_identity(tmp_path):
    """Equivalent interval multisets preserve source ordering and every label."""
    raw = mne.io.RawArray(np.zeros((1, 2000)), mne.create_info(["STI101"], 500., "stim"),
        first_samp=1397000, verbose="ERROR")
    for label in ("BAD_a", "BAD_z"):
        raw.annotations.append(raw.first_time + .333, .123, label)
    mask = {"method": "provisional_automatic", "coverage_status": "selected_components_only",
        "input_binding": {"coordinates": coordinates(raw)}, "detectors": {},
        "annotations": list(reversed(triage.mask_rows(raw)))}
    path = tmp_path / "source_raw.fif"
    raw.save(path, verbose="ERROR")
    saved, footprint, original, changes = adoption.verify_saved({"epoching": {}},
        SimpleNamespace(source_filename="one.fif"), {"mask": mask}, path, [])
    assert list(saved.annotations.description) == ["BAD_a", "BAD_z"]
    assert footprint == original == changes == []
