"""Verify Stage 3 measurement, donor restrictions and immutable handoffs.

Usage: python -m pytest meg/tests/test_visual_timing.py -q

Synthetic native-grid waveforms exercise signed edges, noise and cadence. Tiny
temporary Stage 1/2 fixtures exercise real fingerprints and the future epoching
consumer. No participant signals, shared derivatives or scheduler are touched.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

CODE = Path(__file__).resolve().parents[1] / "code/preprocessing_basic"
sys.path.insert(0, str(CODE))
import trial_preparation as tp
import visual_timing as vt
import visual_timing_utils as vu
from test_trial_preparation import make_case


def pulse(sfreq=1000, lag=19, sign=1):
    """Create a native 250-ms window with one four-frame physical pulse."""
    anchor = round(.075 * sfreq)
    y = np.zeros(round(.250 * sfreq))
    edge = anchor + vu.round_half_away(lag * sfreq / 1000)
    y[edge:edge + round(4 * sfreq / 120)] = .03 * sign
    return y, anchor


@pytest.mark.parametrize("sfreq,sign", [(500, 1), (1000, 1), (2000, -1)])
def test_native_grid_and_fixed_polarity(sfreq, sign):
    """Native offsets/polarity are preserved without resampling or 1-kHz assumptions."""
    y, anchor = pulse(sfreq, sign=sign)
    polarity, evidence = vu.infer_polarity([(y, anchor)] * 3, sfreq)
    assert polarity == sign and evidence["n_controls"] == 3
    row = vu.measure_edge(y, sfreq, anchor, polarity, 4, 120)
    assert row["measurement_status"] == "measured"
    assert row["edge_offset_samples"] == vu.round_half_away(.019 * sfreq)
    assert row["lag_ms"] == row["edge_offset_samples"] * 1000 / sfreq


@pytest.mark.parametrize("kind", ["flat", "one_sample_spike", "nan", "clipped_window"])
def test_no_false_observation_from_invalid_signal(kind):
    """A missing or truncated response cannot become a zero-lag measurement."""
    y, anchor = pulse()
    if kind == "flat": y[:] = .03
    if kind == "one_sample_spike": y[:] = 0; y[94] = .03
    if kind == "nan": y[20] = np.nan
    if kind == "clipped_window": y = y[:120]
    row = vu.measure_edge(y, 1000, anchor, 1, 4, 120)
    assert row["measurement_status"] == "unusable"
    assert row["lag_ms"] is None


def test_signed_negative_lag_is_retained_with_baseline_warning():
    """The diagnostic search exposes negative timing, not forced nonnegative lags."""
    y, anchor = pulse(lag=-10)
    row = vu.measure_edge(y, 1000, anchor, 1, 4, 120)
    assert row["edge_offset_samples"] == -10
    assert row["measurement_status"] == "flagged"
    assert "baseline_contaminated" in row["qc_reasons"]


def test_warmup_previous_high_requires_new_rise():
    """Cadence handling must not mistake a previous plateau for a new onset."""
    y, anchor = pulse()
    y[:anchor+10] = .03
    row = vu.measure_edge(y, 1000, anchor, 1, 4, 120, warmup=True)
    assert row["edge_offset_samples"] == 19
    assert row["baseline_method"] == "warmup_post_anchor_low_trough"


def test_condition_after_warmup_excludes_the_previous_physical_pulse():
    """Known cadence contamination must not assign the preceding warm-up rise."""
    y,anchor=pulse(lag=16)
    y[:anchor+8]=.03
    # Include a preceding warm-up low interval and rise in the bilateral search.
    y[anchor-30:anchor-16]=0
    row=vu.measure_edge(y,1000,anchor,1,4,120,post_warmup=True)
    assert row["measurement_status"]=="measured"
    assert row["edge_offset_samples"]==16
    assert row["baseline_method"]=="post_warmup_low_trough"


def test_low_trough_exception_does_not_hide_arbitrary_bad_baselines():
    """Without the accepted warm-up context, the same waveform stays flagged."""
    y,anchor=pulse(lag=16)
    y[:anchor+8]=.03
    row=vu.measure_edge(y,1000,anchor,1,4,120)
    assert row["measurement_status"]!="measured"


def test_second_response_is_flagged_and_first_edge_kept():
    """Multiple transitions remain inspectable without moving to the later pulse."""
    y, anchor = pulse()
    y[anchor+65:anchor+95] = .03
    row = vu.measure_edge(y, 1000, anchor, 1, 4, 120)
    assert row["edge_offset_samples"] == 19
    assert "multiple_candidate_rises" in row["qc_reasons"]


@pytest.mark.parametrize("value,expected", [(.5, 1), (-.5, -1), (1.5, 2), (-1.5, -2)])
def test_explicit_signed_rounding(value, expected):
    assert vu.round_half_away(value) == expected


def test_decided_donor_directions_and_code_separation():
    """Neither another participant nor a different code/early block can donate."""
    target = {"condition_code": 71, "behavioral_block": 1, "acquisition_run": 1}
    donor = {"subject_bids_id": "sub-11", "measurement_status": "measured", "trigger_code": 71,
             "behavioral_block": 3, "acquisition_run": 3}
    assert vu.donor_allowed("11", target, donor)
    assert not vu.donor_allowed("11", target, dict(donor, behavioral_block=2))
    assert not vu.donor_allowed("11", target, dict(donor, subject_bids_id="sub-12"))
    assert not vu.donor_allowed("11", target, dict(donor, trigger_code=81))
    assert not vu.donor_allowed("11", target, dict(donor, measurement_status="flagged"))
    donor.update(subject_bids_id="sub-24", acquisition_run=5)
    assert vu.donor_allowed("24", target, donor)
    assert not vu.donor_allowed("24", target, dict(donor, acquisition_run=4))
    assert not vu.donor_allowed("17", target, donor)


def test_estimate_balances_blocks_and_exposes_transfer_uncertainty():
    """A long donor block must not silently dominate the location estimate."""
    donors = [dict(source_fif=block, behavioral_run=1, lag_ms=lag, event_key=f"{block}:{i}")
              for block, n, lag in [("b3", 100, 18), ("b4", 10, 20), ("b5", 10, 21)]
              for i in range(n)]
    method = dict(vu.DEFAULT_METHOD, bootstrap_resamples=100)
    first = vu.estimate_lag(donors, method)
    assert first == vu.estimate_lag(donors, method)
    assert first["lag_estimate_ms"] == 20 and first["event_weighted_median_ms"] == 18
    assert first["transfer_uncertainty"].startswith("unknown")
    assert len(first["held_out_block_errors_ms"]) == 3


def test_coverage_uses_first_sample_exclusive_end_and_uncertainty():
    """A plausible shift can change window coverage without changing inclusion."""
    row = {"sfreq_hz": 1000, "raw_sample_frame": 10600, "recording_first_sample_raw": 10000,
           "recording_end_sample_raw_exclusive": 13300, "trial_interval_end_raw_exclusive": 13400}
    c = vu.window_coverage(row, 10620, -.5, 2.667, (10, 40))
    assert c["proposed_window_inside_recording"]
    assert not c["uncertainty_window_inside_recording"]
    assert vu.window_coverage(row, 10633, -.5, 2.667)["proposed_window_inside_recording"] is False


def handoff_case(tmp_path, *, pilot=False):
    """Build a genuine hashed Stage 2 bundle and synthetic reviewed overlay."""
    cfg, _ = make_case(tmp_path)
    tp.prepare_subject(cfg)
    prepared = tp.load_prepared_decisions(cfg)
    rows = []
    for r in prepared:
        rows.append(dict(r, assignment_basis="digital_fallback", assignment_status="proposed",
                         proposed_sample_raw=r["raw_sample_frame"], applied_shift_samples=0,
                         lag_estimate_ms=None, uncertainty_kind="physical_lag_unknown"))
    folder = Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing/test/sub-30/ses-01"
    report = {"schema_version":"1.0", "status":"completed", "mode":"pilot" if pilot else "cohort_proposal",
              "subject_bids_id":"sub-30", "version":"test", "stage2_report_sha256":tp.file_hash(tp.stage2_paths(cfg)["report"]),
              "epoch_window_sec":[-.5,2.667]}
    contents = {"trial_timing_map.csv":tp.csv_text(rows), "event_measurements.csv":"event_key,measurement_status\n",
                "event_timing_proposals.csv":"event_key,assignment_basis\n", "recording_summary.csv":"source_fif\n",
                "estimates.json":"{}", "contexts.json":"{}", "report.md":"Synthetic", "diagnostics.png":b"synthetic"}
    vt.publish(folder, contents, report, cfg["resolved_paths"]["data_root"])
    review = folder.parent/"review.json"
    review.write_text(json.dumps({"schema_version":"1.0","status":"accepted","reviewer":"synthetic test",
        "reviewed_utc":"2026-10-20T00:00:00Z","decision_reference":"fixture",
        "bundle_manifest":"ses-01/manifest.json", "bundle_manifest_sha256":tp.file_hash(folder/"manifest.json"),
        "accepted_assignment_bases":["digital_fallback"]}))
    return cfg, prepared, folder, review


def test_handoff_preserves_original_keys_exclusions_and_unknown_lag(tmp_path):
    cfg, prepared, folder, review = handoff_case(tmp_path)
    before = tp.evidence_fingerprints(tp.evidence_paths(cfg))
    output, gate = vu.load_reviewed_timing(cfg, prepared, review)
    assert [r["trial_key"] for r in output] == [r["trial_key"] for r in prepared]
    assert [r["included_by_trial_policy"] for r in output] == [False, False, True]
    assert output[2]["epoch_event_sample_raw"] == prepared[2]["raw_sample_frame"]
    assert output[2]["timing_lag_ms"] == ""
    assert gate["version"] == "test"
    assert tp.evidence_fingerprints(tp.evidence_paths(cfg)) == before
    with pytest.raises(FileExistsError):
        vt.publish(folder, {}, {}, cfg["resolved_paths"]["data_root"])


def test_pilot_never_satisfies_scientific_handoff(tmp_path):
    cfg, prepared, _, review = handoff_case(tmp_path, pilot=True)
    with pytest.raises(RuntimeError, match="Pilot"):
        vu.load_reviewed_timing(cfg, prepared, review)


def test_edited_map_and_changed_ledger_stop_handoff(tmp_path):
    cfg, prepared, folder, review = handoff_case(tmp_path)
    path = folder/"trial_timing_map.csv"
    path.write_text(path.read_text()+"\n")
    with pytest.raises(RuntimeError, match="fingerprint"):
        vu.load_reviewed_timing(cfg, prepared, review)


def test_stale_stage2_completion_stops_reviewed_handoff(tmp_path):
    cfg, prepared, _, review = handoff_case(tmp_path)
    path = tp.stage2_paths(cfg)["report"]
    path.write_text(path.read_text()+"\n")
    with pytest.raises(RuntimeError, match="Stage 2 completion"):
        vu.load_reviewed_timing(cfg, prepared, review)


def test_full_standalone_proposal_preserves_accepted_products(tmp_path):
    """Exercise publication/resume against real tiny FIF and saved-MAT inputs."""
    from scipy.io import savemat
    cfg, _ = make_case(tmp_path)
    tp.prepare_subject(cfg)
    profile=cfg["subject_profile"]
    profile.update(versions={"runtime":"v25_dotOffsetBlank","stimuli":"V30","schedule":"v23"},
                   completed_behavioral_blocks=[2])
    folder=Path(cfg["resolved_paths"]["subject_data_dir"])/"output_files";folder.mkdir()
    savemat(folder/"MoveDot1_occlusion_v25_dotOffsetBlank_SUB30_BLOCK02.mat",
            {"Conf":{"experimentName":"MoveDot1_occlusion_v25_dotOffsetBlank","refrate":120}})
    authority=Path(cfg["project"]["docs_root"])/"PREPROCESSING_DECISIONS.md"
    authority.write_text(authority.read_text()+"\nStage 2 trial-accounting gate outcome\n")
    before=tp.evidence_fingerprints(tp.evidence_paths(cfg))
    report=vt.run_subject(cfg,"test",vu.DEFAULT_METHOD)
    assert report["n_original_trial_keys"]==3 and report["n_signal_windows_read"]==0
    bundle=Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing/test/sub-30/ses-01"
    assert vu.verify_bundle(bundle)["input_contract_sha256"]==report["input_contract_sha256"]
    assert len(tp.read_csv(bundle/"trial_timing_map.csv"))==3
    assert vt.run_subject(cfg,"test",vu.DEFAULT_METHOD,resume=True)["status"]=="completed"
    with pytest.raises(RuntimeError,match="immutable"):
        vt.run_subject(cfg,"test",dict(vu.DEFAULT_METHOD,minimum_snr=9),resume=True)
    assert tp.evidence_fingerprints(tp.evidence_paths(cfg))==before


def test_actual_epoch_consumer_uses_overlay_without_changing_trial_identity(tmp_path,monkeypatch):
    """Inspect the actual MNE boundary: shifted event, original metadata/key."""
    cfg,prepared,folder,review=handoff_case(tmp_path)
    rows=tp.read_csv(folder/"trial_timing_map.csv")
    row=rows[-1]
    row.update(assignment_basis="measured",proposed_sample_raw="67519",applied_shift_samples="19",
               lag_estimate_ms="19.0",measurement_event_key="edge")
    measurements=[{"event_key":"edge","measurement_status":"measured","source_fif":row["source_fif"],
                   "digital_sample_raw":67500,"trigger_code":71,"diode_edge_sample_raw":67519}]
    for name,content in {"trial_timing_map.csv":tp.csv_text(rows),"event_measurements.csv":tp.csv_text(measurements)}.items():
        (folder/name).write_text(content)
    manifest=json.loads((folder/"manifest.json").read_text())
    for name in ("trial_timing_map.csv","event_measurements.csv"):
        manifest["output_sha256"][name]=tp.file_hash(folder/name)
    (folder/"manifest.json").write_text(json.dumps(manifest))
    accepted=json.loads(review.read_text());accepted.update(bundle_manifest_sha256=tp.file_hash(folder/"manifest.json"),
                                                          accepted_assignment_bases=["measured"])
    review.write_text(json.dumps(accepted))
    spec=importlib.util.spec_from_file_location("stage3_epoch_consumer_test",CODE/"04_epoching.py")
    epoching=importlib.util.module_from_spec(spec);spec.loader.exec_module(epoching)
    source=row["source_fif"]
    raw_path=Path(cfg["resolved_paths"]["subject_data_dir"])/"meg"/source
    runs=[SimpleNamespace(cleaned_path=str(raw_path),source=SimpleNamespace(source_filename=source))]
    monkeypatch.setattr(epoching,"resolve_manual_ica_apply_runs",lambda *args:(runs,{"status":"passed"}))
    class ReachedBoundary(Exception):pass
    def inspect(raw,events,**kwargs):
        assert events[:,0].tolist()==[67519]
        assert kwargs["metadata"]["raw_sample_frame"].tolist()==[67500]
        assert kwargs["metadata"]["trial_key"].tolist()==[prepared[-1]["trial_key"]]
        raise ReachedBoundary
    monkeypatch.setattr(epoching.mne,"Epochs",inspect)
    cfg["epoching"]["require_stage3_timing"]=True
    with pytest.raises(RuntimeError,match="requires --timing-review"):
        epoching.run_epoching(cfg,source_manifest=tp.evidence_paths(cfg)["manifest"],overwrite=False)
    with pytest.raises(ReachedBoundary):
        epoching.run_epoching(cfg,source_manifest=tp.evidence_paths(cfg)["manifest"],overwrite=False,timing_review=review)


def test_pilot_selection_bound_and_contexts():
    """The fixed pilot selects before reading signals and obeys donor coverage."""
    events = [dict(acquisition_run=str(run), trigger_code=str(code), raw_sample_frame=str(run*10000+i),
                   event_index=str(run*100+i)) for run in range(1,6) for i,code in enumerate([71,81,91,111,118,153]*5)]
    for subject in ("08","09","11","24"):
        selected = vt.select_pilot(events, subject)
        assert len(selected) <= 12 and len({r["event_index"] for r in selected}) == len(selected)
    assert {r["acquisition_run"] for r in vt.select_pilot(events,"24")} == {"5"}
    assert {r["acquisition_run"] for r in vt.select_pilot(events,"11")} == {"1","2","3","4","5"}
    with pytest.raises(ValueError): vt.select_pilot(events,"04")


def test_cached_pilot_reassessment_never_reads_raw_signals(monkeypatch):
    """A revised candidate method reuses the exact native cached windows."""
    events=[];cache={}
    source=SimpleNamespace(source_filename="sub08_run01_test.fif",source_path="unused")
    for i,code in enumerate([71,81,91]):
        row=dict(source_fif=source.source_filename,acquisition_run="1",trigger_code=str(code),
                 raw_sample_frame=str(11000+1000*i),event_index=str(i),behavioral_block="1")
        events.append(row)
        y,anchor=pulse()
        sti=np.zeros(len(y));sti[anchor:anchor+2]=code
        cache[vt.event_key("08",row)]=(np.vstack([sti,y]),anchor,1000.)
    class HeaderOnly:
        ch_names=["STI101","MISC008"]
        first_samp=10000;n_times=6000;info={"sfreq":1000.}
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def get_data(self,**kwargs):raise AssertionError("Cached pilot reread raw signals")
    monkeypatch.setattr(vt.mne.io,"read_raw_fif",lambda *args,**kwargs:HeaderOnly())
    measured,waves,reads=vt.measure_events({"subject":{"id":"08"}},events,[source],
        {"1":{"runtime":"v24_allTrialGuard","refrate":120}},vu.DEFAULT_METHOD,pilot=True,cached_windows=cache)
    assert reads==0 and len(waves)==3
    assert all(r["measurement_status"]=="measured" for r in measured)


@pytest.mark.parametrize("value", ["../raw", "a/other", "x..y", ""])
def test_version_traversal_rejected(tmp_path,value):
    cfg={"resolved_paths":{"derivatives_root":str(tmp_path/"derivatives"),"data_root":str(tmp_path/"raw")}}
    with pytest.raises(ValueError):vt.version_root(cfg,value)


@pytest.mark.parametrize("mode", ["subject", "summary"])
def test_stage3_wrapper_executes_only_standalone_timing(tmp_path,mode):
    """Run actual wrappers against harmless environment/Python stubs."""
    repo = Path(__file__).resolve().parents[2]
    cluster = tmp_path/"project/meg/cluster";cluster.mkdir(parents=True)
    wrapper = repo/"meg/cluster/submit_preprocessing_stage3.sh"
    shutil.copy2(wrapper,cluster/wrapper.name)
    binary=tmp_path/"bin";binary.mkdir()
    mamba=binary/"micromamba";mamba.write_text("#!/bin/bash\nprintf 'micromamba() { :; }\\n'\n");mamba.chmod(0o755)
    python=binary/"python";python.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$DAD_TEST_OUTPUT"\n');python.chmod(0o755)
    captured=tmp_path/"calls"
    env=dict(os.environ,PATH=str(binary)+os.pathsep+os.environ["PATH"],DAD_PROJECT_ROOT=str(tmp_path/"project"),
             MICROMAMBA_BIN=str(mamba),TMPDIR=str(tmp_path/"cache"),SLURM_ARRAY_TASK_ID="24",
             DAD_STAGE3_VERSION="test",DAD_TEST_OUTPUT=str(captured))
    for key in ("DAD_STAGE3_MODE","DAD_STAGE3_PILOT","DAD_STAGE3_RESUME","DAD_STAGE3_VALIDATE_ONLY"):env.pop(key,None)
    target=wrapper if mode=="subject" else wrapper.with_name("submit_preprocessing_stage3_summary.sh")
    result=subprocess.run(["bash",str(target)],env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    calls=captured.read_text().splitlines()
    assert len(calls)==1 and "visual_timing.py" in calls[0] and "run_pipeline.py" not in calls[0]
    assert ("--summarize-cohort" in calls[0])==(mode=="summary")
