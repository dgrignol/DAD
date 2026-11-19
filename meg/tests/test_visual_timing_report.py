"""Check Stage 3 report statistics, categorical plots and immutable presentation.

Usage: PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python3 -B -m pytest \
    meg/tests/test_visual_timing_report.py -q

Synthetic rows exercise exact sample/lag arithmetic, missing dispersion, saved
context and donor evidence. Existing-bundle rendering is tested with guarded
temporary derivatives and a blocked raw-reader/detector. No participant signals,
shared outputs or scheduler are touched. Plot checks inspect actual category
positions rather than comparing implementation-generated image snapshots.
"""

from __future__ import annotations

import copy
import base64
import io
import json
import re
from pathlib import Path
import sys

import numpy as np
import pytest

CODE = Path(__file__).resolve().parents[1]/"code/preprocessing_basic"
sys.path.insert(0,str(CODE))
import visual_timing as vt
import visual_timing_report as vr
import visual_timing_utils as vu
from trial_preparation import csv_text,file_hash,read_csv


def observations():
    """Make different trigger factors, rates, QC states and recording strata."""
    rows = []
    for i,(source,code,lag,status) in enumerate([
        ("one.fif",71,18,"measured"),("one.fif",71,20,"measured"),
        ("one.fif",118,8,"measured"),("one.fif",71,-5,"flagged"),
        ("two.fif",71,30,"measured")]):
        sf = 500 if source=="two.fif" else 1000
        rows.append(dict(event_key=f"sub-11:{source}:{1000+i*100}:B{code}",subject_bids_id="sub-11",
            source_fif=source,runtime="v25",trigger_code=code,digital_sample_raw=1000+i*100,
            first_samp=500,sfreq_hz=sf,lag_ms=lag,measurement_status=status,
            qc_reasons="baseline_contaminated" if status=="flagged" else "",behavioral_block=3,
            acquisition_run=3,behavioral_run=1,baseline_method="pre_edge_median",
            post_warmup_context=i==0,saved_trial_onset_guard_retry_count=1 if i==0 else 0,
            amplitude_v=.03,baseline_noise_v=0,threshold_span_ms=2,polarity=1,
            baseline_normalized_v=0,edge_10_offset_samples=lag*sf/1000-1,
            edge_25_offset_samples=lag*sf/1000,edge_50_offset_samples=lag*sf/1000,
            edge_75_offset_samples=lag*sf/1000,edge_90_offset_samples=lag*sf/1000+1))
    return rows


def trial_map():
    """Include an original boundary, an assigned crossing and unresolved timing."""
    return [dict(trial_key=f"trial-{i}",source_fif="one.fif",subject_bids_id="sub-11",
                 raw_sample_frame=sample,condition_code=71,sfreq_hz=1000,
                 recording_first_sample_raw=500,recording_end_sample_raw_exclusive=2000,
                 included_by_trial_policy=True,assignment_basis=basis,proposed_sample_raw=proposed,
                 applied_shift_samples=None if proposed is None else proposed-sample)
            for i,(sample,proposed,basis) in enumerate([(1000,1019,"measured"),(1900,1919,"measured"),(1500,None,"unresolved")])]


def metadata():
    return dict(subject_bids_id="sub-11",version="test-measurement",mode="pilot",epoch_window_sec=[-.1,.09],
                method=copy.deepcopy(vu.DEFAULT_METHOD),stage2_report_sha256="source-stage2-hash")


def fast_graphics(monkeypatch):
    """Keep non-visual contract tests focused; separate tests inspect real axes."""
    tiny = b"not-a-scientific-figure"
    monkeypatch.setattr(vr,"lag_figures",lambda rows:{"delay_violin_with_dots.png":tiny,"lag_by_trigger.png":tiny})
    monkeypatch.setattr(vr,"time_figures",lambda rows,method:{name:tiny for name in ("lag_over_time.png","signal_usability.png","threshold_sensitivity.png")})
    monkeypatch.setattr(vr,"context_figure",lambda rows:(tiny,vr.context_statistics(rows)))
    monkeypatch.setattr(vr,"interval_figure",lambda rows:tiny)
    monkeypatch.setattr(vr,"interval_violin_figure",lambda rows:tiny)
    monkeypatch.setattr(vr,"interval_residual_figures",lambda rows:{name:tiny for name in (
        "interval_residual_vs_duration.png","interval_residual_vs_time.png","interval_residual_by_pair.png")})
    monkeypatch.setattr(vr,"coverage_figure",lambda rows:tiny)
    monkeypatch.setattr(vr,"measurement_coverage_figure",lambda rows:tiny)
    monkeypatch.setattr(vr,"donor_figures",lambda rows,estimates,subject:({"donor_estimates.png":tiny,"donor_validation.png":tiny},[]))
    cache = io.BytesIO(); np.savez_compressed(cache)
    monkeypatch.setattr(vr,"waveform_figures",lambda rows,waveforms:({"diagnostics.png":tiny,"waveform_windows.npz":cache.getvalue(),"diagnostic_index.csv":"array,event_key\n"},[],[]))


def test_statistics_use_sample_sd_and_preserve_missing():
    assert vr.stats([18,20,None,""])["sd"]==pytest.approx(np.sqrt(2))
    assert vr.stats([19])["sd"] is None
    assert vr.stats([19,19])["sd"]==0
    assert vr.stats([])["mean"] is None
    assert vr.number("nan") is None and vr.number("inf") is None


def test_lag_summary_never_pools_codes_recordings_runtimes_or_flags():
    rows=observations()
    result=vr.lag_statistics(rows)
    one=next(r for r in result if r["source_fif"]=="one.fif" and r["trigger_code"]==71)
    assert one["mean_lag_ms"]==19 and one["n_measured"]==2 and one["n_flagged"]==1
    assert one["sd_lag_ms"]==pytest.approx(np.sqrt(2))
    rows.append(dict(rows[0],event_key="new-runtime",runtime="v24",lag_ms=50))
    assert len(vr.lag_statistics(rows))==4


def test_bar_trigger_positions_are_factors_and_sd_is_sample_dispersion(monkeypatch):
    captured=[]
    monkeypatch.setattr(vr,"png",lambda fig:captured.append(fig) or b"figure")
    vr.lag_figures(vr.normalize(observations()))
    bar=captured[1].axes[0]
    assert [t.get_text() for t in bar.get_xticklabels()]==["71","118"]
    assert [p.get_x()+p.get_width()/2 for p in bar.patches]==pytest.approx([0,1])
    assert [p.get_height() for p in bar.patches]==[19,8]
    assert len(bar.containers)>=3
    for figure in captured:vr.pyplot().close(figure)


def test_sparse_constant_and_missing_groups_do_not_invent_violin_density(monkeypatch):
    rows=vr.normalize(observations()[:2])
    rows[1]["lag_ms"]=18
    monkeypatch.setattr(vr.pyplot().Axes,"violinplot",lambda *args,**kwargs:pytest.fail("No density for constant/sparse group"))
    files=vr.lag_figures(rows)
    assert files["lag_by_trigger.png"].startswith(b"\x89PNG")


def test_lag_display_limits_use_group_spread_and_retain_all_finite_tails():
    rows=vr.normalize(observations())
    rows=[r for r in rows if r["lag_ms"]>0]
    low,high=vr.lag_axis_limits(rows)
    assert 0<low<8 and high>30
    assert low<19-2*np.sqrt(2) and high>19+2*np.sqrt(2)
    one=vr.normalize(observations()[:2])
    low,high=vr.lag_axis_limits(one)
    assert low>15 and high<23


def test_recordings_without_finite_lag_do_not_create_blank_panels(monkeypatch):
    captured=[]
    rows=vr.normalize(observations())
    rows.append(dict(rows[0],source_fif="confirmed_unusable.fif",lag_ms=None,measurement_status="confirmed_unusable"))
    monkeypatch.setattr(vr,"png",lambda fig:captured.append(fig) or b"figure")
    vr.lag_figures(rows)
    assert len(captured[0].axes)==2 and "confirmed_unusable.fif" in captured[0].texts[-1].get_text()
    for fig in captured:vr.pyplot().close(fig)


def test_dense_group_draws_a_real_violin_with_dispersion():
    rows=vr.normalize(observations()[:1])
    rows=[dict(rows[0],event_key=f"dense-{i}",lag_ms=18+i%4) for i in range(24)]
    files=vr.lag_figures(rows)
    assert files["delay_violin_with_dots.png"].startswith(b"\x89PNG")


def test_coverage_rechecks_native_bounds_unassigned_and_arithmetic():
    maps=vr.normalize(trial_map())
    coverage,checks=vr.coverage_and_checks(maps,[-.1,.09])
    assert coverage[1]["original_inside"] is True and coverage[1]["proposed_inside"] is False
    assert coverage[2]["proposed_inside"] is None
    assert checks["assigned_arithmetic_mismatches"]==0
    maps[1]["applied_shift_samples"]=99
    assert vr.coverage_and_checks(maps,[-.1,.09])[1]["assigned_arithmetic_mismatches"]==1


def test_collision_and_order_checks_do_not_change_maps():
    maps=vr.normalize(trial_map())
    maps[1].update(proposed_sample_raw=1019,applied_shift_samples=-881)
    before=copy.deepcopy(maps)
    assert len(vr.coverage_and_checks(maps,[-.1,.09])[1]["eligible_proposed_sample_collisions"])==1
    maps[1].update(proposed_sample_raw=1000,applied_shift_samples=-900)
    assert len(vr.coverage_and_checks(maps,[-.1,.09])[1]["eligible_ordering_inversions"])==1
    assert before[0]==maps[0]


def test_context_missing_is_distinct_from_no():
    result=vr.context_statistics(vr.normalize(observations()))
    assert any(r["context_factor"]=="boundary_trial" and r["context_level"]=="unknown" for r in result)
    assert {r["context_level"] for r in result if r["context_factor"]=="guard_retry"}=={"yes","no"}


def test_report_embeds_each_figure_explanation_evidence_and_actual_findings(monkeypatch):
    fast_graphics(monkeypatch)
    rows,maps=observations(),trial_map()
    before=copy.deepcopy((rows,maps))
    files=vr.build_report_artifacts(rows,maps,{},[],[],{},metadata(),{})
    report=files["report.html"]
    summary=json.loads(files["report_summary.json"])
    assert "data:image/png;base64," in report and "download='lag_statistics.csv'" in report
    assert "<script" not in report and "https://" not in report
    for figure in summary["figure_explanations"]:
        assert all(figure[key] for key in ("axes","derived","finding","decision","example"))
    assert "PILOT ONLY" in report and "1.41 ms" in report
    assert "physical visual latency stays unknown" in report
    assert (rows,maps)==before
    assert "1000,71" in files["trial_timing_map.csv"]
    assert "1000.0,71.0" not in files["trial_timing_map.csv"]


def test_p11_p24_and_no_diode_guidance_preserves_constraints():
    assert "starting at block 3" in " ".join(vr.participant_guidance("11",metadata(),{}))
    assert "Only acquisition-5" in " ".join(vr.participant_guidance("24",metadata(),{}))
    assert "historical v20 counterfactual is not a participant calibration" in " ".join(vr.participant_guidance("04",metadata(),{}))
    assert "physical visual latency" in " ".join(vr.participant_guidance("17",metadata(),{}))


def test_revised_eligibility_report_separates_original_and_final_counts(monkeypatch):
    fast_graphics(monkeypatch)
    maps=trial_map()
    maps[0].update(final_analysis_eligible=False)
    maps[1].update(included_by_trial_policy=False,final_analysis_eligible=True)
    maps[2].update(final_analysis_eligible=False)
    files=vr.build_report_artifacts(observations(),maps,{},[],[],{},metadata(),{})
    summary=json.loads(files['report_summary.json'])
    assert summary['original_stage2_eligible_trial_keys']==2 and summary['eligible_trial_keys']==1
    assert 'original Stage 2 eligible: 2' in files['report.html']
    assert 'Applicable analysis eligible: 1' in files['report.html']
    assert sum(summary['eligible_assignment_counts'].values())==1


@pytest.mark.parametrize("subject",["11","24"])
def test_estimated_donor_figures_retain_restricted_sources_and_stored_uncertainty(subject,monkeypatch):
    rows=[]
    for i in range(20):
        rows.append(dict(observations()[0],event_key=f"donor-{i}",subject_bids_id=f"sub-{subject}",
                         source_fif="later-A.fif" if subject=="24" or i<10 else "later-B.fif",
                         acquisition_run=5,behavioral_block=3,behavioral_run=1+i%2,lag_ms=18+i%3))
    estimate=vu.estimate_lag(rows)
    assert estimate["status"]=="estimated"
    before=copy.deepcopy(estimate)
    rows.append(dict(rows[0],event_key="forbidden",source_fif="forbidden.fif",behavioral_block=1,acquisition_run=1,lag_ms=-40))
    captured=[]
    monkeypatch.setattr(vr,"png",lambda fig:captured.append(fig) or b"figure")
    _,summaries=vr.donor_figures(vr.normalize(rows),{f"P{subject}_code71":estimate},subject)
    labels=captured[0].axes[0].get_legend_handles_labels()[1]
    assert "forbidden.fif" not in labels
    assert summaries[0]["lag_estimate_ms"]==estimate["lag_estimate_ms"] and summaries[0]["donor_sd_ms"]==estimate["donor_sd_ms"]
    assert estimate==before
    for fig in captured:vr.pyplot().close(fig)


def test_cached_waveforms_retain_native_grid_and_anchor(tmp_path):
    t=(np.arange(125)-37)*2
    np.savez_compressed(tmp_path/"waveform_windows.npz",one=np.vstack([t,np.zeros(125),np.ones(125)]))
    (tmp_path/"diagnostic_index.csv").write_text("array,event_key\none,event\n")
    result=vr.load_cached_waveforms(tmp_path)["event"]
    assert result[1:]==(37,500.) and result[0].shape==(2,125)
    t[100]+=1
    np.savez_compressed(tmp_path/"waveform_windows.npz",one=np.vstack([t,np.zeros(125),np.ones(125)]))
    with pytest.raises(RuntimeError,match="sample grid"):vr.load_cached_waveforms(tmp_path)


def test_report_only_verifies_hashes_never_opens_raw_and_protects_source(tmp_path,monkeypatch):
    fast_graphics(monkeypatch)
    raw=tmp_path/"Attention/data";raw.mkdir(parents=True)
    bundle=tmp_path/"Attention/derivatives/visual_timing/source/sub-11/ses-01"
    cache=io.BytesIO();np.savez_compressed(cache)
    contents={"event_measurements.csv":csv_text(observations()),"trial_timing_map.csv":csv_text(trial_map()),
              "event_timing_proposals.csv":"event_key,assignment_basis\n", "recording_summary.csv":"source_fif\n",
              "physical_intervals.csv":"left_event_key\n", "estimates.json":"{}","contexts.json":"{}",
              "report.md":"old report","diagnostics.png":b"old", "waveform_windows.npz":cache.getvalue(),"diagnostic_index.csv":"array,event_key\n"}
    vt.publish(bundle,contents,dict(metadata(),schema_version="1.0",status="completed"),raw)
    original={p.name:file_hash(p) for p in bundle.iterdir()}
    monkeypatch.setattr(vt.mne.io,"read_raw_fif",lambda *a,**k:pytest.fail("Report rendering must not read raw FIF"))
    monkeypatch.setattr(vu,"measure_edge",lambda *a,**k:pytest.fail("Report rendering must not remeasure"))
    output=tmp_path/"Attention/derivatives/reports/visual_timing/new/sub-11/ses-01"
    result=vr.render_existing_bundle(bundle,output,raw)
    assert result["source"]["raw_signal_windows_read"]==0 and (output/"report.html").exists()
    assert original=={p.name:file_hash(p) for p in bundle.iterdir()}
    document=(output/"report.html").read_text()
    embedded=re.search(r"download='trial_timing_map.csv' href='data:[^;]+;base64,([^']+)'",document)
    assert base64.b64decode(embedded.group(1))==(bundle/"trial_timing_map.csv").read_bytes()
    assert (output/"source_manifest.json").read_bytes()==(bundle/"manifest.json").read_bytes()
    with pytest.raises(ValueError,match="new report directory"):vr.render_existing_bundle(bundle,output,raw)
    (bundle/"report.md").write_text("edited")
    with pytest.raises(RuntimeError,match="fingerprint"):vr.render_existing_bundle(bundle,output.parent/"other",raw)
