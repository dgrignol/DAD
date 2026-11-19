"""Check the restricted Stage-3 eligibility/donor application and consumer gate.

Usage: PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python3 -B -m pytest \
  -p no:cacheprovider meg/tests/test_visual_timing_application.py -q

Synthetic exact-key ledgers/MAT metadata test saved-class inclusion, independent
exclusions, all eight decided trial identities, P24-only donor removal and
hash-bound review acceptance. No participant data or cohort stages execute.
The epoching boundary uses a synthetic fixture and intercepts MNE before epoch
construction; original and final eligibility must remain separate there.
"""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest

CODE=Path(__file__).resolve().parents[1]/"code/preprocessing_basic"
sys.path.insert(0,str(CODE))
import visual_timing_application as app
import visual_timing_utils as utils
import visual_timing_report as reporting
import visual_timing as timing
from trial_preparation import csv_text,digest,file_hash,read_csv,stage2_paths


# Saved class, source context and independent exclusions are explicit inputs.
def row(sample=1000,*,subject="10",source="one.fif",condition=81,reasons="",block=1,run=1,trial=1,replay=False):
    return dict(trial_key=f"sub-{subject}:{source}:{sample}",subject_bids_id=f"sub-{subject}",source_fif=source,
        raw_sample_frame=sample,sample_start_raw=sample,condition_code=condition,acquisition_run=block,behavioral_block=block,
        behavioral_run=run,behavioral_trial=trial,is_replay_trial=replay,reported_source_csv="saved.csv",
        trial_sequence_identity=4,included_by_trial_policy=not reasons,exclusion_reasons=reasons,
        accepted_exclusion_reason="",sfreq_hz=1000.,recording_first_sample_raw=0,
        recording_end_sample_raw_exclusive=10000,trial_interval_end_raw_exclusive=10000)


def metadata(original,catch=0):
    return dict(source_mat="saved.mat",behavioral_block=original["behavioral_block"],behavioral_run=original["behavioral_run"],
        behavioral_trial=original["behavioral_trial"],is_replay_trial=original["is_replay_trial"],
        saved_condition_code=original["condition_code"],saved_catch_type_code=catch,debug_source_csv="saved.csv")


@pytest.mark.parametrize("catch,condition,reasons,expected",[
    (2,81,"catch_type_2",True),(2,91,"catch_type_2",True),
    (2,81,"catch_type_2,gaze_break",False),(2,91,"catch_type_2,esc_termination",False),
    (2,81,"catch_type_2,stage1_reviewed_interrupted_trial",False),
    (1,71,"catch_type_1",False),(2,71,"catch_type_2",False),(0,81,"catch_type_2",False)])
def test_only_saved_occlusion_type2_can_waive_one_reason(catch,condition,reasons,expected):
    original=row(condition=condition,reasons=reasons);before=copy.deepcopy(original)
    overlay=app.build_eligibility_overlay([original],[metadata(original,catch)],"10")[0]
    assert overlay["final_analysis_eligible"]==expected
    assert overlay["original_exclusion_reasons"]==reasons and original==before
    assert overlay["original_included_by_trial_policy"] is False
    if "gaze_break" in reasons:assert overlay["analysis_exclusion_reasons"]=="gaze_break"


@pytest.mark.parametrize("case",app.EXCLUSIONS)
def test_every_exact_sequence_exclusion_overrides_catch_inclusion(case):
    subject,source,sample,code,acquisition,block,run,trial,replay=case
    original=row(sample,subject=subject,source=source,condition=code,block=block,run=run,trial=trial,replay=replay,reasons="catch_type_2")
    original["acquisition_run"]=acquisition
    overlay=app.build_eligibility_overlay([original],[metadata(original,2)],subject,require_exact_exclusions=False)[0]
    assert not overlay["final_analysis_eligible"] and overlay["sequence_exclusion"]
    assert overlay["analysis_exclusion_reasons"]=="stage3_sequence_exclusion"
    assert app.policy_definition()["donor_excluded_event_keys"]==list(app.P24_DONOR_KEYS)


def test_missing_metadata_and_question_evidence_never_create_inclusion():
    original=row(reasons="catch_type_2,esc_termination")
    original["evidence_trigger_codes"]="81,118,201"
    result=app.build_eligibility_overlay([original],[],"10")[0]
    assert not result["final_analysis_eligible"] and result["saved_catch_type_code"] is None
    assert result["analysis_exclusion_reasons"]==original["exclusion_reasons"]
    assert result["saved_metadata_join"]=="unavailable_original_disposition_preserved"
    assert app.portable_basename('.\\output_files\\saved.csv')=='saved.csv'


def test_wrong_saved_condition_or_duplicate_context_is_rejected():
    original=row(reasons="catch_type_2");saved=metadata(original,2)
    with pytest.raises(RuntimeError,match="condition"):
        app.build_eligibility_overlay([original],[dict(saved,saved_condition_code=91)],"10")
    with pytest.raises(RuntimeError,match="Duplicate saved"):
        app.build_eligibility_overlay([original],[saved,saved],"10")


def test_saved_reader_ignores_empty_planned_runs_and_normalizes_lab_paths(tmp_path):
    from scipy.io import savemat
    folder=tmp_path/'output_files';folder.mkdir()
    path=folder/'saved.mat'
    trial=dict(block_num=1,trial_num=1,is_replay_trial=0,condition_trigger=81,catch_type_code=2,
               sequence_trigger_code=4,trial_completed=1,frame_missed_count=0)
    runs=np.asarray([dict(run_num=1,output=trial,debugTriggerFile='.\\output_files\\saved.csv'),
                     dict(run_num=np.asarray([]),output=np.asarray([]),debugTriggerFile='')],dtype=object)
    savemat(path,{'runResults':runs})
    cfg={'layout':{'subject_bids_id':'sub-10'},'resolved_paths':{'subject_data_dir':str(tmp_path)}}
    result=app.read_saved_metadata(cfg,{'saved_mat_sha256':{'saved.mat':file_hash(path)}},{'1':{'source_mat':'saved.mat'}})
    assert len(result)==1 and result[0]['saved_catch_type_code']==2 and result[0]['debug_source_csv']=='saved.csv'


def test_overlay_preserves_original_flag_and_reporting_uses_final_flag():
    original=row(reasons="catch_type_2");overlay=app.build_eligibility_overlay([original],[metadata(original,2)],"10")
    original.update(assignment_basis="measured",applied_shift_samples=19,proposed_sample_raw=1019)
    mapped=app.attach_eligibility([original],overlay)[0]
    assert mapped["included_by_trial_policy"] is False and mapped["final_analysis_eligible"] is True
    assert utils.analysis_eligible(mapped)
    assert reporting.correction_statistics([mapped])[0]["n_eligible"]==1
    coverage,_=reporting.coverage_and_checks([mapped],[-.5,2.667])
    assert coverage[0]["eligible"] is True


# P24's five explicit measured keys are the only donor removals.
def test_p24_rebuild_retains_measured_evidence_and_all_other_donor_members():
    method=dict(utils.DEFAULT_METHOD,bootstrap_resamples=20)
    measurements=[];estimates={}
    for code in (81,111,112,114,115):
        excluded=next(k for k in app.P24_DONOR_KEYS if k.endswith(f":{code}"))
        n=49 if code==81 else 90
        group=[dict(event_key=excluded,subject_bids_id="sub-24",source_fif="sub24_run05_100926.fif",
                    acquisition_run=5,behavioral_block=5,behavioral_run=1,trigger_code=code,
                    measurement_status="measured",lag_ms=15 if code==81 else 19)]
        group.extend(dict(group[0],event_key=f"keep-{code}-{i}",lag_ms=19,final_analysis_eligible=False) for i in range(n-1))
        measurements.extend(group);estimates[f"P24_code{code}"]=utils.estimate_lag(group,method)
    for code in (116,117):estimates[f"P24_code{code}"]={"status":"insufficient_donors","n_donors":6}
    estimates["P24_code71"]={"status":"estimated","lag_estimate_ms":19,"unchanged":"keep"}
    original=row(subject="24");mapped=dict(original,assignment_basis="estimated",donor_estimate_id="P24_code81",
        applied_shift_samples=19,proposed_sample_raw=1019,lag_estimate_ms=19)
    before=copy.deepcopy((measurements,estimates))
    maps,proposals,new,audit=app.rebuild_p24_donors([mapped],[],estimates,measurements,[original],method,[-.5,2.667])
    assert (measurements,estimates)==before and proposals==[]
    assert new['P24_code81']['n_donors']==48 and new['P24_code81']['lag_estimate_ms']==19
    assert new['P24_code81']['donor_sd_ms']<estimates['P24_code81']['donor_sd_ms']
    for code in (111,112,114,115):assert new[f'P24_code{code}']['n_donors']==89
    assert new['P24_code71']==estimates['P24_code71']
    assert new['P24_code116']==new['P24_code117']=={'status':'insufficient_donors','n_donors':6}
    assert all(not set(e.get('donor_event_keys',[])).intersection(app.P24_DONOR_KEYS) for e in new.values())
    assert len(audit)==5 and maps[0]['applied_shift_samples']==19


def reviewed_application_fixture(tmp_path):
    """Create genuinely hash-bound source/final bundles, with one new catch."""
    raw=tmp_path/'Attention/data';raw.mkdir(parents=True)
    derivatives=raw.parent/'derivatives';subject_data=raw/'sub10';(subject_data/'output_files').mkdir(parents=True)
    cfg={'subject':{'id':'10'},'layout':{'subject_bids_id':'sub-10','session_bids_id':'ses-01'},
         'epoching':{'tmin_sec':-.5,'tmax_sec':2.667},'resolved_paths':{'data_root':str(raw),'subject_data_dir':str(subject_data),
         'derivatives_root':str(derivatives),'derivatives_preprocessing_dir':str(derivatives/'preprocessing/sub-10/ses-01'),
         'derivatives_reports_dir':str(derivatives/'reports/sub-10/ses-01')}}
    original=row(reasons='catch_type_2');prepared=[original]
    paths=stage2_paths(cfg)
    for key in ('ledger','report'):paths[key].parent.mkdir(parents=True,exist_ok=True)
    paths['ledger'].write_text(csv_text(prepared));paths['report'].write_text('{}')
    saved=[metadata(original,2)];overlay=app.build_eligibility_overlay(prepared,saved,'10')
    measurement=dict(event_key='observed',measurement_status='measured',source_fif='one.fif',
                     digital_sample_raw=1000,trigger_code=81,diode_edge_sample_raw=1019)
    mapped=dict(original,assignment_basis='measured',assignment_status='proposed',proposed_sample_raw=1019,
                applied_shift_samples=19,lag_estimate_ms=19,uncertainty_kind='native_grid',measurement_event_key='observed')
    contents={'trial_timing_map.csv':csv_text([mapped]),'event_measurements.csv':csv_text([measurement]),
              'event_timing_proposals.csv':'event_key\n','recording_summary.csv':'source_fif\n','estimates.json':'{}',
              'contexts.json':'{}','report.md':'synthetic','diagnostics.png':b'synthetic'}
    manifest=dict(schema_version='1.0',status='completed',mode='cohort_proposal',subject_bids_id='sub-10',
                  version='s3-v001',stage2_report_sha256=file_hash(paths['report']),epoch_window_sec=[-.5,2.667])
    source=derivatives/'visual_timing/s3-v001/sub-10/ses-01'
    timing.publish(source,contents,manifest,raw)
    saved_text,overlay_text=csv_text(saved),csv_text(overlay)
    application=dict(schema_version='1.0',policy=app.policy_definition(),subject_bids_id='sub-10',target_version='s3-v002',
        source_version='s3-v001',source_manifest_sha256=file_hash(source/'manifest.json'),
        stage2_report_sha256=file_hash(paths['report']),stage2_ledger_sha256=file_hash(paths['ledger']),saved_mat_sha256={},
        saved_metadata_sha256=digest('placeholder'),eligibility_overlay_sha256=digest('placeholder'),donor_excluded_event_keys=[],
        scientific_review_requirements=['final_analysis_eligibility','newly_included_catch_onsets'])
    import hashlib
    application.update(saved_metadata_sha256=hashlib.sha256(saved_text.encode()).hexdigest(),
                       eligibility_overlay_sha256=hashlib.sha256(overlay_text.encode()).hexdigest())
    application_text=json.dumps(application);application_sha=hashlib.sha256(application_text.encode()).hexdigest()
    contents.update({'trial_timing_map.csv':csv_text(app.attach_eligibility([mapped],overlay)),
                     'saved_trial_metadata.csv':saved_text,'eligibility_overlay.csv':overlay_text,'application.json':application_text})
    final=derivatives/'visual_timing/s3-v002/sub-10/ses-01'
    timing.publish(final,contents,dict(manifest,version='s3-v002',task1_staging=False,review_application_sha256=application_sha),raw)
    review=tmp_path/'accepted-synthetic-review.json'
    review.write_text(json.dumps(dict(schema_version='1.0',status='accepted',reviewer='synthetic test only',
        reviewed_utc='2026-10-21T00:00:00Z',decision_reference='synthetic test',bundle_manifest=str(final/'manifest.json'),
        bundle_manifest_sha256=file_hash(final/'manifest.json'),accepted_assignment_bases=['measured'],
        accepted_eligibility_application_sha256=application_sha,accepted_review_scopes=application['scientific_review_requirements'])))
    return cfg,prepared,final,review


def test_later_working_log_update_does_not_change_frozen_application(tmp_path):
    cfg,prepared,folder,review=reviewed_application_fixture(tmp_path)
    live=Path(cfg['resolved_paths']['derivatives_root'])/'visual_timing/s3-v001/review_notes.md'
    live.write_text('Later authorized human reconciliation of the mutable working log.')
    result,_=utils.load_reviewed_timing(cfg,prepared,review)
    assert result[0]['final_analysis_eligible'] is True


def test_hash_bound_consumer_keeps_original_false_and_selects_final_true(tmp_path):
    cfg,prepared,folder,review=reviewed_application_fixture(tmp_path)
    result,gate=utils.load_reviewed_timing(cfg,prepared,review)
    assert result[0]['included_by_trial_policy'] is False
    assert result[0]['final_analysis_eligible'] is True and result[0]['epoch_event_sample_raw']==1019
    assert gate['eligibility_application']['application_sha256']==file_hash(folder/'application.json')
    # Exercise the production outcome boundary with a fake already-selected epoch.
    spec=importlib.util.spec_from_file_location('application_epoch_test',CODE/'04_epoching.py')
    epoching=importlib.util.module_from_spec(spec);spec.loader.exec_module(epoching)
    epoching._mark_mne_outcomes(result,result,SimpleNamespace(selection=[0],drop_log=[()]))
    assert result[0]['included_final'] is True
    removed=dict(result[0],included_by_trial_policy=True,final_analysis_eligible=False)
    epoching._mark_mne_outcomes([removed],[],SimpleNamespace(selection=[],drop_log=[]))
    assert removed['included_final'] is False and removed['final_epoch_inclusion']=='excluded_before_epoching'


def test_mandatory_review_scopes_come_from_overlay_and_participant():
    row={'newly_included':True}
    assert app.mandatory_review_scopes([row],'10')=={
        'final_analysis_eligibility','newly_included_catch_onsets'}
    assert app.mandatory_review_scopes([{'newly_included':False}],'33')=={
        'final_analysis_eligibility','p33_remaining_original_onsets'}


def test_hash_bound_application_cannot_omit_required_scope(tmp_path):
    cfg,prepared,folder,review_path=reviewed_application_fixture(tmp_path)
    application_path=folder/'application.json'
    payload=json.loads(application_path.read_text())
    payload['scientific_review_requirements']=[]
    application_path.write_text(json.dumps(payload))
    manifest_path=folder/'manifest.json';bundle=json.loads(manifest_path.read_text())
    sha=file_hash(application_path)
    bundle['review_application_sha256']=sha
    bundle['output_sha256']['application.json']=sha
    manifest_path.write_text(json.dumps(bundle,sort_keys=True))
    review=json.loads(review_path.read_text())
    review['bundle_manifest_sha256']=file_hash(manifest_path)
    review['accepted_eligibility_application_sha256']=sha
    review['accepted_review_scopes']=['final_analysis_eligibility','newly_included_catch_onsets']
    with pytest.raises(RuntimeError,match='omits a mandatory scientific review scope'):
        app.load_applied_eligibility(cfg,prepared,folder,bundle,review)


def test_production_epoch_event_boundary_uses_final_eligibility(tmp_path,monkeypatch):
    """Intercept the real consumer before MNE construction; no signal is read."""
    cfg,prepared,folder,review=reviewed_application_fixture(tmp_path)
    cfg['epoching'].update(require_stage2_ledger=True,require_stage3_timing=True,
        event_id_map={'occluded_nondeviant':81},baseline_sec=None,detrend=1,reject_by_annotation=True,
        amplitude_reject={},trial_exclusion={})
    cfg['subject_profile']={'identity_range':[1,80],'iti_sec_range':[2,5]}
    cfg['resolved_paths']['derivatives_analysis_dir']=str(tmp_path/'Attention/derivatives/analysis/sub-10/ses-01')
    spec=importlib.util.spec_from_file_location('application_real_event_boundary',CODE/'04_epoching.py')
    epoching=importlib.util.module_from_spec(spec);spec.loader.exec_module(epoching)
    runs=[SimpleNamespace(cleaned_path='unused-synthetic-path',source=SimpleNamespace(source_filename='one.fif'))]
    # Stub the current preparation boundary, preserving the real eligibility
    # application and prospective event/metadata assertions below.
    monkeypatch.setattr(epoching,'resolve_epoch_inputs',lambda *a:(
        runs,{'status':'passed','preparation':{'mode':'saved_ica'}},None))
    monkeypatch.setattr(epoching,'load_trigger_event_rows',lambda *a:[dict(prepared[0],trigger_code=81,event_index=1)])
    monkeypatch.setattr(epoching,'load_prepared_decisions',lambda *a:copy.deepcopy(prepared))
    # The epoch preparation contract binds the canonical table's bytes even
    # when this event-boundary fixture supplies its parsed rows directly.
    events_path=epoching.epoching_paths(cfg)['events_csv']
    events_path.parent.mkdir(parents=True,exist_ok=True)
    events_path.write_text(csv_text([dict(prepared[0],trigger_code=81,event_index=1)]))
    # Native inputs always carry annotations, including an explicit empty set.
    monkeypatch.setattr(epoching.mne.io,'read_raw_fif',lambda *a,**k:SimpleNamespace(first_samp=0,n_times=10000,info={'sfreq':1000.},annotations=epoching.mne.Annotations([],[],[])))
    monkeypatch.setattr(epoching,'prepare_run',lambda cfg,run,*a:(
        epoching.mne.io.read_raw_fif(run.cleaned_path),None))
    class ReachedEventBoundary(Exception):pass
    def inspect(raw,events,**kwargs):
        assert events.tolist()==[[1019,0,81]]
        values=kwargs['metadata'].to_dict(orient='records')[0]
        assert values['included_by_trial_policy'] is False and values['final_analysis_eligible'] is True
        assert values['exclusion_reasons']=='catch_type_2' and values['analysis_exclusion_reasons']==''
        raise ReachedEventBoundary
    monkeypatch.setattr(epoching.mne,'Epochs',inspect)
    with pytest.raises(ReachedEventBoundary):
        epoching.run_epoching(cfg,source_manifest=folder/'source-manifest-unused.json',overwrite=False,timing_review=review)


@pytest.mark.parametrize('problem',['missing_application_acceptance','missing_new_catch_scope','edited_overlay','staging'])
def test_changed_or_unaccepted_application_cannot_feed_consumer(tmp_path,problem):
    cfg,prepared,folder,review=reviewed_application_fixture(tmp_path)
    accepted=json.loads(review.read_text())
    if problem=='missing_application_acceptance':accepted.pop('accepted_eligibility_application_sha256')
    elif problem=='missing_new_catch_scope':accepted['accepted_review_scopes']=['final_analysis_eligibility']
    elif problem=='edited_overlay':(folder/'eligibility_overlay.csv').write_text('changed')
    elif problem=='staging':
        manifest=json.loads((folder/'manifest.json').read_text());manifest['task1_staging']=True
        (folder/'manifest.json').write_text(json.dumps(manifest));accepted['bundle_manifest_sha256']=file_hash(folder/'manifest.json')
    review.write_text(json.dumps(accepted))
    with pytest.raises(RuntimeError):utils.load_reviewed_timing(cfg,prepared,review)

