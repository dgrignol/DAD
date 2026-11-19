"""Numerical checks for an isolated ICA-convergence candidate, not acceptance.

Usage: python -B -m pytest -q meg/tests/test_ica_convergence_reconciliation.py
Known independent non-Gaussian sources exercise real MNE/sklearn whitening and
the fixed-point measure; no participant/model/receipt or native file is touched.
"""
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import mne
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code/preprocessing_basic'))
from reconcile_ica_convergence import fixed_point_residual, model_comparison
from publish_ica_convergence_candidate import require_candidate
from copy import deepcopy
from ica_optimizer_checkpoint import iteration_cap, load_checkpoint, save_checkpoint
from diagnose_ica_convergence import alignment, diagnose_arrays, fitted_basis
from probe_ica_curvature_restart import represented_model, run_restart
from prepare_ica_restart_candidate import portable_summaries


def example():
    """Return a reproducible converged ICA of mixed independent Laplace sources."""
    rng = np.random.RandomState(97)
    sources = rng.laplace(size=(3, 12000))
    mixing = np.array([[1., .2, .4], [.3, 1., .1], [.1, .4, 1.]])
    raw = mne.io.RawArray(mixing @ sources * 1e-12,
        mne.create_info(['MEG0111','MEG0121','MEG0131'],200.,'mag'),verbose='ERROR')
    raw.filter(1.,40.,verbose='ERROR')
    model = mne.preprocessing.ICA(method='fastica',n_components=None,random_state=97,max_iter=2000)
    model.fit(raw,verbose='ERROR')
    return raw,model


def test_fixed_point_check_uses_the_fitted_whitening_and_leaves_data_unchanged():
    raw,model = example()
    before = raw.get_data().copy()
    weights = model.unmixing_matrix_.copy()
    residual = fixed_point_residual(model,raw)
    assert residual['max'] < residual['tol']
    np.testing.assert_array_equal(raw.get_data(),before)
    np.testing.assert_array_equal(model.unmixing_matrix_,weights)


def test_fit_summary_portability_keeps_every_native_interval_and_filter_field():
    cluster=[{'input_path':'/mnt/storage/derivatives/run_ann.fif',
        'intervals':[{'start_raw':12000,'stop_raw_exclusive':15000,'fit_samples':600}],
        'filter':{'n_taps':3301},'fit_samples':600}]
    local=deepcopy(cluster);local[0]['input_path']='/Volumes/MORWUR/derivatives/run_ann.fif'
    assert portable_summaries(cluster)==portable_summaries(local)
    local[0]['intervals'][0]['start_raw']+=1
    assert portable_summaries(cluster)!=portable_summaries(local)
    assert cluster[0]['input_path'].startswith('/mnt/')


def test_model_comparison_matches_identical_components_without_exclusion_choices(tmp_path):
    raw,model = example()
    comparison = model_comparison(model,model,raw,tmp_path)
    assert comparison['minimum_activation_match'] > 1.-1e-12
    assert all(row['original_IC']==row['candidate_IC'] for row in comparison['matches_for_comparison_only'])
    assert model.exclude == [] and (tmp_path/'model-comparison.png').is_file()
    plt.close('all')


def test_fixed_point_check_refuses_a_different_optimizer():
    raw,model = example()
    model.fit_params['algorithm']='deflation'
    with pytest.raises(RuntimeError,match='parallel/logcosh'):
        fixed_point_residual(model,raw)


def test_direction_diagnosis_distinguishes_label_exchange_from_rotation():
    original=np.eye(3)
    exchanged=original[[1,0,2]]
    assert alignment(exchanged,original)['rowwise_max_residual']==1.
    assert alignment(exchanged,original)['permutation_matched_max_residual']==0.
    angle=np.pi/4
    rotation=np.array([[np.cos(angle),-np.sin(angle),0.],
                       [np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])
    assert alignment(rotation,original)['permutation_matched_max_residual']>.2
    # A true two-cycle returns to the original directions, independently of signs.
    np.testing.assert_allclose(alignment(original,original)['rowwise_max_residual'],0.)


def test_converged_real_ica_has_small_recurrence_and_unchanged_inputs():
    raw,model=example()
    data=raw.get_data().copy();weights=model.unmixing_matrix_.copy()
    white,matrix=fitted_basis(model,raw)
    result=diagnose_arrays(white,matrix,model.fit_params,steps=3)
    assert result['whitening_covariance_max_error']<1e-10
    assert result['unmixing_orthogonality_max_error']<1e-10
    assert result['recurrence'][0]['one_step']['rowwise_max_residual']<result['tol']
    assert len(result['component_descriptors'])==3
    assert result['scientific_acceptance'] is False and result['component_decisions'] is None
    np.testing.assert_array_equal(raw.get_data(),data)
    np.testing.assert_array_equal(model.unmixing_matrix_,weights)


def test_restart_observer_matches_stock_optimizer_and_mne_representation(tmp_path):
    """Verify actual optimizer output, unchanged inputs and the PCA normalization."""
    from sklearn.decomposition import FastICA
    raw,original=example()
    before=original.unmixing_matrix_.copy();data=raw.get_data().copy()
    white,matrix=fitted_basis(original,raw)
    parameters=original.fit_params
    estimator,initial,eigenvalues,history,caught,elapsed=run_restart(
        white,matrix,parameters,max_iter=1000,progress_dir=tmp_path)
    stock=FastICA(n_components=None,whiten=False,random_state=97,max_iter=1000,w_init=initial,
        **{k:v for k,v in parameters.items() if k not in ('max_iter','w_init')})
    stock.fit(white.T)
    np.testing.assert_array_equal(estimator.components_,stock.components_)
    assert estimator.n_iter_==stock.n_iter_<1000
    assert not any(row['convergence_warning'] for row in caught)
    candidate=represented_model(original,estimator.components_,estimator.n_iter_,1000)
    np.testing.assert_allclose(candidate.get_sources(raw).get_data(),estimator.components_@white,atol=1e-12)
    assert fixed_point_residual(candidate,raw)['max']<1e-4
    reconstructed=raw.copy();candidate.apply(reconstructed,exclude=[],verbose='ERROR')
    np.testing.assert_allclose(reconstructed.get_data(),data,atol=1e-25)
    np.testing.assert_array_equal(original.unmixing_matrix_,before)
    np.testing.assert_array_equal(raw.get_data(),data)
    assert candidate.n_components is None and candidate.exclude==[]


def test_publication_gate_allows_only_a_proven_iteration_cap_change():
    raw,original = example()
    candidate = deepcopy(original)
    candidate.max_iter=10000; candidate.fit_params['max_iter']=10000
    old_fit={'method':'fastica','max_iter':'auto','random_state':97,'n_components':None}
    new_fit={**old_fit,'max_iter':10000}
    payload={'fingerprint':{'fit':old_fit,'inputs':[]},'fit_execution':{'n_samples_fitted':original.n_samples_}}
    metadata={'fingerprint':{'fit':new_fit,'inputs':[]}}
    result={'numeric_candidate_ok':True,'convergence_warnings':[],
        'candidate_residual':fixed_point_residual(candidate,raw),
        'candidate_iterations':candidate.n_iter_,'no_removal_reconstruction_relative_error':0.}
    require_candidate(result,metadata,candidate,original,payload)
    rejected=deepcopy(result); rejected['convergence_warnings']=['did not converge']
    with pytest.raises(RuntimeError,match='numerical convergence'):
        require_candidate(rejected,metadata,candidate,original,payload)
    metadata['fingerprint']['fit']['random_state']=98
    with pytest.raises(RuntimeError,match='Only the iteration cap'):
        require_candidate(result,metadata,candidate,original,payload)


def test_checkpoint_survives_source_failure_without_acceptance_or_refitting(tmp_path,monkeypatch):
    """Retain actual weights while keeping all source acceptance explicitly pending."""
    from continuous_cleaning import file_hash
    raw,model = example()
    context = {'subject':19,'original_model_sha256':'original-fixture',
               'candidate_settings':{'random_state':97,'max_iter':2000}}
    directory = tmp_path/'new-local-checkpoint'
    warnings = [{'category':'ConvergenceWarning','message':'fixture warning retained',
                 'convergence_warning':True}]
    metadata = save_checkpoint(directory,model,context,warnings,1.5,tmp_path/'raw')
    assert metadata['source_revalidation_required'] is True
    assert metadata['scientific_acceptance'] is False
    before = file_hash(directory/'checkpoint_ica.fif')
    # Simulate a later network-source failure; the helper has no permission to
    # accept sources and needs none of them to preserve already computed bytes.
    with pytest.raises(OSError):
        raise OSError('shared source disconnected after computation')
    monkeypatch.setattr(mne.preprocessing.ICA,'fit',lambda *a,**kw:pytest.fail('Resume must not refit ICA'))
    restored,record = load_checkpoint(directory,context)
    assert iteration_cap(restored) == 2000
    np.testing.assert_array_equal(restored.unmixing_matrix_,model.unmixing_matrix_)
    assert record['warnings'] == warnings and record['source_revalidation_required']
    assert file_hash(directory/'checkpoint_ica.fif') == before
    changed = deepcopy(context); changed['candidate_settings']['random_state'] = 98
    with pytest.raises(RuntimeError,match='context changed'):
        load_checkpoint(directory,changed)
    with pytest.raises(FileExistsError):
        save_checkpoint(directory,model,context,[],1.5,tmp_path/'raw')
    with (directory/'checkpoint_ica.fif').open('ab') as handle:handle.write(b'changed')
    with pytest.raises(RuntimeError,match='bytes changed'):
        load_checkpoint(directory,context)


def test_restart_publication_requires_exact_pca_ancestry_and_full_source_validation():
    """A real converged solution is eligible only with explicit restart history."""
    raw,original=example();candidate=deepcopy(original)
    fingerprint={'fit':{'method':'fastica','max_iter':'auto','random_state':97,'n_components':None},'inputs':[]}
    history={'kind':'deterministic_curvature_restart_v1','fresh_random_draw':False,
        'random_state_argument':97,'requested_n_components':None,
        'restart_phase_iterations':candidate.n_iter_,'original_phase_iterations':original.n_iter_,
        'iteration_cap':iteration_cap(original),'initializer_sha256':'fixture-init',
        'diagnosis_sha256':'fixture-diag'}
    payload={'fingerprint':fingerprint,'fit_execution':{'n_samples_fitted':original.n_samples_}}
    metadata={'fingerprint':deepcopy(fingerprint),'optimizer_restart':deepcopy(history)}
    result={'repair_kind':'deterministic_curvature_restart_v1','optimizer_restart':history,
        'numeric_candidate_ok':True,'convergence_warnings':[],
        'candidate_residual':fixed_point_residual(candidate,raw),'candidate_iterations':candidate.n_iter_,
        'no_removal_reconstruction_relative_error':0.,'native_auxiliary_samples_equal':True,
        'source_validation':'passed_before_and_after_candidate_validation'}
    require_candidate(result,metadata,candidate,original,payload)
    broken=deepcopy(result);broken['source_validation']='pending'
    with pytest.raises(RuntimeError,match='Restart provenance'):
        require_candidate(broken,metadata,candidate,original,payload)
    changed=deepcopy(metadata);changed['fingerprint']['fit']['random_state']=98
    with pytest.raises(RuntimeError,match='Restart provenance'):
        require_candidate(result,changed,candidate,original,payload)
    candidate.pca_mean_[0]+=1e-8
    with pytest.raises(RuntimeError,match='PCA/whitening'):
        require_candidate(result,metadata,candidate,original,payload)


# Publication is a filesystem transaction: an interrupted render must preserve
# old scientific evidence and unrelated work, even after the candidate install.
@pytest.mark.parametrize('restart',[False,True])
def test_failed_publication_restores_model_evidence_and_preserves_other_notes(tmp_path, monkeypatch,restart):
    import json
    import yaml
    import scipy
    import sklearn
    import publish_ica_convergence_candidate as publisher
    import review_final_ica as review
    from continuous_cleaning import file_hash
    from test_final_ica_review import saved_review

    args,cfg,source,model_path,metadata_path,tuning,dirs = saved_review(tmp_path,monkeypatch)
    root = Path(cfg['resolved_paths']['derivatives_root'])
    args.candidate_dir = root/'ica_convergence/sub-04/test-candidate'
    args.candidate_dir.mkdir(parents=True)
    args.archive_label = 'failed-render'
    args.config = 'unused-test-config'; args.site_config = None
    args.subject_config_dir = str(tmp_path)
    payload = json.loads((dirs['reports']/'03_C_ICA_review_report.json').read_text())
    original = mne.preprocessing.read_ica(model_path,verbose='ERROR')
    candidate = deepcopy(original)
    if not restart:
        candidate.max_iter = 10000; candidate.fit_params['max_iter'] = 10000
    candidate.save(args.candidate_dir/model_path.name,verbose='ERROR')
    stage = review.load_stage()
    candidate_cfg = deepcopy(cfg['ica'])
    if not restart:candidate_cfg['max_iter'] = 10000
    metadata = {'fingerprint':stage.build_ica_fingerprint(input_paths=[source],ica_cfg=candidate_cfg),
                'model_sha256':file_hash(args.candidate_dir/model_path.name)}
    (args.candidate_dir/metadata_path.name).write_text(json.dumps(metadata))
    result = {'subject':4,'original_model_sha256':file_hash(model_path),
        'candidate_model_sha256':metadata['model_sha256'],
        'changed_setting_only':{'max_iter':{'original':1000,'candidate':10000}},
        'runtime':{'mne':mne.__version__,'numpy':np.__version__,
                   'scipy':scipy.__version__,'sklearn':sklearn.__version__}}
    if restart:
        result.update(repair_kind='deterministic_curvature_restart_v1',supporting_file_hashes={},
            optimizer_restart={'seed_97_ancestry_original_model_sha256':file_hash(model_path),
                'candidate_model_sha256':metadata['model_sha256'],'candidate_directory':str(args.candidate_dir)})
    (args.candidate_dir/'result.json').write_text(json.dumps(result))
    native = dirs['figures']/'03_C_ICA_review_native_diagnostics.json'
    native.write_text(json.dumps({'retained_trial_metadata':[{'trial_key':'temporary-fixture'}]}))
    manifest = dirs['reports']/'source_manifest.json'
    old_stage = {'status':'original-unreviewed-fixture'}
    manifest.write_text(json.dumps({'derivative_lineage':{'stages':{'ica_manual':old_stage}}}))
    protected = [model_path,metadata_path,source,native,dirs['reports']/'03_C_ICA_review_report.json',
                 *dirs['figures'].glob('*.png')]
    before = {path:file_hash(path) for path in protected}
    old_ica_tuning = yaml.safe_load(tuning.read_text()).get('ica')
    current = review.triage.subject_config
    def configuration(*a):
        updated = current(*a)
        updated['ica'].update(yaml.safe_load(tuning.read_text()).get('ica',{}))
        return updated
    monkeypatch.setattr(review.triage,'subject_config',configuration)
    monkeypatch.setattr(publisher,'load_stage',lambda:stage)
    # Numerical eligibility has separate real-optimizer tests above. This test
    # intentionally isolates rollback after a fully installed candidate cache.
    monkeypatch.setattr(publisher,'require_candidate',lambda *a:None)
    def interrupted_render(argv):
        assert file_hash(model_path) == metadata['model_sha256']
        notes = yaml.safe_load(tuning.read_text()); notes['review_notes'] = {'concurrent':'preserved'}
        tuning.write_text(yaml.safe_dump(notes))
        manifest.write_text(json.dumps({'derivative_lineage':{'stages':{
            'ica_manual':{'status':'partial-new-render'},'other_stage':{'note':'preserved'}}}}))
        raise RuntimeError('injected rendering interruption')
    monkeypatch.setattr(stage,'main',interrupted_render)
    with pytest.raises(RuntimeError,match='injected rendering interruption'):
        publisher.publish(args)
    assert before == {path:file_hash(path) for path in protected}
    notes = yaml.safe_load(tuning.read_text())
    if restart:
        assert notes.get('ica') == old_ica_tuning
    else:
        assert notes['ica']['max_iter'] == cfg['ica']['max_iter']
    assert notes['review_notes']['concurrent'] == 'preserved'
    restored = json.loads(manifest.read_text())['derivative_lineage']['stages']
    assert restored['ica_manual'] == old_stage and restored['other_stage']['note'] == 'preserved'
    archive = dirs['reports']/'history/ica_review/failed-render'
    assert (archive/'originals'/model_path.relative_to(root)).is_file()
    assert not (args.candidate_dir/'publication.json').exists()
    # Retrying the same unpublished transaction must reuse an intact archive,
    # preserve unrelated notes/lineage, and still reject a changed figure stamp.
    args.resume_archive=True
    with pytest.raises(RuntimeError,match='injected rendering interruption'):
        publisher.publish(args)
    assert before == {path:file_hash(path) for path in protected}
    import os
    image=next(dirs['figures'].glob('*.png'));stat=image.stat()
    os.utime(image,ns=(stat.st_atime_ns,stat.st_mtime_ns+1_000_000_000))
    with pytest.raises(RuntimeError,match='Original evidence changed since archive'):
        publisher.publish(args)


def test_reconciled_cache_cannot_silently_revert_to_original_initialization(tmp_path,monkeypatch):
    """Changed cache policy must stop, preserving explicitly reconciled weights."""
    import json
    import review_final_ica as review
    from test_final_ica_review import saved_review
    _,cfg,source,model_path,metadata_path,_,_=saved_review(tmp_path,monkeypatch)
    metadata=json.loads(metadata_path.read_text());metadata['optimizer_restart']={'kind':'fixture-restart'}
    metadata_path.write_text(json.dumps(metadata))
    stage=review.load_stage()
    monkeypatch.setattr(stage.ICA,'fit',lambda *a,**kw:pytest.fail('Must not silently refit a restart'))
    settings=deepcopy(cfg['ica']);settings['reuse_prefit']=False
    with pytest.raises(RuntimeError,match='cannot be silently refitted'):
        stage.fit_or_load_ica(input_paths=[source],ica_cfg=settings,prefit_path=model_path,
            metadata_path=metadata_path,raw_root=cfg['resolved_paths']['data_root'],
            overwrite=True,require_fit_data=False)
