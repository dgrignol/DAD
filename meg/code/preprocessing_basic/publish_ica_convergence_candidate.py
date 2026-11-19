#!/usr/bin/env python3
"""Publish a validated numerical ICA repair and regenerate its unreviewed evidence.

Usage after an explicitly authorized numerical repair candidate:
    python publish_ica_convergence_candidate.py --subject 19 \
        --candidate-dir /derivatives/ica_convergence/sub-19/LABEL \
        --archive-label p19-convergence-20261110-v1

Requires the same MNE version and exact accepted sources/preparation as the
candidate. Rejects nonconverged/nonfinite candidates, changed PCA rank/samples,
existing final exclusions/acceptance, or source/config/model changes. Accepts
either a cap-only candidate or a fully documented deterministic restart of the
saved seed-97 solution with identical configuration/PCA/whitening. Archives the old model, metadata, reports, review
figures, tuning and manifest before any replacement. Numerical checks are not
scientific component or cleaned-output acceptance.

Cap-only publication records the participant-only cap in tuning; restart
publication leaves tuning unchanged and records phase-specific history. Install the
new exact model and calls the ordinary 4C REVIEW path. A guard forbids any ICA
refitting while rendering. All properties/native/task views are regenerated;
no old component indices are transferred. Native trial identities must match.
On failure restore the archived C evidence/model and the original cap, retaining
other manual/timing sections. No annotations, SSS, ICA application or epochs are
written. This single-participant dependency is ordered; bound BLAS threads with
OMP_NUM_THREADS/OPENBLAS_NUM_THREADS or the existing scheduler allocation.
"""
from copy import deepcopy
from contextlib import nullcontext
import argparse
import json
import os
from pathlib import Path
import shutil
import traceback

import mne
import numpy as np
import scipy
import sklearn

from continuous_cleaning import bind_report, digest, file_hash, file_stamp, require_review
from preprocessing_stage_utils import stage_output_dirs, update_tuning_file, utc_now, write_json
from render_stage03_mne_reports import render_stage_report
from review_final_ica import PARAMS, load_stage, settle_review_outputs, validate_review
import annotation_triage as triage
from ica_optimizer_checkpoint import iteration_cap


# -----------------------------------------------------------------------------
# Reject unproven candidates and preserve the currently unreviewed model
# -----------------------------------------------------------------------------

def require_candidate(result, metadata, candidate, original, payload):
    """Recheck numeric gates and the exact model hash/rank before publication."""
    if not result.get('numeric_candidate_ok') or result.get('convergence_warnings'):
        raise RuntimeError('Candidate has not passed the numerical convergence checks.')
    residual = result['candidate_residual']
    if residual['max'] >= residual['tol'] or result['no_removal_reconstruction_relative_error'] >= 1e-10:
        raise RuntimeError('Candidate fixed-point/reconstruction check failed.')
    if candidate.n_iter_ != result['candidate_iterations'] or candidate.n_iter_ >= iteration_cap(candidate):
        raise RuntimeError('Candidate iteration accounting does not establish convergence.')
    restart = result.get('repair_kind') == 'deterministic_curvature_restart_v1'
    if candidate.method != 'fastica' or (not restart and iteration_cap(candidate) != metadata['fingerprint']['fit']['max_iter']):
        raise RuntimeError('Candidate algorithm/cap differs from its evaluated fingerprint.')
    if {k:v for k,v in candidate.fit_params.items() if k!='max_iter'} != {k:v for k,v in original.fit_params.items() if k!='max_iter'}:
        raise RuntimeError('Candidate changed optimizer parameters beyond the iteration cap.')
    if candidate.exclude or candidate.n_components_ != original.n_components_ or candidate.ch_names != original.ch_names:
        raise RuntimeError('Candidate changed component count/channels or contains exclusions.')
    if candidate.n_samples_ != payload['fit_execution']['n_samples_fitted']:
        raise RuntimeError('Candidate changed fitted sample coverage.')
    if not np.allclose(candidate.pca_explained_variance_[:candidate.n_components_],
                       original.pca_explained_variance_[:original.n_components_], rtol=1e-6, atol=1e-9):
        raise RuntimeError('Candidate PCA spectrum changed unexpectedly.')
    if not all(np.isfinite(getattr(candidate, k)).all() for k in ('unmixing_matrix_', 'mixing_matrix_', 'pca_components_')):
        raise RuntimeError('Candidate contains nonfinite weights.')
    old_fit, new_fit = deepcopy(payload['fingerprint']['fit']), deepcopy(metadata['fingerprint']['fit'])
    if restart:
        # The objective/configuration is unchanged; numerical initialization is
        # recorded explicitly as optimizer history, never as a fresh RNG fit.
        history = result.get('optimizer_restart', {})
        if (metadata['fingerprint'] != payload['fingerprint'] or
            metadata.get('optimizer_restart') != history or
            history.get('kind') != 'deterministic_curvature_restart_v1' or
            history.get('fresh_random_draw') is not False or
            history.get('random_state_argument') != 97 or old_fit.get('random_state') != 97 or
            history.get('requested_n_components') is not None or candidate.n_components is not None or
            history.get('restart_phase_iterations') != candidate.n_iter_ or
            history.get('original_phase_iterations') != original.n_iter_ or
            history.get('iteration_cap') != iteration_cap(original) or
            iteration_cap(candidate) != iteration_cap(original) or
            result.get('source_validation') != 'passed_before_and_after_candidate_validation' or
            result.get('native_auxiliary_samples_equal') is not True or
            not history.get('initializer_sha256') or not history.get('diagnosis_sha256')):
            raise RuntimeError('Restart provenance/configuration/source validation is incomplete or changed.')
        for name in ('pca_components_','pca_mean_','pca_explained_variance_','pre_whitener_'):
            if not np.array_equal(getattr(candidate,name),getattr(original,name)):
                raise RuntimeError('Restart changed original PCA/whitening.')
        return
    old_fit.pop('max_iter'); new_fit.pop('max_iter')
    if old_fit != new_fit or metadata['fingerprint']['inputs'] != payload['fingerprint']['inputs']:
        raise RuntimeError('Only the iteration cap may change; source/preparation/runtime must match.')


def main():
    """Archive, publish and fully revalidate one authorized unreviewed repair."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default=str(PARAMS/'pipeline_default.yaml'))
    p.add_argument('--site-config')
    p.add_argument('--subject-config-dir',default=str(PARAMS))
    p.add_argument('--subject',type=int,default=19)
    p.add_argument('--candidate-dir',type=Path,required=True)
    p.add_argument('--archive-label',required=True)
    p.add_argument('--resume-archive',action='store_true',
        help='Reuse a complete same-candidate archive after failed publication; all original file stamps must still match.')
    p.add_argument('--read-cache-certificate', type=Path,
        help='Optional byte-verified local I/O cache; original source checks remain required.')
    args=p.parse_args()
    if not args.archive_label.replace('-','').replace('_','').isalnum():
        p.error('Use a new simple archive label.')
    certificate = None
    if args.read_cache_certificate:
        from verified_read_cache import verified_fif_reads
        certificate = json.loads(args.read_cache_certificate.read_text())
    # The cache changes only where identical bytes are read. Existing validation
    # and report code continue receiving/binding the original accepted paths.
    with verified_fif_reads(certificate) if certificate else nullcontext():
        publish(args, certificate)


def publish(args, certificate=None):
    """Run ordered validation/archive/publication under optional verified reads."""
    cfg,payload,records,inputs,model_path,metadata_path=validate_review(args,args.subject)
    dirs=stage_output_dirs(cfg); root=Path(cfg['resolved_paths']['derivatives_root'])
    raw_root=Path(cfg['resolved_paths']['data_root'])
    if cfg['ica'].get('reviewed_model_sha256') or payload.get('manual_exclude') is not None or cfg.get('stage4_review',{}).get('03_C_ICA',{}).get('status')=='accepted':
        raise RuntimeError('Existing final-model choices/acceptance must be preserved; automated publication refused.')
    if not args.candidate_dir.resolve().is_relative_to((root/f'ica_convergence/sub-{args.subject:02}').resolve()):
        raise RuntimeError('Candidate must be in this participant\'s isolated reconciliation directory.')
    result_path=args.candidate_dir/'result.json'
    result=json.loads(result_path.read_text())
    meta=json.loads((args.candidate_dir/metadata_path.name).read_text())
    candidate_path=args.candidate_dir/model_path.name
    candidate_sha=file_hash(candidate_path)
    if result['subject']!=args.subject or result['original_model_sha256']!=file_hash(model_path):
        raise RuntimeError('Candidate baseline no longer matches the active model.')
    if candidate_sha!=result['candidate_model_sha256'] or candidate_sha!=meta['model_sha256']:
        raise RuntimeError('Candidate model contents changed since evaluation.')
    if mne.__version__!=meta['fingerprint']['fit']['mne_version']:
        raise RuntimeError('Use the exact candidate MNE runtime.')
    actual_runtime = {'mne':mne.__version__, 'numpy':np.__version__,
                      'scipy':scipy.__version__, 'sklearn':sklearn.__version__}
    if any(result['runtime'].get(key) != value for key,value in actual_runtime.items()):
        raise RuntimeError('Use the numerical package versions recorded by the evaluated candidate.')
    candidate=mne.preprocessing.read_ica(candidate_path,verbose='ERROR')
    original=mne.preprocessing.read_ica(model_path,verbose='ERROR')
    require_candidate(result,meta,candidate,original,payload)
    stage=load_stage(); candidate_cfg=deepcopy(cfg['ica'])
    restart = result.get('repair_kind') == 'deterministic_curvature_restart_v1'
    if restart:
        history = result['optimizer_restart']
        if (history['seed_97_ancestry_original_model_sha256'] != result['original_model_sha256'] or
            history['candidate_model_sha256'] != candidate_sha or
            history['candidate_directory'] != str(args.candidate_dir)):
            raise RuntimeError('Restart ancestry/model/directory binding changed.')
        for name, checksum in result['supporting_file_hashes'].items():
            path = args.candidate_dir/name
            if Path(name).name != name or file_hash(path) != checksum:
                raise RuntimeError('Restart supporting evidence changed.')
    else:
        candidate_cfg['max_iter']=result['changed_setting_only']['max_iter']['candidate']
    if stage.build_ica_fingerprint(input_paths=inputs,ica_cfg=candidate_cfg)!=meta['fingerprint']:
        raise RuntimeError('Current source/fit fingerprint does not reproduce the candidate.')
    tuning=dirs['reports']/'03_tuning_choices.yaml'
    source_manifest=dirs['reports']/'source_manifest.json'
    old_manifest=json.loads(source_manifest.read_text()) if source_manifest.exists() else None
    native_manifest=dirs['figures']/'03_C_ICA_review_native_diagnostics.json'
    old_trials=json.loads(native_manifest.read_text())['retained_trial_metadata']
    archive=dirs['reports']/'history/ica_review'/args.archive_label
    backup=[]
    paths = [model_path,metadata_path,tuning,*[dirs['reports']/f'03_C_ICA_review_{x}'
            for x in ('report.json','report.html','components.html')],source_manifest,
            *sorted(dirs['figures'].glob('03_C_ICA_review*'))]
    resume = getattr(args,'resume_archive',False)
    if resume:
        saved = json.loads((archive/'archive.json').read_text())
        if (saved['original_model_sha256'] != file_hash(model_path) or
            saved['candidate_model_sha256'] != candidate_sha or
            saved['candidate_result'] != str(result_path) or
            (args.candidate_dir/'publication.json').exists()):
            raise RuntimeError('Existing archive is not an unpublished same-candidate transaction.')
        expected_paths = [str(p.relative_to(root)) for p in paths if p.is_file()]
        if saved['files'] != expected_paths:
            raise RuntimeError('Existing archive inventory differs from current original evidence.')
        for path in paths:
            if not path.is_file():continue
            relative=path.relative_to(root)
            target=archive/'originals'/relative
            if path==source_manifest:
                # Rollback restores only this stage, preserving unrelated
                # concurrently added lineage. Its JSON mtime may therefore differ.
                active=json.loads(path.read_text()).get('derivative_lineage',{}).get('stages',{}).get('ica_manual')
                saved_stage=json.loads(target.read_text()).get('derivative_lineage',{}).get('stages',{}).get('ica_manual')
                if active!=saved_stage:raise RuntimeError('Original ICA manifest lineage changed since archive.')
            elif path==tuning:
                # Source/model/config/review validation above uses CURRENT
                # tuning. Preserve unrelated review notes added after backup.
                pass
            elif file_stamp(path)!=file_stamp(target):
                raise RuntimeError(f'Original evidence changed since archive: {path.name}')
            backup.append((path,target))
        print(f'Reusing complete original archive after exact inventory/stamp checks: {len(backup)} files',flush=True)
    else:
        archive.mkdir(parents=True,exist_ok=False)
        for path in paths:
            if path.is_file():
                relative=path.relative_to(root)
                target=archive/'originals'/relative
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
                backup.append((path,target))
                if len(backup)%25==0:print(f'Archived {len(backup)} original P{args.subject} evidence files',flush=True)
        write_json(archive/'archive.json',{'generated_utc':utc_now(),'original_model_sha256':file_hash(model_path),
            'candidate_model_sha256':candidate_sha,'candidate_result':str(result_path),
            'files':[str(path.relative_to(root)) for path,_ in backup]},raw_root=raw_root,overwrite=False)
    old_cap=cfg['ica']['max_iter']
    dirty = {model_path,metadata_path}
    try:
        print('Publishing numerical candidate; final component choices remain unreviewed.',flush=True)
        if not restart:
            update_tuning_file(path=tuning,raw_root=raw_root,subject_id=f'{args.subject:02}',
                section='ica',values={'max_iter':candidate_cfg['max_iter']})
        shutil.copy2(candidate_path,model_path)
        meta.update(model_path=str(model_path),convergence_reconciliation_result=str(result_path))
        write_json(metadata_path,meta,raw_root=raw_root)
        # Explicitly finish timestamps on newly written DERIVATIVE outputs.
        # In-place SMB replacement can expose zero/stale mtimes after close.
        # Original sources/archives are never touched; exact bytes and all
        # source/model validation still apply, followed by stable positive stamps.
        for path in (model_path,metadata_path):
            os.utime(path,None)
            print(f'Installed output stamp: {file_stamp(path)}',flush=True)
        if file_hash(model_path)!=candidate_sha:
            raise RuntimeError('Installed candidate bytes changed during output publication.')
        settle_review_outputs([model_path,metadata_path])
        new_cfg=triage.subject_config(args,args.subject)
        expected=stage.build_ica_fingerprint(input_paths=inputs,ica_cfg=new_cfg['ica'])
        if not stage.cache_is_compatible(metadata_path,expected):
            raise RuntimeError('Published cache is not compatible; report rendering cannot proceed.')
        # Rendering must reuse this evaluated model. Any hidden refit is an error.
        original_fit=stage.ICA.fit
        def no_refit(*a,**kw):
            raise RuntimeError('Unexpected ICA fitting during candidate publication.')
        stage.ICA.fit=no_refit
        try:
            command=['--config',args.config,'--subject-config',str(Path(args.subject_config_dir)/f'sub-{args.subject:02}.yaml'),
                '--tuning-config',str(tuning),'--source-manifest',str(source_manifest),
                '--processing-pass','review','--overwrite','--no-display','--print-summary']
            if args.site_config:command+=['--site-config',args.site_config]
            # Only this call may write report/figure/manifest outputs. Earlier
            # failures must not recopy hundreds of untouched shared figures.
            dirty.update(path for path,_ in backup if path!=tuning)
            if stage.main(command)!=0:raise RuntimeError('Candidate review-report regeneration failed.')
        finally:
            stage.ICA.fit=original_fit
        cfg2,new_payload,_,_,new_model,_=validate_review(args,args.subject)
        if file_hash(new_model)!=candidate_sha or not new_payload['component_review_coverage']['complete'] or new_payload['fit_execution']['iteration_limit_reached']:
            raise RuntimeError('Published model/diagnostics did not pass final verification.')
        require_review(cfg2,'03_A_MaxFilt',dirs['reports']);require_review(cfg2,'03_B_ArtefactAnnotation',dirs['reports'])
        if digest(json.loads(native_manifest.read_text())['retained_trial_metadata'])!=digest(old_trials):
            raise RuntimeError('Diagnostic trial selection changed during an optimizer-only repair.')
        report=dirs['reports']/'03_C_ICA_review_report.json'
        new_payload['convergence_reconciliation']={'original_model_sha256':result['original_model_sha256'],
            'candidate_model_sha256':candidate_sha,'archive':str(archive),'result':str(result_path),
            'repair_kind':result.get('repair_kind','iteration_cap_only'),
            'numeric_checks_passed':True,'component_decisions_transferred':False,'scientific_acceptance':False}
        paths=[root/row['relative_path'] for row in new_payload['evidence_files']]+[result_path]
        bind_report(new_payload,paths,{'fingerprint':new_payload['fingerprint'],'model_sha256':candidate_sha,
            'convergence_result_sha256':file_hash(result_path)},dirs['reports'],root)
        write_json(report,new_payload,raw_root=raw_root)
        render_stage_report(derivatives_root=root,raw_root=raw_root,subject_id=f'sub-{args.subject:02}',
            processing_pass='review',stage_name='ica_manual',overwrite=True,primary_report=True)
        write_json(args.candidate_dir/'publication.json',{'generated_utc':utc_now(),'published':True,
            'old_model_sha256':result['original_model_sha256'],'new_model_sha256':candidate_sha,
            'archive':str(archive),'iteration_cap':candidate_cfg['max_iter'],'scientific_acceptance':False,
            'repair_kind':result.get('repair_kind','iteration_cap_only'),
            'optimizer_restart':result.get('optimizer_restart'), 'tuning_changed':not restart,
            'component_review_required':True,'diagnostic_trials_unchanged':True,
            'verified_read_cache':certificate},raw_root=raw_root,overwrite=False)
        print('PUBLISHED numerical repair; all component decisions remain explicit.',flush=True)
    except BaseException:
        traceback.print_exc()
        print('Publication failed; restoring the archived C model/evidence.',flush=True)
        for path,target in backup:
            if path not in dirty or path in (tuning,source_manifest):continue
            shutil.copy2(target,path)
        if not restart:
            update_tuning_file(path=tuning,raw_root=raw_root,subject_id=f'{args.subject:02}',section='ica',values={'max_iter':old_cap})
        if old_manifest is not None and source_manifest in dirty:
            current=json.loads(source_manifest.read_text())
            old_stage=old_manifest.get('derivative_lineage',{}).get('stages',{}).get('ica_manual')
            if old_stage is None:
                current.get('derivative_lineage',{}).get('stages',{}).pop('ica_manual',None)
            else:
                current['derivative_lineage']['stages']['ica_manual']=old_stage
            write_json(source_manifest,current,raw_root=raw_root)
        if file_hash(model_path)!=result['original_model_sha256']:
            raise RuntimeError('Rollback did not restore the original model contents.')
        raise


if __name__=='__main__':
    main()
