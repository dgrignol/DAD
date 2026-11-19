#!/usr/bin/env python3
"""Test a deterministic restart of a cycling FastICA optimizer, never publish.

Usage in the original pinned numerical environment::

    python probe_ica_curvature_restart.py --subject 19 --threads 1 \
        --diagnosis-dir meg/sandbox/VALIDATED-DIAGNOSIS-WITH-NUMERIC-CACHE \
        --output-dir meg/sandbox/NEW-RESTART-PROBE \
        --read-cache-certificate /local/scratch/read-cache/certificate.json

Inputs: an exact currently unreviewed model, ordinary accepted source/policy
bindings, and hash-verified whitened arrays from diagnose_ica_convergence.py.
The original MNE fit selected its PCA dimension with n_components=None. Keep
that entire basis, every fitted sample, parallel FastICA, logcosh and the original
stopping tolerance. Re-express the saved seed-97 optimizer directions in the
eigenbasis of its symmetric update-curvature matrix, then run STANDARD sklearn
FastICA with random_state=97 and this explicit deterministic initialization.
This is NOT a fresh random draw. Its origin, initializer bytes and phase-specific
iterations are recorded so seed/provenance cannot be misrepresented.

Outputs: a NEW LOCAL probe directory, immutable numerical progress checkpoints,
an explicitly UNVALIDATED MNE optimizer checkpoint, numeric checks and candidate
metadata only. The MNE checkpoint retains the original PCA/whitener/channel basis;
only the ICA matrix and its pseudoinverse are replaced with actual optimizer
output. Verify next-update residual, dimensions, finite weights and native
no-removal reconstruction, then repeat ordinary full source/model validation.
No production model, tuning, manual review, component exclusion, cleaned data or
epoch is changed. Successful numerical convergence is not scientific acceptance.
One participant is an ordered dependency; --threads bounds numerical libraries.
"""
from argparse import ArgumentParser
from contextlib import nullcontext
from copy import deepcopy
import inspect
import json
import os
from pathlib import Path
import time
import warnings

import mne
import numpy as np
import scipy
import sklearn
from sklearn.decomposition import FastICA
import sklearn.decomposition._fastica as optimizer
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

from continuous_cleaning import file_hash
from ica_optimizer_checkpoint import iteration_cap,save_checkpoint
from preprocessing_stage_utils import guarded_path,utc_now,write_json
from review_final_ica import PARAMS,validate_review


# -----------------------------------------------------------------------------
# Deterministic initialization; the subsequent optimizer is unchanged sklearn
# -----------------------------------------------------------------------------

def curvature_initializer(white,matrix,parameters):
    """Rotate the saved complete basis into its symmetric update eigenvectors."""
    gx,derivative=optimizer._logcosh(matrix@white,parameters.get('fun_args') or {})
    proposal=gx@white.T/white.shape[1]-derivative[:,None]*matrix
    curvature=proposal@matrix.T
    eigenvalues,vectors=scipy.linalg.eigh((curvature+curvature.T)/2)
    return vectors.T@matrix,eigenvalues


def run_restart(white,matrix,parameters,*,max_iter,seed=97,progress_dir=None):
    """Run real parallel/logcosh FastICA; observation hooks alter no update."""
    initial,eigenvalues=curvature_initializer(white,matrix,parameters)
    base=optimizer._sym_decorrelation
    previous=[None];history=[]
    def observe(values):
        updated=base(values)
        frame=inspect.currentframe().f_back
        try:
            step=frame.f_locals.get('ii') if frame.f_code.co_name=='_ica_par' else None
        finally:del frame
        if step is not None and previous[0] is not None:
            residual=float(np.max(np.abs(np.abs(np.einsum('ij,ij->i',updated,previous[0]))-1.)))
            history.append({'iteration':int(step+1),'residual':residual})
            if step%10==0:print(f'Curvature restart iteration {step+1}/{max_iter}: residual {residual:.6g}',flush=True)
            if progress_dir and (step+1)%50==0:
                np.save(Path(progress_dir)/f'iteration-{step+1:06}-whitened-weights.npy',updated)
        previous[0]=updated.copy()
        return updated
    params={k:v for k,v in parameters.items() if k not in ('max_iter','w_init')}
    estimator=FastICA(n_components=None,whiten=False,random_state=seed,max_iter=max_iter,
                      w_init=initial,**params)
    started=time.monotonic();optimizer._sym_decorrelation=observe
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always');estimator.fit(white.T)
    finally:optimizer._sym_decorrelation=base
    records=[{'category':w.category.__name__,'message':str(w.message),
              'convergence_warning':issubclass(w.category,ConvergenceWarning)} for w in caught]
    return estimator,initial,eigenvalues,history,records,time.monotonic()-started


def represented_model(original,whitened_matrix,n_iterations,max_iter,seed=97):
    """Pack actual optimizer output in the unchanged original MNE PCA basis."""
    model=original.copy()
    norms=np.sqrt(model.pca_explained_variance_[:model.n_components_])
    model.unmixing_matrix_=whitened_matrix/norms[None,:]
    model._update_mixing_matrix()
    model.n_iter_=int(n_iterations);model.max_iter=int(max_iter)
    model.fit_params=deepcopy(original.fit_params);model.fit_params['max_iter']=int(max_iter)
    model.random_state=seed;model.exclude=[]
    return model


# -----------------------------------------------------------------------------
# Current-source validation, local checkpoints and explicit unaccepted result
# -----------------------------------------------------------------------------

def main():
    """Probe one exact current model and preserve all actual numerical evidence."""
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--config',default=str(PARAMS/'pipeline_default.yaml'))
    parser.add_argument('--site-config');parser.add_argument('--subject-config-dir',default=str(PARAMS))
    parser.add_argument('--subject',type=int,default=19);parser.add_argument('--expected-model')
    parser.add_argument('--diagnosis-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--read-cache-certificate',type=Path)
    parser.add_argument('--threads',type=int,default=1)
    args=parser.parse_args()
    allocated=int(os.environ.get('SLURM_CPUS_PER_TASK',os.cpu_count() or 1))
    if not 1<=args.threads<=allocated or args.output_dir.exists():
        parser.error('Use available numerical threads and a NEW local output directory.')
    certificate=None
    if args.read_cache_certificate:
        from verified_read_cache import verified_fif_reads
        certificate=json.loads(args.read_cache_certificate.read_text())
    with threadpool_limits(args.threads),(verified_fif_reads(certificate) if certificate else nullcontext()):
        cfg,payload,records,inputs,model_path,_=validate_review(args,args.subject)
        if payload.get('manual_exclude') is not None or cfg['ica'].get('reviewed_model_sha256'):
            raise RuntimeError('Preserve existing final-model choices; automated restart probe refused.')
        if (mne.__version__,np.__version__,scipy.__version__,sklearn.__version__)!=('1.10.2','1.26.4','1.13.1','1.5.2'):
            raise RuntimeError('Use the original pinned numerical stack.')
        diagnosis=json.loads((args.diagnosis_dir/'result.json').read_text())
        before=file_hash(model_path)
        if (diagnosis['source_validation']!='passed_before_and_after_diagnosis' or
            diagnosis['model_sha256']!=before or diagnosis['fingerprint']!=payload['fingerprint']):
            raise RuntimeError('Diagnostic context differs from the exact current model/sources.')
        for name,checksum in diagnosis['local_diagnostic_cache'].items():
            if file_hash(args.diagnosis_dir/name)!=checksum:raise RuntimeError('Diagnostic numerical cache changed.')
        white=np.load(args.diagnosis_dir/'whitened-fit.npy',mmap_mode='r')
        matrix=np.load(args.diagnosis_dir/'original-whitened-unmixing.npy')
        original=mne.preprocessing.read_ica(model_path,verbose='ERROR')
        if original.method!='fastica' or cfg['ica']['random_state']!=97 or cfg['ica']['n_components'] is not None:
            raise RuntimeError('Unexpected frozen scientific ICA configuration.')
        if original.fit_params['algorithm']!='parallel' or original.fit_params['fun']!='logcosh':
            raise RuntimeError('This probe preserves the original parallel/logcosh algorithm.')
        expected=original.unmixing_matrix_*np.sqrt(original.pca_explained_variance_[:original.n_components_])[None,:]
        if white.shape!=(original.n_components_,original.n_samples_) or not np.array_equal(matrix,expected):
            raise RuntimeError('Cached ICA basis or fitted sample count differs.')
        root=cfg['resolved_paths']['data_root']
        output=guarded_path(args.output_dir,root);output.mkdir(parents=True,exist_ok=False)
        context={'subject':args.subject,'original_model_sha256':before,
            'original_fingerprint':payload['fingerprint'],'diagnosis_sha256':file_hash(args.diagnosis_dir/'result.json'),
            'initialization':'symmetric update-curvature eigenbasis of saved seed-97 solution',
            'fresh_random_draw':False,'random_state_argument':97,'requested_n_components':None}
        write_json(output/'baseline.json',context,raw_root=root,overwrite=False)
        estimator,initial,eigenvalues,history,caught,elapsed=run_restart(white,matrix,original.fit_params,
            max_iter=iteration_cap(original),seed=97,progress_dir=output)
        candidate=represented_model(original,estimator.components_,estimator.n_iter_,iteration_cap(original))
        np.save(output/'initializer-whitened.npy',initial)
        np.save(output/'final-whitened-unmixing.npy',estimator.components_)
        save_checkpoint(output/'unvalidated-checkpoint',candidate,context,caught,elapsed,root)
        gx,derivative=optimizer._logcosh(estimator.components_@white,original.fit_params.get('fun_args') or {})
        next_matrix=optimizer._sym_decorrelation(gx@white.T/white.shape[1]-derivative[:,None]*estimator.components_)
        residual=np.abs(np.abs(np.einsum('ij,ij->i',next_matrix,estimator.components_))-1.)
        raw=mne.io.read_raw_fif(inputs[0],preload=False,verbose='ERROR').crop(0.,10.).load_data()
        meg_before=raw.get_data(picks=original.ch_names)
        aux=mne.pick_types(raw.info,meg=False,stim=True,misc=True,exclude=[])
        aux_before=raw.get_data(picks=aux)
        candidate.apply(raw,exclude=[],verbose='ERROR')
        error=float(np.linalg.norm(raw.get_data(picks=original.ch_names)-meg_before)/np.linalg.norm(meg_before))
        auxiliary_equal=np.array_equal(raw.get_data(picks=aux),aux_before,equal_nan=True)
        valid=bool(estimator.n_iter_<iteration_cap(original) and not any(w['convergence_warning'] for w in caught)
            and residual.max()<estimator.tol and error<1e-10 and auxiliary_equal
            and np.isfinite(estimator.components_).all())
        result={**context,'generated_utc':utc_now(),'purpose':'isolated_optimizer_restart_probe',
            'numeric_checks_passed':valid,'scientific_acceptance':False,'production_model_changed':False,
            'method':'fastica','algorithm':'parallel','contrast':'logcosh','tol':estimator.tol,
            'restart_phase_iterations':int(estimator.n_iter_),'original_phase_iterations':int(original.n_iter_),
            'iteration_cap':iteration_cap(original),'elapsed_restart_sec':elapsed,'history':history,
            'curvature_eigenvalues':eigenvalues.tolist(),'next_update_residual':float(residual.max()),
            'per_component_next_update_residual':residual.tolist(),
            'no_removal_native_MEG_relative_error':error,'native_auxiliary_samples_equal':auxiliary_equal,
            'initializer_sha256':file_hash(output/'initializer-whitened.npy'),
            'candidate_model_sha256':file_hash(output/'unvalidated-checkpoint/checkpoint_ica.fif'),
            'source_validation':'pending_final_source_model_check','warnings':caught,
            'component_decisions_transferred':False}
        write_json(output/'pending-result.json',result,raw_root=root,overwrite=False)
        _,after,_,_,current,_=validate_review(args,args.subject)
        if file_hash(current)!=before or after['fingerprint']!=payload['fingerprint']:
            raise RuntimeError('Source/model context changed during restart probe.')
        result['source_validation']='passed_before_and_after_probe'
        write_json(output/'result.json',result,raw_root=root,overwrite=False)
        print('RESTART RESULT',valid,'iterations',estimator.n_iter_,'next residual',float(residual.max()),
              'native no-removal error',error,flush=True)


if __name__=='__main__':
    main()
