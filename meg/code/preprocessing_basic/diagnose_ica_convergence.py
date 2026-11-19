#!/usr/bin/env python3
"""Inspect FastICA cycling/contrast without refitting or changing a final model.

Usage in the pinned cluster numerical environment::

    python diagnose_ica_convergence.py --subject 19 --steps 8 --threads 1 \
        --output-dir /local/scratch/NEW-P19-diagnosis \
        --read-cache-certificate /local/scratch/read-cache/certificate.json

Inputs: current unreviewed final model, unchanged accepted 4A/4B and Stage 2/3
evidence, exact protected fit copies. Validate these with the ordinary review
guard before AND after diagnosis. Cache use retains original source checks.

Data flow: original stored PCA whitening -> original unmixing directions ->
several unchanged sklearn parallel/logcosh updates. Compare one- and two-step
row alignments, including permutation-invariant alignment. Measure whitening,
non-Gaussianity descriptors and squared-activation dependence. These distinguish
slow convergence, near-cycling and label exchange; they classify no IC as neural
or artifact and select no exclusion, rank, tolerance, seed or mask change.

Outputs: a NEW local directory with an exact original-model copy, baseline and
pending numerical JSON, then final result.json only after repeat source/model
validation. Pending results remain explicitly unvalidated if a volume vanishes.
No production model, tuning, receipt, native source, annotation or epoch changes.
One participant is an ordered dependency; --threads bounds numerical libraries.
Print each source and numerical update. No new ICA optimizer fit is performed.
"""
from argparse import ArgumentParser
from contextlib import nullcontext
import json
import os
from pathlib import Path
import shutil

import mne
import numpy as np
import scipy
from scipy import signal
from scipy.optimize import linear_sum_assignment
import sklearn
from threadpoolctl import threadpool_limits

from continuous_cleaning import file_hash
from preprocessing_stage_utils import guarded_path, utc_now, write_json
from review_final_ica import PARAMS, load_stage, validate_review


# -----------------------------------------------------------------------------
# Exact FastICA basis/update and sign/permutation-aware recurrence measurements
# -----------------------------------------------------------------------------

def fitted_basis(model, fit_data):
    """Reconstruct its stored PCA-whitened samples without mutating either input."""
    if model.fit_params.get('algorithm')!='parallel' or model.fit_params.get('fun')!='logcosh':
        raise RuntimeError('This diagnostic is for the original parallel/logcosh FastICA only.')
    values=model._pre_whiten(fit_data.get_data(picks=model.ch_names))
    values-=model.pca_mean_[:,None]
    n=model.n_components_
    norms=np.sqrt(model.pca_explained_variance_[:n])
    white=(model.pca_components_[:n]@values)/norms[:,None]
    matrix=model.unmixing_matrix_*norms[None,:]
    return white,matrix


def alignment(new,old):
    """Separate sign changes, label permutations and actual direction changes."""
    overlap=new@old.T
    residual=np.abs(np.abs(np.diag(overlap))-1.)
    i,j=linear_sum_assignment(np.abs(overlap),maximize=True)
    return {'rowwise_max_residual':float(residual.max()),
        'rowwise_median_residual':float(np.median(residual)),
        'permutation_matched_max_residual':float(np.max(np.abs(np.abs(overlap[i,j])-1.))),
        'per_component_rowwise_residual':residual.tolist()}


def piece_spectra(values,sfreq,piece_lengths,nperseg=1024):
    """Estimate diagnostic spectra using windows wholly inside retained pieces."""
    windows=[];offset=0
    for size in piece_lengths:
        windows.extend(range(offset,offset+size-nperseg+1,nperseg//2))
        offset+=size
    if offset!=values.shape[1]:raise RuntimeError('Spectral piece map differs from fit samples.')
    if not windows:return {'status':'no_complete_within_piece_windows'}
    data=np.stack([values[:,start:start+nperseg] for start in windows],axis=1)
    frequencies,psd=signal.welch(data,fs=sfreq,nperseg=nperseg,noverlap=0,
                               detrend='constant',axis=-1)
    mean=psd.mean(axis=1);band=(frequencies>=1)&(frequencies<=40)
    result={'status':'within_good_fit_pieces_only','windows':len(windows),
        'window_samples':nperseg,'sfreq':sfreq,'frequencies':frequencies[band].tolist(),
        'power':mean[:,band].tolist(),
        'peak_hz':[float(frequencies[band][np.argmax(row[band])]) for row in mean]}
    if len(values)==2:
        _,cross=signal.csd(data[0],data[1],fs=sfreq,nperseg=nperseg,noverlap=0,
                           detrend='constant',axis=-1)
        cross=cross.mean(axis=0)
        coherence=np.abs(cross)**2/(mean[0]*mean[1])
        peak=np.flatnonzero(band)[np.argmax(np.sqrt(mean[0,band]*mean[1,band]))]
        result['pair']={'joint_peak_hz':float(frequencies[peak]),
            'coherence_at_joint_peak':float(coherence[peak]),
            'phase_degrees_at_joint_peak':float(np.angle(cross[peak],deg=True))}
    return result


def diagnose_arrays(white,matrix,parameters,steps=8,*,sfreq=None,piece_lengths=None):
    """Describe recurrence and contrast strength; establish no scientific acceptance."""
    if sklearn.__version__ not in ('1.5.2','1.7.2'):
        raise RuntimeError('Use audited sklearn update primitives (1.5.2/1.7.2).')
    from sklearn.decomposition._fastica import _logcosh,_sym_decorrelation
    n_samples=white.shape[1]
    covariance=white@white.T/(n_samples-1)
    values=matrix@white
    moments={power:np.mean(values**power,axis=1) for power in (2,3,4)}
    gx,derivative=_logcosh(values.copy(),parameters.get('fun_args') or {})
    beta=np.mean(values*gx,axis=1)-derivative
    proposal=gx@white.T/n_samples-derivative[:,None]*matrix
    singular=scipy.linalg.svdvals(proposal)
    # Express the update in the current IC coordinates. Mixed curvature signs
    # can make its polar factor a reflection, producing near two-step recurrence.
    curvature=proposal@matrix.T
    eigenvalues,eigenvectors=scipy.linalg.eigh((curvature+curvature.T)/2)
    positive=np.flatnonzero(eigenvalues>0)
    virtual=eigenvectors[:,positive].T@values
    directions=[]
    for column,index in enumerate(positive):
        z=virtual[column];variance=np.mean(z**2)
        directions.append({'direction':int(column),'curvature_eigenvalue':float(eigenvalues[index]),
            'current_IC_weights':eigenvectors[:,index].tolist(),
            'skewness':float(np.mean(z**3)/variance**1.5),
            'excess_kurtosis':float(np.mean(z**4)/variance**2-3)})
    spectra=piece_spectra(virtual,sfreq,piece_lengths) if len(positive) and sfreq and piece_lengths else None
    # Power correlations are descriptive nonlinear dependence, not proof of
    # source independence or a neural/artifact label.
    power=values**2;power-=power.mean(axis=1,keepdims=True)
    norms=np.linalg.norm(power,axis=1)
    power_correlation=(power@power.T)/(norms[:,None]*norms[None,:])
    descriptors=[{'IC':int(i),'skewness':float(moments[3][i]/moments[2][i]**1.5),
        'excess_kurtosis':float(moments[4][i]/moments[2][i]**2-3),
        'logcosh_curvature_beta':float(beta[i]),'update_row_norm':float(np.linalg.norm(proposal[i]))}
        for i in range(len(matrix))]
    del values,gx,power
    history=[matrix.copy()];recurrence=[]
    for step in range(1,steps+1):
        gx,derivative=_logcosh(history[-1]@white,parameters.get('fun_args') or {})
        update=_sym_decorrelation(gx@white.T/n_samples-derivative[:,None]*history[-1])
        row={'step':step,'one_step':alignment(update,history[-1])}
        if step>=2:row['two_step']=alignment(update,history[-2])
        recurrence.append(row);history.append(update)
        print(f'FastICA diagnostic update {step}/{steps}: one-step residual {row["one_step"]["rowwise_max_residual"]:.6g}',flush=True)
    return {'scientific_acceptance':False,'component_decisions':None,
        'n_components':len(matrix),'fit_samples':n_samples,
        'whitening_covariance_max_error':float(np.max(np.abs(covariance-np.eye(len(matrix))))),
        'unmixing_orthogonality_max_error':float(np.max(np.abs(matrix@matrix.T-np.eye(len(matrix))))),
        'update_singular_values':singular.tolist(),'component_descriptors':descriptors,
        'curvature_matrix_in_IC_coordinates':curvature.tolist(),
        'curvature_symmetry_relative_error':float(np.linalg.norm(curvature-curvature.T)/np.linalg.norm(curvature)),
        'curvature_eigenvalues':eigenvalues.tolist(),'positive_curvature_directions':directions,
        'positive_curvature_direction_spectra':spectra,
        'squared_activation_correlation':power_correlation.tolist(),
        'one_step_relative_directions':(history[1]@history[0].T).tolist(),
        'two_step_relative_directions':(history[2]@history[0].T).tolist(),
        'recurrence':recurrence,'tol':float(parameters.get('tol',1e-4)),
        'interpretation_limit':'Descriptors and short recurrence do not establish Gaussianity, source type, usefulness of removal, or final-model acceptance.'}


# -----------------------------------------------------------------------------
# Source-bound diagnosis, local preservation and repeat validation before result
# -----------------------------------------------------------------------------

def main():
    """Inspect one exact current model; preserve pending computation locally."""
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--config',default=str(PARAMS/'pipeline_default.yaml'))
    parser.add_argument('--site-config');parser.add_argument('--subject-config-dir',default=str(PARAMS))
    parser.add_argument('--subject',type=int,default=19);parser.add_argument('--expected-model')
    parser.add_argument('--read-cache-certificate',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--steps',type=int,default=8);parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--save-numeric-cache',action='store_true',
        help='Save exact local whitened arrays for subsequent diagnostic probes; never an operative model.')
    args=parser.parse_args()
    allocated=int(os.environ.get('SLURM_CPUS_PER_TASK',os.cpu_count() or 1))
    if args.steps<2 or not 1<=args.threads<=allocated or args.output_dir.exists():
        parser.error('Use at least two steps, positive threads and a NEW output directory.')
    certificate=None
    if args.read_cache_certificate:
        from verified_read_cache import verified_fif_reads
        certificate=json.loads(args.read_cache_certificate.read_text())
    with threadpool_limits(args.threads),(verified_fif_reads(certificate) if certificate else nullcontext()):
        cfg,payload,records,inputs,model_path,_=validate_review(args,args.subject)
        if (mne.__version__,np.__version__,scipy.__version__,sklearn.__version__)!=('1.10.2','1.26.4','1.13.1','1.5.2'):
            raise RuntimeError('Use the pinned original cluster numerical stack.')
        raw_root=cfg['resolved_paths']['data_root']
        output=guarded_path(args.output_dir,raw_root);output.mkdir(parents=True,exist_ok=False)
        before=file_hash(model_path)
        shutil.copy2(model_path,output/'original-model-ica.fif')
        if file_hash(output/'original-model-ica.fif')!=before:raise RuntimeError('Original-model copy differs.')
        baseline={'generated_utc':utc_now(),'subject':args.subject,'model_sha256':before,
            'fingerprint':payload['fingerprint'],'scientific_acceptance':False}
        write_json(output/'baseline.json',baseline,raw_root=raw_root,overwrite=False)
        model=mne.preprocessing.read_ica(model_path,verbose='ERROR')
        fit,summaries=load_stage().prepare_ica_fit_data(inputs,cfg['ica'])
        if fit.n_times!=model.n_samples_ or [r['fit_samples'] for r in summaries]!=[r['fit_samples'] for r in payload['fit_input_summaries']]:
            raise RuntimeError('Protected fit preparation differs from fitted sample accounting.')
        white,matrix=fitted_basis(model,fit)
        cache=None
        if args.save_numeric_cache:
            np.save(output/'whitened-fit.npy',white)
            np.save(output/'original-whitened-unmixing.npy',matrix)
            cache={name:file_hash(output/name) for name in ('whitened-fit.npy','original-whitened-unmixing.npy')}
        lengths=[row['fit_samples'] for summary in summaries for row in summary['intervals'] if row['status']=='included']
        result={**baseline,**diagnose_arrays(white,matrix,model.fit_params,args.steps,
            sfreq=float(fit.info['sfreq']),piece_lengths=lengths),
            'pca_explained_variance':model.pca_explained_variance_.tolist(),
            'local_diagnostic_cache':cache,'fit_piece_lengths':lengths,
            'fit_run_samples':[row['fit_samples'] for row in summaries],
            'source_validation':'pending_final_source_and_model_check'}
        write_json(output/'pending-result.json',result,raw_root=raw_root,overwrite=False)
        # The completed calculations survive a disconnected share; they never
        # become validated just because the arithmetic completed successfully.
        _,after,_,_,current,_=validate_review(args,args.subject)
        if file_hash(current)!=before or after['fingerprint']!=payload['fingerprint']:
            raise RuntimeError('Source/model context changed during diagnosis.')
        result['source_validation']='passed_before_and_after_diagnosis'
        write_json(output/'result.json',result,raw_root=raw_root,overwrite=False)
        print('DIAGNOSTIC RESULT',output/'result.json',flush=True)


if __name__=='__main__':
    main()
