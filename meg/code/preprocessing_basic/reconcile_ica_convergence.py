#!/usr/bin/env python3
"""Evaluate an isolated longer FastICA fit without overwriting a production model.

Usage (same MNE version as the saved model, bounded numerical threads)
-------------------------------------------------------------------
    python reconcile_ica_convergence.py --subject 19 --max-iter 10000 \
        --label NEW-LABEL --threads 4 --checkpoint-dir /tmp/NEW-checkpoint

If shared storage disappears after fitting, retain that explicitly UNVALIDATED
local checkpoint. After reconnecting, use --resume-checkpoint /tmp/NEW-checkpoint
with a NEW output label and the same subject/cap/runtime. Full source/model/
policy validation and numerical checks run again; no ICA refitting occurs.

Inputs: current unreviewed final model/report, unchanged accepted 4A/4B inputs,
protected 1-40 Hz/200 Hz preparation and original final ICA configuration.
Only the numerical iteration cap changes. Seed 97, FastICA contrast/algorithm,
n_components=None, native masks, channels and preparation remain unchanged.
Use the cluster numerical stack: MNE 1.10.2, NumPy 1.26.4, SciPy 1.13.1 and
scikit-learn 1.5.2. The original job log records MNE/NumPy/sklearn; SciPy follows
the pinned cluster requirements. Platform/BLAS details are recorded separately.

Outputs: a NEW derivatives/ica_convergence/sub-XX/LABEL candidate directory with
model/metadata, original-versus-candidate comparison, fit warnings/runtime and
fixed-point/PCA/reconstruction checks. No source, production model, tuning,
component decision, cleaned output or epoch is changed. Promotion is separate.

This participant-specific dependency cannot parallelize its own input fitting;
BLAS threads are bounded by --threads (and any SLURM allocation). Preparation
prints source progress; fitting prints its cap and numerical checks on completion.
Successful numerical convergence establishes no scientific component acceptance.
"""
from __future__ import annotations
import argparse
from contextlib import nullcontext
from copy import deepcopy
import json
import os
import platform
from pathlib import Path
import time
import warnings

import matplotlib.pyplot as plt
import mne
import numpy as np
import scipy
import sklearn
from scipy.optimize import linear_sum_assignment
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_info, threadpool_limits

from continuous_cleaning import file_hash
from ica_review_plots import continuous_display_copy
from preprocessing_stage_utils import guarded_path, utc_now, write_json
from review_final_ica import PARAMS, load_stage, validate_review
from ica_optimizer_checkpoint import iteration_cap, load_checkpoint, save_checkpoint


# -----------------------------------------------------------------------------
# Exact optimizer checks, not an artifact classifier or acceptance shortcut
# -----------------------------------------------------------------------------

def fixed_point_residual(model, fit_data):
    """Measure one parallel/logcosh FastICA update in the fitted PCA basis.

    Uses version-guarded MNE whitening and sklearn's original update primitives
    ONLY for diagnostics. It does not mutate the saved matrix or change the
    fitting algorithm. Residual is sklearn's sign-invariant stopping measure.
    """
    if sklearn.__version__ not in ('1.5.2', '1.7.2'):
        raise RuntimeError('Fixed-point diagnostic requires audited sklearn 1.5.2/1.7.2 update primitives.')
    from sklearn.decomposition._fastica import _logcosh, _sym_decorrelation
    if model.fit_params.get('algorithm') != 'parallel' or model.fit_params.get('fun') != 'logcosh':
        raise RuntimeError('Fixed-point check is defined only for the original parallel/logcosh fit.')
    values = model._pre_whiten(fit_data.get_data(picks=model.ch_names))
    values -= model.pca_mean_[:, None]
    n = model.n_components_
    norms = np.sqrt(model.pca_explained_variance_[:n])
    white = (model.pca_components_[:n] @ values) / norms[:, None]
    del values
    matrix = model.unmixing_matrix_ * norms[None, :]
    gx, derivative = _logcosh(matrix @ white, model.fit_params.get('fun_args') or {})
    updated = _sym_decorrelation(gx @ white.T / white.shape[1] - derivative[:, None] * matrix)
    residuals = np.abs(np.abs(np.einsum('ij,ij->i', updated, matrix)) - 1.)
    return {'max': float(residuals.max()), 'per_component': residuals.tolist(),
        'tol': float(model.fit_params.get('tol', 1e-4))}


def model_comparison(original, candidate, fit_data, output):
    """Quantify changed maps/activations without transferring exclusion indices."""
    display = continuous_display_copy(fit_data)
    a = original.get_sources(display).get_data()
    b = candidate.get_sources(display).get_data()
    a -= a.mean(axis=1, keepdims=True); b -= b.mean(axis=1, keepdims=True)
    correlation = (a @ b.T) / (np.linalg.norm(a, axis=1)[:, None] * np.linalg.norm(b, axis=1)[None, :])
    del a, b, display
    maps_a, maps_b = original.get_components(), candidate.get_components()
    cosine = (maps_a.T @ maps_b) / (np.linalg.norm(maps_a, axis=0)[:, None] * np.linalg.norm(maps_b, axis=0)[None, :])
    rows, columns = linear_sum_assignment(np.abs(correlation), maximize=True)
    matched = [{'original_IC': int(i), 'candidate_IC': int(j), 'activation_abs_correlation': float(abs(correlation[i,j])),
        'map_abs_cosine': float(abs(cosine[i,j]))} for i,j in zip(rows, columns, strict=True)]
    figure, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, values, title in zip(axes, (abs(correlation), abs(cosine)),
            ('Fit activation |correlation|', 'Mixing map |cosine|'), strict=True):
        artist = ax.imshow(values, origin='lower', vmin=0, vmax=1, cmap='viridis')
        ax.set_title(title); ax.set_xlabel('Candidate IC'); ax.set_ylabel('Original IC')
        figure.colorbar(artist, ax=ax)
    figure.suptitle('Model comparison only; no indices/exclusions transferred')
    figure.tight_layout(); figure.savefig(output / 'model-comparison.png', dpi=130); plt.close(figure)
    return {'matches_for_comparison_only': matched,
        'minimum_activation_match': min(row['activation_abs_correlation'] for row in matched),
        'median_activation_match': float(np.median([row['activation_abs_correlation'] for row in matched]))}


# -----------------------------------------------------------------------------
# Isolated candidate fit, exact source/rank checks and immutable result evidence
# -----------------------------------------------------------------------------

def main():
    """Fit/evaluate one candidate; never promote or accept it implicitly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(PARAMS / 'pipeline_default.yaml'))
    parser.add_argument('--site-config')
    parser.add_argument('--subject-config-dir', default=str(PARAMS))
    parser.add_argument('--subject', type=int, default=19)
    parser.add_argument('--max-iter', type=int, default=10000)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--label', required=True)
    parser.add_argument('--read-cache-certificate', type=Path,
        help='Optional certified byte-identical local reads; original validation remains required.')
    parser.add_argument('--checkpoint-dir',type=Path,
        help='NEW local/scratch directory; preserve unvalidated weights before final shared-source checks.')
    parser.add_argument('--resume-checkpoint',type=Path,
        help='Resume exact matching unvalidated weights after ordinary source/model/policy validation.')
    args = parser.parse_args()
    if args.threads < 1 or args.threads > int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1)):
        parser.error('Threads must fit the available CPU allocation.')
    if args.max_iter <= 1000 or not args.label.replace('-', '').replace('_', '').isalnum():
        parser.error('Use a new simple label and an iteration cap above 1000.')
    if args.checkpoint_dir and args.checkpoint_dir.exists():
        parser.error('Checkpoint output must be a NEW directory; preserve earlier computation.')
    certificate = None
    if args.read_cache_certificate:
        from verified_read_cache import verified_fif_reads
        certificate = json.loads(args.read_cache_certificate.read_text())
    with threadpool_limits(limits=args.threads), (verified_fif_reads(certificate) if certificate else nullcontext()):
        cfg, payload, records, inputs, original_path, original_metadata = validate_review(args, args.subject)
        if cfg['ica'].get('reviewed_model_sha256') or payload.get('manual_exclude') is not None or cfg.get('stage4_review',{}).get('03_C_ICA',{}).get('status')=='accepted':
            raise RuntimeError('Preserve existing final-model decisions; this automated candidate is for an unreviewed model only.')
        if mne.__version__ != payload['fingerprint']['fit']['mne_version']:
            raise RuntimeError('Use the exact MNE version of the original fit; no runtime fingerprint bypass.')
        if (mne.__version__, np.__version__, scipy.__version__, sklearn.__version__) != ('1.10.2','1.26.4','1.13.1','1.5.2'):
            raise RuntimeError('Use the pinned cluster numerical stack for an iteration-cap-only comparison.')
        original = mne.preprocessing.read_ica(original_path, verbose='ERROR')
        if original.method != 'fastica' or original.exclude or cfg['ica']['random_state'] != 97:
            raise RuntimeError('Unexpected algorithm, saved exclusions or seed; reconcile explicitly.')
        if original.n_iter_ < iteration_cap(original) or args.max_iter <= iteration_cap(original):
            raise RuntimeError('This repair requires a cap-limited original and a strictly larger cap.')
        before = file_hash(original_path)
        output = guarded_path(Path(cfg['resolved_paths']['derivatives_root']) /
            f"ica_convergence/sub-{args.subject:02}/{args.label}", cfg['resolved_paths']['data_root'])
        output.mkdir(parents=True, exist_ok=False)
        shutil_inputs = {original_path.name: before, original_metadata.name: file_hash(original_metadata)}
        write_json(output / 'baseline.json', {'generated_utc': utc_now(), 'model_files': shutil_inputs,
            'report_sha256': file_hash(Path(cfg['resolved_paths']['derivatives_root']) /
                f"reports/sub-{args.subject:02}/ses-01/03_C_ICA_review_report.json"),
            'fingerprint': payload['fingerprint'], 'model_settings': original.fit_params,
            'verified_read_cache': certificate,
            'authorization': '2026-11-10 user requested autonomous investigation/repair of P19 nonconvergence'},
            raw_root=cfg['resolved_paths']['data_root'])
        candidate_cfg = deepcopy(cfg['ica']); candidate_cfg['max_iter'] = args.max_iter
        core_runtime = {'mne':mne.__version__,'numpy':np.__version__,
                        'scipy':scipy.__version__,'sklearn':sklearn.__version__}
        checkpoint_context = {'subject':args.subject,'original_model_sha256':before,
            'original_fingerprint':payload['fingerprint'],'candidate_settings':candidate_cfg,
            'core_runtime':core_runtime}
        stage = load_stage()
        fit_data, summaries = stage.prepare_ica_fit_data(inputs, candidate_cfg)
        if fit_data.n_times != payload['fit_execution']['n_samples_fitted'] or [r['fit_samples'] for r in summaries] != [r['fit_samples'] for r in payload['fit_input_summaries']]:
            raise RuntimeError('Protected preparation does not reproduce exact fitted sample accounting.')
        print(f'ORIGINAL fixed-point check: {fit_data.n_times} samples', flush=True)
        old_residual = fixed_point_residual(original, fit_data)
        print('ORIGINAL residual:', old_residual['max'], 'tol:', old_residual['tol'], flush=True)
        parameters = {k:v for k,v in original.fit_params.items() if k != 'max_iter'}
        if args.resume_checkpoint:
            candidate,checkpoint = load_checkpoint(args.resume_checkpoint,checkpoint_context)
            elapsed = checkpoint['elapsed_fit_sec']; warning_records = checkpoint['warnings']
            if iteration_cap(candidate) != args.max_iter or {k:v for k,v in candidate.fit_params.items() if k!='max_iter'} != parameters:
                raise RuntimeError('Checkpoint optimizer parameters differ from the requested repair.')
            print('Reusing exact UNVALIDATED checkpoint; no ICA fitting. All numerical/source checks follow.',flush=True)
        else:
            candidate = mne.preprocessing.ICA(method='fastica', n_components=candidate_cfg['n_components'],
                random_state=97, max_iter=args.max_iter, fit_params=parameters, verbose='WARNING')
            print(f'CANDIDATE fitting seed=97, n_components=None, max_iter={args.max_iter}, threads={args.threads}', flush=True)
            started = time.monotonic()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                candidate.fit(fit_data, reject_by_annotation=True, verbose='WARNING')
            elapsed = time.monotonic() - started
            warning_records = [{'category':w.category.__name__,'message':str(w.message),
                'convergence_warning':issubclass(w.category,ConvergenceWarning)} for w in caught]
        convergence_warnings = [w['message'] for w in warning_records
            if w.get('convergence_warning',w['category']=='ConvergenceWarning')]
        print(f'CANDIDATE fit complete: {candidate.n_iter_} iterations, {elapsed:.1f} s, warnings={convergence_warnings}', flush=True)
        # Preserve computation BEFORE any network-dependent final validation.
        # Checkpoint metadata remains unvalidated even if convergence looks good.
        if args.checkpoint_dir:
            save_checkpoint(args.checkpoint_dir,candidate,checkpoint_context,warning_records,
                            elapsed,cfg['resolved_paths']['data_root'])
        new_residual = fixed_point_residual(candidate, fit_data)
        rank_same = candidate.n_components_ == original.n_components_
        pca_same = rank_same and np.allclose(candidate.pca_explained_variance_[:candidate.n_components_],
            original.pca_explained_variance_[:original.n_components_], rtol=1e-6, atol=1e-9)
        finite = all(np.isfinite(getattr(candidate, name)).all() for name in
            ('unmixing_matrix_', 'mixing_matrix_', 'pca_components_', 'pre_whitener_'))
        sample = fit_data.copy().crop(0., min(10., fit_data.times[-1])).load_data()
        before_samples = sample.get_data()
        candidate.apply(sample, exclude=[], verbose='ERROR')
        error = float(np.linalg.norm(sample.get_data()-before_samples) / np.linalg.norm(before_samples))
        del sample, before_samples
        fingerprint = stage.build_ica_fingerprint(input_paths=inputs, ica_cfg=candidate_cfg)
        if fingerprint['inputs'] != payload['fingerprint']['inputs'] or file_hash(original_path) != before:
            raise RuntimeError('Original inputs/model changed during candidate evaluation.')
        path = output / original_path.name
        candidate.save(path, overwrite=False, verbose='ERROR')
        metadata = {'generated_utc':utc_now(), 'fingerprint':fingerprint, 'model_path':str(path),
            'model_sha256':file_hash(path), 'n_components_fitted':int(candidate.n_components_),
            'n_samples_fitted':int(candidate.n_samples_), 'n_iterations':int(candidate.n_iter_),
            'fit_input_summaries':summaries, 'pca_explained_variance':candidate.pca_explained_variance_.tolist()}
        write_json(output / original_metadata.name, metadata, raw_root=cfg['resolved_paths']['data_root'])
        comparison = model_comparison(original, candidate, fit_data, output)
        numeric_ok = (not convergence_warnings and candidate.n_iter_ < args.max_iter and
            new_residual['max'] < new_residual['tol'] and rank_same and pca_same and finite and error < 1e-10 and
            candidate.n_samples_ == fit_data.n_times and candidate.ch_names == original.ch_names)
        result = {'generated_utc':utc_now(), 'subject':args.subject, 'numeric_candidate_ok':bool(numeric_ok),
            'scientific_acceptance':False, 'production_model_changed':False,
            'original_model_sha256':before, 'candidate_model_sha256':metadata['model_sha256'],
            'changed_setting_only':{'max_iter':{'original':iteration_cap(original), 'candidate':args.max_iter}},
            'runtime':{'mne':mne.__version__,'sklearn':sklearn.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
                'python':platform.python_version(),'platform':platform.platform(),'threadpools':threadpool_info()},
            'elapsed_fit_sec':elapsed,'original_residual':old_residual,'candidate_residual':new_residual,
            'verified_read_cache':certificate,
            'candidate_iterations':int(candidate.n_iter_), 'convergence_warnings':convergence_warnings,
            'warnings':warning_records,
            'resumed_checkpoint':str(args.resume_checkpoint) if args.resume_checkpoint else None,
            'rank_unchanged':rank_same,'leading_pca_spectrum_unchanged':pca_same,'weights_finite':finite,
            'no_removal_reconstruction_relative_error':error,'comparison':comparison}
        write_json(output / 'result.json', result, raw_root=cfg['resolved_paths']['data_root'], overwrite=False)
        print('CANDIDATE RESULT:', output / 'result.json', 'numerically acceptable:', numeric_ok, flush=True)


if __name__ == '__main__':
    main()
