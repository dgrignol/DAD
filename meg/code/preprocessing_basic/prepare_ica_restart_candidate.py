#!/usr/bin/env python3
"""Validate a completed deterministic FastICA restart against rebuilt fit data.

Example (the original pinned numerical environment, one participant dependency)::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python prepare_ica_restart_candidate.py \
      --subject 19 --probe-dir meg/sandbox/RESTART --diagnosis-dir meg/sandbox/DIAG \
      --label NEW-LABEL --threads 1 --read-cache-certificate /scratch/certificate.json

Read the unchanged accepted model/sources and an isolated restart checkpoint.
Rebuild the protected full fit preparation, reproduce the deterministic
initializer, check the serialized candidate's actual fixed-point residual and
native no-removal reconstruction, then repeat ordinary source validation.
All 69 P19 directions selected originally by n_components=None are retained;
FastICA parallel/logcosh, tolerance, seed-97 ancestry and cap remain explicit.
This is a restart of saved optimizer state, NOT a fresh random MNE fit.

Write a NEW shared candidate directory containing exact model bytes, metadata,
initializer/output matrices, diagnosis/probe evidence, source-code snapshots,
model comparison and numerical result. Never replace the operative model,
write choices/acceptance, apply cleaning or create production epochs. Publication
is separate. --threads bounds numerical libraries; print source/check progress.
"""
import argparse
from contextlib import nullcontext
from copy import deepcopy
import json
import os
from pathlib import Path
import platform
import shutil

import mne
import numpy as np
import scipy
import sklearn
from threadpoolctl import threadpool_info, threadpool_limits

from continuous_cleaning import file_hash
from diagnose_ica_convergence import fitted_basis
from ica_optimizer_checkpoint import iteration_cap, load_checkpoint
from preprocessing_stage_utils import guarded_path, utc_now, write_json
from probe_ica_curvature_restart import curvature_initializer
from reconcile_ica_convergence import fixed_point_residual, model_comparison
from review_final_ica import PARAMS, load_stage, validate_review


# -----------------------------------------------------------------------------
# Full preparation and serialized-weight checks precede candidate publication
# -----------------------------------------------------------------------------

def portable_summaries(rows):
    """Compare exact interval maps across the documented cluster/Mac roots.

    Only the input mount prefix is presentation-specific. Every coordinate,
    channel, filter-support field, interval and sample count remains exact.
    Ordinary source/annotation validation separately binds the resolved files.
    """
    result = deepcopy(rows)
    for row in result:
        row['input_path'] = Path(row['input_path']).name
    return result


def main():
    """Create immutable, fully source-validated evidence for one restart."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(PARAMS/'pipeline_default.yaml'))
    parser.add_argument('--site-config')
    parser.add_argument('--subject-config-dir', default=str(PARAMS))
    parser.add_argument('--subject', type=int, default=19)
    parser.add_argument('--probe-dir', type=Path, required=True)
    parser.add_argument('--diagnosis-dir', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--read-cache-certificate', type=Path)
    args = parser.parse_args()
    if not 1 <= args.threads <= int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 1)):
        parser.error('Threads exceed the available CPU allocation.')
    if not args.label.replace('-', '').replace('_', '').isalnum():
        parser.error('Use a new simple output label.')
    certificate = None
    if args.read_cache_certificate:
        from verified_read_cache import verified_fif_reads
        certificate = json.loads(args.read_cache_certificate.read_text())
    with threadpool_limits(args.threads), (verified_fif_reads(certificate) if certificate else nullcontext()):
        cfg, payload, _, inputs, original_path, metadata_path = validate_review(args, args.subject)
        if payload.get('manual_exclude') is not None or cfg['ica'].get('reviewed_model_sha256'):
            raise RuntimeError('Preserve existing final-model decisions; restart adoption refused.')
        runtime = {'mne': mne.__version__, 'numpy': np.__version__,
                   'scipy': scipy.__version__, 'sklearn': sklearn.__version__}
        if runtime != {'mne':'1.10.2','numpy':'1.26.4','scipy':'1.13.1','sklearn':'1.5.2'}:
            raise RuntimeError('Use the original pinned numerical stack.')
        baseline = json.loads((args.probe_dir/'baseline.json').read_text())
        probe = json.loads((args.probe_dir/'result.json').read_text())
        diagnosis = json.loads((args.diagnosis_dir/'result.json').read_text())
        old_sha = file_hash(original_path)
        if (baseline['subject'] != args.subject or baseline['original_model_sha256'] != old_sha
                or baseline['original_fingerprint'] != payload['fingerprint']
                or baseline['diagnosis_sha256'] != file_hash(args.diagnosis_dir/'result.json')
                or diagnosis['model_sha256'] != old_sha
                or diagnosis['source_validation'] != 'passed_before_and_after_diagnosis'
                or not probe['numeric_checks_passed']
                or probe['source_validation'] != 'passed_before_and_after_probe'
                or probe['fresh_random_draw'] or probe['random_state_argument'] != 97):
            raise RuntimeError('Probe/diagnosis do not validate the exact current context.')
        candidate, checkpoint = load_checkpoint(args.probe_dir/'unvalidated-checkpoint', baseline)
        if file_hash(args.probe_dir/'unvalidated-checkpoint/checkpoint_ica.fif') != probe['candidate_model_sha256']:
            raise RuntimeError('Evaluated probe model changed.')
        original = mne.preprocessing.read_ica(original_path, verbose='ERROR')
        for name in ('pca_components_', 'pca_mean_', 'pca_explained_variance_', 'pre_whitener_'):
            if not np.array_equal(getattr(candidate, name), getattr(original, name)):
                raise RuntimeError(f'Restart changed the original PCA/whitening: {name}')
        if (candidate.n_components is not None or candidate.n_components_ != original.n_components_
                or candidate.fit_params != original.fit_params or candidate.method != 'fastica'
                or cfg['ica']['random_state'] != 97 or cfg['ica']['n_components'] is not None):
            raise RuntimeError('Restart changed frozen dimension/optimizer configuration.')

        print('Rebuilding every protected fit sample; no ICA refitting.', flush=True)
        fit_data, summaries = load_stage().prepare_ica_fit_data(inputs, cfg['ica'])
        if portable_summaries(summaries) != portable_summaries(payload['fit_input_summaries']) or fit_data.n_times != candidate.n_samples_:
            raise RuntimeError('Rebuilt fit intervals/sample accounting changed.')
        white, original_matrix = fitted_basis(original, fit_data)
        initial, _ = curvature_initializer(white, original_matrix, original.fit_params)
        saved_initial = np.load(args.probe_dir/'initializer-whitened.npy')
        final_matrix = np.load(args.probe_dir/'final-whitened-unmixing.npy')
        if (file_hash(args.probe_dir/'initializer-whitened.npy') != probe['initializer_sha256']
                or not np.allclose(initial, saved_initial, rtol=0., atol=1e-12)
                or not np.allclose(candidate.unmixing_matrix_ * np.sqrt(candidate.pca_explained_variance_[:candidate.n_components_])[None,:],
                                   final_matrix, rtol=0., atol=1e-12)):
            raise RuntimeError('Initializer or actual serialized optimizer output cannot be reproduced.')
        del white, initial, original_matrix, saved_initial, final_matrix
        residual = fixed_point_residual(candidate, fit_data)
        original_residual = fixed_point_residual(original, fit_data)
        raw = mne.io.read_raw_fif(inputs[0], preload=False, verbose='ERROR').crop(0.,10.).load_data()
        before = raw.get_data(picks=original.ch_names)
        aux = mne.pick_types(raw.info, meg=False, stim=True, misc=True, exclude=[])
        aux_before = raw.get_data(picks=aux)
        candidate.apply(raw, exclude=[], verbose='ERROR')
        error = float(np.linalg.norm(raw.get_data(picks=original.ch_names)-before)/np.linalg.norm(before))
        auxiliary_equal = np.array_equal(raw.get_data(picks=aux), aux_before, equal_nan=True)
        ok = (not checkpoint['warnings'] and candidate.n_iter_ < iteration_cap(candidate)
              and residual['max'] < residual['tol'] and error < 1e-10 and auxiliary_equal
              and candidate.ch_names == original.ch_names and not candidate.exclude
              and all(np.isfinite(getattr(candidate,k)).all() for k in
                      ('unmixing_matrix_','mixing_matrix_','pca_components_','pre_whitener_')))
        if not ok:
            raise RuntimeError('Serialized restart failed numerical/source-preservation requirements.')

        # Keep the numerical history and script bytes alongside the exact model.
        root = Path(cfg['resolved_paths']['derivatives_root'])
        raw_root = cfg['resolved_paths']['data_root']
        output = guarded_path(root/f'ica_convergence/sub-{args.subject:02}/{args.label}', raw_root)
        output.mkdir(parents=True, exist_ok=False)
        pairs = [(original_path, 'original-model-ica.fif'),
                 (args.probe_dir/'unvalidated-checkpoint/checkpoint_ica.fif', original_path.name),
                 (args.probe_dir/'unvalidated-checkpoint/checkpoint.json', 'probe-checkpoint.json'),
                 (args.probe_dir/'result.json', 'probe-result.json'),
                 (args.diagnosis_dir/'result.json', 'diagnosis-result.json'),
                 (args.probe_dir/'initializer-whitened.npy', 'initializer-whitened.npy'),
                 (args.probe_dir/'final-whitened-unmixing.npy', 'final-whitened-unmixing.npy')]
        for name in ('probe_ica_curvature_restart.py','prepare_ica_restart_candidate.py',
                     'diagnose_ica_convergence.py','reconcile_ica_convergence.py'):
            pairs.append((Path(__file__).with_name(name), name))
        for source, name in pairs:
            shutil.copy2(source, output/name)
        supporting = {name:file_hash(output/name) for _,name in pairs}
        provenance = {'kind':'deterministic_curvature_restart_v1',
            'initialization':baseline['initialization'], 'fresh_random_draw':False,
            'seed_97_ancestry_original_model_sha256':old_sha, 'random_state_argument':97,
            'original_phase_iterations':int(original.n_iter_),
            'restart_phase_iterations':int(candidate.n_iter_), 'iteration_cap':iteration_cap(candidate),
            'requested_n_components':None, 'effective_n_components':int(candidate.n_components_),
            'n_samples':int(candidate.n_samples_), 'algorithm':'parallel','contrast':'logcosh',
            'tol':residual['tol'], 'initializer_sha256':supporting['initializer-whitened.npy'],
            'diagnosis_sha256':supporting['diagnosis-result.json'],
            'candidate_model_sha256':file_hash(output/original_path.name),
            'candidate_directory':str(output), 'scientific_acceptance':False}
        meta = deepcopy(json.loads(metadata_path.read_text()))
        meta.update(generated_utc=utc_now(), model_path=str(output/original_path.name),
                    model_sha256=provenance['candidate_model_sha256'], n_iterations=int(candidate.n_iter_),
                    optimizer_restart=provenance)
        write_json(output/metadata_path.name, meta, raw_root=raw_root, overwrite=False)
        comparison = model_comparison(original, candidate, fit_data, output)
        result = {'generated_utc':utc_now(), 'subject':args.subject,
            'repair_kind':'deterministic_curvature_restart_v1','optimizer_restart':provenance,
            'numeric_candidate_ok':True,'source_validation':'pending_final_source_model_check',
            'original_model_sha256':old_sha,'candidate_model_sha256':provenance['candidate_model_sha256'],
            'candidate_iterations':int(candidate.n_iter_), 'candidate_residual':residual,
            'original_residual':original_residual,'convergence_warnings':[],
            'no_removal_reconstruction_relative_error':error,'native_auxiliary_samples_equal':auxiliary_equal,
            'rank_unchanged':True,'pca_and_whitener_exactly_unchanged':True,
            'configuration_fingerprint_unchanged':True,'supporting_file_hashes':supporting,
            'runtime':{**runtime,'python':platform.python_version(),'platform':platform.platform(),
                       'threadpools':threadpool_info()},'verified_read_cache':certificate,
            'comparison':comparison,'scientific_acceptance':False,'production_model_changed':False,
            'component_decisions_transferred':False}
        write_json(output/'pending-result.json',result,raw_root=raw_root,overwrite=False)
        _, after, _, _, current, _ = validate_review(args,args.subject)
        if file_hash(current)!=old_sha or after['fingerprint']!=payload['fingerprint']:
            raise RuntimeError('Source/model context changed during candidate validation.')
        result['source_validation']='passed_before_and_after_candidate_validation'
        write_json(output/'result.json',result,raw_root=raw_root,overwrite=False)
        print(f'VALIDATED RESTART CANDIDATE {output}: {candidate.n_iter_} restart iterations; residual {residual["max"]}',flush=True)


if __name__ == '__main__':
    main()
