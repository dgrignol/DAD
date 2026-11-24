"""Prepare fixed, ROI-specific multivariate noise normalization for Stage 07.

Usage: prepare_noise(handoff, output_dir, cfg) from 07_run_matlab_drsa.py.
Input: accepted native Stage-04 epochs and immutable Stage-06 tensors/sidecars.
Output: noise_transforms.mat and a JSON provenance/QC sidecar in the selected
analysis namespace. No input or legacy manifest is modified. One pooled noise
model is shared across all conditions, identities, times and stimulus models.

The named v1 estimator matches the approved sandbox: negative-time epochs are
resampled to 120 Hz using edge padding, optionally Gaussian-smoothed to match
the target MATLAB analysis, and -0.4..-0.2 s is retained. Diagonal noise units
precede covariance shrinkage; an effective relative rank cutoff of 1e-6 avoids
amplifying near-null SSS directions. Shrinkage is chosen by acquisition-run
holdout likelihood, never by task RDM/model fit. Pearson after symmetric
whitening defines a fixed coordinate convention, not a rotation-invariant angle.
"""
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat
from scipy.ndimage import convolve1d
from continuous_cleaning import digest, file_hash, file_stamp
from drsa_analysis_variant import source_config

ALPHAS = np.array([0., 1e-6, 1e-5, 1e-4, .001, .01, .05, .1, .2, .5, 1.])
RANK_RTOL = 1e-6


def gaussian_smooth(data, width):
    """Match MATLAB smoothdata Gaussian weights and truncated edge support."""
    if width <= 1:
        return data
    offsets = np.arange(width) - (width - 1) / 2
    kernel = np.exp(-.5 * (offsets / (width / 5)) ** 2)
    numerator = convolve1d(data, kernel, axis=-1, mode='constant', cval=0.)
    denominator = convolve1d(np.ones(data.shape[-1]), kernel, mode='constant', cval=0.)
    return numerator / denominator


def fit_noise(base, rank_rtol=RANK_RTOL):
    """Estimate unit-invariant covariance and retained standardized support."""
    b = base.transpose(0, 2, 1).reshape(-1, base.shape[1]).copy()
    b -= b.mean(axis=0)
    sigma = np.sqrt(np.mean(b*b, axis=0))
    if not np.all(np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError('Noise SD must be positive and finite for every retained channel.')
    z = b / sigma
    cov = z.T @ z / len(z)
    eigenvalues, vectors = np.linalg.eigh(cov)
    keep = eigenvalues > eigenvalues[-1] * rank_rtol
    return sigma, eigenvalues[keep], vectors[:, keep], cov


def fit_transform(base, runs):
    """Select shrinkage on held-out runs, then fit one symmetric transform."""
    folds = []; heldout = []
    if len(np.unique(runs)) < 2:
        raise ValueError('MNN requires at least two acquisition runs for noise validation.')
    for run in np.unique(runs):
        sigma, eig, u, _ = fit_noise(base[runs != run])
        z = base[runs == run].transpose(0, 2, 1).reshape(-1, len(sigma)) / sigma
        z -= z.mean(axis=0)
        c = z.T @ z / len(z)
        pc = u.T @ c @ u
        losses = [float(np.mean(np.log(ev) + np.diag(pc)/ev))
                  for alpha in ALPHAS for ev in [(1-alpha)*eig + alpha*eig.mean()]]
        folds.append({'run': int(run), 'rank': len(eig), 'losses': losses})
        heldout.append((eig, pc, c))
    alpha = float(ALPHAS[np.argmin(np.mean([f['losses'] for f in folds], axis=0))])
    # These holdouts also select alpha: QC is not an unbiased nested-CV score.
    for fold, (fold_eig, pc, c) in zip(folds, heldout):
        fold_ev = (1-alpha)*fold_eig + alpha*fold_eig.mean()
        test = pc / np.sqrt(fold_ev[:, None]*fold_ev[None, :])
        correlation = test / np.sqrt(np.diag(test)[:, None]*np.diag(test)[None, :])
        diagonal_correlation = c / np.sqrt(np.diag(c)[:, None]*np.diag(c)[None, :])
        np.fill_diagonal(correlation, 0); np.fill_diagonal(diagonal_correlation, 0)
        fold['variance_median'] = float(np.median(np.diag(test)))
        fold['correlation_rms_retained_basis'] = float(np.sqrt(np.mean(correlation**2)))
        fold['correlation_rms_sensor_basis'] = float(np.sqrt(np.mean(diagonal_correlation**2)))
    sigma, eig, u, cov = fit_noise(base)
    ev = (1-alpha)*eig + alpha*eig.mean()
    white = (u / np.sqrt(ev)) @ u.T
    return {'white': white, 'sigma': sigma, 'rank': len(eig), 'alpha': alpha,
            'cov': cov, 'eigenvalues': eig}, folds


def prepare_noise(handoff, output_dir, cfg):
    """Validate and cache transforms for every ROI and temporal profile used."""
    import mne
    source = source_config(cfg)
    src = Path(source['resolved_paths']['derivatives_analysis_dir'])
    subject = cfg['layout']['subject_bids_id']; session = cfg['layout']['session_bids_id']
    epoch_path = src / 'meg' / f'{subject}_{session}_task-Attention_run-all_epo.fif'
    group_file = Path(handoff['project_root']) / 'code/analysis/megin_sensor_groups.json'
    groups = json.loads(group_file.read_text())
    # Common available sensor axes ensure every condition uses the same geometry.
    labels = list(handoff['channel_labels']); good = np.ones(len(labels), dtype=bool)
    for condition in handoff['condition_outputs']:
        tensor = loadmat(condition['neural_mat'], variable_names=['neural'])['neural']
        if tensor.shape[1] != len(labels):
            raise ValueError('Neural channel count disagrees with handoff.')
        bad = np.all(np.isnan(tensor), axis=(0, 2))
        if not np.all(np.isfinite(tensor[:, ~bad, :])):
            raise ValueError('Partially missing neural channels cannot enter MNN.')
        good &= ~bad
    labels = [label for label, keep in zip(labels, good) if keep]
    profiles = {}
    if 'real' in handoff['selected_analyses']:
        profiles['real'] = 0
    if 'alt' in handoff['selected_analyses']:
        settings = cfg.get('matlab_drsa', {})
        profiles['alt'] = int(settings.get('smooth_window_samples', 15)) if settings.get('smooth_neural', True) else 0
    if not profiles:
        return {'mode': 'mnn', 'transforms_mat': '', 'profiles': {}}
    settings = {'version': 'sensor_noise_mnn_v1', 'profiles': profiles,
                'baseline_sec': [-.4, -.2], 'resample_pad': 'edge',
                'rank_rtol': RANK_RTOL, 'alphas': ALPHAS.tolist(), 'labels': labels,
                'physical_units': {label: ('T' if label.endswith('1') else 'T/m') for label in labels},
                'groups_sha256': file_hash(group_file), 'helper_sha256': file_hash(__file__),
                'epochs': [file_stamp(p) for p in sorted(epoch_path.parent.glob(epoch_path.stem+'*.fif'))],
                'exports': [file_hash(c['sidecar_json']) for c in handoff['condition_outputs']]}
    contract = digest(settings)
    output_dir = Path(output_dir); mat_path = output_dir/'noise_transforms.mat'; json_path = output_dir/'noise_transforms.json'
    if json_path.exists() or mat_path.exists():
        if not (json_path.exists() and mat_path.exists()):
            raise RuntimeError('Incomplete noise artifacts; use a fresh analysis label.')
        cached = json.loads(json_path.read_text())
        if cached['contract_sha256'] != contract or cached['mat_sha256'] != file_hash(mat_path):
            raise RuntimeError('Noise inputs/settings/artifact changed; use a fresh analysis label.')
        return {'mode': 'mnn', 'transforms_mat': str(mat_path), 'contract_sha256': contract, 'profiles': profiles}
    epochs = mne.read_epochs(epoch_path, preload=False, proj=False, verbose='error')
    if epochs.metadata is None or 'acquisition_run' not in epochs.metadata:
        raise ValueError('Accepted epochs lack acquisition-run metadata.')
    runs = epochs.metadata['acquisition_run'].to_numpy(dtype=int)
    # Public MNE reads preserve channel calibration and split-FIF semantics.
    print(f'MNN {subject}: reading {len(epochs)} retained baseline epochs, {len(labels)} channels', flush=True)
    native = epochs.get_data(picks=labels, tmin=epochs.tmin, tmax=-1/epochs.info['sfreq'])
    baseline = mne.filter.resample(native, up=120, down=epochs.info['sfreq'], npad='auto', pad='edge', axis=-1, verbose=False)
    del native
    times = epochs.tmin + np.arange(baseline.shape[-1])/120
    mask = (times >= -.4-1e-9) & (times <= -.2+1e-9)
    if mask.sum() != 25:
        raise ValueError('Unexpected noise sample support; expected 25 at 120 Hz.')
    transforms = {}; qc = {}
    for profile, width in profiles.items():
        data = gaussian_smooth(baseline, width)[..., mask]
        transforms[profile] = {}; qc[profile] = {}
        profile_groups = {'wholebrain': labels} if profile == 'real' else groups
        for group, members in profile_groups.items():
            indices = [i for i, label in enumerate(labels) if label in members]
            if not indices:
                raise ValueError(f'No retained channels for {group}.')
            result, folds = fit_transform(data[:, indices], runs)
            result['labels'] = np.array([labels[i] for i in indices], dtype=object)
            transforms[profile][group] = result
            qc[profile][group] = {'rank': result['rank'], 'alpha': result['alpha'],
                                  'n_channels': len(indices), 'folds': folds}
            print(f'MNN {subject} {profile}/{group}: rank {result["rank"]}/{len(indices)}, alpha {result["alpha"]:g}', flush=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    savemat(mat_path, {'transforms': transforms}, do_compression=True)
    json_path.write_text(json.dumps({'contract_sha256': contract, 'settings': settings,
                         'n_epochs': len(epochs), 'software': {'mne': mne.__version__, 'numpy': np.__version__},
                         'mat_sha256': file_hash(mat_path), 'qc': qc}, indent=2)+'\n')
    return {'mode': 'mnn', 'transforms_mat': str(mat_path), 'contract_sha256': contract, 'profiles': profiles}
