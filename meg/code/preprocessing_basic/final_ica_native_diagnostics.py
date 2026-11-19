#!/usr/bin/env python3
"""Default continuous AND task-locked diagnostics for an existing final ICA.

Usage from the report builder (no fitting/application or production epochs)
----------------------------------------------------------------------------
    images, manifest = generate_native_diagnostics(cfg, model, records, inputs,
        figure_dir, model_sha256)
    data = extract_native_diagnostics(cfg, model, records, inputs, components=[0])
    fig, meta = task_figure_from_data(cfg, data, component=0)

Read each native annotated acquisition once in bounded chunks, project the
unchanged final model, and derive all components' diagnostic epochs together.
Keep recording boundaries and absolute sample coordinates. Accepted Stage 2/3
eligibility/timing gates supply motion-onset events; configured BAD rejection,
detrend and baseline apply. Sensor amplitude thresholds do not apply to IC AU.
No Raw/Epochs FIF, model, tuning, mask or review receipt is saved.

Outputs are two PNGs/JSON sidecars per final IC plus a shared JSON manifest with
exact retained trial identities, drop reasons, acquisition windows and policies.
Continuous examples show one deterministic 20 s native window per acquisition,
anchored two seconds before its first ocular BAD (fixed fallback if absent).
They include BAD periods; task diagnostics ordinarily reject BAD-overlap trials.
Window means are subtracted ONLY from the displayed continuous traces; no
temporal filter/resampling or component removal is applied to native signals.

Resource controls belong to the participant process/SLURM task. Projection is
chunked (30 s default) and a run's full IC buffer is released after extraction;
all components reuse the same source read and task selection. Participants run
in disjoint processes; shared publication waits for their join. Reproduction
commands may request a single IC/run, using identical extraction/plot routines.
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import mne
import numpy as np
import pandas as pd

from ica_review_plots import attach_annotations, continuous_display_copy, context_selection
from preprocessing_stage_utils import guarded_path, write_json


# -----------------------------------------------------------------------------
# Native sample projection and accepted task anchors, shared across every IC
# -----------------------------------------------------------------------------

def extract_native_diagnostics(cfg, model, records, inputs, *, components=None,
                               include_bad=False, with_epochs=True, chunk_sec=30.):
    """Return transient IC epochs/windows without changing native files or choices.

    Components are final-model zero-based indices. Task values have shape
    (retained trials, requested ICs, native samples); acquisition windows retain
    source IDs, absolute sample bounds and all overlapping BAD descriptions.
    """
    from trial_preparation import load_prepared_decisions
    from visual_timing_utils import analysis_eligible, load_reviewed_timing
    components = list(range(model.n_components_)) if components is None else list(components)
    if not components or any(not 0 <= c < model.n_components_ for c in components):
        raise ValueError('Requested component is outside the exact final model.')
    rows, gate = [], None
    if with_epochs:
        root = Path(cfg['resolved_paths']['derivatives_root'])
        review = cfg['epoching'].get('timing_review_path') or root / f"visual_timing/s3-v002/reviews/{cfg['layout']['subject_bids_id']}_review.json"
        rows, gate = load_reviewed_timing(cfg, load_prepared_decisions(cfg), review)
    pieces, kept, dropped, windows, run_counts = [], [], {}, [], []
    times = None
    for record, path in zip(records, inputs, strict=True):
        print(f"Native diagnostics: {record.source_filename}; {len(components)} final ICs; no cleaning", flush=True)
        raw = mne.io.read_raw_fif(path, preload=False, verbose='ERROR')
        sfreq = raw.info['sfreq']
        projected = np.empty((len(components), raw.n_times))
        width = max(1, int(round(chunk_sec * sfreq)))
        for start in range(0, raw.n_times, width):
            stop = min(start + width, raw.n_times)
            print(f"  {record.source_filename}: samples {start}-{stop}/{raw.n_times}", flush=True)
            chunk = raw.copy().crop(start / sfreq, (stop - 1) / sfreq)
            sources = model.get_sources(continuous_display_copy(chunk))
            projected[:, start:stop] = sources.get_data(picks=components)
        del chunk, sources

        # A fixed, component-independent example per run, with native timing.
        selection, _ = context_selection(raw, record)
        start_sec = min(selection['window_start_recording_sec'], max(0., raw.times[-1] - 20.))
        left = int(round(start_sec * sfreq))
        right = min(raw.n_times, left + int(round(20. * sfreq)) + 1)
        start_native, stop_native = raw.first_samp + left, raw.first_samp + right
        bad = [{'description': str(a['description']), 'onset_acquisition_sec': float(a['onset']),
                'duration_sec': float(a['duration'])} for a in raw.annotations
               if a['description'].lower().startswith('bad') and a['onset'] < stop_native / sfreq
               and a['onset'] + a['duration'] >= start_native / sfreq]
        window = {**selection, 'window_start_sample_raw': start_native,
            'window_stop_sample_raw_exclusive': stop_native,
            'window_start_acquisition_sec': start_native / sfreq,
            'window_stop_acquisition_sec': stop_native / sfreq,
            'window_start_recording_sec': left / sfreq, 'bad_annotations': bad,
            'duration_policy_sec': 20., 'sfreq_hz': sfreq,
            'component_display': 'window_mean_subtracted_only; original IC AU; no variance normalization'}
        windows.append({'metadata': window, 'data': projected[:, left:right].copy()})
        if not with_epochs:
            continue
        selected = [row for row in rows if row['source_fif'] == record.source_filename and analysis_eligible(row)]
        if not selected:
            run_counts.append({'source': record.source_filename, 'eligible': 0, 'retained': 0})
            continue
        source = mne.io.RawArray(projected,
            mne.create_info([f'ICA{c:03}' for c in components], sfreq, 'misc'),
            first_samp=raw.first_samp, verbose='ERROR')
        source.set_meas_date(raw.info['meas_date'])
        attach_annotations(source, raw.annotations, globalize=True)
        events = np.array([[int(row['epoch_event_sample_raw']), 0, int(row['condition_code'])] for row in selected])
        epochs = mne.Epochs(source, events, event_id=cfg['epoching']['event_id_map'],
            tmin=cfg['epoching']['tmin_sec'], tmax=cfg['epoching']['tmax_sec'],
            baseline=None, detrend=None, picks='all', proj=False, reject=None,
            reject_by_annotation=bool(cfg['epoching']['reject_by_annotation']) and not include_bad,
            metadata=pd.DataFrame(selected), preload=True, on_missing='ignore', verbose='ERROR')
        if times is not None and not np.array_equal(times, epochs.times):
            raise RuntimeError('Native epoch grids differ between acquisitions; do not silently resample.')
        times = epochs.times
        values = epochs.get_data()
        # MNE does not detrend misc channels. Apply configured operators to IC AU
        # explicitly, without falsely relabelling them as EEG/EOG sensor data.
        order = cfg['epoching']['detrend']
        if len(values) and order is not None:
            from scipy.signal import detrend
            values = detrend(values, axis=-1, type='linear' if order == 1 else 'constant')
        baseline = cfg['epoching']['baseline_sec']
        if len(values) and baseline is not None:
            left_b = times[0] if baseline[0] is None else baseline[0]
            right_b = times[-1] if baseline[1] is None else baseline[1]
            picks = (times >= left_b) & (times <= right_b)
            if not picks.any():
                raise ValueError('Configured diagnostic baseline contains no native samples.')
            values -= values[..., picks].mean(axis=-1, keepdims=True)
        pieces.append(values)
        kept.extend(epochs.metadata.to_dict('records'))
        for index, reasons in enumerate(epochs.drop_log):
            if reasons:
                dropped[f'{record.source_filename}:{index}'] = list(reasons)
        run_counts.append({'source': record.source_filename, 'eligible': len(selected), 'retained': len(values)})
        del projected, source, epochs, raw
    if with_epochs and (not pieces or not any(len(piece) for piece in pieces)):
        raise RuntimeError('No diagnostic task epochs survive the requested selection.')
    values = np.concatenate(pieces) if with_epochs else None
    reasons = {}
    for row in dropped.values():
        for reason in row:
            reasons[reason] = reasons.get(reason, 0) + 1
    metadata = {'purpose': 'diagnostic_only_before_cleaning', 'timing_gate': gate,
        'drop_reasons': reasons, 'drop_log_by_source_event': dropped, 'run_counts': run_counts,
        'retained_trial_metadata': kept, 'components': components,
        'tmin_sec': float(times[0]) if times is not None else None,
        'tmax_sec': float(times[-1]) if times is not None else None,
        'epoch_samples': len(times) if times is not None else None,
        'reject_by_annotation': bool(cfg['epoching']['reject_by_annotation']) and not include_bad,
        'detrend': cfg['epoching']['detrend'], 'baseline': cfg['epoching']['baseline_sec'],
        'IC_channel_type': 'misc; configured detrend/baseline explicitly applied to activation arrays',
        'BAD_projection': 'whole operative BAD intervals; native sensor annotation assignments remain unchanged',
        'mean_band': 'pointwise normal approximation: mean +/- 1.96 SEM; no significance claim',
        'amplitude_rejection': 'not_applied_to_arbitrary_IC_units', 'model_components_removed': [],
        'continuous_windows': [w['metadata'] for w in windows]}
    return {'values': values, 'times': times, 'metadata': metadata, 'windows': windows,
        'components': components, 'kept': kept}


# -----------------------------------------------------------------------------
# Two complementary visible figures per component, from the SAME extraction
# -----------------------------------------------------------------------------

def task_figure_from_data(cfg, data, component):
    """Plot condition-specific native trial images and means; preserve exact times."""
    offset = data['components'].index(component)
    values, times = data['values'][:, offset], data['times']
    codes = np.array([int(row['condition_code']) for row in data['kept']])
    conditions = list(cfg['epoching']['event_id_map'].items())
    fig, axes = plt.subplots(2, len(conditions), figsize=(15, 7), sharex=True, squeeze=False,
        gridspec_kw={'height_ratios': [2, 1]})
    limit = max(float(np.percentile(np.abs(values), 99.5)), np.finfo(float).tiny)
    counts = {}
    for column, (name, code) in enumerate(conditions):
        selected = values[codes == int(code)]
        counts[name] = len(selected)
        axes[0, column].set_title(f'{name}; n={len(selected)}', fontsize=10)
        if len(selected):
            axes[0, column].imshow(selected, origin='lower', aspect='auto', cmap='RdBu_r', vmin=-limit, vmax=limit,
                extent=[times[0], times[-1], .5, len(selected) + .5])
            mean = selected.mean(axis=0)
            ci = 1.96 * selected.std(axis=0, ddof=1) / np.sqrt(len(selected)) if len(selected) > 1 else np.zeros_like(mean)
            axes[1, column].plot(times, mean, color='k')
            axes[1, column].fill_between(times, mean - ci, mean + ci, color='0.7', alpha=.5)
        for ax in axes[:, column]:
            ax.axvline(0, color='tab:green', linestyle='--')
            ax.set_xlim(times[0], times[-1])
        axes[1, column].set_xticks([times[0], 0., .5, 1., 1.5, 2., times[-1]])
        axes[1, column].set_xlabel('Seconds from accepted motion onset')
    axes[0, 0].set_ylabel('Diagnostic trial')
    axes[1, 0].set_ylabel('IC activation (AU)')
    meta = {**data['metadata'], 'kind': 'task_locked_native_IC', 'component': component, 'counts': counts,
        'display_clim': '99.5_percentile_absolute_activation; no trials removed by this display scale'}
    fig.suptitle(f"{cfg['layout']['subject_bids_id']} ICA{component:03}: TASK-LOCKED DIAGNOSTIC, before cleaning\n"
        f"BAD rejection={meta['reject_by_annotation']}; detrend={meta['detrend']}; baseline={meta['baseline']}; no production epochs")
    fig.tight_layout()
    return fig, meta


def continuous_figure_from_data(cfg, data, component):
    """Show one deterministic native continuous example per run, with BAD labels."""
    offset = data['components'].index(component)
    fig, axes = plt.subplots(len(data['windows']), 1, figsize=(15, max(3, 2.1 * len(data['windows']))), squeeze=False)
    traces = [w['data'][offset] - w['data'][offset].mean() for w in data['windows']]
    limit = max(max(float(np.max(np.abs(v))) for v in traces) * 1.05, np.finfo(float).eps)
    for ax, window, trace in zip(axes[:, 0], data['windows'], traces, strict=True):
        meta = window['metadata']
        times = (np.arange(len(trace)) + meta['window_start_sample_raw']) / meta['sfreq_hz']
        ax.plot(times, trace, linewidth=.6, color='k')
        labels = sorted({a['description'] for a in meta['bad_annotations']})
        colors = {label: plt.get_cmap('tab10')(i % 10) for i, label in enumerate(labels)}
        for a in meta['bad_annotations']:
            ax.axvspan(a['onset_acquisition_sec'], a['onset_acquisition_sec'] + a['duration_sec'],
                color=colors[a['description']], alpha=.2)
        if labels:
            ax.legend(handles=[Patch(facecolor=colors[label], alpha=.4, label=label) for label in labels],
                fontsize=6, loc='upper right')
        ax.set_xlim(times[0], times[-1]); ax.set_ylim(-limit, limit)
        ax.set_ylabel('IC (AU)')
        ax.set_title(meta['source_filename'] + ' | first ocular example / fixed fallback', fontsize=9)
        ax.set_xlabel('Acquisition seconds (absolute native sample / sfreq)')
    fig.suptitle(f"{cfg['layout']['subject_bids_id']} ICA{component:03}: NATIVE CONTINUOUS examples, before cleaning\n"
        '20 s per run; BAD periods included; window mean subtracted for display; shared y-scale; not task-locked', fontsize=11)
    fig.tight_layout()
    return fig, {'kind': 'native_continuous_IC_examples', 'component': component,
        'windows': data['metadata']['continuous_windows'], 'components_removed': [],
        'signal_preparation': 'native annotated SSS input; no new temporal filter or resampling',
        'display': 'window mean subtracted only; shared per-component y-scale across runs; IC AU'}


def generate_native_diagnostics(cfg, model, records, inputs, figure_dir, model_sha256):
    """Render every final IC's continuous/task views from a single native pass."""
    raw_root = cfg['resolved_paths']['data_root']
    figure_dir = guarded_path(figure_dir, raw_root)
    figure_dir.mkdir(parents=True, exist_ok=True)
    data = extract_native_diagnostics(cfg, model, records, inputs)
    manifest = figure_dir / '03_C_ICA_review_native_diagnostics.json'
    shared = {**data['metadata'], 'model_sha256': model_sha256,
        'continuous_and_task_default': True, 'components_complete': data['components']}
    write_json(manifest, shared, raw_root=raw_root)
    images = []
    for component in data['components']:
        print(f"Native review figures: ICA{component:03}/{model.n_components_ - 1:03}", flush=True)
        for name, renderer in (('continuous', continuous_figure_from_data), ('task', task_figure_from_data)):
            fig, meta = renderer(cfg, data, component)
            path = figure_dir / f'03_C_ICA_review_ICA{component:03}_{name}_01.png'
            fig.savefig(path, dpi=120); plt.close(fig)
            # Keep shared trial identities once, avoiding duplicate large JSONs.
            sidecar = {k: v for k, v in meta.items() if k not in
                ('retained_trial_metadata', 'drop_log_by_source_event', 'continuous_windows')}
            write_json(path.with_suffix('.json'), {**sidecar, 'model_sha256': model_sha256,
                'shared_manifest': manifest.name}, raw_root=raw_root)
            images.append(path)
    return images, shared
