"""Synthetic storage-saving epoch tests; no participant data or jobs are used.

Usage: python -m pytest meg/tests/test_inline_epoch_preparation.py -q

Compare continuous-in-memory preparation with a reference continuous-FIF route,
including actual spectrum fitting, real saved ICA reuse, annotations, acquisition
coordinates, auxiliary channels, trial identities and float32 serialization.
Guard tests cover missing reviews, model tampering, immutable run manifests,
split-file binding and changed settings. Only temporary fixtures are written.
"""
from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import mne
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code/preprocessing_basic'))
import epoch_input_preparation as prep
from continuous_cleaning import annotation_rows, check_preservation, file_hash, filter_continuous
from test_epoch_export_smoke import EPOCHING, EXPORT, _event_row
from test_stage4_continuous import raw_example


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    """Create annotated native data, trial events and minimal accepted-guard seam."""
    source_dir = tmp_path / 'raw/sub01_date/meg'; source_dir.mkdir(parents=True)
    derivatives = tmp_path / 'derivatives'
    upstream_reports = derivatives / 'reports/sub-01/ses-01'; upstream_reports.mkdir(parents=True)
    preprocessing = derivatives / 'preprocessing/sub-01/ses-01'
    sss = preprocessing / 'sss'; sss.mkdir(parents=True)
    events_dir = preprocessing / 'events'; events_dir.mkdir()
    filename = 'sub01_run01_block01_date.fif'
    raw = raw_example(duration=60.)
    raw.set_annotations(mne.Annotations([20., 35.], [2., 0.], ['BAD_blink', 'EDGE boundary']))
    times = np.arange(raw.n_times) / raw.info['sfreq']
    raw._data[:3] += 1e-11 * (np.sin(2*np.pi*50*times) + np.sin(2*np.pi*100*times))
    raw.save(source_dir / filename, overwrite=True, verbose='ERROR')
    ann = sss / (Path(filename).stem + '_ann.fif'); raw.save(ann, overwrite=True, verbose='ERROR')
    rows = []
    for trial, seconds in enumerate((10., 20., 30., 40.), 1):
        sample = raw.first_samp + int(seconds*raw.info['sfreq'])
        rows.extend([_event_row(code=31, sample=sample, label='trial_condition_onset_always_visible', trial=trial, source_name=filename),
                     _event_row(code=1 if trial < 3 else 2, sample=sample+10, label='trial_sequence_identity', trial=trial, source_name=filename)])
        if trial == 4:
            rows.append(_event_row(code=150, sample=sample+20, label='gaze_break', trial=trial, source_name=filename))
    pd.DataFrame(rows).to_csv(events_dir / '02_trigger_events.csv', index=False)
    manifest = upstream_reports / 'source_manifest.json'
    manifest.write_text(json.dumps({'derivative_lineage': {'schema_version': '1.0', 'stages': {'upstream': {'status': 'completed'}}}}))
    trajectory = source_dir.parent / 'input_files'; trajectory.mkdir()
    (trajectory / 'MovDot_Sub01_Vtest_rescueTraject.mat').write_bytes(b'synthetic trajectory')
    channels = tmp_path / 'canonical.txt'; channels.write_text('MEG0111\nMEG0121\nMEG0311\n')
    run_root = derivatives / 'processing_runs/noica-fixture'
    cfg = {'subject': {'id': '01', 'versions': {'stimuli': 'Vtest'}},
        'subject_profile': {'identity_range': [1, 2], 'iti_sec_range': [1.8, 2.]},
        'layout': {'subject_bids_id': 'sub-01', 'session_bids_id': 'ses-01'},
        'resolved_paths': {'project_root': str(tmp_path), 'data_root': str(tmp_path/'raw'),
            'subject_data_dir': str(source_dir.parent), 'subject_input_files_dir': str(trajectory),
            'derivatives_root': str(derivatives), 'derivatives_preprocessing_dir': str(preprocessing),
            'processing_run_root': str(run_root),
            'upstream_derivatives_reports_dir': str(upstream_reports),
            'upstream_derivatives_analysis_dir': str(derivatives/'analysis/sub-01/ses-01'),
            'derivatives_reports_dir': str(run_root/'reports/sub-01/ses-01'),
            'derivatives_analysis_dir': str(run_root/'analysis/sub-01/ses-01')},
        'data_discovery': {'meg_fif_glob': 'meg/sub01_run*.fif', 'trajectory_mat_glob': 'input_files/MovDot_Sub01_V*_rescueTraject.mat'},
        'stage4_review': {'03_A_MaxFilt': {'evidence_sha256': 'sss'}, '03_B_ArtefactAnnotation': {'evidence_sha256': 'ann'}},
        'preprocessing': {'highpass_hz': None, 'lowpass_hz': None, 'resample_sfreq_hz': None, 'notch_freqs_hz': [50., 100.],
            'notch_method': 'spectrum_fit', 'spectrum_fit_window': '10s', 'filter_edge_margin_sec': None},
        'ica': {'manual_exclude': [1], 'reuse_prefit': True},
        'epoching': {'preparation_mode': 'inline_no_ica', 'processing_run': 'noica-fixture',
            'event_id_map': {'always_visible': 31}, 'tmin_sec': -.5, 'tmax_sec': 1., 'baseline_sec': None,
            'detrend': 1, 'reject_by_annotation': True, 'amplitude_reject': {'mag': 1e-10},
            'trial_exclusion': {'gaze_break_codes': [150], 'catch_type_1_codes': [116, 117], 'catch_type_2_codes': [113, 118, 201, 202],
                'esc_termination_codes': [152], 'exclude_gaze_break': True, 'exclude_catch_type_1': True,
                'exclude_catch_type_2': True, 'exclude_esc_termination': True, 'exclude_replay_trials': False, 'exclude_unresolved_transitions': True}},
        'export': {'canonical_meg_channels_file': str(channels), 'tmin_sec': 0., 'tmax_sec': 1., 'sampling_rate_hz': 50.,
            'missing_identity_policy': 'error', 'bad_or_missing_channel_policy': 'nan'}}
    # Unit fixtures replace only acceptance receipts; production review/source
    # validation is exercised separately by the all-participant preflight.
    calls = []
    monkeypatch.setattr(prep, 'require_review', lambda cfg, stage, reports: calls.append(stage))
    return cfg, manifest, ann, calls


def saved_model(cfg, ann, monkeypatch):
    """Write a real tiny model; isolate only the upstream acceptance/cache seam."""
    model = mne.preprocessing.ICA(n_components=3, random_state=97, max_iter=1000)
    raw = mne.io.read_raw_fif(ann, preload=True, verbose='ERROR')
    model.fit(raw.copy().filter(1., 40., verbose='ERROR'), picks='mag', verbose='ERROR')
    root = Path(cfg['resolved_paths']['derivatives_preprocessing_dir'])/'ica_models'; root.mkdir()
    path = root/'fixture-ica.fif'; model.save(path, overwrite=True, verbose='ERROR')
    sha = file_hash(path)
    path.with_name(path.stem+'_metadata.json').write_text(json.dumps({'model_sha256': sha}))
    report = {'prefit_model': str(path), 'prefit_model_sha256': sha, 'configuration_sha256': 'bound', 'n_components_fitted': 3}
    (Path(cfg['resolved_paths']['upstream_derivatives_reports_dir'])/'03_C_ICA_review_report.json').write_text(json.dumps(report))
    cfg['ica']['reviewed_model_sha256'] = sha
    helpers = prep.ica_stage()
    monkeypatch.setattr(helpers, 'build_ica_fingerprint', lambda **kwargs: {'sources': 'bound'})
    monkeypatch.setattr(helpers, 'cache_is_compatible', lambda metadata, expected: True)
    monkeypatch.setattr(prep, 'stage_configuration', lambda *args, **kwargs: 'bound')
    monkeypatch.setattr(mne.preprocessing.ICA, 'fit', lambda *args, **kwargs: pytest.fail('Production must never fit ICA'))
    return path


@pytest.mark.parametrize('mode', ['inline_no_ica', 'inline_ica'])
@pytest.mark.parametrize('tiny_fragment', [False, True])
def test_real_inline_epoch_export_matches_reference_and_preserves_inputs(inputs, monkeypatch, mode, tiny_fragment):
    """Continuous ordering, native preservation, trial identity and FIF parity."""
    cfg, manifest, ann, calls = inputs
    cfg['epoching']['preparation_mode'] = mode
    if tiny_fragment:
        # A five-sample island outside candidate epochs reproduces the real
        # cohort failure without changing any trial's scientific eligibility.
        native = mne.io.read_raw_fif(ann, preload=True, verbose='ERROR')
        gap_sec = 5. / native.info['sfreq']
        native.annotations.append([native.first_time+45., native.first_time+45.5+gap_sec],
            [.5, .5-gap_sec], ['BAD_blink', 'BAD_muscle'])
        native.save(ann, overwrite=True, verbose='ERROR')
    if mode == 'inline_ica':
        saved_model(cfg, ann, monkeypatch)
    original_manifest = manifest.read_bytes(); original_ann = prep.file_stamp(ann)
    runs, gate, model = prep.resolve_epoch_inputs(cfg, manifest)
    assert ('03_C_ICA' in calls) == (mode == 'inline_ica')
    assert cfg['ica']['manual_exclude'] == [1]  # no-ICA must not rewrite choices
    before = mne.io.read_raw_fif(ann, preload=True, verbose='ERROR')
    reference = before.copy()
    if model is not None:
        model.apply(reference, exclude=[1], verbose='ERROR')
    if tiny_fragment:
        # The old public Raw filter supplies an independent reference when its
        # disposable mask explicitly omits just the unsupported island.
        reference.annotations.append(reference.first_time+45.5, 5./reference.info['sfreq'], 'BAD_reference_fragment')
    filter_continuous(reference, cfg['preprocessing'])
    reference._annotations = before.annotations.copy()
    prepared, qc = prep.prepare_run(cfg, runs[0], gate['preparation'], model, Path(cfg['resolved_paths']['derivatives_reports_dir'])/'qc')
    np.testing.assert_array_equal(prepared.get_data(), reference.get_data())
    check_preservation(before, prepared)
    assert annotation_rows(before) == annotation_rows(prepared)
    assert qc['preservation']['auxiliary_samples_verified'] and not qc['continuous_file_saved']
    assert len(qc['final_filter']['spectrum_fit']['skipped_fragments']) == int(tiny_fragment)
    if tiny_fragment:
        assert qc['final_filter']['spectrum_fit']['skipped_fragments'][0]['n_samples'] == 5
    report = EPOCHING.run_epoching(cfg, source_manifest=manifest, overwrite=False)
    epochs = mne.read_epochs(report['outputs']['epochs_fif'], preload=True, verbose='ERROR')
    assert report['counts']['n_final_epochs'] == 2
    assert report['counts']['trial_policy_exclusions'] == {'gaze_break': 1}
    assert report['counts']['mne_drop_reasons'] == {'BAD_blink': 1}
    assert epochs.metadata['trial_sequence_identity'].tolist() == [1, 2]
    decisions = pd.read_csv(report['outputs']['trial_decisions_csv'])
    candidates = decisions[decisions.included_by_trial_policy]
    events = np.array([[int(row.raw_sample_frame), 0, 31] for row in candidates.itertuples()])
    # Reference route adds the historical continuous float32 FIF save/read.
    refpath = ann.parent/'reference_ica.fif'; reference.save(refpath, overwrite=True, verbose='ERROR')
    ref = mne.Epochs(mne.io.read_raw_fif(refpath, preload=False, verbose='ERROR'), events,
        event_id={'always_visible': 31}, tmin=-.5, tmax=1., baseline=None, detrend=1,
        proj=False, picks='all', reject={'mag': 1e-10}, reject_by_annotation=True, preload=True, verbose='ERROR')
    assert ref.selection.tolist() == [0, 2]
    # Different placement of float32 serialization accounts for roundoff only.
    np.testing.assert_allclose(epochs.get_data()[:,:3], ref.get_data()[:,:3], rtol=3e-6, atol=2e-18)
    assert prep.verify_epoch_artifacts(cfg)['status'] == 'completed'
    exported = EXPORT.run_export(cfg, source_manifest=manifest, overwrite=False)
    assert prep.verify_export_artifacts(cfg)['export_contract_sha256'] == exported['export_contract_sha256']
    from audit_ica_application import audit_exports
    audits = audit_exports(cfg, {'reports': Path(cfg['resolved_paths']['derivatives_reports_dir'])}, epochs)
    assert audits[0]['counts_by_identity'] == {1: 1, 2: 1}
    assert manifest.read_bytes() == original_manifest and prep.file_stamp(ann) == original_ann
    assert not list(ann.parent.glob('sub*_ica.fif'))
    assert prep.processing_manifest(cfg, manifest).is_file()
    with pytest.raises(RuntimeError, match='already exist'):
        EPOCHING.run_epoching(cfg, source_manifest=manifest, overwrite=True)
    changed = deepcopy(cfg); changed['epoching']['detrend'] = 0
    with pytest.raises(RuntimeError, match='settings differ'):
        prep.verify_epoch_artifacts(changed)
    # Every split part is checked, not just the first FIF file.
    split = Path(report['outputs']['epochs_fif']).with_name('fake_epo-1.fif'); split.write_bytes(b'split')
    payload = json.loads(Path(report['outputs']['report_json']).read_text())
    payload['artifact_bindings']['epochs_fif']['parts'].append(prep.file_stamp(split))
    Path(report['outputs']['report_json']).write_text(json.dumps(payload)); split.write_bytes(b'changed split')
    with pytest.raises(RuntimeError, match='split part changed'):
        prep.verify_epoch_artifacts(cfg)


def test_inline_guards_and_contract_changes(inputs, monkeypatch):
    cfg, manifest, ann, calls = inputs
    first = prep.resolve_epoch_inputs(cfg, manifest)[1]
    cfg['ica']['manual_exclude'] = [0, 2]
    second = prep.resolve_epoch_inputs(cfg, manifest)[1]
    assert first == second  # dormant proposals are not an applied-ICA choice
    cfg['preprocessing']['highpass_hz'] = .01
    with pytest.raises(RuntimeError, match='established native'):
        prep.resolve_epoch_inputs(cfg, manifest)
    cfg['preprocessing']['highpass_hz'] = None
    target = prep.processing_manifest(cfg, manifest, create=True)
    assert target != manifest
    manifest.write_text(manifest.read_text()+'\n')
    with pytest.raises(RuntimeError, match='ancestry changed'):
        prep.processing_manifest(cfg, manifest)
    cfg['epoching']['processing_run'] = None
    with pytest.raises(RuntimeError, match='explicit mode'):
        prep.resolve_epoch_inputs(cfg, manifest)


def test_inline_ica_blocks_changed_exact_model(inputs, monkeypatch):
    cfg, manifest, ann, calls = inputs
    cfg['epoching']['preparation_mode'] = 'inline_ica'
    path = saved_model(cfg, ann, monkeypatch)
    authorization = Path(cfg['resolved_paths']['derivatives_root'])/'authorization.json'
    authorization.write_text(json.dumps({'participants': [{'subject': 1,
        'model_sha256': file_hash(path), 'component_count': 3, 'exclusions': [1],
        'deferred_components_retained': [2]}]}))
    cfg['stage4_review']['03_C_ICA'] = {'acceptance_kind': 'user_provisional_component_processing_authorization',
        'authorization_relative_path': authorization.name, 'authorization_sha256': file_hash(authorization)}
    prep.resolve_epoch_inputs(cfg, manifest)
    cfg['ica']['manual_exclude'] = [0]
    with pytest.raises(RuntimeError, match='recorded exact-model authorization'):
        prep.resolve_epoch_inputs(cfg, manifest)
    cfg['ica']['manual_exclude'] = [1]
    with path.open('ab') as handle:
        handle.write(b'changed model')
    with pytest.raises(RuntimeError, match='source/model/configuration'):
        prep.resolve_epoch_inputs(cfg, manifest)


def test_processing_namespace_and_environment_are_shared_across_consumers(tmp_path, monkeypatch):
    """Config overrides isolate only downstream paths and reject unsafe labels."""
    import pipeline_config
    import yaml
    from test_pipeline_config import _minimal_config
    cfg = _minimal_config(tmp_path)
    (tmp_path/'raw/sub01_date').mkdir(parents=True)
    path = tmp_path/'defaults.yaml'; path.write_text(yaml.safe_dump(cfg))
    monkeypatch.setenv('DAD_EPOCH_PREPARATION', 'inline_no_ica')
    monkeypatch.setenv('DAD_PROCESSING_RUN', 'noica-test')
    resolved, *_ = pipeline_config.load_and_resolve_config(path)
    paths = resolved['resolved_paths']
    assert paths['derivatives_reports_dir'].endswith('processing_runs/noica-test/reports/sub-01/ses-01')
    assert paths['upstream_derivatives_reports_dir'].endswith('outputs/reports/sub-01/ses-01')
    assert 'processing_runs' not in paths['derivatives_preprocessing_dir']
    assert prep.upstream_config(resolved)['resolved_paths']['derivatives_reports_dir'] == paths['upstream_derivatives_reports_dir']
    monkeypatch.setenv('DAD_PROCESSING_RUN', '../../escape')
    with pytest.raises(ValueError, match='label'):
        pipeline_config.load_and_resolve_config(path)


def test_erf_namespace_and_resume_bind_cohort_and_settings(inputs):
    """ERF variants cannot reuse another epoch cohort, baseline or permutation set."""
    import importlib.util
    path = Path(__file__).resolve().parents[1]/'code/sanity_checks/05_dot_erf_analysis.py'
    spec = importlib.util.spec_from_file_location('dad_inline_erf_test', path)
    erf = importlib.util.module_from_spec(spec); spec.loader.exec_module(erf)
    cfg, *_ = inputs
    output = erf.group_output_paths(cfg, 'mag', 'occipital')
    assert str(output['out_dir']).startswith(cfg['resolved_paths']['processing_run_root'])
    assert str(output['report_dir']).startswith(cfg['resolved_paths']['processing_run_root'])
    args = SimpleNamespace(baseline=[-.3, -.1], analysis_sfreq=200., n_permutations=10000, random_seed=42, resume=False)
    headers = [{'subject': 'sub-01', 'stage04_validation': {'preparation_contract_sha256': 'noica'}}]
    first = erf.erf_resume_contract(args, headers, ['MEG0111'], 'occipital')
    args.resume = True
    assert first == erf.erf_resume_contract(args, headers, ['MEG0111'], 'occipital')
    headers[0]['stage04_validation']['preparation_contract_sha256'] = 'with-ica'
    assert first != erf.erf_resume_contract(args, headers, ['MEG0111'], 'occipital')
    headers[0]['stage04_validation']['preparation_contract_sha256'] = 'noica'
    args.n_permutations = 9999
    assert first != erf.erf_resume_contract(args, headers, ['MEG0111'], 'occipital')
