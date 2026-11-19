#!/usr/bin/env python3
"""Preserve an ICA optimizer result before shared-source revalidation.

Usage from reconcile_ica_convergence.py (use local/allocated scratch storage)::

    save_checkpoint(new_directory, model, context, warnings, elapsed, raw_root)
    model, metadata = load_checkpoint(directory, expected_context)

Inputs: an actual fitted MNE ICA object, the frozen original model/source/fit
context, original warning records and elapsed fitting time. Outputs: a NEW
checkpoint_ica.fif and checkpoint.json. The checkpoint is always explicitly
UNVALIDATED and cannot establish source, component or cleaned-data acceptance.
Its purpose is to retain completed computation when a network volume disappears.

The caller must perform its ordinary full source/model/accepted-policy validation
before loading, then reconstruct the same fit copies and repeat all numerical
checks and final source validation. Context/hash/model-summary mismatches fail.
This helper never writes a production model, tuning, receipt, annotation or epoch.
The numerical fit remains one participant dependency; participant parallelism
and thread/resource controls belong to the calling workflow.
"""
from pathlib import Path
import json

import mne

from continuous_cleaning import file_hash
from preprocessing_stage_utils import guarded_path, utc_now, write_json


# -----------------------------------------------------------------------------
# Checkpoints describe actual serialized state, not an acceptance claim
# -----------------------------------------------------------------------------

def iteration_cap(model):
    """Read the actual optimizer cap, including MNE's FIF read-back behavior.

    MNE 1.10.2 read_ica reconstructs the max_iter attribute with its default
    while retaining the fitted cap in fit_params. That persisted optimizer
    parameter, which was passed to sklearn, governs convergence accounting.
    """
    value = model.fit_params.get('max_iter', model.max_iter)
    if not isinstance(value,int) or isinstance(value,bool) or value<1:
        raise RuntimeError('ICA optimizer has no valid recorded integer iteration cap.')
    return value


def model_summary(model):
    """Return exact optimizer/accounting fields that must survive MNE I/O."""
    return {'method':model.method, 'n_components':int(model.n_components_),
        'n_samples':int(model.n_samples_), 'n_iterations':int(model.n_iter_),
        'max_iter':iteration_cap(model), 'channels':model.ch_names,
        'fit_params':model.fit_params, 'exclude':model.exclude}


def save_checkpoint(directory, model, context, warning_records, elapsed, raw_root):
    """Write a NEW unvalidated checkpoint outside operative derivative outputs."""
    directory = guarded_path(directory, raw_root)
    directory.mkdir(parents=True, exist_ok=False)
    path = directory/'checkpoint_ica.fif'
    model.save(path, overwrite=False, verbose='ERROR')
    metadata = {'schema_version':1, 'generated_utc':utc_now(),
        'kind':'unvalidated_optimizer_checkpoint', 'source_revalidation_required':True,
        'scientific_acceptance':False, 'production_model_changed':False,
        'model_sha256':file_hash(path), 'model_summary':model_summary(model),
        'context':context, 'warnings':warning_records, 'elapsed_fit_sec':float(elapsed)}
    write_json(directory/'checkpoint.json', metadata, raw_root=raw_root, overwrite=False)
    print(f'Preserved UNVALIDATED optimizer checkpoint: {directory}', flush=True)
    return metadata


# -----------------------------------------------------------------------------
# Resume only the same frozen computation; all ordinary gates stay with caller
# -----------------------------------------------------------------------------

def load_checkpoint(directory, expected_context):
    """Load matching actual weights; this does not validate current sources."""
    directory = Path(directory)
    metadata = json.loads((directory/'checkpoint.json').read_text())
    if (metadata.get('kind') != 'unvalidated_optimizer_checkpoint' or
        metadata.get('source_revalidation_required') is not True or
        metadata.get('scientific_acceptance') is not False or
        metadata.get('production_model_changed') is not False):
        raise RuntimeError('Checkpoint is not an explicitly unvalidated optimizer result.')
    if metadata.get('context') != expected_context:
        raise RuntimeError('Checkpoint source/model/configuration/runtime context changed.')
    path = directory/'checkpoint_ica.fif'
    if file_hash(path) != metadata['model_sha256']:
        raise RuntimeError('Checkpoint model bytes changed.')
    model = mne.preprocessing.read_ica(path, verbose='ERROR')
    if model_summary(model) != metadata['model_summary']:
        raise RuntimeError('Checkpoint model accounting/optimizer settings differ from metadata.')
    if model.exclude:
        raise RuntimeError('Numerical repair checkpoint must contain no component decisions.')
    return model,metadata
