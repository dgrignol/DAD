#!/usr/bin/env python3
"""Audit diagnostic IC epochs against the accepted operative trial footprint.

Usage (bounded participant processes, no fitting/application/acceptance)::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_ica_task_diagnostics.py \
        --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
        --label native-task-retention-audit-20261110-v1 --jobs 2

Inputs: the immutable current-mask plan and joined policy receipt, current
model/report/native-diagnostic manifests for P04-P14/P16-P39. Apply the receipt's
explicit FIF serialization trial deltas to the plan's original footprint; then
compare exact retained trial identities, condition/run/identity counts, accepted
Stage 2/3 timing gate and native motion-onset window bounds. Check current model
content hashes and bound input/manifest file stamps. This supplements the full
source/header validation performed by the report workflow; it cannot replace
those gates or establish final epoch/cleaned-output acceptance.

Outputs: one NEW reports/ica_review_updates/LABEL/summary.json after every
participant joins. It records metadata agreement, not task-response validity or
scientific component acceptance. Native sources, models, decisions, annotations,
accepted receipts and all earlier snapshots remain unchanged. --jobs controls
separate participant processes; numerical-library threads are limited by the
shared participant runner. Print one progress line per completed participant.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from participant_parallel import map_participants
from preprocessing_stage_utils import utc_now, write_json

COHORT = list(range(4,15)) + list(range(16,40))


# -----------------------------------------------------------------------------
# Portable bindings and exact retained identities (including serialization deltas)
# -----------------------------------------------------------------------------

def checksum(path):
    """Hash a model or small manifest, never an entire cohort's raw signals."""
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            result.update(block)
    return result.hexdigest()


def portable(value):
    """Use the accepted plan's path-independent timing-gate convention."""
    if isinstance(value,dict):
        return {k:portable(v) for k,v in value.items() if not k.endswith('_path')}
    if isinstance(value,list):
        return [portable(v) for v in value]
    return value


def audit_participant(item):
    """Return exact diagnostic-retention agreement for one disjoint participant."""
    root,planned,accepted = item
    subject = planned['subject']
    reports = root/f'reports/sub-{subject:02}/ses-01'
    directory = root/f'preprocessing/sub-{subject:02}/ses-01'
    report = json.loads((reports/'03_C_ICA_review_report.json').read_text())
    manifest_path = directory/'figures/03_C_ICA_review_native_diagnostics.json'
    native = json.loads(manifest_path.read_text())
    model = directory/'ica_models'/Path(report['prefit_model']).name
    if checksum(model) != report['prefit_model_sha256'] or native['model_sha256'] != report['prefit_model_sha256']:
        raise RuntimeError(f'P{subject}: current report/native/model identities disagree.')
    bound = next(row['stamp'] for row in report['evidence_files']
                 if row['relative_path'] == str(manifest_path.relative_to(root)))
    stat = manifest_path.stat()
    if bound != {'name':manifest_path.name,'size':stat.st_size,'mtime_ns':stat.st_mtime_ns}:
        raise RuntimeError(f'P{subject}: native manifest changed since report publication.')
    for row in report['fingerprint']['inputs']:
        stat = (directory/'sss'/row['source_filename']).stat()
        if (stat.st_size,stat.st_mtime_ns) != (row['size_bytes'],row['mtime_ns']):
            raise RuntimeError(f'P{subject}: an annotated source changed.')
    expected = {r['trial_key']:r for r in planned['expected_trials']}
    for delta in accepted['serialization_outcome_changes']:
        if expected[delta['trial_key']] != delta['before_serialization']:
            raise RuntimeError(f'P{subject}: serialization delta has a different starting trial.')
        expected[delta['trial_key']] = delta['operative']
    expected_keys = {key for key,row in expected.items() if row['analysis_eligible'] and row['annotation_survives']}
    rows = native['retained_trial_metadata']
    actual_keys = {row['trial_key'] for row in rows}
    if len(actual_keys) != len(rows) or actual_keys != expected_keys:
        raise RuntimeError(f'P{subject}: diagnostic retained trial identities disagree with accepted operative masks.')
    if portable(native['timing_gate']) != planned['timing']['gate']:
        raise RuntimeError(f'P{subject}: diagnostic timing gate differs from accepted Stage 2/3 evidence.')

    # Compare every cell, including cells with zero surviving trials. IC AU are
    # not subject to sensor amplitude thresholds; this is annotation-only parity.
    labels = {row['condition_code']:row['condition'] for row in accepted['condition_retention']}
    conditions = Counter(int(row['condition_code']) for row in rows)
    runs = Counter((row['source_fif'],labels[int(row['condition_code'])]) for row in rows)
    identities = Counter((labels[int(row['condition_code'])],int(float(row['trial_sequence_identity']))) for row in rows)
    assert all(conditions[row['condition_code']] == row['annotation_survivors'] for row in accepted['condition_retention'])
    assert all(runs[row['source_fif'],row['condition']] == row['annotation_survivors'] for row in accepted['run_condition_retention'])
    assert all(identities[row['condition'],int(row['identity'])] == row['annotation_survivors'] for row in accepted['identity_retention'])
    assert len(rows) == accepted['annotation_survivors'] == sum(row['retained'] for row in native['run_counts'])
    assert native['reject_by_annotation'] and native['model_components_removed'] == []
    assert native['tmin_sec'] == -.5 and native['tmax_sec'] == 2.667
    for row in rows:
        anchor = int(row['epoch_event_sample_raw'])
        sfreq = float(row['sfreq_hz'])
        first = anchor + round(native['tmin_sec']*sfreq)
        last = anchor + round(native['tmax_sec']*sfreq)
        assert int(row['proposed_epoch_first_sample_raw']) == first
        assert int(row['proposed_epoch_last_sample_raw']) == last
        assert last-first+1 == native['epoch_samples']
        assert first >= int(row['recording_first_sample_raw'])
        assert last < int(row['recording_end_sample_raw_exclusive'])
    print(f'sub-{subject:02}: exact {len(rows)} retained trials; condition/run/identity cells and accepted native anchors agree',flush=True)
    return {'subject':subject,'retained_trials':len(rows),'condition_counts':dict(conditions),
        'model_sha256':native['model_sha256'],'epoch_samples':native['epoch_samples'],
        'retained_trial_keys_sha256':hashlib.sha256(json.dumps(sorted(actual_keys)).encode()).hexdigest(),
        'native_manifest_sha256':checksum(manifest_path),'status':'diagnostic_metadata_agreement_not_final_epochs'}


# -----------------------------------------------------------------------------
# Join every participant before publishing a new cohort evidence snapshot
# -----------------------------------------------------------------------------

def main():
    """Audit the entire working cohort with explicit resources and immutable label."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--derivatives-root',type=Path,required=True)
    parser.add_argument('--plan-label',default='current-masks-20261109-v1')
    parser.add_argument('--label',required=True)
    parser.add_argument('--jobs',type=int,default=2)
    args = parser.parse_args()
    if not args.label.replace('-','').replace('_','').isalnum():
        parser.error('Use a new simple output label.')
    root = args.derivatives_root.resolve()
    plan_path = root/'stage4b_adoption'/args.plan_label/'plan.json'
    receipt_path = plan_path.with_name('policy_acceptance.json')
    plan = json.loads(plan_path.read_text()); receipt = json.loads(receipt_path.read_text())
    assert plan['plan_sha256'] == receipt['adoption_plan_sha256']
    assert receipt['status'] == 'accepted_for_final_ica_fit'
    assert sorted(row['subject'] for row in receipt['participants']) == COHORT
    accepted = {row['subject']:row for row in receipt['participants']}
    rows = map_participants(audit_participant,[(root,row,accepted[row['subject']])
        for row in plan['participants'] if row['subject'] in COHORT],args.jobs)
    assert sorted(row['subject'] for row in rows) == COHORT
    output = root/'reports/ica_review_updates'/args.label/'summary.json'
    write_json(output,{'generated_utc':utc_now(),'participants':rows,'participant_count':len(rows),
        'retained_trials':sum(row['retained_trials'] for row in rows),
        'resources':{'participant_processes':args.jobs,'numerical_threads_per_worker':1},
        'policy_acceptance_sha256':checksum(receipt_path),'plan_sha256':checksum(plan_path),
        'scope':'exact operative retained trial identities/cells, accepted timing gate and native metadata bounds; current model hashes and bound manifest/input stamps; no repeat raw-header audit',
        'scientific_acceptance':False,'production_epochs_created':False},
        raw_root=root.parent/'data',overwrite=False)
    print('JOINED',output,flush=True)


if __name__ == '__main__':
    main()
