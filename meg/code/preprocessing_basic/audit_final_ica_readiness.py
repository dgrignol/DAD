#!/usr/bin/env python3
"""Join current final-ICA report readiness and exact diagnostic trial retention.

Usage after all participant report/reconciliation work has completed::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_final_ica_readiness.py \
      --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
      --label NEW-LABEL --jobs 4

Read all 35 working participants independently in bounded separate processes.
Compare the 34 unchanged model hashes to the immutable array-11197 audit and
P19 to its fully published numerical-repair result. Check actual model hashes,
complete property/continuous/task coverage, sidecars, bound evidence stamps,
unreviewed decisions and iteration accounting. Reuse the exact accepted Stage
4B/Stage 2/3 diagnostic-retention audit, including FIF serialization deltas.

Only after all participants join, write ONE new immutable summary. This is
software/numerical readiness for human component review, not scientific ICA or
cleaned-output acceptance. No raw-header reread, full PNG/HTML hash reread,
model fitting, signal application, decision or production epoch is performed.
--jobs exposes process allocation; participant_parallel limits library threads.
"""
import argparse
import json
from pathlib import Path

from audit_ica_task_diagnostics import COHORT, audit_participant, checksum
from continuous_cleaning import file_stamp
from participant_parallel import map_participants
from preprocessing_stage_utils import utc_now, write_json


# -----------------------------------------------------------------------------
# Each worker has read-only, disjoint participant evidence; shared join is last
# -----------------------------------------------------------------------------

def audit(item):
    """Check current model/report panels plus exact operative native trial rows."""
    root, planned, accepted, baseline = item
    subject = planned['subject']
    directory = root/f'reports/sub-{subject:02}/ses-01'
    payload = json.loads((directory/'03_C_ICA_review_report.json').read_text())
    assert payload['manual_exclude'] is None
    assert payload['review_interface_version'] == 3 and payload['diagnostic_scope'] == 'full_fit'
    assert payload['component_review_coverage']['complete']
    assert not payload['fit_execution']['iteration_limit_reached']
    assert payload['fit_execution']['n_samples_fitted'] == baseline['fit_samples']
    assert payload['n_components_fitted'] == baseline['components']
    if subject == 19:
        record = payload['convergence_reconciliation']
        candidate = root/'ica_convergence/sub-19'/Path(record['result']).parent.name
        result = json.loads((candidate/'result.json').read_text())
        publication = json.loads((candidate/'publication.json').read_text())
        assert publication['published'] and publication['diagnostic_trials_unchanged']
        assert result['numeric_candidate_ok'] and not result['convergence_warnings']
        assert result['original_model_sha256'] == baseline['model_sha256']
        assert result['candidate_model_sha256'] == payload['prefit_model_sha256']
        assert publication['new_model_sha256'] == payload['prefit_model_sha256']
        assert payload['optimizer_restart'] == result['optimizer_restart']
    else:
        assert payload['prefit_model_sha256'] == baseline['model_sha256']
    figures = root/f'preprocessing/sub-{subject:02}/ses-01/figures'
    for kind in ('properties','continuous','task'):
        paths = list(figures.glob(f'03_C_ICA_review_ICA*_{kind}_*.png'))
        assert len(paths) == payload['n_components_fitted']
        for path in paths:
            sidecar = json.loads(path.with_suffix('.json').read_text())
            if kind == 'properties':
                assert sidecar['scope'] == 'full_fit' and sidecar['stimulus_locked'] is False
                assert sidecar['display_samples'] == payload['fit_execution']['n_samples_fitted']
            else:
                assert sidecar['model_sha256'] == payload['prefit_model_sha256']
    for row in payload['evidence_files']:
        assert file_stamp(root/row['relative_path']) == row['stamp'], (subject,row['relative_path'])
    # Every displayed PNG must have independent local/cluster commands pinned
    # to this exact decomposition, including after a numerical replacement.
    pictured = {Path(row['relative_path']).name for row in payload['evidence_files']
                if row['relative_path'].endswith('.png')}
    assert set(payload['plot_recipes']) == pictured
    for recipe in payload['plot_recipes'].values():
        for location in ('local','cluster'):
            assert f"--expected-model {payload['prefit_model_sha256']}" in recipe['commands'][location]
    retention = audit_participant((root,planned,accepted))
    print(f'sub-{subject:02}: model, all panels, evidence stamps, convergence and retained trials verified',flush=True)
    return {**retention,'components':payload['n_components_fitted'],
        'runs':len(payload['fingerprint']['inputs']),
        'n_iterations':payload['fit_execution']['n_iterations'],
        'figures_with_model_pinned_commands':len(pictured),
        'iteration_cap':payload['fit_execution']['max_iter'],
        'optimizer_restart':payload.get('optimizer_restart'),
        'status':'ready_for_explicit_component_review_not_scientifically_accepted'}


def main():
    """Validate cohort/receipts, parallelize participant reads, publish after join."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--derivatives-root',type=Path,required=True)
    parser.add_argument('--label',required=True)
    parser.add_argument('--jobs',type=int,default=4)
    args = parser.parse_args()
    if not args.label.replace('-','').replace('_','').isalnum():
        parser.error('Use a new simple label.')
    root = args.derivatives_root.resolve()
    output = root/'reports/ica_review_updates'/args.label/'summary.json'
    if output.parent.exists():
        parser.error('Readiness output must be NEW; preserve earlier joined snapshots.')
    plan_path = root/'stage4b_adoption/current-masks-20261109-v1/plan.json'
    receipt_path = plan_path.with_name('policy_acceptance.json')
    previous_path = root/'reports/ica_review_updates/array-11197-audit-20261110-v1/summary.json'
    plan = json.loads(plan_path.read_text());receipt = json.loads(receipt_path.read_text())
    previous = json.loads(previous_path.read_text())
    assert plan['plan_sha256'] == receipt['adoption_plan_sha256']
    assert receipt['status'] == 'accepted_for_final_ica_fit'
    assert sorted(row['subject'] for row in receipt['participants']) == COHORT
    accepted = {row['subject']:row for row in receipt['participants']}
    baseline = {row['subject']:row for row in previous['participants']}
    rows = map_participants(audit,[(root,row,accepted[row['subject']],baseline[row['subject']])
        for row in plan['participants'] if row['subject'] in COHORT],args.jobs)
    assert sorted(row['subject'] for row in rows) == COHORT
    write_json(output,{'generated_utc':utc_now(),'participants':rows,'participant_count':len(rows),
        'run_count':sum(row['runs'] for row in rows),'component_count':sum(row['components'] for row in rows),
        'retained_trials':sum(row['retained_trials'] for row in rows),
        'resources':{'participant_processes':args.jobs,'numerical_threads_per_worker':1},
        'policy_acceptance_sha256':checksum(receipt_path),'plan_sha256':checksum(plan_path),
        'baseline_audit_sha256':checksum(previous_path),
        'audit_implementation_sha256':checksum(Path(__file__)),
        'verification_scope':'actual model hashes; all panel counts/sidecars and bound evidence stamps; every PNG has model-pinned local/cluster commands; exact operative trial identities/cells and native corrected-motion-onset bounds; no repeated raw-header or large HTML/PNG hash audit',
        'scientific_acceptance':False,'production_epochs_created':False,
        'status':'ready_for_explicit_component_review'},raw_root=root.parent/'data',overwrite=False)
    print('JOINED',output,flush=True)


if __name__ == '__main__':
    main()
