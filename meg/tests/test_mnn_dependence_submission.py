"""Test the producer/preflight/full-array/group diagnostic submission order.

Usage: python -m pytest meg/tests/test_mnn_dependence_submission.py -q
All scheduler calls are mocked; no cluster jobs or production writes occur.
"""
import argparse
import importlib.util
import json
from pathlib import Path


def test_submit_joins_complete_cohort_after_producer(tmp_path,monkeypatch):
    path=Path(__file__).resolve().parents[1]/'cluster/submit_mnn_rdm_dependence.py'
    spec=importlib.util.spec_from_file_location('submit_mnn_dependence',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    calls=[]
    def scheduler(cmd,**kwargs):
        calls.append(cmd);return str(len(calls))+'\n'
    monkeypatch.setattr(module.subprocess,'check_output',scheduler)
    brain=tmp_path/'mnn';brain.mkdir();production=tmp_path/'production';production.mkdir()
    args=argparse.Namespace(run_dir=tmp_path/'run',brain_results_root=brain,
        diagnostics_root=brain/'dependence',production_run_dir=production,dependency='13012',
        preflight_workers=4,group_workers=4,participant_cpus=2,max_concurrent=30,
        preflight_memory='12G',participant_memory='24G',group_memory='16G',hours=4)
    assert module.submit(args)=={'preflight':'1','participant':'2','group':'3'}
    assert '--dependency=afterok:13012' in calls[0]
    assert '--dependency=afterok:1' in calls[1]
    assert '--dependency=afterok:2' in calls[2]
    array=next(x for x in calls[1] if x.startswith('--array='))
    assert array=='--array='+','.join(str(int(s)) for s in module.COHORT)+'%30'
    assert len(module.COHORT)==30
    launch=json.loads((args.run_dir/'launch.json').read_text())
    assert '--subject-config-dir' in launch['common']
    assert not any(x.startswith('--array=') for x in calls[2])
