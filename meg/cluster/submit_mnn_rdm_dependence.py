#!/usr/bin/env python3
"""Submit expanded RDM dependence diagnostics after the MNN producer array.

Example on DIPSEN from an isolated, frozen diagnostic code deployment:
  python meg/cluster/submit_mnn_rdm_dependence.py \
    --run-dir /home/USER/runs/mnn-dependence \
    --production-run-dir /home/USER/runs/mnn-production \
    --brain-results-root /derivatives/processing_runs/SOURCE-drsa-mnn-v1 \
    --diagnostics-root /derivatives/processing_runs/SOURCE-drsa-mnn-v1/rdm_dependence_mnn_DATE \
    --dependency 13012

This submits a read-only whole-cohort preflight, a 30-participant array, then a
joined group render. Dependencies are producer -> preflight -> array -> group.
Only fresh diagnostic outputs and this run's logs/launch records are written.
The original producer snapshot/jobs and existing reports are never modified.
Worker resource controls are explicit. Shared rendering follows the full join.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

COHORT = '04 05 06 07 08 09 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 31 32 33 36 37 38 39'.split()


def worker(run, role):
    """Run one scheduler action with private caches and bounded threads."""
    launch=json.loads((run/'launch.json').read_text())
    env=os.environ.copy()
    cache=Path(os.environ.get('TMPDIR','/tmp'))/f'dad_mnn_dependence_{os.environ.get("SLURM_JOB_ID","local")}_{os.getpid()}'
    for folder in ['matlab','matplotlib','cache']:(cache/folder).mkdir(parents=True,exist_ok=True)
    env.update(MATLAB_PREFDIR=str(cache/'matlab'),MPLCONFIGDIR=str(cache/'matplotlib'),
               XDG_CACHE_HOME=str(cache/'cache'),PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
    cmd=launch['common']
    if role=='preflight':
        cmd+=['--mode','validate','--workers',str(launch['preflight_workers']),'--threads','1','--subjects',*COHORT]
    elif role=='participant':
        subject=f'{int(os.environ["SLURM_ARRAY_TASK_ID"]):02d}'
        if subject not in COHORT:raise ValueError('Participant outside approved cohort.')
        cmd+=['--mode','participants','--subjects',subject,'--workers','1','--threads',str(launch['participant_cpus'])]
    else:
        cmd+=['--mode','group','--workers',str(launch['group_workers']),'--threads','1','--subjects',*COHORT]
    print(shlex.join(cmd),flush=True)
    subprocess.run(cmd,env=env,check=True)


def submit(args):
    """Reserve one fresh output namespace and persist every submitted job ID."""
    run=args.run_dir.resolve();brain=args.brain_results_root.resolve(strict=True);output=args.diagnostics_root.resolve()
    if brain not in output.parents or any(p in (brain/'analysis',brain/'reports') for p in (output,*output.parents)):
        raise ValueError('Diagnostics must be a separate child of the brain-result root.')
    if (run/'submission.json').exists() or output.exists():
        raise FileExistsError('A submission or diagnostic folder already exists; use a fresh diagnostic root.')
    runner=Path(__file__).resolve().parents[1]/'code/preprocessing_basic/run_drsa_dependence_cohort.py'
    production=args.production_run_dir.resolve(strict=True)
    common=[sys.executable,str(runner),'--source-run','noica-20261117-v2',
            '--brain-results-root',str(brain),'--diagnostics-root',str(output),
            '--subject-config-dir',str(production/'configs'),'--site-config',str(production/'site.yaml'),
            '--matlab-bin','/opt/ohpc/pub/apps/matlab/R2020a/bin/matlab']
    launch={'common':common,'subjects':COHORT,'preflight_workers':args.preflight_workers,
            'group_workers':args.group_workers,'participant_cpus':args.participant_cpus,
            'producer_dependency':args.dependency,'maximum_concurrency':args.max_concurrent}
    run.mkdir(parents=True,exist_ok=True);(run/'logs').mkdir(exist_ok=True);output.mkdir()
    (run/'launch.json').write_text(json.dumps(launch,indent=2)+'\n')
    jobs={}
    def job(role,dependency,cpus,memory,array=False):
        cmd=['sbatch','--parsable','--job-name=DAD_mnn_rdm_'+role,'--chdir='+str(run),
             '--output='+str(run/'logs'/f'{role}_%A_%a.log'),'--error='+str(run/'logs'/f'{role}_%A_%a.err'),
             '--cpus-per-task='+str(cpus),'--mem='+memory,'--time='+str(args.hours)+':00:00',
             '--dependency=afterok:'+dependency]
        if array:cmd+=['--array='+','.join(str(int(s)) for s in COHORT)+f'%{args.max_concurrent}']
        cmd+=['--wrap',shlex.join([sys.executable,str(Path(__file__).resolve()),'--run-dir',str(run),'--worker',role])]
        identity=subprocess.check_output(cmd,text=True).strip().split(';')[0];jobs[role]=identity
        record={'jobs':jobs,'producer':args.dependency,'subjects':COHORT,'output_root':str(output),'launch':launch}
        (run/'submission.json').write_text(json.dumps(record,indent=2)+'\n')
        (output/'submission.json').write_text(json.dumps(record,indent=2)+'\n')
        print('SUBMITTED',role,identity,flush=True);return identity
    first=job('preflight',args.dependency,args.preflight_workers,args.preflight_memory)
    second=job('participant',first,args.participant_cpus,args.participant_memory,True)
    job('group',second,args.group_workers,args.group_memory)
    return jobs


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--worker',choices=['preflight','participant','group'])
    parser.add_argument('--production-run-dir',type=Path)
    parser.add_argument('--brain-results-root',type=Path)
    parser.add_argument('--diagnostics-root',type=Path)
    parser.add_argument('--dependency')
    parser.add_argument('--preflight-workers',type=int,default=4);parser.add_argument('--group-workers',type=int,default=4)
    parser.add_argument('--participant-cpus',type=int,default=2);parser.add_argument('--max-concurrent',type=int,default=30)
    parser.add_argument('--preflight-memory',default='12G');parser.add_argument('--participant-memory',default='24G')
    parser.add_argument('--group-memory',default='16G');parser.add_argument('--hours',type=int,default=4)
    args=parser.parse_args()
    if args.worker:worker(args.run_dir.resolve(),args.worker);return
    if not all([args.production_run_dir,args.brain_results_root,args.diagnostics_root,args.dependency]):
        parser.error('Submission requires production/input/output roots and producer dependency.')
    if min(args.preflight_workers,args.group_workers,args.participant_cpus,args.max_concurrent,args.hours)<1:
        parser.error('Resource controls must be positive.')
    if not args.dependency.isdigit():parser.error('Producer dependency must be a Slurm job ID.')
    submit(args)

if __name__=='__main__':main()
