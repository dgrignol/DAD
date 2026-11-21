#!/usr/bin/env python3
"""Prepare, validate and submit the approved 30-participant MNN expansion.

Usage on DIPSEN, from an isolated deployment of this repository:
  python meg/cluster/run_drsa_mnn_cohort.py prepare --run-dir /absolute/run
  python meg/cluster/run_drsa_mnn_cohort.py preflight --run-dir /absolute/run --workers 4
  python meg/cluster/run_drsa_mnn_cohort.py submit --run-dir /absolute/run --max-concurrent 30

Inputs are the accepted noica-20261117-v2 epochs/exports and their original
contracts. Per-subject selectors reproduce the saved epoch/filter settings,
including existing recovery decisions; they do not regenerate any input.
Outputs use noica-20261117-v2-drsa-mnn-v1. The exact user-approved cohort is
frozen below: 30 participant report folders and one group report folder.

The participant array runs all five ROIs and all nine valid model sources.
Subject summaries run in a second array; group publication joins that array,
then full matrix/line statistics run. Thread/resource controls are explicit;
array tasks have disjoint outputs. No legacy deployment or result is changed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import yaml

COHORT = '04 05 06 07 08 09 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 31 32 33 36 37 38 39'.split()
GROUPS = ['frontal','temporal','parietal','occipital','wholebrain']
MODELS = ['observed','predicted','prediction_error','pe_prior_unconstrained_50','pe_prior_constrained_50',
          'pe_precision_weighted','pe_kalman_innovation','pe_latent_curvature','pe_belief_update']
SOURCE_LABEL = 'noica-20261117-v2'
OUTPUT_LABEL = SOURCE_LABEL+'-drsa-mnn-v1'
DERIVATIVES = Path('/mnt/storage/tier2/morwur/Projects/DAMIANO/Attention/derivatives')
CODE = Path(__file__).resolve().parents[2]
PYTHON = Path(sys.executable).resolve()


def environment():
    """Select immutable source run and isolated analysis outputs consistently."""
    env=os.environ.copy()
    cache=Path(os.environ.get('TMPDIR','/tmp')) / f'dad_mnn_{os.environ.get("SLURM_JOB_ID","local")}_{os.environ.get("SLURM_ARRAY_TASK_ID","group")}_{os.getpid()}'
    for folder in ['matlab','matplotlib','cache']:(cache/folder).mkdir(parents=True,exist_ok=True)
    env.update(MATLAB_PREFDIR=str(cache/'matlab'),MPLCONFIGDIR=str(cache/'matplotlib'),XDG_CACHE_HOME=str(cache/'cache'))
    env.update(DAD_EPOCH_PREPARATION='inline_no_ica', DAD_PROCESSING_RUN=SOURCE_LABEL,
               DAD_DRSA_NORMALIZATION='mnn', DAD_DRSA_ANALYSIS_LABEL=OUTPUT_LABEL,
               DAD_MATLAB_BIN='/opt/ohpc/pub/apps/matlab/R2020a/bin/matlab',
               DAD_MATLAB_MAX_THREADS=os.environ.get('SLURM_CPUS_PER_TASK','4'),
               PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        env[key]='1'
    return env


def command(stage, run, subject='04'):
    """Construct an argument vector without shell interpolation of content."""
    names={'07':'07_run_matlab_drsa.py','08':'08_plot_drsa_line_summaries.py',
           '09':'09_plot_drsa_sliding_lag_surfaces.py','10':'10_drsa_stats.py'}
    return [str(PYTHON),str(CODE/'meg/code/preprocessing_basic'/names[stage]),
            '--config',str(CODE/'meg/params/pipeline_default.yaml'),'--site-config',str(run/'site.yaml'),
            '--subject-config',str(run/f'configs/sub-{subject}.yaml'),'--neural-normalization','mnn']


def execute(cmd, *, log=None):
    """Run one dependent operation and preserve its complete output on failure."""
    print('RUN',shlex.join(cmd),flush=True)
    if log:
        with Path(log).open('w') as stream:
            subprocess.run(cmd,env=environment(),stdout=stream,stderr=subprocess.STDOUT,check=True)
    else:
        subprocess.run(cmd,env=environment(),check=True)


def stage07(run, subject):
    return command('07',run,subject)+['--source-manifest',str(DERIVATIVES/'processing_runs'/SOURCE_LABEL/f'reports/sub-{subject}/ses-01/source_manifest.json'),
           '--analyses','alt','--sensor-groups',*GROUPS,'--model-sources',*MODELS,
           '--smoothing-plots','--plot-smoothing-window','9']


def prepare(run):
    """Read the whole cohort first, then write disjoint analysis selectors."""
    payloads={}
    for subject in COHORT:
        original=yaml.safe_load((CODE/f'meg/params/sub-{subject}.yaml').read_text())
        report=json.loads((DERIVATIVES/'processing_runs'/SOURCE_LABEL/f'reports/sub-{subject}/ses-01/04_epoching_report.json').read_text())
        if report['status']!='completed':raise RuntimeError(f'Incomplete accepted epochs: {subject}')
        definition=report['preparation_contract']['definition']
        original['epoching']=definition['epoching']
        original['preprocessing']=definition['preparation']['continuous_filter']
        payloads[subject]=original
    (run/'configs').mkdir(parents=True,exist_ok=True);(run/'logs').mkdir(exist_ok=True)
    for subject,payload in payloads.items():
        (run/f'configs/sub-{subject}.yaml').write_text(yaml.safe_dump(payload,sort_keys=False))
    site=yaml.safe_load((CODE/'meg/params/pipeline_cluster.yaml').read_text())
    site['project']['root']=str(CODE/'meg')
    (run/'site.yaml').write_text(yaml.safe_dump(site,sort_keys=False))
    (run/'cohort.json').write_text(json.dumps({'subjects':COHORT,'groups':GROUPS,'models':MODELS,
              'source_label':SOURCE_LABEL,'output_label':OUTPUT_LABEL,'n_expected_results_per_subject':55},indent=2)+'\n')


def preflight(run, workers):
    """Validate all accepted inputs in parallel without running MATLAB."""
    def one(subject):
        try:
            execute(stage07(run,subject)+['--validate-only'],log=run/f'logs/preflight-{subject}.log')
            print(f'PREFLIGHT {subject}: passed',flush=True);return {'subject':subject,'passed':True}
        except subprocess.CalledProcessError:
            print(f'PREFLIGHT {subject}: FAILED',flush=True);return {'subject':subject,'passed':False}
    with ThreadPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(one,COHORT))
    (run/'preflight.json').write_text(json.dumps(rows,indent=2)+'\n')
    if not all(x['passed'] for x in rows):raise RuntimeError('Cohort preflight failed; no jobs submitted.')


def verify_join():
    """Require exactly the frozen cohort and all 55 model/ROI matrices each."""
    root=DERIVATIVES/'processing_runs'/OUTPUT_LABEL
    expected={'sub-'+s for s in COHORT}
    actual={p.name for p in (root/'reports').glob('sub-*')}
    if actual!=expected:raise RuntimeError('Output cohort differs from the approved 30 participants.')
    from scipy.io import loadmat
    for subject in COHORT:
        report=json.loads((root/f'reports/sub-{subject}/ses-01/07_matlab_drsa_report.json').read_text())
        if report['status']!='completed':raise RuntimeError('Participant not completed: '+subject)
        mats=list((root/f'analysis/sub-{subject}/ses-01/drsa/matlab_real_neural/A_dRSA_alt').glob('**/results/*.mat'))
        if len(mats)!=55:raise RuntimeError(f'{subject}: expected 55 results, found {len(mats)}')
        for path in mats:
            info=loadmat(path,variable_names=['normalization_info'],squeeze_me=True,struct_as_record=False)['normalization_info']
            if info.mode!='mnn':raise RuntimeError('Mixed normalization at '+str(path))


def submit(run, concurrent, cpus, memory, hours):
    """Submit cohort arrays together and order only shared final publication."""
    rows=json.loads((run/'preflight.json').read_text())
    if [x['subject'] for x in rows]!=COHORT or not all(x['passed'] for x in rows):
        raise RuntimeError('A successful whole-cohort preflight is required.')
    if (run/'submission.json').exists():
        raise RuntimeError('Submission record already exists; inspect existing jobs before resubmission.')
    array=','.join(str(int(s)) for s in COHORT)+f'%{concurrent}'
    script=Path(__file__).resolve();jobs={}
    def job(name,action,dependency=None,is_array=False):
        cmd=['sbatch','--parsable','--job-name='+name,'--chdir='+str(run),
             '--output='+str(run/'logs'/f'{name}_%A_%a.log'),'--error='+str(run/'logs'/f'{name}_%A_%a.err'),
             '--cpus-per-task='+str(cpus),'--mem='+memory,'--time='+str(hours)+':00:00']
        if dependency:cmd+=['--dependency=afterok:'+dependency]
        if is_array:cmd+=['--array='+array]
        wrapped=shlex.join([str(PYTHON),str(script),action,'--run-dir',str(run)])
        if is_array:wrapped+=' --subject "$SLURM_ARRAY_TASK_ID"'
        cmd+=['--wrap',wrapped]
        result=subprocess.check_output(cmd,text=True).strip().split(';')[0]
        jobs[name]=result
        (run/'submission.json').write_text(json.dumps({'jobs':jobs,'subjects':COHORT,'state':'submitting'},indent=2)+'\n')
        print('SUBMITTED',name,result,flush=True);return result
    jobs['participant']=job('DAD_mnn_drsa','participant',is_array=True)
    jobs['subject_reports']=job('DAD_mnn_reports','subject-reports',jobs['participant'],True)
    jobs['group']=job('DAD_mnn_group','group',jobs['subject_reports'])
    jobs['statistics']=job('DAD_mnn_stats','stats',jobs['group'])
    (run/'submission.json').write_text(json.dumps({'jobs':jobs,'subjects':COHORT,'max_concurrent':concurrent,
             'cpus_per_task':cpus,'memory':memory,'hours':hours,'code_root':str(CODE)},indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','preflight','submit','participant','subject-reports','group','stats'])
    parser.add_argument('--run-dir',type=Path,required=True);parser.add_argument('--subject')
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--max-concurrent',type=int,default=30)
    parser.add_argument('--cpus',type=int,default=4);parser.add_argument('--memory',default='24G');parser.add_argument('--hours',type=int,default=24)
    args=parser.parse_args();run=args.run_dir.resolve()
    if min(args.workers,args.max_concurrent,args.cpus,args.hours)<1:parser.error('Resource controls must be positive.')
    subject=f'{int(args.subject):02d}' if args.subject else '04'
    if subject not in COHORT:parser.error('Subject is outside the approved cohort.')
    if args.action=='prepare':prepare(run)
    elif args.action=='preflight':preflight(run,args.workers)
    elif args.action=='submit':submit(run,args.max_concurrent,args.cpus,args.memory,args.hours)
    elif args.action=='participant':execute(stage07(run,subject)+['--resume'])
    elif args.action=='subject-reports':
        execute(command('08',run,subject)+['--subjects-only','--include-smoothed','--resume'])
        for group in GROUPS:
            execute(command('09',run,subject)+['--subjects-only','--sensor-group',group,'--resume'])
    elif args.action=='group':
        verify_join()
        execute(command('07',run)+['--group-average-only','--group-subjects',*COHORT,'--smoothing-plots','--resume'])
        execute(command('08',run)+['--subjects',*COHORT,'--group-only','--include-smoothed','--resume'])
        for group in GROUPS:
            execute(command('09',run)+['--subjects',*COHORT,'--group-only','--sensor-group',group,'--resume'])
    elif args.action=='stats':
        verify_join()
        execute(command('10',run)+['--subjects',*COHORT,'--analysis-families','A_dRSA_alt','--sensor-groups',*GROUPS,
                '--model-sources',*MODELS,'--resume'])

if __name__=='__main__':main()
