#!/usr/bin/env python3
"""Run the complete MNN RDM diagnostic cohort inside one Slurm allocation.

QoS fallback example (after the complete producer array):
  sbatch --dependency=afterok:13012 --cpus-per-task=60 --mem=192G --time=04:00:00 \
    --wrap='python meg/cluster/run_mnn_rdm_allocation.py --run-dir /absolute/run'

Read launch.json from the separately frozen diagnostic deployment. Validate
all 30 participants, run independent participant processes in parallel, then
join before group validation/rendering. This avoids consuming 30 additional
scheduler submission slots; it does not process the cohort sequentially.
Each participant uses a private cache and two numerical/MATLAB threads. The
192-GB allocation budgets the measured approximately 4-GB participant footprint
across 30 workers plus headroom. Shared plotting/publication is ordered and
uses four input/summary workers. No source or producer output is changed.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
from submit_mnn_rdm_dependence import COHORT


def run_allocation(args):
    """Execute dependent phases with explicit, bounded participant resources."""
    allocation=int(os.environ.get('SLURM_CPUS_PER_TASK',os.cpu_count() or 1))
    if min(args.participant_workers,args.participant_threads,args.validation_workers,args.group_workers)<1:
        raise ValueError('Resource controls must be positive.')
    if max(args.participant_workers*args.participant_threads,args.validation_workers,args.group_workers)>allocation:
        raise ValueError('Workers and numerical threads exceed the allocated CPUs.')
    run=args.run_dir.resolve();launch=json.loads((run/'launch.json').read_text())
    if launch['subjects']!=COHORT:raise ValueError('Launch cohort differs from the approved 30 participants.')
    phases=[('validate',args.validation_workers,1),('participants',args.participant_workers,args.participant_threads),
            ('group',args.group_workers,1)]
    progress=[]
    for mode,workers,threads in phases:
        command=launch['common']+['--mode',mode,'--workers',str(workers),'--threads',str(threads),'--subjects',*COHORT]
        print(f'MNN RDM allocation | {mode}: {workers} workers x {threads} threads',flush=True)
        subprocess.run(command,check=True)
        progress.append({'mode':mode,'status':'completed','workers':workers,'threads':threads})
        (run/'allocation_progress.json').write_text(json.dumps(progress,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--participant-workers',type=int,default=30)
    parser.add_argument('--participant-threads',type=int,default=2)
    parser.add_argument('--validation-workers',type=int,default=4)
    parser.add_argument('--group-workers',type=int,default=4)
    run_allocation(parser.parse_args())

if __name__=='__main__':main()
