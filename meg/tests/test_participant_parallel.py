"""Real process-concurrency and validation-barrier checks, using temporary data.

Usage:
    python -B -m pytest -p no:cacheprovider meg/tests/test_participant_parallel.py -q

Two spawned workers must rendezvous before either can finish: a serial executor
would fail. Check deterministic ordering, independent PIDs, one BLAS/OpenMP
thread per child, restored parent environment and allocation limits. Mocked
import stages verify that a failure in ANY participant prevents all save workers
from starting. No real cohort input/output or scheduler job is used.
"""

from argparse import Namespace
import os
from pathlib import Path
import sys
import time

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code/preprocessing_basic"))
import participant_parallel as parallel
import import_ocular_component_choices as importer


def rendezvous(task):
    """Require simultaneous workers; return resource evidence in task order."""
    folder,number=task
    (folder/str(number)).write_text(str(os.getpid()))
    deadline=time.monotonic()+15
    while len(list(folder.iterdir()))<2:
        if time.monotonic()>deadline:raise RuntimeError("Workers did not execute concurrently")
        time.sleep(.02)
    return number,os.getpid(),{k:os.environ[k] for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")}


def test_actual_spawned_workers_run_concurrently_and_restore_parent(tmp_path,monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS","7")
    monkeypatch.delenv("OPENBLAS_NUM_THREADS",raising=False)
    output=parallel.map_participants(rendezvous,[(tmp_path,2),(tmp_path,1)],2)
    assert [item[0] for item in output]==[2,1]
    assert len({item[1] for item in output})==2
    assert all(item[1]!=os.getpid() and set(item[2].values())=={"1"} for item in output)
    assert os.environ["OMP_NUM_THREADS"]=="7"
    assert "OPENBLAS_NUM_THREADS" not in os.environ


def test_worker_defaults_respect_allocation_and_explicit_limits(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK","4")
    assert parallel.participant_workers(None,36)==4
    assert parallel.participant_workers(None,2)==2
    assert parallel.participant_workers(1,36)==1
    with pytest.raises(ValueError,match="exceeds"):parallel.participant_workers(5,36)
    for value in (0,-1,"bad",True):
        with pytest.raises(ValueError):parallel.participant_workers(value,36)


@pytest.mark.parametrize("phase",["preflight","recheck"])
def test_any_participant_validation_failure_prevents_every_save(monkeypatch,phase):
    """Exercise both whole-table barriers across two independent participants."""
    args=Namespace(apply=True,decisions="unused",jobs=2)
    rows=[{"subject":4,"line":2},{"subject":5,"line":3}]
    monkeypatch.setattr(importer,"parse_table",lambda _:rows)
    # The actual process/concurrency mechanism is exercised above. Execute these
    # injected stage functions in-process so their event order remains inspectable.
    monkeypatch.setattr(importer,"map_participants",lambda fn,tasks,workers:[fn(t) for t in tasks])
    events=[]
    def preflight(task):
        subject=task[1][0]["subject"];events.append(("preflight",subject))
        if phase=="preflight" and subject==5:raise RuntimeError("stale model")
        return [{"subject":subject,"line":subject,"previous":None}]
    def recheck(task):
        subject=task[1][0]["subject"];events.append(("recheck",subject))
        if phase=="recheck" and subject==5:raise RuntimeError("new manual save")
    monkeypatch.setattr(importer,"prepare_import_subject",preflight)
    monkeypatch.setattr(importer,"recheck_import_subject",recheck)
    monkeypatch.setattr(importer,"apply_import_subject",lambda task:events.append(("SAVE",task[1][0]["subject"])))
    with pytest.raises(RuntimeError):importer.run_import(args)
    assert not any(event[0]=="SAVE" for event in events)
    assert events[:2]==[("preflight",4),("preflight",5)]
