#!/usr/bin/env python3
"""Resource-bounded process workers for independent participant operations.

Usage from a cohort script (not a standalone scientific command)
---------------------------------------------------------------
    add_parallel_argument(parser)
    workers = participant_workers(args.jobs, len(tasks))
    results = map_participants(top_level_worker, tasks, workers)

Inputs: a picklable top-level worker and participant-specific task values.
Outputs: results in original task order, after every submitted task completes.
Workers use separate spawned Python processes so MNE/Matplotlib state is not
shared. Only disjoint participant outputs may be written by these workers.
Callers retain shared publication and scientific validation barriers in the
parent. A failed validation phase raises before a later write phase can begin.

Default concurrency uses the allocated SLURM CPUs, or local CPU availability,
capped by the number of participants. Explicit --jobs may reduce concurrency;
it cannot exceed a known SLURM CPU allocation. BLAS/OpenMP threads are limited
to one per child to prevent nested oversubscription. No scheduler job or agent
is created here; processes run inside the user's existing compute allocation.
"""

from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import multiprocessing
import os


# -----------------------------------------------------------------------------
# Resource accounting is explicit; more threads do not replace participant jobs
# -----------------------------------------------------------------------------

def positive_jobs(value):
    """Accept positive integer worker counts, including CLI/environment strings."""
    text = str(value)
    if not text.isdigit() or int(text) < 1:
        raise ValueError("Participant workers must be a positive integer.")
    return int(text)


def add_parallel_argument(parser):
    """Expose the same concurrency control only on independently parallel tools."""
    parser.add_argument("--jobs", type=positive_jobs,
        help="Participant workers; default allocated SLURM CPUs or local CPUs, capped by participants.")


def participant_workers(requested, count):
    """Keep all independent participants eligible while respecting allocated CPUs."""
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    allocation = positive_jobs(allocated) if allocated else None
    available = allocation or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()) or 1
    requested = positive_jobs(requested) if requested is not None else available
    if allocation is not None and requested > allocation:
        raise ValueError(f"--jobs={requested} exceeds the {allocation} allocated SLURM CPUs.")
    return max(1, min(requested, count))


@contextmanager
def child_thread_limits():
    """Set limits before spawning/importing numerical libraries; restore the parent."""
    keys = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
    previous = {key: os.environ.get(key) for key in keys}
    try:
        for key in keys:
            os.environ[key] = "1"
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


# -----------------------------------------------------------------------------
# Independent workers finish before their caller enters a shared write barrier
# -----------------------------------------------------------------------------

def map_participants(worker, tasks, workers):
    """Run disjoint participant tasks concurrently and return deterministic order."""
    tasks = list(tasks)
    if not tasks:
        return []
    if workers == 1:
        return [worker(task) for task in tasks]
    with child_thread_limits():
        with ProcessPoolExecutor(max_workers=min(workers, len(tasks)),
                mp_context=multiprocessing.get_context("spawn")) as pool:
            return list(pool.map(worker, tasks, chunksize=1))
