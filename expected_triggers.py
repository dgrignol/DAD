#!/usr/bin/env python3
"""
Compute the expected trigger timeline for MoveDot1.

This script reconstructs the deterministic, experiment-controlled trigger
sequence for each trial/block using the same ordering and timing rules as
`experiment/MoveDot1_experiment_vX.m` and the trigger IDs documented in
`experiment/trigger_codes.md`. It does not include subject-driven events
(responses, gaze breaks, replays) and ignores any replay/abort delays.

Inputs
- experiment/input_files/MovDot_SubXX.mat (stimulus sequences + fps)
- experiment/input_files/SubXX_TrialStruct.mat (BlockOrder, TrialOrder, TrialStruct)

Outputs
- One CSV per block with columns: trigger,time_s
- Files are written under derivatives/triggers/subXX by default

Run selection
- If --run is provided, generate only that run.
- If --run is omitted, generate all runs available in TrialStruct.
  Each run reuses the same RNG seed as a standalone invocation so outputs
  match per-run calls.

Usage examples
  python expected_triggers.py --subject 4 --run 1
  python expected_triggers.py --subject 4
  python expected_triggers.py --subject 4 --seed 123
  python expected_triggers.py --subject 4 --input-dir /path/to/input_files
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
from scipy.io import loadmat


def build_cond_matrix(xy_seqs: np.ndarray) -> np.ndarray:
    """Recreate the 4 x nSeq condMatrix used in the MATLAB code."""
    n_seq = xy_seqs.size
    cond_matrix = np.zeros((4, n_seq), dtype=float)
    # xy_seqs(:) in MATLAB uses column-major order; mirror that with order='F'
    # so the condMatrix aligns with MoveDot1_experiment_vX.m.
    for idx, obj in enumerate(xy_seqs.ravel(order="F")):
        cond_matrix[0, idx] = obj.condition
        cond_matrix[1, idx] = obj.sequence
        cond_matrix[2, idx] = idx + 1  # MATLAB 1-based iSeq
        cond_matrix[3, idx] = obj.PredictionRange
    return cond_matrix


def block_events(
    block_cond: int,
    trial_order: np.ndarray,
    cond_matrix: np.ndarray,
    trial_struct: np.ndarray,
    fps: int,
    n_frames: int,
    rng: np.random.Generator,
) -> List[Tuple[float, int]]:
    """
    Build (time_s, trigger) events for one block.

    Times are relative to the block start (StartTimeX in the MATLAB code).
    """
    # ITI jitter matches the experiment: 0.5 + 0.04 * rand.
    iti = 0.5 + 0.04 * rng.random(len(trial_order))  # matches MATLAB rand usage
    # TrialOrder is 1-based in MATLAB; convert to 0-based indexing for NumPy.
    cond_shuffled = cond_matrix[:, trial_order - 1]  # convert order to 0-based

    # Accumulate time from block start to each trigger event.
    t = 0.0  # seconds from block start
    events: List[Tuple[float, int]] = []

    for trial_idx in range(cond_shuffled.shape[1]):
        i_seq = int(cond_shuffled[2, trial_idx])  # original iSeq (1-based)
        ts = trial_struct[i_seq - 1]  # trial info for this sequence

        # Trial-level trigger for first frame, following the experiment's mapping
        # in experiment/trigger_codes.md and MoveDot1_experiment_vX.m:
        # - Non-catch trials use sequence + condition/block offsets.
        # - Catch trials always use 102.
        num_catch = np.atleast_1d(ts.Start).size
        if num_catch > 0:
            stim_trigger = 102
        else:
            stim_trigger = cond_shuffled[1, trial_idx]
            if cond_shuffled[0, trial_idx] == 45:
                stim_trigger += 20
            if block_cond == 2:
                stim_trigger += 40

        # Apply ITI, then emit first-frame trigger at the new trial onset.
        t += iti[trial_idx]
        events.append((t, int(stim_trigger)))

        # Catch triggers (if any): each catch adds a 100 (start) and 101 (end).
        if num_catch > 0:
            starts = np.atleast_1d(ts.Start)
            ends = np.atleast_1d(ts.End)
            for start_frame, end_frame in zip(starts, ends):
                events.append((t + (start_frame - 1) / fps, 100))
                events.append((t + (end_frame - 1) / fps, 101))
            last_trigger = 103
        else:
            last_trigger = 81

        # Last frame of the video: 81 for no-catch, 103 for catch trials.
        events.append((t + (n_frames - 1) / fps, last_trigger))

        # Advance to next trial (full video duration).
        t += n_frames / fps

    return events


def write_csv(path: Path, events: Iterable[Tuple[float, int]]) -> None:
    """Write two-column CSV: trigger,time_s."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["trigger", "time_s"])
        # Sort by time to keep rows ordered even if events were added out of order.
        for time_s, trigger in sorted(((t, tr) for t, tr in events)):
            writer.writerow([trigger, f"{time_s:.6f}"])


def resolve_runs(run_arg: int | None, trial_struct: np.ndarray) -> List[int]:
    """
    Return 1-based run numbers to generate.

    If run_arg is None, all runs in TrialStruct are returned.
    """
    trial_struct_arr = np.atleast_1d(trial_struct)
    num_runs = trial_struct_arr.size
    if run_arg is None:
        return list(range(1, num_runs + 1))
    return [run_arg]


def main() -> None:
    # Parse CLI arguments to determine subject, run selection, and I/O roots.
    parser = argparse.ArgumentParser(
        description="Construct expected trigger timeline for MoveDot1."
    )
    parser.add_argument(
        "--subject", "-s", type=int, default=99, help="Subject number (e.g., 99)."
    )
    parser.add_argument(
        "--run",
        "-r",
        type=int,
        default=None,
        help=(
            "Run number (1-based). If omitted, generate all runs "
            "available for the subject."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="RNG seed for ITI jitter (default: subject number).",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("experiment/input_files"),
        help="Directory containing MovDot_SubXX.mat and SubXX_TrialStruct.mat.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Where to write CSVs (one per block). Defaults to "
            "derivatives/triggers/subXX."
        ),
    )
    args = parser.parse_args()

    # Resolve output location for subject-specific CSVs.
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = Path("derivatives/triggers") / f"sub{args.subject:02d}"

    # Resolve file paths and load stimulus + trial-structure MATLAB files.
    stim_path = args.input_dir / f"MovDot_Sub{args.subject:02d}.mat"
    trial_path = args.input_dir / f"Sub{args.subject:02d}_TrialStruct.mat"

    stim_mat = loadmat(stim_path, squeeze_me=True, struct_as_record=False)
    trial_mat = loadmat(trial_path, squeeze_me=True, struct_as_record=False)

    # Extract core timing/ordering inputs shared across runs.
    cfg = stim_mat["Cfg"]
    xy_seqs = stim_mat["xySeqs"]
    fps = int(cfg.fps)
    n_frames = xy_seqs.flat[0].xy.shape[0]

    # Trial structure arrays follow the MATLAB layout for this subject.
    block_order = trial_mat["BlockOrder"]
    trial_order = trial_mat["TrialOrder"]
    trial_struct_all = trial_mat["TrialStruct"]

    # Match experiment ITI randomness: default seed is subject ID.
    seed_value = args.seed if args.seed is not None else args.subject
    cond_matrix = build_cond_matrix(xy_seqs)
    runs = resolve_runs(args.run, trial_struct_all)

    # Emit a CSV per block, per run (BlockOrder controls which dot is attended).
    num_blocks = block_order.shape[1]
    for run in runs:
        run_idx = run - 1  # convert to 0-based
        trial_struct = np.atleast_1d(trial_struct_all)[run_idx]
        # Re-seed per run so outputs match standalone per-run invocations.
        rng = np.random.default_rng(seed_value)

        for block_idx in range(num_blocks):
            block_cond = int(block_order[run_idx, block_idx])
            block_events_list = block_events(
                block_cond,
                trial_order[block_idx, :, run_idx],
                cond_matrix,
                trial_struct,
                fps,
                n_frames,
                rng,
            )

            out_path = output_dir / (
                f"expected_triggers_sub{args.subject:02d}"
                f"_run{run:02d}_block{block_idx + 1}.csv"
            )
            write_csv(out_path, block_events_list)
            print(
                f"Wrote {len(block_events_list)} events for run {run} "
                f"block {block_idx + 1} ({out_path})"
            )


if __name__ == "__main__":
    main()
