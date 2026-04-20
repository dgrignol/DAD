#!/usr/bin/env python3
"""Helpers for detecting trigger events with overlapping pulses."""

from __future__ import annotations

import numpy as np
import mne
from mne.utils.config import _get_stim_channel


def _merge_samples(min_duration: float, sfreq: float) -> int:
    # Convert the minimum duration into a sample-count merge threshold
    # compatible with MNE's internal "merge" semantics.
    min_samples = min_duration * sfreq
    if min_samples > 0:
        merge = int(min_samples // 1)
        if merge == min_samples:
            merge -= 1
        return merge
    return 0


def _find_stim_steps_1d(
    data: np.ndarray,
    first_samp: int,
    *,
    pad_stop: int | None = 0,
    merge: int = 0,
) -> np.ndarray:
    # Detect step changes in a 1D stim channel, preserving sample indices.
    data = np.asarray(data)
    if data.size == 0:
        return np.empty((0, 3), dtype=np.int64)
    changed = np.diff(data) != 0
    idx = np.where(changed)[0]
    if idx.size == 0:
        steps = np.empty((0, 3), dtype=np.int64)
    else:
        pre_step = data[idx]
        idx = idx + 1
        post_step = data[idx]
        idx = idx + first_samp
        steps = np.c_[idx, pre_step, post_step]

    if pad_stop is not None:
        v = data[-1]
        if v != pad_stop:
            last_idx = len(data) + first_samp
            steps = np.append(steps, [[last_idx, v, pad_stop]], axis=0)

    if merge != 0 and steps.size:
        # Merge steps closer than the merge threshold, emulating MNE behavior.
        diff = np.diff(steps[:, 0])
        idx_merge = diff <= abs(merge)
        if np.any(idx_merge):
            where = np.where(idx_merge)[0]
            keep = np.logical_not(idx_merge)
            if merge > 0:
                steps[where + 1, 1] = steps[where, 1]
                keep = np.append(keep, True)
            else:
                steps[where, 2] = steps[where + 1, 2]
                keep = np.insert(keep, 0, True)
            is_step = steps[:, 1] != steps[:, 2]
            keep = np.logical_and(keep, is_step)
            steps = steps[keep]

    return steps


def _mask_steps(steps: np.ndarray, mask: int | None, mask_type: str) -> np.ndarray:
    # Apply bitmask filtering to the pre/post step values.
    if steps.size == 0 or mask is None:
        return steps
    if mask_type == "not_and":
        mask = np.bitwise_not(mask)
    elif mask_type != "and":
        raise ValueError(
            "'mask_type' should be either 'and' or 'not_and', "
            f"got {mask_type!r}"
        )
    steps = steps.copy()
    steps[:, 1:] = np.bitwise_and(steps[:, 1:], mask)
    steps = steps[steps[:, 1] != steps[:, 2]]
    return steps


def _as_structured(arr: np.ndarray) -> np.ndarray:
    # Use a structured view to enable fast row-wise comparisons via np.isin.
    if arr.size == 0:
        return arr
    return np.ascontiguousarray(arr).view(
        np.dtype((np.void, arr.dtype.itemsize * arr.shape[1]))
    ).ravel()


def find_events_with_overlaps(
    raw,
    stim_channel=None,
    output="onset",
    consecutive="increasing",
    min_duration=0,
    shortest_event=2,
    mask=None,
    uint_cast=False,
    mask_type="and",
    initial_event=False,
    verbose=None,
):
    """Return find_events output plus overlap steps (non-zero to non-zero)."""
    # Start with standard MNE event detection.
    events = mne.find_events(
        raw,
        stim_channel=stim_channel,
        output=output,
        consecutive=consecutive,
        min_duration=min_duration,
        shortest_event=shortest_event,
        mask=mask,
        uint_cast=uint_cast,
        mask_type=mask_type,
        initial_event=initial_event,
        verbose=verbose,
    )
    if output != "step":
        return events

    # For step output, augment with overlapping non-zero steps.
    stim_channels = _get_stim_channel(stim_channel, raw.info)
    picks = mne.pick_channels(raw.info["ch_names"], include=stim_channels, ordered=False)
    if len(picks) == 0:
        return events

    data = raw.get_data(picks=picks)
    merge = _merge_samples(min_duration, float(raw.info["sfreq"]))
    extra_steps = []

    for channel_data in data:
        if channel_data.size == 0:
            continue
        # Normalize to signed int, optionally unsigned-cast, and remove negatives.
        data_int = channel_data.astype(np.int64, copy=False)
        if uint_cast:
            data_int = data_int.astype(np.uint16).astype(np.int64)
        if data_int.min() < 0:
            data_int = np.abs(data_int)
        steps = _find_stim_steps_1d(
            data_int, raw.first_samp, pad_stop=0, merge=merge
        )
        steps = _mask_steps(steps, mask, mask_type)
        if steps.size == 0:
            continue
        # Keep only transitions where both the pre- and post- values are non-zero.
        non_zero = (steps[:, 1] != 0) & (steps[:, 2] != 0)
        steps = steps[non_zero]
        if steps.size:
            extra_steps.append(steps)

    if not extra_steps:
        return events

    # Remove steps already present in the base event set.
    extra_steps = np.vstack(extra_steps)
    if events.size:
        existing = _as_structured(events)
        candidate = _as_structured(extra_steps)
        missing_mask = ~np.isin(candidate, existing)
        extra_steps = extra_steps[missing_mask]

    if extra_steps.size == 0:
        return events

    # Combine, unique, and sort by sample index for a stable event table.
    combined = (
        np.vstack([events, extra_steps]) if events.size else extra_steps.copy()
    )
    combined = np.unique(combined, axis=0)
    combined = combined[np.argsort(combined[:, 0])]
    return combined.astype(events.dtype, copy=False)
