#!/usr/bin/env python3
"""
05_ERF_manual.py

Purpose
Create a minimal manual/IPython workflow for inspecting event-related fields
(ERFs) from the DAD MEG epochs. This file intentionally behaves like an
interactive notebook: run it in IPython, edit variables, re-run cells, and keep
objects such as `epochs`, `evoked`, and `evokeds` available at the prompt.

Usage
1) Start IPython with Qt plotting from the project root:
   cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg
   ../.venv/bin/ipython --matplotlib=qt

2) Run this file into the live IPython namespace:
   %run -i code/preprocessing_basic/05_ERF_manual.py

3) Inspect variables and re-run cells after edits:
   %whos
   condition = "occluded_deviant"
   %run -i code/preprocessing_basic/05_ERF_manual.py

Inputs
- Epoched FIF from step 04:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/derivatives/analysis/sub-01/ses-01/meg/sub-01_ses-01_task-Attention_run-all_epo.fif

Outputs
- No files are written. All outputs are interactive MNE/matplotlib figures and
  in-memory variables.

Workflow adapted from FLUX EventRelatedFields.ipynb
- Read epochs.
- Select one condition.
- Average trials with `Epochs.average(method="mean")`.
- Low-pass at 30 Hz and crop to -100..400 ms.
- Plot trial image, magnetometer/gradiometer topo traces, and topomaps.

Assumptions
- `ipython`, `mne`, `matplotlib`, `mne-qt-browser`, and `pyqt6` are installed in
  the active virtual environment.
- Step 04 has already produced the epochs file.
- DAD condition names are:
  `always_visible`, `occluded_nondeviant`, and `occluded_deviant`.
"""

# %%
# -------------------------------------------------------------------------
# Imports and interactive plotting setup
# -------------------------------------------------------------------------
from pathlib import Path

import matplotlib.pyplot as plt
import mne

# Keep plots interactive in an IPython session launched with:
# ../.venv/bin/ipython --matplotlib=qt
plt.ion()
mne.viz.set_browser_backend("qt")

# %%
# -------------------------------------------------------------------------
# Load DAD epochs
# -------------------------------------------------------------------------
PROJECT_ROOT = Path("/Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg")
EPOCHS_FPATH = (
    PROJECT_ROOT
    / "derivatives"
    / "analysis"
    / "sub-01"
    / "ses-01"
    / "meg"
    / "sub-01_ses-01_task-Attention_run-all_epo.fif"
)

# The epochs file is large. In IPython, reuse it on repeated `%run -i` calls
# unless the path changes or you delete `epochs` from the namespace.
if globals().get("loaded_epochs_fpath") == str(EPOCHS_FPATH) and "epochs" in globals():
    print(f"Reusing already-loaded epochs: {EPOCHS_FPATH}")
else:
    epochs = mne.read_epochs(EPOCHS_FPATH, proj=False, preload=True, verbose=True)
    loaded_epochs_fpath = str(EPOCHS_FPATH)

print(f"Loaded: {EPOCHS_FPATH}")
print(f"Conditions: {epochs.event_id}")
print(f"Epochs: {len(epochs)}")
print(f"Time window: {epochs.tmin:.3f}..{epochs.tmax:.3f} s")
print(f"Channels: {len(epochs.ch_names)}")

# %%
# -------------------------------------------------------------------------
# Manual choices: edit these variables, then re-run cells below
# -------------------------------------------------------------------------
# These defaults respect variables already set at the IPython prompt. For
# example, run `condition = "occluded_deviant"` before `%run -i` to switch.
condition = globals().get("condition", "always_visible")

# FLUX uses -100..400 ms and a 30 Hz low-pass for cognitive ERFs.
l_freq = globals().get("l_freq", 0.0)
h_freq = globals().get("h_freq", 30.0)
crop_tmin = globals().get("crop_tmin", -0.1)
crop_tmax = globals().get("crop_tmax", 0.4)
baseline = globals().get("baseline", (-0.1, 0.0))
topomap_time = globals().get("topomap_time", 0.110)

# Real DAD channels from the loaded file.
image_pick = globals().get("image_pick", "MEG0111")
compare_pick = globals().get("compare_pick", "MEG0111")

# Set these to numbers if you want fixed FLUX-like scaling, e.g. -500, 500.
image_vmin = globals().get("image_vmin", None)
image_vmax = globals().get("image_vmax", None)

print(f"Selected condition: {condition}")
print(f"Selected image/compare channel: {image_pick}")

# %%
# -------------------------------------------------------------------------
# Average trials for one condition
# -------------------------------------------------------------------------
epochs_condition = epochs[condition].copy()

evoked = (
    epochs_condition
    .average(method="mean")
    .filter(l_freq, h_freq)
    .crop(crop_tmin, crop_tmax)
)

evoked_baseline = evoked.copy().apply_baseline(baseline=baseline)
print(evoked_baseline)

# %%
# -------------------------------------------------------------------------
# Plot single-trial image for one channel, matching the FLUX plot_image step
# -------------------------------------------------------------------------
fig_image = (
    epochs_condition
    .copy()
    .filter(l_freq, h_freq)
    .crop(crop_tmin, crop_tmax)
    .plot_image(
        picks=[image_pick],
        vmin=image_vmin,
        vmax=image_vmax,
    )
)

# %%
# -------------------------------------------------------------------------
# Plot average ERF traces over the sensor layout
# -------------------------------------------------------------------------
fig_mag_topo = evoked_baseline.copy().pick("mag").plot_topo(
    title=f"{condition} | Magnetometers"
)

fig_grad_topo = evoked_baseline.copy().pick("grad").plot_topo(
    title=f"{condition} | Gradiometers"
)

fig_grad_topo_merged = evoked_baseline.copy().pick("grad").plot_topo(
    title=f"{condition} | Gradiometers merged",
    merge_grads=True,
)

# %%
# -------------------------------------------------------------------------
# Plot topographies at one latency
# -------------------------------------------------------------------------
fig_mag_topomap = evoked_baseline.plot_topomap(
    topomap_time,
    ch_type="mag",
    time_unit="s",
)

fig_grad_topomap = evoked_baseline.plot_topomap(
    topomap_time,
    ch_type="grad",
    time_unit="s",
)

# %%
# -------------------------------------------------------------------------
# Optional: compare all DAD conditions at one channel
# -------------------------------------------------------------------------
condition_names = globals().get(
    "condition_names",
    ["always_visible", "occluded_nondeviant", "occluded_deviant"],
)

evokeds = {
    name: (
        epochs[name]
        .copy()
        .average(method="mean")
        .filter(l_freq, h_freq)
        .crop(crop_tmin, crop_tmax)
        .apply_baseline(baseline=baseline)
    )
    for name in condition_names
}

fig_compare = mne.viz.plot_compare_evokeds(
    evokeds,
    picks=compare_pick,
    combine=None,
    title=f"Condition comparison | {compare_pick}",
)
