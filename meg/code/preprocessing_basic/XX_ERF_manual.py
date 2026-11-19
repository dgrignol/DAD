#!/usr/bin/env python3
"""
XX_ERF_manual.py

Purpose
Provide an optional interactive ERF inspection workflow that uses the same
resolved subject, epoch path, and condition definitions as production
epoching. This script is not a required scientific stage and writes no files.

Usage examples
1) Start IPython with an interactive matplotlib backend:

   cd /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg
   .venv/bin/ipython --matplotlib=qt

2) Inspect subject 01 with the default condition:

   %run -i code/preprocessing_basic/XX_ERF_manual.py -- \
     --config params/pipeline_default.yaml \
     --subject-config params/sub-01.yaml

3) Select another condition before re-running:

   condition = "occluded_deviant"
   %run -i code/preprocessing_basic/XX_ERF_manual.py -- \
     --config params/pipeline_default.yaml \
     --subject-config params/sub-04.yaml

Inputs
- The configured step-04 epochs FIF.
- The resolved subject profile and epoch event map.

Outputs
- In-memory ``epochs``, ``epochs_condition``, ``evoked``, ``evokeds`` objects.
- Interactive trial-image, sensor-layout, topomap, and condition-comparison
  figures.

Scientific assumptions
- The default ERF baseline is -0.1 to 0 s. It lies inside the configured
  -0.5 s epoch and is compatible with every verified ITI.
- Step 04 performs no baseline correction; baseline is applied here only to
  copied evoked/inspection objects.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import mne

from pipeline_config import load_and_resolve_config


# -----------------------------------------------------------------------------
# Configuration and shared epoch path resolution
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the optional inspection parser."""
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Interactive DAD ERF inspection.")
    parser.add_argument(
        "--config",
        default=str(project_root / "params" / "pipeline_default.yaml"),
    )
    parser.add_argument("--site-config", default=None)
    parser.add_argument(
        "--subject-config",
        default=str(project_root / "params" / "sub-04.yaml"),
    )
    parser.add_argument("--tuning-config", default=None)
    return parser


args = build_arg_parser().parse_args()
cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
    args.config,
    args.subject_config,
    args.tuning_config,
    site_config_path=args.site_config,
    required_sections=("epoching",),
)
subject = cfg["layout"]["subject_bids_id"]
session = cfg["layout"]["session_bids_id"]
epochs_fpath = (
    Path(cfg["resolved_paths"]["derivatives_analysis_dir"])
    / "meg"
    / f"{subject}_{session}_task-Attention_run-all_epo.fif"
)
if not epochs_fpath.is_file():
    raise FileNotFoundError(f"Configured epochs FIF not found: {epochs_fpath}")


# -----------------------------------------------------------------------------
# Interactive data loading and user-editable choices
# -----------------------------------------------------------------------------


plt.ion()
if globals().get("loaded_epochs_fpath") == str(epochs_fpath) and "epochs" in globals():
    print(f"Reusing already-loaded epochs: {epochs_fpath}")
else:
    epochs = mne.read_epochs(epochs_fpath, proj=False, preload=True, verbose=True)
    loaded_epochs_fpath = str(epochs_fpath)

condition_names = list(cfg["epoching"]["event_id_map"])
condition = globals().get("condition", condition_names[0])
l_freq = globals().get("l_freq", 0.0)
h_freq = globals().get("h_freq", 30.0)
crop_tmin = globals().get("crop_tmin", -0.1)
crop_tmax = globals().get("crop_tmax", 0.4)
baseline = globals().get("baseline", (-0.1, 0.0))
topomap_time = globals().get("topomap_time", 0.110)
image_pick = globals().get("image_pick", "MEG0111")
compare_pick = globals().get("compare_pick", "MEG0111")
image_vmin = globals().get("image_vmin", None)
image_vmax = globals().get("image_vmax", None)

if condition not in epochs.event_id:
    raise ValueError(
        f"Condition {condition!r} is absent; configured conditions are "
        f"{condition_names}."
    )
if baseline[0] < epochs.tmin or baseline[1] > epochs.tmax:
    raise ValueError(
        f"ERF baseline {baseline} lies outside epoch range "
        f"[{epochs.tmin}, {epochs.tmax}]."
    )

print(f"Loaded: {epochs_fpath}")
print(f"Configuration warnings: {warnings}")
print(f"Conditions: {epochs.event_id}")
print(f"Selected condition: {condition}")
print(f"Epoch range: {epochs.tmin:.3f}..{epochs.tmax:.3f} s")
print(f"Inspection baseline: {baseline} s")


# -----------------------------------------------------------------------------
# Condition average and interactive diagnostic figures
# -----------------------------------------------------------------------------


epochs_condition = epochs[condition].copy()
evoked = (
    epochs_condition.average(method="mean")
    .filter(l_freq, h_freq)
    .crop(crop_tmin, crop_tmax)
)
evoked_baseline = evoked.copy().apply_baseline(baseline=baseline)

fig_image = (
    epochs_condition.copy()
    .filter(l_freq, h_freq)
    .crop(crop_tmin, crop_tmax)
    .plot_image(
        picks=[image_pick],
        vmin=image_vmin,
        vmax=image_vmax,
    )
)
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


# -----------------------------------------------------------------------------
# Optional comparison across all configured conditions
# -----------------------------------------------------------------------------


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
