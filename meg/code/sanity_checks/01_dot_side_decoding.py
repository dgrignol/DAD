#!/usr/bin/env python3
"""
01_dot_side_decoding.py

Purpose
Run a MEG sanity check that decodes dot horizontal side (left vs right half of
screen) at each epoch sample for one condition, and generate a 3-panel figure:
1) baseline-corrected contra/ipsi ERFs from posterior ROIs (auto-selected or
   explicitly overridden by user-provided magnetometer channel names),
2) ipsi-minus-contra waveform,
3) time-resolved decoding accuracy curves for unbalanced and class-balanced
   (undersampled) decoding.

This script is intended for the DAD MEG sanity-check workflow and writes all
outputs into a dedicated sanity-check derivatives folder.

Usage examples
1) Default run for sub-01, ses-01, always_visible:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/sanity_checks/01_dot_side_decoding.py

2) Run another subject/session:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/sanity_checks/01_dot_side_decoding.py \
     --subject 02 --session 01 --condition always_visible

3) Explicit input/output overrides:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/sanity_checks/01_dot_side_decoding.py \
     --epochs-fif /abs/path/to/sub-01_ses-01_task-Attention_run-all_epo.fif \
     --model-mat /abs/path/to/MovDot_Sub01_V28_rescueTraject.mat \
     --out-dir /abs/path/to/output/folder

4) Force specific ERF ROIs:
   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/MEG/code/sanity_checks/01_dot_side_decoding.py \
     --left-roi-mag-channels MEG1731 \
     --right-roi-mag-channels MEG2241

Inputs
- Epoched FIF containing the requested condition and metadata column
  `trial_sequence_identity`.
- Model MAT with `xySeqs` entries containing `condition_label`, `sequence`,
  and `xy`, plus `Cfg.rectSize` for center-threshold definition.

Outputs
- Figure PNG: three-panel sanity-check plot.
- Scores CSV: per-timepoint unbalanced and balanced decoding accuracies.
- Summary JSON: settings, channel sets, balancing strategy, and metrics.
- MNE HTML report: summary text + embedded three-panel figure + whole-sensor
  ERF layout figures for magnetometers and RMS-merged planar gradiometers.

Key assumptions
- `trial_sequence_identity` spans 1..20 and maps to model `sequence`.
- Dot side label is `left` when model x < `Cfg.rectSize(1)/2`, else `right`.
- Balanced decoding uses random undersampling (majority -> minority count)
  inside each CV train/test fold; unbalanced decoding uses original fold data.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
import scipy.io as sio
from mne.io.constants import FIFF
from mne.transforms import apply_trans


# -----------------------------------------------------------------------------
# CLI and path resolution
# -----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    """Create CLI parser for dot-side decoding sanity checks."""
    parser = argparse.ArgumentParser(description="Dot-side decoding sanity check")
    parser.add_argument("--subject", type=str, default="01", help="Subject id without 'sub-' prefix.")
    parser.add_argument("--session", type=str, default="01", help="Session id without 'ses-' prefix.")
    parser.add_argument(
        "--condition",
        type=str,
        default="always_visible",
        help="Condition name in epochs.event_id and model condition_label.",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="MEG project root (default inferred from script location).",
    )
    parser.add_argument(
        "--dad-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="DAD root containing experiment input files.",
    )
    parser.add_argument(
        "--epochs-fif",
        type=Path,
        default=None,
        help="Optional explicit epochs FIF path. If omitted, derived from subject/session.",
    )
    parser.add_argument(
        "--model-mat",
        type=Path,
        default=None,
        help="Optional explicit observed trajectory MAT path. If omitted, derived from subject.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Optional explicit output directory. If omitted, uses derivatives/sanity_checks.",
    )
    parser.add_argument(
        "--decode-sfreq",
        type=float,
        default=120.0,
        help="Sampling frequency used for decoding after resampling.",
    )
    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Stratified CV folds for time-resolved decoding.",
    )
    parser.add_argument(
        "--left-roi-mag-channels",
        type=str,
        default=None,
        help=(
            "Comma-separated left ROI magnetometer channels for ERFs "
            "(example: MEG1731 or MEG1731,MEG1711). "
            "If omitted, channels are auto-selected."
        ),
    )
    parser.add_argument(
        "--right-roi-mag-channels",
        type=str,
        default=None,
        help=(
            "Comma-separated right ROI magnetometer channels for ERFs "
            "(example: MEG2241 or MEG2241,MEG2231). "
            "If omitted, channels are auto-selected."
        ),
    )
    return parser


def resolve_input_output_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """
    Resolve and validate epochs/model inputs and output directory.

    Data flow
    - Build default paths from subject/session identifiers.
    - Respect explicit CLI path overrides when provided.
    - Validate required input files exist before any analysis runs.
    """
    subject = str(args.subject)
    session = str(args.session)
    project_root = Path(args.project_root).expanduser().resolve()
    dad_root = Path(args.dad_root).expanduser().resolve()

    default_epochs = (
        project_root
        / "derivatives"
        / "analysis"
        / f"sub-{subject}"
        / f"ses-{session}"
        / "meg"
        / f"sub-{subject}_ses-{session}_task-Attention_run-all_epo.fif"
    )
    default_model = (
        dad_root
        / "experiment"
        / "input_files"
        / f"sub{subject}"
        / f"MovDot_Sub{subject}_V28_rescueTraject.mat"
    )
    default_out = (
        project_root
        / "derivatives"
        / "sanity_checks"
        / f"sub-{subject}"
        / f"ses-{session}"
        / "dot_side_decoding"
        / str(args.condition)
    )

    epochs_fif = Path(args.epochs_fif).expanduser().resolve() if args.epochs_fif else default_epochs
    model_mat = Path(args.model_mat).expanduser().resolve() if args.model_mat else default_model
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else default_out

    if not epochs_fif.exists():
        raise FileNotFoundError(f"Epochs FIF not found: {epochs_fif}")
    if not model_mat.exists():
        raise FileNotFoundError(f"Model MAT not found: {model_mat}")

    out_dir.mkdir(parents=True, exist_ok=True)
    return epochs_fif, model_mat, out_dir


# -----------------------------------------------------------------------------
# Geometry and channel selection
# -----------------------------------------------------------------------------

def _meg_head_coords(info: mne.Info, *, meg_kind: str) -> list[dict[str, float | int | str]]:
    """
    Return per-channel head-coordinate metadata for MEG channels.

    Data flow
    - Read MEG sensor positions from channel `loc` vectors.
    - Convert device coordinates to head coordinates with `dev_head_t`.
    - Keep only selected MEG kind (`mag` or `grad`).
    """
    if meg_kind not in {"mag", "grad"}:
        raise ValueError(f"Unsupported meg_kind={meg_kind!r}; expected 'mag' or 'grad'.")

    picks = mne.pick_types(info, meg=meg_kind, eeg=False, eog=False, stim=False, exclude=[])
    pick_set = set(int(i) for i in picks)
    dev_head_t = info.get("dev_head_t")

    rows: list[dict[str, float | int | str]] = []
    for idx, ch in enumerate(info["chs"]):
        if idx not in pick_set:
            continue
        xyz = np.array(ch["loc"][:3], dtype=float)
        if ch["coord_frame"] == FIFF.FIFFV_COORD_DEVICE and dev_head_t is not None:
            xyz = apply_trans(dev_head_t, xyz)
        rows.append(
            {
                "index": idx,
                "name": ch["ch_name"],
                "x": float(xyz[0]),
                "y": float(xyz[1]),
                "z": float(xyz[2]),
            }
        )
    return rows


def select_posterior_decoder_channels(info: mne.Info) -> list[str]:
    """
    Select posterior gradiometers used as decoder features.

    Rule
    - Keep grad channels with head-y below the 45th percentile to focus on
      occipito-parietal sensors while keeping enough channels for stability.
    """
    grad_rows = _meg_head_coords(info, meg_kind="grad")
    y_vals = np.array([row["y"] for row in grad_rows], dtype=float)
    y_threshold = float(np.quantile(y_vals, 0.45))
    selected = [row["name"] for row in grad_rows if row["y"] <= y_threshold]
    if len(selected) < 10:
        raise RuntimeError("Too few decoder channels selected; inspect sensor geometry.")
    return selected


def select_left_right_mag_rois(info: mne.Info) -> tuple[list[str], list[str]]:
    """
    Select left/right posterior 3-channel magnetometer ROIs for ERFs.

    Rule
    - Keep only posterior magnetometers (lowest 20% head-y).
    - Split by hemisphere using sign of head-x.
    - Pick the 3 most posterior channels per side.
    """
    mag_rows = _meg_head_coords(info, meg_kind="mag")
    y_vals = np.array([row["y"] for row in mag_rows], dtype=float)
    y_threshold = float(np.quantile(y_vals, 0.20))

    left_candidates = [row for row in mag_rows if row["y"] <= y_threshold and row["x"] < 0]
    right_candidates = [row for row in mag_rows if row["y"] <= y_threshold and row["x"] > 0]
    left_sorted = sorted(left_candidates, key=lambda row: row["y"])
    right_sorted = sorted(right_candidates, key=lambda row: row["y"])

    left_roi = [row["name"] for row in left_sorted[:3]]
    right_roi = [row["name"] for row in right_sorted[:3]]

    if len(left_roi) < 3 or len(right_roi) < 3:
        raise RuntimeError("Could not identify 3-channel left/right posterior mag ROIs.")
    return left_roi, right_roi


def parse_mag_roi_override(raw_value: str | None) -> list[str] | None:
    """
    Parse comma-separated ROI override text into channel names.

    Data flow
    - Split on commas.
    - Strip surrounding whitespace from each token.
    - Drop empty entries and preserve user ordering.
    """
    if raw_value is None:
        return None
    channels = [token.strip() for token in str(raw_value).split(",")]
    channels = [ch for ch in channels if ch]
    return channels or None


def validate_mag_roi_override(info: mne.Info, channels: list[str], side_label: str) -> list[str]:
    """
    Validate a user-supplied ERF ROI channel list.

    Inputs
    - info: epochs info structure that defines available channel names/types.
    - channels: requested ROI channels in desired averaging order.
    - side_label: text label used in error messages ("left" or "right").

    Output
    - The same channel list if all names exist and are magnetometers.
    """
    if not channels:
        raise ValueError(f"{side_label} ROI override is empty after parsing.")

    existing_names = set(info["ch_names"])
    missing = [ch for ch in channels if ch not in existing_names]
    if missing:
        raise ValueError(
            f"{side_label} ROI contains missing channels: {missing}. "
            "Use channel names present in the loaded epochs file."
        )

    mag_picks = mne.pick_types(info, meg="mag", eeg=False, eog=False, stim=False, exclude=[])
    mag_names = {info["ch_names"][idx] for idx in mag_picks}
    non_mag = [ch for ch in channels if ch not in mag_names]
    if non_mag:
        raise ValueError(
            f"{side_label} ROI contains non-magnetometer channels: {non_mag}. "
            "ERF ROI overrides must use MEG magnetometer channels."
        )

    return channels


def resolve_erf_mag_rois(
    info: mne.Info,
    *,
    left_override: list[str] | None,
    right_override: list[str] | None,
) -> tuple[list[str], list[str]]:
    """
    Resolve ERF ROI channels from either auto-selection or explicit overrides.

    Rule
    - If no overrides are provided, use posterior 3-channel auto-selection per
      hemisphere.
    - If one or both overrides are provided, validate each supplied list and
      use it as-is for that hemisphere.
    """
    left_auto, right_auto = select_left_right_mag_rois(info)
    left_roi = left_auto if left_override is None else validate_mag_roi_override(info, left_override, "left")
    right_roi = right_auto if right_override is None else validate_mag_roi_override(info, right_override, "right")
    return left_roi, right_roi


# -----------------------------------------------------------------------------
# Label construction from model trajectories
# -----------------------------------------------------------------------------

def build_dot_side_labels(
    *,
    model_mat: Path,
    condition: str,
    trial_identities: np.ndarray,
    n_times_target: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Build per-trial per-time binary labels from model x trajectories.

    Outputs
    - labels: shape (n_trials, n_times), values 0=right, 1=left.
    - x_by_identity: shape (20, n_times), model x for identities 1..20.
    - center_x: scalar x threshold defined as `Cfg.rectSize(1)/2`.
    """
    model = sio.loadmat(model_mat, simplify_cells=True)
    xyseqs = model["xySeqs"]
    cfg = model["Cfg"]

    trials = [row for row in xyseqs if str(row.get("condition_label", "")).lower() == condition.lower()]
    if len(trials) != 20:
        raise RuntimeError(f"Expected 20 model trials for {condition!r}, found {len(trials)}.")

    by_identity: dict[int, np.ndarray] = {}
    for row in trials:
        by_identity[int(row["sequence"])] = np.asarray(row["xy"], dtype=float)

    expected = set(range(1, 21))
    found = set(by_identity.keys())
    if found != expected:
        raise RuntimeError(f"Model identity set mismatch: found={sorted(found)}, expected={sorted(expected)}")

    x_by_identity_full = np.stack([by_identity[idx][:, 0] for idx in range(1, 21)], axis=0)
    n_times = min(int(n_times_target), int(x_by_identity_full.shape[1]))
    x_by_identity = x_by_identity_full[:, :n_times]

    rect_size = np.asarray(cfg["rectSize"], dtype=float).reshape(-1)
    if rect_size.size == 0:
        raise RuntimeError("Cfg.rectSize missing or empty in model MAT.")
    center_x = float(rect_size[0] / 2.0)

    labels_by_identity = (x_by_identity < center_x).astype(int)
    if np.min(trial_identities) < 1 or np.max(trial_identities) > 20:
        raise RuntimeError(
            "trial_sequence_identity outside expected range 1..20: "
            f"min={int(np.min(trial_identities))}, max={int(np.max(trial_identities))}"
        )

    labels = labels_by_identity[trial_identities - 1, :]
    return labels, x_by_identity, center_x


# -----------------------------------------------------------------------------
# Decoding helpers
# -----------------------------------------------------------------------------

def decode_per_timepoint(
    data: np.ndarray,
    labels: np.ndarray,
    *,
    n_splits: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Decode labels independently at each timepoint with two strategies.

    Outputs
    - acc_unbalanced: 5-fold CV mean accuracy using original fold data.
    - acc_balanced: 5-fold CV mean accuracy after random undersampling of
      majority class to minority count in both train and test fold subsets.
    """
    n_trials, _n_channels, n_times = data.shape
    if labels.shape != (n_trials, n_times):
        raise ValueError(f"labels shape mismatch: got {labels.shape}, expected {(n_trials, n_times)}")

    acc_unbalanced = np.full(n_times, np.nan, dtype=float)
    acc_balanced = np.full(n_times, np.nan, dtype=float)

    # Decode each sample independently because class labels evolve over time.
    for t_idx in range(n_times):
        y_t = labels[:, t_idx]
        if np.unique(y_t).size < 2:
            continue

        x_t = data[:, :, t_idx]
        split_indices = stratified_kfold_indices(y_t, n_splits=n_splits, random_state=42)
        fold_acc_unbalanced: list[float] = []
        fold_acc_balanced: list[float] = []

        for fold_idx, (train_idx, test_idx) in enumerate(split_indices):
            x_train = x_t[train_idx]
            y_train = y_t[train_idx]
            x_test = x_t[test_idx]
            y_test = y_t[test_idx]

            # Strategy 1: original class proportions.
            y_pred_unbalanced = lda_predict_regularized(x_train, y_train, x_test)
            fold_acc_unbalanced.append(float(accuracy(y_test, y_pred_unbalanced)))

            # Strategy 2: random undersampling to parity in train and test.
            rng = np.random.default_rng(100000 + (t_idx * 997) + fold_idx)
            train_bal_idx = undersample_majority_to_minority_indices(y_train, rng=rng)
            test_bal_idx = undersample_majority_to_minority_indices(y_test, rng=rng)
            if train_bal_idx.size == 0 or test_bal_idx.size == 0:
                continue
            y_pred_balanced = lda_predict_regularized(
                x_train[train_bal_idx],
                y_train[train_bal_idx],
                x_test[test_bal_idx],
            )
            fold_acc_balanced.append(float(accuracy(y_test[test_bal_idx], y_pred_balanced)))

        if fold_acc_unbalanced:
            acc_unbalanced[t_idx] = float(np.mean(fold_acc_unbalanced))
        if fold_acc_balanced:
            acc_balanced[t_idx] = float(np.mean(fold_acc_balanced))

    return acc_unbalanced, acc_balanced


def stratified_kfold_indices(
    y: np.ndarray,
    *,
    n_splits: int,
    random_state: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    Build stratified K-fold train/test indices with NumPy only.

    Data flow
    - Shuffle each class index list.
    - Split class indices into `n_splits` chunks.
    - Combine class chunks per fold to preserve class ratios in test sets.
    """
    y = np.asarray(y, dtype=int)
    unique = np.unique(y)
    if unique.size != 2:
        raise ValueError("Stratified split expects binary labels.")

    rng = np.random.default_rng(random_state)
    per_class_folds: dict[int, list[np.ndarray]] = {}
    for cls in unique:
        cls_idx = np.flatnonzero(y == cls)
        rng.shuffle(cls_idx)
        per_class_folds[int(cls)] = list(np.array_split(cls_idx, n_splits))

    all_idx = np.arange(y.size, dtype=int)
    splits: list[tuple[np.ndarray, np.ndarray]] = []
    for fold in range(n_splits):
        test_idx = np.concatenate([per_class_folds[int(cls)][fold] for cls in unique])
        test_idx = np.sort(test_idx)
        train_mask = np.ones(y.size, dtype=bool)
        train_mask[test_idx] = False
        train_idx = all_idx[train_mask]
        splits.append((train_idx, test_idx))
    return splits


def accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute standard classification accuracy."""
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    if y_true.size == 0:
        return np.nan
    return float(np.mean(y_true == y_pred))


def undersample_majority_to_minority_indices(
    y: np.ndarray,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Balance binary labels by random undersampling.

    Strategy details
    - Find minority class count.
    - Randomly sample that many rows from both classes without replacement.
    - Shuffle selected indices for downstream model fitting.
    """
    y = np.asarray(y, dtype=int)
    classes = np.unique(y)
    if classes.size != 2:
        return np.array([], dtype=int)
    idx_a = np.flatnonzero(y == classes[0])
    idx_b = np.flatnonzero(y == classes[1])
    if idx_a.size == 0 or idx_b.size == 0:
        return np.array([], dtype=int)
    n_keep = min(idx_a.size, idx_b.size)
    sel_a = rng.choice(idx_a, size=n_keep, replace=False)
    sel_b = rng.choice(idx_b, size=n_keep, replace=False)
    balanced_idx = np.concatenate([sel_a, sel_b])
    rng.shuffle(balanced_idx)
    return balanced_idx


def lda_predict_regularized(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    *,
    reg: float = 1e-2,
) -> np.ndarray:
    """
    Predict binary labels with regularized LDA.

    Non-obvious implementation detail
    - Features are z-scored with train-only statistics.
    - Pooled covariance gets ridge regularization before inversion.
    """
    y_train = np.asarray(y_train, dtype=int)
    x_train = np.asarray(x_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)

    mu = x_train.mean(axis=0, keepdims=True)
    sigma = x_train.std(axis=0, keepdims=True)
    sigma[sigma == 0] = 1.0
    x_tr = (x_train - mu) / sigma
    x_te = (x_test - mu) / sigma

    x0 = x_tr[y_train == 0]
    x1 = x_tr[y_train == 1]
    m0 = x0.mean(axis=0)
    m1 = x1.mean(axis=0)

    c0 = np.cov(x0, rowvar=False, bias=False)
    c1 = np.cov(x1, rowvar=False, bias=False)
    if c0.ndim == 0:
        c0 = np.atleast_2d(c0)
    if c1.ndim == 0:
        c1 = np.atleast_2d(c1)
    pooled = ((x0.shape[0] - 1) * c0 + (x1.shape[0] - 1) * c1) / max(x0.shape[0] + x1.shape[0] - 2, 1)
    pooled = pooled + reg * np.eye(pooled.shape[0], dtype=float)

    w = np.linalg.solve(pooled, (m1 - m0))
    b = -0.5 * float((m1 + m0) @ w)
    score = x_te @ w + b
    return (score > 0.0).astype(int)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def build_three_panel_plot(
    *,
    times: np.ndarray,
    contra_erf_ft: np.ndarray,
    ipsi_erf_ft: np.ndarray,
    ipsi_minus_contra_ft: np.ndarray,
    accuracy_unbalanced: np.ndarray,
    accuracy_balanced: np.ndarray,
    left_roi: list[str],
    right_roi: list[str],
    n_decoder_channels: int,
    erf_baseline_window: tuple[float, float],
) -> plt.Figure:
    """Create 3-panel figure: ERFs, ipsi-contralateral, and decoding accuracy."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
    ax_top, ax_mid, ax_bottom = axes

    # Top panel: baseline-corrected contra/ipsi ERFs.
    roi_count_label = f"L={len(left_roi)} mag, R={len(right_roi)} mag"
    ax_top.plot(times, contra_erf_ft, color="#1f77b4", lw=2, label=f"Contra posterior ROI ({roi_count_label})")
    ax_top.plot(times, ipsi_erf_ft, color="#d62728", lw=2, label=f"Ipsi posterior ROI ({roi_count_label})")
    ax_top.axvline(0.0, color="k", lw=1, ls="--", alpha=0.7)
    ax_top.set_ylabel("Amplitude (fT)")
    ax_top.set_title(
        "Always_visible ERFs (contra vs ipsi, posterior ROI means, "
        f"baseline {erf_baseline_window[0]:.1f}..{erf_baseline_window[1]:.1f} s)"
    )
    ax_top.legend(loc="upper right", frameon=True)

    # Middle panel: signed lateralization summary.
    ax_mid.plot(times, ipsi_minus_contra_ft, color="#9467bd", lw=2, label="Ipsi - Contra")
    ax_mid.axhline(0.0, color="k", lw=1, ls=":", alpha=0.8)
    ax_mid.axvline(0.0, color="k", lw=1, ls="--", alpha=0.7)
    ax_mid.set_ylabel("Amplitude (fT)")
    ax_mid.set_title("Ipsi - Contra")
    ax_mid.legend(loc="upper right", frameon=True)

    # Bottom panel: decoding accuracy comparison.
    ax_bottom.plot(
        times,
        accuracy_unbalanced,
        color="#2ca02c",
        lw=2,
        label="Accuracy (unbalanced train, 5-fold CV)",
    )
    ax_bottom.plot(
        times,
        accuracy_balanced,
        color="#ff7f0e",
        lw=2,
        label="Accuracy (balanced via random undersampling, 5-fold CV)",
    )
    ax_bottom.axhline(0.5, color="k", lw=1, ls=":", alpha=0.8, label="Chance = 0.5")
    ax_bottom.axvline(0.0, color="k", lw=1, ls="--", alpha=0.7)
    ax_bottom.set_xlabel("Time from condition onset (s)")
    ax_bottom.set_ylabel("Accuracy")
    ax_bottom.set_ylim(0.35, 1.0)
    ax_bottom.set_title(f"Dot side decoding accuracy from posterior gradiometers (n={n_decoder_channels})")
    ax_bottom.legend(loc="upper right", frameon=True)

    info_text = f"Left ROI: {', '.join(left_roi)}\nRight ROI: {', '.join(right_roi)}"
    ax_top.text(
        0.01,
        0.01,
        info_text,
        transform=ax_top.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"},
    )

    return fig


def build_sensor_layout_topo_figures(
    *,
    epochs_erf: mne.Epochs,
    labels: np.ndarray,
    condition_label: str,
) -> tuple[plt.Figure, plt.Figure]:
    """
    Build whole-sensor ipsi-minus-contra ERF layout plots.

    Data flow
    - Compute per-sensor ipsi-minus-contra time courses from baseline-corrected
      epochs using dynamic dot-side labels at each sample.
    - Build one evoked object for magnetometers and one for planar
      gradiometers.
    - Plot magnetometers with `plot_topo`.
    - Plot planar gradiometers with `merge_grads=True`, which combines each
      orthogonal grad pair using RMS magnitude (as in FLUX-style usage).

    Outputs
    - fig_mag_topo: ipsi-minus-contra ERFs over the full magnetometer layout.
    - fig_grad_topo_rms: ipsi-minus-contra ERFs over planar grad layout with
      RMS pair merge.
    """
    def _ipsi_minus_contra_evoked_for_kind(meg_kind: str) -> mne.EvokedArray:
        """
        Compute ipsi-minus-contra evoked traces for one MEG sensor kind.

        Non-obvious logic
        - "Ipsi" and "contra" are defined per channel hemisphere and per sample:
          left-hemisphere channels are ipsi when the dot is on the left, while
          right-hemisphere channels are ipsi when the dot is on the right.
        """
        rows = _meg_head_coords(epochs_erf.info, meg_kind=meg_kind)
        if not rows:
            raise RuntimeError(f"No {meg_kind} channels found for ipsi-contra topography.")

        channel_names = [str(row["name"]) for row in rows]
        x_coords = np.array([float(row["x"]) for row in rows], dtype=float)

        # Drop potential midline channels to avoid ambiguous ipsi/contra mapping.
        keep_mask = ~np.isclose(x_coords, 0.0, atol=1e-9)
        if not np.any(keep_mask):
            raise RuntimeError(f"No lateralized {meg_kind} channels remain after excluding x==0.")
        channel_names = [name for name, keep in zip(channel_names, keep_mask, strict=True) if keep]
        x_coords = x_coords[keep_mask]

        data = epochs_erf.get_data(picks=channel_names)  # (n_trials, n_channels, n_times)
        if labels.shape != (data.shape[0], data.shape[2]):
            raise RuntimeError(
                "Label matrix shape does not match ERF data shape for ipsi-contra topography."
            )
        dot_left_mask = labels.astype(bool)[:, np.newaxis, :]  # (n_trials, 1, n_times)
        channel_is_left = (x_coords < 0.0)[np.newaxis, :, np.newaxis]  # (1, n_channels, 1)

        # Per channel/time sample, choose ipsi and contra trial subsets.
        ipsi_mask = np.where(channel_is_left, dot_left_mask, ~dot_left_mask)
        contra_mask = ~ipsi_mask

        ipsi_count = ipsi_mask.sum(axis=0, dtype=float)
        contra_count = contra_mask.sum(axis=0, dtype=float)
        if np.any(ipsi_count == 0.0) or np.any(contra_count == 0.0):
            raise RuntimeError("Encountered empty ipsi/contra bin while building topography evokeds.")

        ipsi_sum = np.where(ipsi_mask, data, 0.0).sum(axis=0)
        contra_sum = np.where(contra_mask, data, 0.0).sum(axis=0)
        ipsi_mean = ipsi_sum / ipsi_count
        contra_mean = contra_sum / contra_count
        ipsi_minus_contra = ipsi_mean - contra_mean

        picks = mne.pick_channels(epochs_erf.ch_names, include=channel_names, ordered=True)
        picked_info = mne.pick_info(epochs_erf.info, picks)
        return mne.EvokedArray(
            ipsi_minus_contra,
            picked_info,
            tmin=float(epochs_erf.times[0]),
            nave=int(len(epochs_erf)),
            comment=f"ipsi_minus_contra_{meg_kind}",
        )

    evoked_mag_diff = _ipsi_minus_contra_evoked_for_kind("mag")
    evoked_grad_diff = _ipsi_minus_contra_evoked_for_kind("grad")

    # Magnetometer layout ERFs across all available mag channels.
    fig_mag_topo = evoked_mag_diff.plot_topo(
        title=f"{condition_label} | Magnetometers (ipsi - contra)",
        show=False,
    )

    # RMS-combined planar gradiometer layout (orthogonal pair merge).
    fig_grad_topo_rms = evoked_grad_diff.plot_topo(
        title=f"{condition_label} | Planar gradiometers (ipsi - contra, RMS pair merge)",
        merge_grads=True,
        show=False,
    )

    return fig_mag_topo, fig_grad_topo_rms


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

def main() -> None:
    """Run the end-to-end sanity-check decoding workflow."""
    args = build_arg_parser().parse_args()
    epochs_fif, model_mat, out_dir = resolve_input_output_paths(args)
    left_roi_override = parse_mag_roi_override(args.left_roi_mag_channels)
    right_roi_override = parse_mag_roi_override(args.right_roi_mag_channels)

    # -------------------------------------------------------------------------
    # Load epochs and prepare separate decode/ERF branches.
    # -------------------------------------------------------------------------
    epochs_all = mne.read_epochs(epochs_fif, preload=True, verbose=False)
    if args.condition not in epochs_all.event_id:
        raise RuntimeError(f"Condition {args.condition!r} not in epochs.event_id={epochs_all.event_id}.")

    epochs_condition = epochs_all[args.condition].copy()
    if epochs_condition.metadata is None or "trial_sequence_identity" not in epochs_condition.metadata:
        raise RuntimeError("Epoch metadata missing required column 'trial_sequence_identity'.")

    decode_window_tmax = 2.667
    epochs_decode = epochs_condition.copy().crop(tmin=0.0, tmax=decode_window_tmax).resample(args.decode_sfreq)
    erf_baseline_window = (-0.4, -0.2)
    epochs_erf = epochs_condition.copy().resample(args.decode_sfreq).apply_baseline(erf_baseline_window)
    epochs_erf = epochs_erf.crop(tmin=0.0, tmax=min(decode_window_tmax, float(epochs_erf.tmax)))

    trial_identities = epochs_decode.metadata["trial_sequence_identity"].to_numpy(dtype=int)

    # -------------------------------------------------------------------------
    # Build model-derived dot-side labels and align time axes.
    # -------------------------------------------------------------------------
    labels, _x_by_identity, center_x = build_dot_side_labels(
        model_mat=model_mat,
        condition=args.condition,
        trial_identities=trial_identities,
        n_times_target=len(epochs_decode.times),
    )
    n_times_aligned = labels.shape[1]
    if n_times_aligned < len(epochs_decode.times):
        target_tmax = float(epochs_decode.times[n_times_aligned - 1])
        epochs_decode = epochs_decode.copy().crop(tmin=float(epochs_decode.times[0]), tmax=target_tmax)
        epochs_erf = epochs_erf.copy().crop(tmin=float(epochs_erf.times[0]), tmax=target_tmax)
    times = epochs_decode.times.copy()

    # -------------------------------------------------------------------------
    # Select channels and build contra/ipsi waveforms.
    # -------------------------------------------------------------------------
    decoder_channels = select_posterior_decoder_channels(epochs_decode.info)
    left_roi, right_roi = resolve_erf_mag_rois(
        epochs_erf.info,
        left_override=left_roi_override,
        right_override=right_roi_override,
    )

    decode_data = epochs_decode.get_data(picks=decoder_channels)
    left_data = epochs_erf.get_data(picks=left_roi)
    right_data = epochs_erf.get_data(picks=right_roi)

    # Dynamic contra/ipsi definition by dot side at each sample:
    # label=1 (dot left): contra=right ROI, ipsi=left ROI
    # label=0 (dot right): contra=left ROI, ipsi=right ROI
    left_roi_mean = left_data.mean(axis=1)
    right_roi_mean = right_data.mean(axis=1)
    dot_left_mask = labels.astype(bool)
    contra_samples = np.where(dot_left_mask, right_roi_mean, left_roi_mean)
    ipsi_samples = np.where(dot_left_mask, left_roi_mean, right_roi_mean)
    contra_erf_ft = contra_samples.mean(axis=0) * 1e15
    ipsi_erf_ft = ipsi_samples.mean(axis=0) * 1e15
    ipsi_minus_contra_ft = ipsi_erf_ft - contra_erf_ft

    # -------------------------------------------------------------------------
    # Decode timepoint-by-timepoint with unbalanced and balanced strategies.
    # -------------------------------------------------------------------------
    accuracy_unbalanced, accuracy_balanced = decode_per_timepoint(
        decode_data,
        labels,
        n_splits=args.n_splits,
    )

    valid_unbalanced = np.isfinite(accuracy_unbalanced)
    valid_balanced = np.isfinite(accuracy_balanced)
    if not np.any(valid_unbalanced) or not np.any(valid_balanced):
        raise RuntimeError("No valid decoding accuracies were produced.")

    peak_idx_unbalanced = int(np.nanargmax(accuracy_unbalanced))
    peak_idx_balanced = int(np.nanargmax(accuracy_balanced))
    peak_time_unbalanced = float(times[peak_idx_unbalanced])
    peak_time_balanced = float(times[peak_idx_balanced])
    peak_accuracy_unbalanced = float(accuracy_unbalanced[peak_idx_unbalanced])
    peak_accuracy_balanced = float(accuracy_balanced[peak_idx_balanced])
    mean_accuracy_unbalanced = float(np.nanmean(accuracy_unbalanced))
    mean_accuracy_balanced = float(np.nanmean(accuracy_balanced))

    # -------------------------------------------------------------------------
    # Export figure and tabular metrics.
    # -------------------------------------------------------------------------
    fig = build_three_panel_plot(
        times=times,
        contra_erf_ft=contra_erf_ft,
        ipsi_erf_ft=ipsi_erf_ft,
        ipsi_minus_contra_ft=ipsi_minus_contra_ft,
        accuracy_unbalanced=accuracy_unbalanced,
        accuracy_balanced=accuracy_balanced,
        left_roi=left_roi,
        right_roi=right_roi,
        n_decoder_channels=len(decoder_channels),
        erf_baseline_window=erf_baseline_window,
    )

    fig_path = out_dir / f"{args.condition}_dot_side_decoding_3panel.png"
    fig.savefig(fig_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    # Build full-layout ERF topographies for report diagnostics.
    fig_mag_topo, fig_grad_topo_rms = build_sensor_layout_topo_figures(
        epochs_erf=epochs_erf,
        labels=labels,
        condition_label=args.condition,
    )

    scores_csv = out_dir / f"{args.condition}_dot_side_decoding_scores.csv"
    pd.DataFrame(
        {
            "time_sec": times,
            "accuracy_unbalanced_train": accuracy_unbalanced,
            "accuracy_balanced_undersampled_sets": accuracy_balanced,
            "left_fraction": labels.mean(axis=0),
        }
    ).to_csv(scores_csv, index=False)

    summary = {
        "condition": args.condition,
        "subject": str(args.subject),
        "session": str(args.session),
        "n_trials": int(len(epochs_decode)),
        "n_times": int(len(times)),
        "sfreq_hz": float(epochs_decode.info["sfreq"]),
        "n_decoder_channels": int(len(decoder_channels)),
        "decoder_channels": decoder_channels,
        "left_roi_mag_channels": left_roi,
        "right_roi_mag_channels": right_roi,
        "erf_baseline_window_sec": [float(erf_baseline_window[0]), float(erf_baseline_window[1])],
        "balancing_strategy": (
            "Random undersampling of majority class to minority count inside each "
            "CV train and test fold. No oversampling is used."
        ),
        "center_x": center_x,
        "peak_accuracy_unbalanced_train": peak_accuracy_unbalanced,
        "peak_time_sec_unbalanced_train": peak_time_unbalanced,
        "mean_accuracy_unbalanced_train": mean_accuracy_unbalanced,
        "peak_accuracy_balanced_train": peak_accuracy_balanced,
        "peak_time_sec_balanced_train": peak_time_balanced,
        "mean_accuracy_balanced_train": mean_accuracy_balanced,
        "figure_path": str(fig_path),
        "scores_csv": str(scores_csv),
        "epochs_fif": str(epochs_fif),
        "model_mat": str(model_mat),
    }

    summary_json = out_dir / f"{args.condition}_dot_side_decoding_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # -------------------------------------------------------------------------
    # Export MNE report with metrics and figure.
    # -------------------------------------------------------------------------
    report = mne.Report(
        title=(
            f"Dot-side decoding sanity check | sub-{args.subject} ses-{args.session} "
            f"| {args.condition}"
        )
    )

    report.add_html(
        "<pre>"
        + html.escape(
            "Dot-side decoding summary\n"
            f"condition={summary['condition']}\n"
            f"subject=sub-{summary['subject']}\n"
            f"session=ses-{summary['session']}\n"
            f"n_trials={summary['n_trials']}\n"
            f"n_times={summary['n_times']} @ {summary['sfreq_hz']} Hz\n"
            f"n_decoder_channels={summary['n_decoder_channels']}\n"
            f"left_roi={summary['left_roi_mag_channels']}\n"
            f"right_roi={summary['right_roi_mag_channels']}\n"
            f"erf_baseline_window_sec={summary['erf_baseline_window_sec']}\n"
            f"balancing_strategy={summary['balancing_strategy']}\n"
            f"peak_accuracy_unbalanced_train={summary['peak_accuracy_unbalanced_train']:.4f} "
            f"at t={summary['peak_time_sec_unbalanced_train']:.3f} s\n"
            f"mean_accuracy_unbalanced_train={summary['mean_accuracy_unbalanced_train']:.4f}\n"
            f"peak_accuracy_balanced_train={summary['peak_accuracy_balanced_train']:.4f} "
            f"at t={summary['peak_time_sec_balanced_train']:.3f} s\n"
            f"mean_accuracy_balanced_train={summary['mean_accuracy_balanced_train']:.4f}\n"
            f"scores_csv={summary['scores_csv']}\n"
            f"figure_path={summary['figure_path']}"
        )
        + "</pre>",
        title="Decoding summary",
        section="sanity_check",
        tags=("summary", "decoding", "sanity_check"),
    )

    fig_for_report = plt.figure(figsize=(12, 8))
    img = plt.imread(fig_path)
    plt.imshow(img)
    plt.axis("off")
    report.add_figure(
        fig_for_report,
        title="ERF + lateralization + decoding accuracy",
        section="sanity_check",
        tags=("figure", "decoding", "sanity_check"),
    )
    plt.close(fig_for_report)
    report.add_figure(
        fig_mag_topo,
        title="ERF Topography | Magnetometers (ipsi - contra)",
        section="sanity_check",
        tags=("figure", "erf", "magnetometers", "sanity_check"),
    )
    plt.close(fig_mag_topo)
    report.add_figure(
        fig_grad_topo_rms,
        title="ERF Topography | Planar Gradiometers (ipsi - contra, RMS pair merge)",
        section="sanity_check",
        tags=("figure", "erf", "gradiometers", "rms", "sanity_check"),
    )
    plt.close(fig_grad_topo_rms)

    report_path = out_dir / f"{args.condition}_dot_side_decoding_report.html"
    report.save(report_path, overwrite=True, open_browser=False)

    print("Sanity-check decoding complete.")
    print(f"Figure: {fig_path}")
    print(f"Scores CSV: {scores_csv}")
    print(f"Summary JSON: {summary_json}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
