#!/usr/bin/env python3
"""
generate_expected_figures.py

Purpose:
    Generate reproducible figure assets for the MoveDot1 MEG meeting
    presentation. The figure set combines:
      1) real Sub80 trajectory-derived visuals, and
      2) deterministic synthetic dRSA expectations that are explicitly plotted
         as neural-time (x) by model-time (y) matrices.

Why this script exists:
    - Keep slide figures tied to repository data and assumptions.
    - Make expected-result visualizations easy to iterate and rebuild.
    - Provide matrix-style dRSA panels that directly match the hypotheses.

Usage examples (from repo root):
    python control-center/presentation/generate_expected_figures.py
    python control-center/presentation/generate_expected_figures.py --subject 80
    python control-center/presentation/generate_expected_figures.py \
        --subject 80 \
        --outdir control-center/presentation

Inputs:
    - experiment/input_files/MovDot_SubXX.mat
    - experiment/input_files/MovDot_SubXX_predicted.mat
    - simulations/input/MovDot_SubXX_nondeviant.mat
    - simulations/input/MovDot_SubXX_deviant.mat
    - simulations/input/MovDot_SubXX_predicted_deviant.mat

Outputs (written under <outdir>/figures):
    - fig01_paths_sub80_realistic.png
    - fig02_deviant_predicted_divergence.png
    - fig03_mock_predrift_lag.png
    - fig04_mock_postdeviant_switch.png
    - fig05_mock_drsa_matrix.png

Key assumptions:
    - Simulation input files already exist for the requested subject.
    - Path arrays follow shape: trials x 2 x time.
    - Synthetic expected-result panels are illustrative (not empirical dRSA).
    - Synthetic dRSA matrix duration is fixed to 2.66 s to match full trial time.
    - Dependencies are limited to numpy, scipy, matplotlib.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


# Deterministic seed so expected-pattern figures are reproducible across runs.
SYNTHETIC_SEED = 20260224
# Fixed synthetic matrix duration (s) aligned with full trial timing.
SYNTHETIC_TRIAL_DURATION_S = 2.66


@dataclass
class FigureInputs:
    """Container for loaded trajectory arrays and timing metadata."""

    subject: int
    fps: float
    trial_duration_s: float
    deviant_onset_s: float
    nondeviant_dot1: np.ndarray
    nondeviant_dot2: np.ndarray
    deviant_dot1: np.ndarray
    deviant_dot2: np.ndarray
    predicted_deviant_dot1: np.ndarray
    predicted_deviant_dot2: np.ndarray


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for subject selection and output directory."""

    parser = argparse.ArgumentParser(
        description=(
            "Generate realistic and synthetic dRSA-style figures for the MEG presentation."
        )
    )
    parser.add_argument(
        "--subject",
        type=int,
        default=80,
        help="Subject identifier used to locate MovDot_SubXX input files (default: 80).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Presentation directory where the figures/ folder will be created.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Data loading helpers
# Data flow: CLI subject -> repo file paths -> validated arrays + timing info.
# ---------------------------------------------------------------------------

def _require_file(path: Path) -> None:
    """Raise an explicit error if an expected file is missing."""

    if not path.is_file():
        raise FileNotFoundError(f"Required input file not found: {path}")


def _loadmat(path: Path) -> dict:
    """Load MATLAB .mat file as a plain dictionary with struct expansion."""

    return sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)


def _validate_paths_array(name: str, arr: np.ndarray) -> np.ndarray:
    """Validate expected path tensor format and return float ndarray."""

    arr = np.asarray(arr, dtype=float)
    if arr.ndim != 3 or arr.shape[1] != 2:
        raise ValueError(f"{name} must have shape (trials, 2, time); got {arr.shape}")
    return arr


def _compute_divergence(observed: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    """Compute frame-wise Euclidean divergence for each trial."""

    obs = np.transpose(observed, (0, 2, 1))
    pred = np.transpose(predicted, (0, 2, 1))
    return np.linalg.norm(obs - pred, axis=2)


def _infer_deviant_onset_s(
    deviant_dot1: np.ndarray,
    deviant_dot2: np.ndarray,
    predicted_dot1: np.ndarray,
    predicted_dot2: np.ndarray,
    fps: float,
) -> float:
    """Infer deviant onset as first frame where observed and predicted diverge."""

    dot1_div = _compute_divergence(deviant_dot1, predicted_dot1)
    dot2_div = _compute_divergence(deviant_dot2, predicted_dot2)
    mean_div = np.concatenate([dot1_div, dot2_div], axis=0).mean(axis=0)

    nonzero = np.where(mean_div > 1e-9)[0]
    onset_frame = int(nonzero[0]) if nonzero.size else int(0.5 * mean_div.size)
    return onset_frame / fps


def load_figure_inputs(subject: int, repo_root: Path) -> FigureInputs:
    """Load required real-data assets used across all figures."""

    experiment_dir = repo_root / "experiment" / "input_files"
    sim_dir = repo_root / "simulations" / "input"

    movdot_path = experiment_dir / f"MovDot_Sub{subject:02d}.mat"
    movdot_predicted_path = experiment_dir / f"MovDot_Sub{subject:02d}_predicted.mat"
    nondeviant_path = sim_dir / f"MovDot_Sub{subject:02d}_nondeviant.mat"
    deviant_path = sim_dir / f"MovDot_Sub{subject:02d}_deviant.mat"
    predicted_deviant_path = sim_dir / f"MovDot_Sub{subject:02d}_predicted_deviant.mat"

    for required_path in (
        movdot_path,
        movdot_predicted_path,
        nondeviant_path,
        deviant_path,
        predicted_deviant_path,
    ):
        _require_file(required_path)

    movdot_data = _loadmat(movdot_path)
    cfg = movdot_data.get("Cfg")
    if cfg is None or not hasattr(cfg, "fps"):
        raise KeyError(f"Cfg.fps not found in {movdot_path}")
    fps = float(cfg.fps)

    nondeviant_data = _loadmat(nondeviant_path)
    deviant_data = _loadmat(deviant_path)
    predicted_deviant_data = _loadmat(predicted_deviant_path)

    nondeviant_dot1 = _validate_paths_array(
        "dot1GreenPathsCenterRelative", nondeviant_data["dot1GreenPathsCenterRelative"]
    )
    nondeviant_dot2 = _validate_paths_array(
        "dot2YellowPathsCenterRelative", nondeviant_data["dot2YellowPathsCenterRelative"]
    )
    deviant_dot1 = _validate_paths_array(
        "dot1GreenPathsCenterRelative", deviant_data["dot1GreenPathsCenterRelative"]
    )
    deviant_dot2 = _validate_paths_array(
        "dot2YellowPathsCenterRelative", deviant_data["dot2YellowPathsCenterRelative"]
    )
    predicted_deviant_dot1 = _validate_paths_array(
        "dot1GreenPathsCenterRelative",
        predicted_deviant_data["dot1GreenPathsCenterRelative"],
    )
    predicted_deviant_dot2 = _validate_paths_array(
        "dot2YellowPathsCenterRelative",
        predicted_deviant_data["dot2YellowPathsCenterRelative"],
    )

    # Keep synthetic matrix axes consistent across scripts (V1/V2/V3), rather
    # than deriving duration from frame count/fps.
    trial_duration_s = SYNTHETIC_TRIAL_DURATION_S
    deviant_onset_s = _infer_deviant_onset_s(
        deviant_dot1,
        deviant_dot2,
        predicted_deviant_dot1,
        predicted_deviant_dot2,
        fps,
    )

    return FigureInputs(
        subject=subject,
        fps=fps,
        trial_duration_s=trial_duration_s,
        deviant_onset_s=deviant_onset_s,
        nondeviant_dot1=nondeviant_dot1,
        nondeviant_dot2=nondeviant_dot2,
        deviant_dot1=deviant_dot1,
        deviant_dot2=deviant_dot2,
        predicted_deviant_dot1=predicted_deviant_dot1,
        predicted_deviant_dot2=predicted_deviant_dot2,
    )


# ---------------------------------------------------------------------------
# Figure 1: realistic path samples
# Data flow: real nondeviant/deviant arrays -> sampled trajectory panel.
# ---------------------------------------------------------------------------

def plot_realistic_paths(data: FigureInputs, out_path: Path) -> None:
    """Plot sampled trajectory bundles for nondeviant vs deviant conditions."""

    trial_count = data.nondeviant_dot1.shape[0]
    sample_count = min(18, trial_count)
    sample_idx = np.linspace(0, trial_count - 1, sample_count, dtype=int)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), sharex=True, sharey=True)
    panel_defs = [
        ("Dot 1 trajectories", data.nondeviant_dot1, data.deviant_dot1),
        ("Dot 2 trajectories", data.nondeviant_dot2, data.deviant_dot2),
    ]

    for axis, (title, nondeviant_paths, deviant_paths) in zip(axes, panel_defs):
        for idx in sample_idx:
            axis.plot(
                nondeviant_paths[idx, 0, :],
                nondeviant_paths[idx, 1, :],
                color="#1f77b4",
                alpha=0.22,
                linewidth=1.0,
            )
            axis.plot(
                deviant_paths[idx, 0, :],
                deviant_paths[idx, 1, :],
                color="#d62728",
                alpha=0.22,
                linewidth=1.0,
            )

        axis.set_title(title)
        axis.set_xlabel("Horizontal position (deg, center-relative)")
        axis.set_ylabel("Vertical position (deg, center-relative)")
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.2)

    fig.suptitle(
        (
            f"Sub{data.subject:02d} realistic path samples: "
            "nondeviant (blue) vs deviant (red)"
        ),
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0.0, 1, 0.95])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: observed-vs-predicted divergence
# Data flow: deviant paths + predicted-deviant baseline -> divergence curve.
# ---------------------------------------------------------------------------

def plot_deviant_predicted_divergence(data: FigureInputs, out_path: Path) -> None:
    """Plot mean divergence over time between deviant and predicted trajectories."""

    dot1_div = _compute_divergence(data.deviant_dot1, data.predicted_deviant_dot1)
    dot2_div = _compute_divergence(data.deviant_dot2, data.predicted_deviant_dot2)
    combined_div = np.concatenate([dot1_div, dot2_div], axis=0)

    mean_div = np.nanmean(combined_div, axis=0)
    sem_div = np.nanstd(combined_div, axis=0, ddof=1) / np.sqrt(combined_div.shape[0])
    time_s = np.arange(mean_div.size) / data.fps

    fig, axis = plt.subplots(figsize=(10.5, 5.5))
    axis.plot(time_s, mean_div, color="#111111", linewidth=2.2, label="Mean divergence")
    axis.fill_between(
        time_s,
        mean_div - 1.96 * sem_div,
        mean_div + 1.96 * sem_div,
        color="#888888",
        alpha=0.25,
        label="95% CI",
    )
    axis.axvline(
        data.deviant_onset_s,
        color="#d62728",
        linestyle="--",
        linewidth=1.8,
        label=f"Empirical deviant onset (~{data.deviant_onset_s:.3f} s)",
    )

    axis.set_title(
        f"Sub{data.subject:02d} deviant vs predicted-deviant path divergence"
    )
    axis.set_xlabel("Time within trial (s)")
    axis.set_ylabel("Euclidean distance (deg)")
    axis.grid(alpha=0.25)
    axis.legend(loc="upper left")

    pre_mask = time_s < data.deviant_onset_s
    post_mask = time_s >= data.deviant_onset_s
    pre_mean = float(np.nanmean(mean_div[pre_mask]))
    post_mean = float(np.nanmean(mean_div[post_mask]))
    axis.text(
        0.995,
        0.02,
        f"pre-onset mean: {pre_mean:.5f} deg | post-onset mean: {post_mean:.5f} deg",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------------------
# dRSA-matrix helpers (synthetic)
# Data flow: model-time dynamics assumptions -> neural x model dRSA matrices.
# ---------------------------------------------------------------------------

def _zscore_matrix(matrix: np.ndarray) -> np.ndarray:
    """Z-score matrix values for consistent visual scaling."""

    std = float(np.nanstd(matrix))
    if std < 1e-12:
        return matrix - np.nanmean(matrix)
    return (matrix - np.nanmean(matrix)) / std


def _build_drsa_matrix(
    time_axis: np.ndarray,
    lag_fn: Callable[[np.ndarray], np.ndarray],
    amp_fn: Callable[[np.ndarray], np.ndarray],
    rng: np.random.Generator,
    ridge_sigma_s: float = 0.055,
    noise_sd: float = 0.04,
) -> np.ndarray:
    """Construct a synthetic dRSA matrix (x=neural time, y=model time)."""

    neural_grid, model_grid = np.meshgrid(time_axis, time_axis, indexing="xy")

    lag_s = lag_fn(model_grid)
    amp = np.clip(amp_fn(model_grid), 0.0, None)

    # Main ridge: where neural time matches model time plus lag profile.
    ridge = np.exp(
        -((neural_grid - (model_grid + lag_s)) ** 2) / (2.0 * ridge_sigma_s**2)
    )

    # Weak offset anti-ridge creates contrast around the main ridge.
    anti = np.exp(
        -((neural_grid - (model_grid + lag_s + 0.18)) ** 2)
        / (2.0 * (ridge_sigma_s * 1.6) ** 2)
    )

    matrix = amp * (ridge - 0.22 * anti)
    matrix += noise_sd * rng.standard_normal(size=matrix.shape)
    return _zscore_matrix(matrix)


def _draw_drsa_panel(
    axis: plt.Axes,
    matrix: np.ndarray,
    time_axis: np.ndarray,
    title: str,
    onset_s: float | None,
) -> plt.AxesImage:
    """Render a single dRSA panel with standard axis semantics."""

    image = axis.imshow(
        matrix,
        origin="lower",
        aspect="equal",
        extent=[time_axis[0], time_axis[-1], time_axis[0], time_axis[-1]],
        cmap="coolwarm",
        vmin=-2.4,
        vmax=2.4,
    )
    axis.set_title(title, fontsize=10)
    axis.set_xlabel("Neural time (s)")
    axis.set_ylabel("Model time (s)")
    if onset_s is not None:
        axis.axvline(onset_s, color="black", linestyle="--", linewidth=1.0)
        axis.axhline(onset_s, color="black", linestyle="--", linewidth=1.0)
    return image


# ---------------------------------------------------------------------------
# Figure 3: pre-deviant hypotheses as dRSA matrices
# Data flow: hypothesis-specific lag/amplitude rules -> 3 matrix panels.
# ---------------------------------------------------------------------------

def plot_mock_predrift_lag(data: FigureInputs, out_path: Path, rng: np.random.Generator) -> None:
    """Plot H1/H1a/H1b as synthetic dRSA matrices before deviant onset."""

    pre_time = np.linspace(0.0, data.deviant_onset_s, 150)
    pre_end = pre_time[-1]

    # H1: predictive drift (lag decreases over model time and may cross zero).
    h1_matrix = _build_drsa_matrix(
        pre_time,
        lag_fn=lambda m: 0.12 - 0.18 * (m / pre_end),
        amp_fn=lambda m: 0.45 + 0.55 * (m / pre_end),
        rng=rng,
        ridge_sigma_s=0.050,
        noise_sd=0.035,
    )

    # H1a: sharpening (lag stays sensory-lagged but representation strengthens).
    h1a_matrix = _build_drsa_matrix(
        pre_time,
        lag_fn=lambda m: np.full_like(m, 0.11),
        amp_fn=lambda m: 0.30 + 0.85 * (m / pre_end),
        rng=rng,
        ridge_sigma_s=0.050,
        noise_sd=0.035,
    )

    # H1b: stable lagged representation (roughly constant lag and amplitude).
    h1b_matrix = _build_drsa_matrix(
        pre_time,
        lag_fn=lambda m: np.full_like(m, 0.10),
        amp_fn=lambda m: np.full_like(m, 0.65),
        rng=rng,
        ridge_sigma_s=0.050,
        noise_sd=0.035,
    )

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True, sharey=True)
    images = []
    images.append(_draw_drsa_panel(axes[0], h1_matrix, pre_time, "H1: predictive drift", None))
    images.append(_draw_drsa_panel(axes[1], h1a_matrix, pre_time, "H1a: sharpening", None))
    images.append(_draw_drsa_panel(axes[2], h1b_matrix, pre_time, "H1b: stable lag", None))

    # Shared colorbar enforces direct visual comparability across hypotheses.
    cbar = fig.colorbar(images[-1], ax=axes.ravel().tolist(), shrink=0.86, pad=0.02)
    cbar.set_label("Synthetic dRSA strength (z)")

    fig.suptitle(
        "Expected result A (synthetic): pre-deviant dRSA matrices (x=neural time, y=model time)",
        fontsize=12,
    )
    fig.text(0.01, 0.01, "Synthetic expected pattern", color="#a00000", weight="bold")
    fig.subplots_adjust(left=0.05, right=0.93, bottom=0.10, top=0.86, wspace=0.28)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: post-deviant model competition as dRSA matrices
# Data flow: model-specific post-deviant dynamics -> predicted/deviant/PE panels.
# ---------------------------------------------------------------------------

def plot_mock_postdeviant_switch(
    data: FigureInputs, out_path: Path, rng: np.random.Generator
) -> None:
    """Plot deviant-trial expected dynamics for predicted/deviant/PE models."""

    rel_time = np.linspace(-0.25, 1.0, 180)

    # Predicted model: strong pre-deviant, then fast collapse after deviance.
    predicted_matrix = _build_drsa_matrix(
        rel_time,
        lag_fn=lambda m: np.full_like(m, 0.08),
        amp_fn=lambda m: 0.08 + 0.90 * (1.0 / (1.0 + np.exp((m - 0.05) / 0.06))),
        rng=rng,
        ridge_sigma_s=0.060,
        noise_sd=0.040,
    )

    # Deviant model: weak pre-deviant, then gradual post-deviant emergence.
    deviant_matrix = _build_drsa_matrix(
        rel_time,
        lag_fn=lambda m: 0.14 - 0.06 * (1.0 / (1.0 + np.exp(-(m - 0.24) / 0.09))),
        amp_fn=lambda m: 0.05 + 0.95 * (1.0 / (1.0 + np.exp(-(m - 0.18) / 0.08))),
        rng=rng,
        ridge_sigma_s=0.065,
        noise_sd=0.040,
    )

    # PE model: transient post-deviant burst with short positive lag.
    pe_matrix = _build_drsa_matrix(
        rel_time,
        lag_fn=lambda m: np.full_like(m, 0.05),
        amp_fn=lambda m: 0.92 * np.exp(-((m - 0.12) ** 2) / (2.0 * 0.09**2)),
        rng=rng,
        ridge_sigma_s=0.060,
        noise_sd=0.040,
    )

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True, sharey=True)
    images = []
    images.append(
        _draw_drsa_panel(
            axes[0],
            predicted_matrix,
            rel_time,
            "Deviant trials vs predicted model",
            onset_s=0.0,
        )
    )
    images.append(
        _draw_drsa_panel(
            axes[1],
            deviant_matrix,
            rel_time,
            "Deviant trials vs deviant model",
            onset_s=0.0,
        )
    )
    images.append(
        _draw_drsa_panel(
            axes[2],
            pe_matrix,
            rel_time,
            "Deviant trials vs PE model",
            onset_s=0.0,
        )
    )

    cbar = fig.colorbar(images[-1], ax=axes.ravel().tolist(), shrink=0.86, pad=0.02)
    cbar.set_label("Synthetic dRSA strength (z)")

    fig.suptitle(
        (
            "Expected result B (synthetic): post-deviant competition "
            "(x=neural time, y=model time, t=0 at deviant onset)"
        ),
        fontsize=12,
    )
    fig.text(0.01, 0.01, "Synthetic expected pattern", color="#a00000", weight="bold")
    fig.subplots_adjust(left=0.05, right=0.93, bottom=0.10, top=0.86, wspace=0.28)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: full-trial model-by-trial-type matrix grid
# Data flow: nondeviant/deviant trial assumptions -> 2x3 dRSA matrix grid.
# ---------------------------------------------------------------------------

def plot_mock_drsa_matrix(data: FigureInputs, out_path: Path, rng: np.random.Generator) -> None:
    """Plot full-trial expected dRSA matrices for trial-type/model combinations."""

    full_time = np.linspace(0.0, data.trial_duration_s, 190)
    onset = data.deviant_onset_s

    def nondev_nondev_amp(m: np.ndarray) -> np.ndarray:
        """Piecewise prior profile for nondeviant trials vs nondeviant model."""

        early = np.clip((m - 0.00) / 0.30, 0.0, 1.0)  # 0-300 ms buildup
        mid = 1.0 - 0.18 * np.clip((m - 1.00) / max(onset - 1.00, 1e-6), 0.0, 1.0)
        late = 0.82 + 0.18 * np.clip((m - onset) / max(data.trial_duration_s - onset, 1e-6), 0.0, 1.0)
        out = 0.20 + 0.75 * early
        out = np.where(m >= 0.30, mid, out)
        out = np.where(m >= onset, late, out)
        return out

    # Row 1: nondeviant trials.
    nd_vs_nd = _build_drsa_matrix(
        full_time,
        lag_fn=lambda m: 0.11 - 0.05 * np.clip((m - onset) / max(data.trial_duration_s - onset, 1e-6), 0.0, 1.0),
        amp_fn=nondev_nondev_amp,
        rng=rng,
        ridge_sigma_s=0.060,
        noise_sd=0.040,
    )
    nd_vs_dev = _build_drsa_matrix(
        full_time,
        lag_fn=lambda m: np.full_like(m, 0.10),
        amp_fn=lambda m: np.full_like(m, 0.08),
        rng=rng,
        ridge_sigma_s=0.065,
        noise_sd=0.040,
    )
    nd_vs_pe = _build_drsa_matrix(
        full_time,
        lag_fn=lambda m: np.full_like(m, 0.06),
        amp_fn=lambda m: 0.08 * np.exp(-((m - onset) ** 2) / (2.0 * 0.10**2)),
        rng=rng,
        ridge_sigma_s=0.060,
        noise_sd=0.040,
    )

    # Row 2: deviant trials.
    d_vs_nd = _build_drsa_matrix(
        full_time,
        lag_fn=lambda m: np.full_like(m, 0.08),
        amp_fn=lambda m: 0.85 * (1.0 / (1.0 + np.exp((m - (onset + 0.05)) / 0.07))) + 0.06,
        rng=rng,
        ridge_sigma_s=0.060,
        noise_sd=0.040,
    )
    d_vs_dev = _build_drsa_matrix(
        full_time,
        lag_fn=lambda m: 0.14 - 0.07 * (1.0 / (1.0 + np.exp(-(m - (onset + 0.22)) / 0.12))),
        amp_fn=lambda m: 0.05 + 0.95 * (1.0 / (1.0 + np.exp(-(m - (onset + 0.18)) / 0.10))),
        rng=rng,
        ridge_sigma_s=0.065,
        noise_sd=0.040,
    )
    d_vs_pe = _build_drsa_matrix(
        full_time,
        lag_fn=lambda m: np.full_like(m, 0.05),
        amp_fn=lambda m: 1.05 * np.exp(-((m - (onset + 0.12)) ** 2) / (2.0 * 0.12**2)),
        rng=rng,
        ridge_sigma_s=0.060,
        noise_sd=0.040,
    )

    grid = [
        [nd_vs_nd, nd_vs_dev, nd_vs_pe],
        [d_vs_nd, d_vs_dev, d_vs_pe],
    ]
    titles = [
        "Neural(nondev) vs model(nondev)",
        "Neural(nondev) vs model(deviant)",
        "Neural(nondev) vs model(PE)",
        "Neural(deviant) vs model(nondev)",
        "Neural(deviant) vs model(deviant)",
        "Neural(deviant) vs model(PE)",
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.4), sharex=True, sharey=True)
    idx = 0
    image = None
    for r in range(2):
        for c in range(3):
            image = _draw_drsa_panel(
                axes[r, c],
                grid[r][c],
                full_time,
                titles[idx],
                onset_s=onset,
            )
            if c > 0:
                axes[r, c].set_ylabel("")
            idx += 1

    # Row labels clarify which neural stream is correlated to each model.
    axes[0, 0].text(
        -0.34,
        0.5,
        "Nondeviant trial neural stream",
        transform=axes[0, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=10,
        weight="bold",
    )
    axes[1, 0].text(
        -0.34,
        0.5,
        "Deviant trial neural stream",
        transform=axes[1, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=10,
        weight="bold",
    )

    cbar = fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.88, pad=0.015)
    cbar.set_label("Synthetic dRSA strength (z)")

    fig.suptitle(
        (
            "Expected result C (synthetic): full-trial dRSA matrix grid "
            "(x=neural time, y=model time)"
        ),
        fontsize=12,
    )
    fig.text(0.01, 0.01, "Synthetic expected pattern", color="#a00000", weight="bold")
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.90, wspace=0.28, hspace=0.28)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main orchestration
# Data flow: parse args -> load data -> generate figures -> print outputs.
# ---------------------------------------------------------------------------

def main() -> None:
    """Generate all required figure assets."""

    args = parse_args()
    outdir = args.outdir.resolve()
    figures_dir = outdir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[2]
    inputs = load_figure_inputs(subject=args.subject, repo_root=repo_root)
    synthetic_rng = np.random.default_rng(SYNTHETIC_SEED)

    outputs = {
        "fig01_paths_sub80_realistic.png": lambda p: plot_realistic_paths(inputs, p),
        "fig02_deviant_predicted_divergence.png": lambda p: plot_deviant_predicted_divergence(
            inputs, p
        ),
        "fig03_mock_predrift_lag.png": lambda p: plot_mock_predrift_lag(
            inputs, p, synthetic_rng
        ),
        "fig04_mock_postdeviant_switch.png": lambda p: plot_mock_postdeviant_switch(
            inputs, p, synthetic_rng
        ),
        "fig05_mock_drsa_matrix.png": lambda p: plot_mock_drsa_matrix(
            inputs, p, synthetic_rng
        ),
    }

    for filename, builder in outputs.items():
        target = figures_dir / filename
        builder(target)

    print("Generated figure assets:")
    for filename in outputs:
        print(f"- {figures_dir / filename}")


if __name__ == "__main__":
    main()
