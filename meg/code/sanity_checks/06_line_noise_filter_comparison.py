#!/usr/bin/env python3
"""
06_line_noise_filter_comparison.py

Purpose
Create a reproducible, presentation-ready comparison of residual 50-Hz line
noise in DAD MEG data after three processing choices:

1. No additional cleaning after Maxwell filtering.
2. MNE's default zero-phase FIR notch at 50 Hz.
3. MNE's multitaper sinusoid fitting (``method="spectrum_fit"``) at 50 Hz.

The script reads one continuous post-Maxwell FIF file, crops a fixed-duration
segment, applies each method to an independent copy, and estimates Welch power
spectral densities for all good MEG sensors. To summarize magnetometers and
gradiometers without mixing their physical units, every sensor spectrum is
expressed in dB relative to that sensor's local 50-Hz sidebands. The plotted
trace is the channel median; the zoom panel additionally shows the interquartile
range across channels.

Usage examples
1) Reproduce the default subject-01/block-01 comparison from the repository
   root:

   python meg/code/sanity_checks/06_line_noise_filter_comparison.py

2) Compare a different file and longer time window:

   python meg/code/sanity_checks/06_line_noise_filter_comparison.py \
     --input meg/derivatives/preprocessing/sub-02/ses-01/sss/SUB02_block01_raw_sss.fif \
     --duration 120 --spectrum-fit-window 10

3) Write the figure and supporting values to a separate directory:

   python meg/code/sanity_checks/06_line_noise_filter_comparison.py \
     --output-dir /tmp/dad_line_noise_comparison --overwrite

Inputs
- A continuous post-Maxwell MNE Raw FIF file. The default is subject 01,
  block 01 under ``meg/derivatives/preprocessing``.

Outputs
- ``line_noise_method_comparison.png``: high-resolution raster figure.
- ``line_noise_method_comparison.pdf``: vector figure for slides/manuscripts.
- ``line_noise_method_comparison_values.csv``: plotted line-excess values.
- ``line_noise_method_comparison_metadata.json``: inputs and method details.

Metric definition and assumptions
- The 50-Hz peak is the maximum Welch PSD within +/-0.25 Hz of 50 Hz.
- The local background is the median PSD from 0.5 to 2.0 Hz away on both
  sides of 50 Hz.
- Line excess is ``10*log10(peak/background)`` per sensor and is summarized by
  the median across good MEG sensors. A negative value means the 50-Hz bin is
  below the neighboring spectral floor; for the FIR method this represents a
  spectral hole rather than retained neural information.
- This is a descriptive comparison of one continuous segment, not a
  participant-level statistical analysis and not a ground-truth recovery test.
- Cleaning is deliberately targeted only at 50 Hz so the figure matches the
  reported three-method comparison. Other harmonics are left unchanged.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import mne  # noqa: E402
import numpy as np  # noqa: E402
from scipy import signal  # noqa: E402


# Resolve repository-relative defaults without depending on the launch folder.
REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "meg"
    / "derivatives"
    / "preprocessing"
    / "sub-01"
    / "ses-01"
    / "sss"
    / "SUB01_block01_raw_sss.fif"
)
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT
    / "meg"
    / "derivatives"
    / "reports"
    / "sub-01"
    / "ses-01"
    / "line_noise_comparison"
)

LINE_FREQUENCY_HZ = 50.0
PEAK_HALF_WIDTH_HZ = 0.25
BACKGROUND_INNER_HZ = 0.5
BACKGROUND_OUTER_HZ = 2.0
N_FFT = 8192
N_OVERLAP = 4096

METHOD_ORDER = ("No additional cleaning", "Default FIR notch", "MNE spectrum_fit")
METHOD_STYLE = {
    "No additional cleaning": {"color": "#4B5563", "linestyle": "-"},
    "Default FIR notch": {"color": "#D97706", "linestyle": "--"},
    "MNE spectrum_fit": {"color": "#2563EB", "linestyle": "-."},
}


def parse_args() -> argparse.Namespace:
    """Parse command-line options for source selection and figure generation."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Post-Maxwell continuous raw FIF file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for PNG, PDF, CSV, and metadata outputs.",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=60.0,
        help="Seconds read from the beginning of the recording (default: 60).",
    )
    parser.add_argument(
        "--spectrum-fit-window",
        type=float,
        default=10.0,
        help="Effective MNE spectrum-fit window in seconds (default: 10).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing comparison outputs.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Reject missing inputs and parameter values that cannot support the PSD."""

    if not args.input.is_file():
        raise FileNotFoundError(f"Input FIF does not exist: {args.input}")
    if args.duration <= 10.0:
        raise ValueError("--duration must exceed 10 seconds for a stable comparison.")
    if args.spectrum_fit_window <= 0.0:
        raise ValueError("--spectrum-fit-window must be positive.")


def load_segment(path: Path, duration: float) -> mne.io.BaseRaw:
    """Load at most ``duration`` seconds while preserving Raw metadata."""

    raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
    available_duration = raw.times[-1]
    if available_duration < 10.0:
        raise ValueError(f"Recording is too short for this comparison: {available_duration:.2f} s")

    # Crop before preloading so only the requested interval enters memory.
    raw.crop(tmin=0.0, tmax=min(duration, available_duration), include_tmax=False)
    raw.load_data(verbose="ERROR")
    return raw


def build_method_copies(
    raw: mne.io.BaseRaw,
    picks: np.ndarray,
    spectrum_fit_window: float,
) -> dict[str, mne.io.BaseRaw]:
    """Return independent Raw objects for the three compared processing paths."""

    methods = {"No additional cleaning": raw.copy()}

    # This reproduces MNE's ordinary zero-phase FIR band-stop implementation.
    fir = raw.copy()
    fir.notch_filter(
        freqs=[LINE_FREQUENCY_HZ],
        picks=picks,
        method="fir",
        verbose="ERROR",
    )
    methods["Default FIR notch"] = fir

    # Spectrum fitting estimates and subtracts a sinusoid instead of removing a band.
    spectrum_fit = raw.copy()
    spectrum_fit.notch_filter(
        freqs=[LINE_FREQUENCY_HZ],
        picks=picks,
        method="spectrum_fit",
        filter_length=f"{spectrum_fit_window:g}s",
        verbose="ERROR",
    )
    methods["MNE spectrum_fit"] = spectrum_fit
    return methods


def welch_psd(raw: mne.io.BaseRaw, picks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Estimate channel-wise Welch PSD with settings shared by all methods."""

    data = raw.get_data(picks=picks)
    n_per_seg = min(N_FFT, data.shape[-1])
    n_overlap = min(N_OVERLAP, n_per_seg // 2)
    frequencies, psd = signal.welch(
        data,
        fs=float(raw.info["sfreq"]),
        window="hamming",
        nperseg=n_per_seg,
        noverlap=n_overlap,
        nfft=N_FFT,
        detrend="constant",
        return_onesided=True,
        scaling="density",
        axis=-1,
        average="mean",
    )
    return frequencies, psd


def summarize_spectrum(
    frequencies: np.ndarray,
    psd: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Normalize channel PSDs to local sidebands and compute line-excess summaries."""

    distance = np.abs(frequencies - LINE_FREQUENCY_HZ)
    peak_mask = distance <= PEAK_HALF_WIDTH_HZ
    background_mask = (distance >= BACKGROUND_INNER_HZ) & (distance <= BACKGROUND_OUTER_HZ)
    if not np.any(peak_mask) or not np.any(background_mask):
        raise RuntimeError("PSD frequency grid does not cover the 50-Hz metric windows.")

    # Work in linear power for peak/background ratios, then convert those ratios to dB.
    background = np.median(psd[:, background_mask], axis=1)
    peak = np.max(psd[:, peak_mask], axis=1)
    line_excess_by_channel = 10.0 * np.log10(peak / background)

    normalized_db = 10.0 * np.log10(psd / background[:, np.newaxis])
    median = np.median(normalized_db, axis=0)
    q25, q75 = np.percentile(normalized_db, [25.0, 75.0], axis=0)
    return (
        median,
        q25,
        q75,
        line_excess_by_channel,
        float(np.median(line_excess_by_channel)),
        float(np.percentile(line_excess_by_channel, 90.0)),
    )


def plot_comparison(
    frequencies: np.ndarray,
    summaries: dict[str, dict[str, Any]],
    input_label: str,
    duration: float,
    output_png: Path,
    output_pdf: Path,
) -> None:
    """Render the shared-scale spectra, 50-Hz zoom, and exact excess estimates."""

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.0,
            "axes.titlesize": 12.0,
            "axes.labelsize": 10.0,
            "legend.fontsize": 9.5,
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "axes.edgecolor": "#374151",
            "axes.linewidth": 0.8,
            "savefig.facecolor": "white",
        }
    )

    fig = plt.figure(figsize=(12.2, 6.6), constrained_layout=False)
    grid = fig.add_gridspec(
        1,
        3,
        width_ratios=(1.30, 1.30, 0.82),
        left=0.075,
        right=0.975,
        bottom=0.19,
        top=0.77,
        wspace=0.32,
    )
    ax_broad = fig.add_subplot(grid[0, 0])
    ax_zoom = fig.add_subplot(grid[0, 1])
    ax_metric = fig.add_subplot(grid[0, 2])

    # Draw one consistent encoding in the broadband and zoomed spectral panels.
    for method in METHOD_ORDER:
        style = METHOD_STYLE[method]
        summary = summaries[method]
        for ax in (ax_broad, ax_zoom):
            ax.plot(
                frequencies,
                summary["median"],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.0,
                label=method,
                zorder=3,
            )
        ax_zoom.fill_between(
            frequencies,
            summary["q25"],
            summary["q75"],
            color=style["color"],
            alpha=0.10,
            linewidth=0,
            zorder=1,
        )

    for ax in (ax_broad, ax_zoom):
        ax.axvline(LINE_FREQUENCY_HZ, color="#9CA3AF", linewidth=1.0, zorder=0)
        ax.axhline(0.0, color="#6B7280", linewidth=0.9, linestyle=":", zorder=0)
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Frequency (Hz)")

    ax_broad.set_xlim(35.0, 110.0)
    ax_broad.set_xticks([40, 50, 60, 80, 100])
    ax_broad.set_ylabel("Power relative to local 50-Hz sidebands (dB)")
    ax_broad.set_title("A  Broadband context", loc="left", fontweight="bold")

    ax_zoom.set_xlim(46.0, 54.0)
    ax_zoom.set_xticks([46, 48, 50, 52, 54])
    ax_zoom.set_title("B  Detail around line frequency", loc="left", fontweight="bold")
    ax_zoom.text(
        50.0,
        0.98,
        "50 Hz",
        transform=ax_zoom.get_xaxis_transform(),
        ha="center",
        va="top",
        color="#6B7280",
        fontsize=9.0,
    )

    # Use aligned dots and exact labels so the table values remain easy to quote.
    y_positions = np.arange(len(METHOD_ORDER))[::-1]
    values = [summaries[method]["median_excess_db"] for method in METHOD_ORDER]
    for y_pos, method, value in zip(y_positions, METHOD_ORDER, values):
        style = METHOD_STYLE[method]
        ax_metric.plot(
            [0.0, value],
            [y_pos, y_pos],
            color="#D1D5DB",
            linewidth=2.0,
            zorder=1,
        )
        ax_metric.scatter(
            value,
            y_pos,
            s=90,
            color=style["color"],
            edgecolor="white",
            linewidth=1.0,
            zorder=3,
        )
        ax_metric.text(
            value + 1.2,
            y_pos,
            f"{value:+.1f} dB",
            ha="left",
            va="center",
            color="#111827",
            fontweight="bold",
        )

    ax_metric.axvline(0.0, color="#374151", linewidth=1.0)
    ax_metric.set_xlim(-23.0, 16.0)
    ax_metric.set_ylim(-0.6, len(METHOD_ORDER) - 0.4)
    ax_metric.set_yticks(y_positions, ["No cleaning", "FIR notch", "spectrum_fit"])
    ax_metric.set_xlabel("Median 50-Hz excess (dB)")
    ax_metric.set_title("C  Residual peak", loc="left", fontweight="bold")
    ax_metric.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax_metric.spines[["top", "right", "left"]].set_visible(False)
    ax_metric.tick_params(axis="y", length=0)

    fig.suptitle(
        "50-Hz line-noise cleaning in DAD MEG",
        x=0.075,
        y=0.965,
        ha="left",
        fontsize=18.0,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.075,
        0.905,
        f"{input_label} | first {duration:g} s after Maxwell filtering | good MEG sensors",
        ha="left",
        va="top",
        fontsize=10.5,
        color="#4B5563",
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=METHOD_STYLE[method]["color"],
            linestyle=METHOD_STYLE[method]["linestyle"],
            linewidth=2.2,
            label=method,
        )
        for method in METHOD_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.075, 0.865),
        frameon=False,
        ncol=3,
        borderaxespad=0.0,
        handlelength=2.7,
        columnspacing=1.8,
    )

    fig.text(
        0.075,
        0.075,
        "Welch PSD; channel median (lines), interquartile range (shading in B). "
        "Each sensor is normalized to its own 48-49.5 and 50.5-52 Hz background.\n"
        "Line excess = maximum power within 50 +/- 0.25 Hz relative to that background. "
        "A negative FIR value indicates a spectral hole, not recovered signal.",
        ha="left",
        va="bottom",
        fontsize=8.7,
        color="#4B5563",
        linespacing=1.35,
    )

    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def ensure_outputs_available(paths: list[Path], overwrite: bool) -> None:
    """Protect existing report artifacts unless explicit replacement is requested."""

    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        rendered = "\n".join(f"- {path}" for path in existing)
        raise FileExistsError(f"Outputs already exist; pass --overwrite to replace them:\n{rendered}")


def main() -> None:
    """Run filtering, PSD estimation, plotting, and provenance export."""

    args = parse_args()
    validate_args(args)
    mne.set_log_level("WARNING")

    output_dir = args.output_dir.resolve()
    output_png = output_dir / "line_noise_method_comparison.png"
    output_pdf = output_dir / "line_noise_method_comparison.pdf"
    output_csv = output_dir / "line_noise_method_comparison_values.csv"
    output_json = output_dir / "line_noise_method_comparison_metadata.json"
    outputs = [output_png, output_pdf, output_csv, output_json]
    ensure_outputs_available(outputs, args.overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load one shared source segment, then derive all method-specific copies from it.
    raw = load_segment(args.input.resolve(), args.duration)
    picks = mne.pick_types(raw.info, meg=True, ref_meg=False, exclude="bads")
    if len(picks) == 0:
        raise RuntimeError("No good MEG channels were found in the input file.")
    methods = build_method_copies(raw, picks, args.spectrum_fit_window)

    # Estimate spectra and retain both plotted summaries and per-channel diagnostics.
    summaries: dict[str, dict[str, Any]] = {}
    frequencies: np.ndarray | None = None
    for method in METHOD_ORDER:
        method_frequencies, psd = welch_psd(methods[method], picks)
        if frequencies is None:
            frequencies = method_frequencies
        elif not np.array_equal(frequencies, method_frequencies):
            raise RuntimeError("PSD frequency grids differ unexpectedly between methods.")

        median, q25, q75, per_channel, median_excess, p90_excess = summarize_spectrum(
            method_frequencies,
            psd,
        )
        summaries[method] = {
            "median": median,
            "q25": q25,
            "q75": q75,
            "line_excess_by_channel_db": per_channel,
            "median_excess_db": median_excess,
            "p90_excess_db": p90_excess,
        }

    assert frequencies is not None
    input_label = "Subject 01, block 01" if args.input.resolve() == DEFAULT_INPUT else args.input.name
    plot_comparison(
        frequencies,
        summaries,
        input_label,
        raw.n_times / float(raw.info["sfreq"]),
        output_png,
        output_pdf,
    )

    # Save the exact summary values used by panel C for review and reuse.
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("method", "median_50hz_excess_db", "p90_50hz_excess_db", "n_meg_channels"),
        )
        writer.writeheader()
        for method in METHOD_ORDER:
            writer.writerow(
                {
                    "method": method,
                    "median_50hz_excess_db": f"{summaries[method]['median_excess_db']:.6f}",
                    "p90_50hz_excess_db": f"{summaries[method]['p90_excess_db']:.6f}",
                    "n_meg_channels": len(picks),
                }
            )

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_fif": str(args.input.resolve()),
        "segment_start_seconds": 0.0,
        "segment_duration_requested_seconds": args.duration,
        "segment_samples": int(raw.n_times),
        "sampling_frequency_hz": float(raw.info["sfreq"]),
        "n_good_meg_channels": int(len(picks)),
        "line_frequency_hz": LINE_FREQUENCY_HZ,
        "spectrum_fit_window_seconds": args.spectrum_fit_window,
        "welch_n_fft": N_FFT,
        "welch_n_overlap": N_OVERLAP,
        "peak_half_width_hz": PEAK_HALF_WIDTH_HZ,
        "background_inner_offset_hz": BACKGROUND_INNER_HZ,
        "background_outer_offset_hz": BACKGROUND_OUTER_HZ,
        "mne_version": mne.__version__,
        "numpy_version": np.__version__,
        "methods": {
            method: {
                "median_50hz_excess_db": summaries[method]["median_excess_db"],
                "p90_50hz_excess_db": summaries[method]["p90_excess_db"],
            }
            for method in METHOD_ORDER
        },
        "interpretation_limit": (
            "Descriptive comparison of one subject/block/segment; not a ground-truth "
            "signal-recovery benchmark or participant-level inference."
        ),
    }
    output_json.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote: {output_png}")
    print(f"Wrote: {output_pdf}")
    print(f"Wrote: {output_csv}")
    for method in METHOD_ORDER:
        print(f"{method}: {summaries[method]['median_excess_db']:+.3f} dB")


if __name__ == "__main__":
    main()
