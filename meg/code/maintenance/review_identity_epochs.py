#!/usr/bin/env python3
"""Render native MEG evidence for every occurrence of one condition/identity.

Usage::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
      meg/code/maintenance/review_identity_epochs.py --subject 21 --identity 31 \
      --condition occluded_deviant --workers 4 \
      --output /tmp/p21-identity31-review

Inputs are immutable annotated SSS FIFs, the selected processing run's trial
decisions, and its recorded epoch settings. Outputs are per-acquisition PNGs,
an HTML contact sheet and a JSON evidence table in a new review directory.
No annotations, thresholds, trial decisions or scientific derivatives change.

For each epoch, read all MEG channels and reproduce spectrum-fit separately
over the COMPLETE native BAD/EDGE-delimited stretches touching the display.
This avoids a short-window filtering approximation. Then compute the same
linear-detrended peak-to-peak check on the exact full epoch, excluding marked
bad sensors. This utility currently supports inline_no_ica, spectrum-fit-only
epochs with baseline=None and detrend=1; other operators fail explicitly.
Acquisitions run in independent processes with bounded workers. Numerical
threads must be limited by the caller. Rendering joins before the HTML index.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import html
import json
import multiprocessing
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
from mne.annotations import _annotations_starts_stops
import numpy as np
from scipy.signal import detrend

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "preprocessing_basic"))
from continuous_cleaning import annotation_rows, digest, file_stamp

FRONTAL = ["MEG0311", "MEG0341", "MEG0511", "MEG0911", "MEG0921", "MEG1211"]


# -----------------------------------------------------------------------------
# Reconstruct exactly the existing operators in a disposable native slice.
# -----------------------------------------------------------------------------

def review_run(task):
    """Return evidence for one acquisition and write its review figures only."""
    source, trials, report, derivatives, output = task
    path = Path(derivatives) / "preprocessing" / report["subject_bids_id"] / "ses-01/sss" / (Path(source).stem + "_ann.fif")
    expected = next(x for x in report["preparation_contract"]["definition"]["preparation"]["inputs"] if x["source_filename"] == source)
    print(f"Reading {source}: {len(trials)} trial(s)", flush=True)
    raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
    if digest(annotation_rows(raw)) != expected["annotations_sha256"]:
        raise RuntimeError(f"Native annotations changed: {source}")
    if [file_stamp(p) for p in raw.filenames] != expected["parts"]:
        raise RuntimeError(f"Native input changed: {source}")
    sfreq = float(raw.info["sfreq"])
    picks = mne.pick_types(raw.info, meg=True, exclude=[])
    names = [raw.ch_names[i] for i in picks]
    types = [raw.get_channel_types()[i] for i in picks]
    starts, stops = _annotations_starts_stops(raw, ("bad", "edge"), invert=True)
    results = []
    for trial in trials:
        anchor = int(float(trial["epoch_event_sample_raw"]))
        left = max(0, anchor - raw.first_samp - round(2.5 * sfreq))
        stop = min(raw.n_times, anchor - raw.first_samp + round(4.5 * sfreq) + 1)
        data = raw.get_data(picks=picks, start=left, stop=stop)
        # Work on complete good stretches: notch fitting depends on stretch
        # boundaries, including its overlap/add placement and final short tail.
        for a, b in zip(starts, stops):
            if a >= stop or b <= left or b - a < 9:
                continue
            segment = raw.get_data(picks=picks, start=int(a), stop=int(b))
            mne.filter.notch_filter(segment, sfreq, [50., 100.], method="spectrum_fit",
                filter_length="10s", copy=False, verbose="ERROR")
            lo, hi = max(int(a), left), min(int(b), stop)
            data[:, lo-left:hi-left] = segment[:, lo-int(a):hi-int(a)]
        epoch_a = anchor - raw.first_samp + round(report["epoch_window_sec"][0] * sfreq)
        epoch_b = anchor - raw.first_samp + round(report["epoch_window_sec"][1] * sfreq) + 1
        epoch = data[:, epoch_a-left:epoch_b-left]
        residual = detrend(epoch, axis=-1, type="linear")
        ptp = np.ptp(residual, axis=-1)
        exceeded = [name for name, kind, value in zip(names, types, ptp)
            if name not in raw.info["bads"] and value > report["amplitude_reject"][kind]]
        recorded = [n for n in trial["mne_drop_reasons"].split(",") if n.startswith("MEG")]
        if recorded and set(recorded) != set(exceeded):
            raise RuntimeError(f"Reproduced amplitude reasons differ: {source}, {anchor}: {exceeded}")
        times = (np.arange(left, stop) + raw.first_samp - anchor) / sfreq
        epoch_times = np.arange(round(report["epoch_window_sec"][0]*sfreq),
            round(report["epoch_window_sec"][1]*sfreq)+1) / sfreq
        ann = [{"onset_relative_sec": float(a["onset"] - anchor/sfreq),
                "duration_sec": float(a["duration"]), "label": str(a["description"])}
            for a in raw.annotations if a["onset"] < (raw.first_samp+stop)/sfreq
            and a["onset"]+a["duration"] > (raw.first_samp+left)/sfreq
            and str(a["description"]).lower().startswith("bad")]

        # Left: native continuous context, centered only for legible offsets.
        # Right: the exact epoch operator; centering the min/max midpoint makes
        # the +/-2 pT band a correct display of a 4 pT peak-to-peak criterion.
        fig, axes = plt.subplots(len(FRONTAL), 2, figsize=(15, 12), sharex="col")
        for row, name in enumerate(FRONTAL):
            index = names.index(name)
            y = data[index] * 1e12
            y = y - np.median(y[epoch_a-left:epoch_b-left])
            axes[row, 0].plot(times, y, color="#26384a", lw=.65)
            axes[row, 0].axvspan(*report["epoch_window_sec"], color="#7db7d8", alpha=.16)
            z = residual[index] * 1e12
            z -= .5 * (z.max()+z.min())
            axes[row, 1].axhspan(-2, 2, color="#a8d5ba", alpha=.2)
            axes[row, 1].plot(epoch_times, z, color="#b34b22" if ptp[index] > 4e-12 else "#26384a", lw=.7)
            axes[row, 1].axhline(-2, color="#647669", ls="--", lw=.7)
            axes[row, 1].axhline(2, color="#647669", ls="--", lw=.7)
            axes[row, 1].set_title(f"{name}: peak-to-peak {ptp[index]*1e12:.2f} pT", fontsize=10, loc="left")
            axes[row, 0].set_ylabel(name + "\npT", fontsize=9)
            axes[row, 1].set_ylabel("pT", fontsize=9)
            for col in (0, 1):
                ax = axes[row, col]
                for a in ann:
                    ax.axvspan(a["onset_relative_sec"], a["onset_relative_sec"]+a["duration_sec"], color="#db6b74", alpha=.24)
                ax.axvline(0, color="#426c8e", lw=.8)
                ax.grid(alpha=.12)
                ax.tick_params(labelsize=9)
                ax.set_xlim((-2.5, 4.5) if col == 0 else report["epoch_window_sec"])
        reason = trial["mne_drop_reasons"] or trial.get("analysis_exclusion_reasons", "") or "retained"
        run = int(trial["acquisition_run"])
        fig.suptitle(f"P{int(report['subject_bids_id'].split('-')[1])} | {trial['condition_name']} | identity {int(float(trial['trial_sequence_identity']))} | acquisition {run}\n"
            f"Behavioral run {trial['behavioral_run']}, trial {trial['behavioral_trial']} | recorded rejection: {reason}", fontsize=12)
        fig.text(.26, .915, "Continuous signal: blue = full epoch; red = BAD annotation", ha="center", fontsize=10)
        fig.text(.75, .915, "Exact epoch after linear detrend; dashed band spans 4 pT", ha="center", fontsize=10)
        axes[-1, 0].set_xlabel("Seconds from corrected stimulus onset")
        axes[-1, 1].set_xlabel("Seconds from corrected stimulus onset")
        fig.tight_layout(rect=(0, .015, 1, .905))
        image = Path(output) / f"acquisition-{run:02d}_sample-{anchor}.png"
        fig.savefig(image, dpi=145); plt.close(fig)
        results.append({"trial_key": trial["trial_key"], "acquisition_run": run,
            "behavioral_run": trial["behavioral_run"], "behavioral_trial": trial["behavioral_trial"],
            "epoch_event_sample_raw": anchor, "recorded_rejection": reason,
            "reproduced_amplitude_channels": exceeded,
            "ptp_pT_by_magnetometer": {n: float(v*1e12) for n,k,v in zip(names,types,ptp) if k == "mag"},
            "annotation_context": ann, "figure": image.name})
        print(f"Finished {source}, anchor {anchor}: {len(exceeded)} amplitude channels", flush=True)
    raw.close()
    return results


# -----------------------------------------------------------------------------
# Select exact recorded trials, join processes, and publish a review-only index.
# -----------------------------------------------------------------------------

def main():
    """Select one identity and publish review evidence without changing decisions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", type=int, required=True)
    parser.add_argument("--identity", type=int, required=True)
    parser.add_argument("--condition", default="occluded_deviant")
    parser.add_argument("--processing-run", default="noica-20261117-v2")
    parser.add_argument("--derivatives", type=Path, default=Path("/Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1: parser.error("--workers must be positive")
    root = args.derivatives / "processing_runs" / args.processing_run
    subject = f"sub-{args.subject:02d}"
    report = json.loads((root / "reports" / subject / "ses-01/04_epoching_report.json").read_text())
    prep = report["preparation_contract"]["definition"]["preparation"]
    if prep["mode"] != "inline_no_ica" or report["baseline_sec"] is not None or report["detrend"] != 1:
        raise RuntimeError("Review supports the recorded no-ICA/no-baseline/linear-detrend variant only")
    if any(prep["continuous_filter"].get(k) is not None for k in ["highpass_hz", "lowpass_hz", "resample_sfreq_hz"]):
        raise RuntimeError("Unsupported continuous operators")
    with (root / "analysis" / subject / f"ses-01/meg/{subject}_ses-01_task-Attention_run-all_trial_decisions.csv").open() as f:
        trials = [r for r in csv.DictReader(f) if r["condition_name"] == args.condition
            and r["trial_sequence_identity"] and int(float(r["trial_sequence_identity"])) == args.identity
            and r["epoch_event_sample_raw"]]
    if not trials: raise RuntimeError("No trials with accepted corrected onset samples")
    args.output.mkdir(parents=True, exist_ok=False)
    groups = {name: [t for t in trials if t["source_fif"] == name] for name in sorted({t["source_fif"] for t in trials})}
    results = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(review_run, (source, rows, report, str(args.derivatives), str(args.output))) for source, rows in groups.items()]
        for future in as_completed(futures): results.extend(future.result())
    results.sort(key=lambda r: (r["acquisition_run"], r["epoch_event_sample_raw"]))
    payload = {"purpose": "review_only_no_decisions_changed", "subject": subject, "identity": args.identity,
        "condition": args.condition, "processing_run": args.processing_run, "epoch_window_sec": report["epoch_window_sec"],
        "amplitude_reject": report["amplitude_reject"], "trials": results}
    (args.output / "evidence.json").write_text(json.dumps(payload, indent=2))
    cards = "\n".join(f'<h2>Acquisition {r["acquisition_run"]}, behavioral run {r["behavioral_run"]}, trial {r["behavioral_trial"]}</h2><img src="{r["figure"]}" alt="Native and detrended frontal signals">' for r in results)
    (args.output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>P21 identity review</title><style>body{font:17px system-ui;max-width:1500px;margin:32px auto;padding:0 20px}img{width:100%}p{max-width:1050px;line-height:1.5}</style>'
        f'<h1>{html.escape(subject)}: {html.escape(args.condition)}, identity {args.identity}</h1>'
        '<p>Review only: original signals, annotations and trial decisions are unchanged. Blue shading marks the full epoch. Red marks BAD annotations. The left column shows continuous line-cleaned data with a constant visual offset. The right shows the exact epoch after linear detrending; its min/max midpoint is centered for display, so the dashed +/-2 pT band spans the 4 pT peak-to-peak threshold. All good MEG sensors were checked, including sensors not displayed.</p>' + cards)
    print(f"Review saved: {args.output / 'index.html'}", flush=True)


if __name__ == "__main__":
    main()
