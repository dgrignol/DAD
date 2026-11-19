#!/usr/bin/env python3
"""Audit P24 acquisition differences and propose visual-onset timing.

Purpose
-------
Deferred reference: this standalone adapter is not part of Stage 1. P24
timing must be reviewed in the later cohort-wide visual-timing stage before
this approach is adopted or new timing derivatives are generated.

P24 runs 1-4 were recorded with a 0.1 Hz online high-pass and without the
photodiode channel. Run 5 has MISC008 and corrected acquisition settings.
This separate adapter measures run-5 trial-onset STI101-to-MISC008 lag at 50%
of local pulse amplitude. It labels earlier-run offsets as estimates derived
solely from accepted run-5 measurements. Nothing here edits FIFs, epochs, or
the generic trigger events; the calibration requires human review.

Usage
-----
Reference invocation for later visual-timing development (not a Stage-1
execution instruction)::

    python meg/code/preprocessing_basic/p24_acquisition_adapter.py \
      --config meg/params/pipeline_default.yaml \
      --site-config meg/params/pipeline_cluster.yaml \
      --subject-config meg/params/sub-24.yaml

Inputs are the read-only P24 FIF headers and run-5 MISC008 samples plus the
Stage-02 trigger CSV. Outputs are a CSV and JSON in the P24 reports directory.
The per-event CSV preserves the digital sample and the source acquisition run.
For runs 1-4, proposed visual samples are explicitly marked estimated. The
0.1 Hz online high-pass cannot be undone and is reported as an acquisition
difference for the later filtering/analysis decision.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import mne
import numpy as np

from pipeline_config import assert_output_path_outside_raw, load_and_resolve_config


# Trial-onset markers are the only events used to estimate the lag. Warm-up and
# mid-trial markers can have different pulse timing and are excluded.
TRIAL_CODES = {71, 81, 91}
THRESHOLDS = (0.10, 0.25, 0.50, 0.75, 0.90)


def measure_rising_edge(signal: np.ndarray, sample: int, sfreq: float) -> dict[str, Any]:
    """Measure one event-local light rise and return sensitivity/QC values.

    `sample` is zero based relative to the first FIF sample. A short pre-edge
    baseline avoids the earlier high plateau of the four-frame trial pulse.
    The voltage and timing limits follow the existing photodiode audit and are
    review flags, not a claim that every accepted edge is physiologically exact.
    """
    pre = max(8, round(0.008 * sfreq))
    after = max(60, round(0.060 * sfreq))
    if sample < pre or sample + after >= signal.size:
        return {"qc_status": "flag", "qc_reason": "event_near_file_boundary"}
    baseline = float(np.median(signal[sample - pre:sample]))
    window = signal[sample:sample + after]
    peak = float(np.quantile(window[:max(30, round(0.030 * sfreq))], 0.95))
    amplitude = peak - baseline
    result: dict[str, Any] = {"baseline_v": baseline, "peak_v": peak,
                              "amplitude_v": amplitude}
    for fraction in THRESHOLDS:
        key = f"lag_{round(fraction * 100):02d}_ms"
        hits = np.flatnonzero(window >= baseline + fraction * amplitude)
        result[key] = float(hits[0] * 1000 / sfreq) if len(hits) else math.nan
    lags = np.asarray([result[f"lag_{round(f * 100):02d}_ms"]
                       for f in THRESHOLDS], dtype=float)
    reasons = []
    if amplitude < 0.005:
        reasons.append("low_amplitude")
    if not np.isfinite(lags).all():
        reasons.append("missing_onset")
    elif np.ptp(lags) > 3:
        reasons.append("threshold_sensitive")
    if np.isfinite(result["lag_50_ms"]) and not 0 <= result["lag_50_ms"] <= 50:
        reasons.append("lag_outside_window")
    result["qc_status"] = "pass" if not reasons else "flag"
    result["qc_reason"] = ";".join(reasons)
    return result


def read_trigger_events(path: Path) -> list[dict[str, str]]:
    """Load Stage-02 events while preserving its digital sample and metadata."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_timing_proposal(events: list[dict[str, str]], signal: np.ndarray,
                          sfreq: float, first_samp: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Measure run 5, then estimate earlier runs from its accepted median.

    The accepted run-5 lag is applied only to aligned trial-onset codes. A
    flagged run-5 event remains uncorrected and must be inspected separately.
    This returns a proposal, not events ready for epoching.
    """
    rows: list[dict[str, Any]] = []
    for event in events:
        if int(event["trigger_code"]) not in TRIAL_CODES:
            continue
        run = int(event["acquisition_run"])
        sample = int(event["raw_sample_frame"])
        row: dict[str, Any] = {
            "event_index": event["event_index"],
            "acquisition_run": run,
            "behavioral_block": event.get("behavioral_block", ""),
            "trigger_code": int(event["trigger_code"]),
            "digital_sample": sample,
            "source_fif": event["source_fif"],
            "reported_alignment_status": event.get("reported_alignment_status", ""),
            "timing_status": "review_required",
            "proposed_visual_sample": "",
        }
        if run == 5:
            row.update(measure_rising_edge(signal, sample - first_samp, sfreq))
            row["timing_status"] = "measured_run5" if row["qc_status"] == "pass" else "flagged_run5"
            if row["qc_status"] == "pass":
                row["proposed_visual_sample"] = sample + round(row["lag_50_ms"] * sfreq / 1000)
        else:
            row.update({"qc_status": "not_recorded", "qc_reason": "MISC008_absent"})
        rows.append(row)
    run5_aligned = [row for row in rows if row["acquisition_run"] == 5
                    and row["reported_alignment_status"] in {"MATCH", "MATCH_WARNING"}]
    passing = [float(row["lag_50_ms"]) for row in rows
               if row["acquisition_run"] == 5 and row["timing_status"] == "measured_run5"
               and row["reported_alignment_status"] in {"MATCH", "MATCH_WARNING"}]
    if len(passing) < 20 or len(passing) < 0.8 * len(run5_aligned):
        raise RuntimeError(
            f"Only {len(passing)}/{len(run5_aligned)} aligned run-5 diode "
            "trial onsets passed QC; no calibration proposed."
        )
    lag_ms = float(np.median(passing))
    lag_samples = round(lag_ms * sfreq / 1000)
    for row in rows:
        if row["acquisition_run"] < 5 and row["reported_alignment_status"] in {"MATCH", "MATCH_WARNING"}:
            row["timing_status"] = "estimated_from_run5_review_required"
            row["lag_50_ms"] = lag_ms
            row["proposed_visual_sample"] = row["digital_sample"] + lag_samples
        elif row["acquisition_run"] < 5:
            row["timing_status"] = "unaligned_review_required"
    summary = {
        "calibration_source": "P24 run 5 only; aligned trial-onset MISC008 rises passing QC",
        "accepted_run5_events": len(passing),
        "aligned_run5_trial_events": len(run5_aligned),
        "run5_median_lag_ms": lag_ms,
        "run5_lag_p05_ms": float(np.quantile(passing, 0.05)),
        "run5_lag_p95_ms": float(np.quantile(passing, 0.95)),
        "run5_lag_mad_ms": float(np.median(np.abs(np.asarray(passing) - lag_ms))),
        "early_runs_status": "estimated; no photodiode was recorded; review before application",
        "proposal_status": "review_required",
    }
    return rows, summary


def main() -> int:
    """Read P24 inputs, guard output location, and write the timing review."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config")
    parser.add_argument("--subject-config", required=True)
    args = parser.parse_args()
    cfg, _, _, _ = load_and_resolve_config(
        args.config, args.subject_config, site_config_path=args.site_config
    )
    if cfg["subject"]["id"] != "24":
        raise ValueError("This acquisition adapter is restricted to P24.")
    resolved = cfg["resolved_paths"]
    reports = Path(resolved["derivatives_reports_dir"])
    event_path = Path(resolved["derivatives_preprocessing_dir"]) / "events" / "02_trigger_events.csv"
    output_csv = reports / "p24_visual_timing_proposal.csv"
    output_json = reports / "p24_acquisition_audit.json"
    for path in (output_csv, output_json):
        assert_output_path_outside_raw(path, resolved["data_root"])
        if path.exists():
            raise FileExistsError(f"Preserving existing review derivative: {path}")
    if not event_path.is_file():
        raise FileNotFoundError(f"Run Stage 02 first: {event_path}")

    # FIF headers are read lazily, preserving the acquisition filter facts. The
    # run-5 photodiode is the only continuous channel loaded into memory.
    fif_paths = sorted(Path(resolved["subject_data_dir"]).glob("meg/sub24_run*.fif"))
    if len(fif_paths) != 5:
        raise RuntimeError(f"Expected P24's five acquisition FIFs; found {len(fif_paths)}")
    headers = []
    run5 = None
    for index, path in enumerate(fif_paths, start=1):
        raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
        has_diode = "MISC008" in raw.ch_names
        headers.append({"run": index, "fif": str(path),
                        "online_highpass_hz": float(raw.info["highpass"]),
                        "online_lowpass_hz": float(raw.info["lowpass"]),
                        "sfreq_hz": float(raw.info["sfreq"]),
                        "n_channels": int(raw.info["nchan"]),
                        "misc_channels": [name for name, kind in zip(
                            raw.ch_names, raw.get_channel_types(), strict=True
                        ) if kind == "misc"],
                        "misc008_present": has_diode})
        if index < 5 and has_diode:
            raise RuntimeError(f"Unexpected MISC008 in early run: {path}")
        if index == 5:
            if not has_diode:
                raise RuntimeError(f"Run 5 has no MISC008: {path}")
            run5 = raw
    assert run5 is not None
    signal = run5.get_data(picks=["MISC008"])[0]
    events = read_trigger_events(event_path)
    rows, summary = build_timing_proposal(
        events, signal, float(run5.info["sfreq"]), int(run5.first_samp)
    )
    summary["acquisition_headers"] = headers
    summary["online_filter_warning"] = (
        "Runs 1-4 have an irreversible 0.1 Hz acquisition high-pass; "
        "do not infer 0.01 Hz information or silently equate low-frequency spectra."
    )
    summary["input_events_csv"] = str(event_path)
    summary["n_proposed_trial_events"] = len(rows)
    reports.mkdir(parents=True, exist_ok=True)
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    output_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"P24 timing proposal requires review: {output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
