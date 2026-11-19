#!/usr/bin/env python3
"""Build explained, offline Stage 3 participant review reports from observations.

Usage from the DAD root to present an existing immutable measurement bundle::

    python meg/code/preprocessing_basic/visual_timing_report.py \
      --bundle /path/to/derivatives/visual_timing/VERSION/sub-11/ses-01 \
      --output /path/to/derivatives/reports/visual_timing/REPORT_VERSION/sub-11/ses-01 \
      --raw-root /path/to/Attention/data

The measurement CLI calls ``build_report_artifacts`` automatically for new
bundles. The report-only CLI verifies every source output hash, reads existing
tables and native waveform caches, and publishes a NEW presentation directory.
It never opens a raw FIF, estimates a new lag, changes a timing assignment, or
authors scientific acceptance. Existing measurements and reports are read-only.

Data flow: verified measured/estimated/fallback proposals -> descriptive group
statistics and sanity checks -> categorical plots -> explanations, findings,
decision examples and exact evidence tables -> self-contained report.html plus
Markdown, PNG and CSV companions. HTML embeds figures and downloadable tables;
it needs neither a web server nor external scripts, styles or network access.

Trigger codes are factors everywhere. Bars show measured mean +/- sample SD;
medians and observed percentiles accompany them. Donor location confidence
intervals remain distinct from observed dispersion and unknown transfer error.
Pilot evidence is explicitly incomplete. A descriptive group center is never
silently substituted for individual measured onsets or a restricted estimate.

Within-trial review retains the original interval/bar figure and adds a violin
with every measured pair, residual versus expected duration, residual versus
recording time, and categorical endpoint-pair distributions. Exact endpoint
joins and source/frame/sample arithmetic are checked before these presentation
views are built. Several interval observations may share one trial/onset; their
sample SD is descriptive spread, not independent-trial or physical uncertainty.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter, defaultdict
import hashlib
import html
import io
import json
from pathlib import Path

import numpy as np

from downstream_stage_utils import csv_bool
from preprocessing_stage_utils import utc_now
from trial_preparation import csv_text, file_hash, read_csv
from visual_timing_utils import DEFAULT_METHOD, analysis_eligible, donor_allowed, round_half_away, verify_bundle


REPORT_SCHEMA_VERSION = "1.1"
TRIGGER_NAMES = {
    71: "Always-visible onset", 81: "Occluded nondeviant onset",
    91: "Occluded deviant onset", 111: "Occlusion start",
    112: "Occlusion complete", 114: "Occlusion end start",
    115: "Occlusion end complete", 116: "Catch disappear",
    117: "Catch reappear", 118: "Question onset", 153: "Warm-up (diagnostic)",
}
TRIGGER_COLORS = dict(zip(TRIGGER_NAMES, [
    "#2563a1", "#c46a1b", "#627b32", "#8653a3", "#b64876", "#267c7b",
    "#9b7333", "#526d9c", "#8c5849", "#52525b", "#a48412",
]))
STATUS_COLORS = {"measured": "#2563a1", "flagged": "#bd6b22", "unusable": "#8653a3",
                 "confirmed_unusable": "#8653a3", "not_recorded": "#777777"}
INTERVAL_PAIR_ORDER = [(71,116),(71,117),(81,111),(81,112),(81,114),(81,115),
                       (91,111),(91,112),(91,114),(91,115),
                       (111,112),(112,114),(114,115),(116,117)]
INTERVAL_PAIR_COLORS = ["#2563a1","#c46a1b","#627b32","#8653a3","#b64876","#267c7b","#9b7333",
                        "#536f99","#8c5849","#52525b","#a48412","#499883","#d38e73","#ab8ebd"]
RECORDING_COLORS = ["#2563a1","#c46a1b","#627b32","#8653a3","#267c7b","#9b7333"]
NUMERIC_FIELDS = {
    "trigger_code", "digital_sample_raw", "first_samp", "sfreq_hz", "lag_ms",
    "amplitude_v", "baseline_noise_v", "threshold_span_ms", "polarity",
    "pulse_duration_ms", "expected_pulse_duration_ms", "baseline_normalized_v",
    "raw_sample_frame", "condition_code", "recording_first_sample_raw",
    "recording_end_sample_raw_exclusive", "proposed_sample_raw", "applied_shift_samples",
    "trial_interval_end_raw_exclusive", "saved_trial_onset_guard_retry_count",
    "diode_edge_sample_raw", "edge_10_offset_samples", "edge_25_offset_samples",
    "edge_50_offset_samples", "edge_75_offset_samples", "edge_90_offset_samples",
    "observed_diode_interval_ms", "expected_frame_interval_ms", "interval_error_ms",
}


# -----------------------------------------------------------------------------
# Normalize CSV evidence and compute transparent descriptive statistics
# -----------------------------------------------------------------------------


def number(value):
    """Return a finite numeric value, preserving missing/undefined values."""
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def normalize(rows):
    """Copy table rows; numeric conversion never changes the source tables."""
    return [{k: number(v) if k in NUMERIC_FIELDS else v for k, v in row.items()} for row in rows]


def fmt(value, digits=2):
    """Use an explicit unavailable marker rather than inventing a zero."""
    value = number(value)
    return "unavailable" if value is None else f"{value:.{digits}f}"


def stats(values):
    """Describe observations; sample SD needs at least two finite observations."""
    a = np.asarray([n for value in values if (n := number(value)) is not None], float)
    if not len(a):
        return dict(n=0, mean=None, median=None, sd=None, p05=None, p95=None, minimum=None, maximum=None)
    return dict(n=len(a), mean=float(a.mean()), median=float(np.median(a)),
                sd=float(a.std(ddof=1)) if len(a) > 1 else None,
                p05=float(np.quantile(a, .05)), p95=float(np.quantile(a, .95)),
                minimum=float(a.min()), maximum=float(a.max()))


def code_label(code):
    """Keep the recorded numeric code as a categorical label with its meaning."""
    code = int(code)
    return f"{code} {TRIGGER_NAMES.get(code, 'Other recorded marker')}"


def strata(rows):
    """Never pool distinct source FIFs or saved runtimes into one distribution."""
    return sorted({(r["source_fif"], r.get("runtime", "unknown")) for r in rows})


def lag_statistics(rows):
    """Per recording/runtime/code summaries separate usable and flagged edges."""
    groups = defaultdict(list)
    for row in rows:
        groups[(row["source_fif"], row.get("runtime", "unknown"), int(row["trigger_code"]))].append(row)
    output = []
    for (source, runtime, code), group in sorted(groups.items()):
        measured = [r for r in group if r["measurement_status"] == "measured"]
        summary = stats(r.get("lag_ms") for r in measured)
        output.append(dict(source_fif=source, runtime=runtime, trigger_code=code,
                           trigger_meaning=TRIGGER_NAMES.get(code, "Other marker"),
                           n_assessed=len(group), n_measured=summary.pop("n"),
                           n_flagged=sum(r["measurement_status"] == "flagged" for r in group),
                           **{f"{k}_lag_ms": v for k, v in summary.items()}))
    return output


def coverage_and_checks(maps, window):
    """Recompute original/proposed coverage and arithmetic without changing maps.

    Inclusive epoch endpoints are checked against exclusive recording ends.
    Unassigned timing and unknown uncertainty remain explicit. Ordering and
    collisions use eligible assigned trials within the same exact source FIF.
    """
    output, assigned = [], defaultdict(list)
    for row in maps:
        sf = row["sfreq_hz"]
        first, end = row["recording_first_sample_raw"], row["recording_end_sample_raw_exclusive"]
        original, proposed = row["raw_sample_frame"], row.get("proposed_sample_raw")
        before, after = [round_half_away(float(v) * sf) for v in window]
        original_inside = first <= original + before and original + after < end
        proposed_inside = None if proposed is None else first <= proposed + before and proposed + after < end
        arithmetic = None if proposed is None else bool(row.get("applied_shift_samples") is not None
                     and np.isclose(proposed-original, row["applied_shift_samples"], atol=1e-8, rtol=0))
        eligible = analysis_eligible(row)
        output.append(dict(trial_key=row["trial_key"], source_fif=row["source_fif"],
                           condition_code=row["condition_code"], eligible=eligible,
                           assignment_basis=row["assignment_basis"], original_inside=original_inside,
                           proposed_inside=proposed_inside, arithmetic_matches=arithmetic,
                           uncertainty_inside=row.get("uncertainty_window_inside_recording", ""),
                           overlaps_trial_end=row.get("epoch_overlaps_trial_interval_end", "")))
        if eligible and proposed is not None:
            assigned[row["source_fif"]].append((original, proposed, row["trial_key"]))
    collisions, inversions = [], []
    for source, group in assigned.items():
        seen = {}
        previous = None
        for original, proposed, key in sorted(group):
            if proposed in seen:
                collisions.append(dict(source_fif=source, trial_key=key, other_key=seen[proposed], sample=proposed))
            seen[proposed] = key
            if previous and proposed < previous[1]:
                inversions.append(dict(source_fif=source, trial_key=key, other_key=previous[2]))
            previous = (original, proposed, key)
    checks = dict(original_keys=len(maps), unique_original_keys=len({r["trial_key"] for r in maps}),
                  assigned_arithmetic_mismatches=sum(r["arithmetic_matches"] is False for r in output),
                  eligible_proposed_sample_collisions=collisions, eligible_ordering_inversions=inversions)
    return output, checks


def correction_statistics(maps):
    """Show the actual proposed shifts; descriptive centers are not corrections."""
    groups = defaultdict(list)
    for row in maps:
        if analysis_eligible(row):
            groups[(row["source_fif"], int(row["condition_code"]), row["assignment_basis"])].append(row)
    result = []
    for (source, code, basis), group in sorted(groups.items()):
        values = [r["applied_shift_samples"] * 1000 / r["sfreq_hz"] for r in group
                  if r.get("applied_shift_samples") is not None]
        s = stats(values)
        result.append(dict(source_fif=source, trigger_code=code, assignment_basis=basis,
                           n_eligible=len(group), n_assigned=s.pop("n"),
                           **{f"{k}_shift_ms": v for k, v in s.items()}))
    return result


# -----------------------------------------------------------------------------
# Reproducible static figures; categories never use numeric trigger distances
# -----------------------------------------------------------------------------


def pyplot():
    """Load the existing headless cluster renderer only when figures are needed."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def png(fig):
    """Export a standalone figure with the same content embedded in HTML."""
    data = io.BytesIO()
    fig.savefig(data, format="png", dpi=125, bbox_inches="tight", facecolor="white")
    pyplot().close(fig)
    return data.getvalue()


def empty_figure(message, title):
    """Unavailable evidence gets an explanation, never fabricated observations."""
    fig, ax = pyplot().subplots(figsize=(11, 2.3))
    ax.axis("off")
    ax.text(.02, .65, title, fontsize=14, weight="bold", transform=ax.transAxes)
    ax.text(.02, .28, message, fontsize=11, wrap=True, transform=ax.transAxes)
    return png(fig)


def lag_axis_limits(rows):
    """Focus on observed centers/dispersion without clipping any finite tail.

    The envelope includes each usable group's mean +/- two sample SD and all
    finite measured/flagged points, then adds a small native-grid-aware margin.
    These are display limits, not an outlier rule or an uncertainty interval.
    """
    values = [r["lag_ms"] for r in rows if r.get("lag_ms") is not None]
    if not values:return (-1,1)
    bounds = list(values)
    for group in lag_statistics(rows):
        if group["mean_lag_ms"] is not None:
            mean,sd = group["mean_lag_ms"],group["sd_lag_ms"] or 0
            bounds.extend([mean-2*sd,mean+2*sd])
    low,high = min(bounds),max(bounds)
    native = min(1000/r["sfreq_hz"] for r in rows if r.get("sfreq_hz"))
    margin = max((high-low)*.06,native*.25,.15)
    return low-margin,high+margin


def lag_figures(rows):
    """Compare trigger-factor distributions and mean +/- sample SD per stratum."""
    plt = pyplot()
    all_groups = strata(rows)
    groups = [(source,runtime) for source,runtime in all_groups if any(
        r["source_fif"]==source and r.get("runtime","unknown")==runtime and r.get("lag_ms") is not None for r in rows)]
    missing = [source for source,runtime in all_groups if (source,runtime) not in groups]
    if not groups:
        image = empty_figure("No finite lag is available in this scope. Consult the measured/failed counts and acquisition facts.", "No lag evidence")
        return {"delay_violin_with_dots.png": image, "lag_by_trigger.png": image}
    codes = sorted({int(r["trigger_code"]) for r in rows})
    limits = lag_axis_limits(rows)
    fig_v, axes_v = plt.subplots(len(groups), 1, figsize=(13, max(4, len(groups)*(1.6+.36*len(codes)))), squeeze=False, sharex=True)
    fig_b, axes_b = plt.subplots(len(groups), 1, figsize=(13, max(4, len(groups)*3.6)), squeeze=False, sharey=True)
    for (source, runtime), axv, axb in zip(groups, axes_v.flat, axes_b.flat):
        group = [r for r in rows if (r["source_fif"], r.get("runtime", "unknown")) == (source, runtime)]
        for position, code in enumerate(codes):
            subset = [r for r in group if int(r["trigger_code"]) == code]
            good = [r for r in subset if r["measurement_status"] == "measured" and r.get("lag_ms") is not None]
            values = [r["lag_ms"] for r in good]
            if len(values) >= 3 and np.ptp(values) > 0:
                violin = axv.violinplot([values], positions=[position], vert=False,
                                       widths=.75, showextrema=False, bw_method="scott")
                for body in violin["bodies"]:
                    body.set_facecolor(TRIGGER_COLORS.get(code, "#555555")); body.set_alpha(.23)
            jitter = np.linspace(-.18, .18, len(values)) if len(values)>1 else np.zeros(len(values))
            axv.scatter(values, position+jitter, s=12, color=TRIGGER_COLORS.get(code, "#555555"), alpha=.65)
            bad = [r["lag_ms"] for r in subset if r["measurement_status"] != "measured" and r.get("lag_ms") is not None]
            axv.scatter(bad, [position]*len(bad), s=30, marker="x", color="#222222")
            if values:
                s = stats(values)
                axv.plot([s["p05"], s["p95"]], [position, position], color="#222222", lw=1)
                axv.scatter([s["median"]], [position], marker="D", facecolors="none", edgecolors="#222222", zorder=5)
                axb.bar(position, s["mean"], width=.65, color=TRIGGER_COLORS.get(code, "#555555"), alpha=.75)
                if s["sd"] is not None:
                    axb.errorbar(position, s["mean"], yerr=s["sd"], fmt="none", ecolor="#222222", capsize=4)
                axb.scatter([position], [s["median"]], marker="D", facecolors="white", edgecolors="#222222", zorder=5)
                axb.annotate(f"n={len(values)}", (position, max(0,s["mean"])+(s["sd"] or 0)),
                             xytext=(0,5), textcoords="offset points", ha="center", fontsize=8)
            elif subset:
                axb.text(position, .04, "no usable\nedge", transform=axb.get_xaxis_transform(), ha="center", fontsize=8)
        labels = []
        for c in codes:
            subset = [r for r in group if int(r["trigger_code"])==c]
            n = sum(r["measurement_status"]=="measured" for r in subset)
            labels.append(f"{code_label(c)} (n={n}/{len(subset)})")
        axv.set(yticks=range(len(codes)), yticklabels=labels,xlim=limits,
                xlabel="Signed lag from STI101 to candidate 50% diode edge (ms)", title=f"{source} | {runtime}")
        axv.invert_yaxis(); axv.grid(axis="x", alpha=.18)
        axv.set_ylim(len(codes)-.5,-.5)
        if limits[0]<=0<=limits[1]:axv.axvline(0,color="#777777",lw=.7)
        if not any(r.get("lag_ms") is not None for r in group):
            axv.text(.03,.5,"No usable or finite candidate lag in this assessed sample",transform=axv.transAxes,fontsize=10)
        axb.set(xticks=range(len(codes)), xticklabels=[str(c) for c in codes],
                xlim=(-.6,len(codes)-.4),
                xlabel="Trigger code (categorical; see code meanings in the table)", ylabel="Measured lag mean +/- SD (ms)",
                title=f"{source} | {runtime}")
        axb.axhline(0, color="#555555", lw=.8); axb.grid(axis="y", alpha=.18); axb.margins(y=.25)
    fig_v.suptitle("Are trigger lags consistent, and where are their centers?\nOne coloured dot = one measured event; diamond = group median; line = observed 5-95%; x = flagged", fontsize=12)
    fig_b.suptitle("Lag by trigger category: bar = mean; whisker = sample SD; diamond = median\nSD is observed dispersion, not confidence in a correction; unavailable for n=1", fontsize=13)
    footer = "No finite lag (see QC/coverage): "+", ".join(missing) if missing else "Groups retain recording/runtime identity; no finite lag points are clipped."
    fig_v.text(.02,.008,footer,fontsize=8);fig_b.text(.02,.008,footer,fontsize=8)
    fig_v.tight_layout(rect=(0,.03,1,1)); fig_b.tight_layout(rect=(0,.03,1,1))
    return {"delay_violin_with_dots.png": png(fig_v), "lag_by_trigger.png": png(fig_b)}


def time_figures(rows, method):
    """Separate recordings and give each signal panel a run-specific y scale."""
    plt = pyplot()
    groups = strata(rows)
    if not groups:
        image = empty_figure("No assessed visual events.", "Timing and signal evidence")
        return {"lag_over_time.png": image, "signal_usability.png": image, "threshold_sensitivity.png": image}
    fig_t, axes_t = plt.subplots(len(groups), 1, figsize=(13, 3.1*len(groups)), squeeze=False, sharey=True)
    # Signal metrics have different ranges across recordings; keep each panel's
    # view limits local to its exact source/runtime rather than sharing by column.
    fig_s, axes_s = plt.subplots(len(groups), 3, figsize=(15, 3*len(groups)), squeeze=False)
    fig_q, axes_q = plt.subplots(len(groups), 2, figsize=(14, 3.1*len(groups)), squeeze=False, sharey="col")
    for index, (source, runtime) in enumerate(groups):
        group = [r for r in rows if (r["source_fif"], r.get("runtime", "unknown")) == (source, runtime)]
        at, signal, quality = axes_t[index,0], axes_s[index], axes_q[index]
        signal_minima, signal_maxima = [0.0]*3, [0.0]*3
        for code in sorted({int(r["trigger_code"]) for r in group}):
            subset = [r for r in group if int(r["trigger_code"]) == code]
            lagged = [r for r in subset if r.get("lag_ms") is not None]
            color = TRIGGER_COLORS.get(code, "#555555")
            for status in sorted({r["measurement_status"] for r in lagged}):
                g = [r for r in lagged if r["measurement_status"] == status]
                at.scatter([(r["digital_sample_raw"]-r["first_samp"])/r["sfreq_hz"] for r in g],
                           [r["lag_ms"] for r in g], s=16, color=color,
                           marker="o" if status=="measured" else "x", label=f"{code} {status}")
            for column, (axis, field) in enumerate(zip(signal, ("amplitude_v", "baseline_noise_v", "threshold_span_ms"))):
                g = [r for r in subset if r.get(field) is not None]
                scale = 1000 if field.endswith("_v") else 1
                y_values = [r[field]*scale for r in g]
                if y_values:
                    signal_minima[column] = min(signal_minima[column], min(y_values))
                    signal_maxima[column] = max(signal_maxima[column], max(y_values))
                axis.scatter([(r["digital_sample_raw"]-r["first_samp"])/r["sfreq_hz"] for r in g],
                             y_values, s=12, color=color, label=str(code))
            good = [r for r in lagged if r["measurement_status"]=="measured"]
            summaries = [stats(r.get(f"edge_{level}_offset_samples")*1000/r["sfreq_hz"]
                              for r in good if r.get(f"edge_{level}_offset_samples") is not None)
                         for level in (10,25,50,75,90)]
            points = [(i,s) for i,s in enumerate(summaries) if s["n"]]
            if points:
                quality[0].plot([i for i,s in points],[s["mean"] for i,s in points],marker="o",color=color,label=str(code))
                with_sd = [(i,s) for i,s in points if s["sd"] is not None]
                if with_sd:
                    quality[0].errorbar([i for i,s in with_sd],[s["mean"] for i,s in with_sd],
                                       yerr=[s["sd"] for i,s in with_sd],fmt="none",ecolor=color,capsize=3)
            g = [r for r in lagged if r.get("threshold_span_ms") is not None]
            quality[1].scatter([r["lag_ms"] for r in g], [r["threshold_span_ms"] for r in g], s=15, color=color, label=str(code))
        at.set(title=f"{source} | {runtime}", xlabel="Time since this recording's first sample (s)", ylabel="Signed lag (ms)")
        if any(r.get("lag_ms") is not None for r in group): at.legend(fontsize=8, ncol=4, loc="upper right")
        else: at.text(.03,.45,"No usable/candidate lag in this assessed sample", transform=at.transAxes)
        at.grid(alpha=.18)
        for column, (axis, label) in enumerate(zip(signal, ("Response amplitude (mV)", "Robust baseline noise (mV)", "10-90% edge span (ms)"))):
            axis.set(xlabel="Time within this recording (s)", ylabel=label); axis.grid(alpha=.18)
            if column == 0:
                reference = method["minimum_amplitude_v"]*1000
                signal_maxima[column] = max(signal_maxima[column], reference)
                signal_minima[column] = min(signal_minima[column], reference)
            # Add a small margin around this run's finite values (and the
            # amplitude reference line); preserve the zero floor for these
            # nonnegative signal metrics while still including any negatives.
            span = signal_maxima[column] - signal_minima[column]
            margin = .08*span if span > 0 else (.08*abs(signal_maxima[column]) if signal_maxima[column] else 0.0)
            lower = min(0.0, signal_minima[column])
            upper = signal_maxima[column] + margin
            if upper <= lower:
                upper = lower + 1.0
            axis.set_ylim(lower, upper)
        signal[0].set_title(f"{source} | {runtime}")
        signal[0].axhline(method["minimum_amplitude_v"]*1000, color="#777777", ls="--", lw=.8)
        if any(r.get("amplitude_v") is not None for r in group): signal[0].legend(title="Trigger", fontsize=7, ncol=4)
        quality[0].set(title=source, xticks=range(5), xticklabels=["10%","25%","50%","75%","90%"],
                       xlabel="Response fraction defining the candidate edge", ylabel="Lag mean +/- sample SD (ms)")
        quality[1].set(xlabel="Candidate 50% lag (ms)", ylabel="10-90% span (ms)")
        quality[1].axhline(method["maximum_sensitivity_ms"], color="#777777", ls="--", lw=.8)
        for axis in quality:
            axis.grid(alpha=.18)
            if any(r.get("lag_ms") is not None for r in group): axis.legend(title="Trigger", fontsize=7, ncol=4)
    fig_t.suptitle("Inspect within-recording lag stability and genuine timing tails", fontsize=14)
    fig_s.suptitle("Signal usability: amplitude, baseline noise and edge sensitivity by recording", fontsize=14)
    fig_q.suptitle("Detector sensitivity: threshold-specific means +/- SD and individual edge spans", fontsize=14)
    for fig in (fig_t,fig_s,fig_q): fig.tight_layout(rect=(0,0,1,.97))
    return {"lag_over_time.png":png(fig_t), "signal_usability.png":png(fig_s), "threshold_sensitivity.png":png(fig_q)}


def measurement_coverage_figure(rows):
    """Show detection failures with the assessed source/code denominator visible."""
    if not rows:return empty_figure("No visual observations were assessed.","Measurement coverage")
    groups = strata(rows)
    codes = sorted({int(r["trigger_code"]) for r in rows})
    statuses = sorted({r["measurement_status"] for r in rows})
    fig,axes = pyplot().subplots(len(groups),1,figsize=(13,2.8*len(groups)),squeeze=False)
    for (source,runtime),axis in zip(groups,axes.flat):
        bottom = np.zeros(len(codes))
        for status in statuses:
            counts = [sum(r["source_fif"]==source and r.get("runtime","unknown")==runtime
                          and int(r["trigger_code"])==code and r["measurement_status"]==status for r in rows) for code in codes]
            axis.bar(range(len(codes)),counts,bottom=bottom,label=status,color=STATUS_COLORS.get(status,"#a1a1aa"))
            bottom += counts
        axis.set(title=f"{source} | {runtime}",xticks=range(len(codes)),xticklabels=[str(c) for c in codes],
                 xlabel="Trigger code (categorical)",ylabel="Assessed event count")
        axis.legend(fontsize=8,ncol=3);axis.set_ylim(bottom=0);axis.margins(y=.25)
    fig.suptitle("Which assessed events have usable measured timing? Counts retain failed/missing observations",fontsize=13)
    fig.tight_layout(rect=(0,0,1,.96))
    return png(fig)


def select_waveforms(rows, waveforms):
    """Choose reproducible per-recording controls, families, tails and QC cases.

    All windows are already in memory/cache. At most 16 per recording are
    displayed and cached; no extra raw reads occur. Pilots retain every window.
    The index and report state the displayed/available counts explicitly.
    """
    candidates = [r for r in rows if r["event_key"] in waveforms]
    if len(candidates) <= 12:
        return candidates
    output = []
    for source in sorted({r["source_fif"] for r in candidates}):
        group = [r for r in candidates if r["source_fif"]==source]
        chosen = []
        def add(row):
            if row not in chosen and len(chosen)<16: chosen.append(row)
        # Reserve one control per factor first, so later code numbers (including
        # question/warm-up families) cannot be crowded out by earlier codes.
        for code in sorted({int(r["trigger_code"]) for r in group}):
            subset = [r for r in group if int(r["trigger_code"])==code]
            add(subset[len(subset)//2])
        for row in [r for r in group if csv_bool(r.get("post_warmup_context",False))][:1]:add(row)
        for row in sorted(group, key=lambda r:(r["measurement_status"]=="measured", r["digital_sample_raw"]))[:3]: add(row)
        for row in sorted(group, key=lambda r:r.get("lag_ms") if r.get("lag_ms") is not None else -999)[-2:]: add(row)
        output.extend(chosen)
    return output


def waveform_figures(rows, waveforms):
    """Annotate cached native voltages, observed STI101 and baseline/edge choice."""
    plt = pyplot()
    selected = select_waveforms(rows, waveforms)
    files, index, arrays, pages = {}, [], {}, []
    by_key = {r["event_key"]:r for r in selected}
    for page, source in enumerate(sorted({r["source_fif"] for r in selected}),1):
        group = [r for r in selected if r["source_fif"]==source]
        fig, axes = plt.subplots(int(np.ceil(len(group)/3)),3,figsize=(15,3.1*int(np.ceil(len(group)/3))),squeeze=False)
        for ax,row in zip(axes.flat,group):
            values,anchor,sf = waveforms[row["event_key"]]
            t = (np.arange(values.shape[1])-anchor)*1000/sf
            ax.plot(t,values[1]*1000,color="#2563a1",lw=1)
            ax.axvline(0,color="#333333",ls="--",lw=1,label="Original STI101")
            if row.get("lag_ms") is not None: ax.axvline(row["lag_ms"],color="#c46a1b",ls=":",label="50% edge")
            polarity = row.get("polarity")
            if polarity and row.get("baseline_normalized_v") is not None and row.get("amplitude_v") is not None:
                level = (row["baseline_normalized_v"]+.5*row["amplitude_v"])*polarity*1000
                ax.axhline(level,color="#777777",ls=":",lw=.7)
            reason = str(row.get("qc_reasons") or "no candidate QC flags")
            ax.set(title=f"{int(row['trigger_code'])}: {row['measurement_status']} | sample {int(row['digital_sample_raw'])}\n{reason[:85]}",
                   xlabel="Time relative to original STI101 (ms)",ylabel="Raw diode (mV)")
            ax.title.set_fontsize(8); ax.grid(alpha=.15)
            name = f"window_{len(index):03d}"
            arrays[name] = np.vstack([t,values])
            index.append(dict(array=name,event_key=row["event_key"],source_fif=source,
                              selection="cached_representative_or_pilot",baseline_method=row.get("baseline_method","unavailable")))
        for ax in list(axes.flat)[len(group):]: ax.axis("off")
        fig.suptitle(f"{source}: native waveform controls, tails and QC cases",fontsize=13)
        fig.tight_layout(rect=(0,0,1,.96))
        name = f"waveforms_recording_{page:02d}.png"
        files[name] = png(fig); pages.append((name,source,len(group)))
    if not pages:
        files["diagnostics.png"] = empty_figure("No raw waveform windows were retained for this evidence scope.","Waveform evidence unavailable")
    else: files["diagnostics.png"] = files[pages[0][0]]
    cache = io.BytesIO(); np.savez_compressed(cache,**arrays)
    files.update({"waveform_windows.npz":cache.getvalue(),"diagnostic_index.csv":csv_text(index,("array","event_key"))})
    return files,pages,index


def context_statistics(rows):
    """Stratify condition onsets by saved context; no PTB time becomes an edge."""
    groups = defaultdict(list)
    for row in rows:
        if int(row["trigger_code"]) not in (71,81,91) or row["measurement_status"]!="measured":
            continue
        contexts = {"after_warmup": "yes" if csv_bool(row.get("post_warmup_context",False)) else "no"}
        retry = number(row.get("saved_trial_onset_guard_retry_count"))
        contexts["guard_retry"] = "unknown" if retry is None else "yes" if retry>0 else "no"
        boundary = row.get("saved_trial_onset_boundary_trial")
        contexts["boundary_trial"] = "unknown" if boundary in (None,"") else "yes" if csv_bool(boundary) else "no"
        for factor,level in contexts.items():
            groups[(row["source_fif"],int(row["trigger_code"]),factor,level)].append(row["lag_ms"])
    return [dict(source_fif=source,trigger_code=code,context_factor=factor,context_level=level,
                 **stats(values)) for (source,code,factor,level),values in sorted(groups.items())]


def context_figure(rows):
    """Plot saved condition context separately by source and categorical code."""
    summaries = context_statistics(rows)
    if not summaries:
        return empty_figure("No measured condition-onset context comparisons are available.","Context comparisons"),summaries
    sources = sorted({r["source_fif"] for r in summaries})
    fig, axes = pyplot().subplots(len(sources),3,figsize=(15,3*len(sources)),squeeze=False,sharey=True)
    for source,axes_row in zip(sources,axes):
        for factor,axis in zip(("after_warmup","guard_retry","boundary_trial"),axes_row):
            for offset,code in zip((-.18,0,.18),(71,81,91)):
                group = [r for r in summaries if r["source_fif"]==source and r["context_factor"]==factor and r["trigger_code"]==code]
                for row in group:
                    x = ("no","yes","unknown").index(row["context_level"])+offset
                    axis.errorbar(x,row["mean"],yerr=row["sd"],fmt="o",capsize=3,color=TRIGGER_COLORS[code])
                    axis.annotate(f"n={row['n']}",(x,row["mean"]+(row["sd"] or 0)),xytext=(0,4),textcoords="offset points",ha="center",fontsize=7)
                axis.plot([],[],"o",color=TRIGGER_COLORS[code],label=str(code))
            axis.set(xticks=range(3),xticklabels=["No","Yes","Unknown"],xlabel=factor.replace("_"," "),
                     ylabel="Measured lag mean +/- SD (ms)",title=source)
            axis.grid(axis="y",alpha=.18); axis.legend(title="Trigger",fontsize=7,ncol=3); axis.margins(y=.3)
    fig.suptitle("Saved context and measured condition-onset lag: association, not a causal correction",fontsize=13)
    fig.tight_layout(rect=(0,0,1,.97))
    return png(fig),summaries


def interval_figure(intervals):
    """Compare only available measured same-trial endpoints to motion frames."""
    paired = [r for r in intervals if r.get("expected_frame_interval_ms") is not None
              and r.get("observed_diode_interval_ms") is not None]
    if not paired:
        return empty_figure("No pair of measured endpoints with comparable motion-frame indices is available in this scope.","Within-trial interval check")
    fig, axes = pyplot().subplots(1,2,figsize=(13,4))
    sources = sorted({str(r["left_event_key"]).split(":")[1] for r in paired})
    palette = ["#2563a1","#c46a1b","#627b32","#8653a3","#267c7b","#9b7333"]
    for i,source in enumerate(sources):
        group = [r for r in paired if str(r["left_event_key"]).split(":")[1]==source]
        axes[0].scatter([r["expected_frame_interval_ms"] for r in group],[r["observed_diode_interval_ms"] for r in group],s=12,label=source,color=palette[i%len(palette)])
        error = [r["observed_diode_interval_ms"]-r["expected_frame_interval_ms"] for r in group]
        s = stats(error)
        axes[1].bar(i,s["mean"],color=palette[i%len(palette)])
        if s["sd"] is not None: axes[1].errorbar(i,s["mean"],yerr=s["sd"],fmt="none",ecolor="#222222",capsize=4)
        axes[1].annotate(f"n={s['n']}",(i,s["mean"]+(s["sd"] or 0)),xytext=(0,4),textcoords="offset points",ha="center",fontsize=8)
    lower = min(r["expected_frame_interval_ms"] for r in paired)
    upper = max(r["expected_frame_interval_ms"] for r in paired)
    axes[0].plot([lower,upper],[lower,upper],color="#555555",ls="--",lw=.8)
    axes[0].set(xlabel="Saved motion-frame interval / refresh rate (ms)",ylabel="Measured diode-to-diode interval (ms)")
    axes[0].legend(fontsize=7)
    axes[1].set(xticks=range(len(sources)),xticklabels=[f"Recording {i+1}" for i in range(len(sources))],
                ylabel="Observed minus frame-based interval: mean +/- SD (ms)",xlabel="Recordings (same order as legend)")
    axes[1].axhline(0,color="#555555",lw=.8)
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle("Physical sequence sanity check: no missing endpoint is reconstructed")
    fig.tight_layout(rect=(0,0,1,.93))
    return png(fig)


# -----------------------------------------------------------------------------
# Additional interval views preserve the original bar figure and exact pairs
# -----------------------------------------------------------------------------


def interval_observations(intervals, rows, contexts):
    """Join existing comparable pairs to recorded events, verifying arithmetic.

    Inputs are normalized existing interval/event rows and saved block context.
    Output has one row per measured comparable pair, with exact keys, source,
    trigger factors, native recording time and saved trial metadata. Missing
    expected intervals remain absent. Contradictory joins or sample/frame
    arithmetic fail instead of being silently repaired or remeasured.
    """
    lookup = {r["event_key"]:r for r in rows}
    if len(lookup)!=len(rows):raise ValueError("Duplicate event keys in interval evidence.")
    identity = ("source_fif","behavioral_block","behavioral_run","behavioral_trial",
                "is_replay_trial","reported_source_csv")
    result = []
    for pair in intervals:
        expected,observed = number(pair.get("expected_frame_interval_ms")),number(pair.get("observed_diode_interval_ms"))
        if expected is None or observed is None:continue
        try:
            left,right = lookup[pair["left_event_key"]],lookup[pair["right_event_key"]]
        except KeyError as error:
            raise ValueError("Interval endpoint is absent from event evidence.") from error
        left_context = tuple(str(left.get(field,"")) for field in identity)
        right_context = tuple(str(right.get(field,"")) for field in identity)
        if left_context!=right_context or left["measurement_status"]!="measured" or right["measurement_status"]!="measured":
            raise ValueError("Interval endpoints must be measured in one exact source/trial context.")
        sf = number(left.get("sfreq_hz"))
        if sf is None or sf<=0 or number(right.get("sfreq_hz"))!=sf:
            raise ValueError("Interval endpoints have incompatible native sampling rates.")
        measured = (right["diode_edge_sample_raw"]-left["diode_edge_sample_raw"])*1000/sf
        residual = observed-expected
        if not np.isclose(measured,observed,atol=1e-8,rtol=0) or not np.isclose(residual,pair.get("interval_error_ms"),atol=1e-8,rtol=0):
            raise ValueError("Interval sample/residual arithmetic contradicts saved evidence.")
        # The saved nominal motion-frame clock is the existing interval contract;
        # measuredRefreshHz and PTB timestamps do not replace its expectation.
        block = str(left["behavioral_block"])
        refresh = number(contexts.get(block,{}).get("refrate"))
        left_frame,right_frame = number(left.get("behavioral_frame")),number(right.get("behavioral_frame"))
        if refresh is None or refresh<=0 or left_frame is None or right_frame is None:
            raise ValueError("Comparable interval lacks saved frame/refresh context.")
        frame_expected = (right_frame-left_frame)*1000/refresh
        if not np.isclose(frame_expected,expected,atol=1e-8,rtol=0):
            raise ValueError("Interval expectation contradicts saved frame/refresh context.")
        a,b = int(left["trigger_code"]),int(right["trigger_code"])
        result.append(dict(pair,source_fif=left["source_fif"],runtime=left.get("runtime","unknown"),
            left_code=a,right_code=b,endpoint_pair=f"{a}->{b}",
            left_meaning=TRIGGER_NAMES.get(a,"Other marker"),right_meaning=TRIGGER_NAMES.get(b,"Other marker"),
            expected_ms=expected,observed_ms=observed,residual_ms=residual,
            recording_time_sec=(left["digital_sample_raw"]-left["first_samp"])/sf,
            trial_context="|".join(left_context),behavioral_block=left["behavioral_block"],
            behavioral_run=left.get("behavioral_run",""),behavioral_trial=left.get("behavioral_trial",""),
            is_replay_trial=left.get("is_replay_trial",""),
            left_lag_ms=left.get("lag_ms"),right_lag_ms=right.get("lag_ms"),
            saved_frame_missed_count=left.get("saved_frame_missed_count","")))
    keys = [(r["left_event_key"],r["right_event_key"]) for r in result]
    if len(set(keys))!=len(keys):raise ValueError("Duplicate comparable interval pairs.")
    return result


def interval_jitter(row, width=.28):
    """Separate overlapping category positions deterministically; never move ms."""
    key = str(row["left_event_key"])+"|"+str(row["right_event_key"])
    unit = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8],"big")/2**64
    return (unit-.5)*width


def interval_statistics(pairs, by_transition=False):
    """Describe source/pair groups, retaining correlated-pair and trial counts."""
    groups = defaultdict(list)
    for row in pairs:
        groups[(row["source_fif"],row["endpoint_pair"] if by_transition else "all")].append(row)
    return [dict(source_fif=source,endpoint_pair=transition,
                 n_trial_contexts=len({r["trial_context"] for r in group}),
                 **stats(r["residual_ms"] for r in group))
            for (source,transition),group in sorted(groups.items())]


def interval_violin_figure(intervals):
    """Add a full-range violin/point view beside the retained original bar view.

    Source groups, their order and left-panel coordinates match interval_figure.
    Every finite comparable pair is a dot, even for sparse/constant groups.
    Diamonds summarize the median; offset squares show mean +/- sample SD of
    those same complete observations. No central zoom or exclusion is applied.
    """
    paired = [r for r in intervals if r.get("expected_frame_interval_ms") is not None
              and r.get("observed_diode_interval_ms") is not None]
    if not paired:
        return empty_figure("No measured endpoint pair with comparable motion frames is available.","Within-trial interval distributions")
    fig,axes = pyplot().subplots(1,2,figsize=(14,4.8))
    sources = sorted({str(r["left_event_key"]).split(":")[1] for r in paired})
    for i,source in enumerate(sources):
        group = [r for r in paired if str(r["left_event_key"]).split(":")[1]==source]
        color = RECORDING_COLORS[i%len(RECORDING_COLORS)]
        axes[0].scatter([r["expected_frame_interval_ms"] for r in group],
                        [r["observed_diode_interval_ms"] for r in group],s=12,label=source,color=color)
        values = np.asarray([r["observed_diode_interval_ms"]-r["expected_frame_interval_ms"] for r in group])
        if len(values)>=3 and np.ptp(values)>1e-9:
            violin = axes[1].violinplot([values],positions=[i],widths=.65,showextrema=False)
            for body in violin["bodies"]:body.set_facecolor(color);body.set_alpha(.22)
        axes[1].scatter([i+interval_jitter(r) for r in group],values,s=13,alpha=.5,color=color,
                        edgecolors="none",zorder=3,label="One measured interval" if i==0 else None)
        summary = stats(values)
        axes[1].scatter([i],[summary["median"]],marker="D",s=35,facecolors="white",edgecolors="#222222",zorder=5,
                        label="Median" if i==0 else None)
        axes[1].errorbar(i+.30,summary["mean"],yerr=summary["sd"],fmt="s",markersize=4,
                         color="#222222",capsize=4,zorder=5,label="Mean +/- sample SD" if i==0 else None)
        axes[1].annotate(f"n={summary['n']} pairs",(i,1),xycoords=axes[1].get_xaxis_transform(),
                         xytext=(0,5),textcoords="offset points",ha="center",fontsize=8)
    low = min(r["expected_frame_interval_ms"] for r in paired)
    high = max(r["expected_frame_interval_ms"] for r in paired)
    axes[0].plot([low,high],[low,high],color="#555555",ls="--",lw=.8)
    axes[0].set(xlabel="Saved motion-frame interval / refresh rate (ms)",ylabel="Measured diode-to-diode interval (ms)")
    axes[0].legend(fontsize=7)
    axes[1].set(xticks=range(len(sources)),xticklabels=[f"Recording {i+1}" for i in range(len(sources))],
                xlim=(-.55,len(sources)-.45),ylabel="Observed - expected interval (ms)",
                xlabel="Recordings (same groups and order as left legend)")
    axes[1].axhline(0,color="#555555",lw=.8)
    axes[1].legend(fontsize=7,loc="lower center",bbox_to_anchor=(.5,1.08),ncol=3)
    for axis in axes:axis.grid(alpha=.15)
    fig.suptitle("Physical within-trial interval consistency: distributions and every measured pair\nRight: violin + individual intervals; sparse/constant groups retain their points",fontsize=12)
    fig.text(.02,.01,"Dots are endpoint-pair intervals; several pairs may share a trial. Horizontal jitter only separates points. No intervals are reconstructed or omitted.",fontsize=8)
    fig.subplots_adjust(left=.06,right=.99,bottom=.17,top=.80,wspace=.19)
    return png(fig)


def interval_residual_figures(pairs):
    """Draw duration, recording-time and categorical-transition residual panels.

    Exact source rows use shared full-range residual limits. Trigger-pair colours
    are stable across figures; numerical codes never become a continuous axis.
    Scatter panels retain original coordinates, so coincidence can hide dots.
    The transition view adds category-only jitter, medians and mean/sample SD.
    """
    names = ("interval_residual_vs_duration.png","interval_residual_vs_time.png","interval_residual_by_pair.png")
    if not pairs:
        return {name:empty_figure("No measured endpoint pair with comparable motion frames is available.","Within-trial interval residuals") for name in names}
    sources = sorted({r["source_fif"] for r in pairs})
    observed = {(r["left_code"],r["right_code"]) for r in pairs}
    order = [p for p in INTERVAL_PAIR_ORDER if p in observed]+sorted(observed-set(INTERVAL_PAIR_ORDER))
    colours = {p:INTERVAL_PAIR_COLORS[INTERVAL_PAIR_ORDER.index(p)] if p in INTERVAL_PAIR_ORDER
               else RECORDING_COLORS[i%len(RECORDING_COLORS)] for i,p in enumerate(order)}
    result = {}
    plt = pyplot()
    for name,field,xlabel,title in (
        (names[0],"expected_ms","Saved motion-frame expected interval (ms)","Does interval discrepancy depend on expected duration?"),
        (names[1],"recording_time_sec","Left endpoint time since this recording's first sample (s)","When within each recording do discrepant intervals occur?")):
        fig,axes = plt.subplots(len(sources),1,figsize=(14,2.9*len(sources)),squeeze=False,sharey=True)
        for source,axis in zip(sources,axes.flat):
            group = [r for r in pairs if r["source_fif"]==source]
            for a,b in order:
                subset = [r for r in group if (r["left_code"],r["right_code"])==(a,b)]
                axis.scatter([r[field] for r in subset],[r["residual_ms"] for r in subset],s=12,alpha=.55,
                             color=colours[(a,b)],label=f"{a}->{b}")
            axis.set(title=f"{source} | n={len(group)} pairs",xlabel=xlabel,ylabel="Observed - expected (ms)")
            axis.axhline(0,color="#333333",lw=.8);axis.grid(alpha=.15)
        handles,labels = axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,title="Endpoint trigger pair (categorical; same colours in every recording)",
                   loc="lower center",ncol=min(7,len(order)),fontsize=8)
        fig.suptitle(title+"\nAll finite residuals retained; shared y scale across recordings",fontsize=13)
        fig.tight_layout(rect=(0,.08,1,.96))
        result[name] = png(fig)
    fig,axes = plt.subplots(len(sources),1,figsize=(15,3.2*len(sources)),squeeze=False,sharey=True)
    for source,axis in zip(sources,axes.flat):
        for i,pair in enumerate(order):
            group = [r for r in pairs if r["source_fif"]==source and (r["left_code"],r["right_code"])==pair]
            values = np.asarray([r["residual_ms"] for r in group])
            if not len(values):continue
            if len(values)>=3 and np.ptp(values)>1e-9:
                violin = axis.violinplot([values],positions=[i],widths=.65,showextrema=False)
                for body in violin["bodies"]:body.set_facecolor(colours[pair]);body.set_alpha(.18)
            axis.scatter([i+interval_jitter(r,.30) for r in group],values,s=9,alpha=.5,color=colours[pair],zorder=3)
            summary = stats(values)
            axis.errorbar(i+.31,summary["mean"],yerr=summary["sd"],fmt="s",color="#222222",markersize=3,capsize=2,zorder=5)
            axis.scatter([i],[summary["median"]],marker="D",facecolors="none",edgecolors="#222222",s=24,zorder=5)
        axis.set(title=source,xticks=range(len(order)),xticklabels=[f"{a}->{b}" for a,b in order],
                 xlabel="Endpoint trigger pair (categorical; colours match duration/time panels)",
                 ylabel="Observed - expected (ms)",xlim=(-.6,len(order)-.4))
        axis.tick_params(axis="x",rotation=55);axis.axhline(0,color="#333333",lw=.8);axis.grid(axis="y",alpha=.15)
    fig.suptitle("Which endpoint transitions carry the discrepancy?\nViolin + every pair; open diamond = median; offset square/whisker = mean +/- sample SD",fontsize=13)
    fig.tight_layout(rect=(0,0,1,.96))
    result[names[2]] = png(fig)
    return result


def donor_figures(rows, estimates, subject):
    """Show only decided own-participant donors and existing estimator results."""
    if not estimates:
        return {},[]
    estimates={name:value for name,value in estimates.items() if value.get("estimator")!="cohort_event_median"}
    if not estimates:return {},[]
    plt = pyplot()
    names = sorted(estimates)
    fig,axes = plt.subplots(len(names),1,figsize=(13,3*len(names)),squeeze=False)
    fig_h,axes_h = plt.subplots(len(names),1,figsize=(13,3*len(names)),squeeze=False)
    output = []
    for name,axis,holdout in zip(names,axes.flat,axes_h.flat):
        estimate = estimates[name]
        code = int(name.split("code")[-1])
        target = {"condition_code":code,"behavioral_block":1,"acquisition_run":1}
        donors = [r for r in rows if donor_allowed(subject,target,r)]
        selected_keys = estimate.get("donor_event_keys")
        if selected_keys is not None:
            donors = [r for r in donors if r["event_key"] in selected_keys]
        sources = sorted({r["source_fif"] for r in donors})
        palette = ["#2563a1","#c46a1b","#627b32","#8653a3","#267c7b","#9b7333"]
        for i,source in enumerate(sources):
            group = [r for r in donors if r["source_fif"]==source]
            axis.scatter([r["lag_ms"] for r in group],np.full(len(group),i),s=16,label=source,color=palette[i%len(palette)])
        if estimate["status"]=="estimated":
            center = estimate["lag_estimate_ms"]
            ci = estimate["location_ci95_ms"]
            spread = estimate["donor_observed_p05_p95_ms"]
            y = len(sources)
            # Lines retain asymmetric bounds even if the estimate lies outside them.
            axis.plot(spread,[y+.22]*2,color="#999999",lw=5,label="Observed donor 5-95%")
            axis.plot(ci,[y]*2,color="#222222",lw=2,label="Location 95% bootstrap CI")
            axis.scatter([center],[y],marker="D",s=50,color="#222222",label="Proposed donor location")
            axis.set_yticks(range(len(sources)+1),[f"Donor recording {i+1}" for i in range(len(sources))]+["Estimate"])
        else:
            axis.set_yticks(range(len(sources)),[f"Donor recording {i+1}" for i in range(len(sources))])
            axis.text(.01,.95,"Insufficient donors: no numerical correction estimated",transform=axis.transAxes,va="top",fontsize=10)
        axis.set(title=f"{name}: {estimate['status']} (n={estimate.get('n_donors',0)})",xlabel="Measured donor lag / proposed location (ms)")
        if sources:axis.legend(fontsize=7,ncol=2,loc="lower right")
        axis.grid(axis="x",alpha=.18); axis.margins(y=.5)
        had_errors = False
        for i,(key,label) in enumerate((("held_out_block_errors_ms","Held-out recording"),("held_out_behavioral_run_errors_ms","Held-out behavioral run"))):
            values = estimate.get(key,[])
            if not values:continue
            had_errors = True
            s = stats(values)
            holdout.bar(i,s["mean"],color="#2563a1",alpha=.6)
            if s["sd"] is not None:holdout.errorbar(i,s["mean"],yerr=s["sd"],fmt="none",ecolor="#222222",capsize=4)
            holdout.scatter(np.full(len(values),i),values,color="#222222",s=18,zorder=4)
            holdout.annotate(f"n={len(values)}",(i,s["mean"]+(s["sd"] or 0)),xytext=(0,5),textcoords="offset points",ha="center",fontsize=8)
        holdout.set(title=name,xticks=[0,1],xticklabels=["Held-out recording","Held-out behavioral run"],
                    ylabel="Held-out median minus other donor location (ms)")
        holdout.axhline(0,color="#555555",lw=.8)
        if not had_errors:holdout.text(.03,.55,"No held-out error comparison available in this scope",transform=holdout.transAxes)
        output.append(dict(estimate_id=name,trigger_code=code,status=estimate["status"],n_donors=estimate.get("n_donors",0),
                           lag_estimate_ms=estimate.get("lag_estimate_ms"),donor_sd_ms=estimate.get("donor_sd_ms"),
                           event_weighted_median_ms=estimate.get("event_weighted_median_ms"),
                           location_ci95_ms=json.dumps(estimate["location_ci95_ms"]) if estimate.get("location_ci95_ms") else None,
                           donor_observed_p05_p95_ms=json.dumps(estimate["donor_observed_p05_p95_ms"]) if estimate.get("donor_observed_p05_p95_ms") else None,
                           transfer_uncertainty=estimate.get("transfer_uncertainty","unknown"),
                           target_context_limitation=estimate.get("target_context_limitation","")))
    fig.suptitle("Restricted donor evidence: location confidence, observed dispersion and transfer differ",fontsize=13)
    fig_h.suptitle("Donor validation: individual held-out errors, mean +/- SD; transfer remains unknown",fontsize=13)
    for f in (fig,fig_h):f.tight_layout(rect=(0,0,1,.95))
    return {"donor_estimates.png":png(fig),"donor_validation.png":png(fig_h)},output


def coverage_figure(coverage):
    """Compare eligible trials with unchanged denominators, retaining unassigned."""
    rows = [r for r in coverage if r["eligible"]]
    sources = sorted({r["source_fif"] for r in rows})
    if not sources:return empty_figure("No eligible Stage 2 trials in this report.","Epoch-window coverage")
    fig,ax = pyplot().subplots(figsize=(13,4.5))
    for offset,kind in ((-.2,"original"),(.2,"proposed")):
        bottoms = np.zeros(len(sources))
        for state,color in (("inside","#2563a1"),("outside","#c46a1b"),("unassigned","#a1a1aa")):
            counts = []
            for source in sources:
                values = [r[f"{kind}_inside"] for r in rows if r["source_fif"]==source]
                counts.append(sum(v is True if state=="inside" else v is False if state=="outside" else v is None for v in values))
            ax.bar(np.arange(len(sources))+offset,counts,bottom=bottoms,width=.36,color=color,
                   label=state if kind=="original" else None)
            bottoms += counts
        for i,total in enumerate(bottoms):ax.text(i+offset,total,kind,ha="center",va="bottom",rotation=60,fontsize=8)
    ax.set(xticks=range(len(sources)),xticklabels=[f"Recording {i+1}" for i in range(len(sources))],
           ylabel="Eligible trials (same denominator in each pair)",xlabel="Recording order matches the coverage table")
    ax.legend(title="Epoch-window status",ncol=3);ax.margins(y=.2)
    fig.suptitle("Does proposed timing change recording coverage? Original versus proposed anchors")
    fig.tight_layout(rect=(0,0,1,.93))
    return png(fig)


def pooled_lag_figure(donors):
    """Compare categorical trigger centers with eventwise observed dispersion."""
    if not donors:return empty_figure("No pooled measured donor events.","Cohort donor lag by trigger")
    plt=pyplot();fig,ax=plt.subplots(figsize=(9,5))
    codes=[71,81,91];means=[];sds=[];medians=[];ns=[]
    for code in codes:
        vals=np.asarray([float(r["lag_ms"]) for r in donors if int(r["trigger_code"])==code],float)
        means.append(float(vals.mean()) if len(vals) else np.nan)
        sds.append(float(vals.std(ddof=1)) if len(vals)>1 else 0.)
        medians.append(float(np.median(vals)) if len(vals) else np.nan);ns.append(len(vals))
    x=np.arange(3);bars=ax.bar(x,means,yerr=sds,capsize=5,color=[TRIGGER_COLORS[c] for c in codes],alpha=.72,label="Mean +/- sample SD")
    ax.scatter(x,medians,marker="D",s=58,facecolor="white",edgecolor="#222",zorder=4,label="Median used for target estimate")
    for i,(bar,n) in enumerate(zip(bars,ns)):ax.text(i,bar.get_height()+sds[i]+.03,f"n={n}",ha="center",va="bottom",fontsize=9)
    ax.set(xticks=x,xticklabels=[f"{c}\n{TRIGGER_NAMES[c]}" for c in codes],ylabel="Signed diode-edge lag from STI101 (ms)",
           xlabel="Trigger code (categorical factor)")
    ax.set_title("Pooled measured onset lag: center and eventwise spread");ax.grid(axis="y",alpha=.2);ax.legend()
    fig.tight_layout();return png(fig)


# -----------------------------------------------------------------------------
# Explained figures, exact evidence tables and participant-specific decisions
# -----------------------------------------------------------------------------


def table_html(rows, columns=None, limit=80):
    """Keep exact evidence readable; complete rows remain in embedded CSVs."""
    if not rows:return "<p>No observations available for this table.</p>"
    columns = columns or list(rows[0])
    escape = lambda value: html.escape(str(value))
    out = ["<div class='table-scroll'><table><thead><tr>"]
    out.extend(f"<th>{escape(k.replace('_',' '))}</th>" for k in columns)
    out.append("</tr></thead><tbody>")
    for row in rows[:limit]:
        out.append("<tr>")
        for key in columns:
            value = row.get(key)
            text = "unavailable" if value is None or value=="" else f"{value:.4g}" if isinstance(value,float) else str(value)
            out.append(f"<td>{escape(text)}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    if len(rows)>limit:out.append(f"<p>Showing {limit} of {len(rows)} rows. The embedded CSV download contains every row.</p>")
    return "".join(out)


def participant_guidance(subject, manifest, estimates):
    """Explain decided restrictions and the decisions that remain genuinely open."""
    general = [
        "For measured trials, the proposal uses each trial's own usable 50% edge. The mean or median in a plot summarizes the evidence; it is not automatically applied as one correction to all trials.",
        "Compare trigger centers AND their SD/5-95% spread within recordings, then compare the same triggers across recordings. A similar center with large dispersion does not establish precise trial timing. A shift between groups needs waveform/context review before pooling.",
        "There is no accepted numerical consistency cutoff. Judge acceptable spread and differences against the intended temporal analysis, native sampling interval and detector sensitivity. Record the rationale and unresolved uncertainty with the exact timing bundle.",
        "A hypothetical +19 ms lag at 1,000 Hz moves the proposed anchor 19 samples later. This example illustrates sign/units; use the actual per-event proposals or restricted donor estimate shown here.",
        "Acceptance requires independent reconciliation of method, measurements, QC, estimates/fallbacks and uncertainty. This report does not author acceptance or final participant inclusion.",
    ]
    specific = {
        "11": ["P11 blocks 1-2 have confirmed unusable diode. Only P11's own usable later blocks, starting at block 3, can donate; early-block assignments must remain estimated.",
               "Inspect condition-specific donor centers/spread, block differences, estimator weighting and held-out errors. The patch moved from [193,0] to [245,208]; its spatial/scanout transfer error remains unknown, even if a bootstrap interval is narrow."],
        "24": ["P24 acquisitions 1-4 have no diode and irreversible 0.1 Hz online high-pass. Only acquisition-5 measured lag can donate; earlier timing remains estimated.",
               "One donor acquisition permits within-acquisition validation, not validation of transfer to acquisitions 1-4. Keep the acquisition-setting and unknown transfer limitations alongside any proposed correction."],
        "09": ["P09 acquisition 2 has a session-note detachment lead. Compare its temporal detection coverage, amplitude/noise and waveforms with neighboring acquisitions; sampled measurable edges alone do not establish whole-recording usability."],
        "08": ["P08 block 1 used v24; later blocks used v25. Inspect these saved-runtime groups separately before deciding whether their measured timing supports the same interpretation."],
        "17": ["P17 has confirmed no diode and irreversible 0.1 Hz online high-pass throughout. The proposed digital fallback keeps the anchor unchanged; physical visual latency and its uncertainty remain unknown."],
    }.get(subject, [])
    if subject in ("04","05","06","07"):
        specific = ["This participant has confirmed no diode and used v20 with a pre-draw trigger and unscheduled first flip. The proposed digital fallback leaves the anchor unchanged; the historical v20 counterfactual is not a participant calibration."]
    pooled=manifest.get("onset_estimation_method",{})
    if subject in ("04","05","06","07","17") and pooled.get("method_id")=="ordinary_median_of_unique_qc_usable_measured_cohort_events_v1":
        specific=["The onset estimate uses the cohort median of measured diode lags and is labelled estimated. The cohort spread describes donor events; transfer uncertainty into this participant's runtime/acquisition remains unknown."]
    insufficient = [name for name,value in estimates.items() if value.get("status")=="insufficient_donors"]
    if insufficient:specific.append("Insufficient donor evidence in this scope for "+", ".join(insufficient)+". No numerical estimate is available for those conditions; do not infer one from a few illustrative points.")
    if manifest.get("mode")=="pilot":specific.insert(0,"PILOT ONLY: the sampled windows and incomplete timing map cannot establish full-recording coverage, a full donor calibration or a production epoching handoff.")
    return specific+general


def build_report_artifacts(measurements, maps, estimates, recordings, intervals, contexts,
                           manifest, waveforms, *, source_identity=None, proposals=None, evidence_bytes=None,
                           donor_measurements=None):
    """Return new report/plot/table bytes without performing I/O or measurements.

    All inputs are copied before CSV normalization. Waveforms are already-read
    ``(STI101/diode array, anchor index, native sfreq)`` tuples. Descriptive
    statistics and report findings are derived from these exact observations;
    no estimator, exclusion or timing-map assignment is changed here.
    """
    original_measurements, original_maps = measurements, maps
    original_recordings, original_intervals = recordings, intervals
    rows, maps = normalize(measurements), normalize(maps)
    recordings, intervals = normalize(recordings), normalize(intervals)
    subject = str(manifest["subject_bids_id"]).replace("sub-","")
    method = manifest.get("method",DEFAULT_METHOD)
    window = manifest.get("epoch_window_sec",[-.5,2.667])
    lag_rows = lag_statistics(rows)
    corrections = correction_statistics(maps)
    coverage, checks = coverage_and_checks(maps,window)
    checks["stored_anchor_code_mismatches"] = sum(
        number(r.get("anchor_recorded_STI101_value")) is not None
        and int(float(r["anchor_recorded_STI101_value"]))!=int(r["trigger_code"]) for r in rows)
    checks["measured_lag_sample_mismatches"] = sum(
        not np.isclose((r["diode_edge_sample_raw"]-r["digital_sample_raw"])*1000/r["sfreq_hz"],r["lag_ms"],atol=1e-8,rtol=0)
        for r in rows if r["measurement_status"]=="measured" and r.get("diode_edge_sample_raw") is not None and r.get("lag_ms") is not None)
    qc = [r for r in rows if r["measurement_status"]!="measured"]
    counts = dict(Counter(r["measurement_status"] for r in rows))
    eligible = [r for r in maps if analysis_eligible(r)]
    original_eligible = sum(csv_bool(r.get("included_by_trial_policy",False)) for r in maps)
    bases = dict(Counter(r["assignment_basis"] for r in eligible))
    coverage_by_source = []
    for source in sorted({r["source_fif"] for r in coverage}):
        group = [r for r in coverage if r["source_fif"]==source and r["eligible"]]
        coverage_by_source.append(dict(source_fif=source,n_eligible=len(group),
            original_inside=sum(r["original_inside"] is True for r in group),
            proposed_inside=sum(r["proposed_inside"] is True for r in group),
            proposed_outside=sum(r["proposed_inside"] is False for r in group),
            timing_unassigned=sum(r["proposed_inside"] is None for r in group),
            newly_outside=sum(r["original_inside"] is True and r["proposed_inside"] is False for r in group),
            uncertainty_inside=sum(csv_bool(r["uncertainty_inside"]) for r in group if r["uncertainty_inside"] not in (None,"")),
            uncertainty_unknown=sum(r["uncertainty_inside"] in (None,"") for r in group),
            trial_end_overlap=sum(csv_bool(r["overlaps_trial_end"]) for r in group if r["overlaps_trial_end"] not in (None,""))))
    signal_rows = []
    for source,runtime in strata(rows):
        group = [r for r in rows if (r["source_fif"],r.get("runtime","unknown"))==(source,runtime)]
        for code in sorted({int(r["trigger_code"]) for r in group}):
            subset = [r for r in group if int(r["trigger_code"])==code]
            for field,scale,unit in (("amplitude_v",1000,"mV"),("baseline_noise_v",1000,"mV"),("threshold_span_ms",1,"ms")):
                s = stats(r[field]*scale for r in subset if r.get(field) is not None)
                signal_rows.append(dict(source_fif=source,runtime=runtime,trigger_code=code,measure=field,unit=unit,**s))
    # Complete tail inventory retains extrema within strata, rather than using
    # a global central-lag cutoff to discard genuine early/late measurements.
    tails = []
    for summary in lag_rows:
        group = sorted([r for r in rows if r["source_fif"]==summary["source_fif"]
                        and r.get("runtime","unknown")==summary["runtime"]
                        and int(r["trigger_code"])==summary["trigger_code"]
                        and r["measurement_status"]=="measured" and r.get("lag_ms") is not None],key=lambda r:r["lag_ms"])
        seen = set()
        for r in group[:2]+group[-2:]:
            if r["event_key"] not in seen:tails.append(r);seen.add(r["event_key"])
    interval_details = interval_observations(intervals,rows,contexts)
    interval_source_stats = interval_statistics(interval_details)
    interval_pair_stats = interval_statistics(interval_details,by_transition=True)
    # This is an inspection inventory, not a new outlier or trial-rejection rule.
    interval_tails = sorted(interval_details,key=lambda r:abs(r["residual_ms"]),reverse=True)[:20]
    figures = {}
    figures.update(lag_figures(rows)); figures.update(time_figures(rows,method))
    figures["measurement_coverage.png"] = measurement_coverage_figure(rows)
    waveform_files,pages,waveform_index = waveform_figures(rows,waveforms)
    figures.update({k:v for k,v in waveform_files.items() if k.endswith(".png")})
    context_image,context_rows = context_figure(rows)
    figures["context_comparisons.png"] = context_image
    figures["physical_intervals.png"] = interval_figure(intervals)
    figures["physical_intervals_violin.png"] = interval_violin_figure(intervals)
    figures.update(interval_residual_figures(interval_details))
    figures["epoch_coverage.png"] = coverage_figure(coverage)
    donor_images,donor_rows = donor_figures(rows,estimates,subject)
    figures.update(donor_images)
    pooled_rows=donor_measurements or []
    pooled_image=pooled_lag_figure(pooled_rows) if pooled_rows else None
    if pooled_image is not None:
        figures["cohort_pooled_trigger_lag.png"]=pooled_image
        if "donor_estimates.png" not in figures:figures["donor_estimates.png"]=pooled_image
        if "donor_validation.png" not in figures:
            figures["donor_validation.png"]=empty_figure("No held-out target-transfer validation was performed.","Transfer error remains unknown")
    # Keep the legacy filename useful while the report uses descriptive names.
    figures["lag_sensitivity.png"] = figures["lag_by_trigger.png"]
    data = {
        "event_measurements.csv":csv_text(original_measurements,("event_key","measurement_status")),
        "trial_timing_map.csv":csv_text(original_maps),
        "recording_summary.csv":csv_text(original_recordings),
        "physical_intervals.csv":csv_text(original_intervals,("left_event_key","right_event_key","observed_diode_interval_ms","expected_frame_interval_ms","interval_error_ms")),
        "lag_statistics.csv":csv_text(lag_rows), "correction_summary.csv":csv_text(corrections),
        "signal_statistics.csv":csv_text(signal_rows), "context_statistics.csv":csv_text(context_rows),
        "interval_details.csv":csv_text(interval_details,("left_event_key","right_event_key","source_fif","endpoint_pair","expected_ms","observed_ms","residual_ms")),
        "interval_statistics.csv":csv_text(interval_source_stats,("source_fif","endpoint_pair","n_trial_contexts","n","mean","sd","median")),
        "interval_transition_statistics.csv":csv_text(interval_pair_stats,("source_fif","endpoint_pair","n_trial_contexts","n","mean","sd","median")),
        "interval_tail_pairs.csv":csv_text(interval_tails,("left_event_key","right_event_key","source_fif","endpoint_pair","expected_ms","observed_ms","residual_ms")),
        "event_qc_review.csv":csv_text(qc,("event_key","measurement_status","qc_reasons")),
        "tail_events.csv":csv_text(tails,("event_key","lag_ms")),
        "coverage_review.csv":csv_text(coverage), "coverage_summary.csv":csv_text(coverage_by_source),
        "donor_summary.csv":csv_text(donor_rows,("estimate_id","status","n_donors")),
        "estimates.json":json.dumps(estimates,indent=2,sort_keys=True,allow_nan=False),
        "contexts.json":json.dumps(contexts,indent=2,sort_keys=True,allow_nan=False),
        "timing_sanity_checks.json":json.dumps(checks,indent=2,sort_keys=True,allow_nan=False),
    }
    if proposals is not None:
        data["event_timing_proposals.csv"] = csv_text(proposals,("event_key","assignment_basis"))
    if pooled_rows:
        data["cohort_pooled_donor_measurements.csv"]=csv_text(pooled_rows)
        data["cohort_pooled_donor_provenance.json"]=json.dumps(manifest.get("cohort_pooled_donor_provenance",{}),indent=2,sort_keys=True)
    if evidence_bytes:
        allowed = {"event_measurements.csv","event_timing_proposals.csv","trial_timing_map.csv",
                   "recording_summary.csv","physical_intervals.csv","estimates.json","contexts.json","source_manifest.json",
                   "cohort_pooled_donor_measurements.csv","cohort_pooled_donor_provenance.json"}
        if not set(evidence_bytes).issubset(allowed):raise ValueError("Unexpected report evidence file.")
        data.update(evidence_bytes)
    data.update({k:v for k,v in waveform_files.items() if not k.endswith(".png")})
    guidance = participant_guidance(subject,manifest,estimates)
    spread_rows = [r for r in lag_rows if r["sd_lag_ms"] is not None]
    widest = max(spread_rows,key=lambda r:r["sd_lag_ms"]) if spread_rows else None
    spread_text = (f"Largest observed within-recording/code SD is {fmt(widest['sd_lag_ms'])} ms for {int(widest['trigger_code'])} in {widest['source_fif']} (n={widest['n_measured']})."
                   if widest else "No stratum has two usable measurements; sample SD is unavailable.")
    measured = [r for r in rows if r["measurement_status"]=="measured" and r.get("lag_ms") is not None]
    measured_summary = stats(r["lag_ms"] for r in measured)
    zero_noise = sum(r.get("baseline_noise_v")==0 for r in rows)
    comparable = [r for r in intervals if r.get("interval_error_ms") is not None]
    interval_summary = stats(r["interval_error_ms"] for r in comparable)
    interval_finding = (f"{len(interval_details)} comparable measured pairs from "
        f"{len({r['trial_context'] for r in interval_details})} distinct recorded trial contexts. "
        f"Residual mean {fmt(interval_summary['mean'])} ms, sample SD {fmt(interval_summary['sd'])} ms, "
        f"median {fmt(interval_summary['median'])} ms; range {fmt(interval_summary['minimum'])} to {fmt(interval_summary['maximum'])} ms.")
    interval_example = "No comparable measured interval is available; do not replace missing pairs with zero residuals."
    if interval_tails:
        worst = interval_tails[0]
        shared = [r for r in interval_details if r["left_event_key"]==worst["left_event_key"]]
        interval_finding += (f" Largest absolute residual: {fmt(worst['residual_ms'],3)} ms for {worst['endpoint_pair']} "
            f"in {worst['source_fif']}, at {fmt(worst['recording_time_sec'],3)} s; "
            f"block {worst['behavioral_block']}/run {worst['behavioral_run']}/trial {worst['behavioral_trial']}, "
            f"replay={worst['is_replay_trial']}. {len(shared)} pairs share this left endpoint; they are not independent trials.")
        interval_example = (f"Actual example: {worst['endpoint_pair']} in {worst['source_fif']} has a saved expectation "
            f"of {fmt(worst['expected_ms'],3)} ms and a measured interval of {fmt(worst['observed_ms'],3)} ms. "
            f"Subtracting gives {fmt(worst['residual_ms'],3)} ms. Locate both exact keys in interval_details.csv "
            "and compare other pairs from that trial before deciding what explains the difference.")
    interval_derivation = (
        "Observed interval = (right diode_edge_sample_raw - left diode_edge_sample_raw) * 1000 / native sfreq_hz. "
        "Expected interval = (right behavioral_frame - left behavioral_frame) * 1000 / saved Conf.refrate. "
        "Residual = observed - expected: a positive value means a longer measured interval. Both endpoints must be measured "
        "in the same exact source/block/run/trial/replay/source-CSV context; these joins and sample/frame arithmetic are verified. "
        "The expected clock is nominal saved motion refresh, not measuredRefreshHz or PTB time. Question-UI frames are not comparable. "
        "At 1000 Hz the endpoint grid is 1 ms; a 120-Hz expectation uses 8.333... ms frame increments. Fractional residuals and smoothing add no measurement precision. "
        "No absent endpoint, unrecorded trial or dot-offset marker is reconstructed. These measured events are not filtered by a later final trial-selection policy.")
    interval_marks = (
        "One coloured dot is one measured endpoint-pair interval; several pairs may share one trial or endpoint. "
        "A median is the middle sorted residual, averaging the middle two for even n. Mean = sum(residuals)/n; "
        "sample SD = sqrt(sum((residual-mean)^2)/(n-1)), unavailable at n=1. Whiskers span mean-SD to mean+SD: "
        "observed pair dispersion, not a confidence interval, independent-trial uncertainty or an acceptance cutoff. "
        "Violin width describes smoothed observed density; sparse/effectively constant groups retain dots without a density. "
        "Horizontal category jitter only separates points and never changes milliseconds. Every finite comparable pair is retained.")
    new_boundary = sum(r["newly_outside"] for r in coverage_by_source)
    cards = []
    def card(identifier,title,file,axes,derived,finding,decision,example):
        cards.append(dict(id=identifier,title=title,file=file,axes=axes,derived=derived,
                          finding=finding,decision=decision,example=example))
    card("violin","Are the trigger lags consistent, and what shift does the evidence support?","delay_violin_with_dots.png",
         "X: signed lag in ms, positive when the diode edge follows STI101. All panels share limits based on group means +/- two SD and every finite observed point, with a small margin; no tail is clipped and zero is not forced into a positive-only view. Y: trigger code and meaning as categories; n is usable/assessed. Each panel is one exact recording/runtime. Recordings without any finite lag are listed in a compact notice instead of empty panels.",
         "Each coloured dot is ONE optically marked event's measured lag: find its native 50% diode edge, subtract that event's recorded STI101 sample, and multiply by 1000/sfreq_hz. Colour identifies the trigger category. Within each recording/runtime/trigger group, sort the usable event lags: the open diamond is the middle value, or the average of the two middle values for even n. A diamond is a descriptive median, not another measurement; for n=1 it surrounds the single event's dot at the same lag. Lines give observed 5-95% percentiles. Only measurement_status=measured enters the violin/median/range; flagged finite candidates are black x marks. These are measured marked events; accepted trial eligibility is reported separately. Jitter only separates dots vertically and encodes neither timing nor trial order. A violin requires at least three nonidentical lags; smoothing and a median between samples add no measurement precision.",
         f"{len(measured)} usable measured edges across {len(strata(rows))} recording/runtime strata; observed lag range {fmt(measured_summary['minimum'])} to {fmt(measured_summary['maximum'])} ms. {spread_text}",
         "Use centers and dispersion to judge lag consistency across triggers within each recording, then compare the same trigger across recordings. Inspect tails and separated centers before considering any pooling. Use the correction table for the actual assigned shift; the plot alone does not select a correction.",
         "If hypothetical code-71 and code-81 distributions both center near 19 ms with small SD, that supports a common latency description. If one centers near 13 ms or has a wide tail, inspect its raw edge, context and threshold sensitivity before accepting or pooling it.")
    card("bars","How large is each trigger's measured lag and dispersion?","lag_by_trigger.png",
         "X: equally spaced trigger categories, not a numeric trigger scale. Y: mean measured lag in ms. Bar height is the mean; whiskers are +/- one sample SD (ddof=1); diamonds are medians. Labels give usable n. Separate panels retain recording/runtime identity.",
         "For each source_fif/runtime/trigger_code group, average measured lag_ms and compute SD using denominator n-1. An n=1 group has no SD whisker; no usable group has no bar. Warm-up 153 remains a diagnostic category. SD is not a confidence interval or a calibrated error bound.",
         spread_text,
         "Read the size of lag and its event-to-event spread together. Compare means with medians and the individual dots in the violin; a mean pulled by a real tail is not a reason to delete that event. Measured trials retain their own edges; donor corrections use their separately reviewed estimator.",
         "A hypothetical mean of 19 ms with SD 0.5 ms describes a tighter observed group than mean 19 ms with SD 4 ms. Neither result, by itself, establishes the unknown patch-to-dot transfer error.")
    if pooled_rows:
        ps=manifest.get("onset_estimation_method",{}).get("summary",{})
        card("cohort-pooled","What measured cohort lag supports the estimated onset?","cohort_pooled_trigger_lag.png",
             "X: trigger code as a categorical factor. Y: signed difference (ms) from original STI101 sample to measured diode edge; positive means the diode edge follows the trigger. Bars show eventwise mean +/- sample SD; white diamonds show ordinary medians used for target shifts. Labels give unique measured event counts.",
             "Keep each unique event whose status is measured, qc_reasons is empty and native-sample arithmetic reproduces lag_ms. Exclude only the five explicit P24 donor keys. Compute ordinary medians independently for 71, 81 and 91, with no bootstrap or block weighting. Exact event identities and source hashes are in the donor CSV and provenance JSON.",
             "; ".join(f"{code}: median {fmt(ps.get(str(code),{}).get('median_ms'))} ms (n={ps.get(str(code),{}).get('n','unavailable')})" for code in (71,81,91))+". Observed donor SD does not quantify transfer error into an unmeasured participant.",
             "Compare each trigger's median with its SD and observed range. Similar centers and modest spread support a common descriptive pattern; a trigger difference or broad spread calls for reviewing donor QC and event family. Keep target transfer uncertainty unknown.",
             "If all three medians are 19 ms with narrow spread, apply the native-rounded 19-ms shift as an estimate to the designated participants. This estimates trigger-to-patch lag; it does not prove exact dot luminance onset.")
    card("time","Is timing stable within each recording?","lag_over_time.png",
         "X: seconds since that exact FIF's first recorded sample. Y: signed candidate lag in ms. Trigger identity is categorical colour/legend; dots are measured and x marks other finite candidates. Panels never overlay clocks from different recordings.",
         "Recording time = (digital_sample_raw - first_samp) / sfreq_hz. Lag comes from the native 50% diode edge. Missing/unusable events have no invented y value; their counts appear in coverage/QC tables.",
         f"{len(measured)} usable points are plotted with {sum(r['measurement_status']!='measured' and r.get('lag_ms') is not None for r in rows)} separately identified finite non-usable candidates. A pilot samples a few times and cannot establish continuous stability.",
         "Look for drift, abrupt shifts, changing spread and isolated tails in the same trigger. Reconcile them with waveforms, acquisition boundaries and saved context before deciding whether a recording or segment supports measured timing or donor transfer.",
         "If a hypothetical recording changes from 18 to 23 ms halfway through, do not average it into one stable correction. Verify the shift and its context first; the report does not automatically fit a drift correction.")
    card("signal","Does the diode signal remain usable?","signal_usability.png",
         "Each row is one recording/runtime. X: seconds within that recording. Y columns: response amplitude in mV, robust baseline noise in mV, and 10-90% threshold span in ms. Every panel has limits computed from that row's own finite values; the amplitude panel also includes its candidate-minimum reference line. Colours/legend identify trigger categories. The reference is not an accepted QC decision.",
         "Amplitude and baseline_noise_v are detector outputs multiplied by 1,000. Robust noise is 1.4826 times the baseline median absolute deviation. Edge span is the native-grid time between 10% and 90% crossings. Per-recording/code mean, SD and counts are in signal_statistics.csv.",
         f"{sum(r.get('amplitude_v') is not None for r in rows)} assessed events have amplitude values. {zero_noise} have zero robust baseline-noise estimates, which may reflect repeated/quantized baseline samples and do not establish physically noiseless recording.",
         "Inspect attenuation/dropouts with detection coverage and actual traces. A channel name or a large amplitude does not prove an unambiguous edge. A small/zero robust-noise estimate is not sufficient evidence of high timing precision.",
         "For P09-like attenuation, measurable isolated pulses support those observations, while an interval of missing/ambiguous edges still requires temporal usability review. Do not classify a whole recording from a few clean pilot points.")
    card("measurement-coverage","How much assessed timing is measured, flagged or unavailable?","measurement_coverage.png",
         "X: equally spaced trigger categories. Y: assessed event counts, stacked by measurement status. Panels retain exact recording/runtime; colours/legend describe measured, flagged, unusable or not-recorded states. Counts do not have an SD error bar.",
         "Count every event_measurements.csv row by source_fif/runtime/trigger_code/measurement_status. The denominator is assessed marked events in that stratum, including failures. A pilot assesses selected windows only; it does not establish whole-recording detection coverage. Stage 2 trial eligibility is a separate denominator.",
         "; ".join(f"{name}: {count}" for name,count in sorted(counts.items())),
         "Read this beside lag distributions: a tight violin among a few successful detections can coexist with substantial missing or failed measurements. Review the failure reasons and temporal pattern before accepting recording usability or donor coverage.",
         "If hypothetical code 71 has 100 assessed events but only 20 measured edges, the 20-edge mean and SD do not establish reliable timing for the other 80. Keep their statuses explicit instead of imputing the mean automatically.")
    card("threshold","How sensitive is lag to the edge definition?","threshold_sensitivity.png",
         "Left X: response fractions 10/25/50/75/90% as ordered threshold categories; Y: corresponding mean lag +/- sample SD in ms, separately by trigger and recording. Right X: candidate 50% lag (ms); Y: individual 10-90% span (ms). The dashed span reference is a candidate flag threshold.",
         "Each threshold lag is edge_LEVEL_offset_samples * 1000 / sfreq_hz. Fraction-specific means/SD use measured events. The scatter retains all finite candidates with a span. Neither threshold span nor SD is a calibrated confidence interval for luminance onset.",
         f"{sum(r.get('threshold_span_ms') is not None for r in rows)} events have stored threshold spans. The proposal's sensitivity flag is {fmt(method['maximum_sensitivity_ms'])} ms; its scientific acceptance remains open.",
         "A shallow/ambiguous rise can make edge choice matter. Inspect traces and polarity when thresholds disagree or extra rises appear; do not hide a problematic threshold by accepting the mean lag alone.",
         "A hypothetical 10-90% span of 1 ms supports less edge-definition sensitivity than a span of 6 ms. It does not remove sampling-grid or patch-to-dot limitations.")
    card("context","Do warm-up, guard and boundary contexts differ?","context_comparisons.png",
         "Rows: exact recording. Columns: preceding warm-up, saved guard retry, saved boundary-trial flag. X: No/Yes/Unknown categories. Y: measured condition-onset lag mean +/- sample SD (ms); colour/legend is trigger 71/81/91. Labels give subgroup n.",
         "Only usable condition onsets enter this comparison. Canonical previous-153 context identifies post-warm-up onsets; saved MAT flags provide guard/boundary groups. PTB values provide context, never substitute optical observations. SD is unavailable for n=1 and a zero-length drawn whisker is not a known zero error.",
         f"{len(context_rows)} observed source/code/context subgroups are available. Missing groups have no fabricated comparison or numerical difference.",
         "Check whether apparent tails concentrate in a saved context, then inspect the matching waveform and trial key. An association is not a causal correction or permission to repair arbitrary baselines.",
         "If a low-lag first onset follows warm-up, verify that the chosen edge belongs to this trial rather than the preceding pulse before accepting its lag. Keep a genuine early edge when the waveform supports it.")
    card("interval","Are physical within-trial intervals consistent with the saved sequence?","physical_intervals.png",
         "Left X: comparable saved motion-frame interval in ms; Y: measured diode-to-diode interval in ms, with the equality diagonal. Right X: recording categories; Y: observed minus frame-based interval, mean +/- SD (ms). Colours identify recordings.",
         interval_derivation+" The original mean/sample-SD bar figure is retained. Sample SD uses n-1; n counts pairs, not independent trials.",
         interval_finding,
         "Use mismatches to investigate edge pairing, event families and frame context. Equal common latency can cancel in an interval, so a good interval check does not determine the absolute lag correction.",
         "If two measured edges each occur 19 ms after their anchors, their difference can still match the saved 500-ms interval. That supports sequence consistency while leaving the absolute onset correction to the lag evidence.")
    card("interval-violin","How are the within-trial interval residuals distributed?","physical_intervals_violin.png",
         "Left X: saved expected duration (ms); Y: measured diode-to-diode duration (ms), with equality diagonal. Right X: the same recording factors/order as the bar figure; Y: residual (ms). Colours identify sources in the left legend. Dots are individual pairs, open diamonds medians, offset black squares/whiskers mean +/- sample SD; n counts pairs.",
         interval_derivation+" "+interval_marks,
         interval_finding,
         "Compare the individual observations, center and spread across recordings. Sparse tails can drive the mean/SD while most pairs agree closely. Inspect exact trial relationships instead of deleting tails as independent failures. Common latency can cancel between endpoints, so this supplements the separate trigger-to-diode lag review; diode-patch edges do not establish exact dot luminance onset.",
         interval_example)
    card("interval-duration","Does interval discrepancy depend on expected interval length?","interval_residual_vs_duration.png",
         "X: saved motion-frame expected interval length (ms). Y: measured minus expected interval (ms). Rows are exact source recordings with a shared full-range Y scale. Colour/legend identifies categorical left->right trigger pairs. Each dot is one existing measured pair; coincident coordinates can overlap. Zero means interval agreement.",
         interval_derivation+" The same points as the primary figure are replotted with residual on Y to expose small discrepancies. No fit, resampling, duration-based correction or point filtering is applied. Pair-specific mean, sample SD, medians and counts are in interval_transition_statistics.csv.",
         interval_finding,
         "Look for duration-related differences within the same transition/recording context. Duration also differs between event families, so an apparent slope alone does not establish clock drift. An isolated onset shared by several long residuals calls for checking its exact sequence and endpoints.",
         "Hypothetical example: residuals growing from 1 ms at 200 ms to 10 ms at 1600 ms within a comparable transition may motivate a frame/clock-scale investigation. Several displaced pairs sharing one onset instead call for a focused trial review. "+interval_example)
    card("interval-recording-time","When within each recording do discrepant intervals occur?","interval_residual_vs_time.png",
         "X: the left endpoint's digital time since this exact recording's first sample (seconds). Y: residual (ms), shared across recording rows. Colours/legend identify categorical endpoint trigger pairs. Several intervals from the same left onset have the same X coordinate; clocks from different sources are never overlaid.",
         interval_derivation+" Recording time = (left digital_sample_raw - first_samp) / sfreq_hz. This anchors context to recorded time, without replacing measured diode samples with behavioral/PTB timestamps. Exact-pair values and source/pair descriptive dispersion are in the interval CSV companions.",
         interval_finding,
         "Distinguish isolated trial context from drift, abrupt changes or boundary effects within the same pair factor. Count distinct trials before describing frequency: several dots at one time can share one onset. Saved missed-frame counters are contextual leads, not causal proof or measured corrections.",
         "Hypothetical example: a gradual increase in the same transition over a recording merits drift investigation; four high pairs at one onset merit inspection of that trial. Neither pattern automatically changes the timing map. "+interval_example)
    card("interval-transition","Which endpoint trigger pairs carry interval discrepancies?","interval_residual_by_pair.png",
         "X: categorical endpoint trigger pairs, equally spaced in fixed sequence order (left code -> right code). Rows: exact source recordings. Y: residual (ms), shared across rows. Colours match duration/time panels. Dots are individual pairs, open diamonds medians, offset squares/whiskers mean +/- sample SD. Numerical code distance has no physical meaning.",
         interval_derivation+" Exact endpoint event keys supply the two trigger factors. "+interval_marks,
         interval_finding,
         "Compare onset-to-midtrial pairs with adjacent midtrial intervals. If several onset-to-event pairs are unusual but adjacent midtrial pairs are ordinary, investigate the shared start/sequence relationship before assuming independent delays at every later transition. Read the trigger definitions and exact trial tables below; nominal pairing still needs scientific review.",
         "Hypothetical example: 81->111 and 81->115 tails in one trial with tight 111->112 and 114->115 intervals support focused review of the shared onset/context. They do not prove an incorrect edge or justify a correction. "+interval_example)
    if estimates:
        card("donors","What correction and uncertainty do the permitted donors support?","donor_estimates.png",
             "X: measured donor lag or existing estimated location in ms. Y: donor recording categories plus an estimate row when available; panels are condition-specific. Dots are donor observations, a diamond the location, a thin dark line its 95% bootstrap CI, and a separate thick light line the observed donor 5-95% range.",
             "Only existing P11-own-later-block or P24-acquisition-5 donor measurements are displayed. The existing estimate uses the median of donor-recording medians; event-weighted median, SD, hierarchical-bootstrap location CI and exact donor keys remain inspectable. Insufficient donors produce no numerical estimate.",
             "; ".join(f"{name}: {e['status']}, n={e.get('n_donors',0)}" for name,e in sorted(estimates.items())),
             "Review donor QC, between-recording consistency, weighting and held-out results before accepting the estimated early-block correction. CI precision, observed donor jitter and unknown cross-context transfer are different quantities. No CI includes unknown patch/settings transfer error.",
             "For P11, similar usable block-3/4/5 locations may support their restricted estimate for blocks 1-2 after review. A narrow location CI cannot demonstrate that the patch-position change has zero transfer error.")
        card("holdout","How well do donors predict another donor recording or run?","donor_validation.png",
             "X: held-out recording or held-out behavioral-run categories. Y: held-out median minus location from other donors in ms. Dots show individual errors; bar/whisker show mean +/- sample SD of available errors; n counts held-out units.",
             "These are the stored held_out_block_errors_ms and held_out_behavioral_run_errors_ms from the existing estimator. They are not new estimates and do not validate transfer into an unmeasured target recording.",
             f"{sum(len(e.get('held_out_block_errors_ms',[])) for e in estimates.values())} recording-level and {sum(len(e.get('held_out_behavioral_run_errors_ms',[])) for e in estimates.values())} behavioral-run error observations are available.",
             "Large or systematic held-out errors weaken pooling of donor recordings. With one donor acquisition, explicitly limit conclusions to within-acquisition validation. Empty comparisons remain unavailable.",
             "For P24, withholding behavioral runs inside acquisition 5 can check within-run stability, but cannot test transfer to acquisitions 1-4.")
    if pooled_rows:
        summary=manifest.get("onset_estimation_method",{}).get("summary",{})
        finding="; ".join(f"{code}: median {fmt(summary.get(str(code),{}).get('median_ms'))} ms, n={summary.get(str(code),{}).get('n','unavailable')}" for code in (71,81,91))
        card("cohort-pooled","What measured cohort lag supports the estimated onset?","cohort_pooled_trigger_lag.png",
             "X: trigger code as a categorical factor. Y: signed difference (ms) from original STI101 sample to the measured diode edge; positive means the diode edge follows the trigger. Bars are the eventwise mean +/- sample SD; white diamonds are ordinary medians, which supply the five target participants' proposed shifts. Labels show the number of distinct measured event keys.",
             "One row per unique event key from P04-P39 is retained when status is measured, qc_reasons is empty and diode_sample - STI101_sample divided by native sfreq matches lag_ms. The only removed keys are the five explicit P24 donor exclusions. The ordinary median is computed separately for 71, 81 and 91, without bootstrap or block weighting. Exact donor rows and source hashes are downloadable as cohort_pooled_donor_measurements.csv and cohort_pooled_donor_provenance.json.",
             finding+". The donor SD describes observed across-event spread; it does not quantify transfer error into a participant/runtime with no diode.",
             "Compare each trigger's center and spread. Similar medians with tight spread support a shared descriptive lag pattern; a trigger-specific shift or wide spread calls for reviewing event family, runtime and donor QC before accepting the target assignment. Keep target transfer uncertainty unknown.",
             "If all three medians are 19 ms and each spread is narrow, the target proposal adds the native-sample-rounded equivalent of 19 ms to its original onset. That is an estimated trigger-to-patch lag, not proof that the dot's luminance begins exactly then.")
    card("coverage","Will proposed timing change epoch-window coverage?","epoch_coverage.png",
         "X: recording categories, each with original and proposed timing bars. Y: number of analysis-eligible trials, with the same denominator in each pair. A hash-bound final eligibility application is used when present; original Stage 2 eligibility is retained separately. Stacks distinguish inside recording, outside recording and timing unassigned. These are counts; a dispersion bar is not applicable.",
         f"The hashed epoch window is {window[0]} to {window[1]} s. Native rounded endpoints are checked inclusively against the exact recording's first sample and exclusive end. Original and proposed anchors are checked separately; uncertainty coverage refers only to supplied components, not unknown transfer error.",
         f"{len(eligible)} eligible original trials; {new_boundary} newly outside recording after assigned shifts; {sum(r['timing_unassigned'] for r in coverage_by_source)} with timing unassigned. Arithmetic mismatches: {checks['assigned_arithmetic_mismatches']}; eligible sample collisions: {len(checks['eligible_proposed_sample_collisions'])}; ordering inversions: {len(checks['eligible_ordering_inversions'])}.",
         "Review changed coverage, trial-end overlap and unassigned timing before epoching. Preserve Stage 2 exclusions and original keys; explain affected future epoch selection separately. Unknown uncertainty is not an inside-coverage result.",
         "A hypothetical +19-ms shift may move an otherwise eligible epoch beyond a recording's last sample. Flag that trial for the later epoch-window decision rather than silently rewriting its accepted Stage 2 eligibility.")
    for i,(name,source,n) in enumerate(pages,1):
        card(f"waveform-{i}",f"Inspect native edges: {source}",name,
             "X: milliseconds relative to the original recorded STI101 anchor (dashed line at zero). Y: raw diode voltage in mV. The dotted vertical line is the stored candidate 50% edge; the horizontal dotted line, when available, is its polarity-adjusted baseline plus half response amplitude.",
             "Native cached MISC008 samples are multiplied by 1,000; time = (window sample - anchor index) * 1000 / sfreq. No smoothing, filtering or resampling is applied. Per-event baseline method and exact source/sample keys are listed in diagnostic_index.csv and event_measurements.csv.",
             f"Showing {n} retained windows from {source}. Across the report, {len(waveform_index)} windows are shown from {len(waveforms)} available cached/in-memory windows. Selection favors recording/code controls, tails and QC cases; it is not a random population sample.",
             "Check polarity, the intended rise, multiple pulses, flat/attenuated signal, baseline context and threshold position. Match suspicious points to their exact event key before accepting or questioning a timing assignment.",
             "If a marked edge sits on the previous warm-up pulse, the lag is not evidence of this trial's visual onset. If an early edge belongs to a clear current pulse and is stable across thresholds, retain it for scientific review rather than removing it as an outlier.")
    if not pages:
        card("waveform-unavailable","Waveform evidence in this scope","diagnostics.png",
             "No waveform axes or observations are drawn because this bundle retains no native signal windows.",
             "Only existing waveform_windows.npz entries may be displayed; report generation does not open raw FIFs or fabricate a diode response.",
             "No waveform cache entries are available for this report.",
             "Use the acquisition facts and explicit unresolved/fallback status. Missing signal is not evidence of zero physical lag. Additional signal review, when justified and available, requires separately generated evidence.",
             "For a confirmed no-diode recording, a zero applied digital fallback does not represent an observed zero-lag pulse.")
    acquisition = [dict(behavioral_block=block,runtime=value.get("runtime"),
                        refresh_hz=value.get("measuredRefreshHz") or value.get("refrate"),
                        patch_position_px=json.dumps(value.get("photodiodePatchPositionPx")),
                        source_mat=value.get("source_mat")) for block,value in sorted(contexts.items())]
    provenance = dict(source_identity or {}, measurement_version=manifest.get("version"),
                      measurement_mode=manifest.get("mode"), measured_generated_utc=manifest.get("generated_utc"),
                      stage2_report_sha256=manifest.get("stage2_report_sha256"),
                      decisions_record_sha256=manifest.get("decisions_record_sha256"),
                      method=method,software_versions=manifest.get("software_versions",{}),
                      reporting_code_sha256=file_hash(Path(__file__)))
    summary = dict(report_schema_version=REPORT_SCHEMA_VERSION,subject_bids_id=manifest["subject_bids_id"],
                   scientific_acceptance="pending_independent_review",measurement_status_counts=counts,
                   eligible_assignment_counts=bases,original_trial_keys=len(maps),eligible_trial_keys=len(eligible),
                   original_stage2_eligible_trial_keys=original_eligible,
                   lag_statistics=lag_rows,coverage_summary=coverage_by_source,sanity_checks=checks,
                   interval_statistics=interval_source_stats,interval_transition_statistics=interval_pair_stats,
                   provenance=provenance,figure_explanations=cards)
    summary["figure_aliases"] = {"lag_sensitivity.png":"lag_by_trigger.png"}
    if pages:summary["figure_aliases"]["diagnostics.png"] = pages[0][0]
    data["report_summary.json"] = json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)
    # Every table needed for review is embedded as a downloadable exact file.
    # Reports remain usable if copied alone; PNG/CSV companions are convenient.
    title = f"Stage 3 visual timing review: P{subject}"
    lead = (f"{manifest.get('mode','unknown')} evidence: {counts.get('measured',0)} usable measured edges out of {len(rows)} assessed visual events. "
            f"Eligible trial timing: "+", ".join(f"{k}={v}" for k,v in sorted(bases.items()))+". Scientific acceptance remains pending.")
    sections = []
    def section(identifier,heading,body):
        sections.append((identifier,heading,body))
    section("review","What needs a decision", "<p>"+html.escape(lead)+f"</p><p>Original accepted trial keys: {len(maps)}; original Stage 2 eligible: {original_eligible}; original excluded: {len(maps)-original_eligible}. Applicable analysis eligible: {len(eligible)}; analysis excluded: {len(maps)-len(eligible)}. Original anchors and original-policy eligibility/reasons are preserved separately from any final eligibility application. Scientific acceptance of newly included catches does not follow automatically from earlier eligible-only review.</p><ul>"+"".join("<li>"+html.escape(p)+"</li>" for p in guidance)+"</ul>")
    section("assignments","Actual proposed corrections", "<p>Measured shifts are event-specific. Estimated, digital fallback and unresolved assignments remain distinct. A digital fallback applies zero shift while its physical visual latency stays unknown. Dispersion here describes proposed applied shifts; zero SD for a constant estimate/fallback does not establish zero physical uncertainty.</p>"+table_html(corrections,["source_fif","trigger_code","assignment_basis","n_eligible","n_assigned","mean_shift_ms","median_shift_ms","sd_shift_ms","p05_shift_ms","p95_shift_ms"]))
    section("acquisition","Acquisition and saved runtime",table_html(acquisition)+"<p>Native rate, first/end samples and irreversible online filters are retained in recording_summary.csv. Channel presence and patch enablement alone do not establish usable signal. Original Stage 2 flags and reasons remain recorded separately from the narrow final analysis-eligibility application.</p>"+table_html(recordings,["source_fif","acquisition_run","sfreq_hz","first_sample_raw","end_sample_raw_exclusive","online_highpass_hz","diode_channel_present","trigger_code","n_assessed","n_measured","n_flagged"]))
    for item in cards:
        if item["id"]=="interval":
            codes = sorted({r["left_code"] for r in interval_details}|{r["right_code"] for r in interval_details})
            definitions = [dict(trigger_code=code,meaning=TRIGGER_NAMES.get(code,"Other marker")) for code in codes]
            section("interval-reading","Reading the within-trial interval figures",
                "<p>The original bar figure and the additional violin view summarize the same measured interval pairs. Three further panels examine expected length, recording time and endpoint trigger pairs. Recording numbers follow sorted exact source FIF order among sources with comparable pairs; they do not rename acquisition runs.</p>"
                "<p>A pair such as 81-&gt;111 is the separation from that trial's occluded nondeviant onset marker to its occlusion-start marker. Numerical differences between trigger codes carry no timing meaning. Each pair has two recorded optical patch edges; patch timing is not automatically the dot's exact luminance timing.</p>"+table_html(definitions)+
                "<p>Several pairs may share one trial/onset. SD describes pair spread and is not independent-trial uncertainty. These views check internal sequence consistency; a common absolute trigger-to-diode lag cancels between endpoints, so use the separate lag/correction evidence to decide the absolute visual shift.</p>")
        encoded = base64.b64encode(figures[item["file"]]).decode("ascii")
        # A CSS-only native-resolution view keeps scientific labels readable in
        # narrow panels without duplicating image data or requiring JavaScript.
        body = (f"<figure><input type='checkbox' class='figure-zoom' id='zoom-{item['id']}'>"
                f"<label for='zoom-{item['id']}'>View at full resolution (scroll horizontally)</label>"
                f"<div class='figure-scroll'><img src='data:image/png;base64,{encoded}' alt='{html.escape(item['title'],quote=True)}'></div></figure>")
        for key,label in (("axes","Axes and marks"),("derived","How it is derived"),("finding","What these data show"),("decision","How to use it"),("example","Decision example")):
            body += f"<p><strong>{label}.</strong> {html.escape(item[key])}</p>"
        section(item["id"],item["title"],body)
    section("statistics","Exact lag and dispersion statistics", "<p>Usable measured edges only, separated by exact recording, runtime and trigger factor. Sample SD needs n >= 2; observed percentiles describe spread and are not a confidence interval.</p>"+table_html(lag_rows))
    interval_columns = ["source_fif","endpoint_pair","expected_ms","observed_ms","residual_ms","recording_time_sec",
                        "behavioral_block","behavioral_run","behavioral_trial","is_replay_trial",
                        "left_lag_ms","right_lag_ms","saved_frame_missed_count","left_event_key","right_event_key"]
    section("interval-statistics","Exact within-trial interval statistics and event keys",
        "<p>Recording summaries use all comparable measured pairs. Transition summaries retain the endpoint-code factor. n counts pairs; n_trial_contexts counts distinct recorded trial contexts, which can contribute correlated pairs. Complete CSVs are embedded below. No central zoom, new pair-coverage panel or additional rejection threshold is applied.</p>"
        "<h3>Recording residuals (ms)</h3>"+table_html(interval_source_stats)+
        "<details><summary>Endpoint-pair residual means, SD, medians and counts</summary>"+table_html(interval_pair_stats,limit=len(interval_pair_stats))+"</details>"
        "<h3>20 largest absolute residuals: inspection inventory</h3><p>These are interval observations, not 20 independent failed or excluded trials. Read both keys and shared trial context. A saved missed-frame count is a contextual lead, not proof of a cause.</p>"+table_html(interval_tails,interval_columns,limit=20))
    section("signal-statistics","Signal amplitude, noise and sensitivity statistics", "<p>Observed detector quantities, including finite flagged candidates, separated by recording/runtime/trigger factor. Units are explicit. Mean, sample SD and observed percentiles summarize each quantity; they are not calibrated uncertainty bounds for physical onset.</p>"+table_html(signal_rows))
    event_columns = ["event_key","trigger_code","measurement_status","lag_ms","qc_reasons","baseline_method","saved_trial_onset_guard_retry_count"]
    section("event-review","Flagged cases and retained timing tails", "<p>QC cases and stratum extrema are investigation leads. Measured tails remain in the measured statistics; no global central-lag cutoff is used. Full event tables are embedded below.</p><h3>QC and unavailable observations</h3>"+table_html(qc,event_columns)+"<h3>Lowest/highest usable lags within each recording/runtime/code</h3>"+table_html(tails,event_columns))
    section("coverage-tables","Coverage and assignment sanity checks",table_html(coverage_by_source)+"<pre>"+html.escape(json.dumps(checks,indent=2))+"</pre><p>These presentation checks do not revalidate the live raw tree or reopen accepted decisions. The source measurement manifest and Stage 2 hashes establish which evidence was used.</p>")
    if estimates:section("donor-tables","Exact donor estimates and restrictions",table_html(donor_rows)+"<details><summary>Estimator results and exact donor event keys</summary><pre>"+html.escape(json.dumps(estimates,indent=2,sort_keys=True))+"</pre></details>")
    downloads = []
    for name,content in sorted(data.items()):
        raw = content if isinstance(content,bytes) else content.encode()
        mime = "application/json" if name.endswith(".json") else "text/csv" if name.endswith(".csv") else "application/octet-stream"
        downloads.append(f"<li><a download='{html.escape(name,quote=True)}' href='data:{mime};base64,{base64.b64encode(raw).decode()}'>{html.escape(name)}</a></li>")
    section("evidence","Embedded evidence and provenance","<p>Downloads work from this HTML alone. Original measured/estimated timing remains in its immutable source bundle; report-only CSVs are presentation copies.</p><ul>"+"".join(downloads)+"</ul><details><summary>Method, software and source fingerprints</summary><pre>"+html.escape(json.dumps(provenance,indent=2,sort_keys=True))+"</pre></details>")
    checklist = ["Reconcile measured edges, polarity, threshold sensitivity and failed/missing observations.",
                 "Review consistency across trigger factors, recordings/runtimes, drift, context and tails; justify pooling or its limits.",
                 "Reconcile donor QC, estimator/weighting, held-out errors, physical uncertainty and permitted donor restrictions.",
                 "Decide unresolved fallbacks and estimates explicitly; preserve original anchors, exclusions and recording boundaries.",
                 "Review epoch-window changes, collisions/order and uncertainty coverage before later epoching.",
                 "Record the accepted exact measurement bundle and permitted assignment bases separately in the authorized decisions/review records."]
    section("checklist","Independent review checklist","<ul>"+"".join("<li><input type='checkbox'> "+html.escape(p)+"</li>" for p in checklist)+"</ul><p>Checkboxes are a local reading aid; they do not persist decisions or constitute scientific acceptance. Cleaning, ICA, epochs and final inclusion remain later work. Diode-patch onset is not automatically the dot's exact luminance onset.</p>")
    css = """body{font:16px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;color:#202a35;margin:0;background:#f5f6f8}header,main,nav{max-width:1260px;margin:auto;padding:22px 32px}header{background:#fff;border-bottom:4px solid #2563a1}h1{font-size:30px;margin:0 0 10px}h2{font-size:25px;line-height:1.25;margin:0 0 20px}h3{font-size:18px}nav{display:flex;gap:8px 18px;flex-wrap:wrap;font-size:14px}a{color:#155c96}section{background:#fff;margin:0 0 28px;padding:30px;border-top:1px solid #dbe1e8}p{max-width:1100px}li{margin:9px 0}figure{margin:0 0 25px}img{width:100%;height:auto;display:block}.table-scroll{overflow:auto;margin:18px 0}table{border-collapse:collapse;font-size:13px;width:100%}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #dbe1e8;white-space:nowrap}th{background:#eef2f6}pre{overflow:auto;font-size:12px;white-space:pre-wrap;overflow-wrap:anywhere}summary{cursor:pointer;font-weight:600}.scope{font-weight:600;color:#734c13}@media(max-width:650px){header,main,nav{padding:16px}section{padding:18px}h1{font-size:25px}h2{font-size:22px}}@media print{body{background:white}nav{display:none}section{padding:12px;break-before:auto}figure,img{break-inside:avoid}a{color:inherit}}"""
    navigation = "".join(f"<a href='#{identifier}'>{html.escape(heading)}</a>" for identifier,heading,body in sections)
    content = "".join(f"<section id='{identifier}'><h2>{html.escape(heading)}</h2>{body}</section>" for identifier,heading,body in sections)
    css += ".contents{max-width:1260px;margin:0 auto;padding:16px 32px}.contents summary{color:#155c96}.contents nav{padding:16px 0}.figure-scroll{overflow:auto;max-width:100%}.figure-zoom:checked~.figure-scroll img{width:auto;max-width:none}.figure-zoom+label{font-size:13px;color:#155c96;cursor:pointer}@media print{.figure-zoom,.figure-zoom+label{display:none}.figure-zoom:checked~.figure-scroll img{width:100%}}"
    document = f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{html.escape(title)}</title><style>{css}</style></head><body><header><h1>{html.escape(title)}</h1><p class='scope'>{html.escape(str(manifest.get('mode')))} | measurement version {html.escape(str(manifest.get('version')))} | independent review pending</p></header><details class='contents'><summary>Report contents and figure index</summary><nav aria-label='Report sections'>{navigation}</nav></details><main>{content}</main></body></html>\n"
    markdown = [f"# {title}","",lead,"","Open `report.html` for the complete offline report with embedded figures, evidence downloads and detailed tables.","","## Decisions to review",""]
    markdown.extend(f"- {p}" for p in guidance)
    for item in cards:
        markdown.extend(["",f"## {item['title']}","",f"![{item['title']}]({item['file']})",""])
        markdown.extend(f"**{label}.** {item[key]}\n" for key,label in (("axes","Axes and marks"),("derived","How it is derived"),("finding","What these data show"),("decision","How to use it"),("example","Decision example")))
    markdown.extend(["","## Evidence","","See `lag_statistics.csv`, `correction_summary.csv`, `coverage_summary.csv`, `event_qc_review.csv`, `tail_events.csv`, `interval_details.csv`, `interval_statistics.csv`, `interval_transition_statistics.csv`, `interval_tail_pairs.csv`, `estimates.json`, `timing_sanity_checks.json` and `report_summary.json`.","","Scientific acceptance remains pending independent reconciliation.",""])
    return {**data,**figures,"report.html":document,"report.md":"\n".join(markdown)}


# -----------------------------------------------------------------------------
# Immutable report-only publication from a verified existing measurement bundle
# -----------------------------------------------------------------------------


def load_cached_waveforms(folder):
    """Decode only retained native windows; no raw FIF or detector is opened."""
    result = {}
    with np.load(Path(folder)/"waveform_windows.npz",allow_pickle=False) as cache:
        for entry in read_csv(Path(folder)/"diagnostic_index.csv"):
            array = cache[entry["array"]]
            time = array[0]
            zero = np.flatnonzero(np.isclose(time,0,atol=1e-9))
            if len(zero)!=1 or array.shape[0]!=3 or len(time)<2:
                raise RuntimeError("Invalid native waveform cache coordinates.")
            sf = 1000/(time[1]-time[0])
            if sf<=0 or not np.allclose(np.diff(time),1000/sf):
                raise RuntimeError("Invalid native sample grid in waveform cache.")
            result[entry["event_key"]] = (array[1:].copy(),int(zero[0]),float(sf))
    return result


def render_existing_bundle(bundle, output, raw_root):
    """Publish a new presentation with exact source fingerprints and no rereads."""
    from visual_timing import publish
    bundle,output,raw_root = Path(bundle).resolve(),Path(output).resolve(),Path(raw_root).resolve()
    allowed = raw_root.parent/"derivatives/reports/visual_timing"
    if not output.is_relative_to(allowed) or output==allowed or output.exists():
        raise ValueError("Choose a new report directory below derivatives/reports/visual_timing.")
    source_hash = file_hash(bundle/"manifest.json")
    code_hash = file_hash(Path(__file__))
    manifest = verify_bundle(bundle)
    identity = dict(bundle_manifest=str(bundle/"manifest.json"),bundle_manifest_sha256=source_hash,
                    report_generated_utc=utc_now(),raw_signal_windows_read=0)
    inputs = {name:read_csv(bundle/name) for name in ("event_measurements.csv","trial_timing_map.csv","recording_summary.csv","physical_intervals.csv")}
    evidence = {name:(bundle/name).read_bytes() for name in ("event_measurements.csv","event_timing_proposals.csv",
        "trial_timing_map.csv","recording_summary.csv","physical_intervals.csv","estimates.json","contexts.json")}
    evidence["source_manifest.json"] = (bundle/"manifest.json").read_bytes()
    artifacts = build_report_artifacts(inputs["event_measurements.csv"],inputs["trial_timing_map.csv"],
        json.loads((bundle/"estimates.json").read_text()),inputs["recording_summary.csv"],inputs["physical_intervals.csv"],
        json.loads((bundle/"contexts.json").read_text()),manifest,load_cached_waveforms(bundle),source_identity=identity,
        proposals=read_csv(bundle/"event_timing_proposals.csv"),evidence_bytes=evidence)
    # A source edited while rendering cannot acquire a plausible finished report.
    verify_bundle(bundle)
    if file_hash(bundle/"manifest.json")!=source_hash:raise RuntimeError("Measurement source changed during report rendering.")
    if file_hash(Path(__file__))!=code_hash:raise RuntimeError("Report implementation changed during rendering.")
    # Preserve the exact scientific table bytes when copying presentation
    # companions; the report never becomes a replacement measurement bundle.
    for name in ("event_measurements.csv","event_timing_proposals.csv","trial_timing_map.csv",
                 "recording_summary.csv","physical_intervals.csv","estimates.json","contexts.json"):
        artifacts[name] = (bundle/name).read_bytes()
    artifacts["code_visual_timing_report.py.txt"] = Path(__file__).read_bytes()
    report = dict(schema_version=REPORT_SCHEMA_VERSION,status="completed",mode="presentation_report",
                  scientific_acceptance="pending_independent_review",source=identity,
                  subject_bids_id=manifest["subject_bids_id"],reporting_code_sha256=file_hash(Path(__file__)))
    publish(output,artifacts,report,raw_root)
    return report


def main():
    """Build one new, portable participant report from existing hashed evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--raw-root",required=True,type=Path)
    args = parser.parse_args()
    result = render_existing_bundle(args.bundle,args.output,args.raw_root)
    print(json.dumps({"status":result["status"],"source":result["source"],"report":str(args.output/"report.html")},sort_keys=True))


if __name__=="__main__":
    main()
