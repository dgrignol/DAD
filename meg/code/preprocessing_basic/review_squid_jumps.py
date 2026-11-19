#!/usr/bin/env python3
"""Prepare/review Stage 4A SQUID candidates before SSS using accepted channels.

Usage (from the DAD code root)
------------------------------
    python meg/code/preprocessing_basic/review_squid_jumps.py --mode prepare \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml

    python meg/code/preprocessing_basic/review_squid_jumps.py --mode review \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml --source sub04_run01_220526.fif

Read-only inspection while a separate editable review remains open::

    python meg/code/preprocessing_basic/review_squid_jumps.py --mode inspect \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml --source sub04_run01_220526.fif \
      --candidate 20

Inspection reads existing candidates only and writes no file, tuning entry or
review decision. Each flagged sensor has separately scaled native context,
an edge close-up and median-filtered first differences in physical units.
--context also shows available colocated Vectorview sensors as labeled controls.

Inputs are accepted Stage 4A channel evidence, native raw FIFs, and immutable
candidate JSONs. Preparation reads existing accepted lists; it runs no channel
detector, SSS, muscle detector or ICA, and writes no tuning/annotation files.
For a new participant, 03_A_MaxFilt.py prepares the first channel basis and
candidates itself; this helper can then use that separate channel-basis file.

Review uses native unfiltered data. MNE centers each channel's displayed window
and sensor-type scales are estimated only from that window, so remote drift or
large acquired DC offsets do not flatten active traces;
this display-only mean subtraction preserves jumps within the window and never
modifies raw samples or detector inputs.
The browser defaults to the selected candidate's flagged sensors only. ALL
shows other screened sensors for comparison; VIEW/NEXT/PREV returns to focus.
DETAIL opens independent per-sensor context/edge/first-difference plots.
Delete/resize possible_squid_jumps, add BAD_squid_jumps, and click channel names
to mark additional bad channels. VIEW reopens the edited object without reseeding
proposals. SAVE explicitly confirms retained windows, channel scope and review
observations, including zero. DISCARD leaves saved work unchanged. Union is the
channel default; run-only exceptions require an explicit selection and can make
the final ICA channel bases incompatible. The displayed +/-0.5 s proposals are
starting points, not measured artifact recovery durations.

Only SQUID candidate/review sidecars and this source's squid_jumps.reviews tuning
entry are outputs. Saved ocular/muscle annotations, their receipts, active MNE
windows, accepted timing/eligibility and prior SSS files remain untouched.
Scientific Stage 4A acceptance and user-submitted rebuild jobs are separate.
After changed detector settings/channel evidence, append --reconcile to reopen
the previous immutable decisions as unreviewed starting intervals. Native
coordinates must still match. New candidates are listed but not automatically
seeded over prior edits/deletions; inspect them and add intervals deliberately.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mne
import numpy as np

from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import discover_raw_run_records
from squid_jumps import (CANDIDATE_LABEL, LABEL, accepted_report, attach_intervals,
    channel_basis_configuration, file_hash, file_stamp, load_review, parameters,
    prepare, root, save_review, read_channel_basis, remember_channel_basis)


# -----------------------------------------------------------------------------
# Read-only accepted channel basis; fresh preparation is kept outside reports
# -----------------------------------------------------------------------------

def channel_lists(cfg, records):
    """Reuse accepted report lists, or validated fresh pre-SSS channel evidence."""
    report_path = Path(cfg["resolved_paths"]["derivatives_reports_dir"]) / "03_A_MaxFilt_review_report.json"
    if report_path.is_file():
        try:
            report = accepted_report(cfg)
        except RuntimeError:
            # Channel decisions remain known while rebuilt SSS awaits its new
            # acceptance. This helper never grants downstream SSS acceptance.
            return read_channel_basis(cfg, records)["channels"]
        rows = {row["source_filename"]: row.get("squid_jumps", {}).get("screening_bad_channels",
                  row["applied_bad_channels"]) for row in report["runs"]}
        remember_channel_basis(cfg, records, rows, reviewed=True, report=report)
    else:
        basis = read_channel_basis(cfg, records)
        rows = basis["channels"]
    if any(r.source_filename not in rows for r in records):
        raise RuntimeError("SQUID channel evidence does not cover all exact sources.")
    return rows


# -----------------------------------------------------------------------------
# Read-only waveform diagnostics and local physical display scaling
# -----------------------------------------------------------------------------

def candidate_order(raw, proposal, index):
    """Put the chosen event's sensors first, retaining every screened channel.

    MNE's default type sorting can place a candidate's affected sensor outside
    the first visible page. Explicit native channel order avoids that trap.
    """
    names = proposal["screened_channels"]
    affected = proposal["candidates"][index]["channels"] if proposal["candidates"] else []
    order = [name for name in affected if name in names]
    order += [name for name in names if name not in order]
    return [raw.ch_names.index(name) for name in order]


def window_scalings(raw, names, start, duration):
    """Set common physical scales by type from this displayed window only.

    Whole-recording pooled automatic scales can be dominated by slow drift and
    channel DC differences. Use per-channel mean removal and the largest 99th
    percentile amplitude of visible sensors of each type. A flat channel uses
    the same scale as its active peers; it is never inflated to look active.
    """
    sfreq = float(raw.info["sfreq"])
    begin = max(0, min(raw.n_times-1, round(start*sfreq)))
    stop = min(raw.n_times, max(begin+1, round((start+duration)*sfreq)))
    data = raw.get_data(picks=names, start=begin, stop=stop)
    result = {"mag": 1e-12, "grad": 4e-11}
    types = raw.get_channel_types(picks=names)
    for kind in result:
        amplitudes = []
        for row, ch_type in zip(data, types):
            finite = row[np.isfinite(row)]
            if ch_type == kind and finite.size:
                amplitudes.append(float(np.percentile(np.abs(finite-np.mean(finite)), 99)))
        if amplitudes and max(amplitudes) > 0:
            result[kind] = max(amplitudes)*1.1
    return result


def inspection_figures(raw, proposal, index, *, context=False):
    """Build read-only per-sensor context/edge/difference pages for one flag.

    Every axis carries physical units and its own y scale; display centering
    subtracts one constant per sensor, never a trend or filter. Median-filtered
    differences show the detector's edge representation, not a new score or
    hardware diagnosis. Colocated controls need not be flagged. Candidate JSON
    v2 has event-level scores; channel-specific scores are intentionally not
    inferred from the short displayed window. No source annotations are edited.
    """
    import matplotlib.pyplot as plt
    from scipy.ndimage import median_filter
    event = proposal["candidates"][index]
    sfreq = float(raw.info["sfreq"])
    peak = int(event["peak_sample"]-raw.first_samp)
    begin, stop = max(0, peak-round(3*sfreq)), min(raw.n_times, peak+round(3*sfreq)+1)
    flagged = list(event["channels"])
    names = list(flagged)
    if context:
        # Vectorview names end in 1/2/3 for a colocated sensor triplet. Do not
        # invent controls for systems/names that do not have those channels.
        for name in flagged:
            for suffix in "123":
                candidate = name[:-1]+suffix
                if candidate in raw.ch_names and candidate not in names:
                    names.append(candidate)
    types = raw.get_channel_types(picks=names)
    data = raw.get_data(picks=names, start=begin, stop=stop)
    times = np.arange(begin, stop)/sfreq
    width = max(1, round(proposal["identity"]["parameters"]["median_window_sec"]*sfreq))
    width += 1-width % 2
    smooth = median_filter(data, size=(1, width), mode="nearest")
    differences = np.diff(smooth, axis=1)
    edge_time = (event["peak_sample"]-raw.first_samp)/sfreq
    span = [(event[k]-raw.first_samp)/sfreq for k in ("start_sample", "stop_sample")]
    pages, metrics = [], []
    for page_start in range(0, len(names), 6):
        page_names = names[page_start:page_start+6]
        fig, axes = plt.subplots(len(page_names), 3, figsize=(16, max(3.4, len(page_names)*2.4)), squeeze=False)
        for local, name in enumerate(page_names):
            row = page_start+local
            kind = types[row]
            factor, unit = (1e15, "fT") if kind == "mag" else ((1e13, "fT/cm") if kind == "grad" else (1., "SI"))
            trace = data[row]
            finite = trace[:min(round(.2*sfreq), len(trace))]
            center = float(np.median(finite[np.isfinite(finite)])) if np.any(np.isfinite(finite)) else 0.
            original, filtered = (trace-center)*factor, (smooth[row]-center)*factor
            role = "FLAGGED" if name in flagged else ("accepted BAD control" if name in proposal["excluded_channels"] else "context only")
            for col in (0, 1):
                axis = axes[local,col]
                axis.plot(times, original, color="#375a7f", lw=.65, label="Native raw")
                if col == 1:
                    axis.plot(times, filtered, color="#bd531a", lw=1., label=f"Median {width/sfreq*1000:g} ms")
                axis.axvline(edge_time, color="#a51d2d", lw=1., linestyle="--")
                axis.set_xlim((times[0], times[-1]) if col == 0 else
                              (max(times[0],edge_time-.15), min(times[-1],edge_time+.15)))
                if col == 0:
                    axis.axvspan(max(times[0],span[0]), min(times[-1],span[1]), color="#d8b365", alpha=.12)
                else:
                    zoom = (times >= edge_time-.15) & (times <= edge_time+.15)
                    values = np.r_[original[zoom], filtered[zoom]]
                    values = values[np.isfinite(values)]
                    if values.size:
                        margin = max(float(np.ptp(values))*.08, np.finfo(float).eps)
                        axis.set_ylim(float(values.min())-margin, float(values.max())+margin)
                axis.set_ylabel(unit+" (constant centered)")
                axis.grid(alpha=.18)
            axis = axes[local,2]
            axis.plot(times[1:], differences[row]*factor, color="#4b8b3b", lw=.85)
            axis.axvline(edge_time, color="#a51d2d", lw=1., linestyle="--")
            axis.set_xlim(max(times[0],edge_time-.15),min(times[-1],edge_time+.15))
            zoom = (times[1:] >= edge_time-.15) & (times[1:] <= edge_time+.15)
            values = differences[row,zoom]*factor
            values = values[np.isfinite(values)]
            if values.size:
                margin = max(float(np.ptp(values))*.08,np.finfo(float).eps)
                axis.set_ylim(float(values.min())-margin,float(values.max())+margin)
            axis.set_ylabel(unit+" / native sample")
            axis.grid(alpha=.18)
            labels = ("6 s context", "300 ms edge close-up", "Median-filtered first difference")
            for col, label in enumerate(labels):
                axes[local,col].set_title(f"{name} | {role} | {label}",fontsize=9)
                axes[local,col].set_xlabel("Time from recording start (s)")
            offset = peak-begin
            step = float(trace[offset]-trace[offset-1]) if offset > 0 else None
            metrics.append({"channel": name, "role": role, "type": kind, "unit": unit,
                "context_peak_to_peak": float(np.ptp(trace))*factor,
                "raw_sample_change_at_event_peak": step*factor if step is not None else None})
        axes[0,1].legend(fontsize=7)
        described = ', '.join(flagged) if len(flagged) <= 6 else f"{len(flagged)} sensors (complete list in terminal)"
        fig.suptitle(f"{proposal['source_filename']} | candidate {event['id']} at {edge_time:.3f} s | score {event['score']:.2f}\n"
            f"Flagged: {described} | independent y scales, labeled physical units | READ ONLY",
            fontsize=10, wrap=True)
        fig.tight_layout(rect=(0,0,1,.91 if len(page_names)<3 else .95))
        pages.append(fig)
    return pages, metrics


def inspect_candidate(raw, proposal, index, *, context=False):
    """Show diagnostic figures without editing or saving a review decision."""
    import matplotlib.pyplot as plt
    pages, metrics = inspection_figures(raw, proposal, index, context=context)
    event = proposal["candidates"][index]
    print(f"READ ONLY candidate {event['id']} at {event['time_sec']:.3f}s; flagged: {', '.join(event['channels'])}")
    print("Each sensor has its own labeled y scale; trace heights across rows are not comparable.")
    for row in metrics:
        print(f"  {row['channel']} [{row['role']}]: context peak-to-peak {row['context_peak_to_peak']:.1f} {row['unit']}")
    print("Closing these figures saves nothing. Return to the separate editable review to make decisions.")
    plt.show(block=True)
    for fig in pages:
        plt.close(fig)


# -----------------------------------------------------------------------------
# Native browser state survives navigation; only explicit SAVE stores decisions
# -----------------------------------------------------------------------------

def review(cfg, record, excluded, proposal_file, proposal, *, reconcile=False, candidate=None):
    """Review one source without touching existing manual_eye files or receipts."""
    raw = mne.io.read_raw_fif(record.source_path, preload=False, verbose="ERROR")
    baseline_annotations = raw.annotations.copy()
    # Original acquisition labels are protected context. Only new labels in
    # this disposable object become this Stage 4A review's operative intervals.
    original_jump = [(a["onset"], a["duration"]) for a in baseline_annotations if a["description"] == LABEL]
    raw.info["bads"] = list(excluded)
    previous = cfg.get("squid_jumps", {}).get("reviews", {}).get(record.source_filename)
    old = None
    if previous:
        # Changed settings/source need deliberate reconciliation, not a reset.
        try:
            old = load_review(cfg, record, excluded)
        except (RuntimeError, FileNotFoundError):
            if not reconcile:
                raise RuntimeError("Saved SQUID review changed; use --reconcile to retain/recheck previous decisions.")
            name = previous["review_file"]
            if Path(name).name != name or file_hash(root(cfg)/name) != previous["review_sha256"]:
                raise RuntimeError("Previous immutable SQUID review changed; reconcile its provenance first.")
            old = json.loads((root(cfg)/name).read_text())
            prior_proposal = root(cfg) / f"{record.source_stem}-{old['fingerprint']}-candidates.json"
            if file_hash(prior_proposal) != old["proposal_sha256"]:
                raise RuntimeError("Previous proposal changed; cannot recover its reviewed decisions.")
            from continuous_cleaning import coordinates
            if json.loads(prior_proposal.read_text())["identity"]["coordinates"] != coordinates(raw):
                raise RuntimeError("Native coordinates changed; previous intervals cannot be transferred.")
            old["application"]["bad_channels"] = [ch for ch in old["application"]["bad_channels"] if ch not in excluded]
            print("Previous decisions are UNREVIEWED starting points. Inspect the NEW candidate list and entire affected source.")
        attach_intervals(raw, old["application"]["intervals"])
        raw.info["bads"] += old["application"]["bad_channels"]
        print("Reopened saved decisions; dismissed candidates were not reseeded.")
    else:
        # Show overlapping proposal windows once. Repeated transparent fills
        # can obscure the traces and force the operator to delete the same
        # muscle burst many times. Original candidate IDs/evidence stay intact.
        for start_sample,stop_sample in proposal_windows(proposal):
            raw.annotations.append(start_sample / raw.info["sfreq"],
                (stop_sample-start_sample) / raw.info["sfreq"], CANDIDATE_LABEL)
    print(f"{record.source_filename}: {len(proposal['candidates'])} candidates; {len(excluded)} accepted bad channels skipped.")
    for row in proposal["candidates"]:
        print(f"  {row['id']:4d}  {row['time_sec']:9.3f}s  score={row['score']:.1f}  {', '.join(row['channels'])}")
    if proposal["unscored_channels"]:
        print("Unscored channels (inspect separately):", proposal["unscored_channels"])
    if not proposal["screened_channels"]:
        print("No remaining MEG channels; inspect the channel decision before saving reviewed zero.")
    selected_candidate = candidate-1 if candidate is not None else 0
    start = max(0., proposal["candidates"][selected_candidate]["time_sec"]-3.) if proposal["candidates"] else 0.
    show_all = False
    while True:
        if proposal["screened_channels"]:
            names = proposal["screened_channels"]
            event = proposal["candidates"][selected_candidate] if proposal["candidates"] else None
            if event and not show_all:
                names = event["channels"]
            duration = min(6., raw.times[-1])
            title = f"Stage 4A SQUID review: {record.source_filename}"
            if event:
                title += f" | candidate {event['id']} @ {event['time_sec']:.3f}s | "
                title += "ALL screened sensors" if show_all else "FLAGGED sensors only"
                print(f"VIEW candidate {event['id']} @ {event['time_sec']:.3f}s; flagged sensors: {', '.join(event['channels'])}")
            raw.plot(picks=names, start=start, duration=duration,
                     order=[i for i in candidate_order(raw,proposal,selected_candidate) if raw.ch_names[i] in names],
                     n_channels=min(24,len(names)), remove_dc=True,
                     scalings=window_scalings(raw,names,start,duration), block=True, title=title)
        choice = input("VIEW [number] / DETAIL [number] / ALL / NEXT / PREV / SAVE / DISCARD: ").strip().split()
        if not choice:
            continue
        if choice[0].upper() == "DISCARD":
            print("Session discarded; existing reviews and annotations were not changed.")
            return None
        if choice[0].upper() in ("VIEW", "DETAIL", "NEXT", "PREV", "ALL"):
            if len(choice) > 1 and (len(choice) != 2 or not choice[1].isdigit() or
                    not 1 <= int(choice[1]) <= len(proposal["candidates"])):
                print("Use a listed candidate number; no edits changed.")
                continue
            if len(choice) == 2 and choice[1].isdigit() and 1 <= int(choice[1]) <= len(proposal["candidates"]):
                selected_candidate = int(choice[1])-1
                start = max(0., proposal["candidates"][selected_candidate]["time_sec"]-3.)
            elif choice[0].upper() in ("NEXT", "PREV") and proposal["candidates"]:
                direction = 1 if choice[0].upper() == "NEXT" else -1
                selected_candidate = min(max(0,selected_candidate+direction),len(proposal["candidates"])-1)
                start = max(0.,proposal["candidates"][selected_candidate]["time_sec"]-3.)
            show_all = choice[0].upper() == "ALL"
            if choice[0].upper() == "DETAIL" and proposal["candidates"]:
                inspect_candidate(raw,proposal,selected_candidate)
            continue
        if choice[0].upper() != "SAVE":
            print("Choose VIEW, DETAIL, ALL, NEXT, PREV, SAVE or DISCARD.")
            continue
        added_channels = sorted(set(raw.info["bads"])-set(excluded))
        scope = old["application"]["channel_scope"] if old else "union"
        if added_channels:
            entered = input(f"New bad channels {added_channels}; scope union/run [{scope}]: ").strip().lower()
            if entered:
                if entered not in ("union", "run"):
                    print("Choose union or run explicitly.")
                    continue
                scope = entered
        reviewer = input("Reviewer: ").strip()
        reason = input("Review observations (include limitations/unscored channels): ").strip()
        if not reviewer or not reason:
            print("Reviewer and observations are required; nothing saved.")
            continue
        intervals = []
        for ann in raw.annotations:
            if ann["description"] not in (LABEL, CANDIDATE_LABEL):
                continue
            if ann["description"] == LABEL and (ann["onset"], ann["duration"]) in original_jump:
                continue
            a = round(ann["onset"]*raw.info["sfreq"])
            b = round((ann["onset"]+ann["duration"])*raw.info["sfreq"])
            # Discard a seeded interval only when all contributing channels
            # have instead been excluded. Manually adjusted/added spans survive.
            seeded = (a,b) in [tuple(pair) for pair in proposal_windows(proposal)]
            contributors = {ch for e in proposal["candidates"] if e["start_sample"] < b and e["stop_sample"] > a
                            for ch in e["channels"]} if seeded else set()
            if seeded and contributors and contributors.issubset(added_channels):
                continue
            if b > a:
                intervals.append([int(a), int(b)])
        intervals = merge_intervals(intervals)
        dismissed = [e["id"] for e in proposal["candidates"]
                     if not any(a <= e["peak_sample"] < b for a, b in intervals)
                     and not set(e["channels"]).issubset(added_channels)]
        action = {"intervals": intervals, "bad_channels": added_channels, "channel_scope": scope}
        print(f"SAVE: {len(intervals)} {LABEL} intervals, {len(added_channels)} channel additions ({scope}).")
        if input("Confirm this source's SQUID review? [yes/no]: ").strip().lower() != "yes":
            continue
        saved = save_review(cfg, record, excluded, proposal, proposal_file, action,
                            reviewer=reviewer, reason=reason, dismissed=dismissed)
        print("Saved SQUID review only. Existing SSS and manual BAD annotations were not changed.")
        return saved


def merge_intervals(intervals):
    """Union reviewed overlaps without inventing extra rejection padding."""
    merged = []
    for start, stop in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], stop)
        else:
            merged.append([start, stop])
    return merged


def proposal_windows(proposal):
    """Union overlapping display proposals without changing detection evidence."""
    return merge_intervals([[e["start_sample"],e["stop_sample"]] for e in proposal["candidates"]])


def main(argv=None):
    """Prepare native proposals or explicitly review exactly one source."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("prepare", "review", "inspect"), required=True)
    for flag in ("config", "subject-config", "tuning-config"):
        parser.add_argument("--"+flag, required=True)
    parser.add_argument("--site-config")
    parser.add_argument("--source")
    parser.add_argument("--candidate", type=int, help="One-based candidate ID to inspect or open first in a review.")
    parser.add_argument("--context", action="store_true", help="Include available colocated sensors in read-only inspection.")
    parser.add_argument("--reconcile", action="store_true",
        help="Retain/recheck prior SQUID decisions after changed settings/basis; never reset manual work.")
    args = parser.parse_args(argv)
    cfg, *_ = load_and_resolve_config(args.config, args.subject_config, args.tuning_config,
                                    site_config_path=args.site_config, required_sections=("maxwell_filter",))
    parameters(cfg)
    cfg.setdefault("squid_jumps", {})["tuning_path"] = str(Path(args.tuning_config).resolve())
    records = discover_raw_run_records(cfg)
    if args.mode in ("review", "inspect") and not args.source:
        parser.error("Review/inspection requires --source with an exact acquisition filename.")
    if args.mode == "inspect" and args.candidate is None:
        parser.error("Inspection requires --candidate with a listed candidate ID.")
    channels = channel_lists(cfg, records) if args.mode != "inspect" else None
    selected = [r for r in records if not args.source or r.source_filename == args.source]
    if not selected:
        parser.error("Source is not an exact configured acquisition filename.")
    for record in selected:
        if args.mode == "inspect":
            # Inspection is strictly read-only: never call preparation or the
            # channel-basis writer, even while another GUI holds unsaved edits.
            from squid_jumps import ALGORITHM
            candidates = []
            for path in root(cfg).glob(f"{record.source_stem}-*-candidates.json"):
                value = json.loads(path.read_text())
                if value["identity"]["algorithm"] == ALGORITHM and value["identity"]["source"] == file_stamp(record.source_path):
                    from continuous_cleaning import digest
                    if value["content_sha256"] != digest({k:v for k,v in value.items() if k != "content_sha256"}):
                        raise RuntimeError("Candidate file changed; cannot inspect it as intact evidence.")
                    candidates.append((path,value))
            if len(candidates) != 1:
                parser.error("Inspection requires one current proposal; prepare first or reconcile multiple detector settings.")
            path, proposal = candidates[0]
        else:
            path, proposal = prepare(cfg, record, channels[record.source_filename])
        print(f"{record.source_filename}: {len(proposal['candidates'])} candidates; {path}")
        if args.candidate is not None and not 1 <= args.candidate <= len(proposal["candidates"]):
            parser.error("Candidate ID is outside the listed proposal range.")
        if args.mode == "inspect":
            raw = mne.io.read_raw_fif(record.source_path,preload=False,verbose="ERROR")
            inspect_candidate(raw,proposal,args.candidate-1,context=args.context)
        if args.mode == "review":
            review(cfg, record, channels[record.source_filename], path, proposal, reconcile=args.reconcile,candidate=args.candidate)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
