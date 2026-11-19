#!/usr/bin/env python3
"""Render self-contained HTML books from saved annotation-retention tables.

Usage from the workspace (no MEG/ICA loading or detector execution)
----------------------------------------------------------------
    python meg/code/preprocessing_basic/render_annotation_triage_report.py \
      --report-dir /derivatives/annotation_triage/reports/retention-v1 \
      --label retention-v1
    python meg/code/preprocessing_basic/render_annotation_triage_report.py \
      --report-dir /derivatives/annotation_triage/reports/retention-v2 \
      --compare-to /derivatives/annotation_triage/reports/retention-v1 \
      --label retention-v2

After changing annotations, FIRST rerun annotation_triage.py with a new snapshot
label. Rendering an old snapshot cannot incorporate newer annotation edits. The
companion submit_annotation_retention_report.sh can recount and render in one job.

Execution: independent participant books render in spawned processes (--jobs;
default allocated CPUs). Cohort figures, source rechecks and immutable publication
remain ordered in the parent after every participant finishes.

Inputs: report.json and its trial/condition/identity/run CSVs. Counts reconcile
against the trial ledger before plotting. Input-error participants are explicitly
unranked and never plotted as zero. Every complete participant has all 40 identity
cells per condition, including real zero survivors. Baseline comparisons require
the same trial identities, eligibility, conditions and native epoch windows for
each compared participant; changed accounting is not called annotation recovery.

Data flow: small read-only tables -> validated summaries -> standard Matplotlib
PNG figures -> embedded data-URI HTML -> a new immutable visual report directory.
Outputs: annotation_triage/visual_reports/LABEL/index.html, one sub-XX.html per
participant, reusable PNGs, summary CSVs and a SHA-256 provenance manifest. Every
HTML embeds its own figures and needs no image directory, server, JS or network.
Separate participant books are linked from the cohort index. Inputs, previous
snapshots, manual annotations, selections and scientific acceptance stay unchanged.
Counts precede final ICA/amplitude rejection and do not establish clean epochs.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import html
import io
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from participant_parallel import add_parallel_argument, map_participants, participant_workers

SCHEMA = "annotation-triage-visual-report-v1"
TABLES = ("trial_retention", "condition_retention", "identity_retention", "run_condition_retention")
KEEP, LOST = "#247a9d", "#e9b7a8"
CAUTION = ("Preliminary annotation-only counts, before final ICA and amplitude rejection. "
           "Surviving windows are not established clean epochs. No participant replacement threshold is applied.")
COVERAGE = ("Unscored ocular coverage includes disabled detectors (explicit NONE), missing/deferred choices, "
            "and native gaps or short stretches. It does not mean that an artifact is present or absent.")


# -----------------------------------------------------------------------------
# Snapshot validation: errors stay unranked; known zeros remain actual zeros
# -----------------------------------------------------------------------------

def sha(path):
    """Hash a small report artifact, never an original MEG recording."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integer(value):
    """Reject malformed/fractional/negative counts instead of silently coercing."""
    text = str(value)
    if not re.fullmatch(r"\d+", text):
        raise ValueError(f"Invalid nonnegative integer: {value!r}")
    return int(text)


def boolean(value):
    """Require the explicit booleans emitted by the triage CSV writer."""
    if value not in ("True", "False"):
        raise ValueError(f"Invalid CSV boolean: {value!r}")
    return value == "True"


def read_table(path):
    """Read an empty diagnostic table or a complete CSV without changing it."""
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def load_snapshot(path):
    """Load complete published tables, reconcile counts, and retain input hashes."""
    folder = Path(path).resolve(strict=True)
    if folder.parent.name != "reports" or folder.parent.parent.name != "annotation_triage":
        raise ValueError("Input must be an annotation_triage/reports/SNAPSHOT directory.")
    files = [folder / "report.json", *(folder / f"{name}.csv" for name in TABLES)]
    hashes = {p.name: sha(p) for p in files}
    payload = json.loads(files[0].read_text())
    if payload.get("schema") != "annotation-triage-v1" or payload.get("purpose") != "preliminary_annotation_retention_only":
        raise ValueError("Not a preliminary annotation-retention snapshot.")
    summaries = {r["subject_bids_id"]: r for r in payload["subjects"]}
    if len(summaries) != len(payload["subjects"]):
        raise ValueError("Duplicate participant summaries.")
    if any(not re.fullmatch(r"sub-\d{2}", s) for s in summaries):
        raise ValueError("Invalid participant label.")
    if any(r["status"] not in ("counted", "incomplete_inputs") for r in summaries.values()):
        raise ValueError("Unsupported participant input status.")
    tables = {name: read_table(folder / f"{name}.csv") for name in TABLES}
    for name, rows in tables.items():
        for row in rows:
            if row["subject_bids_id"] not in summaries:
                raise ValueError(f"Unknown participant in {name}.")
            if name == "trial_retention":
                for key in ("analysis_eligible", "annotation_survives", "incomplete_window"):
                    row[key] = boolean(row[key])
                row["condition_code"] = integer(row["condition_code"])
                row["trial_sequence_identity"] = integer(row["trial_sequence_identity"]) if row["trial_sequence_identity"] else None
                expected = (row["analysis_eligible"] and not row["incomplete_window"] and not row["bad_categories"])
                if row["annotation_survives"] != bool(expected):
                    raise ValueError("Trial survivor flag disagrees with eligibility/BAD/window evidence.")
                outcome = ("excluded_by_accepted_trial_policy" if not row["analysis_eligible"] else
                    "incomplete_window" if row["incomplete_window"] else
                    "overlaps_BAD" if row["bad_categories"] else "annotation_survivor")
                if row["triage_outcome"] != outcome:
                    raise ValueError("Exclusive trial outcome disagrees with its evidence.")
            else:
                for key in ("eligible", "annotation_survivors"):
                    row[key] = integer(row[key])
                if row["annotation_survivors"] > row["eligible"]:
                    raise ValueError("Survivors exceed eligible trials.")
                if name == "identity_retention":
                    row["identity"] = integer(row["identity"])
    snapshot = {"folder": folder, "root": folder.parent.parent, "report": payload,
                "summaries": summaries, "hashes": hashes, **tables, "participants": {}}
    for subject, summary in summaries.items():
        if summary["status"] != "counted":
            continue  # Partial diagnostic tables never masquerade as complete data.
        part = {name: [r for r in rows if r["subject_bids_id"] == subject] for name, rows in tables.items()}
        conditions = sorted(part["condition_retention"], key=lambda r: integer(r["condition_code"]))
        if len(conditions) != 3 or len({r["condition"] for r in conditions}) != 3:
            raise ValueError(f"{subject}: expected three distinct conditions.")
        trial_keys = [(r["source_fif"], r["trial_key"]) for r in part["trial_retention"]]
        if len(trial_keys) != len(set(trial_keys)):
            raise ValueError(f"{subject}: duplicate source/trial identity.")
        codes = {integer(c["condition_code"]) for c in conditions}
        if len(codes) != 3:
            raise ValueError("Condition codes must be distinct.")
        if any(t["condition_code"] not in codes for t in part["trial_retention"]):
            raise ValueError(f"{subject}: unexpected trial condition.")
        for condition in conditions:
            name, code = condition["condition"], integer(condition["condition_code"])
            trials = [r for r in part["trial_retention"] if r["condition_code"] == code]
            ids = sorted([r for r in part["identity_retention"] if r["condition"] == name], key=lambda r: r["identity"])
            if [r["identity"] for r in ids] != list(range(1, 41)):
                raise ValueError(f"{subject} {name}: require all 40 unique identity cells, including zeros.")
            if any(t["analysis_eligible"] and t["trial_sequence_identity"] not in range(1, 41) for t in trials):
                raise ValueError("An eligible trial has no valid identity; do not invent its cell.")
            reconcile(condition, trials)
            for cell in ids:
                reconcile(cell, [t for t in trials if t["trial_sequence_identity"] == cell["identity"]])
            runs = [r for r in part["run_condition_retention"] if r["condition"] == name]
            if len({r["source_fif"] for r in runs}) != len(runs):
                raise ValueError("Duplicate run/condition summary.")
            if {r["source_fif"] for r in part["trial_retention"]} != {r["source_fif"] for r in runs}:
                raise ValueError("Run table and trial source coverage differ.")
            for run in runs:
                reconcile(run, [t for t in trials if t["source_fif"] == run["source_fif"]])
        if len(part["identity_retention"]) != 120 or len(part["run_condition_retention"]) != 3 * len({r["source_fif"] for r in part["trial_retention"]}):
            raise ValueError("Unexpected identity/run condition cells.")
        zeros = sum(r["annotation_survivors"] == 0 for r in part["identity_retention"])
        if (summary["minimum_condition_survivors"] != min(c["annotation_survivors"] for c in conditions)
                or summary["total_annotation_survivors"] != sum(c["annotation_survivors"] for c in conditions)
                or summary["missing_identity_cells"] != zeros):
            raise ValueError("Priority summary disagrees with trial/identity counts.")
        snapshot["participants"][subject] = {**part, "conditions": conditions, "summary": summary}
    snapshot["ranked"] = sorted(snapshot["participants"], key=lambda s: (
        summaries[s]["minimum_condition_survivors"], -summaries[s]["missing_identity_cells"], s))
    return snapshot


def reconcile(cell, trials):
    """Cross-check each plotted table cell against exact trial-level outcomes."""
    if (cell["eligible"] != sum(t["analysis_eligible"] for t in trials)
            or cell["annotation_survivors"] != sum(t["annotation_survives"] for t in trials)):
        raise ValueError("Summary counts disagree with the trial ledger.")


def compare_snapshots(current, previous):
    """Compare only complete common participants with unchanged trial accounting."""
    comparable = set(current["participants"]) & set(previous["participants"])
    contract = ("condition_code", "trial_sequence_identity", "analysis_eligible",
                "sfreq_hz", "recording_first_sample_raw", "recording_end_sample_raw_exclusive",
                "epoch_event_sample_raw", "epoch_start_raw", "epoch_stop_raw_exclusive")
    for subject in comparable:
        def keys(snapshot):
            return {(t["source_fif"], t["trial_key"]): tuple(t.get(k, "") for k in contract)
                    for t in snapshot["participants"][subject]["trial_retention"]}
        if keys(current) != keys(previous):
            raise ValueError(f"{subject}: trial identities/eligibility/epoch windows changed; annotation-only comparison is invalid.")
        a = [(c["condition"], c["condition_code"]) for c in current["participants"][subject]["conditions"]]
        b = [(c["condition"], c["condition_code"]) for c in previous["participants"][subject]["conditions"]]
        if a != b:
            raise ValueError("Condition definitions changed between snapshots.")
    return comparable


# -----------------------------------------------------------------------------
# Scientific plots: shared identity scales, explicit denominators and overlaps
# -----------------------------------------------------------------------------

def figure(width, height):
    """Use an isolated Agg canvas; never change the user's GUI backend."""
    fig = Figure(figsize=(width, height), layout="constrained")
    FigureCanvasAgg(fig)
    return fig


def identity_figure(part, previous=None):
    """Show identities 1-40 per condition; stacked totals equal eligible counts."""
    fig = figure(14, 9)
    axes = fig.subplots(3, 1, sharex=True, sharey=True)
    maximum = max([r["eligible"] for r in part["identity_retention"]] + [1])
    if previous:
        maximum = max(maximum, max(r["annotation_survivors"] for r in previous["identity_retention"]))
    for ax, condition in zip(axes, part["conditions"]):
        cells = sorted([r for r in part["identity_retention"] if r["condition"] == condition["condition"]], key=lambda r: r["identity"])
        x = np.arange(1, 41)
        kept = np.array([r["annotation_survivors"] for r in cells])
        lost = np.array([r["eligible"] for r in cells]) - kept
        ax.bar(x, kept, width=.8, color=KEEP, label="Survives annotations")
        ax.bar(x, lost, bottom=kept, width=.8, color=LOST, label="Lost to BAD/window boundary")
        if previous:
            old = {r["identity"]: r["annotation_survivors"] for r in previous["identity_retention"] if r["condition"] == condition["condition"]}
            ax.scatter(x, [old[i] for i in x], s=18, facecolors="none", edgecolors="#202b33", zorder=3, label="Previous survivors")
        ax.set_title(f'{condition["condition"]} | kept {int(kept.sum())}/{condition["eligible"]} eligible | {int((kept == 0).sum())}/40 empty identity cells', loc="left", fontsize=11)
        ax.set(ylabel="Trials", ylim=(0, maximum * 1.15 + .3), xlim=(.3, 40.7))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
        ax.grid(axis="y", alpha=.18); ax.set_axisbelow(True)
    axes[0].legend(loc="upper right", fontsize=8, ncol=3 if previous else 2)
    axes[-1].set(xticks=range(1, 41), xlabel="Stimulus identity")
    fig.suptitle(f'{part["summary"]["subject_bids_id"]}: preliminary trials per identity and condition', fontsize=14)
    return fig


def heatmap(values, row_names, column_names, title, *, percent=False, delta=False, annotations=None, width=10):
    """Render unknown values as grey; true zeros stay labelled numeric zeros."""
    values = np.array(values, dtype=float)
    fig = figure(width, max(3, .31 * len(row_names) + 1.6))
    ax = fig.subplots()
    cmap = __import__("matplotlib").colormaps["RdBu" if delta else "viridis"].copy()
    cmap.set_bad("#e6e8eb")
    scale = max(1, float(np.nanmax(np.abs(values)))) if np.isfinite(values).any() else 1
    artist = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=cmap,
        vmin=-scale if delta else 0, vmax=scale if delta else 100 if percent else scale)
    ax.set(yticks=range(len(row_names)), yticklabels=row_names,
        xticks=range(len(column_names)), xticklabels=column_names, title=title)
    if len(column_names) > 10:
        ax.tick_params(axis="x", labelsize=6)
    else:
        for i, row in enumerate(values):
            for j, value in enumerate(row):
                text = "n/a" if not np.isfinite(value) else annotations[i][j] if annotations else f"{value:+.0f}" if delta else f"{value:.1f}%" if percent else str(int(value))
                white = np.isfinite(value) and (.22 < float(artist.norm(value)) < .72 if delta else float(artist.norm(value)) < .55)
                ax.text(j, i, text, ha="center", va="center", color="white" if white else "#1b2530", fontsize=8)
    colorbar = fig.colorbar(artist, ax=ax, label="Change in survivors" if delta else "Retention (%)" if percent else "Surviving trials", shrink=.8)
    if not percent:
        colorbar.locator = MaxNLocator(integer=True, nbins=6)
        colorbar.update_ticks()
    return fig


def reason_counts(trials):
    """Count each reason once per eligible trial; causes deliberately overlap."""
    counts = Counter()
    for trial in trials:
        if trial["analysis_eligible"]:
            counts.update(set(filter(None, trial["bad_categories"].split(";"))))
            if trial["incomplete_window"]:
                counts["incomplete_window"] += 1
    return counts


def reason_figure(trials, title):
    """Do not stack overlapping BAD reasons or imply they sum to unique loss."""
    counts = reason_counts(trials)
    fig = figure(10, max(3, len(counts) * .36 + 1.7))
    ax = fig.subplots()
    names = sorted(counts, key=lambda n: (-counts[n], n))
    if names:
        ax.barh(names, [counts[n] for n in names], color="#b76551")
        ax.invert_yaxis()
        ax.bar_label(ax.containers[0], padding=3, fontsize=9)
        ax.margins(x=.15)
    else:
        ax.text(.5, .5, "No BAD overlaps or incomplete windows", ha="center", transform=ax.transAxes)
    ax.set(title=title + "\nOverlapping causes: do not add these counts", xlabel="Affected eligible trial windows")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    return fig


# -----------------------------------------------------------------------------
# Portable HTML and derivative summaries; every image is embedded in its book
# -----------------------------------------------------------------------------

def table(headers, rows):
    """Escape source identifiers and prose before embedding them in HTML."""
    esc = lambda value: html.escape(str(value))
    return '<div class="table-wrap"><table><thead><tr>' + ''.join(f'<th>{esc(h)}</th>' for h in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join(f'<td>{esc(v)}</td>' for v in row) + '</tr>' for row in rows) + '</tbody></table></div>'


def document(title, content):
    """Create a standalone readable book with print styles and no remote assets."""
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title><style>
body{{font:16px/1.55 system-ui,sans-serif;color:#25313d;background:#f3f5f7;margin:0}}
main{{max-width:1500px;margin:auto;padding:28px}}h1,h2,h3{{line-height:1.2}}
section,figure,.notice{{background:white;border:1px solid #dce2e8;border-radius:10px;padding:20px;margin:20px 0}}
.notice{{border-left:5px solid #c38539}}img{{display:block;width:100%;height:auto}}figure{{margin:24px 0}}
figcaption,.muted{{color:#566474;font-size:14px}}table{{border-collapse:collapse;width:100%;font-size:14px}}
th,td{{padding:9px 12px;text-align:right;border-bottom:1px solid #dce2e8}}th:first-child,td:first-child{{text-align:left}}
th{{background:#edf2f6}}.table-wrap{{overflow:auto}}a{{color:#176b91}}.links{{display:flex;flex-wrap:wrap;gap:12px}}
@media print{{body{{background:white}}main{{padding:0}}section,figure{{break-inside:avoid;border:0}}}}
</style></head><body><main><h1>{html.escape(title)}</h1>
<div class="notice">{CAUTION}<p>{COVERAGE}</p></div>{content}</main></body></html>'''


def embed(fig, folder, filename, caption):
    """Write an exportable PNG and embed identical bytes in the HTML figure."""
    buffer = io.BytesIO(); fig.savefig(buffer, format="png", dpi=135)
    data = buffer.getvalue(); fig.clear()
    (folder / "figures" / filename).write_bytes(data)
    encoded = base64.b64encode(data).decode("ascii")
    return f'<figure><img alt="{html.escape(caption, quote=True)}" src="data:image/png;base64,{encoded}"><figcaption>{html.escape(caption)}</figcaption></figure>'


def percent(kept, eligible):
    """An absent denominator is unknown, never 0% retention."""
    return 100 * kept / eligible if eligible else np.nan


def summary_rows(snapshot):
    """Produce priority-ordered condition summaries with identity and coverage gaps."""
    result = []
    for rank, subject in enumerate(snapshot["ranked"], 1):
        part = snapshot["participants"][subject]
        for condition in part["conditions"]:
            cells = [r["annotation_survivors"] for r in part["identity_retention"] if r["condition"] == condition["condition"]]
            trials = [r for r in part["trial_retention"] if r["condition_code"] == integer(condition["condition_code"])]
            result.append({"priority": rank, "participant": subject, "condition": condition["condition"],
                "eligible": condition["eligible"], "survivors": condition["annotation_survivors"],
                "lost": condition["eligible"] - condition["annotation_survivors"],
                "retention_percent": percent(condition["annotation_survivors"], condition["eligible"]),
                "empty_identities": cells.count(0), "minimum_per_identity": min(cells),
                "median_per_identity": float(np.median(cells)), "maximum_per_identity": max(cells),
                "survivors_unscored_blink": sum(t["annotation_survives"] and "blink" in t["unscored_groups"].split(";") for t in trials),
                "survivors_unscored_eye_movement": sum(t["annotation_survives"] and "eye_movement" in t["unscored_groups"].split(";") for t in trials)})
    return result


def participant_book(subject, part, rows, folder, previous=None):
    """One participant book: identity bars, run retention, causes and exact counts."""
    content = '<p><a href="index.html">Cohort overview</a></p>'
    content += '<section><h2>Condition summary</h2>' + table(
        ["Condition", "Eligible", "Survive", "Lost", "Retention", "Empty identities", "Min/median/max per identity", "Survivors: blink unscored", "Survivors: eye unscored"],
        [[r["condition"], r["eligible"], r["survivors"], r["lost"], f'{r["retention_percent"]:.1f}%' if np.isfinite(r["retention_percent"]) else "n/a", f'{r["empty_identities"]}/40',
          f'{r["minimum_per_identity"]}/{r["median_per_identity"]:g}/{r["maximum_per_identity"]}', r["survivors_unscored_blink"], r["survivors_unscored_eye_movement"]] for r in rows]) + '</section>'
    content += embed(identity_figure(part, previous), folder, f"{subject}_identities.png", "Bar heights show eligible trials: blue survives; pale red is lost to BAD intervals or incomplete windows. All three panels share a y scale." + (" Hollow points show the previous snapshot's survivors." if previous else ""))
    conditions = [r["condition"] for r in part["conditions"]]
    sources = sorted({r["source_fif"] for r in part["run_condition_retention"]})
    run_cells = {(r["source_fif"], r["condition"]): r for r in part["run_condition_retention"]}
    values = [[percent(run_cells[s, c]["annotation_survivors"], run_cells[s, c]["eligible"]) for c in conditions] for s in sources]
    labels = [[f'{run_cells[s,c]["annotation_survivors"]}/{run_cells[s,c]["eligible"]}' for c in conditions] for s in sources]
    content += embed(heatmap(values, sources, conditions, f"{subject}: run retention (survive / eligible)", percent=True, annotations=labels, width=13), folder, f"{subject}_runs.png", "Run colours show retention percentage; cell labels show surviving/eligible trials. Grey cells have no eligible denominator.")
    content += embed(reason_figure(part["trial_retention"], f"{subject}: rejection reasons"), folder, f"{subject}_reasons.png", "A trial can overlap multiple BAD categories and a recording boundary. Reason counts overlap; unique loss is shown in the condition summary.")
    outcomes = Counter((t["condition_code"], t["triage_outcome"]) for t in part["trial_retention"])
    content += '<section><h2>Trial accounting</h2>' + table(["Condition", "Policy excluded", "Incomplete window", "BAD overlap only", "Survives"],
        [[c["condition"], *(outcomes[integer(c["condition_code"]), key] for key in ("excluded_by_accepted_trial_policy", "incomplete_window", "overlaps_BAD", "annotation_survivor"))] for c in part["conditions"]]) + '<p>These outcome columns are mutually exclusive. Incomplete-window trials may also overlap BAD intervals.</p></section>'
    ids = {(r["condition"], r["identity"]): r for r in part["identity_retention"]}
    content += '<section><h2>Exact identity counts (survive / eligible)</h2>' + table(["Identity", *conditions], [[i, *(f'{ids[c,i]["annotation_survivors"]}/{ids[c,i]["eligible"]}' for c in conditions)] for i in range(1,41)]) + '</section>'
    methods = Counter(m["method"] for m in part["summary"].get("masks", []))
    content += '<section><h2>Annotation basis</h2>' + table(["Mask method", "Acquisitions"], sorted(methods.items())) + '<p>Manual masks preserve saved review decisions. Automatic masks use provisional component choices and detector settings. EDF/SQUID proposals are outside this triage pass.</p></section>'
    (folder / f"{subject}.html").write_text(document(f"{subject}: preliminary annotation retention", content))


def cohort_book(snapshot, folder, rows, previous, comparable):
    """Show review priority and cohort coverage without hiding incomplete inputs."""
    ranked = snapshot["ranked"]
    content = f'<p>Snapshot: <strong>{html.escape(snapshot["folder"].name)}</strong> | counted: {len(ranked)}/{len(snapshot["summaries"])} participants.</p>'
    content += '<p class="muted">This report is a saved snapshot. After changing annotations, recompute counts before rendering a new report.</p>'
    if ranked:
        conditions = [r["condition"] for r in snapshot["participants"][ranked[0]]["conditions"]]
        if any([c["condition"] for c in snapshot["participants"][s]["conditions"]] != conditions for s in ranked):
            raise ValueError("Cohort condition definitions differ.")
        cells = {(r["participant"], r["condition"]): r for r in rows}
        for key, title, is_percent in (("survivors", "Surviving trials by condition", False), ("retention_percent", "Annotation retention by condition", True)):
            values = [[cells[s,c][key] for c in conditions] for s in ranked]
            content += embed(heatmap(values, ranked, conditions, title, percent=is_percent, width=11), folder, f"cohort_{key}.png", "Participants are ordered by the smallest surviving condition count, then empty identity cells. Input-error participants are unranked.")
        values = [[r["annotation_survivors"] for c in conditions for r in sorted([r for r in snapshot["participants"][s]["identity_retention"] if r["condition"] == c], key=lambda r:r["identity"])] for s in ranked]
        fig = heatmap(values, ranked, [str(i) for _ in conditions for i in range(1,41)], "Identity coverage: all 120 condition-by-identity cells", width=19)
        ax = fig.axes[0]
        for boundary in (39.5,79.5): ax.axvline(boundary, color="white", lw=2)
        # A real secondary axis participates in constrained layout; free text
        # above the image can be clipped, especially for one-participant reports.
        top = ax.secondary_xaxis("top")
        top.set_xticks([19.5,59.5,99.5], conditions)
        top.tick_params(length=0, labelsize=9, pad=4)
        top.spines["top"].set_visible(False)
        ax.set_title("Identity coverage: all 120 condition-by-identity cells", pad=28)
        content += embed(fig, folder, "cohort_identities.png", "Each block contains identities 1-40 for one condition. Dark cells are real low/zero survivor counts; every cell is included.")
        all_trials = [t for s in ranked for t in snapshot["participants"][s]["trial_retention"]]
        content += embed(reason_figure(all_trials, "Cohort: rejection reasons"), folder, "cohort_reasons.png", "Affected eligible windows across fully counted participants. Categories overlap and must not be added.")
        content += '<section><h2>Review priority</h2><p>Start with the smallest surviving condition count. Retention rates and identity coverage provide context; this is not an exclusion rule.</p>' + table(
            ["Priority", "Participant", "Weakest condition", "Total survivors", "Empty cells /120", "Missing/deferred selections", "Survivors with unscored ocular coverage"],
            [[i,s,snapshot["summaries"][s]["minimum_condition_survivors"],snapshot["summaries"][s]["total_annotation_survivors"],snapshot["summaries"][s]["missing_identity_cells"],snapshot["summaries"][s]["runs_missing_selections"],snapshot["summaries"][s]["survivors_with_unscored_ocular"]] for i,s in enumerate(ranked,1)]) + '</section>'
    content += '<section><h2>Participant reports</h2><div class="links">' + ''.join(f'<a href="{s}.html">{s}</a>' for s in sorted(snapshot["summaries"])) + '</div></section>'
    if previous:
        content += f'<section><h2>Comparison with {html.escape(previous["folder"].name)}</h2><p>Only fully counted participants with identical trial identities, eligibility and native epoch windows are compared. Positive changes mean more surviving windows.</p>'
        if comparable:
            names = [s for s in ranked if s in comparable]
            values = [[snapshot["participants"][s]["conditions"][i]["annotation_survivors"] - previous["participants"][s]["conditions"][i]["annotation_survivors"] for i in range(3)] for s in names]
            content += embed(heatmap(values, names, conditions, "Change in surviving trials", delta=True, width=11), folder, "cohort_changes.png", "Current minus previous survivors. Annotation revisions can recover trials or reject additional contaminated trials.")
        excluded = set(snapshot["summaries"]) | set(previous["summaries"])
        content += '<p>Not comparable (absent or incomplete inputs): ' + html.escape(', '.join(sorted(excluded-comparable)) or 'none') + '</p></section>'
    errors = snapshot["report"].get("errors", [])
    content += '<section><h2>Input errors</h2>' + (table(["Participant","Source","Error"], [[e["subject_bids_id"],e["source_fif"],e["error"]] for e in errors]) if errors else '<p>No input errors reported.</p>') + '</section>'
    (folder / "index.html").write_text(document("Preliminary annotation retention", content))


# -----------------------------------------------------------------------------
# Immutable publication, input rechecks and the reusable command-line interface
# -----------------------------------------------------------------------------

def write_csv(path, rows):
    """Export derived summaries for reproducible quantitative comparisons."""
    with path.open("w", newline="") as handle:
        if rows:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)


def render_participant(task):
    """Write one participant's embedded book and PNGs in an isolated process."""
    subject, summary, part, rows, old, folder, errors = task
    print(f"{subject}: rendering HTML and embedded figures ({summary['status']})",flush=True)
    if summary["status"] != "counted":
        errors = [e for e in errors if e["subject_bids_id"] == subject]
        content = '<p><a href="index.html">Cohort overview</a></p><section><h2>Incomplete inputs: unranked</h2><p>Diagnostic partial counts are not displayed as participant retention or zero survivors.</p>' + table(["Source","Error"],[[e["source_fif"],e["error"]] for e in errors]) + '</section>'
        (folder/f"{subject}.html").write_text(document(subject+": incomplete inputs",content))
        return
    participant_book(subject,part,rows,folder,old)


def render_report(report_dir, *, label=None, compare_to=None, jobs=None):
    """Publish a new book directory only after validating every source and output."""
    current = load_snapshot(report_dir)
    previous = load_snapshot(compare_to) if compare_to else None
    comparable = compare_snapshots(current, previous) if previous else set()
    label = label or current["folder"].name + "-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",label) or label in (".",".."):
        raise ValueError("--label must be a simple folder name.")
    parent = current["root"] / "visual_reports"
    if parent.resolve() != parent:
        raise ValueError("Visual output must not redirect through a symlink.")
    parent.mkdir(exist_ok=True)
    target = parent / label
    if target.exists():
        raise FileExistsError(f"Visual report already exists; choose a new label: {target}")
    folder = Path(tempfile.mkdtemp(prefix=".building-"+label+"-", dir=parent))
    try:
        (folder / "figures").mkdir()
        rows = summary_rows(current)
        workers = participant_workers(jobs, len(current["summaries"]))
        print(f"Rendering {len(current['summaries'])} participants with {workers} participant worker(s)", flush=True)
        tasks = [(subject, summary, current["participants"].get(subject),
            [r for r in rows if r["participant"] == subject],
            previous["participants"][subject] if subject in comparable else None,
            folder, current["report"].get("errors", []))
            for subject, summary in sorted(current["summaries"].items())]
        map_participants(render_participant, tasks, workers)
        # The shared cohort index and immutable manifest have exactly one writer.
        cohort_book(current,folder,rows,previous,comparable)
        write_csv(folder/"condition_summary.csv", rows)
        changes=[]
        for subject in sorted(comparable):
            for now,old in zip(current["participants"][subject]["conditions"],previous["participants"][subject]["conditions"]):
                changes.append({"participant":subject,"condition":now["condition"],"previous_survivors":old["annotation_survivors"],"current_survivors":now["annotation_survivors"],"change":now["annotation_survivors"]-old["annotation_survivors"]})
        write_csv(folder/"changes.csv", changes)
        inputs=[]
        for snapshot in (current,previous):
            if snapshot:
                for name,expected in snapshot["hashes"].items():
                    if sha(snapshot["folder"]/name)!=expected:
                        raise RuntimeError("Source snapshot changed during rendering; no report published.")
                inputs.append({"folder":str(snapshot["folder"]),"sha256":snapshot["hashes"]})
        outputs={str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob("*")) if p.is_file()}
        manifest={"schema":SCHEMA,"purpose":"preliminary_annotation_retention_only","generated_utc":datetime.now(timezone.utc).isoformat(),"execution":{"participant_workers":workers},"inputs":inputs,"counted_participants":current["ranked"],"compared_participants":sorted(comparable),"output_sha256":outputs}
        (folder/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        os.rename(folder,target)
    except BaseException:
        shutil.rmtree(folder)
        raise
    print(f"Saved self-contained HTML report: {target/'index.html'}",flush=True)
    return target


def main(argv=None):
    """Render saved counts; detector recomputation is an explicit separate step."""
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    add_parallel_argument(parser)
    parser.add_argument("--report-dir",required=True)
    parser.add_argument("--compare-to",help="Previous annotation snapshot, for compatible before/after counts.")
    parser.add_argument("--label",help="New immutable directory under annotation_triage/visual_reports.")
    args=parser.parse_args(argv)
    render_report(args.report_dir,label=args.label,compare_to=args.compare_to,jobs=args.jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
