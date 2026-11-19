#!/usr/bin/env python3
"""Build offline, navigable ocular-component reports, ONE HTML per participant.

Usage from the DAD workspace
---------------------------
    python meg/code/preprocessing_basic/build_ocular_component_reports.py \
      --subjects 14 --label ocular-20261030
    python meg/code/preprocessing_basic/build_ocular_component_reports.py \
      --subjects 14 15 --label ocular-20261030 --resume
    # Repair only the navigation index after report jobs have finished:
    python meg/code/preprocessing_basic/build_ocular_component_reports.py \
      --subjects 14 --label ocular-20261030 --index-only
    # Cluster: add --site-config meg/params/pipeline_cluster.yaml.
    # The companion submission wrapper exports remaining participants in batch.

Inputs: existing accepted 4A evidence, report-bound muscle context and compatible
preliminary frontal ICA caches. Completed manual work and saved component choices
are preserved and listed, not regenerated or modified. Every pending acquisition
gets maps and a compact peak-preserving source view. All pending runs are embedded
in one participant HTML with a run selector, time slider, zoom and gain controls.
The file is offline/self-contained: viewing reads no FIF, server or external JS.

Data flow: read-only source/cache validation -> mag/grad maps and diagnostic trace
copies -> content-bound internal run assets -> participant HTML/manifest -> cohort
index and blank decision template. Run assets only support generation/resume;
users open the participant HTML. Existing cached ICA is never refitted.

Trace reduction chooses the largest absolute NATIVE sample in each display bin,
storing its exact within-bin offset. Any native NaN/nonfinite gap makes that bin
unscored. Normalize by each native source's maximum absolute finite amplitude.
These are view-only copies in arbitrary units, not inputs to detection. All
model/source identities remain bound; the importer revalidates native caches
before accepting choices. No annotations, choices, tuning or acceptance are
written here. New report products stay under derivatives/annotation_triage.
"""

from __future__ import annotations

import argparse
import base64
import fcntl
import hashlib
import html
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
from uuid import uuid4

import mne
import numpy as np

from annotation_triage import (common_arguments, component_binding, csv_content,
    has_manual_work, selection_path, subject_config, triage_path, triage_root)
from blink_annotation import cache_paths, load_cached_aid
from continuous_cleaning import digest, file_hash, file_stamp, require_review
from preprocessing_stage_utils import (derivative_path_for_run, discover_raw_run_records,
    stage_output_dirs, utc_now)

REPORT_SCHEMA = "ocular-component-reports-v1"
REPORT_PURPOSE = "preliminary_ocular_component_review_only"
VIEWER = Path(__file__).with_name("ocular_component_report.js")
TEMPLATE = Path(__file__).with_name("ocular_component_report.html")


# -----------------------------------------------------------------------------
# Fingerprints and atomic report-only outputs; never create scientific decisions
# -----------------------------------------------------------------------------

def input_fingerprints(cfg, record):
    """Bind all cache-validation dependencies without repeatedly opening large FIFs."""
    paths = cache_paths(cfg, record)
    dirs = stage_output_dirs(cfg)
    roles = {**{f"config_{k}": (Path(v) if v else None)
                for k, v in cfg.get("_config_sources", {}).items()},
        "cache_metadata": paths["metadata"], "cache_model": paths["model"],
        "cache_sources": paths["sources"], "raw_source": Path(record.source_path),
        "sss": derivative_path_for_run(record, dirs["sss"], "_raw_sss.fif"),
        "annotated": derivative_path_for_run(record, dirs["sss"], "_ann.fif"),
        "artifact_report": dirs["reports"] / "03_B_ArtefactAnnotation_review_report.json",
        "sss_report": dirs["reports"] / "03_A_MaxFilt_review_report.json",
        "manual": Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) /
            "manual_eye" / f"{record.source_stem}-annot.fif"}
    result = {}
    for role, path in roles.items():
        if path is None or not path.exists():
            result[role] = None
        elif role in ("sss", "annotated", "cache_sources", "raw_source"):
            result[role] = file_stamp(path)
        else:
            result[role] = {"stamp": file_stamp(path), "sha256": file_hash(path)}
    return result


def atomic_bytes(cfg, path, content, *, immutable=False, locked_index_fallback=False):
    """Prefer atomic publication; permit a narrowly scoped locked index repair.

    Some shared-storage clients can write an existing index but cannot replace
    its directory entry across SMB/NFS ownership or sharing boundaries. Only
    update_index, while holding the batch lock, may fall back to an in-place
    update of index.html/decisions_template.csv. Immutable books, manifests and
    assets never use this fallback. No permission/ownership changes are made.
    """
    path = triage_path(cfg, path)
    if locked_index_fallback and (immutable or path.name not in ("index.html", "decisions_template.csv")):
        raise ValueError("In-place fallback is restricted to locked mutable index products.")
    if immutable and path.exists():
        if path.read_bytes() != content:
            raise RuntimeError(f"Existing report differs; use a new label: {path}")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        try:
            os.replace(temporary, path)
        except PermissionError:
            if not locked_index_fallback or not path.is_file():
                raise
            # Preserve the existing file's ownership/ACL. The caller's lock
            # serializes publishers; a reader may briefly see an updating index.
            # Completed participant books remain independently immutable.
            with path.open("r+b") as existing:
                existing.write(content)
                existing.truncate()
                existing.flush()
                os.fsync(existing.fileno())
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return path


def json_bytes(payload):
    """Use the same deterministic JSON representation for hashes and publication."""
    return (json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def load_report_manifest(path):
    """Validate the manifest and exact self-contained HTML that was reviewed."""
    path = Path(path)
    payload = json.loads(path.read_text())
    expected = payload.get("manifest_sha256")
    if (payload.get("schema") != REPORT_SCHEMA or payload.get("purpose") != REPORT_PURPOSE
            or digest({k: v for k, v in payload.items() if k != "manifest_sha256"}) != expected):
        raise RuntimeError(f"Invalid component-report manifest: {path}")
    name = payload.get("report_html", "")
    if Path(name).name != name or not name.endswith(".html"):
        raise RuntimeError("Component-report HTML must be adjacent to its manifest.")
    if file_hash(path.parent / name) != payload["report_sha256"]:
        raise RuntimeError("Component-report HTML changed after generation.")
    return payload


# -----------------------------------------------------------------------------
# Diagnostic copies: exact native peak times, visible gaps, unchanged ICA caches
# -----------------------------------------------------------------------------

def pack_trace_view(sources, target_rate=200.):
    """Return compact float32 samples and uint16 native offsets for browser display.

    Each complete/partial bin contributes its largest absolute native sample.
    The bin is NaN if ANY acquired sample is nonfinite; padding is not acquired
    data and cannot invalidate the final partial bin. Offsets retain exact peak
    timing. No gap is interpolated and no original source sample is changed.
    """
    sfreq = float(sources.info["sfreq"])
    if not math.isfinite(target_rate) or target_rate <= 0:
        raise ValueError("Display rate must be finite and positive.")
    step = max(1, int(math.floor(sfreq / target_rate)))
    if step > np.iinfo(np.uint16).max:
        raise ValueError("Display-bin width exceeds the native offset representation.")
    data = sources.get_data()
    n_bins = int(math.ceil(sources.n_times / step))
    values = np.empty((len(data), n_bins), dtype="<f4")
    offsets = np.empty((len(data), n_bins), dtype="<u2")
    scales, unscored = [], []
    for component, trace in enumerate(data):
        finite = np.isfinite(trace)
        scale = float(np.max(np.abs(trace[finite]))) if finite.any() else 1.
        scale = max(scale, np.finfo(float).eps)
        padded = np.zeros(n_bins * step)
        padded[:len(trace)] = trace
        blocks = padded.reshape(n_bins, step)
        magnitudes = np.where(np.isfinite(blocks), np.abs(blocks), -np.inf)
        if len(trace) < len(padded):
            magnitudes[-1, len(trace) % step:] = -np.inf
        positions = np.argmax(magnitudes, axis=1)
        selected = blocks[np.arange(n_bins), positions] / scale
        invalid = ~np.isfinite(blocks).all(axis=1)
        selected[invalid] = np.nan
        values[component] = selected
        offsets[component] = positions
        scales.append(scale)
        unscored.append(int(invalid.sum()))
    return {"encoding": "float32-le-and-uint16-le-base64", "method": "peak-abs-native-v1", "step": step, "n_bins": n_bins,
        "display_rate_hz": sfreq / step, "native_sfreq_hz": sfreq,
        "native_n_times": int(sources.n_times), "native_first_samp": int(sources.first_samp),
        "normalization_scales": scales, "unscored_bins_by_component": unscored,
        "values": base64.b64encode(values.tobytes()).decode("ascii"),
        "offsets": base64.b64encode(offsets.tobytes()).decode("ascii")}


def component_maps(model):
    """Render all fitted components for each available MEG sensor type, without fitting."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    images = {}
    for kind in ("mag", "grad"):
        if not len(mne.pick_types(model.info, meg=kind, exclude=[])):
            continue
        figures = model.plot_components(ch_type=kind, show=False, ncols=5, size=1.5)
        figures = figures if isinstance(figures, list) else [figures]
        if len(figures) != 1:
            for figure in figures:
                plt.close(figure)
            raise RuntimeError("Unexpected multi-page preliminary component maps.")
        try:
            buffer = io.BytesIO()
            figures[0].savefig(buffer, format="png", dpi=140, bbox_inches="tight")
            images[kind] = "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")
        finally:
            plt.close(figures[0])
    if not images:
        raise RuntimeError("Cached preliminary ICA has no MEG sensor maps.")
    return images


def build_run_asset(cfg, record, folder, *, display_rate=200.):
    """Validate once and cache complete report data for one exact acquisition."""
    before = input_fingerprints(cfg, record)
    identity = {"source": record.source_filename, "inputs": before, "display_rate": display_rate,
                "format": REPORT_SCHEMA}
    path = triage_path(cfg, folder / "_cache" / cfg["layout"]["subject_bids_id"] /
                       f"{record.source_stem}-{digest(identity)}.json")
    if path.exists():
        asset = json.loads(path.read_text())
        if asset.get("identity") != identity or digest({k: v for k, v in asset.items() if k != "asset_sha256"}) != asset.get("asset_sha256"):
            raise RuntimeError(f"Cached report asset changed: {path}")
        return asset, path
    print(f"  {record.source_filename}: validating cached ICA", flush=True)
    cache = load_cached_aid(cfg, record)
    print(f"  {record.source_filename}: loading diagnostic traces and rendering maps", flush=True)
    sources = cache["native_sources"]
    run = {"status": "ready", "source_filename": record.source_filename,
        "acquisition_run": record.acquisition_run_label, "binding": component_binding(record, cache),
        "n_components": int(cache["model"].n_components_), "sfreq": float(sources.info["sfreq"]),
        "duration_sec": float(sources.n_times / sources.info["sfreq"]),
        "converged": bool(cache["metadata"].get("converged", False)),
        "maps": component_maps(cache["model"]), "display": pack_trace_view(sources, display_rate),
        "summary": "Pending component review"}
    if input_fingerprints(cfg, record) != before or has_manual_work(cfg, record):
        raise RuntimeError("Inputs/manual work changed during export; rerun this participant.")
    asset = {"identity": identity, "run": run}
    asset["asset_sha256"] = digest(asset)
    atomic_bytes(cfg, path, json_bytes(asset), immutable=True)
    return asset, path


# -----------------------------------------------------------------------------
# One offline participant book, plus a batch index and reusable note templates
# -----------------------------------------------------------------------------

def render_html(payload):
    """Embed data and viewer without external URLs or script-injection boundaries."""
    encoded = json.dumps(payload, separators=(",", ":"), allow_nan=False)
    encoded = encoded.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    return TEMPLATE.read_text().replace("__REPORT_VIEWER__", VIEWER.read_text()).replace("__REPORT_DATA__", encoded)


def decision_rows(manifest):
    """Emit only report-ready runs; blank choices remain unreviewed, never NONE."""
    return [{"participant": manifest["participant"], "run": row["acquisition_run"],
        "source_filename": row["source_filename"], "blink": "", "eye_movement": "", "notes": ""}
        for row in manifest["runs"] if row["status"] == "ready"]


def build_participant(cfg, records, folder, *, resume=False, display_rate=200.):
    """Preserve saved work and create one immutable participant report and manifest."""
    participant = "P" + str(cfg["subject"]["id"]).zfill(2)
    stem = cfg["layout"]["subject_bids_id"]
    manifest_path = folder / f"{stem}_manifest.json"
    require_review(cfg, "03_A_MaxFilt", stage_output_dirs(cfg)["reports"])
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Report exists; use --resume or a new --label: {manifest_path}")
        old = load_report_manifest(manifest_path)
        if old["viewer_sha256"] != file_hash(VIEWER) or old["template_sha256"] != file_hash(TEMPLATE):
            raise RuntimeError("Report viewer changed; use a new label to preserve viewed reports.")
        if old["display_rate"] != display_rate:
            raise RuntimeError("Display rate changed; use a new report label.")
        if {r["source_filename"] for r in old["runs"]} != {r.source_filename for r in records}:
            raise RuntimeError("Report source coverage changed; use a new label.")
        for record in records:
            bound = next((r for r in old["runs"] if r["source_filename"] == record.source_filename), None)
            if bound is None or bound.get("inputs") != input_fingerprints(cfg, record):
                raise RuntimeError("Report inputs changed; use a new label for this participant.")
        print(f"{participant}: existing report verified and reused", flush=True)
        return old
    runs, manifest_runs = [], []
    for record in records:
        inputs = input_fingerprints(cfg, record)
        base = {"source_filename": record.source_filename, "acquisition_run": record.acquisition_run_label}
        if has_manual_work(cfg, record):
            run = {**base, "status": "manual_work_preserved", "summary": "Reopen the existing manual review normally."}
        else:
            chosen = selection_path(cfg, record)
            choice = json.loads(chosen.read_text()) if chosen.exists() else None
            if choice and choice.get("status") == "selected":
                run = {**base, "status": "choices_already_saved",
                    "summary": "Blinks: " + str(choice.get("blink_components")) +
                        "; eyes: " + str(choice.get("eye_movement_components"))}
            else:
                asset, asset_path = build_run_asset(cfg, record, folder, display_rate=display_rate)
                run = asset["run"]
                if asset["identity"]["inputs"] != inputs:
                    raise RuntimeError("Inputs changed before report generation; retry this participant.")
        runs.append(run)
        manifest_runs.append({**{k: v for k, v in run.items() if k not in ("maps", "display")},
                              "inputs": inputs})
    for record, bound in zip(records, manifest_runs):
        if input_fingerprints(cfg, record) != bound["inputs"]:
            raise RuntimeError("Inputs/choices changed while building this participant; retry.")
    identity = {"participant": participant, "runs": manifest_runs,
                "display_rate": display_rate, "viewer_sha256": file_hash(VIEWER),
                "template_sha256": file_hash(TEMPLATE)}
    report_id = digest(identity)
    report = render_html({"participant": participant, "report_id": report_id, "runs": runs})
    report_path = folder / f"{stem}.html"
    atomic_bytes(cfg, report_path, report.encode("utf-8"), immutable=True)
    manifest = {"schema": REPORT_SCHEMA, "purpose": REPORT_PURPOSE, "participant": participant,
        "subject_bids_id": stem, "report_id": report_id, "report_html": report_path.name,
        "report_sha256": hashlib.sha256(report.encode("utf-8")).hexdigest(),
        "generated_utc": utc_now(), **identity}
    manifest["manifest_sha256"] = digest(manifest)
    atomic_bytes(cfg, manifest_path, json_bytes(manifest), immutable=True)
    atomic_bytes(cfg, folder / f"{stem}_decisions_template.csv", csv_content(decision_rows(manifest)).encode("utf-8"), immutable=True)
    print(f"{participant}: report ready: {report_path}", flush=True)
    return manifest


def update_index(cfg, folder):
    """Refresh navigation under a lock, with a shared-storage replacement fallback."""
    folder = triage_path(cfg, folder)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / ".index.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        manifests = [json.loads(p.read_text()) for p in sorted(folder.glob("sub-*_manifest.json"))]
        links, decisions = [], []
        for manifest in manifests:
            pending = len(decision_rows(manifest))
            links.append('<tr><td><a href="' + html.escape(manifest["report_html"], quote=True) + '">' +
                html.escape(manifest["participant"]) + '</a></td><td>' + str(pending) +
                '</td><td><a href="' + html.escape(manifest["subject_bids_id"], quote=True) +
                '_decisions_template.csv">Decision template</a></td></tr>')
            decisions.extend(decision_rows(manifest))
        index = ('<!doctype html><html lang="en"><meta charset="utf-8"><meta http-equiv="refresh" content="30"><title>Ocular component reports</title>' +
            '<style>body{font:16px system-ui;max-width:900px;margin:40px auto;padding:20px}' +
            'td,th{padding:12px;text-align:left;border-bottom:1px solid #ddd}table{width:100%}</style>' +
            '<h1>Ocular component reports</h1><p>One report per participant. Each contains a run selector, ' +
            'component maps and navigable traces. This index reloads every 30 seconds; participant reports do not auto-reload.</p>' +
            '<p>Record choices per run using zero-based indices, NONE or DEFER. Reports do not save decisions. ' +
            'Work on a copy of the template or send your notes in chat.</p>' +
            '<table><tr><th>Participant</th><th>Runs to review</th><th>Notes</th></tr>' +
            ''.join(links) + '</table></html>')
        atomic_bytes(cfg, folder / "index.html", index.encode("utf-8"), locked_index_fallback=True)
        atomic_bytes(cfg, folder / "decisions_template.csv", csv_content(decisions).encode("utf-8"), locked_index_fallback=True)


def run_export(args):
    """Generate selected participants without scheduling jobs or opening GUI windows."""
    label = args.label or "components-" + utc_now().replace(":", "").replace("-", "") + "-" + uuid4().hex[:6]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label) or label in (".", ".."):
        raise ValueError("--label must be a simple folder name.")
    if getattr(args, "index_only", False):
        # A navigation repair reads small completed manifests only. It does not
        # revisit source FIFs, fit caches, saved choices or immutable report data.
        cfg = subject_config(args, args.subjects[0])
        folder = triage_path(cfg, triage_root(cfg) / "component_reports" / label)
        if not folder.is_dir() or not any(folder.glob("sub-*_manifest.json")):
            raise FileNotFoundError("No completed participant manifests in this report batch.")
        update_index(cfg, folder)
        print(f"Refreshed navigation only: {folder / 'index.html'}", flush=True)
        return 0
    errors = []
    for subject in sorted(set(args.subjects)):
        try:
            cfg = subject_config(args, subject)
            folder = triage_path(cfg, triage_root(cfg) / "component_reports" / label)
            print(f"P{subject:02}: preparing one participant report", flush=True)
            records = discover_raw_run_records(cfg)
            build_participant(cfg, records, folder, resume=args.resume, display_rate=args.display_rate)
            update_index(cfg, folder)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            errors.append(f"P{subject:02}: {error}")
            print(errors[-1], flush=True)
    print(f"Finished with {len(errors)} input issue(s). Existing decisions and scientific derivatives were not modified.", flush=True)
    return 1 if errors else 0


def main(argv=None):
    """Run batch report generation; --resume verifies existing complete books/assets."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common_arguments(parser)
    parser.add_argument("--label", help="Shared new batch folder name, also used across array tasks.")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--index-only", action="store_true", help="Refresh the existing batch index only; no report/FIF regeneration.")
    parser.add_argument("--display-rate", type=float, default=200., help="Peak-preserving display bins/second; native data unchanged.")
    args = parser.parse_args(argv)
    if not math.isfinite(args.display_rate) or args.display_rate <= 0:
        parser.error("--display-rate must be finite and positive")
    if (args.resume or args.index_only) and not args.label:
        parser.error("--resume/--index-only requires --label")
    return run_export(args)


if __name__ == "__main__":
    raise SystemExit(main())
