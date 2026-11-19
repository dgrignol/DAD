#!/usr/bin/env python3
"""Review cached ocular ICA components across the DAD cohort, with resumable saves.

Usage (from the DAD root, in the local dRSA environment)
-----------------------------------------------------
    python meg/code/preprocessing_basic/review_ocular_components.py --check-only
    python meg/code/preprocessing_basic/review_ocular_components.py
    python meg/code/preprocessing_basic/review_ocular_components.py --subjects 06 \
      --source sub06_run04_260526.fif
    # --revisit reopens saved choices, but never completed manual interval work.
    # --browser-backend matplotlib selects the alternate source browser.

Inputs: current subject tuning, accepted SSS, report-bound muscle evidence, and
the existing frontal ICA model/native-source caches. No ICA is fitted locally.
The viewer shows component maps plus scrollable, display-normalized source
traces. Zero-based indices/ICA names may be selected independently for blinks
and eye movements; the same component may serve both. NONE disables one detector
without establishing absence of artifacts. DEFER and QUIT never become empty
selections. Every SAVE is bound to the exact source/cache/model and persisted
atomically under derivatives/annotation_triage/component_selections, alongside
an exportable per-subject CSV. These records are triage choices, not production
tuning, whole-source interval reviews or ICA component-removal decisions.

Data flow: read-only discovery/validation -> skip existing manual work and valid
saved choices -> maps/source inspection -> explicit choices -> guarded sidecar
save -> next recording. Existing manual files/receipts, including P04/P05 and
P06 runs 1-3, always win. New manual saves are checked again immediately before
component choices are persisted. P06 runs 4-5 are ordinary unreviewed inputs.
--check-only reads headers/models and prints availability without any writes.
Missing/stale caches produce targeted preparation advice, never a silent refit.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import tempfile

import mne
import numpy as np

from annotation_triage import (SCHEMA, common_arguments, component_binding, csv_content,
    has_manual_work, load_context, selection_path, subject_config, triage_path, validate_selection)
from blink_annotation import load_cached_aid
from continuous_cleaning import require_review
from preprocessing_stage_utils import (discover_raw_run_records, guarded_path,
    stage_output_dirs, utc_now, write_text)


# -----------------------------------------------------------------------------
# Explicit user choices and safe resumable persistence
# -----------------------------------------------------------------------------

def parse_components(text, count):
    """Parse zero-based components or NONE; DEFER/QUIT are control actions."""
    text = text.strip().upper()
    if text in ("DEFER", "QUIT"):
        return text
    if text == "NONE":
        return []
    tokens = [re.fullmatch(r"(?:ICA)?(\d+)", part.strip()) for part in text.split(",")]
    if not all(tokens):
        raise ValueError("Use 0 / ICA000, comma-separated components, NONE, DEFER or QUIT.")
    values = sorted({int(token.group(1)) for token in tokens})
    if not all(value < count for value in values):
        raise ValueError("A component is outside the displayed zero-based range.")
    return values


def ask_components(kind, count):
    """Retry invalid answers without accepting a blank as reviewed zero."""
    while True:
        try:
            return parse_components(input(f"{kind} components (0/ICA000, NONE, DEFER, QUIT): "), count)
        except ValueError as error:
            print(error)


def save_selection(cfg, record, payload, *, previous=None):
    """Atomically save only triage choices; reject another writer's intervening edit."""
    path = guarded_path(selection_path(cfg, record), cfg["resolved_paths"]["data_root"])
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = triage_path(cfg, path.with_name(f".{path.name}.lock"))
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        current = json.loads(path.read_text()) if path.exists() else None
        if current != previous:
            raise RuntimeError("Component choices changed in another session; reopen this recording.")
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
                temporary = Path(handle.name)
                json.dump(payload, handle, indent=2, sort_keys=True)
                handle.write("\n")
            os.replace(temporary, path)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()
    # JSONs are authoritative; the CSV is a disposable convenience export.
    rows = []
    for saved in sorted(path.parent.glob("*.json")):
        row = json.loads(saved.read_text())
        rows.append({"source_filename": row["binding"]["source_filename"],
            "acquisition_run": row["binding"]["acquisition_run"], "status": row["status"],
            "blink_components": ",".join(map(str, row.get("blink_components", []))) if row["status"] == "selected" else "",
            "eye_movement_components": ",".join(map(str, row.get("eye_movement_components", []))) if row["status"] == "selected" else "",
            "model_sha256": row["binding"]["model_sha256"], "observations": row.get("observations", "")})
    write_text(triage_path(cfg, path.parent / "component_selections.csv"), csv_content(rows), raw_root=cfg["resolved_paths"]["data_root"])
    return path


# -----------------------------------------------------------------------------
# Component-only display; no annotation browser, EDF or final ICA application
# -----------------------------------------------------------------------------

def show_components(cache, title, backend):
    """Display maps and native cached sources, preserving amplitudes used by detection."""
    import matplotlib.pyplot as plt
    mne.viz.set_browser_backend(backend)
    figures = cache["model"].plot_components(show=False)
    figures = figures if isinstance(figures, list) else [figures]
    for figure in figures:
        figure.suptitle(title + " | preliminary ocular aid (0 = ICA000)")
        figure.show()
    view = cache["native_sources"].copy().load_data()
    def normalize(values):
        finite = np.isfinite(values)
        scale = float(np.max(np.abs(values[finite]))) if finite.any() else 1.
        return values / max(scale, np.finfo(float).eps)
    view.apply_function(normalize, picks="all", verbose="ERROR")
    try:
        view.plot(duration=20, n_channels=len(view.ch_names), scalings={"misc": 1.},
                  title=title + " | close source browser to enter choices", block=True)
    finally:
        for figure in figures:
            plt.close(figure)


def run_review(args):
    """Walk exact sources deterministically, resuming only compatible selections."""
    errors = []
    matched = 0
    for subject in sorted(set(args.subjects)):
        try:
            print(f"P{subject:02}: validating participant inputs", flush=True)
            cfg = subject_config(args, subject)
            require_review(cfg, "03_A_MaxFilt", stage_output_dirs(cfg)["reports"])
            records = discover_raw_run_records(cfg)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            errors.append(f"P{subject:02}: {error}")
            print(errors[-1], flush=True)
            continue
        for record in records:
            if args.source and record.source_filename != args.source:
                continue
            matched += 1
            title = f"P{subject:02} run {record.acquisition_run_label:02}: {record.source_filename}"
            try:
                cfg = subject_config(args, subject)
                if has_manual_work(cfg, record):
                    load_context(cfg, record)  # Validate, then preserve; never reseed manual decisions.
                    print(f"{title}: preserved manual review", flush=True)
                    continue
                cache = load_cached_aid(cfg, record)
                path = selection_path(cfg, record)
                previous = json.loads(path.read_text()) if path.exists() else None
                if previous and previous.get("status") == "selected":
                    if not args.revisit:
                        validate_selection(previous, record, cache)
                        print(f"{title}: compatible choices saved; skipped", flush=True)
                        continue
                if args.check_only:
                    print(f"{title}: cache valid; {cache['model'].n_components_} components; choices needed", flush=True)
                    continue
                print(f"\n{title}\n0 = ICA000 (first component). NONE disables a detector; gaps remain unscored.", flush=True)
                print("Cache convergence:", cache["metadata"].get("converged"), "Previous choices:",
                      {key: previous.get(key) for key in ("status", "blink_components", "eye_movement_components")} if previous else None)
                while True:
                    show_components(cache, title, args.browser_backend)
                    blink = ask_components("Blink", cache["model"].n_components_)
                    eye = ask_components("Eye movement", cache["model"].n_components_) if isinstance(blink, list) else blink
                    if blink == "QUIT" or eye == "QUIT":
                        print("Stopped; earlier explicit saves are retained.")
                        return 1 if errors else 0
                    deferred = blink == "DEFER" or eye == "DEFER"
                    if not deferred:
                        print("Blink:", blink, "Eye movement:", eye)
                        action = input("SAVE choices, VIEW again, DEFER recording or QUIT: ").strip().upper()
                        if action == "VIEW":
                            continue
                        if action == "QUIT":
                            return 1 if errors else 0
                        if action not in ("SAVE", "DEFER"):
                            print("No changes saved; inspect/enter choices again.")
                            continue
                        deferred = action == "DEFER"
                    observations = input("Optional observations: ").strip()
                    latest = subject_config(args, subject)
                    if has_manual_work(latest, record):
                        raise RuntimeError("Manual work was saved during component inspection; it takes precedence.")
                    current_cache = load_cached_aid(latest, record)
                    if component_binding(record, current_cache) != component_binding(record, cache):
                        raise RuntimeError("Cached model changed during inspection; reopen this recording.")
                    payload = {"schema": SCHEMA, "purpose": "ocular_component_selection_only",
                        "status": "deferred" if deferred else "selected", "binding": component_binding(record, cache),
                        "reviewer": args.reviewer, "saved_utc": utc_now(), "observations": observations}
                    if not deferred:
                        payload.update(blink_components=blink, eye_movement_components=eye)
                    print("Saved:", save_selection(latest, record, payload, previous=previous), flush=True)
                    break
            except (OSError, ValueError, RuntimeError, KeyError) as error:
                errors.append(f"{title}: {error}")
                print(errors[-1], flush=True)
                print("Reconcile the stated input issue first. For a stale/missing ICA cache, the targeted cluster command is:")
                print(f"DAD_BLINK_AID_ONLY=1 DAD_REBUILD=1 DAD_SOURCE={record.source_filename} sbatch --array={subject} meg/cluster/submit_preprocessing_review.sh")
    if args.source and not matched:
        raise ValueError("--source did not match a configured recording in the selected subjects.")
    print(f"Finished. {len(errors)} input issue(s); no production files or scientific review receipts changed.")
    return 1 if errors else 0


def main(argv=None):
    """Run the local component-only viewer, or a strictly read-only cache preflight."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common_arguments(parser)
    parser.add_argument("--source", help="Optional exact source filename.")
    parser.add_argument("--check-only", action="store_true", help="Read-only validation; no GUI or saved choices.")
    parser.add_argument("--revisit", action="store_true", help="Reopen compatible choices; manual work still takes precedence.")
    parser.add_argument("--reviewer", default="User")
    parser.add_argument("--browser-backend", choices=("qt", "matplotlib"), default="qt")
    return run_review(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
