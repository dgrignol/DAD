#!/usr/bin/env python3
"""Review ocular, muscle and other BAD intervals for one Stage 4 acquisition.

Usage
-----
    python meg/code/preprocessing_basic/manual_eye_annotations.py \
      --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
      --tuning-config /path/to/03_tuning_choices.yaml \
      --source sub04_run01_220526.fif --basis misc_and_frontal

    # Add --preliminary-ica for supplementary component traces.
    # Recover a failed annotation save with --overwrite --recover-annotations.
    # Add --browser-backend matplotlib for its thin default trace strokes.
    # Add --onset-linewidth 3 for more prominent corrected onset lines (default 2).
    # Add --onset-color "#0057B8" for blue markers; default is opaque black.
    # Add --high-pass-view 0.01 to filter only the browser's frontal MEG copy.
    # Add --no-timing-marks to hide corrected onset/epoch guidance entirely.
    # Targeted triage review: add --triage-proposals /path/to/exact-mask.json.
    # This loads the exact provisional intervals counted in the retention report,
    # disables EDF guidance for this session and requires the usual full SAVE.
    # Existing manual work is reopened normally, never replaced by triage masks.

The annotation selector always offers BAD_blink (manual blinks),
BAD_eye_movement and BAD_segment (unclear/non-prototypical artifacts). These
descriptions are registered in the browser UI only, without dummy intervals.
Unused descriptions are absent from the saved annotation file; every browser
opening offers them again, including after BACK/REDETECT. BAD_blink_auto remains
separate, so manually added blinks can keep their own description.

BAD_squid_jumps belongs to Stage 4A and is protected SSS context. It is displayed
but excluded from editable eye/muscle snapshots; changes must return through
review_squid_jumps.py so SSS and its acceptance are reconciled. Current source
SQUID labels are preserved once during 4B rebuilding. Existing saved manual
intervals/receipts remain intact; after changed SSS they are starting points for
explicit revalidation rather than automatically accepted replacement reviews.

By default, the browser adds corrected motion-onset event lines and green
EPOCH_preview annotations for the accepted Stage 3 window (-0.5 to +2.667 s).
The existing Stage 2/3 loaders verify their accepted hashes and eligibility;
eligible trials use green guides. Gaze-break trials use orange
EPOCH_invalid_preview guides around their stored corrected onset, showing the
intended window of an already excluded trial. Native code-150 markers are orange,
opaque and above traces/annotations, sharing --onset-linewidth. They receive no
visual timing shift and create no BAD interval or eligibility decision. Missing
usable timing or ambiguous marker ownership is reported, never guessed.
Eligible trials can still be rejected later for incomplete epochs/artifacts. The
default review is visual_timing/s3-v002/reviews/sub-XX_review.json below the
derivative root; --timing-review or epoching.timing_review_path can select an
explicit review instead. No new lag is estimated or applied to continuous data.
Press 'a', then Qt 'Select Visible' (or the Matplotlib label checkboxes) to
toggle either preview description. Both are temporary display guides: edits are
discarded and accepted coordinates restored whenever the browser reopens.
They are removed before saving artifacts, and never affect epoch rejection.

Optional --preliminary-ica LOADS a source-bound cluster cache prepared by
blink_annotation.py. It never fits ICA locally. First inspect maps and cached
source traces, then select BLINK components followed by EYE-MOVEMENT components.
Each prompt accepts comma-separated zero-based numbers (0 = ICA000, first) or
ICA names, or explicit NONE. Blank/invalid selections retry. The same component
may serve both groups. Peaks from blink components propose BAD_blink_auto;
sustained median-smoothed excursions from eye components independently propose
BAD_eye_movement_auto. See blink_annotation.py for the duration method, native
sample/NaN-gap handling and the 1-40 Hz cache's long-plateau limitation.
NONE skips only its detector, without claiming the source is artifact-free.
Missing eye fields in older receipts mean NONE in memory, not a corrupt receipt
or a completed whole-source review. Loading never migrates saved files.

The combined disposable browser includes selected AID traces and separate
BLINK_score / EYE_MOVEMENT_score channels and overviews (1 = threshold, capped
at 2 only in the browser). Automatic BAD intervals remain provisional until
SAVE. BACK reopens selection on the same fit, keeping all edits; unchanged
groups are not regenerated. Changed nonempty groups ask before regeneration;
NONE removes only untouched session proposals for that group. Declining keeps
the prior score/proposal provenance distinct from the pending selection.
REDETECT asks BLINK or EYE, then changes only that detector's threshold/results.
Manual edits, deletions, saved intervals and each group's vetoes are retained.
MISC usefulness and
observations are asked AFTER inspection. SAVE confirms whole-source review;
finish those prompts and wait for "Saved reviewed ocular annotations" to confirm
both the FIF and tuning receipt were written. Different runs may be reviewed in
overlapping processes: SAVE merges only this source's receipt under a short
write lock, retaining reviews saved since the browser opened. Avoid two editors
of the same source; annotation edits are not merged across processes.
DISCARDALL exits without saving this session. Blank/typo actions repeat the
prompt. This never removes components or supplies a final model. Failure to
find a component is not evidence of no artifacts.

Inputs: accepted 4A SSS FIF and available MISC/frontal MEG channels. Outputs:
one source-specific -annot.fif plus review metadata in the existing tuning YAML.
SAVE retains ALL BAD categories, including edited BAD_muscle intervals. When
current provisional 4B muscle evidence was loaded, SAVE records its reviewed
muscle set as a replacement: deletions (including deleting every interval),
resizing and additions survive reopening and the later 4B rebuild. Legacy
ocular-only receipts never imply a reviewed-empty muscle set. Changed muscle
settings invalidate replacement reuse and require explicit reconciliation.
Muscle changes can invalidate the supplementary ICA fit mask; rebuild only the
affected aid cache before reusing it. Final 4B acceptance is still separate.
The operator must inspect the entire recording and explicitly confirm saving,
including when zero annotations were marked. MISC meanings/usability must be
checked; no EOG or EDF synchronization study is performed. Current 4B candidate
sidecars seed editable, non-BAD blink intervals. Optional EDF gaze/pupil traces
are displayed by default on each trial's marker map with NaN gaps and coverage.
These normalized EDF MISC channels are temporary; original MISC/MEG channels are intact.
Use --no-edf-traces for candidates alone. Explicit SAVE confirms complete-source
review and converts remaining (including long) candidates to BAD_blink; delete
or adjust uncertain candidates before saving. Missing EDF guidance is not proof
of no artifacts. Unchanged reviewed proposals are not seeded again on reopening.
If a previous save wrote annotations but failed before its receipt, explicitly
use --overwrite --recover-annotations. The existing file becomes an unreviewed
starting point, without reseeding deleted candidates; recheck the current EDF
guidance and the complete source before SAVE. DISCARDALL leaves existing files
and receipts unchanged; it discards only the current session's edits.

Optional --high-pass-view HZ applies a high-pass filter to the selected frontal
MEG channels in the disposable browser copy only. MISC channels remain
unfiltered. It does not modify the source FIF, annotation coordinates, saved
annotations or cached ICA source traces; auxiliary EDF and score traces retain
their existing display values.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from contextlib import nullcontext
import inspect
import json
from pathlib import Path
import re

import mne
import numpy as np

from continuous_cleaning import (annotation_rows, coordinates, digest, file_hash,
    file_stamp, muscle_parameters, require_review, reviewed_muscle_override)
from blink_annotation import (AUTO_LABEL, EYE_AUTO_LABEL, EYE_METHOD, annotation_key,
    blink_settings, detect_blinks, detect_eye_movements, eye_movement_settings,
    frontal_channel_groups, load_cached_aid, remember_deletions, replace_proposals)
from edf_blink_candidates import (CANDIDATE_LABELS, load_candidates, add_edf_display,
                                  promote_candidates_after_review)
from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import (discover_raw_run_records, derivative_path_for_run,
    stage_output_dirs, prepare_output_path, update_tuning_file, utc_now)
from downstream_stage_utils import csv_bool, load_trigger_event_rows
from trial_preparation import evidence_paths, load_prepared_decisions, read_csv
from visual_timing_utils import analysis_eligible, load_reviewed_timing, round_half_away


EYE_LABELS = ("BAD_blink", "BAD_eye_movement")
REVIEW_DESCRIPTIONS = EYE_LABELS + ("BAD_segment",)
EPOCH_LABEL = "EPOCH_preview"
EPOCH_COLOR = "#238b45"
INVALID_EPOCH_LABEL = "EPOCH_invalid_preview"
GAZE_BREAK_CODE = 150
GAZE_BREAK_COLOR = "#ff8c00"
PREVIEW_COLORS = {EPOCH_LABEL: EPOCH_COLOR, INVALID_EPOCH_LABEL: GAZE_BREAK_COLOR}
DEFAULT_ONSET_COLOR = "#000000"


# -----------------------------------------------------------------------------
# Exact-source annotation storage; an empty reviewed set is still a real file
# -----------------------------------------------------------------------------

def manual_path(cfg, record):
    """Resolve the annotation data path without introducing another config file."""
    return Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]) / "manual_eye" / f"{record.source_stem}-annot.fif"


def eye_annotations(raw):
    """Select the reviewed BAD snapshot, including muscle, in acquisition time.

    The historical function/file names are retained for compatibility. A separate
    receipt field distinguishes a full muscle review from legacy eye-only files.
    Non-BAD guides such as EPOCH_preview and EDF_no_guidance are never saved.
    Stage 4A-owned BAD_squid_jumps also remains source context, not an editable
    4B snapshot category; this avoids duplicate or obsolete copied intervals.
    """
    from squid_jumps import without_owned_context
    editable = without_owned_context(raw.annotations)
    keep = [i for i, label in enumerate(editable.description)
            if label.lower().startswith("bad")]
    return editable[keep].copy()


def load_manual(cfg, record, raw, input_path):
    """Require a reviewed file with matching source, coordinates and SSS input."""
    receipt = cfg.get("artifact_annotation", {}).get("manual_eye_reviews", {}).get(record.source_filename)
    path = manual_path(cfg, record)
    if not receipt or receipt.get("status") != "reviewed" or not path.is_file():
        raise RuntimeError(f"Ocular review missing for {record.source_filename}; reviewed zero is explicit.")
    if receipt.get("coordinates") != coordinates(raw) or receipt.get("sss_input") != file_stamp(input_path):
        raise RuntimeError(f"Ocular review source/SSS coordinates changed: {record.source_filename}")
    if receipt.get("annotation_sha256") != file_hash(path):
        raise RuntimeError(f"Ocular annotation file changed since review: {path}")
    ann = mne.read_annotations(path)
    if len(ann) != receipt.get("annotation_count"):
        raise RuntimeError("Ocular review count does not match the saved annotations.")
    muscle_review = receipt.get("muscle_review")
    if muscle_review and sum(label.lower() == "bad_muscle" for label in ann.description) != muscle_review.get("annotation_count"):
        raise RuntimeError("Reviewed muscle count does not match the saved annotations.")
    from squid_jumps import without_owned_context
    return without_owned_context(ann), receipt


def load_muscle_context(cfg, record, raw, input_path, dirs):
    """Load only current, report-bound provisional muscle annotations for QC.

    Returns (annotations, basis) or (None, None) when no provisional evidence
    exists. Stale evidence fails before signal loading. The existing annotated
    FIF signature binds every detector setting and its used manual receipt;
    no new report/cache system is introduced. The basis becomes a reviewed
    replacement only on explicit SAVE, even if zero muscle intervals remain.
    """
    import json
    report_path = dirs["reports"] / "03_B_ArtefactAnnotation_review_report.json"
    if not report_path.is_file():
        return None, None
    report = json.loads(report_path.read_text())
    rows = [row for row in report["runs"] if row["source_filename"] == record.source_filename]
    candidate = derivative_path_for_run(record, dirs["sss"], "_ann.fif")
    if len(rows) != 1 or not candidate.is_file() or rows[0].get("input_stamp") != file_stamp(input_path):
        raise RuntimeError("Missing/stale source-specific 4B muscle evidence; refresh the affected report before review.")
    settings = cfg["artifact_annotation"]
    by_run = settings.get("muscle_threshold_by_run", {})
    run = record.acquisition_run_label
    threshold = next((float(by_run[key]) for key in (run, str(run), f"{run:02d}", f"run{run:02d}")
                      if key in by_run), float(settings["muscle_threshold"]))
    parameters = muscle_parameters(settings, threshold)
    evidence = next((entry for entry in report.get("evidence_files", [])
                     if Path(entry["relative_path"]).name == candidate.name), None)
    if not evidence or evidence["stamp"] != file_stamp(candidate):
        raise RuntimeError("Provisional muscle annotation evidence changed; refresh its report before review.")
    header = mne.io.read_raw_fif(candidate, preload=False, verbose="ERROR")
    used_review = rows[0].get("manual_eye_review", {})
    used_review = used_review if used_review.get("status") == "reviewed" else None
    signature = digest({"input": file_stamp(input_path), **parameters, "eye_review": used_review})
    if (float(rows[0]["muscle_threshold"]) != threshold or
            f"stage4_input_sha256={signature}" not in str(header.info.get("description"))):
        raise RuntimeError("Muscle settings changed; reconcile existing muscle reviews and refresh affected 4B evidence first.")
    if coordinates(header) != coordinates(raw) or header.annotations.orig_time != raw.annotations.orig_time:
        raise RuntimeError("Provisional muscle annotations have different recording coordinates.")
    muscle = header.annotations[[i for i, label in enumerate(header.annotations.description)
                                if label.lower() == "bad_muscle"]]
    return muscle, {"threshold": threshold, "settings_sha256": digest(parameters), "basis_stage4_input_sha256": signature,
        "basis_annotated_input": file_stamp(candidate), "basis_muscle_annotation_count": len(muscle),
        "basis_annotation_sha256": digest(annotation_rows(header))}


# -----------------------------------------------------------------------------
# Accepted timing guides: native event samples and disposable epoch annotations
# -----------------------------------------------------------------------------

def gaze_break_guidance(cfg, record, rows, gate):
    """Match native code-150 markers to exact ledger trials and stored timing.

    The Stage 2/3 loaders have already validated source evidence and acceptance.
    Excluded trials have no epoch-event field in that consumer interface, so
    re-read their existing timing-map entries with the same review/manifest/map
    hashes. This does not estimate timing or change eligibility. Ambiguous or
    unmatched markers remain visible and are reported, without guessing an owner.
    """
    events = load_trigger_event_rows(evidence_paths(cfg)["events"])
    markers = [event for event in events if event["source_fif"] == record.source_filename
               and int(event["trigger_code"]) == GAZE_BREAK_CODE]
    owners, issues = {}, []
    for marker in markers:
        key = (str(GAZE_BREAK_CODE), str(marker["event_index"]), str(marker["raw_sample_frame"]))
        candidates = []
        for row in rows:
            # Match a single accepted evidence tuple, including event identity,
            # rather than assigning to the nearest onset across replay attempts.
            evidence = zip(row.get("evidence_trigger_codes", "").split(","),
                           row.get("evidence_event_indices", "").split(","),
                           row.get("evidence_samples_raw", "").split(","), strict=True)
            if (key in evidence and int(row["raw_sample_frame"]) <= int(marker["raw_sample_frame"])
                    < int(row["trial_interval_end_raw_exclusive"])):
                candidates.append(row)
        if csv_bool(marker.get("unresolved_transition")) or len(candidates) != 1:
            issues.append({"sample_raw": int(marker["raw_sample_frame"]), "event_index": marker["event_index"],
                           "reason": "unresolved_or_ambiguous_trial_match", "candidate_count": len(candidates)})
            continue
        owner = candidates[0]
        if not csv_bool(owner.get("gaze_break_evidence")) or analysis_eligible(owner):
            raise RuntimeError("Code 150 contradicts accepted gaze-break eligibility; reconcile the source trial before display.")
        owners[owner["trial_key"]] = owner
    if not owners:
        return {"markers": markers, "onsets": {}, "issues": issues}

    review_path = Path(gate["review_path"])
    if file_hash(review_path) != gate["review_sha256"]:
        raise RuntimeError("Accepted timing review changed while preparing gaze-break guidance.")
    review = json.loads(review_path.read_text())
    manifest_path = Path(review["bundle_manifest"])
    if not manifest_path.is_absolute():
        manifest_path = review_path.parent / manifest_path
    if file_hash(manifest_path) != gate["bundle_manifest_sha256"]:
        raise RuntimeError("Accepted timing manifest changed while preparing gaze-break guidance.")
    manifest = json.loads(manifest_path.read_text())
    map_path = manifest_path.parent / "trial_timing_map.csv"
    if file_hash(map_path) != manifest["output_sha256"][map_path.name]:
        raise RuntimeError("Accepted timing map changed while preparing gaze-break guidance.")
    timing = {entry["trial_key"]: entry for entry in read_csv(map_path)}
    onsets = {}
    for trial_key, owner in owners.items():
        entry = timing[trial_key]
        basis = entry.get("assignment_basis")
        if (entry.get("assignment_status") != "proposed" or not entry.get("proposed_sample_raw")
                or basis not in review["accepted_assignment_bases"]):
            issues.append({"trial_key": trial_key, "reason": "no_usable_reviewed_timing_for_invalid_preview"})
            continue
        protected = ("source_fif", "raw_sample_frame", "condition_code", "sfreq_hz",
                     "recording_first_sample_raw", "recording_end_sample_raw_exclusive")
        if any(str(entry[name]) != str(owner[name]) for name in protected):
            raise RuntimeError("Stored gaze-break timing changed an original trial coordinate.")
        sample, shift = int(entry["proposed_sample_raw"]), int(entry["applied_shift_samples"])
        if sample != int(owner["raw_sample_frame"]) + shift:
            raise RuntimeError("Stored gaze-break timing has inconsistent signed-shift arithmetic.")
        if basis == "digital_fallback":
            if shift or entry.get("lag_estimate_ms"):
                raise RuntimeError("Gaze-break digital fallback incorrectly implies a measured lag.")
        elif shift != round_half_away(float(entry["lag_estimate_ms"]) * float(owner["sfreq_hz"]) / 1000.):
            raise RuntimeError("Stored gaze-break timing has inconsistent native-sample rounding.")
        onsets[trial_key] = sample
    return {"markers": markers, "onsets": onsets, "issues": issues}


def load_timing_preview(cfg, record, raw, review_path=None):
    """Validate accepted timing and return source-specific display guides.

    Reuse the production read-only Stage 2/3 validators, without requiring 4C
    or creating epochs. Original eligibility/anchors remain unchanged. Rows
    without usable stored timing are not silently given digital timing. Gaze-
    break trials receive a separate invalid guide using their stored assignment;
    code 150 itself remains at its unshifted native sample.
    Annotation onsets use acquisition seconds, including first_samp; event
    arrays retain absolute native samples. Truncated windows are shown clipped
    at the recording boundary, reported as incomplete, and never join runs.
    """
    epoch = cfg["epoching"]
    review_path = (review_path or epoch.get("timing_review_path") or
        Path(cfg["resolved_paths"]["derivatives_root"]) / "visual_timing" /
        "s3-v002" / "reviews" / f'{cfg["layout"]["subject_bids_id"]}_review.json')
    prepared = load_prepared_decisions(cfg)
    rows, gate = load_reviewed_timing(cfg, prepared, review_path)
    rows = [row for row in rows if row["source_fif"] == record.source_filename]
    if not rows:
        raise RuntimeError(f"No accepted timing rows for {record.source_filename}; use --no-timing-marks to review without guides.")
    gaze = gaze_break_guidance(cfg, record, rows, gate)
    sfreq = float(raw.info["sfreq"])
    first, end = int(raw.first_samp), int(raw.first_samp + raw.n_times)
    tmin, tmax = float(epoch["tmin_sec"]), float(epoch["tmax_sec"])
    events, onsets, durations, labels = [], [], [], []
    incomplete = valid_count = invalid_count = 0
    for row in rows:
        if (float(row["sfreq_hz"]) != sfreq or int(row["recording_first_sample_raw"]) != first
                or int(row["recording_end_sample_raw_exclusive"]) != end):
            raise RuntimeError("Accepted timing and current continuous recording coordinates differ.")
        if analysis_eligible(row):
            sample, label = int(row["epoch_event_sample_raw"]), EPOCH_LABEL
            valid_count += 1
        elif row.get("trial_key") in gaze["onsets"]:
            sample, label = gaze["onsets"][row["trial_key"]], INVALID_EPOCH_LABEL
            invalid_count += 1
        else:
            continue
        if not first <= sample < end:
            raise RuntimeError("Corrected onset leaves its exact continuous recording.")
        events.append([sample, 0, int(row["condition_code"])])
        # Match the reviewed window's inclusive sample endpoints. Clipping is
        # visual only; later incomplete-epoch accounting still uses the full span.
        left, right = sample + round(tmin * sfreq), sample + round(tmax * sfreq)
        incomplete += int(left < first or right >= end)
        onsets.append(max(first, left) / sfreq)
        durations.append((min(end - 1, right) - max(first, left)) / sfreq)
        labels.append(label)
    for marker in gaze["markers"]:
        sample = int(marker["raw_sample_frame"])
        if not first <= sample < end:
            raise RuntimeError("Code 150 leaves its exact continuous recording.")
        events.append([sample, 0, GAZE_BREAK_CODE])
    events.sort(key=lambda event: (event[0], event[2]))
    annotations = mne.Annotations(onsets, durations, labels,
                                  orig_time=raw.annotations.orig_time)
    metadata = {"enabled": True, **gate, "source_filename": record.source_filename,
        "window_sec": [tmin, tmax], "corrected_onset_count": valid_count,
        "invalid_preview_count": invalid_count, "gaze_break_marker_count": len(gaze["markers"]),
        "gaze_break_review_issues": gaze["issues"],
        "without_accepted_onset_count": len(rows) - valid_count,
        "incomplete_window_count": incomplete, "purpose": "temporary_display_only"}
    print(f"Accepted {gate['version']} timing: {valid_count} valid and {invalid_count} gaze-break invalid previews; "
          f"{len(gaze['markers'])} native code-150 markers; {incomplete} incomplete windows clipped for display.")
    if gaze["issues"]:
        print("Gaze-break guidance needs inspection:", gaze["issues"])
    return {"events": np.asarray(events, dtype=int).reshape(-1, 3),
            "annotations": annotations, "metadata": metadata}


def epoch_color_context(annotations):
    """Use a scoped public palette for green/orange guides on older MNE.

    MNE 1.10/1.11 assign non-BAD label colors in sorted order from this palette;
    BAD/EDGE labels remain red. Supplying the palette only while the browser is
    open preserves its built-in visibility/editing controls and avoids private
    browser hooks or a dependency upgrade. Newer MNE uses annotation_colors
    directly instead (see plot_review_browser).
    """
    import matplotlib as mpl
    from matplotlib.colors import to_rgb
    labels = sorted({str(label) for label in annotations.description
                     if not label.lower().startswith(("bad", "edge"))})
    # MNE removes red shades before assigning non-BAD colors. Filter them here
    # too so the preview's palette position stays exact on both old backends.
    reds = [to_rgb(color) for color in ("#fa8174", "#d62728", "#ff0000")]
    palette = [color for color in mpl.rcParams["axes.prop_cycle"].by_key()["color"]
               if not any(np.allclose(to_rgb(color), red, atol=1.5 / 255.) for red in reds)]
    palette = palette or ["#1f77b4", "#ff7f0e", "#9467bd"]
    colors = [PREVIEW_COLORS.get(label, palette[i % len(palette)])
              for i, label in enumerate(labels)]
    return mpl.rc_context({"axes.prop_cycle": mpl.cycler(color=colors)})


# -----------------------------------------------------------------------------
# Supplementary frontal ICA: compute sources for display, never apply removal
# -----------------------------------------------------------------------------

def choose_aid_components(n_components, kind="blink"):
    """Return one group's sorted zero-based indices; require explicit NONE.

    Names and numbers may be mixed/deduplicated. Blank, negative, out-of-range
    and mixed NONE/index answers retry without changing the selection.
    """
    print("Component indices start at 0: ICA000 = 0 (first), ICA001 = 1 (second).")
    print("Available:", ", ".join(f"ICA{i:03} [{i}]" for i in range(n_components)))
    while True:
        text = input(f"Clear {kind} components (e.g. 0 or ICA000; comma separated), or NONE: ").strip()
        if text == "NONE":
            print(f"{kind}: detector disabled; this does not establish absence of artifacts.")
            return []
        tokens = [re.fullmatch(r"(?:ICA)?(\d+)", value.strip()) for value in text.split(",")]
        if all(tokens):
            components = sorted({int(token.group(1)) for token in tokens})
            if all(i < n_components for i in components):
                print("Selected:", ", ".join(f"ICA{i:03} -> AID{i:03}" for i in components))
                return components
        print("Invalid selection. Use displayed zero-based numbers/ICA names, or exactly NONE; no selection changed.")


def add_preliminary_sources(display, cfg, *, state=None):
    """Review a loaded frontal model and replace only temporary source traces.

    ``state`` holds a validated model/native-grid source cache. No fitting occurs
    here. Per-trace display normalization keeps AID traces visible beside MISC;
    detection always uses the original cached source amplitudes. BACK preserves
    annotations and the same component numbering. It never applies ICA removal.
    """
    if not state or "model" not in state:
        raise RuntimeError("Load the prepared blink-aid cache before component review.")
    model = state["model"]
    metadata = state["metadata"]
    if not metadata["converged"]:
        print("WARNING: supplementary ICA did not converge; inspect carefully and use NONE for any unclear group.")
    model.plot_components(show=True)
    # Rare, strong blink peaks must not turn into apparent missing-data gaps
    # through browser clipping. Normalize only a temporary view, never detection.
    source_view = state["native_sources"].copy().load_data()
    def normalize_source(trace):
        """Keep every finite source peak visible; ICA units are arbitrary."""
        finite = np.isfinite(trace)
        scale = float(np.max(np.abs(trace[finite]))) if finite.any() else 1.
        return trace / max(scale, np.finfo(float).eps)
    source_view.apply_function(normalize_source, picks="all", verbose="ERROR")
    source_view.plot(duration=20, n_channels=int(model.n_components_), scalings={"misc": 1.}, block=True)
    state["components"] = choose_aid_components(int(model.n_components_), "blink")
    # Blink uses the historical state/receipt keys. Eye state is independent;
    # missing selection in a legacy receipt means NONE only in this session.
    eye = state.setdefault("eye_movement", {})
    eye["components"] = choose_aid_components(int(model.n_components_), "eye-movement")
    components = sorted(set(state["components"] + eye["components"]))
    previous = [ch for ch in display.ch_names if ch.startswith("AID")]
    if previous:
        display.drop_channels(previous)
    if components:
        traces = source_view.copy().pick([f"ICA{i:03}" for i in components])
        traces.apply_function(lambda trace: trace * state.get("display_gain", 1.), picks="all", verbose="ERROR")
        traces.rename_channels({name: f"AID{i:03}" for name, i in zip(traces.ch_names, components)})
        traces.set_channel_types({name: "misc" for name in traces.ch_names})
        display.load_data().add_channels([traces], force_update_info=True)
    return {"purpose": "supplementary_annotation_only", "frontal_channels": metadata["identity"]["channels"],
            "cache_fingerprint": metadata["fingerprint"], "model_sha256": metadata["model_sha256"],
            "sensor_selection": "all good MNE Left-frontal/Right-frontal magnetometers and gradiometers",
            "selected_components": list(state["components"]),
            "selected_eye_movement_components": list(eye["components"]),
            "n_components_fitted": int(model.n_components_),
            "fit_duration_sec": float(metadata["fit_summary"]["fit_duration_sec"]),
            "converged": metadata["converged"], "source_display": "per-component maximum-absolute normalization, arbitrary units"}


def ask_yes_no(prompt):
    """Require an explicit response before regenerating editable proposals."""
    while True:
        answer = input(prompt + " [YES/NO]: ").strip()
        if answer in ("YES", "NO"):
            return answer == "YES"
        print("Type exactly YES or NO; annotations have not changed.")


def aid_groups(state):
    """Return the two small detector states; the model/source cache is shared."""
    return (("blink_auto", state, AUTO_LABEL, "Blink"),
            ("eye_movement_auto", state["eye_movement"], EYE_AUTO_LABEL, "Eye movement"))


def initialize_aid_review(display, state, cfg, previous_review):
    """Restore per-group settings/vetoes without mutating an old receipt.

    Called before component inspection. Legacy selected_components refers only
    to blinks. Absent eye fields imply no component selection, never a completed
    review, and are not added to the loaded configuration or written to disk.
    """
    state["eye_movement"] = {key: state[key] for key in ("native_sources", "metadata", "display_gain")}
    for key, group, label, title in aid_groups(state):
        previous = deepcopy(previous_review.get(key, {}))
        group["previous_receipt"] = previous
        group["components"] = list(previous.get("selected_components", previous.get("components", [])))
        group["removed_spans_raw"] = deepcopy(previous.get("removed_spans_raw", []))
        group["previous_auto"] = [annotation_key(display, a) for a in display.annotations if a["description"] == label]
        group["settings"] = blink_settings(cfg) if key == "blink_auto" else eye_movement_settings(cfg)
        if previous.get("model_sha256") == state["metadata"]["model_sha256"]:
            group["settings"].update(previous.get("settings", {}))


def proposal_fingerprint(state, settings, *, eye_movement=False):
    """Keep the historical blink digest exactly; version only the new method."""
    identity = {"model": state["metadata"]["model_sha256"],
                "components": state["components"], "settings": settings}
    if eye_movement:
        identity["method"] = EYE_METHOD
    return digest(identity)


def update_blink_proposals(display, state, settings, *, replace=True):
    """Compatibility entry point for the unchanged blink peak detector."""
    update_aid_proposals(display, state, settings, replace=replace)


def update_aid_proposals(display, state, settings, *, eye_movement=False, replace=True):
    """Update only the chosen group's proposals, diagnostics and score display.

    NONE skips detector execution. Saved intervals and edited/deleted session
    proposals remain protected by the existing replacement helper. Fingerprint
    equality controls reseeding independently for each group, including a
    previously reviewed-empty set. Score channels are disposable display only.
    """
    import matplotlib.pyplot as plt
    label = EYE_AUTO_LABEL if eye_movement else AUTO_LABEL
    title = "Eye movement" if eye_movement else "Blink"
    score_name = "EYE_MOVEMENT_score" if eye_movement else "BLINK_score"
    if state["components"]:
        detector = detect_eye_movements if eye_movement else detect_blinks
        intervals, scores, result = detector(state["native_sources"], state["components"], settings)
    else:
        intervals, scores = [], np.full(display.n_times, np.nan)
        result = {"status": "no_component_selected", "candidate_count": 0,
                  "components": [], "settings": dict(settings), "unscored_samples": display.n_times}
        if eye_movement:
            result["method"] = EYE_METHOD
    result.setdefault("components", list(state["components"]))
    result.setdefault("settings", dict(settings))
    result["fingerprint"] = proposal_fingerprint(state, settings, eye_movement=eye_movement)
    result["provisional_intervals_raw"] = [[int(a + display.first_samp), int(b + display.first_samp)]
                                           for a, b in intervals]
    if replace:
        replace_proposals(display, intervals, state, label=label)
    state["detector"] = result
    if state.get("overview") is not None:
        plt.close(state.pop("overview"))
    previous = [ch for ch in display.ch_names if ch == score_name]
    if previous:
        display.drop_channels(previous)
    if not state["components"]:
        print(f"{title} detector disabled (NONE); saved/edited intervals remain. Inspect the complete source manually.")
        return
    info = mne.create_info([score_name], display.info["sfreq"], "misc")
    # One browser unit equals the detection threshold; the score is visual
    # guidance only and never enters a production FIF or the final ICA model.
    trace = mne.io.RawArray((np.minimum(scores / settings["threshold_std"], 2.) * state.get("display_gain", 1.))[None],
                           info, first_samp=display.first_samp, verbose="ERROR")
    trace.set_meas_date(display.info["meas_date"])
    display.load_data().add_channels([trace], force_update_info=True)
    fig, ax = plt.subplots(figsize=(12, 3))
    ax.plot(display.times, scores, linewidth=.5)
    ax.axhline(settings["threshold_std"], color="red", label="candidate threshold")
    broad = [p for p in result.get("peaks", []) if p["long_event"]]
    if broad:
        ax.scatter([(p["peak_raw"] - display.first_samp) / display.info["sfreq"] for p in broad],
                   [p["score"] for p in broad], color="orange", label="broad peak: inspect/resize")
    ax.set(xlabel="Seconds from recording start", ylabel="Robust amplitude score",
           title=f'{title} proposals: {result["candidate_count"]}; gaps are unscored, not artifact-free')
    ax.legend(); fig.tight_layout(); plt.show(block=False)
    state["overview"] = fig
    print(f'{title} detector: {result["status"]}; {result["candidate_count"]} merged proposals. Review flagged AND unflagged intervals.')
    if broad:
        print("Broad peak times from recording start (seconds):",
              [round((p["peak_raw"] - display.first_samp) / display.info["sfreq"], 3) for p in broad])
    print(f"{score_name}: 1 = candidate threshold; display capped at 2. Overview/detection scores are uncapped.")


def choose_redetection_group():
    """Ask which detector to regenerate, retrying empty/unknown responses."""
    while True:
        answer = input("Detector to update [BLINK/EYE]: ").strip()
        if answer in ("BLINK", "EYE"):
            return answer == "EYE"
        print("Type exactly BLINK or EYE; neither group's proposals changed.")


def adjust_threshold(settings, title):
    """Return a copied settings map with one validated positive threshold."""
    while True:
        text = input(f'{title} threshold in robust standard deviations (Enter keeps {settings["threshold_std"]:g}): ').strip()
        try:
            value = float(text) if text else settings["threshold_std"]
            if np.isfinite(value) and value > 0:
                return {**settings, "threshold_std": value}
        except ValueError:
            pass
        print("Enter a finite positive threshold or press Enter to keep it.")


def detector_receipt(group, label, annotations):
    """Bind one group's actual detection separately from a pending selection.

    BACK may display a new selection while regeneration is declined. components
    and fingerprint still describe the scores/proposals actually computed;
    selected_components records the current choice. Keep compact earlier method
    provenance when changing a detector while retaining previously saved work.
    """
    result = deepcopy(group["detector"])
    previous = group["previous_receipt"]
    history = deepcopy(previous.get("prior_detections", []))
    if previous.get("fingerprint") and previous["fingerprint"] != result["fingerprint"]:
        old = {key: previous[key] for key in ("fingerprint", "model_sha256", "components", "settings",
               "method", "detection_status", "provisional_intervals_raw") if key in previous}
        if not any(row["fingerprint"] == old["fingerprint"] for row in history):
            history.append(old)
    return {**result, "model_sha256": group["metadata"]["model_sha256"],
        "cache_fingerprint": group["metadata"]["fingerprint"],
        "selected_components": list(group["components"]),
        "selection_applied": group["components"] == result.get("components", []),
        "selection_status": "selected" if group["components"] else "no_component_selected",
        "detection_status": result["status"], "removed_spans_raw": deepcopy(group["removed_spans_raw"]),
        "retained_auto_count": int(sum(value == label for value in annotations.description)),
        "prior_detections": history, "status": "reviewed",
        "acceptance": "Retained automatic intervals explicitly reviewed on complete-source SAVE; NONE is not artifact absence"}


def misc_observations(available):
    """Collect source-specific MISC usefulness after viewing, retrying mistakes."""
    print("Available original MISC channels:", available)
    while True:
        text = input("Usable eye MISC channels (comma separated), or Enter for none: ").strip()
        usable = [ch.strip() for ch in text.split(",") if ch.strip()]
        if all(ch in available for ch in usable):
            break
        print("A declared channel is unavailable; enter only names shown above.")
    while True:
        note = input("Channel meanings/usability and review observations (required): ").strip()
        if note:
            return usable, note
        print("Record a brief source-specific observation, including limitations.")


def review_action(*, preliminary_ica, remaining):
    """Require an explicit action; blank/typo input never discards annotations."""
    print(f"{remaining} retained EDF candidates, including long ones, will become BAD_blink only on SAVE.")
    print("SAVE = confirm entire-source review and save; DISCARDALL = discard this session's edits.")
    print("After SAVE, finish MISC/observations prompts and wait for the final 'Saved reviewed ocular annotations' message.")
    print("SAVE retains edited BAD categories, including the complete reviewed muscle set when muscle QC is available.")
    choices = {"SAVE", "DISCARDALL"}
    if preliminary_ica:
        choices.add("BACK")
        choices.add("REDETECT")
        print("BACK = go back and select aid components again, keeping all annotation edits in memory.")
        print("REDETECT = choose BLINK or EYE, adjust its threshold and regenerate only that group's untouched proposals.")
    while True:
        answer = input("Action [" + "/".join(sorted(choices)) + "]: ").strip()
        if answer in choices:
            return answer
        print("No action taken. Type exactly " + ", ".join(sorted(choices)) + ".")


def onset_color(value):
    """Parse an opaque named/hex marker color before loading recording data.

    Return a canonical RGB hex string accepted by both rendering backends.
    Transparent colors are rejected: opacity is fixed at one, not configurable.
    """
    from matplotlib.colors import to_hex, to_rgba
    try:
        rgba = to_rgba(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("Use a named color or RGB hex such as black or #0057B8.") from error
    if rgba[3] != 1.:
        raise argparse.ArgumentTypeError("Onset colors must be fully opaque; transparency is not supported.")
    return to_hex(rgba, keep_alpha=False)


def style_onset_lines(browser, width, color=DEFAULT_ONSET_COLOR):
    """Style opaque onset/gaze markers with the same width, including on scroll.

    Installed MNE has no public event-line width parameter. This small display
    compatibility helper reads its existing event artists and uses their public
    styling methods. It never patches a browser method or changes event samples.
    Matplotlib recreates its LineCollection on scroll, so a draw callback styles
    each replacement; Qt keeps its existing InfiniteLine objects. Qt initially
    draws events and annotation fills at the same level, which lets fills tint
    markers. Ordinary onsets stay above fills and below traces; orange code 150
    is drawn above all plot layers. Matplotlib groups events in one collection,
    so a separate display-only collection raises just its gaze-break segments.
    No extra events or annotations are added to the data. Qt can raise individual
    event artists directly. Opacity is always one. Backend changes fail explicitly.
    """
    from matplotlib.colors import to_rgba
    from matplotlib.figure import Figure
    rgba = to_rgba(color, alpha=1.)
    gaze_rgba = to_rgba(GAZE_BREAK_COLOR, alpha=1.)
    front_order = 1_000_000  # Above MNE's trace, annotation, cursor and scale layers.
    if not hasattr(browser, "mne") or not hasattr(browser.mne, "event_lines"):
        raise RuntimeError("This MNE browser does not expose the expected onset artists.")
    if isinstance(browser, Figure):
        from matplotlib.collections import LineCollection
        gaze_lines = None
        def update_style(event=None):
            nonlocal gaze_lines
            lines = browser.mne.event_lines
            if lines is None:
                return
            visible = ((browser.mne.event_times >= browser.mne.times[0]) &
                       (browser.mne.event_times <= browser.mne.times[-1]))
            codes = browser.mne.event_nums[visible]
            colors = np.asarray([gaze_rgba if code == GAZE_BREAK_CODE else rgba for code in codes]).reshape(-1, 4)
            changed = False
            if (not np.allclose(lines.get_linewidths(), width) or lines.get_alpha() != 1.
                    or not np.allclose(lines.get_colors(), colors)
                    or lines.get_zorder() != browser.mne.zorder["events"]):
                lines.set_linewidth(width)
                lines.set_color(colors)
                lines.set_alpha(1.)
                lines.set_zorder(browser.mne.zorder["events"])
                changed = True
            segments = np.asarray(lines.get_segments()).reshape(-1, 2, 2)
            gaze_segments = segments[codes == GAZE_BREAK_CODE]
            if len(gaze_segments) and gaze_lines is None:
                gaze_lines = LineCollection([], colors=[gaze_rgba], linewidths=width,
                    alpha=1., zorder=front_order, label="DAD gaze-break markers")
                browser.mne.ax_main.add_collection(gaze_lines)
            if gaze_lines is not None and not np.array_equal(
                    np.asarray(gaze_lines.get_segments()).reshape(-1, 2, 2), gaze_segments):
                gaze_lines.set_segments(gaze_segments)
                changed = True
            # Queue another paint only when an event artist was recreated or
            # the visible gaze segments moved. Scrolling out clears the overlay.
            if changed:
                browser.canvas.draw_idle()
        browser.canvas.mpl_connect("draw_event", update_style)
        update_style()
    else:
        from qtpy.QtGui import QColor, QPen
        for line, code in zip(browser.mne.event_lines, browser.mne.event_nums, strict=True):
            pen = QPen(line.pen)
            pen.setWidthF(width)
            pen.setColor(QColor.fromRgbF(*(gaze_rgba if code == GAZE_BREAK_CODE else rgba)))
            line.setPen(pen)
            line.setOpacity(1.)
            line.setZValue(front_order if code == GAZE_BREAK_CODE else .5)


def prepare_annotation_descriptions(browser):
    """Offer standard manual labels without adding any annotation intervals.

    Reuse the installed browser's Add Description machinery in a small backend
    compatibility helper. Existing/custom labels remain, duplicates are skipped,
    and Qt starts with BAD_blink so hand-marked blinks need not be mislabeled as
    automatic. Matplotlib builds its usual sorted selector when 'a' is pressed;
    Qt already has an annotation dock. No data/annotation object is modified.
    """
    from matplotlib.figure import Figure
    existing = set(browser._get_annotation_labels())
    missing = [label for label in REVIEW_DESCRIPTIONS if label not in existing]
    if isinstance(browser, Figure):
        if missing:
            browser.mne.new_annotation_labels.extend(missing)
            browser._setup_annotation_colors()
            if browser.mne.fig_annotation is not None:
                browser._update_annotation_fig()
    else:
        controls = browser.mne.fig_annotation
        for label in missing:
            controls._add_description(label)
        controls.description_cmbx.setCurrentText("BAD_blink")


def plot_review_browser(display, usable_misc, edf_display=None, timing=None):
    """Put eye guidance and declared usable MISC beside the first frontal MEG traces.

    Explicit public channel ordering prevents type grouping from hiding the
    temporary EDF channels below dozens of MEG channels in the opening view.
    All remaining original channels are still available by scrolling. Accepted
    epoch annotations exist only during this browser; their edits are discarded
    on close, while artifact edits remain. The public event interface displays
    corrected onsets as vertical lines, never shifted STI101 samples.
    Onset width is independent of signal-trace thickness; the small artist
    compatibility helper runs before the usual MNE blocking browser display.
    Standard manual descriptions are available even when unused, without
    creating placeholder segments or changing what SAVE persists.
    """
    first = [ch for ch in display.ch_names if ch.startswith(("EDF_", "AID", "BLINK_", "EYE_MOVEMENT_"))]
    first.extend(ch for ch in usable_misc if ch in display.ch_names)
    display.reorder_channels(first + [ch for ch in display.ch_names if ch not in first])
    options = {}
    color_context = nullcontext()
    if timing is not None:
        preview = timing["annotations"]
        display.annotations.append(preview.onset, preview.duration, preview.description)
        onset = timing.get("onset_color", DEFAULT_ONSET_COLOR)
        codes = np.unique(timing["events"][:, 2])
        colors = ({int(code): GAZE_BREAK_COLOR if code == GAZE_BREAK_CODE else onset for code in codes}
                  if GAZE_BREAK_CODE in codes else onset)
        options.update(events=timing["events"], event_color=colors,
            event_id={"150 gaze_break" if code == GAZE_BREAK_CODE else f"Corrected onset {code}": int(code) for code in codes})
        if len(preview):
            if "annotation_colors" in inspect.signature(display.plot).parameters:
                options["annotation_colors"] = {label: color for label, color in PREVIEW_COLORS.items()
                                                 if label in preview.description}
            else:
                color_context = epoch_color_context(display.annotations)
    try:
        with color_context:
            # Construct before blocking so the event pens can be styled. Reuse
            # MNE's existing backend-aware display helper for its GUI lifecycle.
            from mne.viz.utils import _show_browser
            browser = display.plot(duration=20, n_channels=min(30, len(display.ch_names)),
                order=np.arange(len(display.ch_names)),
                scalings={"mag": "auto", "grad": "auto",
                          "misc": edf_display["display_gain"] if edf_display else "auto"},
                show=False, block=False, **options)
            prepare_annotation_descriptions(browser)
            if timing is not None:
                style_onset_lines(browser, timing.get("onset_linewidth", 2.), timing.get("onset_color", DEFAULT_ONSET_COLOR))
            _show_browser(show=True, block=True, fig=browser)
            return browser
    finally:
        # A moved/resized/deleted preview cannot become a scientific decision.
        # Its non-BAD label is also excluded by eye_annotations at SAVE.
        display.annotations.delete([i for i, label in enumerate(display.annotations.description)
                                    if label in PREVIEW_COLORS])


def main(argv=None):
    """Review one native source, optionally starting from a verified triage mask.

    Triage masks cannot replace existing manual work. They omit EDF proposals
    and retain exact detector/deletion provenance only after whole-source SAVE.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for flag in ("config", "subject-config", "tuning-config", "source", "basis"):
        parser.add_argument("--" + flag, required=True)
    parser.add_argument("--site-config")
    parser.add_argument("--preliminary-ica", action="store_true", help="Load frontal ICA; select blink then eye-movement components, or NONE independently. Never fit locally.")
    parser.add_argument("--triage-proposals", help="Exact provisional mask from annotation_triage.py; no EDF or detector reseeding.")
    parser.add_argument("--no-edf-traces", action="store_true", help="Show EDF candidate labels without loading gaze/pupil samples.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--recover-annotations", action="store_true",
                        help="Recheck an existing annotation file without reseeding EDF candidates; requires --overwrite.")
    parser.add_argument("--browser-backend", choices=("qt", "matplotlib"), default="qt")
    parser.add_argument("--no-timing-marks", action="store_true",
                        help="Omit onset/code-150 lines and temporary valid/invalid epoch previews (shown by default).")
    parser.add_argument("--timing-review", help="Explicit accepted Stage 3 review JSON; otherwise configured path or accepted s3-v002.")
    parser.add_argument("--onset-linewidth", type=float, default=2.,
                        help="Corrected onset-line width, independent of traces (default 2; Matplotlib points / Qt pixels).")
    parser.add_argument("--onset-color", type=onset_color, default=DEFAULT_ONSET_COLOR,
                        help="Opaque onset marker/label color: name or RGB hex (default black).")
    parser.add_argument("--high-pass-view", type=float, metavar="HZ",
                        help="High-pass only frontal MEG in the browser copy; MISC and source files are unchanged.")
    args = parser.parse_args(argv)
    if args.triage_proposals and (args.preliminary_ica or args.recover_annotations):
        parser.error("--triage-proposals cannot combine with --preliminary-ica or --recover-annotations")
    if args.recover_annotations and not args.overwrite:
        parser.error("--recover-annotations requires --overwrite")
    if not np.isfinite(args.onset_linewidth) or args.onset_linewidth <= 0:
        parser.error("--onset-linewidth must be finite and greater than zero")
    if args.high_pass_view is not None and (not np.isfinite(args.high_pass_view) or args.high_pass_view <= 0):
        parser.error("--high-pass-view must be finite and greater than zero")
    cfg, *_ = load_and_resolve_config(args.config, args.subject_config, args.tuning_config,
                                      site_config_path=args.site_config, required_sections=("ica", "artifact_annotation") +
                                      (() if args.no_timing_marks else ("epoching",)))
    dirs = stage_output_dirs(cfg)
    require_review(cfg, "03_A_MaxFilt", dirs["reports"])
    matches = [r for r in discover_raw_run_records(cfg) if r.source_filename == args.source]
    if len(matches) != 1:
        raise ValueError("--source must identify exactly one configured acquisition filename.")
    record = matches[0]
    source_path = derivative_path_for_run(record, dirs["sss"], "_raw_sss.fif")
    raw = mne.io.read_raw_fif(source_path, preload=False)
    # Validate triage provenance before creating any output directory. Loading
    # never promotes a proposal or changes an existing manual review receipt.
    triage = None
    if args.triage_proposals:
        from annotation_triage import load_review_mask
        triage = load_review_mask(cfg, record, args.triage_proposals)
    output = prepare_output_path(manual_path(cfg, record), raw_root=cfg["resolved_paths"]["data_root"], overwrite=args.overwrite)
    if args.recover_annotations and not output.is_file():
        raise ValueError("No existing annotation file to recover.")
    # Start from current muscle evidence. A reviewed snapshot replaces this set
    # on reopening, so deleted automatic intervals cannot reappear.
    muscle, muscle_basis = load_muscle_context(cfg, record, raw, source_path, dirs)
    if muscle is not None:
        raw.annotations.delete([i for i, label in enumerate(raw.annotations.description) if label.lower() == "bad_muscle"])
        raw.annotations.append(muscle.onset, muscle.duration, muscle.description)
    # Existing eye work may be reopened after a changed SSS decision; it is a
    # candidate until the operator rechecks it and saves a newly bound receipt.
    if output.exists():
        previous = mne.read_annotations(output)
        if previous.orig_time != raw.annotations.orig_time:
            raise ValueError("Saved eye annotations have a different acquisition origin.")
        previous_review = cfg["artifact_annotation"].get("manual_eye_reviews", {}).get(record.source_filename, {})
        # An orphan saved file can already contain manually revised muscle
        # intervals. Recovery displays that set, without duplicating the older
        # provisional detector intervals. An orphan without muscle labels cannot
        # establish reviewed zero; retain provisional context in that case.
        recovered_muscle = args.recover_annotations and any(
            str(label).lower() == "bad_muscle" for label in previous.description)
        if muscle_basis is not None and (recovered_muscle or
                reviewed_muscle_override(previous_review, cfg["artifact_annotation"],
                                        threshold=muscle_basis["threshold"])):
            raw.annotations.delete([i for i, label in enumerate(raw.annotations.description) if label.lower() == "bad_muscle"])
        # Stage 4A owns SQUID windows. Old snapshots may contain display copies;
        # retain current SSS context once and revalidate only editable 4B work.
        from squid_jumps import without_owned_context
        previous = without_owned_context(previous)
        raw.annotations.append(previous.onset, previous.duration, previous.description)
        if args.recover_annotations:
            print(f"Recovering {len(previous)} saved annotations as UNREVIEWED work; deleted EDF candidates will not be reseeded.")
            print("Recheck the current source/guidance before SAVE; recovery alone records no acceptance.")
    if triage is not None:
        from annotation_triage import seed_review_mask
        seed_review_mask(raw, triage)
        print("Loaded exact triage intervals. Review flagged, unflagged and unscored periods before SAVE.")
    proposal = None
    existing_review = cfg["artifact_annotation"].get("manual_eye_reviews", {}).get(record.source_filename, {})
    skip_edf = triage is not None or existing_review.get("edf_guidance_policy") == "not_used_annotation_triage"
    if not skip_edf and cfg["artifact_annotation"].get("edf_blink_candidates", {}).get("enabled", False):
        proposal, candidates = load_candidates(cfg, record, raw, source_path)
        previous_review = cfg["artifact_annotation"].get("manual_eye_reviews", {}).get(record.source_filename, {})
        if not args.recover_annotations and previous_review.get("edf_candidates_fingerprint") != proposal["fingerprint"]:
            raw.annotations.append(candidates.onset, candidates.duration, candidates.description)
        print(f'EDF guidance: {proposal["status"]}; {proposal["candidate_count"]} candidate bouts, '
              f'{proposal["long_candidate_count"]} long. Unmatched trials/gaps require manual inspection.')
    frontal = {ch for group in frontal_channel_groups(raw).values() for ch in group}
    channels = [ch for ch, kind in zip(raw.ch_names, raw.get_channel_types())
                if kind == "misc" or ch in frontal]
    # Validate the small cache/headers before the expensive EDF/display read.
    # A missing/stale cache must fail without loading the full frontal signals.
    loaded_aid = load_cached_aid(cfg, record) if args.preliminary_ica else None
    timing = None if args.no_timing_marks else load_timing_preview(cfg, record, raw, args.timing_review)
    if timing is not None:
        timing["onset_linewidth"] = args.onset_linewidth
        timing["onset_color"] = args.onset_color
    display = raw.copy().pick(channels)
    if args.high_pass_view is not None:
        # Filter only this disposable browser copy; raw remains the source of
        # truth for annotations. Explicit picks keep MISC and auxiliary traces
        # out of the filter while preserving them in the browser unchanged.
        frontal_meg_channels = [ch for ch in channels if ch in frontal]
        if frontal_meg_channels:
            print(f"High-pass filtering browser frontal MEG at {args.high_pass_view:g} Hz; MISC stays unfiltered.")
            display.load_data().filter(l_freq=args.high_pass_view, h_freq=None, picks=frontal_meg_channels)
        else:
            print("No frontal MEG channels are available to high-pass; MISC stays unfiltered.")
    edf_display = add_edf_display(display, cfg, proposal) if proposal and not args.no_edf_traces else None
    mne.viz.set_browser_backend(args.browser_backend)
    available_misc = [ch for ch in channels if ch.startswith("MISC")]
    print("Available analog channels:", available_misc)
    if "MISC008" in available_misc:
        print("MISC008 is photodiode context; its presence does not make it an eye signal.")
    # All original MISC traces are initially available. Ask usefulness only
    # after the user has actually inspected them, immediately before saving.
    usable = available_misc
    aid = None
    aid_state = loaded_aid if loaded_aid is not None else {}
    previous_review = cfg["artifact_annotation"].get("manual_eye_reviews", {}).get(record.source_filename, {})
    manual_auto_states = {key: {"previous_auto": [annotation_key(display, a) for a in display.annotations
        if a["description"] == label], "removed_spans_raw": deepcopy(previous_review.get(key, {}).get("removed_spans_raw", []))}
        for key, label in (("blink_auto", AUTO_LABEL), ("eye_movement_auto", EYE_AUTO_LABEL))}
    if triage is not None:
        # Keep detector/model provenance and future deletion vetoes when this
        # exact provisional set becomes a genuine complete-source manual SAVE.
        previous_review = deepcopy(previous_review)
        for key, detector_key in (("blink_auto", "blink"), ("eye_movement_auto", "eye_movement")):
            if detector_key in triage["detectors"]:
                previous_review[key] = deepcopy(triage["detectors"][detector_key])
    if args.preliminary_ica:
        aid_state["display_gain"] = edf_display["display_gain"] if edf_display else 1.
        initialize_aid_review(display, aid_state, cfg, previous_review)
        aid = add_preliminary_sources(display, cfg, state=aid_state)
        for key, group, label, title in aid_groups(aid_state):
            eye = key == "eye_movement_auto"
            same_proposals = group["previous_receipt"].get("fingerprint") == proposal_fingerprint(
                group, group["settings"], eye_movement=eye)
            update_aid_proposals(display, group, group["settings"], eye_movement=eye,
                                 replace=not same_proposals and not args.recover_annotations)
            if same_proposals:
                print(f"{title}: unchanged reviewed proposals were not reseeded; prior edits/deletions remain.")
        if args.recover_annotations:
            print("Recovery kept saved ocular edits/deletions; diagnostics added no proposals. REDETECT is explicit for each group.")
    print("Inspect MISC availability, meaning, flatness, clipping and dropouts; no EOG assumptions.")
    if not skip_edf:
        print("Edit/delete EDF candidates; add BAD_blink and BAD_eye_movement where needed. Review the entire source, including EDF gaps.")
    else:
        print("Edit the triage BAD intervals and add missing artifacts. EDF and SQUID proposals are outside this session.")
    print("BAD_blink, BAD_eye_movement and BAD_segment are already available in the annotation selector; unused labels save no intervals.")
    print("BAD_squid_jumps is Stage 4A context: browser edits to it are not saved; revise it through SQUID review.")
    if muscle_basis is not None:
        print("BAD_muscle can also be moved, resized, deleted or added. SAVE makes this full set the reviewed replacement, including zero.")
    print("EDF_*_norm traces are normalized visual guidance; EDF_coverage=0 and EDF_no_guidance denote unavailable coverage.")
    if timing is not None:
        print("Green EPOCH_preview = eligible trial; orange EPOCH_invalid_preview = already excluded gaze-break trial's intended window.")
        print("Orange code 150 stays at its recorded sample, above traces/annotations; no continuous coordinates are shifted.")
        print("Press a, then Select Visible (Qt) or label checkboxes (Matplotlib) to toggle either preview description.")
        print("Keep both preview labels unchanged; their edits are discarded at browser close and neither is saved.")
    # Reopen through the same browser lifecycle; only onset artist width is styled.
    while True:
        plot_review_browser(display, usable, edf_display, timing)
        if args.preliminary_ica:
            import matplotlib.pyplot as plt
            for key, group, label, title in aid_groups(aid_state):
                if group.get("overview") is not None:
                    plt.close(group.pop("overview"))
        remaining = sum(label in CANDIDATE_LABELS for label in display.annotations.description)
        action = review_action(preliminary_ica=args.preliminary_ica, remaining=remaining)
        if action == "BACK":
            # Capture each group's edits BEFORE changing its selection. Groups
            # with unchanged selections keep their current diagnostics/proposals.
            old_selections = {}
            for key, group, label, title in aid_groups(aid_state):
                remember_deletions(display, group, label=label)
                old_selections[key] = list(group["components"])
            aid = add_preliminary_sources(display, cfg, state=aid_state)
            for key, group, label, title in aid_groups(aid_state):
                if group["components"] == old_selections[key]:
                    continue
                if not group["components"] or ask_yes_no(
                        f"Regenerate {title.lower()} proposals from this selection, preserving edits/deletions?"):
                    update_aid_proposals(display, group, group["settings"], eye_movement=key == "eye_movement_auto")
                else:
                    print(f"{title}: new traces selected; existing scores/proposals retain their previous detector provenance. Use REDETECT to apply the selection.")
        elif action == "REDETECT":
            eye = choose_redetection_group()
            group = aid_state["eye_movement"] if eye else aid_state
            if group["components"]:
                group["settings"] = adjust_threshold(group["settings"], "Eye movement" if eye else "Blink")
            update_aid_proposals(display, group, group["settings"], eye_movement=eye)
        elif action == "DISCARDALL":
            print("Session edits discarded. Existing annotation files and tuning receipts were not changed.")
            return 0
        else:
            break
    usable, mapping_note = misc_observations(available_misc)
    if args.preliminary_ica:
        for key, group, label, title in aid_groups(aid_state):
            remember_deletions(display, group, label=label)
    promote_candidates_after_review(display, confirmed=True)
    ann = eye_annotations(display)
    ann.save(output, overwrite=args.overwrite)
    receipt = {"status": "reviewed", "reviewed_utc": utc_now(), "basis": args.basis,
               "coordinates": coordinates(raw), "sss_input": file_stamp(source_path),
               "annotation_count": len(ann), "annotation_sha256": file_hash(output),
               "edf_candidates_fingerprint": proposal["fingerprint"] if proposal else None,
               "edf_display": edf_display,
               "muscle_review": ({**muscle_basis, "status": "reviewed", "annotation_count":
                   int(sum(label.lower() == "bad_muscle" for label in ann.description))} if muscle_basis is not None else None),
               "timing_preview": timing["metadata"] if timing is not None else {"enabled": False},
               "preliminary_ica": aid, "usable_misc_channels": usable, "channel_meanings_and_limits": mapping_note}
    if args.preliminary_ica:
        for key, group, label, title in aid_groups(aid_state):
            receipt[key] = detector_receipt(group, label, ann)
    else:
        # Manual-only review must not erase earlier detector/model provenance.
        # Do not invent an eye field for an old blink-only receipt.
        receipt["preliminary_ica"] = deepcopy(previous_review.get("preliminary_ica"))
        for key, label in (("blink_auto", AUTO_LABEL), ("eye_movement_auto", EYE_AUTO_LABEL)):
            if key in previous_review:
                remember_deletions(display, manual_auto_states[key], label=label)
                receipt[key] = {**deepcopy(previous_review[key]),
                    "removed_spans_raw": manual_auto_states[key]["removed_spans_raw"],
                    "retained_auto_count": int(sum(value == label for value in ann.description))}
    if triage is not None:
        receipt["triage_origin"] = {"payload_sha256": triage["payload_sha256"],
            "source_mask": str(Path(args.triage_proposals).resolve()), "purpose": "manually_reviewed_triage_starting_intervals"}
        receipt["edf_guidance_policy"] = "not_used_annotation_triage"
    else:
        for key in ("triage_origin", "edf_guidance_policy"):
            if key in existing_review:
                receipt[key] = deepcopy(existing_review[key])
    # Send only this exact-source receipt. The helper reloads the current map
    # under its write lock, so an older GUI cannot erase another run's new SAVE.
    try:
        update_tuning_file(path=args.tuning_config, raw_root=cfg["resolved_paths"]["data_root"],
                           subject_id=cfg["subject"]["id"], section="artifact_annotation",
                           mapping_key="manual_eye_reviews", values={record.source_filename: receipt})
    except Exception:
        print(f"Annotation FIF saved at {output}, but its tuning receipt failed. Preserve the file; use --overwrite --recover-annotations to recover it.")
        raise
    print(f"Saved reviewed ocular annotations ({len(ann)}): {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
