"""Small Stage 4 helpers for coordinates, decisions and protected filtering.

Usage
-----
    fit_raw, intervals = prepare_fit_copy(annotated_raw, ica_config)
    summary = filter_continuous(cleaned_native_raw, preprocessing_config)
    # Muscle overrides use the same detector contract as the reviewed browser.
    use_reviewed_set = reviewed_muscle_override(receipt, artifact_config, threshold)

Inputs are MNE Raw objects and the existing YAML sections. Fit preparation
returns a disposable MEG-only object and exact native-source interval records.
Final filtering modifies only MEG samples in place, preserving native time and
auxiliary signals. Historical exact-source exceptions protect only a future
additional high-pass; spectrum-fit runs for every source, including P38 run 1.
No function discovers, writes or processes a cohort itself.

Every BAD-prefixed annotation is excluded, case-insensitively. Filtering and
resampling operate on separate admissible stretches. The 1 Hz fit filter is
independent of final processing, which defaults to spectrum-fit only. No FIR
margin is needed without an additional high-pass. Historical filter-edge
exposure remains information, never an epoch-exclusion rule.
Large FIFs use portable stat fingerprints;
small decision, annotation and ICA-model files use content hashes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mne
import numpy as np


# Only these downstream settings are removed from the SSS review identity.
# The frozen former values reproduce the 2026-10-23 identity for read-only
# compatibility; no historical receipt, reason or report is rewritten.
LEGACY_FINAL_FILTER = {"highpass_hz": 0.01, "lowpass_hz": None,
    "notch_freqs_hz": [50.0, 100.0], "notch_method": "spectrum_fit",
    "spectrum_fit_window": "10s", "resample_sfreq_hz": None}
FINAL_FILTER_KEYS = (*LEGACY_FINAL_FILTER, "filter_edge_margin_sec", "final_filter_skip_sources")


# -----------------------------------------------------------------------------
# Portable evidence and exact acquisition coordinates
# -----------------------------------------------------------------------------

def digest(value):
    """Hash JSON-compatible scientific settings, independent of mount paths."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def file_hash(path):
    """Hash a small resource/model/annotation file without loading it all at once."""
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def file_stamp(path):
    """Identify a large shared FIF across local and cluster mount points."""
    path = Path(path)
    stat = path.stat()
    return {"name": path.name, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def annotation_rows(raw):
    """Serialize acquisition-time annotations, including their time origin."""
    return {"orig_time": str(raw.annotations.orig_time),
            "rows": [[float(a["onset"]), float(a["duration"]), str(a["description"]),
                      list(a.get("ch_names", ())), dict(a.get("extras", {}))] for a in raw.annotations]}


def coordinates(raw):
    """Return the native coordinates and channel basis that must be preserved."""
    return {"sfreq": float(raw.info["sfreq"]), "n_times": int(raw.n_times),
            "first_samp": int(raw.first_samp), "meas_date": str(raw.info["meas_date"]),
            "channels": raw.ch_names, "types": raw.get_channel_types()}


def check_preservation(source, output, *, annotations=True, samples=True):
    """Check native coordinates and STI/MISC values, allowing FIF single rounding.

    Annotations may be added by 03_B/final edge handling; original rows must
    still be present. The tolerance is only for float32 FIF serialization,
    never for changes to integer STI101 values or sample coordinates.
    """
    if coordinates(source) != coordinates(output):
        raise RuntimeError("Continuous recording coordinates/channel identities changed.")
    if source.info["bads"] != output.info["bads"]:
        raise RuntimeError("Bad-channel metadata changed unexpectedly.")
    if annotations and annotation_rows(source) != annotation_rows(output):
        raise RuntimeError("Annotations changed unexpectedly.")
    if not annotations:
        for original in source.annotations:
            matches = [a for a in output.annotations
                       if a["description"] == original["description"]
                       and np.isclose(a["onset"], original["onset"], atol=1e-6, rtol=0)
                       and np.isclose(a["duration"], original["duration"], atol=1e-6, rtol=0)
                       and tuple(a.get("ch_names", ())) == tuple(original.get("ch_names", ()))
                       and dict(a.get("extras", {})) == dict(original.get("extras", {}))]
            if not matches:
                raise RuntimeError("An original annotation was lost or shifted.")
    if samples:
        picks = [i for i, kind in enumerate(source.get_channel_types()) if kind in ("stim", "misc")]
        for start in range(0, source.n_times, 100000):
            before = source.get_data(picks=picks, start=start, stop=start + 100000)
            after = output.get_data(picks=picks, start=start, stop=start + 100000)
            for row, pick in enumerate(picks):
                if source.get_channel_types()[pick] == "stim":
                    valid = np.array_equal(before[row], after[row])
                else:
                    valid = np.allclose(before[row], after[row], rtol=2e-7, atol=1e-12, equal_nan=True)
                if not valid:
                    raise RuntimeError(f"Auxiliary samples changed: {source.ch_names[pick]}")
    return {**coordinates(output), "auxiliary_samples_verified": bool(samples)}


def require_review(cfg, stage, report_dir):
    """Require a human-accepted evidence fingerprint in the existing tuning YAML."""
    receipt = cfg.get("stage4_review", {}).get(stage, {})
    report = Path(report_dir) / f"{stage}_review_report.json"
    if not report.is_file():
        raise RuntimeError(f"Missing prerequisite review report: {report}")
    evidence = json.loads(report.read_text()).get("evidence_sha256")
    if not evidence or receipt.get("evidence_sha256") != evidence or receipt.get("status") != "accepted":
        raise RuntimeError(f"{stage} is unreviewed or changed; reconcile and accept its exact evidence first.")
    # The report alone is insufficient if files changed after acceptance.
    payload = json.loads(report.read_text())
    for row in payload.get("evidence_files", []):
        path = Path(report_dir).parents[2] / row["relative_path"]
        if file_stamp(path) != row["stamp"]:
            raise RuntimeError(f"Accepted Stage 4 evidence changed: {path}")
    for row in payload.get("source_files", []):
        path = Path(cfg["resolved_paths"]["data_root"]) / row["relative_path"]
        if file_stamp(path) != row["stamp"]:
            raise RuntimeError(f"Accepted raw source changed: {path}")
    if stage == "03_B_ArtefactAnnotation" and payload.get("source_association_sha256"):
        # Additive stage lineage may change this manifest; only the accepted
        # source-association portion defines the EDF trial correspondence.
        manifest = json.loads((Path(report_dir) / "source_manifest.json").read_text())
        associations = manifest.get("trigger_derived_associations", {}).get("meg_acquisition_run_to_behavioral_segments", [])
        if digest(associations) != payload["source_association_sha256"]:
            raise RuntimeError("Accepted source associations changed; regenerate EDF guidance and review it.")
    if stage == "03_B_ArtefactAnnotation":
        # Explicit adoption is independent of manual SAVE provenance. Recheck
        # frozen source/manual/model/timing bindings whenever it is consumed.
        from stage4b_adoption import require_operative_adoption
        require_operative_adoption(cfg, payload, require_publication=
            receipt.get("acceptance_kind") == "explicit_annotation_policy_for_processing")
    accepted_config = receipt.get("configuration_sha256")
    if stage == "03_A_MaxFilt" and accepted_config != payload.get("configuration_sha256"):
        raise RuntimeError("Accepted SSS report/configuration binding changed.")
    compatible = accepted_config == stage_configuration(cfg, stage)
    if stage == "03_A_MaxFilt" and not compatible:
        # Reconstruct the former final-filter portion only. Current Maxwell,
        # calibration contents and every other preprocessing field still have
        # to reproduce the accepted fingerprint; all file checks above remain.
        compatible = accepted_config == legacy_stage4a_configuration(cfg)
    if not compatible:
        raise RuntimeError(f"Accepted {stage} configuration changed; review its affected evidence again.")
    return receipt


def _sss_configuration(cfg):
    """Bind SSS settings and nonempty SQUID actions, preserving empty reviews.

    Detector/UI settings do not alter SSS. Actual interval/channel decisions
    do; they cannot use the final-temporal-only legacy compatibility exception.
    """
    from squid_jumps import application
    settings = dict(cfg["maxwell_filter"])
    for key in ("calibration_file", "crosstalk_file"):
        settings[key] = file_hash(settings[key])
    filters = {k: v for k, v in cfg["preprocessing"].items() if k not in FINAL_FILTER_KEYS}
    result = {"maxwell": settings, "continuous_filter": filters}
    squid = application(cfg)
    if squid:
        result["squid_jumps"] = squid
    return result


def legacy_stage4a_configuration(cfg):
    """Reproduce only the known historical 4A projection, without accepting it.

    Used solely by the prerequisite checker after report/receipt and evidence
    validation. Unknown non-temporal keys and changed SSS settings still fail.
    The superseded temporal policy is never used to process a recording.
    """
    settings = _sss_configuration(cfg)
    settings["continuous_filter"].update(LEGACY_FINAL_FILTER)
    return digest(settings)


def final_filter_configuration(settings, *, include_edge=True):
    """Project final operations independently of annotation and ICA fit settings."""
    result = {k: settings.get(k) for k in FINAL_FILTER_KEYS}
    # A leftover historical margin has no operational meaning without a FIR.
    if not include_edge or settings.get("highpass_hz") is None:
        result["filter_edge_margin_sec"] = None
    return result


def stage_configuration(cfg, stage, *, include_edge=True):
    """Bind SSS, annotation/final policy and fit identity to their own stages."""
    if stage == "03_A_MaxFilt":
        return digest(_sss_configuration(cfg))
    if stage == "03_B_ArtefactAnnotation":
        settings = dict(cfg["artifact_annotation"])
        default = settings.get("muscle_threshold")
        settings["muscle_threshold_by_run"] = {str(k): v for k, v in settings.get("muscle_threshold_by_run", {}).items() if v != default}
        return digest({"artifact": settings,
                       "final_filter": final_filter_configuration(cfg["preprocessing"], include_edge=include_edge)})
    # Component exclusions are a subsequent decision; fit settings define model identity.
    settings = {k: v for k, v in cfg["ica"].items() if k not in ("manual_exclude", "reviewed_model_sha256")}
    return digest(settings)


def bind_report(payload, paths, settings, report_dir, derivatives_root):
    """Add a small portable evidence record to an existing stage report."""
    root = Path(derivatives_root).resolve()
    rows = []
    for path in paths:
        path = Path(path).resolve()
        rows.append({"relative_path": str(path.relative_to(root)), "stamp": file_stamp(path)})
    payload["evidence_files"] = rows
    payload["evidence_sha256"] = digest({"files": rows, "settings": settings})


# -----------------------------------------------------------------------------
# BAD/EDGE masks and isolated, fully mapped fit stretches
# -----------------------------------------------------------------------------

def muscle_parameters(settings, threshold):
    """Project the existing muscle-detection settings for source-bound review.

    Use the same contract for provisional detection, local edits and rebuilds.
    Per-run thresholds are resolved by the caller; unrelated display/eye settings
    do not invalidate a reviewed muscle replacement.
    """
    return {"threshold": float(threshold), "muscle_settings": {
        key: value for key, value in settings.items()
        if key.startswith("muscle_") and key != "muscle_threshold_by_run"}}


def reviewed_muscle_override(receipt, settings, threshold):
    """Recognize a reviewed full muscle set, including an explicitly empty set.

    Legacy eye-only receipts lack this field and keep automatic muscle labels.
    A changed detector contract requires explicit reconciliation before reuse;
    the existing receipt separately binds SSS identity and native coordinates.
    """
    review = (receipt or {}).get("muscle_review")
    if not review:
        return False
    if (review.get("status") != "reviewed" or
            review.get("settings_sha256") != digest(muscle_parameters(settings, threshold))):
        raise RuntimeError("Muscle settings changed after manual review; reconcile/reset this source's muscle review before rebuilding.")
    return True


def good_stretches(raw):
    """Return half-open native indices outside all BAD intervals and EDGE cuts.

    Raw.set_annotations has already converted onsets to the raw acquisition
    origin, whether meas_date exists or not. Subtract first_time exactly once.
    Zero-duration BAD/EDGE boundaries split filtering without deleting samples.
    """
    good = np.ones(raw.n_times, dtype=bool)
    cuts = {0, raw.n_times}
    sfreq = raw.info["sfreq"]
    for ann in raw.annotations:
        label = ann["description"].lower()
        if not label.startswith(("bad", "edge")):
            continue
        start = int(np.floor((ann["onset"] - raw.first_time) * sfreq + 1e-6))
        stop = int(np.ceil((ann["onset"] + ann["duration"] - raw.first_time) * sfreq - 1e-6))
        start, stop = max(0, min(raw.n_times, start)), max(0, min(raw.n_times, stop))
        cuts.update((start, stop))
        if label.startswith("bad"):
            good[start:stop] = False
    transitions = np.flatnonzero(np.diff(np.r_[False, good, False])).tolist()
    cuts.update(transitions)
    ordered = sorted(cuts)
    return [(a, b) for a, b in zip(ordered[:-1], ordered[1:]) if b > a and good[a:b].all()]


def filter_design(sfreq, l_freq, h_freq=None):
    """Report the actual default zero-phase FIR support, not a guessed margin."""
    coefficients = mne.filter.create_filter(None, sfreq, l_freq, h_freq, verbose="ERROR")
    return {"method": "fir", "phase": "zero", "fir_design": "firwin",
            "n_taps": len(coefficients), "half_support_sec": (len(coefficients) - 1) / (2 * sfreq)}


def prepare_fit_copy(raw, ica_cfg):
    """Filter/resample independent good stretches and return mapped fit data.

    Full FIR half-support is removed at both ends; a further 0.1 s per end is
    removed after polyphase resampling to avoid relying on its boundary padding.
    These are fit-only exclusions and never rewrite trial eligibility.
    """
    sfreq = float(raw.info["sfreq"])
    design = filter_design(sfreq, ica_cfg["fit_l_freq_hz"], ica_cfg["fit_h_freq_hz"])
    guard = int(np.ceil(design["half_support_sec"] * sfreq))
    target = float(ica_cfg["fit_resample_sfreq_hz"])
    resample_guard = int(np.ceil(float(ica_cfg.get("fit_resample_edge_sec", 0.1)) * target))
    minimum = float(ica_cfg.get("fit_min_good_sec", 1.0))
    copies, intervals = [], []
    for start, stop in good_stretches(raw):
        if (stop - start - 2 * guard) / sfreq < minimum + 2 * resample_guard / target:
            intervals.append({"start_raw": start + raw.first_samp, "stop_raw_exclusive": stop + raw.first_samp,
                              "status": "too_short_after_filter_edges", "fit_samples": 0})
            continue
        segment = raw.copy().crop(start / sfreq, (stop - 1) / sfreq).pick("meg").load_data()
        segment.filter(ica_cfg["fit_l_freq_hz"], ica_cfg["fit_h_freq_hz"],
                       skip_by_annotation=("edge", "bad"), verbose="ERROR")
        segment.crop(guard / sfreq, (segment.n_times - guard - 1) / sfreq)
        segment.resample(target, method="polyphase", verbose="ERROR")
        if resample_guard:
            segment.crop(resample_guard / target, (segment.n_times - resample_guard - 1) / target)
        # Only disposable fit copies lose source-relative annotations; no BAD
        # sample was selected. Exact provenance is retained in the interval map.
        segment.set_annotations(None)
        copies.append(segment)
        intervals.append({"start_raw": start + raw.first_samp, "stop_raw_exclusive": stop + raw.first_samp,
                          "retained_start_raw": start + guard + raw.first_samp,
                          "retained_stop_raw_exclusive": stop - guard + raw.first_samp,
                          "resampled_first_samp": int(segment.first_samp),
                          "resampled_stop_exclusive": int(segment.first_samp + segment.n_times),
                          "resample_edge_sec": resample_guard / target,
                          "status": "included", "fit_samples": segment.n_times})
    if not copies:
        raise RuntimeError("No ICA fit stretches remain after BAD intervals and filter-edge exclusions.")
    fitted = mne.concatenate_raws(copies, on_mismatch="raise", verbose="ERROR")
    expected = sum(row["fit_samples"] for row in intervals)
    if fitted.n_times != expected:
        raise RuntimeError("ICA fit interval/sample reconciliation failed.")
    return fitted, {"filter": design, "native_coordinates": coordinates(raw),
                    "intervals": intervals, "fit_samples": expected,
                    "fit_duration_sec": expected / target}


# -----------------------------------------------------------------------------
# Final native-resolution filtering and informational edge exposure
# -----------------------------------------------------------------------------

def final_filter_summary(raw, settings, *, source_filename=None):
    """Describe independent high-pass and line-removal operations for review.

    A null high-pass has no FIR design or required margin. The historical
    ``final_filter_skip_sources`` map now protects only an additional high-pass:
    the user clarified on 2026-10-28 that even P38 run 1 receives spectrum-fit.
    Preserve its exact-source warning/reason for any future high-pass choice.
    Durations expose spectrum-fit boundary/short-stretch concerns without
    inventing a FIR margin or an epoch-exclusion policy for that operation.
    """
    exceptions = settings.get("final_filter_skip_sources", {})
    if not isinstance(exceptions, dict) or any(
            Path(name).name != name or not name.endswith(".fif")
            or not isinstance(reason, str) or not reason.strip()
            for name, reason in exceptions.items()):
        raise ValueError("Final-filter exceptions require exact FIF basenames and nonempty reasons.")
    if exceptions and source_filename is None:
        raise ValueError("An exact source_filename is required for final-filter exceptions.")
    highpass = settings.get("highpass_hz")
    skip_highpass = source_filename in exceptions
    design = filter_design(raw.info["sfreq"], highpass) if highpass is not None and not skip_highpass else None
    stretches = good_stretches(raw)
    margin = settings.get("filter_edge_margin_sec") if design is not None else None
    candidate = (design["half_support_sec"] if margin is None else float(margin)) if design is not None else None
    sfreq = raw.info["sfreq"]
    edge = int(np.ceil(candidate * sfreq)) if candidate is not None else 0
    usable = sum(max(0, b - a - 2 * edge) for a, b in stretches)
    good_seconds = sum(b - a for a, b in stretches) / sfreq
    return {"status": ("spectrum_fit_only_requires_review" if design is None else
                       "filter_design_review_required" if usable == 0 else "edge_exposure_requires_review"),
            "source_filename": source_filename,
            "highpass_status": ("skipped_by_reviewed_source_decision" if skip_highpass and highpass is not None
                                else "not_selected" if highpass is None else "selected"),
            "source_highpass_warning": ({"historical_reason": exceptions[source_filename],
                "current_policy": "Spectrum-fit still runs. Revisit this short-source warning before any future additional high-pass."}
                if skip_highpass else None),
            "spectrum_fit": {"method": "spectrum_fit", "freqs_hz": settings.get("notch_freqs_hz", [50., 100.]),
                "window": settings.get("spectrum_fit_window", "10s"), "status": "selected_for_every_source",
                "boundary_policy": "separate_BAD_EDGE_stretches; inspect boundaries and short stretches"},
            "filter": design, "reviewed_margin_sec": margin,
            "candidate_margin_sec": candidate, "good_stretch_count": len(stretches),
            "good_stretch_durations_sec": [(b - a) / sfreq for a, b in stretches],
            "seconds_outside_edge_exposure": usable / sfreq,
            "good_seconds_before_edge_exposure": good_seconds,
            "edge_exposed_good_seconds": good_seconds - usable / sfreq,
            "edge_policy": "informational_only_keep_eligible_epochs",
            "epoch_exclusions_from_filter_edges": 0,
            "settings": dict(settings)}


def make_filter_edges_informational(raw):
    """Keep legacy filter-edge labels from automatically rejecting epochs.

    Rename only the exact historical BAD_filter_edge category, in memory.
    Sensor faults, muscle/ocular artifacts and all other BAD labels retain their
    usual rejection behavior. New continuous outputs use the informational name.
    """
    mapping = {str(name): "filter_edge_exposure" for name in raw.annotations.description
               if str(name).casefold() == "bad_filter_edge"}
    if mapping:
        raw.annotations.rename(mapping)
    return len(mapping)


def filter_continuous(raw, settings, *, source_filename=None):
    """Apply spectrum-fit to every source and an optional high-pass to MEG only.

    All BAD/EDGE boundaries split both operations. Skipped BAD samples remain
    present and annotated. Only an applied high-pass needs a reviewed FIR margin.
    Its annotations are informational and cannot reject otherwise eligible
    epochs. Historical exact-source exceptions protect only a future high-pass;
    they never skip independent spectrum-fit. ICA fit copies remain independent.
    Spectrum-fit boundary/short-stretch concerns still need review in 4B.
    """
    if settings.get("lowpass_hz") is not None or settings.get("resample_sfreq_hz") is not None:
        raise ValueError("Final Stage 4 preserves native sampling and uses no extra low-pass.")
    if settings.get("notch_method", "spectrum_fit") != "spectrum_fit":
        raise ValueError("Stage 4 line removal uses spectrum_fit only.")
    summary = final_filter_summary(raw, settings, source_filename=source_filename)
    margin = settings.get("filter_edge_margin_sec")
    if summary["filter"] is not None and (margin is None or not np.isfinite(float(margin))
            or float(margin) < summary["filter"]["half_support_sec"]):
        raise RuntimeError("Review final-filter edge exposure first; margin must cover the actual half-support.")
    stretches = good_stretches(raw)
    picks = mne.pick_types(raw.info, meg=True, exclude=[])
    if summary["filter"] is not None:
        raw.filter(settings["highpass_hz"], None, picks=picks,
                   skip_by_annotation=("edge", "bad"), verbose="WARNING")
    # This operation is independent of the high-pass and its historical P38
    # exception. Never return early merely because there is no FIR to apply.
    raw.notch_filter(settings.get("notch_freqs_hz", [50., 100.]), picks=picks, method="spectrum_fit",
                     filter_length=settings.get("spectrum_fit_window", "10s"),
                     skip_by_annotation=("edge", "bad"), verbose="WARNING")
    if summary["filter"] is None:
        return summary  # No new FIR-exposure annotations or fabricated margin.
    # Padding exposure is supplementary QC. The name deliberately starts with
    # neither BAD nor EDGE, so MNE retains overlapping otherwise eligible epochs.
    onset, duration = [], []
    sfreq = raw.info["sfreq"]
    edge = int(np.ceil(float(margin) * sfreq))
    for start, stop in stretches:
        intervals = [(start, stop)] if stop - start <= 2 * edge else [(start, start + edge), (stop - edge, stop)]
        for a, b in intervals:
            onset.append(a / sfreq + raw.first_time)
            duration.append((b - a) / sfreq)
    edges = mne.Annotations(onset, duration, ["filter_edge_exposure"] * len(onset),
                            orig_time=raw.annotations.orig_time)
    # Append directly in acquisition coordinates; set_annotations would add
    # first_time again when meas_date is absent.
    raw.annotations.append(edges.onset, edges.duration, edges.description)
    return summary
