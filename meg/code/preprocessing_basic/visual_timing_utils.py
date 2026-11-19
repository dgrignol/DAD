#!/usr/bin/env python3
"""Stage 3 visual-timing measurement and immutable epoching handoff helpers.

Usage::

    from visual_timing_utils import measure_edge, DEFAULT_METHOD
    result = measure_edge(window, 1000., 75, 1, 4, 120., DEFAULT_METHOD)

Inputs are short, unfiltered diode windows on the native FIF sample grid.
``anchor_index`` is relative to that window; all returned offsets are relative
to the original STI101 anchor. Polarity is a recording/segment decision, never
optimized separately for each event. The candidate 50% edge is not the exact
moving-dot luminance onset. QC thresholds are proposals, not accepted science.

Pure measurement -> QC -> restricted donor estimates -> timing overlay. The
reviewed-map loader validates a separately authored review record and original
Stage 2 keys before adding ``epoch_event_sample_raw`` to copies of ledger rows.
No function changes Stage 1/2 evidence or writes a raw FIF.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from downstream_stage_utils import csv_bool
from trial_preparation import file_hash, read_csv, stage2_paths


SCHEMA_VERSION = "1.0"
VISUAL_CODES = {71, 81, 91, 111, 112, 114, 115, 116, 117, 118, 153}
CONDITION_CODES = {71, 81, 91}
DEFAULT_METHOD = {
    "method_id": "fractional_first_response_candidate_v1",
    "window_ms": [-75.0, 175.0],
    "baseline_ms": [-50.0, -5.0],
    "search_ms": [-20.0, 80.0],
    "post_warmup_trough_ms": [0.0, 30.0],
    "post_warmup_max_gap_ms": 150.0,
    "threshold_fractions": [0.10, 0.25, 0.50, 0.75, 0.90],
    "primary_fraction": 0.50,
    "first_response_peak_ms": 8.0,
    "sustained_ms": 2.0,
    "minimum_amplitude_v": 0.005,
    "minimum_snr": 8.0,
    "maximum_sensitivity_ms": 3.0,
    "pulse_tolerance_ms": 4.0,
    "minimum_donor_events_per_code": 10,
    "bootstrap_resamples": 2000,
    "random_seed": 20261020,
    "estimator": "condition_specific_block_balanced_median",
    "rounding": "nearest_integer_half_away_from_zero",
    "approval": "proposal_only",
}


# -----------------------------------------------------------------------------
# Native-grid edge measurements: polarity, low/high transitions and explicit QC
# -----------------------------------------------------------------------------


def round_half_away(value: float) -> int:
    """Convert a finite signed sample shift without Python's ties-to-even rule."""
    if not math.isfinite(value):
        raise ValueError("A sample shift must be finite.")
    return int(math.copysign(math.floor(abs(value) + 0.5), value))


def event_family(code: int) -> str:
    """Keep scheduled, unscheduled and diagnostic paths distinct."""
    return ("condition_onset" if code in CONDITION_CODES else
            "question_onset" if code == 118 else
            "warmup" if code == 153 else
            "mid_trial" if code in VISUAL_CODES else "not_optically_marked")


def analysis_eligible(row: Mapping) -> bool:
    """Select using a verified final overlay while retaining original policy.

    Legacy rows use included_by_trial_policy. Revised rows retain that field
    and add final_analysis_eligible; load_reviewed_timing verifies the complete
    hash-bound application and explicit human review before epoching uses it.
    """
    return csv_bool(row["final_analysis_eligible"] if "final_analysis_eligible" in row else row.get("included_by_trial_policy",False))


def sample_slice(ms, sfreq, anchor, size):
    """Return a clipped native-grid interval plus whether it was truncated."""
    start, stop = [anchor + round_half_away(x * sfreq / 1000) for x in ms]
    return slice(max(0, start), min(size, stop)), start < 0 or stop > size


def infer_polarity(windows, sfreq, method=DEFAULT_METHOD) -> tuple[int | None, dict]:
    """Infer one sign from ordinary-event controls; abstain on weak evidence.

    Each item is (signal, anchor_index). Warm-up must not be a polarity control
    because the previous pulse can still be high at the digital anchor.
    """
    scores = []
    for values, anchor in windows:
        y = np.asarray(values, float)
        pre, clipped = sample_slice(method["baseline_ms"], sfreq, anchor, len(y))
        post, clipped_post = sample_slice([0, 60], sfreq, anchor, len(y))
        if clipped or clipped_post or not np.isfinite(y).all():
            continue
        baseline = float(np.median(y[pre]))
        scores.append((max(0., float(np.quantile(y[post], .90)) - baseline),
                       max(0., baseline - float(np.quantile(y[post], .10)))))
    if not scores:
        return None, {"reason": "no_polarity_controls"}
    positive, negative = np.median(scores, axis=0)
    dominant, other = max(positive, negative), min(positive, negative)
    valid = dominant >= method["minimum_amplitude_v"] and dominant > 2 * other
    return ((1 if positive > negative else -1) if valid else None,
            {"positive_amplitude_v": float(positive), "negative_amplitude_v": float(negative),
             "n_controls": len(scores), "reason": "recording_controls" if valid else "ambiguous_polarity"})


def measure_edge(values, sfreq: float, anchor_index: int, polarity: int,
                 expected_frames: int, refresh_hz: float, method=DEFAULT_METHOD,
                 *, warmup: bool = False, post_warmup: bool = False) -> dict[str, Any]:
    """Find the first noise-verified response and fractional edge candidates.

    A provisional sustained transition finds the response; fractional crossings
    need not themselves persist, avoiding a jump to later plateau modulation.
    Warm-up uses a post-anchor low trough and requires a subsequent rise. An
    ordinary contaminated baseline is flagged rather than silently repaired.
    The documented post-warm-up exception uses a low trough only when the
    accepted preceding marker establishes a nearby warm-up context. Its search
    starts at the digital anchor to exclude the preceding warm-up light rise.
    Signed edges are retained. Missing/ambiguous responses are not zero lags.
    """
    if sfreq <= 0 or refresh_hz <= 0 or polarity not in (-1, 1):
        raise ValueError("Positive sampling/refresh rates and fixed polarity are required.")
    y = np.asarray(values, float) * polarity
    result = {"measurement_status": "unusable", "qc_reasons": "", "edge_offset_samples": None,
              "lag_ms": None, "polarity": polarity, "baseline_method": "pre_edge_median"}
    if y.ndim != 1 or len(y) < 4 or not np.isfinite(y).all():
        return dict(result, qc_reasons="nonfinite_or_short_signal")
    pre, clipped = sample_slice(method["baseline_ms"], sfreq, anchor_index, len(y))
    search, clipped_search = sample_slice([0, 80] if warmup or post_warmup else method["search_ms"],
                                          sfreq, anchor_index, len(y))
    if clipped or clipped_search or pre.stop - pre.start < 2:
        return dict(result, qc_reasons="truncated_detection_window")
    baseline = float(np.median(y[pre]))
    noise = float(1.4826 * np.median(np.abs(y[pre] - baseline)))
    if warmup or post_warmup:
        trough, _ = sample_slice(method["post_warmup_trough_ms"], sfreq, anchor_index, len(y))
        baseline = float(np.quantile(y[trough], .10))
        result["baseline_method"] = "warmup_post_anchor_low_trough" if warmup else "post_warmup_low_trough"
        result["qc_notes"] = "cadence_context_low_trough;previous_high_is_not_baseline"
        # A previous high plateau is not baseline noise for cadence pulses.
        low = y[trough][y[trough] <= baseline + method["minimum_amplitude_v"] / 2]
        noise = float(1.4826 * np.median(np.abs(low - np.median(low)))) if len(low) else noise
    candidate_peak = float(np.quantile(y[search], .95))
    amplitude = candidate_peak - baseline
    result.update(baseline_normalized_v=baseline, baseline_noise_v=noise,
                  amplitude_v=amplitude, rounded_levels=len(np.unique(np.round(y, 8))),
                  dynamic_range_v=float(np.ptp(y)))
    if amplitude < method["minimum_amplitude_v"] or np.ptp(y) < method["minimum_amplitude_v"]:
        return dict(result, qc_reasons="flat_or_low_amplitude")
    sustain = max(1, round_half_away(method["sustained_ms"] * sfreq / 1000))
    provisional = baseline + .5 * amplitude
    hits = []
    for i in range(search.start, search.stop - sustain + 1):
        # Require an observed low-to-high transition, including for warm-up.
        if i and y[i - 1] < provisional and np.all(y[i:i + sustain] >= provisional):
            hits.append(i)
    if not hits:
        return dict(result, qc_reasons="no_verified_rise_or_contaminated_baseline")
    first = hits[0]
    peak_stop = min(len(y), first + max(1, round_half_away(method["first_response_peak_ms"] * sfreq / 1000)) + 1)
    amplitude = float(np.max(y[first:peak_stop]) - baseline)
    local_start = first
    while local_start > search.start and y[local_start - 1] > baseline + .05 * amplitude:
        local_start -= 1
    crossings = {}
    for fraction in method["threshold_fractions"]:
        threshold = baseline + fraction * amplitude
        found = np.flatnonzero(y[local_start:peak_stop] >= threshold)
        crossings[f"edge_{round(fraction * 100):02d}_offset_samples"] = (
            int(local_start + found[0] - anchor_index) if len(found) else None)
    primary = crossings["edge_50_offset_samples"]
    if primary is None or any(value is None for value in crossings.values()):
        return dict(result, **crossings, qc_reasons="missing_fractional_crossing")
    edge = anchor_index + primary
    high_threshold = baseline + .5 * amplitude
    falls = [i for i in range(edge + 1, len(y) - sustain + 1)
             if np.all(y[i:i + sustain] < high_threshold)]
    width = (falls[0] - edge) * 1000 / sfreq if falls else None
    span = (max(crossings.values()) - min(crossings.values())) * 1000 / sfreq
    reasons = []
    snr = amplitude / max(noise, np.finfo(float).eps)
    if amplitude < method["minimum_amplitude_v"]:
        reasons.append("low_first_response_amplitude")
    if snr < method["minimum_snr"]:
        reasons.append("noisy_baseline")
    if not (warmup or post_warmup) and np.std(y[pre]) > .20 * amplitude:
        reasons.append("baseline_contaminated")
    if width is None:
        reasons.append("pulse_fall_unavailable")
    elif abs(width - expected_frames * 1000 / refresh_hz) > method["pulse_tolerance_ms"]:
        reasons.append("pulse_duration")
    if span > method["maximum_sensitivity_ms"]:
        reasons.append("threshold_sensitive")
    if len(hits) > 1 and not warmup:
        reasons.append("multiple_candidate_rises")
    result.update(crossings, amplitude_v=amplitude, snr=snr,
                  edge_offset_samples=int(primary), lag_ms=primary * 1000 / sfreq,
                  threshold_span_ms=span, pulse_duration_ms=width,
                  expected_pulse_duration_ms=expected_frames * 1000 / refresh_hz,
                  candidate_rise_count=len(hits),
                  measurement_status="measured" if not reasons else "flagged",
                  qc_reasons=";".join(reasons),
                  sampling_interval_ms=1000 / sfreq,
                  uncertainty_kind="threshold_span_and_native_grid;patch_to_dot_unknown")
    return result


# -----------------------------------------------------------------------------
# Donor restrictions, robust estimates and separately identified uncertainty
# -----------------------------------------------------------------------------


def donor_allowed(subject: str, target: Mapping, donor: Mapping) -> bool:
    """Enforce decided P11/P24 directions; never borrow another participant."""
    if donor.get("subject_bids_id") != f"sub-{subject}" or donor.get("measurement_status") != "measured":
        return False
    if int(donor["trigger_code"]) != int(target["condition_code"]):
        return False
    if subject == "11" and int(target["behavioral_block"]) in (1, 2):
        return int(donor.get("behavioral_block") or 0) >= 3
    if subject == "24" and int(target["acquisition_run"]) in (1, 2, 3, 4):
        return int(donor["acquisition_run"]) == 5
    return False


def estimate_lag(donors: list[dict], method=DEFAULT_METHOD) -> dict:
    """Estimate a condition's block-balanced location and diagnostic intervals.

    Bootstrap samples blocks, then behavioral-run clusters, then events. One
    block supports only within-block resampling. Percentiles describe observed
    donor variation, not a fully calibrated target prediction interval.
    Unknown between-context/spatial transfer uncertainty remains explicit.
    """
    if len(donors) < method["minimum_donor_events_per_code"]:
        return {"status": "insufficient_donors", "n_donors": len(donors)}
    blocks = defaultdict(lambda: defaultdict(list))
    for row in donors:
        blocks[row["source_fif"]][str(row.get("behavioral_run", ""))].append(float(row["lag_ms"]))
    arrays = {b: {r: np.asarray(v) for r, v in runs.items()} for b, runs in blocks.items()}
    block_values = {b: np.concatenate(list(runs.values())) for b, runs in arrays.items()}
    block_medians = [float(np.median(v)) for v in block_values.values()]
    location = float(np.median(block_medians))
    rng = np.random.default_rng(method["random_seed"])
    bootstrap = []
    names = list(arrays)
    for _ in range(method["bootstrap_resamples"]):
        centers = []
        for b in rng.choice(names, len(names), replace=True):
            run_names = list(arrays[b])
            values = [rng.choice(arrays[b][r], len(arrays[b][r]), replace=True)
                      for r in rng.choice(run_names, len(run_names), replace=True)]
            centers.append(np.median(np.concatenate(values)))
        bootstrap.append(np.median(centers))
    all_values = np.asarray([float(r["lag_ms"]) for r in donors])
    return {"status": "estimated", "lag_estimate_ms": location, "n_donors": len(donors),
            "n_donor_blocks": len(arrays), "estimator": method["estimator"],
            "event_weighted_median_ms": float(np.median(all_values)),
            "block_medians_ms": block_medians,
            "location_ci95_ms": np.quantile(bootstrap, [.025, .975]).tolist(),
            "donor_observed_p05_p95_ms": np.quantile(all_values, [.05, .95]).tolist(),
            "donor_sd_ms": float(np.std(all_values, ddof=1)),
            "transfer_uncertainty": "unknown;not_included_in_intervals",
            "donor_event_keys": [r["event_key"] for r in donors],
            "held_out_block_errors_ms": [float(np.median(v) - np.median(
                [np.median(x) for key, x in block_values.items() if key != b]))
                for b, v in block_values.items()] if len(arrays) > 1 else [],
            "held_out_behavioral_run_errors_ms": [float(np.median(v) - np.median(np.concatenate(
                [other for block, runs in arrays.items() for run, other in runs.items() if (block, run) != (b, r)])))
                for b, runs in arrays.items() for r, v in runs.items()]
                if sum(len(runs) for runs in arrays.values()) > 1 else [],
            "validation_scope": "donor_block_holdout" if len(arrays) > 1 else "within_one_acquisition_only"}


def window_coverage(row: Mapping, proposed: int, tmin: float, tmax: float,
                    uncertainty_ms: tuple[float, float] | None = None) -> dict:
    """Check inclusive epoch endpoints against exclusive native recording end.

    Uncertainty describes only the supplied component, not unknown transfer
    error. Behavioral-boundary overlap is reported without changing eligibility.
    """
    sf = float(row["sfreq_hz"])
    first, end = int(row["recording_first_sample_raw"]), int(row["recording_end_sample_raw_exclusive"])
    a, b = proposed + round(tmin * sf), proposed + round(tmax * sf)
    result = {"proposed_epoch_first_sample_raw": a, "proposed_epoch_last_sample_raw": b,
              "proposed_window_inside_recording": first <= a and b < end,
              "epoch_overlaps_trial_interval_end": b >= int(row["trial_interval_end_raw_exclusive"]),
              "uncertainty_window_inside_recording": None}
    if uncertainty_ms is not None:
        lo, hi = uncertainty_ms
        digital = int(row["raw_sample_frame"])
        earliest = digital + math.floor(lo * sf / 1000) + round(tmin * sf)
        latest = digital + math.ceil(hi * sf / 1000) + round(tmax * sf)
        result["uncertainty_window_inside_recording"] = first <= earliest and latest < end
    return result


# -----------------------------------------------------------------------------
# Immutable bundle validation and the reviewed Stage 3 -> Stage 5 map interface
# -----------------------------------------------------------------------------


def verify_bundle(folder: Path) -> dict:
    """Verify all declared files; reject incomplete, unsafe or edited bundles."""
    folder = Path(folder).resolve()
    report = json.loads((folder / "manifest.json").read_text())
    required = {"event_measurements.csv", "event_timing_proposals.csv", "trial_timing_map.csv", "recording_summary.csv",
                "estimates.json", "contexts.json", "report.md", "diagnostics.png"}
    if (report.get("schema_version") != SCHEMA_VERSION or report.get("status") != "completed"
            or not required.issubset(report.get("output_sha256", {}))):
        raise RuntimeError("Incomplete or unsupported Stage 3 bundle.")
    for name, expected in report["output_sha256"].items():
        if Path(name).name != name or not (folder / name).is_file() or file_hash(folder / name) != expected:
            raise RuntimeError(f"Stage 3 output fingerprint mismatch: {name}")
    return report


def validate_cohort_donor_provenance(cfg, bundle_folder, bundle):
    """Recompute the three pooled medians from hash-bound original measured rows.

    The donor CSV is a convenience artifact, not its own authority. Every source
    manifest and event table is rechecked against s3-v001 and donor identities,
    QC status, native arithmetic, explicit exclusions and median are recomputed.
    """
    import hashlib
    folder=Path(bundle_folder)
    donor_path=folder/"cohort_pooled_donor_measurements.csv"
    provenance_path=folder/"cohort_pooled_donor_provenance.json"
    expected=bundle.get("onset_estimation_method",{})
    if not donor_path.is_file() or not provenance_path.is_file():
        raise RuntimeError("Pooled assignment lacks its donor provenance artifacts.")
    if file_hash(donor_path)!=expected.get("donor_measurements_sha256") or file_hash(provenance_path)!=expected.get("donor_provenance_sha256"):
        raise RuntimeError("Pooled donor evidence is not bound to the final Stage 3 manifest.")
    provenance=json.loads(provenance_path.read_text());donors=read_csv(donor_path)
    if provenance.get("method")!="ordinary_median_of_unique_qc_usable_measured_cohort_events_v1" or provenance.get("analysis_eligibility_filter") is not False:
        raise RuntimeError("Unsupported pooled donor method or eligibility filtering.")
    excluded=set(provenance.get("donor_excluded_event_keys",[]));seen=set();by_code=defaultdict(list);by_subject=defaultdict(set)
    for row in donors:
        key=row["event_key"];code=int(row["trigger_code"])
        if code not in (71,81,91) or key in seen or key in excluded or row.get("qc_reasons",""):
            raise RuntimeError("Pooled donor CSV has an invalid/duplicate/excluded event.")
        seen.add(key);sf=float(row["sfreq_hz"])
        lag=(int(row["diode_edge_sample_raw"])-int(row["digital_sample_raw"]))*1000/sf
        if not math.isclose(lag,float(row["lag_ms"]),rel_tol=0,abs_tol=1e-8):raise RuntimeError("Pooled donor lag arithmetic mismatch.")
        by_code[code].append(lag);by_subject[row["subject_bids_id"]].add(key)
    source_root=Path(cfg["resolved_paths"]["derivatives_root"])/"visual_timing"/"s3-v001"
    for sub,declared in provenance.get("donor_sources",{}).items():
        source=source_root/sub/"ses-01";manifest_path=source/"manifest.json"
        if file_hash(manifest_path)!=declared["source_manifest_sha256"]:raise RuntimeError("Pooled donor source manifest changed.")
        manifest=json.loads(manifest_path.read_text());event_path=source/"event_measurements.csv"
        if file_hash(event_path)!=declared["event_measurements_sha256"] or manifest["output_sha256"]["event_measurements.csv"]!=declared["event_measurements_sha256"]:
            raise RuntimeError("Pooled donor source measurements changed.")
        source_rows=read_csv(event_path)
        expected_keys={r["event_key"] for r in source_rows if int(r["trigger_code"]) in (71,81,91)
            and r["measurement_status"]=="measured" and not r.get("qc_reasons","").strip() and r["event_key"] not in excluded}
        if expected_keys!=by_subject.get(sub,set()):raise RuntimeError("Pooled donor membership differs from complete measured source rows.")
    if not set(by_subject).issubset(provenance.get("donor_sources",{})):raise RuntimeError("Pooled donor source inventory is incomplete.")
    result={}
    for code,expected_n in ((71,7341),(81,7422),(91,7357)):
        values=np.asarray(by_code[code],float)
        if len(values)!=expected_n:raise RuntimeError(f"Pooled trigger {code} donor count mismatch.")
        result[code]=float(np.median(values))
        if result[code]!=float(expected["summary"][str(code)]["median_ms"]):raise RuntimeError("Pooled median differs from manifest summary.")
    return result


def load_reviewed_timing(cfg: Mapping, prepared: list[dict], review_path: str | Path,
                         application_path: str | Path | None = None) -> tuple[list[dict], dict]:
    """Load a human-authored accepted review, join original keys, validate shifts.

    Review JSON fields: schema_version, status='accepted', reviewer,
    reviewed_utc, decision_reference, bundle_manifest (relative path supported),
    bundle_manifest_sha256 and accepted_assignment_bases. Revised eligibility
    additionally requires the exact application hash and accepted review scopes.
    Pilot/staging bundles cannot be accepted. Original Stage 2 eligibility and
    coordinates remain separate from final_analysis_eligible. This function
    never edits an acceptance record or ledger.
    """
    review_path = Path(review_path).resolve()
    review = json.loads(review_path.read_text())
    if (review.get("schema_version") != SCHEMA_VERSION or review.get("status") != "accepted"
            or any(not review.get(k) for k in ("reviewer", "reviewed_utc", "decision_reference"))):
        raise RuntimeError("Stage 3 requires an explicit completed independent review record.")
    manifest = Path(review["bundle_manifest"])
    manifest = manifest if manifest.is_absolute() else review_path.parent / manifest
    if manifest.name != "manifest.json" or file_hash(manifest) != review["bundle_manifest_sha256"]:
        raise RuntimeError("Reviewed Stage 3 manifest hash mismatch.")
    bundle = verify_bundle(manifest.parent)
    if bundle["mode"] != "cohort_proposal" or bundle.get("task1_staging") or bundle["subject_bids_id"] != cfg["layout"]["subject_bids_id"]:
        raise RuntimeError("Pilot or wrong-participant Stage 3 map cannot feed epoching.")
    if bundle["stage2_report_sha256"] != file_hash(stage2_paths(cfg)["report"]):
        raise RuntimeError("Stage 2 completion record changed after timing.")
    if bundle["epoch_window_sec"] != [cfg["epoching"]["tmin_sec"], cfg["epoching"]["tmax_sec"]]:
        raise RuntimeError("Reviewed timing-map epoch window differs; explicit contract review required.")
    selections,application_gate = None,None
    if bundle.get("review_application_sha256"):
        from visual_timing_application import load_applied_eligibility
        selections,application_gate = load_applied_eligibility(cfg,prepared,manifest.parent,bundle,review,application_path)
    elif application_path or review.get("accepted_eligibility_application_sha256"):
        raise RuntimeError("An eligibility application cannot be attached to an unbound timing bundle.")
    maps = read_csv(manifest.parent / "trial_timing_map.csv")
    measurements = read_csv(manifest.parent / "event_measurements.csv")
    measurements_by_key = {r["event_key"]: r for r in measurements}
    estimates = json.loads((manifest.parent / "estimates.json").read_text())
    if len(measurements_by_key) != len(measurements):
        raise RuntimeError("Duplicate Stage 3 measurement keys.")
    by_key = {r["trial_key"]: r for r in maps}
    if len(by_key) != len(maps) or set(by_key) != {r["trial_key"] for r in prepared}:
        raise RuntimeError("Stage 3 trial keys are incomplete or duplicated.")
    output, used = [], defaultdict(set)
    pooled_centers=None
    allowed = set(review.get("accepted_assignment_bases", []))
    if not allowed.issubset({"measured", "estimated", "digital_fallback"}):
        raise RuntimeError("Unknown accepted timing assignment basis.")
    protected = ("source_fif", "raw_sample_frame", "condition_code", "sfreq_hz",
                 "recording_first_sample_raw", "recording_end_sample_raw_exclusive")
    for original in prepared:
        timing = by_key[original["trial_key"]]
        if any(str(timing[k]) != str(original[k]) for k in protected):
            raise RuntimeError("Stage 3 changed an original trial coordinate or condition.")
        if csv_bool(timing["included_by_trial_policy"]) != original["included_by_trial_policy"]:
            raise RuntimeError("Stage 3 changed trial eligibility.")
        row = dict(original)
        if selections is not None:
            selected=selections[original["trial_key"]]
            if (csv_bool(timing.get("final_analysis_eligible"))!=selected["final_analysis_eligible"]
                    or timing.get("analysis_exclusion_reasons","")!=selected["analysis_exclusion_reasons"]):
                raise RuntimeError("Timing map and accepted final eligibility disagree.")
            row.update({key:selected[key] for key in ("final_analysis_eligible","analysis_exclusion_reasons",
                "original_included_by_trial_policy","original_exclusion_reasons","saved_catch_type_code",
                "newly_included","newly_excluded","sequence_exclusion")})
            row["eligibility_application_sha256"]=application_gate["application_sha256"]
        elif "final_analysis_eligible" in timing:
            raise RuntimeError("Final eligibility lacks an accepted hash-bound application.")
        basis = timing["assignment_basis"]
        if analysis_eligible(row):
            if basis not in allowed or timing["assignment_status"] != "proposed":
                raise RuntimeError(f"Unreviewed timing assignment: {original['trial_key']}")
            digital = int(original["raw_sample_frame"])
            proposed, shift = int(timing["proposed_sample_raw"]), int(timing["applied_shift_samples"])
            if proposed != digital + shift:
                raise RuntimeError("Timing-map signed-shift arithmetic mismatch.")
            if basis == "digital_fallback" and (shift or timing["lag_estimate_ms"]):
                raise RuntimeError("Digital fallback must keep physical lag unknown.")
            if basis != "digital_fallback":
                lag = float(timing["lag_estimate_ms"])
                if shift != round_half_away(lag * float(original["sfreq_hz"]) / 1000):
                    raise RuntimeError("Timing-map rounding mismatch.")
            if basis == "measured":
                measurement = measurements_by_key.get(timing.get("measurement_event_key"), {})
                if (measurement.get("measurement_status") != "measured"
                        or measurement.get("source_fif") != original["source_fif"]
                        or int(measurement.get("digital_sample_raw", -1)) != digital
                        or int(measurement.get("trigger_code", -1)) != int(original["condition_code"])
                        or int(measurement.get("diode_edge_sample_raw", -1)) != proposed):
                    raise RuntimeError("Measured assignment lacks its exact usable observed edge.")
            if basis == "estimated":
                estimate = estimates.get(timing.get("donor_estimate_id"), {})
                donor_keys = estimate.get("donor_event_keys", [])
                if str(timing.get("donor_estimate_id","")).startswith("COHORT_code"):
                    if pooled_centers is None:pooled_centers=validate_cohort_donor_provenance(cfg,manifest.parent,bundle)
                    code=int(str(timing["donor_estimate_id"]).split("code")[-1])
                    if (code!=int(original["condition_code"]) or code not in pooled_centers
                            or estimate.get("status")!="estimated" or float(estimate["lag_estimate_ms"])!=pooled_centers[code]
                            or float(timing["lag_estimate_ms"])!=pooled_centers[code]
                            or timing.get("uncertainty_kind")!="observed_cohort_dispersion;target_transfer_unknown"):
                        raise RuntimeError("Cross-participant estimated assignment violates its measured donor provenance.")
                else:
                    donors = [measurements_by_key.get(key, {}) for key in donor_keys]
                    if (not donors or len(set(donor_keys)) != len(donors)
                            or not all(donor_allowed(cfg["subject"]["id"], original, donor) for donor in donors)
                            or bool(set(donor_keys).intersection(application_gate["donor_excluded_event_keys"] if application_gate else []))
                            or estimate.get("status") != "estimated"
                            or float(estimate["lag_estimate_ms"]) != float(timing["lag_estimate_ms"])):
                        raise RuntimeError("Estimated assignment violates its decided donor/provenance restriction.")
            if not int(original["recording_first_sample_raw"]) <= proposed < int(original["recording_end_sample_raw_exclusive"]):
                raise RuntimeError("Proposed epoch onset leaves its original source recording.")
            coverage = window_coverage(original, proposed, *bundle["epoch_window_sec"])
            # Coverage changes require later rejection accounting, not ledger edits.
            row.update(coverage)
            if proposed in used[original["source_fif"]]:
                raise RuntimeError("Timing-map shifted event collision.")
            used[original["source_fif"]].add(proposed)
            row["epoch_event_sample_raw"] = proposed
        row.update(timing_assignment_basis=basis, timing_review_sha256=file_hash(review_path),
                   timing_bundle_manifest_sha256=review["bundle_manifest_sha256"],
                   timing_lag_ms=timing["lag_estimate_ms"], timing_uncertainty=timing["uncertainty_kind"])
        output.append(row)
    for source in used:
        included = sorted((r for r in output if r["source_fif"] == source and analysis_eligible(r)),
                          key=lambda r: int(r["raw_sample_frame"]))
        if any(a["epoch_event_sample_raw"] >= b["epoch_event_sample_raw"] for a, b in zip(included, included[1:])):
            raise RuntimeError("Timing map reverses the original event order.")
    return output, {"review_path": str(review_path), "review_sha256": file_hash(review_path),
                    "bundle_manifest_sha256": review["bundle_manifest_sha256"], "version": bundle["version"],
                    **({"eligibility_application":application_gate} if application_gate else {})}
