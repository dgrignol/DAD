#!/usr/bin/env python3
"""Audit existing ICA application and final epoch identities without accepting them.

Usage (existing MNE environment; run AFTER user-submitted application finishes)::

    python audit_ica_application.py --label NEW-LABEL --subjects 4 5 6 --jobs 3
    # Add --task-signal for matched native / ICA-only / final task-mean evidence.
    # Use --phase epochs after epoch/export; --phase inventory is read-only preflight.
    # --phase inputs validates accepted inline inputs and Stage 2/3 timing only:
    # DAD_EPOCH_PREPARATION=inline_no_ica DAD_PROCESSING_RUN=noica-NEW \
    #   python audit_ica_application.py --phase inputs --label inputs-NEW --jobs 4
    # Cluster: --site-config meg/params/pipeline_cluster.yaml --jobs ALLOCATED_CPUS.

Inputs: current guarded Stage 4 receipts, exact saved model/review/apply reports,
native annotated and applied FIFs, accepted diagnostic trial identities; epoch
mode additionally reads final FIF metadata and complete trial-decision CSVs.
Default cohort is P04-P14/P16-P39. A pending participant is reported explicitly;
select the actual authorized batch with --subjects when auditing a partial batch.

Outputs: immutable participant JSONs and a joined summary under ICA review updates.
Optional task evidence uses ALL annotation-surviving diagnostic trial identities,
the accepted -0.5/+2.667 s native bounds, configured detrend/baseline, and no sensor
amplitude rejection. Matched means compare native, reconstructed ICA-only and
saved final signals; they are diagnostic arrays, not production epochs or ERFs.
Differences/correlations are descriptive, with no invented acceptance threshold.
Native coordinate/annotation and full STI101/MISC preservation checks are reused.
Final spectrum-fit figures and short-stretch evidence remain human review items.

Independent participants run in separate bounded processes with one numerical
thread each; task-signal mode needs about 24 GB per worker (three native arrays).
Use an appropriate allocated --jobs. No signal, model, choice, annotation,
timing, scientific acceptance or earlier output is overwritten by this script.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import traceback

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from scipy.signal import detrend
from scipy.io import loadmat

from annotation_triage import subject_config
from epoch_input_preparation import epoch_stage, verify_epoch_artifacts, verify_export_artifacts
from continuous_cleaning import check_preservation, digest, file_hash, file_stamp, require_review
from participant_parallel import map_participants, participant_workers
from pipeline_config import WORKING_SUBJECT_IDS
from preprocessing_stage_utils import derivative_path_for_run, discover_raw_run_records, stage_output_dirs, utc_now
from record_provisional_ica_choices import require
from review_final_ica import validate_review

PROJECT = Path(__file__).resolve().parents[3]


# -----------------------------------------------------------------------------
# Matched trial signals reuse accepted sample identities and configured display
# operators. Numeric change is evidence to inspect, never proof of preservation.
# -----------------------------------------------------------------------------

def trial_signal_sums(raws, rows, cfg, picks):
    """Accumulate the same native samples in every view with no new trial rule."""
    sums, counts = {}, Counter()
    sfreq = raws[0].info["sfreq"]
    for row in rows:
        first = int(row["proposed_epoch_first_sample_raw"]) - raws[0].first_samp
        stop = int(row["proposed_epoch_last_sample_raw"]) - raws[0].first_samp + 1
        require(first >= 0 and stop <= raws[0].n_times, "Diagnostic trial crosses a native boundary.")
        condition = str(row["condition_code"])
        for view, raw in zip(("native", "ica_only", "final"), raws):
            values = raw.get_data(picks=picks, start=first, stop=stop)
            order = cfg["epoching"]["detrend"]
            if order is not None:
                values = detrend(values, axis=-1, type="linear" if order == 1 else "constant")
            baseline = cfg["epoching"]["baseline_sec"]
            if baseline is not None:
                times = np.arange(values.shape[-1]) / sfreq + cfg["epoching"]["tmin_sec"]
                left = times[0] if baseline[0] is None else baseline[0]
                right = times[-1] if baseline[1] is None else baseline[1]
                mask = (times >= left) & (times <= right)
                require(mask.any(), "Diagnostic baseline has no samples.")
                values -= values[:, mask].mean(axis=-1, keepdims=True)
            key = f"condition_{condition}_{view}"
            sums[key] = sums.get(key, 0) + values
        counts[condition] += 1
    return sums, counts


def save_task_evidence(sums, counts, names, kinds, cfg, output, sfreq):
    """Save matched condition means, sensor-type changes and a compact RMS plot."""
    means = {key: value / counts[key.split("_")[1]] for key, value in sums.items()}
    metrics = []
    fig, axes = plt.subplots(len(counts), 2, figsize=(13, 3 * len(counts)), squeeze=False)
    for row, condition in enumerate(sorted(counts)):
        for col, kind in enumerate(("mag", "grad")):
            selected = np.flatnonzero(np.asarray(kinds) == kind)
            if not len(selected):
                continue
            before = means[f"condition_{condition}_native"][selected]
            denominator = np.linalg.norm(before)
            for view in ("native", "ica_only", "final"):
                values = means[f"condition_{condition}_{view}"][selected]
                corr = (float(np.corrcoef(before.ravel(), values.ravel())[0, 1])
                        if np.std(before) and np.std(values) else None)
                metrics.append({"condition": condition, "sensor_type": kind, "view": view,
                    "trials": counts[condition], "correlation_with_native": corr,
                    "relative_change_from_native": float(np.linalg.norm(values-before) / denominator) if denominator else None})
                times = np.arange(values.shape[-1]) / sfreq + cfg["epoching"]["tmin_sec"]
                scale = 1e15 if kind == "mag" else 1e13
                axes[row, col].plot(times, np.sqrt(np.mean(values**2, axis=0)) * scale, label=view)
            axes[row, col].set(title=f"Condition {condition}; {kind}; n={counts[condition]}",
                xlabel="Seconds from accepted motion onset", ylabel="RMS condition mean (fT)" if kind == "mag" else "RMS condition mean (fT/cm)")
            axes[row, col].legend()
    fig.tight_layout(); fig.savefig(output / "matched_task_means.png", dpi=130); plt.close(fig)
    np.savez_compressed(output / "matched_task_means.npz", **means, channel_names=np.asarray(names), channel_types=np.asarray(kinds))
    return metrics


# -----------------------------------------------------------------------------
# Full application checks preserve source/model guards and all auxiliary samples.
# -----------------------------------------------------------------------------

def audit_application(args, subject, output):
    """Check exact-model application, filter contract and native preservation."""
    cfg, review, records, inputs, model_path, _ = validate_review(args, subject)
    dirs = stage_output_dirs(cfg)
    require_review(cfg, "03_C_ICA", dirs["reports"])
    applied = json.loads((dirs["reports"] / "03_C_ICA_apply_report.json").read_text())
    require(applied["status"] == "completed" and applied["processing_pass"] == "apply", "Application is incomplete.")
    require(applied["cache_status"] == "reused_compatible", "Application did not reuse the saved model.")
    require(applied["prefit_model_sha256"] == review["prefit_model_sha256"] == cfg["ica"]["reviewed_model_sha256"], "Applied model changed.")
    require(applied["manual_exclude"] == cfg["ica"]["manual_exclude"], "Applied exclusions disagree with receipt.")
    require(applied["fingerprint"] == review["fingerprint"], "Application source/fit identity changed.")
    for row in applied["evidence_files"]:
        require(file_stamp(Path(cfg["resolved_paths"]["derivatives_root"]) / row["relative_path"]) == row["stamp"], "Application input evidence changed.")
    native = json.loads((dirs["figures"] / "03_C_ICA_review_native_diagnostics.json").read_text())
    require(native["model_sha256"] == review["prefit_model_sha256"], "Task reference belongs to another model.")
    by_source = {row["source_filename"]: row for row in applied["runs"]}
    require(set(by_source) == {record.source_filename for record in records}, "Applied acquisition set differs.")
    model = mne.preprocessing.read_ica(model_path, verbose="ERROR") if args.task_signal else None
    runs, total_sums, total_counts = [], {}, Counter()
    channel_names, channel_kinds = None, None
    for record, input_path in zip(records, inputs):
        row = by_source[record.source_filename]
        final_path = derivative_path_for_run(record, dirs["sss"], "_ica.fif")
        source = mne.io.read_raw_fif(input_path, preload=args.task_signal, verbose="ERROR")
        final = mne.io.read_raw_fif(final_path, preload=args.task_signal, verbose="ERROR")
        expected = digest({"model": review["prefit_model_sha256"], "input": file_stamp(input_path),
            "exclude": applied["manual_exclude"], "final_filter": cfg["preprocessing"],
            "final_filter_contract": "independent_spectrum_fit_2026-10-28"})
        require(f"application_sha256={expected}" in str(final.info["description"]), "Saved application identity differs.")
        preservation = check_preservation(source, final, annotations=False)
        require(len(source.annotations) == len(final.annotations), "Unexpected annotation addition after spectrum-fit-only application.")
        for field in ("dev_head_t", "ctf_head_t", "dev_ctf_t"):
            before_t, after_t = source.info.get(field), final.info.get(field)
            require((before_t is None) == (after_t is None), f"Coordinate transform presence changed: {field}")
            if before_t is not None:
                require(before_t["from"] == after_t["from"] and before_t["to"] == after_t["to"] and
                        np.allclose(before_t["trans"], after_t["trans"], rtol=2e-7, atol=1e-8),
                        f"Coordinate transform changed: {field}")
        filtering = row["final_filter"]
        require(filtering["highpass_status"] == "not_selected" and filtering["filter"] is None, "Unexpected final high-pass.")
        require(filtering["spectrum_fit"]["freqs_hz"] == [50., 100.] and filtering["spectrum_fit"]["window"] == "10s", "Line-removal settings changed.")
        figures = [dirs["figures"] / Path(row[key]).name for key in ("before_after_ica_figure", "final_filter_psd_figure")]
        require(all(path.is_file() and path.stat().st_size > 0 for path in figures), "Matched ICA-only/final-filter figures missing.")
        if args.task_signal:
            corrected = source.copy()
            model.apply(corrected, exclude=applied["manual_exclude"], verbose="ERROR")
            picks = mne.pick_types(source.info, meg=True, exclude="bads")
            names = [source.ch_names[index] for index in picks]
            kinds = [source.get_channel_types()[index] for index in picks]
            require(channel_names is None or names == channel_names, "Task sensor basis changed across acquisitions.")
            channel_names, channel_kinds = names, kinds
            trials = [trial for trial in native["retained_trial_metadata"] if trial["source_fif"] == record.source_filename]
            sums, counts = trial_signal_sums((source, corrected, final), trials, cfg, picks)
            for key, value in sums.items(): total_sums[key] = total_sums.get(key, 0) + value
            total_counts.update(counts)
            del corrected
        runs.append({"source": record.source_filename, "preservation": preservation,
            "output_stamp": file_stamp(final_path), "spectrum_fit": filtering["spectrum_fit"],
            "good_stretch_durations_sec": filtering["good_stretch_durations_sec"],
            "figures": [{"path": str(path), "sha256": file_hash(path)} for path in figures],
            "p10_slow_contamination_qc": row.get("p10_slow_contamination_qc")})
        print(f"P{subject:02}: preserved {record.source_filename}", flush=True)
        del source, final
    metrics = None
    if args.task_signal:
        require(sum(total_counts.values()) == len(native["retained_trial_metadata"]), "Matched task identities/counts changed.")
        metrics = save_task_evidence(total_sums, total_counts, channel_names, channel_kinds, cfg, output, native["retained_trial_metadata"][0]["sfreq_hz"])
    return {"model_sha256": review["prefit_model_sha256"], "exclusions": applied["manual_exclude"],
        "runs": runs, "task_signal_metrics": metrics, "task_counts": dict(total_counts),
        "review_needed": "Inspect ICA-only examples, final spectra/BAD boundaries/short stretches and matched task changes; P10 retention remains open."}


def audit_epochs(args, subject, output):
    """Reconcile epoch identities/drops and exported tensor axes/counts/NaN masks."""
    cfg = subject_config(args, subject); dirs = stage_output_dirs(cfg)
    verified = verify_epoch_artifacts(cfg)
    verify_export_artifacts(cfg)
    root = Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / "meg"
    stem = f"sub-{subject:02}_ses-01_task-Attention_run-all"
    epochs = mne.read_epochs(root / f"{stem}_epo.fif", preload=False, verbose="ERROR")
    require(epochs.metadata is not None and "trial_key" in epochs.metadata, "Final epoch identities are missing.")
    native = json.loads((dirs["figures"] / "03_C_ICA_review_native_diagnostics.json").read_text())
    expected = {row["trial_key"]: row for row in native["retained_trial_metadata"]}
    actual = epochs.metadata.to_dict("records")
    stored_metadata = pd.read_csv(root / f"{stem}_epoch_metadata.csv")
    require(stored_metadata["trial_key"].tolist() == epochs.metadata["trial_key"].tolist(), "FIF and CSV final trial orders differ.")
    require(np.isclose(epochs.tmin, -.5, atol=1e-12) and np.isclose(epochs.tmax, 2.667, atol=1e-12) and
            epochs.info["sfreq"] == float(native["retained_trial_metadata"][0]["sfreq_hz"]), "Actual native epoch grid changed.")
    keys = {row["trial_key"] for row in actual}
    require(len(keys) == len(actual) and keys.issubset(expected), "Final epoch identity added/duplicated outside accepted footprint.")
    for row in actual:
        for field in ("source_fif", "raw_sample_frame", "epoch_event_sample_raw", "condition_code", "trial_sequence_identity"):
            reference = expected[row["trial_key"]][field]
            equal = row[field] == reference if field == "source_fif" else float(row[field]) == float(reference)
            require(equal, f"Final trial field changed: {field}")
    decisions = pd.read_csv(root / f"{stem}_trial_decisions.csv").fillna("")
    dropped = decisions[decisions["trial_key"].isin(set(expected)-keys)]
    require(len(dropped) == len(set(expected)-keys) and all(dropped["mne_drop_reasons"].astype(str).str.len()), "Unexplained final trial loss.")
    report = json.loads((dirs["reports"] / "04_epoching_report.json").read_text())
    if cfg["epoching"].get("preparation_mode", "saved_ica") != "saved_ica":
        require(all(row["continuous_preparation_qc"]["preservation"]["auxiliary_samples_verified"]
            and not row["continuous_preparation_qc"]["continuous_file_saved"] for row in verified["per_run"]),
            "Inline native preservation evidence is missing.")
    require(report["counts"]["n_final_epochs"] == len(actual) and report["epoch_window_sec"] == [-.5, 2.667], "Final epoch count/window changed.")
    require(report["baseline_sec"] is None and report["detrend"] == 1 and report["reject_by_annotation"], "Epoch operators changed.")
    require(report["stage3_timing_gate"] is not None, "Accepted corrected-onset timing is absent.")
    cells = Counter((str(row["source_fif"]), str(row["condition_code"]), str(row["trial_sequence_identity"])) for row in actual)
    dropped.to_csv(output / "extra_trial_drops.csv", index=False)
    exports = audit_exports(cfg, dirs, epochs)
    return {"native_annotation_survivors": len(expected), "final_epochs": len(actual), "extra_mne_drops": len(dropped),
        "export_audit": exports,
        "condition_counts": dict(Counter(str(row["condition_code"]) for row in actual)),
        "run_condition_identity_counts": [{"source": key[0], "condition": key[1], "identity": key[2], "count": value} for key, value in sorted(cells.items())]}


def audit_exports(cfg, dirs, epochs):
    """Read actual exported MAT/JSONs on the fixed axis and reconcile final counts.

    No averaging or scientific selection is repeated: every cell count must
    match the final FIF's metadata, and stored tensor labels/times/NaN channels
    must agree with the repository axis and condition sidecars.
    """
    path = dirs["reports"] / "06_export_neural_for_drsa_report.json"
    report = json.loads(path.read_text())
    require(report["status"] == "completed", "Neural export is incomplete.")
    channel_path = Path(cfg["resolved_paths"]["project_root"]) / cfg["export"]["canonical_meg_channels_file"]
    canonical = [line.strip() for line in channel_path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
    require(report["channel_axis"]["ordered_channel_labels"] == canonical, "Export canonical channel order changed.")
    export_root = Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / "drsa"
    result = []
    require(len(report["conditions"]) == len(cfg["epoching"]["event_id_map"]), "Export condition set changed.")
    for condition in report["conditions"]:
        name = condition["condition_name"]
        rows = epochs.metadata[epochs.metadata["condition_name"] == name]
        counts = Counter(int(float(value)) for value in rows["trial_sequence_identity"])
        saved_counts = {int(key): int(value) for key, value in condition["counts_by_identity"].items()}
        require(dict(counts) == saved_counts and all(saved_counts.values()), "Export identity counts differ from final epochs.")
        mat_path = export_root / Path(condition["mat_path"]).name
        sidecar = json.loads((export_root / Path(condition["json_path"]).name).read_text())
        require(sidecar == condition, "Condition sidecar and aggregate export disagree.")
        mat = loadmat(mat_path, variable_names=["neural", "channel_labels", "identity_values", "time_axis_sec", "sampling_rate_hz"], squeeze_me=True)
        tensor = mat["neural"]
        require(list(tensor.shape) == condition["tensor_shape"] and
                [str(value) for value in np.asarray(mat["channel_labels"]).ravel()] == canonical, "MAT tensor shape/channel labels changed.")
        first, last = cfg["subject_profile"]["identity_range"]
        require(np.array_equal(np.asarray(mat["identity_values"]).ravel(), np.arange(first, last+1)), "MAT identity axis changed.")
        require(float(mat["sampling_rate_hz"]) == float(cfg["export"]["sampling_rate_hz"]) and
                np.allclose(np.asarray(mat["time_axis_sec"]).ravel(), condition["time_axis_sec"], rtol=0, atol=1e-12), "MAT export time axis/rate changed.")
        nan_channels = set(condition["nan_channel_labels"])
        for index, channel in enumerate(canonical):
            require(np.isnan(tensor[:, index]).all() if channel in nan_channels else np.isfinite(tensor[:, index]).all(),
                    f"Export NaN channel policy violated: {channel}")
        result.append({"condition": name, "tensor_shape": list(tensor.shape), "counts_by_identity": saved_counts,
            "nan_channel_labels": sorted(nan_channels), "mat_stamp": file_stamp(mat_path)})
    return result


# -----------------------------------------------------------------------------
# Participant outputs are disjoint; shared publication occurs only after join.
# -----------------------------------------------------------------------------

def worker(task):
    args, subject, output = task
    directory = output / f"P{subject:02}"; directory.mkdir()
    print(f"P{subject:02}: starting {args.phase}", flush=True)
    try:
        cfg = subject_config(args, subject); dirs = stage_output_dirs(cfg)
        if args.phase == "inventory":
            records = discover_raw_run_records(cfg)
            result = {"application_report_exists": (dirs["reports"] / "03_C_ICA_apply_report.json").exists(),
                "native_runs": len(records), "cleaned_runs_present": sum(derivative_path_for_run(record, dirs["sss"], "_ica.fif").exists() for record in records)}
        elif args.phase == "inputs":
            review = Path(cfg["resolved_paths"]["derivatives_root"]) / "visual_timing/s3-v002/reviews" / f"sub-{subject:02}_review.json"
            manifest = Path(cfg["resolved_paths"].get("upstream_derivatives_reports_dir",
                cfg["resolved_paths"]["derivatives_reports_dir"])) / "source_manifest.json"
            result = epoch_stage().run_epoching(cfg, source_manifest=manifest, overwrite=False,
                validate_only=True, timing_review=review)
        else:
            result = (audit_application if args.phase == "application" else audit_epochs)(args, subject, directory)
        result.update(subject=subject, ok=True)
    except Exception as error:
        result = {"subject": subject, "ok": False, "error": str(error), "traceback": traceback.format_exc()}
    (directory / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"P{subject:02}: {args.phase} {'checked' if result['ok'] else result['error']}", flush=True)
    return result


def main():
    """Expose resources, select the explicit cohort and preserve prior evidence."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=PROJECT / "meg/params/pipeline_default.yaml")
    p.add_argument("--site-config", type=Path)
    p.add_argument("--subject-config-dir", type=Path, default=PROJECT / "meg/params")
    p.add_argument("--subjects", nargs="+", type=int, default=list(map(int, WORKING_SUBJECT_IDS)))
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--label", required=True)
    p.add_argument("--phase", choices=("inventory", "inputs", "application", "epochs"), default="application")
    p.add_argument("--task-signal", action="store_true")
    args = p.parse_args(); args.expected_model = None
    require(len(set(args.subjects)) == len(args.subjects) and set(args.subjects).issubset(set(map(int, WORKING_SUBJECT_IDS))), "Invalid working-cohort selection.")
    require(args.label.replace("-", "").replace("_", "").isalnum(), "Use a new simple label.")
    cfg = subject_config(args, args.subjects[0])
    output = Path(cfg["resolved_paths"]["derivatives_root"]) / "reports/ica_review_updates" / args.label
    output.mkdir(exist_ok=False)
    jobs = participant_workers(args.jobs, len(args.subjects))
    results = map_participants(worker, [(args, subject, output) for subject in args.subjects], jobs)
    summary = {"generated_utc": utc_now(), "phase": args.phase, "participants": results,
        "all_checks_passed": all(row["ok"] for row in results), "resources": {"participant_processes": jobs, "numerical_threads_per_worker": 1},
        "scientific_acceptance": False, "cleaned_output_acceptance": False}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    require(summary["all_checks_passed"], "One or more audits failed; inspect participant evidence.")
    print("JOINED", output, flush=True)


if __name__ == "__main__": main()
