#!/usr/bin/env python3
"""Guarded continuous preparation for epoching without permanent cleaned FIFs.

Usage through the ordinary epoch/export wrapper or stage CLI::

    DAD_EPOCH_PREPARATION=inline_no_ica DAD_PROCESSING_RUN=noica-v1 \
        python 04_epoching.py --config ... --subject-config ... \
        --tuning-config ... --source-manifest ORIGINAL/source_manifest.json \
        --timing-review ACCEPTED/sub-XX_review.json

Modes: saved_ica preserves the existing completed-file route; inline_no_ica reads
accepted annotated native inputs and omits ICA; inline_ica first applies the
exact source/model-bound Stage 4C choices. Both inline modes then reuse the
unchanged spectrum-fit 50/100 Hz, 10 s, BAD/EDGE-isolated continuous operation
before epoching. No additional HP/LP/resampling, fitting or model save exists.

Spectrum-fit bypasses isolated segments of at most eight native samples,
including empty segments. These cannot support the pinned MNE 1.10.2 / SciPy
1.13 multitaper defaults. Their samples and annotations remain untouched by
line removal; ordinary segments retain the same 50/100 Hz, 10 s operation.
Each bypass is recorded in per-acquisition QC. The 2026-11-17 cohort check found
233 such fragments (0.894 s total) with no overlap with any of the 19,988
annotation-surviving epoch windows across all 180 acquisitions.

Outputs: prepared native Raw in memory, preservation/QC records and small
figures in the selected processing run. No new continuous FIF is written.
An explicitly selected P24 auxiliary reconstruction may then replace/add its
six eye channels and omit its run-5 photodiode, after the native preservation
checks. That separate plan binds all reconstructed arrays and keeps every MEG
sample, annotation, trigger and acquisition setting unchanged.
Accepted upstream configuration, timing, annotations and canonical manifest
remain read-only. Each run has its own derived manifest, reports and analysis
outputs; source/model/filter/epoch/timing fingerprints guard resumed results.
Participant concurrency belongs to the array wrapper; acquisitions remain
ordered within a participant to bound memory. Use about 24 GB per participant.

This module intentionally imports the original continuous_cleaning functions
without editing that model-bound protected source file. Skipping ICA is an
explicit processing operation, not an accepted empty component review.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.annotations import _annotations_starts_stops

from continuous_cleaning import (annotation_rows, coordinates, check_preservation, digest,
    file_hash, file_stamp, filter_continuous, final_filter_summary, require_review, stage_configuration)
from downstream_stage_utils import ManualIcaApplyRun, resolve_manual_ica_apply_runs
from preprocessing_stage_utils import (derivative_path_for_run, discover_raw_run_records,
    guarded_path, stage_output_dirs, write_json)

# The pinned multitaper defaults request eight DPSS tapers with half-bandwidth
# four. Nine samples is the first supported size; one sample also breaks COLA.
MIN_SPECTRUM_FIT_SAMPLES = 9


# -----------------------------------------------------------------------------
# Upstream provenance stays canonical; output manifests belong to one run label.
# -----------------------------------------------------------------------------

def upstream_config(cfg):
    """Return the accepted upstream path view without changing scientific settings."""
    result = deepcopy(cfg)
    for key in ("derivatives_reports_dir", "derivatives_analysis_dir"):
        result["resolved_paths"][key] = cfg["resolved_paths"].get("upstream_" + key, cfg["resolved_paths"][key])
    return result


def processing_manifest(cfg, supplied, *, create=False):
    """Use a derived manifest for scoped runs, retaining immutable upstream ancestry."""
    if not cfg["resolved_paths"].get("processing_run_root"):
        return Path(supplied).resolve()
    source = Path(upstream_config(cfg)["resolved_paths"]["derivatives_reports_dir"]) / "source_manifest.json"
    target = Path(cfg["resolved_paths"]["derivatives_reports_dir"]) / "source_manifest.json"
    if Path(supplied).resolve() not in (source.resolve(), target.resolve()):
        raise RuntimeError("Processing run must derive from its canonical participant manifest.")
    stamp = file_hash(source)
    if target.exists():
        existing = json.loads(target.read_text()).get("processing_run", {})
        if existing != {"label": cfg["epoching"]["processing_run"], "upstream_manifest_sha256": stamp}:
            raise RuntimeError("Processing-run manifest ancestry changed; use a new run label.")
    elif create:
        payload = json.loads(source.read_text())
        # Inherit accepted upstream history only. No old downstream completion
        # can make a newly labelled epoch/export/analysis appear already complete.
        stages = payload.get("derivative_lineage", {}).get("stages", {})
        for name in ("epoching", "export_drsa", "matlab_drsa"):
            stages.pop(name, None)
        payload["processing_run"] = {"label": cfg["epoching"]["processing_run"], "upstream_manifest_sha256": stamp}
        write_json(target, payload, raw_root=cfg["resolved_paths"]["data_root"], overwrite=False)
    return target


def portable(value):
    """Remove mount prefixes from runtime evidence, retaining file/identity fields."""
    if isinstance(value, dict):
        return {key: portable(item) for key, item in value.items() if not key.endswith("_path")}
    if isinstance(value, (list, tuple)):
        return [portable(item) for item in value]
    return value


def native_input_evidence(raw, source_filename):
    """Bind every native split part, annotation and acquisition coordinate."""
    return {"source_filename": source_filename, "parts": [file_stamp(path) for path in raw.filenames],
        "coordinates": coordinates(raw), "annotations_sha256": digest(annotation_rows(raw))}


def validate_recorded_exclusions(cfg, report, exclusions):
    """Keep a recorded authorization bound to its exact model and chosen list.

    Ordinary manual reviews bind the list through reviewed_model_sha256. The
    frozen batch authorization additionally binds its content, deferred choices
    and participant row; editing tuning alone cannot change that authorization.
    Review provenance remains in ICA documentation/receipts, not new labels.
    """
    receipt = cfg.get("stage4_review", {}).get("03_C_ICA", {})
    if receipt.get("acceptance_kind") != "user_provisional_component_processing_authorization":
        return
    root = Path(cfg["resolved_paths"]["derivatives_root"]).resolve()
    authorization = (root / receipt["authorization_relative_path"]).resolve()
    authorization.relative_to(root)
    if file_hash(authorization) != receipt["authorization_sha256"]:
        raise RuntimeError("Recorded ICA authorization changed.")
    rows = [row for row in json.loads(authorization.read_text())["participants"]
            if int(row["subject"]) == int(cfg["subject"]["id"])]
    if (len(rows) != 1 or rows[0]["model_sha256"] != report["prefit_model_sha256"]
            or rows[0]["component_count"] != report["n_components_fitted"]
            or rows[0]["exclusions"] != sorted(exclusions)
            or set(exclusions).intersection(rows[0]["deferred_components_retained"])):
        raise RuntimeError("ICA choices differ from their recorded exact-model authorization.")


@lru_cache(maxsize=1)
def ica_stage():
    """Import existing model compatibility/QC helpers without running their CLI."""
    spec = importlib.util.spec_from_file_location("dad_epoch_ica_helpers", Path(__file__).with_name("03_C_ICA.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=1)
def epoch_stage():
    """Load the guarded epoch CLI for cohort preflight without executing it."""
    spec = importlib.util.spec_from_file_location("dad_epoch_stage", Path(__file__).with_name("04_epoching.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# -----------------------------------------------------------------------------
# Whole-participant validation precedes signal preparation and output writes.
# -----------------------------------------------------------------------------

def resolve_epoch_inputs(cfg, source_manifest):
    """Return exact native input records, honest preparation evidence and optional ICA."""
    mode = cfg["epoching"].get("preparation_mode", "saved_ica")
    base = upstream_config(cfg)
    model = None
    if mode == "saved_ica":
        runs, gate = resolve_manual_ica_apply_runs(base, source_manifest)
        evidence = {"mode": mode, "ica": {"exclusions": gate["manual_exclude"],
            "model_sha256": cfg.get("ica", {}).get("reviewed_model_sha256")}}
    else:
        if mode not in ("inline_no_ica", "inline_ica") or not cfg["epoching"].get("processing_run"):
            raise RuntimeError("Inline preparation needs an explicit mode and processing-run label.")
        dirs = stage_output_dirs(base)
        for stage in ("03_A_MaxFilt", "03_B_ArtefactAnnotation"):
            require_review(base, stage, dirs["reports"])
        records = discover_raw_run_records(base)
        runs = [ManualIcaApplyRun(source=record,
            cleaned_path=str(derivative_path_for_run(record, dirs["sss"], "_ann.fif")),
            manifest_output_path="", manifest_status="accepted_annotated_input") for record in records]
        settings = cfg["preprocessing"]
        if (any(settings.get(key) is not None for key in ("highpass_hz", "lowpass_hz", "resample_sfreq_hz")) or
            settings["notch_freqs_hz"] != [50., 100.] or settings["notch_method"] != "spectrum_fit" or
            settings["spectrum_fit_window"] != "10s"):
            raise RuntimeError("Inline preparation requires the established native spectrum-fit-only settings.")
        evidence = {"mode": mode, "ica": {"applied": False, "exclusions": None},
            "accepted_4A": base["stage4_review"]["03_A_MaxFilt"],
            "accepted_4B": base["stage4_review"]["03_B_ArtefactAnnotation"],
            "continuous_filter": settings,
            "spectrum_fit_minimum_samples": MIN_SPECTRUM_FIT_SAMPLES}
        if mode == "inline_ica":
            require_review(base, "03_C_ICA", dirs["reports"])
            report = json.loads((dirs["reports"] / "03_C_ICA_review_report.json").read_text())
            model_path = dirs["ica_models"] / Path(report["prefit_model"]).name
            metadata = model_path.with_name(model_path.stem + "_metadata.json")
            expected = ica_stage().build_ica_fingerprint(input_paths=[Path(run.cleaned_path) for run in runs], ica_cfg=base["ica"])
            if (not base["ica"].get("reuse_prefit") or not ica_stage().cache_is_compatible(metadata, expected) or
                file_hash(model_path) != report["prefit_model_sha256"] or
                base["ica"].get("reviewed_model_sha256") != report["prefit_model_sha256"] or
                report["configuration_sha256"] != stage_configuration(base, "03_C_ICA", include_edge=False)):
                raise RuntimeError("Inline ICA source/model/configuration binding changed; refitting is forbidden.")
            saved = json.loads(metadata.read_text())
            if saved["model_sha256"] != report["prefit_model_sha256"]:
                raise RuntimeError("Inline ICA metadata/model identity changed.")
            exclusions = base["ica"]["manual_exclude"]
            if (not isinstance(exclusions, list) or any(type(i) is not int or not 0 <= i < report["n_components_fitted"] for i in exclusions)
                    or len(set(exclusions)) != len(exclusions)):
                raise RuntimeError("Inline ICA requires explicit in-range model-bound exclusions.")
            validate_recorded_exclusions(base, report, exclusions)
            model = mne.preprocessing.read_ica(model_path, verbose="ERROR")
            if model.n_components_ != report["n_components_fitted"]:
                raise RuntimeError("Saved model component count differs from its review.")
            # MNE unions an explicit exclude argument with model.exclude.
            # Reset the disposable object's field so only the bound list acts.
            model.exclude = list(exclusions)
            evidence["ica"] = {"applied": True, "model_sha256": report["prefit_model_sha256"],
                "metadata_sha256": file_hash(metadata), "exclusions": sorted(set(exclusions))}
        gate = {"status": "passed", "processing_pass": "inline", "input_kind": "accepted_annotated_native",
            "n_raw_runs": len(runs), "ica_applied": mode == "inline_ica"}
    evidence["inputs"] = []
    trim_plans = []
    for run in runs:
        # Reading headers resolves every split part without loading signal data.
        # Annotation contents and native coordinates bind manual edits as well.
        print(f"Validating epoch input header: {run.source.source_filename}", flush=True)
        with mne.io.read_raw_fif(run.cleaned_path, preload=False, verbose="ERROR") as raw:
            evidence["inputs"].append(native_input_evidence(raw, run.source.source_filename))
            if cfg["epoching"].get("annotation_end_trims"):
                from epoch_annotation_trims import plan_trims
                trim_plans.extend(plan_trims(raw, cfg, run.source.source_filename))
    if cfg["epoching"].get("annotation_end_trims"):
        if len(trim_plans) != len(cfg["epoching"]["annotation_end_trims"]):
            raise RuntimeError("Not every approved annotation trim matched a native input")
        evidence["annotation_end_trims"] = {"plans": trim_plans,
            "implementation_sha256": file_hash(Path(__file__).with_name("epoch_annotation_trims.py"))}
    # Auxiliary reconstruction is an explicit new variant, never an implicit
    # repair of old outputs. Bind the whole participant plan before any work.
    if cfg["epoching"].get("auxiliary_reconstruction_manifest") is not None:
        from p24_auxiliary_reconstruction import validate_plan
        evidence["auxiliary_reconstruction"] = validate_plan(cfg, evidence["inputs"])
    gate["preparation"] = evidence
    return runs, gate, model


def epoch_contract(cfg, gate, timing_gate, *, events_csv=None):
    """Fingerprint data-changing inputs/operators; run labels identify storage only."""
    settings = {key: value for key, value in cfg["epoching"].items()
                if key not in ("processing_run", "timing_review_path")}
    payload = {"schema": "epoch-preparation-v1", "preparation": gate["preparation"],
        "epoching": settings, "timing_gate": portable(timing_gate),
        "trigger_events_sha256": file_hash(events_csv) if events_csv else None,
        "implementation_sha256": {name: file_hash(Path(__file__).with_name(name))
            for name in ("continuous_cleaning.py", "epoch_input_preparation.py", "04_epoching.py")}}
    return {"sha256": digest(payload), "definition": payload}


def artifact_binding(paths):
    """Bind all split FIF parts by stamp and small sidecars by content hash."""
    from review_final_ica import settle_review_outputs
    selected = {key: Path(path) for key, path in paths.items()
                if key in ("epochs_fif", "epoch_metadata_csv", "trial_decisions_csv")}
    parts = sorted(selected["epochs_fif"].parent.glob(selected["epochs_fif"].stem + "*.fif"))
    pattern = re.compile(re.escape(selected["epochs_fif"].stem) + r"(?:-\d+)?\.fif")
    parts = [path for path in parts if pattern.fullmatch(path.name)]
    settle_review_outputs([*selected.values(), *parts])
    result = {key: {"stamp": file_stamp(path), "sha256": file_hash(path)}
              for key, path in selected.items() if key != "epochs_fif"}
    result["epochs_fif"] = {"parts": [file_stamp(path) for path in parts]}
    return result


def verify_epoch_artifacts(cfg):
    """Require the completed epoch contract and unchanged stored outputs downstream."""
    reports = Path(cfg["resolved_paths"]["derivatives_reports_dir"])
    report = json.loads((reports / "04_epoching_report.json").read_text())
    if not report.get("preparation_contract"):
        if cfg["epoching"].get("processing_run"):
            raise RuntimeError("Processing-run epochs lack a preparation contract.")
        return report  # Historical unlabelled contracts remain readable.
    analysis = Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / "meg"
    for key, binding in report["artifact_bindings"].items():
        if key == "epochs_fif":
            if not binding["parts"] or any(file_stamp(analysis / part["name"]) != part for part in binding["parts"]):
                raise RuntimeError("Stored epoch FIF or split part changed.")
            continue
        path = analysis / Path(report["outputs"][key]).name
        if file_stamp(path) != binding["stamp"] or binding["sha256"] and file_hash(path) != binding["sha256"]:
            raise RuntimeError(f"Stored epoch output changed: {key}")
    definition = report["preparation_contract"]["definition"]
    if digest(definition) != report["preparation_contract"]["sha256"]:
        raise RuntimeError("Stored epoch preparation contract changed.")
    current = {key: value for key, value in cfg["epoching"].items() if key not in ("processing_run", "timing_review_path")}
    if definition["epoching"] != current:
        raise RuntimeError("Requested epoch settings differ from saved epochs; use a new run label.")
    if definition["preparation"]["mode"] != cfg["epoching"].get("preparation_mode", "saved_ica"):
        raise RuntimeError("Requested processing mode differs from saved epochs.")
    if definition["preparation"]["mode"] != "saved_ica":
        if definition["preparation"]["continuous_filter"] != cfg["preprocessing"]:
            raise RuntimeError("Continuous preparation settings changed; use a new run label.")
    if cfg["epoching"].get("auxiliary_reconstruction_manifest") is not None:
        from p24_auxiliary_reconstruction import validate_plan
        auxiliary = validate_plan(cfg, definition["preparation"]["inputs"])
        if auxiliary != definition["preparation"].get("auxiliary_reconstruction"):
            raise RuntimeError("Auxiliary reconstruction inputs changed; use a new run label.")
    return report


def export_contract(cfg, epoch_report, canonical_channels, trajectory_models):
    """Bind export operators and small trajectory files to their epoch variant."""
    # Trajectory discovery returns nested metadata; actual source MATs are
    # included by content below when the values identify existing files.
    def files(value):
        if isinstance(value, dict):
            return [item for child in value.values() for item in files(child)]
        if isinstance(value, (list, tuple)):
            return [item for child in value for item in files(child)]
        if isinstance(value, (str, Path)) and str(value).endswith(".mat") and Path(value).is_file():
            return [{"name": Path(value).name, "sha256": file_hash(value)}]
        return []
    return digest({"epoch_contract": epoch_report.get("preparation_contract"),
        "epoch_artifacts": epoch_report.get("artifact_bindings"), "export": cfg["export"],
        "canonical_channels": canonical_channels, "trajectory_files": files(trajectory_models)})


def verify_export_artifacts(cfg):
    """Require unchanged epoch and export outputs before MATLAB or resume."""
    epoch = verify_epoch_artifacts(cfg)
    report_path = Path(cfg["resolved_paths"]["derivatives_reports_dir"]) / "06_export_neural_for_drsa_report.json"
    report = json.loads(report_path.read_text())
    if cfg["epoching"].get("processing_run"):
        if report.get("status") != "completed" or not report.get("export_contract_sha256"):
            raise RuntimeError("Processing-run export lacks a completed bound contract.")
        if report.get("epoch_contract_sha256") != epoch["preparation_contract"]["sha256"]:
            raise RuntimeError("Export belongs to different epochs.")
        from downstream_stage_utils import discover_trajectory_models
        canonical = Path(cfg["resolved_paths"]["project_root"]) / cfg["export"]["canonical_meg_channels_file"]
        channels = [line.strip() for line in canonical.read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith("#")]
        if report["export_contract_sha256"] != export_contract(cfg, epoch, channels, discover_trajectory_models(cfg)):
            raise RuntimeError("Export inputs/settings changed; use a new processing-run label.")
        root = Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / "drsa"
        for binding in report["artifact_bindings"]:
            path = root / binding["stamp"]["name"]
            if file_stamp(path) != binding["stamp"] or binding.get("sha256") and file_hash(path) != binding["sha256"]:
                raise RuntimeError(f"Stored neural export changed: {path.name}")
    return report


# -----------------------------------------------------------------------------
# Native signal operations happen once per recording, in the intended order.
# -----------------------------------------------------------------------------

def filter_continuous_for_epoching(raw, settings, *, source_filename=None):
    """Run native line removal while leaving unsupported tiny fragments intact.

    Input is a preloaded native Raw with accepted annotations and spectrum-fit-
    only settings. Use MNE's exact annotation/sample segmentation, including its
    rounding and zero-duration boundaries. On recordings with short fragments,
    call the same public array filter that Raw.notch_filter uses, in place on
    supported MEG spans only. The other recordings use the original helper.

    Return the existing filter summary plus absolute sample ranges and total
    duration of skipped fragments. No annotation, sample position, auxiliary
    channel or epoch selection is changed. ICA, when selected, has already been
    applied; this bypass concerns only line-noise removal. The protected
    continuous_cleaning module and saved-model fingerprints remain unchanged.
    """
    if (any(settings.get(key) is not None for key in
            ("highpass_hz", "lowpass_hz", "resample_sfreq_hz"))
            or settings.get("notch_method", "spectrum_fit") != "spectrum_fit"):
        raise ValueError("Inline preparation requires native spectrum-fit-only settings.")
    if not raw.preload:
        raise RuntimeError("Inline continuous preparation requires preloaded native data.")
    # Match Raw.notch_filter's exact mask, including its sample rounding.
    starts, stops = _annotations_starts_stops(raw, ("bad", "edge"), invert=True)
    short = (stops - starts) < MIN_SPECTRUM_FIT_SAMPLES
    skipped = [{"first_sample_raw": int(a) + int(raw.first_samp),
        "stop_sample_raw_exclusive": int(b) + int(raw.first_samp),
        "n_samples": int(b - a)} for a, b in zip(starts[short], stops[short])]
    if skipped:
        summary = final_filter_summary(raw, settings, source_filename=source_filename)
        picks = mne.pick_types(raw.info, meg=True, exclude=[])
        for start, stop in zip(starts[~short], stops[~short]):
            # Match the public array call made by Raw.notch_filter. Slicing
            # retains native time and cannot bridge a BAD/EDGE boundary.
            mne.filter.notch_filter(raw._data[:, start:stop], raw.info["sfreq"],
                settings.get("notch_freqs_hz", [50., 100.]), picks=picks,
                method="spectrum_fit", filter_length=settings.get("spectrum_fit_window", "10s"),
                copy=False, verbose="WARNING")
        print(f"{source_filename}: line removal skipped {len(skipped)} tiny fragments "
              f"({sum(row['n_samples'] for row in skipped)} native samples)", flush=True)
    else:
        summary = filter_continuous(raw, settings, source_filename=source_filename)
    summary["spectrum_fit"].update({
        "minimum_segment_samples": MIN_SPECTRUM_FIT_SAMPLES,
        "short_fragment_policy": "leave_samples_and_annotations_unchanged",
        "skipped_fragments": skipped,
        "skipped_seconds": sum(row["n_samples"] for row in skipped) / raw.info["sfreq"],
        "processed_segment_count": int(np.count_nonzero(~short))})
    return summary


def prepare_run(cfg, run, evidence, model, figure_dir):
    """Return prepared Raw/QC, including an explicitly bound P24 auxiliary plan.

    Native preservation checks precede the optional auxiliary-only replacement.
    The separate adapter verifies all retained samples, acquisition settings and
    annotations after replacement. Original FIFs and old variants stay intact.
    """
    raw = mne.io.read_raw_fif(run.cleaned_path, preload=evidence["mode"] != "saved_ica", verbose="ERROR")
    expected = next(row for row in evidence["inputs"] if row["source_filename"] == run.source.source_filename)
    if native_input_evidence(raw, run.source.source_filename) != expected:
        raise RuntimeError("Native input changed after validation; no epochs may be produced.")
    if evidence["mode"] == "saved_ica":
        return raw, None
    # Apply the two source/trial-bound approvals before line removal. Preserve
    # the original native-input binding above; subsequent checks use this exact
    # operative mask and no original FIF or accepted mask file is rewritten.
    trims = []
    if "annotation_end_trims" in evidence:
        from epoch_annotation_trims import apply_trims
        if file_hash(Path(__file__).with_name("epoch_annotation_trims.py")) != evidence["annotation_end_trims"]["implementation_sha256"]:
            raise RuntimeError("Annotation trim implementation changed after validation")
        trims = apply_trims(raw, cfg, run.source.source_filename)
    before = raw.copy()
    print(f"{run.source.source_filename}: {evidence['mode']} then spectrum-fit 50/100 Hz", flush=True)
    if model is not None:
        model.apply(raw, exclude=evidence["ica"]["exclusions"], verbose="ERROR")
    figure_dir = guarded_path(figure_dir, cfg["resolved_paths"]["data_root"])
    figure_dir.mkdir(parents=True, exist_ok=True)
    # Exact existing helpers supply comparable native windows and the same
    # bounded spectrum snapshots; no fit copy or new temporal filter is used.
    comparison = None
    if model is not None:
        comparison = str(ica_stage().save_native_comparison(before, raw, run.source, figure_dir))
    spectral_before = ica_stage().spectral_snapshot(raw)
    filtering = filter_continuous_for_epoching(raw, cfg["preprocessing"], source_filename=run.source.source_filename)
    spectral_after = ica_stage().spectral_snapshot(raw, spectral_before["window"])
    figure, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for axis, kind in zip(axes, ("mag", "grad")):
        for snapshot, label in ((spectral_before, "Before line removal"), (spectral_after, "After 50/100 Hz removal")):
            if kind in snapshot["types"]:
                row = snapshot["types"][kind]
                axis.plot(row["frequencies"], 10*np.log10(np.maximum(row["mean_psd"], np.finfo(float).tiny)), label=label)
        axis.set_ylabel(kind + " PSD (dB)")
        if kind in spectral_before["types"]: axis.legend()
    axes[-1].set_xlabel("Frequency (Hz)"); figure.tight_layout()
    path = figure_dir / f"{run.source.source_stem}_line_removal.png"
    figure.savefig(path, dpi=120); plt.close(figure)
    preservation = check_preservation(before, raw, annotations=True)
    if [file_stamp(path) for path in raw.filenames] != expected["parts"]:
        raise RuntimeError("Native input changed during continuous preparation.")
    for key in ("dev_head_t", "dev_ctf_t", "ctf_head_t"):
        left, right = before.info.get(key), raw.info.get(key)
        if (left is None) != (right is None) or left is not None and not np.array_equal(left["trans"], right["trans"]):
            raise RuntimeError(f"Native coordinate transform changed: {key}")
    result = {"preservation": preservation, "final_filter": filtering,
        "line_removal_figure": str(path), "ica_only_figure": comparison,
        "continuous_file_saved": False, "input_annotations_sha256": digest(annotation_rows(before))}
    if trims:
        result["input_annotations_sha256"] = expected["annotations_sha256"]
        result["operative_annotations_sha256"] = digest(annotation_rows(before))
        result["annotation_end_trims"] = trims
    # P10 remains an open concern: show native post-preparation sub-Hz spectra
    # and long-window traces without adding high-pass or display detrending.
    if run.source.source_filename.startswith("sub10_"):
        picks = mne.pick_types(raw.info, meg="grad", exclude="bads")[:4]
        width = min(raw.n_times, round(60*raw.info["sfreq"]))
        starts = sorted(set(int(v) for v in np.linspace(0, raw.n_times-width, 3)))
        fig, axes = plt.subplots(len(starts), 2, figsize=(13, 3*len(starts)), squeeze=False)
        for index, start in enumerate(starts):
            values = raw.get_data(picks=picks, start=start, stop=start+width)
            times = (np.arange(start, start+width)+raw.first_samp)/raw.info["sfreq"]
            axes[index, 0].plot(times, values.T*1e13)
            from scipy.signal import welch
            frequencies, power = welch(values, raw.info["sfreq"], nperseg=width)
            selected = (frequencies >= .1) & (frequencies <= 2.)
            axes[index, 1].plot(frequencies[selected], np.mean(power, axis=0)[selected])
            axes[index, 1].axvline(1/3, linestyle=":", color="red")
        fig.suptitle(f"{run.source.source_filename}: residual slow signal after preparation; retention OPEN")
        fig.tight_layout(); slow_path = figure_dir / f"{run.source.source_stem}_slow_qc.png"
        fig.savefig(slow_path, dpi=120); plt.close(fig); result["p10_slow_qc"] = str(slow_path)
    if "auxiliary_reconstruction" in evidence:
        from p24_auxiliary_reconstruction import apply_plan
        result["auxiliary_reconstruction"] = apply_plan(raw, cfg, run.source.source_filename,
            evidence["auxiliary_reconstruction"])
    return raw, result
