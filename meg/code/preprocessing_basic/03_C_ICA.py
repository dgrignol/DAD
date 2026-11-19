#!/usr/bin/env python3
"""Stage 4C: fit, review and apply source-bound manual ICA.

Usage after independent 4A/4B acceptance:
  python meg/code/preprocessing_basic/03_C_ICA.py \
    --config meg/params/pipeline_default.yaml --subject-config meg/params/sub-04.yaml \
    --tuning-config /path/to/03_tuning_choices.yaml --processing-pass review
Then record the reconciled component decision with review_cleaning.py --stage
4C --exclude [indices]; explicit empty is valid. Apply with the same arguments
and --processing-pass apply --overwrite, or --resume after an interruption.

Data flow: annotated native acquisitions -> separate admissible MEG stretches ->
1-40 Hz FIR + full edge guards -> isolated 200 Hz resampling -> compatible
concatenated fit copies -> MNE ICA(method='fastica', n_components=None, seed 97)
-> all-component diagnostics -> exact-model exclusions -> each original native
run -> ICA-only before/after plots -> spectrum-fit 50/100 Hz removal, without
additional high-pass/low-pass -> final continuous FIF and preservation evidence.

Every BAD category, including routine ocular intervals, is omitted from final
fitting. Preliminary frontal ICA is only a separate annotation aid and never a
production cache. Source coordinate maps and fit sample counts are recorded.
Common transforms/channel/bad-channel bases must agree; mismatches are errors.
Reports embed PNGs and explain their interpretation and review limits.
The component book is the primary human review interface; the standard report
records processing/evidence. Every figure has a model-bound reproduction command.
Properties are two-second non-task segments; separate diagnostic task plots use
accepted Stage 3 motion-onsets. Full-length display copies work around the MNE
1.10 concatenated-source metadata issue without changing or refitting ICA.
An explicitly reconciled deterministic restart retains its initializer ancestry
and separate phase counts in metadata/reports. Such a cache requires explicit
reconciliation if incompatible; it cannot silently revert to a default fit.
No EOG workflow executes. Spectrum-fit boundaries/short stretches need review;
no FIR margin is required without a high-pass. Historical exact-source warnings
protect only a future high-pass, including P38 run 1; line removal always runs.
A changed fit input requires a new model/component review; apply never silently
refits. All native coordinates, auxiliary channels and original annotations are
preserved. No epochs, eligibility edits, exports or analyses are performed.
An explicit provisional component authorization is documented in ICA application
reports. It records processing permission and the need for later manual inspection.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.preprocessing import ICA, read_ica

from continuous_cleaning import (stage_configuration, annotation_rows, bind_report, check_preservation, coordinates,
    digest, file_hash, file_stamp, filter_continuous, final_filter_summary, prepare_fit_copy, require_review)

from pipeline_config import load_and_resolve_config
from preprocessing_stage_utils import (
    AcquisitionRunRecord,
    derivative_path_for_run,
    discover_raw_run_records,
    ensure_output_directory,
    guarded_path,
    lineage_record,
    prepare_output_path,
    stage_output_dirs,
    update_source_manifest_stage,
    update_tuning_file,
    utc_now,
    write_json,
)

from render_stage03_mne_reports import render_stage_report
from final_ica_review_book import write_component_book
from ica_review_plots import continuous_display_copy, describe_properties, native_context
from ica_optimizer_checkpoint import iteration_cap


STAGE_NAME = "ica_manual"
INPUT_SUFFIX = "_ann.fif"
OUTPUT_SUFFIX = "_ica.fif"


# -----------------------------------------------------------------------------
# Cache identity and fit-data preparation
# -----------------------------------------------------------------------------


def build_ica_fingerprint(
    *,
    input_paths: list[Path],
    ica_cfg: Mapping[str, Any],
) -> dict[str, Any]:
    """
    Build a deterministic cache fingerprint from inputs and fit parameters.

    Stable derivative filenames, file size, and nanosecond modification time
    detect changed inputs without hashing multi-gigabyte FIF contents. Absolute
    paths are deliberately excluded so the same shared file remains compatible
    through the macOS ``/Volumes`` and cluster ``/mnt/storage`` mount points.
    """
    payload = {
        "schema_version": 3,
        "inputs": [
            {
                "source_filename": path.name,
                "size_bytes": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
                "annotations": digest(annotation_rows(mne.io.read_raw_fif(path, preload=False, verbose="ERROR"))),
            }
            for path in input_paths
        ],
        "fit": {
            "method": ica_cfg["method"],
            "n_components": ica_cfg["n_components"],
            "random_state": ica_cfg["random_state"],
            "max_iter": ica_cfg["max_iter"],
            "fit_resample_sfreq_hz": ica_cfg["fit_resample_sfreq_hz"],
            "fit_l_freq_hz": ica_cfg["fit_l_freq_hz"],
            "fit_h_freq_hz": ica_cfg["fit_h_freq_hz"],
            "fit_resample_edge_sec": ica_cfg.get("fit_resample_edge_sec", 0.1),
            "fit_min_good_sec": ica_cfg.get("fit_min_good_sec", 1.0),
            "mne_version": mne.__version__,
            "protected_preparation_sha256": file_hash(Path(__file__).with_name("continuous_cleaning.py")),
            "fit_policy": "omit_all_BAD_and_filter_edges_v1",
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["sha256"] = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def _portable_fingerprint_sha256(
    fingerprint: Mapping[str, Any],
) -> str | None:
    """
    Return a mount-independent hash for current or legacy fingerprints.

    Legacy metadata stored absolute paths. Reducing those paths to basenames
    permits reuse of a model fitted from the same shared derivatives before
    site-config support was introduced.
    """
    inputs = fingerprint.get("inputs")
    fit = fingerprint.get("fit")
    if not isinstance(inputs, list) or not isinstance(fit, Mapping):
        return None
    normalized_inputs = []
    for item in inputs:
        if not isinstance(item, Mapping):
            return None
        filename = item.get("source_filename")
        if filename is None and item.get("path") is not None:
            filename = Path(str(item["path"])).name
        if filename is None:
            return None
        normalized_inputs.append(
            {
                "source_filename": str(filename),
                "size_bytes": item.get("size_bytes"),
                "mtime_ns": item.get("mtime_ns"),
                "annotations": item.get("annotations"),
            }
        )
    portable = {
        "schema_version": 2,
        "inputs": normalized_inputs,
        "fit": dict(fit),
    }
    canonical = json.dumps(portable, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def cache_is_compatible(
    metadata_path: Path,
    expected_fingerprint: Mapping[str, Any],
) -> bool:
    """Return whether an ICA cache matches across local/cluster mount paths."""
    if not metadata_path.is_file():
        return False
    try:
        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return _portable_fingerprint_sha256(
        payload.get("fingerprint", {})
    ) == _portable_fingerprint_sha256(
        expected_fingerprint
    )


def prepare_ica_fit_data(input_paths, ica_cfg):
    """Build protected per-source fit stretches, then concatenate compatible data.

    No filtering or resampling crosses any BAD interval. Final spectrum-fit-only
    processing is independent of these disposable 1-40 Hz fit copies.
    """
    prepared, summaries = [], []
    for input_path in input_paths:
        print(f"Final ICA fit preparation: {input_path.name}", flush=True)
        raw = mne.io.read_raw_fif(input_path, preload=False, verbose="WARNING")
        fit_raw, summary = prepare_fit_copy(raw, ica_cfg)
        summary.update(input_path=str(input_path), n_channels=len(fit_raw.ch_names),
                       n_times=fit_raw.n_times, sfreq_hz=fit_raw.info["sfreq"],
                       bad_channels=fit_raw.info["bads"])
        prepared.append(fit_raw)
        summaries.append(summary)
    result = mne.concatenate_raws(prepared, on_mismatch="raise", verbose="WARNING")
    if result.n_times != sum(row["fit_samples"] for row in summaries):
        raise RuntimeError("Final ICA sample count disagrees with native interval maps.")
    return result, summaries


def fit_or_load_ica(
    *,
    input_paths: list[Path],
    ica_cfg: Mapping[str, Any],
    prefit_path: Path,
    metadata_path: Path,
    raw_root: str | Path,
    overwrite: bool,
    require_fit_data: bool,
) -> tuple[ICA, mne.io.BaseRaw | None, list[dict[str, Any]], str, dict[str, Any]]:
    """
    Load a compatible cache or fit a new ICA model.

    ``require_fit_data`` is true for review diagnostics. Apply runs can reuse a
    compatible model without concatenating all fit copies again. An explicitly
    reconciled restart cannot silently fall back to the original initialization.
    """
    fingerprint = build_ica_fingerprint(
        input_paths=input_paths,
        ica_cfg=ica_cfg,
    )
    if str(ica_cfg["method"]).lower() == "fastica" and (
        importlib.util.find_spec("sklearn") is None
    ):
        raise RuntimeError(
            "scikit-learn is required for the production FastICA stage. "
            "Install meg/requirements.txt before running ICA."
        )
    can_reuse = (
        bool(ica_cfg["reuse_prefit"])
        and prefit_path.is_file()
        and cache_is_compatible(metadata_path, fingerprint)
    )
    fit_data: mne.io.BaseRaw | None = None
    fit_summaries: list[dict[str, Any]] = []
    if not can_reuse and metadata_path.is_file():
        saved_metadata = json.loads(metadata_path.read_text())
        if saved_metadata.get("optimizer_restart"):
            raise RuntimeError("Reconciled ICA restart cannot be silently refitted. Reconcile changed input/configuration/model lineage explicitly.")
    if can_reuse:
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("model_sha256") != file_hash(prefit_path):
            raise RuntimeError("ICA cache model contents changed since fitting.")
        fit_summaries = metadata.get("fit_input_summaries", [])
        model = read_ica(prefit_path, verbose="WARNING")
        cache_status = "reused_compatible"
        if require_fit_data:
            fit_data, fit_summaries = prepare_ica_fit_data(input_paths, ica_cfg)
        return model, fit_data, fit_summaries, cache_status, fingerprint

    fit_data, fit_summaries = prepare_ica_fit_data(input_paths, ica_cfg)
    print(f"Final ICA: fitting {fit_data.n_times} samples across {len(input_paths)} acquisitions", flush=True)
    model = ICA(
        method=str(ica_cfg["method"]),
        n_components=ica_cfg["n_components"],
        random_state=int(ica_cfg["random_state"]),
        max_iter=ica_cfg["max_iter"],
        verbose="WARNING",
    )
    model.fit(fit_data, reject_by_annotation=True, verbose="WARNING")
    if int(model.n_samples_) != fit_data.n_times:
        raise RuntimeError("ICA omitted samples outside the reconciled fit interval map.")
    prepared_prefit = prepare_output_path(
        prefit_path,
        raw_root=raw_root,
        overwrite=overwrite,
    )
    model.save(prepared_prefit, overwrite=overwrite, verbose="WARNING")
    write_json(
        metadata_path,
        {
            "generated_utc": utc_now(),
            "fingerprint": fingerprint,
            "model_path": str(prefit_path),
            "n_components_fitted": int(model.n_components_),
            "model_sha256": file_hash(prepared_prefit),
            "fit_input_summaries": fit_summaries,
            "n_samples_fitted": int(model.n_samples_),
            "n_iterations": int(model.n_iter_),
            "pca_explained_variance": model.pca_explained_variance_.tolist(),
        },
        raw_root=raw_root,
        overwrite=True,
    )
    return model, fit_data, fit_summaries, "fitted_new", fingerprint


# -----------------------------------------------------------------------------
# Diagnostics and interactive review
# -----------------------------------------------------------------------------


def summarize_eog_scores(
    *,
    model: ICA,
    input_paths: list[Path],
    configured_channels: list[str],
) -> dict[str, Any]:
    """
    Compute diagnostic component-EOG rankings only when channels are present.

    Failures are recorded per file/channel and never trigger automatic
    exclusion or a fallback cleaning method.
    """
    if not configured_channels:
        return {
            "status": "skipped_no_configured_eog_channels",
            "per_file": [],
            "aggregate_ranking": [],
        }
    component_scores: dict[int, list[float]] = defaultdict(list)
    per_file: list[dict[str, Any]] = []
    for input_path in input_paths:
        raw = mne.io.read_raw_fif(input_path, preload=True, verbose="ERROR")
        available = [
            channel for channel in configured_channels if channel in raw.ch_names
        ]
        file_record: dict[str, Any] = {
            "input_path": str(input_path),
            "available_eog_channels": available,
            "channels": [],
        }
        if not available:
            file_record["status"] = "skipped_no_eog_channels"
            per_file.append(file_record)
            continue
        file_record["status"] = "completed"
        for channel in available:
            try:
                suggested, scores = model.find_bads_eog(
                    raw,
                    ch_name=channel,
                    verbose=False,
                )
                top_indices = np.argsort(np.abs(scores))[::-1][:10]
                top = []
                for index in top_indices:
                    value = float(scores[index])
                    component_scores[int(index)].append(abs(value))
                    top.append(
                        {
                            "component": int(index),
                            "score": value,
                            "absolute_score": abs(value),
                        }
                    )
                file_record["channels"].append(
                    {
                        "channel": channel,
                        "status": "completed",
                        "suggested_components": [
                            int(value) for value in suggested
                        ],
                        "top_scores": top,
                    }
                )
            except Exception as exc:  # pylint: disable=broad-except
                file_record["channels"].append(
                    {
                        "channel": channel,
                        "status": "failed",
                        "error": str(exc),
                    }
                )
        per_file.append(file_record)

    ranking = sorted(
        (
            {
                "component": component,
                "mean_absolute_score": float(np.mean(values)),
                "n_scores": len(values),
            }
            for component, values in component_scores.items()
        ),
        key=lambda item: item["mean_absolute_score"],
        reverse=True,
    )
    completed = any(
        channel.get("status") == "completed"
        for item in per_file
        for channel in item["channels"]
    )
    return {
        "status": "completed" if completed else "skipped_or_failed",
        "per_file": per_file,
        "aggregate_ranking": ranking,
    }


def _save_figure(
    fig: Any,
    path: Path,
    *,
    show: bool,
) -> Path:
    """Save one MNE/matplotlib figure and optionally display it."""
    fig.savefig(path, dpi=150, bbox_inches="tight")
    if show:
        plt.show(block=True)
    plt.close(fig)
    return path


def generate_component_diagnostics(
    *,
    model: ICA,
    fit_data: mne.io.BaseRaw,
    figure_dir: Path,
    subject_label: str,
    processing_pass: str,
    raw_root: str | Path,
    show: bool,
    include_properties: bool,
) -> list[Path]:
    """
    Save source, topography, and optional per-component property diagnostics.

    Topography pages may be displayed interactively. Per-component property
    figures are saved for detailed inspection without forcing 30 sequential GUI
    windows. MNE's source browser uses a backend separate from matplotlib's
    normal backend; noninteractive runs force its matplotlib implementation so
    headless cluster jobs never attempt to connect to Qt/X11.
    """
    images: list[Path] = []
    # A single display RawArray avoids MNE 1.10 get_sources() truncating the
    # metadata of concatenated good stretches to their first piece. No fit changes.
    display_data = continuous_display_copy(fit_data)
    if show:
        source_fig = model.plot_sources(
            display_data,
            title=f"{subject_label} ICA sources",
            show=True,
        )
    else:
        with mne.viz.use_browser_backend("matplotlib"):
            source_fig = model.plot_sources(
                display_data,
                title=f"{subject_label} ICA sources",
                show=False,
            )
    source_path = prepare_output_path(
        figure_dir / f"03_C_ICA_{processing_pass}_sources.png",
        raw_root=raw_root,
        overwrite=True,
    )
    images.append(_save_figure(source_fig, source_path, show=False))

    component_figures = model.plot_components(show=False)
    if not isinstance(component_figures, list):
        component_figures = [component_figures]
    for page, figure in enumerate(component_figures, start=1):
        path = prepare_output_path(
            figure_dir
            / f"03_C_ICA_{processing_pass}_topographies_{page:02d}.png",
            raw_root=raw_root,
            overwrite=True,
        )
        images.append(_save_figure(figure, path, show=show))

    if include_properties:
        # Batch components so MNE projects/segments once per batch rather than
        # once per IC, while bounding simultaneously open figures and arrays.
        for first in range(0, int(model.n_components_), 8):
            picks = list(range(first, min(first + 8, int(model.n_components_))))
            print(f"Final ICA diagnostics: components {picks[0]}-{picks[-1]} (full fit display)", flush=True)
            property_figures = model.plot_properties(display_data, picks=picks,
                psd_args={"fmax": float(display_data.info["lowpass"])}, show=False, verbose=False)
            for component, figure in zip(picks, property_figures, strict=True):
                metadata = describe_properties(figure, display_data.info["sfreq"])
                metadata.update(scope="full_fit", display_samples=int(display_data.n_times),
                    mne_display_version=mne.__version__)
                path = prepare_output_path(figure_dir / f"03_C_ICA_{processing_pass}_ICA{component:03d}_properties_01.png",
                    raw_root=raw_root, overwrite=True)
                write_json(path.with_suffix(".json"), metadata, raw_root=raw_root)
                images.append(_save_figure(figure, path, show=False))
    return images


def parse_component_indices(text: str) -> list[int]:
    """Parse comma/space separated non-negative component indices."""
    if not text.strip():
        return []
    values = []
    for token in text.replace(",", " ").split():
        value = int(token)
        if value < 0:
            raise ValueError("ICA component indices must be non-negative.")
        values.append(value)
    return sorted(set(values))


def interactive_component_review(
    *,
    n_components: int,
    initial_exclusions: list[int] | None,
    diagnostic_paths: list[Path],
    eog_diagnostics: Mapping[str, Any],
) -> list[int]:
    """Prompt repeatedly until the operator accepts an exclusion list."""
    print("\nICA diagnostic files:")
    for path in diagnostic_paths:
        print(f"  {path}")
    ranking = eog_diagnostics.get("aggregate_ranking", [])
    if ranking:
        print("\nEOG diagnostic ranking (not auto-applied):")
        for item in ranking[:10]:
            print(
                f"  ICA{item['component']:03d}: "
                f"mean|score|={item['mean_absolute_score']:.3f} "
                f"(n={item['n_scores']})"
            )
    exclusions = list(initial_exclusions or [])
    while True:
        default_text = ",".join(str(value) for value in exclusions)
        answer = input(
            "ICA components to reject as comma/space separated indices "
            f"(blank means none; current [{default_text}]): "
        ).strip()
        if answer:
            try:
                exclusions = parse_component_indices(answer)
            except (TypeError, ValueError) as exc:
                print(f"Invalid component list: {exc}")
                continue
        else:
            exclusions = []
        invalid = [value for value in exclusions if value >= n_components]
        if invalid:
            print(
                f"Indices outside available range 0-{n_components - 1}: {invalid}"
            )
            continue
        confirmation = input(
            f"Accept ICA exclusions {exclusions}? [Y/n]: "
        ).strip().lower()
        if confirmation in ("", "y", "yes"):
            return exclusions


# -----------------------------------------------------------------------------
# ICA application
# -----------------------------------------------------------------------------


def _validate_resumable_ica_output(
    *,
    record: AcquisitionRunRecord,
    input_path: Path,
    output_path: Path,
    exclusions: list[int],
    raw_root: str | Path,
) -> dict[str, Any]:
    """Validate one existing cleaned acquisition run before resume skips it."""
    resolved = guarded_path(output_path, raw_root)
    try:
        source = mne.io.read_raw_fif(
            input_path,
            preload=False,
            verbose="ERROR",
        )
        existing = mne.io.read_raw_fif(
            resolved,
            preload=False,
            verbose="ERROR",
        )
    except Exception as exc:  # pylint: disable=broad-except
        raise RuntimeError(
            f"Cannot resume from unreadable ICA derivative {resolved}: {exc}"
        ) from exc
    description = str(existing.info.get("description") or "")
    required_tokens = (
        f"source={record.source_filename}",
        f"acquisition_run={record.acquisition_run_label}",
        "stage=03_C_ICA",
        f"excluded={exclusions}",
    )
    missing_tokens = [token for token in required_tokens if token not in description]
    if missing_tokens:
        raise RuntimeError(
            "Cannot resume from ICA derivative with mismatched lineage or "
            f"exclusions: {resolved}; missing={missing_tokens}"
        )
    if existing.ch_names != source.ch_names:
        raise RuntimeError(
            f"Cannot resume because ICA changed channel identity: {resolved}"
        )
    if list(existing.info.get("bads", [])) != list(source.info.get("bads", [])):
        raise RuntimeError(
            f"Cannot resume because ICA changed bad-channel metadata: {resolved}"
        )
    if existing.n_times != source.n_times:
        raise RuntimeError(
            "Cannot resume because ICA changed the acquisition sample count: "
            f"{resolved}; input={source.n_times}, output={existing.n_times}"
        )
    if float(existing.info["sfreq"]) != float(source.info["sfreq"]):
        raise RuntimeError(
            f"Cannot resume because ICA changed sampling frequency: {resolved}"
        )
    check_preservation(source, existing, annotations=False)
    misc_channels = [
        channel
        for channel, kind in zip(
            existing.ch_names,
            existing.get_channel_types(),
            strict=True,
        )
        if kind == "misc"
    ]
    source_misc = [
        channel
        for channel, kind in zip(
            source.ch_names,
            source.get_channel_types(),
            strict=True,
        )
        if kind == "misc"
    ]
    if misc_channels != source_misc:
        raise RuntimeError(
            "Cannot resume because ICA changed MISC channel identity: "
            f"{resolved}; input={source_misc}, output={misc_channels}"
        )
    return {
        "source_filename": record.source_filename,
        "acquisition_run_label": record.acquisition_run_label,
        "execution_status": "resumed_existing",
        "input_path": str(input_path),
        "output_path": str(resolved),
        "excluded_components": exclusions,
        "bad_channels_preserved": list(existing.info.get("bads", [])),
        "misc_channels_preserved": misc_channels,
    }


def apply_ica_to_runs(
    *,
    model: ICA,
    exclusions: list[int],
    records: list[AcquisitionRunRecord],
    input_paths: list[Path],
    sss_dir: Path,
    raw_root: str | Path,
    overwrite: bool,
    resume: bool = False,
    final_settings=None,
    model_sha256=None,
    figure_dir=None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Apply reviewed exclusions to each original annotated acquisition run."""
    model.exclude = list(exclusions)
    run_reports: list[dict[str, Any]] = []
    lineage_records: list[dict[str, Any]] = []
    for record, input_path in zip(records, input_paths, strict=True):
        requested_output = derivative_path_for_run(
            record,
            sss_dir,
            OUTPUT_SUFFIX,
        )
        if resume and requested_output.exists():
            run_report = _validate_resumable_ica_output(
                record=record,
                input_path=input_path,
                output_path=requested_output,
                exclusions=exclusions,
                raw_root=raw_root,
            )
            header = mne.io.read_raw_fif(requested_output, preload=False, verbose="ERROR")
            if model_sha256:
                application_sha = digest({"model": model_sha256, "input": file_stamp(input_path),
                                          "exclude": exclusions, "final_filter": final_settings,
                                          "final_filter_contract": "independent_spectrum_fit_2026-10-28"})
                if f"application_sha256={application_sha}" not in str(header.info.get("description")):
                    raise RuntimeError("Cleaned output belongs to different input/model/filter settings.")
            if final_settings is not None:
                run_report["final_filter"] = final_filter_summary(header, final_settings, source_filename=record.source_filename)
            run_reports.append(run_report)
            lineage_records.append(
                lineage_record(
                    source=record,
                    input_path=input_path,
                    output_path=run_report["output_path"],
                    status="resumed_existing",
                    metadata={
                        "excluded_components": exclusions,
                        "bad_channels_preserved": run_report[
                            "bad_channels_preserved"
                        ],
                        "misc_channels_preserved": run_report[
                            "misc_channels_preserved"
                        ],
                        "execution_status": "resumed_existing",
                    },
                )
            )
            continue
        output_path = prepare_output_path(
            requested_output,
            raw_root=raw_root,
            overwrite=overwrite,
        )
        raw = mne.io.read_raw_fif(
            input_path,
            preload=True,
            verbose="WARNING",
        )
        source = raw.copy()
        input_bads = list(raw.info.get("bads", []))
        misc_before = [
            channel
            for channel, kind in zip(
                raw.ch_names,
                raw.get_channel_types(),
                strict=True,
            )
            if kind == "misc"
        ]
        model.apply(raw, exclude=exclusions, verbose="WARNING")
        raw.info["bads"] = input_bads
        comparison = None
        slow_qc = None
        if figure_dir is not None:
            comparison = save_native_comparison(source, raw, record, Path(figure_dir))
            if record.source_filename.startswith("sub10_"):
                slow_qc = save_p10_slow_qc(source, raw, record, Path(figure_dir))
        spectral_before = spectral_snapshot(raw) if final_settings is not None else None
        final_summary = filter_continuous(raw, final_settings, source_filename=record.source_filename) if final_settings is not None else None
        final_figure = save_final_spectra(spectral_before, raw, record, Path(figure_dir)) if spectral_before is not None and figure_dir is not None else None
        check_preservation(source, raw, annotations=final_settings is None)
        raw.info["description"] = (
            f"DAD source={record.source_filename}; "
            f"acquisition_run={record.acquisition_run_label}; "
            f"stage=03_C_ICA; excluded={exclusions}"
        )
        misc_after = [
            channel
            for channel, kind in zip(
                raw.ch_names,
                raw.get_channel_types(),
                strict=True,
            )
            if kind == "misc"
        ]
        if misc_after != misc_before:
            raise RuntimeError(
                "MISC channel identity changed during ICA application: "
                f"before={misc_before}, after={misc_after}"
            )
        if model_sha256:
            application_sha = digest({"model": model_sha256, "input": file_stamp(input_path),
                                      "exclude": exclusions, "final_filter": final_settings,
                                      "final_filter_contract": "independent_spectrum_fit_2026-10-28"})
            raw.info["description"] += f"; application_sha256={application_sha}"
        raw.save(output_path, overwrite=overwrite, verbose="WARNING")
        saved = mne.io.read_raw_fif(output_path, preload=False, verbose="ERROR")
        preservation = check_preservation(source, saved, annotations=final_settings is None)
        run_report = {
            "source_filename": record.source_filename,
            "acquisition_run_label": record.acquisition_run_label,
            "execution_status": "completed",
            "input_path": str(input_path),
            "output_path": str(output_path),
            "excluded_components": exclusions,
            "preservation": preservation, "final_filter": final_summary,
            "before_after_ica_figure": str(comparison) if comparison else None,
            "p10_slow_contamination_qc": slow_qc,
            "final_filter_psd_figure": str(final_figure) if final_figure else None,
            "bad_channels_preserved": input_bads,
            "misc_channels_preserved": misc_after,
        }
        run_reports.append(run_report)
        lineage_records.append(
            lineage_record(
                source=record,
                input_path=input_path,
                output_path=output_path,
                status="completed",
                metadata={
                    "excluded_components": exclusions,
                    "bad_channels_preserved": input_bads,
                    "misc_channels_preserved": misc_after,
                    "execution_status": "completed",
                },
            )
        )
    return run_reports, lineage_records


def save_p10_slow_qc(before, after, record, figure_dir):
    """Expose P10's open slow-contamination concern on unfiltered native signals.

    Inspect up to three evenly spaced 60 s source windows in every acquisition,
    including annotated activity. The matched samples are display/QC only; no
    detector, detrending, high-pass, selection threshold or retention decision
    is introduced. Full-window Welch spectra resolve approximately 0.33 Hz.
    BAD overlap and exact native bounds remain explicit for interpretation.
    """
    sfreq = before.info["sfreq"]
    width = min(before.n_times, int(round(60 * sfreq)))
    starts = sorted(set(int(value) for value in np.linspace(0, before.n_times - width, 3)))
    names = [before.ch_names[i] for i in mne.pick_types(before.info, meg="mag", exclude="bads")]
    if not names:
        return {"status": "missing_good_magnetometers", "decision": "open"}
    frontal = [ch for ch in names if ch.startswith(("MEG011", "MEG012", "MEG031", "MEG041"))][:4] or names[:4]
    fig, axes = plt.subplots(len(starts), 2, figsize=(14, 4 * len(starts)), squeeze=False)
    windows = []
    for row, start in enumerate(starts):
        stop = start + width
        times = (np.arange(start, stop) + before.first_samp) / sfreq
        for raw, color, label in ((before, "0.55", "before"), (after, "tab:blue", "after ICA, before line removal")):
            for index, ch in enumerate(frontal):
                values = raw.get_data(picks=[ch], start=start, stop=stop)[0] * 1e12
                axes[row, 0].plot(times, values + index * 10, color=color, alpha=.7,
                                  label=label if index == 0 else None)
            preview = raw.copy().crop(start / sfreq, (stop - 1) / sfreq)
            spectrum = preview.compute_psd(method="welch", picks=names, fmin=.1, fmax=2.,
                n_fft=width, n_per_seg=width, reject_by_annotation=False, verbose="ERROR")
            power, frequencies = spectrum.get_data(return_freqs=True)
            axes[row, 1].plot(frequencies, 10 * np.log10(np.maximum(np.mean(power, axis=0), np.finfo(float).tiny)),
                              color=color, label=label)
        overlapping = [str(a["description"]) for a in before.annotations
            if a["description"].lower().startswith("bad") and a["onset"] < times[-1] and a["onset"] + a["duration"] > times[0]]
        windows.append({"start_sample_raw": before.first_samp + start,
            "stop_sample_raw_exclusive": before.first_samp + stop, "bad_categories": sorted(set(overlapping))})
        axes[row, 0].set(xlabel="Native acquisition seconds", ylabel="Frontal mag (pT), offsets 10 pT")
        axes[row, 1].set(xlabel="Frequency (Hz)", ylabel="Mean good-mag PSD (dB)")
        axes[row, 1].axvline(1 / 3, color="tab:red", linestyle=":", label="approximately 0.33 Hz")
        for ax in axes[row]:
            ax.legend(fontsize=7)
    fig.suptitle(record.source_filename + " | P10 slow contamination: matched native before/after ICA; decision OPEN")
    fig.tight_layout()
    path = figure_dir / f"03_C_ICA_apply_{record.source_stem}_p10_slow_qc.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return {"status": "evidence_for_explicit_review", "decision": "open", "figure": str(path),
        "windows": windows, "interpretation": "Matched native signals; no high-pass/detrending. Assess residual contamination and task-signal preservation."}


def save_native_comparison(before, after, record, figure_dir):
    """Plot identical native windows before/after ICA, with MISC context.

    This comparison is made before final high-pass/line processing so a change
    is attributable to the selected ICA exclusions. Ocular annotations guide
    the artifact example; the longest ordinary stretch supplies a control.
    """
    from continuous_cleaning import good_stretches
    sfreq = before.info["sfreq"]
    ocular = [a for a in before.annotations if a["description"].startswith(("BAD_blink", "BAD_eye_movement"))]
    spans = good_stretches(before)
    centre = sum(max(spans, key=lambda p: p[1] - p[0])) / (2 * sfreq) if spans else before.times[-1] / 2
    starts = [max(0., min(centre - 2.5, before.times[-1] - 5))]
    if ocular:
        starts.append(max(0., min(ocular[0]["onset"] - before.first_time - 2., before.times[-1] - 5)))
    picks = [ch for ch, kind in zip(before.ch_names, before.get_channel_types())
             if kind == "mag" and ch.startswith(("MEG011", "MEG012", "MEG031", "MEG041"))][:4]
    if not picks:
        picks = [before.ch_names[i] for i in mne.pick_types(before.info, meg="mag", exclude="bads")[:4]]
    misc = [ch for ch, kind in zip(before.ch_names, before.get_channel_types()) if kind == "misc" and ch != "MISC008"]
    fig, axes = plt.subplots(len(starts), 2, squeeze=False, figsize=(13, 4 * len(starts)))
    for row, start in enumerate(starts):
        a, b = int(start * sfreq), min(before.n_times, int((start + 5) * sfreq))
        times = (np.arange(a, b) + before.first_samp) / sfreq
        for index, ch in enumerate(picks):
            for raw, color, label in [(before, "0.6", "before"), (after, "tab:blue", "after ICA")]:
                axes[row, 0].plot(times, raw.get_data(picks=[ch], start=a, stop=b)[0] * 1e12 + index * 5,
                                  color=color, alpha=.8, label=label if index == 0 else None)
        axes[row, 0].set_ylabel("Frontal mag (pT), offset 5 pT"); axes[row, 0].legend()
        for ch in misc:
            axes[row, 1].plot(times, before.get_data(picks=[ch], start=a, stop=b)[0], label=ch)
        axes[row, 1].set_ylabel("MISC acquisition values; assess usability");
        if misc: axes[row, 1].legend(fontsize=7)
        for ax in axes[row]: ax.set_xlabel("Acquisition seconds (absolute sample / sfreq)")
    fig.suptitle(record.source_filename + " | ICA-only comparison, before final temporal filtering")
    fig.tight_layout()
    path = figure_dir / f"03_C_ICA_apply_{record.source_stem}_before_after.png"
    fig.savefig(path, dpi=130); plt.close(fig)
    return path


def save_component_context(model, raw, record, figure_dir, *, raw_root):
    """Save a labelled native example and its deterministic window-selection record.

    Final ICA is projected without removal or temporal filtering. BAD categories,
    exact native bounds, first-ocular/fallback selection and per-window z-score
    limitations are persisted separately from manual/ICA decision provenance.
    """
    path, metadata = native_context(model, raw, record, figure_dir)
    write_json(path.with_suffix(".json"), metadata, raw_root=raw_root)
    return path


def spectral_snapshot(raw, window=None):
    """Summarize mag/grad PSDs on one exact good window, at most 60 seconds."""
    from continuous_cleaning import good_stretches
    sfreq = raw.info["sfreq"]
    if window is None:
        spans = good_stretches(raw)
        if not spans:
            raise RuntimeError("No ordinary window for final-filter spectral inspection.")
        a, b = max(spans, key=lambda pair: pair[1] - pair[0])
        width = min(b - a, int(60 * sfreq))
        a = (a + b - width) // 2
        window = (a, a + width)
    a, b = window
    preview = raw.copy().crop(a / sfreq, (b - 1) / sfreq)
    result = {"window": window, "types": {}}
    for kind in ("mag", "grad"):
        names = [raw.ch_names[i] for i in mne.pick_types(raw.info, meg=kind, exclude="bads")]
        if not names:
            continue
        spectrum = preview.compute_psd(method="welch", picks=names, fmin=0.1,
            fmax=min(150., sfreq / 2), n_fft=min(8192, b - a), reject_by_annotation=False, verbose="ERROR")
        values, frequencies = spectrum.get_data(return_freqs=True)
        result["types"][kind] = {"names": names, "frequencies": frequencies,
                                  "mean_psd": np.mean(values, axis=0)}
    return result


def save_final_spectra(before, raw, record, figure_dir):
    """Compare final temporal processing on identical sensors/native samples.

    This bounded PSD resolves 50/100 Hz treatment, not the 0.01 Hz transition;
    the latter also requires filter-response and edge-loss review.
    """
    after = spectral_snapshot(raw, before["window"])
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for ax, kind in zip(axes, ("mag", "grad")):
        for snapshot, label in [(before, "ICA, before temporal filters"), (after, "Final continuous")]:
            if kind in snapshot["types"]:
                row = snapshot["types"][kind]
                ax.plot(row["frequencies"], 10 * np.log10(np.maximum(row["mean_psd"], np.finfo(float).tiny)), label=label)
        ax.set_ylabel(kind + " PSD (dB)")
        if kind in before["types"]:
            ax.legend(fontsize=8)
        for frequency in (50, 100):
            ax.axvline(frequency, color="0.7", linestyle=":")
    axes[-1].set_xlabel("Frequency (Hz)")
    fig.suptitle(record.source_filename + " | final filtering; identical native window " + str(before["window"]))
    fig.tight_layout()
    path = figure_dir / f"03_C_ICA_apply_{record.source_stem}_final_filter_psd.png"
    fig.savefig(path, dpi=130); plt.close(fig)
    return path


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the production manual ICA CLI."""
    parser = argparse.ArgumentParser(
        description="Review or apply manual ICA for DAD MEG."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--source-manifest", default=None)
    parser.add_argument(
        "--processing-pass",
        choices=("review", "apply"),
        default="review",
    )
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--no-display", action="store_true")
    parser.add_argument("--decisions-out", default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Validate and skip existing apply outputs; fail on any mismatch. "
            "Mutually exclusive with --overwrite."
        ),
    )
    parser.add_argument(
        "--skip-component-properties",
        action="store_true",
        help="Skip expensive per-component property plots during smoke tests.",
    )
    parser.add_argument("--print-summary", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run ICA review or explicit apply mode."""
    args = build_arg_parser().parse_args(argv)
    try:
        cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
            args.config,
            args.subject_config,
            args.tuning_config,
            site_config_path=args.site_config,
            required_sections=("ica",),
        )
        ica_cfg = cfg["ica"]
        if not ica_cfg["enabled"]:
            raise RuntimeError("Manual ICA is disabled in configuration.")
        raw_root = cfg["resolved_paths"]["data_root"]
        overwrite = bool(args.overwrite or cfg["runtime"]["overwrite_existing"])
        if args.resume and overwrite:
            raise ValueError("--resume and --overwrite are mutually exclusive.")
        output_dirs = stage_output_dirs(cfg)
        sss_dir = ensure_output_directory(output_dirs["sss"], raw_root)
        figure_dir = ensure_output_directory(output_dirs["figures"], raw_root)
        model_dir = ensure_output_directory(output_dirs["ica_models"], raw_root)
        report_dir = ensure_output_directory(output_dirs["reports"], raw_root)
        subject_label = cfg["layout"]["subject_bids_id"]
        session_label = cfg["layout"]["session_bids_id"]
        decisions_out = (
            Path(args.decisions_out).expanduser().resolve()
            if args.decisions_out
            else report_dir / "03_tuning_choices.yaml"
        )

        records = discover_raw_run_records(cfg)
        require_review(cfg, "03_A_MaxFilt", report_dir)
        require_review(cfg, "03_B_ArtefactAnnotation", report_dir)
        if args.processing_pass == "apply":
            require_review(cfg, "03_C_ICA", report_dir)
        receipt = cfg.get("stage4_review", {}).get("03_C_ICA", {})
        processing_provenance = (dict(receipt) if receipt.get("acceptance_kind") ==
            "user_provisional_component_processing_authorization" else None)
        input_paths = [
            derivative_path_for_run(record, sss_dir, INPUT_SUFFIX)
            for record in records
        ]
        missing_inputs = [path for path in input_paths if not path.is_file()]
        if missing_inputs:
            raise FileNotFoundError(
                "Expected lineage-matched annotated derivatives are missing: "
                f"{missing_inputs}. Run 03_B_ArtefactAnnotation.py first."
            )

        prefit_path = (
            model_dir / f"{subject_label}_{session_label}_ica_prefit.fif"
        )
        prefit_metadata_path = (
            model_dir
            / f"{subject_label}_{session_label}_ica_prefit_metadata.json"
        )
        postreject_path = (
            model_dir / f"{subject_label}_{session_label}_ica_postreject.fif"
        )
        require_fit_data = args.processing_pass == "review"
        if args.processing_pass == "apply" and not cache_is_compatible(
                prefit_metadata_path, build_ica_fingerprint(input_paths=input_paths, ica_cfg=ica_cfg)):
            raise RuntimeError("Apply cannot refit ICA: regenerate and review a model for changed inputs first.")
        model, fit_data, fit_summaries, cache_status, fingerprint = (
            fit_or_load_ica(
                input_paths=input_paths,
                ica_cfg=ica_cfg,
                prefit_path=prefit_path,
                metadata_path=prefit_metadata_path,
                raw_root=raw_root,
                overwrite=overwrite,
                require_fit_data=require_fit_data,
            )
        )

        eog_diagnostics = {"status": "legacy_not_in_stage4", "per_file": [], "aggregate_ranking": []}

        diagnostic_images: list[Path] = []
        if args.processing_pass == "review":
            if fit_data is None:
                raise RuntimeError("ICA review requires prepared fit data.")
            diagnostic_images = generate_component_diagnostics(
                model=model,
                fit_data=fit_data,
                figure_dir=figure_dir,
                subject_label=subject_label,
                processing_pass=args.processing_pass,
                raw_root=raw_root,
                show=bool(args.interactive and not args.no_display),
                include_properties=not args.skip_component_properties,
            )

        if args.processing_pass == "review":
            for record, input_path in zip(records, input_paths):
                native = mne.io.read_raw_fif(input_path, preload=False, verbose="ERROR")
                diagnostic_images.append(save_component_context(model, native, record, figure_dir, raw_root=raw_root))

        exclusions = ica_cfg["manual_exclude"]
        if args.processing_pass == "review" and ica_cfg.get("reviewed_model_sha256") != file_hash(prefit_path):
            exclusions = None  # Old component indices never accept a new decomposition.
        if exclusions is not None:
            exclusions = sorted(set(int(value) for value in exclusions))
        if args.interactive:
            exclusions = interactive_component_review(
                n_components=int(model.n_components_),
                initial_exclusions=exclusions,
                diagnostic_paths=diagnostic_images,
                eog_diagnostics=eog_diagnostics,
            )
            update_tuning_file(
                path=decisions_out,
                raw_root=raw_root,
                subject_id=cfg["subject"]["id"],
                section="ica",
                values={"manual_exclude": exclusions, "reviewed_model_sha256": file_hash(prefit_path)},
            )

        run_reports: list[dict[str, Any]] = []
        lineage_records: list[dict[str, Any]] = []
        if args.processing_pass == "review":
            status = (
                "review_required" if exclusions is None else "review_complete"
            )
            lineage_records = [
                lineage_record(
                    source=record,
                    input_path=input_path,
                    output_path=None,
                    status=status,
                    metadata={
                        "prefit_ica_model": str(prefit_path),
                        "manual_exclude": exclusions,
                    },
                )
                for record, input_path in zip(records, input_paths, strict=True)
            ]
        else:
            if exclusions is None:
                raise RuntimeError(
                    "ICA apply pass requires an explicit `ica.manual_exclude` "
                    "list. Use [] after reviewing when no components should be "
                    "rejected."
                )
            invalid = [
                value
                for value in exclusions
                if value < 0 or value >= int(model.n_components_)
            ]
            if invalid:
                raise ValueError(
                    "Configured ICA exclusions are outside the fitted component "
                    f"range 0-{int(model.n_components_) - 1}: {invalid}"
                )
            if ica_cfg.get("reviewed_model_sha256") != file_hash(prefit_path):
                raise RuntimeError("ICA exclusions are unbound or belong to a different prefit model.")
            model.exclude = list(exclusions)
            if args.resume and postreject_path.is_file():
                prepared_postreject = guarded_path(postreject_path, raw_root)
                existing_postreject = read_ica(
                    prepared_postreject,
                    verbose="WARNING",
                )
                same_unmixing = np.allclose(
                    existing_postreject.unmixing_matrix_,
                    model.unmixing_matrix_,
                )
                if (
                    list(existing_postreject.exclude) != exclusions
                    or int(existing_postreject.n_components_)
                    != int(model.n_components_)
                    or not same_unmixing
                    or not np.allclose(existing_postreject.pca_components_, model.pca_components_)
                    or not np.allclose(existing_postreject.pre_whitener_, model.pre_whitener_)
                    or not np.allclose(existing_postreject.pca_mean_, model.pca_mean_)
                ):
                    raise RuntimeError(
                        "Cannot resume from a post-rejection ICA model that "
                        "does not match the current prefit model and exclusions: "
                        f"{prepared_postreject}"
                    )
            else:
                prepared_postreject = prepare_output_path(
                    postreject_path,
                    raw_root=raw_root,
                    overwrite=overwrite,
                )
                model.save(
                    prepared_postreject,
                    overwrite=overwrite,
                    verbose="WARNING",
                )
            run_reports, lineage_records = apply_ica_to_runs(
                model=model,
                exclusions=exclusions,
                records=records,
                input_paths=input_paths,
                sss_dir=sss_dir,
                raw_root=raw_root,
                overwrite=overwrite,
                resume=bool(args.resume), final_settings=cfg["preprocessing"],
                model_sha256=file_hash(prefit_path), figure_dir=figure_dir,
            )
            status = "completed"

        report_payload = {
            "stage": STAGE_NAME,
            "generated_utc": utc_now(),
            "processing_pass": args.processing_pass,
            "interactive": bool(args.interactive),
            "resume": bool(args.resume),
            "subject_bids_id": subject_label,
            "status": status,
            "warnings": warnings,
            "cache_status": cache_status,
            "fingerprint": fingerprint,
            "fit_input_summaries": fit_summaries,
            "n_components_fitted": int(model.n_components_),
            "fit_execution": {"n_samples_fitted": int(model.n_samples_),
                "n_iterations": int(model.n_iter_), "max_iter": iteration_cap(model),
                "iteration_limit_reached": bool(model.n_iter_ >= iteration_cap(model)),
                "mapped_fit_samples": sum(row["fit_samples"] for row in fit_summaries),
                "pca_explained_variance": model.pca_explained_variance_.tolist()},
            "manual_exclude": exclusions,
            "scientific_acceptance": "separate_explicit_stage4_review",
            "component_properties_requested": not args.skip_component_properties,
            "review_interface_version": 2,
            "diagnostic_scope": "full_fit",
            "prefit_model": str(prefit_path),
            "prefit_model_sha256": file_hash(prefit_path),
            "postreject_model": (
                str(postreject_path)
                if args.processing_pass == "apply"
                else None
            ),
            "runs": run_reports,
            "blink_annotation_policy": "No ICA-derived blink annotations added.",
        }
        report_payload["configuration_sha256"] = stage_configuration(cfg, "03_C_ICA", include_edge=False)
        if processing_provenance:
            report_payload["processing_provenance"] = processing_provenance
        # Persist the real optimization trajectory, including a deterministic
        # restart's separate phase count. It is not a fresh seed-97 random fit.
        optimizer_restart = json.loads(prefit_metadata_path.read_text()).get("optimizer_restart")
        if optimizer_restart:
            report_payload["optimizer_restart"] = optimizer_restart
            report_payload["fit_execution"]["iteration_count_scope"] = "restart_phase_only"
        native_manifest_paths = []
        if args.processing_pass == "review":
            from review_final_ica import make_plot_recipes, settle_review_outputs
            if not args.skip_component_properties:
                # Reuse one native projection/epoch extraction for ALL final
                # components, keeping continuous and task evidence visible.
                from final_ica_native_diagnostics import generate_native_diagnostics
                native_images, native_manifest = generate_native_diagnostics(cfg, model,
                    records, input_paths, figure_dir, file_hash(prefit_path))
                diagnostic_images.extend(native_images)
                report_payload["review_interface_version"] = 3
                report_payload["native_component_diagnostics"] = {
                    "manifest": "03_C_ICA_review_native_diagnostics.json",
                    "components_complete": native_manifest["components_complete"],
                    "continuous_and_task_default": True}
                native_manifest_paths = [figure_dir / "03_C_ICA_review_native_diagnostics.json"]
            report_payload["plot_recipes"] = make_plot_recipes(cfg, report_payload, diagnostic_images)
            book = report_dir / "03_C_ICA_review_components.html"
            report_payload["component_review_coverage"] = write_component_book(
                report_payload, diagnostic_images, book, raw_root)
            report_payload["component_review_book"] = str(book)
            settle_review_outputs(diagnostic_images + native_manifest_paths + [book])
            report_payload["component_book_sha256"] = file_hash(book)
        diagnostic_evidence = diagnostic_images + native_manifest_paths + [book] if args.processing_pass == "review" else []
        bind_report(report_payload, input_paths + [prefit_path, prefit_metadata_path] + diagnostic_evidence,
                    {"fingerprint": fingerprint, "model_sha256": file_hash(prefit_path)},
                    report_dir, cfg["resolved_paths"]["derivatives_root"])
        report_json = write_json(
            report_dir / f"03_C_ICA_{args.processing_pass}_report.json",
            report_payload,
            raw_root=raw_root,
        )
        # Render saved JSON/PNGs into the primary self-contained report. This
        # presentation step opens no FIFs and makes no scientific decisions.
        report_html = render_stage_report(
            derivatives_root=Path(cfg["resolved_paths"]["derivatives_root"]),
            raw_root=raw_root,
            subject_id=cfg["layout"]["subject_bids_id"],
            processing_pass=args.processing_pass,
            stage_name=STAGE_NAME,
            overwrite=True,
            primary_report=True,
        )
        manifest_status = update_source_manifest_stage(
            manifest_path=args.source_manifest,
            raw_root=raw_root,
            stage_name=STAGE_NAME,
            processing_pass=args.processing_pass,
            records=lineage_records,
            status=status,
            qc={
                "report_json": str(report_json),
                "report_html": str(report_html),
                "manual_exclude": exclusions,
                "cache_status": cache_status,
            },
        )
        if args.print_summary:
            print(f"Manual ICA: {status}")
            print(f"Cache: {cache_status}")
            print(f"Components: {int(model.n_components_)}")
            print(f"Manual exclusions: {exclusions}")
            print(f"JSON report: {report_json}")
            print(f"HTML report: {report_html}")
            print(f"Source manifest: {manifest_status}")
            if args.interactive:
                print(f"Accepted tuning: {decisions_out}")
        return 0
    except Exception as exc:  # pylint: disable=broad-except
        print(f"03_C ICA error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
