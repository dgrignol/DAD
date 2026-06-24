#!/usr/bin/env python3
"""
03_C_ICA.py

Purpose
Fit, review, and apply the production manual ICA branch for DAD MEG data.
Component rejection is always an explicit subject-level decision. Batch mode
never prompts and never invents exclusions from EOG scores or component index.

Usage examples
1) First-pass review using FLUX-aligned fit defaults:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/03_C_ICA.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-01.yaml \
     --processing-pass review --print-summary

   This fits or reuses a compatible ICA cache, saves topographies and component
   property plots, records EOG diagnostics when channels exist, and exits
   successfully with ``review_required`` when ``ica.manual_exclude`` is null.

2) Interactive component review:

   python .../03_C_ICA.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --processing-pass review \
     --interactive --decisions-out .../03_tuning_choices.yaml

   Component diagnostics are saved and topography pages are optionally
   displayed. The operator enters component indices, reviews the proposed list,
   and repeats until accepting it.

3) Apply reviewed exclusions:

   python .../03_C_ICA.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --tuning-config .../03_tuning_choices.yaml \
     --processing-pass apply --overwrite

4) Resume an interrupted apply pass. Existing cleaned runs and the
   post-rejection model are validated before reuse:

   python .../03_C_ICA.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --tuning-config .../03_tuning_choices.yaml \
     --processing-pass apply --resume

Inputs
- Layered configuration; ``--site-config`` overrides machine paths without
  changing scientific settings.
- ``<exact-source-stem>_ann.fif`` from 03_B for every acquisition run.

Outputs
- Review: compatible prefit ICA model, metadata fingerprint, component PNGs,
  and pass-specific JSON/HTML QC reports.
- Apply: post-rejection ICA model and
  ``<exact-source-stem>_ica.fif`` cleaned acquisition runs.
- Optional accepted tuning YAML and source-manifest derivative lineage.

Data flow
annotated acquisition runs -> MEG-only 200 Hz, 1-40 Hz fit copies ->
concatenated FastICA decomposition -> component/EOG diagnostics -> explicit
manual exclusion decision -> projection onto each original annotated run.

Scientific assumptions
- FLUX defaults are retained where applicable: FastICA, 30 components,
  ``random_state=97``, 200 Hz resampling, and 1-40 Hz fit filtering.
- Different ``dev_head_t`` values across acquisition runs are accepted with an
  explicit warning during concatenation, matching the documented local
  deviation.
- EOG scores are diagnostics only. Missing EOG channels skip cleanly.
- No ICA000 blink workaround is used.
- A null exclusion list means unreviewed. ``[]`` means reviewed and rejecting
  no components.

Safety
Every output passes the shared raw-tree guard. Existing models and cleaned FIFs
require explicit overwrite permission unless a prefit cache fingerprint proves
that the existing cache is compatible and can be reused read-only. Resume also
requires matching source/stage lineage, exclusion list, channel identity, bad
channels, sample count, and MISC channels.
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
    write_simple_html_report,
)


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
        "schema_version": 2,
        "inputs": [
            {
                "source_filename": path.name,
                "size_bytes": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
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


def prepare_ica_fit_data(
    input_paths: list[Path],
    ica_cfg: Mapping[str, Any],
) -> tuple[mne.io.BaseRaw, list[dict[str, Any]]]:
    """
    Build and concatenate the MEG-only ICA fit copies.

    MISC, EOG, and stimulation channels remain in the original annotated files
    and are not part of the ICA decomposition.
    """
    prepared: list[mne.io.BaseRaw] = []
    summaries: list[dict[str, Any]] = []
    for input_path in input_paths:
        raw = mne.io.read_raw_fif(
            input_path,
            preload=True,
            verbose="WARNING",
        )
        fit_raw = raw.copy().pick("meg")
        fit_raw.resample(float(ica_cfg["fit_resample_sfreq_hz"]))
        fit_raw.filter(
            float(ica_cfg["fit_l_freq_hz"]),
            float(ica_cfg["fit_h_freq_hz"]),
            verbose="WARNING",
        )
        summaries.append(
            {
                "input_path": str(input_path),
                "n_channels": fit_raw.info["nchan"],
                "n_times": fit_raw.n_times,
                "sfreq_hz": float(fit_raw.info["sfreq"]),
                "bad_channels": list(fit_raw.info.get("bads", [])),
            }
        )
        prepared.append(fit_raw)
    concatenated = mne.concatenate_raws(
        prepared,
        on_mismatch="warn",
        verbose="WARNING",
    )
    return concatenated, summaries


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
    compatible model without concatenating all fit copies again.
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
    if can_reuse:
        model = read_ica(prefit_path, verbose="WARNING")
        cache_status = "reused_compatible"
        if require_fit_data:
            fit_data, fit_summaries = prepare_ica_fit_data(input_paths, ica_cfg)
        return model, fit_data, fit_summaries, cache_status, fingerprint

    fit_data, fit_summaries = prepare_ica_fit_data(input_paths, ica_cfg)
    model = ICA(
        method=str(ica_cfg["method"]),
        n_components=ica_cfg["n_components"],
        random_state=int(ica_cfg["random_state"]),
        max_iter=ica_cfg["max_iter"],
        verbose="WARNING",
    )
    model.fit(fit_data, verbose="WARNING")
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
    if show:
        source_fig = model.plot_sources(
            fit_data,
            title=f"{subject_label} ICA sources",
            show=True,
        )
    else:
        with mne.viz.use_browser_backend("matplotlib"):
            source_fig = model.plot_sources(
                fit_data,
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
        for component in range(int(model.n_components_)):
            try:
                property_figures = model.plot_properties(
                    fit_data,
                    picks=component,
                    psd_args={"fmax": float(fit_data.info["lowpass"])},
                    show=False,
                    verbose=False,
                )
                if not isinstance(property_figures, list):
                    property_figures = [property_figures]
                for page, figure in enumerate(property_figures, start=1):
                    path = prepare_output_path(
                        figure_dir
                        / (
                            f"03_C_ICA_{processing_pass}_ICA{component:03d}_"
                            f"properties_{page:02d}.png"
                        ),
                        raw_root=raw_root,
                        overwrite=True,
                    )
                    images.append(_save_figure(figure, path, show=False))
            except Exception as exc:  # pylint: disable=broad-except
                print(
                    f"WARNING: ICA{component:03d} property plot failed: {exc}",
                    file=sys.stderr,
                )
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
        raw.save(output_path, overwrite=overwrite, verbose="WARNING")
        run_report = {
            "source_filename": record.source_filename,
            "acquisition_run_label": record.acquisition_run_label,
            "execution_status": "completed",
            "input_path": str(input_path),
            "output_path": str(output_path),
            "excluded_components": exclusions,
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

        eog_diagnostics = (
            summarize_eog_scores(
                model=model,
                input_paths=input_paths,
                configured_channels=list(ica_cfg["eog_channel_names"]),
            )
            if ica_cfg["eog_diagnostics"]
            else {
                "status": "disabled",
                "per_file": [],
                "aggregate_ranking": [],
            }
        )

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

        exclusions = ica_cfg["manual_exclude"]
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
                values={"manual_exclude": exclusions},
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
                resume=bool(args.resume),
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
            "manual_exclude": exclusions,
            "eog_diagnostics": eog_diagnostics,
            "prefit_model": str(prefit_path),
            "postreject_model": (
                str(postreject_path)
                if args.processing_pass == "apply"
                else None
            ),
            "runs": run_reports,
            "blink_annotation_policy": "No ICA-derived blink annotations added.",
        }
        report_json = write_json(
            report_dir / f"03_C_ICA_{args.processing_pass}_report.json",
            report_payload,
            raw_root=raw_root,
        )
        summary_lines = [
            f"status: {status}",
            f"processing_pass: {args.processing_pass}",
            f"subject: {subject_label}",
            f"cache_status: {cache_status}",
            f"n_components_fitted: {int(model.n_components_)}",
            f"manual_exclude: {exclusions}",
            f"EOG diagnostics: {eog_diagnostics['status']}",
            "ICA-derived blink annotations: disabled",
            f"prefit_model: {prefit_path}",
        ]
        if args.processing_pass == "apply":
            summary_lines.append(f"postreject_model: {postreject_path}")
            summary_lines.extend(
                (
                    f"run{item['acquisition_run_label']:02d} "
                    f"{item['source_filename']} "
                    f"({item['execution_status']}) -> {item['output_path']}"
                )
                for item in run_reports
            )
        report_html = write_simple_html_report(
            report_dir / f"03_C_ICA_{args.processing_pass}_report.html",
            title=f"{subject_label} | 03_C Manual ICA ({args.processing_pass})",
            summary_lines=summary_lines,
            image_paths=diagnostic_images,
            raw_root=raw_root,
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
