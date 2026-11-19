#!/usr/bin/env python3
"""
07_run_matlab_drsa.py

Processing runs
Export the same DAD_EPOCH_PREPARATION and DAD_PROCESSING_RUN values used for
Stage 04. Analysis/reports use the matching processing_runs/LABEL directory.
Labelled resume verifies input/settings contracts and artifact bindings.

Purpose
Run the optional MATLAB dRSA analysis layer from the same resolved
configuration and lineage contracts used by the Python MEG pipeline. This
stage consumes the completed Stage-06 identity-averaged neural exports and
creates a small MATLAB handoff JSON so the MATLAB scripts do not reconstruct
paths, subject versions, identity axes, model paths, or time axes themselves.

Usage examples
1) Validate handoff inputs for one local subject without launching MATLAB:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/07_run_matlab_drsa.py \
     --config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/pipeline_default.yaml \
     --subject-config /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/params/sub-04.yaml \
     --source-manifest /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/reports/sub-04/ses-01/source_manifest.json \
     --validate-only --print-summary

2) Run the default real-neural dRSA MATLAB script:

   python .../07_run_matlab_drsa.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --source-manifest .../source_manifest.json \
     --analyses real --resume --print-summary

3) Run sensor-group dRSA and also write smoothed display figures:

   python .../07_run_matlab_drsa.py --config .../pipeline_default.yaml \
     --subject-config .../sub-04.yaml --source-manifest .../source_manifest.json \
     --analyses alt --sensor-groups occipital --smoothing-plots --resume

4) Run observed, expected/predicted, and PE occluded-deviant alt stimulus
   models:

   python .../07_run_matlab_drsa.py --config .../pipeline_default.yaml \
     --subject-config .../sub-04.yaml --source-manifest .../source_manifest.json \
     --analyses alt --conditions occluded_deviant \
     --sensor-groups occipital \
     --model-sources observed predicted prediction_error pe_prior_constrained_50 \
     --overwrite

5) Build a numeric group-average figure report after subject dRSA has run:

   python .../07_run_matlab_drsa.py --config .../pipeline_default.yaml \
     --subject-config .../sub-01.yaml --group-average-only \
     --group-subjects {04..39} --smoothing-plots

6) Run from the cluster through site paths and an explicit MATLAB binary:

   DAD_MATLAB_BIN=/usr/local/MATLAB/R2020a/bin/matlab DAD_MATLAB_MAX_THREADS=4 \
   python .../07_run_matlab_drsa.py --config .../pipeline_default.yaml \
     --site-config .../pipeline_cluster.yaml --subject-config .../sub-04.yaml \
     --source-manifest .../source_manifest.json --analyses real --resume

Cluster wrappers set DAD_MATLAB_MAX_THREADS to the allocated CPU count; generated
batches enforce it with maxNumCompThreads before analysis. An explicit limit
must be positive and cannot exceed SLURM_CPUS_PER_TASK. This affects resources,
not the dRSA model, identity, smoothing or statistical settings.

Inputs
- Completed Stage-06 condition MAT files and JSON sidecars under the configured
  derivative root.
- Completed source_manifest.json containing export_drsa lineage.
- One subject-specific observed trajectory MAT discovered from the resolved
  subject input_files directory. Its matching ``occluded_nondeviant`` trials
  provide the expected continuation for ``occluded_deviant`` and the baseline
  for ``prediction_error`` and explicit ``pe_*`` model-source runs.
- MATLAB scripts in meg/code/analysis plus dRSA helper functions in
  simulations/functions.

Outputs
- reports/sub-XX/ses-01/07_matlab_drsa_handoff.json
- reports/sub-XX/ses-01/07_matlab_drsa_report.json
- reports/sub-XX/ses-01/07_matlab_drsa_report.html
- reports/sub-XX/ses-01/07_matlab_drsa_report_embedded.html
- MATLAB PNG/MAT outputs below analysis/sub-XX/ses-01/drsa/
- Optional group-average PNG/JSON/HTML outputs below analysis/group/ses-01/drsa/
  and reports/group/ses-01/. Group HTML reports also receive sibling
  *_embedded.html copies with all image bytes inlined for sharing.
- export_drsa remains unchanged; this stage records matlab_drsa lineage only.

Safety and assumptions
- The raw data tree is never modified. All generated paths pass the shared
  raw-tree guard.
- MATLAB is optional and independently rerunnable. With --resume, existing
  expected outputs skip recomputation for the selected analyses.
- MATLAB receives a handoff JSON with already validated paths and metadata.
  The MATLAB scripts remain responsible for validating tensors and model trial
  order before computing dRSA.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import mimetypes
import os
import re
import shlex
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from downstream_stage_utils import discover_trajectory_models
from continuous_cleaning import digest, file_hash, file_stamp
from epoch_input_preparation import processing_manifest, verify_export_artifacts
from drsa_result_utils import (
    analysis_family_from_path,
    model_source_from_path,
    model_time_axis_label,
    smooth_matrix_for_plot as shared_smooth_matrix_for_plot,
    stimulus_model_label,
)
from pipeline_config import (
    assert_output_path_outside_raw,
    load_and_resolve_config,
)
from preprocessing_stage_utils import (
    guarded_path,
    update_source_manifest_stage,
    write_json,
    write_simple_html_report,
)


CONDITIONS = ("always_visible", "occluded_nondeviant", "occluded_deviant")
SUPPORTED_ANALYSES = ("real", "alt", "simulations")
PE_MODEL_SOURCES = (
    "prediction_error",
    "pe_prior_unconstrained_50",
    "pe_prior_constrained_50",
    "pe_precision_weighted",
    "pe_kalman_innovation",
    "pe_latent_curvature",
    "pe_belief_update",
)
MODEL_SOURCES = ("observed", "predicted", *PE_MODEL_SOURCES)
DRSA_MATRIX_VARIABLES = (
    "dRSA_neural_position",
    "dRSA_neural_correlation",
    "dRSA_neural_self",
    "dRSA_model_pos_pos",
    "dRSA_model_pos_dir",
)


def utc_now() -> str:
    """Return an ISO8601 UTC timestamp with an explicit Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the optional MATLAB dRSA command-line interface."""
    parser = argparse.ArgumentParser(
        description="Run optional MATLAB dRSA from Stage-06 neural exports."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--site-config", default=None)
    parser.add_argument("--subject-config", default=None)
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument(
        "--source-manifest",
        default=None,
        help=(
            "Subject source_manifest.json. Required for subject MATLAB runs; "
            "not required with --group-average-only."
        ),
    )
    parser.add_argument(
        "--analyses",
        nargs="+",
        choices=SUPPORTED_ANALYSES,
        default=None,
        help="MATLAB analysis families to run.",
    )
    parser.add_argument(
        "--drsa-type",
        choices=("corr", "PCR"),
        default=None,
        help="Override matlab_drsa.drsa_type.",
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        choices=CONDITIONS,
        default=None,
        help="Override matlab_drsa.conditions_to_run.",
    )
    parser.add_argument(
        "--sensor-groups",
        nargs="+",
        default=None,
        help="Override matlab_drsa.sensor_groups_to_run for alt/simulations.",
    )
    parser.add_argument(
        "--model-source",
        choices=MODEL_SOURCES,
        default=None,
        help="Override matlab_drsa.model_source.",
    )
    parser.add_argument(
        "--model-sources",
        nargs="+",
        choices=MODEL_SOURCES,
        default=None,
        help=(
            "Run multiple stimulus model sources. Predicted and "
            "PE model sources are valid only for occluded_deviant. "
            "PE model sources are implemented for alt outputs; simulations "
            "continue to use observed/predicted sources."
        ),
    )
    parser.add_argument(
        "--matlab-bin",
        default=None,
        help=(
            "Absolute MATLAB executable. Overrides DAD_MATLAB_BIN and "
            "matlab_drsa.matlab_bin."
        ),
    )
    parser.add_argument(
        "--use-arch-x86-64",
        action="store_true",
        help="Launch MATLAB through `arch -x86_64` for local R2020a.",
    )
    parser.add_argument(
        "--smoothing-plots",
        dest="smoothing_plots",
        action="store_true",
        default=None,
        help="Also export smoothed display-only PNG variants of dRSA matrices.",
    )
    parser.add_argument(
        "--no-smoothing-plots",
        dest="smoothing_plots",
        action="store_false",
        help="Disable smoothed display PNGs even if enabled in config/env.",
    )
    parser.add_argument(
        "--plot-smoothing-window",
        type=int,
        default=None,
        help="Gaussian smoothing window, in samples, for display PNGs.",
    )
    parser.add_argument(
        "--group-average",
        action="store_true",
        help="After the subject run, build group-average dRSA figures and report.",
    )
    parser.add_argument(
        "--group-average-only",
        action="store_true",
        help="Only build group-average dRSA figures/report from existing subject MAT outputs.",
    )
    parser.add_argument(
        "--group-subjects",
        nargs="+",
        default=None,
        help=(
            "Subjects to include in group averaging, e.g. 01 02 sub-03. "
            "Defaults to config/env or all discovered sub-* folders."
        ),
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--print-summary", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser


def _analysis_paths(cfg: Mapping[str, Any]) -> dict[str, Path]:
    """Return stage input/output paths from the resolved configuration."""
    resolved = cfg["resolved_paths"]
    project_root = Path(resolved["project_root"])
    analysis_dir = Path(resolved["derivatives_analysis_dir"])
    reports_dir = Path(resolved["derivatives_reports_dir"])
    subject = cfg["layout"]["subject_bids_id"]
    session = cfg["layout"]["session_bids_id"]
    group_analysis_dir = analysis_dir.parent.parent / "group" / session
    group_reports_dir = reports_dir.parent.parent / "group" / session
    return {
        "project_root": project_root,
        "dad_root": project_root.parent,
        "analysis_dir": analysis_dir,
        "reports_dir": reports_dir,
        "drsa_dir": analysis_dir / "drsa",
        "group_analysis_dir": group_analysis_dir,
        "group_reports_dir": group_reports_dir,
        "group_drsa_dir": group_analysis_dir / "drsa",
        "group_report_json": group_reports_dir / "07_matlab_drsa_group_report.json",
        "group_report_html": group_reports_dir / "07_matlab_drsa_group_report.html",
        "group_report_smoothed_json": (
            group_reports_dir / "07_matlab_drsa_group_report_smoothed.json"
        ),
        "group_report_smoothed_html": (
            group_reports_dir / "07_matlab_drsa_group_report_smoothed.html"
        ),
        "handoff_json": reports_dir / "07_matlab_drsa_handoff.json",
        "report_json": reports_dir / "07_matlab_drsa_report.json",
        "report_html": reports_dir / "07_matlab_drsa_report.html",
        "real_script": project_root / "code" / "analysis" / "A_dRSA.m",
        "alt_script": project_root / "code" / "analysis" / "A_dRSA_alt.m",
        "simulations_script": (
            project_root / "code" / "analysis" / "dRSA_simulations.m"
        ),
        "drsa_functions_dir": project_root.parent / "simulations" / "functions",
        "epoch_report_json": reports_dir / "04_epoching_report.json",
        "export_report_json": reports_dir / "06_export_neural_for_drsa_report.json",
        "source_manifest": reports_dir / "source_manifest.json",
        "subject": Path(subject),
        "session": Path(session),
    }


def _condition_output_paths(
    cfg: Mapping[str, Any],
    drsa_dir: Path,
    condition: str,
) -> tuple[Path, Path]:
    """Build the Stage-06 MAT/JSON path pair for one condition."""
    subject = cfg["layout"]["subject_bids_id"]
    session = cfg["layout"]["session_bids_id"]
    stem = (
        f"{subject}_{session}_task-Attention_condition-{condition}_"
        "desc-identityAveraged_neural"
    )
    return drsa_dir / f"{stem}.mat", drsa_dir / f"{stem}.json"


def _load_json(path: Path) -> dict[str, Any]:
    """Read one JSON object with a clear error on missing or wrong shape."""
    if not path.is_file():
        raise FileNotFoundError(f"Required JSON file is missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"Expected JSON object in {path}.")
    return payload


def write_embedded_html_copy(
    report_path: str | Path,
    *,
    raw_root: str | Path,
) -> Path:
    """
    Write a self-contained HTML copy with image sources embedded as data URIs.

    Data flow:
    existing report HTML -> resolve absolute/relative image paths -> base64
    inline copies -> sibling ``*_embedded.html`` shareable report.
    """
    report = Path(report_path).expanduser().resolve()
    output = report.with_name(f"{report.stem}_embedded{report.suffix}")
    guarded_path(output, raw_root)
    html_text = report.read_text(encoding="utf-8")

    def replace_src(match: re.Match[str]) -> str:
        quote = match.group(1)
        src = match.group(2)
        if src.startswith("data:"):
            return match.group(0)
        image_path = Path(src).expanduser()
        if not image_path.is_absolute():
            image_path = (report.parent / image_path).resolve()
        if not image_path.is_file():
            return match.group(0)
        mime_type = mimetypes.guess_type(image_path.name)[0] or "image/png"
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        return f"src={quote}data:{mime_type};base64,{encoded}{quote}"

    embedded = re.sub(r"src=(['\"])([^'\"]+)\1", replace_src, html_text)
    output.write_text(embedded, encoding="utf-8")
    return output


def _as_list_of_ints(values: Iterable[Any], label: str) -> list[int]:
    """Normalize JSON numeric values to integers."""
    normalized = [int(value) for value in values]
    if not normalized:
        raise RuntimeError(f"{label} must not be empty.")
    return normalized


def _env_flag(name: str) -> bool | None:
    """Parse an optional boolean environment variable."""
    value = os.environ.get(name)
    if value is None or value == "":
        return None
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y", "on"}:
        return True
    if normalized in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"{name} must be a boolean-like value, got {value!r}.")


def _csv_env_list(name: str) -> list[str] | None:
    """Read a comma-separated environment variable as a string list."""
    value = os.environ.get(name)
    if value is None or value.strip() == "":
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def _normalize_model_sources(values: Iterable[str]) -> list[str]:
    """Normalize model-source selectors while preserving user order."""
    normalized: list[str] = []
    for value in values:
        source = str(value).strip().lower()
        if source == "pe":
            source = "prediction_error"
        if source not in MODEL_SOURCES:
            raise ValueError(
                "model source must be observed, predicted, prediction_error, "
                "or a registered pe_* variant, "
                f"got {value!r}."
            )
        if source not in normalized:
            normalized.append(source)
    if not normalized:
        raise ValueError("At least one model source is required.")
    return normalized


def _model_sources_from_config_or_args(
    args: argparse.Namespace,
    matlab_cfg: Mapping[str, Any],
) -> list[str]:
    """
    Resolve requested observed/predicted/prediction-error stimulus sources.

    Data flow:
    --model-sources or DAD_DRSA_MODEL_SOURCES -> legacy --model-source ->
    matlab_drsa.model_sources/model_source -> normalized ordered source list.
    """
    configured = matlab_cfg.get("model_sources")
    if configured is None:
        configured = matlab_cfg.get("model_source", "observed")
    configured_values: Iterable[str]
    if isinstance(configured, str):
        configured_values = [configured]
    else:
        configured_values = configured
    values = (
        args.model_sources
        or _csv_env_list("DAD_DRSA_MODEL_SOURCES")
        or ([args.model_source] if args.model_source else None)
        or configured_values
    )
    return _normalize_model_sources(values)


def _conditions_for_model_source(
    conditions: Iterable[str],
    model_source: str,
) -> list[str]:
    """
    Return the scientifically valid conditions for one model source.

    The expected and PE models are defined for the
    occluded-deviant question by using matching occluded-nondeviant trials.
    Running them for other conditions would create misleading duplicated
    outputs, so we exclude those here.
    """
    selected = [str(condition) for condition in conditions]
    if model_source == "predicted" or model_source in PE_MODEL_SOURCES:
        return [condition for condition in selected if condition == "occluded_deviant"]
    return selected


def _model_sources_for_analysis(
    analysis: str,
    model_sources: Iterable[str],
) -> list[str]:
    """
    Return model sources supported by one MATLAB analysis family.

    A_dRSA_alt can construct the single-model PE tensors.
    dRSA_simulations still delegates model selection to drsa_stage_utils, which
    only returns observed or expected/predicted trajectory trials.
    """
    selected = list(model_sources)
    if analysis == "simulations":
        return [source for source in selected if source not in PE_MODEL_SOURCES]
    return selected


def _resolve_plot_smoothing(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> tuple[bool, int]:
    """Resolve display-smoothing settings from CLI, env, then config."""
    matlab_cfg = cfg.get("matlab_drsa", {})
    env_enabled = _env_flag("DAD_DRSA_SMOOTHING_PLOTS")
    enabled = (
        args.smoothing_plots
        if args.smoothing_plots is not None
        else (
            env_enabled
            if env_enabled is not None
            else bool(matlab_cfg.get("smoothing_plots", False))
        )
    )
    env_window = os.environ.get("DAD_DRSA_PLOT_SMOOTHING_WINDOW")
    window = int(
        args.plot_smoothing_window
        or (int(env_window) if env_window else matlab_cfg.get("plot_smoothing_window", 9))
    )
    if window < 1:
        raise ValueError("--plot-smoothing-window must be >= 1.")
    return bool(enabled), int(window)


def _normalize_subject_label(value: str) -> str:
    """Normalize 01/sub-01 style input to sub-XX."""
    text = str(value).strip()
    if not text:
        raise ValueError("Group-average subject labels must not be empty.")
    if text.startswith("sub-"):
        return text
    return f"sub-{int(text):02d}" if text.isdigit() else f"sub-{text}"


def _resolve_group_subjects(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> list[str]:
    """Resolve group-average subjects from CLI, env, config, or disk scan."""
    paths = _analysis_paths(cfg)
    session = cfg["layout"]["session_bids_id"]
    matlab_cfg = cfg.get("matlab_drsa", {})
    configured = (
        args.group_subjects
        or _csv_env_list("DAD_DRSA_GROUP_SUBJECTS")
        or matlab_cfg.get("group_average_subjects")
    )
    if configured:
        return [_normalize_subject_label(str(value)) for value in configured]

    analysis_root = paths["analysis_dir"].parent.parent
    discovered = []
    for subject_dir in sorted(analysis_root.glob("sub-*")):
        if (subject_dir / session / "drsa").is_dir():
            discovered.append(subject_dir.name)
    if not discovered:
        raise RuntimeError(
            f"No subject dRSA folders were found under {analysis_root}."
        )
    return discovered


def _occlusion_timing_from_config(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """
    Resolve occlusion boundary timing for Stage-07 figure overlays.

    The stimulus generation code uses 1-based frame indices, while dRSA plots
    label the first sample as t=0. Therefore frame-to-second conversion uses
    ``(frame - 1) / fps``.
    """
    timing_cfg = cfg.get("matlab_drsa", {}).get("occlusion_timing", {})
    fps = float(
        timing_cfg.get(
            "frame_frequency_hz",
            cfg.get("export", {}).get("sampling_rate_hz", 120.0),
        )
    )
    first_complete_frame = int(
        timing_cfg.get("first_complete_occlusion_frame", 130)
    )
    reappearance_frame = int(timing_cfg.get("reappearance_frame", 191))
    return {
        "frame_frequency_hz": fps,
        "first_complete_occlusion_frame": first_complete_frame,
        "reappearance_frame": reappearance_frame,
        "first_complete_occlusion_sec": (first_complete_frame - 1) / fps,
        "reappearance_sec": (reappearance_frame - 1) / fps,
        "labels": {
            "first_complete_occlusion": "first complete occlusion",
            "reappearance": "first reappearance",
        },
    }


def _occlusion_line_times(
    condition: str,
    occlusion_timing: Mapping[str, Any] | None,
) -> list[float]:
    """Return horizontal occlusion boundary times for occluded conditions."""
    if "occluded" not in str(condition).lower() or not occlusion_timing:
        return []
    lines: list[float] = []
    for key in ("first_complete_occlusion_sec", "reappearance_sec"):
        value = occlusion_timing.get(key)
        if value is not None:
            lines.append(float(value))
    return lines


def build_handoff(
    cfg: Mapping[str, Any],
    *,
    source_manifest: str | Path,
    selected_analyses: Iterable[str],
) -> dict[str, Any]:
    """
    Validate Stage-06 outputs and assemble one MATLAB handoff payload.

    Data flow:
    resolved config -> expected Stage-06 files -> sidecar validation ->
    subject trajectory discovery -> handoff JSON consumed by MATLAB.
    """
    paths = _analysis_paths(cfg)
    raw_root = Path(cfg["resolved_paths"]["data_root"])
    manifest_path = processing_manifest(cfg, source_manifest)
    if cfg.get("epoching", {}).get("processing_run"):
        verify_export_artifacts(cfg)
    manifest = _load_json(manifest_path)
    export_lineage = (
        manifest.get("derivative_lineage", {})
        .get("stages", {})
        .get("export_drsa")
    )
    if not isinstance(export_lineage, Mapping) or export_lineage.get("status") != "completed":
        raise RuntimeError(
            "Source manifest must contain completed export_drsa lineage before "
            "MATLAB dRSA can run."
        )

    export_report = _load_json(paths["export_report_json"])
    if export_report.get("status") != "completed":
        raise RuntimeError(
            f"Stage-06 export report is not completed: {paths['export_report_json']}"
        )

    trajectory_models = discover_trajectory_models(cfg)
    condition_payloads: list[dict[str, Any]] = []
    identity_values_ref: list[int] | None = None
    time_axis_ref: list[float] | None = None
    sampling_rate_ref: float | None = None
    channel_labels_ref: list[str] | None = None
    for condition in CONDITIONS:
        neural_mat, sidecar_json = _condition_output_paths(
            cfg,
            paths["drsa_dir"],
            condition,
        )
        guarded_path(neural_mat, raw_root)
        guarded_path(sidecar_json, raw_root)
        if not neural_mat.is_file():
            raise FileNotFoundError(f"Stage-06 neural MAT is missing: {neural_mat}")
        sidecar = _load_json(sidecar_json)
        if sidecar.get("condition_name") != condition:
            raise RuntimeError(
                f"Condition sidecar mismatch: expected {condition}, "
                f"found {sidecar.get('condition_name')!r} in {sidecar_json}"
            )
        identity_values = _as_list_of_ints(
            sidecar.get("identity_values", []),
            f"{condition} identity_values",
        )
        time_axis = [float(value) for value in sidecar.get("time_axis_sec", [])]
        if not time_axis:
            raise RuntimeError(f"{sidecar_json} is missing time_axis_sec.")
        sampling_rate_hz = float(sidecar["sampling_rate_hz"])
        channel_labels = [str(value) for value in sidecar["ordered_channel_labels"]]
        if identity_values_ref is None:
            identity_values_ref = identity_values
            time_axis_ref = time_axis
            sampling_rate_ref = sampling_rate_hz
            channel_labels_ref = channel_labels
        elif (
            identity_values != identity_values_ref
            or time_axis != time_axis_ref
            or sampling_rate_hz != sampling_rate_ref
            or channel_labels != channel_labels_ref
        ):
            raise RuntimeError(
                "Stage-06 condition exports disagree on identity, time, "
                f"sampling, or channel axis for condition {condition}."
            )
        condition_payloads.append(
            {
                "name": condition,
                "code": int(sidecar["condition_code"]),
                "neural_mat": str(neural_mat),
                "sidecar_json": str(sidecar_json),
                "tensor_shape": [int(value) for value in sidecar["tensor_shape"]],
                "identity_values": identity_values,
                "counts_by_identity": sidecar["counts_by_identity"],
                "missing_identities": sidecar["missing_identities"],
            }
        )

    if identity_values_ref is None or time_axis_ref is None or sampling_rate_ref is None:
        raise RuntimeError("No Stage-06 condition outputs were validated.")

    for output_path in (
        paths["handoff_json"],
        paths["report_json"],
        paths["report_html"],
    ):
        guarded_path(output_path, raw_root)

    return {
        "version": 2,
        "generated_utc": utc_now(),
        "stage": "matlab_drsa",
        "status": "validated",
        "selected_analyses": list(selected_analyses),
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "subject_id": cfg["subject"]["id"],
        "project_root": str(paths["project_root"]),
        "dad_root": str(paths["dad_root"]),
        "data_root": str(raw_root),
        "derivatives_analysis_dir": str(paths["analysis_dir"]),
        "derivatives_reports_dir": str(paths["reports_dir"]),
        "drsa_dir": str(paths["drsa_dir"]),
        "drsa_functions_dir": str(paths["drsa_functions_dir"]),
        "source_manifest": str(manifest_path),
        "export_report_json": str(paths["export_report_json"]),
        "trajectory_models": trajectory_models,
        "condition_outputs": condition_payloads,
        "identity_values": identity_values_ref,
        "identity_count": len(identity_values_ref),
        "time_axis_sec": time_axis_ref,
        "sampling_rate_hz": sampling_rate_ref,
        "channel_labels": channel_labels_ref,
        "channel_count": len(channel_labels_ref or []),
        "occlusion_timing": _occlusion_timing_from_config(cfg),
    }


def _matlab_quote(value: str | Path) -> str:
    """Quote one value for a MATLAB single-quoted string literal."""
    return "'" + str(value).replace("'", "''") + "'"


def _matlab_cellstr(values: Iterable[str]) -> str:
    """Build a MATLAB cellstr literal from Python strings."""
    return "{" + ",".join(_matlab_quote(value) for value in values) + "}"


def _matlab_bin_from_config(args: argparse.Namespace, cfg: Mapping[str, Any]) -> Path:
    """Resolve and validate the absolute MATLAB executable path."""
    configured = (
        args.matlab_bin
        or os.environ.get("DAD_MATLAB_BIN")
        or cfg.get("matlab_drsa", {}).get("matlab_bin")
    )
    if not configured:
        raise RuntimeError(
            "No MATLAB executable configured. Set --matlab-bin, DAD_MATLAB_BIN, "
            "or matlab_drsa.matlab_bin."
        )
    matlab_bin = Path(str(configured)).expanduser()
    if not matlab_bin.is_absolute():
        raise RuntimeError(f"MATLAB executable must be an absolute path: {configured}")
    if not matlab_bin.is_file():
        raise FileNotFoundError(f"MATLAB executable not found: {matlab_bin}")
    return matlab_bin


def expected_outputs_for_analysis(
    cfg: Mapping[str, Any],
    analyses: Iterable[str],
    conditions: Iterable[str],
    sensor_groups: Iterable[str],
    drsa_type: str,
    model_source: str | Iterable[str],
    *,
    smoothing_plots: bool = False,
) -> list[Path]:
    """Return representative output files used for --resume decisions."""
    paths = _analysis_paths(cfg)
    drsa_dir = paths["drsa_dir"]
    model_sources = _normalize_model_sources(
        [model_source] if isinstance(model_source, str) else model_source
    )
    outputs: list[Path] = []
    for analysis in analyses:
        if analysis == "real":
            for condition in CONDITIONS:
                root = drsa_dir / "matlab_real_neural" / "07_dRSA" / condition / drsa_type.lower()
                outputs.extend(
                    [
                        root / "matrices" / f"dRSA_{condition}_{drsa_type.lower()}_matrices.png",
                        root / "results" / f"dRSA_{condition}_{drsa_type.lower()}_results.mat",
                    ]
                )
                if smoothing_plots:
                    outputs.append(
                        root
                        / "matrices"
                        / f"dRSA_{condition}_{drsa_type.lower()}_matrices_smoothed.png"
                    )
        elif analysis == "alt":
            for source in _model_sources_for_analysis(analysis, model_sources):
                for condition in _conditions_for_model_source(conditions, source):
                    for sensor_group in sensor_groups:
                        source_parts = (
                            [condition, source, drsa_type.lower()]
                            if condition == "occluded_deviant"
                            else [condition, drsa_type.lower()]
                        )
                        root = drsa_dir / "matlab_real_neural" / "A_dRSA_alt"
                        for part in source_parts:
                            root = root / part
                        tag = (
                            f"A_dRSA_alt_{condition}_{source}_"
                            f"{drsa_type.lower()}_{sensor_group.lower()}"
                            if condition == "occluded_deviant"
                            else f"A_dRSA_alt_{condition}_{drsa_type.lower()}_{sensor_group.lower()}"
                        )
                        outputs.extend(
                            [
                                root / "matrices" / f"{tag}.png",
                                root / "results" / f"{tag}.mat",
                            ]
                        )
                        if smoothing_plots:
                            outputs.append(root / "matrices" / f"{tag}_smoothed.png")
        elif analysis == "simulations":
            simulation_sources = _model_sources_for_analysis(analysis, model_sources)
            for source in simulation_sources:
                for condition in _conditions_for_model_source(conditions, source):
                    source_parts = (
                        [condition, source, drsa_type.lower()]
                        if condition == "occluded_deviant"
                        else [condition, drsa_type.lower()]
                    )
                    model_root = drsa_dir / "simulations" / "model_vs_model"
                    for part in source_parts:
                        model_root = model_root / part
                    model_tag = (
                        f"dRSA_simulations_{condition}_{source}_{drsa_type.lower()}"
                        if condition == "occluded_deviant"
                        else f"dRSA_simulations_{condition}_{drsa_type.lower()}"
                    )
                    outputs.extend(
                        [
                            model_root / "matrices" / f"{model_tag}_position_with_position.png",
                            model_root / "matrices" / f"{model_tag}_position_with_direction.png",
                            model_root / "results" / f"{model_tag}_model_vs_model.mat",
                        ]
                    )
                    if smoothing_plots:
                        outputs.extend(
                            [
                                model_root
                                / "matrices"
                                / f"{model_tag}_position_with_position_smoothed.png",
                                model_root
                                / "matrices"
                                / f"{model_tag}_position_with_direction_smoothed.png",
                            ]
                        )
            if simulation_sources:
                for condition in conditions:
                    neural_root = (
                        drsa_dir
                        / "simulations"
                        / "neural_vs_neural"
                        / condition
                        / "corr"
                    )
                    for sensor_group in sensor_groups:
                        tag = (
                            f"dRSA_simulations_{condition}_corr_neural_with_neural_"
                            f"{sensor_group.lower()}"
                        )
                        outputs.extend(
                            [
                                neural_root / "matrices" / f"{tag}.png",
                                neural_root / "results" / f"{tag}.mat",
                            ]
                        )
                        if smoothing_plots:
                            outputs.append(neural_root / "matrices" / f"{tag}_smoothed.png")
    return outputs


def _build_matlab_batches(
    cfg: Mapping[str, Any],
    *,
    handoff_path: Path,
    analyses: Iterable[str],
    drsa_type: str,
    conditions: Iterable[str],
    sensor_groups: Iterable[str],
    model_source: str | Iterable[str],
    smoothing_plots: bool = False,
    plot_smoothing_window: int = 9,
) -> list[tuple[str, str]]:
    """Build MATLAB batch source for selected analysis families."""
    paths = _analysis_paths(cfg)
    matlab_cfg = cfg.get("matlab_drsa", {})
    condition_list = [str(condition) for condition in conditions]
    model_sources = _normalize_model_sources(
        [model_source] if isinstance(model_source, str) else model_source
    )
    handoff_literal = _matlab_quote(handoff_path)
    drsa_type_literal = _matlab_quote(drsa_type)
    sensor_literal = _matlab_cellstr(sensor_groups)
    pcr_strategy_literal = _matlab_quote(
        matlab_cfg.get("pcr_strategy", "baseline_pcr_border")
    )
    smooth_neural_literal = (
        "true" if bool(matlab_cfg.get("smooth_neural", True)) else "false"
    )
    smooth_window_literal = str(int(matlab_cfg.get("smooth_window_samples", 15)))
    smoothing_plots_literal = "true" if smoothing_plots else "false"
    plot_smoothing_window_literal = str(int(plot_smoothing_window))
    shared_runtime = (
        f"drsa_handoff_json={handoff_literal};"
        f"dRSAtypeToRun={drsa_type_literal};"
        f"pcrStrategy={pcr_strategy_literal};"
        f"doSmoothNeural={smooth_neural_literal};"
        f"smoothWindowSamples={smooth_window_literal};"
        f"doSmoothingPlots={smoothing_plots_literal};"
        f"plotSmoothingWindow={plot_smoothing_window_literal};"
    )
    # One participant MATLAB process must respect its scheduler CPU allocation;
    # library environment variables alone do not cap MATLAB's compute threads.
    thread_limit = os.environ.get("DAD_MATLAB_MAX_THREADS")
    if thread_limit is not None:
        if not thread_limit.isdigit() or int(thread_limit) < 1:
            raise ValueError("DAD_MATLAB_MAX_THREADS must be a positive integer.")
        allocated = os.environ.get("SLURM_CPUS_PER_TASK")
        if allocated and int(thread_limit) > int(allocated):
            raise ValueError("MATLAB thread limit exceeds allocated SLURM CPUs.")
        shared_runtime = f"maxNumCompThreads({int(thread_limit)});" + shared_runtime
    batches: list[tuple[str, str]] = []

    if "real" in analyses:
        script = _matlab_quote(paths["real_script"])
        batches.append(
            (
                "real",
                (
                    shared_runtime +
                    f"run({script});"
                ),
            )
        )
    if "alt" in analyses:
        script = _matlab_quote(paths["alt_script"])
        for source in _model_sources_for_analysis("alt", model_sources):
            source_conditions = _conditions_for_model_source(condition_list, source)
            if not source_conditions:
                continue
            condition_literal = _matlab_cellstr(source_conditions)
            model_source_literal = _matlab_quote(source)
            batches.append(
                (
                    f"alt:{source}",
                    (
                        shared_runtime +
                        f"conditions_to_run={condition_literal};"
                        f"sensor_groups_to_run={sensor_literal};"
                        f"path_type={model_source_literal};"
                        f"run({script});"
                    ),
                )
            )
    if "simulations" in analyses:
        script = _matlab_quote(paths["simulations_script"])
        for source in _model_sources_for_analysis("simulations", model_sources):
            source_conditions = _conditions_for_model_source(condition_list, source)
            if not source_conditions:
                continue
            condition_literal = _matlab_cellstr(source_conditions)
            model_source_literal = _matlab_quote(source)
            batches.append(
                (
                    f"simulations:{source}",
                    (
                        shared_runtime +
                        f"conditions_to_run={condition_literal};"
                        f"sensor_groups_to_run={sensor_literal};"
                        f"model_path_source={model_source_literal};"
                        f"run({script});"
                    ),
                )
            )
    return batches


def _matlab_batch_command(
    matlab_bin: Path,
    batch_code: str,
    *,
    use_arch_x86_64: bool,
    platform_name: str | None = None,
) -> list[str]:
    """
    Build the noninteractive MATLAB command used for one dRSA batch.

    Parameters
    ----------
    matlab_bin
        Absolute MATLAB executable selected by configuration or the cluster
        wrapper.
    batch_code
        MATLAB source expression generated from the validated Stage-07
        handoff.
    use_arch_x86_64
        Prefix the command with ``arch -x86_64`` for local R2020a execution on
        Apple Silicon.
    platform_name
        Optional platform override for deterministic tests. Runtime calls use
        ``sys.platform``.

    Returns
    -------
    list of str
        Argument vector suitable for ``subprocess.run``.

    Notes
    -----
    The JVM stays enabled because R2020a disables ``figure`` under ``-nojvm``.
    Linux cluster jobs use software OpenGL because compute nodes have no
    display or physical graphics context; this exact mode was validated with
    ``figure`` and ``exportgraphics`` on a CIMeC SGE compute node.
    """
    runtime_platform = sys.platform if platform_name is None else platform_name
    command = [
        str(matlab_bin),
        "-nodesktop",
        "-nosplash",
        "-nodisplay",
    ]
    if runtime_platform.startswith("linux"):
        command.append("-softwareopengl")
    command.extend(["-batch", batch_code])
    if use_arch_x86_64:
        command = ["arch", "-x86_64", *command]
    return command


def run_matlab_batches(
    batches: Iterable[tuple[str, str]],
    *,
    matlab_bin: Path,
    use_arch_x86_64: bool,
    print_summary: bool,
) -> list[dict[str, Any]]:
    """
    Run MATLAB batches sequentially and return subprocess provenance records.

    Each generated analysis expression flows through the shared headless
    command builder. Processing stops after the first nonzero MATLAB return
    code so later batches cannot obscure the original failure.
    """
    records: list[dict[str, Any]] = []
    for analysis, batch_code in batches:
        command = _matlab_batch_command(
            matlab_bin,
            batch_code,
            use_arch_x86_64=use_arch_x86_64,
        )
        if print_summary:
            print(f"[matlab_drsa:{analysis}] {shlex.join(command)}")
        start = datetime.now(timezone.utc)
        completed = subprocess.run(command, check=False)
        finish = datetime.now(timezone.utc)
        records.append(
            {
                "analysis": analysis,
                "command": command,
                "command_display": shlex.join(command),
                "started_utc": start.isoformat().replace("+00:00", "Z"),
                "finished_utc": finish.isoformat().replace("+00:00", "Z"),
                "return_code": int(completed.returncode),
                "status": "completed" if completed.returncode == 0 else "failed",
            }
        )
        if completed.returncode != 0:
            break
    return records


def scan_drsa_outputs(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """Scan MATLAB output folders for reportable PNG/MAT products."""
    paths = _analysis_paths(cfg)
    drsa_dir = paths["drsa_dir"]
    pngs = sorted(path for path in drsa_dir.rglob("*.png") if path.is_file())
    mats = sorted(
        path
        for path in drsa_dir.rglob("*.mat")
        if path.is_file() and "desc-identityAveraged_neural" not in path.name
    )
    return {
        "drsa_dir": str(drsa_dir),
        "png_count": len(pngs),
        "mat_result_count": len(mats),
        "pngs": [str(path) for path in pngs],
        "mat_results": [str(path) for path in mats],
    }


def write_report(
    cfg: Mapping[str, Any],
    *,
    handoff: Mapping[str, Any],
    handoff_path: Path,
    matlab_records: list[dict[str, Any]],
    selected_analyses: Iterable[str],
    model_sources: Iterable[str],
    status: str,
    raw_root: Path,
    overwrite: bool,
    smoothing_plots: bool,
    plot_smoothing_window: int,
) -> dict[str, Any]:
    """Write JSON and HTML reports for the MATLAB dRSA layer."""
    paths = _analysis_paths(cfg)
    scan = scan_drsa_outputs(cfg)
    model_source_list = list(model_sources)
    report = {
        "stage": "matlab_drsa",
        "status": status,
        "generated_utc": utc_now(),
        "subject_bids_id": cfg["layout"]["subject_bids_id"],
        "session_bids_id": cfg["layout"]["session_bids_id"],
        "selected_analyses": list(selected_analyses),
        "model_sources": model_source_list,
        "smoothing_plots": smoothing_plots,
        "plot_smoothing_window": plot_smoothing_window,
        "handoff_json": str(handoff_path),
        "matlab_records": matlab_records,
        "handoff_summary": {
            "identity_count": handoff["identity_count"],
            "conditions": [item["name"] for item in handoff["condition_outputs"]],
            "sampling_rate_hz": handoff["sampling_rate_hz"],
            "channel_count": handoff["channel_count"],
            "trajectory_models": handoff["trajectory_models"],
        },
        "outputs": scan,
    }
    write_json(paths["report_json"], report, raw_root=raw_root, overwrite=overwrite)
    report_html = write_simple_html_report(
        paths["report_html"],
        title=f"{cfg['layout']['subject_bids_id']} | 07 MATLAB dRSA",
        summary_lines=[
            f"Stage status: {status}",
            f"Selected analyses: {', '.join(selected_analyses)}",
            f"Model sources: {', '.join(model_source_list)}",
            f"Smoothing plots: {smoothing_plots}",
            f"Plot smoothing window: {plot_smoothing_window}",
            f"PNG figures: {scan['png_count']}",
            f"MAT results: {scan['mat_result_count']}",
            f"Identity count: {handoff['identity_count']}",
            f"Output directory: {scan['drsa_dir']}",
            f"Machine report: {paths['report_json']}",
        ],
        image_paths=[Path(path) for path in scan["pngs"]],
        raw_root=raw_root,
    )
    embedded_html = write_embedded_html_copy(report_html, raw_root=raw_root)
    report["report_embedded_html"] = str(embedded_html)
    write_json(paths["report_json"], report, raw_root=raw_root, overwrite=True)
    return report


def _load_hdf5_mat_arrays(path: Path) -> dict[str, Any]:
    """Load numeric arrays from a MATLAB v7.3 file when h5py is available."""
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Cannot read MATLAB -v7.3 result MAT files because h5py is not "
            f"installed. Re-run Stage 07 after this update to write v7 result "
            f"files, or install h5py in the active Python environment. File: {path}"
        ) from exc

    payload: dict[str, Any] = {}
    with h5py.File(path, "r") as handle:
        for name in (*DRSA_MATRIX_VARIABLES, "time_axis_sec"):
            if name not in handle:
                continue
            array = handle[name][()]
            # MATLAB v7.3 stores arrays in HDF5 dimension order; reversing axes
            # restores the usual MATLAB/scipy [time x time x model] layout.
            if getattr(array, "ndim", 0) > 1:
                array = array.transpose(tuple(reversed(range(array.ndim))))
            payload[name] = array
    return payload


def _is_matlab_v73_file(path: Path) -> bool:
    """Detect MATLAB v7.3/HDF5 MAT files from the cheap header bytes."""
    try:
        with path.open("rb") as file_obj:
            header = file_obj.read(128)
    except OSError:
        return False
    return b"MATLAB 7.3" in header or header.startswith(b"\x89HDF")


def _load_mat_arrays(path: Path) -> dict[str, Any]:
    """Load candidate numeric matrices from one MATLAB result file."""
    if _is_matlab_v73_file(path):
        return _load_hdf5_mat_arrays(path)

    try:
        from scipy.io import loadmat
    except ImportError as exc:
        raise RuntimeError(
            "Group averaging requires scipy to read Stage-07 MAT result files."
        ) from exc

    try:
        return loadmat(
            path,
            variable_names=[*DRSA_MATRIX_VARIABLES, "time_axis_sec"],
            squeeze_me=True,
            struct_as_record=False,
        )
    except (NotImplementedError, ValueError) as exc:
        message = str(exc).lower()
        if "7.3" in message or "hdf" in message or "unknown mat file type" in message:
            return _load_hdf5_mat_arrays(path)
        raise


def _as_numeric_matrix(value: Any, *, label: str, path: Path):
    """Normalize one loaded MAT variable to a finite 2-D/3-D float array."""
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Group averaging requires numpy.") from exc

    array = np.asarray(value, dtype=float).squeeze()
    if array.ndim not in {2, 3}:
        raise RuntimeError(
            f"{label} in {path} must be a 2-D or 3-D numeric matrix, "
            f"found shape {array.shape}."
        )
    if array.ndim == 3 and array.shape[0] != array.shape[1]:
        if array.shape[1] == array.shape[2]:
            array = np.moveaxis(array, 0, -1)
        else:
            raise RuntimeError(
                f"{label} in {path} has unsupported 3-D shape {array.shape}."
            )
    return array


def _read_drsa_result_matrices(path: Path) -> tuple[dict[str, Any], Any | None]:
    """Read all group-averageable dRSA matrices from one result MAT file."""
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Group averaging requires numpy.") from exc

    payload = _load_mat_arrays(path)
    matrices: dict[str, Any] = {}
    for variable in DRSA_MATRIX_VARIABLES:
        if variable in payload:
            matrices[variable] = _as_numeric_matrix(
                payload[variable],
                label=variable,
                path=path,
            )
    time_axis = None
    if "time_axis_sec" in payload:
        candidate = np.asarray(payload["time_axis_sec"], dtype=float).squeeze()
        if candidate.ndim == 1 and candidate.size > 1:
            time_axis = candidate
    return matrices, time_axis


def _smooth_matrix_for_plot_python(matrix_2d: Any, window_samples: int):
    """
    Apply the shared display smoother used by Stage 08 line summaries.

    Keeping this thin wrapper lets the Stage-07 group-average plotting path
    preserve the same undefined matrix support as the exploratory summaries.
    """
    return shared_smooth_matrix_for_plot(matrix_2d, window_samples)


def _smooth_array_for_plot(matrix: Any, window_samples: int):
    """Smooth every 2-D matrix slice in a 2-D or 3-D dRSA result array."""
    import numpy as np

    array = np.asarray(matrix, dtype=float)
    if array.ndim == 2:
        return _smooth_matrix_for_plot_python(array, window_samples)
    slices = [
        _smooth_matrix_for_plot_python(array[:, :, i_model], window_samples)
        for i_model in range(array.shape[2])
    ]
    return np.stack(slices, axis=2)


def _iter_matrix_slices(matrix: Any) -> Iterable[tuple[str, Any]]:
    """Yield named 2-D slices from one averaged dRSA matrix array."""
    import numpy as np

    array = np.asarray(matrix, dtype=float)
    if array.ndim == 2:
        yield "", array
        return
    for i_model in range(array.shape[2]):
        yield f"_model{i_model + 1:02d}", array[:, :, i_model]


def _time_axis_for_matrix(
    matrix_2d: Any,
    *,
    time_axis: Any | None,
    sampling_rate_hz: float,
) -> Any:
    """Return a plotting time axis that matches one dRSA matrix."""
    import numpy as np

    n_time = int(matrix_2d.shape[0])
    if time_axis is not None and len(time_axis) == n_time:
        return np.asarray(time_axis, dtype=float)
    return np.arange(n_time, dtype=float) / float(sampling_rate_hz)


def _compact_plot_title(parts: Iterable[Any], *, width: int = 74) -> str:
    """
    Return a wrapped, compact Matplotlib title from plot-identifying fields.

    Data flow:
    group plot metadata -> short PNG title. Full model and source details stay
    in the HTML figure caption so long labels cannot crowd the axes.
    """
    text = " | ".join(str(part).strip() for part in parts if str(part).strip())
    if not text:
        return ""
    return "\n".join(
        textwrap.wrap(
            text,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )


def _plot_group_matrix(
    matrix_2d: Any,
    *,
    out_file: Path,
    title: str,
    model_axis_label: str,
    condition: str,
    occlusion_timing: Mapping[str, Any] | None,
    time_axis: Any | None,
    sampling_rate_hz: float,
) -> None:
    """Write one group-average dRSA matrix PNG using matplotlib Agg."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    axis_sec = _time_axis_for_matrix(
        matrix_2d,
        time_axis=time_axis,
        sampling_rate_hz=sampling_rate_hz,
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    image = ax.imshow(
        matrix_2d,
        origin="lower",
        aspect="equal",
        extent=[axis_sec[0], axis_sec[-1], axis_sec[0], axis_sec[-1]],
    )
    ax.plot(axis_sec, axis_sec, color="white", linewidth=1.0)
    for line_time in _occlusion_line_times(condition, occlusion_timing):
        if float(axis_sec[0]) <= line_time <= float(axis_sec[-1]):
            ax.axhline(
                line_time,
                color="white",
                linestyle="--",
                linewidth=1.3,
            )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("neural time (s)")
    ax.set_ylabel(model_axis_label)
    fig.colorbar(image, ax=ax)
    fig.savefig(out_file, dpi=300, facecolor="white")
    plt.close(fig)


def _group_output_parent(relative_result_path: Path) -> Path:
    """Map a subject result MAT relative path to a group matrices folder."""
    parts = list(relative_result_path.parts[:-1])
    if "results" in parts:
        parts[parts.index("results")] = "matrices"
    return Path(*parts) if parts else Path()


def _sampling_rate_for_group_plots(cfg: Mapping[str, Any]) -> float:
    """Resolve the group plot sampling rate when result MAT lacks a time axis."""
    return float(
        cfg.get("export", {}).get(
            "sampling_rate_hz",
            cfg.get("matlab_drsa", {}).get("plot_sample_rate_hz", 120.0),
        )
    )


def _figure_variant(path: Path) -> str:
    """Classify a figure path as raw or smoothed from its filename."""
    return "smoothed" if path.stem.endswith("_smoothed") else "raw"


def _condition_from_path(path: Path) -> str:
    """Infer the dRSA condition from a relative result or figure path."""
    parts = [part.lower() for part in path.parts]
    name = path.name.lower()
    for condition in CONDITIONS:
        if condition.lower() in parts or condition.lower() in name:
            return condition
    return "unknown"


def _analysis_family_from_path(path: Path) -> str:
    """Infer a readable analysis family label from a relative dRSA path."""
    text = str(path).lower()
    if "a_drsa_alt" in text:
        return "A_dRSA_alt"
    if "07_drsa" in text:
        return "A_dRSA"
    if "model_vs_model" in text:
        return "simulations model-vs-model"
    if "neural_vs_neural" in text:
        return "simulations neural-vs-neural"
    return "dRSA"


def discover_participant_figures(
    cfg: Mapping[str, Any],
    subjects: Iterable[str],
) -> list[dict[str, str]]:
    """Discover subject-level dRSA PNGs for inclusion in group reports."""
    paths = _analysis_paths(cfg)
    session = cfg["layout"]["session_bids_id"]
    analysis_root = paths["analysis_dir"].parent.parent
    figures: list[dict[str, str]] = []

    for subject in subjects:
        subject_drsa_dir = analysis_root / subject / session / "drsa"
        if not subject_drsa_dir.is_dir():
            continue
        for figure_path in sorted(subject_drsa_dir.rglob("*.png")):
            relative = figure_path.relative_to(subject_drsa_dir)
            condition = _condition_from_path(relative)
            variant = _figure_variant(figure_path)
            figures.append(
                {
                    "subject": subject,
                    "condition": condition,
                    "variant": variant,
                    "analysis_family": _analysis_family_from_path(relative),
                    "figure": str(figure_path),
                    "relative_drsa_path": str(relative),
                    "caption": (
                        f"{subject} | {condition} | "
                        f"{_analysis_family_from_path(relative)} | {relative.name}"
                    ),
                }
            )

    order = {condition: index for index, condition in enumerate(CONDITIONS)}
    figures.sort(
        key=lambda item: (
            item["subject"],
            order.get(item["condition"], len(order)),
            item["variant"],
            item["analysis_family"],
            item["relative_drsa_path"],
        )
    )
    return figures


def _html_id(value: str) -> str:
    """Create a simple stable HTML anchor from a label."""
    safe = "".join(
        char.lower() if char.isalnum() else "-"
        for char in str(value)
    ).strip("-")
    while "--" in safe:
        safe = safe.replace("--", "-")
    return safe or "section"


def _image_href(report_path: Path, image_path: str | Path) -> str:
    """Return an HTML image href relative to the report folder when possible."""
    resolved_image = Path(image_path).expanduser().resolve()
    try:
        return str(resolved_image.relative_to(report_path.parent))
    except ValueError:
        return str(resolved_image)


def _figure_block(report_path: Path, item: Mapping[str, Any]) -> str:
    """Build one HTML figure block with a portable image reference."""
    href = _image_href(report_path, item["figure"])
    caption = str(item.get("caption") or item.get("figure"))
    return (
        "<figure>"
        f"<img src=\"{html.escape(href)}\"/>"
        f"<figcaption>{html.escape(caption)}</figcaption>"
        "</figure>"
    )


def _group_output_caption(item: Mapping[str, Any], variant_label: str) -> str:
    """Build a readable caption for one group-average figure."""
    subjects = ", ".join(str(value) for value in item.get("subjects", []))
    condition = str(item.get("condition", "unknown"))
    model_label = str(item.get("stimulus_model_label", "stimulus model"))
    return (
        f"group {variant_label} | {condition} | {model_label} | "
        f"{item.get('relative_result', 'unknown')} | {item.get('variable', 'matrix')} | "
        f"n={item.get('n_subjects', '?')} | "
        f"subjects: {subjects}"
    )


def _condition_grouped_figure_blocks(
    report_path: Path,
    figures: Iterable[Mapping[str, Any]],
) -> str:
    """Render figure blocks grouped by dRSA condition."""
    figure_list = [dict(item) for item in figures]
    blocks: list[str] = []
    for condition in CONDITIONS:
        condition_figures = [
            item for item in figure_list
            if item.get("condition") == condition
        ]
        if not condition_figures:
            continue
        blocks.append(
            f"<section><h3>{html.escape(condition)}</h3>"
            + "\n".join(_figure_block(report_path, item) for item in condition_figures)
            + "</section>"
        )
    unknown_figures = [
        item for item in figure_list
        if item.get("condition") not in CONDITIONS
    ]
    if unknown_figures:
        blocks.append(
            "<section><h3>unknown condition</h3>"
            + "\n".join(_figure_block(report_path, item) for item in unknown_figures)
            + "</section>"
        )
    return "\n".join(blocks)


def write_group_drsa_html_report(
    path: Path,
    *,
    title: str,
    variant_label: str,
    summary_lines: Iterable[str],
    subjects: Iterable[str],
    group_outputs: Iterable[Mapping[str, Any]],
    participant_figures: Iterable[Mapping[str, Any]],
    skipped_inputs: Iterable[Mapping[str, str]],
    raw_root: Path,
) -> Path:
    """
    Write a navigable group dRSA report with group and participant sections.

    The report is intentionally static HTML rather than an MNE Report so it can
    be opened from the cluster derivative folder without Python dependencies.
    """
    report_path = guarded_path(path, raw_root)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    subject_list = list(subjects)
    group_items = [dict(item) for item in group_outputs]
    participant_items = [dict(item) for item in participant_figures]
    skipped = list(skipped_inputs)

    for item in group_items:
        item["caption"] = _group_output_caption(item, variant_label)

    participant_by_subject = {
        subject: [
            item for item in participant_items
            if item.get("subject") == subject
        ]
        for subject in subject_list
    }

    toc_subject_links = "\n".join(
        f"<li><a href=\"#{_html_id('subject-' + subject)}\">{html.escape(subject)}</a></li>"
        for subject in subject_list
    )
    group_blocks = _condition_grouped_figure_blocks(report_path, group_items)
    if not group_blocks:
        group_blocks = "<p>No group-average figures were available for this report variant.</p>"

    subject_sections: list[str] = []
    for subject in subject_list:
        subject_id = _html_id("subject-" + subject)
        condition_blocks: list[str] = []
        subject_figures = participant_by_subject[subject]
        for condition in CONDITIONS:
            condition_figures = [
                item for item in subject_figures
                if item.get("condition") == condition
            ]
            blocks = "\n".join(
                _figure_block(report_path, item) for item in condition_figures
            ) or "<p>No figures for this condition.</p>"
            condition_blocks.append(
                f"<section><h3>{html.escape(condition)}</h3>{blocks}</section>"
            )
        unknown_figures = [
            item for item in subject_figures
            if item.get("condition") not in CONDITIONS
        ]
        if unknown_figures:
            condition_blocks.append(
                "<section><h3>unknown condition</h3>"
                + "\n".join(_figure_block(report_path, item) for item in unknown_figures)
                + "</section>"
            )
        subject_sections.append(
            f"<section id=\"{subject_id}\"><h2>{html.escape(subject)}</h2>"
            + "\n".join(condition_blocks)
            + "</section>"
        )

    skipped_block = ""
    if skipped:
        skipped_rows = "\n".join(
            "<li>"
            + html.escape(
                " | ".join(f"{key}: {value}" for key, value in sorted(item.items()))
            )
            + "</li>"
            for item in skipped
        )
        skipped_block = f"<details><summary>Skipped inputs</summary><ul>{skipped_rows}</ul></details>"

    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\">",
            f"<title>{html.escape(title)}</title>",
            "<style>",
            "body{font-family:system-ui,sans-serif;max-width:1600px;margin:2rem auto;padding:0 1rem;}",
            "nav,pre{background:#f5f5f5;padding:1rem;border-radius:6px;}",
            "figure{margin:1.5rem 0;padding-bottom:1rem;border-bottom:1px solid #ddd;}",
            "img{max-width:100%;height:auto;border:1px solid #ddd;}",
            "figcaption{font-size:.9rem;color:#333;margin-top:.4rem;}",
            "h2{margin-top:3rem;border-top:3px solid #222;padding-top:1rem;}",
            "h3{margin-top:2rem;color:#333;}",
            "</style>",
            "</head><body>",
            f"<h1>{html.escape(title)}</h1>",
            f"<p><strong>Report variant:</strong> {html.escape(variant_label)}</p>",
            f"<pre>{html.escape(chr(10).join(summary_lines))}</pre>",
            "<nav><strong>Subjects</strong><ul>",
            toc_subject_links,
            "</ul></nav>",
            skipped_block,
            "<section id=\"group-average\"><h2>Group average figures</h2>",
            group_blocks,
            "</section>",
            "<section id=\"participants\"><h2>Participant figures</h2>",
            "\n".join(subject_sections),
            "</section>",
            "</body></html>",
        ]
    )
    report_path.write_text(body, encoding="utf-8")
    write_embedded_html_copy(report_path, raw_root=raw_root)
    return report_path


def build_group_average_report(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """
    Average Stage-07 result matrices numerically across subjects.

    Data flow:
    subject result MAT files -> variable/path grouping -> raw mean matrices;
    when smoothing plots are enabled, each subject matrix is smoothed first and
    those smoothed matrices are averaged into a separate display-only figure.
    """
    import numpy as np

    raw_root = Path(cfg["resolved_paths"]["data_root"])
    paths = _analysis_paths(cfg)
    session = cfg["layout"]["session_bids_id"]
    subjects = _resolve_group_subjects(cfg, args)
    smoothing_plots, smoothing_window = _resolve_plot_smoothing(cfg, args)
    analysis_root = paths["analysis_dir"].parent.parent
    sampling_rate_hz = _sampling_rate_for_group_plots(cfg)
    occlusion_timing = _occlusion_timing_from_config(cfg)

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    skipped_inputs: list[dict[str, str]] = []

    for subject in subjects:
        subject_drsa_dir = analysis_root / subject / session / "drsa"
        if not subject_drsa_dir.is_dir():
            skipped_inputs.append(
                {"subject": subject, "reason": f"missing dRSA folder: {subject_drsa_dir}"}
            )
            continue
        for mat_path in sorted(subject_drsa_dir.rglob("*.mat")):
            if "desc-identityAveraged_neural" in mat_path.name:
                continue
            try:
                matrices, time_axis = _read_drsa_result_matrices(mat_path)
            except RuntimeError as exc:
                skipped_inputs.append(
                    {
                        "subject": subject,
                        "path": str(mat_path),
                        "reason": str(exc),
                    }
                )
                continue
            relative_path = mat_path.relative_to(subject_drsa_dir)
            for variable, matrix in matrices.items():
                grouped.setdefault((str(relative_path), variable), []).append(
                    {
                        "subject": subject,
                        "path": str(mat_path),
                        "matrix": matrix,
                        "time_axis": time_axis,
                    }
                )

    outputs: list[dict[str, Any]] = []
    for (relative_text, variable), items in sorted(grouped.items()):
        shapes = {tuple(item["matrix"].shape) for item in items}
        if len(shapes) != 1:
            skipped_inputs.append(
                {
                    "path": relative_text,
                    "reason": f"shape mismatch across subjects: {sorted(shapes)}",
                }
            )
            continue

        stack = np.stack([item["matrix"] for item in items], axis=0)
        averaged_raw = np.nanmean(stack, axis=0)
        first_time_axis = next(
            (item["time_axis"] for item in items if item["time_axis"] is not None),
            None,
        )
        relative_path = Path(relative_text)
        condition = _condition_from_path(relative_path)
        analysis_family = analysis_family_from_path(relative_path)
        path_model_source = model_source_from_path(relative_path)
        output_parent = _group_output_parent(relative_path)
        figure_dir = paths["group_drsa_dir"] / "group_average" / output_parent
        guarded_path(figure_dir, raw_root).mkdir(parents=True, exist_ok=True)

        base_stem = f"{relative_path.stem}_{variable}_groupAverage"
        for suffix, matrix_slice in _iter_matrix_slices(averaged_raw):
            slice_label = suffix.lstrip("_") or "matrix"
            model_label = stimulus_model_label(
                condition,
                path_model_source,
                variable,
                slice_label=slice_label,
                analysis_family=analysis_family,
            )
            model_axis_label = model_time_axis_label(
                condition,
                path_model_source,
                variable,
                slice_label=slice_label,
                analysis_family=analysis_family,
            )
            figure_path = figure_dir / f"{base_stem}{suffix}.png"
            guarded_path(figure_path, raw_root)
            _plot_group_matrix(
                matrix_slice,
                out_file=figure_path,
                title=_compact_plot_title(
                    [
                        "Group average raw",
                        condition,
                        analysis_family,
                        variable + suffix,
                        f"n={len(items)}",
                    ]
                ),
                model_axis_label=model_axis_label,
                condition=condition,
                occlusion_timing=occlusion_timing,
                time_axis=first_time_axis,
                sampling_rate_hz=sampling_rate_hz,
            )
            outputs.append(
                {
                    "kind": "raw_average",
                    "condition": condition,
                    "analysis_family": analysis_family,
                    "model_source": path_model_source,
                    "stimulus_model_label": model_label,
                    "model_time_axis_label": model_axis_label,
                    "variable": variable,
                    "slice_label": slice_label,
                    "relative_result": relative_text,
                    "subjects": [item["subject"] for item in items],
                    "n_subjects": len(items),
                    "figure": str(figure_path),
                }
            )

        if smoothing_plots:
            smoothed_stack = np.stack(
                [
                    _smooth_array_for_plot(item["matrix"], smoothing_window)
                    for item in items
                ],
                axis=0,
            )
            averaged_smoothed = np.nanmean(smoothed_stack, axis=0)
            for suffix, matrix_slice in _iter_matrix_slices(averaged_smoothed):
                slice_label = suffix.lstrip("_") or "matrix"
                model_label = stimulus_model_label(
                    condition,
                    path_model_source,
                    variable,
                    slice_label=slice_label,
                    analysis_family=analysis_family,
                )
                model_axis_label = model_time_axis_label(
                    condition,
                    path_model_source,
                    variable,
                    slice_label=slice_label,
                    analysis_family=analysis_family,
                )
                figure_path = figure_dir / f"{base_stem}{suffix}_smoothed.png"
                guarded_path(figure_path, raw_root)
                _plot_group_matrix(
                    matrix_slice,
                    out_file=figure_path,
                    title=_compact_plot_title(
                        [
                            "Group average smoothed",
                            condition,
                            analysis_family,
                            variable + suffix,
                            f"n={len(items)}",
                        ]
                    ),
                    model_axis_label=model_axis_label,
                    condition=condition,
                    occlusion_timing=occlusion_timing,
                    time_axis=first_time_axis,
                    sampling_rate_hz=sampling_rate_hz,
                )
                outputs.append(
                    {
                        "kind": "average_of_smoothed",
                        "condition": condition,
                        "analysis_family": analysis_family,
                        "model_source": path_model_source,
                        "stimulus_model_label": model_label,
                        "model_time_axis_label": model_axis_label,
                        "variable": variable,
                        "slice_label": slice_label,
                        "relative_result": relative_text,
                        "subjects": [item["subject"] for item in items],
                        "n_subjects": len(items),
                        "smoothing_window_samples": smoothing_window,
                        "figure": str(figure_path),
                    }
                )

    if not outputs:
        reasons = "; ".join(item.get("reason", "unknown") for item in skipped_inputs[:5])
        raise RuntimeError(
            "No group-average dRSA figures were produced. "
            f"First skipped-input reasons: {reasons or 'none'}"
        )

    participant_figures = discover_participant_figures(cfg, subjects)
    raw_group_outputs = [
        item for item in outputs if item.get("kind") == "raw_average"
    ]
    smoothed_group_outputs = [
        item for item in outputs if item.get("kind") == "average_of_smoothed"
    ]
    raw_participant_figures = [
        item for item in participant_figures if item.get("variant") == "raw"
    ]
    smoothed_participant_figures = [
        item for item in participant_figures if item.get("variant") == "smoothed"
    ]

    common_summary = [
        f"Subjects requested: {', '.join(subjects)}",
        f"All group-average figure count: {len(outputs)}",
        f"Participant figure count: {len(participant_figures)}",
        f"Smoothing plots requested: {smoothing_plots}",
        f"Plot smoothing window: {smoothing_window}",
        f"Skipped inputs: {len(skipped_inputs)}",
        f"Group dRSA directory: {paths['group_drsa_dir']}",
    ]
    raw_report = {
        "stage": "matlab_drsa_group_average",
        "status": "completed",
        "generated_utc": utc_now(),
        "variant": "raw",
        "session_bids_id": session,
        "subjects_requested": subjects,
        "smoothing_plots": smoothing_plots,
        "plot_smoothing_window": smoothing_window,
        "group_drsa_dir": str(paths["group_drsa_dir"]),
        "outputs": raw_group_outputs,
        "participant_figures": raw_participant_figures,
        "skipped_inputs": skipped_inputs,
        "notes": [
            "Raw group figures average the unsmoothed subject dRSA matrices.",
        ],
    }
    write_json(paths["group_report_json"], raw_report, raw_root=raw_root, overwrite=True)
    write_group_drsa_html_report(
        paths["group_report_html"],
        title=f"Group | {session} | 07 MATLAB dRSA averages | raw",
        variant_label="raw / non-smoothed",
        summary_lines=[
            "Stage status: completed",
            *common_summary,
            f"Raw group-average figures: {len(raw_group_outputs)}",
            f"Raw participant figures: {len(raw_participant_figures)}",
            f"Machine report: {paths['group_report_json']}",
        ],
        subjects=subjects,
        group_outputs=raw_group_outputs,
        participant_figures=raw_participant_figures,
        skipped_inputs=skipped_inputs,
        raw_root=raw_root,
    )

    report_files = {
        "raw_json": str(paths["group_report_json"]),
        "raw_html": str(paths["group_report_html"]),
        "raw_embedded_html": str(
            paths["group_report_html"].with_name(
                f"{paths['group_report_html'].stem}_embedded.html"
            )
        ),
    }
    smoothed_report = None
    if smoothing_plots:
        smoothed_report = {
            "stage": "matlab_drsa_group_average",
            "status": "completed",
            "generated_utc": utc_now(),
            "variant": "smoothed",
            "session_bids_id": session,
            "subjects_requested": subjects,
            "smoothing_plots": smoothing_plots,
            "plot_smoothing_window": smoothing_window,
            "group_drsa_dir": str(paths["group_drsa_dir"]),
            "outputs": smoothed_group_outputs,
            "participant_figures": smoothed_participant_figures,
            "skipped_inputs": skipped_inputs,
            "notes": [
                (
                    "Smoothed group figures smooth each subject matrix first "
                    "and then average those display matrices."
                ),
                (
                    "Participant sections include only *_smoothed.png figures "
                    "for this report variant."
                ),
            ],
        }
        write_json(
            paths["group_report_smoothed_json"],
            smoothed_report,
            raw_root=raw_root,
            overwrite=True,
        )
        write_group_drsa_html_report(
            paths["group_report_smoothed_html"],
            title=f"Group | {session} | 07 MATLAB dRSA averages | smoothed",
            variant_label="smoothed display",
            summary_lines=[
                "Stage status: completed",
                *common_summary,
                f"Smoothed group-average figures: {len(smoothed_group_outputs)}",
                f"Smoothed participant figures: {len(smoothed_participant_figures)}",
                f"Machine report: {paths['group_report_smoothed_json']}",
            ],
            subjects=subjects,
            group_outputs=smoothed_group_outputs,
            participant_figures=smoothed_participant_figures,
            skipped_inputs=skipped_inputs,
            raw_root=raw_root,
        )
        report_files.update(
            {
                "smoothed_json": str(paths["group_report_smoothed_json"]),
                "smoothed_html": str(paths["group_report_smoothed_html"]),
                "smoothed_embedded_html": str(
                    paths["group_report_smoothed_html"].with_name(
                        f"{paths['group_report_smoothed_html'].stem}_embedded.html"
                    )
                ),
            }
        )

    return {
        "stage": "matlab_drsa_group_average",
        "status": "completed",
        "generated_utc": utc_now(),
        "session_bids_id": session,
        "subjects_requested": subjects,
        "smoothing_plots": smoothing_plots,
        "plot_smoothing_window": smoothing_window,
        "group_drsa_dir": str(paths["group_drsa_dir"]),
        "outputs": outputs,
        "participant_figures": participant_figures,
        "raw_report": raw_report,
        "smoothed_report": smoothed_report,
        "report_files": report_files,
        "skipped_inputs": skipped_inputs,
    }


def run_stage(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Validate inputs, optionally launch MATLAB, and record outputs."""
    raw_root = Path(cfg["resolved_paths"]["data_root"])
    paths = _analysis_paths(cfg)
    matlab_cfg = cfg.get("matlab_drsa", {})
    if args.group_average_only:
        return build_group_average_report(cfg, args)
    if not args.source_manifest:
        raise RuntimeError("--source-manifest is required unless --group-average-only is used.")
    env_analyses = os.environ.get("DAD_DRSA_ANALYSES")
    env_conditions = os.environ.get("DAD_DRSA_CONDITIONS")
    env_sensor_groups = os.environ.get("DAD_DRSA_SENSOR_GROUPS")
    smoothing_plots, plot_smoothing_window = _resolve_plot_smoothing(cfg, args)
    analyses = list(
        args.analyses
        or (env_analyses.split(",") if env_analyses else None)
        or matlab_cfg.get("analyses", ["real"])
    )
    drsa_type = str(args.drsa_type or matlab_cfg.get("drsa_type", "corr"))
    conditions = list(
        args.conditions
        or (env_conditions.split(",") if env_conditions else None)
        or matlab_cfg.get("conditions_to_run", CONDITIONS)
    )
    sensor_groups = list(
        args.sensor_groups
        or (env_sensor_groups.split(",") if env_sensor_groups else None)
        or matlab_cfg.get("sensor_groups_to_run", ["wholeBrain"])
    )
    model_sources = _model_sources_from_config_or_args(args, matlab_cfg)
    if drsa_type.lower() == "pcr":
        drsa_type = "PCR"
    else:
        drsa_type = "corr"
    conditions = [str(condition) for condition in conditions]
    sensor_groups = [str(group) for group in sensor_groups]

    for analysis in analyses:
        if analysis not in SUPPORTED_ANALYSES:
            raise ValueError(f"Unsupported MATLAB dRSA analysis: {analysis}")
    for condition in conditions:
        if condition not in CONDITIONS:
            raise ValueError(f"Unsupported dRSA condition: {condition}")

    for path_key in (
        "real_script",
        "alt_script",
        "simulations_script",
        "drsa_functions_dir",
    ):
        value = paths[path_key]
        if not value.exists():
            raise FileNotFoundError(f"Required MATLAB analysis path is missing: {value}")

    handoff = build_handoff(
        cfg,
        source_manifest=args.source_manifest,
        selected_analyses=analyses,
    )
    if args.print_summary:
        print(
            json.dumps(
                {
                    "stage": "matlab_drsa",
                    "status": "validated",
                    "subject_bids_id": handoff["subject_bids_id"],
                    "smoothing_plots": smoothing_plots,
                    "plot_smoothing_window": plot_smoothing_window,
                    "model_sources": model_sources,
                    "identity_count": handoff["identity_count"],
                    "conditions": [
                        item["name"] for item in handoff["condition_outputs"]
                    ],
                    "sampling_rate_hz": handoff["sampling_rate_hz"],
                    "channel_count": handoff["channel_count"],
                    "drsa_dir": handoff["drsa_dir"],
                    "trajectory_models": handoff["trajectory_models"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    if args.validate_only:
        return {
            "stage": "matlab_drsa",
            "status": "validated",
            "handoff": handoff,
        }

    paths["reports_dir"].mkdir(parents=True, exist_ok=True)
    # Bind resume to the neural variant, requested settings and MATLAB sources.
    # Output MATs use file stamps; large raw recordings are never hashed here.
    analysis_contract = None
    previous = None
    if cfg.get("epoching", {}).get("processing_run"):
        export_evidence = verify_export_artifacts(cfg)
        analysis_contract = digest({"export": export_evidence["export_contract_sha256"],
            "settings": [analyses, conditions, sensor_groups, drsa_type, model_sources,
                smoothing_plots, plot_smoothing_window, cfg.get("matlab_drsa", {})],
            "scripts": {key: file_hash(paths[key]) for key in ("real_script", "alt_script", "simulations_script")}})
        if paths["report_json"].is_file():
            previous = _load_json(paths["report_json"])
            if previous.get("analysis_contract_sha256") != analysis_contract:
                raise RuntimeError("MATLAB inputs/settings changed; use a new processing-run label.")
            for binding in previous.get("artifact_bindings", []):
                if file_stamp(Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / binding["relative_path"]) != binding["stamp"]:
                    raise RuntimeError("Stored MATLAB output changed.")
    write_json(
        paths["handoff_json"],
        handoff,
        raw_root=raw_root,
        overwrite=True,
    )

    expected_outputs = expected_outputs_for_analysis(
        cfg,
        analyses,
        conditions,
        sensor_groups,
        drsa_type,
        model_sources,
        smoothing_plots=smoothing_plots,
    )
    if analysis_contract is not None and previous is None and any(path.exists() for path in expected_outputs):
        raise RuntimeError("Unbound MATLAB outputs already exist; preserve them and use a new processing-run label.")
    if (args.resume and expected_outputs and all(path.is_file() for path in expected_outputs)
            and (analysis_contract is None or previous is not None and previous.get("status") == "completed")):
        matlab_records: list[dict[str, Any]] = [
            {
                "analysis": ",".join(analyses),
                "status": "skipped",
                "reason": "resume_outputs_exist",
                "expected_outputs": [str(path) for path in expected_outputs],
            }
        ]
        status = "completed"
    else:
        matlab_bin = _matlab_bin_from_config(args, cfg)
        batches = _build_matlab_batches(
            cfg,
            handoff_path=paths["handoff_json"],
            analyses=analyses,
            drsa_type=drsa_type,
            conditions=conditions,
            sensor_groups=sensor_groups,
            model_source=model_sources,
            smoothing_plots=smoothing_plots,
            plot_smoothing_window=plot_smoothing_window,
        )
        if not batches:
            raise RuntimeError(
                "No MATLAB dRSA batches were generated. Predicted and "
                "PE model paths are only valid for occluded_deviant, and "
                "PE variants are available for alt outputs only; include "
                "observed or add occluded_deviant to --conditions."
            )
        matlab_records = run_matlab_batches(
            batches,
            matlab_bin=matlab_bin,
            use_arch_x86_64=bool(
                args.use_arch_x86_64 or matlab_cfg.get("use_arch_x86_64", False)
            ),
            print_summary=args.print_summary,
        )
        status = (
            "completed"
            if matlab_records and all(item["status"] == "completed" for item in matlab_records)
            else "failed"
        )

    report = write_report(
        cfg,
        handoff=handoff,
        handoff_path=paths["handoff_json"],
        matlab_records=matlab_records,
        selected_analyses=analyses,
        model_sources=model_sources,
        status=status,
        raw_root=raw_root,
        overwrite=True,
        smoothing_plots=smoothing_plots,
        plot_smoothing_window=plot_smoothing_window,
    )
    if analysis_contract is not None:
        from review_final_ica import settle_review_outputs
        present = [path for path in expected_outputs if path.is_file()]
        settle_review_outputs(present)
        analysis_root = Path(cfg["resolved_paths"]["derivatives_analysis_dir"])
        report["analysis_contract_sha256"] = analysis_contract
        report["artifact_bindings"] = [{"relative_path": str(path.relative_to(analysis_root)),
            "stamp": file_stamp(path)} for path in present]
        write_json(paths["report_json"], report, raw_root=raw_root, overwrite=True)
    if status == "completed":
        records = [
            {
                "input_path": str(Path(item["neural_mat"]).expanduser().resolve()),
                "output_path": str(paths["report_json"].resolve()),
                "status": "completed",
                "metadata": {
                    "condition": item["name"],
                    "identity_count": len(item["identity_values"]),
                },
            }
            for item in handoff["condition_outputs"]
        ]
        update_source_manifest_stage(
            manifest_path=processing_manifest(cfg, args.source_manifest),
            raw_root=raw_root,
            stage_name="matlab_drsa",
            processing_pass="apply",
            records=records,
            status="completed",
            qc={
                "report_json": str(paths["report_json"]),
                "report_html": str(paths["report_html"]),
                "handoff_json": str(paths["handoff_json"]),
                "png_count": report["outputs"]["png_count"],
                "mat_result_count": report["outputs"]["mat_result_count"],
                "smoothing_plots": smoothing_plots,
                "plot_smoothing_window": plot_smoothing_window,
                "model_sources": model_sources,
            },
        )
    if args.group_average and status == "completed":
        group_report = build_group_average_report(cfg, args)
        report["group_average_report"] = {
            "report_files": group_report.get("report_files", {}),
            "figure_count": len(group_report.get("outputs", [])),
            "participant_figure_count": len(
                group_report.get("participant_figures", [])
            ),
        }
    return report


def main() -> int:
    """Resolve config and execute the optional MATLAB dRSA stage."""
    args = build_arg_parser().parse_args()
    if args.resume and args.overwrite:
        raise SystemExit("--resume and --overwrite are mutually exclusive.")
    cfg, warnings, _config_path, _subject_path = load_and_resolve_config(
        args.config,
        args.subject_config,
        args.tuning_config,
        site_config_path=args.site_config,
        required_sections=("epoching", "export", "matlab_drsa"),
    )
    if warnings and args.print_summary:
        print(json.dumps({"configuration_warnings": warnings}, indent=2))
    try:
        report = run_stage(cfg, args)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"MATLAB dRSA error: {exc}")
        return 2
    if args.print_summary and report.get("stage") == "matlab_drsa_group_average":
        print(
            json.dumps(
                {
                    "status": report["status"],
                    "figure_count": len(report.get("outputs", [])),
                    "skipped_inputs": len(report.get("skipped_inputs", [])),
                    "group_drsa_dir": report.get("group_drsa_dir"),
                },
                indent=2,
                sort_keys=True,
            )
        )
    elif args.print_summary and isinstance(report.get("outputs"), dict):
        print(
            json.dumps(
                {
                    "status": report["status"],
                    "png_count": report["outputs"]["png_count"],
                    "mat_result_count": report["outputs"]["mat_result_count"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0 if report.get("status") in {"validated", "completed"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
