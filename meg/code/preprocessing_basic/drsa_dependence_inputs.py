"""Validate separately located A_dRSA_alt inputs for descriptive diagnostics.

Usage (normally called by Stage 07, not as a standalone command)::
    evidence = validate_brain_inputs(run_root, handoff, diagnostic_report,
                                     conditions, sources, sensor_groups)

The run root contains analysis/ and reports/. Source handoffs continue to name
validated Stage-06 exports in the source run. Historical v0 reports have no
analysis-contract/artifact hashes: their successful execution records and saved
handoffs are checked explicitly, and current input bytes are bound by SHA-256.
No files are written. Mount prefixes may differ between DiPSEN and macOS; only
paths below the known processing_runs/ or data/ anchors are compared portably.
"""
from __future__ import annotations

import json
import hashlib
import io
import re
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from continuous_cleaning import file_hash, file_stamp
from drsa_dependence import DEVIANT, MODEL_LABELS


def brain_result_path(analysis_root, subject, session, condition, source, group):
    """Resolve the existing A_dRSA_alt correlation MAT for a known selection."""
    folder = analysis_root / subject / session / "drsa" / "matlab_real_neural" / "A_dRSA_alt" / condition
    if condition == DEVIANT:
        folder = folder / source
    prefix = f"{condition}_{source}" if condition == DEVIANT else condition
    return folder / "corr" / "results" / f"A_dRSA_alt_{prefix}_corr_{group}.mat"


def portable_path(value):
    """Compare saved source locations across known cluster/Mac mount prefixes."""
    value = str(value)
    for anchor in ("/processing_runs/", "/data/"):
        if anchor in value:
            return anchor + value.split(anchor, 1)[1]
    return value


def validate_handoff(saved, current):
    """Require the same participant, ordered axes, timing and source files."""
    keys = ("subject_bids_id", "session_bids_id", "identity_values", "time_axis_sec",
            "sampling_rate_hz", "channel_labels", "occlusion_timing")
    if any(saved[key] != current[key] for key in keys):
        raise ValueError("Brain-result handoff differs from the validated source axes or timing.")
    expected = current.get("neural_normalization", {"mode": "none"})
    actual = saved.get("neural_normalization", {"mode": "none"})
    if expected["mode"] != actual["mode"] or (expected["mode"] == "mnn" and
            (expected["contract_sha256"] != actual.get("contract_sha256") or
             portable_path(expected["transforms_mat"]) != portable_path(actual.get("transforms_mat")))):
        raise ValueError("Brain and diagnostic MNN normalization contracts differ.")
    for key in ("observed",):
        if portable_path(saved["trajectory_models"][key]) != portable_path(current["trajectory_models"][key]):
            raise ValueError("Brain-result trajectory source differs from Stage 06.")
    for key in ("source_manifest", "export_report_json"):
        if portable_path(saved[key]) != portable_path(current[key]):
            raise ValueError("Brain results reference a different Stage-06 source run.")
    old = {entry["name"]: entry for entry in saved["condition_outputs"]}
    new = {entry["name"]: entry for entry in current["condition_outputs"]}
    if old.keys() != new.keys() or any(
        portable_path(old[c][key]) != portable_path(new[c][key])
        for c in new for key in ("neural_mat", "sidecar_json")
    ):
        raise ValueError("Brain results reference different neural exports.")


def read_brain_matrix(path, subject, condition, source, group, axis, identities,
                      selected_channels=None, handoff=None, content=None):
    """Validate model rows / MEG columns and preserve undefined Pearson values."""
    if not path.is_file():
        raise FileNotFoundError(f"Required brain-model MAT is missing: {path}")
    data = loadmat(io.BytesIO(content) if content is not None else path, simplify_cells=True, variable_names=[
        "condition_str", "dRSA_neural_position", "time_axis_sec", "handoff",
        "model_source", "model_feature", "sensor_group_tag", "selected_channel_labels", "params", "normalization_info",
    ])
    matrix = np.asarray(data["dRSA_neural_position"])
    channels = data.get("selected_channel_labels")
    # R2020a saves some string arrays as opaque MCOS. Their exact channel
    # equality is unavailable to scipy; the sensor tag + handoff channel axis
    # are checked instead, and this limitation is recorded in the report.
    comparable = channels is not None and np.asarray(channels).dtype.names is None
    expected_feature = {"pe_latent_curvature": "latent_curvature", "pe_belief_update": "belief_update"}.get(source, "position")
    if (
        data["condition_str"] != condition or data["handoff"]["subject_bids_id"] != subject
        or matrix.shape != (axis.size, axis.size)
        or not np.array_equal(axis, data["time_axis_sec"])
        or not np.array_equal(identities, data["handoff"]["identity_values"])
        or data["model_source"] != source or data["sensor_group_tag"] != group
        or data.get("model_feature", "position") != expected_feature
        or np.isinf(matrix).any() or (np.abs(matrix) > 1 + 1e-12).any()
        or (selected_channels is not None and comparable and not np.array_equal(
            np.asarray(selected_channels).reshape(-1), np.asarray(channels).reshape(-1)))
    ):
        raise ValueError(f"Brain-model MAT does not match diagnostic input axes: {path}.")
    if handoff is not None:
        # scipy squeezes MATLAB singleton arrays. Reuse the JSON-shaped source
        # handoff after normalizing numeric/cell fields back to Python lists.
        saved = data["handoff"]
        for key in ("identity_values", "time_axis_sec", "channel_labels"):
            saved[key] = np.asarray(saved[key]).reshape(-1).tolist()
        entries = saved["condition_outputs"]
        saved["condition_outputs"] = [entries] if isinstance(entries, dict) else list(entries)
        validate_handoff(saved, handoff)
        expected = handoff.get("neural_normalization", {"mode": "none"})
        actual = data.get("normalization_info", {"mode": "none"})
        if expected["mode"] != actual["mode"] or (expected["mode"] == "mnn" and
                expected["contract_sha256"] != actual.get("contract_sha256")):
            raise ValueError(f"Brain-model MAT normalization differs from diagnostics: {path}.")
        params = data["params"]
        if (params["dRSAtype"] != "corr" or params["neuralDistMeasure"] != "correlation"
                or str(params["modelDistMeasure"]) != "euclidean" or "normalization" in params):
            raise ValueError(f"Unexpected RDM definition in {path}.")
    return matrix


def validate_brain_inputs(run_root, handoff, diagnostic, conditions, sources, groups,
                          matrix_cache=None, selected_channels=None):
    """Preflight all MATs; optionally retain validated arrays for group reuse.

    Reading each MAT's bytes once supplies both the digest and scipy reader.
    The optional caller-owned cache avoids rereading mounted files during
    plotting; it is not part of the serializable provenance evidence.
    """
    root = Path(run_root).resolve(strict=True)
    if not (root / "analysis").is_dir() or not (root / "reports").is_dir():
        raise ValueError("--brain-results-root must contain analysis/ and reports/.")
    subject, session = handoff["subject_bids_id"], handoff["session_bids_id"]
    directory = root / "reports" / subject / session
    report_path = directory / "07_matlab_drsa_report.json"
    handoff_path = directory / "07_matlab_drsa_handoff.json"
    report = json.loads(report_path.read_text())
    validate_handoff(json.loads(handoff_path.read_text()), handoff)
    normalization = handoff.get("neural_normalization", {"mode": "none"})
    if normalization["mode"] == "mnn" and normalization.get("mat_sha256") != file_hash(normalization["transforms_mat"]):
        raise ValueError("Diagnostic production MNN artifact has changed.")
    if handoff.get("neural_normalization", {}).get("mode") == "mnn" and not report.get("analysis_contract_sha256"):
        raise ValueError("MNN diagnostics require a bound production analysis contract.")
    if (report["status"] != "completed" or report["subject_bids_id"] != subject
            or report["session_bids_id"] != session or "alt" not in report["selected_analyses"]
            or not set(sources) <= set(report["model_sources"])):
        raise ValueError(f"Incomplete or mismatched expanded Stage-07 report: {report_path}.")
    if "source_export_contract_sha256" in report and report["source_export_contract_sha256"] != diagnostic["inputs"]["source_export_contract_sha256"]:
        raise ValueError("Expansion and diagnostic source export contracts disagree.")

    # Historical v0 stores smoothing in the executed MATLAB batch, not in a
    # settings object. Check every used model source, including observed.
    records = {record["analysis"]: record for record in report["matlab_records"]}
    for source in sources:
        record = records.get(f"alt:{source}", {})
        batch = record.get("command", [""])[-1]
        flag = re.search(r"(?:^|;)doSmoothNeural=(true|false);", batch)
        window = re.search(r"(?:^|;)smoothWindowSamples=(\d+);", batch)
        if (record.get("status") != "completed" or record.get("return_code") != 0
                or not flag or not window or (flag[1] == "true") != diagnostic["smooth_neural"]
                or int(window[1]) != diagnostic["smooth_window_samples"]
                or f"path_type='{source}';" not in batch or "dRSAtypeToRun='corr';" not in batch):
            raise ValueError(f"Missing successful matching execution/smoothing record for {subject}:{source}.")
    inventory = {portable_path(path) for path in report["outputs"]["mat_results"]}
    recorded_bindings = {entry["relative_path"]: entry["stamp"] for entry in report.get("artifact_bindings", [])}
    bindings = {}
    for condition in conditions:
        selected = sources if condition == DEVIANT else ["observed"]
        for source in selected:
            if source not in MODEL_LABELS:
                raise ValueError(f"Unknown model source {source}.")
            for group in groups:
                path = brain_result_path(root / "analysis", subject, session, condition, source, group)
                relative = str(path.relative_to(root / "analysis" / subject / session))
                if recorded_bindings and file_stamp(path) != recorded_bindings.get(relative):
                    raise ValueError(f"Expanded artifact differs from its saved binding: {path}.")
                # Use the inventory tail to support the historical sibling alias
                # even when the physical expansion folder is selected on macOS.
                tail = str(path.relative_to(root / "analysis"))
                if not any(p.endswith("/analysis/" + tail) for p in inventory):
                    raise ValueError(f"Requested MAT absent from completed inventory: {path}.")
                content = path.read_bytes()
                matrix = read_brain_matrix(path, subject, condition, source, group,
                    np.asarray(handoff["time_axis_sec"]), handoff["identity_values"],
                    selected_channels=selected_channels[condition][group] if selected_channels else None,
                    handoff=handoff, content=content)
                if matrix_cache is not None:
                    matrix_cache[path] = matrix
                bindings[tail] = {"stamp": file_stamp(path), "sha256": hashlib.sha256(content).hexdigest()}
    return {
        "root": str(root), "report_sha256": file_hash(report_path),
        "handoff_sha256": file_hash(handoff_path), "matrices": bindings,
        "provenance_level": ("contract" if "source_export_contract_sha256" in report else
                             "analysis_contract_and_handoff" if report.get("analysis_contract_sha256") else
                             "historical_execution_and_handoff"),
        "analysis_contract_sha256": report.get("analysis_contract_sha256"),
        "channel_validation": "Exact labels when readable; MCOS uses sensor tag and source channel axis",
    }
