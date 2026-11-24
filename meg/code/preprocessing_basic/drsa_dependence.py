"""Stage-07 descriptive RDM correlations and their isolated result contract.

Example from the repository root (the standalone runner supplies Stage-07 APIs):
    python meg/code/preprocessing_basic/run_drsa_dependence_cohort.py \
        --mode participants --subjects 04 --workers 1 --threads 1 \
        --brain-results-root /path/to/MNN-run \
        --diagnostics-root /path/to/MNN-run/fresh-diagnostics \
        --subject-config-dir /path/to/accepted-selectors

MATLAB consumes the validated Stage-06 handoff. Python validates and publishes
a complete participant bundle with its own provenance, independent of the
existing brain-model analysis. The same root is used for the group report.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from continuous_cleaning import digest, file_hash, file_stamp
from drsa_result_utils import normalize_sensor_group
from preprocessing_stage_utils import guarded_path

SENSOR_GROUPS = {"wholebrain", "occipital", "parietal", "frontal", "temporal"}
DEVIANT = "occluded_deviant"
MODEL_AXES = {
    "observed_self": ("Observed position", "Observed position"),
    "predicted_self": ("Predicted position", "Predicted position"),
    "observed_vs_predicted": ("Observed position", "Predicted position"),
}
MODEL_LABELS = {
    "observed": "Observed position", "predicted": "Predicted position",
    "prediction_error": "Prediction error position",
    "pe_prior_unconstrained_50": "50/50 unconstrained-prior PE",
    "pe_prior_constrained_50": "50/50 constrained-prior PE",
    "pe_precision_weighted": "Precision-weighted PE",
    "pe_kalman_innovation": "Kalman innovation PE",
    "pe_latent_curvature": "Latent curvature PE",
    "pe_belief_update": "Belief-update PE",
}
LEGACY_SOURCES = ["observed", "predicted"]


def model_sources(args):
    """Keep legacy defaults; require both reference models for focused PE panels."""
    selected = (getattr(args, "model_sources", None)
                or ([args.model_source] if getattr(args, "model_source", None) else None)
                or os.environ.get("DAD_DRSA_MODEL_SOURCES", "observed,predicted").split(","))
    selected = list(dict.fromkeys(selected))
    if set(selected) - MODEL_LABELS.keys() or not set(LEGACY_SOURCES) <= set(selected):
        raise ValueError("Diagnostics require observed and predicted references and registered PE sources.")
    return selected


def model_axes(condition, sources=LEGACY_SOURCES):
    """Map each stored comparison to explicit row/column feature labels."""
    axes = dict(MODEL_AXES) if condition == DEVIANT else {"observed_self": MODEL_AXES["observed_self"]}
    if condition == DEVIANT:
        for source in sources:
            if source in LEGACY_SOURCES:
                continue
            label = MODEL_LABELS[source]
            axes[f"{source}_self"] = (label, label)
            axes[f"observed_vs_{source}"] = (MODEL_LABELS["observed"], label)
            axes[f"predicted_vs_{source}"] = (MODEL_LABELS["predicted"], label)
    return axes


def diagnostics_root(value, cfg, paths):
    """Require a separate tree outside both raw and ordinary dRSA outputs."""
    root = guarded_path(value, cfg["resolved_paths"]["data_root"])
    for directory in (paths["analysis_dir"], paths["reports_dir"]):
        production_tree = directory.parent.parent.resolve()
        if root == production_tree or production_tree in root.parents:
            raise ValueError("Place diagnostics outside the ordinary analysis and reports trees.")
    return root


def sensor_tags(groups):
    """Canonicalize requested groups before building filenames and MATLAB inputs."""
    tags = list(dict.fromkeys(normalize_sensor_group(group) for group in groups))
    tags = ["parietal" if tag == "pariatel" else tag for tag in tags]
    if not tags or set(tags) - SENSOR_GROUPS:
        raise ValueError(f"Unknown diagnostic sensor groups: {tags}.")
    return list(dict.fromkeys(tags))


def matrix_keys(condition, groups, sources=LEGACY_SOURCES):
    """Return the exact requested matrix set for one condition."""
    models = list(model_axes(condition, sources))
    return [*models, *(f"neural_{tag}" for tag in groups)]


def read_condition(path, condition, groups, handoff, sources=LEGACY_SOURCES):
    """Read and validate identity/time axes, correlation support and shared prefix."""
    data = loadmat(path, simplify_cells=True)
    axis = np.asarray(data["time_axis_sec"])
    if (
        data["condition"] != condition
        or not np.array_equal(axis, handoff["time_axis_sec"])
        or not np.array_equal(data["identities"], handoff["identity_values"])
    ):
        raise ValueError(f"Invalid condition, time or identity axis in {path}.")
    normalization = handoff.get("neural_normalization", {"mode": "none"})
    if normalization["mode"] == "mnn":
        for group in groups:
            saved = data.get("normalization_info", {}).get(group, {})
            if saved.get("mode") != "mnn" or saved.get("contract_sha256") != normalization["contract_sha256"]:
                raise ValueError(f"Diagnostic MNN transform does not match production: {path}:{group}.")
    matrices = data["matrices"]
    if set(matrices) != set(matrix_keys(condition, groups, sources)):
        raise ValueError(f"Unexpected set of diagnostic matrices in {path}.")
    if condition == DEVIANT and set(sources) - set(LEGACY_SOURCES):
        metadata = data.get("model_metadata", {})
        if set(metadata) != set(sources):
            raise ValueError(f"Missing model feature metadata in {path}.")
        for source, item in metadata.items():
            feature = {"pe_latent_curvature": "latent_curvature", "pe_belief_update": "belief_update"}.get(source, "position")
            count = {"pe_latent_curvature": 1, "pe_belief_update": 2,
                     "pe_prior_unconstrained_50": 3, "pe_prior_constrained_50": 3}.get(source, 2)
            if item["feature"] != feature or item["feature_count"] != count or item["rdm_distance"] != "euclidean":
                raise ValueError(f"Unexpected PE feature/RDM definition in {path}:{source}.")

    # Undefined Pearson values are represented by NaN, including on the diagonal.
    for key, matrix in matrices.items():
        if (
            matrix.shape != (axis.size, axis.size)
            or np.isinf(matrix).any()
            or (np.abs(matrix) > 1 + 1e-12).any()
        ):
            raise ValueError(f"Invalid correlation matrix {path}:{key}.")
        if "_vs_" not in key:
            diagonal = np.diag(matrix)
            if not np.allclose(matrix, matrix.T, equal_nan=True, rtol=0, atol=1e-12):
                raise ValueError(f"Nonsymmetric autocorrelation in {path}:{key}.")
            if not np.allclose(diagonal[np.isfinite(diagonal)], 1, rtol=0, atol=1e-12):
                raise ValueError(f"Non-unit autocorrelation diagonal in {path}:{key}.")

    if condition == DEVIANT:
        prefix = int(data["shared_prefix_frames"])
        if not 1 <= prefix <= axis.size:
            raise ValueError(f"Invalid shared prefix length in {path}.")
        cross_prefix = matrices["observed_vs_predicted"][:prefix, :prefix]
        for key in ("observed_self", "predicted_self"):
            auto_prefix = matrices[key][:prefix, :prefix]
            if not np.allclose(cross_prefix, auto_prefix, equal_nan=True, rtol=0, atol=1e-12):
                raise ValueError(f"Shared prefix disagrees with {key} in {path}.")
    return data


def read_bundle(directory, *, contract=None):
    """Verify the bytes of a published participant bundle before using it."""
    report = json.loads((directory / "report.json").read_text())
    if report["status"] != "completed" or (
        contract is not None and report["contract_sha256"] != contract
    ):
        raise ValueError(f"Different diagnostic inputs or settings in {directory}. Use a fresh root.")
    if set(report["sha256"]) != {f"{name}.mat" for name in report["conditions"]} | {"handoff.json"}:
        raise ValueError(f"Incomplete diagnostic bundle in {directory}.")
    for name, expected in report["sha256"].items():
        if file_hash(directory / name) != expected:
            raise ValueError(f"Diagnostic artifact changed: {directory / name}.")
    return report


def run_participant(cfg, args, stage, root, conditions, groups):
    """Validate inputs, run a single MATLAB process and publish atomically."""
    if args.overwrite:
        raise ValueError("Diagnostic bundles are immutable. Use a fresh diagnostic root.")
    handoff = stage.build_handoff(
        cfg, source_manifest=args.source_manifest, selected_analyses=["rdm_dependence"]
    )
    # Use the producer's exact saved transform. Never refit a diagnostic-specific
    # metric, nor silently combine MNN brain matrices with raw MEG self panels.
    if cfg.get("matlab_drsa", {}).get("neural_normalization", "mnn") == "mnn":
        production = Path(args.brain_results_root) / "reports" / handoff["subject_bids_id"] / handoff["session_bids_id"]
        saved_handoff = stage._load_json(production / "07_matlab_drsa_handoff.json")
        normalization = saved_handoff.get("neural_normalization", {})
        if normalization.get("mode") != "mnn":
            raise ValueError("Requested MNN diagnostics require MNN production inputs.")
        transform = Path(normalization["transforms_mat"])
        evidence = stage._load_json(transform.with_suffix(".json"))
        if (evidence["contract_sha256"] != normalization["contract_sha256"]
                or evidence["mat_sha256"] != file_hash(transform)):
            raise ValueError("Production MNN artifact has changed.")
        expected_width = int(cfg["matlab_drsa"].get("smooth_window_samples", 15)) if cfg["matlab_drsa"].get("smooth_neural", True) else 0
        if normalization.get("profiles", {}).get("alt") != expected_width:
            raise ValueError("Diagnostic temporal preprocessing differs from the production MNN profile.")
        normalization = dict(normalization, mat_sha256=file_hash(transform))
        handoff["neural_normalization"] = normalization
    else:
        handoff["neural_normalization"] = {"mode": "none"}
    repository_path = Path(__file__).resolve().parents[3]
    analysis_code = repository_path / "meg" / "code" / "analysis"
    functions = repository_path / "simulations" / "functions"
    handoff["drsa_functions_dir"] = str(functions)
    matlab = cfg["matlab_drsa"]
    smooth_neural = bool(matlab.get("smooth_neural", True))
    smooth_window = int(matlab.get("smooth_window_samples", 15))
    sources_selected = model_sources(args)
    if DEVIANT not in conditions and sources_selected != LEGACY_SOURCES:
        raise ValueError("PE diagnostics require occluded_deviant.")
    if smooth_window < 1:
        raise ValueError("The neural smoothing window must be positive.")

    # Use the checkout's code even when the site config points to a shared deployment.
    sources = [
        Path(__file__), Path(stage.__file__), analysis_code / "dRSA_dependence.m",
        analysis_code / "drsa_stage_utils.m", analysis_code / "drsa_build_alt_model.m",
        analysis_code / "megin_sensor_groups.json",
        Path(__file__).with_name("run_drsa_dependence_cohort.py"),
        Path(__file__).with_name("drsa_analysis_variant.py"),
        *sorted(functions.glob("dRSA_*.m")),
    ]
    settings = {
        "conditions": list(dict.fromkeys(conditions)),
        "sensor_groups": groups,
        "smooth_neural": smooth_neural,
        "smooth_window_samples": smooth_window,
        "neural_rdm_distance": "correlation",
        "position_rdm_distance": "euclidean",
        "rdm_similarity": "Pearson",
        "model_sources": sources_selected,
        "model_rdm_distance": "euclidean",
        "neural_normalization": handoff["neural_normalization"],
    }
    inputs = {
        "normalization": handoff["neural_normalization"],
        "observed_trajectory_sha256": file_hash(handoff["trajectory_models"]["observed"]),
        "neural_exports": {
            entry["name"]: {
                "mat": file_stamp(entry["neural_mat"]),
                "sidecar_sha256": file_hash(entry["sidecar_json"]),
            }
            for entry in handoff["condition_outputs"]
            if entry["name"] in settings["conditions"]
        },
        "identity_values": handoff["identity_values"],
        "time_axis_sec": handoff["time_axis_sec"],
        "occlusion_timing": handoff["occlusion_timing"],
        "processing_run": cfg["epoching"].get("processing_run"),
        "subject": handoff["subject_bids_id"],
        "session": handoff["session_bids_id"],
        "source_export_contract_sha256": stage._load_json(Path(handoff["export_report_json"]))["export_contract_sha256"],
    }
    source_hashes = {p.name: file_hash(p) for p in sources}
    contract = digest({"inputs": inputs, "settings": settings, "code": source_hashes})
    final_dir = (root / inputs["subject"] / inputs["session"]).resolve()
    if root not in final_dir.parents:
        raise ValueError("The participant output path must stay inside the diagnostic root.")
    report = {
        "stage": "rdm_dependence",
        "status": "validated",
        "contract_sha256": contract,
        "directory": str(final_dir),
        "inputs": inputs,
        "code_sha256": source_hashes,
        **settings,
    }
    if getattr(args, "brain_results_root", None):
        from drsa_dependence_inputs import validate_brain_inputs

        report["brain_provenance"] = validate_brain_inputs(
            args.brain_results_root, handoff, report, settings["conditions"], sources_selected, groups,
        )
    if args.validate_only:
        return report
    if final_dir.exists():
        if not args.resume:
            raise FileExistsError(f"Diagnostic bundle exists: {final_dir}. Use --resume.")
        print(f"RDM dependence | {inputs['subject']} | verified existing bundle")
        return read_bundle(final_dir, contract=contract)

    # The final directory appears only after all requested MATs pass validation.
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    quote = stage._matlab_quote
    flags = "true" if smooth_neural else "false"
    with tempfile.TemporaryDirectory(prefix=".rdm-", dir=final_dir.parent) as staging:
        working = Path(staging)
        handoff_path = working / "handoff.json"
        handoff_path.write_text(json.dumps(handoff, indent=2) + "\n")
        batch = (
            f"maxNumCompThreads({int(os.environ.get('DAD_MATLAB_MAX_THREADS', '1'))});"
            + f"addpath({quote(analysis_code)});"
            + f"dRSA_dependence({quote(handoff_path)},{quote(working)},"
            + f"{stage._matlab_cellstr(settings['conditions'])},{stage._matlab_cellstr(groups)},"
            + f"{flags},{smooth_window},{stage._matlab_cellstr(sources_selected)});"
        )
        runs = stage.run_matlab_batches(
            [("rdm_dependence", batch)], matlab_bin=stage._matlab_bin_from_config(args, cfg),
            use_arch_x86_64=args.use_arch_x86_64 or bool(matlab.get("use_arch_x86_64", False)),
            print_summary=args.print_summary,
        )
        if runs[0]["status"] != "completed":
            raise RuntimeError(f"Diagnostic MATLAB process exited with {runs[0]['return_code']}.")
        matrix_metadata = {}
        for condition in settings["conditions"]:
            data = read_condition(working / f"{condition}.mat", condition, groups, handoff, sources_selected)
            axes = model_axes(condition, sources_selected)
            axes.update({f"neural_{tag}": (f"MEG ({tag})", f"MEG ({tag})") for tag in groups})
            matrix_metadata[condition] = {
                key: {"row": axes[key][0], "column": axes[key][1]}
                for key in data["matrices"]
            }
            if condition == DEVIANT:
                report["shared_prefix_frames"] = int(data["shared_prefix_frames"])
        report.update(
            status="completed", matlab=runs[0], matrix_axes=matrix_metadata,
            sha256={path.name: file_hash(path) for path in sorted(working.iterdir())},
        )
        (working / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        os.rename(working, final_dir)
    return report
