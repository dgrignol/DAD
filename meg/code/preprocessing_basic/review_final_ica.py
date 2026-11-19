#!/usr/bin/env python3
"""Independently plot, explain or refresh existing FINAL ICA review evidence.

Usage from any directory (activate the existing MNE environment first)
-------------------------------------------------------------------
    python review_final_ica.py --subject 08 --plot sources --source sub08_run01_210826.fif
    python review_final_ica.py --subject 08 --plot properties --component 0
    python review_final_ica.py --subject 08 --plot context --source sub08_run01_210826.fif
    python review_final_ica.py --subject 08 --plot task --component 0
    python review_final_ica.py --subject 08 --plot continuous --component 0
    python review_final_ica.py --refresh --subject 08 --native-only --label native-task-v1
    python review_final_ica.py --refresh --subjects 04 05 08 --label explained-v2 --jobs 3
    # Add --regenerate for corrected full-fit properties/context; no ICA fitting.
    # Cluster: add --site-config /absolute/DAD/meg/params/pipeline_cluster.yaml.

Inputs: existing final prefit model/metadata, exact review report, native annotated
inputs and accepted 4A/4B. Diagnostic task views additionally validate accepted
Stage 2/3 eligibility, onset samples and timing/application hashes. --expected-model
pins commands copied from HTML to the displayed model, including after a refit.

Outputs: GUI plots, optional --save PNGs, or archived/refreshed HTML/JSON/PNGs.
No ICA.fit, model save, production apply, tuning/review receipt, BAD edit, cleaned
FIF or production epoch export is available. Transient task epochs are diagnostic
projections BEFORE cleaning; they reuse final motion-onset anchors/window but are
not accepted Stage 5 data. Plot-browser clicks affect only a disposable ICA object.

Full-fit property displays rebuild the documented fit copies read-only and flatten
their display bookkeeping to avoid MNE 1.10's concatenated-source length bug.
legacy-first reproduces the OLD panels' first retained fit stretch explicitly.
Both use arbitrary 2 s segments, never stimulus locking. All inputs/model hashes
remain checked; a diagnostic runtime difference does not permit refitting/reuse
for production. Protected preparation source and fitted parameters must match.

Refresh validates the whole selected cohort before disjoint participant writes,
joins workers before shared publication, archives prior evidence, and refuses to
replace reports with existing final-model decisions. --jobs exposes CPU controls;
participant workers use one numerical-library thread. Later component decisions,
application QC and scientific cleaned-output acceptance remain explicit.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from functools import lru_cache
import os
from pathlib import Path
import re
import shlex
import shutil
import time

import mne
import numpy as np

from continuous_cleaning import (annotation_rows, bind_report, digest, file_hash,
    file_stamp, prepare_fit_copy, require_review, stage_configuration)
from participant_parallel import add_parallel_argument, map_participants, participant_workers
from pipeline_config import WORKING_SUBJECT_IDS, read_yaml_mapping
from preprocessing_stage_utils import (derivative_path_for_run, discover_raw_run_records,
    guarded_path, stage_output_dirs, utc_now, write_json)
import annotation_triage as triage

PROJECT = Path(__file__).resolve().parents[3]
PARAMS = PROJECT / "meg/params"


# -----------------------------------------------------------------------------
# Exact existing model/input validation; no fitting or decision persistence
# -----------------------------------------------------------------------------

def settle_review_outputs(paths, *, timeout=30.):
    """Require positive, stable generated-file stamps before binding shared evidence.

    SMB can briefly return mtime_ns=0 just after closing a newly written HTML.
    Two matching positive stamps prevent recording that transient state. This
    applies only to generated displays; existing source/model validation stays
    unchanged and is never relaxed to accept a zero or mismatching timestamp.
    """
    paths = list(map(Path, paths))
    previous = None
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        current = {str(path): file_stamp(path) for path in paths}
        if current == previous and all(row["mtime_ns"] > 0 and row["size"] > 0 for row in current.values()):
            return current
        previous = current
        time.sleep(.2)
    raise RuntimeError("Generated review files have not reached stable positive timestamps; evidence was not bound. "
        f"Last observed generated-output stamps: {previous}")


def load_stage():
    """Load existing plot/fit-copy helpers without calling the production entry point."""
    spec = importlib.util.spec_from_file_location("dad_existing_final_ica", Path(__file__).with_name("03_C_ICA.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_review(args, subject):
    """Return current config/report/paths after validating source and model bindings."""
    cfg = triage.subject_config(args, int(subject))
    print(f"sub-{int(subject):02}: checking accepted inputs and exact final model", flush=True)
    dirs = stage_output_dirs(cfg)
    path = dirs["reports"] / "03_C_ICA_review_report.json"
    payload = json.loads(path.read_text())
    if payload["subject_bids_id"] != cfg["layout"]["subject_bids_id"]:
        raise RuntimeError("Final ICA report belongs to another participant.")
    if payload["configuration_sha256"] != stage_configuration(cfg, "03_C_ICA", include_edge=False):
        raise RuntimeError("Final ICA fit configuration changed; reconcile the model before plotting.")
    require_review(cfg, "03_A_MaxFilt", dirs["reports"])
    require_review(cfg, "03_B_ArtefactAnnotation", dirs["reports"])
    model_path = dirs["ica_models"] / Path(payload["prefit_model"]).name
    metadata_path = model_path.with_name(model_path.stem + "_metadata.json")
    metadata = json.loads(metadata_path.read_text())
    model_sha = file_hash(model_path)
    if (model_sha != payload["prefit_model_sha256"] or model_sha != metadata["model_sha256"] or
        metadata["fingerprint"] != payload["fingerprint"]):
        raise RuntimeError("Final ICA model/metadata/report identity changed.")
    expected = getattr(args, "expected_model", None)
    if expected and expected != model_sha:
        raise RuntimeError("This command belongs to a different final model; reopen the current report.")
    fingerprint = payload["fingerprint"]
    if fingerprint["fit"]["protected_preparation_sha256"] != file_hash(Path(__file__).with_name("continuous_cleaning.py")):
        raise RuntimeError("Protected fit-preparation implementation changed; no diagnostic reconstruction permitted.")
    records = discover_raw_run_records(cfg)
    inputs = [derivative_path_for_run(record, dirs["sss"], "_ann.fif") for record in records]
    current = []
    for source in inputs:
        print(f"sub-{int(subject):02}: checking {source.name}", flush=True)
        current.append({"source_filename": source.name, "size_bytes": source.stat().st_size,
            "mtime_ns": source.stat().st_mtime_ns, "annotations": digest(annotation_rows(
                mne.io.read_raw_fif(source, preload=False, verbose="ERROR")))})
    if current != fingerprint["inputs"]:
        raise RuntimeError("Annotated inputs differ from the exact fitted model.")
    return cfg, payload, records, inputs, model_path, metadata_path


def fit_display(cfg, payload, inputs, scope):
    """Reconstruct full fit display or exactly the first legacy fit stretch, read-only."""
    from ica_review_plots import continuous_display_copy
    if scope == "legacy-first":
        summary = payload["fit_input_summaries"][0]
        interval = next(row for row in summary["intervals"] if row["status"] == "included")
        raw = mne.io.read_raw_fif(inputs[0], preload=False, verbose="ERROR")
        sfreq = raw.info["sfreq"]
        piece = raw.copy().crop((interval["start_raw"] - raw.first_samp) / sfreq,
            (interval["stop_raw_exclusive"] - raw.first_samp - 1) / sfreq)
        data, _ = prepare_fit_copy(piece, cfg["ica"])
        expected = interval["fit_samples"]
    else:
        data, summaries = load_stage().prepare_ica_fit_data(inputs, cfg["ica"])
        expected = payload["fit_execution"]["n_samples_fitted"]
        if [row["fit_samples"] for row in summaries] != [row["fit_samples"] for row in payload["fit_input_summaries"]]:
            raise RuntimeError("Reconstructed fit sample counts differ from fitted source maps.")
    if data.n_times != expected:
        raise RuntimeError("Reconstructed display duration differs from fitted evidence.")
    return continuous_display_copy(data)


def validate_full_property_evidence(cfg, payload):
    """Reuse fit pictures only when their current stamps and full-sample scope match."""
    dirs = stage_output_dirs(cfg)
    root = Path(cfg['resolved_paths']['derivatives_root'])
    expected = {row['relative_path']: row['stamp'] for row in payload.get('evidence_files', [])}
    if payload.get('diagnostic_scope') != 'full_fit' or not payload.get('component_properties_requested'):
        raise RuntimeError('Native-only refresh requires complete full-fit properties; use --regenerate first.')
    for component in range(payload['n_components_fitted']):
        paths = list(dirs['figures'].glob(f'03_C_ICA_review_ICA{component:03}_properties_*.png'))
        if len(paths) != 1 or file_stamp(paths[0]) != expected.get(str(paths[0].relative_to(root))):
            raise RuntimeError(f'Existing ICA{component:03} property evidence changed; use --regenerate.')
        sidecar = json.loads(paths[0].with_suffix('.json').read_text())
        if sidecar.get('scope') != 'full_fit' or sidecar.get('display_samples') != payload['fit_execution']['n_samples_fitted']:
            raise RuntimeError(f'Existing ICA{component:03} property scope does not match the fitted data.')


# -----------------------------------------------------------------------------
# Per-figure explanations and model-bound local/cluster reproduction commands
# -----------------------------------------------------------------------------

@lru_cache(maxsize=1)
def command_roots():
    """Read the two documented project roots once per report worker."""
    local = Path(read_yaml_mapping(PARAMS / "pipeline_default.yaml")["project"]["root"]).parent
    cluster = Path(read_yaml_mapping(PARAMS / "pipeline_cluster.yaml")["project"]["root"]).parent
    return local, cluster


def command_pair(cfg, payload, flags):
    """Quote complete independent commands for both documented mount points."""
    local, cluster = command_roots()
    common = ["--subject", str(cfg["subject"]["id"]), "--expected-model", payload["prefit_model_sha256"], *flags]
    return {"local": shlex.join(["python3", str(local / "meg/code/preprocessing_basic/review_final_ica.py"), *common]),
        "cluster": shlex.join(["python", str(cluster / "meg/code/preprocessing_basic/review_final_ica.py"),
            "--site-config", str(cluster / "meg/params/pipeline_cluster.yaml"), *common])}


def make_plot_recipes(cfg, payload, images):
    """Describe every current review figure using exact model and source/window metadata."""
    from ica_review_plots import context_selection
    records = {record.source_stem: record for record in discover_raw_run_records(cfg)}
    contexts = {}
    result = {}
    for path in map(Path, images):
        metadata = json.loads(path.with_suffix(".json").read_text()) if path.with_suffix(".json").is_file() else {}
        if "_continuous_" in path.name or "_task_" in path.name:
            component = int(re.search(r"_ICA(\d+)_", path.name).group(1))
            kind = "continuous" if "_continuous_" in path.name else "task"
            flags = ["--plot", kind, "--component", str(component)]
            description = ("Native continuous final-IC examples: 20 s per acquisition, including labelled BAD periods. "
                "First ocular annotation / fixed fallback; native acquisition seconds. Window mean subtracted only for display; shared y-scale across runs. "
                "Examples are not representative statistics; scroll the full acquisitions using the commands below."
                if kind == "continuous" else
                "Accepted Stage 2/3 motion-onset-locked native IC trials and condition means, before cleaning. "
                "Configured -0.5/+2.667 s, BAD rejection, detrend and baseline; no sensor amplitude rejection on IC AU. "
                "Grey = mean +/- 1.96 SEM; not a significance test. Full trial identities/drop reasons are in the shared manifest. "
                "Both brain and artifacts may be task-locked; absence of a mean response does not exclude neural origin.")
        elif "_properties_" in path.name:
            component = int(re.search(r"_ICA(\d+)_", path.name).group(1))
            full = metadata.get("scope") == "full_fit" or payload.get("diagnostic_scope") == "full_fit"
            scope = "full" if full else "legacy-first"
            first = next(i for i in payload["fit_input_summaries"][0]["intervals"] if i["status"] == "included")
            count = payload["fit_execution"]["n_samples_fitted"] if full else first["fit_samples"]
            description = ("Arbitrary 2 s segments of the full protected fit copy; not stimulus-locked. " if full else
                "LEGACY: MNE concatenation metadata limited this panel to the FIRST retained fit stretch, not the whole fit. ")
            description += (f"Display input {count} samples at 200 Hz. Segment samples span 0-1.995 s, not -0.5/+2.667 s. "
                "Black curve = segment mean; grey = 95% parametric CI. PSD grey band = 1 SD. "
                "The old Dropped segments percentage omits annotation/boundary rejection and is not trial retention.")
            flags = ["--plot", "properties", "--component", str(component), "--scope", scope]
            metadata.update(scope="full_fit" if full else "legacy_first_fit_stretch", display_samples=count)
        elif "_topographies_" in path.name:
            page = int(re.search(r"_topographies_(\d+)", path.name).group(1))
            flags = ["--plot", "topographies", "--page", str(page)]
            description = "Final-model spatial mixing maps. IC sign/scale are arbitrary; inspect mag and grad maps and combine with temporal evidence."
        elif "_native_context" in path.name:
            stem = path.name.split("03_C_ICA_review_", 1)[1].rsplit("_native_context.png", 1)[0]
            record = records[stem]
            if not metadata:
                raw = mne.io.read_raw_fif(derivative_path_for_run(record, stage_output_dirs(cfg)["sss"], "_ann.fif"),
                    preload=False, verbose="ERROR")
                metadata, _ = context_selection(raw, record)
            contexts[record.source_filename] = metadata
            flags = ["--plot", "context", "--source", record.source_filename]
            labels = sorted({row["description"] for row in metadata["bad_annotations"]})
            description = (f"Native unfiltered final-IC projections. Selected first ocular example: {metadata['anchor_category']}; "
                "start 2 s before that annotation, show 5 s, clamped to recording; absent eyes use fixed fallback. "
                f"BAD categories in this window: {', '.join(labels) or 'none'}. "
                "Each IC is standardized within this window: colours do not rank strength across ICs/windows. "
                "This is qualitative association, not a BAD-enrichment statistic. MISC008 is excluded from eye context.")
            if metadata.get("frontal_magnetometer_display", "").startswith("per_channel_window_mean_subtracted"):
                flags += ["--center-frontal"]
                description += " Frontal traces have only their window mean subtracted for display (pT); no temporal filtering."
        elif "_sources" in path.name:
            flags = ["--plot", "sources", "--data", "fit"]
            description = ("Scrollable final sources on concatenated protected fit copies. Time is compressed fit time, not native acquisition or stimulus time. "
                "The static screenshot shows the initial browser window; use native run sources to inspect BAD periods.")
        else:
            flags = ["--plot", "existing", "--figure", path.name]
            description = "Reopen the exact saved figure. This command views existing evidence without recomputing its processing."
        result[path.name] = {"description": description, "metadata": metadata,
            "commands": command_pair(cfg, payload, flags)}
        if "_continuous_" in path.name:
            result[path.name]["native_browser_commands"] = [
                {"source": record.source_filename, **command_pair(cfg, payload,
                    ["--plot", "sources", "--source", record.source_filename, "--component", str(component)])}
                for record in records.values()]
        if "_properties_" in path.name:
            result[path.name]["task_commands"] = command_pair(cfg, payload,
                ["--plot", "task", "--component", str(component)])
            if not full:
                result[path.name]["full_commands"] = command_pair(cfg, payload,
                    ["--plot", "properties", "--component", str(component), "--scope", "full"])
    payload["native_context_windows"] = contexts
    payload["review_commands"] = {
        "native_sources": command_pair(cfg, payload, ["--plot", "sources", "--source", next(iter(records.values())).source_filename]),
        "task_example": command_pair(cfg, payload, ["--plot", "task", "--component", "0"]),
        "full_properties_example": command_pair(cfg, payload, ["--plot", "properties", "--component", "0", "--scope", "full"])}
    return result


# -----------------------------------------------------------------------------
# Transient task-locked IC diagnostics; never accepted or exported Stage 5 data
# -----------------------------------------------------------------------------

def task_figure(cfg, model, records, inputs, component, include_bad=False):
    """Use the same batch extraction/plotting as default report task panels."""
    from final_ica_native_diagnostics import extract_native_diagnostics, task_figure_from_data
    data = extract_native_diagnostics(cfg, model, records, inputs,
        components=[component], include_bad=include_bad)
    return task_figure_from_data(cfg, data, component)


# -----------------------------------------------------------------------------
# Archived report-only refresh, with parallel participant barriers
# -----------------------------------------------------------------------------

def refresh_participant(task):
    """Archive one unreviewed participant report and refresh its human interface."""
    args, subject, phase = task
    print(f"sub-{int(subject):02}: {phase} existing final ICA evidence", flush=True)
    cfg, payload, records, inputs, model_path, metadata_path = validate_review(args, subject)
    if cfg.get("stage4_review", {}).get("03_C_ICA", {}).get("status") == "accepted" or cfg["ica"].get("reviewed_model_sha256"):
        raise RuntimeError("Saved final-model decisions exist; preserve their bound reports and use a separate diagnostic view.")
    dirs = stage_output_dirs(cfg)
    archive = dirs["reports"] / "history/ica_review" / args.label
    if archive.exists():
        raise FileExistsError(f"Immutable report archive exists: {archive}")
    if getattr(args, 'native_only', False):
        validate_full_property_evidence(cfg, payload)
    if phase == "validate":
        return {"subject": int(subject), "model_sha256": payload["prefit_model_sha256"]}
    from final_ica_review_book import write_component_book
    from render_stage03_mne_reports import render_stage_report
    images = sorted(dirs["figures"].glob("03_C_ICA_review_sources*.png"))
    images += sorted(dirs["figures"].glob("03_C_ICA_review_topographies_*.png"))
    images += sorted(dirs["figures"].glob("03_C_ICA_review_ICA*_properties_*.png"))
    images += sorted(dirs["figures"].glob("03_C_ICA_review_*_native_context.png"))
    images += sorted(dirs["figures"].glob("03_C_ICA_review_ICA*_continuous_*.png"))
    images += sorted(dirs["figures"].glob("03_C_ICA_review_ICA*_task_*.png"))
    archive.mkdir(parents=True)
    previous = {}
    for path in [dirs["reports"] / f"03_C_ICA_review_{suffix}" for suffix in ("report.json", "report.html", "components.html")]:
        if path.exists():
            shutil.copy2(path, archive / path.name)
            previous[path.name] = file_hash(path)
    native_refresh = args.regenerate or getattr(args, "native_only", False)
    if native_refresh:
        (archive / "figures").mkdir()
        for path in images:
            if not args.regenerate and not any(tag in path.name for tag in ("_continuous_", "_task_")):
                continue  # Unchanged fit/context PNGs remain embedded in the archived book.
            shutil.copy2(path, archive / "figures" / path.name)
            if path.with_suffix(".json").is_file():
                shutil.copy2(path.with_suffix(".json"), archive / "figures" / path.with_suffix(".json").name)
        manifest = dirs["figures"] / "03_C_ICA_review_native_diagnostics.json"
        if manifest.exists():
            shutil.copy2(manifest, archive / "figures" / manifest.name)
        model = mne.preprocessing.read_ica(model_path, verbose="ERROR")
    if args.regenerate:
        display = fit_display(cfg, payload, inputs, "full")
        stage = load_stage()
        images = stage.generate_component_diagnostics(model=model, fit_data=display,
            figure_dir=dirs["figures"], subject_label=cfg["layout"]["subject_bids_id"], processing_pass="review",
            raw_root=cfg["resolved_paths"]["data_root"], show=False, include_properties=True)
        for record, path in zip(records, inputs, strict=True):
            raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
            images.append(stage.save_component_context(model, raw, record, dirs["figures"], raw_root=cfg["resolved_paths"]["data_root"]))
        payload["diagnostic_scope"] = "full_fit"
        payload["component_properties_requested"] = True
        del display
    if native_refresh:
        from final_ica_native_diagnostics import generate_native_diagnostics
        images = [p for p in images if not any(tag in p.name for tag in ("_continuous_", "_task_"))]
        native_images, native_manifest = generate_native_diagnostics(cfg, model, records, inputs,
            dirs["figures"], payload["prefit_model_sha256"])
        images.extend(native_images)
        payload["native_component_diagnostics"] = {"manifest": "03_C_ICA_review_native_diagnostics.json",
            "components_complete": native_manifest["components_complete"], "continuous_and_task_default": True}
    payload["review_interface_version"] = 3 if native_refresh else payload.get("review_interface_version", 2)
    payload["diagnostic_scope"] = payload.get("diagnostic_scope", "legacy_first_fit_stretch")
    payload["plot_recipes"] = make_plot_recipes(cfg, payload, images)
    payload["diagnostic_refresh"] = {"generated_utc": utc_now(), "label": args.label,
        "model_refitted": False, "model_sha256": payload["prefit_model_sha256"],
        "mne_display_version": mne.__version__, "mne_fit_version": payload["fingerprint"]["fit"]["mne_version"],
        "previous_reports_sha256": previous, "figures_regenerated": bool(args.regenerate),
        "native_diagnostics_regenerated": native_refresh}
    book = dirs["reports"] / "03_C_ICA_review_components.html"
    payload["component_review_coverage"] = write_component_book(payload, images, book, cfg["resolved_paths"]["data_root"])
    extra = [dirs["figures"] / payload["native_component_diagnostics"]["manifest"]] if payload.get("native_component_diagnostics") else []
    settle_review_outputs(images + extra + [book])
    payload["component_book_sha256"] = file_hash(book)
    bind_report(payload, inputs + [model_path, metadata_path] + images + extra + [book],
        {"fingerprint": payload["fingerprint"], "model_sha256": payload["prefit_model_sha256"],
         "plot_recipes": payload["plot_recipes"], "component_book_sha256": payload["component_book_sha256"]},
        dirs["reports"], cfg["resolved_paths"]["derivatives_root"])
    write_json(dirs["reports"] / "03_C_ICA_review_report.json", payload, raw_root=cfg["resolved_paths"]["data_root"])
    render_stage_report(derivatives_root=Path(cfg["resolved_paths"]["derivatives_root"]),
        raw_root=cfg["resolved_paths"]["data_root"], subject_id=cfg["layout"]["subject_bids_id"],
        processing_pass="review", stage_name="ica_manual", overwrite=True, primary_report=True)
    if file_hash(model_path) != payload["prefit_model_sha256"]:
        raise RuntimeError("Final model changed during report refresh.")
    return {"subject": int(subject), "model_sha256": payload["prefit_model_sha256"], "book": str(book),
        "diagnostic_scope": payload["diagnostic_scope"], "status": "reports_updated_not_accepted"}


def main(argv=None):
    """Run an independent display or an explicitly selected report-only cohort refresh."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", default=str(PARAMS / "pipeline_default.yaml"))
    p.add_argument("--site-config")
    p.add_argument("--subject-config-dir", default=str(PARAMS))
    p.add_argument("--subject", type=int)
    p.add_argument("--subjects", nargs="+", type=int)
    p.add_argument("--expected-model")
    p.add_argument("--plot", choices=("sources", "raw", "properties", "topographies", "context", "continuous", "task", "existing"))
    p.add_argument("--component", type=int)
    p.add_argument("--source")
    p.add_argument("--data", choices=("native", "fit"), default="native")
    p.add_argument("--scope", choices=("full", "legacy-first"), default="full")
    p.add_argument("--page", type=int, default=1)
    p.add_argument("--category")
    p.add_argument("--center-frontal", action="store_true",
        help="Context only: subtract each frontal channel's window mean for display. Default reproduces legacy DC levels.")
    p.add_argument("--include-bad", action="store_true")
    p.add_argument("--figure")
    p.add_argument("--save")
    p.add_argument("--refresh", action="store_true")
    p.add_argument("--regenerate", action="store_true")
    p.add_argument("--native-only", action="store_true", help="Add/refresh default continuous and task views, retaining validated full-fit properties.")
    p.add_argument("--label")
    add_parallel_argument(p)
    args = p.parse_args(argv)
    import matplotlib
    if args.refresh or args.save:
        os.environ["MPLBACKEND"] = "Agg"
        matplotlib.use("Agg", force=True)
    elif matplotlib.get_backend().casefold() == "agg":
        import sys
        matplotlib.use("MacOSX" if sys.platform == "darwin" else "QtAgg", force=True)
    import matplotlib.pyplot as plt
    if args.refresh:
        if args.regenerate and args.native_only:
            p.error("Choose --regenerate or --native-only, not both.")
        if not args.label or not re.fullmatch(r"[A-Za-z0-9_-]+", args.label):
            p.error("refresh requires a new simple --label")
        subjects = sorted(set(args.subjects or ([args.subject] if args.subject else map(int, WORKING_SUBJECT_IDS))))
        if any(str(subject).zfill(2) not in WORKING_SUBJECT_IDS for subject in subjects):
            p.error("Current report refresh excludes P15 and historical P01-P03.")
        workers = participant_workers(args.jobs, len(subjects))
        first = triage.subject_config(args, subjects[0])
        output = Path(first["resolved_paths"]["derivatives_root"]) / "reports/ica_review_updates" / args.label / "summary.json"
        if output.exists():
            raise FileExistsError("Immutable report refresh label already exists.")
        print(f"Report-only refresh: {len(subjects)} participants, {workers} workers; no fitting or decisions", flush=True)
        map_participants(refresh_participant, [(args, subject, "validate") for subject in subjects], workers)
        results = map_participants(refresh_participant, [(args, subject, "render") for subject in subjects], workers)
        write_json(output, {"generated_utc": utc_now(), "purpose": "report_presentation_only", "participants": results},
                   raw_root=first["resolved_paths"]["data_root"], overwrite=False)
        print(f"Joined report refresh: {output}", flush=True)
        return 0
    if args.subject is None or args.plot is None:
        p.error("Independent plotting requires --subject and --plot.")
    cfg, payload, records, inputs, model_path, _ = validate_review(args, args.subject)
    model = mne.preprocessing.read_ica(model_path, verbose="ERROR")
    if args.component is not None and not 0 <= args.component < model.n_components_:
        p.error("Component outside the final fitted range.")
    if args.plot in ("properties", "task", "continuous") and args.component is None:
        p.error("This plot requires --component (zero-based final-model index).")
    from ica_review_plots import continuous_display_copy, describe_properties, native_context
    metadata = {"model_sha256": payload["prefit_model_sha256"], "mne_display_version": mne.__version__,
        "mne_fit_version": payload["fingerprint"]["fit"]["mne_version"]}
    if args.save:
        output = guarded_path(args.save, cfg["resolved_paths"]["data_root"])
        if output.suffix.casefold() != ".png" or output.is_relative_to(Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]).resolve()):
            p.error("Independent plots save PNG diagnostics outside production preprocessing outputs only.")
        if output.exists() or output.with_suffix(".json").exists():
            p.error("Choose a new diagnostic name: existing PNGs and JSON sidecars are preserved.")
    if args.plot in ("sources", "raw", "context"):
        matches = [i for i, record in enumerate(records) if record.source_filename == args.source] if args.source else [0]
        if len(matches) != 1:
            p.error("--source must name one exact acquisition from this participant.")
        i = matches[0]
        native = mne.io.read_raw_fif(inputs[i], preload=False, verbose="ERROR")
    if args.plot == "properties":
        display = fit_display(cfg, payload, inputs, args.scope)
        figs = model.plot_properties(display, picks=[args.component],
            psd_args={"fmax": float(display.info["lowpass"])}, show=False)
        fig = figs[0]
        metadata.update(describe_properties(fig, display.info["sfreq"]), scope=args.scope, display_samples=int(display.n_times))
    elif args.plot == "topographies":
        figs = model.plot_components(show=False)
        figs = figs if isinstance(figs, list) else [figs]
        if not 1 <= args.page <= len(figs):
            p.error("Topography page outside available pages.")
        fig = figs[args.page - 1]
        for other in figs:
            if other is not fig:
                plt.close(other)
    elif args.plot == "sources":
        # MNE's browser directly transforms the complete Raw; it does not use
        # the truncated get_sources metadata path. Keep native input lazy and
        # avoid allocating an additional full sensor-data display copy.
        display = fit_display(cfg, payload, inputs, "full") if args.data == "fit" else native
        with mne.viz.use_browser_backend("matplotlib"):
            fig = model.plot_sources(display, picks=[args.component] if args.component is not None else None,
                show=False, block=False, show_first_samp=args.data == "native",
                title="Final ICA: disposable browser; decisions are not saved")
    elif args.plot == "raw":
        picks = [ch for ch, kind in zip(native.ch_names, native.get_channel_types()) if
            (kind == "misc" and ch != "MISC008") or (kind == "mag" and ch.startswith(("MEG011", "MEG012", "MEG031", "MEG041")) and ch not in native.info["bads"])]
        with mne.viz.use_browser_backend("matplotlib"):
            fig = native.plot(picks=picks, duration=20., show=False, show_first_samp=True,
                title="Native frontal/MISC context; no edits are saved")
    elif args.plot == "context":
        output = Path(args.save) if args.save else Path(__import__("tempfile").gettempdir()) / f"dad-{args.subject:02}-native-context.png"
        guarded_path(output, cfg["resolved_paths"]["data_root"])
        _, info = native_context(model, native, records[i], output.parent, category=args.category,
                                  output=output, show=not bool(args.save), center_frontal=args.center_frontal)
        metadata.update(info)
        print(json.dumps(metadata, indent=2, default=str))
        if args.save:
            write_json(output.with_suffix(".json"), metadata, raw_root=cfg["resolved_paths"]["data_root"])
        return 0
    elif args.plot in ("task", "continuous"):
        from final_ica_native_diagnostics import extract_native_diagnostics, continuous_figure_from_data
        if args.source:
            pairs = [(r, p) for r, p in zip(records, inputs, strict=True) if r.source_filename == args.source]
            if len(pairs) != 1:
                p.error("--source must name one exact acquisition.")
            records, inputs = map(list, zip(*pairs))
        if args.plot == "task":
            fig, info = task_figure(cfg, model, records, inputs, args.component, args.include_bad)
        else:
            data = extract_native_diagnostics(cfg, model, records, inputs,
                components=[args.component], with_epochs=False)
            fig, info = continuous_figure_from_data(cfg, data, args.component)
        metadata.update(info)
    else:
        if not args.figure or Path(args.figure).name != args.figure:
            p.error("existing requires a figure basename.")
        path = stage_output_dirs(cfg)["figures"] / args.figure
        if not path.is_file():
            p.error("Saved figure is missing.")
        fig, ax = plt.subplots(figsize=(13, 10))
        ax.imshow(plt.imread(path)); ax.axis("off")
    if args.save:
        output = guarded_path(args.save, cfg["resolved_paths"]["data_root"])
        if output.suffix.casefold() != ".png" or output.is_relative_to(Path(cfg["resolved_paths"]["derivatives_preprocessing_dir"]).resolve()):
            p.error("Independent plots may save PNG diagnostics outside production preprocessing outputs only.")
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150, bbox_inches="tight")
        write_json(output.with_suffix(".json"), metadata, raw_root=cfg["resolved_paths"]["data_root"])
        plt.close(fig)
    else:
        plt.show(block=True)
    if args.plot == "sources":
        metadata["temporary_browser_marks_not_saved"] = list(model.exclude)
    print(json.dumps(metadata, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
