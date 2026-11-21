"""Run expanded descriptive diagnostics with an explicit source/result split.

Mac example from a reviewed repository root (all 30 participants)::
    python meg/code/preprocessing_basic/run_drsa_dependence_cohort.py \
      --brain-results-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/processing_runs/noica-20261117-v2-drsa-mnn-v1 \
      --diagnostics-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives/processing_runs/noica-20261117-v2-drsa-mnn-v1/rdm_dependence_expanded_mnn_DATE \
      --mode all --workers 4 --threads 1 \
      --matlab-bin /Applications/MATLAB_R2020a.app/bin/matlab --use-arch-x86-64

Inputs: accepted source Stage-06 exports, MNN Stage-07 MATs/reports and exact
production noise transforms. Add --subject-config-dir with the accepted
selectors frozen for the MNN run. Model-only RDMs are not noise-normalized.
Defaults select the accepted 30-person cohort, all five groups, nine model
sources, three conditions and +/-2 s lags. Outputs are immutable participant
bundles and a joined group report under a fresh expansion subfolder. No source
manifest or Stage-07 output is updated. Validation and independent participant
processes use a bounded process pool. Group input validation and summaries also
use --workers; shared plotting/publication is serial after the full join. Modes validate,
participants and group also support Slurm preflight/array/afterok workflows.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import importlib
import multiprocessing
import os
from pathlib import Path
import time

COHORT = "04,05,06,07,08,09,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27,28,31,32,33,36,37,38,39"
SOURCES = "observed,predicted,prediction_error,pe_prior_unconstrained_50,pe_prior_constrained_50,pe_precision_weighted,pe_kalman_innovation,pe_latent_curvature,pe_belief_update"
GROUPS = "frontal,temporal,parietal,occipital,wholebrain"


def execute(subject, options, mode):
    """Resolve reviewed tuning, validate the source run and execute one action."""
    stage = importlib.import_module("07_run_matlab_drsa")
    project = Path(__file__).resolve().parents[3]
    params = project / "meg/params"
    # Reuse the accepted selectors frozen for the production MNN run, including
    # its recovery choices. Never regenerate epochs or refit normalization.
    selectors = Path(options["subject_config_dir"]) if options["subject_config_dir"] else params
    cfg, warnings, _, _ = stage.load_and_resolve_config(
        params / "pipeline_default.yaml", selectors / f"sub-{subject}.yaml", None,
        site_config_path=options["site_config"], neural_normalization=options["neural_normalization"],
        required_sections=("epoching", "export", "matlab_drsa"),
    )
    for warning in warnings:
        print(warning, flush=True)
    cfg["resolved_paths"]["project_root"] = str(project / "meg")
    paths = stage._analysis_paths(cfg)
    from drsa_analysis_variant import source_config
    source_paths = source_config(cfg)["resolved_paths"]
    argv = ["--config", str(params / "pipeline_default.yaml"),
            "--source-manifest", str(Path(source_paths["derivatives_reports_dir"]) / "source_manifest.json"),
            "--sensor-groups", *options["groups"].split(","),
            "--model-sources", *options["sources"].split(","),
            "--conditions", "always_visible", "occluded_nondeviant", "occluded_deviant",
            "--drsa-type", "corr", "--resume"]
    if options["matlab_bin"]:
        argv += ["--matlab-bin", options["matlab_bin"]]
    if options["use_arch_x86_64"]:
        argv += ["--use-arch-x86-64"]
    if mode == "validate":
        argv += ["--validate-only"]
    if mode == "group":
        os.environ["DAD_RDM_GROUP_WORKERS"] = str(options["workers"])
        argv += ["--group-average-only", "--group-subjects", *options["subjects"]]
    started = time.monotonic()
    print(f"RDM cohort | sub-{subject} | {mode} started", flush=True)
    from drsa_dependence import diagnostics_root, run_participant, sensor_tags
    args = stage.build_arg_parser().parse_args(argv)
    args.brain_results_root = options["brain_results_root"]
    args.lag_window_sec = 2.0
    root = diagnostics_root(options["diagnostics_root"], cfg, paths)
    groups = sensor_tags(options["groups"].split(","))
    conditions = ["always_visible", "occluded_nondeviant", "occluded_deviant"]
    if mode == "group":
        from plot_drsa_dependence import run_group
        report = run_group(cfg, args, paths, root, ["sub-"+s for s in options["subjects"]], conditions, groups)
    else:
        report = run_participant(cfg, args, stage, root, conditions, groups)
    print(f"RDM cohort | sub-{subject} | {mode} {report['status']} | {time.monotonic()-started:.1f} s", flush=True)
    return report["directory"]


def parallel_subjects(options, mode):
    """Join all selected participant processes; propagate every failure."""
    with ProcessPoolExecutor(max_workers=options["workers"], mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = {pool.submit(execute, subject, options, mode): subject for subject in options["subjects"]}
        for index, future in enumerate(as_completed(futures), 1):
            future.result()
            print(f"RDM cohort | {mode} joined {index}/{len(futures)} | sub-{futures[future]}", flush=True)


def main():
    """Validate resource controls before importing numerical libraries or writing."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=("validate", "participants", "group", "all"), default="validate")
    parser.add_argument("--subjects", nargs="+", default=COHORT.split(","))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--source-run", default="noica-20261117-v2")
    parser.add_argument("--neural-normalization", choices=("mnn", "none"), default="mnn")
    parser.add_argument("--sources", default=SOURCES)
    parser.add_argument("--groups", default=GROUPS)
    parser.add_argument("--site-config")
    parser.add_argument("--subject-config-dir", help="Accepted selectors from the production MNN run.")
    parser.add_argument("--brain-results-root", required=True)
    parser.add_argument("--diagnostics-root", required=True)
    parser.add_argument("--matlab-bin")
    parser.add_argument("--use-arch-x86-64", action="store_true")
    options = vars(parser.parse_args())
    if (not options["subjects"] or len(set(options["subjects"])) != len(options["subjects"])
            or any(subject not in COHORT.split(",") for subject in options["subjects"])):
        parser.error("Choose unique two-digit participants from the accepted 30-person cohort.")
    allocation = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    if min(options["workers"], options["threads"]) < 1 or options["workers"] * options["threads"] > allocation:
        parser.error("workers * threads must fit the available CPU allocation.")
    brain = Path(options["brain_results_root"]).resolve(strict=True)
    output = Path(options["diagnostics_root"]).resolve()
    if brain not in output.parents or any(p in (brain / "analysis", brain / "reports") for p in (output, *output.parents)):
        parser.error("Diagnostics must use a separate subfolder of the expansion root.")
    os.environ.update(DAD_PROCESSING_RUN=options["source_run"], DAD_EPOCH_PREPARATION="inline_no_ica",
                      DAD_MATLAB_MAX_THREADS=str(options["threads"]),
                      DAD_DRSA_NORMALIZATION=options["neural_normalization"], DAD_DRSA_ANALYSIS_LABEL=brain.name)
    # Separate processes each own at most --threads numerical-library threads.
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[variable] = str(options["threads"])
    print(f"RDM cohort | {len(options['subjects'])} participants | {options['workers']} workers x {options['threads']} threads", flush=True)
    if options["mode"] in ("validate", "participants", "all"):
        parallel_subjects(options, "validate")
    if options["mode"] in ("participants", "all"):
        parallel_subjects(options, "participants")
    if options["mode"] in ("group", "all"):
        print(execute(options["subjects"][0], options, "group"), flush=True)


if __name__ == "__main__":
    main()
