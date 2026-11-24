"""Render the descriptive group report beside matched Stage-07 brain-model data.

Example from the repository root:
    python meg/code/preprocessing_basic/run_drsa_dependence_cohort.py \
        --mode group --workers 4 --threads 1 \
        --brain-results-root /path/to/MNN-run \
        --diagnostics-root /path/to/MNN-run/fresh-diagnostics \
        --subject-config-dir /path/to/accepted-selectors

Participant bundles and ordinary dRSA MATs are read-only inputs. The report
uses Stage-08 interval and Stage-09 sliding-window definitions. All summaries
give equal weight to participants with defined values and save their counts.
Interval lag curves also save and plot the across-participant SEM at each lag.
DAD_RDM_GROUP_WORKERS controls independent input/summary processes (defaults to
SLURM_CPUS_PER_TASK, or one outside Slurm); set numerical thread limits to one
before starting Python. All workers join before serial rendering/publication.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import multiprocessing
import html
import importlib
import json
import os
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from continuous_cleaning import digest, file_hash, file_stamp
from drsa_dependence import (DEVIANT, LEGACY_SOURCES, MODEL_LABELS, model_axes, model_sources, read_bundle, read_condition)
from drsa_dependence_inputs import brain_result_path, read_brain_matrix, validate_brain_inputs
from drsa_result_utils import occlusion_intervals, sampling_step_seconds

LAG_WINDOW_SEC = 1.0
ROW_WINDOW_SEC = 0.200
ROW_STEP_SEC = 0.025
MIN_WINDOW_SAMPLES = 3

INTERVAL_NAMES = ("Pre-occlusion", "During occlusion", "Post-occlusion")


def mean_and_support(values):
    """Average finite values over participants, preserving empty cells as NaN."""
    stack = np.asarray(values, dtype=float)
    finite = np.isfinite(stack)
    support = finite.sum(axis=0)
    mean = np.full(support.shape, np.nan)
    np.divide(np.sum(stack, axis=0, where=finite), support, out=mean, where=support > 0)
    return mean, support


def summarize(matrices, axis, intervals, lag_window_sec=LAG_WINDOW_SEC):
    """Summarize subject-level curves and surfaces at lags within the half-window."""
    stage08 = importlib.import_module("08_plot_drsa_line_summaries")
    stage09 = importlib.import_module("09_plot_drsa_sliding_lag_surfaces")
    dt = sampling_step_seconds(axis, 120.0)
    lag_samples = stage08._lag_samples(axis, lag_window_sec)
    curves, surfaces = [], []
    for matrix in matrices:
        curves.append([
            stage08.compute_interval_lag_curve(matrix, lag_samples, mask)
            for mask in intervals.values()
        ])
        _, lag_axis, row_axis, surface, _, _, _ = stage09.compute_sliding_lag_surface(
            matrix, axis, lag_window_sec=lag_window_sec, model_window_sec=ROW_WINDOW_SEC,
            model_step_sec=ROW_STEP_SEC, min_valid_samples=MIN_WINDOW_SAMPLES,
        )
        surfaces.append(surface)

    matrix, matrix_n = mean_and_support(matrices)
    curve_stack = np.asarray(curves)
    curve, curve_n = mean_and_support(curve_stack)
    # Estimate uncertainty across subjects separately for each interval and lag.
    curve_sem = np.stack([
        stage08.compute_across_subject_sem(interval_curves)
        for interval_curves in curve_stack.swapaxes(0, 1)
    ])
    surface, surface_n = mean_and_support(surfaces)
    return {
        "matrix": matrix, "matrix_n": matrix_n,
        "curves": curve, "curves_n": curve_n, "curves_sem": curve_sem, "lag_sec": lag_samples * dt,
        "surface": surface, "surface_n": surface_n,
        "row_window_sec": row_axis, "surface_lag_sec": lag_axis,
    }


@contextmanager
def group_pool(workers):
    """Spawn group workers with one numerical thread, including from old wrappers.

    Set the child environment before spawning, so NumPy/SciPy initialize with
    bounded libraries. Restore the parent environment after every worker joins.
    """
    keys = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
    previous = {key: os.environ.get(key) for key in keys}
    try:
        os.environ.update({key: "1" for key in keys})
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            yield pool
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def load_group_subject(directory, conditions, groups, sources, brain_root):
    """Read and validate one participant without writes or shared mutable state.

    Used in spawned processes for independent participant IO. Return validated
    diagnostic arrays, brain matrices and byte provenance to the joining parent.
    """
    report = read_bundle(directory)
    handoff = json.loads((directory / "handoff.json").read_text())
    data = {condition: read_condition(directory / f"{condition}.mat", condition,
            report["sensor_groups"], handoff, report.get("model_sources", LEGACY_SOURCES))
            for condition in conditions}
    cache = {}
    evidence = validate_brain_inputs(brain_root, handoff, report, conditions, sources, groups,
        cache, {c: data[c]["selected_channel_labels"] for c in conditions}) if brain_root else {}
    return report, handoff, data, cache, evidence


def summarize_panels(stacks, axis, intervals, lag_window_sec, workers=1):
    """Join independent stacks while preserving each participant reduction order.

    Numerical summaries use the unchanged Stage-08/09 functions. Workers return
    arrays; only the parent renders and publishes shared outputs. Worker count
    affects resources only, and never the participant order inside a stack.
    """
    if workers == 1:
        return {key: summarize(value, axis, intervals, lag_window_sec) for key, value in stacks.items()}
    results = {}
    with group_pool(workers) as pool:
        futures = {pool.submit(summarize, value, axis, intervals, lag_window_sec): key
                   for key, value in stacks.items()}
        for index, future in enumerate(as_completed(futures), 1):
            key = futures[future]
            results[key] = future.result()
            print(f"Group RDM dependence | summary joined {index}/{len(futures)} | {key}", flush=True)
    return results


def count_range(counts):
    """Format visible cohort support while retaining exact counts in the NPZ."""
    positive = counts[counts > 0]
    if not positive.size:
        return "n = 0"
    low, high = int(positive.min()), int(positive.max())
    return f"n = {low}" if low == high else f"n = {low}–{high}"


def plot_condition(path, condition, group, axis, timing, panels, lag_window_sec):
    """Put brain-model, model-model and neural-neural comparisons in columns."""
    pe_source = next((key.removeprefix("brain_") for key in panels
                      if key.startswith("brain_") and key not in ("brain_observed", "brain_predicted")), None)
    pe_label = MODEL_LABELS.get(pe_source)

    def axis_label(label):
        # Keep long variant names in the figure header and metadata, so axis
        # labels and neighboring titles remain readable at normal page width.
        return "PE model" if label == pe_label else label

    def title_label(label):
        return "PE" if label == pe_label else label.replace(" position", "").split(" (")[0]

    brain_limit = max(
        (
            float(np.max(np.abs(summary["matrix"]), where=np.isfinite(summary["matrix"]), initial=0))
            for key, (_, summary) in panels.items() if key.startswith("brain_")
        ),
        default=0,
    ) or 1.0
    # Keep comparable y-scales and include the full shaded bands.
    curve_limits = {True: brain_limit, False: 1.0}
    for key, (_, summary) in panels.items():
        is_brain = key.startswith("brain_")
        extent = np.abs(summary["curves"]) + np.nan_to_num(summary["curves_sem"])
        curve_limits[is_brain] = max(
            curve_limits[is_brain], np.max(extent, where=np.isfinite(extent), initial=0),
        )
    fig, axes = plt.subplots(3, len(panels), figsize=(3.45 * len(panels), 9.5), squeeze=False,
                             layout="constrained")
    boundaries = (
        [timing["first_complete_occlusion_sec"], timing["reappearance_sec"]]
        if condition != "always_visible" else []
    )
    for column, (key, ((row_label, col_label), summary)) in enumerate(panels.items()):
        is_brain = key.startswith("brain_")
        limit = brain_limit if is_brain else 1.0
        curve_limit = curve_limits[is_brain]
        heatmap, curves, surface = axes[:, column]
        heatmap.set_title(f"{title_label(row_label)} × {title_label(col_label)}", fontweight="bold", fontsize=11)
        image = heatmap.pcolormesh(
            axis, axis, summary["matrix"], cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest",
        )
        heatmap.plot(axis, axis, color="0.4", linewidth=0.7, alpha=0.6)
        heatmap.set(xlabel=f"{axis_label(col_label)} time (s)", ylabel=f"{axis_label(row_label)} time (s)", aspect="equal")
        fig.colorbar(image, ax=heatmap, shrink=0.8, label="Pearson r")
        for time in boundaries:
            heatmap.axhline(time, color="0.4", linestyle=":", lw=0.8)
            heatmap.axvline(time, color="0.4", linestyle=":", lw=0.8)

        for label, color, curve, sem, n in zip(
            INTERVAL_NAMES, ("#0072b2", "#d55e00", "#009e73"),
            summary["curves"], summary["curves_sem"], summary["curves_n"],
        ):
            curves.plot(summary["lag_sec"], curve, label=f"{label} ({count_range(n)})", color=color)
            curves.fill_between(
                summary["lag_sec"], curve - sem, curve + sem, color=color, alpha=0.20, linewidth=0,
            )
        curves.axvline(0, color="0.4", lw=0.8, linestyle=":")
        curves.set(xlim=(-lag_window_sec, lag_window_sec), ylim=(-curve_limit, curve_limit),
                   xlabel="Lag (s)", ylabel="Mean Pearson r ± SEM")
        curves.legend(loc="lower right", fontsize=7)

        surface.pcolormesh(
            summary["row_window_sec"], summary["surface_lag_sec"], summary["surface"].T,
            cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest",
        )
        surface.axhline(0, color="0.4", lw=0.8, linestyle=":")
        for time in boundaries:
            surface.axvline(time, color="0.4", lw=0.8, linestyle=":")
        surface.set(
            xlabel=f"{axis_label(row_label)} window centre (s)", ylabel="Lag (s)",
            ylim=(-lag_window_sec, lag_window_sec),
            title=f"Sliding lag · {count_range(summary['surface_n'])}",
        )

    subtitle = f"PE = {pe_label} · " if pe_label else ""
    fig.suptitle(f"{condition.replace('_', ' ').capitalize()} · {group} MEG\n{subtitle}lag = column time − row time",
                 fontweight="bold")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_group(cfg, args, paths, root, subjects, conditions, groups):
    """Join complete aligned diagnostic and brain-model bundles before publishing."""
    if args.overwrite:
        raise ValueError("Diagnostic reports are immutable. Use a fresh diagnostic root.")
    workers = int(os.environ.get("DAD_RDM_GROUP_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    allocation = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    if not 1 <= workers <= allocation:
        raise ValueError("DAD_RDM_GROUP_WORKERS must fit the CPU allocation.")
    subjects = list(dict.fromkeys(subjects))
    conditions = list(dict.fromkeys(conditions))
    session = cfg["layout"]["session_bids_id"]
    sources = model_sources(args)
    brain_root = getattr(args, "brain_results_root", None)
    analysis_root = Path(brain_root).resolve(strict=True) / "analysis" if brain_root else paths["analysis_dir"].parent.parent
    lag_window_sec = LAG_WINDOW_SEC
    group_dir = root / "group"
    if args.lag_window_sec is not None:
        lag_window_sec = args.lag_window_sec
        if not np.isfinite(lag_window_sec) or lag_window_sec <= 0:
            raise ValueError("--lag-window-sec must be positive and finite.")
        # Keep an explicitly sized report separate so the same MATs can be replotted.
        group_dir /= f"lag-{str(lag_window_sec).removesuffix('.0')}s"
    final_dir = (group_dir / session).resolve()
    if root not in final_dir.parents:
        raise ValueError("The group output path must stay inside the diagnostic root.")

    # Enforce a shared participant set and sampling/time axes for every panel.
    bundles = {subject: root / subject / session for subject in subjects}
    if any(source not in LEGACY_SOURCES for source in sources) and not brain_root:
        raise ValueError("PE reports require --brain-results-root for explicit expanded inputs.")
    reports, handoffs, diagnostic_data, brain_cache, brain_evidence = {}, {}, {}, {}, {}
    with group_pool(workers) as pool:
        futures = {pool.submit(load_group_subject, directory, conditions, groups, sources, brain_root): subject
                   for subject, directory in bundles.items()}
        for index, future in enumerate(as_completed(futures), 1):
            subject = futures[future]
            report, handoff, data, cache, evidence = future.result()
            reports[subject], handoffs[subject] = report, handoff
            diagnostic_data.update({(subject, condition): value for condition, value in data.items()})
            brain_cache.update(cache)
            if brain_root:
                brain_evidence[subject] = evidence
            print(f"Group RDM dependence | inputs joined {index}/{len(futures)} | {subject}", flush=True)
    reference = handoffs[subjects[0]]
    axis = np.asarray(reference["time_axis_sec"], dtype=float)
    timing = reference["occlusion_timing"]
    settings = {
        "lag_window_sec": lag_window_sec, "row_window_sec": ROW_WINDOW_SEC,
        "row_step_sec": ROW_STEP_SEC, "minimum_window_samples": MIN_WINDOW_SAMPLES,
        "curve_error_band": "across-subject SEM",
    }
    smoothing = {
        (report["smooth_neural"], report["smooth_window_samples"]) for report in reports.values()
    }
    normalization_modes = {report.get("neural_normalization", {"mode": "none"})["mode"] for report in reports.values()}
    if len(normalization_modes) != 1:
        raise ValueError("Participant diagnostics mix neural normalization methods.")
    normalization_mode = normalization_modes.pop()
    settings["neural_normalization"] = normalization_mode
    if len(smoothing) != 1:
        raise ValueError("Participant diagnostics disagree on neural smoothing.")
    if len({digest(report.get("code_sha256", {})) for report in reports.values()}) != 1:
        raise ValueError("Participant diagnostics were generated by different code revisions.")
    for subject, report in reports.items():
        handoff = handoffs[subject]
        if (
            handoff["subject_bids_id"] != subject or handoff["session_bids_id"] != session
            or handoff["identity_values"] != reference["identity_values"]
            or handoff["sampling_rate_hz"] != reference["sampling_rate_hz"]
            or handoff["time_axis_sec"] != reference["time_axis_sec"]
            or handoff["occlusion_timing"] != timing
            or report["inputs"]["processing_run"] != cfg["epoching"].get("processing_run")
            or not set(groups) <= set(report["sensor_groups"])
            or not set(conditions) <= set(report["conditions"])
            or not set(sources) <= set(report.get("model_sources", LEGACY_SOURCES))
        ):
            raise ValueError(f"Diagnostic inputs do not match the requested group: {subject}.")
        for condition in conditions:
            expected_axes = model_axes(condition, report.get("model_sources", LEGACY_SOURCES))
            expected_axes.update({f"neural_{group}": (f"MEG ({group})", f"MEG ({group})")
                                  for group in report["sensor_groups"]})
            if report["matrix_axes"][condition] != {
                key: {"row": labels[0], "column": labels[1]} for key, labels in expected_axes.items()
            }:
                raise ValueError(f"Diagnostic panel axes disagree with model definitions: {subject}.")

    brain_paths = {
        (subject, condition, source, group): brain_result_path(
            analysis_root, subject, session, condition, source, group
        )
        for subject in subjects for condition in conditions for group in groups
        for source in (sources if condition == DEVIANT else ("observed",))
    }
    if not brain_root:
        # Legacy inputs still receive full axes/range/channel validation before
        # --validate-only can succeed or the shared staging directory is made.
        for (subject, condition, source, group), path in brain_paths.items():
            data = diagnostic_data[(subject, condition)]
            brain_cache[path] = read_brain_matrix(path, subject, condition, source, group, axis,
                data["identities"], data["selected_channel_labels"][group])
    # All participant validators have joined before contract creation or writes.
    brain_bindings = {
        str(path.relative_to(analysis_root)): file_stamp(path) for path in brain_paths.values()
    }
    renderer_code = {name: file_hash(Path(__file__).with_name(name)) for name in (
        Path(__file__).name, "drsa_dependence.py", "drsa_dependence_inputs.py", "08_plot_drsa_line_summaries.py",
        "09_plot_drsa_sliding_lag_surfaces.py", "drsa_result_utils.py",
    )}
    contract = digest({
        "subjects": subjects, "conditions": conditions, "groups": groups, "model_sources": sources,
        "brain_provenance": brain_evidence,
        "diagnostic_bundles": {s: r["sha256"] for s, r in reports.items()},
        "brain_results": brain_bindings, "plot_settings": settings,
        "code": renderer_code,
    })
    report = {
        "stage": "rdm_dependence_group", "status": "validated", "directory": str(final_dir),
        "subjects": subjects, "conditions": conditions, "sensor_groups": groups,
        "processing_run": cfg["epoching"].get("processing_run"),
        "plot_settings": settings,
        "contract_sha256": contract,
        "brain_results": brain_bindings, "brain_provenance": brain_evidence,
        "brain_results_root": str(analysis_root.parent), "model_sources": sources,
        "code_sha256": renderer_code,
        "diagnostic_bundle_contracts": {s: r["contract_sha256"] for s, r in reports.items()},
        "panel_axes": {}, "execution": {"group_workers": workers, "shared_publication": "serial after full join"},
    }
    if args.validate_only:
        return report
    if final_dir.exists():
        if not args.resume:
            raise FileExistsError(f"Diagnostic group report exists: {final_dir}. Use --resume.")
        previous = json.loads((final_dir / "report.json").read_text())
        if previous["contract_sha256"] != contract or any(
            file_hash(final_dir / name) != sha for name, sha in previous["sha256"].items()
        ):
            raise ValueError("The group report has changed inputs or outputs. Use a fresh root.")
        return previous

    # Serial publication follows participant validation. Model summaries are
    # reused across sensor groups, while each brain/MEG summary is group-specific.
    intervals = occlusion_intervals(axis, timing)
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".group-", dir=final_dir.parent) as staging:
        working = Path(staging)
        pages = []
        for condition in conditions:
            print(f"Group RDM dependence | {condition} | {len(subjects)} participants", flush=True)
            data = {subject: diagnostic_data[(subject, condition)] for subject in subjects}
            axes_meta = reports[subjects[0]]["matrix_axes"][condition]
            model_keys = list(model_axes(condition, sources))
            # Submit independent model/MEG/brain stacks together. Each stack
            # keeps the cohort in its original order for deterministic SEM.
            stacks = {key: [data[subject]["matrices"][key] for subject in subjects] for key in model_keys}
            for group in groups:
                stacks[f"neural_{group}"] = [data[subject]["matrices"][f"neural_{group}"] for subject in subjects]
                for source in (sources if condition == DEVIANT else ("observed",)):
                    stacks[f"brain_{source}_{group}"] = [
                        brain_cache[brain_paths[(subject, condition, source, group)]]
                        for subject in subjects]
            summaries = summarize_panels(stacks, axis, intervals, lag_window_sec, workers)
            models = {key: ((axes_meta[key]["row"], axes_meta[key]["column"]), summaries[key]) for key in model_keys}
            for group in groups:
                meg_key, meg_label = f"neural_{group}", f"MEG ({group}; {normalization_mode.upper()})"
                brain_panels = {f"brain_{source}": ((MODEL_LABELS[source], meg_label), summaries[f"brain_{source}_{group}"])
                                for source in (sources if condition == DEVIANT else ("observed",))}
                neural_panel = ((meg_label, meg_label), summaries[meg_key])
                # Preserve the original overview and its NPZ keys. Each PE gets
                # a focused page instead of a single unreadable 35-column plot.
                selections = {"overview": {
                    **{key: value for key, value in brain_panels.items() if key in ("brain_observed", "brain_predicted")},
                    **{key: models[key] for key in model_axes(condition)}, meg_key: neural_panel,
                }}
                if condition == DEVIANT:
                    for source in sources:
                        if source in LEGACY_SOURCES:
                            continue
                        selections[source] = {
                            f"brain_{source}": brain_panels[f"brain_{source}"],
                            **{key: models[key] for key in (f"{source}_self", f"observed_vs_{source}", f"predicted_vs_{source}")},
                            meg_key: neural_panel,
                        }
                for selection, panels in selections.items():
                    stem = f"{condition}_{group}" + (f"_{selection}" if selection != "overview" else "")
                    print(f"Group RDM dependence | {group} | {selection}", flush=True)
                    plot_condition(working / f"{stem}.png", condition, group, axis, timing, panels, lag_window_sec)
                    report["panel_axes"][stem] = {
                        key: {"row": labels[0], "column": labels[1]} for key, (labels, _) in panels.items()
                    }
                    np.savez_compressed(
                        working / f"{stem}.npz", time_axis_sec=axis, interval_names=list(intervals),
                        **{f"{key}_{name}": value for key, (_, summary) in panels.items()
                           for name, value in summary.items()},
                    )
                    label = html.escape(f"{condition} | {group} | {MODEL_LABELS.get(selection, selection)}")
                    pages.append(f'<details><summary>{label}</summary><figure>'
                                 f'<a href="{stem}.png"><img loading="lazy" src="{stem}.png" alt="{label}"></a>'
                                 f'<figcaption><a href="{stem}.npz">Data, SEM and participant counts</a></figcaption>'
                                 '</figure></details>')

        intro = (
            '<!doctype html><html lang="en"><meta charset="utf-8">'
            '<title>RDM correlations</title><style>body{font:16px system-ui;margin:2rem}'
            'figure{margin:2rem 0}img{max-width:100%;height:auto}</style>'
            '<h1>Descriptive RDM correlations</h1>'
            f'<p>Neural RDM normalization: <strong>{html.escape(normalization_mode.upper())}</strong>. Model RDMs are unchanged.</p>'
            '<p>PE pages compare each generated model with itself, observed position and predicted position, '
            'beside its saved brain-model result. Features and Euclidean distances follow A_dRSA_alt; '
            'constant RDMs remain NaN. These descriptive panels do not rank models or test significance.</p>'
            f'<p>{len(subjects)} participants · {html.escape(", ".join(subjects))}</p>'
            f'<p>Source exports: <code>{html.escape(str(cfg["epoching"].get("processing_run")))}</code><br>'
            f'Brain-model inputs: <code>{html.escape(str(analysis_root.parent))}</code><br>'
            f'Sensor groups: {html.escape(", ".join(groups))}<br>'
            f'Model sources: {html.escape(", ".join(sources))}<br>'
            '<a href="report.json">Input provenance, axes, settings and output hashes</a></p>'
            '<p>Lag is column time minus row time. Interval and sliding-window summaries anchor '
            'the row axis. Always-visible curves use matched time windows. Each figure uses a separate '
            'brain–model correlation scale; consult colorbars across PE pages. Shading shows ±1 '
            'across-participant SEM at lags with at least two participants. The linked data '
            'include SEM and participant counts for individual lags and cells.</p>'
        )
        (working / "index.html").write_text(intro + "\n".join(pages) + "</html>\n", encoding="utf-8")
        report.update(status="completed", sha256={
            path.name: file_hash(path) for path in sorted(working.iterdir())
        })
        (working / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        os.rename(working, final_dir)
    return report
