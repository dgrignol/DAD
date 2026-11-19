#!/usr/bin/env python3
"""
render_stage03_mne_reports.py

Purpose
Render self-contained MNE HTML reports from the JSON metadata and PNG figures
already produced by the DAD stage-03 review or apply pass. This is a report-only
utility: it does not open raw FIF files, rerun Maxwell filtering, recompute
artifact annotations, fit ICA, or modify scientific derivatives.

Usage examples
1) Render all completed review reports from the locally mounted derivative
   tree:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/render_stage03_mne_reports.py \
     --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
     --subjects {04..14} {16..39} \
     --processing-pass review --jobs 8

2) Render only the ICA report for one subject:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/render_stage03_mne_reports.py \
     --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
     --subjects 02 --stages ica_manual --processing-pass review

3) Replace MNE reports created by an earlier invocation:

   python .../render_stage03_mne_reports.py \
     --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
     --subjects {04..39} \
     --processing-pass review --overwrite

4) Repair the primary reports after a completed 4A job, without reprocessing:

   python .../render_stage03_mne_reports.py \
     --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
     --subjects {04..39} --stages maxwell_filter \
     --primary-report --overwrite

Inputs
- Pass-specific JSON reports under
  ``<derivatives>/reports/sub-XX/ses-01``.
- Existing stage PNGs under
  ``<derivatives>/preprocessing/sub-XX/ses-01/figures``.
- Optional ``03_tuning_choices.yaml`` beside the JSON reports.

Outputs
One self-contained MNE report per requested stage, named with the suffix
``_mne_report.html`` by default. ``--primary-report`` selects ``_report.html``;
the processing scripts use this option automatically. Both versions embed PNGs
and explain the processing sequence, figure limitations and review criteria.
When tuning YAML exists, its full contents and source path are embedded in each
report so the displayed decision state is explicit.

Data flow
stage JSON + deterministic stage PNG discovery + optional tuning YAML
-> structured summary/tables -> embedded MNE report HTML.
Participant reports run in separate bounded processes (--jobs), with one
numerical-library thread per worker. Stages for the same participant remain
ordered. The parent joins every participant before reporting cohort completion.

Safety assumptions
- The raw tree is read-only and is never inspected by this utility.
- Every report output passes
  ``pipeline_config.assert_output_path_outside_raw`` before writing.
- Missing required JSON or PNG inputs are errors in legacy decision reports.
  Interface v2 processing views may explicitly document missing property panels;
  the component book remains incomplete and cannot establish acceptance.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import mne
import yaml

from participant_parallel import add_parallel_argument, map_participants, participant_workers

from pipeline_config import (
    WORKING_SUBJECT_IDS,
    assert_output_path_outside_raw,
    normalize_subject_id,
)


STAGE_SPECS = {
    "maxwell_filter": {
        "json_stem": "03_A_MaxFilt",
        "output_stem": "03_A_MaxFilt",
        "title": "03_A Maxwell filtering",
        "tuning_key": "maxwell_filter",
    },
    "artifact_annotation": {
        "json_stem": "03_B_ArtefactAnnotation",
        "output_stem": "03_B_ArtefactAnnotation",
        "title": "03_B Physiological artifact annotation",
        "tuning_key": "artifact_annotation",
    },
    "ica_manual": {
        "json_stem": "03_C_ICA",
        "output_stem": "03_C_ICA",
        "title": "03_C Manual ICA",
        "tuning_key": "ica",
    },
}

DEFAULT_RAW_ROOT = Path("/Volumes/MORWUR/Projects/DAMIANO/Attention/data")


# -----------------------------------------------------------------------------
# Generic report formatting and guarded input/output handling
# -----------------------------------------------------------------------------


def _read_json_mapping(path: Path) -> dict[str, Any]:
    """Read one required JSON report and require a mapping at its root."""
    if not path.is_file():
        raise FileNotFoundError(f"Stage JSON report not found: {path}")
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError(f"Stage JSON report must contain a mapping: {path}")
    return parsed


def _read_optional_tuning(path: Path) -> dict[str, Any] | None:
    """Read optional tuning YAML, preserving the distinction from no file."""
    if not path.is_file():
        return None
    parsed = yaml.safe_load(path.read_text(encoding="utf-8"))
    if parsed is None:
        return {}
    if not isinstance(parsed, dict):
        raise ValueError(f"Tuning YAML must contain a mapping: {path}")
    return parsed


def _format_value(value: Any) -> str:
    """Format nested report values compactly without losing null decisions."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, Mapping):
        return yaml.safe_dump(
            dict(value),
            default_flow_style=True,
            sort_keys=False,
        ).strip()
    if isinstance(value, (list, tuple, set)):
        return ", ".join(str(item) for item in value) if value else "(none)"
    return str(value)


def _html_table(
    rows: Iterable[Mapping[str, Any]],
    columns: Sequence[tuple[str, str]],
) -> str:
    """Build a compact escaped HTML table for insertion into an MNE report."""
    body_rows = []
    for row in rows:
        cells = "".join(
            f"<td>{html.escape(_format_value(row.get(key)))}</td>"
            for key, _label in columns
        )
        body_rows.append(f"<tr>{cells}</tr>")
    headers = "".join(
        f"<th>{html.escape(label)}</th>" for _key, label in columns
    )
    return (
        "<div style='overflow-x:auto'>"
        "<table class='table table-striped table-hover'>"
        f"<thead><tr>{headers}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table></div>"
    )


def _summary_table(rows: Sequence[tuple[str, Any]]) -> str:
    """Render key/value stage metadata using the same table styling."""
    return _html_table(
        [{"field": key, "value": value} for key, value in rows],
        (("field", "Field"), ("value", "Value")),
    )


def _require_images(paths: Iterable[Path], *, context: str) -> list[Path]:
    """Require all expected report images and return their resolved paths."""
    resolved = [path.expanduser().resolve() for path in paths]
    missing = [path for path in resolved if not path.is_file()]
    if missing:
        missing_text = "\n".join(f"- {path}" for path in missing)
        raise FileNotFoundError(
            f"Missing required {context} image(s):\n{missing_text}"
        )
    return resolved


def _prepare_report_output(
    path: Path,
    *,
    raw_root: Path,
    overwrite: bool,
) -> Path:
    """Apply the shared raw-tree guard and explicit overwrite policy."""
    guarded = assert_output_path_outside_raw(path, raw_root)
    if guarded.exists() and not overwrite:
        raise FileExistsError(
            f"MNE report already exists; pass --overwrite to replace it: {guarded}"
        )
    guarded.parent.mkdir(parents=True, exist_ok=True)
    return guarded


def _add_tuning_section(
    report: mne.Report,
    *,
    tuning: Mapping[str, Any] | None,
    tuning_path: Path,
    stage_key: str,
) -> None:
    """Embed complete and stage-specific reviewed choices when available."""
    if tuning is None:
        return
    stage_choices = tuning.get(stage_key)
    full_yaml = yaml.safe_dump(dict(tuning), sort_keys=False).strip()
    relevant_yaml = yaml.safe_dump(
        {
            "subject": tuning.get("subject"),
            stage_key: stage_choices,
        },
        sort_keys=False,
    ).strip()
    content = (
        f"<p><strong>Source:</strong> {html.escape(str(tuning_path))}</p>"
        "<h4>Choices relevant to this stage</h4>"
        f"<pre>{html.escape(relevant_yaml)}</pre>"
        "<h4>Complete tuning file</h4>"
        f"<pre>{html.escape(full_yaml)}</pre>"
    )
    report.add_html(
        content,
        title="03 tuning choices",
        section="Tuning and review records",
        tags=("tuning", "reviewed-decisions"),
    )


def _add_images(
    report: mne.Report,
    images: Iterable[tuple[Path, str, str]],
) -> None:
    """Add PNGs to the report as embedded images grouped by section."""
    for image_path, title, section in images:
        report.add_image(
            str(image_path),
            title=title,
            caption=image_path.name,
            section=section,
            tags=("qc",),
        )


# -----------------------------------------------------------------------------
# Stage-specific summaries and deterministic figure discovery
# -----------------------------------------------------------------------------


def _review_guidance(stage_name: str, payload: Mapping[str, Any]) -> str:
    """Explain Stage 4 evidence without interpreting it as scientific acceptance.

    These paragraphs describe the approved workflow, not new measurements.
    Numeric results and recorded decisions remain in the tables and tuning
    section. Historical reports retain their original scientific context.
    """
    if not payload.get("scientific_acceptance"):
        return (
            "<p>This historical report displays the saved processing metadata "
            "and figures. Computational completion does not establish scientific "
            "acceptance; interpret it using the decisions for its acquisition.</p>"
        )
    common = (
        "<p><strong>Review evidence, not automatic acceptance.</strong> "
        "The status describes software execution or missing review inputs. "
        "The user and Codex inspect independently, reconcile findings, and "
        "record an explicit decision before the next Stage 4 gate.</p>"
    )
    if stage_name == "maxwell_filter":
        return common + (
            "<h4>What was generated and what comes later</h4>"
            "<p>4A detects bad channels independently in each acquisition and "
            "performs Maxwell filtering with a common measured destination. "
            "The approved method is SSS only (st_duration=None), using static "
            "per-run HPI measurements; continuous cHPI is unavailable. "
            "Final spectrum-fit removal at 50/100 Hz belongs to 4C after "
            "annotation and ICA review; no additional final high-pass or "
            "low-pass is selected. The 1-40 Hz filter "
            "is restricted to ICA-fit copies. These are not final cleaned FIFs.</p>"
            "<h4>Bad-channel decisions</h4>"
            "<p>Detected bads combine raw-header bads and automatic noisy/flat "
            "detections. Applied bads follow the recorded harmonization policy "
            "and explicit manual additions/removals. With union, a channel "
            "detected in any run is excluded from SSS estimation in every run. "
            "Missing/extra columns compare detected lists with the reference "
            "run; they are review clues, not further rejection decisions. "
            "The channel-list reference and head-destination run have different "
            "roles and need not be the same run.</p>"
            "<h4>How to read the figures</h4>"
            "<p>The first figures show raw PSDs with one curve per sensor and "
            "location-color insets, following the pilot layout. Detected channels "
            "are overlaid with named red dashed curves; applied-only channels "
            "are amber dash-dot curves. All raw channels are included here. "
            "The final PSD figures instead exclude the applied bad-channel list "
            "on both sides and show retained channels individually, before/after "
            "SSS, with matched y scales within each sensor type.</p>"
            "<p>The heatmaps show every channel's saved noisy/flat scores over "
            "the evaluated run windows, beside scores that crossed their own "
            "limits. Red [D] outlines/names identify detected channels in this "
            "run; amber [A] identifies channels applied only through "
            "harmonization/manual choices. Header-excluded/unscored values "
            "remain missing, not zero. Crossings alone do not redefine the "
            "recorded decision: detection also depends on the required count "
            "of problematic windows. A noisy score measures how unusually large "
            "a sensor's difference is between its measured signal and the SSS "
            "reconstruction; values above that channel's limit are crossings. "
            "A flat score is signal variability; values below MNE's fixed flat "
            "threshold indicate very little variation, as with a stuck or "
            "near-zero sensor. Here noisy scores use 5-second windows, a limit "
            "of 7, and min_count=5; a single crossed cell is not by itself an "
            "automatic bad-channel decision.</p>"
            "<p>The raw trace view includes every MEG sensor in one reproducible "
            "random window (20 seconds by default), with detected channels individually named "
            "in red-family colors, applied-only channels amber and others gray. "
            "Mag/grad have separate physical scales. Only constant medians are "
            "removed for display; there is no temporal filtering, resampling or "
            "per-channel amplitude normalization. Exact recording-relative "
            "times, native sample bounds and seed are shown. A random window "
            "can miss intermittent faults; use the heatmaps to guide further "
            "inspection. Older reports may still contain the detection-average "
            "PSD until refreshed with maxwell_qc.py.</p>"
            "<p>All PSD panels use up to the first 120 seconds and 1-60 Hz, "
            "with mag/grad separate and matching Welch settings. The same "
            "native samples and retained channel names are used before/after. "
            "Check spectral changes and unexpected "
            "broad attenuation or new peaks; a reduction alone cannot establish "
            "that brain signals were preserved. A residual 50 Hz peak can "
            "remain before the later spectrum-fit step.</p>"
            "<p>These PSDs cannot verify every transient fault or the whole "
            "recording. Inspect saved per-channel/window bad-channel scores "
            "(NPZ) and justified flagged/unflagged raw/SSS signal intervals. "
            "The figures do not cover 100 Hz, resolve the 0.01 Hz transition, "
            "or establish within-run head stability.</p>"
            "<h4>Decision before 4B</h4>"
            "<p>Reconcile channel detections and harmonization, destination, "
            "SSS settings and the planned final filtering. Use session notes "
            "as context. P17 and P24 runs 1-4 retain their irreversible 0.1 Hz "
            "acquisition high-pass limitation. Record changes explicitly and "
            "rebuild affected 4A evidence before accepting it; a changed 4A "
            "invalidates downstream annotations, ICA models and outputs.</p>"
        )
    if stage_name == "artifact_annotation":
        if payload.get("annotation_review_policy") == "current_annotations_defer_detailed_review_2026-11-07":
            return common + (
                "<p>The explicit 2026-11-07 policy adopts the current validated "
                "automatic ocular/muscle masks and preserves genuine saved manual "
                "reviews, including replacements, reviewed-empty muscle and deletion "
                "provenance. Further exhaustive interval inspection is deferred. "
                "Automatic adoption is separate from manual SAVE; untouched runs "
                "remain manually unreviewed.</p>"
                "<p>Reconcile the operative FIF/CSV, exact mask identities, native "
                "preservation and all trial outcomes across the whole working cohort "
                "before recording policy acceptance for final fitting. Source/model "
                "and accepted Stage 2/3 checks remain mandatory. This pass generates "
                "no EDF guidance or SQUID proposals and reuses accepted SSS.</p>"
                "<p>Every operative BAD category, including routine ocular intervals, "
                "is excluded from protected final ICA fit copies. Preliminary "
                "frontal ICA indices never transfer to final exclusions. Final "
                "component decisions, application QC and clean-data acceptance remain "
                "later explicit steps. Spectrum-fit 50/100 Hz removal includes P38 "
                "run 1; no additional final high-pass or low-pass is selected.</p>")
        return common + (
            "<p>4B combines automatic muscle candidates with source-specific "
            "manual BAD_blink and BAD_eye_movement intervals. Inspect muscle "
            "scores above and below the threshold with corresponding signal "
            "windows before accepting a threshold. Session labels are context, "
            "not automatic exclusions.</p>"
            "<p>Ocular review uses usable MISC eye-tracker signals and frontal "
            "MEG. Optional preliminary frontal ICA is supplementary annotation "
            "evidence; it is neither the final model nor a component-removal "
            "decision. No clear component does not establish absence of ocular "
            "artifacts. EOG is outside this workflow. Missing ocular review is "
            "unreviewed; a reviewed file with zero intervals is valid.</p>"
            "<p>Final spectrum-fit processing needs boundary/short-stretch QC. "
            "No additional high-pass is selected and no FIR exposure margin is "
            "required. All BAD categories, including "
            "routine ocular intervals, are excluded from final ICA fitting. "
            "Annotation-aware filtering splits at every BAD interval and edge; "
            "fit copies discard their own filter and resampling edges. Verify "
            "the resulting fit intervals before proceeding.</p>"
            "<p>Complete and accept 4B before final 4C fitting. Changing an "
            "annotation or threshold requires rebuilding the affected annotated "
            "data and final ICA model, repeating component review, and "
            "reapplying cleaning. This does not rewrite trial eligibility.</p>"
        )
    return common + (
        "<p>4C fits MNE ICA(method='fastica', n_components=None) on admissible "
        "1-40 Hz, resampled MEG copies after accepted 4B annotation. The number "
        "of components is selected using MNE's variance/rank safeguards: there "
        "is no dataset-specific evidence that a fixed smaller number is better. "
        "The exact compatible model is applied to each native-resolution run.</p>"
        "<p>Review topographies, component spectra/properties, source traces "
        "and native MISC/frontal context together. Routine ocular BAD intervals "
        "were excluded from fitting; projection on native data can still show "
        "ocular responses. Annotation and component decisions serve different "
        "roles. A reviewed empty exclusion list is valid; a missing decision "
        "is unreviewed. Example indices from another pipeline are not DAD "
        "exclusions.</p>"
        "<p>For application, inspect matched ICA-only before/after intervals "
        "and final-filter PSDs, plus preservation and filter-edge evidence. "
        "PSD windows alone cannot validate the very-low-frequency transition. "
        "Native sampling, sample count, first_samp, recording boundaries and "
        "STI/MISC samples must remain intact. BAD intervals remain present in "
        "continuous data; later epoch rejection uses overlaps with BAD "
        "annotations and the separately accepted Stage 5 policy.</p>"
        "<p>A changed fit input or model requires new component review; a "
        "changed exclusion list requires reapplication and before/after review. "
        "Cleaning completes only after applied outputs are independently "
        "reviewed, reconciled and accepted, with a documented Stage 5 handoff.</p>"
    )


def _maxwell_content(
    payload: Mapping[str, Any],
    figure_dir: Path,
    processing_pass: str,
    plot_metadata: Mapping[str, Any] | None = None,
) -> tuple[list[tuple[str, str]], list[tuple[Path, str, str]]]:
    """Build per-run detector plots plus the retained matched SSS PSD comparison.

    New plot metadata must agree with the reported detected/applied lists.
    Historical inputs without it retain their old detection PSD until refreshed.
    """
    destination = payload.get("destination") or {}
    summary = _summary_table(
        [
            ("Status", payload.get("status")),
            ("Processing pass", payload.get("processing_pass")),
            ("Generated UTC", payload.get("generated_utc")),
            ("MNE version", payload.get("mne_version")),
            ("Bad-channel policy", payload.get("bad_channel_policy")),
            (
                "Reference acquisition run",
                payload.get("reference_acquisition_run"),
            ),
            ("Destination policy", destination.get("policy")),
            (
                "Destination acquisition run",
                destination.get("selected_acquisition_run"),
            ),
            (
                "Centrality scores (mm)",
                destination.get("score_by_run_mm"),
            ),
            (
                "Destination transform SHA256",
                destination.get("transform_sha256"),
            ),
            ("Manual bad channels", payload.get("manual_bad_channels")),
            ("Calibration file", payload.get("calibration_file")),
            ("Cross-talk file", payload.get("crosstalk_file")),
            ("Warnings", payload.get("warnings")),
        ]
    )
    runs = list(payload.get("runs") or [])
    run_table = _html_table(
        runs,
        (
            ("acquisition_run_label", "Run"),
            ("source_filename", "Source acquisition"),
            ("execution_status", "Execution"),
            ("input_header_bads", "Raw-header bads"),
            ("auto_noisy_channels", "Automatically noisy"),
            ("auto_flat_channels", "Automatically flat"),
            ("detected_bad_channels", "Detected bads"),
            ("missing_vs_reference", "Missing vs reference"),
            ("extra_vs_reference", "Extra vs reference"),
            ("applied_bad_channels", "Applied bads"),
            ("squid_jumps", "Stage 4A SQUID review/application"),
            ("destination_acquisition_run", "Destination run"),
        ),
    )
    figures = []
    expected_paths = []
    plot_runs = {run["source_filename"]: run for run in (plot_metadata or {}).get("runs", [])}
    for run in runs:
        source_stem = Path(str(run["source_filename"])).stem
        detection_path = (
            figure_dir
            / f"03_A_MaxFilt_{processing_pass}_{source_stem}_detection_psd.png"
        )
        comparison_path = (
            figure_dir
            / f"03_A_MaxFilt_{processing_pass}_{source_stem}_before_after_psd.png"
        )
        qc = plot_runs.get(run["source_filename"])
        if plot_metadata is not None:
            if qc is None:
                raise ValueError(f"Missing plot metadata for {run['source_filename']}")
            for key in ("detected_bad_channels", "applied_bad_channels"):
                if set(qc[key]) != set(run[key]):
                    raise ValueError(f"Stale QC {key}: {run['source_filename']}; refresh maxwell_qc.py")
            for plot in qc["plots"]:
                # Require a basename; portable metadata cannot redirect image
                # discovery to another subject or outside the figure directory.
                if Path(plot["filename"]).name != plot["filename"]:
                    raise ValueError("QC plot filename must be a basename.")
                path = figure_dir / plot["filename"]
                expected_paths.append(path)
                figures.append((path, f"Run {int(run['acquisition_run_label']):02d}: {plot['title']}", plot["section"]))
        else:
            expected_paths.append(detection_path)
            figures.append((detection_path,
                f"Run {int(run['acquisition_run_label']):02d}: legacy detection-average PSD",
                "Legacy detection average (refresh for channel plots)"))
        if not qc or not any(plot["filename"] == comparison_path.name for plot in qc["plots"]):
            expected_paths.append(comparison_path)
            figures.append((comparison_path,
                f"Run {int(run['acquisition_run_label']):02d}: before/after Maxwell PSD",
                "Before and after Maxwell filtering"))
    # All raw overviews appear before channel diagnostics. Matched SSS
    # comparisons form the last figure section, with no duplicated PSD images.
    order = {"Raw PSD and sensor locations": 0, "Bad-channel score heatmaps": 1,
             "Bad-channel raw traces": 2, "Before and after Maxwell filtering": 3}
    figures.sort(key=lambda item: order.get(item[2], 1))
    _require_images(expected_paths, context="Maxwell QC")
    sections = [
        ("Stage summary", summary),
        ("Per-run bad-channel decisions", run_table),
    ]
    if plot_metadata is not None:
        sections.append(("Random raw trace windows", _html_table(
            [{"source": r["source_filename"], **r["window"]} for r in plot_runs.values()],
            (("source", "Exact source"), ("start_sec", "Start (s from recording start)"),
             ("stop_sec", "Stop (exclusive, s)"), ("native_start_sample", "First native sample"),
             ("native_stop_sample_exclusive", "Native stop (exclusive)"), ("seed", "Seed")))))
    if destination.get("metric"):
        # The sum is the selection score; individual pairwise distances are
        # separate measurements and must not be described as mean head motion.
        details = _summary_table([
            ("Distance metric", destination.get("metric")),
            ("Rationale", destination.get("rationale")),
            ("Digitization points", destination.get("digitization_point_count")),
            ("Tie candidates", destination.get("tie_candidates")),
            ("Tie tolerance (mm)", destination.get("tie_tolerance_mm")),
            ("Selected transform", destination.get("transform")),
            ("References", destination.get("references")),
        ])
        pairs = [{"run": run, **distances} for run, distances in
                 (destination.get("pairwise_rms_mm") or {}).items()]
        labels = list((destination.get("pairwise_rms_mm") or {}).keys())
        details += _html_table(pairs, [("run", "Run")] + [(run, run) for run in labels])
        sections.append(("Destination geometry and rationale",
            "<p>central_run chooses the actual run with the lowest sum of "
            "pairwise RMS digitization-point displacements in device space. "
            "It includes translation and rotation. The scores are sums, not "
            "averages or continuously measured movement. Each matrix entry "
            "below is a pairwise RMS distance in mm. The lowest acquisition "
            "number resolves ties within the recorded tolerance; invalid or "
            "missing transforms stop processing rather than trigger an "
            "unreported fallback.</p>" + details))
    if payload.get("final_filter_design"):
        sections.append(("Planned final-filter support (not applied in 4A)",
                         _summary_table(list(payload["final_filter_design"].items()))))
    return sections, figures


def _artifact_content(
    payload: Mapping[str, Any],
    figure_dir: Path,
    processing_pass: str,
) -> tuple[list[tuple[str, str]], list[tuple[Path, str, str]]]:
    """Build current ocular/muscle evidence, retaining historical EOG report support."""
    summary = _summary_table(
        [
            ("Status", payload.get("status")),
            ("Processing pass", payload.get("processing_pass")),
            ("Generated UTC", payload.get("generated_utc")),
            ("Default muscle threshold", payload.get("default_muscle_threshold")),
            (
                "Reported thresholds (scientific acceptance is separate)",
                payload.get("proposed_muscle_threshold_by_run", payload.get("accepted_muscle_threshold_by_run")),
            ),
            ("Ocular support policy", payload.get("eog_policy")),
            ("Trigger-150 policy", payload.get("trigger_150_policy")),
            ("Annotation processing policy", payload.get("annotation_review_policy")),
            ("Adoption plan identity", payload.get("adoption_plan_sha256")),
            ("Warnings", payload.get("warnings")),
        ]
    )
    runs = list(payload.get("runs") or [])
    # Reviewed muscle replacements need not match threshold proposals. Keep the
    # applied count and its source visible beside automatic diagnostic counts.
    columns = [("acquisition_run_label", "Run"), ("source_filename", "Source acquisition"),
               ("execution_status", "Execution"), ("muscle_threshold", "Muscle threshold"),
               ("muscle_annotation_count", "Muscle segments"), ("muscle_review_mode", "Muscle label source"),
               ("automatic_muscle_annotation_count", "Automatic proposed segments"),
               ("misc_channels_preserved", "Preserved MISC")]
    if not payload.get("scientific_acceptance"):
        columns.extend([("eog_status", "Historical EOG status"), ("detected_eog_event_count", "EOG events")])
    run_table = _html_table(runs, columns)
    figures = []
    expected_paths = []
    for run in runs:
        run_label = int(run["acquisition_run_label"])
        source_stem = Path(str(run["source_filename"])).stem
        prefix = f"03_B_ArtefactAnnotation_{processing_pass}_{source_stem}"
        muscle_path = figure_dir / f"{prefix}_muscle_zscore.png"
        expected_paths.append(muscle_path)
        figures.append(
            (
                muscle_path,
                f"Run {run_label:02d}: muscle z-score",
                "Muscle diagnostics",
            )
        )
        if run.get("eog_status") not in (None, "skipped_no_eog_channels", "disabled", "legacy_not_in_stage4"):
            eog_path = figure_dir / f"{prefix}_eog_diagnostics.png"
            expected_paths.append(eog_path)
            figures.append(
                (
                    eog_path,
                    f"Run {run_label:02d}: EOG diagnostics",
                    "EOG diagnostics",
                )
            )
    _require_images(expected_paths, context="artifact-annotation QC")
    sections = [("Stage summary", summary), ("Per-run physiological artifact results", run_table)]
    if payload.get("scientific_acceptance"):
        current = [{"source": r.get("source_filename"), "ocular": r.get("manual_eye_review"), "edf": r.get("edf_blink_candidates"),
                    "adoption": r.get("annotation_adoption"), "serialization": r.get("serialization_outcome_changes"),
                    "filter_exposure": r.get("final_filter_preview")} for r in runs]
        sections.append(("Source-specific ocular review and final-filter exposure", _html_table(current,
                         (("source", "Source"), ("ocular", "Manual review evidence"), ("adoption", "Separate processing adoption"),
                          ("serialization", "FIF endpoint representation: exact trial deltas"), ("edf", "EDF guidance (not BAD)"), ("filter_exposure", "Filter edge exposure")))))
    return sections, figures


def _ica_content(
    payload: Mapping[str, Any],
    figure_dir: Path,
    processing_pass: str,
) -> tuple[list[tuple[str, str]], list[tuple[Path, str, str]]]:
    """Build fit evidence, delegating v2 human inspection to the component book.

    Legacy reports require every requested property image. A v2 processing view
    can document missing diagnostics explicitly; the component book keeps its
    incomplete coverage flag and final component acceptance remains blocked.
    ICA application documents prominently retain a provisional choice reminder
    when an explicit processing authorization is present in the saved report.
    """
    separate_book = payload.get("review_interface_version", 1) >= 2 and processing_pass == "review"
    summary = _summary_table(
        [
            ("Status", payload.get("status")),
            ("Processing pass", payload.get("processing_pass")),
            ("Generated UTC", payload.get("generated_utc")),
            ("Cache status", payload.get("cache_status")),
            ("Components fitted", payload.get("n_components_fitted")),
            ("Fit samples/convergence", payload.get("fit_execution")),
            ("Optimizer restart provenance (when present)", payload.get("optimizer_restart")),
            ("Current-model diagnostic coverage", payload.get("component_review_coverage")),
            ("Manual exclusions", payload.get("manual_exclude")),
            ("Prefit model", payload.get("prefit_model")),
            ("Post-rejection model", payload.get("postreject_model")),
            ("Blink annotation policy", payload.get("blink_annotation_policy")),
            ("Warnings", payload.get("warnings")),
        ]
    )
    if payload.get("processing_provenance"):
        summary = ("<p><strong>Provisional ICA choices.</strong> Damiano has not yet "
            "manually inspected these component proposals and needs to return to "
            "inspect them. The exclusions authorize processing; cleaned-output "
            "and final scientific acceptance remain separate.</p>" + summary)
    fit_table = _html_table(
        payload.get("fit_input_summaries") or [],
        (
            ("input_path", "Fit input"),
            ("n_channels", "Channels"),
            ("n_times", "Samples"),
            ("sfreq_hz", "Sampling frequency (Hz)"),
            ("bad_channels", "Bad channels"),
        ),
    )
    eog = payload.get("eog_diagnostics") or {}
    eog_table = _summary_table(
        [
            ("Status", eog.get("status")),
            ("Aggregate ranking", eog.get("aggregate_ranking")),
            ("Per-file diagnostics", eog.get("per_file")),
        ]
    )

    # Final application reuses the exact reviewed model; its component figures
    # remain review evidence. New applied plots are separate native comparisons.
    figure_pass = "review" if processing_pass == "apply" else processing_pass
    prefix = f"03_C_ICA_{figure_pass}"
    source_paths = sorted(figure_dir.glob(f"{prefix}_sources*.png"))
    topography_paths = sorted(figure_dir.glob(f"{prefix}_topographies_*.png"))
    property_paths = sorted(figure_dir.glob(f"{prefix}_ICA*_properties_*.png"))
    n_components = int(payload.get("n_components_fitted") or 0)
    if not source_paths:
        raise FileNotFoundError(f"No ICA source figure found in {figure_dir}")
    if not topography_paths:
        raise FileNotFoundError(f"No ICA topography figures found in {figure_dir}")
    if len(property_paths) != n_components and payload.get("component_properties_requested", True) and not separate_book:
        raise FileNotFoundError(
            "Expected one ICA properties figure per fitted component; "
            f"found {len(property_paths)} for {n_components} components in "
            f"{figure_dir}."
        )

    figures = [
        (path, f"ICA source traces {index}", "ICA sources")
        for index, path in enumerate(source_paths, start=1)
    ]
    figures.extend(
        (path, f"ICA topographies {index}", "ICA topographies")
        for index, path in enumerate(topography_paths, start=1)
    )
    figures.extend(
        (
            path,
            path.name.split("_properties_", maxsplit=1)[0].rsplit("_", 1)[-1],
            "ICA component properties",
        )
        for path in property_paths
    )
    for path in sorted(figure_dir.glob("03_C_ICA_review_*_native_context.png")):
        figures.append((path, path.stem, "Native MISC/frontal/component context"))
    if processing_pass == "apply":
        for path in sorted(figure_dir.glob("03_C_ICA_apply_*_before_after.png")):
            figures.append((path, path.stem, "ICA-only before/after"))
        for path in sorted(figure_dir.glob("03_C_ICA_apply_*_final_filter_psd.png")):
            figures.append((path, path.stem, "Final temporal filtering"))
        for path in sorted(figure_dir.glob("03_C_ICA_apply_*_p10_slow_qc.png")):
            figures.append((path, path.stem, "P10 native slow-contamination QC: decision open"))
    sections = [("Stage summary", summary), ("ICA fit inputs", fit_table)]
    if len(property_paths) != n_components:
        sections.append(("Incomplete component diagnostics",
            "<p><strong>Component property diagnostics are incomplete.</strong> "
            "Regenerate all current-model properties before final component review; "
            "this processing view does not accept missing diagnostics.</p>"))
    if payload.get("eog_diagnostics") and not payload.get("scientific_acceptance"):
        sections.append(("Historical EOG diagnostics", eog_table))
    if separate_book:
        # The navigator owns human component inspection. Keep this report for
        # processing, fit accounting, provenance and tuning, without a second
        # partially overlapping set of the same component figures.
        sections.insert(0, ("Start component review", '<p>Use <a href="03_C_ICA_review_components.html">'
            'the component-review book</a> for the guide, every diagnostic and independent plotting commands. '
            'This processing report records fit accounting, convergence, source/model identity and tuning. '
            'Both refer to the same exact prefit model; neither accepts component decisions.</p>'))
        return sections, []
    if processing_pass == "apply":
        sections.append(("Applied continuous preservation/filtering", _html_table(payload.get("runs", []),
            (("source_filename", "Source"), ("preservation", "Preservation"), ("final_filter", "Final filtering")))))
    return sections, figures


# -----------------------------------------------------------------------------
# Public rendering workflow and CLI
# -----------------------------------------------------------------------------


def render_stage_report(
    *,
    derivatives_root: Path,
    raw_root: Path,
    subject_id: str,
    processing_pass: str,
    stage_name: str,
    overwrite: bool,
    primary_report: bool = False,
) -> Path:
    """
    Render one stage report from existing metadata and figures, without FIF I/O.

    Set primary_report to use the processing scripts' normal _report.html name.
    The default preserves the separate report-only utility's _mne_report.html
    convention. Replacing either existing HTML requires overwrite=True.

    Returns
    -------
    pathlib.Path
        Absolute path to the generated self-contained MNE HTML report.
    """
    subject_id = normalize_subject_id(subject_id)
    spec = STAGE_SPECS[stage_name]
    subject_label = f"sub-{subject_id}"
    report_dir = derivatives_root / "reports" / subject_label / "ses-01"
    figure_dir = (
        derivatives_root
        / "preprocessing"
        / subject_label
        / "ses-01"
        / "figures"
    )
    json_path = (
        report_dir
        / f"{spec['json_stem']}_{processing_pass}_report.json"
    )
    tuning_path = report_dir / "03_tuning_choices.yaml"
    suffix = "report" if primary_report else "mne_report"
    output_path = _prepare_report_output(
        report_dir
        / f"{spec['output_stem']}_{processing_pass}_{suffix}.html",
        raw_root=raw_root,
        overwrite=overwrite,
    )

    payload = _read_json_mapping(json_path)
    if payload.get("subject_bids_id") != subject_label:
        raise ValueError(
            f"Subject mismatch in {json_path}: "
            f"{payload.get('subject_bids_id')!r} != {subject_label!r}"
        )
    if payload.get("processing_pass") != processing_pass:
        raise ValueError(
            f"Processing-pass mismatch in {json_path}: "
            f"{payload.get('processing_pass')!r} != {processing_pass!r}"
        )

    if stage_name == "maxwell_filter":
        plot_path = report_dir / f"03_A_MaxFilt_{processing_pass}_plot_metadata.json"
        plot_metadata = _read_json_mapping(plot_path) if plot_path.exists() else None
        if plot_metadata is not None and plot_metadata.get("stage_report_sha256") != hashlib.sha256(json_path.read_bytes()).hexdigest():
            raise ValueError("Maxwell plot metadata belongs to a changed stage report; refresh maxwell_qc.py.")
        sections, images = _maxwell_content(
            payload,
            figure_dir,
            processing_pass,
            plot_metadata,
        )
    elif stage_name == "artifact_annotation":
        sections, images = _artifact_content(
            payload,
            figure_dir,
            processing_pass,
        )
    else:
        sections, images = _ica_content(
            payload,
            figure_dir,
            processing_pass,
        )

    report = mne.Report(
        title=f"{subject_label} | {spec['title']} ({processing_pass})",
        verbose=False,
    )
    report.add_html(
        _review_guidance(stage_name, payload),
        title="Processing and review guide",
        section="Start here",
        tags=("summary", "review"),
    )
    for title, section_html in sections:
        report.add_html(
            section_html,
            title=title,
            section="Processing summary",
            tags=("summary",),
        )
    _add_tuning_section(
        report,
        tuning=_read_optional_tuning(tuning_path),
        tuning_path=tuning_path,
        stage_key=str(spec["tuning_key"]),
    )
    _add_images(report, images)
    report.save(
        str(output_path),
        open_browser=False,
        overwrite=overwrite,
        sort_content=False,
        verbose=False,
    )
    return output_path


def build_parser() -> argparse.ArgumentParser:
    """Build the report-only command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Render self-contained MNE reports from completed stage-03 JSON "
            "metadata and PNG figures."
        )
    )
    parser.add_argument(
        "--derivatives-root",
        required=True,
        type=Path,
        help="Derivative tree containing reports/ and preprocessing/.",
    )
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=DEFAULT_RAW_ROOT,
        help=(
            "Read-only raw root used only by the shared output guard; defaults "
            f"to {DEFAULT_RAW_ROOT}."
        ),
    )
    parser.add_argument(
        "--subjects",
        nargs="+",
        default=list(WORKING_SUBJECT_IDS),
        help="Current cohort selectors 04 through 39.",
    )
    parser.add_argument(
        "--processing-pass",
        choices=("review", "apply"),
        default="review",
        help="Pass-specific JSON and figure set to render.",
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=tuple(STAGE_SPECS),
        default=list(STAGE_SPECS),
        help="Stage reports to render; defaults to all production stage-03 steps.",
    )
    parser.add_argument(
        "--primary-report",
        action="store_true",
        help="Write the normal *_report.html path instead of *_mne_report.html.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing HTML reports at the selected output paths.",
    )
    add_parallel_argument(parser)
    return parser


def _render_participant_reports(task):
    """Render disjoint participant HTML; keep its stages ordered and print progress."""
    args, subject_id = task
    outputs = []
    for stage_name in args.stages:
        print(f"sub-{subject_id}: rendering {stage_name} {args.processing_pass} HTML", flush=True)
        output = render_stage_report(derivatives_root=args.derivatives_root,
            raw_root=args.raw_root, subject_id=subject_id, processing_pass=args.processing_pass,
            stage_name=stage_name, overwrite=args.overwrite, primary_report=args.primary_report)
        outputs.append(str(output))
    return outputs


def main(argv: Sequence[str] | None = None) -> int:
    """Join bounded participant renderers before reporting cohort completion."""
    args = build_parser().parse_args(argv)
    derivatives_root = args.derivatives_root.expanduser().resolve()
    raw_root = args.raw_root.expanduser().resolve()
    if not derivatives_root.is_dir():
        print(
            f"Report rendering error: derivative root not found: "
            f"{derivatives_root}",
            file=sys.stderr,
        )
        return 2
    try:
        subjects = sorted(set(normalize_subject_id(subject) for subject in args.subjects))
        args.derivatives_root, args.raw_root = derivatives_root, raw_root
        workers = participant_workers(args.jobs, len(subjects))
        print(f"Report rendering: {len(subjects)} participants; {workers} processes", flush=True)
        results = map_participants(_render_participant_reports, [(args, subject) for subject in subjects], workers)
        for outputs in results:
            for output in outputs:
                print(f"Rendered: {output}")
        return 0
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Report rendering error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
