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
     --subjects 01 02 03 04 05 06 07 \
     --processing-pass review

2) Render only the ICA report for one subject:

   python /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/render_stage03_mne_reports.py \
     --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
     --subjects 02 --stages ica_manual --processing-pass review

3) Replace MNE reports created by an earlier invocation:

   python .../render_stage03_mne_reports.py \
     --derivatives-root /Volumes/MORWUR/Projects/DAMIANO/Attention/derivatives \
     --subjects 01 02 03 04 05 06 07 \
     --processing-pass review --overwrite

Inputs
- Pass-specific JSON reports under
  ``<derivatives>/reports/sub-XX/ses-01``.
- Existing stage PNGs under
  ``<derivatives>/preprocessing/sub-XX/ses-01/figures``.
- Optional ``03_tuning_choices.yaml`` beside the JSON reports.

Outputs
One self-contained MNE report per requested stage, named with the suffix
``_mne_report.html``. Existing lightweight stage HTML files are left unchanged.
When tuning YAML exists, its full contents and source path are embedded in each
report so the displayed decision state is explicit.

Data flow
stage JSON + deterministic stage PNG discovery + optional tuning YAML
-> structured summary/tables -> embedded MNE report HTML.

Safety assumptions
- The raw tree is read-only and is never inspected by this utility.
- Every report output passes
  ``pipeline_config.assert_output_path_outside_raw`` before writing.
- Missing required JSON or PNG inputs are errors rather than silently producing
  an incomplete decision report.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import mne
import yaml

from pipeline_config import (
    SUPPORTED_SUBJECT_IDS,
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
        section="Reviewed decisions",
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


def _maxwell_content(
    payload: Mapping[str, Any],
    figure_dir: Path,
    processing_pass: str,
) -> tuple[list[tuple[str, str]], list[tuple[Path, str, str]]]:
    """Build Maxwell summaries plus detection and before/after PSD figures."""
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
            ("detected_bad_channels", "Detected bads"),
            ("missing_vs_reference", "Missing vs reference"),
            ("extra_vs_reference", "Extra vs reference"),
            ("applied_bad_channels", "Applied bads"),
            ("destination_acquisition_run", "Destination run"),
        ),
    )
    figures = []
    expected_paths = []
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
        expected_paths.extend([detection_path, comparison_path])
        figures.append(
            (
                detection_path,
                f"Run {int(run['acquisition_run_label']):02d}: bad-channel PSD",
                "Bad-channel detection",
            )
        )
        figures.append(
            (
                comparison_path,
                (
                    f"Run {int(run['acquisition_run_label']):02d}: "
                    "before/after Maxwell PSD"
                ),
                "Before and after Maxwell filtering",
            )
        )
    _require_images(expected_paths, context="Maxwell QC")
    return [
        ("Stage summary", summary),
        ("Per-run bad-channel decisions", run_table),
    ], figures


def _artifact_content(
    payload: Mapping[str, Any],
    figure_dir: Path,
    processing_pass: str,
) -> tuple[list[tuple[str, str]], list[tuple[Path, str, str]]]:
    """Build artifact summaries plus muscle and available EOG diagnostics."""
    summary = _summary_table(
        [
            ("Status", payload.get("status")),
            ("Processing pass", payload.get("processing_pass")),
            ("Generated UTC", payload.get("generated_utc")),
            ("Default muscle threshold", payload.get("default_muscle_threshold")),
            (
                "Accepted thresholds",
                payload.get("accepted_muscle_threshold_by_run"),
            ),
            ("EOG policy", payload.get("eog_policy")),
            ("Trigger-150 policy", payload.get("trigger_150_policy")),
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
            ("muscle_threshold", "Muscle threshold"),
            ("muscle_annotation_count", "Muscle segments"),
            ("eog_status", "EOG status"),
            ("detected_eog_event_count", "EOG events"),
            ("misc_channels_preserved", "Preserved MISC"),
        ),
    )
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
        if run.get("eog_status") not in (None, "skipped_no_eog_channels"):
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
    return [
        ("Stage summary", summary),
        ("Per-run physiological artifact results", run_table),
    ], figures


def _ica_content(
    payload: Mapping[str, Any],
    figure_dir: Path,
    processing_pass: str,
) -> tuple[list[tuple[str, str]], list[tuple[Path, str, str]]]:
    """Build manual-ICA summaries, EOG ranking, and component diagnostics."""
    summary = _summary_table(
        [
            ("Status", payload.get("status")),
            ("Processing pass", payload.get("processing_pass")),
            ("Generated UTC", payload.get("generated_utc")),
            ("Cache status", payload.get("cache_status")),
            ("Components fitted", payload.get("n_components_fitted")),
            ("Manual exclusions", payload.get("manual_exclude")),
            (
                "EOG diagnostic status",
                (payload.get("eog_diagnostics") or {}).get("status"),
            ),
            ("Prefit model", payload.get("prefit_model")),
            ("Post-rejection model", payload.get("postreject_model")),
            ("Blink annotation policy", payload.get("blink_annotation_policy")),
            ("Warnings", payload.get("warnings")),
        ]
    )
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

    prefix = f"03_C_ICA_{processing_pass}"
    source_paths = sorted(figure_dir.glob(f"{prefix}_sources*.png"))
    topography_paths = sorted(figure_dir.glob(f"{prefix}_topographies_*.png"))
    property_paths = sorted(figure_dir.glob(f"{prefix}_ICA*_properties_*.png"))
    n_components = int(payload.get("n_components_fitted") or 0)
    if not source_paths:
        raise FileNotFoundError(f"No ICA source figure found in {figure_dir}")
    if not topography_paths:
        raise FileNotFoundError(f"No ICA topography figures found in {figure_dir}")
    if len(property_paths) != n_components:
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
    return [
        ("Stage summary", summary),
        ("ICA fit inputs", fit_table),
        ("EOG diagnostics", eog_table),
    ], figures


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
) -> Path:
    """
    Render one stage report from existing metadata and figures.

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
    output_path = _prepare_report_output(
        report_dir
        / f"{spec['output_stem']}_{processing_pass}_mne_report.html",
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
        sections, images = _maxwell_content(
            payload,
            figure_dir,
            processing_pass,
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
        default=list(SUPPORTED_SUBJECT_IDS),
        help="Subject selectors such as 01, 02, or sub-03.",
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
        "--overwrite",
        action="store_true",
        help="Replace existing *_mne_report.html outputs.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Render every requested subject/stage report and print output paths."""
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
        for subject in args.subjects:
            subject_id = normalize_subject_id(subject)
            for stage_name in args.stages:
                output = render_stage_report(
                    derivatives_root=derivatives_root,
                    raw_root=raw_root,
                    subject_id=subject_id,
                    processing_pass=args.processing_pass,
                    stage_name=stage_name,
                    overwrite=args.overwrite,
                )
                print(f"Rendered: {output}")
        return 0
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Report rendering error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
