"""
test_render_stage03_mne_reports.py

Purpose
Verify that the report-only stage-03 renderer creates portable MNE HTML from
synthetic JSON/PNG inputs and includes reviewed tuning choices when available.
No raw data or expensive preprocessing is used.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_render_stage03_mne_reports.py -q
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import yaml


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
SCRIPT_PATH = CODE_DIR / "render_stage03_mne_reports.py"
sys.path.insert(0, str(CODE_DIR))


def _load_renderer():
    """Import the report renderer under a stable test module name."""
    spec = importlib.util.spec_from_file_location("dad_stage03_reports", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


RENDERER = _load_renderer()


def _write_png(path: Path) -> None:
    """Write one tiny valid PNG for MNE's image embedding path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(2, 1))
    axis.plot([0, 1], [0, 1])
    axis.set_title(path.stem)
    figure.savefig(path, dpi=50)
    plt.close(figure)


def _write_json(path: Path, payload: dict) -> None:
    """Write deterministic synthetic stage metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _build_synthetic_inputs(root: Path) -> tuple[Path, Path]:
    """Create one subject's minimum complete stage-03 report artifacts."""
    derivatives_root = root / "derivatives"
    raw_root = root / "raw"
    raw_root.mkdir(parents=True)
    report_dir = derivatives_root / "reports" / "sub-02" / "ses-01"
    figure_dir = (
        derivatives_root
        / "preprocessing"
        / "sub-02"
        / "ses-01"
        / "figures"
    )
    common = {
        "subject_bids_id": "sub-02",
        "processing_pass": "review",
        "generated_utc": "2026-06-11T00:00:00Z",
        "warnings": [],
    }
    run = {
        "acquisition_run_label": 1,
        "source_filename": "sub02_run01_130526.fif",
        "execution_status": "completed",
    }

    _write_json(
        report_dir / "03_A_MaxFilt_review_report.json",
        {
            **common,
            "stage": "maxwell_filter",
            "status": "completed",
            "bad_channel_policy": "union",
            "reference_acquisition_run": 1,
            "manual_bad_channels": {"add": [], "remove": [], "by_run": {}},
            "destination": {
                "policy": "central_run",
                "selected_acquisition_run": 1,
                "score_by_run_mm": {"run01": 0.0},
                "transform_sha256": "synthetic",
            },
            "runs": [
                {
                    **run,
                    "detected_bad_channels": ["MEG0111"],
                    "missing_vs_reference": [],
                    "extra_vs_reference": [],
                    "applied_bad_channels": ["MEG0111"],
                    "destination_acquisition_run": 1,
                }
            ],
        },
    )
    _write_png(
        figure_dir
        / "03_A_MaxFilt_review_sub02_run01_130526_detection_psd.png"
    )
    _write_png(
        figure_dir
        / "03_A_MaxFilt_review_sub02_run01_130526_before_after_psd.png"
    )

    _write_json(
        report_dir / "03_B_ArtefactAnnotation_review_report.json",
        {
            **common,
            "stage": "artifact_annotation",
            "status": "completed",
            "default_muscle_threshold": 10.0,
            "accepted_muscle_threshold_by_run": {"run01": 9.0},
            "runs": [
                {
                    **run,
                    "muscle_threshold": 9.0,
                    "muscle_annotation_count": 2,
                    "eog_status": "skipped_no_eog_channels",
                    "detected_eog_event_count": 0,
                    "misc_channels_preserved": ["MISC001"],
                }
            ],
        },
    )
    _write_png(
        figure_dir
        / "03_B_ArtefactAnnotation_review_sub02_run01_130526_muscle_zscore.png"
    )

    _write_json(
        report_dir / "03_C_ICA_review_report.json",
        {
            **common,
            "stage": "ica_manual",
            "status": "review_complete",
            "cache_status": "reused_compatible",
            "n_components_fitted": 1,
            "manual_exclude": [0],
            "eog_diagnostics": {
                "status": "skipped_no_configured_eog_channels",
                "aggregate_ranking": [],
                "per_file": [],
            },
            "fit_input_summaries": [],
        },
    )
    _write_png(figure_dir / "03_C_ICA_review_sources.png")
    _write_png(figure_dir / "03_C_ICA_review_topographies_01.png")
    _write_png(figure_dir / "03_C_ICA_review_ICA000_properties_01.png")

    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "03_tuning_choices.yaml").write_text(
        yaml.safe_dump(
            {
                "subject": {"id": "02"},
                "maxwell_filter": {"bad_channel_policy": "reference_run"},
                "artifact_annotation": {
                    "muscle_threshold_by_run": {"run01": 9.0}
                },
                "ica": {"manual_exclude": [0]},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return derivatives_root, raw_root


def test_renderer_embeds_images_and_tuning_without_cluster_paths(tmp_path):
    """Generated reports must be portable and expose accepted decisions."""
    derivatives_root, raw_root = _build_synthetic_inputs(tmp_path)
    outputs = [
        RENDERER.render_stage_report(
            derivatives_root=derivatives_root,
            raw_root=raw_root,
            subject_id="02",
            processing_pass="review",
            stage_name=stage,
            overwrite=False,
        )
        for stage in RENDERER.STAGE_SPECS
    ]

    assert len(outputs) == 3
    for output in outputs:
        report_html = output.read_text(encoding="utf-8")
        assert "data:image/png;base64," in report_html
        assert "03 tuning choices" in report_html
        assert "manual_exclude" in report_html
        assert "/mnt/storage/" not in report_html


def test_renderer_refuses_raw_tree_output(tmp_path):
    """The shared guard must reject a derivative root nested under raw data."""
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    unsafe_root = raw_root / "derivatives"
    unsafe_root.mkdir(parents=True)

    try:
        RENDERER.render_stage_report(
            derivatives_root=unsafe_root,
            raw_root=raw_root,
            subject_id="02",
            processing_pass="review",
            stage_name="maxwell_filter",
            overwrite=False,
        )
    except ValueError as exc:
        assert "raw" in str(exc).lower()
    else:
        raise AssertionError("Expected the shared raw-tree guard to reject output")


def test_stage4_eyes_do_not_require_legacy_eog_images(tmp_path):
    """Current 4B rendering needs muscle/manual evidence, never an EOG PNG."""
    path = tmp_path / "03_B_ArtefactAnnotation_review_source_muscle_zscore.png"
    _write_png(path)
    payload = {"scientific_acceptance": "separate_explicit_stage4_review", "runs": [{
        "source_filename": "source.fif", "acquisition_run_label": 1,
        "eog_status": "legacy_not_in_stage4", "manual_eye_review": {"status": "reviewed"},
        "final_filter_preview": {"reviewed_margin_sec": None}}]}
    sections, figures = RENDERER._artifact_content(payload, tmp_path, "review")
    assert len(figures) == 1
    assert any("ocular review" in title for title, body in sections)


def test_stage4_apply_reuses_review_component_images(tmp_path):
    """Applied reports embed exact-model review evidence and new native QC."""
    for name in ["03_C_ICA_review_sources.png", "03_C_ICA_review_topographies_01.png",
                 "03_C_ICA_review_ICA000_properties_01.png", "03_C_ICA_apply_source_before_after.png",
                 "03_C_ICA_apply_source_final_filter_psd.png"]:
        _write_png(tmp_path / name)
    sections, figures = RENDERER._ica_content({"n_components_fitted": 1, "runs": []}, tmp_path, "apply")
    assert len(figures) == 5
    assert any("Applied continuous" in title for title, body in sections)


def test_primary_report_repair_embeds_figures_and_explains_limits(tmp_path):
    """Repair cluster-style HTML from saved evidence without changing its JSON.

    The raw input tree is empty. Requiring embedded data URLs rather than a
    machine-specific image path catches the reported cross-machine failure.
    """
    import re
    derivatives_root, raw_root = _build_synthetic_inputs(tmp_path)
    report_dir = derivatives_root / "reports/sub-02/ses-01"
    json_path = report_dir / "03_A_MaxFilt_review_report.json"
    payload = json.loads(json_path.read_text())
    payload["scientific_acceptance"] = "separate_explicit_stage4_review"
    payload["destination"].update(
        metric="summed_pairwise_rms_digitization_displacement_mm",
        rationale="Synthetic medoid test",
        pairwise_rms_mm={"run01": {"run01": 0.}},
        tie_candidates=[1], tie_tolerance_mm=1e-6,
    )
    _write_json(json_path, payload)
    before = json_path.read_bytes()
    output = report_dir / "03_A_MaxFilt_review_report.html"
    output.write_text('<img src="/mnt/storage/cluster/figure.png">')
    assert RENDERER.main([
        "--derivatives-root", str(derivatives_root), "--raw-root", str(raw_root),
        "--subjects", "02", "--stages", "maxwell_filter", "--primary-report", "--overwrite",
    ]) == 0
    rendered = output.read_text()
    image_sources = re.findall(r'<img[^>]+src="([^"]+)"', rendered)
    assert len(image_sources) >= 2
    assert all(source.startswith("data:image/") for source in image_sources)
    for phrase in ["first 120 seconds", "1-60 Hz", "cannot verify every transient fault",
                   "not automatic acceptance", "union", "translation and rotation",
                   "Synthetic medoid test"]:
        assert phrase in rendered
    assert json_path.read_bytes() == before
    assert not list(raw_root.iterdir())


def test_primary_renderer_is_deployed_with_stage_scripts():
    """Automatic rendering must also work in the minimal cluster deployment."""
    deployment = (CODE_DIR.parents[1] / "cluster/deploy_review_to_cluster.sh").read_text()
    assert '"meg/code/preprocessing_basic/render_stage03_mne_reports.py"' in deployment
