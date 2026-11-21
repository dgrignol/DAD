"""
test_matlab_drsa_stage.py

Purpose
Validate the optional Stage-07 MATLAB dRSA handoff/orchestration contract
without launching MATLAB. The tests use tiny fake sidecars and trajectory-file
placeholders so they exercise path construction, metadata validation, and
expected output naming.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_matlab_drsa_stage.py -q
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.io import savemat


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))


def _load_stage07_module():
    """Load the numbered Stage-07 module for direct function tests."""
    module_path = CODE_DIR / "07_run_matlab_drsa.py"
    spec = importlib.util.spec_from_file_location("stage07_matlab_drsa", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


STAGE07 = _load_stage07_module()


def _fake_cfg(tmp_path: Path) -> dict:
    """Create a minimal resolved configuration for Stage-07 tests."""
    project_root = tmp_path / "DAD" / "meg"
    analysis_dir = tmp_path / "derivatives" / "analysis" / "sub-04" / "ses-01"
    reports_dir = tmp_path / "derivatives" / "reports" / "sub-04" / "ses-01"
    input_dir = tmp_path / "data" / "sub04_220526" / "input_files"
    for path in (project_root / "code" / "analysis", analysis_dir / "drsa", reports_dir, input_dir):
        path.mkdir(parents=True, exist_ok=True)
    (input_dir / "MovDot_Sub04_V30_rescueTraject.mat").touch()
    return {
        "subject": {"id": "04"},
        "layout": {"subject_bids_id": "sub-04", "session_bids_id": "ses-01"},
        "data_discovery": {
            "trajectory_mat_glob": "input_files/MovDot_Sub04_V*_rescueTraject.mat",
        },
        "subject_profile": {"versions": {"stimuli": "V30"}},
        "resolved_paths": {
            "project_root": str(project_root),
            "data_root": str(tmp_path / "data"),
            "derivatives_analysis_dir": str(analysis_dir),
            "derivatives_reports_dir": str(reports_dir),
            "subject_input_files_dir": str(input_dir),
        },
        "matlab_drsa": {
            "neural_normalization": "none",
            "occlusion_timing": {
                "frame_frequency_hz": 120.0,
                "first_complete_occlusion_frame": 130,
                "reappearance_frame": 191,
            }
        },
    }


def _write_stage06_fixture(cfg: dict):
    """Write tiny Stage-06-like files required by build_handoff."""
    drsa_dir = Path(cfg["resolved_paths"]["derivatives_analysis_dir"]) / "drsa"
    reports_dir = Path(cfg["resolved_paths"]["derivatives_reports_dir"])
    manifest = reports_dir / "source_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "derivative_lineage": {
                    "stages": {"export_drsa": {"status": "completed"}}
                }
            }
        ),
        encoding="utf-8",
    )
    (reports_dir / "06_export_neural_for_drsa_report.json").write_text(
        json.dumps({"status": "completed"}),
        encoding="utf-8",
    )
    for condition, code in {
        "always_visible": 71,
        "occluded_nondeviant": 81,
        "occluded_deviant": 91,
    }.items():
        stem = (
            f"sub-04_ses-01_task-Attention_condition-{condition}_"
            "desc-identityAveraged_neural"
        )
        (drsa_dir / f"{stem}.mat").touch()
        (drsa_dir / f"{stem}.json").write_text(
            json.dumps(
                {
                    "condition_name": condition,
                    "condition_code": code,
                    "identity_values": [1, 2, 3, 4],
                    "time_axis_sec": [0.0, 0.01, 0.02],
                    "sampling_rate_hz": 100.0,
                    "ordered_channel_labels": ["MEG0111", "MEG0112"],
                    "tensor_shape": [4, 2, 3],
                    "counts_by_identity": {"1": 1, "2": 1, "3": 1, "4": 1},
                    "missing_identities": [],
                }
            ),
            encoding="utf-8",
        )
    return manifest


def test_build_handoff_validates_stage06_sidecars(tmp_path):
    """The MATLAB handoff is generated from resolved paths and Stage-06 metadata."""
    cfg = _fake_cfg(tmp_path)
    manifest = _write_stage06_fixture(cfg)
    handoff = STAGE07.build_handoff(
        cfg,
        source_manifest=manifest,
        selected_analyses=["real"],
    )
    assert handoff["identity_values"] == [1, 2, 3, 4]
    assert handoff["identity_count"] == 4
    assert len(handoff["condition_outputs"]) == 3
    assert handoff["trajectory_models"]["stimuli_version"] == "V30"
    assert handoff["trajectory_models"]["observed"].endswith(
        "MovDot_Sub04_V30_rescueTraject.mat"
    )
    assert handoff["trajectory_models"]["expected_source_condition"] == (
        "occluded_nondeviant"
    )
    assert handoff["version"] == 2
    assert "predicted" not in handoff["trajectory_models"]
    assert handoff["occlusion_timing"]["first_complete_occlusion_frame"] == 130
    assert handoff["occlusion_timing"]["first_complete_occlusion_sec"] == 129 / 120


def test_expected_output_paths_cover_real_alt_and_simulations(tmp_path):
    """Resume checks know the expected PNG/MAT locations for each analysis family."""
    cfg = _fake_cfg(tmp_path)
    outputs = STAGE07.expected_outputs_for_analysis(
        cfg,
        ["real", "alt", "simulations"],
        ["always_visible"],
        ["wholebrain"],
        "corr",
        "observed",
    )
    as_text = "\n".join(str(path) for path in outputs)
    assert "matlab_real_neural/07_dRSA/always_visible/corr" in as_text
    assert "A_dRSA_alt/always_visible/corr" in as_text
    assert "simulations/model_vs_model/always_visible/corr" in as_text
    assert "simulations/neural_vs_neural/always_visible/corr" in as_text

    smoothed_outputs = STAGE07.expected_outputs_for_analysis(
        cfg,
        ["alt"],
        ["always_visible"],
        ["occipital"],
        "corr",
        "observed",
        smoothing_plots=True,
    )
    assert any(str(path).endswith("_smoothed.png") for path in smoothed_outputs)

    occluded_outputs = STAGE07.expected_outputs_for_analysis(
        cfg,
        ["alt", "simulations"],
        ["always_visible", "occluded_deviant"],
        ["occipital"],
        "corr",
        ["observed", "predicted"],
    )
    occluded_text = "\n".join(str(path) for path in occluded_outputs)
    assert "A_dRSA_alt/occluded_deviant/observed/corr" in occluded_text
    assert "A_dRSA_alt/occluded_deviant/predicted/corr" in occluded_text
    assert "simulations/model_vs_model/occluded_deviant/observed/corr" in occluded_text
    assert "simulations/model_vs_model/occluded_deviant/predicted/corr" in occluded_text
    assert "A_dRSA_alt/always_visible/predicted" not in occluded_text

    pe_outputs = STAGE07.expected_outputs_for_analysis(
        cfg,
        ["alt", "simulations"],
        ["always_visible", "occluded_deviant"],
        ["occipital"],
        "corr",
        ["pe_prior_constrained_50"],
    )
    pe_text = "\n".join(str(path) for path in pe_outputs)
    assert "A_dRSA_alt/occluded_deviant/pe_prior_constrained_50/corr" in pe_text
    assert "A_dRSA_alt/always_visible/pe_prior_constrained_50" not in pe_text
    assert "simulations/model_vs_model/occluded_deviant/pe_prior_constrained_50" not in pe_text


def test_python_group_smoothing_preserves_undefined_matrix_support():
    """
    Stage-07 group-display smoothing must not fill analysis-masked regions.

    Prediction-error matrices can be undefined before occlusion. The Python
    group-average plotting path should smooth only inside the originally finite
    support so downstream Stage-08 summaries cannot see artificial values.
    """
    matrix = np.full((8, 8), np.nan)
    matrix[4:, :] = np.arange(32, dtype=float).reshape(4, 8)

    smoothed = STAGE07._smooth_matrix_for_plot_python(  # pylint: disable=protected-access
        matrix,
        5,
    )

    assert not np.isfinite(smoothed[:4, :]).any()
    assert np.isfinite(smoothed[4:, :]).any()


def test_matlab_batch_uses_handoff_and_selected_script(tmp_path):
    """MATLAB commands pass the handoff explicitly instead of rebuilding paths."""
    cfg = _fake_cfg(tmp_path)
    batches = STAGE07._build_matlab_batches(  # pylint: disable=protected-access
        cfg,
        handoff_path=tmp_path / "handoff.json",
        analyses=["real"],
        drsa_type="corr",
        conditions=["always_visible"],
        sensor_groups=["wholeBrain"],
        model_source="observed",
        smoothing_plots=True,
        plot_smoothing_window=7,
    )
    assert len(batches) == 1
    _analysis, source = batches[0]
    assert "drsa_handoff_json=" in source
    assert "A_dRSA.m" in source
    assert "doSmoothingPlots=true" in source
    assert "plotSmoothingWindow=7" in source


def test_matlab_batches_obey_allocated_compute_threads(tmp_path, monkeypatch):
    """A participant cannot silently use more MATLAB threads than allocated CPUs."""
    cfg = _fake_cfg(tmp_path)
    kwargs = dict(handoff_path=tmp_path / "handoff.json", analyses=["real"],
        drsa_type="corr", conditions=["always_visible"], sensor_groups=["wholeBrain"],
        model_source="observed", smoothing_plots=False)
    monkeypatch.setenv("DAD_MATLAB_MAX_THREADS", "2")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    batches = STAGE07._build_matlab_batches(cfg, **kwargs)
    assert "maxNumCompThreads(2);" in batches[0][1]
    monkeypatch.setenv("DAD_MATLAB_MAX_THREADS", "3")
    with pytest.raises(ValueError, match="exceeds"):
        STAGE07._build_matlab_batches(cfg, **kwargs)
    monkeypatch.setenv("DAD_MATLAB_MAX_THREADS", "0")
    with pytest.raises(ValueError, match="positive"):
        STAGE07._build_matlab_batches(cfg, **kwargs)


def test_matlab_batch_command_uses_validated_headless_flags():
    """MATLAB keeps its JVM while Linux uses the validated software renderer."""
    matlab_bin = Path("/state/partition1/MATLAB/R2020a/bin/matlab")
    batch_code = "disp('probe')"

    linux_command = STAGE07._matlab_batch_command(  # pylint: disable=protected-access
        matlab_bin,
        batch_code,
        use_arch_x86_64=False,
        platform_name="linux",
    )
    assert linux_command == [
        str(matlab_bin),
        "-nodesktop",
        "-nosplash",
        "-nodisplay",
        "-softwareopengl",
        "-batch",
        batch_code,
    ]
    assert "-nojvm" not in linux_command

    mac_command = STAGE07._matlab_batch_command(  # pylint: disable=protected-access
        matlab_bin,
        batch_code,
        use_arch_x86_64=True,
        platform_name="darwin",
    )
    assert mac_command[:3] == ["arch", "-x86_64", str(matlab_bin)]
    assert "-softwareopengl" not in mac_command
    assert mac_command[-2:] == ["-batch", batch_code]


def test_alt_batches_run_observed_and_predicted_occluded_deviant(tmp_path):
    """Observed and predicted occluded-deviant runs are separate MATLAB batches."""
    cfg = _fake_cfg(tmp_path)
    batches = STAGE07._build_matlab_batches(  # pylint: disable=protected-access
        cfg,
        handoff_path=tmp_path / "handoff.json",
        analyses=["alt"],
        drsa_type="corr",
        conditions=["always_visible", "occluded_deviant"],
        sensor_groups=["occipital"],
        model_source=["observed", "predicted"],
    )

    assert [analysis for analysis, _source in batches] == ["alt:observed", "alt:predicted"]
    assert "conditions_to_run={'always_visible','occluded_deviant'}" in batches[0][1]
    assert "path_type='observed'" in batches[0][1]
    assert "conditions_to_run={'occluded_deviant'}" in batches[1][1]
    assert "path_type='predicted'" in batches[1][1]


def test_alt_batches_run_explicit_pe_variant_only_for_occluded_deviant(tmp_path):
    """Explicit PE variants are Stage-07 alt-only occluded-deviant sources."""
    cfg = _fake_cfg(tmp_path)
    batches = STAGE07._build_matlab_batches(  # pylint: disable=protected-access
        cfg,
        handoff_path=tmp_path / "handoff.json",
        analyses=["alt", "simulations"],
        drsa_type="corr",
        conditions=["always_visible", "occluded_deviant"],
        sensor_groups=["occipital"],
        model_source=["pe_prior_constrained_50"],
    )

    assert [analysis for analysis, _source in batches] == [
        "alt:pe_prior_constrained_50"
    ]
    assert "conditions_to_run={'occluded_deviant'}" in batches[0][1]
    assert "path_type='pe_prior_constrained_50'" in batches[0][1]


def test_group_average_report_averages_result_mats(tmp_path):
    """Group-average reporting uses numeric MAT arrays, not rendered PNG pixels."""
    cfg = _fake_cfg(tmp_path)
    result_rel = (
        Path("matlab_real_neural")
        / "A_dRSA_alt"
        / "always_visible"
        / "corr"
        / "results"
        / "A_dRSA_alt_always_visible_corr_occipital.mat"
    )
    for index, subject in enumerate(("sub-01", "sub-02"), start=1):
        subject_drsa = (
            tmp_path
            / "derivatives"
            / "analysis"
            / subject
            / "ses-01"
            / "drsa"
        )
        result_path = (
            subject_drsa / result_rel
        )
        result_path.parent.mkdir(parents=True, exist_ok=True)
        savemat(
            result_path,
            {
                "dRSA_neural_position": np.full((4, 4, 1), float(index)),
                "time_axis_sec": np.array([0.0, 0.1, 0.2, 0.3]),
            },
            do_compression=True,
        )
        figure_dir = result_path.parents[1] / "matrices"
        figure_dir.mkdir(parents=True, exist_ok=True)
        for suffix in ("", "_smoothed"):
            (
                figure_dir
                / f"A_dRSA_alt_always_visible_corr_occipital{suffix}.png"
            ).write_bytes(b"fake png bytes")

    args = argparse.Namespace(
        group_subjects=["01", "02"],
        smoothing_plots=True,
        plot_smoothing_window=3,
    )
    report = STAGE07.build_group_average_report(cfg, args)
    figure_paths = [Path(item["figure"]) for item in report["outputs"]]

    assert report["status"] == "completed"
    assert len(figure_paths) == 2
    assert any(path.name.endswith("_groupAverage.png") for path in figure_paths)
    assert any(path.name.endswith("_groupAverage_smoothed.png") for path in figure_paths)
    assert all(path.is_file() for path in figure_paths)
    assert {item["condition"] for item in report["outputs"]} == {"always_visible"}
    raw_html = Path(
        tmp_path
        / "derivatives"
        / "reports"
        / "group"
        / "ses-01"
        / "07_matlab_drsa_group_report.html"
    )
    smoothed_html = raw_html.with_name("07_matlab_drsa_group_report_smoothed.html")
    raw_embedded_html = raw_html.with_name(
        "07_matlab_drsa_group_report_embedded.html"
    )
    smoothed_embedded_html = raw_html.with_name(
        "07_matlab_drsa_group_report_smoothed_embedded.html"
    )
    assert raw_html.is_file()
    assert smoothed_html.is_file()
    assert raw_embedded_html.is_file()
    assert smoothed_embedded_html.is_file()
    assert "sub-01" in raw_html.read_text(encoding="utf-8")
    assert "always_visible" in raw_html.read_text(encoding="utf-8")
    assert "group raw / non-smoothed | always_visible" in raw_html.read_text(
        encoding="utf-8"
    )
    assert "actual visible position" in raw_html.read_text(encoding="utf-8")
    assert "A_dRSA_alt_always_visible_corr_occipital.png" in raw_html.read_text(
        encoding="utf-8"
    )
    assert "data:image/png;base64," in raw_embedded_html.read_text(encoding="utf-8")
    assert "data:image/png;base64," in smoothed_embedded_html.read_text(
        encoding="utf-8"
    )
    assert report["report_files"]["raw_embedded_html"] == str(raw_embedded_html)
    assert report["report_files"]["smoothed_embedded_html"] == str(
        smoothed_embedded_html
    )
    assert "group smoothed display | always_visible" in smoothed_html.read_text(
        encoding="utf-8"
    )
    assert "A_dRSA_alt_always_visible_corr_occipital_smoothed.png" in smoothed_html.read_text(
        encoding="utf-8"
    )
    assert report["report_files"]["raw_html"] == str(raw_html)
    assert report["report_files"]["smoothed_html"] == str(smoothed_html)
    assert len(report["participant_figures"]) == 4


def test_variant_manifest_accepts_handoff_strings_and_preserves_source(tmp_path):
    """The completed-result provenance step must not mutate input ancestry."""
    source=tmp_path/'accepted'/'source_manifest.json';source.parent.mkdir()
    source.write_text(json.dumps({'derivative_lineage':{'stages':{'export_drsa':{'status':'completed'}}}}))
    before=source.read_bytes();target=tmp_path/'new'/'source_manifest.json'
    cfg={'drsa_input_resolved_paths':{},'drsa_analysis_label':'mnn-v1'}
    handoff={'source_manifest':str(source),'neural_normalization':{'mode':'mnn'}}
    result=STAGE07._completed_analysis_manifest(cfg,{'source_manifest':target},handoff,tmp_path/'raw',str(source))
    assert result==target and source.read_bytes()==before
    assert json.loads(target.read_text())['drsa_analysis']['normalization']['mode']=='mnn'
