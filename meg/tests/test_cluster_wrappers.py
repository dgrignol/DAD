"""
test_cluster_wrappers.py

Purpose
Validate the lightweight Sun Grid Engine submission wrappers without
submitting jobs or running scientific preprocessing. In particular, verify
that the apply wrapper remains valid when Grid Engine executes a temporary
spooled copy rather than the repository copy.

Usage
Run only these tests:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_cluster_wrappers.py -q
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]
APPLY_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_preprocessing_apply.sh"
DRSA_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_drsa_analysis.sh"
DRSA_LINE_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_drsa_line_summaries.sh"
DOT_SIDE_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_dot_side_decoding.sh"
DOT_SIDE_GROUP_WRAPPER = REPO_ROOT / "meg" / "cluster" / "submit_dot_side_group_report.sh"
DEPLOY_WRAPPER = REPO_ROOT / "meg" / "cluster" / "deploy_review_to_cluster.sh"


def test_apply_wrapper_resolves_review_wrapper_from_project_root(tmp_path):
    """
    A spooled apply script must delegate through DAD_PROJECT_ROOT.

    Grid Engine changes BASH_SOURCE to its temporary job-script path. This
    test copies the apply wrapper into an unrelated spool directory and uses a
    harmless review-wrapper stub to record the delegated environment.
    """
    project_root = tmp_path / "project"
    cluster_dir = project_root / "meg" / "cluster"
    cluster_dir.mkdir(parents=True)
    output_path = tmp_path / "delegated.txt"

    review_wrapper = cluster_dir / "submit_preprocessing_review.sh"
    review_wrapper.write_text(
        "#!/bin/bash\n"
        'printf "%s|%s\\n" "$DAD_PROCESSING_PASS" "$DAD_REBUILD" '
        '>"$DAD_TEST_OUTPUT"\n',
        encoding="utf-8",
    )
    review_wrapper.chmod(0o755)

    spool_dir = tmp_path / "gridengine" / "job_scripts"
    spool_dir.mkdir(parents=True)
    spooled_apply = spool_dir / "91761"
    shutil.copy2(APPLY_WRAPPER, spooled_apply)
    spooled_apply.chmod(0o755)

    env = os.environ.copy()
    env.update(
        {
            "DAD_PROJECT_ROOT": str(project_root),
            "DAD_TEST_OUTPUT": str(output_path),
        }
    )
    completed = subprocess.run(
        [str(spooled_apply)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert completed.returncode == 0, completed.stderr
    assert output_path.read_text(encoding="utf-8") == "apply|1\n"


def test_drsa_and_deploy_wrappers_have_valid_bash_syntax():
    """Cluster submission/deployment shell scripts must parse before upload."""
    for script in (
        DRSA_WRAPPER,
        DRSA_LINE_WRAPPER,
        DOT_SIDE_WRAPPER,
        DOT_SIDE_GROUP_WRAPPER,
        DEPLOY_WRAPPER,
    ):
        completed = subprocess.run(
            ["bash", "-n", str(script)],
            check=False,
            capture_output=True,
            text=True,
        )
        assert completed.returncode == 0, completed.stderr


def test_deploy_payload_includes_matlab_analysis_layer():
    """Stage-07 deployment must include MATLAB scripts and helper functions."""
    text = DEPLOY_WRAPPER.read_text(encoding="utf-8")
    assert '"meg/code/analysis/"' in text
    assert '"meg/code/preprocessing_basic/07_run_matlab_drsa.py"' in text
    assert '"meg/code/preprocessing_basic/08_plot_drsa_line_summaries.py"' in text
    assert '"meg/code/preprocessing_basic/drsa_result_utils.py"' in text
    assert '"meg/cluster/submit_drsa_analysis.sh"' in text
    assert '"meg/cluster/submit_drsa_line_summaries.sh"' in text
    assert '"simulations/functions/"' in text


def test_deploy_payload_includes_dot_side_sanity_check():
    """Optional dot-side cluster jobs need their script and wrapper deployed."""
    text = DEPLOY_WRAPPER.read_text(encoding="utf-8")
    assert '"meg/code/sanity_checks/"' in text
    assert '"meg/cluster/submit_dot_side_decoding.sh"' in text
    assert '"meg/cluster/submit_dot_side_group_report.sh"' in text


def test_dot_side_wrapper_exposes_decoding_controls():
    """The dot-side wrapper should expose the requested runtime knobs."""
    text = DOT_SIDE_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DOT_SIDE_ANALYSES" in text
    assert "DAD_DOT_SIDE_CONDITION" in text
    assert "DAD_DOT_SIDE_MODEL_SOURCE" in text
    assert "occluded_deviant_predicted" in text
    assert "DAD_DECODE_WINDOW_SAMPLES" in text
    assert "DAD_TRIALS_PER_AVERAGE" in text
    assert "DAD_RANDOM_SEED" in text
    assert "--validate-only" in text


def test_dot_side_group_wrapper_exposes_cluster_controls():
    """The group wrapper should expose CI/SEM and cluster-threshold settings."""
    text = DOT_SIDE_GROUP_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DOT_SIDE_ANALYSES" in text
    assert "--analyses" in text
    assert "DAD_DOT_SIDE_INTERVAL" in text
    assert "DAD_CLUSTER_FORMING_P" in text
    assert "DAD_CLUSTER_ALPHA" in text
    assert "DAD_DOT_SIDE_GROUP_SUBJECTS" in text
    assert "submit_dot_side_group_report.sh" in text


def test_drsa_line_summary_wrapper_exposes_stage08_controls():
    """The Stage-08 wrapper should expose group line-summary controls."""
    text = DRSA_LINE_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DRSA_LINE_SUBJECTS" in text
    assert "DAD_DRSA_LINE_LAG_WINDOW_SEC" in text
    assert "DAD_DRSA_LINE_INCLUDE_SMOOTHED" in text
    assert "DAD_DRSA_LINE_SHADING" in text
    assert "DAD_DRSA_LINE_INTERVAL_LAG_MODE" in text
    assert "DAD_DRSA_LINE_MODEL_TIME_WINDOW_SAMPLES" not in text
    assert "--model-time-window-samples" not in text
    assert "--shading" in text
    assert "--interval-lag-mode" in text
    assert "08_plot_drsa_line_summaries.py" in text


def test_drsa_wrapper_exposes_model_source_controls():
    """Stage-07 cluster submissions should expose observed/predicted model runs."""
    text = DRSA_WRAPPER.read_text(encoding="utf-8")
    assert "DAD_DRSA_MODEL_SOURCES" in text
    assert "DAD_DRSA_CONDITIONS" in text


def test_drsa_wrapper_uses_validated_cimec_matlab_resources():
    """Stage 07 should discover the compute-node MATLAB and retain JVM memory."""
    text = DRSA_WRAPPER.read_text(encoding="utf-8")
    assert "/state/partition1/MATLAB/R2020a/bin/matlab" in text
    assert "#$ -l h_vmem=24G" in text
