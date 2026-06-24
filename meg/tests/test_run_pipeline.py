"""
test_run_pipeline.py

Purpose
Validate orchestration command contracts for stable preprocessing and
downstream scientific stages.

Coverage
- the numbered trigger entry point is registered,
- the durable source manifest is passed to the stage,
- optional inventory reports are not scientific stage inputs,
- the deep-report opt-out flag is forwarded only when explicitly requested,
- production 03_A/03_B/manual-ICA stages receive review/apply controls,
- epoching/export receive the manifest without Stage-03 interactive controls,
- the deprecated automatic ICA branch is not registered.

Usage
Run only these tests:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_run_pipeline.py -q
"""

from __future__ import annotations

from pathlib import Path
import sys


CODE_DIR = (
    Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
)
sys.path.insert(0, str(CODE_DIR))

from run_pipeline import (  # noqa: E402
    build_stage_commands,
    preserve_durable_manifest_sections,
)


def test_trigger_command_receives_manifest_but_not_inventory(tmp_path):
    """The durable manifest is provenance; inventory JSON is never an input."""
    config_path = tmp_path / "pipeline_default.yaml"
    subject_config_path = tmp_path / "sub-04.yaml"
    manifest_path = tmp_path / "source_manifest.json"

    command = build_stage_commands(
        config_path=config_path,
        subject_config_path=subject_config_path,
        source_manifest_path=manifest_path,
        disable_deep_trigger_report=False,
    )["trigger_extraction"]

    assert command[1].endswith("02_trigger_extraction_custom.py")
    assert command[command.index("--source-manifest") + 1] == str(manifest_path)
    assert command[command.index("--subject-config") + 1] == str(
        subject_config_path
    )
    assert not any("inventory" in argument for argument in command)


def test_site_config_is_forwarded_to_every_stage(tmp_path):
    """All subprocesses must resolve the same cluster path override."""
    site_path = tmp_path / "pipeline_cluster.yaml"
    commands = build_stage_commands(
        config_path=tmp_path / "pipeline_default.yaml",
        site_config_path=site_path,
        subject_config_path=tmp_path / "sub-02.yaml",
        source_manifest_path=tmp_path / "source_manifest.json",
        disable_deep_trigger_report=False,
    )
    for command in commands.values():
        assert command[command.index("--site-config") + 1] == str(site_path)


def test_deep_report_opt_out_is_forwarded_only_when_requested(tmp_path):
    """Deep reporting remains enabled in the default stage command."""
    common = {
        "config_path": tmp_path / "pipeline_default.yaml",
        "subject_config_path": None,
        "source_manifest_path": tmp_path / "source_manifest.json",
    }
    default_command = build_stage_commands(
        **common,
        disable_deep_trigger_report=False,
    )["trigger_extraction"]
    compatibility_command = build_stage_commands(
        **common,
        disable_deep_trigger_report=True,
    )["trigger_extraction"]

    assert "--disable-deep-trigger-report" not in default_command
    assert compatibility_command[-1] == "--disable-deep-trigger-report"


def test_production_preprocessing_commands_receive_review_controls(tmp_path):
    """Stable stages share tuning, pass, interactivity, and overwrite flags."""
    commands = build_stage_commands(
        config_path=tmp_path / "pipeline_default.yaml",
        subject_config_path=tmp_path / "sub-04.yaml",
        source_manifest_path=tmp_path / "source_manifest.json",
        disable_deep_trigger_report=False,
        tuning_config_path=tmp_path / "03_tuning_choices.yaml",
        processing_pass="apply",
        interactive=True,
        no_display=True,
        decisions_out=tmp_path / "accepted.yaml",
        overwrite=True,
    )

    assert set(commands) == {
        "trigger_extraction",
        "maxwell_filter",
        "artifact_annotation",
        "ica_manual",
        "epoching",
        "export_drsa",
        "matlab_drsa",
    }
    assert "ica_auto" not in commands
    for stage in ("maxwell_filter", "artifact_annotation", "ica_manual"):
        command = commands[stage]
        assert command[command.index("--processing-pass") + 1] == "apply"
        assert "--interactive" in command
        assert "--no-display" in command
        assert "--overwrite" in command
        assert command[command.index("--tuning-config") + 1].endswith(
            "03_tuning_choices.yaml"
        )


def test_downstream_commands_receive_lineage_and_noninteractive_controls(tmp_path):
    """Downstream stages use the durable manifest and never inherit review flags."""
    manifest_path = tmp_path / "source_manifest.json"
    commands = build_stage_commands(
        config_path=tmp_path / "pipeline_default.yaml",
        site_config_path=tmp_path / "pipeline_cluster.yaml",
        subject_config_path=tmp_path / "sub-04.yaml",
        source_manifest_path=manifest_path,
        disable_deep_trigger_report=False,
        tuning_config_path=tmp_path / "03_tuning_choices.yaml",
        processing_pass="apply",
        interactive=True,
        overwrite=True,
        resume=False,
    )
    assert commands["epoching"][1].endswith("04_epoching.py")
    assert commands["export_drsa"][1].endswith("06_export_neural_for_drsa.py")
    assert commands["matlab_drsa"][1].endswith("07_run_matlab_drsa.py")
    for stage in ("epoching", "export_drsa", "matlab_drsa"):
        command = commands[stage]
        assert command[command.index("--source-manifest") + 1] == str(
            manifest_path
        )
        assert "--processing-pass" not in command
        assert "--interactive" not in command
        assert "--overwrite" in command


def test_resume_is_forwarded_to_supported_stages(tmp_path):
    """Resume is forwarded only to stages that validate existing outputs."""
    commands = build_stage_commands(
        config_path=tmp_path / "pipeline_default.yaml",
        subject_config_path=tmp_path / "sub-01.yaml",
        source_manifest_path=tmp_path / "source_manifest.json",
        disable_deep_trigger_report=False,
        resume=True,
    )
    assert "--resume" in commands["maxwell_filter"]
    assert "--resume" in commands["artifact_annotation"]
    assert "--resume" in commands["ica_manual"]
    assert "--resume" in commands["epoching"]
    assert "--resume" in commands["export_drsa"]
    assert "--resume" in commands["matlab_drsa"]


def test_manifest_refresh_preserves_evidence_derived_sections(tmp_path):
    """Dry runs must not erase trigger associations or derivative lineage."""
    manifest_path = tmp_path / "source_manifest.json"
    manifest_path.write_text(
        """
        {
          "trigger_derived_associations": {"status": "resolved"},
          "derivative_lineage": {"stages": {"maxwell_filter": {"status": "completed"}}}
        }
        """,
        encoding="utf-8",
    )
    refreshed = preserve_durable_manifest_sections(
        {"sources": {"task_meg_acquisition_runs": []}},
        manifest_path,
    )
    assert refreshed["trigger_derived_associations"]["status"] == "resolved"
    assert (
        refreshed["derivative_lineage"]["stages"]["maxwell_filter"]["status"]
        == "completed"
    )
