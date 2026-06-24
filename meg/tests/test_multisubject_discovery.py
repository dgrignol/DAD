"""
test_multisubject_discovery.py

Purpose
Exercise the external multisubject discovery and source-lineage contracts
without reading or writing the mounted raw dataset.

Usage examples
1) Run this focused suite:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_multisubject_discovery.py -q

2) Run the complete MEG test suite:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests -q

Coverage
- exact subject-folder matching for participants 01-07,
- zero and ambiguous folder failures,
- profile-derived trigger and identity settings,
- sub01 irregular acquisition filenames,
- sub07 date-label handling and empty-room separation,
- missing required scientific inputs,
- manifest preservation without positional MEG/block association.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pytest
import yaml


MEG_ROOT = Path(__file__).resolve().parents[1]
PREPROCESSING_DIR = MEG_ROOT / "code" / "preprocessing_basic"
DEFAULT_CONFIG = MEG_ROOT / "params" / "pipeline_default.yaml"
sys.path.insert(0, str(PREPROCESSING_DIR))

import dataset_inventory  # noqa: E402
import pipeline_config  # noqa: E402


SUBJECT_DATES = {
    "01": "150426",
    "02": "130526",
    "03": "200526",
    "04": "220526",
    "05": "260526",
    "06": "260526",
    "07": "270526",
}


def _load_subject_config(tmp_path: Path, subject_id: str):
    """
    Load the real profile/config schema against a synthetic external tree.

    Only project roots and the subject selector are overridden. This keeps the
    tests coupled to production profile derivation while avoiding the mounted
    read-only volume.
    """
    data_root = tmp_path / "raw"
    subject_folder = data_root / f"sub{subject_id}_{SUBJECT_DATES[subject_id]}"
    subject_folder.mkdir(parents=True, exist_ok=True)
    override_path = tmp_path / f"sub-{subject_id}.yaml"
    override_path.write_text(
        yaml.safe_dump(
            {
                "project": {
                    "data_root": str(data_root),
                    "derivatives_root": str(tmp_path / "derivatives"),
                },
                "subject": {
                    "id": subject_id,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    cfg, warnings, default_path, subject_path = (
        pipeline_config.load_and_resolve_config(
            DEFAULT_CONFIG,
            override_path,
        )
    )
    return cfg, warnings, default_path, subject_path


def _create_source_file(subject_root: Path, relative_path: str) -> Path:
    """Create one small placeholder source file in a synthetic raw tree."""
    path = subject_root / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"source")
    return path


@pytest.mark.parametrize("subject_id", pipeline_config.SUPPORTED_SUBJECT_IDS)
def test_all_subject_folder_matches_are_exact(tmp_path, subject_id):
    """Each supported selector resolves its unique external subNN_* folder."""
    expected = tmp_path / f"sub{subject_id}_{SUBJECT_DATES[subject_id]}"
    expected.mkdir()

    discovered = pipeline_config.discover_external_subject_folder(
        tmp_path,
        subject_id,
    )

    assert discovered == expected.resolve()


def test_subject_folder_discovery_fails_on_zero_matches(tmp_path):
    """A missing participant folder must fail before downstream discovery."""
    with pytest.raises(FileNotFoundError, match="No external subject folder"):
        pipeline_config.discover_external_subject_folder(tmp_path, "03")


def test_subject_folder_discovery_fails_on_ambiguous_matches(tmp_path):
    """Two subNN_* directories are unsafe and must not be guessed."""
    (tmp_path / "sub03_200526").mkdir()
    (tmp_path / "sub03_duplicate").mkdir()

    with pytest.raises(ValueError, match="Ambiguous external subject folders"):
        pipeline_config.discover_external_subject_folder(tmp_path, "03")


@pytest.mark.parametrize(
    ("subject_id", "condition_codes", "identity_max"),
    [
        ("01", [31, 41, 51], 20),
        ("04", [71, 81, 91], 40),
        ("07", [71, 81, 91], 40),
    ],
)
def test_stage_configuration_is_derived_from_subject_profile(
    tmp_path,
    subject_id,
    condition_codes,
    identity_max,
):
    """Condition and identity settings come from the selected profile."""
    cfg, _warnings, _default_path, _subject_path = _load_subject_config(
        tmp_path,
        subject_id,
    )

    assert list(cfg["epoching"]["event_id_map"].values()) == condition_codes
    assert cfg["trigger_parser"]["expected_sequence_code_min"] == 1
    assert cfg["trigger_parser"]["expected_sequence_code_max"] == identity_max
    assert cfg["trigger_parser"]["required_trigger_codes"][:3] == condition_codes
    runtime = cfg["subject_profile"]["versions"]["runtime"]
    assert cfg["trigger_parser"]["runtime_version"] == runtime
    assert f"occlusion_{runtime}_rescueTraject" in cfg["data_discovery"][
        "reported_trigger_csv_glob"
    ]
    assert cfg["trigger_parser"]["gaze_break_code"] == 150
    assert cfg["trigger_parser"]["replay_start_code"] == 151
    assert cfg["trigger_parser"]["esc_termination_code"] == 152


def test_sub01_irregular_run_names_are_preserved(tmp_path):
    """Sub01 discovery keeps all seven acquisition filenames, including run08."""
    cfg, _warnings, _default_path, _subject_path = _load_subject_config(
        tmp_path,
        "01",
    )
    subject_root = Path(cfg["resolved_paths"]["subject_data_dir"])
    names = [
        "sub01_run01_block01_150426.fif",
        "sub01_run02_block02_03_150426.fif",
        "sub01_run03_block04_05_150426.fif",
        "sub01_run04_block06_07_150426.fif",
        "sub01_run05_block08_09_150426.fif",
        "sub01_run06_block10_150426.fif",
        "sub01_run08_block11_150426.fif",
    ]
    for name in names:
        _create_source_file(subject_root, f"meg/{name}")

    sources = dataset_inventory.discover_subject_sources(cfg)

    assert [path.name for path in sources["task_meg_fif"]] == names
    assert sources["empty_room_fif"] == []


def test_sub07_date_mismatch_uses_raw_folder_date(tmp_path):
    """The raw 270526 label wins over the prior 260527 migration-note label."""
    cfg, warnings, _default_path, _subject_path = _load_subject_config(
        tmp_path,
        "07",
    )

    assert warnings == []
    assert cfg["resolved_paths"]["subject_folder_name"] == "sub07_270526"
    assert cfg["resolved_paths"]["subject_folder_date"] == "270526"
    assert (
        cfg["subject_profile"]["external_folder"]["prior_migration_label"]
        == "sub07_260527"
    )
    assert cfg["data_discovery"]["empty_room_fif_glob"] == "empty_room_270526.fif"


def test_sub07_empty_room_is_excluded_from_task_runs(tmp_path):
    """The root-level empty-room FIF is discovered separately from task FIFs."""
    cfg, _warnings, _default_path, _subject_path = _load_subject_config(
        tmp_path,
        "07",
    )
    subject_root = Path(cfg["resolved_paths"]["subject_data_dir"])
    data_root = Path(cfg["resolved_paths"]["data_root"])
    task_fif = _create_source_file(subject_root, "meg/sub07_run01_270526.fif")
    empty_room = data_root / "empty_room_270526.fif"
    empty_room.write_bytes(b"empty-room")

    sources = dataset_inventory.discover_subject_sources(cfg)

    assert sources["task_meg_fif"] == [task_fif.resolve()]
    assert sources["empty_room_fif"] == [empty_room.resolve()]
    assert empty_room.resolve() not in sources["task_meg_fif"]


def test_inventory_reports_missing_required_inputs(tmp_path, monkeypatch):
    """Missing trajectory and TrialStruct files produce explicit fatal errors."""
    cfg, warnings, default_path, subject_path = _load_subject_config(
        tmp_path,
        "07",
    )
    subject_root = Path(cfg["resolved_paths"]["subject_data_dir"])
    _create_source_file(subject_root, "meg/sub07_run01_270526.fif")
    _create_source_file(
        subject_root,
        "output_files/MoveDot1_occlusion_v20_rescueTraject_SUB07_BLOCK01.mat",
    )
    _create_source_file(
        subject_root,
        "output_files/debug_actual_triggers_occlusion_v20_rescueTraject_"
        "sub07_block01_run01.csv",
    )
    _create_source_file(subject_root, "eyetracker/MB070100.edf")

    # This test targets source validation, so avoid parsing the placeholder FIF.
    monkeypatch.setattr(
        dataset_inventory,
        "summarize_fif_file",
        lambda path, stim_channel: {
            "path": str(path),
            "name": Path(path).name,
            "size_bytes": Path(path).stat().st_size,
            "sfreq_hz": 1000.0,
            "n_channels": 1,
            "n_samples": 1,
            "duration_sec": 0.001,
            "stim_channel_expected": stim_channel,
            "stim_channel_present": True,
            "channel_type_counts": {"stim": 1},
            "eog_channel_names": [],
            "misc_channel_names": [],
        },
    )

    payload = dataset_inventory.build_inventory_payload(
        cfg,
        default_path,
        subject_path,
        warnings,
    )

    assert payload["validation"]["is_valid"] is False
    assert any("presented trajectory MAT" in error for error in payload["validation"]["errors"])
    assert any("TrialStruct MAT" in error for error in payload["validation"]["errors"])


def test_source_manifest_keeps_modalities_separate(tmp_path):
    """Manifest lineage preserves filenames and leaves cross-list mapping empty."""
    cfg, warnings, default_path, subject_path = _load_subject_config(
        tmp_path,
        "01",
    )
    subject_root = Path(cfg["resolved_paths"]["subject_data_dir"])
    _create_source_file(
        subject_root,
        "meg/sub01_run02_block02_03_150426.fif",
    )
    _create_source_file(
        subject_root,
        "input_files/MovDot_Sub01_V28_rescueTraject.mat",
    )
    _create_source_file(
        subject_root,
        "input_files/Sub01_TrialStruct_v21_rescueTraject.mat",
    )
    _create_source_file(
        subject_root,
        "output_files/MoveDot1_occlusion_v18_rescueTraject_SUB01_BLOCK02.mat",
    )
    _create_source_file(
        subject_root,
        "output_files/debug_actual_triggers_occlusion_v18_rescueTraject_"
        "sub01_block02_run01.csv",
    )

    manifest = dataset_inventory.build_source_manifest_payload(
        cfg,
        default_path,
        subject_path,
        warnings,
    )

    task_record = manifest["sources"]["task_meg_acquisition_runs"][0]
    assert task_record["source_filename"] == "sub01_run02_block02_03_150426.fif"
    assert task_record["acquisition_run_label"] == 2
    assert "behavioral_block_label" not in task_record
    assert manifest["trigger_derived_associations"]["status"] == "pending"
    assert (
        manifest["trigger_derived_associations"][
            "meg_acquisition_run_to_behavioral_segments"
        ]
        == []
    )
