"""
test_preprocessing_stage_policies.py

Purpose
Exercise pure policy decisions from Maxwell, artifact annotation, and manual
ICA without running expensive participant preprocessing.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_preprocessing_stage_policies.py -q
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import mne
import numpy as np
import pytest


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))


def _load_script(name: str, filename: str):
    """Import a numbered stage under a test-safe module name."""
    spec = importlib.util.spec_from_file_location(name, CODE_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MAXWELL = _load_script("dad_maxwell_stage", "03_A_MaxFilt.py")
ARTIFACT = _load_script("dad_artifact_stage", "03_B_ArtefactAnnotation.py")
ICA_STAGE = _load_script("dad_ica_stage", "03_C_ICA.py")


def test_reference_run_policy_reports_same_applied_list():
    """Legacy FLUX policy propagates the selected reference list explicitly."""
    selected, reference = MAXWELL.select_bad_channels_by_policy(
        {1: ["MEG0011"], 2: ["MEG0021"]},
        policy="reference_run",
        reference_acquisition_run=1,
        manual_cfg={"add": [], "remove": [], "by_run": {}},
    )
    assert reference == 1
    assert selected == {1: ["MEG0011"], 2: ["MEG0011"]}


def test_maxwell_auto_origin_reports_missing_headshape_as_version_issue():
    """Auto origin must not silently fall back when head-shape points vanish."""
    info = mne.create_info(["MEG0011"], sfreq=1000.0, ch_types=["mag"])
    with pytest.raises(RuntimeError, match=r"MNE >= 1\.10\.0"):
        MAXWELL.validate_auto_origin_digitization(
            info,
            origin="auto",
            source_filename="sub03_run01_200526.fif",
        )


def test_maxwell_numeric_origin_does_not_require_headshape_points():
    """An explicitly reviewed numerical origin remains usable without dig."""
    info = mne.create_info(["MEG0011"], sfreq=1000.0, ch_types=["mag"])
    summary = MAXWELL.validate_auto_origin_digitization(
        info,
        origin=(0.0, 0.0, 0.04),
        source_filename="synthetic_raw.fif",
    )
    assert summary["total"] == 0


def test_union_and_per_run_manual_overrides():
    """Alternative policies and manual adjustments remain deterministic."""
    union, _ = MAXWELL.select_bad_channels_by_policy(
        {1: ["MEG0011"], 2: ["MEG0021"]},
        policy="union",
        reference_acquisition_run=None,
        manual_cfg={
            "add": ["MEG0031"],
            "remove": ["MEG0011"],
            "by_run": {"run02": {"add": ["MEG0041"], "remove": ["MEG0021"]}},
        },
    )
    assert union[1] == ["MEG0021", "MEG0031"]
    assert union[2] == ["MEG0031", "MEG0041"]

    per_run, _ = MAXWELL.select_bad_channels_by_policy(
        {1: ["MEG0011"], 2: ["MEG0021"]},
        policy="per_run",
        reference_acquisition_run=None,
        manual_cfg={"add": [], "remove": [], "by_run": {}},
    )
    assert per_run == {1: ["MEG0011"], 2: ["MEG0021"]}


def test_central_destination_uses_actual_transform_medoid():
    """The selected destination must be an actual run, not an average matrix."""
    transforms = {}
    for run, translation_mm in ((1, 0.0), (2, 1.0), (3, 10.0)):
        matrix = np.eye(4)
        matrix[0, 3] = translation_mm / 1000.0
        transforms[run] = mne.transforms.Transform("meg", "head", matrix)
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.09, 0.0, 0.0],
            [0.0, 0.09, 0.0],
            [0.0, 0.0, 0.09],
        ]
    )

    central_run, details = MAXWELL.select_central_destination_run(
        transforms,
        head_points=points,
    )
    assert central_run == 2
    assert details["metric"].startswith("summed_pairwise_rms")
    assert details["score_by_run_mm"]["run02"] < details["score_by_run_mm"]["run03"]


def test_destination_configuration_requires_automatic_central_run():
    """A fixed destination must not bypass transform-medoid selection."""
    transforms = {
        1: mne.transforms.Transform("meg", "head", np.eye(4)),
        2: mne.transforms.Transform("meg", "head", np.eye(4)),
    }
    points = np.zeros((4, 3))
    with pytest.raises(ValueError, match="central_run"):
        MAXWELL.resolve_destination(
            destination_cfg={"policy": "central_run", "acquisition_run": 1},
            dev_head_t_by_run=transforms,
            head_points=points,
        )


def test_maxwell_resume_validates_lineage_bads_and_destination(tmp_path):
    """Resume accepts only the current bad list and destination transform."""
    info = mne.create_info(["MEG0011"], sfreq=1000.0, ch_types=["mag"])
    raw = mne.io.RawArray(np.zeros((1, 1000)), info, verbose="ERROR")
    destination = mne.transforms.Transform("meg", "head", np.eye(4))
    raw.info["dev_head_t"] = destination
    raw.info["bads"] = ["MEG0011"]
    fingerprint = MAXWELL.destination_transform_fingerprint(destination)
    raw.info["description"] = (
        "DAD source=sub01_run01_150426.fif; "
        "acquisition_run=1; stage=03_A_MaxFilt; destination_run=2; "
        f"destination_transform_sha256={fingerprint}"
    )
    output_path = tmp_path / "sub01_run01_150426_raw_sss.fif"
    raw.save(output_path, overwrite=True, verbose="ERROR")
    record = MAXWELL.AcquisitionRunRecord(
        source_path="/raw/sub01_run01_150426.fif",
        source_filename="sub01_run01_150426.fif",
        source_stem="sub01_run01_150426",
        acquisition_run_label=1,
        raw_relative_path="sub01/meg/sub01_run01_150426.fif",
    )

    result = MAXWELL._validate_resumable_output(
        path=output_path,
        record=record,
        applied_bads=["MEG0011"],
        destination_run=2,
        destination_transform=destination,
        destination_fingerprint=fingerprint,
        raw_root=tmp_path / "raw",
    )
    assert result["n_times"] == 1000

    with pytest.raises(RuntimeError, match="different bad channels"):
        MAXWELL._validate_resumable_output(
            path=output_path,
            record=record,
            applied_bads=[],
            destination_run=2,
            destination_transform=destination,
            destination_fingerprint=fingerprint,
            raw_root=tmp_path / "raw",
        )

    shifted_matrix = np.eye(4)
    shifted_matrix[1, 3] = 0.01
    shifted = mne.transforms.Transform("meg", "head", shifted_matrix)
    shifted_fingerprint = MAXWELL.destination_transform_fingerprint(shifted)
    with pytest.raises(RuntimeError, match="mismatched lineage"):
        MAXWELL._validate_resumable_output(
            path=output_path,
            record=record,
            applied_bads=["MEG0011"],
            destination_run=3,
            destination_transform=shifted,
            destination_fingerprint=shifted_fingerprint,
            raw_root=tmp_path / "raw",
        )


def test_eog_auto_skips_when_channels_are_absent():
    """Subjects 02-07 must never enter EOG-only processing."""
    info = mne.create_info(
        ["MEG0011", "MISC001"],
        sfreq=1000.0,
        ch_types=["mag", "misc"],
    )
    raw = mne.io.RawArray(np.zeros((2, 1000)), info, verbose="ERROR")
    mode, available, missing = ARTIFACT.resolve_eog_mode(
        "auto",
        ["EOG001", "EOG002"],
        raw,
    )
    assert mode == "disabled"
    assert available == []
    assert missing == ["EOG001", "EOG002"]


def test_eog_auto_is_diagnostics_only_when_present():
    """Sub01 EOG remains diagnostic rather than an automatic rejection source."""
    info = mne.create_info(
        ["MEG0011", "EOG001", "EOG002"],
        sfreq=1000.0,
        ch_types=["mag", "eog", "eog"],
    )
    raw = mne.io.RawArray(np.zeros((3, 1000)), info, verbose="ERROR")
    mode, available, missing = ARTIFACT.resolve_eog_mode(
        "auto",
        ["EOG001", "EOG002"],
        raw,
    )
    assert mode == "diagnostics_only"
    assert available == ["EOG001", "EOG002"]
    assert missing == []


def test_artifact_resume_lineage_uses_paths_from_common_run_report(tmp_path):
    """Resumed runs must not depend on new-output branch-local variables."""
    record = ARTIFACT.AcquisitionRunRecord(
        source_path="/raw/sub01_run01_block01_150426.fif",
        source_filename="sub01_run01_block01_150426.fif",
        source_stem="sub01_run01_block01_150426",
        acquisition_run_label=1,
        raw_relative_path="sub01_150426/meg/sub01_run01_block01_150426.fif",
    )
    output_fif = tmp_path / "sub01_run01_block01_150426_ann.fif"
    output_csv = tmp_path / "sub01_run01_block01_150426_ann.csv"
    lineage = ARTIFACT.build_annotation_lineage_record(
        record=record,
        input_path=tmp_path / "sub01_run01_block01_150426_raw_sss.fif",
        run_report={
            "output_fif": str(output_fif),
            "output_csv": str(output_csv),
            "execution_status": "resumed_existing",
            "muscle_annotation_count": 3,
            "eog_status": "resumed_diagnostics_only",
            "misc_channels_preserved": [],
        },
        accepted_threshold=10.0,
    )
    assert lineage["status"] == "resumed_existing"
    assert lineage["output_path"] == str(output_fif.resolve())
    assert lineage["metadata"]["annotation_csv"] == str(output_csv)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("", []),
        ("0, 5", [0, 5]),
        ("5 0 5", [0, 5]),
    ],
)
def test_manual_ica_component_parser(text, expected):
    """Manual exclusions are explicit, sorted, and de-duplicated."""
    assert ICA_STAGE.parse_component_indices(text) == expected


def test_manual_ica_component_parser_rejects_negative():
    """Negative component indices are never silently accepted."""
    with pytest.raises(ValueError, match="non-negative"):
        ICA_STAGE.parse_component_indices("-1")


def test_ica_cache_fingerprint_is_mount_independent(tmp_path):
    """Legacy local paths remain compatible with cluster-mounted inputs."""
    input_path = tmp_path / "sub02_run01_130526_ann.fif"
    import numpy as np
    raw = mne.io.RawArray(np.zeros((2, 500)), mne.create_info(["MEG0111", "MEG0121"], 100., ["mag", "mag"]), verbose="ERROR")
    raw.save(input_path, overwrite=True, verbose="ERROR")
    ica_cfg = {
        "method": "fastica",
        "n_components": 30,
        "random_state": 97,
        "max_iter": "auto",
        "fit_resample_sfreq_hz": 200.0,
        "fit_l_freq_hz": 1.0,
        "fit_h_freq_hz": 40.0,
    }
    expected = ICA_STAGE.build_ica_fingerprint(
        input_paths=[input_path],
        ica_cfg=ica_cfg,
    )
    legacy = {
        "inputs": [
            {
                "path": (
                    "/Volumes/MORWUR/elsewhere/"
                    "sub02_run01_130526_ann.fif"
                ),
                "size_bytes": input_path.stat().st_size,
                "mtime_ns": input_path.stat().st_mtime_ns,
                "annotations": expected["inputs"][0]["annotations"],
            }
        ],
        "fit": expected["fit"],
        "sha256": "legacy-absolute-path-hash",
    }
    metadata_path = tmp_path / "ica_metadata.json"
    metadata_path.write_text(
        json.dumps({"fingerprint": legacy}),
        encoding="utf-8",
    )
    assert ICA_STAGE.cache_is_compatible(metadata_path, expected)


def test_auto_ica_is_deprecated_before_optional_dependency_import():
    """Production execution stops with a clear status, not a MEGNet import error."""
    completed = subprocess.run(
        [sys.executable, str(CODE_DIR / "03_C_ICA_auto.py")],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 2
    assert "deprecated experimental code" in completed.stderr
    assert "mne-icalabel" not in completed.stderr


def test_central_run_reports_ties_and_rejects_missing_or_invalid_transforms():
    """Retained medoid selection is deterministic and fails on unusable headers."""
    import numpy as np
    transform = mne.transforms.Transform("meg", "head", np.eye(4))
    points = np.vstack([np.eye(3), np.zeros(3)])
    run, details = MAXWELL.select_central_destination_run({2: transform, 1: transform}, head_points=points)
    assert run == 1 and details["tie_candidates"] == [1, 2]
    with pytest.raises(RuntimeError, match="no static"):
        MAXWELL.select_central_destination_run({1: None}, head_points=points)
    invalid = {"trans": np.full((4, 4), np.nan)}
    with pytest.raises(RuntimeError, match="invalid"):
        MAXWELL.select_central_destination_run({1: invalid}, head_points=points)
