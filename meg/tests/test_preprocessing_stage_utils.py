"""
test_preprocessing_stage_utils.py

Purpose
Validate shared stage safety, acquisition-run lineage, and tuning persistence.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_preprocessing_stage_utils.py -q
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest
import yaml


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))

from preprocessing_stage_utils import (  # noqa: E402
    AcquisitionRunRecord,
    derivative_path_for_run,
    parse_acquisition_run_label,
    update_source_manifest_stage,
    update_tuning_file,
    write_json,
)


def test_run_parser_ignores_behavioral_block_tokens():
    """Only the acquisition run token determines run identity."""
    name = "sub01_run02_block02_03_150426.fif"
    assert parse_acquisition_run_label(name) == 2


def test_derivative_name_preserves_complete_source_stem(tmp_path):
    """Irregular sub01 source identity survives every derivative name."""
    record = AcquisitionRunRecord(
        source_path="/raw/sub01_run08_block11_150426.fif",
        source_filename="sub01_run08_block11_150426.fif",
        source_stem="sub01_run08_block11_150426",
        acquisition_run_label=8,
        raw_relative_path="sub01/meg/sub01_run08_block11_150426.fif",
    )
    assert derivative_path_for_run(record, tmp_path, "_ica.fif").name == (
        "sub01_run08_block11_150426_ica.fif"
    )


def test_tuning_file_merges_stage_decisions(tmp_path):
    """Interactive stages update one reusable subject tuning YAML."""
    raw_root = tmp_path / "raw"
    tuning_path = tmp_path / "derivatives" / "03_tuning_choices.yaml"
    update_tuning_file(
        path=tuning_path,
        raw_root=raw_root,
        subject_id="04",
        section="artifact_annotation",
        values={"muscle_threshold_by_run": {"run01": 8.0}},
    )
    update_tuning_file(
        path=tuning_path,
        raw_root=raw_root,
        subject_id="04",
        section="ica",
        values={"manual_exclude": [1, 4]},
    )
    payload = yaml.safe_load(tuning_path.read_text(encoding="utf-8"))
    assert payload["subject"]["id"] == "04"
    assert payload["artifact_annotation"]["muscle_threshold_by_run"] == {
        "run01": 8.0
    }
    assert payload["ica"]["manual_exclude"] == [1, 4]


def test_tuning_file_rejects_raw_tree_output(tmp_path):
    """Accepted decisions are outputs and must never enter the raw tree."""
    raw_root = tmp_path / "raw"
    with pytest.raises(ValueError, match="read-only raw-data tree"):
        update_tuning_file(
            path=raw_root / "sub01" / "choices.yaml",
            raw_root=raw_root,
            subject_id="01",
            section="ica",
            values={"manual_exclude": []},
        )


def test_json_writer_normalizes_numpy_values(tmp_path):
    """Scientific report payloads serialize NumPy scalars and arrays."""
    output = tmp_path / "derivatives" / "report.json"
    write_json(
        output,
        {
            "count": np.int64(5),
            "score": np.float64(0.25),
            "indices": np.array([1, 3], dtype=np.int64),
        },
        raw_root=tmp_path / "raw",
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload == {"count": 5, "indices": [1, 3], "score": 0.25}


def test_manifest_stage_update_preserves_existing_sections(tmp_path):
    """Derivative lineage augments rather than replaces source provenance."""
    raw_root = tmp_path / "raw"
    manifest_path = tmp_path / "derivatives" / "source_manifest.json"
    manifest_path.parent.mkdir(parents=True)
    manifest_path.write_text(
        json.dumps(
            {
                "raw_dataset": {"data_root": str(raw_root)},
                "sources": {"task_meg_acquisition_runs": [{"name": "run01"}]},
            }
        ),
        encoding="utf-8",
    )
    status = update_source_manifest_stage(
        manifest_path=manifest_path,
        raw_root=raw_root,
        stage_name="maxwell_filter",
        processing_pass="review",
        records=[{"acquisition_run_label": 1}],
        status="completed",
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert status == "updated"
    assert payload["sources"]["task_meg_acquisition_runs"] == [{"name": "run01"}]
    assert (
        payload["derivative_lineage"]["stages"]["maxwell_filter"]["records"][0][
            "acquisition_run_label"
        ]
        == 1
    )
