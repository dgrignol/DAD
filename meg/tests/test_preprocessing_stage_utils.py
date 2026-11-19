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
import subprocess
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


def test_tuning_file_normalizes_nested_ica_metadata_and_preserves_other_decisions(tmp_path):
    """A real NumPy ICA-duration value must not break the final ocular SAVE."""
    tuning = tmp_path / "tuning.yaml"
    tuning.write_text("stage4_review:\n  03_A_MaxFilt:\n    status: accepted\noriginal_integer_keys:\n  1: retained\n")
    update_tuning_file(path=tuning, raw_root=tmp_path / "raw", subject_id="04",
        section="artifact_annotation", values={"manual_eye_reviews": {"source.fif": {
            "preliminary_ica": {"fit_duration_sec": np.float64(800.79),
                "n_components_fitted": np.int64(12), "selected_components": np.array([0, 2])}}}})
    saved = yaml.safe_load(tuning.read_text())
    assert saved["stage4_review"]["03_A_MaxFilt"]["status"] == "accepted"
    assert saved["original_integer_keys"] == {1: "retained"}
    assert saved["artifact_annotation"]["manual_eye_reviews"]["source.fif"]["preliminary_ica"] == {
        "fit_duration_sec": 800.79, "n_components_fitted": 12, "selected_components": [0, 2]}


def test_source_receipts_merge_latest_file_and_replace_one_complete_entry(tmp_path):
    """Saving a stale GUI must retain other runs and remove old same-run fields."""
    path = tmp_path / "tuning.yaml"
    path.write_text(yaml.safe_dump({"stage4_review": {"03_A_MaxFilt": {"status": "accepted"}},
        "artifact_annotation": {"manual_eye_reviews": {
            "run01.fif": {"obsolete": True}, "run02.fif": {"status": "reviewed"}}}}))
    update_tuning_file(path=path, raw_root=tmp_path / "raw", subject_id="04",
        section="artifact_annotation", mapping_key="manual_eye_reviews",
        values={"run01.fif": {"status": "reviewed", "annotation_count": np.int64(0)}})
    saved = yaml.safe_load(path.read_text())
    assert saved["artifact_annotation"]["manual_eye_reviews"] == {
        "run01.fif": {"status": "reviewed", "annotation_count": 0},
        "run02.fif": {"status": "reviewed"}}
    assert saved["stage4_review"]["03_A_MaxFilt"]["status"] == "accepted"


def test_simultaneous_tuning_writers_preserve_all_source_receipts(tmp_path):
    """Real independent processes must serialize read/merge/write, not lose runs."""
    path = tmp_path / "tuning.yaml"
    script = """
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from preprocessing_stage_utils import update_tuning_file
for index in range(5):
    update_tuning_file(path=sys.argv[2], raw_root=Path(sys.argv[2]).parent / 'raw',
        subject_id='04', section='artifact_annotation', mapping_key='manual_eye_reviews',
        values={f'{sys.argv[3]}_{index}.fif': {'status': 'reviewed'}})
"""
    processes = [subprocess.Popen([sys.executable, "-c", script, str(CODE_DIR), str(path), str(worker)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) for worker in range(4)]
    try:
        for process in processes:
            stdout, stderr = process.communicate(timeout=30)
            assert process.returncode == 0, stdout + stderr
    finally:
        for process in processes:
            if process.poll() is None:
                process.kill()
                process.wait()
    saved = yaml.safe_load(path.read_text())["artifact_annotation"]["manual_eye_reviews"]
    assert set(saved) == {f"{worker}_{index}.fif" for worker in range(4) for index in range(5)}


def test_failed_atomic_tuning_replacement_keeps_previous_file(tmp_path, monkeypatch):
    """A failed save must leave readable previous acceptance and no temp debris."""
    import preprocessing_stage_utils as utils
    path = tmp_path / "tuning.yaml"
    path.write_text("subject:\n  id: '04'\nstage4_review:\n  retained: true\n")
    before = path.read_bytes()

    def fail_replace(*args):
        raise OSError("Synthetic replacement failure")

    monkeypatch.setattr(utils.os, "replace", fail_replace)
    with pytest.raises(OSError, match="replacement failure"):
        update_tuning_file(path=path, raw_root=tmp_path / "raw", subject_id="04",
            section="ica", values={"manual_exclude": []})
    assert path.read_bytes() == before
    assert not list(tmp_path.glob("*.tmp"))


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
