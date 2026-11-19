"""Synthetic regression checks for participant books and bulk note import.

Usage
-----
    python -B -m pytest -p no:cacheprovider meg/tests/test_ocular_component_reports.py -q

Only temporary fixtures are used. Check exact native peak coordinates, display
gap preservation, offline one-book/multiple-run structure, resumed generation,
manifest/HTML tamper detection, unfilled versus NONE/DEFER decisions, fitted
component bounds, current-model identity, whole-table preflight and preservation
of existing choices. Scientific detector samples and source files stay unchanged.
"""

from copy import deepcopy
import argparse
import base64
import csv
import json
from pathlib import Path
import sys
import shutil
import subprocess
from types import SimpleNamespace

import mne
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import annotation_triage as triage
import build_ocular_component_reports as reports
import import_ocular_component_choices as importer
from continuous_cleaning import digest
from preprocessing_stage_utils import AcquisitionRunRecord


# -----------------------------------------------------------------------------
# Deterministic peak/gap fixtures and lightweight book-generation dependencies
# -----------------------------------------------------------------------------

def example(tmp_path):
    """Make two native component traces, two independent runs and a stable model identity."""
    cfg = {"subject": {"id": "14"}, "layout": {"subject_bids_id": "sub-14"},
        "resolved_paths": {"derivatives_root": str(tmp_path / "derivatives"), "data_root": str(tmp_path / "raw"),
            "derivatives_preprocessing_dir": str(tmp_path / "derivatives/preprocessing/sub-14/ses-01"),
            "derivatives_reports_dir": str(tmp_path / "derivatives/reports/sub-14/ses-01")}}
    records = [AcquisitionRunRecord("source", f"sub14_run{n:02}_test.fif", f"sub14_run{n:02}_test", n, "source") for n in (2, 3)]
    times = np.arange(10000) / 1000.
    signal = np.sin(times * 2 * np.pi) * .1
    signal += 4 * np.exp(-.5 * ((times - 2.) / .08) ** 2)
    raw = mne.io.RawArray(np.vstack([signal, -signal]), mne.create_info(["ICA000", "ICA001"], 1000., "misc"),
                         first_samp=1397000, verbose="ERROR")
    raw._data[:, 3051:3053] = np.nan
    cache = {"model": SimpleNamespace(n_components_=2), "native_sources": raw,
        "metadata": {"fingerprint": "cache", "model_sha256": "model", "sources_stamp": {"name": "sources"}, "converged": True}}
    return cfg, records, raw, cache


def dependencies(tmp_path, monkeypatch):
    """Mock only validated file/model dependencies; exercise real encoding, HTML and publication."""
    cfg, records, raw, cache = example(tmp_path)
    snapshots = {r.source_filename: {"source": r.source_filename, "choices": None} for r in records}
    monkeypatch.setattr(reports, "input_fingerprints", lambda cfg, r: snapshots[r.source_filename])
    monkeypatch.setattr(reports, "require_review", lambda *a: None)
    monkeypatch.setattr(reports, "has_manual_work", lambda *a: False)
    monkeypatch.setattr(reports, "load_cached_aid", lambda *a: cache)
    # A real, tiny PNG keeps HTML self-contained without requiring MEG geometry in a trace fixture.
    tiny = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+jL1sAAAAASUVORK5CYII="
    monkeypatch.setattr(reports, "component_maps", lambda *a: {"mag": tiny})
    folder = triage.triage_root(cfg) / "component_reports/test-batch"
    return cfg, records, raw, cache, snapshots, folder


def test_peak_view_keeps_native_peak_offsets_and_short_nan_gaps(tmp_path):
    _, _, raw, _ = example(tmp_path)
    original = raw.get_data().copy()
    view = reports.pack_trace_view(raw)
    values = np.frombuffer(base64.b64decode(view["values"]), dtype="<f4").reshape(2, -1)
    offsets = np.frombuffer(base64.b64decode(view["offsets"]), dtype="<u2").reshape(2, -1)
    assert view["step"] == 5 and view["native_first_samp"] == raw.first_samp
    assert np.isnan(values[:, 610]).all()  # Gap does not coincide with a simple every-fifth sample.
    peak = np.argmax(np.abs(original[0, :3000]))
    block = peak // 5
    assert block * 5 + offsets[0, block] == peak
    assert values[0, block] == pytest.approx(1.)
    assert values[1, block] == pytest.approx(-1.)
    np.testing.assert_array_equal(raw.get_data(), original)


def test_new_component_choices_do_not_invalidate_read_only_trace_exports(tmp_path):
    cfg, records, _, _ = example(tmp_path)
    before = reports.input_fingerprints(cfg, records[0])
    path = triage.selection_path(cfg, records[0]); path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"status": "selected", "blink_components": [0]}))
    assert reports.input_fingerprints(cfg, records[0]) == before


def test_final_partial_bin_is_not_mistaken_for_missing_acquired_samples(tmp_path):
    _, _, raw, _ = example(tmp_path)
    raw.crop(tmax=9.997)
    raw._data[:, -1] = 10.
    packed = reports.pack_trace_view(raw)
    values = np.frombuffer(base64.b64decode(packed["values"]), dtype="<f4").reshape(2, -1)
    offsets = np.frombuffer(base64.b64decode(packed["offsets"]), dtype="<u2").reshape(2, -1)
    assert values[0, -1] == pytest.approx(1.)
    assert (packed["n_bins"] - 1) * 5 + offsets[0, -1] == raw.n_times - 1


def test_one_participant_html_contains_both_runs_and_resumes_without_rereading_fif(tmp_path, monkeypatch):
    cfg, records, _, _, snapshots, folder = dependencies(tmp_path, monkeypatch)
    manifest = reports.build_participant(cfg, records, folder)
    book = folder / "sub-14.html"
    assert list(folder.glob("*.html")) == [book]
    text = book.read_text()
    assert all(r.source_filename in text for r in records)
    assert 'id="run-select"' in text and 'id="time-slider"' in text
    assert "fetch(" not in text and "https://" not in text
    assert reports.load_report_manifest(folder / "sub-14_manifest.json") == manifest
    monkeypatch.setattr(reports, "load_cached_aid", lambda *a: pytest.fail("Resume should not reread FIFs"))
    assert reports.build_participant(cfg, records, folder, resume=True) == manifest
    snapshots[records[0].source_filename] = {"source": "changed"}
    with pytest.raises(RuntimeError, match="inputs changed"):
        reports.build_participant(cfg, records, folder, resume=True)
    # Existing report remains byte-for-byte intact after rejected resume.
    assert book.read_text() == text


@pytest.mark.parametrize("empty", [False, True])
def test_viewer_controls_in_hermetic_dom_without_browser_or_network(tmp_path, empty):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is optional for hermetic viewer control checks.")
    _, records, raw, _ = example(tmp_path)
    payload = {"participant": "P14", "report_id": "synthetic-only", "runs": []}
    if not empty:
        for record in records:
            payload["runs"].append({"status": "ready", "source_filename": record.source_filename,
                "acquisition_run": record.acquisition_run_label, "n_components": 2,
                "sfreq": 1000., "duration_sec": 10., "converged": True,
                "display": reports.pack_trace_view(raw), "maps": {"mag": "synthetic-map"}})
    path = tmp_path / "synthetic.json"; path.write_text(json.dumps(payload))
    harness = Path(__file__).with_name("ocular_report_dom_harness.js")
    result = subprocess.run([node, str(harness), str(reports.VIEWER), str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_manifest_rejects_changed_html_and_script_like_notes_are_escaped(tmp_path, monkeypatch):
    cfg, records, _, _, _, folder = dependencies(tmp_path, monkeypatch)
    reports.build_participant(cfg, records, folder)
    path = folder / "sub-14_manifest.json"
    with (folder / "sub-14.html").open("a") as handle:
        handle.write("changed")
    with pytest.raises(RuntimeError, match="HTML changed"):
        reports.load_report_manifest(path)
    text = reports.render_html({"participant": "</script><script>alert(1)</script>", "report_id": "test", "runs": []})
    assert "</script><script>alert(1)" not in text
    assert "\\u003c/script" in text


def test_saved_choices_and_manual_work_are_listed_but_never_rebuilt(tmp_path, monkeypatch):
    cfg, records, _, _, _, folder = dependencies(tmp_path, monkeypatch)
    monkeypatch.setattr(reports, "has_manual_work", lambda cfg, record: record == records[0])
    path = triage.selection_path(cfg, records[1]); path.parent.mkdir(parents=True)
    saved = {"status": "selected", "blink_components": [0], "eye_movement_components": []}
    path.write_text(json.dumps(saved)); before = path.read_bytes()
    monkeypatch.setattr(reports, "load_cached_aid", lambda *a: pytest.fail("Existing choices must be preserved"))
    manifest = reports.build_participant(cfg, records, folder)
    assert [r["status"] for r in manifest["runs"]] == ["manual_work_preserved", "choices_already_saved"]
    assert reports.decision_rows(manifest) == []
    assert path.read_bytes() == before
    reports.update_index(cfg, folder)
    assert "sub-14.html" in (folder / "index.html").read_text()


def test_index_repair_handles_denied_replace_without_rewriting_participant_book(tmp_path, monkeypatch):
    """Reproduce the cluster/NAS index-only error without relaxing immutable files."""
    cfg, records, _, _, _, folder = dependencies(tmp_path, monkeypatch)
    reports.build_participant(cfg, records, folder)
    reports.update_index(cfg, folder)
    book = folder / "sub-14.html"; before = book.read_bytes()
    index = folder / "index.html"; index.write_text("stale index")
    original_inode = index.stat().st_ino
    replace = reports.os.replace
    def denied(source, destination):
        if Path(destination).name in ("index.html", "decisions_template.csv"):
            raise PermissionError("shared storage denies replacement")
        return replace(source, destination)
    monkeypatch.setattr(reports.os, "replace", denied)
    monkeypatch.setattr(reports, "load_cached_aid", lambda *a: pytest.fail("Index repair must not read ICA"))
    monkeypatch.setattr(reports, "subject_config", lambda *a: cfg)
    args = argparse.Namespace(subjects=[14], label="test-batch", index_only=True)
    assert reports.run_export(args) == 0
    assert "sub-14.html" in index.read_text()
    assert '<meta http-equiv="refresh" content="30">' in index.read_text()
    assert 'http-equiv="refresh"' not in book.read_text()
    assert index.stat().st_ino == original_inode
    assert book.read_bytes() == before
    with pytest.raises(ValueError, match="restricted"):
        reports.atomic_bytes(cfg, book, b"must not write", immutable=True, locked_index_fallback=True)
    assert book.read_bytes() == before


def test_permission_fallback_does_not_hide_unwritable_index(tmp_path, monkeypatch):
    cfg, _, _, _, _, folder = dependencies(tmp_path, monkeypatch)
    target = folder / "index.html"; target.parent.mkdir(parents=True); target.write_text("unchanged")
    monkeypatch.setattr(reports.os, "replace", lambda *a: (_ for _ in ()).throw(PermissionError("replacement denied")))
    original_open = Path.open
    def denied_open(path, *args, **kwargs):
        if path == target and args and args[0] == "r+b":
            raise PermissionError("writing also denied")
        return original_open(path, *args, **kwargs)
    monkeypatch.setattr(Path, "open", denied_open)
    with pytest.raises(PermissionError, match="writing also denied"):
        reports.atomic_bytes(cfg, target, b"new", locked_index_fallback=True)
    assert target.read_text() == "unchanged"


# -----------------------------------------------------------------------------
# Human notes: explicit syntax and preflight before any new choice is persisted
# -----------------------------------------------------------------------------

def write_notes(path, rows):
    """Use real CSV quoting for multi-component fields."""
    fields = ["participant", "run", "source_filename", "blink", "eye_movement", "notes"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def note(run=2, blink="0", eye="NONE", notes="inspected"):
    return {"participant": "P14", "run": f"R{run:02}", "source_filename": "", "blink": blink,
            "eye_movement": eye, "notes": notes}


@pytest.mark.parametrize("row", [note(blink="", eye="NONE"), note(blink="NONE", eye="DEFER"),
                                    note(blink="QUIT"), note(blink="not sure")])
def test_note_parser_rejects_partial_or_implicit_choices(tmp_path, row):
    path = tmp_path / "notes.csv"; write_notes(path, [row])
    with pytest.raises(ValueError):
        importer.parse_table(path)


def test_note_parser_retains_none_defer_and_ignores_unfinished_rows(tmp_path):
    path = tmp_path / "notes.csv"
    write_notes(path, [note(blink="ICA000,1"), note(3, "DEFER", "DEFER"), note(4, "", "")])
    rows = importer.parse_table(path)
    assert len(rows) == 2 and rows[0]["blink"] == "ICA000,1" and rows[1]["blink"] == "DEFER"
    write_notes(path, [note(), note()])
    with pytest.raises(ValueError, match="duplicate"):
        importer.parse_table(path)


def import_dependencies(tmp_path, monkeypatch):
    """Generate real manifests and preserve real atomic choice writes in the importer."""
    cfg, records, _, cache, snapshots, folder = dependencies(tmp_path, monkeypatch)
    reports.build_participant(cfg, records, folder)
    monkeypatch.setattr(importer, "subject_config", lambda *a: cfg)
    monkeypatch.setattr(importer, "require_review", lambda *a: None)
    monkeypatch.setattr(importer, "discover_raw_run_records", lambda *a: records)
    monkeypatch.setattr(importer, "has_manual_work", lambda *a: False)
    monkeypatch.setattr(importer, "load_cached_aid", lambda *a: cache)
    monkeypatch.setattr(importer, "input_fingerprints", lambda cfg, r: snapshots[r.source_filename])
    path = tmp_path / "notes.csv"
    args = argparse.Namespace(reports_dir=str(folder), decisions=str(path), reviewer="User", apply=False)
    return cfg, records, cache, snapshots, args, path


def test_import_validates_then_writes_existing_choice_format_and_is_idempotent(tmp_path, monkeypatch):
    cfg, records, _, _, args, path = import_dependencies(tmp_path, monkeypatch)
    write_notes(path, [note(blink="0,1")])
    assert importer.run_import(args) == 0
    target = triage.selection_path(cfg, records[0])
    assert not target.exists()
    args.apply = True
    assert importer.run_import(args) == 0
    before = target.read_bytes()
    saved = json.loads(before)
    assert saved["blink_components"] == [0, 1] and saved["eye_movement_components"] == []
    assert saved["purpose"] == "ocular_component_selection_only"
    assert "report_review" in saved
    assert importer.run_import(args) == 0
    assert target.read_bytes() == before
    write_notes(path, [note(blink="1")])
    with pytest.raises(RuntimeError, match="conflict"):
        importer.run_import(args)
    assert target.read_bytes() == before


def test_invalid_later_row_or_changed_model_prevents_every_new_save(tmp_path, monkeypatch):
    cfg, records, cache, _, args, path = import_dependencies(tmp_path, monkeypatch)
    args.apply = True
    write_notes(path, [note(), note(3, "2")])
    with pytest.raises(ValueError):
        importer.run_import(args)
    assert not triage.selection_path(cfg, records[0]).exists()
    write_notes(path, [note()])
    stale = deepcopy(cache); stale["metadata"]["model_sha256"] = "changed"
    monkeypatch.setattr(importer, "load_cached_aid", lambda *a: stale)
    with pytest.raises(RuntimeError, match="binding changed"):
        importer.run_import(args)
    assert not triage.selection_path(cfg, records[0]).exists()


def test_manual_review_or_source_mismatch_blocks_import(tmp_path, monkeypatch):
    cfg, records, _, _, args, path = import_dependencies(tmp_path, monkeypatch)
    row = note(); row["source_filename"] = "another_source.fif"
    write_notes(path, [row])
    with pytest.raises(ValueError, match="filename"):
        importer.run_import(args)
    write_notes(path, [note()])
    monkeypatch.setattr(importer, "has_manual_work", lambda *a: True)
    with pytest.raises(RuntimeError, match="Manual work"):
        importer.run_import(args)
    assert not triage.selection_path(cfg, records[0]).exists()


def test_cluster_site_text_can_differ_while_current_model_identity_is_revalidated(tmp_path, monkeypatch):
    cfg, records, _, snapshots, args, path = import_dependencies(tmp_path, monkeypatch)
    write_notes(path, [note()])
    snapshots[records[0].source_filename] = {**snapshots[records[0].source_filename],
        "config_site": {"sha256": "local-site-paths"}}
    assert importer.run_import(args) == 0
    assert not triage.selection_path(cfg, records[0]).exists()
    snapshots[records[0].source_filename]["source"] = "changed-native-evidence"
    with pytest.raises(RuntimeError, match="Inputs changed"):
        importer.run_import(args)
