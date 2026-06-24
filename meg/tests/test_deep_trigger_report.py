"""
test_deep_trigger_report.py

Purpose
Validate legacy-format deep trigger report generation from canonical step-02
events and alignment data.

Coverage
- legacy bundle filenames and Markdown sections,
- full alignment and anomaly CSV creation,
- anomaly and transition PNG rendering from a bounded synthetic STI101 FIF,
- preservation of acquisition-run and raw-sample provenance,
- rejection of any report output inside the configured raw tree.

Usage
Run only this file:
  python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_deep_trigger_report.py -q
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import mne
import numpy as np
import pytest


CODE_DIR = (
    Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
)
sys.path.insert(0, str(CODE_DIR))

from deep_trigger_report import generate_legacy_deep_report  # noqa: E402


def _write_synthetic_stim_fif(path: Path) -> None:
    """Write a tiny read-only-style source FIF with two STI101 pulses."""
    sfreq = 1000.0
    data = np.zeros((1, 400), dtype=float)
    data[0, 100:102] = 31
    data[0, 200:202] = 152
    info = mne.create_info(["STI101"], sfreq=sfreq, ch_types=["stim"])
    raw = mne.io.RawArray(data, info, verbose="ERROR")
    path.parent.mkdir(parents=True, exist_ok=True)
    raw.save(path, overwrite=True, verbose="ERROR")


def _canonical_fixture(fif_path: Path):
    """Build canonical events/alignment/warnings for deep-report tests."""
    events = [
        {
            "event_index": 1,
            "trigger_code": 31,
            "time_start_sec": 0.100,
            "time_end_sec": 0.102,
            "width_samples": 2,
            "source_fif": fif_path.name,
            "source_fif_path": str(fif_path),
            "acquisition_run": 1,
            "raw_sample_frame": 100,
        },
        {
            "event_index": 2,
            "trigger_code": 152,
            "time_start_sec": 0.200,
            "time_end_sec": 0.202,
            "width_samples": 2,
            "source_fif": fif_path.name,
            "source_fif_path": str(fif_path),
            "acquisition_run": 1,
            "raw_sample_frame": 200,
        },
    ]
    alignment = [
        {
            "alignment_index": 1,
            "status": "MATCH_WARNING",
            "behavioral_block": 1,
            "behavioral_run": 1,
            "behavioral_trial": 1,
            "reported_trigger_code": 31,
            "reported_label": "trial_condition_onset",
            "reported_frame": 1,
            "reported_seconds": 1.0,
            "reported_source_csv": "/tmp/debug_block01_run01.csv",
            "reported_row_in_csv": 1,
            "reported_replay_status": "original",
            "observed_event_index": 1,
            "observed_trigger_code": 31,
            "source_fif": fif_path.name,
            "acquisition_run": 1,
            "raw_sample_frame": 100,
            "behavioral_context_source": "reported_alignment",
        },
        {
            "alignment_index": 2,
            "status": "OBSERVED_ONLY",
            "behavioral_block": 1,
            "behavioral_run": 1,
            "behavioral_trial": None,
            "reported_trigger_code": None,
            "reported_label": None,
            "reported_frame": None,
            "reported_seconds": None,
            "reported_source_csv": None,
            "reported_row_in_csv": None,
            "reported_replay_status": None,
            "observed_event_index": 2,
            "observed_trigger_code": 152,
            "source_fif": fif_path.name,
            "acquisition_run": 1,
            "raw_sample_frame": 200,
            "behavioral_context_source": "inferred_esc_from_preceding_segment",
        },
    ]
    warnings = [
        {
            "warning_index": 1,
            "source_fif": fif_path.name,
            "source_fif_path": str(fif_path),
            "acquisition_run": 1,
            "order_pattern": "9_then_31",
            "dropped_code": 9,
            "kept_code": 31,
            "expected_trigger_code": 31,
            "decision_source": "expected_trigger",
            "unresolved_transition": False,
            "transition_resolution_confidence": 1.0,
            "dropped_start_sample_raw": 100,
            "dropped_end_sample_raw": 101,
            "dropped_width_samples": 1,
            "kept_start_sample_raw": 101,
            "kept_end_sample_raw": 102,
            "kept_width_samples": 1,
        }
    ]
    segment_summaries = [
        {
            "behavioral_block": 1,
            "behavioral_run": 1,
            "alignment_score": 3,
        }
    ]
    return events, alignment, warnings, segment_summaries


def test_generate_legacy_deep_report_bundle(tmp_path):
    """Generate the complete legacy bundle from canonical step-02 evidence."""
    raw_root = tmp_path / "raw"
    fif_path = raw_root / "sub04_220526" / "meg" / "sub04_run01_220526.fif"
    _write_synthetic_stim_fif(fif_path)
    events, alignment, warnings, segment_summaries = _canonical_fixture(
        fif_path
    )
    output_dir = (
        tmp_path
        / "derivatives"
        / "reports"
        / "sub-04"
        / "ses-01"
        / "02_triggers_extraction_custum_report_deeper"
    )

    result = generate_legacy_deep_report(
        subject_id="04",
        raw_root=raw_root,
        output_dir=output_dir,
        stim_channel="STI101",
        events=events,
        alignment_rows=alignment,
        transition_warnings=warnings,
        segments=[object()],
        segment_summaries=segment_summaries,
    )

    assert result["status"] == "ok"
    assert result["alignment_status_counts"] == {
        "MATCH_WARNING": 1,
        "OBSERVED_ONLY": 1,
    }
    assert result["anomaly_count"] == 1
    assert result["anomaly_plot_count"] == 1
    assert result["transition_plot_count"] == 1

    report_path = Path(result["report_md"])
    assert report_path.name == "sub04_all_blocks_reported_vs_observed.md"
    report = report_path.read_text(encoding="utf-8")
    assert "# Subject Trigger Confrontation Report" in report
    assert "## Block Starts" in report
    assert "## Run Starts" in report
    assert "## Anomaly Plots" in report
    assert "## Transition Warnings" in report
    assert "## Full Alignment Table" in report
    assert "OBSERVED_ONLY" in report
    assert "sub04_run01_220526.fif" in report

    alignment_path = Path(result["alignment_csv"])
    with alignment_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["reported_frame"] == "1"
    assert rows[0]["acquisition_run"] == "1"
    assert rows[1]["raw_sample_frame"] == "200"
    assert (
        rows[1]["behavioral_context_source"]
        == "inferred_esc_from_preceding_segment"
    )

    assert Path(result["anomalies_csv"]).exists()
    assert list(Path(result["plots_dir"]).glob("*.png"))
    assert list((Path(result["plots_dir"]) / "transition_pulses").glob("*.png"))


def test_deep_report_rejects_output_inside_raw_tree(tmp_path):
    """The shared guard must reject deep-report directories below raw data."""
    raw_root = tmp_path / "raw"
    fif_path = raw_root / "sub04_220526" / "meg" / "sub04_run01_220526.fif"
    _write_synthetic_stim_fif(fif_path)
    events, alignment, warnings, segment_summaries = _canonical_fixture(
        fif_path
    )

    with pytest.raises(ValueError, match="read-only raw-data tree"):
        generate_legacy_deep_report(
            subject_id="04",
            raw_root=raw_root,
            output_dir=raw_root / "forbidden-report",
            stim_channel="STI101",
            events=events,
            alignment_rows=alignment,
            transition_warnings=warnings,
            segments=[object()],
            segment_summaries=segment_summaries,
        )
