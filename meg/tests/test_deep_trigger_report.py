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
- distinct aligned binary labels for reviewed pairs that retain both pulses,
- conservative next-frame gaze-break adjudication without changing raw status,
- confirmed unavailable-tail adjudication without creating observed events,
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

from deep_trigger_report import (  # noqa: E402
    _add_transition_alignment_references,
    _annotate_gaze_break_exclusions,
    _render_transition_plots,
    build_anomaly_rows,
    build_legacy_alignment_rows,
    generate_legacy_deep_report,
)


def _gaze_break_adjudication_rows(
    *,
    missing_code=24,
    missing_reason="",
    partner_status="MATCH",
    partner_frame=13,
    partner_trial=4,
    include_condition_onset=True,
):
    """Build a tiny reported trial with a missing event and a later gaze break."""
    rows = []
    if include_condition_onset:
        rows.append(
            {
                "align_row": 1,
                "status": "MATCH",
                "reported_idx": 1,
                "reported_block": 1,
                "reported_run": 1,
                "reported_trial": 4,
                "reported_replay_status": "original",
                "reported_source_csv": "block01_run01.csv",
                "reported_trigger": 71,
                "reported_label": "trial_condition_onset",
                "reported_frame": 1,
                "reported_seconds": 1.0,
                "observed_trigger": 71,
            }
        )
    rows.extend(
        [
            {
                "align_row": 2,
                "status": "REPORTED_ONLY",
                "reported_idx": 2 if include_condition_onset else 1,
                "reported_block": 1,
                "reported_run": 1,
                "reported_trial": 4,
                "reported_replay_status": "original",
                "reported_source_csv": "block01_run01.csv",
                "reported_trigger": missing_code,
                "reported_label": "trial_sequence_identity",
                "reported_frame": 12,
                "reported_seconds": 1.100,
                "reported_only_reason": missing_reason,
                "observed_trigger": "",
                "stage1_adjudication": "REVIEW_REQUIRED",
                "adjudication_partner_align_row": "",
                "adjudication_gap_ms": "",
                "adjudication_detail": "",
            },
            {
                "align_row": 3,
                "status": partner_status,
                "reported_idx": 3 if include_condition_onset else 2,
                "reported_block": 1,
                "reported_run": 1,
                "reported_trial": partner_trial,
                "reported_replay_status": "original",
                "reported_source_csv": "block01_run01.csv",
                "reported_trigger": 150,
                "reported_label": "gaze_break",
                "reported_frame": partner_frame,
                "reported_seconds": 1.1082,
                "observed_trigger": 150,
            },
        ]
    )
    return rows


def test_unavailable_tail_plot_context_stays_in_truncated_acquisition(tmp_path):
    """Missing run-1 data must not be illustrated using run-2 pulses."""
    events, alignment, _, _ = _canonical_fixture(tmp_path / "run01_raw.fif")
    events[1]["source_fif_path"] = str(tmp_path / "run02_raw.fif")
    events[1]["source_fif"] = "run02_raw.fif"
    events[1]["acquisition_run"] = 2
    alignment[1]["acquisition_run"] = 2
    missing = dict(alignment[0])
    missing.update({"alignment_index": 2, "status": "REPORTED_ONLY",
                    "observed_event_index": None, "observed_trigger_code": None,
                    "source_fif": None, "acquisition_run": None, "raw_sample_frame": None,
                    "constrained_acquisition_run": 1,
                    "reported_only_reason": "recording_tail_unavailable"})
    alignment[1]["alignment_index"] = 3
    alignment.insert(1, missing)
    anomalies = build_anomaly_rows(build_legacy_alignment_rows(alignment, events), events)
    tail = next(row for row in anomalies if row["reported_only_reason"])
    assert tail["anchor_source_fif_path"] == str(tmp_path / "run01_raw.fif")
    assert tail["anchor_source"] == "recorded_cutoff_context"
    assert "no recorded waveform" in tail["detail"]
    assert tail["stage1_adjudication"] == "RESOLVED_RECORDING_TAIL_UNAVAILABLE"


@pytest.mark.parametrize(
    "reason, constrained_run, expected",
    [
        (None, None, "REVIEW_REQUIRED"),
        ("recording_tail_unavailable", None, "REVIEW_REQUIRED"),
        ("unexplained", 1, "REVIEW_REQUIRED"),
        ("recording_tail_unavailable", 1, "RESOLVED_RECORDING_TAIL_UNAVAILABLE"),
    ],
)
def test_recording_tail_resolution_requires_canonical_boundary_evidence(
    tmp_path, reason, constrained_run, expected
):
    """Resolve only confirmed unavailable tails; preserve missing observations."""
    events, alignment, _, _ = _canonical_fixture(tmp_path / "run01_raw.fif")
    missing = dict(alignment[0])
    missing.update({
        "status": "REPORTED_ONLY", "observed_event_index": None,
        "observed_trigger_code": None, "source_fif": None,
        "acquisition_run": None, "raw_sample_frame": None,
        "reported_only_reason": reason,
        "constrained_acquisition_run": constrained_run,
    })
    before = dict(missing)
    result = build_legacy_alignment_rows([missing], events)[0]
    assert result["stage1_adjudication"] == expected
    assert result["status"] == "REPORTED_ONLY"
    assert result["observed_idx"] is None
    assert result["raw_sample_frame"] is None
    assert missing == before


def test_next_frame_matched_gaze_break_resolves_trial_only_for_review():
    """Keep raw REPORTED_ONLY while marking its excluded trial adjudicated."""
    rows = _gaze_break_adjudication_rows()

    _annotate_gaze_break_exclusions(
        rows,
        gaze_break_codes=[150],
        gaze_break_exclusion_enabled=True,
    )

    assert rows[1]["status"] == "REPORTED_ONLY"
    assert rows[1]["stage1_adjudication"] == (
        "RESOLVED_TRIAL_EXCLUDED_GAZE_BREAK"
    )
    assert rows[1]["adjudication_partner_align_row"] == 3
    assert rows[1]["adjudication_gap_ms"] == "8.200"


@pytest.mark.parametrize(
    "row_options, exclusion_enabled",
    [
        ({"missing_code": 150}, True),
        ({"partner_frame": 14}, True),
        ({"partner_trial": 5}, True),
        ({"partner_status": "REPORTED_ONLY"}, True),
        ({"missing_reason": "recording_tail_unavailable"}, True),
        ({"include_condition_onset": False}, True),
        ({}, False),
    ],
)
def test_gaze_break_adjudication_leaves_unproven_rows_for_review(
    row_options, exclusion_enabled
):
    """Do not auto-resolve missing breaks, tails, or trials without full evidence."""
    rows = _gaze_break_adjudication_rows(**row_options)

    _annotate_gaze_break_exclusions(
        rows,
        gaze_break_codes=[150],
        gaze_break_exclusion_enabled=exclusion_enabled,
    )

    missing = next(row for row in rows if row["status"] == "REPORTED_ONLY")
    assert missing["stage1_adjudication"] == "REVIEW_REQUIRED"


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
            "sample_start_raw": 101,
            "raw_sample_frame": 101,
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
            "raw_sample_frame": 101,
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
    assert "Alignment rows: kept 1" in report
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


def test_transition_warning_references_adjacent_reported_only_dropped_code(tmp_path):
    """Link both rows when a collapsed pulse matches the next reported code."""
    source = str(tmp_path / "sub37_run01_raw.fif")
    events = [
        {
            "event_index": 1,
            "trigger_code": 20,
            "source_fif_path": source,
            "sample_start_raw": 100,
        }
    ]
    alignment = [
        {
            "alignment_index": 134,
            "reported_trigger_code": 20,
            "observed_event_index": 1,
        },
        {
            "alignment_index": 135,
            "reported_trigger_code": 150,
            "observed_event_index": None,
        },
    ]
    warning = {
        "warning_index": 1,
        "source_fif_path": source,
        "kept_start_sample_raw": 100,
        "kept_code": 20,
        "dropped_code": 150,
    }

    referenced = _add_transition_alignment_references([warning], events, alignment)

    assert referenced[0]["kept_align_row"] == 134
    assert referenced[0]["dropped_align_row"] == 135


def test_reviewed_transition_reports_both_retained_alignment_rows(tmp_path):
    """Show separate row references when a reviewed pair retains both codes."""
    source = str(tmp_path / "sub19_run04_raw.fif")
    events = [
        {
            "event_index": 1,
            "trigger_code": 16,
            "source_fif_path": source,
            "sample_start_raw": 25103,
        },
        {
            "event_index": 2,
            "trigger_code": 150,
            "source_fif_path": source,
            "sample_start_raw": 25104,
        },
    ]
    alignment = [
        {
            "alignment_index": 2163,
            "reported_trigger_code": 16,
            "observed_event_index": 1,
        },
        {
            "alignment_index": 2164,
            "reported_trigger_code": 150,
            "observed_event_index": 2,
        },
    ]
    warning = {
        "warning_index": 15,
        "source_fif_path": source,
        "kept_start_sample_raw": 25103,
        "kept_code": 16,
        "dropped_start_sample_raw": 25104,
        "dropped_code": 150,
        "retained_both": True,
    }

    referenced = _add_transition_alignment_references([warning], events, alignment)

    assert referenced[0]["kept_align_row"] == 2163
    assert referenced[0]["dropped_align_row"] == 2164


def test_transition_plot_title_shows_both_alignment_rows(monkeypatch, tmp_path):
    """Identify both retained pulses in the title and aligned binary panel."""
    captured = {}

    class _ClosableRaw:
        """Minimal cached FIF handle for the plot-renderer test."""

        @staticmethod
        def close():
            return None

    monkeypatch.setattr(
        "deep_trigger_report._save_stim_plot",
        lambda **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "deep_trigger_report.mne.io.read_raw_fif",
        lambda *args, **kwargs: _ClosableRaw(),
    )
    warning = {
        "warning_index": 2,
        "source_fif": "sub37_run01_051026.fif",
        "source_fif_path": str(tmp_path / "sub37_run01_051026.fif"),
        "order_pattern": "20_then_150",
        "kept_code": 20,
        "dropped_code": 150,
        "retained_both": True,
        "retained_codes": [20, 150],
        "kept_align_row": 134,
        "dropped_align_row": 135,
        "dropped_start_sample_raw": 285666,
        "dropped_end_sample_raw": 285668,
    }

    _render_transition_plots(
        [warning],
        plots_dir=tmp_path / "plots",
        stim_channel="STI101",
        half_window_ms=3.0,
    )

    assert "keep1 134; keep2 135" in captured["title"]
    assert "retained both codes 20, 150" in captured["title"]
    assert captured["right_annotation"] == (
        "keep1: 00010100 (20)\nkeep2: 10010110 (150)"
    )
    assert "_both_retained.png" in str(captured["output_path"])


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
