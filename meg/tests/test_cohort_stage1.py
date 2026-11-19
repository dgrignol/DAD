"""Exercise new cohort timing and provenance rules without touching raw data.

Usage: python -m pytest meg/tests/test_cohort_stage1.py -q

Synthetic CSVs and diode samples cover P08's runtime switch, P24's separate
run-5 calibration, and the nonauthoritative session-note extraction. No large
FIF processing or mounted-data writes are part of this test.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"))

from dataset_inventory import read_session_review_clues  # noqa: E402
from p24_acquisition_adapter import build_timing_proposal  # noqa: E402
from trigger_extraction import (  # noqa: E402
    ExpectedTrigger, ObservedWindow, ReportedSegment,
    load_reported_trigger_segments, process_trigger_windows,
)


def _trigger_csv(path: Path) -> Path:
    """Create one reported segment with the production minimum columns."""
    path.write_text("trigger,trial,frame,seconds,label\n71,1,1,1.0,trial\n", encoding="utf-8")
    return path


def test_p08_runtime_changes_only_after_block_one(tmp_path: Path):
    """Accept observed v24/v25 names only in their configured blocks."""
    first = _trigger_csv(tmp_path / "debug_actual_triggers_occlusion_v24_allTrialGuard_sub08_block01_run01.csv")
    second = _trigger_csv(tmp_path / "debug_actual_triggers_occlusion_v25_dotOffsetBlank_sub08_block02_run01.csv")
    selected = {1: "v24_allTrialGuard", 2: "v25_dotOffsetBlank"}
    segments = load_reported_trigger_segments([first, second], subject_id="08",
                                               runtime_version="v25_dotOffsetBlank",
                                               runtime_by_block=selected)
    assert [(item.block, item.run) for item in segments] == [(1, 1), (2, 1)]
    with pytest.raises(ValueError, match="runtime differs"):
        load_reported_trigger_segments([second], subject_id="08",
                                       runtime_version="v25_dotOffsetBlank",
                                       runtime_by_block={2: "v24_allTrialGuard"})


def test_p24_only_run_five_measures_photodiode():
    """The earlier-run estimate is visibly distinct from measured run 5."""
    sfreq = 1000.0
    signal = np.zeros(3000, dtype=float)
    events = []
    for index in range(25):
        sample = 100 + index * 100
        signal[sample + 19:sample + 39] = 0.1
        events.append({"event_index": str(index + 1), "trigger_code": "71",
                       "acquisition_run": "5", "raw_sample_frame": str(sample),
                       "source_fif": "run05.fif", "reported_alignment_status": "MATCH"})
    events.insert(0, {"event_index": "0", "trigger_code": "81",
                      "acquisition_run": "1", "raw_sample_frame": "20",
                      "source_fif": "run01.fif", "reported_alignment_status": "MATCH"})
    rows, summary = build_timing_proposal(events, signal, sfreq, 0)
    assert summary["accepted_run5_events"] == 25
    assert summary["run5_median_lag_ms"] == pytest.approx(19.0)
    assert rows[0]["timing_status"] == "estimated_from_run5_review_required"
    assert rows[0]["proposed_visual_sample"] == 39
    assert rows[1]["timing_status"] == "measured_run5"


def test_session_notes_include_only_review_fields(tmp_path: Path):
    """Do not copy demographics or compensation into derivative reports."""
    notes = tmp_path / "Notes_from_session"
    notes.mkdir()
    with (notes / "Demographics & Sessions log - Session Log.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Participant_ID", "Compensation", "Technical_Issues"])
        writer.writeheader()
        writer.writerow({"Participant_ID": "P24", "Compensation": "secret",
                         "Technical_Issues": "check settings"})
    (notes / "Demographics & Sessions log - Runs Details.csv").write_text(
        "Participant 24,Run 1,Run 2\nheadcoil,poor,okay\nParticipant 25,Run 1\n",
        encoding="utf-8",
    )
    clues = read_session_review_clues(tmp_path, "24")
    assert clues["session"] == {"Technical_Issues": "check settings"}
    assert clues["run_details"] == [{"topic": "headcoil", "runs": {"run_1": "poor", "run_2": "okay"}}]


def test_truncated_acquisition_does_not_consume_next_block():
    """Reset both transition and alignment cursors after missing recorded data.

    The first acquisition ends before its second trial. The next acquisition
    begins with a one-sample line transition in either order; its expected
    code must come from block 2 rather than block 1's missing tail.
    """
    profile = {"identity_range": [1, 40],
               "condition_triggers": {"visible": 71, "nondeviant": 81, "deviant": 91},
               "acquisition_block_map": {1: [1], 2: [2]},
               "truncated_acquisition_runs": [1]}
    settings = {"gaze_break_code": 150, "esc_termination_code": 152}
    codes_by_block = {1: [71, 1, 111, 112, 114, 115, 81, 2, 111, 112, 114, 115],
                      2: [91, 3, 111, 112, 114, 115]}
    segments = []
    for block, codes in codes_by_block.items():
        path = Path(f"/tmp/block{block}_run01.csv")
        segment = ReportedSegment(block, block, 1, path)
        segment.events = [ExpectedTrigger(index, code, str(path), block, 1, index)
                          for index, code in enumerate(codes, start=1)]
        segments.append(segment)
    for pair in ((80, 91), (91, 80)):
        windows = []
        for run, codes in ((1, codes_by_block[1][:6]), (2, codes_by_block[2])):
            for index, code in enumerate(codes, start=1):
                start = index * 100
                parts = [(code, start, start + 2)]
                if run == 2 and index == 1:
                    parts = [(pair[0], start, start + 1), (pair[1], start + 1, start + 2)]
                for trig, first, last in parts:
                    windows.append(ObservedWindow(
                        len(windows) + 1, trig, first, last, first / 1000, last / 1000,
                        first, last, first / 1000, last / 1000, last - first,
                        f"run{run}.fif", f"/tmp/run{run}.fif", run,
                    ))
        events, warnings, _, alignment, summary = process_trigger_windows(
            windows, segments, trigger_config=settings, subject_profile=profile,
        )
        assert summary["status_counts"].get("MISMATCH", 0) == 0
        assert summary["status_counts"]["REPORTED_ONLY"] == 6
        assert summary["acquisition_boundaries"][0]["recording_tail_unavailable_count"] == 6
        assert all(row["behavioral_block"] == 2 for row in events if row["acquisition_run"] == 2)
        assert len(events) == 12  # The six unavailable events were not synthesized.
        assert [row["event_index"] for row in events] == list(range(1, 13))
        assert warnings[0]["kept_code"] == 91
        assert not warnings[0]["unresolved_transition"]
        assert all(row["raw_sample_frame"] is None for row in alignment if row["status"] == "REPORTED_ONLY")


def test_current_cohort_defaults_exclude_early_experiment_versions():
    """Working defaults exclude P15 while accepted Stage 1-3 stay readable."""
    import yaml
    from pipeline_config import ANALYSIS_SUBJECT_IDS, WORKING_SUBJECT_IDS
    config = yaml.safe_load((Path(__file__).resolve().parents[1] / "params/pipeline_default.yaml").read_text())
    assert list(WORKING_SUBJECT_IDS) == config["analysis_cohort"]["included_subjects"]
    assert config["matlab_drsa"]["group_average_subjects"] == list(WORKING_SUBJECT_IDS)
    assert not {"01", "02", "03", "15"} & set(WORKING_SUBJECT_IDS)
    # P15's later exclusion must not rewrite the accepted accounting roster.
    assert "15" in ANALYSIS_SUBJECT_IDS and len(ANALYSIS_SUBJECT_IDS) == 36


def test_stage1_cli_writes_deep_report_for_truncated_recording(tmp_path: Path):
    """Exercise preflight, FIF reconstruction, boundaries, and deep plots end to end."""
    import mne
    import yaml
    meg_root = Path(__file__).resolve().parents[1]
    raw_root = tmp_path / "raw"
    subject = raw_root / "sub37_051026"
    for folder in ("meg", "input_files", "output_files", "eyetracker"):
        (subject / folder).mkdir(parents=True)
    (subject / "eyetracker" / "saved.edf").write_bytes(b"placeholder")
    for name in ("MovDot_Sub37_V30_rescueTraject.mat", "Sub37_TrialStruct_v23_rescueTraject.mat"):
        (subject / "input_files" / name).write_bytes(b"placeholder")
    expected = {1: [71, 1, 111, 112, 114, 115, 81, 2, 111, 112, 114, 115],
                2: [91, 3, 111, 112, 114, 115]}
    for run, codes in expected.items():
        actual = codes[:6] if run == 1 else codes
        data = np.zeros((1, 1000))
        for index, code in enumerate(actual, start=1):
            data[0, index * 100:index * 100 + 2] = code
        if run == 2:
            data[0, 100] = 80
        raw = mne.io.RawArray(data, mne.create_info(["STI101"], 1000.0, ["stim"]), verbose="ERROR")
        raw.save(subject / "meg" / f"sub37_run{run:02d}_051026.fif", overwrite=False, verbose="ERROR")
        (subject / "output_files" / f"MoveDot1_occlusion_v25_dotOffsetBlank_SUB37_BLOCK{run:02d}.mat").write_bytes(b"placeholder")
        path = subject / "output_files" / f"debug_actual_triggers_occlusion_v25_dotOffsetBlank_sub37_block{run:02d}_run01.csv"
        path.write_text("trigger,trial,frame,seconds,label\n" + "".join(
            f"{code},1,1,{index / 10},test\n" for index, code in enumerate(codes, start=1)
        ))
    override = tmp_path / "sub-37.yaml"
    override.write_text(yaml.safe_dump({
        "subject": {"id": "37"},
        "project": {"data_root": str(raw_root), "derivatives_root": str(tmp_path / "derivatives")},
    }))
    # Profile maps merge by key across config layers. Build a reduced default
    # fixture so the real run-3..6 constraints do not enter this two-file test.
    defaults = yaml.safe_load((meg_root / "params/pipeline_default.yaml").read_text())
    defaults["subject_profiles"]["37"].update({
        "completed_behavioral_blocks": [1, 2], "acquisition_block_map": {1: [1], 2: [2]},
    })
    default_path = tmp_path / "pipeline_default.yaml"
    default_path.write_text(yaml.safe_dump(defaults))
    env = dict(os.environ, MPLBACKEND="Agg", MPLCONFIGDIR=str(tmp_path / "mpl"), PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run([
        sys.executable, str(meg_root / "code/preprocessing_basic/run_pipeline.py"),
        "--config", str(default_path),
        "--subject-config", str(override), "--stages", "trigger_extraction", "--write-inventory",
    ], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    reports = tmp_path / "derivatives/reports/sub-37/ses-01"
    report = json.loads((reports / "02_trigger_extraction_custom_report.json").read_text())
    assert report["checks"]["alignment_status_counts"].get("MISMATCH", 0) == 0
    assert report["stats"]["acquisition_boundaries"][0]["recording_tail_unavailable_count"] == 6
    deep = report["files"]["deep_trigger_report"]
    assert deep["status"] == "ok"
    # Accepted tail adjudication retains all six log-only events but shares one
    # recording-cutoff context plot; there is no waveform for unrecorded pulses.
    assert deep["anomaly_count"] == 6
    assert deep["reported_only_resolved_recording_tail"] == 6
    assert deep["reported_only_review_required"] == 0
    assert deep["anomaly_plot_count"] == 1
    assert deep["plot_errors"] == []
    assert deep["transition_plot_count"] == 1
    assert Path(deep["report_md"]).is_file()
