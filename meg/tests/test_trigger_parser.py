"""
test_trigger_parser.py

Purpose
Validate transition resolution, runtime-specific CSV parsing, replay metadata,
and behavioral-segment alignment for
`code/preprocessing_basic/02_trigger_extraction_custom.py`.

Coverage
- expected-trigger-first symmetric transition handling,
- v18 and v20 runtime-specific CSV filename/schema parsing,
- identity code 40,
- gaze break 150, replay start 151, and ESC 152,
- one acquisition FIF containing multiple behavioral blocks,
- an ESC-terminated partial final block/run segment.
- an observed-only ESC inheriting unambiguous same-FIF segment context.

Why this matters
The MEG trigger stream can show short line-update transitions in either order:
`0 -> A -> B -> 0` or `0 -> B -> A -> 0`.
The parser must keep the intended trigger from expected sequence context in
both orders and drop the non-intended trigger.

Usage examples
1) Run this test file only:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_trigger_parser.py -q

2) Run all tests:
   python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests -q

Assumptions
- Pytest is available in the environment.
- Step-02 module exists at:
  /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/code/preprocessing_basic/02_trigger_extraction_custom.py
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


def _load_step02_module():
    """Load the step-02 script dynamically so tests can access its functions/classes."""
    module_path = (
        Path(__file__).resolve().parents[1]
        / "code"
        / "preprocessing_basic"
        / "02_trigger_extraction_custom.py"
    )
    sys.path.insert(0, str(module_path.parent))
    spec = importlib.util.spec_from_file_location("step02_trigger_module", str(module_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module spec from {module_path}")
    module = importlib.util.module_from_spec(spec)
    # Register module before execution so dataclass/type resolution works reliably.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STEP02 = _load_step02_module()


def _window(idx: int, trig: int, start: int, end: int):
    """Build one observed window with deterministic metadata for parser tests."""
    return STEP02.ObservedWindow(
        idx_global=idx,
        trigger=trig,
        start_sample=start,
        end_sample=end,
        start_sec=start / 1000.0,
        end_sec=end / 1000.0,
        start_sample_raw=start,
        end_sample_raw=end,
        start_sec_raw=start / 1000.0,
        end_sec_raw=end / 1000.0,
        width_samples=end - start,
        source_fif="SUB01_block01_meg.fif",
        source_fif_path="/tmp/SUB01_block01_meg.fif",
    )


def _expected(idx: int, trig: int):
    """Build one expected trigger row for expected-trigger-first decisions."""
    return STEP02.ExpectedTrigger(
        idx_global=idx,
        trigger=trig,
        source_csv="debug_actual_triggers_occlusion_v18_rescueTraject_sub01_block01_run01.csv",
        block=1,
        run=1,
        row_in_csv=idx,
    )


def _write_debug_csv(
    tmp_path: Path,
    *,
    runtime: str,
    subject_id: str,
    block: int,
    run: int,
    rows: list[tuple[int, int, int, float, str]],
) -> Path:
    """Write one runtime-shaped synthetic debug CSV with actual saved rows."""
    path = (
        tmp_path
        / (
            f"debug_actual_triggers_occlusion_{runtime}_rescueTraject_"
            f"sub{subject_id}_block{block:02d}_run{run:02d}.csv"
        )
    )
    lines = ["trigger,trial,frame,seconds,label"]
    lines.extend(
        f"{trigger},{trial},{frame},{seconds},{label}"
        for trigger, trial, frame, seconds, label in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _resolved_events(
    codes: list[int],
    *,
    source_fif: str,
    acquisition_run: int,
    seq_max: int,
    condition_triggers: dict[str, int],
    expected: list,
):
    """Resolve separated synthetic windows into authoritative event rows."""
    windows = []
    for index, code in enumerate(codes, start=1):
        start = index * 10
        window = _window(index, code, start, start + 2)
        window.source_fif = source_fif
        window.source_fif_path = f"/tmp/{source_fif}"
        window.acquisition_run = acquisition_run
        windows.append(window)
    events, warnings, _decisions = STEP02.resolve_observed_windows_to_events(
        observed_windows=windows,
        expected_sequence=expected,
        transition_max_samples=1,
        seq_min=1,
        seq_max=seq_max,
        expected_lookahead=6,
        condition_triggers=condition_triggers,
    )
    assert warnings == []
    return events


def test_expected_trigger_first_is_symmetric_for_transition_orders():
    """
    Keep expected code for both contiguous transition orders.

    Construct two transition pairs:
    - `10 -> 31` (keep 31)
    - `31 -> 10` (keep 31)
    plus one non-transition 41 event.
    """
    observed_windows = [
        # Pair 1: unexpected->expected
        _window(1, 10, 100, 101),
        _window(2, 31, 101, 103),
        # Pair 2: expected->unexpected
        _window(3, 31, 200, 202),
        _window(4, 10, 202, 203),
        # Non-transition event
        _window(5, 41, 300, 302),
    ]
    expected_sequence = [
        _expected(1, 31),
        _expected(2, 31),
        _expected(3, 41),
    ]

    events, transition_warnings, decision_counts = STEP02.resolve_observed_windows_to_events(
        observed_windows=observed_windows,
        expected_sequence=expected_sequence,
        transition_max_samples=1,
        seq_min=1,
        seq_max=20,
        expected_lookahead=6,
    )

    assert [row["trigger_code"] for row in events] == [31, 31, 41]
    assert len(transition_warnings) == 2
    assert all(row["decision_source"] == "expected_trigger" for row in transition_warnings)
    assert all(not bool(row["unresolved_transition"]) for row in transition_warnings)
    assert [row["order_pattern"] for row in transition_warnings] == [
        "10_then_31",
        "31_then_10",
    ]
    assert [row["dropped_code"] for row in transition_warnings] == [10, 10]
    assert decision_counts.get("expected_trigger", 0) == 2


def test_fallback_resolution_uses_width_then_unresolved_deterministic():
    """
    Validate fallback path when expected identity is unavailable.

    Case 1 (width fallback):
    - contiguous pair `77 -> 88`
    - first pulse is short (1 sample), second is long (4 samples)
    - expected sequence is empty, so resolver must keep 88 via width_fallback.

    Case 2 (deterministic unresolved fallback):
    - contiguous pair `55 -> 66`
    - both pulses are short (1 sample each)
    - no expected sequence, no asymmetric width evidence
    - resolver keeps first by deterministic fallback and marks unresolved.
    """
    observed_windows = [
        # Width-fallback pair.
        _window(1, 77, 100, 101),  # short
        _window(2, 88, 101, 105),  # long
        # Deterministic unresolved pair (both short).
        _window(3, 55, 200, 201),  # short
        _window(4, 66, 201, 202),  # short
    ]

    events, transition_warnings, decision_counts = STEP02.resolve_observed_windows_to_events(
        observed_windows=observed_windows,
        expected_sequence=[],
        transition_max_samples=1,
        seq_min=1,
        seq_max=20,
        expected_lookahead=6,
    )

    assert len(events) == 2
    assert [row["trigger_code"] for row in events] == [88, 55]

    assert len(transition_warnings) == 2

    width_warning = transition_warnings[0]
    assert width_warning["order_pattern"] == "77_then_88"
    assert width_warning["decision_source"] == "width_fallback"
    assert width_warning["dropped_code"] == 77
    assert width_warning["kept_code"] == 88
    assert not bool(width_warning["unresolved_transition"])
    assert width_warning["transition_resolution_confidence"] == 0.7

    unresolved_warning = transition_warnings[1]
    assert unresolved_warning["order_pattern"] == "55_then_66"
    assert unresolved_warning["decision_source"] == "deterministic_fallback"
    assert unresolved_warning["dropped_code"] == 66
    assert unresolved_warning["kept_code"] == 55
    assert bool(unresolved_warning["unresolved_transition"])
    assert unresolved_warning["transition_resolution_confidence"] == 0.0

    assert decision_counts.get("width_fallback", 0) == 1
    assert decision_counts.get("deterministic_fallback", 0) == 1


@pytest.mark.parametrize(
    ("runtime", "subject_id", "condition_code", "identity_code"),
    [
        ("v18", "01", 31, 20),
        ("v20", "04", 71, 40),
    ],
)
def test_runtime_specific_csv_parser_supports_both_generations(
    tmp_path,
    runtime,
    subject_id,
    condition_code,
    identity_code,
):
    """The profile runtime selects the exact CSV generation and identity range."""
    csv_path = _write_debug_csv(
        tmp_path,
        runtime=runtime,
        subject_id=subject_id,
        block=1,
        run=1,
        rows=[
            (condition_code, 1, 1, 1.0, "trial_condition_onset"),
            (identity_code, 1, 12, 1.1, "trial_sequence_identity"),
        ],
    )

    segments = STEP02.load_reported_trigger_segments(
        [csv_path],
        subject_id=subject_id,
        runtime_version=runtime,
    )

    assert len(segments) == 1
    assert [event.trigger for event in segments[0].events] == [
        condition_code,
        identity_code,
    ]
    assert segments[0].block == 1
    assert segments[0].run == 1


def test_v20_trigger_label_treats_code_31_as_identity():
    """Profile-specific v20 conditions prevent the legacy code-31 label collision."""
    assert STEP02.trigger_label(
        31,
        1,
        40,
        {
            "always_visible": 71,
            "occluded_nondeviant": 81,
            "occluded_deviant": 91,
        },
    ) == "trial_sequence_identity"


def test_replay_gaze_and_escape_metadata_are_preserved(tmp_path):
    """Triggers 150/151/152 remain events and replay phase follows trigger 151."""
    csv_path = _write_debug_csv(
        tmp_path,
        runtime="v20",
        subject_id="04",
        block=2,
        run=3,
        rows=[
            (71, 1, 1, 1.0, "trial_condition_onset"),
            (40, 1, 12, 1.1, "trial_sequence_identity"),
            (150, 1, 100, 1.8, "gaze_break"),
            (151, 2, 0, 2.0, "replay_start"),
            (71, 2, 1, 2.1, "trial_condition_onset"),
            (40, 2, 12, 2.2, "trial_sequence_identity"),
            (152, 2, 0, 2.8, "escape_termination"),
        ],
    )
    segments = STEP02.load_reported_trigger_segments(
        [csv_path],
        subject_id="04",
        runtime_version="v20",
    )
    segment = segments[0]

    assert segment.termination_status == "esc_terminated"
    assert [event.replay_status for event in segment.events] == [
        "original",
        "original",
        "original",
        "replay_start",
        "replay",
        "replay",
        "replay",
    ]

    events = _resolved_events(
        [event.trigger for event in segment.events],
        source_fif="sub04_run03_220526.fif",
        acquisition_run=3,
        seq_max=40,
        condition_triggers={
            "always_visible": 71,
            "occluded_nondeviant": 81,
            "occluded_deviant": 91,
        },
        expected=segment.events,
    )
    events, _alignment, _summary = STEP02.align_reported_segments_to_events(
        events,
        segments,
    )

    by_code = {event["trigger_code"]: event for event in events}
    assert by_code[150]["is_gaze_break"] is True
    assert by_code[151]["replay_status"] == "replay_start"
    assert by_code[152]["is_escape_termination"] is True
    replay_condition = [
        event
        for event in events
        if event["trigger_code"] == 71 and event["behavioral_trial"] == 2
    ][0]
    assert replay_condition["is_replay_trial"] is True


def test_alignment_maps_multi_block_fif_and_partial_final_segment(tmp_path):
    """
    One acquisition FIF can contain multiple behavioral blocks.

    The second segment intentionally contains only one trial plus ESC. The
    parser/alignment uses those saved rows and never expects a prepared count.
    """
    csv_block1 = _write_debug_csv(
        tmp_path,
        runtime="v20",
        subject_id="04",
        block=1,
        run=3,
        rows=[
            (71, 1, 1, 1.0, "trial_condition_onset"),
            (40, 1, 12, 1.1, "trial_sequence_identity"),
            (111, 1, 114, 1.8, "occlusion_start"),
        ],
    )
    csv_partial = _write_debug_csv(
        tmp_path,
        runtime="v20",
        subject_id="04",
        block=2,
        run=3,
        rows=[
            (81, 1, 1, 2.0, "trial_condition_onset"),
            (39, 1, 12, 2.1, "trial_sequence_identity"),
            (152, 1, 0, 2.4, "escape_termination"),
        ],
    )
    segments = STEP02.load_reported_trigger_segments(
        [csv_block1, csv_partial],
        subject_id="04",
        runtime_version="v20",
    )
    expected = STEP02.flatten_reported_segments(segments)
    events = _resolved_events(
        [71, 40, 111, 81, 39, 152],
        source_fif="sub04_run03_220526.fif",
        acquisition_run=3,
        seq_max=40,
        condition_triggers={
            "always_visible": 71,
            "occluded_nondeviant": 81,
            "occluded_deviant": 91,
        },
        expected=expected,
    )
    events, alignment, summary = STEP02.align_reported_segments_to_events(
        events,
        segments,
    )
    associations = STEP02.build_source_associations(events, summary)

    assert [event["behavioral_block"] for event in events] == [1, 1, 1, 2, 2, 2]
    assert all(event["acquisition_run"] == 3 for event in events)
    assert all(
        event["raw_sample_frame"] == event["sample_start_raw"]
        for event in events
    )
    assert all(
        event["source_fif"] == "sub04_run03_220526.fif"
        for event in events
    )
    assert summary["status_counts"]["MATCH"] == 6
    assert not any(row["status"] == "REPORTED_ONLY" for row in alignment)
    source = associations["meg_acquisition_run_to_behavioral_segments"][0]
    assert source["source_fif"] == "sub04_run03_220526.fif"
    assert source["behavioral_segments"] == [
        {
            "behavioral_block": 1,
            "behavioral_run": 3,
            "aligned_event_count": 3,
        },
        {
            "behavioral_block": 2,
            "behavioral_run": 3,
            "aligned_event_count": 3,
        },
    ]
    assert segments[1].termination_status == "esc_terminated"


def test_source_associations_require_exact_trigger_matches(tmp_path):
    """A mismatch remains QC metadata but cannot establish a source mapping."""
    csv_path = _write_debug_csv(
        tmp_path,
        runtime="v18",
        subject_id="01",
        block=1,
        run=1,
        rows=[
            (31, 1, 1, 1.0, "trial_condition_onset"),
            (9, 1, 12, 1.1, "trial_sequence_identity"),
            (111, 1, 114, 1.8, "occlusion_start"),
        ],
    )
    segments = STEP02.load_reported_trigger_segments(
        [csv_path],
        subject_id="01",
        runtime_version="v18",
    )
    events = _resolved_events(
        [31, 8, 111],
        source_fif="SUB01_block01_meg.fif",
        acquisition_run=1,
        seq_max=20,
        condition_triggers={
            "always_visible": 31,
            "occluded_nondeviant": 41,
            "occluded_deviant": 51,
        },
        expected=segments[0].events,
    )
    events, _alignment, summary = STEP02.align_reported_segments_to_events(
        events,
        segments,
    )
    associations = STEP02.build_source_associations(events, summary)

    assert summary["status_counts"] == {"MATCH": 2, "MISMATCH": 1}
    source = associations["meg_acquisition_run_to_behavioral_segments"][0]
    assert source["aligned_event_count"] == 2
    assert source["behavioral_segments"] == [
        {
            "behavioral_block": 1,
            "behavioral_run": 1,
            "aligned_event_count": 2,
        }
    ]


def test_observed_only_escape_inherits_preceding_segment_context(tmp_path):
    """
    Preserve an STI101 ESC when the runtime did not save a matching CSV row.

    This models sub04 block02/run03: the saved CSV ends after the next trial
    onset, while trigger 152 remains visible in the same acquisition FIF.
    """
    csv_path = _write_debug_csv(
        tmp_path,
        runtime="v20",
        subject_id="04",
        block=2,
        run=3,
        rows=[
            (81, 18, 1, 1.0, "trial_condition_onset"),
            (150, 18, 5, 1.1, "gaze_break"),
            (81, 19, 1, 2.0, "trial_condition_onset"),
        ],
    )
    segments = STEP02.load_reported_trigger_segments(
        [csv_path],
        subject_id="04",
        runtime_version="v20",
    )
    events = _resolved_events(
        [81, 150, 81, 152],
        source_fif="sub04_run02_220526.fif",
        acquisition_run=2,
        seq_max=40,
        condition_triggers={
            "always_visible": 71,
            "occluded_nondeviant": 81,
            "occluded_deviant": 91,
        },
        expected=segments[0].events,
    )
    events, alignment, summary = STEP02.align_reported_segments_to_events(
        events,
        segments,
    )
    associations = STEP02.build_source_associations(events, summary)

    escape = events[-1]
    assert escape["trigger_code"] == 152
    assert escape["reported_alignment_status"] == "UNALIGNED_OBSERVED"
    assert escape["behavioral_block"] == 2
    assert escape["behavioral_run"] == 3
    assert (
        escape["behavioral_context_source"]
        == "inferred_esc_from_preceding_segment"
    )
    assert alignment[-1]["status"] == "OBSERVED_ONLY"
    assert alignment[-1]["behavioral_block"] == 2
    assert alignment[-1]["behavioral_run"] == 3
    assert summary["segments"][0]["termination_status"] == (
        "observed_esc_after_saved_rows"
    )
    assert summary["segments"][0]["observed_only_count"] == 1

    source = associations["meg_acquisition_run_to_behavioral_segments"][0]
    assert source["aligned_event_count"] == 3
