"""
test_trigger_parser.py

Purpose
Validate core transition-resolution behavior for
`code/preprocessing_basic/02_trigger_extraction_custom.py`.

This test covers the highest-risk trigger rule:
- expected-trigger-first symmetric handling of contiguous transition pairs.

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
