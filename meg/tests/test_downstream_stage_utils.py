"""
test_downstream_stage_utils.py

Purpose
Validate Chat-4 lineage, run-boundary, trial-exclusion, transition, replay,
identity-axis, trajectory-discovery, and guarded-path contracts without running
full participant epoching.

Usage
python -m pytest /Users/damiano/Documents/UniTn/Dynamo/Attention/DAD/meg/tests/test_downstream_stage_utils.py -q
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest


CODE_DIR = Path(__file__).resolve().parents[1] / "code" / "preprocessing_basic"
sys.path.insert(0, str(CODE_DIR))

import downstream_stage_utils as dsu  # noqa: E402
from preprocessing_stage_utils import AcquisitionRunRecord  # noqa: E402


def _row(
    code: int,
    sample: int,
    *,
    source: str = "sub04_run03_220526.fif",
    label: str | None = None,
    unresolved: bool = False,
    transition: bool = False,
    replay: bool = False,
) -> dict[str, str]:
    """Build one canonical-trigger-like row for trial-contract tests."""
    condition_labels = {
        31: "trial_condition_onset_always_visible",
        41: "trial_condition_onset_occluded_nondeviant",
        51: "trial_condition_onset_occluded_deviant",
        71: "trial_condition_onset_always_visible",
        81: "trial_condition_onset_occluded_nondeviant",
        91: "trial_condition_onset_occluded_deviant",
    }
    return {
        "trigger_code": str(code),
        "event_label": label or condition_labels.get(
            code,
            "trial_sequence_identity" if 1 <= code <= 40 else f"trigger_{code}",
        ),
        "sample_start_raw": str(sample),
        "sample_end_raw": str(sample + 1),
        "raw_sample_frame": str(sample),
        "source_fif": source,
        "source_fif_path": f"/raw/{source}",
        "acquisition_run": "3",
        "is_transition_warning": str(transition),
        "transition_partner_code": "72" if transition else "",
        "expected_trigger_code": str(code),
        "transition_decision_source": "reported" if transition else "",
        "transition_resolution_confidence": "1.0" if transition else "",
        "unresolved_transition": str(unresolved),
        "behavioral_block": "2",
        "behavioral_run": "3",
        "behavioral_trial": "18",
        "behavioral_context_source": "reported_alignment",
        "replay_status": "replay" if replay else "original",
        "is_replay_trial": str(replay),
        "is_gaze_break": str(code == 150),
        "is_escape_termination": str(code == 152),
    }


def _policy() -> dict:
    """Return the production trial-exclusion policy used by synthetic tests."""
    return {
        "gaze_break_codes": [150],
        "catch_type_1_codes": [116, 117],
        "catch_type_2_codes": [113, 118, 201, 202],
        "esc_termination_codes": [152],
        "exclude_gaze_break": True,
        "exclude_catch_type_1": True,
        "exclude_catch_type_2": True,
        "exclude_esc_termination": True,
        "exclude_replay_trials": False,
        "exclude_unresolved_transitions": True,
    }


def test_trajectory_discovery_requires_only_observed_mat(tmp_path):
    """Expected paths come from nondeviant trials in the observed stimulus MAT."""
    input_dir = tmp_path / "input_files"
    input_dir.mkdir()
    observed = input_dir / "MovDot_Sub04_V30_rescueTraject.mat"
    observed.touch()
    cfg = {
        "resolved_paths": {"subject_input_files_dir": str(input_dir)},
        "data_discovery": {
            "trajectory_mat_glob": "input_files/MovDot_Sub04_V*_rescueTraject.mat"
        },
        "subject_profile": {"versions": {"stimuli": "V30"}},
    }

    result = dsu.discover_trajectory_models(cfg)

    assert result == {
        "observed": str(observed.resolve()),
        "expected_source_condition": "occluded_nondeviant",
        "stimuli_version": "V30",
    }


@pytest.mark.parametrize(
    ("event_map", "identity"),
    [
        ({"always_visible": 31, "occluded_nondeviant": 41, "occluded_deviant": 51}, 20),
        ({"always_visible": 71, "occluded_nondeviant": 81, "occluded_deviant": 91}, 40),
    ],
)
def test_trial_builder_supports_both_trigger_generations(event_map, identity):
    """Condition and identity axes are profile-driven for 20/40 subjects."""
    condition = event_map["always_visible"]
    decisions = dsu.build_trial_decisions(
        [_row(condition, 100), _row(identity, 110)],
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=200,
        event_id_map=event_map,
        identity_range=range(1, identity + 1),
        exclusion_policy=_policy(),
    )
    assert decisions[0]["trial_sequence_identity"] == identity
    assert decisions[0]["included_by_trial_policy"] is True


def test_gaze_catch_types_and_esc_are_distinct_exclusion_reasons():
    """Trigger-150 and both catch families remain scientifically distinguishable."""
    rows = [
        _row(71, 100),
        _row(8, 110),
        _row(150, 120),
        _row(71, 200),
        _row(9, 210),
        _row(116, 220),
        _row(71, 300),
        _row(10, 310),
        _row(118, 320),
        _row(201, 330),
        _row(71, 400),
        _row(11, 410),
        _row(152, 420),
    ]
    decisions = dsu.build_trial_decisions(
        rows,
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=430,
        event_id_map={"always_visible": 71},
        identity_range=range(1, 41),
        exclusion_policy=_policy(),
    )
    assert [item["exclusion_reasons"] for item in decisions] == [
        "gaze_break",
        "catch_type_1",
        "catch_type_2",
        "esc_termination",
    ]


def test_excluded_gaze_trial_may_terminate_before_identity():
    """A gaze break before the identity frame excludes, rather than aborts, Stage 04."""
    decision = dsu.build_trial_decisions(
        [_row(31, 100), _row(150, 120)],
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=200,
        event_id_map={"always_visible": 31},
        identity_range=range(1, 21),
        exclusion_policy=_policy(),
    )[0]

    assert decision["trial_sequence_identity"] is None
    assert decision["identity_sample_start_raw"] is None
    assert decision["identity_mapping_status"] == "missing_excluded_before_identity"
    assert decision["exclusion_reasons"] == "gaze_break"
    assert decision["included_by_trial_policy"] is False


def test_v20_identity_31_is_selected_even_with_legacy_condition_label():
    """Numeric profile identity wins over a stale legacy label in canonical CSVs."""
    decision = dsu.build_trial_decisions(
        [
            _row(71, 100),
            _row(31, 110, label="trial_condition_onset_always_visible"),
        ],
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=200,
        event_id_map={"always_visible": 71},
        identity_range=range(1, 41),
        exclusion_policy=_policy(),
    )[0]

    assert decision["trial_sequence_identity"] == 31
    assert decision["identity_mapping_status"] == "resolved_numeric_profile_range"
    assert decision["included_by_trial_policy"] is True


def test_included_trial_without_identity_still_fails():
    """Only explicitly excluded partial trials may omit the identity trigger."""
    with pytest.raises(RuntimeError, match="Expected one resolved trial identity"):
        dsu.build_trial_decisions(
            [_row(71, 100)],
            source_fif="sub04_run03_220526.fif",
            run_end_sample_exclusive=200,
            event_id_map={"always_visible": 71},
            identity_range=range(1, 41),
            exclusion_policy=_policy(),
        )


def test_replay_is_included_and_resolved_transition_condition_is_retained():
    """Presented replay trials survive while resolved transition rows stay usable."""
    rows = [
        _row(81, 100, transition=True, replay=True),
        _row(40, 110, replay=True),
    ]
    decision = dsu.build_trial_decisions(
        rows,
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=200,
        event_id_map={"occluded_nondeviant": 81},
        identity_range=range(1, 41),
        exclusion_policy=_policy(),
    )[0]
    assert decision["is_replay_trial"] is True
    assert decision["is_transition_warning"] is True
    assert decision["included_by_trial_policy"] is True


def test_last_trial_interval_stops_at_run_boundary():
    """Evidence after the acquisition run cannot contaminate its last trial."""
    rows = [_row(71, 100), _row(8, 110), _row(150, 205)]
    decision = dsu.build_trial_decisions(
        rows,
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=200,
        event_id_map={"always_visible": 71},
        identity_range=range(1, 41),
        exclusion_policy=_policy(),
    )[0]
    assert decision["trial_interval_end_raw_exclusive"] == 200
    assert decision["gaze_break_evidence"] is False


def test_partial_esc_trial_is_excluded_without_expected_count_assumption():
    """Sub04-like partial termination uses saved events, not prepared totals."""
    rows = [_row(91, 100), _row(18, 110), _row(152, 150)]
    decisions = dsu.build_trial_decisions(
        rows,
        source_fif="sub04_run03_220526.fif",
        run_end_sample_exclusive=160,
        event_id_map={"occluded_deviant": 91},
        identity_range=range(1, 41),
        exclusion_policy=_policy(),
    )
    assert len(decisions) == 1
    assert decisions[0]["behavioral_trial"] == 18
    assert decisions[0]["exclusion_reasons"] == "esc_termination"


def test_identity_coverage_errors_or_retains_nan_rows_explicitly():
    """Missing identities never cause silent synthesis or truncation."""
    with pytest.raises(RuntimeError, match="no surviving epochs"):
        dsu.identity_coverage([1, 1, 2], range(1, 4), missing_policy="error")
    counts, missing = dsu.identity_coverage(
        [1, 1, 2],
        range(1, 4),
        missing_policy="nan",
    )
    assert counts == {1: 2, 2: 1, 3: 0}
    assert missing == [3]


def test_irregular_sub01_lineage_uses_source_stem_not_sorted_position(
    tmp_path,
    monkeypatch,
):
    """The run08/block11 source resolves only its exact `_ica.fif` output."""
    source = AcquisitionRunRecord(
        source_path="/raw/sub01_run08_block11_150426.fif",
        source_filename="sub01_run08_block11_150426.fif",
        source_stem="sub01_run08_block11_150426",
        acquisition_run_label=8,
        raw_relative_path="sub01_150426/meg/sub01_run08_block11_150426.fif",
    )
    sss = tmp_path / "derivatives" / "sss"
    sss.mkdir(parents=True)
    cleaned = sss / "sub01_run08_block11_150426_ica.fif"
    cleaned.touch()
    manifest = tmp_path / "source_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "derivative_lineage": {
                    "stages": {
                        "ica_manual": {
                            "processing_pass": "apply",
                            "status": "completed",
                            "qc": {"manual_exclude": []},
                            "records": [
                                {
                                    "source_stem": source.source_stem,
                                    "status": "completed",
                                    "output_path": (
                                        "/cluster/derivatives/sss/"
                                        f"{cleaned.name}"
                                    ),
                                }
                            ],
                        }
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(dsu, "discover_raw_run_records", lambda _cfg: [source])
    monkeypatch.setattr(dsu, "stage_output_dirs", lambda _cfg: {"sss": sss})
    runs, gate = dsu.resolve_manual_ica_apply_runs({}, manifest)
    assert runs[0].cleaned_path == str(cleaned)
    assert gate["manual_exclude"] == []


def test_manual_ica_gate_rejects_review_only_lineage(tmp_path, monkeypatch):
    """Downstream stages cannot bypass the manual-ICA apply prerequisite."""
    manifest = tmp_path / "source_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "derivative_lineage": {
                    "stages": {
                        "ica_manual": {
                            "processing_pass": "review",
                            "status": "review_complete",
                            "qc": {"manual_exclude": [1]},
                            "records": [],
                        }
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(dsu, "discover_raw_run_records", lambda _cfg: [])
    with pytest.raises(RuntimeError, match="processing_pass='apply'"):
        dsu.resolve_manual_ica_apply_runs({}, manifest)
