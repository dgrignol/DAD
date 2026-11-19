"""Test Stage 2 scientific accounting and protection of accepted inputs.

Usage: python -m pytest meg/tests/test_trial_preparation.py -q

All FIFs, Stage 1 tables, decisions, and outputs are synthetic temporary files.
Fixtures reproduce absolute sample offsets, P30's two prescribed exclusions,
P37's interrupted trial/unavailable tail, and the Stage-04 consumer contract.
No participant data, scheduler, or shared derivative directory is accessed.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path
from types import SimpleNamespace
import sys

import mne
import numpy as np
import pytest

CODE = Path(__file__).resolve().parents[1] / "code/preprocessing_basic"
sys.path.insert(0, str(CODE))
import trial_preparation as tp
from downstream_stage_utils import REQUIRED_TRIGGER_COLUMNS, build_trial_decisions


def make_case(tmp_path, subject="30"):
    """Create accepted-like evidence with nonzero raw first_samp coordinates."""
    raw_root = tmp_path / "raw"
    folder = raw_root / f"sub{subject}_date" / "meg"
    folder.mkdir(parents=True)
    root = tmp_path / "derivatives"
    reports = root / f"reports/sub-{subject}/ses-01"
    pre = root / f"preprocessing/sub-{subject}/ses-01"
    reports.mkdir(parents=True)
    pre.mkdir(parents=True)
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "PREPROCESSING_DECISIONS.md").write_text("Accepted on 2026-10-19 for Stage 2 trial preparation")
    source = "sub30_run02_210926.fif" if subject == "30" else "sub37_run01_051026.fif"
    acquisition = 2 if subject == "30" else 1
    first = 59000 if subject == "30" else 550000
    n_times = 12000 if subject == "30" else 13000
    raw = mne.io.RawArray(np.zeros((1, n_times)), mne.create_info(["STI101"], 1000, ["stim"]),
                          first_samp=first, verbose="ERROR")
    raw.save(folder / source, overwrite=True, verbose="ERROR")
    profile = {"condition_triggers": {"always_visible": 71, "occluded_deviant": 91},
               "identity_range": [1, 40], "iti_sec_range": [1.0, 1.5],
               "known_exceptions": [], "truncated_acquisition_runs": [1] if subject == "37" else []}
    policy = {"gaze_break_codes": [150], "catch_type_1_codes": [116, 117],
              "catch_type_2_codes": [113, 118, 201, 202], "esc_termination_codes": [152],
              "exclude_gaze_break": True, "exclude_catch_type_1": True,
              "exclude_catch_type_2": True, "exclude_esc_termination": True,
              "exclude_replay_trials": False, "exclude_unresolved_transitions": True}
    if subject == "30":
        codes = [(71, 60099, 1), (15, 60191, 1), (71, 64694, 2), (18, 64786, 2),
                 (152, 65667, 2), (71, 67500, 3), (3, 67592, 3)]
        exclusions = [{"source_fif": source, "raw_sample_frame": sample, "condition_code": 71,
                       "reason": "stage1_unreported_incomplete_sequence", "decision_reference": "accepted/P30"}
                      for sample in (60099, 64694)]
    else:
        codes = [(91, 556000, 2), (6, 556092, 2), (91, 562193, 3), (28, 562284, 3)]
        exclusions = [{"source_fif": source, "raw_sample_frame": 562193, "condition_code": 91,
                       "reason": "stage1_incomplete_recording_tail", "decision_reference": "accepted/P37",
                       "require_recording_tail_evidence": True}]
    cfg = {"subject": {"id": subject}, "subject_profile": profile,
           "layout": {"subject_bids_id": f"sub-{subject}", "session_bids_id": "ses-01"},
           "project": {"docs_root": str(docs), "params_root": str(tmp_path / "params")},
           "data_discovery": {"meg_fif_glob": "meg/*.fif"},
           "_config_sources": {"default": "fixture.yaml"},
           "resolved_paths": {"data_root": str(raw_root), "subject_data_dir": str(folder.parent),
                              "derivatives_root": str(root), "derivatives_reports_dir": str(reports),
                              "derivatives_preprocessing_dir": str(pre),
                              "derivatives_analysis_dir": str(root / f"analysis/sub-{subject}/ses-01")},
           "epoching": {"event_id_map": profile["condition_triggers"], "trial_exclusion": policy,
                        "tmin_sec": -0.5, "tmax_sec": 2.667, "accepted_trial_exclusions": exclusions,
                        "require_stage2_ledger": True, "baseline_sec": None, "detrend": 1,
                        "reject_by_annotation": True, "amplitude_reject": {}}}
    events, alignment, deep = [], [], []
    for index, (code, sample, trial) in enumerate(codes, 1):
        event = dict.fromkeys(REQUIRED_TRIGGER_COLUMNS, "")
        event.update(trigger_code=str(code), sample_start_raw=str(sample), sample_end_raw=str(sample + 1),
                     raw_sample_frame=str(sample), source_fif=source, source_fif_path=str(folder / source),
                     acquisition_run=str(acquisition), event_index=str(index), behavioral_block=str(acquisition),
                     behavioral_run="3" if subject == "37" else "1", behavioral_trial=str(trial),
                     reported_source_csv=f"/original/triggers_{subject}.csv", replay_status="original",
                     reported_alignment_status="MATCH", unresolved_transition="False")
        events.append(event)
        alignment.append(dict(event, alignment_index=str(index), observed_event_index=str(index),
                              reported_replay_status="original", status="MATCH"))
        deep.append({"align_row": str(index), "status": "MATCH", "stage1_adjudication": "",
                     "adjudication_partner_align_row": ""})
    if subject == "37":
        for trial, code in [(3, 111), (3, 112), (4, 91), (4, 10)]:
            index = str(len(alignment) + 1)
            alignment.append(dict(events[-1], alignment_index=index, observed_event_index="", status="REPORTED_ONLY",
                                  source_fif="", raw_sample_frame="", reported_trigger_code=str(code),
                                  behavioral_trial=str(trial), reported_replay_status="original",
                                  constrained_acquisition_run="1", reported_only_reason=tp.TAIL_REASON))
            deep.append({"align_row": index, "status": "REPORTED_ONLY",
                         "stage1_adjudication": "RESOLVED_RECORDING_TAIL_UNAVAILABLE",
                         "adjudication_partner_align_row": ""})
    association = {"source_fif": source, "acquisition_run": acquisition, "association_status": "aligned",
                   "behavioral_segments": [{"behavioral_block": acquisition, "behavioral_run": 1}]}
    manifest = {"subject_bids_id": f"sub-{subject}", "subject_profile": profile, "warnings": [],
                "sources": {"task_meg_acquisition_runs": [{"source_filename": source,
                            "size_bytes": (folder / source).stat().st_size}]},
                "trigger_derived_associations": {"meg_acquisition_run_to_behavioral_segments": [association]}}
    inventory = {"subject_bids_id": f"sub-{subject}", "subject_profile": profile, "status": "ok",
                 "validation": {}, "files": {"meg_fif": [{"name": source, "n_samples": n_times,
                 "sfreq_hz": 1000, "online_highpass_hz": 0, "online_lowpass_hz": 500}]}}
    paths = tp.evidence_paths(cfg)
    for key, value in {"manifest": manifest, "inventory": inventory, "pipeline_run": {"status": "completed"},
                       "associations": association}.items():
        paths[key].parent.mkdir(parents=True, exist_ok=True)
        paths[key].write_text(json.dumps(value))
    for key, rows in {"events": events, "alignment": alignment, "deep": deep, "warnings": []}.items():
        paths[key].parent.mkdir(parents=True, exist_ok=True)
        paths[key].write_text(tp.csv_text(rows, ("unresolved_transition",)))
    return cfg, events


def test_p30_exact_exclusions_and_all_stage1_inputs_preserved(tmp_path):
    """Both incomplete candidates remain visible; a later complete trial survives."""
    cfg, _ = make_case(tmp_path)
    before = tp.evidence_fingerprints(tp.evidence_paths(cfg))
    report = tp.prepare_subject(cfg)
    rows = tp.load_prepared_decisions(cfg)
    assert [r["raw_sample_frame"] for r in rows if r["included_by_trial_policy"]] == [67500]
    assert [r["accepted_exclusion_reason"] for r in rows[:2]] == ["stage1_unreported_incomplete_sequence"] * 2
    assert report["accepted_exclusions_matched"] == 2
    assert tp.evidence_fingerprints(tp.evidence_paths(cfg)) == before
    assert rows[0]["recording_first_sample_raw"] == "59000"
    assert len(tp.read_csv(tp.stage2_paths(cfg)["identities"])) == 80
    with pytest.raises(FileExistsError):
        tp.prepare_subject(cfg)


def test_accepted_exclusion_must_match_exact_source_sample_and_code(tmp_path):
    """An edited key stops preparation, even if a neighboring trial is valid."""
    cfg, _ = make_case(tmp_path)
    cfg["epoching"]["accepted_trial_exclusions"][0]["raw_sample_frame"] += 1
    with pytest.raises(RuntimeError, match="exclusion onset missing"):
        tp.prepare_subject(cfg)
    assert not tp.stage2_paths(cfg)["ledger"].exists()


def test_reviewed_partial_boundary_is_independent_of_general_esc_policy(tmp_path):
    """An accepted interrupted trial stays excluded if the generic ESC rule changes."""
    cfg, _ = make_case(tmp_path)
    rule = cfg["epoching"]["accepted_trial_exclusions"][1]
    rule["required_boundary_marker"] = {"code": 152, "raw_sample_frame": 65667}
    cfg["epoching"]["trial_exclusion"]["exclude_esc_termination"] = False
    tp.prepare_subject(cfg)
    assert tp.load_prepared_decisions(cfg)[1]["included_by_trial_policy"] is False
    rule["required_boundary_marker"]["raw_sample_frame"] += 1
    with pytest.raises(RuntimeError, match="boundary marker missing"):
        tp.prepare_subject(cfg, validate_only=True)


def test_p37_partial_trial_and_unrecorded_tail_are_distinct(tmp_path):
    """Only observed onsets become candidates, with a reviewed final exclusion."""
    cfg, _ = make_case(tmp_path, "37")
    report = tp.prepare_subject(cfg)
    rows = tp.load_prepared_decisions(cfg)
    assert len(rows) == 2
    assert rows[-1]["accepted_exclusion_reason"] == "stage1_incomplete_recording_tail"
    assert rows[-1]["included_by_trial_policy"] is False
    assert rows[-1]["recording_end_sample_raw_exclusive"] == "563000"
    assert report["n_unavailable_logged_events"] == 4
    unavailable = tp.read_csv(tp.stage2_paths(cfg)["unavailable"])
    assert [r["recorded_candidate_present"] for r in unavailable] == ["True", "False"]
    assert all(not r["raw_condition_onset_sample"] for r in unavailable)


def test_p37_reviewed_partial_requires_canonical_tail_evidence(tmp_path):
    """A generic missing trigger cannot justify the reviewed cutoff exclusion."""
    cfg, _ = make_case(tmp_path, "37")
    path = tp.evidence_paths(cfg)["alignment"]
    rows = [r for r in tp.read_csv(path) if r["status"] != "REPORTED_ONLY"]
    path.write_text(tp.csv_text(rows))
    with pytest.raises(RuntimeError, match="lost its unavailable-tail evidence"):
        tp.prepare_subject(cfg)


def test_ledger_loader_allows_later_lineage_but_rejects_changed_evidence(tmp_path):
    """Later orchestration metadata must not invalidate immutable trigger evidence."""
    cfg, _ = make_case(tmp_path)
    tp.prepare_subject(cfg)
    path = tp.evidence_paths(cfg)["manifest"]
    data = json.loads(path.read_text())
    data.update(derivative_lineage={"later_stage": "completed"}, generated_utc="later")
    path.write_text(json.dumps(data))
    tp.evidence_paths(cfg)["pipeline_run"].write_text('{"status": "running"}')
    assert len(tp.load_prepared_decisions(cfg)) == 3
    path = tp.evidence_paths(cfg)["events"]
    path.write_text(path.read_text() + "\n")
    with pytest.raises(RuntimeError, match="Stage 1 evidence changed"):
        tp.load_prepared_decisions(cfg)


def test_validate_only_writes_nothing_and_edited_ledger_is_rejected(tmp_path):
    """Dry checks do not publish; hand editing a published decision is detected."""
    cfg, _ = make_case(tmp_path)
    tp.prepare_subject(cfg, validate_only=True)
    assert not tp.stage2_paths(cfg)["report"].parent.exists()
    tp.prepare_subject(cfg)
    ledger = tp.stage2_paths(cfg)["ledger"]
    ledger.write_text(ledger.read_text() + "\n")
    with pytest.raises(RuntimeError, match="output changed"):
        tp.load_prepared_decisions(cfg)


def test_unresolved_onset_and_missing_identity_remain_reviewable(tmp_path):
    """Unresolved onsets cannot disappear or lend identities to earlier trials."""
    cfg, events = make_case(tmp_path)
    events[-2]["unresolved_transition"] = "True"
    events.pop()
    rows = build_trial_decisions(events, source_fif=events[0]["source_fif"], run_end_sample_exclusive=71000,
        event_id_map=cfg["epoching"]["event_id_map"], identity_range=range(1, 41),
        exclusion_policy=cfg["epoching"]["trial_exclusion"],
        accepted_exclusions=cfg["epoching"]["accepted_trial_exclusions"], strict=False)
    assert len(rows) == 3
    assert rows[-1]["included_by_trial_policy"] is False
    assert rows[-1]["exclusion_reasons"] == "unresolved_transition"


def test_production_epoching_never_sends_p30_exclusions_to_mne(tmp_path, monkeypatch):
    """Exercise the real Stage-04 consumer through its MNE event-array boundary."""
    cfg, _ = make_case(tmp_path)
    tp.prepare_subject(cfg)
    spec = importlib.util.spec_from_file_location("stage2_epoch_consumer_test", CODE / "04_epoching.py")
    epoching = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(epoching)
    source = cfg["epoching"]["accepted_trial_exclusions"][0]["source_fif"]
    raw_path = Path(cfg["resolved_paths"]["subject_data_dir"]) / "meg" / source
    runs = [SimpleNamespace(cleaned_path=str(raw_path), source=SimpleNamespace(source_filename=source))]
    monkeypatch.setattr(epoching, "resolve_manual_ica_apply_runs", lambda *args: (runs, {"status": "passed"}))

    class ReachedEpochBoundary(Exception):
        """Stop after checking the actual prospective MNE events."""

    def check_events(raw, events, **kwargs):
        assert events[:, 0].tolist() == [67500]
        assert kwargs["metadata"]["raw_sample_frame"].tolist() == [67500]
        raise ReachedEpochBoundary

    monkeypatch.setattr(epoching.mne, "Epochs", check_events)
    with pytest.raises(ReachedEpochBoundary):
        epoching.run_epoching(cfg, source_manifest=tp.evidence_paths(cfg)["manifest"], overwrite=False)


def test_symlink_cannot_redirect_explicit_overwrite_to_stage1(tmp_path):
    """Output overwrite permission never extends to a linked Stage 1 table."""
    cfg, _ = make_case(tmp_path)
    path = tp.stage2_paths(cfg)["ledger"]
    path.parent.mkdir(parents=True)
    path.symlink_to(tp.evidence_paths(cfg)["events"])
    before = tp.evidence_fingerprints(tp.evidence_paths(cfg))
    with pytest.raises(ValueError, match="symlink"):
        tp.prepare_subject(cfg, overwrite=True)
    assert tp.evidence_fingerprints(tp.evidence_paths(cfg)) == before


def test_pre_warmup_esc_creates_no_trial_and_does_not_exclude_the_next_one(tmp_path):
    """P29 B6's unpaired early marker stays outside trial evidence."""
    cfg, events = make_case(tmp_path)
    rows = []
    for code, sample in [(152, 44761), (153, 78863), (71, 80000), (5, 80092)]:
        rows.append(dict(events[0], trigger_code=str(code), sample_start_raw=str(sample),
                         source_fif="sub29_run06_150926.fif", behavioral_block="6"))
    decisions = build_trial_decisions(rows, source_fif="sub29_run06_150926.fif",
        run_end_sample_exclusive=90000, event_id_map=cfg["epoching"]["event_id_map"],
        identity_range=range(1, 41), exclusion_policy=cfg["epoching"]["trial_exclusion"])
    assert len(decisions) == 1
    assert decisions[0]["included_by_trial_policy"] is True
    assert decisions[0]["evidence_trigger_codes"] == "71,5"


def test_behavioral_run_boundary_cannot_lend_gaze_evidence(tmp_path):
    """A gaze marker in a later run cannot reject the previous run's trial."""
    cfg, events = make_case(tmp_path)
    rows = [dict(events[0], trigger_code=str(code), sample_start_raw=str(sample), behavioral_run=str(run))
            for code, sample, run in [(71, 100, 1), (5, 110, 1), (150, 180, 2), (71, 200, 2), (6, 210, 2)]]
    decisions = build_trial_decisions(rows, source_fif=events[0]["source_fif"], run_end_sample_exclusive=300,
        event_id_map=cfg["epoching"]["event_id_map"], identity_range=range(1, 41),
        exclusion_policy=cfg["epoching"]["trial_exclusion"])
    assert decisions[0]["trial_interval_end_raw_exclusive"] == 180
    assert all(r["included_by_trial_policy"] for r in decisions)


def test_cohort_summary_requires_matching_complete_subject_bundles(tmp_path, monkeypatch):
    """Aggregation checks per-subject policies, counts, and every output hash."""
    cfg30, _ = make_case(tmp_path / "case30")
    cfg37, _ = make_case(tmp_path / "case37", "37")
    tp.prepare_subject(cfg30)
    tp.prepare_subject(cfg37)
    root30 = Path(cfg30["resolved_paths"]["derivatives_root"])
    root37 = Path(cfg37["resolved_paths"]["derivatives_root"])
    for branch in ("reports", "preprocessing"):
        shutil.copytree(root37 / branch / "sub-37", root30 / branch / "sub-37")
    monkeypatch.setattr(tp, "ANALYSIS_SUBJECT_IDS", ("30", "37"))
    selectors = {"30": cfg30, "37": cfg37}
    monkeypatch.setattr(tp, "load_and_resolve_config", lambda config, selector, **kwargs:
                        (selectors[Path(selector).stem[-2:]], [], None, None))
    report = tp.summarize_cohort(cfg30)
    assert report["n_subjects"] == 2
    assert report["eligibility_counts"] == {"eligible_under_current_policy": 2, "excluded": 3}
    group = root30 / "reports/group/trial_preparation/stage2_cohort_summary.json"
    before = group.read_bytes()
    (root30 / "reports/sub-37/ses-01/trial_preparation/stage2_trial_counts.csv").unlink()
    with pytest.raises(FileNotFoundError):
        tp.summarize_cohort(cfg30, overwrite=True)
    assert group.read_bytes() == before
