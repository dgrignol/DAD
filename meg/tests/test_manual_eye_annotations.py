"""Synthetic regression checks for the interactive ocular-review workflow.

Usage: python -B -m pytest -p no:cacheprovider meg/tests/test_manual_eye_annotations.py -q

No cohort data are opened. Temporary recordings and mocked GUI windows exercise
component numbering, bilateral picks, BACK with edited annotations, strict save
actions, independent NONE/redetection, legacy blink fingerprints, two-category
reviewed-empty decisions, failed-save recovery and tuning persistence. The real ICA projection
is checked separately in test_stage4_continuous.py.
"""

from datetime import datetime, timezone
import csv
from copy import deepcopy
import json
from pathlib import Path
import sys

import mne
import numpy as np
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code/preprocessing_basic"))
import manual_eye_annotations as review
from continuous_cleaning import file_hash
from preprocessing_stage_utils import AcquisitionRunRecord


# -----------------------------------------------------------------------------
# Two small synthetic detectors isolate interactive state/persistence behavior
# -----------------------------------------------------------------------------

def setup_dual_aid(tmp_path, monkeypatch, *, blink_components=(0,), eye_components=(1,)):
    """Reuse the native CLI fixture; replace only detector work and GUI selection.

    Both groups propose the same interval intentionally. The detectors count
    calls and move their proposal when the threshold changes, letting tests
    detect unintended regeneration of the other group. All files are temporary.
    """
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch)
    cfg = review.load_and_resolve_config(None)[0]
    cache = review.load_cached_aid(cfg)
    selection = {"blink": list(blink_components), "eye": list(eye_components)}
    calls = {"blink": [], "eye": []}

    def select(display, cfg, *, state):
        state["components"] = list(selection["blink"])
        state["eye_movement"]["components"] = list(selection["eye"])
        return {"model_sha256": cache["metadata"]["model_sha256"],
                "selected_components": list(selection["blink"]),
                "selected_eye_movement_components": list(selection["eye"])}

    def detector(kind):
        def run(sources, components, settings):
            assert components  # NONE must bypass the detector entirely.
            assert sources is cache["native_sources"]
            calls[kind].append(list(components))
            start = 500 if settings["threshold_std"] == 3. else 400
            return [[start, start + 50]], np.full(raw.n_times, settings["threshold_std"]), {
                "candidate_count": 1, "status": "proposals_for_review", "peaks": [],
                "components": list(components), "settings": dict(settings)}
        return run

    monkeypatch.setattr(review, "load_cached_aid", lambda *args: cache)
    monkeypatch.setattr(review, "add_preliminary_sources", select)
    monkeypatch.setattr(review, "detect_blinks", detector("blink"))
    monkeypatch.setattr(review, "detect_eye_movements", detector("eye"))
    return {"raw": raw, "output": output, "tuning": tuning, "args": args,
            "cfg": cfg, "cache": cache, "selection": selection, "calls": calls}


def saved_dual_review(session):
    """Read this temporary source's normal tuning receipt after explicit SAVE."""
    return yaml.safe_load(session["tuning"].read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]


@pytest.mark.parametrize("blink_components,eye_components", [([0], [1]), ([0], [0]),
    ([0], []), ([], [1]), ([], [])])
def test_each_none_group_is_disabled_and_overlap_stays_separate(tmp_path, monkeypatch, blink_components, eye_components):
    session = setup_dual_aid(tmp_path, monkeypatch, blink_components=blink_components, eye_components=eye_components)

    def inspect(display, *unused):
        assert not session["output"].exists() and not session["tuning"].exists()
        for components, label, score in ((blink_components, review.AUTO_LABEL, "BLINK_score"),
                (eye_components, review.EYE_AUTO_LABEL, "EYE_MOVEMENT_score")):
            assert (label in display.annotations.description) == bool(components)
            assert (score in display.ch_names) == bool(components)

    monkeypatch.setattr(review, "plot_review_browser", inspect)
    answers = iter(["", "SAFE", "SAVE", "MISC001", "Both selected groups and complete source reviewed"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    receipt = saved_dual_review(session)
    for kind, key, components in (("blink", "blink_auto", blink_components), ("eye", "eye_movement_auto", eye_components)):
        assert session["calls"][kind] == ([components] if components else [])
        assert receipt[key]["selected_components"] == components
        assert receipt[key]["detection_status"] == ("proposals_for_review" if components else "no_component_selected")
        assert receipt[key]["status"] == "reviewed"  # Only whole-source SAVE establishes this.
    if blink_components and eye_components:
        assert receipt["blink_auto"]["provisional_intervals_raw"] == receipt["eye_movement_auto"]["provisional_intervals_raw"]
        assert receipt["blink_auto"]["fingerprint"] != receipt["eye_movement_auto"]["fingerprint"]


@pytest.mark.parametrize("chosen", ["BLINK", "EYE"])
def test_redetect_changes_only_chosen_group_and_keeps_other_edits(tmp_path, monkeypatch, chosen):
    session = setup_dual_aid(tmp_path, monkeypatch)
    other_label = review.EYE_AUTO_LABEL if chosen == "BLINK" else review.AUTO_LABEL
    other_score = "EYE_MOVEMENT_score" if chosen == "BLINK" else "BLINK_score"
    windows = []

    def inspect(display, *unused):
        if not windows:
            display.annotations.duration[list(display.annotations.description).index(other_label)] = .7
            display.annotations.append(display.first_time + .2, .1, "BAD_segment")
            windows.append(([review.annotation_key(display, a) for a in display.annotations if a["description"] == other_label],
                            display.get_data(picks=other_score).copy()))
        else:
            assert [review.annotation_key(display, a) for a in display.annotations if a["description"] == other_label] == windows[0][0]
            np.testing.assert_array_equal(display.get_data(picks=other_score), windows[0][1])
            assert "BAD_segment" in display.annotations.description
            windows.append(True)

    monkeypatch.setattr(review, "plot_review_browser", inspect)
    answers = iter(["REDETECT", "", "both", chosen, "nan", "0", "3", "SAVE", "", "Reviewed both categories"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    assert len(windows) == 2
    assert len(session["calls"]["blink"]) == (2 if chosen == "BLINK" else 1)
    assert len(session["calls"]["eye"]) == (2 if chosen == "EYE" else 1)
    receipt = saved_dual_review(session)
    assert receipt["blink_auto"]["settings"]["threshold_std"] == (3. if chosen == "BLINK" else 5.)
    assert receipt["eye_movement_auto"]["settings"]["threshold_std"] == (3. if chosen == "EYE" else 4.)


@pytest.mark.parametrize("regenerate", ["YES", "NO"])
def test_back_reselects_eye_only_and_keeps_actual_provenance_when_declined(tmp_path, monkeypatch, regenerate):
    session = setup_dual_aid(tmp_path, monkeypatch)
    visits = []

    def inspect(display, *unused):
        visits.append(True)
        session["selection"]["eye"] = [0]

    monkeypatch.setattr(review, "plot_review_browser", inspect)
    answers = iter(["BACK", "bad", regenerate, "SAVE", "", "Reselected eye sources after inspection"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    assert len(visits) == 2 and session["calls"]["blink"] == [[0]]
    assert session["calls"]["eye"] == ([[1], [0]] if regenerate == "YES" else [[1]])
    eye = saved_dual_review(session)["eye_movement_auto"]
    assert eye["selected_components"] == [0]
    assert eye["components"] == ([0] if regenerate == "YES" else [1])
    assert eye["selection_applied"] == (regenerate == "YES")


def test_legacy_p04_blink_fingerprint_and_reviewed_empty_are_reused_without_migration(tmp_path, monkeypatch):
    session = setup_dual_aid(tmp_path, monkeypatch, eye_components=[])
    model_sha = "98d5306c5cd561be3c0f06aec548a6b8a975911e461e30e3b00d23bebb599082"
    session["cache"]["metadata"]["model_sha256"] = model_sha
    settings = review.blink_settings(session["cfg"])
    # These are P04's historical model/settings digest, not recording samples.
    fingerprint = "28cb7b7b17fda4edb08ebb8595386c5fae1d34f0d9eaadefceb65dd0f1923694"
    assert review.digest({"model": model_sha, "components": [0], "settings": settings}) == fingerprint
    session["output"].parent.mkdir(parents=True)
    mne.Annotations([], [], [], orig_time=session["raw"].annotations.orig_time).save(session["output"])
    previous = {"status": "reviewed", "annotation_count": 0, "edf_candidates_fingerprint": "current",
        "blink_auto": {"status": "reviewed", "fingerprint": fingerprint, "components": [0],
            "model_sha256": model_sha, "settings": settings, "retained_auto_count": 0,
            "removed_spans_raw": [[session["raw"].first_samp + 400, session["raw"].first_samp + 450]]}}
    session["cfg"]["artifact_annotation"]["manual_eye_reviews"] = {"sub99_run01_test.fif": previous}
    session["tuning"].write_text(yaml.safe_dump({"artifact_annotation": {"manual_eye_reviews": {"sub99_run01_test.fif": previous}}}))
    before_cfg = deepcopy(session["cfg"])
    before = session["output"].read_bytes(), session["tuning"].read_bytes()
    monkeypatch.setattr(review, "plot_review_browser", lambda display, *unused: (
        pytest.fail("Reviewed-empty proposals resurrected") if len(display.annotations) else None))
    answers = iter(["", "discard", "DISCARDALL"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica", "--overwrite"]) == 0
    assert session["calls"]["eye"] == [] and session["cfg"] == before_cfg
    assert (session["output"].read_bytes(), session["tuning"].read_bytes()) == before
    assert "eye_movement_auto" not in previous


def test_back_none_clears_only_untouched_eye_proposals_without_running_detector(tmp_path, monkeypatch):
    session = setup_dual_aid(tmp_path, monkeypatch)
    visits = []

    def inspect(display, *unused):
        visits.append(True)
        if len(visits) == 1:
            session["selection"]["eye"] = []
        else:
            assert review.AUTO_LABEL in display.annotations.description
            assert review.EYE_AUTO_LABEL not in display.annotations.description
            assert "BLINK_score" in display.ch_names and "EYE_MOVEMENT_score" not in display.ch_names

    monkeypatch.setattr(review, "plot_review_browser", inspect)
    answers = iter(["BACK", "SAVE", "", "Eye aid disabled after inspection"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    assert session["calls"] == {"blink": [[0]], "eye": [[1]]}
    eye = saved_dual_review(session)["eye_movement_auto"]
    assert eye["selection_status"] == "no_component_selected" and eye["retained_auto_count"] == 0


def test_manual_only_deletion_preserves_detector_fingerprint_and_remembers_veto(tmp_path, monkeypatch):
    session = setup_dual_aid(tmp_path, monkeypatch)
    monkeypatch.setattr(review, "plot_review_browser", lambda *a: None)
    answers = iter(["SAVE", "", "Retained both automatic categories"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    previous = saved_dual_review(session)
    session["cfg"]["artifact_annotation"]["manual_eye_reviews"] = {"sub99_run01_test.fif": previous}

    def edit(display, *unused):
        display.annotations.delete([i for i, label in enumerate(display.annotations.description) if label == review.AUTO_LABEL])

    monkeypatch.setattr(review, "plot_review_browser", edit)
    answers = iter(["SAVE", "", "Deleted blink proposal in manual-only review"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--overwrite"]) == 0
    saved = saved_dual_review(session)
    assert saved["blink_auto"]["fingerprint"] == previous["blink_auto"]["fingerprint"]
    assert saved["blink_auto"]["removed_spans_raw"] == previous["blink_auto"]["provisional_intervals_raw"]
    assert saved["blink_auto"]["retained_auto_count"] == 0
    assert saved["eye_movement_auto"] == previous["eye_movement_auto"]


@pytest.mark.parametrize("chosen", ["BLINK", "EYE"])
def test_redetect_none_skips_only_that_detector_without_threshold_prompt(tmp_path, monkeypatch, chosen):
    session = setup_dual_aid(tmp_path, monkeypatch,
        blink_components=[] if chosen == "BLINK" else [0], eye_components=[] if chosen == "EYE" else [1])
    monkeypatch.setattr(review, "plot_review_browser", lambda *a: None)
    answers = iter(["REDETECT", chosen, "SAVE", "", "Reviewed disabled detector group manually"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    assert session["calls"] == {"blink": [] if chosen == "BLINK" else [[0]],
                                "eye": [] if chosen == "EYE" else [[1]]}


def test_none_after_saved_review_retains_intervals_and_previous_detector_provenance(tmp_path, monkeypatch):
    session = setup_dual_aid(tmp_path, monkeypatch)
    monkeypatch.setattr(review, "plot_review_browser", lambda *a: None)
    answers = iter(["SAVE", "", "Retained both categories"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    previous = saved_dual_review(session)
    session["cfg"]["artifact_annotation"]["manual_eye_reviews"] = {"sub99_run01_test.fif": previous}
    session["selection"]["eye"] = []
    answers = iter(["SAVE", "", "Disabled new eye proposals; retained previously reviewed interval"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica", "--overwrite"]) == 0
    eye = saved_dual_review(session)["eye_movement_auto"]
    assert eye["retained_auto_count"] == 1 and eye["selection_status"] == "no_component_selected"
    assert eye["prior_detections"][0]["fingerprint"] == previous["eye_movement_auto"]["fingerprint"]
    assert eye["prior_detections"][0]["components"] == [1]
    assert session["calls"]["eye"] == [[1]]


def test_reviewed_empty_for_both_groups_survives_reopening_and_manual_only_save(tmp_path, monkeypatch):
    session = setup_dual_aid(tmp_path, monkeypatch)
    monkeypatch.setattr(review, "plot_review_browser", lambda display, *unused: display.annotations.delete(list(range(len(display.annotations)))))
    answers = iter(["SAVE", "", "Reviewed and deleted every candidate"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--preliminary-ica"]) == 0
    previous = saved_dual_review(session)
    for key in ("blink_auto", "eye_movement_auto"):
        assert previous[key]["retained_auto_count"] == 0 and previous[key]["removed_spans_raw"]
    session["cfg"]["artifact_annotation"]["manual_eye_reviews"] = {"sub99_run01_test.fif": previous}
    before = session["output"].read_bytes(), session["tuning"].read_bytes()
    monkeypatch.setattr(review, "plot_review_browser", lambda display, *unused: pytest.fail("reseeded") if len(display.annotations) else None)
    monkeypatch.setattr("builtins.input", lambda prompt: "DISCARDALL")
    assert review.main(session["args"] + ["--preliminary-ica", "--overwrite"]) == 0
    assert (session["output"].read_bytes(), session["tuning"].read_bytes()) == before
    # A subsequent manual-only SAVE retains both groups' detector provenance.
    answers = iter(["SAVE", "", "Reviewed empty sets without ICA display"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(session["args"] + ["--overwrite"]) == 0
    saved = saved_dual_review(session)
    assert saved["blink_auto"] == previous["blink_auto"]
    assert saved["eye_movement_auto"] == previous["eye_movement_auto"]


def test_maps_and_sources_precede_blink_then_eye_prompt_and_display_shared_component(monkeypatch):
    from types import SimpleNamespace
    steps = []
    sources = mne.io.RawArray(np.ones((2, 1000)), mne.create_info(["ICA000", "ICA001"], 100., "misc"),
                             first_samp=13000, verbose="ERROR")
    display = mne.io.RawArray(np.zeros((1, 1000)), mne.create_info(["MISC001"], 100., "misc"),
                             first_samp=13000, verbose="ERROR")
    state = {"model": SimpleNamespace(n_components_=2, plot_components=lambda **kw: steps.append("maps")),
        "native_sources": sources, "metadata": {"identity": {"channels": ["MEG0121"]},
            "model_sha256": "model", "fingerprint": "cache", "converged": True, "fit_summary": {"fit_duration_sec": 10.}}}
    monkeypatch.setattr(mne.io.BaseRaw, "plot", lambda *a, **k: steps.append("sources"))

    def answer(prompt):
        assert steps[:2] == ["maps", "sources"]
        steps.append("eye" if "eye-movement" in prompt else "blink")
        return "ICA000"

    monkeypatch.setattr("builtins.input", answer)
    receipt = review.add_preliminary_sources(display, {}, state=state)
    assert steps == ["maps", "sources", "blink", "eye"]
    assert display.ch_names == ["MISC001", "AID000"]  # Shared selection displayed once.
    assert receipt["selected_components"] == receipt["selected_eye_movement_components"] == [0]
    assert display.first_samp == sources.first_samp and display.n_times == sources.n_times


# -----------------------------------------------------------------------------
# Terminal choices must match displayed component identities and explicit actions
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("text, expected", [("0", [0]), ("ICA000", [0]), ("1", [1]),
    ("ICA001,0,ICA001", [0, 1]), ("NONE", [])])
def test_component_names_and_numbers_are_zero_based(monkeypatch, text, expected):
    monkeypatch.setattr("builtins.input", lambda prompt: text)
    assert review.choose_aid_components(3) == expected


def test_invalid_component_input_retries(monkeypatch):
    answers = iter(["", "none", "NONE,0", "-1", "ICA100", "eye", "0,", "ICA000"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.choose_aid_components(3) == [0]


@pytest.mark.parametrize("action", ["SAVE", "DISCARDALL", "BACK", "REDETECT"])
def test_blank_and_typo_actions_never_discard(monkeypatch, action):
    answers = iter(["", "SAFE", "save", action])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.review_action(preliminary_ica=True, remaining=0) == action


def test_back_requires_an_existing_preliminary_fit(monkeypatch):
    answers = iter(["BACK", "DISCARDALL"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.review_action(preliminary_ica=False, remaining=2) == "DISCARDALL"


def test_frontal_regions_include_both_sides_and_no_posterior_sensors():
    names = ["MEG0121", "MEG0312", "MEG0811", "MEG0913", "MEG2111", "MISC001"]
    raw = mne.io.RawArray(np.zeros((6, 10)), mne.create_info(names, 100.,
        ["mag", "grad", "mag", "grad", "mag", "misc"]), verbose="ERROR")
    groups = review.frontal_channel_groups(raw)
    assert set(groups["Left-frontal"]) == {"MEG0121", "MEG0312"}
    assert set(groups["Right-frontal"]) == {"MEG0811", "MEG0913"}


# -----------------------------------------------------------------------------
# Exercise the CLI on an exact temporary source, replacing only external/GUI work
# -----------------------------------------------------------------------------

def setup_session(tmp_path, monkeypatch, *, recover=False):
    """Create a native-coordinate source and current non-BAD EDF proposals."""
    dirs = {key: tmp_path / key for key in ("sss", "reports")}
    for path in dirs.values():
        path.mkdir()
    cfg = {"subject": {"id": "99"}, "ica": {}, "artifact_annotation": {
        "edf_blink_candidates": {"enabled": True}}, "resolved_paths": {
        "data_root": str(tmp_path / "raw"), "derivatives_preprocessing_dir": str(tmp_path / "preprocessing")}}
    record = AcquisitionRunRecord("source", "sub99_run01_test.fif", "sub99_run01_test", 1, "source")
    raw = mne.io.RawArray(np.zeros((3, 800)), mne.create_info(
        ["MEG0121", "MEG0811", "MISC001"], 100., ["mag", "mag", "misc"]),
        first_samp=1000, verbose="ERROR")
    raw.set_meas_date(datetime(2026, 1, 1, tzinfo=timezone.utc))
    raw.save(dirs["sss"] / "sub99_run01_test_raw_sss.fif", overwrite=True, verbose="ERROR")
    output = review.manual_path(cfg, record)
    if recover:
        output.parent.mkdir(parents=True)
        # A save failure may leave an explicitly empty annotation file but no receipt.
        mne.Annotations([], [], [], orig_time=raw.annotations.orig_time).save(output)
    proposal = {"fingerprint": "current", "status": "candidates_for_manual_review",
                "candidate_count": 2, "long_candidate_count": 1}
    candidates = mne.Annotations([raw.first_time + 1., raw.first_time + 3.], [.1, .7],
        ["EDF_blink_candidate", "EDF_blink_long_candidate"], orig_time=raw.annotations.orig_time)
    monkeypatch.setattr(review, "load_and_resolve_config", lambda *a, **k: (cfg,))
    monkeypatch.setattr(review, "stage_output_dirs", lambda cfg: dirs)
    monkeypatch.setattr(review, "require_review", lambda *a, **k: None)
    monkeypatch.setattr(review, "discover_raw_run_records", lambda cfg: [record])
    monkeypatch.setattr(review, "load_candidates", lambda *a: (proposal, candidates))
    monkeypatch.setattr(review, "add_edf_display", lambda *a: {"display_gain": 1.})
    monkeypatch.setattr(review, "load_timing_preview", lambda *a: None)
    sources = mne.io.RawArray(np.zeros((2, raw.n_times)), mne.create_info(["ICA000", "ICA001"], 100., "misc"),
                             first_samp=raw.first_samp, verbose="ERROR")
    sources.set_meas_date(raw.info["meas_date"])
    monkeypatch.setattr(review, "load_cached_aid", lambda *a: {"model": None, "native_sources": sources,
        "metadata": {"model_sha256": "synthetic_model", "fingerprint": "synthetic_cache"}})
    monkeypatch.setattr(mne.viz, "set_browser_backend", lambda *a: None)
    tuning = tmp_path / "tuning.yaml"
    args = ["--config", "config", "--subject-config", "subject", "--tuning-config", str(tuning),
            "--source", record.source_filename, "--basis", "synthetic_review"]
    if recover:
        args += ["--overwrite", "--recover-annotations"]
    return raw, output, tuning, args


def test_targeted_triage_save_uses_exact_intervals_without_edf_and_keeps_deletions(tmp_path, monkeypatch):
    """Import a counted mask, delete one auto interval, SAVE and reopen without reseeding."""
    import annotation_triage
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch)
    proposals = {"payload_sha256": "counted-mask", "annotations": [
        {"onset": raw.first_time + 1., "duration": .2, "description": "BAD_blink_auto"},
        {"onset": raw.first_time + 3., "duration": .2, "description": "BAD_eye_movement_auto"}],
        "detectors": {"blink": {"model_sha256": "model", "selected_components": [0],
            "fingerprint": "blink-detector", "settings": {}, "removed_spans_raw": []},
            "eye_movement": {"model_sha256": "model", "selected_components": [1],
            "fingerprint": "eye-detector", "settings": {}, "removed_spans_raw": []}}}
    monkeypatch.setattr(annotation_triage, "load_review_mask", lambda *a: proposals)
    monkeypatch.setattr(review, "load_candidates", lambda *a: pytest.fail("Triage review must not load EDF candidates"))
    monkeypatch.setattr(review, "add_edf_display", lambda *a: pytest.fail("Triage review must not load EDF traces"))
    def browser(display, *unused):
        assert list(display.annotations.description) == ["BAD_blink_auto", "BAD_eye_movement_auto"]
        display.annotations.delete(0)
    monkeypatch.setattr(review, "plot_review_browser", browser)
    answers = iter(["SAVE", "MISC001", "Whole source reviewed from exact triage mask"])
    monkeypatch.setattr("builtins.input", lambda *a: next(answers))
    assert review.main(args + ["--triage-proposals", "mask.json"]) == 0
    assert list(mne.read_annotations(output).description) == ["BAD_eye_movement_auto"]
    receipt = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert receipt["edf_guidance_policy"] == "not_used_annotation_triage"
    assert receipt["triage_origin"]["payload_sha256"] == "counted-mask"
    assert receipt["blink_auto"]["removed_spans_raw"] == [[1100, 1120]]
    assert receipt["eye_movement_auto"]["retained_auto_count"] == 1
    cfg = review.load_and_resolve_config(None)[0]
    cfg["artifact_annotation"]["manual_eye_reviews"] = {"sub99_run01_test.fif": receipt}
    monkeypatch.setattr(review, "plot_review_browser", lambda display, *a:
        None if list(display.annotations.description) == ["BAD_eye_movement_auto"] else pytest.fail("Deleted interval was reseeded"))
    answers = iter(["SAVE", "MISC001", "Reopened reviewed intervals"])
    monkeypatch.setattr("builtins.input", lambda *a: next(answers))
    assert review.main(args + ["--overwrite"]) == 0
    again = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert again["triage_origin"] == receipt["triage_origin"]
    assert again["edf_guidance_policy"] == receipt["edf_guidance_policy"]


def test_back_preserves_candidate_deletions_and_new_annotations_before_save(tmp_path, monkeypatch):
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch)
    windows = []
    selections = []

    def browser(display, *unused):
        windows.append(1)
        if len(windows) == 1:
            display.annotations.append(display.first_time + .1, .1, "BAD_sensor_fault")
            display.annotations.delete([i for i, label in enumerate(display.annotations.description) if label.startswith("EDF_")])
            display.annotations.append(display.first_time + 5., .2, "BAD_blink")
        else:
            assert set(display.annotations.description) == {"BAD_sensor_fault", "BAD_blink"}

    def select(display, cfg, *, state):
        selections.append(state)
        state["components"] = [len(selections) - 1]
        return {"purpose": "supplementary_annotation_only", "selected_components": [len(selections) - 1],
                "fit_duration_sec": np.float64(800.79), "model_sha256": "synthetic_model"}

    monkeypatch.setattr(review, "plot_review_browser", browser)
    monkeypatch.setattr(review, "add_preliminary_sources", select)
    answers = iter(["SAFE", "BACK", "YES", "SAVE", "MISC001", "Meaning checked after inspection"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args + ["--preliminary-ica"]) == 0
    assert len(windows) == 2 and len(selections) == 2 and selections[0] is selections[1]
    ann = mne.read_annotations(output)
    assert set(ann.description) == {"BAD_sensor_fault", "BAD_blink"}
    assert ann.orig_time == raw.annotations.orig_time
    receipt = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert receipt["status"] == "reviewed"
    assert receipt["preliminary_ica"]["fit_duration_sec"] == 800.79
    assert receipt["annotation_sha256"] == file_hash(output)
    assert receipt["edf_candidates_fingerprint"] == "current"


def test_recovery_keeps_deleted_candidates_and_discard_leaves_files_untouched(tmp_path, monkeypatch):
    _, output, tuning, args = setup_session(tmp_path, monkeypatch, recover=True)
    before = output.read_bytes()
    tuning.write_text("stage4_review:\n  03_A_MaxFilt:\n    status: accepted\n")
    before_tuning = tuning.read_bytes()

    def browser(display, *unused):
        assert len(display.annotations) == 0  # Do not resurrect the deleted proposals.
        display.annotations.append(display.first_time + 1., .1, "BAD_blink")

    monkeypatch.setattr(review, "plot_review_browser", browser)
    answers = iter(["", "SAFE", "DISCARDALL"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args) == 0
    assert output.read_bytes() == before
    assert tuning.read_bytes() == before_tuning


def test_recovered_empty_file_can_be_explicitly_reviewed_and_saved(tmp_path, monkeypatch):
    _, output, tuning, args = setup_session(tmp_path, monkeypatch, recover=True)
    monkeypatch.setattr(review, "plot_review_browser", lambda *a: None)
    answers = iter(["SAVE", "", "No usable MISC"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args) == 0
    assert len(mne.read_annotations(output)) == 0
    receipt = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert receipt["status"] == "reviewed" and receipt["annotation_count"] == 0


def test_save_retains_another_run_reviewed_since_this_browser_opened(tmp_path, monkeypatch):
    """Reproduce the P04 overlap: startup cfg lacks a later completed run."""
    _, _, tuning, args = setup_session(tmp_path, monkeypatch)
    cfg = review.load_and_resolve_config(None)[0]
    assert not cfg["artifact_annotation"].get("manual_eye_reviews")

    def other_run_finishes(display, *unused):
        tuning.write_text(yaml.safe_dump({"stage4_review": {"03_A_MaxFilt": {"status": "accepted"}},
            "artifact_annotation": {"manual_eye_reviews": {"sub99_run02_test.fif": {
                "status": "reviewed", "annotation_sha256": "other_run_hash"}}}}))

    monkeypatch.setattr(review, "plot_review_browser", other_run_finishes)
    answers = iter(["SAVE", "MISC001", "Reviewed current source"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args) == 0
    saved = yaml.safe_load(tuning.read_text())
    assert set(saved["artifact_annotation"]["manual_eye_reviews"]) == {
        "sub99_run01_test.fif", "sub99_run02_test.fif"}
    assert saved["artifact_annotation"]["manual_eye_reviews"]["sub99_run02_test.fif"]["annotation_sha256"] == "other_run_hash"
    assert saved["stage4_review"]["03_A_MaxFilt"]["status"] == "accepted"


def test_orphan_recovery_keeps_saved_muscle_and_does_not_reseed_blinks(tmp_path, monkeypatch):
    """The saved edited set must survive recovery with diagnostics and no receipt."""
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch, recover=True)
    previous = mne.Annotations([raw.first_time + 2., raw.first_time + 5.], [.4, .2],
        ["BAD_muscle", "BAD_blink_auto"], orig_time=raw.annotations.orig_time)
    previous.save(output, overwrite=True)
    provisional = mne.Annotations([raw.first_time + 1.], [.1], ["BAD_muscle"],
        orig_time=raw.annotations.orig_time)
    monkeypatch.setattr(review, "load_muscle_context", lambda *args: (provisional,
        {"threshold": 10., "settings_sha256": "synthetic"}))

    def select(display, cfg, *, state):
        state["components"] = [0]
        return {"model_sha256": "synthetic_model", "selected_components": [0]}

    monkeypatch.setattr(review, "add_preliminary_sources", select)
    monkeypatch.setattr(review, "detect_blinks", lambda sources, components, settings: (
        [(600, 650)], np.zeros(raw.n_times), {"candidate_count": 1,
        "status": "candidates_for_review", "peaks": [], "settings": settings}))

    def check_existing(display, *unused):
        assert list(display.annotations.description) == ["BAD_muscle", "BAD_blink_auto"]
        np.testing.assert_allclose(display.annotations.onset, previous.onset, atol=2e-6)

    monkeypatch.setattr(review, "plot_review_browser", check_existing)
    answers = iter(["SAVE", "MISC001", "Confirmed recovered saved intervals"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args + ["--preliminary-ica"]) == 0
    saved = mne.read_annotations(output)
    assert list(saved.description) == list(previous.description)
    np.testing.assert_allclose(saved.onset, previous.onset, atol=2e-6)
    receipt = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert receipt["muscle_review"]["annotation_count"] == 1
    assert receipt["blink_auto"]["retained_auto_count"] == 1


def test_auto_proposals_are_not_saved_before_review_or_reseeded_after_reviewed_deletions(tmp_path, monkeypatch):
    """Exercise real detection and two CLI sessions, with GUI actions mocked."""
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch)
    cfg = review.load_and_resolve_config(None)[0]
    times = raw.times
    values = np.random.default_rng(4).normal(scale=.01, size=raw.n_times)
    values += 5. * np.exp(-.5 * ((times - 4.) / .08) ** 2)
    sources = mne.io.RawArray(values[None], mne.create_info(["ICA000"], 100., "misc"),
                             first_samp=raw.first_samp, verbose="ERROR")
    sources.set_meas_date(raw.info["meas_date"])
    monkeypatch.setattr(review, "load_cached_aid", lambda *a: {"native_sources": sources,
        "metadata": {"model_sha256": "synthetic_model", "fingerprint": "synthetic_cache"}})

    def select(display, cfg, *, state):
        state["components"] = [0]
        return {"model_sha256": "synthetic_model", "selected_components": [0]}

    def delete_after_inspecting(display, *unused):
        assert not output.exists()  # Automatic BAD proposals have no persisted artifact file.
        assert "BAD_blink_auto" in display.annotations.description
        display.annotations.delete(list(range(len(display.annotations))))

    monkeypatch.setattr(review, "add_preliminary_sources", select)
    monkeypatch.setattr(review, "plot_review_browser", delete_after_inspecting)
    answers = iter(["SAFE", "SAVE", "MISC001", "Reviewed entire source; deleted false proposals"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args + ["--preliminary-ica"]) == 0
    assert len(mne.read_annotations(output)) == 0
    payload = yaml.safe_load(tuning.read_text())
    saved = payload["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert saved["blink_auto"]["candidate_count"] >= 1 and saved["blink_auto"]["retained_auto_count"] == 0
    assert saved["blink_auto"]["status"] == "reviewed" and saved["blink_auto"]["removed_spans_raw"]
    cfg["artifact_annotation"]["manual_eye_reviews"] = payload["artifact_annotation"]["manual_eye_reviews"]
    before = output.read_bytes(), tuning.read_bytes()

    def unchanged_reopen(display, *unused):
        assert len(display.annotations) == 0  # Neither EDF nor ICA proposals resurrect.

    monkeypatch.setattr(review, "plot_review_browser", unchanged_reopen)
    monkeypatch.setattr("builtins.input", lambda prompt: "DISCARDALL")
    assert review.main(args + ["--preliminary-ica", "--overwrite"]) == 0
    assert (output.read_bytes(), tuning.read_bytes()) == before


def test_missing_cache_stops_before_loading_edf_or_frontal_signal_data(tmp_path, monkeypatch):
    """A missing cluster preparation must fail quickly, preserving saved work."""
    _, output, tuning, args = setup_session(tmp_path, monkeypatch)

    def missing(*unused):
        raise RuntimeError("Blink-aid cache is missing")

    monkeypatch.setattr(review, "load_cached_aid", missing)
    monkeypatch.setattr(review, "add_edf_display", lambda *a: pytest.fail("No full signal read before cache validation"))
    with pytest.raises(RuntimeError, match="cache is missing"):
        review.main(args + ["--preliminary-ica"])
    assert not output.exists() and not tuning.exists()


# -----------------------------------------------------------------------------
# Accepted timing is visual context: native coordinates, green, never persisted
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("meas_date", [None, datetime(2026, 1, 1, tzinfo=timezone.utc)])
def test_epoch_preview_uses_corrected_native_samples_and_reports_incomplete_windows(monkeypatch, meas_date):
    """An original trigger must not substitute for its accepted shifted onset."""
    raw = mne.io.RawArray(np.zeros((1, 800)), mne.create_info(["MISC001"], 100., "misc"),
                         first_samp=13000, verbose="ERROR")
    raw.set_meas_date(meas_date)
    record = AcquisitionRunRecord("source", "exact.fif", "exact", 1, "source")
    cfg = {"epoching": {"tmin_sec": -.5, "tmax_sec": 2.667},
           "resolved_paths": {"derivatives_root": "/temporary/derivatives"},
           "layout": {"subject_bids_id": "sub-99"}}
    base = {"source_fif": "exact.fif", "sfreq_hz": "100.0", "recording_first_sample_raw": "13000",
            "recording_end_sample_raw_exclusive": "13800", "condition_code": 71,
            "raw_sample_frame": 13010, "final_analysis_eligible": True}
    rows = [dict(base, epoch_event_sample_raw=13020), dict(base, epoch_event_sample_raw=13700),
            dict(base, final_analysis_eligible=False), dict(base, source_fif="another.fif", epoch_event_sample_raw=13100)]
    monkeypatch.setattr(review, "load_prepared_decisions", lambda cfg: ["unchanged ledger"])

    def load(cfg, prepared, path):
        assert prepared == ["unchanged ledger"]
        assert str(path) == "/temporary/derivatives/visual_timing/s3-v002/reviews/sub-99_review.json"
        return rows, {"version": "s3-v002", "review_sha256": "accepted_hash"}

    monkeypatch.setattr(review, "load_reviewed_timing", load)
    monkeypatch.setattr(review, "gaze_break_guidance", lambda *a: {"markers": [], "onsets": {}, "issues": []})
    result = review.load_timing_preview(cfg, record, raw)
    assert result["events"].tolist() == [[13020, 0, 71], [13700, 0, 71]]
    np.testing.assert_allclose(result["annotations"].onset, [130., 136.5])
    np.testing.assert_allclose(result["annotations"].duration, [2.87, 1.49])
    assert result["annotations"].orig_time == raw.annotations.orig_time
    assert result["metadata"]["incomplete_window_count"] == 2
    assert result["metadata"]["without_accepted_onset_count"] == 1
    assert rows[0]["raw_sample_frame"] == 13010
    rows[0]["sfreq_hz"] = "200.0"
    with pytest.raises(RuntimeError, match="coordinates differ"):
        review.load_timing_preview(cfg, record, raw)


def test_timing_preview_propagates_acceptance_failure_and_honors_explicit_path(monkeypatch):
    """Never silently replace a failed accepted map with uncorrected triggers."""
    cfg = {"epoching": {}, "resolved_paths": {}, "layout": {}}
    monkeypatch.setattr(review, "load_prepared_decisions", lambda cfg: [])

    def stale(cfg, rows, path):
        assert path == "explicit_review.json"
        raise RuntimeError("Reviewed Stage 3 manifest hash mismatch")

    monkeypatch.setattr(review, "load_reviewed_timing", stale)
    with pytest.raises(RuntimeError, match="hash mismatch"):
        review.load_timing_preview(cfg, None, None, "explicit_review.json")


def gaze_preview_example(tmp_path, monkeypatch, *, missing_timing=False, unmatched=False):
    """Create real canonical events and hash-bound timing for a replay example."""
    from downstream_stage_utils import REQUIRED_TRIGGER_COLUMNS
    from continuous_cleaning import file_hash
    raw = mne.io.RawArray(np.zeros((1, 1000)), mne.create_info(["MISC001"], 100., "misc"),
                         first_samp=13000, verbose="ERROR")
    record = AcquisitionRunRecord("source", "exact.fif", "exact", 1, "source")
    pre = tmp_path / "preprocessing"; (pre / "events").mkdir(parents=True)
    cfg = {"subject": {"id": "99"}, "layout": {"subject_bids_id": "sub-99"},
        "epoching": {"tmin_sec": -.5, "tmax_sec": 2.667},
        "resolved_paths": {"derivatives_root": str(tmp_path), "derivatives_preprocessing_dir": str(pre),
                           "derivatives_reports_dir": str(tmp_path / "reports")}}
    base = {"source_fif": "exact.fif", "sfreq_hz": "100.0", "recording_first_sample_raw": "13000",
            "recording_end_sample_raw_exclusive": "14000", "condition_code": 71, "gaze_break_evidence": False,
            "final_analysis_eligible": True, "behavioral_block": 1, "behavioral_run": 1, "behavioral_trial": 2}
    rows = [dict(base, trial_key="first", raw_sample_frame=13100, trial_interval_end_raw_exclusive=13300,
                 epoch_event_sample_raw=13102, evidence_trigger_codes="71", evidence_event_indices="1", evidence_samples_raw="13100"),
            dict(base, trial_key="interrupted", raw_sample_frame=13300, trial_interval_end_raw_exclusive=13600,
                 final_analysis_eligible=False, gaze_break_evidence=True, evidence_trigger_codes="71,150",
                 evidence_event_indices="3,4", evidence_samples_raw="13300,13380"),
            dict(base, trial_key="replay", raw_sample_frame=13600, trial_interval_end_raw_exclusive=14000,
                 epoch_event_sample_raw=13602, is_replay_trial=True, evidence_trigger_codes="71",
                 evidence_event_indices="6", evidence_samples_raw="13600")]
    marker = {name: "" for name in REQUIRED_TRIGGER_COLUMNS | {"event_index"}}
    marker.update(source_fif="exact.fif", trigger_code="150", event_index="4", raw_sample_frame="13380",
                  unresolved_transition="False", is_gaze_break="True")
    markers = [marker, dict(marker, source_fif="different_acquisition.fif")]
    if unmatched:
        markers.append(dict(marker, event_index="99", raw_sample_frame="13800"))
    with (pre / "events/02_trigger_events.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(marker)); writer.writeheader(); writer.writerows(markers)
    entries = [dict(row, assignment_status="proposed", assignment_basis="estimated",
                    proposed_sample_raw=row["raw_sample_frame"] + 2, applied_shift_samples=2, lag_estimate_ms=20.) for row in rows]
    if missing_timing:
        entries[1]["proposed_sample_raw"] = ""
    map_path = tmp_path / "trial_timing_map.csv"
    fields = sorted(set().union(*(entry.keys() for entry in entries)))
    with map_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(entries)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"output_sha256": {map_path.name: file_hash(map_path)}}))
    accepted = tmp_path / "review.json"
    accepted.write_text(json.dumps({"bundle_manifest": manifest.name, "accepted_assignment_bases": ["estimated"]}))
    gate = {"review_path": str(accepted), "review_sha256": file_hash(accepted),
            "bundle_manifest_sha256": file_hash(manifest), "version": "synthetic-reviewed"}
    monkeypatch.setattr(review, "load_prepared_decisions", lambda cfg: rows)
    monkeypatch.setattr(review, "load_reviewed_timing", lambda *a: (rows, gate))
    return cfg, record, raw, rows, gate


@pytest.mark.parametrize("missing_timing, unmatched", [(False, False), (True, False), (False, True)])
def test_gaze_break_uses_exact_trial_and_keeps_native_marker_unshifted(tmp_path, monkeypatch, missing_timing, unmatched):
    cfg, record, raw, rows, gate = gaze_preview_example(tmp_path, monkeypatch,
        missing_timing=missing_timing, unmatched=unmatched)
    result = review.load_timing_preview(cfg, record, raw)
    valid = result["annotations"][result["annotations"].description == review.EPOCH_LABEL]
    invalid = result["annotations"][result["annotations"].description == review.INVALID_EPOCH_LABEL]
    np.testing.assert_allclose(valid.onset, [130.52, 135.52])
    assert len(invalid) == (0 if missing_timing else 1)
    if not missing_timing:
        np.testing.assert_allclose(invalid.onset, [132.52])
        np.testing.assert_allclose(invalid.duration, [3.17])  # Inclusive native endpoints, same reviewed rounding.
    events = result["events"]
    assert [13380, 0, 150] in events.tolist() and [13382, 0, 150] not in events.tolist()
    assert not rows[1]["final_analysis_eligible"] and "epoch_event_sample_raw" not in rows[1]
    assert rows[2]["final_analysis_eligible"]  # The replay keeps its own green guide.
    assert result["metadata"]["gaze_break_marker_count"] == 1 + unmatched
    assert result["metadata"]["corrected_onset_count"] == 2
    assert bool(result["metadata"]["gaze_break_review_issues"]) == (missing_timing or unmatched)


def test_gaze_break_rejects_changed_timing_and_reports_ambiguous_ownership(tmp_path, monkeypatch):
    cfg, record, raw, rows, gate = gaze_preview_example(tmp_path, monkeypatch)
    ambiguous = rows + [dict(rows[1], trial_key="ambiguous_duplicate")]
    guidance = review.gaze_break_guidance(cfg, record, ambiguous, gate)
    assert not guidance["onsets"] and guidance["issues"][0]["candidate_count"] == 2
    rows[1]["final_analysis_eligible"] = True
    with pytest.raises(RuntimeError, match="contradicts"):
        review.gaze_break_guidance(cfg, record, rows, gate)
    rows[1]["final_analysis_eligible"] = False
    path = tmp_path / "trial_timing_map.csv"
    path.write_text(path.read_text().replace("13302", "13303"))
    with pytest.raises(RuntimeError, match="timing map changed"):
        review.gaze_break_guidance(cfg, record, rows, gate)


def test_orange_gaze_marker_is_on_top_and_invalid_preview_never_saves(tmp_path, monkeypatch):
    """Exercise native rendering, visibility toggles, scrolling and persistence."""
    import matplotlib.pyplot as plt
    from matplotlib.backend_bases import KeyEvent
    from matplotlib.colors import to_rgba
    from mne.viz import utils
    mne.viz.set_browser_backend("matplotlib")
    raw = mne.io.RawArray(np.sin(np.arange(1200)[None] / 20),
        mne.create_info(["MISC001"], 100., "misc"), first_samp=13000, verbose="ERROR")
    raw.annotations.append(135., .1, "BAD_blink")
    timing = {"events": np.array([[13100, 0, 71], [13200, 0, 81], [13250, 0, 150]]),
        "annotations": mne.Annotations([130.5, 131.5], [3.167] * 2, [review.EPOCH_LABEL, review.INVALID_EPOCH_LABEL]),
        "onset_linewidth": 3., "onset_color": "#0057b8"}
    original_plot = mne.io.BaseRaw.plot
    figures = []

    def narrow_view(self, **kwargs):
        kwargs.update(duration=4.)
        return original_plot(self, **kwargs)

    def inspect_browser(*unused, **kwargs):
        fig = kwargs.get("fig")
        if fig is None or not kwargs.get("block"):
            return
        figures.append(fig)
        assert fig.mne.annotation_segment_colors[review.EPOCH_LABEL] == review.EPOCH_COLOR
        assert fig.mne.annotation_segment_colors[review.INVALID_EPOCH_LABEL] == review.GAZE_BREAK_COLOR
        np.testing.assert_allclose(fig.mne.event_lines.get_colors(),
            [to_rgba("#0057b8"), to_rgba("#0057b8"), to_rgba(review.GAZE_BREAK_COLOR)])
        overlay = next(artist for artist in fig.mne.ax_main.collections if artist.get_label() == "DAD gaze-break markers")
        assert len(overlay.get_segments()) == 1 and overlay.get_alpha() == 1.
        np.testing.assert_allclose(overlay.get_linewidths(), 3.)
        assert overlay.get_zorder() > max(fig.mne.zorder.values())
        assert overlay.get_segments()[0][0, 0] == 2.5  # No lag added to code 150.
        fig.canvas.callbacks.process("key_press_event", KeyEvent("key_press_event", fig.canvas, key="a"))
        checkboxes = fig.mne.show_hide_annotation_checkboxes
        names = [label.get_text() for label in checkboxes.labels]
        checkboxes.set_active(names.index(review.INVALID_EPOCH_LABEL))
        assert not fig.mne.visible_annotations[review.INVALID_EPOCH_LABEL]
        assert fig.mne.visible_annotations[review.EPOCH_LABEL]
        for key, expected_count in (("shift+right", 0), ("shift+left", 1)):
            fig.canvas.callbacks.process("key_press_event", KeyEvent("key_press_event", fig.canvas, key=key))
            fig.canvas.draw()
            assert len(overlay.get_segments()) == expected_count
        raw.annotations.append(137., .1, "BAD_segment")

    monkeypatch.setattr(mne.io.BaseRaw, "plot", narrow_view)
    monkeypatch.setattr(utils, "_show_browser", inspect_browser)
    try:
        review.plot_review_browser(raw, ["MISC001"], timing=timing)
        assert not set(raw.annotations.description).intersection(review.PREVIEW_COLORS)
        path = tmp_path / "manual-annot.fif"
        review.eye_annotations(raw).save(path)
        assert set(mne.read_annotations(path).description) == {"BAD_blink", "BAD_segment"}
    finally:
        for fig in figures:
            if fig.mne.fig_annotation is not None:
                plt.close(fig.mne.fig_annotation)
            plt.close(fig)


@pytest.mark.parametrize("marker_color", [None, "#0057b8"])
def test_epoch_guides_are_green_and_removed_after_browser_artifact_edits(monkeypatch, marker_color):
    """Exercise the actual Matplotlib browser and color assignment on old MNE."""
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    mne.viz.set_browser_backend("matplotlib")
    raw = mne.io.RawArray(np.sin(np.arange(800)[None] / 20),
        mne.create_info(["MISC001"], 100., "misc"), first_samp=13000, verbose="ERROR")
    raw.annotations.append(131., .2, "EDF_blink_candidate")
    timing = {"events": np.array([[13200, 0, 71]]), "annotations":
        mne.Annotations([131.5], [3.167], [review.EPOCH_LABEL], orig_time=raw.annotations.orig_time)}
    expected_color = marker_color or review.DEFAULT_ONSET_COLOR
    if marker_color:
        timing["onset_color"] = marker_color
    original_plot = mne.io.BaseRaw.plot
    original_palette = rcParams["axes.prop_cycle"]
    figures = []

    def nonblocking_plot(self, **kwargs):
        assert kwargs["events"].tolist() == [[13200, 0, 71]]
        assert kwargs["event_color"] == expected_color
        kwargs.update(block=False, show=False, duration=3)
        fig = original_plot(self, **kwargs)
        figures.append(fig)
        assert fig.mne.annotation_segment_colors[review.EPOCH_LABEL] == review.EPOCH_COLOR
        assert fig.mne.visible_annotations[review.EPOCH_LABEL]
        assert all(region.get_alpha() == .3 for region in fig.mne.annotations)
        # Simulate resizing a guide and adding a true artifact in the same view.
        idx = list(self.annotations.description).index(review.EPOCH_LABEL)
        self.annotations.onset[idx] += .1
        self.annotations.append(132., .1, "BAD_blink")
        return fig

    monkeypatch.setattr(mne.io.BaseRaw, "plot", nonblocking_plot)
    try:
        review.plot_review_browser(raw, ["MISC001"], timing=timing)
        assert review.EPOCH_LABEL not in raw.annotations.description
        assert list(review.eye_annotations(raw).description) == ["BAD_blink"]
        assert rcParams["axes.prop_cycle"] == original_palette
        assert timing["annotations"].onset[0] == 131.5  # Fresh opening restores accepted coordinates.
        fig = figures[0]
        np.testing.assert_allclose(fig.mne.event_lines.get_linewidths(), 2.)
        from matplotlib.colors import to_rgba
        np.testing.assert_allclose(fig.mne.event_lines.get_colors(), [to_rgba(expected_color)])
        assert fig.mne.event_lines.get_alpha() == 1.
        assert fig.mne.event_lines.get_zorder() > max(region.get_zorder() for region in fig.mne.annotations)
        trace_widths = [line.get_linewidth() for line in fig.mne.traces]
        # The normal scroll key recreates event artists; styling must persist.
        # Our mock appended BAD_blink directly; perform the label registration
        # that MNE's actual interactive annotation controls normally do.
        fig._setup_annotation_colors()
        from matplotlib.backend_bases import KeyEvent
        fig.canvas.callbacks.process("key_press_event", KeyEvent("key_press_event", fig.canvas, key="right"))
        fig.canvas.draw()
        np.testing.assert_allclose(fig.mne.event_lines.get_linewidths(), 2.)
        np.testing.assert_allclose(fig.mne.event_lines.get_colors(), [to_rgba(expected_color)])
        assert fig.mne.event_lines.get_alpha() == 1.
        assert [line.get_linewidth() for line in fig.mne.traces] == trace_widths
    finally:
        for fig in figures:
            plt.close(fig)


def test_no_timing_marks_skips_the_loader_and_saved_files_never_contain_previews(tmp_path, monkeypatch):
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch)
    calls = []
    preview = {"events": np.array([[raw.first_samp + 200, 0, 71]]),
        "annotations": mne.Annotations([raw.first_time + 1.5], [3.167], [review.EPOCH_LABEL],
                                        orig_time=raw.annotations.orig_time),
        "metadata": {"enabled": True, "review_sha256": "accepted_hash"}}
    monkeypatch.setattr(review, "load_timing_preview", lambda *a: calls.append(1) or preview)

    def browser(display, misc, edf, timing):
        assert timing is preview
        # Even if a browser failed to remove its non-BAD guides, SAVE excludes them.
        display.annotations.append(preview["annotations"].onset, preview["annotations"].duration, preview["annotations"].description)
        display.annotations.append(display.first_time + 1., .1, "BAD_blink")

    monkeypatch.setattr(review, "plot_review_browser", browser)
    answers = iter(["SAVE", "MISC001", "Whole source inspected with timing guides"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args) == 0
    assert review.EPOCH_LABEL not in mne.read_annotations(output).description
    saved = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]["sub99_run01_test.fif"]
    assert saved["timing_preview"]["review_sha256"] == "accepted_hash"
    monkeypatch.setattr(review, "plot_review_browser", lambda display, misc, edf, timing: pytest.fail("unexpected") if timing is not None else None)
    monkeypatch.setattr("builtins.input", lambda prompt: "DISCARDALL")
    assert review.main(args + ["--overwrite", "--no-timing-marks"]) == 0
    assert calls == [1]


@pytest.mark.parametrize("width", ["0", "-1", "nan", "inf"])
def test_invalid_onset_width_is_rejected_before_reading_data(width):
    with pytest.raises(SystemExit) as error:
        review.main(["--config", "unused", "--subject-config", "unused", "--tuning-config", "unused",
                     "--source", "unused", "--basis", "unused", "--onset-linewidth", width])
    assert error.value.code == 2


@pytest.mark.parametrize("color", ["none", "#00000080", "#00000000", "not-a-color", "#12345"])
def test_invalid_or_transparent_onset_colors_are_rejected_before_loading(color):
    with pytest.raises(SystemExit) as error:
        review.main(["--config", "unused", "--subject-config", "unused", "--tuning-config", "unused",
                     "--source", "unused", "--basis", "unused", "--onset-color", color])
    assert error.value.code == 2


@pytest.mark.parametrize("value, expected", [("black", "#000000"), ("tab:blue", "#1f77b4"),
                                           ("#0057B8", "#0057b8"), ("#000000ff", "#000000")])
def test_named_and_hex_onset_colors_are_portable_opaque_rgb(value, expected):
    assert review.onset_color(value) == expected


@pytest.mark.parametrize("existing_labels", [False, True])
def test_default_descriptions_are_available_but_unused_labels_are_not_saved(tmp_path, monkeypatch, existing_labels):
    """Exercise native selector creation with and without existing BAD labels."""
    import matplotlib.pyplot as plt
    from matplotlib.backend_bases import KeyEvent
    from mne.viz import utils
    from continuous_cleaning import annotation_rows, good_stretches
    mne.viz.set_browser_backend("matplotlib")
    raw = mne.io.RawArray(np.zeros((1, 800)), mne.create_info(["MISC001"], 100., "misc"),
                         first_samp=13000, verbose="ERROR")
    if existing_labels:
        raw.annotations.append([131., 133.], [.1, .1], ["BAD_blink", "BAD_blink_auto"])
    original = annotation_rows(raw)
    stretches = good_stretches(raw)
    figures = []

    def inspect_browser(*unused, **kwargs):
        fig = kwargs.get("fig")
        if fig is None or not kwargs.get("block"):
            return
        figures.append(fig)
        assert set(review.REVIEW_DESCRIPTIONS).issubset(fig._get_annotation_labels())
        assert annotation_rows(raw) == original and good_stretches(raw) == stretches
        # Re-registering must neither duplicate choices nor add dummy intervals.
        review.prepare_annotation_descriptions(fig)
        assert len(fig.mne.new_annotation_labels) == len(set(fig.mne.new_annotation_labels))
        fig.canvas.callbacks.process("key_press_event", KeyEvent("key_press_event", fig.canvas, key="a"))
        labels = [label.get_text() for label in fig.mne.fig_annotation.mne.radio_ax.buttons.labels]
        assert set(review.REVIEW_DESCRIPTIONS).issubset(labels)
        assert len(labels) == len(set(labels))
        unused_path = tmp_path / "unused-annot.fif"
        review.eye_annotations(raw).save(unused_path)
        saved = mne.read_annotations(unused_path)
        assert "BAD_segment" not in saved.description and "BAD_eye_movement" not in saved.description
        assert len(saved) == len(original["rows"])
        # Once a genuine interval is marked, the same selector label is saved.
        raw.annotations.append(raw.first_time + 5., .2, "BAD_segment")
        used_path = tmp_path / "used-annot.fif"
        review.eye_annotations(raw).save(used_path)
        used = mne.read_annotations(used_path)
        assert sum(used.description == "BAD_segment") == 1
        assert "BAD_eye_movement" not in used.description

    monkeypatch.setattr(utils, "_show_browser", inspect_browser)
    try:
        review.plot_review_browser(raw, ["MISC001"])
        assert len(figures) == 1
    finally:
        for fig in figures:
            if fig.mne.fig_annotation is not None:
                plt.close(fig.mne.fig_annotation)
            plt.close(fig)


@pytest.mark.parametrize("empty", [False, True])
def test_muscle_edits_and_reviewed_zero_survive_save_and_reopening(tmp_path, monkeypatch, empty):
    """Use real provisional report/FIF evidence and real local persistence."""
    from continuous_cleaning import digest, file_stamp, muscle_parameters
    raw, output, tuning, args = setup_session(tmp_path, monkeypatch)
    cfg = review.load_and_resolve_config(None)[0]
    cfg["artifact_annotation"].update(muscle_threshold=10., muscle_channel_type="mag",
        muscle_filter_freq_hz=[110., 140.], muscle_min_length_good_sec=.2)
    dirs = review.stage_output_dirs(cfg)
    source = dirs["sss"] / "sub99_run01_test_raw_sss.fif"
    signature = digest({"input": file_stamp(source), **muscle_parameters(cfg["artifact_annotation"], 10.), "eye_review": None})
    provisional = raw.copy()
    provisional.annotations.append([raw.first_time + 1., raw.first_time + 3.], [.1, .2], ["BAD_muscle"] * 2)
    provisional.info["description"] = f"stage4_input_sha256={signature}"
    candidate = dirs["sss"] / "sub99_run01_test_ann.fif"
    provisional.save(candidate, overwrite=True, verbose="ERROR")
    report = {"runs": [{"source_filename": "sub99_run01_test.fif", "input_stamp": file_stamp(source), "muscle_threshold": 10.}],
              "evidence_files": [{"relative_path": "sss/" + candidate.name, "stamp": file_stamp(candidate)}]}
    (dirs["reports"] / "03_B_ArtefactAnnotation_review_report.json").write_text(json.dumps(report))

    def edit(display, *unused):
        assert sum(label == "BAD_muscle" for label in display.annotations.description) == 2
        display.annotations.delete([i for i, label in enumerate(display.annotations.description) if label == "BAD_muscle"])
        if not empty:
            display.annotations.append([display.first_time + 1.5, display.first_time + 6.], [1.2, .3], ["BAD_muscle"] * 2)

    monkeypatch.setattr(review, "plot_review_browser", edit)
    answers = iter(["SAVE", "MISC001", "Reviewed muscle edits as well as eyes"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    assert review.main(args) == 0
    saved = mne.read_annotations(output)
    expected = [] if empty else [raw.first_time + 1.5, raw.first_time + 6.]
    np.testing.assert_allclose(saved.onset[saved.description == "BAD_muscle"], expected, atol=2e-6)
    reviews = yaml.safe_load(tuning.read_text())["artifact_annotation"]["manual_eye_reviews"]
    receipt = reviews["sub99_run01_test.fif"]
    assert receipt["muscle_review"]["status"] == "reviewed"
    assert receipt["muscle_review"]["annotation_count"] == len(expected)
    cfg["artifact_annotation"]["manual_eye_reviews"] = reviews
    record = review.discover_raw_run_records(cfg)[0]
    loaded, _ = review.load_manual(cfg, record, raw, source)
    assert len(loaded) == len(saved)
    receipt["muscle_review"]["annotation_count"] += 1
    with pytest.raises(RuntimeError, match="muscle count"):
        review.load_manual(cfg, record, raw, source)
    receipt["muscle_review"]["annotation_count"] -= 1

    def reopen(display, *unused):
        # The old 1.0/3.0 second automatic intervals must not be reintroduced.
        np.testing.assert_allclose(display.annotations.onset[display.annotations.description == "BAD_muscle"], expected, atol=2e-6)

    monkeypatch.setattr(review, "plot_review_browser", reopen)
    monkeypatch.setattr("builtins.input", lambda prompt: "DISCARDALL")
    before = output.read_bytes(), tuning.read_bytes()
    assert review.main(args + ["--overwrite"]) == 0
    assert before == (output.read_bytes(), tuning.read_bytes())
    cfg["artifact_annotation"]["muscle_channel_type"] = "grad"
    with pytest.raises(RuntimeError, match="settings changed"):
        review.main(args + ["--overwrite"])
